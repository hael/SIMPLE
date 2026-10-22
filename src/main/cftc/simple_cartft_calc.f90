!@descr: Cartesian Fourier calculator of continuous pose refinement: references, prepared particles, objectives
! The Cartesian counterpart of polarft_calc: state-by-half references and the prepared particle
! slots of the current batch. Evaluation methods take the calculator intent(in); set_ref runs
! serially, set_ptcl writes only its own slot. A slot is prepared for one objective: OBJFUN_CC
! (1 - cc, uniform pixel weight, no sigma2) or OBJFUN_EUCLID (sigma2-weighted loss L, score exp(-L)).
! The objectives and the reference file are specified in refine3D_policy.md (Cartesian section).
module simple_cartft_calc
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
use, intrinsic :: iso_fortran_env, only: int64
use simple_core_module_api,   only: dp, sp, PI, KBWINSZ, KBALPHA, OSMPL_PAD_FAC, CTFFLAG_NO, CTFFLAG_YES, CTFFLAG_FLIP, &
    &OBJFUN_CC, OBJFUN_EUCLID, ctfparams, ctfvars, kbinterpol, cyci_1d, simple_exception, string, fopen, fclose, &
    &fileiochk, file_exists, logfhandle
use simple_image,             only: image
use simple_ctf,               only: ctf
use simple_gridding,          only: kb_stencil_envelope_1d
use simple_cartesian_fourier, only: center_embed_real3d, gather_packed_bspline3_window_grad, bspline3_stencil_envelope_1d
implicit none
private
public :: cartft_calc, CART_REFVOLS_FORMAT_VERSION

#include "simple_local_flags.inc"

real(dp), parameter :: NUMERIC_FLOOR = epsilon(1._dp)**2
! range of the regression slope of data on model that scales the Gauss-Newton matrices
real(dp), parameter :: CURVATURE_SCALE_MIN = 0.05_dp, CURVATURE_SCALE_MAX = 20._dp
integer,  parameter :: CART_REFVOLS_FORMAT_VERSION = 1 !< format_version of the reference-volume files (6.6)
integer,  parameter :: CART_REFVOLS_NHEADER = 5        !< [format_version, box_crop, nstates, kfrom, kto]

!> One reference: the padded Fourier transform of a prepared physical volume, weighted by
!! ref_weight, times padf**3 (matcher side), or the staged prepared volume itself (materializer
!! side, written by write).
type :: cartft_ref
    complex, allocatable :: cmat(:,:,:)
    real,    allocatable :: vol(:,:,:)
end type cartft_ref

!> One prepared particle for one objective: the observation X and the shift-free transfer
!! abs(C) over its shells, whitened by 1/sqrt(sigma2) of its half under OBJFUN_EUCLID,
!! unweighted under OBJFUN_CC.
type :: cartft_ptcl
    complex, allocatable :: observed(:,:)
    real,    allocatable :: transfer(:,:)
    real(dp) :: power   = 0._dp                 !< sum |X|^2 (whitened under OBJFUN_EUCLID)
    integer  :: kfromto(2) = 0
    integer  :: npix    = 0                     !< pixels of its shells in the full disk
    integer  :: half    = 0                     !< column of its half in sigma2 (OBJFUN_EUCLID only)
    integer  :: objfun  = 0                     !< OBJFUN_CC or OBJFUN_EUCLID once prepared
    logical  :: valid   = .false.
end type cartft_ptcl

type :: cartft_calc
    private
    integer :: nstates = 0                  !< number of states
    integer :: nptcls  = 0                  !< number of particle slots
    integer :: box     = 0                  !< (cropped) box of particles and references
    integer :: boxpd   = 0                  !< padded box of the references
    integer :: padf    = 1                  !< padding factor
    integer :: iwinsz  = 0                  !< integer half-width of the gather stencil
    integer :: lims2(2,2) = 0               !< full-disk 2D Fourier limits
    real    :: padsc   = 1.                 !< padf**3, the native Fourier scale folded into the references
    type(kbinterpol)               :: kbwin !< KB window of the polar branch (particle taper, ref_weight)
    real,              allocatable :: ref_weight(:)  !< (boxpd) per-axis weight of the padded reference (set_ref)
    integer,           allocatable :: wrap(:)     !< periodic wrap table of the padded grid
    type(cartft_ref),  allocatable :: refs(:,:)   !< (2,nstates): 1 even, 2 odd
    type(cartft_ptcl), allocatable :: ptcls(:)    !< (nptcls) prepared particles
    integer,           allocatable :: slots(:)    !< particle index -> slot of the current batch (0 when absent)
    integer,           allocatable :: shells(:,:) !< (lims2) shell nint(sqrt(h*h+k*k)) of a pixel inside the Nyquist disk, -1 outside
    integer,           allocatable :: shell_npix(:) !< (0:box/2) pixels of each shell in the full disk
    real,              allocatable :: sigma2(:,:)       !< (shell, half) the two global spectra (OBJFUN_EUCLID)
    real,              allocatable :: rsqrt_sigma2(:,:) !< 1/sqrt(sigma2) where valid
    logical,           allocatable :: sigma2_ok(:,:)    !< finite and positive
    logical :: exists = .false.
contains
    ! lifecycle
    procedure          :: new
    procedure          :: kill
    ! references
    procedure          :: set_ref
    procedure          :: ref_exists
    procedure          :: has_refs
    procedure          :: get_ref_bytes
    ! reference-volume file (6.6)
    procedure          :: set_refvol
    procedure          :: write
    procedure          :: read
    procedure          :: cart_refvols_header_compatible
    ! particles
    procedure          :: set_sigma2
    procedure, private :: set_ptcl_cc
    procedure, private :: set_ptcl_euclid
    generic            :: set_ptcl => set_ptcl_cc, set_ptcl_euclid
    procedure          :: ptcl_is_valid
    procedure          :: get_ptcl_kfromto
    procedure          :: get_ptcl_npix
    procedure          :: get_ptcl_objfun
    procedure          :: get_ptcl_taper
    procedure          :: get_box
    procedure          :: new_ptcls
    procedure          :: get_nptcls
    procedure          :: set_ptcl_inds
    procedure          :: get_ptcl_slot
    ! evaluation
    procedure          :: predict
    procedure          :: objective_gradient
    procedure          :: shift_normal_terms
    procedure          :: pose_normal_terms
    procedure          :: shift_normal_terms_full_disk
    procedure          :: pose_normal_terms_full_disk
    procedure          :: sigma_contribution
    procedure          :: sigma_from_shell_residual
    procedure          :: score
    procedure, private :: fill_ptcl
    procedure, private :: gather
    procedure, private :: shift_phase_tables
    procedure, private :: ref_index
    procedure, private :: check_eval_args
end type cartft_calc

contains

    ! LIFECYCLE

    !> Allocate nstates even/odd reference slots and nptcls particle slots for an even box.
    !! A previous instance is killed first. THROW_HARD for nstates < 1, nptcls < 1 or a box
    !! that is not even and >= 2.
    subroutine new( self, nstates, box, nptcls )
        class(cartft_calc), intent(inout) :: self
        integer,            intent(in)    :: nstates, box, nptcls
        real, allocatable :: kb_env(:), bspline_env(:)
        integer :: h, k
        call self%kill
        if( nstates < 1 ) THROW_HARD('cartft_calc requires at least one state')
        if( nptcls  < 1 ) THROW_HARD('cartft_calc requires at least one particle slot')
        if( box < 2 .or. mod(box,2) /= 0 ) THROW_HARD('cartft_calc requires an even box')
        self%nstates = nstates
        self%nptcls  = nptcls
        self%box     = box
        self%padf    = OSMPL_PAD_FAC
        self%boxpd   = self%padf*box
        self%padsc   = real(self%padf)**3
        self%kbwin   = kbinterpol(KBWINSZ, KBALPHA)
        self%iwinsz  = ceiling(self%kbwin%get_winsz() - 0.5)
        ! the KB stencil's real-space response over the B-spline gather's (on-grid samples, padded period)
        call kb_stencil_envelope_1d(self%kbwin, self%boxpd, kb_env)
        call bspline3_stencil_envelope_1d(self%boxpd, bspline_env)
        self%ref_weight = kb_env/bspline_env
        self%lims2(1,:) = [-box/2, box/2]
        self%lims2(2,:) = [-box/2, box/2]
        allocate(self%refs(2,nstates), self%ptcls(nptcls))
        ! the shell of every pixel, once for the box: the pixel loops look it up instead of a square root.
        ! A pixel belongs to shell nint(|(h,k)|), the polar ring membership, only inside the Nyquist circle
        ! |(h,k)| <= box/2: the outer half of the Nyquist shell is beyond the band of the padded reference
        ! transform (and beyond its wrap table), where the polar ring of radius box/2 does not sample either
        allocate(self%shells(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2)), source=-1)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k <= (box/2)**2 ) self%shells(h,k) = nint(sqrt(real(h*h + k*k)))
            end do
        end do
        allocate(self%shell_npix(0:box/2), source=0)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( self%shells(h,k) >= 0 ) self%shell_npix(self%shells(h,k)) = self%shell_npix(self%shells(h,k)) + 1
            end do
        end do
        self%exists = .true.
    end subroutine new

    !> Release every reference and particle; safe on a killed or never-built instance.
    subroutine kill( self )
        class(cartft_calc), intent(inout) :: self
        if( allocated(self%refs)  ) deallocate(self%refs)
        if( allocated(self%ptcls) ) deallocate(self%ptcls)
        if( allocated(self%wrap)  ) deallocate(self%wrap)
        if( allocated(self%slots) ) deallocate(self%slots)
        if( allocated(self%shells)) deallocate(self%shells)
        if( allocated(self%shell_npix)   ) deallocate(self%shell_npix)
        if( allocated(self%ref_weight)   ) deallocate(self%ref_weight)
        if( allocated(self%sigma2)       ) deallocate(self%sigma2)
        if( allocated(self%rsqrt_sigma2) ) deallocate(self%rsqrt_sigma2)
        if( allocated(self%sigma2_ok)    ) deallocate(self%sigma2_ok)
        self%nstates = 0
        self%nptcls  = 0
        self%box     = 0
        self%boxpd   = 0
        self%padf    = 1
        self%iwinsz  = 0
        self%lims2   = 0
        self%padsc   = 1.
        self%exists  = .false.
    end subroutine kill

    ! REFERENCES

    !> Build the reference of one state and half from a prepared physical volume (box**3, the
    !! volume the polar branch pads and reprojects): pad by OSMPL_PAD_FAC, weight by ref_weight
    !! and Fourier transform once. THROW_HARD for a state outside 1..nstates or not box**3.
    subroutine set_ref( self, state, iseven, volume )
        class(cartft_calc), intent(inout) :: self
        integer,            intent(in)    :: state
        logical,            intent(in)    :: iseven
        real,               intent(in)    :: volume(:,:,:)
        type(image)       :: padded_image
        real, allocatable :: padded(:,:,:)
        integer :: lims3(3,2), wlims(2), lo, hi, i, j, k, ihalf
        if( .not. self%exists ) THROW_HARD('cartft_calc reference set before new')
        if( state < 1 .or. state > self%nstates ) THROW_HARD('cartft_calc reference state out of range')
        if( any(shape(volume) /= self%box) ) THROW_HARD('cartft_calc reference volume does not match the box')
        ihalf = self%ref_index(iseven)
        call padded_image%new([self%boxpd, self%boxpd, self%boxpd], 1.0)
        ! on grid, the B-spline gather of this transform equals the KB gather of the unweighted one:
        ! the reference carries the polar branch's stencil response, as the particle taper does (C10, O5)
        padded = center_embed_real3d(volume, self%boxpd)
        do k = 1, self%boxpd
            do j = 1, self%boxpd
                padded(:,j,k) = padded(:,j,k)*self%ref_weight*(self%ref_weight(j)*self%ref_weight(k))
            end do
        end do
        call padded_image%set_rmat(padded, .false.)
        deallocate(padded)
        call padded_image%fft()
        if( .not. allocated(self%wrap) )then
            ! the periodic wrap table depends on the padded box only
            lims3 = padded_image%loop_lims(3)
            wlims = lims3(2,:)
            lo    = wlims(1) - self%iwinsz - 1
            hi    = wlims(2) + self%iwinsz + 1
            allocate(self%wrap(lo:hi))
            do i = lo, hi
                self%wrap(i) = cyci_1d(wlims, i)
            end do
        endif
        ! padf**3 folded in: a power-of-two scale is exact, and the gather saves the multiply
        self%refs(ihalf,state)%cmat = self%padsc*padded_image%get_cmat()
        call padded_image%kill
    end subroutine set_ref

    !> True when the reference of this state and half has been set.
    pure logical function ref_exists( self, state, iseven )
        class(cartft_calc), intent(in) :: self
        integer,            intent(in) :: state
        logical,            intent(in) :: iseven
        ref_exists = .false.
        if( .not. self%exists ) return
        if( state < 1 .or. state > self%nstates ) return
        ref_exists = allocated(self%refs(self%ref_index(iseven),state)%cmat)
    end function ref_exists

    ! REFERENCE-VOLUME FILE (plan section 6.6)

    !> Stage the prepared physical volume (box**3) of one state and half for write, without a
    !! transform (the materializer needs no padded reference). THROW_HARD for a state outside
    !! 1..nstates or a volume that is not box**3.
    subroutine set_refvol( self, state, iseven, volume )
        class(cartft_calc), intent(inout) :: self
        integer,            intent(in)    :: state
        logical,            intent(in)    :: iseven
        real,               intent(in)    :: volume(:,:,:)
        if( .not. self%exists ) THROW_HARD('cartft_calc reference volume staged before new')
        if( state < 1 .or. state > self%nstates ) THROW_HARD('cartft_calc reference volume state out of range')
        if( any(shape(volume) /= self%box) ) THROW_HARD('cartft_calc reference volume does not match the box')
        self%refs(self%ref_index(iseven),state)%vol = volume
    end subroutine set_refvol

    !> Write the staged volumes of every state, one file per half (fname_even, fname_odd),
    !! replacing existing files. Stream binary: five default integers [format_version,
    !! box_crop, nstates, kfrom, kto], one default real smpd_crop, then nstates volumes of
    !! box_crop**3 default reals, state 1 first. THROW_HARD for a volume not staged or a header
    !! the reader would reject (cart_refvols_header_compatible).
    subroutine write( self, fname_even, fname_odd, kfromto, smpd )
        class(cartft_calc), intent(in) :: self
        class(string),      intent(in) :: fname_even, fname_odd
        integer,            intent(in) :: kfromto(2)
        real,               intent(in) :: smpd
        character(len=32) :: field
        integer :: header(CART_REFVOLS_NHEADER), ihalf, state, funit, io_stat
        if( .not. self%exists ) THROW_HARD('cartft_calc reference volumes written before new')
        header = [CART_REFVOLS_FORMAT_VERSION, self%box, self%nstates, kfromto(1), kfromto(2)]
        if( .not. self%cart_refvols_header_compatible(header, smpd, header, smpd, smpd, field) )then
            write(logfhandle,*) 'cart_refvols header field: ', trim(field)
            THROW_HARD('cartft_calc refuses to write an incompatible reference-volume header')
        endif
        do ihalf = 1, 2
            do state = 1, self%nstates
                if( .not. allocated(self%refs(ihalf,state)%vol) ) THROW_HARD('cartft_calc reference volume not staged')
            end do
            if( ihalf == 1 )then
                call fopen(funit, fname_even, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
                call fileiochk('cartft_calc write: '//fname_even%to_char(), io_stat)
            else
                call fopen(funit, fname_odd, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
                call fileiochk('cartft_calc write: '//fname_odd%to_char(), io_stat)
            endif
            write(unit=funit, pos=1) header, smpd
            do state = 1, self%nstates
                write(unit=funit) self%refs(ihalf,state)%vol
            end do
            call fclose(funit)
        end do
    end subroutine write

    !> Build the calculator (new(nstates, box, nptcls)) and its references from the two
    !! reference-volume files: both headers are read and checked against the run's box,
    !! state count and sampling (cart_refvols_header_compatible), then every volume is padded
    !! and transformed once (set_ref). kfromto returns the band limit of the header, which the
    !! caller adopts. THROW_HARD for a missing file or an incompatible header, naming the file
    !! and the field.
    subroutine read( self, fname_even, fname_odd, nstates, box, nptcls, smpd, kfromto )
        class(cartft_calc), intent(inout) :: self
        class(string),      intent(in)    :: fname_even, fname_odd
        integer,            intent(in)    :: nstates, box, nptcls
        real,               intent(in)    :: smpd
        integer,            intent(out)   :: kfromto(2)
        character(len=32) :: field
        real, allocatable :: volume(:,:,:)
        integer :: header_even(CART_REFVOLS_NHEADER), header_odd(CART_REFVOLS_NHEADER)
        integer :: funit_even, funit_odd, io_stat, state
        real    :: smpd_even, smpd_odd
        call self%new(nstates, box, nptcls)
        if( .not. file_exists(fname_even) ) THROW_HARD('missing Cartesian reference volumes: '//fname_even%to_char())
        if( .not. file_exists(fname_odd)  ) THROW_HARD('missing Cartesian reference volumes: '//fname_odd%to_char())
        call fopen(funit_even, fname_even, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('cartft_calc read: '//fname_even%to_char(), io_stat)
        call fopen(funit_odd, fname_odd, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('cartft_calc read: '//fname_odd%to_char(), io_stat)
        read(unit=funit_even, pos=1, iostat=io_stat) header_even, smpd_even
        call fileiochk('cartft_calc read header: '//fname_even%to_char(), io_stat)
        read(unit=funit_odd, pos=1, iostat=io_stat) header_odd, smpd_odd
        call fileiochk('cartft_calc read header: '//fname_odd%to_char(), io_stat)
        if( .not. self%cart_refvols_header_compatible(header_even, smpd_even, header_odd, smpd_odd, smpd, field) )then
            write(logfhandle,*) 'incompatible field: ', trim(field), '; files: ', fname_even%to_char(), ' ', fname_odd%to_char()
            write(logfhandle,*) 'even header: ', header_even, smpd_even
            write(logfhandle,*) 'odd header:  ', header_odd, smpd_odd
            write(logfhandle,*) 'expected box, nstates, smpd: ', box, nstates, smpd
            THROW_HARD('incompatible Cartesian reference-volume header')
        endif
        allocate(volume(box,box,box))
        do state = 1, nstates
            read(unit=funit_even, iostat=io_stat) volume
            call fileiochk('cartft_calc read volume: '//fname_even%to_char(), io_stat)
            call self%set_ref(state, .true., volume)
            read(unit=funit_odd, iostat=io_stat) volume
            call fileiochk('cartft_calc read volume: '//fname_odd%to_char(), io_stat)
            call self%set_ref(state, .false., volume)
        end do
        call fclose(funit_even)
        call fclose(funit_odd)
        kfromto = header_even(4:5)
    end subroutine read

    !> The reader's compatibility rule for the two headers of the reference-volume files:
    !! format_version CART_REFVOLS_FORMAT_VERSION, box_crop and nstates those of this
    !! calculator, smpd_crop equal to smpd within single-precision tolerance,
    !! 1 <= kfrom <= kto <= fdim(box_crop) - 1, and identical even and odd headers. field names
    !! the first field that fails ('' when compatible). The reader turns false into a stop.
    logical function cart_refvols_header_compatible( self, header_even, smpd_even, header_odd, smpd_odd, smpd, field ) &
        &result( l_compatible )
        class(cartft_calc), intent(in)  :: self
        integer,            intent(in)  :: header_even(CART_REFVOLS_NHEADER), header_odd(CART_REFVOLS_NHEADER)
        real,               intent(in)  :: smpd_even, smpd_odd, smpd
        character(len=*),   intent(out) :: field
        l_compatible = .false.
        field = ''
        if( any(header_even /= header_odd) .or. smpd_even /= smpd_odd )then
            field = 'even/odd headers differ'
        else if( header_even(1) /= CART_REFVOLS_FORMAT_VERSION )then
            field = 'format_version'
        else if( header_even(2) /= self%box )then
            field = 'box_crop'
        else if( header_even(3) /= self%nstates )then
            field = 'nstates'
        else if( abs(smpd_even - smpd) > 10.*epsilon(smpd)*max(abs(smpd), 1.) )then
            field = 'smpd_crop'
        else if( header_even(4) < 1 .or. header_even(5) < header_even(4) .or. header_even(5) > self%box/2 )then
            field = 'kfrom/kto'
        else
            l_compatible = .true.
        endif
    end function cart_refvols_header_compatible

    !> True when the calculator holds the references of every state and half (a pass has read them).
    pure logical function has_refs( self )
        class(cartft_calc), intent(in) :: self
        integer :: state
        has_refs = self%exists
        if( .not. has_refs ) return
        do state = 1, self%nstates
            has_refs = has_refs .and. self%ref_exists(state, .true.) .and. self%ref_exists(state, .false.)
        end do
    end function has_refs

    !> Bytes held by the reference transforms of every state and half.
    pure integer(int64) function get_ref_bytes( self )
        class(cartft_calc), intent(in) :: self
        integer :: ihalf, state
        get_ref_bytes = 0_int64
        if( .not. self%exists ) return
        do state = 1, self%nstates
            do ihalf = 1, 2
                if( .not. allocated(self%refs(ihalf,state)%cmat) ) cycle
                get_ref_bytes = get_ref_bytes + size(self%refs(ihalf,state)%cmat, kind=int64) * &
                    &int(storage_size(self%refs(ihalf,state)%cmat)/8, int64)
            end do
        end do
    end function get_ref_bytes

    ! PARTICLES

    !> Prepare particle slot iptcl for OBJFUN_CC from its full-disk observation X (prepared as
    !! the polar particle up to the mask, times get_ptcl_taper): X and the shift-free transfer
    !! abs(C) (1 without CTF) over the shells kfromto intersected with the Cartesian Nyquist
    !! limit. No sigma2 is taken, read or held (C5). The objective applies the candidate shift
    !! itself, so including a shift here would apply it twice. An empty range or an observation
    !! without power over it leaves the slot invalid (ptcl_is_valid false): a per-particle
    !! failure, not a stop. THROW_HARD for iptcl outside 1..nptcls or an unsupported CTF flag.
    subroutine set_ptcl_cc( self, iptcl, observed, ctfparms, kfromto )
        class(cartft_calc), intent(inout) :: self
        integer,            intent(in)    :: iptcl
        complex,            intent(in)    :: observed(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2))
        type(ctfparams),    intent(in)    :: ctfparms !< at the cropped sampling
        integer,            intent(in)    :: kfromto(2)
        integer :: lower_shell, upper_shell
        if( .not. self%exists ) THROW_HARD('cartft_calc particle set before new')
        if( iptcl < 1 .or. iptcl > self%nptcls ) THROW_HARD('cartft_calc particle slot out of range')
        lower_shell = max(0, kfromto(1))
        upper_shell = min(kfromto(2), self%box/2)
        call self%fill_ptcl(iptcl, OBJFUN_CC, observed, ctfparms, lower_shell, upper_shell)
    end subroutine set_ptcl_cc

    !> The two global spectra of OBJFUN_EUCLID (shell kfrom.., half), held once for every slot
    !! with their inverse square roots and the shells where they are usable (finite, positive)
    subroutine set_sigma2( self, kfrom, sigma2 )
        class(cartft_calc), intent(inout) :: self
        integer,            intent(in)    :: kfrom
        real,               intent(in)    :: sigma2(kfrom:,:)
        integer :: kto, k, h
        if( size(sigma2,2) /= 2 ) THROW_HARD('cartft_calc sigma2 must hold one spectrum per even/odd half')
        kto = ubound(sigma2,1)
        if( allocated(self%sigma2)       ) deallocate(self%sigma2)
        if( allocated(self%rsqrt_sigma2) ) deallocate(self%rsqrt_sigma2)
        if( allocated(self%sigma2_ok)    ) deallocate(self%sigma2_ok)
        allocate(self%sigma2(kfrom:kto,2), source=sigma2)
        allocate(self%rsqrt_sigma2(kfrom:kto,2), source=0.)
        allocate(self%sigma2_ok(kfrom:kto,2))
        do h = 1, 2
            do k = kfrom, kto
                self%sigma2_ok(k,h) = ieee_is_finite(self%sigma2(k,h)) .and. self%sigma2(k,h) > 0.
                if( self%sigma2_ok(k,h) ) self%rsqrt_sigma2(k,h) = 1./sqrt(self%sigma2(k,h))
            end do
        end do
    end subroutine set_sigma2

    !> Prepare particle slot iptcl for OBJFUN_EUCLID with the spectrum of its half (set_sigma2):
    !! X/sqrt(sigma2) and abs(C)/sqrt(sigma2) of the pixel's shell over kfromto intersected with
    !! the sigma2 shells and the Cartesian Nyquist limit, and the whitened particle power that
    !! normalizes the loss. An empty range, a non-finite or non-positive sigma2 in it or an
    !! observation without power leaves the slot invalid. THROW_HARD for iptcl outside
    !! 1..nptcls, a half other than 1 or 2, no spectra or an unsupported CTF flag.
    subroutine set_ptcl_euclid( self, iptcl, observed, ctfparms, half, kfromto )
        class(cartft_calc), intent(inout) :: self
        integer,            intent(in)    :: iptcl
        complex,            intent(in)    :: observed(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2))
        type(ctfparams),    intent(in)    :: ctfparms !< at the cropped sampling
        integer,            intent(in)    :: half
        integer,            intent(in)    :: kfromto(2)
        integer :: lower_shell, upper_shell
        if( .not. self%exists ) THROW_HARD('cartft_calc particle set before new')
        if( iptcl < 1 .or. iptcl > self%nptcls ) THROW_HARD('cartft_calc particle slot out of range')
        if( .not. allocated(self%sigma2) ) THROW_HARD('cartft_calc euclid particle set before its sigma2')
        if( half < 1 .or. half > 2 ) THROW_HARD('cartft_calc euclid particle without an even/odd half')
        lower_shell = max(lbound(self%sigma2,1), kfromto(1))
        upper_shell = min(kfromto(2), ubound(self%sigma2,1), self%box/2)
        call self%fill_ptcl(iptcl, OBJFUN_EUCLID, observed, ctfparms, lower_shell, upper_shell, half)
    end subroutine set_ptcl_euclid

    !> True when particle slot iptcl holds a valid preparation.
    pure logical function ptcl_is_valid( self, iptcl )
        class(cartft_calc), intent(in) :: self
        integer,            intent(in) :: iptcl
        ptcl_is_valid = .false.
        if( .not. self%exists ) return
        if( iptcl < 1 .or. iptcl > self%nptcls ) return
        ptcl_is_valid = self%ptcls(iptcl)%valid
    end function ptcl_is_valid

    !> Shell range of particle slot iptcl: the requested range intersected with the
    !! available sigma2 shells and the Cartesian Nyquist limit.
    pure function get_ptcl_kfromto( self, iptcl ) result( kfromto )
        class(cartft_calc), intent(in) :: self
        integer,            intent(in) :: iptcl
        integer :: kfromto(2)
        kfromto = 0
        if( .not. self%exists ) return
        if( iptcl < 1 .or. iptcl > self%nptcls ) return
        kfromto = self%ptcls(iptcl)%kfromto
    end function get_ptcl_kfromto

    !> Pixels of slot iptcl's shells in the full disk, the work unit of one evaluation (0 when unprepared).
    pure integer function get_ptcl_npix( self, iptcl )
        class(cartft_calc), intent(in) :: self
        integer,            intent(in) :: iptcl
        get_ptcl_npix = 0
        if( .not. self%exists ) return
        if( iptcl < 1 .or. iptcl > self%nptcls ) return
        get_ptcl_npix = self%ptcls(iptcl)%npix
    end function get_ptcl_npix

    !> Objective the particle slot iptcl was prepared for (0 when never prepared).
    pure integer function get_ptcl_objfun( self, iptcl )
        class(cartft_calc), intent(in) :: self
        integer,            intent(in) :: iptcl
        get_ptcl_objfun = 0
        if( .not. self%exists ) return
        if( iptcl < 1 .or. iptcl > self%nptcls ) return
        get_ptcl_objfun = self%ptcls(iptcl)%objfun
    end function get_ptcl_objfun

    !> Separable real-space taper of the gather stencil along one axis of the box, centre at
    !! box/2+1 (O5 (a)). The polar particle is sampled from its OSMPL_PAD_FAC-padded transform
    !! by the normalized KB stencil (polarize_oversamp); at an on-grid sample that equals the
    !! native transform of the masked image times taper(i)*taper(j), the stencil's response
    !! at period padf*box cropped to the box. The reference carries the same taper through
    !! the gather, so the observation is multiplied by it before its transform.
    function get_ptcl_taper( self ) result( taper )
        class(cartft_calc), intent(in) :: self
        real, allocatable :: taper(:), padded_env(:)
        integer :: offset
        if( .not. self%exists ) THROW_HARD('cartft_calc taper requested before new')
        call kb_stencil_envelope_1d(self%kbwin, self%boxpd, padded_env)
        offset = (self%boxpd - self%box)/2
        taper  = padded_env(offset+1:offset+self%box)
        deallocate(padded_env)
    end function get_ptcl_taper

    !> (Re)allocate nptcls empty particle slots, the slots of one particle batch, keeping the
    !! references. THROW_HARD before new or for nptcls < 1.
    subroutine new_ptcls( self, nptcls )
        class(cartft_calc), intent(inout) :: self
        integer,            intent(in)    :: nptcls
        if( .not. self%exists ) THROW_HARD('cartft_calc particle slots sized before new')
        if( nptcls < 1 ) THROW_HARD('cartft_calc requires at least one particle slot')
        if( allocated(self%ptcls) ) deallocate(self%ptcls)
        if( allocated(self%slots) ) deallocate(self%slots)
        allocate(self%ptcls(nptcls))
        self%nptcls = nptcls
    end subroutine new_ptcls

    !> Record which particle each slot of the current batch holds: slot i holds particle
    !! pinds(i), as the polar calculator maps particle indices (its pinds). THROW_HARD for more
    !! particles than slots, a non-positive or a repeated index.
    subroutine set_ptcl_inds( self, pinds )
        class(cartft_calc), intent(inout) :: self
        integer,            intent(in)    :: pinds(:)
        integer :: i
        if( size(pinds) > self%nptcls ) THROW_HARD('cartft_calc batch has more particles than slots')
        if( size(pinds) < 1 ) THROW_HARD('cartft_calc batch has no particle')
        if( minval(pinds) < 1 ) THROW_HARD('cartft_calc particle index must be positive')
        if( allocated(self%slots) ) deallocate(self%slots)
        allocate(self%slots(minval(pinds):maxval(pinds)), source=0)
        do i = 1, size(pinds)
            if( self%slots(pinds(i)) /= 0 ) THROW_HARD('cartft_calc batch repeats a particle')
            self%slots(pinds(i)) = i
        end do
    end subroutine set_ptcl_inds

    !> The slot of particle iptcl in the current batch, 0 when the batch does not hold it.
    pure integer function get_ptcl_slot( self, iptcl )
        class(cartft_calc), intent(in) :: self
        integer,            intent(in) :: iptcl
        get_ptcl_slot = 0
        if( .not. allocated(self%slots) ) return
        if( iptcl < lbound(self%slots,1) .or. iptcl > ubound(self%slots,1) ) return
        get_ptcl_slot = self%slots(iptcl)
    end function get_ptcl_slot

    !> The number of particle slots.
    pure integer function get_nptcls( self )
        class(cartft_calc), intent(in) :: self
        get_nptcls = self%nptcls
    end function get_nptcls

    !> The box of particles and references.
    pure integer function get_box( self )
        class(cartft_calc), intent(in) :: self
        get_box = self%box
    end function get_box

    ! EVALUATION (calculator intent(in); safe to call concurrently on different particles)

    !> Unweighted full-disk prediction S(shift) G(rotmat) V of one reference, the forward
    !! model before CTF and whitening, over the pixels of the objectives' membership for shells
    !! 0..box/2 (the Nyquist disk). rotmat is the particle rotation, shift in pixels of box.
    subroutine predict( self, state, iseven, rotmat, shift, prediction )
        class(cartft_calc), intent(in)  :: self
        integer,            intent(in)  :: state
        logical,            intent(in)  :: iseven
        real(dp),           intent(in)  :: rotmat(3,3), shift(2)
        complex,            intent(out) :: prediction(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2))
        complex  :: val, dval(3), phase
        real(sp) :: loc(3)
        real(dp) :: l_arg
        integer  :: h, k, ihalf
        logical  :: inside
        if( .not. self%ref_exists(state, iseven) ) THROW_HARD('cartft_calc prediction from an absent reference')
        ihalf = self%ref_index(iseven)
        prediction = cmplx(0.,0.)
        associate( cmat => self%refs(ihalf,state)%cmat )
            do k = self%lims2(2,1), self%lims2(2,2)
                do h = self%lims2(1,1), self%lims2(1,2)
                    if( self%shells(h,k) < 0 ) cycle
                    loc = real(self%padf, sp)*real(matmul(real([h, k, 0], dp), rotmat), sp)
                    call self%gather(cmat, loc, val, dval, inside)
                    if( .not. inside ) THROW_HARD('cartft_calc gather location outside the wrap table')
                    l_arg = 2._dp*real(PI, dp)*(real(h, dp)*shift(1) + real(k, dp)*shift(2))/real(self%box, dp)
                    phase = cmplx(cos(l_arg), sin(l_arg), kind=sp)
                    prediction(h,k) = phase*val
                end do
            end do
        end associate
    end subroutine predict

    !> Objective of particle slot iptcl against one reference and its five derivatives (three
    !! right-increment rotation components in radians, two shifts in pixels), for the
    !! objective the slot was prepared for: 1 - cc (OBJFUN_CC) or the normalized loss L
    !! (OBJFUN_EUCLID). THROW_HARD for an invalid slot or an absent reference.
    subroutine objective_gradient( self, state, iseven, iptcl, rotmat, shift, objective, gradient )
        class(cartft_calc), intent(in)  :: self
        integer,            intent(in)  :: state, iptcl
        logical,            intent(in)  :: iseven
        real(dp),           intent(in)  :: rotmat(3,3), shift(2)
        real(dp),           intent(out) :: objective, gradient(5)
        real(dp) :: hessian(5,5)
        call self%pose_normal_terms(state, iseven, iptcl, rotmat, shift, objective, gradient, hessian)
    end subroutine objective_gradient

    !> Objective, two-vector shift gradient and 2x2 Gauss-Newton block at a fixed rotation,
    !! for the optimizer's shift stage, over the half plane of the slot's band (every summand is
    !! Hermitian between Friedel mates). NaN objective when a correlation is undefined.
    subroutine shift_normal_terms( self, state, iseven, iptcl, rotmat, shift, objective, gradient, hessian )
        class(cartft_calc), intent(in)  :: self
        integer,            intent(in)  :: state, iptcl
        logical,            intent(in)  :: iseven
        real(dp),           intent(in)  :: rotmat(3,3), shift(2)
        real(dp),           intent(out) :: objective, gradient(2), hessian(2,2)
        complex(dp) :: eh(0:self%ptcls(iptcl)%kfromto(2)), ek(-self%ptcls(iptcl)%kfromto(2):self%ptcls(iptcl)%kfromto(2))
        complex     :: val, dval(3)
        complex(dp) :: model, particle, residual, xm
        real(sp)    :: loc(3)
        real(dp)    :: row1(3), row2(3), fh, fk, fw, m2, q
        real(dp)    :: prediction_power, cross_real, weighted_power, weighted_cross, model_gradient(2)
        integer     :: h, k, kto, ihalf, shell
        logical     :: inside, l_euclid
        call self%check_eval_args(state, iseven, iptcl)
        ihalf = self%ref_index(iseven)
        associate( p => self%ptcls(iptcl), cmat => self%refs(ihalf,state)%cmat )
            kto      = p%kfromto(2)
            l_euclid = p%objfun == OBJFUN_EUCLID
            call self%shift_phase_tables(shift, eh, ek)
            row1 = real(self%padf, dp)*rotmat(1,:)
            row2 = real(self%padf, dp)*rotmat(2,:)
            objective        = 0._dp
            gradient         = 0._dp
            hessian          = 0._dp
            prediction_power = 0._dp
            cross_real       = 0._dp
            weighted_power   = 0._dp
            weighted_cross   = 0._dp
            model_gradient   = 0._dp
            do k = -kto, kto
                fk = 2._dp*real(PI, dp)*real(k, dp)/real(self%box, dp)
                ! the half plane: h > 0, and h = 0 for k > 0; the DC pixel is added once below
                do h = merge(0, 1, k > 0), kto
                    shell = self%shells(h,k)
                    if( shell < p%kfromto(1) .or. shell > kto ) cycle
                    loc = real(real(h, dp)*row1 + real(k, dp)*row2, sp)
                    call self%gather(cmat, loc, val, dval, inside)
                    if( .not. inside ) THROW_HARD('cartft_calc gather location outside the wrap table')
                    model    = eh(h)*ek(k)*real(p%transfer(h,k), dp)*cmplx(val, kind=dp)
                    particle = cmplx(p%observed(h,k), kind=dp)
                    fh = 2._dp*real(PI, dp)*real(h, dp)/real(self%box, dp)
                    fw = fh*fh + fk*fk
                    m2 = real(conjg(model)*model, dp)
                    xm = conjg(particle)*model
                    q  = aimag(xm)
                    weighted_power = weighted_power + fw*m2
                    weighted_cross = weighted_cross + fw*real(xm, dp)
                    ! a shift Jacobian is i f model: its gradient terms are -+ f Im(conj(X) M)
                    if( l_euclid )then
                        residual  = model - particle
                        objective = objective + real(conjg(residual)*residual, dp)
                        gradient  = gradient + [fh, fk]*q
                    else
                        prediction_power = prediction_power + m2
                        cross_real       = cross_real + real(xm, dp)
                        gradient         = gradient - [fh, fk]*q
                    endif
                    hessian(1,1) = hessian(1,1) + fh*fh*m2
                    hessian(1,2) = hessian(1,2) + fh*fk*m2
                    hessian(2,2) = hessian(2,2) + fk*fk*m2
                end do
            end do
            ! each half-plane pixel stands for its Friedel mate as well
            objective        = 2._dp*objective
            gradient         = 2._dp*gradient
            hessian          = 2._dp*hessian
            prediction_power = 2._dp*prediction_power
            cross_real       = 2._dp*cross_real
            weighted_power   = 2._dp*weighted_power
            weighted_cross   = 2._dp*weighted_cross
            if( p%kfromto(1) == 0 ) call add_dc_pixel(self, p, cmat, l_euclid, objective, prediction_power, cross_real)
            hessian(2,1) = hessian(1,2)
            if( l_euclid )then
                call finalize_euclid_normal_terms(p%power, weighted_power, weighted_cross, objective, gradient, hessian)
            else
                call finalize_cc_normal_terms(p%power, prediction_power, cross_real, weighted_power, &
                    &weighted_cross, model_gradient, objective, gradient, hessian)
            endif
        end associate
    end subroutine shift_normal_terms

    !> Objective, five-vector gradient and 5x5 Gauss-Newton block (rotations 1:3, shifts 4:5), for
    !! the optimizer's joint stage, over the half plane of the slot's band (every summand is
    !! Hermitian between Friedel mates); with shell_residual the whitened |r|^2 per shell of a
    !! euclid slot (zero for cc), for sigma_from_shell_residual. NaN objective when a correlation is undefined.
    subroutine pose_normal_terms( self, state, iseven, iptcl, rotmat, shift, objective, gradient, hessian, shell_residual )
        class(cartft_calc), intent(in)  :: self
        integer,            intent(in)  :: state, iptcl
        logical,            intent(in)  :: iseven
        real(dp),           intent(in)  :: rotmat(3,3), shift(2)
        real(dp),           intent(out) :: objective, gradient(5), hessian(5,5)
        real(dp), optional, intent(out) :: shell_residual(self%ptcls(iptcl)%kfromto(1):self%ptcls(iptcl)%kfromto(2))
        complex(dp) :: eh(0:self%ptcls(iptcl)%kfromto(2)), ek(-self%ptcls(iptcl)%kfromto(2):self%ptcls(iptcl)%kfromto(2))
        complex     :: val, dval(3)
        complex(dp) :: wp, model, particle, residual, xm, jacobian(3)
        real(sp)    :: loc(3)
        real(dp)    :: row1(3), row2(3), fh, fk, fw, m2, q, u(3), l(3), r2, dc_residual
        real(dp)    :: prediction_power, cross_real, weighted_power, weighted_cross, model_gradient(5)
        integer     :: a, b, h, k, kto, ihalf, shell
        logical     :: inside, l_euclid, l_shells
        call self%check_eval_args(state, iseven, iptcl)
        ihalf = self%ref_index(iseven)
        associate( p => self%ptcls(iptcl), cmat => self%refs(ihalf,state)%cmat )
            kto      = p%kfromto(2)
            l_euclid = p%objfun == OBJFUN_EUCLID
            l_shells = present(shell_residual) .and. l_euclid
            if( present(shell_residual) ) shell_residual = 0._dp
            call self%shift_phase_tables(shift, eh, ek)
            ! the padded rotation rows: pixel (h,k) samples h*row1 + k*row2
            row1 = real(self%padf, dp)*rotmat(1,:)
            row2 = real(self%padf, dp)*rotmat(2,:)
            objective        = 0._dp
            gradient         = 0._dp
            hessian          = 0._dp
            prediction_power = 0._dp
            cross_real       = 0._dp
            weighted_power   = 0._dp
            weighted_cross   = 0._dp
            model_gradient   = 0._dp
            do k = -kto, kto
                fk = 2._dp*real(PI, dp)*real(k, dp)/real(self%box, dp)
                ! the half plane: h > 0, and h = 0 for k > 0; the DC pixel is added once below
                do h = merge(0, 1, k > 0), kto
                    shell = self%shells(h,k)
                    if( shell < p%kfromto(1) .or. shell > kto ) cycle
                    loc = real(real(h, dp)*row1 + real(k, dp)*row2, sp)
                    call self%gather(cmat, loc, val, dval, inside)
                    if( .not. inside ) THROW_HARD('cartft_calc gather location outside the wrap table')
                    wp       = eh(h)*ek(k)*real(p%transfer(h,k), dp)
                    model    = wp*cmplx(val, kind=dp)
                    particle = cmplx(p%observed(h,k), kind=dp)
                    ! the rotation Jacobians, the derivative along loc x e_a: wp (dV x loc)
                    l = real(loc, dp)
                    jacobian(1) = wp*(cmplx(dval(2), kind=dp)*l(3) - cmplx(dval(3), kind=dp)*l(2))
                    jacobian(2) = wp*(cmplx(dval(3), kind=dp)*l(1) - cmplx(dval(1), kind=dp)*l(3))
                    jacobian(3) = wp*(cmplx(dval(1), kind=dp)*l(2) - cmplx(dval(2), kind=dp)*l(1))
                    fh = 2._dp*real(PI, dp)*real(h, dp)/real(self%box, dp)
                    fw = fh*fh + fk*fk
                    m2 = real(conjg(model)*model, dp)
                    xm = conjg(particle)*model
                    q  = aimag(xm)
                    do a = 1, 3
                        u(a) = aimag(conjg(jacobian(a))*model)
                    end do
                    weighted_power = weighted_power + fw*m2
                    weighted_cross = weighted_cross + fw*real(xm, dp)
                    ! a shift Jacobian is i f model: its gradient terms are -+ f Im(conj(X) M), and
                    ! model_gradient has no shift component
                    if( l_euclid )then
                        residual  = model - particle
                        r2        = real(conjg(residual)*residual, dp)
                        objective = objective + r2
                        if( l_shells ) shell_residual(shell) = shell_residual(shell) + r2
                        do a = 1, 3
                            gradient(a) = gradient(a) + real(conjg(jacobian(a))*residual, dp)
                        end do
                        gradient(4) = gradient(4) + fh*q
                        gradient(5) = gradient(5) + fk*q
                    else
                        prediction_power = prediction_power + m2
                        cross_real       = cross_real + real(xm, dp)
                        do a = 1, 3
                            gradient(a)       = gradient(a)       + real(conjg(particle)*jacobian(a), dp)
                            model_gradient(a) = model_gradient(a) + real(conjg(model)*jacobian(a), dp)
                        end do
                        gradient(4) = gradient(4) - fh*q
                        gradient(5) = gradient(5) - fk*q
                    endif
                    ! the 15 distinct entries; rotation-shift entries are -f Im(conj(J_a) M)
                    do a = 1, 3
                        do b = a, 3
                            hessian(a,b) = hessian(a,b) + real(conjg(jacobian(a))*jacobian(b), dp)
                        end do
                        hessian(a,4) = hessian(a,4) - fh*u(a)
                        hessian(a,5) = hessian(a,5) - fk*u(a)
                    end do
                    hessian(4,4) = hessian(4,4) + fh*fh*m2
                    hessian(4,5) = hessian(4,5) + fh*fk*m2
                    hessian(5,5) = hessian(5,5) + fk*fk*m2
                end do
            end do
            ! each half-plane pixel stands for its Friedel mate as well
            objective        = 2._dp*objective
            gradient         = 2._dp*gradient
            hessian          = 2._dp*hessian
            prediction_power = 2._dp*prediction_power
            cross_real       = 2._dp*cross_real
            weighted_power   = 2._dp*weighted_power
            weighted_cross   = 2._dp*weighted_cross
            model_gradient   = 2._dp*model_gradient
            if( l_shells ) shell_residual = 2._dp*shell_residual
            if( p%kfromto(1) == 0 )then
                call add_dc_pixel(self, p, cmat, l_euclid, objective, prediction_power, cross_real, dc_residual)
                if( l_shells ) shell_residual(0) = shell_residual(0) + dc_residual
            endif
            do a = 2, 5
                do b = 1, a - 1
                    hessian(a,b) = hessian(b,a)
                end do
            end do
            if( l_euclid )then
                call finalize_euclid_normal_terms(p%power, weighted_power, weighted_cross, objective, gradient, hessian)
            else
                call finalize_cc_normal_terms(p%power, prediction_power, cross_real, weighted_power, &
                    &weighted_cross, model_gradient, objective, gradient, hessian)
            endif
        end associate
    end subroutine pose_normal_terms

    !> shift_normal_terms as it was before the half-plane kernel (cont_lm_acceleration.md, stage 6):
    !! the full disk, kept as the reference the tests compare the kernel with.
    subroutine shift_normal_terms_full_disk( self, state, iseven, iptcl, rotmat, shift, objective, gradient, hessian )
        class(cartft_calc), intent(in)  :: self
        integer,            intent(in)  :: state, iptcl
        logical,            intent(in)  :: iseven
        real(dp),           intent(in)  :: rotmat(3,3), shift(2)
        real(dp),           intent(out) :: objective, gradient(2), hessian(2,2)
        complex     :: val, dval(3), phase
        complex(dp) :: model, particle, residual, jacobian(2)
        real(sp)    :: loc(3)
        real(dp)    :: l_arg, frequency(2), particle_power, prediction_power, cross_real
        real(dp)    :: frequency_weight, weighted_power, weighted_cross
        real(dp)    :: model_gradient(2)
        integer     :: axis, h, jaxis, k, ihalf, objfun, shell
        logical     :: inside
        call self%check_eval_args(state, iseven, iptcl)
        ihalf = self%ref_index(iseven)
        associate( p => self%ptcls(iptcl), cmat => self%refs(ihalf,state)%cmat )
            objfun           = p%objfun
            objective        = 0._dp
            gradient         = 0._dp
            hessian          = 0._dp
            particle_power   = 0._dp
            prediction_power = 0._dp
            cross_real       = 0._dp
            weighted_power   = 0._dp
            weighted_cross   = 0._dp
            model_gradient   = 0._dp
            do k = self%lims2(2,1), self%lims2(2,2)
                do h = self%lims2(1,1), self%lims2(1,2)
                    shell = self%shells(h,k)
                    if( shell < p%kfromto(1) .or. shell > p%kfromto(2) ) cycle
                    loc = real(self%padf, sp)*real(matmul(real([h, k, 0], dp), rotmat), sp)
                    call self%gather(cmat, loc, val, dval, inside)
                    if( .not. inside ) THROW_HARD('cartft_calc gather location outside the wrap table')
                    l_arg = 2._dp*real(PI, dp)*(real(h, dp)*shift(1) + real(k, dp)*shift(2))/real(self%box, dp)
                    phase = cmplx(cos(l_arg), sin(l_arg), kind=sp)
                    model = cmplx(phase*val, kind=dp)
                    model = model*cmplx(p%transfer(h,k), kind=dp)
                    particle  = cmplx(p%observed(h,k), kind=dp)
                    residual  = model - particle
                    frequency = 2._dp*real(PI, dp)*real([h, k], dp)/real(self%box, dp)
                    frequency_weight = dot_product(frequency, frequency)
                    jacobian  = cmplx(0._dp, frequency, kind=dp)*model
                    select case(objfun)
                        case(OBJFUN_EUCLID)
                            objective = objective + real(conjg(residual)*residual, dp)
                            weighted_power = weighted_power + frequency_weight*real(conjg(model)*model, dp)
                            weighted_cross = weighted_cross + frequency_weight*real(conjg(particle)*model, dp)
                            do axis = 1, 2
                                gradient(axis) = gradient(axis) + real(conjg(jacobian(axis))*residual, dp)
                                do jaxis = 1, 2
                                    hessian(axis,jaxis) = hessian(axis,jaxis) + real(conjg(jacobian(axis))*jacobian(jaxis), dp)
                                end do
                            end do
                        case(OBJFUN_CC)
                            particle_power   = particle_power   + real(conjg(particle)*particle, dp)
                            prediction_power = prediction_power + real(conjg(model)*model, dp)
                            cross_real       = cross_real       + real(conjg(particle)*model, dp)
                            weighted_power   = weighted_power   + frequency_weight*real(conjg(model)*model, dp)
                            weighted_cross   = weighted_cross   + frequency_weight*real(conjg(particle)*model, dp)
                            do axis = 1, 2
                                gradient(axis)       = gradient(axis)       + real(conjg(particle)*jacobian(axis), dp)
                                model_gradient(axis) = model_gradient(axis) + real(conjg(model)*jacobian(axis), dp)
                                do jaxis = 1, 2
                                    hessian(axis,jaxis) = hessian(axis,jaxis) + real(conjg(jacobian(axis))*jacobian(jaxis), dp)
                                end do
                            end do
                    end select
                end do
            end do
            select case(objfun)
                case(OBJFUN_EUCLID)
                    call finalize_euclid_normal_terms(p%power, weighted_power, weighted_cross, objective, gradient, hessian)
                case(OBJFUN_CC)
                    call finalize_cc_normal_terms(particle_power, prediction_power, cross_real, weighted_power, &
                        &weighted_cross, model_gradient, objective, gradient, hessian)
            end select
        end associate
    end subroutine shift_normal_terms_full_disk

    !> pose_normal_terms as it was before the half-plane kernel (cont_lm_acceleration.md, stage 6):
    !! the full disk, kept as the reference the tests compare the kernel with.
    subroutine pose_normal_terms_full_disk( self, state, iseven, iptcl, rotmat, shift, objective, gradient, hessian )
        class(cartft_calc), intent(in)  :: self
        integer,            intent(in)  :: state, iptcl
        logical,            intent(in)  :: iseven
        real(dp),           intent(in)  :: rotmat(3,3), shift(2)
        real(dp),           intent(out) :: objective, gradient(5), hessian(5,5)
        complex     :: val, dval(3), phase
        complex(dp) :: weighted_phase, model, particle, residual, jacobian(5)
        real(sp)    :: loc(3)
        real(dp)    :: l_arg, dloc(3,3), frequency(2)
        real(dp)    :: particle_power, prediction_power, cross_real
        real(dp)    :: frequency_weight, weighted_power, weighted_cross
        real(dp)    :: model_gradient(5)
        integer     :: axis, h, jaxis, k, ihalf, objfun, shell
        logical     :: inside
        call self%check_eval_args(state, iseven, iptcl)
        ihalf = self%ref_index(iseven)
        associate( p => self%ptcls(iptcl), cmat => self%refs(ihalf,state)%cmat )
            objfun           = p%objfun
            objective        = 0._dp
            gradient         = 0._dp
            hessian          = 0._dp
            particle_power   = 0._dp
            prediction_power = 0._dp
            cross_real       = 0._dp
            weighted_power   = 0._dp
            weighted_cross   = 0._dp
            model_gradient   = 0._dp
            do k = self%lims2(2,1), self%lims2(2,2)
                do h = self%lims2(1,1), self%lims2(1,2)
                    shell = self%shells(h,k)
                    if( shell < p%kfromto(1) .or. shell > p%kfromto(2) ) cycle
                    loc = real(self%padf, sp)*real(matmul(real([h, k, 0], dp), rotmat), sp)
                    call self%gather(cmat, loc, val, dval, inside)
                    if( .not. inside ) THROW_HARD('cartft_calc gather location outside the wrap table')
                    ! columns are loc x e1, loc x e2 and loc x e3
                    dloc(:,1) = [0._dp, real(loc(3), dp), -real(loc(2), dp)]
                    dloc(:,2) = [-real(loc(3), dp), 0._dp, real(loc(1), dp)]
                    dloc(:,3) = [real(loc(2), dp), -real(loc(1), dp), 0._dp]
                    l_arg = 2._dp*real(PI, dp)*(real(h, dp)*shift(1) + real(k, dp)*shift(2))/real(self%box, dp)
                    phase = cmplx(cos(l_arg), sin(l_arg), kind=sp)
                    weighted_phase = cmplx(phase, kind=dp)
                    weighted_phase = weighted_phase*cmplx(p%transfer(h,k), kind=dp)
                    model    = weighted_phase*cmplx(val, kind=dp)
                    particle = cmplx(p%observed(h,k), kind=dp)
                    residual = model - particle
                    do axis = 1, 3
                        jacobian(axis) = weighted_phase*sum(cmplx(dval, kind=dp)*dloc(:,axis))
                    end do
                    frequency     = 2._dp*real(PI, dp)*real([h, k], dp)/real(self%box, dp)
                    frequency_weight = dot_product(frequency, frequency)
                    jacobian(4:5) = cmplx(0._dp, frequency, kind=dp)*model
                    select case(objfun)
                        case(OBJFUN_EUCLID)
                            objective = objective + real(conjg(residual)*residual, dp)
                            weighted_power = weighted_power + frequency_weight*real(conjg(model)*model, dp)
                            weighted_cross = weighted_cross + frequency_weight*real(conjg(particle)*model, dp)
                            do axis = 1, 5
                                gradient(axis) = gradient(axis) + real(conjg(jacobian(axis))*residual, dp)
                                do jaxis = 1, 5
                                    hessian(axis,jaxis) = hessian(axis,jaxis) + real(conjg(jacobian(axis))*jacobian(jaxis), dp)
                                end do
                            end do
                        case(OBJFUN_CC)
                            particle_power   = particle_power   + real(conjg(particle)*particle, dp)
                            prediction_power = prediction_power + real(conjg(model)*model, dp)
                            cross_real       = cross_real       + real(conjg(particle)*model, dp)
                            weighted_power   = weighted_power   + frequency_weight*real(conjg(model)*model, dp)
                            weighted_cross   = weighted_cross   + frequency_weight*real(conjg(particle)*model, dp)
                            do axis = 1, 5
                                gradient(axis)       = gradient(axis)       + real(conjg(particle)*jacobian(axis), dp)
                                model_gradient(axis) = model_gradient(axis) + real(conjg(model)*jacobian(axis), dp)
                                do jaxis = 1, 5
                                    hessian(axis,jaxis) = hessian(axis,jaxis) + real(conjg(jacobian(axis))*jacobian(jaxis), dp)
                                end do
                            end do
                    end select
                end do
            end do
            select case(objfun)
                case(OBJFUN_EUCLID)
                    call finalize_euclid_normal_terms(p%power, weighted_power, weighted_cross, objective, gradient, hessian)
                case(OBJFUN_CC)
                    call finalize_cc_normal_terms(particle_power, prediction_power, cross_real, weighted_power, &
                        &weighted_cross, model_gradient, objective, gradient, hessian)
            end select
        end associate
    end subroutine pose_normal_terms_full_disk

    !> Unwhitened per-shell accounting at one terminal pose of an OBJFUN_EUCLID particle slot:
    !! residual power over two per sample (the sigma2 contribution) over the particle's shells.
    !! The optional diagnostics, used as test instruments, add reference and particle power per
    !! shell and v = the normalized loss L (-1 when undefined); production omits them and pays
    !! for neither their storage nor their arithmetic. An invalid slot returns v = -1 and
    !! unallocated arrays. THROW_HARD for an OBJFUN_CC slot: a correlation pass writes no
    !! sigma2 (C5).
    subroutine sigma_contribution( self, state, iseven, iptcl, rotmat, shift, sigma_contrib, ref_pow, ptcl_pow, v )
        class(cartft_calc),          intent(in)  :: self
        integer,                     intent(in)  :: state, iptcl
        logical,                     intent(in)  :: iseven
        real(dp),                    intent(in)  :: rotmat(3,3), shift(2)
        real, allocatable,           intent(out) :: sigma_contrib(:)
        real, allocatable, optional, intent(out) :: ref_pow(:), ptcl_pow(:)
        real,              optional, intent(out) :: v
        complex(dp), allocatable :: eh(:), ek(:)
        complex     :: val, dval(3)
        complex(dp) :: model, raw_observed, residual
        real(dp), allocatable :: sigma_sum(:), ref_sum(:), ptcl_sum(:), root_sigma2(:)
        real(dp)    :: row1(3), row2(3), root_sigma, vnum, vden, weight
        real(sp)    :: loc(3)
        integer, allocatable :: counts(:)
        integer     :: h, k, shell, lower_shell, upper_shell, ihalf, multiplicity
        logical     :: inside, l_diag
        if( .not. self%ref_exists(state, iseven) ) THROW_HARD('cartft_calc sigma contribution from an absent reference')
        if( iptcl < 1 .or. iptcl > self%nptcls )   THROW_HARD('cartft_calc particle slot out of range')
        l_diag = present(ref_pow) .or. present(ptcl_pow) .or. present(v)
        if( present(v) ) v = -1.
        if( .not. self%ptcls(iptcl)%valid ) return
        if( self%ptcls(iptcl)%objfun /= OBJFUN_EUCLID ) THROW_HARD('cartft_calc sigma contribution requires a euclid particle slot')
        ihalf = self%ref_index(iseven)
        associate( p => self%ptcls(iptcl), cmat => self%refs(ihalf,state)%cmat )
            lower_shell = p%kfromto(1)
            upper_shell = p%kfromto(2)
            ! the per-shell sums are tiny; only the diagnostic arithmetic is conditional
            allocate(sigma_sum(lower_shell:upper_shell), source=0._dp)
            allocate(ref_sum(lower_shell:upper_shell),   source=0._dp)
            allocate(ptcl_sum(lower_shell:upper_shell),  source=0._dp)
            allocate(counts(lower_shell:upper_shell),    source=0)
            ! the square root of the half's spectrum per shell, once per call
            allocate(root_sigma2(lower_shell:upper_shell))
            root_sigma2 = sqrt(real(self%sigma2(lower_shell:upper_shell,p%half), dp))
            vnum = 0._dp
            vden = 0._dp
            allocate(eh(0:upper_shell), ek(-upper_shell:upper_shell))
            call self%shift_phase_tables(shift, eh, ek)
            row1 = real(self%padf, dp)*rotmat(1,:)
            row2 = real(self%padf, dp)*rotmat(2,:)
            ! the half plane, h > 0 and h = 0 for k >= 0: a pixel counts twice for its Friedel mate,
            ! the DC pixel once
            do k = -upper_shell, upper_shell
                do h = merge(0, 1, k >= 0), upper_shell
                    shell = self%shells(h,k)
                    if( shell < lower_shell .or. shell > upper_shell ) cycle
                    multiplicity = merge(1, 2, h == 0 .and. k == 0)
                    weight       = real(multiplicity, dp)
                    loc   = real(real(h, dp)*row1 + real(k, dp)*row2, sp)
                    call self%gather(cmat, loc, val, dval, inside)
                    if( .not. inside ) THROW_HARD('cartft_calc gather location outside the wrap table')
                    root_sigma   = root_sigma2(shell)
                    model        = eh(h)*ek(k)*real(p%transfer(h,k), dp)*cmplx(val, kind=dp)*root_sigma
                    raw_observed = cmplx(p%observed(h,k), kind=dp)*root_sigma
                    residual     = raw_observed - model
                    sigma_sum(shell) = sigma_sum(shell) + weight*real(conjg(residual)*residual, dp)
                    counts(shell)    = counts(shell) + multiplicity
                    if( l_diag )then
                        ref_sum(shell)  = ref_sum(shell)  + weight*real(conjg(model)*model, dp)
                        ptcl_sum(shell) = ptcl_sum(shell) + weight*real(conjg(raw_observed)*raw_observed, dp)
                        vnum = vnum + weight*real(conjg(residual)*residual, dp)/real(self%sigma2(shell,p%half), dp)
                        vden = vden + weight*real(conjg(raw_observed)*raw_observed, dp)/real(self%sigma2(shell,p%half), dp)
                    endif
                end do
            end do
        end associate
        if( any(counts == 0) ) THROW_HARD('cartft_calc sigma contribution found an empty active shell')
        allocate(sigma_contrib(lower_shell:upper_shell))
        sigma_contrib = real(sigma_sum/(2._dp*real(counts, dp)), sp)
        if( present(ref_pow) )then
            allocate(ref_pow(lower_shell:upper_shell))
            ref_pow = real(ref_sum/real(counts, dp), sp)
        endif
        if( present(ptcl_pow) )then
            allocate(ptcl_pow(lower_shell:upper_shell))
            ptcl_pow = real(ptcl_sum/real(counts, dp), sp)
        endif
        if( present(v) )then
            if( vden > 0._dp )then
                v = real(vnum/vden, sp)
            else
                v = -1.
            endif
        endif
    end subroutine sigma_contribution

    !> The sigma2 contribution of euclid particle slot iptcl from the whitened |r|^2 per shell of an
    !! evaluation (pose_normal_terms): residual power over two per sample, unwhitened by the half's
    !! spectrum, as sigma_contribution sums it at that pose without a pass of its own.
    subroutine sigma_from_shell_residual( self, iptcl, shell_residual, sigma_contrib )
        class(cartft_calc), intent(in)  :: self
        integer,            intent(in)  :: iptcl
        real(dp),           intent(in)  :: shell_residual(:)
        real, allocatable,  intent(out) :: sigma_contrib(:)
        integer :: shell, lower_shell, upper_shell
        if( iptcl < 1 .or. iptcl > self%nptcls ) THROW_HARD('cartft_calc particle slot out of range')
        if( .not. self%ptcls(iptcl)%valid ) THROW_HARD('cartft_calc sigma contribution requires a valid particle slot')
        if( self%ptcls(iptcl)%objfun /= OBJFUN_EUCLID ) THROW_HARD('cartft_calc sigma contribution requires a euclid particle slot')
        lower_shell = self%ptcls(iptcl)%kfromto(1)
        upper_shell = self%ptcls(iptcl)%kfromto(2)
        if( size(shell_residual) /= upper_shell - lower_shell + 1 ) THROW_HARD('cartft_calc shell residual does not span the slot band')
        if( any(self%shell_npix(lower_shell:upper_shell) == 0) ) THROW_HARD('cartft_calc sigma contribution found an empty active shell')
        allocate(sigma_contrib(lower_shell:upper_shell))
        do shell = lower_shell, upper_shell
            sigma_contrib(shell) = real(real(self%sigma2(shell,self%ptcls(iptcl)%half), dp)*&
                &shell_residual(shell - lower_shell + 1)/(2._dp*real(self%shell_npix(shell), dp)), sp)
        end do
    end subroutine sigma_from_shell_residual

    !> The stored score of an objective value of particle slot iptcl, as the polar branch stores
    !! it: cc = 1 - objective under OBJFUN_CC, exp(-L) under OBJFUN_EUCLID (C6). THROW_HARD for a
    !! slot that was never prepared.
    real(dp) function score( self, iptcl, objective )
        class(cartft_calc), intent(in) :: self
        integer,            intent(in) :: iptcl
        real(dp),           intent(in) :: objective
        select case(self%get_ptcl_objfun(iptcl))
            case(OBJFUN_CC)
                score = 1._dp - objective
            case(OBJFUN_EUCLID)
                score = exp(-objective)
            case DEFAULT
                score = 0._dp
                THROW_HARD('cartft_calc score of an unprepared particle slot')
        end select
    end function score

    ! PRIVATE

    !> Fill particle slot iptcl for objfun over the shells lower_shell..upper_shell; half is
    !! present exactly for OBJFUN_EUCLID (set_ptcl_euclid): its spectrum whitens observation and
    !! transfer.
    subroutine fill_ptcl( self, iptcl, objfun, observed, ctfparms, lower_shell, upper_shell, half )
        class(cartft_calc), intent(inout) :: self
        integer,            intent(in)    :: iptcl, objfun, lower_shell, upper_shell
        complex,            intent(in)    :: observed(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2))
        type(ctfparams),    intent(in)    :: ctfparms
        integer, optional,  intent(in)    :: half
        type(ctf)     :: tfun
        type(ctfvars) :: ctfvals
        real     :: cval, weight
        real(dp) :: half_power, dc_power
        integer  :: h, k, shell
        logical  :: use_ctf
        select case(ctfparms%ctfflag)
            case(CTFFLAG_NO, CTFFLAG_YES, CTFFLAG_FLIP)
            case DEFAULT
                THROW_HARD('cartft_calc particle with an unsupported CTF flag')
        end select
        associate( p => self%ptcls(iptcl) )
            p%valid   = .false.
            p%objfun  = objfun
            p%power   = 0._dp
            p%kfromto = [lower_shell, upper_shell]
            p%half    = 0
            p%npix    = 0
            if( upper_shell < lower_shell ) return
            if( present(half) )then
                if( .not. all(self%sigma2_ok(lower_shell:upper_shell,half)) ) return
                p%half = half
            endif
            ! the slot's arrays are kept across fills (the box is fixed for the calculator's lifetime);
            ! every reader visits only the pixels of the slot's shells, which the fill below writes
            if( .not. allocated(p%observed) ) &
                &allocate(p%observed(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2)), source=cmplx(0.,0.))
            if( .not. allocated(p%transfer) ) &
                &allocate(p%transfer(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2)), source=0.)
            use_ctf = ctfparms%ctfflag /= CTFFLAG_NO
            if( use_ctf )then
                tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
                ctfvals = tfun%get_ctfvars(ctfparms%phshift)
            endif
            half_power = 0._dp
            dc_power   = 0._dp
            do k = self%lims2(2,1), self%lims2(2,2)
                do h = self%lims2(1,1), self%lims2(1,2)
                    shell = self%shells(h,k)
                    if( shell < p%kfromto(1) .or. shell > p%kfromto(2) ) cycle
                    weight = 1.0
                    if( present(half) ) weight = self%rsqrt_sigma2(shell,half)
                    cval = 1.0
                    ! SIMPLE's CTF object without the process-global Fourier maps of the matcher;
                    ! abs(CTF): the observation is phase-flipped (CTFFLAG_YES, O4) or stored flipped
                    if( use_ctf ) cval = abs(ctf_value(tfun, h, k, self%box, ctfvals%phshift))
                    p%transfer(h,k) = cval*weight
                    p%observed(h,k) = observed(h,k)*weight
                    p%npix = p%npix + 1
                    ! the power the kernels see: the half plane twice, in their order, and DC once
                    if( h == 0 .and. k == 0 )then
                        dc_power = abs(cmplx(p%observed(h,k), kind=dp))**2
                    else if( h > 0 .or. (h == 0 .and. k > 0) )then
                        half_power = half_power + abs(cmplx(p%observed(h,k), kind=dp))**2
                    endif
                end do
            end do
            p%power = 2._dp*half_power + dc_power
            ! both objectives normalize by the particle power: an observation without power
            ! has no defined objective
            p%valid = p%power > NUMERIC_FLOOR
        end associate
    end subroutine fill_ptcl

    !> Continuous (B-spline) gather of one reference and its three spatial derivatives at one
    !! oversampled lattice location, in native Fourier scale (folded into the reference). Pure: a
    !! location outside the wrap table returns inside = .false. and zeros; the non-pure caller stops.
    pure subroutine gather( self, cmat, loc, val, dval, inside )
        class(cartft_calc), intent(in)  :: self
        complex,            intent(in)  :: cmat(:,:,:)
        real(sp),           intent(in)  :: loc(3)
        complex,            intent(out) :: val, dval(3)
        logical,            intent(out) :: inside
        call gather_packed_bspline3_window_grad(cmat, lbound(self%wrap,1), self%wrap, loc, val, dval, inside)
    end subroutine gather

    pure integer function ref_index( self, iseven )
        class(cartft_calc), intent(in) :: self
        logical,            intent(in) :: iseven
        ref_index = merge(1, 2, iseven)
    end function ref_index

    !> Caller contract of the objective methods: an existing reference and a valid particle
    !! slot whose shells lie in the native Fourier disk.
    subroutine check_eval_args( self, state, iseven, iptcl )
        class(cartft_calc), intent(in) :: self
        integer,            intent(in) :: state, iptcl
        logical,            intent(in) :: iseven
        if( .not. self%ref_exists(state, iseven) ) THROW_HARD('cartft_calc objective from an absent reference')
        if( .not. self%ptcl_is_valid(iptcl) )      THROW_HARD('cartft_calc objective requires a valid particle slot')
        associate( kfromto => self%ptcls(iptcl)%kfromto )
            if( kfromto(1) < 0 .or. kfromto(2) > self%box/2 .or. kfromto(2) < kfromto(1) ) &
                &THROW_HARD('cartft_calc shell range lies outside the native Fourier disk')
        end associate
    end subroutine check_eval_args

    !> The separable shift phase of one evaluation, exp(2 pi i (h s1 + k s2)/box) = eh(h) ek(k),
    !! over h = 0..kto and k = -kto..kto
    pure subroutine shift_phase_tables( self, shift, eh, ek )
        class(cartft_calc), intent(in)  :: self
        real(dp),           intent(in)  :: shift(2)
        complex(dp),        intent(out) :: eh(0:), ek(-(size(eh)-1):)
        real(dp) :: arg
        integer  :: i
        do i = 0, size(eh) - 1
            arg   = 2._dp*real(PI, dp)*real(i, dp)*shift(1)/real(self%box, dp)
            eh(i) = cmplx(cos(arg), sin(arg), kind=dp)
        end do
        do i = -(size(eh) - 1), size(eh) - 1
            arg   = 2._dp*real(PI, dp)*real(i, dp)*shift(2)/real(self%box, dp)
            ek(i) = cmplx(cos(arg), sin(arg), kind=dp)
        end do
    end subroutine shift_phase_tables

    !> The DC pixel of a band from shell 0, its own Friedel mate: no Jacobian, its objective and
    !! power terms once
    subroutine add_dc_pixel( self, p, cmat, l_euclid, objective, prediction_power, cross_real, dc_residual )
        class(cartft_calc), intent(in)    :: self
        type(cartft_ptcl),  intent(in)    :: p
        complex,            intent(in)    :: cmat(:,:,:)
        logical,            intent(in)    :: l_euclid
        real(dp),           intent(inout) :: objective, prediction_power, cross_real
        real(dp), optional, intent(out)   :: dc_residual
        complex     :: val, dval(3)
        complex(dp) :: model, particle, residual
        logical     :: inside
        if( present(dc_residual) ) dc_residual = 0._dp
        call self%gather(cmat, [0.,0.,0.], val, dval, inside)
        model    = real(p%transfer(0,0), dp)*cmplx(val, kind=dp)
        particle = cmplx(p%observed(0,0), kind=dp)
        if( l_euclid )then
            residual  = model - particle
            objective = objective + real(conjg(residual)*residual, dp)
            if( present(dc_residual) ) dc_residual = real(conjg(residual)*residual, dp)
        else
            prediction_power = prediction_power + real(conjg(model)*model, dp)
            cross_real       = cross_real + real(conjg(particle)*model, dp)
        endif
    end subroutine add_dc_pixel

    !> Convert the whitened residual power sum |r|^2 and its Gauss-Newton sums into the
    !! normalized loss L = sum |r|^2/power, its gradient and Gauss-Newton block, the block scaled
    !! by the |f|^2-weighted slope of data on model (curvature_scale). NaN when the particle power vanishes.
    subroutine finalize_euclid_normal_terms( power, weighted_power, weighted_cross, objective, gradient, hessian )
        real(dp), intent(in)    :: power, weighted_power, weighted_cross
        real(dp), intent(inout) :: objective, gradient(:), hessian(:,:)
        if( power <= NUMERIC_FLOOR )then
            objective = ieee_value(0._dp, ieee_quiet_nan)
            gradient  = objective
            hessian   = objective
            return
        endif
        objective = objective/power
        gradient  = 2._dp*gradient/power
        hessian   = 2._dp*hessian/power
        if( weighted_power > NUMERIC_FLOOR ) hessian = curvature_scale(weighted_cross/weighted_power)*hessian
    end subroutine finalize_euclid_normal_terms

    !> Convert correlation sufficient statistics into the 1-correlation gradient and the
    !! Gauss-Newton block of the normalized model vector, scaled by the |f|^2-weighted slope of the
    !! normalized data on the normalized model (curvature_scale).
    subroutine finalize_cc_normal_terms( particle_power, prediction_power, cross_real, weighted_power, &
        &weighted_cross, model_gradient, objective, gradient, hessian )
        real(dp), intent(in)    :: particle_power, prediction_power, cross_real, weighted_power, weighted_cross
        real(dp), intent(in)    :: model_gradient(:)
        real(dp), intent(out)   :: objective
        real(dp), intent(inout) :: gradient(:), hessian(:,:)
        real(dp) :: normalizer
        integer  :: axis, jaxis
        if( particle_power <= NUMERIC_FLOOR .or. prediction_power <= NUMERIC_FLOOR )then
            objective = ieee_value(0._dp, ieee_quiet_nan)
            gradient  = objective
            hessian   = objective
            return
        endif
        normalizer = sqrt(particle_power*prediction_power)
        objective  = 1._dp - cross_real/normalizer
        gradient   = -gradient/normalizer + cross_real*model_gradient/(normalizer*prediction_power)
        do axis = 1, size(gradient)
            do jaxis = 1, size(gradient)
                hessian(axis,jaxis) = hessian(axis,jaxis)/prediction_power - &
                    &model_gradient(axis)*model_gradient(jaxis)/(prediction_power*prediction_power)
            end do
        end do
        if( weighted_power > NUMERIC_FLOOR ) hessian = &
            &curvature_scale(weighted_cross/weighted_power*sqrt(prediction_power/particle_power))*hessian
    end subroutine finalize_cc_normal_terms

    !> The expected curvature of a noisy objective is the Gauss-Newton matrix times the regression slope
    !! of data on model, weighted by |f|^2 as the matrix is: a reference filtered at high frequency
    !! needs more curvature than its overall amplitude suggests (cont_lm_acceleration.md, section 3.1).
    pure real(dp) function curvature_scale( slope )
        real(dp), intent(in) :: slope
        curvature_scale = min(CURVATURE_SCALE_MAX, max(CURVATURE_SCALE_MIN, slope))
    end function curvature_scale

    !> SIMPLE's CTF object at one signed full-disk Fourier coordinate. Unlike the memoized
    !! hot-loop kernel it owns no process-global Fourier maps, so it is valid in standalone tests.
    real function ctf_value( tfun, h, k, box, phshift ) result( cval )
        type(ctf), intent(in) :: tfun
        integer,   intent(in) :: h, k, box
        real,      intent(in) :: phshift
        real :: angle, spatial_frequency_squared
        spatial_frequency_squared = (real(h)*real(h) + real(k)*real(k))/(real(box)*real(box))
        angle = 0.
        if( h /= 0 .or. k /= 0 ) angle = atan2(real(k), real(h))
        cval = tfun%eval_canonical(spatial_frequency_squared, angle, phshift)
    end function ctf_value

end module simple_cartft_calc
