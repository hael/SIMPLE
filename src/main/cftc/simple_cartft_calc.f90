!@descr: Cartesian Fourier calculator of continuous pose refinement: references, prepared particles, objectives
!
! The Cartesian counterpart of polarft_calc (plan section 6.2). It holds the state-by-half
! references and the prepared particles of the current batch. The evaluation methods
! (predict, objective_gradient, shift_normal_terms, pose_normal_terms, sigma_contribution)
! take the calculator intent(in), so concurrent evaluation on different particles cannot
! write shared state. set_ref runs serially before the particle loop; set_ptcl writes only
! the slot it is given, so one thread per slot may run it concurrently.
!
! Objectives (Phase 3, plan section 5; C4-C6). A particle slot is prepared for one objective
! and every evaluation on it uses that objective, so the objective follows the caller's
! params%cc_objfun at preparation:
!  - OBJFUN_CC, set_ptcl without sigma2 (C5: sigma2 cannot enter): 1 - cc with
!    cc = sum Re(conjg(X) C M) / sqrt(sum |X|^2 sum |C M|^2), a uniform weight per Cartesian
!    pixel, the weighting of the continuous polar stages;
!  - OBJFUN_EUCLID, set_ptcl with sigma2: the polar loss L = sum |X - C M|^2/sigma2 /
!    sum |X|^2/sigma2 (per-pixel weight 1/sigma2 of the pixel's shell), score exp(-L) (C6).
! X is the observation as prepimg4align prepares it up to the mask (centred on the stored
! shift, phase-flipped for CTFFLAG_YES, masked) times the stencil's taper (O4, O5 (a),
! get_ptcl_taper); C = abs(CTF) for CTFFLAG_YES and CTFFLAG_FLIP, 1 without CTF; M = S(t) G(R) V
! the gathered reference with the shift phase S(t). A pixel belongs to shell nint(|(h,k)|)
! and is used when that shell lies in the slot's range, the polar ring membership, and the
! pixel lies inside the Nyquist circle.
! References reach the matcher as a file of prepared real-space volumes, one per half (plan
! section 6.6, C22): the reference materializer stages the volumes it has prepared
! (set_refvol) and writes them (write); every matcher reads them back (read), which pads and
! transforms each once (set_ref). The header carries the band limit of the iteration, which
! the reader returns for the caller to adopt.
! The particle slots are those of the current batch (new_ptcls sizes them; Phase 6): the batch
! preparation of the matcher (prep_cart_batch) fills slot i with the i-th particle of the batch,
! one particle per iteration of its parallel loop, and the sigma owner reads the residual of a
! slot (sigma_contribution).
module simple_cartft_calc
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
use simple_core_module_api,   only: dp, sp, PI, KBWINSZ, KBALPHA, OSMPL_PAD_FAC, CTFFLAG_NO, CTFFLAG_YES, CTFFLAG_FLIP, &
    &OBJFUN_CC, OBJFUN_EUCLID, ctfparams, ctfvars, kbinterpol, cyci_1d, simple_exception, string, fopen, fclose, &
    &fileiochk, file_exists, logfhandle
use simple_image,             only: image
use simple_ctf,               only: ctf
use simple_gridding,          only: kb_stencil_envelope_1d
use simple_cartesian_fourier, only: center_embed_real3d, gather_packed_window_grad
implicit none
private
public :: cartft_calc, CART_REFVOLS_FORMAT_VERSION

#include "simple_local_flags.inc"

real(dp), parameter :: NUMERIC_FLOOR = epsilon(1._dp)**2
integer,  parameter :: CART_REFVOLS_FORMAT_VERSION = 1 !< format_version of the reference-volume files (6.6)
integer,  parameter :: CART_REFVOLS_NHEADER = 5        !< [format_version, box_crop, nstates, kfrom, kto]

!> One reference: the padded Fourier transform of a prepared physical volume (matcher side),
!! or the staged prepared volume itself (materializer side, written by write).
type :: cartft_ref
    complex, allocatable :: cmat(:,:,:)
    real,    allocatable :: vol(:,:,:)
end type cartft_ref

!> One prepared particle for one objective: the observation X and the shift-free transfer
!! abs(C) over its shells, whitened by 1/sqrt(sigma2) under OBJFUN_EUCLID, unweighted under
!! OBJFUN_CC (no sigma2 is held).
type :: cartft_ptcl
    complex, allocatable :: observed(:,:), transfer(:,:)
    real,    allocatable :: sigma2(:)           !< OBJFUN_EUCLID only
    real(dp) :: power   = 0._dp                 !< sum |X|^2 (whitened under OBJFUN_EUCLID)
    integer  :: kfromto(2) = 0
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
    integer :: iwinsz  = 0                  !< integer half-width of the KB window
    integer :: wdim    = 0                  !< KB stencil width
    integer :: lims2(2,2) = 0               !< full-disk 2D Fourier limits
    real    :: padsc   = 1.                 !< padf**3, native Fourier scaling of a gather
    type(kbinterpol)               :: kbwin !< KB window of the gather
    integer,           allocatable :: wrap(:)     !< periodic wrap table of the padded grid
    type(cartft_ref),  allocatable :: refs(:,:)   !< (2,nstates): 1 even, 2 odd
    type(cartft_ptcl), allocatable :: ptcls(:)    !< (nptcls) prepared particles
    integer,           allocatable :: slots(:)    !< particle index -> slot of the current batch (0 when absent)
    logical :: exists = .false.
contains
    ! lifecycle
    procedure          :: new
    procedure          :: kill
    ! references
    procedure          :: set_ref
    procedure          :: ref_exists
    procedure          :: has_refs
    ! reference-volume file (6.6)
    procedure          :: set_refvol
    procedure          :: write
    procedure          :: read
    procedure          :: cart_refvols_header_compatible
    ! particles
    procedure, private :: set_ptcl_cc
    procedure, private :: set_ptcl_euclid
    generic            :: set_ptcl => set_ptcl_cc, set_ptcl_euclid
    procedure          :: ptcl_is_valid
    procedure          :: get_ptcl_kfromto
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
    procedure          :: sigma_contribution
    procedure          :: score
    procedure, private :: fill_ptcl
    procedure, private :: gather
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
        self%wdim    = 2*self%iwinsz + 1
        self%lims2(1,:) = [-box/2, box/2]
        self%lims2(2,:) = [-box/2, box/2]
        allocate(self%refs(2,nstates), self%ptcls(nptcls))
        self%exists = .true.
    end subroutine new

    !> Release every reference and particle; safe on a killed or never-built instance.
    subroutine kill( self )
        class(cartft_calc), intent(inout) :: self
        if( allocated(self%refs)  ) deallocate(self%refs)
        if( allocated(self%ptcls) ) deallocate(self%ptcls)
        if( allocated(self%wrap)  ) deallocate(self%wrap)
        if( allocated(self%slots) ) deallocate(self%slots)
        self%nstates = 0
        self%nptcls  = 0
        self%box     = 0
        self%boxpd   = 0
        self%padf    = 1
        self%iwinsz  = 0
        self%wdim    = 0
        self%lims2   = 0
        self%padsc   = 1.
        self%exists  = .false.
    end subroutine kill

    ! REFERENCES

    !> Build the reference of one state and half from a prepared physical volume
    !! (box**3, the volume the polar branch pads and reprojects): pad by OSMPL_PAD_FAC in
    !! real space and Fourier transform once. No inverse KB envelope is applied (C10).
    !! THROW_HARD for a state outside 1..nstates or a volume that is not box**3.
    subroutine set_ref( self, state, iseven, volume )
        class(cartft_calc), intent(inout) :: self
        integer,            intent(in)    :: state
        logical,            intent(in)    :: iseven
        real,               intent(in)    :: volume(:,:,:)
        type(image) :: padded_image
        integer     :: lims3(3,2), wlims(2), lo, hi, i, ihalf
        if( .not. self%exists ) THROW_HARD('cartft_calc reference set before new')
        if( state < 1 .or. state > self%nstates ) THROW_HARD('cartft_calc reference state out of range')
        if( any(shape(volume) /= self%box) ) THROW_HARD('cartft_calc reference volume does not match the box')
        ihalf = self%ref_index(iseven)
        call padded_image%new([self%boxpd, self%boxpd, self%boxpd], 1.0)
        call padded_image%set_rmat(center_embed_real3d(volume, self%boxpd), .false.)
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
        self%refs(ihalf,state)%cmat = padded_image%get_cmat()
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

    !> Prepare particle slot iptcl for OBJFUN_EUCLID: X/sqrt(sigma2) and abs(C)/sqrt(sigma2) of
    !! the pixel's shell over kfromto intersected with the sigma2 shells and the Cartesian
    !! Nyquist limit, and the whitened particle power that normalizes the loss. An empty range,
    !! a non-finite or non-positive sigma2 in it or an observation without power leaves the
    !! slot invalid. THROW_HARD for iptcl outside 1..nptcls or an unsupported CTF flag.
    subroutine set_ptcl_euclid( self, iptcl, observed, ctfparms, sigma2, kfromto )
        class(cartft_calc), intent(inout) :: self
        integer,            intent(in)    :: iptcl
        complex,            intent(in)    :: observed(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2))
        type(ctfparams),    intent(in)    :: ctfparms !< at the cropped sampling
        real,               intent(in)    :: sigma2(0:)
        integer,            intent(in)    :: kfromto(2)
        integer :: lower_shell, upper_shell
        if( .not. self%exists ) THROW_HARD('cartft_calc particle set before new')
        if( iptcl < 1 .or. iptcl > self%nptcls ) THROW_HARD('cartft_calc particle slot out of range')
        lower_shell = max(0, kfromto(1))
        upper_shell = min(kfromto(2), ubound(sigma2,1), self%box/2)
        call self%fill_ptcl(iptcl, OBJFUN_EUCLID, observed, ctfparms, lower_shell, upper_shell, sigma2)
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
                    if( .not. in_shells(h, k, [0, self%box/2], self%box) ) cycle
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
    !! for the optimizer's shift stage; it avoids derivative planes and a masked five-by-five
    !! solve when only the two image shifts are active. NaN objective when a correlation is
    !! undefined.
    subroutine shift_normal_terms( self, state, iseven, iptcl, rotmat, shift, objective, gradient, hessian )
        class(cartft_calc), intent(in)  :: self
        integer,            intent(in)  :: state, iptcl
        logical,            intent(in)  :: iseven
        real(dp),           intent(in)  :: rotmat(3,3), shift(2)
        real(dp),           intent(out) :: objective, gradient(2), hessian(2,2)
        complex     :: val, dval(3), phase
        complex(dp) :: model, particle, residual, jacobian(2)
        real(sp)    :: loc(3)
        real(dp)    :: l_arg, frequency(2), particle_power, prediction_power, cross_real
        real(dp)    :: model_gradient(2)
        integer     :: axis, h, jaxis, k, ihalf, objfun
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
            model_gradient   = 0._dp
            do k = self%lims2(2,1), self%lims2(2,2)
                do h = self%lims2(1,1), self%lims2(1,2)
                    if( .not. in_shells(h, k, p%kfromto, self%box) ) cycle
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
                    jacobian  = cmplx(0._dp, frequency, kind=dp)*model
                    select case(objfun)
                        case(OBJFUN_EUCLID)
                            objective = objective + real(conjg(residual)*residual, dp)
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
                    call finalize_euclid_normal_terms(p%power, objective, gradient, hessian)
                case(OBJFUN_CC)
                    call finalize_cc_normal_terms(particle_power, prediction_power, cross_real, &
                        &model_gradient, objective, gradient, hessian)
            end select
        end associate
    end subroutine shift_normal_terms

    !> Objective, five-vector gradient and 5x5 Gauss-Newton block (rotations 1:3, shifts
    !! 4:5), for the optimizer's joint stage. NaN objective when a correlation is undefined.
    subroutine pose_normal_terms( self, state, iseven, iptcl, rotmat, shift, objective, gradient, hessian )
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
        real(dp)    :: model_gradient(5)
        integer     :: axis, h, jaxis, k, ihalf, objfun
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
            model_gradient   = 0._dp
            do k = self%lims2(2,1), self%lims2(2,2)
                do h = self%lims2(1,1), self%lims2(1,2)
                    if( .not. in_shells(h, k, p%kfromto, self%box) ) cycle
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
                    jacobian(4:5) = cmplx(0._dp, frequency, kind=dp)*model
                    select case(objfun)
                        case(OBJFUN_EUCLID)
                            objective = objective + real(conjg(residual)*residual, dp)
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
                    call finalize_euclid_normal_terms(p%power, objective, gradient, hessian)
                case(OBJFUN_CC)
                    call finalize_cc_normal_terms(particle_power, prediction_power, cross_real, &
                        &model_gradient, objective, gradient, hessian)
            end select
        end associate
    end subroutine pose_normal_terms

    !> Unwhitened per-shell accounting at one terminal pose of an OBJFUN_EUCLID particle slot:
    !! residual power over two per sample (the sigma2 contribution), reference and particle
    !! power, over the particle's shells, and v = the normalized loss L (-1 when undefined).
    !! An invalid slot returns v = -1 and unallocated arrays. THROW_HARD for an OBJFUN_CC slot:
    !! a correlation pass writes no sigma2 (C5).
    subroutine sigma_contribution( self, state, iseven, iptcl, rotmat, shift, sigma_contrib, ref_pow, ptcl_pow, v )
        class(cartft_calc), intent(in)  :: self
        integer,            intent(in)  :: state, iptcl
        logical,            intent(in)  :: iseven
        real(dp),           intent(in)  :: rotmat(3,3), shift(2)
        real, allocatable,  intent(out) :: sigma_contrib(:), ref_pow(:), ptcl_pow(:)
        real,               intent(out) :: v
        complex     :: val, dval(3), phase
        complex(dp) :: model, raw_observed, residual
        real(dp), allocatable :: sigma_sum(:), ref_sum(:), ptcl_sum(:)
        real(dp)    :: l_arg, root_sigma, vnum, vden
        real(sp)    :: loc(3)
        integer, allocatable :: counts(:)
        integer     :: h, k, shell, lower_shell, upper_shell, ihalf
        logical     :: inside
        if( .not. self%ref_exists(state, iseven) ) THROW_HARD('cartft_calc sigma contribution from an absent reference')
        if( iptcl < 1 .or. iptcl > self%nptcls )   THROW_HARD('cartft_calc particle slot out of range')
        v = -1.
        if( .not. self%ptcls(iptcl)%valid ) return
        if( self%ptcls(iptcl)%objfun /= OBJFUN_EUCLID ) THROW_HARD('cartft_calc sigma contribution requires a euclid particle slot')
        ihalf = self%ref_index(iseven)
        associate( p => self%ptcls(iptcl), cmat => self%refs(ihalf,state)%cmat )
            lower_shell = p%kfromto(1)
            upper_shell = p%kfromto(2)
            allocate(sigma_sum(lower_shell:upper_shell), source=0._dp)
            allocate(ref_sum(lower_shell:upper_shell),   source=0._dp)
            allocate(ptcl_sum(lower_shell:upper_shell),  source=0._dp)
            allocate(counts(lower_shell:upper_shell),    source=0)
            vnum = 0._dp
            vden = 0._dp
            do k = self%lims2(2,1), self%lims2(2,2)
                do h = self%lims2(1,1), self%lims2(1,2)
                    if( .not. in_shells(h, k, p%kfromto, self%box) ) cycle
                    shell = nint(sqrt(real(h*h + k*k)))
                    loc   = real(self%padf, sp)*real(matmul(real([h, k, 0], dp), rotmat), sp)
                    call self%gather(cmat, loc, val, dval, inside)
                    if( .not. inside ) THROW_HARD('cartft_calc gather location outside the wrap table')
                    l_arg = 2._dp*real(PI, dp)*(real(h, dp)*shift(1) + real(k, dp)*shift(2))/real(self%box, dp)
                    phase = cmplx(cos(l_arg), sin(l_arg), kind=sp)
                    root_sigma   = sqrt(real(p%sigma2(shell), dp))
                    model        = cmplx(phase, kind=dp)*cmplx(p%transfer(h,k), kind=dp)*cmplx(val, kind=dp)*root_sigma
                    raw_observed = cmplx(p%observed(h,k), kind=dp)*root_sigma
                    residual     = raw_observed - model
                    sigma_sum(shell) = sigma_sum(shell) + real(conjg(residual)*residual, dp)
                    ref_sum(shell)   = ref_sum(shell)   + real(conjg(model)*model, dp)
                    ptcl_sum(shell)  = ptcl_sum(shell)  + real(conjg(raw_observed)*raw_observed, dp)
                    counts(shell)    = counts(shell) + 1
                    vnum = vnum + real(conjg(residual)*residual, dp)/real(p%sigma2(shell), dp)
                    vden = vden + real(conjg(raw_observed)*raw_observed, dp)/real(p%sigma2(shell), dp)
                end do
            end do
        end associate
        if( any(counts == 0) ) THROW_HARD('cartft_calc sigma contribution found an empty active shell')
        allocate(sigma_contrib(lower_shell:upper_shell), ref_pow(lower_shell:upper_shell), ptcl_pow(lower_shell:upper_shell))
        sigma_contrib = real(sigma_sum/(2._dp*real(counts, dp)), sp)
        ref_pow       = real(ref_sum/real(counts, dp), sp)
        ptcl_pow      = real(ptcl_sum/real(counts, dp), sp)
        if( vden > 0._dp )then
            v = real(vnum/vden, sp)
        else
            v = -1.
        endif
    end subroutine sigma_contribution

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

    !> Fill particle slot iptcl for objfun over the shells lower_shell..upper_shell; sigma2 is
    !! present exactly for OBJFUN_EUCLID (set_ptcl_euclid) and whitens observation and transfer.
    subroutine fill_ptcl( self, iptcl, objfun, observed, ctfparms, lower_shell, upper_shell, sigma2 )
        class(cartft_calc), intent(inout) :: self
        integer,            intent(in)    :: iptcl, objfun, lower_shell, upper_shell
        complex,            intent(in)    :: observed(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2))
        type(ctfparams),    intent(in)    :: ctfparms
        real, optional,     intent(in)    :: sigma2(0:)
        type(ctf)     :: tfun
        type(ctfvars) :: ctfvals
        real     :: cval, weight
        integer  :: h, k, shell
        logical  :: use_ctf
        select case(ctfparms%ctfflag)
            case(CTFFLAG_NO, CTFFLAG_YES, CTFFLAG_FLIP)
            case DEFAULT
                THROW_HARD('cartft_calc particle with an unsupported CTF flag')
        end select
        associate( p => self%ptcls(iptcl) )
            if( allocated(p%observed) ) deallocate(p%observed)
            if( allocated(p%transfer) ) deallocate(p%transfer)
            if( allocated(p%sigma2)   ) deallocate(p%sigma2)
            p%valid   = .false.
            p%objfun  = objfun
            p%power   = 0._dp
            p%kfromto = [lower_shell, upper_shell]
            if( upper_shell < lower_shell ) return
            if( present(sigma2) )then
                do shell = lower_shell, upper_shell
                    if( .not. ieee_is_finite(sigma2(shell)) .or. sigma2(shell) <= 0.0 ) return
                end do
                allocate(p%sigma2(lower_shell:upper_shell), source=sigma2(lower_shell:upper_shell))
            endif
            allocate(p%observed(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2)), source=cmplx(0.,0.))
            allocate(p%transfer(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2)), source=cmplx(0.,0.))
            use_ctf = ctfparms%ctfflag /= CTFFLAG_NO
            if( use_ctf )then
                tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
                ctfvals = tfun%get_ctfvars(ctfparms%phshift)
            endif
            do k = self%lims2(2,1), self%lims2(2,2)
                do h = self%lims2(1,1), self%lims2(1,2)
                    if( .not. in_shells(h, k, p%kfromto, self%box) ) cycle
                    shell  = nint(sqrt(real(h*h + k*k)))
                    weight = 1.0
                    if( present(sigma2) ) weight = 1.0/sqrt(sigma2(shell))
                    cval = 1.0
                    ! SIMPLE's CTF object without the process-global Fourier maps of the matcher;
                    ! abs(CTF): the observation is phase-flipped (CTFFLAG_YES, O4) or stored flipped
                    if( use_ctf ) cval = abs(ctf_value(tfun, h, k, self%box, ctfvals%phshift))
                    p%transfer(h,k) = cval*weight
                    p%observed(h,k) = observed(h,k)*weight
                    p%power = p%power + real(conjg(p%observed(h,k))*p%observed(h,k), dp)
                end do
            end do
            ! both objectives normalize by the particle power: an observation without power
            ! has no defined objective
            p%valid = p%power > NUMERIC_FLOOR
        end associate
    end subroutine fill_ptcl

    !> KB gather of one reference and its three spatial derivatives at one oversampled
    !! lattice location, times padf**3. Pure: a location outside the wrap table returns
    !! inside = .false. and zeros, and the non-pure caller stops (formerly an error stop here).
    pure subroutine gather( self, cmat, loc, val, dval, inside )
        class(cartft_calc), intent(in)  :: self
        complex,            intent(in)  :: cmat(:,:,:)
        real(sp),           intent(in)  :: loc(3)
        complex,            intent(out) :: val, dval(3)
        logical,            intent(out) :: inside
        real(sp) :: w(self%wdim,self%wdim,self%wdim), dw(self%wdim,self%wdim,self%wdim,3), switch_margin(3)
        integer  :: i0(3)
        ! w and dw/dloc on the same fixed interpolation stencil
        call self%kbwin%apod_mat_3d_fast_grad(loc, self%iwinsz, self%wdim, i0, switch_margin, w, dw)
        inside = .not.( any(i0 < lbound(self%wrap,1)) .or. any(i0 + self%wdim - 1 > ubound(self%wrap,1)) )
        if( .not. inside )then
            val  = cmplx(0.,0.)
            dval = cmplx(0.,0.)
            return
        endif
        call gather_packed_window_grad(cmat, lbound(self%wrap,1), self%wrap, i0, w, dw, val, dval)
        ! native Fourier scaling of the value and all derivatives
        val  = self%padsc*val
        dval = self%padsc*dval
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

    !> True when pixel (h,k) belongs to a shell nint(|(h,k)|) inside kfromto, the polar ring
    !! membership, and lies inside the Nyquist circle |(h,k)| <= box/2: the outer half of the
    !! Nyquist shell is beyond the band of the padded reference transform (and beyond its wrap
    !! table), where the polar ring of radius box/2 does not sample either.
    pure logical function in_shells( h, k, kfromto, box )
        integer, intent(in) :: h, k, kfromto(2), box
        integer :: shell
        shell     = nint(sqrt(real(h*h + k*k)))
        in_shells = shell >= kfromto(1) .and. shell <= kfromto(2) .and. h*h + k*k <= (box/2)**2
    end function in_shells

    !> Convert the whitened residual power sum |r|^2 and its Gauss-Newton sums into the
    !! normalized loss L = sum |r|^2/power, its gradient and Gauss-Newton block. NaN when the
    !! particle power vanishes.
    subroutine finalize_euclid_normal_terms( power, objective, gradient, hessian )
        real(dp), intent(in)    :: power
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
    end subroutine finalize_euclid_normal_terms

    !> Convert correlation sufficient statistics into the 1-correlation gradient and the
    !! Gauss-Newton block of the normalized model vector.
    subroutine finalize_cc_normal_terms( particle_power, prediction_power, cross_real, &
        &model_gradient, objective, gradient, hessian )
        real(dp), intent(in)    :: particle_power, prediction_power, cross_real
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
    end subroutine finalize_cc_normal_terms

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
