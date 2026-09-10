!@descr: the abstract data type for sigma2 used when objfun=euclid
module simple_euclid_sigma2
use, intrinsic :: iso_fortran_env, only: int64, real32
use simple_core_module_api
use simple_polarft_calc,   only: polarft_calc
use simple_parameters,     only: parameters
use simple_sigma2_state,   only: sigma2_state_candidate_path, sigma2_state_range_path, sigma2_state_next_generation
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_read_header, &
    &sigma2_state_read_groups, sigma2_state_read_particles, sigma2_state_write_local_range
use simple_starfile_wrappers
implicit none

public :: euclid_sigma2, sigma2_group_iter
! grouped-STAR I/O for the explicit sigma2_convert boundary only
public :: write_groups_starfile, read_sigma2_groups_file
private
#include "simple_local_flags.inc"

integer, parameter :: LENSTR = 48
! euclid scale diagnostics (doc/implementation_notes/drop_legacy_box_division.md, plan step 1):
! the search band is split into NDIAG_BANDS contiguous bands; per particle we keep the
! reference/particle amplitude ratio per band and the euclid objective value v at the
! assigned orientation, and report quantiles once per iteration
integer, parameter :: NDIAG_BANDS = 4

type euclid_sigma2
    private
    class(parameters),    pointer :: p_ptr => null()
    real,    allocatable, public  :: sigma2_noise(:,:)      !< the sigmas for alignment & reconstruction (from groups)
    real,    allocatable          :: sigma2_part(:,:)       !< the actual sigmas per particle (this part only)
    real,    allocatable          :: sigma2_groups(:,:,:)   !< sigmas for groups
    integer, allocatable          :: pinds(:)
    integer, allocatable          :: micinds(:)
    real,    allocatable          :: diag_ratio(:,:)       !< ref/ptcl amplitude ratio per band & particle (this part only)
    real,    allocatable          :: diag_v(:)             !< euclid objective value at assigned orientation (this part only)
    integer                       :: fromp
    integer                       :: top
    integer                       :: kfromto(2) = 0
    type(string)                  :: binfname
    logical                       :: exists     = .false.
contains
    ! constructor
    procedure          :: new
    procedure, private :: init_from_group_header
    ! utils
    procedure          :: write_info
    procedure          :: get_kfromto
    procedure          :: set_kfromto
    ! I/O
    procedure          :: read_part
    procedure          :: read_groups
    procedure          :: allocate_ptcls
    procedure          :: calc_sigma2
    procedure          :: write_sigma2
    procedure          :: report_euclid_diag
    procedure, private :: read_sigma2_groups
    ! destructor
    procedure          :: kill
end type euclid_sigma2

contains

    pure integer function sigma2_group_iter( matcher_iter, matcher_completed ) result( group_iter )
        integer, intent(in) :: matcher_iter
        logical, intent(in) :: matcher_completed
        group_iter = matcher_iter
        if( matcher_completed ) group_iter = group_iter + 1
    end function sigma2_group_iter

    subroutine new( self, params, pftc, binfname, box )
        ! read individual sigmas from binary file, to be modified at the end of the iteration
        ! read group sigmas from starfile, to be used for alignment and volume reconstruction
        ! set up fields for fast access to sigmas
        class(euclid_sigma2), target, intent(inout) :: self
        class(parameters),    target, intent(in)    :: params
        class(polarft_calc),          intent(inout) :: pftc
        class(string),                intent(in)    :: binfname
        integer,                      intent(in)    :: box
        call self%kill
        self%p_ptr => params
        self%kfromto = [1, fdim(box)-1]
        allocate( self%sigma2_noise(self%kfromto(1):self%kfromto(2),self%p_ptr%fromp:self%p_ptr%top))
        call pftc%assign_sigma2_noise(self%sigma2_noise)
        self%binfname     =  binfname
        self%fromp        =  self%p_ptr%fromp
        self%top          =  self%p_ptr%top
        self%sigma2_noise =  0.
        ! scale diagnostics, filled by calc_sigma2, reported & reset by write_sigma2 (euclid_diag=yes)
        if( self%p_ptr%l_euclid_diag )then
            allocate(self%diag_ratio(NDIAG_BANDS,self%fromp:self%top), self%diag_v(self%fromp:self%top), source=-1.)
        endif
        self%exists       =  .true.
    end subroutine new

    !>  This is a minimal constructor to allow I/O of groups
    subroutine init_from_group_header( self, fname )
        class(euclid_sigma2), target, intent(inout) :: self
        class(string),                intent(in)    :: fname
        type(string), allocatable  :: names(:)
        type(starfile_table_type)  :: table
        integer                    :: kfromto(2), ngroups
        logical                    :: l
        if (.not. file_exists(fname)) then
            THROW_HARD('euclid_sigma2_starfile: init_from_group_header; file does not exists: ' // fname%to_char())
        end if
        call starfile_table__new(table)
        call starfile_table__getnames(table, fname, names)
        call starfile_table__read( table, fname, names(1)%to_char() )
        l = starfile_table__getValue_int(table, EMDL_MLMODEL_NR_GROUPS, ngroups)
        l = starfile_table__getValue_int(table, EMDL_SPECTRAL_IDX,  kfromto(1))
        l = starfile_table__getValue_int(table, EMDL_SPECTRAL_IDX2, kfromto(2))
        self%kfromto = kfromto
        call starfile_table__delete(table)
    end subroutine init_from_group_header

    subroutine write_info(self)
        class(euclid_sigma2), intent(in) :: self
        write(logfhandle,*) 'kfromto: ',self%kfromto
        write(logfhandle,*) 'fromp:   ',self%fromp
        write(logfhandle,*) 'top:     ',self%top
    end subroutine write_info

    pure function get_kfromto( self )result( kfromto )
        class(euclid_sigma2), intent(in) :: self
        integer :: kfromto(2)
        kfromto = self%kfromto
    end function get_kfromto

    !>  Set the reconstruction/alignment band. Used by workflows (e.g. flex_pca)
    !!  that reconstruct at a cropped box without a loaded sigma2 group table, so the
    !!  reconstruction backend uses the correct box_crop spectral range.
    subroutine set_kfromto( self, kfromto )
        class(euclid_sigma2), intent(inout) :: self
        integer,              intent(in)    :: kfromto(2)
        self%kfromto = kfromto
    end subroutine set_kfromto

    ! I/O

    subroutine read_part( self, os )
        class(euclid_sigma2), intent(inout) :: self
        class(oris),          intent(inout) :: os
        real(real32), allocatable :: state_part(:,:)
        integer :: status
        character(len=STDLEN) :: message
        call sigma2_state_read_particles(self%binfname%to_char(), self%fromp, self%top, &
            &state_part, status, message)
        if( status /= 0 ) THROW_HARD(trim(message))
        if( size(state_part,1) /= self%kfromto(2)-self%kfromto(1)+1 ) &
            &THROW_HARD('canonical particle sigma2 has incompatible shell bounds')
        allocate(self%sigma2_part(self%kfromto(1):self%kfromto(2),self%fromp:self%top))
        self%sigma2_part = real(state_part)
        deallocate(state_part)
    end subroutine read_part

    subroutine read_groups( self, os )
        class(euclid_sigma2), intent(inout) :: self
        class(oris),          intent(inout) :: os
        integer                             :: iptcl, igroup, ngroups, eo
        real(real32), allocatable           :: state_groups(:,:,:)
        integer                             :: status, ishell
        character(len=STDLEN)               :: message
        if( .not.associated(self%p_ptr) )then
            THROW_HARD('euclid_sigma2: params pointer is not set')
        endif
        call sigma2_state_read_groups(self%binfname%to_char(), state_groups, status, message)
        if( status /= 0 ) THROW_HARD(trim(message))
        ngroups = size(state_groups,3)
        if( size(state_groups,1) /= self%kfromto(2)-self%kfromto(1)+1 ) &
            &THROW_HARD('canonical grouped sigma2 has incompatible shell bounds')
        allocate(self%sigma2_groups(2,ngroups,self%kfromto(1):self%kfromto(2)))
        do igroup = 1, ngroups
            do eo = 1, 2
                do ishell = self%kfromto(1), self%kfromto(2)
                    self%sigma2_groups(eo,igroup,ishell) = &
                        &real(state_groups(ishell-self%kfromto(1)+1,eo,igroup))
                enddo
            enddo
        enddo
        deallocate(state_groups)
        if( self%p_ptr%l_sigma_glob )then
            if( ngroups /= 1 ) THROW_HARD('ngroups must be 1 when global sigma is estimated (p_ptr%l_sigma_glob == .true.)')
            ! copy global sigma to particles
            !$omp parallel do default(shared) private(iptcl,eo) proc_bind(close) schedule(static)
            do iptcl = self%p_ptr%fromp, self%p_ptr%top
                if(os%get_state(iptcl) == 0 ) cycle
                eo = nint(os%get(iptcl, 'eo')) ! 0/1
                self%sigma2_noise(:,iptcl) = self%sigma2_groups(eo+1,1,:)
            end do
            !$omp end parallel do
        else
            ! copy group sigmas to particles
            !$omp parallel do default(shared) private(iptcl,eo,igroup) proc_bind(close) schedule(static)
            do iptcl = self%p_ptr%fromp, self%p_ptr%top
                if(os%get_state(iptcl) == 0 ) cycle
                igroup = os%get_int(iptcl, 'stkind')
                eo     = os%get_eo(iptcl)  ! 0/1
                self%sigma2_noise(:,iptcl) = self%sigma2_groups(eo+1,igroup,:)
            end do
            !$omp end parallel do
        endif
    end subroutine read_groups

    !>  allocate sigma2_part
    subroutine allocate_ptcls( self )
        class(euclid_sigma2), intent(inout) :: self
        if( .not.self%exists ) THROW_HARD('euclid_sigma2 has not been instanciated! allocate_ptcls_from_groups')
        allocate(self%sigma2_part(self%kfromto(1):self%kfromto(2),self%fromp:self%top),&
            &source=self%sigma2_noise)
    end subroutine allocate_ptcls

    !>  Calculates and updates sigma2 within search resolution range
    subroutine calc_sigma2( self, pftc, iptcl, o, refkind )
        class(euclid_sigma2), intent(inout) :: self
        class(polarft_calc),  intent(inout) :: pftc
        integer,              intent(in)    :: iptcl
        class(ori),           intent(in)    :: o
        character(len=*),     intent(in)    :: refkind ! 'proj' or 'class'
        integer :: iref, irot, kfromto(2), nk, ib, klo, khi
        real, allocatable :: sigma_contrib(:), ref_pow(:), ptcl_pow(:)
        real    :: shvec(2), v, rsum, psum
        if( .not.associated(self%p_ptr) )then
            THROW_HARD('euclid_sigma2: params pointer is not set')
        endif
        if ( o%isstatezero() ) return
        kfromto = self%p_ptr%kfromto
        allocate(sigma_contrib(kfromto(1):kfromto(2)), source=0.)
        shvec = o%get_2Dshift()
        iref  = nint(o%get(trim(refkind)))
        if( trim(refkind) == 'proj' .and. self%p_ptr%nstates > 1 )then
            iref = (o%get_state() - 1) * self%p_ptr%nspace + iref
        endif
        irot  = pftc%get_roind(360. - o%e3get())
        if( allocated(self%diag_v) )then
            allocate(ref_pow(kfromto(1):kfromto(2)), ptcl_pow(kfromto(1):kfromto(2)), source=0.)
            call pftc%gen_sigma_contrib(iref, iptcl, shvec, irot, sigma_contrib, ref_pow, ptcl_pow, v)
        else
            call pftc%gen_sigma_contrib(iref, iptcl, shvec, irot, sigma_contrib)
        endif
        self%sigma2_part(kfromto(1):kfromto(2),iptcl) = sigma_contrib
        ! scale diagnostics
        if( allocated(self%diag_v) )then
            self%diag_v(iptcl) = v
            nk = kfromto(2) - kfromto(1) + 1
            do ib = 1, NDIAG_BANDS
                klo  = kfromto(1) + nint(real(ib-1) * real(nk) / real(NDIAG_BANDS))
                khi  = kfromto(1) + nint(real(ib)   * real(nk) / real(NDIAG_BANDS)) - 1
                khi  = min(khi, kfromto(2))
                if( khi < klo ) cycle
                rsum = sum(ref_pow(klo:khi))
                psum = sum(ptcl_pow(klo:khi))
                if( psum > 0. ) self%diag_ratio(ib,iptcl) = sqrt(rsum / psum)
            enddo
        endif
        deallocate(sigma_contrib)
        if( allocated(ref_pow)  ) deallocate(ref_pow)
        if( allocated(ptcl_pow) ) deallocate(ptcl_pow)
    end subroutine calc_sigma2

    subroutine write_sigma2( self )
        class(euclid_sigma2), intent(inout) :: self
        type(sigma2_state_header) :: header
        type(string) :: candidate_path, range_path
        integer(int64) :: next_gen
        integer :: status
        character(len=STDLEN) :: message
        ! transaction-scoped names: the candidate and this range carry the
        ! generation the master's update will commit (2026-09-07)
        call sigma2_state_next_generation(self%binfname%to_char(), next_gen, status, message)
        if( status /= 0 ) THROW_HARD(trim(message))
        candidate_path = sigma2_state_candidate_path(self%binfname%to_char(), next_gen)
        range_path = sigma2_state_range_path(self%binfname%to_char(), next_gen, self%p_ptr%part, self%p_ptr%numlen)
        call sigma2_state_read_header(candidate_path%to_char(), header, status, message)
        if( status /= 0 ) THROW_HARD(trim(message))
        if( header%generation /= next_gen ) THROW_HARD('canonical sigma2 candidate belongs to another transaction')
        call sigma2_state_write_local_range(range_path%to_char(), header%generation, header%layout_digest, &
            &self%fromp, real(self%sigma2_part,real32), self%kfromto(1), self%kfromto(2), status, message)
        if( status /= 0 ) THROW_HARD(trim(message))
        call candidate_path%kill
        call range_path%kill
        call self%report_euclid_diag
        ! reset so that the next iteration's report covers only the particles it updates
        if( allocated(self%diag_ratio) ) self%diag_ratio = -1.
        if( allocated(self%diag_v)     ) self%diag_v     = -1.
    end subroutine write_sigma2

    !>  Once-per-iteration report of the reference/particle amplitude ratio per band and
    !>  the quantiles of the euclid objective value v at the assigned orientations.
    !>  v = sum_k (k/sigma2_k) sum_p |ptcl - CTF*ref|^2 / sum_k (k/sigma2_k) sum_p |ptcl|^2, so a
    !>  reference that explains particle variance gives v < 1; v ~ 1.000 throughout means the
    !>  reference barely enters the residual (refs << ptcls), v > 1 means it adds more power than
    !>  it explains. Healthy target (drop_legacy_box_division.md S3): ratios ~0.1-0.5 falling
    !>  with resolution; v clearly below 1, never ~1.000 throughout, never near the threshold.
    subroutine report_euclid_diag( self )
        class(euclid_sigma2), intent(in) :: self
        real, allocatable :: vals(:)
        real    :: q(NDIAG_BANDS), vq(3), vmax, vthres
        integer :: kfromto(2), nk, ib, klo, khi, n, ninvalid
        character(len=:), allocatable :: str
        if( .not.allocated(self%diag_v) ) return
        if( self%p_ptr%part /= 1 ) return   ! one report per iteration in distributed execution
        kfromto = self%p_ptr%kfromto
        nk      = kfromto(2) - kfromto(1) + 1
        vthres  = real(-log(real(TINY,dp)), kind=kind(vthres))
        ! bands
        str = ''
        do ib = 1, NDIAG_BANDS
            klo = kfromto(1) + nint(real(ib-1) * real(nk) / real(NDIAG_BANDS))
            khi = min(kfromto(1) + nint(real(ib) * real(nk) / real(NDIAG_BANDS)) - 1, kfromto(2))
            call valid_vals(self%diag_ratio(ib,:), vals, n)
            q(ib) = quantile(vals, n, 0.5)
            str = str//' k['//int2str(klo)//'-'//int2str(khi)//']: '//real2str_diag(q(ib))
        enddo
        call valid_vals(self%diag_v, vals, n)
        if( n == 0 ) return
        vq(1)    = quantile(vals, n, 0.05)
        vq(2)    = quantile(vals, n, 0.50)
        vq(3)    = quantile(vals, n, 0.95)
        vmax     = vals(n)
        ninvalid = count(vals(1:n) > vthres)
        if( self%p_ptr%nparts > 1 ) str = str//' [PART 1/'//int2str(self%p_ptr%nparts)//' ONLY]'
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A)') '>>> EUCLID DIAG ITER ', self%p_ptr%which_iter, &
            &' NPTCLS ', n, ' KFROMTO ', kfromto(1), '-', kfromto(2), ' REF/PTCL AMP (q50)'//str
        write(logfhandle,'(A,I0,A,F0.4,A,F0.4,A,F0.4,A,F0.4,A,F0.2,A,I0)') '>>> EUCLID DIAG ITER ', &
            &self%p_ptr%which_iter, ' V q05: ', vq(1), ' q50: ', vq(2), ' q95: ', vq(3), ' max: ', vmax, &
            &' THRES: ', vthres, ' NINVALID: ', ninvalid
        ! maps written before 2026-08 carry a 1/box amplitude convention; such a starting
        ! reference reprojects ~box times below the particle signal (v ~ 1.000 throughout)
        if( q(1) >= 0. .and. q(1) < 0.02 .and. vq(2) > 0.99 )then
            str = 'EUCLID DIAG: reference amplitudes ~'//real2str_diag(q(1))//' of the particle signal; '//&
                &'an old-convention (1/box) or otherwise mis-scaled reference volume escaped the automatic '//&
                &'rescaling in reference preparation; alignment against it is unreliable -- scale the input '//&
                &'volume by the box size and restart'
            THROW_WARN(str)
        endif
        if( allocated(vals) ) deallocate(vals)

        contains

            !> copies the non-negative entries into a sorted array
            subroutine valid_vals( arr, vals, n )
                real,              intent(in)    :: arr(:)
                real, allocatable, intent(inout) :: vals(:)
                integer,           intent(out)   :: n
                integer :: i
                if( allocated(vals) ) deallocate(vals)
                n = count(arr >= 0.)
                allocate(vals(max(1,n)), source=0.)
                if( n == 0 ) return
                n = 0
                do i = 1, size(arr)
                    if( arr(i) >= 0. )then
                        n = n + 1
                        vals(n) = arr(i)
                    endif
                enddo
                call hpsort(vals(1:n))
            end subroutine valid_vals

            real function quantile( vals, n, frac )
                real,    intent(in) :: vals(:)
                integer, intent(in) :: n
                real,    intent(in) :: frac
                if( n == 0 )then
                    quantile = -1.
                else
                    quantile = vals(max(1, min(n, nint(frac * real(n) + 0.5))))
                endif
            end function quantile

            function real2str_diag( r ) result( str )
                real, intent(in) :: r
                character(len=:), allocatable :: str
                character(len=32) :: buf
                if( r < 0. )then
                    str = 'n/a'
                else
                    write(buf,'(ES9.3)') r
                    str = trim(adjustl(buf))
                endif
            end function real2str_diag

    end subroutine report_euclid_diag

    subroutine write_groups_starfile( fname, group_pspecs, ngroups )
        class(string),     intent(in) :: fname
        real,              intent(in) :: group_pspecs(:,:,:)
        integer,           intent(in) :: ngroups
        type(string)                  :: stmp
        integer                       :: kfromto(2), eo, igroup, idx
        type(starfile_table_type)     :: ostar
        call starfile_table__new(ostar)
        call starfile_table__open_ofile(ostar, fname%to_char())
        ! global fields
        kfromto(1) = lbound(group_pspecs,3)
        kfromto(2) = ubound(group_pspecs,3)
        call starfile_table__addObject(ostar)
        call starfile_table__setIsList(ostar, .true.)
        call starfile_table__setname(ostar, "general")
        call starfile_table__setValue_int(ostar, EMDL_MLMODEL_NR_GROUPS, ngroups)
        call starfile_table__setValue_int(ostar, EMDL_SPECTRAL_IDX,  kfromto(1))
        call starfile_table__setValue_int(ostar, EMDL_SPECTRAL_IDX2, kfromto(2))
        call starfile_table__write_ofile(ostar)
        ! values
        do eo = 1, 2
            if( eo == 1 )then
                stmp = 'even'
            else
                stmp = 'odd'
            end if
            do igroup = 1, ngroups
                call starfile_table__clear(ostar)
                call starfile_table__setComment(ostar, stmp%to_char() // ', group ' // trim(int2str(igroup)) )
                call starfile_table__setName(ostar, trim(int2str(eo)) // '_group_' // trim(int2str(igroup)) )
                call starfile_table__setIsList(ostar, .false.)
                do idx = kfromto(1), kfromto(2)
                    call starfile_table__addObject(ostar)
                    call starfile_table__setValue_int(ostar,    EMDL_SPECTRAL_IDX, idx)
                    call starfile_table__setValue_double(ostar, EMDL_MLMODEL_SIGMA2_NOISE,&
                        real(group_pspecs(eo,igroup,idx),dp) )
                end do
                call starfile_table__write_ofile(ostar)
            end do
        end do
        call starfile_table__close_ofile(ostar)
        call starfile_table__delete(ostar)
    end subroutine write_groups_starfile

    ! Hard-coded reader of the grouped sigma2 STAR written by write_groups_starfile
    subroutine read_sigma2_groups( self, fname, pspecs, ngroups )
        class(euclid_sigma2),          intent(inout) :: self
        class(string),                 intent(in)    :: fname
        real,             allocatable, intent(out)   :: pspecs(:,:,:)
        integer,                       intent(out)   :: ngroups
        character(len=LENSTR), allocatable :: strings(:)
        character(len=LENSTR) :: line, string
        real(dp) :: dval
        integer  :: kfromto(2), i, l, funit, iostat, group, eo, idx, igroup, ieo
        if(.not.file_exists(fname))then
            THROW_HARD('File: '//fname%to_char()//' Does not exists; read_sigma2_groups')
        endif
        call fopen(funit, fname, action='READ',status='OLD', form='FORMATTED', iostat=iostat)
        call fileiochk('read_sigma2_groups: '//fname%to_char(), iostat)
        ! read header
        ! # of groups
        read(funit,fmt='(A)') line
        read(funit,fmt='(A)') line
        if( trim(line).ne.'data_general' )then
            THROW_HARD('Unrecognized formatting: '//fname%to_char()//'; read_sigma2_groups')
        endif
        read(funit,fmt='(A)') line
        read(funit,fmt='(A)') line
        call parse_key_int_pair(line, '_rlnNrGroups', ngroups)
        ! sprectral range
        read(funit,fmt='(A)') line
        call parse_key_int_pair(line, '_rlnSpectralIndex', kfromto(1))
        read(funit,fmt='(A)') line
        call parse_key_int_pair(line, '_rlnSpectralIndex2', kfromto(2))
        if( any(kfromto-self%kfromto < 0) )then
            print *,kfromto,self%kfromto
            THROW_HARD('Incorrect resolution range: read_sigma2_groups')
        endif
        ! parse data
        allocate(strings(kfromto(1):kfromto(2)), pspecs(2,ngroups,self%kfromto(1):self%kfromto(2)))
        do eo = 1,2
            do group = 1,ngroups
                do
                    read(funit,fmt='(A)') line
                    if( line(1:5).eq.'data_' ) exit
                enddo
                l = len_trim(line)
                i = index(line(6:l), '_')
                ieo = str2int(line(6:6+i-2), iostat)
                call fileiochk('Invalid formatting: '//trim(line)//'; read_sigma2_groups', iostat)
                if( ieo /= eo ) THROW_HARD('Invalid formatting: '//trim(line)//'; read_sigma2_groups')
                if( line(6+i-1:6+i+5).ne.'_group_')THROW_HARD('Invalid formatting: '//trim(line)//'; read_sigma2_groups')
                i = index(line, '_', back=.true.)
                igroup = str2int(line(i+1:l), iostat)
                call fileiochk('Invalid formatting: '//trim(line)//'; read_sigma2_groups', iostat)
                if( igroup /= group ) THROW_HARD('Invalid formatting: '//trim(line)//'; read_sigma2_groups')
                do
                    read(funit,fmt='(A)') line
                    if( trim(line).eq.'loop_' ) exit
                enddo
                read(funit,fmt='(A)') line
                call parse_key_string_pair(line, '_rlnSpectralIndex', string)
                if( trim(string).ne.'#1' ) THROW_HARD('Invalid formatting: '//trim(line)//'; read_sigma2_groups')
                read(funit,fmt='(A)') line
                call parse_key_string_pair(line, '_rlnSigma2Noise', string)
                if( trim(string).ne.'#2' ) THROW_HARD('Invalid formatting: '//trim(line)//'; read_sigma2_groups')
                ! read all frequencies
                do i = kfromto(1),kfromto(2)
                    read(funit,fmt='(A)') strings(i)
                enddo
                ! parse only relevent frequencies
                do i = self%kfromto(1),self%kfromto(2)
                    line = strings(i)
                    call split(line,' ',string)
                    read(string,*,iostat=iostat) idx
                    call fileiochk('Unrecognized formatting: '//trim(line)//'; read_sigma2_groups', iostat)
                    read(line,*,iostat=iostat) dval
                    call fileiochk('Unrecognized formatting: '//trim(line)//'; read_sigma2_groups', iostat)
                    if( i /= idx ) THROW_HARD('Invalid formatting: '//trim(line)//'; read_sigma2_groups')
                    pspecs(eo,group,idx) = real(dval)
                enddo
            enddo
        enddo
        deallocate(strings)
        call fclose(funit)
        contains
        
            subroutine parse_key_int_pair(string, key, val)
                character(len=*), intent(in)  :: string, key
                integer,          intent(out) :: val
                character(len=LENSTR) :: tmp1, tmp2
                tmp1 = trim(string)
                call split(tmp1,' ',tmp2)
                if( trim(tmp2).ne.trim(key) ) THROW_HARD('Unrecognized formatting: '//trim(string)//'; parse_key_int_pair')
                val = str2int(tmp1, iostat)
                call fileiochk('Unrecognized formatting: '//trim(string)//'; parse_key_int_pair', iostat)
            end subroutine parse_key_int_pair
        
            subroutine parse_key_string_pair(string, key, val)
                character(len=*),   intent(in)  :: string, key
                character(len=LENSTR), intent(out) :: val
                character(len=LENSTR) :: tmp2
                val = trim(string)
                call split(val,' ',tmp2)
                if( trim(tmp2).ne.trim(key) ) THROW_HARD('Unrecognized formatting: '//trim(string)//'; parse_key_string_pair')
            end subroutine parse_key_string_pair

    end subroutine read_sigma2_groups

    subroutine read_sigma2_groups_file( fname, group_pspecs, kfromto, ngroups )
        class(string),             intent(in)  :: fname
        real, allocatable,         intent(out) :: group_pspecs(:,:,:)
        integer,                   intent(out) :: kfromto(2), ngroups
        type(euclid_sigma2) :: sigma
        call sigma%init_from_group_header(fname)
        kfromto = sigma%kfromto
        call sigma%read_sigma2_groups(fname, group_pspecs, ngroups)
        call sigma%kill
    end subroutine read_sigma2_groups_file

    ! Destructor

    subroutine kill( self )
        class(euclid_sigma2), intent(inout) :: self
        if( self%exists )then
            if(allocated(self%micinds))       deallocate(self%micinds)
            if(allocated(self%sigma2_groups)) deallocate(self%sigma2_groups)
            if(allocated(self%sigma2_noise))  deallocate(self%sigma2_noise)
            if( allocated(self%sigma2_part) ) deallocate(self%sigma2_part)
            if( allocated(self%diag_ratio) )  deallocate(self%diag_ratio)
            if( allocated(self%diag_v) )      deallocate(self%diag_v)
            self%kfromto     = 0
            self%fromp       = -1
            self%top         = -1
            self%exists      = .false.
        endif
        self%p_ptr => null()
    end subroutine kill

end module simple_euclid_sigma2
