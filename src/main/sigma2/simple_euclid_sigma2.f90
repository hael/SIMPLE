!@descr: the abstract data type for sigma2 used when objfun=euclid
module simple_euclid_sigma2
use, intrinsic :: iso_fortran_env, only: int64, real32
use simple_core_module_api
use simple_polarft_calc,   only: polarft_calc
use simple_cartft_calc,    only: cartft_calc
use simple_parameters,     only: parameters
use simple_syslib,         only: simple_atomic_replace
use simple_sigma2_state,   only: sigma2_state_candidate_path, sigma2_state_range_path, sigma2_state_next_generation
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_read_header, &
    &sigma2_state_read_groups, sigma2_state_read_particles, sigma2_state_write_local_range, &
    &sigma2_state_read_local_range
use simple_starfile_wrappers
implicit none
private

#include "simple_local_flags.inc"

public :: euclid_sigma2, sigma2_group_iter
! grouped-STAR I/O for the explicit sigma2_convert boundary only
public :: write_groups_starfile, read_sigma2_groups_file

integer, parameter :: LENSTR = 48

type euclid_sigma2
    private
    class(parameters),    pointer :: p_ptr => null()
    real,    allocatable, public  :: sigma2_noise(:,:)      !< the sigmas for alignment & reconstruction (from groups)
    real,    allocatable          :: sigma2_part(:,:)       !< the actual sigmas per particle (this part only)
    real,    allocatable          :: sigma2_groups(:,:,:)   !< sigmas for groups
    integer, allocatable          :: pinds(:)
    integer, allocatable          :: micinds(:)
    integer                       :: fromp
    integer                       :: top
    integer                       :: kfromto(2) = 0
    type(string)                  :: binfname
    logical                       :: exists     = .false.
    logical                       :: replaces_range = .false. !< started from the open transaction's range (read_pending_range)
contains
    ! constructor
    procedure, private :: new_pftc
    procedure, private :: new_cart
    generic            :: new => new_pftc, new_cart
    procedure, private :: init_from_group_header
    ! utils
    procedure          :: write_info
    procedure          :: get_kfromto
    procedure          :: set_kfromto
    ! I/O
    procedure          :: read_part
    procedure          :: read_pending_range
    procedure          :: read_groups
    procedure          :: allocate_ptcls
    procedure, private :: calc_sigma2_pftc
    procedure, private :: calc_sigma2_cart
    generic            :: calc_sigma2 => calc_sigma2_pftc, calc_sigma2_cart
    procedure, private :: store_contribution
    procedure          :: get_sigma2_part
    procedure          :: write_sigma2
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

    !> Polar pass: the noise table is shared with the polar calculator.
    subroutine new_pftc( self, params, pftc, binfname, box )
        ! read individual sigmas from binary file, to be modified at the end of the iteration
        ! read group sigmas from starfile, to be used for alignment and volume reconstruction
        ! set up fields for fast access to sigmas
        class(euclid_sigma2), target, intent(inout) :: self
        class(parameters),    target, intent(in)    :: params
        class(polarft_calc),          intent(inout) :: pftc
        class(string),                intent(in)    :: binfname
        integer,                      intent(in)    :: box
        call self%new_cart(params, binfname, box)
        call pftc%assign_sigma2_noise(self%sigma2_noise)
    end subroutine new_pftc

    !> The noise table of a pass: a Cartesian pass, which constructs no polar calculator, uses
    !! it alone; the polar constructor shares it with the calculator.
    subroutine new_cart( self, params, binfname, box )
        class(euclid_sigma2), target, intent(inout) :: self
        class(parameters),    target, intent(in)    :: params
        class(string),                intent(in)    :: binfname
        integer,                      intent(in)    :: box

        call self%kill
        self%p_ptr => params
        self%kfromto = [1, fdim(box)-1]
        allocate(self%sigma2_noise(self%kfromto(1):self%kfromto(2), &
            &self%p_ptr%fromp:self%p_ptr%top), source=0.)
        self%binfname = binfname
        self%fromp = self%p_ptr%fromp
        self%top = self%p_ptr%top
        self%exists = .true.
    end subroutine new_cart

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

    !> Replace this part's rows by the range an earlier pass of the open transaction wrote for it
    !! (the discrete pass a polish follows), so a particle the later pass leaves untouched keeps
    !! the residual of its current pose rather than the committed generation's.
    subroutine read_pending_range( self )
        class(euclid_sigma2), intent(inout) :: self
        type(string) :: range_path
        real(real32), allocatable :: spectra(:,:)
        integer(int64) :: next_gen, generation, layout_digest
        integer :: first_row, last_row, kfrom, kto, status
        character(len=STDLEN) :: message
        if( .not. allocated(self%sigma2_part) ) THROW_HARD('particle sigma2 storage is not allocated')
        call sigma2_state_next_generation(self%binfname%to_char(), next_gen, status, message)
        if( status /= 0 ) THROW_HARD(trim(message))
        range_path = sigma2_state_range_path(self%binfname%to_char(), next_gen, self%p_ptr%part, self%p_ptr%numlen)
        call sigma2_state_read_local_range(range_path%to_char(), generation, layout_digest, first_row, last_row, &
            &kfrom, kto, spectra, status, message)
        if( status /= 0 ) THROW_HARD('no sigma2 range of the open transaction for this part: '//trim(message))
        if( generation /= next_gen ) THROW_HARD('pending sigma2 range belongs to another transaction')
        if( first_row /= self%fromp .or. last_row /= self%top ) THROW_HARD('pending sigma2 range does not cover this part')
        if( kfrom /= self%kfromto(1) .or. kto /= self%kfromto(2) ) THROW_HARD('pending sigma2 range has incompatible shell bounds')
        self%sigma2_part(:,:) = real(spectra)
        self%replaces_range   = .true.
        deallocate(spectra)
        call range_path%kill
    end subroutine read_pending_range

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
    !> The sigma2 contribution of particle iptcl at its assigned polar orientation o.
    subroutine calc_sigma2_pftc( self, pftc, iptcl, o, refkind )
        class(euclid_sigma2), intent(inout) :: self
        class(polarft_calc),  intent(inout) :: pftc
        integer,              intent(in)    :: iptcl
        class(ori),           intent(in)    :: o
        character(len=*),     intent(in)    :: refkind ! 'proj' or 'class'
        integer :: iref, irot, kfromto(2)
        real, allocatable :: sigma_contrib(:)
        real    :: shvec(2)
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
        call pftc%gen_sigma_contrib(iref, iptcl, shvec, irot, sigma_contrib)
        call self%store_contribution(iptcl, sigma_contrib)
        deallocate(sigma_contrib)
    end subroutine calc_sigma2_pftc

    !> The sigma2 contribution of particle iptcl at the committed Cartesian pose, from particle
    !! slot islot of the Cartesian calculator (prepared for objfun=euclid; sigma2 belongs to the
    !! Euclidean objective only, C5): the mean squared residual over two per sample per shell
    !! of the slot's shell range (cartft_calc%sigma_contribution). rotmat and shift (cropped
    !! pixels, the model shift relative to the slot's observation) are the committed pose;
    !! state and iseven select the reference. THROW_HARD when the slot's shell range is not the
    !! pass's band.
    subroutine calc_sigma2_cart( self, cftc, islot, iptcl, state, iseven, rotmat, shift )
        class(euclid_sigma2), intent(inout) :: self
        class(cartft_calc),   intent(in)    :: cftc
        integer,              intent(in)    :: islot, iptcl, state
        logical,              intent(in)    :: iseven
        real(dp),             intent(in)    :: rotmat(3,3), shift(2)
        real, allocatable :: sigma_contrib(:)
        if( .not.associated(self%p_ptr) )then
            THROW_HARD('euclid_sigma2: params pointer is not set')
        endif
        if( .not. cftc%ptcl_is_valid(islot) ) return
        ! residual only: the diagnostic outputs of sigma_contribution are test instruments
        call cftc%sigma_contribution(state, iseven, islot, rotmat, shift, sigma_contrib)
        if( lbound(sigma_contrib,1) /= self%p_ptr%kfromto(1) .or. ubound(sigma_contrib,1) /= self%p_ptr%kfromto(2) ) &
            &THROW_HARD('Cartesian sigma contribution does not span the band of the pass')
        call self%store_contribution(iptcl, sigma_contrib)
    end subroutine calc_sigma2_cart

    !> The particle's sigma2 contribution of this part over the pass's band (as calc_sigma2 stored it).
    function get_sigma2_part( self, iptcl ) result( sigma2 )
        class(euclid_sigma2), intent(in) :: self
        integer,              intent(in) :: iptcl
        real, allocatable :: sigma2(:)
        if( .not. allocated(self%sigma2_part) ) THROW_HARD('particle sigma2 storage is not allocated')
        if( iptcl < lbound(self%sigma2_part,2) .or. iptcl > ubound(self%sigma2_part,2) ) &
            &THROW_HARD('particle index is outside sigma2 storage')
        sigma2 = self%sigma2_part(self%p_ptr%kfromto(1):self%p_ptr%kfromto(2),iptcl)
    end function get_sigma2_part

    !> Store one particle's per-shell contribution over the pass's band.
    subroutine store_contribution( self, iptcl, sigma_contrib )
        class(euclid_sigma2), intent(inout) :: self
        integer,              intent(in)    :: iptcl
        real,                 intent(in)    :: sigma_contrib(:)
        integer :: kfromto(2)
        if( .not. allocated(self%sigma2_part) ) THROW_HARD('particle sigma2 storage is not allocated')
        if( iptcl < lbound(self%sigma2_part,2) .or. iptcl > ubound(self%sigma2_part,2) ) &
            &THROW_HARD('particle index is outside sigma2 storage')
        kfromto = self%p_ptr%kfromto
        self%sigma2_part(kfromto(1):kfromto(2),iptcl) = sigma_contrib
    end subroutine store_contribution

    subroutine write_sigma2( self )
        class(euclid_sigma2), intent(inout) :: self
        type(sigma2_state_header) :: header
        type(string) :: candidate_path, range_path, staged_path
        integer(int64) :: next_gen
        integer :: status
        character(len=STDLEN) :: message
        ! transaction-scoped names: the candidate and this range carry the
        ! generation the master's update will commit
        call sigma2_state_next_generation(self%binfname%to_char(), next_gen, status, message)
        if( status /= 0 ) THROW_HARD(trim(message))
        candidate_path = sigma2_state_candidate_path(self%binfname%to_char(), next_gen)
        range_path = sigma2_state_range_path(self%binfname%to_char(), next_gen, self%p_ptr%part, self%p_ptr%numlen)
        call sigma2_state_read_header(candidate_path%to_char(), header, status, message)
        if( status /= 0 ) THROW_HARD(trim(message))
        if( header%generation /= next_gen ) THROW_HARD('canonical sigma2 candidate belongs to another transaction')
        if( self%replaces_range )then
            ! the range this pass started from stays until its successor is complete
            staged_path = range_path%to_char()//'.tmp'
            call del_file(staged_path)
            call sigma2_state_write_local_range(staged_path%to_char(), header%generation, header%layout_digest, &
                &self%fromp, real(self%sigma2_part,real32), self%kfromto(1), self%kfromto(2), status, message)
            if( status /= 0 ) THROW_HARD(trim(message))
            call simple_atomic_replace(staged_path, range_path, status)
            if( status /= 0 ) THROW_HARD('cannot replace the pending sigma2 range')
            call staged_path%kill
        else
            call sigma2_state_write_local_range(range_path%to_char(), header%generation, header%layout_digest, &
                &self%fromp, real(self%sigma2_part,real32), self%kfromto(1), self%kfromto(2), status, message)
            if( status /= 0 ) THROW_HARD(trim(message))
        endif
        call candidate_path%kill
        call range_path%kill
    end subroutine write_sigma2

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
            self%kfromto     = 0
            self%fromp       = -1
            self%top         = -1
            self%exists      = .false.
        endif
        self%replaces_range = .false.
        self%p_ptr => null()
    end subroutine kill

end module simple_euclid_sigma2
