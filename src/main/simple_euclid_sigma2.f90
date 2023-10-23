module simple_euclid_sigma2
include 'simple_lib.f08'
use simple_parameters,       only: params_glob
use simple_cartft_corrcalc,  only: cartft_corrcalc, cftcc_glob
use simple_polarft_corrcalc, only: polarft_corrcalc, pftcc_glob
use simple_sigma2_binfile,   only: sigma2_binfile
use simple_starfile_wrappers
implicit none

public :: euclid_sigma2, eucl_sigma2_glob, write_groups_starfile
public :: split_sigma2_into_groups, consolidate_sigma2_groups
public :: sigma2_star_from_iter, fill_sigma2_before_nyq, test_unit
private
#include "simple_local_flags.inc"

integer, parameter :: LENSTR = 48

type euclid_sigma2
    private
    real,    allocatable, public  :: sigma2_noise(:,:)      !< the sigmas for alignment & reconstruction (from groups)
    real,    allocatable          :: sigma2_part(:,:)       !< the actual sigmas per particle (this part only)
    real,    allocatable          :: sigma2_groups(:,:,:)   !< sigmas for groups
    integer, allocatable          :: pinds(:)
    integer, allocatable          :: micinds(:)
    integer                       :: fromp
    integer                       :: top
    integer                       :: kfromto(2) = 0
    character(len=:), allocatable :: binfname
    logical                       :: exists     = .false.
contains
    ! constructor
    procedure          :: new
    procedure, private :: init_from_group_header
    ! utils
    procedure          :: write_info
    ! I/O
    procedure          :: read_part
    procedure          :: read_groups
    procedure          :: allocate_ptcls
    procedure, private :: calc_sigma2_1, calc_sigma2_2
    generic            :: calc_sigma2 => calc_sigma2_1, calc_sigma2_2
    procedure, private :: update_sigma2_1, update_sigma2_2
    generic            :: update_sigma2 => update_sigma2_1, update_sigma2_2
    procedure          :: write_sigma2
    procedure, private :: read_groups_starfile, read_sigma2_groups
    ! destructor
    procedure          :: kill
end type euclid_sigma2

class(euclid_sigma2), pointer :: eucl_sigma2_glob => null()

contains

    subroutine new( self, binfname, box )
        ! read individual sigmas from binary file, to be modified at the end of the iteration
        ! read group sigmas from starfile, to be used for alignment and volume reconstruction
        ! set up fields for fast access to sigmas
        class(euclid_sigma2), target, intent(inout) :: self
        character(len=*),             intent(in)    :: binfname
        integer,                      intent(in)    :: box
        call self%kill
        self%kfromto = [1, fdim(box)-1]
        allocate( self%sigma2_noise(self%kfromto(1):self%kfromto(2),params_glob%fromp:params_glob%top),&
                  self%pinds(params_glob%fromp:params_glob%top) )
        if( associated(pftcc_glob) )then
            call pftcc_glob%assign_sigma2_noise(self%sigma2_noise)
            call pftcc_glob%assign_pinds(self%pinds)
        endif
        if( associated(cftcc_glob) )then
            call cftcc_glob%assign_sigma2_noise(self%sigma2_noise)
            call cftcc_glob%assign_pinds(self%pinds)
        endif
        self%binfname     =  trim(binfname)
        self%fromp        =  params_glob%fromp
        self%top          =  params_glob%top
        self%sigma2_noise =  0.
        self%exists       =  .true.
        eucl_sigma2_glob  => self
    end subroutine new

    !>  This is a minimal constructor to allow I/O of groups
    subroutine init_from_group_header( self, fname )
        class(euclid_sigma2), target, intent(inout) :: self
        character(len=*),             intent(in)    :: fname
        type(str4arr), allocatable :: names(:)
        type(starfile_table_type)  :: table
        integer                    :: kfromto(2), ngroups
        logical                    :: l
        if (.not. file_exists(fname)) then
            THROW_HARD('euclid_sigma2_starfile: read_groups_pspecs; file does not exists: ' // trim(fname))
        end if
        call starfile_table__new(table)
        call starfile_table__getnames(table, trim(fname), names)
        call starfile_table__read( table, trim(fname), names(1)%str )
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

    ! I/O

    subroutine read_part( self, os, ptcl_mask )
        class(euclid_sigma2), intent(inout) :: self
        class(oris),          intent(inout) :: os
        logical,              intent(in)    :: ptcl_mask(params_glob%fromp:params_glob%top)
        type(sigma2_binfile) :: binfile
        call binfile%new_from_file(self%binfname)
        call binfile%read(self%sigma2_part)
    end subroutine read_part

    subroutine read_groups( self, os, ptcl_mask )
        class(euclid_sigma2), intent(inout) :: self
        class(oris),          intent(inout) :: os
        logical,              intent(in)    :: ptcl_mask(params_glob%fromp:params_glob%top)
        integer                             :: iptcl, igroup, ngroups, eo
        if( associated(pftcc_glob) ) call pftcc_glob%assign_pinds(self%pinds)
        if( associated(cftcc_glob) ) call cftcc_glob%assign_pinds(self%pinds)
        call self%read_sigma2_groups( params_glob%which_iter, self%sigma2_groups, ngroups )
        if( params_glob%l_sigma_glob )then
            if( ngroups /= 1 ) THROW_HARD('ngroups must be 1 when global sigma is estimated (params_glob%l_sigma_glob == .true.)')
            ! copy global sigma to particles
            !$omp parallel do default(shared) private(iptcl,eo) proc_bind(close) schedule(static)
            do iptcl = params_glob%fromp, params_glob%top
                if(os%get_state(iptcl) == 0 ) cycle
                eo = nint(os%get(iptcl, 'eo')) ! 0/1
                self%sigma2_noise(:,iptcl) = self%sigma2_groups(eo+1,1,:)
            end do
            !$omp end parallel do
        else
            ! copy group sigmas to particles
            !$omp parallel do default(shared) private(iptcl,eo,igroup) proc_bind(close) schedule(static)
            do iptcl = params_glob%fromp, params_glob%top
                if(os%get_state(iptcl) == 0 ) cycle
                igroup  = nint(os%get(iptcl, 'stkind'))
                eo      = nint(os%get(iptcl, 'eo'    )) ! 0/1
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
    subroutine calc_sigma2_1( self, pftcc, iptcl, o, refkind )
        class(euclid_sigma2),    intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(ori),              intent(in)    :: o
        character(len=*),        intent(in)    :: refkind ! 'proj' or 'class'
        integer :: iref, irot
        real    :: sigma_contrib(params_glob%kfromto(1):params_glob%kfromto(2))
        real    :: shvec(2)
        if ( o%isstatezero() ) return
        shvec = o%get_2Dshift()
        iref  = nint(o%get(trim(refkind)))
        irot  = pftcc_glob%get_roind(360. - o%e3get())
        call pftcc%gencorr_sigma_contrib(iref, iptcl, shvec, irot, sigma_contrib)
        self%sigma2_part(params_glob%kfromto(1):params_glob%kfromto(2),iptcl) = sigma_contrib
    end subroutine calc_sigma2_1

    !>  Calculates and updates sigma2 within search resolution range
    subroutine calc_sigma2_2( self, cftcc, iptcl, o, refkind )
        class(euclid_sigma2),   intent(inout) :: self
        class(cartft_corrcalc), intent(inout) :: cftcc
        integer,                intent(in)    :: iptcl
        class(ori),             intent(in)    :: o
        character(len=*),       intent(in)    :: refkind ! 'proj' or 'class'
        real :: sigma_contrib(params_glob%kfromto(1):params_glob%kfromto(2))
        real :: shvec(2)
        if ( o%isstatezero() ) return
        shvec = o%get_2Dshift()
        call cftcc%calc_sigma_contrib(iptcl, o, shvec, sigma_contrib)
        self%sigma2_part(params_glob%kfromto(1):params_glob%kfromto(2),iptcl) = sigma_contrib
    end subroutine calc_sigma2_2

    !>  Calculates and updates sigma2 within search resolution range
    subroutine update_sigma2_1( self, pftcc, iptcl, o, refkind )
        class(euclid_sigma2),    intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(ori),              intent(in)    :: o
        character(len=*),        intent(in)    :: refkind ! 'proj' or 'class'
        integer :: iref, irot
        real    :: shvec(2)
        if ( o%isstatezero() ) return
        shvec = o%get_2Dshift()
        iref  = nint(o%get(trim(refkind)))
        irot  = pftcc_glob%get_roind(360. - o%e3get())
        call pftcc%update_sigma(iref, iptcl, shvec, irot)
    end subroutine update_sigma2_1

    !>  Calculates and updates sigma2 within search resolution range
    subroutine update_sigma2_2( self, cftcc, iptcl, o, refkind )
        class(euclid_sigma2),   intent(inout) :: self
        class(cartft_corrcalc), intent(inout) :: cftcc
        integer,                intent(in)    :: iptcl
        class(ori),             intent(in)    :: o
        character(len=*),       intent(in)    :: refkind ! 'proj' or 'class'
        real :: shvec(2)
        if ( o%isstatezero() ) return
        shvec = o%get_2Dshift()
        call cftcc%update_sigma(iptcl, o, shvec)
    end subroutine update_sigma2_2

    subroutine write_sigma2( self )
        class(euclid_sigma2), intent(in) :: self
        type(sigma2_binfile) :: binfile
        if( file_exists(self%binfname) )then
            call binfile%new_from_file(self%binfname)
        else
            call binfile%new(self%binfname, self%fromp, self%top, self%kfromto)
        endif
        call binfile%write(self%sigma2_part)
        call binfile%kill
    end subroutine write_sigma2

    subroutine write_groups_starfile( fname, group_pspecs, ngroups )
        character(len=*),  intent(in) :: fname
        real, allocatable, intent(in) :: group_pspecs(:,:,:)
        integer,           intent(in) :: ngroups
        character(len=:), allocatable :: stmp
        integer                       :: kfromto(2), eo, igroup, idx
        type(starfile_table_type)     :: ostar
        call starfile_table__new(ostar)
        call starfile_table__open_ofile(ostar, trim(fname))
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
                call starfile_table__setComment(ostar, stmp // ', group ' // trim(int2str(igroup)) )
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

    ! Is deprecated, kept for compatibility. cf read_sigma2_groups below
    subroutine read_groups_starfile( self, iter, group_pspecs, ngroups, fname )
        class(euclid_sigma2),          intent(inout) :: self
        integer,                       intent(in)    :: iter
        real,             allocatable, intent(out)   :: group_pspecs(:,:,:)
        integer,                       intent(out)   :: ngroups
        character(len=*), optional,    intent(in)    :: fname
        type(str4arr),    allocatable :: names(:)
        type(starfile_table_type)     :: istarfile
        character(len=:), allocatable :: starfile_fname
        character                     :: eo_char
        real(dp)                      :: val
        integer(C_long)               :: num_objs, object_id
        integer                       :: kfromto(2), stat, spec_idx, eo, igroup, idx
        logical                       :: l
        if( present(fname) )then
            starfile_fname = trim(fname)
        else
            starfile_fname = sigma2_star_from_iter(iter)
        endif
        if (.not. file_exists(starfile_fname)) then
            THROW_HARD('euclid_sigma2: read_groups_starfile; file does not exists: ' // starfile_fname)
        end if
        call starfile_table__new(istarfile)
        ! read header
        call starfile_table__getnames(istarfile, starfile_fname, names)
        call starfile_table__read( istarfile, starfile_fname, names(1)%str )
        l = starfile_table__getValue_int(istarfile, EMDL_MLMODEL_NR_GROUPS, ngroups)
        l = starfile_table__getValue_int(istarfile, EMDL_SPECTRAL_IDX, kfromto(1))
        l = starfile_table__getValue_int(istarfile, EMDL_SPECTRAL_IDX2, kfromto(2))
        if( any(kfromto-self%kfromto < 0) )then
            print *,kfromto,self%kfromto
            THROW_HARD('Incorrect resolution range: read_groups_starfile')
        endif
        ! read values
        allocate(group_pspecs(2,ngroups,self%kfromto(1):self%kfromto(2)))
        call starfile_table__getnames(istarfile, starfile_fname, names)
        do idx = 1, size(names)
            if( len(names(idx)%str) < len('1_group_')+1 )cycle
            if( names(idx)%str(2:8) .ne. '_group_' ) cycle
            eo_char = names(idx)%str(1:1)
            if ((eo_char .ne. '1').and.(eo_char .ne. '2')) cycle
            eo = 1
            if (eo_char == '2') eo = 2
            call str2int( names(idx)%str(9:len_trim(names(idx)%str)), stat, igroup )
            if( stat > 0 ) cycle
            if( (igroup < 1).or.(igroup>ngroups) ) cycle
            call starfile_table__read( istarfile, starfile_fname, names(idx)%str )
            object_id = starfile_table__firstobject(istarfile)
            num_objs  = starfile_table__numberofobjects(istarfile)
            do while( (object_id < num_objs) .and. (object_id >= 0) )
                l = starfile_table__getValue_int(istarfile, EMDL_SPECTRAL_IDX, spec_idx)
                if( l ) then
                    l = starfile_table__getValue_double(istarfile, EMDL_MLMODEL_SIGMA2_NOISE, val)
                    if( l ) then
                        if( (spec_idx >= self%kfromto(1)).and.(spec_idx <= self%kfromto(2)) ) then
                            group_pspecs(eo,igroup,spec_idx) = real(val)
                        end if
                    end if
                end if
                object_id = starfile_table__nextobject(istarfile)
            end do
        end do
        call starfile_table__delete(istarfile)
        deallocate(names)
    end subroutine read_groups_starfile

    ! Is a much faster hard-coded fortran substitute to %read_groups_starfile
    ! that employs the starfile library
    subroutine read_sigma2_groups( self, iter, pspecs, ngroups, filename )
        class(euclid_sigma2),          intent(inout) :: self
        integer,                       intent(in)    :: iter
        real,             allocatable, intent(out)   :: pspecs(:,:,:)
        integer,                       intent(out)   :: ngroups
        character(len=*), optional,    intent(in)    :: filename
        character(len=LENSTR), allocatable :: strings(:)
        character(len=:),      allocatable :: fname
        character(len=LENSTR) :: line, string
        real(dp) :: dval
        integer  :: kfromto(2), i, l, funit, iostat, group, eo, idx, igroup, ieo
        if( present(filename) )then
            fname = trim(filename)
        else
            fname = sigma2_star_from_iter(iter)
        endif
        if(.not.file_exists(fname))then
            THROW_HARD('File: '//trim(fname)//' Does not exists; read_sigma2_groups')
        endif
        call fopen(funit, trim(fname), action='READ',status='OLD', form='FORMATTED', iostat=iostat)
        call fileiochk('read_sigma2_groups: '//trim(fname), iostat)
        ! read header
        ! # of groups
        read(funit,fmt='(A)') line
        read(funit,fmt='(A)') line
        if( trim(line).ne.'data_general' )then
            THROW_HARD('Unrecognized formatting: '//trim(fname)//'; read_sigma2_groups')
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
            THROW_HARD('Incorrect resolution range: read_groups_starfile')
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
                call str2int(line(6:6+i-2), iostat, ieo)
                call fileiochk('Invalid formatting: '//trim(line)//'; read_sigma2_groups', iostat)
                if( ieo /= eo ) THROW_HARD('Invalid formatting: '//trim(line)//'; read_sigma2_groups')
                if( line(6+i-1:6+i+5).ne.'_group_')THROW_HARD('Invalid formatting: '//trim(line)//'; read_sigma2_groups')
                i = index(line, '_', back=.true.)
                call str2int(line(i+1:l), iostat, igroup)
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
                call str2int(tmp1, iostat, val)
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

    ! Public modifiers

    ! Updates the lowest resolution info of the file with most frequencies with the other & overwrites it
    subroutine fill_sigma2_before_nyq( fname1, fname2 )
        character(len=*), intent(in)  :: fname1, fname2
        type(euclid_sigma2)           :: sigma2_1, sigma2_2
        integer                       :: ngroups1, ngroups2, k0
        call sigma2_1%init_from_group_header(fname1)
        call sigma2_2%init_from_group_header(fname2)
        call sigma2_1%read_sigma2_groups(0, sigma2_1%sigma2_groups, ngroups1, filename=fname1)
        call sigma2_2%read_sigma2_groups(0, sigma2_2%sigma2_groups, ngroups2, filename=fname2)
        if( ngroups1 /= ngroups2 ) THROW_HARD('Inconsistent dimensions; fill_sigma2_beyond_nyq')
        if( sigma2_1%kfromto(1) /= sigma2_2%kfromto(1) ) THROW_HARD('Inconsistent fourier dimensions 1; fill_sigma2_beyond_nyq')
        k0 = sigma2_1%kfromto(1)
        if( sigma2_1%kfromto(2) > sigma2_2%kfromto(2) )then
            sigma2_1%sigma2_groups(:,:,k0:sigma2_2%kfromto(2)) = sigma2_2%sigma2_groups
            call write_groups_starfile( fname2, sigma2_1%sigma2_groups, ngroups1 )
        else if(sigma2_1%kfromto(2) == sigma2_2%kfromto(2))then
            ! nothing to do?
        else
            sigma2_2%sigma2_groups(:,:,k0:sigma2_1%kfromto(2)) = sigma2_1%sigma2_groups
            call write_groups_starfile( fname1, sigma2_2%sigma2_groups, ngroups1 )
        endif
        call sigma2_1%kill
        call sigma2_2%kill
    end subroutine fill_sigma2_before_nyq

    !> Split a sigma2 doc into individual docs
    subroutine split_sigma2_into_groups( fname, fnames )
        character(len=*),                        intent(in) :: fname
        character(len=LONGSTRLEN), allocatable,  intent(in) :: fnames(:)
        type(euclid_sigma2)           :: euclidsigma2
        real,             allocatable :: sigma2_group(:,:,:)
        integer                       :: igroup, ngroups
        call euclidsigma2%init_from_group_header(fname)
        call euclidsigma2%read_sigma2_groups(0, euclidsigma2%sigma2_groups, ngroups, filename=fname)
        if( ngroups /= size(fnames) ) THROW_HARD('Inconsistent number of groups & stacks! split_group_sigma2')
        do igroup = 1,ngroups
            allocate(sigma2_group(2,1,euclidsigma2%kfromto(1):euclidsigma2%kfromto(2)),&
            &source=euclidsigma2%sigma2_groups(:,igroup:igroup,:))
            call write_groups_starfile( fnames(igroup), sigma2_group, 1 )
            deallocate(sigma2_group)
        enddo
        call euclidsigma2%kill
    end subroutine split_sigma2_into_groups

    ! the reverse of split_sigma2_into_groups
    subroutine consolidate_sigma2_groups( fname, fnames )
        character(len=*),                        intent(in) :: fname
        character(len=LONGSTRLEN), allocatable,  intent(in) :: fnames(:)
        type(euclid_sigma2)           :: euclidsigma2
        real,             allocatable :: sigma2_group(:,:,:)
        integer                       :: igroup, ngroups, n
        call euclidsigma2%init_from_group_header(fnames(1))
        call euclidsigma2%read_sigma2_groups(1, euclidsigma2%sigma2_groups, n, filename=fnames(1))
        ngroups = size(fnames)
        allocate(sigma2_group(2,ngroups,euclidsigma2%kfromto(1):euclidsigma2%kfromto(2)),source=0.0)
        sigma2_group(:,1,:) = euclidsigma2%sigma2_groups(:,1,:)
        call euclidsigma2%kill
        do igroup = 2,ngroups
            call euclidsigma2%read_sigma2_groups(1, euclidsigma2%sigma2_groups, n, filename=fnames(igroup))
            sigma2_group(:,igroup,:) = euclidsigma2%sigma2_groups(:,1,:)
            call euclidsigma2%kill
        enddo
        call write_groups_starfile(fname, sigma2_group, ngroups)
        deallocate(sigma2_group)
    end subroutine consolidate_sigma2_groups

    function sigma2_star_from_iter( iter )
        integer, intent(in) :: iter
        character(len=:), allocatable :: sigma2_star_from_iter
        sigma2_star_from_iter = trim(SIGMA2_GROUP_FBODY) // trim(int2str(iter)) // '.star'
    end function sigma2_star_from_iter

    ! Destructor

    subroutine kill( self )
        class(euclid_sigma2), intent(inout) :: self
        if( self%exists )then
            if(allocated(self%pinds)  )       deallocate(self%pinds)
            if(allocated(self%micinds))       deallocate(self%micinds)
            if(allocated(self%sigma2_groups)) deallocate(self%sigma2_groups)
            if(allocated(self%sigma2_noise))  deallocate(self%sigma2_noise)
            if( allocated(self%sigma2_part) ) deallocate(self%sigma2_part)
            self%kfromto     = 0
            self%fromp       = -1
            self%top         = -1
            self%exists      = .false.
            eucl_sigma2_glob => null()
        endif
    end subroutine kill

    subroutine test_unit
        integer, parameter :: ngroups    = 2000
        integer, parameter :: kfromto(2) = [1,128]
        integer, parameter :: iter       = 7
        real,    parameter :: scale      = 0.3
        type(euclid_sigma2)           :: euclidsigma2
        character(len=STDLEN)         :: fname, fname1, fname2
        character(len=LONGSTRLEN), allocatable :: fnames(:)
        real,                      allocatable :: sigma2(:,:,:)
        integer :: igroup, ng
        logical :: l_err
        l_err = .false.
        ! testing bookkeeping
        allocate(fnames(ngroups))
        do igroup = 1,ngroups
            fnames(igroup) = 'test_'//int2str(igroup)//'.star'
        enddo
        allocate(sigma2(2,ngroups,kfromto(1):kfromto(2)),source=1.0)
        ! call seed_rnd()
        ! call random_number(sigma2)
        ! sigma2 = sigma2 / 1.e6
        fname  = sigma2_star_from_iter(iter)
        call write_groups_starfile( fname, sigma2, ngroups )
        call split_sigma2_into_groups( fname, fnames )
        fname  = sigma2_star_from_iter(666)
        call consolidate_sigma2_groups( fname, fnames)
        call split_sigma2_into_groups( fname, fnames)
        do igroup = 1,ngroups
            fname = fnames(igroup)
            if(.not.file_exists(fname))then
                l_err = .true.
                THROW_WARN('File does not exists for group: '//int2str(igroup))
            else
                call euclidsigma2%init_from_group_header(fname)
                call euclidsigma2%read_sigma2_groups(iter, euclidsigma2%sigma2_groups, ng, filename=fname )
                if( ng /= 1 )then
                    l_err = .true.
                    THROW_WARN('Erroneous group number for group: '//int2str(igroup))
                endif
                if( any(abs(euclidsigma2%sigma2_groups-scale) >  0.000001) )then
                    l_err = .true.
                    THROW_WARN('Scaling failed for group: '//int2str(igroup))
                endif
                if( size(euclidsigma2%sigma2_groups,dim=3) /= (kfromto(2)-kfromto(1)+1) )then
                    l_err = .true.
                    THROW_WARN('Incorrect number of frequencies for group: '//int2str(igroup))
                endif
                call euclidsigma2%kill
            endif
        enddo
        if( l_err )then
            write(*,'(A)')'>>> EUCLID_SIGMA2 UNIT TEST 1 FAILED'
        else
            write(*,'(A)')'>>> EUCLID_SIGMA2 UNIT TEST 1 PASSED'
        endif
        deallocate(sigma2)
        allocate(sigma2(2,ngroups,kfromto(1):kfromto(2)),source=1.0)
        fname1 = sigma2_star_from_iter(1)
        call write_groups_starfile( fname1, sigma2, ngroups )
        deallocate(sigma2)
        allocate(sigma2(2,ngroups,kfromto(1):2*kfromto(2)),source=2.0)
        fname2 = sigma2_star_from_iter(2)
        call write_groups_starfile( fname2, sigma2, ngroups )
        deallocate(sigma2)
        call fill_sigma2_before_nyq(fname1, fname2)
    end subroutine test_unit

end module simple_euclid_sigma2
