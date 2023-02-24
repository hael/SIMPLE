module simple_euclid_sigma2
include 'simple_lib.f08'
use simple_parameters,       only: params_glob
use simple_cartft_corrcalc,  only: cartft_corrcalc, cftcc_glob
use simple_polarft_corrcalc, only: polarft_corrcalc, pftcc_glob
use simple_sigma2_binfile,   only: sigma2_binfile
use simple_starfile_wrappers
implicit none

public :: euclid_sigma2, eucl_sigma2_glob, write_groups_starfile
public :: split_group_sigma2
public :: fill_sigma2_before_nyq, test_unit
private
#include "simple_local_flags.inc"

character(len=STDLEN), parameter :: SIGMA2_ONE_GROUP = 'sigma2_group_'

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
    procedure, private :: calc_sigma2_1, calc_sigma2_2
    generic            :: calc_sigma2 => calc_sigma2_1, calc_sigma2_2
    procedure, private :: update_sigma2_1, update_sigma2_2
    generic            :: update_sigma2 => update_sigma2_1, update_sigma2_2
    procedure          :: write_sigma2
    procedure, private :: read_groups_starfile
    ! destructor
    procedure          :: kill_ptclsigma2
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
        if( params_glob%l_sigma_glob )then
            call self%read_groups_starfile( params_glob%which_iter, self%sigma2_groups, ngroups )
            if( ngroups /= 1 ) THROW_HARD('ngroups must be 1 when global sigma is estimated (params_glob%l_sigma_glob == .true.)')
            ! copy global sigma to particles
            do iptcl = params_glob%fromp, params_glob%top
                eo = nint(os%get(iptcl, 'eo')) ! 0/1
                self%sigma2_noise(:,iptcl) = self%sigma2_groups(eo+1,1,:)
            end do
        else
            call self%read_groups_starfile( params_glob%which_iter, self%sigma2_groups, ngroups )
            ! copy group sigmas to particles
            do iptcl = params_glob%fromp, params_glob%top
                igroup  = nint(os%get(iptcl, 'stkind'))
                eo      = nint(os%get(iptcl, 'eo'    )) ! 0/1
                self%sigma2_noise(:,iptcl) = self%sigma2_groups(eo+1,igroup,:)
            end do
        endif
    end subroutine read_groups

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
        class(euclid_sigma2), intent(inout) :: self
        type(sigma2_binfile)                :: binfile
        call binfile%new_from_file(self%binfname)
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
            starfile_fname = trim(SIGMA2_GROUP_FBODY) // trim(int2str(iter)) // '.star'
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
    end subroutine read_groups_starfile

    ! Public modifiers

    ! Updates the lowest resolution info of the file with most frequencies with the other & overwrites it
    subroutine fill_sigma2_before_nyq( fname1, fname2 )
        character(len=*), intent(in)  :: fname1, fname2
        type(euclid_sigma2)           :: sigma2_1, sigma2_2
        integer                       :: ngroups1, ngroups2, k0
        call sigma2_1%init_from_group_header(fname1)
        call sigma2_2%init_from_group_header(fname2)
        call sigma2_1%read_groups_starfile(0, sigma2_1%sigma2_groups, ngroups1, fname=fname1)
        call sigma2_2%read_groups_starfile(0, sigma2_2%sigma2_groups, ngroups2, fname=fname2)
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
            print *,trim(fname1)
            call write_groups_starfile( fname1, sigma2_2%sigma2_groups, ngroups1 )
        endif
        call sigma2_1%kill
        call sigma2_2%kill
    end subroutine fill_sigma2_before_nyq

    subroutine split_group_sigma2( iter )
        integer,          intent(in)  :: iter
        type(euclid_sigma2)           :: euclidsigma2
        character(len=:), allocatable :: fname
        real,             allocatable :: sigma2_group(:,:,:)
        integer                       :: igroup, ngroups
        fname  = trim(SIGMA2_GROUP_FBODY) // trim(int2str(iter)) // '.star'
        call euclidsigma2%init_from_group_header(fname)
        call euclidsigma2%read_groups_starfile(iter, euclidsigma2%sigma2_groups, ngroups)
        do igroup = 1,ngroups
            fname  = trim(SIGMA2_ONE_GROUP)//int2str(igroup)//'.star'
            allocate(sigma2_group(2,1,euclidsigma2%kfromto(1):euclidsigma2%kfromto(2)),&
            &source=euclidsigma2%sigma2_groups(:,igroup:igroup,:))
            call write_groups_starfile( fname, sigma2_group, 1 )
            deallocate(sigma2_group)
        enddo
        call euclidsigma2%kill
    end subroutine split_group_sigma2

    ! Destructor

    subroutine kill_ptclsigma2( self )
        class(euclid_sigma2), intent(inout) :: self
        if( allocated(self%sigma2_noise) ) deallocate(self%sigma2_noise)
    end subroutine kill_ptclsigma2

    subroutine kill( self )
        class(euclid_sigma2), intent(inout) :: self
        if( self%exists )then
            call self%kill_ptclsigma2
            if(allocated(self%pinds)  )       deallocate(self%pinds)
            if(allocated(self%micinds))       deallocate(self%micinds)
            if(allocated(self%sigma2_groups)) deallocate(self%sigma2_groups)
            if(allocated(self%sigma2_noise))  deallocate(self%sigma2_noise)
            self%kfromto     = 0
            self%fromp       = -1
            self%top         = -1
            self%exists      = .false.
            eucl_sigma2_glob => null()
        endif
    end subroutine kill

    subroutine test_unit
        integer, parameter :: ngroups    = 19
        integer, parameter :: kfromto(2) = [1,64]
        integer, parameter :: iter       = 7
        real,    parameter :: scale      = 0.3
        type(euclid_sigma2)           :: euclidsigma2
        character(len=STDLEN)         :: fname, fname1, fname2
        real,             allocatable :: sigma2(:,:,:)
        integer :: igroup, ng
        logical :: l_err
        l_err = .false.
        ! testing bookkeeping
        allocate(sigma2(2,ngroups,kfromto(1):kfromto(2)),source=1.0)
        fname  = trim(SIGMA2_GROUP_FBODY) // trim(int2str(iter)) // '.star'
        call write_groups_starfile( fname, sigma2, ngroups )
        call split_group_sigma2( iter )
        do igroup = 1,ngroups
            fname = trim(SIGMA2_ONE_GROUP)//int2str(igroup)//'.star'
            if(.not.file_exists(fname))then
                l_err = .true.
                THROW_WARN('File does not exists for group: '//int2str(igroup))
            else
                call euclidsigma2%init_from_group_header(fname)
                call euclidsigma2%read_groups_starfile(iter, euclidsigma2%sigma2_groups, ng, fname=fname )
                if( ng /= 1 )then
                    l_err = .true.
                    THROW_WARN('Erroneous group number for group: '//int2str(igroup))
                endif
                if( any(abs(euclidsigma2%sigma2_groups-scale) >  0.001) )then
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
        fname1 = trim(SIGMA2_GROUP_FBODY) // '1.star'
        call write_groups_starfile( fname1, sigma2, ngroups )
        deallocate(sigma2)
        allocate(sigma2(2,ngroups,kfromto(1):2*kfromto(2)),source=2.0)
        fname2 = trim(SIGMA2_GROUP_FBODY) // '2.star'
        call write_groups_starfile( fname2, sigma2, ngroups )
        deallocate(sigma2)
        call fill_sigma2_before_nyq(fname1, fname2)
    end subroutine test_unit

end module simple_euclid_sigma2
