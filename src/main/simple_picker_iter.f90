! particle picker iterator
module simple_picker_iter
include 'simple_lib.f08'
use simple_picker
use simple_parameters
use simple_picker_utils, only: picker_utils
use simple_cmdline,      only: cmdline
use simple_image,        only: image
implicit none

public :: picker_iter

#include "simple_local_flags.inc"

private

type :: picker_iter
    type(image), allocatable :: pickrefs(:)
    logical                  :: l_pickrefs_exist = .false.
  contains
    procedure          :: iterate
    procedure, private :: read_pickrefs
    procedure          :: kill
end type picker_iter

contains

    subroutine iterate( self, cline, smpd, moviename_intg, boxfile, nptcls_out, dir_out )
        class(picker_iter),        intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        real,                      intent(in)    :: smpd
        character(len=*),          intent(in)    :: moviename_intg
        character(len=LONGSTRLEN), intent(out)   :: boxfile
        integer,                   intent(out)   :: nptcls_out
        character(len=*),          intent(in)    :: dir_out
        type(picker_utils) :: putils
        integer            :: ldim(3), ifoo
        if( .not. file_exists(moviename_intg) ) write(logfhandle,*) 'inputted micrograph does not exist: ', trim(adjustl(moviename_intg))
        write(logfhandle,'(a,1x,a)') '>>> PICKING MICROGRAPH:', trim(adjustl(moviename_intg))
        select case(trim(params_glob%picker))
            case('old')
                if( .not. cline%defined('pickrefs') ) THROW_HARD('Old picker requires pickrefs (2D picking references) input')
                call self%read_pickrefs(params_glob%pickrefs)
                if( cline%defined('thres') )then
                    call init_picker(moviename_intg, self%pickrefs, smpd, lp_in=params_glob%lp,&
                    &distthr_in=params_glob%thres, ndev_in=params_glob%ndev, dir_out=dir_out)
                else
                    call init_picker(moviename_intg, self%pickrefs, smpd, lp_in=params_glob%lp, &
                                                  &ndev_in=params_glob%ndev, dir_out=dir_out)
                endif
                call exec_picker(boxfile, nptcls_out)
                call kill_picker
            case('new')
                ! read micrograph
                call find_ldim_nptcls(moviename_intg, ldim, ifoo)
                call putils%new(moviename_intg, params_glob%pcontrast, smpd, params_glob%moldiam)
                if( cline%defined('pickrefs') )then
                    if( .not. cline%defined('mskdiam') ) THROW_HARD('New picker requires mask diameter (in A) in conjunction with pickrefs')  
                    call self%read_pickrefs(params_glob%pickrefs)
                    call putils%set_refs(self%pickrefs, params_glob%mskdiam)
                else if( cline%defined('moldiam') )then
                    ! at least moldiam is required
                else
                    THROW_HARD('New picker requires 2D references (pickrefs) or moldiam')
                endif
                call putils%exec_picker(boxfile, nptcls_out, dir_out=dir_out)
                call putils%kill
        end select
    end subroutine iterate

    subroutine read_pickrefs( self, pickrefs_fname )
        class(picker_iter), intent(inout) :: self
        character(len=*),   intent(in)    :: pickrefs_fname
        real    :: smpd
        integer :: ldim(3), nrefs, iref
        if( self%l_pickrefs_exist ) return
        if( file_exists(pickrefs_fname) )then
            call find_ldim_nptcls(trim(pickrefs_fname), ldim, nrefs, smpd=smpd)
            ldim(3) = 1
            allocate(self%pickrefs(nrefs))
            do iref = 1,nrefs
                call self%pickrefs(iref)%new(ldim, smpd)
                call self%pickrefs(iref)%read(trim(pickrefs_fname), iref)
            end do
        else
            THROW_HARD('file '//trim(pickrefs_fname)//' does not exist')
        endif
        self%l_pickrefs_exist = .true.
    end subroutine read_pickrefs

    subroutine kill( self )
        class(picker_iter), intent(inout) :: self
        integer :: iref
        if( self%l_pickrefs_exist )then
            do iref = 1,size(self%pickrefs)
                call self%pickrefs(iref)%kill
            enddo
            deallocate(self%pickrefs)
        endif
        self%l_pickrefs_exist = .false.
    end subroutine kill

end module simple_picker_iter
