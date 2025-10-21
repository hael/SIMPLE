! particle picker iterator
module simple_picker_iter
include 'simple_lib.f08'
use simple_picker
use simple_parameters
use simple_picker_utils
use simple_cmdline, only: cmdline
use simple_image,   only: image
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

    subroutine iterate( self, cline, smpd, moviename_intg, dir_out, boxfile, nptcls_out, moldiam_opt )
        class(picker_iter),        intent(inout) :: self
        class(cmdline),            intent(in)    :: cline
        real,                      intent(in)    :: smpd
        character(len=*),          intent(in)    :: moviename_intg
        character(len=*),          intent(in)    :: dir_out
        character(len=LONGSTRLEN), intent(out)   :: boxfile
        integer,                   intent(out)   :: nptcls_out
        real, optional,            intent(out)   :: moldiam_opt
        logical                                  :: l_append
        l_append = .false.
        if( .not. file_exists(moviename_intg) ) write(logfhandle,*) 'inputted micrograph does not exist: ', trim(adjustl(moviename_intg))
        write(logfhandle,'(a,1x,a)') '>>> PICKING MICROGRAPH:', trim(adjustl(moviename_intg))
        if(params_glob%append == 'yes') l_append = .true.
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
                if( cline%defined('pickrefs') )then 
                    call self%read_pickrefs(params_glob%pickrefs)
                    if( cline%defined('nboxes_max') )then
                        call exec_gaupick(moviename_intg, boxfile, smpd, nptcls_out, self%pickrefs, dir_out=dir_out, nboxes_max=params_glob%nboxes_max)
                    else
                        call exec_gaupick(moviename_intg, boxfile, smpd, nptcls_out, self%pickrefs, dir_out=dir_out)
                    endif
                else if( cline%defined('moldiam') .or. cline%defined('multi_moldiams')  )then
                    ! at least moldiam is required
                    call exec_gaupick(moviename_intg, boxfile, smpd, nptcls_out, moldiam_opt=moldiam_opt, dir_out=dir_out, append=l_append)
                else
                    THROW_HARD('New picker requires 2D references (pickrefs) or moldiam')
                endif
            case('seg')
                if( .not. cline%defined('lp') )then
                    THROW_HARD('Segmentation-based picker requires lp (low-pass limit) for filtering')
                endif
                if( cline%defined('moldiam') )then
                    call exec_segpick(moviename_intg, boxfile, nptcls_out, dir_out=dir_out, moldiam=params_glob%moldiam)
                elseif( cline%defined('winsz') )then
                    call exec_segpick(moviename_intg, boxfile, nptcls_out, dir_out=dir_out, winsz=int(params_glob%winsz))
                else
                    call exec_segpick(moviename_intg, boxfile, nptcls_out, dir_out=dir_out)
                endif
            case('segdiam')
                if( cline%defined('moldiam_max') )then
                    call exec_segdiampick(moviename_intg, boxfile, nptcls_out,params_glob%moldiam_max, dir_out=dir_out)
                endif
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
