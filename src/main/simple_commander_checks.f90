! concrete commander: checking routines
module simple_commander_checks
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
implicit none

public :: check_box_commander
public :: check_nptcls_commander
public :: info_image_commander
public :: info_stktab_commander
private

type, extends(commander_base) :: check_box_commander
  contains
    procedure :: execute      => exec_check_box
end type check_box_commander
type, extends(commander_base) :: check_nptcls_commander
  contains
    procedure :: execute      => exec_check_nptcls
end type check_nptcls_commander
type, extends(commander_base) :: info_image_commander
 contains
   procedure :: execute      => exec_info_image
end type info_image_commander
type, extends(commander_base) :: info_stktab_commander
 contains
   procedure :: execute      => exec_info_stktab
end type info_stktab_commander

contains

    !> for checking the image dimensions of MRC and SPIDER stacks and volumes
    subroutine exec_check_box( self, cline )
        class(check_box_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(params) :: p
        if( cline%defined('stk') .or. cline%defined('vol1') )then
            ! all ok
        else
            stop 'Either stack (stk) or volume (vol1) needs to be defined on command line!'
        endif
        p = params(cline) ! parameters generated
        call cline%set('box', real(p%box))
        write(*,'(A,1X,I7)') '>>> BOX:', p%box
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK_BOX NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check_box

    !> for checking the number of images in MRC and SPIDER stacks
    subroutine exec_check_nptcls( self, cline )
        class(check_nptcls_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params) :: p
        p = params(cline) ! parameters generated
        call cline%set('nptcls', real(p%nptcls))
        write(*,'(A,1X,I7)') '>>> NPTCLS:', p%nptcls
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK_NPTCLS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check_nptcls

    !> for printing header information in MRC and SPIDER stacks and volumes
    subroutine exec_info_image( self, cline)
        use simple_image,   only: image
        use simple_imgfile, only: imgfile
        class(info_image_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params)      :: p
        type(image)       :: img
        integer           :: ldim(3), maxim, i, n_nans
        real              :: sdev, ave, minv, maxv
        p = params(cline) ! constants & derived constants produced
        if( cline%defined('fname') )then
            call find_ldim_nptcls(p%fname, ldim, maxim, doprint=.true.)
        endif
        p%box = ldim(1)
        call img%new([ldim(1),ldim(2),1],p%smpd)
        if( p%vis .eq. 'yes' .or. p%stats .ne. 'no' )then
            do i=1,maxim
                call img%read(p%fname, i)
                if( p%stats .ne. 'no' )then
                    call img%cure(maxv, minv, ave, sdev, n_nans)
                    if( p%stats .eq. 'print' .or. n_nans > 0 )then
                        write(*,*) '*********IMAGE********', i, '*******'
                        write(*,*) 'maxv = ',   maxv
                        write(*,*) 'minv = ',   minv
                        write(*,*) 'ave = ',    ave
                        write(*,*) 'sdev = ',   sdev
                        write(*,*) 'n_nans = ', n_nans
                    endif
                endif
                if( p%vis .eq. 'yes' ) call img%vis()
            end do
        endif
        call simple_end('**** SIMPLE_INFO_IMAGE NORMAL STOP ****')
    end subroutine exec_info_image

    !> for printing information about stktab
    subroutine exec_info_stktab( self, cline )
        class(info_stktab_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params) :: p
        p = params(cline) ! constants & derived constants produced
        if( .not. file_exists(p%stktab) )then
            write(*,*) 'file: ', trim(p%stktab), ' not in cwd'
            stop 'commander_checks :: exec_info_stktab'
        endif
        write(*,*) '# micrograps: ', p%nmics
        write(*,*) '# particles : ', p%nptcls
        write(*,*) 'ldim        : ', p%ldim
        write(*,*) 'box size    : ', p%box
        call simple_end('**** SIMPLE_INFO_STKTAB NORMAL STOP ****')
    end subroutine exec_info_stktab

end module simple_commander_checks
