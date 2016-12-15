!==Class simple_commander_checks
!
! This class contains the set of concrete checks commanders of the SIMPLE library (used for checking image stack/volume parameters). 
! This class provides the glue between the reciver (main reciever is simple_exec program) and the abstract action, which is simply 
! execute (defined by the base class: simple_commander_base). Later we can use the composite pattern to create MacroCommanders (or workflows)
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Frederic Bonnet, Cyril Reboul & Hans Elmlund 2016
!
module simple_commander_checks
use simple_defs            ! singleton
use simple_jiffys          ! singleton
use simple_timing          ! singleton
use simple_cuda            ! singleton
use simple_cuda_defs       ! singleton
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
implicit none

public :: check_box_commander
public :: check_nptcls_commander
public :: iminfo_commander
private

type, extends(commander_base) :: check_box_commander
  contains
    procedure :: execute      => exec_check_box
end type check_box_commander
type, extends(commander_base) :: check_nptcls_commander
  contains
    procedure :: execute      => exec_check_nptcls
end type check_nptcls_commander
type, extends(commander_base) :: iminfo_commander
 contains
   procedure :: execute      => exec_iminfo
end type iminfo_commander

contains
    
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
        call simple_end('**** SIMPLE_CHECK_BOX NORMAL STOP ****')
    end subroutine exec_check_box

    subroutine exec_check_nptcls( self, cline )
        class(check_nptcls_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params) :: p
        p = params(cline) ! parameters generated
        call cline%set('nptcls', real(p%nptcls))
        write(*,'(A,1X,I7)') '>>> NPTCLS:', p%nptcls
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK_NPTCLS NORMAL STOP ****')
    end subroutine exec_check_nptcls
    
    subroutine exec_iminfo( self, cline)
        use simple_image,   only: image
        use simple_imgfile, only: imgfile
        class(iminfo_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params)      :: p
        type(image)       :: img
        integer           :: ldim(3), maxim, i, n_nans
        real              :: smpd, sdev, ave, minv, maxv
        p = params(cline) ! constants & derived constants produced
        if( cline%defined('fname') )then
            call find_ldim_nptcls(p%fname, ldim, maxim, doprint=.true.)
        endif
        p%box  = ldim(1)
        p%smpd = smpd          !! UNDEFINED
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
                if( p%vis .eq. 'yes' ) call img%vis
            end do
        endif
        call simple_end('**** SIMPLE_IMINFO NORMAL STOP ****')
    end subroutine exec_iminfo

end module simple_commander_checks
