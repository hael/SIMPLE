!==Class simple_commander_mask
!
! This class contains the set of concrete mask commanders of the SIMPLE library. This class provides the glue between the reciver 
! (main reciever is simple_exec program) and the abstract action, which is simply execute (defined by the base class: simple_commander_base). 
! Later we can use the composite pattern to create MacroCommanders (or workflows)
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Cyril Reboul & Hans Elmlund 2016
!
module simple_commander_mask
use simple_defs
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_filehandling    ! use all in there
use simple_jiffys          ! use all in there
implicit none

public :: mask_commander
public :: automask2D_commander
public :: automask3D_commander
private

type, extends(commander_base) :: mask_commander
 contains
   procedure :: execute      => exec_mask
end type mask_commander
type, extends(commander_base) :: automask2D_commander
  contains
    procedure :: execute      => exec_automask2D
end type automask2D_commander
type, extends(commander_base) :: automask3D_commander
  contains
    procedure :: execute      => exec_automask3D
end type automask3D_commander

contains
  
    subroutine exec_mask( self, cline )
        use simple_procimgfile,   only: mask_imgfile
        class(mask_commander), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        type(build)                :: b
        type(params)               :: p
        type(automask2D_commander) :: automask2D
        type(automask3D_commander) :: automask3D
        logical                    :: here
        p = params(cline)                   ! parameters generated
        p%boxmatch = p%box                  !!!!!!!!!!!!!!!!!! 4 NOW
        call b%build_general_tbox(p, cline) ! general objects built
        if( cline%defined('stk') .and. cline%defined('vol1') )stop 'Cannot operate on images AND volume at once'
        if( p%automsk.eq.'yes' .and..not.cline%defined('mw')  )stop 'Missing mw argument for automasking'
        if( p%msktype.ne.'soft' .and. p%msktype.ne.'hard' .and. p%msktype.ne.'cavg' )stop 'Invalid mask type'
        call b%vol%new([p%box,p%box,p%box], p%smpd) ! reallocate vol (boxmatch issue)
        if( cline%defined('stk') )then
            ! 2D
            if( p%automsk.eq.'yes' )then
                ! auto
                if( .not. cline%defined('amsklp') )call cline%set('amsklp', 25.)
                if( .not. cline%defined('edge')   )call cline%set('edge', 20.)
                call exec_automask2D( automask2D, cline )
            else if( cline%defined('msk') .or. cline%defined('inner') )then
                ! spherical
                if( cline%defined('inner') )then
                    if( cline%defined('width') )then
                        call mask_imgfile(p%stk, p%outstk, p%msk, inner=p%inner, width=p%width, which=p%msktype)
                    else
                        call mask_imgfile(p%stk, p%outstk, p%msk, inner=p%inner, which=p%msktype)
                    endif
                else
                    call mask_imgfile(p%stk, p%outstk, p%msk, which=p%msktype)
                endif
            else
                stop 'Nothing to do!'
            endif
        else if( cline%defined('vol1') )then
            ! 3D
            here = .false.
            inquire(FILE=p%vols(1), EXIST=here)
            if( .not.here )stop 'Cannot find input volume'
            call b%vol%read(p%vols(1))
            if( p%automsk.eq.'yes' )then
                ! auto
                call exec_automask3D( automask3D, cline )
            else if( cline%defined('msk') )then
                ! spherical
                if( cline%defined('inner') )then
                    if( cline%defined('width') )then
                        call b%vol%mask(p%msk, p%msktype, inner=p%inner, width=p%width)
                    else
                        call b%vol%mask(p%msk, p%msktype, inner=p%inner)
                    endif
                else
                    call b%vol%mask(p%msk, p%msktype)
                endif
                if( p%outvol .ne. '' )call b%vol%write(p%outvol, del_if_exists=.true.)
            else
                stop 'Nothing to do!'
            endif
        else
            stop 'No input images(s) or volume provided'
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_MASK NORMAL STOP ****')
    end subroutine exec_mask

    subroutine exec_automask2D( self, cline )
        use simple_masker, only: automask2D
        class(automask2D_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer      :: iptcl
        p = params(cline)                   ! parameters generated
        p%boxmatch = p%box                   !!!!!!!!!!!!!!!!!! 4 NOW
        call b%build_general_tbox(p, cline) ! general objects built
        write(*,'(A,F8.2,A)') '>>> AUTOMASK LOW-PASS:',        p%amsklp, ' ANGSTROMS'
        write(*,'(A,I3,A)')   '>>> AUTOMASK SOFT EDGE WIDTH:', p%edge,   ' PIXELS'
        p%outstk = add2fbody(p%stk, p%ext, 'msk') 
        do iptcl=1,p%nptcls
            call b%img%read(p%stk, iptcl)
            call automask2D(b%img, p, b%img_msk)
            call b%img_msk%write('automasks2D'//p%ext, iptcl)
            call b%img%write(p%outstk, iptcl)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_AUTOMASK2D NORMAL STOP ****')
    end subroutine exec_automask2D
    
    subroutine exec_automask3D( self, cline )
        use simple_masker,  only: automask
        use simple_strings, only: int2str_pad
        class(automask3D_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer      :: istate
        p = params(cline)                   ! parameters generated
        p%boxmatch = p%box                  !!!!!!!!!!!!!!!!!! 4 NOW
        call b%build_general_tbox(p, cline) ! general objects built
        write(*,'(A,F14.1,A)') '>>> AUTOMASK LOW-PASS:',            p%amsklp,  ' ANGSTROMS'
        write(*,'(A,I7,A)')    '>>> AUTOMASK SOFT EDGE WIDTH:',     p%edge,    ' PIXEL(S)'
        write(*,'(A,I3,A)')    '>>> AUTOMASK BINARY LAYERS WIDTH:', p%binwidth,' PIXEL(S)'
        do istate=1,p%nstates
            p%masks(istate)    = 'automask_state'//int2str_pad(istate,2)//p%ext
            p%vols_msk(istate) = add2fbody(p%vols(istate), p%ext, 'msk')
            call b%vol%read(p%vols(istate))
            call automask(b, p, cline, b%vol,  b%mskvol, p%vols_msk(istate), p%masks(istate))
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_AUTOMASK3D NORMAL STOP ****')
    end subroutine exec_automask3D

end module simple_commander_mask
