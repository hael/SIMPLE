! concrete commander: checking routines
module simple_commander_checks
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: parameters
implicit none

public :: check_box_commander
public :: check_nptcls_commander
public :: info_image_commander
public :: info_stktab_commander
private
#include "simple_local_flags.inc"

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
        type(parameters) :: params
        if( cline%defined('stk') .or. cline%defined('vol1') )then
            ! all ok
        else
            THROW_HARD('Either stack (stk) or volume (vol1) needs to be defined on command line!')
        endif
        call params%new(cline)
        call cline%set('box', real(params%box))
        write(logfhandle,'(A,1X,I7)') '>>> BOX:', params%box
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK_BOX NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check_box

    !> for checking the number of images in MRC and SPIDER stacks
    subroutine exec_check_nptcls( self, cline )
        class(check_nptcls_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters) :: params
        call params%new(cline)
        call cline%set('nptcls', real(params%nptcls))
        write(logfhandle,'(A,1X,I7)') '>>> NPTCLS:', params%nptcls
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK_NPTCLS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check_nptcls

    !> for printing header information in MRC and SPIDER stacks and volumes
    subroutine exec_info_image( self, cline)
        use simple_image,   only: image
        class(info_image_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters) :: params
        type(image)      :: img
        integer          :: ldim(3), maxim, i, n_nans
        real             :: sdev, ave, minv, maxv
        call params%new(cline)
        if( cline%defined('fname') )then
            call find_ldim_nptcls(params%fname, ldim, maxim, doprint=.true.)
        endif
        if( params%vis .eq. 'yes' .or. params%stats .ne. 'no' )then
            params%box = ldim(1)
            call img%new([ldim(1),ldim(2),1],params%smpd)
            do i=1,maxim
                call img%read(params%fname, i)
                if( params%stats .ne. 'no' )then
                    call img%cure(maxv, minv, ave, sdev, n_nans)
                    if( params%stats .eq. 'print' .or. n_nans > 0 )then
                        write(logfhandle,*) '*********IMAGE********', i, '*******'
                        write(logfhandle,*) 'maxv = ',   maxv
                        write(logfhandle,*) 'minv = ',   minv
                        write(logfhandle,*) 'ave = ',    ave
                        write(logfhandle,*) 'sdev = ',   sdev
                        write(logfhandle,*) 'n_nans = ', n_nans
                    endif
                endif
                if( params%vis .eq. 'yes' ) call img%vis()
            end do
        endif
        call simple_end('**** SIMPLE_INFO_IMAGE NORMAL STOP ****')
    end subroutine exec_info_image

    !> for printing information about stktab
    subroutine exec_info_stktab( self, cline )
        class(info_stktab_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        character(LONGSTRLEN), allocatable :: stkfnames(:)
        type(parameters) :: params
        integer          :: ldim(3), box, istk, np_stk
        call params%new(cline)
        if( .not. file_exists(params%stktab) )then
            THROW_HARD('file: '//trim(params%stktab)//' not in cwd; exec_info_stktab')
        endif
        call read_filetable(params%stktab, stkfnames)
        params%nmics  = size(stkfnames)
        if( params%nmics < 1 ) THROW_HARD('Nonsensical # of stacks: '//int2str(params%nmics))
        params%nptcls = 0
        do istk = 1, params%nmics
            ! get dimension of stack and # ptcls
            call find_ldim_nptcls(trim(stkfnames(istk)), ldim, np_stk)
            if( istk == 1 )then
                box         = ldim(1)
                params%ldim = ldim
            else
                if( box /= ldim(1) )then
                    THROW_HARD('Box size of stk #'//int2str(istk)//' is '//int2str(ldim(1))//' which is inconsistent with box='//int2str(box)//' of stk #1')
                endif
            endif
            ! update nptcls
            params%nptcls = params%nptcls + np_stk
        end do
        params%box = box
        write(logfhandle,*) '# micrograps   : ', params%nmics
        write(logfhandle,*) '# particles    : ', params%nptcls
        write(logfhandle,*) 'ldim of stk #1 : ', params%ldim
        write(logfhandle,*) 'box size       : ', params%box
        call simple_end('**** SIMPLE_INFO_STKTAB NORMAL STOP ****')
    end subroutine exec_info_stktab

end module simple_commander_checks
