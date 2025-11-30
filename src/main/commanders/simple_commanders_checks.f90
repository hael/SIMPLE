! concrete commander: checking routines
module simple_commanders_checks
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: parameters
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_check_box
  contains
    procedure :: execute      => exec_check_box
end type commander_check_box

type, extends(commander_base) :: commander_check_nptcls
  contains
    procedure :: execute      => exec_check_nptcls
end type commander_check_nptcls

type, extends(commander_base) :: commander_check_stoch_update
  contains
    procedure :: execute      => exec_check_stoch_update
end type commander_check_stoch_update

type, extends(commander_base) :: commander_check_update_frac
  contains
    procedure :: execute      => exec_check_update_frac
end type commander_check_update_frac

type, extends(commander_base) :: commander_info_image
 contains
   procedure :: execute      => exec_info_image
end type commander_info_image

type, extends(commander_base) :: commander_info_stktab
 contains
   procedure :: execute      => exec_info_stktab
end type commander_info_stktab

contains

    !> for checking the image dimensions of MRC and SPIDER stacks and volumes
    subroutine exec_check_box( self, cline )
        class(commander_check_box), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters) :: params
        if( cline%defined('stk') .or. cline%defined('vol1') )then
            ! all ok
        else
            THROW_HARD('Either stack (stk) or volume (vol1) needs to be defined on command line!')
        endif
        call params%new(cline)
        call cline%set('box', params%box)
        write(logfhandle,'(A,1X,I7)') '>>> BOX:', params%box
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK_BOX NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check_box

    !> for checking the number of images in MRC and SPIDER stacks
    subroutine exec_check_nptcls( self, cline )
        class(commander_check_nptcls), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters) :: params
        call params%new(cline)
        call cline%set('nptcls', params%nptcls)
        write(logfhandle,'(A,1X,I7)') '>>> NPTCLS:', params%nptcls
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK_NPTCLS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check_nptcls

    subroutine exec_check_stoch_update( self, cline )
        use simple_decay_funs
        class(commander_check_stoch_update), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters) :: params
        integer          :: i, nsampl_fromto(2), nsampl
        call params%new(cline)
        nsampl_fromto = calc_nsampl_fromto(params%nptcls, NSAMPLE_MINMAX_DEFAULT)
        write(logfhandle,'(A,1X,I7,1X,I7)') '>>> NSAMPLE, FROM/TO:', nsampl_fromto(1), nsampl_fromto(2)
        do i = 1, params%maxits
            nsampl = inv_nsampl_decay(i, params%maxits, params%nptcls, NSAMPLE_MINMAX_DEFAULT)
            write(logfhandle,'(A,1X,I7,1X,I7)') 'ITER/NSAMPLE:', i, nsampl
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK_STOCH_UPDATE NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check_stoch_update

    subroutine exec_check_update_frac( self, cline )
        use simple_decay_funs
        class(commander_check_update_frac), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters) :: params
        integer          :: minmax(2), nsampl, i
        real             :: update_frac
        call params%new(cline)
        minmax = NSAMPLE_MINMAX_DEFAULT
        if( cline%defined('nsample_start') .and. cline%defined('nsample_stop') )then
            do i = 1, params%maxits
                update_frac = calc_update_frac_dyn(params%nptcls, 1, [params%nsample_start,params%nsample_stop], i, params%maxits)
                nsampl      = update_frac * real(params%nptcls)
                write(logfhandle,'(A,1X,I7,1X,A,1X,I7)') 'ITER:', i, 'NSAMPLE:', nsampl
            enddo
        else
            if( cline%defined('nsample_max') ) minmax(2) = params%nsample_max
            update_frac = calc_update_frac(params%nptcls, 1, minmax)
            nsampl      = update_frac * real(params%nptcls)
            write(logfhandle,'(A,1X,I7)') 'NSAMPLE:', nsampl
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK_UPDATE_FRAC NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check_update_frac

    !> for printing header information in MRC and SPIDER stacks and volumes
    subroutine exec_info_image( self, cline)
        use simple_image,   only: image
        class(commander_info_image), intent(inout) :: self
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
        class(commander_info_stktab), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(string), allocatable :: stkfnames(:)
        type(parameters) :: params
        integer          :: ldim(3), box, istk, np_stk
        call params%new(cline)
        if( .not. file_exists(params%stktab) )then
            THROW_HARD('file: '//params%stktab%to_char()//' not in cwd; exec_info_stktab')
        endif
        call read_filetable(params%stktab, stkfnames)
        params%nmics  = size(stkfnames)
        if( params%nmics < 1 ) THROW_HARD('Nonsensical # of stacks: '//int2str(params%nmics))
        params%nptcls = 0
        do istk = 1, params%nmics
            ! get dimension of stack and # ptcls
            call find_ldim_nptcls(stkfnames(istk), ldim, np_stk)
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

end module simple_commanders_checks
