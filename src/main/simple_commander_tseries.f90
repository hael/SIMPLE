! concrete commander: time-series analysis
module simple_commander_tseries
use simple_defs
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_strings,        only: int2str, int2str_pad
use simple_filehandling    ! use all in there
use simple_jiffys          ! use all in there
implicit none

public :: tseries_extract_commander
public :: tseries_track_commander
public :: tseries_split_commander
private
#include "simple_local_flags.inc"
type, extends(commander_base) :: tseries_extract_commander
  contains
    procedure :: execute      => exec_tseries_extract
end type tseries_extract_commander
type, extends(commander_base) :: tseries_track_commander
  contains
    procedure :: execute      => exec_tseries_track
end type tseries_track_commander
type, extends(commander_base) :: tseries_split_commander
  contains
    procedure :: execute      => exec_tseries_split
end type tseries_split_commander

contains

    subroutine exec_tseries_extract( self, cline )
        use simple_image, only: image
        class(tseries_extract_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(params) :: p
        character(len=STDLEN), allocatable :: filenames(:)
        character(len=STDLEN)              :: outfname
        integer      :: ldim(3), nframes, frame_from, frame_to, numlen, cnt
        integer      :: iframe, jframe, nfiles
        type(image)  :: frame_img
        p = params(cline) ! parameters generated
        if( cline%defined('filetab') )then
            call read_filetable(p%filetab, filenames)
            nfiles = size(filenames)
            numlen = len(int2str(nfiles))
        else
            stop 'need filetab input, listing all the individual frames of&
            &the time series; simple_commander_tseries :: exec_tseries_extract'
        endif
        call find_ldim_nptcls(filenames(1),ldim,nframes)
        if( nframes == 1 .and. ldim(3) == 1 )then
            ! all ok
            call frame_img%new(ldim, p%smpd)
        else
            write(*,*) 'ldim(3): ', ldim(3)
            write(*,*) 'nframes: ', nframes
            stop 'simple_commander_imgproc :: exec_tseries_extract assumes one frame per file' 
        endif
        if( cline%defined('nframesgrp') )then
            if( p%nframesgrp < 3 )then
                stop 'nframesgrp integer (nr of frames to average) needs to be >= 3; &
                &simple_commander_imgproc :: exec_tseries_extract'
            endif
        else
            stop 'need nframesgrp integer input = nr of frames to average; &
            &simple_commander_imgproc :: exec_tseries_extract'
        endif
        do iframe=1,nfiles - p%nframesgrp + 1
            if( cline%defined('fbody') )then
                outfname = 'tseries_frames'//int2str_pad(iframe,numlen)//p%ext
            else
                outfname = trim(p%fbody)//'tseries_frames'//int2str_pad(iframe,numlen)//p%ext
            endif
            frame_from = iframe
            frame_to   = iframe + p%nframesgrp - 1
            cnt = 0
            do jframe=frame_from,frame_to
                cnt = cnt + 1
                call frame_img%read(filenames(jframe),1)
                call frame_img%write(outfname,cnt)
            end do
        end do
        call frame_img%kill
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_EXTRACT NORMAL STOP ****')
    end subroutine exec_tseries_extract

    subroutine exec_tseries_track( self, cline )
        use simple_tseries_tracker
        use simple_nrtxtfile, only: nrtxtfile
        use simple_qsys_funs, only: qsys_job_finished
        class(tseries_track_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(params)      :: p
        type(nrtxtfile)   :: boxfile
        integer           :: ndatlines, alloc_stat, j, orig_box, numlen
        real, allocatable :: boxdata(:,:)
        p = params(cline, checkdistr=.false.) ! parameters generated
        numlen = 5 ! default value
        orig_box = p%box
        ! check file inout existence and read filetables
        if( .not. file_exists(p%filetab)  ) stop 'inputted filetab does not exist in cwd'
        if( cline%defined('boxfile') )then
            if( .not. file_exists(p%boxfile)  ) stop 'inputted boxfile does not exist in cwd'
            if( nlines(p%boxfile) > 0 )then
                call boxfile%new(p%boxfile, 1)
                ndatlines = boxfile%get_ndatalines()
                numlen    = len(int2str(ndatlines))
                allocate( boxdata(ndatlines,boxfile%get_nrecs_per_line()), stat=alloc_stat)
                call alloc_err('In: simple_commander_tseries :: exec_tseries_track', alloc_stat)
                do j=1,ndatlines
                    call boxfile%readNextDataLine(boxdata(j,:))
                    orig_box = nint(boxdata(j,3))
                    if( nint(boxdata(j,3)) /= nint(boxdata(j,4)) )then
                        stop 'Only square windows are currently allowed!'
                    endif
                end do
            else
                stop 'inputted boxfile is empty; simple_commander_tseries :: exec_tseries_track'
            endif
        else if( cline%defined('xcoord') .and. cline%defined('ycoord') )then
            if( .not. cline%defined('box') ) stop 'need box to be part of command linefor this mode of&
            &execution; simple_commander_tseries :: exec_tseries_track'
            allocate( boxdata(1,2) )
            boxdata(1,1) = real(p%xcoord)
            boxdata(1,2) = real(p%ycoord)
            ndatlines = 1
        else
            stop 'need either boxfile or xcoord/ycoord to be part of command line&
            &; simple_commander_tseries :: exec_tseries_track'
        endif
        do j=1,ndatlines
            call init_tracker(p, nint(boxdata(j,1:2)))
            call track_particle
            if( cline%defined('ind') )then
                if( .not. cline%defined('numlen') ) stop 'need numlen to be part of command line if ind is&
                &; simple_commander_tseries :: exec_tseries_track'
                call write_tracked_series(trim(p%fbody)//int2str_pad(p%ind,p%numlen))
            else
                call write_tracked_series(trim(p%fbody)//int2str_pad(j,numlen))
            endif
            call kill_tracker
        end do
        if( p%l_distr_exec ) call qsys_job_finished(p, 'simple_commander_tseries :: exec_tseries_track')
    end subroutine exec_tseries_track

    subroutine exec_tseries_split( self, cline )
        use simple_oris,     only: oris
        use simple_ori,      only: ori
        use simple_syscalls, only: exec_cmdline
        class(tseries_split_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        type(ori)    :: o
        type(oris)   :: os
        character(len=:), allocatable :: stkname, oriname
        integer :: i, iptcl, numlen, istart, istop, cnt, cnt2, ntot
        p = params(cline) ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        numlen = len(int2str(p%nptcls))
        ! count the number of chunks
        ntot = 0
        do i=1,p%nptcls,p%stepsz
            istart = i
            istop  = i + p%chunksz - 1
            if( istop > p%nptcls ) exit
            ntot = ntot + 1
        end do
        ! split
        do i=1,p%nptcls,p%stepsz
            istart = i
            istop  = i + p%chunksz - 1
            if( istop > p%nptcls ) exit
            cnt    = cnt + 1
            call progress(cnt,ntot)
            call os%new(p%chunksz)
            call exec_cmdline('mkdir -p '//'tseries_chunk'//int2str_pad(cnt,numlen))
            allocate( stkname, source='./tseries_chunk'//int2str_pad(cnt,numlen)//'/imgs'//p%ext )
            allocate( oriname, source='tseries_chunk'//int2str_pad(cnt,numlen)//'/oris'//'.txt' )
            call del_file( stkname )
            call del_file( oriname )
            cnt2 = 0
            do iptcl=istart,istop
                cnt2 = cnt2 + 1
                o = b%a%get_ori(iptcl)
                call os%set_ori(cnt2, o)
                call b%img%read(p%stk, iptcl)
                call b%img%write(stkname,cnt2)
            end do
            call os%write(oriname)
            call os%kill
            deallocate(stkname, oriname)
        end do
    end subroutine exec_tseries_split

end module simple_commander_tseries
