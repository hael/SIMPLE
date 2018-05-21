! concrete commander: time-series analysis
module simple_commander_tseries
include 'simple_lib.f08'
use simple_parameters,     only: parameters
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
implicit none

public :: tseries_extract_commander
public :: tseries_track_commander
public :: tseries_backgr_subtr_commander
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
type, extends(commander_base) :: tseries_backgr_subtr_commander
  contains
    procedure :: execute      => exec_tseries_backgr_subtr
end type tseries_backgr_subtr_commander

contains

    subroutine exec_tseries_extract( self, cline )
        use simple_image,   only: image
        class(tseries_extract_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        character(len=LONGSTRLEN), allocatable :: filenames(:)
        character(len=LONGSTRLEN)              :: outfname
        type(parameters) :: params
        integer          :: ldim(3), nframes, frame_from, frame_to, numlen, cnt
        integer          :: iframe, jframe, nfiles, endit
        type(image)      :: frame_img
        params = parameters(cline)
        if( cline%defined('filetab') )then
            call read_filetable(params%filetab, filenames)
            nfiles = size(filenames)
            numlen = len(int2str(nfiles))
        else
            stop 'need filetab input, listing all the individual frames of&
            &the time series; simple_commander_tseries :: exec_tseries_extract'
        endif
        call find_ldim_nptcls(filenames(1),ldim,nframes)
        if( nframes == 1 .and. ldim(3) == 1 )then
            ! all ok
            call frame_img%new(ldim, params%smpd)
        else
            write(*,*) 'ldim(3): ', ldim(3)
            write(*,*) 'nframes: ', nframes
            stop 'simple_commander_imgproc :: exec_tseries_extract assumes one frame per file'
        endif
        if( cline%defined('nframesgrp') )then
            if( params%nframesgrp < 3 )then
                stop 'nframesgrp integer (nr of frames to average) needs to be >= 3; &
                &simple_commander_imgproc :: exec_tseries_extract'
            endif
        else
            stop 'need nframesgrp integer input = nr of frames to average; &
            &simple_commander_imgproc :: exec_tseries_extract'
        endif
        endit = nfiles - params%nframesgrp + 1
        do iframe=1,endit
            call progress(iframe, endit)
            if( cline%defined('fbody') )then
                outfname = 'tseries_frames'//int2str_pad(iframe,numlen)//params%ext
            else
                outfname = trim(params%fbody)//'tseries_frames'//int2str_pad(iframe,numlen)//params%ext
            endif
            frame_from = iframe
            frame_to   = iframe + params%nframesgrp - 1
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
        use simple_qsys_funs, only: qsys_job_finished
        class(tseries_track_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)  :: params
        type(nrtxtfile)   :: boxfile
        integer           :: ndatlines, j, orig_box, numlen
        real, allocatable :: boxdata(:,:)
        params = parameters(cline)
        numlen = 5 ! default value
        orig_box = params%box
        ! check file inout existence and read filetables
        if( .not. file_exists(params%filetab)  ) stop 'inputted filetab does not exist in cwd'
        if( cline%defined('boxfile') )then
            if( .not. file_exists(params%boxfile)  ) stop 'inputted boxfile does not exist in cwd'
            if( nlines(params%boxfile) > 0 )then
                call boxfile%new(params%boxfile, 1)
                ndatlines = boxfile%get_ndatalines()
                numlen    = len(int2str(ndatlines))
                allocate( boxdata(ndatlines,boxfile%get_nrecs_per_line()), stat=alloc_stat)
                if(alloc_stat.ne.0)call allocchk('In: simple_commander_tseries :: exec_tseries_track boxdata')
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
            boxdata(1,1) = real(params%xcoord)
            boxdata(1,2) = real(params%ycoord)
            ndatlines = 1
        else
            stop 'need either boxfile or xcoord/ycoord to be part of command line&
            &; simple_commander_tseries :: exec_tseries_track'
        endif
        do j=1,ndatlines
            call init_tracker( nint(boxdata(j,1:2)))
            call track_particle
            if( cline%defined('ind') )then
                if( .not. cline%defined('numlen') ) stop 'need numlen to be part of command line if ind is&
                &; simple_commander_tseries :: exec_tseries_track'
                call write_tracked_series(trim(params%fbody)//int2str_pad(params%ind,params%numlen))
            else
                call write_tracked_series(trim(params%fbody)//int2str_pad(j,numlen))
            endif
            call kill_tracker
        end do
        if( params%l_distr_exec ) call qsys_job_finished( 'simple_commander_tseries :: exec_tseries_track')
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_TRACK NORMAL STOP ****')
    end subroutine exec_tseries_track

    subroutine exec_tseries_backgr_subtr( self, cline )
        use simple_image, only: image
        use simple_ctf,   only: ctf
        class(tseries_backgr_subtr_commander), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(image)      :: img_backgr, img_backgr_wctf
        type(ctf)        :: tfun
        logical          :: params_present(4)
        real             :: dfx, dfy, angast
        integer          :: iptcl
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! get background image
        call img_backgr%new([params%box,params%box,1], params%smpd)
        call img_backgr_wctf%new([params%box,params%box,1], params%smpd)
        call img_backgr%read(params%stk_backgr, 1)
        if( cline%defined('deftab') )then
            params_present(1) = build%spproj_field%isthere('kv')
            params_present(2) = build%spproj_field%isthere('cs')
            params_present(3) = build%spproj_field%isthere('fraca')
            params_present(4) = build%spproj_field%isthere('dfx')
            if( all(params_present) )then
                ! alles ok
            else
                if( .not. params_present(1) ) write(*,*) 'ERROR! input deftab lacks kv'
                if( .not. params_present(2) ) write(*,*) 'ERROR! input deftab lacks cs'
                if( .not. params_present(3) ) write(*,*) 'ERROR! input deftab lacks fraca'
                if( .not. params_present(4) ) write(*,*) 'ERROR! input deftab lacks dfx'
                stop
            endif
        endif
        do iptcl=1,params%nptcls
            call progress(iptcl,params%nptcls)
            ! read particle image
            call build%img%read(params%stk, iptcl)
            ! copy background image
            img_backgr_wctf = img_backgr
            ! CTF logics
            if( cline%defined('deftab') )then
                tfun = ctf(params%smpd, build%spproj_field%get(iptcl,'kv'), build%spproj_field%get(iptcl,'cs'), build%spproj_field%get(iptcl,'fraca'))
                dfx = build%spproj_field%get(iptcl, 'dfx')
                dfy    = build%spproj_field%get(iptcl, 'dfy'   )
                angast = build%spproj_field%get(iptcl, 'angast')
                call tfun%apply(build%img, dfx, 'flip', dfy=dfy, angast=angast)
                call tfun%apply(img_backgr_wctf, dfx, 'flip', dfy=dfy, angast=angast)
            endif
            ! fwd ft
            call build%img%fft()
            call img_backgr_wctf%fft()
            ! filter background image
            call img_backgr_wctf%bp(params%hp,params%lp,width=params%width)
            ! subtract background
            call build%img%subtr(img_backgr_wctf)
            ! bwd ft
            call build%img%ifft()
            ! normalise
            call build%img%norm()
            ! output corrected image
            call build%img%write(params%outstk, iptcl)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_BACKGR_SUBTR NORMAL STOP ****')
    end subroutine exec_tseries_backgr_subtr

    subroutine exec_tseries_split( self, cline )
        use simple_oris,     only: oris
        class(tseries_split_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(oris)       :: os
        character(len=:), allocatable :: stkname, oriname
        integer :: i, iptcl, numlen, istart, istop, cnt, cnt2, ntot
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        numlen = len(int2str(params%nptcls))
        ! count the number of chunks
        ntot = 0
        do i=1,params%nptcls,params%stepsz
            istart = i
            istop  = i + params%chunksz - 1
            if( istop > params%nptcls ) exit
            ntot = ntot + 1
        end do
        ! split
        do i=1,params%nptcls,params%stepsz
            istart = i
            istop  = i + params%chunksz - 1
            if( istop > params%nptcls ) exit
            cnt    = cnt + 1
            call progress(cnt,ntot)
            call os%new(params%chunksz)
            call simple_mkdir('tseries_chunk'//int2str_pad(cnt,numlen))
            allocate( stkname, source='./tseries_chunk'//int2str_pad(cnt,numlen)//'/imgs'//params%ext)
            allocate( oriname, source='tseries_chunk'//int2str_pad(cnt,numlen)//'/oris'//trim(TXT_EXT))
            call del_file( stkname )
            call del_file( oriname )
            cnt2 = 0
            do iptcl=istart,istop
                cnt2 = cnt2 + 1
                !o = build%spproj_field%get_ori(iptcl)
                call os%set_ori(cnt2, build%spproj_field%get_ori(iptcl) )
                call build%img%read(params%stk, iptcl)
                call build%img%write(stkname,cnt2)
            end do
            call os%write(oriname)
            call os%kill
            deallocate(stkname, oriname)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_SPLIT NORMAL STOP ****')
    end subroutine exec_tseries_split

end module simple_commander_tseries
