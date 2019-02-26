! concrete commander: time-series analysis
module simple_commander_tseries
include 'simple_lib.f08'
use simple_parameters,     only: parameters
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
implicit none

public :: tseries_import_commander
public :: tseries_track_commander
public :: tseries_backgr_subtr_commander
public :: tseries_split_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: tseries_import_commander
  contains
    procedure :: execute      => exec_tseries_import
end type tseries_import_commander
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

    !> for creating a new project
    subroutine exec_tseries_import( self, cline )
        use simple_sp_project,        only: sp_project
        use simple_image,             only: image
        class(tseries_import_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)            :: params
        type(sp_project)            :: spproj
        type(image)                 :: frame_img
        type(ctfparams)             :: ctfvars
        character(len=LONGSTRLEN), allocatable :: filenames(:), movfnames(:)
        character(len=LONGSTRLEN)              :: outfname
        integer          :: ldim(3), nframes, frame_from, frame_to, numlen, cnt
        integer          :: iframe, jframe, nfiles, endit
        call cline%set('oritype','mic')
        ! check input
        if( cline%defined('nframesgrp') )then
            if( nint(cline%get_rarg('nframesgrp')) < 3 )then
                THROW_HARD('nframesgrp integer (nr of frames to average) needs to be >= 3; exec_tseries_extract')
            endif
        else
            THROW_HARD('need nframesgrp integer input = nr of frames to average; exec_tseries_extract')
        endif
        select case(trim(params%ctf))
            case('yes','no','flip')
                ! fine
            case DEFAULT
                THROW_HARD('Invalid CTF flag')
        end select
        ! build parameters
        call params%new(cline)
        call read_filetable(params%filetab, filenames)
        nfiles = size(filenames)
        numlen = len(int2str(nfiles))
        if( params%mkdir.eq.'yes' )then
            do iframe = 1,nfiles
                if(filenames(iframe)(1:1).ne.'/') filenames(iframe) = '../'//trim(filenames(iframe))
            enddo
        endif
        call find_ldim_nptcls(filenames(1),ldim,nframes)
        if( nframes == 1 .and. ldim(3) == 1 )then
            ! all ok
            call frame_img%new(ldim, params%smpd)
        else
            write(logfhandle,*) 'ldim(3): ', ldim(3)
            write(logfhandle,*) 'nframes: ', nframes
            THROW_HARD('exec_tseries_import assumes one frame per file')
        endif
        ! generates 'movies' to import
        call spproj%read(params%projfile)
        endit = nfiles - params%nframesgrp + 1
        allocate(movfnames(endit))
        do iframe=1,endit
            call progress(iframe, endit)
            if( cline%defined('fbody') )then
                outfname = 'tseries_frames'//int2str_pad(iframe,numlen)//params%ext
            else
                outfname = trim(params%fbody)//'tseries_frames'//int2str_pad(iframe,numlen)//params%ext
            endif
            movfnames(iframe) = trim(outfname)
            frame_from = iframe
            frame_to   = iframe + params%nframesgrp - 1
            cnt = 0
            do jframe=frame_from,frame_to
                cnt = cnt + 1
                call frame_img%read(filenames(jframe),1)
                call frame_img%write(movfnames(iframe),cnt)
            end do
        end do
        call frame_img%kill
        ! import
        select case(trim(params%ctf))
            case('yes')
                ctfvars%ctfflag = CTFFLAG_YES
            case('no')
                ctfvars%ctfflag = CTFFLAG_YES
            case('flip')
                ctfvars%ctfflag = CTFFLAG_FLIP
            case DEFAULT
        end select
        ctfvars%smpd    = params%smpd
        ctfvars%cs      = params%cs
        ctfvars%kv      = params%kv
        ctfvars%fraca   = params%fraca
        ctfvars%l_phaseplate = .false.
        ctfvars%dfx     = 0.
        ctfvars%dfy     = 0.
        ctfvars%angast  = 0.
        ctfvars%phshift = 0.
        call spproj%add_movies( movfnames, ctfvars )
        !final write
        call spproj%write
        ! end gracefully
        call simple_end('**** TSERIES_IMPORT NORMAL STOP ****')
    end subroutine exec_tseries_import

    subroutine exec_tseries_track( self, cline )
        use simple_tseries_tracker
        use simple_qsys_funs,  only: qsys_job_finished
        use simple_sp_project, only: sp_project
        class(tseries_track_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(sp_project)                       :: spproj
        type(parameters)                       :: params
        type(nrtxtfile)                        :: boxfile
        character(len=:),          allocatable :: dir
        character(len=LONGSTRLEN), allocatable :: framenames(:)
        real,                      allocatable :: boxdata(:,:)
        integer :: i, iframe, ndatlines, j, orig_box, numlen, nmovs
        call params%new(cline)
        numlen   = 5 ! default value
        orig_box = params%box
        ! boxfile input
        if( cline%defined('boxfile') )then
            if( .not. file_exists(params%boxfile) ) THROW_HARD('inputted boxfile does not exist in cwd')
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
                        THROW_HARD('Only square windows allowed!')
                    endif
                end do
            else
                THROW_HARD('inputted boxfile is empty; exec_tseries_track')
            endif
        else if( cline%defined('xcoord') .and. cline%defined('ycoord') )then
            if( .not. cline%defined('box') ) THROW_HARD('need box to be part of command line for this mode of execution; exec_tseries_track')
            allocate( boxdata(1,2) )
            boxdata(1,1) = real(params%xcoord)
            boxdata(1,2) = real(params%ycoord)
            ndatlines = 1
        else
            THROW_HARD('need either boxfile or xcoord/ycoord to be part of command line; exec_tseries_track')
        endif
        ! frames input
        call spproj%read(params%projfile)
        nmovs = spproj%get_nintgs()
        allocate(framenames(nmovs))
        iframe = 0
        do i = 1,spproj%os_mic%get_noris()
            if( spproj%os_mic%isthere(i,'intg') )then
                iframe = iframe + 1
                framenames(iframe) = trim(spproj%os_mic%get_static(i,'intg'))
            endif
        enddo
        ! actual tracking
        do j=1,ndatlines
            call init_tracker( nint(boxdata(j,1:2)), framenames)
            call track_particle
            if( cline%defined('ind') )then
                dir = trim(params%fbody)//int2str(params%ind)
                call simple_mkdir(dir)
                if( .not. cline%defined('numlen') ) THROW_HARD('need numlen to be part of command line if ind is; exec_tseries_track')
                call write_tracked_series(dir,trim(params%fbody)//int2str_pad(params%ind,params%numlen))
            else
                dir = trim(params%fbody)//int2str(params%ind)
                call simple_mkdir(dir)
                call write_tracked_series(dir,trim(params%fbody)//int2str_pad(j,numlen))
            endif
            call kill_tracker
        end do
        call qsys_job_finished(  'simple_commander_tseries :: exec_tseries_track')
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
                if( .not. params_present(1) ) write(logfhandle,*) 'ERROR! input deftab lacks kv'
                if( .not. params_present(2) ) write(logfhandle,*) 'ERROR! input deftab lacks cs'
                if( .not. params_present(3) ) write(logfhandle,*) 'ERROR! input deftab lacks fraca'
                if( .not. params_present(4) ) write(logfhandle,*) 'ERROR! input deftab lacks dfx'
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
            call simple_mkdir('tseries_chunk'//int2str_pad(cnt,numlen),errmsg="commander_tseries::exec_tseries_split")
            stkname = filepath(PATH_HERE,'tseries_chunk'//int2str_pad(cnt,numlen),'imgs'//params%ext)
            oriname = filepath('tseries_chunk'//int2str_pad(cnt,numlen),'oris'//trim(TXT_EXT))
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
