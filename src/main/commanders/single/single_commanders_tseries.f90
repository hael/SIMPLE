!@descr: commanders operating on the full time-series field of view, used in SINGLE for nanoparticle processing
module single_commanders_tseries
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_tseries_import
  contains
    procedure :: execute      => exec_tseries_import
end type commander_tseries_import

type, extends(commander_base) :: commander_tseries_motion_correct_distr
  contains
    procedure :: execute      => exec_tseries_motion_correct_distr
end type commander_tseries_motion_correct_distr

type, extends(commander_base) :: commander_tseries_motion_correct
  contains
    procedure :: execute      => exec_tseries_motion_correct
end type commander_tseries_motion_correct

type, extends(commander_base) :: commander_tseries_make_pickavg
  contains
    procedure :: execute      => exec_tseries_make_pickavg
end type commander_tseries_make_pickavg

contains

    subroutine exec_tseries_import( self, cline )
        class(commander_tseries_import), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        type(ctfparams)  :: ctfvars
        type(string), allocatable :: filenames(:)
        integer :: iframe, nfiles, numlen
        call cline%set('oritype','mic')
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        ! read files
        call read_filetable(params%filetab, filenames)
        nfiles = size(filenames)
        numlen = len(int2str(nfiles))
        if( params%mkdir.eq.'yes' )then
            do iframe = 1,nfiles
                if(filenames(iframe)%to_char([1,1]).ne.'/') filenames(iframe) = '../'//filenames(iframe)%to_char()
            enddo
        endif
        ! set CTF parameters
        ctfvars%smpd  = params%smpd
        ctfvars%cs    = params%cs
        ctfvars%kv    = params%kv
        ctfvars%fraca = params%fraca
        ! import the individual frames
        call spproj%read(params%projfile)
        call spproj%add_movies(filenames, ctfvars, singleframe=.true.)
        call spproj%write
        ! end gracefully
        call simple_end('**** TSERIES_IMPORT NORMAL STOP ****')
    end subroutine exec_tseries_import

    subroutine exec_tseries_motion_correct_distr( self, cline )
        class(commander_tseries_motion_correct_distr), intent(inout) :: self
        class(cmdline),                                intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        integer          :: nframes
        call cline%set('oritype',    'mic')
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('nframesgrp') ) call cline%set('nframesgrp',     5)
        if( .not. cline%defined('mcpatch')    ) call cline%set('mcpatch',    'yes')
        if( .not. cline%defined('nxpatch')    ) call cline%set('nxpatch',       10)
        if( .not. cline%defined('nypatch')    ) call cline%set('nypatch',       10)
        if( .not. cline%defined('trs')        ) call cline%set('trs',          10.)
        if( .not. cline%defined('lpstart')    ) call cline%set('lpstart',       5.)
        if( .not. cline%defined('lpstop')     ) call cline%set('lpstop',        3.)
        if( .not. cline%defined('bfac')       ) call cline%set('bfac',          5.)
        if( .not. cline%defined('wcrit')      ) call cline%set('wcrit',  'softmax')
        if( .not. cline%defined('algorithm')  ) call cline%set('algorithm','patch')
        if( .not. cline%defined('downscale')  ) call cline%set('downscale',   'no')
        call params%new(cline)
        call cline%set('numlen', real(params%numlen))
        if( cline%defined('boxfile') )then
            if( .not.file_exists(params%boxfile) ) THROW_HARD('BOXFILE not found!')
        endif
        ! sanity check
        call spproj%read_segment(params%oritype, params%projfile)
        nframes = spproj%get_nframes()
        if( nframes ==0 ) THROW_HARD('no movie frames to process! exec_tseries_motion_correct_distr')
        call spproj%kill
        ! setup the environment for distributed execution
        call qenv%new(params, params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=string(ALGN_FBODY), array=L_USE_SLURM_ARR, extra_params=params)
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(nframes, params%nparts, 'mic', ALGN_FBODY)
        call spproj%kill
        ! clean
        call qsys_cleanup(params)
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_TSERIES_MOTION_CORRECT NORMAL STOP ****')
    end subroutine exec_tseries_motion_correct_distr

    subroutine exec_tseries_motion_correct( self, cline )
        use simple_motion_correct_iter, only: motion_correct_iter
        class(commander_tseries_motion_correct), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        type(string), allocatable :: framenames(:)
        type(string)              :: frames2align
        type(image)               :: img
        type(sp_project)          :: spproj
        type(parameters)          :: params
        type(cmdline)             :: cline_mcorr
        type(motion_correct_iter) :: mciter
        type(ctfparams)           :: ctfvars
        type(ori)                 :: o
        integer :: i, iframe, nframes, frame_counter, ldim(3), fromto(2), nframesgrp
        integer :: numlen_nframes, cnt
        call cline%set('mkdir',       'no') ! shared-memory workflow, dir making in driver
        if( .not. cline%defined('nframesgrp') ) call cline%set('nframesgrp',     5)
        if( .not. cline%defined('mcpatch')    ) call cline%set('mcpatch',    'yes')
        if( .not. cline%defined('nxpatch')    ) call cline%set('nxpatch',       10)
        if( .not. cline%defined('nypatch')    ) call cline%set('nypatch',       10)
        if( .not. cline%defined('trs')        ) call cline%set('trs',          10.)
        if( .not. cline%defined('lpstart')    ) call cline%set('lpstart',       5.)
        if( .not. cline%defined('lpstop')     ) call cline%set('lpstop',        3.)
        if( .not. cline%defined('bfac')       ) call cline%set('bfac',          5.)
        if( .not. cline%defined('wcrit')      ) call cline%set('wcrit',  'softmax')
        if( .not. cline%defined('downscale')  ) call cline%set('downscale',   'no')
        call params%new(cline)
        call spproj%read(params%projfile)
        nframes  = spproj%get_nframes()
        numlen_nframes = len(int2str(nframes))
        allocate(framenames(nframes))
        do i = 1,nframes
            if( spproj%os_mic%isthere(i,'frame') )then
                framenames(i) = spproj%os_mic%get_str(i,'frame')
                params%smpd   = spproj%os_mic%get(i,'smpd')
            endif
        enddo
        call cline%set('smpd', params%smpd)
        frames2align = 'frames2align_part'//int2str_pad(params%part,params%numlen)//'.mrc'
        ! prepare 4 motion_correct
        ldim = [nint(spproj%os_mic%get(1,'xdim')), nint(spproj%os_mic%get(1,'ydim')), 1]
        call img%new(ldim, ctfvars%smpd)
        cline_mcorr = cline
        call cline_mcorr%delete('nframesgrp')
        nframesgrp = params%nframesgrp
        params%nframesgrp = 0
        call cline_mcorr%set('prg', 'motion_correct')
        call cline_mcorr%set('mkdir', 'no')
        ctfvars%smpd = params%smpd
        do iframe=params%fromp,params%top
            call spproj%os_mic%get_ori(iframe, o)
            ! set time window
            fromto(1) = iframe - (nframesgrp-1)/2
            fromto(2) = fromto(1) + nframesgrp - 1
            ! shift the window if it's outside the time-series
            do while(fromto(1) < 1)
                fromto = fromto + 1
            end do
            do while(fromto(2) > nframes)
                fromto = fromto - 1
            end do
            ! make stack
            cnt = 0
            do i = fromto(1),fromto(2)
                cnt = cnt + 1
                call img%read(framenames(i))
                call img%write(frames2align, cnt)
            enddo
            ! motion corr
            frame_counter = 0
            if( cline%defined('gainref') )then
                call mciter%iterate(params, cline_mcorr, ctfvars, o, string('tseries_win')//int2str_pad(iframe,numlen_nframes),&
                &frame_counter, frames2align, string('./'), tseries='yes', gainref_fname=params%gainref)
            else
                call mciter%iterate(params, cline_mcorr, ctfvars, o, string('tseries_win')//int2str_pad(iframe,numlen_nframes),&
                &frame_counter, frames2align, string('./'), tseries='yes')
            endif
            call spproj%os_mic%set_ori(iframe, o)
        end do
        ! output
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, [params%fromp,params%top], isegment=MIC_SEG)
        ! done!
        call del_file(frames2align)
        call img%kill
        call o%kill
        call qsys_job_finished(params, string('single_commanders_tseries :: exec_tseries_motion_correct'))
        call simple_end('**** SIMPLE_TSERIES_MOTION_CORRECT NORMAL STOP ****')
    end subroutine exec_tseries_motion_correct

    subroutine exec_tseries_make_pickavg( self, cline )
        use simple_motion_correct_iter, only: motion_correct_iter
        use simple_tvfilter,            only: tvfilter
        class(commander_tseries_make_pickavg), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        real, parameter :: LAM_TV = 1.5
        type(string), allocatable :: framenames(:)
        type(sp_project)          :: spproj
        type(parameters)          :: params, params_mc
        type(cmdline)             :: cline_mcorr
        type(motion_correct_iter) :: mciter
        type(ctfparams)           :: ctfvars
        type(ori)                 :: o
        type(tvfilter)            :: tvfilt
        type(image)               :: img_intg, img_frame
        integer                   :: istart, istop, period, win_start, win_stop, nwin, win_nframes
        integer                   :: i, nframes, frame_counter, ldim(3), ifoo, numlen_nframes, iwin
        real                      :: smpd_win
        logical                   :: l_periodic
        type(string)              :: stkname, fname_intg, fname_denoised, fbody_mc
        if( .not. cline%defined('nframesgrp') ) call cline%set('nframesgrp',    100)
        if( .not. cline%defined('mcpatch')    ) call cline%set('mcpatch',     'yes')
        if( .not. cline%defined('nxpatch')    ) call cline%set('nxpatch',         3)
        if( .not. cline%defined('nypatch')    ) call cline%set('nypatch',         3)
        if( .not. cline%defined('trs')        ) call cline%set('trs',           20.)
        if( .not. cline%defined('lpstart')    ) call cline%set('lpstart',        5.)
        if( .not. cline%defined('lpstop')     ) call cline%set('lpstop',         3.)
        if( .not. cline%defined('bfac')       ) call cline%set('bfac',           5.)
        if( .not. cline%defined('wcrit')      ) call cline%set('wcrit',   'softmax')
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',       'yes')
        call cline%set('mcconvention','relion') ! ensures alignment to first frame
        call params%new(cline)
        call spproj%read(params%projfile)
        nframes  = spproj%get_nframes()
        if( params%nframesgrp < 1 ) THROW_HARD('nframesgrp must be >= 1; exec_tseries_make_pickavg')
        numlen_nframes = len(int2str(nframes))
        istart = 1
        if( cline%defined('fromf') ) istart = params%fromf
        if( istart < 1 .or. istart > nframes ) THROW_HARD('fromf is outside available frame range; exec_tseries_make_pickavg')
        period = 0
        if( cline%defined('period') ) period = params%period
        if( period < 0 ) THROW_HARD('period must be >= 0; exec_tseries_make_pickavg')
        l_periodic = period > 0
        nwin = 0
        if( l_periodic )then
            do win_start = istart, nframes, period
                win_stop = win_start + params%nframesgrp - 1
                if( win_stop > nframes ) exit
                nwin = nwin + 1
                win_nframes = win_stop - win_start + 1
                allocate(framenames(win_nframes))
                do i = win_start,win_stop
                    if( spproj%os_mic%isthere(i,'frame') )then
                        iwin = i - win_start + 1
                        framenames(iwin) = spproj%os_mic%get_str(i,'frame')
                        smpd_win         = spproj%os_mic%get(i,'smpd')
                    else
                        THROW_HARD('Missing frame entry in periodic make_pickavg window')
                    endif
                enddo
                call cline%set('smpd', smpd_win)
                stkname  = 'frames2align_from'//int2str_pad(win_start,numlen_nframes)//'.mrc'
                fbody_mc = 'frames2align_from'//int2str_pad(win_start,numlen_nframes)
                ! build temporary stack directly from selected frames
                call find_ldim_nptcls(framenames(1),ldim,ifoo)
                call img_frame%new([ldim(1),ldim(2),1], smpd_win)
                do i=1,win_nframes
                    call img_frame%read(framenames(i))
                    call img_frame%write(stkname, i)
                enddo
                call img_frame%kill
                ! prepare and execute motion correction
                cline_mcorr = cline
                call cline_mcorr%delete('nframesgrp')
                params_mc = params
                params_mc%nframesgrp = 0
                params_mc%smpd = smpd_win
                call cline_mcorr%set('prg', 'motion_correct')
                call cline_mcorr%set('mkdir', 'no')
                call o%new(is_ptcl=.false.)
                ctfvars%smpd  = smpd_win
                frame_counter = 0
                if( cline%defined('gainref') )then
                    call mciter%iterate(params_mc, cline_mcorr, ctfvars, o, fbody_mc, frame_counter, stkname, string('./'),&
                        &gainref_fname=params%gainref, tseries='yes')
                else
                    call mciter%iterate(params_mc, cline_mcorr, ctfvars, o, fbody_mc, frame_counter, stkname, string('./'), tseries='yes')
                endif
                call o%kill
                ! apply TV filter for de-noising
                fname_intg = fbody_mc//'_intg.mrc'
                call find_ldim_nptcls(fname_intg,ldim,ifoo)
                call img_intg%new(ldim, smpd_win)
                call img_intg%read(fname_intg,1)
                call tvfilt%new
                call tvfilt%apply_filter(img_intg, LAM_TV)
                call tvfilt%kill
                fname_denoised = fbody_mc//'_intg_denoised.mrc'
                call img_intg%write(fname_denoised)
                call img_intg%kill
                call del_file(stkname)
                deallocate(framenames)
            enddo
            if( nwin == 0 ) THROW_HARD('No complete window found for periodic make_pickavg')
        else
            istop  = min(nframes, istart + params%nframesgrp - 1)
            win_nframes = istop - istart + 1
            allocate(framenames(win_nframes))
            do i = istart,istop
                if( spproj%os_mic%isthere(i,'frame') )then
                    iwin = i - istart + 1
                    framenames(iwin) = spproj%os_mic%get_str(i,'frame')
                    smpd_win        = spproj%os_mic%get(i,'smpd')
                else
                    THROW_HARD('Missing frame entry in make_pickavg window')
                endif
            enddo
            call cline%set('smpd', smpd_win)
            stkname  = 'frames2align.mrc'
            fbody_mc = 'frames2align'
            ! build temporary stack directly from selected frames
            call find_ldim_nptcls(framenames(1),ldim,ifoo)
            call img_frame%new([ldim(1),ldim(2),1], smpd_win)
            do i=1,win_nframes
                call img_frame%read(framenames(i))
                call img_frame%write(stkname, i)
            enddo
            call img_frame%kill
            ! prepare and execute motion correction
            cline_mcorr = cline
            call cline_mcorr%delete('nframesgrp')
            params_mc = params
            params_mc%nframesgrp = 0
            params_mc%smpd = smpd_win
            call cline_mcorr%set('prg', 'motion_correct')
            call cline_mcorr%set('mkdir', 'no')
            call o%new(is_ptcl=.false.)
            ctfvars%smpd  = smpd_win
            frame_counter = 0
            if( cline%defined('gainref') )then
                call mciter%iterate(params_mc, cline_mcorr, ctfvars, o, fbody_mc, frame_counter, stkname, string('./'),&
                    &gainref_fname=params%gainref, tseries='yes')
            else
                call mciter%iterate(params_mc, cline_mcorr, ctfvars, o, fbody_mc, frame_counter, stkname, string('./'), tseries='yes')
            endif
            call o%kill
            ! apply TV filter for de-noising
            call find_ldim_nptcls(string('frames2align_intg.mrc'),ldim,ifoo)
            call img_intg%new(ldim, smpd_win)
            call img_intg%read(string('frames2align_intg.mrc'),1)
            call tvfilt%new
            call tvfilt%apply_filter(img_intg, LAM_TV)
            call tvfilt%kill
            call img_intg%write(string('frames2align_intg_denoised.mrc'))
            call img_intg%kill
            call del_file(stkname)
            deallocate(framenames)
        endif
        call simple_end('**** SIMPLE_TSERIES_MAKE_PICKAVG NORMAL STOP ****')
    end subroutine exec_tseries_make_pickavg

end module single_commanders_tseries
