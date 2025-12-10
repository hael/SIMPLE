! concrete commander: time-series analysis
module simple_commanders_tseries
include 'simple_lib.f08'
use simple_builder,           only: builder
use simple_cmdline,           only: cmdline
use simple_commander_base,    only: commander_base
use simple_commanders_oris,   only: commander_vizoris
use simple_commanders_rec,    only: commander_reconstruct3D
use simple_commanders_volops, only: commander_reproject
use simple_exec_helpers,      only: set_shmem_flag
use simple_image,             only: image
use simple_image_bin,         only: image_bin
use simple_parameters,        only: parameters, params_glob
use simple_qsys_env,          only: qsys_env
use simple_sp_project,        only: sp_project
use simple_stack_io,          only: stack_io
use simple_binoris_io
use simple_commanders_cluster2D
use simple_nanoparticle
use simple_nice
use simple_qsys_funs
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_tseries_import
  contains
    procedure :: execute      => exec_tseries_import
end type commander_tseries_import

type, extends(commander_base) :: commander_tseries_import_particles
  contains
    procedure :: execute      => exec_tseries_import_particles
end type commander_tseries_import_particles

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

type, extends(commander_base) :: commander_tseries_track_particles_distr
  contains
    procedure :: execute      => exec_tseries_track_particles_distr
end type commander_tseries_track_particles_distr

type, extends(commander_base) :: commander_tseries_track_particles
  contains
    procedure :: execute      => exec_tseries_track_particles
end type commander_tseries_track_particles

type, extends(commander_base) :: commander_analysis2D_nano
  contains
    procedure :: execute      => exec_analysis2D_nano
end type commander_analysis2D_nano

type, extends(commander_base) :: commander_center2D_nano
  contains
    procedure :: execute      => exec_center2D_nano
end type commander_center2D_nano

type, extends(commander_base) :: commander_cluster2D_nano
  contains
    procedure :: execute      => exec_cluster2D_nano
end type commander_cluster2D_nano

type, extends(commander_base) :: commander_tseries_backgr_subtr
  contains
    procedure :: execute      => exec_tseries_backgr_subtr
end type commander_tseries_backgr_subtr

type, extends(commander_base) :: commander_tseries_ctf_estimate
  contains
    procedure :: execute      => exec_tseries_ctf_estimate
end type commander_tseries_ctf_estimate

type, extends(commander_base) :: commander_extract_substk
contains
    procedure :: execute      => exec_extract_substk
end type commander_extract_substk

type, extends(commander_base) :: commander_autorefine3D_nano
  contains
    procedure :: execute      => exec_autorefine3D_nano
end type commander_autorefine3D_nano

type, extends(commander_base) :: commander_cavgsproc_nano
  contains
    procedure :: execute      => exec_cavgsproc_nano
end type commander_cavgsproc_nano

type, extends(commander_base) :: commander_cavgseoproc_nano
  contains
    procedure :: execute      => exec_cavgseoproc_nano
end type commander_cavgseoproc_nano

type, extends(commander_base) :: commander_ptclsproc_nano
  contains
    procedure :: execute      => exec_ptclsproc_nano
end type commander_ptclsproc_nano

type, extends(commander_base) :: commander_refine3D_nano
  contains
    procedure :: execute      => exec_refine3D_nano
end type commander_refine3D_nano

type, extends(commander_base) :: commander_graphene_subtr
  contains
    procedure :: execute      => exec_graphene_subtr
end type commander_graphene_subtr

type, extends(commander_base) :: commander_denoise_trajectory
  contains
    procedure :: execute      => exec_denoise_trajectory
end type commander_denoise_trajectory

type, extends(commander_base) :: commander_tseries_swap_stack
  contains
    procedure :: execute      => exec_tseries_swap_stack
end type commander_tseries_swap_stack

type, extends(commander_base) :: commander_tseries_reconstruct3D_distr
  contains
    procedure :: execute      => exec_commander_tseries_reconstruct3D_distr
end type commander_tseries_reconstruct3D_distr

type, extends(commander_base) :: commander_tseries_core_finder
  contains
    procedure :: execute      => exec_tseries_core_finder
end type commander_tseries_core_finder

type, extends(commander_base) :: commander_tseries_make_projavgs
  contains
    procedure :: execute      => exec_tseries_make_projavgs
end type commander_tseries_make_projavgs

contains

    subroutine exec_tseries_import( self, cline )
        use simple_sp_project, only: sp_project
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

    subroutine exec_tseries_import_particles( self, cline )
        use simple_sp_project,       only: sp_project
        use simple_ctf_estimate_fit, only: ctf_estimate_fit
        class(commander_tseries_import_particles), intent(inout) :: self
        class(cmdline),                            intent(inout) :: cline
        type(parameters)       :: params
        type(sp_project)       :: spproj
        type(ctfparams)        :: ctfvars
        type(ctf_estimate_fit) :: ctffit
        type(oris)             :: os
        integer :: iframe, lfoo(3)
        call cline%set('oritype','mic')
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call spproj%read(params%projfile)
        ! # of particles
        call find_ldim_nptcls(params%stk, lfoo, params%nptcls)
        call os%new(params%nptcls, is_ptcl=.true.)
        ! CTF parameters
        ctfvars%smpd    = params%smpd
        ctfvars%phshift = 0.
        ctfvars%dfx     = 0.
        ctfvars%dfy     = 0.
        ctfvars%angast  = 0.
        ctfvars%l_phaseplate = .false.
        if( cline%defined('deftab') )then
            call ctffit%read_doc(params%deftab)
            call ctffit%get_parms(ctfvars)
            if( .not.is_equal(ctfvars%smpd,params%smpd) )then
                THROW_HARD('Iconsistent sampling distance; exec_tseries_import_particles')
            endif
            ctfvars%ctfflag = CTFFLAG_YES
            do iframe = 1,spproj%os_mic%get_noris()
                call spproj%os_mic%set(iframe,'ctf','yes')
            enddo
            call os%set_all2single('dfx', ctfvars%dfx)
            call os%set_all2single('dfy', ctfvars%dfy)
            call os%set_all2single('angast', ctfvars%angast)
        else
            ctfvars%ctfflag = CTFFLAG_NO
            do iframe = 1,spproj%os_mic%get_noris()
                call spproj%os_mic%set(iframe,'ctf','no')
            enddo
        endif
        ! import stack
        call os%set_all2single('state', 1.0)
        call spproj%add_single_stk(params%stk, ctfvars, os)
        call spproj%write
        ! end gracefully
        call simple_end('**** TSERIES_IMPORT_PARTICLES NORMAL STOP ****')
    end subroutine exec_tseries_import_particles

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
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=string(ALGN_FBODY), array=L_USE_SLURM_ARR)
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(nframes, params%nparts, 'mic', ALGN_FBODY)
        call spproj%kill
        ! clean
        call qsys_cleanup
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
                call mciter%iterate(cline_mcorr, ctfvars, o, string('tseries_win')//int2str_pad(iframe,numlen_nframes),&
                &frame_counter, frames2align, string('./'), tseries='yes', gainref_fname=params%gainref)
            else
                call mciter%iterate(cline_mcorr, ctfvars, o, string('tseries_win')//int2str_pad(iframe,numlen_nframes),&
                frame_counter, frames2align, string('./'), tseries='yes')
            endif
            call spproj%os_mic%set_ori(iframe, o)
        end do
        ! output
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, [params%fromp,params%top], isegment=MIC_SEG)
        ! done!
        call del_file(frames2align)
        call img%kill
        call o%kill
        call qsys_job_finished(string('simple_commanders_tseries :: exec_tseries_motion_correct'))
        call simple_end('**** SIMPLE_TSERIES_MOTION_CORRECT NORMAL STOP ****')
    end subroutine exec_tseries_motion_correct

    subroutine exec_tseries_make_pickavg( self, cline )
        use simple_commanders_imgproc,   only: commander_stack
        use simple_motion_correct_iter, only: motion_correct_iter
        use simple_tvfilter,            only: tvfilter
        class(commander_tseries_make_pickavg), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        real, parameter :: LAM_TV = 1.5
        type(string), allocatable :: framenames(:)
        type(string)              :: filetabname
        type(sp_project)          :: spproj
        type(parameters)          :: params
        type(cmdline)             :: cline_stack, cline_mcorr
        type(commander_stack)     :: xstack
        type(motion_correct_iter) :: mciter
        type(ctfparams)           :: ctfvars
        type(ori)                 :: o
        type(tvfilter)            :: tvfilt
        type(image)               :: img_intg
        integer                   :: istart, istop
        integer :: i, nframes, frame_counter, ldim(3), ifoo, cnt
        if( .not. cline%defined('nframesgrp') ) call cline%set('nframesgrp',     10)
        if( .not. cline%defined('mcpatch')    ) call cline%set('mcpatch',     'yes')
        if( .not. cline%defined('nxpatch')    ) call cline%set('nxpatch',         3)
        if( .not. cline%defined('nypatch')    ) call cline%set('nypatch',         3)
        if( .not. cline%defined('trs')        ) call cline%set('trs',           10.)
        if( .not. cline%defined('lpstart')    ) call cline%set('lpstart',        5.)
        if( .not. cline%defined('lpstop')     ) call cline%set('lpstop',         3.)
        if( .not. cline%defined('bfac')       ) call cline%set('bfac',           5.)
        if( .not. cline%defined('wcrit')      ) call cline%set('wcrit',   'softmax')
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',       'yes')
        call cline%set('mcconvention','relion') ! ensures alignment to first frame
        call params%new(cline)
        call spproj%read(params%projfile)
        nframes  = spproj%get_nframes()
        allocate(framenames(params%nframesgrp))
        istart = 1
        if( cline%defined('fromf') ) istart = params%fromf
        istop  = istart + params%nframesgrp - 1
        cnt    = 0
        do i = istart,istop
            if( spproj%os_mic%isthere(i,'frame') )then
                cnt = cnt + 1
                framenames(cnt) = spproj%os_mic%get_str(i,'frame')
                params%smpd     = spproj%os_mic%get(i,'smpd')
            endif
        enddo
        call cline%set('smpd', params%smpd)
        filetabname = 'filetab_first'//int2str(params%nframesgrp)//'frames.txt'
        call write_filetable(filetabname, framenames)
        ! prepare stack command
        call cline_stack%set('mkdir',   'no')
        call cline_stack%set('filetab', filetabname)
        call cline_stack%set('outstk',  'frames2align.mrc')
        call cline_stack%set('smpd',    params%smpd)
        call cline_stack%set('nthr',    1)
        ! execute stack commander
        call xstack%execute(cline_stack)
        ! prepare 4 motion_correct
        cline_mcorr = cline
        call cline_mcorr%delete('nframesgrp')
        params%nframesgrp = 0
        call cline_mcorr%set('prg', 'motion_correct')
        call cline_mcorr%set('mkdir', 'no')
        call o%new(is_ptcl=.false.)
        ctfvars%smpd  = params%smpd
        frame_counter = 0
        ! motion corr
        if( cline%defined('gainref') )then
            call mciter%iterate(cline_mcorr, ctfvars, o, string('frames2align'), frame_counter,&
                &string('frames2align.mrc'), string('./'), gainref_fname=params%gainref, tseries='yes')
        else
            call mciter%iterate(cline_mcorr, ctfvars, o, string('frames2align'), frame_counter,&
                &string('frames2align.mrc'), string('./'), tseries='yes')
        endif
        call o%kill
        ! apply TV filter for de-noising
        call find_ldim_nptcls(string('frames2align_intg.mrc'),ldim,ifoo)
        call img_intg%new(ldim, params%smpd)
        call img_intg%read(string('frames2align_intg.mrc'),1)
        call tvfilt%new
        call tvfilt%apply_filter(img_intg, LAM_TV)
        call tvfilt%kill
        call img_intg%write(string('frames2align_intg_denoised.mrc'))
        call img_intg%kill
        call simple_end('**** SIMPLE_TSERIES_MAKE_PICKAVG NORMAL STOP ****')
    end subroutine exec_tseries_make_pickavg

    subroutine exec_tseries_track_particles_distr( self, cline )
        class(commander_tseries_track_particles_distr), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters)              :: params
        type(qsys_env)                :: qenv
        type(chash)                   :: job_descr
        type(nrtxtfile)               :: boxfile
        real,        allocatable      :: boxdata(:,:)
        type(chash), allocatable      :: part_params(:)
        integer :: ndatlines, numlen, j, orig_box, ipart
        if( .not. cline%defined('neg')       ) call cline%set('neg',      'yes')
        if( .not. cline%defined('lp')        ) call cline%set('lp',         2.3)
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      5.0)
        if( .not. cline%defined('nframesgrp')) call cline%set('nframesgrp',  30)
        if( .not. cline%defined('filter'))     call cline%set('filter',    'tv')
        if( .not. cline%defined('offset'))     call cline%set('offset',     10.)
        call cline%set('oritype','mic')
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        if( .not. file_exists(params%boxfile) ) THROW_HARD('inputted boxfile does not exist in cwd')
        if( nlines(params%boxfile) > 0 )then
            call boxfile%new(params%boxfile, 1)
            ndatlines = boxfile%get_ndatalines()
            numlen    = len(int2str(ndatlines))
            allocate( boxdata(ndatlines,boxfile%get_nrecs_per_line()) )
            do j=1,ndatlines
                call boxfile%readNextDataLine(boxdata(j,:))
                orig_box = nint(boxdata(j,3))
                if( nint(boxdata(j,3)) /= nint(boxdata(j,4)) )then
                    THROW_HARD('Only square windows allowed!')
                endif
            end do
        else
            THROW_HARD('inputted boxfile is empty; exec_tseries_track_particles')
        endif
        call boxfile%kill
        params%nptcls  = ndatlines
        params%nparts  = params%nptcls
        params%ncunits = params%nparts
        ! box and numlen need to be part of command line
        if( .not. cline%defined('hp') ) call cline%set('hp', real(orig_box) )
        call cline%set('box',    real(orig_box))
        call cline%set('numlen', real(numlen)  )
        call cline%delete('fbody')
        ! prepare part-dependent parameters
        allocate(part_params(params%nparts))
        do ipart=1,params%nparts
            call part_params(ipart)%new(4)
            call part_params(ipart)%set('xcoord', real2str(boxdata(ipart,1)))
            call part_params(ipart)%set('ycoord', real2str(boxdata(ipart,2)))
            call part_params(ipart)%set('box',    real2str(boxdata(ipart,3)))
            if( params%nparts > 1 )then
                call part_params(ipart)%set('fbody', params%fbody%to_char()//'_'//int2str_pad(ipart,numlen))
            else
                call part_params(ipart)%set('fbody', params%fbody%to_char())
            endif
        end do
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! schedule & clean
        call cline%gen_job_descr(job_descr)
        call qenv%gen_scripts_and_schedule_jobs(job_descr, part_params=part_params, array=L_USE_SLURM_ARR)
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_TRACK_PARTICLES NORMAL STOP ****')
    end subroutine exec_tseries_track_particles_distr

    subroutine exec_tseries_track_particles( self, cline )
        use simple_tseries_track_particles
        use simple_qsys_funs,  only: qsys_job_finished
        use simple_sp_project, only: sp_project
        class(commander_tseries_track_particles), intent(inout) :: self
        class(cmdline),                           intent(inout) :: cline
        type(sp_project)          :: spproj
        type(parameters)          :: params
        type(string)              :: dir, forctf
        type(string), allocatable :: intg_names(:), frame_names(:)
        real,         allocatable :: boxdata(:,:)
        logical,      allocatable :: frames_are_there(:), intgs_are_there(:)
        integer :: i, orig_box, nframes
        call cline%set('oritype','mic')
        call params%new(cline)
        orig_box = params%box
        ! coordinates input
        if( cline%defined('xcoord') .and. cline%defined('ycoord') )then
            if( .not. cline%defined('box') ) THROW_HARD('need box to be part of command line for this mode of execution; exec_tseries_track_particles')
            allocate( boxdata(1,2) )
            boxdata(1,1) = real(params%xcoord)
            boxdata(1,2) = real(params%ycoord)
        else
            THROW_HARD('need xcoord/ycoord to be part of command line; exec_tseries_track_particles')
        endif
        ! frames input
        call spproj%read(params%projfile)
        nframes = spproj%get_nframes()
        allocate(frames_are_there(nframes), intgs_are_there(nframes), source=.false.)
        do i = 1,nframes
            frames_are_there(i) = spproj%os_mic%isthere(i,'frame')
            intgs_are_there(i)  = spproj%os_mic%isthere(i,'intg')
        enddo
        if( all(frames_are_there) )then
            allocate(frame_names(nframes))
            do i = 1,nframes
                frame_names(i) = spproj%os_mic%get_str(i,'frame')
            enddo
        endif
        if( all(intgs_are_there) )then
            allocate(intg_names(nframes))
            do i = 1,nframes
                intg_names(i) = spproj%os_mic%get_str(i,'intg')
            enddo
        endif
        ! actual tracking
        dir = params%fbody
        call simple_mkdir(dir)
        if( allocated(frame_names) )then
            if( allocated(intg_names) )then
                call init_tracker( nint(boxdata(1,1:2)), intg_names, frame_names, dir, params%fbody)
            else
                call init_tracker( nint(boxdata(1,1:2)), frame_names, frame_names, dir, params%fbody)
            endif
        endif
        if( cline%defined('fromf') )then
            call track_particle( forctf, params%fromf )
        else
            call track_particle( forctf )
        endif
        ! clean tracker
        call kill_tracker
        ! end gracefully
        call qsys_job_finished(string('simple_commanders_tseries :: exec_tseries_track_particles'))
        call spproj%kill
        call simple_end('**** SIMPLE_TSERIES_TRACK_PARTICLES NORMAL STOP ****')
    end subroutine exec_tseries_track_particles

    subroutine exec_analysis2D_nano( self, cline )
        use simple_commanders_imgproc, only: commander_estimate_diam
        use simple_commanders_sim,     only: commander_simulate_atoms
        class(commander_analysis2D_nano), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        ! commanders
        type(commander_center2D_nano)  :: xcenter2D
        type(commander_cluster2D_nano) :: xcluster2D
        type(commander_estimate_diam)  :: xest_diam
        type(commander_simulate_atoms) :: xsim_atms
        ! other variables
        type(simple_nice_communicator) :: nice_communicator
        type(parameters)               :: params
        type(sp_project)               :: spproj
        type(cmdline)                  :: cline_est_diam, cline_sim_atms, cline_copy
        type(string)                   :: stkname
        character(len=*), parameter    :: STARTVOL    = 'startvol.mrc'
        real,             parameter    :: LP_EST_DIAM = 3.
        integer :: ncls, nptcls, ldim(3)
        real    :: smpd, diam_min, diam_max, mskdiam
        call cline%set('dir_exec', 'analysis2D_nano')
        if( .not. cline%defined('objfun')  ) call cline%set('objfun', 'cc') ! best objfun
        if( .not. cline%defined('ml_reg')  ) call cline%set('ml_reg', 'no') ! ml_reg=yes -> too few atoms 
        call params%new(cline)
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        call nice_communicator%cycle()
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        cline_copy = cline
        ! centering
        call cline%delete('nptcls_per_cls')
        if( cline%defined('center') .and. trim(params%center).eq.'yes' )then
            call cline%set('center', 'yes')
        else
            call cline%set('center', 'no')
        endif
        call xcenter2D%execute(cline)
        ! prep for diameter estimation
        call spproj%read(params%projfile)
        call spproj%get_cavgs_stk(stkname, ncls, smpd)
        call cline_est_diam%set('stk',     stkname)
        call cline_est_diam%set('smpd',    smpd)
        call cline_est_diam%set('mskdiam', 0.)
        call cline_est_diam%set('nthr',    params%nthr)
        call cline_est_diam%set('lp',      LP_EST_DIAM)
        ! estimate diameter
        call xest_diam%execute(cline_est_diam)
        diam_min = cline_est_diam%get_rarg('min_diam')
        diam_max = cline_est_diam%get_rarg('max_diam')
        mskdiam  = diam_max * 1.5
        call cline%set('mskdiam', mskdiam)
        write(logfhandle,'(A,2F6.1)') '>>> MASK DIAMETER (MSKDIAM) (IN A): ', mskdiam
        ! make a starting volume for initialization of 3D refinement
        call find_ldim_nptcls(stkname, ldim, nptcls)
        call cline_sim_atms%set('outvol',  STARTVOL)
        call cline_sim_atms%set('smpd',    params%smpd)
        call cline_sim_atms%set('element', params%element)
        call cline_sim_atms%set('moldiam', diam_min)
        call cline_sim_atms%set('box',     ldim(1))
        call cline_sim_atms%set('nthr',    params%nthr)
        call xsim_atms%execute(cline_sim_atms)
        ! run final 2D analysis
        cline = cline_copy
        call exec_cmdline('rm -rf cavgs* clusters2D*star *_FINISHED start2Drefs* frcs*')
        call cline%set('center', 'no')
        call xcluster2D%execute(cline)
        ! end gracefully
        call nice_communicator%terminate()
        call simple_end('**** SIMPLE_ANALYSIS2D_NANO NORMAL STOP ****')
    end subroutine exec_analysis2D_nano

    subroutine exec_center2D_nano( self, cline )
        class(commander_center2D_nano), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        ! commanders
        type(commander_cluster2D_nano)   :: xcluster2D_nano ! shared-memory by default
        type(commander_make_cavgs_distr) :: xmake_cavgs
        ! constants
        integer, parameter               :: NCLS_CEN_NANO = 10
        ! other variables
        type(parameters)                 :: params
        type(sp_project)                 :: spproj
        type(cmdline)                    :: cline_make_cavgs, cline_cluster2D_nano
        type(string)                     :: orig_projfile
        type(string)                     :: finalcavgs
        integer :: last_iter_stage2, nptcls
        call cline%set('dir_exec', 'center2D_nano')
        if( .not. cline%defined('center') ) call cline%set('center', 'no')
        ! master parameters
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! read project file
        call spproj%read(params%projfile)
        orig_projfile = params%projfile
        ! sanity checks
        nptcls = spproj%get_nptcls()
        if( nptcls == 0 )then
            THROW_HARD('No particles found in project file: '//params%projfile%to_char()//'; exec_center2D_nano')
        endif
        ! delete any previous solution
        if( .not. spproj%is_virgin_field(params%oritype) )then
            ! removes previous cluster2D solution (states are preserved)
            call spproj%os_ptcl2D%delete_2Dclustering
            call spproj%write_segment_inside(params%oritype)
        endif
        ! make initial class averages (averages of time-chunks)
        cline_make_cavgs = cline
        call cline_make_cavgs%set('prg',     'make_cavgs')
        if( .not. cline%defined('ncls'))then
            call cline_make_cavgs%set('ncls', NCLS_CEN_NANO)
        endif
        call cline_make_cavgs%set('tseries', 'yes')
        call cline_make_cavgs%set('nparts',   1)
        call cline_make_cavgs%set('refs',    'start2Drefs'//params%ext%to_char())
        call cline_make_cavgs%set('projfile', params%projfile)
        call xmake_cavgs%execute_safe(cline_make_cavgs)
        ! do centering
        cline_cluster2D_nano = cline
        call cline_cluster2D_nano%set('prg',     'cluster2D_nano')
        call cline_cluster2D_nano%set('mskdiam',  0.)
        call cline_cluster2D_nano%set('refine',  'inpl')
        call cline_cluster2D_nano%set('projfile', params%projfile)
        call xcluster2D_nano%execute_safe(cline_cluster2D_nano)        
        last_iter_stage2 = cline_cluster2D_nano%get_iarg('endit')
        finalcavgs       = CAVGS_ITER_FBODY//int2str_pad(last_iter_stage2,3)//params%ext%to_char()
        ! adding cavgs & FRCs to project
        params%projfile = orig_projfile
        call spproj%read( params%projfile )
        call spproj%add_frcs2os_out(string(FRCS_FILE), 'frc2D')
        call spproj%add_cavgs2os_out(finalcavgs, spproj%get_smpd(), imgkind='cavg')
        ! transfer 2D shifts to 3D field
        call spproj%map2Dshifts23D
        call spproj%write
        call spproj%kill
        ! cleanup
        call del_file('start2Drefs'//params%ext%to_char())
        ! end gracefully
        call simple_end('**** SIMPLE_CENTER2D_NANO NORMAL STOP ****')
    end subroutine exec_center2D_nano

    subroutine exec_cluster2D_nano( self, cline )
        class(commander_cluster2D_nano), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        ! commander
        type(commander_cluster2D) :: xcluster2D ! shared-memory
        type(string) :: str_refine
        ! static parameters
        call cline%delete('nparts') ! always shared-memory
        call cline%set('prg',           'cluster2D')
        call cline%set('dir_exec', 'cluster2D_nano')
        call cline%set('autoscale',            'no')
        if( .not. cline%defined('tseries') ) call cline%set('tseries', 'yes')
        if( .not. cline%defined('refine')  ) call cline%set('refine','greedy')
        str_refine = cline%get_carg('refine')
        select case(str_refine%to_char())
            case('no','greedy')
                call cline%set('refine','greedy')
                if( .not. cline%defined('nptcls_per_cls') ) call cline%set('nptcls_per_cls', 20)
                if( .not. cline%defined('maxits')         ) call cline%set('maxits',         15)
            case('inpl')
                if( .not. cline%defined('maxits')         ) call cline%set('maxits',         10)
            case DEFAULT
                THROW_HARD('Unsupported refinement mode!')
        end select
        if( .not. cline%defined('center')         ) call cline%set('center',       'yes')
        if( .not. cline%defined('graphene_filt')  ) call cline%set('graphene_filt', 'no')
        if( .not. cline%defined('hp')             ) call cline%set('hp',             3.0)
        if( .not. cline%defined('lp')             ) call cline%set('lp',             1.0)
        if( .not. cline%defined('cenlp')          ) call cline%set('cenlp',           5.)
        if( .not. cline%defined('trs')            ) call cline%set('trs',             5.)
        if( .not. cline%defined('objfun')         ) call cline%set('objfun',        'cc') ! best objfun
        if( .not. cline%defined('ml_reg')         ) call cline%set('ml_reg',        'no') ! ml_reg=yes -> too few atoms 
        if( .not. cline%defined('oritype')        ) call cline%set('oritype',   'ptcl2D')
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        call xcluster2D%execute_safe(cline)
        call str_refine%kill
        call simple_end('**** SIMPLE_CLUSTER2D_NANO NORMAL STOP ****')
    end subroutine exec_cluster2D_nano

    subroutine exec_tseries_backgr_subtr( self, cline )
        ! for background subtraction in time-series data. The goal is to subtract the two graphene
        ! peaks @ 2.14 A and @ 1.23 A. This is done by band-pass filtering the background image,
        ! recommended (and default settings) are hp=5.0 lp=1.1 and width=5.0.
        use simple_ctf,   only: ctf
        class(commander_tseries_backgr_subtr), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(image)      :: img_backgr, img_backgr_wctf, ave_img
        type(ctf)        :: tfun
        type(ctfparams)  :: ctfvars
        type(string)     :: ext,imgname
        real             :: ave,sdev,minv,maxv
        integer          :: iptcl, nptcls, ldim(3)
        logical          :: do_flip, err
        call cline%set('oritype','mic')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! dimensions
        call find_ldim_nptcls(params%stk,ldim,nptcls)
        ! CTF logics
        do_flip = .false.
        if( (build%spproj_field%isthere('ctf')) .and. (cline%defined('ctf').and.params%ctf.eq.'flip') )then
            if( build%spproj%get_nintgs() /= nptcls ) THROW_HARD('Incompatible # of images and micrographs!')
            if( .not.build%spproj_field%isthere('dfx') ) THROW_HARD('Missing CTF parameters')
            do_flip = .true.
        endif
        ! get background image
        call img_backgr%new([params%box,params%box,1], params%smpd)
        call img_backgr_wctf%new([params%box,params%box,1], params%smpd)
        call img_backgr%read(params%stk_backgr, 1)
        ! background image is skipped if sdev is small because this neighbour was
        ! outside the micrograph
        call img_backgr%stats(ave, sdev, maxv, minv, errout=err)
        if( sdev>TINY .and. .not.err )then
            call ave_img%new([params%box,params%box,1], params%smpd)
            ave_img = 0.
            ! main loop
            do iptcl=1,nptcls
                call progress(iptcl,nptcls)
                ! read particle image
                call build%img%read(params%stk, iptcl)
                ! copy background image & CTF
                img_backgr_wctf = img_backgr
                ! fwd ft
                call build%img%fft()
                call img_backgr_wctf%fft()
                if( do_flip )then
                    ctfvars = build%spproj_field%get_ctfvars(iptcl)
                    tfun = ctf(ctfvars%smpd, ctfvars%kv, ctfvars%cs, ctfvars%fraca)
                    call tfun%apply_serial(build%img, 'flip', ctfvars)
                    call tfun%apply_serial(img_backgr_wctf, 'flip', ctfvars)
                endif
                ! filter background image
                call img_backgr_wctf%bp(params%hp,params%lp,width=params%width)
                ! subtract background
                call build%img%subtr(img_backgr_wctf)
                ! bwd ft
                call build%img%ifft()
                ! normalise
                call build%img%norm()
                call ave_img%add(build%img)
                ! output corrected image
                call build%img%write(params%outstk, iptcl)
            end do
            ! generates average
            call ave_img%div(real(nptcls))
            ext     = fname2ext(params%outstk)
            imgname = get_fbody(params%outstk,ext,.true.)//'_ave.'//ext%to_char()
            call ave_img%write(imgname)
        endif
        call ave_img%kill
        call img_backgr%kill
        call img_backgr_wctf%kill
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_BACKGR_SUBTR NORMAL STOP ****')
    end subroutine exec_tseries_backgr_subtr

    subroutine exec_tseries_ctf_estimate( self, cline )
        use simple_ctf_estimate_fit, only: ctf_estimate_fit
        class(commander_tseries_ctf_estimate), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        character(len=*), parameter :: pspec_fname = 'tseries_ctf_estimate_pspec.mrc'
        character(len=*), parameter :: diag_fname  = 'tseries_ctf_estimate_diag'//JPG_EXT
        integer,          parameter :: nmics4ctf    = 10
        type(parameters)       :: params
        type(builder)          :: build
        type(ctf_estimate_fit) :: ctffit
        type(ctfparams)        :: ctfvars
        type(string)           :: fname_diag, tmpl_fname, docname
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('hp')      ) call cline%set('hp', 5.)
        if( .not. cline%defined('lp')      ) call cline%set('lp', 1.)
        if( .not. cline%defined('dfmin')   ) call cline%set('dfmin', -0.05)
        if( .not. cline%defined('dfmax')   ) call cline%set('dfmax',  0.05)
        if( .not. cline%defined('astigtol')) call cline%set('astigtol', 0.001)
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        ! prep
        tmpl_fname = get_fbody(basename(params%stk), params%ext, separator=.false.)
        fname_diag = filepath('./',tmpl_fname//'_ctf_estimate_diag'//JPG_EXT)
        docname    = filepath('./',tmpl_fname//'_ctf'//TXT_EXT)
        if( build%spproj%os_mic%get_noris() /= 0 )then
            ctfvars = build%spproj%os_mic%get_ctfvars(1)
        else if( build%spproj%os_stk%get_noris() /= 0 )then
            ctfvars = build%spproj%os_stk%get_ctfvars(1)
        else
            THROW_HARD('Insufficient information found in the project')
        endif
        ! command-line override
        if( cline%defined('cs').or.cline%defined('kv').or.cline%defined('fraca') )then
            if( .not.(cline%defined('cs').and.cline%defined('kv').and.cline%defined('fraca')) )then
                THROW_HARD('Insufficient number of CTF parameters')
            endif
            ctfvars%cs    = params%cs
            ctfvars%kv    = params%kv
            ctfvars%fraca = params%fraca
        endif
        ctfvars%ctfflag = CTFFLAG_YES
        ! fitting
        call ctffit%fit_nano(params%stk, params%box, ctfvars, [params%dfmin,params%dfmax], [params%hp,params%lp], params%astigtol)
        ! output
        call ctffit%write_doc(params%stk, docname)
        call ctffit%write_diagnostic(fname_diag, nano=.true.)
        ! cleanup
        call ctffit%kill
        call build%spproj%kill
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_CTF_ESTIMATE NORMAL STOP ****')
    end subroutine exec_tseries_ctf_estimate

    subroutine exec_refine3D_nano( self, cline )
        use simple_commanders_refine3D, only: commander_refine3D_distr
        class(commander_refine3D_nano), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        ! commander
        type(commander_refine3D_distr) :: xrefine3D_distr
        ! static parameters
        call cline%set('prg',           'refine3D')
        call cline%set('dir_exec', 'refine3D_nano')
        ! dynamic parameters
        if( .not. cline%defined('cenlp')          ) call cline%set('cenlp',            5.)
        if( .not. cline%defined('graphene_filt')  ) call cline%set('graphene_filt',  'no') ! since Graphene subtraction is part of the workflow
        if( .not. cline%defined('keepvol')        ) call cline%set('keepvol',       'yes')
        if( .not. cline%defined('hp')             ) call cline%set('hp',              3.0)
        if( .not. cline%defined('lp')             ) call cline%set('lp',              1.0)
        if( .not. cline%defined('lpstart_nonuni') ) call cline%set('lpstart_nonuni',  2.5)
        if( .not. cline%defined('lpstop')         ) call cline%set('lpstop',          0.5)
        if( .not. cline%defined('maxits')         ) call cline%set('maxits',           30)
        if( .not. cline%defined('refine')         ) call cline%set('refine',      'neigh')
        if( .not. cline%defined('oritype')        ) call cline%set('oritype',    'ptcl3D')
        if( .not. cline%defined('trs')            ) call cline%set('trs',             5.0)
        if( .not. cline%defined('objfun')         ) call cline%set('objfun',         'cc') ! best objfun for this kind of data
        if( .not. cline%defined('ml_reg')         ) call cline%set('ml_reg',         'no') ! ml_reg=yes -> too few atoms 
        if( .not. cline%defined('sigma_est')      ) call cline%set('sigma_est',  'global') ! only sensible option for this kind of data
        if( .not. cline%defined('icm')            ) call cline%set('icm',           'yes') ! ICM regualrization works 
        if( .not. cline%defined('lambda')         ) call cline%set('lambda',          0.1) ! this is an empirically determined regularization parameter
        call xrefine3D_distr%execute_safe(cline)
    end subroutine exec_refine3D_nano

    subroutine exec_extract_substk( self, cline )
        use simple_image, only: image
        class(commander_extract_substk), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        call cline%set('mkdir', 'no')
        ! init params
        call params%new(cline)
        ! read the project file
        call spproj%read(params%projfile)
        call spproj%write_segment_inside('projinfo')
        call spproj%write_substk([params%fromp,params%top], params%outstk)
        ! end gracefully
        call simple_end('**** SINGLE_EXTRACT_SUBSTK NORMAL STOP ****')
    end subroutine exec_extract_substk

    subroutine exec_autorefine3D_nano( self, cline )
        use simple_commanders_atoms, only: commander_detect_atoms
        class(commander_autorefine3D_nano), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters)              :: params
        type(commander_refine3D_nano) :: xrefine3D_nano
        type(commander_detect_atoms)  :: xdetect_atms
        type(commander_reproject)     :: xreproject
        type(commander_vizoris)       :: xvizoris
        type(commander_make_cavgs)    :: xmake_cavgs
        type(cmdline)                 :: cline_refine3D_nano, cline_detect_atms, cline_reproject
        type(cmdline)                 :: cline_vizoris, cline_make_cavgs
        type(image), allocatable      :: imgs(:)
        type(sp_project)              :: spproj
        character(len=*), parameter   :: RECVOL     = 'recvol_state01.mrc'
        character(len=*), parameter   :: EVEN       = 'recvol_state01_even.mrc'
        character(len=*), parameter   :: ODD        = 'recvol_state01_odd.mrc'
        character(len=*), parameter   :: SIMVOL     = 'recvol_state01_SIM.mrc'
        character(len=*), parameter   :: ATOMS      = 'recvol_state01_ATMS.pdb'
        character(len=*), parameter   :: BINARY     = 'recvol_state01_BIN.mrc'
        character(len=*), parameter   :: CCS        = 'recvol_state01_CC.mrc'
        character(len=*), parameter   :: SPLITTED   = 'split_ccs.mrc'
        character(len=*), parameter   :: FINAL_MAPS = './final_results/'
        character(len=*), parameter   :: TAG        = 'xxx' ! for checking command lines
        integer,          parameter   :: NSPACE_CLS3D = 500
        type(string)                  :: iter_dir, cavgs_stk, fname
        real,             allocatable :: rstates(:), corrs(:)
        logical,          allocatable :: state_mask(:)
        type(string) :: fbody, fbody_split, fname_reprojs, fname_reprojs_sim, fname_cvags_vs_reprojs
        integer      :: i, iter, cnt, cnt2, funit, io_stat, endit
        real         :: smpd
        logical      :: fall_over
        fbody       = get_fbody(RECVOL,   'mrc')
        fbody_split = get_fbody(SPLITTED, 'mrc')
        if(       cline%defined('nparts')         ) call cline%delete('nparts') ! shared-memory workflow
        if( .not. cline%defined('maxits')         ) call cline%set('maxits',          5)
        if( .not. cline%defined('maxits_between') ) call cline%set('maxits_between', 10)
        if( .not. cline%defined('overlap')        ) call cline%set('overlap',      0.98)
        if( .not. cline%defined('fracsrch')       ) call cline%set('fracsrch',      0.9)
        if( .not. cline%defined('objfun')         ) call cline%set('objfun',       'cc') ! needs to be here to avoid ERROR! file sigma2_it_10.star does not exist; simple_fileio.f90; line:   932
        if( .not. cline%defined('trail_rec')      ) call cline%set('trail_rec',   'yes') 
        if( .not. cline%defined('ufrac_trec')     ) call cline%set('ufrac_trec',    0.5)
        call cline%set('mkdir', 'yes') ! because we want to create the directory X_autorefine3D_nano & copy the project file
        call params%new(cline)         ! because the parameters class manages directory creation and project file copying, mkdir = yes
        params%mkdir = 'no'            ! to prevent the input vol to be appended with ../
        call cline%set('mkdir', 'no')  ! because we do not want a nested directory structure in the execution directory
        ! read the project file
        call spproj%read(params%projfile)
        call spproj%write_segment_inside('projinfo')
        ! sanity checks
        rstates = spproj%os_ptcl2D%get_all('state')
        fall_over = .false.
        if( any(rstates < 0.5 ) ) fall_over = .true.
        deallocate(rstates)
        rstates = spproj%os_ptcl3D%get_all('state')
        if( any(rstates < 0.5 ) ) fall_over = .true.
        if( fall_over ) THROW_HARD('There are state=0s in the ptcl2D/3D fields of the project, which is not allowed. Use simple_exec prg=prune_project before executing autorefine3D_nano')
        ! copy the input command line as templates for the refine3D_nano/detect_atoms command line
        cline_refine3D_nano = cline
        cline_detect_atms   = cline
        ! then update cline_refine3D_nano accordingly
        call cline_refine3D_nano%set('prg',     'refine3D_nano')
        call cline_refine3D_nano%set('projfile', params%projfile) ! since we are not making directories (non-standard execution) we need to keep track of project file
        call cline_refine3D_nano%set('keepvol',  'yes')
        call cline_refine3D_nano%set('maxits',   params%maxits_between) ! turn maxits_between into maxits (max # iterations between model building)
        call cline_refine3D_nano%delete('maxits_between')
        ! then update cline_detect_atoms accordingly
        call cline_detect_atms%set('prg', 'detect_atoms')
        call cline_detect_atms%set('vol1', RECVOL)               ! this is ALWYAS going to be the input volume to detect_atoms
        iter = 0
        do i = 1, params%maxits
            ! first refinement pass on the initial volume uses the low-pass limit defined by the user
            call xrefine3D_nano%execute_safe(cline_refine3D_nano)
            call cline_refine3D_nano%set('vol1', SIMVOL)         ! the reference volume is ALWAYS SIMVOL
            call cline_refine3D_nano%delete('lp')                ! uses the default 1.0 A low-pass limit
            endit = cline_refine3D_nano%get_iarg('endit')        ! last iteration executed by refine3D_nano
            call cline_refine3D_nano%delete('endit')             ! used internally but not technically allowed
            call cline_refine3D_nano%set('prg', 'refine3D_nano') ! because the command line is modified refine3D_nano -> refine3D internally
            ! model building
            call xdetect_atms%execute_safe(cline_detect_atms)
            ! copy critical output
            iter_dir = 'iteration_'//int2str_pad(i,2)//'/'
            call simple_mkdir(iter_dir)
            call simple_copy_file(string(RECVOL),   iter_dir//fbody//'_iter'//int2str_pad(i,3)//'.mrc')
            call simple_copy_file(string(EVEN),     iter_dir//fbody//'_iter'//int2str_pad(i,3)//'_even.mrc')
            call simple_copy_file(string(ODD),      iter_dir//fbody//'_iter'//int2str_pad(i,3)//'_odd.mrc')
            call simple_copy_file(string(SIMVOL),   iter_dir//fbody//'_iter'//int2str_pad(i,3)//'_SIM.mrc')
            call simple_copy_file(string(ATOMS),    iter_dir//fbody//'_iter'//int2str_pad(i,3)//'_ATMS.pdb')
            call simple_copy_file(string(BINARY),   iter_dir//fbody//'_iter'//int2str_pad(i,3)//'_BIN.mrc')
            call simple_copy_file(string(CCS),      iter_dir//fbody//'_iter'//int2str_pad(i,3)//'_CC.mrc')
            call simple_copy_file(string(SPLITTED), iter_dir//fbody_split//'_iter'//int2str_pad(i,3)//'.mrc')
            if( params%l_needs_sigma )then
                call simple_copy_file(string(SIGMA2_GROUP_FBODY)//int2str(endit)//STAR_EXT,&
                    &iter_dir//SIGMA2_GROUP_FBODY//int2str_pad(i,3)//STAR_EXT)
            endif
            ! clean
            call exec_cmdline('rm -f recvol_state01_iter*')
            if( params%l_needs_sigma ) call exec_cmdline('rm -f '//SIGMA2_GROUP_FBODY//'*'//STAR_EXT)
            call del_file(ATOMS)
            call del_file(BINARY)
            call del_file(CCS)
            call del_file(SPLITTED)
            iter = iter + 1
        end do
        call xdetect_atms%execute_safe(cline_detect_atms)
        call simple_mkdir(FINAL_MAPS)
        call simple_copy_file(string(RECVOL),   string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'.mrc')
        call simple_copy_file(string(EVEN),     string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'_even.mrc')
        call simple_copy_file(string(ODD),      string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'_odd.mrc')
        call simple_copy_file(string(SIMVOL),   string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'_SIM.mrc')
        call simple_copy_file(string(ATOMS),    string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'_ATMS.pdb')
        call simple_copy_file(string(BINARY),   string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'_BIN.mrc')
        call simple_copy_file(string(CCS),      string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'_CC.mrc')
        call simple_copy_file(string(SPLITTED), string(FINAL_MAPS)//fbody_split//'_iter'//int2str_pad(iter,3)//'.mrc')
        ! clean
        call del_file(SIMVOL)
        call del_file(ATOMS)
        call del_file(BINARY)
        call del_file(CCS)
        call del_file(SPLITTED)
        ! generate 3d classes
        cavgs_stk = 'cavgs3D.mrc'
        call cline_make_cavgs%set('prg',      'make_cavgs')
        call cline_make_cavgs%set('nspace',   NSPACE_CLS3D)
        call cline_make_cavgs%set('pgrp',     params%pgrp)
        call cline_make_cavgs%set('projfile', params%projfile)
        call cline_make_cavgs%set('nthr',     params%nthr)
        call cline_make_cavgs%set('mkdir',    'no')
        call cline_make_cavgs%set('refs',     cavgs_stk)
        call cline_make_cavgs%set('outfile',  'cavgs_oris.txt')
        call cline_make_cavgs%set('ml_reg',   'no')
        call xmake_cavgs%execute_safe(cline_make_cavgs)
        call spproj%os_cls3D%new(NSPACE_CLS3D, is_ptcl=.false.)
        call spproj%os_cls3D%read(string('cavgs_oris.txt')) ! will not be written as part of document
        if( allocated(rstates) ) deallocate(rstates)
        rstates = spproj%os_cls3D%get_all('state')
        ! prepare for re-projection
        call cline_reproject%set('vol1',   string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'.mrc')
        call cline_reproject%set('outstk', 'reprojs_recvol.mrc')
        call cline_reproject%set('smpd',   params%smpd)
        call cline_reproject%set('oritab', 'cavgs_oris.txt')
        call cline_reproject%set('pgrp',   params%pgrp)
        call cline_reproject%set('nthr',   params%nthr)
        call xreproject%execute_safe(cline_reproject)
        call cline_reproject%set('vol1',   string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'_SIM.mrc')
        call cline_reproject%set('outstk', 'reprojs_SIM.mrc')
        ! re-project
        call xreproject%execute_safe(cline_reproject)
        ! write cavgs & reprojections in triplets
        allocate(imgs(3), state_mask(NSPACE_CLS3D))
        call imgs(1)%new([params%box,params%box,1], smpd)
        call imgs(2)%new([params%box,params%box,1], smpd)
        call imgs(3)%new([params%box,params%box,1], smpd)
        cnt  = 0
        cnt2 = 1
        fname_reprojs          = 'reprojs_recvol.mrc'
        fname_reprojs_sim      = 'reprojs_SIM.mrc'
        fname_cvags_vs_reprojs = 'cavgs_vs_reprojections_rec_and_sim.mrc'
        do i = 1,3*NSPACE_CLS3D,3
            cnt = cnt + 1
            if( rstates(cnt) > 0.5 )then
                state_mask(cnt) = .true.
                call imgs(1)%read(cavgs_stk,         cnt)
                call imgs(2)%read(fname_reprojs,     cnt)
                call imgs(3)%read(fname_reprojs_sim, cnt)
                call imgs(1)%norm
                call imgs(2)%norm
                call imgs(3)%norm
                call imgs(1)%write(fname_cvags_vs_reprojs, cnt2    )
                call imgs(2)%write(fname_cvags_vs_reprojs, cnt2 + 1)
                call imgs(3)%write(fname_cvags_vs_reprojs, cnt2 + 2)
                cnt2 = cnt2 + 3
            else
                state_mask(cnt) = .false.
            endif
        end do
        call imgs(1)%kill
        call imgs(2)%kill
        call imgs(3)%kill
        deallocate(imgs)
        call exec_cmdline('rm -rf fsc* fft* recvol* RES* reprojs_recvol* reprojs* cavgs3D*mrc reproject_oris.txt cavgs_oris.txt stderrout')
        ! visualization of particle orientations
        ! read the ptcl3D segment first to make sure that we are using the latest information
        call spproj%read_segment('ptcl3D', params%projfile)
        ! extract ptcls oritab
        call spproj%os_ptcl3D%write(string('ptcls_oris.txt'))
        call cline_vizoris%set('oritab', 'ptcls_oris.txt')
        call cline_vizoris%set('pgrp',        params%pgrp)
        call cline_vizoris%set('nspace',     NSPACE_CLS3D)
        call cline_vizoris%set('tseries',           'yes')
        call xvizoris%execute_safe(cline_vizoris)
        ! print CSV file of correlation vs particle number
        corrs = spproj%os_ptcl3D%get_all('corr')
        fname = 'ptcls_vs_reprojs_corrs.csv'
        call fopen(funit, fname, 'replace', 'unknown', iostat=io_stat, form='formatted')
        call fileiochk('autorefine3D_nano fopen failed'//fname%to_char(), io_stat)
        write(funit,*) 'PTCL_INDEX'//CSV_DELIM//'CORR'
        do i = 1,size(corrs)
            write(funit,*) int2str(i)//CSV_DELIM//real2str(corrs(i))
        end do
        call fclose(funit)
        ! deallocate
        call iter_dir%kill
        call cavgs_stk%kill
        ! end gracefully
        call simple_end('**** AUTOREFINE3D_NANO NORMAL STOP ****')
    end subroutine exec_autorefine3D_nano

    subroutine exec_ptclsproc_nano( self, cline )
        use simple_strategy2D3D_common, only: read_imgbatch, prepimgbatch, discrete_read_imgbatch
        class(commander_ptclsproc_nano), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        integer,       allocatable :: pinds(:)
        type(parameters)     :: params
        type(builder)        :: build
        type(image)          :: cavg_img
        type(sp_project)     :: spproj
        type(string)         :: command_plot
        type(string)         :: cavgs_stk
        integer              :: cnt, iptcl, icavgs, ncavgs, nptcls
        real                 :: smpd, corr
        call cline%set('mkdir', 'yes') ! because we want to create the directory X_ptclsproc_nano & copy the project file
        call params%new(cline)         ! because the parameters class manages directory creation and project file copying, mkdir = yes
        params%mkdir = 'no'            ! to prevent the input vol to be appended with ../
        call cline%set('mkdir', 'no')  ! because we do not want a nested directory structure in the execution directory
        ! read the project file
        call spproj%read(params%projfile)
        call build%init_params_and_build_general_tbox(cline, params)
        ! retrieve number of cavgs stack
        call spproj%get_cavgs_stk(cavgs_stk, ncavgs, smpd, fail=.false.)
        !sz_cls2D = spproj%os_cls2D%get_noris()
        print *, "Number of cavgs ", ncavgs
        do icavgs = 1, ncavgs
            cnt = 0
            call spproj%os_ptcl2D%get_pinds(icavgs, 'class', pinds)
            if( .not. allocated(pinds) ) cycle
            nptcls = size(pinds)
            call cavg_img%new([params%box,params%box,1], smpd)
            call cavg_img%read(cavgs_stk,     icavgs)
            write(logfhandle,'(2(A,I3))')'>>> PROCESSING CLASS ', icavgs, ' NUMBER OF PARTICLES ', nptcls
            if( .not. allocated(pinds) ) cycle
            open(unit=25, file="ptcls_cc_analysis.csv")
            do iptcl = 1, nptcls
                cnt = cnt + 1
                ! compute cross-correlation between particle image and the class average image
                call prepimgbatch(nptcls)
                call read_imgbatch([1, nptcls])
                call discrete_read_imgbatch(nptcls, pinds(:), [1,nptcls] )
                corr = build%imgbatch(iptcl)%real_corr(cavg_img)
                write(25,'(2i6, 2F18.6)') icavgs, pinds(cnt), real(cnt)/real(nptcls), corr
            enddo
            write(25, '(A)') "          "
            call cavg_img%kill
        enddo
        command_plot = "gnuplot -e "//'"'//"set view map; set zrange[-.4:1]; splot " //"'"//"ptcls_cc_analysis.csv"//"'"// &
            " u 1:2:4 ; set term png; set xlabel " //"'"//"Class average"//"'"// "; set ylabel " //"'"//"Time"// &
            "'"// "; set title " //"'"//"Time Particles Class Analysis"//"'"// "; set nokey; set output 'ptcl_analysis.png'; replot" //'"'
        call execute_command_line(command_plot%to_char())
        call exec_cmdline('rm -rf fsc* fft* stderrout')
        ! deallocate
        call cavgs_stk%kill
        ! end gracefully
        call simple_end('**** PTCLSPROC_NANO NORMAL STOP ****')
    end subroutine exec_ptclsproc_nano

    subroutine exec_cavgseoproc_nano( self, cline )
        class(commander_cavgseoproc_nano), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters)     :: params
        type(image)          :: cavg_even, cavg_odd, img_w
        type(sp_project)     :: spproj
        type(string)         :: cavgs_stk, cavgs_stk_even, cavgs_stk_odd
        type(string)         :: command_plot
        real,    allocatable :: rstates(:), rad_cc(:,:), rad_dists(:,:)
        logical, allocatable :: state_mask(:)
        integer, allocatable :: pinds(:)
        integer :: ncavgs, i, j, cnt, tmax, tmin, tstamp
        real    :: smpd
        call cline%set('mkdir', 'yes') ! because we want to create the directory X_cavgseoproc_nano & copy the project file
        call params%new(cline)         ! because the parameters class manages directory creation and project file copying, mkdir = yes
        params%mkdir = 'no'            ! to prevent the input vol to be appended with ../
        call cline%set('mkdir', 'no')  ! because we do not want a nested directory structure in the execution directory
        ! read the project file
        call spproj%read(params%projfile)
        ! retrieve even and odd cavgs stack
        call spproj%get_cavgs_stk(cavgs_stk, ncavgs, smpd, fail=.false.)
        cavgs_stk_even = add2fbody(cavgs_stk,'.mrc','_even')
        cavgs_stk_odd  = add2fbody(cavgs_stk,'.mrc','_odd')
        if( allocated(rstates) ) deallocate(rstates)
        rstates = spproj%os_cls3D%get_all('state')
        ! compute radial cross-correlation between the even and odd cavgs
        print *,"params", params%box
        allocate(rad_cc(ncavgs,params%box/2), rad_dists(ncavgs,params%box/2), state_mask(ncavgs))
        cnt = 0
        open(unit=25, file="evenodd_radial_analysis.csv")
        print *, "ncavgs", ncavgs
        do i = 1, ncavgs
            cnt = cnt + 1
            if( rstates(cnt) > 0.5 )then
                print *,">>>PROCESSING CLASS ", i
                call cavg_odd%new([params%box,params%box,1], smpd)
                call cavg_odd%new([params%box,params%box,1], smpd)
                call cavg_even%new([params%box,params%box,1], smpd)
                call img_w%new([params%box,params%box,1], smpd)
                call cavg_odd%read(cavgs_stk_even,     i)
                call cavg_even%read(cavgs_stk_odd,     i)
                call spproj%os_ptcl2D%get_pinds(cnt, 'class', pinds)
                if( .not. allocated(pinds) ) cycle
                tmax   = maxval(pinds)
                tmin   = minval(pinds)
                tstamp = tmin + (tmax-tmin)/2
                call cavg_odd%radial_cc(cavg_even, img_w, smpd, rad_cc(cnt,:), rad_dists(cnt,:)) 
                do j = 1, size(rad_dists,dim=2)
                    write(25,'(i6, i6, 2F18.6, i6, i6)') cnt, tstamp, rad_dists(cnt,j), rad_cc(cnt,j), tmax-tmin, size(pinds)
                enddo
            else
                state_mask(cnt) = .false.
            endif
            write(25, '(A)') "          "
        enddo
        command_plot = "gnuplot -e "//'"'//"set pm3d map; set zrange[-.4:1]; splot " //"'"//"evenodd_radial_analysis.csv"//"'"// &
            " u 2:3:4 ; set term png; set xlabel " //"'"//"Time"//"'"// "; set ylabel " //"'"//"Radius({\305})"// &
            "'"// "; set title " //"'"//"Even-Odd Radial Cross-correlation"//"'"// "; set nokey; set output 'radial_analysis.png'; replot" //'"'
        call execute_command_line(command_plot%to_char())
        close(25)
        ! deallocate
        call cavgs_stk%kill
        ! end gracefully
        call simple_end('**** CAVGSEOPROC_NANO NORMAL STOP ****')
    end subroutine exec_cavgseoproc_nano

    subroutine exec_cavgsproc_nano( self, cline )
        class(commander_cavgsproc_nano), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)              :: params
        type(commander_refine3D_nano) :: xrefine3D_nano
        type(commander_reproject)     :: xreproject
        type(cmdline)                 :: cline_refine3D_cavgs, cline_reproject
        type(image),      allocatable :: imgs(:)
        type(image)                   :: img_w
        type(sp_project)              :: spproj
        type(string)                  :: cavgs_stk, command_plot, fname
        real,             allocatable :: rstates(:), rad_cc(:,:), rad_dists(:,:)
        logical,          allocatable :: state_mask(:)
        integer,          allocatable :: pinds(:)
        integer :: ncavgs, i, j, cnt, cnt2, tmax, tmin, tstamp
        real    :: smpd
        call cline%set('mkdir', 'yes') ! because we want to create the directory X_cavgsproc_nano & copy the project file
        call params%new(cline)         ! because the parameters class manages directory creation and project file copying, mkdir = yes
        params%mkdir = 'no'            ! to prevent the input vol to be appended with ../
        call cline%set('mkdir', 'no')  ! because we do not want a nested directory structure in the execution directory
        ! read the project file
        call spproj%read(params%projfile)
        ! retrieve cavgs stack
        call spproj%get_cavgs_stk(cavgs_stk, ncavgs, smpd, fail=.false.)
        if( ncavgs /= 0 )then
            ! update cline_refine3D_cavgs accordingly
            call cline_refine3D_cavgs%set('prg',      'refine3D_nano')
            call cline_refine3D_cavgs%set('vol1',      params%vols(1))
            call cline_refine3D_cavgs%set('pgrp',         params%pgrp)
            call cline_refine3D_cavgs%set('mskdiam',   params%mskdiam)
            call cline_refine3D_cavgs%set('nthr',         params%nthr)
            call cline_refine3D_cavgs%set('mkdir',               'no')
            call cline_refine3D_cavgs%set('maxits',                 1)
            call cline_refine3D_cavgs%set('projfile', params%projfile)
            call cline_refine3D_cavgs%set('oritype',          'cls3D')
            call cline_refine3D_cavgs%set('objfun',              'cc')
            call cline_refine3D_cavgs%set('lp',                   1.0)
            call cline_refine3D_cavgs%set('trs',                  5.0)
            call cline_refine3D_cavgs%set('nspace',             10000)
            call cline_refine3D_cavgs%set('center',              'no')
            ! convention for executing shared-memory workflows from within another workflow with a parameters object declared
            call xrefine3D_nano%execute_safe(cline_refine3D_cavgs)
            ! align cavgs
            call spproj%read_segment('cls3D', params%projfile) ! now the newly generated cls3D field will be read...
            ! ...so write out its content
            call spproj%os_cls3D%write(string('cavgs_oris.txt'))
            ! ...and get the state flags
            if( allocated(rstates) ) deallocate(rstates)
            rstates = spproj%os_cls3D%get_all('state')
            ! prepare for re-projection
            call cline_reproject%set('vol1',      params%vols(1))
            call cline_reproject%set('outstk',     'reprojs.mrc')
            call cline_reproject%set('smpd',                smpd)
            call cline_reproject%set('oritab',  'cavgs_oris.txt')
            call cline_reproject%set('pgrp',         params%pgrp)
            call cline_reproject%set('nthr',         params%nthr)
            call xreproject%execute_safe(cline_reproject)
            ! compute radial cross-correlation between cavgs and reproj
            allocate(rad_cc(ncavgs,params%box/2), rad_dists(ncavgs,params%box/2))
            ! write cavgs & reprojections
            allocate(imgs(2*ncavgs), state_mask(ncavgs))
            cnt = 0
            open(unit=25, file="radial_analysis.csv")
            do i = 1,2*ncavgs,2
                cnt = cnt + 1
                if( rstates(cnt) > 0.5 )then
                    call imgs(i    )%new([params%box,params%box,1], smpd)
                    call imgs(i + 1)%new([params%box,params%box,1], smpd)
                    call imgs(i    )%read(cavgs_stk,     cnt)
                    call imgs(i + 1)%read(string('reprojs.mrc'), cnt)
                    call img_w%new([params%box,params%box,1], smpd)
                    ! cavgs images are weighted using radial cross-correlation
                    call spproj%os_ptcl2D%get_pinds(cnt, 'class', pinds)
                    if( .not. allocated(pinds) ) cycle
                    tmax   = maxval(pinds)
                    tmin   = minval(pinds)
                    tstamp = tmin + (tmax-tmin)/2
                    call imgs(i)%radial_cc(imgs(i+1), img_w, smpd, rad_cc(cnt,:), rad_dists(cnt,:)) 
                    do j = 1, size(rad_dists,dim=2)
                        write(25,'(i6, i6, 2F18.6, i6, i6)') cnt, tstamp, rad_dists(cnt,j), rad_cc(cnt,j), tmax-tmin, size(pinds)
                    enddo
                    ! filter out cavgs
                    call imgs(i    )%norm
                    call imgs(i + 1)%norm
                    state_mask(cnt) = .true.
                else
                    state_mask(cnt) = .false.
                endif
                write(25, '(A)') "          "

            enddo
            command_plot = "gnuplot -e "//'"'//"set pm3d map; set zrange[-.4:1]; splot " //"'"//"radial_analysis.csv"//"'"// &
            " u 2:3:4 ; set term png; set xlabel " //"'"//"Time"//"'"// "; set ylabel " //"'"//"Radius({\305})"// &
            "'"// "; set title " //"'"//"Radial Cross-correlation"//"'"// "; set nokey; set output 'radial_analysis.png'; replot" //'"'
            call execute_command_line(command_plot%to_char())
            close(25)
            cnt   = 0
            cnt2  = 1 ! needed because we have omissions
            fname = 'cavgs_vs_reprojections.mrc'
            do i = 1,2*ncavgs,2
                cnt = cnt + 1
                if( state_mask(cnt) )then
                    call imgs(i    )%write(fname, cnt2    )
                    call imgs(i + 1)%write(fname, cnt2 + 1)
                    call imgs(i    )%kill
                    call imgs(i + 1)%kill
                    cnt2 = cnt2 + 2
                endif
            enddo
            deallocate(imgs)
        endif ! end of class average-based validation
        call exec_cmdline('rm -rf fsc* fft* recvol* RES* reprojs_recvol* reprojs* reproject_oris.txt stderrout')
        ! deallocate
        call cavgs_stk%kill
        ! end gracefully
        call simple_end('**** CAVGSPROC_NANO NORMAL STOP ****')
    end subroutine exec_cavgsproc_nano

    subroutine exec_graphene_subtr( self, cline )
        use simple_tseries_graphene_subtr
        class(commander_graphene_subtr), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)   :: params
        type(builder)      :: build
        type(image)        :: ave_pre, ave_post, img_tmp, img_tmp2, img_tmp3
        real,  allocatable :: angles1(:), angles2(:)
        real               :: smpd, ave,var,sdev
        integer            :: iptcl, ldim_ptcl(3), ldim(3), n, nptcls
        logical            :: err
        call cline%set('objfun','cc')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        ! sanity checks & dimensions
        call find_ldim_nptcls(params%stk,ldim_ptcl,nptcls,smpd=smpd)
        if( .not.cline%defined('smpd') )then
            if( smpd < 1.e-4 ) THROW_HARD('Please provide SMPD!')
            params%smpd = smpd
        endif
        ldim_ptcl(3) = 1
        call find_ldim_nptcls(params%stk2,ldim,n)
        ldim(3) = 1
        if( any(ldim-ldim_ptcl/=0) )THROW_HARD('Inconsistent dimensions between stacks')
        if( n /= nptcls )THROW_HARD('Inconsistent number of images between stacks')
        if( .not.cline%defined('outstk') )then
            params%outstk = add2fbody(basename(params%stk),params%ext,'_subtr')
        endif
        ! initialize subtracter
        allocate(angles1(nptcls),angles2(nptcls),source=0.)
        call init_graphene_subtr(params%box,params%smpd)
        ! init images
        call build%img%new(ldim,params%smpd)       ! particle
        call img_tmp%new(ldim,params%smpd)
        call img_tmp2%new(ldim,params%smpd)   ! nn background
        call img_tmp3%new(ldim,params%smpd)
        call ave_pre%new(ldim,params%smpd)
        call ave_post%new(ldim,params%smpd)
        ave_pre  = 0.
        ave_post = 0.
        ! read, subtract & write
        do iptcl = 1,nptcls
            if( mod(iptcl,50)==0 ) call progress(iptcl,nptcls)
            ! particle
            call build%img%read(params%stk,iptcl)
            call build%img%norm()
            ! neighbours background
            call img_tmp2%read(params%stk2,iptcl)
            ! detection
            call calc_peaks(img_tmp2, angles1(iptcl), angles2(iptcl))
            ! pre-subtraction average
            call img_tmp3%copy(build%img)
            call img_tmp3%zero_edgeavg
            call img_tmp3%fft()
            call img_tmp3%ft2img('sqrt', img_tmp)
            call ave_pre%add(img_tmp)
            ! subtraction
            call remove_lattices(build%img, angles1(iptcl), angles2(iptcl))
            call build%img%norm()
            call build%img%write(params%outstk, iptcl)
            ! graphene subtracted average
            call build%img%zero_edgeavg
            call build%img%fft()
            call build%img%ft2img('sqrt', img_tmp)
            call ave_post%add(img_tmp)
        enddo
        call progress(iptcl,nptcls)
        call ave_pre%div(real(nptcls))
        call ave_post%div(real(nptcls))
        call ave_pre%write(string('pre_subtr_ave_pspec.mrc'))
        call ave_post%write(string('subtr_ave_pspec.mrc'))
        ! stats
        call moment(angles1, ave, sdev, var, err)
        write(logfhandle,'(A,F6.2,A2,F6.2,A1)')'>>> POSITION GRAPHENE SHEET 1 (DEGREES): ',ave,' (',sdev,')'
        call moment(angles2, ave, sdev, var, err)
        write(logfhandle,'(A,F6.2,A2,F6.2,A1)')'>>> POSITION GRAPHENE SHEET 2 (DEGREES): ',ave,' (',sdev,')'
        angles1 = angles2 - angles1
        where(angles1>30.) angles1 = -(angles1 - 60.)
        call moment(angles1, ave, sdev, var, err)
        write(logfhandle,'(A,F6.2,A2,F6.2,A1)')'>>> RELATIVE ROTATION (DEGREES): ',ave,' (',sdev,')'
        ! cleanup
        call build%kill_general_tbox
        call kill_graphene_subtr
        call img_tmp%kill
        call img_tmp2%kill
        call img_tmp3%kill
        call ave_pre%kill
        call ave_post%kill
        ! end gracefully
        call simple_end('**** SIMPLE_GRAPHENE_SUBTR NORMAL STOP ****')
    end subroutine exec_graphene_subtr

    subroutine exec_denoise_trajectory( self, cline )
        use simple_commanders_imgproc, only: commander_ppca_denoise
        class(commander_denoise_trajectory), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(commander_ppca_denoise) :: xkpca_den
        if( .not. cline%defined('neigs')    ) call cline%set('neigs', 500)
        if( .not. cline%defined('pca_mode') ) call cline%set('pca_mode', 'kpca')
        call xkpca_den%execute(cline)
    end subroutine exec_denoise_trajectory

    subroutine exec_tseries_swap_stack( self, cline )
        use simple_commanders_project
        class(commander_tseries_swap_stack), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(sp_project) :: spproj, spproj_tmp
        type(parameters) :: params
        type(ctfparams)  :: ctfparms
        integer :: ldim(3), nimgs, nstks
        call cline%set('oritype','stk')
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call spproj%read(params%projfile)
        nstks = spproj%os_stk%get_noris()
        if( nstks < 1 ) THROW_HARD('No stack could be detected in the project!')
        call find_ldim_nptcls(params%stk, ldim, nimgs)
        ldim(3) = 1
        if( spproj%get_box() /= ldim(1) .or. spproj%get_box() /= ldim(2))then
            THROW_HARD('Incompatible dimensions between stacks')
        endif
        if( nimgs /= spproj%os_ptcl2D%get_noris() ) THROW_HARD('Incompatible number of images and orientation parameters!')
        ctfparms = spproj%get_ctfparams(params%oritype, 1)
        call spproj_tmp%read(params%projfile)
        call spproj%os_stk%kill
        call spproj%os_ptcl2D%kill
        call spproj%os_ptcl3D%kill
        call spproj%add_stk(params%stk, ctfparms)
        spproj%os_ptcl2D = spproj_tmp%os_ptcl2D
        spproj%os_ptcl3D = spproj_tmp%os_ptcl3D
        call spproj_tmp%kill
        if( nstks > 1 )then
            call spproj%os_ptcl2D%set_all2single('stkind',1)
            call spproj%os_ptcl3D%set_all2single('stkind',1)
        endif
        call spproj%write(params%projfile)
        call simple_end('**** SINGLE_TSERIES_SWAP_STACK NORMAL STOP ****')
    end subroutine exec_tseries_swap_stack

    subroutine exec_commander_tseries_reconstruct3D_distr( self, cline )
        use simple_commanders_rec, only: commander_volassemble
        real, parameter :: LP_LIST(4) = [1.5,2.0,2.5,3.0]
        real, parameter :: HP_LIM = 5.0 ! no information at lower res for these kind of data
        class(commander_tseries_reconstruct3D_distr), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(string),          allocatable :: vol_fnames(:)
        real,                  allocatable :: ccs(:,:,:), fsc(:), rstates(:), rad_cc(:), rad_dists(:)
        integer,               allocatable :: parts(:,:)
        type(string)                  :: recname, fname, str_state, recname_even, res_fname
        type(string)                  :: recname_odd, vol_fname_even, vol_fname_odd
        type(commander_reconstruct3D) :: xrec3D_shmem
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(cmdline)                 :: cline_rec
        type(image)                   :: vol1, vol2
        integer :: state, ipart, istate, nptcls, frame_start, frame_end
        integer :: funit, nparts, i, ind, nlps, ilp, iostat, hp_ind, lifetime
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('trs')     ) call cline%set('trs',           5.) ! to assure that shifts are being used
        if( .not. cline%defined('stepsz')  ) call cline%set('stepsz',      500.)
        if( .not. cline%defined('objfun')  ) call cline%set('objfun',      'cc') ! best objfun
        if( .not. cline%defined('ml_reg')  ) call cline%set('ml_reg',      'no') ! ml_reg=yes -> too few atoms 
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call cline%delete('refine')
        call params%new(cline)
        call spproj%read(params%projfile)
        if( cline%defined('fromp') .and. cline%defined('top') )then
            call cline%delete('nparts')   ! shared-memory implementation
            call cline%set('mkdir', 'no') ! to avoid nested directory structure
            call xrec3D_shmem%execute(cline)
            return
        endif
        ! state exception
        rstates = spproj%os_ptcl3D%get_all('state')
        if( any(rstates < 0.5) ) THROW_HARD('state=0 entries not allowed, prune project beforehand')
        ! states/stepz
        nptcls = size(rstates)
        if( cline%defined('nparts') )then
            nparts = params%nparts
        else
            nparts = nint(real(nptcls)/real(params%stepsz))
        endif
        parts = split_nobjs_even(nptcls, nparts)
        allocate(vol_fnames(nparts), rad_cc(params%box/2), rad_dists(params%box/2))
        recname      = VOL_FBODY//int2str_pad(1,2)//params%ext%to_char()
        recname_even = VOL_FBODY//int2str_pad(1,2)//'_even'//params%ext%to_char()
        recname_odd  = VOL_FBODY//int2str_pad(1,2)//'_odd'//params%ext%to_char()
        fname        = 'lifetimes.csv'
        call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=iostat)
        write(funit,*) 'PARTITION, ', 'FRAME_START, ', 'FRAME_END, ', 'LIFETIME'
        do ipart = 1,nparts
            str_state         = int2str_pad(ipart,2)
            vol_fnames(ipart) = 'partvol'//str_state%to_char()//params%ext%to_char()
            vol_fname_even    = 'partvol'//str_state%to_char()//'_even'//params%ext%to_char()
            vol_fname_odd     = 'partvol'//str_state%to_char()//'_odd'//params%ext%to_char()
            res_fname         = 'RESOLUTION_STATE'//str_state%to_char()
            ! prep 3D rec command line
            cline_rec = cline
            call cline_rec%delete('nparts') ! shared-memory implementation
            call cline_rec%set('fromp', parts(ipart,1))
            call cline_rec%set('top',   parts(ipart,2))
            frame_start = spproj%os_ptcl3D%get(parts(ipart,1), 'pind')
            frame_end   = spproj%os_ptcl3D%get(parts(ipart,2), 'pind')
            lifetime    = frame_end - frame_start + 1
            write(funit,'(I6,I6,I6,I6)') ipart, frame_start, frame_end, lifetime
            call cline_rec%set('mkdir', 'no')
            ! rec
            call xrec3D_shmem%execute(cline_rec)
            ! rename volumes and resolution files
            call simple_rename(recname,      vol_fnames(ipart))
            call simple_rename(recname_even, vol_fname_even)
            call simple_rename(recname_odd,  vol_fname_odd)
            if( ipart == 1 )then
                call simple_rename('RESOLUTION_STATE01', 'RESOLUTION_FILE_FIRST')
            else
                call simple_rename('RESOLUTION_STATE01',  res_fname)
            endif
        end do
        call fclose(funit)
        call simple_rename('RESOLUTION_FILE_FIRST', 'RESOLUTION_STATE01')
        ! Calculate correlation matrices
        nlps   = size(LP_LIST)
        hp_ind = calc_fourier_index(HP_LIM, params%box, params_glob%smpd)
        allocate(fsc(fdim(params%box)-1),ccs(nlps,nparts,nparts))
        call vol1%new([params%box,params%box,params%box],params%smpd)
        call vol2%new([params%box,params%box,params%box],params%smpd)
        ccs = 1.
        do state = 1, nparts - 1
            call vol1%zero_and_unflag_ft
            call vol1%read(vol_fnames(state))
            call vol1%fft
            do istate = state + 1, nparts
                call vol2%zero_and_unflag_ft
                call vol2%read(vol_fnames(istate))
                call vol2%fft
                call vol1%fsc(vol2,fsc)
                do ilp = 1, nlps
                    ind = calc_fourier_index(LP_LIST(ilp), params%box, params_glob%smpd)
                    ccs(ilp,state,istate) = sum(fsc(hp_ind:ind)) / real(ind - hp_ind + 1)
                    ccs(ilp,istate,state) = ccs(ilp,state,istate)
                enddo
            enddo
        enddo
        ! replace the diagonal elements with the maximum corr value in the columns
        ! for improved plotting / graphing
        do ilp = 1, nlps
            do istate = 1, nparts
                ccs(ilp,istate,istate) = maxval(ccs(ilp,istate,:), mask=ccs(ilp,istate,:) < 0.99)
            end do
        end do
        do ilp = 1,nlps
            fname = 'ccmat_lp'//trim(real2str(LP_LIST(ilp)))//'.csv'
            call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=iostat)
            do istate = 1, nparts
                do i = 1, nparts - 1
                    write(funit,'(F8.3,A2)', advance='no') ccs(ilp,istate,i), ', '
                enddo
                write(funit,'(F8.3)', advance='yes') ccs(ilp,istate,nparts)
            enddo
            call fclose(funit)
            fname = 'lp'//trim(real2str(LP_LIST(ilp)))//'.txt'
            fname = 'ccneigh_lp'//trim(real2str(LP_LIST(ilp)))//'.csv'
            call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=iostat)
            do istate = 1, nparts - 1
                write(funit,'(I3,A3,F8.3)',advance='yes') istate, ', ', ccs(ilp,istate,istate + 1)
            end do
            call fclose(funit)
        enddo
        call simple_end('**** SIMPLE_TSERIES_RECONSTRUCT3D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_commander_tseries_reconstruct3D_distr

    subroutine exec_tseries_core_finder( self, cline )
        use simple_opt_mask,  only: estimate_spher_mask
        use simple_image_msk, only: image_msk
        class(commander_tseries_core_finder), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline !< command line input
        real, parameter   :: RAD_LB = 5.0 ! 5.0 A radial boundary for averaging
        type(parameters)  :: params
        integer           :: ivol, nvols, ldim(3), ifoo, loc(1), best_msk, mskfromto(2)
        real, allocatable :: msk_in_pix(:)
        type(image_msk)   :: mskvol
        type(image)       :: vol, vol_ref, mskimg, mskimg2, mskimg2inv, vol_sum, vol_sum2, vol_w
        type(string), allocatable :: volnames(:)
        call params%new(cline)
        call read_filetable(params%filetab, volnames)
        nvols = size(volnames)
        if( params%mkdir.eq.'yes' )then
            do ivol = 1,nvols
                if(volnames(ivol)%to_char([1,1]).ne.'/') volnames(ivol) = '../'//volnames(ivol)%to_char()
            enddo
        endif
        call find_ldim_nptcls(volnames(1), ldim, ifoo)
        call vol%new(ldim, params%smpd)
        allocate(msk_in_pix(nvols), source=0.)
        do ivol = 1,nvols
            call vol%read(volnames(ivol))
            call mskvol%estimate_spher_mask_diam(vol, AMSKLP_NANO, msk_in_pix(ivol))
            write(logfhandle,*) ivol, 'mask diameter in A: ', 2. * msk_in_pix(ivol) * params%smpd
            call mskvol%kill
        end do
        loc = minloc(msk_in_pix) ! use the smallest NP to drive the radial averaging
        call vol_ref%new(ldim, params%smpd)
        call vol_ref%read(volnames(loc(1)))
        call vol_sum%copy(vol_ref)
        call vol_w%new(ldim, params%smpd)
        vol_w = 1.
        call mskimg%disc(ldim, params%smpd, msk_in_pix(loc(1)))
        mskfromto(1) = int(RAD_LB/params%smpd)
        mskfromto(2) = int(msk_in_pix(loc(1)))
        do ivol = 1,nvols
            if( ivol == loc(1) ) cycle
            call vol%zero_and_unflag_ft
            call vol%read(volnames(ivol))
            call estimate_spher_mask(vol_ref, vol, mskimg, mskfromto, best_msk)
            if( best_msk > mskfromto(1) )then ! agreement beyond the inputted lower radial limit
                call mskimg2%disc(ldim, params%smpd, real(best_msk))
                call vol%mul(mskimg2)
                call vol_sum%add(vol)
                call vol_w%add(mskimg2)
            else
                best_msk = 0
            endif
        enddo
        call vol_sum%div(vol_w)
        call vol_sum%write(string('radial_average_vol_')// int2str_pad(loc(1),2)//'.mrc')
        ! now do the reverse
        call vol_sum2%new(ldim, params%smpd)
        do ivol = 1,nvols
            vol_ref = vol_sum
            call vol%zero_and_unflag_ft
            call vol%read(volnames(ivol))
            call estimate_spher_mask(vol_ref, vol, mskimg, mskfromto, best_msk)
            call mskimg2%disc(ldim, params%smpd, real(best_msk))
            call mskimg2inv%copy(mskimg2)
            call mskimg2inv%bin_inv
            call vol_ref%mul(mskimg2)
            call vol%mul(mskimg2inv)
            call vol_sum2%zero_and_unflag_ft
            call vol_sum2%add(vol)
            call vol_sum2%add(vol_ref)
            call vol_sum2%write(string('core_inserted_vol_')// int2str_pad(ivol,2)//'.mrc')
        end do
    end subroutine exec_tseries_core_finder

    subroutine exec_tseries_make_projavgs( self, cline )
        use simple_strategy2D3D_common
        class(commander_tseries_make_projavgs), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        type(parameters)               :: params
        type(builder)                  :: build
        type(image),       allocatable :: pavgs(:), rot_imgs(:)
        real,              allocatable :: sumw(:), ref_weights(:,:)
        integer,           allocatable :: pinds(:), batches(:,:)
        logical,           allocatable :: ptcl_mask(:)
        real    :: euls_ref(3), euls(3), x, y, sdev_noise, dist, dist_threshold, w
        real    :: spiral_step
        integer :: nbatches, batchsz_max, batch_start, batch_end, batchsz, cnt
        integer :: iptcl, iref, ibatch, nptcls2update, i
        logical :: fall_over
        call cline%set('tseries', 'yes')
        call cline%set('pgrp',    'c1')
        call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('nspace') ) call cline%set('nspace',  300)
        if( .not. cline%defined('athres') ) call cline%set('athres',  10.)
        call build%init_params_and_build_strategy3D_tbox(cline, params)
        ! sanity check
        fall_over = .false.
        select case(trim(params%oritype))
            case('ptcl3D')
                fall_over = build%spproj%get_nptcls() == 0
            case DEFAULT
                THROW_HARD('unsupported ORITYPE')
        end select
        if( fall_over ) THROW_HARD('No images found!')
        call build%eulspace%set_all2single('state',1.)
        ! allocations
        allocate(pavgs(params%nspace),sumw(params%nspace))
        sumw = 0.0
        !$omp parallel do default(shared) private(iref) proc_bind(close) schedule(static)
        do iref = 1,params%nspace
            call pavgs(iref)%new([params%box,params%box,1],params%smpd,wthreads=.false.)
            call pavgs(iref)%zero_and_unflag_ft
        enddo
        !$omp end parallel do
        ! particles mask and indices
        allocate(ptcl_mask(params%nptcls),source=.true.)
        !$omp parallel do default(shared) private(iptcl) proc_bind(close) schedule(static)
        do iptcl = 1,params%nptcls
            ptcl_mask(iptcl) = build%spproj_field%get_state(iptcl) > 0
        enddo
        !$omp end parallel do
        nptcls2update = count(ptcl_mask)
        allocate(pinds(nptcls2update))
        i = 0
        do iptcl = 1,params%nptcls
            if( ptcl_mask(iptcl) )then
                i        = i+1
                pinds(i) = iptcl
            endif
        enddo
        ! batch prep
        batchsz_max = min(nptcls2update,params_glob%nthr*BATCHTHRSZ)
        nbatches    = ceiling(real(nptcls2update)/real(batchsz_max))
        batches     = split_nobjs_even(nptcls2update, nbatches)
        batchsz_max = maxval(batches(:,2)-batches(:,1)+1)
        call prepimgbatch(batchsz_max)
        allocate(ref_weights(batchsz_max,params%nspace),rot_imgs(batchsz_max))
        !$omp parallel do default(shared) private(iref) proc_bind(close) schedule(static)
        do i = 1,batchsz_max
            call rot_imgs(i)%new([params%box,params%box,1],params%smpd,wthreads=.false.)
        enddo
        !$omp end parallel do
        ! angular threshold
        euls_ref       = 0.
        euls           = [0.,params%athres,0.]
        dist_threshold = geodesic_frobdev(euls_ref,euls) / (2.*sqrt(2.))
        spiral_step    = rad2deg(3.809/sqrt(real(params%nspace)))
        do ibatch=1,nbatches
            call progress_gfortran(ibatch,nbatches)
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            ! Weights
            ref_weights = -1.0
            do i = 1,batchsz
                iptcl   = pinds(batch_start+i-1)
                euls    = build%spproj_field%get_euler(iptcl)
                euls(3) = 0.
                do iref = 1,params%nspace
                    if( build%eulspace%get_state(iref) == 0 ) cycle
                    euls_ref    = build%eulspace%get_euler(iref)
                    euls_ref(3) = 0.
                    dist        = geodesic_frobdev(euls_ref,euls) / (2.*sqrt(2.)) ! => [0;1]
                    if( dist > dist_threshold ) cycle
                    w = 1. / (1. + 999.0*dist)
                    ref_weights(i,iref) = w
                enddo
            enddo
            ! Images prep
            call discrete_read_imgbatch(batchsz, pinds(batch_start:batch_end), [1,batchsz])
            !$omp parallel do default(shared) private(i,iptcl,x,y,sdev_noise) proc_bind(close) schedule(static)
            do i = 1, batchsz
                iptcl = pinds(batch_start+i-1)
                call build%imgbatch(i)%norm_noise(build%lmsk, sdev_noise)
                x = build%spproj_field%get(iptcl, 'x')
                y = build%spproj_field%get(iptcl, 'y')
                call build%imgbatch(i)%fft
                call build%imgbatch(i)%shift2Dserial(-[x,y])
                call build%imgbatch(i)%ifft
                call rot_imgs(i)%zero_and_flag_ft
                call build%imgbatch(i)%rtsq(-build%spproj_field%e3get(iptcl),0.,0.,rot_imgs(i))
            enddo
            !$omp end parallel do
            ! Projection direction weighted sum
            !$omp parallel do default(shared) private(i,iref,w) proc_bind(close) schedule(static)
            do iref = 1,params%nspace
                if( build%eulspace%get_state(iref) == 0 ) cycle
                do i = 1,batchsz
                    w = ref_weights(i,iref)
                    if( w < TINY ) cycle
                    sumw(iref) = sumw(iref) + w
                    call pavgs(iref)%add(rot_imgs(i),w=w)
                enddo
            enddo
            !$omp end parallel do
        enddo
        ! Weights normalization
        !$omp parallel do default(shared) private(iref,w) proc_bind(close) schedule(static)
        do iref = 1,params%nspace
            if( build%eulspace%get_state(iref) == 0 ) cycle
            if( sumw(iref) > 0.001 )then
                call pavgs(iref)%div(sumw(iref))
            else
                call pavgs(iref)%zero_and_unflag_ft
            endif
        enddo
        !$omp end parallel do
        ! write
        cnt = 0
        do iref = 1,params%nspace
            if( build%eulspace%get_state(iref) == 0 ) cycle
            cnt = cnt + 1
            call pavgs(iref)%write(params%outstk,cnt)
        enddo
        call build%eulspace%write(string('projdirs.txt'))
        ! end gracefully
        call simple_end('**** SIMPLE_MAKE_PROJAVGS NORMAL STOP ****')
    end subroutine exec_tseries_make_projavgs

  end module simple_commanders_tseries
