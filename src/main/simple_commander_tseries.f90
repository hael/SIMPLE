! concrete commander: time-series analysis
module simple_commander_tseries
include 'simple_lib.f08'
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_oris,           only: oris
use simple_parameters,     only: parameters
use simple_sp_project,     only: sp_project
use simple_image,          only: image
use simple_binimage,       only: binimage
use simple_qsys_env,       only: qsys_env
use simple_nanoparticle
use simple_qsys_funs
implicit none

public :: tseries_import_commander
public :: tseries_import_particles_commander
public :: tseries_make_pickavg_commander
public :: tseries_motion_correct_commander_distr
public :: tseries_motion_correct_commander
public :: tseries_track_particles_commander_distr
public :: tseries_track_particles_commander
public :: center2D_nano_commander_distr
public :: cluster2D_nano_commander_hlev
public :: tseries_ctf_estimate_commander
public :: refine3D_nano_commander_distr
public :: graphene_subtr_commander
public :: initial_3Dmodel_nano_commander_distr
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: tseries_import_commander
  contains
    procedure :: execute      => exec_tseries_import
end type tseries_import_commander
type, extends(commander_base) :: tseries_import_particles_commander
  contains
    procedure :: execute      => exec_tseries_import_particles
end type tseries_import_particles_commander
type, extends(commander_base) :: tseries_motion_correct_commander_distr
  contains
    procedure :: execute      => exec_tseries_motion_correct_distr
end type tseries_motion_correct_commander_distr
type, extends(commander_base) :: tseries_motion_correct_commander
  contains
    procedure :: execute      => exec_tseries_motion_correct
end type tseries_motion_correct_commander
type, extends(commander_base) :: tseries_make_pickavg_commander
  contains
    procedure :: execute      => exec_tseries_make_pickavg
end type tseries_make_pickavg_commander
type, extends(commander_base) :: tseries_track_particles_commander_distr
  contains
    procedure :: execute      => exec_tseries_track_particles_distr
end type tseries_track_particles_commander_distr
type, extends(commander_base) :: tseries_track_particles_commander
  contains
    procedure :: execute      => exec_tseries_track_particles
end type tseries_track_particles_commander
type, extends(commander_base) :: center2D_nano_commander_distr
  contains
    procedure :: execute      => exec_center2D_nano_distr
end type center2D_nano_commander_distr
type, extends(commander_base) :: cluster2D_nano_commander_hlev
  contains
    procedure :: execute      => exec_cluster2D_nano_hlev
end type cluster2D_nano_commander_hlev
type, extends(commander_base) :: tseries_backgr_subtr_commander
  contains
    procedure :: execute      => exec_tseries_backgr_subtr
end type tseries_backgr_subtr_commander
type, extends(commander_base) :: tseries_ctf_estimate_commander
  contains
    procedure :: execute      => exec_tseries_ctf_estimate
end type tseries_ctf_estimate_commander
type, extends(commander_base) :: refine3D_nano_commander_distr
  contains
    procedure :: execute      => exec_refine3D_nano_distr
end type refine3D_nano_commander_distr
type, extends(commander_base) :: graphene_subtr_commander
  contains
    procedure :: execute      => exec_graphene_subtr
end type graphene_subtr_commander
type, extends(commander_base) :: initial_3Dmodel_nano_commander_distr
  contains
    procedure :: execute      => exec_initial_3Dmodel_nano
end type initial_3Dmodel_nano_commander_distr
contains

    subroutine exec_tseries_import( self, cline )
        use simple_sp_project, only: sp_project
        class(tseries_import_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        type(ctfparams)  :: ctfvars
        character(len=LONGSTRLEN), allocatable :: filenames(:)
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
                if(filenames(iframe)(1:1).ne.'/') filenames(iframe) = '../'//trim(filenames(iframe))
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
        class(tseries_import_particles_commander), intent(inout) :: self
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
        call os%new(params%nptcls)
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
        call spproj%add_single_stk(params%stk, ctfvars, os)
        call spproj%write
        ! end gracefully
        call simple_end('**** TSERIES_IMPORT_PARTICLES NORMAL STOP ****')
    end subroutine exec_tseries_import_particles

    subroutine exec_tseries_motion_correct_distr( self, cline )
        class(tseries_motion_correct_commander_distr), intent(inout) :: self
        class(cmdline),                                intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        integer          :: nframes
        call cline%set('oritype',    'mic')
        call cline%set('groupframes', 'no')
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('nframesgrp') ) call cline%set('nframesgrp',    5.)
        if( .not. cline%defined('mcpatch')    ) call cline%set('mcpatch',    'yes')
        if( .not. cline%defined('nxpatch')    ) call cline%set('nxpatch',       3.)
        if( .not. cline%defined('nypatch')    ) call cline%set('nypatch',       3.)
        if( .not. cline%defined('trs')        ) call cline%set('trs',          10.)
        if( .not. cline%defined('lpstart')    ) call cline%set('lpstart',       5.)
        if( .not. cline%defined('lpstop')     ) call cline%set('lpstop',        3.)
        if( .not. cline%defined('bfac')       ) call cline%set('bfac',          5.)
        if( .not. cline%defined('wcrit')      ) call cline%set('wcrit',  'softmax')
        call params%new(cline)
        call cline%set('numlen', real(params%numlen))
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
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
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
        use simple_binoris_io,          only: binwrite_oritab
        use simple_motion_correct_iter, only: motion_correct_iter
        use simple_ori,                 only: ori
        class(tseries_motion_correct_commander), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        character(len=LONGSTRLEN), allocatable :: framenames(:)
        character(len=:),          allocatable :: frames2align
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
        call cline%set('groupframes', 'no')
        if( .not. cline%defined('nframesgrp') ) call cline%set('nframesgrp',    5.)
        if( .not. cline%defined('mcpatch')    ) call cline%set('mcpatch',    'yes')
        if( .not. cline%defined('nxpatch')    ) call cline%set('nxpatch',       3.)
        if( .not. cline%defined('nypatch')    ) call cline%set('nypatch',       3.)
        if( .not. cline%defined('trs')        ) call cline%set('trs',          10.)
        if( .not. cline%defined('lpstart')    ) call cline%set('lpstart',       5.)
        if( .not. cline%defined('lpstop')     ) call cline%set('lpstop',        3.)
        if( .not. cline%defined('bfac')       ) call cline%set('bfac',          5.)
        if( .not. cline%defined('wcrit')      ) call cline%set('wcrit',  'softmax')
        call params%new(cline)
        call spproj%read(params%projfile)
        nframes  = spproj%get_nframes()
        numlen_nframes = len(int2str(nframes))
        allocate(framenames(nframes))
        do i = 1,nframes
            if( spproj%os_mic%isthere(i,'frame') )then
                framenames(i) = trim(spproj%os_mic%get_static(i,'frame'))
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
                call mciter%iterate(cline_mcorr, ctfvars, o, 'tseries_win'//int2str_pad(iframe,numlen_nframes), frame_counter, frames2align, './', tseries='yes', gainref_fname=params%gainref)
            else
                call mciter%iterate(cline_mcorr, ctfvars, o, 'tseries_win'//int2str_pad(iframe,numlen_nframes), frame_counter, frames2align, './', tseries='yes')
            endif
            call spproj%os_mic%set_ori(iframe, o)
        end do
        ! output
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, [params%fromp,params%top], isegment=MIC_SEG)
        ! done!
        call del_file(frames2align)
        call img%kill
        call o%kill
        call qsys_job_finished('simple_commander_tseries :: exec_tseries_motion_correct' )
        call simple_end('**** SIMPLE_TSERIES_MOTION_CORRECT NORMAL STOP ****')
    end subroutine exec_tseries_motion_correct

    subroutine exec_tseries_make_pickavg( self, cline )
        use simple_commander_imgproc,   only: stack_commander
        use simple_motion_correct_iter, only: motion_correct_iter
        use simple_ori,                 only: ori
        use simple_tvfilter,            only: tvfilter
        class(tseries_make_pickavg_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        real, parameter :: LAM_TV = 1.5
        character(len=LONGSTRLEN), allocatable :: framenames(:)
        character(len=:),          allocatable :: filetabname
        type(sp_project)          :: spproj
        type(parameters)          :: params
        type(cmdline)             :: cline_stack, cline_mcorr
        type(stack_commander)     :: xstack
        type(motion_correct_iter) :: mciter
        type(ctfparams)           :: ctfvars
        type(ori)                 :: o
        type(tvfilter)            :: tvfilt
        type(image)               :: img_intg
        integer :: i, nframes, frame_counter, ldim(3), ifoo
        if( .not. cline%defined('nframesgrp') ) call cline%set('nframesgrp',    10.)
        if( .not. cline%defined('mcpatch')     ) call cline%set('mcpatch',    'yes')
        if( .not. cline%defined('nxpatch')    ) call cline%set('nxpatch',        3.)
        if( .not. cline%defined('nypatch')    ) call cline%set('nypatch',        3.)
        if( .not. cline%defined('trs')        ) call cline%set('trs',           10.)
        if( .not. cline%defined('lpstart')    ) call cline%set('lpstart',        5.)
        if( .not. cline%defined('lpstop')     ) call cline%set('lpstop',         3.)
        if( .not. cline%defined('bfac')       ) call cline%set('bfac',           5.)
        if( .not. cline%defined('groupframes')) call cline%set('groupframes',  'no')
        if( .not. cline%defined('wcrit')      ) call cline%set('wcrit',   'softmax')
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',       'yes')
        call cline%set('mcconvention','relion') ! ensures alignment to first frame
        call params%new(cline)
        call spproj%read(params%projfile)
        nframes  = spproj%get_nframes()
        allocate(framenames(params%nframesgrp))
        do i = 1,params%nframesgrp
            if( spproj%os_mic%isthere(i,'frame') )then
                framenames(i) = trim(spproj%os_mic%get_static(i,'frame'))
                params%smpd   = spproj%os_mic%get(i,'smpd')
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
        call cline_stack%set('nthr',    1.0)
        ! execute stack commander
        call xstack%execute(cline_stack)
        ! prepare 4 motion_correct
        cline_mcorr = cline
        call cline_mcorr%delete('nframesgrp')
        params%nframesgrp = 0
        call cline_mcorr%set('prg', 'motion_correct')
        call cline_mcorr%set('mkdir', 'no')
        call o%new
        ctfvars%smpd  = params%smpd
        frame_counter = 0
        ! motion corr
        if( cline%defined('gainref') )then
            call mciter%iterate(cline_mcorr, ctfvars, o, 'frames2align', frame_counter,&
                &'frames2align.mrc', './', gainref_fname=params%gainref, tseries='yes')
        else
            call mciter%iterate(cline_mcorr, ctfvars, o, 'frames2align', frame_counter,&
                &'frames2align.mrc', './', tseries='yes')
        endif
        call o%kill
        ! apply TV filter for de-noising
        call find_ldim_nptcls('frames2align_intg.mrc',ldim,ifoo)
        call img_intg%new(ldim, params%smpd)
        call img_intg%read('frames2align_intg.mrc',1)
        call tvfilt%new
        call tvfilt%apply_filter(img_intg, LAM_TV)
        call tvfilt%kill
        call img_intg%write('frames2align_intg_denoised.mrc')
        call img_intg%kill
        call simple_end('**** SIMPLE_GEN_INI_TSERIES_AVG NORMAL STOP ****')
    end subroutine exec_tseries_make_pickavg

    subroutine exec_tseries_track_particles_distr( self, cline )
        use simple_nrtxtfile, only: nrtxtfile
        class(tseries_track_particles_commander_distr), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters)              :: params
        type(qsys_env)                :: qenv
        type(chash)                   :: job_descr
        type(nrtxtfile)               :: boxfile
        real,        allocatable      :: boxdata(:,:)
        type(chash), allocatable      :: part_params(:)
        integer :: ndatlines, numlen, alloc_stat, j, orig_box, ipart
        if( .not. cline%defined('neg')       ) call cline%set('neg',      'yes')
        if( .not. cline%defined('lp')        ) call cline%set('lp',         2.3)
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      5.0)
        if( .not. cline%defined('nframesgrp')) call cline%set('nframesgrp', 30.)
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
            allocate( boxdata(ndatlines,boxfile%get_nrecs_per_line()), stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('In: simple_commander_tseries :: exec_tseries_track_particles', alloc_stat)
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
                call part_params(ipart)%set('fbody', trim(params%fbody)//'_'//int2str_pad(ipart,numlen))
            else
                call part_params(ipart)%set('fbody', trim(params%fbody))
            endif
        end do
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! schedule & clean
        call cline%gen_job_descr(job_descr)
        call qenv%gen_scripts_and_schedule_jobs( job_descr, part_params=part_params)
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_tseries_track_particles NORMAL STOP ****')
    end subroutine exec_tseries_track_particles_distr

    subroutine exec_tseries_track_particles( self, cline )
        use simple_tseries_track_particles
        use simple_qsys_funs,        only: qsys_job_finished
        use simple_sp_project,       only: sp_project
        class(tseries_track_particles_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(sp_project)                       :: spproj
        type(parameters)                       :: params
        character(len=:),          allocatable :: dir, forctf
        character(len=LONGSTRLEN), allocatable :: framenames(:)
        real,                      allocatable :: boxdata(:,:)
        integer :: i, iframe, orig_box, nframes
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
        allocate(framenames(nframes))
        iframe = 0
        do i = 1,spproj%os_mic%get_noris()
            if( spproj%os_mic%isthere(i,'frame') )then
                iframe = iframe + 1
                if( spproj%os_mic%isthere(i,'intg') )then
                    framenames(iframe) = trim(spproj%os_mic%get_static(i,'intg'))
                else
                    framenames(iframe) = trim(spproj%os_mic%get_static(i,'frame'))
                endif
            endif
        enddo
        ! actual tracking
        dir = trim(params%fbody)
        call simple_mkdir(dir)
        call init_tracker( nint(boxdata(1,1:2)), framenames, dir, params%fbody)
        call track_particle( forctf )
        ! clean tracker
        call kill_tracker
        ! end gracefully
        call qsys_job_finished('simple_commander_tseries :: exec_tseries_track_particles')
        call spproj%kill
        call simple_end('**** SIMPLE_TSERIES_TRACK_PARTICLES NORMAL STOP ****')
    end subroutine exec_tseries_track_particles

    subroutine exec_center2D_nano_distr( self, cline )
        use simple_commander_cluster2D, only: make_cavgs_commander_distr,cluster2D_commander_distr
        use simple_commander_imgproc,   only: pspec_int_rank_commander
        class(center2D_nano_commander_distr), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        ! commanders
        type(cluster2D_commander_distr) :: xcluster2D_distr
        type(pspec_int_rank_commander)  :: xpspec_rank
        ! other variables
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(cmdline)                 :: cline_pspec_rank
        character(len=:), allocatable :: orig_projfile
        character(len=STDLEN)         :: prev_ctfflag
        character(len=LONGSTRLEN)     :: finalcavgs
        integer  :: nparts, last_iter_stage2, nptcls
        call cline%set('dir_exec', 'center2D_nano')
        call cline%set('match_filt',          'no')
        call cline%set('ptclw',               'no')
        call cline%set('center',             'yes')
        call cline%set('autoscale',           'no')
        call cline%set('refine',          'greedy')
        call cline%set('tseries',            'yes')
        if( .not. cline%defined('graphene_filt') ) call cline%set('graphene_filt', 'yes')
        if( .not. cline%defined('lp')      )       call cline%set('lp',            1.)
        if( .not. cline%defined('ncls')    )       call cline%set('ncls',         20.)
        if( .not. cline%defined('cenlp')   )       call cline%set('cenlp',         5.)
        if( .not. cline%defined('trs')     )       call cline%set('trs',           5.)
        if( .not. cline%defined('maxits')  )       call cline%set('maxits',       15.)
        if( .not. cline%defined('oritype') )       call cline%set('oritype', 'ptcl2D')
        call params%new(cline)
        nparts = params%nparts
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! read project file
        call spproj%read(params%projfile)
        orig_projfile = trim(params%projfile)
        ! sanity checks
        nptcls = spproj%get_nptcls()
        if( nptcls == 0 )then
            THROW_HARD('No particles found in project file: '//trim(params%projfile)//'; exec_center2D_nano')
        endif
        ! delete any previous solution
        if( .not. spproj%is_virgin_field(params%oritype) )then
            ! removes previous cluster2D solution (states are preserved)
            call spproj%os_ptcl2D%delete_2Dclustering
            call spproj%write_segment_inside(params%oritype)
        endif
        ! de-activating CTF correction for this application
        prev_ctfflag = spproj%get_ctfflag(params%oritype)
        ! chronological initialisation
        params%nptcls_per_cls = ceiling(real(nptcls)/real(params%ncls))
        call cline%set('nptcls_per_cls', real(params%nptcls_per_cls))
        ! splitting
        call spproj%split_stk(params%nparts, dir=PATH_PARENT)
        ! no auto-scaling
        call cline%set('prg', 'cluster2D')
        call xcluster2D_distr%execute(cline)
        last_iter_stage2 = nint(cline%get_rarg('endit'))
        finalcavgs       = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//params%ext
        ! adding cavgs & FRCs to project
        params%projfile = trim(orig_projfile)
        call spproj%read( params%projfile )
        call spproj%add_frcs2os_out( trim(FRCS_FILE), 'frc2D')
        call spproj%add_cavgs2os_out(trim(finalcavgs), spproj%get_smpd(), imgkind='cavg')
        ! rank based on maximum of power spectrum
        call cline_pspec_rank%set('mkdir',   'no')
        call cline_pspec_rank%set('moldiam', params%moldiam)
        call cline_pspec_rank%set('nthr',    real(params%nthr))
        call cline_pspec_rank%set('smpd',    params%smpd)
        call cline_pspec_rank%set('stk',     finalcavgs)
        if( cline%defined('lp_backgr') ) call cline_pspec_rank%set('lp_backgr', params%lp_backgr)
        call xpspec_rank%execute(cline_pspec_rank)
        ! transfer 2D shifts to 3D field
        call spproj%map2Dshifts23D
        call spproj%write
        call spproj%kill
        ! cleanup
        call del_file('start2Drefs'//params%ext)
        ! end gracefully
        call simple_end('**** SIMPLE_CENTER2D_NANO NORMAL STOP ****')
    end subroutine exec_center2D_nano_distr

    subroutine exec_cluster2D_nano_hlev( self, cline )
        use simple_commander_cluster2D, only: cluster2D_autoscale_commander_hlev
        class(cluster2D_nano_commander_hlev), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        ! commander
        type(cluster2D_autoscale_commander_hlev) :: xcluster2D_distr
        ! static parameters
        call cline%set('prg',      'cluster2D')
        call cline%set('dir_exec', 'cluster2D_nano')
        call cline%set('match_filt',     'no')
        call cline%set('ptclw',          'no')
        call cline%set('center',        'yes')
        call cline%set('autoscale',      'no')
        call cline%set('tseries',       'yes')
        call cline%set('refine',     'greedy')
        ! dynamic parameters
        if( .not. cline%defined('graphene_filt')  ) call cline%set('graphene_filt','yes')
        if( .not. cline%defined('maxits')         ) call cline%set('maxits',        15.0)
        if( .not. cline%defined('lpstart')        ) call cline%set('lpstart',        1.0)
        if( .not. cline%defined('lpstop')         ) call cline%set('lpstop',         1.0)
        if( .not. cline%defined('lp')             ) call cline%set('lp',             1.0)
        if( .not. cline%defined('nptcls_per_cls') ) call cline%set('nptcls_per_cls', 35.)
        if( .not. cline%defined('winsz')          ) call cline%set('winsz',           3.)
        if( .not. cline%defined('cenlp')          ) call cline%set('cenlp',           5.)
        if( .not. cline%defined('trs')            ) call cline%set('trs',             5.)
        if( .not. cline%defined('oritype')        ) call cline%set('oritype',   'ptcl2D')
        call xcluster2D_distr%execute(cline)
        call simple_end('**** SIMPLE_CLUSTER2D_NANO NORMAL STOP ****')
    end subroutine exec_cluster2D_nano_hlev

    subroutine exec_tseries_backgr_subtr( self, cline )
        ! for background subtraction in time-series data. The goal is to subtract the two graphene
        ! peaks @ 2.14 A and @ 1.23 A. This is done by band-pass filtering the background image,
        ! recommended (and default settings) are hp=5.0 lp=1.1 and width=5.0.
        use simple_ctf,   only: ctf
        class(tseries_backgr_subtr_commander), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(image)      :: img_backgr, img_backgr_wctf, ave_img
        type(ctf)        :: tfun
        type(ctfparams)  :: ctfvars
        character(len=:), allocatable :: ext,imgname
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
            imgname = trim(get_fbody(params%outstk,ext,.true.))//'_ave.'//trim(ext)
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
        class(tseries_ctf_estimate_commander), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        character(len=LONGSTRLEN), parameter :: pspec_fname  = 'tseries_ctf_estimate_pspec.mrc'
        character(len=LONGSTRLEN), parameter :: diag_fname   = 'tseries_ctf_estimate_diag'//JPG_EXT
        integer,                   parameter :: nmics4ctf    = 10
        type(parameters)              :: params
        type(builder)                 :: build
        type(ctf_estimate_fit)        :: ctffit
        type(ctfparams)               :: ctfvars
        character(len=:), allocatable :: fname_diag, tmpl_fname, docname
        ! call cline%set('oritype','mic')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('hp')      ) call cline%set('hp', 5.)
        if( .not. cline%defined('lp')      ) call cline%set('lp', 1.)
        if( .not. cline%defined('dfmin')   ) call cline%set('dfmin', -0.05)
        if( .not. cline%defined('dfmax')   ) call cline%set('dfmax',  0.05)
        if( .not. cline%defined('astigtol')) call cline%set('astigtol', 0.001)
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        ! prep
        tmpl_fname = trim(get_fbody(basename(trim(params%stk)), params%ext, separator=.false.))
        fname_diag = filepath('./',trim(adjustl(tmpl_fname))//'_ctf_estimate_diag'//trim(JPG_EXT))
        docname    = filepath('./',trim(adjustl(tmpl_fname))//'_ctf'//trim(TXT_EXT))
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

    subroutine exec_refine3D_nano_distr( self, cline )
        use simple_commander_refine3D, only: refine3D_commander_distr
        class(refine3D_nano_commander_distr), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        ! commander
        type(refine3D_commander_distr) :: xrefine3D_distr
        ! static parameters
        call cline%set('prg',      'refine3D')
        call cline%set('match_filt',     'no')
        call cline%set('keepvol',       'yes')
        call cline%set('ninplpeaks',      1.0)
        ! dynamic parameters
        if( .not. cline%defined('graphene_filt') ) call cline%set('graphene_filt', 'yes')
        if( .not. cline%defined('ptclw')         ) call cline%set('ptclw',          'no')
        if( .not. cline%defined('nspace')        ) call cline%set('nspace',       10000.)
        if( .not. cline%defined('shcfrac')       ) call cline%set('shcfrac',         10.)
        if( .not. cline%defined('wcrit')         ) call cline%set('wcrit',         'inv')
        if( .not. cline%defined('trs')           ) call cline%set('trs',             5.0)
        if( .not. cline%defined('lp')            ) call cline%set('lp',              1.0)
        if( .not. cline%defined('cenlp')         ) call cline%set('cenlp',            5.)
        if( .not. cline%defined('maxits')        ) call cline%set('maxits',          15.)
        if( .not. cline%defined('oritype')       ) call cline%set('oritype',    'ptcl3D')
        call xrefine3D_distr%execute(cline)
    end subroutine exec_refine3D_nano_distr

    subroutine exec_graphene_subtr( self, cline )
        use simple_nano_utils, only: remove_graphene_peaks
        class(graphene_subtr_commander), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(image)      :: ave_pre, ave_post, tmp
        real             :: smpd, w, ang
        integer          :: iptcl, ldim_ptcl(3), ldim(3), n, nptcls
        call cline%set('oritype','mic')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        ! snaity checks
        call find_ldim_nptcls(params%stk,ldim_ptcl,nptcls,smpd=smpd)
        ldim_ptcl(3) = 1
        if( abs(smpd-params%smpd) > 1.e-4 )then
            THROW_HARD('Inconsistent sampling distances between project & input')
        endif
        call find_ldim_nptcls(params%stk2,ldim,n,smpd=smpd)
        ldim(3) = 1
        if( abs(smpd-params%smpd) > 1.e-4 )then
            THROW_HARD('Inconsistent sampling distances between project & input')
        endif
        if( any(ldim-ldim_ptcl/=0) )THROW_HARD('Inconsistent dimensions between stacks')
        if( n /= nptcls )THROW_HARD('Inconsistent number of images between stacks')
        ! init images
        ldim(3) = 1
        w = 1./real(nptcls)
        call build%img%new(ldim,smpd)       ! particle
        call build%img_tmp%new(ldim,smpd)   ! nn background
        call build%img_copy%new(ldim,smpd)
        call tmp%new(ldim,smpd)
        call ave_pre%new(ldim,smpd)
        call ave_post%new(ldim,smpd)
        ave_pre  = 0.
        ave_post = 0.
        ! read, subtract & write
        do iptcl = 1,nptcls
            call build%img%read(params%stk,iptcl)
            call build%img_tmp%read(params%stk2,iptcl)
            ! pre-subtraction average
            build%img_copy = build%img
            call build%img_copy%norm()
            call build%img_copy%zero_edgeavg
            call build%img_copy%fft()
            call build%img_copy%ft2img('sqrt', tmp)
            call ave_pre%add(tmp)
            ! subtraction
            call remove_graphene_peaks(build%img, build%img_tmp, ang)
            call build%img%write(params%outstk,iptcl)
            ! graphene subtracted average
            call build%img%norm()
            call build%img%zero_edgeavg
            call build%img%fft()
            call build%img%ft2img('sqrt', tmp)
            call ave_post%add(tmp)
        enddo
        call ave_pre%div(real(nptcls))
        call ave_post%div(real(nptcls))
        call ave_pre%write('pre_subtr_ave_pspec.mrc')
        call ave_post%write('subtr_ave_pspec.mrc')
        ! cleanup
        call build%kill_general_tbox
        call build%spproj%kill
        call tmp%kill
        call ave_pre%kill
        call ave_post%kill
        ! end gracefully
        call simple_end('**** SIMPLE_GRAPHENE_SUBTR NORMAL STOP ****')
    end subroutine exec_graphene_subtr

    !> for generation of an initial 3d model from class averages
    subroutine exec_initial_3Dmodel_nano( self, cline )
        use simple_commander_refine3D, only: refine3D_commander_distr
        use simple_commander_sim,      only: simulate_atoms_commander
        use simple_image,              only: image
        use simple_commander_volops,   only: reproject_commander
        use simple_parameters,         only: params_glob
        use simple_commander_quant,    only: detect_atoms_commander
        use simple_commander_rec,      only: reconstruct3D_commander_distr
        class(initial_3Dmodel_nano_commander_distr), intent(inout) :: self
        class(cmdline),                              intent(inout) :: cline
        ! constants
        real,                  parameter :: CENLP  = 5. !< consistency with refine3D
        integer,               parameter :: MAXITS = 15, NSPACE = 1000
        character(len=STDLEN), parameter :: WORK_PROJFILE      = 'initial_3Dmodel_tmpproj.simple'   ! temporary project
        ! initialization
        character(len=STDLEN), parameter :: START_SIM_VOL      = 'start_sim_vol.mrc'                ! starting volume
        ! stage 1
        character(len=STDLEN), parameter :: STAGE1_FBODY       = 'stage1'
        character(len=STDLEN), parameter :: STAGE1_VOL         = trim(STAGE1_FBODY)//'_vol'
        character(len=STDLEN), parameter :: STAGE1_VOL_PPROC   = trim(STAGE1_VOL)//trim(PPROC_SUFFIX)
        character(len=STDLEN), parameter :: STAGE1_ORIS        = trim(STAGE1_FBODY)//'_oris.txt'
        character(len=STDLEN), parameter :: STAGE1_PROJS       = trim(STAGE1_FBODY)//'_reprojs.mrc'
        character(len=STDLEN), parameter :: STAGE1_ATOMS_PDB   = trim(STAGE1_VOL)//'_atom_centers.pdb'
        character(len=STDLEN), parameter :: STAGE1_SIM_VOL     = trim(STAGE1_FBODY)//'sim_vol.mrc'
        character(len=STDLEN), parameter :: STAGE1_SIM_PROJS   = trim(STAGE1_FBODY)//'sim_reprojs.mrc'
        ! ranked reconstruction
        character(len=STDLEN), parameter :: RANKED_FBODY       = 'ranked'
        character(len=STDLEN), parameter :: RANKED_VOL         = trim(RANKED_FBODY)//'_vol.mrc'
        character(len=STDLEN), parameter :: RANKED_ATOMS_PDB   = trim(RANKED_FBODY)//'_vol_atom_centers.pdb'
        character(len=STDLEN), parameter :: RANKED_SIM_VOL     = trim(RANKED_FBODY)//'sim_vol.mrc'
        character(len=STDLEN), parameter :: RANKED_PROJS       = trim(RANKED_FBODY)//'sim_reprojs.mrc'
        ! stage 2
        character(len=STDLEN), parameter :: STAGE2_FBODY       = 'final'
        character(len=STDLEN), parameter :: STAGE2_VOL         = trim(STAGE2_FBODY)//'_vol'
        character(len=STDLEN), parameter :: STAGE2_VOL_PPROC   = trim(STAGE2_VOL)//trim(PPROC_SUFFIX)
        character(len=STDLEN), parameter :: STAGE2_ORIS        = trim(STAGE2_FBODY)//'_oris.txt'
        character(len=STDLEN), parameter :: STAGE2_PROJS       = trim(STAGE2_FBODY)//'_reprojs.mrc'
        character(len=STDLEN), parameter :: STAGE2_ATOMS_PDB   = trim(STAGE2_VOL)//'_atom_centers.pdb'
        character(len=STDLEN), parameter :: STAGE2_SIM_VOL     = trim(STAGE2_FBODY)//'_sim_vol.mrc'
        character(len=STDLEN), parameter :: STAGE2_SIM_PROJS   = trim(STAGE2_FBODY)//'_sim_reprojs.mrc'
        ! final re-projections
        character(len=STDLEN), parameter :: FINAL_REPROJS      = 'tseries_reprojs.mrc'
        ! distributed commanders
        type(refine3D_commander_distr)      :: xrefine3D_distr
        type(reconstruct3D_commander_distr) :: xreconstruct3D_distr
        ! shared-mem commanders
        type(reproject_commander)      :: xreproject
        type(simulate_atoms_commander) :: xsimatoms
        type(detect_atoms_commander)   :: xdetect_atoms
        ! command lines
        type(cmdline) :: cline_refine3D1, cline_refine3D2, cline_reconstruct3D
        type(cmdline) :: cline_detect_atoms, cline_simatoms, cline_reproject
        ! other variables
        character(len=:), allocatable :: stk
        integer,          allocatable :: states(:)
        real,             allocatable :: real_corrs(:), corrs_msk(:), tmp(:)
        logical,          allocatable :: lmsk(:,:,:), msk(:)
        character(len=2)              :: str_state
        type(parameters)      :: params
        type(ctfparams)       :: ctfvars ! ctf=yes by default
        type(sp_project)      :: spproj, work_proj
        type(image)           :: img, img2, img_msk
        character(len=STDLEN) :: vol_iter
        real                  :: iter, smpd, ave, sdev, var, corr_thresh
        integer               :: icls, ncavgs, status, cnt
        logical :: err
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',     'yes')
        ! static parameters
        call cline%set('match_filt', 'no')
        call cline%set('ninplpeaks',  1.0)
        call cline%set('objfun',     'cc')
        ! dynamic parameters
        if( .not. cline%defined('graphene_filt') ) call cline%set('graphene_filt', 'yes')
        if( .not. cline%defined('ptclw')         ) call cline%set('ptclw',          'no')
        if( .not. cline%defined('nspace')        ) call cline%set('nspace',        1000.)
        if( .not. cline%defined('wcrit')         ) call cline%set('wcrit',         'inv')
        if( .not. cline%defined('trs')           ) call cline%set('trs',             5.0)
        if( .not. cline%defined('lp')            ) call cline%set('lp',              1.5)
        if( .not. cline%defined('cenlp')         ) call cline%set('cenlp',         CENLP)
        if( .not. cline%defined('oritype')       ) call cline%set('oritype',    'ptcl3D')
        if( .not. cline%defined('element')       ) call cline%set('element',        'PT')
        ! hard set oritype
        call cline%set('oritype', 'out') ! because cavgs are part of out segment
        ! class averages, so no CTF!!
        ctfvars%ctfflag = CTFFLAG_NO
        ! make master parameters
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! from now on we are in the ptcl3D segment, final report is in the cls3D segment
        call cline%set('oritype', 'ptcl3D')
        ! state string
        str_state = int2str_pad(1,2)
        ! read project & update sampling distance
        call spproj%read(params%projfile)
        ! retrieve cavgs stack
        call spproj%get_cavgs_stk(stk, ncavgs, smpd)
        if( .not.spproj%os_cls2D%isthere('state') )then
            allocate(states(ncavgs), source=1) ! start from import
        else
            states = nint(spproj%os_cls2D%get_all('state')) ! start from previous 2D
        endif
        if( count(states==0) .eq. ncavgs )then
            THROW_HARD('no class averages detected in project file: '//trim(params%projfile)//'; initial_3Dmodel')
        endif
        ctfvars%smpd = params%smpd
        ! SANITY CHECKS
        ! prepare a temporary project file for the class average processing
        call del_file(WORK_PROJFILE)
        work_proj%projinfo = spproj%projinfo
        work_proj%compenv  = spproj%compenv
        if( spproj%jobproc%get_noris()  > 0 ) work_proj%jobproc = spproj%jobproc
        call work_proj%add_stk(trim(stk), ctfvars)
        call work_proj%os_ptcl3D%set_all('state', real(states)) ! takes care of states
        ! name change
        call work_proj%projinfo%delete_entry('projname')
        call work_proj%projinfo%delete_entry('projfile')
        call cline%set('projfile', trim(WORK_PROJFILE))
        call cline%set('projname', trim(get_fbody(trim(WORK_PROJFILE),trim('simple'))))
        call work_proj%update_projinfo(cline)
        params%box = work_proj%get_box()
        ! split
        if(params%nparts == 1 )then
            call work_proj%write()
        else
            call work_proj%split_stk(params%nparts)
        endif
        ! prepare command lines from prototype
        ! projects names are subject to change and updated individually
        call cline%delete('projname')
        call cline%delete('projfile')
        call cline%set('objfun', 'cc') ! always standard cross-correlation,
        cline_refine3D1     = cline
        cline_refine3D2     = cline
        cline_reconstruct3D = cline
        cline_reproject     = cline
        cline_simatoms      = cline
        ! reproject & simulate_atoms are not distributed executions
        call cline_reproject%delete('nparts')
        call cline_simatoms%delete('nparts')
        ! initialise command line parameters
        ! (1) INITIALIZATION, SIMULATED LATTICE
        call cline_simatoms%set('prg',    'simulate_atoms')
        call cline_simatoms%set('outvol', START_SIM_VOL)
        call cline_simatoms%set('smpd',   params%smpd)
        call cline_simatoms%set('box',    real(params%box))
        call cline_simatoms%delete('lp')
        ! STAGE 1
        ! (2) REFINE3D
        call cline_refine3D1%set('prg',      'refine3D')
        call cline_refine3D1%set('projfile', trim(WORK_PROJFILE))
        call cline_refine3D1%set('box',      real(params%box))
        call cline_refine3D1%set('maxits',   real(MAXITS))
        call cline_refine3D1%set('vol1',     START_SIM_VOL)
        call cline_refine3D1%set('refine',   'single')
        call cline_refine3D1%set('lp',       params%lp)
        ! (3) RE-PROJECT VOLUME, produces REC_PROJS
        call cline_reproject%set('prg',    'reproject')
        call cline_reproject%set('outstk', STAGE1_PROJS)
        call cline_reproject%set('smpd',   params%smpd)
        call cline_reproject%set('msk',    params%msk)
        call cline_reproject%set('box',    real(params%box))
        call cline_reproject%set('vol1',   trim(STAGE1_VOL)//params%ext)
        call cline_reproject%set('oritab', STAGE1_ORIS)
        call cline_reproject%delete('lp')
        ! (4) DETECT ATOMS, produces DETECTED_ATOMS_PDB
        call cline_detect_atoms%set('prg',      'detect_atoms')
        call cline_detect_atoms%set('smpd',     real(params%smpd))
        call cline_detect_atoms%set('cs_thres', 6.)
        call cline_detect_atoms%set('vol1',     trim(STAGE1_VOL)//params%ext)
        call cline_detect_atoms%set('nthr',     real(params%nthr))
        call cline_detect_atoms%set('msk',      params%msk)
        call cline_detect_atoms%set('element',  trim(params%element))
        call cline_detect_atoms%set('mkdir',    'no')
        call cline_detect_atoms%delete('lp')
        ! (5) SIMULATED DETECTED ATOMS, produces DETECTED_SIM_VOL
        ! (6) REPROJECTION OF SIMULATED DETECTED ATOMS
        ! RANKED STAGE
        ! (7) RECONSTRUCT WITH BEST RANKING IMAGES
        call cline_reconstruct3D%set('prg',      'reconstruct3D')
        call cline_reconstruct3D%set('msk',      params%msk)
        call cline_reconstruct3D%set('box',      real(params%box))
        call cline_reconstruct3D%set('projfile', trim(WORK_PROJFILE))
        ! (8) DETECT ATOMS
        ! (9) SIMULATED DETECTED ATOMS
        ! STAGE2
        ! (10) REFINE3D
        call cline_refine3D2%set('prg',      'refine3D')
        call cline_refine3D2%set('projfile', trim(WORK_PROJFILE))
        call cline_refine3D2%set('box',      real(params%box))
        call cline_refine3D2%set('maxits',   real(MAXITS))
        call cline_refine3D2%set('vol1',     RANKED_SIM_VOL)
        call cline_refine3D2%set('refine',  'greedy_single')
        ! (11) RE-PROJECTION OF FINAL VOLUME
        ! (12) DETECT ATOMS, produces RANKED_ATOMS_PDB
        ! (13) SIMULATE VOLUME
        ! (14) RE-PROJECTION OF FINAL SIMULATED VOLUME
        ! EXECUTIONS
        ! lattice simulation
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> SIMULATED SPHERICAL LATTICE'
        write(logfhandle,'(A)') '>>>'
        call xsimatoms%execute(cline_simatoms)
        call cline_simatoms%delete('vol1')
        call cline_simatoms%delete('element')
        call cline_simatoms%delete('moldiam')
        ! refine3D
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> INITIAL 3D MODEL GENERATION WITH REFINE3D - STAGE 1'
        write(logfhandle,'(A)') '>>>'
        call xrefine3D_distr%execute(cline_refine3D1)
        iter     = cline_refine3D1%get_rarg('endit')
        vol_iter = trim(VOL_FBODY)//trim(str_state)//params%ext
        call del_files(O_PEAKS_FBODY, params_glob%nparts, ext=BIN_EXT)
        ! deals with final volume
        status = simple_rename(vol_iter, trim(STAGE1_VOL)//params%ext)
        status = simple_rename(add2fbody(vol_iter,params%ext,PPROC_SUFFIX),trim(STAGE1_VOL_PPROC)//params%ext)
        call work_proj%read(WORK_PROJFILE)
        call work_proj%os_ptcl3D%write(STAGE1_ORIS)
        ! reprojections
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RE-PROJECTION OF STAGE 1 VOLUME'
        write(logfhandle,'(A)') '>>>'
        call xreproject%execute(cline_reproject)
        ! Atoms detection
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> ATOMS POSITIONS DETECTION WITH QUANT'
        write(logfhandle,'(A)') '>>>'
        call xdetect_atoms%execute(cline_detect_atoms)
        ! Simulation of detected atoms
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> SIMULATION OF DETECTED ATOMS'
        write(logfhandle,'(A)') '>>>'
        call cline_simatoms%set('pdbfile', STAGE1_ATOMS_PDB)
        call cline_simatoms%set('outvol',  STAGE1_SIM_VOL)
        call xsimatoms%execute(cline_simatoms)
        ! Re-projection of simulated detected atoms
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RE-PROJECTION OF SIMULATED DETECTED ATOMS'
        write(logfhandle,'(A)') '>>>'
        call cline_reproject%set('outstk', STAGE1_SIM_PROJS)
        call cline_reproject%set('vol1',   STAGE1_SIM_VOL)
        call cline_reproject%set('oritab', STAGE1_ORIS)
        call xreproject%execute(cline_reproject)
        ! write alternated stack
        allocate(real_corrs(ncavgs),source=0.)
        msk = states > 0
        call img%new([params%box,params%box,1], smpd)
        call img2%new([params%box,params%box,1],smpd)
        call img_msk%disc([params%box,params%box,1], params%smpd, params%msk, lmsk)
        call img_msk%kill
        cnt = -2
        do icls=1,ncavgs
            cnt = cnt + 3
            ! class-averages
            call img%read(stk, icls)
            call img%norm
            call img%write(FINAL_REPROJS, cnt)
            ! volume re-projections
            call img%read(STAGE1_PROJS, icls)
            call img%norm
            call img%write(FINAL_REPROJS, cnt+1)
            ! simulated atoms volume reprojections
            call img2%read(STAGE1_SIM_PROJS,icls)
            call img2%norm
            call img2%write(FINAL_REPROJS,cnt+2)
            ! real space correlations
            if( msk(icls) ) real_corrs(icls) = img%real_corr(img2, lmsk)
        enddo
        corrs_msk = pack(real_corrs, mask=msk)
        call moment(corrs_msk, ave, sdev, var, err)
        call otsu(corrs_msk, corr_thresh)
        corr_thresh = min(corr_thresh, max(minval(corrs_msk), ave-sdev))
        write(logfhandle,'(A)') '>>>'
        do icls=1,ncavgs
            if( states(icls) == 0 ) cycle
            if( real_corrs(icls) < corr_thresh )then
                states(icls) = 0 ! preserved in msk
                write(logfhandle,'(A,I4)') '>>> POTENTIALLY MISALIGNED CLASS AVERAGE: ',icls
            endif
        enddo
        ! Reconstruction with best ranking images
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RECONSTRUCTION FROM RANKED TIME SERIES AVERAGES'
        write(logfhandle,'(A)') '>>>'
        call work_proj%os_ptcl3D%set_all('state',real(states))
        call work_proj%write
        call xreconstruct3D_distr%execute(cline_reconstruct3D)
        status = simple_rename('recvol_state01.mrc', RANKED_VOL)
        where( msk ) states = 1
        ! Atoms detection
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> ATOMS POSITIONS DETECTION WITH QUANT'
        write(logfhandle,'(A)') '>>>'
        call cline_detect_atoms%set('vol1', RANKED_VOL)
        call xdetect_atoms%execute(cline_detect_atoms)
        ! Simulation
        call cline_simatoms%set('pdbfile', RANKED_ATOMS_PDB)
        call cline_simatoms%set('outvol',  RANKED_SIM_VOL)
        call xsimatoms%execute(cline_simatoms)
        ! final refine3D
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> INITIAL 3D MODEL GENERATION WITH REFINE3D - STAGE 2'
        write(logfhandle,'(A)') '>>>'
        call xrefine3D_distr%execute(cline_refine3D2)
        iter     = cline_refine3D2%get_rarg('endit')
        vol_iter = trim(VOL_FBODY)//trim(str_state)//params%ext
        call del_files(O_PEAKS_FBODY, params_glob%nparts, ext=BIN_EXT)
        ! deals with final volume
        status = simple_rename(vol_iter, trim(STAGE2_VOL)//params%ext)
        status = simple_rename(add2fbody(vol_iter,params%ext,PPROC_SUFFIX),trim(STAGE2_VOL_PPROC)//params%ext)
        ! updates original cls3D segment
        call work_proj%read_segment('ptcl3D', WORK_PROJFILE)
        params_glob%nptcls = ncavgs
        spproj%os_cls3D    = work_proj%os_ptcl3D
        call spproj%os_cls3D%delete_entry('eo')
        call spproj%os_cls3D%set_all2single('stkind',1.) ! revert splitting
        ! map the orientation parameters obtained for the clusters back to the particles
        call spproj%map2ptcls
        ! add rec_final to os_out
        call spproj%add_vol2os_out(trim(STAGE2_VOL)//params%ext, params%smpd, 1, 'vol_cavg')
        ! write results (this needs to be a full write as multiple segments are updated)
        call spproj%write()
        ! reprojections
        call spproj%os_cls3D%write(STAGE2_ORIS)
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RE-PROJECTION OF FINAL VOLUME'
        write(logfhandle,'(A)') '>>>'
        call cline_reproject%set('outstk', STAGE2_PROJS)
        call cline_reproject%set('vol1',   trim(STAGE2_VOL)//params%ext)
        call cline_reproject%set('oritab', STAGE2_ORIS)
        call xreproject%execute(cline_reproject)
        ! Atoms detection
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> ATOMS POSITIONS DETECTION WITH QUANT'
        write(logfhandle,'(A)') '>>>'
        call cline_detect_atoms%set('vol1', trim(STAGE2_VOL)//params%ext)
        call xdetect_atoms%execute(cline_detect_atoms)
        ! Simulation of detected atoms
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> SIMULATION OF DETECTED ATOMS'
        write(logfhandle,'(A)') '>>>'
        call cline_simatoms%set('pdbfile', STAGE2_ATOMS_PDB)
        call cline_simatoms%set('outvol',  STAGE2_SIM_VOL)
        call xsimatoms%execute(cline_simatoms)
        ! Re-projection of simulated detected atoms
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RE-PROJECTION OF SIMULATED DETECTED ATOMS'
        write(logfhandle,'(A)') '>>>'
        call cline_reproject%set('outstk', STAGE2_SIM_PROJS)
        call cline_reproject%set('vol1',   STAGE2_SIM_VOL)
        call cline_reproject%set('oritab', STAGE2_ORIS)
        call xreproject%execute(cline_reproject)
        ! write alternated stack
        call img_msk%disc([params%box,params%box,1], params%smpd, params%msk, lmsk)
        call img_msk%kill
        cnt = -2
        do icls=1,ncavgs
            cnt = cnt + 3
            ! volume re-projections
            call img%read(STAGE2_PROJS, icls)
            call img%norm
            call img%write(FINAL_REPROJS, cnt+1)
            ! simulated atoms volume reprojections
            call img2%read(STAGE2_SIM_PROJS,icls)
            call img2%norm
            call img2%write(FINAL_REPROJS,cnt+2)
        enddo
        ! end gracefully
        call img%kill
        call img2%kill
        call work_proj%kill
        call spproj%kill
        call del_file(WORK_PROJFILE)
        call simple_rmdir(STKPARTSDIR)
        call simple_end('**** SIMPLE_INITIAL_3DMODEL_NANO NORMAL STOP ****')
    end subroutine exec_initial_3Dmodel_nano

  end module simple_commander_tseries
