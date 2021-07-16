! concrete commander: time-series analysis
module simple_commander_tseries
include 'simple_lib.f08'
use simple_builder,          only: builder
use simple_cmdline,          only: cmdline
use simple_commander_base,   only: commander_base
use simple_oris,             only: oris
use simple_parameters,       only: parameters, params_glob
use simple_sp_project,       only: sp_project
use simple_image,            only: image
use simple_binimage,         only: binimage
use simple_qsys_env,         only: qsys_env
use simple_commander_volops, only: reproject_commander
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
public :: validate_nano_commander
public :: tseries_swap_stack_commander
public :: tseries_reconstruct3D_distr
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
type, extends(commander_base) :: validate_nano_commander
  contains
    procedure :: execute      => exec_validate_nano
end type validate_nano_commander
type, extends(commander_base) :: tseries_swap_stack_commander
  contains
    procedure :: execute      => exec_tseries_swap_stack
end type tseries_swap_stack_commander
type, extends(commander_base) :: tseries_reconstruct3D_distr
  contains
    procedure :: execute      => exec_tseries_reconstruct3D_distr
end type tseries_reconstruct3D_distr

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
        if( .not. cline%defined('algorithm')  ) call cline%set('algorithm','patch')
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
        call o%new(is_ptcl=.false.)
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
        call simple_end('**** SIMPLE_TSERIES_TRACK_PARTICLES NORMAL STOP ****')
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
        character(len=LONGSTRLEN), allocatable :: intg_names(:), frame_names(:)
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
        allocate(intg_names(nframes),frame_names(nframes))
        iframe = 0
        do i = 1,spproj%os_mic%get_noris()
            if( spproj%os_mic%isthere(i,'frame') )then
                iframe = iframe + 1
                intg_names(iframe)  = trim(spproj%os_mic%get_static(i,'intg'))
                frame_names(iframe) = trim(spproj%os_mic%get_static(i,'frame'))
            endif
        enddo
        ! actual tracking
        dir = trim(params%fbody)
        call simple_mkdir(dir)
        call init_tracker( nint(boxdata(1,1:2)), intg_names, frame_names, dir, params%fbody)
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
        ! dynamic parameters
        if( .not. cline%defined('refine') )then
            call cline%set('refine','greedy')
        endif
        select case(trim(cline%get_carg('refine')))
        case('no','greedy')
            call cline%set('refine','greedy')
            if( .not. cline%defined('nptcls_per_cls') ) call cline%set('nptcls_per_cls', 35.)
            if( .not. cline%defined('maxits')         ) call cline%set('maxits',        15.0)
        case('inpl')
            call cline%set('center','no')
            if( .not. cline%defined('maxits')         ) call cline%set('maxits',         5.0)
        case DEFAULT
            THROW_HARD('Unsupported refinement mode!')
        end select
        if( .not. cline%defined('center')               ) call cline%set('center',       'yes')
        if( .not. cline%defined('graphene_filt')  ) call cline%set('graphene_filt','yes')
        if( .not. cline%defined('lpstart')        ) call cline%set('lpstart',        1.0)
        if( .not. cline%defined('lpstop')         ) call cline%set('lpstop',         1.0)
        if( .not. cline%defined('lp')             ) call cline%set('lp',             1.0)
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
        use simple_tseries_graphene_subtr
        class(graphene_subtr_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)   :: params
        type(builder)      :: build
        type(image)        :: ave_pre, ave_post, tmp
        real,  allocatable :: angles1(:), angles2(:)
        real               :: smpd, ave,var,sdev
        integer            :: iptcl, ldim_ptcl(3), ldim(3), n, nptcls
        logical            :: err
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
        call build%img_tmp%new(ldim,params%smpd)   ! nn background
        call build%img_copy%new(ldim,params%smpd)
        call tmp%new(ldim,params%smpd)
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
            call build%img_tmp%read(params%stk2,iptcl)
            ! detection
            call calc_peaks(build%img_tmp, angles1(iptcl), angles2(iptcl))
            ! pre-subtraction average
            call build%img_copy%copy(build%img)
            call build%img_copy%zero_edgeavg
            call build%img_copy%fft()
            call build%img_copy%ft2img('sqrt', tmp)
            call ave_pre%add(tmp)
            ! subtraction
            call remove_lattices(build%img, angles1(iptcl), angles2(iptcl))
            call build%img%norm()
            call build%img%write(params%outstk, iptcl)
            ! graphene subtracted average
            call build%img%zero_edgeavg
            call build%img%fft()
            call build%img%ft2img('sqrt', tmp)
            call ave_post%add(tmp)
        enddo
        call progress(iptcl,nptcls)
        call ave_pre%div(real(nptcls))
        call ave_post%div(real(nptcls))
        call ave_pre%write('pre_subtr_ave_pspec.mrc')
        call ave_post%write('subtr_ave_pspec.mrc')
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
        call tmp%kill
        call ave_pre%kill
        call ave_post%kill
        ! end gracefully
        call simple_end('**** SIMPLE_GRAPHENE_SUBTR NORMAL STOP ****')
    end subroutine exec_graphene_subtr

    subroutine exec_validate_nano( self, cline )
        use simple_commander_project
        class(validate_nano_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(new_project_commander)         :: xnew_project
        type(import_particles_commander)    :: ximport_ptcls
        type(refine3D_nano_commander_distr) :: xrefine3D
        type(reproject_commander)           :: xreproject
        type(sp_project)                    :: spproj
        character(len=:), allocatable       :: vol_fname
        real    :: smpd
        integer :: box
        call cline%set('mkdir', 'no')
        ! make new project called validation
        call cline%set('projname', 'validation')
        call xnew_project%execute(cline)
        ! import particles
        call cline%set('ctf', 'no')
        call cline%set('projfile', 'validation.simple')
        call ximport_ptcls%execute(cline)
        ! refine3D_nano
        params_glob => null()
        call cline%set('maxits', 5.)
        call xrefine3D%execute(cline)
        ! print_project_field
        call exec_cmdline('simple_exec prg=print_project_field oritype=ptcl3D > oris.txt')
        ! reproject
        params_glob => null()
        call spproj%read('validation.simple')
        call spproj%get_vol( 'vol', 1, vol_fname, smpd, box)
        call cline%delete('msk')
        call cline%set('vol1', vol_fname)
        call cline%set('outstk', 'reprojections.mrc')
        call cline%set('oritab', 'oris.txt')
        call xreproject%execute(cline)
        call exec_cmdline('rm -f frcs.bin fsc_state01.bin&
        & nohup.out oris.txt recvol_state01_even.spi recvol_state01_odd.spi&
        & *pproc* *part* RESOLUTION* *filelist* reproject_oris.txt')
    end subroutine exec_validate_nano

    subroutine exec_tseries_swap_stack( self, cline )
        use simple_commander_project
        class(tseries_swap_stack_commander), intent(inout) :: self
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
            call spproj%os_ptcl2D%set_all2single('stkind',1.)
            call spproj%os_ptcl3D%set_all2single('stkind',1.)
        endif
        call spproj%write(params%projfile)
        call simple_end('**** SINGLE_TSERIES_SWAP_STACK NORMAL STOP ****')
    end subroutine exec_tseries_swap_stack

    subroutine exec_tseries_reconstruct3D_distr( self, cline )
        use gnufor2
        real, parameter :: LP_LIST(4) = [1.5,2.0,2.5,3.0]
        real, parameter :: HP_LIM = 5.0 ! no information at lower res for these kind of data
        class(tseries_reconstruct3D_distr), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        character(len=LONGSTRLEN), allocatable :: list(:)
        character(len=:),          allocatable :: target_name
        character(len=STDLEN),     allocatable :: state_assemble_finished(:), vol_fnames(:)
        real,                      allocatable :: ccs(:,:,:), fsc(:)
        integer,                   allocatable :: parts(:,:)
        character(len=LONGSTRLEN) :: fname
        character(len=STDLEN)     :: volassemble_output, str_state, fsc_file, optlp_file
        type(parameters) :: params
        type(builder)    :: build
        type(qsys_env)   :: qenv
        type(cmdline)    :: cline_volassemble
        type(chash)      :: job_descr
        type(image)      :: vol1, vol2
        real             :: w, sumw
        integer          :: state, ipart, sz_list, istate, iptcl, cnt, nptcls, nptcls_per_state
        integer          :: funit, nparts, i, ind, nlps, ilp, iostat, hp_ind
        logical          :: fall_over
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('ptclw')   ) call cline%set('ptclw',       'no')
        if( .not. cline%defined('trs')     ) call cline%set('trs',           5.) ! to assure that shifts are being used
        if( .not. cline%defined('stepsz')  ) call cline%set('stepsz',      500.)
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call cline%delete('refine')
        call build%init_params_and_build_spproj(cline, params)
        ! sanity check
        fall_over = .false.
        select case(trim(params%oritype))
            case('ptcl3D')
                fall_over = build%spproj%get_nptcls() == 0
            case('cls3D')
                fall_over = build%spproj%os_out%get_noris() == 0
        case DEFAULT
            THROW_HARD('unsupported ORITYPE')
        end select
        if( fall_over ) THROW_HARD('No images found!')
        ! states/stepz
        nptcls = build%spproj_field%get_noris(consider_state=.true.)
        nparts = ceiling(real(nptcls)/real(params%stepsz))
        parts  = split_nobjs_even(nptcls, nparts)
        istate = 1
        cnt    = 0
        nptcls_per_state = parts(istate,2) - parts(istate,1) + 1
        do iptcl = 1,build%spproj_field%get_noris()
            if( build%spproj_field%get_state(iptcl) == 0 ) cycle
            cnt = cnt+1
            if( cnt > nptcls_per_state )then
                istate = istate + 1
                if( istate > nparts ) exit
                nptcls_per_state = parts(istate,2) - parts(istate,1) + 1
                cnt = 1
            endif
            call build%spproj_field%set(iptcl, 'state', real(istate))
            call build%spproj_field%set(iptcl, 'w',     1.)
        enddo
        call build%spproj%write_segment_inside(params%oritype)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        call cline%gen_job_descr(job_descr)
        call job_descr%set('prg','reconstruct3D')
        call job_descr%set('nstates',int2str(nparts))
        ! splitting
        if( trim(params%oritype).eq.'ptcl3D' )then
            call build%spproj%split_stk(params%nparts, dir=PATH_PARENT)
        endif
        ! eo partitioning
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype)
        endif
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr)
        ! assemble volumes
        allocate(state_assemble_finished(nparts), vol_fnames(nparts))
        do state = 1, nparts
            state_assemble_finished(state) = 'VOLASSEMBLE_FINISHED_STATE'//int2str_pad(state,2)
        enddo
        cline_volassemble = cline
        call cline_volassemble%set('prg', 'volassemble')
        do state = 1,nparts,params%nparts
            do istate = state,min(nparts,state+params%nparts-1)
                str_state = int2str_pad(istate,2)
                volassemble_output = 'RESOLUTION_STATE'//trim(str_state)
                call cline_volassemble%set( 'state', real(istate) )
                call cline_volassemble%set( 'nstates', real(nparts) )
                if( nparts>1 )call cline_volassemble%set('part', real(istate))
                call qenv%exec_simple_prg_in_queue_async(cline_volassemble,'simple_script_state'//trim(str_state), trim(volassemble_output))
                vol_fnames(istate) = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
            enddo
            call qsys_watcher(state_assemble_finished(state:min(nparts,state+params%nparts-1)))
            do istate = state,min(nparts,state+params%nparts-1)
                str_state = int2str_pad(istate,2)
                do ipart = 1,params%nparts
                    call del_file(trim(VOL_FBODY)//trim(str_state)//'_part'//trim(int2str(ipart))//'_even'//trim(params%ext))
                    call del_file(trim(VOL_FBODY)//trim(str_state)//'_part'//trim(int2str(ipart))//'_odd'//trim(params%ext))
                    call del_file('rho_'//trim(VOL_FBODY)//trim(str_state)//'_part'//trim(int2str(ipart))//'_even'//trim(params%ext))
                    call del_file('rho_'//trim(VOL_FBODY)//trim(str_state)//'_part'//trim(int2str(ipart))//'_odd'//trim(params%ext))
                enddo
            enddo
        end do
        ! Assemble states
        call vol1%new([params%box,params%box,params%box],params%smpd)
        call vol2%new([params%box,params%box,params%box],params%smpd)
        do state = 1,nparts
            call vol1%zero_and_unflag_ft
            sumw = 0.
            do istate = max(1,state-6),min(nparts,state+6)
                w     = exp(-real(istate-state)**2. / 16.0)
                sumw  = sumw + w
                call vol2%zero_and_unflag_ft
                call vol2%read(vol_fnames(istate))
                call vol1%add(vol2,w)
            enddo
            call vol1%div(sumw)
            call vol1%fft
            call vol1%bp(0.,params%lp)
            call vol1%ifft
            call vol1%mask(params%msk,'soft')
            call vol1%write('state_'// int2str_pad(state,2)//'.mrc')
        enddo
        ! Calculate correlation matrices
        nlps   = size(LP_LIST)
        hp_ind = calc_fourier_index(HP_LIM, params%box, params_glob%smpd)
        allocate(fsc(fdim(params%box)-1),ccs(nlps,nparts,nparts))
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
            call gnufor_image(ccs(ilp,:,:),palette='gray',filename=fname,persist='persist')
            fname = 'ccneigh_lp'//trim(real2str(LP_LIST(ilp)))//'.csv'
            call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=iostat)
            do istate = 1, nparts - 1
                write(funit,'(I3,A3,F8.3)',advance='yes') istate, ', ', ccs(ilp,istate,istate + 1)
            end do
            call fclose(funit)
        enddo
        ! termination
        call qsys_cleanup
        call build%spproj_field%kill
        call simple_end('**** SIMPLE_TSERIES_RECONSTRUCT3D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_tseries_reconstruct3D_distr

  end module simple_commander_tseries
