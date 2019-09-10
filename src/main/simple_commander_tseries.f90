! concrete commander: time-series analysis
module simple_commander_tseries
include 'simple_lib.f08'
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_oris,           only: oris
use simple_parameters,     only: parameters
use simple_sp_project,     only: sp_project
use simple_qsys_funs
implicit none

public :: tseries_import_commander
public :: tseries_track_commander_distr
public :: tseries_track_commander
public :: cleanup2D_nano_commander
public :: tseries_average_commander
public :: tseries_ctf_estimate_commander
public :: refine3D_nano_commander_distr
public :: tseries_split_commander
public :: compare_nano_commander
public :: detect_atoms_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: tseries_import_commander
  contains
    procedure :: execute      => exec_tseries_import
end type tseries_import_commander
type, extends(commander_base) :: tseries_track_commander_distr
  contains
    procedure :: execute      => exec_tseries_track_distr
end type tseries_track_commander_distr
type, extends(commander_base) :: tseries_track_commander
  contains
    procedure :: execute      => exec_tseries_track
end type tseries_track_commander
type, extends(commander_base) :: cleanup2D_nano_commander
  contains
    procedure :: execute      => exec_cleanup2D_nano
end type cleanup2D_nano_commander
type, extends(commander_base) :: tseries_average_commander
  contains
    procedure :: execute      => exec_tseries_average
end type tseries_average_commander
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
type, extends(commander_base) :: tseries_split_commander
  contains
    procedure :: execute      => exec_tseries_split
end type tseries_split_commander
type, extends(commander_base) :: compare_nano_commander
  contains
    procedure :: execute      => exec_compare_nano
end type compare_nano_commander
type, extends(commander_base) :: detect_atoms_commander
  contains
    procedure :: execute      => exec_detect_atoms
end type detect_atoms_commander

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
        do iframe=1,nfiles
            call spproj%add_single_movie_frame(filenames(iframe), ctfvars)
        end do
        call spproj%write
        ! end gracefully
        call simple_end('**** tseries_import NORMAL STOP ****')
    end subroutine exec_tseries_import

    subroutine exec_tseries_track_distr( self, cline )
        use simple_nrtxtfile, only: nrtxtfile
        use simple_qsys_env,  only: qsys_env
        class(tseries_track_commander_distr), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters)              :: params
        type(qsys_env)                :: qenv
        type(chash)                   :: job_descr
        type(nrtxtfile)               :: boxfile
        real,        allocatable      :: boxdata(:,:)
        type(chash), allocatable      :: part_params(:)
        integer :: ndatlines, numlen, alloc_stat, j, orig_box, ipart
        call cline%set('nthr', 1.0)
        if( .not. cline%defined('neg')       ) call cline%set('neg',      'yes')
        if( .not. cline%defined('lp')        ) call cline%set('lp',         2.0)
        if( .not. cline%defined('width')     ) call cline%set('width',      1.1)
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      5.0)
        if( .not. cline%defined('nframesgrp')) call cline%set('nframesgrp', 10.)
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        if( .not. file_exists(params%boxfile) ) THROW_HARD('inputted boxfile does not exist in cwd')
        if( nlines(params%boxfile) > 0 )then
            call boxfile%new(params%boxfile, 1)
            ndatlines = boxfile%get_ndatalines()
            numlen    = len(int2str(ndatlines))
            allocate( boxdata(ndatlines,boxfile%get_nrecs_per_line()), stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('In: simple_commander_tseries :: exec_tseries_track', alloc_stat)
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
        call boxfile%kill
        call cline%delete('boxfile')
        params%nptcls = ndatlines
        params%nparts = params%nptcls
        ! box and numlen need to be part of command line
        call cline%set('box',    real(orig_box))
        call cline%set('numlen', real(numlen)  )
        ! prepare part-dependent parameters
        allocate(part_params(params%nparts))
        do ipart=1,params%nparts
            call part_params(ipart)%new(3)
            call part_params(ipart)%set('xcoord', real2str(boxdata(ipart,1)))
            call part_params(ipart)%set('ycoord', real2str(boxdata(ipart,2)))
            call part_params(ipart)%set('ind',    int2str(ipart))
        end do
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! schedule & clean
        call cline%gen_job_descr(job_descr)
        call qenv%gen_scripts_and_schedule_jobs( job_descr, part_params=part_params)
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_TRACK NORMAL STOP ****')
    end subroutine exec_tseries_track_distr

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
        integer :: i, iframe, ndatlines, j, orig_box, numlen, nframes
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
                    if( nint(boxdata(j,3)) /= nint(boxdata(j,4)) )THROW_HARD('Only square windows allowed!')
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
        nframes = spproj%get_nframes()
        allocate(framenames(nframes))
        iframe = 0
        do i = 1,spproj%os_mic%get_noris()
            if( spproj%os_mic%isthere(i,'frame') )then
                iframe = iframe + 1
                framenames(iframe) = trim(spproj%os_mic%get_static(i,'frame'))
            endif
        enddo
        ! actual tracking
        do j=1,ndatlines
            if( cline%defined('ind') )then
                dir = trim(params%fbody)//int2str(params%ind)
                call simple_mkdir(dir)
                if( .not. cline%defined('numlen') ) THROW_HARD('need numlen to be part of command line if ind is; exec_tseries_track')
            else
                dir = trim(params%fbody)//int2str(params%ind)
                call simple_mkdir(dir)
            endif
            call init_tracker( nint(boxdata(j,1:2)), framenames, dir, trim(params%fbody)//int2str_pad(j,numlen))
            call track_particle
            ! clean
            call kill_tracker
        end do
        call qsys_job_finished('simple_commander_tseries :: exec_tseries_track')
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_TRACK NORMAL STOP ****')
    end subroutine exec_tseries_track

    !> for distributed cleanup2D optimized for time-series of nanoparticles
    subroutine exec_cleanup2D_nano( self, cline )
        use simple_commander_cluster2D, only: make_cavgs_commander_distr,cluster2D_commander_distr
        class(cleanup2D_nano_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        ! commanders
        type(cluster2D_commander_distr) :: xcluster2D_distr
        ! other variables
        type(parameters)              :: params
        type(sp_project)              :: spproj
        character(len=:), allocatable :: orig_projfile
        character(len=LONGSTRLEN)     :: finalcavgs
        integer  :: nparts, last_iter_stage2
        call cline%set('center',    'yes')
        call cline%set('autoscale', 'no')
        call cline%set('refine',    'greedy')
        call cline%set('tseries',   'yes')
        if( .not. cline%defined('lp')      ) call cline%set('lp',     1.)
        if( .not. cline%defined('ncls')    ) call cline%set('ncls',   20.)
        if( .not. cline%defined('cenlp')   ) call cline%set('cenlp',  5.)
        if( .not. cline%defined('maxits')  ) call cline%set('maxits', 15.)
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        call params%new(cline)
        nparts = params%nparts
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! read project file
        call spproj%read(params%projfile)
        orig_projfile = trim(params%projfile)
        ! sanity checks
        if( spproj%get_nptcls() == 0 )then
            THROW_HARD('No particles found in project file: '//trim(params%projfile)//'; exec_cleanup2D_nano')
        endif
        ! delete any previous solution
        if( .not. spproj%is_virgin_field(params%oritype) )then
            ! removes previous cluster2D solution (states are preserved)
            call spproj%os_ptcl2D%delete_2Dclustering
            call spproj%write_segment_inside(params%oritype)
        endif
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
        ! transfer 2D shifts to 3D field
        call spproj%map2Dshifts23D
        call spproj%write
        call spproj%kill
        ! cleanup
        call del_file('start2Drefs'//params%ext)
        ! end gracefully
        call simple_end('**** SIMPLE_CLEANUP2D_NANO NORMAL STOP ****')
    end subroutine exec_cleanup2D_nano

    subroutine exec_tseries_average( self, cline )
        use simple_tseries_averager
        class(tseries_average_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(parameters) :: params
        if( .not. cline%defined('nframesgrp') ) call cline%set('nframesgrp',  10.)
        if( .not. cline%defined('corrw')      ) call cline%set('corrw', 'softmax')
        if( .not. cline%defined('rankw')      ) call cline%set('rankw',      'no')
        if( .not. cline%defined('outstk')     ) call cline%set('outstk', 'time_window_wavgs.mrcs')
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',     'yes')
        call params%new(cline)
        call init_tseries_averager
        call tseries_average
        call kill_tseries_averager
        call simple_end('**** SIMPLE_TSERIES_AVERAGE NORMAL STOP ****')
    end subroutine exec_tseries_average

    subroutine exec_tseries_backgr_subtr( self, cline )
        ! for background subtraction in time-series data. The goal is to subtract the two graphene
        ! peaks @ 2.14 A and @ 1.23 A. This is done by band-pass filtering the background image,
        ! recommended (and default settings) are hp=5.0 lp=1.1 and width=5.0.
        use simple_image, only: image
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
        use simple_image,             only: image
        use simple_ctf_estimate_fit,  only: ctf_estimate_fit
        class(tseries_ctf_estimate_commander), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        character(len=LONGSTRLEN), parameter :: pspec_fname  = 'tseries_ctf_estimate_pspec.mrc'
        character(len=LONGSTRLEN), parameter :: diag_fname   = 'tseries_ctf_estimate_diag'//JPG_EXT
        integer,                   parameter :: nmics4ctf    = 10
        type(parameters)              :: params
        type(builder)                 :: build
        type(ctf_estimate_fit)        :: ctffit
        type(image)                   :: mic_img, pspec_avg, pspec
        type(ctfparams)               :: ctfvars
        character(len=:), allocatable :: micname
        integer :: nmics, imic, dims(3), nselmics
        call cline%set('oritype','mic')
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('hp') )    call cline%set('hp', 10.)
        if( .not. cline%defined('lp') )    call cline%set('lp', 2.5)
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        ! generates ave
        dims(1) = nint(build%spproj_field%get(1,'xdim'))
        dims(2) = nint(build%spproj_field%get(1,'ydim'))
        dims(3) = 1
        call mic_img%new(dims, params%smpd)
        call pspec%new(    [params%pspecsz,params%pspecsz,1], params%smpd)
        call pspec_avg%new([params%pspecsz,params%pspecsz,1], params%smpd)
        nmics = build%spproj%get_nintgs()
        if( nmics /= build%spproj%os_mic%get_noris() ) THROW_HARD('Invalid number of micrographs')
        write(logfhandle,'(A)')'>>> BUILDING AVERAGE SPECTRUM'
        ! accumulate spectrum
        nselmics = 0
        do imic = 1,nmics,10
            call progress(imic,nmics)
            nselmics = nselmics + 1
            call build%spproj_field%getter(imic,'forctf',micname)
            call mic_img%read(micname)
            pspec = 0.
            call mic_img%mic2spec(params%pspecsz, 'sqrt', params%hp, pspec, postproc=.false.)
            call pspec_avg%add(pspec, w=0.01)
        enddo
        call progress(1,1)
        call mic_img%kill
        call pspec%kill
        ! average spectrum
        call pspec_avg%div(0.01*real(nselmics))
        ! post-process
        call pspec_avg%dampen_pspec_central_cross
        call pspec_avg%subtr_backgr(params%hp)
        ! ctf estimation
        write(logfhandle,'(A)')'>>> CTF ESTIMATION'
        ctfvars = build%spproj_field%get_ctfvars(1)
        call ctffit%new(pspec_avg, params%pspecsz, ctfvars, [params%dfmin,params%dfmax],&
            &[params%hp,params%lp],params%astigtol)
        call ctffit%fit(ctfvars, pspec_avg)
        call ctffit%write_diagnostic(diag_fname)
        do imic = 1,nmics
            call build%spproj_field%set_ctfvars(imic, ctfvars)
            call build%spproj_field%set(imic, 'ctf_estimate_cc', ctffit%get_ccfit())
            call build%spproj_field%set(imic, 'ctfscore', ctffit%get_ctfscore())
        enddo
        call ctffit%kill
        ! final write
        call build%spproj%write(params%projfile)
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
        call cline%set('graphene_filt', 'yes')
        call cline%set('ptclw',          'no')
        call cline%set('keepvol',       'yes')
        call cline%set('ninplpeaks',      1.0)
        call cline%set('wscheme',       'loc')
        ! dynamic parameters
        if( .not. cline%defined('nspace')      ) call cline%set('nspace',    10000.)
        if( .not. cline%defined('shcfrac')     ) call cline%set('shcfrac',      10.)
        if( .not. cline%defined('rankw')       ) call cline%set('rankw',      'exp')
        if( .not. cline%defined('trs')         ) call cline%set('trs',          2.0)
        if( .not. cline%defined('update_frac') ) call cline%set('update_frac',  0.2)
        if( .not. cline%defined('lp')          ) call cline%set('lp',            1.)
        if( .not. cline%defined('cenlp')       ) call cline%set('cenlp',         5.)
        if( .not. cline%defined('maxits')      ) call cline%set('maxits',       15.)
        if( .not. cline%defined('oritype')     ) call cline%set('oritype', 'ptcl3D')
        call xrefine3D_distr%execute(cline)
    end subroutine exec_refine3D_nano_distr

    subroutine exec_tseries_split( self, cline )
        use simple_ori, only: ori
        class(tseries_split_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(oris)       :: os
        type(ori)        :: o_tmp
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
                call build%spproj_field%get_ori(iptcl, o_tmp)
                call os%set_ori(cnt2, o_tmp )
                call build%img%read(params%stk, iptcl)
                call build%img%write(stkname,cnt2)
            end do
            call os%write(oriname)
            call os%kill
            deallocate(stkname, oriname)
        end do
        call o_tmp%kill
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_SPLIT NORMAL STOP ****')
    end subroutine exec_tseries_split

    ! for comparison of atomic models of 2nanoparticles
    subroutine exec_compare_nano( self, cline )
        use simple_nanoparticles_mod
        use simple_image, only : image
        class(compare_nano_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline !< command line input
        type(parameters)   :: params
        type(nanoparticle) :: nano1, nano2
        integer :: nptcls
        integer :: ldim1(3), ldim2(3)
        real    :: smpd1,smpd2
        call params%new(cline)
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_compare_nano')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! vol1 needs to be present; exec_compare_nano')
        endif
        if( .not. cline%defined('vol2') )then
            THROW_HARD('ERROR! vol2 needs to be present; exec_compare_nano')
        endif
        ! COMPARING
        call nano1%new(params%vols(1), params%smpd)
        call nano2%new(params%vols(2), params%smpd)
        call find_ldim_nptcls (params%vols(1), ldim1, nptcls, smpd1)
        call find_ldim_nptcls (params%vols(2), ldim2, nptcls, smpd2)
        if(any(ldim1 .ne. ldim2))   THROW_HARD('Non compatible dimensions of the particles to compare; compare_atomic_models')
        if(abs(smpd1-smpd2) > TINY) THROW_HARD('Non compatible particles, different smpd; compare_atomic_models')
        ! Nanoparticle binarization
        call nano1%compare_atomic_models(nano2)
        ! kill
        call nano1%kill
        call nano2%kill
        ! end gracefully
        call simple_end('**** SIMPLE_COMPARE_NANO NORMAL STOP ****')
    end subroutine exec_compare_nano

    ! for binarizing a nanoparticle and identiying its atomic positions
    subroutine exec_detect_atoms( self, cline )
        use simple_nanoparticles_mod
        use simple_image, only : image
        class(detect_atoms_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline !< command line input
        type(parameters)   :: params
        type(nanoparticle) :: nano
        integer :: ldim(3), nptcls
        real    :: smpd
        call params%new(cline)
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_detect_atoms')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! vol1 needs to be present; exec_detect_atoms')
        endif
        if( cline%defined('element') )then
            call nano%new(params%vols(1), params%smpd,params%element)
        else
            call nano%new(params%vols(1), params%smpd) !default pt
        endif
        call find_ldim_nptcls (params%vols(1), ldim, nptcls, smpd)
        ! execute
        call nano%detect_atoms()
        ! kill
        call nano%kill
        ! end gracefully
        call simple_end('**** SIMPLE_DETECT_ATOMS NORMAL STOP ****')
    end subroutine exec_detect_atoms

end module simple_commander_tseries
