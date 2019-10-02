! concrete commander: cluster2D for simultanous 2D alignment and clustering of single-particle images
module simple_commander_cluster2D
include 'simple_lib.f08'
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: parameters, params_glob
use simple_sp_project,     only: sp_project
use simple_qsys_env,       only: qsys_env
use simple_image,          only: image
use simple_qsys_funs
implicit none

public :: cleanup2D_commander
public :: cluster2D_commander_stream
public :: cluster2D_autoscale_commander
public :: cluster2D_commander_distr
public :: cluster2D_commander
public :: make_cavgs_commander_distr
public :: make_cavgs_commander
public :: cavgassemble_commander
public :: check_2Dconv_commander
public :: rank_cavgs_commander
public :: cluster_cavgs_commander
public :: write_classes_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: cleanup2D_commander
  contains
    procedure :: execute      => exec_cleanup2D
end type cleanup2D_commander
type, extends(commander_base) :: cluster2D_commander_stream
  contains
    procedure :: execute      => exec_cluster2D_stream
end type cluster2D_commander_stream
type, extends(commander_base) :: cluster2D_autoscale_commander
  contains
    procedure :: execute      => exec_cluster2D_autoscale
end type cluster2D_autoscale_commander
type, extends(commander_base) :: cluster2D_commander_distr
  contains
    procedure :: execute      => exec_cluster2D_distr
end type cluster2D_commander_distr
type, extends(commander_base) :: cluster2D_commander
  contains
    procedure :: execute      => exec_cluster2D
end type cluster2D_commander
type, extends(commander_base) :: make_cavgs_commander_distr
  contains
    procedure :: execute      => exec_make_cavgs_distr
end type make_cavgs_commander_distr
type, extends(commander_base) :: make_cavgs_commander
  contains
    procedure :: execute      => exec_make_cavgs
end type make_cavgs_commander
type, extends(commander_base) :: cavgassemble_commander
  contains
    procedure :: execute      => exec_cavgassemble
end type cavgassemble_commander
type, extends(commander_base) :: check_2Dconv_commander
  contains
    procedure :: execute      => exec_check_2Dconv
end type check_2Dconv_commander
type, extends(commander_base) :: rank_cavgs_commander
  contains
    procedure :: execute      => exec_rank_cavgs
end type rank_cavgs_commander
type, extends(commander_base) :: cluster_cavgs_commander
  contains
    procedure :: execute      => exec_cluster_cavgs
end type cluster_cavgs_commander
type, extends(commander_base) :: write_classes_commander
  contains
    procedure :: execute      => exec_write_classes
end type write_classes_commander

contains

    subroutine exec_make_cavgs_distr( self, cline )
        class(make_cavgs_commander_distr), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters) :: params
        type(cmdline)    :: cline_cavgassemble
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! prepare command lines from prototype master
        cline_cavgassemble = cline
        call cline_cavgassemble%set('prg', 'cavgassemble')
        call cline_cavgassemble%set('nthr', 0.) ! to ensure the use of all resources in assembly
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr)
        ! assemble class averages
        call qenv%exec_simple_prg_in_queue(cline_cavgassemble, 'CAVGASSEMBLE_FINISHED')
        call qsys_cleanup
        call simple_end('**** SIMPLE_DISTR_MAKE_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_make_cavgs_distr

    subroutine exec_make_cavgs( self, cline )
        use simple_classaverager
        class(make_cavgs_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer :: ncls_here
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_strategy2D_tbox(cline, params)
        write(logfhandle,'(a)') '>>> GENERATING CLUSTER CENTERS'
        ! deal with the orientations
        ncls_here = build%spproj_field%get_n('class')
        if( .not. cline%defined('ncls') ) params%ncls = build%spproj_field%get_n('class')
        if( params%l_remap_cls )then
            call build%spproj_field%remap_cls()
            if( cline%defined('ncls') )then
                if( params%ncls < ncls_here ) THROW_HARD('inputted ncls < ncls_in_oritab not allowed!')
                if( params%ncls > ncls_here )then
                    call build%spproj_field%expand_classes(params%ncls)
                endif
            endif
        else if( params%tseries .eq. 'yes' )then
            if( .not. cline%defined('ncls') )then
                THROW_HARD('# class averages (ncls) need to be part of command line when tseries=yes')
            endif
            call build%spproj_field%ini_tseries(params%ncls, 'class')
            call build%spproj_field%partition_eo(tseries=.true.)
        endif
        ! shift multiplication
        if( params%mul > 1. )then
            call build%spproj_field%mul_shifts(params%mul)
        endif
        ! setup weights
        if( params%l_ptclw )then
            call build%spproj_field%calc_soft_weights2D
        else
            call build%spproj_field%calc_hard_weights2D(params%frac, params%ncls)
        endif
        ! even/odd partitioning
        if( build%spproj_field%get_nevenodd() == 0 ) call build%spproj_field%partition_eo
        ! write
        if( nint(cline%get_rarg('part')) .eq. 1 )then
            call build%spproj%write_segment_inside(params%oritype)
        endif
        ! create class averager
        call cavger_new('class')
        ! transfer ori data to object
        call cavger_transf_oridat(build%spproj)
        ! standard cavg assembly
        call cavger_assemble_sums( .false. )
        ! write sums
        call cavger_readwrite_partial_sums('write')
        call qsys_job_finished(  'simple_commander_cluster2D :: exec_make_cavgs' )
        call cavger_kill
        ! end gracefully
        call build%kill_strategy2D_tbox
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_MAKE_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_make_cavgs

    subroutine exec_cleanup2D( self, cline )
        use simple_commander_project, only: scale_project_commander_distr
        use simple_procimgfile,       only: random_selection_from_imgfile, random_cls_from_imgfile
        class(cleanup2D_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        ! commanders
        type(cluster2D_commander_distr)     :: xcluster2D_distr
        type(scale_project_commander_distr) :: xscale_distr
        type(rank_cavgs_commander)          :: xrank_cavgs
        ! command lines
        type(cmdline)                       :: cline_cluster2D1, cline_cluster2D2
        type(cmdline)                       :: cline_rank_cavgs, cline_scale
        ! other variables
        type(parameters)                    :: params
        type(sp_project)                    :: spproj, spproj_sc
        character(len=:),       allocatable :: projfile, orig_projfile
        character(len=LONGSTRLEN)           :: finalcavgs, finalcavgs_ranked, cavgs
        real                                :: scale_factor, smpd, msk, ring2, lp1, lp2
        integer                             :: last_iter, box, status
        logical                             :: do_scaling
        ! parameters
        character(len=STDLEN) :: orig_projfile_bak = 'orig_bak.simple'
        integer, parameter    :: MINBOX      = 92
        real,    parameter    :: TARGET_LP   = 15.
        real,    parameter    :: MINITS      = 5.
        real,    parameter    :: MAXITS      = 15.
        real                  :: SMPD_TARGET = 4.
        if( .not. cline%defined('lp')        ) call cline%set('lp',         15. )
        if( .not. cline%defined('ncls')      ) call cline%set('ncls',      200. )
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      20. )
        if( .not. cline%defined('center')    ) call cline%set('center',     'no')
        if( .not. cline%defined('maxits')    ) call cline%set('maxits',     15. )
        if( .not. cline%defined('center')    ) call cline%set('center',    'no' )
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('oritype')   ) call cline%set('oritype', 'ptcl2D')
        call params%new(cline)
        orig_projfile = params%projfile
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! read project file
        call spproj%read(params%projfile)
        ! sanity checks
        if( spproj%get_nptcls() == 0 )then
            THROW_HARD('No particles found in project file: '//trim(params%projfile)//'; exec_cleanup2D_autoscale')
        endif
        ! delete any previous solution
        call spproj%os_ptcl2D%delete_2Dclustering
        call spproj%write_segment_inside(params%oritype)
        ! splitting
        call spproj%split_stk(params%nparts, dir=PATH_PARENT)
        ! first stage
        ! down-scaling for fast execution, greedy optimisation, no match filter, bi-linear interpolation,
        ! no incremental learning, objective function is standard cross-correlation (cc), no centering
        cline_cluster2D1 = cline
        cline_cluster2D2 = cline
        cline_scale      = cline
        call cline_cluster2D1%set('prg',        'cluster2D')
        call cline_cluster2D1%set('refine',     'greedy')
        call cline_cluster2D1%set('maxits',     MINITS)
        call cline_cluster2D1%set('objfun',     'cc')
        call cline_cluster2D1%set('match_filt', 'no')
        call cline_cluster2D1%set('center',     'no')
        call cline_cluster2D1%set('autoscale',  'no')
        call cline_cluster2D1%set('ptclw',      'no')
        call cline_cluster2D1%delete('update_frac')
        ! second stage
        ! down-scaling for fast execution, greedy optimisation, no match filter, bi-linear interpolation,
        ! objective function default is standard cross-correlation (cc)
        call cline_cluster2D2%set('prg',        'cluster2D')
        call cline_cluster2D2%set('refine',     'greedy')
        call cline_cluster2D2%set('match_filt', 'no')
        call cline_cluster2D2%set('autoscale',  'no')
        call cline_cluster2D2%set('trs',         MINSHIFT)
        call cline_cluster2D2%set('objfun',     'cc')
        if( .not.cline%defined('maxits') ) call cline_cluster2D2%set('maxits', MAXITS)
        if( cline%defined('update_frac') )call cline_cluster2D2%set('update_frac',params%update_frac)
        ! Scaling
        do_scaling = .true.
        if( params%box < MINBOX .or. params%autoscale.eq.'no')then
            do_scaling   = .false.
            smpd         = params%smpd
            scale_factor = 1.
            box          = params%box
            projfile     = trim(params%projfile)
        else
            call autoscale(params%box, params%smpd, SMPD_TARGET, box, smpd, scale_factor)
            if( box < MINBOX ) SMPD_TARGET = params%smpd * real(params%box) / real(MINBOX)
            call spproj%scale_projfile(SMPD_TARGET, projfile, cline_cluster2D1, cline_scale, dir=trim(STKPARTSDIR))
            call spproj%kill
            scale_factor = cline_scale%get_rarg('scale')
            smpd         = cline_scale%get_rarg('smpd')
            box          = nint(cline_scale%get_rarg('newbox'))
            call cline_scale%set('state',1.)
            call cline_scale%delete('smpd') !!
            call simple_mkdir(trim(STKPARTSDIR),errmsg="commander_hlev_wflows :: exec_cluster2D_autoscale;  ")
            call xscale_distr%execute( cline_scale )
            ! rename scaled projfile and stash original project file
            ! such that the scaled project file has the same name as the original and can be followed from the GUI
            call simple_copy_file(orig_projfile, orig_projfile_bak)
            call spproj%read_non_data_segments(projfile)
            call spproj%projinfo%set(1,'projname',get_fbody(orig_projfile,METADATA_EXT,separator=.false.))
            call spproj%projinfo%set(1,'projfile',orig_projfile)
            call spproj%write_non_data_segments(projfile)
            call spproj%kill
            status   = simple_rename(projfile,orig_projfile)
            projfile = trim(orig_projfile)
        endif
        if( cline%defined('msk') )then
            msk = params%msk*scale_factor
        else
            msk = real(box/2)-COSMSKHALFWIDTH
        endif
        ring2 = 0.8*msk
        lp1   = max(2.*smpd, max(params%lp,TARGET_LP))
        lp2   = max(2.*smpd, params%lp)
        ! execute initialiser
        params%refs = 'start2Drefs' // params%ext
        call spproj%read(projfile)
        if( params%avg.eq.'yes' )then
            call random_cls_from_imgfile(spproj, params%refs, params%ncls)
        else
            call random_selection_from_imgfile(spproj, params%refs, box, params%ncls)
        endif
        call spproj%kill
        ! updates command-lines
        call cline_cluster2D1%set('refs',   params%refs)
        call cline_cluster2D1%set('msk',    msk)
        call cline_cluster2D1%set('ring2',  ring2)
        call cline_cluster2D1%set('lp',     lp1)
        call cline_cluster2D2%set('msk',    msk)
        call cline_cluster2D2%set('lp',     lp2)
        ! execution 1
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A,F6.1)') '>>> STAGE 1, LOW-PASS LIMIT: ',lp1
        write(logfhandle,'(A)') '>>>'
        call cline_cluster2D1%set('projfile', trim(projfile))
        call xcluster2D_distr%execute(cline_cluster2D1)
        last_iter  = nint(cline_cluster2D1%get_rarg('endit'))
        finalcavgs = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter,3)//params%ext
        ! execution 2
        if( cline%defined('maxits') )then
            if( last_iter < params%maxits )then
                write(logfhandle,'(A)') '>>>'
                write(logfhandle,'(A,F6.1)') '>>> STAGE 2, LOW-PASS LIMIT: ',lp2
                write(logfhandle,'(A)') '>>>'
                call cline_cluster2D2%set('projfile', trim(projfile))
                call cline_cluster2D2%set('startit',  real(last_iter+1))
                call cline_cluster2D2%set('refs',     trim(finalcavgs))
                call xcluster2D_distr%execute(cline_cluster2D2)
                last_iter  = nint(cline_cluster2D2%get_rarg('endit'))
                finalcavgs = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter,3)//params%ext
            endif
        endif
        ! restores project file name
        params%projfile = trim(orig_projfile)
        ! update original project
        if( do_scaling )then
            call spproj_sc%read(projfile)
            call spproj%read(orig_projfile_bak)
            call spproj_sc%os_ptcl2D%mul_shifts(1./scale_factor)
            call rescale_cavgs(finalcavgs)
            cavgs = add2fbody(finalcavgs,params%ext,'_even')
            call rescale_cavgs(cavgs)
            cavgs = add2fbody(finalcavgs,params%ext,'_odd')
            call rescale_cavgs(cavgs)
            call spproj%add_cavgs2os_out(trim(finalcavgs), params%smpd, imgkind='cavg')
            spproj%os_ptcl2D = spproj_sc%os_ptcl2D
            spproj%os_cls2D  = spproj_sc%os_cls2D
            ! restores original project and deletes backup & scaled
            call spproj%write(params%projfile)
            call del_file(orig_projfile_bak)
        else
            call spproj%read_segment('out', params%projfile)
            call spproj%add_cavgs2os_out(trim(finalcavgs), params%smpd, imgkind='cavg')
            call spproj%add_frcs2os_out( trim(FRCS_FILE), 'frc2D')
            call spproj%write_segment_inside('out', params%projfile)
        endif
        call spproj_sc%kill
        call spproj%kill
        ! ranking
        finalcavgs_ranked = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter,3)//'_ranked'//params%ext
        call cline_rank_cavgs%set('projfile', trim(params%projfile))
        call cline_rank_cavgs%set('stk',      trim(finalcavgs))
        call cline_rank_cavgs%set('outstk',   trim(finalcavgs_ranked))
        call xrank_cavgs%execute(cline_rank_cavgs)
        ! cleanup
        if( do_scaling ) call simple_rmdir(STKPARTSDIR)
        ! end gracefully
        call simple_end('**** SIMPLE_CLEANUP2D NORMAL STOP ****')
        contains

            subroutine rescale_cavgs(cavgs)
                character(len=*), intent(in) :: cavgs
                type(image)                  :: img, img_pad
                integer                      :: icls, iostat
                call img%new([box,box,1],smpd)
                call img_pad%new([params%box,params%box,1],params%smpd)
                do icls = 1,params%ncls
                    call img%read(cavgs,icls)
                    call img%fft
                    call img%pad(img_pad, backgr=0.)
                    call img_pad%ifft
                    call img_pad%write('tmp_cavgs.mrc',icls)
                enddo
                iostat = simple_rename('tmp_cavgs.mrc',cavgs)
                call img%kill
                call img_pad%kill
            end subroutine

    end subroutine exec_cleanup2D

    subroutine exec_cluster2D_stream( self, cline )
        use simple_projection_frcs, only: projection_frcs
        use simple_oris,            only: oris
        use simple_ori,             only: ori
        class(cluster2D_commander_stream), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        integer,               parameter   :: WAIT_WATCHER        = 30    ! seconds prior to new stack detection
        integer,               parameter   :: ORIGPROJ_WRITEFREQ  = 600   ! 10mins, Frequency at which the original project file should be updated
        integer,               parameter   :: MINBOXSZ            = 72   ! minimum boxsize for scaling
        ! dev settings
        ! integer,               parameter   :: WAIT_WATCHER        = 30    ! seconds prior to new stack detection
        ! integer,               parameter   :: ORIGPROJ_WRITEFREQ  = 60    ! 10mins, Frequency at which the original project file should be updated
        real,                 parameter    :: GREEDY_TARGET_LP   = 15.
        real                               :: SMPD_TARGET         = 4.    ! target sampling distance
        character(len=STDLEN), parameter   :: MICS_SELECTION_FILE = 'stream2D_selection.txt'
        character(len=STDLEN), parameter   :: PROJFILE_BUFFER     = 'buffer.simple'
        character(len=STDLEN), parameter   :: PROJFILE_POOL       = 'pool.simple'
        character(len=STDLEN), parameter   :: PROJFILE2D          = 'cluster2D.simple'
        character(len=STDLEN), parameter   :: SCALE_DIR           = './scaled_stks/'
        logical,               parameter   :: debug_here = .false.
        type(parameters)                   :: params
        type(make_cavgs_commander_distr)   :: xmake_cavgs
        type(rank_cavgs_commander)         :: xrank_cavgs
        type(cmdline)                      :: cline_cluster2D, cline_cluster2D_buffer
        type(cmdline)                      :: cline_make_cavgs, cline_rank_cavgs
        type(sp_project)                   :: orig_proj, stream_proj, buffer_proj, pool_proj
        type(ctfparams)                    :: ctfvars
        type(oris)                         :: os_stk
        type(ori)                          :: o_stk
        character(LONGSTRLEN), allocatable :: spproj_list(:), stk_list(:)
        character(len=:),      allocatable :: spproj_list_fname, orig_projfile, stk
        character(len=STDLEN)              :: str_iter, refs_glob, refs_glob_ranked
        character(len=LONGSTRLEN)          :: buffer_dir
        real    :: orig_smpd, msk, scale_factor, orig_msk, smpd, large_msk, lp_greedy
        integer :: iter, orig_box, box, nptcls_glob, iproj, ncls_glob, n_transfers, large_ring2, iptcl
        integer :: nptcls_glob_prev, n_spprojs, orig_nparts, last_injection, nparts, box_coords(2)
        integer :: origproj_time, max_ncls, nptcls_per_buffer, buffer_ptcls_range(2), pool_iter, i
        logical :: do_autoscale, l_maxed, buffer_exists, do_wait, l_greedy
        if( cline%defined('refine') )then
            if( trim(cline%get_carg('refine')).ne.'greedy' )then
                if( .not.cline%defined('msk') ) THROW_HARD('MSK must be defined!')
            endif
        else
            if( .not.cline%defined('msk') ) THROW_HARD('MSK must be defined!')
        endif
        if( .not. cline%defined('lp')        ) call cline%set('lp',          15.)
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',       20.)
        if( .not. cline%defined('center')    ) call cline%set('center',     'no')
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('lpthresh')  ) call cline%set('lpthresh',    30.)
        if( .not. cline%defined('ndev')      ) call cline%set('ndev',        1.5)
        if( .not. cline%defined('oritype')   ) call cline%set('oritype', 'ptcl2D')
        call cline%set('stream','yes') ! only for parameters determination
        call seed_rnd
        call params%new(cline)
        ! sanity
        if( .not.file_exists(params%projfile) )then
            THROW_HARD('project file: '//trim(params%projfile)//' does not exist!')
        endif
        orig_projfile = trim(params%projfile)
        if( .not.file_exists(params%dir_target) )then
            THROW_HARD('folder: '//trim(params%dir_target)//' does not exist!')
        endif
        call cline%set('stream','no') ! was only for parameters determination
        call cline%set('mkdir','no')
        ! init
        do_autoscale      = params%autoscale.eq.'yes'
        orig_nparts       = params%nparts
        max_ncls          = floor(real(params%ncls)/real(params%ncls_start))*params%ncls_start ! effective maximum # of classes
        nptcls_per_buffer = params%nptcls_per_cls*params%ncls_start         ! # of particles in each buffer
        buffer_exists     = .false.                                         ! whether the buffer exists
        l_maxed           = .false.                                         ! whether all chunks have been merged
        do_wait           = .true.
        spproj_list_fname = filepath(trim(params%dir_target), trim(STREAM_SPPROJFILES))
        ncls_glob         = 0
        n_transfers       = 0
        l_greedy          = trim(params%refine).eq.'greedy'
        lp_greedy         = GREEDY_TARGET_LP
        if( cline%defined('lp') ) lp_greedy = params%lp
        ! for microscopes that don't work too good, automatically turned off after 2 hours
        if(.not.cline%defined('time_inactive'))params%time_inactive = 2*3600
        ! init command-lines
        call cline%delete('lp')
        call cline%delete('refine')
        cline_cluster2D         = cline
        cline_cluster2D_buffer  = cline
        cline_make_cavgs        = cline
        ! buffer classification
        ! down-scaling for fast execution, aggressive stochastic optimisation, no match filter, bi-linear interpolation,
        ! no incremental learning, objective function default is standard cross-correlation (cc), no centering
        call cline_cluster2D_buffer%set('prg',       'cluster2D')
        call cline_cluster2D_buffer%set('projfile',  trim(PROJFILE_BUFFER))
        call cline_cluster2D_buffer%set('projname',  trim(get_fbody(trim(PROJFILE_BUFFER),trim('simple'))))
        call cline_cluster2D_buffer%set('objfun',    'cc')
        call cline_cluster2D_buffer%set('center',    'no')
        call cline_cluster2D_buffer%set('match_filt','no')
        call cline_cluster2D_buffer%set('autoscale', 'no')
        call cline_cluster2D_buffer%set('refine',    'snhc')
        call cline_cluster2D_buffer%set('ptclw',     'no')
        call cline_cluster2D_buffer%delete('update_frac')
        if( l_greedy )then
            call cline_cluster2D_buffer%set('maxits', 10.)
            call cline_cluster2D_buffer%set('refine', 'greedy')
            call cline_cluster2D_buffer%set('lp',     GREEDY_TARGET_LP)
        else
            call cline_cluster2D_buffer%set('maxits',    12.) ! guaranties 3 iterations with withdrawal
        endif
        ! pool classification
        ! down-scaling for fast execution, stochastic optimisation, optional match filter, bi-linear interpolation,
        ! no incremental learning, objective function is standard cross-correlation (cc)
        call cline_cluster2D%set('prg',       'cluster2D')
        call cline_cluster2D%set('autoscale', 'no')
        call cline_cluster2D%set('extr_iter', 100.)
        call cline_cluster2D%set('trs',       MINSHIFT)
        call cline_cluster2D%set('projfile',  trim(PROJFILE_POOL))
        call cline_cluster2D%set('projname',  trim(get_fbody(trim(PROJFILE_POOL),trim('simple'))))
        call cline_cluster2D%set('objfun',    'cc')
        if( .not.cline%defined('match_filt') ) call cline_cluster2D%set('match_filt','no')
        if( l_greedy )then
            call cline_cluster2D%set('center', 'no')
            call cline_cluster2D%set('refine', 'greedy')
            call cline_cluster2D%set('lp',     lp_greedy)
        endif
        call cline_cluster2D%delete('update_frac')
        ! final averaging
        call cline_make_cavgs%set('prg', 'make_cavgs')
        call cline_make_cavgs%delete('autoscale')
        call cline_make_cavgs%delete('remap_cls')
        ! WAIT FOR FIRST STACKS
        nptcls_glob = 0
        do
            if( file_exists(spproj_list_fname) )then
                if( .not.is_file_open(spproj_list_fname) )then
                    call read_mics
                    write(logfhandle,'(A,I8,A,A)')'>>> # OF PARTICLES: ', nptcls_glob, ' : ',cast_time_char(simple_gettime())
                    call flush(6)
                    if( nptcls_glob > nptcls_per_buffer )then
                        exit ! Enough particles to initiate cluster2D
                    endif
                endif
            endif
            call simple_sleep(WAIT_WATCHER)
        enddo
        ! transfer projects info, rename & rewrite
        call orig_proj%read(params%projfile)
        pool_proj%projinfo = orig_proj%projinfo
        pool_proj%compenv  = orig_proj%compenv
        if( orig_proj%jobproc%get_noris()>0 ) pool_proj%jobproc = orig_proj%jobproc
        call pool_proj%projinfo%delete_entry('projname')
        call pool_proj%projinfo%delete_entry('projfile')
        call pool_proj%update_projinfo(cline_cluster2D)
        ! getting general parameters from the first sp_project
        call stream_proj%read(trim(spproj_list(1)))
        orig_box  = stream_proj%get_box()
        orig_smpd = stream_proj%get_smpd()
        orig_msk  = params%msk
        call stream_proj%kill
        params%smpd_targets2D(1) = max(orig_smpd, params%lp*LP2SMPDFAC)
        if( do_autoscale )then
            if( orig_box < MINBOXSZ )then
                do_autoscale = .false.
            else
                call autoscale(orig_box, orig_smpd, SMPD_TARGET, box, smpd, scale_factor)
                if( box < MINBOXSZ )then
                    SMPD_TARGET = orig_smpd * real(orig_box) / real(MINBOXSZ)
                    call autoscale(orig_box, orig_smpd, SMPD_TARGET, box, smpd, scale_factor)
                endif
                if( box == orig_box ) do_autoscale = .false.
            endif
        endif
        if( do_autoscale )then
            msk = orig_msk * scale_factor
        else
            smpd = orig_smpd
            box  = orig_box
            msk  = orig_msk
            scale_factor = 1.
        endif
        large_msk   = real(box/2)-COSMSKHALFWIDTH
        large_ring2 = round2even(0.8*large_msk)
        ! FIRST IMPORT
        allocate(stk_list(n_spprojs))
        call os_stk%new(n_spprojs)
        do iproj=1,n_spprojs
            call stream_proj%read(spproj_list(iproj))
            call stream_proj%os_stk%get_ori(1, o_stk)
            ctfvars         = stream_proj%get_ctfparams('ptcl2D', 1)
            ctfvars%smpd    = smpd
            stk             = stream_proj%get_stkname(1)
            stk_list(iproj) = trim(stk)
            call o_stk%set_ctfvars(ctfvars)
            call os_stk%set_ori(iproj, o_stk)
        enddo
        call stream_proj%kill
        ! updates & scales stacks
        call scale_stks( stk_list ) ! must come first as names updated
        call pool_proj%add_stktab(stk_list, os_stk)
        nptcls_glob = pool_proj%get_nptcls()
        if( file_exists(MICS_SELECTION_FILE) )then
            call flag_selection
        else
            do iproj=1,pool_proj%os_mic%get_noris()
                call pool_proj%os_mic%set(iproj,'state',1.)
            enddo
            do iproj=1,pool_proj%os_stk%get_noris()
                call pool_proj%os_stk%set(iproj,'state',1.)
            enddo
        endif
        ! transfer picking coordinates
        iptcl = 0
        do iproj=1,n_spprojs
            call stream_proj%read(spproj_list(iproj))
            do i = 1,stream_proj%os_ptcl2D%get_noris()
                iptcl = iptcl+1
                if( iptcl > nptcls_glob ) THROW_HARD('picking coordinates indexing error 1')
                if( stream_proj%has_boxcoords(i) )then
                    call stream_proj%get_boxcoords(i,box_coords)
                    call pool_proj%set_boxcoords(iptcl, box_coords)
                endif
            enddo
            call stream_proj%kill
        enddo
        if( iptcl /= nptcls_glob ) THROW_HARD('picking coordinates indexing error 2')
        ! first write
        call pool_proj%write
        ! generates buffer project
        buffer_ptcls_range = 0
        call gen_buffer_from_pool
        ! MAIN LOOP
        last_injection = simple_gettime()
        origproj_time  = last_injection
        pool_iter      = 1
        do iter = 1,9999
            str_iter  = int2str_pad(iter,3)
            pool_iter = min(999,pool_iter)
            write(logfhandle,'(A,I3)')'>>> WAIT CYCLE ',iter
            if( is_timeout(simple_gettime()) )exit
            if( buffer_exists )then
                ! ten iterations of the buffer
                call classify_buffer
                ! particles selection
                call reject_from_buffer
                ! book keeping
                call transfer_buffer_to_pool
                buffer_exists = .false.
                ! cleanup
                call buffer_proj%kill
                call simple_rmdir(buffer_dir)
            endif
            ! one iteration of the whole dataset when not converged,
            ! or after transfer from buffer or each 5 iterations
            if( .not.cline_cluster2D%defined('converged') .or. mod(iter,5)==0 )then
                call classify_pool
                pool_iter = pool_iter+1
                ! rejection each 5 iterations
                if( mod(iter,5)==0 )then
                    call reject_from_pool
                    call pool_proj%write
                endif
                do_wait = .false.
            else if( mod(iter,5) /= 0 )then
                do_wait = .true.
            endif
            ! termination and/or pause
            do while( file_exists(trim(PAUSE_STREAM)) )
                if( file_exists(trim(TERM_STREAM)) ) exit
                call write_singlelineoftext(PAUSE_STREAM, 'PAUSED')
                write(logfhandle,'(A,A)')'>>> CLUSTER2D STREAM PAUSED ',cast_time_char(simple_gettime())
                call simple_sleep(WAIT_WATCHER)
            enddo
            if( file_exists(TERM_STREAM) )then
                write(logfhandle,'(A,A)')'>>> TERMINATING CLUSTER2D STREAM ',cast_time_char(simple_gettime())
                exit
            endif
            if( do_wait ) call simple_sleep(WAIT_WATCHER)
            ! detect new project files
            call append_new_mics
            ! optionally builds new buffer
            call gen_buffer_from_pool
            ! update original project
            if( simple_gettime()-origproj_time > ORIGPROJ_WRITEFREQ )then
                call update_orig_proj
                origproj_time = simple_gettime()
            endif
        enddo
        call qsys_cleanup
        ! updates original project
        call update_orig_proj
        ! class averages at original sampling
        if ( do_autoscale )then
            nptcls_glob = orig_proj%get_nptcls()
            nparts      = calc_nparts(orig_proj, nptcls_glob)
            call orig_proj%kill
            call cline_make_cavgs%set('nparts', real(nparts))
            call cline_make_cavgs%set('ncls',   real(ncls_glob))
            call cline_make_cavgs%set('refs',   refs_glob)
            if( .not.cline%defined('msk') )then
                large_msk = real(orig_box/2)-COSMSKHALFWIDTH
                call cline_make_cavgs%set('msk',large_msk)
            endif
            call xmake_cavgs%execute(cline_make_cavgs)
            call orig_proj%read(orig_projfile)
        endif
        call orig_proj%add_cavgs2os_out(refs_glob, orig_smpd)
        call orig_proj%add_frcs2os_out(trim(FRCS_FILE), 'frc2D')
        call orig_proj%write_segment_inside('out',orig_projfile)
        call orig_proj%kill()
        ! ranking
        refs_glob_ranked = add2fbody(refs_glob,params%ext,'_ranked')
        call cline_rank_cavgs%set('projfile', orig_projfile)
        call cline_rank_cavgs%set('stk',      refs_glob)
        call cline_rank_cavgs%set('outstk',   trim(refs_glob_ranked))
        call xrank_cavgs%execute(cline_rank_cavgs)
        ! cleanup
        call qsys_cleanup
        call simple_rmdir(SCALE_DIR)
        call del_file(PROJFILE_POOL)
        call del_file(PROJFILE_BUFFER)
        call del_file(PROJFILE2D)
        call o_stk%kill
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER2D_STREAM NORMAL STOP ****')
        contains

            subroutine reject_from_buffer
                type(image)          :: img
                logical, allocatable :: cls_mask(:)
                character(STDLEN)    :: refs_buffer
                real                 :: ave, sdev, maxv, minv
                integer              :: nptcls_rejected, ncls_rejected, iptcl
                integer              :: boxmatch, icls, endit, ncls_here, cnt
                if( debug_here ) print *,'in reject from_buffer'; call flush(6)
                ncls_rejected   = 0
                nptcls_rejected = 0
                ncls_here       = buffer_proj%os_cls2D%get_noris()
                boxmatch        = find_boxmatch(box, msk)
                endit           = nint(cline_cluster2D_buffer%get_rarg('endit'))
                refs_buffer     = trim(buffer_dir)//'/'//trim(CAVGS_ITER_FBODY)//int2str_pad(endit,3)//params%ext
                allocate(cls_mask(ncls_here), source=.true.)
                if( debug_here )then
                    ! variance
                    call img%new([box,box,1],smpd)
                    do icls=1,ncls_here
                        call img%read(refs_buffer,icls)
                        call img%stats(ave, sdev, maxv, minv)
                        call buffer_proj%os_cls2D%set(icls,'sdev',sdev)
                    enddo
                endif
                ! resolution and correlation
                call buffer_proj%os_cls2D%find_best_classes(boxmatch,smpd,params%lpthresh,cls_mask,params%ndev)
                if( debug_here ) call buffer_proj%os_cls2D%write('buffer_'//trim(str_iter)//'.txt')
                if( any(cls_mask) )then
                    ncls_rejected = count(.not.cls_mask)
                    do iptcl=1,buffer_proj%os_ptcl2D%get_noris()
                        if( buffer_proj%os_ptcl2D%get_state(iptcl) == 0 )cycle
                        icls = nint(buffer_proj%os_ptcl2D%get(iptcl,'class'))
                        if( cls_mask(icls) ) cycle
                        nptcls_rejected = nptcls_rejected+1
                        call buffer_proj%os_ptcl2D%set(iptcl,'state',0.)
                    enddo
                    do icls=1,ncls_here
                        if( .not.cls_mask(icls) )then
                            call buffer_proj%os_cls2D%set(icls,'pop',0.)
                            call buffer_proj%os_cls2D%set(icls,'corr',-1.)
                        endif
                    enddo
                    call img%new([box,box,1],smpd)
                    cnt = 0
                    do icls=1,ncls_here
                        if( cls_mask(icls) ) cycle
                        cnt = cnt+1
                        call img%read(refs_buffer,icls)
                        call img%write('rejected_'//int2str(pool_iter)//'.mrc',cnt)
                        img = 0.
                        call img%write(refs_buffer,icls)
                    enddo
                    call img%read(refs_buffer, params%ncls_start)
                    call img%write(refs_buffer,params%ncls_start)
                    call img%kill
                    deallocate(cls_mask)
                    write(logfhandle,'(A,I4,A,I6,A)')'>>> REJECTED FROM BUFFER: ',nptcls_rejected,' PARTICLES IN ',ncls_rejected,' CLUSTERS'
                endif
                if( debug_here ) print *,'end reject from_buffer'; call flush(6)
            end subroutine reject_from_buffer

            subroutine reject_from_pool
                type(image)          :: img
                logical, allocatable :: cls_mask(:)
                real                 :: ave, sdev, minv, maxv, ndev_here
                integer              :: nptcls_rejected, ncls_rejected, iptcl
                integer              :: boxmatch, icls, cnt
                if( debug_here ) print *,'in reject from_pool'; call flush(6)
                ncls_rejected   = 0
                nptcls_rejected = 0
                boxmatch        = find_boxmatch(box, msk)
                allocate(cls_mask(ncls_glob), source=.true.)
                if( debug_here )then
                    call img%new([box,box,1],smpd)
                    do icls=1,ncls_glob
                        if( pool_proj%os_cls2D%get_state(icls)==0 .or. pool_proj%os_cls2D%get(icls,'pop')<1. ) cycle
                        call img%read(refs_glob,icls)
                        call img%stats(ave, sdev, maxv, minv)
                        call pool_proj%os_cls2D%set(icls,'sdev',sdev)
                    enddo
                endif
                if( debug_here )call pool_proj%os_cls2D%write('pool_'//int2str(pool_iter)//'.txt')
                ! correlation & resolution
                ndev_here = 1.5*params%ndev ! less stringent rejection
                call pool_proj%os_cls2D%find_best_classes(boxmatch,smpd,params%lpthresh,cls_mask,ndev_here)
                if( debug_here )call pool_proj%os_cls2D%write('pool_aftersel_'//int2str(pool_iter)//'.txt')
                if( .not.all(cls_mask) )then
                    ncls_rejected = count(.not.cls_mask)
                    do iptcl=1,pool_proj%os_ptcl2D%get_noris()
                        if( pool_proj%os_ptcl2D%get_state(iptcl) == 0 )cycle
                        icls = nint(pool_proj%os_ptcl2D%get(iptcl,'class'))
                        if( cls_mask(icls) ) cycle
                        nptcls_rejected = nptcls_rejected+1
                        call pool_proj%os_ptcl2D%set(iptcl,'state',0.)
                    enddo
                    if( nptcls_rejected > 0 )then
                        do icls=1,ncls_glob
                            if( .not.cls_mask(icls) )then
                                call pool_proj%os_cls2D%set(icls,'pop',0.)
                                call pool_proj%os_cls2D%set(icls,'corr',-1.)
                            endif
                        enddo
                        cnt = 0
                        call img%new([box,box,1],smpd)
                        do icls=1,ncls_glob
                            if( cls_mask(icls) ) cycle
                            cnt = cnt+1
                            !if( debug_here )then
                                call img%read(refs_glob,icls)
                                call img%write('rejected_pool_'//int2str(pool_iter)//'.mrc',cnt)
                            !endif
                            img = 0.
                            call img%write(refs_glob,icls)
                        enddo
                        call img%read(refs_glob, ncls_glob)
                        call img%write(refs_glob, ncls_glob)
                        call img%kill
                        deallocate(cls_mask)
                        write(logfhandle,'(A,I4,A,I6,A)')'>>> REJECTED FROM POOL: ',nptcls_rejected,' PARTICLES IN ',ncls_rejected,' CLUSTERS'
                    endif
                else
                    write(logfhandle,'(A,I4,A,I6,A)')'>>> NO PARTICLES FLAGGED FOR REJECTION FROM POOL'
                endif
                if( debug_here ) print *,'end reject from_pool'; call flush(6)
            end subroutine reject_from_pool

            subroutine append_new_mics
                integer :: i, n_spprojs_prev, n_new_spprojs, cnt, iptcl, n_new_ptcls
                if( debug_here ) print *,'in append_new_mics'; call flush(6)
                if( .not.is_file_open(spproj_list_fname) )then
                    do_wait        = .false.
                    n_spprojs_prev = n_spprojs
                    n_spprojs      = nlines(spproj_list_fname)
                    n_new_spprojs  = n_spprojs - n_spprojs_prev
                    if( n_new_spprojs > 0 )then
                        ! fetch new stacks
                        n_new_ptcls = 0
                        if(allocated(spproj_list))deallocate(spproj_list)
                        if(allocated(stk_list))   deallocate(stk_list)
                        allocate(stk_list(n_new_spprojs))
                        call read_filetable(spproj_list_fname, spproj_list)
                        cnt = 0
                        do iproj=n_spprojs_prev+1,n_spprojs
                            cnt = cnt + 1
                            call stream_proj%read(spproj_list(iproj))
                            n_new_ptcls   = n_new_ptcls + stream_proj%get_nptcls()
                            stk           = stream_proj%get_stkname(1)
                            stk_list(cnt) = trim(stk)
                        enddo
                        if( n_new_ptcls > 0 )then
                            call scale_stks( stk_list )
                            ! update project with new images
                            cnt   = 0
                            iptcl = nptcls_glob
                            do iproj=n_spprojs_prev+1,n_spprojs
                                cnt = cnt + 1
                                call stream_proj%read(spproj_list(iproj))
                                ctfvars      = stream_proj%get_ctfparams('ptcl2D', 1)
                                ctfvars%smpd = smpd
                                call pool_proj%add_stk(stk_list(cnt), ctfvars)
                                ! transfer picking coordinates
                                do i = 1,stream_proj%os_ptcl2D%get_noris()
                                    iptcl = iptcl+1
                                    if( stream_proj%has_boxcoords(i) )then
                                        call stream_proj%get_boxcoords(i,box_coords)
                                        call pool_proj%set_boxcoords(iptcl, box_coords)
                                    endif
                                enddo
                            enddo
                            ! global number of particles update
                            nptcls_glob = pool_proj%get_nptcls()
                            ! deactivate new particles by default
                            do iptcl=nptcls_glob-n_new_ptcls+1,nptcls_glob
                                call pool_proj%os_ptcl2D%set(iptcl,'state',0.)
                            enddo
                            write(logfhandle,'(A,I8,A,A)')'>>> # OF PARTICLES: ', nptcls_glob, ' ; ',cast_time_char(simple_gettime())
                            last_injection = simple_gettime()
                            call cline_cluster2D%delete('converged') ! reactivates pool classification
                        endif
                    endif
                    ! updates pool with selection
                    if( file_exists(MICS_SELECTION_FILE) ) call flag_selection
                else
                    do_wait = .true.
                endif
                if( debug_here ) print *,'end append_new_mics'; call flush(6)
            end subroutine append_new_mics

            !>  flags mics and stacks segements with selection
            subroutine flag_selection
                character(len=:), allocatable :: mic_name, mic_name_from_proj
                type(oris)                    :: mics_sel
                integer                       :: nmics, imic, iproj
                logical                       :: included
                nmics = nlines(MICS_SELECTION_FILE)
                call mics_sel%new(nmics)
                call mics_sel%read(MICS_SELECTION_FILE)
                do iproj=1,pool_proj%os_mic%get_noris()
                    if(.not.pool_proj%os_mic%isthere('intg'))cycle
                    call pool_proj%os_mic%getter(iproj,'intg',mic_name_from_proj)
                    ! check whether corresponding mic is selected
                    included = .true.
                    do imic=1,nmics
                        call mics_sel%getter(imic,'intg',mic_name)
                        if( trim(mic_name).eq.trim(mic_name_from_proj) )then
                            included  = mics_sel%get_state(imic) == 1
                            exit
                        endif
                    enddo
                    if( included )then
                        call pool_proj%os_mic%set(iproj,'state',1.)
                        call pool_proj%os_stk%set(iproj,'state',1.)
                    else
                        write(logfhandle,'(A,A)')'>>> DESELECTING MICROGRAPH: ',trim(mic_name)
                        call pool_proj%os_mic%set(iproj,'state',0.)
                        call pool_proj%os_stk%set(iproj,'state',0.)
                    endif
                enddo
                call mics_sel%kill
            end subroutine flag_selection

            integer function calc_nparts(spproj, nptcls_in)
                type(sp_project), intent(in) :: spproj
                integer,          intent(in) :: nptcls_in
                integer :: i, tailsz
                ! adjust number of parts
                tailsz = 0
                do i=nptcls_in,1,-1
                    if( spproj%os_ptcl2D%get_state(i) == 0 )then
                        tailsz = tailsz + 1
                    else
                        exit
                    endif
                enddo
                do calc_nparts=orig_nparts,1,-1
                    if(real(nptcls_in)/real(calc_nparts) > real(tailsz) )exit
                enddo
            end function calc_nparts

            !>  runs iterations of cluster2D for the buffer in a separate folder
            subroutine classify_buffer
                type(cluster2D_commander_distr) :: xcluster2D_distr
                integer :: nptcls, nparts
                if( debug_here ) print *,'in classify_buffer'; call flush(6)
                write(logfhandle,'(A)')'>>> 2D CLASSIFICATION OF NEW BUFFER'
                ! directory structure
                call chdir('buffer2D')
                ! init
                nptcls = buffer_proj%get_nptcls()
                nparts = calc_nparts(buffer_proj, nptcls)
                call buffer_proj%kill
                ! cluster2d execution
                params_glob%projfile = trim('./'//trim(PROJFILE_BUFFER))
                call cline_cluster2D_buffer%set('startit', 1.)
                call cline_cluster2D_buffer%set('box',     real(box))
                call cline_cluster2D_buffer%set('ncls',    real(params%ncls_start))
                call cline_cluster2D_buffer%set('nparts',  real(nparts))
                call cline_cluster2D_buffer%set('msk',     large_msk)
                call cline_cluster2D_buffer%set('ring2',   real(large_ring2))
                call cline_cluster2D_buffer%set('stream',  'yes')
                if( l_greedy )then
                    ! done
                else
                    ! stochastic setting
                    call cline_cluster2D_buffer%set('extr_iter', real(MAX_EXTRLIM2D-2))
                endif
                call cline_cluster2D_buffer%delete('trs')
                call cline_cluster2D_buffer%delete('endit')
                call cline_cluster2D_buffer%delete('converged')
                params_glob%nptcls = nptcls
                call xcluster2D_distr%execute(cline_cluster2D_buffer)
                call getcwd(buffer_dir)
                call rename(PROJFILE_BUFFER, '../'//trim(PROJFILE_BUFFER))
                call chdir('..')
                call buffer_proj%read(PROJFILE_BUFFER)
                params_glob%nparts   = orig_nparts
                params_glob%projfile = trim(orig_projfile)
                if( debug_here ) print *,'end classify_buffer'; call flush(6)
            end subroutine classify_buffer

            !>  runs one iteration of cluster2D for the merged buffers in the cwd
            subroutine classify_pool
                type(cluster2D_commander_distr)  :: xcluster2D_distr
                type(sp_project)                 :: spproj2D
                type(ori)                        :: o_tmp
                character(len=:), allocatable    :: prev_refs
                integer :: iptcl, nptcls, nptcls_sel, istk, nstks, nptcls_glob
                ! transfer from pool
                call cline_cluster2D%set('projfile',  trim(PROJFILE2D))
                call cline_cluster2D%set('projname',  trim(get_fbody(trim(PROJFILE2D),trim('simple'))))
                spproj2D%projinfo = pool_proj%projinfo
                spproj2D%compenv  = pool_proj%compenv
                if( pool_proj%jobproc%get_noris()>0 ) spproj2D%jobproc = pool_proj%jobproc
                call spproj2D%projinfo%delete_entry('projname')
                call spproj2D%projinfo%delete_entry('projfile')
                call spproj2D%update_projinfo(cline_cluster2D)
                nptcls_glob = pool_proj%get_nptcls()
                do iptcl=nptcls_glob,1,-1
                    if( pool_proj%os_ptcl2D%get_state(iptcl) == 1 )exit
                enddo
                nstks  = nint(pool_proj%os_ptcl2D%get(iptcl, 'stkind'))
                nptcls = nint(pool_proj%os_stk%get(nstks, 'top'))
                call spproj2D%os_stk%new(nstks)
                call spproj2D%os_ptcl2D%new(nptcls)
                do iptcl=1,nptcls
                    call pool_proj%os_ptcl2D%get_ori(iptcl, o_tmp)
                    call spproj2D%os_ptcl2D%set_ori(iptcl, o_tmp)
                enddo
                do istk=1,nstks
                    call pool_proj%os_stk%get_ori(istk, o_tmp)
                    call spproj2D%os_stk%set_ori(istk, o_tmp)
                enddo
                spproj2D%os_cls2D = pool_proj%os_cls2D
                spproj2D%os_out   = pool_proj%os_out
                call spproj2D%write(PROJFILE2D)
                nptcls_sel = spproj2D%os_ptcl2D%get_noris(consider_state=.true.)
                write(logfhandle,'(A,I8,A,I4,A)')'>>> 2D CLASSIFICATION OF POOL: ',nptcls_sel,' PARTICLES IN ',ncls_glob,' CLUSTERS'
                ! cluster2d execution
                params_glob%projfile = trim(PROJFILE_POOL)
                call cline_cluster2D%set('startit', real(pool_iter))
                call cline_cluster2D%set('maxits',  real(pool_iter))
                call cline_cluster2D%set('ncls',    real(ncls_glob))
                call cline_cluster2D%set('box',     real(box))
                call cline_cluster2D%set('refs',    trim(refs_glob))
                call cline_cluster2D%set('frcs',    trim(FRCS_FILE))
                call cline_cluster2D%set('msk',     msk)
                call cline_cluster2D%set('stream',  'yes')
                if( l_greedy .and. .not.cline%defined('msk') )then
                    call cline_cluster2D%set('msk',   large_msk)
                    call cline_cluster2D%set('ring2', real(large_ring2))
                endif
                call cline_cluster2D%delete('endit')
                call cline_cluster2D%delete('converged')
                params_glob%nptcls = nptcls
                call xcluster2D_distr%execute(cline_cluster2D)
                params_glob%projfile = trim(orig_projfile)
                refs_glob            = trim(CAVGS_ITER_FBODY)//trim(int2str_pad(pool_iter,3))//trim(params%ext)
                params_glob%nptcls   = nptcls_glob
                ! transfer back to pool
                call spproj2D%read(PROJFILE2D)
                do iptcl=1,nptcls
                    call spproj2D%os_ptcl2D%get_ori(iptcl, o_tmp)
                    call pool_proj%os_ptcl2D%set_ori(iptcl, o_tmp)
                enddo
                do istk=1,nstks
                    call spproj2D%os_stk%get_ori(istk, o_tmp)
                    call pool_proj%os_stk%set_ori(istk, o_tmp)
                enddo
                pool_proj%os_cls2D = spproj2D%os_cls2D
                pool_proj%os_out   = spproj2D%os_out
                ! removes previous references
                if(pool_iter>1)then
                    prev_refs = trim(CAVGS_ITER_FBODY)//trim(int2str_pad(pool_iter-1,3))//trim(params%ext)
                    call del_file(prev_refs)
                    if( l_greedy )then
                        ! no even/odd stored
                    else
                        call del_file(add2fbody(trim(prev_refs),params%ext,'_even'))
                        call del_file(add2fbody(trim(prev_refs),params%ext,'_odd'))
                    endif
                endif
                ! cleanup
                call o_tmp%kill
                call spproj2D%kill
            end subroutine classify_pool

            !>  append the classified buffer to the pool
            subroutine transfer_buffer_to_pool
                type(projection_frcs)         :: frcs_glob, frcs_buffer, frcs_prev
                type(image)                   :: img
                type(ori)                     :: o
                character(len=:), allocatable :: refs_buffer, stkout, stkin
                integer,          allocatable :: cls_pop(:), cls_buffer_pop(:), pinds(:)
                integer                       :: endit, iptcl, ind, state, ncls_here, icls, i, cnt, stat
                real                          :: stkind
                n_transfers = n_transfers+1
                write(logfhandle,'(A,I4)')'>>> TRANSFER BUFFER PARTICLES CLASSIFICATION TO POOL #',n_transfers
                ! max # of classes reached ?
                l_maxed = ncls_glob >= max_ncls
                ! updates # of classes
                if( .not.l_maxed ) ncls_here = ncls_glob+params%ncls_start
                ! transfer class parameters to pool
                if( l_maxed )then
                    cls_pop = nint(pool_proj%os_cls2D%get_all('pop'))
                else
                    if( ncls_glob == 0 )then
                        pool_proj%os_cls2D = buffer_proj%os_cls2D
                    else
                        call pool_proj%os_cls2D%reallocate(ncls_here)
                        do icls=1,params%ncls_start
                            ind = ncls_glob+icls
                            call buffer_proj%os_cls2D%get_ori(icls, o)
                            call o%set('class',real(ind))
                            call pool_proj%os_cls2D%set_ori(ind,o)
                        enddo
                    endif
                endif
                ! transfer particles parameters to pool
                do iptcl=1,buffer_proj%os_ptcl2D%get_noris()
                    call buffer_proj%os_ptcl2D%get_ori(iptcl, o)
                    ! greedy search
                    call o%set('updatecnt', 1.)
                    ! updates class
                    if( l_maxed )then
                        icls = irnd_uni(ncls_glob)
                        do while(cls_pop(icls)==0)
                            icls = irnd_uni(ncls_glob)
                        enddo
                        call o%set('class',real(icls))
                    else
                        icls = nint(o%get('class'))
                        call o%set('class',real(ncls_glob+icls))
                    endif
                    ind    = buffer_ptcls_range(1)+iptcl-1
                    stkind = pool_proj%os_ptcl2D%get(ind,'stkind')
                    ! preserves stack index
                    call o%set('stkind',stkind)
                    ! keep picking coordinates
                    if( pool_proj%has_boxcoords(ind) )then
                        call pool_proj%get_boxcoords(ind, box_coords)
                        call o%set('xpos',real(box_coords(1)))
                        call o%set('ypos',real(box_coords(2)))
                    endif
                    ! updates orientation
                    call pool_proj%os_ptcl2D%set_ori(ind,o)
                enddo
                do iptcl=ind+1,pool_proj%os_ptcl2D%get_noris()
                    call pool_proj%os_ptcl2D%set(iptcl,'state',0.)
                enddo
                ! transfer references to pool
                endit       = nint(cline_cluster2D_buffer%get_rarg('endit'))
                refs_buffer = trim(buffer_dir)//'/'//trim(CAVGS_ITER_FBODY)//int2str_pad(endit,3)//params%ext
                if( .not.l_maxed )then
                    if( ncls_glob == 0 )then
                        ! first time
                        refs_glob = 'start_cavgs'//params%ext
                        call simple_copy_file(refs_buffer,refs_glob)
                        if( l_greedy )then
                            ! no even/odd stored
                        else
                            stkin  = add2fbody(trim(refs_buffer),params%ext,'_even')
                            stkout = add2fbody(trim(refs_glob),params%ext,'_even')
                            call simple_copy_file(stkin,stkout)
                            stkin  = add2fbody(trim(refs_buffer),params%ext,'_odd')
                            stkout = add2fbody(trim(refs_glob),params%ext,'_odd')
                            call simple_copy_file(stkin,stkout)
                        endif
                        call simple_copy_file(trim(buffer_dir)//'/'//trim(FRCS_FILE),FRCS_FILE)
                    else
                        ! class averages &
                        call img%new([box,box,1],smpd)
                        do icls=1,params%ncls_start
                            ind = ncls_glob+icls
                            call img%read(refs_buffer, icls)
                            call img%write(refs_glob, ind)
                        enddo
                        if( l_greedy )then
                            ! no even/odd stored
                        else
                            stkin  = add2fbody(trim(refs_buffer),params%ext,'_even')
                            stkout = add2fbody(trim(refs_glob),params%ext,'_even')
                            do icls=1,params%ncls_start
                                ind = ncls_glob+icls
                                call img%read(stkin, icls)
                                call img%write(stkout, ind)
                            enddo
                            stkin  = add2fbody(trim(refs_buffer),params%ext,'_odd')
                            stkout = add2fbody(trim(refs_glob),params%ext,'_odd')
                            do icls=1,params%ncls_start
                                ind = ncls_glob+icls
                                call img%read(stkin, icls)
                                call img%write(stkout, ind)
                            enddo
                        endif
                        ! FRCs
                        state = 1
                        call frcs_prev%new(ncls_glob, box, smpd, state)
                        call frcs_buffer%new(params%ncls_start, box, smpd, state)
                        call frcs_glob%new(ncls_here, box, smpd, state)
                        call frcs_prev%read(FRCS_FILE)
                        call frcs_buffer%read(trim(buffer_dir)//'/'//trim(FRCS_FILE))
                        do icls=1,ncls_glob
                            call frcs_glob%set_frc(icls,frcs_prev%get_frc(icls, box, state), state)
                        enddo
                        ind = 0
                        do icls=ncls_glob+1,ncls_here
                            ind = ind + 1
                            call frcs_glob%set_frc(icls,frcs_buffer%get_frc(ind, box, state), state)
                        enddo
                        call frcs_glob%write(FRCS_FILE)
                        ! cleanup
                        call frcs_prev%kill
                        call frcs_glob%kill
                        call frcs_buffer%kill
                        call img%kill
                    endif
                    ! global parameters
                    ncls_glob = ncls_here
                endif
                ! remapping
                if( l_maxed .and. trim(params%remap_cls).eq.'yes' )then
                    if( any(cls_pop==0) )then
                        state = 1
                        cls_buffer_pop = nint(buffer_proj%os_cls2D%get_all('pop'))
                        call frcs_buffer%new(params%ncls_start, box, smpd, state)
                        call frcs_glob%new(ncls_glob, box, smpd, state)
                        call frcs_glob%read(FRCS_FILE)
                        call frcs_buffer%read(trim(buffer_dir)//'/'//trim(FRCS_FILE))
                        call img%new([box,box,1],smpd)
                        do icls=1,ncls_glob
                            if( cls_pop(icls)>0 )cycle
                            cnt = 0
                            ind = irnd_uni(params%ncls_start)
                            do while( cls_buffer_pop(ind)==0 )
                                ind = irnd_uni(params%ncls_start)
                                 ! pathological cases and insufficient # of classes
                                cnt = cnt + 1
                                if( cnt>params%ncls_start )exit
                            enddo
                            if( cnt>params%ncls_start )cycle
                            cls_buffer_pop(ind) = 0 ! excludes from being picked again
                            ! classes
                            call img%read(refs_buffer, ind)
                            call img%write(refs_glob, icls)
                            if( l_greedy )then
                                ! no even/odd stored
                            else
                                stkin  = add2fbody(trim(refs_buffer),params%ext,'_even')
                                stkout = add2fbody(trim(refs_glob),params%ext,'_even')
                                call img%read(stkin, ind)
                                call img%write(stkout, icls)
                                stkin  = add2fbody(trim(refs_buffer),params%ext,'_odd')
                                stkout = add2fbody(trim(refs_glob),params%ext,'_odd')
                                call img%read(stkin, ind)
                                call img%write(stkout, icls)
                            endif
                            call buffer_proj%os_cls2D%get_ori(ind, o)
                            call o%set('class',real(icls))
                            call pool_proj%os_cls2D%set_ori(icls,o)
                            ! frcs
                            call frcs_glob%set_frc(icls,frcs_buffer%get_frc(ind, box, state), state)
                            ! assignments
                            call buffer_proj%os_ptcl2D%get_pinds(ind,'class',pinds,consider_w=.false.)
                            pinds = pinds + buffer_ptcls_range(1)
                            do i=1,size(pinds)
                                iptcl = pinds(i)
                                call pool_proj%os_ptcl2D%set(iptcl,'class',real(icls))
                            enddo
                        enddo
                        if( l_greedy )then
                            ! no even/odd stored
                        else
                            call img%read(refs_glob, ncls_glob)
                            call img%write(refs_glob, ncls_glob)
                            stkout = add2fbody(trim(refs_glob),params%ext,'_even')
                            call img%read(stkout, ncls_glob)
                            call img%write(stkout, ncls_glob)
                            stkout = add2fbody(trim(refs_glob),params%ext,'_odd')
                            call img%read(stkout, ncls_glob)
                            call img%write(stkout, ncls_glob)
                        endif
                        call frcs_glob%write(FRCS_FILE)
                        call frcs_glob%kill
                        call frcs_buffer%kill
                        call img%kill
                    endif
                endif
                ! preserve buffer
                stat = rename(refs_buffer,'cavgs_buffer'//int2str_pad(n_transfers,4)//params%ext)
                ! cleanup
                call o%kill
                if(allocated(cls_pop))deallocate(cls_pop)
                if(allocated(cls_buffer_pop))deallocate(cls_buffer_pop)
                if( debug_here )print *,'end transfer_buffer_to_pool'; call flush(6)
            end subroutine transfer_buffer_to_pool

            subroutine read_mics
                character(len=:), allocatable :: mic_name, mic_name_from_proj
                type(oris)                    :: mics_sel
                integer                       :: nptcls, nmics, imic
                logical                       :: included, do_selection
                if( debug_here )print *,'in read_mics'; call flush(6)
                do_selection = file_exists(MICS_SELECTION_FILE)
                call read_filetable(spproj_list_fname, spproj_list)
                if( do_selection )then
                    nmics = nlines(MICS_SELECTION_FILE)
                    call mics_sel%new(nmics)
                    call mics_sel%read(MICS_SELECTION_FILE)
                endif
                ! determine number of particles
                nptcls_glob_prev = nptcls_glob
                nptcls_glob      = 0
                if( allocated(spproj_list) )then
                    n_spprojs = size(spproj_list)
                    do iproj = 1,n_spprojs
                        call stream_proj%read(spproj_list(iproj))
                        nptcls   = stream_proj%get_nptcls()
                        included = .true. ! included by default
                        if( do_selection )then
                            call stream_proj%os_mic%getter(1,'intg',mic_name_from_proj)
                            ! check whether corresponding mic is selected
                            do imic=1,nmics
                                call mics_sel%getter(imic,'intg',mic_name)
                                if( trim(mic_name).eq.trim(mic_name_from_proj) )then
                                    included  = mics_sel%get_state(imic) == 1
                                    exit
                                endif
                            enddo
                        endif
                        if( included ) nptcls_glob = nptcls_glob + nptcls
                    enddo
                    call stream_proj%kill
                else
                    n_spprojs = 0
                endif
                if(do_selection)call mics_sel%kill
                if( debug_here )print *,'end read_mics'; call flush(6)
            end subroutine read_mics

            !> builds and writes new buffer project from the pool
            subroutine gen_buffer_from_pool
                type(ctfparams) :: ctfparms
                character(len=:), allocatable :: stkname
                integer :: iptcl, stkind, ind_in_stk, imic, micind, nptcls_here
                integer :: fromp, top, istk, state, nmics_here, nptcls_tot
                if( debug_here )print *,'in gen_buffer_from_pool'; call flush(6)
                buffer_exists = .false.
                call buffer_proj%kill
                ! to account for directory structure
                call simple_rmdir('buffer2D') ! clean previous round
                call simple_mkdir('buffer2D')
                call simple_chdir('buffer2D')
                ! determines whether there are enough particles for a buffer
                if( buffer_ptcls_range(2) > 0 )then
                    call pool_proj%map_ptcl_ind2stk_ind('ptcl2D', buffer_ptcls_range(2), stkind, ind_in_stk)
                else
                    stkind = 0
                endif
                micind = stkind + 1
                nptcls_here = 0
                do istk=micind,pool_proj%os_stk%get_noris()
                    if(.not.pool_proj%os_stk%isthere(istk,'state'))then
                        write(logfhandle,*)'error: missing state flag; gen_buffer_from_pool'
                        stop
                    endif
                    state = pool_proj%os_stk%get_state(istk)
                    if( state == 0 )cycle
                    fromp       = nint(pool_proj%os_stk%get(istk,'fromp'))
                    top         = nint(pool_proj%os_stk%get(istk,'top'))
                    nptcls_here = nptcls_here + top-fromp+1
                enddo
                if( nptcls_here < nptcls_per_buffer )then
                    call simple_chdir('..')
                    return
                endif
                ! builds buffer if enough particles
                buffer_proj%projinfo = pool_proj%projinfo
                buffer_proj%compenv  = pool_proj%compenv
                if( pool_proj%jobproc%get_noris()>0 ) buffer_proj%jobproc = pool_proj%jobproc
                call buffer_proj%projinfo%delete_entry('projname')
                call buffer_proj%projinfo%delete_entry('projfile')
                call buffer_proj%update_projinfo(cline_cluster2D_buffer)
                istk        = 0
                nptcls_here = 0
                nmics_here  = 0
                nptcls_tot  = 0
                do imic=micind,pool_proj%os_stk%get_noris()
                    state = pool_proj%os_stk%get_state(imic)
                    istk  = istk+1
                    call pool_proj%os_stk%getter(imic,'stk',stkname)
                    ! to account for directory structure
                    if( stkname(1:1) /= '/') stkname = '../'//trim(stkname)
                    ! import
                    ctfparms = pool_proj%os_stk%get_ctfvars(imic)
                    call buffer_proj%add_stk(stkname, ctfparms)
                    fromp      = nint(buffer_proj%os_stk%get(istk,'fromp'))
                    top        = nint(buffer_proj%os_stk%get(istk,'top'))
                    call buffer_proj%os_stk%set(istk,'state',real(state))
                    do iptcl=fromp,top
                        call buffer_proj%os_ptcl2D%set(iptcl,'state',real(state))
                    enddo
                    if( state == 0 ) cycle
                    nmics_here  = nmics_here+1
                    nptcls_here = nptcls_here + top-fromp+1
                    if( nptcls_here > nptcls_per_buffer )exit
                enddo
                buffer_ptcls_range(1) = buffer_ptcls_range(2)+1
                buffer_ptcls_range(2) = buffer_ptcls_range(1)+top-1
                call buffer_proj%write
                call simple_chdir('..')
                buffer_exists = .true.
                write(logfhandle,'(A,I4,A,I6,A)')'>>> BUILT NEW BUFFER WITH ', nmics_here, ' MICROGRAPHS, ',nptcls_here, ' PARTICLES'
            end subroutine gen_buffer_from_pool

            logical function is_timeout( time_now )
                integer, intent(in) :: time_now
                is_timeout = .false.
                if(time_now-last_injection > params%time_inactive)then
                    write(logfhandle,'(A,A)')'>>> TIME LIMIT WITHOUT NEW IMAGES REACHED: ',cast_time_char(time_now)
                    is_timeout = .true.
                    if( .not.cline_cluster2D%defined('converged') )is_timeout = .false.
                else if(time_now-last_injection > 3600)then
                    write(logfhandle,'(A,A)')'>>> OVER ONE HOUR WITHOUT NEW PARTICLES: ',cast_time_char(time_now)
                    call flush(6)
                endif
                return
            end function is_timeout

            subroutine update_orig_proj
                type(ori)                          :: o
                type(oris)                         :: o_stks
                type(ctfparams)                    :: ctfparms
                character(LONGSTRLEN), allocatable :: stk_list(:), mic_list(:)
                real                               :: rstate
                integer                            :: nprev, n2append, cnt
                logical                            :: has_mics
                write(logfhandle,'(A,A)')'>>> UPDATING PROJECT AT: ',cast_time_char(simple_gettime())
                call orig_proj%read(orig_projfile)
                nprev    = orig_proj%os_stk%get_noris()
                n2append = n_spprojs-nprev
                has_mics = .false.
                if( n2append > 0 )then
                    allocate(stk_list(n2append),mic_list(n2append))
                    call o_stks%new(n2append)
                    cnt   = 0
                    do iproj=nprev+1,n_spprojs
                        cnt = cnt + 1
                        ! stk
                        call stream_proj%read(spproj_list(iproj))
                        call stream_proj%os_stk%get_ori(1, o)
                        stk_list(cnt) = trim(stream_proj%get_stkname(1))
                        ctfparms      = stream_proj%get_ctfparams('ptcl2D', 1)
                        call o%set_ctfvars(ctfparms)
                        call o_stks%set_ori(cnt,o)
                        ! mic
                        if( stream_proj%get_nintgs() == 1 )then
                            has_mics      = .true.
                            ctfparms      = stream_proj%get_micparams(1)
                            mic_list(cnt) = trim(stream_proj%os_mic%get_static(1,'intg'))
                        endif
                    enddo
                    call orig_proj%add_stktab(stk_list, o_stks)
                    if( has_mics )then
                        call orig_proj%add_movies(mic_list, ctfparms)
                    endif
                    ! updates 2D, os_out not updated as not correct scale
                    do iproj=nprev+1,n_spprojs
                        rstate = pool_proj%os_stk%get(iproj,'state')
                        call orig_proj%os_stk%set(iproj,'state', rstate)
                        if( has_mics )call orig_proj%os_mic%set(iproj,'state', rstate)
                    enddo
                    orig_proj%os_cls2D  = pool_proj%os_cls2D
                    orig_proj%os_ptcl2D = pool_proj%os_ptcl2D
                    if( do_autoscale )call orig_proj%os_ptcl2D%mul_shifts( 1./scale_factor )
                    ! write
                    call orig_proj%write
                endif
                ! cleanup
                call o%kill
                call o_stks%kill
                if(allocated(mic_list))     deallocate(mic_list)
                if(allocated(stk_list))     deallocate(stk_list)
                if( debug_here )print *,'end update_orig_proj'; call flush(6)
            end subroutine update_orig_proj

            subroutine scale_stks( stk_fnames )
                character(len=*), allocatable, intent(inout) :: stk_fnames(:)
                character(len=*), parameter   :: SCALE_FILETAB = 'stkscale.txt'
                character(len=:), allocatable :: fname
                type(qsys_env) :: qenv
                type(cmdline)  :: cline_scale
                integer        :: istk
                if( .not.do_autoscale )return
                if( .not.allocated(stk_fnames) )return
                call simple_mkdir(SCALE_DIR, errmsg= "commander_stream_wflows:: cluster2D_stream scale_stks")
                call qenv%new(params%nparts)
                call cline_scale%set('prg',        'scale')
                call cline_scale%set('smpd',       orig_smpd)
                call cline_scale%set('box',        real(orig_box))
                call cline_scale%set('newbox',     real(box))
                call cline_scale%set('filetab',    trim(SCALE_FILETAB))
                call cline_scale%set('nthr',       real(params%nthr))
                call cline_scale%set('dir_target', trim(SCALE_DIR))
                call cline_scale%set('stream',     'yes')
                call write_filetable(trim(SCALE_FILETAB), stk_fnames)
                call qenv%exec_simple_prg_in_queue(cline_scale, 'JOB_FINISHED_1')
                call qsys_cleanup
                do istk=1,size(stk_fnames)
                    fname            = add2fbody(stk_fnames(istk), params%ext, SCALE_SUFFIX)
                    stk_fnames(istk) = filepath(trim(SCALE_DIR), basename(fname))
                enddo
                call del_file(SCALE_FILETAB)
                call qenv%kill
            end subroutine scale_stks

    end subroutine exec_cluster2D_stream

    subroutine exec_cluster2D_autoscale( self, cline )
        use simple_commander_project, only: scale_project_commander_distr, prune_project_commander_distr
        use simple_commander_imgproc, only: scale_commander
        class(cluster2D_autoscale_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        ! constants
        integer,               parameter :: MAXITS_STAGE1      = 10
        integer,               parameter :: MAXITS_STAGE1_EXTR = 15
        character(len=STDLEN), parameter :: orig_projfile_bak  = 'orig_bak.simple'
        ! commanders
        type(prune_project_commander_distr) :: xprune_project
        type(make_cavgs_commander_distr)    :: xmake_cavgs
        type(cluster2D_commander_distr)     :: xcluster2D_distr
        type(rank_cavgs_commander)          :: xrank_cavgs
        type(scale_commander)               :: xscale
        type(scale_project_commander_distr) :: xscale_distr
        ! command lines
        type(cmdline) :: cline_cluster2D_stage1
        type(cmdline) :: cline_cluster2D_stage2
        type(cmdline) :: cline_scalerefs, cline_scale1, cline_scale2
        type(cmdline) :: cline_make_cavgs, cline_rank_cavgs, cline_prune_project
        ! other variables
        type(parameters)              :: params
        type(sp_project)              :: spproj, spproj_sc
        character(len=:), allocatable :: projfile_sc, orig_projfile
        character(len=LONGSTRLEN)     :: finalcavgs, finalcavgs_ranked, refs_sc
        real     :: scale_stage1, scale_stage2, trs_stage2
        integer  :: nparts, last_iter_stage1, last_iter_stage2, status
        integer  :: nptcls_sel
        logical  :: scaling
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        if( .not. cline%defined('lpstart')   ) call cline%set('lpstart',    15. )
        if( .not. cline%defined('lpstop')    ) call cline%set('lpstop',      8. )
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      30. )
        if( .not. cline%defined('maxits')    ) call cline%set('maxits',     30. )
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        call cline%delete('clip')
        ! master parameters
        call params%new(cline)
        nparts = params%nparts
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! read project file
        call spproj%read(params%projfile)
        orig_projfile = trim(params%projfile)
        ! sanity checks
        if( spproj%get_nptcls() == 0 )then
            THROW_HARD('No particles found in project file: '//trim(params%projfile)//'; exec_cluster2D_autoscale')
        endif
        ! delete any previous solution
        if( .not. spproj%is_virgin_field(params%oritype) )then
            ! removes previous cluster2D solution (states are preserved)
            call spproj%os_ptcl2D%delete_2Dclustering
            call spproj%write_segment_inside(params%oritype)
        endif
        ! automated pruning
        nptcls_sel = spproj%os_ptcl2D%get_noris(consider_state=.true.)
        if( nptcls_sel < nint(PRUNE_FRAC*real(spproj%os_ptcl2D%get_noris())) )then
            write(logfhandle,'(A)')'>>> AUTO-PRUNING PROJECT FILE'
            call spproj%kill
            cline_prune_project = cline
            call xprune_project%execute(cline_prune_project)
            call spproj%read(params%projfile)
            params%nptcls = nptcls_sel
        endif
        ! refinement flag
        if(.not.cline%defined('refine')) call cline%set('refine','snhc')
        ! splitting
        call spproj%split_stk(params%nparts, dir=PATH_PARENT)
        ! general options planning
        if( params%l_autoscale )then
            ! this workflow executes two stages of CLUSTER2D
            ! Stage 1: high down-scaling for fast execution, hybrid extremal/SHC optimisation for
            !          improved population distribution of clusters, no incremental learning,
            !          objective function is standard cross-correlation (cc)
            cline_cluster2D_stage1 = cline
            call cline_cluster2D_stage1%set('objfun',     'cc')
            call cline_cluster2D_stage1%set('match_filt', 'no')
            if( params%l_frac_update )then
                call cline_cluster2D_stage1%delete('update_frac') ! no incremental learning in stage 1
                call cline_cluster2D_stage1%set('maxits', real(MAXITS_STAGE1_EXTR))
            else
                call cline_cluster2D_stage1%set('maxits', real(MAXITS_STAGE1))
            endif
            ! Scaling
            call spproj%scale_projfile(params%smpd_targets2D(1), projfile_sc,&
                &cline_cluster2D_stage1, cline_scale1, dir=trim(STKPARTSDIR))
            call spproj%kill
            scale_stage1 = cline_scale1%get_rarg('scale')
            scaling      = basename(projfile_sc) /= basename(orig_projfile)
            if( scaling )then
                call cline_scale1%delete('smpd') !!
                call cline_scale1%set('state',1.)
                call simple_mkdir(trim(STKPARTSDIR),errmsg="commander_hlev_wflows :: exec_cluster2D_autoscale;  ")
                call xscale_distr%execute( cline_scale1 )
                ! rename scaled projfile and stash original project file
                ! such that the scaled project file has the same name as the original and can be followed from the GUI
                call simple_copy_file(orig_projfile, orig_projfile_bak)
                call spproj%read_non_data_segments(projfile_sc)
                call spproj%projinfo%set(1,'projname',get_fbody(orig_projfile,METADATA_EXT,separator=.false.))
                call spproj%projinfo%set(1,'projfile',orig_projfile)
                call spproj%write_non_data_segments(projfile_sc)
                call spproj%kill
                status = simple_rename(projfile_sc,orig_projfile)
                deallocate(projfile_sc)
                ! scale references
                if( cline%defined('refs') )then
                    call cline_scalerefs%set('stk', trim(params%refs))
                    refs_sc = 'refs'//trim(SCALE_SUFFIX)//params%ext
                    call cline_scalerefs%set('outstk', trim(refs_sc))
                    call cline_scalerefs%set('smpd', params%smpd)
                    call cline_scalerefs%set('newbox', cline_scale1%get_rarg('newbox'))
                    call xscale%execute(cline_scalerefs)
                    call cline_cluster2D_stage1%set('refs',trim(refs_sc))
                endif
            endif
            ! execution
            call cline_cluster2D_stage1%set('projfile', trim(orig_projfile))
            call xcluster2D_distr%execute(cline_cluster2D_stage1)
            last_iter_stage1 = nint(cline_cluster2D_stage1%get_rarg('endit'))
            ! update original project backup and copy to original project file
            if( scaling )then
                call spproj_sc%read_segment('ptcl2D', orig_projfile)
                call spproj_sc%os_ptcl2D%mul_shifts(1./scale_stage1)
                call spproj%read(orig_projfile_bak)
                spproj%os_ptcl2D = spproj_sc%os_ptcl2D
                call spproj%write_segment_inside('ptcl2D',fname=orig_projfile_bak)
                call spproj%kill()
                call simple_copy_file(orig_projfile_bak, orig_projfile)
                ! clean stacks
                call simple_rmdir(STKPARTSDIR)
            endif
            ! Stage 2: refinement stage, less down-scaling, no extremal updates, incremental
            !          learning for acceleration
            cline_cluster2D_stage2 = cline
            call cline_cluster2D_stage2%delete('refs')
            call cline_cluster2D_stage2%set('startit', real(last_iter_stage1 + 1))
            if( cline%defined('update_frac') )then
                call cline_cluster2D_stage2%set('update_frac', params%update_frac)
            endif
            ! Scaling
            call spproj%read(orig_projfile)
            call spproj%scale_projfile( params%smpd_targets2D(2), projfile_sc,&
                &cline_cluster2D_stage2, cline_scale2, dir=trim(STKPARTSDIR))
            call spproj%kill
            scale_stage2 = cline_scale2%get_rarg('scale')
            scaling      = basename(projfile_sc) /= basename(orig_projfile)
            if( scaling )then
                call cline_scale2%delete('smpd') !!
                call cline_scale2%set('state',1.)
                call xscale_distr%execute( cline_scale2 )
                ! rename scaled projfile and stash original project file
                ! such that the scaled project file has the same name as the original and can be followed from the GUI
                call spproj%read_non_data_segments(projfile_sc)
                call spproj%projinfo%set(1,'projname',get_fbody(orig_projfile,METADATA_EXT,separator=.false.))
                call spproj%projinfo%set(1,'projfile',orig_projfile)
                call spproj%write_non_data_segments(projfile_sc)
                call spproj%kill
                status = simple_rename(projfile_sc,orig_projfile)
                deallocate(projfile_sc)
            endif
            trs_stage2 = MSK_FRAC*cline_cluster2D_stage2%get_rarg('msk')
            trs_stage2 = min(MAXSHIFT,max(MINSHIFT,trs_stage2))
            call cline_cluster2D_stage2%set('trs', trs_stage2)
            ! execution
            call cline_cluster2D_stage2%set('projfile', trim(orig_projfile))
            call xcluster2D_distr%execute(cline_cluster2D_stage2)
            last_iter_stage2 = nint(cline_cluster2D_stage2%get_rarg('endit'))
            finalcavgs       = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//params%ext
            ! Updates project and references
            if( scaling )then
                ! shift modulation
                call spproj_sc%read_segment('ptcl2D', orig_projfile)
                call spproj_sc%os_ptcl2D%mul_shifts(1./scale_stage2)
                call spproj%read(orig_projfile_bak)
                spproj%os_ptcl2D = spproj_sc%os_ptcl2D
                call spproj%write_segment_inside('ptcl2D',fname=orig_projfile_bak)
                call spproj%kill()
                call spproj_sc%kill()
                status = simple_rename(orig_projfile_bak,orig_projfile)
                ! clean stacks
                call simple_rmdir(STKPARTSDIR)
                ! original scale references
                cline_make_cavgs = cline ! ncls is transferred here
                call cline_make_cavgs%delete('autoscale')
                call cline_make_cavgs%delete('balance')
                call cline_make_cavgs%set('prg',      'make_cavgs')
                call cline_make_cavgs%set('projfile', orig_projfile)
                call cline_make_cavgs%set('nparts',   real(nparts))
                call cline_make_cavgs%set('refs',     trim(finalcavgs))
                call xmake_cavgs%execute(cline_make_cavgs)
            endif
        else
            ! no auto-scaling
            call xcluster2D_distr%execute(cline)
            last_iter_stage2 = nint(cline%get_rarg('endit'))
            finalcavgs       = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//params%ext
        endif
        ! adding cavgs & FRCs to project
        params%projfile = trim(orig_projfile)
        call spproj%read( params%projfile )
        call spproj%add_frcs2os_out( trim(FRCS_FILE), 'frc2D')
        call spproj%add_cavgs2os_out(trim(finalcavgs), spproj%get_smpd(), imgkind='cavg')
        call spproj%write_segment_inside('out')
        call spproj%kill()
        ! ranking
        finalcavgs_ranked = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//'_ranked'//params%ext
        call cline_rank_cavgs%set('projfile', trim(params%projfile))
        call cline_rank_cavgs%set('stk',      trim(finalcavgs))
        call cline_rank_cavgs%set('outstk',   trim(finalcavgs_ranked))
        call xrank_cavgs%execute( cline_rank_cavgs )
        ! cleanup
        call del_file('start2Drefs'//params%ext)
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER2D NORMAL STOP ****')
    end subroutine exec_cluster2D_autoscale

    subroutine exec_cluster2D_distr( self, cline )
        use simple_procimgfile
        class(cluster2D_commander_distr), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        ! commanders
        type(check_2Dconv_commander)     :: xcheck_2Dconv
        type(make_cavgs_commander_distr) :: xmake_cavgs
        ! command lines
        type(cmdline) :: cline_check_2Dconv
        type(cmdline) :: cline_cavgassemble
        type(cmdline) :: cline_make_cavgs
        ! other variables
        type(parameters)          :: params
        type(builder)             :: build
        type(qsys_env)            :: qenv
        character(len=LONGSTRLEN) :: refs, refs_even, refs_odd, str, str_iter
        integer                   :: iter
        type(chash)               :: job_descr
        real                      :: frac_srch_space
        logical                   :: l_stream
        if( .not. cline%defined('lpstart')   ) call cline%set('lpstart',    15. )
        if( .not. cline%defined('lpstop')    ) call cline%set('lpstop',      8. )
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      30. )
        if( .not. cline%defined('maxits')    ) call cline%set('maxits',     30. )
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('oritype')   ) call cline%set('oritype', 'ptcl2D')
        ! streaming gets its own logics because it is an exception to rules in parameters
        l_stream = .false.
        if( cline%defined('stream') ) l_stream = trim(cline%get_carg('stream')).eq.'yes'
        call cline%set('stream','no')
        ! builder & params
        call build%init_params_and_build_spproj(cline, params)
        ! sanity check
        if( build%spproj%get_nptcls() == 0 )then
            THROW_HARD('no particles found! exec_cluster2D_distr')
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        if( l_stream )then
            call job_descr%set('stream','yes')
            ! need to be explicity defined for streaming
            call job_descr%set('box',   int2str(params%box))
            call job_descr%set('smpd',  real2str(params%smpd))
            call job_descr%set('nptcls',int2str(params%nptcls))
        endif
        ! splitting
        call build%spproj%split_stk(params%nparts)
        ! prepare command lines from prototype master
        cline_check_2Dconv = cline
        cline_cavgassemble = cline
        cline_make_cavgs   = cline ! ncls is transferred here
        ! initialise static command line parameters and static job description parameters
        call cline_cavgassemble%set('prg', 'cavgassemble')
        call cline_cavgassemble%set('nthr', 0.) ! to ensure use of all resources in assembly
        call cline_make_cavgs%set('prg',   'make_cavgs')
        if( l_stream )call cline_check_2Dconv%set('stream','yes')
        ! execute initialiser
        if( .not. cline%defined('refs') )then
            refs             = 'start2Drefs' // params%ext
            params%refs      = trim(refs)
            params%refs_even = 'start2Drefs_even'//params%ext
            params%refs_odd  = 'start2Drefs_odd'//params%ext
            if( build%spproj%is_virgin_field('ptcl2D') .or. params%startit == 1 )then
                if( params%tseries .eq. 'yes' )then
                    call selection_from_tseries_imgfile(build%spproj, params%refs, params%box, params%ncls)
                else
                    call random_selection_from_imgfile(build%spproj, params%refs, params%box, params%ncls)
                endif
                call copy_imgfile(trim(params%refs), trim(params%refs_even), params%smpd, [1,params%ncls])
                call copy_imgfile(trim(params%refs), trim(params%refs_odd),  params%smpd, [1,params%ncls])
            else
                call cline_make_cavgs%set('refs', params%refs)
                call xmake_cavgs%execute(cline_make_cavgs)
            endif
        else
            refs = trim(params%refs)
        endif
        ! variable neighbourhood size
        if( cline%defined('extr_iter') )then
            params%extr_iter = params%extr_iter - 1
        else
            params%extr_iter = params%startit - 1
        endif
        ! deal with eo partitioning
        if( build%spproj_field%get_nevenodd() == 0 )then
            if( params%tseries .eq. 'yes' )then
                call build%spproj_field%partition_eo(tseries=.true.)
            else
                call build%spproj_field%partition_eo
            endif
            call build%spproj%write_segment_inside(params%oritype)
        endif
        ! main loop
        iter = params%startit - 1
        do
            iter = iter + 1
            str_iter = int2str_pad(iter,3)
            write(logfhandle,'(A)')   '>>>'
            write(logfhandle,'(A,I6)')'>>> ITERATION ', iter
            write(logfhandle,'(A)')   '>>>'
            ! cooling of the randomization rate
            params%extr_iter = params%extr_iter + 1
            call job_descr%set('extr_iter', trim(int2str(params%extr_iter)))
            call cline%set('extr_iter', real(params%extr_iter))
            ! updates
            call job_descr%set('refs', trim(refs))
            call job_descr%set('startit', int2str(iter))
            ! the only FRC we have is from the previous iteration, hence the iter - 1
            call job_descr%set('frcs', trim(FRCS_FILE))
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=trim(ALGN_FBODY))
            ! merge orientation documents
            call build%spproj%merge_algndocs(params%nptcls, params%nparts, 'ptcl2D', ALGN_FBODY)
            ! assemble class averages
            refs      = trim(CAVGS_ITER_FBODY) // trim(str_iter)            // params%ext
            refs_even = trim(CAVGS_ITER_FBODY) // trim(str_iter) // '_even' // params%ext
            refs_odd  = trim(CAVGS_ITER_FBODY) // trim(str_iter) // '_odd'  // params%ext
            call cline_cavgassemble%set('refs', trim(refs))
            call qenv%exec_simple_prg_in_queue(cline_cavgassemble, 'CAVGASSEMBLE_FINISHED')
            ! check convergence
            call xcheck_2Dconv%execute(cline_check_2Dconv)
            frac_srch_space = 0.
            if( iter > 1 ) frac_srch_space = cline_check_2Dconv%get_rarg('frac_srch')
            ! the below activates shifting & automasking
            if( iter > 3 .and. (frac_srch_space >= FRAC_SH_LIM .or. cline_check_2Dconv%defined('trs')) )then
                if( .not.job_descr%isthere('trs') )then
                    ! activates shift search
                    str = real2str(cline_check_2Dconv%get_rarg('trs'))
                    call job_descr%set('trs', trim(str) )
                endif
            endif
            if( cline_check_2Dconv%get_carg('converged').eq.'yes' .or. iter==params%maxits )then
                if( cline_check_2Dconv%get_carg('converged').eq.'yes' )call cline%set('converged','yes')
                exit
            endif
        end do
        call qsys_cleanup
        ! report the last iteration on exit
        call cline%delete( 'startit' )
        call cline%set('endit', real(iter))
        ! end gracefully
        call build%spproj_field%kill
        call simple_end('**** SIMPLE_DISTR_CLUSTER2D NORMAL STOP ****')
    end subroutine exec_cluster2D_distr

    subroutine exec_cluster2D( self, cline )
        use simple_strategy2D_matcher, only: cluster2D_exec
        class(cluster2D_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer :: startit, ncls_from_refs, lfoo(3)
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_strategy2D_tbox(cline, params)
        if( cline%defined('refs') )then
            call find_ldim_nptcls(params%refs, lfoo, ncls_from_refs)
            ! consistency check
            if( params%ncls /=  ncls_from_refs ) THROW_HARD('nrefs /= inputted ncls')
        endif
        startit = 1
        if( cline%defined('startit') )startit = params%startit
        if( startit == 1 )call build%spproj_field%clean_updatecnt
        ! execute
        if( .not. cline%defined('outfile') ) THROW_HARD('need unique output file for parallel jobs')
        call cluster2D_exec( cline, startit ) ! partition or not, depending on 'part'
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER2D NORMAL STOP ****')
        call qsys_job_finished('simple_commander_cluster2D :: exec_cluster2D')
    end subroutine exec_cluster2D

    subroutine exec_cavgassemble( self, cline )
        use simple_classaverager
        class(cavgassemble_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_strategy2D_tbox(cline, params)
        call cavger_new( 'class')
        call cavger_transf_oridat( build%spproj )
        call cavger_assemble_sums_from_parts()
        if( cline%defined('which_iter') )then
            params%refs      = trim(CAVGS_ITER_FBODY)//int2str_pad(params%which_iter,3)//params%ext
            params%refs_even = trim(CAVGS_ITER_FBODY)//int2str_pad(params%which_iter,3)//'_even'//params%ext
            params%refs_odd  = trim(CAVGS_ITER_FBODY)//int2str_pad(params%which_iter,3)//'_odd'//params%ext
        else if( .not. cline%defined('refs') )then
            params%refs      = 'start2Drefs'//params%ext
            params%refs_even = 'start2Drefs_even'//params%ext
            params%refs_odd  = 'start2Drefs_odd'//params%ext
        endif
        call cavger_calc_and_write_frcs_and_eoavg(params%frcs)
        ! classdoc gen needs to be after calc of FRCs
        call cavger_gen2Dclassdoc(build%spproj)
        ! write references
        call cavger_write(trim(params%refs),      'merged')
        call cavger_write(trim(params%refs_even), 'even'  )
        call cavger_write(trim(params%refs_odd),  'odd'   )
        call cavger_kill()
        ! write project
        call build%spproj%write_segment_inside('cls2D', params%projfile)
        ! end gracefully
        call simple_end('**** SIMPLE_CAVGASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        call simple_touch('CAVGASSEMBLE_FINISHED', errmsg='In: commander_rec :: eo_cavgassemble ')
    end subroutine exec_cavgassemble

    subroutine exec_check_2Dconv( self, cline )
        use simple_convergence, only: convergence
        use simple_parameters,  only: params_glob
        class(check_2Dconv_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)  :: params
        type(builder)     :: build
        type(convergence) :: conv
        logical :: converged
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        ! convergence check
        converged = conv%check_conv2D(cline, build%spproj_field%get_n('class'), params%msk)
        call cline%set('frac_srch', conv%get('frac_srch'))
        if( params_glob%l_doshift )then
            ! activates shift search
            call cline%set('trs', params_glob%trs)
        endif
        if( converged )then
            call cline%set('converged', 'yes')
        else
            call cline%set('converged', 'no')
        endif
        ! end gracefully
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_CHECK_2DCONV NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check_2Dconv

    subroutine exec_rank_cavgs( self, cline )
        use simple_oris, only: oris
        class(rank_cavgs_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters)     :: params
        type(builder)        :: build
        integer, allocatable :: order(:)
        real,    allocatable :: res(:)
        integer    :: ldim(3), ncls, iclass
        type(oris) :: clsdoc_ranked
        call cline%set('oritype', 'cls2D')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        call find_ldim_nptcls(params%stk, ldim, ncls)
        params%ncls = ncls
        if( build%spproj_field%get_noris() == params%ncls )then
            ! all we need to do is fetch from classdoc in projfile &
            ! order according to resolution
            call clsdoc_ranked%new(params%ncls)
            res = build%spproj%os_cls2D%get_all('res')
            allocate(order(params%ncls))
            order = (/(iclass,iclass=1,params%ncls)/)
            call hpsort(res, order)
            do iclass=1,params%ncls
                call clsdoc_ranked%set(iclass, 'class',     real(order(iclass)))
                call clsdoc_ranked%set(iclass, 'rank',      real(iclass))
                call clsdoc_ranked%set(iclass, 'pop',       build%spproj_field%get(order(iclass),  'pop'))
                call clsdoc_ranked%set(iclass, 'res',       build%spproj_field%get(order(iclass),  'res'))
                call clsdoc_ranked%set(iclass, 'corr',      build%spproj_field%get(order(iclass), 'corr'))
                call clsdoc_ranked%set(iclass, 'w',         build%spproj_field%get(order(iclass),    'w'))
                call clsdoc_ranked%set(iclass, 'specscore', build%spproj_field%get(order(iclass), 'specscore'))
                write(logfhandle,'(a,1x,i5,1x,a,1x,i5,1x,a,i5,1x,a,1x,f6.2)') 'CLASS:', order(iclass),&
                    &'RANK:', iclass ,'POP:', nint(build%spproj_field%get(order(iclass), 'pop')),&
                    &'RES:', build%spproj_field%get(order(iclass), 'res')
                call build%img%read(params%stk, order(iclass))
                call build%img%write(params%outstk, iclass)
            end do
            call clsdoc_ranked%write('classdoc_ranked.txt')
        else
            ! nothing to do
        endif
        ! end gracefully
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_RANK_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_rank_cavgs

    subroutine exec_cluster_cavgs( self, cline )
        use simple_cluster_cavgs, only: cluster_cavgs_exec
        class(cluster_cavgs_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)     :: params
        type(builder)        :: build
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir', 'yes')
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        call cluster_cavgs_exec( )
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_cluster_cavgs

    subroutine exec_write_classes( self, cline )
        class(write_classes_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        type(image)      :: img_cavg
        character(len=:),   allocatable :: cavgsstk, stkname, classname
        type(image),        allocatable :: imgs_class(:)
        real,               allocatable :: states(:), inpls(:,:)
        integer,            allocatable :: pops(:), pinds(:)
        real(kind=c_float), allocatable :: rmat_rot(:,:,:)
        integer :: ncls, n, ldim(3), icls, pop_max, ind_in_stk, i, cnt
        real    :: smpd
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call spproj%read(params%projfile)
        ! get class average stack
        call spproj%get_cavgs_stk(cavgsstk, ncls, smpd)
        call find_ldim_nptcls(cavgsstk, ldim, n)
        ldim(3) = 1
        if( n /= ncls ) THROW_HARD('Incosistent # classes in project file vs cavgs stack; exec_write_classes')
        ! get state flag array
        states = spproj%os_cls2D%get_all('state')
        ! get ncls from ptcl2D field
        n = spproj%os_ptcl2D%get_n('class')
        if( n /= ncls ) THROW_HARD('Incosistent # classes in ptcl2D field of spproj vs cavgs stack; exec_write_classes')
        ! find out maximum population and allocate image arrays accordingly
        call spproj%os_ptcl2D%get_pops(pops, 'class')
        pop_max = maxval(pops)
        write(logfhandle,'(A,I5)') '>>> MAXIMUM CLASS POPULATION: ', pop_max
        allocate(imgs_class(pop_max), inpls(pop_max,3), rmat_rot(ldim(1),ldim(2),1))
        rmat_rot = 0.
        inpls    = 0.
        do i=1,pop_max
            call imgs_class(i)%new(ldim, smpd)
        end do
        call img_cavg%new(ldim, smpd)
        ! loop over classes
        do icls=1,ncls
            if( states(icls) < 0.5 ) cycle
            ! get particle indices of class
            call spproj%os_ptcl2D%get_pinds(icls, 'class', pinds)
            if( .not. allocated(pinds) ) cycle
            ! read the class average
            call img_cavg%read(cavgsstk, icls)
            ! read the images and get the in-plane parameters
            do i=1,size(pinds)
                ! read
                call spproj%get_stkname_and_ind('ptcl2D', pinds(i), stkname, ind_in_stk)
                call imgs_class(i)%read(stkname, ind_in_stk)
                ! get params
                inpls(i,1)  = spproj%os_ptcl2D%e3get(pinds(i))
                inpls(i,2:) = spproj%os_ptcl2D%get_2Dshift(pinds(i))
            end do
            ! rotate the images (in parallel)
            !$omp parallel do default(shared) private(i,rmat_rot) schedule(static) proc_bind(close)
            do i=1,size(pinds)
                call imgs_class(i)%rtsq_serial(-inpls(i,1), -inpls(i,2), -inpls(i,3), rmat_rot)
                call imgs_class(i)%set_rmat(rmat_rot)
            end do
            !$omp end parallel do
            ! make a filename for the class
            classname = 'class'//int2str_pad(icls,5)//'.mrcs'
            ! write the class average first, followed by the rotated and shifted particles
            call img_cavg%write(classname, 1)
            cnt = 1
            do i=1,size(pinds)
                cnt = cnt + 1
                call imgs_class(i)%write(classname, cnt)
            end do
        end do
        ! destruct
        call spproj%kill
        call img_cavg%kill
        do i=1,size(imgs_class)
            call imgs_class(i)%kill
        end do
        deallocate(imgs_class, inpls, rmat_rot)
        if( allocated(states) ) deallocate(states)
        if( allocated(pops)   ) deallocate(pops)
        if( allocated(pinds)  ) deallocate(pinds)
        ! end gracefully
        call simple_end('**** SIMPLE_WRITE_CLASSES NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_write_classes

end module simple_commander_cluster2D
