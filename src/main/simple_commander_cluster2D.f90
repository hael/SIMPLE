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

public :: cleanup2D_commander_hlev
public :: cluster2D_commander_stream
public :: cluster2D_autoscale_commander_hlev
public :: cluster2D_commander_distr
public :: cluster2D_commander
public :: make_cavgs_commander_distr
public :: make_cavgs_commander
public :: cavgassemble_commander
public :: rank_cavgs_commander
public :: cluster_cavgs_commander
public :: write_classes_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: cleanup2D_commander_hlev
  contains
    procedure :: execute      => exec_cleanup2D
end type cleanup2D_commander_hlev
type, extends(commander_base) :: cluster2D_commander_stream
  contains
    procedure :: execute      => exec_cluster2D_stream
end type cluster2D_commander_stream
type, extends(commander_base) :: cluster2D_autoscale_commander_hlev
  contains
    procedure :: execute      => exec_cluster2D_autoscale
end type cluster2D_autoscale_commander_hlev
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
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',      'yes')
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
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
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
            call build%spproj_field%partition_eo
        else if( params%proj_is_class.eq.'yes' )then
            call build%spproj_field%proj2class
        endif
        ! shift multiplication
        if( params%mul > 1. )then
            call build%spproj_field%mul_shifts(params%mul)
        endif
        ! setup weights
        if( trim(params%ptclw) .eq. 'yes' )then
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
        call qsys_job_finished('simple_commander_cluster2D :: exec_make_cavgs')
        call cavger_kill
        ! end gracefully
        call build%kill_strategy2D_tbox
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_MAKE_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_make_cavgs

    subroutine exec_cleanup2D( self, cline )
        use simple_commander_project, only: scale_project_commander_distr
        use simple_procimgfile,       only: random_selection_from_imgfile, random_cls_from_imgfile
        use simple_commander_imgproc, only: scale_commander
        use simple_projection_frcs,   only: projection_frcs
        class(cleanup2D_commander_hlev), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        ! commanders
        type(cluster2D_commander_distr)     :: xcluster2D_distr
        type(scale_project_commander_distr) :: xscale_distr
        type(scale_commander)               :: xscale
        type(rank_cavgs_commander)          :: xrank_cavgs
        ! command lines
        type(cmdline)                       :: cline_cluster2D1, cline_cluster2D2
        type(cmdline)                       :: cline_rank_cavgs, cline_scale, cline_scalerefs
        ! other variables
        type(parameters)                    :: params
        type(sp_project)                    :: spproj, spproj_sc
        type(projection_frcs)               :: frcs, frcs_sc
        character(len=:),       allocatable :: projfile, orig_projfile
        character(len=LONGSTRLEN)           :: finalcavgs, finalcavgs_ranked, cavgs, refs_sc
        real                                :: scale_factor, smpd, msk, ring2, lp1, lp2
        integer                             :: last_iter, box, status
        logical                             :: do_scaling, do_ranking
        ! parameters
        character(len=STDLEN) :: orig_projfile_bak = 'orig_bak.simple'
        integer, parameter    :: MINBOX      = 92
        real,    parameter    :: TARGET_LP   = 15.
        real,    parameter    :: MINITS      =  5., MINITS_FAST =  9.
        real,    parameter    :: MAXITS      = 15., MAXITS_FAST = 18.
        real                  :: SMPD_TARGET = 4.
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',     'yes')
        if( .not. cline%defined('lp')        ) call cline%set('lp',         15. )
        if( .not. cline%defined('ncls')      ) call cline%set('ncls',      200. )
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      20. )
        if( .not. cline%defined('center')    ) call cline%set('center',     'no')
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('refine')    ) call cline%set('refine', 'greedy')
        if( .not. cline%defined('oritype')   ) call cline%set('oritype', 'ptcl2D')
        do_ranking = .true.
        if( cline%defined('stream') )then
            if( cline%get_carg('stream').eq.'yes' ) do_ranking = .false.
        endif
        call cline%set('stream', 'no')
        call params%new(cline)
        orig_projfile = params%projfile
        if( .not. cline%defined('maxits') )then
            if( params%refine.eq.'fast' )then
                params%maxits = nint(MAXITS_FAST)
                call cline%set('maxits', MAXITS_FAST)
            else
                params%maxits = nint(MAXITS)
                call cline%set('maxits', MAXITS)
            endif
        endif
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
        ! down-scaling for fast execution, greedy optimisation, no match filter
        ! no incremental learning, objective function is standard cross-correlation (cc), no centering
        cline_cluster2D1 = cline
        cline_cluster2D2 = cline
        cline_scale      = cline
        call cline_cluster2D1%set('prg', 'cluster2D')
        if( params%refine.eq.'fast' )then
            call cline_cluster2D1%set('maxits', MINITS_FAST)
        else
            call cline_cluster2D1%set('maxits', MINITS)
        endif
        call cline_cluster2D1%set('objfun',     'cc')
        call cline_cluster2D1%set('match_filt', 'no')
        call cline_cluster2D1%set('center',     'no')
        call cline_cluster2D1%set('autoscale',  'no')
        call cline_cluster2D1%set('ptclw',      'no')
        call cline_cluster2D1%delete('update_frac')
        ! second stage
        ! down-scaling for fast execution, greedy optimisation, no match filter
        ! objective function default is standard cross-correlation (cc)
        call cline_cluster2D2%set('prg', 'cluster2D')
        call cline_cluster2D2%set('match_filt', 'no')
        call cline_cluster2D2%set('autoscale',  'no')
        call cline_cluster2D2%set('trs',         MINSHIFT)
        call cline_cluster2D2%set('objfun',     'cc')
        if( .not.cline%defined('maxits') )then
            if( params%refine.eq.'fast' )then
                call cline_cluster2D2%set('maxits', MAXITS_FAST)
            else
                call cline_cluster2D2%set('maxits', MAXITS)
            endif
        endif
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
        if( cline%defined('refs') )then
            if( params%autoscale.eq.'yes')then
                call cline_scalerefs%set('stk', params%refs)
                refs_sc = 'refs'//trim(SCALE_SUFFIX)//'.mrc'
                call cline_scalerefs%set('outstk', trim(refs_sc))
                call cline_scalerefs%set('smpd',   smpd)
                call cline_scalerefs%set('newbox', real(box))
                call xscale%execute(cline_scalerefs)
                params%refs = trim(refs_sc)
            else
                ! all good
            endif
        else
            params%refs = 'start2Drefs' // params%ext
            call spproj%read(projfile)
            if( params%avg.eq.'yes' )then
                call random_cls_from_imgfile(spproj, params%refs, params%ncls)
            else
                call random_selection_from_imgfile(spproj, params%refs, box, params%ncls)
            endif
            call spproj%kill
        endif
        ! updates command-lines
        call cline_cluster2D1%set('refs',   params%refs)
        call cline_cluster2D1%set('msk',    msk)
        call cline_cluster2D1%set('ring2',  ring2)
        call cline_cluster2D1%set('lp',     lp1)
        call cline_cluster2D2%set('msk',    msk)
        call cline_cluster2D2%set('lp',     lp2)
        if( trim(params%pssnr).eq.'yes' )then
            call cline_cluster2D2%set('match_filt', 'yes')
            call cline_cluster2D2%set('pssnr', 'yes')
            call cline_cluster2D2%delete('lp')
        endif
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
        call cline%set('endit',real(last_iter))
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
            call frcs_sc%read(FRCS_FILE)
            call frcs_sc%upsample(params%smpd, params%box, frcs)
            call frcs%write(FRCS_FILE)
            call spproj%add_frcs2os_out(FRCS_FILE, 'frc2D')
            call frcs%kill
            call frcs_sc%kill
            spproj%os_ptcl2D = spproj_sc%os_ptcl2D
            spproj%os_cls2D  = spproj_sc%os_cls2D
            spproj%os_cls3D  = spproj_sc%os_cls3D
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
        if( do_ranking )then
            finalcavgs_ranked = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter,3)//'_ranked'//params%ext
            call cline_rank_cavgs%set('projfile', trim(params%projfile))
            call cline_rank_cavgs%set('stk',      trim(finalcavgs))
            call cline_rank_cavgs%set('outstk',   trim(finalcavgs_ranked))
            call xrank_cavgs%execute(cline_rank_cavgs)
        endif
        ! cleanup
        if( do_scaling ) call simple_rmdir(STKPARTSDIR)
        ! end gracefully
        call simple_end('**** SIMPLE_CLEANUP2D NORMAL STOP ****', print_simple=do_ranking)
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
        class(cluster2D_commander_stream), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        integer,               parameter   :: WAIT_WATCHER       = 30    ! seconds prior to new stack detection
        integer,               parameter   :: ORIGPROJ_WRITEFREQ = 7200  ! Frequency at which the original project file should be updated
        integer,               parameter   :: MINBOXSZ           = 72    ! minimum boxsize for scaling
        ! integer,               parameter   :: ORIGPROJ_WRITEFREQ = 1800 ! dev settings
        real,                  parameter   :: GREEDY_TARGET_LP   = 15.
        real                               :: SMPD_TARGET        = 4.    ! target sampling distance
        character(len=STDLEN), parameter   :: USER_PARAMS        = 'stream2D_user_params.txt'
        character(len=STDLEN), parameter   :: SPPROJ_SNAPSHOT    = 'SIMPLE_PROJECT_SNAPSHOT'
        character(len=STDLEN), parameter   :: PROJFILE_BUFFER    = 'buffer.simple'
        character(len=STDLEN), parameter   :: PROJFILE_POOL        = 'cluster2D.simple'
        character(len=STDLEN), parameter   :: SCALE_DIR          = './scaled_stks/'
        logical,               parameter   :: debug_here = .false.
        type(parameters)                   :: params
        type(make_cavgs_commander_distr)   :: xmake_cavgs
        type(rank_cavgs_commander)         :: xrank_cavgs
        type(cmdline)                      :: cline_cluster2D, cline_cluster2D_buffer, cline_rank_cavgs, cline_make_cavgs
        type(sp_project)                   :: orig_proj, stream_proj, buffer_proj, pool_proj
        character(LONGSTRLEN), allocatable :: spproj_list(:), new_buffer_spprojs(:), imported_spprojs(:), imported_stks(:)
        character(len=:),      allocatable :: spproj_list_fname, orig_projfile
        character(len=STDLEN)              :: str_iter, refs_glob, refs_glob_ranked
        character(len=LONGSTRLEN)          :: buffer_dir, prev_snapshot_frcs, prev_snapshot_cavgs
        real    :: orig_smpd, msk, scale_factor, orig_msk, smpd, large_msk, lp_greedy, lpstart_stoch
        integer :: nptcls_per_buffer, nmics_imported, nstks_imported ! target number of particles per buffer, # of projects & stacks in pool
        integer :: n_buffer_spprojs, n_buffer_stks, n_buffer_ptcls   ! # of projects, stacks & particles, in new buffer
        integer :: n_rejected, n_transfers, ncls_glob, last_injection, max_ncls, buffer_id
        integer :: iter, orig_box, box, boxpd, n_spprojs, pool_iter, origproj_time, i
        logical :: do_autoscale, l_maxed, buffer_exists, do_wait, l_greedy, l_forced_exit
        if( cline%defined('refine') )then
            if( trim(cline%get_carg('refine')).ne.'greedy' )then
                if( .not.cline%defined('msk') ) THROW_HARD('MSK must be defined!')
            endif
        else
            if( .not.cline%defined('msk') ) THROW_HARD('MSK must be defined!')
        endif
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',     'yes')
        if( .not. cline%defined('lp')        ) call cline%set('lp',          15.)
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',       20.)
        if( .not. cline%defined('center')    ) call cline%set('center',     'no')
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('lpthresh')  ) call cline%set('lpthresh',    30.)
        if( .not. cline%defined('ndev')      ) call cline%set('ndev',        1.5)
        if( .not. cline%defined('oritype')   ) call cline%set('oritype','ptcl2D')
        if( .not. cline%defined('ptclw')     ) call cline%set('ptclw',      'no')
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
        do_autoscale      = params%autoscale .eq. 'yes'
        max_ncls          = floor(real(params%ncls)/real(params%ncls_start))*params%ncls_start ! effective maximum # of classes
        nptcls_per_buffer = params%nptcls_per_cls*params%ncls_start         ! # of particles in each buffer
        buffer_exists     = .false.                                         ! whether the buffer exists
        l_maxed           = .false.                                         ! whether all chunks have been merged
        do_wait           = .true.
        l_forced_exit     = .false.
        spproj_list_fname = filepath(trim(params%dir_target), trim(STREAM_SPPROJFILES))
        ncls_glob         = 0
        n_transfers       = 0
        nmics_imported    = 0
        nstks_imported    = 0
        n_rejected        = 0
        buffer_id         = 0
        l_greedy          = trim(params%refine).eq.'greedy'
        lp_greedy         = GREEDY_TARGET_LP
        if( cline%defined('lp') ) lp_greedy = params%lp
        SMPD_TARGET       = lp_greedy/2.
        lpstart_stoch     = (GREEDY_TARGET_LP+lp_greedy)/2.
        prev_snapshot_cavgs = ''
        call write_user_params
        ! Automated exit after 2 hours without new particles
        if(.not.cline%defined('time_inactive')) params%time_inactive = 2*60
        ! init command-lines
        call cline%delete('lp')
        call cline%delete('refine')
        cline_cluster2D         = cline
        cline_cluster2D_buffer  = cline
        ! buffer classification
        ! Greedy optimisation, no match filter, no centering
        ! no incremental learning, objective function default is standard cross-correlation (cc)
        call cline_cluster2D_buffer%set('prg',       'cluster2D')
        call cline_cluster2D_buffer%set('projfile',  trim(PROJFILE_BUFFER))
        call cline_cluster2D_buffer%set('projname',  trim(get_fbody(trim(PROJFILE_BUFFER),trim('simple'))))
        call cline_cluster2D_buffer%set('objfun',    'cc')
        call cline_cluster2D_buffer%set('center',    'no')
        call cline_cluster2D_buffer%set('match_filt','no')
        call cline_cluster2D_buffer%set('autoscale', 'no')
        call cline_cluster2D_buffer%set('refine',    'snhc')
        call cline_cluster2D_buffer%set('ptclw',     'no')
        call cline_cluster2D_buffer%delete('stream')
        call cline_cluster2D_buffer%delete('update_frac')
        if( l_greedy )then
            call cline_cluster2D_buffer%set('maxits', 10.)
            call cline_cluster2D_buffer%set('refine', 'greedy')
            call cline_cluster2D_buffer%set('lp',     GREEDY_TARGET_LP)
        else
            call cline_cluster2D_buffer%set('maxits',  12.) ! guaranties 3 iterations with withdrawal
            call cline_cluster2D_buffer%set('lpstart', GREEDY_TARGET_LP)
            call cline_cluster2D_buffer%set('lpstop',  lp_greedy)
        endif
        ! pool classification: optional stochastic optimisation, optional match filter
        ! automated incremental learning, objective function is standard cross-correlation (cc)
        call cline_cluster2D%set('prg',       'cluster2D')
        call cline_cluster2D%set('autoscale', 'no')
        call cline_cluster2D%set('trs',       MINSHIFT)
        call cline_cluster2D%set('projfile',  trim(PROJFILE_POOL))
        call cline_cluster2D%set('projname',  trim(get_fbody(trim(PROJFILE_POOL),trim('simple'))))
        call cline_cluster2D%set('objfun',    'cc')
        if( .not.cline%defined('match_filt') ) call cline_cluster2D%set('match_filt','no')
        if( l_greedy )then
            call cline_cluster2D%set('center', 'no')
            call cline_cluster2D%set('refine', 'greedy')
            call cline_cluster2D%set('lp',     lp_greedy)
        else
            call cline_cluster2D%set('lpstart', lpstart_stoch)
        endif
        ! WAIT FOR FIRST STACKS
        n_buffer_spprojs = 0
        n_buffer_stks    = 0
        do
            if( file_exists(spproj_list_fname) )then
                if( .not.is_file_open(spproj_list_fname) )then
                    call read_mics(new_buffer_spprojs, n_buffer_spprojs, n_buffer_stks, n_buffer_ptcls)
                    if( n_buffer_spprojs > 0 )then
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
        call orig_proj%kill
        ! getting general parameters from the first sp_project
        orig_msk = params%msk
        orig_box = 0
        do i=1,size(spproj_list)
            call stream_proj%read_segment('mic',spproj_list(i))
            if( stream_proj%os_mic%get(1,'nptcls') > 0.5 )then
                call stream_proj%read_segment('stk',spproj_list(i))
                orig_box  = stream_proj%get_box()
                orig_smpd = stream_proj%get_smpd()
                exit
            endif
        enddo
        if( orig_box == 0 ) THROW_HARD('FATAL ERROR')
        call stream_proj%kill
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
        boxpd     = 2*round2even(params%alpha*real(box/2)) ! logics from parameters
        large_msk = max(msk, real(box/2)-COSMSKHALFWIDTH)
        if( l_greedy )then
            if( .not.cline%defined('msk') )then
                call cline_cluster2D%set('msk',large_msk)
            else
                call cline_cluster2D%set('msk',msk)
            endif
        else
            call cline_cluster2D%set('msk',    msk)
        endif
        ! First import
        call update_pool
        ! generates buffer project
        call gen_buffer_from_pool
        ! MAIN LOOP
        last_injection = simple_gettime()
        origproj_time  = last_injection
        pool_iter      = 1
        do iter = 1,9999
            str_iter  = int2str_pad(iter,3)
            pool_iter = min(999,pool_iter)
            write(logfhandle,'(A,I4)')'>>> WAIT CYCLE ',iter
            if( is_timeout(simple_gettime()) )exit
            if( buffer_exists )then
                call classify_buffer
                ! particles selection
                call reject_from_buffer
                ! book keeping
                call transfer_buffer_to_pool
                ! cleanup
                call buffer_proj%kill
                call simple_rmdir(buffer_dir)
                ! updates global variables
                buffer_exists  = .false.
                last_injection = simple_gettime()
            endif
            ! one iteration of the whole dataset when not converged or each 5 iterations
            if( .not.cline_cluster2D%defined('converged') .or. mod(iter,5)==0 )then
                call classify_pool
                pool_iter = pool_iter+1
                ! rejection each 5 iterations
                if( mod(iter,5)==0 ) call reject_from_pool
                do_wait = .false.
            else if( mod(iter,5) /= 0 )then
                do_wait = .true.
            endif
            call write_snapshot( .false., .true.)
            ! termination and/or pause
            do while( file_exists(trim(PAUSE_STREAM)) )
                if( file_exists(trim(TERM_STREAM)) ) exit
                call write_singlelineoftext(PAUSE_STREAM, 'PAUSED')
                write(logfhandle,'(A,A)')'>>> CLUSTER2D STREAM PAUSED ',cast_time_char(simple_gettime())
                call simple_sleep(WAIT_WATCHER)
            enddo
            if( file_exists(TERM_STREAM) )then
                write(logfhandle,'(A,A)')'>>> TERMINATING CLUSTER2D STREAM ',cast_time_char(simple_gettime())
                l_forced_exit = .true.
                exit
            endif
            ! detect new project files
            call read_mics(new_buffer_spprojs, n_buffer_spprojs, n_buffer_stks, n_buffer_ptcls)
            do_wait = do_wait .and. (n_buffer_spprojs == 0)
            if( do_wait ) call simple_sleep(WAIT_WATCHER)
            ! update original project
            if( simple_gettime()-origproj_time > ORIGPROJ_WRITEFREQ )then
                call write_snapshot( .true., .true.)
                origproj_time = simple_gettime()
            endif
            ! updates pool with new projects
            call update_pool
            ! optionally builds new buffer
            call gen_buffer_from_pool
        enddo
        call qsys_cleanup(keep2D=.false.)
        ! updates original project
        call write_snapshot( .true., .false.)
        ! whether to generate class averages
        if( .not.l_forced_exit )then
            call cline_make_cavgs%set('prg',     'make_cavgs')
            call cline_make_cavgs%set('mkdir',   'no')
            call cline_make_cavgs%set('ncls',    real(ncls_glob))
            call cline_make_cavgs%set('nparts',  real(params%nparts))
            call cline_make_cavgs%set('nthr',    real(params%nthr))
            call cline_make_cavgs%set('projfile',orig_projfile)
            call cline_make_cavgs%set('msk',     orig_msk)
            call cline_make_cavgs%set('ptclw',   'no')
            call cline_make_cavgs%set('refs',    refs_glob)
            call xmake_cavgs%execute(cline_make_cavgs)
        endif
        ! ranking
        refs_glob_ranked = add2fbody(refs_glob,params%ext,'_ranked')
        call cline_rank_cavgs%set('projfile', orig_projfile)
        call cline_rank_cavgs%set('stk',      refs_glob)
        call cline_rank_cavgs%set('outstk',   trim(refs_glob_ranked))
        call xrank_cavgs%execute(cline_rank_cavgs)
        ! cleanup
        if( .not.debug_here )then
            call qsys_cleanup
            call simple_rmdir(SCALE_DIR)
            call del_file(PROJFILE_POOL)
            call del_file(PROJFILE_BUFFER)
            call simple_rmdir(buffer_dir)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER2D_STREAM NORMAL STOP ****')
        contains

            ! optimized 2D classification commander
            subroutine exec_cluster2D_stream_distr( cline, spproj )
                use simple_timer
                class(cmdline),    intent(inout) :: cline
                class(sp_project), intent(inout) :: spproj
                logical, parameter :: L_BENCH = .false.
                ! command lines
                type(cmdline) :: cline_check_2Dconv
                type(cmdline) :: cline_cavgassemble
                ! other variables
                type(parameters)          :: params
                type(qsys_env)            :: qenv
                type(chash)               :: job_descr
                type(sp_project)          :: transfer_spproj
                logical, allocatable      :: transfer_mask(:)
                integer, allocatable      :: update_cnts(:), eo_pops(:,:)
                character(len=LONGSTRLEN) :: refs, str_iter
                real                      :: srch_frac
                integer                   :: iter, cnt, iptcl, nptcls, nptcls_old, n2update, ntot, indinstk, stkind, eo, icls
                integer(timer_int_kind)   :: t_init,  t,  t_tot
                if( L_BENCH )then
                    t_init = tic()
                    t_tot  = t_init
                endif
                ! sanity checks
                if( spproj%get_nptcls() == 0 )    THROW_HARD('no particles found! exec_cluster2D_stream_distr')
                if( .not. cline%defined('refs') ) THROW_HARD('refs is undefined!')
                if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',   30.)
                call cline%set('oritype', 'ptcl2D')
                call cline%set('mkdir', 'no')
                call cline%set('stream','yes')
                call cline%set('extr_iter', 100.) ! cooling of the randomization rate
                call cline%delete('update_frac')
                ! taking care of particles history
                ntot = spproj%os_ptcl2D%get_noris()
                allocate(transfer_mask(ntot))
                allocate(update_cnts(ntot), source=-1)
                nptcls_old = 0
                do iptcl = 1,ntot
                    transfer_mask(iptcl) = spproj%os_ptcl2D%get(iptcl,'state') > 0.5
                    if( transfer_mask(iptcl) ) update_cnts(iptcl) = nint(spproj%os_ptcl2D%get(iptcl,'updatecnt'))
                enddo
                nptcls     = count(transfer_mask)                   ! state=1
                nptcls_old = count(update_cnts >= STREAM_SRCHLIM)   ! old state=1
                if( nptcls_old > 0 )then
                    call cline%set('center', 'no') ! deactivates centering
                    srch_frac = STREAM_SRCHFRAC
                    if( nint(real(nptcls_old)*STREAM_SRCHFRAC) > MAX_STREAM_NPTCLS )then
                        ! cap reached
                        srch_frac = max(0.05,(real(MAX_STREAM_NPTCLS) / real(nptcls_old)))
                    endif
                    ! randomly deactivating 'old' particles & calculate deselected class populations for cavg_assemble
                    allocate(eo_pops(ncls_glob,2),source=0)
                    do iptcl = 1,ntot
                        if( transfer_mask(iptcl) )then
                            if( update_cnts(iptcl) >= STREAM_SRCHLIM )then
                                transfer_mask(iptcl) = ran3() < srch_frac
                                if( .not.transfer_mask(iptcl) )then
                                    icls = nint(spproj%os_ptcl2D%get(iptcl,'class'))
                                    eo   = nint(spproj%os_ptcl2D%get(iptcl,'eo')) + 1
                                    eo_pops(icls,eo) = eo_pops(icls,eo) + 1
                                endif
                            endif
                        endif
                    end do
                    ! temporary project for selected particles
                    n2update = count(transfer_mask)
                    transfer_spproj%projinfo = spproj%projinfo
                    transfer_spproj%compenv  = spproj%compenv
                    if( spproj%jobproc%get_noris()>0 ) transfer_spproj%jobproc = spproj%jobproc
                    call transfer_spproj%projinfo%delete_entry('projname')
                    call transfer_spproj%projinfo%delete_entry('projfile')
                    call transfer_spproj%update_projinfo( cline )
                    transfer_spproj%os_stk   = spproj%os_stk
                    transfer_spproj%os_cls2D = spproj%os_cls2D
                    do icls = 1,ncls_glob
                        call transfer_spproj%os_cls2D%set(icls,'pop_even',real(eo_pops(icls,1)))
                        call transfer_spproj%os_cls2D%set(icls,'pop_odd', real(eo_pops(icls,2)))
                    enddo
                    deallocate(eo_pops)
                    call transfer_spproj%os_ptcl2D%new(n2update)
                    call transfer_spproj%os_ptcl3D%new(n2update)
                    cnt = 0
                    do iptcl = 1,ntot
                        if( transfer_mask(iptcl) )then
                            cnt = cnt+1
                            call transfer_spproj%os_ptcl2D%transfer_ori(cnt, spproj%os_ptcl2D, iptcl)
                            call spproj%map_ptcl_ind2stk_ind('ptcl2D', iptcl, stkind, indinstk)
                            call transfer_spproj%os_ptcl2D%set(cnt,'indstk',real(indinstk)) ! to bypass stack indexing convention
                            call transfer_spproj%os_ptcl3D%set(cnt,'state',1.)
                        endif
                    enddo
                    call transfer_spproj%write(PROJFILE_POOL)
                    call transfer_spproj%os_ptcl3D%kill
                    call transfer_spproj%os_stk%kill
                    call transfer_spproj%os_ptcl2D%new(n2update)
                else
                    call spproj%write(PROJFILE_POOL)
                    n2update = ntot
                endif
                deallocate(update_cnts)
                ! params, project is passed, not read
                call cline%set('stream',      'yes')
                call cline%set('box',         real(box))
                call cline%set('smpd',        smpd)
                call cline%set('nptcls',      real(n2update))
                call params%new(cline)
                params%update_frac = 1.
                if( nptcls_old > 0 )then
                    params_glob%nptcls = n2update
                    params%update_frac = max(0.1,real(params%nptcls) / real(nptcls))
                endif
                call cline%set('update_frac', params%update_frac)
                ! prepare command lines from prototype
                cline_check_2Dconv = cline
                cline_cavgassemble = cline
                ! setup the environment for distributed execution
                call qenv%new(params%nparts)
                ! prepare job description
                call cline%gen_job_descr(job_descr)
                ! initialise static command line parameters and static job description parameters
                call cline_cavgassemble%set('prg',   'cavgassemble')
                ! single iteration
                iter = params%startit
                str_iter = int2str_pad(iter,3)
                write(logfhandle,'(A)')   '>>>'
                write(logfhandle,'(A,I6)')'>>> ITERATION ', iter
                write(logfhandle,'(A)')   '>>>'
                ! updates
                call job_descr%set('refs', trim(params%refs))
                call job_descr%set('startit', int2str(iter))
                ! the only FRC we have is from the previous iteration, hence the iter - 1
                call job_descr%set('frcs', trim(FRCS_FILE))
                if( L_BENCH )then
                    print *,'timer init        : ',toc(t_init)
                    t = tic()
                endif
                ! schedule
                call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=trim(ALGN_FBODY))
                if( L_BENCH )then
                    print *,'timer iter        : ',toc(t)
                    t = tic()
                endif
                ! merge orientation documents and transfers to memory
                if( nptcls_old > 0 )then
                    call transfer_spproj%merge_algndocs(params%nptcls, params%nparts, 'ptcl2D', ALGN_FBODY)
                    ! transfers alignment parameters back
                    cnt = 0
                    do iptcl = 1,ntot
                        if( transfer_mask(iptcl) )then
                            cnt = cnt+1
                            call spproj%os_ptcl2D%transfer_2Dparams(iptcl, transfer_spproj%os_ptcl2D, cnt)
                            call spproj%os_ptcl2D%set(iptcl, 'dist_inpl',  transfer_spproj%os_ptcl2D%get(cnt, 'dist_inpl'))
                            call spproj%os_ptcl2D%set(iptcl, 'mi_class',   transfer_spproj%os_ptcl2D%get(cnt, 'mi_class'))
                        endif
                    enddo
                    call transfer_spproj%kill
                else
                    call spproj%merge_algndocs(params%nptcls, params%nparts, 'ptcl2D', ALGN_FBODY)
                endif
                if( L_BENCH )then
                    print *,'timer merge       : ',toc(t)
                    t = tic()
                endif
                ! assemble class averages & update in memory
                refs = trim(CAVGS_ITER_FBODY) // trim(str_iter) // params%ext
                call cline_cavgassemble%set('refs', trim(refs))
                call qenv%exec_simple_prg_in_queue(cline_cavgassemble, 'CAVGASSEMBLE_FINISHED')
                call spproj%read_segment('cls2D', params_glob%projfile)
                call spproj%read_segment('out', params_glob%projfile)
                if( L_BENCH )then
                    print *,'timer cavgassemble: ',toc(t)
                    t = tic()
                endif
                ! check convergence
                call check_2Dconv(cline_check_2Dconv, spproj%os_ptcl2D)
                if( iter == 1 )then
                    call cline%set('converged','no')
                else
                    if( cline_check_2Dconv%get_carg('converged').eq.'yes' )then
                        if( cline_check_2Dconv%get_carg('converged').eq.'yes' )call cline%set('converged','yes')
                    endif
                endif
                if( L_BENCH )then
                    print *,'timer convergence : ',toc(t)
                endif
                ! report the last iteration on exit
                call cline%delete( 'startit' )
                call cline%set('endit', real(iter))
                ! restores global parameters
                params_glob%nptcls = ntot
                ! end gracefully
                call qsys_cleanup
                call simple_end('**** SIMPLE_DISTR_CLUSTER2D NORMAL STOP ****',print_simple=.false.)
                if( L_BENCH ) print *,'timer tot         : ',toc(t_tot)
            end subroutine exec_cluster2D_stream_distr

            subroutine rescale_cavgs( src, dest )
                character(len=*), intent(in) :: src, dest
                type(image) :: img, img_pad
                integer, allocatable :: cls_pop(:)
                integer     :: icls, iostat
                call img%new([box,box,1],smpd)
                call img_pad%new([orig_box,orig_box,1],params%smpd)
                cls_pop = nint(pool_proj%os_cls2D%get_all('pop'))
                do icls = 1,ncls_glob
                    if( cls_pop(icls) > 0 )then
                        call img%zero_and_unflag_ft
                        call img%read(src,icls)
                        call img%fft
                        call img%pad(img_pad, backgr=0.)
                        call img_pad%ifft
                    else
                        img_pad = 0.
                    endif
                    call img_pad%write('tmp_cavgs.mrc',icls)
                enddo
                iostat = simple_rename('tmp_cavgs.mrc',dest)
                call img%kill
                call img_pad%kill
            end subroutine rescale_cavgs

            subroutine reject_from_buffer
                type(image)          :: img
                logical, allocatable :: cls_mask(:)
                character(STDLEN)    :: refs_buffer
                integer              :: nptcls_rejected, ncls_rejected, iptcl
                integer              :: boxmatch, icls, endit, ncls_here, cnt
                if( debug_here ) print *,'in reject from_buffer'; call flush(6)
                call update_user_params
                ncls_rejected   = 0
                nptcls_rejected = 0
                ncls_here       = buffer_proj%os_cls2D%get_noris()
                boxmatch        = find_boxmatch(box, msk)
                endit           = nint(cline_cluster2D_buffer%get_rarg('endit'))
                refs_buffer     = trim(buffer_dir)//'/'//trim(CAVGS_ITER_FBODY)//int2str_pad(endit,3)//params%ext
                allocate(cls_mask(ncls_here), source=.true.)
                ! print out prior to rejection
                if( debug_here ) call buffer_proj%os_cls2D%write('classdoc_buffer_'//int2str_pad(iter,3)//'.txt')
                ! resolution and correlation
                call buffer_proj%os_cls2D%find_best_classes(boxmatch,smpd,params%lpthresh,cls_mask,params%ndev)
                if( any(cls_mask) )then
                    ncls_rejected = count(.not.cls_mask)
                    ! rejects particles 2D/3D
                    do iptcl=1,buffer_proj%os_ptcl2D%get_noris()
                        if( buffer_proj%os_ptcl2D%get_state(iptcl) == 0 )cycle
                        icls = nint(buffer_proj%os_ptcl2D%get(iptcl,'class'))
                        if( cls_mask(icls) ) cycle
                        nptcls_rejected = nptcls_rejected+1
                        call buffer_proj%os_ptcl2D%set(iptcl,'state',0.)
                        call buffer_proj%os_ptcl3D%set(iptcl,'state',0.)
                    enddo
                    ! updates cls2D field
                    do icls=1,ncls_here
                        if( .not.cls_mask(icls) )then
                            call buffer_proj%os_cls2D%set(icls,'pop',  0.)
                            call buffer_proj%os_cls2D%set(icls,'corr',-1.)
                        endif
                    enddo
                    ! updates class averages
                    call img%new([box,box,1],smpd)
                    cnt = 0
                    do icls=1,ncls_here
                        if( cls_mask(icls) ) cycle
                        cnt = cnt+1
                        if( debug_here )then
                            call img%read(refs_buffer,icls)
                            call img%write('rejected_'//int2str(pool_iter)//'.mrc',cnt)
                        endif
                        img = 0.
                        call img%write(refs_buffer,icls)
                    enddo
                    call img%read(refs_buffer, params%ncls_start)
                    call img%write(refs_buffer,params%ncls_start)
                    n_rejected = n_rejected + nptcls_rejected
                    write(logfhandle,'(A,I6,A,I6,A)')'>>> REJECTED FROM BUFFER: ',nptcls_rejected,' PARTICLES IN ',ncls_rejected,' CLUSTERS'
                    if( debug_here ) call buffer_proj%os_cls2D%write('classdoc_buffer_aftersel'//trim(str_iter)//'.txt')
                endif
                call img%kill
                if( debug_here ) print *,'end reject from_buffer'; call flush(6)
            end subroutine reject_from_buffer

            subroutine reject_from_pool
                type(image)          :: img
                logical, allocatable :: cls_mask(:)
                real                 :: ndev_here
                integer              :: nptcls_rejected, ncls_rejected, iptcl
                integer              :: boxmatch, icls, cnt
                if( debug_here ) print *,'in reject from_pool'; call flush(6)
                call update_user_params
                ncls_rejected   = 0
                nptcls_rejected = 0
                boxmatch        = find_boxmatch(box, msk)
                allocate(cls_mask(ncls_glob), source=.true.)
                if( debug_here )call pool_proj%os_cls2D%write('classdoc_pool_beforesel_'//int2str(pool_iter)//'.txt')
                ! correlation & resolution
                ndev_here = 1.5*params%ndev ! less stringent rejection
                call pool_proj%os_cls2D%find_best_classes(boxmatch,smpd,params%lpthresh,cls_mask,ndev_here)
                if( .not.all(cls_mask) )then
                    ncls_rejected = 0
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
                                ncls_rejected = ncls_rejected+1
                            endif
                        enddo
                        cnt = 0
                        call img%new([box,box,1],smpd)
                        do icls=1,ncls_glob
                            if( cls_mask(icls) ) cycle
                            cnt = cnt+1
                            if( debug_here )then
                                call img%read(refs_glob,icls)
                                call img%write('rejected_pool_'//int2str(pool_iter)//'.mrc',cnt)
                            endif
                            img = 0.
                            call img%write(refs_glob,icls)
                        enddo
                        call img%read(refs_glob, ncls_glob)
                        call img%write(refs_glob, ncls_glob)
                        call img%kill
                        deallocate(cls_mask)
                        n_rejected = n_rejected + nptcls_rejected
                        write(logfhandle,'(A,I4,A,I6,A)')'>>> REJECTED FROM POOL: ',nptcls_rejected,' PARTICLES IN ',ncls_rejected,' CLUSTERS'
                        if( debug_here )call pool_proj%os_cls2D%write('classdoc_pool_aftersel_'//int2str(pool_iter)//'.txt')
                    endif
                else
                    write(logfhandle,'(A,I4,A,I6,A)')'>>> NO PARTICLES FLAGGED FOR REJECTION FROM POOL'
                endif
                call img%kill
                if( debug_here ) print *,'end reject from_pool'; call flush(6)
            end subroutine reject_from_pool

            ! import new micrographs, stacks & particles to pool project
            subroutine update_pool
                type(sp_project) :: stream_proj
                type(oris)       :: os_stk
                character(LONGSTRLEN), allocatable :: stk_list(:), tmp(:)
                integer,               allocatable :: nptcls(:)
                real    :: dfx,dfy,angast,phshift
                integer :: box_coords(2),iproj,istk,i,n_newstks,iptcl,n_prevptcls
                logical :: transfer_ctf_params
                if( debug_here ) print *,'start update_pool'; call flush(6)
                if( n_buffer_spprojs == 0 )return
                n_prevptcls    = pool_proj%os_ptcl2D%get_noris()
                nmics_imported = pool_proj%os_mic%get_noris()
                ! add movies
                if( nmics_imported == 0 )then
                    call pool_proj%os_mic%new(n_buffer_spprojs)
                else
                    call pool_proj%os_mic%reallocate(nmics_imported+n_buffer_spprojs)
                endif
                allocate(nptcls(n_buffer_spprojs),source=0)
                do iproj=1,n_buffer_spprojs
                    call stream_proj%read_segment('mic', new_buffer_spprojs(iproj))
                    call pool_proj%os_mic%transfer_ori(nmics_imported+iproj, stream_proj%os_mic, 1)
                    nptcls(iproj) = nint(stream_proj%os_mic%get(1,'nptcls'))
                enddo
                call stream_proj%kill
                n_newstks = count(nptcls>0)
                if( n_newstks /= n_buffer_stks )THROW_HARD('incompatible # of stacks; update_pool')
                nmics_imported = pool_proj%os_mic%get_noris()
                ! fetch stack list
                allocate(stk_list(n_buffer_stks))
                call os_stk%new(n_buffer_stks)
                istk = 0
                do iproj=1,n_buffer_spprojs
                    if( nptcls(iproj) == 0 )cycle
                    istk = istk + 1
                    call stream_proj%read_segment('stk', new_buffer_spprojs(iproj))
                    stk_list(istk) = trim(stream_proj%get_stkname(1))
                    call os_stk%transfer_ori(istk, stream_proj%os_stk, 1)
                enddo
                call stream_proj%kill
                ! updates stacks list & scales stacks
                if( allocated(imported_stks) )then
                    tmp = imported_stks(:)
                    deallocate(imported_stks); allocate(imported_stks(nstks_imported+n_buffer_stks))
                    imported_stks(1:nstks_imported) = tmp(:)
                    imported_stks(nstks_imported+1:nstks_imported+n_buffer_stks) = stk_list(:)
                else
                    imported_stks  = stk_list
                endif
                call scale_stks( stk_list ) ! must come first as names updated
                if( do_autoscale) call os_stk%set_all2single('smpd',smpd)
                call pool_proj%add_stktab(stk_list, os_stk)
                nstks_imported = pool_proj%os_stk%get_noris()
                ! transfer picking coordinates & ctf parameters; needs to be included in add_stktab ??
                iptcl = n_prevptcls
                do istk=1,n_buffer_stks
                    call stream_proj%read_segment('ptcl2D', new_buffer_spprojs(istk))
                    transfer_ctf_params = stream_proj%os_ptcl2D%isthere('dfx')
                    do i = 1,stream_proj%os_ptcl2D%get_noris()
                        iptcl = iptcl+1
                        if( stream_proj%has_boxcoords(i) )then
                            call stream_proj%get_boxcoords(i,box_coords)
                            call pool_proj%set_boxcoords(iptcl, box_coords)
                        endif
                        if( transfer_ctf_params )then
                            dfx = stream_proj%os_ptcl2D%get(i,'dfx')
                            dfy = stream_proj%os_ptcl2D%get(i,'dfy')
                            angast = stream_proj%os_ptcl2D%get(i,'angast')
                            call pool_proj%os_ptcl2D%set(iptcl, 'dfx',   dfx)
                            call pool_proj%os_ptcl2D%set(iptcl, 'dfy',   dfy)
                            call pool_proj%os_ptcl2D%set(iptcl, 'angast',angast)
                            call pool_proj%os_ptcl3D%set(iptcl, 'dfx',   dfx)
                            call pool_proj%os_ptcl3D%set(iptcl, 'dfy',   dfy)
                            call pool_proj%os_ptcl3D%set(iptcl, 'angast',angast)
                            if( stream_proj%os_ptcl2D%isthere(i,'phshift') )then
                                phshift = stream_proj%os_ptcl2D%get(i,'phshift')
                                call pool_proj%os_ptcl2D%set(iptcl,'phshift',phshift)
                                call pool_proj%os_ptcl3D%set(iptcl,'phshift',phshift)
                            endif
                        endif
                    enddo
                enddo
                call stream_proj%kill
                call os_stk%kill
                if( debug_here ) print *,'end update_pool'; call flush(6)
            end subroutine update_pool

            subroutine classify_buffer
                type(cluster2D_commander_distr) :: xcluster2D_distr
                if( debug_here ) print *,'in classify_buffer'; call flush(6)
                buffer_id = buffer_id + 1
                write(logfhandle,'(A,I6)')'>>> 2D CLASSIFICATION OF NEW BUFFER: ',buffer_id
                ! directory structure
                call chdir('buffer2D')
                ! cluster2d execution
                params_glob%projfile = trim('./'//trim(PROJFILE_BUFFER))
                call cline_cluster2D_buffer%set('startit',  1.)
                call cline_cluster2D_buffer%set('box',      real(box))
                call cline_cluster2D_buffer%set('ncls',     real(params%ncls_start))
                call cline_cluster2D_buffer%set('nparts',   real(params%nparts))
                call cline_cluster2D_buffer%set('msk',      large_msk)
                if( l_greedy )then
                    ! done
                else
                    ! stochastic setting
                    call cline_cluster2D_buffer%set('extr_iter', real(MAX_EXTRLIM2D-2))
                endif
                call cline_cluster2D_buffer%delete('trs')
                call cline_cluster2D_buffer%delete('endit')
                call cline_cluster2D_buffer%delete('converged')
                params_glob%nptcls = buffer_proj%get_nptcls()
                call buffer_proj%kill
                call xcluster2D_distr%execute(cline_cluster2D_buffer)
                call getcwd(buffer_dir)
                call rename(PROJFILE_BUFFER, '../'//trim(PROJFILE_BUFFER))
                call chdir('..')
                call simple_getcwd(cwd_glob)
                call buffer_proj%read(PROJFILE_BUFFER)
                params_glob%projfile = trim(orig_projfile)
                params_glob%nptcls   = pool_proj%os_ptcl2D%get_noris()
                if( debug_here ) print *,'end classify_buffer'; call flush(6)
            end subroutine classify_buffer

            !>  runs one iteration of cluster2D for the merged buffers in the cwd
            subroutine classify_pool
                type(oris)                      :: os_cls2D
                character(len=:), allocatable   :: prev_refs
                integer :: nptcls_sel
                nptcls_sel = pool_proj%os_ptcl2D%get_noris() - n_rejected
                write(logfhandle,'(A,I8,A,I4,A)')'>>> 2D CLASSIFICATION OF POOL: ',nptcls_sel,' PARTICLES IN ',ncls_glob,' CLUSTERS'
                ! cluster2d execution
                params_glob%projfile = trim(PROJFILE_POOL)
                params_glob%nptcls   = pool_proj%os_ptcl2D%get_noris()
                call cline_cluster2D%set('projfile',trim(PROJFILE_POOL))
                call cline_cluster2D%set('startit', real(pool_iter))
                call cline_cluster2D%set('maxits',  real(pool_iter))
                call cline_cluster2D%set('ncls',    real(ncls_glob))
                call cline_cluster2D%set('box',     real(box))
                call cline_cluster2D%set('refs',    trim(refs_glob))
                call cline_cluster2D%set('frcs',    trim(FRCS_FILE))
                call cline_cluster2D%delete('endit')
                call cline_cluster2D%delete('converged')
                call exec_cluster2D_stream_distr(cline_cluster2D, pool_proj)
                refs_glob = trim(CAVGS_ITER_FBODY)//trim(int2str_pad(pool_iter,3))//trim(params%ext)
                params_glob%projfile = trim(orig_projfile)
                call pool_proj%write_segment_inside('cls2D', params_glob%projfile) ! for gui
                call pool_proj%write_segment_inside('out', params_glob%projfile)   ! for gui
                ! removes previous references
                if(pool_iter>1)then
                    prev_refs = trim(CAVGS_ITER_FBODY)//trim(int2str_pad(pool_iter-1,3))//trim(params%ext)
                    ! call del_file(prev_refs) ! for now
                    call del_file(add2fbody(trim(prev_refs),params%ext,'_even'))
                    call del_file(add2fbody(trim(prev_refs),params%ext,'_odd'))
                endif
                !! DEBUG
                call pool_proj%os_cls2D%write('classdoc_pool_current.txt')
                !! DEBUG
                if( debug_here ) call pool_proj%os_cls2D%write('classdoc_pool_'//int2str_pad(n_transfers,4)//'.txt')
                call os_cls2D%kill
            end subroutine classify_pool

            !>  append the classified buffer to the pool
            subroutine transfer_buffer_to_pool
                type(projection_frcs)         :: frcs_glob, frcs_buffer, frcs_prev
                character(len=:), allocatable :: refs_buffer
                integer,          allocatable :: cls_pop(:), cls_buffer_pop(:), pinds(:), pool_inds(:), states(:)
                integer :: endit, iptcl, ind, ncls_here, icls, i, stat, poolind, n_remap, pop
                if( debug_here )print *,'in transfer_buffer_to_pool'; call flush(6)
                n_transfers = n_transfers+1
                write(logfhandle,'(A,I4)')'>>> TRANSFER BUFFER PARTICLES CLASSIFICATION TO POOL #',n_transfers
                ! max # of classes reached ?
                l_maxed = ncls_glob >= max_ncls
                ! updates # of classes
                ncls_here = ncls_glob
                if( .not.l_maxed ) ncls_here = ncls_glob+params%ncls_start
                ! buffer to pool & states mapping
                pool_inds = nint(buffer_proj%os_ptcl2D%get_all('poolind'))
                if( any(pool_inds==0) ) THROW_HARD('pool index indexing error! transfer_buffer_to_pool')
                states = nint(buffer_proj%os_ptcl2D%get_all('state'))
                ! class averages & FRCs
                endit       = nint(cline_cluster2D_buffer%get_rarg('endit'))
                refs_buffer = trim(buffer_dir)//'/'//trim(CAVGS_ITER_FBODY)//int2str_pad(endit,3)//params%ext
                call frcs_buffer%new(params%ncls_start, box, smpd, nstates=1)
                call frcs_buffer%read(trim(buffer_dir)//'/'//trim(FRCS_FILE))
                if( l_maxed )then
                    ! transfer all others
                    cls_pop = nint(pool_proj%os_cls2D%get_all('pop'))
                    n_remap = 0
                    if( any(cls_pop==0) )then
                        if( all(cls_pop==0) ) THROW_HARD('Empty os_cls2D!')
                        ! remapping
                        cls_buffer_pop = nint(buffer_proj%os_cls2D%get_all('pop'))
                        call frcs_glob%new(ncls_glob, box, smpd, nstates=1)
                        call frcs_glob%read(FRCS_FILE)
                        do icls=1,ncls_glob
                            if( cls_pop(icls)>0 ) cycle          ! class already filled
                            if( all(cls_buffer_pop == 0 ) ) exit ! no more buffer class available
                            ind = irnd_uni(params%ncls_start)    ! selects buffer class stochastically
                            do while( cls_buffer_pop(ind) == 0 )
                                ind = irnd_uni(params%ncls_start)
                            enddo
                            cls_buffer_pop(ind) = 0             ! excludes from being picked again
                            n_remap = n_remap+1
                            ! class average
                            call transfer_cavg(refs_buffer, ind, refs_glob, icls)
                            ! frcs
                            call frcs_glob%set_frc(icls,frcs_buffer%get_frc(ind, box, 1), 1)
                            ! class parameters transfer
                            call pool_proj%os_cls2D%transfer_ori(icls, buffer_proj%os_cls2D, ind)
                            call pool_proj%os_cls2D%set(icls, 'class', real(icls))
                            ! particles
                            call buffer_proj%os_ptcl2D%get_pinds(ind,'class',pinds,consider_w=.false.)
                            pop = size(pinds)
                            do i=1,pop
                                iptcl   = pinds(i)
                                poolind = pool_inds(iptcl)
                                call pool_proj%os_ptcl2D%transfer_2Dparams(poolind, buffer_proj%os_ptcl2D, iptcl)
                                call pool_proj%os_ptcl2D%set(poolind,'class',real(icls))
                            enddo
                            cls_pop(icls) = cls_pop(icls) + pop ! updates class populations
                        enddo
                        ! now transfer particles that were not remapped
                        do icls = 1,params%ncls_start
                            if( cls_buffer_pop(icls) == 0 ) cycle
                            ! particles 2D
                            call buffer_proj%os_ptcl2D%get_pinds(icls,'class',pinds,consider_w=.false.)
                            pop = size(pinds)
                            do i = 1,pop
                                iptcl   = pinds(i)
                                poolind = pool_inds(iptcl)
                                ind     = irnd_uni(ncls_glob)  ! stochastic labelling
                                do while( cls_pop(ind) == 0 )  ! but followed by greedy search
                                    ind = irnd_uni(ncls_glob)
                                enddo
                                call pool_proj%os_ptcl2D%transfer_2Dparams(poolind,buffer_proj%os_ptcl2D,iptcl)
                                call pool_proj%os_ptcl2D%set(poolind,'class',real(ind))
                                cls_pop(ind) = cls_pop(ind) + 1 ! updates class populations
                            enddo
                        enddo
                        if( n_remap > 0 )then
                            ! finalize
                            call transfer_cavg(refs_glob,ncls_glob,refs_glob,ncls_glob, self_transfer=.true.)
                            call frcs_glob%write(FRCS_FILE)
                            write(logfhandle,'(A,I4)')'>>> # OF RE-MAPPED CLASS AVERAGES: ',n_remap
                        endif
                    else
                        ! no remapping, just transfer particles & updates 2D population
                        do iptcl=1,buffer_proj%os_ptcl2D%get_noris()
                            if( states(iptcl) /= 0 )then
                                icls          = irnd_uni(ncls_glob) ! stochastic labelling, but followed by greedy search
                                cls_pop(icls) = cls_pop(icls) + 1   ! updates populations
                                poolind       = pool_inds(iptcl)
                                call pool_proj%os_ptcl2D%transfer_2Dparams(poolind, buffer_proj%os_ptcl2D, iptcl)
                                call pool_proj%os_ptcl2D%set(poolind, 'class', real(icls))
                            endif
                        enddo
                    endif
                    ! updates class populations
                    call pool_proj%os_cls2D%set_all('pop',real(cls_pop))
                else
                    if( ncls_glob == 0 )then
                        ! first transfer
                        ! class averages
                        refs_glob = 'start_cavgs'//params%ext
                        do icls= 1,params%ncls_start
                            call transfer_cavg(refs_buffer,icls,refs_glob,icls)
                        enddo
                        ! FRCs
                        call simple_copy_file(trim(buffer_dir)//'/'//trim(FRCS_FILE),FRCS_FILE)
                        ! class parameters
                        pool_proj%os_cls2D = buffer_proj%os_cls2D
                    else
                        ! append new classes
                        ! class averages
                        do icls=1,params%ncls_start
                            call transfer_cavg(refs_buffer, icls, refs_glob, ncls_glob+icls)
                        enddo
                        ! FRCs
                        call frcs_glob%new(ncls_here, box, smpd, nstates=1)
                        call frcs_prev%new(ncls_glob, box, smpd, nstates=1)
                        call frcs_prev%read(FRCS_FILE)
                        do icls=1,ncls_glob
                            call frcs_glob%set_frc(icls,frcs_prev%get_frc(icls, box, 1), 1)
                        enddo
                        do icls=1,params%ncls_start
                            call frcs_glob%set_frc(ncls_glob+icls,frcs_buffer%get_frc(icls, box, 1), 1)
                        enddo
                        call frcs_glob%write(FRCS_FILE)
                        ! class parameters
                        call pool_proj%os_cls2D%reallocate(ncls_here)
                        do icls=1,params%ncls_start
                            ind = ncls_glob+icls
                            call pool_proj%os_cls2D%transfer_ori(ind, buffer_proj%os_cls2D, icls)
                            call pool_proj%os_cls2D%set(ind, 'class', real(ind))
                        enddo
                    endif
                    ! particles 2D
                    do iptcl=1,buffer_proj%os_ptcl2D%get_noris()
                        if( states(iptcl) /= 0 )then
                            call buffer_proj%os_ptcl2D%set(iptcl, 'class', real(ncls_glob+nint(buffer_proj%os_ptcl2D%get(iptcl,'class'))))
                            call pool_proj%os_ptcl2D%transfer_2Dparams(pool_inds(iptcl), buffer_proj%os_ptcl2D, iptcl)
                        endif
                    enddo
                    ! global # of classes
                    ncls_glob = ncls_here
                endif
                ! Finally, deals with states and particles 3D, and search strategy
                do iptcl = 1,buffer_proj%os_ptcl2D%get_noris()
                    poolind = pool_inds(iptcl)
                    if( states(iptcl) == 0 )then
                        ! maps de-selected particles
                        call pool_proj%os_ptcl2D%set(poolind, 'state', 0.)
                        call pool_proj%os_ptcl3D%set(poolind, 'state', 0.)
                    else
                        call pool_proj%os_ptcl2D%delete_entry(poolind, 'poolind') ! declutter
                        ! greedy search at next iteration for all newly added particles!
                        call pool_proj%os_ptcl2D%set(poolind, 'updatecnt', 0.)
                        ! 3D eo flags
                        call pool_proj%os_ptcl3D%set(poolind, 'eo', pool_proj%os_ptcl2D%get(poolind,'eo'))
                    endif
                enddo
                ! preserve buffer for analysis. To delete
                stat = rename(refs_buffer,'cavgs_buffer'//int2str_pad(n_transfers,4)//params%ext)
                ! cleanup
                call frcs_prev%kill
                call frcs_glob%kill
                call frcs_buffer%kill
                if( debug_here )print *,'end transfer_buffer_to_pool'; call flush(6)
            end subroutine transfer_buffer_to_pool

            !>  Convenience function for transfer_buffer_to_pool
            subroutine transfer_cavg( refs_in, indin, refs_out, indout, self_transfer )
                character(len=*),  intent(in) :: refs_in, refs_out
                integer,           intent(in) :: indin, indout
                logical, optional, intent(in) :: self_transfer
                type(image)                   :: img
                character(len=:), allocatable :: stkout, stkin
                integer :: ipart
                logical :: l_self
                l_self = .false.
                if( present(self_transfer) ) l_self = self_transfer
                call img%new([box,box,1],smpd)
                call img%read( refs_in, indin)
                call img%write(refs_out,indout)
                stkin  = add2fbody(refs_in, params%ext,'_even')
                stkout = add2fbody(refs_out,params%ext,'_even')
                call img%read( stkin, indin)
                call img%write(stkout,indout)
                stkin  = add2fbody(refs_in,params%ext,'_odd')
                stkout = add2fbody(refs_out,params%ext,'_odd')
                call img%read( stkin, indin)
                call img%write(stkout,indout)
                ! temporary matrices
                call img%new([boxpd,boxpd,1],smpd)
                call img%zero_and_flag_ft
                do ipart = 1,params%nparts
                    stkin  = trim(buffer_dir)//'/cavgs_even_part'//int2str_pad(ipart,params%numlen)//params%ext
                    stkout = 'cavgs_even_part'//int2str_pad(ipart,params%numlen)//params%ext
                    if( l_self ) stkin = stkout
                    call img%read(stkin, indin)
                    call img%write(stkout,indout)
                    stkin  = trim(buffer_dir)//'/cavgs_odd_part'//int2str_pad(ipart,params%numlen)//params%ext
                    stkout = 'cavgs_odd_part'//int2str_pad(ipart,params%numlen)//params%ext
                    if( l_self ) stkin = stkout
                    call img%read(stkin, indin)
                    call img%write(stkout,indout)
                    stkin  = trim(buffer_dir)//'/ctfsqsums_even_part'//int2str_pad(ipart,params%numlen)//params%ext
                    stkout = 'ctfsqsums_even_part'//int2str_pad(ipart,params%numlen)//params%ext
                    if( l_self ) stkin = stkout
                    call img%read(stkin, indin)
                    call img%write(stkout,indout)
                    stkin  = trim(buffer_dir)//'/ctfsqsums_odd_part'//int2str_pad(ipart,params%numlen)//params%ext
                    stkout = 'ctfsqsums_odd_part'//int2str_pad(ipart,params%numlen)//params%ext
                    if( l_self ) stkin = stkout
                    call img%read(stkin, indin)
                    call img%write(stkout,indout)
                enddo
                ! cleanup
                call img%kill
            end subroutine transfer_cavg

            !> returns the list of projects necessary for creating a new buffer
            subroutine read_mics( new_buffer_spprojs, n_buffer_spprojs, n_buffer_stks, n_buffer_ptcls )
                character(len=LONGSTRLEN), allocatable, intent(inout) :: new_buffer_spprojs(:)
                integer,                                intent(inout) :: n_buffer_spprojs, n_buffer_stks, n_buffer_ptcls
                type(sp_project)                       :: stream_proj
                character(len=LONGSTRLEN), allocatable :: tmp(:)
                logical,                   allocatable :: spproj_mask(:)
                integer :: nptcls, iproj, jproj, cnt, nmics_imported, n, cnt2
                logical :: isnew
                if( debug_here ) print *,'in read_mics'; call flush(6)
                call read_filetable(spproj_list_fname, spproj_list)
                if(allocated(new_buffer_spprojs)) deallocate(new_buffer_spprojs)
                n_buffer_spprojs = 0
                n_buffer_stks    = 0
                n_buffer_ptcls   = 0
                nmics_imported   = merge(size(imported_spprojs),0,allocated(imported_spprojs))
                if( allocated(spproj_list) )then
                    n_spprojs = size(spproj_list)
                    allocate(spproj_mask(n_spprojs),source=.false.)
                    nptcls = 0
                    cnt    = 0
                    cnt2   = 0
                    do iproj = 1,n_spprojs
                        ! identifies whether this is a new project
                        isnew = .true.
                        do jproj = 1,nmics_imported
                            if( trim(spproj_list(iproj)) .eq. trim(imported_spprojs(jproj)) )then
                                isnew = .false.
                                exit
                            endif
                        enddo
                        spproj_mask(iproj) = isnew
                        if( .not. spproj_mask(iproj) ) cycle
                        ! enough for a new buffer?
                        call stream_proj%read_segment('mic',spproj_list(iproj))
                        n = nint(stream_proj%os_mic%get(1,'nptcls'))
                        cnt    = cnt + 1          ! micrographs
                        nptcls = nptcls + n       ! particles
                        if( n > 0 ) cnt2 = cnt2+1 ! stacks
                        if( nptcls > nptcls_per_buffer )then
                            n_buffer_spprojs = cnt
                            n_buffer_stks    = cnt2
                            n_buffer_ptcls   = nptcls
                            exit
                        endif
                    enddo
                    call stream_proj%kill
                    if( n_buffer_spprojs > 0 )then
                        ! new projects list for buffer
                        allocate(new_buffer_spprojs(n_buffer_spprojs))
                        cnt = 0
                        do iproj = 1,n_spprojs
                            if( spproj_mask(iproj) )then
                                cnt = cnt + 1
                                new_buffer_spprojs(cnt) = trim(spproj_list(iproj))
                                if( debug_here )print *,'new mics ',cnt,trim(spproj_list(iproj)); call flush(6)
                            endif
                        enddo
                        ! updates microgrpahs imported list
                        if( nmics_imported == 0 )then
                            imported_spprojs = new_buffer_spprojs(:)
                        else
                            tmp = imported_spprojs(:)
                            deallocate(imported_spprojs)
                            allocate(imported_spprojs(nmics_imported+n_buffer_spprojs))
                            imported_spprojs(1:nmics_imported) = tmp(:)
                            imported_spprojs(nmics_imported+1:nmics_imported+n_buffer_spprojs) = new_buffer_spprojs(:)
                        endif
                    endif
                endif
                call stream_proj%kill
                if( debug_here )then
                    print *,'n_buffer_spprojs ',n_buffer_spprojs
                    print *,'n_buffer_stks    ',n_buffer_stks
                    print *,'n_buffer_ptcls   ',n_buffer_ptcls
                    print *,'end read_mics'; call flush(6)
                endif
            end subroutine read_mics

            !> builds and writes new buffer project from the pool
            subroutine gen_buffer_from_pool
                use simple_ori, only: ori
                type(ori)                     :: ostk, optcl
                character(len=:), allocatable :: stkname
                integer :: i, iptcl, fromp_pool, fromp, top, istk, cnt, nptcls, fromstk, tostk
                if( debug_here )print *,'in gen_buffer_from_pool'; call flush(6)
                buffer_exists = .false.
                if( n_buffer_spprojs == 0 )return
                call buffer_proj%kill
                call simple_mkdir('buffer2D')
                call simple_chdir('buffer2D')
                buffer_proj%projinfo = pool_proj%projinfo
                buffer_proj%compenv  = pool_proj%compenv
                if( pool_proj%jobproc%get_noris()>0 ) buffer_proj%jobproc = pool_proj%jobproc
                call buffer_proj%projinfo%delete_entry('projname')
                call buffer_proj%projinfo%delete_entry('projfile')
                call buffer_proj%update_projinfo(cline_cluster2D_buffer)
                call buffer_proj%os_stk%new(n_buffer_stks)
                call buffer_proj%os_ptcl2D%new(n_buffer_ptcls)
                call buffer_proj%os_ptcl3D%new(n_buffer_ptcls)
                fromstk = max(1,nstks_imported-n_buffer_stks+1)
                tostk   = nstks_imported
                cnt   = 0
                fromp = 1
                do istk = fromstk,tostk
                    cnt = cnt+1
                    call pool_proj%os_stk%get_ori(istk, ostk)
                    ! transfer particles
                    nptcls     = nint(ostk%get('nptcls'))
                    fromp_pool = nint(ostk%get('fromp'))
                    top = fromp + nptcls-1
                    i   = fromp_pool-1
                    do iptcl=fromp,top
                        i = i+1
                        call pool_proj%os_ptcl2D%get_ori(i, optcl)
                        call optcl%set('stkind',real(cnt))
                        call buffer_proj%os_ptcl3D%set_ori(iptcl, optcl)
                        call optcl%set('poolind', real(i)) ! Original index in pool, critical!
                        call buffer_proj%os_ptcl2D%set_ori(iptcl, optcl)
                    enddo
                    ! stack transfer
                    call ostk%getter('stk',stkname)
                    call ostk%set('stk',   '../'//trim(stkname))
                    call ostk%set('fromp', real(fromp))
                    call ostk%set('top',   real(top))
                    call buffer_proj%os_stk%set_ori(cnt, ostk)
                    fromp = top+1
                enddo
                if( debug_here )then
                    if( top /= n_buffer_ptcls )then
                        print *,'top            ',iptcl
                        print *,'n_buffer_ptcls ',n_buffer_ptcls
                        THROW_HARD('Inconsistent # of particles; gen_buffer_from_pool')
                    endif
                endif
                call buffer_proj%write(PROJFILE_BUFFER)
                call simple_chdir('..')
                call simple_getcwd(cwd_glob)
                buffer_exists = .true.
                write(logfhandle,'(A,I4,A,I6,A)')'>>> BUILT NEW BUFFER WITH ',n_buffer_spprojs,' MICROGRAPHS, ',n_buffer_ptcls,' PARTICLES'
                ! cleanup
                call ostk%kill
                call optcl%kill
            end subroutine gen_buffer_from_pool

            logical function is_timeout( time_now )
                integer, intent(in) :: time_now
                is_timeout = .false.
                if( nint(real(time_now-last_injection)/60.) > params%time_inactive)then
                    write(logfhandle,'(A,A)')'>>> TIME LIMIT WITHOUT NEW IMAGES REACHED: ',cast_time_char(time_now)
                    is_timeout = .true.
                    if( .not.cline_cluster2D%defined('converged') )is_timeout = .false.
                    if( .not.is_timeout) write(logfhandle,'(A,A)')'>>> WAITING FOR CONVERGENCE'
                else if(time_now-last_injection > 3600)then
                    write(logfhandle,'(A,A)')'>>> OVER ONE HOUR WITHOUT NEW PARTICLES: ',cast_time_char(time_now)
                    call flush(6)
                endif
            end function is_timeout

            !> produces consolidated project at original scale
            subroutine write_snapshot( force, add_suffix )
                logical,           intent(in) :: force, add_suffix
                type(projection_frcs)         :: frcs, frcs_sc
                type(oris)                    :: os_backup1, os_backup2, os_backup3
                character(len=:), allocatable :: projfname, cavgsfname, suffix, frcsfname, src, dest
                integer :: istk
                if( debug_here ) print *,'in write_snapshot'; call flush(6)
                if( .not.force )then
                    if( .not.file_exists(SPPROJ_SNAPSHOT) )return
                    call del_file(SPPROJ_SNAPSHOT)
                endif
                projfname  = get_fbody(orig_projfile, METADATA_EXT, separator=.false.)
                cavgsfname = get_fbody(refs_glob, params%ext, separator=.false.)
                frcsfname  = get_fbody(FRCS_FILE, BIN_EXT, separator=.false.)
                if( add_suffix )then
                    suffix     = '_'//trim(int2str(pool_proj%os_ptcl2D%get_noris()))//'nptcls'
                    cavgsfname = trim(cavgsfname)//trim(suffix)
                    frcsfname  = trim(frcsfname)//trim(suffix)
                endif
                call pool_proj%projinfo%set(1,'projname', projfname)
                projfname  = trim(projfname)//trim(METADATA_EXT)
                cavgsfname = trim(cavgsfname)//trim(params%ext)
                frcsfname  = trim(frcsfname)//trim(BIN_EXT)
                call pool_proj%projinfo%set(1,'projfile', projfname)
                if( add_suffix )then
                    if( trim(prev_snapshot_cavgs) /= '' )then
                        call del_file(prev_snapshot_frcs)
                        call del_file(prev_snapshot_cavgs)
                        src = add2fbody(prev_snapshot_cavgs, params%ext,'_even')
                        call del_file(src)
                        src = add2fbody(prev_snapshot_cavgs, params%ext,'_odd')
                        call del_file(src)
                    endif
                endif
                write(logfhandle,'(A,A,A,A)')'>>> PRODUCING PROJECT SNAPSHOT ',trim(projfname), ' AT: ',cast_time_char(simple_gettime())
                os_backup3 = pool_proj%os_cls2D
                if( do_autoscale )then
                    os_backup1 = pool_proj%os_ptcl2D
                    os_backup2 = pool_proj%os_stk
                    ! rescale classes
                    call rescale_cavgs(refs_glob, cavgsfname)
                    src  = add2fbody(refs_glob, params%ext,'_even')
                    dest = add2fbody(cavgsfname,params%ext,'_even')
                    call rescale_cavgs(src, dest)
                    src  = add2fbody(refs_glob, params%ext,'_odd')
                    dest = add2fbody(cavgsfname,params%ext,'_odd')
                    call rescale_cavgs(src, dest)
                    call pool_proj%os_out%kill
                    call pool_proj%add_cavgs2os_out(cavgsfname, orig_smpd, 'cavg')
                    pool_proj%os_cls2D = os_backup3
                    ! rescale frcs
                    call frcs_sc%read(FRCS_FILE)
                    call frcs_sc%upsample(orig_smpd, orig_box, frcs)
                    call frcs%write(frcsfname)
                    call frcs%kill
                    call frcs_sc%kill
                    call pool_proj%add_frcs2os_out(frcsfname, 'frc2D')
                    ! project updates to original scale
                    call pool_proj%os_stk%set_all2single('box', real(orig_box))
                    call pool_proj%os_stk%set_all2single('smpd',orig_smpd)
                    call pool_proj%os_ptcl2D%mul_shifts( 1./scale_factor )
                    do istk = 1,pool_proj%os_stk%get_noris()
                        call pool_proj%os_stk%set(istk,'stk',imported_stks(istk))
                    enddo
                else
                    if( add_suffix )then
                        call simple_copy_file(FRCS_FILE, frcsfname)
                        call simple_copy_file(refs_glob, cavgsfname)
                        src  = add2fbody(refs_glob, params%ext,'_even')
                        dest = add2fbody(cavgsfname,params%ext,'_even')
                        call simple_copy_file(src, dest)
                        src  = add2fbody(refs_glob, params%ext,'_odd')
                        dest = add2fbody(cavgsfname,params%ext,'_odd')
                        call simple_copy_file(src, dest)
                    endif
                    call pool_proj%os_out%kill
                    call pool_proj%add_cavgs2os_out(cavgsfname, orig_smpd, 'cavg')
                    pool_proj%os_cls2D = os_backup3
                    call pool_proj%add_frcs2os_out(frcsfname, 'frc2D')
                endif
                call os_backup3%kill
                call pool_proj%write(projfname)
                if( do_autoscale )then
                    ! preserve down-scaling
                    pool_proj%os_ptcl2D = os_backup1
                    call os_backup1%kill
                    pool_proj%os_stk    = os_backup2
                    call os_backup2%kill
                endif
                ! cleanup previous snapshot
                prev_snapshot_frcs  = trim(frcsfname)
                prev_snapshot_cavgs = trim(cavgsfname)
                if( debug_here ) print *,'end write_snapshot'; call flush(6)
            end subroutine write_snapshot

            !> scales stacks and returns updated list of names
            subroutine scale_stks( stk_fnames )
                use simple_commander_project,  only: scale_project_commander_distr
                character(len=*), allocatable, intent(inout) :: stk_fnames(:)
                character(len=:), allocatable       :: fname
                type(scale_project_commander_distr) :: xscale_distr
                type(sp_project) :: dummy_proj
                type(oris)       :: os
                type(cmdline)    :: cline_scale
                integer          :: istk, n
                logical          :: err
                if( .not.do_autoscale )return
                if( .not.allocated(stk_fnames) )return
                err = .false.
                call simple_mkdir(SCALE_DIR, errmsg= "commander_stream_wflows:: cluster2D_stream scale_stks")
                n = size(stk_fnames)
                ! command-line
                call cline_scale%set('prg',        'scale_project')
                call cline_scale%set('projfile',   'forscale.simple')
                call cline_scale%set('projname',   'forscale')
                call cline_scale%set('smpd',       orig_smpd)
                call cline_scale%set('box',        real(orig_box))
                call cline_scale%set('newbox',     real(box))
                call cline_scale%set('nparts',     real(params%nparts))
                call cline_scale%set('nthr',       real(params%nthr))
                call cline_scale%set('mkdir',      'no')
                call cline_scale%set('dir_target', trim(SCALE_DIR))
                ! dummy project & import
                dummy_proj%projinfo = pool_proj%projinfo
                dummy_proj%compenv  = pool_proj%compenv
                call dummy_proj%projinfo%delete_entry('projname')
                call dummy_proj%projinfo%delete_entry('projfile')
                call dummy_proj%update_projinfo(cline_scale)
                call os%new(n)
                call os%set_all2single('state', 1.)
                call os%set_all2single('ctf',   'no')
                call os%set_all2single('smpd',  orig_smpd)
                call os%set_all2single('kv',    300.)
                call os%set_all2single('cs',    2.7)
                call os%set_all2single('fraca', 0.1)
                call os%set_all2single('dfx',   1.)
                call os%set_all2single('dfy',   1.)
                call os%set_all2single('angast',0.)
                call os%set_all2single('phaseplate','no')
                call dummy_proj%add_stktab(stk_fnames, os)
                call dummy_proj%write('forscale.simple')
                ! execution
                call xscale_distr%execute(cline_scale)
                do istk = 1,dummy_proj%os_stk%get_noris()
                    fname = add2fbody(stk_fnames(istk), params%ext, SCALE_SUFFIX)
                    stk_fnames(istk) = filepath(trim(SCALE_DIR), basename(fname))
                    if( .not.file_exists(stk_fnames(istk)))then
                        write(logfhandle,*)'stack does not exists: ',n,istk,trim(stk_fnames(istk))
                        err = .true.
                    endif
                enddo
                call os%kill
                call dummy_proj%kill
                if( err )then
                    THROW_HARD('CWD: '//trim(CWD_GLOB))
                else
                    call qsys_cleanup
                    call del_file('forscale.simple')
                endif
            end subroutine scale_stks

            !> for initial write of set of user adjustable parameters
            subroutine write_user_params
                type(oris) :: os
                call os%new(1)
                call os%set(1,'lpthresh',params%lpthresh)
                call os%set(1,'ndev',    params%ndev)
                call os%write(USER_PARAMS)
                call os%kill
            end subroutine write_user_params

            !> updates current parameters with user input
            subroutine update_user_params
                type(oris) :: os
                real       :: lpthresh, ndev
                if( .not.file_exists(USER_PARAMS) ) return ! use of default/last update
                lpthresh = params%lpthresh
                ndev     = params%ndev
                call os%new(1)
                call os%read(USER_PARAMS)
                if( os%isthere(1,'lpthresh') )then
                    lpthresh = os%get(1,'lpthresh')
                    if( abs(lpthresh-params%lpthresh) > 0.001 )then
                        params%lpthresh = lpthresh
                        write(logfhandle,'(A,F8.2)')'>>> REJECTION LPTHESH UPDATED TO: ',params%lpthresh
                    endif
                endif
                if( os%isthere(1,'ndev') )then
                    ndev = os%get(1,'ndev')
                    if( abs(ndev-params%ndev) > 0.001 )then
                        params%ndev = ndev
                        write(logfhandle,'(A,F8.2)')'>>> REJECTION NDEV    UPDATED TO: ',params%ndev
                    endif
                endif
                call os%kill
            end subroutine update_user_params

    end subroutine exec_cluster2D_stream

    subroutine exec_cluster2D_autoscale( self, cline )
        use simple_commander_project, only: scale_project_commander_distr, prune_project_commander_distr
        use simple_commander_imgproc, only: scale_commander, pspec_int_rank_commander
        class(cluster2D_autoscale_commander_hlev), intent(inout) :: self
        class(cmdline),                            intent(inout) :: cline
        ! constants
        integer,               parameter :: MAXITS_STAGE1      = 10
        integer,               parameter :: MAXITS_STAGE1_EXTR = 15
        character(len=STDLEN), parameter :: orig_projfile_bak  = 'orig_bak.simple'
        ! commanders
        type(prune_project_commander_distr) :: xprune_project
        type(make_cavgs_commander_distr)    :: xmake_cavgs
        type(cluster2D_commander_distr)     :: xcluster2D_distr
        type(pspec_int_rank_commander)      :: xpspec_rank
        type(rank_cavgs_commander)          :: xrank_cavgs
        type(scale_commander)               :: xscale
        type(scale_project_commander_distr) :: xscale_distr
        ! command lines
        type(cmdline) :: cline_cluster2D_stage1
        type(cmdline) :: cline_cluster2D_stage2
        type(cmdline) :: cline_scalerefs, cline_scale1, cline_scale2
        type(cmdline) :: cline_make_cavgs, cline_rank_cavgs, cline_pspec_rank
        type(cmdline) :: cline_prune_project
        ! other variables
        type(parameters)              :: params
        type(sp_project)              :: spproj, spproj_sc
        character(len=:), allocatable :: projfile_sc, orig_projfile
        character(len=LONGSTRLEN)     :: finalcavgs, finalcavgs_ranked, refs_sc
        real     :: scale_stage1, scale_stage2, trs_stage2
        integer  :: nparts, last_iter_stage1, last_iter_stage2, status
        integer  :: nptcls_sel
        logical  :: scaling
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('oritype')   ) call cline%set('oritype', 'ptcl2D')
        if( .not. cline%defined('lpstart')   ) call cline%set('lpstart',     15. )
        if( .not. cline%defined('lpstop')    ) call cline%set('lpstop',       8. )
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',       30. )
        if( .not. cline%defined('maxits')    ) call cline%set('maxits',      30. )
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale',  'yes')
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
        ! nptcls_sel = spproj%os_ptcl2D%get_noris(consider_state=.true.)
        ! if( nptcls_sel < nint(PRUNE_FRAC*real(spproj%os_ptcl2D%get_noris())) )then
        !     write(logfhandle,'(A)')'>>> AUTO-PRUNING PROJECT FILE'
        !     call spproj%kill
        !     cline_prune_project = cline
        !     call xprune_project%execute(cline_prune_project)
        !     call spproj%read(params%projfile)
        !     params%nptcls = nptcls_sel
        ! endif
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
        call spproj%write_segment_inside('out', params%projfile)
        ! re-activating CTF correction for time series
        if( trim(params%tseries).eq.'yes' )then
        endif
        ! clean
        call spproj%kill()
        ! ranking
        if( trim(params%tseries).eq.'yes' )then
            ! rank based on maximum of power spectrum
            call cline_pspec_rank%set('mkdir',   'no')
            call cline_pspec_rank%set('moldiam', params%moldiam)
            call cline_pspec_rank%set('nthr',    real(params%nthr))
            call cline_pspec_rank%set('smpd',    params%smpd)
            call cline_pspec_rank%set('stk',     finalcavgs)
            if( cline%defined('lp_backgr') ) call cline_pspec_rank%set('lp_backgr', params%lp_backgr)
            call xpspec_rank%execute(cline_pspec_rank)
        else
            ! rank baseed on gold-standard resolution estimates
            finalcavgs_ranked = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//'_ranked'//params%ext
            call cline_rank_cavgs%set('projfile', trim(params%projfile))
            call cline_rank_cavgs%set('stk',      trim(finalcavgs))
            call cline_rank_cavgs%set('outstk',   trim(finalcavgs_ranked))
            call xrank_cavgs%execute( cline_rank_cavgs )
        endif
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
        integer                   :: iter, cnt, iptcl, ptclind
        type(chash)               :: job_descr
        real                      :: frac_srch_space
        if( .not. cline%defined('lpstart')   ) call cline%set('lpstart',    15. )
        if( .not. cline%defined('lpstop')    ) call cline%set('lpstop',      8. )
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      30. )
        if( .not. cline%defined('maxits')    ) call cline%set('maxits',     30. )
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('oritype')   ) call cline%set('oritype', 'ptcl2D')
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
        ! splitting
        call build%spproj%split_stk(params%nparts)
        ! prepare command lines from prototype master
        cline_check_2Dconv = cline
        cline_cavgassemble = cline
        cline_make_cavgs   = cline ! ncls is transferred here
        ! initialise static command line parameters and static job description parameters
        call cline_cavgassemble%set('prg', 'cavgassemble')
        call cline_make_cavgs%set('prg',   'make_cavgs')
        ! execute initialiser
        if( .not. cline%defined('refs') )then
            refs             = 'start2Drefs' // params%ext
            params%refs      = trim(refs)
            params%refs_even = 'start2Drefs_even'//params%ext
            params%refs_odd  = 'start2Drefs_odd'//params%ext
            if( build%spproj%is_virgin_field('ptcl2D') .or. params%startit == 1 )then
                if( params%tseries .eq. 'yes' )then
                    if( cline%defined('nptcls_per_cls') )then
                        if( build%spproj%os_ptcl2D%any_state_zero() )then
                            THROW_HARD('cluster2D_nano does not allow state=0 particles, prune project before execution; exec_cluster2D_distr')
                        endif
                        cnt = 0
                        do iptcl=1,params%nptcls,params%nptcls_per_cls
                            cnt = cnt + 1
                            params%ncls = cnt
                            do ptclind=iptcl,min(params%nptcls, iptcl + params%nptcls_per_cls - 1)
                                call build%spproj%os_ptcl2D%set(ptclind, 'class', real(cnt))
                            end do
                        end do
                        call job_descr%set('ncls',int2str(params%ncls))
                        call cline%set('ncls', real(params%ncls))
                        call cline_make_cavgs%set('ncls', real(params%ncls))
                        call cline_cavgassemble%set('ncls', real(params%ncls))
                        call cline_make_cavgs%set('refs', params%refs)
                        call xmake_cavgs%execute(cline_make_cavgs)
                    else
                        call selection_from_tseries_imgfile(build%spproj, params%refs, params%box, params%ncls)
                    endif
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
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype,params%projfile)
        endif
        ! main loop
        iter = params%startit - 1
        do
            iter = iter + 1
            params_glob%which_iter = iter
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
            ! merge orientation documents and transfers to memory
            call build%spproj%merge_algndocs(params%nptcls, params%nparts, 'ptcl2D', ALGN_FBODY)
            ! assemble class averages
            refs      = trim(CAVGS_ITER_FBODY) // trim(str_iter)            // params%ext
            refs_even = trim(CAVGS_ITER_FBODY) // trim(str_iter) // '_even' // params%ext
            refs_odd  = trim(CAVGS_ITER_FBODY) // trim(str_iter) // '_odd'  // params%ext
            call cline_cavgassemble%set('refs', trim(refs))
            call qenv%exec_simple_prg_in_queue(cline_cavgassemble, 'CAVGASSEMBLE_FINISHED')
            ! check convergence
            call check_2Dconv(cline_check_2Dconv, build%spproj_field)
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
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
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
        type(parameters)  :: params
        type(builder)     :: build
        real, allocatable :: states(:)
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
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
        ! write project: cls2D and state congruent cls3D
        call build%spproj%os_cls3D%new(params%ncls)
        states = build%spproj%os_cls2D%get_all('state')
        call build%spproj%os_cls3D%set_all('state',states)
        call build%spproj%write_segment_inside('cls2D', params%projfile)
        call build%spproj%write_segment_inside('cls3D', params%projfile)
        ! end gracefully
        call simple_end('**** SIMPLE_CAVGASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        call simple_touch('CAVGASSEMBLE_FINISHED', errmsg='In: commander_rec :: eo_cavgassemble ')
    end subroutine exec_cavgassemble

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
        use simple_polarizer, only: polarizer
        use simple_cluster_cavgs
        class(cluster_cavgs_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)             :: params
        type(sp_project)             :: spproj
        type(polarizer), allocatable :: cavg_imgs(:)
        character(len=:),allocatable :: cavgsstk, classname
        integer,         allocatable :: centers(:), labels(:), cntarr(:)
        real,            allocatable :: states(:), rtmparr(:), orig_cls_inds(:)
        integer :: ncls, n, ldim(3), ncls_sel, i, icls, ncls_aff_prop
        real    :: smpd
        ! defaults
        call cline%set('match_filt', 'no')
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir', 'yes')
        call cline%set('oritype', 'ptcl2D')
        call params%new(cline)
        call spproj%read(params%projfile)
        ! need to turn off CTF since it may be on based on the imported images and now we operate on class averages
        params%ctf = 'no'
        ! get class average stack
        call spproj%get_cavgs_stk(cavgsstk, ncls, smpd)
        call find_ldim_nptcls(cavgsstk, ldim, n)
        ldim(3) = 1
        if( n /= ncls ) THROW_HARD('Incosistent # classes in project file vs cavgs stack; exec_cluster_cavgs')
        ! ensure correct smpd in params class
        params%smpd = smpd
        ! get state flag array
        states = spproj%os_cls2D%get_all('state')
        ! get ncls from ptcl2D field
        n = spproj%os_ptcl2D%get_n('class')
        if( n /= ncls ) THROW_HARD('Incosistent # classes in ptcl2D field of spproj vs cavgs stack; exec_cluster_cavgs')
        ! find out how many selected class averages
        ncls_sel = count(states > 0.5)
        ! make class index array for the original classes excluding state=0 ones
        rtmparr = spproj%os_cls2D%get_all('class')
        if( allocated(rtmparr) )then
            orig_cls_inds = pack(rtmparr, mask=states > 0.5)
            deallocate(rtmparr)
        endif
        ! allocate polarizer images and read
        allocate(cavg_imgs(ncls_sel))
        do i=1,ncls_sel
            call cavg_imgs(i)%new(ldim, smpd, wthreads=.false.)
            icls = nint(orig_cls_inds(i))
            call cavg_imgs(i)%read(cavgsstk, icls)
        end do
        ! rotationally invariant clustering of class averages with affinity propagation
        call cluster_cavgs(cavg_imgs, centers, labels)
        ncls_aff_prop = size(centers)
        write(logfhandle,'(A,I3)') '>>> # CLUSTERS FOUND BY AFFINITY PROPAGATION: ', ncls_aff_prop
        allocate(cntarr(ncls_aff_prop), source=0)
        ! read back the original (unprocessed) images
        do i=1,ncls_sel
            call cavg_imgs(i)%new(ldim, smpd, wthreads=.false.)
            icls = nint(orig_cls_inds(i))
            call cavg_imgs(i)%read(cavgsstk, icls)
        end do
        ! write the classes
        do icls=1,ncls_aff_prop
            ! make a filename for the class
            do i=1,ncls_sel
                if( labels(i) == icls )then
                    classname = 'class'//int2str_pad(icls,5)//'.mrcs'
                    cntarr(labels(i)) = cntarr(labels(i)) + 1
                    call cavg_imgs(i)%write(classname, cntarr(labels(i)))
                endif
            end do
        end do
        ! destruct
        call spproj%kill
        do icls=1,ncls_sel
            call cavg_imgs(icls)%kill
        end do
        deallocate(cavg_imgs)
        if( allocated(centers)       ) deallocate(centers)
        if( allocated(labels)        ) deallocate(labels)
        if( allocated(states)        ) deallocate(states)
        if( allocated(orig_cls_inds) ) deallocate(orig_cls_inds)
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER_CAVGS NORMAL STOP ****')
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
            call imgs_class(i)%new(ldim, smpd, wthreads=.false.)
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
                call imgs_class(i)%fft
                call imgs_class(i)%shift2Dserial([-inpls(i,2),-inpls(i,3)])
                call imgs_class(i)%ifft
                call imgs_class(i)%rtsq_serial(inpls(i,1), 0., 0., rmat_rot)
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
        call simple_end('**** SIMPLE_WRITE_CLASSES NORMAL STOP ****')
    end subroutine exec_write_classes

    ! UTILITIES

    subroutine check_2Dconv( cline, os )
        use simple_convergence, only: convergence
        use simple_parameters,  only: params_glob
        use simple_oris,        only: oris
        class(cmdline), intent(inout) :: cline
        class(oris),    intent(inout) :: os
        type(parameters)  :: params
        type(convergence) :: conv
        logical :: converged
        call cline%set('oritype', 'ptcl2D')
        call params%new(cline)
        ! convergence check
        converged = conv%check_conv2D(cline, os, os%get_n('class'), params%msk)
        call cline%set('frac_srch', conv%get('frac_srch'))
        ! activates shift search
        if( params_glob%l_doshift ) call cline%set('trs', params_glob%trs)
        if( converged )then
            call cline%set('converged', 'yes')
        else
            call cline%set('converged', 'no')
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK_2DCONV NORMAL STOP ****', print_simple=.false.)
    end subroutine check_2Dconv

end module simple_commander_cluster2D
