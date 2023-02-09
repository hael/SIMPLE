! concrete commander: cluster2D for simultanous 2D alignment and clustering of single-particle images
module simple_commander_cluster2D
include 'simple_lib.f08'
use simple_builder,          only: builder, build_glob
use simple_cmdline,          only: cmdline
use simple_commander_base,   only: commander_base
use simple_parameters,       only: parameters, params_glob
use simple_sp_project,       only: sp_project
use simple_qsys_env,         only: qsys_env
use simple_image,            only: image
use simple_stack_io,         only: stack_io
use simple_starproject,      only: starproject
use simple_euclid_sigma2
use simple_commander_euclid
use simple_qsys_funs
use simple_procimgstk
use simple_progress
use simple_default_clines
implicit none

public :: cleanup2D_commander_hlev
public :: cluster2D_autoscale_commander
public :: cluster2D_commander_distr
public :: cluster2D_commander
public :: make_cavgs_commander_distr
public :: make_cavgs_commander
public :: cavgassemble_commander
public :: rank_cavgs_commander
public :: cluster_cavgs_commander
public :: write_classes_commander
public :: check_2dconv
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: cleanup2D_commander_hlev
  contains
    procedure :: execute      => exec_cleanup2D
end type cleanup2D_commander_hlev

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
        if( .not. cline%defined('ptclw')   ) call cline%set('ptclw',       'no')
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
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR)
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
        logical :: l_shmem
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        ! set shared-memory flag
        if( cline%defined('nparts') )then
            if( nint(cline%get_rarg('nparts')) == 1 )then
                l_shmem = .true.
                call cline%delete('nparts')
            else
                l_shmem = .false.
            endif
        else
            l_shmem = .true.
        endif
        if( l_shmem .and. .not. cline%defined('refs') ) THROW_HARD('need input refs (filename) for shared-memory execution')
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
        if( L_VERBOSE_GLOB ) write(logfhandle,'(a)') '>>> GENERATING CLUSTER CENTERS'
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
        if( params%mul > 1. ) call build%spproj_field%mul_shifts(params%mul)
        ! setup weights
        if( trim(params%ptclw).eq.'yes' )then
            ! weights are set at search time, so nothing to do here.
        else
            call build%spproj_field%calc_hard_weights2D(params%frac, params%ncls)
        endif
        ! even/odd partitioning
        if( build%spproj_field%get_nevenodd() == 0 ) call build%spproj_field%partition_eo
        ! write
        if( l_shmem )then
            call build%spproj%write_segment_inside(params%oritype)
        else
            if( params%part .eq. 1 ) call build%spproj%write_segment_inside(params%oritype)
        endif
        ! create class averager
        call cavger_new
        ! transfer ori data to object
        call cavger_transf_oridat(build%spproj)
        call cavger_read_euclid_sigma2
        ! standard cavg assembly
        call cavger_assemble_sums( .false. )
        if( l_shmem )then
            call cavger_merge_eos_and_norm
            call cavger_calc_and_write_frcs_and_eoavg(params%frcs, params%which_iter)
            ! classdoc gen needs to be after calc of FRCs
            call cavger_gen2Dclassdoc(build%spproj)
            ! write references
            call cavger_write(trim(params%refs),      'merged')
            call cavger_write(trim(params%refs_even), 'even'  )
            call cavger_write(trim(params%refs_odd),  'odd'   )
        else
            ! write partital sums
            call cavger_readwrite_partial_sums('write')
            call qsys_job_finished('simple_commander_cluster2D :: exec_make_cavgs')
        endif
        call cavger_kill
        ! end gracefully
        call build%kill_strategy2D_tbox
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_MAKE_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_make_cavgs

    subroutine exec_cleanup2D( self, cline )
        use simple_commander_project, only: scale_project_commander_distr
        use simple_procimgstk,        only: random_selection_from_imgfile, random_cls_from_imgfile
        use simple_commander_imgproc, only: scale_commander
        use simple_class_frcs,        only: class_frcs
        class(cleanup2D_commander_hlev), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        ! commanders
        type(cluster2D_commander_distr)     :: xcluster2D_distr
        type(cluster2D_commander)           :: xcluster2D ! shared-memory implementation
        type(scale_project_commander_distr) :: xscale_proj
        type(scale_commander)               :: xscale
        type(rank_cavgs_commander)          :: xrank_cavgs
        type(calc_pspec_commander_distr)    :: xcalc_pspec_distr
        ! command lines
        type(cmdline)                       :: cline_cluster2D1, cline_cluster2D2
        type(cmdline)                       :: cline_rank_cavgs, cline_scale, cline_scalerefs
        type(cmdline)                       :: cline_calc_pspec_distr
        ! other variables
        class(parameters), pointer          :: params_ptr => null()
        type(parameters)                    :: params
        type(sp_project)                    :: spproj, spproj_sc
        type(class_frcs)                    :: frcs, frcs_sc
        character(len=:), allocatable       :: projfile, orig_projfile
        character(len=LONGSTRLEN)           :: finalcavgs, finalcavgs_ranked, cavgs, refs_sc
        real                                :: scale_factor, smpd, lp1, lp2
        integer                             :: last_iter, box
        logical                             :: do_scaling, l_shmem, l_euclid
        ! parameters
        character(len=STDLEN) :: orig_projfile_bak = 'orig_bak.simple'
        integer, parameter    :: MINBOX      = 128
        real,    parameter    :: TARGET_LP   = 15.
        real,    parameter    :: MINITS      =  5., MINITS_FAST =  9.
        real,    parameter    :: MAXITS      = 15., MAXITS_FAST = 18.
        real                  :: SMPD_TARGET = 4.
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('lp')        ) call cline%set('lp',           15.)
        if( .not. cline%defined('ncls')      ) call cline%set('ncls',        200.)
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',        20.)
        if( .not. cline%defined('center')    ) call cline%set('center',      'no')
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale',  'yes')
        if( .not. cline%defined('refine')    ) call cline%set('refine',  'greedy')
        if( .not. cline%defined('oritype')   ) call cline%set('oritype', 'ptcl2D')
        if( .not. cline%defined('wiener')    ) call cline%set('wiener',    'full')
        if( .not. cline%defined('ptclw')     ) call cline%set('ptclw',       'no')
        call cline%set('stream', 'no')
        ! set shared-memory flag
        if( cline%defined('nparts') )then
            if( nint(cline%get_rarg('nparts')) == 1 )then
                l_shmem = .true.
                call cline%delete('nparts')
            else
                l_shmem = .false.
            endif
        else
            l_shmem = .true.
        endif
        ! parse parameters
        call params%new(cline)
        orig_projfile = params%projfile
        if( .not. cline%defined('maxits') )then
            params%maxits = nint(MAXITS)
            call cline%set('maxits', MAXITS)
        endif
        if( .not. cline%defined('mskdiam') )then
            call cline%set('mskdiam', (real(params%box) * 0.8) * params%smpd)
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! read project file
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo')
        ! sanity checks
        if( spproj%get_nptcls() == 0 )then
            THROW_HARD('No particles found in project file: '//trim(params%projfile)//'; exec_cleanup2D_autoscale')
        endif
        ! delete any previous solution
        call spproj%os_ptcl2D%delete_2Dclustering
        call spproj%write_segment_inside(params%oritype)
        ! splitting
        call spproj%split_stk(params%nparts, dir=PATH_PARENT)
        ! deal with eo partitioning
        if( spproj%os_ptcl2D%get_nevenodd() == 0 )then
            call spproj%os_ptcl2D%partition_eo
            call spproj%write_segment_inside(params%oritype,params%projfile)
        endif
        ! objfun = euclid
        l_euclid = .false.
        if( cline%defined('objfun') )then
            l_euclid = ( trim(cline%get_carg('objfun')).eq.'euclid' .or. trim(cline%get_carg('objfun')).eq.'prob' .or. trim(cline%get_carg('objfun')).eq.'test' )
            if( l_euclid )then
                cline_calc_pspec_distr  = cline
                call cline_calc_pspec_distr%set( 'prg', 'calc_pspec' )
                call spproj%os_ptcl2D%set_all2single('w', 1.0)
                call spproj%write_segment_inside(params%oritype)
            endif
        endif
        ! first stage
        ! down-scaling for fast execution, greedy optimisation, no match filter
        ! no incremental learning, no centering
        cline_cluster2D1 = cline
        cline_cluster2D2 = cline
        cline_scale      = cline
        call cline_cluster2D1%set('prg', 'cluster2D')
        call cline_cluster2D1%set('maxits',   MINITS)
        call cline_cluster2D1%set('center',     'no')
        call cline_cluster2D1%set('autoscale',  'no')
        call cline_cluster2D1%set('ptclw',      'no')
        call cline_cluster2D1%set('objfun',     'cc')
        call cline_cluster2D1%delete('update_frac')
        ! second stage
        ! down-scaling for fast execution, greedy optimisation, no match filter
        call cline_cluster2D2%set('prg', 'cluster2D')
        call cline_cluster2D2%set('autoscale',  'no')
        call cline_cluster2D2%set('trs',    MINSHIFT)
        if( .not.cline%defined('maxits') )then
            call cline_cluster2D2%set('maxits', MAXITS)
        endif
        if( l_euclid )then
            call cline_cluster2D2%set('objfun', trim(cline%get_carg('objfun')))
            call cline_cluster2D2%set('needs_sigma', 'yes')
        else
            call cline_cluster2D1%set('objfun', 'cc')
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
            if( l_euclid ) call xcalc_pspec_distr%execute( cline_calc_pspec_distr )
        else
            call autoscale(params%box, params%smpd, SMPD_TARGET, box, smpd, scale_factor)
            if( box < MINBOX ) SMPD_TARGET = params%smpd * real(params%box) / real(MINBOX)
            call spproj%scale_projfile(SMPD_TARGET, projfile, cline_cluster2D1, cline_scale, dir=trim(STKPARTSDIR))
            call spproj%kill
            scale_factor = cline_scale%get_rarg('scale')
            smpd         = cline_scale%get_rarg('smpd')
            box          = nint(cline_scale%get_rarg('newbox'))
            ! noise power estimates for objfun = euclid at scaled sampling
            if( l_euclid )then
                call cline_calc_pspec_distr%set('scale', 1./scale_factor**2)
                call xcalc_pspec_distr%execute( cline_calc_pspec_distr )
            endif
            call cline_scale%set('state',1.)
            call cline_scale%delete('smpd') !!
            call simple_mkdir(trim(STKPARTSDIR),errmsg="commander_hlev_wflows :: exec_cluster2D_autoscale;  ")
            call xscale_proj%execute( cline_scale )
            ! rename scaled projfile and stash original project file
            ! such that the scaled project file has the same name as the original and can be followed from the GUI
            call simple_copy_file(orig_projfile, orig_projfile_bak)
            call spproj%read_non_data_segments(projfile)
            call spproj%projinfo%set(1,'projname',get_fbody(orig_projfile,METADATA_EXT,separator=.false.))
            call spproj%projinfo%set(1,'projfile',orig_projfile)
            call spproj%write_non_data_segments(projfile)
            call spproj%kill
            call simple_rename(projfile,orig_projfile)
            projfile = trim(orig_projfile)
        endif
        lp1 = max(2.*smpd, max(params%lp,TARGET_LP))
        lp2 = max(2.*smpd, params%lp)
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
        call cline_cluster2D1%set('refs', params%refs)
        call cline_cluster2D1%set('lp',   lp1)
        call cline_cluster2D2%set('lp',   lp2)
        ! initialise progress monitor
        call progressfile_init()
        ! execution 1
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A,F6.1)') '>>> STAGE 1, LOW-PASS LIMIT: ',lp1
        write(logfhandle,'(A)') '>>>'
        call cline_cluster2D1%set('projfile', trim(projfile))
        if( l_shmem )then
            params_ptr  => params_glob
            params_glob => null()
            call xcluster2D%execute(cline_cluster2D1)
            params_glob => params_ptr
            params_ptr  => null()
        else
            call xcluster2D_distr%execute(cline_cluster2D1)
        endif
        last_iter  = nint(cline_cluster2D1%get_rarg('endit'))
        finalcavgs = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter,3)//params%ext
        if( .not. file_exists(trim(finalcavgs)) ) THROW_HARD('File '//trim(finalcavgs)//' does not exist')
        ! execution 2
        if( cline%defined('maxits') )then
            if( last_iter < params%maxits )then
                write(logfhandle,'(A)') '>>>'
                write(logfhandle,'(A,F6.1)') '>>> STAGE 2, LOW-PASS LIMIT: ',lp2
                write(logfhandle,'(A)') '>>>'
                call cline_cluster2D2%set('projfile', trim(projfile))
                call cline_cluster2D2%set('startit',  real(last_iter+1))
                call cline_cluster2D2%set('refs',     trim(finalcavgs))
                if( l_shmem )then
                    params_ptr  => params_glob
                    params_glob => null()
                    call xcluster2D%execute(cline_cluster2D2)
                    params_glob => params_ptr
                    params_ptr  => null()
                else
                    call xcluster2D_distr%execute(cline_cluster2D2)
                endif
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
                type(image)    :: img, img_pad
                type(stack_io) :: stkio_w
                integer        :: icls
                call img%new([box,box,1],smpd)
                call img_pad%new([params%box,params%box,1],params%smpd)
                call stkio_w%open('tmp_cavgs.mrc', params%smpd, 'write', is_ft=.false., box=params%box)
                do icls = 1,params%ncls
                    call img%read(cavgs,icls)
                    call img%fft
                    call img%pad(img_pad, backgr=0.)
                    call img_pad%ifft
                    call stkio_w%write(icls, img_pad)
                enddo
                call stkio_w%close
                call simple_rename('tmp_cavgs.mrc',cavgs)
                call img%kill
                call img_pad%kill
            end subroutine rescale_cavgs

    end subroutine exec_cleanup2D

    subroutine exec_cluster2D_autoscale( self, cline )
        use simple_commander_project, only: scale_project_commander_distr
        use simple_commander_imgproc, only: scale_commander, pspec_int_rank_commander
        class(cluster2D_autoscale_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        ! constants
        integer,               parameter :: MAXITS_STAGE1      = 10
        integer,               parameter :: MAXITS_STAGE1_EXTR = 15
        character(len=STDLEN), parameter :: orig_projfile_bak  = 'orig_bak.simple'
        ! commanders
        type(make_cavgs_commander_distr)    :: xmake_cavgs_distr
        type(make_cavgs_commander)          :: xmake_cavgs
        type(cluster2D_commander_distr)     :: xcluster2D_distr
        type(cluster2D_commander)           :: xcluster2D
        type(pspec_int_rank_commander)      :: xpspec_rank
        type(rank_cavgs_commander)          :: xrank_cavgs
        type(scale_commander)               :: xscale
        type(scale_project_commander_distr) :: xscale_proj
        type(calc_pspec_commander_distr)    :: xcalc_pspec_distr
        ! command lines
        type(cmdline) :: cline_cluster2D_stage1, cline_cluster2D_stage2
        type(cmdline) :: cline_scalerefs, cline_scale
        type(cmdline) :: cline_make_cavgs, cline_rank_cavgs, cline_pspec_rank
        type(cmdline) :: cline_calc_pspec_distr
        ! other variables
        class(parameters), pointer    :: params_ptr => null()
        type(parameters)              :: params
        type(sp_project)              :: spproj, spproj_sc
        character(len=:), allocatable :: projfile_sc, orig_projfile
        character(len=LONGSTRLEN)     :: finalcavgs, finalcavgs_ranked, cavgs, refs_sc
        real     :: scale, trs_stage2
        integer  :: last_iter_stage1, last_iter_stage2
        logical  :: scaling, l_shmem, l_euclid
        call set_cluster2D_defaults( cline )
        call cline%delete('clip')
        ! set shared-memory flag
        if( cline%defined('nparts') )then
            if( nint(cline%get_rarg('nparts')) == 1 )then
                l_shmem = .true.
                call cline%delete('nparts')
            else
                l_shmem = .false.
            endif
        else
            l_shmem = .true.
        endif
        ! master parameters
        call params%new(cline)
        ! report limits used
        write(logfhandle,'(A,F5.1)') '>>> DID SET STARTING  LOW-PASS LIMIT (IN A) TO: ', params%lpstart
        write(logfhandle,'(A,F5.1)') '>>> DID SET HARD      LOW-PASS LIMIT (IN A) TO: ', params%lpstop
        write(logfhandle,'(A,F5.1)') '>>> DID SET CENTERING LOW-PASS LIMIT (IN A) TO: ', params%cenlp
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! read project file
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo')
        orig_projfile = trim(params%projfile)
        ! sanity checks
        if( spproj%get_nptcls() == 0 )then
            THROW_HARD('No particles found in project file: '//trim(params%projfile)//'; exec_cluster2D_autoscale')
        endif
        ! delete any previous solution
        if( .not. spproj%is_virgin_field(params%oritype) .and. trim(params%refine).ne.'inpl' )then
            ! removes previous cluster2D solution (states are preserved)
            call spproj%os_ptcl2D%delete_2Dclustering
            call spproj%write_segment_inside(params%oritype)
        endif
        ! splitting
        if( .not. l_shmem ) call spproj%split_stk(params%nparts, dir=PATH_PARENT)
        ! deal with eo partitioning
        if( spproj%os_ptcl2D%get_nevenodd() == 0 )then
            call spproj%os_ptcl2D%partition_eo
            call spproj%write_segment_inside(params%oritype,params%projfile)
        endif
        ! noise power estimates for objfun = euclid at original sampling
        l_euclid = .false.
        if( cline%defined('objfun') )then
            l_euclid = ( trim(cline%get_carg('objfun')).eq.'euclid' .or. trim(cline%get_carg('objfun')).eq.'prob' .or. trim(cline%get_carg('objfun')).eq.'test' )
            if( l_euclid )then
                cline_calc_pspec_distr  = cline
                call cline_calc_pspec_distr%delete('scale')
                call cline_calc_pspec_distr%set( 'prg', 'calc_pspec' )
                call spproj%os_ptcl2D%set_all2single('w', 1.0)
            endif
        endif
        ! general options planning
        if( params%l_autoscale )then
            ! this workflow executes two stages of CLUSTER2D
            ! Stage 1: down-scaling for fast execution, hybrid extremal/SHC optimisation for
            !          improved population distribution of clusters, no incremental learning,
            cline_cluster2D_stage1 = cline
            call cline_cluster2D_stage1%set('lpstop',     params%lpstart)
            call cline_cluster2D_stage1%set('ptclw','no')
            call cline_cluster2D_stage1%set('objfun','cc') ! cc-based search in first phase
            if( params%l_frac_update )then
                call cline_cluster2D_stage1%delete('update_frac') ! no incremental learning in stage 1
                call cline_cluster2D_stage1%set('maxits', real(MAXITS_STAGE1_EXTR))
            else
                call cline_cluster2D_stage1%set('maxits', real(MAXITS_STAGE1))
            endif
            ! Scaling
            call spproj%scale_projfile(params%smpd_targets2D(2), projfile_sc,&
                &cline_cluster2D_stage1, cline_scale, dir=trim(STKPARTSDIR))
            call spproj%kill
            scale   = cline_scale%get_rarg('scale')
            scaling = basename(projfile_sc) /= basename(orig_projfile)
            if( scaling )then
                if( l_euclid )then
                    ! sigma2 need to be scaled, which will be undone prior
                    ! to execution of make_cavgs at original sampling
                    call cline_calc_pspec_distr%set('scale', 1./scale**2)
                    call xcalc_pspec_distr%execute( cline_calc_pspec_distr )
                endif
                call cline_scale%delete('smpd') !!
                call cline_scale%set('state',1.)
                call simple_mkdir(trim(STKPARTSDIR),errmsg="commander_hlev_wflows :: exec_cluster2D_autoscale;  ")
                call xscale_proj%execute( cline_scale )
                ! rename scaled projfile and stash original project file
                ! such that the scaled project file has the same name as the original and can be followed from the GUI
                call simple_copy_file(orig_projfile, orig_projfile_bak)
                call spproj%read_non_data_segments(projfile_sc)
                call spproj%projinfo%set(1,'projname',get_fbody(orig_projfile,METADATA_EXT,separator=.false.))
                call spproj%projinfo%set(1,'projfile',orig_projfile)
                call spproj%write_non_data_segments(projfile_sc)
                call spproj%kill
                call simple_rename(projfile_sc,orig_projfile)
                deallocate(projfile_sc)
                ! scale references
                if( cline%defined('refs') )then
                    call cline_scalerefs%set('stk', trim(params%refs))
                    refs_sc = 'refs'//trim(SCALE_SUFFIX)//params%ext
                    call cline_scalerefs%set('outstk', trim(refs_sc))
                    call cline_scalerefs%set('smpd', params%smpd)
                    call cline_scalerefs%set('newbox', cline_scale%get_rarg('newbox'))
                    call xscale%execute(cline_scalerefs)
                    call cline_cluster2D_stage1%set('refs',trim(refs_sc))
                endif
            else
                if( l_euclid ) call xcalc_pspec_distr%execute( cline_calc_pspec_distr )
            endif
            ! execution
            call cline_cluster2D_stage1%set('projfile', trim(orig_projfile))
            if( l_shmem )then
                params_ptr  => params_glob
                params_glob => null()
                call xcluster2D%execute(cline_cluster2D_stage1)
                params_glob => params_ptr
                params_ptr  => null()
            else
                call xcluster2D_distr%execute(cline_cluster2D_stage1)
            endif
            last_iter_stage1 = nint(cline_cluster2D_stage1%get_rarg('endit'))
            cavgs            = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage1,3)//params%ext
            ! Stage 2: refinement stage, down-scaling, no extremal updates, optional incremental
            !          learning for acceleration
            cline_cluster2D_stage2 = cline
            call cline_cluster2D_stage2%set('refs', cavgs)
            call cline_cluster2D_stage2%set('startit', real(last_iter_stage1 + 1))
            if( cline%defined('update_frac') )then
                call cline_cluster2D_stage2%set('update_frac', params%update_frac)
            endif
            if( l_euclid )then
                call cline_cluster2D_stage2%set('objfun',      trim(cline%get_carg('objfun')))
                call cline_cluster2D_stage2%set('needs_sigma', 'yes')
            endif
            trs_stage2 = MSK_FRAC * params%mskdiam / (2 * params%smpd_targets2D(2))
            trs_stage2 = min(MAXSHIFT,max(MINSHIFT,trs_stage2))
            call cline_cluster2D_stage2%set('trs', trs_stage2)
            ! execution
            call cline_cluster2D_stage2%set('projfile', trim(orig_projfile))
            if( l_shmem )then
                params_ptr  => params_glob
                params_glob => null()
                call xcluster2D%execute(cline_cluster2D_stage2)
                params_glob => params_ptr
                params_ptr  => null()
            else
                call xcluster2D_distr%execute(cline_cluster2D_stage2)
            endif
            last_iter_stage2 = nint(cline_cluster2D_stage2%get_rarg('endit'))
            finalcavgs       = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//params%ext
            ! Updates project and references
            if( scaling )then
                ! shift modulation
                call spproj_sc%read(orig_projfile)
                call spproj_sc%os_ptcl2D%mul_shifts(1./scale)
                call spproj%read(orig_projfile_bak)
                spproj%os_ptcl2D = spproj_sc%os_ptcl2D
                spproj%os_cls2D  = spproj_sc%os_cls2D
                spproj%os_cls3D  = spproj_sc%os_cls3D
                call spproj%write(orig_projfile_bak)
                call spproj%kill()
                call spproj_sc%kill()
                call simple_rename(orig_projfile_bak,orig_projfile)
                if( l_euclid )then
                    ! adjusts sigma2
                    call scale_group_sigma2_magnitude(last_iter_stage2, scale**2)
                endif
                ! clean stacks
                call simple_rmdir(STKPARTSDIR)
                ! original scale references
                cline_make_cavgs = cline ! ncls is transferred here
                call cline_make_cavgs%delete('autoscale')
                call cline_make_cavgs%delete('balance')
                call cline_make_cavgs%set('prg',      'make_cavgs')
                call cline_make_cavgs%set('projfile', orig_projfile)
                call cline_make_cavgs%set('nparts',   real(params%nparts))
                call cline_make_cavgs%set('refs',     trim(finalcavgs))
                call cline_make_cavgs%delete('wiener') ! to ensure that full Wiener restoration is done for the final cavgs
                call cline_make_cavgs%set('which_iter', real(last_iter_stage2)) ! to ensure masks are generated and used
                if( l_shmem )then
                    params_ptr  => params_glob
                    params_glob => null()
                    call xmake_cavgs%execute(cline_make_cavgs)
                    params_glob => params_ptr
                    params_ptr  => null()
                else
                    call xmake_cavgs_distr%execute(cline_make_cavgs)
                endif
            endif
        else
            ! no auto-scaling
            if( l_shmem )then
                params_ptr  => params_glob
                params_glob => null()
                call xcluster2D%execute(cline)
                params_glob => params_ptr
                params_ptr  => null()
            else
                call xcluster2D_distr%execute(cline)
            endif
            last_iter_stage2 = nint(cline%get_rarg('endit'))
            finalcavgs       = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//params%ext
        endif
        ! adding cavgs & FRCs to project
        params%projfile = trim(orig_projfile)
        call spproj%read( params%projfile )
        call spproj%add_frcs2os_out( trim(FRCS_FILE), 'frc2D')
        call spproj%add_cavgs2os_out(trim(finalcavgs), spproj%get_smpd(), imgkind='cavg')
        call spproj%write_segment_inside('out', params%projfile)
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
            ! rank based on gold-standard resolution estimates
            finalcavgs_ranked = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//'_ranked'//params%ext
            call cline_rank_cavgs%set('projfile', trim(params%projfile))
            call cline_rank_cavgs%set('stk',      trim(finalcavgs))
            call cline_rank_cavgs%set('outstk',   trim(finalcavgs_ranked))
            call xrank_cavgs%execute( cline_rank_cavgs )
        endif
        ! cleanup
        call del_file('start2Drefs'//params%ext)
        call del_file('start2Drefs_even'//params%ext)
        call del_file('start2Drefs_odd'//params%ext)
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER2D NORMAL STOP ****')
    end subroutine exec_cluster2D_autoscale

    subroutine exec_cluster2D_distr( self, cline )
        class(cluster2D_commander_distr), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        ! commanders
        type(make_cavgs_commander_distr) :: xmake_cavgs
        ! command lines
        type(cmdline) :: cline_check_2Dconv
        type(cmdline) :: cline_cavgassemble
        type(cmdline) :: cline_make_cavgs
        type(cmdline) :: cline_calc_sigma
        integer(timer_int_kind)   :: t_init,   t_scheduled,  t_merge_algndocs,  t_cavgassemble,  t_tot
        real(timer_int_kind)      :: rt_init, rt_scheduled, rt_merge_algndocs, rt_cavgassemble, rt_tot
        character(len=STDLEN)     :: benchfname, orig_objfun
        ! other variables
        type(parameters)          :: params
        type(builder)             :: build
        type(qsys_env)            :: qenv
        character(len=LONGSTRLEN) :: refs, refs_even, refs_odd, str, str_iter, finalcavgs
        integer                   :: iter, cnt, iptcl, ptclind, fnr, iter_switch2euclid
        type(chash)               :: job_descr
        real                      :: frac_srch_space
        logical                   :: l_stream, l_switch2euclid, l_ptclw, l_griddingset
        call cline%set('prg','cluster2D')
        call set_cluster2D_defaults( cline )
        ! streaming
        l_stream = .false.
        if( cline%defined('stream') )then
            l_stream = trim(cline%get_carg('stream'))=='yes'
        endif
        call cline%set('stream','no') ! for parameters determination
        ! objfun=euclid logics, part 1
        l_griddingset   = cline%defined('gridding')
        l_switch2euclid = .false.
        if( cline%defined('objfun') )then
            if( trim(cline%get_carg('objfun')).eq.'euclid' .or. trim(cline%get_carg('objfun')).eq.'prob' .or. trim(cline%get_carg('objfun')).eq.'test' )then
                orig_objfun = trim(cline%get_carg('objfun'))
                l_ptclw     = trim(cline%get_carg('ptclw')).eq.'yes'
                if( cline%defined('needs_sigma') )then
                    if(trim(cline%get_carg('needs_sigma')).eq.'yes')then
                        ! were are already doing ML, no need to switch from cc nor calculate noise power
                        if( .not.l_griddingset ) call cline%set('gridding','yes')
                    else
                        l_switch2euclid = .true.
                        call cline%set('objfun','cc')
                        call cline%set('needs_sigma','yes')
                    endif
                else
                    l_switch2euclid = .true.
                    call cline%set('objfun','cc')
                    call cline%set('needs_sigma','yes')
                endif
            endif
        endif
        ! builder & params
        call build%init_params_and_build_spproj(cline, params)
        if( l_stream ) call cline%set('stream','yes')
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
        ! objfun=euclid logics, part 2
        iter_switch2euclid = -1
        if( l_switch2euclid )then
            iter_switch2euclid = params%startit
            if( cline%defined('update_frac') ) iter_switch2euclid = ceiling(1./(params%update_frac+0.001))
        endif
        ! prepare command lines from prototype master
        cline_check_2Dconv      = cline
        cline_cavgassemble      = cline
        cline_make_cavgs        = cline ! ncls is transferred here
        cline_calc_sigma        = cline
        ! initialise static command line parameters and static job description parameters
        call cline_cavgassemble%set(     'prg', 'cavgassemble')
        call cline_make_cavgs%set(       'prg', 'make_cavgs')
        call cline_calc_sigma%set('prg', 'calc_group_sigmas' ) ! required for local call
        ! execute initialiser
        if( .not. cline%defined('refs') )then
            refs             = 'start2Drefs'//params%ext
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
                        if( trim(params%refine).eq.'inpl' )then
                            params%ncls = build%spproj%os_ptcl2D%get_n('class')
                            call job_descr%set('ncls',int2str(params%ncls))
                            call cline%set('ncls', real(params%ncls))
                            call cline_make_cavgs%set('ncls', real(params%ncls))
                            call cline_make_cavgs%delete('tseries')
                            call cline_cavgassemble%set('ncls', real(params%ncls))
                            call cline_make_cavgs%set('refs', params%refs)
                            call xmake_cavgs%execute(cline_make_cavgs)
                        else
                            call selection_from_tseries_imgfile(build%spproj, params%refs, params%box, params%ncls)
                        endif
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
        ! initialise progress monitor
        if(.not. l_stream) call progressfile_init()
        ! main loop
        iter = params%startit - 1
        do
            if( L_BENCH_GLOB )then
                t_init = tic()
                t_tot  = t_init
            endif
            iter = iter + 1
            params%which_iter = iter
            call cline%set(     'which_iter', trim(int2str(params%which_iter)))
            call job_descr%set( 'which_iter', trim(int2str(params%which_iter)))
            str_iter = int2str_pad(iter,3)
            write(logfhandle,'(A)')   '>>>'
            write(logfhandle,'(A,I6)')'>>> ITERATION ', params%which_iter
            write(logfhandle,'(A)')   '>>>'
            ! noise power
            if( trim(params%objfun).eq.'euclid' .or. trim(params%objfun).eq.'prob' .or. trim(params%objfun).eq.'test' .or. l_switch2euclid )then
                call cline_calc_sigma%set('which_iter',real(params%which_iter))
                call qenv%exec_simple_prg_in_queue(cline_calc_sigma, 'CALC_GROUP_SIGMAS_FINISHED')
            endif
            ! cooling of the randomization rate
            params%extr_iter = params%extr_iter + 1
            call job_descr%set('extr_iter', trim(int2str(params%extr_iter)))
            call cline%set('extr_iter', real(params%extr_iter))
            ! updates
            call job_descr%set('refs', trim(refs))
            call job_descr%set('startit', trim(int2str(iter)))
            ! the only FRC we have is from the previous iteration, hence the iter - 1
            call job_descr%set('frcs', trim(FRCS_FILE))
            ! schedule
            if( L_BENCH_GLOB )then
                rt_init = toc(t_init)
                t_scheduled = tic()
            endif
            call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=trim(ALGN_FBODY), array=L_USE_SLURM_ARR)
            call terminate_stream('SIMPLE_DISTR_CLUSTER2D HARD STOP 1')
            ! assemble alignment docs
            if( L_BENCH_GLOB )then
                rt_scheduled = toc(t_scheduled)
                t_merge_algndocs = tic()
            endif
            call build%spproj%merge_algndocs(params%nptcls, params%nparts, 'ptcl2D', ALGN_FBODY)
            if( L_BENCH_GLOB )then
                rt_merge_algndocs = toc(t_merge_algndocs)
                t_cavgassemble = tic()
            endif
            ! assemble class averages
            refs      = trim(CAVGS_ITER_FBODY) // trim(str_iter)            // params%ext
            refs_even = trim(CAVGS_ITER_FBODY) // trim(str_iter) // '_even' // params%ext
            refs_odd  = trim(CAVGS_ITER_FBODY) // trim(str_iter) // '_odd'  // params%ext
            call cline_cavgassemble%set('refs', trim(refs))
            call terminate_stream('SIMPLE_DISTR_CLUSTER2D HARD STOP 2')
            call qenv%exec_simple_prg_in_queue(cline_cavgassemble, 'CAVGASSEMBLE_FINISHED')
            if( L_BENCH_GLOB ) rt_cavgassemble = toc(t_cavgassemble)
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
            ! objfun=euclid, part 3: actual switch
            if( l_switch2euclid .and. (iter==iter_switch2euclid) )then
                write(logfhandle,'(A)')'>>>'
                write(logfhandle,'(A)')'>>> SWITCHING TO OBJFUN=EUCLID'
                call cline%set('objfun', orig_objfun)
                if(.not.l_griddingset )then
                    call cline%set('gridding',     'yes')
                    call job_descr%set('gridding', 'yes')
                endif
                call job_descr%set('objfun', orig_objfun)
                call cline_cavgassemble%set('objfun', orig_objfun)
                if( l_ptclw )then
                    call cline%set('ptclw',    'yes')
                    call job_descr%set('ptclw','yes')
                endif
                params%objfun = orig_objfun
                if( params%objfun .eq. 'euclid' )then
                    params%cc_objfun = OBJFUN_EUCLID
                elseif( params%objfun .eq. 'prob' )then
                    params%cc_objfun = OBJFUN_PROB
                elseif( params%objfun .eq. 'test' )then
                    params%cc_objfun = OBJFUN_TEST
                endif
                l_switch2euclid = .false.
            endif
            if( L_BENCH_GLOB )then
                rt_tot  = toc(t_init)
                benchfname = 'CLUSTER2D_DISTR_BENCH_ITER'//int2str_pad(iter,3)//'.txt'
                call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
                write(fnr,'(a)') '*** TIMINGS (s) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation  : ', rt_init
                write(fnr,'(a,1x,f9.2)') 'scheduled jobs  : ', rt_scheduled
                write(fnr,'(a,1x,f9.2)') 'merge_algndocs  : ', rt_merge_algndocs
                write(fnr,'(a,1x,f9.2)') 'cavgassemble    : ', rt_cavgassemble
                write(fnr,'(a,1x,f9.2)') 'total time      : ', rt_tot
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation  : ', (rt_init/rt_tot)           * 100.
                write(fnr,'(a,1x,f9.2)') 'scheduled jobs  : ', (rt_scheduled/rt_tot)      * 100.
                write(fnr,'(a,1x,f9.2)') 'merge_algndocs  : ', (rt_merge_algndocs/rt_tot) * 100.
                write(fnr,'(a,1x,f9.2)') 'cavgassemble    : ', (rt_cavgassemble/rt_tot)   * 100.
                write(fnr,'(a,1x,f9.2)') '% accounted for : ',&
                    &((rt_init+rt_scheduled+rt_merge_algndocs+rt_cavgassemble)/rt_tot) * 100.
                call fclose(fnr)
            endif
        end do
        ! updates os_out
        finalcavgs = trim(CAVGS_ITER_FBODY)//int2str_pad(iter,3)//params%ext
        call build%spproj%add_cavgs2os_out(trim(finalcavgs), build%spproj%get_smpd(), imgkind='cavg')
        call build%spproj%write_segment_inside('out', params%projfile)
        call qsys_cleanup
        ! report the last iteration on exit
        call cline%delete( 'startit' )
        call cline%set('endit', real(iter))
        ! end gracefully
        call build%spproj_field%kill
        call simple_touch(CLUSTER2D_FINISHED)
        call simple_end('**** SIMPLE_DISTR_CLUSTER2D NORMAL STOP ****')
    end subroutine exec_cluster2D_distr

    subroutine exec_cluster2D( self, cline )
        use simple_strategy2D_matcher, only: cluster2D_exec
        use simple_classaverager
        class(cluster2D_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(make_cavgs_commander)        :: xmake_cavgs
        type(calc_group_sigmas_commander) :: xcalc_group_sigmas
        type(cmdline)              :: cline_make_cavgs
        type(parameters)           :: params
        type(builder), target      :: build
        type(starproject)          :: starproj
        character(len=LONGSTRLEN)  :: finalcavgs, orig_objfun
        logical                    :: converged, l_stream, l_switch2euclid
        integer                    :: startit, ncls_from_refs, lfoo(3), i, cnt, iptcl, ptclind
        call cline%set('oritype', 'ptcl2D')
        if( .not. cline%defined('maxits') ) call cline%set('maxits', 30.)
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
        if( cline%defined('refs') )then
            call find_ldim_nptcls(params%refs, lfoo, ncls_from_refs)
            ! consistency check
            if( params%ncls /=  ncls_from_refs ) THROW_HARD('nrefs /= inputted ncls')
        endif
        l_stream = .false.
        if( cline%defined('stream') )then
            l_stream = trim(cline%get_carg('stream'))=='yes'
        endif
        startit = 1
        if( cline%defined('startit') )startit = params%startit
        if( startit == 1 )call build%spproj_field%clean_updatecnt
        if( params%l_distr_exec )then
            if( .not. cline%defined('outfile') ) THROW_HARD('need unique output file for parallel jobs')
            call cluster2D_exec( cline, startit, converged )
            ! end gracefully
            call simple_end('**** SIMPLE_CLUSTER2D NORMAL STOP ****')
            call qsys_job_finished('simple_commander_cluster2D :: exec_cluster2D')
        else
            if( .not. cline%defined('refs') )then
                cline_make_cavgs = cline ! ncls is transferred here
                params%refs      = 'start2Drefs'//params%ext
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
                            call cline%set('ncls', real(params%ncls))
                            call cline_make_cavgs%set('ncls', real(params%ncls))
                            call cline_make_cavgs%set('refs', params%refs)
                            call xmake_cavgs%execute(cline_make_cavgs)
                        else
                            if( trim(params%refine).eq.'inpl' )then
                                params%ncls = build%spproj%os_ptcl2D%get_n('class')
                                call cline%set('ncls', real(params%ncls))
                                call cline_make_cavgs%set('ncls', real(params%ncls))
                                call cline_make_cavgs%delete('tseries')
                                call cline_make_cavgs%set('refs', params%refs)
                                call xmake_cavgs%execute(cline_make_cavgs)
                            else
                                call selection_from_tseries_imgfile(build%spproj, params%refs, params%box, params%ncls)
                            endif
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
                call cline%set('refs', params%refs)
                ! put back the pointer to builder (no idea why this is needed but it bugs out otherwise)
                build_glob => build
            endif
            params%startit = startit
            params%outfile = 'algndoc'//METADATA_EXT
            ! variable neighbourhood size
            if( cline%defined('extr_iter') )then
                params%extr_iter = params%extr_iter - 1
            else
                params%extr_iter = params%startit - 1
            endif
            ! ML sigmas
            l_switch2euclid = ( params%cc_objfun.eq.OBJFUN_EUCLID .or. params%cc_objfun.eq.OBJFUN_PROB .or. params%cc_objfun.eq.OBJFUN_TEST )
            orig_objfun     = trim(cline%get_carg('objfun'))
            if( cline%defined('needs_sigma') .and. params%l_needs_sigma )then
                ! we are continuing from an ML iterartion
                l_switch2euclid = .false.
            endif
            if( l_switch2euclid )then
                call cline%set('objfun','cc')
                params%objfun        = 'cc'
                params%cc_objfun     = OBJFUN_CC
                params%l_needs_sigma = .true.
            endif
            if( params%l_needs_sigma )then
                params%needs_sigma = 'yes'
                call cline%set('gridding','yes')
                call cline%set('needs_sigma', 'yes')
            endif
            ! initialise progress monitor
            if(.not. l_stream) call progressfile_init()
            do i = 1, params%maxits
                params%which_iter = params%startit
                ! sigmas2
                if( params%l_needs_sigma )then
                    call xcalc_group_sigmas%execute(cline)
                endif
                write(logfhandle,'(A)')   '>>>'
                write(logfhandle,'(A,I6)')'>>> ITERATION ', params%which_iter
                write(logfhandle,'(A)')   '>>>'
                call cline%set('which_iter', real(params%which_iter))
                call cluster2D_exec( cline, params%startit, converged )
                ! objfun=euclid
                if( l_switch2euclid )then
                    params%objfun    = orig_objfun
                    if( params%objfun .eq. 'euclid' )then
                        params%cc_objfun = OBJFUN_EUCLID
                    elseif( params%objfun .eq. 'prob' )then
                        params%cc_objfun = OBJFUN_PROB
                    elseif( params%objfun .eq. 'test' )then
                        params%cc_objfun = OBJFUN_TEST
                    endif
                    l_switch2euclid  = .false.
                endif
                ! cooling of the randomization rate
                params%extr_iter = params%extr_iter + 1
                ! update project with the new orientations (this needs to be here rather than after convergence for the current implementation to work)
                call build%spproj%write_segment_inside(params%oritype)
                ! write cavgs starfile for iteration
                call starproj%export_cls2D(build%spproj, params%which_iter)
                if( converged .or. i == params%maxits )then
                    ! report the last iteration on exit
                    call cline%delete( 'startit' )
                    call cline%set('endit', real(params%startit))
                    call del_file(params%outfile)
                    ! update os_out
                    finalcavgs = trim(CAVGS_ITER_FBODY)//int2str_pad(params%startit,3)//params%ext
                    call build%spproj%add_cavgs2os_out(trim(finalcavgs), build%spproj%get_smpd(), imgkind='cavg')
                    call build%spproj%write_segment_inside('out', params%projfile)
                    exit
                endif
                ! update iteration counter
                params%startit = params%startit + 1
            end do
            if( params%l_needs_sigma .and. (i > 1) )then
                params%which_iter = params%which_iter + 1
                call xcalc_group_sigmas%execute(cline)
                params%which_iter = params%which_iter - 1
            endif
            ! end gracefully
            call simple_touch(CLUSTER2D_FINISHED)
            call simple_end('**** SIMPLE_CLUSTER2D NORMAL STOP ****')
        endif
    end subroutine exec_cluster2D

    subroutine exec_cavgassemble( self, cline )
        use simple_classaverager
        class(cavgassemble_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)  :: params
        type(builder)     :: build
        type(starproject) :: starproj
        real, allocatable :: states(:)
        logical           :: l_stream
        integer           :: iterstr_start, iterstr_end, iter, io_stat
        call cline%set('oritype', 'ptcl2D')
        l_stream = .false.
        if( cline%defined('stream') )then
            l_stream = trim(cline%get_carg('stream'))=='yes'
            call cline%set('stream','no')
        endif
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
        if( l_stream )then
            params_glob%stream = 'yes'
        endif
        call cavger_new
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
        call terminate_stream('SIMPLE_CAVGASSEMBLE HARD STOP 1')
        call cavger_calc_and_write_frcs_and_eoavg(params%frcs, params%which_iter)
        ! classdoc gen needs to be after calc of FRCs
        call cavger_gen2Dclassdoc(build%spproj)
        ! get iteration from which_iter else from refs filename and write cavgs starfile
        if( cline%defined('which_iter') ) then
            call starproj%export_cls2D(build%spproj, params%which_iter)
        else if( cline%defined('refs') .and. index(params%refs, trim(CAVGS_ITER_FBODY)) > 0 ) then
            iterstr_start = index(params%refs, trim(CAVGS_ITER_FBODY)) + 10
            iterstr_end = index(params%refs, trim(params%ext)) - 1
            call str2int(params%refs(iterstr_start:iterstr_end), io_stat, iter)
            call starproj%export_cls2D(build%spproj, iter)
        end if
        ! write references
        call terminate_stream('SIMPLE_CAVGASSEMBLE HARD STOP 2')
        call cavger_write(trim(params%refs),      'merged')
        call cavger_write(trim(params%refs_even), 'even'  )
        call cavger_write(trim(params%refs_odd),  'odd'   )
        call cavger_kill()
        ! write project: cls2D and state congruent cls3D
        call build%spproj%os_cls3D%new(params%ncls, is_ptcl=.false.)
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
        class(rank_cavgs_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters)     :: params
        type(sp_project)     :: spproj
        type(oris)           :: clsdoc_ranked
        type(stack_io)       :: stkio_r, stkio_w
        type(image)          :: img
        integer, allocatable :: order(:)
        real,    allocatable :: res(:)
        integer              :: ldim(3), ncls, icls
        call cline%set('oritype', 'cls2D')
        call params%new(cline)
        call spproj%read_segment(params%oritype, params%projfile)
        call find_ldim_nptcls(params%stk, ldim, ncls)
        params%ncls = ncls
        if( spproj%os_cls2D%get_noris() == params%ncls )then
            ! all we need to do is fetch from classdoc in projfile &
            ! order according to resolution
            call img%new([params%box,params%box,1], params%smpd)
            call clsdoc_ranked%new(params%ncls, is_ptcl=.false.)
            res = spproj%os_cls2D%get_all('res')
            allocate(order(params%ncls))
            order = (/(icls,icls=1,params%ncls)/)
            call hpsort(res, order)
            call stkio_r%open(params%stk, params%smpd, 'read', bufsz=params%ncls)
            call stkio_r%read_whole ! because need asynchronous access
            call stkio_w%open(params%outstk, params%smpd, 'write', box=ldim(1), bufsz=params%ncls)
            do icls=1,params%ncls
                call clsdoc_ranked%set(icls, 'class',     real(order(icls)))
                call clsdoc_ranked%set(icls, 'rank',      real(icls))
                call clsdoc_ranked%set(icls, 'pop',       spproj%os_cls2D%get(order(icls),  'pop'))
                call clsdoc_ranked%set(icls, 'res',       spproj%os_cls2D%get(order(icls),  'res'))
                call clsdoc_ranked%set(icls, 'corr',      spproj%os_cls2D%get(order(icls), 'corr'))
                call clsdoc_ranked%set(icls, 'w',         spproj%os_cls2D%get(order(icls),    'w'))
                call clsdoc_ranked%set(icls, 'specscore', spproj%os_cls2D%get(order(icls), 'specscore'))
                write(logfhandle,'(a,1x,i5,1x,a,1x,i5,1x,a,i5,1x,a,1x,f6.2)') 'CLASS:', order(icls),&
                    &'RANK:', icls ,'POP:', nint(spproj%os_cls2D%get(order(icls), 'pop')),&
                    &'RES:', spproj%os_cls2D%get(order(icls), 'res')
                call flush(logfhandle)
                call stkio_r%get_image(order(icls), img)
                call stkio_w%write(icls, img)
            end do
            call stkio_r%close
            call stkio_w%close
            call clsdoc_ranked%write('classdoc_ranked.txt')
        else
            ! nothing to do
        endif
        ! end gracefully
        call clsdoc_ranked%kill
        call img%kill
        call spproj%kill
        call simple_end('**** SIMPLE_RANK_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_rank_cavgs

    subroutine exec_cluster_cavgs( self, cline )
        use simple_polarizer,        only: polarizer
        use simple_class_frcs,       only: class_frcs
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_aff_prop,         only: aff_prop
        class(cluster_cavgs_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(class_frcs)              :: clsfrcs
        type(image)                   :: img_msk
        type(polarft_corrcalc)        :: pftcc
        type(aff_prop)                :: aprop
        type(polarizer),  allocatable :: cavg_imgs(:), cavg_imgs_good(:)
        character(len=:), allocatable :: cavgsstk, cavgsstk_shifted, classname, frcs_fname
        real,             allocatable :: states(:), orig_cls_inds(:), frc(:), filter(:), clspops(:), clsres(:)
        real,             allocatable :: corrs(:), corrmat_comlin(:,:), corrs_top_ranking(:)
        logical,          allocatable :: l_msk(:,:,:), mask_top_ranking(:), mask_otsu(:), mask_icls(:)
        integer,          allocatable :: order(:), nloc(:), centers(:), labels(:), cntarr(:), clsinds(:), pops(:)
        integer,          allocatable :: labels_mapped(:), classmapping(:)
        logical,          parameter   :: DEBUG = .false.
        integer :: ncls, n, ldim(3), ncls_sel, i, j, icls, cnt, filtsz, pop1, pop2, nsel, ncls_aff_prop, icen, jcen, loc(1)
        real    :: smpd, sdev_noise, simsum, cmin, cmax, pref, corr_icls
        logical :: l_apply_optlp, use_shifted
        ! defaults
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir', 'yes')
        call cline%set('oritype', 'ptcl2D')
        call params%new(cline)
        call spproj%read(params%projfile)
        ! need to turn off CTF since it may be on based on the imported images and now we operate on class averages
        params%ctf = 'no'
        ! get class average stack
        call spproj%get_cavgs_stk(cavgsstk_shifted, ncls, smpd, imgkind='cavg_shifted', fail=.false.)
        use_shifted = .true.
        if( cavgsstk_shifted .eq. NIL ) use_shifted = .false.
        call spproj%get_cavgs_stk(cavgsstk, ncls, smpd)
        call find_ldim_nptcls(cavgsstk, ldim, n)
        write(logfhandle,'(A,I3)') '# classes in stack ', n
        ldim(3) = 1
        if( n /= ncls ) THROW_HARD('Inconsistent # classes in project file vs cavgs stack; exec_cluster_cavgs')
        ! ensure correct smpd/box in params class
        params%smpd = smpd
        params%box  = ldim(1)
        ! threshold based on states/population/resolution
        states = spproj%os_cls2D%get_all('state')
        if( .not. DEBUG )then
            clspops = spproj%os_cls2D%get_all('pop')
            clsres  = spproj%os_cls2D%get_all('res')
            where( clsres  >= params%lpthres     ) states = 0.
            where( clspops <  real(MINCLSPOPLIM) ) states = 0.
        endif
        ! find out how many selected class averages initially
        ncls_sel = count(states > 0.5)
        write(logfhandle,'(A,I3)') '# classes left after standard rejection ', ncls_sel
        ! keep track of the original class indices
        allocate(clsinds(ncls))
        clsinds = (/(i,i=1,ncls)/)
        clsinds = pack(clsinds, mask=states > 0.5)
        ! get FRCs
        l_apply_optlp = .false.
        call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
        if( file_exists(frcs_fname) )then
            call clsfrcs%read(frcs_fname)
            l_apply_optlp = .true.
        endif
        ! create the stuff needed in the loop
        filtsz = clsfrcs%get_filtsz()
        allocate(cavg_imgs(ncls_sel), frc(filtsz), filter(filtsz))
        ! prep mask
        call img_msk%new([params%box,params%box,1], params%smpd)
        img_msk = 1.
        call img_msk%mask(params%msk, 'hard')
        l_msk = img_msk%bin2logical()
        call img_msk%kill
        write(logfhandle,'(A)') '>>> PREPARING CLASS AVERAGES'
        cnt = 0
        do i = 1 , ncls
            if( states(i) > 0.5 )then
                cnt = cnt + 1
            else
                cycle
            endif
            call cavg_imgs(cnt)%new(ldim, smpd, wthreads=.false.)
            if( use_shifted )then
                call cavg_imgs(cnt)%read(cavgsstk_shifted, i)
            else
                call cavg_imgs(cnt)%read(cavgsstk, i)
                if( DEBUG ) call cavg_imgs(cnt)%write('1_read_cavgs.mrc', cnt)
            endif
            ! FRC-based filter
            if( l_apply_optlp )then
                call clsfrcs%frc_getter(i, params%hpind_fsc, params%l_phaseplate, frc)
                if( any(frc > 0.143) )then
                    call fsc2optlp_sub(clsfrcs%get_filtsz(), frc, filter)
                    where( filter > TINY ) filter = sqrt(filter) ! because the filter is applied to the average not the even or odd
                    call cavg_imgs(cnt)%fft()
                    call cavg_imgs(cnt)%apply_filter_serial(filter)
                    call cavg_imgs(cnt)%ifft()
                endif
            endif
            ! normalization
            call cavg_imgs(cnt)%norm_noise(l_msk, sdev_noise)
            if( DEBUG ) call cavg_imgs(cnt)%write('3_norm_noisealized.mrc', cnt)
            ! mask
            call cavg_imgs(cnt)%mask(params%msk, 'soft')
            if( DEBUG ) call cavg_imgs(cnt)%write('4_masked.mrc', cnt)
        end do
        ! create the polarft_corrcalc object
        params%kfromto(1) = max(2, calc_fourier_index(params%hp, params%box, params%smpd))
        params%kfromto(2) =        calc_fourier_index(params%lp, params%box, params%smpd)
        call pftcc%new(ncls_sel, [1,1], params%kfromto)
        if( trim(params%bin_cls).ne.'no' .and. .not. DEBUG )then
            ! initialize polarizer for the first image, then copy it to the rest
            call cavg_imgs(1)%init_polarizer(pftcc, params%alpha)
            !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
            do icls = 1, ncls_sel
                call cavg_imgs(icls)%fft()
                call cavg_imgs(1)%polarize(pftcc, cavg_imgs(icls), icls, isptcl=.false., iseven=.true.) ! 2 polar coords
                call cavg_imgs(icls)%ifft()
            end do
            !$omp end parallel do
            ! some allocations needed
            nsel = ceiling(0.1 * real(ncls_sel))
            allocate(nloc(nsel), order(ncls_sel),                                    source=0)
            allocate(corrmat_comlin(ncls_sel,ncls_sel), corrs_top_ranking(ncls_sel), source=-1.)
            allocate(mask_top_ranking(ncls_sel),                                     source=.true.)
            write(logfhandle,'(A)') '>>> CALCULATING COMMON-LINE CORRELATION MATRIX'
            !$omp parallel do default(shared) private(i,j) schedule(dynamic) proc_bind(close)
            do i = 1, ncls_sel - 1
                do j = i + 1, ncls_sel
                    ! Common line correlation
                    corrmat_comlin(i,j) = pftcc%genmaxcorr_comlin(i,j)
                    corrmat_comlin(j,i) = corrmat_comlin(i,j)
                enddo
                corrmat_comlin(i,i) = 1.
            enddo
            !$omp end parallel do
            ! take care of the last diagonal element
            corrmat_comlin(ncls_sel,ncls_sel) = 1.
            write(logfhandle,'(A)') '>>> GOOD/BAD CLASSIFICATION AND RANKING OF CLASS AVERAGES'
            ! average the nsel best correlations (excluding self) to create a scoring function for garbage removal
            do i = 1, ncls_sel
                mask_top_ranking    = .true.
                mask_top_ranking(i) = .false. ! to remove the diagonal element
                corrs = pack(corrmat_comlin(i,:), mask_top_ranking)
                nloc = maxnloc(corrs, nsel)
                corrs_top_ranking(i) = 0.
                do j = 1, nsel
                    corrs_top_ranking(i) = corrs_top_ranking(i) + corrs(nloc(j))
                end do
                corrs_top_ranking(i) = corrs_top_ranking(i) / real(nsel)
            end do
            ! use Otsu's algorithm to remove the junk
            call otsu(ncls_sel, corrs_top_ranking, mask_otsu)
            pop1 = count(      mask_otsu)
            pop2 = count(.not. mask_otsu)
            write(logfhandle,*) 'average corr cluster 1: ', sum(corrs_top_ranking, mask=      mask_otsu) / real(pop1), ' pop ', pop1
            write(logfhandle,*) 'average corr cluster 2: ', sum(corrs_top_ranking, mask=.not. mask_otsu) / real(pop2), ' pop ', pop2
            pop1 = 0
            pop2 = 0
            do i = 1, ncls_sel
                if( mask_otsu(i) )then
                    pop1 = pop1 + 1
                    call cavg_imgs(i)%write('good.mrc', pop1)
                else
                    pop2 = pop2 + 1
                    call cavg_imgs(i)%write('bad.mrc',  pop2)
                endif
            end do
            ! rank
            order = (/(i,i=1,ncls_sel)/)
            call hpsort(corrs_top_ranking, order)
            do i = 1, ncls_sel
                call cavg_imgs(order(i))%write('ranked.mrc', i)
            end do
            ! update the states array according to the good/bad classification
            do i = 1, ncls_sel
                if( .not. mask_otsu(i) ) states(clsinds(i)) = 0.
            end do
            ! update the class indices (bookkeeping is getting hairy)
            deallocate(clsinds)
            allocate(clsinds(ncls))
            clsinds = (/(i,i=1,ncls)/)
            clsinds = pack(clsinds, mask=states > 0.5)
            ! toss out the bad ones
            allocate(cavg_imgs_good(pop1))
            pop1 = 0
            do i = 1, ncls_sel
                if( mask_otsu(i) )then
                    pop1 = pop1 + 1
                    call cavg_imgs_good(pop1)%copy(cavg_imgs(i))
                endif
            end do
            ncls_sel = pop1
            write(logfhandle,'(A,I3)') '# classes left after binary classification rejection ', ncls_sel
        else
            allocate(cavg_imgs_good(ncls_sel))
            do i = 1, ncls_sel
                call cavg_imgs_good(i)%copy(cavg_imgs(i))
                if( DEBUG ) call cavg_imgs_good(i)%write('5_good_cavgs.mrc', i)
            end do
        endif
        if( trim(params%bin_cls).ne.'only' )then
            ! re-calculate the similarity matrix
            write(logfhandle,'(A)') '>>> CALCULATING COMMON-LINE CORRELATION MATRIX'
            ! re-create the polarft_corrcalc object
            call pftcc%kill
            call pftcc%new(ncls_sel, [1,1], params%kfromto)
            ! initialize polarizer for the first image, then copy it to the rest
            call cavg_imgs_good(1)%init_polarizer(pftcc, params%alpha)
            !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
            do icls = 1, ncls_sel
                call cavg_imgs_good(icls)%fft()
                call cavg_imgs_good(1)%polarize(pftcc, cavg_imgs_good(icls), icls, isptcl=.false., iseven=.true.) ! 2 polar coords
                call cavg_imgs_good(icls)%ifft()
            end do
            !$omp end parallel do
            if( allocated(corrmat_comlin) ) deallocate(corrmat_comlin)
            allocate(corrmat_comlin(ncls_sel,ncls_sel), source=-1.)
            !$omp parallel do default(shared) private(i,j) schedule(dynamic) proc_bind(close)
            do i = 1, ncls_sel - 1
                do j = i + 1, ncls_sel
                    ! Common line orrelations
                    corrmat_comlin(i,j) = pftcc%genmaxcorr_comlin(i,j)
                    corrmat_comlin(j,i) = corrmat_comlin(i,j)
                enddo
                corrmat_comlin(i,i) = 1.
            enddo
            !$omp end parallel do
            ! take care of the last diagonal element
            corrmat_comlin(ncls_sel,ncls_sel) = 1.
            ! calculate a preference that generates a small number of clusters
            call analyze_smat(corrmat_comlin, .false., cmin, cmax)
            pref = cmin - (cmax - cmin)
            write(logfhandle,'(A)') '>>> CLUSTERING CLASS AVERAGES WITH AFFINITY PROPAGATION'
            call aprop%new(ncls_sel, corrmat_comlin, pref=pref)
            call aprop%propagate(centers, labels, simsum)
            ncls_aff_prop = size(centers)
            write(logfhandle,'(A,I3)') '>>> # CLUSTERS FOUND BY AFFINITY PROPAGATION (AP): ', ncls_aff_prop
            allocate(cntarr(ncls_aff_prop), source=0)
            if( cline%defined('ncls') )then
                write(logfhandle,'(A,I3)') '>>> RE-MAPPING THE HIGHEST POPULATED AP CLUSTERS TO INPUT NCLS: ', params%ncls
                ! identify params%ncls highest populated clusters
                if( allocated(mask_icls) ) deallocate(mask_icls)
                allocate( mask_icls(ncls_sel), pops(ncls_aff_prop) )
                do icls = 1, ncls_aff_prop
                    ! mask out the cluster
                    where( labels == icls )
                        mask_icls = .true.
                    elsewhere
                        mask_icls = .false.
                    end where
                    pops(icls) = count(mask_icls)
                end do
                if( allocated(nloc) ) deallocate(nloc)
                allocate(nloc(params%ncls), source=0)
                nloc = maxnloc(real(pops), params%ncls)
                ! remap the AP clusters
                if( allocated(corrs) ) deallocate(corrs)
                allocate(corrs(params%ncls), classmapping(ncls_aff_prop))
                ! put back the diagonal elements of the corrmat
                do i = 1, ncls_sel
                    corrmat_comlin(i,i) = 1.
                end do
                ! make a per cluster mapping
                do icls = 1, ncls_aff_prop
                    do i = 1,params%ncls
                        corrs(i) = corrmat_comlin(nloc(i),icls)
                    end do
                    loc = maxloc(corrs)
                    classmapping(icls) = nloc(loc(1))
                    write(logfhandle,*) 'mapped cluster ', icls, ' to ', classmapping(icls)
                end do
                ! do the mapping
                allocate(labels_mapped(ncls_sel), source = 0)
                do i = 1, ncls_sel
                    do icls = 1, ncls_aff_prop
                        if( labels(i) == icls ) labels_mapped(i) = classmapping(icls)
                    end do
                end do
                if( any(labels_mapped == 0) ) THROW_HARD('class mapping did not work, class=0 not allowed')
                ! re-map agaian
                cnt = 0
                do i = 1, ncls_aff_prop
                    if( any(labels_mapped == i) )then
                        cnt = cnt + 1
                        where( labels_mapped == classmapping(i) ) labels = cnt
                    endif
                end do
                deallocate(labels_mapped)
                ! update # AP clusters
                ncls_aff_prop = params%ncls
            endif
            ! put back the original (unprocessed) images
            cnt = 0
            do i = 1 , ncls
                if( states(i) > 0.5 )then
                    cnt = cnt + 1
                else
                    cycle
                endif
                call cavg_imgs(cnt)%read(cavgsstk, i)
            end do
            ! write the classes
            do icls = 1, ncls_aff_prop
                ! make a filename for the class
                do i=1,ncls_sel
                    if( labels(i) == icls )then
                        classname = 'class'//int2str_pad(icls,3)//trim(STK_EXT)
                        cntarr(labels(i)) = cntarr(labels(i)) + 1
                        call cavg_imgs(i)%write(classname, cntarr(labels(i)))
                    endif
                end do
            end do
            ! calculate average common line correlation of the AP clusters
            if( allocated(mask_icls) ) deallocate(mask_icls)
            allocate( mask_icls(ncls_sel) )
            do icls = 1, ncls_aff_prop
                ! mask out the cluster
                where( labels == icls )
                    mask_icls = .true.
                elsewhere
                    mask_icls = .false.
                end where
                ! calculate the average common line correlation of the cluster
                corr_icls = 0.
                cnt       = 0
                do i = 1, ncls_sel - 1
                    do j = i + 1, ncls_sel
                        if( mask_icls(i) .and. mask_icls(j) )then
                            corr_icls = corr_icls + corrmat_comlin(i,j)
                            cnt       = cnt + 1
                        endif
                    end do
                end do
                corr_icls = corr_icls / real(cnt)
                write(logfhandle,*) 'average corr AP cluster '//int2str_pad(icls,3)//': ', corr_icls, ' pop ', count(mask_icls)
            end do
            ! print the correlations between the cluster centers
            write(logfhandle,'(A)') '>>> COMMON LINE CORRELATIONS BETEEN CLUSTER CENTERS'
            do i = 1,ncls_aff_prop - 1
                icen = centers(i)
                do j = i + 1,ncls_aff_prop
                    jcen = centers(j)
                    write(logfhandle,*) 'corr btw cluster  '//int2str_pad(i,3)//' & '//int2str_pad(j,3)//': ', corrmat_comlin(icen,jcen)
                end do
            end do
        endif
        ! destruct
        call spproj%kill
        call clsfrcs%kill
        call pftcc%kill
        call aprop%kill
        do icls=1,ncls_sel
            call cavg_imgs(icls)%kill
        end do
        do icls=1,ncls_sel
            call cavg_imgs_good(icls)%kill
        end do
        ! deallocate
        deallocate(cavg_imgs, cavg_imgs_good)
        ! only need to deallocate stuff that was not explicitly allocated in this scope
        if( allocated(states)        ) deallocate(states)
        if( allocated(orig_cls_inds) ) deallocate(orig_cls_inds)
        if( allocated(clspops)       ) deallocate(clspops)
        if( allocated(clsres)        ) deallocate(clsres)
        if( allocated(corrs)         ) deallocate(corrs)
        if( allocated(l_msk)         ) deallocate(l_msk)
        if( allocated(mask_otsu)     ) deallocate(mask_otsu)
        if( allocated(centers)       ) deallocate(centers)
        if( allocated(labels)        ) deallocate(labels)
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
                call imgs_class(i)%set_rmat(rmat_rot,.false.)
            end do
            !$omp end parallel do
            ! make a filename for the class
            classname = 'class'//int2str_pad(icls,5)//trim(STK_EXT)
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
        class(cmdline), intent(inout) :: cline
        class(oris),    intent(inout) :: os
        type(parameters)  :: params
        type(convergence) :: conv
        logical :: converged, l_stream
        call cline%set('oritype', 'ptcl2D')
        call params%new(cline)
        l_stream = .false.
        if( cline%defined('stream') )then
            l_stream = trim(cline%get_carg('stream'))=='yes'
        endif
        ! convergence check
        converged = conv%check_conv2D(cline, os, os%get_n('class'), params%msk)
        ! Update progress file
        if(.not. l_stream) call progressfile_update(conv%get('progress'))
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

    subroutine terminate_stream( msg )
        character(len=*), intent(in) :: msg
        if(trim(params_glob%async).eq.'yes')then
            if( file_exists(TERM_STREAM) )then
                call simple_end('**** '//trim(msg)//' ****', print_simple=.false.)
            endif
        endif
    end subroutine terminate_stream

end module simple_commander_cluster2D
