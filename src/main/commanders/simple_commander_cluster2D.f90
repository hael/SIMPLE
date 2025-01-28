! concrete commander: cluster2D for simultanous 2D alignment and clustering of single-particle images
module simple_commander_cluster2D
include 'simple_lib.f08'
use simple_builder,           only: builder, build_glob
use simple_cmdline,           only: cmdline
use simple_commander_base,    only: commander_base
use simple_parameters,        only: parameters, params_glob
use simple_sp_project,        only: sp_project
use simple_qsys_env,          only: qsys_env
use simple_image,             only: image
use simple_stack_io,          only: stack_io
use simple_starproject,       only: starproject
use simple_commander_imgproc, only: scale_commander
use simple_exec_helpers,      only: set_shmem_flag, set_master_num_threads
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
public :: prob_tab2D_commander
public :: prob_tab2D_commander_distr
public :: make_cavgs_commander_distr
public :: make_cavgs_commander
public :: cavgassemble_commander
public :: rank_cavgs_commander
public :: cluster_cavgs_commander
public :: write_classes_commander
public :: ppca_denoise_classes_commander
public :: partition_cavgs_commander
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

type, extends(commander_base) :: prob_tab2D_commander
  contains
    procedure :: execute      => exec_prob_tab2D
end type prob_tab2D_commander

type, extends(commander_base) :: prob_tab2D_commander_distr
  contains
    procedure :: execute      => exec_prob_tab2D_distr
end type prob_tab2D_commander_distr

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

type, extends(commander_base) :: ppca_denoise_classes_commander
  contains
    procedure :: execute      => exec_ppca_denoise_classes
end type ppca_denoise_classes_commander

type, extends(commander_base) :: partition_cavgs_commander
  contains
    procedure :: execute      => exec_partition_cavgs
end type partition_cavgs_commander

contains

    subroutine exec_make_cavgs_distr( self, cline )
        class(make_cavgs_commander_distr), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters)             :: params
        type(builder)                :: build
        type(cmdline)                :: cline_cavgassemble
        type(qsys_env)               :: qenv
        type(chash)                  :: job_descr
        type(make_cavgs_commander)   :: xmk_cavgs_shmem
        type(cavgassemble_commander) :: xcavgassemble
        integer :: ncls_here, nthr_here
        logical :: l_shmem
        call cline%set('wiener', 'full')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('ml_reg')  ) call cline%set('ml_reg',      'no')
        if( (cline%defined('ncls')).and. cline%defined('nspace') )then
            THROW_HARD('NCLS and NSPACE cannot be both defined!')
        endif
        if( cline%defined('nspace') )then
            if( trim(cline%get_carg('oritype')).eq.'ptcl2D' )then
                THROW_HARD('NSPACE & PTCL2D are incompatible!')
            endif
            call cline%set('oritype', 'ptcl3D')
        else
            call cline%set('oritype', 'ptcl2D')
        endif
        if( cline%defined('nparts') )then
            l_shmem = cline%get_iarg('nparts') == 1
        else
            l_shmem = .true.
        endif
        ! deal with # threads for the master process
        if( .not.l_shmem ) call set_master_num_threads(nthr_here, 'CLUSTER2D')
        ! parse parameters & project
        call build%init_params_and_build_spproj(cline, params)
        if( cline%defined('nspace') )then
            ! handled in exec_make_cavgs
        else
            ncls_here = build%spproj_field%get_n('class')
            if( .not. cline%defined('ncls') )then
                call cline%set('ncls', ncls_here)
                params%ncls = ncls_here
            endif
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        if( l_shmem  )then
            call xmk_cavgs_shmem%execute_safe(cline)
            return
        endif
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! prepare command lines from prototype master
        cline_cavgassemble = cline
        call cline_cavgassemble%set('prg',  'cavgassemble')
        call cline_cavgassemble%set('nthr',  nthr_here)
        if( trim(params%oritype).eq.'ptcl3D' )then
            call cline_cavgassemble%set('ncls', params%nspace)
        endif
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR)
        ! assemble class averages
        call xcavgassemble%execute_safe(cline_cavgassemble)
        ! end
        call qsys_cleanup
        call simple_end('**** SIMPLE_DISTR_MAKE_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_make_cavgs_distr

    subroutine exec_make_cavgs( self, cline )
        use simple_classaverager
        class(make_cavgs_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer :: ncls_here, icls
        logical :: l_shmem
        if( (cline%defined('ncls')).and. cline%defined('nspace') )then
            THROW_HARD('NCLS and NSPACE cannot be both defined!')
        endif
        if( cline%defined('nspace') )then
            if( trim(cline%get_carg('oritype')).eq.'ptcl2D' )then
                THROW_HARD('NSPACE & PTCL2D are incompatible!')
            endif
            call cline%set('oritype', 'ptcl3D')
            call cline%set('ncls',    cline%get_iarg('nspace'))
        else
            call cline%set('oritype', 'ptcl2D')
        endif
        ! set shared-memory flag
        l_shmem = set_shmem_flag( cline )
        if( l_shmem .and. .not. cline%defined('refs') ) THROW_HARD('need input refs (filename) for shared-memory execution')
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
        if( L_VERBOSE_GLOB ) write(logfhandle,'(a)') '>>> GENERATING CLUSTER CENTERS'
        ! deal with the orientations
        if( trim(params%oritype).eq.'ptcl3D' )then
            ! 3D class averages
            call build%eulspace%new(params%nspace, is_ptcl=.false.)
            call build%pgrpsyms%build_refspiral(build%eulspace)
            call build%spproj%os_ptcl3D%set_projs(build%eulspace)
            call build%spproj%os_ptcl3D%proj2class
        else
            ! 2D
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
        endif
        ! setup weights
        call build%spproj_field%calc_hard_weights2D(params%frac, params%ncls)
        ! even/odd partitioning
        if( build%spproj_field%get_nevenodd() == 0 ) call build%spproj_field%partition_eo
        ! write
        if( l_shmem )then
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        else
            if( params%part .eq. 1 ) call build%spproj%write_segment_inside(params%oritype, params%projfile)
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
            select case(trim(params%oritype))
                case('ptcl2D')
                    call build%spproj%add_frcs2os_out( trim(FRCS_FILE), 'frc2D')
                    call build%spproj%add_cavgs2os_out(trim(params%refs), build%spproj%get_smpd(), imgkind='cavg')
                    call build%spproj%write_segment_inside('out', params%projfile)
                case('ptcl3D')
                    do icls = 1,params%nspace
                        call build%spproj%os_cls3D%set_euler(icls, build%eulspace%get_euler(icls))
                    enddo
                    if( cline%defined('outfile') )then
                        call build%spproj%os_cls3D%write(params%outfile)
                    else
                        call build%spproj%os_cls3D%write('cls3D_oris.txt')
                    endif
                    call build%spproj%write_segment_inside('cls3D', params%projfile)
                    call build%spproj%add_frcs2os_out( trim(FRCS_FILE), 'frc3D')
                    call build%spproj%add_cavgs2os_out(trim(params%refs), build%spproj%get_smpd(), imgkind='cavg3D')
                    call build%spproj%write_segment_inside('out', params%projfile)
                case DEFAULT
                    THROW_HARD('Unsupported ORITYPE: '//trim(params%oritype))
            end select
        else
            ! write partial sums
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
        use simple_class_frcs,        only: class_frcs
        class(cleanup2D_commander_hlev), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        ! commanders
        type(cluster2D_commander_distr)  :: xcluster2D_distr
        type(cluster2D_commander)        :: xcluster2D ! shared-memory implementation
        type(scale_commander)            :: xscale
        type(rank_cavgs_commander)       :: xrank_cavgs
        type(calc_pspec_commander_distr) :: xcalc_pspec_distr
        ! command lines
        type(cmdline)                    :: cline_cluster2D1, cline_cluster2D2
        type(cmdline)                    :: cline_rank_cavgs, cline_scalerefs
        type(cmdline)                    :: cline_calc_pspec_distr
        ! other variables
        type(parameters)                 :: params
        type(sp_project)                 :: spproj
        type(class_frcs)                 :: frcs
        character(len=LONGSTRLEN)        :: finalcavgs, finalcavgs_ranked, cavgs, refs_sc
        real                             :: scale_factor, lp1, lp2, cenlp, smpd_target
        integer                          :: last_iter
        logical                          :: l_scaling, l_shmem, l_euclid
        ! parameters
        integer, parameter    :: MINBOX      = 128
        real,    parameter    :: MINITS      = 5.
        real,    parameter    :: MAXITS      = 15.
        real                  :: SMPD_TARGET_DEFAULT = 3.
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('ncls')      ) call cline%set('ncls',         200)
        if( .not. cline%defined('center')    ) call cline%set('center',      'no')
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale',  'yes')
        if( .not. cline%defined('refine')    ) call cline%set('refine',  'greedy_smpl')
        if( .not. cline%defined('oritype')   ) call cline%set('oritype', 'ptcl2D')
        if( .not. cline%defined('wiener')    ) call cline%set('wiener',    'full')
        if( .not. cline%defined('objfun')    ) call cline%set('objfun',  'euclid')
        if( .not. cline%defined('ml_reg')    ) call cline%set('ml_reg',      'no')
        if( .not. cline%defined('masscen')   ) call cline%set('masscen',    'yes')
        if( .not. cline%defined('sh_first')  ) call cline%set('sh_first',    'no')
        call cline%set('stream', 'no')
        ! set shared-memory flag
        l_shmem = set_shmem_flag( cline )
        ! parse parameters
        call params%new(cline)
        if( .not. cline%defined('maxits') )then
            params%maxits = nint(MAXITS)
            call cline%set('maxits', MAXITS)
        endif
        if( .not. cline%defined('mskdiam') )then
            params%mskdiam = (real(params%box) * 0.8) * params%smpd
            call cline%set('mskdiam', params%mskdiam)
            write(logfhandle,'(A,F6.1)')'>>> MASK DIAMETER (Angstroms): ',params%mskdiam
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
        ! eo partitioning
        if( spproj%os_ptcl2D%get_nevenodd() == 0 )then
            call spproj%os_ptcl2D%partition_eo
            call spproj%write_segment_inside(params%oritype,params%projfile)
        endif
        ! resolutions limits
        call mskdiam2lplimits(params%mskdiam, params%lpstart, params%lpstop, cenlp)
        if( cline%defined('lp') )then
            lp1 = max(params%lp, params%lpstart)
            lp2 = params%lp
        else
            lp1 = params%lpstart
            lp2 = (params%lpstart+params%lpstop)/2.
        endif
        smpd_target = min(lp2/2. * (2./3.), SMPD_TARGET_DEFAULT)
        ! Cropped dimensions
        scale_factor     = 1.0
        params%smpd_crop = params%smpd
        params%box_crop  = params%box
        params%msk_crop  = params%msk
        l_scaling        = .false.
        if( params%l_autoscale .and. params%box >= MINBOX )then
            call autoscale(params%box, params%smpd, smpd_target, params%box_crop, params%smpd_crop, scale_factor, minbox=MINBOX)
            l_scaling       = params%box_crop < params%box
            if( l_scaling )then
                params%msk_crop = round2even(params%msk * scale_factor)
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',params%box_crop
            endif
        endif
        ! objfun = euclid
        l_euclid = .false.
        if( cline%defined('objfun') )then
            l_euclid = trim(cline%get_carg('objfun')).eq.'euclid'
            if( l_euclid )then
                cline_calc_pspec_distr  = cline
                call cline_calc_pspec_distr%set( 'prg', 'calc_pspec' )
                call spproj%os_ptcl2D%set_all2single('w', 1.0)
                call spproj%write_segment_inside(params%oritype)
                call xcalc_pspec_distr%execute( cline_calc_pspec_distr )
            endif
        endif
        ! Clustering command lines & cropping
        cline_cluster2D1 = cline
        cline_cluster2D2 = cline
        call cline_cluster2D1%set('smpd_crop', params%smpd_crop)
        call cline_cluster2D1%set('box_crop',  params%box_crop)
        call cline_cluster2D2%set('smpd_crop', params%smpd_crop)
        call cline_cluster2D2%set('box_crop',  params%box_crop)
        ! resolution limits
        call cline_cluster2D1%set('lp', lp1)
        call cline_cluster2D2%set('lp', lp2)
        if( .not.cline%defined('cenlp') )then
            call cline_cluster2D1%set('cenlp', cenlp)
            call cline_cluster2D2%set('cenlp', cenlp)
        endif
        ! first stage
        ! down-scaling for fast execution, greedy optimisation, no incremental learning, no centering
        call cline_cluster2D1%set('prg', 'cluster2D')
        call cline_cluster2D1%set('maxits',   MINITS)
        call cline_cluster2D1%set('minits',   MINITS)
        call cline_cluster2D1%set('center',     'no')
        call cline_cluster2D1%set('autoscale',  'no')
        call cline_cluster2D1%set('objfun',     'cc')
        call cline_cluster2D1%set('ml_reg','no')
        if( l_euclid ) call cline_cluster2D1%set('cc_iters', MINITS)
        call cline_cluster2D1%delete('update_frac')
        ! second stage
        ! down-scaling for fast execution, greedy optimisation
        call cline_cluster2D2%set('prg', 'cluster2D')
        call cline_cluster2D2%set('autoscale',  'no')
        call cline_cluster2D2%set('trs',    MINSHIFT)
        if( .not.cline%defined('maxits') )then
            call cline_cluster2D2%set('maxits', MAXITS)
        endif
        call cline_cluster2D2%set('minits', min(MINITS+3,MAXITS))
        if( l_euclid )then
            call cline_cluster2D2%set('objfun',   trim(cline%get_carg('objfun')))
            call cline_cluster2D2%set('cc_iters', 0)
        else
            call cline_cluster2D2%set('objfun', 'cc')
        endif
        ! optional non-uniform filtering
        if( cline%defined('update_frac') )call cline_cluster2D2%set('update_frac',params%update_frac)
        ! scale references
        if( l_scaling )then
            if( cline%defined('refs') )then
                refs_sc = 'refs'//trim(SCALE_SUFFIX)//params%ext
                call cline_scalerefs%set('stk',    params%refs)
                call cline_scalerefs%set('outstk', refs_sc)
                call cline_scalerefs%set('smpd',   params%smpd)
                call cline_scalerefs%set('newbox', params%box_crop)
                call xscale%execute_safe(cline_scalerefs)
                call cline_cluster2D1%set('refs',  refs_sc)
            endif
        endif
        ! initialise progress monitor
        call progressfile_init()
        ! execution stage 1
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A,F6.1)') '>>> STAGE 1, LOW-PASS LIMIT: ',lp1
        write(logfhandle,'(A)') '>>>'
        if( l_shmem )then
            call xcluster2D%execute_safe(cline_cluster2D1)
        else
            call xcluster2D_distr%execute(cline_cluster2D1)
        endif
        last_iter  = cline_cluster2D1%get_iarg('endit')
        finalcavgs = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter,3)//params%ext
        if( .not. file_exists(trim(finalcavgs)) ) THROW_HARD('File '//trim(finalcavgs)//' does not exist')
        ! execution stage 2
        if( cline%defined('maxits') )then
            if( last_iter < params%maxits )then
                write(logfhandle,'(A)') '>>>'
                write(logfhandle,'(A,F6.1)') '>>> STAGE 2, LOW-PASS LIMIT: ',lp2
                write(logfhandle,'(A)') '>>>'
                call cline_cluster2D2%set('startit',  last_iter+1)
                call cline_cluster2D2%set('refs',     trim(finalcavgs))
                if( l_shmem )then
                    call xcluster2D%execute_safe(cline_cluster2D2)
                else
                    call xcluster2D_distr%execute(cline_cluster2D2)
                endif
                last_iter  = cline_cluster2D2%get_iarg('endit')
                finalcavgs = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter,3)//params%ext
            endif
        endif
        ! update project
        call cline%set('endit', last_iter)
        if( l_scaling )then
            call rescale_cavgs(finalcavgs)
            cavgs = add2fbody(finalcavgs,params%ext,'_even')
            call rescale_cavgs(cavgs)
            cavgs = add2fbody(finalcavgs,params%ext,'_odd')
            call rescale_cavgs(cavgs)
            call spproj%read_segment('out', params%projfile)
            call spproj%read_segment('cls2D', params%projfile)
            call spproj%read_segment('cls3D', params%projfile)
            call spproj%add_cavgs2os_out(trim(finalcavgs), params%smpd, imgkind='cavg')
            call frcs%read(FRCS_FILE)
            call frcs%pad(params%smpd, params%box)
            call frcs%write(FRCS_FILE)
            call spproj%add_frcs2os_out(FRCS_FILE, 'frc2D')
            call frcs%kill
            ! transfer 2D shift parameters to 3D
            call spproj%read_segment('ptcl2D', params%projfile)
            call spproj%read_segment('ptcl3D', params%projfile)
            call spproj%os_ptcl3D%transfer_2Dshifts(spproj%os_ptcl2D)
            call spproj%write(params%projfile)
        else
            call spproj%read_segment('ptcl2D', params%projfile)
            call spproj%read_segment('ptcl3D', params%projfile)
            call spproj%os_ptcl3D%transfer_2Dshifts(spproj%os_ptcl2D)
            call spproj%write_segment_inside('ptcl3D', params%projfile)
        endif
        call spproj%kill
        ! ranking
        finalcavgs_ranked = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter,3)//'_ranked'//params%ext
        call cline_rank_cavgs%set('projfile', params%projfile)
        call cline_rank_cavgs%set('stk',      finalcavgs)
        call cline_rank_cavgs%set('outstk',   finalcavgs_ranked)
        call xrank_cavgs%execute_safe(cline_rank_cavgs)
        ! end gracefully
        call simple_end('**** SIMPLE_CLEANUP2D NORMAL STOP ****')

        contains

            subroutine rescale_cavgs(cavgs)
                character(len=*), intent(in) :: cavgs
                type(image)    :: img, img_pad
                type(stack_io) :: stkio_w
                integer        :: icls
                call img%new([params%box_crop,params%box_crop,1],params%smpd_crop)
                call img_pad%new([params%box,params%box,1],params%smpd)
                call stkio_w%open('tmp_cavgs.mrc', params%smpd, 'write', is_ft=.false., box=params%box)
                do icls = 1,params%ncls
                    call img%read(cavgs,icls)
                    call img%fft
                    call img%pad(img_pad, backgr=0., antialiasing=.false.)
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
        use simple_commander_imgproc, only: pspec_int_rank_commander
        class(cluster2D_autoscale_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        ! constants
        integer, parameter :: MAXITS_STAGE1      = 10
        integer, parameter :: MAXITS_STAGE1_EXTR = 15
        integer, parameter :: MINBOX             = 88
        ! commanders
        type(make_cavgs_commander_distr)    :: xmake_cavgs_distr
        type(make_cavgs_commander)          :: xmake_cavgs
        type(cluster2D_commander_distr)     :: xcluster2D_distr
        type(cluster2D_commander)           :: xcluster2D
        type(pspec_int_rank_commander)      :: xpspec_rank
        type(rank_cavgs_commander)          :: xrank_cavgs
        type(calc_pspec_commander_distr)    :: xcalc_pspec_distr
        type(scale_commander)               :: xscale
        ! command lines
        type(cmdline) :: cline_cluster2D_stage1, cline_cluster2D_stage2
        type(cmdline) :: cline_scalerefs, cline_calc_pspec_distr
        type(cmdline) :: cline_make_cavgs, cline_rank_cavgs, cline_pspec_rank
        ! other variables
        class(parameters), pointer    :: params_ptr => null()
        type(parameters)              :: params
        type(sp_project)              :: spproj
        character(len=LONGSTRLEN)     :: finalcavgs, finalcavgs_ranked, cavgs, refs_sc
        real     :: scale, trs_stage2, smpd_target
        integer  :: last_iter_stage1, last_iter_stage2
        logical  :: l_scaling, l_shmem, l_euclid, l_cc_iters
        call set_cluster2D_defaults( cline )
        if( .not.cline%defined('objfun') ) call cline%set('objfun', 'cc')
        call cline%delete('clip')
        l_cc_iters = cline%defined('cc_iters')
        ! shared memory flag
        l_shmem = set_shmem_flag( cline )
        ! master parameters
        call params%new(cline)
        ! report limits used
        write(logfhandle,'(A,F5.1)') '>>> DID SET STARTING  LOW-PASS LIMIT (IN A) TO: ', params%lpstart
        write(logfhandle,'(A,F5.1)') '>>> DID SET HARD      LOW-PASS LIMIT (IN A) TO: ', params%lpstop
        write(logfhandle,'(A,F5.1)') '>>> DID SET CENTERING LOW-PASS LIMIT (IN A) TO: ', params%cenlp
        smpd_target = max(2.*params%smpd_crop,min(params%smpd_targets2D(2),params%lpstop))
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! read project file
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo')
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
        ! Cropped dimensions
        l_scaling = .false.
        scale     = 1.0
        params%smpd_crop = params%smpd
        params%box_crop  = params%box
        params%msk_crop  = params%msk
        if( params%l_autoscale .and. (params%box > MINBOX) )then
            call autoscale(params%box, params%smpd, smpd_target, params%box_crop, params%smpd_crop, scale, minbox=MINBOX)
            l_scaling = params%box_crop < params%box
            if( l_scaling )then
                params%msk_crop = round2even(params%msk * scale)
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',params%box_crop
            endif
        endif
        ! noise power estimates for objfun = euclid at original sampling
        l_euclid = .false.
        if( cline%defined('objfun') )then
            l_euclid = trim(cline%get_carg('objfun')).eq.'euclid'
            if( l_euclid )then
                cline_calc_pspec_distr  = cline
                call cline_calc_pspec_distr%set( 'prg', 'calc_pspec' )
                call spproj%os_ptcl2D%set_all2single('w', 1.0)
                call spproj%write_segment_inside(params%oritype)
                call xcalc_pspec_distr%execute( cline_calc_pspec_distr )
            endif
        endif
        ! Clustering command lines
        cline_cluster2D_stage1 = cline
        cline_cluster2D_stage2 = cline
        call cline_cluster2D_stage1%set('smpd_crop', params%smpd_crop)
        call cline_cluster2D_stage1%set('box_crop',  params%box_crop)
        call cline_cluster2D_stage2%set('smpd_crop', params%smpd_crop)
        call cline_cluster2D_stage2%set('box_crop',  params%box_crop)
        if( l_scaling )then
            ! scale references
            if( cline%defined('refs') )then
                refs_sc = 'refs'//trim(SCALE_SUFFIX)//params%ext
                call cline_scalerefs%set('stk',    params%refs)
                call cline_scalerefs%set('outstk', refs_sc)
                call cline_scalerefs%set('smpd',   params%smpd)
                call cline_scalerefs%set('newbox', params%box_crop)
                call xscale%execute_safe(cline_scalerefs)
                call cline_cluster2D_stage1%set('refs', refs_sc)
            endif
        endif
        ! this workflow executes two stages of CLUSTER2D
        ! Stage 1: down-scaling for fast execution, hybrid extremal/SHC optimisation for
        !          improved population distribution of clusters, no incremental learning
        if( l_euclid )then
            call cline_cluster2D_stage1%set('objfun', trim(cline%get_carg('objfun')))
        else
            call cline_cluster2D_stage1%set('objfun', 'cc') ! cc-based search in first phase
        endif
        call cline_cluster2D_stage1%set('lpstop',     params%lpstart)
        call cline_cluster2D_stage1%set('ml_reg',     'no')
        if( params%l_update_frac )then
            call cline_cluster2D_stage1%delete('update_frac') ! no incremental learning in stage 1
            call cline_cluster2D_stage1%set('maxits', MAXITS_STAGE1_EXTR)
            if( l_euclid )then
                if( l_cc_iters )then
                    params%cc_iters = min(params%cc_iters,MAXITS_STAGE1_EXTR)
                    call cline_cluster2D_stage1%set('cc_iters', params%cc_iters)
                else
                    call cline_cluster2D_stage1%set('cc_iters', MAXITS_STAGE1)
                endif
            endif
        else
            call cline_cluster2D_stage1%set('maxits', MAXITS_STAGE1)
            if( l_euclid )then
                if( l_cc_iters )then
                    params%cc_iters = min(params%cc_iters,MAXITS_STAGE1)
                    call cline_cluster2D_stage1%set('cc_iters', params%cc_iters)
                else
                    call cline_cluster2D_stage1%set('cc_iters', MAXITS_STAGE1)
                endif
            endif
        endif
        ! execution
        if( l_shmem )then
            params_ptr  => params_glob
            params_glob => null()
            call xcluster2D%execute(cline_cluster2D_stage1)
            params_glob => params_ptr
            params_ptr  => null()
        else
            call xcluster2D_distr%execute(cline_cluster2D_stage1)
        endif
        last_iter_stage1 = cline_cluster2D_stage1%get_iarg('endit')
        cavgs            = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage1,3)//params%ext
        ! Stage 2: refinement stage, little extremal updates, optional incremental
        !          learning for acceleration
        call cline_cluster2D_stage2%delete('cc_iters')
        call cline_cluster2D_stage2%set('refs',    cavgs)
        call cline_cluster2D_stage2%set('startit', last_iter_stage1+1)
        if( params%l_update_frac )then
            call cline_cluster2D_stage2%set('update_frac', params%update_frac)
        endif
        if( l_euclid )then
            call cline_cluster2D_stage2%set('objfun',   cline%get_carg('objfun'))
            call cline_cluster2D_stage2%set('cc_iters', 0.)
        endif
        trs_stage2 = MSK_FRAC * params%mskdiam / (2. * params%smpd_targets2D(2))
        trs_stage2 = min(MAXSHIFT,max(MINSHIFT,trs_stage2))
        call cline_cluster2D_stage2%set('trs', trs_stage2)
        ! for testing
        if( cline%defined('extr_iter') )then
            call cline_cluster2D_stage2%set('extr_iter', cline_cluster2D_stage1%get_iarg('extr_iter'))
        endif
        ! execution
        if( l_shmem )then
            params_ptr  => params_glob
            params_glob => null()
            call xcluster2D%execute(cline_cluster2D_stage2)
            params_glob => params_ptr
            params_ptr  => null()
        else
            call xcluster2D_distr%execute(cline_cluster2D_stage2)
        endif
        last_iter_stage2 = cline_cluster2D_stage2%get_iarg('endit')
        finalcavgs       = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//params%ext
        ! Updates project and references
        if( l_scaling )then
            ! original scale references
            cline_make_cavgs = cline ! ncls is transferred here
            call cline_make_cavgs%delete('autoscale')
            call cline_make_cavgs%delete('balance')
            call cline_make_cavgs%set('prg',      'make_cavgs')
            call cline_make_cavgs%set('nparts',   params%nparts)
            call cline_make_cavgs%set('refs',     finalcavgs)
            call cline_make_cavgs%delete('wiener') ! to ensure that full Wiener restoration is done for the final cavgs
            call cline_make_cavgs%set('which_iter', last_iter_stage2) ! to ensure masks are generated and used
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
        ! adding cavgs & FRCs to project
        call spproj%read( params%projfile )
        call spproj%add_frcs2os_out( trim(FRCS_FILE), 'frc2D')
        call spproj%add_cavgs2os_out(trim(finalcavgs), params%smpd, imgkind='cavg')
        call spproj%write_segment_inside('out', params%projfile)
        call spproj%os_ptcl3D%transfer_2Dshifts(spproj%os_ptcl2D)
        call spproj%write_segment_inside('ptcl3D', params%projfile)
        ! clean
        call spproj%kill()
        ! ranking
        if( trim(params%tseries).eq.'yes' )then
            ! rank based on maximum of power spectrum
            call cline_pspec_rank%set('mkdir',   'no')
            call cline_pspec_rank%set('moldiam', params%moldiam)
            call cline_pspec_rank%set('nthr',    params%nthr)
            call cline_pspec_rank%set('smpd',    params%smpd)
            call cline_pspec_rank%set('stk',     finalcavgs)
            if( cline%defined('lp_backgr') ) call cline_pspec_rank%set('lp_backgr', params%lp_backgr)
            call xpspec_rank%execute_safe(cline_pspec_rank)
        else
            ! rank based on gold-standard resolution estimates
            finalcavgs_ranked = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//'_ranked'//params%ext
            call cline_rank_cavgs%set('projfile', params%projfile)
            call cline_rank_cavgs%set('stk',      finalcavgs)
            call cline_rank_cavgs%set('outstk',   finalcavgs_ranked)
            call xrank_cavgs%execute_safe( cline_rank_cavgs )
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
        type(make_cavgs_commander_distr)  :: xmake_cavgs
        type(scale_commander)             :: xscale
        type(calc_group_sigmas_commander) :: xcalc_group_sigmas
        type(cavgassemble_commander)      :: xcavgassemble
        type(prob_tab2D_commander_distr)  :: xprob_tab2D_distr
        ! command lines
        type(cmdline) :: cline_check_2Dconv
        type(cmdline) :: cline_cavgassemble
        type(cmdline) :: cline_make_cavgs
        type(cmdline) :: cline_calc_sigma
        type(cmdline) :: cline_scalerefs
        type(cmdline) :: cline_prob_tab2D_distr
        integer(timer_int_kind)   :: t_init,   t_scheduled,  t_merge_algndocs,  t_cavgassemble,  t_tot
        real(timer_int_kind)      :: rt_init, rt_scheduled, rt_merge_algndocs, rt_cavgassemble, rt_tot
        character(len=STDLEN)     :: benchfname, orig_objfun
        ! other variables
        type(parameters)          :: params
        type(builder)             :: build
        type(qsys_env)            :: qenv
        character(len=LONGSTRLEN) :: refs, refs_even, refs_odd, str, str_iter, finalcavgs, refs_sc
        integer                   :: iter, cnt, iptcl, ptclind, fnr, iter_switch2euclid
        type(chash)               :: job_descr
        real                      :: frac_srch_space
        integer                   :: nthr_here
        logical                   :: l_stream, l_switch2euclid, l_griddingset, l_converged, l_ml_reg, l_scale_inirefs
        call cline%set('prg','cluster2D')
        call set_cluster2D_defaults( cline )
        ! streaming
        l_stream = .false.
        if( cline%defined('stream') )then
            l_stream = trim(cline%get_carg('stream'))=='yes'
        endif
        call cline%set('stream','no') ! for parameters determination
        ! objective functions part 1
        l_griddingset   = cline%defined('gridding')
        l_switch2euclid = .false.
        if( cline%defined('objfun') )then
            if( trim(cline%get_carg('objfun')).eq.'euclid' )then
                orig_objfun     = trim(cline%get_carg('objfun'))
                l_switch2euclid = .true.
            endif
        endif
        ! deal with # threads for the master process
        call set_master_num_threads(nthr_here, 'CLUSTER2D')
        ! builder & params
        call build%init_params_and_build_spproj(cline, params)
        if( l_stream ) call cline%set('stream','yes')
        ! objective functions part 2: scheduling
        l_ml_reg = params%l_ml_reg
        iter_switch2euclid = -1
        if( l_switch2euclid )then
            if( params%cc_iters < 1 )then
                ! already performing euclidian-based optimization, no switch
                iter_switch2euclid = params%startit
                l_switch2euclid    = .false.
                if( .not.l_griddingset ) call cline%set('gridding','yes')
            else
                ! switching objective function from cc_iters+1
                call cline%set('objfun','cc')
                call cline%set('ml_reg','no')
                params%cc_objfun   = OBJFUN_CC
                params%objfun      = 'cc'
                iter_switch2euclid = params%startit
                if( cline%defined('cc_iters') ) iter_switch2euclid = params%cc_iters
                params%l_needs_sigma = .false.
                params%needs_sigma   = 'no'
                params%ml_reg        = 'no'
                params%l_ml_reg      = .false.
            endif
        else
            ! Correlation only, but sigma2 calculated from iter_switch2euclid
            params%l_needs_sigma = .false.
            params%needs_sigma   = 'no'
            params%ml_reg        = 'no'
            params%l_ml_reg      = .false.
            call cline%set('ml_reg','no')
            if( cline%defined('cc_iters') ) iter_switch2euclid = params%cc_iters
        endif
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
        cline_check_2Dconv      = cline
        cline_cavgassemble      = cline
        cline_make_cavgs        = cline ! ncls is transferred here
        cline_calc_sigma        = cline
        ! initialise static command line parameters and static job description parameters
        call cline_cavgassemble%set('prg', 'cavgassemble')
        call cline_make_cavgs%set(  'prg', 'make_cavgs')
        call cline_calc_sigma%set(  'prg', 'calc_group_sigmas')
        ! execute initialiser
        if( .not. cline%defined('refs') )then
            refs             = 'start2Drefs'//params%ext
            params%refs      = trim(refs)
            params%refs_even = 'start2Drefs_even'//params%ext
            params%refs_odd  = 'start2Drefs_odd'//params%ext
            l_scale_inirefs  = .false.
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
                                call build%spproj%os_ptcl2D%set(ptclind, 'class', cnt)
                            end do
                        end do
                        call job_descr%set('ncls',int2str(params%ncls))
                        call cline%set('ncls', params%ncls)
                        call cline_make_cavgs%set('ncls', params%ncls)
                        call cline_cavgassemble%set('ncls', params%ncls)
                        call cline_make_cavgs%set('refs', params%refs)
                        call xmake_cavgs%execute(cline_make_cavgs)
                        l_scale_inirefs = .false.
                    else
                        if( trim(params%refine).eq.'inpl' )then
                            params%ncls = build%spproj%os_ptcl2D%get_n('class')
                            call job_descr%set('ncls',int2str(params%ncls))
                            call cline%set('ncls', params%ncls)
                            call cline_make_cavgs%set('ncls', params%ncls)
                            call cline_make_cavgs%delete('tseries')
                            call cline_cavgassemble%set('ncls', params%ncls)
                            call cline_make_cavgs%set('refs', params%refs)
                            call xmake_cavgs%execute(cline_make_cavgs)
                            l_scale_inirefs = .false.
                        else
                            call selection_from_tseries_imgfile(build%spproj, params%refs, params%box, params%ncls)
                            l_scale_inirefs = .true.
                        endif
                    endif
                else
                    select case(trim(params%cls_init))
                    case('ptcl')
                        ! initialization from raw images
                        call random_selection_from_imgfile(build%spproj, params%refs, params%box, params%ncls)
                        l_scale_inirefs = .true.
                    case('rand')
                        ! from noise
                        call noise_imgfile(params%refs, params%ncls, params%box_crop, params%smpd_crop)
                        l_scale_inirefs = .false.
                    case('randcls')
                        if(.not.cline%defined('ncls')) THROW_HARD('NCLS must be provide with CLS_INIT=RANDCLS')
                        ! initialization from random classes
                        do iptcl=1,params%nptcls
                            if( build%spproj_field%get_state(iptcl) == 0 ) cycle
                            call build%spproj_field%set(iptcl, 'class', irnd_uni(params%ncls))
                            call build%spproj_field%set(iptcl, 'w',     1.0)
                            call build%spproj_field%e3set(iptcl,ran3()*360.0)
                        end do
                        call build%spproj%write_segment_inside(params%oritype, params%projfile)
                        call cline_make_cavgs%set('refs', params%refs)
                        call xmake_cavgs%execute(cline_make_cavgs)
                        l_scale_inirefs = .false.
                    case DEFAULT
                        THROW_HARD('Unsupported mode of initial class generation CLS_INIT='//trim(params%cls_init))
                    end select
                endif
            else
                call cline_make_cavgs%set('refs', params%refs)
                call xmake_cavgs%execute_safe(cline_make_cavgs)
                l_scale_inirefs = .false.
            endif
            ! scale references to box_crop
            if( l_scale_inirefs )then
                refs_sc = 'refs'//trim(SCALE_SUFFIX)//params%ext
                call cline_scalerefs%set('stk',    trim(params%refs))
                call cline_scalerefs%set('outstk', trim(refs_sc))
                call cline_scalerefs%set('smpd',   params%smpd)
                call cline_scalerefs%set('newbox', params%box_crop)
                call cline_scalerefs%set('nthr',   nthr_here)
                call xscale%execute(cline_scalerefs)
                call simple_rename(refs_sc, params%refs)
            endif
            call copy_imgfile(trim(params%refs), trim(params%refs_even), params%smpd_crop, [1,params%ncls])
            call copy_imgfile(trim(params%refs), trim(params%refs_odd),  params%smpd_crop, [1,params%ncls])
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
            ! cooling of the randomization rate
            params%extr_iter = params%extr_iter + 1
            call job_descr%set('extr_iter', trim(int2str(params%extr_iter)))
            call cline%set('extr_iter', params%extr_iter)
            ! objfun function part 3: activate sigma2 calculation
            if( iter==iter_switch2euclid )then
                call cline%set('needs_sigma','yes')
                call job_descr%set('needs_sigma','yes')
                params%needs_sigma   = 'yes'
                params%l_needs_sigma = .true.
            endif
            ! build probability table
            if( str_has_substr(params%refine, 'prob') )then
                cline_prob_tab2D_distr = cline
                call cline_prob_tab2D_distr%set('refs',      refs)
                call cline_prob_tab2D_distr%set('frcs',      FRCS_FILE)
                call cline_prob_tab2D_distr%set('startit',   iter)
                call xprob_tab2D_distr%execute_safe(cline_prob_tab2D_distr)
            endif
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
            call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=trim(ALGN_FBODY), array=L_USE_SLURM_ARR, extra_params=params)
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
            call cline_cavgassemble%set('refs', refs)
            call cline_cavgassemble%set('nthr', nthr_here)
            call terminate_stream('SIMPLE_DISTR_CLUSTER2D HARD STOP 2')
            call xcavgassemble%execute_safe(cline_cavgassemble)
            if( L_BENCH_GLOB ) rt_cavgassemble = toc(t_cavgassemble)
            ! objfun=euclid, part 4: sigma2 consolidation
            if( params%l_needs_sigma )then
                call cline_calc_sigma%set('which_iter', params%which_iter+1)
                call cline_calc_sigma%set('nthr',       nthr_here)
                call xcalc_group_sigmas%execute_safe(cline_calc_sigma)
            endif
            ! print out particle parameters per iteration
            if( trim(params%print_corrs).eq.'yes' )then
                call build%spproj_field%write('ptcl2D_'//int2str_pad(params%which_iter,2)//'.txt')
            endif
            ! check convergence
            call check_2Dconv(cline_check_2Dconv, build%spproj_field)
            frac_srch_space = 0.
            if( iter > 1 ) frac_srch_space = cline_check_2Dconv%get_rarg('frac_srch')
            ! the below activates shifting
            if( iter > 3 .and. (frac_srch_space >= FRAC_SH_LIM .or. cline_check_2Dconv%defined('trs')) )then
                if( .not.job_descr%isthere('trs') )then
                    ! activates shift search
                    str = real2str(cline_check_2Dconv%get_rarg('trs'))
                    call job_descr%set('trs', trim(str) )
                endif
            endif
            l_converged = (iter >= params%minits) .and. (cline_check_2Dconv%get_carg('converged').eq.'yes')
            if( l_converged .or. iter==params%maxits ) then
                exit
            else
                call cline_check_2Dconv%delete('converged')
            endif
            ! objfun=euclid, part 5:  actual switch
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
                params%objfun = trim(orig_objfun)
                if( params%objfun == 'euclid' ) params%cc_objfun = OBJFUN_EUCLID
                l_switch2euclid = .false.
                if( l_ml_reg )then
                    call cline%set('ml_reg',     'yes')
                    call job_descr%set('ml_reg', 'yes')
                    params%ml_reg   = 'yes'
                    params%l_ml_reg = .true.
                endif
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
        call cline%set('endit', iter)
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
        type(scale_commander)             :: xscale
        type(prob_tab2D_commander_distr)  :: xprob_tab2D_distr
        type(cmdline)              :: cline_make_cavgs, cline_scalerefs, cline_prob_tab2D
        type(parameters)           :: params
        type(builder), target      :: build
        type(starproject)          :: starproj
        character(len=LONGSTRLEN)  :: finalcavgs, orig_objfun, refs_sc, fname
        integer                    :: startit, ncls_from_refs, lfoo(3), i, cnt, iptcl, ptclind
        integer                    :: iter_switch2euclid, j, io_stat, funit, class_ind, class_max
        logical                    :: converged, l_stream, l_switch2euclid, l_griddingset, l_ml_reg, l_scale_inirefs
        real,    allocatable       :: corrs(:), corrs_all(:), class_all(:)
        integer, allocatable       :: order(:), class_cnt(:)
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
        if( (startit == 1) .and. (.not.str_has_substr(params%refine,'prob')) )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
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
                l_scale_inirefs  = .false.
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
                                    call build%spproj%os_ptcl2D%set(ptclind, 'class', cnt)
                                end do
                            end do
                            call cline%set('ncls', params%ncls)
                            call cline_make_cavgs%set('ncls', params%ncls)
                            call cline_make_cavgs%set('refs', params%refs)
                            call xmake_cavgs%execute(cline_make_cavgs)
                            l_scale_inirefs  = .false.
                        else
                            if( trim(params%refine).eq.'inpl' )then
                                params%ncls = build%spproj%os_ptcl2D%get_n('class')
                                call cline%set('ncls', params%ncls)
                                call cline_make_cavgs%set('ncls', params%ncls)
                                call cline_make_cavgs%delete('tseries')
                                call cline_make_cavgs%set('refs', params%refs)
                                call xmake_cavgs%execute(cline_make_cavgs)
                                l_scale_inirefs  = .false.
                            else
                                call selection_from_tseries_imgfile(build%spproj, params%refs, params%box, params%ncls)
                                l_scale_inirefs  = .true.
                            endif
                        endif
                    else
                        select case(trim(params%cls_init))
                        case('ptcl')
                            ! initialization from raw images
                            call random_selection_from_imgfile(build%spproj, params%refs, params%box, params%ncls)
                            l_scale_inirefs  = .true.
                        case('rand')
                            ! from noise
                            call noise_imgfile(params%refs, params%ncls, params%box_crop, params%smpd_crop)
                            l_scale_inirefs  = .false.
                        case('randcls')
                            if(.not.cline%defined('ncls')) THROW_HARD('NCLS must be provide with CLS_INIT=RANDCLS ')
                            ! initialization from random classes
                            do iptcl=1,params%nptcls
                                if( build%spproj_field%get_state(iptcl) == 0 ) cycle
                                call build%spproj_field%set(iptcl, 'class', irnd_uni(params%ncls))
                                call build%spproj_field%set(iptcl, 'w',     1.0)
                                call build%spproj_field%e3set(iptcl,ran3()*360.0)
                            end do
                            call build%spproj%write_segment_inside(params%oritype, params%projfile)
                            call cline_make_cavgs%set('refs', params%refs)
                            call xmake_cavgs%execute(cline_make_cavgs)
                            l_scale_inirefs  = .false.
                        case DEFAULT
                            THROW_HARD('Unsupported mode of initial class generation CLS_INIT='//trim(params%cls_init))
                        end select
                    endif
                    ! scale references to box_crop
                    if( l_scale_inirefs )then
                        refs_sc = 'refs'//trim(SCALE_SUFFIX)//params%ext
                        call cline_scalerefs%set('stk',    trim(params%refs))
                        call cline_scalerefs%set('outstk', trim(refs_sc))
                        call cline_scalerefs%set('smpd',   params%smpd)
                        call cline_scalerefs%set('newbox', params%box_crop)
                        call xscale%execute(cline_scalerefs)
                        call simple_rename(refs_sc, params%refs)
                    endif
                    call copy_imgfile(trim(params%refs), trim(params%refs_even), params%smpd_crop, [1,params%ncls])
                    call copy_imgfile(trim(params%refs), trim(params%refs_odd),  params%smpd_crop, [1,params%ncls])
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
                ! all good
            else
                params%extr_iter = params%startit
            endif
            ! objective functions
            l_switch2euclid = params%cc_objfun==OBJFUN_EUCLID
            orig_objfun     = trim(cline%get_carg('objfun'))
            l_griddingset   = cline%defined('gridding')
            l_ml_reg        = params%l_ml_reg
            iter_switch2euclid = -1
            if( l_switch2euclid )then
                if( params%cc_iters < 1 )then
                    ! already performing euclidian-based optimization, no switch
                    if( .not.file_exists(sigma2_star_from_iter(params%startit)) )then
                        THROW_HARD('Sigma2 file does not exists: '//sigma2_star_from_iter(params%startit))
                    endif
                    iter_switch2euclid = params%startit
                    l_switch2euclid    = .false.
                    if( .not.l_griddingset ) call cline%set('gridding','yes')
                else
                    ! switching objective function from cc_iters+1
                    call cline%set('objfun','cc')
                    call cline%set('ml_reg','no')
                    params%cc_objfun   = OBJFUN_CC
                    params%objfun      = 'cc'
                    iter_switch2euclid = params%startit
                    if( cline%defined('cc_iters') ) iter_switch2euclid = params%cc_iters
                    params%l_needs_sigma = .false.
                    params%needs_sigma   = 'no'
                    params%ml_reg        = 'no'
                    params%l_ml_reg      = .false.
                endif
            else
                ! Correlation only, but sigma2 calculated from iter_switch2euclid
                params%l_needs_sigma = .false.
                params%needs_sigma   = 'no'
                if( cline%defined('cc_iters') ) iter_switch2euclid = params%cc_iters
                call cline%set('ml_reg','no')
                params%ml_reg        = 'no'
                params%l_ml_reg      = .false.
            endif
            ! initialise progress monitor
            if(.not. l_stream) call progressfile_init()
            ! Main loop
            do i = 1, params%maxits
                params%which_iter = params%startit
                write(logfhandle,'(A)')   '>>>'
                write(logfhandle,'(A,I6)')'>>> ITERATION ', params%which_iter
                write(logfhandle,'(A)')   '>>>'
                call cline%set('which_iter', params%which_iter)
                if( params%which_iter == iter_switch2euclid )then
                    call cline%set('needs_sigma','yes')
                    params%needs_sigma   = 'yes'
                    params%l_needs_sigma = .true.
                endif
                ! refine=prob
                if( str_has_substr(params%refine, 'prob') )then
                    cline_prob_tab2D = cline
                    call cline_prob_tab2D%set('prg',        'prob_tab2D')
                    call cline_prob_tab2D%set('which_iter', params%which_iter)
                    call cline_prob_tab2D%set('startit',    params%which_iter)
                    call cline_prob_tab2D%set('refs',       params%refs)
                    call cline_prob_tab2D%set('frcs',       FRCS_FILE)
                    call xprob_tab2D_distr%execute( cline_prob_tab2D )
                endif
                ! stochastic search
                call cluster2D_exec( cline, params%startit, converged )
                ! objective functions
                if( params%l_needs_sigma )then
                    params%which_iter = params%which_iter + 1
                    call xcalc_group_sigmas%execute(cline)
                    params%which_iter = params%which_iter - 1
                endif
                if( l_switch2euclid .and. params%which_iter==iter_switch2euclid )then
                    params%objfun = trim(orig_objfun)
                    if( params%objfun == 'euclid' ) params%cc_objfun = OBJFUN_EUCLID
                    l_switch2euclid = .false.
                    if( l_ml_reg )then
                        call cline%set('ml_reg','yes')
                        params%ml_reg        = 'yes'
                        params%l_ml_reg      = .true.
                    endif
                endif
                ! cooling of the randomization rate
                params%extr_iter = params%extr_iter + 1
                ! update project with the new orientations (this needs to be here rather than after convergence for the current implementation to work)
                call build%spproj%write_segment_inside(params%oritype)
                ! write cavgs starfile for iteration
                call starproj%export_cls2D(build%spproj, params%which_iter)
                ! exit condition
                converged = converged .and. (params%which_iter >= params%minits)
                if( converged .or. params%which_iter >= params%maxits )then
                    ! report the last iteration on exit
                    call cline%delete( 'startit' )
                    call cline%set('endit', params%startit)
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
            ! print CSV file of correlation vs particle number
            if( trim(params_glob%print_corrs).eq.'yes' )then
                class_all = build%spproj%os_ptcl2D%get_all('class')
                class_max = int(maxval(class_all))
                ! counting the number of particles in each class
                allocate(class_cnt(class_max))
                do class_ind = 1, class_max
                    class_cnt(class_ind) = 0
                    do j = 1, size(class_all)
                        if( class_all(j) == class_ind ) class_cnt(class_ind) = class_cnt(class_ind) + 1
                    enddo
                enddo
                ! print all sorted corrs
                corrs_all = build%spproj%os_ptcl2D%get_all('corr')
                order     = (/(j,j=1,size(corrs_all))/)
                call hpsort(corrs_all, order)
                fname = 'ptcls_vs_cavgs_corrs_iter'// int2str(params%which_iter) //'.csv'
                call fopen(funit, trim(fname), 'replace', 'unknown', iostat=io_stat, form='formatted')
                call fileiochk('cluster2D fopen failed: '//trim(fname), io_stat)
                write(funit,*) 'PTCL_INDEX'//CSV_DELIM//'CORR'
                do j = 1,size(corrs_all)
                    write(funit,*) int2str(order(j))//CSV_DELIM//real2str(corrs_all(j))
                end do
                call fclose(funit)
                ! print sorted corrs for each class
                do class_ind = 1, class_max
                    if( allocated(corrs) ) deallocate(corrs)
                    allocate(corrs(class_cnt(class_ind)), source=0.)
                    cnt = 0
                    do j = 1, size(corrs_all)
                        if( class_all(j) == class_ind )then
                            cnt = cnt + 1
                            corrs(cnt) = corrs_all(j)
                        endif
                    enddo
                    if( allocated(order) ) deallocate(order)
                    order = (/(j,j=1,class_cnt(class_ind))/)
                    call hpsort(corrs, order)
                    fname = 'ptcls_vs_cavgs_corrs_cls_'// int2str(class_ind) //'.csv'
                    call fopen(funit, trim(fname), 'replace', 'unknown', iostat=io_stat, form='formatted')
                    call fileiochk('cluster2D fopen failed: '//trim(fname), io_stat)
                    write(funit,*) 'PTCL_INDEX'//CSV_DELIM//'CORR'
                    do j = 1,size(corrs)
                        write(funit,*) int2str(order(j))//CSV_DELIM//real2str(corrs(j))
                    end do
                    call fclose(funit)
                enddo
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
        integer           :: iterstr_start, iterstr_end, iter, io_stat, icls
        if( .not.cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
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
        call cavger_gen2Dclassdoc(build%spproj) ! populates the cls2D field in project
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
        ! updates project
        select case(trim(params%oritype))
            case('ptcl2D')
                ! cls2D and state congruent cls3D
                call build%spproj%os_cls3D%new(params%ncls, is_ptcl=.false.)
                states = build%spproj%os_cls2D%get_all('state')
                call build%spproj%os_cls3D%set_all('state',states)
                deallocate(states)
                call build%spproj%add_frcs2os_out( trim(FRCS_FILE), 'frc2D')
                call build%spproj%add_cavgs2os_out(trim(params%refs), build%spproj%get_smpd(), imgkind='cavg')
                ! multiple fields updated, do a full write
                call build%spproj%write(params%projfile)
            case('ptcl3D')
                call build%eulspace%new(params%nspace, is_ptcl=.false.)
                call build%pgrpsyms%build_refspiral(build%eulspace)
                do icls = 1,params%ncls
                    call build%spproj%os_cls3D%set_euler(icls, build%eulspace%get_euler(icls))
                enddo
                call build%spproj%write_segment_inside('cls3D', params%projfile)
                if( cline%defined('outfile') )then
                    call build%spproj%os_cls3D%write(params%outfile)
                else
                    call build%spproj%os_cls3D%write('cls3D_oris.txt')
                endif
                call build%spproj%add_frcs2os_out( trim(FRCS_FILE), 'frc3D')
                call build%spproj%add_cavgs2os_out(trim(params%refs), build%spproj%get_smpd(), imgkind='cavg3D')
                call build%spproj%write_segment_inside('out', params%projfile)
            case DEFAULT
                THROW_HARD('Unsupported ORITYPE: '//trim(params%oritype))
        end select
        ! end gracefully
        call starproj%kill
        call build%spproj%kill
        call build%kill_general_tbox
        call build%kill_strategy2D_tbox
        call simple_end('**** SIMPLE_CAVGASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        call simple_touch('CAVGASSEMBLE_FINISHED', errmsg='In: commander_rec :: eo_cavgassemble ')
    end subroutine exec_cavgassemble

    subroutine exec_prob_tab2D_distr( self, cline )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_eul_prob_tab2D,      only: eul_prob_tab2D
        use simple_strategy2D3D_common, only: sample_ptcls4update
        class(prob_tab2D_commander_distr), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        integer,            allocatable :: pinds(:)
        character(len=:),   allocatable :: fname
        type(builder)                   :: build
        type(parameters)                :: params
        type(prob_tab2D_commander)      :: xprob_tab2D
        type(eul_prob_tab2D)            :: eulprob
        type(cmdline)                   :: cline_prob_tab2D
        type(qsys_env)                  :: qenv
        type(chash)                     :: job_descr
        integer :: nptcls, ipart
        if( associated(build_glob) )then
            if( .not.associated(params_glob) )then
                THROW_HARD('Builder & parameters must be associated for shared memory execution!')
            endif
        else
            call cline%set('mkdir',   'no')
            call cline%set('stream',  'no')
            call cline%set('oritype', 'ptcl2D')
            call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        endif
        if( params_glob%startit == 1 )then
            if( cline%defined('updatecnt_ini') )then
                call build_glob%spproj_field%set_nonzero_updatecnt(params%updatecnt_ini)
                call build_glob%spproj_field%clean_entry('sampled')
            else
                call build_glob%spproj_field%clean_entry('updatecnt', 'sampled')
            endif
        endif
        ! sampled incremented
        call sample_ptcls4update([1,params_glob%nptcls], .true., nptcls, pinds)
        ! communicate to project file
        call build_glob%spproj%write_segment_inside(params_glob%oritype, params_glob%projfile)
        ! more prep
        call eulprob%new(pinds)
        ! generating all scores
        cline_prob_tab2D = cline
        call cline_prob_tab2D%set('prg', 'prob_tab2D' )
        ! execution
        ! if( .not.cline_prob_tab2D%defined('nparts') )then
        !     call xprob_tab2D%execute_safe(cline_prob_tab)
        ! else
            ! setup the environment for distributed execution
            call qenv%new(params_glob%nparts, nptcls=params_glob%nptcls)
            call cline_prob_tab2D%gen_job_descr(job_descr)
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        ! endif
        ! reading scores from all parts
        do ipart = 1, params_glob%nparts
            fname = trim(DIST_FBODY)//int2str_pad(ipart,params_glob%numlen)//'.dat'
            call eulprob%read_table_to_glob(fname)
        enddo
        ! perform assignment
        select case(trim(params%refine))
        case('prob')
            if( params%which_iter == 1 )then
                call eulprob%assign_cls_greedy
            else
                call eulprob%normalize_table
                call eulprob%assign_cls_stoch(build_glob%spproj_field)
            endif
        case('prob_greedy')
            call eulprob%assign_cls_greedy
        end select
        ! write
        fname = trim(ASSIGNMENT_FBODY)//'.dat'
        call eulprob%write_assignment(fname)
        ! cleanup
        call eulprob%kill
        call cline_prob_tab2D%kill
        call qenv%kill
        call job_descr%kill
        call qsys_job_finished('simple_commander_cluster2D :: exec_prob_tab2D_distr')
        call qsys_cleanup
        call simple_end('**** SIMPLE_PROB_TAB2D_DISTR NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab2D_distr

    subroutine exec_prob_tab2D( self, cline )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_strategy2D3D_common
        use simple_classaverager
        use simple_strategy2D_matcher
        use simple_polarft_corrcalc,  only: polarft_corrcalc
        use simple_eul_prob_tab2D,    only: eul_prob_tab2D
        class(prob_tab2D_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        integer,          allocatable :: pinds(:)
        character(len=:), allocatable :: fname
        type(polarft_corrcalc)        :: pftcc
        type(builder)                 :: build
        type(parameters)              :: params
        type(eul_prob_tab2D)          :: eulprob
        type(euclid_sigma2)           :: eucl_sigma_glob
        real    :: frac_srch_space
        integer :: nptcls
        logical :: l_ctf
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! The policy here ought to be that nothing is done with regards to sampling other than reproducing
        ! what was generated in the driver (prob_align, below). Sampling is delegated to prob_align (below)
        ! and merely reproduced here
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab requires prior particle sampling (in exec_prob_align)')
        endif
        ! Resolution range
        frac_srch_space = build%spproj_field%get_avg('frac')
        if( file_exists(params_glob%frcs) ) call build%clsfrcs%read(params_glob%frcs)
        call set_bp_range2D( cline, params%which_iter, frac_srch_space )
        ! Read references
        call cavger_new(pinds)
        call cavger_read(params%refs, 'merged')
        if( file_exists(params%refs_even) )then
            call cavger_read(params%refs_even, 'even')
        else
            call cavger_read(params%refs, 'even')
        endif
        if( file_exists(params%refs_odd) )then
            call cavger_read(params%refs_odd, 'odd')
        else
            call cavger_read(params%refs, 'odd')
        endif
        ! init scorer & prep references
        call preppftcc4align2D(pftcc, nptcls, params%which_iter)
        ! prep particles
        l_ctf = build%spproj%get_ctfflag('ptcl2D',iptcl=params%fromp).ne.'no'
        call prep_batch_particles2D(nptcls)
        call build_batch_particles2D(pftcc, nptcls, pinds, l_ctf)
        ! init prob table
        call eulprob%new(pinds)
        fname = trim(DIST_FBODY)//int2str_pad(params%part,params%numlen)//'.dat'
        ! algorithm
        select case(trim(params%refine))
        case('prob')
            call eulprob%fill_table_stoch_inpl
        case('prob_greedy')
            call eulprob%fill_table_greedy_inpl
        end select
        ! write
        call eulprob%write_table(fname)
        ! clean & end
        if( associated(eucl_sigma2_glob) ) call eucl_sigma2_glob%kill
        call clean_batch_particles2D
        call eulprob%kill
        call pftcc%kill
        call build%kill_general_tbox
        call qsys_job_finished('simple_commander_cluster2D :: exec_prob_tab')
        call simple_end('**** SIMPLE_PROB_TAB2D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab2D

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
                call clsdoc_ranked%set(icls, 'class',     order(icls))
                call clsdoc_ranked%set(icls, 'rank',      icls)
                call clsdoc_ranked%set(icls, 'pop',       spproj%os_cls2D%get(order(icls),  'pop'))
                call clsdoc_ranked%set(icls, 'res',       spproj%os_cls2D%get(order(icls),  'res'))
                call clsdoc_ranked%set(icls, 'corr',      spproj%os_cls2D%get(order(icls), 'corr'))
                call clsdoc_ranked%set(icls, 'w',         spproj%os_cls2D%get(order(icls),    'w'))
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
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_polarizer,        only: polarizer
        use simple_class_frcs,       only: class_frcs
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_aff_prop,         only: aff_prop
        use simple_pftcc_shsrch_fm
        class(cluster_cavgs_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(pftcc_shsrch_fm), allocatable :: fm_correlators(:)
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(class_frcs)              :: clsfrcs
        type(image)                   :: img_msk
        type(polarft_corrcalc)        :: pftcc
        type(aff_prop)                :: aprop
        type(image),      allocatable :: cavg_imgs(:), ccimgs(:,:)
        type(polarizer)               :: polartransform
        character(len=:), allocatable :: cavgsstk, cavgsstk_shifted, classname, frcs_fname
        real,             allocatable :: states(:), frc(:), filter(:), clspops(:), clsres(:)
        real,             allocatable :: corrs(:), corrmat(:,:), corrs_top_ranking(:)
        logical,          allocatable :: l_msk(:,:,:), mask_top_ranking(:), mask_otsu(:), mask_icls(:)
        integer,          allocatable :: order(:), nloc(:), centers(:), labels(:), cntarr(:), clsinds(:), pops(:)
        integer,          allocatable :: labels_mapped(:), classmapping(:)
        logical,          parameter   :: DEBUG = .true.
        integer :: ldim(3), loc(1), ncls, n, ncls_sel, i, j, icls, cnt, filtsz, pop1, pop2, nsel
        integer :: ithr, ncls_aff_prop, icen, jcen
        real    :: smpd, simsum, cmin, cmax, pref, corr_icls, cc,ccm
        logical :: l_apply_optlp, use_shifted
        ! defaults
        call cline%set('oritype', 'cls2D')
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'cc')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('trs')     ) call cline%set('trs',       10.)
        if( .not. cline%defined('kweight') ) call cline%set('kweight', 'all')
        ! master parameters
        call params%new(cline)
        ! get class average stack
        l_apply_optlp = .false.
        if( cline%defined('stk') )then
            cavgsstk      = trim(params%stk)
            use_shifted   = .false.
            l_apply_optlp = .false.
            call find_ldim_nptcls(cavgsstk, ldim, ncls, smpd=smpd)
            write(logfhandle,'(A,I3)') '# classes in stack ', ncls
            ldim(3) = 1
            allocate(states(ncls),source=1.)
            filtsz = fdim(params%box)-1
        else
            call spproj%read_segment(params%oritype, params%projfile)
            call spproj%get_cavgs_stk(cavgsstk_shifted, ncls, smpd, imgkind='cavg_shifted', fail=.false.)
            use_shifted = .true.
            if( cavgsstk_shifted .eq. NIL ) use_shifted = .false.
            call spproj%get_cavgs_stk(cavgsstk, ncls, smpd)
            call find_ldim_nptcls(cavgsstk, ldim, n)
            write(logfhandle,'(A,I3)') '# classes in stack ', n
            ldim(3) = 1
            if( n /= ncls ) THROW_HARD('Inconsistent # classes in project file vs cavgs stack; exec_cluster_cavgs')
            ! threshold based on states/population/resolution
            states = spproj%os_cls2D%get_all('state')
            if( .not. DEBUG )then
                clspops = spproj%os_cls2D%get_all('pop')
                clsres  = spproj%os_cls2D%get_all('res')
                where( clsres  >= params%lpthres     ) states = 0.
                where( clspops <  real(MINCLSPOPLIM) ) states = 0.
            endif
            ! get FRCs
            call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
            if( file_exists(frcs_fname) )then
                call clsfrcs%read(frcs_fname)
                l_apply_optlp = .true.
            endif
            filtsz = clsfrcs%get_filtsz()
        endif
        ! ensure correct smpd/box in params class
        params%smpd = smpd
        params%box  = ldim(1)
        params%msk  = min(real(params%box/2)-COSMSKHALFWIDTH-1., 0.5*params%mskdiam /params%smpd)
        ! find out how many selected class averages initially
        ncls_sel = count(states > 0.5)
        write(logfhandle,'(A,I3)') '# classes left after standard rejection ', ncls_sel
        ! keep track of the original class indices
        clsinds = pack((/(i,i=1,ncls)/), mask=states > 0.5)
        ! create the stuff needed in the loop
        allocate(cavg_imgs(ncls_sel), frc(filtsz), filter(filtsz))
        ! prep mask
        call img_msk%new([params%box,params%box,1], params%smpd)
        img_msk = 1.
        call img_msk%mask(params%msk, 'hard')
        l_msk = img_msk%bin2logical()
        call img_msk%kill
        ! polarizer
        call polartransform%new([params%box,params%box,1], params%smpd)
        write(logfhandle,'(A)') '>>> PREPARING CLASS AVERAGES'
        ! read images
        do i = 1, ncls_sel
            call cavg_imgs(i)%new(ldim, smpd, wthreads=.false.)
            j = clsinds(i)
            if( use_shifted )then
                call cavg_imgs(i)%read(cavgsstk_shifted, j)
            else
                call cavg_imgs(i)%read(cavgsstk, j)
            endif
        enddo
        !$omp parallel do default(shared) private(i,j,frc,filter) schedule(static) proc_bind(close)
        do i = 1, ncls_sel
            j = clsinds(i)
            ! FRC-based filter
            if( l_apply_optlp )then
                call clsfrcs%frc_getter(j, frc)
                if( any(frc > 0.143) )then
                    call fsc2optlp_sub(clsfrcs%get_filtsz(), frc, filter)
                    where( filter > TINY ) filter = sqrt(filter) ! because the filter is applied to the average not the even or odd
                    call cavg_imgs(i)%fft()
                    call cavg_imgs(i)%apply_filter_serial(filter)
                    call cavg_imgs(i)%ifft()
                endif
            endif
            ! normalization
            call cavg_imgs(i)%norm_within(l_msk)
            ! mask
            call cavg_imgs(i)%mask(params%msk, 'soft', backgr=0.)
        end do
        !$omp end parallel do
        if( DEBUG )then
            do i = 1, ncls_sel
                call cavg_imgs(i)%write('cavgs_prepped.mrc', i)
            enddo
        endif
        ! resolution limits
        params%kfromto(1) = max(2, calc_fourier_index(params%hp, params%box, params%smpd))
        params%kfromto(2) =        calc_fourier_index(params%lp, params%box, params%smpd)
        ! Fourier-Mellin transform
        write(logfhandle,'(A)') '>>> CALCULATING CORRELATION OF MAGNITUDES MATRIX'
        allocate(corrmat(ncls_sel,ncls_sel), source=-1.)
        ! polarizer, pftcc
        params%sh_inv = 'yes'
        call pftcc%new(ncls_sel, [1,ncls_sel], params%kfromto)
        call polartransform%init_polarizer(pftcc, params%alpha)
        !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
        do icls = 1, ncls_sel
            call cavg_imgs(icls)%fft()
            ! particles & even refs are the same, odd refs are even mirrored
            call polartransform%polarize(pftcc, cavg_imgs(icls), icls, isptcl=.false., iseven=.true.)
            call pftcc%cp_even_ref2ptcl(icls,icls)
            call pftcc%mirror_pft(pftcc%pfts_refs_even(:,:,icls), pftcc%pfts_refs_odd(:,:,icls))            
        end do
        !$omp end parallel do
        call pftcc%memoize_refs
        call pftcc%memoize_ptcls
        allocate(fm_correlators(nthr_glob),ccimgs(nthr_glob,2))
        do i = 1,nthr_glob
            call fm_correlators(i)%new(params%trs,1.,opt_angle=.false.)
            call ccimgs(i,1)%new(ldim, smpd, wthreads=.false.)
            call ccimgs(i,2)%new(ldim, smpd, wthreads=.false.)
        enddo
        !$omp parallel do default(shared) private(i,j,cc,ccm,ithr)&
        !$omp schedule(dynamic) proc_bind(close)
        do i = 1, ncls_sel - 1
            ithr = omp_get_thread_num()+1
            corrmat(i,i) = 1.
            do j = i, ncls_sel
                ! correlation of reference to particle
                call pftcc%set_eo(i,.true.)
                call fm_correlators(ithr)%calc_phasecorr(j, i, cavg_imgs(j), cavg_imgs(i),&
                    &ccimgs(ithr,1), ccimgs(ithr,2), cc)
                ! correlation of mirrored reference to particle
                ! call pftcc%set_eo(i,.false.)
                ! call fm_correlators(ithr)%calc_phasecorr(j, i, cavg_imgs(j), cavg_imgs(i),&
                !     &ccimgs(ithr,1), ccimgs(ithr,2), ccm, mirror=.true.)
                ! cc = max(cc,ccm)
                corrmat(i,j) = cc
                corrmat(j,i) = cc
            enddo
        enddo
        !$omp end parallel do
        corrmat(ncls_sel,ncls_sel) = 1.
        do i = 1,nthr_glob
            call fm_correlators(i)%kill
            call ccimgs(i,1)%kill
            call ccimgs(i,2)%kill
        enddo
        deallocate(fm_correlators,ccimgs)
        if( trim(params%bin_cls).eq.'yes' )then
            write(logfhandle,'(A)') '>>> GOOD/BAD CLASSIFICATION AND RANKING OF CLASS AVERAGES'
            nsel = ceiling(0.1 * real(ncls_sel))
            allocate(mask_top_ranking(ncls_sel),corrs_top_ranking(ncls_sel),nloc(nsel))
            ! average the nsel best correlations (excluding self) to create a scoring function for garbage removal
            do i = 1, ncls_sel
                mask_top_ranking    = .true.
                mask_top_ranking(i) = .false. ! to remove the diagonal element
                corrs = pack(corrmat(i,:), mask_top_ranking)
                nloc  = maxnloc(corrs, nsel)
                corrs_top_ranking(i) = sum(corrs(nloc(1:nsel))) / real(nsel)
            end do
            ! use Otsu's algorithm to remove the junk
            allocate(mask_otsu(ncls_sel))
            call otsu(ncls_sel, corrs_top_ranking, mask_otsu)
            pop1 = count(      mask_otsu)
            pop2 = count(.not. mask_otsu)
            write(logfhandle,*) 'average corr cluster 1: ', sum(corrs_top_ranking, mask=      mask_otsu) / real(pop1), ' pop ', pop1
            write(logfhandle,*) 'average corr cluster 2: ', sum(corrs_top_ranking, mask=.not. mask_otsu) / real(pop2), ' pop ', pop2
            pop1 = 0
            pop2 = 0
            do i = 1, ncls_sel
                call cavg_imgs(i)%ifft()
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
        else
            write(logfhandle,'(A)') '>>> CLUSTERING CLASS AVERAGES WITH AFFINITY PROPAGATION'
            ! calculate a preference that generates a small number of clusters
            call analyze_smat(corrmat, .false., cmin, cmax)
            pref = cmin - (cmax - cmin)
            call aprop%new(ncls_sel, corrmat, pref=pref)
            call aprop%propagate(centers, labels, simsum)
            call aprop%kill
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
                    corrmat(i,i) = 1.
                end do
                ! make a per cluster mapping
                do icls = 1, ncls_aff_prop
                    do i = 1,params%ncls
                        corrs(i) = corrmat(nloc(i),icls)
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
            if( DEBUG )then
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
                cntarr = 0
                do icls = 1, ncls_aff_prop
                    do i=1,ncls_sel
                        if( labels(i) == icls )then
                            classname = 'class'//int2str_pad(icls,3)//trim(STK_EXT)
                            cntarr(labels(i)) = cntarr(labels(i)) + 1
                            call cavg_imgs(i)%ifft
                            call cavg_imgs(i)%write(classname, cntarr(labels(i)))
                        endif
                    end do
                end do
            endif
            ! calculate average correlation of the AP clusters
            if( allocated(mask_icls) ) deallocate(mask_icls)
            allocate( mask_icls(ncls_sel) )
            do icls = 1, ncls_aff_prop
                ! mask out the cluster
                where( labels == icls )
                    mask_icls = .true.
                elsewhere
                    mask_icls = .false.
                end where
                ! calculate the average correlation of the cluster
                if( count(mask_icls) == 1 )then
                    corr_icls = 1.
                else
                    corr_icls = 0.
                    cnt       = 0
                    do i = 1, ncls_sel - 1
                        do j = i + 1, ncls_sel
                            if( mask_icls(i) .and. mask_icls(j) )then
                                corr_icls = corr_icls + corrmat(i,j)
                                cnt       = cnt + 1
                            endif
                        end do
                    end do
                    corr_icls = corr_icls / real(cnt)
                endif
                write(logfhandle,*) 'average corr AP cluster '//int2str_pad(icls,3)//': ', corr_icls, ' pop ', count(mask_icls)
            end do
            ! print the correlations between the cluster centers
            write(logfhandle,'(A)') '>>> CORRELATIONS BETWEEN CLUSTER CENTERS'
            do i = 1,ncls_aff_prop - 1
                icen = centers(i)
                do j = i + 1,ncls_aff_prop
                    jcen = centers(j)
                    write(logfhandle,*) 'corr btw cluster  '//int2str_pad(i,3)//' & '//int2str_pad(j,3)//': ', corrmat(icen,jcen)
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
        deallocate(cavg_imgs)
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

    subroutine exec_ppca_denoise_classes( self, cline )
        use simple_imgproc,       only: make_pcavecs
        use simple_image,         only: image
        use simple_classaverager, only: transform_ptcls
        use simple_pca,           only: pca
        use simple_pca_svd,       only: pca_svd
        use simple_kpca_svd,      only: kpca_svd
        use simple_ppca_inmem,    only: ppca_inmem
        class(ppca_denoise_classes_commander), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        integer,          parameter   :: MAXPCAITS = 15
        class(pca),       pointer     :: pca_ptr  => null()
        type(parameters)              :: params
        type(builder)                 :: build
        type(image),      allocatable :: imgs(:), imgs_ori(:)
        type(image)                   :: cavg
        type(oris)                    :: os
        type(sp_project), target      :: spproj
        character(len=:), allocatable :: label, fname, fname_denoised, fname_cavgs, fname_cavgs_denoised, fname_oris
        integer,          allocatable :: cls_inds(:), pinds(:), cls_pops(:)
        real,             allocatable :: avg(:), avg_pix(:), pcavecs(:,:), tmpvec(:)
        real             :: std
        integer          :: npix, i, j, ncls, nptcls, cnt1, cnt2, neigs
        logical          :: l_phflip, l_transp_pca, l_pre_norm ! pixel-wise learning
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        if( .not. cline%defined('neigs')   ) call cline%set('neigs',    4)
        call build%init_params_and_build_general_tbox(cline, params, do3d=(trim(params%oritype) .eq. 'ptcl3D'))
        call spproj%read(params%projfile)
        select case(trim(params%oritype))
            case('ptcl2D')
                label = 'class'
            case('ptcl3D')
                label = 'proj'
                call build%spproj_field%proj2class
            case DEFAULT
                THROW_HARD('ORITYPE not supported!')
        end select
        l_transp_pca = (trim(params%transp_pca) .eq. 'yes')
        l_pre_norm   = (trim(params%pre_norm)   .eq. 'yes')
        l_phflip     = .false.
        select case( spproj%get_ctfflag_type(params%oritype) )
            case(CTFFLAG_NO)
                THROW_WARN('No CTF information could be found, phase flipping is deactivated')
            case(CTFFLAG_FLIP)
                THROW_WARN('Images have already been phase-flipped, phase flipping is deactivated')
            case(CTFFLAG_YES)
                l_phflip = .true.
            case DEFAULT
                THROW_HARD('UNSUPPORTED CTF FLAG')
        end select
        cls_inds = build%spproj_field%get_label_inds(label)
        if( cline%defined('class') .and. cline%defined('ncls') )then
            THROW_HARD('EITHER class OR cls CAN BE DEFINED')
        endif
        if( cline%defined('class') )then
            cls_inds = pack(cls_inds, mask=(cls_inds == params%class))
        endif
        ncls = size(cls_inds)
        if( cline%defined('ncls') )then
            ncls     = params%ncls
            cls_inds = cls_inds(1:ncls)
        endif
        allocate(cls_pops(ncls), source=0)
        do i = 1, ncls
            call build%spproj_field%get_pinds(cls_inds(i), label, pinds)
            if( allocated(pinds) )then
                cls_pops(i) = size(pinds)
                nptcls = nptcls + cls_pops(i)
                deallocate(pinds)
            endif
        end do
        cls_inds = pack(cls_inds, mask=cls_pops > 2)
        nptcls   = sum(cls_pops,  mask=cls_pops > 2)
        ncls     = size(cls_inds)
        call os%new(nptcls, is_ptcl=.true.)
        fname                = 'ptcls.mrcs'
        fname_denoised       = 'ptcls_denoised.mrcs'
        fname_cavgs          = 'cavgs.mrcs'
        fname_cavgs_denoised = 'cavgs_denoised.mrcs'
        fname_oris           = 'oris_denoised.txt'
        cnt1 = 0
        cnt2 = 0
        ! pca allocation
        select case(trim(params_glob%pca_mode))
            case('ppca')
                allocate(ppca_inmem :: pca_ptr)
            case('pca_svd')
                allocate(pca_svd    :: pca_ptr)
            case('kpca')
                allocate(kpca_svd   :: pca_ptr)
        end select
        do i = 1, ncls
            call progress_gfortran(i,ncls)
            if( trim(params%pca_img_ori) .eq. 'yes' )then
                call transform_ptcls(spproj, params%oritype, cls_inds(i), imgs, pinds, phflip=l_phflip, cavg=cavg, imgs_ori=imgs_ori)
                do j = 1, size(imgs)
                    call imgs(j)%copy_fast(imgs_ori(j))
                enddo
            else
                call transform_ptcls(spproj, params%oritype, cls_inds(i), imgs, pinds, phflip=l_phflip, cavg=cavg)
            endif
            nptcls = size(imgs)
            neigs  = params%neigs
            if( neigs >= nptcls )then
                THROW_WARN('neigs is greater than the number of particles within this class. All eigens are used now!')
                neigs = nptcls - 1
            endif
            if( l_pre_norm )then
                do j = 1, nptcls
                    call imgs(j)%norm
                end do
            endif
            do j = 1, nptcls
                cnt1 = cnt1 + 1
                call imgs(j)%write(fname, cnt1)
            end do
            call cavg%write(fname_cavgs, i)
            ! performs ppca
            if( trim(params%projstats).eq.'yes' )then
                call make_pcavecs(imgs, npix, avg, pcavecs, transp=l_transp_pca, avg_pix=avg_pix)
            else
                call make_pcavecs(imgs, npix, avg, pcavecs, transp=l_transp_pca)
            endif
            if( allocated(tmpvec) ) deallocate(tmpvec)
            if( l_transp_pca )then
                call pca_ptr%new(npix, nptcls, neigs)
                call pca_ptr%master(pcavecs, MAXPCAITS)
                allocate(tmpvec(nptcls))
                !$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
                do j = 1, npix
                    call pca_ptr%generate(j, avg, tmpvec)
                    pcavecs(:,j) = tmpvec
                end do
                !$omp end parallel do
                pcavecs = transpose(pcavecs)
            else
                call pca_ptr%new(nptcls, npix, neigs)
                call pca_ptr%master(pcavecs, MAXPCAITS)
                allocate(tmpvec(npix))
                !$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
                do j = 1, nptcls
                    call pca_ptr%generate(j, avg, tmpvec)
                    pcavecs(:,j) = tmpvec
                end do
                !$omp end parallel do
            endif
            if( trim(params%projstats).eq.'yes' )then
                call cavg%unserialize(avg_pix)
                call cavg%write('cavgs_unserialized.mrcs', i)
                !$omp parallel do private(j,std) default(shared) proc_bind(close) schedule(static)
                do j = 1,nptcls
                    std = sqrt(sum((pcavecs(:,j)-avg_pix)**2) / real(npix))
                enddo
                !$omp end parallel do
            endif
            ! output
            call cavg%zero_and_unflag_ft
            if( trim(params%pca_img_ori) .eq. 'yes' )then
                do j = 1, nptcls
                    cnt2 = cnt2 + 1
                    call imgs_ori(j)%unserialize(pcavecs(:,j))
                    call os%transfer_ori(cnt2, build%spproj_field, pinds(j))
                    call imgs_ori(j)%write(fname_denoised, cnt2)
                end do
                call transform_ptcls(spproj, params%oritype, cls_inds(i), imgs, pinds, phflip=l_phflip, cavg=cavg, imgs_ori=imgs_ori, just_transf=.true.)
            else
                do j = 1, nptcls
                    cnt2 = cnt2 + 1
                    call imgs(j)%unserialize(pcavecs(:,j))
                    call cavg%add(imgs(j))
                    call os%transfer_ori(cnt2, build%spproj_field, pinds(j))
                    call imgs(j)%write(fname_denoised, cnt2)
                    call imgs(j)%kill
                end do
                call cavg%div(real(nptcls))
            endif
            call cavg%write(fname_cavgs_denoised, i)
        end do
        call os%zero_inpl
        call os%write(fname_oris)
        if( trim(params%projstats).eq.'yes' ) call build%spproj_field%write('ptcl_field.txt')
        ! cleanup
        deallocate(imgs)
        call build%kill_general_tbox
        call os%kill
        ! end gracefully
        call simple_end('**** SIMPLE_PPCA_DENOISE_CLASSES NORMAL STOP ****')
    end subroutine exec_ppca_denoise_classes

    subroutine exec_partition_cavgs( self, cline )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_polarizer,           only: polarizer
        use simple_class_frcs,          only: class_frcs
        use simple_polarft_corrcalc,    only: polarft_corrcalc
        use simple_aff_prop,            only: aff_prop
        use simple_spectral_clustering, only: spec_clust
        use simple_pftcc_shsrch_fm,     only: pftcc_shsrch_fm
        class(partition_cavgs_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        character(len=STDLEN),   parameter :: SIMMAT_FNAME = 'simmat.bin'
        character(len=STDLEN),   parameter :: LABELS_FNAME = 'labels.txt'
        real,                    parameter :: LP_DEFAULT   = 6.0 ! Angs
        logical,                 parameter :: DEBUG = .true.
        type(pftcc_shsrch_fm), allocatable :: fm_correlators(:)
        type(image),           allocatable :: cavg_imgs(:), ccimgs(:,:)
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(class_frcs)              :: clsfrcs
        type(image)                   :: img_msk
        type(polarft_corrcalc)        :: pftcc
        type(aff_prop)                :: aprop
        type(spec_clust)              :: specclust
        type(polarizer)               :: polartransform
        character(len=:), allocatable :: cavgsstk, frcs_fname, fmt, header1, header2, header3
        real,             allocatable :: states(:), frc(:), filter(:), clspops(:), clsres(:)
        real,             allocatable :: corrmat(:,:), S(:,:), rotmat(:,:), xoffset(:,:), yoffset(:,:)
        integer,          allocatable :: centers(:), labels(:), clsinds(:), multi_labels(:,:),tmp(:)
        logical,          allocatable :: l_msk(:,:,:)
        character(len=8)              :: tmpstr
        real    :: offset(2),minmax(2),smpd,simsum,simmin,simmax,simmed,pref,cc,ccm,dunn,dunnmax
        integer :: ldim(3), n, ncls_sel, i, j, k, icls, filtsz, ind, ithr, funit, io_stat
        logical :: l_apply_optlp, l_mirr
        ! defaults
        call cline%set('oritype', 'out')
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'cc')
        call cline%set('sh_inv',  'yes')
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',     'yes')
        if( .not. cline%defined('kweight')   ) call cline%set('kweight',   'all')
        if( .not. cline%defined('mirr')      ) call cline%set('mirr',      'no')
        if( .not. cline%defined('algorithm') ) call cline%set('algorithm', 'affprop')
        if( .not. cline%defined('nsearch')   ) call cline%set('nsearch',   7)
        ! parse parameters
        call params%new(cline)
        l_mirr = trim(params%mirr).eq.'yes'
        select case(trim(params%algorithm))
        case('affprop','spc') ! supported methods
        case DEFAULT
            THROW_HARD('Unsupported clustering algorithm')
        end select
        ! fetch cavgs stack
        call spproj%read(params%projfile)
        call spproj%get_cavgs_stk(cavgsstk, params%ncls, smpd)
        if( trim(params%mkdir).eq.'no' )then
            ! to deal with relative path, the project will be updated in place
            if( .not.file_exists(cavgsstk) )then
                if( cavgsstk(1:3) == '../' ) cavgsstk = cavgsstk(4:len_trim(cavgsstk))
            endif
        endif
        params%stk = trim(cavgsstk)
        call find_ldim_nptcls(params%stk, ldim, n)
        ldim(3) = 1
        if( n /= params%ncls ) THROW_HARD('Inconsistent # classes in project file vs cavgs stack; partition_cavgs')
        ! threshold based on states/population/resolution
        states = spproj%os_cls2D%get_all('state')
        if( spproj%os_cls2D%isthere('pop') )then
            clspops = spproj%os_cls2D%get_all('pop')
            where( clspops <  real(MINCLSPOPLIM) ) states = 0.
        endif
        if( spproj%os_cls2D%isthere('res') )then
            if( cline%defined('lpthres') )then
                clsres  = spproj%os_cls2D%get_all('res')
                where( clsres  >= params%lpthres ) states = 0.
            endif
        endif
        ! number of classes to partition
        ncls_sel = count(states>0.5)
        ! keep track of the original class indices
        clsinds = pack((/(i,i=1,params%ncls)/), mask=states > 0.5)
        ! get FRCs
        call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
        if( trim(params%mkdir).eq.'no' )then
            ! to deal with relative path, the project will be updated in place
            if( .not.file_exists(frcs_fname) )then
                if( frcs_fname(1:3) == '../' ) frcs_fname = frcs_fname(4:len_trim(frcs_fname))
            endif
        endif
        ! ensure correct smpd/box in params class
        params%smpd = smpd
        params%box  = ldim(1)
        if( cline%defined('mskdiam') )then
            params%msk = min(real(params%box/2)-COSMSKHALFWIDTH-1., 0.5*params%mskdiam /params%smpd)
        else
            params%msk = 0.9*real(params%box/2)-COSMSKHALFWIDTH-1.
        endif
        write(logfhandle,'(A)') '>>> PREPARING CLASS AVERAGES'
        ! read images
        allocate(cavg_imgs(ncls_sel))
        do i = 1, ncls_sel
            call cavg_imgs(i)%new(ldim, smpd, wthreads=.false.)
            j = clsinds(i)
            call cavg_imgs(i)%read(params%stk, j)
            minmax = cavg_imgs(i)%minmax()
            if( is_equal(minmax(1),minmax(2)) )then
                ! empty class
                clsinds(i) = 0
                states(j)  = 0
            endif
        enddo
        if( ncls_sel /= count(clsinds>0) )then
            ! updating for empty classes
            cavg_imgs(:) = pack(cavg_imgs(:), mask=clsinds(:)>0)
            ncls_sel     = count(clsinds>0)
            tmp          = pack(clsinds(:), mask=clsinds(:)>0)
            clsinds      = tmp(:)
            deallocate(tmp)
        endif
        ! resolution limits
        l_apply_optlp = .false.
        if( file_exists(frcs_fname) )then
            call clsfrcs%read(frcs_fname)
            ! whether to filter
            l_apply_optlp = .true.
            ! resolution limit
            if( .not.cline%defined('lp') )then
                ! FRC=0.5 criterion
                params%lp = clsfrcs%estimate_lp_for_align(crit0143=.false.)
            endif
            filtsz = clsfrcs%get_filtsz()
            allocate(frc(filtsz), filter(filtsz))
        else
            THROW_WARN('LP or frcs.bin must be present to determine resolution, default value used')
            params%lp = LP_DEFAULT
        endif
        write(logfhandle,'(A,F6.1)')'>>> RESOLUTION LIMIT: ',params%lp
        ! prep mask
        call img_msk%new([params%box,params%box,1], params%smpd)
        img_msk = 1.
        call img_msk%mask(params%msk, 'hard')
        l_msk = img_msk%bin2logical()
        call img_msk%kill
        ! prep images
        !$omp parallel do default(shared) private(i,j,frc,filter) schedule(static) proc_bind(close)
        do i = 1, ncls_sel
            j = clsinds(i)
            if( l_apply_optlp )then                                 ! FRC-based filter
                call clsfrcs%frc_getter(j, frc)
                if( any(frc > 0.143) )then
                    call fsc2optlp_sub(clsfrcs%get_filtsz(), frc, filter)
                    where( filter > TINY ) filter = sqrt(filter)
                    call cavg_imgs(i)%fft()
                    call cavg_imgs(i)%apply_filter_serial(filter)
                    call cavg_imgs(i)%ifft()
                endif
            endif
            call cavg_imgs(i)%norm_within(l_msk)                    ! normalization
            call cavg_imgs(i)%mask(params%msk, 'soft', backgr=0.)   ! mask
        end do
        !$omp end parallel do
        if( DEBUG )then
            do i = 1, ncls_sel
                call cavg_imgs(i)%write('cavgs_prepped.mrc', i)
            enddo
        endif
        ! Resolution limits
        params%kfromto(1) = max(2, calc_fourier_index(params%hp, params%box, params%smpd))
        params%kfromto(2) =        calc_fourier_index(params%lp, params%box, params%smpd)
        allocate(corrmat(ncls_sel,ncls_sel), source=-1.)
        ! FOURIER-MELLIN TRANSFORM
        write(logfhandle,'(A)') '>>> CALCULATING CORRELATION OF MAGNITUDES MATRIX'
        ! Shift boundaries
        if( .not. cline%defined('trs') ) params%trs = real(params%box)/6.
        ! pftcc init
        call pftcc%new(ncls_sel, [1,ncls_sel], params%kfromto)
        ! PFT transform
        call polartransform%new([params%box,params%box,1], params%smpd)
        call polartransform%init_polarizer(pftcc, params%alpha)
        !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
        do icls = 1, ncls_sel
            call cavg_imgs(icls)%fft()
            ! particles & even refs are the same, odd refs are even mirrored
            call polartransform%polarize(pftcc, cavg_imgs(icls), icls, isptcl=.false., iseven=.true.)
            call pftcc%cp_even_ref2ptcl(icls,icls)
            if( l_mirr ) call pftcc%mirror_pft(pftcc%pfts_refs_even(:,:,icls), pftcc%pfts_refs_odd(:,:,icls))
        end do
        !$omp end parallel do
        call pftcc%memoize_refs
        call pftcc%memoize_ptcls
        ! correlation matrix calculation
        allocate(fm_correlators(nthr_glob),ccimgs(nthr_glob,2),rotmat(ncls_sel,ncls_sel),&
            &xoffset(ncls_sel,ncls_sel),yoffset(ncls_sel,ncls_sel))
        do i = 1,nthr_glob
            call fm_correlators(i)%new(params%trs,1.,opt_angle=.false.)
            call ccimgs(i,1)%new(ldim, smpd, wthreads=.false.)
            call ccimgs(i,2)%new(ldim, smpd, wthreads=.false.)
        enddo
        rotmat = 0.
        !$omp parallel do default(shared) private(i,j,cc,ccm,ithr,offset)&
        !$omp schedule(dynamic) proc_bind(close)
        do i = 1, ncls_sel - 1
            ithr = omp_get_thread_num()+1
            corrmat(i,i) = 1.
            do j = i + 1, ncls_sel
                ! reference to particle
                call pftcc%set_eo(i,.true.)
                call fm_correlators(ithr)%calc_phasecorr(j, i, cavg_imgs(j), cavg_imgs(i),&
                    &ccimgs(ithr,1), ccimgs(ithr,2), cc, rotang=rotmat(i,j), shift=offset)
                rotmat(j,i)  = rotmat(i,j)
                xoffset(i,j) = offset(1); xoffset(j,i) = offset(1)
                yoffset(i,j) = offset(2); yoffset(j,i) = offset(2)
                ! mirrored reference to particle
                if( l_mirr )then
                    call pftcc%set_eo(i,.false.)
                    call fm_correlators(ithr)%calc_phasecorr(j, i, cavg_imgs(j), cavg_imgs(i),&
                    &ccimgs(ithr,1), ccimgs(ithr,2), ccm, mirror=.true.)
                    cc = max(cc,ccm)
                endif
                corrmat(i,j) = cc
                corrmat(j,i) = cc
            enddo
        enddo
        !$omp end parallel do
        corrmat(ncls_sel,ncls_sel) = 1.
        ! write similarity matrix
        call rmat2file(corrmat, SIMMAT_FNAME)
        ! tidy
        call pftcc%kill
        call polartransform%kill_polarizer
        call polartransform%kill
        do i = 1,nthr_glob
            call fm_correlators(i)%kill
            call ccimgs(i,1)%kill
            call ccimgs(i,2)%kill
        enddo
        deallocate(ccimgs,fm_correlators)
        ! CLUSTERING
        allocate(multi_labels(params%ncls,params%nsearch),source=0)
        ! correlation matrix preprocessing
        ! correlation to squared euclidian distance, values within [0.;4.]
        corrmat = 2.*(1.-corrmat)
        where( corrmat < 0. ) corrmat = 0.
        where( corrmat > 4. ) corrmat = 4.
        select case(trim(params%algorithm))
        case('affprop')
            ! Affinity Propagation: negative squared euclidian distance as similarity
            ! to circumvent possible correlation sign change. Values within [-4.;0]
            corrmat = -corrmat
        case('spc')
            ! Spectral Clustering: euclidian distance as similarity
            where( corrmat < 0. ) corrmat =  0.
            corrmat = sqrt(corrmat)
        end select
        ! for labellings output
        header1 = '# K     '
        select case(trim(params%algorithm))
        case('affprop')
            header2 = '# PREF  '
            header3 = '# SIM   '
            ! Preference bounds
            call analyze_smat(corrmat, .false., simmin, simmax)
            simmed = median(pack(corrmat,.true.))  ! taken as upper bound
            ! clustering nsearch times with preference in [min,median]
            write(logfhandle,'(A)') '>>> CLUSTERING CLASS AVERAGES WITH AFFINITY PROPAGATION'
            do i = 1,params%nsearch
                ! preference
                pref = simmin + (simmed-simmin)*real(i-1)/real(params%nsearch-1)
                ! clustering
                S = corrmat ! because corrmat is modified
                call aprop%new(ncls_sel, S, pref=pref)
                call aprop%propagate(centers, labels, simsum)
                call aprop%kill
                ! update labels & header
                write(tmpstr,'(I8)') maxval(labels)
                header1 = header1//tmpstr
                write(tmpstr,'(F8.3)') pref
                header2 = header2//tmpstr
                write(tmpstr,'(F8.3)') simsum
                header3 = header3//tmpstr
                multi_labels(clsinds(:),i) = labels(:)
                k = maxval(labels)
                write(*,'(A16,I3,2F9.3)') '>>> K,PREF,SIM: ',k,pref,simsum
                if( DEBUG ) call write_partition('part'//int2str_pad(i,3))
                ! cleanup
                deallocate(centers, labels,S)
            enddo
            ! update & write project with labelling with smallest preference
            call update_project(multi_labels(1,:))
        case('spc')
            dunnmax = -1.
            ind     = 0
            header2 = '# Dunn  '
            do i = 1,params%nsearch
                ! clustering
                k = i+1
                call specclust%new(ncls_sel,k,corrmat,algorithm='cpqr')
                call specclust%cluster
                call specclust%get_labels(labels)
                dunn = specclust%dunnindex()
                call specclust%kill
                if( dunn > dunnmax )then
                    dunnmax = dunn
                    ind     = i
                endif
                write(*,'(A12,I3,F9.3)') '>>> K,Dunn: ',k,dunn
                ! update labels & header
                write(tmpstr,'(I8)') k
                header1 = header1//tmpstr
                write(tmpstr,'(F8.3)') dunn
                header2 = header2//tmpstr
                multi_labels(clsinds(:),i) = labels(:)
                if( DEBUG )call write_partition('part'//int2str_pad(k,3))
            enddo
            ! update & write project with labelling with highest Dunn index
            call update_project(multi_labels(ind,:))
        end select
        ! write nsearch labellings
        call fopen(funit,LABELS_FNAME, 'replace', 'unknown', iostat=io_stat, form='formatted')
        call fileiochk("fopen failed "//trim(LABELS_FNAME),io_stat)
        write(funit,'(A2,I6,I8)') '# ',params%ncls,params%nsearch
        write(funit,'(A)') header1
        write(funit,'(A)') header2
        if( trim(params%algorithm)=='affprop' ) write(funit,'(A)') header3
        fmt = '('//int2str(params%nsearch+1)//'I8)'
        do i = 1,params%ncls
            write(funit,fmt) i, multi_labels(i,:)
        end do
        call fclose(funit)
        ! cleanup
        call spproj%kill
        call clsfrcs%kill
        if( allocated(cavg_imgs) )then
            do i = 1, ncls_sel
                call cavg_imgs(i)%kill
            enddo
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_PARTITION_CAVGS NORMAL STOP ****')
        contains

            ! write labels to project
            subroutine update_project( labels )
                integer :: labels(ncls_sel)
                integer :: i, icls
                ! setting all to zero for the case where rejection has occured
                do icls = 1,params%ncls
                    if( spproj%os_cls2D%get_state(icls) /= 0 )then
                        call spproj%os_cls2D%set(icls, 'cluster', 0)
                    endif
                enddo
                ! reporting label to cls2D
                do i = 1,ncls_sel
                    icls = clsinds(i)
                    call spproj%os_cls2D%set(icls, 'cluster', labels(i))
                enddo
                ! to particles
                call spproj%map_cls2D_flag_to_ptcls('cluster')
                call spproj%write(params%projfile)
            end subroutine update_project

            ! write the classes for debugging purpose
            subroutine write_partition( flag )
                character(len=*), intent(in) :: flag
                type(image)                   :: img1,img2
                character(len=:), allocatable :: classname
                integer,          allocatable :: cntarr(:)
                integer :: ncls_here, icls, i, ifirst
                logical :: first
                ncls_here = size(labels)
                allocate(cntarr(ncls_here), source=0)
                call img1%copy(cavg_imgs(1))
                do icls = 1, ncls_here
                    classname  = trim(flag)//'_class'//int2str_pad(icls,3)//trim(STK_EXT)
                    first = .true.
                    do i = 1,ncls_sel
                        if( labels(i) == icls )then
                            cntarr(labels(i)) = cntarr(labels(i)) + 1
                            call cavg_imgs(i)%ifft
                            if( first )then
                                ifirst = i
                                first  = .false.
                                call img2%copy(cavg_imgs(i))
                            else
                                call img1%copy_fast(cavg_imgs(i))
                                call img1%fft
                                call img1%shift2Dserial([xoffset(ifirst,i), yoffset(ifirst,i)])
                                call img1%ifft
                                call img1%rtsq(rotmat(ifirst,i), 0., 0.,img2)
                            endif
                            call img2%write(classname, cntarr(labels(i)))
                        endif
                    end do
                end do
                deallocate(cntarr)
                call img1%kill
                call img2%kill
            end subroutine write_partition

    end subroutine exec_partition_cavgs

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
