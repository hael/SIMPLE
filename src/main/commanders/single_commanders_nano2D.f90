!@descr: 2D analysis commanders used in SINGLE for nanoparticle processing
module single_commanders_nano2D
use simple_commander_module_api
use simple_commanders_cluster2D, only: commander_cluster2D
use simple_commanders_mkcavgs,   only: commander_make_cavgs_distr
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_center2D_nano
  contains
    procedure :: execute      => exec_center2D_nano
end type commander_center2D_nano

type, extends(commander_base) :: commander_cluster2D_nano
  contains
    procedure :: execute      => exec_cluster2D_nano
end type commander_cluster2D_nano

type, extends(commander_base) :: commander_analysis2D_nano
  contains
    procedure :: execute      => exec_analysis2D_nano
end type commander_analysis2D_nano

contains

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

end module single_commanders_nano2D
