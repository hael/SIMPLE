!@descr: solve3D_pyramid commander: multistate 3D by class-average assignment in one pass, stream p07's method on a 2D-classified project, with a hard exit at every failed step
!==============================================================================
! MODULE: simple_commanders_solve3D_pyramid
!
! PURPOSE:
!   Does offline, in one pass, what stream p07 does session-long (doc/policies/3D/solve3D_pyramid_policy.md):
!   the states' reference volumes are made once from class averages of the first nptcls_refs selected
!   particles, every particle is assigned to a state by matching sub-class averages against them, each
!   state's sub-class averages are refined alone by solve3D_cavgs (its own frame and resolution limits),
!   and one refine3D refine=cont iteration refines the poses with the states fixed. The method's rules, constants and
!   job command lines are p07's, shared through simple_state_assignment; none of p07's stream machinery
!   (watching, publications, retries, fallbacks, GUI). Every step runs in-process in a folder of the
!   execution directory (its commander on a child command line, then back), and every failed step is a hard
!   exit after its log and files are on disk, with the report so far.
!
! STEPS (folders of the execution directory):
!   refs_subset/                 the first set: whole stacks in arrival order (cut_whole_stacks), in a
!                                copy of the project holding its particles only (prune_particles)
!   refs_subset/solve2D/         solve2D into NCLS2D_FIRST classes
!   refs_subset/quality_selection/  the pool model's selection of them (use_model=yes)
!   refs_subset/cls_expansion/   cls_expansion into expansion_budget sub-classes
!   refs_subset/solve3D_cavgs/   solve3D_cavgs, nstates
!   state_refs/                  the check (check_state_refs): the references, their record, class averages
!   assign/quality_selection/    the pool model on the project's classes, logged only
!   assign/cls_expansion/        cls_expansion of every selected particle's classes
!   assign/                      assign_cavgs_states, the rejection, flags and mapping (map_assignment)
!   assign/state_<NN>/           solve3D_cavgs nstates=1 on state NN's sub-class averages, its volume docked
!                                onto the state's reference, its poses mapped onto the state's particles
!   refine3D/                    refine3D refine=cont maxits=1, the states fixed
!   The result is the execution directory's project, with a report (REPORT_FNAME) beside it.
!==============================================================================
module simple_commanders_solve3D_pyramid
use simple_defs,                only: logfhandle, STDLEN, CWD_GLOB
use simple_defs_fname,          only: METADATA_EXT, TASK_FINISHED
use simple_string_utils,        only: int2str_pad
use simple_error,               only: simple_exception
use simple_string,              only: string
use simple_fileio,              only: file_exists, file2rarr, fopen, fclose
use simple_syslib,              only: simple_abspath, simple_chdir, simple_getcwd, simple_mkdir, simple_touch
use simple_jiffys,              only: simple_end
use simple_math_ft,             only: get_resarr
use simple_estimate_ssnr,       only: get_resolution
use simple_chash,               only: chash
use simple_cmdline,             only: cmdline
use simple_parameters,          only: parameters
use simple_sp_project,          only: sp_project
use simple_commander_base,      only: commander_base
use simple_commanders_solve2D,  only: commander_solve2D
use simple_commanders_denoise,  only: commander_cls_expansion
use simple_commanders_solve3D,  only: commander_solve3D_cavgs
use simple_commanders_cavgs,    only: commander_assign_cavgs_states
use simple_commanders_refine3D, only: commander_refine3D
use simple_cavg_quality_model,  only: CAVG_QUALITY_MODEL_POOL_DEFAULT
use simple_state_assignment,    only: MIN_PTCLS_PER_STATE, NSTATES_MIN, NCLS2D_FIRST, SOLVE2D_DIR, CLSEXP_DIR,&
                                     &CAVGS3D_DIR, ASSIGN_LABEL, REFINE_DIR, REFS_DIR, QUALITY_DIR, refs_check,&
                                     &assign_stats, cut_whole_stacks, usable_job_2D, count_active_classes,&
                                     &select_state_classes, map_state_poses, docked_cc_common_grid,&
                                     &set_state_cavgs3D_cline,&
                                     &expansion_budget, fit_mskdiam, select_route_classes, score_cavgs_with_model,&
                                     &model_rejected_rows, check_state_refs, read_refs_record, map_assignment,&
                                     &write_matched_cavgs, set_refine_seeds, write_refine_project, map_refine_result,&
                                     &set_solve2D_route_cline, set_clsexp_cline, set_cavgs3D_cline, set_assign_cline,&
                                     &set_refine_cont_cline
implicit none

public :: commander_solve3D_pyramid
private
#include "simple_local_flags.inc"

integer,          parameter :: NPTCLS_REFS     = 100000 ! nptcls_refs' default: the first set's selected particles at most
character(len=*), parameter :: REFS_SUBSET_DIR = 'refs_subset'
character(len=*), parameter :: ASSIGN_DIR      = 'assign'
character(len=*), parameter :: REPORT_FNAME    = 'solve3D_pyramid_report.txt'
integer,          parameter :: NSTEPS_MAX      = 32     ! the jobs timed: solve2D, two cls_expansion, the solve3D_cavgs runs, matching, refine3D
! each state's own solve3D_cavgs (provisional, until a validation run sets them)
integer,          parameter :: NCAVGS_STATE_MIN = 20    ! sub-class averages a state needs for its own run
real,             parameter :: STATE_VOL_CC_MIN = 0.80  ! its volume docked onto the state's reference, at least

type, extends(commander_base) :: commander_solve3D_pyramid
  contains
    procedure :: execute => exec_solve3D_pyramid
end type commander_solve3D_pyramid

!> what the report holds, filled as the steps finish (write_report after each)
type :: pyramid_report
    character(len=:), allocatable :: command, source, status
    integer :: nstates     = 0
    integer :: nselected   = 0      !< selected particles of the project
    integer :: nptcls_refs = 0      !< the first set's count asked for
    integer :: nfirst      = 0      !< the first set's selected particles (whole stacks)
    logical :: l_model     = .true. !< the pool model selected the first set's classes
    integer :: ncls_kept   = -1     !< the route's classes kept, of NCLS2D_FIRST (-1: not yet)
    integer :: nptcls_kept = 0      !< and their particles
    integer :: nsub_route  = 0      !< the route's sub-classes
    logical :: l_chk       = .false.
    type(refs_check) :: chk
    integer :: ncls_scored = -1     !< the project's classes the model scored (-1: not scored) and kept
    integer :: ncls_model_kept = 0
    integer :: nsub_assign = 0      !< the assignment's sub-classes
    logical :: l_assigned  = .false.
    type(assign_stats) :: stats
    integer, allocatable :: state_cavgs(:)       !< sub-class averages assigned per state
    real,    allocatable :: state_vol_cc(:)      !< each state's own volume docked onto its reference (-1: not run)
    real,    allocatable :: res0143(:), res05(:) !< the refined states' resolution (A; 0: none)
    integer :: nsteps = 0
    character(len=32) :: step_names(NSTEPS_MAX) = ''
    real              :: step_secs(NSTEPS_MAX)  = 0.
end type pyramid_report

contains

    subroutine exec_solve3D_pyramid( self, cline )
        class(commander_solve3D_pyramid), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(commander_solve2D)             :: xsolve2D
        type(commander_cls_expansion)       :: xclsexp
        type(commander_solve3D_cavgs)       :: xcavgs3D
        type(commander_assign_cavgs_states) :: xassign
        type(commander_refine3D)            :: xrefine3D
        type(parameters)          :: params
        type(sp_project)          :: spproj, sub, res
        type(cmdline)             :: cline_job
        type(pyramid_report)      :: rep
        type(string)              :: edir, source_in, stk, dir, projfile, report_fname, match_projfile, state_vol
        type(string), allocatable :: ref_vols(:)
        integer,      allocatable :: stkind(:), cls_states(:), states_matched(:)
        logical,      allocatable :: selected(:), l_first(:), model_rej(:), desel(:)
        character(len=STDLEN)     :: msg
        real    :: mskdiam, smpd, box_stk, ref_smpd, smpd_vol, cc, cc_direct, cc_mirror
        integer :: n, i, s, nstates, nparts, nsplit, nsub, nheld, nreg, ref_box, ncls, nsync, nfirst_rows
        integer :: box_vol, nsel, nposes
        ! ---- the command line, normalised before the parameters
        if( .not. cline%defined('projfile') ) THROW_HARD('solve3D_pyramid needs projfile')
        if( .not. cline%defined('nstates')  ) THROW_HARD('solve3D_pyramid needs nstates')
        if( cline%get_iarg('nstates') < NSTATES_MIN ) THROW_HARD('solve3D_pyramid needs nstates >= 2')
        if( cline%defined('mkdir') )then
            if( .not. (cline%get_carg('mkdir') == 'yes') ) THROW_HARD('solve3D_pyramid needs mkdir=yes: its steps live in its execution directory')
        endif
        call cline%set('mkdir', 'yes')
        if( .not. cline%defined('nptcls_refs') ) call cline%set('nptcls_refs', NPTCLS_REFS)
        if( cline%get_iarg('nptcls_refs') < 1 ) THROW_HARD('solve3D_pyramid needs nptcls_refs >= 1')
        if( .not. cline%defined('use_model') ) call cline%set('use_model', 'yes')
        if( .not. (cline%get_carg('use_model') == 'yes' .or. cline%get_carg('use_model') == 'no') )then
            THROW_HARD('use_model is yes or no')
        endif
        source_in = simple_abspath(cline%get_carg('projfile'))
        rep%command = command_text(cline)
        rep%source  = source_in%to_char()
        rep%status  = 'RUNNING'
        ! ---- the execution directory, the input project copied into it
        call params%new(cline)
        call simple_getcwd(edir)
        if( simple_abspath(params%projfile) == source_in )then
            THROW_HARD('solve3D_pyramid did not enter its own execution directory; its steps would overwrite the input project')
        endif
        report_fname = edir//'/'//REPORT_FNAME
        nstates = params%nstates
        nparts  = 1
        if( cline%defined('nparts') ) nparts = max(1, params%nparts)
        if( mod(params%nspace, 2) /= 0 ) THROW_HARD('solve3D_pyramid needs an even nspace (the reference spiral''s)')
        rep%nstates     = nstates
        rep%nptcls_refs = params%nptcls_refs
        rep%l_model     = params%use_model == 'yes'
        ! ---- the input: a 2D-classified project
        call spproj%read(params%projfile)
        n = spproj%os_ptcl2D%get_noris()
        if( n < 1 ) call fail('the project has no particles (ptcl2D)')
        if( spproj%os_ptcl3D%get_noris() /= n ) call fail('ptcl3D and ptcl2D of the project differ in size')
        if( .not. usable_job_2D(spproj, n) )then
            call fail('the project is not 2D-classified: class averages (one per cls2D entry), their halves or the FRCs are missing')
        endif
        ! a row is selected in both particle segments or in neither
        nsync = 0
        do i = 1,n
            if( (spproj%os_ptcl2D%get_state(i) > 0) .eqv. (spproj%os_ptcl3D%get_state(i) > 0) ) cycle
            call spproj%os_ptcl2D%set_state(i, 0)
            call spproj%os_ptcl3D%set_state(i, 0)
            nsync = nsync + 1
        enddo
        if( nsync > 0 ) write(logfhandle,'(A,I8,A)') '>>> ', nsync, ' PARTICLES SELECTED IN ONE PARTICLE SEGMENT ONLY ARE DESELECTED'
        ! the mask diameter: the command line's, else the class averages'
        mskdiam = 0.
        if( cline%defined('mskdiam') )then
            mskdiam = params%mskdiam
        else
            call spproj%get_mskdiam('cavg', mskdiam)
        endif
        if( mskdiam <= 0. ) call fail('no mask diameter: the project''s class averages record none; give mskdiam=<A>')
        write(logfhandle,'(A,F8.2,A)') '>>> SOLVE3D_PYRAMID MASK DIAMETER: ', mskdiam, ' A'
        ! ---- 1. the first set: whole stacks in arrival order
        allocate(selected(n), stkind(n))
        do i = 1,n
            selected(i) = spproj%os_ptcl2D%get_state(i) > 0
            stkind(i)   = spproj%os_ptcl2D%get_int(i, 'stkind')
        enddo
        l_first       = cut_whole_stacks(selected, stkind, params%nptcls_refs)
        rep%nselected = count(selected)
        rep%nfirst    = count(l_first)
        write(logfhandle,'(A,I8,A,I8,A,I8,A)') '>>> SOLVE3D_PYRAMID: ', rep%nselected, ' SELECTED PARTICLES; THE REFERENCES FROM THE FIRST ',&
            &rep%nfirst, ' (nptcls_refs=', params%nptcls_refs, ', WHOLE STACKS IN ARRIVAL ORDER)'
        ! a copy holding the first set's particles only (prune_particles): distributed 2D splits its
        ! parts by rows, so rows left in deselected would give parts without a particle. The first
        ! set's 2D is never mapped back onto the project's rows
        call sub%read(params%projfile)
        do i = 1,n
            if( l_first(i) ) cycle
            call sub%os_ptcl2D%set_state(i, 0)
            call sub%os_ptcl3D%set_state(i, 0)
        enddo
        call sub%prune_particles()
        nfirst_rows = sub%os_ptcl2D%get_noris()
        if( nfirst_rows /= rep%nfirst )then
            ! prune_particles also drops the particles of a stack with state 0
            write(logfhandle,'(A,I8,A)') '>>> WARNING: ', rep%nfirst - nfirst_rows,&
                &' PARTICLES OF THE FIRST SET ARE IN DESELECTED STACKS AND ARE LEFT OUT'
            rep%nfirst = nfirst_rows
        endif
        write(logfhandle,'(A,I8,A,I6,A)') '>>> THE FIRST SET''S PROJECT: ', nfirst_rows, ' PARTICLES IN ', sub%get_nstks(),&
            &' STACKS'
        call write_report(rep, report_fname)
        if( rep%nfirst < MIN_PTCLS_PER_STATE * nstates )then
            write(msg,'(A,I0,A,I0,A)') 'the first set holds ', rep%nfirst, ' selected particles; ', MIN_PTCLS_PER_STATE * nstates,&
                &' needed (5 per state)'
            call fail(trim(msg))
        endif
        call simple_mkdir(REFS_SUBSET_DIR)
        ! ---- 2. the route's solve2D
        dir      = string(REFS_SUBSET_DIR//'/'//SOLVE2D_DIR)
        projfile = dir//'/'//SOLVE2D_DIR//METADATA_EXT
        call simple_mkdir(dir)
        call sub%write(projfile)
        call set_solve2D_route_cline(cline_job, SOLVE2D_DIR//METADATA_EXT, mskdiam)
        call run_job(xsolve2D, 'solve2D (first set)', .true.)
        call take_2D(sub, nfirst_rows, 'the first set''s solve2D')
        ! ---- 3. the pool model's selection of the first set's classes (use_model)
        call select_route_classes(sub, rep%l_model, mskdiam, string(REFS_SUBSET_DIR//'/'//QUALITY_DIR))
        rep%ncls_kept   = count(sub%os_cls2D%get_all_asint('state') > 0)
        rep%nptcls_kept = sub%os_ptcl3D%count_state_gt_zero()
        call write_report(rep, report_fname)
        ! ---- 4. the route's cls_expansion
        dir      = string(REFS_SUBSET_DIR//'/'//CLSEXP_DIR)
        projfile = dir//'/'//CLSEXP_DIR//METADATA_EXT
        nsplit   = count_active_classes(sub%os_ptcl2D%get_all_asint('state'), sub%os_ptcl2D%get_all_asint('class'))
        nsub     = expansion_budget(nsplit)
        write(logfhandle,'(A,I6,A,I6,A)') '>>> CLS_EXPANSION OF THE FIRST SET: ', nsplit, ' CLASSES INTO ', nsub, ' SUBCLASSES IN ALL'
        call simple_mkdir(dir)
        call sub%write(projfile)
        call set_clsexp_cline(cline_job, CLSEXP_DIR//METADATA_EXT, nsub, fit_mskdiam(mskdiam, sub%get_box(), sub%get_smpd()))
        call run_job(xclsexp, 'cls_expansion (first set)', .true.)
        call take_2D(sub, nfirst_rows, 'the first set''s cls_expansion')
        rep%nsub_route = sub%os_cls2D%get_noris()
        ! ---- 5. solve3D_cavgs on the sub-class averages
        dir      = string(REFS_SUBSET_DIR//'/'//CAVGS3D_DIR)
        projfile = dir//'/'//CAVGS3D_DIR//METADATA_EXT
        call simple_mkdir(dir)
        call sub%write(projfile)
        call set_cavgs3D_cline(cline_job, CAVGS3D_DIR//METADATA_EXT, nstates, mskdiam)
        call run_job(xcavgs3D, 'solve3D_cavgs', .true.)
        if( .not. file_exists(projfile) ) call fail('no project from solve3D_cavgs: '//projfile%to_char())
        call res%kill
        call res%read_segment('cls3D', projfile)
        call res%read_segment('out',   projfile)
        ! ---- 6. the check; the references in state_refs/
        call check_state_refs(res, nstates, mskdiam, string(REFS_DIR), 'SOLVE3D_CAVGS', .true., rep%chk)
        call res%kill
        rep%l_chk = .true.
        call write_report(rep, report_fname)
        if( .not. rep%chk%l_vols ) call fail('solve3D_cavgs gave no volume for a state that holds class averages, or no state')
        if( .not. rep%chk%l_pass ) call fail('the reference volumes fail the check; see '//REFS_DIR//'/ and the report')
        if( .not. read_refs_record(string(REFS_DIR), nstates, ref_vols, ref_box, ref_smpd) )then
            call fail('the references recorded in '//REFS_DIR//' cannot be read')
        endif
        write(logfhandle,'(A,I3,A)') '>>> THE REFERENCES: ', size(ref_vols), ' STATES'
        call sub%kill
        ! ---- 7. the pool model on the project's classes: logged beside the fit, not applied
        if( score_cavgs_with_model(spproj, CAVG_QUALITY_MODEL_POOL_DEFAULT, mskdiam, string(ASSIGN_DIR//'/'//QUALITY_DIR),&
            &cls_states) )then
            model_rej           = model_rejected_rows(spproj, cls_states)
            rep%ncls_scored     = size(cls_states)
            rep%ncls_model_kept = count(cls_states > 0)
        endif
        ! ---- 8. cls_expansion of every selected particle's classes
        dir      = string(ASSIGN_DIR//'/'//CLSEXP_DIR)
        projfile = dir//'/'//CLSEXP_DIR//METADATA_EXT
        nsplit   = count_active_classes(spproj%os_ptcl2D%get_all_asint('state'), spproj%os_ptcl2D%get_all_asint('class'))
        nsub     = expansion_budget(nsplit)
        write(logfhandle,'(A,I6,A,I6,A)') '>>> CLS_EXPANSION OF THE PROJECT: ', nsplit, ' CLASSES INTO ', nsub, ' SUBCLASSES IN ALL'
        call simple_mkdir(dir)
        call spproj%write(projfile)
        call set_clsexp_cline(cline_job, CLSEXP_DIR//METADATA_EXT, nsub, fit_mskdiam(mskdiam, spproj%get_box(), spproj%get_smpd()))
        call run_job(xclsexp, 'cls_expansion (assignment)', .true.)
        if( .not. file_exists(projfile) ) call fail('no project from the assignment''s cls_expansion: '//projfile%to_char())
        call res%kill
        call res%read(projfile)
        if( .not. usable_job_2D(res, n) ) call fail('no usable result from the assignment''s cls_expansion: '//projfile%to_char())
        rep%nsub_assign = res%os_cls2D%get_noris()
        ! ---- 9. the matching against the references
        dir      = string(ASSIGN_DIR)
        projfile = dir//'/'//ASSIGN_LABEL//METADATA_EXT
        call res%write(projfile)
        box_stk = 0.
        call res%get_cavgs_stk(stk, ncls, smpd, fail=.false., box=box_stk)
        write(logfhandle,'(A,I6,A,I3,A)') '>>> MATCHING ', ncls, ' SUBCLASS AVERAGES AGAINST THE ', size(ref_vols), ' REFERENCES'
        call set_assign_cline(cline_job, ASSIGN_LABEL//METADATA_EXT, ref_vols, fit_mskdiam(mskdiam, nint(box_stk), smpd))
        call run_job(xassign, 'assign_cavgs_states', .false.)
        if( .not. file_exists(projfile) ) call fail('no project from the matching: '//projfile%to_char())
        call res%kill
        call res%read(projfile)
        ! ---- 10. the rejection, the flags and the mapping onto the particles
        if( res%os_cls3D%get_noris() > 0 )then
            states_matched = res%os_cls3D%get_all_asint('state')
        else
            allocate(states_matched(0))
        endif
        if( .not. map_assignment(spproj, res, size(ref_vols), 'SUBCLASSES', rep%stats, desel, model_rej) )then
            call fail('no usable matching in '//projfile%to_char()//': no candidate sub-class, or other particles')
        endif
        call res%write(projfile)
        call write_matched_cavgs(res, states_matched, size(ref_vols), string(ASSIGN_DIR), 'ASSIGNMENT')
        allocate(rep%state_cavgs(size(ref_vols)), source=0)
        do s = 1,size(ref_vols)
            rep%state_cavgs(s) = count(res%os_cls3D%get_all_asint('state') == s)
        enddo
        call res%kill
        rep%l_assigned = .true.
        call write_report(rep, report_fname)
        ! ---- 11. each state's sub-class averages refined alone: solve3D_cavgs nstates=1 with the limits
        ! it derives from their FRCs, its volume docked onto the state's reference (a run that found
        ! another structure stops here), its poses mapped onto the state's particles in place of the
        ! matching's. Each state is then in its own run's frame, which refine3D allows: it compares a
        ! particle with its own state's volume only
        match_projfile = projfile
        allocate(rep%state_vol_cc(size(ref_vols)), source=-1.)
        do s = 1,size(ref_vols)
            if( rep%state_cavgs(s) < NCAVGS_STATE_MIN )then
                write(msg,'(A,I0,A,I0,A,I0,A)') 'state ', s, ' holds ', rep%state_cavgs(s),&
                    &' sub-class averages; its own solve3D_cavgs needs ', NCAVGS_STATE_MIN, ' at least'
                call fail(trim(msg))
            endif
            dir      = string(ASSIGN_DIR//'/state_')//int2str_pad(s, 2)
            projfile = dir//'/'//CAVGS3D_DIR//METADATA_EXT
            call simple_mkdir(dir)
            call res%read(match_projfile)
            nsel = select_state_classes(res, s)
            call res%write(projfile)
            call res%kill
            write(logfhandle,'(A,I3,A,I6,A)') '>>> STATE ', s, ': SOLVE3D_CAVGS ON ITS ', nsel, ' SUBCLASS AVERAGES'
            call set_state_cavgs3D_cline(cline_job, CAVGS3D_DIR//METADATA_EXT, mskdiam)
            call run_job(xcavgs3D, 'solve3D_cavgs (state '//int2str_pad(s, 2)//')', .true.)
            if( .not. file_exists(projfile) ) call fail('no project from the state''s solve3D_cavgs: '//projfile%to_char())
            call res%read(projfile)
            if( .not. res%isthere_in_osout('vol_cavg', 1) ) call fail('the state''s solve3D_cavgs registered no volume: '//projfile%to_char())
            call res%get_vol('vol_cavg', 1, state_vol, smpd_vol, box_vol)
            if( .not. file_exists(state_vol) ) call fail('the state''s volume is missing: '//state_vol%to_char())
            cc = docked_cc_common_grid(ref_vols(s), state_vol, mskdiam, cc_direct, cc_mirror)
            rep%state_vol_cc(s) = cc
            write(logfhandle,'(A,I3,A,F7.4,A,F7.4,A)') '>>> STATE ', s, ': ITS VOLUME DOCKED ONTO ITS REFERENCE ', cc_direct,&
                &' (', cc_mirror, ' MIRRORED)'
            call write_report(rep, report_fname)
            if( cc < STATE_VOL_CC_MIN )then
                write(msg,'(A,I0,A,F6.3,A,F4.2,A)') 'state ', s, '''s own volume docks onto its reference at ', cc,&
                    &' (at least ', STATE_VOL_CC_MIN, ' needed); see '//dir%to_char()
                call fail(trim(msg))
            endif
            nposes = map_state_poses(spproj, res, s)
            call res%kill
            if( nposes < 1 ) call fail('no particle of the state takes a pose from its solve3D_cavgs: '//projfile%to_char())
            write(logfhandle,'(A,I3,A,I8,A)') '>>> STATE ', s, ': ', nposes, ' PARTICLES TAKE ITS POSES'
        enddo
        ! ---- 12. the local refinement, the states fixed
        call set_refine_seeds(spproj, params%nspace)
        dir      = string(REFINE_DIR)
        projfile = dir//'/'//REFINE_DIR//METADATA_EXT
        call simple_mkdir(dir)
        call write_refine_project(spproj, projfile, nheld)
        if( nheld < 1 ) call fail('no particle has a state after the assignment; no local refinement')
        write(logfhandle,'(A,I8,A,I3,A)') '>>> REFINING ', spproj%os_ptcl3D%count_state_gt_zero(), ' PARTICLES LOCALLY IN ', nheld,&
            &' FIXED STATES'
        call set_refine_cont_cline(cline_job, REFINE_DIR//METADATA_EXT, nheld, params%nspace, mskdiam)
        ! distributed even with one part: the distributed refine3D reconstructs its starting volumes
        call cline_job%set('nparts', nparts)
        call run_job(xrefine3D, 'refine3D refine=cont', .false.)
        if( .not. file_exists(projfile) ) call fail('no project from the local refinement: '//projfile%to_char())
        call res%read_segment('ptcl3D', projfile)
        call res%read_segment('out',    projfile)
        nreg = map_refine_result(spproj, res, nstates)
        call res%kill
        if( nreg == 0 ) call fail('the local refinement registered no volume: '//projfile%to_char())
        ! ---- 13. the result over the execution directory's project, and the report
        call state_resolutions()
        call spproj%write(params%projfile)
        rep%status = 'DONE'
        call write_report(rep, report_fname)
        write(logfhandle,'(A,A)') '>>> SOLVE3D_PYRAMID RESULT: ', params%projfile%to_char()
        write(logfhandle,'(A,A)') '>>> SOLVE3D_PYRAMID REPORT: ', report_fname%to_char()
        call spproj%kill
        call cline_job%kill
        call simple_touch(TASK_FINISHED)
        call simple_end('**** SIMPLE_SOLVE3D_PYRAMID NORMAL STOP ****')

    contains

        !> A job on cline_job (its folder dir, its project written there): the threads, with
        !! @p l_parts the parts when more than one, executed in dir and back, timed as @p label.
        subroutine run_job( x, label, l_parts )
            class(commander_base), intent(inout) :: x
            character(len=*),      intent(in)    :: label
            logical,               intent(in)    :: l_parts
            integer(8) :: t0, t1, rate
            call cline_job%set('nthr', params%nthr)
            if( l_parts .and. nparts > 1 ) call cline_job%set('nparts', nparts)
            write(logfhandle,'(A,A,A,A)') '>>> SOLVE3D_PYRAMID STEP: ', label, ' IN ', dir%to_char()
            call cline_job%printline()
            call simple_chdir(dir)
            call system_clock(t0, rate)
            call x%execute(cline_job)
            call system_clock(t1)
            call simple_chdir(edir)
            CWD_GLOB = edir%to_char()
            if( rep%nsteps < NSTEPS_MAX )then
                rep%nsteps = rep%nsteps + 1
                rep%step_names(rep%nsteps) = label
                rep%step_secs(rep%nsteps)  = real(t1 - t0) / real(rate)
            endif
            call cline_job%kill
        end subroutine run_job

        !> A 2D job's project (projfile: solve2D's, cls_expansion's) as the first set @p proj, of
        !! @p nrows particles, or exit when it is missing or unusable (usable_job_2D).
        subroutine take_2D( proj, nrows, what )
            class(sp_project), intent(inout) :: proj
            integer,           intent(in)    :: nrows
            character(len=*),  intent(in)    :: what
            if( .not. file_exists(projfile) ) call fail('no project from '//what//': '//projfile%to_char())
            call proj%kill
            call proj%read(projfile)
            if( .not. usable_job_2D(proj, nrows) ) call fail('no usable result from '//what//': '//projfile%to_char())
        end subroutine take_2D

        !> Each refined state's FSC=0.143 and FSC=0.5 resolution (A), for the report.
        subroutine state_resolutions()
            type(string)      :: vol, fsc_fname
            real, allocatable :: fsc(:), resarr(:)
            real    :: smpd_vol
            integer :: box_vol, fsc_box, m
            allocate(rep%res0143(size(ref_vols)), rep%res05(size(ref_vols)), source=0.)
            do s = 1,size(ref_vols)
                if( .not. spproj%isthere_in_osout('vol', s) .or. .not. spproj%isthere_in_osout('fsc', s) ) cycle
                call spproj%get_vol('vol', s, vol, smpd_vol, box_vol)
                call spproj%get_fsc(s, fsc_fname, fsc_box)
                if( .not. file_exists(fsc_fname) ) cycle
                fsc    = file2rarr(fsc_fname)
                resarr = get_resarr(fsc_box, smpd_vol)
                m      = min(size(fsc), size(resarr))
                if( m == 0 ) cycle
                call get_resolution(fsc(:m), resarr(:m), rep%res05(s), rep%res0143(s))
                write(logfhandle,'(A,I3,A,F7.2,A,F7.2,A)') '>>>   STATE ', s, ': FSC=0.143 ', rep%res0143(s), ' A, FSC=0.5 ',&
                    &rep%res05(s), ' A'
            enddo
        end subroutine state_resolutions

        !> The failed step @p why recorded in the report, then the hard exit.
        subroutine fail( why )
            character(len=*), intent(in) :: why
            rep%status = 'FAILED: '//why
            if( report_fname%strlen() > 0 ) call write_report(rep, report_fname)
            THROW_HARD(why)
        end subroutine fail

    end subroutine exec_solve3D_pyramid

    !> The command line as normalised, for the report.
    function command_text( cline ) result( text )
        class(cmdline), intent(inout) :: cline
        character(len=:), allocatable :: text
        type(chash)  :: job_descr
        type(string) :: str
        call cline%gen_job_descr(job_descr)
        str  = job_descr%chash2str()
        text = 'simple_exec '//str%to_char()
        call job_descr%kill
    end function command_text

    !> The report so far, a text table (REPORT_FNAME): the first set, the route, the check, the
    !! assignment, per state its sub-class averages, particles and resolution, and the steps' times.
    subroutine write_report( rep, fname )
        class(pyramid_report), intent(in) :: rep
        class(string),         intent(in) :: fname
        real    :: vol_cc
        integer :: funit, ios, s, s2, nptcls
        call fopen(funit, file=fname, status='replace', action='write', iostat=ios)
        if( ios /= 0 )then
            THROW_WARN('could not write the solve3D_pyramid report: '//fname%to_char())
            return
        endif
        write(funit,'(A)') 'SOLVE3D_PYRAMID REPORT'
        write(funit,'(A)') 'status: '//rep%status
        write(funit,'(A)') 'source project: '//rep%source
        write(funit,'(A)') 'command: '//rep%command
        write(funit,'(A,I0,A,I0)') 'states: ', rep%nstates, '  selected particles: ', rep%nselected
        write(funit,'(A,I0,A,I0,A)') 'first set: ', rep%nfirst, ' selected particles (nptcls_refs=', rep%nptcls_refs,&
            &'; whole stacks in arrival order)'
        if( rep%ncls_kept >= 0 )then
            write(funit,'(A,I0,A,I0,A,I0,A,A,A)') 'route: solve2D keeps ', rep%ncls_kept, ' of ', NCLS2D_FIRST, ' classes and ',&
                &rep%nptcls_kept, ' particles (use_model=', trim(merge('yes', 'no ', rep%l_model)), ')'
        endif
        if( rep%nsub_route > 0 ) write(funit,'(A,I0,A)') 'route: cls_expansion into ', rep%nsub_route, ' sub-classes'
        if( rep%l_chk )then
            write(funit,'(A)') ''
            write(funit,'(A)') 'check of the reference volumes'
            write(funit,'(A,F8.2,A,F8.2,A)') '  even/odd angular distance (median, halves in one state): ', rep%chk%angdist,&
                &' deg (mean over all ', rep%chk%angdist_all, ')'
            write(funit,'(A,F8.2,A)') '  even/odd state agreement: ', rep%chk%overlap, ' %'
            write(funit,'(A,F8.4)') '  largest docked correlation of two states: ', rep%chk%cc_max
            if( allocated(rep%chk%cc_pairs) )then
                do s = 1,size(rep%chk%cc_pairs, 1) - 1
                    do s2 = s + 1,size(rep%chk%cc_pairs, 2)
                        write(funit,'(A,I0,A,I0,A,F8.4)') '  states ', s, ' and ', s2, ': docked correlation ', rep%chk%cc_pairs(s,s2)
                    enddo
                enddo
            endif
            if( allocated(rep%chk%counts) )then
                do s = 1,size(rep%chk%counts)
                    write(funit,'(A,I0,A,I0,A)') '  state ', s, ': ', rep%chk%counts(s), ' class averages'
                enddo
            endif
            if( rep%chk%l_pass )then
                write(funit,'(A)') '  verdict: PASS'
            else
                write(funit,'(A)') '  verdict: FAIL'
            endif
        endif
        if( rep%l_assigned )then
            write(funit,'(A)') ''
            write(funit,'(A,I0,A)') 'assignment of ', rep%stats%ncls, ' sub-classes'
            if( rep%ncls_scored >= 0 ) write(funit,'(A,I0,A,I0,A)') '  pool model (logged, not applied): keeps ', &
                &rep%ncls_model_kept, ' of ', rep%ncls_scored, ' classes'
            write(funit,'(A,I0,A,I0,A,I0,A,I0)') '  candidates ', rep%stats%ncand, '  rejected by fit ', rep%stats%nrej,&
                &' (at most ', rep%stats%nmax, ')  ambiguous ', rep%stats%namb
            write(funit,'(A,F8.4,A,F8.4,A,F8.4)') '  fit: median best correlation ', rep%stats%med, '  MAD ', rep%stats%mad,&
                &'  threshold ', rep%stats%thres
            write(funit,'(A,I0,A,I0)') '  particles assigned ', rep%stats%nassigned, '  deselected ', rep%stats%ndesel
            if( rep%stats%l_model ) write(funit,'(A,I0,A,I0,A,I0)') '  particles rejected by the fit and the model ',&
                &rep%stats%nboth, '  the model alone ', rep%stats%nmodel, '  the fit alone ', rep%stats%nfit
        endif
        if( allocated(rep%state_cavgs) )then
            write(funit,'(A)') ''
            write(funit,'(A)') 'state  sub-class averages  particles  own volume cc  FSC0.143(A)  FSC0.5(A)'
            do s = 1,size(rep%state_cavgs)
                nptcls = 0
                if( allocated(rep%stats%counts) )then
                    if( s <= size(rep%stats%counts) ) nptcls = rep%stats%counts(s)
                endif
                vol_cc = -1.
                if( allocated(rep%state_vol_cc) )then
                    if( s <= size(rep%state_vol_cc) ) vol_cc = rep%state_vol_cc(s)
                endif
                if( allocated(rep%res0143) )then
                    write(funit,'(I5,2X,I18,2X,I9,2X,F13.4,2X,F11.2,2X,F9.2)') s, rep%state_cavgs(s), nptcls, vol_cc,&
                        &rep%res0143(s), rep%res05(s)
                else
                    write(funit,'(I5,2X,I18,2X,I9,2X,F13.4,2X,A11,2X,A9)') s, rep%state_cavgs(s), nptcls, vol_cc, '-', '-'
                endif
            enddo
            write(funit,'(A)') '(own volume cc: the state''s own solve3D_cavgs volume docked onto its reference; -1 not run)'
        endif
        if( rep%nsteps > 0 )then
            write(funit,'(A)') ''
            write(funit,'(A)') 'step                              wall(s)'
            do s = 1,rep%nsteps
                write(funit,'(A32,2X,F9.1)') rep%step_names(s), rep%step_secs(s)
            enddo
        endif
        call fclose(funit)
    end subroutine write_report

end module simple_commanders_solve3D_pyramid
