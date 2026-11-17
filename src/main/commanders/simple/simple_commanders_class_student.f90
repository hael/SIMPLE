!@descr: the class student's commanders: learn_class_student (the teacher's rows against frozen class averages and the training), predict_classes (the class, the angle, the candidates and the confidence of unseen particles) and the validate_class_student controller
module simple_commanders_class_student
use simple_commander_base,         only: commander_base
use simple_cmdline,                only: cmdline
use simple_strategy_class_student, only: run_learn_class_student, run_class_posterior_mode, run_predict_classes
use simple_defs_hhcnn,             only: KMAX_CLASS, STUDENT_SHIFT_SIGMA, STUDENT_TRS, STUDENT_LP_VALIDATE, CLASS_DENSE_SHIFT_SIGMA, &
    &CLASS_MINPOP_FRC, CLASS_STUDENT_PROB_BETA, CLASS_REFERENCE_FNAME, CLASS_REFERENCE_CAVGS, STUDENT_TRAIN_PROJ, STUDENT_UNSEEN_PROJ, &
    &CLASS_STUDENT_MODEL_FNAME, CLASS_STUDENT_PRED_TABLE, CLASS_STUDENT_PRED_TIMES, CLASS_STUDENT_CAVGS, &
    &VALIDATE_CLASS_STUDENT_REPORT, VALIDATE_CLASS_STUDENT_DETAILS
implicit none
private
public :: commander_learn_class_student, commander_predict_classes, commander_validate_class_student
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_learn_class_student
  contains
    procedure :: execute => exec_learn_class_student
end type commander_learn_class_student

type, extends(commander_base) :: commander_predict_classes
  contains
    procedure :: execute => exec_predict_classes
end type commander_predict_classes

type, extends(commander_base) :: commander_validate_class_student
  contains
    procedure :: execute => exec_validate_class_student
end type commander_validate_class_student


contains

    !> The teacher's rows of the project's selected particles against the class averages refs at the band lp, with
    !! their features, and the student trained on them, written as the model fname with its weights beside it
    !! (simple_private_exec). The dense table's mode and objective are set here: the student has no key of its own
    subroutine exec_learn_class_student( self, cline )
        class(commander_learn_class_student), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'no')
        if( .not. cline%defined('trs')         ) call cline%set('trs',         STUDENT_TRS)
        if( .not. cline%defined('shift_sigma') ) call cline%set('shift_sigma', STUDENT_SHIFT_SIGMA)
        if( .not. cline%defined('fname')       ) call cline%set('fname',       CLASS_STUDENT_MODEL_FNAME)
        call set_ncls_from_refs(cline)
        call cline%set('oritype', 'ptcl2D')
        call cline%set('objfun',  'euclid')
        call cline%set('refine',  'prob')
        call cline%set('center',  'no')
        ! the physical teacher: the posterior a solve2D run ends at, no tempering
        call cline%set('prob_beta', CLASS_STUDENT_PROB_BETA)
        call run_learn_class_student(cline)
    end subroutine exec_learn_class_student

    !> The student's prediction for the project's selected particles under the model infile: per particle the class,
    !! its in-plane angle, the candidates holding 99% of the posterior and the confidence, and, given the class averages
    !! refs, the alignment over the reference pass's disk (simple_private_exec). With nparts it is the master of that
    !! many parts; a part is given its range as the master hands it (part, nparts, fromp, top, outfile)
    subroutine exec_predict_classes( self, cline )
        use simple_error, only: simple_exception
        class(commander_predict_classes), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        if( cline%defined('part') )then
            if( .not. (cline%defined('nparts') .and. cline%defined('fromp') .and. cline%defined('top')) ) &
                &THROW_HARD('a part of predict_classes needs its range: nparts, fromp and top')
            if( .not. cline%defined('outfile') ) THROW_HARD('a part of predict_classes needs outfile')
        endif
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'no')
        if( .not. cline%defined('fname') ) call cline%set('fname', CLASS_STUDENT_MODEL_FNAME)
        call cline%set('oritype', 'ptcl2D')
        call cline%set('objfun',  'euclid')
        if( cline%defined('refs') )then
            call cline%set('refine', 'prob')
            call cline%set('center', 'no')
            ! the alignment's disk and prior are the reference pass's, at the physical temperature
            if( .not. cline%defined('trs')         ) call cline%set('trs',         STUDENT_TRS)
            if( .not. cline%defined('shift_sigma') ) call cline%set('shift_sigma', CLASS_DENSE_SHIFT_SIGMA)
            call cline%set('prob_beta', CLASS_STUDENT_PROB_BETA)
        endif
        call run_predict_classes(cline)
    end subroutine exec_predict_classes

    !> validate_class_student: solve2D on a random nsample (the teacher), the student trained, the unseen predicted and
    !! aligned at their picked positions and searched densely, the class averages, the reads and the timings
    !! (doc/algorithms/class_student.md). Every step but the training runs over nparts when given
    subroutine exec_validate_class_student( self, cline )
        use simple_core_module_api
        use simple_parameters,              only: parameters
        use simple_sp_project,              only: sp_project
        use simple_commanders_solve2D,      only: commander_solve2D
        use simple_commanders_mkcavgs,      only: commander_make_cavgs
        use simple_eul_prob_tab_utils,      only: disk_radius_px
        class(commander_validate_class_student), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        integer,          parameter :: NPHASES = 9    !< the validation's timed phases
        integer,          parameter :: NARM    = 6    !< arm_stats' reads
        integer,          parameter :: NALN    = 6    !< aln_stats' reads
        character(len=*), parameter :: PHASE_NAMES(NPHASES) = [character(len=48) :: &
            &'Drawing the training set', 'Teacher: solve2D on the training set', &
            &'Training the student (learn_class_student)', 'Prediction and alignment (predict_classes)', &
            &'Reference: the physical posterior''s best class', 'Class averages of the student''s alignment', &
            &'Statistics', 'Project copies and bookkeeping', 'Prediction at the reference''s shifts (details)']
        integer, parameter :: PHASE_ORDER(NPHASES) = [1, 2, 3, 4, 5, 9, 6, 7, 8] !< the run's order
        type(commander_solve2D)             :: xsolve2D
        type(commander_make_cavgs)          :: xmkcavgs
        type(commander_learn_class_student) :: xlearn
        type(commander_predict_classes)     :: xpredict
        type(parameters)  :: params
        type(sp_project)  :: spproj, res
        type(cmdline)     :: cline_entry, cl
        type(ran_tabu)    :: rt
        type(string)      :: edir, result_proj, cavgs, sigma2_abs, model_abs, fname, dense_cavgs, student_cavgs
        integer, allocatable :: pinds(:), cls_s(:), cls_d(:), ncand(:), route(:), cands(:,:)
        integer, allocatable :: cls_i(:), ncand_i(:), route_i(:), cands_i(:,:)
        logical, allocatable :: active(:), train(:), unseen(:)
        real,    allocatable :: sh_s(:,:), sh_d(:,:), e3_s(:), e3_d(:)
        real(timer_int_kind) :: rt_ph(NPHASES), rt_all
        integer(timer_int_kind) :: t0, t_all
        real    :: smpd_cavgs, t_pred, t_aln
        integer :: nptcls, nactive, ntrain, nunseen, ncls_t, i, funit, dunit, io_stat
        logical :: l_sigma2, l_times
        t_all = tic()
        rt_ph = 0.
        ! ---- the command line: an execution directory, the training sample's size, the classes
        if( .not. cline%defined('nsample') ) THROW_HARD('validate_class_student needs nsample, the training sample size')
        if( .not. cline%defined('ncls')    ) THROW_HARD('validate_class_student needs ncls, the teacher''s classes')
        if( cline%defined('mkdir') )then
            if( .not. (cline%get_carg('mkdir') == 'yes') ) THROW_HARD('validate_class_student needs mkdir=yes: its steps live in its execution directory')
        endif
        call cline%set('mkdir', 'yes')
        call cline%set('oritype', 'ptcl2D')
        cline_entry = cline
        call params%new(cline)
        call simple_getcwd(edir)
        ! ---- the training sample: a random draw of nsample active particles; the others are unseen
        t0 = tic()
        call spproj%read(params%projfile)
        nptcls = spproj%os_ptcl2D%get_noris()
        if( nptcls < 1 ) THROW_HARD('the project has no particles (ptcl2D)')
        allocate(active(nptcls), train(nptcls), unseen(nptcls), source=.false.)
        do i = 1,nptcls
            active(i) = spproj%os_ptcl2D%get_state(i) > 0
        enddo
        nactive = count(active)
        call spproj%os_ptcl2D%get_pinds(1, 'state', pinds)
        if( .not. allocated(pinds) ) THROW_HARD('no particle in state 1')
        if( size(pinds) /= nactive ) THROW_HARD('validate_class_student takes particles of state 0 and 1 only')
        ntrain = min(params%nsample, nactive)
        if( ntrain >= nactive ) THROW_HARD('nsample leaves no unseen particle')
        rt = ran_tabu(nactive)
        call rt%shuffle(pinds)
        call rt%kill
        train(pinds(1:ntrain)) = .true.
        unseen  = active .and. .not. train
        nunseen = count(unseen)
        write(logfhandle,'(A,I0,A,I0,A)') '>>> VALIDATE_CLASS_STUDENT: ', ntrain, ' training particles, ', nunseen, ' unseen'
        rt_ph(1) = toc(t0)
        t0 = tic()
        call spproj%os_ptcl2D%delete_2Dclustering(keepshifts=.true.)
        call select_rows(spproj, train)
        fname = edir%to_char()//'/'//STUDENT_TRAIN_PROJ
        call spproj%update_projinfo(fname)
        call spproj%write(fname)
        rt_ph(8) = rt_ph(8) + toc(t0)
        ! ---- the teacher: solve2D on every training particle (nsample at the sample's size: full sampling)
        write(logfhandle,'(A)') '>>> VALIDATE_CLASS_STUDENT: SOLVE2D ON THE TRAINING SET'
        t0 = tic()
        cl = cline_entry
        call cl%set('prg',      'solve2D')
        call cl%set('projfile', fname%to_char())
        call cl%set('nsample',  ntrain)
        call cl%set('mkdir',    'yes')
        call xsolve2D%execute(cl)
        call simple_chdir(edir)
        result_proj = simple_abspath(cl%get_carg('projfile'))
        call cl%kill
        rt_ph(2) = toc(t0)
        ! ---- the student: trained on the teacher's rows against its class averages, at the teacher's shifts
        t0 = tic()
        call res%read(result_proj)
        call res%get_cavgs_stk(cavgs, ncls_t, smpd_cavgs)
        cavgs = simple_abspath(cavgs)
        call res%get_sigma2_state_path(sigma2_abs, l_sigma2)
        if( l_sigma2 ) l_sigma2 = file_exists(sigma2_abs)
        if( l_sigma2 ) sigma2_abs = simple_abspath(sigma2_abs)
        allocate(cls_s(nptcls), cls_d(nptcls), ncand(nptcls), route(nptcls), source=0)
        allocate(cls_i(nptcls), ncand_i(nptcls), route_i(nptcls), source=0)
        allocate(cands(KMAX_CLASS,nptcls), cands_i(KMAX_CLASS,nptcls), source=0)
        allocate(sh_s(2,nptcls), sh_d(2,nptcls), source=0.)
        allocate(e3_s(nptcls), e3_d(nptcls), source=0.)
        call enter_subdir('learn')
        call write_copy(res, STUDENT_TRAIN_PROJ)
        rt_ph(8) = rt_ph(8) + toc(t0)
        write(logfhandle,'(A)') '>>> VALIDATE_CLASS_STUDENT: STUDENT TRAINING'
        t0 = tic()
        call cl%set('prg',      'learn_class_student')
        call cl%set('projfile', STUDENT_TRAIN_PROJ)
        call cl%set('refs',     cavgs%to_char())
        call cl%set('lp',       STUDENT_LP_VALIDATE)
        call cl%set('mskdiam',  params%mskdiam)
        call cl%set('nthr',     params%nthr)
        call cl%set('fname',    CLASS_STUDENT_MODEL_FNAME)
        call xlearn%execute(cl)
        call cl%kill
        call simple_chdir(edir)
        model_abs = simple_abspath(string('learn/'//CLASS_STUDENT_MODEL_FNAME))
        rt_ph(3) = toc(t0)
        ! ---- the unseen particles at their picked positions (no centring): shift zero, not yet centred, the same
        !      project for the prediction and the reference
        t0 = tic()
        call res%os_ptcl2D%delete_2Dclustering(keepshifts=.true.)
        call select_rows(res, unseen)
        do i = 1,nptcls
            if( .not. unseen(i) ) cycle
            call res%os_ptcl2D%set_shift(i, [0.,0.])
            call res%os_ptcl2D%set(i, 'shsig',  -1.)
            call res%os_ptcl2D%set(i, 'shincr',  0.)
        enddo
        call enter_subdir('dense')
        call write_copy(res, STUDENT_UNSEEN_PROJ)
        call enter_subdir('predict')
        call write_copy(res, STUDENT_UNSEEN_PROJ)
        rt_ph(8) = rt_ph(8) + toc(t0)
        write(logfhandle,'(A)') '>>> VALIDATE_CLASS_STUDENT: STUDENT PREDICTION AND ALIGNMENT FOR THE UNSEEN PARTICLES'
        t0 = tic()
        call cl%set('prg',      'predict_classes')
        call cl%set('projfile', STUDENT_UNSEEN_PROJ)
        call cl%set('infile',   model_abs%to_char())
        call cl%set('refs',     cavgs%to_char())
        call cl%set('mskdiam',  params%mskdiam)
        call cl%set('nthr',     params%nthr)
        call cl%set('fname',    CLASS_STUDENT_MODEL_FNAME)
        call set_parts(cl)
        call xpredict%execute(cl)
        call cl%kill
        rt_ph(4) = toc(t0)
        t0 = tic()
        call read_pred_table(string(CLASS_STUDENT_PRED_TABLE), ncand, route, cands)
        call read_pred_times(string(CLASS_STUDENT_PRED_TIMES), t_pred, t_aln, l_times)
        call res%kill
        call res%read(string(STUDENT_UNSEEN_PROJ))
        do i = 1,nptcls
            if( .not. unseen(i) ) cycle
            cls_s(i)  = res%os_ptcl2D%get_class(i)
            e3_s(i)   = res%os_ptcl2D%e3get(i)
            sh_s(:,i) = res%os_ptcl2D%get_2Dshift(i)
        enddo
        call simple_chdir(edir)
        rt_ph(8) = rt_ph(8) + toc(t0)
        ! ---- the reference the student is judged against: the physical posterior's best class of the same particles
        !      against the teacher's averages (the class of largest in-plane marginal mass at the physical temperature,
        !      as the hard assignment at the end of a solve2D run takes it), every particle searching the
        !      disk from its picked position, at the student's box and band, with the class averages of that assignment
        call enter_subdir('dense')
        write(logfhandle,'(A)') '>>> VALIDATE_CLASS_STUDENT: THE PHYSICAL POSTERIOR''S BEST CLASS OF THE UNSEEN PARTICLES'
        t0 = tic()
        call cl%set('prg',         'learn_class_student')
        call cl%set('projfile',    STUDENT_UNSEEN_PROJ)
        call cl%set('oritype',     'ptcl2D')
        call cl%set('ncls',        ncls_t)
        call cl%set('refs',        cavgs%to_char())
        call cl%set('refine',      'prob')
        call cl%set('objfun',      'euclid')
        call cl%set('lp',          STUDENT_LP_VALIDATE)
        call cl%set('trs',         STUDENT_TRS)
        call cl%set('shift_sigma', CLASS_DENSE_SHIFT_SIGMA)
        call cl%set('prob_beta',   CLASS_STUDENT_PROB_BETA)
        call cl%set('center',      'no')
        call cl%set('mskdiam',     params%mskdiam)
        call cl%set('nthr',        params%nthr)
        call cl%set('fname',       CLASS_REFERENCE_FNAME)
        call cl%set('mkdir',       'no')
        call run_class_posterior_mode(cl)
        call cl%kill
        call cl%set('prg',      'make_cavgs')
        call cl%set('projfile', STUDENT_UNSEEN_PROJ)
        call cl%set('ncls',     ncls_t)
        call cl%set('refs',     CLASS_REFERENCE_CAVGS)
        call cl%set('mskdiam',  params%mskdiam)
        call cl%set('nthr',     params%nthr)
        call cl%set('mkdir',    'no')
        call set_parts(cl)
        call xmkcavgs%execute(cl)
        call cl%kill
        call res%kill
        call res%read(string(STUDENT_UNSEEN_PROJ))
        do i = 1,nptcls
            if( .not. unseen(i) ) cycle
            cls_d(i)  = res%os_ptcl2D%get_class(i)
            e3_d(i)   = res%os_ptcl2D%e3get(i)
            sh_d(:,i) = res%os_ptcl2D%get_2Dshift(i)
        enddo
        dense_cavgs = simple_abspath(string(CLASS_REFERENCE_CAVGS))
        call simple_chdir(edir)
        rt_ph(5) = toc(t0)
        ! ---- the network at the reference's shifts: the student's prediction of the same particles, no alignment
        t0 = tic()
        call enter_subdir('predict_dense')
        call write_copy(res, STUDENT_UNSEEN_PROJ)
        call cl%set('prg',      'predict_classes')
        call cl%set('projfile', STUDENT_UNSEEN_PROJ)
        call cl%set('infile',   model_abs%to_char())
        call cl%set('mskdiam',  params%mskdiam)
        call cl%set('nthr',     params%nthr)
        call cl%set('fname',    CLASS_STUDENT_MODEL_FNAME)
        call set_parts(cl)
        call xpredict%execute(cl)
        call cl%kill
        call read_pred_table(string(CLASS_STUDENT_PRED_TABLE), ncand_i, route_i, cands_i)
        call simple_chdir(edir)
        rt_ph(9) = toc(t0)
        ! ---- the class averages of the student's assignments (the prediction's aligned records)
        t0 = tic()
        call res%kill
        call res%read(string(edir%to_char()//'/predict/'//STUDENT_UNSEEN_PROJ))
        call enter_subdir('cavgs')
        call write_copy(res, STUDENT_UNSEEN_PROJ)
        call make_cavgs_here(CLASS_STUDENT_CAVGS)
        student_cavgs = simple_abspath(string(CLASS_STUDENT_CAVGS))
        call simple_chdir(edir)
        rt_ph(6) = toc(t0)
        ! ---- the statistics and the report
        t0 = tic()
        call fopen(funit, string(VALIDATE_CLASS_STUDENT_REPORT), status='replace', action='write', iostat=io_stat)
        call fileiochk('validate_class_student; open '//VALIDATE_CLASS_STUDENT_REPORT, io_stat)
        call fopen(dunit, string(VALIDATE_CLASS_STUDENT_DETAILS), status='replace', action='write', iostat=io_stat)
        call fileiochk('validate_class_student; open '//VALIDATE_CLASS_STUDENT_DETAILS, io_stat)
        rt_ph(7) = toc(t0)
        rt_all = toc(t_all)
        call report_statistics(funit, dunit)
        call fclose(funit)
        call fclose(dunit)
        call cline_entry%kill
        call spproj%kill
        call res%kill
        call simple_end('**** VALIDATE_CLASS_STUDENT NORMAL STOP ****')

    contains

        !> a step's command line given the entry's parts (the teacher's solve2D inherits them with the entry's line)
        subroutine set_parts( cl_here )
            type(cmdline), intent(inout) :: cl_here
            if( cline_entry%defined('nparts') ) call cl_here%set('nparts', params%nparts)
        end subroutine set_parts

        subroutine enter_subdir( dir )
            character(len=*), intent(in) :: dir
            call simple_chdir(edir)
            call simple_mkdir(dir)
            call simple_chdir(string(dir))
        end subroutine enter_subdir

        !> the rows of mask in state 1, the others in state 0, in both particle fields
        subroutine select_rows( proj, mask )
            type(sp_project), intent(inout) :: proj
            logical,          intent(in)    :: mask(nptcls)
            integer :: j
            do j = 1,nptcls
                call proj%os_ptcl2D%set_state(j, merge(1, 0, mask(j)))
                if( proj%os_ptcl3D%get_noris() == nptcls ) call proj%os_ptcl3D%set_state(j, merge(1, 0, mask(j)))
            enddo
        end subroutine select_rows

        !> the project written as fname in the working directory with its own copy of the teacher's sigma2 state (a
        !! child may commit to it; the teacher's stays frozen)
        subroutine write_copy( proj, fname_here )
            type(sp_project), intent(inout) :: proj
            character(len=*), intent(in)    :: fname_here
            type(string) :: base
            if( l_sigma2 )then
                base = basename(sigma2_abs)
                call simple_copy_file(sigma2_abs, base)
                call proj%set_sigma2_state_path(base)
            endif
            call proj%update_projinfo(string(fname_here))
            call proj%write(string(fname_here))
        end subroutine write_copy

        !> the class averages of the working directory's unseen project as refs_here
        subroutine make_cavgs_here( refs_here )
            character(len=*), intent(in) :: refs_here
            call cl%set('prg',      'make_cavgs')
            call cl%set('projfile', STUDENT_UNSEEN_PROJ)
            call cl%set('ncls',     ncls_t)
            call cl%set('refs',     refs_here)
            call cl%set('mskdiam',  params%mskdiam)
            call cl%set('nthr',     params%nthr)
            call cl%set('mkdir',    'no')
            call set_parts(cl)
            call xmkcavgs%execute(cl)
            call cl%kill
        end subroutine make_cavgs_here

        !> the candidates, their count and the route of every predicted row
        subroutine read_pred_table( table, ncand_o, route_o, cands_o )
            class(string), intent(in)    :: table
            integer,       intent(inout) :: ncand_o(nptcls), route_o(nptcls), cands_o(KMAX_CLASS,nptcls)
            character(len=4096) :: line
            integer :: u, ios, iptcl, eo, cls, kk, nc, rte, cstored, j
            real    :: pmax, cf, e3, x, y
            call fopen(u, table, status='old', action='read', iostat=ios)
            call fileiochk('validate_class_student; open '//table%to_char(), ios)
            do
                read(u,'(A)', iostat=ios) line
                if( ios /= 0 ) exit
                if( line(1:1) == '#' ) cycle
                read(line,*) iptcl, eo, cls, pmax, kk, nc, cf, rte, e3, x, y, cstored, (cands_o(j,iptcl), j = 1,min(nc,KMAX_CLASS))
                ncand_o(iptcl) = nc
                route_o(iptcl) = rte
            enddo
            call fclose(u)
        end subroutine read_pred_table

        !> the prediction's wall times (a shared-memory run's file; parts write none)
        subroutine read_pred_times( fname_t, tp, ta, l_found )
            class(string), intent(in)  :: fname_t
            real,          intent(out) :: tp, ta
            logical,       intent(out) :: l_found
            character(len=256) :: line
            integer :: u, ios, n
            tp = -1.
            ta = -1.
            l_found = file_exists(fname_t)
            if( .not. l_found ) return
            call fopen(u, fname_t, status='old', action='read', iostat=ios)
            call fileiochk('validate_class_student; open '//fname_t%to_char(), ios)
            do
                read(u,'(A)', iostat=ios) line
                if( ios /= 0 ) exit
                if( line(1:1) == '#' ) cycle
                read(line,*) n, tp, ta
            enddo
            call fclose(u)
        end subroutine read_pred_times

        !> the verdict (can the student replace the full search?) with how to read it and the run time, to unit u and the
        !! log; the diagnostics to unit ud (doc/algorithms/class_student.md, section Validation)
        subroutine report_statistics( u, ud )
            use simple_class_student, only: class_student
            integer, intent(in) :: u, ud
            character(len=96), parameter :: RULE = repeat('-', 96)
            type(class_student) :: stu
            real    :: thr, target, hs, hr, cost_p, cost_a, cost_r, cost_d, rdisk, va(NARM), vi(NARM)
            real    :: av(NALN), frc(3,2), pick(2)
            integer :: iu, units(2), ncmp
            units = [u, logfhandle]
            va = arm_stats(ncand, route, cands)
            vi = arm_stats(ncand_i, route_i, cands_i)
            av = aln_stats(cls_s, e3_s, sh_s)
            pick = picking_error()
            call cavgs_frc(frc, ncmp)
            call stu%read(model_abs)
            call stu%get_threshold(thr, target, hs, hr)
            call stu%kill
            target = 100. * target
            rdisk  = disk_radius_px(CLASS_DENSE_SHIFT_SIGMA, params%smpd, params%box) * params%smpd
            ! the costs per unseen particle: the prediction's wall time with reading and preparation, the alignment's
            ! wall time, and the route (both); the full search's wall time
            cost_d = 1000. * real(rt_ph(5)) / real(max(1,nunseen))
            if( l_times )then
                cost_a = 1000. * t_aln / real(max(1,nunseen))
                cost_p = 1000. * (real(rt_ph(4)) - t_aln) / real(max(1,nunseen))
            else
                cost_a = -1.
                cost_p = -1.
            endif
            cost_r = 1000. * real(rt_ph(4)) / real(max(1,nunseen))
            do iu = 1,2
                associate( uu => units(iu) )
                write(uu,'(a)') 'CLASS STUDENT VALIDATION'
                write(uu,'(a)') repeat('=', 96)
                write(uu,'(a,i0,a,i0,a)') 'The student learns the classes of a solve2D run on ', ntrain, ' particles (', ncls_t, &
                    &' classes). For each'
                write(uu,'(a)') 'new particle, at its picked position (no centring), it ranks the classes by how likely they are.'
                write(uu,'(a)') 'Its candidate list holds the classes making up 99% of that ranking (at most 16). Where it is'
                write(uu,'(a)') 'confident (its calibrated threshold), the list is the particle''s shortlist: the only classes a'
                write(uu,'(a)') 'later search tests. Every particle is then aligned to the teacher''s class averages of its'
                write(uu,'(a,a,a)') 'candidate list: every in-plane rotation and shifts within ', r2s(rdisk, 1), &
                    &' A of its picked position;'
                write(uu,'(a)') 'the best alignment is kept.'
                write(uu,'(a,i0,a)') 'Checked on ', nunseen, ' particles it has not seen, against a full search (the reference).'
                write(uu,'(a)') ''
                write(uu,'(a,a10,3x,a)') lj('', 54), 'result', 'goal'
                write(uu,'(a)') RULE
                write(uu,'(2x,a,2x,f10.1,3x,a)') lj('Particles given a shortlist (%)', 50), va(3), 'as many as possible'
                write(uu,'(2x,a,2x,f10.1,3x,a,a)') lj('Shortlists containing the reference''s class (%)', 50), va(4), &
                    &lj(r2s(target, 1)//' or more', 26), met(va(4) >= target)
                write(uu,'(2x,a,2x,f10.1,3x,a,a)') lj('  same, predicted at the reference''s shifts (%)', 50), vi(4), &
                    &lj(r2s(target, 1)//' or more', 26), met(vi(4) >= target)
                if( va(6) > 0. )then
                    write(uu,'(2x,a,2x,i10,3x,a,i0,a)') lj('Shortlist length (median)', 50), nint(va(6)), 'of ', ncls_t, ' classes'
                else
                    write(uu,'(2x,a,2x,a10,3x,a,i0,a)') lj('Shortlist length (median)', 50), '-', 'of ', ncls_t, ' classes'
                endif
                write(uu,'(2x,a,2x,f10.1,3x,a)') lj('Aligned class = the reference''s class (%)', 50), av(1), &
                    &'higher is better'
                if( ncmp > 0 )then
                    write(uu,'(2x,a,2x,f10.2,3x,a)') lj('Class averages vs the teacher''s, FRC 0.5 (A)', 50), frc(2,1), &
                        &'reference '//r2s(frc(2,2), 2)//'; lower is better'
                else
                    write(uu,'(2x,a,2x,a10,3x,a)') lj('Class averages vs the teacher''s, FRC 0.5 (A)', 50), '-', 'no class compared'
                endif
                write(uu,'(2x,a,2x,f10.3,3x,a,a)') lj('Cost per particle (ms)', 50), cost_r, &
                    &lj('below full search '//r2s(cost_d, 3), 26), met(cost_r < cost_d)
                write(uu,'(a)') ''
                if( va(4) >= target )then
                    write(uu,'(a)') 'The shortlists are safe: they hold the reference''s class as often as the goal asks.'
                else if( vi(4) >= target )then
                    write(uu,'(a)') 'The shortlists miss the reference''s class because of the picking error: predicted at the'
                    write(uu,'(a)') 'reference''s shifts they meet the goal.'
                else
                    write(uu,'(a)') 'The shortlists miss the reference''s class even at the reference''s shifts: the student'
                    write(uu,'(a)') 'itself falls short.'
                endif
                if( cost_r >= cost_d )then
                    write(uu,'(a,a,a,a,a)') 'The student''s route costs more than the full search: prediction ', r2s(cost_p, 3), &
                        &', alignment ', r2s(cost_a, 3), ' ms.'
                endif
                write(uu,'(a)') ''
                write(uu,'(a)') 'The reference is the physical posterior''s best class (the largest in-plane marginal at the'
                write(uu,'(a)') 'physical temperature) over all classes against the teacher''s averages, as a solve2D run''s hard'
                write(uu,'(a)') 'assignment ends, from the picked positions over the same disk.'
                write(uu,'(a)') 'Class averages: every unseen particle, shortlisted or not, at its aligned class, angle and shift.'
                write(uu,'(a)') 'Costs: wall time per particle, reading and preparation included; no centring.'
                write(uu,'(a)') 'Diagnostics: '//VALIDATE_CLASS_STUDENT_DETAILS
                write(uu,'(a)') ''
                call write_timings(uu)
                end associate
            enddo
            ! the diagnostics
            write(ud,'(a)') 'CLASS STUDENT VALIDATION, DIAGNOSTICS'
            write(ud,'(a)') repeat('=', 96)
            write(ud,'(a,2x,a10,2x,a19)') lj('THE NETWORK ON THE UNSEEN PARTICLES', 52), 'picked', 'reference''s shifts'
            write(ud,'(a)') RULE
            write(ud,'(2x,a,2x,i10)') lj('Unseen particles', 50), nunseen
            write(ud,'(2x,a,2x,f10.1,2x,f19.1)') lj('Most likely class = the reference''s class (%)', 50), va(1), vi(1)
            write(ud,'(2x,a,2x,f10.1,2x,f19.1)') lj('Every list holds the reference''s class (%)', 50), va(2), vi(2)
            write(ud,'(2x,a,2x,f10.1,2x,f19.1)') lj('Particles given a shortlist (%)', 50), va(3), vi(3)
            write(ud,'(2x,a,2x,f10.1,2x,f19.1)') lj('Shortlists containing the reference''s class (%)', 50), va(4), vi(4)
            write(ud,'(2x,a,2x,i10,2x,i19)') lj('List length, every particle (median)', 50), nint(va(5)), nint(vi(5))
            write(ud,'(2x,a,2x,i10,2x,i19)') lj('Shortlist length (median)', 50), nint(va(6)), nint(vi(6))
            write(ud,'(a)') ''
            write(ud,'(a)') 'THE ALIGNMENT'
            write(ud,'(a)') RULE
            write(ud,'(2x,a,2x,f10.1)') lj('Aligned class = the reference''s class (%)', 50), av(1)
            write(ud,'(2x,a,2x,f10.3)') lj('Shift to the reference''s, same class (A, median)', 50), av(2)
            write(ud,'(2x,a,2x,f10.3)') lj('Shift to the reference''s, same class (A, 90th pct)', 50), av(3)
            write(ud,'(2x,a,2x,f10.3)') lj('Shift to the reference''s, all (A, median)', 50), av(4)
            write(ud,'(2x,a,2x,f10.3)') lj('Shift to the reference''s, all (A, 90th pct)', 50), av(5)
            write(ud,'(2x,a,2x,f10.3)') lj('Angle to the reference''s, same class (deg, median)', 50), av(6)
            write(ud,'(2x,a,2x,f10.3,2x,f10.3)') lj('Picking error (A, median and 90th pct)', 50), pick
            write(ud,'(a)') ''
            write(ud,'(a,i0,a,i0,a)') 'CLASS AVERAGES AGAINST THE TEACHER''S (', ncmp, ' classes holding ', CLASS_MINPOP_FRC, &
                &' particles in every assignment)'
            write(ud,'(a,2x,a10,2x,a10)') lj('', 52), 'student', 'reference'
            write(ud,'(a)') RULE
            write(ud,'(2x,a,2x,f10.2,2x,f10.2)') lj('FRC 0.5, whole shells (A, median)', 50), frc(1,:)
            write(ud,'(2x,a,2x,f10.2,2x,f10.2)') lj('FRC 0.5, interpolated (A, median)', 50), frc(2,:)
            write(ud,'(2x,a,2x,f10.3,2x,f10.3)') lj('Mean FRC to '//r2s(STUDENT_LP_VALIDATE, 1)//' A', 50), frc(3,:)
            write(ud,'(a)') ''
            write(ud,'(a)') 'COSTS PER PARTICLE (ms)'
            write(ud,'(a)') RULE
            write(ud,'(2x,a,2x,f10.3)') lj('Prediction, reading and preparation included', 50), cost_p
            write(ud,'(2x,a,2x,f10.3)') lj('Alignment', 50),                                     cost_a
            write(ud,'(2x,a,2x,f10.3)') lj('Full search, per particle', 50),                     cost_d
            write(ud,'(a)') ''
            write(ud,'(a)') 'A particle''s list: the student''s likely classes, a shortlist where it is confident. Picked: at the'
            write(ud,'(a)') 'picked positions; reference''s shifts: predicted again at the full search''s shifts, no alignment.'
            write(ud,'(a)') 'Same class: the particles whose aligned class is the reference''s. Picking error: the length of the'
            write(ud,'(a)') 'reference''s shift from the picked position. Costs of -1: a run over parts.'
        end subroutine report_statistics

        !> a goal's verdict
        function met( l_met ) result( str )
            logical, intent(in) :: l_met
            character(len=:), allocatable :: str
            str = 'NOT MET'
            if( l_met ) str = 'met'
        end function met

        !> x with nd decimals (1 to 3) and its leading zero, without blanks (f0.d drops the zero)
        function r2s( x, nd ) result( str )
            real,    intent(in) :: x
            integer, intent(in) :: nd
            character(len=:), allocatable :: str
            character(len=32) :: buf
            select case(nd)
                case(1)
                    write(buf,'(f16.1)') x
                case(2)
                    write(buf,'(f16.2)') x
                case default
                    write(buf,'(f16.3)') x
            end select
            str = trim(adjustl(buf))
        end function r2s

        !> str left-justified in n characters (an aw edit descriptor right-justifies)
        pure function lj( str, n ) result( r )
            character(len=*), intent(in) :: str
            integer,          intent(in) :: n
            character(len=n) :: r
            r = str
        end function lj

        !> one prediction's reads against the reference's classes, in %: the most likely class (the list's first) the
        !! reference's, the reference's class among the candidates (all, seeded), the seeded share; and the median
        !! candidate count (all, seeded; 0 without seeds)
        function arm_stats( ncand_a, route_a, cands_a ) result( v )
            integer, intent(in) :: ncand_a(nptcls), route_a(nptcls), cands_a(KMAX_CLASS,nptcls)
            real    :: v(NARM), nc(nunseen), ncs(nunseen)
            integer :: j, m, nagr, ninc, nseed, nincs
            logical :: l_in
            nagr  = 0
            ninc  = 0
            nseed = 0
            nincs = 0
            m     = 0
            do j = 1,nptcls
                if( .not. unseen(j) ) cycle
                m     = m + 1
                nc(m) = real(ncand_a(j))
                if( cands_a(1,j) == cls_d(j) ) nagr = nagr + 1
                l_in = any(cands_a(1:min(ncand_a(j),KMAX_CLASS),j) == cls_d(j))
                if( l_in ) ninc = ninc + 1
                if( route_a(j) == 1 )then
                    nseed      = nseed + 1
                    ncs(nseed) = real(ncand_a(j))
                    if( l_in ) nincs = nincs + 1
                endif
            enddo
            v(1) = 100. * real(nagr)  / real(max(1,nunseen))
            v(2) = 100. * real(ninc)  / real(max(1,nunseen))
            v(3) = 100. * real(nseed) / real(max(1,nunseen))
            v(4) = 100. * real(nincs) / real(max(1,nseed))
            v(5) = median(nc(1:m))
            v(6) = 0.
            if( nseed > 0 ) v(6) = median(ncs(1:nseed))
        end function arm_stats

        !> one alignment's reads against the reference: the aligned class the reference's (%); the shift's distance to the
        !! reference's (A) on the particles of the same class (median, 90th percentile) and on all (median, 90th
        !! percentile); the in-plane angle's distance on the same class (degrees, median)
        function aln_stats( cls_a, e3_a, sh_a ) result( v )
            integer, intent(in) :: cls_a(nptcls)
            real,    intent(in) :: e3_a(nptcls), sh_a(2,nptcls)
            real    :: v(NALN), dsame(nunseen), dall(nunseen), dang(nunseen)
            integer :: j, m, ns
            v  = 0.
            m  = 0
            ns = 0
            do j = 1,nptcls
                if( .not. unseen(j) ) cycle
                m       = m + 1
                dall(m) = sqrt(sum((sh_a(:,j) - sh_d(:,j))**2)) * params%smpd
                if( cls_a(j) /= cls_d(j) ) cycle
                ns        = ns + 1
                dsame(ns) = dall(m)
                dang(ns)  = abs(e3_a(j) - e3_d(j))
                dang(ns)  = min(dang(ns), 360. - dang(ns))
            enddo
            v(1) = 100. * real(ns) / real(max(1,nunseen))
            if( ns > 0 )then
                v(2) = median(dsame(1:ns))
                v(3) = pct90(dsame(1:ns))
                v(6) = median(dang(1:ns))
            endif
            if( m > 0 )then
                v(4) = median(dall(1:m))
                v(5) = pct90(dall(1:m))
            endif
        end function aln_stats

        !> the picking error: the reference's shift from the picked position (A; median, 90th percentile)
        function picking_error() result( v )
            real    :: v(2), d(nunseen)
            integer :: j, m
            m = 0
            do j = 1,nptcls
                if( .not. unseen(j) ) cycle
                m    = m + 1
                d(m) = sqrt(sum(sh_d(:,j)**2)) * params%smpd
            enddo
            v = 0.
            if( m > 0 ) v = [median(d(1:m)), pct90(d(1:m))]
        end function picking_error

        !> the 90th percentile of v: the sorted value at rank ceiling(0.9 n)
        function pct90( v ) result( p )
            real, intent(in) :: v(:)
            real :: p, w(size(v))
            w = v
            call hpsort(w)
            p = w(max(1, ceiling(0.9 * real(size(w)))))
        end function pct90

        !> the class averages of the student and of the reference against the teacher's, on the classes holding at least
        !! CLASS_MINPOP_FRC unseen particles in both assignments: per assignment the medians of the FRC 0.5 crossing
        !! (whole shells; interpolated between the shells about it) and the mean FRC over the shells to the band's limit
        subroutine cavgs_frc( frc_o, ncmp )
            use simple_image, only: image
            real,    intent(out) :: frc_o(3,2)
            integer, intent(out) :: ncmp
            type(image)  :: t, a
            type(string) :: stks(2)
            real, allocatable :: corrs(:), res_k(:), r05(:,:), r05i(:,:), band(:,:)
            integer :: ldim(3), ldim_a(3), n, ic, k, nf, iside, ires
            real    :: r0143, mskrad
            frc_o = -1.
            ncmp  = 0
            stks  = [student_cavgs, dense_cavgs]
            call find_ldim_nptcls(cavgs, ldim, n)
            do iside = 1,2
                call find_ldim_nptcls(stks(iside), ldim_a, n)
                if( any(ldim(1:2) /= ldim_a(1:2)) )then
                    THROW_WARN('the class averages differ in size; no FRC')
                    return
                endif
            enddo
            ldim(3) = 1
            mskrad  = 0.5 * params%mskdiam / smpd_cavgs
            allocate(r05(ncls_t,2), r05i(ncls_t,2), band(ncls_t,2), source=0.)
            call t%new(ldim, smpd_cavgs)
            call a%new(ldim, smpd_cavgs)
            nf = t%get_filtsz()
            allocate(corrs(nf), res_k(nf))
            do k = 1,nf
                res_k(k) = real(ldim(1)) * smpd_cavgs / real(k)
            enddo
            do ic = 1,ncls_t
                if( count(unseen .and. cls_s == ic) < CLASS_MINPOP_FRC ) cycle
                if( count(unseen .and. cls_d == ic) < CLASS_MINPOP_FRC ) cycle
                ncmp = ncmp + 1
                do iside = 1,2
                    call t%read(cavgs, ic)
                    call a%read(stks(iside), ic)
                    call t%mask2D_soft(mskrad, backgr=0.)
                    call a%mask2D_soft(mskrad, backgr=0.)
                    call t%fft
                    call a%fft
                    call t%fsc(a, corrs)
                    call get_resolution(corrs, res_k, r05(ncmp,iside), r0143)
                    ! the crossing between the last shell at or above 0.5 (get_resolution's) and the next
                    r05i(ncmp,iside) = r05(ncmp,iside)
                    ires = 0
                    if( r05(ncmp,iside) > 0. ) ires = minloc(abs(res_k - r05(ncmp,iside)), dim=1)
                    if( ires == 0 )then
                        r05i(ncmp,iside) = res_k(1)
                    else if( ires < nf )then
                        r05i(ncmp,iside) = res_k(1) / (real(ires) + (corrs(ires) - 0.5) / max(corrs(ires) - corrs(ires+1), 1.e-9))
                    endif
                    band(ncmp,iside) = sum(corrs, mask=res_k >= STUDENT_LP_VALIDATE) / real(max(1, count(res_k >= STUDENT_LP_VALIDATE)))
                enddo
            enddo
            if( ncmp > 0 )then
                do iside = 1,2
                    frc_o(1,iside) = median(r05(1:ncmp,iside))
                    frc_o(2,iside) = median(r05i(1:ncmp,iside))
                    frc_o(3,iside) = sum(band(1:ncmp,iside)) / real(ncmp)
                enddo
            endif
            call t%kill
            call a%kill
        end subroutine cavgs_frc

        !> the phases' wall times and their shares, in the run's order, to unit uu
        subroutine write_timings( uu )
            integer, intent(in) :: uu
            integer :: j, k
            write(uu,'(a,2x,a10,3x,a)') lj('RUN TIME', 52), 'seconds', 'share'
            write(uu,'(a)') repeat('-', 96)
            do k = 1,NPHASES
                j = PHASE_ORDER(k)
                write(uu,'(2x,a,2x,f10.2,3x,f5.1,a)') lj(PHASE_NAMES(j), 50), real(rt_ph(j)), &
                    &100. * real(rt_ph(j)) / max(real(rt_all), 1.e-9), ' %'
            enddo
            write(uu,'(2x,a,2x,f10.2)') lj('Total', 50), real(rt_all)
        end subroutine write_timings

    end subroutine exec_validate_class_student

    !> the classes are the class averages in refs
    subroutine set_ncls_from_refs( cline )
        use simple_error,   only: simple_exception
        use simple_string,  only: string
        use simple_imghead, only: find_ldim_nptcls
        class(cmdline), intent(inout) :: cline
        type(string) :: refs
        integer      :: ldim(3), n
        if( .not. cline%defined('refs') ) THROW_HARD('the class averages are refs; learn_class_student')
        refs = cline%get_carg('refs')
        call find_ldim_nptcls(refs, ldim, n)
        if( cline%defined('ncls') )then
            if( cline%get_iarg('ncls') /= n ) THROW_HARD('ncls differs from the class averages in refs; learn_class_student')
        endif
        call cline%set('ncls', n)
        call refs%kill
    end subroutine set_ncls_from_refs

end module simple_commanders_class_student
