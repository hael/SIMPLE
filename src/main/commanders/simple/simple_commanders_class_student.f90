!@descr: the class student's commanders: learn_class_student (the teacher's rows against frozen class averages and the training), predict_classes (the class, the angle, the candidates and the confidence of unseen particles) and the validate_class_student controller
module simple_commanders_class_student
use simple_commander_base,         only: commander_base
use simple_cmdline,                only: cmdline
use simple_strategy_class_student, only: run_learn_class_student, run_predict_classes
use simple_defs_hhcnn,             only: KMAX_CLASS, STUDENT_SHIFT_SIGMA, STUDENT_TRS, STUDENT_LP_VALIDATE, CLASS_DENSE_SHIFT_SIGMA, &
    &CLASS_DENSE_PROB_BETA, CLASS_MINPOP_FRC, STUDENT_TRAIN_PROJ, STUDENT_ALL_PROJ, STUDENT_UNSEEN_PROJ, &
    &CLASS_STUDENT_MODEL_FNAME, CLASS_STUDENT_PRED_TABLE, CLASS_STUDENT_CAVGS, VALIDATE_CLASS_STUDENT_REPORT, VALIDATE_CLASS_STUDENT_DETAILS
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
        call run_learn_class_student(cline)
    end subroutine exec_learn_class_student

    !> The student's prediction for the project's selected particles under the model infile: per particle the class,
    !! its in-plane angle, the candidates holding 99% of the posterior and the confidence, and, given the class averages
    !! refs, the shift registered against the predicted class (simple_private_exec). With nparts it is the master of
    !! that many parts; a part is given its range as the master hands it (part, nparts, fromp, top, outfile)
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
        endif
        call run_predict_classes(cline)
    end subroutine exec_predict_classes

    !> validate_class_student: solve2D on a random nsample (the teacher), the student trained, center_ptcls on all, the
    !! unseen predicted and searched densely against the teacher's averages, the student's class averages, the reads
    !! and the timings (doc/algorithms/class_student.md). Every step but the training runs over nparts when given
    subroutine exec_validate_class_student( self, cline )
        use simple_core_module_api
        use simple_parameters,              only: parameters
        use simple_sp_project,              only: sp_project
        use simple_commanders_solve2D,      only: commander_solve2D
        use simple_commanders_refine2D,     only: commander_refine2D
        use simple_commanders_center_ptcls, only: commander_center_ptcls
        use simple_commanders_mkcavgs,      only: commander_make_cavgs
        use simple_refine2D_fnames,         only: refine2D_cavgs_fname
        class(commander_validate_class_student), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        integer,          parameter :: NPHASES = 10   !< the validation's timed phases
        integer,          parameter :: NARM    = 6    !< arm_stats' reads
        character(len=*), parameter :: PHASE_NAMES(NPHASES) = [character(len=48) :: &
            &'Drawing the training set', 'Teacher: solve2D on the training set', &
            &'Training the student (learn_class_student)', 'Centring every particle (center_ptcls)', &
            &'Student prediction (predict_classes)', 'Reference: best class of the unseen (refine2D)', &
            &'Class averages of the student''s classes', 'Statistics', 'Project copies and bookkeeping', &
            &'Prediction at ideal shifts (details)']
        integer, parameter :: PHASE_ORDER(NPHASES) = [1, 2, 3, 4, 5, 6, 10, 7, 8, 9] !< the run's order
        type(commander_solve2D)             :: xsolve2D
        type(commander_refine2D)            :: xrefine2D
        type(commander_center_ptcls)        :: xcenter
        type(commander_make_cavgs)          :: xmkcavgs
        type(commander_learn_class_student) :: xlearn
        type(commander_predict_classes)     :: xpredict
        type(parameters)  :: params
        type(sp_project)  :: spproj, res, work
        type(cmdline)     :: cline_entry, cl
        type(ran_tabu)    :: rt
        type(string)      :: edir, result_proj, cavgs, sigma2_abs, model_abs, fname, dense_cavgs
        integer, allocatable :: pinds(:), cls_t(:), cls_s(:), cls_d(:), ncand(:), route(:), cands(:,:)
        integer, allocatable :: cls_s2(:), ncand2(:), route2(:), cands2(:,:)
        logical, allocatable :: active(:), train(:), unseen(:)
        real,    allocatable :: sh_t(:,:), sh_c(:,:), sh_s(:,:), sh_d(:,:), e3_s(:), e3_d(:)
        real(timer_int_kind) :: rt_ph(NPHASES), rt_all
        integer(timer_int_kind) :: t0, t_all
        real    :: smpd_cavgs
        integer :: nptcls, nactive, ntrain, nunseen, ncls_t, i, funit, dunit, io_stat
        logical :: l_sigma2
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
        rt_ph(9) = rt_ph(9) + toc(t0)
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
        ! ---- the student: trained on the teacher's rows against its class averages
        t0 = tic()
        call res%read(result_proj)
        call res%get_cavgs_stk(cavgs, ncls_t, smpd_cavgs)
        cavgs = simple_abspath(cavgs)
        call res%get_sigma2_state_path(sigma2_abs, l_sigma2)
        if( l_sigma2 ) l_sigma2 = file_exists(sigma2_abs)
        if( l_sigma2 ) sigma2_abs = simple_abspath(sigma2_abs)
        allocate(cls_t(nptcls), cls_s(nptcls), cls_d(nptcls), ncand(nptcls), route(nptcls), source=0)
        allocate(cands(KMAX_CLASS,nptcls), source=0)
        allocate(cls_s2(nptcls), ncand2(nptcls), route2(nptcls), source=0)
        allocate(cands2(KMAX_CLASS,nptcls), source=0)
        allocate(sh_t(2,nptcls), sh_c(2,nptcls), sh_s(2,nptcls), sh_d(2,nptcls), source=0.)
        allocate(e3_s(nptcls), e3_d(nptcls), source=0.)
        do i = 1,nptcls
            if( .not. train(i) ) cycle
            cls_t(i)  = res%os_ptcl2D%get_class(i)
            sh_t(:,i) = res%os_ptcl2D%get_2Dshift(i)
        enddo
        call enter_subdir('learn')
        call write_copy(res, STUDENT_TRAIN_PROJ)
        rt_ph(9) = rt_ph(9) + toc(t0)
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
        ! ---- every active particle centred from its stored shift, as raw picks are before labelling (decision 9)
        t0 = tic()
        call enter_subdir('center')
        call select_rows(res, active)
        call write_copy(res, STUDENT_ALL_PROJ)
        rt_ph(9) = rt_ph(9) + toc(t0)
        write(logfhandle,'(A)') '>>> VALIDATE_CLASS_STUDENT: CENTRING EVERY PARTICLE'
        t0 = tic()
        call cl%set('prg',      'center_ptcls')
        call cl%set('projfile', STUDENT_ALL_PROJ)
        call cl%set('mskdiam',  params%mskdiam)
        call cl%set('nthr',     params%nthr)
        call cl%set('mkdir',    'no')
        call set_parts(cl)
        call xcenter%execute(cl)
        call cl%kill
        call work%read(string(STUDENT_ALL_PROJ))
        do i = 1,nptcls
            if( active(i) ) sh_c(:,i) = work%os_ptcl2D%get_2Dshift(i)
        enddo
        call simple_chdir(edir)
        rt_ph(4) = toc(t0)
        ! ---- the unseen particles at their centred shifts: the student's prediction
        t0 = tic()
        call work%os_ptcl2D%delete_2Dclustering(keepshifts=.true.)
        call select_rows(work, unseen)
        call enter_subdir('predict')
        call write_copy(work, STUDENT_UNSEEN_PROJ)
        rt_ph(9) = rt_ph(9) + toc(t0)
        write(logfhandle,'(A)') '>>> VALIDATE_CLASS_STUDENT: STUDENT PREDICTION FOR THE UNSEEN PARTICLES'
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
        call read_pred_table(string(CLASS_STUDENT_PRED_TABLE), ncand, route, cands)
        call res%kill
        call res%read(string(STUDENT_UNSEEN_PROJ))
        do i = 1,nptcls
            if( .not. unseen(i) ) cycle
            cls_s(i)  = res%os_ptcl2D%get_class(i)
            e3_s(i)   = res%os_ptcl2D%e3get(i)
            sh_s(:,i) = res%os_ptcl2D%get_2Dshift(i)
        enddo
        call simple_chdir(edir)
        rt_ph(5) = toc(t0)
        ! ---- the search the student replaces: one dense pass of the same particles against the teacher's averages, as
        !      the first iteration of its own run (every particle searches the disk), at the student's box and band
        t0 = tic()
        call enter_subdir('dense')
        call write_copy(work, STUDENT_UNSEEN_PROJ)
        rt_ph(9) = rt_ph(9) + toc(t0)
        write(logfhandle,'(A)') '>>> VALIDATE_CLASS_STUDENT: DENSE PASS OF THE UNSEEN PARTICLES'
        t0 = tic()
        call cl%set('prg',         'refine2D')
        call cl%set('projfile',    STUDENT_UNSEEN_PROJ)
        call cl%set('oritype',     'ptcl2D')
        call cl%set('ncls',        ncls_t)
        call cl%set('refs',        cavgs%to_char())
        call cl%set('refine',      'prob')
        call cl%set('objfun',      'euclid')
        call cl%set('lp',          STUDENT_LP_VALIDATE)
        call cl%set('shift_sigma', CLASS_DENSE_SHIFT_SIGMA)
        call cl%set('prob_beta',   CLASS_DENSE_PROB_BETA)
        call cl%set('maxits',      1)
        call cl%set('center',      'no')
        call cl%set('mskdiam',     params%mskdiam)
        call cl%set('nthr',        params%nthr)
        call cl%set('mkdir',       'no')
        call set_parts(cl)
        call xrefine2D%execute(cl)
        call cl%kill
        call res%kill
        call res%read(string(STUDENT_UNSEEN_PROJ))
        do i = 1,nptcls
            if( .not. unseen(i) ) cycle
            cls_d(i)  = res%os_ptcl2D%get_class(i)
            e3_d(i)   = res%os_ptcl2D%e3get(i)
            sh_d(:,i) = res%os_ptcl2D%get_2Dshift(i)
        enddo
        dense_cavgs = simple_abspath(refine2D_cavgs_fname(1, params%ext%to_char()))
        call simple_chdir(edir)
        rt_ph(6) = toc(t0)
        ! ---- the centring's ceiling: the student's prediction of the same particles at the dense pass's shifts
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
        call read_pred_table(string(CLASS_STUDENT_PRED_TABLE), ncand2, route2, cands2)
        call res%kill
        call res%read(string(STUDENT_UNSEEN_PROJ))
        do i = 1,nptcls
            if( unseen(i) ) cls_s2(i) = res%os_ptcl2D%get_class(i)
        enddo
        call simple_chdir(edir)
        rt_ph(10) = toc(t0)
        ! ---- the class averages of the student's assignments alone (no search)
        t0 = tic()
        call enter_subdir('cavgs')
        call res%kill
        call res%read(string(edir%to_char()//'/predict/'//STUDENT_UNSEEN_PROJ))
        call write_copy(res, STUDENT_UNSEEN_PROJ)
        call cl%set('prg',      'make_cavgs')
        call cl%set('projfile', STUDENT_UNSEEN_PROJ)
        call cl%set('ncls',     ncls_t)
        call cl%set('refs',     CLASS_STUDENT_CAVGS)
        call cl%set('mskdiam',  params%mskdiam)
        call cl%set('nthr',     params%nthr)
        call cl%set('mkdir',    'no')
        call set_parts(cl)
        call xmkcavgs%execute(cl)
        call cl%kill
        call simple_chdir(edir)
        rt_ph(7) = toc(t0)
        ! ---- the statistics and the report
        t0 = tic()
        call fopen(funit, string(VALIDATE_CLASS_STUDENT_REPORT), status='replace', action='write', iostat=io_stat)
        call fileiochk('validate_class_student; open '//VALIDATE_CLASS_STUDENT_REPORT, io_stat)
        call fopen(dunit, string(VALIDATE_CLASS_STUDENT_DETAILS), status='replace', action='write', iostat=io_stat)
        call fileiochk('validate_class_student; open '//VALIDATE_CLASS_STUDENT_DETAILS, io_stat)
        rt_ph(8) = toc(t0)
        rt_all = toc(t_all)
        call report_statistics(funit, dunit)
        call fclose(funit)
        call fclose(dunit)
        call cline_entry%kill
        call spproj%kill
        call res%kill
        call work%kill
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

        !> the verdict (can the student replace the full search?) with how to read it and the run time, to unit u and the
        !! log; the diagnostics to unit ud (doc/algorithms/class_student.md, section Validation)
        subroutine report_statistics( u, ud )
            use simple_class_student, only: class_student
            integer, intent(in) :: u, ud
            character(len=96), parameter :: RULE = repeat('-', 96)
            type(class_student) :: stu
            real, allocatable :: dang(:), dsh(:), dsh_c(:), dfx(:), dfy(:)
            real    :: off(2), r_s, r_d, thr, target, hs, hr, cost_s, cost_d, cost_c, cost_p, va(NARM), vi(NARM)
            integer :: j, iu, units(2), nagree, n, ntr, ncmp
            character(len=:), allocatable :: flag
            units = [u, logfhandle]
            allocate(dang(nptcls), dsh(nptcls), dsh_c(nptcls), dfx(nptcls), dfy(nptcls), source=0.)
            nagree = 0
            n      = 0
            do j = 1,nptcls
                if( .not. unseen(j) ) cycle
                n = n + 1
                dsh_c(n) = sqrt(sum((sh_c(:,j) - sh_d(:,j))**2)) * params%smpd
                if( cls_s(j) /= cls_d(j) ) cycle
                nagree       = nagree + 1
                dang(nagree) = abs(e3_s(j) - e3_d(j))
                dang(nagree) = min(dang(nagree), 360. - dang(nagree))
                dsh(nagree)  = sqrt(sum((sh_s(:,j) - sh_d(:,j))**2)) * params%smpd
            enddo
            ! the frame offset: center_ptcls's shifts against the teacher's on the training particles, their common
            ! translation (the median per axis) and the residual after it
            ntr = 0
            do j = 1,nptcls
                if( .not. train(j) ) cycle
                ntr = ntr + 1
                dfx(ntr) = (sh_c(1,j) - sh_t(1,j)) * params%smpd
                dfy(ntr) = (sh_c(2,j) - sh_t(2,j)) * params%smpd
            enddo
            off(1) = median(dfx(1:ntr))
            off(2) = median(dfy(1:ntr))
            do j = 1,ntr
                dfx(j) = sqrt((dfx(j) - off(1))**2 + (dfy(j) - off(2))**2)
            enddo
            fname = edir%to_char()//'/cavgs/'//CLASS_STUDENT_CAVGS
            call cavgs_frc(fname, r_s, r_d, ncmp)
            va = arm_stats(cls_s,  ncand,  route,  cands)
            vi = arm_stats(cls_s2, ncand2, route2, cands2)
            call stu%read(model_abs)
            call stu%get_threshold(thr, target, hs, hr)
            call stu%kill
            target = 100. * target
            cost_c = 1000. * real(rt_ph(4)) / real(max(1,nactive))
            cost_p = 1000. * real(rt_ph(5)) / real(max(1,nunseen))
            cost_s = cost_c + cost_p
            cost_d = 1000. * real(rt_ph(6)) / real(max(1,nunseen))
            do iu = 1,2
                associate( uu => units(iu) )
                write(uu,'(a)') 'CLASS STUDENT VALIDATION'
                write(uu,'(a)') repeat('=', 96)
                write(uu,'(a,i0,a,i0,a)') 'The student learns the classes of a solve2D run on ', ntrain, ' particles (', ncls_t, &
                    &' classes). For each'
                write(uu,'(a)') 'new particle it proposes a shortlist of likely classes; where it is confident, alignment tests'
                write(uu,'(a,i0,a)') 'only the shortlist. Checked on ', nunseen, ' particles it has not seen, against a full search.'
                write(uu,'(a)') ''
                write(uu,'(a,a10,3x,a)') lj('', 54), 'result', 'goal'
                write(uu,'(a)') RULE
                write(uu,'(2x,a,2x,f10.1,3x,a)') lj('Particles given a shortlist (%)', 50), va(3), 'as many as possible'
                write(uu,'(2x,a,2x,f10.1,3x,a,a)') lj('Shortlists containing the full search''s class (%)', 50), va(4), &
                    &lj(r2s(target, 1)//' or more', 26), met(va(4) >= target)
                write(uu,'(2x,a,2x,f10.1,3x,a,a)') lj('  same, with the particles perfectly centred (%)', 50), vi(4), &
                    &lj(r2s(target, 1)//' or more', 26), met(vi(4) >= target)
                if( va(6) > 0. )then
                    write(uu,'(2x,a,2x,i10,3x,a,i0,a)') lj('Shortlist length (median)', 50), nint(va(6)), 'of ', ncls_t, ' classes'
                else
                    write(uu,'(2x,a,2x,a10,3x,a,i0,a)') lj('Shortlist length (median)', 50), '-', 'of ', ncls_t, ' classes'
                endif
                if( ncmp > 0 )then
                    write(uu,'(2x,a,2x,f10.2,3x,a)') lj('Class averages vs the teacher''s, FRC 0.5 (A)', 50), r_s, &
                        &'full search '//r2s(r_d, 2)//'; lower is better'
                else
                    write(uu,'(2x,a,2x,a10,3x,a)') lj('Class averages vs the teacher''s, FRC 0.5 (A)', 50), '-', 'no class compared'
                endif
                write(uu,'(2x,a,2x,f10.3,3x,a,a)') lj('Cost per particle (ms)', 50), cost_s, &
                    &lj('below full search '//r2s(cost_d, 3), 26), met(cost_s < cost_d)
                write(uu,'(a)') ''
                if( va(4) >= target )then
                    write(uu,'(a)') 'The shortlists are safe: they hold the full search''s class as often as the goal asks.'
                else if( vi(4) >= target )then
                    write(uu,'(a)') 'The shortlists miss the full search''s class because of the centring: with the particles'
                    write(uu,'(a)') 'perfectly centred they meet the goal.'
                else
                    write(uu,'(a)') 'The shortlists miss the full search''s class even with the particles perfectly centred:'
                    write(uu,'(a)') 'the student itself falls short.'
                endif
                if( cost_s >= cost_d )then
                    write(uu,'(a,a,a,a,a)') 'The student costs more than the full search; centring takes ', r2s(cost_c, 3), &
                        &' of its ', r2s(cost_s, 3), ' ms.'
                endif
                write(uu,'(a)') ''
                write(uu,'(a)') 'Full search: refine2D refine=prob, one iteration over all classes against the teacher''s averages.'
                write(uu,'(a)') 'The student''s cost includes centring (center_ptcls). Diagnostics: '//VALIDATE_CLASS_STUDENT_DETAILS
                write(uu,'(a)') ''
                call write_timings(uu)
                end associate
            enddo
            ! the diagnostics
            write(ud,'(a)') 'CLASS STUDENT VALIDATION, DIAGNOSTICS'
            write(ud,'(a)') repeat('=', 96)
            write(ud,'(a,2x,a10,2x,a19)') lj('THE UNSEEN PARTICLES', 52), 'centred', 'perfectly centred'
            write(ud,'(a)') RULE
            write(ud,'(2x,a,2x,i10)') lj('Unseen particles', 50), nunseen
            write(ud,'(2x,a,2x,f10.1,2x,f19.1)') lj('Student''s class = full search''s class (%)', 50), va(1), vi(1)
            write(ud,'(2x,a,2x,f10.1,2x,f19.1)') lj('Every particle''s list has the full-search class (%)', 50), va(2), vi(2)
            write(ud,'(2x,a,2x,f10.1,2x,f19.1)') lj('Particles given a shortlist (%)', 50), va(3), vi(3)
            write(ud,'(2x,a,2x,f10.1,2x,f19.1)') lj('Shortlists containing the full search''s class (%)', 50), va(4), vi(4)
            write(ud,'(2x,a,2x,i10,2x,i19)') lj('List length, every particle (median)', 50), nint(va(5)), nint(vi(5))
            write(ud,'(2x,a,2x,i10,2x,i19)') lj('Shortlist length (median)', 50), nint(va(6)), nint(vi(6))
            if( nagree > 0 )then
                write(ud,'(2x,a,2x,f10.3)') lj('In-plane angle error, same class (deg, median)', 50), median(dang(1:nagree))
                write(ud,'(2x,a,2x,f10.3)') lj('Registered shift error, same class (A, median)', 50), median(dsh(1:nagree))
            endif
            write(ud,'(2x,a,2x,f10.3)') lj('Centring error before registration (A, median)', 50), median(dsh_c(1:n))
            write(ud,'(a)') ''
            write(ud,'(a)') 'CENTER_PTCLS AGAINST THE TEACHER (training particles)'
            write(ud,'(a)') RULE
            write(ud,'(2x,a,2x,2f10.3)') lj('Common offset, x and y (A)', 50),                 off
            write(ud,'(2x,a,2x,f10.3)') lj('Remaining difference after it (A, median)', 50), median(dfx(1:ntr))
            write(ud,'(a)') ''
            write(ud,'(a)') 'COSTS PER PARTICLE (ms)'
            write(ud,'(a)') RULE
            write(ud,'(2x,a,2x,f10.3)') lj('Centring, per particle', 50),              cost_c
            write(ud,'(2x,a,2x,f10.3)') lj('Student prediction, per particle', 50),    cost_p
            write(ud,'(2x,a,2x,f10.3)') lj('Full search, per particle', 50),           cost_d
            write(ud,'(a)') ''
            write(ud,'(a)') 'A particle''s list: the student''s likely classes, a shortlist where it is confident. Its class:'
            write(ud,'(a)') 'the most probable. Perfectly centred: predicted again at the full search''s shifts.'
            deallocate(dang, dsh, dsh_c, dfx, dfy)
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

        !> one prediction's reads against the full search's classes, in %: the same class, the full-search class among the
        !! candidates (all, seeded), the seeded share; and the median candidate count (all, seeded; 0 without seeds)
        function arm_stats( cls_a, ncand_a, route_a, cands_a ) result( v )
            integer, intent(in) :: cls_a(nptcls), ncand_a(nptcls), route_a(nptcls), cands_a(KMAX_CLASS,nptcls)
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
                if( cls_a(j) == cls_d(j) ) nagr = nagr + 1
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

        !> per class populated by at least CLASS_MINPOP_FRC unseen particles on both routes, the FRC 0.5 resolution of the
        !! student's and of the dense pass's class average against the teacher's; their medians
        subroutine cavgs_frc( student_stk, r_s, r_d, ncmp )
            use simple_image, only: image
            class(string), intent(in)  :: student_stk
            real,          intent(out) :: r_s, r_d
            integer,       intent(out) :: ncmp
            type(image) :: t, a
            real, allocatable :: corrs(:), res_k(:), rs(:), rd(:)
            integer :: ldim(3), ldim_s(3), ldim_d(3), n, ic, k, nf, iside
            real    :: r05, r0143, mskrad
            r_s  = -1.
            r_d  = -1.
            ncmp = 0
            call find_ldim_nptcls(cavgs, ldim, n)
            call find_ldim_nptcls(student_stk, ldim_s, n)
            call find_ldim_nptcls(dense_cavgs, ldim_d, n)
            if( any(ldim(1:2) /= ldim_s(1:2)) .or. any(ldim(1:2) /= ldim_d(1:2)) )then
                THROW_WARN('the class averages differ in size; no FRC')
                return
            endif
            ldim(3) = 1
            mskrad  = 0.5 * params%mskdiam / smpd_cavgs
            allocate(rs(ncls_t), rd(ncls_t), source=0.)
            call t%new(ldim, smpd_cavgs)
            call a%new(ldim, smpd_cavgs)
            nf = t%get_filtsz()
            allocate(corrs(nf), res_k(nf))
            do k = 1,nf
                res_k(k) = real(ldim(1)) * smpd_cavgs / real(k)
            enddo
            do ic = 1,ncls_t
                if( count(unseen .and. cls_s == ic) < CLASS_MINPOP_FRC .or. count(unseen .and. cls_d == ic) < CLASS_MINPOP_FRC ) cycle
                ncmp = ncmp + 1
                do iside = 1,2
                    call t%read(cavgs, ic)
                    if( iside == 1 )then
                        call a%read(student_stk, ic)
                    else
                        call a%read(dense_cavgs, ic)
                    endif
                    call t%mask2D_soft(mskrad, backgr=0.)
                    call a%mask2D_soft(mskrad, backgr=0.)
                    call t%fft
                    call a%fft
                    call t%fsc(a, corrs)
                    call get_resolution(corrs, res_k, r05, r0143)
                    if( iside == 1 )then
                        rs(ncmp) = r05
                    else
                        rd(ncmp) = r05
                    endif
                enddo
            enddo
            if( ncmp > 0 )then
                r_s = median(rs(1:ncmp))
                r_d = median(rd(1:ncmp))
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
