!@descr: the pose student's commanders: learn_pose_student (the teacher's rows and the training), predict_poses (the mode, the cone and the confidence of unseen particles) and the validate_pose_student controller
module simple_commanders_pose_student
use simple_commander_base,        only: commander_base
use simple_cmdline,               only: cmdline
use simple_error,                 only: simple_exception
use simple_strategy_pose_student, only: run_learn_pose_student, run_predict_poses, run_pose_posterior_mode
use simple_strategy_student_common, only: fname_base
use simple_defs_hhcnn,            only: NSPACE_STUDENT, KMAX_STUDENT, STUDENT_TRS, STUDENT_SHIFT_SIGMA, STUDENT_LP_VALIDATE, &
    &POSE_AGREE_DEG, POSE_TWO_STEPS_DEG, STUDENT_TRAIN_PROJ, STUDENT_UNSEEN_PROJ, POSE_STUDENT_MODEL_FNAME, &
    &POSE_STUDENT_PRED_TABLE, POSE_STUDENT_PRED_TIMES, VALIDATE_POSE_STUDENT_REPORT, VALIDATE_POSE_STUDENT_DETAILS, &
    &POSE_STUDENT_PROB_BETA, POSE_DENSE_SHIFT_SIGMA, POSE_REFERENCE_FNAME, &
    &POSE_REFERENCE_MODEL_FNAME
implicit none
private
public :: commander_learn_pose_student, commander_predict_poses, commander_validate_pose_student
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_learn_pose_student
  contains
    procedure :: execute => exec_learn_pose_student
end type commander_learn_pose_student

type, extends(commander_base) :: commander_predict_poses
  contains
    procedure :: execute => exec_predict_poses
end type commander_predict_poses

type, extends(commander_base) :: commander_validate_pose_student
  contains
    procedure :: execute => exec_validate_pose_student
end type commander_validate_pose_student


contains

    !> The teacher's rows of the project's selected particles at the student grid (nspace) and band (lp) against the
    !! references vol1, with their resolvability and their features, and the student trained on them, written as the
    !! model fname with its weights beside it (simple_private_exec). The dense table's mode and objective are set
    !! here: the student has no key of its own
    subroutine exec_learn_pose_student( self, cline )
        class(commander_learn_pose_student), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',  'no')
        if( .not. cline%defined('nspace') ) call cline%set('nspace', NSPACE_STUDENT)
        if( .not. cline%defined('trs')    ) call cline%set('trs',    STUDENT_TRS)
        if( .not. cline%defined('fname')  ) call cline%set('fname',  POSE_STUDENT_MODEL_FNAME)
        call cline%set('oritype', 'ptcl3D')
        call cline%set('objfun',  'euclid')
        call cline%set('refine',  'prob')
        call cline%set('center',  'no')
        ! the teacher is the physical posterior a solve3D run ends at: no tempering
        call cline%set('prob_beta', POSE_STUDENT_PROB_BETA)
        call run_learn_pose_student(cline)
    end subroutine exec_learn_pose_student

    !> The student's prediction for the project's selected particles under the model infile: per particle the mode
    !! direction, the cone holding 99% of the posterior and the confidence, and, given the map vol1, the alignment over
    !! the reference pass's disk, to the text tables beside fname (simple_private_exec). With nparts it is the master of
    !! that many parts; a part is given its range as the master hands it (part, nparts, fromp, top, outfile)
    subroutine exec_predict_poses( self, cline )
        class(commander_predict_poses), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        if( cline%defined('part') )then
            if( .not. (cline%defined('nparts') .and. cline%defined('fromp') .and. cline%defined('top')) ) &
                &THROW_HARD('a part of predict_poses needs its range: nparts, fromp and top')
            if( .not. cline%defined('outfile') ) THROW_HARD('a part of predict_poses needs outfile')
        endif
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'no')
        if( .not. cline%defined('fname') ) call cline%set('fname', POSE_STUDENT_MODEL_FNAME)
        call cline%set('oritype', 'ptcl3D')
        call cline%set('objfun',  'euclid')
        ! the map is the model's frame: never re-centred (centring would move every particle's shift by its offset)
        call cline%set('center',  'no')
        if( cline%defined('vol1') )then
            ! the alignment's disk and prior are the reference pass's, at the physical temperature
            call cline%set('refine', 'prob')
            if( .not. cline%defined('trs')         ) call cline%set('trs',         STUDENT_TRS)
            if( .not. cline%defined('shift_sigma') ) call cline%set('shift_sigma', POSE_DENSE_SHIFT_SIGMA)
            call cline%set('prob_beta', POSE_STUDENT_PROB_BETA)
        endif
        call run_predict_poses(cline)
    end subroutine exec_predict_poses

    !> The physical posterior's best pose of the project's selected particles against vol1 (the validation's reference),
    !! the dense table at the student's grid and the band lp, from the stored shifts over the reference pass's disk
    subroutine pose_posterior_mode( cline )
        class(cmdline), intent(inout) :: cline
        call cline%set('mkdir',       'no')
        call cline%set('nspace',      NSPACE_STUDENT)
        call cline%set('trs',         STUDENT_TRS)
        call cline%set('shift_sigma', POSE_DENSE_SHIFT_SIGMA)
        call cline%set('prob_beta',   POSE_STUDENT_PROB_BETA)
        call cline%set('fname',       POSE_REFERENCE_MODEL_FNAME)
        call cline%set('oritype',     'ptcl3D')
        call cline%set('objfun',      'euclid')
        call cline%set('refine',      'prob')
        call cline%set('center',      'no')
        call run_pose_posterior_mode(cline)
    end subroutine pose_posterior_mode


    !> validate_pose_student: a class-balanced greedy sample of nsample particles of the 2D solution, solve3D on all of
    !! it (the teacher; or the frozen teacher of a solve3D directory dir_refine run on the same sample), the student
    !! trained on its poses and map, the other selected particles predicted at their 2D shifts (centred by the prior
    !! solve2D, as the training particles) and aligned after the prediction, judged against the physical posterior's best
    !! pose (doc/algorithms/pose_student.md)
    subroutine exec_validate_pose_student( self, cline )
        use simple_core_module_api
        use simple_parameters,              only: parameters
        use simple_sp_project,              only: sp_project
        use simple_commanders_solve3D,      only: commander_solve3D
        use simple_view_partition_sampling, only: make_class_samples
        use simple_bench_io,                only: open_bench_report
        use simple_commanders_rec,          only: commander_rec3D
        class(commander_validate_pose_student), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        integer,          parameter :: NPHASES = 8   !< the validation's timed phases
        integer,          parameter :: NRECS   = 2   !< the unseen particles' maps
        character(len=*), parameter :: PHASE_NAMES(NPHASES) = [character(len=48) :: &
            &'class-balanced training sample', 'solve3D on the training set (teacher)', &
            &'student training (learn_pose_student)', 'reference and dense search (the dense table)', &
            &'prediction and alignment (predict_poses)', 'statistics', 'reconstructions from the poses', &
            &'project copies and bookkeeping']
        character(len=*), parameter :: REC_NAMES(NRECS) = [character(len=32) :: 'unseen_aligned_cone', &
            &'unseen_reference']
        character(len=*), parameter :: REC_HELP(NRECS) = [character(len=80) :: &
            &'the unseen particles at the alignment over the cone''s poses and shifts (SK)', &
            &'the unseen particles at the reference''s poses and shifts']
        type(commander_solve3D)            :: xsolve3D
        type(commander_learn_pose_student) :: xlearn
        type(commander_predict_poses)      :: xpredict
        type(commander_rec3D)              :: xrec3D
        type(parameters)                   :: params
        type(sp_project)                   :: spproj, res, tproj
        type(cmdline)                      :: cline_entry, cl
        type(class_sample), allocatable    :: clssmp(:)
        type(string) :: edir, result_proj, vol_fname, model_abs, table_abs, times_abs, fname, teacher_dir
        integer, allocatable :: states(:)
        logical, allocatable :: active(:), train(:), unseen(:), has_ref(:)
        real,    allocatable :: sh_2D(:,:), euls_ref(:,:), sh_ref(:,:)
        real,    allocatable :: euls_sk(:,:), sh_sk(:,:)
        real(dp) :: ref_cost(2)               ! the reference pass per particle (s): the table fill; setup, particles and fill
        real     :: rec_res(4,NRECS)          ! per map: FSC 0.5 and 0.143 against the teacher's map, then of its halves
        real     :: times(2)                  ! the prediction's and the alignment's walls (s)
        integer  :: rec_n(NRECS)
        real(timer_int_kind) :: rt(NPHASES), rt_all
        integer(timer_int_kind) :: t0, t_all
        real    :: smpd_vol
        integer :: nptcls, navail, ntrain, nactive, nunseen, i, box_vol, funit, dunit, fnr, io_stat, ncls2D
        logical :: l_frozen
        t_all = tic()
        rt    = 0.
        ! ---- the command line: one state, an execution directory, the training sample's size
        if( .not. cline%defined('nsample') ) THROW_HARD('validate_pose_student needs nsample, the training sample size')
        if( cline%defined('nstates') )then
            if( cline%get_iarg('nstates') /= 1 ) THROW_HARD('the pose student is single-state: nstates=1')
        endif
        if( cline%defined('mkdir') )then
            if( .not. (cline%get_carg('mkdir') == 'yes') ) THROW_HARD('validate_pose_student needs mkdir=yes: its steps live in its execution directory')
        endif
        call cline%set('mkdir', 'yes')
        l_frozen = cline%defined('dir_refine')
        cline_entry = cline
        call params%new(cline)
        call simple_getcwd(edir)
        if( l_frozen ) teacher_dir = simple_abspath(params%dir_refine)
        ! ---- the training sample: class-balanced and greedy (best 2D scores first) over the selected classes
        t0 = tic()
        call spproj%read(params%projfile)
        nptcls = spproj%os_ptcl2D%get_noris()
        if( nptcls < 1 ) THROW_HARD('the project has no particles (ptcl2D)')
        if( spproj%os_ptcl3D%get_noris() /= nptcls ) THROW_HARD('ptcl2D and ptcl3D differ in size')
        if( spproj%os_cls2D%get_noris() < 1 ) THROW_HARD('validate_pose_student needs a 2D solution (cls2D)')
        allocate(active(nptcls), train(nptcls), unseen(nptcls), has_ref(nptcls), source=.false.)
        do i = 1,nptcls
            active(i) = spproj%os_ptcl2D%get_state(i) > 0
        enddo
        nactive = count(active)
        ncls2D  = spproj%os_cls2D%get_noris()
        call make_class_samples(params, spproj, 'class', params%nthr, .false., clssmp)
        ! the training set takes each class's best 2D scores
        call spproj%os_ptcl2D%rank_class_samples(clssmp, 'corr')
        navail = sum(clssmp(:)%pop)
        ntrain = min(params%nsample, navail)
        if( ntrain < params%nsample )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> VALIDATE_POSE_STUDENT: nsample ', params%nsample, ' exceeds the ', navail, &
                &' particles of the selected classes; all of them train'
        endif
        allocate(states(nptcls), source=0)
        call spproj%os_ptcl2D%sample_balanced(clssmp, ntrain, .true., states, exact=.true.)
        call deallocate_class_samples(clssmp)
        train   = states > 0
        ntrain  = count(train)
        unseen  = active .and. .not. train
        nunseen = count(unseen)
        if( nunseen < 1 ) THROW_HARD('validate_pose_student: no selected particle is left unseen; lower nsample')
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> VALIDATE_POSE_STUDENT: ', ntrain, ' training particles of ', nactive, &
            &' selected, from ', ncls2D, ' classes'
        rt(1) = toc(t0)
        if( l_frozen )then
            ! ---- the frozen teacher: a solve3D directory's result, trained on this sample
            t0 = tic()
            call spproj%kill
            fname = teacher_dir%to_char()//'/'//STUDENT_TRAIN_PROJ
            if( .not. file_exists(fname) ) THROW_HARD('dir_refine holds no training_set.simple: not a teacher''s solve3D directory')
            call tproj%read(fname)
            if( tproj%os_ptcl3D%get_noris() /= nptcls ) THROW_HARD('the frozen teacher holds another particle count')
            do i = 1,nptcls
                if( (tproj%os_ptcl3D%get_state(i) > 0) .neqv. train(i) ) THROW_HARD('the frozen teacher was trained on another sample')
            enddo
            call tproj%kill
            result_proj = simple_abspath(fname)
            write(logfhandle,'(A,A)') '>>> VALIDATE_POSE_STUDENT: THE FROZEN TEACHER ', result_proj%to_char()
            rt(8) = rt(8) + toc(t0)
        else
            ! ---- the training project: the sample alone selected
            t0 = tic()
            do i = 1,nptcls
                if( train(i) ) cycle
                call spproj%os_ptcl2D%set_state(i, 0)
                call spproj%os_ptcl3D%set_state(i, 0)
            enddo
            fname = edir%to_char()//'/'//STUDENT_TRAIN_PROJ
            call spproj%projinfo%set(1, 'projfile', fname%to_char())
            call spproj%projinfo%set(1, 'projname', 'training_set')
            call spproj%write(fname)
            call spproj%kill
            rt(8) = rt(8) + toc(t0)
            ! ---- the teacher: solve3D on every training particle in every iteration (nsample at the sample's size forces
            !      full sampling: no fractional update, no trailing reconstruction)
            write(logfhandle,'(A)') '>>> VALIDATE_POSE_STUDENT: SOLVE3D ON THE TRAINING SET'
            t0 = tic()
            cl = cline_entry
            call cl%set('prg',      'solve3D')
            call cl%set('projfile', fname%to_char())
            call cl%set('nsample',  ntrain)
            call cl%set('mkdir',    'yes')
            call xsolve3D%execute(cl)
            call simple_chdir(edir)
            result_proj = simple_abspath(cl%get_carg('projfile'))
            call cl%kill
            rt(2) = toc(t0)
        endif
        ! ---- the student: trained on the teacher's poses and map
        t0 = tic()
        call res%read(result_proj)
        if( .not. res%isthere_in_osout('vol', 1) ) THROW_HARD('the solve3D result registers no volume')
        call res%get_vol('vol', 1, vol_fname, smpd_vol, box_vol)
        vol_fname = simple_abspath(vol_fname)
        allocate(sh_2D(2,nptcls), euls_ref(3,nptcls), sh_ref(2,nptcls), &
            &euls_sk(3,nptcls), sh_sk(2,nptcls), source=0.)
        do i = 1,nptcls
            sh_2D(:,i) = res%os_ptcl2D%get_2Dshift(i)
        enddo
        call enter_subdir('learn')
        fname = edir%to_char()//'/learn/'//STUDENT_TRAIN_PROJ
        call res%projinfo%set(1, 'projfile', fname%to_char())
        call res%write(fname)
        rt(8) = rt(8) + toc(t0)
        write(logfhandle,'(A)') '>>> VALIDATE_POSE_STUDENT: STUDENT TRAINING'
        t0 = tic()
        call cl%set('prg',      'learn_pose_student')
        call cl%set('projfile', fname%to_char())
        call cl%set('vol1',     vol_fname%to_char())
        call cl%set('lp',       STUDENT_LP_VALIDATE)
        call cl%set('mskdiam',  params%mskdiam)
        call cl%set('pgrp',     trim(params%pgrp))
        call cl%set('nthr',     params%nthr)
        call cl%set('fname',    POSE_STUDENT_MODEL_FNAME)
        call cl%set('shift_sigma', STUDENT_SHIFT_SIGMA)
        call xlearn%execute(cl)
        call cl%kill
        call simple_chdir(edir)
        model_abs = simple_abspath(string('learn/'//POSE_STUDENT_MODEL_FNAME))
        rt(3) = toc(t0)
        ! ---- the reference (the physical posterior's best pose, from the 2D solution's shifts): a frozen teacher's
        !      directory may hold it, else the dense table's pass makes it
        t0 = tic()
        call dense_pass('reference', POSE_REFERENCE_FNAME, has_ref, euls_ref, sh_ref, ref_cost)
        rt(4) = toc(t0)
        ! ---- the prediction and the alignment: the unseen particles at their 2D shifts, the reference's pose stored for
        !      the table's reference columns
        t0 = tic()
        call select_unseen(res)
        do i = 1,nptcls
            if( .not. unseen(i) ) cycle
            call res%os_ptcl3D%set_euler(i, [0.,0.,0.])
            if( has_ref(i) ) call res%os_ptcl3D%set_euler(i, euls_ref(:,i))
        enddo
        call enter_subdir('predict')
        fname = edir%to_char()//'/predict/'//STUDENT_UNSEEN_PROJ
        call res%projinfo%set(1, 'projfile', fname%to_char())
        call res%projinfo%set(1, 'projname', 'unseen_particles')
        call res%write(fname)
        rt(8) = rt(8) + toc(t0)
        write(logfhandle,'(A)') '>>> VALIDATE_POSE_STUDENT: STUDENT PREDICTION AND ALIGNMENT FOR THE UNSEEN PARTICLES'
        t0 = tic()
        call cl%set('prg',      'predict_poses')
        call cl%set('projfile', fname%to_char())
        call cl%set('infile',   model_abs%to_char())
        call cl%set('vol1',     vol_fname%to_char())
        call cl%set('mskdiam',  params%mskdiam)
        call cl%set('nthr',     params%nthr)
        call cl%set('fname',    POSE_STUDENT_MODEL_FNAME)
        call set_parts(cl)
        call xpredict%execute(cl)
        call cl%kill
        call simple_chdir(edir)
        table_abs = simple_abspath(string('predict/'//POSE_STUDENT_PRED_TABLE))
        times_abs = simple_abspath(string('predict/'//POSE_STUDENT_PRED_TIMES))
        rt(5) = toc(t0)
        ! ---- the reads and the report
        t0 = tic()
        call read_times(times_abs, times)
        call fopen(funit, string(VALIDATE_POSE_STUDENT_REPORT), status='replace', action='write', iostat=io_stat)
        call fileiochk('validate_pose_student; open '//VALIDATE_POSE_STUDENT_REPORT, io_stat)
        call fopen(dunit, string(VALIDATE_POSE_STUDENT_DETAILS), status='replace', action='write', iostat=io_stat)
        call fileiochk('validate_pose_student; open '//VALIDATE_POSE_STUDENT_DETAILS, io_stat)
        call write_header(funit)
        call report_reads(funit, dunit, model_abs, table_abs, trim(params%pgrp), params%smpd, nptcls, unseen, has_ref, &
            &euls_ref, sh_ref, sh_2D, ref_cost, times, euls_sk, sh_sk)
        call fclose(dunit)
        rt(6) = toc(t0)
        ! ---- the maps from the aligned poses and from the reference's; against the teacher's map they are
        !      cross-validated (disjoint particles)
        write(logfhandle,'(A)') '>>> VALIDATE_POSE_STUDENT: RECONSTRUCTIONS FROM THE POSES'
        t0 = tic()
        call rec_from_poses(1, unseen, euls_sk, sh_sk)
        call rec_from_poses(2, unseen .and. has_ref, euls_ref, sh_ref)
        call write_recs(funit)
        rt(7) = toc(t0)
        rt_all = toc(t_all)
        call write_timings(funit)
        call fclose(funit)
        if( L_BENCH_GLOB )then
            call open_bench_report(fnr, 'VALIDATE_POSE_STUDENT')
            write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
            write(fnr,'(a,i0)')   'validation selected particles       : ', nactive
            write(fnr,'(a,i0)')   'validation training particles       : ', ntrain
            write(fnr,'(a,i0)')   'validation unseen particles         : ', nunseen
            write(fnr,'(a,i0)')   'validation classes (cls2D)          : ', ncls2D
            write(fnr,'(a,i0)')   'validation threads                  : ', params%nthr
            write(fnr,'(a,f0.2)') 'validation student band (lp, A)     : ', STUDENT_LP_VALIDATE
            write(fnr,'(a,l1)')   'validation frozen teacher           : ', l_frozen
            write(fnr,'(a)') ''
            call write_timings(fnr)
            call fclose(fnr)
        endif
        call cline_entry%kill
        call res%kill
        deallocate(states, active, train, unseen, has_ref, sh_2D, euls_ref, sh_ref, euls_sk, sh_sk)
        call simple_end('**** VALIDATE_POSE_STUDENT NORMAL STOP ****')

    contains

        !> a step's command line given the entry's parts (the teacher's solve3D inherits them with the entry's line)
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

        !> the unseen particles alone selected in proj at the 2D solution's shifts, their shifts not yet refined
        subroutine select_unseen( proj )
            type(sp_project), intent(inout) :: proj
            integer :: j
            do j = 1,nptcls
                if( unseen(j) )then
                    call proj%os_ptcl2D%set_state(j, 1)
                    call proj%os_ptcl3D%set_state(j, 1)
                    call proj%os_ptcl3D%set_euler(j, [0.,0.,0.])
                    call proj%os_ptcl3D%set_shift(j, sh_2D(:,j))
                    call proj%os_ptcl3D%set(j, 'shsig',  -1.)
                    call proj%os_ptcl3D%set(j, 'shincr',  0.)
                else
                    call proj%os_ptcl2D%set_state(j, 0)
                    call proj%os_ptcl3D%set_state(j, 0)
                endif
            enddo
        end subroutine select_unseen

        !> the poses of a pass of the dense table over the unseen particles (the physical posterior's best, from the 2D
        !! solution's shifts): the frozen teacher's file fname when it holds one, else the pass in <tag>/, its poses
        !! copied to fname in the execution directory
        subroutine dense_pass( tag, fname_poses, has, euls, shs, cost )
            character(len=*), intent(in)    :: tag, fname_poses
            logical,          intent(inout) :: has(nptcls)
            real,             intent(inout) :: euls(3,nptcls), shs(2,nptcls)
            real(dp),         intent(out)   :: cost(2)
            type(sp_project) :: proj
            type(string)     :: poses, made
            integer          :: status
            poses = ''
            if( l_frozen )then
                poses = teacher_dir%to_char()//'/'//fname_poses
                if( .not. file_exists(poses) ) poses = ''
            endif
            if( poses%strlen_trim() > 0 )then
                write(logfhandle,'(A,A,A,A)') '>>> VALIDATE_POSE_STUDENT: ', tag, ' FROM THE FROZEN TEACHER ', poses%to_char()
            else
                write(logfhandle,'(A,A)') '>>> VALIDATE_POSE_STUDENT: THE DENSE TABLE AT THE PHYSICAL TEMPERATURE, ', tag
                call proj%read(result_proj)
                call select_unseen(proj)
                call enter_subdir(tag)
                fname = edir%to_char()//'/'//tag//'/'//STUDENT_UNSEEN_PROJ
                call proj%projinfo%set(1, 'projfile', fname%to_char())
                call proj%projinfo%set(1, 'projname', 'unseen_particles')
                call proj%write(fname)
                call proj%kill
                call cl%set('prg',      'learn_pose_student')
                call cl%set('projfile', fname%to_char())
                call cl%set('vol1',     vol_fname%to_char())
                call cl%set('lp',       STUDENT_LP_VALIDATE)
                call cl%set('mskdiam',  params%mskdiam)
                call cl%set('pgrp',     trim(params%pgrp))
                call cl%set('nthr',     params%nthr)
                call pose_posterior_mode(cl)
                call cl%kill
                call simple_chdir(edir)
                made = fname_base(string(POSE_REFERENCE_MODEL_FNAME))
                made = edir%to_char()//'/'//tag//'/'//made%to_char()//'_poses.txt'
                poses = edir%to_char()//'/'//fname_poses
                call simple_copy_file(made, poses, status)
                if( status /= 0 ) THROW_HARD('could not copy a dense pass''s poses')
            endif
            call read_pose_file(poses, nptcls, unseen, has, euls, shs, cost)
            write(logfhandle,'(A,A,A,I0,A,I0,A)') '>>> VALIDATE_POSE_STUDENT: ', tag, ' holds ', count(has), ' of ', &
                &nunseen, ' unseen particles'
        end subroutine dense_pass

        !> the validation's provenance: the teacher, the sample, the model's calibration
        subroutine write_header( u )
            integer, intent(in) :: u
            integer :: iu, units(2)
            units = [u, logfhandle]
            do iu = 1,2
                write(units(iu),'(a)') '*** POSE STUDENT VALIDATION ***'
                write(units(iu),'(a,a)') 'teacher: ', result_proj%to_char()
                write(units(iu),'(a,l1)') 'teacher frozen (dir_refine): ', l_frozen
                write(units(iu),'(a,i0,a,i0,a,i0,a)') 'particles: ', nactive, ' selected, ', ntrain, ' training, ', &
                    &nunseen, ' unseen (predicted at their 2D shifts)'
            enddo
        end subroutine write_header

        !> reconstruct3D of the masked rows at the given poses and shifts in rec_<name>/, and its map's FSC against the
        !! teacher's map and between its halves
        subroutine rec_from_poses( ir, mask, euls, shifts )
            integer, intent(in) :: ir
            logical, intent(in) :: mask(nptcls)
            real,    intent(in) :: euls(3,nptcls), shifts(2,nptcls)
            type(sp_project) :: proj
            type(string)     :: dir, recname, even, odd
            integer :: j
            dir = 'rec_'//trim(REC_NAMES(ir))
            call enter_subdir(dir%to_char())
            call proj%read(string(edir%to_char()//'/predict/'//STUDENT_UNSEEN_PROJ))
            do j = 1,nptcls
                if( mask(j) )then
                    call proj%os_ptcl2D%set_state(j, 1)
                    call proj%os_ptcl3D%set_state(j, 1)
                    call proj%os_ptcl3D%set_euler(j, euls(:,j))
                    call proj%os_ptcl3D%set_shift(j, shifts(:,j))
                else
                    call proj%os_ptcl2D%set_state(j, 0)
                    call proj%os_ptcl3D%set_state(j, 0)
                endif
            enddo
            recname = edir%to_char()//'/'//dir%to_char()//'/rec.simple'
            call proj%projinfo%set(1, 'projfile', recname%to_char())
            call proj%projinfo%set(1, 'projname', 'rec')
            call proj%write(recname)
            call proj%kill
            call cl%set('prg',      'reconstruct3D')
            call cl%set('projfile', recname%to_char())
            call cl%set('mkdir',    'no')
            call cl%set('mskdiam',  params%mskdiam)
            call cl%set('pgrp',     trim(params%pgrp))
            call cl%set('nthr',     params%nthr)
            call set_parts(cl)
            call xrec3D%execute(cl)
            call cl%kill
            call simple_chdir(edir)
            rec_n(ir) = count(mask)
            recname = dir%to_char()//'/recvol_state01.mrc'
            even    = dir%to_char()//'/recvol_state01_even_unfil.mrc'
            odd     = dir%to_char()//'/recvol_state01_odd_unfil.mrc'
            if( .not. file_exists(even) )then
                even = dir%to_char()//'/recvol_state01_even.mrc'
                odd  = dir%to_char()//'/recvol_state01_odd.mrc'
            endif
            call map_fsc(recname, vol_fname, params%box, params%smpd, 0.5 * params%mskdiam / params%smpd, &
                &rec_res(1,ir), rec_res(2,ir))
            call map_fsc(even, odd, params%box, params%smpd, 0.5 * params%mskdiam / params%smpd, rec_res(3,ir), rec_res(4,ir))
        end subroutine rec_from_poses

        subroutine write_recs( u )
            integer, intent(in) :: u
            integer :: j, iu, units(2)
            units = [u, logfhandle]
            do iu = 1,2
                write(units(iu),'(a)') '*** MAPS FROM THE POSES (no further search; FSC resolution in A) ***'
                write(units(iu),'(a)') '# map                              particles | vs teacher map 0.5 0.143 | halves 0.5 0.143'
                do j = 1,NRECS
                    write(units(iu),'(a32,1x,i9,a,2(1x,f6.2),a,2(1x,f6.2))') REC_NAMES(j), rec_n(j), '   |       ', &
                        &rec_res(1,j), rec_res(2,j), '  |  ', rec_res(3,j), rec_res(4,j)
                enddo
                do j = 1,NRECS
                    write(units(iu),'(a,a,a,a)') '  ', trim(REC_NAMES(j)), ': ', trim(REC_HELP(j))
                enddo
                write(units(iu),'(a)') '  an unseen map against the teacher''s is cross-validated: the two share no particle'
            enddo
        end subroutine write_recs

        !> the phases' wall times, their shares, the per-particle costs and the share accounted for, to unit u and the log
        subroutine write_timings( u )
            integer, intent(in) :: u
            character(len=200) :: line
            integer :: j, iu, units(2)
            units = [u, logfhandle]
            do iu = 1,2
                write(units(iu),'(a)') '*** TIMINGS (s) ***'
                do j = 1,NPHASES
                    write(line,'(a48,a,1x,f10.2,1x,a,f6.1,a)') PHASE_NAMES(j), ':', rt(j), '(', &
                        &100. * rt(j) / max(rt_all, real(1.e-9,timer_int_kind)), ' %)'
                    write(units(iu),'(a)') trim(line)
                enddo
                write(units(iu),'(a48,a,1x,f10.2)') 'total time', ':', rt_all
                write(units(iu),'(a48,a,1x,f10.1)') '% accounted for', ':', 100. * sum(rt) / max(rt_all, real(1.e-9,timer_int_kind))
                write(units(iu),'(a48,a,1x,f10.3)') 'solve3D per training particle (ms)', ':', 1000. * rt(2) / real(max(1,ntrain))
                write(units(iu),'(a48,a,1x,f10.3)') 'training per training particle (ms)', ':', 1000. * rt(3) / real(max(1,ntrain))
                write(units(iu),'(a)') 'phase breakdowns: BENCH/ of the solve3D directory, learn/BENCH/POSE_STUDENT_LEARN_BENCH.txt, '//&
                    &'predict/BENCH/POSE_STUDENT_PREDICT_BENCH.txt'
            enddo
        end subroutine write_timings

    end subroutine exec_validate_pose_student

    !> the poses of a dense pass (iptcl e1 e2 e3 x y ...; its costs per particle after 'fill:' in the header); a row of
    !! a particle that is not unseen stops the validation
    subroutine read_pose_file( fname, nptcls, unseen, has, euls, shs, cost )
        use simple_defs,   only: dp
        use simple_string, only: string
        use simple_fileio, only: fopen, fclose, fileiochk
        class(string), intent(in)    :: fname
        integer,       intent(in)    :: nptcls
        logical,       intent(in)    :: unseen(nptcls)
        logical,       intent(inout) :: has(nptcls)
        real,          intent(inout) :: euls(3,nptcls), shs(2,nptcls)
        real(dp),      intent(out)   :: cost(2)
        character(len=1024) :: line
        real    :: v(6)
        integer :: u, ios, iptcl, j
        cost = -1.d0
        call fopen(u, fname, status='old', action='read', iostat=ios)
        call fileiochk('validate_pose_student; open '//fname%to_char(), ios)
        do
            read(u,'(A)', iostat=ios) line
            if( ios /= 0 ) exit
            if( line(1:1) == '#' )then
                j = index(line, 'fill:')
                if( j > 0 ) read(line(j+5:),*, iostat=ios) cost
                cycle
            endif
            read(line,*, iostat=ios) v
            if( ios /= 0 ) cycle
            iptcl = nint(v(1))
            if( iptcl < 1 .or. iptcl > nptcls ) THROW_HARD('a pose of a dense pass names no particle of the project')
            if( .not. unseen(iptcl) ) THROW_HARD('a pose of a dense pass is not an unseen particle''s: another sample')
            euls(:,iptcl) = v(2:4)
            shs(:,iptcl)  = v(5:6)
            has(iptcl)    = .true.
        enddo
        call fclose(u)
    end subroutine read_pose_file

    !> the prediction's and the alignment's wall times (s) of predict_poses' times table
    subroutine read_times( fname, times )
        use simple_string, only: string
        use simple_fileio, only: fopen, fclose, fileiochk, file_exists
        class(string), intent(in)  :: fname
        real,          intent(out) :: times(2)
        character(len=256) :: line
        integer :: u, ios, n
        times = -1.
        n     = 0
        ! a distributed prediction's parts time themselves: no table of the whole run
        if( .not. file_exists(fname) ) return
        call fopen(u, fname, status='old', action='read', iostat=ios)
        call fileiochk('validate_pose_student; open '//fname%to_char(), ios)
        do
            read(u,'(A)', iostat=ios) line
            if( ios /= 0 ) exit
            if( line(1:1) == '#' ) cycle
            read(line,*, iostat=ios) n, times
            if( ios /= 0 ) n = 0
            exit
        enddo
        call fclose(u)
        if( n > 0 )then
            times = times / real(n)
        else
            times = -1.
        endif
    end subroutine read_times

    !> The reads on the unseen particles against the reference, per arm (RO the readout at the stored shift sh_in, SK
    !! the alignment over the cone), with a READS line each, the details, the costs per particle; SK's committed poses
    !! for the map (doc/algorithms/pose_student.md)
    subroutine report_reads( u, dunit, model, table, pgrp, smpd, nptcls, unseen, has_ref, euls_ref, sh_ref, sh_in, &
        &ref_cost, times, euls_sk, sh_sk )
        use simple_defs,          only: logfhandle, dp
        use simple_string,        only: string
        use simple_fileio,        only: fopen, fclose, fileiochk
        use simple_sym,           only: sym
        use simple_srch_sort_loc, only: hpsort
        use simple_ori,           only: ori
        use simple_pose_student,  only: pose_student
        integer,          intent(in)  :: u, dunit, nptcls
        class(string),    intent(in)  :: model, table
        character(len=*), intent(in)  :: pgrp
        real,             intent(in)  :: smpd, times(2)
        logical,          intent(in)  :: unseen(nptcls), has_ref(nptcls)
        real,             intent(in)  :: euls_ref(3,nptcls), sh_ref(2,nptcls), sh_in(2,nptcls)
        real(dp),         intent(in)  :: ref_cost(2)
        real,             intent(out) :: euls_sk(3,nptcls), sh_sk(2,nptcls)
        integer,          parameter   :: NARMS = 2
        character(len=2), parameter   :: ARMS(NARMS) = ['RO', 'SK']
        type(pose_student) :: student
        type(sym) :: symop
        type(ori) :: o_ref, o_st, o_sym
        integer, allocatable :: mode(:), k99(:), route(:), in_cone(:), rank_ref(:)
        real,    allocatable :: conf(:), dist_mode(:), euls_ro(:,:), euls(:,:), sh(:,:), d_dir(:), d_inpl(:), d_sh(:)
        logical, allocatable :: m(:), has(:)
        real    :: vals(22), q(16), conf_threshold, conf_fraction, heldout_seeded, heldout_recovered
        integer :: ia, i, funit, io_stat, iptcl, n, iu, units(2)
        character(len=1024) :: line
        call student%read(model)
        call student%get_threshold(conf_threshold, conf_fraction, heldout_seeded, heldout_recovered)
        call student%kill
        call symop%new(pgrp)
        call o_ref%new(is_ptcl=.true.)
        call o_st%new(is_ptcl=.true.)
        allocate(mode(nptcls), k99(nptcls), route(nptcls), in_cone(nptcls), rank_ref(nptcls), source=0)
        allocate(conf(nptcls), dist_mode(nptcls), d_dir(nptcls), d_inpl(nptcls), d_sh(nptcls), source=-1.)
        allocate(euls_ro(3,nptcls), euls(3,nptcls), sh(2,nptcls), source=0.)
        allocate(m(nptcls), has(nptcls), source=.false.)
        euls_sk = 0.
        sh_sk   = 0.
        ! the prediction table (predict_poses' columns): the committed pose is the cone's alignment, then the readout
        call fopen(funit, table, status='old', action='read', iostat=io_stat)
        call fileiochk('validate_pose_student; open '//table%to_char(), io_stat)
        do
            read(funit,'(A)', iostat=io_stat) line
            if( io_stat /= 0 ) exit
            if( line(1:1) == '#' ) cycle
            read(line,*, iostat=io_stat) vals
            if( io_stat /= 0 ) cycle
            iptcl = nint(vals(1))
            if( iptcl < 1 .or. iptcl > nptcls ) cycle
            mode(iptcl)      = nint(vals(3))
            k99(iptcl)       = nint(vals(5))
            conf(iptcl)      = vals(7)
            dist_mode(iptcl) = vals(10)
            rank_ref(iptcl)  = nint(vals(12))
            in_cone(iptcl)   = nint(vals(13))
            route(iptcl)     = nint(vals(14))
            euls_sk(:,iptcl) = vals(15:17)
            sh_sk(:,iptcl)   = vals(18:19)
            euls_ro(:,iptcl) = vals(20:22)
            has(iptcl)       = mode(iptcl) > 0
        enddo
        call fclose(funit)
        write(dunit,'(A)') '# iptcl arm route in_cone rank_ref dist_mode_deg d_dir_deg d_inpl_deg d_shift_A conf k99'
        units = [u, logfhandle]
        do iu = 1,2
            write(units(iu),'(a,f0.4,a,f0.3,a)') 'model: confidence threshold ', conf_threshold, ' (teacher in the cone for ', &
                &conf_fraction, ' of the held-out rows)'
            write(units(iu),'(a,f0.4,a,f0.4)') 'model: held-out rows seeded ', heldout_seeded, ', recovered ', heldout_recovered
            write(units(iu),'(a)') '*** THE UNSEEN PARTICLES AGAINST THE REFERENCE (the physical posterior''s best pose) ***'
            write(units(iu),'(a,i0,a,i0,a,i0)') 'unseen particles ', count(unseen), ', with a reference ', &
                &count(unseen .and. has_ref), ', predicted ', count(unseen .and. has)
            call quant(pack(sqrt(sh_ref(1,:)**2 + sh_ref(2,:)**2) * smpd, unseen .and. has_ref), q(1:3))
            write(units(iu),'(a,3(1x,f0.2))') 'the reference''s shift from the picked position (A): q50 q90 q99', q(1:3)
            write(units(iu),'(a)') 'arms: RO the readout''s pose (no search), SK aligned over the cone''s directions and '//&
                &'their mirrors (committed)'
        enddo
        do ia = 1,NARMS
            d_dir = -1.
            d_inpl = -1.
            d_sh = -1.
            select case(ARMS(ia))
                case('RO')
                    euls = euls_ro
                    sh   = sh_in
                    m    = unseen .and. has_ref .and. has
                case DEFAULT
                    euls = euls_sk
                    sh   = sh_sk
                    m    = unseen .and. has_ref .and. has
            end select
            do i = 1,nptcls
                if( .not. m(i) ) cycle
                call o_ref%set_euler(euls_ref(:,i))
                call o_st%set_euler(euls(:,i))
                call symop%sym_dists(o_ref, o_st, o_sym, d_dir(i), d_inpl(i))
                d_sh(i) = sqrt(sum((sh(:,i) - sh_ref(:,i))**2)) * smpd
                write(dunit,'(I0,1X,A,1X,I0,1X,I0,1X,I0,4(1X,F8.3),1X,F8.5,1X,I0)') i, ARMS(ia), route(i), in_cone(i), &
                    &rank_ref(i), dist_mode(i), d_dir(i), d_inpl(i), d_sh(i), conf(i), k99(i)
            enddo
            n = count(m)
            if( n < 1 ) cycle
            q = -1.
            q(1)  = real(count(m .and. route == 1)) / real(n)
            q(2)  = real(count(m .and. in_cone == 1)) / real(n)
            if( count(m .and. route == 1) > 0 ) q(3) = real(count(m .and. route == 1 .and. in_cone == 1)) / &
                &real(count(m .and. route == 1))
            q(4)  = real(count(m .and. dist_mode <= POSE_AGREE_DEG)) / real(n)
            q(5)  = real(count(m .and. dist_mode <= POSE_TWO_STEPS_DEG)) / real(n)
            call quant(pack(dist_mode, m), q(8:9))
            q(16) = real(count(m .and. k99 > KMAX_STUDENT)) / real(n)
            q(6) = real(count(m .and. d_dir <= POSE_AGREE_DEG)) / real(n)
            q(7) = real(count(m .and. d_dir <= POSE_TWO_STEPS_DEG)) / real(n)
            call quant(pack(d_dir, m), q(10:11))
            call quant(pack(d_inpl, m .and. d_dir <= POSE_AGREE_DEG), q(12:13))
            call quant(pack(d_sh, m), q(14:15))
            do iu = 1,2
                write(units(iu),'(a,a,a,i0)') 'arm ', ARMS(ia), ': n ', n
                write(units(iu),'(a,f0.4,a,f0.4,a,f0.4,a,f0.4)') '  seeded ', q(1), ', cone holds the reference''s direction ', &
                    &q(2), ' (seeded: ', q(3), '), cone capped ', q(16)
                write(units(iu),'(a,f0.4,a,f0.4,a,f0.2,a,f0.2,a)') '  mode within 6 deg ', q(4), ', 12 deg ', q(5), &
                    &' (median ', q(8), ', q90 ', q(9), ' deg; folded pair)'
                write(units(iu),'(a,f0.4,a,f0.4,a,f0.2,a,f0.2,a)') '  committed direction within 6 deg ', q(6), &
                    &', 12 deg ', q(7), ' (median ', q(10), ', q90 ', q(11), ' deg; mirror not folded)'
                write(units(iu),'(a,f0.2,a,f0.2,a,f0.2,a,f0.2,a)') '  in-plane error (direction within 6 deg) median ', q(12), &
                    &' q90 ', q(13), ' deg; shift to the reference''s median ', q(14), ' q90 ', q(15), ' A'
                write(units(iu),'(a,a,1x,i0,16(1x,f0.5))') 'READS ', ARMS(ia), n, q
            enddo
        enddo
        do iu = 1,2
            write(units(iu),'(a)') '*** COSTS PER PARTICLE (ms, wall) ***'
            write(units(iu),'(a,f0.4)') 'prediction (features and forward pass)          : ', 1000. * times(1)
            write(units(iu),'(a,f0.4)') 'alignment over the cone (SK)                    : ', 1000. * times(2)
            write(units(iu),'(a,f0.4)') 'the reference pass (the full search), fill      : ', 1000. * real(ref_cost(1))
            write(units(iu),'(a,f0.4)') 'the reference pass, setup, particles and fill   : ', 1000. * real(ref_cost(2))
            write(units(iu),'(a,3(1x,f0.5))') 'COSTS predict sk reference ', 1000. * times(1), 1000. * times(2), &
                &1000. * real(ref_cost(1))
        enddo
        call symop%kill
        call o_ref%kill
        call o_st%kill
        call o_sym%kill

    contains

        !> the median and the 90th percentile of x (q50, q90, q99 when three are asked); -1 when x is empty
        subroutine quant( x, qq )
            real, intent(in)  :: x(:)
            real, intent(out) :: qq(:)
            real, allocatable :: y(:)
            integer :: nn
            qq = -1.
            nn = size(x)
            if( nn < 1 ) return
            y = x
            call hpsort(y)
            qq(1) = y(max(1,nint(0.5*nn)))
            qq(2) = y(max(1,nint(0.9*nn)))
            if( size(qq) > 2 ) qq(3) = y(max(1,nint(0.99*nn)))
        end subroutine quant

    end subroutine report_reads

    !> the FSC of two maps on disk, both soft-masked at mskrad (pixels): the resolutions at 0.5 and 0.143 (A)
    subroutine map_fsc( f1, f2, box, smpd, mskrad, r05, r0143 )
        use simple_image,         only: image
        use simple_string,        only: string
        use simple_estimate_ssnr, only: get_resolution
        class(string), intent(in)  :: f1, f2
        integer,       intent(in)  :: box
        real,          intent(in)  :: smpd, mskrad
        real,          intent(out) :: r05, r0143
        type(image) :: a, b
        real, allocatable :: corrs(:), res(:)
        integer :: k, n
        call a%new([box,box,box], smpd)
        call b%new([box,box,box], smpd)
        call a%read(f1)
        call b%read(f2)
        call a%mask3D_soft(mskrad, backgr=0.)
        call b%mask3D_soft(mskrad, backgr=0.)
        call a%fft
        call b%fft
        n = a%get_filtsz()
        allocate(corrs(n), res(n))
        do k = 1,n
            res(k) = real(box) * smpd / real(k)
        enddo
        call a%fsc(b, corrs)
        call get_resolution(corrs, res, r05, r0143)
        call a%kill
        call b%kill
    end subroutine map_fsc

end module simple_commanders_pose_student
