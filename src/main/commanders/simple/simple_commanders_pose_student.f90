!@descr: the pose student's commanders: learn_pose_student (the teacher's rows and the training), predict_poses (the mode, the cone and the confidence of unseen particles) and the validate_pose_student controller
module simple_commanders_pose_student
use simple_commander_base,        only: commander_base
use simple_cmdline,               only: cmdline
use simple_strategy_pose_student, only: run_learn_pose_student, run_predict_poses
use simple_defs_hhcnn,            only: NSPACE_STUDENT, KMAX_STUDENT, STUDENT_TRS, STUDENT_SHIFT_SIGMA, STUDENT_LP_VALIDATE, &
    &POSE_AGREE_DEG, POSE_TWO_STEPS_DEG, STUDENT_TRAIN_PROJ, STUDENT_ALL_PROJ, POSE_STUDENT_MODEL_FNAME, &
    &POSE_STUDENT_PRED_TABLE, VALIDATE_POSE_STUDENT_REPORT
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
        call run_learn_pose_student(cline)
    end subroutine exec_learn_pose_student

    !> The student's prediction for the project's selected particles under the model infile: per particle the mode
    !! direction, the cone holding 99% of the posterior and the confidence, to the text table beside fname
    !! (simple_private_exec). With nparts it is the master of that many parts; a part is given its range as the master
    !! hands it (part, nparts, fromp, top, outfile)
    subroutine exec_predict_poses( self, cline )
        use simple_error, only: simple_exception
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
        call run_predict_poses(cline)
    end subroutine exec_predict_poses

    !> validate_pose_student: a class-balanced greedy sample of nsample particles of the 2D solution, solve3D on all of
    !! it (full sampling), the student trained on its poses and map, its prediction for every selected particle, the
    !! statistics and the phase timings (doc/algorithms/pose_student.md). Every step but the training runs over nparts
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
        integer,          parameter :: NPHASES = 7   !< the validation's timed phases
        integer,          parameter :: NRECS   = 6   !< the reconstructions from the student's poses and the control
        character(len=*), parameter :: PHASE_NAMES(NPHASES) = [character(len=40) :: &
            &'class-balanced training sample', 'solve3D on the training set (teacher)', &
            &'student training (learn_pose_student)', 'student prediction (predict_poses)', &
            &'statistics', 'reconstructions from the poses', 'project copies and bookkeeping']
        character(len=*), parameter :: REC_NAMES(NRECS) = [character(len=32) :: 'unseen_student_2Dshifts', &
            &'unseen_student_registered', 'training_student_2Dshifts', 'training_student_registered', &
            &'training_student_teacher_shifts', 'training_teacher']
        character(len=*), parameter :: REC_HELP(NRECS) = [character(len=80) :: &
            &'unseen particles: the student''s poses, the 2D shifts (no search)', &
            &'unseen particles: the student''s poses, the shifts registered at them', &
            &'training particles: the student''s poses, the 2D shifts (no search)', &
            &'training particles: the student''s poses, the shifts registered at them', &
            &'training particles: the student''s poses, the teacher''s shifts', &
            &'training particles: the teacher''s poses and shifts (the control)']
        type(commander_solve3D)          :: xsolve3D
        type(commander_learn_pose_student) :: xlearn
        type(commander_predict_poses)    :: xpredict
        type(commander_rec3D)            :: xrec3D
        type(parameters)                 :: params
        type(sp_project)                 :: spproj, res
        type(cmdline)                    :: cline_entry, cl
        type(class_sample), allocatable  :: clssmp(:)
        type(string) :: edir, result_proj, vol_fname, model_abs, table_abs, fname
        integer, allocatable :: states(:), cls(:)
        logical, allocatable :: active(:), train(:)
        real,    allocatable :: euls_teacher(:,:), dshift(:), sh_teacher(:,:), euls_student(:,:), sh_2D(:,:), sh_student(:,:)
        real    :: rec_res(4,NRECS)   ! per reconstruction: FSC 0.5 and 0.143 against the teacher's map, then of its halves
        integer :: rec_n(NRECS)
        real(timer_int_kind) :: rt(NPHASES), rt_all
        integer(timer_int_kind) :: t0, t_all
        real    :: smpd_vol
        integer :: nptcls, navail, ntrain, nactive, i, box_vol, funit, fnr, io_stat, ncls2D
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
        cline_entry = cline
        call params%new(cline)
        call simple_getcwd(edir)
        ! ---- the training sample: class-balanced and greedy (best 2D scores first) over the selected classes
        t0 = tic()
        call spproj%read(params%projfile)
        nptcls = spproj%os_ptcl2D%get_noris()
        if( nptcls < 1 ) THROW_HARD('the project has no particles (ptcl2D)')
        if( spproj%os_ptcl3D%get_noris() /= nptcls ) THROW_HARD('ptcl2D and ptcl3D differ in size')
        if( spproj%os_cls2D%get_noris() < 1 ) THROW_HARD('validate_pose_student needs a 2D solution (cls2D)')
        allocate(active(nptcls), train(nptcls), source=.false.)
        allocate(cls(nptcls), source=0)
        do i = 1,nptcls
            active(i) = spproj%os_ptcl2D%get_state(i) > 0
            cls(i)    = spproj%os_ptcl2D%get_class(i)
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
        train = states > 0
        ntrain = count(train)
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> VALIDATE_POSE_STUDENT: ', ntrain, ' training particles of ', nactive, &
            &' selected, from ', ncls2D, ' classes'
        rt(1) = toc(t0)
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
        rt(7) = rt(7) + toc(t0)
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
        ! ---- the student: trained on the teacher's poses and map
        t0 = tic()
        call res%read(result_proj)
        if( .not. res%isthere_in_osout('vol', 1) ) THROW_HARD('the solve3D result registers no volume')
        call res%get_vol('vol', 1, vol_fname, smpd_vol, box_vol)
        vol_fname = simple_abspath(vol_fname)
        ! the teacher's poses, and how far its shifts moved from the 2D solution's (the unseen rows keep the latter)
        allocate(euls_teacher(3,nptcls), sh_teacher(2,nptcls), euls_student(3,nptcls), sh_2D(2,nptcls), &
            &sh_student(2,nptcls), source=0.)
        allocate(dshift(nptcls), source=-1.)
        do i = 1,nptcls
            sh_2D(:,i) = res%os_ptcl2D%get_2Dshift(i)
            if( .not. train(i) ) cycle
            euls_teacher(:,i) = res%os_ptcl3D%get_euler(i)
            sh_teacher(:,i)   = res%os_ptcl3D%get_2Dshift(i)
            dshift(i) = sqrt(sum((res%os_ptcl3D%get_2Dshift(i) - res%os_ptcl2D%get_2Dshift(i))**2))
        enddo
        call enter_subdir('learn')
        fname = edir%to_char()//'/learn/'//STUDENT_TRAIN_PROJ
        call res%projinfo%set(1, 'projfile', fname%to_char())
        call res%write(fname)
        rt(7) = rt(7) + toc(t0)
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
        ! ---- the prediction: every selected particle, the training rows at the teacher's poses, the others at their
        !      2D shifts and no direction
        t0 = tic()
        do i = 1,nptcls
            if( active(i) )then
                call res%os_ptcl2D%set_state(i, 1)
                call res%os_ptcl3D%set_state(i, 1)
            else
                call res%os_ptcl2D%set_state(i, 0)
                call res%os_ptcl3D%set_state(i, 0)
            endif
            if( active(i) .and. .not. train(i) )then
                call res%os_ptcl3D%set_euler(i, [0.,0.,0.])
                call res%os_ptcl3D%set_shift(i, res%os_ptcl2D%get_2Dshift(i))
            endif
        enddo
        call enter_subdir('predict')
        fname = edir%to_char()//'/predict/'//STUDENT_ALL_PROJ
        call res%projinfo%set(1, 'projfile', fname%to_char())
        call res%projinfo%set(1, 'projname', 'all_particles')
        call res%write(fname)
        call res%kill
        rt(7) = rt(7) + toc(t0)
        write(logfhandle,'(A)') '>>> VALIDATE_POSE_STUDENT: STUDENT PREDICTION FOR EVERY SELECTED PARTICLE'
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
        rt(4) = toc(t0)
        ! ---- the statistics and the report
        t0 = tic()
        call fopen(funit, string(VALIDATE_POSE_STUDENT_REPORT), status='replace', action='write', iostat=io_stat)
        call fileiochk('validate_pose_student; open '//VALIDATE_POSE_STUDENT_REPORT, io_stat)
        call report_statistics(funit, model_abs, table_abs, trim(params%pgrp), nptcls, active, train, cls, euls_teacher, &
            &dshift, 3. * STUDENT_SHIFT_SIGMA / params%smpd, sh_teacher, sh_2D, euls_student, sh_student)
        rt(5) = toc(t0)
        ! ---- the maps from the student's poses straight from its readout, no search and no refinement; the unseen
        !      particles' map against the teacher's is cross-validated (disjoint particles)
        write(logfhandle,'(A)') '>>> VALIDATE_POSE_STUDENT: RECONSTRUCTIONS FROM THE STUDENT''S POSES'
        t0 = tic()
        call rec_from_poses(1, active .and. .not. train, euls_student, sh_2D)
        call rec_from_poses(2, active .and. .not. train, euls_student, sh_student)
        call rec_from_poses(3, train, euls_student, sh_2D)
        call rec_from_poses(4, train, euls_student, sh_student)
        call rec_from_poses(5, train, euls_student, sh_teacher)
        call rec_from_poses(6, train, euls_teacher, sh_teacher)
        call write_recs(funit)
        rt(6) = toc(t0)
        rt_all = toc(t_all)
        call write_timings(funit)
        call fclose(funit)
        if( L_BENCH_GLOB )then
            call open_bench_report(fnr, 'VALIDATE_POSE_STUDENT')
            write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
            write(fnr,'(a,i0)')   'validation selected particles       : ', nactive
            write(fnr,'(a,i0)')   'validation training particles       : ', ntrain
            write(fnr,'(a,i0)')   'validation classes (cls2D)          : ', ncls2D
            write(fnr,'(a,i0)')   'validation threads                  : ', params%nthr
            write(fnr,'(a,f0.2)') 'validation student band (lp, A)     : ', STUDENT_LP_VALIDATE
            write(fnr,'(a)') ''
            call write_timings(fnr)
            call fclose(fnr)
        endif
        call cline_entry%kill
        deallocate(states, cls, active, train, euls_teacher, dshift, sh_teacher, euls_student, sh_2D, sh_student)
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
            call proj%read(string(edir%to_char()//'/predict/'//STUDENT_ALL_PROJ))
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
                write(units(iu),'(a)') '*** MAPS FROM THE POSES (no search, no refinement; FSC resolution in A) ***'
                write(units(iu),'(a)') '# map                              particles | vs teacher map 0.5 0.143 | halves 0.5 0.143'
                do j = 1,NRECS
                    write(units(iu),'(a32,1x,i9,a,2(1x,f6.2),a,2(1x,f6.2))') REC_NAMES(j), rec_n(j), '   |       ', &
                        &rec_res(1,j), rec_res(2,j), '  |  ', rec_res(3,j), rec_res(4,j)
                enddo
                do j = 1,NRECS
                    write(units(iu),'(a,a,a,a)') '  ', trim(REC_NAMES(j)), ': ', trim(REC_HELP(j))
                enddo
                write(units(iu),'(a)') '  the unseen map against the teacher''s is cross-validated: the two share no particle'
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
                    write(line,'(a40,a,1x,f10.2,1x,a,f6.1,a)') PHASE_NAMES(j), ':', rt(j), '(', &
                        &100. * rt(j) / max(rt_all, real(1.e-9,timer_int_kind)), ' %)'
                    write(units(iu),'(a)') trim(line)
                enddo
                write(units(iu),'(a40,a,1x,f10.2)') 'total time', ':', rt_all
                write(units(iu),'(a40,a,1x,f10.1)') '% accounted for', ':', 100. * sum(rt) / max(rt_all, real(1.e-9,timer_int_kind))
                write(units(iu),'(a40,a,1x,f10.3)') 'solve3D per training particle (ms)', ':', 1000. * rt(2) / real(max(1,ntrain))
                write(units(iu),'(a40,a,1x,f10.3)') 'training per training particle (ms)', ':', 1000. * rt(3) / real(max(1,ntrain))
                write(units(iu),'(a40,a,1x,f10.3)') 'prediction per selected particle (ms)', ':', 1000. * rt(4) / real(max(1,nactive))
                if( rt(4) > 0. ) write(units(iu),'(a40,a,1x,f10.1)') 'prediction throughput (particles/s)', ':', &
                    &real(nactive) / rt(4)
                write(units(iu),'(a)') 'phase breakdowns: BENCH/ of the solve3D directory, learn/BENCH/POSE_STUDENT_LEARN_BENCH.txt, '//&
                    &'predict/BENCH/POSE_STUDENT_PREDICT_BENCH.txt'
            enddo
        end subroutine write_timings

    end subroutine exec_validate_pose_student

    !> The validation's statistics from the model and its prediction table, to unit u and the log. Distances fold the
    !! point group and the mirror pair. Unseen rows have no teacher but a 2D class: a class's teacher medoid is the
    !! training pose nearest all the class's others. dshift: the training rows' 3D-2D shift distance (px), r_aug the disk;
    !! euls_student, sh_student: the readout's full pose and the registered shift of every predicted row (the table's
    !! last columns), compared with sh_teacher and the 2D solution's sh_2D
    subroutine report_statistics( u, model, table, pgrp, nptcls, active, train, cls, euls_teacher, dshift, r_aug, &
        &sh_teacher, sh_2D, euls_student, sh_student )
        use simple_defs,          only: logfhandle
        use simple_string,        only: string
        use simple_fileio,        only: fopen, fclose, fileiochk
        use simple_sym,           only: sym
        use simple_ori_utils,     only: euler_mirror
        use simple_srch_sort_loc, only: hpsort
        use simple_pose_student,  only: pose_student
        use simple_ori,           only: ori
        integer,          intent(in) :: u, nptcls, cls(nptcls)
        class(string),    intent(in) :: model, table
        character(len=*), intent(in) :: pgrp
        logical,          intent(in) :: active(nptcls), train(nptcls)
        real,             intent(in) :: euls_teacher(3,nptcls), dshift(nptcls), r_aug
        real,             intent(in)  :: sh_teacher(2,nptcls), sh_2D(2,nptcls)
        real,             intent(out) :: euls_student(3,nptcls), sh_student(2,nptcls)
        type(pose_student) :: student
        type(ori)          :: o_st, o_te, o_sym
        type(sym)          :: symop
        integer, allocatable :: mode(:), k99(:), route(:), in_cone(:), medoid(:), members(:)
        real,    allocatable :: conf(:), radius(:), d_teach(:), d_pred_med(:), d_teach_med(:), dmat(:,:)
        real,    allocatable :: d_dir(:), d_inpl(:), e_reg(:), e_2D(:)
        real    :: vals(19), conf_threshold, conf_fraction, heldout_seeded, heldout_recovered
        integer :: i, j, c, ncls, nm, io_stat, funit, iptcl, iu, units(2)
        character(len=512) :: line
        call student%read(model)
        call student%get_threshold(conf_threshold, conf_fraction, heldout_seeded, heldout_recovered)
        call symop%new(pgrp)
        allocate(mode(nptcls), k99(nptcls), route(nptcls), in_cone(nptcls), source=0)
        allocate(conf(nptcls), radius(nptcls), d_teach(nptcls), d_pred_med(nptcls), d_teach_med(nptcls), &
            &d_dir(nptcls), d_inpl(nptcls), e_reg(nptcls), e_2D(nptcls), source=-1.)
        euls_student = 0.
        sh_student   = sh_2D
        ! the prediction table: iptcl eo mode pmax k99 ncone conf radius_deg iproj_ref dist_mode dist_cone rank_ref in_cone route
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
            mode(iptcl)    = nint(vals(3))
            k99(iptcl)     = nint(vals(5))
            conf(iptcl)    = vals(7)
            radius(iptcl)  = vals(8)
            in_cone(iptcl) = nint(vals(13))
            route(iptcl)   = nint(vals(14))
            euls_student(:,iptcl) = vals(15:17)
            sh_student(:,iptcl)   = vals(18:19)
        enddo
        call fclose(funit)
        ! the mode against the teacher's pose (training rows); the readout's full pose against it, unfolded: the
        ! direction with the mirror the readout chose, and the in-plane angle
        call o_st%new(is_ptcl=.true.)
        call o_te%new(is_ptcl=.true.)
        do i = 1,nptcls
            if( .not. (train(i) .and. mode(i) > 0) ) cycle
            d_teach(i) = pair_angle(mode(i), euls_teacher(:,i))
            call o_st%set_euler(euls_student(:,i))
            call o_te%set_euler(euls_teacher(:,i))
            call symop%sym_dists(o_te, o_st, o_sym, d_dir(i), d_inpl(i))
            e_reg(i) = sqrt(sum((sh_student(:,i) - sh_teacher(:,i))**2))
            e_2D(i)  = sqrt(sum((sh_2D(:,i)      - sh_teacher(:,i))**2))
        enddo
        call o_st%kill
        call o_te%kill
        call o_sym%kill
        ! the classes' teacher medoids and every row's distance to its class's
        ncls = maxval(cls)
        allocate(medoid(max(1,ncls)), source=0)
        do c = 1,ncls
            members = pack([(i, i=1,nptcls)], train .and. cls == c)
            nm = size(members)
            if( nm < 2 ) cycle
            allocate(dmat(nm,nm), source=0.)
            do i = 1,nm
                do j = i+1,nm
                    dmat(i,j) = euls_angle(euls_teacher(:,members(i)), euls_teacher(:,members(j)))
                    dmat(j,i) = dmat(i,j)
                enddo
            enddo
            medoid(c) = members(minloc(sum(dmat, dim=2), dim=1))
            deallocate(dmat)
        enddo
        do i = 1,nptcls
            if( .not. active(i) ) cycle
            if( cls(i) < 1 .or. cls(i) > ncls ) cycle
            if( medoid(cls(i)) == 0 ) cycle
            if( mode(i) > 0 ) d_pred_med(i) = pair_angle(mode(i), euls_teacher(:,medoid(cls(i))))
            if( train(i) .and. i /= medoid(cls(i)) ) d_teach_med(i) = &
                &euls_angle(euls_teacher(:,i), euls_teacher(:,medoid(cls(i))))
        enddo
        units = [u, logfhandle]
        do iu = 1,2
            write(units(iu),'(a)') '*** POSE STUDENT VALIDATION ***'
            write(units(iu),'(a,f0.4,a,f0.3,a)') 'model: confidence threshold ', conf_threshold, ' (teacher in the cone for ', &
                &conf_fraction, ' of the held-out rows)'
            write(units(iu),'(a,f0.4,a,f0.4)') 'model: held-out rows seeded ', heldout_seeded, ', recovered ', heldout_recovered
            call shift_stats(units(iu))
            call group_stats(units(iu), 'all selected particles', active, .false.)
            call group_stats(units(iu), 'training particles (in sample)', active .and. train, .true.)
            call group_stats(units(iu), 'unseen particles', active .and. .not. train, .false.)
            write(units(iu),'(a)') 'student mode against the teacher (training particles, folded pair and point group):'
            call dist_stats(units(iu), '  all', d_teach, active .and. train)
            call dist_stats(units(iu), '  seeded', d_teach, active .and. train .and. route == 1)
            write(units(iu),'(a)') 'the readout''s full pose against the teacher''s (training particles, no search; point group folded, mirror not):'
            call dist_stats(units(iu), '  direction', d_dir, active .and. train)
            call dist_stats(units(iu), '  in-plane, direction within 12 deg', d_inpl, active .and. train .and. d_dir <= POSE_TWO_STEPS_DEG)
            call dist_stats(units(iu), '  direction, seeded', d_dir, active .and. train .and. route == 1)
            call dist_stats(units(iu), '  in-plane, seeded and direction within 12 deg', d_inpl, &
                &active .and. train .and. route == 1 .and. d_dir <= POSE_TWO_STEPS_DEG)
            write(units(iu),'(a)') 'shift against the teacher''s (training particles):'
            call shift_err_stats(units(iu), '  2D shift', e_2D, active .and. train)
            call shift_err_stats(units(iu), '  registered at the read-out pose', e_reg, active .and. train)
            call shift_err_stats(units(iu), '  2D shift, direction within 6 deg', e_2D, active .and. train .and. d_dir <= POSE_AGREE_DEG)
            call shift_err_stats(units(iu), '  registered, direction within 6 deg', e_reg, active .and. train .and. d_dir <= POSE_AGREE_DEG)
            write(units(iu),'(a)') 'class consistency: distance to the 2D class''s teacher medoid (classes with two or more training particles):'
            call dist_stats(units(iu), '  teacher, training particles (baseline)', d_teach_med, active .and. train)
            call dist_stats(units(iu), '  student, training particles', d_pred_med, active .and. train)
            call dist_stats(units(iu), '  student, unseen particles', d_pred_med, active .and. .not. train)
            call dist_stats(units(iu), '  student, unseen and seeded', d_pred_med, active .and. .not. train .and. route == 1)
        enddo
        call student%kill
        call symop%kill

    contains

        !> the angle between a grid direction's mirror pair and an orientation (the smaller of the two)
        real function pair_angle( iproj, e )
            integer, intent(in) :: iproj
            real,    intent(in) :: e(3)
            pair_angle = min(student%angle_between(iproj, e), student%angle_between(student%mirror_of(iproj), e))
        end function pair_angle

        !> the angle between two orientations' directions, folded by the point group and the mirror
        real function euls_angle( e1, e2 )
            real, intent(in) :: e1(3), e2(3)
            real :: em(3), es(3), d1, d2, inpl
            call symop%sym_dists(e1, e2, es, d1, inpl)
            call euler_mirror(e1, em)
            call symop%sym_dists(em, e2, es, d2, inpl)
            euls_angle = min(d1, d2)
        end function euls_angle

        !> l_teacher: the rows carry the teacher's poses, so the cone's coverage of the stored direction means something
        subroutine group_stats( uu, label, mask, l_teacher )
            integer,          intent(in) :: uu
            character(len=*), intent(in) :: label
            logical,          intent(in) :: mask(nptcls), l_teacher
            real, allocatable :: r(:), cf(:)
            integer :: n
            n = count(mask .and. mode > 0)
            write(uu,'(a,a,i0,a)') label, ': ', n, ' predicted'
            if( n < 1 ) return
            r  = pack(radius, mask .and. mode > 0)
            cf = pack(conf,   mask .and. mode > 0)
            call hpsort(r)
            call hpsort(cf)
            write(uu,'(a,f0.4,a,f0.4)') '  seeded (neighbourhood search from the cone) ', &
                &real(count(mask .and. route == 1)) / real(n), ', cone capped at KMAX_STUDENT ', &
                &real(count(mask .and. mode > 0 .and. k99 > KMAX_STUDENT)) / real(n)
            if( l_teacher ) write(uu,'(a,f0.4)') '  teacher''s direction in the cone ', real(count(mask .and. in_cone == 1)) / real(n)
            write(uu,'(a,3(1x,f0.4))') '  confidence quartiles  ', cf(max(1,nint(0.25*n))), cf(max(1,nint(0.5*n))), &
                &cf(max(1,nint(0.75*n)))
            write(uu,'(a,3(1x,f0.2))') '  cone radius (deg) q50 q90 max', r(max(1,nint(0.5*n))), r(max(1,nint(0.9*n))), r(n)
        end subroutine group_stats

        !> the training rows' shift change from the 2D solution to the teacher's, against the training's shift disk
        subroutine shift_stats( uu )
            integer, intent(in) :: uu
            real, allocatable :: x(:)
            integer :: n
            x = pack(dshift, train .and. dshift >= 0.)
            n = size(x)
            if( n < 1 ) return
            call hpsort(x)
            write(uu,'(a,f0.2,a,f0.2,a,f0.2,a,f0.4,a,f0.2,a)') 'teacher''s shift minus the 2D shift (training, px): median ', &
                &x(max(1,nint(0.5*n))), ', q90 ', x(max(1,nint(0.9*n))), ', max ', x(n), '; within the training''s disk ', &
                &real(count(x <= r_aug)) / real(n), ' (', r_aug, ' px); unseen particles are predicted at their 2D shifts'
        end subroutine shift_stats

        !> a shift error's quantiles (px) and the fractions within 1 and 3 px
        subroutine shift_err_stats( uu, label, e, mask )
            integer,          intent(in) :: uu
            character(len=*), intent(in) :: label
            real,             intent(in) :: e(nptcls)
            logical,          intent(in) :: mask(nptcls)
            real, allocatable :: x(:)
            integer :: n
            x = pack(e, mask .and. e >= 0.)
            n = size(x)
            if( n < 1 )then
                write(uu,'(a,a)') label, ': none'
                return
            endif
            call hpsort(x)
            write(uu,'(a,a,i0,a,f0.2,a,f0.2,a,f0.4,a,f0.4)') label, ': n ', n, ', median ', x(max(1,nint(0.5*n))), &
                &' px, q90 ', x(max(1,nint(0.9*n))), ' px, within 1 px ', real(count(x <= 1.)) / real(n), &
                &', within 3 px ', real(count(x <= 3.)) / real(n)
        end subroutine shift_err_stats

        subroutine dist_stats( uu, label, d, mask )
            integer,          intent(in) :: uu
            character(len=*), intent(in) :: label
            real,             intent(in) :: d(nptcls)
            logical,          intent(in) :: mask(nptcls)
            real, allocatable :: x(:)
            integer :: n
            x = pack(d, mask .and. d >= 0.)
            n = size(x)
            if( n < 1 )then
                write(uu,'(a,a)') label, ': none'
                return
            endif
            call hpsort(x)
            write(uu,'(a,a,i0,a,f0.2,a,f0.2,a,f0.4,a,f0.4)') label, ': n ', n, ', median ', x(max(1,nint(0.5*n))), &
                &' deg, q90 ', x(max(1,nint(0.9*n))), ' deg, within 6 deg ', real(count(x <= POSE_AGREE_DEG)) / real(n), &
                &', within 12 deg ', real(count(x <= POSE_TWO_STEPS_DEG)) / real(n)
        end subroutine dist_stats

    end subroutine report_statistics

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
