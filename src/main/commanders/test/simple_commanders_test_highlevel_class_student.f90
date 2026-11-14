!@descr: validate_class_student in miniature on the multi-state fixture's beta-galactosidase particles with planted shifts, a gate of simple_commanders_test_highlevel
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_class_student
implicit none
#include "simple_local_flags.inc"

contains

    !> validate_class_student on 2 000 simulated beta-galactosidase particles with planted shifts, half of them the
    !! teacher's: the programs' contract gated, the validation's reads reported to metrics.tsv (the test policy's
    !! entry class_student)
    module subroutine exec_test_class_student( self, cline )
        use simple_test_state_fixture,       only: state_fixture, FIX_SEARCH
        use simple_test_gate,                only: test_gate
        use simple_commanders_class_student, only: commander_validate_class_student
        use simple_class_student,            only: class_student
        use simple_timer,                    only: tic, toc, timer_int_kind
        use, intrinsic :: ieee_arithmetic,   only: ieee_is_finite
        class(commander_test_class_student), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        integer, parameter :: NPTCLS   = 2000
        integer, parameter :: NTRAIN   = 1000
        integer, parameter :: NCLS     = 10
        real,    parameter :: SNR      = 0.5
        real,    parameter :: SHERR_PX = 6.    !< the planted shifts' range (pixels at 2.2 A): the center_ptcls test's
        integer, parameter :: NREADS   = 17
        ! each read: its file (1 the report, 2 the diagnostics), its row's label and the column of its number
        integer,          parameter :: READ_FILE(NREADS) = [1,1,1,1,1,1,1,1,1,1, 2,2,2,2,2,2,2]
        integer,          parameter :: READ_COL(NREADS)  = [1,1,1,1,2,1,2,1,1,1, 1,2,1,1,1,1,1]
        character(len=*), parameter :: READ_LABELS(NREADS) = [character(len=50) :: &
            &'Particles given a shortlist (%)', 'Shortlists containing the full search''s class (%)', &
            &'Shortlist length (median)', 'Class averages vs the teacher''s, FRC 0.5 (A)', &
            &'Class averages vs the teacher''s, FRC 0.5 (A)', 'Cost per particle (ms)', 'Cost per particle (ms)', &
            &'same, with the particles perfectly centred (%)', 'Student prediction (predict_classes)', &
            &'Reference: best class of the unseen (refine2D)', &
            &'Student''s class = full search''s class (%)', 'Student''s class = full search''s class (%)', &
            &'Every particle''s list has the full-search class (%)', 'In-plane angle error, same class (deg, median)', &
            &'Registered shift error, same class (A, median)', 'Centring error before registration (A, median)', &
            &'Remaining difference after it (A, median)']
        character(len=*), parameter :: READ_KEYS(NREADS) = [character(len=40) :: 'shortlisted_pct', &
            &'shortlist_has_class_pct', 'shortlist_length_median', 'student_cavgs_frc05_A', &
            &'full_search_cavgs_frc05_A', 'student_cost_ms', 'full_search_cost_ms', &
            &'perfect_centring_shortlist_has_class_pct', 'predict_wall_s', 'full_search_wall_s', &
            &'same_class_pct', 'perfect_centring_same_class_pct', 'list_has_class_pct', 'inpl_error_deg', &
            &'registered_shift_error_A', 'centred_shift_error_A', 'frame_residual_A']
        class(parameters), allocatable :: params
        type(state_fixture)   :: fix
        type(test_gate)       :: gate
        type(sp_project)      :: proj
        type(class_student)   :: student
        type(commander_validate_class_student) :: xvalidate
        type(cmdline)         :: cl
        type(string)          :: cwd_saved, root, projfile, edir
        integer(timer_int_kind) :: t0
        real    :: wall, vals(NREADS), sh(2), e3, pconf, top1
        integer :: i, status, npred, ntrain_pred, cls, rte
        logical :: l_records, l_found(NREADS), l_mtable, l_mcands, l_mrecords, l_mclean, l_empty
        allocate(params)
        call simple_getcwd(cwd_saved)
        root = filepath(cwd_saved, 'test_class_student_'//int2str(get_process_id()))
        if( dir_exists(root) ) call simple_rmdir(root)
        call simple_mkdir(root)
        call simple_chdir(root, status)
        if( status /= 0 ) THROW_HARD('TEST_CLASS_STUDENT FAILED: could not enter the test directory')
        call params%new(cline)
        call gate%new(filepath(cwd_saved, 'metrics.tsv'))
        ! the fixture: raw picks with planted shifts
        call fix%new('betagal', FIX_SEARCH, 1, NPTCLS, SNR, params%nthr, kind='identical', sherr=SHERR_PX)
        call fix%build
        projfile = simple_abspath(fix%project_fname())
        ! the program in its execution directory
        call cl%set('prg',      'validate_class_student')
        call cl%set('projfile', projfile%to_char())
        call cl%set('nsample',  NTRAIN)
        call cl%set('ncls',     NCLS)
        call cl%set('mskdiam',  fix%mskdiam())
        call cl%set('nthr',     params%nthr)
        call cl%set('mkdir',    'yes')
        t0 = tic()
        call xvalidate%execute(cl)
        wall = real(toc(t0))
        call cl%kill
        call simple_getcwd(edir)
        ! the model reads back with the teacher's classes
        call student%read(string('learn/class_student.txt'))
        call gate%check('model_reads_back_with_the_teachers_classes', student%get_ncls() == NCLS)
        call student%kill
        ! the prediction's records: every unseen particle predicted, no training particle
        call proj%read(string('predict/unseen_particles.simple'))
        l_records   = proj%os_ptcl2D%get_noris() == NPTCLS
        npred       = 0
        ntrain_pred = 0
        do i = 1,proj%os_ptcl2D%get_noris()
            if( proj%os_ptcl2D%get_state(i) == 0 )then
                if( proj%os_ptcl2D%isthere(i, 'route') ) ntrain_pred = ntrain_pred + 1
                cycle
            endif
            npred = npred + 1
            cls   = proj%os_ptcl2D%get_class(i)
            e3    = proj%os_ptcl2D%e3get(i)
            sh    = proj%os_ptcl2D%get_2Dshift(i)
            pconf = proj%os_ptcl2D%get(i, 'pconf')
            rte   = nint(proj%os_ptcl2D%get(i, 'route'))
            if( cls < 1 .or. cls > NCLS )                    l_records = .false.
            if( .not. (e3 >= 0. .and. e3 < 360.) )           l_records = .false.
            if( .not. all(ieee_is_finite(sh)) )              l_records = .false.
            if( rte /= 1 .and. rte /= 2 )                    l_records = .false.
            if( .not. (pconf >= 0. .and. pconf <= 1.0001) )  l_records = .false.
        enddo
        call proj%kill
        call gate%check('every_unseen_particle_carries_a_prediction', l_records .and. npred == NPTCLS - NTRAIN)
        call gate%check('no_training_particle_predicted', ntrain_pred == 0)
        ! the distributed master (two parts) against shared memory on the predicted particles, both registering the
        ! shifts against the teacher's averages: the table, the candidates and the records the same, byte for byte
        call master_as_shared_memory(l_mtable, l_mcands, l_mrecords, l_mclean, l_empty)
        call gate%check('master_table_as_shared_memory',      l_mtable)
        call gate%check('master_candidates_as_shared_memory', l_mcands)
        call gate%check('master_records_as_shared_memory',    l_mrecords)
        call gate%check('master_leaves_no_part_file',         l_mclean)
        call gate%check('part_of_unselected_rows_writes_its_contributions', l_empty)
        ! the reads of the report
        call read_reads(vals, l_found)
        call gate%check('the_report_holds_every_read', all(l_found))
        top1 = vals(11) / 100.
        call gate%check('the_student_agrees_with_the_full_search_above_chance', top1 > 2. / real(NCLS))
        do i = 1,NREADS
            call gate%report(trim(READ_KEYS(i)), vals(i))
        enddo
        call gate%report('validate_class_student_wall_s', wall)
        call fix%kill
        call simple_chdir(cwd_saved, status)
        if( status /= 0 ) THROW_HARD('TEST_CLASS_STUDENT FAILED: could not restore the original directory')
        if( gate%passed() )then
            call gate%kill
            call simple_rmdir(root)
            write(logfhandle,'(A)') 'PASS: validate_class_student in miniature'
            call simple_end('**** SIMPLE_TEST_CLASS_STUDENT NORMAL STOP ****')
        else
            call gate%kill
            THROW_HARD('TEST_CLASS_STUDENT FAILED')
        endif

      contains

        !> predict/'s particles predicted again in parts_check/shm (shared memory) and parts_check/distr (the master over
        !! two parts), each from a copy of the project and its sigma2 store, with the model and the teacher's averages;
        !! then one part whose range holds only an unselected row, in parts_check/empty
        subroutine master_as_shared_memory( l_table, l_cands, l_recs, l_clean, l_empty_part )
            use simple_fileio,                   only: files_identical
            use simple_defs_hhcnn,               only: CLASS_STUDENT_MODEL_FNAME, CLASS_STUDENT_PRED_TABLE, STUDENT_UNSEEN_PROJ
            use simple_commanders_class_student, only: commander_predict_classes
            use simple_strategy_student_common,  only: fname_base
            logical, intent(out) :: l_table, l_cands, l_recs, l_clean, l_empty_part
            character(len=*), parameter :: RUN_DIRS(3) = [character(len=17) :: 'parts_check/shm', 'parts_check/distr', &
                &'parts_check/empty']
            character(len=*), parameter :: PART_FILES(5) = [character(len=32) :: 'class_student_part1_pred.txt', &
                &'class_student_part2_cands.bin', 'algndoc_1.simple', 'JOB_FINISHED_2', 'distr_simple_script_1']
            type(commander_predict_classes) :: xpredict
            class(sp_project), allocatable  :: input, runs(:)
            type(string) :: model_abs, cavgs, sigma2_src, sigma2_base, cands_fname, here
            integer :: irun, ncls_c, irow, u, ios, nlines
            real    :: smpd_c
            character(len=512) :: line
            logical :: l_sigma2
            allocate(input, runs(2))
            call simple_getcwd(here)
            call input%read(string('predict/'//STUDENT_UNSEEN_PROJ))
            call input%get_cavgs_stk(cavgs, ncls_c, smpd_c)
            cavgs       = simple_abspath(cavgs)
            model_abs   = simple_abspath(string('learn/'//CLASS_STUDENT_MODEL_FNAME))
            cands_fname = fname_base(string(CLASS_STUDENT_MODEL_FNAME))
            cands_fname = cands_fname%to_char()//'_cands.bin'
            call input%get_sigma2_state_path(sigma2_src, l_sigma2)
            if( l_sigma2 )then
                sigma2_base = basename(sigma2_src)
                sigma2_src  = simple_abspath(string('predict/'//sigma2_base%to_char()))
                l_sigma2    = file_exists(sigma2_src)
            endif
            call simple_mkdir('parts_check')
            do irun = 1,3
                call simple_mkdir(trim(RUN_DIRS(irun)))
                call simple_chdir(string(trim(RUN_DIRS(irun))), status)
                if( status /= 0 ) THROW_HARD('TEST_CLASS_STUDENT FAILED: could not enter a parts check directory')
                if( l_sigma2 )then
                    call simple_copy_file(sigma2_src, sigma2_base)
                    call input%set_sigma2_state_path(sigma2_base)
                endif
                call input%update_projinfo(string(STUDENT_UNSEEN_PROJ))
                call input%write(string(STUDENT_UNSEEN_PROJ))
                call cl%set('prg',      'predict_classes')
                call cl%set('projfile', simple_abspath(string(STUDENT_UNSEEN_PROJ)))
                call cl%set('infile',   model_abs%to_char())
                call cl%set('refs',     cavgs%to_char())
                call cl%set('mskdiam',  fix%mskdiam())
                call cl%set('nthr',     params%nthr)
                call cl%set('fname',    CLASS_STUDENT_MODEL_FNAME)
                if( irun == 2 ) call cl%set('nparts', 2)
                if( irun == 3 )then
                    irow = 1
                    do while( input%os_ptcl2D%get_state(irow) > 0 )
                        irow = irow + 1
                    enddo
                    call cl%set('part',    1)
                    call cl%set('nparts',  1)
                    call cl%set('fromp',   irow)
                    call cl%set('top',     irow)
                    call cl%set('outfile', 'part1.simple')
                endif
                call xpredict%execute(cl)
                call cl%kill
                if( irun < 3 ) call runs(irun)%read(string(STUDENT_UNSEEN_PROJ))
                call simple_chdir(here, status)
            enddo
            l_table = files_identical(string(trim(RUN_DIRS(1))//'/'//CLASS_STUDENT_PRED_TABLE), &
                &string(trim(RUN_DIRS(2))//'/'//CLASS_STUDENT_PRED_TABLE))
            l_cands = files_identical(string(trim(RUN_DIRS(1))//'/'//cands_fname%to_char()), &
                &string(trim(RUN_DIRS(2))//'/'//cands_fname%to_char()))
            l_clean = .true.
            do irun = 1,size(PART_FILES)
                if( file_exists(string(trim(RUN_DIRS(2))//'/'//trim(PART_FILES(irun)))) ) l_clean = .false.
            enddo
            ! the empty part: its table's header alone, its range's records and its candidates part file
            nlines = 0
            call fopen(u, string(trim(RUN_DIRS(3))//'/class_student_part1_pred.txt'), status='old', action='read', iostat=ios)
            if( ios == 0 )then
                do
                    read(u,'(A)', iostat=ios) line
                    if( ios /= 0 ) exit
                    nlines = nlines + 1
                enddo
                call fclose(u)
            endif
            l_empty_part = nlines == 1 .and. file_exists(string(trim(RUN_DIRS(3))//'/part1.simple')) .and. &
                &file_exists(string(trim(RUN_DIRS(3))//'/class_student_part1_cands.bin'))
            l_recs  = runs(1)%os_ptcl2D%get_noris() == runs(2)%os_ptcl2D%get_noris()
            if( l_recs )then
                do i = 1,runs(1)%os_ptcl2D%get_noris()
                    if( runs(1)%os_ptcl2D%ori2str(i) /= runs(2)%os_ptcl2D%ori2str(i) ) l_recs = .false.
                enddo
            endif
            call input%kill
            call runs(1)%kill
            call runs(2)%kill
            deallocate(input, runs)
        end subroutine master_as_shared_memory

        !> each read's value: the READ_COL-th number after its label, on the first line of its file that starts with the
        !! label (after leading blanks)
        subroutine read_reads( v, found )
            use simple_defs_hhcnn, only: VALIDATE_CLASS_STUDENT_REPORT, VALIDATE_CLASS_STUDENT_DETAILS
            real,    intent(out) :: v(NREADS)
            logical, intent(out) :: found(NREADS)
            character(len=512) :: line, rest
            character(len=64)  :: tok
            real    :: x
            integer :: u, ios, j, ifile, nnum, k, ios_tok
            v     = -1.
            found = .false.
            do ifile = 1,2
                if( ifile == 1 )then
                    call fopen(u, string(VALIDATE_CLASS_STUDENT_REPORT), status='old', action='read', iostat=ios)
                else
                    call fopen(u, string(VALIDATE_CLASS_STUDENT_DETAILS), status='old', action='read', iostat=ios)
                endif
                call fileiochk('test_class_student; open the validation report', ios)
                do
                    read(u,'(A)', iostat=ios) line
                    if( ios /= 0 ) exit
                    line = adjustl(line)
                    do j = 1,NREADS
                        if( found(j) .or. READ_FILE(j) /= ifile ) cycle
                        if( index(line, trim(READ_LABELS(j))) /= 1 ) cycle
                        ! the numbers after the label, in order; words between them are skipped
                        rest = line(len_trim(READ_LABELS(j))+1:)
                        nnum = 0
                        do while( len_trim(rest) > 0 .and. nnum < READ_COL(j) )
                            rest = adjustl(rest)
                            k    = index(rest, ' ')
                            tok  = rest(1:k-1)
                            rest = rest(k:)
                            read(tok, *, iostat=ios_tok) x
                            if( ios_tok /= 0 ) cycle
                            nnum = nnum + 1
                            if( nnum == READ_COL(j) )then
                                v(j)     = x
                                found(j) = .true.
                            endif
                        enddo
                    enddo
                enddo
                call fclose(u)
            enddo
        end subroutine read_reads

    end subroutine exec_test_class_student

end submodule simple_commanders_test_highlevel_class_student
