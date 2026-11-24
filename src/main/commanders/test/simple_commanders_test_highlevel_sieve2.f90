!@descr: particle_sieving_2 in miniature on 2 000 simulated beta-galactosidase particles in 20 stacks, a gate of simple_commanders_test_highlevel
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_sieve2
implicit none
#include "simple_local_flags.inc"

contains

    !> particle_sieving_2 at the miniature's sizes (a seed of five stacks ended by its particle cap, ten classes,
    !! batches of 300, a round every two, the student retrained after round 1) with the CNN fixture's weights at a
    !! frozen threshold: counts gated, reads reported to metrics.tsv; the run directory is kept for test=identity
    module subroutine exec_test_particle_sieving_2( self, cline )
        use simple_test_state_fixture, only: state_fixture, FIX_SEARCH
        use simple_test_gate,          only: test_gate
        use simple_ptcl_sieve2,        only: ptcl_sieve2, sieve2_settings
        use simple_defs_sieve2,        only: SIEVE2_LABELS_FBODY, SIEVE2_ROUND_DIR, SIEVE2_BATCH_DIR, SIEVE2_LEARN_DIR
        use simple_timer,              only: tic, toc, timer_int_kind
        use simple_sigma2_state_file,  only: sigma2_state_read_particles, SIGMA2_STATE_FNAME
        use, intrinsic :: iso_fortran_env, only: real32
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        class(commander_test_particle_sieving_2), intent(inout) :: self
        class(cmdline),                           intent(inout) :: cline
        integer,  parameter :: NPTCLS         = 2000
        integer,  parameter :: NSTKS          = 20     !< one stack per micrograph, 100 particles each
        integer,  parameter :: NSEED_MICS     = 5      !< the seed's stacks, ended by its particle cap
        integer,  parameter :: NMICS          = 10     !< nmics, which the cap comes before
        integer,  parameter :: NCLS           = 10
        integer,  parameter :: BATCH_NPTCLS   = 300    !< three stacks, so that every batch closes at it
        integer,  parameter :: ROUND_NBATCHES = 2
        integer,  parameter :: TRAIN_NPTCLS   = 400
        integer,  parameter :: NSEED_PTCLS    = NSEED_MICS * (NPTCLS / NSTKS)
        integer,  parameter :: SEED_NPTCLS    = NSEED_PTCLS - (NPTCLS / NSTKS) / 2 !< the cap, reached by the fifth stack
        integer,  parameter :: NBATCHES       = (NPTCLS - NSEED_PTCLS) / BATCH_NPTCLS
        integer,  parameter :: NROUNDS        = (NBATCHES + ROUND_NBATCHES - 1) / ROUND_NBATCHES
        real,     parameter :: SNR            = 0.5
        real,     parameter :: SHERR_PX       = 6.  !< the planted shifts' range (pixels at 2.2 A), the class_student test's
        real,     parameter :: HALVES_IMBALANCE = 0.05 !< the even/odd partition's own warning bound (IMBALANCE_WARN)
        !> frozen between the p_junk values the fixture network gives the miniature's classes, so both labels occur
        real(dp), parameter :: JUNK_THRES     = 0.503d0
        character(len=*), parameter :: WEIGHTS = SIMPLE_SOURCE_PATH//'/testdata/cavg_cnn/weights.raw'
        class(parameters),  allocatable :: params, params_run
        class(ptcl_sieve2), allocatable :: sieve
        type(sieve2_settings)   :: settings
        type(state_fixture)     :: fix
        type(test_gate)         :: gate
        type(sp_project)        :: proj
        type(cmdline)           :: cl
        type(string)            :: cwd_saved, root, projfile, projname, edir, fname
        integer(timer_int_kind) :: t0
        integer, allocatable    :: lab_pop(:)
        real(real32), allocatable :: spectra(:,:)
        character(len=STDLEN)   :: msg
        logical, allocatable    :: lab_bad(:)
        real    :: wall
        integer :: status, ib, ir, n, nsum, nshort, iptcl, icls, nacc, nwrong, ir_bad(NROUNDS), neo(2)
        logical :: l_ok, l_active, l_halves, l_found
        allocate(params)
        call simple_getcwd(cwd_saved)
        root = filepath(cwd_saved, 'test_particle_sieving_2_mini_'//int2str(get_process_id()))
        if( dir_exists(root) ) call simple_rmdir(root)
        call simple_mkdir(root)
        call simple_chdir(root, status)
        if( status /= 0 ) THROW_HARD('TEST_PARTICLE_SIEVING_2 FAILED: could not enter the test directory')
        call params%new(cline)
        call gate%new(filepath(cwd_saved, 'metrics.tsv'))
        ! the fixture: raw picks with planted shifts, its one stack cut into NSTKS (split_stk declines a stack that
        ! add_stk recorded as a split one), a micrograph row per stack and no halves, as an imported project: the
        ! program makes the halves by micrograph group
        call fix%new('betagal', FIX_SEARCH, 1, NPTCLS, SNR, params%nthr, kind='identical', sherr=SHERR_PX)
        call fix%build
        projfile = simple_abspath(fix%project_fname())
        projname = basename(projfile)
        call proj%read(projfile)
        call proj%os_stk%set(1, 'stkkind', 'single')
        call proj%split_stk(NSTKS)
        if( proj%os_stk%get_noris() /= NSTKS ) THROW_HARD('TEST_PARTICLE_SIEVING_2 FAILED: the fixture was not cut')
        call proj%os_mic%new(NSTKS, is_ptcl=.false.)
        do ib = 1,NSTKS
            call proj%os_mic%set(ib, 'nptcls', NPTCLS / NSTKS)
            call proj%os_mic%set_state(ib, 1)
        enddo
        call proj%os_stk%delete_entry('eo_byindex')
        call proj%os_stk%delete_entry('eo_group')
        do iptcl = 1,NPTCLS
            call proj%os_ptcl2D%set(iptcl, 'eo', 0.)
            call proj%os_ptcl3D%set(iptcl, 'eo', 0.)
        enddo
        call proj%write(projfile)
        call proj%kill
        ! the program at the miniature's sizes, in the test's directory (simple_test_exec makes no execution directory)
        call cl%set('prg',      'particle_sieving_2')
        call cl%set('projfile', projfile%to_char())
        call cl%set('infile',   WEIGHTS)
        call cl%set('nmics',    NMICS)
        call cl%set('ncls',     NCLS)
        call cl%set('nthr',     params%nthr)
        call cl%set('mkdir',    'yes')
        call cl%set('oritype',  'ptcl2D')
        allocate(params_run)
        call params_run%new(cl)
        settings%nmics          = NMICS
        settings%seed_nptcls    = SEED_NPTCLS
        settings%ncls           = NCLS
        settings%batch_nptcls   = BATCH_NPTCLS
        settings%round_nbatches = ROUND_NBATCHES
        settings%train_nptcls   = TRAIN_NPTCLS
        settings%junk_thres     = JUNK_THRES
        settings%mskdiam        = (real(params_run%box) - COSMSKHALFWIDTH) * params_run%smpd ! the box's mask
        settings%retrain_rounds = [1, NROUNDS] ! no batch follows the last round: retrained after round 1 only
        allocate(sieve)
        t0 = tic()
        call sieve%new(params_run, settings)
        call sieve%run
        call sieve%kill
        wall = real(toc(t0))
        call simple_getcwd(edir)
        ! the batches: each but the tail at least the batch size, every particle after the seed in one
        nsum    = 0
        nshort  = 0
        l_ok    = .true.
        l_halves = .true.
        do ib = 1,NBATCHES
            fname = SIEVE2_BATCH_DIR//int2str_pad(ib, 3)//'/'//projname%to_char()
            if( .not. file_exists(fname) )then
                l_ok     = .false.
                l_halves = .false.
                cycle
            endif
            call count_selected_halves(fname, n, neo)
            nsum = nsum + n
            if( ib < NBATCHES .and. n < BATCH_NPTCLS ) nshort = nshort + 1
            if( any(neo == 0) ) l_halves = .false.
        enddo
        fname = SIEVE2_BATCH_DIR//int2str_pad(NBATCHES + 1, 3)
        call gate%check('the_batches_cut_as_planned', l_ok .and. .not. dir_exists(fname))
        call gate%check('every_batch_but_the_tail_holds_the_batch_size', nshort == 0)
        call gate%check('every_particle_after_the_seed_in_a_batch', nsum == NPTCLS - NSEED_PTCLS)
        call gate%check('every_batch_holds_both_halves', l_halves)
        ! the rounds: one every ROUND_NBATCHES batches and after the last, each on every particle ingested by then
        l_ok     = .true.
        l_active = .true.
        do ir = 1,NROUNDS
            fname = SIEVE2_ROUND_DIR//int2str_pad(ir, 2)//'/'//projname%to_char()
            if( .not. file_exists(fname) )then
                l_ok     = .false.
                l_active = .false.
                cycle
            endif
            call count_selected_halves(fname, n, neo)
            if( n /= min(NPTCLS, NSEED_PTCLS + ir * ROUND_NBATCHES * BATCH_NPTCLS) ) l_active = .false.
        enddo
        fname = SIEVE2_ROUND_DIR//int2str_pad(NROUNDS + 1, 2)
        call gate%check('one_greedy_round_per_round_nbatches_and_at_the_end', l_ok .and. .not. dir_exists(fname))
        call gate%check('every_round_holds_every_ingested_particle_active', l_active)
        ! the student: trained after the seed and again after round 1, the one listed round a batch follows
        l_ok = dir_exists(string(SIEVE2_LEARN_DIR//'00')) .and. dir_exists(string(SIEVE2_LEARN_DIR//'01'))
        do ir = 2,NROUNDS
            if( dir_exists(string(SIEVE2_LEARN_DIR//int2str_pad(ir, 2))) ) l_ok = .false.
        enddo
        call gate%check('the_student_trained_after_the_seed_and_the_listed_rounds', l_ok)
        ! the last round's sigma2 store: a valid record (every shell finite and positive) for every particle, the
        ! batches' particles seeded from the global curves and then updated over the round's band
        fname = SIEVE2_ROUND_DIR//int2str_pad(NROUNDS, 2)//'/'//SIGMA2_STATE_FNAME
        call sigma2_state_read_particles(fname%to_char(), 1, NPTCLS, spectra, status, msg)
        l_ok = status == 0
        if( l_ok ) l_ok = all(ieee_is_finite(spectra))
        if( l_ok ) l_ok = all(spectra > 0._real32)
        call gate%check('every_particle_has_a_valid_sigma2_record', l_ok)
        ! the result: accepted exactly the particles whose class the last labels keep, every particle decided
        do ir = 1,NROUNDS
            call read_labels(SIEVE2_LABELS_FBODY//'round_'//int2str_pad(ir, 2)//'.txt', lab_pop, lab_bad, l_found)
            ir_bad(ir) = -1
            if( l_found ) ir_bad(ir) = count(lab_bad)
        enddo
        call read_labels(SIEVE2_LABELS_FBODY//'round_'//int2str_pad(NROUNDS, 2)//'.txt', lab_pop, lab_bad, l_found)
        call gate%check('the_last_labels_are_written', l_found)
        call proj%read(projname)
        nacc   = 0
        nwrong = 0
        do iptcl = 1,proj%os_ptcl2D%get_noris()
            icls = proj%os_ptcl2D%get_class(iptcl)
            l_ok = .false.
            if( l_found .and. icls >= 1 .and. icls <= NCLS ) l_ok = .not. lab_bad(icls)
            if( proj%os_ptcl2D%get_state(iptcl) == 1 ) nacc = nacc + 1
            if( (proj%os_ptcl2D%get_state(iptcl) == 1) .neqv. l_ok ) nwrong = nwrong + 1
            if( proj%os_ptcl3D%get_state(iptcl) /= proj%os_ptcl2D%get_state(iptcl) ) nwrong = nwrong + 1
        enddo
        call gate%check('the_result_holds_every_particle', proj%os_ptcl2D%get_noris() == NPTCLS)
        neo = 0
        do iptcl = 1,proj%os_ptcl2D%get_noris()
            if( proj%os_ptcl2D%get_eo(iptcl) == 0 ) neo(1) = neo(1) + 1
            if( proj%os_ptcl2D%get_eo(iptcl) == 1 ) neo(2) = neo(2) + 1
        enddo
        call gate%check('the_halves_hold_every_particle_balanced', sum(neo) == NPTCLS .and. &
            &real(abs(neo(1) - neo(2))) <= HALVES_IMBALANCE * real(NPTCLS))
        call gate%check('accepted_exactly_the_particles_of_classes_not_junk', l_found .and. nwrong == 0)
        call gate%check('both_labels_occur_at_the_end', l_found .and. any(lab_bad .and. lab_pop > 0) .and. &
            &any(.not. lab_bad .and. lab_pop > 0))
        call proj%kill
        ! the reads
        call gate%report('accepted_particles', real(nacc))
        do ir = 1,NROUNDS
            call gate%report('junk_classes_round_'//int2str(ir), real(ir_bad(ir)))
        enddo
        call gate%report('particle_sieving_2_wall_s', wall)
        call fix%kill
        call simple_chdir(cwd_saved, status)
        if( status /= 0 ) THROW_HARD('TEST_PARTICLE_SIEVING_2 FAILED: could not restore the original directory')
        if( gate%passed() )then
            call gate%kill
            write(logfhandle,'(A)') 'PASS: particle_sieving_2 in miniature'
            call simple_end('**** SIMPLE_TEST_PARTICLE_SIEVING_2 NORMAL STOP ****')
        else
            call gate%kill
            THROW_HARD('TEST_PARTICLE_SIEVING_2 FAILED')
        endif

      contains

        !> the selected particles of a step's project and how many of them are in each half
        subroutine count_selected_halves( pname, nsel, nhalves )
            class(string), intent(in)  :: pname
            integer,       intent(out) :: nsel, nhalves(2)
            type(sp_project) :: p
            integer :: i, eo
            call p%read_segment('ptcl2D', pname)
            nsel    = 0
            nhalves = 0
            do i = 1,p%os_ptcl2D%get_noris()
                if( p%os_ptcl2D%get_state(i) == 0 ) cycle
                nsel = nsel + 1
                eo   = p%os_ptcl2D%get_eo(i)
                if( eo == 0 .or. eo == 1 ) nhalves(eo+1) = nhalves(eo+1) + 1
            enddo
            call p%kill
        end subroutine count_selected_halves

        !> a labels table's populations and junk flags, by class
        subroutine read_labels( tname, pops, bad, found )
            character(len=*),     intent(in)  :: tname
            integer, allocatable, intent(out) :: pops(:)
            logical, allocatable, intent(out) :: bad(:)
            logical,              intent(out) :: found
            character(len=256) :: line
            real(dp) :: pj, ps, pl
            integer  :: u, ios, c, pop, ie, ibad
            allocate(pops(NCLS), source=0)
            allocate(bad(NCLS), source=.false.)
            found = file_exists(tname)
            if( .not. found ) return
            call fopen(u, string(tname), status='old', action='read', iostat=ios)
            call fileiochk('read_labels; open '//tname, ios)
            do
                read(u,'(A)', iostat=ios) line
                if( ios /= 0 ) exit
                if( line(1:1) == '#' ) cycle
                read(line,*) c, pop, pj, ps, pl, ie, ibad
                if( c < 1 .or. c > NCLS ) cycle
                pops(c) = pop
                bad(c)  = ibad == 1
            enddo
            call fclose(u)
        end subroutine read_labels

    end subroutine exec_test_particle_sieving_2

end submodule simple_commanders_test_highlevel_sieve2
