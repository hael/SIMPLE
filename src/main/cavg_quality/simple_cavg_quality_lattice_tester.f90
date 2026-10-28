!@descr: unit tests of the nano context of the class-average quality backend: the lattice family, its normalization, the gates, the preset, the bank by context, model and table round trips
module simple_cavg_quality_lattice_tester
!$ use omp_lib
use simple_test_utils
use simple_defs,                   only: dp, DPI
use simple_string,                 only: string
use simple_string_utils,           only: int2str
use simple_fileio,                 only: file_exists, del_file
use simple_syslib,                 only: get_process_id
use simple_image,                  only: image
use simple_oris,                   only: oris
use simple_cavg_quality_types,     only: CAVG_QUALITY_NFEATS, CAVG_QUALITY_CONTEXT_NANO, CAVG_QUALITY_CONTEXT_CHUNK, &
    CAVG_RELATIONAL_SCHEMA_CORR_KNN_SIGNAL_V1, CAVG_RELATIONAL_NANO_CORR_HP, CAVG_RELATIONAL_NANO_CORR_LP, &
    CAVG_RELATIONAL_NANO_CORR_TRS, CLIP_Z, cavg_quality_result
use simple_cavg_quality_lattice,   only: cavg_quality_lattice, LAT_NFEATS, LAT_LATTICE_PEAK, LAT_SPOT_COUNT, &
    LAT_SPOT_ARC, LAT_SPOT_DIAM, LAT_BAND_PRESENCE, LAT_PEAK_MIN
use simple_cavg_quality_feats,     only: cavg_quality_bank_size, cavg_quality_feature_name, cavg_quality_feature_family
use simple_cavg_quality_analysis,  only: evaluate_cavg_quality_hard_reject, write_cavg_quality_training_table
use simple_cavg_quality_model,     only: cavg_quality_model, cavg_quality_model_spec, write_cavg_quality_model_builtin_code
use simple_cavg_quality_learn,     only: learn_cavg_quality_model
use simple_cavg_quality_relations, only: cavg_quality_relation_analysis
implicit none
private
public :: run_all_cavg_quality_lattice_tests

integer,  parameter :: BOX  = 120
real,     parameter :: SMPD = 0.358
integer,  parameter :: FH = 18, FK = 12           !< the fringe wave vector (shell 22, 1.99 A at this box)
integer,  parameter :: KIND_NOISE = 0, KIND_LATTICE = 1, KIND_LOWPASS = 2
integer(kind=8)     :: rng_state = 20261028_8

contains

    subroutine run_all_cavg_quality_lattice_tests()
        write(*,'(A)') '**** running all class-average nano context tests ****'
        call test_bank_by_context()
        call test_lattice_against_noise()
        call test_normalization()
        call test_nano_preset()
        call test_nano_gates()
        call test_diameter_calibration()
        call test_spot_arc_straight_bent()
        call test_band_presence_unmasked()
        call test_team_identity()
        call test_model_file_round_trips()
        call test_table_round_trip_learn()
    end subroutine run_all_cavg_quality_lattice_tests

    !> the cryo-EM contexts keep their fourteen features; nano has the six lattice features
    subroutine test_bank_by_context()
        character(len=16), parameter :: NANO_NAMES(LAT_NFEATS) = [character(len=16) :: 'lattice_peak', &
            'spot_count', 'spot_arc', 'spot_diam', 'spot_snr', 'band_presence']
        integer :: i
        write(*,'(A)') 'test_bank_by_context'
        call assert_int(CAVG_QUALITY_NFEATS, cavg_quality_bank_size('chunk'), 'chunk bank size')
        call assert_int(CAVG_QUALITY_NFEATS, cavg_quality_bank_size('pool'),  'pool bank size')
        call assert_int(CAVG_QUALITY_NFEATS, cavg_quality_bank_size('sieve'), 'sieve bank size')
        call assert_int(LAT_NFEATS,          cavg_quality_bank_size('nano'),  'nano bank size')
        call assert_char('log_pop',           trim(cavg_quality_feature_name(1, 'chunk')),  'first chunk feature')
        call assert_char('fuzzy_ball_signal', trim(cavg_quality_feature_name(14, 'pool')), 'last pool feature')
        do i = 1, LAT_NFEATS
            call assert_char(trim(NANO_NAMES(i)), trim(cavg_quality_feature_name(i, 'nano')), 'nano feature name')
            call assert_char('lattice', trim(cavg_quality_feature_family(i, 'nano')), 'nano feature family')
        end do
    end subroutine test_bank_by_context

    !> a lattice disc stands far above noise; the same lattice low-passed (its envelope only) falls to the noise level
    subroutine test_lattice_against_noise()
        type(cavg_quality_lattice) :: lat
        type(image)                :: imgs(3)
        real    :: raw(3,LAT_NFEATS)
        logical :: blank(3)
        write(*,'(A)') 'test_lattice_against_noise'
        call make_cavg(imgs(1), KIND_LATTICE, 35., 0.3)
        call make_cavg(imgs(2), KIND_LOWPASS, 35., 0.3)
        call make_cavg(imgs(3), KIND_NOISE,   35., 0.3)
        call lat%new(BOX, SMPD)
        call lat%extract(imgs, [.true., .true., .true.], raw, blank)
        call assert_false(any(blank), 'no blank class')
        call assert_true(raw(1,LAT_LATTICE_PEAK) > raw(3,LAT_LATTICE_PEAK) + 0.5, 'lattice above noise')
        call assert_true(abs(raw(2,LAT_LATTICE_PEAK) - raw(3,LAT_LATTICE_PEAK)) < 0.3, 'low-passed lattice at noise')
        call lat%kill
        call kill_imgs(imgs)
    end subroutine test_lattice_against_noise

    !> the nano normalization: lattice_peak enters raw (lattice classes above LAT_PEAK_MIN, empty ones below, in a
    !! mixed project as in an all-lattice or an all-empty one, where a floor would not exist), the other features by
    !! median and MAD over the valid classes, hard-rejected rows at -CLIP_Z
    subroutine test_normalization()
        type(cavg_quality_result) :: q
        type(image), allocatable  :: imgs(:)
        integer, parameter :: N = 16, NU = 8
        integer :: kinds(N), ukinds(NU), j
        real    :: med
        write(*,'(A)') 'test_normalization'
        kinds(1:N/2)   = KIND_LATTICE
        kinds(N/2+1:N) = KIND_NOISE
        call project_quality(kinds, q, imgs)
        call assert_int(0, count(q%hard_reject), 'mixed project: no hard reject')
        call assert_true(all(q%features(:,LAT_LATTICE_PEAK) == q%raw(:,LAT_LATTICE_PEAK)), 'mixed: lattice_peak enters raw')
        call assert_true(all(q%raw(1:N/2,LAT_LATTICE_PEAK) > LAT_PEAK_MIN), 'mixed: lattice classes above the lattice grade')
        call assert_true(all(q%raw(N/2+1:N,LAT_LATTICE_PEAK) < LAT_PEAK_MIN), 'mixed: empty classes below the lattice grade')
        call assert_true(all(q%raw(1:N/2,LAT_SPOT_COUNT) >= 1.) .and. all(q%raw(N/2+1:N,LAT_SPOT_COUNT) == 0.), &
            'mixed: spot_count counts lattice-grade maxima only')
        do j = 1, LAT_NFEATS
            if( j == LAT_LATTICE_PEAK ) cycle
            med = median_of(q%features(:,j))
            call assert_true(abs(med) < 1.e-4 .or. all(q%features(:,j) == 0.), 'mixed: median/MAD centres the other features')
        end do
        call q%kill
        call kill_imgs(imgs)
        ukinds = KIND_LATTICE
        call project_quality(ukinds, q, imgs)
        call assert_int(0, count(q%hard_reject), 'all-lattice project: no hard reject')
        call assert_true(all(q%features(:,LAT_LATTICE_PEAK) > LAT_PEAK_MIN), 'all-lattice project: every class above the grade')
        call q%kill
        call kill_imgs(imgs)
        ukinds = KIND_NOISE
        call project_quality(ukinds, q, imgs)
        call assert_int(0, count(q%hard_reject), 'all-empty project: no hard reject')
        call assert_true(all(q%features(:,LAT_LATTICE_PEAK) < LAT_PEAK_MIN), 'all-empty project: every class below the grade')
        call q%kill
        call kill_imgs(imgs)

        contains

            real function median_of( x )
                real, intent(in) :: x(:)
                real    :: xs(size(x)), t
                integer :: a, b, n
                xs = x
                n  = size(xs)
                do a = 2, n
                    t = xs(a)
                    b = a - 1
                    do while( b >= 1 )
                        if( xs(b) <= t ) exit
                        xs(b+1) = xs(b)
                        b = b - 1
                    end do
                    xs(b+1) = t
                end do
                if( mod(n, 2) == 1 )then
                    median_of = xs(n/2 + 1)
                else
                    median_of = 0.5 * (xs(n/2) + xs(n/2 + 1))
                endif
            end function median_of

    end subroutine test_normalization

    !> the nano preset is the lattice gate: context nano, no relational term, a threshold on lattice_peak at
    !! LAT_PEAK_MIN (probability 1/2 there), movable by set_nano_threshold; on a mixed project it keeps the lattice
    !! classes and rejects the empty ones, and it decides the same on an all-lattice and an all-empty project
    subroutine test_nano_preset()
        type(cavg_quality_model)      :: m
        type(cavg_quality_model_spec) :: spec
        type(cavg_quality_result)     :: q
        type(image), allocatable      :: imgs(:)
        integer, parameter :: N = 16, NU = 8
        integer :: kinds(N), ukinds(NU)
        real    :: x(LAT_NFEATS), thr, p_at, p_below, p_above
        write(*,'(A)') 'test_nano_preset'
        call m%init_preset('nano')
        spec = m%get_spec()
        call assert_char(CAVG_QUALITY_CONTEXT_NANO, trim(spec%context), 'nano preset context')
        call assert_char('lattice', trim(spec%feature_policy), 'nano preset policy')
        call assert_false(m%supports_relational(), 'nano preset has no relational term')
        call assert_true(spec%linear_coefficients(LAT_LATTICE_PEAK) > 0. .and. &
            &all(abs(spec%linear_coefficients(2:LAT_NFEATS)) < 1.e-6), 'nano preset reads lattice_peak only')
        thr = -spec%intercept / spec%linear_coefficients(LAT_LATTICE_PEAK)
        call assert_true(abs(thr - LAT_PEAK_MIN) < 1.e-4, 'nano preset threshold: lattice_peak at LAT_PEAK_MIN')
        x = 0.
        x(LAT_LATTICE_PEAK) = thr
        p_at = m_probability(m, x)
        x(LAT_LATTICE_PEAK) = thr - 0.1
        p_below = m_probability(m, x)
        x(LAT_LATTICE_PEAK) = thr + 0.1
        p_above = m_probability(m, x)
        call assert_true(abs(p_at - 0.5) < 1.e-5 .and. p_below < 0.5 .and. p_above > 0.5, &
            'nano preset: probability 1/2 at the threshold, below and above it on either side')
        kinds(1:N/2)   = KIND_LATTICE
        kinds(N/2+1:N) = KIND_NOISE
        call project_quality(kinds, q, imgs)
        call m%classify(q)
        call assert_true(gate_is_decision(q, thr), 'nano preset: the decision is the lattice_peak threshold, class by class')
        call assert_true(all(q%states(1:N/2) == 1), 'nano preset keeps the lattice classes')
        call assert_true(all(q%states(N/2+1:N) == 0), 'nano preset rejects the empty classes')
        ! thres: the threshold moved above the lattice classes rejects them too; back at the default it keeps them
        call m%set_nano_threshold(2.5)
        call m%classify(q)
        call assert_true(gate_is_decision(q, 2.5) .and. all(q%states == 0), 'set_nano_threshold moves the gate')
        call m%init_preset('nano')
        call m%classify(q)
        call assert_true(all(q%states(1:N/2) == 1), 'nano preset restored')
        call q%kill
        call kill_imgs(imgs)
        ukinds = KIND_LATTICE
        call project_quality(ukinds, q, imgs)
        call m%classify(q)
        call assert_true(gate_is_decision(q, thr) .and. all(q%states == 1), 'nano preset keeps an all-lattice project')
        call q%kill
        call kill_imgs(imgs)
        ukinds = KIND_NOISE
        call project_quality(ukinds, q, imgs)
        call m%classify(q)
        call assert_true(gate_is_decision(q, thr) .and. all(q%states == 0), 'nano preset rejects an all-empty project')
        call q%kill
        call kill_imgs(imgs)

        contains

            !> the probability the preset gives one class with bank values xc (lattice_peak raw)
            real function m_probability( model, xc )
                type(cavg_quality_model), intent(in) :: model
                real,                     intent(in) :: xc(:)
                type(cavg_quality_result) :: one
                allocate(one%features(1, size(xc)), source=0.)
                one%features(1,:) = xc
                allocate(one%raw(1, size(xc)), source=0.)
                one%raw(1,:) = xc
                allocate(one%hard_reject(1), source=.false.)
                call model%classify(one)
                m_probability = one%scores(1)
                call one%kill
            end function m_probability

            !> every class is accepted exactly when it passes the validity gates and its lattice_peak >= thr_c
            logical function gate_is_decision( qc, thr_c )
                type(cavg_quality_result), intent(in) :: qc
                real,                      intent(in) :: thr_c
                integer :: ic
                gate_is_decision = .true.
                do ic = 1, size(qc%states)
                    if( qc%hard_reject(ic) )then
                        gate_is_decision = gate_is_decision .and. qc%states(ic) == 0
                    else if( qc%raw(ic,LAT_LATTICE_PEAK) >= thr_c )then
                        gate_is_decision = gate_is_decision .and. qc%states(ic) == 1
                    else
                        gate_is_decision = gate_is_decision .and. qc%states(ic) == 0
                    endif
                end do
            end function gate_is_decision

    end subroutine test_nano_preset

    !> nano keeps the validity gates only: an empty class and a flat class are rejected, a lattice class whose
    !! foreground geometry would fail the cryo-EM gates is not
    subroutine test_nano_gates()
        type(cavg_quality_result) :: q
        type(image), allocatable  :: imgs(:)
        type(oris)                :: os
        integer :: i
        write(*,'(A)') 'test_nano_gates'
        allocate(imgs(6))
        do i = 1, 6
            call make_cavg(imgs(i), merge(KIND_LATTICE, KIND_NOISE, i <= 3), 55., 0.3)
        end do
        imgs(6) = 0.
        call class_oris(os, 6)
        call os%set(5, 'pop', 0.)
        call evaluate_cavg_quality_hard_reject(imgs, os, 0., q, CAVG_QUALITY_CONTEXT_NANO)
        call assert_int(LAT_NFEATS, size(q%features, dim=2), 'nano features')
        call assert_false(any(q%hard_reject(1:4)), 'valid classes pass')
        call assert_true(q%hard_reject(5), 'empty class rejected')
        call assert_true(q%hard_reject(6), 'flat class rejected')
        call assert_true(all(q%features(5:6,:) <= -CLIP_Z + 1.e-6), 'rejected rows at -CLIP_Z')
        call evaluate_cavg_quality_hard_reject(imgs, os, 0., q, CAVG_QUALITY_CONTEXT_CHUNK)
        call assert_int(CAVG_QUALITY_NFEATS, size(q%features, dim=2), 'chunk features')
        call q%kill
        call os%kill
        call kill_imgs(imgs)
    end subroutine test_nano_gates

    !> the projected-sphere fit of the strongest spot recovers the diameter of a lattice-filled sphere and,
    !! less closely, of a lattice-filled disc (25 A)
    subroutine test_diameter_calibration()
        type(cavg_quality_lattice) :: lat
        type(image)                :: imgs(2)
        real    :: raw(2,LAT_NFEATS)
        logical :: blank(2)
        write(*,'(A)') 'test_diameter_calibration'
        call make_cavg(imgs(1), KIND_LATTICE, 35., 0.02, sphere=.true.)
        call make_cavg(imgs(2), KIND_LATTICE, 35., 0.02)
        call lat%new(BOX, SMPD)
        call lat%extract(imgs, [.true., .true.], raw, blank)
        write(*,'(A,2F8.2)') 'spot_diam sphere / disc (A): ', raw(1,LAT_SPOT_DIAM), raw(2,LAT_SPOT_DIAM)
        call assert_real(25., raw(1,LAT_SPOT_DIAM), 2.5, 'sphere diameter within 10%')
        call assert_real(25., raw(2,LAT_SPOT_DIAM), 5.0, 'disc diameter within 20%')
        call lat%kill
        call kill_imgs(imgs)
    end subroutine test_diameter_calibration

    !> straight fringes give a round spot; fringes bent about a centre outside the particle smear it along the ring
    subroutine test_spot_arc_straight_bent()
        type(cavg_quality_lattice) :: lat
        type(image)                :: imgs(2)
        real    :: raw(2,LAT_NFEATS)
        logical :: blank(2)
        write(*,'(A)') 'test_spot_arc_straight_bent'
        call make_cavg(imgs(1), KIND_LATTICE, 35., 0.02)
        call make_cavg(imgs(2), KIND_LATTICE, 35., 0.02, bent=.true.)
        call lat%new(BOX, SMPD)
        call lat%extract(imgs, [.true., .true.], raw, blank)
        write(*,'(A,2F8.3)') 'spot_arc straight / bent: ', raw(1,LAT_SPOT_ARC), raw(2,LAT_SPOT_ARC)
        call assert_true(raw(2,LAT_SPOT_ARC) > raw(1,LAT_SPOT_ARC) + 0.1, 'bent fringes widen the spot along the ring')
        call lat%kill
        call kill_imgs(imgs)
    end subroutine test_spot_arc_straight_bent

    !> band presence is measured on the unmasked image: a lattice at the centre is positive, a lattice ring between
    !! a quarter and a half of the box negative
    subroutine test_band_presence_unmasked()
        type(cavg_quality_lattice) :: lat
        type(image)                :: imgs(2)
        real    :: raw(2,LAT_NFEATS)
        logical :: blank(2)
        write(*,'(A)') 'test_band_presence_unmasked'
        call make_cavg(imgs(1), KIND_LATTICE, 25., 0.05)
        call make_cavg(imgs(2), KIND_LATTICE, 25., 0.05, ring=.true.)
        call lat%new(BOX, SMPD)
        call lat%extract(imgs, [.true., .true.], raw, blank)
        call assert_true(raw(1,LAT_BAND_PRESENCE) > 0.3,  'central lattice: positive presence')
        call assert_true(raw(2,LAT_BAND_PRESENCE) < -0.3, 'lattice ring: negative presence')
        call lat%kill
        call kill_imgs(imgs)
    end subroutine test_band_presence_unmasked

    !> the parallel extraction gives the same bits with a team as with one thread
    subroutine test_team_identity()
        type(cavg_quality_lattice) :: lat
        type(image),  allocatable  :: imgs(:)
        real,         allocatable  :: raw1(:,:), rawn(:,:)
        logical,      allocatable  :: valid(:), blank(:)
        integer, parameter :: N = 12
        integer :: i, nthr
        write(*,'(A)') 'test_team_identity'
        allocate(imgs(N), raw1(N,LAT_NFEATS), rawn(N,LAT_NFEATS), valid(N), blank(N))
        valid = .true.
        do i = 1, N
            call make_cavg(imgs(i), merge(KIND_LATTICE, KIND_NOISE, mod(i,2) == 0), 35., 0.3)
        end do
        call lat%new(BOX, SMPD)
        nthr = 1
        !$ nthr = omp_get_max_threads()
        !$ call omp_set_num_threads(1)
        call lat%extract(imgs, valid, raw1, blank)
        !$ call omp_set_num_threads(nthr)
        call lat%extract(imgs, valid, rawn, blank)
        call assert_true(all(raw1 == rawn), 'raw lattice features: team == one thread')
        call lat%kill
        call kill_imgs(imgs)
    end subroutine test_team_identity

    !> a nano model survives write and read; a file naming a preset without its interactions keeps the preset's
    subroutine test_model_file_round_trips()
        type(cavg_quality_model)      :: m, r, p
        type(cavg_quality_model_spec) :: spec
        type(string) :: fname, snippet
        integer :: i, j, iterm, funit
        write(*,'(A)') 'test_model_file_round_trips'
        fname   = 'nano_model_test_'//int2str(get_process_id())//'.txt'
        snippet = 'nano_model_snippet_'//int2str(get_process_id())//'.txt'
        call m%init_preset('chunk100mics')
        spec = m%get_spec()
        spec%name           = 'nano_test'
        spec%context        = CAVG_QUALITY_CONTEXT_NANO
        spec%feature_policy = 'lattice'
        spec%weights        = 0.
        spec%weights(1:LAT_NFEATS) = 1. / real(LAT_NFEATS)
        spec%linear_coefficients = 0.
        spec%linear_coefficients(1:LAT_NFEATS) = [2.5, 0.25, -0.5, 0.125, 1.5, 0.75]
        spec%intercept      = -1.25
        spec%n_interactions = 0
        spec%interaction_terms        = 0
        spec%interaction_coefficients = 0.
        iterm = 0
        do i = 1, LAT_NFEATS - 1
            do j = i + 1, LAT_NFEATS
                iterm = iterm + 1
                spec%interaction_terms(iterm,:)    = [i, j]
                spec%interaction_coefficients(iterm) = 0.0625 * real(iterm)
            end do
        end do
        spec%n_interactions      = iterm
        spec%relational_knn      = 5
        spec%relational_corr_hp  = CAVG_RELATIONAL_NANO_CORR_HP
        spec%relational_corr_lp  = CAVG_RELATIONAL_NANO_CORR_LP
        spec%relational_corr_trs = CAVG_RELATIONAL_NANO_CORR_TRS
        spec%relational_coefficient = 0.5
        spec%prob_threshold      = 0.4
        call m%init_spec(spec)
        call m%write(fname%to_char())
        call r%read(fname%to_char())
        call assert_char('nano', trim(r%context), 'nano context read back')
        call assert_true(all(r%linear_coefficients(1:LAT_NFEATS) == m%linear_coefficients(1:LAT_NFEATS)), &
            'linear coefficients read back')
        call assert_true(all(r%linear_coefficients(LAT_NFEATS+1:) == 0.), 'no coefficient beyond the nano bank')
        call assert_int(15, r%n_interactions, 'nano interactions read back')
        call assert_true(all(r%interaction_terms(1:15,:) == m%interaction_terms(1:15,:)), 'interaction terms read back')
        call assert_real(m%relational_corr_lp, r%relational_corr_lp, 1.e-6, 'nano relational band read back')
        call write_cavg_quality_model_builtin_code(r, snippet%to_char())
        call assert_true(file_exists(snippet), 'promotion snippet written')
        call del_file(snippet)
        ! a preset without its interaction lines
        open(newunit=funit, file=fname%to_char(), status='replace', action='write')
        write(funit,'(A)') 'model_version=11'
        write(funit,'(A)') 'preset=sieve'
        write(funit,'(A)') 'relational_feature_schema=corr_knn_signal_v1'
        write(funit,'(A)') 'relational_corr_lp=8.0'
        close(funit)
        call p%init_preset('sieve')
        call r%read(fname%to_char())
        call assert_int(p%n_interactions, r%n_interactions, 'preset interactions kept')
        call assert_real(8.0, r%relational_corr_lp, 1.e-6, 'preset override applied')
        call del_file(fname)
    end subroutine test_model_file_round_trips

    !> analysis tables of both bank sizes go through learn: the model has the tables' context and bank
    subroutine test_table_round_trip_learn()
        type(string) :: tabs(2), fmodel
        type(cavg_quality_model) :: learned, r
        character(len=5) :: ctx(2) = ['nano ', 'chunk']
        integer :: ic, it
        write(*,'(A)') 'test_table_round_trip_learn'
        do ic = 1, 2
            do it = 1, 2
                tabs(it) = 'cq_table_'//trim(ctx(ic))//int2str(it)//'_'//int2str(get_process_id())//'.txt'
                call write_synthetic_table(trim(ctx(ic)), it, tabs(it)%to_char())
            end do
            fmodel = 'cq_learned_'//trim(ctx(ic))//'_'//int2str(get_process_id())//'.txt'
            call learn_cavg_quality_model(tabs, learned, fmodel%to_char(), 'cq_learn_report_'//&
                int2str(get_process_id())//'.txt')
            call r%read(fmodel%to_char())
            call assert_char(trim(ctx(ic)), trim(r%context), 'learned context = table context')
            call assert_true(all(r%interaction_terms(1:r%n_interactions,:) <= cavg_quality_bank_size(ctx(ic))), &
                'learned interactions within the bank')
            if( ic == 1 ) call assert_char('lattice', trim(r%feature_policy), 'nano policy')
            call del_file(fmodel)
            call del_file('cq_learn_report_'//int2str(get_process_id())//'.txt')
            do it = 1, 2
                call del_file(tabs(it))
            end do
        end do
    end subroutine test_table_round_trip_learn

    ! a training table of a context: 40 classes, the first feature separates the labels, one hard reject
    subroutine write_synthetic_table( context, iseed, fname )
        character(len=*), intent(in) :: context, fname
        integer,          intent(in) :: iseed
        type(cavg_quality_result) :: q
        type(cavg_quality_model)  :: m
        type(cavg_quality_model_spec) :: spec
        type(cavg_quality_relation_analysis) :: rel
        integer, parameter :: N = 40
        integer :: refs(N), i, nf
        nf = cavg_quality_bank_size(context)
        rng_state = 7_8 + int(iseed, kind=8)
        allocate(q%raw(N,nf), q%features(N,nf), q%scores(N), q%states(N), q%labels(N), q%hard_reject(N), &
            q%reasons(N))
        do i = 1, N
            refs(i) = merge(1, 0, mod(i,2) == 0)
        end do
        do i = 1, N
            q%features(i,:) = 0.5 * real(gauss_vec(nf))
            q%features(i,1) = q%features(i,1) + 1.5 * real(2 * refs(i) - 1)
        end do
        q%raw         = q%features
        q%hard_reject = .false.
        q%hard_reject(1) = .true.
        q%features(1,:) = -CLIP_Z
        q%states = 1
        q%labels = 1
        q%scores = 1.
        q%reasons = 0
        allocate(rel%raw(N), rel%feature(N))
        rel%feature = 0.3 * real(gauss_vec(N))
        rel%raw     = rel%feature
        rel%k_requested = 5
        rel%hp  = 15.
        rel%lp  = 1.5
        rel%trs = 3.
        call m%init_preset('chunk100mics')
        spec = m%get_spec()
        spec%context = trim(context)
        if( trim(context) == CAVG_QUALITY_CONTEXT_NANO )then
            spec%n_interactions = 0
            spec%interaction_terms = 0
            spec%interaction_coefficients = 0.
        endif
        call m%init_spec(spec)
        call write_cavg_quality_training_table(q, refs, m, fname, 'synthetic_'//int2str(iseed), context, rel)
        call q%kill
        call rel%kill
    end subroutine write_synthetic_table


    ! the nano-context quality of a project of synthetic class averages of the given kinds
    subroutine project_quality( kinds, q, imgs )
        integer,                   intent(in)    :: kinds(:)
        type(cavg_quality_result), intent(inout) :: q
        type(image), allocatable,  intent(inout) :: imgs(:)
        type(oris) :: os
        integer :: i
        if( allocated(imgs) ) deallocate(imgs)
        allocate(imgs(size(kinds)))
        do i = 1, size(kinds)
            call make_cavg(imgs(i), kinds(i), 35., 0.3)
        end do
        call class_oris(os, size(kinds))
        call evaluate_cavg_quality_hard_reject(imgs, os, 0., q, CAVG_QUALITY_CONTEXT_NANO)
        call os%kill
    end subroutine project_quality

    subroutine class_oris( os, n )
        type(oris), intent(inout) :: os
        integer,    intent(in)    :: n
        integer :: i
        call os%new(n, is_ptcl=.false.)
        do i = 1, n
            call os%set(i, 'class', real(i))
            call os%set(i, 'pop',   20.)
            call os%set(i, 'res',   12.)
            call os%set(i, 'corr',  0.2)
            call os%set(i, 'state', 1.)
        end do
    end subroutine class_oris

    ! a synthetic class average: a particle of radius rad_px (a disc, or a projected sphere) carrying fringes at
    ! (FH,FK) (straight, bent, or low-passed away), plus white noise of sd noise_sd; ring puts the fringes on an
    ! annulus instead of a disc
    subroutine make_cavg( img, kind, rad_px, noise_sd, sphere, bent, ring )
        type(image),       intent(inout) :: img
        integer,           intent(in)    :: kind
        real,              intent(in)    :: rad_px, noise_sd
        logical, optional, intent(in)    :: sphere, bent, ring
        real    :: rmat(BOX,BOX,1)
        real(dp):: x, y, r, profile, phase, kx, ky, kmag
        logical :: l_sphere, l_bent, l_ring
        integer :: i, j
        l_sphere = .false.; if( present(sphere) ) l_sphere = sphere
        l_bent   = .false.; if( present(bent)   ) l_bent   = bent
        l_ring   = .false.; if( present(ring)   ) l_ring   = ring
        kx   = 2._dp * DPI * real(FH, dp) / real(BOX, dp)
        ky   = 2._dp * DPI * real(FK, dp) / real(BOX, dp)
        kmag = hypot(kx, ky)
        do j = 1, BOX
            do i = 1, BOX
                x = real(i - 1 - BOX/2, dp)
                y = real(j - 1 - BOX/2, dp)
                r = hypot(x, y)
                if( l_ring )then
                    profile = merge(1._dp, 0._dp, r >= real(BOX,dp)/4._dp + 2._dp .and. r < real(BOX,dp)/2._dp - 6._dp)
                else if( l_sphere )then
                    profile = sqrt(max(0._dp, real(rad_px, dp)**2 - r**2)) / real(rad_px, dp)
                else
                    profile = merge(1._dp, 0._dp, r <= real(rad_px, dp))
                endif
                if( l_bent )then
                    ! circular fringes about a centre 1.5 particle radii away, same spacing
                    phase = kmag * hypot(x - 1.5_dp * real(rad_px, dp), y)
                else
                    phase = kx * x + ky * y
                endif
                select case(kind)
                    case(KIND_LATTICE)
                        rmat(i,j,1) = real(profile * (1._dp + 0.5_dp * cos(phase)))
                    case(KIND_LOWPASS, KIND_NOISE)
                        rmat(i,j,1) = real(profile)
                end select
                rmat(i,j,1) = rmat(i,j,1) + noise_sd * real(gauss())
            end do
        end do
        call img%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call img%set_rmat(rmat, .false.)
    end subroutine make_cavg

    subroutine kill_imgs( imgs )
        type(image), intent(inout) :: imgs(:)
        integer :: i
        do i = 1, size(imgs)
            call imgs(i)%kill
        end do
    end subroutine kill_imgs

    ! a deterministic standard normal deviate (Lehmer stream, Box-Muller)
    real(dp) function gauss()
        real(dp) :: u1, u2
        u1 = next_uniform()
        u2 = next_uniform()
        gauss = sqrt(-2._dp * log(u1)) * cos(2._dp * DPI * u2)
    end function gauss

    function gauss_vec( n ) result( v )
        integer, intent(in) :: n
        real(dp) :: v(n)
        integer  :: i
        do i = 1, n
            v(i) = gauss()
        end do
    end function gauss_vec

    real(dp) function next_uniform()
        rng_state = modulo(rng_state * 48271_8, 2147483647_8)
        next_uniform = (real(rng_state, dp) + 0.5_dp) / 2147483647._dp
    end function next_uniform

end module simple_cavg_quality_lattice_tester
