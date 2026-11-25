!@descr: particle_sieving_2: solve2D on the particles of the first micrographs, the class-average CNN's junk labels, and the class student for the rest in batches with a greedy refine2D round every few batches
! A controller of whole programs (solve2D, learn_class_student, predict_classes, refine2D), each run on a copy of
! the work project in a folder of its own; it owns the work project and writes it once, at the end. Rules:
! doc/policies/sieving_and_rejection/ptcl_sieve_policy.md (its section on particle_sieving_2); decisions:
! doc/implementation_notes/planned/particle_sieving_2.md.
module simple_ptcl_sieve2
use simple_defs,         only: dp, logfhandle
use simple_string,       only: string
use simple_string_utils, only: int2str, int2str_pad
use simple_timer,        only: tic, toc, timer_int_kind
use simple_error,        only: simple_exception
use simple_parameters,   only: parameters
use simple_sp_project,   only: sp_project
use simple_defs_sieve2,  only: SIEVE2_NMICS, SIEVE2_SEED_NPTCLS, SIEVE2_NCLS, SIEVE2_BATCH_NPTCLS, &
    &SIEVE2_ROUND_NBATCHES, SIEVE2_ROUND_OLD_FRAC, SIEVE2_TRAIN_NPTCLS, SIEVE2_JUNK_THRES, SIEVE2_LPSTOP, &
    &SIEVE2_FREEZE_CORR, SIEVE2_FREEZE_NROUNDS
implicit none

public :: ptcl_sieve2, sieve2_settings, run_all_ptcl_sieve2_tests
private
#include "simple_local_flags.inc"

!> the workflow's numbers; nmics, ncls, nparts, mskdiam, lpstop and retrain_rounds are the command line's, the
!! others its constants
type :: sieve2_settings
    integer              :: nmics          = SIEVE2_NMICS
    integer              :: seed_nptcls    = SIEVE2_SEED_NPTCLS
    integer              :: ncls           = SIEVE2_NCLS
    integer              :: batch_nptcls   = SIEVE2_BATCH_NPTCLS
    integer              :: round_nbatches = SIEVE2_ROUND_NBATCHES
    integer              :: train_nptcls   = SIEVE2_TRAIN_NPTCLS
    integer              :: freeze_nrounds = SIEVE2_FREEZE_NROUNDS
    real(dp)             :: junk_thres     = SIEVE2_JUNK_THRES
    real                 :: round_old_frac = SIEVE2_ROUND_OLD_FRAC
    real                 :: freeze_corr    = SIEVE2_FREEZE_CORR    !< above 1 the references never freeze
    integer              :: nparts         = 0                     !< the children's parts; 0 runs them in shared memory
    real                 :: mskdiam        = 0.                    !< the mask diameter (A, mskdiam); 0: the seed's solve2D estimates it
    real                 :: lpstop         = SIEVE2_LPSTOP         !< the alignment's finest low-pass limit (A, lpstop)
    logical              :: final_restore  = .false.               !< restore_cavgs: a frozen last round restores the references
    integer, allocatable :: retrain_rounds(:)                      !< greedy rounds after which the student is trained again (retrain)
end type sieve2_settings

type :: ptcl_sieve2
    private
    type(sieve2_settings)         :: s
    type(sp_project), allocatable :: proj                        !< the work project: the input's rows, owned here
    type(string)                  :: edir                        !< the execution directory
    type(string)                  :: projfile                    !< the work project's file
    type(string)                  :: copy_name                   !< its name, the name of every step's copy
    type(string)                  :: weights                     !< the class-average CNN's weights
    type(string)                  :: cavgs, frcs, sigma2, model  !< the current artifacts, absolute
    type(string)                  :: carried                     !< the carried class sums at the full box, from round 1 on
    real                          :: mskdiam  = 0.
    real                          :: smpd     = 0.               !< the particles' and the class averages' sampling (A)
    integer                       :: box      = 0                !< the full box of the class averages the rounds compare
    integer                       :: nthr     = 1
    integer                       :: nptcls   = 0
    integer                       :: nbatches = 0
    integer                       :: nrounds  = 0
    integer                       :: n_new      = 0              !< particles ingested since the last greedy round
    integer                       :: nnew_srch  = 0              !< new particles the rounds searched, for the student's agreement
    integer                       :: nnew_kept  = 0              !< of them, those the search left in the student's class
    integer                       :: nlab_srch  = 0              !< the new particles searched in rounds with labels
    integer                       :: nlab_kept  = 0              !< of them, those the search left in the student's class
    integer                       :: njunk_srch = 0              !< the new particles searched whose predicted class was junk
    integer                       :: njunk_kept = 0              !< of them, those the search left in the student's class
    integer(timer_int_kind)       :: t_run      = 0
    logical,          allocatable :: sel(:)                      !< selected in the input
    logical,          allocatable :: ingested(:)                 !< classified by the seed or a batch
    logical,          allocatable :: bad(:)                      !< the classes the last labelling found junk
    logical,          allocatable :: pending(:)                  !< batch particles whose sigma2 records are not seeded yet
    real(dp),         allocatable :: pjunk(:)                    !< the last labelling's p_junk, by class
    integer,          allocatable :: cls_pred(:)                 !< new particles the student put in a class, over the rounds
    integer,          allocatable :: cls_kept(:)                 !< of them, those their first round left there
    integer,          allocatable :: cls_into(:)                 !< new particles their first round put in a class
    integer,          allocatable :: ptcl_batch(:)               !< 0 the seed, b > 0 batch b, -1 not selected
    logical                       :: exists   = .false.
  contains
    procedure                  :: new
    procedure                  :: run
    procedure                  :: kill
    procedure, private         :: run_seed
    procedure, private         :: label_classes
    procedure, private         :: train_student
    procedure, private         :: predict_batch
    procedure, private         :: greedy_round
    procedure, private         :: write_result
    procedure, private         :: seed_sigma2
    procedure, private         :: enter_step
    procedure, private         :: write_copy
    procedure, private         :: adopt_artifacts
    procedure, private         :: report_step
    procedure, private         :: refs_change
    procedure, private, nopass :: plan_batches
    procedure, private, nopass :: retrains_after
    procedure, private, nopass :: class_labels
    procedure, private, nopass :: greedy_frac
    procedure, private, nopass :: round_lowpass
    procedure, private, nopass :: round_agreement
    procedure, private, nopass :: class_agreement
    procedure, private, nopass :: kept_order
    procedure, private, nopass :: class_pops
    procedure, private, nopass :: freeze_count
    procedure, private, nopass :: newonly_frac
    procedure, private, nopass :: final_states
    procedure, private, nopass :: check_classes_left
    procedure, private, nopass :: check_none_pending
    procedure, private, nopass :: both_halves
end type ptcl_sieve2

interface
    !> white-box tests (submodule simple_ptcl_sieve2_tester, built only with BUILD_TESTS=ON)
    module subroutine run_all_ptcl_sieve2_tests()
    end subroutine run_all_ptcl_sieve2_tests
end interface

contains

    !> the work project read from params%projfile in the execution directory (the working directory), its selected
    !! particles cut into the seed's stacks and the batches
    subroutine new( self, params, settings )
        use simple_syslib, only: simple_getcwd, simple_abspath
        use simple_fileio, only: basename
        class(ptcl_sieve2),    intent(inout) :: self
        class(parameters),     intent(in)    :: params
        type(sieve2_settings), intent(in)    :: settings
        integer,          allocatable :: stk_nsel(:), stk_batch(:)
        character(len=:), allocatable :: rounds_str
        integer                       :: istk, iptcl, nstks, nmics, nrounds_plan, i
        call self%kill
        if( settings%nmics < 1 .or. settings%ncls < 2 .or. settings%batch_nptcls < 1 .or. settings%round_nbatches < 1 .or. &
            &settings%train_nptcls < 1 .or. settings%seed_nptcls < 1 .or. settings%mskdiam < 0. .or. settings%lpstop <= 0. .or. &
            &settings%round_old_frac < 0. .or. settings%round_old_frac > 1. .or. settings%freeze_nrounds < 1 )then
            THROW_HARD('invalid settings; ptcl_sieve2%new')
        endif
        if( allocated(settings%retrain_rounds) )then
            if( any(settings%retrain_rounds < 1) ) THROW_HARD('retrain rounds start at 1; ptcl_sieve2%new')
        endif
        self%s       = settings
        self%nthr    = params%nthr
        self%mskdiam = settings%mskdiam
        self%box     = params%box
        self%smpd    = params%smpd
        call simple_getcwd(self%edir)
        self%projfile  = simple_abspath(params%projfile)
        self%copy_name = basename(self%projfile)
        self%weights   = simple_abspath(params%infile)
        allocate(self%proj)
        call self%proj%read(self%projfile)
        self%nptcls = self%proj%os_ptcl2D%get_noris()
        nstks       = self%proj%os_stk%get_noris()
        nmics       = self%proj%os_mic%get_noris()
        if( self%nptcls < 1 ) THROW_HARD('no particles in the project; ptcl_sieve2%new')
        if( nstks < 1 )       THROW_HARD('no stacks in the project; ptcl_sieve2%new')
        if( nmics > 0 .and. nmics /= nstks ) THROW_HARD('micrographs and stacks differ in number, use prune_project')
        if( self%proj%os_ptcl3D%get_noris() /= self%nptcls ) THROW_HARD('ptcl2D and ptcl3D differ in length; ptcl_sieve2%new')
        allocate(self%sel(self%nptcls), self%ingested(self%nptcls), self%pending(self%nptcls), source=.false.)
        allocate(self%ptcl_batch(self%nptcls), source=-1)
        allocate(self%cls_pred(settings%ncls), self%cls_kept(settings%ncls), self%cls_into(settings%ncls), source=0)
        allocate(stk_nsel(nstks), source=0)
        do istk = 1,nstks
            if( self%proj%os_stk%get_state(istk) == 0 ) cycle
            do iptcl = self%proj%os_stk%get_fromp(istk), self%proj%os_stk%get_top(istk)
                if( self%proj%os_ptcl2D%get_state(iptcl) == 0 ) cycle
                self%sel(iptcl) = .true.
                stk_nsel(istk)  = stk_nsel(istk) + 1
            enddo
        enddo
        if( .not. any(self%sel) ) THROW_HARD('no selected particle; ptcl_sieve2%new')
        call self%plan_batches(stk_nsel, settings%nmics, settings%seed_nptcls, settings%batch_nptcls, stk_batch, &
            &self%nbatches)
        do istk = 1,nstks
            do iptcl = self%proj%os_stk%get_fromp(istk), self%proj%os_stk%get_top(istk)
                if( self%sel(iptcl) ) self%ptcl_batch(iptcl) = stk_batch(istk)
            enddo
        enddo
        write(logfhandle,'(A,I0,A,I0,A)') '>>> PARTICLE_SIEVING_2 SEED: ', count(stk_batch == 0), ' micrographs, ', &
            &count(self%ptcl_batch == 0), ' particles'
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> PARTICLE_SIEVING_2 BATCHES: ', self%nbatches, ' holding ', &
            &count(self%ptcl_batch > 0), ' particles, a greedy round every ', settings%round_nbatches, ' batches'
        ! a retraining serves only the batches after its round
        nrounds_plan = (self%nbatches + settings%round_nbatches - 1) / settings%round_nbatches
        rounds_str   = ''
        if( allocated(settings%retrain_rounds) )then
            do i = 1,size(settings%retrain_rounds)
                if( settings%retrain_rounds(i) < nrounds_plan )then
                    rounds_str = rounds_str//' '//int2str(settings%retrain_rounds(i))
                else
                    THROW_WARN('retrain: no batch follows greedy round '//int2str(settings%retrain_rounds(i))//'; ignored')
                endif
            enddo
        endif
        if( len(rounds_str) > 0 ) write(logfhandle,'(A,A)') '>>> PARTICLE_SIEVING_2 STUDENT: TRAINED AFTER THE SEED AND &
            &AFTER GREEDY ROUNDS', rounds_str
        if( self%mskdiam > 0. )then
            write(logfhandle,'(A,F8.1,A)') '>>> PARTICLE_SIEVING_2 MASK DIAMETER: ', self%mskdiam, ' A'
        else
            write(logfhandle,'(A)') '>>> PARTICLE_SIEVING_2 MASK DIAMETER: the seed''s solve2D estimate'
        endif
        write(logfhandle,'(A,F6.1,A)') '>>> PARTICLE_SIEVING_2 ALIGNMENT LOW-PASS LIMIT: ', self%s%lpstop, ' A (lpstop)'
        deallocate(stk_nsel, stk_batch)
        self%exists = .true.
    end subroutine new

    !> the seed and the student's training, then the batches, a greedy round every round_nbatches batches and after the
    !! last, the student trained again after the rounds in retrain_rounds while batches remain; once freeze_nrounds
    !! rounds running left the class averages at a weighted correlation of freeze_corr or more with the previous ones,
    !! the class-average CNN labels them, once, and the rounds search the new particles only and keep the references
    !! (the last restores them with final_restore); without a freeze the last round's averages are labelled; the result
    subroutine run( self )
        class(ptcl_sieve2), intent(inout) :: self
        real    :: refcorr
        integer :: ib, ncons
        logical :: l_frozen
        if( .not. self%exists ) THROW_HARD('not constructed; ptcl_sieve2%run')
        self%t_run = tic()
        call self%run_seed
        ! no batch, no round: the seed's averages are the final ones
        if( self%nbatches < 1 ) call self%label_classes('seed')
        if( self%nbatches > 0 ) call self%train_student
        ncons = 0
        do ib = 1,self%nbatches
            call self%predict_batch(ib)
            if( mod(ib, self%s%round_nbatches) == 0 .or. ib == self%nbatches )then
                l_frozen = ncons >= self%s%freeze_nrounds
                if( l_frozen .and. ncons == self%s%freeze_nrounds ) write(logfhandle,'(A,I0,A,I0,A,F5.3,A)') &
                    &'>>> PARTICLE_SIEVING_2: REFERENCES FROZEN FROM ROUND ', self%nrounds + 1, ' (', ncons, &
                    &' ROUNDS AT A CORRELATION OF ', self%s%freeze_corr, ' OR MORE)'
                call self%greedy_round(l_frozen, ib == self%nbatches, refcorr)
                ncons = self%freeze_count(ncons, refcorr, self%s%freeze_corr)
                ! the class-average CNN's one selection: on the references once they froze, else on the last averages
                if( .not. allocated(self%bad) )then
                    if( ncons >= self%s%freeze_nrounds .or. ib == self%nbatches ) &
                        &call self%label_classes('round_'//int2str_pad(self%nrounds, 2))
                endif
                if( ib < self%nbatches )then
                    if( self%retrains_after(self%s%retrain_rounds, self%nrounds) ) call self%train_student
                endif
            endif
        enddo
        call self%write_result
    end subroutine run

    subroutine kill( self )
        class(ptcl_sieve2), intent(inout) :: self
        if( allocated(self%proj) )then
            call self%proj%kill
            deallocate(self%proj)
        endif
        if( allocated(self%sel)        ) deallocate(self%sel)
        if( allocated(self%ingested)   ) deallocate(self%ingested)
        if( allocated(self%ptcl_batch) ) deallocate(self%ptcl_batch)
        if( allocated(self%bad)        ) deallocate(self%bad)
        if( allocated(self%pending)    ) deallocate(self%pending)
        if( allocated(self%pjunk)      ) deallocate(self%pjunk)
        if( allocated(self%cls_pred)   ) deallocate(self%cls_pred, self%cls_kept, self%cls_into)
        call self%edir%kill
        call self%projfile%kill
        call self%copy_name%kill
        call self%weights%kill
        call self%cavgs%kill
        call self%frcs%kill
        call self%sigma2%kill
        call self%carried%kill
        call self%model%kill
        self%s        = sieve2_settings()
        self%mskdiam  = 0.
        self%smpd     = 0.
        self%box      = 0
        self%nthr     = 1
        self%nptcls   = 0
        self%nbatches = 0
        self%nrounds  = 0
        self%n_new    = 0
        self%nnew_srch = 0
        self%nnew_kept = 0
        self%nlab_srch = 0
        self%nlab_kept = 0
        self%njunk_srch = 0
        self%njunk_kept = 0
        self%t_run    = 0
        self%exists   = .false.
    end subroutine kill

    ! ---- the steps -------------------------------------------------------------------------------------------------

    !> the halves of every selected particle, then solve2D on the seed's particles, its stages marching to lpstop;
    !! its particle fields, classes and artifacts become the work project's, and each seed particle counts one update,
    !! so that the rounds count greedy searches
    subroutine run_seed( self )
        use simple_cmdline,            only: cmdline
        use simple_commanders_solve2D, only: commander_solve2D
        use simple_defs_sieve2,        only: SIEVE2_SEED_DIR
        use simple_eo_partition,       only: ensure_field_halves
        use simple_syslib,             only: simple_chdir
        class(ptcl_sieve2), intent(inout) :: self
        type(commander_solve2D) :: xsolve2D
        type(cmdline)           :: cl
        type(sp_project)        :: res
        integer(timer_int_kind) :: t0
        integer, allocatable    :: pinds(:)
        integer                 :: i
        t0 = tic()
        ! the halves made here, by micrograph group over every selected particle: the seed's solve2D would label only
        ! its own, and the stacks it does not hold would all fall into the lighter half
        call ensure_field_halves(self%proj, 'ptcl2D', self%projfile)
        write(logfhandle,'(A)') '>>> PARTICLE_SIEVING_2: SOLVE2D ON THE SEED'
        call self%enter_step(SIEVE2_SEED_DIR)
        call self%write_copy(self%ptcl_batch == 0)
        call cl%set('prg',      'solve2D')
        call cl%set('projfile', self%copy_name%to_char())
        call cl%set('ncls',     self%s%ncls)
        call cl%set('mskdiam',  self%mskdiam)
        call cl%set('lpstop',   self%s%lpstop)
        call cl%set('nthr',     self%nthr)
        call cl%set('mkdir',    'no')
        if( self%s%nparts > 0 ) call cl%set('nparts', self%s%nparts)
        call xsolve2D%execute(cl)
        call cl%kill
        call res%read(self%copy_name)
        self%proj%os_stk    = res%os_stk
        self%proj%os_ptcl2D = res%os_ptcl2D
        self%proj%os_ptcl3D = res%os_ptcl3D
        self%proj%os_cls2D  = res%os_cls2D
        self%proj%os_cls3D  = res%os_cls3D
        self%proj%os_out    = res%os_out
        call self%adopt_artifacts(res)
        if( self%mskdiam <= 0. )then
            ! the automatic mask: the one solve2D registered with its final class averages, for every later step
            call res%get_mskdiam('cavg', self%mskdiam)
            if( self%mskdiam <= 0. ) THROW_HARD('the seed''s solve2D registered no mask diameter; ptcl_sieve2')
            write(logfhandle,'(A,F8.1,A)') '>>> PARTICLE_SIEVING_2 MASK DIAMETER: ', self%mskdiam, ' A (the seed''s solve2D)'
        endif
        call res%kill
        call simple_chdir(self%edir)
        self%ingested = self%ptcl_batch == 0
        pinds = pack([(i, i = 1,self%nptcls)], self%ingested)
        call self%proj%os_ptcl2D%set_updatecnt(1, pinds)
        call self%report_step('seed', toc(t0))
    end subroutine run_seed

    !> the CNN's probabilities of the current class averages; a populated class is junk when p_junk exceeds the
    !! threshold. The table labels_<tag>.txt; stops when no populated class is left
    subroutine label_classes( self, tag )
        use simple_cavg_cnn_predict, only: cavg_cnn_class_probs
        use simple_defs_sieve2,      only: SIEVE2_LABELS_FBODY
        use simple_fileio,           only: fopen, fclose, fileiochk
        class(ptcl_sieve2), intent(inout) :: self
        character(len=*),   intent(in)    :: tag
        real(dp), allocatable   :: pclass(:,:)
        logical,  allocatable   :: empty(:)
        integer,  allocatable   :: pops(:)
        type(string)            :: stk, fname
        integer(timer_int_kind) :: t0
        real                    :: smpd, smpd_net
        integer                 :: ncls, icls, iptcl, funit, ios
        t0 = tic()
        call self%proj%get_cavgs_stk(stk, ncls, smpd)
        call cavg_cnn_class_probs(self%weights, self%cavgs, ncls, smpd, pclass, empty, smpd_net)
        allocate(pops(ncls), source=0)
        do iptcl = 1,self%nptcls
            if( .not. self%ingested(iptcl) ) cycle
            icls = self%proj%os_ptcl2D%get_class(iptcl)
            if( icls < 1 .or. icls > ncls ) cycle
            pops(icls) = pops(icls) + 1
        enddo
        call self%class_labels(pclass(1,:), empty, pops, self%s%junk_thres, self%bad)
        self%pjunk = pclass(1,:)
        fname = SIEVE2_LABELS_FBODY//tag//'.txt'
        call fopen(funit, fname, status='replace', action='write', iostat=ios)
        call fileiochk('ptcl_sieve2; open '//fname%to_char(), ios)
        write(funit,'(A)') '# class pop p_junk p_ss p_likely empty bad'
        do icls = 1,ncls
            write(funit,'(I0,1X,I0,3(1X,ES15.8),2(1X,I0))') icls, pops(icls), pclass(:,icls), merge(1, 0, empty(icls)), &
                &merge(1, 0, self%bad(icls))
        enddo
        call fclose(funit)
        write(logfhandle,'(A,A,A,I0,A,I0,A,F6.3,A,F6.3,A)') '>>> PARTICLE_SIEVING_2 LABELS ', tag, ': ', count(self%bad), &
            &' of ', count(pops > 0), ' populated classes junk (p_junk > ', self%s%junk_thres, '; network at ', smpd_net, ' A)'
        call self%check_classes_left(self%bad, pops)
        call self%report_step('labels_'//tag, toc(t0))
    end subroutine label_classes

    !> learn_class_student on a random sample of the ingested particles (after the seed, the seed's), those of junk
    !! classes included, against the current class averages: the model the next batches are predicted with
    subroutine train_student( self )
        use simple_cmdline,                  only: cmdline
        use simple_commanders_class_student, only: commander_learn_class_student
        use simple_defs_hhcnn,               only: STUDENT_LP_VALIDATE, CLASS_STUDENT_MODEL_FNAME
        use simple_defs_sieve2,              only: SIEVE2_LEARN_DIR
        use simple_ran_tabu,                 only: ran_tabu
        use simple_syslib,                   only: simple_chdir, simple_abspath
        class(ptcl_sieve2), intent(inout) :: self
        type(commander_learn_class_student) :: xlearn
        type(cmdline)                       :: cl
        type(ran_tabu)                      :: rt
        integer, allocatable                :: pinds(:)
        logical, allocatable                :: train(:)
        integer(timer_int_kind)             :: t0
        integer                             :: i, n, ntrain
        t0 = tic()
        pinds  = pack([(i, i = 1,self%nptcls)], self%ingested)
        n      = size(pinds)
        ntrain = min(self%s%train_nptcls, n)
        if( ntrain < n )then
            rt = ran_tabu(n)
            call rt%shuffle(pinds)
            call rt%kill
        endif
        allocate(train(self%nptcls), source=.false.)
        train(pinds(1:ntrain)) = .true.
        write(logfhandle,'(A,I0,A,I0,A)') '>>> PARTICLE_SIEVING_2: TRAINING THE CLASS STUDENT ON ', ntrain, ' OF ', n, &
            &' INGESTED PARTICLES'
        call self%enter_step(SIEVE2_LEARN_DIR//int2str_pad(self%nrounds, 2))
        call self%write_copy(train)
        call cl%set('prg',      'learn_class_student')
        call cl%set('projfile', self%copy_name%to_char())
        call cl%set('refs',     self%cavgs%to_char())
        call cl%set('lp',       STUDENT_LP_VALIDATE)
        call cl%set('mskdiam',  self%mskdiam)
        call cl%set('nthr',     self%nthr)
        call cl%set('fname',    CLASS_STUDENT_MODEL_FNAME)
        call cl%set('mkdir',    'no')
        call xlearn%execute(cl)
        call cl%kill
        self%model = simple_abspath(string(CLASS_STUDENT_MODEL_FNAME))
        call simple_chdir(self%edir)
        call self%report_step(SIEVE2_LEARN_DIR//int2str_pad(self%nrounds, 2), toc(t0))
    end subroutine train_student

    !> predict_classes on batch ib against the current class averages: the student's candidates and the alignment
    !! over them give each particle its class, angle and shift (nothing else is searched); then the sigma2 records of
    !! the pending particles are seeded
    subroutine predict_batch( self, ib )
        use simple_cmdline,                  only: cmdline
        use simple_commanders_class_student, only: commander_predict_classes
        use simple_defs_hhcnn,               only: CLASS_STUDENT_MODEL_FNAME
        use simple_defs_sieve2,              only: SIEVE2_BATCH_DIR
        use simple_syslib,                   only: simple_chdir
        class(ptcl_sieve2), intent(inout) :: self
        integer,            intent(in)    :: ib
        type(commander_predict_classes) :: xpredict
        type(cmdline)                   :: cl
        type(sp_project)                :: res
        logical, allocatable            :: mask(:)
        integer(timer_int_kind)         :: t0
        integer                         :: iptcl
        t0   = tic()
        mask = self%ptcl_batch == ib
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> PARTICLE_SIEVING_2: BATCH ', ib, ' OF ', self%nbatches, ', ', &
            &count(mask), ' PARTICLES'
        call self%enter_step(SIEVE2_BATCH_DIR//int2str_pad(ib, 3))
        call self%write_copy(mask)
        call cl%set('prg',      'predict_classes')
        call cl%set('projfile', self%copy_name%to_char())
        call cl%set('infile',   self%model%to_char())
        call cl%set('refs',     self%cavgs%to_char())
        call cl%set('mskdiam',  self%mskdiam)
        call cl%set('nthr',     self%nthr)
        call cl%set('fname',    CLASS_STUDENT_MODEL_FNAME)
        call cl%set('mkdir',    'no')
        if( self%s%nparts > 0 ) call cl%set('nparts', self%s%nparts)
        call xpredict%execute(cl)
        call cl%kill
        call res%read_segment('ptcl2D', self%copy_name)
        do iptcl = 1,self%nptcls
            if( mask(iptcl) ) call self%proj%os_ptcl2D%transfer_ori(iptcl, res%os_ptcl2D, iptcl)
        enddo
        call res%kill
        self%pending = self%pending .or. mask
        call self%seed_sigma2
        call simple_chdir(self%edir)
        self%ingested = self%ingested .or. mask
        self%n_new    = self%n_new + count(mask)
        call self%report_step('batch_'//int2str_pad(ib, 3), toc(t0))
    end subroutine predict_batch

    !> one greedy refine2D iteration at the full box on the ingested particles, junk classes included, at the class
    !! FRCs' limit, never finer than lpstop: the particles since the last round and a share of the earlier ones, lowest
    !! update count first, blended with the carried class sums (the first round, without any, a full update); frozen,
    !! the new particles only, the references kept but in the last round with final_restore. refcorr: how much the
    !! class averages moved
    subroutine greedy_round( self, l_frozen, l_last, refcorr )
        use simple_cmdline,             only: cmdline
        use simple_class_frcs,          only: class_frcs
        use simple_commanders_refine2D, only: commander_refine2D
        use simple_defs_fname,          only: FRCS_FILE, CAVG_STATE_FILE
        use simple_defs_sieve2,         only: SIEVE2_ROUND_DIR, SIEVE2_ROUND_ITER0
        use simple_fileio,              only: simple_copy_file
        use simple_syslib,              only: simple_chdir, simple_abspath, file_exists
        class(ptcl_sieve2), intent(inout) :: self
        logical,            intent(in)    :: l_frozen, l_last
        real,               intent(out)   :: refcorr
        type(commander_refine2D)      :: xrefine2D
        type(class_frcs)              :: clsfrcs
        type(cmdline)                 :: cl
        type(sp_project)              :: res
        type(string)                  :: cavgs_prev
        integer(timer_int_kind)       :: t0
        real                          :: frac, lp
        real                          :: split(2)
        integer                       :: n_active, iptcl, nnew, nnew_kept, nold, nold_kept, njunk, njunk_kept
        logical                       :: l_carried, l_labels, l_keep
        integer,          allocatable :: cls0(:), cnt0(:), cls1(:), cnt1(:)
        logical,          allocatable :: nolabels(:)
        character(len=:), allocatable :: msg
        t0 = tic()
        self%nrounds = self%nrounds + 1
        n_active  = count(self%ingested)
        l_carried = self%carried%is_allocated()
        frac = 1.
        if( l_carried ) frac = self%greedy_frac(self%n_new, n_active, self%s%round_old_frac)
        if( l_carried .and. l_frozen ) frac = self%newonly_frac(self%n_new, n_active)
        ! the limit refine2D would take from the class FRCs, held at lpstop
        call clsfrcs%read(self%frcs)
        lp = self%round_lowpass(clsfrcs%estimate_lp_for_align(), self%s%lpstop)
        call clsfrcs%kill
        write(logfhandle,'(A,I0,A,I0,A,I0,A,F6.3,A,F6.2,A)') '>>> PARTICLE_SIEVING_2: GREEDY ROUND ', self%nrounds, ' ON ', &
            &n_active, ' PARTICLES, ', self%n_new, ' NEW, UPDATE FRACTION ', frac, ', LOW-PASS ', lp, ' A'
        l_keep = l_frozen .and. .not. (l_last .and. self%s%final_restore)
        if( l_frozen )then
            if( l_keep )then
                write(logfhandle,'(A)') '>>> PARTICLE_SIEVING_2: FROZEN: THE NEW PARTICLES ONLY, THE REFERENCES KEPT'
            else
                write(logfhandle,'(A)') '>>> PARTICLE_SIEVING_2: FROZEN: THE NEW PARTICLES ONLY, THE REFERENCES RESTORED (LAST ROUND)'
            endif
        endif
        allocate(cls0(self%nptcls), cnt0(self%nptcls), cls1(self%nptcls), cnt1(self%nptcls))
        do iptcl = 1,self%nptcls
            cls0(iptcl) = self%proj%os_ptcl2D%get_class(iptcl)
            cnt0(iptcl) = self%proj%os_ptcl2D%get_updatecnt(iptcl)
        enddo
        cavgs_prev = self%cavgs
        call self%enter_step(SIEVE2_ROUND_DIR//int2str_pad(self%nrounds, 2))
        call self%write_copy(self%ingested)
        call simple_copy_file(self%frcs, string(FRCS_FILE))
        if( l_carried ) call simple_copy_file(self%carried, string(CAVG_STATE_FILE))
        call cl%set('prg',       'refine2D')
        call cl%set('projfile',  self%copy_name%to_char())
        call cl%set('refs',      self%cavgs%to_char())
        call cl%set('frcs',      FRCS_FILE)
        call cl%set('ncls',      self%s%ncls)
        call cl%set('mskdiam',   self%mskdiam)
        call cl%set('nthr',      self%nthr)
        call cl%set('mkdir',     'no')
        call cl%set('refine',    'greedy')
        call cl%set('autoscale', 'no')
        call cl%set('center',    'no')
        call cl%set('startit',   SIEVE2_ROUND_ITER0 + self%nrounds - 1)
        call cl%set('minits',    1)
        call cl%set('maxits',    1)
        call cl%set('lp',        lp)
        ! solve2D's stage keys: with lp set, a standalone refine2D would otherwise align against merged references
        call cl%set('stage',          1)
        call cl%set('gold_std_stage', 1)
        if( frac < 1. ) call cl%set('update_frac', frac)
        if( l_keep ) call cl%set('restore_cavgs', 'no')
        if( self%s%nparts > 0 ) call cl%set('nparts', self%s%nparts)
        call xrefine2D%execute(cl)
        call cl%kill
        call res%read(self%copy_name)
        self%proj%os_ptcl2D = res%os_ptcl2D
        self%proj%os_cls2D  = res%os_cls2D
        self%proj%os_cls3D  = res%os_cls3D
        self%proj%os_out    = res%os_out
        call self%adopt_artifacts(res)
        call res%kill
        if( .not. file_exists(string(CAVG_STATE_FILE)) ) THROW_HARD('the greedy round left no carried class sums; ptcl_sieve2')
        self%carried = simple_abspath(string(CAVG_STATE_FILE))
        call simple_chdir(self%edir)
        self%n_new = 0
        ! the student's agreement: the new particles' classes against the search's, split by the label their predicted
        ! class had when the batch was predicted (the labels this round's own replace next), beside the search's own
        ! agreement on the earlier particles it searched again
        do iptcl = 1,self%nptcls
            cls1(iptcl) = self%proj%os_ptcl2D%get_class(iptcl)
            cnt1(iptcl) = self%proj%os_ptcl2D%get_updatecnt(iptcl)
        enddo
        ! the split by the predicted class's label needs the labels, which come with the freeze
        l_labels = allocated(self%bad)
        if( l_labels )then
            call self%round_agreement(self%ingested, cls0, cls1, cnt0, cnt1, self%bad, nnew, nnew_kept, njunk, &
                &njunk_kept, nold, nold_kept)
        else
            allocate(nolabels(self%s%ncls), source=.false.)
            call self%round_agreement(self%ingested, cls0, cls1, cnt0, cnt1, nolabels, nnew, nnew_kept, njunk, &
                &njunk_kept, nold, nold_kept)
        endif
        self%nnew_srch = self%nnew_srch + nnew
        self%nnew_kept = self%nnew_kept + nnew_kept
        split = -1.
        if( l_labels )then
            self%nlab_srch  = self%nlab_srch  + nnew
            self%nlab_kept  = self%nlab_kept  + nnew_kept
            self%njunk_srch = self%njunk_srch + njunk
            self%njunk_kept = self%njunk_kept + njunk_kept
            split = [pct(nnew_kept - njunk_kept, nnew - njunk), pct(njunk_kept, njunk)]
        endif
        call self%class_agreement(self%ingested, cls0, cls1, cnt0, cnt1, self%cls_pred, self%cls_kept, self%cls_into)
        msg = '>>> PARTICLE_SIEVING_2 STUDENT: ROUND '//int2str(self%nrounds)//' KEPT '//int2str(nnew_kept)//' OF '// &
            &int2str(nnew)//' NEW PARTICLES IN THE PREDICTED CLASS ('//pct_str(pct(nnew_kept, nnew))//' %'
        if( l_labels ) msg = msg//'; NOT JUNK '//pct_str(split(1))//' %, JUNK '//pct_str(split(2))//' %'
        msg = msg//')'
        if( nold > 0 ) msg = msg//'; EARLIER PARTICLES SEARCHED AGAIN: '//int2str(nold_kept)//' OF '//int2str(nold)// &
            &' KEPT ('//pct_str(pct(nold_kept, nold))//' %)'
        write(logfhandle,'(A)') msg
        ! how much the references moved: each class's average against its previous one, by its new population
        refcorr = self%refs_change(cavgs_prev, self%cavgs, self%ingested .and. cnt0 > 0, cls0, self%ingested, cls1, lp)
        write(logfhandle,'(A,I0,A,F7.4,A,F6.2,A)') '>>> PARTICLE_SIEVING_2: ROUND ', self%nrounds, &
            &' CLASS AVERAGES AGAINST THE PREVIOUS ONES: WEIGHTED CORRELATION ', refcorr, ' (TO ', lp, ' A)'
        call self%report_step('round_'//int2str_pad(self%nrounds, 2), toc(t0), [pct(nnew_kept, nnew), split(1), &
            &split(2), pct(nold_kept, nold)], refcorr)
    end subroutine greedy_round

    !> the states from the last labels (a particle is accepted when its class is not junk), the class states, the 2D
    !! shifts in ptcl3D unless it holds a 3D alignment; no sigma2 store, as particle_sieving's result has none (a later
    !! program rebuilds it); one write
    subroutine write_result( self )
        use simple_defs_fname,  only: METADATA_EXT
        use simple_defs_sieve2, only: SIEVE2_REPORT
        use simple_fileio,      only: add2fbody
        use simple_syslib,      only: simple_rename, del_file, file_exists
        class(ptcl_sieve2), intent(inout) :: self
        integer, allocatable :: states(:), classes(:)
        type(string)         :: tmp
        integer              :: iptcl, icls, nunclassified
        call self%check_none_pending(self%sel, self%ingested)
        allocate(classes(self%nptcls))
        do iptcl = 1,self%nptcls
            classes(iptcl) = self%proj%os_ptcl2D%get_class(iptcl)
        enddo
        call self%final_states(self%sel, self%ingested, classes, self%bad, states)
        nunclassified = count(self%ingested .and. (classes < 1 .or. classes > size(self%bad)))
        if( nunclassified > 0 ) write(logfhandle,'(A,I0,A)') '>>> PARTICLE_SIEVING_2: ', nunclassified, &
            &' ingested particles hold no class and are rejected'
        if( any(self%pending) ) write(logfhandle,'(A,I0,A)') '>>> PARTICLE_SIEVING_2 SIGMA2: ', count(self%pending), &
            &' particles of one half were never seeded (no later batch)'
        do iptcl = 1,self%nptcls
            call self%proj%os_ptcl2D%set_state(iptcl, states(iptcl))
            call self%proj%os_ptcl3D%set_state(iptcl, states(iptcl))
        enddo
        if( self%proj%os_cls2D%get_noris() == size(self%bad) )then
            do icls = 1,size(self%bad)
                if( self%bad(icls) ) call self%proj%os_cls2D%set_state(icls, 0)
            enddo
        endif
        if( self%proj%os_cls3D%get_noris() == size(self%bad) )then
            do icls = 1,size(self%bad)
                if( self%bad(icls) ) call self%proj%os_cls3D%set_state(icls, 0)
            enddo
        endif
        if( .not. self%proj%has_valid_ptcl3D_alignment() ) call self%proj%os_ptcl3D%transfer_2Dshifts(self%proj%os_ptcl2D)
        call self%proj%update_projinfo(self%projfile)
        if( self%proj%projinfo%isthere(1, 'sigma2_state') ) call self%proj%projinfo%delete_entry('sigma2_state')
        tmp = add2fbody(self%projfile, METADATA_EXT, '_tmp_sieve2')
        if( file_exists(tmp) ) call del_file(tmp)
        call self%proj%write(tmp)
        call simple_rename(tmp, self%projfile)
        call self%report_step('result', toc(self%t_run))
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> PARTICLE_SIEVING_2: ', count(states == 1), ' OF ', count(self%sel), &
            &' PARTICLES ACCEPTED, ', count(self%bad), ' CLASSES JUNK'
        write(logfhandle,'(A,I0,A,I0,A)') '>>> PARTICLE_SIEVING_2: ', self%nbatches, ' BATCHES, ', self%nrounds, ' GREEDY ROUNDS'
        if( self%nnew_srch > 0 ) write(logfhandle,'(A,I0,A,I0,A,F5.1,A,F5.1,A,F5.1,A)') '>>> PARTICLE_SIEVING_2 STUDENT: ', &
            &self%nnew_kept, ' OF ', self%nnew_srch, ' NEW PARTICLES KEPT THE PREDICTED CLASS IN THEIR FIRST ROUND (', &
            &pct(self%nnew_kept, self%nnew_srch), ' %; IN THE ROUNDS WITH LABELS: NOT JUNK ', pct(self%nlab_kept - &
            &self%njunk_kept, self%nlab_srch - self%njunk_srch), ' %, JUNK ', pct(self%njunk_kept, self%njunk_srch), ' %)'
        write(logfhandle,'(A)') '>>> PARTICLE_SIEVING_2 REPORT: '//SIEVE2_REPORT
        if( self%nrounds > 0 ) call write_student_table

      contains

        !> the student's agreement by class, ranked by the share kept (classes without predictions last): rank, class,
        !! new particles predicted there, kept, kept %, new particles the rounds put there, the final population, the
        !! last p_junk and junk label
        subroutine write_student_table
            use simple_defs_sieve2, only: SIEVE2_STUDENT_TABLE
            use simple_fileio,      only: fopen, fclose, fileiochk
            integer, allocatable :: order(:), pops(:)
            character(len=8)     :: kept_str
            integer :: funit, ios, k, c
            call self%class_pops(self%ingested, classes, size(self%cls_pred), pops)
            order = self%kept_order(self%cls_pred, self%cls_kept)
            call fopen(funit, string(SIEVE2_STUDENT_TABLE), status='replace', action='write', iostat=ios)
            call fileiochk('ptcl_sieve2; open '//SIEVE2_STUDENT_TABLE, ios)
            write(funit,'(A)') '# rank class predicted kept kept_pct placed population p_junk junk'
            do k = 1,size(order)
                c = order(k)
                kept_str = '-'
                if( self%cls_pred(c) > 0 ) write(kept_str,'(F6.1)') pct(self%cls_kept(c), self%cls_pred(c))
                write(funit,'(I0,3(1X,I0),1X,A,2(1X,I0),1X,F6.4,1X,I0)') k, c, self%cls_pred(c), self%cls_kept(c), &
                    &trim(adjustl(kept_str)), self%cls_into(c), pops(c), self%pjunk(c), merge(1, 0, self%bad(c))
            enddo
            call fclose(funit)
            write(logfhandle,'(A)') '>>> PARTICLE_SIEVING_2 STUDENT BY CLASS: '//SIEVE2_STUDENT_TABLE
        end subroutine write_student_table

    end subroutine write_result

    !> the pending batch particles' sigma2 records started from the store's global curve of their half, the others
    !! kept, published over the folder's store: a round updates only its band, so an unseeded record stays invalid.
    !! A pending set of one half waits for the next batch: the extension takes both halves' means over the seeded rows
    subroutine seed_sigma2( self )
        use, intrinsic :: iso_fortran_env, only: int64
        use simple_defs,              only: STDLEN
        use simple_fileio,            only: basename
        use simple_syslib,            only: simple_abspath
        use simple_sigma2_state,      only: sigma2_state_extend, sigma2_state_candidate_path, sigma2_state_next_generation
        use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_init_header, sigma2_state_read_header, &
            &sigma2_state_publish
        class(ptcl_sieve2), intent(inout) :: self
        type(sigma2_state_header) :: frozen, header
        type(string)              :: store, candidate
        integer, allocatable      :: eo(:)
        character(len=STDLEN)     :: msg
        integer(int64)            :: generation
        integer                   :: iptcl, status
        if( .not. any(self%pending) ) return
        allocate(eo(self%nptcls))
        do iptcl = 1,self%nptcls
            eo(iptcl) = self%proj%os_ptcl2D%get_eo(iptcl)
        enddo
        if( .not. self%both_halves(self%pending, eo) )then
            write(logfhandle,'(A,I0,A)') '>>> PARTICLE_SIEVING_2 SIGMA2: ', count(self%pending), &
                &' particles of one half wait for the next batch''s seeding'
            return
        endif
        store = simple_abspath(basename(self%sigma2))
        call sigma2_state_read_header(store%to_char(), frozen, status, msg)
        if( status /= 0 ) THROW_HARD(trim(msg))
        if( int(frozen%nptcls) /= self%nptcls ) THROW_HARD('the sigma2 store holds another particle count; ptcl_sieve2')
        call sigma2_state_next_generation(store%to_char(), generation, status, msg)
        if( status /= 0 ) THROW_HARD(trim(msg))
        call sigma2_state_init_header(header, int(frozen%kfrom), int(frozen%kto), self%nptcls, int(frozen%box), &
            &real(frozen%smpd), generation, frozen%layout_digest, frozen%provenance)
        candidate = sigma2_state_candidate_path(store%to_char(), generation)
        call sigma2_state_extend(store%to_char(), candidate%to_char(), header, .not. self%pending, self%pending, eo, &
            &status, msg)
        if( status == 0 ) call sigma2_state_publish(candidate%to_char(), store%to_char(), status, msg)
        if( status /= 0 ) THROW_HARD('ptcl_sieve2 sigma2 seeding: '//trim(msg))
        write(logfhandle,'(A,I0,A)') '>>> PARTICLE_SIEVING_2 SIGMA2: ', count(self%pending), &
            &' records seeded from the global curves of their halves'
        self%sigma2  = store
        self%pending = .false.
    end subroutine seed_sigma2

    ! ---- the steps' frame ------------------------------------------------------------------------------------------

    !> the folder of a step, made in the execution directory and entered; the step returns to the execution directory
    subroutine enter_step( self, dir )
        use simple_syslib, only: simple_chdir, simple_mkdir
        class(ptcl_sieve2), intent(inout) :: self
        character(len=*),   intent(in)    :: dir
        call simple_chdir(self%edir)
        call simple_mkdir(dir)
        call simple_chdir(string(dir))
    end subroutine enter_step

    !> the work project under its own file name in the working directory, the particles of mask selected in both fields
    !! and the others not, with its own copy of the current sigma2 store (a child may commit to it)
    subroutine write_copy( self, mask )
        use simple_fileio, only: basename, simple_copy_file
        class(ptcl_sieve2), intent(inout) :: self
        logical,            intent(in)    :: mask(:)
        type(string) :: base
        integer      :: iptcl, state
        do iptcl = 1,self%nptcls
            state = merge(1, 0, mask(iptcl))
            call self%proj%os_ptcl2D%set_state(iptcl, state)
            call self%proj%os_ptcl3D%set_state(iptcl, state)
        enddo
        if( self%sigma2%is_allocated() )then
            base = basename(self%sigma2)
            call simple_copy_file(self%sigma2, base)
            call self%proj%set_sigma2_state_path(base)
        endif
        call self%proj%update_projinfo(self%copy_name)
        call self%proj%write(self%copy_name)
    end subroutine write_copy

    !> the class averages, their FRCs and the sigma2 store a step registered, as absolute paths (in the step's folder)
    subroutine adopt_artifacts( self, res )
        use simple_syslib, only: simple_abspath
        class(ptcl_sieve2), intent(inout) :: self
        class(sp_project),  intent(inout) :: res
        type(string) :: path
        real         :: smpd
        integer      :: ncls
        logical      :: found
        call res%get_cavgs_stk(path, ncls, smpd)
        if( ncls /= self%s%ncls ) THROW_HARD('the class averages hold another class count; ptcl_sieve2')
        self%cavgs = simple_abspath(path)
        call res%get_frcs(path, 'frc2D')
        self%frcs = simple_abspath(path)
        call res%get_sigma2_state_path(path, found)
        if( .not. found ) THROW_HARD('the step registered no sigma2 store; ptcl_sieve2')
        self%sigma2 = simple_abspath(path)
    end subroutine adopt_artifacts

    !> one line of the run's report: the step, its wall time, the particles ingested and accepted (class not junk by
    !! the last labels), the junk classes and, for a round, kept (%): the new particles left in the student's class
    !! (all, then by their predicted class not junk and junk) and the earlier ones searched again left in theirs, and
    !! refcorr, the class averages' weighted correlation with the previous ones
    subroutine report_step( self, step, wall, kept, refcorr )
        use simple_defs_sieve2, only: SIEVE2_REPORT
        use simple_fileio,      only: fopen, fclose, fileiochk, filepath
        use simple_syslib,      only: file_exists
        class(ptcl_sieve2),   intent(inout) :: self
        character(len=*),     intent(in)    :: step
        real(timer_int_kind), intent(in)    :: wall
        real, optional,       intent(in)    :: kept(4), refcorr
        character(len=8) :: kept_str(4), corr_str, acc_str, junk_str
        logical          :: l_labels
        integer          :: i
        type(string)     :: fname
        integer          :: funit, ios, iptcl, icls, naccepted, njunk
        logical          :: l_new
        naccepted = 0
        njunk     = 0
        l_labels  = allocated(self%bad)
        if( l_labels )then
            do iptcl = 1,self%nptcls
                if( .not. self%ingested(iptcl) ) cycle
                icls = self%proj%os_ptcl2D%get_class(iptcl)
                if( icls < 1 .or. icls > size(self%bad) ) cycle
                if( .not. self%bad(icls) ) naccepted = naccepted + 1
            enddo
        endif
        fname = filepath(self%edir, SIEVE2_REPORT)
        l_new = .not. file_exists(fname)
        call fopen(funit, fname, status='unknown', action='write', position='append', iostat=ios)
        call fileiochk('ptcl_sieve2; open '//fname%to_char(), ios)
        if( l_new ) write(funit,'(A)') '# step wall_s ingested accepted junk_classes new_kept_pct new_nonjunk_kept_pct &
            &new_junk_kept_pct searched_kept_pct ref_corr'
        if( allocated(self%bad) ) njunk = count(self%bad)
        kept_str = '-'
        if( present(kept) )then
            do i = 1,4
                if( kept(i) >= 0. ) write(kept_str(i),'(F6.1)') kept(i)
            enddo
        endif
        corr_str = '-'
        if( present(refcorr) )then
            if( refcorr >= -1. ) write(corr_str,'(F7.4)') refcorr
        endif
        ! before the labelling, accepted and junk are unknown
        acc_str  = '-'
        junk_str = '-'
        if( l_labels )then
            acc_str  = int2str(naccepted)
            junk_str = int2str(njunk)
        endif
        write(funit,'(A,1X,F10.1,1X,I0,7(1X,A))') step, wall, count(self%ingested), trim(acc_str), trim(junk_str), &
            &(trim(adjustl(kept_str(i))), i = 1,4), trim(adjustl(corr_str))
        call fclose(funit)
    end subroutine report_step

    !> the class averages' change over a round: each class populated before and after (rows of mask0 by cls0, of mask1
    !! by cls1), its new average against its previous one inside the mask and up to lp, weighted by its new
    !! population; -2 when no class compares
    real function refs_change( self, stk0, stk1, mask0, cls0, mask1, cls1, lp )
        use simple_image, only: image
        class(ptcl_sieve2), intent(in) :: self
        class(string),      intent(in) :: stk0, stk1
        logical,            intent(in) :: mask0(:), mask1(:)
        integer,            intent(in) :: cls0(:), cls1(:)
        real,               intent(in) :: lp
        type(image)          :: img0, img1
        integer, allocatable :: pops0(:), pops1(:)
        real    :: wsum, csum
        integer :: icls, ncls
        ncls = self%s%ncls
        call self%class_pops(mask0, cls0, ncls, pops0)
        call self%class_pops(mask1, cls1, ncls, pops1)
        call img0%new([self%box, self%box, 1], self%smpd)
        call img1%new([self%box, self%box, 1], self%smpd)
        wsum = 0.
        csum = 0.
        do icls = 1,ncls
            if( pops0(icls) < 1 .or. pops1(icls) < 1 ) cycle
            call img0%read(stk0, icls)
            call img1%read(stk1, icls)
            call img0%mask2D_soft(self%mskdiam / (2. * self%smpd))
            call img1%mask2D_soft(self%mskdiam / (2. * self%smpd))
            csum = csum + real(pops1(icls)) * img0%corr(img1, lp_dyn=lp)
            wsum = wsum + real(pops1(icls))
        enddo
        call img0%kill
        call img1%kill
        refs_change = -2.
        if( wsum > 0. ) refs_change = csum / wsum
    end function refs_change

    ! ---- the rules -------------------------------------------------------------------------------------------------

    !> the seed (stk_batch 0): the first stacks with selected particles while it holds fewer than nmics stacks and
    !! fewer than seed_nptcls particles; then batches of whole stacks in project order, each closed once it holds
    !! batch_nptcls selected particles, the tail the last; -1 for a stack without selected particles
    subroutine plan_batches( stk_nsel, nmics, seed_nptcls, batch_nptcls, stk_batch, nbatches )
        integer,              intent(in)               :: stk_nsel(:), nmics, seed_nptcls, batch_nptcls
        integer,              allocatable, intent(out) :: stk_batch(:)
        integer,              intent(out)              :: nbatches
        integer :: istk, nseed, nseed_ptcls, cnt
        allocate(stk_batch(size(stk_nsel)), source=-1)
        nseed       = 0
        nseed_ptcls = 0
        nbatches    = 0
        cnt         = 0
        do istk = 1,size(stk_nsel)
            if( stk_nsel(istk) < 1 ) cycle
            if( nseed < nmics .and. nseed_ptcls < seed_nptcls )then
                nseed           = nseed + 1
                nseed_ptcls     = nseed_ptcls + stk_nsel(istk)
                stk_batch(istk) = 0
                cycle
            endif
            if( cnt == 0 ) nbatches = nbatches + 1
            stk_batch(istk) = nbatches
            cnt = cnt + stk_nsel(istk)
            if( cnt >= batch_nptcls ) cnt = 0
        enddo
    end subroutine plan_batches

    !> whether the student is trained again after greedy round iround: it is one of rounds (none when unallocated)
    logical function retrains_after( rounds, iround )
        integer,              allocatable, intent(in) :: rounds(:)
        integer,              intent(in)              :: iround
        retrains_after = .false.
        if( allocated(rounds) ) retrains_after = any(rounds == iround)
    end function retrains_after

    !> a class is junk when it is populated, not empty and the CNN's p_junk exceeds thres
    subroutine class_labels( pjunk, empty, pops, thres, bad )
        real(dp),             intent(in)               :: pjunk(:)
        logical,              intent(in)               :: empty(:)
        integer,              intent(in)               :: pops(:)
        real(dp),             intent(in)               :: thres
        logical,              allocatable, intent(out) :: bad(:)
        if( size(empty) /= size(pjunk) .or. size(pops) /= size(pjunk) ) THROW_HARD('sizes differ; ptcl_sieve2 class_labels')
        allocate(bad(size(pjunk)))
        bad = pops > 0 .and. .not. empty .and. pjunk > thres
    end subroutine class_labels

    !> the greedy round's update fraction: the new particles and old_frac of the earlier ones, at least as many
    !! earlier ones as new, at most all; a share that grows with the set bounds how long an earlier particle waits
    real function greedy_frac( n_new, n_active, old_frac )
        integer, intent(in) :: n_new, n_active
        real,    intent(in) :: old_frac
        real :: n_old
        greedy_frac = 1.
        if( n_active < 1 ) return
        n_old       = real(n_active - n_new)
        greedy_frac = min(1., (real(n_new) + max(real(n_new), old_frac * n_old)) / real(n_active))
    end function greedy_frac

    !> a round's agreement counts: of the new particles it searched (no update before it) how many kept the class
    !! they came with, the student's, and of those whose class bad marked junk how many; of the earlier particles it
    !! searched again how many kept theirs
    subroutine round_agreement( active, cls0, cls1, cnt0, cnt1, bad, nnew, nnew_kept, njunk, njunk_kept, nold, &
        &nold_kept )
        logical, intent(in)  :: active(:), bad(:)
        integer, intent(in)  :: cls0(:), cls1(:), cnt0(:), cnt1(:)
        integer, intent(out) :: nnew, nnew_kept, njunk, njunk_kept, nold, nold_kept
        logical, allocatable :: srch(:), new(:), same(:), junk(:)
        integer :: i
        srch = active .and. cnt1 > cnt0
        new  = srch .and. cnt0 == 0
        same = cls1 == cls0
        allocate(junk(size(cls0)), source=.false.)
        do i = 1,size(cls0)
            if( cls0(i) >= 1 .and. cls0(i) <= size(bad) ) junk(i) = bad(cls0(i))
        enddo
        nnew       = count(new)
        nnew_kept  = count(new .and. same)
        njunk      = count(new .and. junk)
        njunk_kept = count(new .and. junk .and. same)
        nold       = count(srch .and. .not. new)
        nold_kept  = count(srch .and. .not. new .and. same)
    end subroutine round_agreement

    !> adds a round's new particles (no update before it, searched by it) to the class counts: npred their predicted
    !! class, nkept that class when the round left them there, ninto the class the round put them in
    subroutine class_agreement( active, cls0, cls1, cnt0, cnt1, npred, nkept, ninto )
        logical, intent(in)    :: active(:)
        integer, intent(in)    :: cls0(:), cls1(:), cnt0(:), cnt1(:)
        integer, intent(inout) :: npred(:), nkept(:), ninto(:)
        integer :: i
        do i = 1,size(cls0)
            if( .not. active(i) ) cycle
            if( cnt0(i) /= 0 .or. cnt1(i) <= cnt0(i) ) cycle
            if( cls0(i) >= 1 .and. cls0(i) <= size(npred) )then
                npred(cls0(i)) = npred(cls0(i)) + 1
                if( cls1(i) == cls0(i) ) nkept(cls0(i)) = nkept(cls0(i)) + 1
            endif
            if( cls1(i) >= 1 .and. cls1(i) <= size(ninto) ) ninto(cls1(i)) = ninto(cls1(i)) + 1
        enddo
    end subroutine class_agreement

    !> the classes by the share of their predictions kept, best first, ties in class order; classes without a
    !! prediction last
    function kept_order( npred, nkept ) result( order )
        use simple_srch_sort_loc, only: scores2order_stable
        integer, intent(in)  :: npred(:), nkept(:)
        integer, allocatable :: order(:)
        real(dp), allocatable :: share(:)
        integer :: c
        allocate(share(size(npred)), source=-1._dp)
        do c = 1,size(npred)
            if( npred(c) > 0 ) share(c) = real(nkept(c), dp) / real(npred(c), dp)
        enddo
        order = scores2order_stable(share)
    end function kept_order

    !> the classes' populations over the rows of mask (classes outside 1..ncls not counted)
    subroutine class_pops( mask, classes, ncls, pops )
        logical,              intent(in)  :: mask(:)
        integer,              intent(in)  :: classes(:), ncls
        integer, allocatable, intent(out) :: pops(:)
        integer :: i
        allocate(pops(ncls), source=0)
        do i = 1,size(classes)
            if( .not. mask(i) ) cycle
            if( classes(i) >= 1 .and. classes(i) <= ncls ) pops(classes(i)) = pops(classes(i)) + 1
        enddo
    end subroutine class_pops

    !> the rounds running that left the class averages at a weighted correlation of thres or more, after one with refcorr
    integer function freeze_count( ncons, refcorr, thres )
        integer, intent(in) :: ncons
        real,    intent(in) :: refcorr, thres
        freeze_count = 0
        if( refcorr >= thres ) freeze_count = ncons + 1
    end function freeze_count

    !> a frozen round's update fraction: the new particles, which hold the lowest update count, alone
    real function newonly_frac( n_new, n_active )
        integer, intent(in) :: n_new, n_active
        newonly_frac = 1.
        if( n_active > 0 ) newonly_frac = min(1., real(n_new) / real(n_active))
    end function newonly_frac

    !> a share in percent as text, one decimal; - for none (negative)
    function pct_str( p ) result( str )
        real, intent(in) :: p
        character(len=:), allocatable :: str
        character(len=8) :: buf
        if( p < 0. )then
            str = '-'
        else
            write(buf,'(F6.1)') p
            str = trim(adjustl(buf))
        endif
    end function pct_str

    !> the share n of ntot in percent; -1 for none
    real function pct( n, ntot )
        integer, intent(in) :: n, ntot
        pct = -1.
        if( ntot > 0 ) pct = 100. * real(n) / real(ntot)
    end function pct

    !> a greedy round's low-pass limit (A): the class FRCs' estimate, never finer than lpstop
    real function round_lowpass( lp_frc, lpstop )
        real, intent(in) :: lp_frc, lpstop
        round_lowpass = max(lp_frc, lpstop)
    end function round_lowpass

    !> state 1 for a selected, ingested particle whose class is not junk, 0 otherwise (a particle without a class too)
    subroutine final_states( sel, ingested, classes, bad, states )
        logical,              intent(in)               :: sel(:), ingested(:), bad(:)
        integer,              intent(in)               :: classes(:)
        integer,              allocatable, intent(out) :: states(:)
        integer :: iptcl, icls
        allocate(states(size(sel)), source=0)
        do iptcl = 1,size(sel)
            if( .not. sel(iptcl) )      cycle
            if( .not. ingested(iptcl) ) cycle
            icls = classes(iptcl)
            if( icls < 1 .or. icls > size(bad) ) cycle
            if( .not. bad(icls) ) states(iptcl) = 1
        enddo
    end subroutine final_states

    !> stops when no populated class is left that is not junk
    subroutine check_classes_left( bad, pops )
        logical, intent(in) :: bad(:)
        integer, intent(in) :: pops(:)
        if( count(pops > 0 .and. .not. bad) == 0 ) THROW_HARD('no class left: every populated class is junk; particle_sieving_2')
    end subroutine check_classes_left

    !> rows holds particles of both halves (eo 0 and 1)
    logical function both_halves( rows, eo )
        logical, intent(in) :: rows(:)
        integer, intent(in) :: eo(:)
        both_halves = any(rows .and. eo == 0) .and. any(rows .and. eo == 1)
    end function both_halves

    !> stops when a selected particle was never classified
    subroutine check_none_pending( sel, ingested )
        logical, intent(in) :: sel(:), ingested(:)
        if( any(sel .and. .not. ingested) ) THROW_HARD('a selected particle was never classified; particle_sieving_2')
    end subroutine check_none_pending

end module simple_ptcl_sieve2
