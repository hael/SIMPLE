!@descr: the class student: a teacher's class averages as the classes of the polar student's core, the candidate classes holding most of the posterior, the class and in-plane angle a readout gives, and the model file's class keys
! The core (simple_polar_student) holds the network, the features, the band, the training and the model file's common
! part; this class adds what is 2D: the classes without mirrors (the 2D search tests none), their states (a rejected
! class stays a class, so that a prediction into it is a rejection signal), the candidates and the read-out pose
! (doc/algorithms/class_student.md).
module simple_class_student
use, intrinsic :: iso_fortran_env, only: int64
use simple_defs,          only: sp, dp
use simple_error,         only: simple_exception
use simple_string,        only: string
use simple_polar_feats,   only: polar_feat
use simple_polar_student, only: polar_student, polar_student_header, student_train_profile, polar_student_targets
use simple_math,          only: hpsort
use simple_defs_hhcnn,    only: KMAX_CLASS
implicit none

public :: class_student, class_pred, CLASS_STUDENT_VERSION, run_all_class_student_tests
private
#include "simple_local_flags.inc"

integer,  parameter :: CLASS_STUDENT_VERSION = 1
real(sp), parameter :: CAND_MASS             = 0.99  !< the candidates hold this share of the posterior
character(len=*), parameter :: MODEL_TITLE   = 'class student model'

!> the record of one prediction: the mode, the candidates (the smallest set of classes holding CAND_MASS, at most
!! KMAX_CLASS) and their mass
type :: class_pred
    integer  :: mode  = 0     !< the class of largest mass
    real(sp) :: pmax  = 0.    !< its mass
    integer  :: k99   = 0     !< classes holding CAND_MASS
    integer  :: ncand = 0     !< the candidates (k99 capped)
    real(sp) :: conf  = 0.    !< the mass inside the candidates
    integer  :: cands(KMAX_CLASS) = 0
end type class_pred

type :: class_student
    private
    type(polar_student), allocatable :: core   !< the network, the features, the band, the training, the common model file
    integer :: ncls = 0
    integer,  allocatable :: states(:)          !< (ncls) the classes' states in the teacher's project (0: rejected)
    integer(int64) :: cavgs_hash = 0_int64      !< the teacher's class averages
    logical :: exists = .false.
  contains
    procedure :: new
    procedure :: init_readout
    procedure :: train
    procedure :: prepare_predict
    procedure :: predict_feat
    procedure :: predict_tensor
    procedure :: cands_of
    procedure :: rotation_of
    procedure :: write
    procedure :: read
    procedure :: get_ncls
    procedure :: get_state
    procedure :: get_pftsz
    procedure :: get_kfromto
    procedure :: get_lp
    procedure :: get_box
    procedure :: get_box_crop
    procedure :: get_smpd
    procedure :: get_mskdiam
    procedure :: get_with_ctf
    procedure :: get_weights
    procedure :: get_record
    procedure :: get_train_profile
    procedure :: set_provenance
    procedure :: set_threshold
    procedure :: get_threshold
    procedure :: kill
end type class_student

interface
    module subroutine run_all_class_student_tests()
    end subroutine run_all_class_student_tests
end interface

contains

    !> A student for training: the calculator's geometry, the band's low-pass lp, the sampling, the mask, the shell
    !! weights wk the rows were whitened with (frozen), whether the CTF is applied, the classes and their states, the
    !! network's width and the seed its weights start from
    subroutine new( self, pftsz, kfromto, box, box_crop, smpd, mskdiam, lp, wk, with_ctf, ncls, states, width, seed )
        class(class_student), intent(inout) :: self
        integer,              intent(in)    :: pftsz, kfromto(2), box, box_crop, ncls, states(ncls), width
        real(sp),             intent(in)    :: smpd, mskdiam, lp
        real(sp),             intent(in)    :: wk(kfromto(1):kfromto(2),2)
        logical,              intent(in)    :: with_ctf
        integer(int64),       intent(in)    :: seed
        call self%kill
        self%ncls = ncls
        allocate(self%states(ncls), source=states)
        allocate(self%core)
        ! classes without mirrors: the 2D search tests none, so the readout scores the map alone
        call self%core%new(pftsz, kfromto, box, box_crop, smpd, mskdiam, lp, wk, with_ctf, ncls, .false., width, seed)
        self%exists = .true.
    end subroutine new

    !> the readout's templates from the class averages: xcavgs(pftsz_net, nk, 2, ncls) their tensors through the features
    !! (no CTF: the ring-averaged CTF squared one); done on every new network before its training
    subroutine init_readout( self, xcavgs )
        class(class_student), intent(inout) :: self
        real(sp),             intent(in)    :: xcavgs(:,:,:,:)
        if( .not. self%exists ) THROW_HARD('student not constructed; class_student%init_readout')
        call self%core%init_templates(xcavgs)
    end subroutine init_readout

    !> The training on the teacher's rows of n particles (the core's training): x, ctf2, side the tensors, dists(ncls, n)
    !! the rows' losses and temps(n) their temperatures, whose posterior over the classes is each row's soft target
    subroutine train( self, n, x, ctf2, side, dists, temps, held, prob_beta, shift_aug_px, max_epochs, patience, &
        &nthr, log_unit )
        class(class_student), intent(inout) :: self
        integer,              intent(in)    :: n, max_epochs, patience, nthr, log_unit
        real(sp),             intent(in)    :: x(:,:,:,:), ctf2(:,:), side(n)
        logical,              intent(in)    :: held(n)
        real(sp),             intent(in)    :: dists(self%ncls, n), temps(n), prob_beta, shift_aug_px
        real(sp), allocatable :: targets(:,:)
        integer :: isample
        if( .not. self%exists ) THROW_HARD('student not constructed; class_student%train')
        allocate(targets(self%ncls, n))
        !$omp parallel do default(shared) private(isample) schedule(static) proc_bind(close)
        do isample = 1,n
            call polar_student_targets(self%ncls, dists(:,isample), temps(isample), targets(:,isample))
        enddo
        !$omp end parallel do
        call self%core%train(n, x, ctf2, side, targets, held, prob_beta, shift_aug_px, max_epochs, patience, nthr, log_unit)
        deallocate(targets)
    end subroutine train

    !> the workspaces of a team of nthr threads for predict_feat
    subroutine prepare_predict( self, nthr )
        class(class_student), intent(inout) :: self
        integer,              intent(in)    :: nthr
        if( .not. self%exists ) THROW_HARD('student not constructed; class_student%prepare_predict')
        call self%core%prepare_predict(nthr)
    end subroutine prepare_predict

    !> the posterior over the classes of one feature (built with this student's geometry and weights) on thread ithr
    subroutine predict_feat( self, ithr, feat, prob )
        class(class_student), intent(inout) :: self
        integer,              intent(in)    :: ithr
        type(polar_feat),     intent(in)    :: feat
        real(sp),             intent(out)   :: prob(self%ncls)
        call self%core%predict_tensor(ithr, feat%x, feat%ctf2, feat%logpow, prob)
    end subroutine predict_feat

    !> the posterior over the classes of a packed tensor with its ring-averaged CTF squared and raw side scalar
    subroutine predict_tensor( self, ithr, x, ctf2, side, prob )
        class(class_student), intent(inout) :: self
        integer,              intent(in)    :: ithr
        real(sp),             intent(in)    :: x(:,:,:), ctf2(:), side
        real(sp),             intent(out)   :: prob(self%ncls)
        call self%core%predict_tensor(ithr, x, ctf2, side, prob)
    end subroutine predict_tensor

    !> the mode, the candidates (the smallest set of classes holding CAND_MASS, capped at KMAX_CLASS, in descending
    !! mass) and the mass inside them
    subroutine cands_of( self, prob, pred )
        class(class_student), intent(in)  :: self
        real(sp),             intent(in)  :: prob(self%ncls)
        type(class_pred),     intent(out) :: pred
        real(sp) :: w(self%ncls)
        integer  :: inds(self%ncls), i, j
        real(dp) :: csum
        w    = prob
        inds = [(i, i = 1,self%ncls)]
        call hpsort(w, inds)
        pred%mode = inds(self%ncls)
        pred%pmax = w(self%ncls)
        csum      = 0.d0
        pred%k99  = 0
        do i = self%ncls,1,-1
            csum     = csum + real(w(i),dp)
            pred%k99 = pred%k99 + 1
            if( csum >= real(CAND_MASS,dp) ) exit
        enddo
        pred%ncand = min(pred%k99, KMAX_CLASS)
        pred%cands = 0
        pred%conf  = 0.
        do j = 1,pred%ncand
            pred%cands(j) = inds(self%ncls - j + 1)
            pred%conf     = pred%conf + w(self%ncls - j + 1)
        enddo
    end subroutine cands_of

    !> the in-plane angle the readout of thread ithr's last prediction gives class icls, without any search: the readout
    !! turns the sample's map by ang from the template, so e3 = 360 - ang, as the search's rotation index
    real(sp) function rotation_of( self, ithr, icls ) result( e3 )
        class(class_student), intent(in) :: self
        integer,              intent(in) :: ithr, icls
        real(sp) :: ang
        logical  :: refl
        call self%core%class_rotation(ithr, icls, ang, refl)
        e3 = modulo(360. - ang, 360.)
    end function rotation_of

    ! ===== THE MODEL FILE =====

    !> the header fname (key=value) and the weights beside it (the core's model file), with the classes' keys and the
    !! caller's after the core's
    subroutine write( self, fname, extra_keys, extra_vals )
        class(class_student),       intent(in) :: self
        class(string),              intent(in) :: fname
        character(len=*), optional, intent(in) :: extra_keys(:), extra_vals(:)
        integer, parameter :: NOWN = 5
        character(len=64),   allocatable :: allkeys(:)
        character(len=:),    allocatable :: states_str
        character(len=4096), allocatable :: allvals(:)
        character(len=16) :: tmp
        integer :: nextra, i
        if( .not. self%exists ) THROW_HARD('student not constructed; class_student%write')
        nextra = 0
        if( present(extra_keys) ) nextra = size(extra_keys)
        allocate(allkeys(NOWN + nextra), allvals(NOWN + nextra))
        states_str = ''
        do i = 1,self%ncls
            write(tmp,'(I0)') self%states(i)
            if( i > 1 ) states_str = states_str//','
            states_str = states_str//trim(tmp)
        enddo
        allkeys(1) = 'ncls';       write(allvals(1),'(I0)')   self%ncls
        allkeys(2) = 'kmax_class'; write(allvals(2),'(I0)')   KMAX_CLASS
        allkeys(3) = 'cand_mass';  write(allvals(3),'(F6.3)') CAND_MASS
        allkeys(4) = 'cavgs_hash'; write(allvals(4),'(I0)')   self%cavgs_hash
        allkeys(5) = 'states';     allvals(5) = states_str
        if( nextra > 0 )then
            allkeys(NOWN+1:) = extra_keys
            allvals(NOWN+1:) = extra_vals
        endif
        call self%core%write(fname, MODEL_TITLE, CLASS_STUDENT_VERSION, allkeys, allvals)
        deallocate(allkeys, allvals)
    end subroutine write

    !> the student from its header and the weights the header names beside it; the weights must be the network the
    !! header describes
    subroutine read( self, fname )
        class(class_student), intent(inout) :: self
        class(string),        intent(in)    :: fname
        type(polar_student_header) :: hdr
        character(len=:), allocatable :: key, val
        integer, allocatable :: states(:)
        integer(int64) :: cavgs_hash
        integer :: ncls, i
        logical :: l_states
        call hdr%read(fname)
        if( hdr%version /= CLASS_STUDENT_VERSION ) THROW_HARD('model_version mismatch; class_student%read')
        ncls       = 0
        cavgs_hash = 0_int64
        l_states   = .false.
        do i = 1,hdr%nextra
            key = hdr%keys(i)%to_char()
            val = hdr%vals(i)%to_char()
            if( key == 'ncls' ) read(val,*) ncls
        enddo
        if( ncls < 2 ) THROW_HARD('incomplete classes in the model header; class_student%read')
        allocate(states(ncls), source=1)
        do i = 1,hdr%nextra
            key = hdr%keys(i)%to_char()
            val = hdr%vals(i)%to_char()
            select case(key)
                case('cavgs_hash'); read(val,*) cavgs_hash
                case('states')
                    read(val,*) states
                    l_states = .true.
                case default
            end select
        enddo
        if( .not. l_states ) THROW_HARD('the model header lacks its class states; class_student%read')
        call self%new(hdr%pftsz, hdr%kfromto, hdr%box, hdr%box_crop, hdr%smpd, hdr%mskdiam, hdr%lp, hdr%wk, &
            &hdr%with_ctf, ncls, states, hdr%width, hdr%seed)
        call self%core%load(hdr, fname)
        self%cavgs_hash = cavgs_hash
        call hdr%kill
        deallocate(states)
    end subroutine read

    ! ===== GETTERS AND SETTERS =====

    pure integer function get_ncls( self )
        class(class_student), intent(in) :: self
        get_ncls = self%ncls
    end function get_ncls

    pure integer function get_state( self, icls )
        class(class_student), intent(in) :: self
        integer,              intent(in) :: icls
        get_state = self%states(icls)
    end function get_state

    pure integer function get_pftsz( self )
        class(class_student), intent(in) :: self
        get_pftsz = self%core%get_pftsz()
    end function get_pftsz

    pure function get_kfromto( self ) result( kfromto )
        class(class_student), intent(in) :: self
        integer :: kfromto(2)
        kfromto = self%core%get_kfromto()
    end function get_kfromto

    pure real(sp) function get_lp( self )
        class(class_student), intent(in) :: self
        get_lp = self%core%get_lp()
    end function get_lp

    pure integer function get_box( self )
        class(class_student), intent(in) :: self
        get_box = self%core%get_box()
    end function get_box

    pure integer function get_box_crop( self )
        class(class_student), intent(in) :: self
        get_box_crop = self%core%get_box_crop()
    end function get_box_crop

    pure real(sp) function get_smpd( self )
        class(class_student), intent(in) :: self
        get_smpd = self%core%get_smpd()
    end function get_smpd

    pure real(sp) function get_mskdiam( self )
        class(class_student), intent(in) :: self
        get_mskdiam = self%core%get_mskdiam()
    end function get_mskdiam

    pure logical function get_with_ctf( self )
        class(class_student), intent(in) :: self
        get_with_ctf = self%core%get_with_ctf()
    end function get_with_ctf

    subroutine get_weights( self, wk )
        class(class_student),  intent(in)  :: self
        real(sp), allocatable, intent(out) :: wk(:,:)
        call self%core%get_weights(wk)
    end subroutine get_weights

    !> the training record for the caller's reports
    subroutine get_record( self, nptcls_train, nheldout, epochs_run, best_epoch, loss_train, loss_held, mode_held, &
        &time_train_s, nweights )
        class(class_student), intent(in)  :: self
        integer,              intent(out) :: nptcls_train, nheldout, epochs_run, best_epoch, nweights
        real(sp),             intent(out) :: loss_train, loss_held, mode_held, time_train_s
        call self%core%get_record(nptcls_train, nheldout, epochs_run, best_epoch, loss_train, loss_held, mode_held, &
            &time_train_s, nweights)
    end subroutine get_record

    !> the cost record of the last training
    function get_train_profile( self ) result( prof )
        class(class_student), intent(in) :: self
        type(student_train_profile) :: prof
        prof = self%core%get_train_profile()
    end function get_train_profile

    !> the teacher's provenance, recorded in the model file: the core's and the class averages' hash
    subroutine set_provenance( self, teacher_project, selection_hash, cavgs_hash, holdout_rule )
        class(class_student), intent(inout) :: self
        character(len=*),     intent(in)    :: teacher_project, holdout_rule
        integer(int64),       intent(in)    :: selection_hash, cavgs_hash
        call self%core%set_provenance(teacher_project, selection_hash, holdout_rule)
        self%cavgs_hash = cavgs_hash
    end subroutine set_provenance

    !> the confidence threshold set on the held-out rows: above it the teacher's class lies among the candidates in
    !! conf_fraction of the held-out particles
    subroutine set_threshold( self, conf_threshold, conf_fraction, heldout_seeded, heldout_recovered )
        class(class_student), intent(inout) :: self
        real(sp),             intent(in)    :: conf_threshold, conf_fraction, heldout_seeded, heldout_recovered
        call self%core%set_threshold(conf_threshold, conf_fraction, heldout_seeded, heldout_recovered)
    end subroutine set_threshold

    pure subroutine get_threshold( self, conf_threshold, conf_fraction, heldout_seeded, heldout_recovered )
        class(class_student), intent(in)  :: self
        real(sp),             intent(out) :: conf_threshold, conf_fraction, heldout_seeded, heldout_recovered
        call self%core%get_threshold(conf_threshold, conf_fraction, heldout_seeded, heldout_recovered)
    end subroutine get_threshold

    subroutine kill( self )
        class(class_student), intent(inout) :: self
        if( allocated(self%core) )then
            call self%core%kill
            deallocate(self%core)
        endif
        if( allocated(self%states) ) deallocate(self%states)
        self%ncls       = 0
        self%cavgs_hash = 0_int64
        self%exists     = .false.
    end subroutine kill

end module simple_class_student
