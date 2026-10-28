!@descr: bank-free weighted, regularized pairwise-logistic fit and the dataset-role scoring shared by the quality models
! The caller owns its feature bank: it gives the problem the feature indices of a policy, one row of features per
! observation, the label and the weight, and scores each dataset itself; this module knows nothing of features.
module simple_quality_logistic
use simple_error,       only: simple_exception
use simple_optimizer,   only: optimizer
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
implicit none
private
#include "simple_local_flags.inc"

public :: quality_logistic_problem
public :: dataset_learn_role, logistic_dataset_class_weights, robust_macro_score, stable_sigmoid_8
public :: QUALITY_MAX_FEATS, QUALITY_MAX_INTERACTIONS
public :: LEARN_ROLE_SKIP, LEARN_ROLE_BALANCED, LEARN_ROLE_RECALL_ONLY, LEARN_ROLE_SPECIFICITY_ONLY
public :: LOGISTIC_LAMBDAS, LOGISTIC_THRESHOLDS, LEARN_ROBUST_TAIL_FRAC, LEARN_ROBUST_TAIL_WEIGHT
public :: LEARN_BALANCED_GOOD_LOSS_WEIGHT, LEARN_BALANCED_BAD_LOSS_WEIGHT

integer, parameter :: QUALITY_MAX_FEATS        = 16   !< the largest feature bank of any quality model
integer, parameter :: QUALITY_MAX_INTERACTIONS = (QUALITY_MAX_FEATS * (QUALITY_MAX_FEATS - 1)) / 2
integer, parameter :: LEARN_ROLE_SKIP             = 0
integer, parameter :: LEARN_ROLE_BALANCED         = 1
integer, parameter :: LEARN_ROLE_RECALL_ONLY      = 2
integer, parameter :: LEARN_ROLE_SPECIFICITY_ONLY = 3
real,    parameter :: LEARN_ROBUST_TAIL_FRAC          = 0.25
real,    parameter :: LEARN_ROBUST_TAIL_WEIGHT        = 0.50
real,    parameter :: LEARN_BALANCED_GOOD_LOSS_WEIGHT = 0.65
real,    parameter :: LEARN_BALANCED_BAD_LOSS_WEIGHT  = 1.0 - LEARN_BALANCED_GOOD_LOSS_WEIGHT
real,    parameter :: LOGISTIC_LAMBDAS(7)     = [1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3, &
                                                 1.0e-2, 3.0e-2, 1.0e-1]
real,    parameter :: LOGISTIC_THRESHOLDS(20) = [0.02, 0.05, 0.075, 0.10, 0.15, 0.20, &
                                                 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, &
                                                 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, &
                                                 0.85, 0.90]
real,    parameter :: LOGISTIC_COEFF_ABS_BOUND = 20.0

!> the design: intercept, the policy's features, their pairwise products, optionally one extra (relational) column
type :: quality_logistic_problem
    private
    integer, public :: nobs           = 0
    integer, public :: ndim           = 0
    integer, public :: n_linear       = 0
    integer, public :: n_interactions = 0
    integer, public :: linear_features(QUALITY_MAX_FEATS) = 0
    integer, public :: interaction_terms(QUALITY_MAX_INTERACTIONS,2) = 0
    logical, public :: has_extra      = .false.
    real(kind=8), public :: lambda       = 0.0d0
    real(kind=8), public :: total_weight = 0.0d0
    real(kind=8), allocatable :: design(:,:)
    real(kind=8), allocatable :: labels(:)
    real(kind=8), allocatable :: weights(:)
contains
    procedure :: new
    procedure :: set_row
    procedure :: finalize
    procedure :: fit
    procedure :: fdf
    procedure :: kill
end type quality_logistic_problem

contains

    !> feature indices of the policy (into the caller's bank), the number of observations, the extra column, lambda
    subroutine new( self, linear_features, nobs, has_extra, lambda )
        class(quality_logistic_problem), intent(inout) :: self
        integer,                         intent(in)    :: linear_features(:), nobs
        logical,                         intent(in)    :: has_extra
        real,                            intent(in)    :: lambda
        integer :: ilinear, jlinear
        call self%kill()
        self%n_linear = size(linear_features)
        if( self%n_linear < 1 ) THROW_HARD('empty feature policy; quality_logistic_problem%new')
        if( self%n_linear > QUALITY_MAX_FEATS ) THROW_HARD('too many features; quality_logistic_problem%new')
        if( nobs < 1 ) THROW_HARD('no observations; quality_logistic_problem%new')
        self%linear_features(1:self%n_linear) = linear_features
        do ilinear = 1, self%n_linear - 1
            do jlinear = ilinear + 1, self%n_linear
                self%n_interactions = self%n_interactions + 1
                self%interaction_terms(self%n_interactions,1) = linear_features(ilinear)
                self%interaction_terms(self%n_interactions,2) = linear_features(jlinear)
            end do
        end do
        self%has_extra = has_extra
        self%ndim   = 1 + self%n_linear + self%n_interactions + merge(1, 0, has_extra)
        self%lambda = real(lambda, kind=8)
        self%nobs   = nobs
        allocate(self%design(nobs, self%ndim), source=0.0d0)
        allocate(self%labels(nobs),  source=0.0d0)
        allocate(self%weights(nobs), source=0.0d0)
    end subroutine new

    !> row iobs from the full feature vector of the caller's bank (single precision, as the training tables hold it)
    subroutine set_row( self, iobs, features, good, weight, extra )
        class(quality_logistic_problem), intent(inout) :: self
        integer,                         intent(in)    :: iobs
        real,                            intent(in)    :: features(:)
        logical,                         intent(in)    :: good
        real(kind=8),                    intent(in)    :: weight
        real, optional,                  intent(in)    :: extra
        integer :: ilinear, iterm, ifeat, jfeat
        if( iobs < 1 .or. iobs > self%nobs ) THROW_HARD('row out of range; quality_logistic_problem%set_row')
        self%design(iobs,1) = 1.0d0
        do ilinear = 1, self%n_linear
            ifeat = self%linear_features(ilinear)
            self%design(iobs, 1 + ilinear) = real(features(ifeat), kind=8)
        end do
        do iterm = 1, self%n_interactions
            ifeat = self%interaction_terms(iterm,1)
            jfeat = self%interaction_terms(iterm,2)
            self%design(iobs, 1 + self%n_linear + iterm) = real(features(ifeat) * features(jfeat), kind=8)
        end do
        if( self%has_extra )then
            if( .not. present(extra) ) THROW_HARD('missing extra column; quality_logistic_problem%set_row')
            self%design(iobs,self%ndim) = real(extra, kind=8)
        endif
        if( good )then
            self%labels(iobs) = 1.0d0
        else
            self%labels(iobs) = 0.0d0
        endif
        self%weights(iobs) = weight
    end subroutine set_row

    subroutine finalize( self )
        class(quality_logistic_problem), intent(inout) :: self
        self%total_weight = sum(self%weights)
        if( self%total_weight <= 0.0d0 ) THROW_HARD('nonpositive total weight; quality_logistic_problem%finalize')
    end subroutine finalize

    !> L-BFGS-B from the weighted log-odds intercept, coefficients bounded by LOGISTIC_COEFF_ABS_BOUND
    subroutine fit( self, solution, objective )
        class(quality_logistic_problem), intent(inout) :: self
        real(kind=8), allocatable,       intent(inout) :: solution(:)
        real,                            intent(out)   :: objective
        type(opt_factory)         :: ofac
        type(opt_spec)            :: ospec
        class(optimizer), pointer :: opt_obj
        real, allocatable         :: limits(:,:)
        real(kind=8)              :: pos_weight, neg_weight
        integer :: i
        if( allocated(solution) ) deallocate(solution)
        allocate(solution(self%ndim), source=0.0d0)
        allocate(limits(self%ndim,2), source=0.0)
        limits(:,1) = -LOGISTIC_COEFF_ABS_BOUND
        limits(:,2) =  LOGISTIC_COEFF_ABS_BOUND
        pos_weight = sum(self%weights, mask=self%labels > 0.5d0)
        neg_weight = sum(self%weights, mask=self%labels <= 0.5d0)
        solution(1) = log(max(pos_weight, 1.0d-12) / max(neg_weight, 1.0d-12))
        opt_obj => null()
        call ospec%specify('lbfgsb', self%ndim, ftol=1e-7, gtol=1e-7, maxits=250, &
            factr=1.0d+7, pgtol=1.0d-5, limits=limits)
        call ofac%new(ospec, opt_obj)
        call ospec%set_costfun_8(logistic_cost_wrapper)
        call ospec%set_gcostfun_8(logistic_gradient_wrapper)
        call ospec%set_fdfcostfun_8(logistic_fdf_wrapper)
        ospec%x   = real(solution)
        ospec%x_8 = solution
        call opt_obj%minimize(ospec, self, objective)
        do i = 1, self%ndim
            solution(i) = ospec%x_8(i)
        end do
        call opt_obj%kill()
        deallocate(opt_obj)
        call ospec%kill()
        deallocate(limits)
    end subroutine fit

    !> weighted mean log-loss plus the ridge on every coefficient but the intercept, and its gradient
    subroutine fdf( self, vec, f, grad, D )
        class(quality_logistic_problem), intent(in)  :: self
        integer,                         intent(in)  :: D
        real(kind=8),                    intent(in)  :: vec(D)
        real(kind=8),                    intent(out) :: f, grad(D)
        integer :: iobs, idim
        real(kind=8) :: eta, prob, resid, loss
        if( D /= self%ndim ) THROW_HARD('dimension mismatch; quality_logistic_problem%fdf')
        if( self%total_weight <= 0.0d0 ) THROW_HARD('nonpositive total weight; quality_logistic_problem%fdf')
        f    = 0.0d0
        grad = 0.0d0
        do iobs = 1, self%nobs
            eta  = dot_product(vec, self%design(iobs,:))
            prob = stable_sigmoid_8(eta)
            if( eta >= 0.0d0 )then
                loss = (1.0d0 - self%labels(iobs)) * eta + log(1.0d0 + exp(-min(eta, 80.0d0)))
            else
                loss = -self%labels(iobs) * eta + log(1.0d0 + exp(max(eta, -80.0d0)))
            endif
            f = f + self%weights(iobs) * loss
            resid = self%weights(iobs) * (prob - self%labels(iobs))
            do idim = 1, D
                grad(idim) = grad(idim) + resid * self%design(iobs,idim)
            end do
        end do
        f    = f    / self%total_weight
        grad = grad / self%total_weight
        do idim = 2, D
            f = f + 0.5d0 * self%lambda * vec(idim)**2
            grad(idim) = grad(idim) + self%lambda * vec(idim)
        end do
    end subroutine fdf

    subroutine kill( self )
        class(quality_logistic_problem), intent(inout) :: self
        if( allocated(self%design)  ) deallocate(self%design)
        if( allocated(self%labels)  ) deallocate(self%labels)
        if( allocated(self%weights) ) deallocate(self%weights)
        self%nobs              = 0
        self%ndim              = 0
        self%n_linear          = 0
        self%n_interactions    = 0
        self%linear_features   = 0
        self%interaction_terms = 0
        self%has_extra         = .false.
        self%lambda            = 0.0d0
        self%total_weight      = 0.0d0
    end subroutine kill

    function logistic_cost_wrapper( fun_self, vec, D ) result( cost )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(in)    :: vec(D)
        real(kind=8)                :: cost
        real(kind=8) :: grad(D)
        select type(problem => fun_self)
        class is(quality_logistic_problem)
            call problem%fdf(vec, cost, grad, D)
        class default
            THROW_HARD('invalid optimization context; logistic_cost_wrapper')
        end select
    end function logistic_cost_wrapper

    subroutine logistic_gradient_wrapper( fun_self, vec, grad, D )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: grad(D)
        real(kind=8) :: cost
        select type(problem => fun_self)
        class is(quality_logistic_problem)
            call problem%fdf(vec, cost, grad, D)
        class default
            THROW_HARD('invalid optimization context; logistic_gradient_wrapper')
        end select
    end subroutine logistic_gradient_wrapper

    subroutine logistic_fdf_wrapper( fun_self, vec, f, grad, D )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: f, grad(D)
        select type(problem => fun_self)
        class is(quality_logistic_problem)
            call problem%fdf(vec, f, grad, D)
        class default
            THROW_HARD('invalid optimization context; logistic_fdf_wrapper')
        end select
    end subroutine logistic_fdf_wrapper

    real(kind=8) function stable_sigmoid_8( eta )
        real(kind=8), intent(in) :: eta
        real(kind=8) :: z
        if( eta >= 0.0d0 )then
            z = exp(-min(eta, 80.0d0))
            stable_sigmoid_8 = 1.0d0 / (1.0d0 + z)
        else
            z = exp(max(eta, -80.0d0))
            stable_sigmoid_8 = z / (1.0d0 + z)
        endif
    end function stable_sigmoid_8

    !> the role a dataset plays in learning, from its trainable counts and its hard-rejected good rows
    integer function dataset_learn_role( ngood, nbad, nhard_good )
        integer, intent(in) :: ngood, nbad, nhard_good
        if( ngood > 0 .and. nbad > 0 )then
            dataset_learn_role = LEARN_ROLE_BALANCED
        else if( ngood > 0 )then
            dataset_learn_role = LEARN_ROLE_RECALL_ONLY
        else if( nbad > 0 .and. nhard_good == 0 )then
            dataset_learn_role = LEARN_ROLE_SPECIFICITY_ONLY
        else
            dataset_learn_role = LEARN_ROLE_SKIP
        endif
    end function dataset_learn_role

    !> per-row loss weights of a dataset's good and bad rows: each dataset weighs the same, a balanced one split
    !! between its classes
    subroutine logistic_dataset_class_weights( role, ngood, nbad, good_weight, bad_weight )
        integer,      intent(in)  :: role, ngood, nbad
        real(kind=8), intent(out) :: good_weight, bad_weight
        good_weight = 0.0d0
        bad_weight  = 0.0d0
        select case(role)
        case(LEARN_ROLE_BALANCED)
            if( ngood <= 0 .or. nbad <= 0 ) THROW_HARD('invalid balanced dataset; logistic_dataset_class_weights')
            good_weight = real(LEARN_BALANCED_GOOD_LOSS_WEIGHT, kind=8) / real(ngood, kind=8)
            bad_weight  = real(LEARN_BALANCED_BAD_LOSS_WEIGHT,  kind=8) / real(nbad,  kind=8)
        case(LEARN_ROLE_RECALL_ONLY)
            if( ngood <= 0 ) THROW_HARD('invalid recall-only dataset; logistic_dataset_class_weights')
            good_weight = 1.0d0 / real(ngood, kind=8)
        case(LEARN_ROLE_SPECIFICITY_ONLY)
            if( nbad <= 0 ) THROW_HARD('invalid specificity-only dataset; logistic_dataset_class_weights')
            bad_weight = 1.0d0 / real(nbad, kind=8)
        case default
            THROW_HARD('unsupported dataset role; logistic_dataset_class_weights')
        end select
    end subroutine logistic_dataset_class_weights

    !> the mean of the per-dataset scores blended with the mean of their lower tail
    real function robust_macro_score( scores )
        real, intent(in) :: scores(:)
        real    :: sorted(size(scores)), tmp, mean_score, tail_score
        integer :: nused, tail_n, i, j
        nused = size(scores)
        if( nused == 0 ) THROW_HARD('no scoreable datasets; robust_macro_score')
        sorted = scores
        do i = 1, nused - 1
            do j = i + 1, nused
                if( sorted(j) < sorted(i) )then
                    tmp       = sorted(i)
                    sorted(i) = sorted(j)
                    sorted(j) = tmp
                endif
            end do
        end do
        mean_score = sum(sorted(1:nused)) / real(nused)
        tail_n     = max(1, min(nused, ceiling(LEARN_ROBUST_TAIL_FRAC * real(nused))))
        tail_score = sum(sorted(1:tail_n)) / real(tail_n)
        robust_macro_score = (1.0 - LEARN_ROBUST_TAIL_WEIGHT) * mean_score + LEARN_ROBUST_TAIL_WEIGHT * tail_score
    end function robust_macro_score

end module simple_quality_logistic
