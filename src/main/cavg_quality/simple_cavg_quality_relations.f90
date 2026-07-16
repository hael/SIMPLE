!@descr: one promoted pairwise-neighbour feature for class-average quality
module simple_cavg_quality_relations
use simple_error,              only: simple_exception
use simple_image,              only: image
use simple_imgarr_utils,       only: pack_imgarr, dealloc_imgarr
use simple_parameters,         only: parameters
use simple_srch_sort_loc,      only: hpsort
use simple_stat,               only: median, mad_gau
use simple_strategy2D_utils,   only: calc_cavg_pairwise_algninfo, calc_cavg_sigstats_components, &
                                    cavg_sigstats_matrices
use simple_type_defs,          only: inpl_struct
use simple_cavg_quality_types, only: EPS, CLIP_Z, CAVG_RELATIONAL_SCHEMA_CORR_KNN_SIGNAL_V1, &
                                    cavg_quality_result
implicit none
private
#include "simple_local_flags.inc"

public :: CAVG_RELATIONAL_FEATURE_NAME
public :: cavg_quality_relation_analysis
public :: test_cavg_quality_relations

character(len=*), parameter :: CAVG_RELATIONAL_FEATURE_NAME = 'signal_stats_anchor_topk_mean'

type, public :: cavg_quality_relation_analysis
    integer :: ncls        = 0
    integer :: nfit        = 0
    integer :: k_requested = 0
    real    :: hp          = 0.0
    real    :: lp          = 0.0
    real    :: trs         = 0.0
    real    :: signal_msk  = 0.0
    real    :: signal_oa_minmax(2) = 0.0
    logical, allocatable :: eligible(:)
    integer, allocatable :: class_inds(:)
    real,    allocatable :: raw(:)
    real,    allocatable :: feature(:)
    type(cavg_sigstats_matrices) :: signal
contains
    procedure :: calculate => calculate_relation_analysis
    procedure :: kill      => kill_relation_analysis
end type cavg_quality_relation_analysis

contains

    subroutine calculate_relation_analysis( self, imgs, quality, params, k_requested, hp, lp, trs )
        class(cavg_quality_relation_analysis), intent(inout) :: self
        class(image),                          intent(inout) :: imgs(:)
        type(cavg_quality_result),             intent(in)    :: quality
        type(parameters),                      intent(in)    :: params
        integer,                               intent(in)    :: k_requested
        real,                                  intent(in)    :: hp, lp, trs
        type(image),       allocatable :: fit_imgs(:)
        type(inpl_struct), allocatable :: matches(:,:)
        type(parameters) :: pair_params
        real,              allocatable :: cc(:,:), distance(:,:)
        real :: mm(2)
        integer :: i, j

        call self%kill()
        self%ncls        = size(imgs)
        self%k_requested = k_requested
        self%hp          = hp
        self%lp          = lp
        self%trs         = trs
        self%signal_msk  = params%msk
        if( k_requested < 1 ) THROW_HARD('relation analysis: neighbour count must be positive')
        if( hp <= 0.0 .or. lp <= 0.0 .or. hp < lp ) &
            THROW_HARD('relation analysis: invalid correlation resolution limits')
        if( trs < 0.0 ) THROW_HARD('relation analysis: shift range must be nonnegative')
        if( .not. allocated(quality%hard_reject) ) THROW_HARD('relation analysis: missing hard-reject mask')
        if( size(quality%hard_reject) /= self%ncls ) THROW_HARD('relation analysis: hard-reject size mismatch')

        allocate(self%eligible(self%ncls), source=.not. quality%hard_reject)
        self%nfit = count(self%eligible)
        self%class_inds = pack((/(i, i=1,self%ncls)/), self%eligible)
        allocate(self%raw(self%ncls), source=0.0)
        allocate(self%feature(self%ncls), source=-CLIP_Z)
        if( self%nfit <= 1 )then
            if( self%nfit == 1 ) self%feature(self%class_inds(1)) = 0.0
            return
        end if

        fit_imgs = pack_imgarr(imgs, self%eligible)
        self%signal_oa_minmax = [huge(1.0), -huge(1.0)]
        do i = 1, self%nfit
            mm = fit_imgs(i)%minmax(self%signal_msk)
            self%signal_oa_minmax(1) = min(self%signal_oa_minmax(1), mm(1))
            self%signal_oa_minmax(2) = max(self%signal_oa_minmax(2), mm(2))
        end do
        call calc_cavg_sigstats_components(params, fit_imgs, self%signal_oa_minmax, self%signal)
        pair_params = params
        pair_params%ctf    = 'no'
        pair_params%objfun = 'cc'
        matches = calc_cavg_pairwise_algninfo(pair_params, fit_imgs, hp, lp, trs)
        allocate(cc(self%nfit,self%nfit), source=0.0)
        do i = 1, self%nfit - 1
            do j = i + 1, self%nfit
                cc(i,j) = matches(i,j)%corr
                cc(j,i) = cc(i,j)
            end do
        end do
        distance = self%signal%signal_stats_distance
        call calculate_promoted_feature(cc, distance, self%class_inds, self%k_requested, self%raw)
        call normalize_relation_feature(self%raw, self%eligible, self%feature)
        deallocate(cc, distance, matches)
        call dealloc_imgarr(fit_imgs)
    end subroutine calculate_relation_analysis

    subroutine calculate_promoted_feature( cc, distance, class_inds, k_requested, raw )
        real,    intent(in)    :: cc(:,:), distance(:,:)
        integer, intent(in)    :: class_inds(:), k_requested
        real,    intent(inout) :: raw(:)
        integer, allocatable :: order(:)
        real,    allocatable :: corrs(:)
        integer :: ipacked, icls, k
        if( size(cc,dim=1) /= size(cc,dim=2) .or. any(shape(distance) /= shape(cc)) ) &
            THROW_HARD('calculate_promoted_feature: pairwise matrix size mismatch')
        if( size(class_inds) /= size(cc,dim=1) ) THROW_HARD('calculate_promoted_feature: class map size mismatch')
        do ipacked = 1, size(class_inds)
            call sorted_cc_neighbours(cc, ipacked, order, corrs)
            k = min(k_requested, size(order))
            icls = class_inds(ipacked)
            if( k > 0 ) raw(icls) = sum(distance(ipacked,order(1:k))) / real(k)
            deallocate(order, corrs)
        end do
    end subroutine calculate_promoted_feature

    subroutine sorted_cc_neighbours( cc, ipacked, order, corrs )
        real,    intent(in)  :: cc(:,:)
        integer, intent(in)  :: ipacked
        integer, allocatable, intent(out) :: order(:)
        real,    allocatable, intent(out) :: corrs(:)
        integer, allocatable :: asc_order(:)
        real,    allocatable :: asc_corrs(:)
        integer :: i, cnt, n
        n = size(cc,dim=1) - 1
        allocate(asc_order(n), asc_corrs(n), order(n), corrs(n))
        cnt = 0
        do i = 1, size(cc,dim=1)
            if( i == ipacked ) cycle
            cnt = cnt + 1
            asc_order(cnt) = i
            asc_corrs(cnt) = cc(ipacked,i)
        end do
        call hpsort(asc_corrs, asc_order)
        do i = 1, n
            order(i) = asc_order(n-i+1)
            corrs(i) = asc_corrs(n-i+1)
        end do
        deallocate(asc_order, asc_corrs)
    end subroutine sorted_cc_neighbours

    subroutine normalize_relation_feature( raw, eligible, normalized )
        real,    intent(in)  :: raw(:)
        logical, intent(in)  :: eligible(:)
        real,    intent(out) :: normalized(:)
        real, allocatable :: vals(:)
        real :: med, dev
        normalized = -CLIP_Z
        if( count(eligible) == 0 ) return
        if( count(eligible) == 1 )then
            where( eligible ) normalized = 0.0
            return
        end if
        vals = pack(raw, eligible)
        med = median(vals)
        dev = mad_gau(vals, med)
        where( eligible )
            normalized = 0.0
        end where
        if( dev > EPS )then
            where( eligible ) normalized = max(-CLIP_Z, min(CLIP_Z, (raw - med) / dev))
        end if
        deallocate(vals)
    end subroutine normalize_relation_feature

    subroutine kill_relation_analysis( self )
        class(cavg_quality_relation_analysis), intent(inout) :: self
        if( allocated(self%eligible)   ) deallocate(self%eligible)
        if( allocated(self%class_inds) ) deallocate(self%class_inds)
        if( allocated(self%raw)        ) deallocate(self%raw)
        if( allocated(self%feature)    ) deallocate(self%feature)
        call self%signal%kill()
        self%ncls = 0
        self%nfit = 0
        self%k_requested = 0
        self%hp = 0.0
        self%lp = 0.0
        self%trs = 0.0
        self%signal_msk = 0.0
        self%signal_oa_minmax = 0.0
    end subroutine kill_relation_analysis

    subroutine test_cavg_quality_relations()
        use simple_test_utils, only: assert_real
        real :: cc(4,4), distance(4,4), raw(4)
        integer :: class_inds(4)
        real, parameter :: TOL = 1.0e-6
        cc = 0.0
        distance = 0.0
        raw = 0.0
        class_inds = [1, 2, 3, 4]
        cc(1,2) = 0.9
        cc(1,3) = 0.7
        cc(1,4) = 0.2
        cc(2,3) = 0.8
        cc(2,4) = 0.1
        cc(3,4) = 0.4
        cc = cc + transpose(cc)
        distance(1,2) = 0.10
        distance(1,3) = 0.20
        distance(1,4) = 0.90
        distance(2,3) = 0.30
        distance(2,4) = 0.80
        distance(3,4) = 0.70
        distance = distance + transpose(distance)
        call calculate_promoted_feature(cc, distance, class_inds, 2, raw)
        call assert_real(0.15, raw(1), TOL, 'relational feature: CC-anchor top-k mean')
        call assert_real(0.20, raw(2), TOL, 'relational feature: per-class neighbour ordering')
        if( trim(CAVG_RELATIONAL_SCHEMA_CORR_KNN_SIGNAL_V1) == '' ) &
            THROW_HARD('relational feature schema must be named')
    end subroutine test_cavg_quality_relations

end module simple_cavg_quality_relations
