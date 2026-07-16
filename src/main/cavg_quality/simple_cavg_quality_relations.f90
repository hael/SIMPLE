!@descr: analyze-only pairwise relationship features for class-average quality
module simple_cavg_quality_relations
use simple_error,              only: simple_exception
use simple_image,              only: image
use simple_imgarr_utils,       only: pack_imgarr, dealloc_imgarr
use simple_oris,               only: oris
use simple_parameters,         only: parameters
use simple_srch_sort_loc,      only: hpsort
use simple_stat,               only: median, mad_gau
use simple_strategy2D_utils,   only: calc_cavg_pairwise_algninfo, calc_cavg_sigstats_components, &
                                    cavg_sigstats_matrices
use simple_type_defs,          only: inpl_struct
use simple_cavg_quality_feats, only: I_BP40_100_CENTER_EDGE_VAR, SIEVE_BP_CENTER_EDGE_VAR_HARD_REJECT_MIN
use simple_cavg_quality_types, only: EPS, CLIP_Z, cavg_quality_result
implicit none
private
#include "simple_local_flags.inc"

public :: test_cavg_quality_relations

integer, parameter :: N_RELATIONAL_FEATURES = 8
integer, parameter :: N_ANCHOR_REDUCERS     = 10
integer, parameter :: N_CHANNEL_REDUCERS    = 6
real,    parameter :: BP_CENTER_EDGE_POC_LOG_THRESHOLD = &
    log(SIEVE_BP_CENTER_EDGE_VAR_HARD_REJECT_MIN)

character(len=40), parameter :: RELATIONAL_FEATURE_NAMES(N_RELATIONAL_FEATURES) = [ character(len=40) :: &
    'corr_local_density',             'corr_local_cohesion', &
    'corr_to_local_medoid',           'corr_neighbourhood_margin', &
    'corr_effective_support',         'corr_population_support', &
    'corr_isolation',                 'corr_neighbour_bp_fail_fraction' ]

character(len=32), parameter :: ANCHOR_CHANNELS(N_ANCHOR_REDUCERS) = [ character(len=32) :: &
    'rotmax_cc', 'rotmax_cc', 'rotmax_cc', 'rotmax_cc', 'rotmax_cc', 'rotmax_cc', 'rotmax_cc', &
    'rotmax_cc', 'rotmax_cc', 'rotmax_cc' ]

character(len=40), parameter :: ANCHOR_REDUCERS(N_ANCHOR_REDUCERS) = [ character(len=40) :: &
    'nearest', 'topk_mean', 'topk_median', 'topk_lower_quartile', 'topk_upper_quartile', &
    'topk_spread', 'nearest_next_shell_margin', 'effective_support', 'population_support', &
    'bp_fail_fraction' ]

character(len=40), parameter :: CHANNEL_REDUCERS(N_CHANNEL_REDUCERS) = [ character(len=40) :: &
    'anchor_nearest', 'anchor_topk_mean', 'anchor_topk_median', &
    'anchor_topk_lower_quartile', 'anchor_topk_upper_quartile', 'anchor_topk_spread' ]

type, public :: cavg_quality_relation_analysis
    integer :: ncls        = 0
    integer :: nfit        = 0
    integer :: k_requested = 5
    integer :: k_stored    = 0
    real    :: hp          = 100.0
    real    :: lp          = 15.0
    real    :: trs         = 10.0
    real    :: tau         = 0.05
    real    :: signal_msk  = 0.0
    real    :: signal_oa_minmax(2) = 0.0
    logical,          allocatable :: eligible(:)
    logical,          allocatable :: fit_bp_fail(:)
    integer,          allocatable :: class_inds(:)
    integer,          allocatable :: fit_populations(:)
    integer,          allocatable :: neighbour_count(:)
    integer,          allocatable :: neighbours(:,:)
    real,             allocatable :: neighbour_corr(:,:)
    real,             allocatable :: neighbour_frc_distance(:,:)
    real,             allocatable :: neighbour_weight(:,:)
    real,             allocatable :: cc(:,:)
    real,             allocatable :: frc_distance(:,:)
    real,             allocatable :: raw(:,:)
    real,             allocatable :: features(:,:)
    type(inpl_struct), allocatable :: matches(:,:)
    type(cavg_sigstats_matrices)   :: signal
contains
    procedure :: calculate             => calculate_relation_analysis
    procedure :: write_pair_table      => write_relation_pair_table
    procedure :: write_neighbour_table => write_relation_neighbour_table
    procedure :: write_feature_table   => write_relation_feature_table
    procedure :: write_candidate_table => write_relation_candidate_table
    procedure :: kill                  => kill_relation_analysis
end type cavg_quality_relation_analysis

contains

    subroutine calculate_relation_analysis( self, imgs, cls_oris, quality, params )
        class(cavg_quality_relation_analysis), intent(inout) :: self
        class(image),                          intent(inout) :: imgs(:)
        type(oris),                            intent(in)    :: cls_oris
        type(cavg_quality_result),             intent(in)    :: quality
        type(parameters),                      intent(in)    :: params
        type(image), allocatable :: fit_imgs(:)
        integer,     allocatable :: populations(:)
        real :: mm(2)
        integer :: i, j, icls

        call self%kill()
        self%ncls        = size(imgs)
        self%k_requested = params%relational_knn
        self%hp          = params%relational_corr_hp
        self%lp          = params%relational_corr_lp
        self%trs         = params%relational_corr_trs
        self%tau         = params%relational_weight_tau
        self%signal_msk  = params%msk
        if( cls_oris%get_noris() /= self%ncls ) &
            THROW_HARD('relation analysis: class-average/orientation count mismatch')
        if( .not. allocated(quality%hard_reject) ) &
            THROW_HARD('relation analysis: missing hard-reject mask')
        if( .not. allocated(quality%raw) ) THROW_HARD('relation analysis: missing base raw features')
        if( size(quality%hard_reject) /= self%ncls ) &
            THROW_HARD('relation analysis: hard-reject size mismatch')

        allocate(self%eligible(self%ncls), source=.not. quality%hard_reject)
        self%nfit = count(self%eligible)
        self%class_inds = pack((/(i, i=1,self%ncls)/), self%eligible)
        populations = cls_oris%get_all_asint('pop')
        self%fit_populations = pack(populations, self%eligible)
        deallocate(populations)
        allocate(self%fit_bp_fail(self%nfit), source=.false.)
        do i = 1, self%nfit
            icls = self%class_inds(i)
            self%fit_bp_fail(i) = quality%raw(icls, I_BP40_100_CENTER_EDGE_VAR) < &
                                  BP_CENTER_EDGE_POC_LOG_THRESHOLD
        end do

        allocate(self%raw(self%ncls,N_RELATIONAL_FEATURES), source=0.0)
        allocate(self%features(self%ncls,N_RELATIONAL_FEATURES), source=0.0)
        allocate(self%neighbour_count(self%ncls), source=0)
        self%k_stored = min(self%k_requested, max(0, self%nfit - 1))
        allocate(self%neighbours(self%k_stored,self%ncls), source=0)
        allocate(self%neighbour_corr(self%k_stored,self%ncls), source=0.0)
        allocate(self%neighbour_frc_distance(self%k_stored,self%ncls), source=0.0)
        allocate(self%neighbour_weight(self%k_stored,self%ncls), source=0.0)

        if( self%nfit > 1 )then
            fit_imgs = pack_imgarr(imgs, self%eligible)
            self%signal_oa_minmax = [huge(1.0), -huge(1.0)]
            do i = 1, self%nfit
                mm = fit_imgs(i)%minmax(self%signal_msk)
                self%signal_oa_minmax(1) = min(self%signal_oa_minmax(1), mm(1))
                self%signal_oa_minmax(2) = max(self%signal_oa_minmax(2), mm(2))
            end do
            call calc_cavg_sigstats_components(params, fit_imgs, self%signal_oa_minmax, self%signal)
            self%matches = calc_cavg_pairwise_algninfo(params, fit_imgs, self%hp, self%lp, self%trs)
            call dealloc_imgarr(fit_imgs)
            allocate(self%cc(self%nfit,self%nfit), source=0.0)
            allocate(self%frc_distance(self%nfit,self%nfit), source=0.0)
            do i = 1, self%nfit - 1
                do j = i + 1, self%nfit
                    self%cc(i,j) = self%matches(i,j)%corr
                    self%cc(j,i) = self%cc(i,j)
                    self%frc_distance(i,j) = 1.0 / real(max(1, self%matches(i,j)%find_fsc05))
                    self%frc_distance(j,i) = self%frc_distance(i,j)
                end do
            end do
            do i = 1, self%nfit
                call calculate_wide_row(self, i)
            end do
        else if( self%nfit == 1 )then
            self%raw(self%class_inds(1),7) = 1.0
        end if
        call normalize_relation_features(self%raw, self%eligible, self%features)
    end subroutine calculate_relation_analysis

    subroutine calculate_wide_row( self, ipacked )
        class(cavg_quality_relation_analysis), intent(inout) :: self
        integer,                               intent(in)    :: ipacked
        integer, allocatable :: order(:), local_inds(:)
        real,    allocatable :: corrs(:), weights(:), pair_values(:)
        real :: sumw, support, best_score, score
        integer :: k, nnext, i, j, cnt, candidate, medoid, icls

        call sorted_cc_neighbours(self, ipacked, order, corrs)
        k    = min(self%k_requested, size(order))
        icls = self%class_inds(ipacked)
        if( k == 0 )then
            self%raw(icls,7) = 1.0
            return
        end if
        weights = exp((corrs(1:k) - corrs(1)) / self%tau)
        sumw = sum(weights)
        self%raw(icls,1) = sum(corrs(1:k)) / real(k)
        self%raw(icls,5) = log(1.0 + sumw * sumw / max(sum(weights * weights), EPS))
        support = 0.0
        do i = 1, k
            support = support + weights(i) * real(max(0, self%fit_populations(order(i))))
        end do
        self%raw(icls,6) = log(1.0 + support)
        self%raw(icls,7) = 1.0 - corrs(1)
        self%raw(icls,8) = sum(weights, mask=self%fit_bp_fail(order(1:k))) / max(sumw, EPS)
        nnext = min(k, size(order) - k)
        if( nnext > 0 ) self%raw(icls,4) = self%raw(icls,1) - &
            sum(corrs(k+1:k+nnext)) / real(nnext)

        allocate(local_inds(k+1))
        local_inds = [ipacked, order(1:k)]
        allocate(pair_values((k * (k + 1)) / 2))
        cnt = 0
        do i = 1, k
            do j = i + 1, k + 1
                cnt = cnt + 1
                pair_values(cnt) = self%cc(local_inds(i),local_inds(j))
            end do
        end do
        self%raw(icls,2) = median(pair_values)

        medoid = ipacked
        best_score = -huge(1.0)
        do i = 1, k + 1
            candidate = local_inds(i)
            score = 0.0
            do j = 1, k + 1
                if( j == i ) cycle
                score = score + self%cc(candidate,local_inds(j))
            end do
            score = score / real(k)
            if( score > best_score )then
                best_score = score
                medoid = candidate
            end if
        end do
        if( medoid == ipacked )then
            self%raw(icls,3) = self%raw(icls,1)
        else
            self%raw(icls,3) = self%cc(ipacked,medoid)
        end if

        self%neighbour_count(icls) = k
        do i = 1, k
            self%neighbours(i,icls) = self%class_inds(order(i))
            self%neighbour_corr(i,icls) = corrs(i)
            self%neighbour_frc_distance(i,icls) = self%frc_distance(ipacked,order(i))
            self%neighbour_weight(i,icls) = weights(i) / max(sumw, EPS)
        end do
        deallocate(order, corrs, weights, local_inds, pair_values)
    end subroutine calculate_wide_row

    subroutine sorted_cc_neighbours( self, ipacked, order, corrs )
        class(cavg_quality_relation_analysis), intent(in) :: self
        integer,                               intent(in) :: ipacked
        integer, allocatable,                  intent(out) :: order(:)
        real,    allocatable,                  intent(out) :: corrs(:)
        integer, allocatable :: asc_order(:)
        real,    allocatable :: asc_corrs(:)
        integer :: i, cnt, n
        n = self%nfit - 1
        allocate(asc_order(n), asc_corrs(n), order(n), corrs(n))
        cnt = 0
        do i = 1, self%nfit
            if( i == ipacked ) cycle
            cnt = cnt + 1
            asc_order(cnt) = i
            asc_corrs(cnt) = self%cc(ipacked,i)
        end do
        call hpsort(asc_corrs, asc_order)
        do i = 1, n
            order(i) = asc_order(n-i+1)
            corrs(i) = asc_corrs(n-i+1)
        end do
        deallocate(asc_order, asc_corrs)
    end subroutine sorted_cc_neighbours

    subroutine calculate_anchor_candidate_matrix( self, k_requested, values, valid_counts )
        class(cavg_quality_relation_analysis), intent(in) :: self
        integer,                               intent(in) :: k_requested
        real,    allocatable,                  intent(out) :: values(:,:)
        integer, allocatable,                  intent(out) :: valid_counts(:)
        integer, allocatable :: order(:)
        real,    allocatable :: corrs(:), weights(:)
        real :: sumw
        integer :: ipacked, icls, k, nnext
        allocate(values(self%ncls,N_ANCHOR_REDUCERS), source=0.0)
        allocate(valid_counts(self%ncls), source=0)
        do ipacked = 1, self%nfit
            icls = self%class_inds(ipacked)
            call sorted_cc_neighbours(self, ipacked, order, corrs)
            k = min(k_requested, size(order))
            valid_counts(icls) = k
            if( k == 0 )then
                deallocate(order, corrs)
                cycle
            end if
            weights = exp((corrs(1:k) - corrs(1)) / self%tau)
            sumw = sum(weights)
            values(icls,1)  = corrs(1)
            values(icls,2)  = sum(corrs(1:k)) / real(k)
            values(icls,3)  = median(corrs(1:k))
            values(icls,4)  = sample_quantile(corrs(1:k), 0.25)
            values(icls,5)  = sample_quantile(corrs(1:k), 0.75)
            values(icls,6)  = maxval(corrs(1:k)) - minval(corrs(1:k))
            nnext = min(k, size(order) - k)
            if( nnext > 0 ) values(icls,7) = values(icls,2) - &
                sum(corrs(k+1:k+nnext)) / real(nnext)
            values(icls,8) = log(1.0 + sumw * sumw / max(sum(weights * weights), EPS))
            values(icls,9) = log(1.0 + sum(weights * real(max(0, self%fit_populations(order(1:k))))))
            values(icls,10) = sum(weights, mask=self%fit_bp_fail(order(1:k))) / max(sumw, EPS)
            deallocate(order, corrs, weights)
        end do
    end subroutine calculate_anchor_candidate_matrix

    subroutine calculate_channel_candidate_matrix( self, k_requested, matrix, values, valid_counts )
        class(cavg_quality_relation_analysis), intent(in) :: self
        integer,                               intent(in) :: k_requested
        real,                                  intent(in) :: matrix(:,:)
        real,    allocatable,                  intent(out) :: values(:,:)
        integer, allocatable,                  intent(out) :: valid_counts(:)
        integer, allocatable :: order(:)
        real,    allocatable :: corrs(:), selected(:)
        integer :: ipacked, icls, k, i
        if( size(matrix,dim=1) /= self%nfit .or. size(matrix,dim=2) /= self%nfit ) &
            THROW_HARD('relation analysis: secondary channel size mismatch')
        allocate(values(self%ncls,N_CHANNEL_REDUCERS), source=0.0)
        allocate(valid_counts(self%ncls), source=0)
        do ipacked = 1, self%nfit
            icls = self%class_inds(ipacked)
            call sorted_cc_neighbours(self, ipacked, order, corrs)
            k = min(k_requested, size(order))
            valid_counts(icls) = k
            if( k == 0 )then
                deallocate(order, corrs)
                cycle
            end if
            allocate(selected(k))
            do i = 1, k
                selected(i) = matrix(ipacked,order(i))
            end do
            values(icls,1) = selected(1)
            values(icls,2) = sum(selected) / real(k)
            values(icls,3) = median(selected)
            values(icls,4) = sample_quantile(selected, 0.25)
            values(icls,5) = sample_quantile(selected, 0.75)
            values(icls,6) = maxval(selected) - minval(selected)
            deallocate(order, corrs, selected)
        end do
    end subroutine calculate_channel_candidate_matrix

    subroutine normalize_relation_features( raw, eligible, normalized )
        real,    intent(in)  :: raw(:,:)
        logical, intent(in)  :: eligible(:)
        real,    intent(out) :: normalized(:,:)
        real, allocatable :: vals(:)
        real :: med, dev
        integer :: j
        normalized = 0.0
        do j = 1, size(raw,dim=2)
            if( count(eligible) > 1 )then
                vals = pack(raw(:,j), eligible)
                med = median(vals)
                dev = mad_gau(vals, med)
                if( dev > EPS ) normalized(:,j) = max(-CLIP_Z, min(CLIP_Z, (raw(:,j) - med) / dev))
                deallocate(vals)
            end if
            where( .not. eligible ) normalized(:,j) = -CLIP_Z
        end do
    end subroutine normalize_relation_features

    real function sample_quantile( vals, fraction )
        real, intent(in) :: vals(:), fraction
        real, allocatable :: sorted(:)
        integer :: ind
        if( size(vals) == 0 )then
            sample_quantile = 0.0
            return
        end if
        sorted = vals
        call hpsort(sorted)
        ind = 1 + nint(max(0.0, min(1.0, fraction)) * real(size(sorted)-1))
        sample_quantile = sorted(ind)
        deallocate(sorted)
    end function sample_quantile

    subroutine write_relation_pair_table( self, fname )
        class(cavg_quality_relation_analysis), intent(in) :: self
        character(len=*),                      intent(in) :: fname
        integer :: funit, i, j
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        call write_relation_metadata(self, funit)
        write(funit,'(A)') 'class_i,class_j,rotmax_cc,frc05_index,frc05_distance,'//&
                           'radial_power_distance,histogram_tvd,histogram_jsd,histogram_hellinger,'//&
                           'histogram_distance,signal_stats_distance,best_rotation,best_shift_x,'//&
                           'best_shift_y,best_mirror'
        if( allocated(self%matches) )then
            do i = 1, self%nfit - 1
                do j = i + 1, self%nfit
                    write(funit,'(I0,A,I0,A,ES14.6,A,I0,A,ES14.6)',advance='no') &
                        self%class_inds(i), ',', self%class_inds(j), ',', self%cc(i,j), ',', &
                        self%matches(i,j)%find_fsc05, ',', self%frc_distance(i,j)
                    write(funit,'(A,ES14.6)',advance='no') ',', self%signal%power_distance(i,j)
                    write(funit,'(A,ES14.6)',advance='no') ',', self%signal%histogram_tvd(i,j)
                    write(funit,'(A,ES14.6)',advance='no') ',', self%signal%histogram_jsd(i,j)
                    write(funit,'(A,ES14.6)',advance='no') ',', self%signal%histogram_hellinger(i,j)
                    write(funit,'(A,ES14.6)',advance='no') ',', self%signal%histogram_distance(i,j)
                    write(funit,'(A,ES14.6)',advance='no') ',', self%signal%signal_stats_distance(i,j)
                    write(funit,'(A,ES14.6,A,ES14.6,A,ES14.6,A,L1)') ',', self%matches(i,j)%e3, ',', &
                        self%matches(i,j)%x, ',', self%matches(i,j)%y, ',', self%matches(i,j)%l_mirr
                end do
            end do
        end if
        close(funit)
    end subroutine write_relation_pair_table

    subroutine write_relation_neighbour_table( self, fname )
        class(cavg_quality_relation_analysis), intent(in) :: self
        character(len=*),                      intent(in) :: fname
        integer :: funit, icls, rank, ipacked, jpacked
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        call write_relation_metadata(self, funit)
        write(funit,'(A)') 'class,rank,neighbour,rotmax_cc,frc05_distance,radial_power_distance,'//&
                           'histogram_tvd,histogram_jsd,histogram_hellinger,histogram_distance,'//&
                           'signal_stats_distance,weight,neighbour_bp_fail'
        do icls = 1, self%ncls
            if( self%neighbour_count(icls) == 0 ) cycle
            ipacked = packed_index_for_class(self, icls)
            do rank = 1, self%neighbour_count(icls)
                jpacked = packed_index_for_class(self, self%neighbours(rank,icls))
                write(funit,'(I0,A,I0,A,I0,A,ES14.6,A,ES14.6)',advance='no') &
                    icls, ',', rank, ',', self%neighbours(rank,icls), ',', &
                    self%neighbour_corr(rank,icls), ',', self%neighbour_frc_distance(rank,icls)
                write(funit,'(A,ES14.6)',advance='no') ',', self%signal%power_distance(ipacked,jpacked)
                write(funit,'(A,ES14.6)',advance='no') ',', self%signal%histogram_tvd(ipacked,jpacked)
                write(funit,'(A,ES14.6)',advance='no') ',', self%signal%histogram_jsd(ipacked,jpacked)
                write(funit,'(A,ES14.6)',advance='no') ',', self%signal%histogram_hellinger(ipacked,jpacked)
                write(funit,'(A,ES14.6)',advance='no') ',', self%signal%histogram_distance(ipacked,jpacked)
                write(funit,'(A,ES14.6)',advance='no') ',', self%signal%signal_stats_distance(ipacked,jpacked)
                write(funit,'(A,ES14.6,A,L1)') ',', self%neighbour_weight(rank,icls), ',', &
                    self%fit_bp_fail(packed_index_for_class(self, self%neighbours(rank,icls)))
            end do
        end do
        close(funit)
    end subroutine write_relation_neighbour_table

    subroutine write_relation_feature_table( self, fname )
        class(cavg_quality_relation_analysis), intent(in) :: self
        character(len=*),                      intent(in) :: fname
        integer :: funit, icls, j
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        call write_relation_metadata(self, funit)
        write(funit,'(A)',advance='no') 'class,eligible,neighbour_count'
        do j = 1, N_RELATIONAL_FEATURES
            write(funit,'(A,A)',advance='no') ',raw_', trim(RELATIONAL_FEATURE_NAMES(j))
        end do
        do j = 1, N_RELATIONAL_FEATURES
            write(funit,'(A,A)',advance='no') ',z_', trim(RELATIONAL_FEATURE_NAMES(j))
        end do
        write(funit,*)
        do icls = 1, self%ncls
            write(funit,'(I0,A,L1,A,I0)',advance='no') icls, ',', self%eligible(icls), ',', &
                self%neighbour_count(icls)
            do j = 1, N_RELATIONAL_FEATURES
                write(funit,'(A,ES14.6)',advance='no') ',', self%raw(icls,j)
            end do
            do j = 1, N_RELATIONAL_FEATURES
                write(funit,'(A,ES14.6)',advance='no') ',', self%features(icls,j)
            end do
            write(funit,*)
        end do
        close(funit)
    end subroutine write_relation_feature_table

    subroutine write_relation_candidate_table( self, fname )
        class(cavg_quality_relation_analysis), intent(in) :: self
        character(len=*),                      intent(in) :: fname
        integer, parameter :: NK = 5
        integer :: requested_k(NK), unique_k(NK), nk_unique
        real,    allocatable :: values(:,:), normalized(:,:)
        integer, allocatable :: valid_counts(:)
        integer :: funit, ik, j, icls, k
        requested_k = [1, 3, 5, 8, self%k_requested]
        unique_k = 0
        nk_unique = 0
        do ik = 1, NK
            k = max(1, min(requested_k(ik), max(1, self%nfit-1)))
            if( nk_unique > 0 )then
                if( any(unique_k(1:nk_unique) == k) ) cycle
            end if
            nk_unique = nk_unique + 1
            unique_k(nk_unique) = k
        end do
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        call write_relation_metadata(self, funit)
        write(funit,'(A)') 'class,channel,reducer,k,tau,raw_value,normalized_value,valid_count'
        do ik = 1, nk_unique
            k = unique_k(ik)
            call calculate_anchor_candidate_matrix(self, k, values, valid_counts)
            allocate(normalized(self%ncls,N_ANCHOR_REDUCERS), source=0.0)
            call normalize_relation_features(values, self%eligible, normalized)
            do icls = 1, self%ncls
                if( .not. self%eligible(icls) ) cycle
                do j = 1, N_ANCHOR_REDUCERS
                    write(funit,'(I0,A,A,A,A,A,I0,A,ES14.6,A,ES14.6,A,ES14.6,A,I0)') &
                        icls, ',', trim(ANCHOR_CHANNELS(j)), ',', trim(ANCHOR_REDUCERS(j)), ',', &
                        k, ',', self%tau, ',', values(icls,j), ',', normalized(icls,j), ',', valid_counts(icls)
                end do
            end do
            deallocate(values, normalized, valid_counts)
            if( allocated(self%signal%power_distance) )then
                call write_channel_candidate_rows(self, funit, k, 'rotmax_frc_distance', self%frc_distance)
                call write_channel_candidate_rows(self, funit, k, 'radial_power_distance', &
                    self%signal%power_distance)
                call write_channel_candidate_rows(self, funit, k, 'histogram_tvd', self%signal%histogram_tvd)
                call write_channel_candidate_rows(self, funit, k, 'histogram_jsd', self%signal%histogram_jsd)
                call write_channel_candidate_rows(self, funit, k, 'histogram_hellinger', &
                    self%signal%histogram_hellinger)
                call write_channel_candidate_rows(self, funit, k, 'histogram_distance', &
                    self%signal%histogram_distance)
                call write_channel_candidate_rows(self, funit, k, 'signal_stats_distance', &
                    self%signal%signal_stats_distance)
            end if
        end do
        close(funit)
    end subroutine write_relation_candidate_table

    subroutine write_channel_candidate_rows( self, funit, k, channel, matrix )
        class(cavg_quality_relation_analysis), intent(in) :: self
        integer,                               intent(in) :: funit, k
        character(len=*),                      intent(in) :: channel
        real,                                  intent(in) :: matrix(:,:)
        real,    allocatable :: values(:,:), normalized(:,:)
        integer, allocatable :: valid_counts(:)
        integer :: icls, j
        call calculate_channel_candidate_matrix(self, k, matrix, values, valid_counts)
        allocate(normalized(self%ncls,N_CHANNEL_REDUCERS), source=0.0)
        call normalize_relation_features(values, self%eligible, normalized)
        do icls = 1, self%ncls
            if( .not. self%eligible(icls) ) cycle
            do j = 1, N_CHANNEL_REDUCERS
                write(funit,'(I0,A,A,A,A,A,I0,A,ES14.6,A,ES14.6,A,ES14.6,A,I0)') &
                    icls, ',', trim(channel), ',', trim(CHANNEL_REDUCERS(j)), ',', k, ',', self%tau, ',', &
                    values(icls,j), ',', normalized(icls,j), ',', valid_counts(icls)
            end do
        end do
        deallocate(values, normalized, valid_counts)
    end subroutine write_channel_candidate_rows

    subroutine write_relation_metadata( self, funit )
        class(cavg_quality_relation_analysis), intent(in) :: self
        integer,                               intent(in) :: funit
        write(funit,'(A)') '# relational_analysis_version=2'
        write(funit,'(A)') '# relational_feature_schema=corr_knn_v1'
        write(funit,'(A)') '# relation_channels=rotmax_cc,rotmax_frc_distance,radial_power_distance,'//&
                           'histogram_tvd,histogram_jsd,histogram_hellinger,histogram_distance,'//&
                           'signal_stats_distance'
        write(funit,'(A)') '# relation_anchor=rotmax_cc'
        write(funit,'(A)') '# relation_channel_semantics=rotmax_cc:similarity;all_others:distance'
        write(funit,'(A)') '# signal_stats_provider=cluster_cavgs_signal_stats_v1'
        write(funit,'(A)') '# radial_power_distance_scale=raw_euclidean'
        write(funit,'(A)') '# histogram_components_scale=native'
        write(funit,'(A)') '# histogram_distance_scale=component_minmax_composite'
        write(funit,'(A)') '# signal_stats_distance_scale=power_histogram_minmax_composite'
        write(funit,'(A,I0)') '# signal_stats_histogram_bins=', self%signal%histogram_bins
        write(funit,'(A,F10.4)') '# signal_stats_power_hp=', self%signal%power_hp
        write(funit,'(A,F10.4)') '# signal_stats_power_lp=', self%signal%power_lp
        write(funit,'(A,F10.4)') '# signal_stats_mask_radius_pixels=', self%signal_msk
        write(funit,'(A,F14.6)') '# signal_stats_oa_min=', self%signal_oa_minmax(1)
        write(funit,'(A,F14.6)') '# signal_stats_oa_max=', self%signal_oa_minmax(2)
        write(funit,'(A,F10.4)') '# relational_corr_hp=', self%hp
        write(funit,'(A,F10.4)') '# relational_corr_lp=', self%lp
        write(funit,'(A,F10.4)') '# relational_corr_trs=', self%trs
        write(funit,'(A,I0)') '# relational_knn=', self%k_requested
        write(funit,'(A,F10.6)') '# relational_weight_tau=', self%tau
        write(funit,'(A,F10.4)') '# bp_center_edge_threshold=', &
            SIEVE_BP_CENTER_EDGE_VAR_HARD_REJECT_MIN
        write(funit,'(A,I0)') '# n_classes=', self%ncls
        write(funit,'(A,I0)') '# n_eligible=', self%nfit
    end subroutine write_relation_metadata

    integer function packed_index_for_class( self, icls )
        class(cavg_quality_relation_analysis), intent(in) :: self
        integer,                               intent(in) :: icls
        integer :: i
        packed_index_for_class = 0
        do i = 1, self%nfit
            if( self%class_inds(i) == icls )then
                packed_index_for_class = i
                return
            end if
        end do
        if( packed_index_for_class == 0 ) THROW_HARD('relation analysis: class is not eligible')
    end function packed_index_for_class

    subroutine kill_relation_analysis( self )
        class(cavg_quality_relation_analysis), intent(inout) :: self
        if( allocated(self%eligible)                       ) deallocate(self%eligible)
        if( allocated(self%fit_bp_fail)                    ) deallocate(self%fit_bp_fail)
        if( allocated(self%class_inds)                     ) deallocate(self%class_inds)
        if( allocated(self%fit_populations)                ) deallocate(self%fit_populations)
        if( allocated(self%neighbour_count)                ) deallocate(self%neighbour_count)
        if( allocated(self%neighbours)                     ) deallocate(self%neighbours)
        if( allocated(self%neighbour_corr)                 ) deallocate(self%neighbour_corr)
        if( allocated(self%neighbour_frc_distance)         ) deallocate(self%neighbour_frc_distance)
        if( allocated(self%neighbour_weight)               ) deallocate(self%neighbour_weight)
        if( allocated(self%cc)                             ) deallocate(self%cc)
        if( allocated(self%frc_distance)                   ) deallocate(self%frc_distance)
        if( allocated(self%raw)                            ) deallocate(self%raw)
        if( allocated(self%features)                       ) deallocate(self%features)
        if( allocated(self%matches)                        ) deallocate(self%matches)
        call self%signal%kill()
        self%ncls = 0
        self%nfit = 0
        self%k_stored = 0
        self%signal_msk = 0.0
        self%signal_oa_minmax = 0.0
    end subroutine kill_relation_analysis

    subroutine test_cavg_quality_relations()
        use simple_test_utils, only: assert_int, assert_real, assert_true
        type(cavg_quality_relation_analysis) :: analysis
        real,    allocatable :: values(:,:)
        integer, allocatable :: valid_counts(:), order(:)
        real,    allocatable :: corrs(:)
        character(len=1024) :: line
        integer :: funit, ios
        logical :: found_version, found_pair_channel, found_candidate_channel
        real, parameter :: TOL = 1.0e-5

        analysis%ncls = 4
        analysis%nfit = 4
        analysis%tau  = 0.1
        analysis%class_inds     = [1, 2, 3, 4]
        analysis%fit_populations = [10, 20, 30, 40]
        analysis%fit_bp_fail     = [.false., .true., .false., .true.]
        allocate(analysis%eligible(4), source=.true.)
        allocate(analysis%cc(4,4), source=0.0)
        allocate(analysis%frc_distance(4,4), source=0.0)
        analysis%cc(1,2) = 0.9
        analysis%cc(1,3) = 0.7
        analysis%cc(1,4) = 0.2
        analysis%cc(2,3) = 0.8
        analysis%cc(2,4) = 0.1
        analysis%cc(3,4) = 0.4
        analysis%cc = analysis%cc + transpose(analysis%cc)
        analysis%frc_distance(1,2) = 0.10
        analysis%frc_distance(1,3) = 0.20
        analysis%frc_distance(1,4) = 0.50
        analysis%frc_distance = analysis%frc_distance + transpose(analysis%frc_distance)

        call sorted_cc_neighbours(analysis, 1, order, corrs)
        call assert_int(2, order(1), 'relational neighbours: highest correlation first')
        call assert_int(3, order(2), 'relational neighbours: second correlation')
        call assert_int(4, order(3), 'relational neighbours: lowest correlation last')
        deallocate(order, corrs)

        call calculate_anchor_candidate_matrix(analysis, 2, values, valid_counts)
        call assert_int(2, valid_counts(1), 'relational reducers: requested top-k count')
        call assert_real(0.9, values(1,1), TOL, 'relational reducers: nearest correlation')
        call assert_real(0.8, values(1,2), TOL, 'relational reducers: top-k mean')
        call assert_real(0.7, values(1,4), TOL, 'relational reducers: lower quartile')
        call assert_real(0.9, values(1,5), TOL, 'relational reducers: upper quartile')
        call assert_real(0.6, values(1,7), TOL, 'relational reducers: next-shell margin')
        call assert_real(1.0 / (1.0 + exp(-2.0)), values(1,10), TOL, &
            'relational reducers: weighted neighbour fail fraction')
        deallocate(values, valid_counts)

        call calculate_channel_candidate_matrix(analysis, 2, analysis%frc_distance, values, valid_counts)
        call assert_int(2, valid_counts(1), 'secondary reducers: anchor top-k count')
        call assert_real(0.10, values(1,1), TOL, 'secondary reducers: anchor-nearest value')
        call assert_real(0.15, values(1,2), TOL, 'secondary reducers: anchor top-k mean')
        call assert_real(0.15, values(1,3), TOL, 'secondary reducers: anchor top-k median')
        call assert_real(0.10, values(1,4), TOL, 'secondary reducers: lower quartile')
        call assert_real(0.20, values(1,5), TOL, 'secondary reducers: upper quartile')
        call assert_real(0.10, values(1,6), TOL, 'secondary reducers: spread')
        deallocate(values, valid_counts)

        analysis%signal%power_distance        = analysis%frc_distance
        analysis%signal%histogram_tvd         = 2.0 * analysis%frc_distance
        analysis%signal%histogram_jsd         = 3.0 * analysis%frc_distance
        analysis%signal%histogram_hellinger   = 4.0 * analysis%frc_distance
        analysis%signal%histogram_distance    = 5.0 * analysis%frc_distance
        analysis%signal%signal_stats_distance = 6.0 * analysis%frc_distance
        call analysis%write_pair_table('/tmp/simple_test_relational_pairs.tmp')
        found_version = .false.
        found_pair_channel = .false.
        open(newunit=funit, file='/tmp/simple_test_relational_pairs.tmp', status='old', action='read')
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( index(line, 'relational_analysis_version=2') > 0 ) found_version = .true.
            if( index(line, 'radial_power_distance') > 0 .and. &
                index(line, 'histogram_hellinger') > 0 ) found_pair_channel = .true.
        end do
        close(funit, status='delete')
        call assert_true(found_version, 'relational output: analysis schema version')
        call assert_true(found_pair_channel, 'relational output: signal channel columns')

        call analysis%write_candidate_table('/tmp/simple_test_relational_candidates.tmp')
        found_candidate_channel = .false.
        open(newunit=funit, file='/tmp/simple_test_relational_candidates.tmp', status='old', action='read')
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( index(line, 'histogram_tvd,anchor_topk_mean') > 0 ) found_candidate_channel = .true.
        end do
        close(funit, status='delete')
        call assert_true(found_candidate_channel, 'relational candidates: signal channel reducer')
        call analysis%kill()
        call assert_true(.not. allocated(analysis%signal%power_distance), &
            'relational signal matrices: kill releases components')
    end subroutine test_cavg_quality_relations

end module simple_cavg_quality_relations
