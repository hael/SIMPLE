!@descr: support-regularized 2D even/odd class-average filtering
module simple_support2D_filter
use simple_core_module_api
use simple_image, only: image
implicit none
private
public :: support2D_filter_eo
#include "simple_local_flags.inc"

integer, parameter :: LABEL_BG=0, LABEL_BND=1, LABEL_FG=2
integer, parameter :: MAX_ITS=8, ALPHA_SMOOTH_PASSES=3
real,    parameter :: PROB_FLOOR=1.e-4, SUPPORT_Z0=1.25, SUPPORT_ZWIDTH=0.50
real,    parameter :: BOUNDARY_BIAS=0.15, BETA_FRAC=0.75
real,    parameter :: DIRECT_BG_FG_COST=1.00, STEP_LABEL_COST=0.20
real,    parameter :: MIN_SIGMA=1.e-6

contains

    subroutine support2D_filter_eo( even, odd, mskdiam, lambda, support, labels, avg, n_support, n_changed )
        class(image), intent(inout) :: even, odd
        real,         intent(in)    :: mskdiam, lambda
        class(image), intent(inout) :: support, labels, avg
        integer,      intent(out)   :: n_support, n_changed
        type(image) :: even_lp, odd_lp, merged_lp
        integer :: ldim(3), ldim_odd(3), nx, ny, i, j
        integer :: n_allowed, nedge, ilab
        real    :: smpd, radius_px, guard_px, cx, cy, r2, seg_lp
        real    :: edge_mean, edge_sigma, diff_sigma, zscore, agreement
        real    :: beta, best_e, second_e, beta_sum, lambda_here
        real    :: sumd2
        real, allocatable    :: even_r(:,:,:), odd_r(:,:,:), even_lpr(:,:,:), odd_lpr(:,:,:)
        real, allocatable    :: merged_r(:,:,:), avg_r(:,:,:), support_r(:,:,:), labels_r(:,:,:)
        real, allocatable    :: prob(:,:), alpha(:,:), unary(:,:,:)
        integer, allocatable :: lab(:,:)
        logical, allocatable :: allowed(:,:)
        if( even%is_ft() .or. odd%is_ft() ) THROW_HARD('real-space images required; support2D_filter_eo')
        ldim     = even%get_ldim()
        ldim_odd = odd%get_ldim()
        if( any(ldim /= ldim_odd) ) THROW_HARD('even/odd dimensions differ; support2D_filter_eo')
        if( ldim(3) /= 1 ) THROW_HARD('2D images only; support2D_filter_eo')
        smpd = even%get_smpd()
        if( smpd <= TINY ) THROW_HARD('positive smpd required; support2D_filter_eo')
        if( mskdiam <= TINY ) THROW_HARD('positive mskdiam required; support2D_filter_eo')
        nx = ldim(1)
        ny = ldim(2)
        radius_px = 0.5 * mskdiam / smpd
        if( radius_px < 2. ) THROW_HARD('mskdiam too small for image sampling; support2D_filter_eo')
        guard_px  = 2.
        cx        = (real(nx) + 1.) / 2.
        cy        = (real(ny) + 1.) / 2.
        seg_lp    = max(4. * smpd, min(30., 0.5 * mskdiam))
        call even_lp%copy(even)
        call odd_lp%copy(odd)
        call even_lp%zero_edgeavg
        call odd_lp%zero_edgeavg
        call even_lp%bp(0., seg_lp)
        call odd_lp%bp(0., seg_lp)
        call merged_lp%copy(even_lp)
        call merged_lp%add(odd_lp)
        call merged_lp%mul(0.5)
        allocate(even_r(nx,ny,1), odd_r(nx,ny,1), even_lpr(nx,ny,1), odd_lpr(nx,ny,1))
        allocate(merged_r(nx,ny,1), avg_r(nx,ny,1), support_r(nx,ny,1), labels_r(nx,ny,1))
        allocate(prob(nx,ny), alpha(nx,ny), unary(nx,ny,3), lab(nx,ny), allowed(nx,ny))
        call even%get_rmat_sub(even_r)
        call odd%get_rmat_sub(odd_r)
        call even_lp%get_rmat_sub(even_lpr)
        call odd_lp%get_rmat_sub(odd_lpr)
        call merged_lp%get_rmat_sub(merged_r)
        n_allowed = 0
        do j = 1, ny
            do i = 1, nx
                r2 = (real(i) - cx)**2 + (real(j) - cy)**2
                allowed(i,j) = sqrt(r2) <= radius_px + guard_px
                if( allowed(i,j) ) n_allowed = n_allowed + 1
            end do
        end do
        if( n_allowed == 0 ) THROW_HARD('empty circular support prior; support2D_filter_eo')

        call edge_stats(merged_r(:,:,1), edge_mean, edge_sigma, nedge)
        sumd2 = 0.
        do j = 1, ny
            do i = 1, nx
                sumd2 = sumd2 + (even_lpr(i,j,1) - odd_lpr(i,j,1))**2
            end do
        end do
        diff_sigma = sqrt(sumd2 / real(max(1, nx * ny)))
        diff_sigma = max(diff_sigma, 0.5 * edge_sigma, MIN_SIGMA)
        do j = 1, ny
            do i = 1, nx
                zscore    = abs(merged_r(i,j,1) - edge_mean) / edge_sigma
                agreement = exp(-0.5 * min(50., ((even_lpr(i,j,1) - odd_lpr(i,j,1)) / diff_sigma)**2))
                prob(i,j) = sigmoid((zscore - SUPPORT_Z0) / SUPPORT_ZWIDTH) * sqrt(max(0., agreement))
                prob(i,j) = min(1. - PROB_FLOOR, max(PROB_FLOOR, prob(i,j)))
                if( .not. allowed(i,j) ) prob(i,j) = PROB_FLOOR
                unary(i,j,LABEL_BG+1)  = -log(max(PROB_FLOOR, 1. - prob(i,j)))
                unary(i,j,LABEL_BND+1) = -log(max(PROB_FLOOR, 4. * prob(i,j) * (1. - prob(i,j)))) + BOUNDARY_BIAS
                unary(i,j,LABEL_FG+1)  = -log(max(PROB_FLOOR, prob(i,j)))
                if( .not. allowed(i,j) )then
                    unary(i,j,LABEL_BND+1) = unary(i,j,LABEL_BND+1) + 1000.
                    unary(i,j,LABEL_FG+1)  = unary(i,j,LABEL_FG+1)  + 1000.
                endif
                lab(i,j) = min_unary_label(unary(i,j,:))
            end do
        end do
        beta_sum = 0.
        do j = 1, ny
            do i = 1, nx
                if( .not. allowed(i,j) ) cycle
                best_e   = huge(best_e)
                second_e = huge(second_e)
                do ilab = 1, 3
                    call update_best_second(unary(i,j,ilab), best_e, second_e)
                end do
                if( second_e < huge(second_e) ) beta_sum = beta_sum + max(0., second_e - best_e)
            end do
        end do
        lambda_here = max(0., lambda)
        beta = lambda_here * BETA_FRAC * beta_sum / real(max(1, n_allowed))
        call refine_labels_icm(unary, allowed, lab, beta, n_changed)
        call keep_largest_support_component(lab, n_support)
        call mark_support_boundary(lab)
        alpha = 0.
        where(lab == LABEL_FG)
            alpha = 1.
        elsewhere(lab == LABEL_BND)
            alpha = 0.5
        elsewhere
            alpha = 0.
        end where
        call smooth_alpha(alpha, allowed, ALPHA_SMOOTH_PASSES)
        even_r(:,:,1) = even_r(:,:,1) * alpha
        odd_r(:,:,1)  = odd_r(:,:,1)  * alpha
        avg_r(:,:,1)  = 0.5 * (even_r(:,:,1) + odd_r(:,:,1))
        support_r(:,:,1) = alpha
        labels_r(:,:,1)  = real(lab)
        call even%set_rmat(even_r, .false.)
        call odd%set_rmat(odd_r, .false.)
        call avg%new(ldim, smpd, .false.)
        call avg%set_rmat(avg_r, .false.)
        call support%new(ldim, smpd, .false.)
        call support%set_rmat(support_r, .false.)
        call labels%new(ldim, smpd, .false.)
        call labels%set_rmat(labels_r, .false.)
        call even_lp%kill
        call odd_lp%kill
        call merged_lp%kill
    end subroutine support2D_filter_eo

    subroutine edge_stats( img, mean, sigma, nedge )
        real,    intent(in)  :: img(:,:)
        real,    intent(out) :: mean, sigma
        integer, intent(out) :: nedge
        integer :: nx, ny, edgew, i, j
        real    :: sumv, sumv2, val
        nx = size(img, 1)
        ny = size(img, 2)
        edgew = max(2, min(8, min(nx,ny) / 8))
        sumv  = 0.
        nedge = 0
        do j = 1, ny
            do i = 1, nx
                if( i > edgew .and. i <= nx - edgew .and. j > edgew .and. j <= ny - edgew ) cycle
                sumv  = sumv + img(i,j)
                nedge = nedge + 1
            end do
        end do
        mean = sumv / real(max(1, nedge))
        sumv2 = 0.
        do j = 1, ny
            do i = 1, nx
                if( i > edgew .and. i <= nx - edgew .and. j > edgew .and. j <= ny - edgew ) cycle
                val   = img(i,j) - mean
                sumv2 = sumv2 + val * val
            end do
        end do
        sigma = sqrt(sumv2 / real(max(1, nedge - 1)))
        if( sigma < MIN_SIGMA ) sigma = MIN_SIGMA
    end subroutine edge_stats

    real function sigmoid( x )
        real, intent(in) :: x
        if( x > 40. )then
            sigmoid = 1.
        else if( x < -40. )then
            sigmoid = 0.
        else
            sigmoid = 1. / (1. + exp(-x))
        endif
    end function sigmoid

    integer function min_unary_label( u )
        real, intent(in) :: u(3)
        integer :: ilab
        min_unary_label = LABEL_BG
        do ilab = LABEL_BND, LABEL_FG
            if( u(ilab+1) < u(min_unary_label+1) ) min_unary_label = ilab
        end do
    end function min_unary_label

    subroutine update_best_second( val, best, second )
        real, intent(in)    :: val
        real, intent(inout) :: best, second
        if( val < best )then
            second = best
            best   = val
        else if( val < second )then
            second = val
        endif
    end subroutine update_best_second

    subroutine refine_labels_icm( unary, allowed, lab, beta, n_changed_total )
        real,    intent(in)    :: unary(:,:,:), beta
        logical, intent(in)    :: allowed(:,:)
        integer, intent(inout) :: lab(:,:)
        integer, intent(out)   :: n_changed_total
        integer :: nx, ny, iter, i, j, cand, cur, best_lab, n_changed
        real    :: e, best_e
        nx = size(lab, 1)
        ny = size(lab, 2)
        n_changed_total = 0
        if( beta <= TINY ) return
        do iter = 1, MAX_ITS
            n_changed = 0
            do j = 1, ny
                do i = 1, nx
                    if( .not. allowed(i,j) ) cycle
                    cur      = lab(i,j)
                    best_lab = cur
                    best_e   = label_site_energy(cur, i, j, unary, lab, beta)
                    do cand = LABEL_BG, LABEL_FG
                        if( cand == cur ) cycle
                        e = label_site_energy(cand, i, j, unary, lab, beta)
                        if( e < best_e - 1.e-6 * max(1., abs(best_e)) )then
                            best_e   = e
                            best_lab = cand
                        endif
                    end do
                    if( best_lab /= cur )then
                        lab(i,j) = best_lab
                        n_changed = n_changed + 1
                    endif
                end do
            end do
            n_changed_total = n_changed_total + n_changed
            if( n_changed == 0 ) exit
        end do
    end subroutine refine_labels_icm

    real function label_site_energy( cand, i, j, unary, lab, beta )
        integer, intent(in) :: cand, i, j
        real,    intent(in) :: unary(:,:,:), beta
        integer, intent(in) :: lab(:,:)
        label_site_energy = unary(i,j,cand+1) + beta * neighborhood_cost(cand, i, j, lab)
    end function label_site_energy

    real function neighborhood_cost( cand, i, j, lab )
        integer, intent(in) :: cand, i, j
        integer, intent(in) :: lab(:,:)
        integer :: nx, ny, di, dj, ni, nj, degree, neigh_lab
        nx = size(lab, 1)
        ny = size(lab, 2)
        neighborhood_cost = 0.
        degree = 0
        do dj = -1, 1
            do di = -1, 1
                if( di == 0 .and. dj == 0 ) cycle
                ni = i + di
                nj = j + dj
                neigh_lab = LABEL_BG
                if( ni >= 1 .and. ni <= nx .and. nj >= 1 .and. nj <= ny ) neigh_lab = lab(ni,nj)
                neighborhood_cost = neighborhood_cost + pair_cost(cand, neigh_lab)
                degree = degree + 1
            end do
        end do
        if( degree > 0 ) neighborhood_cost = neighborhood_cost / real(degree)
    end function neighborhood_cost

    real function pair_cost( ilab, jlab )
        integer, intent(in) :: ilab, jlab
        select case(abs(ilab - jlab))
            case(0)
                pair_cost = 0.
            case(1)
                pair_cost = STEP_LABEL_COST
            case default
                pair_cost = DIRECT_BG_FG_COST
        end select
    end function pair_cost

    subroutine keep_largest_support_component( lab, n_support )
        integer, intent(inout) :: lab(:,:)
        integer, intent(out)   :: n_support
        integer :: nx, ny, i, j, comp, max_comp, max_size, cur_size
        integer :: head, tail, ni, nj, di, dj, qi, qj
        integer, allocatable :: visited(:,:), queue_i(:), queue_j(:)
        nx = size(lab, 1)
        ny = size(lab, 2)
        allocate(visited(nx,ny), queue_i(nx*ny), queue_j(nx*ny))
        visited = 0
        queue_i = 0
        queue_j = 0
        comp = 0
        max_comp = 0
        max_size = 0
        do j = 1, ny
            do i = 1, nx
                if( lab(i,j) == LABEL_BG .or. visited(i,j) /= 0 ) cycle
                comp = comp + 1
                head = 1
                tail = 1
                queue_i(tail) = i
                queue_j(tail) = j
                visited(i,j)  = comp
                cur_size = 0
                do while( head <= tail )
                    qi = queue_i(head)
                    qj = queue_j(head)
                    head = head + 1
                    cur_size = cur_size + 1
                    do dj = -1, 1
                        do di = -1, 1
                            if( di == 0 .and. dj == 0 ) cycle
                            ni = qi + di
                            nj = qj + dj
                            if( ni < 1 .or. ni > nx .or. nj < 1 .or. nj > ny ) cycle
                            if( lab(ni,nj) == LABEL_BG .or. visited(ni,nj) /= 0 ) cycle
                            tail = tail + 1
                            queue_i(tail) = ni
                            queue_j(tail) = nj
                            visited(ni,nj) = comp
                        end do
                    end do
                end do
                if( cur_size > max_size )then
                    max_size = cur_size
                    max_comp = comp
                endif
            end do
        end do
        if( max_comp == 0 )then
            lab = LABEL_BG
            n_support = 0
        else
            do j = 1, ny
                do i = 1, nx
                    if( visited(i,j) /= max_comp ) lab(i,j) = LABEL_BG
                end do
            end do
            n_support = count(lab /= LABEL_BG)
        endif
    end subroutine keep_largest_support_component

    subroutine mark_support_boundary( lab )
        integer, intent(inout) :: lab(:,:)
        integer :: nx, ny, i, j, di, dj, ni, nj
        logical, allocatable :: support(:,:)
        nx = size(lab, 1)
        ny = size(lab, 2)
        allocate(support(nx,ny), source=lab /= LABEL_BG)
        do j = 1, ny
            do i = 1, nx
                if( .not. support(i,j) )then
                    lab(i,j) = LABEL_BG
                    cycle
                endif
                lab(i,j) = LABEL_FG
                do dj = -1, 1
                    do di = -1, 1
                        if( di == 0 .and. dj == 0 ) cycle
                        ni = i + di
                        nj = j + dj
                        if( ni < 1 .or. ni > nx .or. nj < 1 .or. nj > ny )then
                            lab(i,j) = LABEL_BND
                        else if( .not. support(ni,nj) )then
                            lab(i,j) = LABEL_BND
                        endif
                    end do
                end do
            end do
        end do
    end subroutine mark_support_boundary

    subroutine smooth_alpha( alpha, allowed, npass )
        real,    intent(inout) :: alpha(:,:)
        logical, intent(in)    :: allowed(:,:)
        integer, intent(in)    :: npass
        integer :: nx, ny, pass, i, j, di, dj, ni, nj, cnt
        real, allocatable :: tmp(:,:)
        nx = size(alpha, 1)
        ny = size(alpha, 2)
        allocate(tmp(nx,ny), source=0.)
        do pass = 1, npass
            tmp = 0.
            do j = 1, ny
                do i = 1, nx
                    if( .not. allowed(i,j) ) cycle
                    cnt = 0
                    do dj = -1, 1
                        do di = -1, 1
                            ni = i + di
                            nj = j + dj
                            if( ni < 1 .or. ni > nx .or. nj < 1 .or. nj > ny ) cycle
                            tmp(i,j) = tmp(i,j) + alpha(ni,nj)
                            cnt = cnt + 1
                        end do
                    end do
                    if( cnt > 0 ) tmp(i,j) = tmp(i,j) / real(cnt)
                end do
            end do
            alpha = tmp
            where(.not. allowed) alpha = 0.
        end do
        where(alpha < 0.02) alpha = 0.
        where(alpha > 0.98) alpha = 1.
    end subroutine smooth_alpha

end module simple_support2D_filter
