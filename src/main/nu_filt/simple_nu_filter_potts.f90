!@descr: simple nu filter potts implementation for volume-domain nonuniform filtering
submodule (simple_nu_filter) simple_nu_filter_potts
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine refine_nu_candidate_map_ordered_labels( candmap, n_candidates )
        integer, intent(inout) :: candmap(:,:,:)
        integer, intent(in)    :: n_candidates
        integer :: iter, color, i, j, k, imask, icand, cur_icand, best_icand, n_full(3,NU_LABEL_SMOOTH_NNEIGH), nsz
        integer :: nchanged
        integer :: n_base, n_aux
        real    :: beta, e, best_e, site_energy
        if( n_candidates < 2 ) return
        if( .not. allocated(candidate_coords) ) THROW_HARD('candidate_coords not allocated; refine_nu_candidate_map_ordered_labels')
        if( size(candidate_coords) /= n_candidates ) &
            &THROW_HARD('candidate_coords size mismatch; refine_nu_candidate_map_ordered_labels')
        n_base = size(cutoff_finds)
        n_aux  = max(0, n_candidates - n_base)
        beta = estimate_nu_label_smooth_beta(n_candidates)
        write(logfhandle,'(A,ES12.4,A,I0,A,I0,A,I0,A,I0)') '>>> NU ordered-label smoothing: beta=', beta, &
            &', max iterations=', NU_LABEL_SMOOTH_MAXITS, ', candidates=', n_candidates, &
            &', auxiliary=', n_aux, ', step tolerance=', NU_LABEL_SMOOTH_STEP_TOL
        write(logfhandle,'(A,I0,A,I0)') '>>> NU ordered-label smoothing neighborhood: ', &
            &NU_LABEL_SMOOTH_NNEIGH, '-connected, color passes=', NU_LABEL_SMOOTH_NCOLORS
        write(logfhandle,'(A,F6.3)') '>>> NU ordered-label smoothing quadratic jump fraction: ', &
            &NU_LABEL_SMOOTH_QUAD_FRAC
        call log_nu_candidate_coords
        if( beta <= TINY )then
            write(logfhandle,'(A)') '>>> NU ordered-label smoothing skipped: beta <= TINY'
            return
        endif
        site_energy = calc_nu_label_smooth_site_energy(candmap, beta)
        write(logfhandle,'(A,F12.5)') '>>> NU ordered-label smoothing initial mean site energy: ', site_energy
        do iter = 1, NU_LABEL_SMOOTH_MAXITS
            nchanged = 0
            do color = 0, NU_LABEL_SMOOTH_NCOLORS - 1
                !$omp parallel do collapse(3) schedule(static) default(shared) &
                !$omp private(i,j,k,imask,icand,cur_icand,best_icand,n_full,nsz,e,best_e) &
                !$omp reduction(+:nchanged) proc_bind(close)
                do k = 1, ldim(3)
                    do j = 1, ldim(2)
                        do i = 1, ldim(1)
                            if( .not.nu_lmask(i,j,k) ) cycle
                            if( nu_label_smooth_color(i,j,k) /= color ) cycle
                            imask = nu_mask_index(i,j,k)
                            call neigh_8_3D(ldim, [i,j,k], n_full, nsz)
                            cur_icand  = candmap(i,j,k)
                            best_icand = cur_icand
                            best_e     = dmats_mask(imask,cur_icand) + beta * &
                                &nu_label_smooth_neighborhood_cost(cur_icand, candmap, n_full, nsz)
                            do icand = 1, n_candidates
                                if( icand == cur_icand ) cycle
                                e = dmats_mask(imask,icand) + beta * &
                                    &nu_label_smooth_neighborhood_cost(icand, candmap, n_full, nsz)
                                if( nu_label_smooth_is_better(e, best_e) )then
                                    best_e     = e
                                    best_icand = icand
                                endif
                            end do
                            if( best_icand /= cur_icand )then
                                nchanged = nchanged + 1
                                candmap(i,j,k) = best_icand
                            endif
                        end do
                    end do
                end do
                !$omp end parallel do
            end do
            site_energy = calc_nu_label_smooth_site_energy(candmap, beta)
            write(logfhandle,'(A,I0,A,I0,A,F12.5)') '>>> NU ordered-label smoothing iteration ', iter, &
                &' changed voxels: ', nchanged, ', mean site energy: ', site_energy
            if( nchanged == 0 ) exit
        end do
    end subroutine refine_nu_candidate_map_ordered_labels

    module real function estimate_nu_label_smooth_beta( n_candidates )
        integer, intent(in) :: n_candidates
        integer :: i, j, k, imask, icand, nvox
        real    :: best_e, second_e, cur_e
        estimate_nu_label_smooth_beta = 0.
        nvox = 0
        if( n_candidates < 2 ) return
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.nu_lmask(i,j,k) ) cycle
                    imask = nu_mask_index(i,j,k)
                    best_e   = huge(best_e)
                    second_e = huge(second_e)
                    do icand = 1, n_candidates
                        cur_e = dmats_mask(imask,icand)
                        if( cur_e < best_e )then
                            second_e = best_e
                            best_e   = cur_e
                        else if( cur_e < second_e )then
                            second_e = cur_e
                        endif
                    end do
                    if( second_e < huge(second_e) )then
                        estimate_nu_label_smooth_beta = estimate_nu_label_smooth_beta + max(0., second_e - best_e)
                        nvox = nvox + 1
                    endif
                end do
            end do
        end do
        if( nvox > 0 ) estimate_nu_label_smooth_beta = &
            &NU_LABEL_SMOOTH_BETA_FRAC * estimate_nu_label_smooth_beta / real(nvox)
    end function estimate_nu_label_smooth_beta

    module real function nu_label_smooth_neighborhood_cost( icand, candmap, neigh, nsz )
        integer, intent(in) :: icand, candmap(:,:,:), neigh(3,NU_LABEL_SMOOTH_NNEIGH), nsz
        integer :: ineigh, ni, nj, nk, degree
        nu_label_smooth_neighborhood_cost = 0.
        degree = 0
        do ineigh = 1, nsz
            ni = neigh(1,ineigh)
            nj = neigh(2,ineigh)
            nk = neigh(3,ineigh)
            if( .not.nu_lmask(ni,nj,nk) ) cycle
            degree = degree + 1
            nu_label_smooth_neighborhood_cost = nu_label_smooth_neighborhood_cost + &
                &nu_label_smooth_pair_cost(icand, candmap(ni,nj,nk))
        end do
        if( degree > 0 ) nu_label_smooth_neighborhood_cost = nu_label_smooth_neighborhood_cost / real(degree)
    end function nu_label_smooth_neighborhood_cost

    module real function nu_label_smooth_pair_cost( icand, jcand )
        integer, intent(in) :: icand, jcand
        nu_label_smooth_pair_cost = nu_label_smooth_coord_pair_cost(candidate_coords(icand), candidate_coords(jcand))
    end function nu_label_smooth_pair_cost

    module real function nu_label_smooth_coord_pair_cost( icoord, jcoord )
        real, intent(in) :: icoord, jcoord
        real :: step_jump
        if( abs(icoord - jcoord) <= TINY )then
            nu_label_smooth_coord_pair_cost = 0.
        else
            step_jump = max(0., abs(icoord - jcoord) - real(NU_LABEL_SMOOTH_STEP_TOL))
            nu_label_smooth_coord_pair_cost = step_jump + NU_LABEL_SMOOTH_QUAD_FRAC * step_jump * step_jump
        endif
    end function nu_label_smooth_coord_pair_cost

    module integer function nu_label_smooth_color( i, j, k )
        integer, intent(in) :: i, j, k
        nu_label_smooth_color = mod(i,2) + 2 * mod(j,2) + 4 * mod(k,2)
    end function nu_label_smooth_color

    module logical function nu_label_smooth_is_better( e, best_e )
        real, intent(in) :: e, best_e
        real :: tol
        tol = NU_LABEL_SMOOTH_TIE_EPS * max(1., abs(best_e))
        nu_label_smooth_is_better = e < best_e - tol
    end function nu_label_smooth_is_better

    module real function calc_nu_label_smooth_site_energy( candmap, beta )
        integer, intent(in) :: candmap(:,:,:)
        real,    intent(in) :: beta
        integer :: i, j, k, imask, n_full(3,NU_LABEL_SMOOTH_NNEIGH), nsz, nvox
        real :: energy_sum
        calc_nu_label_smooth_site_energy = 0.
        energy_sum = 0.
        nvox = 0
        !$omp parallel do collapse(3) schedule(static) default(shared) &
        !$omp private(i,j,k,imask,n_full,nsz) reduction(+:energy_sum,nvox) proc_bind(close)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.nu_lmask(i,j,k) ) cycle
                    imask = nu_mask_index(i,j,k)
                    call neigh_8_3D(ldim, [i,j,k], n_full, nsz)
                    energy_sum = energy_sum + dmats_mask(imask,candmap(i,j,k)) + beta * &
                        &nu_label_smooth_neighborhood_cost(candmap(i,j,k), candmap, n_full, nsz)
                    nvox = nvox + 1
                end do
            end do
        end do
        !$omp end parallel do
        if( nvox > 0 ) calc_nu_label_smooth_site_energy = energy_sum / real(nvox)
    end function calc_nu_label_smooth_site_energy

end submodule simple_nu_filter_potts
