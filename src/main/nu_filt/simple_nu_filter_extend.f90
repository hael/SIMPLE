!@descr: simple nu filter extend implementation for volume-domain nonuniform filtering
submodule (simple_nu_filter) simple_nu_filter_extend
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine extend_nu_filter_highres( vol_even, vol_odd, threshold_pct, accept_pct, new_limit, stats )
        class(image), intent(in) :: vol_even, vol_odd
        real,         intent(in) :: threshold_pct   ! e.g. 10.0
        real,         intent(in) :: accept_pct      ! minimum selected frontier percentage
        real,         intent(in) :: new_limit        ! Angstrom limit for the proposed shell
        type(nu_highres_extension_stats), optional, intent(out) :: stats
        type(image)       :: vol_even_filt_new, vol_odd_filt_new
        type(string)      :: even_cache_fname, odd_cache_fname
        real, allocatable :: dmat_new(:,:,:), dmat_tmp(:,:,:), dmat_finest(:,:,:)
        integer, allocatable :: cutoff_finds_new(:)
        integer           :: new_find, n_finest, n_total, n_extended, sz_old
        real              :: pct_finest, x
        logical, allocatable :: extend_mask(:,:,:), extend_to_new(:,:,:)
        integer, allocatable :: extend_choice(:,:,:)
        type(nu_highres_extension_stats) :: local_stats
        integer           :: i, j, k
        logical           :: l_use_aux_extension
        local_stats%new_limit = new_limit
        if( .not.allocated(filtmap)      ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds first')
        if( .not.allocated(srcmap)       ) THROW_HARD('srcmap not allocated; run optimize_nu_cutoff_finds first')
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated')
        if( .not.allocated(nu_lmask)     ) THROW_HARD('nu_lmask not allocated; run setup_nu_dmats first')
        sz_old    = size(cutoff_finds)
        n_total   = count(nu_lmask)
        local_stats%n_mask = n_total
        if( n_total == 0 )then
            if( present(stats) ) stats = local_stats
            return
        endif
        n_finest  = count(nu_lmask .and. filtmap == sz_old)
        local_stats%n_tested = n_finest
        pct_finest = 100. * real(n_finest) / real(n_total)
        local_stats%pct_tested_mask = pct_finest
        if( n_finest == 0 )then
            if( present(stats) ) stats = local_stats
            return
        endif
        if( pct_finest < threshold_pct )then
            if( present(stats) ) stats = local_stats
            return   ! trigger not met, nothing to do
        endif
        new_find = calc_fourier_index(new_limit, box, smpd)
        if( new_find <= cutoff_finds(sz_old) )then
            if( present(stats) ) stats = local_stats
            return
        endif
        if( any(cutoff_finds == new_find) )then
            if( present(stats) ) stats = local_stats
            return
        endif
        local_stats%attempted = .true.
        ! --- build the new Butterworth filter and cache filtered vols ---
        even_cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), new_find)
        odd_cache_fname  = filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  new_find)
        if( .not.file_exists(even_cache_fname) .or. .not.file_exists(odd_cache_fname) ) then
            ! generate just this one new filtered pair — cheap, one FFT+filter+IFFT pass
            call generate_single_filtered_pair(vol_even, vol_odd, new_find, even_cache_fname, odd_cache_fname)
        end if
        ! --- build the extend mask: voxels currently at the finest limit ---
        allocate(extend_mask(ldim(1),ldim(2),ldim(3)), source=.false.)
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k)
        do k = 1, ldim(3)
        do j = 1, ldim(2)
            do i = 1, ldim(1)
                extend_mask(i,j,k) = (nu_lmask(i,j,k) .and. filtmap(i,j,k) == sz_old)
            end do
        end do
        end do
        !$omp end parallel do
        ! --- evaluate the new objective only within the mask ---
        allocate(dmat_new(ldim(1),ldim(2),ldim(3)), source=huge(x))
        allocate(dmat_tmp(ldim(1),ldim(2),ldim(3)), source=0.)
        call vol_even_filt_new%new(ldim, smpd)
        call vol_odd_filt_new%new(ldim, smpd)
        call vol_even_filt_new%read(even_cache_fname)
        call vol_odd_filt_new%read(odd_cache_fname)
        call vol_even%nu_objective(vol_even_filt_new, vol_odd, vol_odd_filt_new, dmat_new, nu_lmask)
        call smooth_nu_objective(dmat_new, dmat_tmp)
        ! dmat_tmp is only a work buffer; smooth_nu_objective updates dmat_new in place.
        allocate(dmat_finest(ldim(1),ldim(2),ldim(3)), source=huge(x))
        if( allocated(dmat_finest_cached) ) then
            if( all(shape(dmat_finest_cached) == ldim) ) then
                dmat_finest = dmat_finest_cached
            else
                call vol_even_filt_new%read(filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(sz_old)))
                call vol_odd_filt_new%read(filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  cutoff_finds(sz_old)))
                call vol_even%nu_objective(vol_even_filt_new, vol_odd, vol_odd_filt_new, dmat_finest, nu_lmask)
                call smooth_nu_objective(dmat_finest, dmat_tmp)
            end if
        else
            call vol_even_filt_new%read(filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(sz_old)))
            call vol_odd_filt_new%read(filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  cutoff_finds(sz_old)))
            call vol_even%nu_objective(vol_even_filt_new, vol_odd, vol_odd_filt_new, dmat_finest, nu_lmask)
            call smooth_nu_objective(dmat_finest, dmat_tmp)
        end if
        ! --- update filtmap in place for the masked voxels ---
        l_use_aux_extension = l_aux_source_unordered_potts .and. allocated(dmats_aux_mask)
        n_extended = 0
        if( l_use_aux_extension )then
            allocate(extend_choice(ldim(1),ldim(2),ldim(3)), source=0)
            call init_nu_highres_extension_selection_aux(extend_mask, dmat_finest, dmat_new, &
                &extend_choice, n_extended)
            call refine_nu_highres_extension_selection_aux(extend_mask, dmat_finest, dmat_new, &
                &extend_choice, sz_old, n_extended)
        else
            allocate(extend_to_new(ldim(1),ldim(2),ldim(3)), source=.false.)
            call init_nu_highres_extension_selection(extend_mask, dmat_finest, dmat_new, extend_to_new, n_extended)
            call refine_nu_highres_extension_selection(extend_mask, dmat_finest, dmat_new, &
                &extend_to_new, sz_old, n_extended)
        endif
        local_stats%n_extended = n_extended
        if( n_finest > 0 ) local_stats%pct_extended_tested = 100. * real(n_extended) / real(n_finest)
        local_stats%applied = n_extended > 0 .and. local_stats%pct_extended_tested >= accept_pct
        local_stats%promote_next = local_stats%applied
        if( .not. local_stats%applied ) then
            call delete_cached_filtered_pair(new_find)
            call vol_even_filt_new%kill
            call vol_odd_filt_new%kill
            if( allocated(extend_to_new) ) deallocate(extend_to_new)
            if( allocated(extend_choice) ) deallocate(extend_choice)
            deallocate(extend_mask, dmat_new, dmat_finest, dmat_tmp)
            if( present(stats) ) stats = local_stats
            return
        end if
        if( l_use_aux_extension )then
            call apply_nu_highres_extension_selection_aux(extend_choice, sz_old, sz_old + 1)
        else
            call apply_nu_highres_extension_selection(extend_to_new, sz_old + 1)
        endif
        ! --- grow cutoff_finds to include the new level ---
        allocate(cutoff_finds_new(sz_old + 1))
        cutoff_finds_new(:sz_old)  = cutoff_finds
        cutoff_finds_new(sz_old+1) = new_find
        call move_alloc(cutoff_finds_new, cutoff_finds)
        call append_nu_highres_candidate_coord(sz_old, real(sz_old + 1))
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        if( l_use_aux_extension )then
            call cache_nu_highres_extension_frontier_after_aux_selection(dmat_new, sz_old + 1)
        else
            allocate(dmat_finest_cached(ldim(1),ldim(2),ldim(3)), source=dmat_new)
        endif
        write(logfhandle,'(A,I12,A,F8.2,A)') '>>> Extended ', n_extended, ' voxels to ', new_limit, ' A'
        call vol_even_filt_new%kill
        call vol_odd_filt_new%kill
        if( allocated(extend_to_new) ) deallocate(extend_to_new)
        if( allocated(extend_choice) ) deallocate(extend_choice)
        deallocate(extend_mask, dmat_new, dmat_finest, dmat_tmp)
        if( present(stats) ) stats = local_stats
    end subroutine extend_nu_filter_highres

    module subroutine extend_nu_filter_highres_shell_next( vol_even, vol_odd, stats, accept_pct )
        class(image), intent(in) :: vol_even, vol_odd
        type(nu_highres_extension_stats), optional, intent(out) :: stats
        real, optional, intent(in) :: accept_pct
        type(nu_highres_extension_stats) :: local_stats
        integer :: new_find, sz_old
        real    :: new_limit, accept_threshold
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; extend_nu_filter_highres_shell_next')
        accept_threshold = NU_REFINE_EXTENSION_ACCEPT_PCT
        if( present(accept_pct) ) accept_threshold = max(0., accept_pct)
        sz_old   = size(cutoff_finds)
        new_find = cutoff_finds(sz_old) + 1
        do while( new_find <= box/2 )
            if( .not.any(cutoff_finds == new_find) ) exit
            new_find = new_find + 1
        end do
        if( new_find > box/2 )then
            if( present(stats) ) stats = local_stats
            return
        endif
        new_limit = calc_lowpass_lim(new_find, box, smpd)
        call extend_nu_filter_highres(vol_even, vol_odd, NU_HIGHRES_EXTENSION_THRESHOLD_PCT, &
            &accept_threshold, new_limit, stats=local_stats)
        if( present(stats) ) stats = local_stats
    end subroutine extend_nu_filter_highres_shell_next

    module subroutine extend_nu_filter_highres_shells( vol_even, vol_odd, nsteps, accept_pct )
        class(image), intent(in) :: vol_even, vol_odd
        integer, optional, intent(out) :: nsteps
        real, optional, intent(in) :: accept_pct
        type(nu_highres_extension_stats) :: step_stats
        integer :: nsteps_local
        real :: accept_threshold
        nsteps_local = 0
        accept_threshold = NU_POSTPROCESS_EXTENSION_ACCEPT_PCT
        if( present(accept_pct) ) accept_threshold = max(0., accept_pct)
        do
            call extend_nu_filter_highres_shell_next(vol_even, vol_odd, stats=step_stats, &
                &accept_pct=accept_threshold)
            if( .not. step_stats%attempted ) exit
            if( .not. step_stats%applied   ) exit
            nsteps_local = nsteps_local + 1
            if( .not. step_stats%promote_next ) exit
        end do
        if( present(nsteps) ) nsteps = nsteps_local
    end subroutine extend_nu_filter_highres_shells

    module subroutine init_nu_highres_extension_selection( extend_mask, dmat_old, dmat_new, extend_to_new, n_extended )
        logical, intent(in)    :: extend_mask(:,:,:)
        real,    intent(in)    :: dmat_old(:,:,:), dmat_new(:,:,:)
        logical, intent(inout) :: extend_to_new(:,:,:)
        integer, intent(out)   :: n_extended
        integer :: i, j, k
        extend_to_new = .false.
        n_extended = 0
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) reduction(+:n_extended)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.extend_mask(i,j,k) ) cycle
                    if( dmat_new(i,j,k) < dmat_old(i,j,k) )then
                        extend_to_new(i,j,k) = .true.
                        n_extended = n_extended + 1
                    endif
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine init_nu_highres_extension_selection

    module subroutine refine_nu_highres_extension_selection( extend_mask, dmat_old, dmat_new, extend_to_new, old_label, n_extended )
        logical, intent(in)    :: extend_mask(:,:,:)
        real,    intent(in)    :: dmat_old(:,:,:), dmat_new(:,:,:)
        logical, intent(inout) :: extend_to_new(:,:,:)
        integer, intent(in)    :: old_label
        integer, intent(out)   :: n_extended
        integer :: iter, color, i, j, k, n_full(3,NU_LABEL_SMOOTH_NNEIGH), nsz, nchanged
        real    :: beta, old_coord, new_coord, e_old, e_new
        beta = estimate_nu_highres_extension_beta(extend_mask, dmat_old, dmat_new)
        old_coord = nu_candidate_coord_for_label(old_label)
        new_coord = real(old_label + 1)
        if( beta <= TINY )then
            n_extended = count(extend_to_new)
            return
        endif
        do iter = 1, NU_LABEL_SMOOTH_MAXITS
            nchanged = 0
            do color = 0, NU_LABEL_SMOOTH_NCOLORS - 1
                !$omp parallel do collapse(3) schedule(static) default(shared) &
                !$omp private(i,j,k,n_full,nsz,e_old,e_new) reduction(+:nchanged) proc_bind(close)
                do k = 1, ldim(3)
                    do j = 1, ldim(2)
                        do i = 1, ldim(1)
                            if( .not.extend_mask(i,j,k) ) cycle
                            if( nu_label_smooth_color(i,j,k) /= color ) cycle
                            call neigh_8_3D(ldim, [i,j,k], n_full, nsz)
                            e_old = dmat_old(i,j,k) + beta * &
                                &nu_highres_extension_neighborhood_cost(old_coord, extend_mask, &
                                &extend_to_new, old_label, new_coord, n_full, nsz)
                            e_new = dmat_new(i,j,k) + beta * &
                                &nu_highres_extension_neighborhood_cost(new_coord, extend_mask, &
                                &extend_to_new, old_label, new_coord, n_full, nsz)
                            if( extend_to_new(i,j,k) )then
                                if( nu_label_smooth_is_better(e_old, e_new) )then
                                    extend_to_new(i,j,k) = .false.
                                    nchanged = nchanged + 1
                                endif
                            else
                                if( nu_label_smooth_is_better(e_new, e_old) )then
                                    extend_to_new(i,j,k) = .true.
                                    nchanged = nchanged + 1
                                endif
                            endif
                        end do
                    end do
                end do
                !$omp end parallel do
            end do
            if( nchanged == 0 ) exit
        end do
        n_extended = count(extend_to_new)
    end subroutine refine_nu_highres_extension_selection

    module real function estimate_nu_highres_extension_beta( extend_mask, dmat_old, dmat_new )
        logical, intent(in) :: extend_mask(:,:,:)
        real,    intent(in) :: dmat_old(:,:,:), dmat_new(:,:,:)
        integer :: i, j, k, nvox
        estimate_nu_highres_extension_beta = 0.
        nvox = 0
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.extend_mask(i,j,k) ) cycle
                    estimate_nu_highres_extension_beta = estimate_nu_highres_extension_beta + &
                        &abs(dmat_new(i,j,k) - dmat_old(i,j,k))
                    nvox = nvox + 1
                end do
            end do
        end do
        if( nvox > 0 ) estimate_nu_highres_extension_beta = &
            &NU_LABEL_SMOOTH_BETA_FRAC * estimate_nu_highres_extension_beta / real(nvox)
    end function estimate_nu_highres_extension_beta

    module real function estimate_nu_highres_extension_beta_aux( extend_mask, dmat_old, dmat_new )
        logical, intent(in) :: extend_mask(:,:,:)
        real,    intent(in) :: dmat_old(:,:,:), dmat_new(:,:,:)
        integer :: i, j, k, iaux, imask, nvox
        real    :: cur_e, best_e, second_e
        estimate_nu_highres_extension_beta_aux = 0.
        nvox = 0
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.extend_mask(i,j,k) ) cycle
                    imask = nu_mask_index(i,j,k)
                    best_e   = dmat_old(i,j,k)
                    second_e = dmat_new(i,j,k)
                    if( second_e < best_e )then
                        cur_e   = best_e
                        best_e   = second_e
                        second_e = cur_e
                    endif
                    do iaux = 1, size(dmats_aux_mask,2)
                        cur_e = dmats_aux_mask(imask,iaux)
                        if( cur_e < best_e )then
                            second_e = best_e
                            best_e   = cur_e
                        else if( cur_e < second_e )then
                            second_e = cur_e
                        endif
                    end do
                    if( second_e < huge(second_e) )then
                        estimate_nu_highres_extension_beta_aux = estimate_nu_highres_extension_beta_aux + &
                            &max(0., second_e - best_e)
                        nvox = nvox + 1
                    endif
                end do
            end do
        end do
        if( nvox > 0 ) estimate_nu_highres_extension_beta_aux = &
            &NU_LABEL_SMOOTH_BETA_FRAC * estimate_nu_highres_extension_beta_aux / real(nvox)
    end function estimate_nu_highres_extension_beta_aux

    module subroutine init_nu_highres_extension_selection_aux( extend_mask, dmat_old, dmat_new, &
            &extend_choice, n_extended )
        logical, intent(in)    :: extend_mask(:,:,:)
        real,    intent(in)    :: dmat_old(:,:,:), dmat_new(:,:,:)
        integer, intent(inout) :: extend_choice(:,:,:)
        integer, intent(out)   :: n_extended
        integer :: i, j, k, iaux, imask, best_choice
        real    :: best_dmat, cur_dmat
        extend_choice = 0
        n_extended = 0
        !$omp parallel do collapse(3) schedule(static) default(shared) &
        !$omp private(i,j,k,iaux,imask,best_choice,best_dmat,cur_dmat) reduction(+:n_extended)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.extend_mask(i,j,k) ) cycle
                    imask = nu_mask_index(i,j,k)
                    best_choice = 0
                    best_dmat   = dmat_old(i,j,k)
                    if( dmat_new(i,j,k) < best_dmat )then
                        best_choice = 1
                        best_dmat   = dmat_new(i,j,k)
                    endif
                    do iaux = 1, size(dmats_aux_mask,2)
                        cur_dmat = dmats_aux_mask(imask,iaux)
                        if( cur_dmat < best_dmat )then
                            best_choice = iaux + 1
                            best_dmat   = cur_dmat
                        endif
                    end do
                    extend_choice(i,j,k) = best_choice
                    if( best_choice == 1 ) n_extended = n_extended + 1
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine init_nu_highres_extension_selection_aux

    module subroutine refine_nu_highres_extension_selection_aux( extend_mask, dmat_old, dmat_new, &
            &extend_choice, old_label, n_extended )
        logical, intent(in)    :: extend_mask(:,:,:)
        real,    intent(in)    :: dmat_old(:,:,:), dmat_new(:,:,:)
        integer, intent(inout) :: extend_choice(:,:,:)
        integer, intent(in)    :: old_label
        integer, intent(out)   :: n_extended
        integer :: iter, color, i, j, k, n_full(3,NU_LABEL_SMOOTH_NNEIGH), nsz, nchanged
        integer :: choice, cur_choice, best_choice, iaux, imask
        real    :: beta, new_coord, e, best_e
        beta = estimate_nu_highres_extension_beta_aux(extend_mask, dmat_old, dmat_new)
        new_coord = real(old_label + 1)
        if( beta <= TINY )then
            n_extended = count(extend_mask .and. (extend_choice == 1))
            return
        endif
        do iter = 1, NU_LABEL_SMOOTH_MAXITS
            nchanged = 0
            do color = 0, NU_LABEL_SMOOTH_NCOLORS - 1
                !$omp parallel do collapse(3) schedule(static) default(shared) &
                !$omp private(i,j,k,n_full,nsz,choice,cur_choice,best_choice,iaux,imask,e,best_e) &
                !$omp reduction(+:nchanged) proc_bind(close)
                do k = 1, ldim(3)
                    do j = 1, ldim(2)
                        do i = 1, ldim(1)
                            if( .not.extend_mask(i,j,k) ) cycle
                            if( nu_label_smooth_color(i,j,k) /= color ) cycle
                            imask = nu_mask_index(i,j,k)
                            call neigh_8_3D(ldim, [i,j,k], n_full, nsz)
                            cur_choice  = extend_choice(i,j,k)
                            best_choice = cur_choice
                            best_e = nu_highres_extension_choice_unary(i, j, k, imask, cur_choice, &
                                &dmat_old, dmat_new) + beta * &
                                &nu_highres_extension_choice_neighborhood_cost(cur_choice, i, j, k, extend_mask, &
                                &extend_choice, old_label, new_coord, n_full, nsz)
                            do choice = 0, 1
                                if( choice == cur_choice ) cycle
                                e = nu_highres_extension_choice_unary(i, j, k, imask, choice, dmat_old, dmat_new) + &
                                    &beta * nu_highres_extension_choice_neighborhood_cost(choice, i, j, k, &
                                    &extend_mask, extend_choice, old_label, new_coord, n_full, nsz)
                                if( nu_label_smooth_is_better(e, best_e) )then
                                    best_e = e
                                    best_choice = choice
                                endif
                            end do
                            do iaux = 1, size(dmats_aux_mask,2)
                                choice = iaux + 1
                                if( choice == cur_choice ) cycle
                                e = nu_highres_extension_choice_unary(i, j, k, imask, choice, dmat_old, dmat_new) + &
                                    &beta * nu_highres_extension_choice_neighborhood_cost(choice, i, j, k, &
                                    &extend_mask, extend_choice, old_label, new_coord, n_full, nsz)
                                if( nu_label_smooth_is_better(e, best_e) )then
                                    best_e = e
                                    best_choice = choice
                                endif
                            end do
                            if( best_choice /= cur_choice )then
                                extend_choice(i,j,k) = best_choice
                                nchanged = nchanged + 1
                            endif
                        end do
                    end do
                end do
                !$omp end parallel do
            end do
            if( nchanged == 0 ) exit
        end do
        n_extended = count(extend_mask .and. (extend_choice == 1))
    end subroutine refine_nu_highres_extension_selection_aux

    real function nu_highres_extension_choice_unary( i, j, k, imask, choice, dmat_old, dmat_new ) result(e)
        integer, intent(in) :: i, j, k, imask, choice
        real,    intent(in) :: dmat_old(:,:,:), dmat_new(:,:,:)
        if( choice == 0 )then
            e = dmat_old(i,j,k)
        else if( choice == 1 )then
            e = dmat_new(i,j,k)
        else
            e = dmats_aux_mask(imask,choice-1)
        endif
    end function nu_highres_extension_choice_unary

    module real function nu_highres_extension_choice_neighborhood_cost( choice, i, j, k, extend_mask, &
            &extend_choice, old_label, new_coord, neigh, nsz )
        integer, intent(in) :: choice, i, j, k, old_label, neigh(3,NU_LABEL_SMOOTH_NNEIGH), nsz
        logical, intent(in) :: extend_mask(:,:,:)
        integer, intent(in) :: extend_choice(:,:,:)
        real,    intent(in) :: new_coord
        integer :: ineigh, ni, nj, nk, degree, neigh_choice
        real    :: icoord, jcoord
        logical :: l_aux_i, l_aux_j
        call nu_highres_extension_choice_state(i, j, k, choice, old_label, new_coord, icoord, l_aux_i)
        nu_highres_extension_choice_neighborhood_cost = 0.
        degree = 0
        do ineigh = 1, nsz
            ni = neigh(1,ineigh)
            nj = neigh(2,ineigh)
            nk = neigh(3,ineigh)
            if( .not.nu_lmask(ni,nj,nk) ) cycle
            degree = degree + 1
            neigh_choice = 0
            if( extend_mask(ni,nj,nk) ) neigh_choice = extend_choice(ni,nj,nk)
            call nu_highres_extension_choice_state(ni, nj, nk, neigh_choice, old_label, new_coord, jcoord, l_aux_j)
            nu_highres_extension_choice_neighborhood_cost = nu_highres_extension_choice_neighborhood_cost + &
                &nu_label_smooth_source_pair_cost(icoord, jcoord, l_aux_i, l_aux_j)
        end do
        if( degree > 0 ) nu_highres_extension_choice_neighborhood_cost = &
            &nu_highres_extension_choice_neighborhood_cost / real(degree)
    end function nu_highres_extension_choice_neighborhood_cost

    subroutine nu_highres_extension_choice_state( i, j, k, choice, old_label, new_coord, coord, l_aux )
        integer, intent(in)  :: i, j, k, choice, old_label
        real,    intent(in)  :: new_coord
        real,    intent(out) :: coord
        logical, intent(out) :: l_aux
        integer :: iaux
        l_aux = .false.
        if( choice == 1 )then
            coord = new_coord
        else if( choice > 1 )then
            iaux = choice - 1
            coord = nu_candidate_coord_for_label(old_label + iaux)
            l_aux = .true.
        else if( allocated(srcmap) .and. srcmap(i,j,k) > 1 )then
            iaux = srcmap(i,j,k) - 1
            coord = nu_candidate_coord_for_label(old_label + iaux)
            l_aux = .true.
        else
            coord = nu_candidate_coord_for_label(filtmap(i,j,k))
        endif
    end subroutine nu_highres_extension_choice_state

    module real function nu_highres_extension_neighborhood_cost( icoord, extend_mask, extend_to_new, old_label, &
        &new_coord, neigh, nsz )
        real,    intent(in) :: icoord, new_coord
        logical, intent(in) :: extend_mask(:,:,:), extend_to_new(:,:,:)
        integer, intent(in) :: old_label, neigh(3,NU_LABEL_SMOOTH_NNEIGH), nsz
        integer :: ineigh, ni, nj, nk, degree
        real    :: jcoord
        nu_highres_extension_neighborhood_cost = 0.
        degree = 0
        do ineigh = 1, nsz
            ni = neigh(1,ineigh)
            nj = neigh(2,ineigh)
            nk = neigh(3,ineigh)
            if( .not.nu_lmask(ni,nj,nk) ) cycle
            degree = degree + 1
            jcoord = nu_highres_extension_current_coord(ni, nj, nk, extend_mask, extend_to_new, old_label, new_coord)
            nu_highres_extension_neighborhood_cost = nu_highres_extension_neighborhood_cost + &
                &nu_label_smooth_coord_pair_cost(icoord, jcoord)
        end do
        if( degree > 0 ) nu_highres_extension_neighborhood_cost = &
            &nu_highres_extension_neighborhood_cost / real(degree)
    end function nu_highres_extension_neighborhood_cost

    module real function nu_highres_extension_current_coord( i, j, k, extend_mask, extend_to_new, old_label, new_coord )
        integer, intent(in) :: i, j, k, old_label
        logical, intent(in) :: extend_mask(:,:,:), extend_to_new(:,:,:)
        real,    intent(in) :: new_coord
        integer :: aux_label
        if( extend_mask(i,j,k) )then
            if( extend_to_new(i,j,k) )then
                nu_highres_extension_current_coord = new_coord
            else
                nu_highres_extension_current_coord = nu_candidate_coord_for_label(old_label)
            endif
        else
            if( allocated(srcmap) )then
                if( srcmap(i,j,k) > 1 )then
                    aux_label = old_label + srcmap(i,j,k) - 1
                    nu_highres_extension_current_coord = nu_candidate_coord_for_label(aux_label)
                    return
                endif
            endif
            nu_highres_extension_current_coord = nu_candidate_coord_for_label(filtmap(i,j,k))
        endif
    end function nu_highres_extension_current_coord

    module subroutine apply_nu_highres_extension_selection( extend_to_new, new_label )
        logical, intent(in) :: extend_to_new(:,:,:)
        integer, intent(in) :: new_label
        integer :: i, j, k
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.extend_to_new(i,j,k) ) cycle
                    srcmap(i,j,k)  = 1
                    filtmap(i,j,k) = new_label
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine apply_nu_highres_extension_selection

    module subroutine apply_nu_highres_extension_selection_aux( extend_choice, old_label, new_label )
        integer, intent(in) :: extend_choice(:,:,:)
        integer, intent(in) :: old_label, new_label
        integer :: i, j, k, choice, aux_icand
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k,choice,aux_icand)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    choice = extend_choice(i,j,k)
                    if( choice == 0 ) cycle
                    if( choice == 1 )then
                        srcmap(i,j,k)  = 1
                        filtmap(i,j,k) = new_label
                    else
                        srcmap(i,j,k) = choice
                        aux_icand = old_label + choice - 1
                        filtmap(i,j,k) = nu_effective_base_label_for_candidate(aux_icand, old_label)
                    endif
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine apply_nu_highres_extension_selection_aux

    subroutine cache_nu_highres_extension_frontier_after_aux_selection( dmat_new, new_label )
        real,    intent(in) :: dmat_new(:,:,:)
        integer, intent(in) :: new_label
        integer :: i, j, k, imask, iaux
        allocate(dmat_finest_cached(ldim(1),ldim(2),ldim(3)), source=huge(0.))
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k,imask,iaux)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.nu_lmask(i,j,k) ) cycle
                    if( filtmap(i,j,k) /= new_label ) cycle
                    if( srcmap(i,j,k) == 1 )then
                        dmat_finest_cached(i,j,k) = dmat_new(i,j,k)
                    else
                        iaux = srcmap(i,j,k) - 1
                        if( iaux < 1 .or. iaux > size(dmats_aux_mask,2) ) cycle
                        imask = nu_mask_index(i,j,k)
                        dmat_finest_cached(i,j,k) = dmats_aux_mask(imask,iaux)
                    endif
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine cache_nu_highres_extension_frontier_after_aux_selection

    module subroutine append_nu_highres_candidate_coord( old_n_base, new_coord )
        integer, intent(in) :: old_n_base
        real,    intent(in) :: new_coord
        real, allocatable :: new_coords(:)
        integer :: n_aux
        if( .not.allocated(candidate_coords) ) return
        if( size(candidate_coords) < old_n_base ) &
            &THROW_HARD('candidate_coords shorter than base bank in append_nu_highres_candidate_coord')
        n_aux = size(candidate_coords) - old_n_base
        allocate(new_coords(old_n_base + 1 + n_aux), source=0.)
        new_coords(:old_n_base) = candidate_coords(:old_n_base)
        new_coords(old_n_base + 1) = new_coord
        if( n_aux > 0 ) new_coords(old_n_base + 2:) = candidate_coords(old_n_base + 1:)
        call move_alloc(new_coords, candidate_coords)
    end subroutine append_nu_highres_candidate_coord

end submodule simple_nu_filter_extend
