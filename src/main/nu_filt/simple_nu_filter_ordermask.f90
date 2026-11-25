!@descr: density-domain order attenuation of nonuniform matching references
! automsk=order: M = 1/(1 + NU_ORDER_LAMBDA (1 - w_f)), w_f a logistic of the NU evidence margin at the (1-f)
! quantile over the observed dense set, f the smallest that keeps NU_ORDER_CORE_KEEP of the core's cross-half
! signal. Serial on purpose: the choice of f must not depend on a reduction order. Design and rulings:
! doc/implementation_notes/planned/automsk_order.md.
submodule (simple_nu_filter) simple_nu_filter_ordermask
implicit none
#include "simple_local_flags.inc"
contains

    module subroutine select_nu_order_mask( margin, locres, dense, even, odd, res05, factor, st, coarsest_lp )
        real,                 intent(in)                 :: margin(:,:,:), locres(:,:,:), even(:,:,:), odd(:,:,:), res05
        logical,              intent(in)                 :: dense(:,:,:)
        real,                 intent(out)                :: factor(:,:,:)
        type(nu_order_stats), intent(out)                :: st
        real,                 optional,       intent(in) :: coarsest_lp
        real,    allocatable :: vals(:), core_margin(:), core_even(:), core_odd(:)
        logical, allocatable :: core(:,:,:)
        real(dp)             :: baseline, kept, sm, sm2, peri, peri_kept, fine, fine_kept
        real                 :: threshold, width, arg, w, m, coarse_resolution
        integer              :: nd, nc, ifr, it, i, j, k, v
        st = nu_order_stats()
        factor = 1.
        coarse_resolution = lowpass_limits(1)
        if( present(coarsest_lp) ) coarse_resolution = coarsest_lp
        if( any(shape(margin) /= shape(dense)) .or. any(shape(locres) /= shape(dense)) .or. &
            &any(shape(even) /= shape(dense)) .or. any(shape(odd) /= shape(dense)) .or. &
            &any(shape(factor) /= shape(dense)) ) THROW_HARD('order mask input dimensions differ')
        nd = count(dense)
        st%n_dense = nd
        if( res05 <= 0. )then
            st%reason = 'no positive FSC=0.5 resolution'
            return
        endif
        allocate(core(size(dense,1),size(dense,2),size(dense,3)))
        core = dense .and. locres > 0. .and. locres <= NU_ORDER_CORE_RES_FAC * res05
        nc = count(core)
        st%n_core = nc
        if( nd > 0 ) st%core_share = real(nc) / real(nd)
        if( nc < NU_ORDER_CORE_MIN_VOX )then
            st%reason = 'fewer than 1000 core voxels'
            return
        endif
        vals = pack(margin, dense)
        call hpsort(vals)
        core_margin = pack(margin, core)
        core_even = pack(even, core)
        core_odd = pack(odd, core)
        baseline = sum(real(core_even,dp) * real(core_odd,dp))
        if( baseline <= 0._dp )then
            st%reason = 'nonpositive core cross-half signal'
            return
        endif
        do ifr = 1, 99
            it = min(nd, max(2, nint((1. - real(ifr)*NU_ORDER_FRAC_STEP) * real(nd))))
            threshold = vals(it)
            if( threshold <= 0. ) cycle
            sm = sum(real(vals(1:it-1),dp)) / real(it-1,dp)
            sm2 = sum(real(vals(1:it-1),dp)**2) / real(it-1,dp)
            width = max(real(sqrt(max(sm2-sm**2,0._dp))), 1.e-3*threshold)
            kept = 0._dp
            do v = 1, nc
                arg = min(40., max(-40., (core_margin(v)-threshold)/width))
                w = 1. / (1. + exp(-arg))
                m = 1. / (1. + NU_ORDER_LAMBDA * (1.-w))
                kept = kept + real(m,dp)**2 * real(core_even(v),dp) * real(core_odd(v),dp)
            enddo
            if( kept / baseline < real(NU_ORDER_CORE_KEEP,dp) ) cycle
            st%fraction = real(ifr)*NU_ORDER_FRAC_STEP
            st%threshold = threshold
            st%width = width
            st%kept_core = real(kept/baseline)
            st%l_off = .false.
            exit
        enddo
        if( st%l_off )then
            st%reason = 'no fraction satisfies the core guard'
            return
        endif
        peri = 0._dp
        peri_kept = 0._dp
        fine = 0._dp
        fine_kept = 0._dp
        do k=1,size(dense,3)
            do j=1,size(dense,2)
                do i=1,size(dense,1)
                    arg = min(40., max(-40., (margin(i,j,k)-st%threshold)/st%width))
                    w = 1. / (1. + exp(-arg))
                    m = 1. / (1. + NU_ORDER_LAMBDA * (1.-w))
                    factor(i,j,k) = m
                    if( .not.dense(i,j,k) .or. core(i,j,k) ) cycle
                    sm = real(even(i,j,k),dp)*real(odd(i,j,k),dp)
                    peri = peri + sm
                    peri_kept = peri_kept + real(m,dp)**2*sm
                    if( locres(i,j,k) > 0. .and. locres(i,j,k) < coarse_resolution )then
                        fine = fine + sm
                        fine_kept = fine_kept + real(m,dp)**2*sm
                    endif
                enddo
            enddo
        enddo
        if( peri > 0._dp )then
            st%l_peri_signal = .true.
            st%removed_peri = real(1._dp-peri_kept/peri)
        endif
        if( fine > 0._dp )then
            st%l_fine_signal = .true.
            st%kept_fine_peri = real(fine_kept/fine)
        endif
    end subroutine select_nu_order_mask

    module subroutine make_nu_order_mask( core_image, even, odd, res05, lp_smooth, factor, st )
        class(image),         intent(in)    :: core_image, even, odd
        real,                 intent(in)    :: res05, lp_smooth
        type(image),          intent(inout) :: factor
        type(nu_order_stats), intent(out)   :: st
        type(image)                :: evmap, locres
        real(c_float), pointer     :: ev(:,:,:), lr(:,:,:), c(:,:,:), e(:,:,:), o(:,:,:), f(:,:,:)
        logical,       allocatable :: dense(:,:,:)
        integer                    :: imask, i, j, k, omitted
        call get_nu_evidence_map(evmap, lp_smooth, NU_ENVMASK_RELATIVE)
        call get_nu_local_resolution_map(locres)
        call evmap%get_rmat_ptr(ev)
        call locres%get_rmat_ptr(lr)
        call core_image%get_rmat_ptr(c)
        call even%get_rmat_ptr(e)
        call odd%get_rmat_ptr(o)
        call factor%new(ldim, smpd, wthreads=.false.)
        call factor%get_rmat_ptr(f)
        allocate(dense(ldim(1),ldim(2),ldim(3)), source=.false.)
        omitted = 0
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            if( c(i,j,k) < 0.5 ) cycle
            if( nu_observed_mask(imask) )then
                dense(i,j,k) = .true.
            else
                omitted = omitted + 1
            endif
        enddo
        call select_nu_order_mask(ev(:ldim(1),:ldim(2),:ldim(3)), lr(:ldim(1),:ldim(2),:ldim(3)), dense, &
            &e(:ldim(1),:ldim(2),:ldim(3)), o(:ldim(1),:ldim(2),:ldim(3)), res05, &
            &f(:ldim(1),:ldim(2),:ldim(3)), st, coarsest_lp=nu_label_lowpass_limit(1))
        st%n_unobserved = omitted
        nullify(ev, lr, c, e, o, f)
        call evmap%kill
        call locres%kill
    end subroutine make_nu_order_mask

    module subroutine report_nu_order_mask( state, st )
        integer,              intent(in) :: state
        type(nu_order_stats), intent(in) :: st
        character(len=24) :: removed, kept_fine
        if( st%l_off )then
            write(logfhandle,'(A,I0,A,A,A,I0)') '>>> NU ORDER MASK: STATE ', state, &
                &', no attenuation this iteration: ', trim(st%reason), '; unobserved dense voxels ', st%n_unobserved
        else
            removed = 'unavailable'
            kept_fine = 'unavailable'
            if( st%l_peri_signal ) write(removed,'(F10.4)') st%removed_peri
            if( st%l_fine_signal ) write(kept_fine,'(F10.4)') st%kept_fine_peri
            write(logfhandle,'(A,I0,A,F5.2,A,ES11.3,A,ES11.3,A,F7.4,A,A,A,A,A,F7.4,A,I0)') &
                &'>>> NU ORDER MASK: STATE ', state, ', f ', st%fraction, ', threshold ', st%threshold, &
                &', width ', st%width, ', K_core ', st%kept_core, ', removed around core ', trim(adjustl(removed)), &
                &', kept fine outside core ', trim(adjustl(kept_fine)), ', core share ', st%core_share, &
                &', unobserved dense voxels ', st%n_unobserved
        endif
    end subroutine report_nu_order_mask
end submodule simple_nu_filter_ordermask
