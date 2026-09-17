!@descr: flex_pca projection-aware latent model: Fourier projection/backprojection helpers, particle prep, the coupled M-step solve
module simple_flex_reconstructor_latent_ops
use simple_flex_pca_plane_cache, only: plane_cache_fill
use simple_core_module_api
use simple_reconstructor, only: reconstructor
use simple_builder,          only: builder
use simple_image,            only: image
use simple_linalg,           only: eigsrt, jacobi
use simple_memoize_ft_maps,  only: memoize_ft_maps
use simple_parameters,       only: parameters
implicit none

public :: insert_planes_oversamp_multi_scaled_batch
public :: insert_planes_oversamp_coupled_batch_scaled
public :: project_fplane_mean, project_fplanes_mean_basis
!> exported for the polar (G,b) former, which must use the SAME interpolation kernel and the SAME
!! weight normalisation as the Cartesian projector or the two paths cannot be compared
public :: latent_projection_weights, weighted_expanded_cmat, LATENT_WDIM
!> the projection-aware latent model (merged from simple_flex_projected_latent_model)
public :: prep_imgs4projected_model, solve_coupled_basis_exp, projected_model_kfromto
public :: add_invtausq2rho_coupled, pair_index
public :: cap_fplane_for_projected_model, flex_dev_prep_hook
private
!> Device prep hook. The device variant of prep_imgs4projected_model lives in simple_flex_gpu,
!! which itself uses this module; the hook (set by flex_gpu_prep_begin_f, cleared by
!! flex_gpu_prep_free_f) lets the shared prep funnel take the device branch without a
!! module dependency cycle. Unassociated = CPU prep.
abstract interface
    subroutine flex_dev_prep_iface( params, build, nptcls, ptcl_imgs, pinds, fplanes, fetch )
        import :: parameters, builder, image, fplane_type
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls
        class(image),      intent(inout) :: ptcl_imgs(nptcls)
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), intent(inout) :: fplanes(nptcls)
        logical, optional, intent(in)    :: fetch
    end subroutine flex_dev_prep_iface
end interface
procedure(flex_dev_prep_iface), pointer :: flex_dev_prep_hook => null()
#include "simple_local_flags.inc"

integer, parameter :: LATENT_WDIM = 2 * ceiling(KBWINSZ - 0.5) + 1
! cmat_exp stores h>=0 as the independent Friedel half; h<0 is only
! interpolation halo and must not receive independent projection samples.
integer, parameter :: NONREDUNDANT_HMIN = 0
! Source h-lines in one OpenMP colour must map to non-overlapping 3-D
! interpolation windows for every rotation. A separation of LATENT_WDIM in
! the source plane is not sufficient after rotation; sqrt(3)*LATENT_WDIM
! guarantees that at least one target-grid coordinate differs by a full
! window width.
integer, parameter :: LATENT_SAFE_STRIDE = ceiling(sqrt(3.0) * real(LATENT_WDIM))

real(dp), parameter :: COUPLED_MSTEP_RIDGE_REL = 1.0d-8
real(dp), parameter :: COUPLED_DENSITY_FLOOR = 1.0d-6

contains

    ! Batched insert_plane_oversamp_multi_scaled: identical arithmetic, one OpenMP region per batch
    ! instead of one per particle. Same shape as insert_planes_oversamp_coupled_batch_scaled -- every
    ! per-record quantity is derived serially up front (se%apply goes through oris%get_ori, not
    ! guaranteed thread-safe) and only the h-line sweep is threaded, with the stride keeping
    ! concurrent lines off any one interpolation cell, so no privatisation and no atomics. Bit-for-bit
    ! identical to calling the per-particle routine in a loop.
    subroutine insert_planes_oversamp_multi_scaled_batch( recs, se, orientations, fpls, &
        &data_scales, density_scales, valid, nrecords )
        use simple_math, only: ceil_div, floor_div
        type(reconstructor), intent(inout) :: recs(:)
        class(sym),          intent(inout) :: se
        type(ori),           intent(inout) :: orientations(:)
        type(fplane_type),   intent(in)    :: fpls(:)
        real(dp),            intent(in)    :: data_scales(:,:), density_scales(:,:)
        logical,             intent(in)    :: valid(:)
        integer,             intent(in)    :: nrecords
        type(kbinterpol) :: kbwin
        type(ori) :: o_sym
        complex   :: comp_base, cmplx_raw
        real, allocatable :: rotmats(:,:,:,:), dscale(:,:), rscale(:,:)
        integer, allocatable :: fpllims(:,:,:), nyq_disks(:), act_all(:,:), nact_all(:)
        real      :: loc(3), hrow(3), ctfsq_raw
        real      :: wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM), ww
        real      :: r11, r12, r13, r21, r22, r23
        integer   :: win(2,3), h, k, l, nsym, isym, iwinsz, stride, fpllims_pd(3,2)
        integer   :: hp, kp, pf, ix, iy, iz, hx, ky, mz, q, iq, ncomp, i, nact
        integer   :: nyq_eff, h_sq, k_max_h, k_lo, k_hi, exp_lb(3), exp_ub(3)
        real      :: pf2, eps_norm, inv_wdim
        ncomp = size(recs)
        if( ncomp <= 0 .or. nrecords <= 0 ) return
        if( size(orientations) < nrecords .or. size(fpls) < nrecords .or. size(valid) < nrecords )then
            THROW_HARD('record array smaller than batch; insert_planes_oversamp_multi_scaled_batch')
        endif
        if( size(data_scales,1) < ncomp .or. size(data_scales,2) < nrecords .or. &
            &size(density_scales,1) < ncomp .or. size(density_scales,2) < nrecords )then
            THROW_HARD('scale array smaller than batch; insert_planes_oversamp_multi_scaled_batch')
        endif
        if( .not. allocated(recs(1)%cmat_exp) )then
            THROW_HARD('expanded matrix does not exist; insert_planes_oversamp_multi_scaled_batch')
        endif
        kbwin    = kbinterpol(KBWINSZ, KBALPHA)
        iwinsz   = ceiling(KBWINSZ - 0.5)
        stride   = LATENT_SAFE_STRIDE
        exp_lb   = lbound(recs(1)%cmat_exp)
        exp_ub   = ubound(recs(1)%cmat_exp)
        nsym     = se%get_nsym()
        pf       = OSMPL_PAD_FAC
        pf2      = real(pf*pf)
        eps_norm = epsilon(1.0)
        inv_wdim = 1.0 / real(LATENT_WDIM)
        allocate(rotmats(3,3,nsym,nrecords), source=0.)
        allocate(dscale(ncomp,nrecords), rscale(ncomp,nrecords), source=0.)
        allocate(fpllims(3,2,nrecords), nyq_disks(nrecords), nact_all(nrecords), source=0)
        allocate(act_all(ncomp,nrecords), source=0)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            ! live components only -- the kernel state weights have compact support, so a particle
            ! typically sits inside one or two of the nstates targets and the rest are exact zeros
            nact = 0
            do q = 1, ncomp
                dscale(q,i) = real(data_scales(q,i))
                rscale(q,i) = real(max(0.d0, density_scales(q,i)))
                if( dscale(q,i) /= 0. .or. rscale(q,i) /= 0. )then
                    nact = nact + 1
                    act_all(nact,i) = q
                endif
            end do
            nact_all(i) = nact
            if( nact == 0 ) cycle
            rotmats(:,:,1,i) = orientations(i)%get_mat()
            do isym = 2, nsym
                call se%apply(orientations(i), isym, o_sym)
                rotmats(:,:,isym,i) = o_sym%get_mat()
            end do
            fpllims_pd     = fpls(i)%frlims
            fpllims(:,:,i) = fpllims_pd
            fpllims(1,1,i) = ceil_div (fpllims_pd(1,1), pf)
            fpllims(1,2,i) = floor_div(fpllims_pd(1,2), pf)
            fpllims(2,1,i) = ceil_div (fpllims_pd(2,1), pf)
            fpllims(2,2,i) = floor_div(fpllims_pd(2,2), pf)
            nyq_eff = recs(1)%get_lfny(1)
            if( fpls(i)%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpls(i)%nyq / pf))
            nyq_disks(i) = nyq_eff * (nyq_eff + 1)
        end do
        call o_sym%kill
        !$omp parallel default(shared) private(i,h,k,l,h_sq,k_max_h,k_lo,k_hi,cmplx_raw,&
        !$omp& ctfsq_raw,comp_base,wx,wy,wz,ww,win,loc,hrow,hp,kp,r11,r12,r13,r21,r22,r23,&
        !$omp& isym,ix,iy,iz,hx,ky,mz,q,iq,nact) proc_bind(close)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            nact = nact_all(i)
            if( nact == 0 ) cycle
            do isym = 1, nsym
                r11 = rotmats(1,1,isym,i); r12 = rotmats(1,2,isym,i); r13 = rotmats(1,3,isym,i)
                r21 = rotmats(2,1,isym,i); r22 = rotmats(2,2,isym,i); r23 = rotmats(2,3,isym,i)
                do l = 0, stride-1
                    !$omp do schedule(static,1)
                    do h = fpllims(1,1,i)+l, fpllims(1,2,i), stride
                        h_sq = h*h
                        if( h_sq > nyq_disks(i) ) cycle
                        k_max_h = int(sqrt(real(nyq_disks(i) - h_sq)))
                        k_lo    = max(fpllims(2,1,i), -k_max_h)
                        k_hi    = min(fpllims(2,2,i),  k_max_h)
                        hp      = h * pf
                        hrow(1) = real(h) * r11
                        hrow(2) = real(h) * r12
                        hrow(3) = real(h) * r13
                        do k = k_lo, k_hi
                            kp = k * pf
                            if( kp <= 0 )then
                                cmplx_raw = fpls(i)%cmplx_plane(hp,kp)
                                ctfsq_raw = fpls(i)%ctfsq_plane(hp,kp)
                            else
                                cmplx_raw = conjg(fpls(i)%cmplx_plane(-hp,-kp))
                                ctfsq_raw = fpls(i)%ctfsq_plane(-hp,-kp)
                            endif
                            if( abs(real(cmplx_raw)) + abs(aimag(cmplx_raw)) <= TINY .and. &
                                &ctfsq_raw <= TINY ) cycle
                            loc(1) = hrow(1) + real(k) * r21
                            loc(2) = hrow(2) + real(k) * r22
                            loc(3) = hrow(3) + real(k) * r23
                            win(1,:) = nint(loc)
                            win(2,:) = win(1,:) + iwinsz
                            win(1,:) = win(1,:) - iwinsz
                            if( win(2,1) < NONREDUNDANT_HMIN ) cycle
                            if( any(win(1,:) < exp_lb) .or. any(win(2,:) > exp_ub) ) cycle
                            comp_base = pf2 * cmplx_raw
                            call kb_apod_vecs_3d_fast_b(loc, wx, wy, wz)
                            do iz = 1, LATENT_WDIM
                                mz = win(1,3) + iz - 1
                                do iy = 1, LATENT_WDIM
                                    ky = win(1,2) + iy - 1
                                    do ix = 1, LATENT_WDIM
                                        hx = win(1,1) + ix - 1
                                        ww = wx(ix) * (wy(iy) * wz(iz))
                                        do iq = 1, nact
                                            q = act_all(iq,i)
                                            recs(q)%cmat_exp(hx,ky,mz) = recs(q)%cmat_exp(hx,ky,mz) + &
                                                &(dscale(q,i) * comp_base) * ww
                                            recs(q)%rho_exp(hx,ky,mz) = recs(q)%rho_exp(hx,ky,mz) + &
                                                &(rscale(q,i) * ctfsq_raw) * ww
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                    !$omp end do
                end do
            end do
        end do
        !$omp end parallel
        deallocate(rotmats, dscale, rscale, fpllims, nyq_disks, act_all, nact_all)

    contains

        subroutine kb_apod_vecs_3d_fast_b( loc, wx, wy, wz )
            real, intent(in)  :: loc(3)
            real, intent(out) :: wx(:), wy(:), wz(:)
            integer :: i2, win_lo(3)
            real    :: base(3), ww3(3), sx, sy, sz
            win_lo = nint(loc) - iwinsz
            base   = real(win_lo) - loc
            do i2 = 1, LATENT_WDIM
                ww3    = kbwin%apod_fast(base + real(i2-1))
                wx(i2) = ww3(1)
                wy(i2) = ww3(2)
                wz(i2) = ww3(3)
            end do
            sx = sum(wx)
            sy = sum(wy)
            sz = sum(wz)
            if( abs(sx) > eps_norm )then
                wx = wx * (1.0 / sx)
            else
                wx = inv_wdim
            endif
            if( abs(sy) > eps_norm )then
                wy = wy * (1.0 / sy)
            else
                wy = inv_wdim
            endif
            if( abs(sz) > eps_norm )then
                wz = wz * (1.0 / sz)
            else
                wz = inv_wdim
            endif
        end subroutine kb_apod_vecs_3d_fast_b

    end subroutine insert_planes_oversamp_multi_scaled_batch



    subroutine insert_planes_oversamp_coupled_batch_scaled( recs, rho_cross_exp, se, orientations, fpls, &
        &data_scales, density_scales, valid, nrecords )
        use simple_math, only: ceil_div, floor_div
        type(reconstructor), intent(inout) :: recs(:)
        real,                intent(inout) :: rho_cross_exp(:,:,:,:)
        class(sym),          intent(inout) :: se
        type(ori),           intent(inout) :: orientations(:)
        type(fplane_type),   intent(in)    :: fpls(:)
        real(dp),            intent(in)    :: data_scales(:,:), density_scales(:,:,:)
        logical,             intent(in)    :: valid(:)
        integer,             intent(in)    :: nrecords
        type(ori) :: o_sym
        type(kbinterpol) :: kbwin
        complex   :: comp_base, cmplx_raw
        real, allocatable :: rotmats(:,:,:,:), data_scale_sp(:,:), density_scale_packed(:,:)
        integer, allocatable :: fpllims(:,:,:), nyq_disks(:)
        real      :: loc(3), hrow(3), ctfsq_raw
        real      :: wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM), ww
        real      :: r11, r12, r13, r21, r22, r23
        integer   :: win(2,3), h, k, l, nsym, isym, iwinsz, stride, fpllims_pd(3,2)
        integer   :: hp, kp, pf, ix, iy, iz, hx, ky, mz, q, r, i, ncomp, ipair
        integer   :: h_sq, k_max_h, k_lo, k_hi, ih, ik, im, nyq_eff
        integer   :: exp_lb(3), exp_ub(3), exp_shape(3), npairs
        logical   :: shared_density, diagonal_density
        real      :: pf2, eps_norm, inv_wdim
        ncomp = size(recs)
        if( ncomp <= 0 .or. nrecords <= 0 ) return
        npairs = (ncomp * (ncomp + 1)) / 2
        diagonal_density = size(rho_cross_exp,1) == ncomp
        shared_density = size(rho_cross_exp,1) == 1 .and. .not.diagonal_density
        if( size(orientations)<nrecords .or. size(fpls)<nrecords .or. size(valid)<nrecords )then
            THROW_HARD('record array smaller than batch; insert_planes_oversamp_coupled_batch_scaled')
        endif
        if( size(data_scales,1)<ncomp .or. size(data_scales,2)<nrecords .or. &
            &size(density_scales,1)<ncomp .or. size(density_scales,2)<ncomp .or. &
            &size(density_scales,3)<nrecords )then
            THROW_HARD('scale array smaller than batch; insert_planes_oversamp_coupled_batch_scaled')
        endif
        if( .not.allocated(recs(1)%cmat_exp) )then
            THROW_HARD('expanded matrix does not exist; insert_planes_oversamp_coupled_batch_scaled')
        endif
        exp_lb    = lbound(recs(1)%cmat_exp)
        exp_ub    = ubound(recs(1)%cmat_exp)
        exp_shape = shape(recs(1)%cmat_exp)
        if( (.not.shared_density .and. .not.diagonal_density .and. size(rho_cross_exp,1)<npairs) .or. &
            &size(rho_cross_exp,2)<exp_shape(1) .or. &
            &size(rho_cross_exp,3)<exp_shape(2) .or. size(rho_cross_exp,4)<exp_shape(3) )then
            THROW_HARD('cross-density array shape mismatch; insert_planes_oversamp_coupled_batch_scaled')
        endif
        nsym     = se%get_nsym()
        iwinsz   = ceiling(KBWINSZ - 0.5)
        stride   = LATENT_SAFE_STRIDE
        pf       = OSMPL_PAD_FAC
        pf2      = real(pf*pf)
        eps_norm = epsilon(1.0)
        inv_wdim = 1.0 / real(LATENT_WDIM)
        ! Use the reconstructor's own interpolation window.  The previous
        ! hand-inlined approximation is close in raw accumulation, but its
        ! small differences become large after density correction in weakly
        ! sampled Fourier cells.
        kbwin = recs(1)%get_kbwin()
        allocate(rotmats(3,3,nsym,nrecords), data_scale_sp(ncomp,nrecords), source=0.)
        if( .not.shared_density .and. .not.diagonal_density ) allocate(density_scale_packed(npairs,nrecords), source=0.)
        allocate(fpllims(3,2,nrecords), nyq_disks(nrecords), source=0)
        do i = 1, nrecords
            if( .not.valid(i) ) cycle
            if( .not.allocated(fpls(i)%transfer_plane) )then
                THROW_HARD('forward transfer plane does not exist; insert_planes_oversamp_coupled_batch_scaled')
            endif
            rotmats(:,:,1,i) = orientations(i)%get_mat()
            do isym = 2, nsym
                call se%apply(orientations(i), isym, o_sym)
                rotmats(:,:,isym,i) = o_sym%get_mat()
            end do
            fpllims_pd = fpls(i)%frlims
            fpllims(:,:,i) = fpllims_pd
            fpllims(1,1,i) = ceil_div (fpllims_pd(1,1), pf)
            fpllims(1,2,i) = floor_div(fpllims_pd(1,2), pf)
            fpllims(2,1,i) = ceil_div (fpllims_pd(2,1), pf)
            fpllims(2,2,i) = floor_div(fpllims_pd(2,2), pf)
            nyq_eff = recs(1)%get_lfny(1)
            if( fpls(i)%nyq>0 ) nyq_eff = min(nyq_eff, max(1, fpls(i)%nyq/pf))
            nyq_disks(i) = nyq_eff * (nyq_eff + 1)
            do q = 1, ncomp
                data_scale_sp(q,i) = real(data_scales(q,i))
            end do
            if( .not.shared_density .and. .not.diagonal_density )then
                do r = 1, ncomp
                    do q = 1, r
                        density_scale_packed(pair_index(q,r),i) = real(density_scales(q,r,i))
                    end do
                end do
            endif
        end do
        call o_sym%kill
        !$omp parallel default(shared) private(i,h,k,l,h_sq,k_max_h,k_lo,k_hi,cmplx_raw,ctfsq_raw,&
        !$omp& comp_base,wx,wy,wz,ww,win,loc,hrow,hp,kp,r11,r12,r13,r21,r22,r23,isym,&
        !$omp& ix,iy,iz,hx,ky,mz,ih,ik,im,q,ipair) proc_bind(close)
        do i = 1, nrecords
            if( .not.valid(i) ) cycle
            do isym = 1, nsym
                r11 = rotmats(1,1,isym,i); r12 = rotmats(1,2,isym,i); r13 = rotmats(1,3,isym,i)
                r21 = rotmats(2,1,isym,i); r22 = rotmats(2,2,isym,i); r23 = rotmats(2,3,isym,i)
                do l = 0, stride-1
                    !$omp do schedule(static,1)
                    do h = fpllims(1,1,i)+l, fpllims(1,2,i), stride
                        h_sq = h*h
                        if( h_sq>nyq_disks(i) ) cycle
                        k_max_h = int(sqrt(real(nyq_disks(i)-h_sq)))
                        k_lo = max(fpllims(2,1,i),-k_max_h)
                        k_hi = min(fpllims(2,2,i), k_max_h)
                        hp = h*pf
                        hrow = real(h)*[r11,r12,r13]
                        loc = hrow + real(k_lo-1)*[r21,r22,r23]
                        do k = k_lo, k_hi
                            loc = loc + [r21,r22,r23]
                            kp = k*pf
                            if( kp<=0 )then
                                cmplx_raw = conjg(fpls(i)%transfer_plane(hp,kp))*fpls(i)%cmplx_plane(hp,kp)
                                ctfsq_raw = fpls(i)%ctfsq_plane(hp,kp)
                            else
                                cmplx_raw = conjg(conjg(fpls(i)%transfer_plane(-hp,-kp))*fpls(i)%cmplx_plane(-hp,-kp))
                                ctfsq_raw = fpls(i)%ctfsq_plane(-hp,-kp)
                            endif
                            if( abs(real(cmplx_raw))+abs(aimag(cmplx_raw))<=TINY .and. ctfsq_raw<=TINY ) cycle
                            win(1,:) = nint(loc)-iwinsz
                            win(2,:) = nint(loc)+iwinsz
                            if( win(2,1) < NONREDUNDANT_HMIN ) cycle
                            if( any(win(1,:)<exp_lb) .or. any(win(2,:)>exp_ub) ) cycle
                            comp_base = pf2*cmplx_raw
                            call kb_apod_vecs_3d_fast(loc,wx,wy,wz)
                            do iz = 1, LATENT_WDIM
                                mz = win(1,3)+iz-1
                                im = mz-exp_lb(3)+1
                                do iy = 1, LATENT_WDIM
                                    ky = win(1,2)+iy-1
                                    ik = ky-exp_lb(2)+1
                                    do ix = 1, LATENT_WDIM
                                        hx = win(1,1)+ix-1
                                        ih = hx-exp_lb(1)+1
                                        ww = wx(ix)*(wy(iy)*wz(iz))
                                        do q = 1, ncomp
                                            recs(q)%cmat_exp(hx,ky,mz) = recs(q)%cmat_exp(hx,ky,mz) + &
                                                &(data_scale_sp(q,i)*comp_base)*ww
                                        end do
                                        if( shared_density )then
                                            rho_cross_exp(1,ih,ik,im) = rho_cross_exp(1,ih,ik,im) + ctfsq_raw*ww
                                        else if( diagonal_density )then
                                            do q = 1, ncomp
                                                rho_cross_exp(q,ih,ik,im) = rho_cross_exp(q,ih,ik,im) + &
                                                    &real(density_scales(q,q,i))*ctfsq_raw*ww
                                            end do
                                        else
                                            !$omp simd
                                            do ipair = 1, npairs
                                                rho_cross_exp(ipair,ih,ik,im) = rho_cross_exp(ipair,ih,ik,im) + &
                                                    &(density_scale_packed(ipair,i)*ctfsq_raw)*ww
                                            end do
                                        endif
                                    end do
                                end do
                            end do
                        end do
                    end do
                    !$omp end do
                end do
            end do
        end do
        !$omp end parallel
        deallocate(rotmats,data_scale_sp,fpllims,nyq_disks)
        if( allocated(density_scale_packed) ) deallocate(density_scale_packed)

    contains

        integer pure function pair_index( q, r ) result( ipair )
            integer, intent(in) :: q, r
            ipair = (r*(r-1))/2+q
        end function pair_index

        subroutine kb_apod_vecs_3d_fast( loc, wx, wy, wz )
            real, intent(in)  :: loc(3)
            real, intent(out) :: wx(:), wy(:), wz(:)
            integer :: j, win_lo(3)
            real :: base(3), ww3(3), sx, sy, sz
            win_lo = nint(loc)-iwinsz
            base = real(win_lo)-loc
            do j = 1, LATENT_WDIM
                ww3=kbwin%apod_fast(base+real(j-1))
                wx(j)=ww3(1)
                wy(j)=ww3(2)
                wz(j)=ww3(3)
            end do
            sx=sum(wx); sy=sum(wy); sz=sum(wz)
            if( abs(sx)>eps_norm )then; wx=wx/sx; else; wx=inv_wdim; endif
            if( abs(sy)>eps_norm )then; wy=wy/sy; else; wy=inv_wdim; endif
            if( abs(sz)>eps_norm )then; wz=wz/sz; else; wz=inv_wdim; endif
        end subroutine kb_apod_vecs_3d_fast

    end subroutine insert_planes_oversamp_coupled_batch_scaled



    subroutine project_fplane_mean( mean_rec, o, fpl_ref, mean_fpl, apply_ctf_amp )
        type(reconstructor), intent(in)    :: mean_rec
        class(ori),          intent(inout) :: o
        class(fplane_type),  intent(in)    :: fpl_ref
        type(fplane_type),   intent(inout) :: mean_fpl
        logical, optional,   intent(in)    :: apply_ctf_amp
        call mean_rec%project_fplane(o, fpl_ref, mean_fpl, apply_ctf_amp)
    end subroutine project_fplane_mean

    subroutine project_fplanes_mean_basis( mean_rec, basis_recs, o, fpl_ref, mean_fpl, basis_fpls, apply_ctf_amp )
        use simple_math, only: ceil_div, floor_div
        type(reconstructor), intent(in)    :: mean_rec
        type(reconstructor), intent(in)    :: basis_recs(:)
        class(ori),          intent(inout) :: o
        class(fplane_type),  intent(in)    :: fpl_ref
        type(fplane_type),   intent(inout) :: mean_fpl
        type(fplane_type),   intent(inout) :: basis_fpls(:)
        logical, optional,   intent(in)    :: apply_ctf_amp
        type(kbinterpol) :: kbwin
        complex :: transfer, mean_val, basis_val
        real    :: rotmat(3,3), loc(3), hrow(3), ctfamp
        real    :: wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM)
        integer :: fpllims_pd(3,2), fpllims(3,2), h, k, hp, kp, pf, q, ncomp
        integer :: h_sq, k_max_h, k_lo, k_hi, nyq_disk, nyq_eff, win(2,3)
        logical :: l_apply_ctf_amp, l_conjg
        ! per-sample geometry, so the volume loop can be hoisted out of the (h,k) sweep
        integer,     allocatable :: swin(:,:,:), shp(:), skp(:)
        real,        allocatable :: swx(:,:), swy(:,:), swz(:,:)
        complex,     allocatable :: stf(:)
        logical,     allocatable :: scj(:)
        ! in-bounds samples first, so the volume loops carry no per-tap window test
        integer,     allocatable :: jok(:), jbad(:)
        integer :: exp_lb(3), exp_ub(3), ns_ok, ns_bad, jj
        integer :: ns, nsmax, j
        if( .not. allocated(mean_rec%cmat_exp) )then
            THROW_HARD('expanded mean matrix does not exist; project_fplanes_mean_basis')
        endif
        if( .not. allocated(fpl_ref%cmplx_plane) )then
            THROW_HARD('reference Fourier plane does not exist; project_fplanes_mean_basis')
        endif
        ncomp = size(basis_recs)
        if( size(basis_fpls) < ncomp )then
            THROW_HARD('basis output plane array too small; project_fplanes_mean_basis')
        endif
        do q = 1, ncomp
            if( .not. allocated(basis_recs(q)%cmat_exp) )then
                THROW_HARD('expanded basis matrix does not exist; project_fplanes_mean_basis')
            endif
        end do
        l_apply_ctf_amp = .false.
        if( present(apply_ctf_amp) ) l_apply_ctf_amp = apply_ctf_amp
        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        call ensure_latent_projection_plane(fpl_ref, mean_fpl)
        do q = 1, ncomp
            call ensure_latent_projection_plane(fpl_ref, basis_fpls(q))
        end do
        rotmat      = o%get_mat()
        pf          = OSMPL_PAD_FAC
        fpllims_pd  = fpl_ref%frlims
        fpllims     = fpllims_pd
        fpllims(1,1)= ceil_div (fpllims_pd(1,1), pf)
        fpllims(1,2)= floor_div(fpllims_pd(1,2), pf)
        fpllims(2,1)= ceil_div (fpllims_pd(2,1), pf)
        fpllims(2,2)= floor_div(fpllims_pd(2,2), pf)
        nyq_eff = mean_rec%get_lfny(1)
        if( fpl_ref%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpl_ref%nyq / pf))
        nyq_disk = nyq_eff * (nyq_eff + 1)
        ! The sample geometry -- location, KB window, interpolation weights, CTF transfer -- depends on
        ! (h,k) and the orientation ALONE; the ncomp+1 volumes differ only in what is read through it.
        ! Interleaving volumes inside the (h,k) loop leaves none of them resident, so essentially every
        ! gather is a cold miss: build the sample list once, then hoist the volume loop outside it.
        ! Bit-exact: every output element is an independent expression of its own sample and volume.
        nsmax = (fpllims(1,2) - fpllims(1,1) + 1) * (nyq_eff + 1)
        allocate(swin(2,3,nsmax), swx(LATENT_WDIM,nsmax), swy(LATENT_WDIM,nsmax), &
            &swz(LATENT_WDIM,nsmax), stf(nsmax), shp(nsmax), skp(nsmax), scj(nsmax))
        ns = 0
        do h = fpllims(1,1), fpllims(1,2)
            h_sq = h*h
            if( h_sq > nyq_disk ) cycle
            k_max_h = int(sqrt(real(nyq_disk - h_sq)))
            k_lo    = max(fpllims(2,1), -k_max_h)
            k_hi    = min(0, min(fpllims(2,2), k_max_h))
            hp      = h * pf
            hrow(1) = real(h) * rotmat(1,1)
            hrow(2) = real(h) * rotmat(1,2)
            hrow(3) = real(h) * rotmat(1,3)
            do k = k_lo, k_hi
                kp     = k * pf
                loc(1) = hrow(1) + real(k) * rotmat(2,1)
                loc(2) = hrow(2) + real(k) * rotmat(2,2)
                loc(3) = hrow(3) + real(k) * rotmat(2,3)
                l_conjg = loc(1) < 0.
                if( l_conjg ) loc = -loc
                call latent_projection_weights(kbwin, loc, win, wx, wy, wz)
                transfer = cmplx(1., 0.)
                if( l_apply_ctf_amp )then
                    if( allocated(fpl_ref%transfer_plane) )then
                        transfer = fpl_ref%transfer_plane(hp,kp)
                    else
                        ctfamp   = sqrt(max(0., fpl_ref%ctfsq_plane(hp,kp)))
                        transfer = cmplx(ctfamp, 0.)
                    endif
                endif
                ns = ns + 1
                swin(:,:,ns) = win
                swx(:,ns)    = wx
                swy(:,ns)    = wy
                swz(:,ns)    = wz
                stf(ns)      = transfer
                shp(ns)      = hp
                skp(ns)      = kp
                scj(ns)      = l_conjg
            end do
        end do
        ! The window-in-lattice test depends on the SAMPLE alone, not on the volume: every
        ! reconstructor here shares one expanded lattice (the same assumption
        ! insert_planes_oversamp_coupled_batch_scaled makes when it reads recs(1)'s bounds for all).
        ! Testing it inside the volume loop re-evaluates lbound/ubound and two any() temporaries
        ! ncomp+1 times per sample; partitioning it out once is bit-identical -- an out-of-bounds
        ! sample contributed CMPLX_ZERO before and is written as zero below.
        exp_lb = lbound(mean_rec%cmat_exp)
        exp_ub = ubound(mean_rec%cmat_exp)
        allocate(jok(ns), jbad(ns))
        ns_ok = 0; ns_bad = 0
        do j = 1, ns
            if( any(swin(1,:,j) < exp_lb) .or. any(swin(2,:,j) > exp_ub) )then
                ns_bad = ns_bad + 1
                jbad(ns_bad) = j
            else
                ns_ok = ns_ok + 1
                jok(ns_ok) = j
            endif
        end do
        do jj = 1, ns_ok
            j = jok(jj)
            mean_val = weighted_expanded_cmat(mean_rec, swin(:,:,j), swx(:,j), swy(:,j), swz(:,j))
            if( scj(j) ) mean_val = conjg(mean_val)
            mean_fpl%cmplx_plane(shp(j),skp(j)) = stf(j) * mean_val
        end do
        do jj = 1, ns_bad
            j = jbad(jj)
            mean_fpl%cmplx_plane(shp(j),skp(j)) = CMPLX_ZERO
        end do
        do q = 1, ncomp
            do jj = 1, ns_ok
                j = jok(jj)
                basis_val = weighted_expanded_cmat(basis_recs(q), swin(:,:,j), swx(:,j), swy(:,j), swz(:,j))
                if( scj(j) ) basis_val = conjg(basis_val)
                basis_fpls(q)%cmplx_plane(shp(j),skp(j)) = stf(j) * basis_val
            end do
            do jj = 1, ns_bad
                j = jbad(jj)
                basis_fpls(q)%cmplx_plane(shp(j),skp(j)) = CMPLX_ZERO
            end do
        end do
        deallocate(swin, swx, swy, swz, stf, shp, skp, scj, jok, jbad)

    end subroutine project_fplanes_mean_basis

    subroutine ensure_latent_projection_plane( fpl_in, fpl_out )
        type(fplane_type), intent(in)    :: fpl_in
        type(fplane_type), intent(inout) :: fpl_out
        logical :: l_realloc
        l_realloc = .not. allocated(fpl_out%cmplx_plane)
        if( .not. l_realloc )then
            ! nyq is part of the test, not just the bounds: the projection writes only inside the disc
            ! that nyq defines, and every consumer reads a disc derived from the same nyq. While the
            ! geometry is unchanged the written set is identical for every particle, so the out-of-disc
            ! remainder keeps the zeros it was given at allocation, and re-zeroing the whole plane once
            ! per particle per basis volume is pure memory traffic -- at d_tilde=128, the plane, 129
            ! times, for every particle. Worth ~35 % of the projection stage.
            l_realloc = any(lbound(fpl_out%cmplx_plane) /= lbound(fpl_in%cmplx_plane)) .or. &
                &any(ubound(fpl_out%cmplx_plane) /= ubound(fpl_in%cmplx_plane)) .or. &
                &fpl_out%nyq /= fpl_in%nyq
        endif
        if( l_realloc )then
            if( allocated(fpl_out%cmplx_plane) ) deallocate(fpl_out%cmplx_plane)
            allocate(fpl_out%cmplx_plane(lbound(fpl_in%cmplx_plane,1):ubound(fpl_in%cmplx_plane,1), &
                &lbound(fpl_in%cmplx_plane,2):ubound(fpl_in%cmplx_plane,2)))
            fpl_out%cmplx_plane = CMPLX_ZERO
        endif
        if( allocated(fpl_out%ctfsq_plane) ) deallocate(fpl_out%ctfsq_plane)
        if( allocated(fpl_out%transfer_plane) ) deallocate(fpl_out%transfer_plane)
        fpl_out%frlims  = fpl_in%frlims
        fpl_out%shconst = fpl_in%shconst
        fpl_out%nyq     = fpl_in%nyq
    end subroutine ensure_latent_projection_plane

    pure subroutine latent_projection_weights( kbwin, loc, win, wx, wy, wz )
        type(kbinterpol), intent(in)  :: kbwin
        real,             intent(in)  :: loc(3)
        integer,          intent(out) :: win(2,3)
        real,             intent(out) :: wx(:), wy(:), wz(:)
        integer :: i, iwinsz, win_lo(3)
        real    :: base(3), ww3(3), sx, sy, sz, inv_wdim, eps_norm
        iwinsz   = ceiling(KBWINSZ - 0.5)
        win(1,:) = nint(loc)
        win(2,:) = win(1,:) + iwinsz
        win(1,:) = win(1,:) - iwinsz
        win_lo   = win(1,:)
        base     = real(win_lo) - loc
        do i = 1, LATENT_WDIM
            ww3   = kbwin%apod(base + real(i-1))
            wx(i) = ww3(1)
            wy(i) = ww3(2)
            wz(i) = ww3(3)
        end do
        sx        = sum(wx)
        sy        = sum(wy)
        sz        = sum(wz)
        inv_wdim  = 1.0 / real(LATENT_WDIM)
        eps_norm  = epsilon(1.0)
        if( abs(sx) > eps_norm )then
            wx = wx * (1.0 / sx)
        else
            wx = inv_wdim
        endif
        if( abs(sy) > eps_norm )then
            wy = wy * (1.0 / sy)
        else
            wy = inv_wdim
        endif
        if( abs(sz) > eps_norm )then
            wz = wz * (1.0 / sz)
        else
            wz = inv_wdim
        endif
    end subroutine latent_projection_weights

    ! Out-of-lattice windows return zero. Keep this test INSIDE: bounds are per volume, and callers
    ! pass volumes (utilde) that need not share mean_rec's lattice, so hoisting it to the caller and
    ! testing once against mean_rec is an out-of-bounds read.
    pure function weighted_expanded_cmat( rec, win, wx, wy, wz ) result( val )
        type(reconstructor), intent(in) :: rec
        integer,             intent(in) :: win(2,3)
        real,                intent(in) :: wx(:), wy(:), wz(:)
        complex :: val
        integer :: ix, iy, iz, hx, ky, mz
        real    :: wyz
        val = CMPLX_ZERO
        do iz = 1, LATENT_WDIM
            mz = win(1,3) + iz - 1
            do iy = 1, LATENT_WDIM
                ky  = win(1,2) + iy - 1
                wyz = wy(iy) * wz(iz)
                do ix = 1, LATENT_WDIM
                    hx  = win(1,1) + ix - 1
                    val = val + rec%cmat_exp(hx,ky,mz) * (wx(ix) * wyz)
                end do
            end do
        end do
    end function weighted_expanded_cmat




    subroutine solve_coupled_basis_exp( basis_recs, rho_cross_exp, ncomp )
        integer,             intent(in)    :: ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real,                intent(in)    :: rho_cross_exp(:,:,:,:)
        complex(dp) :: rhs(ncomp), sol(ncomp)
        real(dp)    :: amat(ncomp,ncomp)
        real(dp)    :: diag_sum, diag_max, ridge, denom
        integer     :: lb(3), ub(3), h, k, m, ih, ik, im, q, r, flag, shell, nyq
        logical     :: diagonal_density
        ! Same shape convention insert_planes_oversamp_coupled_batch_scaled uses to pick its
        ! accumulation mode: a leading extent of ncomp means only the diagonal of the coupled normal
        ! matrix was accumulated, so the per-voxel system decouples into ncomp scalar divisions.
        diagonal_density = size(rho_cross_exp,1) == ncomp .and. ncomp /= (ncomp*(ncomp+1))/2
        lb = lbound(basis_recs(1)%cmat_exp)
        ub = ubound(basis_recs(1)%cmat_exp)
        nyq = basis_recs(1)%get_lfny(1)
        !$omp parallel do collapse(3) default(shared) schedule(static) &
        !$omp private(h,k,m,ih,ik,im,q,r,amat,rhs,sol,diag_sum,diag_max,ridge,denom,flag,shell) proc_bind(close)
        do m = lb(3), ub(3)
            do k = lb(2), ub(2)
                do h = lb(1), ub(1)
                    ih = h - lb(1) + 1
                    ik = k - lb(2) + 1
                    im = m - lb(3) + 1
                    ! Match reconstructor%sampl_dens_correct: values outside
                    ! the spherical Nyquist support are not reconstructable,
                    ! even though they lie inside the Cartesian FFT cube.
                    shell = nint(sqrt(real(h*h+k*k+m*m)))
                    if( shell > nyq )then
                        do q = 1, ncomp
                            basis_recs(q)%cmat_exp(h,k,m) = CMPLX_ZERO
                        end do
                        cycle
                    endif
                    rhs  = DCMPLX_ZERO
                    diag_sum = 0.d0
                    diag_max = 0.d0
                    do q = 1, ncomp
                        rhs(q) = cmplx(basis_recs(q)%cmat_exp(h,k,m), kind=dp)
                        if( diagonal_density )then
                            denom = max(0.d0,real(rho_cross_exp(q,ih,ik,im),dp))
                        else
                            denom = max(0.d0,real(rho_cross_exp(pair_index(q,q),ih,ik,im),dp))
                        endif
                        diag_sum = diag_sum + denom
                        diag_max = max(diag_max,denom)
                    end do
                    if( diag_max <= COUPLED_DENSITY_FLOOR )then
                        do q = 1, ncomp
                            basis_recs(q)%cmat_exp(h,k,m) = CMPLX_ZERO
                        end do
                        cycle
                    endif
                    ridge = COUPLED_MSTEP_RIDGE_REL * diag_sum / real(max(1,ncomp), dp)
                    if( diagonal_density )then
                        ! No off-diagonals were accumulated, so the normal matrix IS its diagonal and
                        ! the ncomp x ncomp Cholesky collapses to ncomp divisions. Same ridge, same
                        ! floor, same fallback -- only the cross terms between components are dropped.
                        do q = 1, ncomp
                            denom = max(0.d0,real(rho_cross_exp(q,ih,ik,im),dp)) + ridge
                            if( denom > DTINY )then
                                sol(q) = rhs(q) / denom
                            else
                                sol(q) = DCMPLX_ZERO
                            endif
                        end do
                        do q = 1, ncomp
                            basis_recs(q)%cmat_exp(h,k,m) = cmplx(real(sol(q), sp), real(aimag(sol(q)), sp))
                        end do
                        cycle
                    endif
                    amat = 0.d0
                    do q = 1, ncomp
                        do r = q, ncomp
                            amat(q,r) = real(rho_cross_exp(pair_index(q,r),ih,ik,im), dp)
                            amat(r,q) = amat(q,r)
                        end do
                    end do
                    do q = 1, ncomp
                        amat(q,q) = amat(q,q) + ridge
                    end do
                    call solve_real_spd_complex(amat, rhs, sol, ncomp, flag)
                    if( flag /= 0 )then
                        do q = 1, ncomp
                            denom = max(abs(amat(q,q)), ridge)
                            if( denom > DTINY )then
                                sol(q) = rhs(q) / denom
                            else
                                sol(q) = DCMPLX_ZERO
                            endif
                        end do
                    endif
                    do q = 1, ncomp
                        basis_recs(q)%cmat_exp(h,k,m) = cmplx(real(sol(q), sp), real(aimag(sol(q)), sp))
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine solve_coupled_basis_exp

    integer pure function pair_index( q, r ) result( ipair )
        integer, intent(in) :: q, r
        ipair = (r * (r - 1)) / 2 + q
    end function pair_index

    !> The flex analog of simple_reconstructor::add_invtausq2rho (the ml_reg SSNR ridge), for the
    !! packed coupled normal matrix: adds the per-component, per-shell inverse prior variance
    !! invtau2(q,sh) -- from simple_flex_pca_crossfsc::crossfsc_to_invtau2 -- to the DIAGONAL rows
    !! pair_index(q,q) of rho_cross_exp, the exact analog of `self%rho(phys) = self%rho(phys) +
    !! invtau2` (simple_reconstructor.f90:1151). Because rho stays per-voxel and invtau2 is
    !! per-shell, solve_coupled_basis_exp then sees H_q(v) + invtau2_q(sh) in the denominator: the
    !! sampling-aware Gilles-Singer S.11 shrinkage H/(H+R), inherited for free. Off-diagonal rows
    !! are untouched -- the prior imposed is independent per component per shell (a diagonal tau^2
    !! prior), the only form a per-component cross-fit FSC curve can inform.
    !! Deliberately a mutate-rho routine rather than an optional argument threaded into
    !! solve_coupled_basis_exp: it matches the precedent's semantics, keeps the solve signature
    !! stable for the GPU/branch surface, and lets the dead-voxel floor (COUPLED_DENSITY_FLOOR) and
    !! the Cholesky-failure fallback see the regularized diagonal. The tiny relative ridge
    !! (COUPLED_MSTEP_RIDGE_REL) stays: it is a conditioning floor with a different job, invisible
    !! at 1e-8 next to any real invtau2.
    !! Call it AFTER any distributed reduction (part files on disk are never mutated) and BEFORE
    !! the solve. Shells below the conversion's k_lo carry invtau2 = 0, so "no addition at very low
    !! resolution" holds by construction; shells above nyq / the invtau2 band are skipped.
    subroutine add_invtausq2rho_coupled( basis_recs, rho_cross_exp, ncomp, invtau2 )
        integer,             intent(in)    :: ncomp
        type(reconstructor), intent(in)    :: basis_recs(ncomp)  !< lattice geometry reference only
        real,                intent(inout) :: rho_cross_exp(:,:,:,:)
        real,                intent(in)    :: invtau2(:,:)       !< (ncomp, nshells) per component, per shell
        integer :: lb(3), ub(3), h, k, m, ih, ik, im, q, shell, nyq, shmax
        if( size(rho_cross_exp,1) /= (ncomp*(ncomp+1))/2 ) &
            &THROW_HARD('add_invtausq2rho_coupled requires the FULL packed coupled density')
        if( size(invtau2,1) < ncomp ) THROW_HARD('add_invtausq2rho_coupled: invtau2 rank mismatch')
        lb    = lbound(basis_recs(1)%cmat_exp)
        ub    = ubound(basis_recs(1)%cmat_exp)
        nyq   = basis_recs(1)%get_lfny(1)
        shmax = min(nyq, size(invtau2,2))
        !$omp parallel do collapse(3) default(shared) schedule(static) &
        !$omp private(h,k,m,ih,ik,im,q,shell) proc_bind(close)
        do m = lb(3), ub(3)
            do k = lb(2), ub(2)
                do h = lb(1), ub(1)
                    ih = h - lb(1) + 1
                    ik = k - lb(2) + 1
                    im = m - lb(3) + 1
                    ! same shell convention and spherical-Nyquist support rule as the solve
                    shell = nint(sqrt(real(h*h+k*k+m*m)))
                    if( shell < 1 .or. shell > shmax ) cycle
                    do q = 1, ncomp
                        rho_cross_exp(pair_index(q,q),ih,ik,im) = &
                            &rho_cross_exp(pair_index(q,q),ih,ik,im) + invtau2(q,shell)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine add_invtausq2rho_coupled

    subroutine solve_real_spd_complex( amat_in, rhs, sol, n, flag )
        integer,     intent(in)  :: n
        real(dp),    intent(in)  :: amat_in(n,n)
        complex(dp), intent(in)  :: rhs(n)
        complex(dp), intent(out) :: sol(n)
        integer,     intent(out) :: flag
        real(dp) :: chol(n,n), yr(n), yi(n), xr(n), xi(n)
        real(dp) :: sumr, sumi, sumv, tol
        integer  :: i, j, l
        flag = 0
        sol  = DCMPLX_ZERO
        chol = 0.d0
        tol  = max(DTINY, epsilon(1.d0) * max(1.d0, maxval(abs(amat_in))))
        do j = 1, n
            sumv = amat_in(j,j)
            do l = 1, j - 1
                sumv = sumv - chol(j,l) * chol(j,l)
            end do
            if( sumv <= tol )then
                flag = 1
                return
            endif
            chol(j,j) = sqrt(sumv)
            do i = j + 1, n
                sumv = amat_in(i,j)
                do l = 1, j - 1
                    sumv = sumv - chol(i,l) * chol(j,l)
                end do
                chol(i,j) = sumv / chol(j,j)
            end do
        end do
        do i = 1, n
            sumr = real(rhs(i), dp)
            sumi = aimag(rhs(i))
            do l = 1, i - 1
                sumr = sumr - chol(i,l) * yr(l)
                sumi = sumi - chol(i,l) * yi(l)
            end do
            yr(i) = sumr / chol(i,i)
            yi(i) = sumi / chol(i,i)
        end do
        do i = n, 1, -1
            sumr = yr(i)
            sumi = yi(i)
            do l = i + 1, n
                sumr = sumr - chol(l,i) * xr(l)
                sumi = sumi - chol(l,i) * xi(l)
            end do
            xr(i) = sumr / chol(i,i)
            xi(i) = sumi / chol(i,i)
        end do
        do i = 1, n
            sol(i) = cmplx(xr(i), xi(i), kind=dp)
        end do
    end subroutine solve_real_spd_complex

    !!  mskrad (optional, pixels at params%box): when present the particle is soft-masked to
    !!  that radius after noise normalization instead of edge-tapered. This is SIMPLE's
    !!  equivalent of the reference's mask_images_in_H_B/mask_images_in_proj (covariance_estimation
    !!  options, both default True), which masks each image to the projected molecular
    !!  envelope before the covariance accumulation. Solvent outside the particle contributes
    !!  only noise, and at box_crop=64/mskdiam=200 the disc keeps ~21% of the frame, so the
    !!  noise in every per-image inner product drops by roughly the same factor. Left absent
    !!  the behaviour is exactly as before.
    subroutine prep_imgs4projected_model( params, build, nptcls, ptcl_imgs, pinds, fplanes, &
        &mskrad, force_cpu, resident, cached )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls
        class(image),      intent(inout) :: ptcl_imgs(nptcls)
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), intent(inout) :: fplanes(nptcls)
        real, optional,    intent(in)    :: mskrad
        logical, optional, intent(in)    :: force_cpu   !< cross-check reference building
        logical, optional, intent(in)    :: resident    !< device path: leave planes resident only
        logical, optional, intent(in)    :: cached      !< serve reads from the downscaled cache
        type(ctfparams) :: ctfparms(nthr_glob)
        real    :: shift(2), crop_factor
        integer :: iptcl, i, ithr, kfromto(2)
        logical :: l_mask, l_cpu, l_res, l_cached
        l_mask = .false.
        if( present(mskrad) ) l_mask = mskrad > 0.0
        l_cpu = .false.
        if( present(force_cpu) ) l_cpu = force_cpu
        l_res = .false.
        if( present(resident) ) l_res = resident
        l_cached = .false.
        if( present(cached) ) l_cached = cached
        ! A cache entry is the noise-normalised, Fourier-cropped particle at box_crop. That prefix
        ! is equivalent to the full-box path only for the TAPER variant, which is the one
        ! prep_imgs4rec certified: norm_noise_mask_pad_fft has no renorm= switch, so a masked run
        ! would noise-normalise a second time. Refuse rather than change the numerics silently.
        if( l_cached .and. l_mask ) THROW_HARD('particle cache is incompatible with image masking &
            &(COV_MASK_IMAGES); prep_imgs4projected_model')
        ! device path: when a stage driver has begun the GPU prep lifecycle, the whole
        ! taper->norm->pad->FFT->plane chain runs on device and the planes are fetched packed
        ! (taper variant only; the mask variant stays on the CPU)
        if( associated(flex_dev_prep_hook) .and. .not. l_mask .and. .not. l_cpu .and. .not. l_cached )then
            call flex_dev_prep_hook(params, build, nptcls, ptcl_imgs, pinds, &
                &fplanes, fetch=.not. l_res)
            return
        endif
        if( l_res ) THROW_HARD('resident prep requested without the device prep lifecycle')
        ! logical/physical address mapping for padded Fourier planes: a cached particle already
        ! lives on the cropped grid, so the pad heap and the map must both be box_croppd
        if( l_cached )then
            call memoize_ft_maps([params%box_croppd, params%box_croppd, 1], params%smpd_crop)
        else
            call memoize_ft_maps([params%boxpd, params%boxpd, 1], params%smpd)
        endif
        kfromto = projected_model_kfromto(params)
        if( l_cached ) kfromto(2) = min(kfromto(2), params%box_crop/2)
        crop_factor = real(params%box_crop) / real(params%box)
        !$omp parallel do default(shared) private(i,ithr,iptcl,shift) schedule(static) proc_bind(close)
        do i = 1, nptcls
            ithr   = omp_get_thread_num() + 1
            iptcl  = pinds(i)
            if( l_mask )then
                call ptcl_imgs(i)%norm_noise_mask_pad_fft(build%lmsk, mskrad, build%img_pad_heap(ithr))
            else if( l_cached )then
                ! the plane cache holds the full-box prep's padded transform on this grid: load it
                ! as the heap image's transform and continue exactly as the full-box path does
                call plane_cache_fill(i, build%img_pad_heap(ithr))
            else
                call ptcl_imgs(i)%norm_noise_taper_edge_pad_fft(build%lmsk, build%img_pad_heap(ithr))
            endif
            ctfparms(ithr) = build%spproj%get_ctfparams(params%oritype, iptcl)
            shift = build%spproj_field%get_2Dshift(iptcl)
            if( l_cached )then
                ! shconst is in pixels of the padded box the image actually has, and the CTF kernel
                ! reads cycles/pixel of the current grid -- both must move to the cropped grid
                ctfparms(ithr)%smpd = ctfparms(ithr)%smpd / crop_factor   ! = smpd_crop
                shift               = shift * crop_factor
            endif
            if( params%l_ml_reg )then
                if( .not. allocated(build%esig%sigma2_noise) )then
                    THROW_HARD('projected covariance model requested whitening without loaded sigma2 spectra')
                endif
                if( iptcl < lbound(build%esig%sigma2_noise,2) .or. &
                    &iptcl > ubound(build%esig%sigma2_noise,2) )then
                    THROW_HARD('projected covariance particle index is outside the sigma2 table')
                endif
                call build%img_pad_heap(ithr)%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms(ithr), &
                    &shift, fplanes(i), build%esig%sigma2_noise(kfromto(1):kfromto(2),iptcl), &
                    &store_transfer=.true., observation_model=.true.)
            else
                call build%img_pad_heap(ithr)%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms(ithr), &
                    &shift, fplanes(i), store_transfer=.true., observation_model=.true.)
            endif
            call cap_fplane_for_projected_model(fplanes(i), kfromto)
        end do
        !$omp end parallel do
    end subroutine prep_imgs4projected_model

    subroutine cap_fplane_for_projected_model( fpl, kfromto )
        type(fplane_type), intent(inout) :: fpl
        integer,           intent(in)    :: kfromto(2)
        integer :: nyq_eff
        nyq_eff = max(OSMPL_PAD_FAC, OSMPL_PAD_FAC * kfromto(2))
        if( fpl%nyq > 0 ) fpl%nyq = min(fpl%nyq, nyq_eff)
    end subroutine cap_fplane_for_projected_model

    function projected_model_kfromto( params ) result( kfromto )
        class(parameters), intent(in) :: params
        integer :: kfromto(2), kto_full
        real    :: dstep_crop
        kto_full = max(1, fdim(params%box_crop) - 1)
        kfromto(1) = 1
        kfromto(2) = kto_full
        if( params%lp > 2.0 * params%smpd_crop + TINY )then
            dstep_crop = real(max(1, params%box_crop - 1)) * params%smpd_crop
            kfromto(2) = max(1, min(kto_full, int(dstep_crop / params%lp)))
        endif
    end function projected_model_kfromto

    subroutine cleanup_plane( fpl )
        type(fplane_type), intent(inout) :: fpl
        if( allocated(fpl%cmplx_plane) ) deallocate(fpl%cmplx_plane)
        if( allocated(fpl%ctfsq_plane) ) deallocate(fpl%ctfsq_plane)
        if( allocated(fpl%transfer_plane) ) deallocate(fpl%transfer_plane)
        fpl%frlims  = 0
        fpl%shconst = 0.
        fpl%nyq     = 0
    end subroutine cleanup_plane

end module simple_flex_reconstructor_latent_ops
