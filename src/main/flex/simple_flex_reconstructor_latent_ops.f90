!@descr: Latent-volume projection/backprojection helpers for flex_analysis
module simple_flex_reconstructor_latent_ops
use simple_core_module_api
use simple_reconstructor, only: reconstructor
implicit none

public :: insert_plane_oversamp_multi_scaled, insert_plane_oversamp_coupled_scaled
public :: insert_planes_oversamp_coupled_batch_scaled
public :: accumulate_planes_oversamp_coupled_stats_batch
public :: test_coupled_batch_accumulation
public :: test_cartesian_projection_contract
public :: accumulate_plane_oversamp_coupled_stats, project_fplane_mean, project_fplanes_mean_basis
private
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

contains

    subroutine insert_plane_oversamp_multi_scaled( recs, se, o, fpl, data_scales, density_scales )
        use simple_math, only: ceil_div, floor_div
        type(reconstructor), intent(inout) :: recs(:)
        class(sym),          intent(inout) :: se
        class(ori),          intent(inout) :: o
        class(fplane_type),  intent(in)    :: fpl
        real(dp),            intent(in)    :: data_scales(:), density_scales(:)
        type(kbinterpol) :: kbwin
        type(ori) :: o_sym
        complex   :: comp_base, cmplx_raw
        real      :: rotmats(se%get_nsym(),3,3), loc(3), hrow(3), ctfsq_raw
        real      :: wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM), ww
        real      :: data_scale_sp(size(recs)), density_scale_sp(size(recs))
        real      :: r11, r12, r13, r21, r22, r23
        integer   :: win(2, 3), h, k, l, nsym, isym, iwinsz, stride, fpllims_pd(3, 2)
        integer   :: fpllims(3, 2), hp, kp, pf, ix, iy, iz, hx, ky, mz, q, ncomp
        integer   :: nyq_disk, nyq_eff, h_sq, k_max_h, k_lo, k_hi, exp_lb(3), exp_ub(3)
        real      :: pf2, eps_norm, inv_wdim
        ncomp = size(recs)
        if( ncomp <= 0 ) return
        if( size(data_scales) < ncomp .or. size(density_scales) < ncomp )then
            THROW_HARD('scale vector smaller than reconstructor array; insert_plane_oversamp_multi_scaled')
        endif
        if( .not. allocated(recs(1)%cmat_exp) )then
            THROW_HARD('expanded matrix does not exist; insert_plane_oversamp_multi_scaled')
        endif
        kbwin    = kbinterpol(KBWINSZ, KBALPHA)
        iwinsz   = ceiling(KBWINSZ - 0.5)
        stride   = LATENT_SAFE_STRIDE
        exp_lb   = lbound(recs(1)%cmat_exp)
        exp_ub   = ubound(recs(1)%cmat_exp)
        nsym     = se%get_nsym()
        rotmats(1,:,:) = o%get_mat()
        if( nsym > 1 )then
            do isym = 2, nsym
                call se%apply(o, isym, o_sym)
                rotmats(isym,:,:) = o_sym%get_mat()
            end do
        endif
        fpllims_pd   = fpl%frlims
        pf           = OSMPL_PAD_FAC
        pf2          = real(pf*pf)
        fpllims      = fpllims_pd
        fpllims(1,1) = ceil_div (fpllims_pd(1,1), pf)
        fpllims(1,2) = floor_div(fpllims_pd(1,2), pf)
        fpllims(2,1) = ceil_div (fpllims_pd(2,1), pf)
        fpllims(2,2) = floor_div(fpllims_pd(2,2), pf)
        nyq_eff = recs(1)%get_lfny(1)
        if( fpl%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpl%nyq / pf))
        nyq_disk = nyq_eff * (nyq_eff + 1)
        eps_norm = epsilon(1.0)
        inv_wdim = 1.0 / real(LATENT_WDIM)
        do q = 1, ncomp
            data_scale_sp(q)    = real(data_scales(q))
            density_scale_sp(q) = real(max(0.d0, density_scales(q)))
        end do
        !$omp parallel default(shared) private(h,k,l,h_sq,k_max_h,k_lo,k_hi,cmplx_raw,&
        !$omp& ctfsq_raw,comp_base,wx,wy,wz,ww,win,loc,hrow,hp,kp,r11,r12,r13,r21,r22,r23,&
        !$omp& isym,ix,iy,iz,hx,ky,mz,q) proc_bind(close)
        do isym = 1, nsym
            r11 = rotmats(isym,1,1); r12 = rotmats(isym,1,2); r13 = rotmats(isym,1,3)
            r21 = rotmats(isym,2,1); r22 = rotmats(isym,2,2); r23 = rotmats(isym,2,3)
            do l = 0, stride-1
                !$omp do schedule(static,1)
                do h = fpllims(1,1)+l, fpllims(1,2), stride
                    h_sq = h*h
                    if( h_sq > nyq_disk ) cycle
                    k_max_h = int(sqrt(real(nyq_disk - h_sq)))
                    k_lo    = max(fpllims(2,1), -k_max_h)
                    k_hi    = min(fpllims(2,2),  k_max_h)
                    hp      = h * pf
                    hrow(1) = real(h) * r11
                    hrow(2) = real(h) * r12
                    hrow(3) = real(h) * r13
                    do k = k_lo, k_hi
                        kp = k * pf
                        if( kp <= 0 )then
                            cmplx_raw = fpl%cmplx_plane(hp,kp)
                            ctfsq_raw = fpl%ctfsq_plane(hp,kp)
                        else
                            cmplx_raw = conjg(fpl%cmplx_plane(-hp,-kp))
                            ctfsq_raw = fpl%ctfsq_plane(-hp,-kp)
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
                        call kb_apod_vecs_3d_fast(loc, wx, wy, wz)
                        do iz = 1, LATENT_WDIM
                            mz = win(1,3) + iz - 1
                            do iy = 1, LATENT_WDIM
                                ky = win(1,2) + iy - 1
                                do ix = 1, LATENT_WDIM
                                    hx = win(1,1) + ix - 1
                                    ww = wx(ix) * (wy(iy) * wz(iz))
                                    do q = 1, ncomp
                                        recs(q)%cmat_exp(hx,ky,mz) = recs(q)%cmat_exp(hx,ky,mz) + &
                                            &(data_scale_sp(q) * comp_base) * ww
                                        recs(q)%rho_exp(hx,ky,mz) = recs(q)%rho_exp(hx,ky,mz) + &
                                            &(density_scale_sp(q) * ctfsq_raw) * ww
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
                !$omp end do
            end do
        end do
        !$omp end parallel
        call o_sym%kill

    contains

        subroutine kb_apod_vecs_3d_fast( loc, wx, wy, wz )
            real, intent(in)  :: loc(3)
            real, intent(out) :: wx(:), wy(:), wz(:)
            integer :: i, win_lo(3)
            real    :: base(3), ww3(3), sx, sy, sz
            win_lo = nint(loc) - iwinsz
            base   = real(win_lo) - loc
            do i = 1, LATENT_WDIM
                ww3   = kbwin%apod_fast(base + real(i-1))
                wx(i) = ww3(1)
                wy(i) = ww3(2)
                wz(i) = ww3(3)
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
        end subroutine kb_apod_vecs_3d_fast

    end subroutine insert_plane_oversamp_multi_scaled

    subroutine insert_plane_oversamp_coupled_scaled( recs, rho_cross_exp, se, o, fpl, data_scales, density_scales )
        use simple_math, only: ceil_div, floor_div
        type(reconstructor), intent(inout) :: recs(:)
        real,                intent(inout) :: rho_cross_exp(:,:,:,:)
        class(sym),          intent(inout) :: se
        class(ori),          intent(inout) :: o
        class(fplane_type),  intent(in)    :: fpl
        real(dp),            intent(in)    :: data_scales(:), density_scales(:,:)
        type(kbinterpol) :: kbwin
        type(ori) :: o_sym
        complex   :: comp_base, cmplx_raw
        real      :: rotmats(se%get_nsym(),3,3), loc(3), hrow(3), ctfsq_raw
        real      :: wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM), ww
        real      :: data_scale_sp(size(recs)), density_scale_packed((size(recs)*(size(recs)+1))/2)
        real      :: r11, r12, r13, r21, r22, r23
        integer   :: win(2, 3), h, k, l, nsym, isym, iwinsz, stride, fpllims_pd(3, 2)
        integer   :: fpllims(3, 2), hp, kp, pf, ix, iy, iz, hx, ky, mz, q, r, ncomp, ipair
        integer   :: nyq_disk, nyq_eff, h_sq, k_max_h, k_lo, k_hi, ih, ik, im
        integer   :: exp_lb(3), exp_ub(3), exp_shape(3), npairs
        real      :: pf2, eps_norm, inv_wdim
        ncomp = size(recs)
        if( ncomp <= 0 ) return
        npairs = (ncomp * (ncomp + 1)) / 2
        if( size(data_scales) < ncomp .or. size(density_scales,1) < ncomp .or. size(density_scales,2) < ncomp )then
            THROW_HARD('scale array smaller than reconstructor array; insert_plane_oversamp_coupled_scaled')
        endif
        if( .not. allocated(recs(1)%cmat_exp) )then
            THROW_HARD('expanded matrix does not exist; insert_plane_oversamp_coupled_scaled')
        endif
        if( .not. allocated(fpl%transfer_plane) )then
            THROW_HARD('forward transfer plane does not exist; insert_plane_oversamp_coupled_scaled')
        endif
        exp_lb    = lbound(recs(1)%cmat_exp)
        exp_ub    = ubound(recs(1)%cmat_exp)
        exp_shape = shape(recs(1)%cmat_exp)
        if( size(rho_cross_exp,1) < npairs .or. size(rho_cross_exp,2) < exp_shape(1) .or. &
            &size(rho_cross_exp,3) < exp_shape(2) .or. size(rho_cross_exp,4) < exp_shape(3) )then
            THROW_HARD('cross-density array shape mismatch; insert_plane_oversamp_coupled_scaled')
        endif
        kbwin    = kbinterpol(KBWINSZ, KBALPHA)
        iwinsz   = ceiling(KBWINSZ - 0.5)
        stride   = LATENT_SAFE_STRIDE
        nsym     = se%get_nsym()
        rotmats(1,:,:) = o%get_mat()
        if( nsym > 1 )then
            do isym = 2, nsym
                call se%apply(o, isym, o_sym)
                rotmats(isym,:,:) = o_sym%get_mat()
            end do
        endif
        fpllims_pd   = fpl%frlims
        pf           = OSMPL_PAD_FAC
        pf2          = real(pf*pf)
        fpllims      = fpllims_pd
        fpllims(1,1) = ceil_div (fpllims_pd(1,1), pf)
        fpllims(1,2) = floor_div(fpllims_pd(1,2), pf)
        fpllims(2,1) = ceil_div (fpllims_pd(2,1), pf)
        fpllims(2,2) = floor_div(fpllims_pd(2,2), pf)
        nyq_eff = recs(1)%get_lfny(1)
        if( fpl%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpl%nyq / pf))
        nyq_disk = nyq_eff * (nyq_eff + 1)
        eps_norm = epsilon(1.0)
        inv_wdim = 1.0 / real(LATENT_WDIM)
        do q = 1, ncomp
            data_scale_sp(q) = real(data_scales(q))
        end do
        do r = 1, ncomp
            do q = 1, r
                density_scale_packed(pair_index(q,r)) = real(density_scales(q,r))
            end do
        end do
        !$omp parallel default(shared) private(h,k,l,h_sq,k_max_h,k_lo,k_hi,cmplx_raw,&
        !$omp& ctfsq_raw,comp_base,wx,wy,wz,ww,win,loc,hrow,hp,kp,r11,r12,r13,r21,r22,r23,&
        !$omp& isym,ix,iy,iz,hx,ky,mz,ih,ik,im,q,r,ipair) proc_bind(close)
        do isym = 1, nsym
            r11 = rotmats(isym,1,1); r12 = rotmats(isym,1,2); r13 = rotmats(isym,1,3)
            r21 = rotmats(isym,2,1); r22 = rotmats(isym,2,2); r23 = rotmats(isym,2,3)
            do l = 0, stride-1
                !$omp do schedule(static,1)
                do h = fpllims(1,1)+l, fpllims(1,2), stride
                    h_sq = h*h
                    if( h_sq > nyq_disk ) cycle
                    k_max_h = int(sqrt(real(nyq_disk - h_sq)))
                    k_lo    = max(fpllims(2,1), -k_max_h)
                    k_hi    = min(fpllims(2,2),  k_max_h)
                    hp      = h * pf
                    hrow(1) = real(h) * r11
                    hrow(2) = real(h) * r12
                    hrow(3) = real(h) * r13
                    do k = k_lo, k_hi
                        kp = k * pf
                        if( kp <= 0 )then
                            cmplx_raw = conjg(fpl%transfer_plane(hp,kp)) * fpl%cmplx_plane(hp,kp)
                            ctfsq_raw = fpl%ctfsq_plane(hp,kp)
                        else
                            cmplx_raw = conjg(conjg(fpl%transfer_plane(-hp,-kp)) * fpl%cmplx_plane(-hp,-kp))
                            ctfsq_raw = fpl%ctfsq_plane(-hp,-kp)
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
                        call kb_apod_vecs_3d_fast(loc, wx, wy, wz)
                        do iz = 1, LATENT_WDIM
                            mz = win(1,3) + iz - 1
                            im = mz - exp_lb(3) + 1
                            do iy = 1, LATENT_WDIM
                                ky = win(1,2) + iy - 1
                                ik = ky - exp_lb(2) + 1
                                do ix = 1, LATENT_WDIM
                                    hx = win(1,1) + ix - 1
                                    ih = hx - exp_lb(1) + 1
                                    ww = wx(ix) * (wy(iy) * wz(iz))
                                    do q = 1, ncomp
                                        recs(q)%cmat_exp(hx,ky,mz) = recs(q)%cmat_exp(hx,ky,mz) + &
                                            &(data_scale_sp(q) * comp_base) * ww
                                    end do
                                    !$omp simd
                                    do ipair = 1, npairs
                                        rho_cross_exp(ipair,ih,ik,im) = rho_cross_exp(ipair,ih,ik,im) + &
                                            &(density_scale_packed(ipair) * ctfsq_raw) * ww
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
                !$omp end do
            end do
        end do
        !$omp end parallel
        call o_sym%kill

    contains

        integer pure function pair_index( q, r ) result( ipair )
            integer, intent(in) :: q, r
            ipair = (r * (r - 1)) / 2 + q
        end function pair_index

        subroutine kb_apod_vecs_3d_fast( loc, wx, wy, wz )
            real, intent(in)  :: loc(3)
            real, intent(out) :: wx(:), wy(:), wz(:)
            integer :: i, win_lo(3)
            real    :: base(3), ww3(3), sx, sy, sz
            win_lo = nint(loc) - iwinsz
            base   = real(win_lo) - loc
            do i = 1, LATENT_WDIM
                ww3   = kbwin%apod_fast(base + real(i-1))
                wx(i) = ww3(1)
                wy(i) = ww3(2)
                wz(i) = ww3(3)
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
        end subroutine kb_apod_vecs_3d_fast

    end subroutine insert_plane_oversamp_coupled_scaled

    subroutine insert_planes_oversamp_coupled_batch_scaled( recs, rho_cross_exp, se, orientations, fpls, &
        &data_scales, density_scales, valid, nrecords )
        use simple_math, only: ceil_div, floor_div
        use simple_kbinterpol, only: apod_kb15_a2
        type(reconstructor), intent(inout) :: recs(:)
        real,                intent(inout) :: rho_cross_exp(:,:,:,:)
        class(sym),          intent(inout) :: se
        type(ori),           intent(inout) :: orientations(:)
        type(fplane_type),   intent(in)    :: fpls(:)
        real(dp),            intent(in)    :: data_scales(:,:), density_scales(:,:,:)
        logical,             intent(in)    :: valid(:)
        integer,             intent(in)    :: nrecords
        type(ori) :: o_sym
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
        real      :: pf2, eps_norm, inv_wdim
        ncomp = size(recs)
        if( ncomp <= 0 .or. nrecords <= 0 ) return
        npairs = (ncomp * (ncomp + 1)) / 2
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
        if( size(rho_cross_exp,1)<npairs .or. size(rho_cross_exp,2)<exp_shape(1) .or. &
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
        allocate(rotmats(3,3,nsym,nrecords), data_scale_sp(ncomp,nrecords), &
            &density_scale_packed(npairs,nrecords), source=0.)
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
            do r = 1, ncomp
                do q = 1, r
                    density_scale_packed(pair_index(q,r),i) = real(density_scales(q,r,i))
                end do
            end do
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
                                        !$omp simd
                                        do ipair = 1, npairs
                                            rho_cross_exp(ipair,ih,ik,im) = rho_cross_exp(ipair,ih,ik,im) + &
                                                &(density_scale_packed(ipair,i)*ctfsq_raw)*ww
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
        deallocate(rotmats,data_scale_sp,density_scale_packed,fpllims,nyq_disks)

    contains

        integer pure function pair_index( q, r ) result( ipair )
            integer, intent(in) :: q, r
            ipair = (r*(r-1))/2+q
        end function pair_index

        subroutine kb_apod_vecs_3d_fast( loc, wx, wy, wz )
            real, intent(in)  :: loc(3)
            real, intent(out) :: wx(:), wy(:), wz(:)
            integer :: j, win_lo(3)
            real :: base(3), sx, sy, sz
            win_lo = nint(loc)-iwinsz
            base = real(win_lo)-loc
            do j = 1, LATENT_WDIM
                wx(j)=apod_kb15_a2(base(1)+real(j-1))
                wy(j)=apod_kb15_a2(base(2)+real(j-1))
                wz(j)=apod_kb15_a2(base(3)+real(j-1))
            end do
            sx=sum(wx); sy=sum(wy); sz=sum(wz)
            if( abs(sx)>eps_norm )then; wx=wx/sx; else; wx=inv_wdim; endif
            if( abs(sy)>eps_norm )then; wy=wy/sy; else; wy=inv_wdim; endif
            if( abs(sz)>eps_norm )then; wz=wz/sz; else; wz=inv_wdim; endif
        end subroutine kb_apod_vecs_3d_fast

    end subroutine insert_planes_oversamp_coupled_batch_scaled

    subroutine test_coupled_batch_accumulation()
        use simple_ftiter, only: ftiter
        integer, parameter :: TEST_BOX=16, TEST_NCOMP=3, TEST_NRECORDS=2
        type(reconstructor) :: rec_single(TEST_NCOMP),rec_batch(TEST_NCOMP)
        type(fplane_type) :: fpls(TEST_NRECORDS)
        type(ori) :: orientations(TEST_NRECORDS)
        type(sym) :: se
        type(ftiter) :: fit_pd
        real(dp) :: z(TEST_NCOMP,TEST_NRECORDS),zz(TEST_NCOMP,TEST_NCOMP,TEST_NRECORDS)
        complex, allocatable :: rhs_stats_single(:,:,:,:),rhs_stats_batch(:,:,:,:)
        real, allocatable :: rho_single(:,:,:,:),rho_batch(:,:,:,:),rho_stats_single(:,:,:,:),rho_stats_batch(:,:,:,:)
        logical :: valid(TEST_NRECORDS)
        integer :: lims_pd(3,2),exp_lb(3),exp_ub(3),exp_shape(3)
        integer :: i,q,r,h,k,hp,kp,npairs
        real :: rhs_err,rho_err
        call se%new('c1')
        fit_pd=ftiter([OSMPL_PAD_FAC*TEST_BOX,OSMPL_PAD_FAC*TEST_BOX,1],1.)
        lims_pd=fit_pd%loop_lims(3)
        exp_lb=[-TEST_BOX/2-2,-TEST_BOX/2-2,-TEST_BOX/2-2]
        exp_ub=[ TEST_BOX/2+2, TEST_BOX/2+2, TEST_BOX/2+2]
        exp_shape=exp_ub-exp_lb+1
        npairs=(TEST_NCOMP*(TEST_NCOMP+1))/2
        do q=1,TEST_NCOMP
            call rec_single(q)%new([TEST_BOX,TEST_BOX,TEST_BOX],1.)
            call rec_batch(q)%new([TEST_BOX,TEST_BOX,TEST_BOX],1.)
            allocate(rec_single(q)%cmat_exp(exp_lb(1):exp_ub(1),exp_lb(2):exp_ub(2),exp_lb(3):exp_ub(3)), &
                &source=CMPLX_ZERO)
            allocate(rec_batch(q)%cmat_exp(exp_lb(1):exp_ub(1),exp_lb(2):exp_ub(2),exp_lb(3):exp_ub(3)), &
                &source=CMPLX_ZERO)
        end do
        do i=1,TEST_NRECORDS
            call orientations(i)%new(.true.)
            call orientations(i)%set_euler([17.*real(i),31.*real(i),11.*real(i)])
            allocate(fpls(i)%cmplx_plane(lims_pd(1,1):lims_pd(1,2),lims_pd(2,1):0),source=CMPLX_ZERO)
            allocate(fpls(i)%ctfsq_plane(lims_pd(1,1):lims_pd(1,2),lims_pd(2,1):0),source=0.)
            allocate(fpls(i)%transfer_plane(lims_pd(1,1):lims_pd(1,2),lims_pd(2,1):0),source=CMPLX_ZERO)
            fpls(i)%frlims=lims_pd
            fpls(i)%nyq=fit_pd%get_lfny(1)
            do k=-TEST_BOX/2,0
                kp=OSMPL_PAD_FAC*k
                do h=-TEST_BOX/2,TEST_BOX/2
                    if( h*h+k*k>(TEST_BOX/2)*(TEST_BOX/2+1) ) cycle
                    hp=OSMPL_PAD_FAC*h
                    fpls(i)%cmplx_plane(hp,kp)=cmplx(0.01*real(3*h-2*k+i),0.02*real(h+k-i))
                    fpls(i)%transfer_plane(hp,kp)=cmplx(0.7+0.01*real(mod(abs(h+2*k+i),7)),0.)
                    fpls(i)%ctfsq_plane(hp,kp)=0.2+0.03*real(mod(abs(2*h-k+i),9))
                end do
            end do
        end do
        z(:,1)=[0.4d0,-0.2d0,0.7d0]
        z(:,2)=[-0.1d0,0.8d0,0.3d0]
        do i=1,TEST_NRECORDS
            do r=1,TEST_NCOMP
                do q=1,TEST_NCOMP
                    zz(q,r,i)=z(q,i)*z(r,i)
                end do
            end do
        end do
        valid=.true.
        allocate(rho_single(npairs,exp_shape(1),exp_shape(2),exp_shape(3)),source=0.)
        allocate(rho_batch(npairs,exp_shape(1),exp_shape(2),exp_shape(3)),source=0.)
        do i=1,TEST_NRECORDS
            call insert_plane_oversamp_coupled_scaled(rec_single,rho_single,se,orientations(i),fpls(i),z(:,i),zz(:,:,i))
        end do
        call insert_planes_oversamp_coupled_batch_scaled(rec_batch,rho_batch,se,orientations,fpls,z,zz,valid,TEST_NRECORDS)
        rhs_err=0.
        do q=1,TEST_NCOMP
            rhs_err=max(rhs_err,maxval(abs(rec_single(q)%cmat_exp-rec_batch(q)%cmat_exp)))
        end do
        rho_err=maxval(abs(rho_single-rho_batch))
        if( rhs_err>2.e-5 .or. rho_err>2.e-5 ) THROW_HARD('coupled batch accumulation differs from scalar path')
        allocate(rhs_stats_single(exp_shape(1),exp_shape(2),exp_shape(3),TEST_NCOMP),source=CMPLX_ZERO)
        allocate(rhs_stats_batch(exp_shape(1),exp_shape(2),exp_shape(3),TEST_NCOMP),source=CMPLX_ZERO)
        allocate(rho_stats_single(npairs,exp_shape(1),exp_shape(2),exp_shape(3)),source=0.)
        allocate(rho_stats_batch(npairs,exp_shape(1),exp_shape(2),exp_shape(3)),source=0.)
        do i=1,TEST_NRECORDS
            call accumulate_plane_oversamp_coupled_stats(rhs_stats_single,rho_stats_single,exp_lb,TEST_BOX/2, &
                &se,orientations(i),fpls(i),z(:,i),zz(:,:,i))
        end do
        call accumulate_planes_oversamp_coupled_stats_batch(rhs_stats_batch,rho_stats_batch,exp_lb,TEST_BOX/2, &
            &se,orientations,fpls,z,zz,valid,TEST_NRECORDS)
        rhs_err=maxval(abs(rhs_stats_single-rhs_stats_batch))
        rho_err=maxval(abs(rho_stats_single-rho_stats_batch))
        if( rhs_err>2.e-5 .or. rho_err>2.e-5 ) THROW_HARD('coupled statistics batch differs from scalar path')
        do i=1,TEST_NRECORDS
            if( allocated(fpls(i)%cmplx_plane) ) deallocate(fpls(i)%cmplx_plane)
            if( allocated(fpls(i)%ctfsq_plane) ) deallocate(fpls(i)%ctfsq_plane)
            if( allocated(fpls(i)%transfer_plane) ) deallocate(fpls(i)%transfer_plane)
            call orientations(i)%kill
        end do
        do q=1,TEST_NCOMP
            call rec_single(q)%kill
            call rec_batch(q)%kill
        end do
        call se%kill
        deallocate(rho_single,rho_batch,rhs_stats_single,rhs_stats_batch,rho_stats_single,rho_stats_batch)
    end subroutine test_coupled_batch_accumulation

    subroutine test_cartesian_projection_contract()
        use simple_ftiter,     only: ftiter
        use simple_parameters, only: parameters
        use simple_sp_project, only: sp_project
        integer, parameter :: TEST_BOX=16, TEST_NBASIS=2
        type(parameters) :: params
        type(sp_project) :: spproj
        type(reconstructor) :: mean_rec, basis_recs(TEST_NBASIS)
        type(reconstructor) :: splat_std, splat_latent(1)
        type(fplane_type) :: fpl_ref, mean_fpl, mean_fpl_std
        type(fplane_type) :: basis_fpls(TEST_NBASIS), basis_fpls_std(TEST_NBASIS)
        type(ori) :: orientation
        type(sym) :: se
        type(ftiter) :: fit_pd
        type(kbinterpol) :: kbwin
        complex :: x(LATENT_WDIM,LATENT_WDIM,LATENT_WDIM)
        complex :: splat(LATENT_WDIM,LATENT_WDIM,LATENT_WDIM), y, lhs, rhs, gathered
        real, allocatable :: rho_cross(:,:,:,:)
        real(dp) :: unit_scale(1), unit_second(1,1)
        real :: wx(LATENT_WDIM),wy(LATENT_WDIM),wz(LATENT_WDIM),loc(3),err
        integer :: lims_pd(3,2),win(2,3),lb(3),ub(3),h,k,m,hp,kp,q,ix,iy,iz
        params%box       = TEST_BOX
        params%box_crop  = TEST_BOX
        params%smpd      = 1.
        params%smpd_crop = 1.
        params%oritype   = 'ptcl3D'
        call spproj%os_stk%new(1,is_ptcl=.false.)
        call spproj%os_stk%set(1,'ctf','no')
        call spproj%os_ptcl3D%new(1,is_ptcl=.true.)
        call spproj%os_ptcl3D%set(1,'stkind',1)
        call spproj%os_ptcl3D%set(1,'indstk',1)
        call mean_rec%new([TEST_BOX,TEST_BOX,TEST_BOX],1.)
        call mean_rec%alloc_rho(params,spproj,expand=.true.)
        do q=1,TEST_NBASIS
            call basis_recs(q)%new([TEST_BOX,TEST_BOX,TEST_BOX],1.)
            call basis_recs(q)%alloc_rho(params,spproj,expand=.true.)
        end do
        lb=lbound(mean_rec%cmat_exp); ub=ubound(mean_rec%cmat_exp)
        do m=lb(3),ub(3)
            do k=lb(2),ub(2)
                do h=lb(1),ub(1)
                    mean_rec%cmat_exp(h,k,m)=cmplx(0.003*real(2*h-k+3*m),0.002*real(h+2*k-m))
                    do q=1,TEST_NBASIS
                        basis_recs(q)%cmat_exp(h,k,m)=cmplx(0.002*real(q*h-2*k+m), &
                            &0.001*real(h+q*k-3*m))
                    end do
                end do
            end do
        end do
        fit_pd=ftiter([OSMPL_PAD_FAC*TEST_BOX,OSMPL_PAD_FAC*TEST_BOX,1],1.)
        lims_pd=fit_pd%loop_lims(3)
        allocate(fpl_ref%cmplx_plane(lims_pd(1,1):lims_pd(1,2),lims_pd(2,1):0),source=CMPLX_ZERO)
        allocate(fpl_ref%ctfsq_plane(lims_pd(1,1):lims_pd(1,2),lims_pd(2,1):0),source=0.)
        allocate(fpl_ref%transfer_plane(lims_pd(1,1):lims_pd(1,2),lims_pd(2,1):0),source=CMPLX_ZERO)
        fpl_ref%frlims=lims_pd
        fpl_ref%nyq=fit_pd%get_lfny(1)
        do k=-TEST_BOX/2,0
            kp=OSMPL_PAD_FAC*k
            do h=-TEST_BOX/2,TEST_BOX/2
                if( h*h+k*k>(TEST_BOX/2)*(TEST_BOX/2+1) ) cycle
                hp=OSMPL_PAD_FAC*h
                fpl_ref%cmplx_plane(hp,kp)=cmplx(0.01*real(3*h-k+2),0.015*real(h+2*k-1))
                fpl_ref%ctfsq_plane(hp,kp)=0.25+0.01*real(mod(abs(h-2*k),11))
                fpl_ref%transfer_plane(hp,kp)=cmplx(0.5+0.02*real(mod(abs(2*h+k),9)), &
                    &0.01*real(mod(h-k,5)))
            end do
        end do
        call orientation%new(.true.)
        call orientation%set_euler([23.,47.,19.])
        call mean_rec%project_fplane(orientation,fpl_ref,mean_fpl_std,apply_ctf_amp=.true.)
        call project_fplane_mean(mean_rec,orientation,fpl_ref,mean_fpl,apply_ctf_amp=.true.)
        err=maxval(abs(mean_fpl%cmplx_plane-mean_fpl_std%cmplx_plane))
        if( err>3.e-5 ) THROW_HARD('latent mean projector differs from reconstructor projector')
        call project_fplanes_mean_basis(mean_rec,basis_recs,orientation,fpl_ref,mean_fpl,basis_fpls, &
            &apply_ctf_amp=.true.)
        err=maxval(abs(mean_fpl%cmplx_plane-mean_fpl_std%cmplx_plane))
        if( err>3.e-5 ) THROW_HARD('fused latent mean projector differs from reconstructor projector')
        do q=1,TEST_NBASIS
            call basis_recs(q)%project_fplane(orientation,fpl_ref,basis_fpls_std(q),apply_ctf_amp=.true.)
            err=maxval(abs(basis_fpls(q)%cmplx_plane-basis_fpls_std(q)%cmplx_plane))
            if( err>3.e-5 ) THROW_HARD('fused latent basis projector differs from reconstructor projector')
        end do
        call se%new('c1')
        call splat_std%new([TEST_BOX,TEST_BOX,TEST_BOX],1.)
        call splat_std%alloc_rho(params,spproj,expand=.true.)
        call splat_latent(1)%new([TEST_BOX,TEST_BOX,TEST_BOX],1.)
        call splat_latent(1)%alloc_rho(params,spproj,expand=.true.)
        fpl_ref%cmplx_plane = conjg(fpl_ref%transfer_plane) * fpl_ref%cmplx_plane
        fpl_ref%transfer_plane = cmplx(1.,0.)
        allocate(rho_cross(1,size(splat_latent(1)%cmat_exp,1),size(splat_latent(1)%cmat_exp,2), &
            &size(splat_latent(1)%cmat_exp,3)),source=0.)
        unit_scale=1.d0
        unit_second=1.d0
        call splat_std%insert_plane_oversamp(se,orientation,fpl_ref)
        call insert_plane_oversamp_coupled_scaled(splat_latent,rho_cross,se,orientation,fpl_ref, &
            &unit_scale,unit_second)
        err=maxval(abs(splat_std%cmat_exp-splat_latent(1)%cmat_exp))
        if( err>3.e-5 ) write(logfhandle,'(A,ES12.4)') 'latent splat numerator max_abs_error=',err
        if( err>3.e-5 ) THROW_HARD('latent splat numerator differs from reconstructor insertion')
        err=maxval(abs(splat_std%rho_exp-rho_cross(1,:,:,:)))
        if( err>3.e-5 ) THROW_HARD('latent splat density differs from reconstructor insertion')
        kbwin=kbinterpol(KBWINSZ,KBALPHA)
        loc=[1.17,-0.43,2.31]
        call latent_projection_weights(kbwin,loc,win,wx,wy,wz)
        do iz=1,LATENT_WDIM
            do iy=1,LATENT_WDIM
                do ix=1,LATENT_WDIM
                    x(ix,iy,iz)=cmplx(0.07*real(2*ix-iy+iz),0.03*real(ix+iy-2*iz))
                end do
            end do
        end do
        y=cmplx(0.37,-0.21)
        gathered=CMPLX_ZERO
        do iz=1,LATENT_WDIM
            do iy=1,LATENT_WDIM
                do ix=1,LATENT_WDIM
                    gathered=gathered+x(ix,iy,iz)*(wx(ix)*wy(iy)*wz(iz))
                    splat(ix,iy,iz)=y*(wx(ix)*wy(iy)*wz(iz))
                end do
            end do
        end do
        lhs=conjg(gathered)*y
        rhs=sum(conjg(x)*splat)
        if( abs(lhs-rhs)>2.e-6 ) THROW_HARD('Cartesian KB gather/splat adjoint check failed')
        call cleanup_test_plane(fpl_ref)
        call cleanup_test_plane(mean_fpl)
        call cleanup_test_plane(mean_fpl_std)
        do q=1,TEST_NBASIS
            call cleanup_test_plane(basis_fpls(q))
            call cleanup_test_plane(basis_fpls_std(q))
            call basis_recs(q)%dealloc_rho
            call basis_recs(q)%kill
        end do
        call mean_rec%dealloc_rho
        call mean_rec%kill
        call splat_std%dealloc_rho
        call splat_std%kill
        call splat_latent(1)%dealloc_rho
        call splat_latent(1)%kill
        call se%kill
        deallocate(rho_cross)
        call orientation%kill
        call spproj%kill
    contains
        subroutine cleanup_test_plane( fpl )
            type(fplane_type), intent(inout) :: fpl
            if( allocated(fpl%cmplx_plane) ) deallocate(fpl%cmplx_plane)
            if( allocated(fpl%ctfsq_plane) ) deallocate(fpl%ctfsq_plane)
            if( allocated(fpl%transfer_plane) ) deallocate(fpl%transfer_plane)
        end subroutine cleanup_test_plane
    end subroutine test_cartesian_projection_contract

    subroutine accumulate_plane_oversamp_coupled_stats( basis_rhs, rho_cross_exp, exp_lb, model_nyq, &
        &se, o, fpl, data_scales, density_scales )
        use simple_math, only: ceil_div, floor_div
        complex,            intent(inout) :: basis_rhs(:,:,:,:)
        real,               intent(inout) :: rho_cross_exp(:,:,:,:)
        integer,            intent(in)    :: exp_lb(3), model_nyq
        class(sym),         intent(inout) :: se
        class(ori),         intent(inout) :: o
        class(fplane_type), intent(in)    :: fpl
        real(dp),           intent(in)    :: data_scales(:), density_scales(:,:)
        type(kbinterpol) :: kbwin
        type(ori) :: o_sym
        complex   :: comp_base, cmplx_raw
        real      :: rotmats(se%get_nsym(),3,3), loc(3), hrow(3), ctfsq_raw
        real      :: wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM), ww
        real      :: data_scale_sp(size(basis_rhs,4)), &
            &density_scale_packed((size(basis_rhs,4)*(size(basis_rhs,4)+1))/2)
        real      :: r11, r12, r13, r21, r22, r23
        integer   :: win(2,3), h, k, l, nsym, isym, iwinsz, stride, fpllims_pd(3,2)
        integer   :: fpllims(3,2), hp, kp, pf, ix, iy, iz, hx, ky, mz, q, r, ncomp, ipair
        integer   :: nyq_disk, nyq_eff, h_sq, k_max_h, k_lo, k_hi, ih, ik, im
        integer   :: exp_ub(3), exp_shape(3), npairs
        real      :: pf2, eps_norm, inv_wdim
        ncomp = size(basis_rhs,4)
        if( ncomp <= 0 ) return
        npairs   = (ncomp * (ncomp + 1)) / 2
        exp_shape = shape(basis_rhs(:,:,:,1))
        exp_ub    = exp_lb + exp_shape - 1
        if( model_nyq < 1 ) THROW_HARD('invalid model Nyquist; accumulate_plane_oversamp_coupled_stats')
        if( size(data_scales) < ncomp .or. size(density_scales,1) < ncomp .or. size(density_scales,2) < ncomp )then
            THROW_HARD('scale array smaller than statistics component count; accumulate_plane_oversamp_coupled_stats')
        endif
        if( size(rho_cross_exp,1) < npairs .or. size(rho_cross_exp,2) < exp_shape(1) .or. &
            &size(rho_cross_exp,3) < exp_shape(2) .or. size(rho_cross_exp,4) < exp_shape(3) )then
            THROW_HARD('cross-density array shape mismatch; accumulate_plane_oversamp_coupled_stats')
        endif
        if( .not. allocated(fpl%transfer_plane) )then
            THROW_HARD('forward transfer plane does not exist; accumulate_plane_oversamp_coupled_stats')
        endif
        kbwin    = kbinterpol(KBWINSZ, KBALPHA)
        iwinsz   = ceiling(KBWINSZ - 0.5)
        stride   = LATENT_SAFE_STRIDE
        nsym     = se%get_nsym()
        rotmats(1,:,:) = o%get_mat()
        if( nsym > 1 )then
            do isym = 2, nsym
                call se%apply(o, isym, o_sym)
                rotmats(isym,:,:) = o_sym%get_mat()
            end do
        endif
        fpllims_pd   = fpl%frlims
        pf           = OSMPL_PAD_FAC
        pf2          = real(pf*pf)
        fpllims      = fpllims_pd
        fpllims(1,1) = ceil_div (fpllims_pd(1,1), pf)
        fpllims(1,2) = floor_div(fpllims_pd(1,2), pf)
        fpllims(2,1) = ceil_div (fpllims_pd(2,1), pf)
        fpllims(2,2) = floor_div(fpllims_pd(2,2), pf)
        nyq_eff = model_nyq
        if( fpl%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpl%nyq / pf))
        nyq_disk = nyq_eff * (nyq_eff + 1)
        eps_norm = epsilon(1.0)
        inv_wdim = 1.0 / real(LATENT_WDIM)
        do q = 1, ncomp
            data_scale_sp(q) = real(data_scales(q))
        end do
        do r = 1, ncomp
            do q = 1, r
                density_scale_packed(pair_index(q,r)) = real(density_scales(q,r))
            end do
        end do
        !$omp parallel default(shared) private(h,k,l,h_sq,k_max_h,k_lo,k_hi,cmplx_raw,&
        !$omp& ctfsq_raw,comp_base,wx,wy,wz,ww,win,loc,hrow,hp,kp,r11,r12,r13,r21,r22,r23,&
        !$omp& isym,ix,iy,iz,hx,ky,mz,ih,ik,im,q,r,ipair) proc_bind(close)
        do isym = 1, nsym
            r11 = rotmats(isym,1,1); r12 = rotmats(isym,1,2); r13 = rotmats(isym,1,3)
            r21 = rotmats(isym,2,1); r22 = rotmats(isym,2,2); r23 = rotmats(isym,2,3)
            do l = 0, stride-1
                !$omp do schedule(static,1)
                do h = fpllims(1,1)+l, fpllims(1,2), stride
                    h_sq = h*h
                    if( h_sq > nyq_disk ) cycle
                    k_max_h = int(sqrt(real(nyq_disk - h_sq)))
                    k_lo    = max(fpllims(2,1), -k_max_h)
                    k_hi    = min(fpllims(2,2),  k_max_h)
                    hp      = h * pf
                    hrow(1) = real(h) * r11
                    hrow(2) = real(h) * r12
                    hrow(3) = real(h) * r13
                    do k = k_lo, k_hi
                        kp = k * pf
                        if( kp <= 0 )then
                            cmplx_raw = conjg(fpl%transfer_plane(hp,kp)) * fpl%cmplx_plane(hp,kp)
                            ctfsq_raw = fpl%ctfsq_plane(hp,kp)
                        else
                            cmplx_raw = conjg(conjg(fpl%transfer_plane(-hp,-kp)) * fpl%cmplx_plane(-hp,-kp))
                            ctfsq_raw = fpl%ctfsq_plane(-hp,-kp)
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
                        call kb_apod_vecs_3d_fast(loc, wx, wy, wz)
                        do iz = 1, LATENT_WDIM
                            mz = win(1,3) + iz - 1
                            im = mz - exp_lb(3) + 1
                            do iy = 1, LATENT_WDIM
                                ky = win(1,2) + iy - 1
                                ik = ky - exp_lb(2) + 1
                                do ix = 1, LATENT_WDIM
                                    hx = win(1,1) + ix - 1
                                    ih = hx - exp_lb(1) + 1
                                    ww = wx(ix) * (wy(iy) * wz(iz))
                                    do q = 1, ncomp
                                        basis_rhs(ih,ik,im,q) = basis_rhs(ih,ik,im,q) + &
                                            &(data_scale_sp(q) * comp_base) * ww
                                    end do
                                    !$omp simd
                                    do ipair = 1, npairs
                                        rho_cross_exp(ipair,ih,ik,im) = rho_cross_exp(ipair,ih,ik,im) + &
                                            &(density_scale_packed(ipair) * ctfsq_raw) * ww
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
                !$omp end do
            end do
        end do
        !$omp end parallel
        call o_sym%kill

    contains

        integer pure function pair_index( q, r ) result( ipair )
            integer, intent(in) :: q, r
            ipair = (r * (r - 1)) / 2 + q
        end function pair_index

        subroutine kb_apod_vecs_3d_fast( loc, wx, wy, wz )
            real, intent(in)  :: loc(3)
            real, intent(out) :: wx(:), wy(:), wz(:)
            integer :: i, win_lo(3)
            real    :: base(3), ww3(3), sx, sy, sz
            win_lo = nint(loc) - iwinsz
            base   = real(win_lo) - loc
            do i = 1, LATENT_WDIM
                ww3   = kbwin%apod_fast(base + real(i-1))
                wx(i) = ww3(1)
                wy(i) = ww3(2)
                wz(i) = ww3(3)
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
        end subroutine kb_apod_vecs_3d_fast

    end subroutine accumulate_plane_oversamp_coupled_stats

    subroutine accumulate_planes_oversamp_coupled_stats_batch( basis_rhs, rho_cross_exp, exp_lb, model_nyq, &
        &se, orientations, fpls, data_scales, density_scales, valid, nrecords )
        use simple_math, only: ceil_div, floor_div
        use simple_kbinterpol, only: apod_kb15_a2
        complex,            intent(inout) :: basis_rhs(:,:,:,:)
        real,               intent(inout) :: rho_cross_exp(:,:,:,:)
        integer,            intent(in)    :: exp_lb(3), model_nyq, nrecords
        class(sym),         intent(inout) :: se
        type(ori),          intent(inout) :: orientations(:)
        type(fplane_type),  intent(in)    :: fpls(:)
        real(dp),           intent(in)    :: data_scales(:,:), density_scales(:,:,:)
        logical,            intent(in)    :: valid(:)
        type(ori) :: o_sym
        complex :: comp_base,cmplx_raw
        real, allocatable :: rotmats(:,:,:,:),data_scale_sp(:,:),density_scale_packed(:,:)
        integer, allocatable :: fpllims(:,:,:),nyq_disks(:)
        real :: loc(3),hrow(3),ctfsq_raw,wx(LATENT_WDIM),wy(LATENT_WDIM),wz(LATENT_WDIM),ww
        real :: r11,r12,r13,r21,r22,r23,pf2,eps_norm,inv_wdim
        integer :: win(2,3),fpllims_pd(3,2),exp_ub(3),exp_shape(3)
        integer :: h,k,l,isym,nsym,iwinsz,stride,hp,kp,pf,ix,iy,iz,hx,ky,mz
        integer :: q,r,i,ncomp,npairs,ipair,h_sq,k_max_h,k_lo,k_hi,ih,ik,im,nyq_eff
        ncomp=size(basis_rhs,4)
        if( ncomp<=0 .or. nrecords<=0 ) return
        npairs=(ncomp*(ncomp+1))/2
        exp_shape=shape(basis_rhs(:,:,:,1))
        exp_ub=exp_lb+exp_shape-1
        if( model_nyq<1 ) THROW_HARD('invalid model Nyquist; accumulate_planes_oversamp_coupled_stats_batch')
        if( size(orientations)<nrecords .or. size(fpls)<nrecords .or. size(valid)<nrecords )then
            THROW_HARD('record array smaller than batch; accumulate_planes_oversamp_coupled_stats_batch')
        endif
        if( size(data_scales,1)<ncomp .or. size(data_scales,2)<nrecords .or. &
            &size(density_scales,1)<ncomp .or. size(density_scales,2)<ncomp .or. &
            &size(density_scales,3)<nrecords )then
            THROW_HARD('scale array smaller than batch; accumulate_planes_oversamp_coupled_stats_batch')
        endif
        if( size(rho_cross_exp,1)<npairs .or. size(rho_cross_exp,2)<exp_shape(1) .or. &
            &size(rho_cross_exp,3)<exp_shape(2) .or. size(rho_cross_exp,4)<exp_shape(3) )then
            THROW_HARD('cross-density shape mismatch; accumulate_planes_oversamp_coupled_stats_batch')
        endif
        nsym=se%get_nsym()
        iwinsz=ceiling(KBWINSZ-0.5)
        stride=LATENT_SAFE_STRIDE
        pf=OSMPL_PAD_FAC
        pf2=real(pf*pf)
        eps_norm=epsilon(1.0)
        inv_wdim=1.0/real(LATENT_WDIM)
        allocate(rotmats(3,3,nsym,nrecords),data_scale_sp(ncomp,nrecords), &
            &density_scale_packed(npairs,nrecords),source=0.)
        allocate(fpllims(3,2,nrecords),nyq_disks(nrecords),source=0)
        do i=1,nrecords
            if( .not.valid(i) ) cycle
            if( .not.allocated(fpls(i)%transfer_plane) )then
                THROW_HARD('forward transfer plane does not exist; accumulate_planes_oversamp_coupled_stats_batch')
            endif
            rotmats(:,:,1,i)=orientations(i)%get_mat()
            do isym=2,nsym
                call se%apply(orientations(i),isym,o_sym)
                rotmats(:,:,isym,i)=o_sym%get_mat()
            end do
            fpllims_pd=fpls(i)%frlims
            fpllims(:,:,i)=fpllims_pd
            fpllims(1,1,i)=ceil_div (fpllims_pd(1,1),pf)
            fpllims(1,2,i)=floor_div(fpllims_pd(1,2),pf)
            fpllims(2,1,i)=ceil_div (fpllims_pd(2,1),pf)
            fpllims(2,2,i)=floor_div(fpllims_pd(2,2),pf)
            nyq_eff=model_nyq
            if( fpls(i)%nyq>0 ) nyq_eff=min(nyq_eff,max(1,fpls(i)%nyq/pf))
            nyq_disks(i)=nyq_eff*(nyq_eff+1)
            do q=1,ncomp
                data_scale_sp(q,i)=real(data_scales(q,i))
            end do
            do r=1,ncomp
                do q=1,r
                    density_scale_packed(pair_index(q,r),i)=real(density_scales(q,r,i))
                end do
            end do
        end do
        call o_sym%kill
        !$omp parallel default(shared) private(i,h,k,l,h_sq,k_max_h,k_lo,k_hi,cmplx_raw,ctfsq_raw,&
        !$omp& comp_base,wx,wy,wz,ww,win,loc,hrow,hp,kp,r11,r12,r13,r21,r22,r23,isym,&
        !$omp& ix,iy,iz,hx,ky,mz,ih,ik,im,q,ipair) proc_bind(close)
        do i=1,nrecords
            if( .not.valid(i) ) cycle
            do isym=1,nsym
                r11=rotmats(1,1,isym,i); r12=rotmats(1,2,isym,i); r13=rotmats(1,3,isym,i)
                r21=rotmats(2,1,isym,i); r22=rotmats(2,2,isym,i); r23=rotmats(2,3,isym,i)
                do l=0,stride-1
                    !$omp do schedule(static,1)
                    do h=fpllims(1,1,i)+l,fpllims(1,2,i),stride
                        h_sq=h*h
                        if( h_sq>nyq_disks(i) ) cycle
                        k_max_h=int(sqrt(real(nyq_disks(i)-h_sq)))
                        k_lo=max(fpllims(2,1,i),-k_max_h)
                        k_hi=min(fpllims(2,2,i), k_max_h)
                        hp=h*pf
                        hrow=real(h)*[r11,r12,r13]
                        loc=hrow+real(k_lo-1)*[r21,r22,r23]
                        do k=k_lo,k_hi
                            loc=loc+[r21,r22,r23]
                            kp=k*pf
                            if( kp<=0 )then
                                cmplx_raw=conjg(fpls(i)%transfer_plane(hp,kp))*fpls(i)%cmplx_plane(hp,kp)
                                ctfsq_raw=fpls(i)%ctfsq_plane(hp,kp)
                            else
                                cmplx_raw=conjg(conjg(fpls(i)%transfer_plane(-hp,-kp))*fpls(i)%cmplx_plane(-hp,-kp))
                                ctfsq_raw=fpls(i)%ctfsq_plane(-hp,-kp)
                            endif
                            if( abs(real(cmplx_raw))+abs(aimag(cmplx_raw))<=TINY .and. ctfsq_raw<=TINY ) cycle
                            win(1,:)=nint(loc)-iwinsz
                            win(2,:)=nint(loc)+iwinsz
                            if( win(2,1)<NONREDUNDANT_HMIN ) cycle
                            if( any(win(1,:)<exp_lb) .or. any(win(2,:)>exp_ub) ) cycle
                            comp_base=pf2*cmplx_raw
                            call kb_apod_vecs_3d_fast(loc,wx,wy,wz)
                            do iz=1,LATENT_WDIM
                                mz=win(1,3)+iz-1; im=mz-exp_lb(3)+1
                                do iy=1,LATENT_WDIM
                                    ky=win(1,2)+iy-1; ik=ky-exp_lb(2)+1
                                    do ix=1,LATENT_WDIM
                                        hx=win(1,1)+ix-1; ih=hx-exp_lb(1)+1
                                        ww=wx(ix)*(wy(iy)*wz(iz))
                                        do q=1,ncomp
                                            basis_rhs(ih,ik,im,q)=basis_rhs(ih,ik,im,q)+ &
                                                &(data_scale_sp(q,i)*comp_base)*ww
                                        end do
                                        !$omp simd
                                        do ipair=1,npairs
                                            rho_cross_exp(ipair,ih,ik,im)=rho_cross_exp(ipair,ih,ik,im)+ &
                                                &(density_scale_packed(ipair,i)*ctfsq_raw)*ww
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
        deallocate(rotmats,data_scale_sp,density_scale_packed,fpllims,nyq_disks)

    contains

        integer pure function pair_index(q,r) result(ipair)
            integer,intent(in)::q,r
            ipair=(r*(r-1))/2+q
        end function pair_index

        subroutine kb_apod_vecs_3d_fast(loc,wx,wy,wz)
            real,intent(in)::loc(3)
            real,intent(out)::wx(:),wy(:),wz(:)
            integer::j,win_lo(3)
            real::base(3),sx,sy,sz
            win_lo=nint(loc)-iwinsz
            base=real(win_lo)-loc
            do j=1,LATENT_WDIM
                wx(j)=apod_kb15_a2(base(1)+real(j-1))
                wy(j)=apod_kb15_a2(base(2)+real(j-1))
                wz(j)=apod_kb15_a2(base(3)+real(j-1))
            end do
            sx=sum(wx); sy=sum(wy); sz=sum(wz)
            if( abs(sx)>eps_norm )then; wx=wx/sx; else; wx=inv_wdim; endif
            if( abs(sy)>eps_norm )then; wy=wy/sy; else; wy=inv_wdim; endif
            if( abs(sz)>eps_norm )then; wz=wz/sz; else; wz=inv_wdim; endif
        end subroutine kb_apod_vecs_3d_fast

    end subroutine accumulate_planes_oversamp_coupled_stats_batch

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
                mean_val = weighted_expanded_cmat(mean_rec, win, wx, wy, wz)
                if( l_conjg ) mean_val = conjg(mean_val)
                mean_fpl%cmplx_plane(hp,kp) = transfer * mean_val
                do q = 1, ncomp
                    basis_val = weighted_expanded_cmat(basis_recs(q), win, wx, wy, wz)
                    if( l_conjg ) basis_val = conjg(basis_val)
                    basis_fpls(q)%cmplx_plane(hp,kp) = transfer * basis_val
                end do
            end do
        end do

    end subroutine project_fplanes_mean_basis

    subroutine ensure_latent_projection_plane( fpl_in, fpl_out )
        type(fplane_type), intent(in)    :: fpl_in
        type(fplane_type), intent(inout) :: fpl_out
        logical :: l_realloc
        l_realloc = .not. allocated(fpl_out%cmplx_plane)
        if( .not. l_realloc )then
            l_realloc = any(lbound(fpl_out%cmplx_plane) /= lbound(fpl_in%cmplx_plane)) .or. &
                &any(ubound(fpl_out%cmplx_plane) /= ubound(fpl_in%cmplx_plane))
        endif
        if( l_realloc )then
            if( allocated(fpl_out%cmplx_plane) ) deallocate(fpl_out%cmplx_plane)
            allocate(fpl_out%cmplx_plane(lbound(fpl_in%cmplx_plane,1):ubound(fpl_in%cmplx_plane,1), &
                &lbound(fpl_in%cmplx_plane,2):ubound(fpl_in%cmplx_plane,2)))
        endif
        if( allocated(fpl_out%ctfsq_plane) ) deallocate(fpl_out%ctfsq_plane)
        if( allocated(fpl_out%transfer_plane) ) deallocate(fpl_out%transfer_plane)
        fpl_out%frlims  = fpl_in%frlims
        fpl_out%shconst = fpl_in%shconst
        fpl_out%nyq     = fpl_in%nyq
        fpl_out%cmplx_plane = CMPLX_ZERO
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

    pure function weighted_expanded_cmat( rec, win, wx, wy, wz ) result( val )
        type(reconstructor), intent(in) :: rec
        integer,             intent(in) :: win(2,3)
        real,                intent(in) :: wx(:), wy(:), wz(:)
        complex :: val
        integer :: ix, iy, iz, hx, ky, mz, lb(3), ub(3)
        real    :: wyz
        val = CMPLX_ZERO
        lb = lbound(rec%cmat_exp)
        ub = ubound(rec%cmat_exp)
        if( any(win(1,:) < lb) .or. any(win(2,:) > ub) ) return
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

end module simple_flex_reconstructor_latent_ops
