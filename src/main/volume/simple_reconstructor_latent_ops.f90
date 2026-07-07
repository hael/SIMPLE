!@descr: Latent-volume projection/backprojection helpers for flex workflows
module simple_reconstructor_latent_ops
use simple_core_module_api
use simple_reconstructor, only: reconstructor
implicit none

public :: insert_plane_oversamp_multi_scaled, insert_plane_oversamp_coupled_scaled, project_fplanes_mean_basis
private
#include "simple_local_flags.inc"

integer, parameter :: LATENT_WDIM = 2 * ceiling(KBWINSZ - 0.5) + 1

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
        stride   = LATENT_WDIM
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
        real      :: data_scale_sp(size(recs)), density_scale_sp(size(recs),size(recs))
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
        exp_lb    = lbound(recs(1)%cmat_exp)
        exp_ub    = ubound(recs(1)%cmat_exp)
        exp_shape = shape(recs(1)%cmat_exp)
        if( size(rho_cross_exp,1) < npairs .or. size(rho_cross_exp,2) < exp_shape(1) .or. &
            &size(rho_cross_exp,3) < exp_shape(2) .or. size(rho_cross_exp,4) < exp_shape(3) )then
            THROW_HARD('cross-density array shape mismatch; insert_plane_oversamp_coupled_scaled')
        endif
        kbwin    = kbinterpol(KBWINSZ, KBALPHA)
        iwinsz   = ceiling(KBWINSZ - 0.5)
        stride   = LATENT_WDIM
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
            do r = 1, ncomp
                density_scale_sp(q,r) = real(density_scales(q,r))
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
                                        do r = q, ncomp
                                            ipair = pair_index(q, r)
                                            rho_cross_exp(ipair,ih,ik,im) = rho_cross_exp(ipair,ih,ik,im) + &
                                                &(density_scale_sp(q,r) * ctfsq_raw) * ww
                                        end do
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

    subroutine project_fplanes_mean_basis( mean_rec, basis_recs, o, fpl_ref, mean_fpl, basis_fpls, apply_ctf_amp )
        use simple_math, only: ceil_div, floor_div
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), intent(inout) :: basis_recs(:)
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
        logical :: l_apply_ctf_amp
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
        call ensure_projection_plane(fpl_ref, mean_fpl)
        do q = 1, ncomp
            call ensure_projection_plane(fpl_ref, basis_fpls(q))
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
                call projection_weights(loc, win, wx, wy, wz)
                transfer = cmplx(1., 0.)
                if( l_apply_ctf_amp )then
                    if( allocated(fpl_ref%transfer_plane) )then
                        transfer = fpl_ref%transfer_plane(hp,kp)
                    else
                        ctfamp   = sqrt(max(0., fpl_ref%ctfsq_plane(hp,kp)))
                        transfer = cmplx(ctfamp, 0.)
                    endif
                endif
                mean_val = weighted_cmat(mean_rec, win, wx, wy, wz)
                mean_fpl%cmplx_plane(hp,kp) = transfer * mean_val
                do q = 1, ncomp
                    basis_val = weighted_cmat(basis_recs(q), win, wx, wy, wz)
                    basis_fpls(q)%cmplx_plane(hp,kp) = transfer * basis_val
                end do
            end do
        end do

    contains

        subroutine ensure_projection_plane( fpl_in, fpl_out )
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
        end subroutine ensure_projection_plane

        subroutine projection_weights( loc, win, wx, wy, wz )
            real,    intent(in)  :: loc(3)
            integer, intent(out) :: win(2,3)
            real,    intent(out) :: wx(:), wy(:), wz(:)
            integer :: i, iwinsz, win_lo(3)
            real    :: base(3), ww3(3), sx, sy, sz, inv_wdim, eps_norm
            iwinsz   = ceiling(KBWINSZ - 0.5)
            win(1,:) = nint(loc)
            win(2,:) = win(1,:) + iwinsz
            win(1,:) = win(1,:) - iwinsz
            win_lo   = win(1,:)
            base     = real(win_lo) - loc
            do i = 1, LATENT_WDIM
                ww3   = kbwin%apod_fast(base + real(i-1))
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
        end subroutine projection_weights

        function weighted_cmat( rec, win, wx, wy, wz ) result( val )
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
        end function weighted_cmat

    end subroutine project_fplanes_mean_basis

end module simple_reconstructor_latent_ops
