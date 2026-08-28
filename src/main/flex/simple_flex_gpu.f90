!@descr: CUDA-C GPU path for the flex_pca insertion family (P1 of the polar/GPU plan)
!
! Thin host side over src/main/flex/cuda/simple_flex_gpu_kernels.cu. The GPU accumulates
! per-component cmat/rho on-device across batches; at stage end the totals are fetched and
! ADDED into the reconstructor pair, so everything downstream (normalization, e/o handling,
! FFTs) is untouched. Numerics match insert_planes_oversamp_multi_scaled_batch except for
! fp32 atomic summation order — gate statistically, never byte-wise, across this boundary.
!
! Compiled in only under USE_FLEX_CUDA; without it every entry point is a hard error except
! flex_gpu_available(), which is then .false. — call sites guard on it.
module simple_flex_gpu
use simple_core_module_api
use simple_image,         only: image
use simple_reconstructor, only: reconstructor
use simple_kbinterpol,    only: kbinterpol
use iso_c_binding
implicit none

public :: flex_gpu_available
public :: flex_gpu_insert_begin_f, flex_gpu_insert_batch_f, flex_gpu_insert_batch_res_f
public :: flex_gpu_insert_end_f
public :: flex_gpu_coupled_begin_f, flex_gpu_coupled_batch_f, flex_gpu_coupled_batch_raw_f
public :: flex_gpu_coupled_end_f, flex_gpu_coupled_bank_f, flex_gpu_coupled_batch_banked_f
public :: flex_gpu_coupled_bank_free_f
public :: flex_gpu_cols_begin_f, flex_gpu_cols_batch_f, flex_gpu_cols_end_f
public :: flex_gpu_cols_lookup_f, flex_gpu_cols_fused_batch_f, flex_gpu_cols_test_upload_f
public :: flex_gpu_snr_begin_f, flex_gpu_snr_batch_f, flex_gpu_snr_batch_res_f
public :: flex_gpu_snr_end_f
public :: flex_gpu_solve_begin_f, flex_gpu_solve_cache_f, flex_gpu_solve_chunk_f, flex_gpu_solve_end_f
public :: flex_gpu_psample_begin_f, flex_gpu_psample_batch_f, flex_gpu_psample_batch_res_f
public :: flex_gpu_psample_free_f
public :: flex_gpu_psample_diag_begin_f, flex_gpu_psample_diag_batch_f, flex_gpu_psample_diag_fetch_f
public :: flex_gpu_estep_vols_f, flex_gpu_estep_batch_f, flex_gpu_estep_resid_f, flex_gpu_estep_free_f
public :: flex_gpu_coupled_batch_banked_res_f
public :: flex_gpu_prep_begin_f, flex_gpu_prep_batch_f, flex_gpu_prep_free_f, flex_gpu_estep_batch_res_f
public :: flex_gpu_prep_check_f, flex_gpu_prep_fetch_sep_f, flex_gpu_prep_fetch_f
public :: flex_gpu_prep_ready
public :: flex_gpu_poles_begin_f, flex_gpu_poles_bank_f, flex_gpu_poles_batch_f, flex_gpu_poles_free_f
public :: test_flex_gpu_insert, test_flex_gpu_coupled, test_flex_gpu_coupled_banked
public :: test_flex_gpu_psample, test_flex_gpu_estep
private
#include "simple_local_flags.inc"

#ifdef USE_FLEX_CUDA
interface
    integer(c_int) function flex_gpu_device_count() bind(c, name="flex_gpu_device_count")
        import :: c_int
    end function
    integer(c_int) function c_insert_begin(ncomp, nvox) bind(c, name="flex_gpu_insert_begin")
        import :: c_int, c_long
        integer(c_int),  value :: ncomp
        integer(c_long), value :: nvox
    end function
    integer(c_int) function c_insert_batch(pl_c, pl_ct, rotm, dscale, rscale, act, nact, &
        &nyqdisk, nrec, ncomp, nsym, hpd_lo, hpd_hi, kpd_lo, hlo, hhi, klo, khi, pf, &
        &lb, ub, tiny_val) bind(c, name="flex_gpu_insert_batch")
        import :: c_int, c_ptr, c_float
        type(c_ptr),    value :: pl_c, pl_ct, rotm, dscale, rscale, act, nact, nyqdisk
        integer(c_int), value :: nrec, ncomp, nsym, hpd_lo, hpd_hi, kpd_lo
        integer(c_int), value :: hlo, hhi, klo, khi, pf
        integer(c_int)        :: lb(3), ub(3)
        real(c_float),  value :: tiny_val
    end function
    integer(c_int) function c_insert_batch_res(rotm, dscale, rscale, act, nact, nyqdisk, &
        &nrec, ncomp, nsym, lb, ub, tiny_val) bind(c, name="flex_gpu_insert_batch_resident")
        import :: c_int, c_ptr, c_float
        type(c_ptr),    value :: rotm, dscale, rscale, act, nact, nyqdisk
        integer(c_int), value :: nrec, ncomp, nsym
        integer(c_int)        :: lb(3), ub(3)
        real(c_float),  value :: tiny_val
    end function
    integer(c_int) function c_insert_fetch(q, cmat_host, rho_host) bind(c, name="flex_gpu_insert_fetch")
        import :: c_int, c_ptr
        integer(c_int), value :: q
        type(c_ptr),    value :: cmat_host, rho_host
    end function
    integer(c_int) function c_insert_free() bind(c, name="flex_gpu_insert_free")
        import :: c_int
    end function
    integer(c_int) function c_coupled_begin(ncomp, nrho, nvox) bind(c, name="flex_gpu_coupled_begin")
        import :: c_int, c_long
        integer(c_int),  value :: ncomp, nrho
        integer(c_long), value :: nvox
    end function
    integer(c_int) function c_coupled_batch(pl_c, pl_ct, rotm, dscale, rscale, vmask, nyqdisk, &
        &nrec, ncomp, nrho, nsym, hpd_lo, hpd_hi, kpd_lo, hlo, hhi, klo, khi, pf, &
        &lb, ub, tiny_val) bind(c, name="flex_gpu_coupled_batch")
        import :: c_int, c_ptr, c_float
        type(c_ptr),    value :: pl_c, pl_ct, rotm, dscale, rscale, vmask, nyqdisk
        integer(c_int), value :: nrec, ncomp, nrho, nsym, hpd_lo, hpd_hi, kpd_lo
        integer(c_int), value :: hlo, hhi, klo, khi, pf
        integer(c_int)        :: lb(3), ub(3)
        real(c_float),  value :: tiny_val
    end function
    integer(c_int) function c_coupled_bank(rmat, ndir) bind(c, name="flex_gpu_coupled_bank")
        import :: c_int, c_ptr
        type(c_ptr),    value :: rmat
        integer(c_int), value :: ndir
    end function
    integer(c_int) function c_coupled_batch_banked(pl_c, pl_ct, cav, sav, dscale, rscale, &
        &vmask, dirid, nrec, ncomp, nrho, hpd_lo, hpd_hi, kpd_lo, hlo, hhi, klo, khi, pf, &
        &nyq2, lb, ub) bind(c, name="flex_gpu_coupled_batch_banked")
        import :: c_int, c_ptr
        type(c_ptr),    value :: pl_c, pl_ct, cav, sav, dscale, rscale, vmask
        integer(c_int), value :: nrec, ncomp, nrho, hpd_lo, hpd_hi, kpd_lo
        integer(c_int), value :: hlo, hhi, klo, khi, pf, nyq2
        integer(c_int)        :: dirid(*)
        integer(c_int)        :: lb(3), ub(3)
    end function
    integer(c_int) function c_coupled_bank_free() bind(c, name="flex_gpu_coupled_bank_free")
        import :: c_int
    end function
    integer(c_int) function c_coupled_fetch_cmat(q, cmat_host) bind(c, name="flex_gpu_coupled_fetch_cmat")
        import :: c_int, c_ptr
        integer(c_int), value :: q
        type(c_ptr),    value :: cmat_host
    end function
    integer(c_int) function c_coupled_fetch_rho(rho_host) bind(c, name="flex_gpu_coupled_fetch_rho")
        import :: c_int, c_ptr
        type(c_ptr),    value :: rho_host
    end function
    integer(c_int) function c_coupled_free() bind(c, name="flex_gpu_coupled_free")
        import :: c_int
    end function
    integer(c_int) function c_cols_begin(ncol, nvox) bind(c, name="flex_gpu_cols_begin")
        import :: c_int, c_long
        integer(c_int),  value :: ncol
        integer(c_long), value :: nvox
    end function
    integer(c_int) function c_cols_batch(cwin, cwx, cwy, cwz, cpl, cct, ncache, vmask, gcol, &
        &hcolv, qhkl, nb, maxc, ncol, lb, ub, deadskip, debias, pf2) bind(c, name="flex_gpu_cols_batch")
        import :: c_int, c_ptr, c_float
        type(c_ptr),    value :: cwin, cwx, cwy, cwz, cpl, cct, ncache, vmask, gcol, hcolv, qhkl
        integer(c_int), value :: nb, maxc, ncol, deadskip, debias
        integer(c_int)        :: lb(3), ub(3)
        real(c_float),  value :: pf2
    end function
    integer(c_int) function c_cols_fetch(which, B_host, H_host) bind(c, name="flex_gpu_cols_fetch")
        import :: c_int, c_ptr
        integer(c_int), value :: which
        type(c_ptr),    value :: B_host, H_host
    end function
    integer(c_int) function c_cols_free() bind(c, name="flex_gpu_cols_free")
        import :: c_int
    end function
    integer(c_int) function c_snr_begin(nvox, nslab) bind(c, name="flex_gpu_snr_begin")
        import :: c_int, c_long
        integer(c_long), value :: nvox
        integer(c_int),  value :: nslab
    end function
    integer(c_int) function c_snr_batch(pl_c, pl_ct, rotm, nyqdisk, vmask, nrec, nsym, &
        &hpd_lo, hpd_hi, kpd_lo, hlo, hhi, klo, khi, pf, lb, ub, tiny_val, inv_pf2) &
        &bind(c, name="flex_gpu_snr_batch")
        import :: c_int, c_ptr, c_float
        type(c_ptr),    value :: pl_c, pl_ct, rotm, nyqdisk, vmask
        integer(c_int), value :: nrec, nsym, hpd_lo, hpd_hi, kpd_lo, hlo, hhi, klo, khi, pf
        integer(c_int)        :: lb(3), ub(3)
        real(c_float),  value :: tiny_val, inv_pf2
    end function
    integer(c_int) function c_snr_batch_res(rotm, nyqdisk, vmask, nrec, nsym, lb, ub, &
        &tiny_val, inv_pf2) bind(c, name="flex_gpu_snr_batch_resident")
        import :: c_int, c_ptr, c_float
        type(c_ptr),    value :: rotm, nyqdisk, vmask
        integer(c_int), value :: nrec, nsym
        integer(c_int)        :: lb(3), ub(3)
        real(c_float),  value :: tiny_val, inv_pf2
    end function
    integer(c_int) function c_snr_end(var_host, dens_host) bind(c, name="flex_gpu_snr_end")
        import :: c_int, c_ptr
        type(c_ptr), value :: var_host, dens_host
    end function
    integer(c_int) function c_solve_begin(np, nsamp2, nk, dt, npk, maxdir, maxm, unitc, alo, ahi) &
        &bind(c, name="flex_gpu_solve_begin")
        import :: c_int, c_long_long, c_float
        integer(c_long_long), value :: np, maxm
        integer(c_int),       value :: nsamp2, nk, dt, npk, maxdir, unitc
        real(c_float),        value :: alo, ahi
    end function
    integer(c_int) function c_solve_cache(xws, wrs) bind(c, name="flex_gpu_solve_cache")
        import :: c_int, c_ptr
        type(c_ptr), value :: xws, wrs
    end function
    integer(c_int) function c_solve_chunk(Us, Cpk, Cm0, c00, dlist, dcnt, ndc, mtot, wrsum_host) &
        &bind(c, name="flex_gpu_solve_chunk")
        import :: c_int, c_long_long, c_ptr
        type(c_ptr),          value :: Us, Cpk, Cm0, c00, dlist, dcnt, wrsum_host
        integer(c_int),       value :: ndc
        integer(c_long_long), value :: mtot
    end function
    integer(c_int) function c_solve_end(Sbb_host, Mspk_host) bind(c, name="flex_gpu_solve_end")
        import :: c_int, c_ptr
        type(c_ptr), value :: Sbb_host, Mspk_host
    end function
    integer(c_int) function c_psample_begin(rad, cs, sn, sqwq, ring_of, nang, nsamp, nk, &
        &nrad, ncs, nsn, nwq, nsampn) bind(c, name="flex_gpu_psample_begin")
        import :: c_int, c_ptr
        type(c_ptr),    value :: rad, cs, sn, sqwq, nrad, ncs, nsn, nwq
        integer(c_int)        :: ring_of(*), nang(*)
        integer(c_int), value :: nsamp, nk, nsampn
    end function
    integer(c_int) function c_psample_batch(ply, plt, cav, sav, vmask, nrec, want_halves, &
        &hpd_lo, hpd_hi, kpd_lo, xw, xw1, xw2, wr, pw, taz) bind(c, name="flex_gpu_psample_batch")
        import :: c_int, c_ptr
        type(c_ptr),    value :: ply, plt, cav, sav, vmask, xw, xw1, xw2, wr, pw, taz
        integer(c_int), value :: nrec, want_halves, hpd_lo, hpd_hi, kpd_lo
    end function
    integer(c_int) function c_psample_batch_res(cav, sav, vmask, nrec, want_halves, &
        &xw, xw1, xw2, wr, pw, taz) bind(c, name="flex_gpu_psample_batch_resident")
        import :: c_int, c_ptr
        type(c_ptr),    value :: cav, sav, vmask, xw, xw1, xw2, wr, pw, taz
        integer(c_int), value :: nrec, want_halves
    end function
    integer(c_int) function c_psample_diag_begin(nyq) bind(c, name="flex_gpu_psample_diag_begin")
        import :: c_int
        integer(c_int), value :: nyq
    end function
    integer(c_int) function c_psample_diag_batch(nrec, sh_lo) bind(c, name="flex_gpu_psample_diag_batch")
        import :: c_int
        integer(c_int), value :: nrec, sh_lo
    end function
    integer(c_int) function c_psample_diag_fetch(acc) bind(c, name="flex_gpu_psample_diag_fetch")
        import :: c_int, c_ptr
        type(c_ptr), value :: acc
    end function
    integer(c_int) function c_psample_free() bind(c, name="flex_gpu_psample_free")
        import :: c_int
    end function
    integer(c_int) function c_estep_vols(vol, q1, ncomp1, lb, ub) bind(c, name="flex_gpu_estep_vols")
        import :: c_int, c_ptr
        type(c_ptr),    value :: vol
        integer(c_int), value :: q1, ncomp1
        integer(c_int)        :: lb(3), ub(3)
    end function
    integer(c_int) function c_estep_batch(pl_c, pl_ct, rotm, vmask, nrec, hs_lo, hs_hi, &
        &ks_lo, khi, nyq2, stats_host) bind(c, name="flex_gpu_estep_batch")
        import :: c_int, c_ptr
        type(c_ptr),    value :: pl_c, pl_ct, rotm, vmask, stats_host
        integer(c_int), value :: nrec, hs_lo, hs_hi, ks_lo, khi, nyq2
    end function
    integer(c_int) function c_estep_resid(avec, vmask, nrec, hs_lo, hs_hi, ks_lo, khi, &
        &resid_host) bind(c, name="flex_gpu_estep_resid")
        import :: c_int, c_ptr
        type(c_ptr),    value :: avec, vmask, resid_host
        integer(c_int), value :: nrec, hs_lo, hs_hi, ks_lo, khi
    end function
    integer(c_int) function c_estep_free() bind(c, name="flex_gpu_estep_free")
        import :: c_int
    end function
    integer(c_int) function c_cols_lookup(lookup, qhkl, ncol, lb, ub) &
        &bind(c, name="flex_gpu_cols_lookup")
        import :: c_int, c_ptr
        type(c_ptr),    value :: lookup, qhkl
        integer(c_int), value :: ncol
        integer(c_int)        :: lb(3), ub(3)
    end function
    integer(c_int) function c_cols_fused_batch(rotm, vmask, nb, maxc, nyq_disk, iwinsz, &
        &lb, ub, deadskip, debias, pf2, tiny_val) bind(c, name="flex_gpu_cols_fused_batch")
        import :: c_int, c_ptr, c_float
        type(c_ptr),    value :: rotm, vmask
        integer(c_int), value :: nb, maxc, nyq_disk, iwinsz, deadskip, debias
        integer(c_int)        :: lb(3), ub(3)
        real(c_float),  value :: pf2, tiny_val
    end function
    integer(c_int) function c_cols_test_upload(plc, plct, nrec, hs_lo, hs_hi, ks_lo, khi) &
        &bind(c, name="flex_gpu_cols_test_upload")
        import :: c_int, c_ptr
        type(c_ptr),    value :: plc, plct
        integer(c_int), value :: nrec, hs_lo, hs_hi, ks_lo, khi
    end function
    integer(c_int) function c_prep_begin(lmsk, cs2, n1, n1o, maxrec, mskrad, minlen, wbg, &
        &l_ctf) bind(c, name="flex_gpu_prep_begin")
        import :: c_int, c_ptr, c_float
        type(c_ptr),    value :: lmsk, cs2
        integer(c_int), value :: n1, n1o, maxrec, wbg, l_ctf
        real(c_float),  value :: mskrad, minlen
    end function
    integer(c_int) function c_prep_batch(raw, ctfp, psh, vmask, sig2, nsig, nrec, hs_lo, &
        &hs_hi, ks_lo, nyq_pd) bind(c, name="flex_gpu_prep_batch")
        import :: c_int, c_ptr
        type(c_ptr),    value :: raw, ctfp, psh, vmask, sig2
        integer(c_int), value :: nsig, nrec, hs_lo, hs_hi, ks_lo, nyq_pd
    end function
    integer(c_int) function c_prep_free() bind(c, name="flex_gpu_prep_free")
        import :: c_int
    end function
    integer(c_int) function c_prep_fetch(plc, plct, nrec, hs_lo, hs_hi, ks_lo) &
        &bind(c, name="flex_gpu_prep_fetch")
        import :: c_int, c_ptr
        type(c_ptr),    value :: plc, plct
        integer(c_int), value :: nrec, hs_lo, hs_hi, ks_lo
    end function
    integer(c_int) function c_prep_fetch_sep(plcy, plt, plct, nrec, hs_lo, hs_hi, ks_lo) &
        &bind(c, name="flex_gpu_prep_fetch_sep")
        import :: c_int, c_ptr
        type(c_ptr),    value :: plcy, plt, plct
        integer(c_int), value :: nrec, hs_lo, hs_hi, ks_lo
    end function
    integer(c_int) function c_estep_batch_res(rotm, nrec, hs_lo, hs_hi, ks_lo, khi, nyq2, &
        &stats_host) bind(c, name="flex_gpu_estep_batch_resident")
        import :: c_int, c_ptr
        type(c_ptr),    value :: rotm, stats_host
        integer(c_int), value :: nrec, hs_lo, hs_hi, ks_lo, khi, nyq2
    end function
    integer(c_int) function c_coupled_batch_banked_res(cav, sav, dscale, rscale, vmask, &
        &dirid, nrec, ncomp, nrho, lb, ub) bind(c, name="flex_gpu_coupled_batch_banked_res")
        import :: c_int, c_ptr
        type(c_ptr),    value :: cav, sav, dscale, rscale, vmask
        integer(c_int)        :: dirid(*)
        integer(c_int), value :: nrec, ncomp, nrho
        integer(c_int)        :: lb(3), ub(3)
    end function
    integer(c_int) function c_poles_begin(rad, cs, sn, sqwq, ring_of, nang, nsamp, nk) &
        &bind(c, name="flex_gpu_poles_begin")
        import :: c_int, c_ptr
        type(c_ptr),    value :: rad, cs, sn, sqwq
        integer(c_int)        :: ring_of(*), nang(*)
        integer(c_int), value :: nsamp, nk
    end function
    integer(c_int) function c_poles_bank(Us, Cf, Cm0, c00, nc, ndir) &
        &bind(c, name="flex_gpu_poles_bank")
        import :: c_int, c_ptr
        type(c_ptr),    value :: Us, Cf, Cm0, c00
        integer(c_int), value :: nc, ndir
    end function
    integer(c_int) function c_poles_batch(ply, plt, rotm, cav, sav, dirid, vmask, nrec, &
        &hs_lo, hs_hi, ks_lo, nyq2_ex, stats_host) bind(c, name="flex_gpu_poles_batch")
        import :: c_int, c_ptr
        type(c_ptr),    value :: ply, plt, rotm, cav, sav, vmask, stats_host
        integer(c_int)        :: dirid(*)
        integer(c_int), value :: nrec, hs_lo, hs_hi, ks_lo, nyq2_ex
    end function
    integer(c_int) function c_poles_free() bind(c, name="flex_gpu_poles_free")
        import :: c_int
    end function
end interface
#endif

integer    :: g_lb(3) = 0, g_ub(3) = 0, g_ncomp = 0, g_lfny = 0
integer(8) :: g_nvox = 0
integer    :: g2_lb(3) = 0, g2_ub(3) = 0, g2_ncomp = 0, g2_nrho = 0, g2_lfny = 0
integer(8) :: g2_nvox = 0
integer    :: gc_lb(3) = 0, gc_ub(3) = 0, gc_ncol = 0
integer(8) :: gc_nvox = 0
integer    :: gs_lb(3) = 0, gs_ub(3) = 0
integer(8) :: gs_nvox = 0
logical    :: l_prep_dev_ready = .false.   ! device prep begun and not yet freed

contains

    !> true when the device prep lifecycle is active in this process (begin called, not freed);
    !! the shared stage prep funnel branches on this
    logical function flex_gpu_prep_ready() result( ok )
        ok = l_prep_dev_ready
    end function flex_gpu_prep_ready

    logical function flex_gpu_available() result( ok )
#ifdef USE_FLEX_CUDA
        ok = flex_gpu_device_count() > 0
#else
        ok = .false.
#endif
    end function flex_gpu_available

    !> allocate + zero the device accumulators for ncomp components shaped like recs(1)%cmat_exp
    subroutine flex_gpu_insert_begin_f( recs, ncomp )
        type(reconstructor), intent(in) :: recs(:)
        integer,             intent(in) :: ncomp
#ifdef USE_FLEX_CUDA
        integer :: ierr
        if( .not. allocated(recs(1)%cmat_exp) ) THROW_HARD('cmat_exp not allocated; flex_gpu_insert_begin_f')
        g_lb    = lbound(recs(1)%cmat_exp)
        g_ub    = ubound(recs(1)%cmat_exp)
        g_ncomp = ncomp
        g_lfny  = recs(1)%get_lfny(1)
        g_nvox  = int(g_ub(1)-g_lb(1)+1,8) * int(g_ub(2)-g_lb(2)+1,8) * int(g_ub(3)-g_lb(3)+1,8)
        ierr = c_insert_begin(int(ncomp,c_int), int(g_nvox,c_long))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_insert_begin failed')
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_insert_begin_f')
#endif
    end subroutine flex_gpu_insert_begin_f

    !> pack one batch and accumulate it on the device; mirrors the CPU batch front-end exactly
    subroutine flex_gpu_insert_batch_f( se, orientations, fpls, data_scales, density_scales, &
        &valid, nrecords )
        use simple_math, only: ceil_div, floor_div
        class(sym),        intent(inout) :: se
        type(ori),         intent(inout) :: orientations(:)
        type(fplane_type), intent(in)    :: fpls(:)
        real(dp),          intent(in)    :: data_scales(:,:), density_scales(:,:)
        logical,           intent(in)    :: valid(:)
        integer,           intent(in)    :: nrecords
#ifdef USE_FLEX_CUDA
        type(ori) :: o_sym
        complex(sp), allocatable, target :: pl_c(:,:,:)
        real(sp),    allocatable, target :: pl_ct(:,:,:), rotm(:,:,:,:), dsc(:,:), rsc(:,:)
        integer(c_int), allocatable, target :: act(:,:), nact(:), nyqd(:)
        integer :: frl(3,2), i, j, q, na, nsym, isym, nkeep, nyq_eff, ierr
        integer :: hpd_lo, hpd_hi, kpd_lo, hlo, hhi, klo, khi
        integer(c_int) :: lb_c(3), ub_c(3)
        integer, allocatable :: keep(:)
        if( g_ncomp <= 0 ) THROW_HARD('flex_gpu_insert_begin_f not called; flex_gpu_insert_batch_f')
        ! compact to valid records with at least one live component
        allocate(keep(nrecords), source=0)
        nkeep = 0
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            na = 0
            do q = 1, g_ncomp
                if( real(data_scales(q,i)) /= 0. .or. real(max(0.d0,density_scales(q,i))) /= 0. ) na = na + 1
            end do
            if( na == 0 ) cycle
            nkeep = nkeep + 1
            keep(nkeep) = i
        end do
        if( nkeep == 0 ) return
        ! common plane geometry (asserted): all planes in a batch share frlims
        frl = fpls(keep(1))%frlims
        do j = 2, nkeep
            if( any(fpls(keep(j))%frlims /= frl) ) THROW_HARD('heterogeneous frlims in batch; flex_gpu_insert_batch_f')
        end do
        hpd_lo = frl(1,1); hpd_hi = frl(1,2); kpd_lo = frl(2,1)
        hlo = ceil_div (frl(1,1), OSMPL_PAD_FAC); hhi = floor_div(frl(1,2), OSMPL_PAD_FAC)
        klo = ceil_div (frl(2,1), OSMPL_PAD_FAC); khi = floor_div(frl(2,2), OSMPL_PAD_FAC)
        nsym = se%get_nsym()
        ! PACKED upload: the kernels only ever read pf-multiples of the padded storage, so ship
        ! the strided section — 1/pf^2 of the bytes over the 4-worker-contended PCIe link
        allocate(pl_c (hlo:hhi, klo:0, nkeep), source=cmplx(0.,0.))
        allocate(pl_ct(hlo:hhi, klo:0, nkeep), source=0.)
        allocate(rotm(3,3,nsym,nkeep), source=0.)
        allocate(dsc(g_ncomp,nkeep), rsc(g_ncomp,nkeep), source=0.)
        allocate(act(g_ncomp,nkeep), source=0_c_int)
        allocate(nact(nkeep), nyqd(nkeep), source=0_c_int)
        !$omp parallel do default(shared) private(j,i) schedule(static) proc_bind(close)
        do j = 1, nkeep
            i = keep(j)
            pl_c (:,:,j) = fpls(i)%cmplx_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)
            pl_ct(:,:,j) = fpls(i)%ctfsq_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)
        end do
        !$omp end parallel do
        do j = 1, nkeep
            i = keep(j)
            na = 0
            do q = 1, g_ncomp
                dsc(q,j) = real(data_scales(q,i))
                rsc(q,j) = real(max(0.d0, density_scales(q,i)))
                if( dsc(q,j) /= 0. .or. rsc(q,j) /= 0. )then
                    na = na + 1
                    act(na,j) = int(q, c_int)
                endif
            end do
            nact(j) = int(na, c_int)
            rotm(:,:,1,j) = orientations(i)%get_mat()
            do isym = 2, nsym
                call se%apply(orientations(i), isym, o_sym)
                rotm(:,:,isym,j) = o_sym%get_mat()
            end do
            nyq_eff = g_lfny
            if( fpls(i)%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpls(i)%nyq / OSMPL_PAD_FAC))
            nyqd(j) = int(nyq_eff * (nyq_eff + 1), c_int)
        end do
        call o_sym%kill
        lb_c = int(g_lb, c_int); ub_c = int(g_ub, c_int)
        ierr = c_insert_batch(c_loc(pl_c), c_loc(pl_ct), c_loc(rotm), c_loc(dsc), c_loc(rsc), &
            &c_loc(act), c_loc(nact), c_loc(nyqd), int(nkeep,c_int), int(g_ncomp,c_int), &
            &int(nsym,c_int), int(hlo,c_int), int(hhi,c_int), int(klo,c_int), &
            &int(hlo,c_int), int(hhi,c_int), int(klo,c_int), int(khi,c_int), &
            &int(OSMPL_PAD_FAC,c_int), lb_c, ub_c, real(TINY,c_float))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_insert_batch failed')
        deallocate(pl_c, pl_ct, rotm, dsc, rsc, act, nact, nyqd, keep)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_insert_batch_f')
#endif
    end subroutine flex_gpu_insert_batch_f

    !> resident-plane insertion: consumes the packed planes the device prep left in place
    !! (zero plane traffic). Records are NOT compacted -- dead records carry nact=0. The
    !! per-record nyq disk is uniform, from the shared plane band.
    subroutine flex_gpu_insert_batch_res_f( se, orientations, data_scales, density_scales, &
        &valid, nrecords, nyq_unpd )
        class(sym),        intent(inout) :: se
        type(ori),         intent(inout) :: orientations(:)
        real(dp),          intent(in)    :: data_scales(:,:), density_scales(:,:)
        logical,           intent(in)    :: valid(:)
        integer,           intent(in)    :: nrecords, nyq_unpd
#ifdef USE_FLEX_CUDA
        type(ori) :: o_sym
        real(sp),       allocatable, target :: rotm(:,:,:,:), dsc(:,:), rsc(:,:)
        integer(c_int), allocatable, target :: act(:,:), nact(:), nyqd(:)
        integer :: i, q, na, nsym, isym, nyq_eff, ierr
        integer(c_int) :: lb_c(3), ub_c(3)
        if( g_ncomp <= 0 ) THROW_HARD('flex_gpu_insert_begin_f not called; flex_gpu_insert_batch_res_f')
        nsym = se%get_nsym()
        allocate(rotm(3,3,nsym,nrecords), source=0.)
        allocate(dsc(g_ncomp,nrecords), rsc(g_ncomp,nrecords), source=0.)
        allocate(act(g_ncomp,nrecords), source=0_c_int)
        allocate(nact(nrecords), nyqd(nrecords), source=0_c_int)
        nyq_eff = min(g_lfny, max(1, nyq_unpd))
        nyqd    = int(nyq_eff * (nyq_eff + 1), c_int)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            na = 0
            do q = 1, g_ncomp
                dsc(q,i) = real(data_scales(q,i))
                rsc(q,i) = real(max(0.d0, density_scales(q,i)))
                if( dsc(q,i) /= 0. .or. rsc(q,i) /= 0. )then
                    na = na + 1
                    act(na,i) = int(q, c_int)
                endif
            end do
            nact(i) = int(na, c_int)
            if( na == 0 ) cycle
            rotm(:,:,1,i) = orientations(i)%get_mat()
            do isym = 2, nsym
                call se%apply(orientations(i), isym, o_sym)
                rotm(:,:,isym,i) = o_sym%get_mat()
            end do
        end do
        call o_sym%kill
        lb_c = int(g_lb, c_int); ub_c = int(g_ub, c_int)
        ierr = c_insert_batch_res(c_loc(rotm), c_loc(dsc), c_loc(rsc), c_loc(act), &
            &c_loc(nact), c_loc(nyqd), int(nrecords,c_int), int(g_ncomp,c_int), &
            &int(nsym,c_int), lb_c, ub_c, real(TINY,c_float))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_insert_batch_resident failed')
        deallocate(rotm, dsc, rsc, act, nact, nyqd)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_insert_batch_res_f')
#endif
    end subroutine flex_gpu_insert_batch_res_f

    !> fetch the device totals and ADD into the reconstructors, then release the device buffers.
    !! Components 1..size(recs) land in recs; with recs2 present, size(recs)+1.. land in recs2
    !! (the fused even/odd layout of the state stage: one device accumulation, two target arrays).
    subroutine flex_gpu_insert_end_f( recs, recs2 )
        type(reconstructor),           intent(inout) :: recs(:)
        type(reconstructor), optional, intent(inout) :: recs2(:)
#ifdef USE_FLEX_CUDA
        complex(sp), allocatable, target :: cbuf(:,:,:)
        real(sp),    allocatable, target :: rbuf(:,:,:)
        integer :: q, n1, ierr
        if( g_ncomp <= 0 ) return
        n1 = size(recs)
        if( g_ncomp > n1 .and. .not. present(recs2) ) &
            &THROW_HARD('component count exceeds recs and no recs2; flex_gpu_insert_end_f')
        allocate(cbuf(g_lb(1):g_ub(1), g_lb(2):g_ub(2), g_lb(3):g_ub(3)))
        allocate(rbuf(g_lb(1):g_ub(1), g_lb(2):g_ub(2), g_lb(3):g_ub(3)))
        do q = 1, g_ncomp
            ierr = c_insert_fetch(int(q,c_int), c_loc(cbuf), c_loc(rbuf))
            if( ierr /= 0 ) THROW_HARD('flex_gpu_insert_fetch failed')
            if( q <= n1 )then
                recs(q)%cmat_exp = recs(q)%cmat_exp + cbuf
                recs(q)%rho_exp  = recs(q)%rho_exp  + rbuf
            else
                recs2(q-n1)%cmat_exp = recs2(q-n1)%cmat_exp + cbuf
                recs2(q-n1)%rho_exp  = recs2(q-n1)%rho_exp  + rbuf
            endif
        end do
        ierr = c_insert_free()
        deallocate(cbuf, rbuf)
        g_ncomp = 0; g_nvox = 0
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_insert_end_f')
#endif
    end subroutine flex_gpu_insert_end_f

    !> allocate + zero the coupled-stats device accumulators (probe M-step shape):
    !! ncomp cmplx components + nrho density slots (1 shared / ncomp diagonal / npairs packed)
    subroutine flex_gpu_coupled_begin_f( recs, ncomp, nrho )
        type(reconstructor), intent(inout) :: recs(:)   ! inout only because get_kbwin is not intent(in)-safe
        integer,             intent(in)    :: ncomp, nrho
#ifdef USE_FLEX_CUDA
        type(kbinterpol) :: kbw
        integer :: ierr
        if( .not. allocated(recs(1)%cmat_exp) ) THROW_HARD('cmat_exp not allocated; flex_gpu_coupled_begin_f')
        ! the coupled CPU path uses the reconstructor's OWN window; the device polynomial is the
        ! (Whalf=1.5, alpha=2) fast path, so anything else must refuse rather than drift
        kbw = recs(1)%get_kbwin()
        if( abs(kbw%get_winsz() - 1.5) > 1.e-6 .or. abs(kbw%get_alpha() - 2.0) > 1.e-6 ) &
            &THROW_HARD('reconstructor kb window is not (1.5, 2.0); GPU coupled path unsupported')
        g2_lb    = lbound(recs(1)%cmat_exp)
        g2_ub    = ubound(recs(1)%cmat_exp)
        g2_ncomp = ncomp
        g2_nrho  = nrho
        g2_lfny  = recs(1)%get_lfny(1)
        g2_nvox  = int(g2_ub(1)-g2_lb(1)+1,8) * int(g2_ub(2)-g2_lb(2)+1,8) * int(g2_ub(3)-g2_lb(3)+1,8)
        ierr = c_coupled_begin(int(ncomp,c_int), int(nrho,c_int), int(g2_nvox,c_long))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_coupled_begin failed')
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_coupled_begin_f')
#endif
    end subroutine flex_gpu_coupled_begin_f

    !> pack one batch for the coupled kernel: the plane value is conj(transfer)*cmplx built here,
    !! the density scales are the packed upper triangle of dens (the probe's mode)
    subroutine flex_gpu_coupled_batch_f( se, orientations, fpls, data_scales, dens, valid, nrecords )
        class(sym),        intent(inout) :: se
        type(ori),         intent(inout) :: orientations(:)
        type(fplane_type), intent(in)    :: fpls(:)
        real(dp),          intent(in)    :: data_scales(:,:), dens(:,:,:)
        logical,           intent(in)    :: valid(:)
        integer,           intent(in)    :: nrecords
#ifdef USE_FLEX_CUDA
        real(dp), allocatable :: dsc(:,:), rsc(:,:)
        integer :: i, q, r
        if( g2_ncomp <= 0 ) THROW_HARD('flex_gpu_coupled_begin_f not called; flex_gpu_coupled_batch_f')
        allocate(dsc(g2_ncomp,nrecords), rsc(g2_nrho,nrecords), source=0.d0)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            dsc(1:g2_ncomp,i) = data_scales(1:g2_ncomp,i)
            do r = 1, g2_ncomp
                do q = 1, r
                    rsc((r*(r-1))/2 + q, i) = dens(q,r,i)
                end do
            end do
        end do
        call flex_gpu_coupled_batch_raw_f(se, orientations, fpls, dsc, rsc, valid, nrecords)
        deallocate(dsc, rsc)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_coupled_batch_f')
#endif
    end subroutine flex_gpu_coupled_batch_f

    !> as coupled_batch_f but with the per-record component/density scale slots already built by
    !! the caller — the combined even/odd layouts route halfsets through these slots
    subroutine flex_gpu_coupled_batch_raw_f( se, orientations, fpls, dscales, rscales, valid, nrecords )
        use simple_math, only: ceil_div, floor_div
        class(sym),        intent(inout) :: se
        type(ori),         intent(inout) :: orientations(:)
        type(fplane_type), intent(in)    :: fpls(:)
        real(dp),          intent(in)    :: dscales(:,:), rscales(:,:)
        logical,           intent(in)    :: valid(:)
        integer,           intent(in)    :: nrecords
#ifdef USE_FLEX_CUDA
        type(ori) :: o_sym
        complex(sp), allocatable, target :: pl_c(:,:,:)
        real(sp),    allocatable, target :: pl_ct(:,:,:), rotm(:,:,:,:), dsc(:,:), rsc(:,:)
        integer(c_int),    allocatable, target :: nyqd(:)
        integer(c_int8_t), allocatable, target :: vmask(:)
        integer :: frl(3,2), i, nsym, isym, nyq_eff, ierr, iref
        integer :: hpd_lo, hpd_hi, kpd_lo, hlo, hhi, klo, khi
        integer(c_int) :: lb_c(3), ub_c(3)
        if( g2_ncomp <= 0 ) THROW_HARD('flex_gpu_coupled_begin_f not called; flex_gpu_coupled_batch_raw_f')
        iref = 0
        do i = 1, nrecords
            if( valid(i) )then
                if( .not. allocated(fpls(i)%transfer_plane) ) &
                    &THROW_HARD('forward transfer plane does not exist; flex_gpu_coupled_batch_raw_f')
                if( iref == 0 ) iref = i
            endif
        end do
        if( iref == 0 ) return
        frl = fpls(iref)%frlims
        do i = 1, nrecords
            if( valid(i) .and. any(fpls(i)%frlims /= frl) ) &
                &THROW_HARD('heterogeneous frlims in batch; flex_gpu_coupled_batch_raw_f')
        end do
        hpd_lo = frl(1,1); hpd_hi = frl(1,2); kpd_lo = frl(2,1)
        hlo = ceil_div (frl(1,1), OSMPL_PAD_FAC); hhi = floor_div(frl(1,2), OSMPL_PAD_FAC)
        klo = ceil_div (frl(2,1), OSMPL_PAD_FAC); khi = floor_div(frl(2,2), OSMPL_PAD_FAC)
        nsym = se%get_nsym()
        allocate(pl_c (hlo:hhi, klo:0, nrecords), source=cmplx(0.,0.))
        allocate(pl_ct(hlo:hhi, klo:0, nrecords), source=0.)
        allocate(rotm(3,3,nsym,nrecords), source=0.)
        allocate(dsc(g2_ncomp,nrecords), rsc(g2_nrho,nrecords), source=0.)
        allocate(nyqd(nrecords), source=0_c_int)
        allocate(vmask(nrecords), source=0_c_int8_t)
        !$omp parallel do default(shared) private(i,nyq_eff) schedule(static) proc_bind(close)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            vmask(i) = 1_c_int8_t
            ! packed pf-multiple section (see flex_gpu_insert_batch_f)
            pl_c (:,:,i) = conjg(fpls(i)%transfer_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)) * &
                &fpls(i)%cmplx_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)
            pl_ct(:,:,i) = fpls(i)%ctfsq_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)
            dsc(:,i) = real(dscales(1:g2_ncomp,i))
            rsc(:,i) = real(rscales(1:g2_nrho,i))
            rotm(:,:,1,i) = orientations(i)%get_mat()
            nyq_eff = g2_lfny
            if( fpls(i)%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpls(i)%nyq / OSMPL_PAD_FAC))
            nyqd(i) = int(nyq_eff * (nyq_eff + 1), c_int)
        end do
        !$omp end parallel do
        if( nsym > 1 )then
            ! symmetry copies stay serial: se%apply mutates the shared o_sym scratch
            do i = 1, nrecords
                if( .not. valid(i) ) cycle
                do isym = 2, nsym
                    call se%apply(orientations(i), isym, o_sym)
                    rotm(:,:,isym,i) = o_sym%get_mat()
                end do
            end do
        endif
        call o_sym%kill
        lb_c = int(g2_lb, c_int); ub_c = int(g2_ub, c_int)
        ierr = c_coupled_batch(c_loc(pl_c), c_loc(pl_ct), c_loc(rotm), c_loc(dsc), c_loc(rsc), &
            &c_loc(vmask), c_loc(nyqd), int(nrecords,c_int), int(g2_ncomp,c_int), &
            &int(g2_nrho,c_int), int(nsym,c_int), int(hlo,c_int), int(hhi,c_int), &
            &int(klo,c_int), int(hlo,c_int), int(hhi,c_int), int(klo,c_int), int(khi,c_int), &
            &int(OSMPL_PAD_FAC,c_int), lb_c, ub_c, real(TINY,c_float))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_coupled_batch failed')
        deallocate(pl_c, pl_ct, rotm, dsc, rsc, nyqd, vmask)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_coupled_batch_raw_f')
#endif
    end subroutine flex_gpu_coupled_batch_raw_f

    !> upload the direction bank for the banked coupled adjoint (one call per probe stage)
    subroutine flex_gpu_coupled_bank_f( rmat_bank, ndir )
        integer, intent(in) :: ndir
        real,    intent(in) :: rmat_bank(3,3,ndir)
#ifdef USE_FLEX_CUDA
        real(sp), allocatable, target :: rmb(:,:,:)
        integer :: ierr
        allocate(rmb(3,3,ndir), source=real(rmat_bank,sp))
        ierr = c_coupled_bank(c_loc(rmb), int(ndir,c_int))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_coupled_bank failed')
        deallocate(rmb)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_coupled_bank_f')
#endif
    end subroutine flex_gpu_coupled_bank_f

    subroutine flex_gpu_coupled_bank_free_f
#ifdef USE_FLEX_CUDA
        integer :: ierr
        ierr = c_coupled_bank_free()
#endif
    end subroutine flex_gpu_coupled_bank_free_f

    !> banked coupled batch: records MUST be sorted by dirid; geometry comes from the bank
    !! (dirid + in-plane cos/sin per record), so no orientations/symmetry are taken — the
    !! caller gates this path to nsym == 1.
    subroutine flex_gpu_coupled_batch_banked_f( fpls, dscales, rscales, valid, nrecords, &
        &dirid, cav, sav, premultiplied )
        use simple_math, only: ceil_div, floor_div
        type(fplane_type), intent(in) :: fpls(:)
        real(dp),          intent(in) :: dscales(:,:), rscales(:,:)
        logical,           intent(in) :: valid(:)
        integer,           intent(in) :: nrecords
        integer,           intent(in) :: dirid(:)
        real,              intent(in) :: cav(:), sav(:)
        logical, optional, intent(in) :: premultiplied
#ifdef USE_FLEX_CUDA
        complex(sp), allocatable, target :: pl_c(:,:,:)
        real(sp),    allocatable, target :: pl_ct(:,:,:), dsc(:,:), rsc(:,:), ca_c(:), sa_c(:)
        integer(c_int),    allocatable, target :: did(:)
        integer(c_int8_t), allocatable, target :: vmask(:)
        integer :: frl(3,2), i, nyq_eff, ierr, iref
        integer :: hpd_lo, hpd_hi, kpd_lo, hlo, hhi, klo, khi
        integer(c_int) :: lb_c(3), ub_c(3)
        logical :: l_premul
        l_premul = .false.
        if( present(premultiplied) ) l_premul = premultiplied
        if( g2_ncomp <= 0 ) THROW_HARD('flex_gpu_coupled_begin_f not called; flex_gpu_coupled_batch_banked_f')
        iref = 0
        do i = 1, nrecords
            if( valid(i) )then
                if( .not. allocated(fpls(i)%transfer_plane) ) &
                    &THROW_HARD('forward transfer plane does not exist; flex_gpu_coupled_batch_banked_f')
                if( iref == 0 ) iref = i
            endif
        end do
        if( iref == 0 ) return
        frl = fpls(iref)%frlims
        do i = 1, nrecords
            if( valid(i) .and. any(fpls(i)%frlims /= frl) ) &
                &THROW_HARD('heterogeneous frlims in batch; flex_gpu_coupled_batch_banked_f')
        end do
        hpd_lo = frl(1,1); hpd_hi = frl(1,2); kpd_lo = frl(2,1)
        hlo = ceil_div (frl(1,1), OSMPL_PAD_FAC); hhi = floor_div(frl(1,2), OSMPL_PAD_FAC)
        klo = ceil_div (frl(2,1), OSMPL_PAD_FAC); khi = floor_div(frl(2,2), OSMPL_PAD_FAC)
        allocate(pl_c (hlo:hhi, klo:0, nrecords), source=cmplx(0.,0.))
        allocate(pl_ct(hlo:hhi, klo:0, nrecords), source=0.)
        allocate(dsc(g2_ncomp,nrecords), rsc(g2_nrho,nrecords), source=0.)
        allocate(ca_c(nrecords), source=1._sp)
        allocate(sa_c(nrecords), source=0._sp)
        allocate(did(nrecords), source=1_c_int)
        allocate(vmask(nrecords), source=0_c_int8_t)
        nyq_eff = g2_lfny
        if( fpls(iref)%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpls(iref)%nyq / OSMPL_PAD_FAC))
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, nrecords
            did(i) = int(max(1, dirid(i)), c_int)
            if( .not. valid(i) ) cycle
            vmask(i) = 1_c_int8_t
            ! packed pf-multiple section (see flex_gpu_insert_batch_f); with l_premul the
            ! cmplx plane already IS conj(T)y (residual formed on device by the fused E-step)
            if( l_premul )then
                pl_c(:,:,i) = fpls(i)%cmplx_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                    &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)
            else
                pl_c(:,:,i) = conjg(fpls(i)%transfer_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                    &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)) * &
                    &fpls(i)%cmplx_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                    &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)
            endif
            pl_ct(:,:,i) = fpls(i)%ctfsq_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)
            dsc(:,i) = real(dscales(1:g2_ncomp,i))
            rsc(:,i) = real(rscales(1:g2_nrho,i))
            ca_c(i)  = cav(i)
            sa_c(i)  = sav(i)
        end do
        !$omp end parallel do
        lb_c = int(g2_lb, c_int); ub_c = int(g2_ub, c_int)
        ierr = c_coupled_batch_banked(c_loc(pl_c), c_loc(pl_ct), c_loc(ca_c), c_loc(sa_c), &
            &c_loc(dsc), c_loc(rsc), c_loc(vmask), did, int(nrecords,c_int), &
            &int(g2_ncomp,c_int), int(g2_nrho,c_int), int(hlo,c_int), int(hhi,c_int), &
            &int(klo,c_int), int(hlo,c_int), int(hhi,c_int), int(klo,c_int), int(khi,c_int), &
            &int(OSMPL_PAD_FAC,c_int), int(nyq_eff*(nyq_eff+1),c_int), lb_c, ub_c)
        if( ierr /= 0 ) THROW_HARD('flex_gpu_coupled_batch_banked failed')
        deallocate(pl_c, pl_ct, dsc, rsc, ca_c, sa_c, did, vmask)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_coupled_batch_banked_f')
#endif
    end subroutine flex_gpu_coupled_batch_banked_f

    !> fetch coupled totals: cmat ADDS into the reconstructors, rho ADDS into the
    !! (nrho, shape) cross-density array, then release. With recs2/rho_cross2 present the
    !! combined even/odd layout splits at size(recs) components and size(rho_cross,1) slots.
    subroutine flex_gpu_coupled_end_f( recs, rho_cross, recs2, rho_cross2 )
        type(reconstructor),           intent(inout) :: recs(:)
        real(sp),                      intent(inout) :: rho_cross(:,:,:,:)
        type(reconstructor), optional, intent(inout) :: recs2(:)
        real(sp),  optional,           intent(inout) :: rho_cross2(:,:,:,:)
#ifdef USE_FLEX_CUDA
        complex(sp), allocatable, target :: cbuf(:,:,:)
        real(sp),    allocatable, target :: rbuf(:,:,:,:)
        integer :: q, n1, nr1, ierr
        if( g2_ncomp <= 0 ) return
        n1  = size(recs)
        nr1 = size(rho_cross,1)
        if( g2_ncomp > n1 .and. .not. present(recs2) ) &
            &THROW_HARD('component count exceeds recs and no recs2; flex_gpu_coupled_end_f')
        if( g2_nrho > nr1 .and. .not. present(rho_cross2) ) &
            &THROW_HARD('density count exceeds rho_cross and no rho_cross2; flex_gpu_coupled_end_f')
        allocate(cbuf(g2_lb(1):g2_ub(1), g2_lb(2):g2_ub(2), g2_lb(3):g2_ub(3)))
        do q = 1, g2_ncomp
            ierr = c_coupled_fetch_cmat(int(q,c_int), c_loc(cbuf))
            if( ierr /= 0 ) THROW_HARD('flex_gpu_coupled_fetch_cmat failed')
            if( q <= n1 )then
                recs(q)%cmat_exp = recs(q)%cmat_exp + cbuf
            else
                recs2(q-n1)%cmat_exp = recs2(q-n1)%cmat_exp + cbuf
            endif
        end do
        allocate(rbuf(g2_nrho, g2_ub(1)-g2_lb(1)+1, g2_ub(2)-g2_lb(2)+1, g2_ub(3)-g2_lb(3)+1))
        ierr = c_coupled_fetch_rho(c_loc(rbuf))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_coupled_fetch_rho failed')
        rho_cross(1:min(g2_nrho,nr1), 1:size(rbuf,2), 1:size(rbuf,3), 1:size(rbuf,4)) = &
            &rho_cross(1:min(g2_nrho,nr1), 1:size(rbuf,2), 1:size(rbuf,3), 1:size(rbuf,4)) + &
            &rbuf(1:min(g2_nrho,nr1),:,:,:)
        if( g2_nrho > nr1 )then
            rho_cross2(1:g2_nrho-nr1, 1:size(rbuf,2), 1:size(rbuf,3), 1:size(rbuf,4)) = &
                &rho_cross2(1:g2_nrho-nr1, 1:size(rbuf,2), 1:size(rbuf,3), 1:size(rbuf,4)) + &
                &rbuf(nr1+1:g2_nrho,:,:,:)
        endif
        ierr = c_coupled_free()
        deallocate(cbuf, rbuf)
        g2_ncomp = 0; g2_nrho = 0; g2_nvox = 0
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_coupled_end_f')
#endif
    end subroutine flex_gpu_coupled_end_f

    !> columns phase-B device accumulators: ncol lattices per halfset, zeroed
    subroutine flex_gpu_cols_begin_f( ncol, lb, ub )
        integer, intent(in) :: ncol, lb(3), ub(3)
#ifdef USE_FLEX_CUDA
        integer :: ierr
        gc_lb = lb; gc_ub = ub; gc_ncol = ncol
        gc_nvox = int(ub(1)-lb(1)+1,8) * int(ub(2)-lb(2)+1,8) * int(ub(3)-lb(3)+1,8)
        ierr = c_cols_begin(int(ncol,c_int), int(gc_nvox,c_long))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_cols_begin failed')
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_cols_begin_f')
#endif
    end subroutine flex_gpu_cols_begin_f

    !> one phase-B block: the CPU-built sample caches and column values pass through verbatim
    subroutine flex_gpu_cols_batch_f( cwin_b, cwx_b, cwy_b, cwz_b, cpl_b, cct_b, ncache_b, &
        &valid_b, eo_b, gcol_b, hcolv_b, col_hkl, nb, maxc, l_deadskip, l_debias )
        integer, intent(in), target :: cwin_b(3,maxc,nb)
        real,    intent(in), target :: cwx_b(3,maxc,nb), cwy_b(3,maxc,nb), cwz_b(3,maxc,nb)
        complex, intent(in), target :: cpl_b(maxc,nb)
        real,    intent(in), target :: cct_b(maxc,nb)
        integer, intent(in)         :: ncache_b(nb), eo_b(nb)
        logical, intent(in)         :: valid_b(nb)
        complex, intent(in), target :: gcol_b(:,:)
        real,    intent(in), target :: hcolv_b(:,:)
        integer, intent(in), target :: col_hkl(:,:)
        integer, intent(in)         :: nb, maxc
        logical, intent(in)         :: l_deadskip, l_debias
#ifdef USE_FLEX_CUDA
        integer(c_int),    allocatable, target :: nc_c(:)
        integer(c_int8_t), allocatable, target :: vmask(:)
        integer(c_int) :: lb_c(3), ub_c(3)
        integer :: i, ierr
        if( gc_ncol <= 0 ) THROW_HARD('flex_gpu_cols_begin_f not called; flex_gpu_cols_batch_f')
        allocate(nc_c(nb), vmask(nb))
        do i = 1, nb
            nc_c(i)  = int(ncache_b(i), c_int)
            vmask(i) = 0_c_int8_t
            if( valid_b(i) ) vmask(i) = int(merge(1, 2, eo_b(i) == 0), c_int8_t)
        end do
        lb_c = int(gc_lb, c_int); ub_c = int(gc_ub, c_int)
        ierr = c_cols_batch(c_loc(cwin_b), c_loc(cwx_b), c_loc(cwy_b), c_loc(cwz_b), &
            &c_loc(cpl_b), c_loc(cct_b), c_loc(nc_c), c_loc(vmask), c_loc(gcol_b), &
            &c_loc(hcolv_b), c_loc(col_hkl), int(nb,c_int), int(maxc,c_int), &
            &int(gc_ncol,c_int), lb_c, ub_c, merge(1_c_int, 0_c_int, l_deadskip), &
            &merge(1_c_int, 0_c_int, l_debias), real(OSMPL_PAD_FAC*OSMPL_PAD_FAC, c_float))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_cols_batch failed')
        deallocate(nc_c, vmask)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_cols_batch_f')
#endif
    end subroutine flex_gpu_cols_batch_f

    !> fused columns: upload the column lookup volume + coordinates once per stage
    subroutine flex_gpu_cols_lookup_f( col_lookup, col_hkl, ncol )
        integer, intent(in), target, contiguous :: col_lookup(:,:,:)
        integer, intent(in), target, contiguous :: col_hkl(:,:)
        integer, intent(in) :: ncol
#ifdef USE_FLEX_CUDA
        integer(c_int) :: lb_c(3), ub_c(3)
        integer :: ierr
        if( gc_ncol <= 0 ) THROW_HARD('flex_gpu_cols_begin_f not called; flex_gpu_cols_lookup_f')
        if( any(shape(col_lookup) /= gc_ub - gc_lb + 1) ) &
            &THROW_HARD('column lookup shape mismatch; flex_gpu_cols_lookup_f')
        lb_c = int(gc_lb, c_int); ub_c = int(gc_ub, c_int)
        ierr = c_cols_lookup(c_loc(col_lookup), c_loc(col_hkl), int(ncol,c_int), lb_c, ub_c)
        if( ierr /= 0 ) THROW_HARD('flex_gpu_cols_lookup failed')
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_cols_lookup_f')
#endif
    end subroutine flex_gpu_cols_lookup_f

    !> fused columns batch: device cache build + gather + scatter over the RESIDENT adjoint
    !! pair (after the E-step residual). Zero cache traffic; entry order inside a particle's
    !! cache is atomic-append nondeterministic, so gate deliverables statistically.
    subroutine flex_gpu_cols_fused_batch_f( orientations, valid_b, eo_b, nb, nyq_rec, maxc, &
        &l_deadskip, l_debias )
        type(ori),  intent(inout) :: orientations(:)
        logical,    intent(in)    :: valid_b(:)
        integer,    intent(in)    :: eo_b(:)
        integer,    intent(in)    :: nb, nyq_rec, maxc
        logical,    intent(in)    :: l_deadskip, l_debias
#ifdef USE_FLEX_CUDA
        real(sp),          allocatable, target :: rotm(:,:,:)
        integer(c_int8_t), allocatable, target :: vmask(:)
        integer(c_int) :: lb_c(3), ub_c(3)
        integer :: i, ierr, iwinsz
        if( gc_ncol <= 0 ) THROW_HARD('flex_gpu_cols_begin_f not called; flex_gpu_cols_fused_batch_f')
        allocate(rotm(3,3,nb), source=0.0_sp)
        allocate(vmask(nb), source=0_c_int8_t)
        do i = 1, nb
            if( valid_b(i) )then
                vmask(i)    = int(merge(1, 2, eo_b(i) == 0), c_int8_t)
                rotm(:,:,i) = orientations(i)%get_mat()
            endif
        end do
        iwinsz = ceiling(KBWINSZ - 0.5)
        lb_c = int(gc_lb, c_int); ub_c = int(gc_ub, c_int)
        ierr = c_cols_fused_batch(c_loc(rotm), c_loc(vmask), int(nb,c_int), int(maxc,c_int), &
            &int(nyq_rec*(nyq_rec+1),c_int), int(iwinsz,c_int), lb_c, ub_c, &
            &merge(1_c_int,0_c_int,l_deadskip), merge(1_c_int,0_c_int,l_debias), &
            &real(OSMPL_PAD_FAC*OSMPL_PAD_FAC,c_float), real(TINY,c_float))
        if( ierr == 5 ) THROW_HARD('fused columns cache overflow; flex_gpu_cols_fused_batch_f')
        if( ierr /= 0 ) THROW_HARD('flex_gpu_cols_fused_batch failed')
        deallocate(rotm, vmask)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_cols_fused_batch_f')
#endif
    end subroutine flex_gpu_cols_fused_batch_f

    !> unit-test helper: inject synthetic packed adjoint planes as the resident pair
    subroutine flex_gpu_cols_test_upload_f( plc, plct, nrec, hlo, hhi, klo, khi )
        complex(sp), intent(in), target, contiguous :: plc(:,:,:)
        real(sp),    intent(in), target, contiguous :: plct(:,:,:)
        integer,     intent(in) :: nrec, hlo, hhi, klo, khi
#ifdef USE_FLEX_CUDA
        integer :: ierr
        ierr = c_cols_test_upload(c_loc(plc), c_loc(plct), int(nrec,c_int), int(hlo,c_int), &
            &int(hhi,c_int), int(klo,c_int), int(khi,c_int))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_cols_test_upload failed')
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_cols_test_upload_f')
#endif
    end subroutine flex_gpu_cols_test_upload_f

    !> fetch the device column accumulators and ADD into the host arrays, then release
    subroutine flex_gpu_cols_end_f( Bcol_e, Hcol_e, Bcol_o, Hcol_o )
        complex, intent(inout) :: Bcol_e(:,:,:,:), Bcol_o(:,:,:,:)
        real,    intent(inout) :: Hcol_e(:,:,:,:), Hcol_o(:,:,:,:)
#ifdef USE_FLEX_CUDA
        complex(sp), allocatable, target :: Bbuf(:,:,:,:)
        real(sp),    allocatable, target :: Hbuf(:,:,:,:)
        integer :: ierr
        if( gc_ncol <= 0 ) return
        allocate(Bbuf(gc_ub(1)-gc_lb(1)+1, gc_ub(2)-gc_lb(2)+1, gc_ub(3)-gc_lb(3)+1, gc_ncol))
        allocate(Hbuf(gc_ub(1)-gc_lb(1)+1, gc_ub(2)-gc_lb(2)+1, gc_ub(3)-gc_lb(3)+1, gc_ncol))
        ierr = c_cols_fetch(0_c_int, c_loc(Bbuf), c_loc(Hbuf))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_cols_fetch even failed')
        Bcol_e = Bcol_e + Bbuf
        Hcol_e = Hcol_e + Hbuf
        ierr = c_cols_fetch(1_c_int, c_loc(Bbuf), c_loc(Hbuf))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_cols_fetch odd failed')
        Bcol_o = Bcol_o + Bbuf
        Hcol_o = Hcol_o + Hbuf
        ierr = c_cols_free()
        deallocate(Bbuf, Hbuf)
        gc_ncol = 0; gc_nvox = 0
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_cols_end_f')
#endif
    end subroutine flex_gpu_cols_end_f

    !> SNR device state: per-particle slab pool + the var/dens accumulators, zeroed
    subroutine flex_gpu_snr_begin_f( lb, ub, nslab )
        integer, intent(in) :: lb(3), ub(3), nslab
#ifdef USE_FLEX_CUDA
        integer :: ierr
        gs_lb = lb; gs_ub = ub
        gs_nvox = int(ub(1)-lb(1)+1,8) * int(ub(2)-lb(2)+1,8) * int(ub(3)-lb(3)+1,8)
        ierr = c_snr_begin(int(gs_nvox,c_long), int(nslab,c_int))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_snr_begin failed')
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_snr_begin_f')
#endif
    end subroutine flex_gpu_snr_begin_f

    !> one batch of packed adjoint-residual planes: insert-into-slab + square-sweep on device
    subroutine flex_gpu_snr_batch_f( pl_c, pl_ct, rotm, nyqd, vmask, nb, nsym, frl )
        use simple_math, only: ceil_div, floor_div
        complex,    intent(in), target :: pl_c(:,:,:)
        real,       intent(in), target :: pl_ct(:,:,:)
        real,       intent(in), target :: rotm(:,:,:,:)
        integer(1), intent(in), target :: vmask(:)
        integer,    intent(in), target :: nyqd(:)
        integer,    intent(in)         :: nb, nsym, frl(3,2)
#ifdef USE_FLEX_CUDA
        integer(c_int), allocatable, target :: nyqd_c(:)
        complex(sp),    allocatable, target :: plp_c(:,:,:)
        real(sp),       allocatable, target :: plp_ct(:,:,:)
        integer(c_int) :: lb_c(3), ub_c(3)
        integer :: hlo, hhi, klo, khi, ierr, ho, ko
        if( gs_nvox <= 0 ) THROW_HARD('flex_gpu_snr_begin_f not called; flex_gpu_snr_batch_f')
        hlo = ceil_div (frl(1,1), OSMPL_PAD_FAC); hhi = floor_div(frl(1,2), OSMPL_PAD_FAC)
        klo = ceil_div (frl(2,1), OSMPL_PAD_FAC); khi = floor_div(frl(2,2), OSMPL_PAD_FAC)
        allocate(nyqd_c(nb)); nyqd_c = int(nyqd(1:nb), c_int)
        ! packed pf-multiple repack of the caller's padded planes (see flex_gpu_insert_batch_f);
        ! incoming assumed-shape arrays index from 1, offsets translate the frlims origin
        ho = 1 - frl(1,1); ko = 1 - frl(2,1)
        allocate(plp_c (hlo:hhi, klo:0, nb), plp_ct(hlo:hhi, klo:0, nb))
        plp_c  = pl_c (ho+OSMPL_PAD_FAC*hlo:ho+OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
            &ko+OSMPL_PAD_FAC*klo:ko:OSMPL_PAD_FAC, 1:nb)
        plp_ct = pl_ct(ho+OSMPL_PAD_FAC*hlo:ho+OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
            &ko+OSMPL_PAD_FAC*klo:ko:OSMPL_PAD_FAC, 1:nb)
        lb_c = int(gs_lb, c_int); ub_c = int(gs_ub, c_int)
        ierr = c_snr_batch(c_loc(plp_c), c_loc(plp_ct), c_loc(rotm), c_loc(nyqd_c), c_loc(vmask), &
            &int(nb,c_int), int(nsym,c_int), int(hlo,c_int), int(hhi,c_int), &
            &int(klo,c_int), int(hlo,c_int), int(hhi,c_int), int(klo,c_int), int(khi,c_int), &
            &int(OSMPL_PAD_FAC,c_int), lb_c, ub_c, real(TINY,c_float), &
            &real(1.0/real(OSMPL_PAD_FAC*OSMPL_PAD_FAC), c_float))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_snr_batch failed')
        deallocate(nyqd_c, plp_c, plp_ct)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_snr_batch_f')
#endif
    end subroutine flex_gpu_snr_batch_f

    !> resident slab insertion over the device-formed adjoint residual pair
    subroutine flex_gpu_snr_batch_res_f( rotm, nyqd, vmask, nb, nsym )
        real,       intent(in), target :: rotm(:,:,:,:)
        integer(1), intent(in), target :: vmask(:)
        integer,    intent(in), target :: nyqd(:)
        integer,    intent(in)         :: nb, nsym
#ifdef USE_FLEX_CUDA
        integer(c_int), allocatable, target :: nyqd_c(:)
        integer(c_int) :: lb_c(3), ub_c(3)
        integer :: ierr
        if( gs_nvox <= 0 ) THROW_HARD('flex_gpu_snr_begin_f not called; flex_gpu_snr_batch_res_f')
        allocate(nyqd_c(nb)); nyqd_c = int(nyqd(1:nb), c_int)
        lb_c = int(gs_lb, c_int); ub_c = int(gs_ub, c_int)
        ierr = c_snr_batch_res(c_loc(rotm), c_loc(nyqd_c), c_loc(vmask), &
            &int(nb,c_int), int(nsym,c_int), lb_c, ub_c, real(TINY,c_float), &
            &real(1.0/real(OSMPL_PAD_FAC*OSMPL_PAD_FAC), c_float))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_snr_batch_resident failed')
        deallocate(nyqd_c)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_snr_batch_res_f')
#endif
    end subroutine flex_gpu_snr_batch_res_f

    !> fetch var/dens and ADD into the host accumulators, then release
    subroutine flex_gpu_snr_end_f( var_acc, dens_acc )
        real, intent(inout) :: var_acc(:,:,:), dens_acc(:,:,:)
#ifdef USE_FLEX_CUDA
        real(sp), allocatable, target :: vbuf(:,:,:), dbuf(:,:,:)
        integer :: ierr
        if( gs_nvox <= 0 ) return
        allocate(vbuf(gs_ub(1)-gs_lb(1)+1, gs_ub(2)-gs_lb(2)+1, gs_ub(3)-gs_lb(3)+1))
        allocate(dbuf(gs_ub(1)-gs_lb(1)+1, gs_ub(2)-gs_lb(2)+1, gs_ub(3)-gs_lb(3)+1))
        ierr = c_snr_end(c_loc(vbuf), c_loc(dbuf))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_snr_end failed')
        var_acc  = var_acc  + vbuf
        dens_acc = dens_acc + dbuf
        deallocate(vbuf, dbuf)
        gs_nvox = 0
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_snr_end_f')
#endif
    end subroutine flex_gpu_snr_end_f

    !> polar reduced-solve device pipeline (P2+P3): resident particle cache, per-chunk cuBLAS,
    !! Vg device-only with the d^4 update as one resident FP64-emulated GEMM
    subroutine flex_gpu_solve_begin_f( np, nsamp2, nk, dt, npk, maxdir, maxm, l_unitc, alo, ahi )
        integer, intent(in) :: np, nsamp2, nk, dt, npk, maxdir, maxm
        logical, intent(in) :: l_unitc
        real,    intent(in) :: alo, ahi
#ifdef USE_FLEX_CUDA
        integer :: ierr
        ierr = c_solve_begin(int(np,c_long_long), int(nsamp2,c_int), int(nk,c_int), &
            &int(dt,c_int), int(npk,c_int), int(maxdir,c_int), int(maxm,c_long_long), &
            &merge(1_c_int,0_c_int,l_unitc), real(alo,c_float), real(ahi,c_float))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_solve_begin failed')
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_solve_begin_f')
#endif
    end subroutine flex_gpu_solve_begin_f

    subroutine flex_gpu_solve_cache_f( xws, wrs )
        real, intent(in), target :: xws(:,:), wrs(:,:)
#ifdef USE_FLEX_CUDA
        integer :: ierr
        ierr = c_solve_cache(c_loc(xws), c_loc(wrs))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_solve_cache failed')
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_solve_cache_f')
#endif
    end subroutine flex_gpu_solve_cache_f

    subroutine flex_gpu_solve_chunk_f( Us_c, Cpk_c, Cm0_c, c00_c, dlist_sl, mcnt, ndc, mtot, wrsum )
        real,     intent(in), target :: Us_c(:,:,:)
        real(dp), intent(in), target :: Cpk_c(:,:,:), Cm0_c(:,:,:), c00_c(:,:)
        integer,  intent(in), target :: dlist_sl(:), mcnt(:)
        integer,  intent(in)         :: ndc, mtot
        real(dp), intent(out), target :: wrsum(:,:)
#ifdef USE_FLEX_CUDA
        integer :: ierr
        ierr = c_solve_chunk(c_loc(Us_c), c_loc(Cpk_c), c_loc(Cm0_c), c_loc(c00_c), &
            &c_loc(dlist_sl), c_loc(mcnt), int(ndc,c_int), int(mtot,c_long_long), c_loc(wrsum))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_solve_chunk failed')
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_solve_chunk_f')
#endif
    end subroutine flex_gpu_solve_chunk_f

    !> fetch and ADD the device Sbb (upper) and Mspk (full) into the caller's accumulators
    subroutine flex_gpu_solve_end_f( Sbb, Mspk, dt, npk )
        integer,  intent(in)    :: dt, npk
        real(dp), intent(inout) :: Sbb(dt,dt), Mspk(npk,npk)
#ifdef USE_FLEX_CUDA
        real(dp), allocatable, target :: Sb(:,:), Mp(:,:)
        integer :: ierr
        allocate(Sb(dt,dt), Mp(npk,npk))
        ierr = c_solve_end(c_loc(Sb), c_loc(Mp))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_solve_end failed')
        Sbb  = Sbb  + Sb
        Mspk = Mspk + Mp
        deallocate(Sb, Mp)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_solve_end_f')
#endif
    end subroutine flex_gpu_solve_end_f

    !> A/B the GPU insertion against the CPU batch on a synthetic case (mirrors
    !! test_coupled_batch_accumulation's setup). Skips cleanly when no device/build.
    subroutine test_flex_gpu_insert()
        use simple_ftiter, only: ftiter
        use simple_flex_reconstructor_latent_ops, only: insert_planes_oversamp_multi_scaled_batch
        integer, parameter :: TEST_BOX=16, TEST_NCOMP=3, TEST_NRECORDS=2
        type(reconstructor) :: rec_cpu(TEST_NCOMP), rec_gpu(TEST_NCOMP)
        type(fplane_type)   :: fpls(TEST_NRECORDS)
        type(ori)   :: orientations(TEST_NRECORDS)
        type(sym)   :: se
        type(ftiter) :: fit_pd
        real(dp)    :: dsc(TEST_NCOMP,TEST_NRECORDS), rsc(TEST_NCOMP,TEST_NRECORDS)
        logical     :: valid(TEST_NRECORDS)
        integer     :: lims_pd(3,2), exp_lb(3), exp_ub(3), i, q, h, k, hp, kp
        real        :: cerr, rerr
        if( .not. flex_gpu_available() )then
            write(logfhandle,'(A)') '>>> TEST flex_gpu insert SKIPPED (no USE_FLEX_CUDA build or no device)'
            return
        endif
        call se%new('c1')
        fit_pd  = ftiter([OSMPL_PAD_FAC*TEST_BOX, OSMPL_PAD_FAC*TEST_BOX, 1], 1.)
        lims_pd = fit_pd%loop_lims(3)
        exp_lb  = [-TEST_BOX/2-2, -TEST_BOX/2-2, -TEST_BOX/2-2]
        exp_ub  = [ TEST_BOX/2+2,  TEST_BOX/2+2,  TEST_BOX/2+2]
        do q = 1, TEST_NCOMP
            call rec_cpu(q)%new([TEST_BOX,TEST_BOX,TEST_BOX], 1.)
            call rec_gpu(q)%new([TEST_BOX,TEST_BOX,TEST_BOX], 1.)
            allocate(rec_cpu(q)%cmat_exp(exp_lb(1):exp_ub(1),exp_lb(2):exp_ub(2),exp_lb(3):exp_ub(3)), &
                &source=CMPLX_ZERO)
            allocate(rec_gpu(q)%cmat_exp(exp_lb(1):exp_ub(1),exp_lb(2):exp_ub(2),exp_lb(3):exp_ub(3)), &
                &source=CMPLX_ZERO)
            allocate(rec_cpu(q)%rho_exp(exp_lb(1):exp_ub(1),exp_lb(2):exp_ub(2),exp_lb(3):exp_ub(3)), &
                &source=0.)
            allocate(rec_gpu(q)%rho_exp(exp_lb(1):exp_ub(1),exp_lb(2):exp_ub(2),exp_lb(3):exp_ub(3)), &
                &source=0.)
        end do
        do i = 1, TEST_NRECORDS
            call orientations(i)%new(.true.)
            call orientations(i)%set_euler([17.*real(i), 31.*real(i), 11.*real(i)])
            allocate(fpls(i)%cmplx_plane(lims_pd(1,1):lims_pd(1,2), lims_pd(2,1):0), source=CMPLX_ZERO)
            allocate(fpls(i)%ctfsq_plane(lims_pd(1,1):lims_pd(1,2), lims_pd(2,1):0), source=0.)
            fpls(i)%frlims = lims_pd
            fpls(i)%nyq    = fit_pd%get_lfny(1)
            do k = -TEST_BOX/2, 0
                kp = OSMPL_PAD_FAC*k
                do h = -TEST_BOX/2, TEST_BOX/2
                    if( h*h + k*k > (TEST_BOX/2)*(TEST_BOX/2+1) ) cycle
                    hp = OSMPL_PAD_FAC*h
                    fpls(i)%cmplx_plane(hp,kp) = cmplx(0.01*real(3*h-2*k+i), 0.02*real(h+k-i))
                    fpls(i)%ctfsq_plane(hp,kp) = 0.2 + 0.03*real(mod(abs(2*h-k+i),9))
                end do
            end do
        end do
        ! zeros exercise the live-component compaction on both paths
        dsc(:,1) = [0.7d0,  0.0d0, 0.3d0];  rsc(:,1) = [0.9d0, 0.0d0, 0.4d0]
        dsc(:,2) = [0.0d0, -1.2d0, 0.5d0];  rsc(:,2) = [0.0d0, 1.1d0, 0.6d0]
        valid    = .true.
        call insert_planes_oversamp_multi_scaled_batch(rec_cpu, se, orientations, fpls, dsc, rsc, &
            &valid, TEST_NRECORDS)
        call flex_gpu_insert_begin_f(rec_gpu, TEST_NCOMP)
        call flex_gpu_insert_batch_f(se, orientations, fpls, dsc, rsc, valid, TEST_NRECORDS)
        call flex_gpu_insert_end_f(rec_gpu)
        cerr = 0.; rerr = 0.
        do q = 1, TEST_NCOMP
            cerr = max(cerr, maxval(abs(rec_cpu(q)%cmat_exp - rec_gpu(q)%cmat_exp)))
            rerr = max(rerr, maxval(abs(rec_cpu(q)%rho_exp  - rec_gpu(q)%rho_exp)))
        end do
        write(logfhandle,'(A,ES10.2,A,ES10.2)') '>>> TEST flex_gpu insert: max|dcmat|=', cerr, &
            &'  max|drho|=', rerr
        if( cerr > 2.e-5 .or. rerr > 2.e-5 ) THROW_HARD('GPU insertion differs from CPU batch path')
        write(logfhandle,'(A)') '>>>   PASSED (GPU N-accumulator insertion matches CPU batch)'
        do i = 1, TEST_NRECORDS
            if( allocated(fpls(i)%cmplx_plane) ) deallocate(fpls(i)%cmplx_plane)
            if( allocated(fpls(i)%ctfsq_plane) ) deallocate(fpls(i)%ctfsq_plane)
            call orientations(i)%kill
        end do
        do q = 1, TEST_NCOMP
            call rec_cpu(q)%kill
            call rec_gpu(q)%kill
        end do
        call se%kill
    end subroutine test_flex_gpu_insert

    !> A/B the GPU coupled-stats kernel (probe M-step shape) against the CPU coupled batch
    subroutine test_flex_gpu_coupled()
        use simple_ftiter, only: ftiter
        use simple_flex_reconstructor_latent_ops, only: insert_planes_oversamp_coupled_batch_scaled
        integer, parameter :: TEST_BOX=16, TEST_NCOMP=3, TEST_NRECORDS=2
        integer, parameter :: NPAIRS = (TEST_NCOMP*(TEST_NCOMP+1))/2
        type(reconstructor) :: rec_cpu(TEST_NCOMP), rec_gpu(TEST_NCOMP)
        type(fplane_type)   :: fpls(TEST_NRECORDS)
        type(ori)    :: orientations(TEST_NRECORDS)
        type(sym)    :: se
        type(ftiter) :: fit_pd
        real(dp) :: z(TEST_NCOMP,TEST_NRECORDS), zz(TEST_NCOMP,TEST_NCOMP,TEST_NRECORDS)
        real(sp), allocatable :: rho_cpu(:,:,:,:), rho_gpu(:,:,:,:)
        logical  :: valid(TEST_NRECORDS)
        integer  :: lims_pd(3,2), exp_lb(3), exp_ub(3), exp_shape(3), i, q, r, h, k, hp, kp
        real     :: cerr, rerr
        if( .not. flex_gpu_available() )then
            write(logfhandle,'(A)') '>>> TEST flex_gpu coupled SKIPPED (no USE_FLEX_CUDA build or no device)'
            return
        endif
        call se%new('c1')
        fit_pd  = ftiter([OSMPL_PAD_FAC*TEST_BOX, OSMPL_PAD_FAC*TEST_BOX, 1], 1.)
        lims_pd = fit_pd%loop_lims(3)
        exp_lb  = [-TEST_BOX/2-2, -TEST_BOX/2-2, -TEST_BOX/2-2]
        exp_ub  = [ TEST_BOX/2+2,  TEST_BOX/2+2,  TEST_BOX/2+2]
        exp_shape = exp_ub - exp_lb + 1
        do q = 1, TEST_NCOMP
            call rec_cpu(q)%new([TEST_BOX,TEST_BOX,TEST_BOX], 1.)
            call rec_gpu(q)%new([TEST_BOX,TEST_BOX,TEST_BOX], 1.)
            allocate(rec_cpu(q)%cmat_exp(exp_lb(1):exp_ub(1),exp_lb(2):exp_ub(2),exp_lb(3):exp_ub(3)), &
                &source=CMPLX_ZERO)
            allocate(rec_gpu(q)%cmat_exp(exp_lb(1):exp_ub(1),exp_lb(2):exp_ub(2),exp_lb(3):exp_ub(3)), &
                &source=CMPLX_ZERO)
        end do
        allocate(rho_cpu(NPAIRS,exp_shape(1),exp_shape(2),exp_shape(3)), source=0.)
        allocate(rho_gpu(NPAIRS,exp_shape(1),exp_shape(2),exp_shape(3)), source=0.)
        do i = 1, TEST_NRECORDS
            call orientations(i)%new(.true.)
            call orientations(i)%set_euler([17.*real(i), 31.*real(i), 11.*real(i)])
            allocate(fpls(i)%cmplx_plane(lims_pd(1,1):lims_pd(1,2), lims_pd(2,1):0), source=CMPLX_ZERO)
            allocate(fpls(i)%ctfsq_plane(lims_pd(1,1):lims_pd(1,2), lims_pd(2,1):0), source=0.)
            allocate(fpls(i)%transfer_plane(lims_pd(1,1):lims_pd(1,2), lims_pd(2,1):0), source=CMPLX_ZERO)
            fpls(i)%frlims = lims_pd
            fpls(i)%nyq    = fit_pd%get_lfny(1)
            do k = -TEST_BOX/2, 0
                kp = OSMPL_PAD_FAC*k
                do h = -TEST_BOX/2, TEST_BOX/2
                    if( h*h + k*k > (TEST_BOX/2)*(TEST_BOX/2+1) ) cycle
                    hp = OSMPL_PAD_FAC*h
                    fpls(i)%cmplx_plane(hp,kp)    = cmplx(0.01*real(3*h-2*k+i), 0.02*real(h+k-i))
                    fpls(i)%transfer_plane(hp,kp) = cmplx(0.7+0.01*real(mod(abs(h+2*k+i),7)), 0.)
                    fpls(i)%ctfsq_plane(hp,kp)    = 0.2 + 0.03*real(mod(abs(2*h-k+i),9))
                end do
            end do
        end do
        z(:,1) = [0.4d0, -0.2d0, 0.7d0]
        z(:,2) = [-0.1d0, 0.8d0, 0.3d0]
        do i = 1, TEST_NRECORDS
            do r = 1, TEST_NCOMP
                do q = 1, TEST_NCOMP
                    zz(q,r,i) = z(q,i)*z(r,i)
                end do
            end do
        end do
        valid = .true.
        call insert_planes_oversamp_coupled_batch_scaled(rec_cpu, rho_cpu, se, orientations, fpls, &
            &z, zz, valid, TEST_NRECORDS)
        call flex_gpu_coupled_begin_f(rec_gpu, TEST_NCOMP, NPAIRS)
        call flex_gpu_coupled_batch_f(se, orientations, fpls, z, zz, valid, TEST_NRECORDS)
        call flex_gpu_coupled_end_f(rec_gpu, rho_gpu)
        cerr = 0.
        do q = 1, TEST_NCOMP
            cerr = max(cerr, maxval(abs(rec_cpu(q)%cmat_exp - rec_gpu(q)%cmat_exp)))
        end do
        rerr = maxval(abs(rho_cpu - rho_gpu))
        write(logfhandle,'(A,ES10.2,A,ES10.2)') '>>> TEST flex_gpu coupled: max|dcmat|=', cerr, &
            &'  max|drho|=', rerr
        if( cerr > 2.e-5 .or. rerr > 2.e-5 ) THROW_HARD('GPU coupled insertion differs from CPU batch path')
        write(logfhandle,'(A)') '>>>   PASSED (GPU coupled-stats insertion matches CPU batch)'
        do i = 1, TEST_NRECORDS
            if( allocated(fpls(i)%cmplx_plane) )    deallocate(fpls(i)%cmplx_plane)
            if( allocated(fpls(i)%ctfsq_plane) )    deallocate(fpls(i)%ctfsq_plane)
            if( allocated(fpls(i)%transfer_plane) ) deallocate(fpls(i)%transfer_plane)
            call orientations(i)%kill
        end do
        do q = 1, TEST_NCOMP
            call rec_cpu(q)%kill
            call rec_gpu(q)%kill
        end do
        call se%kill
        deallocate(rho_cpu, rho_gpu)
    end subroutine test_flex_gpu_coupled

    !> upload the polar grid tables for the device sampler (once per stage)
    subroutine flex_gpu_psample_begin_f( rad, cs, sn, sqwq, rbeg, rend, nsamp, nk, &
        &nrad, ncs, nsn, nwq, nsampn )
        integer, intent(in) :: nsamp, nk, nsampn
        real,    intent(in), target :: rad(nsamp), cs(nsamp), sn(nsamp), sqwq(nsamp)
        integer, intent(in) :: rbeg(nk), rend(nk)
        real,    intent(in), target :: nrad(*), ncs(*), nsn(*), nwq(*)
#ifdef USE_FLEX_CUDA
        integer(c_int), allocatable :: ring0(:), nang(:)
        integer :: ir, j, ierr
        allocate(ring0(nsamp), nang(nk))
        do ir = 1, nk
            nang(ir) = int(rend(ir) - rbeg(ir) + 1, c_int)
            do j = rbeg(ir), rend(ir)
                ring0(j) = int(ir - 1, c_int)
            end do
        end do
        ierr = c_psample_begin(c_loc(rad), c_loc(cs), c_loc(sn), c_loc(sqwq), ring0, nang, &
            &int(nsamp,c_int), int(nk,c_int), c_loc(nrad), c_loc(ncs), c_loc(nsn), c_loc(nwq), &
            &int(nsampn,c_int))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_psample_begin failed')
        deallocate(ring0, nang)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_psample_begin_f')
#endif
    end subroutine flex_gpu_psample_begin_f

    !> one batch through the device polar sampler. Outputs land DIRECTLY in the caller's
    !! full-size arrays at column col0 (the batch occupies contiguous columns); pwv/tazv are
    !! per-record (pw = sum nwq*|y|^2 on the noise rings; taz = sum over rings of
    !! sqrt(var)/mean of |T|^2 -- divide by nk and reduce on the host).
    subroutine flex_gpu_psample_batch_f( fpls, cav, sav, valid, nrecords, want_halves, &
        &xwall, xw1all, xw2all, wrall, col0, pwv, tazv )
        use simple_math, only: ceil_div, floor_div
        type(fplane_type), intent(in)    :: fpls(:)
        real,              intent(in)    :: cav(:), sav(:)
        logical,           intent(in)    :: valid(:)
        integer,           intent(in)    :: nrecords, col0
        logical,           intent(in)    :: want_halves
        real,              intent(inout), target :: xwall(:,:), xw1all(:,:), xw2all(:,:), wrall(:,:)
        real,              intent(out),   target :: pwv(:), tazv(:)
#ifdef USE_FLEX_CUDA
        complex(sp), allocatable, target :: ply(:,:,:), plt(:,:,:)
        real(sp),    allocatable, target :: ca_c(:), sa_c(:)
        integer(c_int8_t), allocatable, target :: vmask(:)
        integer :: frl(3,2), i, ierr, iref
        integer :: hlo, hhi, klo
        iref = 0
        do i = 1, nrecords
            if( valid(i) )then
                if( .not. allocated(fpls(i)%transfer_plane) ) &
                    &THROW_HARD('transfer plane missing; flex_gpu_psample_batch_f')
                if( iref == 0 ) iref = i
            endif
        end do
        pwv = 0.; tazv = 0.
        if( iref == 0 ) return
        frl = fpls(iref)%frlims
        hlo = ceil_div (frl(1,1), OSMPL_PAD_FAC); hhi = floor_div(frl(1,2), OSMPL_PAD_FAC)
        klo = ceil_div (frl(2,1), OSMPL_PAD_FAC)
        allocate(ply(hlo:hhi, klo:0, nrecords), plt(hlo:hhi, klo:0, nrecords), source=cmplx(0.,0.))
        allocate(ca_c(nrecords), source=1._sp)
        allocate(sa_c(nrecords), source=0._sp)
        allocate(vmask(nrecords), source=0_c_int8_t)
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            vmask(i) = 1_c_int8_t
            ply(:,:,i) = fpls(i)%cmplx_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)
            plt(:,:,i) = fpls(i)%transfer_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)
            ca_c(i) = cav(i)
            sa_c(i) = sav(i)
        end do
        !$omp end parallel do
        ierr = c_psample_batch(c_loc(ply), c_loc(plt), c_loc(ca_c), c_loc(sa_c), c_loc(vmask), &
            &int(nrecords,c_int), merge(1_c_int, 0_c_int, want_halves), &
            &int(hlo,c_int), int(hhi,c_int), int(klo,c_int), &
            &c_loc(xwall(1,col0)), c_loc(xw1all(1,col0)), c_loc(xw2all(1,col0)), &
            &c_loc(wrall(1,col0)), c_loc(pwv), c_loc(tazv))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_psample_batch failed')
        deallocate(ply, plt, ca_c, sa_c, vmask)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_psample_batch_f')
#endif
    end subroutine flex_gpu_psample_batch_f

    !> resident psample: samples straight off the separated device-prep pair with the
    !! registered geometry — no plane packing, no upload
    subroutine flex_gpu_psample_batch_res_f( cav, sav, valid, nrecords, want_halves, &
        &xwall, xw1all, xw2all, wrall, col0, pwv, tazv )
        real,    intent(in)  :: cav(:), sav(:)
        logical, intent(in)  :: valid(:)
        integer, intent(in)  :: nrecords, col0
        logical, intent(in)  :: want_halves
        real, intent(inout), target :: xwall(:,:), xw1all(:,:), xw2all(:,:), wrall(:,:)
        real, intent(out),   target :: pwv(:), tazv(:)
#ifdef USE_FLEX_CUDA
        real(sp),          allocatable, target :: ca_c(:), sa_c(:)
        integer(c_int8_t), allocatable, target :: vmask(:)
        integer :: i, ierr
        pwv = 0.; tazv = 0.
        if( .not. any(valid(:nrecords)) ) return
        allocate(ca_c(nrecords), source=1._sp)
        allocate(sa_c(nrecords), source=0._sp)
        allocate(vmask(nrecords), source=0_c_int8_t)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            vmask(i) = 1_c_int8_t
            ca_c(i)  = cav(i)
            sa_c(i)  = sav(i)
        end do
        ierr = c_psample_batch_res(c_loc(ca_c), c_loc(sa_c), c_loc(vmask), &
            &int(nrecords,c_int), merge(1_c_int, 0_c_int, want_halves), &
            &c_loc(xwall(1,col0)), c_loc(xw1all(1,col0)), c_loc(xw2all(1,col0)), &
            &c_loc(wrall(1,col0)), c_loc(pwv), c_loc(tazv))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_psample_batch_resident failed')
        deallocate(ca_c, sa_c, vmask)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_psample_batch_res_f')
#endif
    end subroutine flex_gpu_psample_batch_res_f

    !> stage-persistent Cartesian noise diagnostics over the resident whitened planes
    subroutine flex_gpu_psample_diag_begin_f( nyq )
        integer, intent(in) :: nyq
#ifdef USE_FLEX_CUDA
        integer :: ierr
        ierr = c_psample_diag_begin(int(nyq,c_int))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_psample_diag_begin failed')
#endif
    end subroutine flex_gpu_psample_diag_begin_f

    subroutine flex_gpu_psample_diag_batch_f( nrecords, sh_lo )
        integer, intent(in) :: nrecords, sh_lo
#ifdef USE_FLEX_CUDA
        integer :: ierr
        ierr = c_psample_diag_batch(int(nrecords,c_int), int(sh_lo,c_int))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_psample_diag_batch failed')
#endif
    end subroutine flex_gpu_psample_diag_batch_f

    !> fetch [pwsh(0:nyq) | cntsh(0:nyq) | hfpw | hfcnt] accumulated across the stage
    subroutine flex_gpu_psample_diag_fetch_f( pwsh, cntsh, hfpw, hfcnt, nyq )
        integer,  intent(in)  :: nyq
        real(dp), intent(out) :: pwsh(0:nyq), cntsh(0:nyq), hfpw, hfcnt
#ifdef USE_FLEX_CUDA
        real(dp), allocatable, target :: acc(:)
        integer :: ierr
        allocate(acc(2*(nyq+1)+2), source=0.d0)
        ierr = c_psample_diag_fetch(c_loc(acc))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_psample_diag_fetch failed')
        pwsh  = acc(1:nyq+1)
        cntsh = acc(nyq+2:2*(nyq+1))
        hfpw  = acc(2*(nyq+1)+1)
        hfcnt = acc(2*(nyq+1)+2)
        deallocate(acc)
#else
        pwsh = 0.d0; cntsh = 0.d0; hfpw = 0.d0; hfcnt = 0.d0
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_psample_diag_fetch_f')
#endif
    end subroutine flex_gpu_psample_diag_fetch_f

    subroutine flex_gpu_psample_free_f
#ifdef USE_FLEX_CUDA
        integer :: ierr
        ierr = c_psample_free()
#endif
    end subroutine flex_gpu_psample_free_f

    !> upload the mean + basis expanded lattices for the fused E-step (per probe iteration)
    subroutine flex_gpu_estep_vols_f( mean_rec, basis_recs, ncomp )
        type(reconstructor), intent(in) :: mean_rec
        type(reconstructor), intent(in) :: basis_recs(ncomp)
        integer,             intent(in) :: ncomp
#ifdef USE_FLEX_CUDA
        complex(sp), allocatable, target :: vbuf(:,:,:)
        integer(c_int) :: lb_c(3), ub_c(3)
        integer :: q, ierr
        lb_c = int(lbound(mean_rec%cmat_exp), c_int)
        ub_c = int(ubound(mean_rec%cmat_exp), c_int)
        allocate(vbuf, mold=mean_rec%cmat_exp)
        vbuf = mean_rec%cmat_exp
        ierr = c_estep_vols(c_loc(vbuf), 1_c_int, int(ncomp+1,c_int), lb_c, ub_c)
        if( ierr /= 0 ) THROW_HARD('flex_gpu_estep_vols failed (mean)')
        do q = 1, ncomp
            vbuf = basis_recs(q)%cmat_exp
            ierr = c_estep_vols(c_loc(vbuf), int(q+1,c_int), &
                &int(ncomp+1,c_int), lb_c, ub_c)
            if( ierr /= 0 ) THROW_HARD('flex_gpu_estep_vols failed (basis)')
        end do
        deallocate(vbuf)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_estep_vols_f')
#endif
    end subroutine flex_gpu_estep_vols_f

    !> one E-step batch: packed (conj(T)y, ctfsq) up, per-record stat vector back.
    !! stats layout per record: [e_mm, myv, b(1..nc), c(1..nc), G upper q<=r]
    subroutine flex_gpu_estep_batch_f( fpls, orientations, valid, nrecords, stats )
        use simple_math, only: ceil_div, floor_div
        type(fplane_type), intent(in)    :: fpls(:)
        type(ori),         intent(inout) :: orientations(:)
        logical,           intent(in)    :: valid(:)
        integer,           intent(in)    :: nrecords
        real(dp),          intent(out), target, contiguous :: stats(:,:)
#ifdef USE_FLEX_CUDA
        complex(sp), allocatable, target :: pl_c(:,:,:)
        real(sp),    allocatable, target :: pl_ct(:,:,:), rotm(:,:,:)
        integer(c_int8_t), allocatable, target :: vmask(:)
        integer :: frl(3,2), i, ierr, iref, nyq_eff
        integer :: hlo, hhi, klo
        iref = 0
        do i = 1, nrecords
            if( valid(i) )then
                if( iref == 0 ) iref = i
            endif
        end do
        stats = 0.d0
        if( iref == 0 ) return
        frl = fpls(iref)%frlims
        hlo = ceil_div (frl(1,1), OSMPL_PAD_FAC); hhi = floor_div(frl(1,2), OSMPL_PAD_FAC)
        klo = ceil_div (frl(2,1), OSMPL_PAD_FAC)
        allocate(pl_c (hlo:hhi, klo:0, nrecords), source=cmplx(0.,0.))
        allocate(pl_ct(hlo:hhi, klo:0, nrecords), source=0.)
        allocate(rotm(3,3,nrecords), source=0.)
        allocate(vmask(nrecords), source=0_c_int8_t)
        nyq_eff = g2_lfny
        if( fpls(iref)%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpls(iref)%nyq / OSMPL_PAD_FAC))
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            vmask(i) = 1_c_int8_t
            pl_c (:,:,i) = conjg(fpls(i)%transfer_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)) * &
                &fpls(i)%cmplx_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)
            pl_ct(:,:,i) = fpls(i)%ctfsq_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)
            rotm(:,:,i) = orientations(i)%get_mat()
        end do
        !$omp end parallel do
        ierr = c_estep_batch(c_loc(pl_c), c_loc(pl_ct), c_loc(rotm), c_loc(vmask), &
            &int(nrecords,c_int), int(hlo,c_int), int(hhi,c_int), int(klo,c_int), 0_c_int, &
            &int(nyq_eff*(nyq_eff+1),c_int), c_loc(stats))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_estep_batch failed')
        deallocate(pl_c, pl_ct, rotm, vmask)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_estep_batch_f')
#endif
    end subroutine flex_gpu_estep_batch_f

    !> device residual conj(T)y - a*ctfsq*proj0, fetched packed and unpacked into the fplanes
    !! (the M-step then consumes them with premultiplied=.true.)
    subroutine flex_gpu_estep_resid_f( fpls, avec, valid, nrecords, fetch )
        use simple_math, only: ceil_div, floor_div
        type(fplane_type), intent(inout) :: fpls(:)
        real,              intent(in), target, contiguous :: avec(:)
        logical,           intent(in)    :: valid(:)
        integer,           intent(in)    :: nrecords
        logical, optional, intent(in)    :: fetch
#ifdef USE_FLEX_CUDA
        complex(sp), allocatable, target :: resid(:,:,:)
        integer(c_int8_t), allocatable, target :: vmask(:)
        integer :: frl(3,2), i, ierr, iref
        integer :: hlo, hhi, klo
        logical :: l_fetch
        iref = 0
        do i = 1, nrecords
            if( valid(i) .and. iref == 0 ) iref = i
        end do
        if( iref == 0 ) return
        frl = fpls(iref)%frlims
        hlo = ceil_div (frl(1,1), OSMPL_PAD_FAC); hhi = floor_div(frl(1,2), OSMPL_PAD_FAC)
        klo = ceil_div (frl(2,1), OSMPL_PAD_FAC)
        l_fetch = .true.
        if( present(fetch) ) l_fetch = fetch
        allocate(vmask(nrecords), source=0_c_int8_t)
        do i = 1, nrecords
            if( valid(i) ) vmask(i) = 1_c_int8_t
        end do
        if( l_fetch )then
            allocate(resid(hlo:hhi, klo:0, nrecords), source=cmplx(0.,0.))
            ierr = c_estep_resid(c_loc(avec), c_loc(vmask), int(nrecords,c_int), &
                &int(hlo,c_int), int(hhi,c_int), int(klo,c_int), 0_c_int, c_loc(resid))
            if( ierr /= 0 ) THROW_HARD('flex_gpu_estep_resid failed')
            !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
            do i = 1, nrecords
                if( .not. valid(i) ) cycle
                fpls(i)%cmplx_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                    &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC) = resid(:,:,i)
            end do
            !$omp end parallel do
            deallocate(resid)
        else
            ! residual stays resident: the M-step consumes it via the _res entry. Sentinel
            ! bounds (1 > 0) make the device use its REGISTERED geometry -- fpls may be
            ! unbuilt on this path, and reading their empty frlims silently disabled the
            ! mean deflation (the accidental no-deflation configuration; see the deflation A/B)
            ierr = c_estep_resid(c_loc(avec), c_loc(vmask), int(nrecords,c_int), &
                &1_c_int, 0_c_int, 0_c_int, 0_c_int, c_null_ptr)
            if( ierr /= 0 ) THROW_HARD('flex_gpu_estep_resid failed')
        endif
        deallocate(vmask)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_estep_resid_f')
#endif
    end subroutine flex_gpu_estep_resid_f

    !> banked M-step over the RESIDENT residual (fused E-step path): no plane re-upload
    subroutine flex_gpu_coupled_batch_banked_res_f( dscales, rscales, valid, nrecords, &
        &dirid, cav, sav )
        real(dp), intent(in) :: dscales(:,:), rscales(:,:)
        logical,  intent(in) :: valid(:)
        integer,  intent(in) :: nrecords
        integer,  intent(in) :: dirid(:)
        real,     intent(in) :: cav(:), sav(:)
#ifdef USE_FLEX_CUDA
        real(sp), allocatable, target :: dsc(:,:), rsc(:,:), ca_c(:), sa_c(:)
        integer(c_int),    allocatable, target :: did(:)
        integer(c_int8_t), allocatable, target :: vmask(:)
        integer(c_int) :: lb_c(3), ub_c(3)
        integer :: i, ierr
        if( g2_ncomp <= 0 ) THROW_HARD('flex_gpu_coupled_begin_f not called; banked_res')
        allocate(dsc(g2_ncomp,nrecords), rsc(g2_nrho,nrecords), source=0.)
        allocate(ca_c(nrecords), source=1._sp)
        allocate(sa_c(nrecords), source=0._sp)
        allocate(did(nrecords), source=1_c_int)
        allocate(vmask(nrecords), source=0_c_int8_t)
        do i = 1, nrecords
            did(i) = int(max(1, dirid(i)), c_int)
            if( .not. valid(i) ) cycle
            vmask(i) = 1_c_int8_t
            dsc(:,i) = real(dscales(1:g2_ncomp,i))
            rsc(:,i) = real(rscales(1:g2_nrho,i))
            ca_c(i)  = cav(i)
            sa_c(i)  = sav(i)
        end do
        lb_c = int(g2_lb, c_int); ub_c = int(g2_ub, c_int)
        ierr = c_coupled_batch_banked_res(c_loc(ca_c), c_loc(sa_c), c_loc(dsc), c_loc(rsc), &
            &c_loc(vmask), did, int(nrecords,c_int), int(g2_ncomp,c_int), int(g2_nrho,c_int), &
            &lb_c, ub_c)
        if( ierr /= 0 ) THROW_HARD('flex_gpu_coupled_batch_banked_res failed')
        deallocate(dsc, rsc, ca_c, sa_c, did, vmask)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_coupled_batch_banked_res_f')
#endif
    end subroutine flex_gpu_coupled_batch_banked_res_f

    subroutine flex_gpu_estep_free_f
#ifdef USE_FLEX_CUDA
        integer :: ierr
        ierr = c_estep_free()
#endif
    end subroutine flex_gpu_estep_free_f

    !> POLAR E-STEP device port (stage 2): per-stage ring geometry upload. Same grid-array
    !! packing convention as the psample wrapper (0-based ring index per sample, per-ring
    !! angle counts from rbeg/rend).
    subroutine flex_gpu_poles_begin_f( rad, cs, sn, sqwq, rbeg, rend, nsamp, nk )
        integer, intent(in) :: nsamp, nk
        real,    intent(in), target :: rad(nsamp), cs(nsamp), sn(nsamp), sqwq(nsamp)
        integer, intent(in) :: rbeg(nk), rend(nk)
#ifdef USE_FLEX_CUDA
        integer(c_int), allocatable :: ring0(:), nang(:)
        integer :: ir, j, ierr
        allocate(ring0(nsamp), nang(nk))
        do ir = 1, nk
            nang(ir) = int(rend(ir) - rbeg(ir) + 1, c_int)
            do j = rbeg(ir), rend(ir)
                ring0(j) = int(ir - 1, c_int)
            end do
        end do
        ierr = c_poles_begin(c_loc(rad), c_loc(cs), c_loc(sn), c_loc(sqwq), ring0, nang, &
            &int(nsamp,c_int), int(nk,c_int))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_poles_begin failed')
        deallocate(ring0, nang)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_poles_begin_f')
#endif
    end subroutine flex_gpu_poles_begin_f

    !> once per EM iteration: the host-built shared-direction bank and ring Gram tables go up
    !! (UsallE(nsamp2,0:ncomp,ndir) float, CfE(ncomp*ncomp,nk,ndir)/Cm0E(ncomp,nk,ndir)/
    !! c00E(nk,ndir) double, all Fortran column-major -- the device indexes them as such)
    subroutine flex_gpu_poles_bank_f( Us, Cf, Cm0, c00, ncomp, ndir )
        real(sp), intent(in), target, contiguous :: Us(:,:,:)
        real(dp), intent(in), target, contiguous :: Cf(:,:,:), Cm0(:,:,:), c00(:,:)
        integer,  intent(in) :: ncomp, ndir
#ifdef USE_FLEX_CUDA
        integer :: ierr
        ierr = c_poles_bank(c_loc(Us), c_loc(Cf), c_loc(Cm0), c_loc(c00), &
            &int(ncomp,c_int), int(ndir,c_int))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_poles_bank failed')
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_poles_bank_f')
#endif
    end subroutine flex_gpu_poles_bank_f

    !> one polar E-step batch: RAW (y, T) plane sections up (the device derives the packed
    !! conj(T)y / |T|^2 payload itself so the ring sampler can keep the product-of-interpolants
    !! convention), per-record stat vectors back in the fused-E-step layout
    !! [e_mm, myv, b(1..nc), c(1..nc), G upper]. nyq2_ex = rhyb*(rhyb+1) or 0 for pure rings.
    subroutine flex_gpu_poles_batch_f( fpls, orientations, cav, sav, dirid, valid, nrecords, &
            &nyq2_ex, stats )
        use simple_math, only: ceil_div, floor_div
        type(fplane_type), intent(in)    :: fpls(:)
        type(ori),         intent(inout) :: orientations(:)
        real,              intent(in)    :: cav(:), sav(:)
        integer,           intent(in)    :: dirid(:)
        logical,           intent(in)    :: valid(:)
        integer,           intent(in)    :: nrecords, nyq2_ex
        real(dp),          intent(out), target, contiguous :: stats(:,:)
#ifdef USE_FLEX_CUDA
        complex(sp), allocatable, target :: pl_y(:,:,:), pl_t(:,:,:)
        real(sp),    allocatable, target :: rotm(:,:,:), ca_c(:), sa_c(:)
        integer(c_int),    allocatable :: did(:)
        integer(c_int8_t), allocatable, target :: vmask(:)
        integer :: frl(3,2), i, ierr, iref
        integer :: hlo, hhi, klo
        stats = 0.d0
        iref = 0
        do i = 1, nrecords
            if( valid(i) )then
                if( .not. allocated(fpls(i)%transfer_plane) ) &
                    &THROW_HARD('forward transfer plane does not exist; flex_gpu_poles_batch_f')
                if( iref == 0 ) iref = i
            endif
        end do
        if( iref == 0 ) return
        frl = fpls(iref)%frlims
        do i = 1, nrecords
            if( valid(i) .and. any(fpls(i)%frlims /= frl) ) &
                &THROW_HARD('heterogeneous frlims in batch; flex_gpu_poles_batch_f')
        end do
        hlo = ceil_div (frl(1,1), OSMPL_PAD_FAC); hhi = floor_div(frl(1,2), OSMPL_PAD_FAC)
        klo = ceil_div (frl(2,1), OSMPL_PAD_FAC)
        allocate(pl_y(hlo:hhi, klo:0, nrecords), pl_t(hlo:hhi, klo:0, nrecords), &
            &source=cmplx(0.,0.))
        allocate(rotm(3,3,nrecords), source=0.)
        allocate(ca_c(nrecords), source=1._sp)
        allocate(sa_c(nrecords), source=0._sp)
        allocate(did(nrecords), source=1_c_int)
        allocate(vmask(nrecords), source=0_c_int8_t)
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, nrecords
            did(i) = int(max(1, dirid(i)), c_int)
            if( .not. valid(i) ) cycle
            vmask(i) = 1_c_int8_t
            ! packed pf-multiple sections (see flex_gpu_insert_batch_f), RAW pair
            pl_y(:,:,i) = fpls(i)%cmplx_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)
            pl_t(:,:,i) = fpls(i)%transfer_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)
            rotm(:,:,i) = orientations(i)%get_mat()
            ca_c(i) = cav(i)
            sa_c(i) = sav(i)
        end do
        !$omp end parallel do
        ierr = c_poles_batch(c_loc(pl_y), c_loc(pl_t), c_loc(rotm), c_loc(ca_c), c_loc(sa_c), &
            &did, c_loc(vmask), int(nrecords,c_int), int(hlo,c_int), int(hhi,c_int), &
            &int(klo,c_int), int(nyq2_ex,c_int), c_loc(stats))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_poles_batch failed')
        deallocate(pl_y, pl_t, rotm, ca_c, sa_c, did, vmask)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_poles_batch_f')
#endif
    end subroutine flex_gpu_poles_batch_f

    subroutine flex_gpu_poles_free_f
#ifdef USE_FLEX_CUDA
        integer :: ierr
        ierr = c_poles_free()
#endif
    end subroutine flex_gpu_poles_free_f

    !> device prep setup: geometry + resolution mask + soft-mask coordinates, once per stage
    subroutine flex_gpu_prep_begin_f( lmsk, n1, n1o, maxrec, mskrad, l_ctf )
        integer, intent(in) :: n1, n1o, maxrec
        logical, intent(in) :: lmsk(:,:,:)
        real,    intent(in) :: mskrad
        logical, intent(in) :: l_ctf
#ifdef USE_FLEX_CUDA
        integer(c_int8_t), allocatable, target :: lm(:,:)
        real(sp),          allocatable, target :: cs2(:)
        real    :: minlen, c
        integer :: i, j, wbg, ierr
        allocate(lm(n1,n1), source=0_c_int8_t)
        do j = 1, n1
            do i = 1, n1
                if( lmsk(i,j,1) ) lm(i,j) = 1_c_int8_t
            end do
        end do
        allocate(cs2(n1))
        do i = 1, n1
            c      = -0.5*real(n1) + real(i-1)
            cs2(i) = c*c
        end do
        if( mskrad > 0. )then
            minlen = real(min(nint(2.0*(mskrad + COSMSKHALFWIDTH)), n1))
        else
            ! taper path: the minlen slot carries the taper width wtap
            minlen = real(min(nint(COSMSKHALFWIDTH), n1/2))
        endif
        wbg    = max(1, nint(COSMSKHALFWIDTH)/2)
        ierr = c_prep_begin(c_loc(lm), c_loc(cs2), int(n1,c_int), int(n1o,c_int), &
            &int(maxrec,c_int), real(mskrad,c_float), real(minlen,c_float), int(wbg,c_int), &
            &merge(1_c_int, 0_c_int, l_ctf))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_prep_begin failed')
        l_prep_dev_ready = .true.
        deallocate(lm, cs2)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_prep_begin_f')
#endif
    end subroutine flex_gpu_prep_begin_f

    !> one prep batch: raw images + per-record CTF scalars and shift phases up; the packed
    !! planes are left RESIDENT for estep_batch_res / the banked M-step. frlims/nyq_pd give
    !! the packed geometry (from the same ftiter convention gen_fplane4rec uses).
    subroutine flex_gpu_prep_batch_f( imgs, ctfparms_arr, shifts, valid, nrecords, n1, &
        &frlims, nyq_pd, sig2_ups )
        use simple_math, only: ceil_div, floor_div
        use simple_ctf,  only: ctf
        class(image),     intent(inout) :: imgs(:)
        type(ctfparams),  intent(in)    :: ctfparms_arr(:)
        real,             intent(in)    :: shifts(:,:)
        logical,          intent(in)    :: valid(:)
        integer,          intent(in)    :: nrecords, n1, frlims(3,2), nyq_pd
        real, optional,   intent(in), target, contiguous :: sig2_ups(:,:)
#ifdef USE_FLEX_CUDA
        integer :: nsig
        type(ctf)     :: tfun
        type(ctfvars) :: cv
        real(sp), allocatable, target :: raw(:,:,:), ctfp(:,:), psh(:,:)
        integer(c_int8_t), allocatable, target :: vmask(:)
        real, pointer :: rmatp(:,:,:)
        real     :: shconst_pd
        integer  :: i, ierr, hlo, hhi, klo, n1o
        hlo = ceil_div (frlims(1,1), OSMPL_PAD_FAC); hhi = floor_div(frlims(1,2), OSMPL_PAD_FAC)
        klo = ceil_div (frlims(2,1), OSMPL_PAD_FAC)
        n1o = 2*(frlims(1,2))    ! padded box from redundant lims: hmax = n1o/2
        allocate(raw(n1,n1,nrecords), source=0.0_sp)
        allocate(ctfp(8,nrecords), psh(2,nrecords), source=0.0_sp)
        allocate(vmask(nrecords), source=0_c_int8_t)
        shconst_pd = PI/real(n1o/2)
        !$omp parallel do default(shared) private(i,rmatp,tfun,cv) schedule(static) proc_bind(close)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            vmask(i) = 1_c_int8_t
            call imgs(i)%get_rmat_ptr(rmatp)
            raw(:,:,i) = rmatp(1:n1,1:n1,1)
            psh(:,i)   = -shifts(:,i) * shconst_pd
            if( ctfparms_arr(i)%ctfflag /= CTFFLAG_NO )then
                tfun = ctf(ctfparms_arr(i)%smpd, ctfparms_arr(i)%kv, ctfparms_arr(i)%cs, &
                    &ctfparms_arr(i)%fraca)
                call tfun%init(ctfparms_arr(i)%dfx, ctfparms_arr(i)%dfy, ctfparms_arr(i)%angast)
                cv = tfun%get_ctfvars(ctfparms_arr(i)%phshift)
                ctfp(1,i) = cv%dfx + cv%dfy
                ctfp(2,i) = cv%dfx - cv%dfy
                ctfp(3,i) = cv%angast
                ctfp(4,i) = cv%phshift
                ctfp(5,i) = cv%amp_contr_const
                ctfp(6,i) = cv%wl
                ctfp(7,i) = 0.5*cv%wl*cv%wl*cv%cs
                ctfp(8,i) = merge(2.0, 1.0, ctfparms_arr(i)%ctfflag == CTFFLAG_FLIP)
            else
                ctfp(8,i) = 0.0
            endif
        end do
        !$omp end parallel do
        nsig = 0
        if( present(sig2_ups) ) nsig = size(sig2_ups,1)
        if( nsig > 0 )then
            ierr = c_prep_batch(c_loc(raw), c_loc(ctfp), c_loc(psh), c_loc(vmask), &
                &c_loc(sig2_ups), int(nsig,c_int), &
                &int(nrecords,c_int), int(hlo,c_int), int(hhi,c_int), int(klo,c_int), &
                &int(nyq_pd,c_int))
        else
            ierr = c_prep_batch(c_loc(raw), c_loc(ctfp), c_loc(psh), c_loc(vmask), &
                &c_loc(raw), 0_c_int, &
                &int(nrecords,c_int), int(hlo,c_int), int(hhi,c_int), int(klo,c_int), &
                &int(nyq_pd,c_int))
        endif
        if( ierr /= 0 ) THROW_HARD('flex_gpu_prep_batch failed')
        deallocate(raw, ctfp, psh, vmask)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_prep_batch_f')
#endif
    end subroutine flex_gpu_prep_batch_f

    !> validation cross-check: compare the device-prepped packed planes against the CPU
    !! prep (conj(T)y and ctfsq) for one batch; prints max deviations
    subroutine flex_gpu_prep_check_f( fpls, valid, nrecords )
        use simple_math, only: ceil_div, floor_div
        type(fplane_type), intent(in) :: fpls(:)
        logical,           intent(in) :: valid(:)
        integer,           intent(in) :: nrecords
#ifdef USE_FLEX_CUDA
        complex(sp), allocatable, target :: plc(:,:,:)
        real(sp),    allocatable, target :: plct(:,:,:)
        complex(sp), allocatable :: ref_c(:,:)
        real(sp),    allocatable :: ref_t(:,:)
        integer :: frl(3,2), i, ierr, iref, hlo, hhi, klo
        real    :: ec, et, em, ej, en, ek, denom
        iref = 0
        do i = 1, nrecords
            if( valid(i) .and. iref == 0 ) iref = i
        end do
        if( iref == 0 ) return
        frl = fpls(iref)%frlims
        hlo = ceil_div (frl(1,1), OSMPL_PAD_FAC); hhi = floor_div(frl(1,2), OSMPL_PAD_FAC)
        klo = ceil_div (frl(2,1), OSMPL_PAD_FAC)
        allocate(plc(hlo:hhi, klo:0, nrecords), source=cmplx(0.,0.))
        allocate(plct(hlo:hhi, klo:0, nrecords), source=0.)
        ierr = c_prep_fetch(c_loc(plc), c_loc(plct), int(nrecords,c_int), int(hlo,c_int), &
            &int(hhi,c_int), int(klo,c_int))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_prep_fetch failed')
        allocate(ref_c(hlo:hhi, klo:0), ref_t(hlo:hhi, klo:0))
        ec = 0.; et = 0.; em = 0.; ej = 0.; en = 0.; ek = 0.; denom = 1.e-12
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            ref_c = conjg(fpls(i)%transfer_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)) * &
                &fpls(i)%cmplx_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)
            ref_t = fpls(i)%ctfsq_plane(OSMPL_PAD_FAC*hlo:OSMPL_PAD_FAC*hhi:OSMPL_PAD_FAC, &
                &OSMPL_PAD_FAC*klo:0:OSMPL_PAD_FAC)
            ec    = max(ec, maxval(abs(ref_c - plc(:,:,i))))
            et    = max(et, maxval(abs(ref_t - plct(:,:,i))))
            em    = max(em, maxval(abs(abs(ref_c) - abs(plc(:,:,i)))))
            ej    = max(ej, maxval(abs(ref_c - conjg(plc(:,:,i)))))
            en    = max(en, maxval(abs(ref_c + plc(:,:,i))))
            ek    = max(ek, maxval(abs(ref_c + conjg(plc(:,:,i)))))
            denom = max(denom, maxval(abs(ref_c)))
        end do
        write(logfhandle,'(A,ES10.2,A,ES10.2,A,ES10.2)') &
            &'>>> FLEX_GPU PREP CHECK: max|d(conjT*y)|=', ec, '  rel=', ec/denom, &
            &'  max|d(ctfsq)|=', et
        write(logfhandle,'(A,ES10.2,A,ES10.2)') &
            &'>>> FLEX_GPU PREP CHECK: max|d MAGNITUDE|=', em, &
            &'  max|d vs CONJUGATE|=', ej
        write(logfhandle,'(A,ES10.2,A,ES10.2)') &
            &'>>> FLEX_GPU PREP CHECK: max|d vs NEGATED|=', en, &
            &'  max|d vs NEG-CONJ|=', ek
        call flush(logfhandle)
        deallocate(plc, plct, ref_c, ref_t)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_prep_check_f')
#endif
    end subroutine flex_gpu_prep_check_f

    !> fetch the premultiplied reconstruction planes (conj(T)y, ctfsq) packed on the
    !! unpadded lattice; the state-reconstruction prep unpacks these into fplane_type
    subroutine flex_gpu_prep_fetch_f( plc, plct, nrecords, hlo, hhi, klo )
        complex(sp), intent(inout), target, contiguous :: plc(:,:,:)
        real(sp),    intent(inout), target, contiguous :: plct(:,:,:)
        integer,     intent(in)    :: nrecords, hlo, hhi, klo
#ifdef USE_FLEX_CUDA
        integer :: ierr
        ierr = c_prep_fetch(c_loc(plc), c_loc(plct), int(nrecords,c_int), &
            &int(hlo,c_int), int(hhi,c_int), int(klo,c_int))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_prep_fetch failed')
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_prep_fetch_f')
#endif
    end subroutine flex_gpu_prep_fetch_f

    !> fetch the separated observation-model planes (image, CTF, ctfsq) packed on the
    !! unpadded lattice; the caller unpacks into fplane_type at the pf-multiple positions
    subroutine flex_gpu_prep_fetch_sep_f( plcy, plt, plct, nrecords, hlo, hhi, klo )
        complex(sp), intent(inout), target, contiguous :: plcy(:,:,:)
        real(sp),    intent(inout), target, contiguous :: plt(:,:,:), plct(:,:,:)
        integer,     intent(in)    :: nrecords, hlo, hhi, klo
#ifdef USE_FLEX_CUDA
        integer :: ierr
        ierr = c_prep_fetch_sep(c_loc(plcy), c_loc(plt), c_loc(plct), int(nrecords,c_int), &
            &int(hlo,c_int), int(hhi,c_int), int(klo,c_int))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_prep_fetch_sep failed')
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_prep_fetch_sep_f')
#endif
    end subroutine flex_gpu_prep_fetch_sep_f

    subroutine flex_gpu_prep_free_f
#ifdef USE_FLEX_CUDA
        integer :: ierr
        ierr = c_prep_free()
#endif
        l_prep_dev_ready = .false.
    end subroutine flex_gpu_prep_free_f

    !> fused E-step over the device-prepped resident planes
    subroutine flex_gpu_estep_batch_res_f( orientations, nrecords, frlims, nyq_unpd, stats )
        use simple_math, only: ceil_div, floor_div
        type(ori), intent(inout) :: orientations(:)
        integer,   intent(in)    :: nrecords, frlims(3,2), nyq_unpd
        real(dp),  intent(out), target, contiguous :: stats(:,:)
#ifdef USE_FLEX_CUDA
        real(sp), allocatable, target :: rotm(:,:,:)
        integer :: i, ierr, hlo, hhi, klo
        hlo = ceil_div (frlims(1,1), OSMPL_PAD_FAC); hhi = floor_div(frlims(1,2), OSMPL_PAD_FAC)
        klo = ceil_div (frlims(2,1), OSMPL_PAD_FAC)
        allocate(rotm(3,3,nrecords), source=0.0_sp)
        do i = 1, nrecords
            rotm(:,:,i) = orientations(i)%get_mat()
        end do
        stats = 0.d0
        ierr = c_estep_batch_res(c_loc(rotm), int(nrecords,c_int), int(hlo,c_int), &
            &int(hhi,c_int), int(klo,c_int), 0_c_int, &
            &int(nyq_unpd*(nyq_unpd+1),c_int), c_loc(stats))
        if( ierr /= 0 ) THROW_HARD('flex_gpu_estep_batch_res failed')
        deallocate(rotm)
#else
        THROW_HARD('SIMPLE was built without USE_FLEX_CUDA; flex_gpu_estep_batch_res_f')
#endif
    end subroutine flex_gpu_estep_batch_res_f

    !> banked coupled adjoint vs a CPU reference of the SAME math (align by in-plane angle with
    !! the polar former's unit-tap 2D KB stencil, per-direction linear combine, one splat per
    !! direction at the bank rotation). Everything after the align is linear, so agreement to
    !! fp32 roundoff is the expectation; the direction-quantization APPROXIMATION itself is
    !! science-gated on full runs, not unit-tested.
    subroutine test_flex_gpu_coupled_banked()
        use simple_ftiter, only: ftiter
        use simple_flex_reconstructor_latent_ops, only: latent_projection_weights
        integer, parameter :: TEST_BOX=16, TEST_NCOMP=3, TEST_NRECORDS=3, TEST_NDIR=2
        integer, parameter :: NPAIRS = (TEST_NCOMP*(TEST_NCOMP+1))/2
        type(reconstructor) :: rec_ref(TEST_NCOMP), rec_gpu(TEST_NCOMP)
        type(fplane_type)   :: fpls(TEST_NRECORDS)
        type(ori)    :: o_b
        type(ftiter) :: fit_pd
        type(kbinterpol) :: kbwin
        real(dp) :: dsc(TEST_NCOMP,TEST_NRECORDS), rsc(NPAIRS,TEST_NRECORDS)
        real     :: rmat_bank(3,3,TEST_NDIR), cav(TEST_NRECORDS), sav(TEST_NRECORDS)
        real     :: ang, ca, sa, hu, ku, w, wtap, lx, ly, lz
        real     :: wx(3), wy(3), wz(3), wwx(3), wwy(3), wwz(3)
        integer  :: win(2,3), dirid(TEST_NRECORDS)
        complex, allocatable :: al_c(:,:,:), cmb_c(:,:,:)
        real,    allocatable :: al_ct(:,:,:), cmb_r(:,:,:)
        real(sp), allocatable :: rho_ref(:,:,:,:), rho_gpu(:,:,:,:)
        logical  :: valid(TEST_NRECORDS)
        integer  :: lims_pd(3,2), exp_lb(3), exp_ub(3), exp_shape(3)
        integer  :: i, q, r, h, k, hp, kp, hx, ky, ix, iy, iz, id, nyq_eff
        integer  :: hlo, hhi, klo, khi, wx0, wy0, wz0, xg, yg, zg
        complex  :: cv, alv
        real     :: ctv, cerr, rerr
        if( .not. flex_gpu_available() )then
            write(logfhandle,'(A)') '>>> TEST flex_gpu banked SKIPPED (no USE_FLEX_CUDA build or no device)'
            return
        endif
        fit_pd  = ftiter([OSMPL_PAD_FAC*TEST_BOX, OSMPL_PAD_FAC*TEST_BOX, 1], 1.)
        lims_pd = fit_pd%loop_lims(3)
        exp_lb  = [-TEST_BOX/2-2, -TEST_BOX/2-2, -TEST_BOX/2-2]
        exp_ub  = [ TEST_BOX/2+2,  TEST_BOX/2+2,  TEST_BOX/2+2]
        exp_shape = exp_ub - exp_lb + 1
        do q = 1, TEST_NCOMP
            call rec_ref(q)%new([TEST_BOX,TEST_BOX,TEST_BOX], 1.)
            call rec_gpu(q)%new([TEST_BOX,TEST_BOX,TEST_BOX], 1.)
            allocate(rec_ref(q)%cmat_exp(exp_lb(1):exp_ub(1),exp_lb(2):exp_ub(2),exp_lb(3):exp_ub(3)), &
                &source=CMPLX_ZERO)
            allocate(rec_gpu(q)%cmat_exp(exp_lb(1):exp_ub(1),exp_lb(2):exp_ub(2),exp_lb(3):exp_ub(3)), &
                &source=CMPLX_ZERO)
        end do
        allocate(rho_ref(NPAIRS,exp_shape(1),exp_shape(2),exp_shape(3)), source=0.)
        allocate(rho_gpu(NPAIRS,exp_shape(1),exp_shape(2),exp_shape(3)), source=0.)
        kbwin = rec_gpu(1)%get_kbwin()
        ! two bank directions; records 1,2 share dir 1 with different in-plane angles
        call o_b%new(.true.)
        do id = 1, TEST_NDIR
            call o_b%set_euler([23.*real(id), 37.*real(id), 0.])
            rmat_bank(:,:,id) = o_b%get_mat()
        end do
        call o_b%kill
        dirid = [1, 1, 2]
        do i = 1, TEST_NRECORDS
            ang    = 0.31*real(i)
            cav(i) = cos(ang); sav(i) = sin(ang)
            allocate(fpls(i)%cmplx_plane(lims_pd(1,1):lims_pd(1,2), lims_pd(2,1):0), source=CMPLX_ZERO)
            allocate(fpls(i)%ctfsq_plane(lims_pd(1,1):lims_pd(1,2), lims_pd(2,1):0), source=0.)
            allocate(fpls(i)%transfer_plane(lims_pd(1,1):lims_pd(1,2), lims_pd(2,1):0), source=CMPLX_ZERO)
            fpls(i)%frlims = lims_pd
            fpls(i)%nyq    = fit_pd%get_lfny(1)
            do k = -TEST_BOX/2, 0
                kp = OSMPL_PAD_FAC*k
                do h = -TEST_BOX/2, TEST_BOX/2
                    if( h*h + k*k > (TEST_BOX/2)*(TEST_BOX/2+1) ) cycle
                    hp = OSMPL_PAD_FAC*h
                    fpls(i)%cmplx_plane(hp,kp)    = cmplx(0.01*real(3*h-2*k+i), 0.02*real(h+k-i))
                    fpls(i)%transfer_plane(hp,kp) = cmplx(0.7+0.01*real(mod(abs(h+2*k+i),7)), 0.)
                    fpls(i)%ctfsq_plane(hp,kp)    = 0.2 + 0.03*real(mod(abs(2*h-k+i),9))
                end do
            end do
            do q = 1, TEST_NCOMP
                dsc(q,i) = 0.1d0*real(q,dp) - 0.25d0*real(i,dp)/real(q,dp)
            end do
            do r = 1, NPAIRS
                rsc(r,i) = 0.05d0*real(r+i,dp)
            end do
        end do
        valid = .true.
        ! ---- GPU leg ----
        call flex_gpu_coupled_begin_f(rec_gpu, TEST_NCOMP, NPAIRS)
        call flex_gpu_coupled_bank_f(rmat_bank, TEST_NDIR)
        call flex_gpu_coupled_batch_banked_f(fpls, dsc, rsc, valid, TEST_NRECORDS, dirid, cav, sav)
        call flex_gpu_coupled_end_f(rec_gpu, rho_gpu)
        call flex_gpu_coupled_bank_free_f
        ! ---- CPU reference ----
        hlo = ceil_div_t(lims_pd(1,1)); hhi = floor_div_t(lims_pd(1,2))
        klo = ceil_div_t(lims_pd(2,1)); khi = floor_div_t(lims_pd(2,2))
        nyq_eff = min(rec_gpu(1)%get_lfny(1), max(1, fpls(1)%nyq / OSMPL_PAD_FAC))
        allocate(al_c (hlo:hhi, klo:khi, TEST_NRECORDS), source=CMPLX_ZERO)
        allocate(al_ct(hlo:hhi, klo:khi, TEST_NRECORDS), source=0.)
        do i = 1, TEST_NRECORDS
            ca = cav(i); sa = sav(i)
            do k = klo, khi
                do h = hlo, hhi
                    if( h*h + k*k > nyq_eff*(nyq_eff+1) ) cycle
                    hu =  h*ca + k*sa
                    ku = -h*sa + k*ca
                    call latent_projection_weights(kbwin, [hu, ku, 0.], win, wx, wy, wz)
                    alv = CMPLX_ZERO; ctv = 0.
                    do iy = 1, 3
                        ky = win(1,2) + iy - 1
                        do ix = 1, 3
                            hx = win(1,1) + ix - 1
                            w  = wx(ix)*wy(iy)
                            kp = OSMPL_PAD_FAC*ky; hp = OSMPL_PAD_FAC*hx
                            if( kp <= 0 )then
                                if( hp < lims_pd(1,1) .or. hp > lims_pd(1,2) .or. kp < lims_pd(2,1) ) cycle
                                cv  = conjg(fpls(i)%transfer_plane(hp,kp))*fpls(i)%cmplx_plane(hp,kp)
                                ctv = ctv + w*fpls(i)%ctfsq_plane(hp,kp)
                            else
                                if( -hp < lims_pd(1,1) .or. -hp > lims_pd(1,2) .or. -kp < lims_pd(2,1) ) cycle
                                cv  = conjg(conjg(fpls(i)%transfer_plane(-hp,-kp))*fpls(i)%cmplx_plane(-hp,-kp))
                                ctv = ctv + w*fpls(i)%ctfsq_plane(-hp,-kp)
                            endif
                            alv = alv + w*cv
                        end do
                    end do
                    al_c (h,k,i) = real(OSMPL_PAD_FAC*OSMPL_PAD_FAC)*alv
                    al_ct(h,k,i) = ctv
                end do
            end do
        end do
        allocate(cmb_c(hlo:hhi, klo:khi, TEST_NCOMP), cmb_r(hlo:hhi, klo:khi, NPAIRS))
        do id = 1, TEST_NDIR
            cmb_c = CMPLX_ZERO; cmb_r = 0.
            do i = 1, TEST_NRECORDS
                if( dirid(i) /= id ) cycle
                do q = 1, TEST_NCOMP
                    cmb_c(:,:,q) = cmb_c(:,:,q) + real(dsc(q,i))*al_c(:,:,i)
                end do
                do r = 1, NPAIRS
                    cmb_r(:,:,r) = cmb_r(:,:,r) + real(rsc(r,i))*al_ct(:,:,i)
                end do
            end do
            do k = klo, khi
                do h = hlo, hhi
                    if( h*h + k*k > nyq_eff*(nyq_eff+1) ) cycle
                    lx = h*rmat_bank(1,1,id) + k*rmat_bank(2,1,id)
                    ly = h*rmat_bank(1,2,id) + k*rmat_bank(2,2,id)
                    lz = h*rmat_bank(1,3,id) + k*rmat_bank(2,3,id)
                    call latent_projection_weights(kbwin, [lx, ly, lz], win, wwx, wwy, wwz)
                    wx0 = win(1,1); wy0 = win(1,2); wz0 = win(1,3)
                    if( wx0 + 2 < 0 ) cycle
                    if( wx0 < exp_lb(1) .or. wy0 < exp_lb(2) .or. wz0 < exp_lb(3) ) cycle
                    if( wx0+2 > exp_ub(1) .or. wy0+2 > exp_ub(2) .or. wz0+2 > exp_ub(3) ) cycle
                    do iz = 1, 3
                        zg = wz0 + iz - 1
                        do iy = 1, 3
                            yg = wy0 + iy - 1
                            do ix = 1, 3
                                xg   = wx0 + ix - 1
                                wtap = wwx(ix)*wwy(iy)*wwz(iz)
                                do q = 1, TEST_NCOMP
                                    rec_ref(q)%cmat_exp(xg,yg,zg) = rec_ref(q)%cmat_exp(xg,yg,zg) &
                                        &+ wtap*cmb_c(h,k,q)
                                end do
                                do r = 1, NPAIRS
                                    rho_ref(r, xg-exp_lb(1)+1, yg-exp_lb(2)+1, zg-exp_lb(3)+1) = &
                                        &rho_ref(r, xg-exp_lb(1)+1, yg-exp_lb(2)+1, zg-exp_lb(3)+1) &
                                        &+ wtap*cmb_r(h,k,r)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
        cerr = 0.
        do q = 1, TEST_NCOMP
            cerr = max(cerr, maxval(abs(rec_ref(q)%cmat_exp - rec_gpu(q)%cmat_exp)))
        end do
        rerr = maxval(abs(rho_ref - rho_gpu))
        write(logfhandle,'(A,ES10.2,A,ES10.2)') '>>> TEST flex_gpu banked: max|dcmat|=', cerr, &
            &'  max|drho|=', rerr
        if( cerr > 2.e-5 .or. rerr > 2.e-5 ) THROW_HARD('GPU banked coupled adjoint differs from CPU reference')
        write(logfhandle,'(A)') '>>>   PASSED (banked align+combine+splat matches CPU reference)'
        do i = 1, TEST_NRECORDS
            if( allocated(fpls(i)%cmplx_plane) )    deallocate(fpls(i)%cmplx_plane)
            if( allocated(fpls(i)%ctfsq_plane) )    deallocate(fpls(i)%ctfsq_plane)
            if( allocated(fpls(i)%transfer_plane) ) deallocate(fpls(i)%transfer_plane)
        end do
        do q = 1, TEST_NCOMP
            call rec_ref(q)%kill
            call rec_gpu(q)%kill
        end do
        deallocate(rho_ref, rho_gpu, al_c, al_ct, cmb_c, cmb_r)
    contains
        integer function ceil_div_t( a )
            integer, intent(in) :: a
            ceil_div_t = a/OSMPL_PAD_FAC
            if( mod(a, OSMPL_PAD_FAC) /= 0 .and. a > 0 ) ceil_div_t = ceil_div_t + 1
        end function ceil_div_t
        integer function floor_div_t( a )
            integer, intent(in) :: a
            floor_div_t = a/OSMPL_PAD_FAC
            if( mod(a, OSMPL_PAD_FAC) /= 0 .and. a < 0 ) floor_div_t = floor_div_t - 1
        end function floor_div_t
    end subroutine test_flex_gpu_coupled_banked

    !> device polar sampler vs the CPU polar former on synthetic planes: same grid, same
    !! in-plane angles, xw/halves/wr/noise-power/tazim compared. The device apod is the
    !! 14-term polynomial (apod_fast); the CPU one is exact -- tolerance covers the ~1e-6 gap.
    subroutine test_flex_gpu_psample()
        use simple_ftiter, only: ftiter
        use simple_flex_pca_polar, only: polar_grid_t, polar_grid_build, polar_grid_kill, &
            &polar_sample_particle
        integer, parameter :: TEST_BOX=32, TEST_NREC=3
        type(polar_grid_t) :: pg
        type(fplane_type)  :: fpls(TEST_NREC)
        type(ftiter)       :: fit_pd
        complex, allocatable :: xw(:), xw1(:), xw2(:)
        real,    allocatable :: xws_ref(:,:), xws1_ref(:,:), xws2_ref(:,:), wr_ref(:,:)
        real,    allocatable :: xws_gpu(:,:), xws1_gpu(:,:), xws2_gpu(:,:), wr_gpu(:,:)
        real,    allocatable :: pw_gpu(:), taz_gpu(:)
        real     :: cav(TEST_NREC), sav(TEST_NREC), taz_ref(TEST_NREC), ang, tazim
        real(dp) :: pw_ref(TEST_NREC), cnt_ref
        logical  :: valid(TEST_NREC)
        integer  :: lims_pd(3,2), i, j, h, k, hp, kp, kfrom, kto, knfrom, knto, hlo, hhi, klo
        real     :: e1, e2, e3, e4
        if( .not. flex_gpu_available() )then
            write(logfhandle,'(A)') '>>> TEST flex_gpu psample SKIPPED (no USE_FLEX_CUDA build or no device)'
            return
        endif
        fit_pd  = ftiter([OSMPL_PAD_FAC*TEST_BOX, OSMPL_PAD_FAC*TEST_BOX, 1], 1.)
        lims_pd = fit_pd%loop_lims(3)
        hlo = lims_pd(1,1)/OSMPL_PAD_FAC; hhi = lims_pd(1,2)/OSMPL_PAD_FAC
        klo = lims_pd(2,1)/OSMPL_PAD_FAC
        kfrom = 3; kto = TEST_BOX/2 - 4; knfrom = kto + 1; knto = TEST_BOX/2 - 1
        call polar_grid_build(pg, kfrom, kto, knfrom, knto, hlo, hhi, klo, lims_pd(1,1), lims_pd(2,1))
        do i = 1, TEST_NREC
            ang    = 0.47*real(i) - 0.3
            cav(i) = cos(ang); sav(i) = sin(ang)
            allocate(fpls(i)%cmplx_plane(lims_pd(1,1):lims_pd(1,2), lims_pd(2,1):0), source=CMPLX_ZERO)
            allocate(fpls(i)%ctfsq_plane(lims_pd(1,1):lims_pd(1,2), lims_pd(2,1):0), source=0.)
            allocate(fpls(i)%transfer_plane(lims_pd(1,1):lims_pd(1,2), lims_pd(2,1):0), source=CMPLX_ZERO)
            fpls(i)%frlims = lims_pd
            fpls(i)%nyq    = fit_pd%get_lfny(1)
            do k = -TEST_BOX/2, 0
                kp = OSMPL_PAD_FAC*k
                do h = -TEST_BOX/2, TEST_BOX/2
                    if( h*h + k*k > (TEST_BOX/2)*(TEST_BOX/2+1) ) cycle
                    hp = OSMPL_PAD_FAC*h
                    fpls(i)%cmplx_plane(hp,kp)    = cmplx(0.01*real(2*h-k+i), 0.015*real(h+2*k-i))
                    fpls(i)%transfer_plane(hp,kp) = cmplx(0.6+0.02*real(mod(abs(h-k+i),5)), &
                        &0.05*real(mod(abs(h+k),3)))
                end do
            end do
        end do
        valid = [.true., .true., .true.]
        allocate(xws_ref(2*pg%nsamp,TEST_NREC), xws1_ref(2*pg%nsamp,TEST_NREC), &
            &xws2_ref(2*pg%nsamp,TEST_NREC), wr_ref(pg%nk,TEST_NREC), source=0.)
        allocate(xws_gpu(2*pg%nsamp,TEST_NREC), xws1_gpu(2*pg%nsamp,TEST_NREC), &
            &xws2_gpu(2*pg%nsamp,TEST_NREC), wr_gpu(pg%nk,TEST_NREC), source=0.)
        allocate(pw_gpu(TEST_NREC), taz_gpu(TEST_NREC))
        allocate(xw(pg%nsamp), xw1(pg%nsamp), xw2(pg%nsamp))
        do i = 1, TEST_NREC
            call polar_sample_particle(fpls(i)%cmplx_plane, fpls(i)%transfer_plane, pg, &
                &cav(i), sav(i), xw, wr_ref(:,i), pw_ref(i), cnt_ref, taz_ref(i), xw1, xw2)
            do j = 1, pg%nsamp
                xws_ref (2*j-1,i) = pg%sqwq(j)*real(xw(j));  xws_ref (2*j,i) = pg%sqwq(j)*aimag(xw(j))
                xws1_ref(2*j-1,i) = pg%sqwq(j)*real(xw1(j)); xws1_ref(2*j,i) = pg%sqwq(j)*aimag(xw1(j))
                xws2_ref(2*j-1,i) = pg%sqwq(j)*real(xw2(j)); xws2_ref(2*j,i) = pg%sqwq(j)*aimag(xw2(j))
            end do
        end do
        call flex_gpu_psample_begin_f(pg%rad, pg%cs, pg%sn, pg%sqwq, pg%rbeg, pg%rend, &
            &pg%nsamp, pg%nk, pg%nrad, pg%ncs, pg%nsn, pg%nwq, pg%nsamp_n)
        call flex_gpu_psample_batch_f(fpls, cav, sav, valid, TEST_NREC, .true., &
            &xws_gpu, xws1_gpu, xws2_gpu, wr_gpu, 1, pw_gpu, taz_gpu)
        call flex_gpu_psample_free_f
        e1 = maxval(abs(xws_ref - xws_gpu))
        e2 = max(maxval(abs(xws1_ref - xws1_gpu)), maxval(abs(xws2_ref - xws2_gpu)))
        e3 = maxval(abs(wr_ref - wr_gpu))
        e4 = 0.
        do i = 1, TEST_NREC
            ! CPU tazim comes back ALREADY divided by nk; the device returns the raw ring sum
            e4 = max(e4, abs(real(pw_ref(i)) - pw_gpu(i))/max(abs(real(pw_ref(i))), 1.e-6), &
                &abs(taz_ref(i) - taz_gpu(i)/real(max(1,pg%nk))))
        end do
        write(logfhandle,'(A,ES10.2,A,ES10.2,A,ES10.2,A,ES10.2)') '>>> TEST flex_gpu psample: max|dxw|=', &
            &e1, '  max|dhalves|=', e2, '  max|dwr|=', e3, '  rel|dpw|,|dtaz|=', e4
        if( e1 > 5.e-5 .or. e2 > 5.e-5 .or. e3 > 5.e-5 .or. e4 > 5.e-4 ) &
            &THROW_HARD('GPU polar sampler differs from CPU polar former')
        write(logfhandle,'(A)') '>>>   PASSED (device polar sampler matches the CPU former)'
        do i = 1, TEST_NREC
            deallocate(fpls(i)%cmplx_plane, fpls(i)%ctfsq_plane, fpls(i)%transfer_plane)
        end do
        call polar_grid_kill(pg)
        deallocate(xws_ref, xws1_ref, xws2_ref, wr_ref, xws_gpu, xws1_gpu, xws2_gpu, wr_gpu)
        deallocate(pw_gpu, taz_gpu, xw, xw1, xw2)
    end subroutine test_flex_gpu_psample

    !> fused E-step vs the CPU reference (project_fplanes_mean_basis + cov_herm_inner) on
    !! synthetic volumes/planes at random exact orientations
    subroutine test_flex_gpu_estep()
        use simple_ftiter, only: ftiter
        use simple_flex_reconstructor_latent_ops, only: project_fplanes_mean_basis
        integer, parameter :: TEST_BOX=16, TEST_NCOMP=3, TEST_NREC=3
        type(reconstructor) :: mean_rec, basis_recs(TEST_NCOMP)
        type(fplane_type)   :: fpls(TEST_NREC), mean_fpl, basis_fpls(TEST_NCOMP)
        type(ori)    :: orientations(TEST_NREC)
        type(ftiter) :: fit_pd
        real(dp), allocatable :: stats(:,:)
        real(dp) :: e_ref, m_ref, b_ref(TEST_NCOMP), c_ref(TEST_NCOMP), G_ref(TEST_NCOMP,TEST_NCOMP)
        logical  :: valid(TEST_NREC)
        integer  :: lims_pd(3,2), exp_lb(3), exp_ub(3), i, q, r, h, k, hp, kp, nst, off
        real     :: err
        if( .not. flex_gpu_available() )then
            write(logfhandle,'(A)') '>>> TEST flex_gpu estep SKIPPED (no USE_FLEX_CUDA build or no device)'
            return
        endif
        fit_pd  = ftiter([OSMPL_PAD_FAC*TEST_BOX, OSMPL_PAD_FAC*TEST_BOX, 1], 1.)
        lims_pd = fit_pd%loop_lims(3)
        exp_lb  = [-TEST_BOX/2-2, -TEST_BOX/2-2, -TEST_BOX/2-2]
        exp_ub  = [ TEST_BOX/2+2,  TEST_BOX/2+2,  TEST_BOX/2+2]
        call mean_rec%new([TEST_BOX,TEST_BOX,TEST_BOX], 1.)
        allocate(mean_rec%cmat_exp(exp_lb(1):exp_ub(1),exp_lb(2):exp_ub(2),exp_lb(3):exp_ub(3)), &
            &source=CMPLX_ZERO)
        do q = 1, TEST_NCOMP
            call basis_recs(q)%new([TEST_BOX,TEST_BOX,TEST_BOX], 1.)
            allocate(basis_recs(q)%cmat_exp(exp_lb(1):exp_ub(1),exp_lb(2):exp_ub(2),exp_lb(3):exp_ub(3)), &
                &source=CMPLX_ZERO)
        end do
        do k = exp_lb(3), exp_ub(3)
            do h = exp_lb(1), exp_ub(1)
                do r = exp_lb(2), exp_ub(2)
                    mean_rec%cmat_exp(h,r,k) = cmplx(0.02*real(h-r+k), 0.01*real(h+r))
                    do q = 1, TEST_NCOMP
                        basis_recs(q)%cmat_exp(h,r,k) = cmplx(0.01*real(q*h+r), 0.015*real(k-q*r))
                    end do
                end do
            end do
        end do
        do i = 1, TEST_NREC
            call orientations(i)%new(.true.)
            call orientations(i)%set_euler([25.*real(i)+7., 40.*real(i)+11., 13.*real(i)])
            allocate(fpls(i)%cmplx_plane(lims_pd(1,1):lims_pd(1,2), lims_pd(2,1):0), source=CMPLX_ZERO)
            allocate(fpls(i)%ctfsq_plane(lims_pd(1,1):lims_pd(1,2), lims_pd(2,1):0), source=0.)
            allocate(fpls(i)%transfer_plane(lims_pd(1,1):lims_pd(1,2), lims_pd(2,1):0), source=CMPLX_ZERO)
            fpls(i)%frlims = lims_pd
            fpls(i)%nyq    = fit_pd%get_lfny(1)
            do k = -TEST_BOX/2, 0
                kp = OSMPL_PAD_FAC*k
                do h = -TEST_BOX/2, TEST_BOX/2
                    if( h*h + k*k > (TEST_BOX/2)*(TEST_BOX/2+1) ) cycle
                    hp = OSMPL_PAD_FAC*h
                    fpls(i)%cmplx_plane(hp,kp)    = cmplx(0.01*real(3*h-2*k+i), 0.02*real(h+k-i))
                    fpls(i)%transfer_plane(hp,kp) = cmplx(0.7+0.01*real(mod(abs(h+2*k+i),7)), 0.)
                    fpls(i)%ctfsq_plane(hp,kp)    = real(fpls(i)%transfer_plane(hp,kp))**2
                end do
            end do
        end do
        valid = .true.
        ! GPU leg (g2_lfny is read by the batch wrapper; set it via a coupled begin on the recs)
        call flex_gpu_coupled_begin_f(basis_recs, TEST_NCOMP, TEST_NCOMP)
        call flex_gpu_estep_vols_f(mean_rec, basis_recs, TEST_NCOMP)
        nst = 2 + 2*TEST_NCOMP + (TEST_NCOMP*(TEST_NCOMP+1))/2
        allocate(stats(nst,TEST_NREC))
        call flex_gpu_estep_batch_f(fpls, orientations, valid, TEST_NREC, stats)
        call flex_gpu_estep_free_f
        call flex_gpu_coupled_end_f_noop()
        ! CPU reference
        err = 0.
        do i = 1, TEST_NREC
            call project_fplanes_mean_basis(mean_rec, basis_recs, orientations(i), fpls(i), &
                &mean_fpl, basis_fpls, apply_ctf_amp=.true.)
            e_ref = real(cov_herm_inner_t(mean_fpl, mean_fpl), dp)
            m_ref = real(cov_herm_inner_t(mean_fpl, fpls(i)), dp)
            do q = 1, TEST_NCOMP
                b_ref(q) = real(cov_herm_inner_t(basis_fpls(q), fpls(i)), dp)
                c_ref(q) = real(cov_herm_inner_t(basis_fpls(q), mean_fpl), dp)
                do r = q, TEST_NCOMP
                    G_ref(q,r) = real(cov_herm_inner_t(basis_fpls(q), basis_fpls(r)), dp)
                end do
            end do
            err = max(err, abs(real(e_ref - stats(1,i))), abs(real(m_ref - stats(2,i))))
            do q = 1, TEST_NCOMP
                err = max(err, abs(real(b_ref(q) - stats(2+q,i))))
                err = max(err, abs(real(c_ref(q) - stats(2+TEST_NCOMP+q,i))))
            end do
            off = 2 + 2*TEST_NCOMP
            do q = 1, TEST_NCOMP
                do r = q, TEST_NCOMP
                    off = off + 1
                    err = max(err, abs(real(G_ref(q,r) - stats(off,i))))
                end do
            end do
        end do
        write(logfhandle,'(A,ES10.2)') '>>> TEST flex_gpu estep: max|dstat|=', err
        if( err > 5.e-4 ) THROW_HARD('GPU fused E-step differs from CPU reference')
        write(logfhandle,'(A)') '>>>   PASSED (fused E-step stats match projector + inner products)'
        do i = 1, TEST_NREC
            deallocate(fpls(i)%cmplx_plane, fpls(i)%ctfsq_plane, fpls(i)%transfer_plane)
            call orientations(i)%kill
        end do
        call mean_rec%kill
        do q = 1, TEST_NCOMP
            call basis_recs(q)%kill
        end do
        deallocate(stats)
    contains
        !> local copy of cov_herm_inner's sum convention (the module one lives in columns)
        function cov_herm_inner_t( lhs, rhs ) result( val )
            type(fplane_type), intent(in) :: lhs, rhs
            complex(dp) :: val
            integer :: h2, k2, nyq2, nd2
            val = cmplx(0.d0,0.d0,dp)
            nyq2 = min(lhs%nyq, rhs%nyq)
            nd2  = nyq2*(nyq2+1)
            do k2 = -TEST_BOX/2, 0
                do h2 = -TEST_BOX/2, merge(0, TEST_BOX/2, k2 == 0)
                    if( h2*h2 + k2*k2 > nd2 ) cycle
                    val = val + conjg(cmplx(lhs%cmplx_plane(OSMPL_PAD_FAC*h2,OSMPL_PAD_FAC*k2),kind=dp)) &
                        &* cmplx(rhs%cmplx_plane(OSMPL_PAD_FAC*h2,OSMPL_PAD_FAC*k2),kind=dp)
                end do
            end do
        end function cov_herm_inner_t
        subroutine flex_gpu_coupled_end_f_noop()
#ifdef USE_FLEX_CUDA
            integer :: ierr
            ierr = c_coupled_free()
#endif
        end subroutine
    end subroutine test_flex_gpu_estep

end module simple_flex_gpu
