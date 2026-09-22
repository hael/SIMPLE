!@descr: flex_pca EM: the M-step (coupled solve, FSC-Wiener merge, re-orthonormalisation, deflation, mixture update) per fit and iteration
submodule (simple_flex_pca_em) simple_flex_pca_em_mstep
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_flex_pca_util, only: dilation_template
use simple_matcher_3Drec,   only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch
use simple_flex_reconstructor_latent_ops, only: project_fplane_mean, project_fplanes_mean_basis,&
    &insert_planes_oversamp_coupled_batch_scaled
use simple_flex_reconstructor_latent_ops, only: prep_imgs4projected_model, solve_coupled_basis_exp,&
    &projected_model_kfromto, add_invtausq2rho_coupled
use simple_flex_pca_crossfsc, only: crossfsc_file, crossfsc_record, crossfsc_load, crossfsc_write,&
    &crossfsc_append, crossfsc_latest_upto, crossfsc_kill, crossfsc_kill_record, crossfsc_to_invtau2,&
    &crossfsc_harvest_h, crossfsc_stop_stat, crossfsc_inband_mean, crossfsc_khi_deepest,&
    &crossfsc_assert_paired, COV_XFSC_FNAME
use simple_flex_gpu,        only: flex_gpu_available, flex_gpu_coupled_begin_f,&
    &flex_gpu_coupled_batch_raw_f, flex_gpu_coupled_end_f, flex_gpu_coupled_bank_f,&
    &flex_gpu_coupled_batch_banked_f, flex_gpu_coupled_bank_free_f, flex_gpu_estep_vols_f,&
    &flex_gpu_estep_batch_f, flex_gpu_estep_resid_f, flex_gpu_estep_free_f,&
    &flex_gpu_coupled_batch_banked_res_f, flex_gpu_prep_begin_f, flex_gpu_prep_batch_f,&
    &flex_gpu_prep_free_f, flex_gpu_estep_batch_res_f, flex_gpu_prep_check_f,&
    &flex_gpu_poles_begin_f, flex_gpu_poles_bank_f, flex_gpu_poles_batch_f, flex_gpu_poles_free_f
use simple_flex_pca_polar,  only: polar_grid_build, polar_grid_kill, polar_project_recs,&
    &polar_relative_inplane, polar_assign_directions, polar_sample_particle_fused
implicit none
#include "simple_local_flags.inc"

contains

    !> The master-only tail of one EM iteration for one fit, in the single-fit order: Gamma
    !! accumulate, likelihood log, coupled per-halfset M-step
    !! solve + FSC-Wiener merge + band-limit + mask, BAND/RANK diagnostic (writes khi_fit),
    !! even/odd update-agreement bookkeeping, mean-shaped deflation, re-orthonormalisation,
    !! principal-angle convergence, basis swap + eigenvolume write, Gamma log, and the fit's
    !! per-iteration frees.
    module subroutine fit_iter_finish( params, build, fit, it_eff, nthr )
        class(parameters),  intent(inout) :: params
        type(builder),      intent(inout) :: build
        type(probe_fit_t),  intent(inout) :: fit
        integer,            intent(in)    :: it_eff, nthr
        type(reconstructor), allocatable :: utilde(:)
        type(image), allocatable :: realvols(:), utilde_real(:), eimgs(:), oimgs(:), dfl_basis(:)
        type(image)  :: img_o, mstep_gridcorr, mvol_dfl
        type(string) :: fname
        real,     allocatable :: filt(:), corrs(:), fscq_dg(:,:)
        real(dp), allocatable :: sv_eo(:), Mconv(:,:), sconv(:)
        real,     pointer :: rmatp(:,:,:), rm_dfl(:,:,:), rv_dfl(:,:,:)
        real(dp) :: eo_dim, cos_mean, mu_q, sd_q
        real(dp) :: mm_dfl, mv_dfl, rem_dfl, tot_dfl, mnorm_dfl
        real(dp) :: fmean_dg(512), fbest_dg
        real     :: fc, res_lo, res_hi, res_dg
        integer  :: q, ithr, sh, filtsz, d_new
        integer  :: ndfl, ndfl_sh, idfl, jdfl, nkeep_dfl, kfr_dfl(2)
        integer  :: khi_dg, nsig_dg, ntop_dg, sel_dg(4), tq_dg, bq_dg
        logical  :: l_dfl_bg, l_dfl_dil, l_z_local
        logical  :: l_mstep_wiener
        integer  :: iblk_dfl, nblk_dfl, qlo_dfl, qhi_dfl
        type(flex_pcg_outcome_t) :: pcg_out
        do q = 1, fit%ncomp
            fit%gam_acc(q) = fit%gam_sum(q) / real(max(1,fit%nval),dp)
        end do
        deallocate(fit%gam_sum)
        ! Marginal likelihood (the EM objective):
        !   -2 log p(y_i) = ||y_i - a_i T_i u_0||^2/sig2 - h_i' A_i^-1 h_i + log det A_i + log det Gamma
        ! The first term is fixed across iterations (a_i and mu are held), so only the rest is
        ! accumulated; log det Gamma uses `prior`, the Gamma the E-step just used. The filtered
        ! M-step sits outside the likelihood, so this is a diagnostic, not a monotonicity guarantee.
        if( fit%l_mix_used )then
            ! the mixture marginal's Omega log det, with the Omega this iteration's E-step used
            fit%nll_tot = fit%nll_tot + real(fit%nval,dp)*fit%ldOm_used
        else
            fit%nll_tot = fit%nll_tot - real(fit%nval,dp)*sum(log(max(fit%prior(1:fit%ncomp), DTINY)))
        endif
        write(logfhandle,'(A,I0,A,ES14.6,A,ES12.4)') '>>> FLEX_PCA PROBE ITER ',it_eff, &
            &'  -2logL/N (varying part)=',fit%nll_tot/real(max(1,fit%nval),dp), &
            &'  delta=',(fit%nll_tot/real(max(1,fit%nval),dp)) - fit%nll_prev
        fit%nll_prev = fit%nll_tot/real(max(1,fit%nval),dp)
        call flush(logfhandle)
        ! Cross-FSC sampling profiles H: per-shell means of the rho_e/rho_o diagonal pair rows,
        ! harvested before any invtau2 mutates them. On the single-fit engine h_a/h_b are the
        ! internal even/odd halves; the paired writer sums e+o per fit.
        if( fit%l_xf_harvest )then
            block
                integer :: xf_lb(3), xf_nyq, xf_fsz
                xf_lb  = lbound(fit%Yeven(1)%cmat_exp)
                xf_nyq = fit%Yeven(1)%get_lfny(1)
                xf_fsz = max(1, fdim(params%box_crop) - 1)
                if( allocated(fit%xf_h_e) ) deallocate(fit%xf_h_e)
                if( allocated(fit%xf_h_o) ) deallocate(fit%xf_h_o)
                if( allocated(fit%xf_cnt) ) deallocate(fit%xf_cnt)
                allocate(fit%xf_h_e(xf_fsz,fit%ncomp), fit%xf_h_o(xf_fsz,fit%ncomp), &
                    &fit%xf_cnt(xf_fsz))
                call crossfsc_harvest_h(fit%rho_e, fit%npairs, fit%ncomp, xf_lb, xf_nyq, &
                    &xf_fsz, fit%xf_h_e, fit%xf_cnt)
                call crossfsc_harvest_h(fit%rho_o, fit%npairs, fit%ncomp, xf_lb, xf_nyq, &
                    &xf_fsz, fit%xf_h_o, fit%xf_cnt)
            end block
        endif
        ! Cross-FSC SSNR shrinkage ridge (SIMPLE_COV_XFSC_REG=1|2), the flex analog of
        ! add_invtausq2rho: invtau2 from the previous iteration's record, built by the driver
        ! (xfsc_prep_ridge), added to the diagonal rows of rho_e and rho_o after any distributed
        ! reduction and immediately before the solves. One-shot: the driver rebuilds it each iteration.
        if( allocated(fit%xf_invtau2) )then
            call add_invtausq2rho_coupled(fit%Yeven, fit%rho_e, fit%ncomp, fit%xf_invtau2)
            call add_invtausq2rho_coupled(fit%Yodd,  fit%rho_o, fit%ncomp, fit%xf_invtau2)
            ! the PCG operator carries the same ridge as the preconditioner's density
            if( fit%l_pcg ) call fit%pcg%set_ridge(fit%xf_invtau2)
            deallocate(fit%xf_invtau2)
        endif
        ! Coupled M-step solve: at every grid point the components share one k x k normal matrix
        ! sum_i |CTF|^2 E[z_i z_i'], so the basis volumes are solved together per halfset. A unit
        ! density is then handed to compress_exp (the divide has already happened).
        if( fit%l_pcg )then
            ! rec_backend=pcg: the per-voxel solve is the preconditioner of the coupled normal equations
            ! on the pair Gram kernels; each half is solved from it on the spherical support
            call fit%pcg%finalize(fit%kpk_e)
            call fit%pcg%solve(fit%Yeven, fit%rho_e, fit%rpk_e, params%maxits_pcg, params%rtol, pcg_out, &
                &'FLEX_PCA PCG MSTEP even')
            call log_pcg_outcome('even', pcg_out)
            call fit%pcg%finalize(fit%kpk_o)
            call fit%pcg%solve(fit%Yodd,  fit%rho_o, fit%rpk_o, params%maxits_pcg, params%rtol, pcg_out, &
                &'FLEX_PCA PCG MSTEP odd')
            call log_pcg_outcome('odd', pcg_out)
            call fit%pcg%clear_ridge
        else
            call solve_coupled_basis_exp(fit%Yeven, fit%rho_e, fit%ncomp)
            call solve_coupled_basis_exp(fit%Yodd,  fit%rho_o, fit%ncomp)
        endif
        ! finalize even/odd Y_q, half-set FSC Wiener-merge -> band-limited masked basis volume;
        ! the two half-bases are kept for the update-agreement diagnostic below
        allocate(realvols(fit%ncomp))
        allocate(eimgs(fit%ncomp), oimgs(fit%ncomp))
        filtsz = max(1, fdim(params%box_crop) - 1)
        allocate(filt(filtsz), corrs(filtsz))
        l_mstep_wiener = .true.   ! measured 2026-09-11: shipping the unshrunk basis adds noise, not resolution
        if( .not. l_mstep_wiener .and. it_eff == 1 )then
            write(logfhandle,'(A)') '>>> FLEX_PCA M-STEP: half-set FSC Wiener filter OFF (SIMPLE_COV_MSTEP_WIENER=0); &
                &components band-limited and windowed only'
            call flush(logfhandle)
        endif
        ! native-lattice inverse KB envelope, applied to each half before FSC/merge/mask
        mstep_gridcorr = prep3D_inv_kbenvelope4mul([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        allocate(fscq_dg(filtsz, fit%ncomp), source=0.)
        do q = 1, fit%ncomp
            ! even half -> band-limited real image (unmasked, for an unbiased FSC)
            fit%Yeven(q)%rho_exp = 1.0
            call fit%Yeven(q)%compress_exp; call fit%Yeven(q)%ifft
            call realvols(q)%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call fit%Yeven(q)%get_rmat_ptr(rmatp); call realvols(q)%set_rmat(rmatp, .false.)
            call realvols(q)%mul(mstep_gridcorr)
            ! odd half
            fit%Yodd(q)%rho_exp = 1.0
            call fit%Yodd(q)%compress_exp; call fit%Yodd(q)%ifft
            call img_o%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call fit%Yodd(q)%get_rmat_ptr(rmatp); call img_o%set_rmat(rmatp, .false.)
            call img_o%mul(mstep_gridcorr)
            ! half-set FSC -> Wiener filter 2F/(1+F)
            call realvols(q)%fft; call img_o%fft
            call realvols(q)%fsc(img_o, corrs)
            fscq_dg(:,q) = corrs
            do sh = 1, filtsz
                fc = max(0., min(0.999, corrs(sh)))
                filt(sh) = 2.*fc/(1.+fc)
            end do
            ! merged, FSC-Wiener + band-limit filtered, back to real, masked
            call realvols(q)%add(img_o); call realvols(q)%mul(0.5)
            ! SIMPLE_COV_MSTEP_WIENER=0: ship the unshrunk merged basis (band limit and window only). The
            ! half-set FSC Wiener pins each component at its own cross-half resolution (~20 A on 16.8k PfCRT
            ! particles), below the scale of a helix displacement; without it the basis keeps every shell
            ! to the band, noisier but not blurred, and the state reconstruction from all particles decides
            if( l_mstep_wiener ) call realvols(q)%apply_filter(filt)
            if( fit%lp_it > 2.0*params%smpd_crop + TINY ) call realvols(q)%bp(0., fit%lp_it)
            call realvols(q)%ifft
            call flex_window_apply(realvols(q), params)
            ! stash both halves, filtered exactly as the merged basis is
            call eimgs(q)%copy(realvols(q))
            call img_o%ifft
            if( fit%lp_it > 2.0*params%smpd_crop + TINY )then
                call img_o%fft; call img_o%bp(0., fit%lp_it); call img_o%ifft
            endif
            call flex_window_apply(img_o, params)
            call oimgs(q)%copy(img_o)
            call img_o%kill
        end do
        call mstep_gridcorr%kill
        ! Band/rank diagnostic (report only): heterogeneity resolution, where the leading
        ! components' split-half FSC crosses 0.143, and the count of components with significant
        ! in-band FSC.
        khi_dg = filtsz
        if( params%lp > 2.0*params%smpd_crop + TINY ) &
            &khi_dg = max(2, min(filtsz, int(fit%dstep_ann/params%lp)))
        fit%khi_fit = khi_dg   ! per-fit band record
        nsig_dg = 0
        do q = 1, fit%ncomp
            fmean_dg(q) = sum(fscq_dg(1:khi_dg,q))/real(khi_dg)
            if( fmean_dg(q) > 0.143 ) nsig_dg = nsig_dg + 1
        end do
        ! heterogeneity resolution off the strongest components, ranked by in-band mean FSC:
        ! the EM basis is orthonormalised in fixed order and never variance-sorted
        ntop_dg = min(4, fit%ncomp)
        sel_dg  = 0
        do tq_dg = 1, ntop_dg
            fbest_dg = -huge(1.d0); bq_dg = 1
            do q = 1, fit%ncomp
                if( any(sel_dg(1:tq_dg-1) == q) ) cycle
                if( fmean_dg(q) > fbest_dg )then
                    fbest_dg = fmean_dg(q); bq_dg = q
                endif
            end do
            sel_dg(tq_dg) = bq_dg
        end do
        res_dg  = fit%dstep_ann/real(khi_dg)
        do sh = 2, filtsz
            fc = 0.
            do tq_dg = 1, ntop_dg
                fc = fc + fscq_dg(sh,sel_dg(tq_dg))
            end do
            fc = fc/real(ntop_dg)
            if( fc < 0.143 )then
                res_dg = fit%dstep_ann/real(sh)
                exit
            endif
        end do
        write(logfhandle,'(A,I0,A,F7.2,A,I0,A,I0,A)') '>>> FLEX_PCA BAND/RANK it=',it_eff, &
            &'  het-resolution(FSC 0.143)=',res_dg,' A   components with in-band FSC>0.143: ', &
            &nsig_dg,' of ',fit%ncomp,'   (auto-lp / auto-neigs candidates)'
        call flush(logfhandle)
        ! cross-FSC writer payloads: the internal e/o FSC curves and this iteration's Gamma,
        ! stashed for the driver's writer
        if( fit%l_xf_harvest )then
            if( allocated(fit%xf_fscq) ) deallocate(fit%xf_fscq)
            allocate(fit%xf_fscq(filtsz,fit%ncomp), source=fscq_dg)
            if( allocated(fit%xf_gam) ) deallocate(fit%xf_gam)
            allocate(fit%xf_gam(fit%ncomp), source=fit%gam_acc(1:fit%ncomp))
        endif
        deallocate(fscq_dg)
        deallocate(filt, corrs)
        ! Even/odd update agreement: both half-bases inherit the current basis, so comparing them
        ! directly saturates at once. Each half is instead projected onto the orthogonal complement
        ! of the previous basis, and the principal angles between those residual spans (sum =
        ! soft count of reproducible new directions) say whether the update carries signal.
        if( allocated(fit%prev_real) )then
            call deflate_against_basis(eimgs, fit%ncomp, fit%prev_real, size(fit%prev_real))
            call deflate_against_basis(oimgs, fit%ncomp, fit%prev_real, size(fit%prev_real))
            call cross_half_subspace_angles(eimgs, oimgs, fit%ncomp, sv_eo)
            eo_dim = sum(sv_eo)
        else
            allocate(sv_eo(fit%ncomp), source=0.d0)
            eo_dim = -1.d0            ! first iteration: no previous basis
        endif
        write(logfhandle,'(A,I0,A,F8.3,A,I0,A,I0)') '>>> FLEX_PCA PROBE ITER ',it_eff, &
            &'  update agreement (even|odd, prev basis deflated)=',eo_dim, &
            &'   n>=0.9: ',count(sv_eo >= 0.9d0),' of ',fit%ncomp
        call flush(logfhandle)
        if( eo_dim < 0.d0 )then
            continue                  ! nothing to judge yet
        else if( eo_dim > fit%eo_best + COV_EO_TOL )then
            fit%eo_best = eo_dim; fit%eo_stall = 0
        else
            fit%eo_stall = fit%eo_stall + 1
        endif
        deallocate(sv_eo)
        do q = 1, fit%ncomp
            call eimgs(q)%kill; call oimgs(q)%kill
        end do
        deallocate(eimgs, oimgs)
        ! Mean-shaped deflation (SIMPLE_COV_EM_DEFLATE=n): a_i is a scalar contrast fit, so any
        ! frequency-dependent per-image scale (envelope/B-factor spread) lands in the residual as
        ! a consensus-shaped term coherent across particles and would take a whole component.
        ! n=1 removes the mean direction; n>1 removes n resolution shells of the consensus, which
        ! covers any scale that varies smoothly with frequency. Soft masking breaks the Parseval
        ! orthogonality of disjoint bands, so the shells are explicitly orthonormalised.
        if( fit%l_deflate_mean )then
            ndfl = max(1, fit%vdfl)
            ! the flat-in-mask background template is part of the deflation set by default (the
            ! ice/background mode is not consensus-shaped, so the shells never remove it; the merge
            ! step already had it baked on, the per-fit M-step read an env switch that the unified
            ! recipe never exports); SIMPLE_COV_DEFLATE_BG=0 opts out
            l_dfl_bg = .not. cov_env_int_off('SIMPLE_COV_DEFLATE_BG')
            if( l_dfl_bg ) ndfl = ndfl + 1
            ndfl_sh = ndfl
            if( l_dfl_bg )   ndfl_sh = ndfl_sh - 1
            ! the dilation of the consensus, (x-c).grad rho, is the breathing mode (magnification
            ! and defocus scatter); deflated by default, SIMPLE_COV_DEFLATE_DILATION=0 opts out
            l_dfl_dil = .not. cov_env_int_off('SIMPLE_COV_DEFLATE_DILATION')
            if( l_dfl_dil ) ndfl = ndfl + 1
            ! complement block: each block deflates against templates carried on ITS window (a template
            ! windowed by the other block's support is orthogonal to this block's vectors anyway)
            nblk_dfl = 1
            do iblk_dfl = 1, nblk_dfl
            qlo_dfl = 1; qhi_dfl = fit%ncomp
            allocate(dfl_basis(ndfl))
            call mvol_dfl%read_and_crop(params%vols(1), params%smpd, params%box_crop, params%smpd_crop)
            kfr_dfl = covariance_kfromto(params)
            if( l_dfl_bg )then
                ! after the shells: uniform density under the soft spherical mask
                call dfl_basis(ndfl_sh+1)%copy(mvol_dfl)
                call dfl_basis(ndfl_sh+1)%get_rmat_ptr(rm_dfl)
                rm_dfl = 0.
                rm_dfl(1:params%box_crop, 1:params%box_crop, 1:params%box_crop) = 1.
                call flex_window_apply(dfl_basis(ndfl_sh+1), params)
                write(logfhandle,'(A)') '>>> FLEX_PCA DEFLATE_BG: flat-in-mask background template &
                    &appended to the deflation set'
                call flush(logfhandle)
            endif
            if( l_dfl_dil )then
                call dfl_basis(ndfl)%copy(mvol_dfl)
                call dilation_template(dfl_basis(ndfl), params%box_crop)
                call flex_window_apply(dfl_basis(ndfl), params)
                write(logfhandle,'(A)') '>>> FLEX_PCA DEFLATE_DILATION: consensus dilation (breathing) template &
                    &appended to the deflation set'
                call flush(logfhandle)
            endif
            do idfl = 1, ndfl_sh
                call dfl_basis(idfl)%copy(mvol_dfl)
                if( ndfl_sh > 1 )then
                    ! shell idfl of ndfl in resolution (dstep/shell, covariance_kfromto's
                    ! convention); the first shell is a pure low-pass so no high-pass edge at
                    ! k~1 shaves the shells carrying most of the consensus power
                    if( idfl == 1 )then
                        res_lo = 0.
                    else
                        res_lo = fit%dstep_ann / &
                            &max(1., real(kfr_dfl(1)) + real(idfl-1)*real(max(1,kfr_dfl(2)-kfr_dfl(1)))/real(ndfl_sh))
                    endif
                    res_hi = fit%dstep_ann / &
                        &max(1., real(kfr_dfl(1)) + real(idfl)  *real(max(1,kfr_dfl(2)-kfr_dfl(1)))/real(ndfl_sh))
                    call dfl_basis(idfl)%fft
                    ! width=1: bp's default cosine edge is 10 shells per side, wider than a
                    ! shell of the fitted band
                    call dfl_basis(idfl)%bp(res_lo, res_hi, width=1.0)
                    call dfl_basis(idfl)%ifft
                    ! padding columns are not guaranteed zero after an ifft, and the inner
                    ! products below run over the whole padded array
                    call dfl_basis(idfl)%get_rmat_ptr(rm_dfl)
                    rm_dfl(params%box_crop+1:,:,:) = 0.
                    ! band-passing spreads density outside the particle and the deflated
                    ! volumes are soft-masked, so the shells must be too; not applied at ndfl=1
                    call flex_window_apply(dfl_basis(idfl), params)
                endif
            end do
            call mvol_dfl%get_rmat_ptr(rv_dfl)
            mnorm_dfl = sqrt(sum(real(rv_dfl,dp)**2))
            if( ndfl > 1 )then
                write(logfhandle,'(A)',advance='no') '>>> FLEX_PCA deflation shell norms / |consensus|:'
                do idfl = 1, ndfl
                    call dfl_basis(idfl)%get_rmat_ptr(rm_dfl)
                    write(logfhandle,'(A,ES9.2)',advance='no') ' ', &
                        &sqrt(sum(real(rm_dfl,dp)**2))/max(mnorm_dfl,DTINY)
                end do
                write(logfhandle,*)
            endif
            ! modified Gram-Schmidt; a shell that the band edges or the mask emptied drops out
            nkeep_dfl = 0
            do idfl = 1, ndfl
                call dfl_basis(idfl)%get_rmat_ptr(rm_dfl)
                do jdfl = 1, nkeep_dfl
                    call dfl_basis(jdfl)%get_rmat_ptr(rv_dfl)
                    mv_dfl = sum(real(rm_dfl,dp)*real(rv_dfl,dp))
                    rm_dfl = rm_dfl - real(mv_dfl)*rv_dfl
                end do
                mm_dfl = sqrt(sum(real(rm_dfl,dp)*real(rm_dfl,dp)))
                if( mm_dfl <= 1.d-3*mnorm_dfl ) cycle
                rm_dfl    = rm_dfl / real(mm_dfl)
                nkeep_dfl = nkeep_dfl + 1
                if( nkeep_dfl /= idfl ) call dfl_basis(nkeep_dfl)%copy(dfl_basis(idfl))
            end do
            if( nkeep_dfl > 0 )then
                rem_dfl = 0.d0; tot_dfl = 0.d0
                do q = qlo_dfl, qhi_dfl
                    call realvols(q)%get_rmat_ptr(rv_dfl)
                    tot_dfl = tot_dfl + sum(real(rv_dfl,dp)**2)
                    do idfl = 1, nkeep_dfl
                        call dfl_basis(idfl)%get_rmat_ptr(rm_dfl)
                        mv_dfl  = sum(real(rm_dfl,dp)*real(rv_dfl,dp))   ! unit-norm already
                        rem_dfl = rem_dfl + mv_dfl*mv_dfl
                        rv_dfl  = rv_dfl - real(mv_dfl)*rm_dfl
                    end do
                end do
                write(logfhandle,'(A,I0,A,I0,A,F6.2,A)') '>>> FLEX_PCA PROBE ITER ',it_eff, &
                    &'  mean-shaped deflation (rank ',nkeep_dfl,') removed ', &
                    &100.d0*rem_dfl/max(tot_dfl,DTINY),' % of basis energy'
                call flush(logfhandle)
            endif
            do idfl = 1, ndfl
                call dfl_basis(idfl)%kill
            end do
            deallocate(dfl_basis)
            call mvol_dfl%kill
            end do   ! iblk_dfl
        endif
        ! orthonormalize the probe volumes -> refined basis
        call orthonormalize_representatives(params, build, realvols, fit%ncomp, utilde, utilde_real, d_new)
        ! Convergence: principal angles between successive bases. Both come out of
        ! orthonormalize_representatives orthonormal, so the cross-Gram singular values are
        ! principal-angle cosines (align_basis_to_reference's per-vector normalisation is a no-op).
        fit%l_converged = .false.
        if( allocated(fit%prev_real) )then
            call align_basis_to_reference(fit%prev_real, size(fit%prev_real), utilde_real, d_new, Mconv, sconv)
            cos_mean = sum(sconv) / real(max(1,size(sconv)),dp)
            write(logfhandle,'(A,I0,A,F9.6)') '>>> FLEX_PCA PROBE ITER ',it_eff, &
                &'  mean principal-angle cosine vs previous basis=',cos_mean
            call flush(logfhandle)
            ! reported only: the mean over all components is dominated by the non-reproducing
            ! tail and fires early; it is the criterion only when the even/odd signal is unavailable
            if( fit%ncomp < 2 .and. cos_mean >= fit%conv_thresh ) fit%l_converged = .true.
            deallocate(Mconv, sconv)
            do q = 1, size(fit%prev_real)
                call fit%prev_real(q)%kill
            end do
            deallocate(fit%prev_real)
        endif
        allocate(fit%prev_real(d_new))
        do q = 1, d_new
            call fit%prev_real(q)%copy(utilde_real(q))
        end do
        ! replace basis_recs with the refined (projection-ready) basis; eigvals = latent variances
        do q = 1, size(fit%basis_recs)
            call fit%basis_recs(q)%dealloc_rho; call fit%basis_recs(q)%kill
        end do
        deallocate(fit%basis_recs); allocate(fit%basis_recs(d_new))
        if( allocated(fit%eigvals) ) deallocate(fit%eigvals); allocate(fit%eigvals(d_new))
        do q = 1, d_new
            ! projection-ready basis reconstructor from the clean real basis image (mean_rec idiom)
            call init_basis_reconstructor(params, build, fit%basis_recs(q))
            call fit%basis_recs(q)%set_rmat(utilde_real(q)%get_rmat(), .false.)
            call fit%basis_recs(q)%fft
            call fit%basis_recs(q)%expand_exp
            ! EM Gamma update: the posterior second moment (1/n) sum_i (z_iq^2 + [A_i^-1]_qq); the
            ! MAP point-estimate variance underestimates Gamma and collapses the prior
            fit%eigvals(q) = max(fit%gam_acc(min(q,fit%ncomp)), DTINY)
            ! overwrite the eigenvolume MRC with the refined basis vector, in the fit's namespace
            fname = fit%fprefix//int2str_pad(q,3)//MRC_EXT
            call utilde_real(q)%write(fname, del_if_exists=.true.); call fname%kill
        end do
        ! Gamma is the EM update and is always valid: it is reduced from the per-particle posterior
        ! second moments, on either execution path. The point-estimate columns are NOT: the master of
        ! a distributed round holds no z of its own (the workers do, and ship sufficient statistics),
        ! so they are reported only when this process actually solved the particles.
        ! z is zero-initialised at allocation and only ever written by a process that solved the
        ! particles itself, so an all-zero z IS the distributed master's 'I do not hold it'
        l_z_local = any(fit%z /= 0.d0)
        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PROBE EM Gamma update (n=',fit%nval,' valid particles):'
        do q = 1, min(fit%ncomp,10)
            if( l_z_local )then
                mu_q = sum(fit%z(:,q)) / real(max(1,fit%npp),dp)
                sd_q = sum(fit%z(:,q)**2) / real(max(1,fit%nval),dp)
                if( ieee_is_finite(mu_q) .and. ieee_is_finite(sd_q) )then
                    write(logfhandle,'(A,I3,A,ES11.3,A,ES11.3,A,F7.3,A,ES10.2)') '>>>   z',q, &
                        &'  <z^2>=',sd_q,'  Gamma=',fit%gam_acc(q), &
                        &'  posterior_frac=',real(1.d0 - sd_q/max(fit%gam_acc(q),DTINY)), &
                        &'  mean(z)=',mu_q
                else
                    write(logfhandle,'(A,I3,A,ES11.3,A)') '>>>   z',q,'  Gamma=',fit%gam_acc(q), &
                        &'  (point-estimate moments non-finite; the latent solve did not stay bounded)'
                endif
            else
                write(logfhandle,'(A,I3,A,ES11.3,A)') '>>>   z',q,'  Gamma=',fit%gam_acc(q), &
                    &'  (distributed master: the point-estimate moments live on the workers)'
            endif
        end do
        call flush(logfhandle)
        fit%ncomp = d_new
        do ithr = 1, nthr
            call cleanup_plane(fit%mean_fpl(ithr))
            do q = 1, size(fit%basis_fpls,1); call cleanup_plane(fit%basis_fpls(q,ithr)); end do
        end do
        do q = 1, size(fit%Yeven); call fit%Yeven(q)%dealloc_rho; call fit%Yeven(q)%kill; call fit%Yodd(q)%dealloc_rho; call fit%Yodd(q)%kill; end do
        do q = 1, size(utilde); call utilde(q)%dealloc_rho; call utilde(q)%kill; call utilde_real(q)%kill; end do
        do q = 1, size(realvols); call realvols(q)%kill; end do
        deallocate(utilde, utilde_real, realvols)
        deallocate(fit%Yeven, fit%Yodd, fit%rho_e, fit%rho_o, fit%prior)
        deallocate(fit%Gth, fit%Ath, fit%bth, fit%cth, fit%zth, fit%basis_fpls, fit%mean_fpl, fit%zbatch, fit%dens, fit%valid, fit%valid_e, fit%valid_o)
        deallocate(fit%Ainvth, fit%Acpth, fit%gam_thr, fit%gam_acc, fit%nval_thr, fit%hth, fit%nll_thr)
        ! The mixture state (mix_*) and its work arrays are not freed here: they must survive
        ! from one iteration's M-step to the next E-step, and are freed only by kill_probe_fit
        ! (the resize block at iteration start handles dimension changes).
        if( allocated(fit%gam_dbg) ) deallocate(fit%gam_dbg)
        ! the update agreement is a diagnostic, not a stopping rule: it decays from the first
        ! iteration while the basis keeps improving (non-reproducible update directions can
        ! still carry a reproducible bias)
        if( .false. .and. fit%eo_stall >= fit%eo_patience ) fit%l_converged = .true.

      contains

        subroutine log_pcg_outcome( half, res )
            character(len=*),         intent(in) :: half
            type(flex_pcg_outcome_t), intent(in) :: res
            write(logfhandle,'(A,I0,A,A,A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,A,A,ES10.3,A,ES10.3,A,F7.4,A,ES10.3,A,F8.1)') &
                &'>>> FLEX_PCA PCG MSTEP it=', it_eff, ' ', trim(half), '  iters=', res%iteration_count, &
                &'  init=', res%initial_rel_residual, '  resid=', res%final_rel_residual, '  update=', &
                &res%final_rel_update, '  stop=', trim(res%stop_reason), '  |b|=', res%rhs_norm, '  |x0|=', &
                &res%start_norm, '  corr(b,Bx0)=', res%start_corr, '  scale=', res%start_scale, '  seconds=', res%seconds
            call flush(logfhandle)
        end subroutine log_pcg_outcome

    end subroutine fit_iter_finish

end submodule simple_flex_pca_em_mstep
