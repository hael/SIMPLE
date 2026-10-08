!@descr: flex_pca probe-fit iteration setup and M-step update procedures
submodule (simple_flex_probe_fit) simple_flex_probe_fit_update
use simple_core_module_api, only: dtiny, fdim, get_resarr, get_resolution_at_fsc, int2str_pad, logfhandle, &
    &maximgbatchsz, mrc_ext, tiny
use simple_defs_flex,          only: FLEX_FSC_SIGNAL_THRESHOLD
use simple_image,              only: image
use simple_gridding,           only: prep3D_inv_kbenvelope4mul
use simple_flex_pca_crossfsc,  only: crossfsc_harvest_h
use simple_flex_pca_mstep,     only: init_basis_reconstructor
use simple_flex_pca_basis,     only: covariance_kfromto, orthonormalize_representatives, &
    &align_basis_to_reference, deflate_against_basis, cross_half_subspace_angles
use simple_flex_pca_util,      only: dilation_template
use simple_flex_pca_fit_types, only: cleanup_plane
use , intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
#include "simple_local_flags.inc"

contains

    !> Master tail of one EM iteration for one fit: Gamma and likelihood, cross-FSC ridge (when a record
    !! exists), coupled per-half solve, FSC-Wiener merge, mean-shaped deflation, re-orthonormalisation and
    !! basis swap. Frees per-iteration fields only; mix_* persist (see probe_fit_t).
    module subroutine fit_iter_finish( params, build, fit, it_eff, nthr )
        class(parameters),     intent(inout) :: params
        type(builder),         intent(inout) :: build
        class(flex_probe_fit), intent(inout) :: fit
        integer,               intent(in)    :: it_eff, nthr
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
        real,     allocatable :: fcurve_dg(:), resarr_dg(:)
        integer  :: q, ithr, sh, filtsz, d_new
        integer  :: ndfl, ndfl_sh, idfl, jdfl, nkeep_dfl, kfr_dfl(2)
        integer  :: khi_dg, nsig_dg, ntop_dg, sel_dg(4), tq_dg, bq_dg
        logical  :: l_z_local
        integer  :: iblk_dfl, nblk_dfl, qlo_dfl, qhi_dfl
        do q = 1, fit%model%ncomp
            fit%iter%gam_acc(q) = fit%iter%gam_sum(q) / real(max(1,fit%iter%nval),dp)
        end do
        deallocate(fit%iter%gam_sum)
        ! Marginal likelihood (the EM objective):
        !   -2 log p(y_i) = ||y_i - a_i T_i u_0||^2/sig2 - h_i' A_i^-1 h_i + log det A_i + log det Gamma
        ! The first term is fixed across iterations (a_i and mu are held), so only the rest is
        ! accumulated; log det Gamma uses `prior`, the Gamma the E-step just used. The filtered
        ! M-step sits outside the likelihood, so this is a diagnostic, not a monotonicity guarantee.
        if( fit%history%l_mix_used )then
            ! the mixture marginal's Omega log det, with the Omega this iteration's E-step used
            fit%iter%nll_tot = fit%iter%nll_tot + real(fit%iter%nval,dp)*fit%history%ldOm_used
        else
            fit%iter%nll_tot = fit%iter%nll_tot - real(fit%iter%nval,dp)*sum(log(max(fit%iter%prior(1:fit%model%ncomp), DTINY)))
        endif
        write(logfhandle,'(A,I0,A,ES14.6,A,ES12.4)') '>>> FLEX_PCA PROBE ITER ',it_eff, &
            &'  -2logL/N (varying part)=',fit%iter%nll_tot/real(max(1,fit%iter%nval),dp), &
            &'  delta=',(fit%iter%nll_tot/real(max(1,fit%iter%nval),dp)) - fit%history%nll_prev
        fit%history%nll_prev = fit%iter%nll_tot/real(max(1,fit%iter%nval),dp)
        call flush(logfhandle)
        ! Cross-FSC sampling profiles H: per-shell means of the rho_e/rho_o diagonal pair rows,
        ! harvested before any invtau2 mutates them. On the single-fit engine h_a/h_b are the
        ! internal even/odd halves; the paired writer sums e+o per fit.
        if( fit%diag%l_xf_harvest )then
            block
                integer :: xf_lb(3), xf_nyq, xf_fsz
                xf_lb  = lbound(fit%mstep%Yeven(1)%cmat_exp)
                xf_nyq = fit%mstep%Yeven(1)%get_lfny(1)
                xf_fsz = max(1, fdim(params%box_crop) - 1)
                if( allocated(fit%diag%xf_h_e) ) deallocate(fit%diag%xf_h_e)
                if( allocated(fit%diag%xf_h_o) ) deallocate(fit%diag%xf_h_o)
                allocate(fit%diag%xf_h_e(xf_fsz,fit%model%ncomp), fit%diag%xf_h_o(xf_fsz,fit%model%ncomp))
                call crossfsc_harvest_h(fit%mstep%rho_e, fit%mstep%npairs, fit%model%ncomp, xf_lb, xf_nyq, &
                    &xf_fsz, fit%diag%xf_h_e)
                call crossfsc_harvest_h(fit%mstep%rho_o, fit%mstep%npairs, fit%model%ncomp, xf_lb, xf_nyq, &
                    &xf_fsz, fit%diag%xf_h_o)
            end block
        endif
        ! Cross-FSC SSNR ridge (flex analog of add_invtausq2rho): invtau2 from record t-1, built by
        ! xfsc_prep_iter, added to the rho_e/rho_o diagonals right before the solves; one-shot.
        if( allocated(fit%diag%xf_invtau2) )then
            call fit%mstep%apply_ridge(fit%model%ncomp, fit%diag%xf_invtau2)
            deallocate(fit%diag%xf_invtau2)
        endif
        ! Coupled M-step solve: at every grid point the components share one k x k normal matrix
        ! sum_i |CTF|^2 E[z_i z_i'], so the basis volumes are solved together per halfset. A unit
        ! density is then handed to compress_exp (the divide has already happened).
        call fit%mstep%solve_halves(params, fit%model%ncomp, it_eff)
        ! finalize even/odd Y_q, half-set FSC Wiener-merge -> band-limited masked basis volume;
        ! the two half-bases are kept for the update-agreement diagnostic below
        allocate(realvols(fit%model%ncomp))
        allocate(eimgs(fit%model%ncomp), oimgs(fit%model%ncomp))
        filtsz = max(1, fdim(params%box_crop) - 1)
        allocate(filt(filtsz), corrs(filtsz))
        ! native-lattice inverse KB envelope, applied to each half before FSC/merge/mask
        mstep_gridcorr = prep3D_inv_kbenvelope4mul([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        allocate(fscq_dg(filtsz, fit%model%ncomp), source=0.)
        do q = 1, fit%model%ncomp
            ! even half -> band-limited real image (unmasked, for an unbiased FSC)
            fit%mstep%Yeven(q)%rho_exp = 1.0
            call fit%mstep%Yeven(q)%compress_exp; call fit%mstep%Yeven(q)%ifft
            call realvols(q)%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call fit%mstep%Yeven(q)%get_rmat_ptr(rmatp); call realvols(q)%set_rmat(rmatp, .false.)
            call realvols(q)%mul(mstep_gridcorr)
            ! odd half
            fit%mstep%Yodd(q)%rho_exp = 1.0
            call fit%mstep%Yodd(q)%compress_exp; call fit%mstep%Yodd(q)%ifft
            call img_o%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call fit%mstep%Yodd(q)%get_rmat_ptr(rmatp); call img_o%set_rmat(rmatp, .false.)
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
            call realvols(q)%apply_filter(filt)
            if( fit%spec%lp_it > 2.0*params%smpd_crop + TINY ) call realvols(q)%bp(0., fit%spec%lp_it)
            call realvols(q)%ifft
            call fit%mstep%env%apply(realvols(q), params%box_crop, params%msk_crop)
            ! stash both halves, filtered exactly as the merged basis is
            call eimgs(q)%copy(realvols(q))
            call img_o%ifft
            if( fit%spec%lp_it > 2.0*params%smpd_crop + TINY )then
                call img_o%fft; call img_o%bp(0., fit%spec%lp_it); call img_o%ifft
            endif
            call fit%mstep%env%apply(img_o, params%box_crop, params%msk_crop)
            call oimgs(q)%copy(img_o)
            call img_o%kill
        end do
        call mstep_gridcorr%kill
        ! Band/rank diagnostic (report only): heterogeneity resolution, where the leading
        ! components' split-half FSC crosses 0.143, and the count of components with significant
        ! in-band FSC.
        khi_dg = filtsz
        if( params%lp > 2.0*params%smpd_crop + TINY ) &
            &khi_dg = max(2, min(filtsz, int(fit%spec%dstep_ann/params%lp)))
        nsig_dg = 0
        do q = 1, fit%model%ncomp
            fmean_dg(q) = sum(fscq_dg(1:khi_dg,q))/real(khi_dg)
            if( fmean_dg(q) > real(FLEX_FSC_SIGNAL_THRESHOLD) ) nsig_dg = nsig_dg + 1
        end do
        ! heterogeneity resolution off the strongest components, ranked by in-band mean FSC:
        ! the EM basis is orthonormalised in fixed order and never variance-sorted
        ntop_dg = min(4, fit%model%ncomp)
        sel_dg  = 0
        do tq_dg = 1, ntop_dg
            fbest_dg = -huge(1.d0); bq_dg = 1
            do q = 1, fit%model%ncomp
                if( any(sel_dg(1:tq_dg-1) == q) ) cycle
                if( fmean_dg(q) > fbest_dg )then
                    fbest_dg = fmean_dg(q); bq_dg = q
                endif
            end do
            sel_dg(tq_dg) = bq_dg
        end do
        ! the resolution where the mean curve of the strongest components last passes 0.143
        allocate(fcurve_dg(filtsz), source=0.)
        do tq_dg = 1, ntop_dg
            fcurve_dg = fcurve_dg + fscq_dg(:,sel_dg(tq_dg))
        end do
        fcurve_dg = fcurve_dg/real(ntop_dg)
        resarr_dg = get_resarr(params%box_crop, params%smpd_crop)
        call get_resolution_at_fsc(fcurve_dg, resarr_dg, real(FLEX_FSC_SIGNAL_THRESHOLD), res_dg)
        deallocate(fcurve_dg, resarr_dg)
        write(logfhandle,'(A,I0,A,F7.2,A,I0,A,I0,A)') '>>> FLEX_PCA BAND/RANK it=',it_eff, &
            &'  het-resolution(FSC 0.143)=',res_dg,' A   components with in-band FSC>0.143: ', &
            &nsig_dg,' of ',fit%model%ncomp,'   (auto-lp / auto-neigs candidates)'
        call flush(logfhandle)
        deallocate(fscq_dg)
        deallocate(filt, corrs)
        ! Even/odd update agreement: both half-bases inherit the current basis, so comparing them
        ! directly saturates at once. Each half is instead projected onto the orthogonal complement
        ! of the previous basis, and the principal angles between those residual spans (sum =
        ! soft count of reproducible new directions) say whether the update carries signal.
        if( allocated(fit%history%prev_real) )then
            call deflate_against_basis(eimgs, fit%model%ncomp, fit%history%prev_real, size(fit%history%prev_real))
            call deflate_against_basis(oimgs, fit%model%ncomp, fit%history%prev_real, size(fit%history%prev_real))
            call cross_half_subspace_angles(eimgs, oimgs, fit%model%ncomp, sv_eo)
            eo_dim = sum(sv_eo)
        else
            allocate(sv_eo(fit%model%ncomp), source=0.d0)
            eo_dim = -1.d0            ! first iteration: no previous basis
        endif
        write(logfhandle,'(A,I0,A,F8.3,A,I0,A,I0)') '>>> FLEX_PCA PROBE ITER ',it_eff, &
            &'  update agreement (even|odd, prev basis deflated)=',eo_dim, &
            &'   n>=0.9: ',count(sv_eo >= 0.9d0),' of ',fit%model%ncomp
        call flush(logfhandle)
        deallocate(sv_eo)
        do q = 1, fit%model%ncomp
            call eimgs(q)%kill; call oimgs(q)%kill
        end do
        deallocate(eimgs, oimgs)
        ! Mean-shaped deflation: the scalar contrast a_i leaves any frequency-dependent per-image scale in the
        ! residual as a consensus-shaped term; vdfl resolution shells of the consensus (1 = the mean) remove it.
        if( fit%spec%l_deflate_mean )then
            ndfl_sh = max(1, fit%spec%vdfl)
            ! The flat-in-mask background and consensus dilation are part of the production
            ! deflation recipe in addition to the consensus shells.
            ndfl = ndfl_sh + 2
            ! a single block spanning every component (nblk_dfl = 1, qlo_dfl:qhi_dfl = 1:ncomp)
            nblk_dfl = 1
            do iblk_dfl = 1, nblk_dfl
                qlo_dfl = 1; qhi_dfl = fit%model%ncomp
                allocate(dfl_basis(ndfl))
                call mvol_dfl%read_and_crop(params%vols(1), params%smpd, params%box_crop, params%smpd_crop)
                kfr_dfl = covariance_kfromto(params)
                ! after the shells: uniform density under the soft spherical mask
                call dfl_basis(ndfl_sh+1)%copy(mvol_dfl)
                call dfl_basis(ndfl_sh+1)%get_rmat_ptr(rm_dfl)
                rm_dfl = 0.
                rm_dfl(1:params%box_crop, 1:params%box_crop, 1:params%box_crop) = 1.
                call fit%mstep%env%apply(dfl_basis(ndfl_sh+1), params%box_crop, params%msk_crop)
                write(logfhandle,'(A)') '>>> FLEX_PCA DEFLATE_BG: flat-in-mask background template &
                    &appended to the deflation set'
                call flush(logfhandle)
                call dfl_basis(ndfl)%copy(mvol_dfl)
                call dilation_template(dfl_basis(ndfl), params%box_crop)
                call fit%mstep%env%apply(dfl_basis(ndfl), params%box_crop, params%msk_crop)
                write(logfhandle,'(A)') '>>> FLEX_PCA DEFLATE_DILATION: consensus dilation (breathing) template &
                    &appended to the deflation set'
                call flush(logfhandle)
                do idfl = 1, ndfl_sh
                    call dfl_basis(idfl)%copy(mvol_dfl)
                    if( ndfl_sh > 1 )then
                        ! shell idfl of ndfl in resolution (dstep/shell, covariance_kfromto's
                        ! convention); the first shell is a pure low-pass so no high-pass edge at
                        ! k~1 shaves the shells carrying most of the consensus power
                        if( idfl == 1 )then
                            res_lo = 0.
                        else
                            res_lo = fit%spec%dstep_ann / &
                                &max(1., real(kfr_dfl(1)) + real(idfl-1)*real(max(1,kfr_dfl(2)-kfr_dfl(1)))/real(ndfl_sh))
                        endif
                        res_hi = fit%spec%dstep_ann / &
                            &max(1., real(kfr_dfl(1)) + real(idfl)  *real(max(1,kfr_dfl(2)-kfr_dfl(1)))/real(ndfl_sh))
                        call dfl_basis(idfl)%fft
                        ! width=1: bp's default cosine edge is 10 shells per side, wider than a
                        ! shell of the fitted band
                        call dfl_basis(idfl)%bp(res_lo, res_hi, width=1.0)
                        call dfl_basis(idfl)%ifft
                        ! band-passing spreads density outside the particle and the deflated
                        ! volumes are soft-masked, so the shells must be too; not applied at ndfl=1
                        call fit%mstep%env%apply(dfl_basis(idfl), params%box_crop, params%msk_crop)
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
        call orthonormalize_representatives(params, build, realvols, fit%model%ncomp, utilde, utilde_real, d_new)
        ! Convergence: principal angles between successive bases. Both come out of
        ! orthonormalize_representatives orthonormal, so the cross-Gram singular values are
        ! principal-angle cosines (align_basis_to_reference's per-vector normalisation is a no-op).
        fit%history%l_converged = .false.
        if( allocated(fit%history%prev_real) )then
            call align_basis_to_reference(fit%history%prev_real, size(fit%history%prev_real), utilde_real, d_new, Mconv, sconv)
            cos_mean = sum(sconv) / real(max(1,size(sconv)),dp)
            write(logfhandle,'(A,I0,A,F9.6)') '>>> FLEX_PCA PROBE ITER ',it_eff, &
                &'  mean principal-angle cosine vs previous basis=',cos_mean
            call flush(logfhandle)
            ! a stopping rule for a rank-1 fit only: over several components the non-reproducing tail
            ! dominates the mean and it fires early
            if( fit%model%ncomp < 2 .and. cos_mean >= fit%spec%conv_thresh ) fit%history%l_converged = .true.
            deallocate(Mconv, sconv)
            do q = 1, size(fit%history%prev_real)
                call fit%history%prev_real(q)%kill
            end do
            deallocate(fit%history%prev_real)
        endif
        allocate(fit%history%prev_real(d_new))
        do q = 1, d_new
            call fit%history%prev_real(q)%copy(utilde_real(q))
        end do
        ! replace basis_recs with the refined (projection-ready) basis; eigvals = latent variances
        do q = 1, size(fit%model%basis_recs)
            call fit%model%basis_recs(q)%dealloc_rho; call fit%model%basis_recs(q)%kill
        end do
        deallocate(fit%model%basis_recs); allocate(fit%model%basis_recs(d_new))
        if( allocated(fit%model%eigvals) ) deallocate(fit%model%eigvals); allocate(fit%model%eigvals(d_new))
        do q = 1, d_new
            ! projection-ready basis reconstructor from the clean real basis image (mean_rec idiom)
            call init_basis_reconstructor(params, build, fit%model%basis_recs(q))
            call fit%model%basis_recs(q)%set_rmat(utilde_real(q)%get_rmat(), .false.)
            call fit%model%basis_recs(q)%fft
            call fit%model%basis_recs(q)%expand_exp
            ! EM Gamma update: the posterior second moment (1/n) sum_i (z_iq^2 + [A_i^-1]_qq); the
            ! MAP point-estimate variance underestimates Gamma and collapses the prior
            fit%model%eigvals(q) = max(fit%iter%gam_acc(min(q,fit%model%ncomp)), DTINY)
            ! overwrite the eigenvolume MRC with the refined basis vector, in the fit's namespace
            fname = fit%spec%fprefix//int2str_pad(q,3)//MRC_EXT
            call utilde_real(q)%write(fname, del_if_exists=.true.); call fname%kill
        end do
        ! Gamma comes from reduced posterior moments on either path; the z moments only where z was solved.
        ! z is zeroed at allocation and written only by a solving process: all-zero = distributed master.
        l_z_local = any(fit%model%z /= 0.d0)
        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PROBE EM Gamma update (n=',fit%iter%nval,' valid particles):'
        do q = 1, min(fit%model%ncomp,10)
            if( l_z_local )then
                mu_q = sum(fit%model%z(:,q)) / real(max(1,fit%spec%npp),dp)
                sd_q = sum(fit%model%z(:,q)**2) / real(max(1,fit%iter%nval),dp)
                if( ieee_is_finite(mu_q) .and. ieee_is_finite(sd_q) )then
                    write(logfhandle,'(A,I3,A,ES11.3,A,ES11.3,A,F7.3,A,ES10.2)') '>>>   z',q, &
                        &'  <z^2>=',sd_q,'  Gamma=',fit%iter%gam_acc(q), &
                        &'  posterior_frac=',real(1.d0 - sd_q/max(fit%iter%gam_acc(q),DTINY)), &
                        &'  mean(z)=',mu_q
                else
                    write(logfhandle,'(A,I3,A,ES11.3,A)') '>>>   z',q,'  Gamma=',fit%iter%gam_acc(q), &
                        &'  (point-estimate moments non-finite; the latent solve did not stay bounded)'
                endif
            else
                write(logfhandle,'(A,I3,A,ES11.3,A)') '>>>   z',q,'  Gamma=',fit%iter%gam_acc(q), &
                    &'  (distributed master: the point-estimate moments live on the workers)'
            endif
        end do
        call flush(logfhandle)
        fit%model%ncomp = d_new
        do ithr = 1, nthr
            call cleanup_plane(fit%iter%mean_fpl(ithr))
            call fit%iter%mean_es(ithr)%kill
            do q = 1, size(fit%iter%basis_fpls,1); call cleanup_plane(fit%iter%basis_fpls(q,ithr)); end do
        end do
        call fit%mstep%kill_iteration
        do q = 1, size(utilde); call utilde(q)%dealloc_rho; call utilde(q)%kill; call utilde_real(q)%kill; end do
        do q = 1, size(realvols); call realvols(q)%kill; end do
        deallocate(utilde, utilde_real, realvols)
        deallocate(fit%iter%prior)
        deallocate(fit%iter%Gth, fit%iter%Ath, fit%iter%bth, fit%iter%cth, fit%iter%zth, fit%iter%basis_fpls, fit%iter%mean_fpl, fit%iter%mean_es, fit%iter%zbatch, fit%iter%dens, fit%iter%valid, fit%iter%valid_e, fit%iter%valid_o)
        deallocate(fit%iter%Ainvth, fit%iter%Acpth, fit%iter%gam_thr, fit%iter%gam_acc, fit%iter%nval_thr, fit%iter%hth, fit%iter%nll_thr)
        ! The mixture state (mix_*) and its work arrays are not freed here: they must survive
        ! from one iteration's M-step to the next E-step, and are freed only by kill_probe_fit
        ! (the resize block at iteration start handles dimension changes).
        if( allocated(fit%diag%gam_dbg) ) deallocate(fit%diag%gam_dbg)

    end subroutine fit_iter_finish

    !> Stage begin: the E-step formulation and every per-fit policy value, chosen once per stage
    !! and never inside a batch loop. Under the paired driver this runs once per fit; the values
    !! are run settings, so the two fits agree by construction.
    module subroutine fit_estep_begin_stage( fit, params, nthr )
        class(flex_probe_fit), intent(inout) :: fit
        class(parameters),     intent(inout) :: params
        integer,               intent(in)    :: nthr
        ! ---- POLAR E-STEP: in-batch-loop shared-direction bank (the E-step) ----
        fit%estep%l_pol_hyb  = .false.
        fit%estep%rhyb_es    = 0
        fit%estep%npos_es    = 0
        fit%estep%l_pol_grid    = .false.
        fit%estep%l_pol_bank_it = .false.
        fit%estep%sec_bank      = 0.
        write(logfhandle,'(A)') '>>> FLEX_PCA POLAR E-STEP ON (stage 1): shared-direction bank &
            &supplies G/b/c; data-plane prep and M-step insertion stay Cartesian'
        call flush(logfhandle)
        ! MCFA: one basis, K latent Gaussians (Baek et al., IEEE TPAMI 2010); K=1 is plain PPCA EM.
        fit%spec%kmix   = COV_EM_MIX
        fit%spec%l_mix_req = fit%spec%kmix >= 1
        ! the MCFA accumulators and a bounded latent subsample (mcfa_init's seed) ride in the probe part files
        fit%spec%n_mix_warm = 1
        fit%history%l_mix_active = .false.
        fit%history%l_mix_used   = .false.
        fit%history%ldOm_mix     = 0.d0
        fit%history%ldOm_used    = 0.d0
        if( fit%spec%l_mix_req ) write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA MIX requested: K=', &
            &fit%spec%kmix,'  (plain warm-up: ',fit%spec%n_mix_warm,' iterations)'
        if( .not. allocated(fit%diag%sec_ring_thr) )then
            allocate(fit%diag%sec_bank_thr(nthr), fit%diag%sec_ring_thr(nthr), &
                &fit%diag%sec_exact_thr(nthr), fit%diag%sec_solve_thr(nthr), source=0.d0)
        endif
        fit%history%nll_prev    = 0.d0
        ! mean-shaped deflation: the per-particle scalar contrast cannot absorb a frequency-dependent
        ! scale, which would otherwise take a whole consensus-shaped component
        fit%spec%vdfl           = COV_EM_DEFLATE
        fit%spec%l_deflate_mean = .true.
        ! Per-particle contrast is the clamped mean-only estimate and stays fixed in the posterior.
        fit%spec%kfr_ann   = covariance_kfromto(params)
        fit%spec%dstep_ann = real(max(1, params%box_crop - 1)) * params%smpd_crop
        fit%spec%conv_thresh = COV_PROBE_CONV
    end subroutine fit_estep_begin_stage

    !> Per-iteration begin: band/lp schedule, prior from the current Gamma, MCFA freeze +
    !! rank-change resize, the per-iteration accumulators and thread scratch, and the
    !! projection-ready expansion of the fit's mean + basis.
    module subroutine fit_iter_begin( params, build, fit, mean_rec, it_eff, niters_eff, nthr )
        class(parameters),     intent(inout) :: params
        type(builder),         intent(inout) :: build
        class(flex_probe_fit), intent(inout) :: fit
        type(reconstructor),   intent(inout) :: mean_rec
        integer,               intent(in)    :: it_eff, niters_eff, nthr
        integer :: q
        fit%spec%lp_it  = params%lp
        fit%diag%sec_bank_thr  = 0.d0; fit%diag%sec_ring_thr = 0.d0
        fit%diag%sec_exact_thr = 0.d0; fit%diag%sec_solve_thr = 0.d0
        fit%estep%sec_bank = 0.; fit%estep%l_pol_bank_it = .false.   ! the polar E-step bank is rebuilt every iteration
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA PROBE SUBSPACE ITERATION ',it_eff,' / ',niters_eff, &
            &'  basis dim=',fit%model%ncomp
        call flush(logfhandle)
        allocate(fit%iter%prior(fit%model%ncomp))
        do q = 1, fit%model%ncomp
            fit%iter%prior(q) = 1.d0 / max(fit%model%eigvals(q), DTINY)
        end do
        ! MCFA bookkeeping. l_mix_used freezes which prior THIS iteration's E-step runs
        ! under (the M-step hook below flips l_mix_active mid-iteration); a basis dimension
        ! change invalidates xi/Omega, so the mixture re-warms and re-initialises.
        fit%history%l_mix_used = fit%history%l_mix_active
        fit%history%ldOm_used  = fit%history%ldOm_mix
        if( fit%spec%l_mix_req )then
            if( allocated(fit%history%rhs0th) )then
                if( size(fit%history%rhs0th,1) /= fit%model%ncomp )then
                    deallocate(fit%history%rhs0th, fit%history%mkth, fit%history%lwth, fit%history%rkth, fit%history%mxa_sr, fit%history%mxa_sm, fit%history%mxa_smm, fit%history%mxa_sainv)
                    if( allocated(fit%history%mix_xi) ) deallocate(fit%history%mix_xi, fit%history%mix_Om, fit%history%mix_Ominv, fit%history%mix_pi, &
                        &fit%history%mix_Omxi, fit%history%mix_xiOx, fit%history%mix_lpi)
                    if( fit%history%l_mix_active ) write(logfhandle,'(A)') &
                        &'>>> FLEX_PCA MIX re-initialising: basis dimension changed'
                    fit%history%l_mix_active = .false.
                    fit%history%l_mix_used   = .false.
                endif
            endif
            if( .not. allocated(fit%history%rhs0th) )then
                allocate(fit%history%rhs0th(fit%model%ncomp,nthr), fit%history%mkth(fit%model%ncomp,fit%spec%kmix,nthr), fit%history%lwth(fit%spec%kmix,nthr), &
                    &fit%history%rkth(fit%spec%kmix,nthr), fit%history%mxa_sr(fit%spec%kmix,nthr), fit%history%mxa_sm(fit%model%ncomp,fit%spec%kmix,nthr), &
                    &fit%history%mxa_smm(fit%model%ncomp,fit%model%ncomp,fit%spec%kmix,nthr), fit%history%mxa_sainv(fit%model%ncomp,fit%model%ncomp,nthr))
            endif
            fit%history%mxa_sr = 0.d0; fit%history%mxa_sm = 0.d0; fit%history%mxa_smm = 0.d0; fit%history%mxa_sainv = 0.d0
        endif
        ! the M-step system at this rank: even/odd numerators, the coupled normal matrix and,
        ! under rec_backend=pcg, the solver's lattice, support and packed kernel sums
        call fit%mstep%begin_iteration(params, build, fit%model%ncomp)
        allocate(fit%iter%Gth(fit%model%ncomp,fit%model%ncomp,nthr), fit%iter%Ath(fit%model%ncomp,fit%model%ncomp,nthr), fit%iter%bth(fit%model%ncomp,nthr), fit%iter%cth(fit%model%ncomp,nthr), fit%iter%zth(fit%model%ncomp,nthr))
        allocate(fit%iter%Ainvth(fit%model%ncomp,fit%model%ncomp,nthr), fit%iter%Acpth(fit%model%ncomp,fit%model%ncomp,nthr))
        allocate(fit%iter%basis_fpls(fit%model%ncomp,nthr), fit%iter%mean_fpl(nthr), fit%iter%mean_es(nthr))
        allocate(fit%iter%zbatch(fit%model%ncomp,MAXIMGBATCHSZ), fit%iter%dens(fit%model%ncomp,fit%model%ncomp,MAXIMGBATCHSZ))
        allocate(fit%iter%valid(MAXIMGBATCHSZ), fit%iter%valid_e(MAXIMGBATCHSZ), fit%iter%valid_o(MAXIMGBATCHSZ))
        allocate(fit%iter%gam_thr(fit%model%ncomp,nthr), source=0.d0)
        allocate(fit%iter%gam_acc(fit%model%ncomp), source=0.d0)
        allocate(fit%iter%nval_thr(nthr), source=0)
        allocate(fit%iter%hth(fit%model%ncomp,nthr), source=0.d0)
        allocate(fit%iter%nll_thr(nthr), source=0.d0)
        allocate(fit%diag%gam_dbg(4,nthr), source=0.d0)
        fit%iter%nll_tot = 0.d0
        fit%iter%dens = 0.d0
        ! the Gamma reduction target (filled by fit_iter_reduce or the distributed
        ! master's part reduce; freed by fit_iter_finish at the gam_acc step)
        allocate(fit%iter%gam_sum(fit%model%ncomp), source=0.d0)
        fit%iter%nval = 0
        call mean_rec%expand_exp
        do q = 1, fit%model%ncomp
            call fit%model%basis_recs(q)%expand_exp
        end do
    end subroutine fit_iter_begin

    !> Snapshot one fit's LAST-iteration raw M-step sufficient statistics + entry frame.
    !! Called by the paired driver after the batch loop and reductions, BEFORE fit_iter_finish
    !! (which ridges rho, mutates the numerators in the coupled solve, and frees everything).
    !! Overwrites every iteration: convergence is evaluated after the tails, so any iteration can
    !! turn out to be the last. prev_real at this point holds the
    !! PREVIOUS iteration's delivered basis == the frame the statistics' latents live in.
    module subroutine probe_fit_merge_stash( fit )
        class(flex_probe_fit), intent(inout) :: fit
        call fit%mstep%snapshot_for_pair_merge(fit%model%ncomp, fit%history%prev_real)
    end subroutine probe_fit_merge_stash

end submodule simple_flex_probe_fit_update
