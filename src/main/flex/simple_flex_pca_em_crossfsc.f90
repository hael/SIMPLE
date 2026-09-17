!@descr: flex_pca EM: the cross-fit-FSC driver context (artifact + SSNR ridge)
submodule (simple_flex_pca_em) simple_flex_pca_em_crossfsc
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

!> par.5.2 deepest-crossing criterion for the per-fit internal-FSC bands recorded in the artifact
real, parameter :: XFSC_CRIT = 0.143

contains



    !> Read the SIMPLE_COV_XFSC_REG arm (default 0: internal e/o Wiener only), reload the
    !! artifact when the ridge or the paired writer is live (restart-complete series -- the
    !! load_probe_state idiom) and fail fast on the contract violations. Master-only state:
    !! workers pass l_master=.false. and the whole subsystem stays inert on them (the ridge
    !! and the writer are master M-step / master-loop operations).
    module subroutine xfsc_setup( ctx, params, kfr_ann, l_paired, l_master )
        type(xfsc_ctx_t),  intent(inout) :: ctx
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: kfr_ann(2)
        logical,           intent(in)    :: l_paired, l_master
        ctx%v_reg      = 1     ! the cross-fit FSC ridge (the internal e/o arm was the scaffolding)
        ctx%l_paired   = l_paired
        ! the paired master writes honest paired=1 records EVERY iteration (spec par.2.2: the
        ! artifact is the engine's lasting product); the single-fit engine writes none
        ctx%l_writer   = l_paired .and. l_master
        ctx%pairing_id = 0
        if( l_paired )then
            ! one-time env read (hazard 7: never mid-loop); the driver already validated 1|3
            ctx%pairing_id = 1
            call cov_env_int('SIMPLE_COV_MOD4_PAIRING', ctx%pairing_id)
        endif
        ctx%l_any      = ctx%l_writer .or. ctx%v_reg > 0
        ctx%l_loaded   = .false.
        ctx%filtsz     = max(1, fdim(params%box_crop) - 1)
        ! crossfsc low-resolution exemption index: the reslim_ind analog (spec par.3.1 step 5)
        ctx%klo        = max(6, kfr_ann(1))
        if( .not. l_master )then
            ctx%l_any    = .false.
            ctx%l_writer = .false.
            return
        endif
        if( ctx%l_any )then
            ! restart-complete series reload (the load_probe_state idiom, applied to the artifact)
            call crossfsc_load(ctx%xf, ctx%l_loaded)
            if( ctx%l_loaded )then
                if( ctx%xf%box_crop /= params%box_crop .or. ctx%xf%filtsz /= ctx%filtsz ) &
                    &THROW_HARD('existing '//COV_XFSC_FNAME//' was written on a different lattice; &
                    &delete it to restart the series')
                if( ctx%l_writer .and. ctx%xf%paired /= 1 ) &
                    &THROW_HARD('refusing to append PAIRED records to a scaffolding (paired=0) &
                    &crossfsc artifact; move '//COV_XFSC_FNAME//' aside')
                write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA XFSC artifact reloaded: ',ctx%xf%nrec, &
                    &' record(s)'
                call flush(logfhandle)
            endif
            if( ctx%v_reg > 0 .and. .not. ctx%l_writer .and. .not. ctx%l_loaded )then
                write(logfhandle,'(A)') '>>> FLEX_PCA XFSC ridge requested but no artifact exists &
                    &and no writer is active; it stays inert (the paired engine writes the records)'
                call flush(logfhandle)
            endif
        endif
    end subroutine xfsc_setup


    !> Build one fit's ridge for THIS iteration (fit%xf_invtau2, applied by fit_iter_finish
    !! immediately before the coupled solves), from the latest record stamped <= t-1 -- the
    !! independence timing rule. Iteration 1 (no previous record) and any t-1 record that is
    !! scaffolding (paired=0) silently degrade to arm 0 for that iteration; the record's reg_mode
    !! field logs the degradation. Also sets l_xf_harvest when the writer needs the payloads.
    module subroutine xfsc_prep_iter( ctx, params, fit, it_eff, tag )
        type(xfsc_ctx_t),  intent(inout) :: ctx
        class(parameters), intent(in)    :: params
        type(probe_fit_t), intent(inout) :: fit
        integer,           intent(in)    :: it_eff
        character(len=*),  intent(in)    :: tag   !< '' single-fit; '  fit=A'/'  fit=B' paired
        logical :: l_fitb
        integer :: irec
        fit%l_xf_harvest = ctx%l_writer
        if( allocated(fit%xf_invtau2) ) deallocate(fit%xf_invtau2)
        ctx%reg_active = 0
        if( ctx%v_reg <= 0 ) return
        irec = crossfsc_latest_upto(ctx%xf, it_eff - 1)
        if( irec < 1 )then
            write(logfhandle,'(A,I0,A,A)') '>>> FLEX_PCA XFSC REG it=',it_eff, &
                &'  arm degraded to 0: no record stamped <= t-1',tag
        else if( ctx%xf%paired /= 1 )then
            write(logfhandle,'(A,I0,A,A)') '>>> FLEX_PCA XFSC REG it=',it_eff, &
                &'  arm degraded to 0: previous record is scaffolding (paired=0)',tag
        else
            ! which side of the record is THIS fit: the paired engine's fit%id names it
            l_fitb = fit%id == 2
            allocate(fit%xf_invtau2(fit%ncomp,ctx%filtsz), source=0.0)
            call xfsc_build_invtau2(ctx, params, ctx%xf%recs(irec), ctx%v_reg, l_fitb, it_eff, &
                &fit%xf_invtau2)
            ctx%reg_active = ctx%v_reg
            ! fit-level invtau2, added once per halfset by fit_iter_finish (the gold-standard
            ! precedent adds the shared curve into both halves' rho every iteration)
            write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA XFSC REG it=',it_eff, &
                &'  arm ',ctx%v_reg,' ridge added from record it=',ctx%xf%recs(irec)%it_eff,tag
        endif
        call flush(logfhandle)
    end subroutine xfsc_prep_iter

    !> Assemble the per-component, per-shell inverse prior variances for the crossfsc ridge
    !! (spec par.4.1) from one PAIRED (paired=1) record, mapped into THIS fit's component
    !! indices via the record's signed permutation. Arm 1: the matched cross-fit curve alone.
    !! Arm 2: per shell min(internal, cross) -- the conservative blend; the internal curve
    !! saturates high, so the blend defaults to the honest curve except where matching noise
    !! makes the cross curve spuriously HIGH, which the internal curve then caps. Unmatched
    !! components (spec par.4.4): arm 1 kills them (F=0 across the band -> kill-branch
    !! invtau2, the intended suppression of an unreproducible direction); arm 2 falls back to
    !! the internal curve alone, a finite weaker ridge that caps the damage of a matching
    !! failure -- which is why arm 2 is the recommended first A/B. Components beyond the
    !! record's per-fit rank (rank grew since t-1) get no ridge. Unmatched/unknown components
    !! are logged loudly. Fit X's ridge always uses fit X's OWN H (invariant, spec par.2.3):
    !! the record's per-fit H block for this fit, with the (shared) cross-fit F. The sign of
    !! the match is irrelevant here (FSC of sign-matched volumes is stored sign-resolved).
    subroutine xfsc_build_invtau2( ctx, params, rec, arm, l_fitb, it_eff, invtau2 )
        type(xfsc_ctx_t),      intent(in)  :: ctx
        class(parameters),     intent(in)  :: params
        type(crossfsc_record), intent(in)  :: rec
        integer,               intent(in)  :: arm, it_eff
        logical,               intent(in)  :: l_fitb
        real,                  intent(out) :: invtau2(:,:)   !< (ncomp, filtsz)
        real    :: fcur(ctx%filtsz), hcur(ctx%filtsz)
        integer :: qc, kc, kmatch_q, ncomp_fit, nun
        invtau2 = 0.0
        nun     = 0
        do qc = 1, size(invtau2,1)
            ncomp_fit = rec%ncomp_a
            if( l_fitb ) ncomp_fit = rec%ncomp_b
            if( qc > ncomp_fit )then
                ! component born after t-1: no curve, no H -> no ridge this iteration
                write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA XFSC REG: component ',qc, &
                    &' exceeds the previous record''s rank; no ridge applied to it'
                cycle
            endif
            if( l_fitb )then
                hcur = rec%h_b(:,qc)
            else
                hcur = rec%h_a(:,qc)
            endif
            ! locate this fit's component in the matched pairing
            kmatch_q = 0
            do kc = 1, rec%kmatch
                if( .not. l_fitb .and. rec%match_a(kc) == qc ) kmatch_q = kc
                if(       l_fitb .and. rec%match_b(kc) == qc ) kmatch_q = kc
            end do
            if( kmatch_q > 0 )then
                fcur = rec%fsc_cross(:,kmatch_q)
                if( arm == 2 )then
                    if( l_fitb )then
                        fcur = min(rec%fsc_int_b(:,qc), fcur)
                    else
                        fcur = min(rec%fsc_int_a(:,qc), fcur)
                    endif
                endif
            else
                nun = nun + 1
                if( arm == 2 )then
                    ! internal curve alone: finite, weaker ridge
                    if( l_fitb )then
                        fcur = rec%fsc_int_b(:,qc)
                    else
                        fcur = rec%fsc_int_a(:,qc)
                    endif
                    write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA XFSC UNMATCHED component ',qc, &
                        &': no partner cleared the match floor; arm-2 ridge from the internal curve alone'
                else
                    ! arm 1: F=0 across the band -> kill-branch invtau2 everywhere in band
                    fcur = -1.0
                    write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA XFSC UNMATCHED component ',qc, &
                        &': no partner cleared the match floor; arm-1 ridge SUPPRESSES it (kill branch)'
                endif
            endif
            call crossfsc_to_invtau2(fcur, hcur, params%tau, ctx%klo, invtau2(qc,:))
        end do
        if( nun > 0 )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA XFSC REG it=',it_eff,'  ',nun, &
                &' unmatched component(s) this iteration (see lines above)'
        endif
        call flush(logfhandle)
    end subroutine xfsc_build_invtau2



    !> crossfsc state teardown (the artifact on disk is the persistent series); the per-fit
    !! payloads are freed by kill_probe_fit / scope exit.
    module subroutine xfsc_teardown( ctx )
        type(xfsc_ctx_t), intent(inout) :: ctx
        if( .not. ctx%l_any ) return
        call crossfsc_kill(ctx%xf)
    end subroutine xfsc_teardown

    !> CROSSFSC artifact writer, PAIRED mode (spec par.2.2): one honest paired=1 record per
    !! iteration, appended after BOTH fits' master tails complete -- the impl-map step-3 hook in
    !! the paired loop. Per-fit blocks come from each fit's own stashes (internal e/o FSC curves
    !! + Gamma at the BAND/RANK site; H = the fit's own e+o harvest sum per spec par.2.3).
    !!
    !! MATCHED BLOCKS, v1-minimal (documented deviation from the full par.5 contract): the signed
    !! pairing is GREEDY |cos| matching on the two fits' realized in-memory bases
    !! (fit%prev_real -- the delivered orthonormal basis of this iteration, band-limited to the
    !! working band and soft-masked identically in both fits, so the raw real-space cosine IS the
    !! comparison-band masked cosine), rather than mask-mean-subtracted per-half varimax +
    !! exhaustive/Hungarian matching. fsc_cross is then the honest cross-fit FSC between fit A's
    !! component and sign * fit B's matched component (image%fsc, the em_iter Wiener-site call).
    !! The offline instrument (~/ribo_local/rank_criterion.py) remains the reference for rank
    !! decisions; upgrading the in-engine matching to matched varimax is step-3 follow-up work.
    module subroutine xfsc_paired_record( ctx, params, fits, it_eff )
        type(xfsc_ctx_t),  intent(inout) :: ctx
        class(parameters), intent(in)    :: params
        type(probe_fit_t), intent(inout) :: fits(2)
        integer,           intent(in)    :: it_eff
        type(crossfsc_record) :: rec
        type(image) :: imga, imgb
        real,     pointer     :: pra(:,:,:), prb(:,:,:)
        real(dp), allocatable :: cosm(:,:)
        real,     allocatable :: corrs(:)
        logical,  allocatable :: used_a(:), used_b(:)
        integer  :: nc_a, nc_b, km, q, p, k, best_q, best_p
        real(dp) :: nrm_a, nrm_b, c, best_c
        if( .not. (allocated(fits(1)%xf_fscq) .and. allocated(fits(2)%xf_fscq) .and. &
            &allocated(fits(1)%xf_h_e) .and. allocated(fits(2)%xf_h_e)) ) &
            &THROW_HARD('xfsc_paired_record: writer payloads missing (l_xf_harvest not set?)')
        if( .not. (allocated(fits(1)%prev_real) .and. allocated(fits(2)%prev_real)) ) &
            &THROW_HARD('xfsc_paired_record: realized bases (prev_real) missing')
        ! delivered ranks; the stashes were taken at entry rank >= delivered rank, truncated here
        nc_a = min(fits(1)%ncomp, size(fits(1)%xf_fscq,2))
        nc_b = min(fits(2)%ncomp, size(fits(2)%xf_fscq,2))
        km   = min(nc_a, nc_b)
        if( .not. ctx%l_loaded .and. ctx%xf%nrec == 0 )then
            ! file header, frozen at first record: khi_cmp = khi_full at the production lp
            ctx%xf%paired     = 1
            ctx%xf%pairing_id = ctx%pairing_id
            ctx%xf%box_crop   = params%box_crop
            ctx%xf%filtsz     = ctx%filtsz
            ctx%xf%smpd_crop  = params%smpd_crop
            ctx%xf%khi_full   = fits(1)%khi_full
            ctx%xf%khi_cmp    = fits(1)%khi_full
        endif
        rec%it_eff     = it_eff
        rec%ncomp_a    = nc_a
        rec%ncomp_b    = nc_b
        rec%kmatch     = km
        ! per-fit internal-FSC bands by the par.5.2 criterion (deepest crossing + margin, clamped
        ! to the full band), recorded for offline rank/band diagnostics
        rec%khi_a      = min(fits(1)%khi_full, max(1, crossfsc_khi_deepest(fits(1)%xf_fscq, nc_a, XFSC_CRIT) + 2))
        rec%khi_b      = min(fits(2)%khi_full, max(1, crossfsc_khi_deepest(fits(2)%xf_fscq, nc_b, XFSC_CRIT) + 2))
        rec%khi_shared = fits(1)%khi_full
        rec%reg_mode   = ctx%reg_active
        rec%march_on   = 0
        allocate(rec%match_a(km), rec%match_b(km), rec%match_sign(km))
        allocate(rec%match_cos(km), rec%fsc_cross(ctx%filtsz,km))
        allocate(rec%fsc_int_a(ctx%filtsz,nc_a), rec%fsc_int_b(ctx%filtsz,nc_b))
        allocate(rec%eigvals_a(nc_a), rec%eigvals_b(nc_b))
        allocate(rec%h_a(ctx%filtsz,nc_a), rec%h_b(ctx%filtsz,nc_b), rec%cnt(ctx%filtsz))
        rec%fsc_int_a = fits(1)%xf_fscq(:,1:nc_a)
        rec%fsc_int_b = fits(2)%xf_fscq(:,1:nc_b)
        rec%eigvals_a = fits(1)%xf_gam(1:nc_a)
        rec%eigvals_b = fits(2)%xf_gam(1:nc_b)
        ! a paired-engine fit stores its OWN e+o sampling sum (spec par.2.3): H is the per-shell
        ! mean of (rho_e + rho_o) diagonals = mean_e + mean_o (same voxel counts)
        rec%h_a = fits(1)%xf_h_e(:,1:nc_a) + fits(1)%xf_h_o(:,1:nc_a)
        rec%h_b = fits(2)%xf_h_e(:,1:nc_b) + fits(2)%xf_h_o(:,1:nc_b)
        rec%cnt = fits(1)%xf_cnt            ! shared lattice: identical for both fits
        ! ---- greedy signed |cos| matching on the realized bases ----
        allocate(cosm(nc_a,nc_b), used_a(nc_a), used_b(nc_b))
        used_a = .false.; used_b = .false.
        do p = 1, nc_b
            call fits(2)%prev_real(p)%get_rmat_ptr(prb)
            do q = 1, nc_a
                call fits(1)%prev_real(q)%get_rmat_ptr(pra)
                nrm_a = sqrt(sum(real(pra,dp)**2))
                nrm_b = sqrt(sum(real(prb,dp)**2))
                cosm(q,p) = sum(real(pra,dp)*real(prb,dp)) / max(nrm_a*nrm_b, DTINY)
            end do
        end do
        allocate(corrs(ctx%filtsz))
        do k = 1, km
            best_c = -1.d0; best_q = 0; best_p = 0
            do q = 1, nc_a
                if( used_a(q) ) cycle
                do p = 1, nc_b
                    if( used_b(p) ) cycle
                    c = abs(cosm(q,p))
                    if( c > best_c )then
                        best_c = c; best_q = q; best_p = p
                    endif
                end do
            end do
            used_a(best_q) = .true.
            used_b(best_p) = .true.
            rec%match_a(k)    = best_q
            rec%match_b(k)    = best_p
            rec%match_sign(k) = merge(1, -1, cosm(best_q,best_p) >= 0.d0)
            rec%match_cos(k)  = real(best_c)
            ! honest cross-fit FSC between A's component and sign * B's matched component
            call imga%copy(fits(1)%prev_real(best_q))
            call imgb%copy(fits(2)%prev_real(best_p))
            if( rec%match_sign(k) < 0 ) call imgb%mul(-1.0)
            call imga%fft
            call imgb%fft
            call imga%fsc(imgb, corrs)
            rec%fsc_cross(:,k) = corrs
            call imga%kill
            call imgb%kill
        end do
        deallocate(cosm, used_a, used_b, corrs)
        rec%s_stop = crossfsc_stop_stat(rec, ctx%xf%khi_cmp)
        write(logfhandle,'(A,I0,A,I0,A,F6.3,A,F6.3)') '>>> FLEX_PCA XFSC PAIRED record it=', &
            &it_eff,'  matched pairs=',km,'  |cos| lead=',rec%match_cos(1),'  S(t)=',real(rec%s_stop)
        call crossfsc_append(ctx%xf, rec)
        call crossfsc_write(ctx%xf)
        call crossfsc_kill_record(rec)
    end subroutine xfsc_paired_record

end submodule simple_flex_pca_em_crossfsc
