!@descr: flex_pca probe fit E-step: stage setup, polar bank/former, posterior solve, insertion and reduction
submodule (simple_flex_probe_fit) simple_flex_probe_fit_estep
!$ use omp_lib, only: omp_get_thread_num, omp_get_wtime
use simple_core_module_api, only: cmplx_zero, dtiny, kbalpha, kbinterpol, kbwinsz, logfhandle, &
    &maximgbatchsz, oris, osmpl_pad_fac, tic, timer_int_kind, toc
use simple_math,                          only: ceil_div, floor_div
use simple_flex_pca_polar,                only: polar_grid_build, polar_project_recs, polar_relative_inplane, &
    &polar_assign_directions, polar_sample_particle_fused
use simple_flex_reconstructor_latent_ops, only: latent_projection_weights, &
    &weighted_expanded_cmat, planes_batch_load, LATENT_WDIM
use simple_flex_pca_posterior,            only: probe_solve_plain, probe_solve_mix, mcfa_init, mcfa_condition, mcfa_mstep
use simple_flex_pca_basis,                only: cov_image_mask_radius
use simple_matcher_3Drec,                 only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io,               only: prepimgbatch
implicit none
#include "simple_local_flags.inc"

contains

    !> Polar E-step bank for one fit: grid geometry + pose-fixed direction assignment (once per
    !! stage), then the per-iteration shared-direction bank + ring Gram tables at the fit's
    !! current rank, restricted to the directions the fit's current window touches.
    module subroutine fit_polar_bank_build( build, fit, mean_rec, fpl1, nthr )
        type(builder),        intent(inout) :: build
        type(flex_probe_fit), intent(inout) :: fit
        type(reconstructor),  intent(inout) :: mean_rec
        type(fplane_type),    intent(in)    :: fpl1
        integer,              intent(in)    :: nthr
        type(oris) :: dirs_es
        type(ori)  :: o_es
        real,    allocatable :: rmatp_es(:,:,:), nrmp_es(:,:)
        real     :: ca1, sa1
        real(dp) :: t_bank0, t_bank1
        integer  :: i, q, r, ir, id_es, ithr, jx_es, kx_es
        if( .not. fit%estep%l_pol_grid )then
            ! grid geometry from the first prepped plane: the identical derivation
            ! embed_accumulate_polar uses, so polar embed and polar E-step share
            ! band/quadrature conventions. No noise rings: sig2 arrives as sig2_eff.
            fit%estep%ph0_es  = lbound(fpl1%cmplx_plane,1)
            fit%estep%pk0_es  = lbound(fpl1%cmplx_plane,2)
            fit%estep%hlo_es  = ceil_div (lbound(fpl1%cmplx_plane,1), OSMPL_PAD_FAC)
            fit%estep%hhi_es  = floor_div(ubound(fpl1%cmplx_plane,1), OSMPL_PAD_FAC)
            fit%estep%klo_es  = ceil_div (lbound(fpl1%cmplx_plane,2), OSMPL_PAD_FAC)
            fit%estep%nyqr_es = mean_rec%get_lfny(1)
            fit%estep%nyqb_es = fit%estep%nyqr_es
            if( fpl1%nyq > 0 ) fit%estep%nyqb_es = min(fit%estep%nyqb_es, max(1, fpl1%nyq / OSMPL_PAD_FAC))
            ! hybrid split point: auto at 0.72*band -- the measured knee of the
            ! real-data ladder (10049 gate, band 11: rhyb 6 -> b err 8.9%, min z
            ! corr 0.949, lambda head 153/282/326; rhyb 8 -> b err 3.1%, min z
            ! corr 0.990, lambda 166/273/394 vs Cartesian 166/264/382).
            fit%estep%rhyb_es = nint(0.72*real(fit%estep%nyqb_es))
            fit%estep%rhyb_es   = max(0, min(fit%estep%rhyb_es, fit%estep%nyqb_es-1))
            fit%estep%l_pol_hyb = fit%estep%rhyb_es > 0
            if( fit%estep%l_pol_hyb )then
                call polar_grid_build(fit%estep%pg_es, fit%estep%rhyb_es+1, fit%estep%nyqb_es, &
                    &fit%estep%hlo_es, fit%estep%hhi_es, fit%estep%klo_es, fit%estep%ph0_es, fit%estep%pk0_es, &
                    &gate_lo=fit%estep%rhyb_es*(fit%estep%rhyb_es+1))
                ! exact-part lattice positions, cov_herm_inner's half-plane rule
                ! (k<=0; on the k=0 line only h<=0), shells 0..rhyb by the nint
                ! convention: h^2+k^2 <= rhyb*(rhyb+1). Raster order k-outer.
                fit%estep%npos_es = 0
                do kx_es = -fit%estep%rhyb_es, 0
                    do jx_es = -fit%estep%rhyb_es, merge(0, fit%estep%rhyb_es, kx_es == 0)
                        if( jx_es*jx_es + kx_es*kx_es > fit%estep%rhyb_es*(fit%estep%rhyb_es+1) ) cycle
                        fit%estep%npos_es = fit%estep%npos_es + 1
                    end do
                end do
                allocate(fit%estep%hex_es(fit%estep%npos_es), fit%estep%kex_es(fit%estep%npos_es))
                fit%estep%npos_es = 0
                do kx_es = -fit%estep%rhyb_es, 0
                    do jx_es = -fit%estep%rhyb_es, merge(0, fit%estep%rhyb_es, kx_es == 0)
                        if( jx_es*jx_es + kx_es*kx_es > fit%estep%rhyb_es*(fit%estep%rhyb_es+1) ) cycle
                        fit%estep%npos_es = fit%estep%npos_es + 1
                        fit%estep%hex_es(fit%estep%npos_es) = jx_es
                        fit%estep%kex_es(fit%estep%npos_es) = kx_es
                    end do
                end do
                write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA POLAR HYBRID: exact &
                    &Cartesian statistics for shells 0..',fit%estep%rhyb_es,' (',fit%estep%npos_es,' lattice &
                    &points/particle incl. DC), rings ',fit%estep%rhyb_es+1,'..',fit%estep%nyqb_es
                call flush(logfhandle)
            else
                call polar_grid_build(fit%estep%pg_es, 1, fit%estep%nyqb_es, &
                    &fit%estep%hlo_es, fit%estep%hhi_es, fit%estep%klo_es, &
                    &fit%estep%ph0_es, fit%estep%pk0_es)
            endif
            fit%estep%nsamp_es = fit%estep%pg_es%nsamp; fit%estep%nsamp2_es = 2*fit%estep%nsamp_es; fit%estep%nk_es = fit%estep%pg_es%nk
            ! shared direction table: same refspiral + count derivation as the polar embed
            fit%estep%ndir_es = cov_polar_ndir(fit%spec%npp)
            call dirs_es%new(fit%estep%ndir_es, is_ptcl=.false.)
            call build%pgrpsyms%build_refspiral(dirs_es)
            allocate(fit%estep%rmatb_es(3,3,fit%estep%ndir_es), fit%estep%nrmb_es(3,fit%estep%ndir_es))
            do id_es = 1, fit%estep%ndir_es
                fit%estep%rmatb_es(:,:,id_es) = dirs_es%get_mat(id_es)
                fit%estep%nrmb_es(:,id_es)    = fit%estep%rmatb_es(3,:,id_es)
            end do
            call dirs_es%kill
            ! pose-fixed per-particle (direction, in-plane) assignment, once per stage
            allocate(rmatp_es(3,3,fit%spec%npp), nrmp_es(3,fit%spec%npp), fit%estep%dir_es(fit%spec%npp), fit%estep%cae(fit%spec%npp), fit%estep%sae(fit%spec%npp))
            do i = 1, fit%spec%npp
                call build%spproj_field%get_ori(fit%spec%ppinds(i), o_es)
                rmatp_es(:,:,i) = o_es%get_mat()
                nrmp_es(:,i)    = rmatp_es(3,:,i)
            end do
            call o_es%kill
            call polar_assign_directions(nrmp_es, fit%spec%npp, fit%estep%nrmb_es, fit%estep%ndir_es, fit%estep%dir_es)
            do i = 1, fit%spec%npp
                call polar_relative_inplane(rmatp_es(:,:,i), fit%estep%rmatb_es(:,:,fit%estep%dir_es(i)), ca1, sa1)
                fit%estep%cae(i) = ca1; fit%estep%sae(i) = sa1
            end do
            deallocate(rmatp_es, nrmp_es)
            allocate(fit%estep%dused_es(fit%estep%ndir_es))
            fit%estep%l_pol_grid = .true.
            write(logfhandle,'(A,I0,A,I0,A,I0,A,F8.1,A,F8.1,A)') &
                &'>>> FLEX_PCA POLAR ESTEP BANK: ',fit%model%ncomp+1,' volumes x ',fit%estep%ndir_es, &
                &' directions x ',fit%estep%nsamp_es,' ring samples = ', &
                &4.d0*real(fit%estep%nsamp2_es,dp)*real(fit%model%ncomp+1,dp)*real(fit%estep%ndir_es,dp)/1.d6, &
                &' MB (+ ring tables ', &
                &8.d0*real(fit%model%ncomp*fit%model%ncomp+fit%model%ncomp+1,dp)*real(fit%estep%nk_es,dp)*real(fit%estep%ndir_es,dp)/1.d6,' MB)'
            call flush(logfhandle)
        endif
        ! (re)allocate at this iteration's rank (ncomp can change between iterations)
        if( allocated(fit%estep%UsallE) )then
            if( size(fit%estep%UsallE,2) /= fit%model%ncomp+1 ) deallocate(fit%estep%UsallE, fit%estep%CfE, fit%estep%Cm0E, fit%estep%c00E, &
                &fit%estep%UbankE, fit%estep%CspE, fit%estep%xws_es, fit%estep%wr_es, fit%estep%wrd_es, fit%estep%Reb_es)
        endif
        if( .not. allocated(fit%estep%UsallE) )then
            allocate(fit%estep%UsallE(fit%estep%nsamp2_es,0:fit%model%ncomp,fit%estep%ndir_es), fit%estep%CfE(fit%model%ncomp*fit%model%ncomp,fit%estep%nk_es,fit%estep%ndir_es), &
                &fit%estep%Cm0E(fit%model%ncomp,fit%estep%nk_es,fit%estep%ndir_es), fit%estep%c00E(fit%estep%nk_es,fit%estep%ndir_es))
            allocate(fit%estep%UbankE(fit%estep%nsamp_es,0:fit%model%ncomp,nthr), fit%estep%CspE(0:fit%model%ncomp,0:fit%model%ncomp,nthr))
            allocate(fit%estep%xws_es(fit%estep%nsamp2_es,nthr), fit%estep%wr_es(fit%estep%nk_es,nthr), fit%estep%wrd_es(fit%estep%nk_es,nthr), &
                &fit%estep%Reb_es(0:fit%model%ncomp,nthr))
        endif
        ! only the directions this iteration's window touches
        fit%estep%dused_es = .false.
        do i = 1, fit%spec%npp
            if( fit%estep%dir_es(i) > 0 ) fit%estep%dused_es(fit%estep%dir_es(i)) = .true.
        end do
        !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
        !$omp& private(id_es,ithr,q,r,ir,t_bank0,t_bank1)
        do id_es = 1, fit%estep%ndir_es
            if( .not. fit%estep%dused_es(id_es) ) cycle
            ithr = omp_get_thread_num() + 1
            t_bank0 = omp_get_wtime()
            call polar_project_recs(mean_rec, fit%model%basis_recs, fit%model%ncomp, fit%estep%rmatb_es(:,:,id_es), &
                &fit%estep%pg_es, fit%estep%UbankE(:,:,ithr))
            do q = 0, fit%model%ncomp
                do r = 1, fit%estep%nsamp_es
                    fit%estep%UsallE(2*r-1,q,id_es) = fit%estep%pg_es%sqwq(r)*real (fit%estep%UbankE(r,q,ithr))
                    fit%estep%UsallE(2*r,  q,id_es) = fit%estep%pg_es%sqwq(r)*aimag(fit%estep%UbankE(r,q,ithr))
                end do
            end do
            do ir = 1, fit%estep%nk_es
                call polar_ring_gram(fit%estep%UsallE(1,0,id_es), fit%estep%nsamp2_es, fit%model%ncomp, fit%estep%pg_es%rbeg(ir), &
                    &fit%estep%pg_es%rend(ir)-fit%estep%pg_es%rbeg(ir)+1, fit%estep%CspE(0,0,ithr), fit%estep%CfE(1,ir,id_es), &
                    &fit%estep%Cm0E(1,ir,id_es))
                fit%estep%c00E(ir,id_es) = polar_ring_selfpower(fit%estep%UsallE(1,0,id_es), fit%estep%nsamp2_es, &
                    &fit%estep%pg_es%rbeg(ir), fit%estep%pg_es%rend(ir)-fit%estep%pg_es%rbeg(ir)+1)
            end do
            t_bank1 = omp_get_wtime()
            fit%diag%sec_bank_thr(ithr) = fit%diag%sec_bank_thr(ithr) + (t_bank1 - t_bank0)
        end do
        !$omp end parallel do
    end subroutine fit_polar_bank_build

    !> Production polar shared-direction former for one particle: banded mean projection,
    !! fused polar sampling, bank GEMVs, hybrid low-k exact statistics, contrast fit.
    module subroutine fit_estep_former_polar( fit, mean_rec, o, fpl, row, ithr, a, e_mm, myv )
        type(flex_probe_fit), intent(inout) :: fit
        type(reconstructor),  intent(inout) :: mean_rec
        class(ori),           intent(inout) :: o
        type(fplane_type),    intent(inout) :: fpl
        integer,              intent(in)    :: row, ithr
        real(dp),             intent(out)   :: a, e_mm, myv
        integer  :: idp_es, q
        real     :: taz_es
        real(dp) :: twp0, twp1, twp2
        twp0 = omp_get_wtime()
        idp_es  = fit%estep%dir_es(row)
        ! mean only, in Cartesian: the M-step backprojects y - a*(T mu) and
        ! there is no polar->volume adjoint. Banded variant: identical
        ! interpolation, none of project_fplane's per-call full-plane
        ! zero-fill + ctfsq/transfer copies (measured as the bulk of the
        ! polar project bucket at the native padded plane).
        call project_fplane_mean_banded(mean_rec, o, fpl, &
            &fit%iter%mean_fpl(ithr))
        ! polar-sample the prepped data plane ONCE at (bank direction, relative
        ! in-plane angle): CTF amplitude, shift phase and per-shell whitening all
        ! ride in from the same cmplx/transfer planes the Cartesian former reads.
        ! Fused sampler: one KB geometry per ring sample shared by both plane
        ! gathers, packed output written in place, no per-call allocations --
        ! bit-identical statistics (see polar_sample_particle_fused).
        call polar_sample_particle_fused(fpl%cmplx_plane, fpl%transfer_plane, &
            &fit%estep%pg_es, fit%estep%cae(row), fit%estep%sae(row), fit%estep%xws_es(:,ithr), fit%estep%wr_es(:,ithr), taz_es)
        fit%estep%wrd_es(:,ithr) = real(fit%estep%wr_es(:,ithr), dp)
        ! b and the mean row, exact per-sample CTF: one GEMV against the bank
        call sgemv('T', fit%estep%nsamp2_es, fit%model%ncomp+1, 1.0, fit%estep%UsallE(1,0,idp_es), fit%estep%nsamp2_es, &
            &fit%estep%xws_es(1,ithr), 1, 0.0, fit%estep%Reb_es(0,ithr), 1)
        ! G, c, e_mm by the radial factorisation: ring Grams x per-ring mean |T|^2
        call dgemv('N', fit%model%ncomp*fit%model%ncomp, fit%estep%nk_es, 1.d0, fit%estep%CfE(1,1,idp_es), fit%model%ncomp*fit%model%ncomp, &
            &fit%estep%wrd_es(1,ithr), 1, 0.d0, fit%iter%Gth(1,1,ithr), 1)
        call dgemv('N', fit%model%ncomp, fit%estep%nk_es, 1.d0, fit%estep%Cm0E(1,1,idp_es), fit%model%ncomp, &
            &fit%estep%wrd_es(1,ithr), 1, 0.d0, fit%iter%cth(1,ithr), 1)
        e_mm = dot_product(fit%estep%c00E(:,idp_es), fit%estep%wrd_es(:,ithr))
        myv  = real(fit%estep%Reb_es(0,ithr), dp)
        do q = 1, fit%model%ncomp
            fit%iter%bth(q,ithr) = real(fit%estep%Reb_es(q,ithr), dp)
        end do
        twp1 = omp_get_wtime()
        fit%diag%sec_ring_thr(ithr) = fit%diag%sec_ring_thr(ithr) + (twp1 - twp0)
        ! hybrid: the low-k shells enter as exact Cartesian statistics
        if( fit%estep%l_pol_hyb )then
            call polar_hybrid_exact_accum(mean_rec, fit%model%basis_recs, &
                &fit%model%ncomp, o, fpl, fit%estep%hex_es, fit%estep%kex_es, fit%estep%npos_es, &
                &fit%iter%Gth(:,:,ithr), fit%iter%bth(:,ithr), fit%iter%cth(:,ithr), e_mm, myv)
        endif
        twp2 = omp_get_wtime()
        fit%diag%sec_exact_thr(ithr) = fit%diag%sec_exact_thr(ithr) + (twp2 - twp1)
        a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
    end subroutine fit_estep_former_polar

    !> Shared per-particle tail of every CPU former: posterior solve (plain or mixture), latent
    !! and density batch rows, Gamma/likelihood accumulation, and the in-place mean-subtracted
    !! residual the M-step backprojects.
    module subroutine fit_estep_solve_stats( fit, fpl, i, row, ithr, a )
        type(flex_probe_fit), intent(inout) :: fit
        type(fplane_type),    intent(inout) :: fpl
        integer,              intent(in)    :: i, row, ithr
        real(dp),             intent(in)    :: a
        integer  :: q, r
        logical  :: lok
        real(dp) :: ldA, qml, nll_mix_add, twp1, twp2
        twp1 = omp_get_wtime()
        ! Posterior precision A = (a^2/sig2) G + Gamma^-1. The whole normal system is already
        ! scaled by 1/sig2, so Cov[z|y] = A^-1 exactly -- no further sig2 factor.
        if( fit%history%l_mix_used )then
            ! ---- MCFA E-step via the ONE shared solver (probe_solve_mix) ----
            call probe_solve_mix(fit, ithr, i, a, ldA, lok, nll_mix_add)
            if( lok ) fit%iter%nll_thr(ithr) = fit%iter%nll_thr(ithr) + nll_mix_add
        else
            call probe_solve_plain(fit, ithr, a, ldA, lok, qml)
            if( lok ) fit%iter%nll_thr(ithr) = fit%iter%nll_thr(ithr) + ldA - qml
        endif
        fit%model%z(row,:)          = fit%iter%zth(:,ithr)
        fit%iter%zbatch(:,i)       = fit%iter%zth(:,ithr)
        ! EM sufficient statistic E[z z'|y] = z z' + Cov[z|y]. BOTH the coupled M-step normal
        ! matrix and the Gamma update below need it. Dropping Cov underestimates Gamma, which
        ! tightens the prior, which shrinks z further: the bias compounds across iterations.
        ! Under the mixture E[zz'|y] = A^-1 + sum_k r_k m_k m_k', NOT A^-1 + E[z]E[z]'
        ! -- the between-component spread is real posterior variance; probe_solve_mix
        ! has already written dens(:,:,i) in that case.
        if( .not. fit%history%l_mix_used )then
            do r = 1, fit%model%ncomp
                do q = 1, fit%model%ncomp
                    fit%iter%dens(q,r,i) = fit%iter%zth(q,ithr)*fit%iter%zth(r,ithr) + fit%iter%Ainvth(q,r,ithr)
                end do
            end do
        endif
        do q = 1, fit%model%ncomp
            fit%iter%gam_thr(q,ithr) = fit%iter%gam_thr(q,ithr) + fit%iter%dens(q,q,i)
        end do
        fit%iter%nval_thr(ithr)    = fit%iter%nval_thr(ithr) + 1
        fit%iter%valid(i)          = .true.
        fit%diag%gam_dbg(1,ithr) = fit%diag%gam_dbg(1,ithr) + sum([(fit%iter%Gth(q,q,ithr), q=1,fit%model%ncomp)])
        fit%diag%gam_dbg(2,ithr) = fit%diag%gam_dbg(2,ithr) + dot_product(fit%iter%bth(:,ithr), fit%iter%bth(:,ithr))
        fit%diag%gam_dbg(3,ithr) = fit%diag%gam_dbg(3,ithr) + dot_product(fit%iter%cth(:,ithr), fit%iter%cth(:,ithr))
        fit%diag%gam_dbg(4,ithr) = fit%diag%gam_dbg(4,ithr) + a
        ! residual observation r_i = y - a*(T mu) in place (transfer/ctfsq intact for backprojection)
        call subtract_mean_banded(fpl, fit%iter%mean_fpl(ithr), real(a), fit%estep%nyqr_es)
        twp2 = omp_get_wtime()
        ! This bucket is the complete per-particle tail: posterior solve,
        ! posterior moments, accumulation bookkeeping and residual subtraction.
        fit%diag%sec_solve_thr(ithr) = fit%diag%sec_solve_thr(ithr) + (twp2 - twp1)
    end subroutine fit_estep_solve_stats

    !> Per-iteration bank of a polar fit, built at the first prepped batch of the iteration
    !! (basis_recs change every M-step, so the bank cannot be cached across iterations; the grid
    !! and the pose-fixed direction assignment are built once per stage). A no-op for the
    !! Cartesian formulation and once the iteration's bank exists.
    module subroutine fit_estep_bank_prepare( fit, params, build, mean_rec, fpl1, nthr, it_eff, tag )
        class(flex_probe_fit), intent(inout) :: fit
        class(parameters),     intent(inout) :: params
        type(builder),         intent(inout) :: build
        type(reconstructor),   intent(inout) :: mean_rec
        type(fplane_type),     intent(in)    :: fpl1
        integer,               intent(in)    :: nthr, it_eff
        character(len=*),      intent(in)    :: tag   !< ' ' single fit; '  fit=A' / '  fit=B' paired
        integer(timer_int_kind) :: t_bank
        if( fit%estep%l_pol_bank_it ) return
        t_bank = tic()
        call fit_polar_bank_build(build, fit, mean_rec, fpl1, nthr)
        fit%estep%sec_bank      = fit%estep%sec_bank + real(toc(t_bank))
        fit%estep%l_pol_bank_it = .true.
        write(logfhandle,'(A,I0,A,A,I0,A,I0,A,F7.1)') '>>> FLEX_PCA POLAR ESTEP BANK it=', it_eff, tag, &
            &' directions built=', count(fit%estep%dused_es), ' of ', fit%estep%ndir_es, &
            &'  build seconds=', fit%estep%sec_bank
        call flush(logfhandle)
    end subroutine fit_estep_bank_prepare

    !> Reset the batch rows before a batch is formed.
    module subroutine fit_estep_batch_begin( fit, batchsz )
        class(flex_probe_fit), intent(inout) :: fit
        integer,               intent(in)    :: batchsz
        fit%iter%valid(:batchsz)    = .false.
        fit%iter%zbatch(:,:batchsz) = 0.d0
        fit%iter%dens(:,:,:batchsz) = 0.d0
    end subroutine fit_estep_batch_begin

    !> One particle: the formulation's former (G, b, c, e_mm, myv and the contrast) then the shared
    !! posterior tail (solve, batch rows, Gamma/likelihood accounting, the in-place residual).
    !! Thread-safe over distinct (i, ithr); called from the batch loops.
    module subroutine fit_estep_particle( fit, mean_rec, o, fpl, row, i, ithr )
        class(flex_probe_fit), intent(inout) :: fit
        type(reconstructor),   intent(inout) :: mean_rec
        class(ori),            intent(inout) :: o
        type(fplane_type),     intent(inout) :: fpl
        integer,               intent(in)    :: row, i, ithr
        real(dp) :: a, e_mm, myv
        call fit_estep_former_polar(fit, mean_rec, o, fpl, row, ithr, a, e_mm, myv)
        call fit_estep_solve_stats(fit, fpl, i, row, ithr, a)
    end subroutine fit_estep_particle

    !> Halfset masks + the CPU coupled M-step insertion for one batch of one fit.
    module subroutine fit_batch_insert( build, fit, orientations, fpls, eo, batchsz )
        type(builder),         intent(inout) :: build
        class(flex_probe_fit), intent(inout) :: fit
        type(ori),             intent(inout) :: orientations(:)
        type(fplane_type),     intent(inout) :: fpls(:)
        integer,               intent(in)    :: eo(:), batchsz
        integer :: i
        do i = 1, batchsz
            fit%iter%valid_e(i) = fit%iter%valid(i) .and. eo(i)==0
            fit%iter%valid_o(i) = fit%iter%valid(i) .and. eo(i)==1
        end do
        call fit%mstep%accumulate_batch(build, orientations, fpls, fit%iter%zbatch, fit%iter%dens, &
            &fit%iter%valid_e, fit%iter%valid_o, batchsz)
    end subroutine fit_batch_insert

    !> In-process end-of-iteration reductions (thread sums) and the MCFA mixture M-step or its
    !! initialisation. Runs on the shared-memory master and inside workers; the distributed master
    !! takes its sums from reduce_probe_parts.
    module subroutine fit_iter_reduce( fit, it_eff, nthr , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(flex_probe_fit),  intent(inout) :: fit
        integer,                intent(in)    :: it_eff, nthr
        integer :: q, kk2
        ! reduce the EM Gamma accumulator before ncomp is replaced below; Gamma travels between
        ! parts as a sum and is divided by the reduced nval
        fit%iter%nval = sum(fit%iter%nval_thr)
        ! E-step statistics summary: the polar and Cartesian formers must agree on G, b, c
        if( fit%iter%nval > 0 )then
            write(logfhandle,'(A,I0,A,ES12.4,A,ES12.4,A,ES12.4,A,F7.4)') &
                &'>>> FLEX_PCA PROBE ESTAT it=',it_eff,'  <trG>=',sum(fit%diag%gam_dbg(1,:))/real(fit%iter%nval,dp), &
                &'  <b.b>=',sum(fit%diag%gam_dbg(2,:))/real(fit%iter%nval,dp), &
                &'  <c.c>=',sum(fit%diag%gam_dbg(3,:))/real(fit%iter%nval,dp), &
                &'  <a>=',real(sum(fit%diag%gam_dbg(4,:))/real(fit%iter%nval,dp))
            call flush(logfhandle)
        endif
        ! in-process reduction of the likelihood accumulator (distributed runs sum it in
        ! reduce_probe_parts)
        fit%iter%nll_tot = sum(fit%iter%nll_thr)
        do q = 1, fit%model%ncomp
            fit%iter%gam_sum(q) = sum(fit%iter%gam_thr(q,:))
        end do
        ! MCFA: mixture M-step, or its (re)initialisation from the current latents; runs on the
        ! master with this iteration's accumulators complete
        if( fit%spec%l_mix_req .and. .not. rounds%is_worker() )then
            if( fit%history%l_mix_active )then
                block
                    real(dp), allocatable :: rr_sr(:), rr_sm(:,:), rr_smm(:,:,:), rr_sai(:,:)
                    integer  :: tt, kk3
                    allocate(rr_sr(fit%spec%kmix), rr_sm(fit%model%ncomp,fit%spec%kmix), &
                        &rr_smm(fit%model%ncomp,fit%model%ncomp,fit%spec%kmix), rr_sai(fit%model%ncomp,fit%model%ncomp))
                    rr_sr = 0.d0; rr_sm = 0.d0; rr_smm = 0.d0; rr_sai = 0.d0
                    do tt = 1, nthr
                        rr_sr  = rr_sr  + fit%history%mxa_sr(:,tt)
                        rr_sm  = rr_sm  + fit%history%mxa_sm(:,:,tt)
                        rr_sai = rr_sai + fit%history%mxa_sainv(:,:,tt)
                        do kk3 = 1, fit%spec%kmix
                            rr_smm(:,:,kk3) = rr_smm(:,:,kk3) + fit%history%mxa_smm(:,:,kk3,tt)
                        end do
                    end do
                    if( allocated(fit%history%dm_sr) .and. rounds%nparts() > 1 )then
                        ! distributed: the reduce already summed every part's thread-sums
                        call mcfa_mstep(fit, fit%history%dm_sr, fit%history%dm_sm, fit%history%dm_smm, fit%history%dm_sai)
                    else
                        call mcfa_mstep(fit, rr_sr, rr_sm, rr_smm, rr_sai)
                    endif
                    deallocate(rr_sr, rr_sm, rr_smm, rr_sai)
                end block
            else if( it_eff >= fit%spec%n_mix_warm )then
                allocate(fit%history%mix_xi(fit%model%ncomp,fit%spec%kmix), fit%history%mix_Om(fit%model%ncomp,fit%model%ncomp), fit%history%mix_Ominv(fit%model%ncomp,fit%model%ncomp), &
                    &fit%history%mix_pi(fit%spec%kmix), fit%history%mix_Omxi(fit%model%ncomp,fit%spec%kmix), fit%history%mix_xiOx(fit%spec%kmix), fit%history%mix_lpi(fit%spec%kmix))
                if( rounds%nparts() > 1 .and. fit%history%dm_nz > 0 )then
                    ! distributed: seed from the pooled per-part latent subsample
                    call mcfa_init(fit%history%dm_z(:fit%history%dm_nz,:), fit%history%dm_nz, fit%model%ncomp, fit%spec%kmix, fit%iter%gam_sum, fit%iter%nval, &
                        &fit%history%mix_xi, fit%history%mix_pi, fit%history%mix_Om)
                else
                    call mcfa_init(fit%model%z, size(fit%model%z,1), fit%model%ncomp, fit%spec%kmix, fit%iter%gam_sum, fit%iter%nval, fit%history%mix_xi, fit%history%mix_pi, fit%history%mix_Om)
                endif
                call mcfa_condition(fit%model%ncomp, fit%spec%kmix == 1, fit%history%mix_Om, fit%history%mix_Ominv, fit%history%ldOm_mix)
                fit%history%l_mix_active = .true.
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA MIX initialised: K=',fit%spec%kmix, &
                    &'  at iteration ',it_eff
            endif
            if( fit%history%l_mix_active )then
                do kk2 = 1, fit%spec%kmix
                    fit%history%mix_Omxi(:,kk2) = matmul(fit%history%mix_Ominv, fit%history%mix_xi(:,kk2))
                    fit%history%mix_xiOx(kk2)   = dot_product(fit%history%mix_xi(:,kk2), fit%history%mix_Omxi(:,kk2))
                    fit%history%mix_lpi(kk2)    = log(max(fit%history%mix_pi(kk2), 1.d-12))
                end do
                write(logfhandle,'(A,I0,A,ES11.3,A,F7.4,A,F8.3)') '>>> FLEX_PCA MIX it=',it_eff, &
                    &'  logdetOm=',fit%history%ldOm_mix,'  pi_max=',maxval(fit%history%mix_pi), &
                    &'  |xi|_max=',sqrt(maxval(sum(fit%history%mix_xi**2, dim=1)))
                call flush(logfhandle)
            endif
        endif
        call fit%mstep%reduce_local
    end subroutine fit_iter_reduce

    !> The E-step accumulate pass of one iteration over the fits' merged read list (sorted by
    !! project row; each fit's ppinds is ascending -- hazard 1): the shared-memory master runs it
    !! over the full selection, each distributed worker over its fromp/top shard. Every particle
    !! is formed against exactly the basis of the fit that owns it (per-particle cost unchanged;
    !! what the paired engine doubles is basis memory, the per-fit accumulators, the master's
    !! per-voxel solves and the per-iteration bank build). Owns its read/prep buffers; the caller
    !! owns the fits' accumulators (iter_begin / iter_reduce).
    module subroutine fit_estep_pass( params, build, plane_store, fits, means, nfits, it_eff, nthr )
        class(parameters),    intent(inout) :: params
        type(builder),        intent(inout) :: build
        class(flex_plane_store), intent(inout) :: plane_store
        integer,              intent(in)    :: nfits
        type(flex_probe_fit), intent(inout) :: fits(nfits)
        type(flex_mean_ref),  intent(in)    :: means(nfits)
        integer,              intent(in)    :: it_eff, nthr
        type(fplane_type), allocatable :: fpls(:)
        type(ori),         allocatable :: orientations(:)
        integer,           allocatable :: eo(:)
        integer,           allocatable :: mrg_pinds(:), mrg_owner(:), mrg_row(:)
        character(len=:),  allocatable :: pfx, btag, ftag
        integer  :: f, i, ithr, nmrg, ibatch, batchlims(2), batchsz, row, head(nfits), fmin
        real(timer_int_kind) :: sec_read, sec_prep, sec_estep, sec_ins
        logical  :: l_pcache
        integer(timer_int_kind) :: t_sec
        pfx = 'PROBE'
        if( nfits == 2 ) pfx = 'PAIRED'
        sec_read = 0.; sec_prep = 0.; sec_estep = 0.; sec_ins = 0.
        allocate(orientations(MAXIMGBATCHSZ), eo(MAXIMGBATCHSZ))
        ! ---- merged read list: linear k-way merge of the fits' windows, sorted by project row
        ! (ties go to the lower fit) ----
        nmrg = sum(fits(1:nfits)%spec%npp)
        allocate(mrg_pinds(nmrg), mrg_owner(nmrg), mrg_row(nmrg))
        head = 1
        do i = 1, nmrg
            fmin = 0
            do f = 1, nfits
                if( head(f) > fits(f)%spec%npp ) cycle
                if( fmin == 0 )then
                    fmin = f
                elseif( fits(f)%spec%ppinds(head(f)) < fits(fmin)%spec%ppinds(head(fmin)) )then
                    fmin = f
                endif
            end do
            if( fmin == 0 ) THROW_HARD('merged read list lost rows')
            mrg_pinds(i) = fits(fmin)%spec%ppinds(head(fmin)); mrg_owner(i) = fmin; mrg_row(i) = head(fmin)
            head(fmin) = head(fmin) + 1
        end do
        ! downscaled-particle cache (cache=yes): a cache-served batch is read at box_crop into
        ! cropped planes (init_rec cropped=, prepimgbatch at box_crop) and prepped with
        ! cached=.true. (already noise-normalised when the entry was written)
        l_pcache = plane_store%cache_in_use()
        call init_rec(params, build, MAXIMGBATCHSZ, fpls, cropped=l_pcache)
        if( l_pcache )then
            call prepimgbatch(params, build, MAXIMGBATCHSZ, box=params%box_crop, smpd=params%smpd_crop)
        else
            call prepimgbatch(params, build, MAXIMGBATCHSZ)
        endif
        do ibatch = 1, nmrg, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nmrg, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call planes_batch_load(plane_store, params, build, nmrg, mrg_pinds, batchlims, fpls, &
                &cov_image_mask_radius(params), l_pcache, sec_read, sec_prep)
            do i = 1, batchsz
                call build%spproj_field%get_ori(mrg_pinds(batchlims(1)+i-1), orientations(i))
                eo(i) = build%spproj_field%get_eo(mrg_pinds(batchlims(1)+i-1))
            end do
            ! per-fit polar bank, built once per iteration at the first prepped batch, then the
            ! batch rows reset
            do f = 1, nfits
                btag = ' '
                if( nfits == 2 ) btag = merge('  fit=A','  fit=B', f==1)
                call fits(f)%estep_bank_prepare(params, build, means(f)%p, fpls(1), nthr, it_eff, btag)
                call fits(f)%estep_batch_begin(batchsz)
            end do
            t_sec = tic()
            ! interleaved owners: each particle is formed against the basis of the fit that owns it
            !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) private(i,f,row,ithr)
            do i = 1, batchsz
                if( orientations(i)%isstatezero() ) cycle
                ithr = omp_get_thread_num() + 1
                f    = mrg_owner(batchlims(1)+i-1)
                row  = mrg_row(batchlims(1)+i-1)
                call fits(f)%estep_particle(means(f)%p, orientations(i), fpls(i), row, i, ithr)
            end do
            !$omp end parallel do
            sec_estep = sec_estep + toc(t_sec)
            t_sec = tic()
            ! M-step by halfset: Y_q += sum_i z_iq * backproject(r_i), and the coupled normal matrix
            ! rho(q,r) += sum_i |CTF|^2 E[z_iq z_ir]   (batched KB)
            do f = 1, nfits
                call fits(f)%mstep_insert_batch(build, orientations, fpls, eo, batchsz)
            end do
            sec_ins = sec_ins + toc(t_sec)
            if( mod(batchlims(2), 5*MAXIMGBATCHSZ) == 0 .or. batchlims(2) == nmrg )then
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA '//pfx//' PASS PARTICLES: ', &
                    &batchlims(2), ' / ', nmrg
                call flush(logfhandle)
            endif
        end do
        write(logfhandle,'(A,F7.1,A,F7.1,A,F7.1,A,F7.1)') '>>> FLEX_PCA '//pfx//' E-STEP WALL &
            &(wall-seconds): read=', sec_read, '  prep=', sec_prep, '  project+solve=', sec_estep, &
            &'  insert=', sec_ins
        do f = 1, nfits
            ftag = ''
            if( nfits == 2 ) ftag = merge('  fit=A','  fit=B', f==1)
            write(logfhandle,'(A,A,A,F7.1,A,F7.1,A,F7.1,A,F7.1)') &
                &'>>> FLEX_PCA '//pfx//' E-STEP BUCKETS', ftag, &
                &' (thread-seconds): bank=', sum(fits(f)%diag%sec_bank_thr), &
                &'  ring=', sum(fits(f)%diag%sec_ring_thr), &
                &'  exact_lowk=', sum(fits(f)%diag%sec_exact_thr), &
                &'  solve+moments+residual=', sum(fits(f)%diag%sec_solve_thr)
            write(logfhandle,'(A,I0,A,A,F7.1,A,F7.1)') '>>> FLEX_PCA POLAR ESTEP it=', it_eff, ftag, &
                &'  bank build seconds=', fits(f)%estep%sec_bank, '  estep seconds=', sec_estep
        end do
        call flush(logfhandle)
        do i = 1, size(orientations)
            call orientations(i)%kill
        end do
        call cleanup_rec_buffers(build, fpls)
        deallocate(orientations, eo, mrg_pinds, mrg_owner, mrg_row)
    end subroutine fit_estep_pass

    !> Is the polar (shared-direction) E-step former requested? Always on.
    logical module function cov_polar_enabled()
        cov_polar_enabled = .true.
    end function cov_polar_enabled

    !> Bank directions: ~40 particles per direction, clamped to [1000,4000] and even (build_refspiral).
    !! The bank holds all ndir directions, so memory scales with ndir.
    integer module function cov_polar_ndir( nptcls )
        integer, intent(in) :: nptcls
        integer :: v
        v = min(4000, max(1000, nptcls/40))
        v = 2*((v+1)/2)                                  ! build_refspiral needs an even count
        cov_polar_ndir = v
    end function cov_polar_ndir

    !> One ring's contribution to the Gram of the basis (columns 1..ncomp of Us) and to the mean
    !! cross term (column 0). `Cout` is the full ncomp x ncomp block flattened column-major, `Mout`
    !! the ncomp-vector <U_q, T mu>.
    module subroutine polar_ring_gram( Us, ldu, ncomp, row0, nrow, Csp, Cout, Mout )
        integer,  intent(in)    :: ldu, ncomp, row0, nrow
        real,     intent(in)    :: Us(ldu,0:ncomp)
        real,     intent(inout) :: Csp(0:ncomp,0:ncomp)      !< caller-owned scratch
        real(dp), intent(out)   :: Cout(ncomp*ncomp), Mout(ncomp)
        integer :: q, r, i0, n2
        i0 = 2*row0 - 1
        n2 = 2*nrow
        if( n2 <= 0 )then
            Cout = 0.d0; Mout = 0.d0
            return
        endif
        call ssyrk('U','T', ncomp+1, n2, 1.0, Us(i0,0), ldu, 0.0, Csp, ncomp+1)
        do r = 1, ncomp
            do q = 1, r
                Cout((r-1)*ncomp+q) = real(Csp(q,r), dp)
                Cout((q-1)*ncomp+r) = real(Csp(q,r), dp)
            end do
            Mout(r) = real(Csp(0,r), dp)
        end do
    end subroutine polar_ring_gram

    real(dp) module function polar_ring_selfpower( Us, ldu, row0, nrow )
        integer, intent(in) :: ldu, row0, nrow
        real,    intent(in) :: Us(ldu,0:*)
        integer :: j
        polar_ring_selfpower = 0.d0
        do j = 2*row0-1, 2*(row0+nrow-1)
            polar_ring_selfpower = polar_ring_selfpower + real(Us(j,0),dp)*real(Us(j,0),dp)
        end do
    end function polar_ring_selfpower

    !> project_fplane(apply_ctf_amp=.true.) for the mean, cmplx_plane only, same interpolation; the plane
    !! is zeroed only at (re)allocation and no ctfsq/transfer copies are made, so out-of-disc samples stay
    !! zero only while the disc (frlims/nyq) is fixed. subtract_mean_banded sweeps the same disc.
    module subroutine project_fplane_mean_banded( rec, o, fpl_ref, fpl_out )
        type(reconstructor), intent(in)    :: rec
        class(ori),          intent(inout) :: o
        type(fplane_type),   intent(in)    :: fpl_ref
        type(fplane_type),   intent(inout) :: fpl_out
        type(kbinterpol) :: kbwin
        real    :: rotmat(3,3), loc(3), loc_friedel(3), hrow(3)
        real    :: w3(LATENT_WDIM,LATENT_WDIM,LATENT_WDIM)
        integer :: fpllims_pd(3,2), fpllims(3,2), h, k, hp, kp, pf, iwinsz, win(2,3)
        integer :: h_sq, k_max_h, k_lo, k_hi, nyq_disk, nyq_eff
        logical :: l_conjg, l_realloc
        complex :: comp
        kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        iwinsz = ceiling(KBWINSZ - 0.5)
        fpl_out%frlims  = fpl_ref%frlims
        fpl_out%shconst = fpl_ref%shconst
        fpl_out%nyq     = fpl_ref%nyq
        l_realloc = .not. allocated(fpl_out%cmplx_plane)
        if( .not. l_realloc )then
            l_realloc = any(lbound(fpl_out%cmplx_plane) /= lbound(fpl_ref%cmplx_plane)) .or. &
                &any(ubound(fpl_out%cmplx_plane) /= ubound(fpl_ref%cmplx_plane))
        endif
        if( l_realloc )then
            if( allocated(fpl_out%cmplx_plane) ) deallocate(fpl_out%cmplx_plane)
            allocate(fpl_out%cmplx_plane(lbound(fpl_ref%cmplx_plane,1):ubound(fpl_ref%cmplx_plane,1), &
                &lbound(fpl_ref%cmplx_plane,2):ubound(fpl_ref%cmplx_plane,2)))
            fpl_out%cmplx_plane = CMPLX_ZERO
        endif
        rotmat      = o%get_mat()
        pf          = OSMPL_PAD_FAC
        fpllims_pd  = fpl_ref%frlims
        fpllims     = fpllims_pd
        fpllims(1,1)= ceil_div (fpllims_pd(1,1), pf)
        fpllims(1,2)= floor_div(fpllims_pd(1,2), pf)
        fpllims(2,1)= ceil_div (fpllims_pd(2,1), pf)
        fpllims(2,2)= floor_div(fpllims_pd(2,2), pf)
        nyq_eff = rec%get_lfny(1)
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
                ! interp_cmat_exp, verbatim (it is private to simple_reconstructor)
                l_conjg     = loc(1) < 0.
                loc_friedel = loc
                if( l_conjg ) loc_friedel = -loc_friedel
                win(1,:) = nint(loc_friedel)
                win(2,:) = win(1,:) + iwinsz
                win(1,:) = win(1,:) - iwinsz
                call kbwin%apod_mat_3d(loc_friedel, iwinsz, LATENT_WDIM, w3)
                comp = sum(w3 * rec%cmat_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)))
                if( l_conjg ) comp = conjg(comp)
                ! apply_ctf_amp=.true. semantics of project_fplane
                if( allocated(fpl_ref%transfer_plane) )then
                    fpl_out%cmplx_plane(hp,kp) = fpl_ref%transfer_plane(hp,kp) * comp
                else
                    fpl_out%cmplx_plane(hp,kp) = sqrt(max(0., fpl_ref%ctfsq_plane(hp,kp))) * comp
                endif
            end do
        end do
    end subroutine project_fplane_mean_banded

    !> Exact Cartesian G/b/c/e_mm/myv increments over the hybrid E-step's low-k positions (hex,kex), equal
    !! to what project_fplanes_mean_basis + cov_herm_inner add there, DC included; one KB window per position
    !! serves all ncomp+1 volumes. Rings sample these few, steep shells worst and would bias the latent scale.
    module subroutine polar_hybrid_exact_accum( rec0, recs, ncomp, o, fpl, hex, kex, npos, &
            &Gd, bd, cd, e_mm, myv )
        type(reconstructor), intent(in)    :: rec0
        type(reconstructor), intent(in)    :: recs(ncomp)
        integer,             intent(in)    :: ncomp, npos
        class(ori),          intent(inout) :: o
        type(fplane_type),   intent(in)    :: fpl
        integer,             intent(in)    :: hex(npos), kex(npos)
        real(dp),            intent(inout) :: Gd(ncomp,ncomp), bd(ncomp), cd(ncomp)
        real(dp),            intent(inout) :: e_mm, myv
        type(kbinterpol) :: kbwin
        real        :: rotmat(3,3), loc(3), wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM)
        integer     :: j, q, r, win(2,3), hp, kp, exp_lb(3), exp_ub(3), pf
        logical     :: l_conjg, l_tf
        complex     :: tf, yv, u0, val
        complex     :: uq(ncomp)
        complex(dp) :: u0d, yd
        kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        rotmat = o%get_mat()
        pf     = OSMPL_PAD_FAC
        exp_lb = lbound(rec0%cmat_exp)
        exp_ub = ubound(rec0%cmat_exp)
        l_tf   = allocated(fpl%transfer_plane)
        do j = 1, npos
            loc(1) = real(hex(j))*rotmat(1,1) + real(kex(j))*rotmat(2,1)
            loc(2) = real(hex(j))*rotmat(1,2) + real(kex(j))*rotmat(2,2)
            loc(3) = real(hex(j))*rotmat(1,3) + real(kex(j))*rotmat(2,3)
            l_conjg = loc(1) < 0.
            if( l_conjg ) loc = -loc
            call latent_projection_weights(kbwin, loc, win, wx, wy, wz)
            if( any(win(1,:) < exp_lb) .or. any(win(2,:) > exp_ub) ) cycle
            hp = pf*hex(j)
            kp = pf*kex(j)
            if( l_tf )then
                tf = fpl%transfer_plane(hp,kp)
            else
                tf = cmplx(sqrt(max(0., fpl%ctfsq_plane(hp,kp))), 0.)
            endif
            yv  = fpl%cmplx_plane(hp,kp)
            val = weighted_expanded_cmat(rec0, win, wx, wy, wz)
            if( l_conjg ) val = conjg(val)
            u0 = tf * val
            do q = 1, ncomp
                val = weighted_expanded_cmat(recs(q), win, wx, wy, wz)
                if( l_conjg ) val = conjg(val)
                uq(q) = tf * val
            end do
            u0d  = cmplx(u0, kind=dp)
            yd   = cmplx(yv, kind=dp)
            e_mm = e_mm + real(conjg(u0d)*u0d, dp)
            myv  = myv  + real(conjg(u0d)*yd,  dp)
            do q = 1, ncomp
                bd(q) = bd(q) + real(conjg(cmplx(uq(q),kind=dp))*yd,  dp)
                cd(q) = cd(q) + real(conjg(cmplx(uq(q),kind=dp))*u0d, dp)
                do r = q, ncomp
                    Gd(q,r) = Gd(q,r) + real(conjg(cmplx(uq(q),kind=dp))*cmplx(uq(r),kind=dp), dp)
                end do
            end do
        end do
        ! mirror the accumulated upper triangle (the ring dgemv filled both triangles already;
        ! the exact increments above touched q<=r only)
        do r = 1, ncomp
            do q = r+1, ncomp
                Gd(q,r) = Gd(r,q)
            end do
        end do
    end subroutine polar_hybrid_exact_accum

    !> Banded residual subtraction, fpl = fpl - a*mean over EXACTLY the disc the banded (or any
    !! full-plane) mean projection wrote. Everywhere outside that disc the mean plane is
    !! identically zero, so the full-array statement this replaces only rewrote unchanged values
    !! there -- another few MB of per-particle traffic at the native padded lattice for no effect.
    !! The loop bounds are the same expressions as project_fplane_mean_banded's, so written and
    !! subtracted sample sets coincide by construction.
    module subroutine subtract_mean_banded( fpl, mean_fpl, a, rec_nyq )
        type(fplane_type), intent(inout) :: fpl
        type(fplane_type), intent(in)    :: mean_fpl
        real,              intent(in)    :: a
        integer,           intent(in)    :: rec_nyq
        integer :: fpllims_pd(3,2), fpllims(3,2), h, k, hp, kp, pf
        integer :: h_sq, k_max_h, k_lo, k_hi, nyq_disk, nyq_eff
        pf          = OSMPL_PAD_FAC
        fpllims_pd  = fpl%frlims
        fpllims     = fpllims_pd
        fpllims(1,1)= ceil_div (fpllims_pd(1,1), pf)
        fpllims(1,2)= floor_div(fpllims_pd(1,2), pf)
        fpllims(2,1)= ceil_div (fpllims_pd(2,1), pf)
        fpllims(2,2)= floor_div(fpllims_pd(2,2), pf)
        nyq_eff = rec_nyq
        if( fpl%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpl%nyq / pf))
        nyq_disk = nyq_eff * (nyq_eff + 1)
        do h = fpllims(1,1), fpllims(1,2)
            h_sq = h*h
            if( h_sq > nyq_disk ) cycle
            k_max_h = int(sqrt(real(nyq_disk - h_sq)))
            k_lo    = max(fpllims(2,1), -k_max_h)
            k_hi    = min(0, min(fpllims(2,2), k_max_h))
            hp      = h * pf
            do k = k_lo, k_hi
                kp = k * pf
                fpl%cmplx_plane(hp,kp) = fpl%cmplx_plane(hp,kp) - a*mean_fpl%cmplx_plane(hp,kp)
            end do
        end do
    end subroutine subtract_mean_banded

end submodule simple_flex_probe_fit_estep
