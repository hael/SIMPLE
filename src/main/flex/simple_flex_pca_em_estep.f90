!@descr: flex_pca EM: the E-step passes (polar bank, formers, per-particle solve, batch insert, reduce, paired pass and its part reduce)
submodule (simple_flex_pca_em) simple_flex_pca_em_estep
use simple_matcher_3Drec,   only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch
use simple_flex_pca_plane_cache, only: plane_cache_in_use
use simple_flex_pca_planes, only: planes_batch_load
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

    !> Polar E-step bank for one fit: grid geometry + pose-fixed direction assignment (once per
    !! stage), then the per-iteration shared-direction bank + ring Gram tables at the fit's
    !! current rank, restricted to the directions the fit's current window touches.
    module subroutine fit_polar_bank_build( params, build, fit, mean_rec, fpl1, nthr, l_dev_geom )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(probe_fit_t),   intent(inout) :: fit
        type(reconstructor), intent(inout) :: mean_rec
        type(fplane_type),   intent(in)    :: fpl1
        integer,             intent(in)    :: nthr
        !> device polar E-step requests its per-stage ring geometry (driver-owned selector)
        logical,             intent(in)    :: l_dev_geom
        type(oris) :: dirs_es
        type(ori)  :: o_es
        real,    allocatable :: rmatp_es(:,:,:), nrmp_es(:,:)
        real     :: ca1, sa1
        integer  :: i, q, r, ir, id_es, ithr, jx_es, kx_es
                    if( .not. fit%l_pol_grid )then
                        ! grid geometry from the first prepped plane: the identical derivation
                        ! embed_accumulate_polar uses, so polar embed and polar E-step share
                        ! band/quadrature conventions. No noise rings: sig2 arrives as sig2_eff.
                        fit%ph0_es  = lbound(fpl1%cmplx_plane,1)
                        fit%pk0_es  = lbound(fpl1%cmplx_plane,2)
                        fit%hlo_es  = ceil_div (lbound(fpl1%cmplx_plane,1), OSMPL_PAD_FAC)
                        fit%hhi_es  = floor_div(ubound(fpl1%cmplx_plane,1), OSMPL_PAD_FAC)
                        fit%klo_es  = ceil_div (lbound(fpl1%cmplx_plane,2), OSMPL_PAD_FAC)
                        fit%nyqr_es = mean_rec%get_lfny(1)
                        fit%nyqb_es = fit%nyqr_es
                        if( fpl1%nyq > 0 ) fit%nyqb_es = min(fit%nyqb_es, max(1, fpl1%nyq / OSMPL_PAD_FAC))
                        ! hybrid split point: auto at 0.72*band -- the measured knee of the
                        ! real-data ladder (10049 gate, band 11: rhyb 6 -> b err 8.9%, min z
                        ! corr 0.949, lambda head 153/282/326; rhyb 8 -> b err 3.1%, min z
                        ! corr 0.990, lambda 166/273/394 vs Cartesian 166/264/382). Env
                        ! override; explicit 0 = pure rings (the pre-hybrid baseline).
                        fit%rhyb_es = nint(0.72*real(fit%nyqb_es))
                        if( fit%rhyb_req > 0 ) fit%rhyb_es = fit%rhyb_req
                        if( fit%l_rhyb_off )    fit%rhyb_es = 0
                        fit%rhyb_es   = max(0, min(fit%rhyb_es, fit%nyqb_es-1))
                        fit%l_pol_hyb = fit%rhyb_es > 0
                        if( fit%l_pol_hyb )then
                            call polar_grid_build(fit%pg_es, fit%rhyb_es+1, fit%nyqb_es, fit%nyqb_es+1, fit%nyqb_es, &
                                &fit%hlo_es, fit%hhi_es, fit%klo_es, fit%ph0_es, fit%pk0_es, ang_osamp=fit%osamp_pol, &
                                &gate_lo=fit%rhyb_es*(fit%rhyb_es+1))
                            ! exact-part lattice positions, cov_herm_inner's half-plane rule
                            ! (k<=0; on the k=0 line only h<=0), shells 0..rhyb by the nint
                            ! convention: h^2+k^2 <= rhyb*(rhyb+1). Raster order k-outer.
                            fit%npos_es = 0
                            do kx_es = -fit%rhyb_es, 0
                                do jx_es = -fit%rhyb_es, merge(0, fit%rhyb_es, kx_es == 0)
                                    if( jx_es*jx_es + kx_es*kx_es > fit%rhyb_es*(fit%rhyb_es+1) ) cycle
                                    fit%npos_es = fit%npos_es + 1
                                end do
                            end do
                            allocate(fit%hex_es(fit%npos_es), fit%kex_es(fit%npos_es))
                            fit%npos_es = 0
                            do kx_es = -fit%rhyb_es, 0
                                do jx_es = -fit%rhyb_es, merge(0, fit%rhyb_es, kx_es == 0)
                                    if( jx_es*jx_es + kx_es*kx_es > fit%rhyb_es*(fit%rhyb_es+1) ) cycle
                                    fit%npos_es = fit%npos_es + 1
                                    fit%hex_es(fit%npos_es) = jx_es
                                    fit%kex_es(fit%npos_es) = kx_es
                                end do
                            end do
                            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA POLAR HYBRID: exact &
                                &Cartesian statistics for shells 0..',fit%rhyb_es,' (',fit%npos_es,' lattice &
                                &points/particle incl. DC), rings ',fit%rhyb_es+1,'..',fit%nyqb_es
                            if( fit%osamp_pol > 1 ) write(logfhandle,'(A,I0)') &
                                &'>>> FLEX_PCA POLAR OSAMP: ring angular oversampling x', fit%osamp_pol
                            call flush(logfhandle)
                        else
                            call polar_grid_build(fit%pg_es, 1, fit%nyqb_es, fit%nyqb_es+1, fit%nyqb_es, fit%hlo_es, fit%hhi_es, &
                                &fit%klo_es, fit%ph0_es, fit%pk0_es, ang_osamp=fit%osamp_pol)
                            if( fit%osamp_pol > 1 )then
                                write(logfhandle,'(A,I0)') &
                                    &'>>> FLEX_PCA POLAR OSAMP: ring angular oversampling x', fit%osamp_pol
                                call flush(logfhandle)
                            endif
                        endif
                        fit%nsamp_es = fit%pg_es%nsamp; fit%nsamp2_es = 2*fit%nsamp_es; fit%nk_es = fit%pg_es%nk
                        ! shared direction table: same refspiral + count derivation as the polar embed
                        fit%ndir_es = cov_polar_ndir(fit%npp)
                        call dirs_es%new(fit%ndir_es, is_ptcl=.false.)
                        call build%pgrpsyms%build_refspiral(dirs_es)
                        allocate(fit%rmatb_es(3,3,fit%ndir_es), fit%nrmb_es(3,fit%ndir_es))
                        do id_es = 1, fit%ndir_es
                            fit%rmatb_es(:,:,id_es) = dirs_es%get_mat(id_es)
                            fit%nrmb_es(:,id_es)    = fit%rmatb_es(3,:,id_es)
                        end do
                        call dirs_es%kill
                        ! pose-fixed per-particle (direction, in-plane) assignment, once per stage
                        allocate(rmatp_es(3,3,fit%npp), nrmp_es(3,fit%npp), fit%dir_es(fit%npp), fit%cae(fit%npp), fit%sae(fit%npp))
                        do i = 1, fit%npp
                            call build%spproj_field%get_ori(fit%ppinds(i), o_es)
                            rmatp_es(:,:,i) = o_es%get_mat()
                            nrmp_es(:,i)    = rmatp_es(3,:,i)
                        end do
                        call o_es%kill
                        call polar_assign_directions(nrmp_es, fit%npp, fit%nrmb_es, fit%ndir_es, fit%dir_es)
                        do i = 1, fit%npp
                            call polar_relative_inplane(rmatp_es(:,:,i), fit%rmatb_es(:,:,fit%dir_es(i)), ca1, sa1)
                            fit%cae(i) = ca1; fit%sae(i) = sa1
                        end do
                        deallocate(rmatp_es, nrmp_es)
                        allocate(fit%dused_es(fit%ndir_es))
                        ! device ring geometry, once per stage (the grid is pose- and band-fixed)
                        if( l_dev_geom ) call flex_gpu_poles_begin_f(fit%pg_es%rad, fit%pg_es%cs, &
                            &fit%pg_es%sn, fit%pg_es%sqwq, fit%pg_es%rbeg, fit%pg_es%rend, fit%nsamp_es, fit%nk_es)
                        fit%l_pol_grid = .true.
                        write(logfhandle,'(A,I0,A,I0,A,I0,A,F8.1,A,F8.1,A)') &
                            &'>>> FLEX_PCA POLAR ESTEP BANK: ',fit%ncomp+1,' volumes x ',fit%ndir_es, &
                            &' directions x ',fit%nsamp_es,' ring samples = ', &
                            &4.d0*real(fit%nsamp2_es,dp)*real(fit%ncomp+1,dp)*real(fit%ndir_es,dp)/1.d6, &
                            &' MB (+ ring tables ', &
                            &8.d0*real(fit%ncomp*fit%ncomp+fit%ncomp+1,dp)*real(fit%nk_es,dp)*real(fit%ndir_es,dp)/1.d6,' MB)'
                        call flush(logfhandle)
                    endif
                    ! (re)allocate at this iteration's rank (ncomp can change between iterations)
                    if( allocated(fit%UsallE) )then
                        if( size(fit%UsallE,2) /= fit%ncomp+1 ) deallocate(fit%UsallE, fit%CfE, fit%Cm0E, fit%c00E, &
                            &fit%UbankE, fit%CspE, fit%xws_es, fit%wr_es, fit%wrd_es, fit%Reb_es)
                    endif
                    if( .not. allocated(fit%UsallE) )then
                        allocate(fit%UsallE(fit%nsamp2_es,0:fit%ncomp,fit%ndir_es), fit%CfE(fit%ncomp*fit%ncomp,fit%nk_es,fit%ndir_es), &
                            &fit%Cm0E(fit%ncomp,fit%nk_es,fit%ndir_es), fit%c00E(fit%nk_es,fit%ndir_es))
                        allocate(fit%UbankE(fit%nsamp_es,0:fit%ncomp,nthr), fit%CspE(0:fit%ncomp,0:fit%ncomp,nthr))
                        allocate(fit%xws_es(fit%nsamp2_es,nthr), fit%wr_es(fit%nk_es,nthr), fit%wrd_es(fit%nk_es,nthr), &
                            &fit%Reb_es(0:fit%ncomp,nthr))
                    endif
                    ! only the directions this iteration's window touches
                    fit%dused_es = .false.
                    do i = 1, fit%npp
                        if( fit%dir_es(i) > 0 ) fit%dused_es(fit%dir_es(i)) = .true.
                    end do
                    !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
                    !$omp& private(id_es,ithr,q,r,ir)
                    do id_es = 1, fit%ndir_es
                        if( .not. fit%dused_es(id_es) ) cycle
                        ithr = omp_get_thread_num() + 1
                        call polar_project_recs(mean_rec, fit%basis_recs, fit%ncomp, fit%rmatb_es(:,:,id_es), &
                            &fit%pg_es, fit%UbankE(:,:,ithr))
                        do q = 0, fit%ncomp
                            do r = 1, fit%nsamp_es
                                fit%UsallE(2*r-1,q,id_es) = fit%pg_es%sqwq(r)*real (fit%UbankE(r,q,ithr))
                                fit%UsallE(2*r,  q,id_es) = fit%pg_es%sqwq(r)*aimag(fit%UbankE(r,q,ithr))
                            end do
                        end do
                        do ir = 1, fit%nk_es
                            call polar_ring_gram(fit%UsallE(1,0,id_es), fit%nsamp2_es, fit%ncomp, fit%pg_es%rbeg(ir), &
                                &fit%pg_es%rend(ir)-fit%pg_es%rbeg(ir)+1, fit%CspE(0,0,ithr), fit%CfE(1,ir,id_es), &
                                &fit%Cm0E(1,ir,id_es))
                            fit%c00E(ir,id_es) = polar_ring_selfpower(fit%UsallE(1,0,id_es), fit%nsamp2_es, &
                                &fit%pg_es%rbeg(ir), fit%pg_es%rend(ir)-fit%pg_es%rbeg(ir)+1)
                        end do
                    end do
                    !$omp end parallel do
    end subroutine fit_polar_bank_build

    !> Production polar shared-direction former for one particle: banded mean projection,
    !! fused polar sampling, bank GEMVs, hybrid low-k exact statistics, contrast fit.
    module subroutine fit_estep_former_polar( fit, mean_rec, o, fpl, row, ithr, a, aa, e_mm, myv )
        type(probe_fit_t),   intent(inout) :: fit
        type(reconstructor), intent(inout) :: mean_rec
        class(ori),          intent(inout) :: o
        type(fplane_type),   intent(inout) :: fpl
        integer,             intent(in)    :: row, ithr
        real(dp),            intent(out)   :: a, aa, e_mm, myv
        integer  :: idp_es, q
        real     :: taz_es
        real(dp) :: twp0, twp1
        twp0 = omp_get_wtime()
                        idp_es  = fit%dir_es(row)
                            ! mean only, in Cartesian: the M-step backprojects y - a*(T mu) and
                            ! there is no polar->volume adjoint. Banded variant: identical
                            ! interpolation, none of project_fplane's per-call full-plane
                            ! zero-fill + ctfsq/transfer copies (measured as the bulk of the
                            ! polar project bucket at the native padded plane).
                            call project_fplane_mean_banded(mean_rec, o, fpl, &
                                &fit%mean_fpl(ithr))
                        ! polar-sample the prepped data plane ONCE at (bank direction, relative
                        ! in-plane angle): CTF amplitude, shift phase and per-shell whitening all
                        ! ride in from the same cmplx/transfer planes the Cartesian former reads.
                        ! Fused sampler: one KB geometry per ring sample shared by both plane
                        ! gathers, packed output written in place, no per-call allocations --
                        ! bit-identical statistics (see polar_sample_particle_fused).
                        call polar_sample_particle_fused(fpl%cmplx_plane, fpl%transfer_plane, &
                            &fit%pg_es, fit%cae(row), fit%sae(row), fit%xws_es(:,ithr), fit%wr_es(:,ithr), taz_es)
                        fit%wrd_es(:,ithr) = real(fit%wr_es(:,ithr), dp)
                        ! b and the mean row, exact per-sample CTF: one GEMV against the bank
                        call sgemv('T', fit%nsamp2_es, fit%ncomp+1, 1.0, fit%UsallE(1,0,idp_es), fit%nsamp2_es, &
                            &fit%xws_es(1,ithr), 1, 0.0, fit%Reb_es(0,ithr), 1)
                        ! G, c, e_mm by the radial factorisation: ring Grams x per-ring mean |T|^2
                        call dgemv('N', fit%ncomp*fit%ncomp, fit%nk_es, 1.d0, fit%CfE(1,1,idp_es), fit%ncomp*fit%ncomp, &
                            &fit%wrd_es(1,ithr), 1, 0.d0, fit%Gth(1,1,ithr), 1)
                        call dgemv('N', fit%ncomp, fit%nk_es, 1.d0, fit%Cm0E(1,1,idp_es), fit%ncomp, &
                            &fit%wrd_es(1,ithr), 1, 0.d0, fit%cth(1,ithr), 1)
                        e_mm = dot_product(fit%c00E(:,idp_es), fit%wrd_es(:,ithr))
                        myv  = real(fit%Reb_es(0,ithr), dp)
                        do q = 1, fit%ncomp
                            fit%bth(q,ithr) = real(fit%Reb_es(q,ithr), dp)
                        end do
                        ! hybrid: the low-k shells enter as exact Cartesian statistics
                        if( fit%l_pol_hyb ) call polar_hybrid_exact_accum(mean_rec, fit%basis_recs, &
                            &fit%ncomp, o, fpl, fit%hex_es, fit%kex_es, fit%npos_es, &
                            &fit%Gth(:,:,ithr), fit%bth(:,ithr), fit%cth(:,ithr), e_mm, myv)
                        a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
                        aa   = a*a
                        twp1 = omp_get_wtime()
                        fit%sec_proj_thr(ithr) = fit%sec_proj_thr(ithr) + (twp1 - twp0)
    end subroutine fit_estep_former_polar

    !> Production Cartesian former for one particle: mean + basis central sections, exact
    !! Hermitian inner products.
    module subroutine fit_estep_former_cart( fit, mean_rec, o, fpl, ithr, a, aa, e_mm, myv )
        type(probe_fit_t),   intent(inout) :: fit
        type(reconstructor), intent(inout) :: mean_rec
        class(ori),          intent(inout) :: o
        type(fplane_type),   intent(inout) :: fpl
        integer,             intent(in)    :: ithr
        real(dp),            intent(out)   :: a, aa, e_mm, myv
        integer  :: q, r
        real(dp) :: twp0, twp1
        twp0 = omp_get_wtime()
                    call project_fplanes_mean_basis(mean_rec, fit%basis_recs, o, fpl, &
                        &fit%mean_fpl(ithr), fit%basis_fpls(:,ithr), apply_ctf_amp=.true.)
                    twp1 = omp_get_wtime()
                    fit%sec_proj_thr(ithr) = fit%sec_proj_thr(ithr) + (twp1 - twp0)
                    e_mm = real(cov_herm_inner(fit%mean_fpl(ithr), fit%mean_fpl(ithr)), dp)
                    myv  = real(cov_herm_inner(fit%mean_fpl(ithr), fpl), dp)
                    a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
                    aa   = a*a
                    do q = 1, fit%ncomp
                        fit%bth(q,ithr) = real(cov_herm_inner(fit%basis_fpls(q,ithr), fpl), dp)
                        fit%cth(q,ithr) = real(cov_herm_inner(fit%basis_fpls(q,ithr), fit%mean_fpl(ithr)), dp)
                        do r = q, fit%ncomp
                            fit%Gth(q,r,ithr) = real(cov_herm_inner(fit%basis_fpls(q,ithr), fit%basis_fpls(r,ithr)), dp)
                            fit%Gth(r,q,ithr) = fit%Gth(q,r,ithr)
                        end do
                    end do
    end subroutine fit_estep_former_cart

    !> Shared per-particle tail of every CPU former: posterior solve (plain or mixture), latent
    !! and density batch rows, Gamma/likelihood accumulation, and the in-place mean-subtracted
    !! residual the M-step backprojects.
    module subroutine fit_estep_solve_stats( fit, fpl, i, row, ithr, a, aa, e_mm, myv )
        type(probe_fit_t), intent(inout) :: fit
        type(fplane_type), intent(inout) :: fpl
        integer,           intent(in)    :: i, row, ithr
        real(dp),          intent(inout) :: a, aa, e_mm, myv
        integer  :: q, r
        logical  :: lok
        real(dp) :: ldA, qml, nll_mix_add, twp1, twp2
        twp1 = omp_get_wtime()
                    ! Posterior precision A = (a^2/sig2) G + Gamma^-1. The whole normal system is already
                    ! scaled by 1/sig2, so Cov[z|y] = A^-1 exactly -- no further sig2 factor.
                    if( fit%l_mix_used )then
                        ! ---- MCFA E-step via the ONE shared solver (probe_solve_mix) ----
                        call probe_solve_mix(fit%ncomp, fit%kmix, fit%Gth(:,:,ithr), fit%bth(:,ithr), &
                            &fit%cth(:,ithr), myv, e_mm, fit%nml_plain, a, fit%sig2, fit%mix_Ominv, fit%mix_Omxi, &
                            &fit%mix_lpi, fit%mix_xiOx, fit%zth(:,ithr), fit%Ainvth(:,:,ithr), fit%dens(:,:,i), &
                            &ldA, lok, nll_mix_add, fit%mxa_sr(:,ithr), fit%mxa_sm(:,:,ithr), &
                            &fit%mxa_smm(:,:,:,ithr), fit%mxa_sainv(:,:,ithr))
                        if( lok ) fit%nll_thr(ithr) = fit%nll_thr(ithr) + nll_mix_add
                    else
                        call probe_solve_ecm(fit%ncomp, fit%Gth(:,:,ithr), fit%bth(:,ithr), fit%cth(:,ithr), &
                            &myv, e_mm, fit%prior, fit%sig2, fit%nml_plain, a, fit%zth(:,ithr), &
                            &fit%Ainvth(:,:,ithr), ldA, lok, qml)
                        if( lok ) fit%nll_thr(ithr) = fit%nll_thr(ithr) + ldA - qml
                    endif
                    aa = a*a
                    fit%z(row,:)          = fit%zth(:,ithr)
                    fit%zbatch(:,i)       = fit%zth(:,ithr)
                    ! EM sufficient statistic E[z z'|y] = z z' + Cov[z|y]. BOTH the coupled M-step normal
                    ! matrix and the Gamma update below need it. Dropping Cov underestimates Gamma, which
                    ! tightens the prior, which shrinks z further: the bias compounds across iterations.
                    ! Under the mixture E[zz'|y] = A^-1 + sum_k r_k m_k m_k', NOT A^-1 + E[z]E[z]'
                    ! -- the between-component spread is real posterior variance; probe_solve_mix
                    ! has already written dens(:,:,i) in that case.
                    if( .not. fit%l_mix_used )then
                        do r = 1, fit%ncomp
                            do q = 1, fit%ncomp
                                fit%dens(q,r,i) = fit%zth(q,ithr)*fit%zth(r,ithr) + fit%Ainvth(q,r,ithr)
                            end do
                        end do
                    endif
                    do q = 1, fit%ncomp
                        fit%gam_thr(q,ithr) = fit%gam_thr(q,ithr) + fit%dens(q,q,i)
                    end do
                    if( fit%l_probe_mls )then
                        fit%zbatch(:,i) = a*fit%zbatch(:,i)
                        fit%dens(:,:,i) = aa*fit%dens(:,:,i)
                    endif
                    fit%nval_thr(ithr)    = fit%nval_thr(ithr) + 1
                    fit%valid(i)          = .true.
                    fit%gam_dbg(1,ithr) = fit%gam_dbg(1,ithr) + sum([(fit%Gth(q,q,ithr), q=1,fit%ncomp)])
                    fit%gam_dbg(2,ithr) = fit%gam_dbg(2,ithr) + dot_product(fit%bth(:,ithr), fit%bth(:,ithr))
                    fit%gam_dbg(3,ithr) = fit%gam_dbg(3,ithr) + dot_product(fit%cth(:,ithr), fit%cth(:,ithr))
                    fit%gam_dbg(4,ithr) = fit%gam_dbg(4,ithr) + a
                    ! residual observation r_i = y - a*(T mu) in place (transfer/ctfsq intact for backprojection)
                    if( fit%l_pol_es )then
                        ! banded: the mean plane is zero outside the working disc (both formers
                        ! write the same disc), so the full-array statement only rewrote
                        ! unchanged values there -- at the native padded lattice that traffic
                        ! was most of the polar path's gram+solve bucket
                        call subtract_mean_banded(fpl, fit%mean_fpl(ithr), real(a), fit%nyqr_es)
                    else
                        fpl%cmplx_plane = fpl%cmplx_plane - real(a)*fit%mean_fpl(ithr)%cmplx_plane
                    endif
                    twp2 = omp_get_wtime()
                    fit%sec_gram_thr(ithr) = fit%sec_gram_thr(ithr) + (twp2 - twp1)
    end subroutine fit_estep_solve_stats

    !> Halfset masks + the CPU coupled M-step insertion for one batch of one fit.
    module subroutine fit_batch_insert( build, fit, orientations, fpls, eo, batchsz )
        type(builder),     intent(inout) :: build
        type(probe_fit_t), intent(inout) :: fit
        type(ori),         intent(inout) :: orientations(:)
        type(fplane_type), intent(inout) :: fpls(:)
        integer,           intent(in)    :: eo(:), batchsz
        integer :: i
        do i = 1, batchsz
            fit%valid_e(i) = fit%valid(i) .and. eo(i)==0
            fit%valid_o(i) = fit%valid(i) .and. eo(i)==1
        end do
        call insert_planes_oversamp_coupled_batch_scaled(fit%Yeven, fit%rho_e, build%pgrpsyms, &
            &orientations(:batchsz), fpls(:batchsz), fit%zbatch(:,:batchsz), fit%dens(:,:,:batchsz), &
            &fit%valid_e(:batchsz), batchsz)
        call insert_planes_oversamp_coupled_batch_scaled(fit%Yodd, fit%rho_o, build%pgrpsyms, &
            &orientations(:batchsz), fpls(:batchsz), fit%zbatch(:,:batchsz), fit%dens(:,:,:batchsz), &
            &fit%valid_o(:batchsz), batchsz)
        ! rec_backend=pcg: the pair-weighted Gram kernels of the same batch at doubled coordinates
        ! (accumulators live only on inserting processes; the master reduces the packed sums)
        if( fit%l_pcg )then
            if( .not. allocated(fit%kacc_e) ) call fit%pcg%alloc_accum(fit%kacc_e)
            if( .not. allocated(fit%kacc_o) ) call fit%pcg%alloc_accum(fit%kacc_o)
            if( .not. allocated(fit%racc_e) ) call fit%pcg%alloc_rhs_accum(fit%racc_e)
            if( .not. allocated(fit%racc_o) ) call fit%pcg%alloc_rhs_accum(fit%racc_o)
            call fit%pcg%accumulate(fit%kacc_e, build%pgrpsyms, orientations(:batchsz), fpls(:batchsz), &
                &fit%dens(:,:,:batchsz), fit%valid_e(:batchsz), batchsz)
            call fit%pcg%accumulate(fit%kacc_o, build%pgrpsyms, orientations(:batchsz), fpls(:batchsz), &
                &fit%dens(:,:,:batchsz), fit%valid_o(:batchsz), batchsz)
            ! the same batch's right-hand sides at doubled coordinates (the exact adjoint for the solve)
            call fit%pcg%accumulate_rhs(fit%racc_e, build%pgrpsyms, orientations(:batchsz), fpls(:batchsz), &
                &fit%zbatch(:,:batchsz), fit%valid_e(:batchsz), batchsz)
            call fit%pcg%accumulate_rhs(fit%racc_o, build%pgrpsyms, orientations(:batchsz), fpls(:batchsz), &
                &fit%zbatch(:,:batchsz), fit%valid_o(:batchsz), batchsz)
        endif
    end subroutine fit_batch_insert

    !> In-process end-of-iteration reductions (thread sums) and the MCFA mixture M-step or its
    !! initialisation. Runs on the shared-memory master and inside workers; the distributed master
    !! takes its sums from reduce_probe_parts.
    module subroutine fit_iter_reduce( fit, it_eff, nthr , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        type(probe_fit_t), intent(inout) :: fit
        integer,           intent(in)    :: it_eff, nthr
        integer :: q, kk2
            ! reduce the EM Gamma accumulator before ncomp is replaced below; Gamma travels between
            ! parts as a sum and is divided by the reduced nval
        fit%nval = sum(fit%nval_thr)
        ! E-step statistics summary: the polar and Cartesian formers must agree on G, b, c
        if( fit%nval > 0 )then
            write(logfhandle,'(A,I0,A,ES12.4,A,ES12.4,A,ES12.4,A,F7.4)') &
                &'>>> FLEX_PCA PROBE ESTAT it=',it_eff,'  <trG>=',sum(fit%gam_dbg(1,:))/real(fit%nval,dp), &
                &'  <b.b>=',sum(fit%gam_dbg(2,:))/real(fit%nval,dp), &
                &'  <c.c>=',sum(fit%gam_dbg(3,:))/real(fit%nval,dp), &
                &'  <a>=',real(sum(fit%gam_dbg(4,:))/real(fit%nval,dp))
            call flush(logfhandle)
        endif
        ! in-process reduction of the likelihood accumulator (distributed runs sum it in
        ! reduce_probe_parts)
        fit%nll_tot = sum(fit%nll_thr)
        do q = 1, fit%ncomp
            fit%gam_sum(q) = sum(fit%gam_thr(q,:))
        end do
        ! MCFA: mixture M-step, or its (re)initialisation from the current latents; runs on the
        ! master with this iteration's accumulators complete
        if( fit%l_mix_req .and. .not. rounds%is_worker() )then
            if( fit%l_mix_active )then
                block
                    real(dp), allocatable :: rr_sr(:), rr_sm(:,:), rr_smm(:,:,:), rr_sai(:,:)
                    integer  :: tt, kk3
                    allocate(rr_sr(fit%kmix), rr_sm(fit%ncomp,fit%kmix), &
                        &rr_smm(fit%ncomp,fit%ncomp,fit%kmix), rr_sai(fit%ncomp,fit%ncomp))
                    rr_sr = 0.d0; rr_sm = 0.d0; rr_smm = 0.d0; rr_sai = 0.d0
                    do tt = 1, nthr
                        rr_sr  = rr_sr  + fit%mxa_sr(:,tt)
                        rr_sm  = rr_sm  + fit%mxa_sm(:,:,tt)
                        rr_sai = rr_sai + fit%mxa_sainv(:,:,tt)
                        do kk3 = 1, fit%kmix
                            rr_smm(:,:,kk3) = rr_smm(:,:,kk3) + fit%mxa_smm(:,:,kk3,tt)
                        end do
                    end do
                    if( allocated(fit%dm_sr) .and. rounds%nparts() > 1 )then
                        ! distributed: the reduce already summed every part's thread-sums
                        call mcfa_mstep(fit%ncomp, fit%kmix, fit%nval, fit%dm_sr, fit%dm_sm, fit%dm_smm, fit%dm_sai, &
                            &fit%kmix == 1, fit%mix_pi, fit%mix_xi, fit%mix_Om, fit%mix_Ominv, fit%ldOm_mix)
                    else
                        call mcfa_mstep(fit%ncomp, fit%kmix, fit%nval, rr_sr, rr_sm, rr_smm, rr_sai, &
                            &fit%kmix == 1, fit%mix_pi, fit%mix_xi, fit%mix_Om, fit%mix_Ominv, fit%ldOm_mix)
                    endif
                    deallocate(rr_sr, rr_sm, rr_smm, rr_sai)
                end block
            else if( it_eff >= fit%n_mix_warm )then
                allocate(fit%mix_xi(fit%ncomp,fit%kmix), fit%mix_Om(fit%ncomp,fit%ncomp), fit%mix_Ominv(fit%ncomp,fit%ncomp), &
                    &fit%mix_pi(fit%kmix), fit%mix_Omxi(fit%ncomp,fit%kmix), fit%mix_xiOx(fit%kmix), fit%mix_lpi(fit%kmix))
                if( rounds%nparts() > 1 .and. fit%dm_nz > 0 )then
                    ! distributed: seed from the pooled per-part latent subsample
                    call mcfa_init(fit%dm_z(:fit%dm_nz,:), fit%dm_nz, fit%ncomp, fit%kmix, fit%gam_sum, fit%nval, &
                        &fit%mix_xi, fit%mix_pi, fit%mix_Om)
                else
                    call mcfa_init(fit%z, size(fit%z,1), fit%ncomp, fit%kmix, fit%gam_sum, fit%nval, fit%mix_xi, fit%mix_pi, fit%mix_Om)
                endif
                call mcfa_condition(fit%ncomp, fit%kmix == 1, fit%mix_Om, fit%mix_Ominv, fit%ldOm_mix)
                fit%l_mix_active = .true.
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA MIX initialised: K=',fit%kmix, &
                    &'  at iteration ',it_eff
            endif
            if( fit%l_mix_active )then
                do kk2 = 1, fit%kmix
                    fit%mix_Omxi(:,kk2) = matmul(fit%mix_Ominv, fit%mix_xi(:,kk2))
                    fit%mix_xiOx(kk2)   = dot_product(fit%mix_xi(:,kk2), fit%mix_Omxi(:,kk2))
                    fit%mix_lpi(kk2)    = log(max(fit%mix_pi(kk2), 1.d-12))
                end do
                write(logfhandle,'(A,I0,A,ES11.3,A,F7.4,A,F8.3)') '>>> FLEX_PCA MIX it=',it_eff, &
                    &'  logdetOm=',fit%ldOm_mix,'  pi_max=',maxval(fit%mix_pi), &
                    &'  |xi|_max=',sqrt(maxval(sum(fit%mix_xi**2, dim=1)))
                call flush(logfhandle)
            endif
        endif
        ! rec_backend=pcg: fold this process's kernel accumulators into the packed sums (adds; the
        ! distributed master has no accumulators and receives the parts' packed sums instead)
        if( fit%l_pcg )then
            if( allocated(fit%kacc_e) ) call fit%pcg%fold_accum(fit%kacc_e, fit%kpk_e)
            if( allocated(fit%kacc_o) ) call fit%pcg%fold_accum(fit%kacc_o, fit%kpk_o)
            if( allocated(fit%racc_e) ) call fit%pcg%fold_rhs(fit%racc_e, fit%rpk_e)
            if( allocated(fit%racc_o) ) call fit%pcg%fold_rhs(fit%racc_o, fit%rpk_o)
        endif
    end subroutine fit_iter_reduce

    !> Paired E-step accumulate pass over one iteration's merged read list: the shared-memory master
    !! runs it over the full halves, each distributed worker over its fromp/top shard. Owns its
    !! read/prep buffers; the caller owns the fits' accumulators (fit_iter_begin / fit_iter_reduce).
    module subroutine paired_estep_pass( params, build, fits, it_eff, nthr )
        class(parameters), intent(inout) :: params
        type(builder),     intent(inout) :: build
        type(probe_fit_t), intent(inout) :: fits(2)
        integer,           intent(in)    :: it_eff, nthr
        type(fplane_type), allocatable :: fpls(:)
        type(ori),         allocatable :: orientations(:)
        integer,           allocatable :: eo(:)
        integer,           allocatable :: mrg_pinds(:), mrg_owner(:), mrg_row(:)
        integer  :: f, i, ithr, nmrg, ibatch, batchlims(2), batchsz, ia, ib, row
        real(dp) :: a, aa, e_mm, myv
        real(timer_int_kind) :: sec_read, sec_prep, sec_estep, sec_ins
        logical  :: l_pcache
        integer(timer_int_kind) :: t_sec, t_bank
        sec_read = 0.; sec_prep = 0.; sec_estep = 0.; sec_ins = 0.
        allocate(orientations(MAXIMGBATCHSZ), eo(MAXIMGBATCHSZ))
        ! ---- merged read list: linear two-pointer merge of the two fits' windows, sorted
        ! by project row (each fit's ppinds is ascending) ----
        nmrg = fits(1)%npp + fits(2)%npp
        allocate(mrg_pinds(nmrg), mrg_owner(nmrg), mrg_row(nmrg))
        ia = 1
        ib = 1
        i  = 0
        do while( ia <= fits(1)%npp .or. ib <= fits(2)%npp )
            i = i + 1
            if( ib > fits(2)%npp )then
                mrg_pinds(i) = fits(1)%ppinds(ia); mrg_owner(i) = 1; mrg_row(i) = ia; ia = ia + 1
            elseif( ia > fits(1)%npp )then
                mrg_pinds(i) = fits(2)%ppinds(ib); mrg_owner(i) = 2; mrg_row(i) = ib; ib = ib + 1
            elseif( fits(1)%ppinds(ia) <= fits(2)%ppinds(ib) )then
                mrg_pinds(i) = fits(1)%ppinds(ia); mrg_owner(i) = 1; mrg_row(i) = ia; ia = ia + 1
            else
                mrg_pinds(i) = fits(2)%ppinds(ib); mrg_owner(i) = 2; mrg_row(i) = ib; ib = ib + 1
            endif
        end do
        if( i /= nmrg ) THROW_HARD('paired merged read list lost rows')
        ! downscaled-particle cache (cache=yes): a cache-served batch is read at box_crop into
        ! cropped planes (init_rec cropped=, prepimgbatch at box_crop) and prepped with
        ! cached=.true. (already noise-normalised when the entry was written)
        l_pcache = plane_cache_in_use(params, build)
        call init_rec(params, build, MAXIMGBATCHSZ, fpls, cropped=l_pcache)
        if( l_pcache )then
            call prepimgbatch(params, build, MAXIMGBATCHSZ, box=params%box_crop, smpd=params%smpd_crop)
            write(logfhandle,'(A)') '>>> FLEX_PCA E-STEP: particles served from the downscaled cache'
            call flush(logfhandle)
        else
            call prepimgbatch(params, build, MAXIMGBATCHSZ)
        endif
        do ibatch = 1, nmrg, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nmrg, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call planes_batch_load(params, build, nmrg, mrg_pinds, batchlims, fpls, &
                &cov_image_mask_radius(params), l_pcache, sec_read, sec_prep)
            do i = 1, batchsz
                call build%spproj_field%get_ori(mrg_pinds(batchlims(1)+i-1), orientations(i))
                eo(i) = build%spproj_field%get_eo(mrg_pinds(batchlims(1)+i-1))
            end do
            ! per-fit polar bank, built once per iteration at the first prepped batch
            do f = 1, 2
                if( fits(f)%l_pol_es .and. .not. fits(f)%l_pol_bank_it )then
                    t_bank = tic()
                    call fit_polar_bank_build(params, build, fits(f), fits(f)%mean_rec, &
                        &fpls(1), nthr, .false.)
                    fits(f)%sec_bank      = fits(f)%sec_bank + real(toc(t_bank))
                    fits(f)%l_pol_bank_it = .true.
                    write(logfhandle,'(A,I0,A,A,A,I0,A,I0,A,F7.1)') &
                        &'>>> FLEX_PCA POLAR ESTEP BANK it=', it_eff, '  fit=', &
                        &merge('A','B',f==1), ' directions built=', count(fits(f)%dused_es), &
                        &' of ', fits(f)%ndir_es, '  build seconds=', fits(f)%sec_bank
                    call flush(logfhandle)
                endif
                fits(f)%valid(:batchsz)    = .false.
                fits(f)%zbatch(:,:batchsz) = 0.d0
                fits(f)%dens(:,:,:batchsz) = 0.d0
            end do
            t_sec = tic()
            !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
            !$omp& private(i,f,row,ithr,a,aa,e_mm,myv)
            do i = 1, batchsz
                if( orientations(i)%isstatezero() ) cycle
                ithr = omp_get_thread_num() + 1
                f    = mrg_owner(batchlims(1)+i-1)
                row  = mrg_row(batchlims(1)+i-1)
                if( fits(f)%l_pol_es )then
                    call fit_estep_former_polar(fits(f), fits(f)%mean_rec, orientations(i), &
                        &fpls(i), row, ithr, a, aa, e_mm, myv)
                else
                    call fit_estep_former_cart(fits(f), fits(f)%mean_rec, orientations(i), &
                        &fpls(i), ithr, a, aa, e_mm, myv)
                endif
                call fit_estep_solve_stats(fits(f), fpls(i), i, row, ithr, a, aa, e_mm, myv)
            end do
            !$omp end parallel do
            sec_estep = sec_estep + toc(t_sec)
            t_sec = tic()
            do f = 1, 2
                call fit_batch_insert(build, fits(f), orientations, fpls, eo, batchsz)
            end do
            sec_ins = sec_ins + toc(t_sec)
            if( mod(batchlims(2), 5*MAXIMGBATCHSZ) == 0 .or. batchlims(2) == nmrg )then
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA PAIRED PASS PARTICLES: ', &
                    &batchlims(2), ' / ', nmrg
                call flush(logfhandle)
            endif
        end do
        write(logfhandle,'(A,F7.1,A,F7.1,A,F7.1,A,F7.1)') '>>> FLEX_PCA PAIRED E-STEP SPLIT &
            &(seconds): read=', sec_read, '  prep=', sec_prep, '  project+solve=', sec_estep, &
            &'  insert=', sec_ins
        call flush(logfhandle)
        do i = 1, size(orientations)
            call orientations(i)%kill
        end do
        call cleanup_rec_buffers(build, fpls)
        deallocate(orientations, eo, mrg_pinds, mrg_owner, mrg_row)
    end subroutine paired_estep_pass

    !> Master-side v5 reduce: fold every worker's paired part into BOTH fits' accumulators,
    !! streaming -- one part resident at a time, per-fit blocks in fit order, file deleted after.
    module subroutine paired_reduce_parts_v5( params, fits , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters), intent(in)    :: params
        type(probe_fit_t), intent(inout) :: fits(2)
        !> keep equal to probe_subspace_iteration's MIX_ZSUB_MAX (per-part latent subsample)
        integer, parameter :: MIX_ZSUB_MAX5 = 2000
        complex, allocatable :: cme1(:,:,:,:), cmo1(:,:,:,:), cme2(:,:,:,:), cmo2(:,:,:,:)
        real,    allocatable :: rhe1(:,:,:,:), rhoo1(:,:,:,:), rhe2(:,:,:,:), rhoo2(:,:,:,:)
        type(string) :: fname
        integer :: ipart, funit, f, q
        integer(timer_int_kind) :: t_red
        ! per-thread private accumulators of the concurrent part reads
        t_red = tic()
        associate( fa => fits(1), fb => fits(2) )
        allocate(cme1(fa%es(1),fa%es(2),fa%es(3),fa%ncomp), cmo1(fa%es(1),fa%es(2),fa%es(3),fa%ncomp), source=(0.,0.))
        allocate(rhe1(fa%es(1),fa%es(2),fa%es(3),fa%ncomp), rhoo1(fa%es(1),fa%es(2),fa%es(3),fa%ncomp), source=0.)
        allocate(cme2(fb%es(1),fb%es(2),fb%es(3),fb%ncomp), cmo2(fb%es(1),fb%es(2),fb%es(3),fb%ncomp), source=(0.,0.))
        allocate(rhe2(fb%es(1),fb%es(2),fb%es(3),fb%ncomp), rhoo2(fb%es(1),fb%es(2),fb%es(3),fb%ncomp), source=0.)
        end associate
        do f = 1, 2
            if( fits(f)%l_mix_req )then
                ! a rank change between iterations leaves these at the old ncomp; the part files carry
                ! the new one, and buffers sized from stale arrays misalign every read that follows
                if( allocated(fits(f)%dm_sm) )then
                    if( size(fits(f)%dm_sm,1) /= fits(f)%ncomp .or. size(fits(f)%dm_sr) /= fits(f)%kmix )then
                        deallocate(fits(f)%dm_sr, fits(f)%dm_sm, fits(f)%dm_smm, fits(f)%dm_sai, fits(f)%dm_z)
                    endif
                endif
                if( .not. allocated(fits(f)%dm_sr) )then
                    allocate(fits(f)%dm_sr(fits(f)%kmix), fits(f)%dm_sm(fits(f)%ncomp,fits(f)%kmix), &
                        &fits(f)%dm_smm(fits(f)%ncomp,fits(f)%ncomp,fits(f)%kmix), &
                        &fits(f)%dm_sai(fits(f)%ncomp,fits(f)%ncomp))
                    allocate(fits(f)%dm_z(MIX_ZSUB_MAX5*rounds%nparts(), fits(f)%ncomp))
                endif
                fits(f)%dm_sr = 0.d0; fits(f)%dm_sm = 0.d0; fits(f)%dm_smm = 0.d0
                fits(f)%dm_sai = 0.d0; fits(f)%dm_nz = 0
            endif
        end do
        do ipart = 1, rounds%nparts()
            fname = flex_pca_part_fname('probe', ipart, params%numlen)
            call open_probe_part_v5_read(fname, 2, funit)
            do f = 1, 2
                if( fits(f)%l_mix_req )then
                    if( f == 1 )then
                        call fold_probe_part_v5_fit(funit, cme1, rhe1, cmo1, rhoo1, &
                            &fits(f)%rho_e, fits(f)%rho_o, fits(f)%gam_sum, fits(f)%nll_tot, &
                            &fits(f)%nval, fits(f)%ncomp, fits(f)%kpk_e, fits(f)%kpk_o, fits(f)%rpk_e, fits(f)%rpk_o, fits(f)%dm_sr, fits(f)%dm_sm, &
                            &fits(f)%dm_smm, fits(f)%dm_sai, fits(f)%dm_z, fits(f)%dm_nz)
                    else
                        call fold_probe_part_v5_fit(funit, cme2, rhe2, cmo2, rhoo2, &
                            &fits(f)%rho_e, fits(f)%rho_o, fits(f)%gam_sum, fits(f)%nll_tot, &
                            &fits(f)%nval, fits(f)%ncomp, fits(f)%kpk_e, fits(f)%kpk_o, fits(f)%rpk_e, fits(f)%rpk_o, fits(f)%dm_sr, fits(f)%dm_sm, &
                            &fits(f)%dm_smm, fits(f)%dm_sai, fits(f)%dm_z, fits(f)%dm_nz)
                    endif
                else
                    if( f == 1 )then
                        call fold_probe_part_v5_fit(funit, cme1, rhe1, cmo1, rhoo1, &
                            &fits(f)%rho_e, fits(f)%rho_o, fits(f)%gam_sum, fits(f)%nll_tot, &
                            &fits(f)%nval, fits(f)%ncomp, fits(f)%kpk_e, fits(f)%kpk_o, fits(f)%rpk_e, fits(f)%rpk_o)
                    else
                        call fold_probe_part_v5_fit(funit, cme2, rhe2, cmo2, rhoo2, &
                            &fits(f)%rho_e, fits(f)%rho_o, fits(f)%gam_sum, fits(f)%nll_tot, &
                            &fits(f)%nval, fits(f)%ncomp, fits(f)%kpk_e, fits(f)%kpk_o, fits(f)%rpk_e, fits(f)%rpk_o)
                    endif
                endif
            end do
            call close_probe_part_v5_read(funit, fname)
            call fname%kill
        end do
        do q = 1, fits(1)%ncomp
            fits(1)%Yeven(q)%cmat_exp = cme1(:,:,:,q); fits(1)%Yeven(q)%rho_exp = rhe1(:,:,:,q)
            fits(1)%Yodd(q)%cmat_exp  = cmo1(:,:,:,q); fits(1)%Yodd(q)%rho_exp  = rhoo1(:,:,:,q)
        end do
        do q = 1, fits(2)%ncomp
            fits(2)%Yeven(q)%cmat_exp = cme2(:,:,:,q); fits(2)%Yeven(q)%rho_exp = rhe2(:,:,:,q)
            fits(2)%Yodd(q)%cmat_exp  = cmo2(:,:,:,q); fits(2)%Yodd(q)%rho_exp  = rhoo2(:,:,:,q)
        end do
        deallocate(cme1, cmo1, rhe1, rhoo1, cme2, cmo2, rhe2, rhoo2)
        write(logfhandle,'(A,I0,A,I0,A,I0,A,F8.1)') '>>> FLEX_PCA reduced v5 paired parts=', &
            &rounds%nparts(), '  valid A=', fits(1)%nval, '  valid B=', fits(2)%nval, &
            &'  seconds=', toc(t_red)
        call flush(logfhandle)
    end subroutine paired_reduce_parts_v5

    !> One probe iteration's accumulators from this part's particle range. Payload, in order: per
    !! component the even and odd Y_q reconstructor (cmat then rho), the coupled normal-matrix arrays
    !! rho_e / rho_o (one entry per (q,r) pair), the EM Gamma numerator and the valid-particle count.
    !! Gamma is shipped as a sum (the master divides by the reduced nval). The MCFA accumulators are
    !! additive sufficient statistics; header(9) carries kmix, 0 when absent.
    module subroutine write_probe_part( fname, cmat_e, rho_ex, cmat_o, rho_ox, rho_e, rho_o, &
        &gam_sum, nll_sum, nval, ncomp, kpk_e, kpk_o, rpk_e, rpk_o, mix_sr, mix_sm, mix_smm, mix_sainv, z_sub )
        class(string), intent(in) :: fname
        complex,       intent(in) :: cmat_e(:,:,:,:), cmat_o(:,:,:,:)
        real,          intent(in) :: rho_ex(:,:,:,:), rho_ox(:,:,:,:)
        real,          intent(in) :: rho_e(:,:,:,:),  rho_o(:,:,:,:)
        real, allocatable, intent(in) :: kpk_e(:,:), kpk_o(:,:)
        complex, allocatable, intent(in) :: rpk_e(:,:), rpk_o(:,:)
        real(dp),      intent(in) :: gam_sum(:)
        !> sum over this part's particles of log det A_i - h_i' A_i^-1 h_i, the particle-dependent
        !! part of -2 log p(y_i)
        real(dp),      intent(in) :: nll_sum
        integer,       intent(in) :: nval, ncomp
        real(dp), optional, intent(in) :: mix_sr(:), mix_sm(:,:), mix_smm(:,:,:), mix_sainv(:,:)
        !> bounded latent subsample: mcfa_init seeds k-means from z, so each part ships a slice
        !! and the master pools them
        real(dp), optional, intent(in) :: z_sub(:,:)
        type(string) :: tmp_fname
        integer :: funit, io_stat, header(9), q, kmix_w
        kmix_w = 0
        if( present(mix_sr) ) kmix_w = size(mix_sr)
        header = [FLEX_PCA_PART_MAGIC, PROBE_PART_VERSION, ncomp, &
            &size(cmat_e,1), size(cmat_e,2), size(cmat_e,3), size(rho_e,1), nval, kmix_w]
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_probe_part; open', io_stat)
        write(funit, iostat=io_stat) header
        call fileiochk('write_probe_part; header', io_stat)
        do q = 1, ncomp
            write(funit, iostat=io_stat) cmat_e(:,:,:,q), rho_ex(:,:,:,q), cmat_o(:,:,:,q), rho_ox(:,:,:,q)
            call fileiochk('write_probe_part; basis payload', io_stat)
        end do
        write(funit, iostat=io_stat) rho_e, rho_o, gam_sum, nll_sum
        call fileiochk('write_probe_part; coupled payload', io_stat)
        if( kmix_w > 0 )then
            write(funit, iostat=io_stat) mix_sr, mix_sm, mix_smm, mix_sainv
            call fileiochk('write_probe_part; mixture payload', io_stat)
            if( present(z_sub) )then
                write(funit, iostat=io_stat) size(z_sub,1)
                write(funit, iostat=io_stat) z_sub
            else
                write(funit, iostat=io_stat) 0
            endif
            call fileiochk('write_probe_part; z subsample', io_stat)
        endif
        call write_probe_part_kernels(funit, kpk_e, kpk_o, rpk_e, rpk_o, 'write_probe_part')
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call tmp_fname%kill
    end subroutine write_probe_part

    !> Streaming sum of every part into the master's accumulators; one part resident at a time.
    module subroutine reduce_probe_parts( params, nparts, cmat_e, rho_ex, cmat_o, rho_ox, rho_e, rho_o, &
        &gam_sum, nll_sum, nval, ncomp, kpk_e, kpk_o, rpk_e, rpk_o, mix_sr, mix_sm, mix_smm, mix_sainv, z_pool, nz_pool )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: nparts, ncomp
        complex,           intent(inout) :: cmat_e(:,:,:,:), cmat_o(:,:,:,:)
        real,              intent(inout) :: rho_ex(:,:,:,:), rho_ox(:,:,:,:)
        real,              intent(inout) :: rho_e(:,:,:,:),  rho_o(:,:,:,:)
        real, allocatable, intent(inout) :: kpk_e(:,:), kpk_o(:,:)
        complex, allocatable, intent(inout) :: rpk_e(:,:), rpk_o(:,:)
        real(dp),          intent(inout) :: gam_sum(:)
        real(dp),          intent(inout) :: nll_sum
        integer,           intent(inout) :: nval
        real(dp), optional, intent(inout) :: mix_sr(:), mix_sm(:,:), mix_smm(:,:,:), mix_sainv(:,:)
        real(dp), optional, intent(inout) :: z_pool(:,:)
        integer,  optional, intent(inout) :: nz_pool
        real(dp), allocatable :: zbuf(:,:)
        integer :: nz_part
        real(dp), allocatable :: msr(:), msm(:,:), msmm(:,:,:), msai(:,:)
        complex, allocatable :: cbuf(:,:,:)
        real,    allocatable :: rbuf(:,:,:), rcbuf(:,:,:,:)
        real(dp), allocatable :: gbuf(:)
        real(dp) :: nllbuf
        type(string) :: fname
        integer :: ipart, funit, io_stat, header(9), q
        integer(timer_int_kind) :: t_red
        t_red = tic()
        allocate(cbuf(size(cmat_e,1),size(cmat_e,2),size(cmat_e,3)))
        allocate(rbuf(size(rho_ex,1),size(rho_ex,2),size(rho_ex,3)))
        allocate(rcbuf(size(rho_e,1),size(rho_e,2),size(rho_e,3),size(rho_e,4)))
        allocate(gbuf(size(gam_sum)))
        if( present(mix_sr) ) allocate(msr(size(mix_sr)), msm(size(mix_sm,1),size(mix_sm,2)), &
            &msmm(size(mix_smm,1),size(mix_smm,2),size(mix_smm,3)), &
            &msai(size(mix_sainv,1),size(mix_sainv,2)))
        do ipart = 1, nparts
            fname = flex_pca_part_fname('probe', ipart, params%numlen)
            if( .not. file_exists(fname) ) THROW_HARD('missing probe part: '//fname%to_char())
            call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
            call fileiochk('reduce_probe_parts; open '//fname%to_char(), io_stat)
            read(funit, iostat=io_stat) header
            call fileiochk('reduce_probe_parts; header', io_stat)
            if( header(1) /= FLEX_PCA_PART_MAGIC ) THROW_HARD('bad probe part magic')
            if( header(2) /= PROBE_PART_VERSION  ) THROW_HARD('bad probe part version')
            if( header(3) /= ncomp               ) THROW_HARD('probe part ncomp mismatch')
            if( header(7) /= size(rho_e,1)       ) THROW_HARD('probe part npairs mismatch')
            do q = 1, ncomp
                read(funit, iostat=io_stat) cbuf
                call fileiochk('reduce_probe_parts; cmat_e', io_stat)
                !$omp parallel workshare
                cmat_e(:,:,:,q) = cmat_e(:,:,:,q) + cbuf
                !$omp end parallel workshare
                read(funit, iostat=io_stat) rbuf
                call fileiochk('reduce_probe_parts; rho_ex', io_stat)
                !$omp parallel workshare
                rho_ex(:,:,:,q) = rho_ex(:,:,:,q) + rbuf
                !$omp end parallel workshare
                read(funit, iostat=io_stat) cbuf
                call fileiochk('reduce_probe_parts; cmat_o', io_stat)
                !$omp parallel workshare
                cmat_o(:,:,:,q) = cmat_o(:,:,:,q) + cbuf
                !$omp end parallel workshare
                read(funit, iostat=io_stat) rbuf
                call fileiochk('reduce_probe_parts; rho_ox', io_stat)
                !$omp parallel workshare
                rho_ox(:,:,:,q) = rho_ox(:,:,:,q) + rbuf
                !$omp end parallel workshare
            end do
            read(funit, iostat=io_stat) rcbuf
            call fileiochk('reduce_probe_parts; rho_e', io_stat)
            !$omp parallel workshare
            rho_e = rho_e + rcbuf
            !$omp end parallel workshare
            read(funit, iostat=io_stat) rcbuf
            call fileiochk('reduce_probe_parts; rho_o', io_stat)
            !$omp parallel workshare
            rho_o = rho_o + rcbuf
            !$omp end parallel workshare
            read(funit, iostat=io_stat) gbuf
            call fileiochk('reduce_probe_parts; gamma', io_stat)
            gam_sum = gam_sum + gbuf
            read(funit, iostat=io_stat) nllbuf
            call fileiochk('reduce_probe_parts; loglik', io_stat)
            nll_sum = nll_sum + nllbuf
            nval    = nval + header(8)
            if( present(mix_sr) )then
                if( header(9) /= size(mix_sr) ) THROW_HARD('probe part kmix mismatch')
                read(funit, iostat=io_stat) msr, msm, msmm, msai
                call fileiochk('reduce_probe_parts; mixture', io_stat)
                mix_sr    = mix_sr    + msr
                mix_sm    = mix_sm    + msm
                mix_smm   = mix_smm   + msmm
                mix_sainv = mix_sainv + msai
                read(funit, iostat=io_stat) nz_part
            call fileiochk('fold_probe_part_v5_fit; latent subsample count', io_stat)
            if( nz_part < 0 .or. nz_part > 1000000 ) THROW_HARD('v5 probe part reader: corrupt latent subsample count (stream misaligned)')
                call fileiochk('reduce_probe_parts; z subsample count', io_stat)
                if( nz_part > 0 )then
                    allocate(zbuf(nz_part, size(mix_sm,1)))
                    read(funit, iostat=io_stat) zbuf
                    call fileiochk('reduce_probe_parts; z subsample', io_stat)
                    if( present(z_pool) .and. present(nz_pool) )then
                        do q = 1, nz_part
                            if( nz_pool >= size(z_pool,1) ) exit
                            nz_pool = nz_pool + 1
                            z_pool(nz_pool,:) = zbuf(q,:)
                        end do
                    endif
                    deallocate(zbuf)
                endif
            endif
            call fold_probe_part_kernels(funit, kpk_e, kpk_o, rpk_e, rpk_o, 'reduce_probe_parts')
            call fclose(funit)
            call del_file(fname)
            call fname%kill
        end do
        deallocate(cbuf, rbuf, rcbuf, gbuf)
        write(logfhandle,'(A,I0,A,I0,A,F8.1)') '>>> FLEX_PCA reduced probe parts=',nparts, &
            &'  valid particles=',nval,'  seconds=',toc(t_red)
        call flush(logfhandle)
    end subroutine reduce_probe_parts

    !> Open one v5 part for writing and emit the file header. The caller then writes one
    !! block per fit with write_probe_part_v5_fit, IN FIT ORDER, and closes with
    !! close_probe_part_v5_write (which renames the .tmp so the master only ever sees
    !! complete files -- same contract as write_probe_part).
    module subroutine open_probe_part_v5_write( fname, nfits, funit, tmp_fname )
        class(string), intent(in)  :: fname
        integer,       intent(in)  :: nfits
        integer,       intent(out) :: funit
        type(string),  intent(out) :: tmp_fname
        integer :: io_stat
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('open_probe_part_v5_write; open', io_stat)
        write(funit, iostat=io_stat) FLEX_PCA_PART_MAGIC, PROBE_PART_VERSION5, nfits
        call fileiochk('open_probe_part_v5_write; header', io_stat)
    end subroutine open_probe_part_v5_write

    !> One fit's block: sub-header + the write_probe_part payload.
    module subroutine write_probe_part_v5_fit( funit, cmat_e, rho_ex, cmat_o, rho_ox, rho_e, rho_o, &
        &gam_sum, nll_sum, nval, ncomp, kpk_e, kpk_o, rpk_e, rpk_o, mix_sr, mix_sm, mix_smm, mix_sainv, z_sub )
        integer,       intent(in) :: funit
        complex,       intent(in) :: cmat_e(:,:,:,:), cmat_o(:,:,:,:)
        real,          intent(in) :: rho_ex(:,:,:,:), rho_ox(:,:,:,:)
        real,          intent(in) :: rho_e(:,:,:,:),  rho_o(:,:,:,:)
        real, allocatable, intent(in) :: kpk_e(:,:), kpk_o(:,:)
        complex, allocatable, intent(in) :: rpk_e(:,:), rpk_o(:,:)
        real(dp),      intent(in) :: gam_sum(:)
        real(dp),      intent(in) :: nll_sum
        integer,       intent(in) :: nval, ncomp
        real(dp), optional, intent(in) :: mix_sr(:), mix_sm(:,:), mix_smm(:,:,:), mix_sainv(:,:)
        real(dp), optional, intent(in) :: z_sub(:,:)
        integer :: io_stat, subhdr(7), q, kmix_w, t, pnnz
        logical, allocatable :: pm(:,:,:)
        integer, allocatable :: pii(:), pjj(:), pkk(:)
        complex, allocatable :: gce(:), gco(:)
        real,    allocatable :: gre(:), gro(:), gpe(:,:), gpo(:,:)
        kmix_w = 0
        if( present(mix_sr) ) kmix_w = size(mix_sr)
        subhdr = [ncomp, size(cmat_e,1), size(cmat_e,2), size(cmat_e,3), size(rho_e,1), nval, kmix_w]
        write(funit, iostat=io_stat) subhdr
        call fileiochk('write_probe_part_v5_fit; sub-header', io_stat)
        ! the coupled rho rows must live on the SAME lattice as the basis accumulators, otherwise
        ! one band box cannot describe both and the boxed payload would be mis-scattered
        if( size(rho_e,2) /= size(cmat_e,1) .or. size(rho_e,3) /= size(cmat_e,2) .or. &
            &size(rho_e,4) /= size(cmat_e,3) ) &
            &THROW_HARD('write_probe_part_v5_fit: rho_e lattice differs from the basis lattice')
        ! one nonzero bounding box for the whole crop lattice (see "band boxing" above)
        ! index list of every populated lattice point, unioned over the six shipped arrays
        allocate(pm(size(cmat_e,1),size(cmat_e,2),size(cmat_e,3)), source=.false.)
        call pk_mask_c1(cmat_e, pm)
        call pk_mask_r1(rho_ex, pm)
        call pk_mask_c1(cmat_o, pm)
        call pk_mask_r1(rho_ox, pm)
        call pk_mask_r2(rho_e,  pm)
        call pk_mask_r2(rho_o,  pm)
        call pk_mask_to_idx(pm, pii, pjj, pkk, pnnz)
        deallocate(pm)
        write(funit, iostat=io_stat) pnnz
        call fileiochk('write_probe_part_v5_fit; index count', io_stat)
        write(funit, iostat=io_stat) pii, pjj, pkk
        call fileiochk('write_probe_part_v5_fit; index list', io_stat)
        allocate(gce(pnnz), gco(pnnz), gre(pnnz), gro(pnnz))
        do q = 1, ncomp
            do t = 1, pnnz
                gce(t) = cmat_e(pii(t),pjj(t),pkk(t),q)
                gre(t) = rho_ex(pii(t),pjj(t),pkk(t),q)
                gco(t) = cmat_o(pii(t),pjj(t),pkk(t),q)
                gro(t) = rho_ox(pii(t),pjj(t),pkk(t),q)
            end do
            write(funit, iostat=io_stat) gce, gre, gco, gro
            call fileiochk('write_probe_part_v5_fit; basis payload', io_stat)
        end do
        deallocate(gce, gco, gre, gro)
        allocate(gpe(size(rho_e,1),pnnz), gpo(size(rho_o,1),pnnz))
        do t = 1, pnnz
            gpe(:,t) = rho_e(:,pii(t),pjj(t),pkk(t))
            gpo(:,t) = rho_o(:,pii(t),pjj(t),pkk(t))
        end do
        write(funit, iostat=io_stat) gpe, gpo, gam_sum, nll_sum
        call fileiochk('write_probe_part_v5_fit; coupled payload', io_stat)
        deallocate(gpe, gpo, pii, pjj, pkk)
        if( kmix_w > 0 )then
            write(funit, iostat=io_stat) mix_sr, mix_sm, mix_smm, mix_sainv
            call fileiochk('write_probe_part_v5_fit; mixture payload', io_stat)
            if( present(z_sub) )then
                write(funit, iostat=io_stat) size(z_sub,1)
                write(funit, iostat=io_stat) z_sub
            else
                write(funit, iostat=io_stat) 0
            endif
            call fileiochk('write_probe_part_v5_fit; z subsample', io_stat)
        endif
        call write_probe_part_kernels(funit, kpk_e, kpk_o, rpk_e, rpk_o, 'write_probe_part_v5_fit')
    end subroutine write_probe_part_v5_fit

    ! ---- band boxing -------------------------------------------------------------------------
    ! Every probe-part payload is populated only inside the covariance band: the data planes are
    ! capped at projected_model_kfromto before insertion, so the expanded lattices carry EXACT
    ! zeros outside a small centred region (measured 6.5 % of the crop lattice and 6.9 % of the
    ! doubled PCG lattice at box_crop=128 / lp=8, i.e. ~93 % of every part file was zeros).
    ! We therefore ship only the NONZERO BOUNDING BOX of each lattice. The box is derived from the
    ! data, never from an assumed band radius, so nothing nonzero can be dropped: outside the box
    ! the part contributes exact zeros and the fold is bit-identical to shipping the full array.
    ! A box that is merely a superset of the nonzero region is equally correct (it just ships more
    ! zeros), which is why widening is always safe and narrowing never happens.
    ! box = [i0,i1,j0,j1,k0,k1] over the three LATTICE dimensions. Two memory layouts occur:
    !   layout 1: lattice is dims 1-3, slice index is dim 4  -- cmat_e/rho_ex(nx,ny,nz,ncomp)
    !   layout 2: lattice is dims 2-4, slice index is dim 1  -- rho_e/kpk_e/rpk_e(nslice,nx,ny,nz)
    ! The scan runs on the WRITE side, i.e. on the parallel part workers, never on the master, so
    ! it is off the critical path that the reduce sits on.

    subroutine box_reset( box )
        integer, intent(out) :: box(6)
        box = [huge(1), -huge(1), huge(1), -huge(1), huge(1), -huge(1)]
    end subroutine box_reset

    !> an empty box (nothing nonzero anywhere) collapses onto a single voxel: folding one zero
    !! voxel is still bit-identical and it keeps every extent positive for the record I/O.
    subroutine box_finalize( box, n1, n2, n3 )
        integer, intent(inout) :: box(6)
        integer, intent(in)    :: n1, n2, n3
        if( box(1) > box(2) .or. box(3) > box(4) .or. box(5) > box(6) )then
            box = [1,1,1,1,1,1]
            return
        endif
        box(1) = max(1,box(1)); box(2) = min(n1,box(2))
        box(3) = max(1,box(3)); box(4) = min(n2,box(4))
        box(5) = max(1,box(5)); box(6) = min(n3,box(6))
    end subroutine box_finalize

    !> the box is derived from the union of exactly the arrays that get shipped, so "every nonzero
    !! lies inside the box" holds by construction -- until someone adds a seventh array to the
    !! writer and forgets the union. These verify that invariant directly on the shipped arrays and
    !! are layout-, addressing- and wrap-agnostic. Write side only, so off the reduce's critical path.
    subroutine box_verify_c1( a, box, who )
        complex,          intent(in) :: a(:,:,:,:)
        integer,          intent(in) :: box(6)
        character(len=*), intent(in) :: who
        if( count(a /= (0.,0.)) /= &
            &count(a(box(1):box(2),box(3):box(4),box(5):box(6),:) /= (0.,0.)) ) &
            &THROW_HARD(who//': nonzero content outside the band box (box/union out of sync)')
    end subroutine box_verify_c1

    subroutine box_verify_r1( a, box, who )
        real,             intent(in) :: a(:,:,:,:)
        integer,          intent(in) :: box(6)
        character(len=*), intent(in) :: who
        if( count(a /= 0.) /= &
            &count(a(box(1):box(2),box(3):box(4),box(5):box(6),:) /= 0.) ) &
            &THROW_HARD(who//': nonzero content outside the band box (box/union out of sync)')
    end subroutine box_verify_r1

    subroutine box_verify_r2( a, box, who )
        real,             intent(in) :: a(:,:,:,:)
        integer,          intent(in) :: box(6)
        character(len=*), intent(in) :: who
        if( count(a /= 0.) /= &
            &count(a(:,box(1):box(2),box(3):box(4),box(5):box(6)) /= 0.) ) &
            &THROW_HARD(who//': nonzero content outside the band box (box/union out of sync)')
    end subroutine box_verify_r2

    subroutine box_verify_c2( a, box, who )
        complex,          intent(in) :: a(:,:,:,:)
        integer,          intent(in) :: box(6)
        character(len=*), intent(in) :: who
        if( count(a /= (0.,0.)) /= &
            &count(a(:,box(1):box(2),box(3):box(4),box(5):box(6)) /= (0.,0.)) ) &
            &THROW_HARD(who//': nonzero content outside the band box (box/union out of sync)')
    end subroutine box_verify_c2

    ! ---- index-list packing --------------------------------------------------------------------
    ! The bounding box loses most of its value on the DOUBLED PCG lattice: that payload is in
    ! physical FFT addressing, so real mass sits on opposite faces and an axis-aligned box covers
    ! ~50 % of the lattice even though only ~6.9 % of voxels are nonzero. An explicit index list of
    ! the populated lattice points is addressing-agnostic -- it does not care whether the support is
    ! a centred ball, eight wrapped octants, or anything else -- and recovers the full ~15x.
    ! The list is built from the UNION of the nonzeros of exactly the arrays being shipped, so
    ! nothing nonzero is ever dropped and the fold stays bit-identical (the omitted voxels add 0.0).
    ! One list per lattice per fit, shared by every slice, so its 4 bytes/voxel is amortised over
    ! ncomp (or npairs) x 2 half-sets of payload.

    subroutine pk_mask_r2( a, m )
        real,    intent(in)    :: a(:,:,:,:)
        logical, intent(inout) :: m(:,:,:)
        integer :: i, j, k
        do k = 1, size(a,4)
            do j = 1, size(a,3)
                do i = 1, size(a,2)
                    if( any(a(:,i,j,k) /= 0.) ) m(i,j,k) = .true.
                end do
            end do
        end do
    end subroutine pk_mask_r2

    subroutine pk_mask_c1( a, m )
        complex, intent(in)    :: a(:,:,:,:)
        logical, intent(inout) :: m(:,:,:)
        integer :: i, j, k
        !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static)
        do k = 1, size(a,3)
            do j = 1, size(a,2)
                do i = 1, size(a,1)
                    if( any(a(i,j,k,:) /= (0.,0.)) ) m(i,j,k) = .true.
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine pk_mask_c1

    subroutine pk_mask_r1( a, m )
        real,    intent(in)    :: a(:,:,:,:)
        logical, intent(inout) :: m(:,:,:)
        integer :: i, j, k
        !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static)
        do k = 1, size(a,3)
            do j = 1, size(a,2)
                do i = 1, size(a,1)
                    if( any(a(i,j,k,:) /= 0.) ) m(i,j,k) = .true.
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine pk_mask_r1



    !> flatten the union mask into the three coordinate lists the payload is gathered on.
    !! An empty mask keeps ONE voxel so every extent stays positive; folding one zero is a no-op.
    subroutine pk_mask_to_idx( m, ii, jj, kk, nnz )
        logical,              intent(in)  :: m(:,:,:)
        integer, allocatable, intent(out) :: ii(:), jj(:), kk(:)
        integer,              intent(out) :: nnz
        integer :: i, j, k, t
        nnz = count(m)
        if( nnz == 0 )then
            allocate(ii(1), jj(1), kk(1)); ii = 1; jj = 1; kk = 1; nnz = 1
            return
        endif
        allocate(ii(nnz), jj(nnz), kk(nnz))
        t = 0
        do k = 1, size(m,3)
            do j = 1, size(m,2)
                do i = 1, size(m,1)
                    if( m(i,j,k) )then
                        t = t + 1; ii(t) = i; jj(t) = j; kk(t) = k
                    endif
                end do
            end do
        end do
    end subroutine pk_mask_to_idx

    subroutine box_union_c1( a, box )
        complex, intent(in)    :: a(:,:,:,:)
        integer, intent(inout) :: box(6)
        integer :: i, j, k, i0, i1, j0, j1, k0, k1
        i0 = box(1); i1 = box(2); j0 = box(3); j1 = box(4); k0 = box(5); k1 = box(6)
        !$omp parallel do collapse(2) default(shared) private(i,j,k) schedule(static) &
        !$omp reduction(min:i0,j0,k0) reduction(max:i1,j1,k1)
        do k = 1, size(a,3)
            do j = 1, size(a,2)
                do i = 1, size(a,1)
                    if( any(a(i,j,k,:) /= (0.,0.)) )then
                        i0 = min(i0,i); i1 = max(i1,i)
                        j0 = min(j0,j); j1 = max(j1,j)
                        k0 = min(k0,k); k1 = max(k1,k)
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        box = [i0,i1,j0,j1,k0,k1]
    end subroutine box_union_c1

    subroutine box_union_r1( a, box )
        real,    intent(in)    :: a(:,:,:,:)
        integer, intent(inout) :: box(6)
        integer :: i, j, k, i0, i1, j0, j1, k0, k1
        i0 = box(1); i1 = box(2); j0 = box(3); j1 = box(4); k0 = box(5); k1 = box(6)
        !$omp parallel do collapse(2) default(shared) private(i,j,k) schedule(static) &
        !$omp reduction(min:i0,j0,k0) reduction(max:i1,j1,k1)
        do k = 1, size(a,3)
            do j = 1, size(a,2)
                do i = 1, size(a,1)
                    if( any(a(i,j,k,:) /= 0.) )then
                        i0 = min(i0,i); i1 = max(i1,i)
                        j0 = min(j0,j); j1 = max(j1,j)
                        k0 = min(k0,k); k1 = max(k1,k)
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        box = [i0,i1,j0,j1,k0,k1]
    end subroutine box_union_r1

    subroutine box_union_r2( a, box )
        real,    intent(in)    :: a(:,:,:,:)
        integer, intent(inout) :: box(6)
        integer :: i, j, k, i0, i1, j0, j1, k0, k1
        i0 = box(1); i1 = box(2); j0 = box(3); j1 = box(4); k0 = box(5); k1 = box(6)
        !$omp parallel do collapse(2) default(shared) private(i,j,k) schedule(static) &
        !$omp reduction(min:i0,j0,k0) reduction(max:i1,j1,k1)
        do k = 1, size(a,4)
            do j = 1, size(a,3)
                do i = 1, size(a,2)
                    if( any(a(:,i,j,k) /= 0.) )then
                        i0 = min(i0,i); i1 = max(i1,i)
                        j0 = min(j0,j); j1 = max(j1,j)
                        k0 = min(k0,k); k1 = max(k1,k)
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        box = [i0,i1,j0,j1,k0,k1]
    end subroutine box_union_r2

    subroutine box_union_c2( a, box )
        complex, intent(in)    :: a(:,:,:,:)
        integer, intent(inout) :: box(6)
        integer :: i, j, k, i0, i1, j0, j1, k0, k1
        i0 = box(1); i1 = box(2); j0 = box(3); j1 = box(4); k0 = box(5); k1 = box(6)
        !$omp parallel do collapse(2) default(shared) private(i,j,k) schedule(static) &
        !$omp reduction(min:i0,j0,k0) reduction(max:i1,j1,k1)
        do k = 1, size(a,4)
            do j = 1, size(a,3)
                do i = 1, size(a,2)
                    if( any(a(:,i,j,k) /= (0.,0.)) )then
                        i0 = min(i0,i); i1 = max(i1,i)
                        j0 = min(j0,j); j1 = max(j1,j)
                        k0 = min(k0,k); k1 = max(k1,k)
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        box = [i0,i1,j0,j1,k0,k1]
    end subroutine box_union_c2

    !> trailing PCG kernel block of one fit: the pair count (0 when the fit does not solve by PCG), then
    !! the packed even and odd kernel sums
    !> PCG kernel and rhs payload of a part: the packed arrays on the band list, slot for slot. Every
    !! process derives the same list from (box, band, rank), so the reader adds by slot; the list
    !! length travels in the shape and is checked. (Format 12: no per-part index list.)
    subroutine write_probe_part_kernels( funit, kpk_e, kpk_o, rpk_e, rpk_o, who )
        integer,           intent(in) :: funit
        real, allocatable, intent(in) :: kpk_e(:,:), kpk_o(:,:)
        complex, allocatable, intent(in) :: rpk_e(:,:), rpk_o(:,:)
        character(len=*),  intent(in) :: who
        integer :: io_stat, nkp, nrp
        nkp = 0
        if( allocated(kpk_e) .and. allocated(kpk_o) ) nkp = size(kpk_e,1)
        write(funit, iostat=io_stat) nkp
        call fileiochk(who//'; PCG kernel count', io_stat)
        if( nkp > 0 )then
            write(funit, iostat=io_stat) shape(kpk_e)
            write(funit, iostat=io_stat) kpk_e, kpk_o
            call fileiochk(who//'; PCG kernel payload', io_stat)
        endif
        nrp = 0
        if( allocated(rpk_e) .and. allocated(rpk_o) ) nrp = size(rpk_e,1)
        write(funit, iostat=io_stat) nrp
        call fileiochk(who//'; PCG rhs count', io_stat)
        if( nrp > 0 )then
            write(funit, iostat=io_stat) shape(rpk_e)
            write(funit, iostat=io_stat) rpk_e, rpk_o
            call fileiochk(who//'; PCG rhs payload', io_stat)
        endif
    end subroutine write_probe_part_kernels

    !> read one fit's trailing PCG kernel block and ADD it into the packed sums (a part written by a
    !! non-PCG fit carries a zero count; a PCG part folded into a non-PCG reduce is a mismatch)
    subroutine fold_probe_part_kernels( funit, kpk_e, kpk_o, rpk_e, rpk_o, who )
        integer,           intent(in)    :: funit
        real, allocatable, intent(inout) :: kpk_e(:,:), kpk_o(:,:)
        complex, allocatable, intent(inout) :: rpk_e(:,:), rpk_o(:,:)
        character(len=*),  intent(in)    :: who
        integer :: io_stat, nkp, nrp, kshape(2)
        real,    allocatable :: gke(:,:), gko(:,:)
        complex, allocatable :: gqe(:,:), gqo(:,:)
        read(funit, iostat=io_stat) nkp
        call fileiochk(who//'; PCG kernel count', io_stat)
        if( nkp == 0 )then
            if( allocated(kpk_e) ) THROW_HARD(who//': part carries no PCG kernels but the reduce expects them')
        else
            if( .not. (allocated(kpk_e) .and. allocated(kpk_o)) ) &
                &THROW_HARD(who//': part carries PCG kernels the reduce did not expect')
            read(funit, iostat=io_stat) kshape
            call fileiochk(who//'; PCG kernel shape', io_stat)
            if( any(kshape /= shape(kpk_e)) ) THROW_HARD(who//': PCG kernel band-list shape mismatch (part from another lattice or band)')
            allocate(gke(kshape(1),kshape(2)), gko(kshape(1),kshape(2)))
            read(funit, iostat=io_stat) gke, gko
            call fileiochk(who//'; PCG kernel payload', io_stat)
            !$omp parallel workshare default(shared)
            kpk_e = kpk_e + gke
            kpk_o = kpk_o + gko
            !$omp end parallel workshare
            deallocate(gke, gko)
        endif
        read(funit, iostat=io_stat) nrp
        call fileiochk(who//'; PCG rhs count', io_stat)
        if( nrp == 0 )then
            if( allocated(rpk_e) ) THROW_HARD(who//': part carries no PCG right-hand sides but the reduce expects them')
            return
        endif
        if( .not. (allocated(rpk_e) .and. allocated(rpk_o)) ) &
            &THROW_HARD(who//': part carries PCG right-hand sides the reduce did not expect')
        read(funit, iostat=io_stat) kshape
        call fileiochk(who//'; PCG rhs shape', io_stat)
        if( any(kshape /= shape(rpk_e)) ) THROW_HARD(who//': PCG rhs band-list shape mismatch (part from another lattice or band)')
        allocate(gqe(kshape(1),kshape(2)), gqo(kshape(1),kshape(2)))
        read(funit, iostat=io_stat) gqe, gqo
        call fileiochk(who//'; PCG rhs payload', io_stat)
        !$omp parallel workshare default(shared)
        rpk_e = rpk_e + gqe
        rpk_o = rpk_o + gqo
        !$omp end parallel workshare
        deallocate(gqe, gqo)
    end subroutine fold_probe_part_kernels

    module subroutine close_probe_part_v5_write( funit, tmp_fname, fname )
        integer,       intent(in)    :: funit
        type(string),  intent(inout) :: tmp_fname
        class(string), intent(in)    :: fname
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call tmp_fname%kill
    end subroutine close_probe_part_v5_write

    !> Open one v5 part for reduction and validate the file header. The caller folds one
    !! block per fit with fold_probe_part_v5_fit, IN FIT ORDER (the file is streamed).
    module subroutine open_probe_part_v5_read( fname, nfits, funit )
        class(string), intent(in)  :: fname
        integer,       intent(in)  :: nfits
        integer,       intent(out) :: funit
        integer :: io_stat, magic, ver, nfits_in
        if( .not. file_exists(fname) ) THROW_HARD('missing v5 probe part: '//fname%to_char())
        call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('open_probe_part_v5_read; open '//fname%to_char(), io_stat)
        read(funit, iostat=io_stat) magic, ver, nfits_in
        call fileiochk('open_probe_part_v5_read; header', io_stat)
        if( magic /= FLEX_PCA_PART_MAGIC  ) THROW_HARD('bad v5 probe part magic')
        if( ver   /= PROBE_PART_VERSION5  ) THROW_HARD('bad v5 probe part version')
        if( nfits_in /= nfits             ) THROW_HARD('v5 probe part nfits mismatch')
    end subroutine open_probe_part_v5_read

    !> Fold one fit's block from an open v5 part into that fit's accumulators; ncomp is checked
    !! against that fit's own expectation.
    module subroutine fold_probe_part_v5_fit( funit, cmat_e, rho_ex, cmat_o, rho_ox, rho_e, rho_o, &
        &gam_sum, nll_sum, nval, ncomp, kpk_e, kpk_o, rpk_e, rpk_o, mix_sr, mix_sm, mix_smm, mix_sainv, z_pool, nz_pool )
        integer,           intent(in)    :: funit
        complex,           intent(inout) :: cmat_e(:,:,:,:), cmat_o(:,:,:,:)
        real,              intent(inout) :: rho_ex(:,:,:,:), rho_ox(:,:,:,:)
        real,              intent(inout) :: rho_e(:,:,:,:),  rho_o(:,:,:,:)
        real, allocatable, intent(inout) :: kpk_e(:,:), kpk_o(:,:)
        complex, allocatable, intent(inout) :: rpk_e(:,:), rpk_o(:,:)
        real(dp),          intent(inout) :: gam_sum(:)
        real(dp),          intent(inout) :: nll_sum
        integer,           intent(inout) :: nval
        integer,           intent(in)    :: ncomp
        real(dp), optional, intent(inout) :: mix_sr(:), mix_sm(:,:), mix_smm(:,:,:), mix_sainv(:,:)
        real(dp), optional, intent(inout) :: z_pool(:,:)
        integer,  optional, intent(inout) :: nz_pool
        real(dp), allocatable :: gbuf(:), zbuf(:,:)
        real(dp), allocatable :: msr(:), msm(:,:), msmm(:,:,:), msai(:,:)
        real(dp) :: nllbuf
        integer  :: io_stat, subhdr(7), q, nz_part, t, pnnz
        integer, allocatable :: pii(:), pjj(:), pkk(:)
        complex, allocatable :: gce(:), gco(:)
        real,    allocatable :: gre(:), gro(:), gpe(:,:), gpo(:,:)
        read(funit, iostat=io_stat) subhdr
        call fileiochk('fold_probe_part_v5_fit; sub-header', io_stat)
        if( subhdr(1) /= ncomp          ) THROW_HARD('v5 probe part per-fit ncomp mismatch')
        if( subhdr(5) /= size(rho_e,1)  ) THROW_HARD('v5 probe part per-fit npairs mismatch')
        ! the lattice dims MUST be validated, not just ncomp/npairs: a part written on a different
        ! expanded lattice yields an in-range band box whose voxels land at the wrong addresses,
        ! i.e. silently misplaced mass instead of a loud failure
        if( subhdr(2) /= size(cmat_e,1) .or. subhdr(3) /= size(cmat_e,2) .or. &
            &subhdr(4) /= size(cmat_e,3) ) &
            &THROW_HARD('v5 probe part per-fit lattice mismatch (part written on another lattice)')
        ! index list of the populated lattice points (see "index-list packing"): the part carries
        ! only these voxels; every omitted one was an exact zero on the write side, so the fold is
        ! bit-identical to adding the full array
        read(funit, iostat=io_stat) pnnz
        call fileiochk('fold_probe_part_v5_fit; index count', io_stat)
        if( pnnz < 1 ) THROW_HARD('v5 probe part index count is not positive')
        allocate(pii(pnnz), pjj(pnnz), pkk(pnnz))
        read(funit, iostat=io_stat) pii, pjj, pkk
        call fileiochk('fold_probe_part_v5_fit; index list', io_stat)
        if( minval(pii) < 1 .or. maxval(pii) > size(cmat_e,1) .or. &
            &minval(pjj) < 1 .or. maxval(pjj) > size(cmat_e,2) .or. &
            &minval(pkk) < 1 .or. maxval(pkk) > size(cmat_e,3) ) &
            &THROW_HARD('v5 probe part index list outside the basis lattice')
        allocate(gce(pnnz), gco(pnnz), gre(pnnz), gro(pnnz))
        allocate(gbuf(size(gam_sum)))
        do q = 1, ncomp
            read(funit, iostat=io_stat) gce, gre, gco, gro
            call fileiochk('fold_probe_part_v5_fit; basis payload', io_stat)
            !$omp parallel do default(shared) private(t) schedule(static)
            do t = 1, pnnz
                cmat_e(pii(t),pjj(t),pkk(t),q) = cmat_e(pii(t),pjj(t),pkk(t),q) + gce(t)
                rho_ex(pii(t),pjj(t),pkk(t),q) = rho_ex(pii(t),pjj(t),pkk(t),q) + gre(t)
                cmat_o(pii(t),pjj(t),pkk(t),q) = cmat_o(pii(t),pjj(t),pkk(t),q) + gco(t)
                rho_ox(pii(t),pjj(t),pkk(t),q) = rho_ox(pii(t),pjj(t),pkk(t),q) + gro(t)
            end do
            !$omp end parallel do
        end do
        deallocate(gce, gco, gre, gro)
        allocate(gpe(size(rho_e,1),pnnz), gpo(size(rho_o,1),pnnz))
        read(funit, iostat=io_stat) gpe, gpo
        call fileiochk('fold_probe_part_v5_fit; coupled payload', io_stat)
        !$omp parallel do default(shared) private(t) schedule(static)
        do t = 1, pnnz
            rho_e(:,pii(t),pjj(t),pkk(t)) = rho_e(:,pii(t),pjj(t),pkk(t)) + gpe(:,t)
            rho_o(:,pii(t),pjj(t),pkk(t)) = rho_o(:,pii(t),pjj(t),pkk(t)) + gpo(:,t)
        end do
        !$omp end parallel do
        deallocate(gpe, gpo)
        read(funit, iostat=io_stat) gbuf
        call fileiochk('fold_probe_part_v5_fit; gamma', io_stat)
        gam_sum = gam_sum + gbuf
        read(funit, iostat=io_stat) nllbuf
        call fileiochk('fold_probe_part_v5_fit; loglik', io_stat)
        nll_sum = nll_sum + nllbuf
        nval    = nval + subhdr(6)
        if( present(mix_sr) )then
            if( subhdr(7) /= size(mix_sr) ) THROW_HARD('v5 probe part per-fit kmix mismatch')
            if( size(mix_sm,1) /= ncomp .or. size(mix_smm,1) /= ncomp .or. size(mix_sainv,1) /= ncomp ) &
                &THROW_HARD('v5 probe part reader: mixture accumulators not sized to ncomp (stale after a rank change)')
            allocate(msr(size(mix_sr)), msm(size(mix_sm,1),size(mix_sm,2)), &
                &msmm(size(mix_smm,1),size(mix_smm,2),size(mix_smm,3)), &
                &msai(size(mix_sainv,1),size(mix_sainv,2)))
            read(funit, iostat=io_stat) msr, msm, msmm, msai
            call fileiochk('fold_probe_part_v5_fit; mixture', io_stat)
            mix_sr    = mix_sr    + msr
            mix_sm    = mix_sm    + msm
            mix_smm   = mix_smm   + msmm
            mix_sainv = mix_sainv + msai
            deallocate(msr, msm, msmm, msai)
            read(funit, iostat=io_stat) nz_part
            call fileiochk('fold_probe_part_v5_fit; latent subsample count', io_stat)
            if( nz_part < 0 .or. nz_part > 1000000 ) THROW_HARD('v5 probe part reader: corrupt latent subsample count (stream misaligned)')
            call fileiochk('fold_probe_part_v5_fit; z subsample count', io_stat)
            if( nz_part > 0 )then
                allocate(zbuf(nz_part, size(mix_sm,1)))
                read(funit, iostat=io_stat) zbuf
                call fileiochk('fold_probe_part_v5_fit; z subsample', io_stat)
                if( present(z_pool) .and. present(nz_pool) )then
                    do q = 1, nz_part
                        if( nz_pool >= size(z_pool,1) ) exit
                        nz_pool = nz_pool + 1
                        z_pool(nz_pool,:) = zbuf(q,:)
                    end do
                endif
                deallocate(zbuf)
            endif
        else
            if( subhdr(7) /= 0 ) THROW_HARD('v5 probe part carries a mixture block this fit did not expect')
        endif
        call fold_probe_part_kernels(funit, kpk_e, kpk_o, rpk_e, rpk_o, 'fold_probe_part_v5_fit')
        deallocate(gbuf, pii, pjj, pkk)
    end subroutine fold_probe_part_v5_fit

    module subroutine close_probe_part_v5_read( funit, fname )
        integer,       intent(in) :: funit
        class(string), intent(in) :: fname
        call fclose(funit)
        call del_file(fname)
    end subroutine close_probe_part_v5_read

end submodule simple_flex_pca_em_estep
