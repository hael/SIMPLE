!@descr: covariance-column estimation, reduced covariance solve and latent embedding for flex_pca
module simple_flex_pca_columns
use simple_core_module_api
use simple_builder,          only: builder
use simple_image,            only: image
use simple_parameters,       only: parameters
use simple_reconstructor,    only: reconstructor
use simple_flex_pca_distr,   only: flex_pca_is_master, flex_pca_is_worker, flex_pca_nparts, &
    &flex_pca_run_stage, PCA_STAGE_PROBE, PCA_STAGE_EMBED
use simple_flex_pca_parts,   only: flex_pca_part_fname, write_snr_part, reduce_snr_parts, &
    &write_cols_part, reduce_cols_parts, write_solve_part, reduce_solve_parts, &
    &write_probe_part, reduce_probe_parts, &
    &write_embed_stats_part, reduce_embed_zhalf_parts, read_embed_stats_part, &
    &write_mean_scale, read_mean_scale, &
    &write_columns_hkl, read_columns_hkl
use simple_kbinterpol,       only: kbinterpol
use simple_gridding,         only: prep3D_inv_kbenvelope4mul
use simple_linalg,           only: jacobi, eigsrt, svd_solve
use simple_math,             only: ceil_div, floor_div
use simple_srch_sort_loc,    only: hpsort
use simple_matcher_3Drec,    only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io,  only: discrete_read_imgbatch, prepimgbatch
use simple_flex_reconstructor_latent_ops, only: project_fplane_mean, project_fplanes_mean_basis, &
    &latent_projection_weights, weighted_expanded_cmat, LATENT_WDIM, &
    &insert_planes_oversamp_coupled_batch_scaled
use simple_flex_projected_latent_model,   only: prep_imgs4projected_model, solve_coupled_basis_exp, &
    &projected_model_kfromto
use simple_math_ft, only: resample_sigma2
use simple_flex_gpu, only: flex_gpu_available, flex_gpu_coupled_begin_f, &
    &flex_gpu_coupled_batch_raw_f, flex_gpu_coupled_end_f, &
    &flex_gpu_coupled_bank_f, flex_gpu_coupled_batch_banked_f, flex_gpu_coupled_bank_free_f, &
    &flex_gpu_estep_vols_f, flex_gpu_estep_batch_f, flex_gpu_estep_resid_f, flex_gpu_estep_free_f, &
    &flex_gpu_coupled_batch_banked_res_f, &
    &flex_gpu_prep_begin_f, flex_gpu_prep_batch_f, flex_gpu_prep_free_f, flex_gpu_estep_batch_res_f, &
    &flex_gpu_prep_check_f, flex_gpu_prep_ready, &
    &flex_gpu_psample_begin_f, flex_gpu_psample_batch_f, flex_gpu_psample_free_f, &
    &flex_gpu_psample_batch_res_f, flex_gpu_psample_diag_begin_f, &
    &flex_gpu_psample_diag_batch_f, flex_gpu_psample_diag_fetch_f, &
    &flex_gpu_cols_begin_f, flex_gpu_cols_batch_f, flex_gpu_cols_end_f, &
    &flex_gpu_cols_lookup_f, flex_gpu_cols_fused_batch_f, flex_gpu_cols_test_upload_f, &
    &flex_gpu_snr_begin_f, flex_gpu_snr_batch_f, flex_gpu_snr_batch_res_f, flex_gpu_snr_end_f, &
    &flex_gpu_solve_begin_f, flex_gpu_solve_cache_f, flex_gpu_solve_chunk_f, flex_gpu_solve_end_f, &
    &flex_gpu_poles_begin_f, flex_gpu_poles_bank_f, flex_gpu_poles_batch_f, flex_gpu_poles_free_f
use simple_ftiter,           only: ftiter
use simple_flex_pca_polar,   only: polar_grid_t, polar_grid_build, polar_grid_kill, &
    &polar_project_recs, polar_sample_particle, polar_relative_inplane, polar_assign_directions, &
    &polar_sample_at_pose, polar_apply_shift, polar_dir_neighbours, polar_sample_particle_fused
implicit none
private
#include "simple_local_flags.inc"

public :: build_covariance_eigenbasis, embed_latents_with_contrast, estimate_covariance_mean
public :: probe_subspace_iteration, align_basis_to_reference, probe_external_basis
public :: bag_basis_pool, basis_recs_from_images
public :: cov_env_int_pub, save_probe_state, cov_perturb_project_poses

! Density observability floor, matching simple_image_arith::div_cmat_at_1 and the projected-latent coupled
! solve.
real(dp), parameter :: COV_DENSITY_FLOOR = 1.0d-6
! Resident-volume capacity of the fused device E-step. MUST match the `ncomp1 > 64` guard and
! the u_re/u_im extents in cuda/simple_flex_gpu_kernels.cu; exceeding it returns a hard error.
! Raised 24 -> 64 (2026-08-18): the per-thread arrays are dynamically indexed, hence in local
! memory either way -- 24 was a stack-frame guess, not a register wall. neigs=40 now rides the
! device path; parity/speed vs the CPU E-step is validated per-rank by the phase-12 A/B.
integer,  parameter :: COV_GPU_ESTEP_MAXVOLS = 64
! Relative ridge used ONLY for the covariance diagonal in the S.B SNR proxy, which runs before the S.C
! weights exist (Algorithm 1 precedes Algorithm 2).
real,     parameter :: COV_RIDGE_REL     = 5.0e-2
! Relative eigenvalue floor for retaining direct-column PCA components.
real(dp), parameter :: COV_EIG_REL_FLOOR = 1.0d-6
! Cap on the column-subspace dimension. The accumulation is a batched dsyrk on the Van Loan-Pitsianis
! rearrangement (see unrearrange_kron_selfsum), which needs ONE shared d^4 array regardless of thread count.
integer,  parameter :: COV_MAX_DTILDE    = 320
! Default column-subspace dimension, applied as a min against the memory budget so the rank follows
! the data rather than free RAM. Override with SIMPLE_COV_DTILDE.
integer,  parameter :: COV_DEFAULT_DTILDE = 128
! Total particles the probe / SNR / column-accumulation initialiser fit on, summed across processes.
! 0 = OFF: capping traded a recovered conformation for speed, which is not a trade worth taking.
! See doc/policies/flex_pca_policy.md. Enable per-run with SIMPLE_COV_PROBE_MAX / SIMPLE_COV_BASIS_MAX.
integer,  parameter :: COV_PROBE_MAX_PTCLS = 0
integer,  parameter :: COV_BASIS_MAX_PTCLS = 0
! Particle budget for em_calibrate_noise_prior. That pass estimates exactly TWO global
! scalars (the whitened-noise constant sig2 and the initial prior variance Gamma^0), whose
! precision improves only as 1/sqrt(N) -- 20k particles already give ~0.7% on sig2, far
! tighter than the EM needs. It previously shared SIMPLE_COV_BASIS_MAX, whose default of 0
! means UNCAPPED, so on a full dataset it read every particle on the MASTER ALONE (workers
! idle, no stage guard covers it): measured ~5 min of 1/10-capacity time at 105k/box64
! before the first distributed round. Capping removes the work rather than parallelising it.
integer,  parameter :: COV_CALIB_MAX_PTCLS = 20000
! How far above the spectrum's noise bulk a direction must stand to count as signal. Loose by design:
! keeping a noise direction costs one rank, dropping a real one costs a conformation.
real(dp), parameter :: COV_SIGNAL_FACTOR = 4.0d0
! Samples per free parameter for the rank bound d ~ sqrt(2N/R). REPORT ONLY.
real(dp), parameter :: COV_SAMPLES_PER_PARAM = 10.0d0
!> probe stops when successive bases agree to this mean principal-angle cosine: tight enough that the
!! remaining rotation cannot move a state target, and it fires at the measured knee rather than
!! running the tuned count out.
!> Amplitude of the project-level pose degradation, if any. Recorded so the pose block can report
!! how much of a KNOWN error it removed, which is the only honest way to ask that on a dataset whose
!! poses are ground truth to begin with.
real :: COV_PROJ_PERTURB_ROT = 0.0, COV_PROJ_PERTURB_SH = 0.0, COV_PROJ_PERTURB_DIR = 0.0
!> Viewing directions BEFORE the degradation, kept so the direction search can be scored against
!! truth. A direction is 2 DOF on a sphere, so unlike the in-plane angle it cannot be reconstructed
!! from the hash alone without also knowing the frame it was tilted in.
real, allocatable :: COV_TRUTH_NRM(:,:)
real(dp), parameter :: COV_PROBE_CONV    = 0.97d0
!> Convergence is declared when the REPRODUCIBLE dimension -- the sum of principal-angle cosines
!! between the even and odd half-bases, which the M-step already produces every iteration -- has
!! not improved by more than COV_EO_TOL for COV_EO_PATIENCE consecutive iterations. This is the
!! only criterion here with no dataset-specific constant in it: it asks the data how much of the
!! basis survives a change of particles, and stops when that stops growing.
real(dp), parameter :: COV_EO_TOL        = 0.02d0
integer,  parameter :: COV_EO_PATIENCE   = 3
! Memory budget for the shared A accumulator, in bytes.
real(dp), parameter :: COV_ATHR_BUDGET   = 8.0d9
! Accumulate the columns against the unscaled mean (a==1) rather than the per-particle ML contrast a_i.
! Subtracting a_i*T*mu also deletes the component of the conformational signal parallel to T*mu.
logical,  parameter :: COV_UNIT_CONTRAST  = .true.
! Runtime override of COV_UNIT_CONTRAST (SIMPLE_COV_CONTRAST=1): accumulate deviations against the
! per-particle fitted scale a_i = Re<Tmu,y>/Re<Tmu,Tmu> instead of a==1. JUSTIFICATION (RECOVAR,
! PNAS 2025, Appendix A.4): with per-image contrast a_i, the unscaled covariance "is corrupted by a
! rank-one component" -- a mean-shaped spurious principal component that embeds contrast in the
! latent "along with the heterogeneity, resulting in poor interpretability". That is this dataset's
! measured signature (all-positive PC1; latents tracking fitted contrast). Synthetic data has
! uniform contrast, which is why the defect never shows there. Set ONCE before any parallel region.
logical :: cov_fit_contrast_rt = .false.
! cross-halfset generalizing rank (consumed by auto-neigs)
integer :: cov_d_gen_rt = 0
! Grid-search the per-particle contrast in the embedding instead of using the closed-form estimate.
logical,  parameter :: COV_EMBED_CONTRAST_GRID = .false.
integer,  parameter :: COV_CG_MAXIT = 2000     ! CG iteration cap; convergence is reported, not assumed
real(dp), parameter :: COV_CG_TOL   = 1.d-10   ! relative residual target
integer,  parameter :: GRAM_DIAG_STRIDE = 200   ! subsample for the projected-Gram spectrum
integer,  parameter :: NCONTRAST_GRID = 50
real(dp), parameter :: A_GRID_HI = 2.0d0
! bracket for the fitted per-particle contrast; wide enough not to bind on real amplitude spread,
! tight enough that a particle the mask or the band has emptied cannot drag its latent to infinity
real(dp), parameter :: A_CONTRAST_LO = 0.2d0
real(dp), parameter :: A_CONTRAST_HI = 3.0d0
real(dp), parameter :: COV_PINV_RCOND = 1.0d-6

! Source of the covariance mean mu.
logical, parameter :: COV_MEAN_FROM_DATA = .false.

! Soft-mask each particle image to the projected molecular envelope before the column accumulation, so
! solvent (pure noise) does not enter the inner products.
logical, parameter :: COV_MASK_IMAGES = .false.
! Radial margin on that disc, as a multiple of the model radius. It must cover CTF delocalisation,
! lambda*defocus/d, which reaches ~70 A at 5.5 um defocus and 15 A resolution.
real, parameter :: COV_MASK_MARGIN = 1.4

! Subtract the analytic per-sample noise bias K_R(.,q_s)|T|^2 from the column numerator. Without it the
! bias survives into the half-set column FSC and the Wiener shrinkage deletes the low-frequency band.
logical, parameter :: COV_COLUMN_NOISE_DEBIAS = .true.

! Width of the RIGHT kernel -- the one that reads each image's value AT the column frequency
! (gather_column_values). Zero uses the shared 3-tap KB backprojection stencil for both, whose support is
! |d| <= 1.5 per axis against |d| < 2, so it gathers roughly half as many image samples into each column.
character(len=*), parameter :: COV_UTILDE_FBODY = 'flex_pca_utilde'
character(len=*), parameter :: COV_UTILDE_META  = 'flex_pca_utilde.txt'
!> master -> probe-worker handoff: the basis dimension, its prior variances and the whitened-noise
!! level. The basis volumes themselves are already on disk as flex_pca_pc*.mrc.
character(len=*), parameter :: COV_PROBE_META   = 'flex_pca_probe.txt'
real :: COV_RIGHT_KERNEL_W = 0.0
logical :: cov_rkw_read = .false.
! Half-width of the KB backprojection stencil in grid units, as cov_kb_weights derives it.
integer, parameter :: COV_KB_IWINSZ = ceiling(KBWINSZ - 0.5)

contains


    !> Pearson correlation of two double vectors.
    real(dp) function corr_dp( a, b, n ) result( r )
        integer,  intent(in) :: n
        real(dp), intent(in) :: a(n), b(n)
        real(dp) :: ma, mb, sa, sb, sab
        integer  :: i
        r  = 0.d0
        if( n < 3 ) return
        ma = sum(a)/real(n,dp); mb = sum(b)/real(n,dp)
        sa = 0.d0; sb = 0.d0; sab = 0.d0
        do i = 1, n
            sa  = sa  + (a(i)-ma)**2
            sb  = sb  + (b(i)-mb)**2
            sab = sab + (a(i)-ma)*(b(i)-mb)
        end do
        if( sa <= DTINY .or. sb <= DTINY ) return
        r = sab / sqrt(sa*sb)
    end function corr_dp

    !>  True only when an environment variable is explicitly set to zero (an opt-OUT switch).
    logical function cov_env_int_off( name ) result(off)
        character(len=*), intent(in) :: name
        character(len=32) :: envval
        integer :: stat, ln, ival
        off = .false.
        call get_environment_variable(name, envval, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        read(envval(:ln), *, iostat=stat) ival
        if( stat == 0 ) off = ival == 0
    end function cov_env_int_off

    ! True only when set AND non-zero. Not the complement of cov_env_int_off: an opt-in reads unset as OFF.
    logical function cov_env_int_on( name ) result(on)
        character(len=*), intent(in) :: name
        character(len=32) :: envval
        integer :: stat, ln, ival
        on = .false.
        call get_environment_variable(name, envval, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        read(envval(:ln), *, iostat=stat) ival
        if( stat == 0 ) on = ival /= 0
    end function cov_env_int_on

    ! Rank at which the Gram spectrum enters its noise bulk. Noise level = median of the lower half,
    ! so the leading signal directions cannot inflate it. Scale-free.
    pure integer function cov_signal_rank( eval, n ) result( d )
        integer,  intent(in) :: n
        real(dp), intent(in) :: eval(n)          !< DESCENDING eigenvalues
        real(dp) :: noise
        integer  :: lo, m
        d = 1
        if( n < 4 ) return
        lo    = n/2 + 1
        m     = n - lo + 1
        noise = eval(lo + m/2)
        if( noise <= DTINY )then
            d = n
            return
        endif
        d = 0
        do while( d < n )
            if( eval(d+1) <= COV_SIGNAL_FACTOR*noise ) exit
            d = d + 1
        end do
        d = max(1, min(n, d))
    end function cov_signal_rank

    ! Halfset-safe capped subsample, shared by the column-subspace initialiser and the probe EM.
    ! `eo` alternates strictly by particle index, so a plain stride of 2 selects one halfset entirely
    ! and the even/odd FSC that regularises every M-step is then computed against nothing; stride
    ! WITHIN each halfset instead. `maxtot` is a total across processes, so only a WORKER passes
    ! nparts -- the master holds every particle and dividing there inflates the stride by nparts.
    subroutine cov_stage_subsample( build, pinds, nptcls, nparts, maxtot, env_stride, env_max, &
        &label, spinds, nsel )
        type(builder),        intent(inout) :: build
        integer,              intent(in)    :: pinds(:), nptcls, nparts, maxtot
        character(len=*),     intent(in)    :: env_stride, env_max, label
        integer, allocatable, intent(out)   :: spinds(:)
        integer,              intent(out)   :: nsel
        integer :: stride, nmax_tot, nmax_part, ihalf, i, ii, nkept, n_half, ntgt, vstrat
        logical :: l_stride, l_ostrat
        integer, allocatable :: worder(:)
        real,    allocatable :: okey(:)
        real :: th_os, ph_os
        nmax_tot = maxtot
        call cov_env_int(env_max, nmax_tot)
        ! cap off (the default): hand back every particle, in project order
        if( nmax_tot < 1 )then
            allocate(spinds(nptcls), source=pinds(:nptcls))
            nsel = nptcls
            call hpsort(spinds)
            return
        endif
        nmax_part = max(1, nmax_tot / max(1, nparts))
        ! integer stride quantises the budget badly (100000 against a cap of 80000 wants 1.25, gets 2,
        ! keeps 50000), so default to exact-count Bresenham. env_stride restores it; =1 keeps all.
        stride   = 0
        call cov_env_int(env_stride, stride)
        l_stride = stride > 0
        stride   = max(1, stride)
        ! ---- SIMPLE_COV_OSTRAT=1: walk each halfset in ORIENTATION order instead of project
        ! order before the Bresenham pick. Project order is micrograph order, and orientation
        ! correlates with micrograph (ice, support, charging), so an acquisition-order stride
        ! can over/under-hit view clusters; the sorted walk gives low-discrepancy proportional
        ! coverage of the view sphere -- no colatitude band missed by ordering luck. Key: 48
        ! equal-area colatitude bands (1-cos theta), azimuth spreading within each band.
        ! Measured on the polar branch: +5.4 pts GT basis capture on preferred-orientation
        ! data, nothing on isotropic. Selection happens ONCE per stage; per-iteration
        ! resampling remains measured-negative (this is not an incremental EM).
        vstrat = 0
        call cov_env_int('SIMPLE_COV_OSTRAT', vstrat)
        l_ostrat = vstrat > 0
        allocate(worder(nptcls))
        do i = 1, nptcls
            worder(i) = i
        end do
        if( l_ostrat )then
            allocate(okey(nptcls))
            do i = 1, nptcls
                th_os = build%spproj_field%e2get(pinds(i))
                ph_os = build%spproj_field%e1get(pinds(i))
                okey(i) = real(int((1.0 - cos(th_os*PI/180.0))*24.0))*1000.0 + ph_os + 180.0
            end do
            call hpsort(okey, worder)
            deallocate(okey)
            write(logfhandle,'(A,A,A)') '>>> FLEX_PCA ',trim(label), &
                &' subsample: ORIENTATION-STRATIFIED walk (48 equal-area colatitude bands)'
            call flush(logfhandle)
        endif
        allocate(spinds(nptcls))
        nsel = 0
        do ihalf = 0, 1
            n_half = 0
            do i = 1, nptcls
                if( build%spproj_field%get_eo(pinds(i)) == ihalf ) n_half = n_half + 1
            end do
            if( n_half < 1 ) cycle
            ! split the per-part budget evenly between halfsets, never starving one
            ntgt  = min(n_half, max(1, (nmax_part + 1 - ihalf)/2))
            nkept = 0
            do ii = 1, nptcls
                i = worder(ii)
                if( build%spproj_field%get_eo(pinds(i)) /= ihalf ) cycle
                if( l_stride )then
                    if( mod(nkept, stride) == 0 )then
                        nsel = nsel + 1
                        spinds(nsel) = pinds(i)
                    endif
                else
                    ! real(dp) rather than integer products: nkept*ntgt overflows int32 at these sizes
                    if( int(real(nkept+1,dp)*real(ntgt,dp)/real(n_half,dp)) > &
                       &int(real(nkept,  dp)*real(ntgt,dp)/real(n_half,dp)) )then
                        nsel = nsel + 1
                        spinds(nsel) = pinds(i)
                    endif
                endif
                nkept = nkept + 1
            end do
        end do
        deallocate(worder)
        if( nsel < 2 ) THROW_HARD('stage subsample left too few particles; raise '//trim(env_max))
        if( l_ostrat .and. vstrat >= 2 )then
            ! SIMPLE_COV_OSTRAT=2: KEEP a stratified order instead of restoring project order,
            ! so that under ONLINE windows every contiguous window is view-balanced AND
            ! eo-balanced: deal the selection round-robin over (halfset x colatitude band).
            ! Costs scattered reads; on preferred-orientation data (cFAR ~0.06) the window
            ! composition is what the basis sees each iteration, and acquisition-ordered
            ! windows inherit micrograph-correlated view skew.
            block
                integer, allocatable :: sp2(:), bkey(:), bcnt(:), boff(:), bcur(:), bmem(:)
                integer :: kdeal, nout, ib2, nb2, maxcnt
                real    :: th2
                allocate(sp2(nsel), bkey(nsel))
                do i = 1, nsel
                    th2 = build%spproj_field%e2get(spinds(i))
                    bkey(i) = int((1.0 - cos(th2*PI/180.0))*24.0)               ! 48 equal-area bands
                    bkey(i) = bkey(i)*2 + build%spproj_field%get_eo(spinds(i))  ! x eo
                end do
                nb2 = maxval(bkey) + 1
                allocate(bcnt(0:nb2-1), boff(0:nb2-1), bcur(0:nb2-1), bmem(nsel))
                bcnt = 0
                do i = 1, nsel
                    bcnt(bkey(i)) = bcnt(bkey(i)) + 1
                end do
                boff(0) = 0
                do ib2 = 1, nb2 - 1
                    boff(ib2) = boff(ib2-1) + bcnt(ib2-1)
                end do
                bcur = 0
                do i = 1, nsel
                    ib2 = bkey(i)
                    bcur(ib2) = bcur(ib2) + 1
                    bmem(boff(ib2) + bcur(ib2)) = spinds(i)
                end do
                maxcnt = maxval(bcnt)
                nout = 0
                do kdeal = 1, maxcnt
                    do ib2 = 0, nb2 - 1
                        if( kdeal <= bcnt(ib2) )then
                            nout = nout + 1
                            sp2(nout) = bmem(boff(ib2) + kdeal)
                        endif
                    end do
                end do
                spinds(:nsel) = sp2(:nsel)
                deallocate(sp2, bkey, bcnt, boff, bcur, bmem)
                write(logfhandle,'(A,A,A,I0,A)') '>>> FLEX_PCA ',trim(label), &
                    &' subsample: STRATIFIED ORDER retained (', nb2, ' eo x band buckets, round-robin)'
                call flush(logfhandle)
            end block
        else
            call hpsort(spinds(:nsel))   ! restore project order so batched image reads stay sequential
        endif
        if( nsel < nptcls )then
            write(logfhandle,'(A,A,A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA ',trim(label),' subsample (', &
                &merge(1,0,l_stride),'=stride): using ',nsel,' of ',nptcls,' particles'
            call flush(logfhandle)
        endif
    end subroutine cov_stage_subsample

    !>  Write the two half-data latent solves so the per-particle error model can be calibrated
    !!  against the error actually observed. The halves are disjoint checkerboard subsets of one
    !!  particle's own Fourier samples (see cov_herm_inner's `half` argument), so var(z1 - z2)
    !!  measures that particle's estimation error directly, including the model misspecification the
    !!  analytic posterior covariance cannot express. `prior` is written alongside because the half
    !!  solves are shrunk by it. Gated on SIMPLE_COV_ZHALF.
    subroutine write_zhalf_replicates( zhalf, prior, nptcls, ncomp )
        real(dp), intent(in) :: zhalf(nptcls,ncomp,2), prior(ncomp)
        integer,  intent(in) :: nptcls, ncomp
        integer :: enable, funit, io_stat
        enable = 0
        call cov_env_int('SIMPLE_COV_ZHALF', enable)
        if( enable <= 0 ) return
        call fopen(funit, file=string('flex_pca_zhalf.bin'), access='stream', action='WRITE', &
            &status='REPLACE', iostat=io_stat)
        if( io_stat /= 0 )then
            THROW_WARN('could not open flex_pca_zhalf.bin; skipping half-solve export')
            return
        endif
        write(funit) nptcls, ncomp
        write(funit) zhalf
        write(funit) prior
        call fclose(funit)
        write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA wrote flex_pca_zhalf.bin (nptcls=', &
            &nptcls,' ncomp=',ncomp,')'
        call flush(logfhandle)
    end subroutine write_zhalf_replicates

    !>  Override an integer from the environment, if the variable is set and parses.
    subroutine cov_env_int_pub( name, val )
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: val
        call cov_env_int(name, val)
    end subroutine cov_env_int_pub

    !> Is an environment variable present at all? Used where the DEFAULT must be "behave exactly as
    !! before", not "behave as if the variable were at its documented default".
    logical function cov_env_is_set( name )
        character(len=*), intent(in) :: name
        character(len=32) :: envval
        integer :: stat, ln
        call get_environment_variable(name, envval, ln, stat)
        cov_env_is_set = (stat == 0 .and. ln >= 1)
    end function cov_env_is_set

    subroutine cov_env_int( name, val )
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: val
        character(len=32) :: envval
        integer :: stat, ln, ival
        call get_environment_variable(name, envval, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        read(envval(:ln), *, iostat=stat) ival
        if( stat == 0 .and. ival > 0 )then
            val = ival
            write(logfhandle,'(A,A,A,I0)') '>>> FLEX_PCA ',trim(name),' override: ',ival
            call flush(logfhandle)
        endif
    end subroutine cov_env_int


    !>  Packed accumulation + matrix-free CG is the DEFAULT reduced solve, so an UNSET
    !!  SIMPLE_COV_CGSOLVE selects packed. SIMPLE_COV_CGSOLVE=0 is the documented escape hatch back to
    !!  the dense d^2 x d^2 accumulator + Cholesky; cov_env_int cannot express that (it ignores every
    !!  value <= 0), which is why this goes through the presence-and-zero test in cov_env_int_off.
    !!  EVERY site whose memory model depends on the choice must call this -- if the dimension budget
    !!  and the solve disagree, d_tilde is sized against an accumulator that is never allocated.
    logical function cov_packed_cgsolve() result( packed )
        packed = .not. cov_env_int_off('SIMPLE_COV_CGSOLVE')
    end function cov_packed_cgsolve

    !>  Bytes in the reduced solve's ONE shared accumulator at column dimension d.
    pure real(dp) function cov_accum_bytes( d, packed ) result( nbytes )
        integer, intent(in) :: d
        logical, intent(in) :: packed
        real(dp) :: n
        if( packed )then
            n = real(d,dp)*real(d+1,dp)/2.d0   ! Mspk(npk,npk), npk = d(d+1)/2
        else
            n = real(d,dp)**2                  ! A(d^2,d^2)
        endif
        nbytes = 8.d0*n*n
    end function cov_accum_bytes

    !>  Largest d whose accumulator fits COV_ATHR_BUDGET under the model the solve will ACTUALLY use,
    !!  i.e. cov_accum_bytes(d, packed) <= COV_ATHR_BUDGET.
    pure integer function cov_dim_budget( packed ) result( d )
        logical, intent(in) :: packed
        if( packed )then
            ! d(d+1)/2 = sqrt(BUDGET/8)  =>  d = (-1 + sqrt(1 + 8*sqrt(BUDGET/8)))/2
            d = max(1, int((-1.d0 + sqrt(1.d0 + 8.d0*sqrt(COV_ATHR_BUDGET/8.d0)))/2.d0))
        else
            d = max(1, int((COV_ATHR_BUDGET/8.d0)**0.25d0))
        endif
    end function cov_dim_budget

    !>  One-shot read of the SIMPLE_COV_RKW override for COV_RIGHT_KERNEL_W.

    !> begin the device prep lifecycle for a stage batch loop when the GPU is present and
    !! enabled; the shared prep funnel (prep_imgs4projected_model) then takes its device
    !! branch for every batch. No-op (l_on=.false.) when an outer lifecycle already owns it.
    subroutine cov_dev_prep_start( params, build, l_on )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        logical,           intent(out)   :: l_on
        integer :: vprep
        l_on = .false.
        if( .not. flex_gpu_available() ) return
        if( flex_gpu_prep_ready() )      return   ! outer owner
        ! OPT-IN (SIMPLE_COV_GPU_PREP_STAGES=1): measured on the 4-worker Ribosembly arm, the
        ! fetch+unpack funnel is slower than the 6.4 s threaded CPU prep at these stages (8.4 s
        ! under device contention), and its 1e-6-level numerics nudged the probe convergence
        ! rule from 3 to 5 rounds (+31 s wall). The probe's own resident prep (no fetch) and
        ! the STATEREC resident hand-off are the configurations that pay.
        vprep = 0
        call cov_env_int('SIMPLE_COV_GPU_PREP_STAGES', vprep)
        if( vprep <= 0 ) return
        if( params%l_ml_reg )then
            ! whitening path needs the loaded sigma2 spectra (CPU prep THROWs without them)
            if( .not. allocated(build%esig%sigma2_noise) ) return
        endif
        if( cov_image_mask_radius(params) > 0. ) return   ! mask variant stays on the CPU
        call flex_gpu_prep_begin_f(build%lmsk, params%box, params%boxpd, MAXIMGBATCHSZ, &
            &0.0, .true.)
        l_on = .true.
    end subroutine cov_dev_prep_start

    subroutine cov_dev_prep_stop( l_on )
        logical, intent(in) :: l_on
        if( l_on ) call flex_gpu_prep_free_f
    end subroutine cov_dev_prep_stop

    subroutine cov_init_right_kernel_width
        character(len=32) :: envval
        integer :: stat, ln
        real    :: rval
        if( cov_rkw_read ) return
        cov_rkw_read = .true.
        call get_environment_variable('SIMPLE_COV_RKW', envval, ln, stat)
        if( stat == 0 .and. ln > 0 )then
            read(envval(:ln), *, iostat=stat) rval
            if( stat == 0 ) COV_RIGHT_KERNEL_W = rval
        endif
        if( COV_RIGHT_KERNEL_W > 0. )then
            write(logfhandle,'(A,F6.3)') '>>> FLEX_PCA right (column-gather) kernel: triangular, width ', &
                &COV_RIGHT_KERNEL_W
        else
            write(logfhandle,'(A)') '>>> FLEX_PCA right (column-gather) kernel: shared KB backprojection stencil'
        endif
        call flush(logfhandle)
    end subroutine cov_init_right_kernel_width

    !> Sampling precision of the MAP latent estimate, Q = A*Gtil^+*A with A = Gtil + diag(prior). This is
    !! the precision of the ESTIMATOR z_hat, not the posterior precision A, so distances measured with it
    !! reflect how well each component was actually determined for the particle.
    subroutine map_sampling_precision( Gtil, prior, n, Qout )
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: Gtil(n,n), prior(n)
        real(dp), intent(out) :: Qout(n,n)
        real(dp) :: Amat(n,n), Gpinv(n,n), Vmat(n,n), Awork(n,n), ev(n), thresh
        integer  :: ii, jj, kk, nrot
        Amat = Gtil
        do ii = 1, n
            Amat(ii,ii) = Amat(ii,ii) + prior(ii)
        end do
        Awork = Gtil
        call jacobi(Awork, n, n, ev, Vmat, nrot)   ! symmetric eigendecomposition (LAPACK dsyev)
        thresh = COV_PINV_RCOND * maxval(abs(ev))
        Gpinv  = 0.d0
        do kk = 1, n
            if( abs(ev(kk)) <= thresh ) cycle      ! drop the null space, as pinv does
            do jj = 1, n
                do ii = 1, n
                    Gpinv(ii,jj) = Gpinv(ii,jj) + Vmat(ii,kk)*Vmat(jj,kk)/ev(kk)
                end do
            end do
        end do
        Qout = matmul(Amat, matmul(Gpinv, Amat))
        Qout = 0.5d0*(Qout + transpose(Qout))      ! symmetrise away round-off
    end subroutine map_sampling_precision

    !> Full column-covariance eigenbasis pipeline.
    subroutine build_covariance_eigenbasis( params, build, mean_rec, pinds, nptcls, &
        &col_sep, neigs_req, basis_recs, eigvals, ncomp_out, sig2_out, basis_imgs, fprefix )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, col_sep, neigs_req
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        real(dp),            allocatable, intent(out) :: eigvals(:)
        integer,             intent(out)   :: ncomp_out
        real(dp),            intent(out)   :: sig2_out
        !> optional clean real-space eigenvolumes + output-name prefix, used by the
        !! held-out (cross-halfset) embedding to align two independently fitted bases
        type(image), allocatable, optional, intent(out) :: basis_imgs(:)
        character(len=*),        optional, intent(in)  :: fprefix
        integer,   allocatable :: col_hkl(:,:), col_lookup(:,:,:)
        complex,   allocatable :: Bcol_e(:,:,:,:), Bcol_o(:,:,:,:), colvol(:,:,:,:)
        real,      allocatable :: Hcol_e(:,:,:,:), Hcol_o(:,:,:,:), col_fsc(:)
        type(image), allocatable :: realvols(:), utilde_real(:)
        type(reconstructor) :: work
        type(reconstructor), allocatable :: utilde(:)
        !> capped, halfset-balanced particle subset the column-subspace initialiser runs on. Equals
        !! pinds when the cap is off, which is the default.
        integer, allocatable :: bpinds(:)
        integer :: nbp, nparts_sub
        !> reproducibility-weighted column pruning (SIMPLE_COV_COLFSC_W / _MIN)
        integer :: vcolw, vcolmin, s_keep, s_src, ngood_col
        real    :: colfsc_thresh, colfsc_w
        !> probe-worker handoff, read back from the master's flex_pca_probe.txt
        real(dp),            allocatable :: eig_probe(:)
        real(dp),            allocatable :: zw(:,:), contrastw(:), precw(:,:,:), rew(:), rmew(:)
        real(dp) :: sig2_probe
        integer  :: ncomp_probe
        real(dp), allocatable :: vred(:,:)
        integer :: ncol, nreal, s, lb(3), ub(3), nyq_rec, d_tilde, q, directsvd
        !> optional solve-specific subset (SIMPLE_COV_SOLVE_MAX/_STRIDE); defaults to bpinds
        integer, allocatable :: spinds_solve(:)
        integer :: nsolve
        integer(timer_int_kind) :: t_blk
        logical :: l_cols_ok
        real(dp), allocatable :: svals(:)
        ! one work reconstructor defines the expanded lattice / Nyquist / grid correction
        call init_column_reconstructor(params, build, work)
        lb      = lbound(work%cmat_exp)
        ub      = ubound(work%cmat_exp)
        nyq_rec = work%get_lfny(1)
        ! NOTE: do not restrict the column sample sweep to the lp band. Everything above the lp shell is
        ! hard-zeroed by realize_hermitian_volume, yet it still moves the eigenvalues by up to 1.2e-4
        ! relative for only ~1.09x -- some path out of band escapes that argument. Find it first.
        ! column selection. A WORKER in the COLS round must not re-run the SNR-greedy selection: its
        ! own estimate_snr_volume call only produced this part's contribution, so the selection would
        ! see a fraction of the variance. It reads the master's choice instead. A worker in the SNR
        ! round runs the selection path only to reach estimate_snr_volume, and exits below.
        ! PROBE WORKER. Like the SOLVE worker, nothing upstream is needed: the master refreshed
        ! flex_pca_pc*.mrc and flex_pca_probe.txt before scheduling this round, so the worker rebuilds
        ! the current basis from disk and contributes one EM half-pass over its own particle range.
        ! niters=1 -- the iteration loop lives on the master, one qsys round per iteration, because
        ! the basis the E-step projects against changes every iteration.
        if( flex_pca_is_worker() .and. params%stage == PCA_STAGE_PROBE )then
            call load_probe_state(ncomp_probe, eig_probe, sig2_probe)
            call load_probe_basis(params, build, ncomp_probe, basis_recs)
            call probe_subspace_iteration(params, build, mean_rec, basis_recs, eig_probe, sig2_probe, &
                &pinds, nptcls, ncomp_probe, 1)
            do s = 1, size(basis_recs)
                call basis_recs(s)%dealloc_rho; call basis_recs(s)%kill
            end do
            deallocate(basis_recs)
            if( allocated(eig_probe) ) deallocate(eig_probe)
            if( .not. allocated(eigvals) ) allocate(eigvals(0))
            allocate(basis_recs(0))
            ncomp_out = 0
            sig2_out  = 0._dp
            call work%dealloc_rho; call work%kill
            return
        endif
        ! embed worker: same handoff as the probe worker (basis on disk as flex_pca_pc*.mrc,
        ! dimension/prior variances/noise level in flex_pca_probe.txt) but only one round, since the
        ! basis is final by now and does not change under the workers
        if( flex_pca_is_worker() .and. params%stage == PCA_STAGE_EMBED )then
            call load_probe_state(ncomp_probe, eig_probe, sig2_probe)
            call load_probe_basis(params, build, ncomp_probe, basis_recs)
            allocate(zw(nptcls,ncomp_probe), contrastw(nptcls), precw(ncomp_probe,ncomp_probe,nptcls))
            allocate(rew(nptcls), rmew(nptcls))
            call embed_latents_with_contrast(params, build, mean_rec, basis_recs, ncomp_probe, &
                &eig_probe, sig2_probe, pinds, nptcls, zw, contrastw, precw, rew, rmew, &
                &stats_only=.true.)
            deallocate(zw, contrastw, precw, rew, rmew)
            do s = 1, size(basis_recs)
                call basis_recs(s)%dealloc_rho; call basis_recs(s)%kill
            end do
            deallocate(basis_recs)
            if( allocated(eig_probe) ) deallocate(eig_probe)
            if( .not. allocated(eigvals) ) allocate(eigvals(0))
            allocate(basis_recs(0))
            ncomp_out = 0
            sig2_out  = 0._dp
            call work%dealloc_rho; call work%kill
            return
        endif
        ! ---- EM PATH (the only estimator) ----
        ! Hand probe_subspace_iteration a data-free basis. The moment/covariance estimator that used
        ! to live here -- SNR volume, column selection and accumulation, reduced solve -- has been
        ! removed; probe_subspace_iteration is itself a PPCA EM and supersedes it. Only reachable on
        ! the master: a worker in a PROBE or EMBED round has already returned above with the basis it
        ! read off disk, so nothing here runs twice.
        ! basis_imgs is the PRE-refinement basis for the held-out/bagged arms. Under EM that basis is
        ! a data-free placeholder and returning it would silently compare initialisers instead of
        ! fits, so refuse rather than mislead. Two-fit reproducibility is measured by running two
        ! jobs and comparing their written eigenvolumes.
        if( present(basis_imgs) ) THROW_HARD('the held-out/bagged basis arms need a pre-refinement &
            &basis, which the EM initialiser cannot provide; run two jobs and compare flex_pca_pc*.mrc')
        call init_basis_datafree(params, build, mean_rec, pinds, nptcls, col_sep, neigs_req, &
            &basis_recs, eigvals, ncomp_out, sig2_out)
        call work%dealloc_rho; call work%kill
    end subroutine build_covariance_eigenbasis

    !> DATA-FREE EM INITIALISER (SIMPLE_COV_EM=1).
    !!
    !! Replaces the whole moment estimator -- SNR volume, column selection, column accumulation,
    !! merge, reduced solve -- as the thing that hands probe_subspace_iteration its starting basis.
    !! probe_subspace_iteration is ALREADY a PPCA EM (posterior precision A = (a^2/sig2)G + Gamma^-1,
    !! E[zz'] = zz' + A^-1, coupled per-voxel M-step, Gamma from the posterior second moment). What it
    !! has never been given is a start that is not the covariance eigenbasis, which is why running it
    !! today is evidence for neither the EM formulation nor against it: if the columns are noise, the
    !! probe starts inside a noise subspace and its Gamma^(0) are noise variances.
    !!
    !! The subspace here is the lowest-frequency Fourier modes the band admits, realised as cos/sin
    !! Friedel pairs through exactly the same band-limit + mask + deapodisation the column path uses,
    !! then orthonormalised. NO data moment is formed anywhere: select_covariance_columns_lowfreq is
    !! pure geometry (greedy smallest-|xi|, col_sep-separated) and the "columns" fed to it are unit
    !! impulses, not estimates. Deterministic, so two arms are comparable run to run.
    subroutine init_basis_datafree( params, build, mean_rec, pinds, nptcls, col_sep, neigs_req, &
        &basis_recs, eigvals, ncomp_out, sig2_out )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, col_sep, neigs_req
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        real(dp),            allocatable, intent(out) :: eigvals(:)
        integer,             intent(out)   :: ncomp_out
        real(dp),            intent(out)   :: sig2_out
        type(reconstructor) :: work
        type(reconstructor), allocatable :: utilde(:)
        type(image),         allocatable :: realvols(:), utilde_real(:)
        integer,             allocatable :: col_hkl(:,:)
        complex,             allocatable :: colvol(:,:,:,:)
        real(dp),            allocatable :: svals(:)
        integer  :: lb(3), ub(3), ncols_req, ncol, nreal, d_tilde, s, q, h, k, l
        integer, allocatable :: cpinds(:)
        integer  :: ncal
        real(dp) :: gam0
        type(string) :: fname
        call init_column_reconstructor(params, build, work)
        lb = lbound(work%cmat_exp)
        ub = ubound(work%cmat_exp)
        ! each impulse yields a cos AND a sin volume, so half as many impulses as components are
        ! needed; a small margin covers the pairs orthonormalisation drops at the energy floor
        ncols_req = max(1, (neigs_req + 1)/2 + 2)
        call select_covariance_columns_lowfreq(params, ncols_req, col_sep, col_hkl, ncol)
        allocate(colvol(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol), source=cmplx(0.,0.))
        do s = 1, ncol
            h = col_hkl(1,s); k = col_hkl(2,s); l = col_hkl(3,s)
            if( h < lb(1) .or. h > ub(1) .or. k < lb(2) .or. k > ub(2) &
                &.or. l < lb(3) .or. l > ub(3) ) cycle
            colvol(h,k,l,s) = cmplx(1.,0.)
        end do
        call columns_to_real_representatives(params, work, colvol, ncol, lb, ub, realvols, nreal)
        deallocate(colvol, col_hkl)
        if( nreal < 1 ) THROW_HARD('flex_pca EM initialiser produced no basis representatives')
        call orthonormalize_representatives(params, build, realvols, nreal, utilde, utilde_real, &
            &d_tilde, svals, nptcls_basis=nptcls)
        do s = 1, nreal
            call realvols(s)%kill
        end do
        deallocate(realvols)
        ncomp_out = max(1, min(neigs_req, d_tilde))
        call basis_recs_from_images(params, build, utilde_real(1:ncomp_out), ncomp_out, basis_recs)
        ! The eigenvolume MRCs are the master->worker handoff for every distributed probe round; on
        ! the moment path form_eigenbasis_from_reduced writes them, and nothing else does, so the EM
        ! path has to write its own initial basis or the first PROBE round finds no flex_pca_pc001.mrc.
        if( flex_pca_is_master() .or. .not. flex_pca_is_worker() )then
            do q = 1, ncomp_out
                fname = string('flex_pca_pc')//int2str_pad(q,3)//MRC_EXT
                call utilde_real(q)%write(fname, del_if_exists=.true.); call fname%kill
            end do
        endif
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA EM INIT (data-free): impulses=',ncol, &
            &'  representatives=',nreal,'  rank=',ncomp_out
        call flush(logfhandle)
        ! the ONLY data pass: the noise convention constant and the initial prior variance.
        ! Two scalars do not need a million particles, and at box_crop=96 on the full dataset an
        ! uncapped pass costs more than several EM iterations, so it shares the initialiser's
        ! budget SIMPLE_COV_CALIB_MAX (default COV_CALIB_MAX_PTCLS; SIMPLE_COV_BASIS_MAX used to
        ! be shared here but its default of 0 left this pass uncapped). Always master here -- a worker returned
        ! long before this point -- so nparts is 1 and the budget is not divided twice.
        call cov_stage_subsample(build, pinds, nptcls, 1, COV_CALIB_MAX_PTCLS, &
            &'SIMPLE_COV_CALIB_STRIDE', 'SIMPLE_COV_CALIB_MAX', 'EM CALIBRATION', cpinds, ncal)
        call em_calibrate_noise_prior(params, build, mean_rec, basis_recs, ncomp_out, cpinds, ncal, &
            &sig2_out, gam0)
        deallocate(cpinds)
        allocate(eigvals(ncomp_out), source=gam0)
        do s = 1, d_tilde
            call utilde(s)%dealloc_rho; call utilde(s)%kill
            call utilde_real(s)%kill
        end do
        deallocate(utilde, utilde_real)
        if( allocated(svals) ) deallocate(svals)
        call work%dealloc_rho; call work%kill
    end subroutine init_basis_datafree

    !> One pass over the particles to fix the two scalars the EM cannot invent: the whitened-noise
    !! convention constant sig2 (measured on the high-frequency shells, where conformational signal is
    !! negligible) and the initial prior variance Gamma^(0).
    !!
    !! Gamma^(0) is deliberately set from the mean-deflated residual with only the noise floor removed,
    !! i.e. an OVER-estimate of the conformational signal, which makes the iteration-1 prior WEAK. The
    !! asymmetry is on purpose: a prior that starts too tight shrinks z, which shrinks the Gamma update,
    !! which tightens the prior -- the self-reinforcing collapse that is the single most plausible cause
    !! of the previous EM's failure. Starting loose is corrected by the very first Gamma M-step; starting
    !! tight is not corrected by anything.
    subroutine em_calibrate_noise_prior( params, build, mean_rec, basis_recs, ncomp, pinds, nptcls, &
        &sig2_out, gam0 )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), intent(inout) :: basis_recs(:)
        integer,             intent(in)    :: ncomp, pinds(:), nptcls
        real(dp),            intent(out)   :: sig2_out, gam0
        type(fplane_type), allocatable :: fpls(:), basis_fpls(:,:), mean_fpl(:)
        type(ori),         allocatable :: orientations(:)
        real(dp), allocatable :: res_thr(:), trg_thr(:), hfp_thr(:), hfc_thr(:), aa_thr(:)
        integer,  allocatable :: nval_thr(:)
        integer  :: nthr, ithr, i, q, ibatch, batchlims(2), batchsz, nyq_rec, nval
        real(dp) :: a, aa, e_mm, e_yy, myv, res, trg, pw, cnt, wcnt
        real(dp) :: res_sum, trg_sum, aa_sum, hfpw, hfcnt
        nthr    = nthr_glob
        nyq_rec = mean_rec%get_lfny(1)
        call mean_rec%expand_exp
        do q = 1, ncomp
            call basis_recs(q)%expand_exp
        end do
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        allocate(mean_fpl(nthr), basis_fpls(ncomp,nthr), orientations(MAXIMGBATCHSZ))
        allocate(res_thr(nthr), trg_thr(nthr), hfp_thr(nthr), hfc_thr(nthr), aa_thr(nthr), source=0.d0)
        allocate(nval_thr(nthr), source=0)
        wcnt = 0.d0
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz), mskrad=cov_image_mask_radius(params))
            do i = 1, batchsz
                call build%spproj_field%get_ori(pinds(batchlims(1)+i-1), orientations(i))
            end do
            if( wcnt <= 0.d0 ) call cov_herm_selfpower(fpls(1), pw, wcnt)
            !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
            !$omp& private(i,ithr,q,a,aa,e_mm,e_yy,myv,res,trg,pw,cnt)
            do i = 1, batchsz
                if( orientations(i)%isstatezero() ) cycle
                ithr = omp_get_thread_num() + 1
                call project_fplanes_mean_basis(mean_rec, basis_recs, orientations(i), fpls(i), &
                    &mean_fpl(ithr), basis_fpls(:,ithr), apply_ctf_amp=.true.)
                e_mm = real(cov_herm_inner(mean_fpl(ithr), mean_fpl(ithr)), dp)
                myv  = real(cov_herm_inner(mean_fpl(ithr), fpls(i)), dp)
                e_yy = real(cov_herm_inner(fpls(i), fpls(i)), dp)
                a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
                aa   = a*a
                res  = max(e_yy - 2.d0*a*myv + aa*e_mm, 0.d0)
                trg  = 0.d0
                do q = 1, ncomp
                    trg = trg + real(cov_herm_inner(basis_fpls(q,ithr), basis_fpls(q,ithr)), dp)
                end do
                call plane_hf_power(fpls(i), nyq_rec, 0.7, pw, cnt)
                res_thr(ithr)  = res_thr(ithr)  + res
                trg_thr(ithr)  = trg_thr(ithr)  + trg
                aa_thr(ithr)   = aa_thr(ithr)   + aa
                hfp_thr(ithr)  = hfp_thr(ithr)  + pw
                hfc_thr(ithr)  = hfc_thr(ithr)  + cnt
                nval_thr(ithr) = nval_thr(ithr) + 1
            end do
            !$omp end parallel do
        end do
        res_sum = sum(res_thr); trg_sum = sum(trg_thr); aa_sum = sum(aa_thr)
        hfpw    = sum(hfp_thr); hfcnt   = sum(hfc_thr); nval = sum(nval_thr)
        if( nval < 1 ) THROW_HARD('flex_pca EM calibration saw no valid particles')
        sig2_out = max(hfpw / max(hfcnt, 1.d0), DTINY)
        ! signal power = mean-deflated residual with the noise floor removed; floor it at a tenth of
        ! the residual so a mis-measured sig2 cannot drive the prior to zero
        res_sum = max(res_sum/real(nval,dp) - sig2_out*wcnt, 0.1d0*res_sum/real(nval,dp))
        trg_sum = max(trg_sum/real(nval,dp), DTINY)
        aa_sum  = max(aa_sum /real(nval,dp), DTINY)
        gam0    = max(res_sum / (aa_sum * trg_sum), DTINY)
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,I0)') '>>> FLEX_PCA EM CALIBRATION: sig2=',sig2_out, &
            &'  gamma0=',gam0,'  particles=',nval
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,F6.3)') '>>>   mean deflated signal=',res_sum, &
            &'  mean tr(G)=',trg_sum,'  mean a^2=',real(aa_sum)
        call flush(logfhandle)
        do ithr = 1, nthr
            call cleanup_plane(mean_fpl(ithr))
            do q = 1, ncomp
                call cleanup_plane(basis_fpls(q,ithr))
            end do
        end do
        do i = 1, size(orientations)
            call orientations(i)%kill
        end do
        call cleanup_rec_buffers(build, fpls)
        deallocate(mean_fpl, basis_fpls, orientations)
        deallocate(res_thr, trg_thr, hfp_thr, hfc_thr, aa_thr, nval_thr)
    end subroutine em_calibrate_noise_prior





    !>  Master -> probe-worker handoff: basis dimension, prior variances, whitened-noise level.
    subroutine save_probe_state( ncomp, eigvals, sig2_eff )
        integer,  intent(in) :: ncomp
        real(dp), intent(in) :: eigvals(:), sig2_eff
        integer :: funit, io_stat, q
        call fopen(funit, file=string(COV_PROBE_META), action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('save_probe_state', io_stat)
        write(funit,*) ncomp
        write(funit,*) sig2_eff
        do q = 1, ncomp
            write(funit,*) eigvals(q)
        end do
        call fclose(funit)
    end subroutine save_probe_state

    subroutine load_probe_state( ncomp, eigvals, sig2_eff )
        integer,               intent(out) :: ncomp
        real(dp), allocatable, intent(out) :: eigvals(:)
        real(dp),              intent(out) :: sig2_eff
        integer :: funit, io_stat, q
        if( .not. file_exists(string(COV_PROBE_META)) ) &
            &THROW_HARD('flex_pca probe worker found no '//COV_PROBE_META//' from the master')
        call fopen(funit, file=string(COV_PROBE_META), action='READ', status='OLD', iostat=io_stat)
        call fileiochk('load_probe_state', io_stat)
        read(funit,*) ncomp
        read(funit,*) sig2_eff
        if( ncomp < 1 ) THROW_HARD('invalid cached probe basis dimension')
        allocate(eigvals(ncomp))
        do q = 1, ncomp
            read(funit,*) eigvals(q)
        end do
        call fclose(funit)
    end subroutine load_probe_state

    !>  Rebuild the projection-ready basis a probe worker needs from the master's flex_pca_pc*.mrc.
    !!  Same idiom as load_utilde_stack and probe_external_basis: set_rmat then fft then expand_exp,
    !!  never add(), which would leave the reconstructor flagged Fourier and propagate an
    !!  untransformed grid.
    subroutine load_probe_basis( params, build, ncomp, basis_recs )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        integer,             intent(in)    :: ncomp
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        type(image)  :: vol
        type(string) :: fname
        integer      :: q
        allocate(basis_recs(ncomp))
        call vol%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        do q = 1, ncomp
            fname = string('flex_pca_pc')//int2str_pad(q,3)//MRC_EXT
            if( .not. file_exists(fname) ) &
                &THROW_HARD('flex_pca probe worker found no '//fname%to_char()//' from the master')
            call vol%read(fname)
            call init_column_reconstructor(params, build, basis_recs(q))
            call basis_recs(q)%set_rmat(vol%get_rmat(), .false.)
            call basis_recs(q)%fft
            call basis_recs(q)%expand_exp
            call fname%kill
        end do
        call vol%kill
    end subroutine load_probe_basis



    !> Deterministic low-frequency column selection: repeatedly take the lowest-|xi| candidate in the
    !! canonical Hermitian half that is at least col_sep away from every already-chosen column.
    subroutine select_covariance_columns_lowfreq( params, ncols_req, col_sep, col_hkl, ncol )
        class(parameters),    intent(in)  :: params
        integer,              intent(in)  :: ncols_req, col_sep
        integer, allocatable, intent(out) :: col_hkl(:,:)
        integer,              intent(out) :: ncol
        integer, allocatable :: cand(:,:)
        integer :: kfromto(2), kmax, kmax_sq, h, k, l, sep, r_sq, ncand, i, target
        kfromto = covariance_kfromto(params)
        kmax    = max(2, kfromto(2))
        kmax_sq = kmax*(kmax+1)
        sep     = max(1, col_sep)
        target  = max(1, ncols_req)
        allocate(cand(3, (2*kmax+1)**3))
        ncand = 0
        do h = 0, kmax
            do k = -kmax, kmax
                do l = -kmax, kmax
                    r_sq = h*h + k*k + l*l
                    if( r_sq == 0 .or. r_sq > kmax_sq ) cycle
                    if( h == 0 )then
                        if( k < 0 ) cycle
                        if( k == 0 .and. l < 0 ) cycle
                    endif
                    ncand = ncand + 1
                    cand(:,ncand) = [h,k,l]
                end do
            end do
        end do
        allocate(col_hkl(3, target))
        ncol = 0
        do
            call pick_next_lowfreq(cand, ncand, col_hkl, ncol, sep, i)
            if( i == 0 ) exit
            ncol = ncol + 1
            col_hkl(:,ncol) = cand(:,i)
            cand(:,i) = huge(1)
            if( ncol >= target ) exit
        end do
        if( ncol < 1 ) THROW_HARD('flex_pca laid out no impulse columns for the EM initialiser; increase lp or neigs')
        deallocate(cand)
    end subroutine select_covariance_columns_lowfreq

    subroutine pick_next_lowfreq( cand, ncand, chosen, nchosen, sep, best )
        integer, intent(in)  :: cand(:,:), ncand, chosen(:,:), nchosen, sep
        integer, intent(out) :: best
        integer :: i, j, r_sq, best_r, d(3)
        logical :: ok
        best   = 0
        best_r = huge(1)
        do i = 1, ncand
            if( cand(1,i) == huge(1) ) cycle
            r_sq = sum(cand(:,i)**2)
            if( r_sq >= best_r ) cycle
            ok = .true.
            do j = 1, nchosen
                d = cand(:,i) - chosen(:,j)
                if( sum(d**2) < sep*sep )then
                    ok = .false.; exit
                endif
            end do
            if( ok )then
                best   = i
                best_r = r_sq
            endif
        end do
    end subroutine pick_next_lowfreq

    function covariance_kfromto( params ) result( kfromto )
        class(parameters), intent(in) :: params
        integer :: kfromto(2), kto_full
        real    :: dstep_crop
        kto_full   = max(1, fdim(params%box_crop) - 1)
        kfromto(1) = 1
        kfromto(2) = kto_full
        if( params%lp > 2.0 * params%smpd_crop + TINY )then
            dstep_crop = real(max(1, params%box_crop - 1)) * params%smpd_crop
            kfromto(2) = max(1, min(kto_full, int(dstep_crop / params%lp)))
        endif
    end function covariance_kfromto


    !> EXTERNAL-BASIS PROBE: embed the particles in a basis read from disk (the run's own eigenvolumes,
    !! optionally with extra probe volumes appended) and write the coefficients, so the embedding stage can
    !! be exercised against a known basis without re-fitting the covariance.
    subroutine probe_external_basis( params, build, mean_rec, pinds, nptcls, eigdir, neigs, eigvals, &
        &sig2_eff, probe_prefix, nprobe )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, neigs, nprobe
        character(len=*),    intent(in)    :: eigdir, probe_prefix
        real(dp),            intent(in)    :: eigvals(:), sig2_eff
        type(reconstructor), allocatable :: basis_recs(:)
        type(image)  :: vol
        type(string) :: fname
        real(dp), allocatable :: ev(:), z(:,:), contrast(:), precision(:,:,:)
        real(dp), allocatable :: resid_energy(:), resid_mean_energy(:), sorted(:)
        real     :: dummy
        integer  :: ncomb, k, u, i, q
        real(dp) :: evmed
        ! nprobe = 0 is legitimate: an appended probe volume can dominate the projected Gram spectrum, in
        ! which case the reported conditioning describes the probe rather than the basis.
        if( nprobe < 0 ) return
        ncomb = neigs + nprobe
        write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA EXTERNAL-BASIS PROBE: ', neigs, &
            &' eigenvolumes + ', nprobe, ' probe volumes in one joint solve'
        if( ncomb < 2 ) THROW_HARD('probe_external_basis: need at least 2 basis volumes')
        call flush(logfhandle)
        allocate(basis_recs(ncomb), ev(ncomb))
        call vol%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        do k = 1, ncomb
            if( k <= neigs )then
                fname = trim(eigdir)//'flex_pca_pc'//int2str_pad(k,3)//MRC_EXT
            else
                fname = trim(probe_prefix)//int2str_pad(k-neigs,3)//MRC_EXT
            endif
            if( .not. file_exists(fname%to_char()) )then
                write(logfhandle,'(A,A)') '>>> FLEX_PCA probe basis: missing ', fname%to_char()
                call fname%kill; call vol%kill
                do i = 1, k-1
                    call basis_recs(i)%dealloc_rho; call basis_recs(i)%kill
                end do
                deallocate(basis_recs, ev)
                return
            endif
            call vol%read(fname)
            call fname%kill
            if( params%msk_crop > TINY ) call vol%mask3D_soft(params%msk_crop, backgr=0.)
            ! the set_rmat + fft + expand_exp idiom -- NEVER add(), which silently yields a
            ! zero projected basis when the reconstructor is left flagged as Fourier
            call init_column_reconstructor(params, build, basis_recs(k))
            call basis_recs(k)%set_rmat(vol%get_rmat(), .false.)
            call basis_recs(k)%fft
            call basis_recs(k)%expand_exp
        end do
        call vol%kill
        allocate(sorted(neigs), source=eigvals(1:neigs))
        evmed = sorted(max(1,neigs/2))
        do k = 1, ncomb
            if( k <= neigs )then
                ev(k) = max(eigvals(k), DTINY)
            else
                ev(k) = max(evmed, DTINY)
            endif
        end do
        write(logfhandle,'(A,ES12.4)') '>>> FLEX_PCA probe prior variance (median eigenvalue): ', evmed
        allocate(z(nptcls,ncomb), contrast(nptcls), precision(ncomb,ncomb,nptcls), &
            &resid_energy(nptcls), resid_mean_energy(nptcls))
        call embed_latents_with_contrast(params, build, mean_rec, basis_recs, ncomb, ev, sig2_eff, &
            &pinds, nptcls, z, contrast, precision, resid_energy, resid_mean_energy)
        call del_file('flex_pca_probe_coordinates.txt')
        open(newunit=u, file='flex_pca_probe_coordinates.txt', status='replace', action='write')
        write(u,'(A)',advance='no') '# particle'
        do q = 1, neigs
            write(u,'(A,I0)',advance='no') ' pc', q
        end do
        do q = 1, nprobe
            write(u,'(A,I0)',advance='no') ' probe', q
        end do
        write(u,*)
        do i = 1, nptcls
            write(u,'(I10)',advance='no') pinds(i)
            do q = 1, ncomb
                write(u,'(1X,ES16.8)',advance='no') z(i,q)
            end do
            write(u,*)
        end do
        close(u)
        write(logfhandle,'(A)') '>>> FLEX_PCA probe coefficients -> flex_pca_probe_coordinates.txt'
        call flush(logfhandle)
        do k = 1, ncomb
            call basis_recs(k)%dealloc_rho; call basis_recs(k)%kill
        end do
        deallocate(basis_recs, ev, z, contrast, precision, resid_energy, resid_mean_energy, sorted)
        dummy = 0.
    end subroutine probe_external_basis

    subroutine init_column_reconstructor( params, build, rec )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: rec
        call rec%new([params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
        call rec%alloc_rho(params,build%spproj,expand=.true.)
        call rec%reset
        call rec%reset_exp
    end subroutine init_column_reconstructor

    !> Particle-image mask radius in pixels at params%box, or 0 to disable.
    function cov_image_mask_radius( params ) result( r )
        class(parameters), intent(in) :: params
        real :: r
        integer :: vmi
        ! Runtime override (SIMPLE_COV_MASK_IMAGES=1). The compile-time default is OFF, which is
        ! safe ONLY when the solvent is pure noise -- true for synthetic data, FALSE for real data,
        ! where the region outside the envelope carries ice-thickness gradients, neighbouring
        ! particles and carbon edges. Those are low-frequency AND reproducible between halfsets, so
        ! they enter the covariance as apparent signal and can dominate the leading eigenvectors.
        r = 0.
        vmi = 0
        call cov_env_int('SIMPLE_COV_MASK_IMAGES', vmi)
        if( .not. (COV_MASK_IMAGES .or. vmi > 0) ) return
        if( params%msk_crop <= 0. .or. params%box_crop <= 0 ) return
        r = COV_MASK_MARGIN * params%msk_crop * real(params%box) / real(params%box_crop)
        r = min(r, 0.5*real(params%box) - COSMSKHALFWIDTH - 1.)
    end function cov_image_mask_radius

    !> SIMPLE_COV_DUMP_ACC=1 dumps the raw stage accumulators (SNR var/dens, column B/H) and the
    !! column selection to the run dir, for byte-level A/B validation of threading restructures.
    !! Reuses the distributed part-file writers, so the dumps carry the versioned headers.
    logical function cov_dump_acc_on() result( on )
        character(len=32) :: envval
        integer :: stat, ln, ival
        on = .false.
        call get_environment_variable('SIMPLE_COV_DUMP_ACC', envval, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        read(envval(:ln), *, iostat=stat) ival
        if( stat == 0 ) on = ival /= 0
    end function cov_dump_acc_on

    !> Single entry point for the covariance mean.
    subroutine estimate_covariance_mean( params, build, mean_rec, pinds, nptcls )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
        integer :: vfc
        ! one-time runtime init of the fitted-contrast override (single-threaded here; every stage
        ! of every process passes through the mean before touching particles)
        vfc = 0
        call cov_env_int('SIMPLE_COV_CONTRAST', vfc)
        cov_fit_contrast_rt = vfc > 0
        if( cov_fit_contrast_rt )then
            write(logfhandle,'(A)') '>>> FLEX_PCA PER-PARTICLE CONTRAST ON (deviations against &
                &fitted a_i; kills the rank-one contrast mode -- RECOVAR A.4)'
            call flush(logfhandle)
        endif
        if( COV_MEAN_FROM_DATA )then
            call estimate_mean_from_data(params, build, mean_rec, pinds, nptcls)
        else
            call init_mean_reconstructor(params, build, mean_rec)
            if( flex_pca_is_worker() )then
                call apply_cached_mean_scale(params, mean_rec)
            else
                call estimate_mean_scale(params, build, mean_rec, pinds, nptcls)
            endif
        endif
    end subroutine estimate_covariance_mean

    !> Kernel-regression consensus mean (eq. S.1) estimated from the particles themselves, as an
    !! alternative to reading the supplied consensus volume. Selected by COV_MEAN_FROM_DATA.
    subroutine estimate_mean_from_data( params, build, mean_rec, pinds, nptcls )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
        type(fplane_type), allocatable :: fpls(:)
        type(fplane_type) :: num_fpl
        type(ori)    :: orientation
        type(image)  :: gridcorr_img
        type(string) :: fname
        integer :: batchlims(2), batchsz, ibatch, i, iptcl, used
        logical :: l_devprep
        integer(timer_int_kind) :: t_phase
        call init_column_reconstructor(params, build, mean_rec)
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        call cov_dev_prep_start(params, build, l_devprep)
        used    = 0
        t_phase = tic()
        write(logfhandle,'(A)') '>>> FLEX_PCA MEAN ESTIMATION (eq. S.1 kernel regression)'
        call flush(logfhandle)
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz), mskrad=cov_image_mask_radius(params))
            do i = 1, batchsz
                iptcl = pinds(batchlims(1)+i-1)
                call build%spproj_field%get_ori(iptcl, orientation)
                if( orientation%isstatezero() ) cycle
                call form_reconstruction_plane(fpls(i), num_fpl)
                call mean_rec%insert_plane_oversamp(build%pgrpsyms, orientation, num_fpl)
                used = used + 1
            end do
        end do
        call cov_dev_prep_stop(l_devprep)
        call orientation%kill
        call cleanup_plane(num_fpl)
        call cleanup_rec_buffers(build, fpls)
        if( used < 1 ) THROW_HARD('flex_pca mean estimation found no valid particles')
        ! canonical gridding finalization, identical to reconstruct3D_reference
        call mean_rec%compress_exp
        call mean_rec%sampl_dens_correct
        call mean_rec%ifft
        gridcorr_img = prep3D_inv_kbenvelope4mul([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        call mean_rec%mul(gridcorr_img)
        call gridcorr_img%kill
        fname = 'flex_pca_mean'//MRC_EXT
        call mean_rec%write(fname, del_if_exists=.true.)
        call fname%kill
        ! back to the projectable expanded-Fourier state
        call mean_rec%fft
        call mean_rec%expand_exp
        write(logfhandle,'(A,I0,A,F8.1)') '>>> FLEX_PCA mean estimated from ',used, &
            &' particles, seconds=',toc(t_phase)
        call flush(logfhandle)
    end subroutine estimate_mean_from_data

    !> Worker-side mean scaling: apply the radial scale the MASTER fitted, rather than re-fitting it
    !! from this part's particles (which would use a different stride and hence a different subset).
    subroutine apply_cached_mean_scale( params, mean_rec )
        class(parameters),   intent(inout) :: params
        type(reconstructor), intent(inout) :: mean_rec
        real, allocatable :: filt(:)
        integer :: nyq
        logical :: ok
        nyq = max(1, fdim(params%box_crop) - 1)
        allocate(filt(nyq))
        call read_mean_scale(nyq, filt, ok)
        if( .not. ok ) THROW_HARD('flex_pca worker found no flex_pca_mean_scale.bin from the master')
        call mean_rec%apply_filter(filt)
        call mean_rec%expand_exp
        deallocate(filt)
    end subroutine apply_cached_mean_scale

    !> Self-estimate the amplitude scale of the consensus mean map relative to the whitened data, which
    !! carry SIMPLE's non-unitary gridding convention. A smoothed, clamped per-shell scale is applied to
    !! the mean so that y - T*mu is a residual rather than a difference of two amplitude conventions.
    subroutine estimate_mean_scale( params, build, mean_rec, pinds, nptcls )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
        integer, parameter :: NSAMPLE = 4000
        type(fplane_type), allocatable :: fpls(:)
        type(fplane_type), allocatable :: mean_fpl_t(:)
        type(ori),         allocatable :: ori_t(:)
        integer  :: batchlims(2), batchsz, ibatch, i, iptcl, stride, used, nyq, sh, j, nsub
        integer  :: nthr_here, ithr, t
        integer, allocatable :: sub_pinds(:), used_t(:)
        real(dp) :: s_my, s_mm, s, sm
        real(dp), allocatable :: smy_sh(:), smm_sh(:), sprof(:)
        real(dp), allocatable :: s_my_t(:), s_mm_t(:), smy_sh_t(:,:), smm_sh_t(:,:)
        real,     allocatable :: filt(:)
        logical  :: l_devprep
        nyq = max(1, fdim(params%box_crop) - 1)
        allocate(smy_sh(0:nyq), smm_sh(0:nyq), source=0.d0)
        stride = max(1, nptcls / NSAMPLE)
        s_my = 0.d0; s_mm = 0.d0; used = 0
        ! THREADED OVER PARTICLES: per-thread partial sums, folded in fixed thread order below.
        ! Reproducible at a given nthr; changing nthr moves the fitted scale only at rounding level.
        !$ call omp_set_max_active_levels(1)
        nthr_here = max(1, params%nthr)
        allocate(mean_fpl_t(nthr_here), ori_t(nthr_here), used_t(nthr_here))
        allocate(s_my_t(nthr_here), s_mm_t(nthr_here), source=0.d0)
        allocate(smy_sh_t(0:nyq,nthr_here), smm_sh_t(0:nyq,nthr_here), source=0.d0)
        used_t = 0
        call mean_rec%expand_exp
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        call cov_dev_prep_start(params, build, l_devprep)
        ! Select the strided sample UP FRONT, not inside the batch loop -- otherwise every particle is
        ! read, normalised, padded, FFT'd and CTF-evaluated before ~(1 - 1/stride) of that is discarded.
        ! Same particles in the same order as the serial code; only the summation grouping is per-thread.
        nsub = 0
        do j = 1, nptcls, stride
            nsub = nsub + 1
        end do
        allocate(sub_pinds(nsub))
        nsub = 0
        do j = 1, nptcls, stride
            nsub = nsub + 1
            sub_pinds(nsub) = pinds(j)
        end do
        do ibatch = 1, nsub, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nsub, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call discrete_read_imgbatch(params, build, nsub, sub_pinds, batchlims)
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                &sub_pinds(batchlims(1):batchlims(2)), fpls(:batchsz), mskrad=cov_image_mask_radius(params))
            !$omp parallel do default(shared) private(i,iptcl,ithr) schedule(static) proc_bind(close)
            do i = 1, batchsz
                ithr  = omp_get_thread_num() + 1
                iptcl = sub_pinds(batchlims(1)+i-1)
                call build%spproj_field%get_ori(iptcl, ori_t(ithr))
                if( ori_t(ithr)%isstatezero() ) cycle
                call project_fplane_mean(mean_rec, ori_t(ithr), fpls(i), mean_fpl_t(ithr), &
                    &apply_ctf_amp=.true.)
                s_my_t(ithr) = s_my_t(ithr) + real(cov_herm_inner(mean_fpl_t(ithr), fpls(i)), dp)
                s_mm_t(ithr) = s_mm_t(ithr) + real(cov_herm_inner(mean_fpl_t(ithr), mean_fpl_t(ithr)), dp)
                call plane_shell_cross_accum(mean_fpl_t(ithr), fpls(i), nyq, smy_sh_t(:,ithr), smm_sh_t(:,ithr))
                used_t(ithr) = used_t(ithr) + 1
            end do
            !$omp end parallel do
        end do
        call cov_dev_prep_stop(l_devprep)
        deallocate(sub_pinds)
        do t = 1, nthr_here
            s_my   = s_my + s_my_t(t)
            s_mm   = s_mm + s_mm_t(t)
            smy_sh = smy_sh + smy_sh_t(:,t)
            smm_sh = smm_sh + smm_sh_t(:,t)
            used   = used + used_t(t)
            call ori_t(t)%kill
            call cleanup_plane(mean_fpl_t(t))
        end do
        deallocate(mean_fpl_t, ori_t, used_t, s_my_t, s_mm_t, smy_sh_t, smm_sh_t)
        call cleanup_rec_buffers(build, fpls)
        if( s_mm > DTINY )then
            s = s_my / s_mm
        else
            s = 1.d0
        endif
        if( s <= 0.d0 ) s = 1.d0
        write(logfhandle,'(A,ES12.4,A,I0,A)') '>>> FLEX_PCA mean amplitude self-scale s=',s, &
            &' (from ',used,' particles)'
        ! Per-shell mean/data amplitude scale
        allocate(sprof(0:nyq), filt(nyq))
        write(logfhandle,'(A)') '>>> FLEX_PCA per-shell mean scale s(sh) and s(sh)/s_global (D5):'
        do sh = 0, nyq
            if( smm_sh(sh) > DTINY )then
                sprof(sh) = smy_sh(sh) / smm_sh(sh)
            else
                sprof(sh) = s
            endif
            if( sh <= min(nyq,20) ) write(logfhandle,'(A,I3,A,ES11.3,A,F7.3)') '>>>   sh=',sh, &
                &'  s=',sprof(sh),'  ratio=',sprof(sh)/s
        end do
        call flush(logfhandle)
        ! 3-point-smoothed, clamped radial scale applied to the mean (FT state), then re-expand
        do sh = 1, nyq
            if( sh == 1 )then
                sm = 0.5d0*sprof(1) + 0.5d0*sprof(min(2,nyq))
            else if( sh == nyq )then
                sm = 0.5d0*sprof(nyq) + 0.5d0*sprof(nyq-1)
            else
                sm = 0.25d0*sprof(sh-1) + 0.5d0*sprof(sh) + 0.25d0*sprof(sh+1)
            endif
            filt(sh) = real(min(2.d0*s, max(0.5d0*s, sm)))
        end do
        call mean_rec%apply_filter(filt)
        call mean_rec%expand_exp
        if( flex_pca_is_master() ) call write_mean_scale(nyq, filt)
        deallocate(smy_sh, smm_sh, sprof, filt)
    end subroutine estimate_mean_scale
    !>  Accumulate per-shell mean/data cross power Re<T mu, y> and mean auto power |T mu|^2 over the
    !!  native k<=0 half. The per-shell ratio s(sh)=sum my_sh/sum mm_sh is the ML mean amplitude scale
    !!  at each shell; the k=0 double-count cancels in the ratio so no weighting is needed.
    subroutine plane_shell_cross_accum( mean_fpl, fpl, nyq, my_sh, mm_sh )
        type(fplane_type), intent(in)    :: mean_fpl, fpl
        integer,           intent(in)    :: nyq
        real(dp),          intent(inout) :: my_sh(0:), mm_sh(0:)
        integer     :: pf, h, k, hmin, hmax, kmin, kmax, sh
        complex(dp) :: m, y
        pf   = OSMPL_PAD_FAC
        hmin = pf*ceil_div(lbound(fpl%cmplx_plane,1),pf); hmax = pf*floor_div(ubound(fpl%cmplx_plane,1),pf)
        kmin = pf*ceil_div(lbound(fpl%cmplx_plane,2),pf); kmax = min(0, pf*floor_div(ubound(fpl%cmplx_plane,2),pf))
        do k = kmin, kmax, pf
            do h = hmin, hmax, pf
                sh = nint(sqrt(real((h/pf)**2 + (k/pf)**2)))
                if( sh > nyq ) cycle
                m = cmplx(mean_fpl%cmplx_plane(h,k), kind=dp)
                y = cmplx(fpl%cmplx_plane(h,k),      kind=dp)
                my_sh(sh) = my_sh(sh) + real(conjg(m)*y, dp)
                mm_sh(sh) = mm_sh(sh) + real(conjg(m)*m, dp)
            end do
        end do
    end subroutine plane_shell_cross_accum
    subroutine init_mean_reconstructor( params, build, mean_rec )
        class(parameters),  intent(inout) :: params
        type(builder),      intent(inout) :: build
        type(reconstructor),intent(inout) :: mean_rec
        type(image) :: meanvol
        ! alloc_rho() ends with reset(), which zeros the reconstructor's cmat (and, since rmat/cmat
        ! share the in-place FFT buffer, the real map too).
        call mean_rec%read_and_crop(params%vols(1),params%smpd,params%box_crop,params%smpd_crop)
        call mean_rec%alloc_rho(params,build%spproj,expand=.true.)
        call meanvol%read_and_crop(params%vols(1),params%smpd,params%box_crop,params%smpd_crop)
        call mean_rec%set_rmat(meanvol%get_rmat(), .false.)
        call meanvol%kill
        call mean_rec%fft
        call mean_rec%expand_exp
    end subroutine init_mean_reconstructor
    !> Reconstruction-mode plane from a whitened observation-model plane: numerator T*y and density |T|^2.
    subroutine form_reconstruction_plane( fpl, num )
        type(fplane_type), intent(in)    :: fpl
        type(fplane_type), intent(inout) :: num
        integer :: h1, h2, k1, k2
        h1 = lbound(fpl%cmplx_plane,1); h2 = ubound(fpl%cmplx_plane,1)
        k1 = lbound(fpl%cmplx_plane,2); k2 = ubound(fpl%cmplx_plane,2)
        if( .not. allocated(num%cmplx_plane) ) allocate(num%cmplx_plane(h1:h2,k1:k2))
        if( .not. allocated(num%ctfsq_plane) ) allocate(num%ctfsq_plane(h1:h2,k1:k2))
        num%cmplx_plane = conjg(fpl%transfer_plane) * fpl%cmplx_plane
        num%ctfsq_plane = fpl%ctfsq_plane
        num%frlims  = fpl%frlims
        num%nyq     = fpl%nyq
        num%shconst = fpl%shconst
    end subroutine form_reconstruction_plane









    !> Convert each merged complex column C_q into its two real spatial representatives Re(ifft
    !! C_q)=Sigma*cos_q and Im(ifft C_q)=Sigma*sin_q.
    subroutine columns_to_real_representatives( params, work, colvol, ncol, lb, ub, realvols, nreal )
        class(parameters),   intent(inout) :: params
        type(reconstructor), intent(inout) :: work
        complex,             intent(in)    :: colvol(:,:,:,:)
        integer,             intent(in)    :: ncol, lb(3), ub(3)
        type(image), allocatable, intent(out) :: realvols(:)
        integer,                  intent(out) :: nreal
        type(image)  :: gridcorr_img
        complex, allocatable :: vr(:,:,:), vi(:,:,:)
        integer :: s, i1, i2, i3, n1, n2, n3, hn, kn, ln, h, k, l
        real    :: energy
        n1 = ub(1)-lb(1)+1; n2 = ub(2)-lb(2)+1; n3 = ub(3)-lb(3)+1
        allocate(vr(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
        allocate(vi(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
        gridcorr_img = prep3D_inv_kbenvelope4mul([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        allocate(realvols(2*ncol))
        nreal = 0
        do s = 1, ncol
            do i3 = 1, n3
                l = lb(3)+i3-1; ln = -l
                do i2 = 1, n2
                    k = lb(2)+i2-1; kn = -k
                    do i1 = 1, n1
                        h = lb(1)+i1-1; hn = -h
                        if( hn < lb(1) .or. hn > ub(1) .or. kn < lb(2) .or. kn > ub(2) &
                            &.or. ln < lb(3) .or. ln > ub(3) )then
                            vr(h,k,l) = 0.5*colvol(i1,i2,i3,s)
                            vi(h,k,l) = cmplx(0.,-0.5)*colvol(i1,i2,i3,s)
                        else
                            vr(h,k,l) = 0.5*(colvol(i1,i2,i3,s) + conjg(colvol(n1-i1+1,n2-i2+1,n3-i3+1,s)))
                            vi(h,k,l) = cmplx(0.,-0.5)*(colvol(i1,i2,i3,s) - conjg(colvol(n1-i1+1,n2-i2+1,n3-i3+1,s)))
                        endif
                    end do
                end do
            end do
            call realize_hermitian_volume(params, work, vr, gridcorr_img, energy)
            if( energy > 0. )then
                nreal = nreal + 1
                call realvols(nreal)%copy(work)
            endif
            call realize_hermitian_volume(params, work, vi, gridcorr_img, energy)
            if( energy > 0. )then
                nreal = nreal + 1
                call realvols(nreal)%copy(work)
            endif
        end do
        call gridcorr_img%kill
        deallocate(vr, vi)
    end subroutine columns_to_real_representatives

    !>  Load a Hermitian expanded Fourier volume into the work reconstructor, fold to
    !!  compressed storage, inverse-FFT to a real volume, deapodize, low-pass and mask.
    subroutine realize_hermitian_volume( params, work, vherm, gridcorr_img, energy )
        class(parameters),   intent(in)    :: params
        type(reconstructor), intent(inout) :: work
        complex,             intent(in)    :: vherm(:,:,:)
        type(image),         intent(inout) :: gridcorr_img
        real,                intent(out)   :: energy
        real, pointer :: rmat(:,:,:)
        integer       :: ldim_work(3)
        call work%reset
        call work%reset_exp
        work%cmat_exp = vherm
        call work%compress_exp
        ! Band-limit the covariance column to the signal band FIRST, in Fourier space, before any real-space
        ! operation, so out-of-band shells never enter the masking or the Gram products downstream.
        if( params%lp > 2.0*params%smpd_crop + TINY ) call work%bp(0., params%lp)
        call work%ifft
        ! deapodize on the native lattice (same correction as production gridding)
        call work%mul(gridcorr_img)
        if( params%msk_crop > TINY ) call work%mask3D_soft(params%msk_crop, backgr=0.)
        if( work%is_ft() ) call work%ifft
        call work%get_rmat_ptr(rmat)
        ldim_work = work%get_ldim()
        energy = sum(rmat(1:ldim_work(1),1:ldim_work(2),1:ldim_work(3))**2)
    end subroutine realize_hermitian_volume

    !> Orthonormalize the real column representatives into the column subspace Utilde by Gram
    !! eigendecomposition, keeping every direction above a relative energy floor.
    subroutine orthonormalize_representatives( params, build, realvols, nreal, utilde, utilde_real, d_tilde, svals, &
        &nptcls_basis )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(image),         intent(inout) :: realvols(:)
        integer,             intent(in)    :: nreal
        type(reconstructor), allocatable, intent(out) :: utilde(:)
        type(image),         allocatable, intent(out) :: utilde_real(:)
        integer,             intent(out)   :: d_tilde
        !> squared singular values of the representative set.
        real(dp), allocatable, optional, intent(out) :: svals(:)
        !> particles that actually reached the basis stages. Only used to REPORT the
        !! samples-per-parameter rank bound; it selects nothing.
        integer, optional, intent(in) :: nptcls_basis
        real(dp), allocatable :: gram(:,:), evec(:,:), eval(:)
        real, pointer :: rmat_i(:,:,:), rmat_j(:,:,:)
        integer :: i, q, nrot, keep, d_budget, d_cap, d_signal, d_samples, nbasis
        real(dp) :: lam_max, nrm
        logical  :: l_packed
        character(len=9) :: accum_model
        if( nreal < 1 ) THROW_HARD('flex_pca produced no covariance column representatives')
        allocate(gram(nreal,nreal), evec(nreal,nreal), eval(nreal))
        do i = 1, nreal
            call realvols(i)%get_rmat_ptr(rmat_i)
            do q = i, nreal
                call realvols(q)%get_rmat_ptr(rmat_j)
                gram(i,q) = sum(real(rmat_i,dp)*real(rmat_j,dp))
                gram(q,i) = gram(i,q)
            end do
        end do
        call jacobi(gram, nreal, nreal, eval, evec, nrot)
        call eigsrt(eval, evec, nreal, nreal)              ! descending
        lam_max = max(eval(1), DTINY)
        keep = 0
        do q = 1, nreal
            if( eval(q) > COV_EIG_REL_FLOOR*lam_max ) keep = keep + 1
        end do
        ! Largest d the reduced solve's accumulator can afford: one shared d^2 x d^2 array -> 8*d^4 bytes,
        ! against 8*[d(d+1)/2]^2 ~ 2*d^4 bytes for the packed CG path, which never forms the operator.
        ! Size it against the model the solve will ACTUALLY use -- sizing d for the dense accumulator and
        ! then solving packed spends a quarter of the budget and caps the column subspace 41 % below what
        ! the data supports (at 8 GB: d 177 instead of 250).
        l_packed = cov_packed_cgsolve()
        d_budget = cov_dim_budget(l_packed)
        ! data-driven rank, REPORT ONLY: the energy floor and the memory budget never ask how many
        ! directions are real, so log what the data would support and let the discrepancy show
        d_signal = cov_signal_rank(eval, nreal)
        nbasis   = 0
        if( present(nptcls_basis) ) nbasis = nptcls_basis
        d_samples = 0
        if( nbasis > 0 ) d_samples = &
            &max(1, int((-1.d0 + sqrt(1.d0 + 8.d0*real(nbasis,dp)/COV_SAMPLES_PER_PARAM))/2.d0))
        ! memory budget is a GUARD, not the chooser. SIMPLE_COV_DTILDE replaces both, so a fixed-d A/B
        ! is never silently clamped by the box's RAM.
        d_cap = min(d_budget, COV_DEFAULT_DTILDE)
        call cov_env_int('SIMPLE_COV_DTILDE', d_cap)
        d_tilde  = max(1, min(keep, COV_MAX_DTILDE, d_cap))
        if( cov_env_int_on('SIMPLE_COV_DSIGNAL') ) d_tilde = max(1, min(d_tilde, d_signal))
        if( d_samples > 0 )then
            write(logfhandle,'(A,I0,A,I0,A,F4.1,A)') '>>> FLEX_PCA d_samples=',d_samples, &
                &'  (samples-per-parameter bound from N=',nbasis,' at R=',COV_SAMPLES_PER_PARAM, &
                &') -- REPORT ONLY'
        endif
        write(logfhandle,'(A,I0,A,A,A)') '>>> FLEX_PCA d_signal=',d_signal, &
            &' (spectrum noise-bulk estimate; ', &
            &trim(merge('BINDING ','reported',cov_env_int_on('SIMPLE_COV_DSIGNAL'))), &
            &' -- set SIMPLE_COV_DSIGNAL=1 to bind it)'
        if( l_packed )then
            accum_model = 'packed+CG'
        else
            accum_model = 'dense'
        endif
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA d_tilde=',d_tilde, &
            &'  (above energy floor=',keep,', memory cap=',d_budget,', rank cap=',COV_MAX_DTILDE, &
            &', default=',COV_DEFAULT_DTILDE,')'
        write(logfhandle,'(A,A,A,F8.3,A,F6.3,A)') '>>> FLEX_PCA reduced-solve accumulator model: ', &
            &trim(accum_model),', ',cov_accum_bytes(d_tilde, l_packed)/1.d9, &
            &' GB at this d_tilde (budget ',COV_ATHR_BUDGET/1.d9,' GB)'
        if( d_tilde == d_budget .and. keep > d_budget )then
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA NOTE: the column subspace is limited by the &
                &reduced-solve memory budget, not by the data; ',keep,' directions cleared the energy floor.'
        else if( d_tilde == COV_DEFAULT_DTILDE .and. d_budget > COV_DEFAULT_DTILDE )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA NOTE: d_tilde is the measured default; the &
                &memory budget would have allowed ',d_budget,'. Override with SIMPLE_COV_DTILDE=',d_budget, &
                &' to restore the budget-chosen rank.'
        endif
        call flush(logfhandle)
        allocate(utilde(d_tilde), utilde_real(d_tilde))
        do q = 1, d_tilde
            nrm = sqrt(max(eval(q), DTINY))
            ! Build the unit-norm orthonormal basis vector as a CLEAN plain image via image arithmetic on
            ! the (verified band-limited) representatives, rather than raw rmat-pointer math on padded
            ! reconstructor buffers.
            call utilde_real(q)%copy(realvols(1))
            call utilde_real(q)%zero_and_unflag_ft
            do i = 1, nreal
                call utilde_real(q)%add(realvols(i), real(evec(i,q)/nrm))
            end do
            ! set_rmat(...,.false.) then fft then expand_exp -- NEVER add(), which leaves the reconstructor
            ! flagged as Fourier and silently propagates an untransformed grid (see form_eigenbasis_from_reduced)
            call init_column_reconstructor(params, build, utilde(q))
            call utilde(q)%set_rmat(utilde_real(q)%get_rmat(), .false.)
            call utilde(q)%fft
            call utilde(q)%expand_exp
        end do
        if( present(svals) )then
            allocate(svals(d_tilde))
            do q = 1, d_tilde
                svals(q) = max(eval(q), DTINY)
            end do
        endif
        deallocate(gram, evec, eval)
    end subroutine orthonormalize_representatives



    !> Is the polar (shared-direction) former requested for the reduced solve?
    logical function cov_polar_enabled()
        integer :: v
        v = 0
        call cov_env_int('SIMPLE_COV_POLAR', v)
        cov_polar_enabled = v > 0
    end function cov_polar_enabled

    !> ... and for the embedding? Defaults to following SIMPLE_COV_POLAR, but is separable so the two
    !! stages can be A/B'd one at a time -- changing both at once makes a regression unattributable.
    logical function cov_polar_embed_enabled()
        integer :: v
        v = 0
        call cov_env_int('SIMPLE_COV_POLAR_EMBED', v)
        if( cov_env_is_set('SIMPLE_COV_POLAR_EMBED') )then
            cov_polar_embed_enabled = v > 0
        else
            cov_polar_embed_enabled = cov_polar_enabled()
        endif
    end function cov_polar_embed_enabled

    !> Number of bank directions.
    !!
    !! MEASURED (IgG-RL 10k, box_crop 64, lp 15, d_tilde 64): the largest reduced eigenvalue is
    !! 1194.2 / 1195.5 / 1196.2 / 1198.2 at ndir = 1000 / 2000 / 8000 / 32000 and 1198.3 with the
    !! grid removed entirely (SIMPLE_COV_POLAR_EXACT=1), and ground-truth basis capture is 0.5828
    !! at ndir=2000 against 0.5819 exact and 0.5848 Cartesian. A 6 degree direction grid is
    !! indistinguishable from no discretisation at all here, because this stage never uses a
    !! per-particle b on its own -- it accumulates Sbb and sum_i G_i (x) G_i over 10^4-10^5
    !! particles, and the Gram is additionally a sum over ~10^3 plane samples, so direction error
    !! enters suppressed by 1/sqrt(nsamp) rather than as a per-sample decorrelation.
    !!
    !! So the default targets AMORTISATION (~40 particles per direction) rather than resolution,
    !! with a floor so small datasets still get a reasonable grid. Raise it with
    !! SIMPLE_COV_POLAR_NDIR if a dataset ever shows direction sensitivity -- the bank is streamed
    !! direction by direction, so ndir costs no memory, only bank-build time.
    integer function cov_polar_ndir( nptcls )
        integer, intent(in) :: nptcls
        integer :: v
        v = min(4000, max(1000, nptcls/40))
        v = 2*((v+1)/2)                                  ! build_refspiral needs an even count
        call cov_env_int('SIMPLE_COV_POLAR_NDIR', v)
        v = 2*((v+1)/2)
        cov_polar_ndir = v
    end function cov_polar_ndir




    !> ---------------------------------------------------------------------------------------------
    !> POLAR E-STEP FOR THE EMBEDDING
    !>
    !> Fills exactly the sufficient statistics the Cartesian batch loop of
    !> embed_latents_with_contrast produces -- Gcache, bcache, ccache, contrast, the two residual
    !> energies and the split-half solves zhalf -- so the reliability prior, the re-solve, the
    !> precision matrices and every downstream consumer are untouched.
    !>
    !> Same two-pass shape as the solve's former. The one thing the solve did not need is the
    !> SPLIT-HALF statistics `rho` is built from. Those come out of the grid's ring half-split
    !> (see polar_grid_t%rmid): each ring stores its even angles then its odd ones, so both halves
    !> are contiguous and each half's Gram and b are still a single BLAS call.
    !> ---------------------------------------------------------------------------------------------
    subroutine embed_accumulate_polar( params, build, mean_rec, basis_recs, ncomp, eigvals, sig2, &
        &pinds, nptcls, Gcache, bcache, ccache, zhalf, contrast, resid_energy, resid_mean_energy, &
        &prior, nvalid )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        integer,             intent(in)    :: ncomp, pinds(:), nptcls
        real(dp),            intent(in)    :: eigvals(ncomp), sig2, prior(ncomp)
        real(dp),            intent(inout) :: Gcache(ncomp,ncomp,nptcls), bcache(ncomp,nptcls)
        real(dp),            intent(inout) :: ccache(ncomp,nptcls), zhalf(nptcls,ncomp,2)
        real(dp),            intent(inout) :: contrast(nptcls), resid_energy(nptcls)
        real(dp),            intent(inout) :: resid_mean_energy(nptcls)
        integer,             intent(out)   :: nvalid
        type(polar_grid_t)             :: pg
        type(oris)                     :: dirs
        type(ori)                      :: o
        type(fplane_type), allocatable :: fpls(:)
        real,    allocatable :: rmat_p(:,:,:), nrm_p(:,:), nrm_b(:,:), rmat_b(:,:,:), cav(:), sav(:)
        real,    allocatable :: eul_b(:,:)
        logical :: l_exact_embed
        integer :: j_exact
        integer, allocatable :: dir_of(:), dcnt(:), dptr(:), dlist(:)
        logical, allocatable :: l_zero(:)
        real,    allocatable :: xws(:,:), wrs(:,:), xws1(:,:), xws2(:,:)
        real(dp),allocatable :: eyy(:)
        integer, allocatable :: dirs_in_chunk(:), nvalid_thr(:)
        ! per-thread pass-B workspace; index 0 of the half dimension is the FULL ring
        complex, allocatable :: Ubank(:,:,:)
        real,    allocatable :: Us(:,:,:), Xs(:,:,:), Xs1(:,:,:), Xs2(:,:,:), Reb(:,:,:,:), Csp(:,:,:)
        real(dp),allocatable :: Cf(:,:,:,:), Wmat(:,:,:), Gall(:,:,:,:), Corr(:,:,:,:), den_v(:,:)
        real(dp),allocatable :: c00(:,:), Cm0(:,:,:,:)
        real(dp),allocatable :: Ath(:,:,:), zth(:,:)
        ! pose-refinement workspace: the bank for EVERY direction, kept while the images stream past
        logical :: l_pose
        integer :: ipose_test, i_tmp
        real    :: pose_rot_amp, pose_sh_amp, pose_rot_step, pose_sh_step
        real(dp):: bankgb
        real,    allocatable :: Usall(:,:,:), pose_rot(:), pose_shx(:), pose_shy(:)
        integer, allocatable :: dnn(:,:), ndmove_thr(:)
        real(dp),allocatable :: dangle_thr(:), dq0_thr(:), dq1_thr(:)
        logical :: l_p2
        integer :: nnn_dir, jdacc
        real(dp),allocatable :: Cfall(:,:,:), Cm0all(:,:,:), c00all(:,:)
        real(dp),allocatable :: pose_e0(:), pose_e1(:)
        ! pose band guard
        real     :: sh_cap_A, sh_cap_px, rot_cap_d, smpd_c, mskrad_A, shmag
        real(dp) :: acc_sh_A, acc_rot_A
        integer  :: icap, nclamp_sh, nclamp_rot, npose_w
        real(dp), allocatable :: pose_s0(:), pose_s1(:), pose_ea(:), pose_sa(:)
        integer, allocatable :: pose_n(:), pose_t(:), pose_r(:)
        logical :: l_guard
        integer :: nthr, ithr, i, j, q, r, ir, ibatch, batchlims(2), batchsz, row
        integer :: ndir, nk, nsamp, nsamp2, kto, kfrom, knfrom, knto, pf, nc1
        integer :: hlo, hhi, klo, ph0, pk0, nyq_band, nyq_rec, mmax, m, id, ic, ndchunk
        integer :: jc, idir, ip, ii, i0, nsl, ih, nrb, ncc
        real    :: rmb(3,3), ca, sa, tazim, e3p, shp(2)
        type(string) :: fn_pose
        real(dp):: pw1, cnt1, cc, aa, e_yy, e_mm, myv, res
        integer(timer_int_kind) :: t_a, t_b, t_es
        real :: sec_eread, sec_eprep, sec_esamp
        ! device polar sampler (SIMPLE_COV_GPU=1); pose refinement forces the CPU path
        logical  :: l_gpu_eps, l_eps_res, l_devprep
        integer  :: vps_e
        real,    allocatable :: pwv_e(:), tazv_e(:)
        logical, allocatable :: vmask_e(:)
        real,   parameter :: A_LO_C = 0.1, A_HI_C = 5.0
        nthr   = nthr_glob
        nvalid = 0
        pf     = OSMPL_PAD_FAC
        call mean_rec%expand_exp
        do q = 1, ncomp
            call basis_recs(q)%expand_exp
        end do
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        call cov_dev_prep_start(params, build, l_devprep)
        batchlims = [1, min(nptcls, MAXIMGBATCHSZ)]
        call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
        call prep_imgs4projected_model(params, build, batchlims(2), build%imgbatch(:batchlims(2)), &
            &pinds(1:batchlims(2)), fpls(:batchlims(2)), mskrad=cov_image_mask_radius(params))
        ph0 = lbound(fpls(1)%cmplx_plane,1); pk0 = lbound(fpls(1)%cmplx_plane,2)
        hlo = ceil_div (lbound(fpls(1)%cmplx_plane,1), pf)
        hhi = floor_div(ubound(fpls(1)%cmplx_plane,1), pf)
        klo = ceil_div (lbound(fpls(1)%cmplx_plane,2), pf)
        nyq_rec  = mean_rec%get_lfny(1)
        nyq_band = nyq_rec
        if( fpls(1)%nyq > 0 ) nyq_band = min(nyq_band, max(1, fpls(1)%nyq / pf))
        kfrom  = 1
        kto    = nyq_band
        knfrom = max(kto+1, nint(0.7*real(nyq_rec)))
        knto   = max(knfrom, nyq_rec - 2)
        call polar_grid_build(pg, kfrom, kto, knfrom, knto, hlo, hhi, klo, ph0, pk0)
        nsamp = pg%nsamp; nsamp2 = 2*nsamp; nk = pg%nk
        ! EXACT-DIRECTION MODE. SIMPLE_COV_POLAR_EXACT was only ever honoured by
        ! reduced_solve_accumulate_polar; setting it on this former was bit-identical, i.e. a dead
        ! knob. It is the control that separates the two things the banked path confounds -- the
        ! polar QUADRATURE from the direction QUANTISATION -- so it is worth having here. Every
        ! particle becomes its own bank direction, which costs the whole shared-direction
        ! amortisation (C_qr(k) then serves one particle instead of many) and is therefore a
        ! diagnostic, not an operating point.
        j_exact = 0
        call cov_env_int('SIMPLE_COV_POLAR_EXACT', j_exact)
        l_exact_embed = j_exact > 0
        ndir  = cov_polar_ndir(nptcls)
        if( l_exact_embed ) ndir = nptcls
        ! ---------------- orientations and direction assignment
        allocate(rmat_p(3,3,nptcls), nrm_p(3,nptcls), dir_of(nptcls), cav(nptcls), sav(nptcls))
        allocate(l_zero(nptcls), source=.false.)
        do i = 1, nptcls
            call build%spproj_field%get_ori(pinds(i), o)
            rmat_p(:,:,i) = o%get_mat()
            nrm_p(:,i)    = rmat_p(3,:,i)
            l_zero(i)     = o%isstatezero()
        end do
        call o%kill
        allocate(nrm_b(3,ndir), rmat_b(3,3,ndir), eul_b(3,ndir))
        if( l_exact_embed )then
            ! each particle is its own direction: the relative in-plane angle is the identity, so
            ! the polar sampler reads the particle plane at its own orientation exactly
            rmat_b = rmat_p
            do i = 1, nptcls
                nrm_b(:,i) = rmat_b(3,:,i)
                call build%spproj_field%get_ori(pinds(i), o)
                eul_b(:,i) = o%get_euler()
                dir_of(i)  = i
                cav(i)     = 1.0
                sav(i)     = 0.0
                if( l_zero(i) ) dir_of(i) = 0
            end do
            call o%kill
        else
        call dirs%new(ndir, is_ptcl=.false.)
        call build%pgrpsyms%build_refspiral(dirs)
        do j = 1, ndir
            rmat_b(:,:,j) = dirs%get_mat(j)
            nrm_b(:,j)    = rmat_b(3,:,j)
            eul_b(:,j)    = dirs%get_euler(j)
        end do
        call dirs%kill
        call polar_assign_directions(nrm_p, nptcls, nrm_b, ndir, dir_of)
        do i = 1, nptcls
            call polar_relative_inplane(rmat_p(:,:,i), rmat_b(:,:,dir_of(i)), ca, sa)
            cav(i) = ca; sav(i) = sa
            if( l_zero(i) ) dir_of(i) = 0
        end do
        endif
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA POLAR EMBED: samples=',nsamp, &
            &'  rings=',nk,'  directions=',ndir
        call flush(logfhandle)
        ! bank + ring-Gram scratch: allocated here because the pose block below needs them
        ! before pass A, and pass B reuses the same buffers
        allocate(Ubank(nsamp,0:ncomp,nthr), Us(nsamp2,0:ncomp,nthr))
        allocate(Csp(0:ncomp,0:ncomp,nthr))
        ! ---------------- POSE REFINEMENT SETUP (P1: in-plane angle + shift)
        !
        ! Pose search needs the model in hand while the IMAGE is in hand, so when it is on the bank
        ! is built for every direction up front rather than direction-major in pass B. At embedding
        ! rank (ncomp ~ 16-20, not d_tilde ~ 128-250) that bank is a few hundred MB, which is exactly
        ! why this is affordable here and was not in the solve.
        l_pose = cov_pose_mode() > 0
        ipose_test = 0
        call cov_env_int('SIMPLE_COV_POSE_TEST', ipose_test)
        l_p2 = .false.
        if( l_pose ) l_p2 = cov_pose_mode() >= 2
        ! injected-perturbation amplitudes and search steps, all in 0.1 units so they can be swept
        ! from the environment without a rebuild
        i_tmp = ipose_test; call cov_env_int('SIMPLE_COV_POSE_TEST_ROT', i_tmp)
        pose_rot_amp  = 0.1*real(i_tmp)
        i_tmp = ipose_test; call cov_env_int('SIMPLE_COV_POSE_TEST_SH',  i_tmp)
        pose_sh_amp   = 0.1*real(i_tmp)
        i_tmp = 20; call cov_env_int('SIMPLE_COV_POSE_ROTSTEP', i_tmp)
        pose_rot_step = 0.1*real(i_tmp)
        i_tmp = 10; call cov_env_int('SIMPLE_COV_POSE_SHSTEP',  i_tmp)
        pose_sh_step  = 0.1*real(i_tmp)
        ! The split-half acceptance test is ON for P1 and OFF for P2, because that is what the
        ! measurements say. P1's signal (dominated by shift) is strong enough that the guard rejects
        ! mostly noise: it halves the damage to correct shifts at no cost to recovery. P2's direction
        ! signal is weak, and the same guard throws away half the REAL improvements -- direction
        ! recovery 3.11 -> 2.77 deg with it against 3.11 -> 2.25 without, and the refit ends up worse
        ! than not refining at all (capture 0.546 vs 0.573). SIMPLE_COV_POSE_GUARD overrides.
        l_guard = .not. l_p2
        if( cov_env_is_set('SIMPLE_COV_POSE_GUARD') ) &
            &l_guard = .not. cov_env_int_off('SIMPLE_COV_POSE_GUARD')
        if( l_pose )then
            bankgb = 4.d0*real(nsamp2,dp)*real(ncomp+1,dp)*real(ndir,dp)/1.d9
            write(logfhandle,'(A,F7.2,A)') '>>> FLEX_PCA POLAR POSE: full bank ',bankgb,' GB'
            if( bankgb > 6.d0 )then
                write(logfhandle,'(A)') '>>> FLEX_PCA POLAR POSE DISABLED: bank too large; &
                    &lower SIMPLE_COV_POLAR_NDIR'
                l_pose = .false.
            endif
        endif
        nnn_dir = 1
        if( l_p2 )then
            nnn_dir = 12
            call cov_env_int('SIMPLE_COV_POSE_NNN', nnn_dir)
            nnn_dir = max(2, min(nnn_dir, ndir-1))
        endif
        allocate(dnn(nnn_dir,ndir))
        if( l_pose )then
            allocate(Usall(nsamp2,0:ncomp,ndir), Cfall(ncomp*ncomp,nk,ndir))
            allocate(Cm0all(ncomp,nk,ndir), c00all(nk,ndir))
            if( l_p2 )then
                call polar_dir_neighbours(ndir, nnn_dir, nrm_b, dnn)
                write(logfhandle,'(A,I0,A,F6.2,A)') '>>> FLEX_PCA POLAR POSE P2: direction search &
                    &over ',nnn_dir,' neighbours (grid spacing ~', &
                    &2.0*180.0/(PI*sqrt(real(ndir))),' deg)'
            else
                do id = 1, ndir
                    dnn(:,id) = id
                end do
            endif
            !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
            !$omp& private(id,ithr,rmb,q,j,ir)
            do id = 1, ndir
                ithr = omp_get_thread_num() + 1
                rmb  = rmat_b(:,:,id)
                call polar_project_recs(mean_rec, basis_recs, ncomp, rmb, pg, Ubank(:,:,ithr))
                do q = 0, ncomp
                    do j = 1, nsamp
                        Usall(2*j-1,q,id) = pg%sqwq(j)*real(Ubank(j,q,ithr))
                        Usall(2*j,  q,id) = pg%sqwq(j)*aimag(Ubank(j,q,ithr))
                    end do
                end do
                do ir = 1, nk
                    call polar_ring_gram(Usall(1,0,id), nsamp2, ncomp, pg%rbeg(ir), &
                        &pg%rend(ir)-pg%rbeg(ir)+1, Csp(0,0,ithr), Cfall(1,ir,id), Cm0all(1,ir,id))
                    c00all(ir,id) = polar_ring_selfpower(Usall(1,0,id), nsamp2, pg%rbeg(ir), &
                        &pg%rend(ir)-pg%rbeg(ir)+1)
                end do
            end do
            !$omp end parallel do
        endif
        ! ---------------- PASS A
        t_a = tic()
        allocate(xws(nsamp2,nptcls), source=0.0)
        ! the two lattice-parity half-fields. Storing them costs 2 x nsamp2 x nptcls sp, which is the
        ! price of a reliability estimate that is not contaminated by shared interpolation support.
        allocate(xws1(nsamp2,nptcls), xws2(nsamp2,nptcls), source=0.0)
        allocate(wrs(nk,nptcls), source=0.0)
        allocate(eyy(nptcls), source=0.d0)
        allocate(pose_rot(nptcls), pose_shx(nptcls), pose_shy(nptcls), source=0.0)
        allocate(pose_e0(nthr), pose_e1(nthr), pose_s0(nthr), pose_s1(nthr), source=0.d0)
        allocate(pose_ea(nthr), pose_sa(nthr), source=0.d0)
        allocate(pose_n(nthr), pose_t(nthr), pose_r(nthr), ndmove_thr(nthr), source=0)
        allocate(dangle_thr(nthr), dq0_thr(nthr), dq1_thr(nthr), source=0.d0)
        sec_eread = 0.; sec_eprep = 0.; sec_esamp = 0.
        vps_e = 0
        call cov_env_int('SIMPLE_COV_GPU_PSAMPLE', vps_e)
        l_gpu_eps = vps_e > 0 .and. flex_gpu_available() .and. .not. l_pose
        ! resident sampler default (see the solve pass-A note); pose refinement forces CPU
        l_eps_res = (.not. l_gpu_eps) .and. flex_gpu_available() .and. .not. l_pose .and. &
            &.not. cov_env_int_off('SIMPLE_COV_GPU_PSAMPLE') .and. &
            &cov_image_mask_radius(params) <= 0.
        if( l_eps_res .and. params%l_ml_reg ) l_eps_res = allocated(build%esig%sigma2_noise)
        if( l_gpu_eps .or. l_eps_res )then
            call flex_gpu_psample_begin_f(pg%rad, pg%cs, pg%sn, pg%sqwq, pg%rbeg, pg%rend, &
                &pg%nsamp, pg%nk, pg%nrad, pg%ncs, pg%nsn, pg%nwq, pg%nsamp_n)
            allocate(pwv_e(MAXIMGBATCHSZ), tazv_e(MAXIMGBATCHSZ), vmask_e(MAXIMGBATCHSZ))
            if( l_eps_res )then
                call flex_gpu_prep_begin_f(build%lmsk, params%box, params%boxpd, &
                    &MAXIMGBATCHSZ, 0.0, .true.)
                write(logfhandle,'(A)') '>>> FLEX_PCA POLAR EMBED PASS-A RESIDENT device sampler ON'
            else
                write(logfhandle,'(A)') '>>> FLEX_PCA POLAR EMBED PASS-A device sampler ON'
            endif
            call flush(logfhandle)
        endif
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            t_es = tic()
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
            sec_eread = sec_eread + toc(t_es)
            t_es = tic()
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz), &
                &mskrad=cov_image_mask_radius(params), resident=l_eps_res)
            sec_eprep = sec_eprep + toc(t_es)
            t_es = tic()
            if( l_gpu_eps .or. l_eps_res )then
                do i = 1, batchsz
                    vmask_e(i) = dir_of(batchlims(1)+i-1) > 0
                end do
                if( l_eps_res )then
                    call flex_gpu_psample_batch_res_f(cav(batchlims(1):batchlims(2)), &
                        &sav(batchlims(1):batchlims(2)), vmask_e(:batchsz), batchsz, .true., &
                        &xws, xws1, xws2, wrs, batchlims(1), pwv_e(:batchsz), tazv_e(:batchsz))
                else
                    call flex_gpu_psample_batch_f(fpls(:batchsz), cav(batchlims(1):batchlims(2)), &
                        &sav(batchlims(1):batchlims(2)), vmask_e(:batchsz), batchsz, .true., &
                        &xws, xws1, xws2, wrs, batchlims(1), pwv_e(:batchsz), tazv_e(:batchsz))
                endif
                !$omp parallel do default(shared) schedule(static) proc_bind(close) private(i,row)
                do i = 1, batchsz
                    row = batchlims(1) + i - 1
                    if( dir_of(row) <= 0 ) cycle
                    eyy(row) = polar_self_energy(xws(:,row), wrs(:,row), pg)
                end do
                !$omp end parallel do
                sec_esamp = sec_esamp + toc(t_es)
                cycle
            endif
            !$omp parallel do default(shared) schedule(static) proc_bind(close) &
            !$omp& private(i,row,ithr,tazim,pw1,cnt1,jdacc)
            do i = 1, batchsz
                row = batchlims(1) + i - 1
                if( dir_of(row) <= 0 ) cycle
                ithr = omp_get_thread_num() + 1
                call polar_sample_particle_packed(fpls(i), pg, cav(row), sav(row), &
                    &xws(:,row), wrs(:,row), pw1, cnt1, tazim, xws1(:,row), xws2(:,row))
                ! <y,y> in the polar measure; the Cartesian path takes it with cov_herm_inner
                eyy(row) = polar_self_energy(xws(:,row), wrs(:,row), pg)
                if( l_pose )then
                    jdacc = dir_of(row)
                    call polar_pose_refine_one(pg, fpls(i), Usall, Cfall, Cm0all, c00all, &
                        &nsamp2, ncomp, nk, ndir, dir_of(row), row, wrs(:,row), prior, sig2, &
                        &ipose_test, pose_rot_amp, pose_sh_amp, pose_rot_step, pose_sh_step, &
                        &l_guard, l_p2, nnn_dir, dnn(:,dir_of(row)), rmat_p(:,:,row), rmat_b, &
                        &jdacc, cav(row), sav(row), pose_rot(row), pose_shx(row), pose_shy(row), &
                        &xws(:,row), pose_e0(ithr), pose_e1(ithr), pose_s0(ithr), pose_s1(ithr), &
                        &pose_ea(ithr), pose_sa(ithr), pose_n(ithr), pose_t(ithr), pose_r(ithr))
                    ! a direction move changes which bank pass B must use for this particle
                    if( jdacc /= dir_of(row) )then
                        ndmove_thr(ithr) = ndmove_thr(ithr) + 1
                        dangle_thr(ithr) = dangle_thr(ithr) + real(acos(max(-1.,min(1., &
                            &dot_product(rmat_b(3,:,jdacc), rmat_b(3,:,dir_of(row))))))*180./PI, dp)
                    endif
                    ! recovery of a KNOWN direction degradation, if one was injected: the residual is
                    ! measured for BOTH the incoming and the accepted grid direction, so the grid's
                    ! own quantisation floor is present in both and cancels from the comparison
                    if( allocated(COV_TRUTH_NRM) )then
                        dq0_thr(ithr) = dq0_thr(ithr) + real(acos(max(-1.,min(1., &
                            &dot_product(rmat_b(3,:,dir_of(row)), COV_TRUTH_NRM(:,row)))))*180./PI, dp)**2
                        dq1_thr(ithr) = dq1_thr(ithr) + real(acos(max(-1.,min(1., &
                            &dot_product(rmat_b(3,:,jdacc), COV_TRUTH_NRM(:,row)))))*180./PI, dp)**2
                    endif
                    dir_of(row) = jdacc
                endif
            end do
            !$omp end parallel do
            sec_esamp = sec_esamp + toc(t_es)
        end do
        ! ---- POSE BAND GUARD, applied before anything is written back.
        !
        ! The pose block sees ONLY frequencies below `params%lp`, so it cannot justify a correction
        ! larger than that band can resolve. Nothing enforced this, and nothing reported the
        ! correction in physical units, which is how the following got through on EMPIAR-10330
        ! (PfCRT, published at 3.2 A): the block ran at the commander's default `box_crop=64` on a
        ! 300-box particle, i.e. `smpd_crop = 4.85 A/px`, and applied a 1.84 px shift -- which reads
        ! as small and is a **8.9 A** displacement. Consensus resolution went 3.83 A -> 7.06 A.
        !
        ! The same 1.84 px is 2.4 A at `smpd_crop = 1.3`. Pixels hide the failure; Angstroms do not.
        ! Everything below is therefore in Angstroms, and the cap is a fraction of the band: a
        ! displacement of `d` puts a phase error of `2*pi*d/lp` on the band edge, so `lp/4` is a
        ! quarter turn there and is already generous.
        !
        ! Rotation is converted to the displacement it causes at the mask radius, so one budget
        ! covers both degrees of freedom -- which is the right way to compare them
        ! (see the shift-vs-rotation result: they are one identifiability problem in displacement).
        if( l_pose )then
            sh_cap_A = 0.25 * params%lp
            icap = 0
            call cov_env_int('SIMPLE_COV_POSE_MAXSH', icap)      ! tenths of an Angstrom
            if( icap > 0 ) sh_cap_A = 0.1 * real(icap)
            smpd_c   = params%smpd_crop
            mskrad_A = 0.5 * params%msk_crop * smpd_c
            if( mskrad_A <= 0. ) mskrad_A = 0.25 * real(params%box_crop) * smpd_c
            sh_cap_px = sh_cap_A / max(smpd_c, TINY)
            rot_cap_d = sh_cap_A / max(mskrad_A, TINY) * 180. / PI
            nclamp_sh = 0; nclamp_rot = 0
            acc_sh_A  = 0.d0; acc_rot_A = 0.d0; npose_w = 0
            do i = 1, nptcls
                if( dir_of(i) <= 0 ) cycle
                npose_w  = npose_w + 1
                shmag    = sqrt(pose_shx(i)**2 + pose_shy(i)**2)
                acc_sh_A = acc_sh_A + real(shmag*smpd_c, dp)**2
                acc_rot_A= acc_rot_A + real(abs(pose_rot(i))*PI/180.*mskrad_A, dp)**2
                if( shmag > sh_cap_px )then
                    pose_shx(i) = pose_shx(i) * (sh_cap_px/shmag)
                    pose_shy(i) = pose_shy(i) * (sh_cap_px/shmag)
                    nclamp_sh   = nclamp_sh + 1
                endif
                if( abs(pose_rot(i)) > rot_cap_d )then
                    pose_rot(i) = sign(rot_cap_d, pose_rot(i))
                    nclamp_rot  = nclamp_rot + 1
                endif
            end do
            write(logfhandle,'(A,F7.2,A,F6.3,A,F7.2,A,F6.2,A)') &
                &'>>> FLEX_PCA POSE BAND GUARD: lp=',params%lp,' A  smpd_crop=',smpd_c, &
                &' A/px  cap=',sh_cap_A,' A (',sh_cap_px,' px)'
            write(logfhandle,'(A,F7.2,A,F7.2,A)') '>>>   accepted correction rms BEFORE cap: shift=', &
                &sqrt(sum_dp_safe(acc_sh_A, npose_w)),' A   in-plane edge displacement=', &
                &sqrt(sum_dp_safe(acc_rot_A, npose_w)),' A'
            write(logfhandle,'(A,I0,A,I0,A,I0,A,F6.1,A)') '>>>   clamped: shift ',nclamp_sh, &
                &'  rotation ',nclamp_rot,'  of ',npose_w,' particles (', &
                &100.*real(nclamp_sh)/real(max(1,npose_w)),' % of shifts)'
            if( npose_w > 0 .and. real(nclamp_sh)/real(npose_w) > 0.2 )then
                write(logfhandle,'(A,F6.1,A)') '>>> FLEX_PCA POSE WARNING: ', &
                    &100.*real(nclamp_sh)/real(npose_w),' % of shift corrections exceeded the band cap.'
                write(logfhandle,'(A)') '>>>   The pose block is proposing moves its own band cannot &
                    &justify. Either the incoming poses are much better than lp (do NOT refine), or &
                    &box_crop/lp are too coarse for this data. Compare a consensus reconstruction at &
                    &the refined poses against one at the incoming poses BEFORE trusting the result.'
            endif
            call flush(logfhandle)
        endif
        ! ---- WRITE THE REFINED POSES BACK.
        !
        ! The pose block lives in the EMBEDDING, which runs after the mean, the covariance columns
        ! and the basis. So refining here cannot change the delivered eigenvolumes at all -- measured:
        ! GT basis capture was identical to 4 decimal places with refinement on and off. To reach the
        ! basis the refined poses have to go back into the project and the fit has to be repeated;
        ! that outer loop is the joint refinement the design doc's section 4 describes, and this is
        ! the hook for it. Serial on purpose: the pose search runs inside an OpenMP region and the
        ! project field is shared.
        if( l_pose .and. cov_env_is_set('SIMPLE_COV_POSE_WRITE') )then
            do i = 1, nptcls
                if( dir_of(i) <= 0 ) cycle
                if( l_p2 )then
                    ! P2 may have moved the particle to a different bank direction. The accepted
                    ! orientation is that direction's (e1,e2) with the accepted in-plane angle --
                    ! taken from the bank's own Euler triplet rather than by inverting a rotation
                    ! matrix, so no m2euler branch/convention can silently corrupt it.
                    call build%spproj_field%set_euler(pinds(i), &
                        &[eul_b(1,dir_of(i)), eul_b(2,dir_of(i)), &
                        & eul_b(3,dir_of(i)) + polar_inplane_deg(cav(i), sav(i))])
                else
                    e3p = build%spproj_field%e3get(pinds(i))
                    call build%spproj_field%e3set(pinds(i), e3p + pose_rot(i))
                endif
                shp = build%spproj_field%get_2Dshift(pinds(i))
                call build%spproj_field%set_shift(pinds(i), shp + [pose_shx(i), pose_shy(i)])
            end do
            fn_pose = 'flex_pca_refined_poses.simple'
            call build%spproj%write(fn_pose)
            call fn_pose%kill
            write(logfhandle,'(A)') '>>> FLEX_PCA REFINED POSES WRITTEN: flex_pca_refined_poses.simple &
                &-- feed it back as projfile= for the next outer iteration'
            call flush(logfhandle)
        endif
        call cov_dev_prep_stop(l_devprep)
        if( l_eps_res ) call flex_gpu_prep_free_f
        write(logfhandle,'(A,F8.1)') '>>> FLEX_PCA POLAR EMBED PASS-A SECONDS: ',toc(t_a)
        write(logfhandle,'(A,F7.1,A,F7.1,A,F7.1)') '>>> FLEX_PCA POLAR EMBED PASS-A SPLIT (seconds): read=', &
            &sec_eread,'  prep=',sec_eprep,'  sampling=',sec_esamp
        if( l_gpu_eps .or. l_eps_res )then
            call flex_gpu_psample_free_f
            deallocate(pwv_e, tazv_e, vmask_e)
        endif
        if( l_pose .and. sum(pose_n) > 0 )then
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA POSE REFINEMENT (P1: in-plane + shift) on ', &
                &sum(pose_n),' particles'
            if( ipose_test > 0 )then
                write(logfhandle,'(A,F8.4,A,F8.4)') '>>>   PERTURB-RECOVER angle (deg): injected rms=', &
                    &sqrt(sum(pose_e0)/real(sum(pose_n),dp)),'  residual rms=', &
                    &sqrt(sum(pose_e1)/real(sum(pose_n),dp))
                write(logfhandle,'(A,F8.4,A,F8.4)') '>>>   PERTURB-RECOVER shift (px):  injected rms=', &
                    &sqrt(sum(pose_s0)/real(sum(pose_n),dp)),'  residual rms=', &
                    &sqrt(sum(pose_s1)/real(sum(pose_n),dp))
                write(logfhandle,'(A,F7.3,A)') '>>>   the TRUE pose outscores the search answer for ', &
                    &100.d0*real(sum(pose_t),dp)/real(sum(pose_n),dp), &
                    &' % of particles (0 % = search converged, ~50 % = objective uninformative)'
            else if( COV_PROJ_PERTURB_ROT > 0. .or. COV_PROJ_PERTURB_SH > 0. )then
                write(logfhandle,'(A,F8.4,A,F8.4,A,F8.4)') &
                    &'>>>   RECOVERY vs the PROJECT degradation, angle (deg): injected rms=', &
                    &sqrt(sum(pose_e0)/real(sum(pose_n),dp)),'  residual(+)=', &
                    &sqrt(sum(pose_e1)/real(sum(pose_n),dp)),'  residual(-)=', &
                    &sqrt(sum(pose_ea)/real(sum(pose_n),dp))
                write(logfhandle,'(A,F8.4,A,F8.4,A,F8.4)') &
                    &'>>>   RECOVERY vs the PROJECT degradation, shift (px):  injected rms=', &
                    &sqrt(sum(pose_s0)/real(sum(pose_n),dp)),'  residual(+)=', &
                    &sqrt(sum(pose_s1)/real(sum(pose_n),dp)),'  residual(-)=', &
                    &sqrt(sum(pose_sa)/real(sum(pose_n),dp))
            else
                write(logfhandle,'(A,F8.4,A,F8.4)') '>>>   accepted move rms: angle (deg)=', &
                    &sqrt(sum(pose_e1)/real(sum(pose_n),dp)),'  shift (px)=', &
                    &sqrt(sum(pose_s1)/real(sum(pose_n),dp))
                write(logfhandle,'(A,F7.3,A)') '>>>   the score improved for ', &
                    &100.d0*real(sum(pose_t),dp)/real(sum(pose_n),dp),' % of particles'
            endif
            write(logfhandle,'(A,F7.3,A,L2)') '>>>   split-half guard rejected ', &
                &100.d0*real(sum(pose_r),dp)/real(sum(pose_n),dp),' % of proposed moves; guard=',l_guard
            if( l_p2 ) write(logfhandle,'(A,F7.3,A,F7.3,A)') '>>>   direction changed for ', &
                &100.d0*real(sum(ndmove_thr),dp)/real(sum(pose_n),dp),' % of particles, mean move ', &
                &sum(dangle_thr)/max(real(sum(ndmove_thr),dp),1.d0),' deg'
            if( allocated(COV_TRUTH_NRM) ) write(logfhandle,'(A,F8.4,A,F8.4)') &
                &'>>>   DIRECTION recovery vs truth (deg): before rms=', &
                &sqrt(sum(dq0_thr)/real(sum(pose_n),dp)),'  after rms=', &
                &sqrt(sum(dq1_thr)/real(sum(pose_n),dp))
        endif
        call flush(logfhandle)
        call cleanup_rec_buffers(build, fpls)
        deallocate(nrm_p, nrm_b)
        ! ---------------- direction -> particle CSR
        allocate(dcnt(ndir), source=0)
        do i = 1, nptcls
            if( dir_of(i) > 0 ) dcnt(dir_of(i)) = dcnt(dir_of(i)) + 1
        end do
        allocate(dptr(ndir+1)); dptr(1) = 1
        do j = 1, ndir
            dptr(j+1) = dptr(j) + dcnt(j)
        end do
        allocate(dlist(max(1,dptr(ndir+1)-1)))
        dcnt = 0
        do i = 1, nptcls
            j = dir_of(i)
            if( j <= 0 ) cycle
            dlist(dptr(j)+dcnt(j)) = i
            dcnt(j) = dcnt(j) + 1
        end do
        mmax = 0
        do j = 1, ndir
            mmax = max(mmax, dcnt(j))
        end do
        mmax = min(max(64, mmax), 256)
        ! ---------------- PASS B
        t_b  = tic()
        ncc  = ncomp + 1                    ! bank slot 0 is the mean
        ! half index 0 = the full ring, 1 and 2 = the two interleaved half-sets
        allocate(Cf(ncomp*ncomp,nk,0:2,nthr), c00(nk,nthr), Cm0(ncomp,nk,0:2,nthr))
        allocate(Xs(nsamp2,mmax,nthr), Reb(0:ncomp,mmax,0:2,nthr), Wmat(nk,mmax,nthr))
        allocate(Xs1(nsamp2,mmax,nthr), Xs2(nsamp2,mmax,nthr))
        allocate(Gall(ncomp*ncomp,mmax,0:2,nthr), Corr(ncomp,mmax,0:2,nthr), den_v(mmax,nthr))
        allocate(Ath(ncomp,ncomp,nthr), zth(ncomp,nthr))
        allocate(nvalid_thr(nthr), source=0)
        ndchunk = max(4*nthr, 64)
        allocate(dirs_in_chunk(ndchunk))
        idir = 1
        do while( idir <= ndir )
            ic = 0
            do while( idir <= ndir .and. ic < ndchunk )
                if( dcnt(idir) > 0 )then
                    ic = ic + 1
                    dirs_in_chunk(ic) = idir
                endif
                idir = idir + 1
            end do
            if( ic == 0 ) cycle
            !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
            !$omp& private(jc,id,m,ithr,rmb,j,q,r,ir,ip,ii,i0,nsl,ih,nrb,nc1,row) &
            !$omp& private(cc,aa,e_yy,e_mm,myv,res)
            do jc = 1, ic
                id   = dirs_in_chunk(jc)
                m    = dcnt(id)
                ithr = omp_get_thread_num() + 1
                rmb  = rmat_b(:,:,id)
                call polar_project_recs(mean_rec, basis_recs, ncomp, rmb, pg, Ubank(:,:,ithr))
                do q = 0, ncomp
                    do j = 1, nsamp
                        Us(2*j-1,q,ithr) = pg%sqwq(j)*real(Ubank(j,q,ithr))
                        Us(2*j,  q,ithr) = pg%sqwq(j)*aimag(Ubank(j,q,ithr))
                    end do
                end do
                ! Ring Gram tables. The half-set G is HALF the full one, not a subset sum over
                ! polar samples: the split is now by Cartesian lattice parity (see the b-halves
                ! below), each parity carries half the lattice measure, and the basis is noise-free
                ! and smooth so both parities see the same signal. That keeps A_half consistent with
                ! b_half, which is what the reliability solve needs.
                do ir = 1, nk
                    call polar_ring_gram(Us(1,0,ithr), nsamp2, ncomp, pg%rbeg(ir), &
                        &pg%rend(ir) - pg%rbeg(ir) + 1, &
                        &Csp(0,0,ithr), Cf(1,ir,0,ithr), Cm0(1,ir,0,ithr))
                    Cf(:,ir,1,ithr)  = 0.5d0*Cf(:,ir,0,ithr)
                    Cf(:,ir,2,ithr)  = Cf(:,ir,1,ithr)
                    Cm0(:,ir,1,ithr) = 0.5d0*Cm0(:,ir,0,ithr)
                    Cm0(:,ir,2,ithr) = Cm0(:,ir,1,ithr)
                    ! <T mu, T mu> per ring, from the mean's own diagonal entry
                    c00(ir,ithr) = polar_ring_selfpower(Us(1,0,ithr), nsamp2, pg%rbeg(ir), &
                        &pg%rend(ir) - pg%rbeg(ir) + 1)
                end do
                do i0 = 1, m, mmax
                    nsl = min(mmax, m - i0 + 1)
                    do ii = 1, nsl
                        ip = dlist(dptr(id)+i0+ii-2)
                        Xs(:,ii,ithr)   = xws(:,ip)
                        Xs1(:,ii,ithr)  = xws1(:,ip)
                        Xs2(:,ii,ithr)  = xws2(:,ip)
                        Wmat(:,ii,ithr) = real(wrs(:,ip), dp)
                    end do
                    ! b for the full plane and for each half. The halves are not contiguous ACROSS
                    ! rings, so each is accumulated ring by ring; the full one is one GEMM.
                    call sgemm('T','N', ncomp+1, nsl, nsamp2, 1.0, Us(1,0,ithr), nsamp2, &
                        &Xs(1,1,ithr), nsamp2, 0.0, Reb(0,1,0,ithr), ncomp+1)
                    ! Half-b from the two LATTICE-PARITY fields, over the whole sample list. The old
                    ! ring-interleaved split alternated polar samples one lattice unit apart, which
                    ! the 3-tap kernel makes share most of their support: measured, that pushed every
                    ! rho to 0.95-1.00 on EMPIAR-10076 and flattened the reliability prior.
                    call sgemm('T','N', ncomp+1, nsl, nsamp2, 1.0, Us(1,0,ithr), nsamp2, &
                        &Xs1(1,1,ithr), nsamp2, 0.0, Reb(0,1,1,ithr), ncomp+1)
                    call sgemm('T','N', ncomp+1, nsl, nsamp2, 1.0, Us(1,0,ithr), nsamp2, &
                        &Xs2(1,1,ithr), nsamp2, 0.0, Reb(0,1,2,ithr), ncomp+1)
                    ! G and the mean cross term, full and per half, all particles at once
                    do ih = 0, 2
                        call dgemm('N','N', ncomp*ncomp, nsl, nk, 1.d0, Cf(1,1,ih,ithr), ncomp*ncomp, &
                            &Wmat(1,1,ithr), nk, 0.d0, Gall(1,1,ih,ithr), ncomp*ncomp)
                        call dgemm('N','N', ncomp, nsl, nk, 1.d0, Cm0(1,1,ih,ithr), ncomp, &
                            &Wmat(1,1,ithr), nk, 0.d0, Corr(1,1,ih,ithr), ncomp)
                    end do
                    call dgemv('T', nk, nsl, 1.d0, Wmat(1,1,ithr), nk, c00(1,ithr), 1, 0.d0, &
                        &den_v(1,ithr), 1)
                    do ii = 1, nsl
                        row  = dlist(dptr(id)+i0+ii-2)
                        e_mm = den_v(ii,ithr)
                        myv  = real(Reb(0,ii,0,ithr),dp)
                        e_yy = eyy(row)
                        do q = 1, ncomp
                            bcache(q,row) = real(Reb(q,ii,0,ithr),dp)
                            ccache(q,row) = Corr(q,ii,0,ithr)
                            do r = 1, ncomp
                                Gcache(q,r,row) = Gall((r-1)*ncomp+q,ii,0,ithr)
                            end do
                        end do
                        resid_mean_energy(row) = e_yy - 2.d0*myv + e_mm
                        if( COV_UNIT_CONTRAST .and. .not. cov_fit_contrast_rt )then
                            cc = 1.d0
                        else
                            cc = max(real(A_LO_C,dp), min(real(A_HI_C,dp), myv/max(e_mm,DTINY)))
                        endif
                        contrast(row) = cc
                        aa = cc*cc
                        ! MAP solve at the chosen contrast, identical arithmetic to the Cartesian body
                        Ath(:,:,ithr) = (aa/sig2)*Gcache(:,:,row)
                        do q = 1, ncomp
                            Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                            zth(q,ithr)   = (cc*bcache(q,row) - aa*ccache(q,row))/sig2
                        end do
                        call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), ncomp)
                        res = e_yy + aa*e_mm - 2.d0*cc*myv + quad_form(Gcache(:,:,row), zth(:,ithr), ncomp)*aa
                        do q = 1, ncomp
                            res = res + 2.d0*aa*zth(q,ithr)*ccache(q,row) - 2.d0*cc*zth(q,ithr)*bcache(q,row)
                        end do
                        res = res/sig2
                        do q = 1, ncomp
                            res = res + prior(q)*zth(q,ithr)*zth(q,ithr)
                        end do
                        resid_energy(row) = res
                        ! the two half-data solves the reliability prior is built from
                        do ih = 1, 2
                            do q = 1, ncomp
                                do r = 1, ncomp
                                    Ath(q,r,ithr) = (aa/sig2)*Gall((r-1)*ncomp+q,ii,ih,ithr)
                                end do
                            end do
                            do q = 1, ncomp
                                Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                                zth(q,ithr)   = (cc*real(Reb(q,ii,ih,ithr),dp) &
                                    &- aa*Corr(q,ii,ih,ithr))/sig2
                            end do
                            call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), ncomp)
                            zhalf(row,:,ih) = zth(:,ithr)
                        end do
                    end do
                end do
                nvalid_thr(ithr) = nvalid_thr(ithr) + m
            end do
            !$omp end parallel do
        end do
        nvalid = sum(nvalid_thr)
        write(logfhandle,'(A,F8.1,A,I0)') '>>> FLEX_PCA POLAR EMBED PASS-B SECONDS: ',toc(t_b), &
            &'   particles=',nvalid
        call flush(logfhandle)
        call polar_grid_kill(pg)
        deallocate(xws, xws1, xws2, wrs, eyy, dir_of, cav, sav, dcnt, dptr, dlist, l_zero, rmat_b, rmat_p, eul_b)
        deallocate(pose_rot, pose_shx, pose_shy, pose_e0, pose_e1, pose_s0, pose_s1, pose_n, pose_t, pose_r)
        deallocate(pose_ea, pose_sa, ndmove_thr, dangle_thr, dq0_thr, dq1_thr)
        if( allocated(dnn) ) deallocate(dnn)
        if( allocated(Usall) ) deallocate(Usall, Cfall, Cm0all, c00all)
        deallocate(Ubank, Us, Cf, c00, Cm0, Xs, Reb, Wmat, Gall, Corr, den_v, Ath, zth, Csp)
        deallocate(nvalid_thr, dirs_in_chunk)
    end subroutine embed_accumulate_polar

    !> Recover the in-plane angle in degrees from the (cos,sin) pair the pose block carries.
    pure real function polar_inplane_deg( ca, sa ) result( d )
        real, intent(in) :: ca, sa
        d = atan2(sa, ca) * 180.0 / PI
    end function polar_inplane_deg

    integer function cov_pose_mode()
        integer :: v
        v = 0
        call cov_env_int('SIMPLE_COV_POSE', v)
        cov_pose_mode = v
    end function cov_pose_mode

    !> ---------------------------------------------------------------------------------------------
    !> P1 POSE REFINEMENT FOR ONE PARTICLE: in-plane angle and shift, scored by the MARGINAL
    !> likelihood with the conformation integrated out.
    !>
    !> Why the marginal and not the fitted residual: the failure mode of joint pose+heterogeneity
    !> refinement is that pose absorbs conformation -- a small rotation mimics a low-order
    !> eigenvolume and an alternating optimiser takes it, because the loss does not care which knob
    !> explains the residual. Integrating z out under its prior removes the knob: a pose is preferred
    !> only if it explains the data better AFTER the whole ensemble has been given its best shot at
    !> every competing pose.
    !>
    !> With  A = (a^2/sig2) G + Lambda^-1  and  u = (a b - a^2 c)/sig2,
    !>     -2 log p(y | pose)  =  [e_yy + a^2 e_mm - 2 a my]/sig2  -  u^T A^-1 u  +  log det A + const
    !> and for IN-PLANE and SHIFT moves G is invariant (the noise weight is radial), so `A` and the
    !> log-det are constant across trials. One Cholesky per particle serves the whole search, and
    !> each trial costs one GEMV plus one triangular solve.
    !>
    !> Two further consequences of working in polar coordinates: the in-plane angle is just a
    !> different set of sampling angles, and the shift is a phase multiply that needs no resampling
    !> at all -- so an entire shift grid is scored from one angular resample.
    !>
    !> `itest > 0` injects a deterministic per-particle perturbation first and reports what fraction
    !> of it the search recovers. On IgG-RL and Ribosembly the project poses ARE the ground truth, so
    !> perturb-and-recover is the only honest way to ask whether the score locates the right pose.
    !> ---------------------------------------------------------------------------------------------
    subroutine polar_pose_refine_one( pg, fpl, Usall, Cfall, Cm0all, c00all, nsamp2, ncomp, nk, &
        &ndir, idir, iptcl, wr, prior, sig2, itest, rot_amp, sh_amp, rot_step0, sh_step0, l_guard, &
        &l_p2, nnn, dnn, rmat_pi, rmat_b, jdir, &
        &ca, sa, drot, dshx, dshy, xws, e0, e1, s0, s1, ealt, salt, npose, ntrue, nrej )
        type(polar_grid_t), intent(in)    :: pg
        type(fplane_type),  intent(in)    :: fpl
        integer,            intent(in)    :: nsamp2, ncomp, nk, ndir, idir, iptcl, itest
        real,               intent(in)    :: Usall(nsamp2,0:ncomp,ndir)
        real(dp),           intent(in)    :: Cfall(ncomp*ncomp,nk,ndir), Cm0all(ncomp,nk,ndir)
        real(dp),           intent(in)    :: c00all(nk,ndir)
        real,               intent(in)    :: wr(nk)
        real(dp),           intent(in)    :: prior(ncomp), sig2
        real,               intent(in)    :: rot_amp, sh_amp, rot_step0, sh_step0
        logical,            intent(in)    :: l_guard, l_p2
        integer,            intent(in)    :: nnn, dnn(nnn)
        real,               intent(in)    :: rmat_pi(3,3), rmat_b(3,3,ndir)
        integer,            intent(out)   :: jdir          !< the ACCEPTED bank direction
        real,               intent(inout) :: ca, sa
        real,               intent(out)   :: drot, dshx, dshy
        real,               intent(inout) :: xws(nsamp2)
        real(dp),           intent(inout) :: e0, e1, s0, s1, ealt, salt
        integer,            intent(inout) :: npose, ntrue, nrej
        integer, parameter :: NROT = 7, NSH = 7, NROUND = 3
        !> in-plane trials per candidate DIRECTION. Coarse on purpose: the direction stage only has
        !! to rank directions; the winner then gets the full multi-scale refinement.
        integer, parameter :: NROT_D = 3
        ! Coarsest step. Rotation and shift are NOT interchangeable units: a 1 degree in-plane
        ! rotation displaces the particle edge by mskrad*delta pixels, so for a 200 A particle at
        ! 6 A/px it is worth ~0.3 px of shift. The rotation grid therefore starts coarser.
        real,    parameter :: ROT_STEP0_DEF = 2.0, SH_STEP0_DEF = 1.0
        complex, allocatable :: xw0(:), xw1(:)
        real,    allocatable :: xtrial(:), bq(:)
        real(dp),allocatable :: Ach(:,:), uvec(:), cvec(:)
        real(dp) :: e_mm, aa, acon, best, sc, sc_start, sc_true, shc
        real     :: ca0, sa0, cad, sad, cdel, sdel, px, py, brot, bshx, bshy, rstep, sstep
        real     :: rot_in, shx_in, shy_in, cur_rot, cur_shx, cur_shy, brot0, bshx0, bshy0
        real     :: pr, px_in, py_in
        real(dp) :: h1s, h2s, h1f, h2f, hfull, best2, logdet, logdet2, bestd, e_mm2
        real(dp), allocatable :: Ach2(:,:), cvec2(:)
        real     :: drot2, dshx2, dshy2, caj, saj
        integer  :: q, r, ir, info, it, ish, jsh, iround, ic, jd, jdir0
        allocate(xw0(pg%nsamp), xw1(pg%nsamp), xtrial(nsamp2), bq(0:ncomp))
        allocate(Ach(ncomp,ncomp), uvec(ncomp), cvec(ncomp))
        allocate(Ach2(ncomp,ncomp), cvec2(ncomp))
        shc  = real(fpl%shconst(1), dp)
        acon = 1.d0                                        ! contrast; COV_UNIT_CONTRAST is the default
        aa   = acon*acon
        jdir = idir
        call polar_dir_tables(Cfall, Cm0all, c00all, ncomp, nk, ndir, jdir, wr, prior, sig2, &
            &acon, Ach, cvec, e_mm, logdet, info)
        if( info /= 0 )then
            drot = 0.; dshx = 0.; dshy = 0.
            deallocate(xw0, xw1, xtrial, bq, Ach, uvec, cvec, Ach2, cvec2)
            return
        endif
        ca0 = ca; sa0 = sa
        ! optional injected perturbation -- deterministic in the particle index so the test is
        ! reproducible and carries no RNG state
        rot_in = 0.; shx_in = 0.; shy_in = 0.
        if( itest > 0 )then
            rot_in = rot_amp*(2.0*pose_hash(iptcl, 1) - 1.0)
            shx_in = sh_amp *(2.0*pose_hash(iptcl, 2) - 1.0)
            shy_in = sh_amp *(2.0*pose_hash(iptcl, 3) - 1.0)
            cdel = cos(rot_in*PI/180.); sdel = sin(rot_in*PI/180.)
            ca0  = ca*cdel - sa*sdel
            sa0  = sa*cdel + ca*sdel
        endif
        ! The trial shift is expressed in UNPADDED lattice units, but shconst is defined against the
        ! PADDED grid gen_fplane4rec writes on, so the phase carries the pad factor. Missing it makes
        ! every trial shift half of what it claims to be.
        shc  = shc * real(OSMPL_PAD_FAC, dp)
        ! the search starts from the (possibly perturbed) incoming pose
        brot = 0.; bshx = shx_in; bshy = shy_in
        drot = brot; dshx = bshx; dshy = bshy
        brot0 = brot; bshx0 = bshx; bshy0 = bshy
        jdir0 = jdir
        ! ---- P2: OUT-OF-PLANE DIRECTION SEARCH over the neighbours of the current direction.
        !
        ! Each candidate has its own bank, hence its own G, Cholesky and log-determinant. The score
        ! is the full marginal INCLUDING log det A -- see polar_dir_tables for why that term cannot
        ! be dropped here even though P1 was free to ignore it. Candidates are scored with a coarse
        ! in-plane scan at the incoming shift; the winner then gets the full P1 refinement below.
        if( l_p2 )then
            bestd = huge(1.d0)
            do ic = 1, nnn
                jd = dnn(ic)
                call polar_dir_tables(Cfall, Cm0all, c00all, ncomp, nk, ndir, jd, wr, prior, sig2, &
                    &acon, Ach2, cvec2, e_mm2, logdet2, info)
                if( info /= 0 ) cycle
                ! the relative in-plane angle is different at every candidate direction
                call polar_relative_inplane(rmat_pi, rmat_b(:,:,jd), caj, saj)
                if( itest > 0 )then
                    cdel = cos(rot_in*PI/180.); sdel = sin(rot_in*PI/180.)
                    cad  = caj*cdel - saj*sdel
                    sad  = saj*cdel + caj*sdel
                    caj  = cad; saj = sad
                endif
                do it = 1, NROT_D
                    cur_rot = real(it - (NROT_D+1)/2)*rot_step0
                    cdel = cos(cur_rot*PI/180.); sdel = sin(cur_rot*PI/180.)
                    cad  = caj*cdel - saj*sdel
                    sad  = saj*cdel + caj*sdel
                    call polar_sample_at_pose(fpl%cmplx_plane, fpl%transfer_plane, pg, cad, sad, &
                        &-bshx*real(shc), -bshy*real(shc), xw0)
                    call polar_pose_score_halves(pg, xw0, Usall(1,0,jd), nsamp2, ncomp, Ach2, &
                        &cvec2, sig2, acon, e_mm2, xtrial, uvec, h1f, h2f, hfull)
                    sc = hfull + logdet2
                    if( sc < bestd )then
                        bestd = sc; jdir = jd; ca0 = caj; sa0 = saj; brot = cur_rot
                        Ach = Ach2; cvec = cvec2; e_mm = e_mm2; logdet = logdet2
                    endif
                end do
            end do
            drot = brot
        endif
        ! The trial shift is expressed in UNPADDED lattice units, but shconst is defined against the
        ! PADDED grid gen_fplane4rec writes on, so the phase carries the pad factor. Missing it makes
        ! every trial shift half of what it claims to be.
        ! score at the starting pose, for the "did the search find anything" diagnostic
        call polar_sample_at_pose(fpl%cmplx_plane, fpl%transfer_plane, pg, ca0, sa0, &
            &-bshx*real(shc), -bshy*real(shc), xw0)
        call polar_pose_score_halves(pg, xw0, Usall(1,0,jdir), nsamp2, ncomp, Ach, cvec, &
            &sig2, acon, e_mm, xtrial, uvec, h1s, h2s, hfull)
        sc_start = hfull
        best     = sc_start
        best2    = h2s
        drot2    = brot; dshx2 = bshx; dshy2 = bshy
        ! multi-scale coordinate descent: angle grid, then shift grid, halving the step each round
        rstep = rot_step0
        sstep = sh_step0
        do iround = 1, NROUND
            ! ---- angle
            do it = 1, NROT
                cur_rot = brot + real(it - (NROT+1)/2)*rstep
                cdel = cos(cur_rot*PI/180.); sdel = sin(cur_rot*PI/180.)
                cad  = ca0*cdel - sa0*sdel
                sad  = sa0*cdel + ca0*sdel
                call polar_sample_at_pose(fpl%cmplx_plane, fpl%transfer_plane, pg, cad, sad, &
                    &-bshx*real(shc), -bshy*real(shc), xw0)
                call polar_pose_score_halves(pg, xw0, Usall(1,0,jdir), nsamp2, ncomp, Ach, cvec, &
                    &sig2, acon, e_mm, xtrial, uvec, h1f, h2f, hfull)
                sc = hfull                             ! selection always uses the full plane
                if( sc < best )then
                    best = sc; drot = cur_rot; dshx = bshx; dshy = bshy
                endif
                if( h2f < best2 )then
                    best2 = h2f; drot2 = cur_rot; dshx2 = bshx; dshy2 = bshy
                endif
            end do
            brot = drot
            ! ---- shift, at the accepted angle: ONE resample, the whole grid from phase multiplies
            cdel = cos(brot*PI/180.); sdel = sin(brot*PI/180.)
            cad  = ca0*cdel - sa0*sdel
            sad  = sa0*cdel + ca0*sdel
            call polar_sample_at_pose(fpl%cmplx_plane, fpl%transfer_plane, pg, cad, sad, 0.0, 0.0, xw0)
            do ish = 1, NSH
                do jsh = 1, NSH
                    cur_shx = bshx + real(ish - (NSH+1)/2)*sstep
                    cur_shy = bshy + real(jsh - (NSH+1)/2)*sstep
                    px = -cur_shx*real(shc); py = -cur_shy*real(shc)
                    call polar_apply_shift(pg, cad, sad, px, py, xw0, xw1)
                    call polar_pose_score_halves(pg, xw1, Usall(1,0,jdir), nsamp2, ncomp, Ach, &
                        &cvec, sig2, acon, e_mm, xtrial, uvec, h1f, h2f, hfull)
                    sc = hfull
                    if( sc < best )then
                        best = sc; drot = brot; dshx = cur_shx; dshy = cur_shy
                    endif
                    if( h2f < best2 )then
                        best2 = h2f; drot2 = brot; dshx2 = cur_shx; dshy2 = cur_shy
                    endif
                end do
            end do
            bshx = dshx; bshy = dshy
            rstep = 0.5*rstep; sstep = 0.5*sstep
        end do
        ! ---- SPLIT-HALF ACCEPTANCE.
        !
        ! The pose is CHOSEN on the full plane -- restricting selection to one half costs more
        ! recovery than the extra independence buys (measured: 3.45 deg -> 1.69 rather than 1.36).
        ! It is then accepted only if BOTH interleaved halves independently score it better than the
        ! incoming pose. Measured on IgG-10k this halves the damage done to already-correct poses
        ! (0.71 -> 0.35 px) at no cost to recovery from real error (1.6 px -> 0.41 either way).
        !
        ! It does NOT fix the angular case, and cannot: the ~1.0 deg of residual movement is below
        ! the 1.2 deg identifiability floor, so no test built on this data can tell it from signal.
        if( l_guard .and. (abs(drot-brot0) > TINY .or. abs(dshx-bshx0) > TINY .or. &
            &abs(dshy-bshy0) > TINY .or. jdir /= jdir0) )then
            cdel = cos(drot*PI/180.); sdel = sin(drot*PI/180.)
            cad  = ca0*cdel - sa0*sdel
            sad  = sa0*cdel + ca0*sdel
            call polar_sample_at_pose(fpl%cmplx_plane, fpl%transfer_plane, pg, cad, sad, &
                &-dshx*real(shc), -dshy*real(shc), xw0)
            call polar_pose_score_halves(pg, xw0, Usall(1,0,jdir), nsamp2, ncomp, Ach, cvec, &
                &sig2, acon, e_mm, xtrial, uvec, h1f, h2f, hfull)
            if( .not. (h1f < h1s .and. h2f < h2s) )then
                drot = brot0; dshx = bshx0; dshy = bshy0     ! rejected: keep the incoming pose
                jdir = jdir0                                 ! ... including the direction
                nrej = nrej + 1
            endif
        endif
        ! In test mode, score the TRUE pose too. If truth does not beat what the search found, the
        ! score itself does not identify pose and no amount of search will fix that -- this is the
        ! number that separates "search too coarse" from "objective uninformative".
        sc_true = 0.d0
        if( itest > 0 )then
            cdel = cos(-rot_in*PI/180.); sdel = sin(-rot_in*PI/180.)
            cad  = ca0*cdel - sa0*sdel
            sad  = sa0*cdel + ca0*sdel
            call polar_sample_at_pose(fpl%cmplx_plane, fpl%transfer_plane, pg, cad, sad, 0.0, 0.0, xw0)
            sc_true = polar_pose_score(pg, xw0, Usall(1,0,jdir), nsamp2, ncomp, Ach, cvec, prior, &
                &sig2, acon, e_mm, xtrial, bq, uvec)
        endif
        ! adopt the accepted pose: resample there and hand the samples back
        cdel = cos(drot*PI/180.); sdel = sin(drot*PI/180.)
        cad  = ca0*cdel - sa0*sdel
        sad  = sa0*cdel + ca0*sdel
        px   = -dshx*real(shc); py = -dshy*real(shc)
        call polar_sample_at_pose(fpl%cmplx_plane, fpl%transfer_plane, pg, cad, sad, px, py, xw0)
        do q = 1, pg%nsamp
            xws(2*q-1) = pg%sqwq(q)*real(xw0(q))
            xws(2*q)   = pg%sqwq(q)*aimag(xw0(q))
        end do
        ca = cad; sa = sad
        ! statistics
        npose = npose + 1
        if( itest > 0 )then
            e0 = e0 + real(rot_in,dp)**2
            e1 = e1 + real(rot_in + drot, dp)**2
            s0 = s0 + real(shx_in,dp)**2 + real(shy_in,dp)**2
            s1 = s1 + real(dshx,dp)**2 + real(dshy,dp)**2
            ! ntrue counts particles where the TRUE pose scores better than the search's answer
            if( sc_true < best ) ntrue = ntrue + 1
        else if( COV_PROJ_PERTURB_ROT > 0. .or. COV_PROJ_PERTURB_SH > 0. )then
            ! the project's poses were degraded by a KNOWN amount; report what is left of it. The
            ! sign of the in-plane convention is not assumed -- both are accumulated and the log
            ! prints whichever the search actually produced.
            pr = COV_PROJ_PERTURB_ROT*(2.0*pose_hash(iptcl, 1) - 1.0)
            px_in = COV_PROJ_PERTURB_SH*(2.0*pose_hash(iptcl, 2) - 1.0)
            py_in = COV_PROJ_PERTURB_SH*(2.0*pose_hash(iptcl, 3) - 1.0)
            e0 = e0 + real(pr,dp)**2
            e1 = e1 + real(pr + drot, dp)**2
            s0 = s0 + real(px_in,dp)**2 + real(py_in,dp)**2
            s1 = s1 + real(px_in + dshx, dp)**2 + real(py_in + dshy, dp)**2
            ealt = ealt + real(pr - drot, dp)**2
            salt = salt + real(px_in - dshx, dp)**2 + real(py_in - dshy, dp)**2
            ntrue = ntrue + 1
        else
            e1 = e1 + real(drot,dp)**2
            s1 = s1 + real(dshx,dp)**2 + real(dshy,dp)**2
            if( best < sc_start ) ntrue = ntrue + 1
        endif
        deallocate(xw0, xw1, xtrial, bq, Ach, uvec, cvec, Ach2, cvec2)
    end subroutine polar_pose_refine_one

    !> Everything in the pose objective that depends on the candidate DIRECTION: the Gram, its
    !! Cholesky, the log-determinant, the mean cross term and <T mu, T mu>.
    !!
    !! For P1 this was computed once per particle because `G` is invariant under in-plane rotation
    !! and shift. It is NOT invariant across directions, so P2 pays this per candidate -- and, more
    !! importantly, must KEEP the log-det term. It cancels between in-plane trials and does not
    !! cancel between directions; dropping it biases the search toward directions where the basis
    !! happens to be more expressive, which is the term a naive implementation silently loses.
    subroutine polar_dir_tables( Cfall, Cm0all, c00all, ncomp, nk, ndir, jd, wr, prior, sig2, &
        &acon, Ach, cvec, e_mm, logdet, info )
        integer,  intent(in)  :: ncomp, nk, ndir, jd
        real(dp), intent(in)  :: Cfall(ncomp*ncomp,nk,ndir), Cm0all(ncomp,nk,ndir), c00all(nk,ndir)
        real,     intent(in)  :: wr(nk)
        real(dp), intent(in)  :: prior(ncomp), sig2, acon
        real(dp), intent(out) :: Ach(ncomp,ncomp), cvec(ncomp), e_mm, logdet
        integer,  intent(out) :: info
        real(dp) :: aa, w
        integer  :: q, r, ir
        aa   = acon*acon
        e_mm = 0.d0
        cvec = 0.d0
        Ach  = 0.d0
        do ir = 1, nk
            w    = real(wr(ir),dp)
            e_mm = e_mm + w*c00all(ir,jd)
            do q = 1, ncomp
                cvec(q) = cvec(q) + w*Cm0all(q,ir,jd)
                do r = 1, ncomp
                    Ach(q,r) = Ach(q,r) + w*Cfall((r-1)*ncomp+q,ir,jd)
                end do
            end do
        end do
        Ach = (aa/sig2)*Ach
        do q = 1, ncomp
            Ach(q,q) = Ach(q,q) + prior(q)
        end do
        call dpotrf('U', ncomp, Ach, ncomp, info)
        logdet = 0.d0
        if( info == 0 )then
            do q = 1, ncomp
                logdet = logdet + 2.d0*log(max(Ach(q,q), DTINY))
            end do
        endif
    end subroutine polar_dir_tables

    !> The same objective evaluated on each interleaved half of the sample list. Used ONLY to decide
    !! whether to accept a move, never to choose it: a pose picked as the argmax over ~170 trials is
    !! biased to fit noise, and the measured per-particle floor (~0.4 px of displacement at SNR 0.01)
    !! is large enough that refining an already-correct pose makes it WORSE. Requiring both halves to
    !! independently prefer the move is a threshold-free guard against exactly that.
    subroutine polar_pose_score_halves( pg, xw, Us, nsamp2, ncomp, Ach, cvec, sig2, acon, e_mm, &
        &xtrial, uvec, sc1, sc2, scfull )
        type(polar_grid_t), intent(in)    :: pg
        complex,            intent(in)    :: xw(:)
        integer,            intent(in)    :: nsamp2, ncomp
        real,               intent(in)    :: Us(nsamp2,0:ncomp)
        real(dp),           intent(in)    :: Ach(ncomp,ncomp), cvec(ncomp), sig2, acon, e_mm
        real,               intent(inout) :: xtrial(nsamp2)
        real(dp),           intent(inout) :: uvec(ncomp)
        real(dp),           intent(out)   :: sc1, sc2, scfull
        real(dp) :: b1(0:ncomp), b2(0:ncomp), aa, quad
        integer  :: j, q
        do j = 1, pg%nsamp
            xtrial(2*j-1) = pg%sqwq(j)*real(xw(j))
            xtrial(2*j)   = pg%sqwq(j)*aimag(xw(j))
        end do
        b1 = 0.d0; b2 = 0.d0
        do q = 0, ncomp
            do j = 1, nsamp2
                if( pg%hrow(j) == 1 )then
                    b1(q) = b1(q) + real(Us(j,q),dp)*real(xtrial(j),dp)
                else
                    b2(q) = b2(q) + real(Us(j,q),dp)*real(xtrial(j),dp)
                endif
            end do
        end do
        aa = acon*acon
        do q = 1, ncomp
            uvec(q) = (acon*b1(q) - 0.5d0*aa*cvec(q))/sig2
        end do
        call dtrsv('U', 'T', 'N', ncomp, Ach, ncomp, uvec, 1)
        quad = 0.d0
        do q = 1, ncomp
            quad = quad + uvec(q)*uvec(q)
        end do
        sc1 = (0.5d0*aa*e_mm - 2.d0*acon*b1(0))/sig2 - quad
        do q = 1, ncomp
            uvec(q) = (acon*b2(q) - 0.5d0*aa*cvec(q))/sig2
        end do
        call dtrsv('U', 'T', 'N', ncomp, Ach, ncomp, uvec, 1)
        quad = 0.d0
        do q = 1, ncomp
            quad = quad + uvec(q)*uvec(q)
        end do
        sc2 = (0.5d0*aa*e_mm - 2.d0*acon*b2(0))/sig2 - quad
        ! The FULL-plane marginal is NOT sc1+sc2: the quadratic form is (u1+u2)^T A^-1 (u1+u2), and
        ! dropping the cross term halves it relative to the linear term. Selection must use this one;
        ! sc1/sc2 exist only for the split-half acceptance test.
        do q = 1, ncomp
            uvec(q) = (acon*(b1(q) + b2(q)) - aa*cvec(q))/sig2
        end do
        call dtrsv('U', 'T', 'N', ncomp, Ach, ncomp, uvec, 1)
        quad = 0.d0
        do q = 1, ncomp
            quad = quad + uvec(q)*uvec(q)
        end do
        scfull = (aa*e_mm - 2.d0*acon*(b1(0) + b2(0)))/sig2 - quad
    end subroutine polar_pose_score_halves

    !> -2 log p(y | pose) up to terms constant across trials of the same particle.
    real(dp) function polar_pose_score( pg, xw, Us, nsamp2, ncomp, Ach, cvec, prior, sig2, acon, &
        &e_mm, xtrial, bq, uvec ) result( sc )
        type(polar_grid_t), intent(in)    :: pg
        complex,            intent(in)    :: xw(:)
        integer,            intent(in)    :: nsamp2, ncomp
        real,               intent(in)    :: Us(nsamp2,0:ncomp)
        real(dp),           intent(in)    :: Ach(ncomp,ncomp), cvec(ncomp), prior(ncomp), sig2
        real(dp),           intent(in)    :: acon, e_mm
        real,               intent(inout) :: xtrial(nsamp2), bq(0:ncomp)
        real(dp),           intent(inout) :: uvec(ncomp)
        integer  :: j, q
        real(dp) :: aa, myv, quad
        do j = 1, pg%nsamp
            xtrial(2*j-1) = pg%sqwq(j)*real(xw(j))
            xtrial(2*j)   = pg%sqwq(j)*aimag(xw(j))
        end do
        call sgemv('T', nsamp2, ncomp+1, 1.0, Us, nsamp2, xtrial, 1, 0.0, bq, 1)
        aa  = acon*acon
        myv = real(bq(0), dp)
        do q = 1, ncomp
            uvec(q) = (acon*real(bq(q),dp) - aa*cvec(q))/sig2
        end do
        call dtrsv('U', 'T', 'N', ncomp, Ach, ncomp, uvec, 1)
        quad = 0.d0
        do q = 1, ncomp
            quad = quad + uvec(q)*uvec(q)
        end do
        sc = (aa*e_mm - 2.d0*acon*myv)/sig2 - quad
    end function polar_pose_score

    !> ---------------------------------------------------------------------------------------------
    !> Degrade the PROJECT's in-plane angles and shifts before anything reads them.
    !>
    !> SIMPLE_COV_POSE_TEST perturbs only inside the embedding, which measures whether the score can
    !> find a pose but leaves the model itself built from correct poses. That is the easy half of the
    !> question. This perturbs `build%spproj_field` at the very start of the run, so the mean, the
    !> covariance columns, the basis AND the embedding are all built from degraded poses -- the
    !> feedback loop that a real ab-initio run has and the embed-only test does not.
    !>
    !> The perturbation is hashed from the particle index (same hash as the embed-only test), so it
    !> is exactly reproducible and its magnitude is known, which is what makes the delivered basis
    !> comparable against an unperturbed run.
    !> ---------------------------------------------------------------------------------------------
    subroutine cov_perturb_project_poses( build, pinds, nptcls )
        type(builder), intent(inout) :: build
        integer,       intent(in)    :: pinds(:), nptcls
        type(ori) :: o
        integer :: irot, ish, idir, i, j, iptcl
        real    :: ramp, samp, e3, dr, dx, dy, sh(2), dd, psi, cps, sps, cdd, sdd
        real    :: rm(3,3), rm2(3,3), ax(3), v(3), cr(3)
        real(dp):: acc_r, acc_s, acc_d
        irot = 0; call cov_env_int('SIMPLE_COV_POSE_PERTURB_ROT', irot)
        ish  = 0; call cov_env_int('SIMPLE_COV_POSE_PERTURB_SH',  ish)
        idir = 0; call cov_env_int('SIMPLE_COV_POSE_PERTURB_DIR', idir)
        if( irot <= 0 .and. ish <= 0 .and. idir <= 0 ) return
        COV_PROJ_PERTURB_ROT = 0.1*real(irot)
        COV_PROJ_PERTURB_SH  = 0.1*real(ish)
        COV_PROJ_PERTURB_DIR = 0.1*real(idir)
        ramp = 0.1*real(irot)
        samp = 0.1*real(ish)
        acc_r = 0.d0; acc_s = 0.d0; acc_d = 0.d0
        if( COV_PROJ_PERTURB_DIR > 0. )then
            if( allocated(COV_TRUTH_NRM) ) deallocate(COV_TRUTH_NRM)
            allocate(COV_TRUTH_NRM(3,nptcls))
        endif
        do i = 1, nptcls
            iptcl = pinds(i)
            dr = ramp*(2.0*pose_hash(i, 1) - 1.0)
            dx = samp*(2.0*pose_hash(i, 2) - 1.0)
            dy = samp*(2.0*pose_hash(i, 3) - 1.0)
            if( COV_PROJ_PERTURB_DIR > 0. )then
                ! tilt the VIEWING DIRECTION by dd about an axis lying in the projection plane, so
                ! the normal moves by exactly dd. Rodrigues on each row of the frame.
                call build%spproj_field%get_ori(iptcl, o)
                rm = o%get_mat()
                COV_TRUTH_NRM(:,i) = rm(3,:)
                dd  = COV_PROJ_PERTURB_DIR*(2.0*pose_hash(i, 4) - 1.0)
                psi = 360.0*pose_hash(i, 5)
                cps = cos(psi*PI/180.); sps = sin(psi*PI/180.)
                ax  = cps*rm(1,:) + sps*rm(2,:)                    ! axis in the plane
                cdd = cos(dd*PI/180.); sdd = sin(dd*PI/180.)
                do j = 1, 3
                    v          = rm(j,:)
                    cr(1)      = ax(2)*v(3) - ax(3)*v(2)
                    cr(2)      = ax(3)*v(1) - ax(1)*v(3)
                    cr(3)      = ax(1)*v(2) - ax(2)*v(1)
                    rm2(j,:)   = v*cdd + cr*sdd + ax*dot_product(ax,v)*(1.0 - cdd)
                end do
                call o%ori_from_rotmat(rm2, .true.)
                call build%spproj_field%set_euler(iptcl, o%get_euler())
                acc_d = acc_d + real(dd,dp)**2
                call o%kill
            endif
            e3 = build%spproj_field%e3get(iptcl)
            call build%spproj_field%e3set(iptcl, e3 + dr)
            sh = build%spproj_field%get_2Dshift(iptcl)
            call build%spproj_field%set_shift(iptcl, sh + [dx, dy])
            acc_r = acc_r + real(dr,dp)**2
            acc_s = acc_s + real(dx,dp)**2 + real(dy,dp)**2
        end do
        write(logfhandle,'(A,F8.4,A,F8.4,A,F8.4)') '>>> FLEX_PCA PROJECT POSES DEGRADED: in-plane rms=', &
            &sqrt(acc_r/real(nptcls,dp)),' deg   shift rms=',sqrt(acc_s/real(nptcls,dp)), &
            &' px   direction rms=',sqrt(acc_d/real(nptcls,dp))
        write(logfhandle,'(A)') '>>> FLEX_PCA   (every stage below now sees the degraded poses)'
        call flush(logfhandle)
    end subroutine cov_perturb_project_poses

    !> Deterministic per-particle uniform-ish value in [0,1) for the perturbation test. A hash, not
    !! an RNG, so the injected perturbation is identical across runs and across thread counts.
    pure real function pose_hash( i, k ) result( h )
        integer, intent(in) :: i, k
        integer(kind=8) :: x
        x = int(i, 8)*2654435761_8 + int(k, 8)*40503_8
        x = iand(ieor(x, ishft(x, -13)), 2147483647_8)
        h = real(modulo(x, 100000_8)) / 100000.0
    end function pose_hash

    !> One ring's contribution to the Gram of the basis (columns 1..ncomp of Us) and to the mean
    !! cross term (column 0). `Cout` is the full ncomp x ncomp block flattened column-major, `Mout`
    !! the ncomp-vector <U_q, T mu>.
    subroutine polar_ring_gram( Us, ldu, ncomp, row0, nrow, Csp, Cout, Mout )
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

    real(dp) function polar_ring_selfpower( Us, ldu, row0, nrow )
        integer, intent(in) :: ldu, row0, nrow
        real,    intent(in) :: Us(ldu,0:*)
        integer :: j
        polar_ring_selfpower = 0.d0
        do j = 2*row0-1, 2*(row0+nrow-1)
            polar_ring_selfpower = polar_ring_selfpower + real(Us(j,0),dp)*real(Us(j,0),dp)
        end do
    end function polar_ring_selfpower

    !> <y,y> in the polar measure. xws already carries sqrt(wq) and the CTF adjoint, so the CTF has
    !! to be divided back out ring by ring to recover the observation's own energy.
    real(dp) function polar_self_energy( xws, wr, pg ) result( e )
        real,               intent(in) :: xws(:), wr(:)
        type(polar_grid_t), intent(in) :: pg
        integer  :: ir, j
        real(dp) :: acc
        e = 0.d0
        do ir = 1, pg%nk
            if( real(wr(ir),dp) <= DTINY ) cycle
            acc = 0.d0
            do j = 2*pg%rbeg(ir)-1, 2*pg%rend(ir)
                acc = acc + real(xws(j),dp)*real(xws(j),dp)
            end do
            e = e + acc / real(wr(ir),dp)
        end do
    end function polar_self_energy

    pure real(dp) function sum_dp_safe( acc, n ) result( v )
        real(dp), intent(in) :: acc
        integer,  intent(in) :: n
        v = acc / real(max(1,n), dp)
    end function sum_dp_safe

    !> thin wrapper so the OpenMP body stays readable; folds sqrt(wq) into the stored samples
    !> resample a stored half-plane by an in-plane rotation into the bank frame: unit-tap 2D KB
    !! at pf-multiples with per-tap Friedel, per-axis normalized weights, OOB taps dropped --
    !! the polar former's interpolation scheme (polar_interp_plane) on the Cartesian lattice.
    !! Positions outside the nyq disk come back ZERO (nothing downstream reads them).
    subroutine align_halfplane_inplane( frlims, nyq_eff, src, ca, sa, dst )
        integer, intent(in)  :: frlims(3,2), nyq_eff
        complex, intent(in)  :: src(frlims(1,1):frlims(1,2), frlims(2,1):0)
        real,    intent(in)  :: ca, sa
        complex, intent(out) :: dst(frlims(1,1):frlims(1,2), frlims(2,1):0)
        type(kbinterpol) :: kbwin
        real    :: hu, ku, w, wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM)
        integer :: win(2,3), h, k, hlo2, hhi2, klo2, hx, ky, ix, iy, nd, pf
        complex :: acc, cv
        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        pf    = OSMPL_PAD_FAC
        hlo2  = ceil_div (frlims(1,1), pf); hhi2 = floor_div(frlims(1,2), pf)
        klo2  = ceil_div (frlims(2,1), pf)
        nd    = nyq_eff*(nyq_eff+1)
        dst   = CMPLX_ZERO
        do k = klo2, 0
            do h = hlo2, hhi2
                if( h*h + k*k > nd ) cycle
                hu =  h*ca + k*sa
                ku = -h*sa + k*ca
                call latent_projection_weights(kbwin, [hu, ku, 0.], win, wx, wy, wz)
                acc = CMPLX_ZERO
                do iy = 1, LATENT_WDIM
                    ky = win(1,2) + iy - 1
                    do ix = 1, LATENT_WDIM
                        hx = win(1,1) + ix - 1
                        w  = wx(ix)*wy(iy)
                        if( pf*ky <= 0 )then
                            if( pf*hx < frlims(1,1) .or. pf*hx > frlims(1,2) .or. pf*ky < frlims(2,1) ) cycle
                            cv = src(pf*hx, pf*ky)
                        else
                            if( -pf*hx < frlims(1,1) .or. -pf*hx > frlims(1,2) .or. -pf*ky < frlims(2,1) ) cycle
                            cv = conjg(src(-pf*hx, -pf*ky))
                        endif
                        acc = acc + w*cv
                    end do
                end do
                dst(pf*h, pf*k) = acc
            end do
        end do
    end subroutine align_halfplane_inplane

    subroutine polar_sample_particle_packed( fpl, pg, ca, sa, xws, wr, hfpw, hfcnt, tazim, xws1, xws2 )
        type(fplane_type),  intent(in)    :: fpl
        type(polar_grid_t), intent(in)    :: pg
        real,               intent(in)    :: ca, sa
        real,               intent(out)   :: xws(:)
        real,               intent(out)   :: wr(:)
        real(dp),           intent(inout) :: hfpw, hfcnt
        real,               intent(out)   :: tazim
        !> packed forms of the two lattice-parity half-fields, for the reliability prior
        real, optional,     intent(out)   :: xws1(:), xws2(:)
        complex, allocatable :: xw(:), xw1(:), xw2(:)
        real(dp) :: pw, cnt
        integer  :: j
        allocate(xw(pg%nsamp), xw1(pg%nsamp), xw2(pg%nsamp))
        if( present(xws1) )then
            call polar_sample_particle(fpl%cmplx_plane, fpl%transfer_plane, pg, ca, sa, xw, wr, &
                &pw, cnt, tazim, xw1, xw2)
        else
            call polar_sample_particle(fpl%cmplx_plane, fpl%transfer_plane, pg, ca, sa, xw, wr, &
                &pw, cnt, tazim)
        endif
        do j = 1, pg%nsamp
            xws(2*j-1) = pg%sqwq(j)*real(xw(j))
            xws(2*j)   = pg%sqwq(j)*aimag(xw(j))
        end do
        if( present(xws1) )then
            do j = 1, pg%nsamp
                xws1(2*j-1) = pg%sqwq(j)*real(xw1(j))
                xws1(2*j)   = pg%sqwq(j)*aimag(xw1(j))
                xws2(2*j-1) = pg%sqwq(j)*real(xw2(j))
                xws2(2*j)   = pg%sqwq(j)*aimag(xw2(j))
            end do
        endif
        hfpw  = hfpw  + pw
        hfcnt = hfcnt + cnt
        deallocate(xw, xw1, xw2)
    end subroutine polar_sample_particle_packed

    !> Banded mean projection for the polar E-step: reconstructor%project_fplane's numerics --
    !! the SAME banded (h,k) sweep, apod_mat_3d interpolation weights (including their final
    !! global renormalization), per-sample Friedel conjugation and transfer multiply -- with the
    !! per-call full-plane work removed. project_fplane zero-fills the whole PADDED plane and
    !! copies the reference ctfsq and transfer planes into the output on EVERY call; at the
    !! native padded lattice that is several MB of memory traffic per particle, which measured
    !! as ~80% of the polar E-step's project bucket, all spent on values the polar branch never
    !! reads (only mean_fpl%cmplx_plane is consumed, by the residual subtraction). Here the
    !! plane is zero-filled once at (re)allocation; every call rewrites exactly the in-band disc
    !! samples, and out-of-disc positions stay zero -- the same invariant the Cartesian former's
    !! ensure_latent_projection_plane relies on. The interpolated values are bit-identical to
    !! project_fplane's (same expressions, same kbwin), so the residual planes the M-step
    !! consumes are unchanged wherever the mean is nonzero and unchanged-because-zero elsewhere.
    subroutine project_fplane_mean_banded( rec, o, fpl_ref, fpl_out )
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

    !> Exact Cartesian statistics of the low-k shells for the HYBRID polar E-step, added on top
    !! of the ring statistics. Per lattice position the KB window geometry is computed once and
    !! all ncomp+1 volumes are gathered through it (the Cartesian former's hoist); the data value,
    !! CTF/whitening transfer and quadrature weight (1 per lattice point) are exactly the
    !! Cartesian former's, so the shells this covers contribute to G/b/c/e_mm/myv precisely what
    !! project_fplanes_mean_basis + cov_herm_inner would contribute for them -- including the DC
    !! sample. This is what removes the ring quadrature's multiplicative posterior-variance bias:
    !! after whitening the low-k shells still anchor the latent scale, and rings sample them
    !! worst (few samples, steep integrand).
    subroutine polar_hybrid_exact_accum( rec0, recs, ncomp, o, fpl, hex, kex, npos, &
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
    subroutine subtract_mean_banded( fpl, mean_fpl, a, rec_nyq )
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



    !>  Cross-halfset basis alignment for the held-out embedding. Both eigenbases are
    !!  first normalized to unit real-space norm, then M(i,j) = <U_ref_i, U_tgt_j>.
    !!  A latent expressed in the TARGET basis is mapped into the REFERENCE frame by
    !!  z_ref = M z_tgt (since x ~ U_tgt z_tgt and z_ref = U_ref^T x). The singular
    !!  values of M are the cosines of the principal angles between the two subspaces,
    !!  i.e. a gold-standard measure of how many latent dimensions actually reproduce
    !!  across independent halves -- unlike per-component FSC, this cannot be fooled by
    !!  a shared basis, because the two bases are estimated from disjoint particles.
    subroutine align_basis_to_reference( ref_imgs, nref_c, tgt_imgs, ntgt_c, M, svals )
        integer,     intent(in)    :: nref_c, ntgt_c
        type(image), intent(inout) :: ref_imgs(nref_c), tgt_imgs(ntgt_c)
        real(dp), allocatable, intent(out) :: M(:,:), svals(:)
        real, pointer :: rmat_i(:,:,:), rmat_j(:,:,:)
        real(dp), allocatable :: nrm_r(:), nrm_t(:), Mwork(:,:), V(:,:), ev(:)
        integer  :: i, j, nrot, nsv
        allocate(M(nref_c,ntgt_c), source=0.d0)
        allocate(nrm_r(nref_c), nrm_t(ntgt_c), source=0.d0)
        do i = 1, nref_c
            call ref_imgs(i)%get_rmat_ptr(rmat_i)
            nrm_r(i) = sqrt(max(sum(real(rmat_i,dp)*real(rmat_i,dp)), DTINY))
        end do
        do j = 1, ntgt_c
            call tgt_imgs(j)%get_rmat_ptr(rmat_j)
            nrm_t(j) = sqrt(max(sum(real(rmat_j,dp)*real(rmat_j,dp)), DTINY))
        end do
        do i = 1, nref_c
            call ref_imgs(i)%get_rmat_ptr(rmat_i)
            do j = 1, ntgt_c
                call tgt_imgs(j)%get_rmat_ptr(rmat_j)
                M(i,j) = sum(real(rmat_i,dp)*real(rmat_j,dp)) / (nrm_r(i)*nrm_t(j))
            end do
        end do
        ! principal-angle cosines = singular values of M, via the eigenvalues of M^T M
        nsv = min(nref_c, ntgt_c)
        allocate(Mwork(ntgt_c,ntgt_c), V(ntgt_c,ntgt_c), ev(ntgt_c), svals(nsv))
        Mwork = matmul(transpose(M), M)
        call jacobi(Mwork, ntgt_c, ntgt_c, ev, V, nrot)
        call eigsrt(ev, V, ntgt_c, ntgt_c)
        do i = 1, nsv
            svals(i) = sqrt(max(ev(i), 0.d0))
        end do
        deallocate(nrm_r, nrm_t, Mwork, V, ev)
    end subroutine align_basis_to_reference



    !> Bagged basis: pool two independently fitted eigenbases, keep the top-k singular directions.
    !!
    !!  Nominally equivalent fits of the same data vary 5x in explanatory power, so the variance is
    !!  in the fit, not the sample. Directions the two fits AGREE on add coherently in the pooled
    !!  Gram and survive truncation; directions either invented on its own do not. Ceiling is the
    !!  basis-error share of latent error, ~38 % -- the rest is per-particle estimation noise, which
    !!  the deconvolution in simple_flex_pca_model::gmm_state_weights targets instead.
    !!
    !!  Inputs are compared and combined over the WHOLE rmat, matching align_basis_to_reference.
    subroutine bag_basis_pool( imgs_a, na_c, eig_a, imgs_b, nb_c, eig_b, ncomp_out, pooled, eig_pooled )
        integer,     intent(in)    :: na_c, nb_c, ncomp_out
        type(image), intent(inout) :: imgs_a(na_c), imgs_b(nb_c)
        real(dp),    intent(in)    :: eig_a(na_c), eig_b(nb_c)
        type(image), allocatable, intent(out) :: pooled(:)
        real(dp),    allocatable, intent(out) :: eig_pooled(:)
        real, pointer :: rmat_i(:,:,:), rmat_j(:,:,:), rmat_o(:,:,:)
        real(dp), allocatable :: G(:,:), V(:,:), ev(:), nrm(:), eigin(:)
        integer  :: m, k, i, j, q, nrot, ldim(3)
        real(dp) :: onrm
        real     :: smpd
        m = na_c + nb_c
        k = min(ncomp_out, m)
        if( k < 1 ) THROW_HARD('flex_pca bagging asked for a rank below 1')
        allocate(nrm(m), eigin(m), G(m,m), V(m,m), ev(m))
        do j = 1, m
            call pick(j, rmat_j)
            nrm(j)   = sqrt(max(sum(real(rmat_j,dp)*real(rmat_j,dp)), DTINY))
            eigin(j) = merge(eig_a(min(j,na_c)), eig_b(max(1,j-na_c)), j <= na_c)
        end do
        ! Gram of the unit-normalised union: equal weight per direction is what makes this a bag,
        ! since eigenvalue weighting would reinstate a ranking measured not to track signal.
        do i = 1, m
            call pick(i, rmat_i)
            do j = i, m
                call pick(j, rmat_j)
                G(i,j) = sum(real(rmat_i,dp)*real(rmat_j,dp)) / (nrm(i)*nrm(j))
                G(j,i) = G(i,j)
            end do
        end do
        call jacobi(G, m, m, ev, V, nrot)
        call eigsrt(ev, V, m, m)
        allocate(pooled(k), eig_pooled(k))
        ldim = imgs_a(1)%get_ldim()
        smpd = imgs_a(1)%get_smpd()
        do q = 1, k
            call pooled(q)%new(ldim, smpd)
            call pooled(q)%get_rmat_ptr(rmat_o)
            rmat_o = 0.
            do j = 1, m
                call pick(j, rmat_j)
                rmat_o = rmat_o + real(V(j,q)/nrm(j)) * rmat_j
            end do
            onrm = sqrt(max(sum(real(rmat_o,dp)*real(rmat_o,dp)), DTINY))
            rmat_o = rmat_o / real(onrm)
            ! variance along the pooled direction, propagated from the input spectra: the MAP prior
            ! needs the inputs' units, and the Gram eigenvalue is an overlap, not a variance.
            eig_pooled(q) = sum(V(:,q)*V(:,q)*eigin(:))
        end do
        deallocate(nrm, eigin, G, V, ev)
      contains
        !> leading na_c indices address basis A, the remainder basis B
        subroutine pick( idx, p )
            integer,       intent(in)  :: idx
            real, pointer, intent(out) :: p(:,:,:)
            if( idx <= na_c )then
                call imgs_a(idx)%get_rmat_ptr(p)
            else
                call imgs_b(idx-na_c)%get_rmat_ptr(p)
            endif
        end subroutine pick
    end subroutine bag_basis_pool

    !> Turn a set of real-space basis volumes into embedding-ready column reconstructors.
    subroutine basis_recs_from_images( params, build, imgs, ncomp, basis_recs )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        integer,             intent(in)    :: ncomp
        type(image),         intent(inout) :: imgs(ncomp)
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        integer :: q
        allocate(basis_recs(ncomp))
        do q = 1, ncomp
            call init_column_reconstructor(params, build, basis_recs(q))
            call basis_recs(q)%set_rmat(imgs(q)%get_rmat(), .false.)
            call basis_recs(q)%fft
            call basis_recs(q)%expand_exp
        end do
    end subroutine basis_recs_from_images

    !> Contrast-aware MAP embedding (supplement S.E, eqs S.14-S.15).
    subroutine embed_latents_with_contrast( params, build, mean_rec, basis_recs, ncomp, eigvals, sig2_eff, &
        &pinds, nptcls, z, contrast, precision, resid_energy, resid_mean_energy, rho_out, stats_only, &
        &from_parts )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        integer,             intent(in)    :: ncomp, pinds(:), nptcls
        real(dp),            intent(in)    :: eigvals(ncomp)
        real(dp),            intent(in)    :: sig2_eff
        real(dp),            intent(out)   :: z(nptcls,ncomp), contrast(nptcls)
        real(dp),            intent(out)   :: precision(ncomp,ncomp,nptcls)
        real(dp),            intent(out)   :: resid_energy(nptcls), resid_mean_energy(nptcls)
        ! per-component split-half reliability, exported so the STATE stage can order the latent
        ! path by how well each component is measured rather than by how much variance it carries.
        ! A high-variance, low-rho nuisance component otherwise dominates target placement.
        real(dp), optional,  intent(out)   :: rho_out(ncomp)
        !> worker: run the image pass over THIS part's particles, ship the sufficient statistics, stop
        logical,  optional,  intent(in)    :: stats_only
        !> master: skip the image pass entirely, gather the parts, run the coupled phase
        logical,  optional,  intent(in)    :: from_parts
        real(dp), parameter :: A_LO = 0.1d0, A_HI = 5.0d0
        type(fplane_type), allocatable :: fpls(:)
        type(fplane_type), allocatable :: basis_fpls(:,:), mean_fpl(:), data_fpl(:)
        type(ori), allocatable :: orientations(:)
        real(dp), allocatable :: prior(:), rho(:), Gcache(:,:,:), bcache(:,:), ccache(:,:)
        real(dp), allocatable :: zhalf(:,:,:), Ghf(:,:,:,:), bhf(:,:,:), chf(:,:,:)
        real(dp), allocatable :: Gpart(:,:,:), bpart(:,:), cpart(:,:)
        integer,  allocatable :: prows(:)
        integer :: ipart, pn_part
        real(dp), parameter   :: RHO_FLOOR = 1.d-3
        real(dp) :: rho_max, rrel
        logical :: l_relprior, l_stats_only, l_from_parts, l_polar_embed, l_devprep, l_cache_stats
        integer :: ihf
        integer :: batchlims(2), batchsz, ibatch, i, q, r, ithr, nthr, ia, row
        integer, allocatable :: nzeroG_thr(:), nzeroR_thr(:), nzeroZ_thr(:)   ! dead-basis counters
        real(dp) :: a, a_best, a_keep, a_num, a_den, e_yy, e_mm, best_res, res, aa, sig2
        real(dp), allocatable :: Ainv_th(:,:,:)     ! posterior covariance, for the tr(G A^-1) term
        integer  :: icm, ncm_contrast
        integer(timer_int_kind) :: t_phase
        real(dp), allocatable :: Gth(:,:,:), Ath(:,:,:), zth(:,:), zbest(:,:), cth(:,:), bth(:,:), myth(:)
        real(dp), allocatable :: Gtilth(:,:,:)   ! per-thread noise-whitened projected Gram
        real(dp), allocatable :: gwork(:,:,:), gvec(:,:,:), gev(:,:), gspec_thr(:,:), gspec(:)
        integer,  allocatable :: gcnt_thr(:)
        integer :: nrot_t, gcnt, blk_rp
        real(dp) :: gsum
        nthr = nthr_glob
        sig2 = max(sig2_eff, DTINY)      ! whitened-noise variance for the MAP shrinkage
        allocate(nzeroG_thr(nthr), nzeroR_thr(nthr), nzeroZ_thr(nthr), source=0)   ! dead-basis counters
        ! DEFAULT FLIPPED 2026-08-19: the PLAIN prior wins on AMI at n=3 replications
        ! (+0.024/+0.022/+0.012 across three operating points incl. production nparts=4),
        ! so the reliability prior is now OPT-IN via SIMPLE_COV_RELPRIOR=1. cov_env_int
        ! cannot express "off", hence the explicit >0 read (the recurring env-int trap).
        blk_rp = 0
        call cov_env_int('SIMPLE_COV_RELPRIOR', blk_rp)
        l_relprior = blk_rp > 0
        l_stats_only = .false.
        l_from_parts = .false.
        if( present(stats_only) ) l_stats_only = stats_only
        if( present(from_parts) ) l_from_parts = from_parts
        ! RELPRIOR=0 no longer forces the stage in-process. The distributed flow ships G/b/c
        ! and the split-half solves so the master can re-solve under the reliability-scaled
        ! prior; with the prior PLAIN the same re-solve applies 1/Gamma_q instead, so the
        ! caches are kept and only the rho computation is skipped. Measured motivation: two
        ! nparts=1-matched pairs read +0.037/+0.024 and +0.025/+0.022 (ARI/AMI) for the plain
        ! prior on the EM arm -- the rho^2 rescaling over-shrinks the reproducible directions
        ! 5-8 whose rho sits at 0.46-0.55.
        l_cache_stats = l_relprior .or. l_stats_only .or. l_from_parts
        allocate(prior(ncomp))
        if( l_cache_stats )then
            ! the reducing master never holds Gcache at nptcls: it reads one part's blocks at a time
            ! in the re-solve below
            if( l_from_parts )then
                allocate(Gcache(0,0,0), bcache(0,0), ccache(0,0))
            else
                allocate(Gcache(ncomp,ncomp,nptcls), bcache(ncomp,nptcls), ccache(ncomp,nptcls), source=0.d0)
            endif
            allocate(zhalf(nptcls,ncomp,2), source=0.d0)
            allocate(Ghf(ncomp,ncomp,2,nthr), bhf(ncomp,2,nthr), chf(ncomp,2,nthr), source=0.d0)
        endif
        do q = 1, ncomp
            prior(q) = 1.d0 / max(eigvals(q), DTINY)
        end do
        allocate(Gth(ncomp,ncomp,nthr), Ath(ncomp,ncomp,nthr), zth(ncomp,nthr), zbest(ncomp,nthr), &
            &cth(ncomp,nthr), bth(ncomp,nthr), myth(nthr), Gtilth(ncomp,ncomp,nthr))
        allocate(gwork(ncomp,ncomp,nthr), gvec(ncomp,ncomp,nthr), gev(ncomp,nthr))
        allocate(Ainv_th(ncomp,ncomp,nthr), source=0.d0)
        allocate(gspec_thr(ncomp,nthr), source=0.d0)
        allocate(gcnt_thr(nthr), source=0)
        allocate(basis_fpls(ncomp,nthr), mean_fpl(nthr), data_fpl(nthr), orientations(MAXIMGBATCHSZ))
        call mean_rec%expand_exp
        do q = 1, ncomp
            call basis_recs(q)%expand_exp
        end do
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        z = 0.d0; contrast = 1.d0; precision = 0.d0; resid_energy = 0.d0; resid_mean_energy = 0.d0
        t_phase = tic()
        ! master reducing parts: the workers have already paid the image pass below, so the master
        ! only needs their sufficient statistics to run the coupled phase (rho over every particle,
        ! then the re-solve)
        if( l_from_parts )then
            ! pass 1 only: zhalf and the per-particle scalars are all rho needs, and rho has to exist
            ! before any particle can be solved. The Gram blocks stay on disk until the re-solve
            ! reads them one part at a time.
            call reduce_embed_zhalf_parts(params, flex_pca_nparts(), pinds, contrast, resid_energy, &
                &resid_mean_energy, zhalf, nptcls, ncomp)
            goto 200
        endif
        write(logfhandle,'(A)') '>>> FLEX_PCA CONTRAST-AWARE EMBEDDING'
        call flush(logfhandle)
        ! POLAR PATH: the shared-direction former replaces the per-particle batch loop below. It
        ! fills the same sufficient statistics, so the reliability prior, the re-solve, the precision
        ! matrices and the worker part-write are all untouched. It needs Gcache, which only exists on
        ! the reliability-prior path -- the same constraint the distributed path already carries.
        l_polar_embed = cov_polar_embed_enabled()
        if( l_polar_embed )then
            if( .not. l_relprior ) THROW_HARD('polar embedding requires the reliability prior; &
                &unset SIMPLE_COV_RELPRIOR')
            call embed_accumulate_polar(params, build, mean_rec, basis_recs, ncomp, eigvals, sig2, &
                &pinds, nptcls, Gcache, bcache, ccache, zhalf, contrast, resid_energy, &
                &resid_mean_energy, prior, nzeroG_thr(1))
            nzeroG_thr(1) = 0
        endif
        if( .not. l_polar_embed )then
        call cov_dev_prep_start(params, build, l_devprep)
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz), mskrad=cov_image_mask_radius(params))
            do i = 1, batchsz
                call build%spproj_field%get_ori(pinds(batchlims(1)+i-1), orientations(i))
            end do
            !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
            !$omp& private(i,ithr,q,r,ia,a,a_best,a_keep,a_num,a_den,icm,aa,e_yy,e_mm,best_res,res,row)
            do i = 1, batchsz
                if( orientations(i)%isstatezero() ) cycle
                ithr = omp_get_thread_num() + 1
                row  = batchlims(1) + i - 1
                call project_fplanes_mean_basis(mean_rec, basis_recs, orientations(i), fpls(i), &
                    &mean_fpl(ithr), basis_fpls(:,ithr), apply_ctf_amp=.true.)
                ! data plane = whitened observation (fpls(i)); mean_fpl = T mu ; basis = T U
                e_yy = real(cov_herm_inner(fpls(i), fpls(i)), dp)
                e_mm = real(cov_herm_inner(mean_fpl(ithr), mean_fpl(ithr)), dp)
                myth(ithr) = real(cov_herm_inner(mean_fpl(ithr), fpls(i)), dp)
                do q = 1, ncomp
                    bth(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), fpls(i)), dp)      ! (TU)* y
                    cth(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), mean_fpl(ithr)), dp) ! (TU)* T mu
                    do r = q, ncomp
                        Gth(q,r,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), basis_fpls(r,ithr)), dp)
                        Gth(r,q,ithr) = Gth(q,r,ithr)
                    end do
                end do
                ! split-half sufficient statistics, for the reliability-weighted prior below
                if( l_relprior )then
                    do ihf = 1, 2
                        do q = 1, ncomp
                            bhf(q,ihf,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), fpls(i), ihf), dp)
                            chf(q,ihf,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), mean_fpl(ithr), ihf), dp)
                            do r = q, ncomp
                                Ghf(q,r,ihf,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), basis_fpls(r,ithr), ihf), dp)
                                Ghf(r,q,ihf,ithr) = Ghf(q,r,ihf,ithr)
                            end do
                        end do
                    end do
                endif
                resid_mean_energy(row) = e_yy - 2.d0*myth(ithr) + e_mm                       ! contrast=1 mean residual
                ! Contrast (S.E): for each a on the grid solve the fixed-a MAP
                ! (a^2 G/sig2 + Gamma^-1) z = (a b - a^2 c)/sig2 and keep the a with the lowest residual.
                ! With COV_EMBED_CONTRAST_GRID off this is a single pass at a = 1.
                best_res = huge(1.d0)
                a_best   = 1.d0
                a_keep   = 1.d0     ! a NaN residual would leave this unset and hand garbage downstream
                    aa = a_best*a_best
                    Ath(:,:,ithr) = (aa/sig2)*Gth(:,:,ithr)
                    do q = 1, ncomp
                        Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                        zth(q,ithr)   = (a_best*bth(q,ithr) - aa*cth(q,ithr))/sig2
                    end do
                        call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), ncomp)
                    a = a_best
                    aa  = a*a
                    res = e_yy + aa*e_mm - 2.d0*a*myth(ithr) + quad_form(Gth(:,:,ithr), zth(:,ithr), ncomp)*aa
                    do q = 1, ncomp
                        res = res + 2.d0*aa*zth(q,ithr)*cth(q,ithr) - 2.d0*a*zth(q,ithr)*bth(q,ithr)
                    end do
                    res = res/sig2
                    do q = 1, ncomp
                        res = res + prior(q)*zth(q,ithr)*zth(q,ithr)
                    end do
                    if( res < best_res )then
                        best_res      = res
                        a_keep        = a
                        zbest(:,ithr) = zth(:,ithr)
                    endif
                a_best = a_keep
                ! projected-Gram spectrum on a subsample, for the conditioning report below
                if( mod(row, GRAM_DIAG_STRIDE) == 0 )then
                    gwork(:,:,ithr) = Gth(:,:,ithr)
                    call jacobi(gwork(:,:,ithr), ncomp, ncomp, gev(:,ithr), gvec(:,:,ithr), nrot_t)
                    call eigsrt(gev(:,ithr), gvec(:,:,ithr), ncomp, ncomp)
                    do q = 1, ncomp
                        gspec_thr(q,ithr) = gspec_thr(q,ithr) + max(gev(q,ithr), 0.d0)
                    end do
                    gcnt_thr(ithr) = gcnt_thr(ithr) + 1
                endif
                ! count the particles whose projected basis or rhs came out numerically dead
                if( maxval(abs(Gth(:,:,ithr))) <= 0.d0 ) nzeroG_thr(ithr) = nzeroG_thr(ithr) + 1
                if( maxval(abs(bth(:,ithr) - cth(:,ithr))) <= 0.d0 ) nzeroR_thr(ithr) = nzeroR_thr(ithr) + 1
                if( maxval(abs(zbest(:,ithr))) <= 0.d0 ) nzeroZ_thr(ithr) = nzeroZ_thr(ithr) + 1
                contrast(row)     = a_best
                z(row,:)          = zbest(:,ithr)
                resid_energy(row) = best_res
                aa = contrast(row)*contrast(row)
                Gtilth(:,:,ithr) = (aa/sig2)*Gth(:,:,ithr)
                call map_sampling_precision(Gtilth(:,:,ithr), prior, ncomp, precision(:,:,row))
                if( l_cache_stats )then
                    ! cache the sufficient statistics so the master's re-solve can run in closed
                    ! form with no second pass over the images. Cached whenever stats FLOW
                    ! (reliability prior on, or distributed either side): with RELPRIOR=0 under
                    ! distribution the master re-solves with the PLAIN prior, so skipping the
                    ! cache here would ship all-zero blocks and silently zero every latent.
                    Gcache(:,:,row) = Gth(:,:,ithr)
                    bcache(:,row)   = bth(:,ithr)
                    ccache(:,row)   = cth(:,ithr)
                    ! and the two half-data solves at the chosen contrast
                    do ihf = 1, 2
                        Ath(:,:,ithr) = (aa/sig2)*Ghf(:,:,ihf,ithr)
                        do q = 1, ncomp
                            Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                            zth(q,ithr)   = (contrast(row)*bhf(q,ihf,ithr) - aa*chf(q,ihf,ithr))/sig2
                        end do
                        call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), ncomp)
                        zhalf(row,:,ihf) = zth(:,ithr)
                    end do
                endif
            end do
            !$omp end parallel do
            if( batchlims(2) == nptcls .or. mod(batchlims(2), 5*MAXIMGBATCHSZ) == 0 )then
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA CONTRAST EMBED PARTICLES: ',batchlims(2),' / ',nptcls
                call flush(logfhandle)
            endif
        end do
        call cov_dev_prep_stop(l_devprep)
        endif
        ! the reducing master lands here rather than after the diagnostics, so it still frees the
        ! per-thread Gram workspace; gcnt_thr is all zero when no batch loop ran, so the per-particle
        ! spectrum report below skips itself
200     continue
        allocate(gspec(ncomp), source=0.d0)
        gcnt = sum(gcnt_thr)
        if( gcnt > 0 )then
            do q = 1, ncomp
                gspec(q) = sum(gspec_thr(q,:)) / real(gcnt, dp)
            end do
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA projected-Gram spectrum (mean over ', &
                &gcnt,' particles), largest first:'
            write(logfhandle,'(A,10(1X,ES9.2))') '>>>   ', (gspec(q), q=1,min(10,ncomp))
            write(logfhandle,'(A,ES11.4,A,ES11.4)') '>>>   Gram condition number lam1/lamN = ', &
                &gspec(1)/max(gspec(ncomp),DTINY), '   lam1/lam5 = ', gspec(1)/max(gspec(min(5,ncomp)),DTINY)
            ! the participation ratio is only a rank if the spectrum is normalised first
            gsum = sum(gspec)
            if( gsum > 0.d0 )then
                gspec = gspec / gsum
                write(logfhandle,'(A,F7.3,A,I0)') '>>>   effective rank (participation ratio) = ', &
                    &1.d0 / max(sum(gspec**2), 1.d-300), '  of ', ncomp
            endif
            call flush(logfhandle)
        endif
        deallocate(gwork, gvec, gev, gspec_thr, gcnt_thr, gspec)
        deallocate(Ainv_th)
        ! The distribution of the fitted contrast is the first falsification test for the whole
        ! contrast story, and it should be readable without a script: if SIMPLE's per-particle
        ! normalisation has already absorbed the amplitude then this spread is ~0 and freeing a_i
        ! cannot help anything. The count at each bound says whether the bracket is binding, which
        ! would mean the bracket is doing the fitting rather than the data.
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA D2 zeroG=',sum(nzeroG_thr), &
            &' zeroRHS=',sum(nzeroR_thr),' zeroZ=',sum(nzeroZ_thr),' of nptcls=',nptcls
        write(logfhandle,'(A,F8.1)') '>>> FLEX_PCA CONTRAST EMBED SECONDS: ', toc(t_phase)
        call flush(logfhandle)
        ! worker: the image pass is done for this part's particles and everything below is coupled
        ! across all of them, so ship the sufficient statistics and leave the rest to the master
        if( l_stats_only )then
            call write_embed_stats_part(flex_pca_part_fname('embedstats', params%part, params%numlen), &
                &pinds, contrast, resid_energy, resid_mean_energy, Gcache, bcache, ccache, zhalf, &
                &nptcls, ncomp)
            if( present(rho_out) ) rho_out = 1.d0
            goto 900
        endif
        ! ---- RELIABILITY-WEIGHTED PRIOR ---- The plain prior precision sig2/Gamma_q hands the LARGEST
        ! eigenvalue the WEAKEST prior, so a high-variance but poorly measured component becomes
        ! near-unregularized least squares. Rescale each prior by the component's split-half reliability.
        if( l_cache_stats )then
            if( l_relprior )then
            allocate(rho(ncomp))
            ! optional export of the half-data solves for calibrating the per-particle error model:
            ! rho below collapses the pair to one correlation per component, whereas the variance of
            ! their difference measures the error directly. Off unless SIMPLE_COV_ZHALF is set.
            call write_zhalf_replicates(zhalf, prior, nptcls, ncomp)
            do q = 1, ncomp
                rho(q) = corr_dp(zhalf(:,q,1), zhalf(:,q,2), nptcls)
                rho(q) = max(0.d0, rho(q))
                rho(q) = 2.d0*rho(q) / (1.d0 + rho(q))            ! Spearman-Brown to full length
            end do
            ! Scale rho RELATIVE to the most reliable component, not absolutely: an absolute rho^2 shrinks
            ! the informative components as well and compresses the latent spread state placement needs.
            rho_max = maxval(rho)
            if( rho_max <= DTINY ) rho_max = 1.d0
            do q = 1, ncomp
                rrel = (rho(q)*rho(q)) / (rho_max*rho_max)
                prior(q) = 1.d0 / max(max(rrel, RHO_FLOOR) * eigvals(q), DTINY)
            end do
            ! all components, not the leading 10: rho drives state-target placement via comp_rho, so
            ! its ranking has to be checkable over the whole basis
            write(logfhandle,'(A)') '>>> FLEX_PCA split-half reliability per component (rho, corrected):'
            do q = 1, ncomp
                write(logfhandle,'(A,I3,A,F7.4,A,ES11.3,A,ES11.3)') '>>>   z',q,'  rho=',rho(q), &
                    &'  eigval=',eigvals(q),'  prior_precision=',prior(q)
            end do
            call flush(logfhandle)
            else
                write(logfhandle,'(A)') '>>> FLEX_PCA re-solving with the PLAIN prior &
                    &(RELPRIOR=0, distributed): rho not computed, no rescaling'
                call flush(logfhandle)
            endif
            ! re-solve every particle in closed form from the cached sufficient statistics. Two
            ! routes to the same arithmetic: in process the blocks are already in Gcache, while a
            ! reducing master streams them back one part at a time to keep its footprint flat.
            if( l_from_parts )then
                do ipart = 1, flex_pca_nparts()
                    call read_embed_stats_part(params, ipart, pinds, prows, Gpart, bpart, cpart, &
                        &pn_part, nptcls, ncomp)
                    !$omp parallel do default(shared) private(i,row,q,aa,ithr) schedule(static) proc_bind(close)
                    do i = 1, pn_part
                        ithr = omp_get_thread_num() + 1
                        row  = prows(i)
                        aa   = contrast(row)*contrast(row)
                        Ath(:,:,ithr) = (aa/sig2)*Gpart(:,:,i)
                        do q = 1, ncomp
                            Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                            zth(q,ithr)   = (contrast(row)*bpart(q,i) - aa*cpart(q,i))/sig2
                        end do
                        call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), ncomp)
                        z(row,:) = zth(:,ithr)
                        Gtilth(:,:,ithr) = (aa/sig2)*Gpart(:,:,i)
                        call map_sampling_precision(Gtilth(:,:,ithr), prior, ncomp, precision(:,:,row))
                    end do
                    !$omp end parallel do
                    deallocate(prows, Gpart, bpart, cpart)
                end do
            else
                !$omp parallel do default(shared) private(row,q,aa,ithr) schedule(static) proc_bind(close)
                do row = 1, nptcls
                    ithr = omp_get_thread_num() + 1
                    aa   = contrast(row)*contrast(row)
                    Ath(:,:,ithr) = (aa/sig2)*Gcache(:,:,row)
                    do q = 1, ncomp
                        Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                        zth(q,ithr)   = (contrast(row)*bcache(q,row) - aa*ccache(q,row))/sig2
                    end do
                    call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), ncomp)
                    z(row,:) = zth(:,ithr)
                    Gtilth(:,:,ithr) = (aa/sig2)*Gcache(:,:,row)
                    call map_sampling_precision(Gtilth(:,:,ithr), prior, ncomp, precision(:,:,row))
                end do
                !$omp end parallel do
            endif
            write(logfhandle,'(A)') '>>> FLEX_PCA latents re-solved from the cached statistics'
            call flush(logfhandle)
            if( present(rho_out) )then
                if( l_relprior )then
                    rho_out = rho
                else
                    rho_out = 1.d0
                endif
            endif
            if( allocated(rho) ) deallocate(rho)
            deallocate(Gcache, bcache, ccache, zhalf, Ghf, bhf, chf)
        else
            ! no split-half statistics available; treat every component as equally measured
            if( present(rho_out) ) rho_out = 1.d0
        endif
900     continue
        do i = 1, size(orientations)
            call orientations(i)%kill
        end do
        do ithr = 1, nthr
            call cleanup_plane(mean_fpl(ithr)); call cleanup_plane(data_fpl(ithr))
            do q = 1, ncomp
                call cleanup_plane(basis_fpls(q,ithr))
            end do
        end do
        call cleanup_rec_buffers(build, fpls)
        deallocate(prior, Gth, Ath, zth, zbest, cth, bth, myth, Gtilth, basis_fpls, mean_fpl, data_fpl, orientations)
        deallocate(nzeroG_thr, nzeroR_thr, nzeroZ_thr)
        if( allocated(Gcache) ) deallocate(Gcache, bcache, ccache, zhalf)
        if( allocated(Ghf)    ) deallocate(Ghf, bhf, chf)
        if( allocated(gwork)  ) deallocate(gwork, gvec, gev, gspec_thr, gcnt_thr)
    end subroutine embed_latents_with_contrast

    !> Probe-based subspace iteration: alternate a Wiener E-step (per-particle latents in the current
    !! basis) with a weighted-backprojection M-step (Y_q += sum_i z_iq * backproject(r_i)), then
    !! orthonormalize the refined probe volumes into the next basis.
    subroutine probe_subspace_iteration( params, build, mean_rec, basis_recs, eigvals, sig2_eff, &
        &pinds, nptcls, ncomp, niters )
        use simple_ptcl_cache, only: ptcl_cache_in_use, ptcl_cache_read_batch
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), allocatable, intent(inout) :: basis_recs(:)
        real(dp),            allocatable, intent(inout) :: eigvals(:)
        real(dp),            intent(in)    :: sig2_eff
        integer,             intent(in)    :: pinds(:), nptcls, niters
        integer,             intent(inout) :: ncomp
        type(fplane_type), allocatable :: fpls(:), basis_fpls(:,:), mean_fpl(:)
        type(ori),           allocatable :: orientations(:)
        type(reconstructor), allocatable :: Yeven(:), Yodd(:), utilde(:)
        type(image),         allocatable :: realvols(:), utilde_real(:)
        type(image) :: img_o, mstep_gridcorr
        real,                allocatable :: rho_e(:,:,:,:), rho_o(:,:,:,:), filt(:), corrs(:)
        ! GPU M-step path (SIMPLE_COV_GPU=1): one coupled device accumulation over the combined
        ! [Yeven, Yodd] x [rho_e, rho_o] layout, halfsets routed through the scale slots
        logical  :: l_gpu_probe
        integer  :: vgpu, gq, gr
        real(dp), allocatable :: gdsc(:,:), grsc(:,:)
        ! banked coupled adjoint (SIMPLE_COV_GPU_BANKADJ, default on with the GPU M-step):
        ! probe records sorted by bank direction, in-plane aligned on device, GEMM-combined per
        ! direction, one splat per (direction, slot) instead of per particle
        logical  :: l_bank_adj
        integer  :: vbank, ndir_pb, jdir
        type(oris) :: dirs_pb
        type(ori)  :: o_pb
        real     :: ca1, sa1
        real,    allocatable :: rmat_pb(:,:,:), nrm_pb(:,:), rmat_pp(:,:,:), nrm_pp(:,:)
        real,    allocatable :: cap(:), sap(:)
        integer, allocatable :: dirof_pb(:), cnt_pb(:), perm_pb(:), pptmp(:)
        real(dp),            allocatable :: zbatch(:,:), dens(:,:,:), z(:,:), prior(:)
        real(dp),            allocatable :: Gth(:,:,:), Ath(:,:,:), bth(:,:), cth(:,:), zth(:,:)
        !> per-thread posterior covariance A^-1 and the untouched copy of A it is taken from
        real(dp),            allocatable :: Ainvth(:,:,:), Acpth(:,:,:)
        !> marginal-likelihood scratch: h_i before the solve, and the per-thread accumulator of
        !! log det A_i - h_i' A_i^-1 h_i. The remaining terms of -2 log p(y_i) are the per-particle
        !! data energy (fixed across iterations, since a_i and mu are) and N*log det Gamma, added at
        !! the reduction. See spd_logdet_dp for why resid_energy could never serve this purpose.
        real(dp),            allocatable :: hth(:,:), nll_thr(:), gam_dbg(:,:)
        real(dp) :: ldA, nll_tot, nll_prev
        !> principal-angle convergence threshold. Overridable because the default 0.97 was tuned for a
        !! probe that STARTED from the covariance eigenbasis and only had to clean it: from a data-free
        !! start the subspace stops rotating several iterations before Gamma stops moving (measured on
        !! 10076: cos crosses 0.97 at iteration 5 while max var is still falling ~20 % per iteration),
        !! so a pure-EM fit needs a tighter bar. SIMPLE_COV_PROBE_CONV is in PER MILLE (995 = 0.995).
        real(dp) :: conv_thresh
        integer  :: vconv
        !> even/odd half-basis agreement: the dataset-agnostic convergence signal
        type(image), allocatable :: eimgs(:), oimgs(:)
        real(dp),    allocatable :: sv_eo(:)
        real(dp) :: eo_dim, eo_best
        integer  :: eo_stall, vpat, eo_patience
        !> BAND ANNEALING, SIMPLE_COV_EM_ANNEAL = number of iterations to reach the full band.
        !! Expressed as a k_hi schedule and converted to lp, NOT set through lp directly: on this
        !! geometry covariance_kfromto only narrows the band when lp > 2*smpd_crop, so a schedule
        !! written in lp is silently a no-op at common settings. Rationale is that a from-scratch
        !! basis has nothing to constrain the high shells in the first iterations, so letting them
        !! in early spends rank on noise. NOTE: this schedule is the proposal's own invention -- it
        !! is not from 3DVA, which uses a single fixed filter -- so it is OFF by default and is a
        !! thing to A/B, not a thing to assume.
        integer  :: vanneal, khi_full, khi_it, kfr_ann(2)
        integer, parameter :: MIX_ZSUB_MAX = 2000   ! per-part latent subsample shipped for mcfa_init
        real(dp), allocatable :: dm_sr(:), dm_sm(:,:), dm_smm(:,:,:), dm_sai(:,:), dm_z(:,:)
        integer :: dm_nz = 0
        real     :: lp_it, dstep_ann
        !> POLAR EM E-STEP, SIMPLE_COV_EM_POLAR=1.
        !! probe_subspace_iteration has never had a polar former: its three E-step branches are all
        !! Cartesian (GPU-fused exact, banked-forward, CPU exact), and cov_polar_enabled() is only
        !! consulted by reduced_covariance_solve and embed_latents_with_contrast. This routes the
        !! per-particle sufficient statistics through embed_accumulate_polar -- the SAME polar
        !! former the embedding uses, which fills exactly Gcache/bcache/ccache/contrast -- and then
        !! keeps the Cartesian residual pass for the M-step, because there is no polar->volume
        !! adjoint anywhere in the tree and there is no prospect of one.
        !! The residual pass projects the MEAN ONLY, not the r+1 volumes, so at large r this is
        !! cheaper than the Cartesian E-step it replaces as well as being the right representation.
        logical  :: l_em_polar
        integer  :: vempol, nvalid_pol
        real(dp), allocatable :: pGc(:,:,:), pbc(:,:), pcc(:,:), pzh(:,:,:)
        real(dp), allocatable :: pcon(:), pre(:), prme(:)
        !> POLAR SHARED-DIRECTION E-STEP, SIMPLE_COV_POLAR_ESTEP=1 (stage 1).
        !! IN the batch loop, unlike SIMPLE_COV_EM_POLAR (which routes the whole statistics pass
        !! through embed_accumulate_polar up front and is windows/mixture-incompatible): once per
        !! EM iteration the mean + ncomp basis volumes are projected onto the polar rings of a
        !! shared direction table (the bank); per particle, the bank rings at (direction, psi)
        !! replace the ncomp+1 per-particle Cartesian projections in forming G/b/c/e_mm/myv, and
        !! everything downstream -- the ECM/MCFA solve, mixture accumulators, ONLINE windows,
        !! e/o halfsets, the Cartesian M-step insertion -- is untouched. CTF amplitude, shift
        !! phase and per-shell whitening ride in exactly as in the Cartesian former, because the
        !! SAME prepped cmplx/transfer planes are polar-sampled (the validated embed formulation).
        !! The mean-shaped deflation needs no special handling here: dfl_basis deflates the
        !! refined basis VOLUMES at the end of each M-step, and the bank is rebuilt from those
        !! same basis_recs at the next iteration, so the deflation enters the bank exactly as it
        !! enters the per-particle Cartesian projections.
        logical  :: l_pol_es, l_pol_grid, l_pol_bank_it
        integer  :: vpoles, n_pol_check, ndir_es, nsamp_es, nsamp2_es, nk_es, id_es, idp_es, ir
        integer  :: ph0_es, pk0_es, hlo_es, hhi_es, klo_es, nyqr_es, nyqb_es
        !> HYBRID exact/ring quadrature (accuracy): the shells 0..rhyb_es are accumulated as
        !! EXACT Cartesian lattice statistics per particle (data, CTF and model read exactly as
        !! the Cartesian former reads them -- including the DC sample the rings never had) and
        !! the shared-direction rings cover only the annulus above. Measured on the calibrated
        !! synthetic dissector: the pure ring quadrature carries a MULTIPLICATIVE low bias of
        !! the posterior latents (z-variance ratio ~0.85 vs Cartesian, the eigen-spectrum
        !! collapse mechanism of the first science A/B), which angular/radial oversampling and
        !! the DC sample alone do NOT fix; the hybrid at rhyb ~ 0.55*band restores ratio ~1.00
        !! with G err 0.4% and b err ~3%. SIMPLE_COV_POLAR_RHYB overrides (0=off: pure rings,
        !! tonight's baseline); SIMPLE_COV_POLAR_OSAMP multiplies ring angular sampling.
        logical  :: l_pol_hyb, l_rhyb_off
        integer  :: vpol_rhyb, vpol_osamp, rhyb_es, npos_es, jx_es, kx_es
        integer,  allocatable :: hex_es(:), kex_es(:)
        type(polar_grid_t)    :: pg_es
        type(oris)            :: dirs_es
        type(ori)             :: o_es
        real,     allocatable :: rmatp_es(:,:,:), nrmp_es(:,:), rmatb_es(:,:,:), nrmb_es(:,:)
        real,     allocatable :: cae(:), sae(:)
        integer,  allocatable :: dir_es(:)
        logical,  allocatable :: dused_es(:)
        real,     allocatable :: UsallE(:,:,:)                       !< (nsamp2, 0:ncomp, ndir) the bank
        real(dp), allocatable :: CfE(:,:,:), Cm0E(:,:,:), c00E(:,:)  !< ring Gram tables per direction
        complex,  allocatable :: UbankE(:,:,:)                       !< per-thread projection scratch
        real,     allocatable :: CspE(:,:,:)                         !< per-thread ring-Gram scratch
        real,     allocatable :: xws_es(:,:), wr_es(:,:), Reb_es(:,:)
        real(dp), allocatable :: wrd_es(:,:)
        real(dp) :: pw_es, cnt_es
        real     :: taz_es, sec_bank
        integer(timer_int_kind) :: t_bank
        !> DEVICE polar E-step (stage 2), SIMPLE_COV_POLAR_GPU=1 (requires POLAR_ESTEP=1): the
        !! host-built bank + ring Gram tables upload once per iteration; per batch the RAW
        !! (y, T) plane sections go up and per-particle G/b/c/e_mm/myv come back in the fused
        !! E-step's stat layout -- the ring part via a device sampler + bank/table contractions,
        !! the hybrid low-k exact part via the unmodified estep_stats_kernel with the disc
        !! capped at rhyb*(rhyb+1) against the resident mean+basis volumes. The solve (plain or
        !! mixture), Gamma/nll accumulation, ONLINE windows, e/o bookkeeping, the Cartesian
        !! mean projection/residual and the M-step insertion are the CPU polar path's,
        !! untouched. All device state is per-stage (ring geometry), per-iteration (bank,
        !! volumes) or per-batch (planes), so the path is window-safe under SIMPLE_COV_ONLINE
        !! by construction (the ONLINE x fused-device breakage was in the prep/banked caches,
        !! which this path does not use).
        logical  :: l_pol_gpu, l_pol_gpu_any
        integer  :: vpolg, nst_pg
        real(dp), allocatable :: stats_pg(:,:)
        real,     allocatable :: ca_pg(:), sa_pg(:)
        integer,  allocatable :: dir_pg(:)
        logical,  allocatable :: vld_pg(:)
        !> gate-1 GPU parity arrays: the CPU polar former re-run on the check particles
        real(dp), allocatable :: zchk_g(:,:), gerrg_chk(:), berrg_chk(:)
        real(dp), allocatable :: Gpc_chk(:,:,:), bpc_chk(:,:), cpc_chk(:,:), zpc_chk(:,:)
        real(dp) :: emmg_chk, myvg_chk, ag_chk
        !> gate-1 parity diagnostic, SIMPLE_COV_POLAR_CHECK=N: the first N particles of EM
        !! iteration 1 run BOTH formers; per-component z correlation and the median relative
        !! G/b errors are printed, then the run continues on the polar path.
        logical  :: l_chk_act, l_chk_p, lok_chk
        integer  :: nchk_max, ichk
        real(dp), allocatable :: zchk_c(:,:), zchk_p(:,:), gerr_chk(:), berr_chk(:)
        logical,  allocatable :: lchk(:)
        real(dp), allocatable :: Gc_chk(:,:,:), bc_chk(:,:), cc_chk(:,:), zc_chk(:,:), Ai_chk(:,:,:)
        real(dp) :: emmc_chk, myvc_chk, ac_chk, lda_chk, qml_chk
        real(dp) :: gnum_chk, gden_chk, bnum_chk, bden_chk
        !> exact-direction polar arm of the check: the SAME quadrature at the particle's own
        !! direction (no snap), which attributes any disagreement between direction quantization
        !! (banked vs exact) and the polar formulation itself (exact vs Cartesian). tazim -- the
        !! mean within-ring relative spread of |T|^2 -- bounds the radial-factorisation error in G.
        real(dp), allocatable :: zchk_x(:,:), gerrx_chk(:), berrx_chk(:)
        real,     allocatable :: UsX_chk(:,:,:), xwsX_chk(:,:), wrX_chk(:,:), RebX_chk(:,:)
        real(dp), allocatable :: CfX_chk(:,:,:), Cm0X_chk(:,:,:), c00X_chk(:,:), wrdX_chk(:,:)
        real(dp), allocatable :: GX_chk(:,:,:), bX_chk(:,:), cX_chk(:,:), zx_chk(:,:)
        real(dp), allocatable :: tazsum_thr(:)
        integer,  allocatable :: tazcnt_thr(:)
        real     :: rmx_chk(3,3), tazx_chk
        real(dp) :: pwx_chk, cntx_chk, emmx_chk, myvx_chk, ax_chk
        !> DC-less Cartesian arm of the check: cov_herm_inner includes the (0,0) sample, the polar
        !! quadrature starts at ring 1 and has none. If the exact-direction polar arm agrees with
        !! THIS reference, the polar/Cartesian gap is the DC term, not the ring quadrature.
        real(dp), allocatable :: zchk_d(:,:), gerrd_chk(:), berrd_chk(:)
        real(dp), allocatable :: Gd_chk(:,:,:), bd_chk(:,:), cd_chk(:,:), zd_chk(:,:)
        complex,  allocatable :: dcb_chk(:,:)
        complex  :: dcy_chk, dcm_chk
        real(dp) :: emmd_chk, myvd_chk, ad_chk
        logical  :: l_probe_distr_pre
        !> mean-shaped (contrast) deflation of the refined basis, SIMPLE_COV_EM_DEFLATE
        logical       :: l_deflate_mean
        integer       :: vdfl, ndfl, ndfl_sh, idfl, jdfl, nkeep_dfl, kfr_dfl(2)
        logical       :: l_dfl_bg
        logical       :: l_dfl_pose
        integer       :: ipdfl, ixp, iyp, izp
        real(dp)      :: gpx, gpy, gpz, cp_dfl
        type(image)   :: mvol_dfl
        type(image), allocatable :: dfl_basis(:)
        real, pointer :: rm_dfl(:,:,:), rv_dfl(:,:,:)
        real          :: res_lo, res_hi
        real(dp)      :: mm_dfl, mv_dfl, rem_dfl, tot_dfl, mnorm_dfl
        logical  :: lok
        !> per-thread accumulator for the EM Gamma update, and its reduction
        real(dp),            allocatable :: gam_thr(:,:), gam_acc(:)
        integer,             allocatable :: nval_thr(:)
        logical,             allocatable :: valid(:), valid_e(:), valid_o(:)
        integer,             allocatable :: eo(:)
        real, pointer :: rmatp(:,:,:)
        real     :: fc
        real(dp) :: sig2, a, aa, e_mm, myv, mu_q, sd_q
        integer  :: it, q, r, i, ithr, nthr, batchlims(2), batchsz, ibatch, row, d_new, es(3), filtsz, sh
        integer  :: npairs, nval, npp, nparts_sub
        logical  :: l_probe_distr
        real(dp),            allocatable :: gam_sum(:)
        complex,             allocatable :: cme(:,:,:,:), cmo(:,:,:,:)
        real,                allocatable :: rhe(:,:,:,:), rhoo(:,:,:,:)
        integer,             allocatable :: ppinds(:)
        type(image),         allocatable :: prev_real(:)      ! previous iteration's orthonormal basis
        real(dp),            allocatable :: Mconv(:,:), sconv(:)
        real(dp) :: cos_mean
        logical  :: l_converged
        integer  :: n_probe_cm, nml_plain
        real,     allocatable :: fscq_dg(:,:)
        real(dp) :: fmean_dg(512), fbest_dg
        real     :: res_dg
        integer  :: khi_dg, nsig_dg, ntop_dg, sel_dg(4), tq_dg, bq_dg
        ! ---- ONLINE (minibatch) EM: SIMPLE_COV_ONLINE=1 + SIMPLE_COV_BATCH=B ----
        ! Each iteration E-steps ONE rotating contiguous window of the (full) probe set and the
        ! M-step consumes RUNNING-AVERAGED sufficient statistics S_t = (1-g_t) R(S_{t-1}) + g_t s_t
        ! (Cappe-Moulines), with R the Procrustes rotation between successive basis frames. The
        ! Wiener merge therefore sees an effective sample that GROWS with iterations while each
        ! iteration costs one batch -- "the whole dataset informs the basis without increasing
        ! per-step compute". Naive per-iteration resampling without averaging is measured negative;
        ! full-N at fixed regularisation is measured negative (subsampling is regularisation) --
        ! this is the formulation that addresses both. v1: plain CPU body, nparts=1, contiguous
        ! project-order windows (sequential IO); rho/gamma rotate with R^2 weights (exact only for
        ! the map stats; second-order error near convergence, diluted by g early).
        logical  :: l_online, l_focus_annulus
        integer  :: vonl_pr, nb_online, ob_lo, ob_hi, b_online, b_win, it_win, vfoc_r1, vfoc_ed
        real     :: foc_r1, foc_edge
        real(dp) :: gam_onl
        complex, allocatable :: oh_cme(:,:,:,:), oh_cmo(:,:,:,:)
        real,    allocatable :: oh_rhe(:,:,:,:), oh_rho(:,:,:,:)
        real(dp), allocatable :: oh_gam(:), ut_hist(:,:), r_onl(:,:)
        logical  :: l_oh_live
        logical  :: l_probe_mls
        real(dp) :: qml, nll_mix_add
        ! ---- MCFA state: tied-covariance mixture prior over the latents ----
        integer  :: kmix_pr, n_mix_warm, kk2
        logical  :: l_mix_req, l_mix_active, l_mix_used
        real(dp) :: ldOm_mix, ldOm_used, lwm, wsm
        real(dp), allocatable :: mix_xi(:,:), mix_Om(:,:), mix_Ominv(:,:), mix_pi(:)
        real(dp), allocatable :: mix_Omxi(:,:), mix_xiOx(:), mix_lpi(:)
        real(dp), allocatable :: rhs0th(:,:), mkth(:,:,:), lwth(:,:), rkth(:,:)
        real(dp), allocatable :: mxa_sr(:,:), mxa_sm(:,:,:), mxa_smm(:,:,:,:), mxa_sainv(:,:,:)
        ! unified-mixture state: frame rotation (Gap B), full-N running average of the
        ! reduced mixture statistics (Gap A), starved-component reseeding, final full pass
        logical  :: l_onl_final
        type(string) :: fname
        integer(timer_int_kind) :: t_it, t_sec
        real    :: sec_read, sec_prep, sec_estep, sec_ins
        real(dp) :: twp0, twp1, twp2
        real(dp), allocatable :: sec_proj_thr(:), sec_gram_thr(:)
        ! banked FORWARD: shared raw projections per direction segment; data aligned into the
        ! bank frame once (which also hands the M-step a pre-aligned residual).
        ! MEASURED REFUTED 2026-08-14 (gpu15): science ARI 0.836 vs the 0.95 band -- the E-step
        ! latents cannot take direction-quantized projections (the M-step adjoint can) -- and
        ! the full-array CTF multiplies are memory-traffic-bound (gram bucket 30->200 thread-s).
        ! Kept opt-in behind SIMPLE_COV_GPU_BANKFWD for reference; default is the exact
        ! per-particle forward with the banked adjoint.
        type(fplane_type), allocatable :: mean_fplC(:), basis_fplsC(:,:)
        type(ori),         allocatable :: o_pb_thr(:)
        complex,           allocatable :: tmpc(:,:,:)
        real,              allocatable :: cap1(:), sap0(:)
        integer,           allocatable :: seg_dir_b(:), seg_beg_b(:), seg_cnt_b(:)
        integer  :: iseg, ii, nseg_b, nyqb
        logical  :: l_bank_fwd
        ! fused device E-step (exact per-particle directions, volumes resident)
        logical  :: l_gpu_estep, l_deflate
        integer  :: ves, nst_es, off_es
        real(dp), allocatable :: stats_es(:,:)
        real,     allocatable :: avec_es(:)
        logical,  allocatable :: vld_es(:)
        ! device prep (raw images up; packed planes born on device)
        logical  :: l_gpu_prep, l_pcache
        integer  :: vprep, frlims_pb(3,2), nyqpd_pb, kf_pb(2)
        type(ftiter) :: fit_pb
        type(ctfparams), allocatable :: ctfp_arr(:)
        real,            allocatable :: shf_pb(:,:), sig2_pb(:,:)
        integer :: signyq_pb, nyqfull_pb
        nthr = nthr_glob
        sig2 = max(sig2_eff, DTINY)
        ! ---- optional STRIDE subsample, for the basis refinement only ----
        ! The probe refines ncomp band-limited, FSC-regularised volumes, and every iteration costs a
        ! full pass over the data -- far more particles than that many parameters need. A stride keeps
        ! both halfsets and every state proportionally represented; fromp/top would NOT, because
        ! particles are commonly ordered by state (on Ribosembly a contiguous window selects whole
        ! states). The embedding stage that follows still uses every particle: only the basis
        ! refinement is subsampled.
        ! The stride MUST be applied within each halfset, not across the particle list. `eo` alternates
        ! strictly by particle index (0,1,0,1,...), so a plain stride of 2 selects one halfset entirely
        ! and leaves the other empty -- and every probe M-step is regularised by an even/odd FSC, which
        ! is then computed against nothing. Measured: the Wiener filter kills the basis and the run dies
        ! at the "embedding collapsed" guard. Striding per halfset keeps both populated at any stride.
        ! ---- absolute cap, not a fixed ratio ----
        ! The probe refines ncomp band-limited, FSC-regularised volumes, and that parameter count
        ! does not grow with the dataset, so the particles needed to determine it do not either. A
        ! constant stride would leave the probe scaling linearly and dominating the run; capping the
        ! count makes it O(1) in dataset size.
        !
        ! COV_PROBE_MAX_PTCLS is the total across all processes, so each takes its share -- a worker
        ! sees only its own partition and would otherwise take the whole budget nparts times over.
        ! SIMPLE_COV_PROBE_STRIDE still wins if set.
        ! Only a WORKER divides the total by nparts: it holds one fromp/top partition. The master
        ! holds every particle, so passing nparts there divides twice and inflates the stride by
        ! exactly nparts. See cov_stage_subsample, which the initialiser shares.
        nparts_sub = 1
        if( flex_pca_is_worker() ) nparts_sub = params%nparts
        call cov_stage_subsample(build, pinds, nptcls, nparts_sub, COV_PROBE_MAX_PTCLS, &
            &'SIMPLE_COV_PROBE_STRIDE', 'SIMPLE_COV_PROBE_MAX', 'PROBE', ppinds, npp)
        allocate(z(npp,ncomp))
        ! ---- banked-adjoint setup: pose geometry is fixed across EM iterations, so the bank,
        ! the direction assignment and the direction sort are built ONCE per stage. Sorting
        ! The polar and mixture E-step selectors are read BEFORE the GPU E-step gating below,
        ! which has to stand those paths down. The EM_POLAR read used to sit a hundred lines
        ! AFTER the gate that tests it, so the gate read an UNINITIALISED logical -- usually
        ! false on this stack, which is why it never bit in practice, but it was undefined
        ! behaviour standing in for "off".
        vempol     = 0
        call cov_env_int('SIMPLE_COV_EM_POLAR', vempol)
        l_em_polar = vempol > 0
        ! ---- POLAR E-STEP (SIMPLE_COV_POLAR_ESTEP=1): in-batch-loop shared-direction bank ----
        vpoles   = 0
        call cov_env_int('SIMPLE_COV_POLAR_ESTEP', vpoles)
        l_pol_es = vpoles > 0
        if( l_pol_es .and. l_em_polar ) THROW_HARD('SIMPLE_COV_POLAR_ESTEP + SIMPLE_COV_EM_POLAR &
            &are mutually exclusive: both replace the E-step former')
        n_pol_check = 0
        call cov_env_int('SIMPLE_COV_POLAR_CHECK', n_pol_check)
        if( n_pol_check > 0 .and. .not. l_pol_es )then
            write(logfhandle,'(A)') '>>> FLEX_PCA POLAR CHECK ignored: needs SIMPLE_COV_POLAR_ESTEP=1'
            n_pol_check = 0
        endif
        ! hybrid/oversampling accuracy knobs (see the declaration block for the measurements)
        vpol_rhyb  = 0
        call cov_env_int('SIMPLE_COV_POLAR_RHYB', vpol_rhyb)
        l_rhyb_off = cov_env_int_off('SIMPLE_COV_POLAR_RHYB')
        vpol_osamp = 1
        call cov_env_int('SIMPLE_COV_POLAR_OSAMP', vpol_osamp)
        vpol_osamp = max(1, min(8, vpol_osamp))
        l_pol_hyb  = .false.
        rhyb_es    = 0
        npos_es    = 0
        l_pol_grid    = .false.
        l_pol_bank_it = .false.
        l_chk_act     = .false.
        nchk_max      = 0
        sec_bank      = 0.
        if( l_pol_es )then
            write(logfhandle,'(A)') '>>> FLEX_PCA POLAR E-STEP ON (stage 1): shared-direction bank &
                &supplies G/b/c; data-plane prep and M-step insertion stay Cartesian'
            if( n_pol_check > 0 ) write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA POLAR CHECK: dual &
                &former on the first ',n_pol_check,' particles of iteration 1'
            call flush(logfhandle)
        endif
        ! ---- DEVICE polar E-step (stage 2), SIMPLE_COV_POLAR_GPU=1: opt-in, composes with the
        ! production env set; falls back to the CPU polar path (same answer) when the card or
        ! the CUDA build is absent, or when the basis exceeds the resident-volume cap below.
        l_pol_gpu = .false.
        if( l_pol_es )then
            vpolg = 0
            call cov_env_int('SIMPLE_COV_POLAR_GPU', vpolg)
            if( vpolg > 0 )then
                if( flex_gpu_available() )then
                    l_pol_gpu = .true.
                    write(logfhandle,'(A)') '>>> FLEX_PCA POLAR E-STEP DEVICE ON (stage 2): bank + &
                        &ring tables resident per iteration, per-batch ring/hybrid statistics on &
                        &device; solve, windows, mixture and M-step unchanged'
                else
                    write(logfhandle,'(A)') '>>> FLEX_PCA POLAR GPU requested but no device/CUDA &
                        &build available: using the CPU polar E-step'
                endif
                call flush(logfhandle)
            endif
        endif
        l_pol_gpu_any = l_pol_gpu
        ! ---- MCFA (SIMPLE_COV_EM_MIX=K): tied-covariance K-component mixture latent prior ----
        ! Mixtures of common factor analysers (Baek, McLachlan & Flack, IEEE TPAMI 32:1298,
        ! 2010): ONE shared basis, K Gaussian components in the latent space, fitted in the same
        ! EM as the basis itself. The structural point: a direction earns rank in U only if it
        ! improves the fit under a MULTI-MODAL prior, so purely continuous nuisance variation --
        ! fit equally well by any single component -- gains nothing. This is the one thing a
        ! moment estimator cannot reproduce, because responsibilities are posterior objects.
        ! K=1 pins the component at the origin with a diagonal covariance, which reduces EXACTLY
        ! to the plain PPCA EM -- that is the mandated regression test, not a feature.
        kmix_pr = 0
        call cov_env_int('SIMPLE_COV_EM_MIX', kmix_pr)
        kmix_pr   = max(0, min(64, kmix_pr))
        l_mix_req = kmix_pr >= 1
        ! v2 (2026-08-21): the four MCFA accumulators are additive sufficient statistics and now
        ! ride in the probe part files (PROBE_PART_VERSION 3), together with a bounded latent
        ! subsample so the master can seed mcfa_init. nparts>1 is therefore supported.
        if( l_mix_req .and. l_em_polar ) THROW_HARD('SIMPLE_COV_EM_MIX + SIMPLE_COV_EM_POLAR &
            &are mutually exclusive')
        vonl_pr = 0
        call cov_env_int('SIMPLE_COV_ONLINE', vonl_pr)
        l_online = vonl_pr > 0
        vfoc_r1 = 0
        call cov_env_int('SIMPLE_COV_FOCUS_R1', vfoc_r1)
        l_focus_annulus = vfoc_r1 > 0
        foc_r1 = real(vfoc_r1)
        vfoc_ed = 10
        call cov_env_int('SIMPLE_COV_FOCUS_EDGE', vfoc_ed)
        foc_edge = real(max(1, vfoc_ed))
        if( l_focus_annulus )then
            write(logfhandle,'(A,F6.1,A,F5.1,A)') '>>> FLEX_PCA FOCUSED HETEROGENEITY: basis &
                &restricted to the shell beyond r=', foc_r1, ' A (soft edge ', foc_edge, ' A)'
            call flush(logfhandle)
        endif
        l_oh_live = .false.
        b_online = 40000
        call cov_env_int('SIMPLE_COV_BATCH', b_online)
        ! Under windows-only ONLINE the delivered basis is the LAST window's M-step --
        ! and which window that is falls out of mod(iters, nwindows). ONLINE_FINAL runs
        ! the final iteration over the FULL probe set, so the endpoint estimator is
        ! built from every particle with no history blending involved (the
        ! full-N-under-mixture-prior configuration, measured good at the ceiling).
        l_onl_final = l_online .and. cov_env_int_on('SIMPLE_COV_ONLINE_FINAL')
        ! assigned BEFORE its first read: this guard used to read l_probe_distr_pre ~200 lines
        ! before the assignment below -- undefined behaviour standing in for "off", the same trap
        ! the EM_POLAR selector note above documents
        l_probe_distr_pre = flex_pca_is_master() .and. flex_pca_nparts() > 1
        if( l_online )then
            ! v2 (2026-08-21): ONLINE distributes. The recurrence itself -- the running average
            ! S_t = (1-g_t) R(S_{t-1}) + g_t s_t and its Procrustes frame rotation -- is MASTER-side
            ! and operates on the REDUCED statistics, so workers stay stateless: each contributes
            ! s_t for its slice of the current window. The only thing that has to be made
            ! distribution-aware is the window walk itself (below): each rank walks its own shard
            ! with a per-rank window of b_online/nparts, so the UNION across ranks is one global
            ! window of the intended size b_online, evenly spread and load-balanced. At nparts=1
            ! this reduces identically to the shared-memory walk.
        endif
        n_mix_warm = 3
        call cov_env_int('SIMPLE_COV_EM_MIX_WARM', n_mix_warm)
        n_mix_warm   = max(1, min(20, n_mix_warm))
        l_mix_active = .false.
        l_mix_used   = .false.
        ldOm_mix     = 0.d0
        ldOm_used    = 0.d0
        if( l_mix_req ) write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA MIX requested: K=', &
            &kmix_pr,'  (plain warm-up: ',n_mix_warm,' iterations)'
        ! ppinds is safe: z is only ever consumed through order-insensitive sums, and z/ppinds
        ! stay aligned. Gated to nsym==1 (the bank factorization R = R_dir * Rz(psi) has no
        ! symmetry replication) and to the in-process E-step (the distributed master never runs
        ! the batch loop).
        l_bank_adj = .false.
        l_bank_fwd = .false.
        l_gpu_estep = .false.
        ! initialized HERE, not inside the SIMPLE_COV_GPU block below: the batch loop reads it
        ! on the CPU path too (skips the CPU prep and touches unallocated ctfp_arr if garbage-true)
        l_gpu_prep = .false.
        vgpu = 0;  call cov_env_int('SIMPLE_COV_GPU', vgpu)
        vbank = 1; call cov_env_int('SIMPLE_COV_GPU_BANKADJ', vbank)
        if( cov_env_int_off('SIMPLE_COV_GPU_BANKADJ') ) vbank = 0
        if( vgpu > 0 .and. vbank > 0 .and. flex_gpu_available() .and. &
            &build%pgrpsyms%get_nsym() == 1 .and. &
            &.not. (flex_pca_is_master() .and. flex_pca_nparts() > 1) )then
            l_bank_adj = .true.
            ndir_pb = cov_polar_ndir(npp)
            call dirs_pb%new(ndir_pb, is_ptcl=.false.)
            call build%pgrpsyms%build_refspiral(dirs_pb)
            allocate(rmat_pb(3,3,ndir_pb), nrm_pb(3,ndir_pb))
            do jdir = 1, ndir_pb
                rmat_pb(:,:,jdir) = dirs_pb%get_mat(jdir)
                nrm_pb(:,jdir)    = rmat_pb(3,:,jdir)
            end do
            ! dirs_pb stays alive: the banked forward projects the shared basis at bank
            ! orientations pulled from it per segment
            vbank = 0
            call cov_env_int('SIMPLE_COV_GPU_BANKFWD', vbank)
            l_bank_fwd = vbank > 0
            ves = 1
            call cov_env_int('SIMPLE_COV_GPU_ESTEP', ves)
            ! cov_env_int only assigns for values > 0, so it cannot express "off" and =0 was
            ! silently ignored -- which voided the first MCFA reduction test: the plain arm ran
            ! the FUSED body while the mixture arm ran the CPU body, and the pair compared
            ! former numerics instead of the mixture algebra. Same trap, same fix, as
            ! SIMPLE_COV_EM_DEFLATE=0.
            l_gpu_estep = ves > 0 .and. .not. cov_env_int_off('SIMPLE_COV_GPU_ESTEP')
            ! the polar former supplies G, b, c and the contrast itself, so the device E-step and
            ! the banked forward must both stand down or their statistics would be the ones used
            if( l_em_polar .or. l_pol_es )then
                l_gpu_estep = .false.
                l_bank_fwd  = .false.
            endif
            ! the banked forward still lacks a mixture solve; the fused device E-step COMPOSES
            ! with the mixture since 2026-08-19 -- the device fetches the same G/b/c sufficient
            ! statistics and the host branches to probe_solve_mix (one shared solver, both bodies)
            if( l_mix_req )then
                l_bank_fwd  = .false.
            endif
            ! CAPACITY GATE. The device holds mean+basis as ONE resident allocation and the CUDA
            ! entry point refuses ncomp+1 > COV_GPU_ESTEP_MAXVOLS. Without this check the worker
            ! dies inside flex_gpu_estep_vols_f (THROW_HARD) and, because a dead worker never writes
            ! JOB_FINISHED, the distributed MASTER then waits forever -- an 8-hour silent hang
            ! observed at neigs=40. Fall back to the CPU E-step instead: slower, same answer.
            if( l_gpu_estep .and. (ncomp + 1) > COV_GPU_ESTEP_MAXVOLS )then
                l_gpu_estep = .false.
                write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA PROBE fused device E-STEP OFF: needs ', &
                    &ncomp+1,' resident volumes, device entry point caps at ',COV_GPU_ESTEP_MAXVOLS, &
                    &' -- using the CPU E-step for this basis size'
                call flush(logfhandle)
            endif
            ! ONLINE x device is MEASURED BROKEN (2026-08-19, onlk8_dev: ceiling 0.2087 vs
            ! 0.2742 CPU at the identical config -- the representation itself collapses).
            ! Some device-side state IS keyed to a fixed particle set across iterations
            ! (candidate: the banked direction segments / prep caches; the streaming argument
            ! that justified removing this guard was wrong). CPU body under windows until the
            ! device path is made window-safe and re-validated.
            if( l_online .and. (l_gpu_estep .or. l_bank_adj) )then
                l_gpu_estep = .false.
                l_bank_adj  = .false.
                write(logfhandle,'(A)') '>>> FLEX_PCA ONLINE: plain CPU E-step body (device path &
                    &measured broken under rotating windows; see 2026-08-19 note)'
                call flush(logfhandle)
            endif
            if( l_gpu_estep )then
                allocate(avec_es(MAXIMGBATCHSZ), source=0.0)
                allocate(vld_es(MAXIMGBATCHSZ), source=.false.)
            endif
            ! mean deflation of the M-step residual (SIMPLE_COV_PROBE_DEFLATE=0 to skip):
            ! a geometry defect used to skip it silently on the fused path, and the
            ! no-deflation variant measured BETTER twice (ARI 0.970/0.981 incl. 16/16 GT
            ! states, +2 probe rounds) -- kept as an explicit arm until the A/B settles it.
            ! presence-and-zero test: cov_env_int silently ignores <= 0 (the CGSOLVE trap)
            l_deflate = .not. cov_env_int_off('SIMPLE_COV_PROBE_DEFLATE')
            if( l_gpu_estep .and. .not. l_deflate )then
                write(logfhandle,'(A)') '>>> FLEX_PCA PROBE M-STEP DEFLATION OFF (raw-data subspace iteration)'
                call flush(logfhandle)
            endif
            ! device prep: mask path only, no sigma2 whitening; planes born on device
            l_gpu_prep = .false.
            if( l_gpu_estep )then
                vprep = 1
                call cov_env_int('SIMPLE_COV_GPU_PREP', vprep)
                l_gpu_prep = vprep > 0
                if( l_gpu_prep .and. params%l_ml_reg )then
                    ! whitening path: needs the loaded sigma2 spectra (the CPU prep THROWs
                    ! without them, so this mirrors its requirement)
                    l_gpu_prep = allocated(build%esig%sigma2_noise)
                endif
            endif
            if( l_gpu_prep )then
                kf_pb = projected_model_kfromto(params)
                call fit_pb%new([params%boxpd, params%boxpd, 1], params%smpd_crop)
                frlims_pb = fit_pb%loop_lims(3)
                ! fill to the FULL padded band like the CPU generator (the working-band cap
                ! acts downstream through fpl%nyq, not through the stored values)
                nyqpd_pb  = fit_pb%get_lfny(1)
                call flex_gpu_prep_begin_f(build%lmsk, params%box, params%boxpd, &
                    &MAXIMGBATCHSZ, cov_image_mask_radius(params), .true.)
                allocate(ctfp_arr(MAXIMGBATCHSZ), shf_pb(2,MAXIMGBATCHSZ))
                if( params%l_ml_reg )then
                    block
                        type(ftiter) :: fit_cr
                        call fit_cr%new([params%box_crop, params%box_crop, 1], params%smpd_crop)
                        signyq_pb = fit_cr%get_lfny(1)
                    end block
                    nyqfull_pb = fit_pb%get_lfny(1)
                    allocate(sig2_pb(0:nyqfull_pb, MAXIMGBATCHSZ), source=1.0)
                endif
                write(logfhandle,'(A)') '>>> FLEX_PCA PROBE DEVICE PREP ON (planes born on device)'
                call flush(logfhandle)
            endif
            allocate(o_pb_thr(nthr))
            allocate(cap1(MAXIMGBATCHSZ), source=1.0)
            allocate(sap0(MAXIMGBATCHSZ), source=0.0)
            allocate(seg_dir_b(MAXIMGBATCHSZ), seg_beg_b(MAXIMGBATCHSZ), seg_cnt_b(MAXIMGBATCHSZ))
            allocate(rmat_pp(3,3,npp), nrm_pp(3,npp), dirof_pb(npp), cap(npp), sap(npp))
            do i = 1, npp
                call build%spproj_field%get_ori(ppinds(i), o_pb)
                rmat_pp(:,:,i) = o_pb%get_mat()
                nrm_pp(:,i)    = rmat_pp(3,:,i)
            end do
            call o_pb%kill
            call polar_assign_directions(nrm_pp, npp, nrm_pb, ndir_pb, dirof_pb)
            do i = 1, npp
                call polar_relative_inplane(rmat_pp(:,:,i), rmat_pb(:,:,dirof_pb(i)), ca1, sa1)
                cap(i) = ca1; sap(i) = sa1
            end do
            deallocate(rmat_pp, nrm_pp)
            ! stable counting sort by bank direction: contiguous same-direction runs inside each
            ! batch are what turn per-particle splats into per-direction splats
            allocate(cnt_pb(ndir_pb), source=0)
            do i = 1, npp
                cnt_pb(dirof_pb(i)) = cnt_pb(dirof_pb(i)) + 1
            end do
            allocate(perm_pb(npp))
            block
                integer, allocatable :: pos(:)
                allocate(pos(ndir_pb)); pos(1) = 1
                do jdir = 2, ndir_pb
                    pos(jdir) = pos(jdir-1) + cnt_pb(jdir-1)
                end do
                do i = 1, npp
                    perm_pb(pos(dirof_pb(i))) = i
                    pos(dirof_pb(i)) = pos(dirof_pb(i)) + 1
                end do
                deallocate(pos)
            end block
            allocate(pptmp(npp))
            pptmp = ppinds;   ppinds   = pptmp(perm_pb)
            pptmp = dirof_pb; dirof_pb = pptmp(perm_pb)
            deallocate(pptmp)
            block
                real, allocatable :: rtmp(:)
                allocate(rtmp(npp))
                rtmp = cap; cap = rtmp(perm_pb)
                rtmp = sap; sap = rtmp(perm_pb)
                deallocate(rtmp)
            end block
            deallocate(perm_pb, cnt_pb)
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA PROBE BANKED ADJOINT ON (', &
                &ndir_pb,' bank directions, ',npp,' records)'
            call flush(logfhandle)
        endif
        if( .not. allocated(sec_proj_thr) )then
            allocate(sec_proj_thr(nthr), sec_gram_thr(nthr), source=0.d0)
        endif
        nll_prev    = 0.d0
        eo_best     = -1.d0
        eo_stall    = 0
        eo_patience = COV_EO_PATIENCE
        vpat        = 0
        call cov_env_int('SIMPLE_COV_EO_PATIENCE', vpat)
        if( vpat > 0 ) eo_patience = vpat
        ! ON by default. The consensus-shaped term is one the model provably cannot represent:
        ! the contrast coefficient is pinned at 1, so a particle of true amplitude a_i leaves
        ! (a_i - 1) * T_i P(R_i) mu behind, rank one along the consensus and identical in every
        ! particle. Deflating it costs a projection per basis volume per iteration and was
        ! positive in all four cells measured at K=20 on EMPIAR-10076 -- moment 0.1756 -> 0.1906,
        ! EM 0.1193 -> 0.1764, and positive again with latent whitening layered on both arms.
        ! It removes 34.6% of the basis energy on the first EM iteration and 16.2% at convergence.
        ! Validated on 10076 only so far; SIMPLE_COV_EM_DEFLATE=0 restores the old behaviour.
        vdfl        = 1
        call cov_env_int('SIMPLE_COV_EM_DEFLATE', vdfl)
        ! cov_env_int only accepts values > 0, so it cannot express "off" and =0 would silently
        ! leave the default standing. The explicit off-check is what makes the escape hatch real.
        l_deflate_mean = vdfl > 0 .and. .not. cov_env_int_off('SIMPLE_COV_EM_DEFLATE')
        ! ---- per-particle contrast INSIDE the EM (a_i in the basis loop) ----
        ! The probe has always fit a scale a_i = <m,y>/||m||^2 per particle per iteration and
        ! subtracted a_i*(T mu) from the residual -- the claim that the EM path has no scale
        ! parameter was about a different routine. What the historical path does NOT do:
        !   1. refine a_i against the basis: the projection fit ignores B z, which biases a_i
        !      wherever the basis is not orthogonal to the consensus (deflation removes most of
        !      that, which is one reason the two compose), and
        !   2. carry a_i into the M-step: the residual y - a m ~ a B z + n is inserted with
        !      weight z and density E[zz'], where the ML weighting is a z and a^2 E[zz'].
        ! 3DVA fits its alpha_i jointly in the M-step, which is both of these at once.
        ! SIMPLE_COV_PROBE_CONTRAST=n enables n ECM alternations per particle (fixes 1);
        ! SIMPLE_COV_PROBE_MLSCALE=1 enables the consistent M-step weighting (fixes 2).
        ! Both off by default; the polar statistics path keeps the fixed polar-fit contrast.
        n_probe_cm = 0
        call cov_env_int('SIMPLE_COV_PROBE_CONTRAST', n_probe_cm)
        n_probe_cm = max(0, min(8, n_probe_cm))
        nml_plain = n_probe_cm
        if( l_em_polar ) nml_plain = 0
        l_probe_mls = cov_env_int_on('SIMPLE_COV_PROBE_MLSCALE')
        if( n_probe_cm > 0 ) write(logfhandle,'(A,I0,A)') &
            &'>>> FLEX_PCA PROBE per-particle contrast ECM ON (', n_probe_cm, ' alternations)'
        if( l_probe_mls ) write(logfhandle,'(A)') &
            &'>>> FLEX_PCA PROBE ML-scaled M-step ON (a z insertion, a^2 density)'
        vanneal   = 0
        call cov_env_int('SIMPLE_COV_EM_ANNEAL', vanneal)
        kfr_ann   = covariance_kfromto(params)
        khi_full  = max(1, kfr_ann(2))
        dstep_ann = real(max(1, params%box_crop - 1)) * params%smpd_crop
        conv_thresh = COV_PROBE_CONV
        vconv       = 0
        call cov_env_int('SIMPLE_COV_PROBE_CONV', vconv)
        if( vconv > 0 ) conv_thresh = min(0.999999d0, real(vconv,dp)/1000.d0)
        ! ---- DOWNSCALED-PARTICLE CACHE ----
        ! Every EM iteration re-reads and re-preps the SAME particles, and one qsys round is
        ! launched per iteration, so the workers are fresh processes each time and nothing held in
        ! memory survives them. The on-disk cache does survive: it stores the iteration-independent
        ! prefix (noise normalisation, FFT, crop to box_crop), which at box=360/box_crop=64 is where
        ! the measured 48% of probe time goes. Masking and the device prep both keep the full box.
        l_pcache = ptcl_cache_in_use(params, build) .and. .not. l_gpu_prep .and. &
            &cov_image_mask_radius(params) <= 0.0
        if( l_pcache )then
            write(logfhandle,'(A)') '>>> FLEX_PCA PROBE: particles served from the downscaled cache'
            call flush(logfhandle)
        endif
        do it = 1, niters
            lp_it = params%lp
            if( vanneal > 0 )then
                khi_it = max(1, min(khi_full, 1 + ceiling(real(khi_full - 1) * &
                    &min(1.0, real(it) / real(vanneal)))))
                lp_it  = max(params%lp, dstep_ann / real(khi_it))
                write(logfhandle,'(A,I0,A,I0,A,I0,A,F7.2)') '>>> FLEX_PCA PROBE ITER ',it, &
                    &'  band anneal k_hi=',khi_it,' / ',khi_full,'  lp=',lp_it
            endif
            t_it = tic()
            sec_read = 0.; sec_prep = 0.; sec_estep = 0.; sec_ins = 0.
            sec_proj_thr = 0.d0; sec_gram_thr = 0.d0
            sec_bank = 0.; l_pol_bank_it = .false.   ! the polar E-step bank is rebuilt every iteration
            write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA PROBE SUBSPACE ITERATION ',it,' / ',niters, &
                &'  basis dim=',ncomp
            call flush(logfhandle)
            allocate(prior(ncomp))
            do q = 1, ncomp
                prior(q) = 1.d0 / max(eigvals(q), DTINY)
            end do
            ! MCFA bookkeeping. l_mix_used freezes which prior THIS iteration's E-step runs
            ! under (the M-step hook below flips l_mix_active mid-iteration); a basis dimension
            ! change invalidates xi/Omega, so the mixture re-warms and re-initialises.
            l_mix_used = l_mix_active
            ldOm_used  = ldOm_mix
            if( l_mix_req )then
                if( allocated(rhs0th) )then
                    if( size(rhs0th,1) /= ncomp )then
                        deallocate(rhs0th, mkth, lwth, rkth, mxa_sr, mxa_sm, mxa_smm, mxa_sainv)
                        if( allocated(mix_xi) ) deallocate(mix_xi, mix_Om, mix_Ominv, mix_pi, &
                            &mix_Omxi, mix_xiOx, mix_lpi)
                        if( l_mix_active ) write(logfhandle,'(A)') &
                            &'>>> FLEX_PCA MIX re-initialising: basis dimension changed'
                        l_mix_active = .false.
                        l_mix_used   = .false.
                    endif
                endif
                if( .not. allocated(rhs0th) )then
                    allocate(rhs0th(ncomp,nthr), mkth(ncomp,kmix_pr,nthr), lwth(kmix_pr,nthr), &
                        &rkth(kmix_pr,nthr), mxa_sr(kmix_pr,nthr), mxa_sm(ncomp,kmix_pr,nthr), &
                        &mxa_smm(ncomp,ncomp,kmix_pr,nthr), mxa_sainv(ncomp,ncomp,nthr))
                endif
                mxa_sr = 0.d0; mxa_sm = 0.d0; mxa_smm = 0.d0; mxa_sainv = 0.d0
            endif
            ! even/odd Y_q accumulators (half-set FSC regularization) + the COUPLED latent normal matrix.
            ! rho carries one entry per (q,r) pair, not one shared density: the M-step below solves the
            ! components together at every grid point.
            allocate(Yeven(ncomp), Yodd(ncomp))
            do q = 1, ncomp
                call init_column_reconstructor(params, build, Yeven(q)); call Yeven(q)%reset; call Yeven(q)%reset_exp
                call init_column_reconstructor(params, build, Yodd(q));  call Yodd(q)%reset;  call Yodd(q)%reset_exp
            end do
            es     = shape(Yeven(1)%cmat_exp)
            ! The per-voxel coupled normal matrix is always the FULL ncomp x ncomp SPD system.
            ! A diagonal approximation used to be selectable here; it dropped the cross-component
            ! terms, and on EMPIAR-10028 that cost the 40S rotation outright -- the mode only
            ! appeared once the full solve was restored. It is not a speed/accuracy trade worth
            ! offering, so the option is gone rather than merely defaulted off.
            npairs = (ncomp*(ncomp+1))/2
            allocate(rho_e(npairs,es(1),es(2),es(3)), rho_o(npairs,es(1),es(2),es(3)), source=0.)
            write(logfhandle,'(A,I0,A,F8.2,A)') '>>> FLEX_PCA PROBE coupled normal matrix rows=',npairs, &
                &' (full)  rho even+odd ', 8.d0*real(npairs,dp)*real(es(1),dp)*real(es(2),dp)*real(es(3),dp)/1.d9,' GB'
            call flush(logfhandle)
            allocate(Gth(ncomp,ncomp,nthr), Ath(ncomp,ncomp,nthr), bth(ncomp,nthr), cth(ncomp,nthr), zth(ncomp,nthr))
            allocate(Ainvth(ncomp,ncomp,nthr), Acpth(ncomp,ncomp,nthr))
            allocate(basis_fpls(ncomp,nthr), mean_fpl(nthr), orientations(MAXIMGBATCHSZ))
            if( l_bank_fwd ) allocate(basis_fplsC(ncomp,nthr), mean_fplC(nthr))
            allocate(zbatch(ncomp,MAXIMGBATCHSZ), dens(ncomp,ncomp,MAXIMGBATCHSZ))
            allocate(valid(MAXIMGBATCHSZ), valid_e(MAXIMGBATCHSZ), valid_o(MAXIMGBATCHSZ), eo(MAXIMGBATCHSZ))
            allocate(gam_thr(ncomp,nthr), source=0.d0)
            allocate(gam_acc(ncomp), source=0.d0)
            allocate(nval_thr(nthr), source=0)
            allocate(hth(ncomp,nthr), source=0.d0)
            allocate(nll_thr(nthr), source=0.d0)
            allocate(gam_dbg(4,nthr), source=0.d0)
            nll_tot = 0.d0
            dens = 0.d0
            call mean_rec%expand_exp
            do q = 1, ncomp
                call basis_recs(q)%expand_exp
            end do
            ! ---- POLAR E-STEP (SIMPLE_COV_EM_POLAR=1) ----
            ! Must run BEFORE init_rec: embed_accumulate_polar ends with cleanup_rec_buffers, which
            ! calls killimgbatch / dealloc_imgarr(img_pad_heap) / forget_ft_maps -- shared BUILD
            ! state, not just its own planes -- so running it after init_rec pulls the batch loop's
            ! buffers out from under it and segfaults on the first insertion.
            ! One pass through the shared-direction bank fills every particle's G, b, c and
            ! contrast; the batch loop below then projects the MEAN ONLY, to form the Cartesian
            ! residual the M-step adjoint needs. There is no polar->volume adjoint in the tree, so
            ! the M-step stays Cartesian by necessity.
            if( l_em_polar )then
                if( l_probe_distr_pre ) THROW_HARD('SIMPLE_COV_EM_POLAR is shared-memory only for &
                    &now; the polar former would have to run inside the worker. Use nparts=1.')
                allocate(pGc(ncomp,ncomp,npp), source=0.d0)
                allocate(pbc(ncomp,npp), pcc(ncomp,npp), source=0.d0)
                allocate(pzh(npp,ncomp,2), source=0.d0)
                allocate(pcon(npp), pre(npp), prme(npp), source=0.d0)
                call embed_accumulate_polar(params, build, mean_rec, basis_recs, ncomp, eigvals, &
                    &sig2, ppinds, npp, pGc, pbc, pcc, pzh, pcon, pre, prme, prior, nvalid_pol)
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA PROBE POLAR E-STEP: statistics from &
                    &the shared-direction bank for ',nvalid_pol,' of ',npp
                call flush(logfhandle)
            endif
            call init_rec(params, build, MAXIMGBATCHSZ, fpls, cached=l_pcache)
            ! ---- ACCUMULATE: distributed when the master has parts, in-process otherwise ----
            ! One qsys round per EM iteration, because the basis the E-step projects changes every
            ! iteration: workers are relaunched against the master's refreshed flex_pca_pc*.mrc
            ! rather than looping locally. Everything below the reduction -- the coupled M-step
            ! solve, the FSC merge, the re-orthonormalisation -- stays master-only and unchanged, so
            ! the shared-memory result remains the reference.
            l_probe_distr = flex_pca_is_master() .and. flex_pca_nparts() > 1
            allocate(gam_sum(ncomp), source=0.d0)
            nval = 0
            if( l_probe_distr )then
                call save_probe_state(ncomp, eigvals, sig2_eff)
                call flex_pca_run_stage(PCA_STAGE_PROBE, 'probe iteration')
                allocate(cme(es(1),es(2),es(3),ncomp), cmo(es(1),es(2),es(3),ncomp), source=(0.,0.))
                allocate(rhe(es(1),es(2),es(3),ncomp), rhoo(es(1),es(2),es(3),ncomp), source=0.)
                if( l_mix_req )then
                    if( .not. allocated(dm_sr) )then
                        allocate(dm_sr(kmix_pr), dm_sm(ncomp,kmix_pr), &
                            &dm_smm(ncomp,ncomp,kmix_pr), dm_sai(ncomp,ncomp))
                        allocate(dm_z(MIX_ZSUB_MAX*flex_pca_nparts(), ncomp))
                    endif
                    dm_sr = 0.d0; dm_sm = 0.d0; dm_smm = 0.d0; dm_sai = 0.d0; dm_nz = 0
                    call reduce_probe_parts(params, flex_pca_nparts(), cme, rhe, cmo, rhoo, &
                        &rho_e, rho_o, gam_sum, nll_tot, nval, ncomp, &
                        &dm_sr, dm_sm, dm_smm, dm_sai, dm_z, dm_nz)
                    write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA MIX distributed reduce: kmix=', &
                        &kmix_pr,'  pooled latent subsample rows=',dm_nz
                    call flush(logfhandle)
                else
                    call reduce_probe_parts(params, flex_pca_nparts(), cme, rhe, cmo, rhoo, &
                        &rho_e, rho_o, gam_sum, nll_tot, nval, ncomp)
                endif
                do q = 1, ncomp
                    Yeven(q)%cmat_exp = cme(:,:,:,q); Yeven(q)%rho_exp = rhe(:,:,:,q)
                    Yodd(q)%cmat_exp  = cmo(:,:,:,q); Yodd(q)%rho_exp  = rhoo(:,:,:,q)
                end do
                deallocate(cme, cmo, rhe, rhoo)
            else
            if( l_pcache )then
                call prepimgbatch(params, build, MAXIMGBATCHSZ, box=params%box_crop, smpd=params%smpd_crop)
            else
                call prepimgbatch(params, build, MAXIMGBATCHSZ)
            endif
            vgpu = 0
            call cov_env_int('SIMPLE_COV_GPU', vgpu)
            l_gpu_probe  = vgpu > 0 .and. flex_gpu_available()
            if( l_gpu_probe )then
                write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA PROBE M-step GPU insertion ON (', &
                    &2*ncomp, ' components, ', 2*npairs, ' density rows)'
                call flush(logfhandle)
                call flex_gpu_coupled_begin_f(Yeven, 2*ncomp, 2*npairs)
                allocate(gdsc(2*ncomp,MAXIMGBATCHSZ), grsc(2*npairs,MAXIMGBATCHSZ), source=0.d0)
                if( l_bank_adj ) call flex_gpu_coupled_bank_f(rmat_pb, ndir_pb)
                if( l_gpu_estep .and. l_bank_adj )then
                    call flex_gpu_estep_vols_f(mean_rec, basis_recs, ncomp)
                    nst_es = 2 + 2*ncomp + (ncomp*(ncomp+1))/2
                    allocate(stats_es(nst_es, MAXIMGBATCHSZ), source=0.d0)
                    write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PROBE FUSED E-STEP ON (', &
                        &ncomp+1, ' volumes resident, exact directions)'
                    call flush(logfhandle)
                endif
            endif
            z = 0.d0
            ! parity-check bookkeeping: only iteration 1, only when the polar former runs
            l_chk_act = l_pol_es .and. n_pol_check > 0 .and. it == 1
            if( l_chk_act .and. .not. allocated(zchk_c) )then
                nchk_max = min(n_pol_check, npp)
                allocate(zchk_c(ncomp,nchk_max), zchk_p(ncomp,nchk_max), source=0.d0)
                allocate(gerr_chk(nchk_max), berr_chk(nchk_max), source=0.d0)
                allocate(lchk(nchk_max), source=.false.)
                allocate(Gc_chk(ncomp,ncomp,nthr), bc_chk(ncomp,nthr), cc_chk(ncomp,nthr))
                allocate(zc_chk(ncomp,nthr), Ai_chk(ncomp,ncomp,nthr))
                allocate(zchk_d(ncomp,nchk_max), source=0.d0)
                allocate(gerrd_chk(nchk_max), berrd_chk(nchk_max), source=0.d0)
                allocate(Gd_chk(ncomp,ncomp,nthr), bd_chk(ncomp,nthr), cd_chk(ncomp,nthr))
                allocate(zd_chk(ncomp,nthr), dcb_chk(ncomp,nthr))
            endif
            ob_lo = 1
            ob_hi = npp
            if( l_online )then
                b_win     = max(1, b_online / max(1, params%nparts))   ! params%nparts is worker-visible; flex_pca_nparts() is master-only
                nb_online = max(1, (npp + b_win - 1) / b_win)
                it_win    = mod(it - 1, nb_online)
                ob_lo     = 1 + it_win * b_win
                ob_hi     = min(npp, ob_lo + b_win - 1)
                if( l_onl_final .and. it == niters )then
                    ob_lo = 1
                    ob_hi = npp
                    write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA ONLINE FINAL it=', it, &
                        &'  full-set pass over ', npp, ' particles'
                else
                    write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA ONLINE it=', it, &
                        &'  window ', ob_lo, '..', ob_hi, '  of ', npp
                endif
                call flush(logfhandle)
            endif
            do ibatch = ob_lo, ob_hi, MAXIMGBATCHSZ
                batchlims = [ibatch, min(ob_hi, ibatch + MAXIMGBATCHSZ - 1)]
                batchsz   = batchlims(2) - batchlims(1) + 1
                t_sec = tic()
                if( l_pcache )then
                    call ptcl_cache_read_batch(params, build, npp, ppinds, batchlims)
                else
                    call discrete_read_imgbatch(params, build, npp, ppinds, batchlims)
                endif
                sec_read = sec_read + toc(t_sec)
                t_sec = tic()
                if( .not. l_gpu_prep )then
                    call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                        &ppinds(batchlims(1):batchlims(2)), fpls(:batchsz), &
                        &mskrad=cov_image_mask_radius(params), cached=l_pcache)
                endif
                sec_prep = sec_prep + toc(t_sec)
                do i = 1, batchsz
                    call build%spproj_field%get_ori(ppinds(batchlims(1)+i-1), orientations(i))
                    eo(i) = build%spproj_field%get_eo(ppinds(batchlims(1)+i-1))
                    if( l_gpu_prep )then
                        ctfp_arr(i)  = build%spproj%get_ctfparams(params%oritype, ppinds(batchlims(1)+i-1))
                        shf_pb(:,i)  = build%spproj_field%get_2Dshift(ppinds(batchlims(1)+i-1))
                        if( params%l_ml_reg ) call resample_sigma2(kf_pb(1), signyq_pb, &
                            &build%esig%sigma2_noise(kf_pb(1):kf_pb(2), ppinds(batchlims(1)+i-1)), &
                            &nyqfull_pb, real(signyq_pb)/real(nyqfull_pb), sig2_pb(:,i))
                    endif
                end do
                ! ---- POLAR E-STEP BANK, built once per EM iteration at the first prepped batch
                ! (basis_recs change every M-step, so the bank cannot be cached across iterations;
                ! the grid and the pose-fixed direction assignment are built once per stage) ----
                if( l_pol_es .and. .not. l_pol_bank_it )then
                    t_bank = tic()
                    if( .not. l_pol_grid )then
                        ! grid geometry from the first prepped plane: the identical derivation
                        ! embed_accumulate_polar uses, so polar embed and polar E-step share
                        ! band/quadrature conventions. No noise rings: sig2 arrives as sig2_eff.
                        ph0_es  = lbound(fpls(1)%cmplx_plane,1)
                        pk0_es  = lbound(fpls(1)%cmplx_plane,2)
                        hlo_es  = ceil_div (lbound(fpls(1)%cmplx_plane,1), OSMPL_PAD_FAC)
                        hhi_es  = floor_div(ubound(fpls(1)%cmplx_plane,1), OSMPL_PAD_FAC)
                        klo_es  = ceil_div (lbound(fpls(1)%cmplx_plane,2), OSMPL_PAD_FAC)
                        nyqr_es = mean_rec%get_lfny(1)
                        nyqb_es = nyqr_es
                        if( fpls(1)%nyq > 0 ) nyqb_es = min(nyqb_es, max(1, fpls(1)%nyq / OSMPL_PAD_FAC))
                        ! hybrid split point: auto at 0.72*band -- the measured knee of the
                        ! real-data ladder (10049 gate, band 11: rhyb 6 -> b err 8.9%, min z
                        ! corr 0.949, lambda head 153/282/326; rhyb 8 -> b err 3.1%, min z
                        ! corr 0.990, lambda 166/273/394 vs Cartesian 166/264/382). Env
                        ! override; explicit 0 = pure rings (the pre-hybrid baseline).
                        rhyb_es = nint(0.72*real(nyqb_es))
                        if( vpol_rhyb > 0 ) rhyb_es = vpol_rhyb
                        if( l_rhyb_off )    rhyb_es = 0
                        rhyb_es   = max(0, min(rhyb_es, nyqb_es-1))
                        l_pol_hyb = rhyb_es > 0
                        if( l_pol_hyb )then
                            call polar_grid_build(pg_es, rhyb_es+1, nyqb_es, nyqb_es+1, nyqb_es, &
                                &hlo_es, hhi_es, klo_es, ph0_es, pk0_es, ang_osamp=vpol_osamp, &
                                &gate_lo=rhyb_es*(rhyb_es+1))
                            ! exact-part lattice positions, cov_herm_inner's half-plane rule
                            ! (k<=0; on the k=0 line only h<=0), shells 0..rhyb by the nint
                            ! convention: h^2+k^2 <= rhyb*(rhyb+1). Raster order k-outer.
                            npos_es = 0
                            do kx_es = -rhyb_es, 0
                                do jx_es = -rhyb_es, merge(0, rhyb_es, kx_es == 0)
                                    if( jx_es*jx_es + kx_es*kx_es > rhyb_es*(rhyb_es+1) ) cycle
                                    npos_es = npos_es + 1
                                end do
                            end do
                            allocate(hex_es(npos_es), kex_es(npos_es))
                            npos_es = 0
                            do kx_es = -rhyb_es, 0
                                do jx_es = -rhyb_es, merge(0, rhyb_es, kx_es == 0)
                                    if( jx_es*jx_es + kx_es*kx_es > rhyb_es*(rhyb_es+1) ) cycle
                                    npos_es = npos_es + 1
                                    hex_es(npos_es) = jx_es
                                    kex_es(npos_es) = kx_es
                                end do
                            end do
                            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA POLAR HYBRID: exact &
                                &Cartesian statistics for shells 0..',rhyb_es,' (',npos_es,' lattice &
                                &points/particle incl. DC), rings ',rhyb_es+1,'..',nyqb_es
                            if( vpol_osamp > 1 ) write(logfhandle,'(A,I0)') &
                                &'>>> FLEX_PCA POLAR OSAMP: ring angular oversampling x', vpol_osamp
                            call flush(logfhandle)
                        else
                            call polar_grid_build(pg_es, 1, nyqb_es, nyqb_es+1, nyqb_es, hlo_es, hhi_es, &
                                &klo_es, ph0_es, pk0_es, ang_osamp=vpol_osamp)
                            if( vpol_osamp > 1 )then
                                write(logfhandle,'(A,I0)') &
                                    &'>>> FLEX_PCA POLAR OSAMP: ring angular oversampling x', vpol_osamp
                                call flush(logfhandle)
                            endif
                        endif
                        nsamp_es = pg_es%nsamp; nsamp2_es = 2*nsamp_es; nk_es = pg_es%nk
                        ! shared direction table: same refspiral + count derivation as the polar embed
                        ndir_es = cov_polar_ndir(npp)
                        call dirs_es%new(ndir_es, is_ptcl=.false.)
                        call build%pgrpsyms%build_refspiral(dirs_es)
                        allocate(rmatb_es(3,3,ndir_es), nrmb_es(3,ndir_es))
                        do id_es = 1, ndir_es
                            rmatb_es(:,:,id_es) = dirs_es%get_mat(id_es)
                            nrmb_es(:,id_es)    = rmatb_es(3,:,id_es)
                        end do
                        call dirs_es%kill
                        ! pose-fixed per-particle (direction, in-plane) assignment, once per stage
                        allocate(rmatp_es(3,3,npp), nrmp_es(3,npp), dir_es(npp), cae(npp), sae(npp))
                        do i = 1, npp
                            call build%spproj_field%get_ori(ppinds(i), o_es)
                            rmatp_es(:,:,i) = o_es%get_mat()
                            nrmp_es(:,i)    = rmatp_es(3,:,i)
                        end do
                        call o_es%kill
                        call polar_assign_directions(nrmp_es, npp, nrmb_es, ndir_es, dir_es)
                        do i = 1, npp
                            call polar_relative_inplane(rmatp_es(:,:,i), rmatb_es(:,:,dir_es(i)), ca1, sa1)
                            cae(i) = ca1; sae(i) = sa1
                        end do
                        deallocate(rmatp_es, nrmp_es)
                        allocate(dused_es(ndir_es))
                        ! device ring geometry, once per stage (the grid is pose- and band-fixed)
                        if( l_pol_gpu ) call flex_gpu_poles_begin_f(pg_es%rad, pg_es%cs, &
                            &pg_es%sn, pg_es%sqwq, pg_es%rbeg, pg_es%rend, nsamp_es, nk_es)
                        l_pol_grid = .true.
                        write(logfhandle,'(A,I0,A,I0,A,I0,A,F8.1,A,F8.1,A)') &
                            &'>>> FLEX_PCA POLAR ESTEP BANK: ',ncomp+1,' volumes x ',ndir_es, &
                            &' directions x ',nsamp_es,' ring samples = ', &
                            &4.d0*real(nsamp2_es,dp)*real(ncomp+1,dp)*real(ndir_es,dp)/1.d6, &
                            &' MB (+ ring tables ', &
                            &8.d0*real(ncomp*ncomp+ncomp+1,dp)*real(nk_es,dp)*real(ndir_es,dp)/1.d6,' MB)'
                        call flush(logfhandle)
                    endif
                    ! (re)allocate at this iteration's rank (ncomp can change between iterations)
                    if( allocated(UsallE) )then
                        if( size(UsallE,2) /= ncomp+1 ) deallocate(UsallE, CfE, Cm0E, c00E, &
                            &UbankE, CspE, xws_es, wr_es, wrd_es, Reb_es)
                    endif
                    if( .not. allocated(UsallE) )then
                        allocate(UsallE(nsamp2_es,0:ncomp,ndir_es), CfE(ncomp*ncomp,nk_es,ndir_es), &
                            &Cm0E(ncomp,nk_es,ndir_es), c00E(nk_es,ndir_es))
                        allocate(UbankE(nsamp_es,0:ncomp,nthr), CspE(0:ncomp,0:ncomp,nthr))
                        allocate(xws_es(nsamp2_es,nthr), wr_es(nk_es,nthr), wrd_es(nk_es,nthr), &
                            &Reb_es(0:ncomp,nthr))
                    endif
                    ! GPU parity arrays (gate 1): the CPU polar former re-run on the check
                    ! particles; the exact-direction/DC-less attribution arms stand down under
                    ! the device path (they answer the polar-vs-Cartesian question, already
                    ! measured on this recipe; the device gate is GPU-vs-CPU-polar)
                    if( l_chk_act .and. l_pol_gpu .and. .not. allocated(zchk_g) )then
                        allocate(zchk_g(ncomp,nchk_max), source=0.d0)
                        allocate(gerrg_chk(nchk_max), berrg_chk(nchk_max), source=0.d0)
                        allocate(Gpc_chk(ncomp,ncomp,nthr), bpc_chk(ncomp,nthr), &
                            &cpc_chk(ncomp,nthr), zpc_chk(ncomp,nthr))
                    endif
                    ! exact-direction check arm scratch (grid-shaped, so allocated here)
                    if( l_chk_act .and. .not. l_pol_gpu .and. .not. allocated(UsX_chk) )then
                        allocate(zchk_x(ncomp,nchk_max), source=0.d0)
                        allocate(gerrx_chk(nchk_max), berrx_chk(nchk_max), source=0.d0)
                        allocate(UsX_chk(nsamp2_es,0:ncomp,nthr), xwsX_chk(nsamp2_es,nthr), &
                            &wrX_chk(nk_es,nthr), RebX_chk(0:ncomp,nthr))
                        allocate(CfX_chk(ncomp*ncomp,nk_es,nthr), Cm0X_chk(ncomp,nk_es,nthr), &
                            &c00X_chk(nk_es,nthr), wrdX_chk(nk_es,nthr))
                        allocate(GX_chk(ncomp,ncomp,nthr), bX_chk(ncomp,nthr), cX_chk(ncomp,nthr), &
                            &zx_chk(ncomp,nthr))
                        allocate(tazsum_thr(nthr), source=0.d0)
                        allocate(tazcnt_thr(nthr), source=0)
                    endif
                    ! only the directions this iteration's window touches
                    dused_es = .false.
                    do i = ob_lo, ob_hi
                        if( dir_es(i) > 0 ) dused_es(dir_es(i)) = .true.
                    end do
                    !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
                    !$omp& private(id_es,ithr,q,r,ir)
                    do id_es = 1, ndir_es
                        if( .not. dused_es(id_es) ) cycle
                        ithr = omp_get_thread_num() + 1
                        call polar_project_recs(mean_rec, basis_recs, ncomp, rmatb_es(:,:,id_es), &
                            &pg_es, UbankE(:,:,ithr))
                        do q = 0, ncomp
                            do r = 1, nsamp_es
                                UsallE(2*r-1,q,id_es) = pg_es%sqwq(r)*real (UbankE(r,q,ithr))
                                UsallE(2*r,  q,id_es) = pg_es%sqwq(r)*aimag(UbankE(r,q,ithr))
                            end do
                        end do
                        do ir = 1, nk_es
                            call polar_ring_gram(UsallE(1,0,id_es), nsamp2_es, ncomp, pg_es%rbeg(ir), &
                                &pg_es%rend(ir)-pg_es%rbeg(ir)+1, CspE(0,0,ithr), CfE(1,ir,id_es), &
                                &Cm0E(1,ir,id_es))
                            c00E(ir,id_es) = polar_ring_selfpower(UsallE(1,0,id_es), nsamp2_es, &
                                &pg_es%rbeg(ir), pg_es%rend(ir)-pg_es%rbeg(ir)+1)
                        end do
                    end do
                    !$omp end parallel do
                    ! ---- device upload (stage 2): the bank + ring tables once per iteration, and
                    ! the mean+basis volumes for the hybrid-exact kernel. The rank can change
                    ! between iterations, so the capacity gate is re-checked here; falling back
                    ! mid-run is seamless because every consumer below branches per iteration.
                    if( l_pol_gpu .and. l_pol_hyb .and. (ncomp + 1) > COV_GPU_ESTEP_MAXVOLS )then
                        l_pol_gpu = .false.
                        write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA POLAR E-STEP DEVICE OFF: needs ', &
                            &ncomp+1,' resident volumes, device entry point caps at ', &
                            &COV_GPU_ESTEP_MAXVOLS,' -- using the CPU polar E-step'
                        call flush(logfhandle)
                    endif
                    if( l_pol_gpu )then
                        call flex_gpu_poles_bank_f(UsallE, CfE, Cm0E, c00E, ncomp, ndir_es)
                        if( l_pol_hyb ) call flex_gpu_estep_vols_f(mean_rec, basis_recs, ncomp)
                        nst_pg = 2 + 2*ncomp + (ncomp*(ncomp+1))/2
                        if( allocated(stats_pg) )then
                            if( size(stats_pg,1) /= nst_pg ) deallocate(stats_pg)
                        endif
                        if( .not. allocated(stats_pg) ) allocate(stats_pg(nst_pg, MAXIMGBATCHSZ), &
                            &source=0.d0)
                        if( .not. allocated(vld_pg) ) allocate(vld_pg(MAXIMGBATCHSZ), &
                            &dir_pg(MAXIMGBATCHSZ), ca_pg(MAXIMGBATCHSZ), sa_pg(MAXIMGBATCHSZ))
                    endif
                    sec_bank      = sec_bank + toc(t_bank)
                    l_pol_bank_it = .true.
                    write(logfhandle,'(A,I0,A,I0,A,I0,A,F7.1)') '>>> FLEX_PCA POLAR ESTEP BANK it=', &
                        &it,'  directions built=',count(dused_es),' of ',ndir_es, &
                        &'  build seconds=',sec_bank
                    call flush(logfhandle)
                endif
                valid(:batchsz) = .false.
                zbatch(:,:batchsz) = 0.d0
                dens(:,:,:batchsz) = 0.d0
                t_sec = tic()
                if( l_gpu_estep .and. l_gpu_probe .and. l_bank_adj )then
                    ! ---- fused device E-step: stats on device at EXACT directions, tiny fetch,
                    ! host posterior solves, residual formed on device and fetched packed
                    twp0 = omp_get_wtime()
                    do i = 1, batchsz
                        vld_es(i) = .not. orientations(i)%isstatezero()
                    end do
                    if( l_gpu_prep )then
                        if( params%l_ml_reg )then
                            call flex_gpu_prep_batch_f(build%imgbatch(:batchsz), ctfp_arr(:batchsz), &
                                &shf_pb(:,:batchsz), vld_es(:batchsz), batchsz, params%box, &
                                &frlims_pb, nyqpd_pb, sig2_ups=sig2_pb(:,:batchsz))
                        else
                            call flex_gpu_prep_batch_f(build%imgbatch(:batchsz), ctfp_arr(:batchsz), &
                                &shf_pb(:,:batchsz), vld_es(:batchsz), batchsz, params%box, &
                                &frlims_pb, nyqpd_pb)
                        endif
                        ! one-time in-situ cross-check against the CPU prep on real particles
                        if( it == 1 .and. ibatch == 1 )then
                            vprep = 0
                            call cov_env_int('SIMPLE_COV_PREP_CHECK', vprep)
                            if( vprep > 0 )then
                                call prep_imgs4projected_model(params, build, batchsz, &
                                    &build%imgbatch(:batchsz), ppinds(batchlims(1):batchlims(2)), &
                                    &fpls(:batchsz), mskrad=cov_image_mask_radius(params), &
                                    &force_cpu=.true.)
                                call flex_gpu_prep_check_f(fpls(:batchsz), vld_es(:batchsz), batchsz)
                            endif
                        endif
                        call flex_gpu_estep_batch_res_f(orientations(:batchsz), batchsz, &
                            &frlims_pb, kf_pb(2), stats_es(:,:batchsz))
                    else
                        call flex_gpu_estep_batch_f(fpls(:batchsz), orientations(:batchsz), &
                            &vld_es(:batchsz), batchsz, stats_es(:,:batchsz))
                    endif
                    sec_proj_thr(1) = sec_proj_thr(1) + (omp_get_wtime() - twp0)
                    twp0 = omp_get_wtime()
                    !$omp parallel do default(shared) schedule(static) proc_bind(close) &
                    !$omp& private(i,ithr,row,q,r,off_es,a,aa,e_mm,myv,ldA,lok,qml)
                    do i = 1, batchsz
                        if( .not. vld_es(i) ) cycle
                        ithr = omp_get_thread_num() + 1
                        row  = batchlims(1) + i - 1
                        e_mm = stats_es(1,i)
                        myv  = stats_es(2,i)
                        a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
                        aa   = a*a
                        do q = 1, ncomp
                            bth(q,ithr) = stats_es(2+q,i)
                            cth(q,ithr) = stats_es(2+ncomp+q,i)
                        end do
                        off_es = 2 + 2*ncomp
                        do q = 1, ncomp
                            do r = q, ncomp
                                off_es = off_es + 1
                                Gth(q,r,ithr) = stats_es(off_es,i)
                                Gth(r,q,ithr) = Gth(q,r,ithr)
                            end do
                        end do
                        if( l_mix_used )then
                            ! same shared MCFA solver as the plain body: the device fetched the
                            ! identical G/b/c sufficient statistics, only the host solve differs
                            call probe_solve_mix(ncomp, kmix_pr, Gth(:,:,ithr), bth(:,ithr), &
                                &cth(:,ithr), myv, e_mm, n_probe_cm, a, sig2, mix_Ominv, mix_Omxi, &
                                &mix_lpi, mix_xiOx, zth(:,ithr), Ainvth(:,:,ithr), dens(:,:,i), &
                                &ldA, lok, nll_mix_add, mxa_sr(:,ithr), mxa_sm(:,:,ithr), &
                                &mxa_smm(:,:,:,ithr), mxa_sainv(:,:,ithr))
                            if( lok ) nll_thr(ithr) = nll_thr(ithr) + nll_mix_add
                        else
                            call probe_solve_ecm(ncomp, Gth(:,:,ithr), bth(:,ithr), cth(:,ithr), &
                                &myv, e_mm, prior, sig2, n_probe_cm, a, zth(:,ithr), &
                                &Ainvth(:,:,ithr), ldA, lok, qml)
                            if( lok ) nll_thr(ithr) = nll_thr(ithr) + ldA - qml
                        endif
                        aa = a*a
                        z(row,:)    = zth(:,ithr)
                        zbatch(:,i) = zth(:,ithr)
                        if( .not. l_mix_used )then
                            do r = 1, ncomp
                                do q = 1, ncomp
                                    dens(q,r,i) = zth(q,ithr)*zth(r,ithr) + Ainvth(q,r,ithr)
                                end do
                            end do
                        endif
                        do q = 1, ncomp
                            gam_thr(q,ithr) = gam_thr(q,ithr) + dens(q,q,i)
                        end do
                        ! Gamma above reads the UNSCALED E[zz'] (prior on z, not on a z);
                        ! everything below feeds the M-step, which under the fitted scale
                        ! solves  r ~ a B z + n
                        if( l_probe_mls )then
                            zbatch(:,i) = a*zbatch(:,i)
                            dens(:,:,i) = aa*dens(:,:,i)
                        endif
                        nval_thr(ithr) = nval_thr(ithr) + 1
                        valid(i)       = .true.
                        avec_es(i)     = real(a)
                    end do
                    !$omp end parallel do
                    ! residual stays on device; the M-step consumes it via the resident entry
                    if( l_deflate )then
                        call flex_gpu_estep_resid_f(fpls(:batchsz), avec_es(:batchsz), &
                            &vld_es(:batchsz), batchsz, fetch=.false.)
                    endif
                    sec_gram_thr(1) = sec_gram_thr(1) + (omp_get_wtime() - twp0)
                elseif( l_pol_gpu )then
                    ! ---- POLAR E-STEP ON DEVICE (stage 2): per-batch ring + hybrid-exact
                    ! statistics in the fused stat layout; everything downstream of the
                    ! statistics -- solve, mixture, Gamma, windows, e/o, mean projection,
                    ! residual, insertion -- is the CPU polar path's, verbatim.
                    twp0 = omp_get_wtime()
                    do i = 1, batchsz
                        row = batchlims(1) + i - 1
                        vld_pg(i) = .not. orientations(i)%isstatezero()
                        dir_pg(i) = dir_es(row)
                        ca_pg(i)  = cae(row)
                        sa_pg(i)  = sae(row)
                    end do
                    call flex_gpu_poles_batch_f(fpls(:batchsz), orientations(:batchsz), &
                        &ca_pg(:batchsz), sa_pg(:batchsz), dir_pg(:batchsz), vld_pg(:batchsz), &
                        &batchsz, merge(rhyb_es*(rhyb_es+1), 0, l_pol_hyb), stats_pg(:,:batchsz))
                    sec_proj_thr(1) = sec_proj_thr(1) + (omp_get_wtime() - twp0)
                    twp0 = omp_get_wtime()
                    !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
                    !$omp& private(i,ithr,row,q,r,off_es,a,aa,e_mm,myv,ldA,lok,qml,nll_mix_add) &
                    !$omp& private(idp_es,taz_es,l_chk_p,ichk,emmc_chk,myvc_chk,ac_chk,lda_chk) &
                    !$omp& private(qml_chk,lok_chk,gnum_chk,gden_chk,bnum_chk,bden_chk) &
                    !$omp& private(emmg_chk,myvg_chk,ag_chk)
                    do i = 1, batchsz
                        if( .not. vld_pg(i) ) cycle
                        ithr = omp_get_thread_num() + 1
                        row  = batchlims(1) + i - 1
                        e_mm = stats_pg(1,i)
                        myv  = stats_pg(2,i)
                        a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
                        aa   = a*a
                        do q = 1, ncomp
                            bth(q,ithr) = stats_pg(2+q,i)
                            cth(q,ithr) = stats_pg(2+ncomp+q,i)
                        end do
                        off_es = 2 + 2*ncomp
                        do q = 1, ncomp
                            do r = q, ncomp
                                off_es = off_es + 1
                                Gth(q,r,ithr) = stats_pg(off_es,i)
                                Gth(r,q,ithr) = Gth(q,r,ithr)
                            end do
                        end do
                        l_chk_p = l_chk_act .and. row <= nchk_max
                        if( l_chk_p )then
                            ! Cartesian reference former (its mean plane doubles as the
                            ! residual's mean projection below)
                            call project_fplanes_mean_basis(mean_rec, basis_recs, orientations(i), &
                                &fpls(i), mean_fpl(ithr), basis_fpls(:,ithr), apply_ctf_amp=.true.)
                            emmc_chk = real(cov_herm_inner(mean_fpl(ithr), mean_fpl(ithr)), dp)
                            myvc_chk = real(cov_herm_inner(mean_fpl(ithr), fpls(i)), dp)
                            do q = 1, ncomp
                                bc_chk(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), fpls(i)), dp)
                                cc_chk(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), mean_fpl(ithr)), dp)
                                do r = q, ncomp
                                    Gc_chk(q,r,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), &
                                        &basis_fpls(r,ithr)), dp)
                                    Gc_chk(r,q,ithr) = Gc_chk(q,r,ithr)
                                end do
                            end do
                            ac_chk = max(0.1d0, min(5.0d0, myvc_chk / max(emmc_chk, DTINY)))
                            call probe_solve_ecm(ncomp, Gc_chk(:,:,ithr), bc_chk(:,ithr), &
                                &cc_chk(:,ithr), myvc_chk, emmc_chk, prior, sig2, 0, ac_chk, &
                                &zc_chk(:,ithr), Ai_chk(:,:,ithr), lda_chk, lok_chk, qml_chk)
                            ichk = row
                            zchk_c(:,ichk) = zc_chk(:,ithr)
                            ! GPU-vs-Cartesian statistic errors (the polar-vs-Cartesian bar)
                            gnum_chk = 0.d0; gden_chk = 0.d0; bnum_chk = 0.d0; bden_chk = 0.d0
                            do r = 1, ncomp
                                bnum_chk = bnum_chk + (bth(r,ithr) - bc_chk(r,ithr))**2
                                bden_chk = bden_chk + bc_chk(r,ithr)**2
                                do q = 1, ncomp
                                    gnum_chk = gnum_chk + (Gth(q,r,ithr) - Gc_chk(q,r,ithr))**2
                                    gden_chk = gden_chk + Gc_chk(q,r,ithr)**2
                                end do
                            end do
                            gerr_chk(ichk) = sqrt(gnum_chk / max(gden_chk, DTINY))
                            berr_chk(ichk) = sqrt(bnum_chk / max(bden_chk, DTINY))
                            ! the CPU polar former, the production stage-1 path verbatim,
                            ! for the GPU parity gate (same bank, same tables, same hybrid)
                            idp_es = dir_es(row)
                            call polar_sample_particle_fused(fpls(i)%cmplx_plane, &
                                &fpls(i)%transfer_plane, pg_es, cae(row), sae(row), &
                                &xws_es(:,ithr), wr_es(:,ithr), taz_es)
                            wrd_es(:,ithr) = real(wr_es(:,ithr), dp)
                            call sgemv('T', nsamp2_es, ncomp+1, 1.0, UsallE(1,0,idp_es), &
                                &nsamp2_es, xws_es(1,ithr), 1, 0.0, Reb_es(0,ithr), 1)
                            call dgemv('N', ncomp*ncomp, nk_es, 1.d0, CfE(1,1,idp_es), &
                                &ncomp*ncomp, wrd_es(1,ithr), 1, 0.d0, Gpc_chk(1,1,ithr), 1)
                            call dgemv('N', ncomp, nk_es, 1.d0, Cm0E(1,1,idp_es), ncomp, &
                                &wrd_es(1,ithr), 1, 0.d0, cpc_chk(1,ithr), 1)
                            emmg_chk = dot_product(c00E(:,idp_es), wrd_es(:,ithr))
                            myvg_chk = real(Reb_es(0,ithr), dp)
                            do q = 1, ncomp
                                bpc_chk(q,ithr) = real(Reb_es(q,ithr), dp)
                            end do
                            if( l_pol_hyb ) call polar_hybrid_exact_accum(mean_rec, basis_recs, &
                                &ncomp, orientations(i), fpls(i), hex_es, kex_es, npos_es, &
                                &Gpc_chk(:,:,ithr), bpc_chk(:,ithr), cpc_chk(:,ithr), emmg_chk, &
                                &myvg_chk)
                            ! GPU-vs-CPU-polar statistic errors (the gate-1 numbers)
                            gnum_chk = 0.d0; gden_chk = 0.d0; bnum_chk = 0.d0; bden_chk = 0.d0
                            do r = 1, ncomp
                                bnum_chk = bnum_chk + (bth(r,ithr) - bpc_chk(r,ithr))**2
                                bden_chk = bden_chk + bpc_chk(r,ithr)**2
                                do q = 1, ncomp
                                    gnum_chk = gnum_chk + (Gth(q,r,ithr) - Gpc_chk(q,r,ithr))**2
                                    gden_chk = gden_chk + Gpc_chk(q,r,ithr)**2
                                end do
                            end do
                            gerrg_chk(ichk) = sqrt(gnum_chk / max(gden_chk, DTINY))
                            berrg_chk(ichk) = sqrt(bnum_chk / max(bden_chk, DTINY))
                            ! CPU-polar z through the SAME solver settings as the production
                            ! solve below, so the correlation isolates the statistics
                            ag_chk = max(0.1d0, min(5.0d0, myvg_chk / max(emmg_chk, DTINY)))
                            call probe_solve_ecm(ncomp, Gpc_chk(:,:,ithr), bpc_chk(:,ithr), &
                                &cpc_chk(:,ithr), myvg_chk, emmg_chk, prior, sig2, n_probe_cm, &
                                &ag_chk, zpc_chk(:,ithr), Ai_chk(:,:,ithr), lda_chk, lok_chk, &
                                &qml_chk)
                            zchk_g(:,ichk) = zpc_chk(:,ithr)
                            lchk(ichk)     = .true.
                        endif
                        if( l_mix_used )then
                            ! same shared MCFA solver as the CPU bodies: the device fetched the
                            ! identical G/b/c sufficient statistics, only the former differs
                            call probe_solve_mix(ncomp, kmix_pr, Gth(:,:,ithr), bth(:,ithr), &
                                &cth(:,ithr), myv, e_mm, n_probe_cm, a, sig2, mix_Ominv, mix_Omxi, &
                                &mix_lpi, mix_xiOx, zth(:,ithr), Ainvth(:,:,ithr), dens(:,:,i), &
                                &ldA, lok, nll_mix_add, mxa_sr(:,ithr), mxa_sm(:,:,ithr), &
                                &mxa_smm(:,:,:,ithr), mxa_sainv(:,:,ithr))
                            if( lok ) nll_thr(ithr) = nll_thr(ithr) + nll_mix_add
                        else
                            call probe_solve_ecm(ncomp, Gth(:,:,ithr), bth(:,ithr), cth(:,ithr), &
                                &myv, e_mm, prior, sig2, n_probe_cm, a, zth(:,ithr), &
                                &Ainvth(:,:,ithr), ldA, lok, qml)
                            if( lok ) nll_thr(ithr) = nll_thr(ithr) + ldA - qml
                        endif
                        aa = a*a
                        z(row,:)    = zth(:,ithr)
                        zbatch(:,i) = zth(:,ithr)
                        ! parity check: the device-path z actually used, after the shared solve
                        if( l_chk_act )then
                            if( row <= nchk_max ) zchk_p(:,row) = zth(:,ithr)
                        endif
                        if( .not. l_mix_used )then
                            do r = 1, ncomp
                                do q = 1, ncomp
                                    dens(q,r,i) = zth(q,ithr)*zth(r,ithr) + Ainvth(q,r,ithr)
                                end do
                            end do
                        endif
                        do q = 1, ncomp
                            gam_thr(q,ithr) = gam_thr(q,ithr) + dens(q,q,i)
                        end do
                        if( l_probe_mls )then
                            zbatch(:,i) = a*zbatch(:,i)
                            dens(:,:,i) = aa*dens(:,:,i)
                        endif
                        nval_thr(ithr) = nval_thr(ithr) + 1
                        valid(i)       = .true.
                        gam_dbg(1,ithr) = gam_dbg(1,ithr) + sum([(Gth(q,q,ithr), q=1,ncomp)])
                        gam_dbg(2,ithr) = gam_dbg(2,ithr) + dot_product(bth(:,ithr), bth(:,ithr))
                        gam_dbg(3,ithr) = gam_dbg(3,ithr) + dot_product(cth(:,ithr), cth(:,ithr))
                        gam_dbg(4,ithr) = gam_dbg(4,ithr) + a
                        ! residual r_i = y - a*(T mu) in place, banded, exactly the CPU polar
                        ! path's (the check arm's Cartesian mean plane serves where present)
                        if( .not. l_chk_p )then
                            call project_fplane_mean_banded(mean_rec, orientations(i), fpls(i), &
                                &mean_fpl(ithr))
                        endif
                        call subtract_mean_banded(fpls(i), mean_fpl(ithr), real(a), nyqr_es)
                    end do
                    !$omp end parallel do
                    sec_gram_thr(1) = sec_gram_thr(1) + (omp_get_wtime() - twp0)
                elseif( l_bank_fwd )then
                ! ---- banked forward: one raw (no-CTF) projection set per direction segment;
                ! per particle, y/T are RESAMPLED into the bank frame (the same unit-tap scheme
                ! the polar former and the device align use), the CTF is applied to the shared
                ! projections as an elementwise multiply, and the residual is formed already
                ! aligned -- the banked M-step then runs at unit in-plane rotation.
                if( .not. allocated(tmpc) ) allocate(tmpc(fpls(1)%frlims(1,1):fpls(1)%frlims(1,2), &
                    &fpls(1)%frlims(2,1):0, nthr))
                nseg_b = 0
                do i = 1, batchsz
                    row = batchlims(1) + i - 1
                    if( nseg_b > 0 )then
                        if( dirof_pb(row) == seg_dir_b(nseg_b) )then
                            seg_cnt_b(nseg_b) = seg_cnt_b(nseg_b) + 1
                            cycle
                        endif
                    endif
                    nseg_b = nseg_b + 1
                    seg_dir_b(nseg_b) = dirof_pb(row)
                    seg_beg_b(nseg_b) = i
                    seg_cnt_b(nseg_b) = 1
                end do
                !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
                !$omp& private(iseg,ii,i,ithr,q,r,a,aa,e_mm,myv,row,twp0,twp1,twp2,nyqb,ldA,lok,qml)
                do iseg = 1, nseg_b
                    ithr = omp_get_thread_num() + 1
                    call dirs_pb%get_ori(seg_dir_b(iseg), o_pb_thr(ithr))
                    twp0 = omp_get_wtime()
                    call project_fplanes_mean_basis(mean_rec, basis_recs, o_pb_thr(ithr), &
                        &fpls(seg_beg_b(iseg)), mean_fpl(ithr), basis_fpls(:,ithr), &
                        &apply_ctf_amp=.false.)
                    twp1 = omp_get_wtime()
                    sec_proj_thr(ithr) = sec_proj_thr(ithr) + (twp1 - twp0)
                    do ii = 1, seg_cnt_b(iseg)
                        i = seg_beg_b(iseg) + ii - 1
                        if( orientations(i)%isstatezero() ) cycle
                        row  = batchlims(1) + i - 1
                        twp1 = omp_get_wtime()
                        nyqb = max(1, fpls(i)%nyq / OSMPL_PAD_FAC)
                        call align_halfplane_inplane(fpls(i)%frlims, nyqb, fpls(i)%cmplx_plane, &
                            &cap(row), sap(row), tmpc(:,:,ithr))
                        fpls(i)%cmplx_plane = tmpc(:,:,ithr)
                        call align_halfplane_inplane(fpls(i)%frlims, nyqb, fpls(i)%transfer_plane, &
                            &cap(row), sap(row), tmpc(:,:,ithr))
                        fpls(i)%transfer_plane = tmpc(:,:,ithr)
                        fpls(i)%ctfsq_plane = real(fpls(i)%transfer_plane*conjg(fpls(i)%transfer_plane))
                        call ensure_bank_ctf_planes(fpls(i), mean_fplC(ithr), basis_fplsC(:,ithr), ncomp)
                        mean_fplC(ithr)%cmplx_plane = fpls(i)%transfer_plane * mean_fpl(ithr)%cmplx_plane
                        do q = 1, ncomp
                            basis_fplsC(q,ithr)%cmplx_plane = fpls(i)%transfer_plane * &
                                &basis_fpls(q,ithr)%cmplx_plane
                        end do
                        e_mm = real(cov_herm_inner(mean_fplC(ithr), mean_fplC(ithr)), dp)
                        myv  = real(cov_herm_inner(mean_fplC(ithr), fpls(i)), dp)
                        a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
                        aa   = a*a
                        do q = 1, ncomp
                            bth(q,ithr) = real(cov_herm_inner(basis_fplsC(q,ithr), fpls(i)), dp)
                            cth(q,ithr) = real(cov_herm_inner(basis_fplsC(q,ithr), mean_fplC(ithr)), dp)
                            do r = q, ncomp
                                Gth(q,r,ithr) = real(cov_herm_inner(basis_fplsC(q,ithr), &
                                    &basis_fplsC(r,ithr)), dp)
                                Gth(r,q,ithr) = Gth(q,r,ithr)
                            end do
                        end do
                        call probe_solve_ecm(ncomp, Gth(:,:,ithr), bth(:,ithr), cth(:,ithr), &
                            &myv, e_mm, prior, sig2, n_probe_cm, a, zth(:,ithr), &
                            &Ainvth(:,:,ithr), ldA, lok, qml)
                        aa = a*a
                        if( lok ) nll_thr(ithr) = nll_thr(ithr) + ldA - qml
                        z(row,:)          = zth(:,ithr)
                        zbatch(:,i)       = zth(:,ithr)
                        do r = 1, ncomp
                            do q = 1, ncomp
                                dens(q,r,i) = zth(q,ithr)*zth(r,ithr) + Ainvth(q,r,ithr)
                            end do
                        end do
                        do q = 1, ncomp
                            gam_thr(q,ithr) = gam_thr(q,ithr) + dens(q,q,i)
                        end do
                        if( l_probe_mls )then
                            zbatch(:,i) = a*zbatch(:,i)
                            dens(:,:,i) = aa*dens(:,:,i)
                        endif
                        nval_thr(ithr)    = nval_thr(ithr) + 1
                        valid(i)          = .true.
                        ! residual in the BANK frame (transfer/ctfsq already aligned above)
                        fpls(i)%cmplx_plane = fpls(i)%cmplx_plane - real(a)*mean_fplC(ithr)%cmplx_plane
                        twp2 = omp_get_wtime()
                        sec_gram_thr(ithr) = sec_gram_thr(ithr) + (twp2 - twp1)
                    end do
                end do
                !$omp end parallel do
                else
                !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
                !$omp& private(i,ithr,q,r,a,aa,e_mm,myv,row,twp0,twp1,twp2,ldA,lok,qml,kk2,lwm,wsm) &
                !$omp& private(idp_es,pw_es,cnt_es,taz_es,l_chk_p,ichk,emmc_chk,myvc_chk,ac_chk) &
                !$omp& private(lda_chk,qml_chk,lok_chk,gnum_chk,gden_chk,bnum_chk,bden_chk) &
                !$omp& private(ir,rmx_chk,tazx_chk,pwx_chk,cntx_chk,emmx_chk,myvx_chk,ax_chk) &
                !$omp& private(dcy_chk,dcm_chk,emmd_chk,myvd_chk,ad_chk)
                do i = 1, batchsz
                    if( orientations(i)%isstatezero() ) cycle
                    ithr = omp_get_thread_num() + 1
                    row  = batchlims(1) + i - 1
                    twp0 = omp_get_wtime()
                    if( l_em_polar )then
                        ! statistics already formed in the polar quadrature; the mean plane is still
                        ! needed, in Cartesian, because the M-step backprojects y - a*T*mu
                        call project_fplane_mean(mean_rec, orientations(i), fpls(i), &
                            &mean_fpl(ithr), apply_ctf_amp=.true.)
                        a  = pcon(row)
                        aa = a*a
                        e_mm = 0.d0; myv = 0.d0   ! unused at nml=0; never pass uninitialised
                        Gth(:,:,ithr) = pGc(:,:,row)
                        bth(:,ithr)   = pbc(:,row)
                        cth(:,ithr)   = pcc(:,row)
                        twp1 = omp_get_wtime()
                        sec_proj_thr(ithr) = sec_proj_thr(ithr) + (twp1 - twp0)
                    elseif( l_pol_es )then
                        ! ---- POLAR SHARED-DIRECTION FORMER (stage 1) ----
                        idp_es  = dir_es(row)
                        l_chk_p = l_chk_act .and. row <= nchk_max
                        if( l_chk_p )then
                            ! parity check: the production Cartesian former on this particle; its
                            ! mean plane doubles as the residual's mean projection below
                            call project_fplanes_mean_basis(mean_rec, basis_recs, orientations(i), &
                                &fpls(i), mean_fpl(ithr), basis_fpls(:,ithr), apply_ctf_amp=.true.)
                            emmc_chk = real(cov_herm_inner(mean_fpl(ithr), mean_fpl(ithr)), dp)
                            myvc_chk = real(cov_herm_inner(mean_fpl(ithr), fpls(i)), dp)
                            do q = 1, ncomp
                                bc_chk(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), fpls(i)), dp)
                                cc_chk(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), mean_fpl(ithr)), dp)
                                do r = q, ncomp
                                    Gc_chk(q,r,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), &
                                        &basis_fpls(r,ithr)), dp)
                                    Gc_chk(r,q,ithr) = Gc_chk(q,r,ithr)
                                end do
                            end do
                        else
                            ! mean only, in Cartesian: the M-step backprojects y - a*(T mu) and
                            ! there is no polar->volume adjoint. Banded variant: identical
                            ! interpolation, none of project_fplane's per-call full-plane
                            ! zero-fill + ctfsq/transfer copies (measured as the bulk of the
                            ! polar project bucket at the native padded plane).
                            call project_fplane_mean_banded(mean_rec, orientations(i), fpls(i), &
                                &mean_fpl(ithr))
                        endif
                        ! polar-sample the prepped data plane ONCE at (bank direction, relative
                        ! in-plane angle): CTF amplitude, shift phase and per-shell whitening all
                        ! ride in from the same cmplx/transfer planes the Cartesian former reads.
                        ! Fused sampler: one KB geometry per ring sample shared by both plane
                        ! gathers, packed output written in place, no per-call allocations --
                        ! bit-identical statistics (see polar_sample_particle_fused).
                        call polar_sample_particle_fused(fpls(i)%cmplx_plane, fpls(i)%transfer_plane, &
                            &pg_es, cae(row), sae(row), xws_es(:,ithr), wr_es(:,ithr), taz_es)
                        wrd_es(:,ithr) = real(wr_es(:,ithr), dp)
                        ! b and the mean row, exact per-sample CTF: one GEMV against the bank
                        call sgemv('T', nsamp2_es, ncomp+1, 1.0, UsallE(1,0,idp_es), nsamp2_es, &
                            &xws_es(1,ithr), 1, 0.0, Reb_es(0,ithr), 1)
                        ! G, c, e_mm by the radial factorisation: ring Grams x per-ring mean |T|^2
                        call dgemv('N', ncomp*ncomp, nk_es, 1.d0, CfE(1,1,idp_es), ncomp*ncomp, &
                            &wrd_es(1,ithr), 1, 0.d0, Gth(1,1,ithr), 1)
                        call dgemv('N', ncomp, nk_es, 1.d0, Cm0E(1,1,idp_es), ncomp, &
                            &wrd_es(1,ithr), 1, 0.d0, cth(1,ithr), 1)
                        e_mm = dot_product(c00E(:,idp_es), wrd_es(:,ithr))
                        myv  = real(Reb_es(0,ithr), dp)
                        do q = 1, ncomp
                            bth(q,ithr) = real(Reb_es(q,ithr), dp)
                        end do
                        ! hybrid: the low-k shells enter as exact Cartesian statistics
                        if( l_pol_hyb ) call polar_hybrid_exact_accum(mean_rec, basis_recs, &
                            &ncomp, orientations(i), fpls(i), hex_es, kex_es, npos_es, &
                            &Gth(:,:,ithr), bth(:,ithr), cth(:,ithr), e_mm, myv)
                        a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
                        aa   = a*a
                        twp1 = omp_get_wtime()
                        sec_proj_thr(ithr) = sec_proj_thr(ithr) + (twp1 - twp0)
                        if( l_chk_p )then
                            ! Cartesian z at this particle (plain fixed-a solve; iteration 1 only,
                            ! where the mixture is never active, so both paths solve the same system)
                            ac_chk = max(0.1d0, min(5.0d0, myvc_chk / max(emmc_chk, DTINY)))
                            call probe_solve_ecm(ncomp, Gc_chk(:,:,ithr), bc_chk(:,ithr), &
                                &cc_chk(:,ithr), myvc_chk, emmc_chk, prior, sig2, 0, ac_chk, &
                                &zc_chk(:,ithr), Ai_chk(:,:,ithr), lda_chk, lok_chk, qml_chk)
                            ichk = row
                            zchk_c(:,ichk) = zc_chk(:,ithr)
                            gnum_chk = 0.d0; gden_chk = 0.d0; bnum_chk = 0.d0; bden_chk = 0.d0
                            do r = 1, ncomp
                                bnum_chk = bnum_chk + (bth(r,ithr) - bc_chk(r,ithr))**2
                                bden_chk = bden_chk + bc_chk(r,ithr)**2
                                do q = 1, ncomp
                                    gnum_chk = gnum_chk + (Gth(q,r,ithr) - Gc_chk(q,r,ithr))**2
                                    gden_chk = gden_chk + Gc_chk(q,r,ithr)**2
                                end do
                            end do
                            gerr_chk(ichk) = sqrt(gnum_chk / max(gden_chk, DTINY))
                            berr_chk(ichk) = sqrt(bnum_chk / max(bden_chk, DTINY))
                            lchk(ichk)     = .true.
                            ! exact-direction polar arm: the same quadrature at the particle's OWN
                            ! direction (identity in-plane), separating the direction snap from the
                            ! polar formulation itself
                            rmx_chk = orientations(i)%get_mat()
                            call polar_project_recs(mean_rec, basis_recs, ncomp, rmx_chk, pg_es, &
                                &UbankE(:,:,ithr))
                            do q = 0, ncomp
                                do r = 1, nsamp_es
                                    UsX_chk(2*r-1,q,ithr) = pg_es%sqwq(r)*real (UbankE(r,q,ithr))
                                    UsX_chk(2*r,  q,ithr) = pg_es%sqwq(r)*aimag(UbankE(r,q,ithr))
                                end do
                            end do
                            do ir = 1, nk_es
                                call polar_ring_gram(UsX_chk(1,0,ithr), nsamp2_es, ncomp, &
                                    &pg_es%rbeg(ir), pg_es%rend(ir)-pg_es%rbeg(ir)+1, &
                                    &CspE(0,0,ithr), CfX_chk(1,ir,ithr), Cm0X_chk(1,ir,ithr))
                                c00X_chk(ir,ithr) = polar_ring_selfpower(UsX_chk(1,0,ithr), &
                                    &nsamp2_es, pg_es%rbeg(ir), pg_es%rend(ir)-pg_es%rbeg(ir)+1)
                            end do
                            call polar_sample_particle_fused(fpls(i)%cmplx_plane, &
                                &fpls(i)%transfer_plane, pg_es, 1.0, 0.0, xwsX_chk(:,ithr), &
                                &wrX_chk(:,ithr), tazx_chk)
                            wrdX_chk(:,ithr) = real(wrX_chk(:,ithr), dp)
                            call sgemv('T', nsamp2_es, ncomp+1, 1.0, UsX_chk(1,0,ithr), nsamp2_es, &
                                &xwsX_chk(1,ithr), 1, 0.0, RebX_chk(0,ithr), 1)
                            call dgemv('N', ncomp*ncomp, nk_es, 1.d0, CfX_chk(1,1,ithr), ncomp*ncomp, &
                                &wrdX_chk(1,ithr), 1, 0.d0, GX_chk(1,1,ithr), 1)
                            call dgemv('N', ncomp, nk_es, 1.d0, Cm0X_chk(1,1,ithr), ncomp, &
                                &wrdX_chk(1,ithr), 1, 0.d0, cX_chk(1,ithr), 1)
                            emmx_chk = dot_product(c00X_chk(:,ithr), wrdX_chk(:,ithr))
                            myvx_chk = real(RebX_chk(0,ithr), dp)
                            do q = 1, ncomp
                                bX_chk(q,ithr) = real(RebX_chk(q,ithr), dp)
                            end do
                            if( l_pol_hyb ) call polar_hybrid_exact_accum(mean_rec, basis_recs, &
                                &ncomp, orientations(i), fpls(i), hex_es, kex_es, npos_es, &
                                &GX_chk(:,:,ithr), bX_chk(:,ithr), cX_chk(:,ithr), emmx_chk, myvx_chk)
                            ax_chk = max(0.1d0, min(5.0d0, myvx_chk / max(emmx_chk, DTINY)))
                            call probe_solve_ecm(ncomp, GX_chk(:,:,ithr), bX_chk(:,ithr), &
                                &cX_chk(:,ithr), myvx_chk, emmx_chk, prior, sig2, 0, ax_chk, &
                                &zx_chk(:,ithr), Ai_chk(:,:,ithr), lda_chk, lok_chk, qml_chk)
                            zchk_x(:,ichk) = zx_chk(:,ithr)
                            gnum_chk = 0.d0; gden_chk = 0.d0; bnum_chk = 0.d0; bden_chk = 0.d0
                            do r = 1, ncomp
                                bnum_chk = bnum_chk + (bX_chk(r,ithr) - bc_chk(r,ithr))**2
                                bden_chk = bden_chk + bc_chk(r,ithr)**2
                                do q = 1, ncomp
                                    gnum_chk = gnum_chk + (GX_chk(q,r,ithr) - Gc_chk(q,r,ithr))**2
                                    gden_chk = gden_chk + Gc_chk(q,r,ithr)**2
                                end do
                            end do
                            gerrx_chk(ichk)  = sqrt(gnum_chk / max(gden_chk, DTINY))
                            berrx_chk(ichk)  = sqrt(bnum_chk / max(bden_chk, DTINY))
                            tazsum_thr(ithr) = tazsum_thr(ithr) + real(tazx_chk, dp)
                            tazcnt_thr(ithr) = tazcnt_thr(ithr) + 1
                            ! DC-less Cartesian arm: subtract the single (0,0) sample from every
                            ! Cartesian statistic (the polar quadrature has no r=0 sample). Exact-
                            ! polar agreement with THIS reference convicts or clears the DC term.
                            dcy_chk = fpls(i)%cmplx_plane(0,0)
                            dcm_chk = mean_fpl(ithr)%cmplx_plane(0,0)
                            do q = 1, ncomp
                                dcb_chk(q,ithr) = basis_fpls(q,ithr)%cmplx_plane(0,0)
                            end do
                            emmd_chk = emmc_chk - real(conjg(dcm_chk)*dcm_chk, dp)
                            myvd_chk = myvc_chk - real(conjg(dcm_chk)*dcy_chk, dp)
                            do q = 1, ncomp
                                bd_chk(q,ithr) = bc_chk(q,ithr) - real(conjg(dcb_chk(q,ithr))*dcy_chk, dp)
                                cd_chk(q,ithr) = cc_chk(q,ithr) - real(conjg(dcb_chk(q,ithr))*dcm_chk, dp)
                                do r = q, ncomp
                                    Gd_chk(q,r,ithr) = Gc_chk(q,r,ithr) &
                                        &- real(conjg(dcb_chk(q,ithr))*dcb_chk(r,ithr), dp)
                                    Gd_chk(r,q,ithr) = Gd_chk(q,r,ithr)
                                end do
                            end do
                            ad_chk = max(0.1d0, min(5.0d0, myvd_chk / max(emmd_chk, DTINY)))
                            call probe_solve_ecm(ncomp, Gd_chk(:,:,ithr), bd_chk(:,ithr), &
                                &cd_chk(:,ithr), myvd_chk, emmd_chk, prior, sig2, 0, ad_chk, &
                                &zd_chk(:,ithr), Ai_chk(:,:,ithr), lda_chk, lok_chk, qml_chk)
                            zchk_d(:,ichk) = zd_chk(:,ithr)
                            ! exact-direction polar vs DC-less Cartesian statistic errors
                            gnum_chk = 0.d0; gden_chk = 0.d0; bnum_chk = 0.d0; bden_chk = 0.d0
                            do r = 1, ncomp
                                bnum_chk = bnum_chk + (bX_chk(r,ithr) - bd_chk(r,ithr))**2
                                bden_chk = bden_chk + bd_chk(r,ithr)**2
                                do q = 1, ncomp
                                    gnum_chk = gnum_chk + (GX_chk(q,r,ithr) - Gd_chk(q,r,ithr))**2
                                    gden_chk = gden_chk + Gd_chk(q,r,ithr)**2
                                end do
                            end do
                            gerrd_chk(ichk) = sqrt(gnum_chk / max(gden_chk, DTINY))
                            berrd_chk(ichk) = sqrt(bnum_chk / max(bden_chk, DTINY))
                        endif
                    else
                    call project_fplanes_mean_basis(mean_rec, basis_recs, orientations(i), fpls(i), &
                        &mean_fpl(ithr), basis_fpls(:,ithr), apply_ctf_amp=.true.)
                    twp1 = omp_get_wtime()
                    sec_proj_thr(ithr) = sec_proj_thr(ithr) + (twp1 - twp0)
                    e_mm = real(cov_herm_inner(mean_fpl(ithr), mean_fpl(ithr)), dp)
                    myv  = real(cov_herm_inner(mean_fpl(ithr), fpls(i)), dp)
                    a    = max(0.1d0, min(5.0d0, myv / max(e_mm, DTINY)))
                    aa   = a*a
                    do q = 1, ncomp
                        bth(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), fpls(i)), dp)
                        cth(q,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), mean_fpl(ithr)), dp)
                        do r = q, ncomp
                            Gth(q,r,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), basis_fpls(r,ithr)), dp)
                            Gth(r,q,ithr) = Gth(q,r,ithr)
                        end do
                    end do
                    endif
                    ! Posterior precision A = (a^2/sig2) G + Gamma^-1. The whole normal system is already
                    ! scaled by 1/sig2, so Cov[z|y] = A^-1 exactly -- no further sig2 factor.
                    if( l_mix_used )then
                        ! ---- MCFA E-step via the ONE shared solver (probe_solve_mix) ----
                        call probe_solve_mix(ncomp, kmix_pr, Gth(:,:,ithr), bth(:,ithr), &
                            &cth(:,ithr), myv, e_mm, nml_plain, a, sig2, mix_Ominv, mix_Omxi, &
                            &mix_lpi, mix_xiOx, zth(:,ithr), Ainvth(:,:,ithr), dens(:,:,i), &
                            &ldA, lok, nll_mix_add, mxa_sr(:,ithr), mxa_sm(:,:,ithr), &
                            &mxa_smm(:,:,:,ithr), mxa_sainv(:,:,ithr))
                        if( lok ) nll_thr(ithr) = nll_thr(ithr) + nll_mix_add
                    else
                        call probe_solve_ecm(ncomp, Gth(:,:,ithr), bth(:,ithr), cth(:,ithr), &
                            &myv, e_mm, prior, sig2, nml_plain, a, zth(:,ithr), &
                            &Ainvth(:,:,ithr), ldA, lok, qml)
                        if( lok ) nll_thr(ithr) = nll_thr(ithr) + ldA - qml
                    endif
                    aa = a*a
                    z(row,:)          = zth(:,ithr)
                    zbatch(:,i)       = zth(:,ithr)
                    ! parity check: the polar-path z actually used, captured after the shared solve
                    if( l_chk_act )then
                        if( row <= nchk_max ) zchk_p(:,row) = zth(:,ithr)
                    endif
                    ! EM sufficient statistic E[z z'|y] = z z' + Cov[z|y]. BOTH the coupled M-step normal
                    ! matrix and the Gamma update below need it. Dropping Cov underestimates Gamma, which
                    ! tightens the prior, which shrinks z further: the bias compounds across iterations.
                    ! Under the mixture E[zz'|y] = A^-1 + sum_k r_k m_k m_k', NOT A^-1 + E[z]E[z]'
                    ! -- the between-component spread is real posterior variance; probe_solve_mix
                    ! has already written dens(:,:,i) in that case.
                    if( .not. l_mix_used )then
                        do r = 1, ncomp
                            do q = 1, ncomp
                                dens(q,r,i) = zth(q,ithr)*zth(r,ithr) + Ainvth(q,r,ithr)
                            end do
                        end do
                    endif
                    do q = 1, ncomp
                        gam_thr(q,ithr) = gam_thr(q,ithr) + dens(q,q,i)
                    end do
                    if( l_probe_mls )then
                        zbatch(:,i) = a*zbatch(:,i)
                        dens(:,:,i) = aa*dens(:,:,i)
                    endif
                    nval_thr(ithr)    = nval_thr(ithr) + 1
                    valid(i)          = .true.
                    gam_dbg(1,ithr) = gam_dbg(1,ithr) + sum([(Gth(q,q,ithr), q=1,ncomp)])
                    gam_dbg(2,ithr) = gam_dbg(2,ithr) + dot_product(bth(:,ithr), bth(:,ithr))
                    gam_dbg(3,ithr) = gam_dbg(3,ithr) + dot_product(cth(:,ithr), cth(:,ithr))
                    gam_dbg(4,ithr) = gam_dbg(4,ithr) + a
                    ! residual observation r_i = y - a*(T mu) in place (transfer/ctfsq intact for backprojection)
                    if( l_pol_es )then
                        ! banded: the mean plane is zero outside the working disc (both formers
                        ! write the same disc), so the full-array statement only rewrote
                        ! unchanged values there -- at the native padded lattice that traffic
                        ! was most of the polar path's gram+solve bucket
                        call subtract_mean_banded(fpls(i), mean_fpl(ithr), real(a), nyqr_es)
                    else
                        fpls(i)%cmplx_plane = fpls(i)%cmplx_plane - real(a)*mean_fpl(ithr)%cmplx_plane
                    endif
                    twp2 = omp_get_wtime()
                    sec_gram_thr(ithr) = sec_gram_thr(ithr) + (twp2 - twp1)
                end do
                !$omp end parallel do
                endif
                sec_estep = sec_estep + toc(t_sec)
                t_sec = tic()
                ! M-step by halfset: Y_q += sum_i z_iq * backproject(r_i), and the coupled normal matrix
                ! rho(q,r) += sum_i |CTF|^2 E[z_iq z_ir]   (batched KB)
                do i = 1, batchsz
                    valid_e(i) = valid(i) .and. eo(i)==0
                    valid_o(i) = valid(i) .and. eo(i)==1
                end do
                if( l_gpu_probe )then
                    ! halfset routing through the slots: even records fill 1..ncomp / 1..npairs,
                    ! odd records the upper halves; density slots honor the diag/full mode
                    gdsc(:,:batchsz) = 0.d0
                    grsc(:,:batchsz) = 0.d0
                    do i = 1, batchsz
                        if( .not. (valid_e(i) .or. valid_o(i)) ) cycle
                        gq = merge(0, ncomp,  valid_e(i))
                        gr = merge(0, npairs, valid_e(i))
                        gdsc(gq+1:gq+ncomp,i) = zbatch(:,i)
                        do r = 1, ncomp
                            do q = 1, r
                                grsc(gr+(r*(r-1))/2+q,i) = dens(q,r,i)
                            end do
                        end do
                    end do
                    if( l_bank_adj )then
                        if( l_bank_fwd )then
                            ! residuals were built in the bank frame by the forward -- unit rotation
                            call flex_gpu_coupled_batch_banked_f(fpls(:batchsz), gdsc(:,:batchsz), &
                                &grsc(:,:batchsz), valid_e(:batchsz) .or. valid_o(:batchsz), batchsz, &
                                &dirof_pb(batchlims(1):batchlims(2)), cap1(:batchsz), sap0(:batchsz))
                        else if( l_gpu_estep )then
                            ! fused path: the residual is already resident on device
                            call flex_gpu_coupled_batch_banked_res_f(gdsc(:,:batchsz), &
                                &grsc(:,:batchsz), valid_e(:batchsz) .or. valid_o(:batchsz), batchsz, &
                                &dirof_pb(batchlims(1):batchlims(2)), cap(batchlims(1):batchlims(2)), &
                                &sap(batchlims(1):batchlims(2)))
                        else
                            call flex_gpu_coupled_batch_banked_f(fpls(:batchsz), gdsc(:,:batchsz), &
                                &grsc(:,:batchsz), valid_e(:batchsz) .or. valid_o(:batchsz), batchsz, &
                                &dirof_pb(batchlims(1):batchlims(2)), cap(batchlims(1):batchlims(2)), &
                                &sap(batchlims(1):batchlims(2)))
                        endif
                    else
                        call flex_gpu_coupled_batch_raw_f(build%pgrpsyms, orientations(:batchsz), &
                            &fpls(:batchsz), gdsc(:,:batchsz), grsc(:,:batchsz), &
                            &valid_e(:batchsz) .or. valid_o(:batchsz), batchsz)
                    endif
                else
                    call insert_planes_oversamp_coupled_batch_scaled(Yeven, rho_e, build%pgrpsyms, &
                        &orientations(:batchsz), fpls(:batchsz), zbatch(:,:batchsz), dens(:,:,:batchsz), &
                        &valid_e(:batchsz), batchsz)
                    call insert_planes_oversamp_coupled_batch_scaled(Yodd, rho_o, build%pgrpsyms, &
                        &orientations(:batchsz), fpls(:batchsz), zbatch(:,:batchsz), dens(:,:,:batchsz), &
                        &valid_o(:batchsz), batchsz)
                endif
                sec_ins = sec_ins + toc(t_sec)
                if( batchlims(2)==nptcls .or. mod(batchlims(2), 5*MAXIMGBATCHSZ)==0 )then
                    write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA PROBE PASS PARTICLES: ',batchlims(2),' / ',npp
                    call flush(logfhandle)
                endif
            end do
            if( l_gpu_probe )then
                call flex_gpu_coupled_end_f(Yeven, rho_e, Yodd, rho_o)
                deallocate(gdsc, grsc)
            endif
            write(logfhandle,'(A,F7.1,A,F7.1,A,F7.1,A,F7.1)') '>>> FLEX_PCA PROBE E-STEP SPLIT (seconds): read=', &
                &sec_read,'  prep=',sec_prep,'  project+solve=',sec_estep,'  insert=',sec_ins
            write(logfhandle,'(A,F7.1,A,F7.1)') '>>> FLEX_PCA PROBE E-STEP INNER (thread-seconds): project=', &
                &sum(sec_proj_thr),'  gram+solve=',sum(sec_gram_thr)
            call flush(logfhandle)
            if( l_pol_es )then
                ! the two numbers the polar A/B reads from run.log
                write(logfhandle,'(A,I0,A,F7.1,A,F7.1)') '>>> FLEX_PCA POLAR ESTEP it=',it, &
                    &'  bank build seconds=',sec_bank,'  estep seconds=',sec_estep
                call flush(logfhandle)
            endif
            ! ---- POLAR CHECK SUMMARY (iteration 1): three formers ran on the first N particles:
            ! Cartesian (reference), polar at the bank direction (the production stage-1 path) and
            ! polar at the exact direction. banked-vs-exact isolates the direction quantization;
            ! exact-vs-Cartesian isolates the polar formulation (quadrature + radial |T|^2). ----
            if( l_chk_act )then
                block
                    integer  :: nck, qq2, j2
                    real(dp), allocatable :: va_ck(:), vb_ck(:), vx_ck(:)
                    real,     allocatable :: se_ck(:)
                    real(dp) :: cq_ck(ncomp), cqx_ck(ncomp), cqbx_ck(ncomp)
                    real(dp) :: gmed_ck, bmed_ck, gmedx_ck, bmedx_ck, taz_ck
                    nck = count(lchk)
                    if( nck >= 3 )then
                        allocate(va_ck(nck), vb_ck(nck), vx_ck(nck), se_ck(nck))
                        do qq2 = 1, ncomp
                            j2 = 0
                            do ichk = 1, nchk_max
                                if( .not. lchk(ichk) ) cycle
                                j2 = j2 + 1
                                va_ck(j2) = zchk_c(qq2,ichk)
                                vb_ck(j2) = zchk_p(qq2,ichk)
                                vx_ck(j2) = 0.d0
                                if( allocated(zchk_x) ) vx_ck(j2) = zchk_x(qq2,ichk)
                            end do
                            cq_ck(qq2)   = corr_dp(va_ck, vb_ck, nck)
                            cqx_ck(qq2)  = 0.d0
                            cqbx_ck(qq2) = 0.d0
                            if( allocated(zchk_x) )then
                                cqx_ck(qq2)  = corr_dp(va_ck, vx_ck, nck)
                                cqbx_ck(qq2) = corr_dp(vb_ck, vx_ck, nck)
                            endif
                        end do
                        j2 = 0
                        do ichk = 1, nchk_max
                            if( .not. lchk(ichk) ) cycle
                            j2 = j2 + 1
                            se_ck(j2) = real(gerr_chk(ichk))
                        end do
                        call hpsort(se_ck)
                        gmed_ck = real(se_ck((nck+1)/2), dp)
                        j2 = 0
                        do ichk = 1, nchk_max
                            if( .not. lchk(ichk) ) cycle
                            j2 = j2 + 1
                            se_ck(j2) = real(berr_chk(ichk))
                        end do
                        call hpsort(se_ck)
                        bmed_ck = real(se_ck((nck+1)/2), dp)
                        gmedx_ck = 0.d0
                        bmedx_ck = 0.d0
                        if( allocated(gerrx_chk) )then
                            j2 = 0
                            do ichk = 1, nchk_max
                                if( .not. lchk(ichk) ) cycle
                                j2 = j2 + 1
                                se_ck(j2) = real(gerrx_chk(ichk))
                            end do
                            call hpsort(se_ck)
                            gmedx_ck = real(se_ck((nck+1)/2), dp)
                            j2 = 0
                            do ichk = 1, nchk_max
                                if( .not. lchk(ichk) ) cycle
                                j2 = j2 + 1
                                se_ck(j2) = real(berrx_chk(ichk))
                            end do
                            call hpsort(se_ck)
                            bmedx_ck = real(se_ck((nck+1)/2), dp)
                        endif
                        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA POLAR CHECK: Cartesian vs polar &
                            &E-step statistics on ',nck,' particles (iteration 1)'
                        if( allocated(zchk_g) )then
                            write(logfhandle,'(A)') '>>>   z correlation per component (DEVICE polar vs Cartesian):'
                        else
                            write(logfhandle,'(A)') '>>>   z correlation per component (BANKED polar vs Cartesian):'
                        endif
                        do qq2 = 1, ncomp, 8
                            write(logfhandle,'(A,I3,A,8F8.4)') '>>>     z',qq2,'..', &
                                &cq_ck(qq2:min(ncomp,qq2+7))
                        end do
                        write(logfhandle,'(A,F8.4,A,ES10.2,A,ES10.2)') '>>>   BANKED vs Cartesian:   min z corr=', &
                            &minval(cq_ck),'  median rel |G-Gc|_F/|Gc|_F=',gmed_ck, &
                            &'  median rel |b-bc|/|bc|=',bmed_ck
                        ! ---- GPU parity summary (gate 1): device former vs the CPU polar former,
                        ! identical bank/tables/hybrid and identical solver settings ----
                        if( allocated(zchk_g) )then
                            block
                                real(dp) :: cqg_ck(ncomp), gmedg_ck, bmedg_ck
                                do qq2 = 1, ncomp
                                    j2 = 0
                                    do ichk = 1, nchk_max
                                        if( .not. lchk(ichk) ) cycle
                                        j2 = j2 + 1
                                        va_ck(j2) = zchk_g(qq2,ichk)
                                        vb_ck(j2) = zchk_p(qq2,ichk)
                                    end do
                                    cqg_ck(qq2) = corr_dp(va_ck, vb_ck, nck)
                                end do
                                j2 = 0
                                do ichk = 1, nchk_max
                                    if( .not. lchk(ichk) ) cycle
                                    j2 = j2 + 1
                                    se_ck(j2) = real(gerrg_chk(ichk))
                                end do
                                call hpsort(se_ck)
                                gmedg_ck = real(se_ck((nck+1)/2), dp)
                                j2 = 0
                                do ichk = 1, nchk_max
                                    if( .not. lchk(ichk) ) cycle
                                    j2 = j2 + 1
                                    se_ck(j2) = real(berrg_chk(ichk))
                                end do
                                call hpsort(se_ck)
                                bmedg_ck = real(se_ck((nck+1)/2), dp)
                                write(logfhandle,'(A)') '>>>   z correlation per component (DEVICE &
                                    &polar vs CPU polar):'
                                do qq2 = 1, ncomp, 8
                                    write(logfhandle,'(A,I3,A,8F8.4)') '>>>     z',qq2,'..', &
                                        &cqg_ck(qq2:min(ncomp,qq2+7))
                                end do
                                write(logfhandle,'(A,F8.4,A,F8.4,A,ES10.2,A,ES10.2)') &
                                    &'>>>   DEVICE vs CPU-POLAR:   min z corr=',minval(cqg_ck), &
                                    &'  mean=',sum(cqg_ck)/real(ncomp,dp), &
                                    &'  median rel |G-Gp|_F/|Gp|_F=',gmedg_ck, &
                                    &'  median rel |b-bp|/|bp|=',bmedg_ck
                            end block
                        endif
                        if( allocated(zchk_x) )then
                        write(logfhandle,'(A)') '>>>   z correlation per component (EXACT-DIR polar vs Cartesian):'
                        do qq2 = 1, ncomp, 8
                            write(logfhandle,'(A,I3,A,8F8.4)') '>>>     z',qq2,'..', &
                                &cqx_ck(qq2:min(ncomp,qq2+7))
                        end do
                        write(logfhandle,'(A,F8.4,A,ES10.2,A,ES10.2)') '>>>   EXACT-DIR vs Cartesian: min z corr=', &
                            &minval(cqx_ck),'  median rel |G-Gc|_F/|Gc|_F=',gmedx_ck, &
                            &'  median rel |b-bc|/|bc|=',bmedx_ck
                        write(logfhandle,'(A,F8.4,A,F8.4)') '>>>   BANKED vs EXACT-DIR (direction &
                            &quantization alone): min z corr=',minval(cqbx_ck),'  mean=', &
                            &sum(cqbx_ck)/real(ncomp,dp)
                        taz_ck = sum(tazsum_thr) / real(max(1, sum(tazcnt_thr)), dp)
                        write(logfhandle,'(A,F8.4,A)') '>>>   mean within-ring |T|^2 spread (tazim)=', &
                            &taz_ck,'  (bounds the radial-factorisation error in G)'
                        endif
                        ! DC attribution: exact-direction polar against the DC-less Cartesian
                        ! (CPU attribution arm only; stands down under the device path)
                        if( allocated(zchk_x) )then
                        block
                            real(dp) :: cqd_ck(ncomp), cqcd_ck(ncomp), gdm_ck, bdm_ck
                            do qq2 = 1, ncomp
                                j2 = 0
                                do ichk = 1, nchk_max
                                    if( .not. lchk(ichk) ) cycle
                                    j2 = j2 + 1
                                    va_ck(j2) = zchk_d(qq2,ichk)
                                    vx_ck(j2) = zchk_x(qq2,ichk)
                                    vb_ck(j2) = zchk_c(qq2,ichk)
                                end do
                                cqd_ck(qq2)  = corr_dp(va_ck, vx_ck, nck)   ! exact polar vs DC-less cart
                                cqcd_ck(qq2) = corr_dp(va_ck, vb_ck, nck)   ! DC-less cart vs full cart
                            end do
                            j2 = 0
                            do ichk = 1, nchk_max
                                if( .not. lchk(ichk) ) cycle
                                j2 = j2 + 1
                                se_ck(j2) = real(gerrd_chk(ichk))
                            end do
                            call hpsort(se_ck)
                            gdm_ck = real(se_ck((nck+1)/2), dp)
                            j2 = 0
                            do ichk = 1, nchk_max
                                if( .not. lchk(ichk) ) cycle
                                j2 = j2 + 1
                                se_ck(j2) = real(berrd_chk(ichk))
                            end do
                            call hpsort(se_ck)
                            bdm_ck = real(se_ck((nck+1)/2), dp)
                            write(logfhandle,'(A)') '>>>   z correlation per component (EXACT-DIR polar &
                                &vs DC-LESS Cartesian):'
                            do qq2 = 1, ncomp, 8
                                write(logfhandle,'(A,I3,A,8F8.4)') '>>>     z',qq2,'..', &
                                    &cqd_ck(qq2:min(ncomp,qq2+7))
                            end do
                            write(logfhandle,'(A,F8.4,A,ES10.2,A,ES10.2)') &
                                &'>>>   EXACT-DIR vs DC-LESS Cartesian: min z corr=',minval(cqd_ck), &
                                &'  median rel G err=',gdm_ck,'  median rel b err=',bdm_ck
                            write(logfhandle,'(A,F8.4,A,F8.4)') '>>>   DC-LESS vs FULL Cartesian &
                                &(the DC term alone): min z corr=',minval(cqcd_ck),'  mean=', &
                                &sum(cqcd_ck)/real(ncomp,dp)
                        end block
                        endif
                        call flush(logfhandle)
                        deallocate(va_ck, vb_ck, vx_ck, se_ck)
                    else
                        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA POLAR CHECK: only ',nck, &
                            &' particles checked -- too few for a summary'
                        call flush(logfhandle)
                    endif
                end block
                deallocate(zchk_c, zchk_p, gerr_chk, berr_chk, lchk, Gc_chk, bc_chk, cc_chk, &
                    &zc_chk, Ai_chk)
                deallocate(zchk_d, gerrd_chk, berrd_chk, Gd_chk, bd_chk, cd_chk, zd_chk, dcb_chk)
                if( allocated(zchk_x) ) deallocate(zchk_x, gerrx_chk, berrx_chk, UsX_chk, &
                    &xwsX_chk, wrX_chk, RebX_chk, CfX_chk, Cm0X_chk, c00X_chk, wrdX_chk, &
                    &GX_chk, bX_chk, cX_chk, zx_chk, tazsum_thr, tazcnt_thr)
                if( allocated(zchk_g) ) deallocate(zchk_g, gerrg_chk, berrg_chk, Gpc_chk, &
                    &bpc_chk, cpc_chk, zpc_chk)
                l_chk_act   = .false.
                n_pol_check = 0
            endif
            ! reduce the EM Gamma accumulator BEFORE ncomp is replaced by d_new below.
            ! Gamma travels between parts as a SUM and is divided by the REDUCED nval: dividing per
            ! part would weight a small part equally with a large one.
            nval = sum(nval_thr)
            ! E-STEP STATISTICS SUMMARY. The polar and Cartesian formers are supposed to produce
            ! the SAME G, b, c on the same basis; at iteration 1 the basis is the deterministic
            ! data-free init, so any difference here is the former, not the fit. Printing the
            ! actual moments is the only way to tell a scale mismatch from a genuine difference.
            if( nval > 0 )then
                write(logfhandle,'(A,I0,A,ES12.4,A,ES12.4,A,ES12.4,A,F7.4)') &
                    &'>>> FLEX_PCA PROBE ESTAT it=',it,'  <trG>=',sum(gam_dbg(1,:))/real(nval,dp), &
                    &'  <b.b>=',sum(gam_dbg(2,:))/real(nval,dp), &
                    &'  <c.c>=',sum(gam_dbg(3,:))/real(nval,dp), &
                    &'  <a>=',real(sum(gam_dbg(4,:))/real(nval,dp))
                call flush(logfhandle)
            endif
            ! In-process reduction of the likelihood accumulator. This was MISSING: nll_tot was
            ! summed only in the worker branch, so a shared-memory run (nparts=1) reported just the
            ! N*log det Gamma term and none of the per-particle part -- which silently invalidated
            ! every nparts=1 likelihood comparison. Distributed runs were unaffected, since
            ! reduce_probe_parts sums the workers' contributions.
            nll_tot = sum(nll_thr)
            do q = 1, ncomp
                gam_sum(q) = sum(gam_thr(q,:))
            end do
            ! ---- MCFA: mixture M-step, or its (re)initialisation from the current latents ----
            ! nparts>1 is excluded at the flag read, so this always runs on the in-process
            ! master with this iteration's accumulators complete.
            if( l_mix_req .and. .not. flex_pca_is_worker() )then
                if( l_mix_active )then
                    block
                        real(dp), allocatable :: rr_sr(:), rr_sm(:,:), rr_smm(:,:,:), rr_sai(:,:)
                        real(dp) :: g_mx, pi_floor, best_d, dmin, dsq
                        integer  :: tt, kk3, kkp, ii2, best_i
                        allocate(rr_sr(kmix_pr), rr_sm(ncomp,kmix_pr), &
                            &rr_smm(ncomp,ncomp,kmix_pr), rr_sai(ncomp,ncomp))
                        rr_sr = 0.d0; rr_sm = 0.d0; rr_smm = 0.d0; rr_sai = 0.d0
                        do tt = 1, nthr
                            rr_sr  = rr_sr  + mxa_sr(:,tt)
                            rr_sm  = rr_sm  + mxa_sm(:,:,tt)
                            rr_sai = rr_sai + mxa_sainv(:,:,tt)
                            do kk3 = 1, kmix_pr
                                rr_smm(:,:,kk3) = rr_smm(:,:,kk3) + mxa_smm(:,:,kk3,tt)
                            end do
                        end do
                        if( allocated(dm_sr) .and. flex_pca_nparts() > 1 )then
                            ! distributed: the reduce already summed every part's thread-sums
                            call mcfa_mstep(ncomp, kmix_pr, nval, dm_sr, dm_sm, dm_smm, dm_sai, &
                                &kmix_pr == 1, mix_pi, mix_xi, mix_Om, mix_Ominv, ldOm_mix)
                        else
                            call mcfa_mstep(ncomp, kmix_pr, nval, rr_sr, rr_sm, rr_smm, rr_sai, &
                                &kmix_pr == 1, mix_pi, mix_xi, mix_Om, mix_Ominv, ldOm_mix)
                        endif
                        deallocate(rr_sr, rr_sm, rr_smm, rr_sai)
                    end block
                else if( it >= n_mix_warm )then
                    allocate(mix_xi(ncomp,kmix_pr), mix_Om(ncomp,ncomp), mix_Ominv(ncomp,ncomp), &
                        &mix_pi(kmix_pr), mix_Omxi(ncomp,kmix_pr), mix_xiOx(kmix_pr), mix_lpi(kmix_pr))
                    if( flex_pca_nparts() > 1 .and. dm_nz > 0 )then
                        ! distributed: seed from the pooled per-part latent subsample
                        call mcfa_init(dm_z(:dm_nz,:), dm_nz, ncomp, kmix_pr, gam_sum, nval, &
                            &mix_xi, mix_pi, mix_Om)
                    else
                        call mcfa_init(z, size(z,1), ncomp, kmix_pr, gam_sum, nval, mix_xi, mix_pi, mix_Om)
                    endif
                    call mcfa_condition(ncomp, kmix_pr == 1, mix_Om, mix_Ominv, ldOm_mix)
                    l_mix_active = .true.
                    write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA MIX initialised: K=',kmix_pr, &
                        &'  at iteration ',it
                endif
                if( l_mix_active )then
                    do kk2 = 1, kmix_pr
                        mix_Omxi(:,kk2) = matmul(mix_Ominv, mix_xi(:,kk2))
                        mix_xiOx(kk2)   = dot_product(mix_xi(:,kk2), mix_Omxi(:,kk2))
                        mix_lpi(kk2)    = log(max(mix_pi(kk2), 1.d-12))
                    end do
                    write(logfhandle,'(A,I0,A,ES11.3,A,F7.4,A,F8.3)') '>>> FLEX_PCA MIX it=',it, &
                        &'  logdetOm=',ldOm_mix,'  pi_max=',maxval(mix_pi), &
                        &'  |xi|_max=',sqrt(maxval(sum(mix_xi**2, dim=1)))
                    call flush(logfhandle)
                endif
            endif
            if( flex_pca_is_worker() )then
                allocate(cme(es(1),es(2),es(3),ncomp), cmo(es(1),es(2),es(3),ncomp))
                allocate(rhe(es(1),es(2),es(3),ncomp), rhoo(es(1),es(2),es(3),ncomp))
                do q = 1, ncomp
                    cme(:,:,:,q) = Yeven(q)%cmat_exp; rhe(:,:,:,q)  = Yeven(q)%rho_exp
                    cmo(:,:,:,q) = Yodd(q)%cmat_exp;  rhoo(:,:,:,q) = Yodd(q)%rho_exp
                end do
                nll_tot = sum(nll_thr)
                if( l_mix_req )then
                    block
                        real(dp), allocatable :: w_sr(:), w_sm(:,:), w_smm(:,:,:), w_sai(:,:), w_z(:,:)
                        integer :: tt2, kk4, nzs, izs, istep
                        allocate(w_sr(kmix_pr), w_sm(ncomp,kmix_pr), &
                            &w_smm(ncomp,ncomp,kmix_pr), w_sai(ncomp,ncomp))
                        w_sr = 0.d0; w_sm = 0.d0; w_smm = 0.d0; w_sai = 0.d0
                        do tt2 = 1, nthr
                            w_sr  = w_sr  + mxa_sr(:,tt2)
                            w_sm  = w_sm  + mxa_sm(:,:,tt2)
                            w_sai = w_sai + mxa_sainv(:,:,tt2)
                            do kk4 = 1, kmix_pr
                                w_smm(:,:,kk4) = w_smm(:,:,kk4) + mxa_smm(:,:,kk4,tt2)
                            end do
                        end do
                        ! deterministic stride subsample of this part's latents (bounded)
                        nzs   = min(MIX_ZSUB_MAX, npp)
                        istep = max(1, npp / max(1,nzs))
                        nzs   = min(nzs, (npp + istep - 1)/istep)
                        allocate(w_z(nzs,ncomp))
                        do izs = 1, nzs
                            w_z(izs,:) = z(min(npp, 1 + (izs-1)*istep), :)
                        end do
                        call write_probe_part(flex_pca_part_fname('probe', params%part, params%numlen), &
                            &cme, rhe, cmo, rhoo, rho_e, rho_o, gam_sum, nll_tot, nval, ncomp, &
                            &w_sr, w_sm, w_smm, w_sai, w_z)
                        deallocate(w_sr, w_sm, w_smm, w_sai, w_z)
                    end block
                else
                    call write_probe_part(flex_pca_part_fname('probe', params%part, params%numlen), &
                        &cme, rhe, cmo, rhoo, rho_e, rho_o, gam_sum, nll_tot, nval, ncomp)
                endif
                deallocate(cme, cmo, rhe, rhoo, gam_sum)
                call cleanup_rec_buffers(build, fpls)
                do q = 1, size(Yeven)
                    call Yeven(q)%dealloc_rho; call Yeven(q)%kill
                    call Yodd(q)%dealloc_rho;  call Yodd(q)%kill
                end do
                deallocate(Yeven, Yodd, rho_e, rho_o, prior)
                deallocate(Gth, Ath, bth, cth, zth, basis_fpls, mean_fpl, orientations, zbatch, dens)
                deallocate(valid, valid_e, valid_o, eo, Ainvth, Acpth, gam_thr, gam_acc, nval_thr)
                deallocate(hth, nll_thr)
                if( l_bank_fwd )then
                    do ithr = 1, nthr
                        call cleanup_plane(mean_fplC(ithr))
                        do q = 1, size(basis_fplsC,1); call cleanup_plane(basis_fplsC(q,ithr)); end do
                    end do
                    deallocate(mean_fplC, basis_fplsC)
                endif
                if( l_bank_adj )then
                    call flex_gpu_coupled_bank_free_f
                    call flex_gpu_estep_free_f
                    if( l_gpu_prep )then
                        call flex_gpu_prep_free_f
                        deallocate(ctfp_arr, shf_pb)
                        if( allocated(sig2_pb) ) deallocate(sig2_pb)
                    endif
                    call dirs_pb%kill
                    deallocate(rmat_pb, nrm_pb, dirof_pb, cap, sap, o_pb_thr, cap1, sap0)
                    deallocate(seg_dir_b, seg_beg_b, seg_cnt_b)
                    if( allocated(tmpc) ) deallocate(tmpc)
                    if( allocated(stats_es) ) deallocate(stats_es)
                    if( allocated(avec_es) ) deallocate(avec_es, vld_es)
                endif
                deallocate(z, ppinds)
                return
            endif
            endif
            do q = 1, ncomp
                gam_acc(q) = gam_sum(q) / real(max(1,nval),dp)
            end do
            deallocate(gam_sum)
            ! ---- ONLINE EM: blend this batch's sufficient statistics into the running average,
            ! rotating the history into the current basis frame first (rationale at the flag).
            ! MEASURED 2026-08-19 on 10076: the averaging HURTS (0.1936/0.3278/ceil 0.2537 vs
            ! windows-only 0.2361/0.3509/0.2648) -- stale M-step inputs + double regularization
            ! on top of the Wiener. Windowed rotation without averaging is the winner, so
            ! SIMPLE_COV_ONLINE=1 is windows-only and the blend is opt-in at =2 (research).
            if( l_online .and. vonl_pr >= 2 )then
                block
                    real,     pointer     :: onl_pr(:,:,:)
                    real(dp), allocatable :: Mo(:,:), MtMo(:,:), Vo(:,:), evo(:)
                    real,     allocatable :: ut_now(:,:,:,:)
                    complex,  allocatable :: tcme(:,:,:,:), tcmo(:,:,:,:)
                    real,     allocatable :: trhe(:,:,:,:), trho(:,:,:,:)
                    real(dp), allocatable :: tgam(:)
                    real(dp) :: g_onl, wpq
                    integer  :: p2, q2, nrot_o
                    logical  :: l_rot_ok
                    g_onl = (2.d0 + real(it,dp))**(-0.7d0)
                    if( it == 1 .or. .not. l_oh_live ) g_onl = 1.d0
                    ! the CURRENT E-step frame lives in prev_real (utilde_real is per-iteration
                    ! scratch, killed after propagation -- it is NEVER alive at blend time)
                    l_rot_ok = .false.
                    if( allocated(prev_real) .and. l_oh_live .and. allocated(ut_hist) )then
                        if( size(prev_real) == ncomp )then
                            call prev_real(1)%get_rmat_ptr(onl_pr)
                            if( size(ut_hist,1) == size(onl_pr) ) l_rot_ok = .true.
                        endif
                    endif
                    if( g_onl < 1.d0 .and. .not. l_rot_ok )then
                        ! rank changed or first usable frame: restart the average rather than mix frames
                        g_onl = 1.d0
                        write(logfhandle,'(A)') '>>> FLEX_PCA ONLINE: history reset (frame not alignable)'
                    endif
                    if( g_onl < 1.d0 )then
                        ! Procrustes rotation hist-frame -> current frame, from the basis overlap
                        allocate(Mo(ncomp,ncomp), MtMo(ncomp,ncomp), Vo(ncomp,ncomp), evo(ncomp))
                        block
                            real(dp), allocatable :: vnow(:)
                            allocate(vnow(size(ut_hist,1)))
                            do q2 = 1, ncomp
                                call prev_real(q2)%get_rmat_ptr(onl_pr)
                                vnow = real(reshape(onl_pr, [size(onl_pr)]), dp)
                                do p2 = 1, ncomp
                                    Mo(p2,q2) = sum(ut_hist(:,p2) * vnow)
                                end do
                            end do
                        end block
                        MtMo = matmul(transpose(Mo), Mo)
                        call jacobi(MtMo, ncomp, ncomp, evo, Vo, nrot_o)
                        do q2 = 1, ncomp
                            evo(q2) = 1.d0/sqrt(max(evo(q2), 1.d-12))
                        end do
                        do q2 = 1, ncomp
                            MtMo(:,q2) = Vo(:,q2)*evo(q2)     ! scratch: V * diag(1/sqrt(ev))
                        end do
                        Mo = matmul(Mo, matmul(MtMo, transpose(Vo)))
                        ! rotate history: linear for the map sums, R^2 weights for rho/gamma
                        allocate(tcme, mold=oh_cme); allocate(tcmo, mold=oh_cmo)
                        allocate(trhe, mold=oh_rhe); allocate(trho, mold=oh_rho)
                        allocate(tgam(ncomp))
                        tcme = cmplx(0.,0.); tcmo = cmplx(0.,0.); trhe = 0.; trho = 0.; tgam = 0.d0
                        do q2 = 1, ncomp
                            do p2 = 1, ncomp
                                wpq = Mo(p2,q2)
                                if( abs(wpq) < 1.d-6 ) cycle
                                tcme(:,:,:,q2) = tcme(:,:,:,q2) + real(wpq)*oh_cme(:,:,:,p2)
                                tcmo(:,:,:,q2) = tcmo(:,:,:,q2) + real(wpq)*oh_cmo(:,:,:,p2)
                                trhe(:,:,:,q2) = trhe(:,:,:,q2) + real(wpq*wpq)*oh_rhe(:,:,:,p2)
                                trho(:,:,:,q2) = trho(:,:,:,q2) + real(wpq*wpq)*oh_rho(:,:,:,p2)
                                tgam(q2)       = tgam(q2)       + wpq*wpq*oh_gam(p2)
                            end do
                        end do
                        call move_alloc(tcme, oh_cme); call move_alloc(tcmo, oh_cmo)
                        call move_alloc(trhe, oh_rhe); call move_alloc(trho, oh_rho)
                        oh_gam = tgam
                        deallocate(tgam, Mo, MtMo, Vo, evo)
                    endif
                    if( .not. allocated(oh_cme) )then
                        associate( ce => Yeven(1)%cmat_exp )
                            allocate(oh_cme(size(ce,1),size(ce,2),size(ce,3),ncomp))
                            allocate(oh_cmo(size(ce,1),size(ce,2),size(ce,3),ncomp))
                        end associate
                        associate( re => Yeven(1)%rho_exp )
                            allocate(oh_rhe(size(re,1),size(re,2),size(re,3),ncomp))
                            allocate(oh_rho(size(re,1),size(re,2),size(re,3),ncomp))
                        end associate
                        allocate(oh_gam(ncomp))
                    endif
                    ! blend and write back so the M-step consumes the running average
                    do q2 = 1, ncomp
                        if( g_onl >= 1.d0 )then
                            oh_cme(:,:,:,q2) = Yeven(q2)%cmat_exp
                            oh_cmo(:,:,:,q2) = Yodd(q2)%cmat_exp
                            oh_rhe(:,:,:,q2) = Yeven(q2)%rho_exp
                            oh_rho(:,:,:,q2) = Yodd(q2)%rho_exp
                            oh_gam(q2)       = gam_acc(q2)
                        else
                            oh_cme(:,:,:,q2) = real(1.d0-g_onl)*oh_cme(:,:,:,q2) + real(g_onl)*Yeven(q2)%cmat_exp
                            oh_cmo(:,:,:,q2) = real(1.d0-g_onl)*oh_cmo(:,:,:,q2) + real(g_onl)*Yodd(q2)%cmat_exp
                            oh_rhe(:,:,:,q2) = real(1.d0-g_onl)*oh_rhe(:,:,:,q2) + real(g_onl)*Yeven(q2)%rho_exp
                            oh_rho(:,:,:,q2) = real(1.d0-g_onl)*oh_rho(:,:,:,q2) + real(g_onl)*Yodd(q2)%rho_exp
                            oh_gam(q2)       = (1.d0-g_onl)*oh_gam(q2) + g_onl*gam_acc(q2)
                        endif
                        Yeven(q2)%cmat_exp = oh_cme(:,:,:,q2)
                        Yodd(q2)%cmat_exp  = oh_cmo(:,:,:,q2)
                        Yeven(q2)%rho_exp  = oh_rhe(:,:,:,q2)
                        Yodd(q2)%rho_exp   = oh_rho(:,:,:,q2)
                        gam_acc(q2)        = oh_gam(q2)
                    end do
                    l_oh_live = .true.
                    ! snapshot the CURRENT frame as the one the history now lives in
                    if( allocated(prev_real) )then
                        if( size(prev_real) == ncomp )then
                            call prev_real(1)%get_rmat_ptr(onl_pr)
                            if( .not. allocated(ut_hist) ) allocate(ut_hist(size(onl_pr), ncomp))
                            if( size(ut_hist,1) /= size(onl_pr) )then
                                deallocate(ut_hist); allocate(ut_hist(size(onl_pr), ncomp))
                            endif
                            do q2 = 1, ncomp
                                call prev_real(q2)%get_rmat_ptr(onl_pr)
                                ut_hist(:,q2) = real(reshape(onl_pr, [size(onl_pr)]), dp)
                            end do
                        endif
                    endif
                    write(logfhandle,'(A,I0,A,F6.3,A,F9.0)') '>>> FLEX_PCA ONLINE it=', it, &
                        &'  gamma=', g_onl, '  effective N~', real(nval,dp)/max(g_onl,1.d-3)
                    call flush(logfhandle)
                end block
            endif
            ! ---- MARGINAL LIKELIHOOD (the EM objective; NOT resid_energy) ----
            !   -2 log p(y_i) = ||y_i - a_i T_i u_0||^2/sig2  -  h_i' A_i^-1 h_i
            !                   + log det A_i + log det Gamma + const
            ! The first term is FIXED across iterations -- a_i and mu are both held fixed -- so only
            ! the rest is accumulated, and differences between iterations are exact. log det Gamma is
            ! taken from `prior`, which is the Gamma the E-step just USED, not the one it just
            ! produced. Reported per particle so the number is comparable across N.
            ! The filtered M-step (FSC-Wiener + band-limit + mask + re-orthonormalisation) sits
            ! OUTSIDE the likelihood, so this is a DIAGNOSTIC, not a monotonicity guarantee: a rise
            ! means the regularisation is costing more than it buys, which is worth seeing.
            if( l_mix_used )then
                ! the mixture marginal's Omega log det, with the Omega this iteration's E-step
                ! actually USED (the M-step above has already moved mix_Om on)
                nll_tot = nll_tot + real(nval,dp)*ldOm_used
            else
                nll_tot = nll_tot - real(nval,dp)*sum(log(max(prior(1:ncomp), DTINY)))
            endif
            write(logfhandle,'(A,I0,A,ES14.6,A,ES12.4)') '>>> FLEX_PCA PROBE ITER ',it, &
                &'  -2logL/N (varying part)=',nll_tot/real(max(1,nval),dp), &
                &'  delta=',(nll_tot/real(max(1,nval),dp)) - nll_prev
            nll_prev = nll_tot/real(max(1,nval),dp)
            call flush(logfhandle)
            ! COUPLED M-step solve. At every grid point the components share ONE k x k normal matrix
            ! sum_i |CTF|^2 E[z_i z_i'], so the basis volumes have to be solved together. Dividing each
            ! Y_q by the scalar density sum_i |CTF|^2 instead drops the cross-component coupling and
            ! degrades the update to plain subspace iteration. Solve per halfset, then hand a unit
            ! density to compress_exp and skip sampl_dens_correct -- the divide has already happened.
            call solve_coupled_basis_exp(Yeven, rho_e, ncomp)
            call solve_coupled_basis_exp(Yodd,  rho_o, ncomp)
            ! finalize even/odd Y_q, half-set FSC Wiener-merge -> clean band-limited masked basis volume.
            allocate(realvols(ncomp))
            ! the two half-bases are kept, not discarded: comparing them is the only stopping
            ! criterion here that does not embed a dataset-specific constant (see below)
            allocate(eimgs(ncomp), oimgs(ncomp))
            filtsz = max(1, fdim(params%box_crop) - 1)
            allocate(filt(filtsz), corrs(filtsz))
            ! native-lattice inverse KB envelope, applied to each half before FSC/merge/mask
            mstep_gridcorr = prep3D_inv_kbenvelope4mul([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            allocate(fscq_dg(filtsz, ncomp), source=0.)
            do q = 1, ncomp
                ! even half -> band-limited real image (UNmasked, for an unbiased FSC)
                Yeven(q)%rho_exp = 1.0
                call Yeven(q)%compress_exp; call Yeven(q)%ifft
                call realvols(q)%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
                call Yeven(q)%get_rmat_ptr(rmatp); call realvols(q)%set_rmat(rmatp, .false.)
                call realvols(q)%mul(mstep_gridcorr)
                ! odd half
                Yodd(q)%rho_exp = 1.0
                call Yodd(q)%compress_exp; call Yodd(q)%ifft
                call img_o%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
                call Yodd(q)%get_rmat_ptr(rmatp); call img_o%set_rmat(rmatp, .false.)
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
                if( lp_it > 2.0*params%smpd_crop + TINY ) call realvols(q)%bp(0., lp_it)
                call realvols(q)%ifft
                if( params%msk_crop > TINY ) call realvols(q)%mask3D_soft(params%msk_crop, backgr=0.)
                ! ---- FOCUSED HETEROGENEITY (SIMPLE_COV_FOCUS_R1, Angstrom): soft inner
                ! exclusion so the basis explains ONLY the peripheral shell (arms/NBD on
                ! 10049 — 10-24% of motion variance lives beyond r=60 A but never owns an
                ! axis against core+background variance). The mean stays global; only the
                ! heterogeneity subspace is focused — focused-classification semantics.
                if( l_focus_annulus )then
                    block
                        real, pointer :: fpr(:,:,:)
                        integer :: fi, fj, fk
                        real    :: frad, fw
                        call realvols(q)%get_rmat_ptr(fpr)
                        do fk = 1, params%box_crop
                            do fj = 1, params%box_crop
                                do fi = 1, params%box_crop
                                    frad = sqrt(real(fi-params%box_crop/2-1)**2 + &
                                        &real(fj-params%box_crop/2-1)**2 + &
                                        &real(fk-params%box_crop/2-1)**2)*params%smpd_crop
                                    fw = min(1.0, max(0.0, (frad - foc_r1)/max(foc_edge, 1.0)))
                                    fpr(fi,fj,fk) = fpr(fi,fj,fk)*fw
                                end do
                            end do
                        end do
                    end block
                endif
                ! stash both halves, filtered exactly as the merged basis is, so the comparison
                ! below is between two bases that have had identical treatment
                call eimgs(q)%copy(realvols(q))
                call img_o%ifft
                if( lp_it > 2.0*params%smpd_crop + TINY )then
                    call img_o%fft; call img_o%bp(0., lp_it); call img_o%ifft
                endif
                if( params%msk_crop > TINY ) call img_o%mask3D_soft(params%msk_crop, backgr=0.)
                if( l_focus_annulus )then
                    block
                        real, pointer :: fpr(:,:,:)
                        integer :: fi, fj, fk
                        real    :: frad, fw
                        call img_o%get_rmat_ptr(fpr)
                        do fk = 1, params%box_crop
                            do fj = 1, params%box_crop
                                do fi = 1, params%box_crop
                                    frad = sqrt(real(fi-params%box_crop/2-1)**2 + &
                                        &real(fj-params%box_crop/2-1)**2 + &
                                        &real(fk-params%box_crop/2-1)**2)*params%smpd_crop
                                    fw = min(1.0, max(0.0, (frad - foc_r1)/max(foc_edge, 1.0)))
                                    fpr(fi,fj,fk) = fpr(fi,fj,fk)*fw
                                end do
                            end do
                        end do
                    end block
                endif
                call oimgs(q)%copy(img_o)
                call img_o%kill
            end do
            call mstep_gridcorr%kill
            ! ---- AUTO-lp / AUTO-neigs CANDIDATE DIAGNOSTIC (report only) ----
            ! Both knobs are derivable from this one measurement, which the Wiener merge
            ! computes anyway. The HETEROGENEITY resolution -- where the leading components'
            ! split-half FSC crosses 0.143 -- is what lp should encode (NOT the consensus
            ! resolution: compositional signal dies far coarser than the consensus map).
            ! The count of components with significant in-band FSC is what neigs should
            ! encode. Validation bar for any rule built on these: reproduce the measured
            ! ladders (band: lp16 beats lp20 delivered; rank: delivered peaks near 40).
            khi_dg = filtsz
            if( params%lp > 2.0*params%smpd_crop + TINY ) &
                &khi_dg = max(2, min(filtsz, int(dstep_ann/params%lp)))
            nsig_dg = 0
            do q = 1, ncomp
                fmean_dg(q) = sum(fscq_dg(1:khi_dg,q))/real(khi_dg)
                if( fmean_dg(q) > 0.143 ) nsig_dg = nsig_dg + 1
            end do
            ! heterogeneity resolution off the STRONGEST components. Ranked by in-band mean
            ! FSC, NOT by component index: the EM basis is orthonormalised in fixed order and
            ! never variance-sorted (measured on 10076: the 39th/40th components carried the
            ! 2nd/3rd largest latent sd), so "first 4" would read arbitrary components.
            ntop_dg = min(4, ncomp)
            sel_dg  = 0
            do tq_dg = 1, ntop_dg
                fbest_dg = -huge(1.d0); bq_dg = 1
                do q = 1, ncomp
                    if( any(sel_dg(1:tq_dg-1) == q) ) cycle
                    if( fmean_dg(q) > fbest_dg )then
                        fbest_dg = fmean_dg(q); bq_dg = q
                    endif
                end do
                sel_dg(tq_dg) = bq_dg
            end do
            res_dg  = dstep_ann/real(khi_dg)
            do sh = 2, filtsz
                fc = 0.
                do tq_dg = 1, ntop_dg
                    fc = fc + fscq_dg(sh,sel_dg(tq_dg))
                end do
                fc = fc/real(ntop_dg)
                if( fc < 0.143 )then
                    res_dg = dstep_ann/real(sh)
                    exit
                endif
            end do
            write(logfhandle,'(A,I0,A,F7.2,A,I0,A,I0,A)') '>>> FLEX_PCA BAND/RANK it=',it, &
                &'  het-resolution(FSC 0.143)=',res_dg,' A   components with in-band FSC>0.143: ', &
                &nsig_dg,' of ',ncomp,'   (auto-lp / auto-neigs candidates)'
            call flush(logfhandle)
            deallocate(fscq_dg)
            deallocate(filt, corrs)
            ! ---- DATASET-AGNOSTIC STOPPING RULE ----
            ! Every fixed iteration count is a tuned constant standing in for a convergence test,
            ! and the number that is right for one dataset is wrong for the next. The successive-
            ! basis principal angle is not the fix either: it averages over ALL components, so it
            ! is dominated by the tail directions that never reproduce, and it fires at iteration
            ! 3-5 while the leading directions are still moving.
            !
            ! What the algorithm already has, and throws away, is TWO INDEPENDENT HALF-FITS: Yeven
            ! and Yodd are solved separately every iteration and only then merged. The principal
            ! angles between those two spans measure the REPRODUCIBLE dimension of the current
            ! basis directly, and their sum is a soft count of it. Iterating past the point where
            ! that stops rising is refining structure that does not survive a change of particles
            ! -- measured on 10076: the leading eight directions are converged by iteration ~20
            ! while the likelihood keeps falling to iteration 80, and only ~8 reproduce.
            !
            ! So: stop when the reproducible dimension has not improved for COV_EO_PATIENCE
            ! iterations. No particle count, no band, no dataset in the rule.
            ! Comparing the half-bases DIRECTLY does not work, and the reason is worth recording:
            ! both halves are updates to the SAME current basis, so both inherit it and the
            ! agreement saturates immediately -- measured 15.3 of 16 at iteration 1, which is the
            ! shared basis showing through, not reproducible signal.
            !
            ! The quantity that is actually informative is whether the two halves agree about what
            ! to ADD. Projecting each half-basis onto the orthogonal complement of the previous
            ! basis leaves exactly the new directions each half wants, and the principal angles
            ! between THOSE say whether the current update carries signal or noise. When the
            ! updates stop agreeing, further iterations are fitting the half they happen to see.
            if( allocated(prev_real) )then
                call deflate_against_basis(eimgs, ncomp, prev_real, size(prev_real))
                call deflate_against_basis(oimgs, ncomp, prev_real, size(prev_real))
                call cross_half_subspace_angles(eimgs, oimgs, ncomp, sv_eo)
                eo_dim = sum(sv_eo)
            else
                allocate(sv_eo(ncomp), source=0.d0)
                eo_dim = -1.d0            ! first iteration: no previous basis to deflate against
            endif
            write(logfhandle,'(A,I0,A,F8.3,A,I0,A,I0)') '>>> FLEX_PCA PROBE ITER ',it, &
                &'  update agreement (even|odd, prev basis deflated)=',eo_dim, &
                &'   n>=0.9: ',count(sv_eo >= 0.9d0),' of ',ncomp
            call flush(logfhandle)
            if( eo_dim < 0.d0 )then
                continue                  ! nothing to judge yet
            else if( eo_dim > eo_best + COV_EO_TOL )then
                eo_best = eo_dim; eo_stall = 0
            else
                eo_stall = eo_stall + 1
            endif
            deallocate(sv_eo)
            do q = 1, ncomp
                call eimgs(q)%kill; call oimgs(q)%kill
            end do
            deallocate(eimgs, oimgs)
            ! ---- MEAN-SHAPED (CONTRAST) DEFLATION, SIMPLE_COV_EM_DEFLATE=1 ----
            ! a_i is a scalar ML fit of the contrast against the mean, so any FREQUENCY-dependent
            ! part of the per-image scale (envelope/B-factor spread) survives it and lands in the
            ! residual as a term proportional to the consensus shape, coherent across every particle.
            ! It is therefore the strongest single direction available and the EM spends a whole
            ! component on it: measured on 10076, pc001 is a featureless blob in the EM arm AND in
            ! the moment arm. The column path already deflates it in volume space (SIMPLE_COV_RANK1);
            ! this is the same operation applied to the refined basis each iteration.
            ! Only the GLOBAL mean-shaped component is removed -- a domain-local mass change keeps
            ! whatever part of it is orthogonal to the consensus, which is most of it.
            ! SIMPLE_COV_EM_DEFLATE=n with n>1 deflates an n-dimensional ENVELOPE subspace instead
            ! of the single mean direction, which is the gap the paragraph above names. Splitting
            ! the fitted band into n equal shells of the consensus and removing all of them takes
            ! out any per-particle scale that varies smoothly with frequency -- a B-factor spread,
            ! an envelope difference, an ice-thickness gradient -- not just the one scalar a_i
            ! already absorbs. Disjoint Fourier bands are orthogonal by Parseval, but soft-masking
            ! in real space breaks that, so the shells are explicitly orthonormalised rather than
            ! assumed independent.
            if( l_deflate_mean )then
                ndfl = max(1, vdfl)
                ! SIMPLE_COV_DEFLATE_BG=1 appends a FLAT-IN-MASK background template to the
                ! deflation set: the ice/background mode is NOT consensus-shaped (measured
                ! overlap 0.034 on 10049) so consensus shells never remove it, and it owns
                ! the leading axis (27x) on real data whenever the mask is wide enough to
                ! admit solvent. Deflating an explicit background direction lets the mask
                ! OPEN (arms/NBD variance survives) without the nuisance eating the basis.
                if( cov_env_int_on('SIMPLE_COV_DEFLATE_BG') )then
                    ndfl = ndfl + 1
                    l_dfl_bg = .true.
                else
                    l_dfl_bg = .false.
                endif
                ! SIMPLE_COV_DEFLATE_POSE=1: six global pose-derivative templates (below)
                if( cov_env_int_on('SIMPLE_COV_DEFLATE_POSE') )then
                    ndfl = ndfl + 6
                    l_dfl_pose = .true.
                else
                    l_dfl_pose = .false.
                endif
                ndfl_sh = ndfl
                if( l_dfl_bg )   ndfl_sh = ndfl_sh - 1
                if( l_dfl_pose ) ndfl_sh = ndfl_sh - 6
                allocate(dfl_basis(ndfl))
                call mvol_dfl%read_and_crop(params%vols(1), params%smpd, params%box_crop, params%smpd_crop)
                kfr_dfl = covariance_kfromto(params)
                if( l_dfl_bg )then
                    ! last slot: uniform density under the soft spherical mask
                    call dfl_basis(ndfl)%copy(mvol_dfl)
                    call dfl_basis(ndfl)%get_rmat_ptr(rm_dfl)
                    rm_dfl = 0.
                    rm_dfl(1:params%box_crop, 1:params%box_crop, 1:params%box_crop) = 1.
                    if( params%msk_crop > TINY ) call dfl_basis(ndfl)%mask3D_soft(params%msk_crop, backgr=0.)
                    write(logfhandle,'(A)') '>>> FLEX_PCA DEFLATE_BG: flat-in-mask background template &
                        &appended to the deflation set'
                    call flush(logfhandle)
                endif
                do idfl = 1, ndfl_sh
                    call dfl_basis(idfl)%copy(mvol_dfl)
                    if( ndfl_sh > 1 )then
                        ! shell idfl of ndfl, in resolution not shell index, so the bands carry
                        ! comparable signal rather than the outermost one carrying nearly all of it
                        ! dstep_crop / shell, matching covariance_kfromto's own convention.
                        ! The FIRST shell is a pure low-pass: a high-pass edge at k~1 would shave
                        ! the shells where nearly all the consensus power lives.
                        if( idfl == 1 )then
                            res_lo = 0.
                        else
                            res_lo = dstep_ann / &
                                &max(1., real(kfr_dfl(1)) + real(idfl-1)*real(max(1,kfr_dfl(2)-kfr_dfl(1)))/real(ndfl_sh))
                        endif
                        res_hi = dstep_ann / &
                            &max(1., real(kfr_dfl(1)) + real(idfl)  *real(max(1,kfr_dfl(2)-kfr_dfl(1)))/real(ndfl_sh))
                        call dfl_basis(idfl)%fft
                        ! width=1: bp's DEFAULT cosine edge is 10 shells wide PER SIDE, and at
                        ! ndfl=2 each shell of the fitted band is only ~8 shells, so the default
                        ! obliterated the shells and unit-normalisation amplified the residue into
                        ! directions near-orthogonal to everything: rank-2 deflation removed 0.12%
                        ! of basis energy where rank-1 removed 34.6% (measured, hd_defl2)
                        call dfl_basis(idfl)%bp(res_lo, res_hi, width=1.0)
                        call dfl_basis(idfl)%ifft
                        ! real-space padding columns are not guaranteed zero after an ifft, and
                        ! every inner product below runs over the whole padded array
                        call dfl_basis(idfl)%get_rmat_ptr(rm_dfl)
                        rm_dfl(params%box_crop+1:,:,:) = 0.
                        ! band-passing spreads density outside the particle, and the volumes being
                        ! deflated are soft-masked, so the shells must be too or the projection is
                        ! taken against solvent. Deliberately NOT applied at ndfl=1: that is the
                        ! arm measured at ARI 0.176 and it stays bit-identical.
                        if( params%msk_crop > TINY ) call dfl_basis(idfl)%mask3D_soft(params%msk_crop, backgr=0.)
                    endif
                end do
                ! SIMPLE_COV_DEFLATE_POSE=1 appends the six GLOBAL pose-derivative
                ! templates: translation gradients d(mu)/dx,dy,dz and whole-map rotation
                ! derivatives -(a x r).grad(mu) about x/y/z. Per-particle alignment error
                ! IS a small global rotation/shift, so its variance is nuisance in every
                ! dataset regardless of heterogeneity type; DIFFERENTIAL (subunit-local)
                ! motion is orthogonal to whole-map rotation and survives the deflation.
                ! Measured need: on EMPIAR-10028 the delivered "ratcheting" states were
                ! globally rotated consensus copies (pose jitter), and on 3.42 A cryoSPARC
                ! poses that mode collapses -- so deflate it at source. Slots ndfl_sh+1..+6;
                ! Gram-Schmidt below orthonormalises them against the consensus shells.
                if( l_dfl_pose )then
                    call mvol_dfl%get_rmat_ptr(rv_dfl)
                    cp_dfl = real(params%box_crop/2 + 1, dp)
                    do ipdfl = 1, 6
                        call dfl_basis(ndfl_sh+ipdfl)%copy(mvol_dfl)
                        call dfl_basis(ndfl_sh+ipdfl)%get_rmat_ptr(rm_dfl)
                        rm_dfl = 0.
                        do izp = 2, params%box_crop-1
                            do iyp = 2, params%box_crop-1
                                do ixp = 2, params%box_crop-1
                                    gpx = 0.5d0*(real(rv_dfl(ixp+1,iyp,izp),dp) - real(rv_dfl(ixp-1,iyp,izp),dp))
                                    gpy = 0.5d0*(real(rv_dfl(ixp,iyp+1,izp),dp) - real(rv_dfl(ixp,iyp-1,izp),dp))
                                    gpz = 0.5d0*(real(rv_dfl(ixp,iyp,izp+1),dp) - real(rv_dfl(ixp,iyp,izp-1),dp))
                                    select case(ipdfl)
                                        case(1); rm_dfl(ixp,iyp,izp) = real(gpx)
                                        case(2); rm_dfl(ixp,iyp,izp) = real(gpy)
                                        case(3); rm_dfl(ixp,iyp,izp) = real(gpz)
                                        case(4); rm_dfl(ixp,iyp,izp) = &
                                            &real((real(izp,dp)-cp_dfl)*gpy - (real(iyp,dp)-cp_dfl)*gpz)
                                        case(5); rm_dfl(ixp,iyp,izp) = &
                                            &real((real(ixp,dp)-cp_dfl)*gpz - (real(izp,dp)-cp_dfl)*gpx)
                                        case(6); rm_dfl(ixp,iyp,izp) = &
                                            &real((real(iyp,dp)-cp_dfl)*gpx - (real(ixp,dp)-cp_dfl)*gpy)
                                    end select
                                end do
                            end do
                        end do
                        if( params%msk_crop > TINY ) call dfl_basis(ndfl_sh+ipdfl)%mask3D_soft(params%msk_crop, backgr=0.)
                    end do
                    write(logfhandle,'(A)') '>>> FLEX_PCA DEFLATE_POSE: 6 global pose-derivative templates &
                        &appended to the deflation set'
                    call flush(logfhandle)
                endif
                ! shell health: a shell whose norm is tiny relative to the consensus was emptied
                ! by the band edges or the mask, and unit-normalising it manufactures a spurious
                ! deflation direction out of numerical residue -- which is exactly how rank-2
                ! went inert. Print the norms, and keep only shells above a RELATIVE floor.
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
                    do q = 1, ncomp
                        call realvols(q)%get_rmat_ptr(rv_dfl)
                        tot_dfl = tot_dfl + sum(real(rv_dfl,dp)**2)
                        do idfl = 1, nkeep_dfl
                            call dfl_basis(idfl)%get_rmat_ptr(rm_dfl)
                            mv_dfl  = sum(real(rm_dfl,dp)*real(rv_dfl,dp))   ! unit-norm already
                            rem_dfl = rem_dfl + mv_dfl*mv_dfl
                            rv_dfl  = rv_dfl - real(mv_dfl)*rm_dfl
                        end do
                    end do
                    write(logfhandle,'(A,I0,A,I0,A,F6.2,A)') '>>> FLEX_PCA PROBE ITER ',it, &
                        &'  mean-shaped deflation (rank ',nkeep_dfl,') removed ', &
                        &100.d0*rem_dfl/max(tot_dfl,DTINY),' % of basis energy'
                    call flush(logfhandle)
                endif
                do idfl = 1, ndfl
                    call dfl_basis(idfl)%kill
                end do
                deallocate(dfl_basis)
                call mvol_dfl%kill
            endif
            ! orthonormalize the probe volumes -> refined basis
            call orthonormalize_representatives(params, build, realvols, ncomp, utilde, utilde_real, d_new)
            ! ---- CONVERGENCE: principal angles between successive bases ----
            ! This is what n_probe_iters should have been. The measured ladder on Ribosembly saturates
            ! at 2 (25-NN 0.7712 -> 0.7715 from 2 to 3, ceiling pinned at 14/16) while the max-variance
            ! that gets logged keeps oscillating 1.96e4 / 6.71e3 / 1.40e4 -- so the basis had settled
            ! while the quantity being watched had not, and the iteration count was a tuned constant
            ! standing in for a convergence test. Both bases come out of orthonormalize_representatives
            ! ORTHONORMAL, so the singular values of their cross-Gram really are principal-angle
            ! cosines here; align_basis_to_reference only per-vector normalises, which is a no-op on an
            ! already-orthonormal set (it is NOT safe on arbitrary bases -- see its own caveat).
            l_converged = .false.
            if( allocated(prev_real) )then
                call align_basis_to_reference(prev_real, size(prev_real), utilde_real, d_new, Mconv, sconv)
                cos_mean = sum(sconv) / real(max(1,size(sconv)),dp)
                write(logfhandle,'(A,I0,A,F9.6)') '>>> FLEX_PCA PROBE ITER ',it, &
                    &'  mean principal-angle cosine vs previous basis=',cos_mean
                call flush(logfhandle)
                ! The successive-basis cosine is kept as a REPORTED diagnostic only. It averages
                ! over all components, so the never-reproducing tail dominates it and it fires far
                ! too early -- on 10076 at iteration 3-5, with the leading directions still moving
                ! and the likelihood an order of magnitude from its plateau. It becomes the
                ! criterion again only if the even/odd signal is unavailable.
                if( ncomp < 2 .and. cos_mean >= conv_thresh ) l_converged = .true.
                deallocate(Mconv, sconv)
                do q = 1, size(prev_real)
                    call prev_real(q)%kill
                end do
                deallocate(prev_real)
            endif
            allocate(prev_real(d_new))
            do q = 1, d_new
                call prev_real(q)%copy(utilde_real(q))
            end do
            ! replace basis_recs with the refined (projection-ready) basis; eigvals = latent variances
            do q = 1, size(basis_recs)
                call basis_recs(q)%dealloc_rho; call basis_recs(q)%kill
            end do
            deallocate(basis_recs); allocate(basis_recs(d_new))
            if( allocated(eigvals) ) deallocate(eigvals); allocate(eigvals(d_new))
            do q = 1, d_new
                ! projection-ready basis reconstructor from the clean real basis image (mean_rec idiom)
                call init_column_reconstructor(params, build, basis_recs(q))
                call basis_recs(q)%set_rmat(utilde_real(q)%get_rmat(), .false.)   ! NOT add(): see above
                call basis_recs(q)%fft
                call basis_recs(q)%expand_exp
                ! EM Gamma update: the POSTERIOR second moment (1/n) sum_i (z_iq^2 + [A_i^-1]_qq), not the
                ! sample variance of the MAP point estimates. MAP latents are shrunk toward zero by the
                ! prior, so the point-estimate variance underestimates Gamma, which sets a tighter prior
                ! next iteration, which shrinks them further -- a self-reinforcing collapse.
                eigvals(q) = max(gam_acc(min(q,ncomp)), DTINY)
                ! overwrite the eigenvolume MRC with the refined basis vector
                fname = 'flex_pca_pc'//int2str_pad(q,3)//MRC_EXT
                call utilde_real(q)%write(fname, del_if_exists=.true.); call fname%kill
            end do
            ! posterior_frac is the share of Gamma that the old point-estimate update was throwing away;
            ! mean(z) should sit near zero -- drift there means the consensus mean mu is off.
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PROBE EM Gamma update (n=',nval,' valid particles):'
            do q = 1, min(ncomp,10)
                mu_q = sum(z(:,q)) / real(max(1,npp),dp)
                sd_q = sum(z(:,q)**2) / real(max(1,nval),dp)
                write(logfhandle,'(A,I3,A,ES11.3,A,ES11.3,A,F7.3,A,ES10.2)') '>>>   z',q, &
                    &'  <z^2>=',sd_q,'  Gamma=',gam_acc(q), &
                    &'  posterior_frac=',real(1.d0 - sd_q/max(gam_acc(q),DTINY)), &
                    &'  mean(z)=',mu_q
            end do
            call flush(logfhandle)
            ncomp = d_new
            write(logfhandle,'(A,I0,A,ES12.4,A,ES12.4,A,F8.1)') '>>> FLEX_PCA PROBE ITER ',it, &
                &' refined dim=',real(ncomp),' max var=',maxval(eigvals),' seconds=',toc(t_it)
            call flush(logfhandle)
            ! cleanup per-iteration scratch
            do i = 1, size(orientations); call orientations(i)%kill; end do
            do ithr = 1, nthr
                call cleanup_plane(mean_fpl(ithr))
                do q = 1, size(basis_fpls,1); call cleanup_plane(basis_fpls(q,ithr)); end do
                if( l_bank_fwd )then
                    call cleanup_plane(mean_fplC(ithr))
                    do q = 1, size(basis_fplsC,1); call cleanup_plane(basis_fplsC(q,ithr)); end do
                endif
            end do
            if( l_bank_fwd ) deallocate(mean_fplC, basis_fplsC)
            if( allocated(stats_es) ) deallocate(stats_es)
            call cleanup_rec_buffers(build, fpls)
            do q = 1, size(Yeven); call Yeven(q)%dealloc_rho; call Yeven(q)%kill; call Yodd(q)%dealloc_rho; call Yodd(q)%kill; end do
            do q = 1, size(utilde); call utilde(q)%dealloc_rho; call utilde(q)%kill; call utilde_real(q)%kill; end do
            do q = 1, size(realvols); call realvols(q)%kill; end do
            deallocate(Yeven, Yodd, utilde, utilde_real, realvols, rho_e, rho_o, prior)
            deallocate(Gth, Ath, bth, cth, zth, basis_fpls, mean_fpl, orientations, zbatch, dens, valid, valid_e, valid_o, eo)
            deallocate(Ainvth, Acpth, gam_thr, gam_acc, nval_thr, hth, nll_thr)
            ! The MCFA state (mix_*) and its work arrays (rhs0th etc) are deliberately NOT
            ! deallocated here: this block runs at the end of EVERY iteration, and the mixture
            ! must survive from one iteration's M-step to the next iteration's E-step. Putting
            ! their cleanup here was the crash: iteration 3 initialised the mixture, this block
            ! freed it, l_mix_active stayed true, and iteration 4's E-step walked into an
            ! unallocated mix_Ominv -- a null descriptor, which is why -fcheck=bounds stayed
            ! silent while every thread segfaulted on the same line. They are subroutine-local
            ! allocatables, so Fortran frees them at exit; the resize block at iteration start
            ! handles dimension changes.
            if( allocated(gam_dbg) ) deallocate(gam_dbg)
            if( allocated(pGc) ) deallocate(pGc, pbc, pcc, pzh, pcon, pre, prme)
            ! NOT used to stop. The update agreement decays monotonically from the first
            ! iteration (measured 15.0 -> 13.7 -> 10.6 -> 8.2, strongly-reproducible update
            ! directions 15 -> 2 -> 0), so any "peaked and stalled" rule fires at iteration 4-5.
            ! But the basis demonstrably keeps improving to ~20 on the same data, because an
            ! update whose DIRECTION is not reproducible can still carry a reproducible bias that
            ! accumulates. So this is a genuine diagnostic of update quality and NOT a stopping
            ! rule, and pretending otherwise would just be a new tuned constant wearing a
            ! principled hat. The defensible in-run rule is cross-fit agreement between two
            ! INDEPENDENT half-fits, which costs 2x; see the acceleration note.
            if( .false. .and. eo_stall >= eo_patience ) l_converged = .true.
            if( l_converged )then
                write(logfhandle,'(A,I0,A,F8.3,A,I0,A)') '>>> FLEX_PCA PROBE converged after ',it, &
                    &' iterations: reproducible dim (even|odd) peaked at ',eo_best, &
                    &' and did not improve for ',eo_patience,' iterations'
                call flush(logfhandle)
                exit
            endif
        end do
        if( allocated(prev_real) )then
            do q = 1, size(prev_real)
                call prev_real(q)%kill
            end do
            deallocate(prev_real)
        endif
        if( l_pol_es )then
            call polar_grid_kill(pg_es)
            if( allocated(UsallE) ) deallocate(UsallE, CfE, Cm0E, c00E, UbankE, CspE, &
                &xws_es, wr_es, wrd_es, Reb_es)
            if( allocated(dir_es) ) deallocate(dir_es, cae, sae, dused_es, rmatb_es, nrmb_es)
            if( allocated(hex_es) ) deallocate(hex_es, kex_es)
            if( l_pol_gpu_any )then
                call flex_gpu_poles_free_f
                ! resident hybrid volumes (idempotent; also freed by the banked-adjoint teardown)
                call flex_gpu_estep_free_f
            endif
            if( allocated(stats_pg) ) deallocate(stats_pg)
            if( allocated(vld_pg)   ) deallocate(vld_pg, dir_pg, ca_pg, sa_pg)
        endif
        if( l_bank_adj )then
            call flex_gpu_coupled_bank_free_f
            call flex_gpu_estep_free_f
            if( l_gpu_prep )then
                call flex_gpu_prep_free_f
                deallocate(ctfp_arr, shf_pb)
                if( allocated(sig2_pb) ) deallocate(sig2_pb)
            endif
            call dirs_pb%kill
            deallocate(rmat_pb, nrm_pb, dirof_pb, cap, sap, o_pb_thr, cap1, sap0)
            deallocate(seg_dir_b, seg_beg_b, seg_cnt_b)
            if( allocated(tmpc) ) deallocate(tmpc)
            if( allocated(avec_es) ) deallocate(avec_es, vld_es)
        endif
        deallocate(z, ppinds)

    contains

        subroutine ensure_bank_ctf_planes( ref, mfc, bfc, nc )
            type(fplane_type), intent(in)    :: ref
            type(fplane_type), intent(inout) :: mfc
            integer,           intent(in)    :: nc
            type(fplane_type), intent(inout) :: bfc(nc)
            integer :: qq
            if( .not. allocated(mfc%cmplx_plane) )then
                allocate(mfc%cmplx_plane, mold=ref%cmplx_plane)
                mfc%cmplx_plane = CMPLX_ZERO
            endif
            mfc%frlims = ref%frlims; mfc%nyq = ref%nyq
            do qq = 1, nc
                if( .not. allocated(bfc(qq)%cmplx_plane) )then
                    allocate(bfc(qq)%cmplx_plane, mold=ref%cmplx_plane)
                    bfc(qq)%cmplx_plane = CMPLX_ZERO
                endif
                bfc(qq)%frlims = ref%frlims; bfc(qq)%nyq = ref%nyq
            end do
        end subroutine ensure_bank_ctf_planes

    end subroutine probe_subspace_iteration

    !> log(det A) for symmetric positive-definite A, by Cholesky on a private copy.
    !!
    !! This is the term the probe's `resid_energy` has always been missing. resid_energy is the JOINT
    !! MAP objective ||.||^2/sig2 + z'Gamma^-1 z evaluated at zhat, which decreases by construction and
    !! therefore says nothing about convergence. The MARGINAL likelihood needs log det A_i and
    !! log det Gamma as well, and without them the EM has no objective to watch -- which is why the
    !! iteration count was a tuned constant standing in for a stopping rule.
    pure subroutine spd_logdet_dp( A, n, logdet, ok )
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: A(n,n)
        real(dp), intent(out) :: logdet
        logical,  intent(out) :: ok
        real(dp) :: L(n,n), s
        integer  :: i, j
        L      = 0.d0
        logdet = 0.d0
        ok     = .false.
        do j = 1, n
            s = A(j,j) - sum(L(j,1:j-1)**2)
            if( s <= 0.d0 ) return
            L(j,j) = sqrt(s)
            logdet = logdet + 2.d0*log(L(j,j))
            do i = j+1, n
                L(i,j) = (A(i,j) - sum(L(i,1:j-1)*L(j,1:j-1))) / L(j,j)
            end do
        end do
        ok = .true.
    end subroutine spd_logdet_dp

    !> Project a stack of volumes onto the orthogonal complement of an ORTHONORMAL basis stack.
    !! Used to isolate what an EM update ADDS to the current basis, so two half-updates can be
    !! compared without the basis they share dominating the comparison.
    subroutine deflate_against_basis( imgs, n, basis, nb )
        integer,     intent(in)    :: n, nb
        type(image), intent(inout) :: imgs(n), basis(nb)
        real, pointer :: rv(:,:,:), rb(:,:,:)
        real(dp) :: c, bb
        integer  :: i, k
        do k = 1, nb
            call basis(k)%get_rmat_ptr(rb)
            bb = sum(real(rb,dp)*real(rb,dp))
            if( bb <= DTINY ) cycle
            do i = 1, n
                call imgs(i)%get_rmat_ptr(rv)
                c  = sum(real(rv,dp)*real(rb,dp)) / bb
                rv = rv - real(c)*rb
            end do
        end do
    end subroutine deflate_against_basis

    !> Principal angles between the subspaces spanned by two NON-orthonormal volume stacks.
    !!
    !! align_basis_to_reference only per-vector normalises, which is exact only when both inputs
    !! are already orthonormal -- it says so itself. The even/odd bases coming out of the coupled
    !! solve are not, so this computes the honest quantity: for spans E and O the principal-angle
    !! cosines are the singular values of (E'E)^-1/2 (E'O) (O'O)^-1/2. Everything is done in the
    !! n x n Gram algebra, so the volume work is three symmetric Gram products and nothing else.
    subroutine cross_half_subspace_angles( eimgs, oimgs, n, svals )
        integer,     intent(in)    :: n
        type(image), intent(inout) :: eimgs(n), oimgs(n)
        real(dp), allocatable, intent(out) :: svals(:)
        real, pointer :: ri(:,:,:), rj(:,:,:)
        real(dp), allocatable :: Gee(:,:), Goo(:,:), Geo(:,:), We(:,:), Wo(:,:)
        real(dp), allocatable :: ev(:), vec(:,:), M(:,:), MtM(:,:), V2(:,:), ev2(:)
        integer  :: i, j, nrot
        real(dp) :: lam
        allocate(Gee(n,n), Goo(n,n), Geo(n,n), source=0.d0)
        do i = 1, n
            call eimgs(i)%get_rmat_ptr(ri)
            do j = i, n
                call eimgs(j)%get_rmat_ptr(rj)
                Gee(i,j) = sum(real(ri,dp)*real(rj,dp)); Gee(j,i) = Gee(i,j)
            end do
            do j = 1, n
                call oimgs(j)%get_rmat_ptr(rj)
                Geo(i,j) = sum(real(ri,dp)*real(rj,dp))
            end do
        end do
        do i = 1, n
            call oimgs(i)%get_rmat_ptr(ri)
            do j = i, n
                call oimgs(j)%get_rmat_ptr(rj)
                Goo(i,j) = sum(real(ri,dp)*real(rj,dp)); Goo(j,i) = Goo(i,j)
            end do
        end do
        ! inverse square roots by symmetric eigendecomposition, with a relative floor so a
        ! collapsed half-basis direction cannot blow the whitening up
        allocate(We(n,n), Wo(n,n), source=0.d0)
        call inv_sqrt_sym(Gee, n, We)
        call inv_sqrt_sym(Goo, n, Wo)
        allocate(M(n,n), MtM(n,n), V2(n,n), ev2(n), svals(n))
        M   = matmul(We, matmul(Geo, Wo))
        MtM = matmul(transpose(M), M)
        call jacobi(MtM, n, n, ev2, V2, nrot)
        call eigsrt(ev2, V2, n, n)
        do i = 1, n
            svals(i) = sqrt(max(0.d0, min(1.d0, ev2(i))))
        end do
        deallocate(Gee, Goo, Geo, We, Wo, M, MtM, V2, ev2)
      contains
        subroutine inv_sqrt_sym( A, m, Ainvsq )
            integer,  intent(in)  :: m
            real(dp), intent(in)  :: A(m,m)
            real(dp), intent(out) :: Ainvsq(m,m)
            real(dp), allocatable :: Aw(:,:), Vv(:,:), ee(:)
            integer :: k, nr
            real(dp) :: mx
            allocate(Aw(m,m), Vv(m,m), ee(m))
            Aw = A
            call jacobi(Aw, m, m, ee, Vv, nr)
            mx = maxval(ee)
            Ainvsq = 0.d0
            do k = 1, m
                lam = ee(k)
                if( lam > 1.d-10*max(mx,DTINY) )then
                    Ainvsq = Ainvsq + matmul(reshape(Vv(:,k),[m,1]), &
                        &reshape(Vv(:,k),[1,m])) / sqrt(lam)
                endif
            end do
            deallocate(Aw, Vv, ee)
        end subroutine inv_sqrt_sym
    end subroutine cross_half_subspace_angles

    !>  z' M z for symmetric M.
    pure function quad_form( M, z, n ) result( val )
        integer,  intent(in) :: n
        real(dp), intent(in) :: M(n,n), z(n)
        real(dp) :: val
        integer  :: q, r
        val = 0.d0
        do r = 1, n
            do q = 1, n
                val = val + z(q)*M(q,r)*z(r)
            end do
        end do
    end function quad_form

    !> In-place symmetric positive-definite solve A x = b (b overwritten by x) via Cholesky. A is first
    !! scaled by its mean diagonal, so the retry ridge is RELATIVE: an absolute ridge either swamps a
    !! small-diagonal system or fails to rescue a large one, and the b=0 fallback then collapses the
    !! latents of essentially every particle.
    subroutine spd_solve_dp( A, b, n )
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: A(n,n), b(n)
        real(dp) :: L(n,n), s, y(n), ridge, dscale
        integer  :: i, j, attempt
        dscale = 0.d0
        do i = 1, n
            dscale = dscale + abs(A(i,i))
        end do
        dscale = dscale / real(n,dp)
        if( dscale > 0.d0 )then
            A = A / dscale
            b = b / dscale
        endif
        do attempt = 1, 3
            L = 0.d0
            do j = 1, n
                s = A(j,j) - sum(L(j,1:j-1)**2)
                if( s <= 0.d0 ) exit
                L(j,j) = sqrt(s)
                do i = j+1, n
                    L(i,j) = (A(i,j) - sum(L(i,1:j-1)*L(j,1:j-1))) / L(j,j)
                end do
            end do
            if( j > n )then
                ! forward/back substitution
                do i = 1, n
                    y(i) = (b(i) - sum(L(i,1:i-1)*y(1:i-1))) / L(i,i)
                end do
                do i = n, 1, -1
                    b(i) = (y(i) - sum(L(i+1:n,i)*b(i+1:n))) / L(i,i)
                end do
                return
            endif
            ridge = 1.d-8 * (abs(A(1,1))+1.d0) * (10.d0**(attempt-1))
            do i = 1, n
                A(i,i) = A(i,i) + ridge
            end do
        end do
        b = 0.d0
    end subroutine spd_solve_dp

    !>  Inverse of a symmetric positive-definite matrix by Cholesky, using the same diagonal rescaling
    !!  and ridge-escalation policy as spd_solve_dp so a matrix that solves also inverts. A is DESTROYED.
    !!  A rank-deficient input is rescued by the ridge exactly as in spd_solve_dp, so the result can be
    !!  large; only a matrix that all three attempts fail on returns zeros. The embedding never hits that
    !!  path: A = (a^2/sig2) G + diag(prior) with G PSD, so lambda_min(A) >= min(prior) and every
    !!  [A^-1]_qq is bounded above by max(eigvals) -- a dead particle reproduces the prior, as it should.
    !> One particle's MAP solve at fixed contrast, with optional ECM alternations of the
    !! closed-form scale against the CURRENT basis:
    !!     a <- (m'y + b'z) / (||m||^2 + 2 c'z + z'Gz + tr(G A^-1))
    !! nml = 0 reproduces the historical fixed-a solve exactly, call order and all: the log det
    !! comes off the pristine A, the inverse off a copy, the solve off A (the latter two rescale
    !! their argument in place). The tr(G A^-1) term is the posterior variance of z and is
    !! load-bearing -- dropping it uses E[z]E[z]' for E[zz'] and biases a high (measured in the
    !! embedding-stage ECM). The bracket matches the projection fit's [0.1, 5].
    subroutine probe_solve_ecm( n, G, b, c, myv, e_mm, prior_, sig2, nml, a, z_, Ainv_, ldA, lok, quad )
        integer,  intent(in)    :: n, nml
        real(dp), intent(in)    :: G(n,n), b(n), c(n), myv, e_mm, prior_(n), sig2
        real(dp), intent(inout) :: a
        real(dp), intent(out)   :: z_(n), Ainv_(n,n), ldA
        logical,  intent(out)   :: lok
        real(dp), intent(out)   :: quad
        real(dp) :: Amat(n,n), Acp(n,n), h(n), aa, a_num, a_den
        integer  :: q, icm
        do icm = 0, max(0, nml)
            aa   = a*a
            Amat = (aa/sig2)*G
            do q = 1, n
                Amat(q,q) = Amat(q,q) + prior_(q)
                z_(q)     = (a*b(q) - aa*c(q))/sig2
            end do
            Acp = Amat
            call spd_logdet_dp(Amat, n, ldA, lok)
            h = z_
            call spd_inv_dp(Acp, Ainv_, n)
            call spd_solve_dp(Amat, z_, n)
            quad = dot_product(h, z_)
            if( icm >= max(0, nml) ) exit
            a_num = myv + dot_product(b, z_)
            a_den = e_mm + 2.d0*dot_product(c, z_) + quad_form(G, z_, n) + sum(G*Ainv_)
            if( a_den > DTINY ) a = min(5.0d0, max(0.1d0, a_num/a_den))
        end do
    end subroutine probe_solve_ecm

    !> Per-particle MCFA posterior from fetched sufficient statistics (G, b, c) -- the ONE
    !! source for the mixture E-step, shared by the plain CPU body and the fused device
    !! body (the device computes the same G/b/c on card; only this host solve differs from
    !! the single-Gaussian path). Updates the caller's thread-slice accumulators exactly as
    !! the historical in-line block did.
    subroutine probe_solve_mix( n, kmix, G, b, c, myv, e_mm, nml, a, sig2, Ominv, Omxi, lpi, xiOx, &
            &zbar, Ainv_, dens_, ldA, lok, nll_add, sr_acc, sm_acc, smm_acc, sainv_acc )
        integer,  intent(in)    :: n, kmix, nml
        real(dp), intent(in)    :: G(n,n), b(n), c(n), myv, e_mm, sig2
        real(dp), intent(inout) :: a
        real(dp), intent(in)    :: Ominv(n,n), Omxi(n,kmix), lpi(kmix), xiOx(kmix)
        real(dp), intent(out)   :: zbar(n), Ainv_(n,n), dens_(n,n), ldA, nll_add
        logical,  intent(out)   :: lok
        real(dp), intent(inout) :: sr_acc(kmix), sm_acc(n,kmix), smm_acc(n,n,kmix), sainv_acc(n,n)
        real(dp) :: Amat(n,n), Acp(n,n), rhs0(n), mk(n,kmix), lw(kmix), rk(kmix)
        real(dp) :: aa, lwm, wsm, a_num, a_den
        integer  :: q, r, kk, icm
        do icm = 0, max(0, nml)
            aa   = a*a
            Amat = (aa/sig2)*G + Ominv
            Acp  = Amat
            call spd_logdet_dp(Amat, n, ldA, lok)
            call spd_inv_dp(Acp, Ainv_, n)
            do q = 1, n
                rhs0(q) = (a*b(q) - aa*c(q))/sig2
            end do
            do kk = 1, kmix
                mk(:,kk) = matmul(Ainv_, rhs0 + Omxi(:,kk))
                lw(kk)   = lpi(kk) - 0.5d0*xiOx(kk) + 0.5d0*dot_product(rhs0 + Omxi(:,kk), mk(:,kk))
            end do
            lwm = maxval(lw)
            wsm = 0.d0
            do kk = 1, kmix
                rk(kk) = exp(max(-7.d2, lw(kk) - lwm))
                wsm    = wsm + rk(kk)
            end do
            rk = rk / max(wsm, DTINY)
            zbar = 0.d0
            do kk = 1, kmix
                zbar = zbar + rk(kk)*mk(:,kk)
            end do
            ! E[zz'|y] = A^-1 + sum_k r_k m_k m_k' (between-component spread is real variance)
            dens_ = Ainv_
            do kk = 1, kmix
                do r = 1, n
                    do q = 1, n
                        dens_(q,r) = dens_(q,r) + rk(kk)*mk(q,kk)*mk(r,kk)
                    end do
                end do
            end do
            if( icm >= max(0, nml) ) exit
            ! ECM amplitude update under the mixture: the plain path's a-update with the
            ! mixture posterior moments in place of the single-Gaussian ones. This is the
            ! lever that separates on-particle amplitude (ice) from conformation.
            a_num = myv + dot_product(b, zbar)
            a_den = e_mm + 2.d0*dot_product(c, zbar) + sum(G*dens_)
            if( a_den > DTINY ) a = min(5.0d0, max(0.1d0, a_num/a_den))
        end do
        ! mixture marginal, varying part; the N*logdet(Omega) term is added once globally
        nll_add = ldA - 2.d0*(lwm + log(max(wsm, DTINY)))
        sainv_acc = sainv_acc + Ainv_
        do kk = 1, kmix
            sr_acc(kk)   = sr_acc(kk)   + rk(kk)
            sm_acc(:,kk) = sm_acc(:,kk) + rk(kk)*mk(:,kk)
            do r = 1, n
                smm_acc(:,r,kk) = smm_acc(:,r,kk) + rk(kk)*mk(:,kk)*mk(r,kk)
            end do
        end do
    end subroutine probe_solve_mix

    !> MCFA initialisation. K=1 pins the single component at the origin over the current Gamma
    !! diagonal -- with the diagonal constraint in mcfa_condition this makes the mixture path
    !! reproduce the plain PPCA EM exactly (the regression test). K>=2 seeds by deterministic
    !! farthest-point selection on the current latents and polishes with Lloyd iterations;
    !! nothing here is random, so a rerun reproduces bit-for-bit.
    subroutine mcfa_init( z, nptcls, ncomp, kmix, gam_sum, nval, xi, ppi, Om )
        integer,  intent(in)  :: nptcls, ncomp, kmix, nval
        real(dp), intent(in)  :: z(nptcls,ncomp), gam_sum(ncomp)
        real(dp), intent(out) :: xi(ncomp,kmix), ppi(kmix), Om(ncomp,ncomp)
        integer,  allocatable :: rows(:), lab(:), cnt(:)
        real(dp), allocatable :: d2min(:)
        real(dp) :: best, dd
        integer  :: n, i, j, k, kbest, it2
        if( kmix == 1 )then
            xi  = 0.d0
            ppi = 1.d0
            Om  = 0.d0
            do i = 1, ncomp
                Om(i,i) = max(gam_sum(i)/real(max(1,nval),dp), DTINY)
            end do
            return
        endif
        ! particles the E-step skipped (state zero) keep z = 0 and must not seed a component
        allocate(rows(nptcls))
        n = 0
        do i = 1, nptcls
            if( any(z(i,1:ncomp) /= 0.d0) )then
                n = n + 1
                rows(n) = i
            endif
        end do
        allocate(lab(n), cnt(kmix), d2min(n))
        best = -1.d0
        j = 1
        do i = 1, n
            dd = sum(z(rows(i),1:ncomp)**2)
            if( dd > best )then
                best = dd
                j = i
            endif
        end do
        xi(:,1) = z(rows(j),1:ncomp)
        d2min = huge(0.d0)
        do k = 2, kmix
            best  = -1.d0
            kbest = 1
            do i = 1, n
                dd = sum((z(rows(i),1:ncomp) - xi(:,k-1))**2)
                if( dd < d2min(i) ) d2min(i) = dd
                if( d2min(i) > best )then
                    best  = d2min(i)
                    kbest = i
                endif
            end do
            xi(:,k) = z(rows(kbest),1:ncomp)
        end do
        do it2 = 1, 12
            cnt = 0
            do i = 1, n
                best  = huge(0.d0)
                kbest = 1
                do k = 1, kmix
                    dd = sum((z(rows(i),1:ncomp) - xi(:,k))**2)
                    if( dd < best )then
                        best  = dd
                        kbest = k
                    endif
                end do
                lab(i) = kbest
            end do
            xi = 0.d0
            do i = 1, n
                xi(:,lab(i)) = xi(:,lab(i)) + z(rows(i),1:ncomp)
                cnt(lab(i))  = cnt(lab(i)) + 1
            end do
            do k = 1, kmix
                if( cnt(k) > 0 ) xi(:,k) = xi(:,k)/real(cnt(k),dp)
            end do
        end do
        ! mixing proportions floored so no component is born dead
        do k = 1, kmix
            ppi(k) = max(real(cnt(k),dp)/real(max(1,n),dp), 0.25d0/real(kmix,dp))
        end do
        ppi = ppi / sum(ppi)
        ! tied Omega from the pooled within-cluster scatter of the POINT latents; the M-step
        ! adds the posterior covariance back from its own accumulators, so this only has to
        ! carry a sane scale, not be unbiased
        Om = 0.d0
        do i = 1, n
            do j = 1, ncomp
                Om(:,j) = Om(:,j) + (z(rows(i),1:ncomp) - xi(:,lab(i)))*(z(rows(i),j) - xi(j,lab(i)))
            end do
        end do
        Om = Om / real(max(1,n),dp)
        deallocate(rows, lab, cnt, d2min)
    end subroutine mcfa_init

    !> Symmetrise, eigen-floor and invert the tied prior covariance. diag_only enforces the
    !! K=1 reduction-test constraint (matching the plain path's diagonal Gamma exactly).
    subroutine mcfa_condition( ncomp, diag_only, Om, Ominv, ldOm )
        integer,  intent(in)    :: ncomp
        logical,  intent(in)    :: diag_only
        real(dp), intent(inout) :: Om(ncomp,ncomp)
        real(dp), intent(out)   :: Ominv(ncomp,ncomp), ldOm
        real(dp) :: ev(ncomp), evec(ncomp,ncomp), work(ncomp,ncomp), floor_ev
        integer  :: q, r2, nrot
        if( diag_only )then
            do q = 1, ncomp
                do r2 = 1, ncomp
                    if( q /= r2 ) Om(q,r2) = 0.d0
                end do
            end do
            Ominv = 0.d0
            ldOm  = 0.d0
            do q = 1, ncomp
                Om(q,q)    = max(Om(q,q), DTINY)
                Ominv(q,q) = 1.d0/Om(q,q)
                ldOm       = ldOm + log(Om(q,q))
            end do
            return
        endif
        Om   = 0.5d0*(Om + transpose(Om))
        work = Om
        call jacobi(work, ncomp, ncomp, ev, evec, nrot)
        floor_ev = max(1.d-6*maxval(ev), DTINY)
        ldOm = 0.d0
        do q = 1, ncomp
            ev(q) = max(ev(q), floor_ev)
            ldOm  = ldOm + log(ev(q))
        end do
        do q = 1, ncomp
            do r2 = 1, ncomp
                Om(q,r2)    = sum(evec(q,:)*ev(:)*evec(r2,:))
                Ominv(q,r2) = sum(evec(q,:)/ev(:)*evec(r2,:))
            end do
        end do
    end subroutine mcfa_condition

    !> MCFA M-step from REDUCED sufficient statistics (the thread reduction happens at
    !! the call site, so a rotated running-averaged history can be substituted for the
    !! this-iteration statistics transparently):
    !!   pi_k = Sr_k/N,  xi_k = Sm_k/Sr_k,
    !!   Omega = (1/N)[ sum_i A_i^-1 + sum_k (Smm_k - xi Sm' - Sm xi' + Sr xi xi') ]
    !! pin_origin (K=1) keeps xi at 0 and Omega diagonal -- the plain-PPCA reduction.
    subroutine mcfa_mstep( ncomp, kmix, nval, sr, sm, smm, sai, &
        &pin_origin, ppi, xi, Om, Ominv, ldOm )
        integer,  intent(in)    :: ncomp, kmix, nval
        real(dp), intent(in)    :: sr(kmix), sm(ncomp,kmix)
        real(dp), intent(in)    :: smm(ncomp,ncomp,kmix), sai(ncomp,ncomp)
        logical,  intent(in)    :: pin_origin
        real(dp), intent(inout) :: ppi(kmix), xi(ncomp,kmix), Om(ncomp,ncomp)
        real(dp), intent(out)   :: Ominv(ncomp,ncomp), ldOm
        integer  :: k, r2
        do k = 1, kmix
            ppi(k) = max(sr(k)/real(max(1,nval),dp), 1.d-4)
        end do
        ppi = ppi / sum(ppi)
        if( .not. pin_origin )then
            do k = 1, kmix
                ! a starved component keeps its old mean rather than dividing by ~0
                if( sr(k) > 1.d-8*real(max(1,nval),dp) ) xi(:,k) = sm(:,k)/sr(k)
            end do
        endif
        Om = sai
        do k = 1, kmix
            do r2 = 1, ncomp
                Om(:,r2) = Om(:,r2) + smm(:,r2,k) - xi(:,k)*sm(r2,k) - sm(:,k)*xi(r2,k) &
                    &+ sr(k)*xi(:,k)*xi(r2,k)
            end do
        end do
        Om = Om / real(max(1,nval),dp)
        call mcfa_condition(ncomp, pin_origin, Om, Ominv, ldOm)
    end subroutine mcfa_mstep


    subroutine spd_inv_dp( A, Ainv, n )
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: A(n,n)
        real(dp), intent(out)   :: Ainv(n,n)
        real(dp) :: L(n,n), Linv(n,n), s, ridge, dscale
        integer  :: i, j, attempt
        Ainv   = 0.d0
        dscale = 0.d0
        do i = 1, n
            dscale = dscale + abs(A(i,i))
        end do
        dscale = dscale / real(n,dp)
        if( dscale > 0.d0 ) A = A / dscale
        do attempt = 1, 3
            L = 0.d0
            do j = 1, n
                s = A(j,j) - sum(L(j,1:j-1)**2)
                if( s <= 0.d0 ) exit
                L(j,j) = sqrt(s)
                do i = j+1, n
                    L(i,j) = (A(i,j) - sum(L(i,1:j-1)*L(j,1:j-1))) / L(j,j)
                end do
            end do
            if( j > n )then
                ! Linv = L^-1 by forward substitution on the identity, column by column
                Linv = 0.d0
                do j = 1, n
                    Linv(j,j) = 1.d0 / L(j,j)
                    do i = j+1, n
                        Linv(i,j) = -sum(L(i,j:i-1)*Linv(j:i-1,j)) / L(i,i)
                    end do
                end do
                ! A = L L' => A^-1 = (L^-1)' (L^-1), lower-triangular so the sum starts at max(i,j)
                do i = 1, n
                    do j = 1, i
                        Ainv(i,j) = sum(Linv(i:n,i)*Linv(i:n,j))
                        Ainv(j,i) = Ainv(i,j)
                    end do
                end do
                ! undo the rescaling: A_orig = dscale*A_scaled, so A_orig^-1 = A_scaled^-1 / dscale
                if( dscale > 0.d0 ) Ainv = Ainv / dscale
                return
            endif
            ridge = 1.d-8 * (abs(A(1,1))+1.d0) * (10.d0**(attempt-1))
            do i = 1, n
                A(i,i) = A(i,i) + ridge
            end do
        end do
    end subroutine spd_inv_dp

    !> Complex Hermitian-half plane inner product over the native k<=0 half (the planes are stored
    !! half-plane, k in [kmin,0]). Properness is handled by the caller: reduced_covariance_solve
    !! accumulates Re(b_q)Re(b_r), not |b|^2, and pairs it with the 0.5*sig2 noise scaling, since
    !! E[Re(b_q)Re(b_r)]_noise = 0.5*sig2*G. The optional half selects a checkerboard sub-half.
    !> The (h,k) sample sequence cov_herm_inner visits, in exactly its order, for a plane pair that
    !! shares nyq and bounds with `ref` and is called without the `half` selector. Returns nsamp = -1
    !! if the sequence does not fit `slist`, in which case the caller must use cov_herm_inner directly.
    subroutine cov_herm_sample_list( ref, slist, nsamp )
        type(fplane_type), intent(in)  :: ref
        integer,           intent(out) :: slist(:,:)
        integer,           intent(out) :: nsamp
        integer :: h, k, hmin, hmax, kmin, kmax, nyq_eff, pf, h_hi, nyq_disk, k_sq
        pf      = OSMPL_PAD_FAC
        nyq_eff = ref%nyq
        if( nyq_eff <= 0 ) nyq_eff = ubound(ref%cmplx_plane,1)
        hmin = max(pf*ceil_div(lbound(ref%cmplx_plane,1),pf), pf*ceil_div(-nyq_eff,pf))
        hmax = min(pf*floor_div(ubound(ref%cmplx_plane,1),pf), pf*floor_div(nyq_eff,pf))
        kmin = max(pf*ceil_div(lbound(ref%cmplx_plane,2),pf), pf*ceil_div(-nyq_eff,pf))
        kmax = min(0, pf*floor_div(nyq_eff,pf))
        nyq_disk = nyq_eff * (nyq_eff + 1)
        nsamp    = 0
        do k = kmin, kmax, pf
            h_hi = hmax
            if( k == 0 ) h_hi = 0
            k_sq = k*k
            if( k_sq > nyq_disk ) cycle
            do h = hmin, h_hi, pf
                if( h*h + k_sq > nyq_disk ) cycle
                nsamp = nsamp + 1
                if( nsamp > size(slist,2) )then
                    nsamp = -1
                    return
                endif
                slist(1,nsamp) = h
                slist(2,nsamp) = k
            end do
        end do
    end subroutine cov_herm_sample_list

    function cov_herm_inner( lhs, rhs, half ) result( val )
        type(fplane_type), intent(in) :: lhs, rhs
        integer, optional, intent(in) :: half
        complex(dp) :: val, acc
        integer :: h, k, hmin, hmax, kmin, kmax, nyq_eff, pf, h_hi, hlf, par, nyq_disk, k_sq
        hlf = 0
        if( present(half) ) hlf = half
        acc = cmplx(0.d0,0.d0,dp)
        pf  = OSMPL_PAD_FAC
        nyq_eff = lhs%nyq
        if( rhs%nyq > 0 ) nyq_eff = min(nyq_eff, rhs%nyq)
        if( nyq_eff <= 0 ) nyq_eff = ubound(lhs%cmplx_plane,1)
        hmin = max(pf*ceil_div(lbound(lhs%cmplx_plane,1),pf), pf*ceil_div(-nyq_eff,pf))
        hmax = min(pf*floor_div(ubound(lhs%cmplx_plane,1),pf), pf*floor_div(nyq_eff,pf))
        kmin = max(pf*ceil_div(lbound(lhs%cmplx_plane,2),pf), pf*ceil_div(-nyq_eff,pf))
        kmax = min(0, pf*floor_div(nyq_eff,pf))
        ! Integer form of the shell test below. For integer x >= 0 and integer n >= 0,
        ! nint(sqrt(x)) > n  <=>  sqrt(x) >= n+0.5  <=>  x >= n^2+n+0.25  <=>  x > n*(n+1),
        ! so the disc gate selects exactly the same samples without a square root and a round per
        ! element. This routine is called d_tilde*(d_tilde+1)/2 times per particle in the reduced solve.
        nyq_disk = nyq_eff * (nyq_eff + 1)
        do k = kmin, kmax, pf
            ! the k=0 line is its own Friedel mate, so only h<=0 there, or it is counted twice
            h_hi = hmax
            if( k == 0 ) h_hi = 0
            k_sq = k*k
            if( k_sq > nyq_disk ) cycle
            do h = hmin, h_hi, pf
                if( h*h + k_sq > nyq_disk ) cycle
                if( hlf /= 0 )then
                    par = modulo((h/pf) + (k/pf), 2) + 1
                    if( par /= hlf ) cycle
                endif
                acc = acc + conjg(cmplx(lhs%cmplx_plane(h,k),kind=dp)) * cmplx(rhs%cmplx_plane(h,k),kind=dp)
            end do
        end do
        val = acc
    end function cov_herm_inner

    !> Closed-form per-particle ML contrast a = <T mu, y> / <T mu, T mu>, clamped to a sane range.
    real function particle_contrast( mean_fpl, fpl )
        type(fplane_type), intent(in) :: mean_fpl, fpl
        real(dp) :: emm, emy
        if( COV_UNIT_CONTRAST .and. .not. cov_fit_contrast_rt )then
            particle_contrast = 1.0
            return
        endif
        emm = real(cov_herm_inner(mean_fpl, mean_fpl), dp)
        emy = real(cov_herm_inner(mean_fpl, fpl), dp)
        particle_contrast = real(max(0.1d0, min(5.0d0, emy / max(emm, DTINY))))
    end function particle_contrast

    !> Whitened self-power of a plane in cov_herm_inner's index convention. The reduced-solve debias
    !! assumes E[Re(b_q)Re(b_r)]_noise = 0.5*sig2*G with b and G from cov_herm_inner, so the sig2 fed to
    !! it has to be the noise variance measured in THAT convention; this is what the log compares against.
    subroutine cov_herm_selfpower( fpl, pw, cnt )
        type(fplane_type), intent(in)  :: fpl
        real(dp),          intent(out) :: pw, cnt
        integer     :: h, k, hmin, hmax, kmin, kmax, nyq_eff, pf, h_hi, nyq_disk, k_sq
        complex(dp) :: c
        pf  = OSMPL_PAD_FAC
        pw  = 0.d0; cnt = 0.d0
        nyq_eff = fpl%nyq
        if( nyq_eff <= 0 ) nyq_eff = ubound(fpl%cmplx_plane,1)
        hmin = max(pf*ceil_div(lbound(fpl%cmplx_plane,1),pf), pf*ceil_div(-nyq_eff,pf))
        hmax = min(pf*floor_div(ubound(fpl%cmplx_plane,1),pf), pf*floor_div(nyq_eff,pf))
        kmin = max(pf*ceil_div(lbound(fpl%cmplx_plane,2),pf), pf*ceil_div(-nyq_eff,pf))
        kmax = min(0, pf*floor_div(nyq_eff,pf))
        nyq_disk = nyq_eff * (nyq_eff + 1)   ! see cov_herm_inner for why this replaces nint(sqrt(.))
        do k = kmin, kmax, pf
            h_hi = hmax
            if( k == 0 ) h_hi = 0
            k_sq = k*k
            if( k_sq > nyq_disk ) cycle
            do h = hmin, h_hi, pf
                if( h*h + k_sq > nyq_disk ) cycle
                c  = cmplx(fpl%cmplx_plane(h,k), kind=dp)
                pw = pw + real(c*conjg(c), dp)
                cnt= cnt + 1.d0
            end do
        end do
    end subroutine cov_herm_selfpower

    !> Signal-free whitened-noise variance per coefficient, from the high-frequency shells of a residual
    !! plane (where conformational signal is negligible).
    subroutine plane_hf_power( fpl, nyq, frac, pw, cnt )
        type(fplane_type), intent(in)  :: fpl
        integer,           intent(in)  :: nyq
        real,              intent(in)  :: frac
        real(dp),          intent(out) :: pw, cnt
        integer  :: pf, h, k, hmin, hmax, kmin, kmax, sh_lo, sh
        complex(dp) :: c
        pf   = OSMPL_PAD_FAC
        sh_lo= nint(frac*real(nyq))
        pw   = 0.d0; cnt = 0.d0
        hmin = pf*ceil_div(lbound(fpl%cmplx_plane,1),pf); hmax = pf*floor_div(ubound(fpl%cmplx_plane,1),pf)
        kmin = pf*ceil_div(lbound(fpl%cmplx_plane,2),pf); kmax = min(0, pf*floor_div(nyq,pf))
        do k = kmin, kmax, pf
            do h = hmin, hmax, pf
                sh = nint(sqrt(real((h/pf)**2 + (k/pf)**2)))
                if( sh < sh_lo .or. sh > nyq ) cycle
                c  = cmplx(fpl%cmplx_plane(h,k), kind=dp)
                pw = pw + real(c*conjg(c), dp)
                cnt= cnt + 1.d0
            end do
        end do
    end subroutine plane_hf_power

    !> Accumulate whitened plane power and voxel count per radial shell (native indices) into the passed
    !! shell-indexed arrays.
    subroutine plane_shell_power_accum( fpl, nyq, pw_sh, cnt_sh )
        type(fplane_type), intent(in)    :: fpl
        integer,           intent(in)    :: nyq
        real(dp),          intent(inout) :: pw_sh(0:), cnt_sh(0:)
        integer  :: pf, h, k, hmin, hmax, kmin, kmax, sh
        complex(dp) :: c
        pf   = OSMPL_PAD_FAC
        hmin = pf*ceil_div(lbound(fpl%cmplx_plane,1),pf); hmax = pf*floor_div(ubound(fpl%cmplx_plane,1),pf)
        kmin = pf*ceil_div(lbound(fpl%cmplx_plane,2),pf); kmax = min(0, pf*floor_div(nyq,pf))
        do k = kmin, kmax, pf
            do h = hmin, hmax, pf
                sh = nint(sqrt(real((h/pf)**2 + (k/pf)**2)))
                if( sh > nyq ) cycle
                c = cmplx(fpl%cmplx_plane(h,k), kind=dp)
                pw_sh(sh)  = pw_sh(sh)  + real(c*conjg(c), dp)
                cnt_sh(sh) = cnt_sh(sh) + 1.d0
            end do
        end do
    end subroutine plane_shell_power_accum

    !>  Whitened residual plane r_i = data - transfer * P_i mu (mean_fpl already carries
    !!  the CTF-amplitude-weighted mean projection).
    subroutine form_residual_plane( fpl, mean_fpl, resid, contrast )
        type(fplane_type), intent(in)    :: fpl, mean_fpl
        type(fplane_type), intent(inout) :: resid
        real, optional,    intent(in)    :: contrast
        integer :: h1, h2, k1, k2
        real    :: a
        a = 1.0; if( present(contrast) ) a = contrast
        h1 = lbound(fpl%cmplx_plane,1); h2 = ubound(fpl%cmplx_plane,1)
        k1 = lbound(fpl%cmplx_plane,2); k2 = ubound(fpl%cmplx_plane,2)
        if( .not. allocated(resid%cmplx_plane) ) allocate(resid%cmplx_plane(h1:h2,k1:k2))
        resid%cmplx_plane = fpl%cmplx_plane - a*mean_fpl%cmplx_plane
        resid%frlims = fpl%frlims
        resid%nyq    = fpl%nyq
    end subroutine form_residual_plane


    subroutine cleanup_plane( fpl )
        type(fplane_type), intent(inout) :: fpl
        if( allocated(fpl%cmplx_plane)    ) deallocate(fpl%cmplx_plane)
        if( allocated(fpl%ctfsq_plane)    ) deallocate(fpl%ctfsq_plane)
        if( allocated(fpl%transfer_plane) ) deallocate(fpl%transfer_plane)
    end subroutine cleanup_plane

    ! ============================ SELF-CONTAINED TESTS ============================ Synthetic inputs only.


end module simple_flex_pca_columns
