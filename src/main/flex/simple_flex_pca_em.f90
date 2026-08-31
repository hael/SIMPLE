!@descr: EM estimation of the low-rank flex_pca covariance model: basis fit, latent embedding and mean estimation
module simple_flex_pca_em
use simple_core_module_api
use simple_builder,         only: builder
use simple_image,           only: image
use simple_parameters,      only: parameters
use simple_reconstructor,   only: reconstructor
use simple_kbinterpol,      only: kbinterpol
use simple_gridding,        only: prep3D_inv_kbenvelope4mul
use simple_linalg,          only: jacobi, eigsrt
use simple_math,            only: ceil_div, floor_div
use simple_srch_sort_loc,   only: hpsort
use simple_math_ft,         only: resample_sigma2
use simple_ftiter,          only: ftiter
use simple_flex_pca_polar,  only: polar_grid_t
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


interface

    ! ===== implemented in simple_flex_pca_em_env =====

    real(dp) module function corr_dp( a, b, n ) result( r )
        integer,  intent(in) :: n
        real(dp), intent(in) :: a(n), b(n)
    end function corr_dp

    logical module function cov_env_int_off( name ) result(off)
        character(len=*), intent(in) :: name
    end function cov_env_int_off

    logical module function cov_env_int_on( name ) result(on)
        character(len=*), intent(in) :: name
    end function cov_env_int_on

    pure integer module function cov_signal_rank( eval, n ) result( d )
        integer,  intent(in) :: n
        real(dp), intent(in) :: eval(n)          !< DESCENDING eigenvalues
    end function cov_signal_rank

    module subroutine cov_stage_subsample( build, pinds, nptcls, nparts, maxtot, env_stride, env_max, &
        &label, spinds, nsel )
        type(builder),        intent(inout) :: build
        integer,              intent(in)    :: pinds(:), nptcls, nparts, maxtot
        character(len=*),     intent(in)    :: env_stride, env_max, label
        integer, allocatable, intent(out)   :: spinds(:)
        integer,              intent(out)   :: nsel
    end subroutine cov_stage_subsample

    module subroutine write_zhalf_replicates( zhalf, prior, nptcls, ncomp )
        real(dp), intent(in) :: zhalf(nptcls,ncomp,2), prior(ncomp)
        integer,  intent(in) :: nptcls, ncomp
    end subroutine write_zhalf_replicates

    module subroutine cov_env_int_pub( name, val )
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: val
    end subroutine cov_env_int_pub

    logical module function cov_env_is_set( name )
        character(len=*), intent(in) :: name
    end function cov_env_is_set

    module subroutine cov_env_int( name, val )
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: val
    end subroutine cov_env_int

    logical module function cov_packed_cgsolve() result( packed )
    end function cov_packed_cgsolve

    pure real(dp) module function cov_accum_bytes( d, packed ) result( nbytes )
        integer, intent(in) :: d
        logical, intent(in) :: packed
    end function cov_accum_bytes

    pure integer module function cov_dim_budget( packed ) result( d )
        logical, intent(in) :: packed
    end function cov_dim_budget

    module subroutine cov_dev_prep_start( params, build, l_on )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        logical,           intent(out)   :: l_on
    end subroutine cov_dev_prep_start

    module subroutine cov_dev_prep_stop( l_on )
        logical, intent(in) :: l_on
    end subroutine cov_dev_prep_stop

    module subroutine cov_init_right_kernel_width
    end subroutine cov_init_right_kernel_width

    module subroutine map_sampling_precision( Gtil, prior, n, Qout )
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: Gtil(n,n), prior(n)
        real(dp), intent(out) :: Qout(n,n)
    end subroutine map_sampling_precision

    ! ===== implemented in simple_flex_pca_em_fit =====

    module subroutine build_covariance_eigenbasis( params, build, mean_rec, pinds, nptcls, &
        &col_sep, neigs_req, basis_recs, eigvals, ncomp_out, sig2_out, basis_imgs, fprefix )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, col_sep, neigs_req
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        real(dp),            allocatable, intent(out) :: eigvals(:)
        integer,             intent(out)   :: ncomp_out
        real(dp),            intent(out)   :: sig2_out
        type(image), allocatable, optional, intent(out) :: basis_imgs(:)
        character(len=*),        optional, intent(in)  :: fprefix
    end subroutine build_covariance_eigenbasis

    module subroutine init_basis_datafree( params, build, mean_rec, pinds, nptcls, col_sep, neigs_req, &
        &basis_recs, eigvals, ncomp_out, sig2_out )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, col_sep, neigs_req
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        real(dp),            allocatable, intent(out) :: eigvals(:)
        integer,             intent(out)   :: ncomp_out
        real(dp),            intent(out)   :: sig2_out
    end subroutine init_basis_datafree

    module subroutine em_calibrate_noise_prior( params, build, mean_rec, basis_recs, ncomp, pinds, nptcls, &
        &sig2_out, gam0 )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), intent(inout) :: basis_recs(:)
        integer,             intent(in)    :: ncomp, pinds(:), nptcls
        real(dp),            intent(out)   :: sig2_out, gam0
    end subroutine em_calibrate_noise_prior

    module subroutine save_probe_state( ncomp, eigvals, sig2_eff )
        integer,  intent(in) :: ncomp
        real(dp), intent(in) :: eigvals(:), sig2_eff
    end subroutine save_probe_state

    module subroutine load_probe_state( ncomp, eigvals, sig2_eff )
        integer,               intent(out) :: ncomp
        real(dp), allocatable, intent(out) :: eigvals(:)
        real(dp),              intent(out) :: sig2_eff
    end subroutine load_probe_state

    module subroutine load_probe_basis( params, build, ncomp, basis_recs )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        integer,             intent(in)    :: ncomp
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
    end subroutine load_probe_basis

    module subroutine select_frequencies_lowfreq( params, ncols_req, col_sep, col_hkl, ncol )
        class(parameters),    intent(in)  :: params
        integer,              intent(in)  :: ncols_req, col_sep
        integer, allocatable, intent(out) :: col_hkl(:,:)
        integer,              intent(out) :: ncol
    end subroutine select_frequencies_lowfreq

    module subroutine pick_next_lowfreq( cand, ncand, chosen, nchosen, sep, best )
        integer, intent(in)  :: cand(:,:), ncand, chosen(:,:), nchosen, sep
        integer, intent(out) :: best
    end subroutine pick_next_lowfreq

    module function covariance_kfromto( params ) result( kfromto )
        class(parameters), intent(in) :: params
        integer :: kfromto(2)
    end function covariance_kfromto

    module subroutine probe_external_basis( params, build, mean_rec, pinds, nptcls, eigdir, neigs, eigvals, &
        &sig2_eff, probe_prefix, nprobe )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, neigs, nprobe
        character(len=*),    intent(in)    :: eigdir, probe_prefix
        real(dp),            intent(in)    :: eigvals(:), sig2_eff
    end subroutine probe_external_basis

    module subroutine init_basis_reconstructor( params, build, rec )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: rec
    end subroutine init_basis_reconstructor

    module function cov_image_mask_radius( params ) result( r )
        class(parameters), intent(in) :: params
        real :: r
    end function cov_image_mask_radius

    logical module function cov_dump_acc_on() result( on )
    end function cov_dump_acc_on

    ! ===== implemented in simple_flex_pca_em_mean =====

    module subroutine estimate_covariance_mean( params, build, mean_rec, pinds, nptcls )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
    end subroutine estimate_covariance_mean

    module subroutine estimate_mean_from_data( params, build, mean_rec, pinds, nptcls )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
    end subroutine estimate_mean_from_data

    module subroutine apply_cached_mean_scale( params, mean_rec )
        class(parameters),   intent(inout) :: params
        type(reconstructor), intent(inout) :: mean_rec
    end subroutine apply_cached_mean_scale

    module subroutine estimate_mean_scale( params, build, mean_rec, pinds, nptcls )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
    end subroutine estimate_mean_scale

    module subroutine plane_shell_cross_accum( mean_fpl, fpl, nyq, my_sh, mm_sh )
        type(fplane_type), intent(in)    :: mean_fpl, fpl
        integer,           intent(in)    :: nyq
        real(dp),          intent(inout) :: my_sh(0:), mm_sh(0:)
    end subroutine plane_shell_cross_accum

    module subroutine init_mean_reconstructor( params, build, mean_rec )
        class(parameters),  intent(inout) :: params
        type(builder),      intent(inout) :: build
        type(reconstructor),intent(inout) :: mean_rec
    end subroutine init_mean_reconstructor

    module subroutine form_reconstruction_plane( fpl, num )
        type(fplane_type), intent(in)    :: fpl
        type(fplane_type), intent(inout) :: num
    end subroutine form_reconstruction_plane

    module subroutine cov_herm_sample_list( ref, slist, nsamp )
        type(fplane_type), intent(in)  :: ref
        integer,           intent(out) :: slist(:,:)
        integer,           intent(out) :: nsamp
    end subroutine cov_herm_sample_list

    module function cov_herm_inner( lhs, rhs, half ) result( val )
        type(fplane_type), intent(in) :: lhs, rhs
        integer, optional, intent(in) :: half
        complex(dp) :: val
    end function cov_herm_inner

    real module function particle_contrast( mean_fpl, fpl )
        type(fplane_type), intent(in) :: mean_fpl, fpl
    end function particle_contrast

    module subroutine cov_herm_selfpower( fpl, pw, cnt )
        type(fplane_type), intent(in)  :: fpl
        real(dp),          intent(out) :: pw, cnt
    end subroutine cov_herm_selfpower

    module subroutine plane_hf_power( fpl, nyq, frac, pw, cnt )
        type(fplane_type), intent(in)  :: fpl
        integer,           intent(in)  :: nyq
        real,              intent(in)  :: frac
        real(dp),          intent(out) :: pw, cnt
    end subroutine plane_hf_power

    module subroutine plane_shell_power_accum( fpl, nyq, pw_sh, cnt_sh )
        type(fplane_type), intent(in)    :: fpl
        integer,           intent(in)    :: nyq
        real(dp),          intent(inout) :: pw_sh(0:), cnt_sh(0:)
    end subroutine plane_shell_power_accum

    module subroutine form_residual_plane( fpl, mean_fpl, resid, contrast )
        type(fplane_type), intent(in)    :: fpl, mean_fpl
        type(fplane_type), intent(inout) :: resid
        real, optional,    intent(in)    :: contrast
    end subroutine form_residual_plane

    module subroutine cleanup_plane( fpl )
        type(fplane_type), intent(inout) :: fpl
    end subroutine cleanup_plane

    ! ===== implemented in simple_flex_pca_em_basis =====

    module subroutine basis_to_real_representatives( params, work, colvol, ncol, lb, ub, realvols, nreal )
        class(parameters),   intent(inout) :: params
        type(reconstructor), intent(inout) :: work
        complex,             intent(in)    :: colvol(:,:,:,:)
        integer,             intent(in)    :: ncol, lb(3), ub(3)
        type(image), allocatable, intent(out) :: realvols(:)
        integer,                  intent(out) :: nreal
    end subroutine basis_to_real_representatives

    module subroutine realize_hermitian_volume( params, work, vherm, gridcorr_img, energy )
        class(parameters),   intent(in)    :: params
        type(reconstructor), intent(inout) :: work
        complex,             intent(in)    :: vherm(:,:,:)
        type(image),         intent(inout) :: gridcorr_img
        real,                intent(out)   :: energy
    end subroutine realize_hermitian_volume

    module subroutine orthonormalize_representatives( params, build, realvols, nreal, utilde, utilde_real, d_tilde, svals, &
        &nptcls_basis )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(image),         intent(inout) :: realvols(:)
        integer,             intent(in)    :: nreal
        type(reconstructor), allocatable, intent(out) :: utilde(:)
        type(image),         allocatable, intent(out) :: utilde_real(:)
        integer,             intent(out)   :: d_tilde
        real(dp), allocatable, optional, intent(out) :: svals(:)
        integer, optional, intent(in) :: nptcls_basis
    end subroutine orthonormalize_representatives

    module subroutine align_basis_to_reference( ref_imgs, nref_c, tgt_imgs, ntgt_c, M, svals )
        integer,     intent(in)    :: nref_c, ntgt_c
        type(image), intent(inout) :: ref_imgs(nref_c), tgt_imgs(ntgt_c)
        real(dp), allocatable, intent(out) :: M(:,:), svals(:)
    end subroutine align_basis_to_reference

    module subroutine bag_basis_pool( imgs_a, na_c, eig_a, imgs_b, nb_c, eig_b, ncomp_out, pooled, eig_pooled )
        integer,     intent(in)    :: na_c, nb_c, ncomp_out
        type(image), intent(inout) :: imgs_a(na_c), imgs_b(nb_c)
        real(dp),    intent(in)    :: eig_a(na_c), eig_b(nb_c)
        type(image), allocatable, intent(out) :: pooled(:)
        real(dp),    allocatable, intent(out) :: eig_pooled(:)
    end subroutine bag_basis_pool

    module subroutine basis_recs_from_images( params, build, imgs, ncomp, basis_recs )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        integer,             intent(in)    :: ncomp
        type(image),         intent(inout) :: imgs(ncomp)
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
    end subroutine basis_recs_from_images

    ! ===== implemented in simple_flex_pca_em_pose =====

    logical module function cov_polar_enabled()
    end function cov_polar_enabled

    logical module function cov_polar_embed_enabled()
    end function cov_polar_embed_enabled

    integer module function cov_polar_ndir( nptcls )
        integer, intent(in) :: nptcls
    end function cov_polar_ndir

    module subroutine embed_accumulate_polar( params, build, mean_rec, basis_recs, ncomp, eigvals, sig2, &
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
    end subroutine embed_accumulate_polar

    pure real module function polar_inplane_deg( ca, sa ) result( d )
        real, intent(in) :: ca, sa
    end function polar_inplane_deg

    integer module function cov_pose_mode()
    end function cov_pose_mode

    module subroutine polar_pose_refine_one( pg, fpl, Usall, Cfall, Cm0all, c00all, nsamp2, ncomp, nk, &
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
    end subroutine polar_pose_refine_one

    module subroutine polar_dir_tables( Cfall, Cm0all, c00all, ncomp, nk, ndir, jd, wr, prior, sig2, &
        &acon, Ach, cvec, e_mm, logdet, info )
        integer,  intent(in)  :: ncomp, nk, ndir, jd
        real(dp), intent(in)  :: Cfall(ncomp*ncomp,nk,ndir), Cm0all(ncomp,nk,ndir), c00all(nk,ndir)
        real,     intent(in)  :: wr(nk)
        real(dp), intent(in)  :: prior(ncomp), sig2, acon
        real(dp), intent(out) :: Ach(ncomp,ncomp), cvec(ncomp), e_mm, logdet
        integer,  intent(out) :: info
    end subroutine polar_dir_tables

    module subroutine polar_pose_score_halves( pg, xw, Us, nsamp2, ncomp, Ach, cvec, sig2, acon, e_mm, &
        &xtrial, uvec, sc1, sc2, scfull )
        type(polar_grid_t), intent(in)    :: pg
        complex,            intent(in)    :: xw(:)
        integer,            intent(in)    :: nsamp2, ncomp
        real,               intent(in)    :: Us(nsamp2,0:ncomp)
        real(dp),           intent(in)    :: Ach(ncomp,ncomp), cvec(ncomp), sig2, acon, e_mm
        real,               intent(inout) :: xtrial(nsamp2)
        real(dp),           intent(inout) :: uvec(ncomp)
        real(dp),           intent(out)   :: sc1, sc2, scfull
    end subroutine polar_pose_score_halves

    real(dp) module function polar_pose_score( pg, xw, Us, nsamp2, ncomp, Ach, cvec, prior, sig2, acon, &
        &e_mm, xtrial, bq, uvec ) result( sc )
        type(polar_grid_t), intent(in)    :: pg
        complex,            intent(in)    :: xw(:)
        integer,            intent(in)    :: nsamp2, ncomp
        real,               intent(in)    :: Us(nsamp2,0:ncomp)
        real(dp),           intent(in)    :: Ach(ncomp,ncomp), cvec(ncomp), prior(ncomp), sig2
        real(dp),           intent(in)    :: acon, e_mm
        real,               intent(inout) :: xtrial(nsamp2), bq(0:ncomp)
        real(dp),           intent(inout) :: uvec(ncomp)
    end function polar_pose_score

    module subroutine cov_perturb_project_poses( build, pinds, nptcls )
        type(builder), intent(inout) :: build
        integer,       intent(in)    :: pinds(:), nptcls
    end subroutine cov_perturb_project_poses

    pure real module function pose_hash( i, k ) result( h )
        integer, intent(in) :: i, k
    end function pose_hash

    module subroutine polar_ring_gram( Us, ldu, ncomp, row0, nrow, Csp, Cout, Mout )
        integer,  intent(in)    :: ldu, ncomp, row0, nrow
        real,     intent(in)    :: Us(ldu,0:ncomp)
        real,     intent(inout) :: Csp(0:ncomp,0:ncomp)      !< caller-owned scratch
        real(dp), intent(out)   :: Cout(ncomp*ncomp), Mout(ncomp)
    end subroutine polar_ring_gram

    real(dp) module function polar_ring_selfpower( Us, ldu, row0, nrow )
        integer, intent(in) :: ldu, row0, nrow
        real,    intent(in) :: Us(ldu,0:*)
    end function polar_ring_selfpower

    real(dp) module function polar_self_energy( xws, wr, pg ) result( e )
        real,               intent(in) :: xws(:), wr(:)
        type(polar_grid_t), intent(in) :: pg
    end function polar_self_energy

    pure real(dp) module function sum_dp_safe( acc, n ) result( v )
        real(dp), intent(in) :: acc
        integer,  intent(in) :: n
    end function sum_dp_safe

    module subroutine align_halfplane_inplane( frlims, nyq_eff, src, ca, sa, dst )
        integer, intent(in)  :: frlims(3,2), nyq_eff
        complex, intent(in)  :: src(frlims(1,1):frlims(1,2), frlims(2,1):0)
        real,    intent(in)  :: ca, sa
        complex, intent(out) :: dst(frlims(1,1):frlims(1,2), frlims(2,1):0)
    end subroutine align_halfplane_inplane

    module subroutine polar_sample_particle_packed( fpl, pg, ca, sa, xws, wr, hfpw, hfcnt, tazim, xws1, xws2 )
        type(fplane_type),  intent(in)    :: fpl
        type(polar_grid_t), intent(in)    :: pg
        real,               intent(in)    :: ca, sa
        real,               intent(out)   :: xws(:)
        real,               intent(out)   :: wr(:)
        real(dp),           intent(inout) :: hfpw, hfcnt
        real,               intent(out)   :: tazim
        real, optional,     intent(out)   :: xws1(:), xws2(:)
    end subroutine polar_sample_particle_packed

    module subroutine project_fplane_mean_banded( rec, o, fpl_ref, fpl_out )
        type(reconstructor), intent(in)    :: rec
        class(ori),          intent(inout) :: o
        type(fplane_type),   intent(in)    :: fpl_ref
        type(fplane_type),   intent(inout) :: fpl_out
    end subroutine project_fplane_mean_banded

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
    end subroutine polar_hybrid_exact_accum

    module subroutine subtract_mean_banded( fpl, mean_fpl, a, rec_nyq )
        type(fplane_type), intent(inout) :: fpl
        type(fplane_type), intent(in)    :: mean_fpl
        real,              intent(in)    :: a
        integer,           intent(in)    :: rec_nyq
    end subroutine subtract_mean_banded

    ! ===== implemented in simple_flex_pca_em_embed =====

    module subroutine embed_latents_with_contrast( params, build, mean_rec, basis_recs, ncomp, eigvals, sig2_eff, &
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
        real(dp), optional,  intent(out)   :: rho_out(ncomp)
        logical,  optional,  intent(in)    :: stats_only
        logical,  optional,  intent(in)    :: from_parts
    end subroutine embed_latents_with_contrast

    ! ===== implemented in simple_flex_pca_em_iter =====

    module subroutine probe_subspace_iteration( params, build, mean_rec, basis_recs, eigvals, sig2_eff, &
        &pinds, nptcls, ncomp, niters )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), allocatable, intent(inout) :: basis_recs(:)
        real(dp),            allocatable, intent(inout) :: eigvals(:)
        real(dp),            intent(in)    :: sig2_eff
        integer,             intent(in)    :: pinds(:), nptcls, niters
        integer,             intent(inout) :: ncomp
    end subroutine probe_subspace_iteration

    ! ===== implemented in simple_flex_pca_em_solve =====

    pure module subroutine spd_logdet_dp( A, n, logdet, ok )
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: A(n,n)
        real(dp), intent(out) :: logdet
        logical,  intent(out) :: ok
    end subroutine spd_logdet_dp

    module subroutine deflate_against_basis( imgs, n, basis, nb )
        integer,     intent(in)    :: n, nb
        type(image), intent(inout) :: imgs(n), basis(nb)
    end subroutine deflate_against_basis

    module subroutine cross_half_subspace_angles( eimgs, oimgs, n, svals )
        integer,     intent(in)    :: n
        type(image), intent(inout) :: eimgs(n), oimgs(n)
        real(dp), allocatable, intent(out) :: svals(:)
    end subroutine cross_half_subspace_angles

    pure module function quad_form( M, z, n ) result( val )
        integer,  intent(in) :: n
        real(dp), intent(in) :: M(n,n), z(n)
        real(dp) :: val
    end function quad_form

    module subroutine spd_solve_dp( A, b, n )
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: A(n,n), b(n)
    end subroutine spd_solve_dp

    module subroutine probe_solve_ecm( n, G, b, c, myv, e_mm, prior_, sig2, nml, a, z_, Ainv_, ldA, lok, quad )
        integer,  intent(in)    :: n, nml
        real(dp), intent(in)    :: G(n,n), b(n), c(n), myv, e_mm, prior_(n), sig2
        real(dp), intent(inout) :: a
        real(dp), intent(out)   :: z_(n), Ainv_(n,n), ldA
        logical,  intent(out)   :: lok
        real(dp), intent(out)   :: quad
    end subroutine probe_solve_ecm

    module subroutine probe_solve_mix( n, kmix, G, b, c, myv, e_mm, nml, a, sig2, Ominv, Omxi, lpi, xiOx, &
            &zbar, Ainv_, dens_, ldA, lok, nll_add, sr_acc, sm_acc, smm_acc, sainv_acc )
        integer,  intent(in)    :: n, kmix, nml
        real(dp), intent(in)    :: G(n,n), b(n), c(n), myv, e_mm, sig2
        real(dp), intent(inout) :: a
        real(dp), intent(in)    :: Ominv(n,n), Omxi(n,kmix), lpi(kmix), xiOx(kmix)
        real(dp), intent(out)   :: zbar(n), Ainv_(n,n), dens_(n,n), ldA, nll_add
        logical,  intent(out)   :: lok
        real(dp), intent(inout) :: sr_acc(kmix), sm_acc(n,kmix), smm_acc(n,n,kmix), sainv_acc(n,n)
    end subroutine probe_solve_mix

    module subroutine mcfa_init( z, nptcls, ncomp, kmix, gam_sum, nval, xi, ppi, Om )
        integer,  intent(in)  :: nptcls, ncomp, kmix, nval
        real(dp), intent(in)  :: z(nptcls,ncomp), gam_sum(ncomp)
        real(dp), intent(out) :: xi(ncomp,kmix), ppi(kmix), Om(ncomp,ncomp)
    end subroutine mcfa_init

    module subroutine mcfa_condition( ncomp, diag_only, Om, Ominv, ldOm )
        integer,  intent(in)    :: ncomp
        logical,  intent(in)    :: diag_only
        real(dp), intent(inout) :: Om(ncomp,ncomp)
        real(dp), intent(out)   :: Ominv(ncomp,ncomp), ldOm
    end subroutine mcfa_condition

    module subroutine mcfa_mstep( ncomp, kmix, nval, sr, sm, smm, sai, &
        &pin_origin, ppi, xi, Om, Ominv, ldOm )
        integer,  intent(in)    :: ncomp, kmix, nval
        real(dp), intent(in)    :: sr(kmix), sm(ncomp,kmix)
        real(dp), intent(in)    :: smm(ncomp,ncomp,kmix), sai(ncomp,ncomp)
        logical,  intent(in)    :: pin_origin
        real(dp), intent(inout) :: ppi(kmix), xi(ncomp,kmix), Om(ncomp,ncomp)
        real(dp), intent(out)   :: Ominv(ncomp,ncomp), ldOm
    end subroutine mcfa_mstep

    module subroutine spd_inv_dp( A, Ainv, n )
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: A(n,n)
        real(dp), intent(out)   :: Ainv(n,n)
    end subroutine spd_inv_dp

end interface

end module simple_flex_pca_em
