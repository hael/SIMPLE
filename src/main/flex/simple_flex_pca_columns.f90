!@descr: covariance-column estimation, reduced covariance solve and latent embedding for flex_pca
module simple_flex_pca_columns
use simple_core_module_api
use simple_builder,          only: builder
use simple_image,            only: image
use simple_parameters,       only: parameters
use simple_reconstructor,    only: reconstructor
use simple_flex_pca_distr,   only: flex_pca_is_master, flex_pca_is_worker, flex_pca_nparts, &
    &flex_pca_run_stage, PCA_STAGE_SNR, PCA_STAGE_COLS, PCA_STAGE_SOLVE, PCA_STAGE_PROBE, PCA_STAGE_EMBED
use simple_flex_pca_parts,   only: flex_pca_part_fname, write_snr_part, reduce_snr_parts, &
    &write_cols_part, reduce_cols_parts, write_solve_part, reduce_solve_parts, &
    &write_probe_part, reduce_probe_parts, &
    &write_embed_stats_part, reduce_embed_zhalf_parts, read_embed_stats_part, &
    &write_mean_scale, read_mean_scale, &
    &write_columns_hkl, read_columns_hkl
use simple_kbinterpol,       only: kbinterpol
use simple_gridding,         only: prep3D_inv_instrfun4mul
use simple_linalg,           only: jacobi, eigsrt, svd_solve
use simple_math,             only: ceil_div, floor_div
use simple_srch_sort_loc,    only: hpsort
use simple_matcher_3Drec,    only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io,  only: discrete_read_imgbatch, prepimgbatch
use simple_flex_reconstructor_latent_ops, only: project_fplane_mean, project_fplanes_mean_basis, &
    &insert_planes_oversamp_coupled_batch_scaled
use simple_flex_projected_latent_model,   only: prep_imgs4projected_model, solve_coupled_basis_exp
implicit none
private
#include "simple_local_flags.inc"

public :: build_covariance_eigenbasis, embed_latents_with_contrast, estimate_covariance_mean
public :: svec_idx, apply_packed_A, cg_solve_packed
public :: probe_subspace_iteration, align_basis_to_reference, probe_external_basis
public :: cov_env_int_pub, save_probe_state
public :: test_flex_pca_svec_isometry, test_flex_pca_packed_solve

! Density observability floor, matching simple_image_arith::div_cmat_at_1 and the projected-latent coupled
! solve.
real(dp), parameter :: COV_DENSITY_FLOOR = 1.0d-6
! Relative ridge used ONLY for the covariance diagonal in the S.B SNR proxy, which runs before the S.C
! weights exist (Algorithm 1 precedes Algorithm 2).
real,     parameter :: COV_RIDGE_REL     = 5.0e-2
! Relative eigenvalue floor for retaining direct-column PCA components.
real(dp), parameter :: COV_EIG_REL_FLOOR = 1.0d-6
! Cap on the column-subspace dimension. The accumulation is a batched dsyrk on the Van Loan-Pitsianis
! rearrangement (see unrearrange_kron_selfsum), which needs ONE shared d^4 array regardless of thread count.
integer,  parameter :: COV_MAX_DTILDE    = 320
!> Total particles the probe refines the basis on, summed across processes. The basis has a fixed
!! parameter count, so this is a CAP rather than a fraction -- see the stride derivation below.
integer,  parameter :: COV_PROBE_MAX_PTCLS = 25000
!> probe stops when successive bases agree to this mean principal-angle cosine. 0.999 is tight
!! enough that the remaining rotation cannot move a state target, and on Ribosembly it fires at the
!! measured knee (iteration 2) rather than running the tuned count out.
real(dp), parameter :: COV_PROBE_CONV    = 0.97d0
! Memory budget for the shared A accumulator, in bytes.
real(dp), parameter :: COV_ATHR_BUDGET   = 8.0d9
! Accumulate the columns against the unscaled mean (a==1) rather than the per-particle ML contrast a_i.
! Subtracting a_i*T*mu also deletes the component of the conformational signal parallel to T*mu.
logical,  parameter :: COV_UNIT_CONTRAST  = .true.
! Grid-search the per-particle contrast in the embedding instead of using the closed-form estimate.
logical,  parameter :: COV_EMBED_CONTRAST_GRID = .false.
integer,  parameter :: COV_CG_MAXIT = 2000     ! CG iteration cap; convergence is reported, not assumed
real(dp), parameter :: COV_CG_TOL   = 1.d-10   ! relative residual target
integer,  parameter :: GRAM_DIAG_STRIDE = 200   ! subsample for the projected-Gram spectrum
integer,  parameter :: NCONTRAST_GRID = 50
real(dp), parameter :: A_GRID_HI = 2.0d0
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

    !> HeteroPCA (Zhang, Cai & Wu): delete the noise-contaminated diagonal of S and re-impute it from the
    !! rank-r reconstruction of the off-diagonals, iterating to a fixed point. Whatever the analytic debias
    !! Dmat = Sbb - 0.5*sig2_eff*SG fails to remove piles onto the diagonal, where it inflates the leading
    !! eigenvalue and tilts its eigenvector toward the lowest-frequency mode.
    subroutine heteropca_impute( S, n, r, niter )
        integer,           intent(in)    :: n, r
        real(dp),          intent(inout) :: S(n,n)
        integer, optional, intent(in)    :: niter
        real(dp), allocatable :: Mimp(:,:), evec(:,:), eval(:), work(:,:)
        integer  :: it, nit, i, j, k, nrot, rr
        real(dp) :: dmax
        if( n < 2 ) return
        rr  = max(1, min(r, n-1))
        nit = 30
        if( present(niter) ) nit = max(1, niter)
        allocate(Mimp(n,n), evec(n,n), eval(n), work(n,n))
        Mimp = S
        do i = 1, n
            Mimp(i,i) = 0.d0                      ! delete the noise-contaminated diagonal
        end do
        do it = 1, nit
            work = Mimp
            call jacobi(work, n, n, eval, evec, nrot)
            call eigsrt(eval, evec, n, n)
            dmax = 0.d0
            do i = 1, n
                ! diagonal of the rank-r reconstruction sum_k lambda_k v_k v_k'
                work(i,1) = 0.d0
                do k = 1, rr
                    work(i,1) = work(i,1) + eval(k)*evec(i,k)*evec(i,k)
                end do
                dmax = max(dmax, abs(work(i,1) - Mimp(i,i)))
            end do
            do i = 1, n
                Mimp(i,i) = work(i,1)             ! impute; off-diagonals untouched throughout
            end do
            if( dmax <= 1.d-12 ) exit
        end do
        write(logfhandle,'(A,I0,A,I0,A,ES11.3)') '>>> FLEX_PCA HeteroPCA diagonal imputation: rank=', &
            &rr,' iters=',it,' final diagonal change=',dmax
        call flush(logfhandle)
        S = Mimp
        deallocate(Mimp, evec, eval, work)
    end subroutine heteropca_impute

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

    !>  Write the two half-data latent solves so the per-particle error MODEL can be calibrated
    !!  against the error actually observed. The halves are disjoint checkerboard subsets of ONE
    !!  particle's own Fourier samples (see cov_herm_inner's `half` argument), so var(z1 - z2)
    !!  measures that particle's estimation error directly -- including model misspecification,
    !!  which the analytic posterior covariance cannot express. `prior` is written alongside because
    !!  the half solves are shrunk by it and any calibration has to undo that.
    !!  Gated on SIMPLE_COV_ZHALF; writes nothing and costs nothing when unset.
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
        &ncols_req, col_sep, neigs_req, basis_recs, eigvals, ncomp_out, sig2_out, basis_imgs, fprefix )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, ncols_req, col_sep, neigs_req
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
        !> probe-worker handoff, read back from the master's flex_pca_probe.txt
        real(dp),            allocatable :: eig_probe(:)
        real(dp),            allocatable :: zw(:,:), contrastw(:), precw(:,:,:), rew(:), rmew(:)
        real(dp) :: sig2_probe
        integer  :: ncomp_probe
        real(dp), allocatable :: vred(:,:)
        integer :: ncol, nreal, s, lb(3), ub(3), nyq_rec, d_tilde, q, directsvd
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
        ! EMBED WORKER. Same handoff as the probe worker -- the basis is on disk as flex_pca_pc*.mrc
        ! and flex_pca_probe.txt carries its dimension, prior variances and noise level -- but only ONE
        ! round, because the basis is final by now and does not change under the workers.
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
        if( flex_pca_is_worker() .and. params%stage == PCA_STAGE_SOLVE )then
            ! nothing upstream of the solve is needed: the basis is on disk
            call load_utilde_stack(params, build, utilde, utilde_real, d_tilde)
            call reduced_covariance_solve(params, build, mean_rec, utilde, d_tilde, pinds, nptcls, &
                &neigs_req, vred, eigvals, ncomp_out, sig2_out)
            do s = 1, d_tilde
                call utilde(s)%dealloc_rho; call utilde(s)%kill
                call utilde_real(s)%kill
            end do
            deallocate(utilde, utilde_real)
            if( allocated(vred) ) deallocate(vred)
            if( .not. allocated(eigvals) ) allocate(eigvals(0))
            allocate(basis_recs(0))
            ncomp_out = 0
            sig2_out  = 0._dp
            call work%dealloc_rho; call work%kill
            return
        endif
        if( flex_pca_is_worker() .and. params%stage == PCA_STAGE_COLS )then
            call read_columns_hkl(col_hkl, ncol, l_cols_ok)
            if( .not. l_cols_ok ) THROW_HARD('flex_pca worker found no flex_pca_columns.bin from the master')
        else if( trim(params%column_sampling) == 'lowfreq' )then
            call select_covariance_columns_lowfreq(params, ncols_req, col_sep, col_hkl, ncol)
        else
            call select_covariance_columns_snr(params, build, mean_rec, pinds, nptcls, &
                &lb, ub, nyq_rec, ncols_req, col_sep, col_hkl, ncol)
        endif
        if( flex_pca_is_worker() .and. params%stage == PCA_STAGE_SNR )then
            ! the SNR part file was written inside estimate_snr_volume; this round is done
            ncomp_out = 0
            sig2_out  = 0._dp
            allocate(eigvals(0))
            allocate(basis_recs(0))
            call work%dealloc_rho; call work%kill
            if( allocated(col_hkl) ) deallocate(col_hkl)
            return
        endif
        write(logfhandle,'(A,I0,A,A,A,I0)') '>>> FLEX_PCA selected covariance columns=',ncol, &
            &' sampling=',trim(params%column_sampling),' separation=',col_sep
        call flush(logfhandle)
        if( flex_pca_is_master() ) call write_columns_hkl(col_hkl, ncol)
        call build_column_lookup(col_hkl, ncol, lb, ub, col_lookup)
        allocate(Bcol_e(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol), source=cmplx(0.,0.))
        allocate(Bcol_o(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol), source=cmplx(0.,0.))
        allocate(Hcol_e(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol), source=0.)
        allocate(Hcol_o(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol), source=0.)
        call accumulate_covariance_columns(params, build, mean_rec, pinds, nptcls, &
            &col_hkl, col_lookup, ncol, lb, ub, nyq_rec, Bcol_e, Hcol_e, Bcol_o, Hcol_o)
        if( flex_pca_is_worker() )then
            ! the column part file was written inside accumulate_covariance_columns; this round is done
            ncomp_out = 0
            sig2_out  = 0._dp
            allocate(eigvals(0))
            allocate(basis_recs(0))
            deallocate(Bcol_e, Bcol_o, Hcol_e, Hcol_o, col_lookup, col_hkl)
            call work%dealloc_rho; call work%kill
            return
        endif
        ! regularized merge + half-column FSC diagnostics
        allocate(colvol(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol), col_fsc(ncol))
        t_blk = tic()
        call regularize_and_merge_columns(Bcol_e, Hcol_e, Bcol_o, Hcol_o, ncol, lb, colvol, col_fsc)
        write(logfhandle,'(A,F9.1)') '>>> FLEX_PCA STAGE merge_columns seconds=', toc(t_blk)
        call write_column_diagnostics(col_hkl, col_fsc, ncol)
        deallocate(Bcol_e, Bcol_o, Hcol_e, Hcol_o, col_lookup, col_hkl)
        ! Friedel-symmetrize each complex column into two real (cos/sin) volumes
        t_blk = tic()
        call columns_to_real_representatives(params, work, colvol, ncol, lb, ub, realvols, nreal)
        write(logfhandle,'(A,F9.1)') '>>> FLEX_PCA STAGE real_representatives seconds=', toc(t_blk)
        deallocate(colvol)
        call work%dealloc_rho; call work%kill
        write(logfhandle,'(A,I0)') '>>> FLEX_PCA real column representatives=',nreal
        call flush(logfhandle)
        t_blk = tic()
        call orthonormalize_representatives(params, build, realvols, nreal, utilde, utilde_real, d_tilde, svals)
        write(logfhandle,'(A,F9.1)') '>>> FLEX_PCA STAGE orthonormalize seconds=', toc(t_blk)
        do s = 1, nreal
            call realvols(s)%kill
        end do
        deallocate(realvols, col_fsc)
        write(logfhandle,'(A,I0)') '>>> FLEX_PCA column subspace dimension d_tilde=',d_tilde
        call flush(logfhandle)
        ! Persist the orthonormal column basis for the solve round. orthonormalize_representatives
        ! builds utilde(q) from utilde_real(q) with exactly set_rmat(.,.false.) -> fft -> expand_exp,
        ! so shipping the real-space stack lets a worker reproduce the reconstructors bit-for-bit
        ! rather than re-deriving the whole column pipeline.
        if( flex_pca_is_master() ) call write_utilde_stack(utilde_real, d_tilde)
        ! S.B reduced projected-covariance solve (eqs S.6-S.9):
        call reduced_covariance_solve(params, build, mean_rec, utilde, d_tilde, pinds, nptcls, &
            &neigs_req, vred, eigvals, ncomp_out, sig2_out)
        ! the basis handoff has served its purpose; do not leave d_tilde volumes in the run directory
        if( flex_pca_is_master() ) call cleanup_utilde_stack(d_tilde)
        ! SIMPLE_COV_DIRECTSVD=1 takes the basis straight from the regularised-column SVD, bypassing the
        ! reduced solve. The default path uses the columns only to define a SPAN and re-estimates the
        ! covariance inside it, which conditions better.
        directsvd = 0
        call cov_env_int('SIMPLE_COV_DIRECTSVD', directsvd)
        if( directsvd /= 0 )then
            ncomp_out = max(1, min(neigs_req, d_tilde))
            if( allocated(vred) ) deallocate(vred)
            if( allocated(eigvals) ) deallocate(eigvals)
            allocate(vred(d_tilde, ncomp_out), source=0._dp)
            allocate(eigvals(ncomp_out))
            do q = 1, ncomp_out
                vred(q,q)  = 1._dp
                eigvals(q) = svals(q)
            end do
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA DIRECT-SVD basis: ', ncomp_out, &
                &' components taken straight from the regularised-column SVD (reduced solve bypassed)'
            write(logfhandle,'(A,ES12.4,A,ES12.4,A,F7.3)') '>>>   svals max=', eigvals(1), '  min=', &
                &eigvals(ncomp_out), '  lam1/lam5=', real(eigvals(1)/max(eigvals(min(5,ncomp_out)),DTINY))
            call flush(logfhandle)
        endif
        t_blk = tic()
        call form_eigenbasis_from_reduced(params, build, mean_rec, utilde_real, d_tilde, vred, eigvals, &
            &ncomp_out, basis_recs, basis_imgs=basis_imgs, fprefix=fprefix)
        write(logfhandle,'(A,F9.1)') '>>> FLEX_PCA STAGE form_eigenbasis seconds=', toc(t_blk)
        do s = 1, d_tilde
            call utilde(s)%dealloc_rho; call utilde(s)%kill
            call utilde_real(s)%kill
        end do
        deallocate(utilde, utilde_real, vred)
    end subroutine build_covariance_eigenbasis

    !> Persist / restore the orthonormal column basis across the master->worker boundary for the solve
    !! round. Only the real-space volumes travel; the reconstructors are rebuilt with the same three
    !! calls orthonormalize_representatives uses, so the worker's utilde is identical to the master's.
    subroutine write_utilde_stack( utilde_real, d_tilde )
        type(image), intent(inout) :: utilde_real(:)
        integer,     intent(in)    :: d_tilde
        type(string) :: fname
        integer :: q, funit, io_stat
        ! one file per basis volume: image%write refuses to place a 3D image into a stack
        do q = 1, d_tilde
            fname = string(COV_UTILDE_FBODY)//int2str_pad(q,3)//MRC_EXT
            call utilde_real(q)%write(fname, del_if_exists=.true.)
            call fname%kill
        end do
        call fopen(funit, file=string(COV_UTILDE_META), action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_utilde_stack; meta', io_stat)
        write(funit,'(I0)') d_tilde
        call fclose(funit)
    end subroutine write_utilde_stack

    subroutine cleanup_utilde_stack( d_tilde )
        integer, intent(in) :: d_tilde
        type(string) :: fname
        integer :: q
        do q = 1, d_tilde
            fname = string(COV_UTILDE_FBODY)//int2str_pad(q,3)//MRC_EXT
            call del_file(fname)
            call fname%kill
        end do
        call del_file(string(COV_UTILDE_META))
    end subroutine cleanup_utilde_stack

    subroutine load_utilde_stack( params, build, utilde, utilde_real, d_tilde )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), allocatable, intent(out) :: utilde(:)
        type(image),         allocatable, intent(out) :: utilde_real(:)
        integer,                          intent(out) :: d_tilde
        type(string) :: fname
        integer :: q, funit, io_stat
        if( .not. file_exists(string(COV_UTILDE_META)) )then
            THROW_HARD('flex_pca worker found no '//COV_UTILDE_META//' from the master')
        endif
        call fopen(funit, file=string(COV_UTILDE_META), action='READ', status='OLD', iostat=io_stat)
        call fileiochk('load_utilde_stack; meta', io_stat)
        read(funit,*) d_tilde
        call fclose(funit)
        if( d_tilde < 1 ) THROW_HARD('invalid cached column-subspace dimension')
        allocate(utilde(d_tilde), utilde_real(d_tilde))
        do q = 1, d_tilde
            fname = string(COV_UTILDE_FBODY)//int2str_pad(q,3)//MRC_EXT
            if( .not. file_exists(fname) ) THROW_HARD('flex_pca worker missing '//fname%to_char())
            call utilde_real(q)%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call utilde_real(q)%read(fname)
            call fname%kill
            call init_column_reconstructor(params, build, utilde(q))
            call utilde(q)%set_rmat(utilde_real(q)%get_rmat(), .false.)
            call utilde(q)%fft
            call utilde(q)%expand_exp
        end do
        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA worker restored the column basis, d_tilde=',d_tilde, &
            &' from the master'
        call flush(logfhandle)
    end subroutine load_utilde_stack


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

    !> SNR-greedy column selection (proposal 5.2, Algorithm 1).
    subroutine select_covariance_columns_snr( params, build, mean_rec, pinds, nptcls, &
        &lb, ub, nyq_rec, ncols_req, col_sep, col_hkl, ncol )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, lb(3), ub(3), nyq_rec, ncols_req, col_sep
        integer, allocatable, intent(out)  :: col_hkl(:,:)
        integer,              intent(out)  :: ncol
        real, allocatable :: snr(:,:,:)
        integer :: kfromto(2), kmax, kmax_sq, h, k, l, sep, r_sq, target, hn, kn, ln, bh, bk, bl, s
        real    :: best
        logical :: taken
        call estimate_snr_volume(params, build, mean_rec, pinds, nptcls, lb, ub, nyq_rec, snr)
        kfromto = covariance_kfromto(params)
        kmax    = max(2, kfromto(2))
        kmax_sq = kmax*(kmax+1)
        sep     = max(1, col_sep)
        target  = max(1, ncols_req)
        allocate(col_hkl(3, target))
        ncol = 0
        do
            ! find the highest-SNR unclaimed candidate in the canonical Hermitian half
            best = -huge(1.0); bh = 0; bk = 0; bl = 0
            do l = -kmax, kmax
                do k = -kmax, kmax
                    do h = 0, kmax
                        r_sq = h*h + k*k + l*l
                        if( r_sq == 0 .or. r_sq > kmax_sq ) cycle
                        if( h == 0 )then
                            if( k < 0 ) cycle
                            if( k == 0 .and. l < 0 ) cycle
                        endif
                        if( snr(h,k,l) <= best ) cycle
                        ! reject if within sep of an already-chosen column (or its conjugate)
                        taken = .false.
                        do s = 1, ncol
                            if( (h-col_hkl(1,s))**2+(k-col_hkl(2,s))**2+(l-col_hkl(3,s))**2 < sep*sep )taken=.true.
                            if( (h+col_hkl(1,s))**2+(k+col_hkl(2,s))**2+(l+col_hkl(3,s))**2 < sep*sep )taken=.true.
                            if( taken ) exit
                        end do
                        if( taken ) cycle
                        best = snr(h,k,l); bh = h; bk = k; bl = l
                    end do
                end do
            end do
            if( best <= -huge(1.0)*0.5 ) exit
            ncol = ncol + 1
            col_hkl(:,ncol) = [bh,bk,bl]
            if( ncol >= target ) exit
        end do
        if( ncol < 1 ) THROW_HARD('flex_pca SNR selection produced no columns')
        ! Report the radial placement of the selected columns
        block
            real :: rad, rmin, rmax, rmean
            integer :: s2
            rmin = huge(1.0); rmax = 0.; rmean = 0.
            do s2 = 1, ncol
                rad   = sqrt(real(sum(col_hkl(:,s2)**2)))
                rmin  = min(rmin, rad); rmax = max(rmax, rad); rmean = rmean + rad
            end do
            rmean = rmean / real(ncol)
            write(logfhandle,'(A,I0,A,F6.1,A,F6.1,A,F6.1,A,I0)') '>>> FLEX_PCA selected ',ncol, &
                &' columns, |xi| min=',rmin,' mean=',rmean,' max=',rmax,' of band kmax=',kmax
            if( rmean > 0.5*real(kmax) )then
                write(logfhandle,'(A)') '>>> FLEX_PCA WARNING: the selected columns sit in the OUTER &
                    &half of the covariance band. The conformational covariance is a low-frequency quantity; &
                    &columns there carry noise and the eigenbasis will not span the motion. Suspect the SNR proxy.'
            endif
            call flush(logfhandle)
        end block
        deallocate(snr)
    end subroutine select_covariance_columns_snr

    !> Per-voxel signal-variance (SNR proxy) from the whitened adjoint residuals, backprojected and
    !! normalized by sampling density, with the outer-shell noise floor subtracted.
    subroutine estimate_snr_volume( params, build, mean_rec, pinds, nptcls, lb, ub, nyq_rec, snr )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, lb(3), ub(3), nyq_rec
        real, allocatable,   intent(out)   :: snr(:,:,:)
        type(reconstructor) :: scratch
        type(fplane_type), allocatable :: fpls(:)
        type(fplane_type) :: mean_fpl, adj_fpl
        type(ori) :: orientation
        real, allocatable :: var_acc(:,:,:), dens_acc(:,:,:), hi(:)
        real(dp) :: floor_noise
        integer :: batchlims(2), batchsz, ibatch, i, iptcl, pf, h, k, l, nhi, sh, nyq_use
        real    :: inv_pf2, av
        complex :: cval
        integer(timer_int_kind) :: t_phase
        pf = OSMPL_PAD_FAC; inv_pf2 = 1.0/real(pf*pf)
        allocate(var_acc(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), source=0.)
        allocate(dens_acc(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), source=0.)
        ! DISTRIBUTED: the master fans the particle range out and reduces the two accumulators; the
        ! finalisation below (noise floor, normalisation) then runs once, on the global sums.
        if( flex_pca_is_master() )then
            call flex_pca_run_stage(PCA_STAGE_SNR, 'SNR variance')
            call reduce_snr_parts(params, flex_pca_nparts(), var_acc, dens_acc, lb, ub)
            goto 100
        endif
        call init_column_reconstructor(params, build, scratch)
        call mean_rec%expand_exp
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        t_phase = tic()
        write(logfhandle,'(A)') '>>> FLEX_PCA SNR VARIANCE ESTIMATION'
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
                call project_fplane_mean(mean_rec, orientation, fpls(i), mean_fpl, apply_ctf_amp=.true.)
                call form_adjoint_residual_plane(fpls(i), mean_fpl, adj_fpl, particle_contrast(mean_fpl, fpls(i)))
                call scratch%reset_exp
                call scratch%insert_plane_oversamp(build%pgrpsyms, orientation, adj_fpl)
                ! One threaded sweep instead of two serial whole-array expressions; each cell is
                ! independent, so this is elementwise-identical to the array form it replaces.
                !$omp parallel do default(shared) private(h,k,l,cval) collapse(2) &
                !$omp& schedule(static) proc_bind(close)
                do l = lb(3), ub(3)
                    do k = lb(2), ub(2)
                        do h = lb(1), ub(1)
                            cval = scratch%cmat_exp(h,k,l)
                            var_acc(h,k,l)  = var_acc(h,k,l) + (real(cval)*inv_pf2)**2 &
                                &                            + (aimag(cval)*inv_pf2)**2
                            dens_acc(h,k,l) = dens_acc(h,k,l) + scratch%rho_exp(h,k,l)
                        end do
                    end do
                end do
                !$omp end parallel do
            end do
        end do
        write(logfhandle,'(A,F8.1)') '>>> FLEX_PCA SNR VARIANCE SECONDS: ', toc(t_phase)
        call flush(logfhandle)
        call orientation%kill
        call cleanup_plane(mean_fpl); call cleanup_plane(adj_fpl)
        call scratch%dealloc_rho; call scratch%kill
        call cleanup_rec_buffers(build, fpls)
        if( flex_pca_is_worker() )then
            ! this part's contribution only; the master sums the parts and finalises
            call write_snr_part(flex_pca_part_fname('snr', params%part, params%numlen), &
                &var_acc, dens_acc, lb, ub)
            deallocate(var_acc, dens_acc)
            allocate(snr(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), source=0.)
            return
        endif
100     continue
        ! variance per unit coverage
        allocate(snr(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), source=0.)
        where( dens_acc > real(COV_DENSITY_FLOOR) ) snr = var_acc / dens_acc
        ! robust noise floor from the outer (signal-poor) shells, subtracted so the
        ! ranking reflects excess (conformational) variance
        nyq_use = nyq_rec
        allocate(hi(0))
        nhi = 0
        do l = lb(3), ub(3); do k = lb(2), ub(2); do h = lb(1), ub(1)
            sh = nint(sqrt(real(h*h+k*k+l*l)))
            if( sh >= nint(0.75*nyq_use) .and. sh <= nyq_use .and. dens_acc(h,k,l) > real(COV_DENSITY_FLOOR) )then
                nhi = nhi + 1
            endif
        end do; end do; end do
        deallocate(hi); allocate(hi(max(1,nhi)))
        nhi = 0
        do l = lb(3), ub(3); do k = lb(2), ub(2); do h = lb(1), ub(1)
            sh = nint(sqrt(real(h*h+k*k+l*l)))
            if( sh >= nint(0.75*nyq_use) .and. sh <= nyq_use .and. dens_acc(h,k,l) > real(COV_DENSITY_FLOOR) )then
                nhi = nhi + 1; hi(nhi) = snr(h,k,l)
            endif
        end do; end do; end do
        if( nhi > 0 )then
            call hpsort(hi(1:nhi))
            floor_noise = real(hi((nhi+1)/2), dp)          ! median outer-shell variance
        else
            floor_noise = 0.d0
        endif
        snr = max(0., snr - real(floor_noise))
        ! zero unobserved cells
        where( dens_acc <= real(COV_DENSITY_FLOOR) ) snr = 0.
        deallocate(var_acc, dens_acc, hi)
    end subroutine estimate_snr_volume

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
        if( ncol < 1 ) THROW_HARD('flex_pca selected no covariance columns; increase lp or ncols')
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

    subroutine build_column_lookup( col_hkl, ncol, lb, ub, col_lookup )
        integer,              intent(in)  :: col_hkl(:,:), ncol, lb(3), ub(3)
        integer, allocatable, intent(out) :: col_lookup(:,:,:)
        integer :: s, h, k, l
        allocate(col_lookup(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), source=0)
        do s = 1, ncol
            h = col_hkl(1,s); k = col_hkl(2,s); l = col_hkl(3,s)
            if( h < lb(1) .or. h > ub(1) .or. k < lb(2) .or. k > ub(2) .or. l < lb(3) .or. l > ub(3) )then
                THROW_HARD('flex_pca selected column outside expanded volume bounds')
            endif
            col_lookup(h,k,l) = s
        end do
    end subroutine build_column_lookup

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
    pure function cov_image_mask_radius( params ) result( r )
        class(parameters), intent(in) :: params
        real :: r
        r = 0.
        if( .not. COV_MASK_IMAGES ) return
        if( params%msk_crop <= 0. .or. params%box_crop <= 0 ) return
        r = COV_MASK_MARGIN * params%msk_crop * real(params%box) / real(params%box_crop)
        r = min(r, 0.5*real(params%box) - COSMSKHALFWIDTH - 1.)
    end function cov_image_mask_radius

    !> Single entry point for the covariance mean.
    subroutine estimate_covariance_mean( params, build, mean_rec, pinds, nptcls )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
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
        integer(timer_int_kind) :: t_phase
        call init_column_reconstructor(params, build, mean_rec)
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
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
        call orientation%kill
        call cleanup_plane(num_fpl)
        call cleanup_rec_buffers(build, fpls)
        if( used < 1 ) THROW_HARD('flex_pca mean estimation found no valid particles')
        ! canonical gridding finalization, identical to reconstruct3D_reference
        call mean_rec%compress_exp
        call mean_rec%sampl_dens_correct
        call mean_rec%ifft
        call mean_rec%div(real(params%box))
        gridcorr_img = prep3D_inv_instrfun4mul([params%box_crop,params%box_crop,params%box_crop], &
            &OSMPL_PAD_FAC*[params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
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
        type(fplane_type) :: mean_fpl
        type(ori) :: orientation
        integer  :: batchlims(2), batchsz, ibatch, i, iptcl, stride, used, nyq, sh, j, nsub
        integer, allocatable :: sub_pinds(:)
        real(dp) :: s_my, s_mm, s, sm
        real(dp), allocatable :: smy_sh(:), smm_sh(:), sprof(:)
        real,     allocatable :: filt(:)
        nyq = max(1, fdim(params%box_crop) - 1)
        allocate(smy_sh(0:nyq), smm_sh(0:nyq), source=0.d0)
        stride = max(1, nptcls / NSAMPLE)
        s_my = 0.d0; s_mm = 0.d0; used = 0
        call mean_rec%expand_exp
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        ! Select the strided sample UP FRONT, not inside the batch loop -- otherwise every particle is
        ! read, normalised, padded, FFT'd and CTF-evaluated before ~(1 - 1/stride) of that is discarded.
        ! Bit-exact: same particles in the same order, hence the same accumulation order.
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
            do i = 1, batchsz
                iptcl = sub_pinds(batchlims(1)+i-1)
                call build%spproj_field%get_ori(iptcl, orientation)
                if( orientation%isstatezero() ) cycle
                call project_fplane_mean(mean_rec, orientation, fpls(i), mean_fpl, apply_ctf_amp=.true.)
                s_my = s_my + real(cov_herm_inner(mean_fpl, fpls(i)), dp)
                s_mm = s_mm + real(cov_herm_inner(mean_fpl, mean_fpl), dp)
                call plane_shell_cross_accum(mean_fpl, fpls(i), nyq, smy_sh, smm_sh)
                used = used + 1
            end do
        end do
        deallocate(sub_pinds)
        call orientation%kill
        call cleanup_plane(mean_fpl)
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

    !> Matched-KB selected-column accumulation over the full expanded sphere.
    subroutine accumulate_covariance_columns( params, build, mean_rec, pinds, nptcls, &
        &col_hkl, col_lookup, ncol, lb, ub, nyq_rec, Bcol_e, Hcol_e, Bcol_o, Hcol_o )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls, ncol, lb(3), ub(3), nyq_rec
        integer,             intent(in)    :: col_hkl(:,:)
        integer,             intent(in)    :: col_lookup(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3))
        complex,             intent(inout) :: Bcol_e(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol)
        complex,             intent(inout) :: Bcol_o(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol)
        real,                intent(inout) :: Hcol_e(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol)
        real,                intent(inout) :: Hcol_o(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol)
        type(fplane_type), allocatable :: fpls(:)
        type(fplane_type) :: mean_fpl, adj_fpl
        type(ori) :: orientation
        type(kbinterpol) :: kbwin
        ! per-particle sample cache (full expanded plane)
        integer, allocatable :: cwin(:,:)      ! (3,ncache) window lower corner
        real,    allocatable :: cwx(:,:), cwy(:,:), cwz(:,:)
        real,    allocatable :: cloc(:,:)   ! continuous sample locations, for the right kernel
        complex, allocatable :: cpl(:)
        real,    allocatable :: cct(:)
        complex, allocatable :: gcol(:)
        real,    allocatable :: hcolv(:)
        integer :: batchlims(2), batchsz, ibatch, i, iptcl, eo, ncache, progress_stride, maxcache
        integer(timer_int_kind) :: t_phase
        call cov_init_right_kernel_width
        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        maxcache = (2*nyq_rec + 3)**2
        allocate(cwin(3,maxcache), cwx(3,maxcache), cwy(3,maxcache), cwz(3,maxcache), &
            &cpl(maxcache), cct(maxcache), gcol(ncol), hcolv(ncol), cloc(3,maxcache))
        call mean_rec%expand_exp
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        progress_stride = max(1, 5 * MAXIMGBATCHSZ)
        t_phase = tic()
        write(logfhandle,'(A)') '>>> FLEX_PCA COLUMN ACCUMULATION'
        call flush(logfhandle)
        if( flex_pca_is_master() )then
            call flex_pca_run_stage(PCA_STAGE_COLS, 'column accumulation')
            call reduce_cols_parts(params, flex_pca_nparts(), Bcol_e, Hcol_e, Bcol_o, Hcol_o, lb, ub, ncol)
            call cleanup_plane(mean_fpl); call cleanup_plane(adj_fpl)
            call cleanup_rec_buffers(build, fpls)
            deallocate(cwin, cwx, cwy, cwz, cpl, cct, gcol, hcolv, cloc)
            return
        endif
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
                eo = build%spproj_field%get_eo(iptcl)
                call project_fplane_mean(mean_rec, orientation, fpls(i), mean_fpl, apply_ctf_amp=.true.)
                call form_adjoint_residual_plane(fpls(i), mean_fpl, adj_fpl, particle_contrast(mean_fpl, fpls(i)))
                call build_sample_cache(kbwin, orientation, adj_fpl, nyq_rec, lb, ub, &
                    &cwin, cwx, cwy, cwz, cpl, cct, ncache, cloc)
                call gather_column_values(col_lookup, ncol, lb, ub, cwin, cwx, cwy, cwz, cpl, cct, &
                    &ncache, gcol, hcolv, cloc)
                if( eo == 0 )then
                    call backproject_columns(col_hkl, ncol, lb, ub, cwin, cwx, cwy, cwz, cpl, cct, &
                        &ncache, gcol, hcolv, Bcol_e, Hcol_e)
                else
                    call backproject_columns(col_hkl, ncol, lb, ub, cwin, cwx, cwy, cwz, cpl, cct, &
                        &ncache, gcol, hcolv, Bcol_o, Hcol_o)
                endif
            end do
            if( batchlims(2) == nptcls .or. mod(batchlims(2), progress_stride) == 0 )then
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA COLUMN PARTICLES: ', batchlims(2), ' / ', nptcls
                call flush(logfhandle)
            endif
        end do
        write(logfhandle,'(A,F10.1)') '>>> FLEX_PCA COLUMN ACCUMULATION SECONDS: ', toc(t_phase)
        call flush(logfhandle)
        call orientation%kill
        call cleanup_plane(mean_fpl)
        call cleanup_plane(adj_fpl)
        call cleanup_rec_buffers(build, fpls)
        deallocate(cwin, cwx, cwy, cwz, cpl, cct, gcol, hcolv, cloc)
        if( flex_pca_is_worker() )then
            call write_cols_part(flex_pca_part_fname('cols', params%part, params%numlen), &
                &Bcol_e, Hcol_e, Bcol_o, Hcol_o, lb, ub, ncol)
        endif
    end subroutine accumulate_covariance_columns

    !> T*r plane and |T|^2 plane for a particle.
    subroutine form_adjoint_residual_plane( fpl, mean_fpl, adj, contrast )
        type(fplane_type), intent(in)    :: fpl, mean_fpl
        type(fplane_type), intent(inout) :: adj
        real, optional,    intent(in)    :: contrast
        integer :: h1, h2, k1, k2
        real    :: a
        a = 1.0; if( present(contrast) ) a = contrast
        h1 = lbound(fpl%cmplx_plane,1); h2 = ubound(fpl%cmplx_plane,1)
        k1 = lbound(fpl%cmplx_plane,2); k2 = ubound(fpl%cmplx_plane,2)
        if( .not. allocated(adj%cmplx_plane) ) allocate(adj%cmplx_plane(h1:h2,k1:k2))
        if( .not. allocated(adj%ctfsq_plane) ) allocate(adj%ctfsq_plane(h1:h2,k1:k2))
        adj%cmplx_plane = conjg(fpl%transfer_plane) * (fpl%cmplx_plane - a*mean_fpl%cmplx_plane)
        adj%ctfsq_plane = fpl%ctfsq_plane
        adj%frlims  = fpl%frlims
        adj%nyq     = fpl%nyq
        adj%shconst = fpl%shconst
    end subroutine form_adjoint_residual_plane

    !> Iterate the full redundant Fourier plane (both k signs, Friedel-reading the stored k<=0 half) and
    !! cache each sample's native location, KB stencil, adjoint residual value and |T|^2 value.
    subroutine build_sample_cache( kbwin, o, adj, nyq_rec, lb, ub, cwin, cwx, cwy, cwz, cpl, cct, ncache, cloc )
        type(kbinterpol),  intent(in)    :: kbwin
        class(ori),        intent(inout) :: o
        type(fplane_type), intent(in)    :: adj
        integer,           intent(in)    :: nyq_rec, lb(3), ub(3)
        integer,           intent(out)   :: cwin(:,:)
        real,              intent(out)   :: cwx(:,:), cwy(:,:), cwz(:,:)
        complex,           intent(out)   :: cpl(:)
        real,              intent(out)   :: cct(:)
        integer,           intent(out)   :: ncache
        ! continuous 3D location of each cached sample, for the independent right kernel
        real, optional,    intent(out)   :: cloc(:,:)
        real    :: rotmat(3,3), loc(3), wx(3), wy(3), wz(3), ctfval
        complex :: plval
        integer :: fpllims(3,2), fpllims_pd(3,2), pf, h, k, hp, kp, nyq_disk
        integer :: h_sq, k_max_h, k_lo, k_hi, win(2,3), plb(2), pub(2), win_lo(3)
        rotmat     = o%get_mat()
        pf         = OSMPL_PAD_FAC
        plb        = lbound(adj%cmplx_plane)
        pub        = ubound(adj%cmplx_plane)
        fpllims_pd = adj%frlims
        fpllims(1,1)= ceil_div(fpllims_pd(1,1),pf); fpllims(1,2)= floor_div(fpllims_pd(1,2),pf)
        fpllims(2,1)= ceil_div(fpllims_pd(2,1),pf); fpllims(2,2)= floor_div(fpllims_pd(2,2),pf)
        nyq_disk = nyq_rec * (nyq_rec + 1)
        ncache   = 0
        do h = fpllims(1,1), fpllims(1,2)
            h_sq = h*h
            if( h_sq > nyq_disk ) cycle
            k_max_h = int(sqrt(real(nyq_disk - h_sq)))
            k_lo    = max(fpllims(2,1), -k_max_h)
            k_hi    = min(fpllims(2,2),  k_max_h)
            hp      = h*pf
            do k = k_lo, k_hi
                kp = k*pf
                if( kp <= 0 )then
                    if( hp < plb(1) .or. hp > pub(1) .or. kp < plb(2) .or. kp > pub(2) ) cycle
                    plval  = adj%cmplx_plane(hp,kp)
                    ctfval = adj%ctfsq_plane(hp,kp)
                else
                    if( -hp < plb(1) .or. -hp > pub(1) .or. -kp < plb(2) .or. -kp > pub(2) ) cycle
                    plval  = conjg(adj%cmplx_plane(-hp,-kp))
                    ctfval = adj%ctfsq_plane(-hp,-kp)
                endif
                if( abs(real(plval)) + abs(aimag(plval)) <= TINY .and. ctfval <= TINY ) cycle
                loc(1) = real(h)*rotmat(1,1) + real(k)*rotmat(2,1)
                loc(2) = real(h)*rotmat(1,2) + real(k)*rotmat(2,2)
                loc(3) = real(h)*rotmat(1,3) + real(k)*rotmat(2,3)
                ! Reject out-of-grid samples BEFORE the Bessel evaluations: cov_kb_weights derives its
                ! window from nint(loc) alone, so the same test can be made on the window corner first.
                win_lo = nint(loc) - COV_KB_IWINSZ
                if( any(win_lo < lb) .or. any(win_lo + 2*COV_KB_IWINSZ > ub) ) cycle
                call cov_kb_weights(kbwin, loc, win, wx, wy, wz)
                ncache = ncache + 1
                cwin(:,ncache) = win(1,:)
                cwx(:,ncache)  = wx
                cwy(:,ncache)  = wy
                cwz(:,ncache)  = wz
                cpl(ncache)    = plval
                cct(ncache)    = ctfval
                if( present(cloc) ) cloc(:,ncache) = loc
            end do
        end do
    end subroutine build_sample_cache

    !>  Right-kernel column values g_is=sum K_R (T*r), h_is=sum K_R |T|^2 by scattering
    !!  each cached sample's KB stencil into the selected-voxel lookup.
    subroutine gather_column_values( col_lookup, ncol, lb, ub, cwin, cwx, cwy, cwz, cpl, cct, ncache, gcol, hcolv, cloc )
        integer, intent(in)  :: ncol, lb(3), ub(3), ncache
        integer, intent(in)  :: col_lookup(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), cwin(:,:)
        real,    intent(in)  :: cwx(:,:), cwy(:,:), cwz(:,:), cct(:)
        complex, intent(in)  :: cpl(:)
        complex, intent(out) :: gcol(ncol)
        real,    intent(out) :: hcolv(ncol)
        real, optional, intent(in) :: cloc(:,:)
        integer :: j, ix, iy, iz, hx, ky, mz, s, iw, h0, k0, m0
        real    :: w, rkw, dx, dy, dz
        gcol  = cmplx(0.,0.)
        hcolv = 0.
        rkw   = COV_RIGHT_KERNEL_W
        if( rkw > 0. .and. present(cloc) )then
            ! independent right kernel: separable triangular window of half-width rkw
            iw = ceiling(rkw)
            do j = 1, ncache
                h0 = nint(cloc(1,j)); k0 = nint(cloc(2,j)); m0 = nint(cloc(3,j))
                do mz = max(lb(3), m0-iw), min(ub(3), m0+iw)
                    dz = abs(real(mz) - cloc(3,j)); if( dz >= rkw ) cycle
                    do ky = max(lb(2), k0-iw), min(ub(2), k0+iw)
                        dy = abs(real(ky) - cloc(2,j)); if( dy >= rkw ) cycle
                        do hx = max(lb(1), h0-iw), min(ub(1), h0+iw)
                            s = col_lookup(hx,ky,mz)
                            if( s == 0 ) cycle
                            dx = abs(real(hx) - cloc(1,j)); if( dx >= rkw ) cycle
                            w = (1.-dx/rkw)*(1.-dy/rkw)*(1.-dz/rkw)
                            gcol(s)  = gcol(s)  + w * cpl(j)
                            hcolv(s) = hcolv(s) + w * cct(j)
                        end do
                    end do
                end do
            end do
            return
        endif
        do j = 1, ncache
            do iz = 1, 3
                mz = cwin(3,j) + iz - 1
                do iy = 1, 3
                    ky = cwin(2,j) + iy - 1
                    do ix = 1, 3
                        hx = cwin(1,j) + ix - 1
                        s  = col_lookup(hx,ky,mz)
                        if( s == 0 ) cycle
                        w = cwx(ix,j)*cwy(iy,j)*cwz(iz,j)
                        gcol(s)  = gcol(s)  + w * cpl(j)
                        hcolv(s) = hcolv(s) + w * cct(j)
                    end do
                end do
            end do
        end do
    end subroutine gather_column_values

    !> Backproject (P*) the per-column numerator and density into the full expanded sphere.
    subroutine backproject_columns( col_hkl, ncol, lb, ub, cwin, cwx, cwy, cwz, cpl, cct, ncache, gcol, hcolv, Bcol, Hcol )
        integer, intent(in)    :: col_hkl(:,:), ncol, lb(3), ub(3), cwin(:,:), ncache
        real,    intent(in)    :: cwx(:,:), cwy(:,:), cwz(:,:), cct(:)
        complex, intent(in)    :: cpl(:), gcol(ncol)
        real,    intent(in)    :: hcolv(ncol)
        complex, intent(inout) :: Bcol(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol)
        real,    intent(inout) :: Hcol(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol)
        integer :: s, is, nlive, j, ix, iy, iz, hx, ky, mz, qh, qk, ql, iqx, iqy, iqz
        integer :: live(ncol)
        real    :: pf2, w, wqs, wyz, hs
        complex :: cg, cnum
        logical :: in_win
        pf2 = real(OSMPL_PAD_FAC*OSMPL_PAD_FAC)
        ! A column only receives anything from this particle if some plane sample's KB window contains it,
        ! which is what gather_column_values scores: gcol(s)==0 .and. hcolv(s)==0 means the whole body adds
        ! exact zeros. Most columns are dead for any one particle. Compact rather than `cycle` so the
        ! static schedule sees a dense iteration space. Only valid for the shared-stencil right kernel:
        ! with COV_RIGHT_KERNEL_W > 0 the gather window differs from the KB window wqs is built on, so a
        ! dead column can still carry a nonzero noise-debias term.
        if( COV_RIGHT_KERNEL_W > 0. )then
            nlive = ncol
            do s = 1, ncol
                live(s) = s
            end do
        else
            nlive = 0
            do s = 1, ncol
                if( gcol(s) /= cmplx(0.,0.) .or. hcolv(s) /= 0. )then
                    nlive = nlive + 1
                    live(nlive) = s
                endif
            end do
            if( nlive == 0 ) return
        endif
        !$omp parallel do default(shared) schedule(static) proc_bind(close) &
        !$omp& private(is,s,j,ix,iy,iz,hx,ky,mz,qh,qk,ql,iqx,iqy,iqz,w,wqs,wyz,hs,cg,cnum,in_win)
        do is = 1, nlive
            s  = live(is)
            cg = pf2 * conjg(gcol(s))
            hs = hcolv(s)
            qh = col_hkl(1,s); qk = col_hkl(2,s); ql = col_hkl(3,s)
            do j = 1, ncache
                ! right-kernel weight of this sample to q_s (nonzero only if q_s in window)
                in_win = qh >= cwin(1,j) .and. qh <= cwin(1,j)+2 .and. &
                    &    qk >= cwin(2,j) .and. qk <= cwin(2,j)+2 .and. &
                    &    ql >= cwin(3,j) .and. ql <= cwin(3,j)+2
                wqs = 0.
                if( in_win .and. COV_COLUMN_NOISE_DEBIAS )then
                    iqx = qh - cwin(1,j) + 1
                    iqy = qk - cwin(2,j) + 1
                    iqz = ql - cwin(3,j) + 1
                    wqs = pf2 * cwx(iqx,j)*cwy(iqy,j)*cwz(iqz,j) * cct(j)   ! noise-bias scale
                endif
                cnum = cg * cpl(j)
                do iz = 1, 3
                    mz = cwin(3,j) + iz - 1
                    do iy = 1, 3
                        ky  = cwin(2,j) + iy - 1
                        wyz = cwy(iy,j)*cwz(iz,j)
                        do ix = 1, 3
                            hx = cwin(1,j) + ix - 1
                            w  = cwx(ix,j)*wyz
                            Bcol(hx,ky,mz,s) = Bcol(hx,ky,mz,s) + cnum*w - wqs*w
                            Hcol(hx,ky,mz,s) = Hcol(hx,ky,mz,s) + hs*w*cct(j)
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine backproject_columns

    !> KB separable stencil, replicating the one simple_flex_reconstructor_latent_ops uses, so the column
    !! gather and the backprojection see identical interpolation weights.
    pure subroutine cov_kb_weights( kbwin, loc, win, wx, wy, wz )
        type(kbinterpol), intent(in)  :: kbwin
        real,             intent(in)  :: loc(3)
        integer,          intent(out) :: win(2,3)
        real,             intent(out) :: wx(3), wy(3), wz(3)
        integer :: i, iwinsz, win_lo(3)
        real    :: base(3), ww3(3), sx, sy, sz, eps_norm
        iwinsz   = ceiling(KBWINSZ - 0.5)
        win(1,:) = nint(loc)
        win(2,:) = win(1,:) + iwinsz
        win(1,:) = win(1,:) - iwinsz
        win_lo   = win(1,:)
        base     = real(win_lo) - loc
        do i = 1, 3
            ww3   = kbwin%apod(base + real(i-1))
            wx(i) = ww3(1); wy(i) = ww3(2); wz(i) = ww3(3)
        end do
        sx = sum(wx); sy = sum(wy); sz = sum(wz)
        eps_norm = epsilon(1.0)
        if( abs(sx) > eps_norm )then; wx = wx/sx; else; wx = 1.0/3.0; endif
        if( abs(sy) > eps_norm )then; wy = wy/sy; else; wy = 1.0/3.0; endif
        if( abs(sz) > eps_norm )then; wz = wz/sz; else; wz = 1.0/3.0; endif
    end subroutine cov_kb_weights

    !> Generalized-FSC column regularization (supplement S.C, eqs S.11-S.13).
    subroutine regularize_and_merge_columns( Bcol_e, Hcol_e, Bcol_o, Hcol_o, ncol, lb_exp, colvol, col_fsc )
        complex, intent(in)  :: Bcol_e(:,:,:,:), Bcol_o(:,:,:,:)
        real,    intent(in)  :: Hcol_e(:,:,:,:), Hcol_o(:,:,:,:)
        integer, intent(in)  :: ncol
        !> lower bounds of the expanded reconstructor grid, needed for the frequency origin
        integer, intent(in)  :: lb_exp(3)
        complex, intent(out) :: colvol(:,:,:,:)
        real,    intent(out) :: col_fsc(ncol)
        integer, parameter :: NFSC_ITERS = 20
        real(dp), allocatable :: Rsh(:), Psh(:), fsc(:), num(:), den_e(:), den_o(:), top(:), bot(:), hbar_sum(:)
        real(dp), allocatable :: nvox_sh(:), pri_accum(:), shr_accum(:), den_accum(:)
        real(dp), allocatable :: fsc_accum(:)
        integer,  allocatable :: shell(:,:,:)
        real(dp) :: hbar, de, do_, w2, snr, fscbar, fscden, fscnum
        complex(dp) :: ce, co
        integer  :: s, i1, i2, i3, n1, n2, n3, h, k, l, c1, c2, c3, sh, nsh, it
        n1 = size(Bcol_e,1); n2 = size(Bcol_e,2); n3 = size(Bcol_e,3)
        ! Frequency origin of the EXPANDED reconstructor grid.
        c1 = 1 - lb_exp(1); c2 = 1 - lb_exp(2); c3 = 1 - lb_exp(3)
        nsh = nint(sqrt(real((n1/2)**2 + (n2/2)**2 + (n3/2)**2))) + 1
        allocate(shell(n1,n2,n3))
        allocate(Rsh(0:nsh), Psh(0:nsh), fsc(0:nsh), num(0:nsh), den_e(0:nsh), den_o(0:nsh), &
            &top(0:nsh), bot(0:nsh), hbar_sum(0:nsh), fsc_accum(0:nsh), &
            &nvox_sh(0:nsh), pri_accum(0:nsh), shr_accum(0:nsh), den_accum(0:nsh))
        fsc_accum = 0.d0; pri_accum = 0.d0; shr_accum = 0.d0; den_accum = 0.d0
        do i3 = 1, n3
            l = i3 - c3
            do i2 = 1, n2
                k = i2 - c2
                do i1 = 1, n1
                    h = i1 - c1
                    shell(i1,i2,i3) = min(nsh, nint(sqrt(real(h*h + k*k + l*l))))
                end do
            end do
        end do
        ! ---- Regularizer initialization. The supplement's v v^T prior cannot be formed here: it assumes
        ! the whitened convention sigma^2 = 1, whereas SIMPLE's non-unitary FFT/gridding convention puts
        ! sigma^2 at ~1e-6 and the mean map carries yet another amplitude convention. The prior is instead
        ! seeded at the density floor and driven entirely by the generalized-FSC iterations below.
        do s = 1, ncol
            hbar_sum = 0.d0; num = 0.d0; nvox_sh = 0.d0
            do i3 = 1, n3; do i2 = 1, n2; do i1 = 1, n1
                sh = shell(i1,i2,i3)
                hbar_sum(sh) = hbar_sum(sh) + 0.5d0*(real(Hcol_e(i1,i2,i3,s),dp)+real(Hcol_o(i1,i2,i3,s),dp))
                num(sh) = num(sh) + 1.d0
                nvox_sh(sh) = nvox_sh(sh) + 1.d0
            end do; end do; end do
            ! Initial ridge, at the density floor so the first iteration is essentially unregularized
            do sh = 0, nsh
                Rsh(sh) = COV_DENSITY_FLOOR
                Psh(sh) = 1.d0 / Rsh(sh)
            end do
            ! generalized-FSC iterations (no data pass needed)
            do it = 1, NFSC_ITERS
                num = 0.d0; den_e = 0.d0; den_o = 0.d0; top = 0.d0; bot = 0.d0
                do i3 = 1, n3; do i2 = 1, n2; do i1 = 1, n1
                    sh  = shell(i1,i2,i3)
                    de  = real(Hcol_e(i1,i2,i3,s),dp) + Rsh(sh)
                    do_ = real(Hcol_o(i1,i2,i3,s),dp) + Rsh(sh)
                    if( de <= COV_DENSITY_FLOOR .or. do_ <= COV_DENSITY_FLOOR ) cycle
                    ce = cmplx(Bcol_e(i1,i2,i3,s),kind=dp) / de
                    co = cmplx(Bcol_o(i1,i2,i3,s),kind=dp) / do_
                    num(sh)   = num(sh)   + real(ce*conjg(co), dp)
                    den_e(sh) = den_e(sh) + real(ce*conjg(ce), dp)
                    den_o(sh) = den_o(sh) + real(co*conjg(co), dp)
                    hbar = 0.5d0*(real(Hcol_e(i1,i2,i3,s),dp)+real(Hcol_o(i1,i2,i3,s),dp))
                    w2   = 1.d0 / (hbar + Rsh(sh))**2                 ! (H+R)^-2
                    top(sh) = top(sh) + w2*hbar                       ! sum (H+R)^-2 H
                    bot(sh) = bot(sh) + w2*hbar*hbar                  ! sum (H+R)^-2 H^2
                end do; end do; end do
                do sh = 0, nsh
                    ! Guard at zero, NOT at DTINY: these are sums of |B/(H+R)|^2 in the pipeline's
                    ! non-unitary Fourier convention (sig2_eff ~ 1e-6, not 1), so a DTINY floor zeroes the
                    ! FSC at the lowest shells, and via R = 1/P with P proportional to FSC/(1-FSC) that
                    ! deletes the conformational band from every column.
                    if( den_e(sh) > 0.d0 .and. den_o(sh) > 0.d0 )then
                        fsc(sh) = num(sh) / sqrt(den_e(sh)*den_o(sh))
                    else
                        fsc(sh) = 0.d0
                    endif
                    snr = min(1.d0-1.d-3, max(1.d-3, fsc(sh))) / (1.d0 - min(1.d0-1.d-3, max(1.d-3, fsc(sh))))
                    if( bot(sh) > DTINY )then
                        Psh(sh) = max(COV_DENSITY_FLOOR, snr * top(sh) / bot(sh))
                        Rsh(sh) = 1.d0 / Psh(sh)
                    endif
                end do
            end do
            ! per-shell profiles for the log
            do sh = 0, nsh
                fscbar = max(0.d0, min(0.999d0, fsc(sh)))
                fsc_accum(sh) = fsc_accum(sh) + fscbar     ! raw FSC for the profile log
                ! Converged prior and the Wiener shrinkage it produces at the shell's mean density.
                pri_accum(sh) = pri_accum(sh) + Psh(sh)
                den_accum(sh) = den_accum(sh) + sqrt(max(den_e(sh)*den_o(sh),0.d0))
                hbar = 0.d0
                if( nvox_sh(sh) > 0.d0 ) hbar = hbar_sum(sh)/nvox_sh(sh)
                shr_accum(sh) = shr_accum(sh) + hbar/(hbar + Rsh(sh))
            end do
            ! merged regularized column (Wiener/MAP denominator H + 1/P) and mean column FSC diagnostic
            fscnum = 0.d0; fscden = 0.d0
            do i3 = 1, n3; do i2 = 1, n2; do i1 = 1, n1
                sh = shell(i1,i2,i3)
                ! the halves are SUMMED, not averaged, so the ridge enters once: Rsh, not 2*Rsh
                colvol(i1,i2,i3,s) = (Bcol_e(i1,i2,i3,s) + Bcol_o(i1,i2,i3,s)) &
                    &/ (Hcol_e(i1,i2,i3,s) + Hcol_o(i1,i2,i3,s) + real(Rsh(sh)))
            end do; end do; end do
            do sh = 0, nsh
                if( den_e(sh) > DTINY .and. den_o(sh) > DTINY )then
                    fscbar = num(sh)/sqrt(den_e(sh)*den_o(sh))
                    fscnum = fscnum + fscbar*(den_e(sh)+den_o(sh))
                    fscden = fscden + (den_e(sh)+den_o(sh))
                endif
            end do
            if( fscden > DTINY )then
                col_fsc(s) = real(fscnum/fscden)
            else
                col_fsc(s) = 0.
            endif
        end do
        ! mean per-shell column FSC profile (verifies high-freq suppression)
        write(logfhandle,'(A)') '>>> FLEX_PCA mean column FSC / prior / Wiener shrinkage per shell:'
        write(logfhandle,'(A)') '>>>   (shrink = <H>/(<H>+R) is the fraction of the column that survives the ridge)'
        do sh = 0, min(nsh, 32)
            write(logfhandle,'(A,I3,A,I7,A,F7.4,A,ES10.2,A,ES11.3,A,F7.4)') '>>>   sh=', sh, &
                &'  nvox=', nint(nvox_sh(sh)), &
                &'  fsc=', real(fsc_accum(sh)/real(ncol,dp)), &
                &'  |col|=', real(sqrt(max(den_accum(sh)/real(ncol,dp),0.d0))), &
                &'  prior=', real(pri_accum(sh)/real(ncol,dp)), &
                &'  shrink=', real(shr_accum(sh)/real(ncol,dp))
        end do
        deallocate(shell, Rsh, Psh, fsc, num, den_e, den_o, top, bot, hbar_sum, fsc_accum, &
            &nvox_sh, pri_accum, shr_accum, den_accum)
    end subroutine regularize_and_merge_columns

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
        gridcorr_img = prep3D_inv_instrfun4mul([params%box_crop,params%box_crop,params%box_crop], &
            &OSMPL_PAD_FAC*[params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
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
        call work%reset
        call work%reset_exp
        work%cmat_exp = vherm
        call work%compress_exp
        ! Band-limit the covariance column to the signal band FIRST, in Fourier space, before any real-space
        ! operation, so out-of-band shells never enter the masking or the Gram products downstream.
        if( params%lp > 2.0*params%smpd_crop + TINY ) call work%bp(0., params%lp)
        call work%ifft
        call work%div(real(params%box))
        if( params%msk_crop > TINY ) call work%mask3D_soft(params%msk_crop, backgr=0.)
        if( work%is_ft() ) call work%ifft
        call work%get_rmat_ptr(rmat)
        energy = sum(rmat*rmat)
    end subroutine realize_hermitian_volume

    !> Orthonormalize the real column representatives into the column subspace Utilde by Gram
    !! eigendecomposition, keeping every direction above a relative energy floor.
    subroutine orthonormalize_representatives( params, build, realvols, nreal, utilde, utilde_real, d_tilde, svals )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(image),         intent(inout) :: realvols(:)
        integer,             intent(in)    :: nreal
        type(reconstructor), allocatable, intent(out) :: utilde(:)
        type(image),         allocatable, intent(out) :: utilde_real(:)
        integer,             intent(out)   :: d_tilde
        !> squared singular values of the representative set.
        real(dp), allocatable, optional, intent(out) :: svals(:)
        real(dp), allocatable :: gram(:,:), evec(:,:), eval(:)
        real, pointer :: rmat_i(:,:,:), rmat_j(:,:,:)
        integer :: i, q, nrot, keep, d_budget
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
        ! SIMPLE_COV_DTILDE pins the subspace dimension, so accumulation or solver changes can be compared
        ! at a FIXED d instead of confounding the two
        call cov_env_int('SIMPLE_COV_DTILDE', d_budget)
        d_tilde  = max(1, min(keep, COV_MAX_DTILDE, d_budget))
        if( l_packed )then
            accum_model = 'packed+CG'
        else
            accum_model = 'dense'
        endif
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA d_tilde=',d_tilde, &
            &'  (above energy floor=',keep,', memory cap=',d_budget,', rank cap=',COV_MAX_DTILDE,')'
        write(logfhandle,'(A,A,A,F8.3,A,F6.3,A)') '>>> FLEX_PCA reduced-solve accumulator model: ', &
            &trim(accum_model),', ',cov_accum_bytes(d_tilde, l_packed)/1.d9, &
            &' GB at this d_tilde (budget ',COV_ATHR_BUDGET/1.d9,' GB)'
        if( d_tilde == d_budget .and. keep > d_budget )then
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA NOTE: the column subspace is limited by the &
                &reduced-solve memory budget, not by the data; ',keep,' directions cleared the energy floor.'
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


    !> Packed symmetric (svec) index of the unordered pair {i,j}, upper-triangle column order.
    pure integer function svec_idx( i, j, d ) result( k )
        integer, intent(in) :: i, j, d
        integer :: a, b
        a = min(i,j); b = max(i,j)
        k = (b-1)*b/2 + a
    end function svec_idx

    !> sum_i G_i[a,b] G_i[c,d] read back out of the PACKED accumulator Ms = sum_i svec(G_i) svec(G_i)^T.
    pure real(dp) function packed_T( Ms, np, a, b, c, dd_, d ) result( t )
        integer,  intent(in) :: np, a, b, c, dd_, d
        real(dp), intent(in) :: Ms(np,np)
        real(dp) :: w1, w2
        w1 = merge(1.d0, sqrt(2.d0), a == b)
        w2 = merge(1.d0, sqrt(2.d0), c == dd_)
        t  = Ms(svec_idx(a,b,d), svec_idx(c,dd_,d)) / (w1*w2)
    end function packed_T

    !> Apply the packed operator A_s = P A P^T (A = sum_i G_i (x) G_i, P the svec projector) WITHOUT EVER
    !! FORMING A_s. Ms here is the RAW accumulator, never the output of rearrange_packed_selfsum.
    !! The production solve rearranges once and uses dsymv instead; this remains the definition of
    !! the operator and the reference the rearrangement is tested against.
    subroutine apply_packed_A( Ms, np, d, x, y )
        integer,  intent(in)  :: np, d
        real(dp), intent(in)  :: Ms(np,np), x(np)
        real(dp), intent(out) :: y(np)
        integer  :: al, be, ga, de, k, l
        real(dp) :: acc, term, wk, wl
        !$omp parallel do default(shared) private(al,be,ga,de,k,l,acc,term,wk,wl) schedule(static)
        do k = 1, np
            ! recover (al,be) with al<=be from the packed index k
            be = int((sqrt(8.d0*real(k,dp)-7.d0)+1.d0)/2.d0)
            if( svec_idx(1,be,d) > k ) be = be - 1
            al = k - (be-1)*be/2
            wk = merge(1.d0, sqrt(2.d0), al == be)
            acc = 0.d0
            do de = 1, d
                do ga = 1, de
                    l  = svec_idx(ga,de,d)
                    wl = merge(1.d0, sqrt(2.d0), ga == de)
                    if( al == be .and. ga == de )then
                        term = packed_T(Ms, np, al, ga, al, ga, d)
                    else if( al == be )then
                        term = sqrt(2.d0)*packed_T(Ms, np, al, ga, al, de, d)
                    else if( ga == de )then
                        term = sqrt(2.d0)*packed_T(Ms, np, al, ga, be, ga, d)
                    else
                        term = packed_T(Ms, np, al, ga, be, de, d) + packed_T(Ms, np, al, de, be, ga, d)
                    endif
                    acc = acc + term*x(l)
                end do
            end do
            y(k) = acc
        end do
        !$omp end parallel do
    end subroutine apply_packed_A

    !> Conjugate gradients on the packed normal equations A_s s = rhs, with Jacobi preconditioning.
    !! As must ALREADY have been through rearrange_packed_selfsum: this applies it with a single
    !! dsymv per iteration, which is the entire point of rearranging. Passing the raw accumulator
    !! here solves the wrong system silently.
    subroutine cg_solve_packed( As, np, rhs, ridge, maxit, tol, x, iters, relres )
        integer,  intent(in)    :: np, maxit
        real(dp), intent(in)    :: As(np,np), rhs(np), ridge, tol
        real(dp), intent(out)   :: x(np)
        integer,  intent(out)   :: iters
        real(dp), intent(out)   :: relres
        real(dp), allocatable :: r(:), p(:), Ap(:), z(:), diag(:)
        real(dp) :: rz, rz_new, alpha, beta_cg, pAp, nrm0
        integer  :: it, k
        allocate(r(np), p(np), Ap(np), z(np), diag(np))
        ! Jacobi preconditioner: after the rearrangement the operator diagonal IS the matrix diagonal
        do k = 1, np
            diag(k) = max(As(k,k) + ridge, DTINY)
        end do
        x = 0.d0
        r = rhs
        nrm0 = sqrt(sum(r*r))
        if( nrm0 <= DTINY )then
            iters = 0; relres = 0.d0
            deallocate(r,p,Ap,z,diag); return
        endif
        z  = r/diag
        p  = z
        rz = sum(r*z)
        iters = maxit; relres = 1.d0
        do it = 1, maxit
            call dsymv('U', np, 1.d0, As, np, p, 1, 0.d0, Ap, 1)
            Ap  = Ap + ridge*p
            pAp = sum(p*Ap)
            if( pAp <= 0.d0 )then          ! loss of positive definiteness: stop, do not diverge
                iters = it; exit
            endif
            alpha = rz/pAp
            x = x + alpha*p
            r = r - alpha*Ap
            relres = sqrt(sum(r*r))/nrm0
            if( relres < tol )then
                iters = it; exit
            endif
            z = r/diag
            rz_new = sum(r*z)
            beta_cg = rz_new/rz
            rz = rz_new
            p = z + beta_cg*p
        end do
        deallocate(r,p,Ap,z,diag)
    end subroutine cg_solve_packed

    !> svec weight of packed index p: 1 on a diagonal pair (a,a), sqrt(2) otherwise.
    pure real(dp) function packed_weight( p ) result( w )
        integer, intent(in) :: p
        integer :: b
        b = int((sqrt(8.d0*real(p,dp)-7.d0)+1.d0)/2.d0)
        if( (b-1)*b/2 + 1 > p ) b = b - 1
        w = merge(1.d0, sqrt(2.d0), p - (b-1)*b/2 == b)
    end function packed_weight

    !> IN-PLACE packed counterpart of unrearrange_kron_selfsum: turns the accumulator
    !! Ms = sum_i svec(G_i) svec(G_i)^T into the operator A_s = P (sum_i G_i (x) G_i) P^T, so CG can
    !! apply it with ONE dsymv instead of apply_packed_A's scalar gather per term. Measured at
    !! d=177: the gather runs ~20x below the bandwidth this matrix streams at, and the one-time
    !! rearrangement below costs less than a single gathered apply -- it repays inside the first
    !! CG iteration, and there are typically ~100.
    !!
    !! Every position of Ms is an unordered pair-of-pairs {P,Q} of four indices, and the THREE
    !! pairings of four indices form a closed orbit. Orbits are disjoint and every orbit map is
    !! invertible, so no second npk x npk buffer is needed. Writing T for the accumulator value with
    !! the svec weights divided out, the orbit maps by index multiplicity are
    !!
    !!   4 distinct   3 positions   A1 = T2+T3,          A2 = T1+T3,  A3 = T1+T2
    !!   {a,a,b,c}    2 positions   A(P) = sqrt2*T(Q),   A(Q) = T(P)+T(Q)
    !!   {a,a,b,b}    2 positions   A(P) = T(Q),         A(Q) = T(P)+T(Q)
    !!   {a,a,a,b}    1 position    A = sqrt2*T
    !!   {a,a,a,a}    1 position    A = T
    !!
    !! The last four are exactly the positions whose row pair or column pair is diagonal, i.e. where
    !! two of the three pairings collapse onto one and the generic sum would double-count. This is
    !! verified against apply_packed_A entry by entry in test_flex_pca_packed_solve.
    !!
    !! NOTE the rearrangement is NOT a permutation -- it changes the spectrum (Ms has rank <= nptcls,
    !! A_s is generically full rank), so it must run exactly once, and apply_packed_A must never be
    !! handed the result.
    subroutine rearrange_packed_selfsum( Ms, np, d )
        integer,  intent(in)    :: np, d
        real(dp), intent(inout) :: Ms(np,np)
        integer  :: i, j, k, l, a, b, c, r1, c1, r2, c2, r3, c3
        real(dp) :: T1, T2, T3, A1, A2, A3
        !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
        !$omp& private(i,j,k,l,a,b,c,r1,c1,r2,c2,r3,c3,T1,T2,T3,A1,A2,A3)
        do i = 1, d
            do j = i, d
                do k = j, d
                    do l = k, d
                        if( i == l )then                                    ! {a,a,a,a}
                            r1 = svec_idx(i,i,d)
                            Ms(r1,r1) = packed_t_at(Ms,np,r1,r1)
                        else if( i == k )then                               ! {a,a,a,b}, a=i b=l
                            r1 = svec_idx(i,i,d); c1 = svec_idx(i,l,d)
                            A1 = sqrt(2.d0)*packed_t_at(Ms,np,r1,c1)
                            Ms(r1,c1) = A1; Ms(c1,r1) = A1
                        else if( j == l )then                               ! {a,a,a,b}, a=j b=i
                            r1 = svec_idx(j,j,d); c1 = svec_idx(j,i,d)
                            A1 = sqrt(2.d0)*packed_t_at(Ms,np,r1,c1)
                            Ms(r1,c1) = A1; Ms(c1,r1) = A1
                        else if( i == j .and. k == l )then                  ! {a,a,b,b}
                            r1 = svec_idx(i,i,d); c1 = svec_idx(k,k,d)
                            r2 = svec_idx(i,k,d)
                            T1 = packed_t_at(Ms,np,r1,c1); T2 = packed_t_at(Ms,np,r2,r2)
                            A1 = T2; A2 = T1 + T2
                            Ms(r1,c1) = A1; Ms(c1,r1) = A1
                            Ms(r2,r2) = A2
                        else if( i == j .or. j == k .or. k == l )then       ! {a,a,b,c}
                            if( i == j )then
                                a = i; b = k; c = l
                            else if( j == k )then
                                a = j; b = i; c = l
                            else
                                a = k; b = i; c = j
                            endif
                            r1 = svec_idx(a,a,d); c1 = svec_idx(b,c,d)
                            r2 = svec_idx(a,b,d); c2 = svec_idx(a,c,d)
                            T1 = packed_t_at(Ms,np,r1,c1); T2 = packed_t_at(Ms,np,r2,c2)
                            A1 = sqrt(2.d0)*T2; A2 = T1 + T2
                            Ms(r1,c1) = A1; Ms(c1,r1) = A1
                            Ms(r2,c2) = A2; Ms(c2,r2) = A2
                        else                                                ! four distinct
                            r1 = svec_idx(i,j,d); c1 = svec_idx(k,l,d)
                            r2 = svec_idx(i,k,d); c2 = svec_idx(j,l,d)
                            r3 = svec_idx(i,l,d); c3 = svec_idx(j,k,d)
                            T1 = packed_t_at(Ms,np,r1,c1)
                            T2 = packed_t_at(Ms,np,r2,c2)
                            T3 = packed_t_at(Ms,np,r3,c3)
                            A1 = T2 + T3; A2 = T1 + T3; A3 = T1 + T2
                            Ms(r1,c1) = A1; Ms(c1,r1) = A1
                            Ms(r2,c2) = A2; Ms(c2,r2) = A2
                            Ms(r3,c3) = A3; Ms(c3,r3) = A3
                        endif
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine rearrange_packed_selfsum

    !> T at a packed position: the stored entry with the svec weights divided back out.
    pure real(dp) function packed_t_at( Ms, np, r, c ) result( t )
        integer,  intent(in) :: np, r, c
        real(dp), intent(in) :: Ms(np,np)
        t = Ms(r,c)/(packed_weight(r)*packed_weight(c))
    end function packed_t_at

    subroutine unrearrange_kron_selfsum( A, d )
        integer,  intent(in)    :: d
        real(dp), intent(inout) :: A(d*d,d*d)
        integer  :: alpha, beta, gam, delta, r1, r2, c1, c2
        real(dp) :: tmp
        !$omp parallel do collapse(2) default(shared) schedule(static) &
        !$omp& private(alpha,beta,gam,delta,r1,r2,c1,c2,tmp)
        do delta = 1, d
            do alpha = 1, d
                do beta = 1, d
                    do gam = beta+1, d
                        r1 = (beta-1)*d + alpha; c1 = (delta-1)*d + gam
                        r2 = (gam -1)*d + alpha; c2 = (delta-1)*d + beta
                        tmp     = A(r1,c1)
                        A(r1,c1) = A(r2,c2)
                        A(r2,c2) = tmp
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine unrearrange_kron_selfsum

    subroutine reduced_covariance_solve( params, build, mean_rec, utilde, d_tilde, pinds, nptcls, &
        &neigs_req, vred, gamma_out, ncomp_out, sig2_out )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), intent(inout) :: utilde(d_tilde)
        integer,             intent(in)    :: d_tilde, pinds(:), nptcls, neigs_req
        real(dp), allocatable, intent(out) :: vred(:,:), gamma_out(:)
        integer,               intent(out) :: ncomp_out
        real(dp),              intent(out) :: sig2_out
        type(fplane_type), allocatable :: fpls(:)
        type(fplane_type), allocatable :: basis_fpls(:,:), mean_fpl(:), resid_fpl(:)
        type(ori), allocatable :: orientations(:)
        real(dp), allocatable :: A(:,:), Sbb(:,:), SG(:,:), Dmat(:,:), Sig(:,:), V(:,:), gam(:), xsol(:), bvec(:)
        real(dp), allocatable :: Mspk(:,:), spk(:), rhspk(:)
        integer :: cgsolve, npk, cgit
        real(dp) :: cgres
        real(dp), allocatable :: Vg(:,:), Sbb_thr(:,:,:), SG_thr(:,:,:), Gi_thr(:,:,:)
        ! TEMPORARY sub-stage instrumentation for the solve
        real(dp), allocatable :: t_proj_thr(:), t_gram_thr(:), t_diag_thr(:)
        real(dp) :: tw0, t_syrk
        ! compacted per-particle basis samples for the BLAS-3 Gram, real interleaved (re,im,re,im,...)
        real(dp),    allocatable :: gbuf(:,:,:)
        integer,     allocatable :: slist(:,:,:)
        complex(dp) :: gacc
        integer     :: nsamp, nsamp_max, j
        logical     :: l_gram_fast
        real(dp), allocatable :: solve_scal(:)
        real(dp), allocatable :: hfpw_thr(:), hfcnt_thr(:)
        real(dp), allocatable :: cvpw_thr(:), cvcnt_thr(:)   ! log only (see cov_herm_selfpower)
        real(dp), allocatable :: babs_thr(:), bre_thr(:), gdi_thr(:)   ! log only: |b|^2 vs Re(b)^2 vs G_qq
        real(dp), allocatable :: pwsh_thr(:,:), cntsh_thr(:,:), pwsh(:), cntsh(:)   ! per-shell noise profile
        complex(dp), allocatable :: bc_thr(:,:)
        integer,  allocatable :: nvalid_thr(:)
        integer :: batchlims(2), batchsz, ibatch, i, iptcl, q, r, dd, a1, a2, nrot, keep, ncomp
        integer :: alpha, beta, gam1, delta, nvalid, ithr, nthr, nyq_rec, info
        real(dp) :: ridge, trc, gmax, sig2_eff, tr_bb, tr_g, tr_match, pw, cnt
        logical  :: l_trmatch
        integer(timer_int_kind) :: t_phase
        dd   = d_tilde*d_tilde
        nthr = nthr_glob
        ! PACKED accumulation and matrix-free CG, on by default; SIMPLE_COV_CGSOLVE=0 restores the dense
        ! -- through cov_packed_cgsolve, the SAME predicate that sized d_tilde in
        ! orthonormalize_representatives, so the accumulator below is the one that budget paid for.
        ! d^4 accumulator and direct solve. Packing exploits the symmetry of every G_i, so the rank-k
        ! update runs on npk = d(d+1)/2 rather than d^2 rows -- 4x less memory, which matters more than
        ! the speed: the dense accumulator is 8*d^4 bytes, 2.1 GB at ncols=64 and 34 GB at ncols=128.
        ! CG at COV_CG_TOL shifts eigenvalues by ~1e-5 relative; convergence is reported, not assumed.
        cgsolve = merge(1, 0, cov_packed_cgsolve())
        npk = d_tilde*(d_tilde+1)/2
        if( cgsolve /= 0 )then
            allocate(Mspk(npk,npk), source=0.d0)
            write(logfhandle,'(A,I0,A,F8.3,A,F8.3,A)') '>>> FLEX_PCA PACKED+CG solve: npk=',npk, &
                &'  accumulator ', 8.d0*real(npk,dp)**2/1.d9, ' GB (dense would be ', &
                &8.d0*real(dd,dp)**2/1.d9, ' GB)'
            call flush(logfhandle)
        else
            allocate(A(dd,dd), source=0.d0)
            write(logfhandle,'(A,I0,A,F8.3,A)') '>>> FLEX_PCA DENSE solve: dd=',dd, &
                &'  accumulator ', cov_accum_bytes(d_tilde, .false.)/1.d9, ' GB'
            call flush(logfhandle)
        endif
        allocate(Sbb(d_tilde,d_tilde), source=0.d0)
        allocate(SG(d_tilde,d_tilde), source=0.d0)
        ! per-batch staging for the rank-k update: column i holds vec(G_i)
        if( cgsolve /= 0 )then
            allocate(Vg(npk,MAXIMGBATCHSZ), source=0.d0)
        else
            allocate(Vg(dd,MAXIMGBATCHSZ), source=0.d0)
        endif
        allocate(Sbb_thr(d_tilde,d_tilde,nthr), source=0.d0)
        allocate(SG_thr(d_tilde,d_tilde,nthr), source=0.d0)
        allocate(Gi_thr(d_tilde,d_tilde,nthr), bc_thr(d_tilde,nthr))
        allocate(basis_fpls(d_tilde,nthr), mean_fpl(nthr), resid_fpl(nthr))
        allocate(orientations(MAXIMGBATCHSZ), nvalid_thr(nthr))
        allocate(t_proj_thr(nthr), t_gram_thr(nthr), t_diag_thr(nthr), source=0.d0)
        t_syrk = 0.d0
        allocate(hfpw_thr(nthr), hfcnt_thr(nthr), source=0.d0)
        allocate(cvpw_thr(nthr), cvcnt_thr(nthr), source=0.d0)   ! log only: noise in cov_herm_inner's convention
        allocate(babs_thr(nthr), bre_thr(nthr), gdi_thr(nthr), source=0.d0)   ! log only
        nyq_rec = utilde(1)%get_lfny(1)
        allocate(pwsh_thr(0:nyq_rec,nthr), cntsh_thr(0:nyq_rec,nthr), source=0.d0)
        ! Upper bound on the Gram sample count: cov_herm_inner steps OSMPL_PAD_FAC in both axes over
        ! |h|,|k| <= nyq_eff, and nyq_eff never exceeds the padded plane Nyquist. cov_herm_sample_list
        ! reports -1 if it would overrun, and the pair loops then fall back to cov_herm_inner.
        nsamp_max = (2*OSMPL_PAD_FAC*nyq_rec + 4) * (OSMPL_PAD_FAC*nyq_rec + 4) / (OSMPL_PAD_FAC*OSMPL_PAD_FAC)
        allocate(slist(2, nsamp_max, nthr))
        allocate(gbuf(2*nsamp_max, d_tilde, nthr))
        write(logfhandle,'(A,F8.1,A,I0)') '>>> FLEX_PCA Gram sample buffer: ', &
            &16.d0*real(nsamp_max,dp)*real(d_tilde,dp)*real(nthr,dp)/1.d6,' MB  max_samples=',nsamp_max
        call flush(logfhandle)
        allocate(pwsh(0:nyq_rec), cntsh(0:nyq_rec), source=0.d0)
        nvalid_thr = 0
        call mean_rec%expand_exp
        do q = 1, d_tilde
            call utilde(q)%expand_exp
        end do
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        t_phase = tic()
        write(logfhandle,'(A,I0)') '>>> FLEX_PCA REDUCED COVARIANCE SOLVE d_tilde=',d_tilde
        call flush(logfhandle)
        if( flex_pca_is_master() )then
            call flex_pca_run_stage(PCA_STAGE_SOLVE, 'reduced solve')
            allocate(solve_scal(8), source=0.d0)
            if( cgsolve /= 0 )then
                call reduce_solve_parts(params, flex_pca_nparts(), Mspk, npk, Sbb, SG, d_tilde, &
                    &pwsh, cntsh, nyq_rec, solve_scal, nvalid)
            else
                call reduce_solve_parts(params, flex_pca_nparts(), A, dd, Sbb, SG, d_tilde, &
                    &pwsh, cntsh, nyq_rec, solve_scal, nvalid)
            endif
            hfpw_thr(1)=solve_scal(1); hfcnt_thr(1)=solve_scal(2)
            cvpw_thr(1)=solve_scal(3); cvcnt_thr(1)=solve_scal(4)
            babs_thr(1)=solve_scal(5); bre_thr(1)=solve_scal(6); gdi_thr(1)=solve_scal(7)
            if( nthr > 1 )then
                hfpw_thr(2:)=0.d0; hfcnt_thr(2:)=0.d0; cvpw_thr(2:)=0.d0; cvcnt_thr(2:)=0.d0
                babs_thr(2:)=0.d0; bre_thr(2:)=0.d0;  gdi_thr(2:)=0.d0
            endif
            nvalid_thr    = 0
            nvalid_thr(1) = nvalid
            Sbb_thr = 0.d0; SG_thr = 0.d0     ! already folded into Sbb/SG by the reduce
            pwsh_thr = 0.d0; cntsh_thr = 0.d0 ! ditto for pwsh/cntsh
            deallocate(solve_scal)
            goto 200
        endif
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz), mskrad=cov_image_mask_radius(params))
            ! gather orientations serially (get_ori is not guaranteed thread-safe)
            do i = 1, batchsz
                call build%spproj_field%get_ori(pinds(batchlims(1)+i-1), orientations(i))
            end do
            !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
            !$omp& private(i,ithr,q,r,alpha,beta,gam1,delta,a1,a2,pw,cnt,tw0) &
            !$omp& private(j,nsamp,gacc,l_gram_fast)
            do i = 1, batchsz
                if( orientations(i)%isstatezero() ) cycle
                ithr = omp_get_thread_num() + 1
                tw0 = omp_get_wtime()
                call project_fplanes_mean_basis(mean_rec, utilde, orientations(i), fpls(i), &
                    &mean_fpl(ithr), basis_fpls(:,ithr), apply_ctf_amp=.true.)
                t_proj_thr(ithr) = t_proj_thr(ithr) + (omp_get_wtime() - tw0)
                call form_residual_plane(fpls(i), mean_fpl(ithr), resid_fpl(ithr), &
                    &particle_contrast(mean_fpl(ithr), fpls(i)))
                tw0 = omp_get_wtime()
                ! signal-free noise variance from the outer residual shells
                call plane_hf_power(resid_fpl(ithr), nyq_rec, 0.7, pw, cnt)
                hfpw_thr(ithr)  = hfpw_thr(ithr)  + pw
                hfcnt_thr(ithr) = hfcnt_thr(ithr) + cnt
                ! log only: the same residual in cov_herm_inner's index convention
                call cov_herm_selfpower(resid_fpl(ithr), pw, cnt)
                cvpw_thr(ithr)  = cvpw_thr(ithr)  + pw
                cvcnt_thr(ithr) = cvcnt_thr(ithr) + cnt
                call plane_shell_power_accum(resid_fpl(ithr), nyq_rec, pwsh_thr(:,ithr), cntsh_thr(:,ithr))
                t_diag_thr(ithr) = t_diag_thr(ithr) + (omp_get_wtime() - tw0)
                tw0 = omp_get_wtime()
                do q = 1, d_tilde
                    bc_thr(q,ithr) = cov_herm_inner(basis_fpls(q,ithr), resid_fpl(ithr))
                    ! log only: Re(b)^2 against |b|^2 and G_qq
                    babs_thr(ithr) = babs_thr(ithr) + real(bc_thr(q,ithr)*conjg(bc_thr(q,ithr)),dp)
                    bre_thr(ithr)  = bre_thr(ithr)  + real(bc_thr(q,ithr))**2
                end do
                ! ---- per-particle Gram, d_tilde*(d_tilde+1)/2 inner products ----
                ! Hottest loop in flex_pca. Gather each basis plane's samples once, in cov_herm_inner's
                ! visit order, into a contiguous buffer, then form the whole Gram as one dsyrk per particle
                ! (Re(G_qr) = sum_j re_q re_r + im_q im_r). Exact: each Gi(q,r) is an independent
                ! reduction so pair order is free, and the (h,k) order WITHIN a reduction is preserved by
                ! the sample list. gdi_thr stays in q order -- that one accumulator's order does matter.
                ! Needs a common nyq across basis planes and the sample list to fit; else fall back.
                l_gram_fast = nsamp_max > 0
                if( l_gram_fast )then
                    do q = 2, d_tilde
                        if( basis_fpls(q,ithr)%nyq /= basis_fpls(1,ithr)%nyq )then
                            l_gram_fast = .false.
                            exit
                        endif
                    end do
                endif
                if( l_gram_fast )then
                    call cov_herm_sample_list(basis_fpls(1,ithr), slist(:,:,ithr), nsamp)
                    l_gram_fast = nsamp > 0
                endif
                if( l_gram_fast )then
                    ! Real interleaved layout: Re(G_qr) = sum_j (re_q re_r + im_q im_r), so stacking
                    ! [re_1, im_1, re_2, im_2, ...] as a (2*nsamp x d_tilde) real matrix X makes the whole
                    ! Gram exactly X^T X -- one dsyrk instead of d(d+1)/2 hand-rolled dot products.
                    do q = 1, d_tilde
                        do j = 1, nsamp
                            gacc = cmplx(basis_fpls(q,ithr)%cmplx_plane(slist(1,j,ithr), &
                                &slist(2,j,ithr)), kind=dp)
                            gbuf(2*j-1,q,ithr) = real(gacc)
                            gbuf(2*j,  q,ithr) = aimag(gacc)
                        end do
                    end do
                    call dsyrk('U', 'T', d_tilde, 2*nsamp, 1.d0, gbuf(1,1,ithr), size(gbuf,1), &
                        &0.d0, Gi_thr(1,1,ithr), d_tilde)
                    do q = 1, d_tilde
                        do r = q+1, d_tilde
                            Gi_thr(r,q,ithr) = Gi_thr(q,r,ithr)
                        end do
                    end do
                else
                    do q = 1, d_tilde
                        do r = q, d_tilde
                            Gi_thr(q,r,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), basis_fpls(r,ithr)), dp)
                            Gi_thr(r,q,ithr) = Gi_thr(q,r,ithr)
                        end do
                    end do
                endif
                do q = 1, d_tilde
                    gdi_thr(ithr) = gdi_thr(ithr) + Gi_thr(q,q,ithr)
                end do
                t_gram_thr(ithr) = t_gram_thr(ithr) + (omp_get_wtime() - tw0)
                do q = 1, d_tilde
                    do r = 1, d_tilde
                        ! properness: E[Re(b_q)Re(b_r)]_noise = 0.5*sig2*G, so the RHS uses Re(b), not |b|^2
                        Sbb_thr(q,r,ithr) = Sbb_thr(q,r,ithr) + real(bc_thr(q,ithr))*real(bc_thr(r,ithr))
                        SG_thr(q,r,ithr)  = SG_thr(q,r,ithr) + Gi_thr(q,r,ithr)
                    end do
                end do
                ! Stage vec(G_i) for the batched rank-k update below.
                if( cgsolve /= 0 )then
                    ! svec packing: off-diagonals carry sqrt(2) so the packing is an isometry
                    do gam1 = 1, d_tilde
                        do alpha = 1, gam1
                            Vg(svec_idx(alpha,gam1,d_tilde),i) = &
                                &merge(1.d0, sqrt(2.d0), alpha == gam1)*Gi_thr(alpha,gam1,ithr)
                        end do
                    end do
                else
                    do gam1 = 1, d_tilde
                        a2 = (gam1-1)*d_tilde
                        do alpha = 1, d_tilde
                            Vg(a2+alpha,i) = Gi_thr(alpha,gam1,ithr)
                        end do
                    end do
                endif
                nvalid_thr(ithr) = nvalid_thr(ithr) + 1
            end do
            !$omp end parallel do
            ! A_rearranged += Vg Vg^T over this batch, in ONE BLAS-3 call on the SHARED accumulator.
            tw0 = omp_get_wtime()
            if( cgsolve /= 0 )then
                call dsyrk('U', 'N', npk, batchsz, 1.d0, Vg, npk, 1.d0, Mspk, npk)
            else
                call dsyrk('U', 'N', dd, batchsz, 1.d0, Vg, dd, 1.d0, A, dd)
            endif
            Vg(:,1:batchsz) = 0.d0
            t_syrk = t_syrk + (omp_get_wtime() - tw0)
        end do
        if( flex_pca_is_worker() )then
            ! Fold this part's per-thread accumulators and ship them. The upper triangle alone is
            ! populated by dsyrk; summing upper triangles across parts and mirroring once on the master
            ! is equivalent to mirroring per part and summing, so the mirror stays master-side below.
            allocate(solve_scal(8), source=0.d0)
            do ithr = 1, nthr
                Sbb = Sbb + Sbb_thr(:,:,ithr)
                SG  = SG  + SG_thr(:,:,ithr)
                pwsh  = pwsh  + pwsh_thr(:,ithr)
                cntsh = cntsh + cntsh_thr(:,ithr)
            end do
            solve_scal(1)=sum(hfpw_thr);  solve_scal(2)=sum(hfcnt_thr)
            solve_scal(3)=sum(cvpw_thr);  solve_scal(4)=sum(cvcnt_thr)
            solve_scal(5)=sum(babs_thr);  solve_scal(6)=sum(bre_thr); solve_scal(7)=sum(gdi_thr)
            nvalid = sum(nvalid_thr)
            if( cgsolve /= 0 )then
                call write_solve_part(flex_pca_part_fname('solve', params%part, params%numlen), &
                    &Mspk, npk, Sbb, SG, d_tilde, pwsh, cntsh, nyq_rec, solve_scal, nvalid)
            else
                call write_solve_part(flex_pca_part_fname('solve', params%part, params%numlen), &
                    &A, dd, Sbb, SG, d_tilde, pwsh, cntsh, nyq_rec, solve_scal, nvalid)
            endif
            deallocate(solve_scal)
            ncomp_out = 0
            sig2_out  = 0._dp
            allocate(vred(0,0), gamma_out(0))
            return
        endif
200     continue
        ! dsyrk touched only the upper triangle; mirror it before the un-rearrangement.
        if( cgsolve /= 0 )then
            !$omp parallel do default(shared) private(a1,a2) schedule(static)
            do a2 = 1, npk
                do a1 = a2+1, npk
                    Mspk(a1,a2) = Mspk(a2,a1)
                end do
            end do
            !$omp end parallel do
            ! Rearrange ONCE, here, so CG below is a dsymv per iteration rather than a scalar gather
            ! per term. Mirrors the dense branch, which unrearranges at exactly this point -- so the
            ! ridge and the trace diagnostics downstream now see the OPERATOR's diagonal in both
            ! paths, where the packed path previously scaled its ridge by the raw accumulator's.
            t_phase = tic()
            call rearrange_packed_selfsum(Mspk, npk, d_tilde)
            write(logfhandle,'(A,F8.2)') '>>> FLEX_PCA packed rearrangement seconds=', toc(t_phase)
            call flush(logfhandle)
        else
            !$omp parallel do default(shared) private(a1,a2) schedule(static)
            do a2 = 1, dd
                do a1 = a2+1, dd
                    A(a1,a2) = A(a2,a1)
                end do
            end do
            !$omp end parallel do
            call unrearrange_kron_selfsum(A, d_tilde)
        endif
        ! reduce per-thread accumulators
        do ithr = 1, nthr
            Sbb = Sbb + Sbb_thr(:,:,ithr)
            SG  = SG  + SG_thr(:,:,ithr)
        end do
        nvalid = sum(nvalid_thr)
        write(logfhandle,'(A,F9.1,A,F9.1,A,F9.1,A,F9.1)') &
            &'>>> FLEX_PCA SOLVE SPLIT (thread-seconds): projections=',sum(t_proj_thr), &
            &'  gram=',sum(t_gram_thr),'  diagnostics=',sum(t_diag_thr),'  dsyrk(wall)=',t_syrk
        call flush(logfhandle)
        write(logfhandle,'(A,I0,A,F10.1)') '>>> FLEX_PCA REDUCED SOLVE particles=',nvalid, &
            &' seconds=',toc(t_phase)
        call flush(logfhandle)
        do i = 1, size(orientations)
            call orientations(i)%kill
        end do
        do ithr = 1, nthr
            call cleanup_plane(mean_fpl(ithr)); call cleanup_plane(resid_fpl(ithr))
            do q = 1, d_tilde
                call cleanup_plane(basis_fpls(q,ithr))
            end do
        end do
        call cleanup_rec_buffers(build, fpls)
        deallocate(basis_fpls, mean_fpl, resid_fpl, Gi_thr, bc_thr, Vg, Sbb_thr, SG_thr, orientations, nvalid_thr)
        deallocate(gbuf, slist)
        if( nvalid < 1 ) THROW_HARD('flex_pca reduced covariance solve found no valid particles')
        ! Noise-bias scale sig2_eff = signal-free per-coefficient whitened-noise variance E|n|^2, measured
        ! from the outer residual shells. With the properness-consistent RHS the debias is
        ! D = Sbb - 0.5*sig2_eff*SG (below), which removes the noise bias WITHOUT over-subtracting the
        ! signal, unlike trace matching, which folds the signal into sig2_eff and collapses the rank.
        sig2_eff = sum(hfpw_thr) / max(sum(hfcnt_thr), 1.d0)
        ! per-shell whitened residual variance profile, for the log
        do ithr = 1, nthr
            pwsh  = pwsh  + pwsh_thr(:,ithr)
            cntsh = cntsh + cntsh_thr(:,ithr)
        end do
        write(logfhandle,'(A)') '>>> FLEX_PCA per-shell whitened residual variance (shell: var / sig2_eff):'
        do q = 0, min(nyq_rec,32)
            if( cntsh(q) > 0.d0 ) write(logfhandle,'(A,I3,A,ES11.3,A,F8.3)') '>>>   sh=',q, &
                &'  var=',pwsh(q)/cntsh(q),'  ratio=',(pwsh(q)/cntsh(q))/max(sig2_eff,DTINY)
        end do
        call flush(logfhandle)
        tr_bb = 0.d0; tr_g = 0.d0
        do q = 1, d_tilde
            tr_bb = tr_bb + Sbb(q,q)
            tr_g  = tr_g  + SG(q,q)
        end do
        ! The trace-match estimate tr(Sbb)/tr(SG) is only meaningful when tr(SG) is a real quantity rather
        ! than a rounding residue: with a numerically dead projected basis tr_g underflows any fixed floor
        ! and the ratio then reports the floor, not a noise scale. Gate it on the relative epsilon instead.
        l_trmatch = tr_g > epsilon(1.d0) * max(abs(tr_bb), 1.d0)
        tr_match  = 0.d0
        if( l_trmatch ) tr_match = tr_bb / tr_g
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,ES12.4)') '>>> FLEX_PCA noise sig2_eff(hf)=',sig2_eff, &
            &'  tr(Sbb)=',tr_bb,'  tr(SG)=',tr_g
        if( l_trmatch )then
            write(logfhandle,'(A,ES12.4)') '>>> FLEX_PCA trace-match sig2=tr(Sbb)/tr(SG)=',tr_match
        else
            write(logfhandle,'(A)') '>>> FLEX_PCA trace-match sig2 UNDEFINED: tr(SG) is at the rounding &
                &floor -- the projected basis is numerically dead; suspect a reconstructor loaded with add() &
                &while flagged Fourier (see the D2 zeroG/zeroRHS/zeroZ counters).'
        endif
        ! log only: the debias identity E[Re(b)Re(b)] = 0.5*sig2*G requires sig2 in THIS convention
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,ES12.4)') '>>> FLEX_PCA D1b mean|b|^2=', &
            &sum(babs_thr)/max(real(nvalid*d_tilde,dp),1.d0),'  meanRe(b)^2=', &
            &sum(bre_thr)/max(real(nvalid*d_tilde,dp),1.d0),'  meanG_qq=', &
            &sum(gdi_thr)/max(real(nvalid*d_tilde,dp),1.d0)
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,ES12.4)') '>>> FLEX_PCA D1 sig2_covherm=', &
            &sum(cvpw_thr)/max(sum(cvcnt_thr),1.d0),'  0.5*sig2_covherm=', &
            &0.5d0*sum(cvpw_thr)/max(sum(cvcnt_thr),1.d0),'  measured Sbb/SG=',tr_match
        if( .not. l_trmatch ) write(logfhandle,'(A)') &
            &'>>> FLEX_PCA D1 measured Sbb/SG is reported as 0: tr(SG) is at the rounding floor'
        allocate(Dmat(d_tilde,d_tilde))
        ! SIMPLE's FFT/gridding convention is not unitary, so after whitening the noise covariance is
        ! Lambda_i = sig2_conv * I with sig2_conv a pure convention constant (box, padding, gen_fplane4rec
        ! scaling), NOT the supplement's 1. The 0.5 is exact: it is E[Re(b)Re(b)] for proper complex noise.
        Dmat = Sbb - 0.5d0 * sig2_eff * SG
        ! The same identity fixes the embedding's noise scale: b = G z + n with Cov(n) = 0.5*sig2*G
        sig2_out = 0.5d0 * sig2_eff
        write(logfhandle,'(A,ES12.4,A,ES12.4)') '>>> FLEX_PCA noise convention constant sig2_eff=', &
            &sig2_eff,'  debias coefficient=',0.5d0*sig2_eff
        call flush(logfhandle)
        trc = 0.d0
        if( cgsolve /= 0 )then
            ! ridge scaled to the packed accumulator's diagonal
            do a1 = 1, npk
                trc = trc + Mspk(a1,a1)
            end do
            ridge = 1.d-6 * trc / real(npk,dp)
        else
            do a1 = 1, dd
                trc = trc + A(a1,a1)
            end do
            ridge = 1.d-6 * trc / real(dd,dp)
            do a1 = 1, dd
                A(a1,a1) = A(a1,a1) + ridge
            end do
        endif
        ! solve A x = vec(D)
        allocate(bvec(dd), xsol(dd))
        do q = 1, d_tilde
            do r = 1, d_tilde
                bvec((r-1)*d_tilde + q) = Dmat(q,r)
            end do
        end do
        ! A = sum_i G_i (x) G_i + ridge*I is symmetric POSITIVE DEFINITE by construction (a sum of Kronecker
        ! squares of symmetric matrices is PSD), so Cholesky is the right factorization
        if( cgsolve /= 0 )then
            ! pack the RHS and solve matrix-free
            allocate(rhspk(npk), spk(npk))
            do r = 1, d_tilde
                do q = 1, r
                    rhspk(svec_idx(q,r,d_tilde)) = merge(1.d0, sqrt(2.d0), q == r)*Dmat(q,r)
                end do
            end do
            call cg_solve_packed(Mspk, npk, rhspk, ridge, COV_CG_MAXIT, COV_CG_TOL, &
                &spk, cgit, cgres)
            write(logfhandle,'(A,I0,A,ES10.3)') '>>> FLEX_PCA CG iterations=',cgit, &
                &'  relative residual=',cgres
            if( cgres > COV_CG_TOL )then
                write(logfhandle,'(A)') '>>> FLEX_PCA WARNING: CG did not reach tolerance; &
                    &Sigma is not converged and the eigenvalues below should not be trusted'
            endif
            call flush(logfhandle)
            ! unpack svec -> full symmetric, undoing the sqrt(2) scaling
            do r = 1, d_tilde
                do q = 1, r
                    xsol((r-1)*d_tilde + q) = spk(svec_idx(q,r,d_tilde)) / merge(1.d0, sqrt(2.d0), q == r)
                    xsol((q-1)*d_tilde + r) = xsol((r-1)*d_tilde + q)
                end do
            end do
            deallocate(rhspk, spk)
            info = 0
        else
        xsol = bvec
        call dposv('U', dd, 1, A, dd, xsol, dd, info)
        endif
        if( info /= 0 )then
            ! A is destroyed by a failed factorization, so there is nothing to retry with
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA reduced-solve Cholesky failed at pivot ',info, &
                &'; the projected-covariance system is not positive definite'
            THROW_HARD('flex_pca reduced covariance solve is not positive definite')
        endif
        allocate(Sig(d_tilde,d_tilde))
        do q = 1, d_tilde
            do r = 1, d_tilde
                Sig(q,r) = xsol((r-1)*d_tilde + q)
            end do
        end do
        Sig = 0.5d0 * (Sig + transpose(Sig))          ! symmetrize
        ! HeteroPCA diagonal deletion + imputation, in place of trusting the analytic debias.
        if( .not. cov_env_int_off('SIMPLE_COV_HETPCA') ) &
            &call heteropca_impute(Sig, d_tilde, min(neigs_req, d_tilde))
        ! eigendecompose Sigma = V Gamma V^T
        allocate(V(d_tilde,d_tilde), gam(d_tilde))
        call jacobi(Sig, d_tilde, d_tilde, gam, V, nrot)
        call eigsrt(gam, V, d_tilde, d_tilde)          ! descending
        ! Higham PSD projection: clamp the negative eigenvalues the debias can produce
        do q = 1, d_tilde
            if( gam(q) < 0.d0 ) gam(q) = 0.d0
        end do
        gmax = max(gam(1), DTINY)
        keep = 0
        do q = 1, d_tilde
            if( gam(q) > COV_EIG_REL_FLOOR*gmax ) keep = keep + 1
        end do
        ncomp = max(1, min(neigs_req, keep))
        ncomp_out = ncomp
        allocate(vred(d_tilde,ncomp), gamma_out(ncomp))
        do q = 1, ncomp
            vred(:,q)    = V(:,q)
            gamma_out(q) = max(gam(q), DTINY)
        end do
        write(logfhandle,'(A,ES12.4,A,ES12.4)') '>>> FLEX_PCA reduced covariance eigenvalues: max=', &
            &gamma_out(1),' min=',gamma_out(ncomp)
        call flush(logfhandle)
        if( allocated(A) )    deallocate(A)
        if( allocated(Mspk) ) deallocate(Mspk)
        deallocate(Sbb, SG, Dmat, Sig, V, gam, xsol, bvec)
        deallocate(pwsh_thr, cntsh_thr, pwsh, cntsh)
    end subroutine reduced_covariance_solve

    !> Final eigenbasis WITH the supplement S.D contrast correction (Algorithm 3).
    subroutine form_eigenbasis_from_reduced( params, build, mean_rec, utilde_real, d_tilde, vred, gamma, &
        &ncomp, basis_recs, basis_imgs, fprefix )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(in)    :: mean_rec
        type(image),         intent(inout) :: utilde_real(d_tilde)
        integer,             intent(in)    :: d_tilde, ncomp
        real(dp),            intent(in)    :: vred(d_tilde,ncomp)
        !> covariance eigenvalues; overwritten by the contrast-corrected Gamma = S^2
        real(dp),            intent(inout) :: gamma(ncomp)
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        !> clean real-space eigenvolumes, retained for cross-halfset basis alignment (held-out embedding)
        !! and for any downstream volume-space inner product.
        type(image), allocatable, optional, intent(out) :: basis_imgs(:)
        character(len=*),        optional, intent(in)  :: fprefix
        type(image), allocatable :: cvols(:)
        type(image)  :: eigimg, meanimg
        type(string) :: fname
        real, pointer :: qr(:,:,:), ck(:,:,:), cl(:,:,:)
        real, pointer :: esgn(:,:,:)
        integer :: eloc(3)
        real(dp), allocatable :: gramC(:,:), W(:,:), lam(:)
        real(dp) :: qnrm2, ip, lam_m, max_ov, cnrm
        integer :: k, j, m, nrot
        ! ---- materialize C = Utilde V Gamma^(1/2), one volume per component ----
        allocate(cvols(ncomp))
        do k = 1, ncomp
            call cvols(k)%copy(utilde_real(1))
            call cvols(k)%zero_and_unflag_ft
            do j = 1, d_tilde
                call cvols(k)%add(utilde_real(j), real(vred(j,k)*sqrt(max(gamma(k),DTINY))))
            end do
        end do
        ! ---- project the mean component out of C (S.D, Algorithm 3) ----
        call meanimg%copy(mean_rec)
        if( meanimg%is_ft() ) call meanimg%ifft
        ! the column representatives were soft-masked in realize_hermitian_volume, so the mean
        ! must be masked identically or the two live in different inner-product spaces
        if( params%msk_crop > TINY ) call meanimg%mask3D_soft(params%msk_crop, backgr=0.)
        call meanimg%get_rmat_ptr(qr)
        qnrm2  = sum(real(qr,dp)*real(qr,dp))
        max_ov = 0.d0
        if( qnrm2 > DTINY )then
            do k = 1, ncomp
                call cvols(k)%get_rmat_ptr(ck)
                ip   = sum(real(ck,dp)*real(qr,dp))
                cnrm = sqrt(sum(real(ck,dp)*real(ck,dp)))
                if( cnrm > DTINY ) max_ov = max(max_ov, abs(ip)/(sqrt(qnrm2)*cnrm))
                ck = ck - real(ip/qnrm2)*qr
            end do
            write(logfhandle,'(A,F8.5)') '>>> FLEX_PCA S.D mean-component projection: &
                &max overlap removed=',real(max_ov)
        else
            write(logfhandle,'(A)') '>>> FLEX_PCA S.D mean-component projection skipped: empty mean'
        endif
        call flush(logfhandle)
        call meanimg%kill
        ! ---- svd(C) through the Gram eigendecomposition; Gamma <- S^2 ----
        allocate(gramC(ncomp,ncomp), W(ncomp,ncomp), lam(ncomp))
        do k = 1, ncomp
            call cvols(k)%get_rmat_ptr(ck)
            do j = k, ncomp
                call cvols(j)%get_rmat_ptr(cl)
                gramC(k,j) = sum(real(ck,dp)*real(cl,dp))
                gramC(j,k) = gramC(k,j)
            end do
        end do
        call jacobi(gramC, ncomp, ncomp, lam, W, nrot)
        call eigsrt(lam, W, ncomp, ncomp)          ! descending
        allocate(basis_recs(ncomp))
        if( present(basis_imgs) ) allocate(basis_imgs(ncomp))
        do m = 1, ncomp
            k     = m                              ! output naming/index follows the sorted order
            lam_m = max(lam(m), DTINY)
            ! left singular vector u_m = lambda_m^(-1/2) sum_k W(k,m) c_k
            call eigimg%copy(utilde_real(1))
            call eigimg%zero_and_unflag_ft
            do j = 1, ncomp
                call eigimg%add(cvols(j), real(W(j,m)/sqrt(lam_m)))
            end do
            gamma(m) = lam_m
            if( params%msk_crop > TINY ) call eigimg%mask3D_soft(params%msk_crop, backgr=0.)
            ! CANONICAL SIGN. Eigenvectors are defined up to sign and the solver's choice depends on the
            ! last bits of the matrix, so the same data reduced in a different order (nparts, threads,
            ! BLAS) gives the same subspace with arbitrary signs. Not cosmetic: the state ordering uses
            ! the ordering eigenvolume's sign as the trajectory DIRECTION, so a flip plays the motion
            ! backwards and swaps flex_pca_traj_001 with _00N. Making the dominant voxel positive is
            ! stable -- the leading lobe sits far above the 1e-6 perturbations that flip the solver.
            call eigimg%get_rmat_ptr(esgn)
            eloc = maxloc(abs(esgn))
            if( esgn(eloc(1),eloc(2),eloc(3)) < 0. ) call eigimg%mul(-1.)
            if( present(fprefix) )then
                fname = trim(fprefix)//int2str_pad(k,3)//MRC_EXT
            else
                fname = 'flex_pca_pc'//int2str_pad(k,3)//MRC_EXT
            endif
            call eigimg%write(fname, del_if_exists=.true.)
            call fname%kill
            if( present(basis_imgs) ) call basis_imgs(k)%copy(eigimg)
            ! Reconstructor for the embedding projection, built with the mean_rec idiom: set_rmat, fft,
            ! expand_exp. NEVER add() -- it deposits into the shared rmat/cmat buffer of a reconstructor
            ! still flagged Fourier, the following fft() is then a no-op, expand_exp propagates an
            ! untransformed grid, and the projected basis comes out identically zero.
            call init_column_reconstructor(params, build, basis_recs(k))
            call basis_recs(k)%set_rmat(eigimg%get_rmat(), .false.)
            call basis_recs(k)%fft
            call basis_recs(k)%expand_exp
        end do
        do k = 1, ncomp
            call cvols(k)%kill
        end do
        deallocate(cvols, gramC, W, lam)
        call eigimg%kill
        write(logfhandle,'(A,ES12.4,A,ES12.4)') '>>> FLEX_PCA S.D-corrected eigenvalues: max=', &
            &gamma(1),' min=',gamma(ncomp)
        call flush(logfhandle)
    end subroutine form_eigenbasis_from_reduced

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
        logical :: l_relprior, l_stats_only, l_from_parts
        integer :: ihf
        integer :: batchlims(2), batchsz, ibatch, i, iptcl, q, r, ithr, nthr, ia, row
        integer, allocatable :: nzeroG_thr(:), nzeroR_thr(:), nzeroZ_thr(:)   ! dead-basis counters
        real(dp) :: a, a_best, e_yy, e_mm, best_res, res, aa, sig2
        integer(timer_int_kind) :: t_phase
        real(dp), allocatable :: Gth(:,:,:), Ath(:,:,:), zth(:,:), zbest(:,:), cth(:,:), bth(:,:), myth(:)
        real(dp), allocatable :: Gtilth(:,:,:)   ! per-thread noise-whitened projected Gram
        real(dp), allocatable :: gwork(:,:,:), gvec(:,:,:), gev(:,:), gspec_thr(:,:), gspec(:)
        integer,  allocatable :: gcnt_thr(:)
        integer :: nrot_t, gcnt
        real(dp) :: gsum
        nthr = nthr_glob
        sig2 = max(sig2_eff, DTINY)      ! whitened-noise variance for the MAP shrinkage
        allocate(nzeroG_thr(nthr), nzeroR_thr(nthr), nzeroZ_thr(nthr), source=0)   ! dead-basis counters
        l_relprior = .not. cov_env_int_off('SIMPLE_COV_RELPRIOR')
        l_stats_only = .false.
        l_from_parts = .false.
        if( present(stats_only) ) l_stats_only = stats_only
        if( present(from_parts) ) l_from_parts = from_parts
        ! Distribution rides on the sufficient statistics the reliability prior needs. With
        ! SIMPLE_COV_RELPRIOR=0 those caches are never formed, so there is nothing to ship and the
        ! stage stays in process -- correct, just not parallel.
        if( .not. l_relprior .and. (l_stats_only .or. l_from_parts) ) &
            &THROW_HARD('distributed embedding requires the reliability prior; unset SIMPLE_COV_RELPRIOR')
        allocate(prior(ncomp))
        if( l_relprior )then
            ! The reducing master never holds Gcache at nptcls: it reads one part's blocks at a time in
            ! the re-solve below. At a million particles with ncomp=48 that is the difference between
            ! 1.8 GB and 18.4 GB, on top of the precision array it must hold either way.
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
        ! MASTER REDUCING PARTS. The image pass below is the whole cost of this stage and the workers
        ! have already paid it; the master only needs their sufficient statistics to run the coupled
        ! phase (rho over every particle, then the re-solve). Jump straight there.
        if( l_from_parts )then
            ! PASS 1 ONLY. zhalf and the per-particle scalars are all rho needs, and rho has to exist
            ! before any particle can be solved. The Gram blocks stay on disk until the re-solve reads
            ! them one part at a time -- Gcache below is deliberately never allocated at nptcls.
            call reduce_embed_zhalf_parts(params, flex_pca_nparts(), pinds, contrast, resid_energy, &
                &resid_mean_energy, zhalf, nptcls, ncomp)
            goto 200
        endif
        write(logfhandle,'(A)') '>>> FLEX_PCA CONTRAST-AWARE EMBEDDING'
        call flush(logfhandle)
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
            !$omp& private(i,ithr,q,r,ia,a,a_best,aa,e_yy,e_mm,best_res,res,row)
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
                do ia = 1, NCONTRAST_GRID
                    if( COV_EMBED_CONTRAST_GRID )then
                        a = A_GRID_HI * real(ia,dp) / real(NCONTRAST_GRID,dp)   ! 0.04 .. 2.0
                    else
                        a = 1.d0
                    endif
                    aa = a*a
                    Ath(:,:,ithr) = (aa/sig2)*Gth(:,:,ithr)
                    do q = 1, ncomp
                        Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                        zth(q,ithr)   = (a*bth(q,ithr) - aa*cth(q,ithr))/sig2
                    end do
                    call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), ncomp)
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
                        a_best        = a
                        zbest(:,ithr) = zth(:,ithr)
                    endif
                    if( .not. COV_EMBED_CONTRAST_GRID ) exit
                end do
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
                if( l_relprior )then
                    ! cache the sufficient statistics so the reliability-corrected prior can be
                    ! applied by re-solving in closed form -- no second pass over the images
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
        ! The reducing master lands here rather than after the diagnostics, so it still frees the
        ! per-thread Gram workspace. gcnt_thr is all zero when no batch loop ran, so the spectrum
        ! report below skips itself -- those diagnostics are per-particle and belong to the parts.
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
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA D2 zeroG=',sum(nzeroG_thr), &
            &' zeroRHS=',sum(nzeroR_thr),' zeroZ=',sum(nzeroZ_thr),' of nptcls=',nptcls
        write(logfhandle,'(A,F8.1)') '>>> FLEX_PCA CONTRAST EMBED SECONDS: ', toc(t_phase)
        call flush(logfhandle)
        ! WORKER. The image pass is done for this part's particles; everything below is coupled across
        ! ALL of them, so it belongs to the master. Ship the sufficient statistics and stop.
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
        if( l_relprior )then
            allocate(rho(ncomp))
            ! OPTIONAL EXPORT of the two half-data solves themselves, for calibrating the per-particle
            ! error MODEL against the error actually observed. rho below collapses the pair to one
            ! correlation per component; the variance of their DIFFERENCE is what measures the error
            ! directly, including the model misspecification that precision(:,:,row) cannot see.
            ! Off unless SIMPLE_COV_ZHALF is set, and writes nothing else -- purely additive.
            call write_zhalf_replicates(zhalf, prior, nptcls, ncomp)
            do q = 1, ncomp
                rho(q) = corr_dp(zhalf(:,q,1), zhalf(:,q,2), nptcls)
                rho(q) = max(0.d0, rho(q))
                rho(q) = 2.d0*rho(q) / (1.d0 + rho(q))            ! Spearman-Brown to full length
            end do
            ! Scale rho RELATIVE to the most reliable component, not absolutely: an absolute rho^2 shrinks
            ! the informative components as well and compresses the latent spread the trajectory needs.
            rho_max = maxval(rho)
            if( rho_max <= DTINY ) rho_max = 1.d0
            do q = 1, ncomp
                rrel = (rho(q)*rho(q)) / (rho_max*rho_max)
                prior(q) = 1.d0 / max(max(rrel, RHO_FLOOR) * eigvals(q), DTINY)
            end do
            ! All components, not the leading 10: rho drives state-target placement via comp_rho, so
            ! its ranking has to be checkable against a signal measure over the whole basis.
            write(logfhandle,'(A)') '>>> FLEX_PCA split-half reliability per component (rho, corrected):'
            do q = 1, ncomp
                write(logfhandle,'(A,I3,A,F7.4,A,ES11.3,A,ES11.3)') '>>>   z',q,'  rho=',rho(q), &
                    &'  eigval=',eigvals(q),'  prior_precision=',prior(q)
            end do
            call flush(logfhandle)
            ! re-solve every particle in closed form from the cached sufficient statistics.
            ! Two routes to the SAME arithmetic: in process the blocks are already in Gcache, while a
            ! reducing master streams them back one part at a time so its footprint does not grow with
            ! the dataset. The solve is identical either way -- only where the blocks live differs.
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
            write(logfhandle,'(A)') '>>> FLEX_PCA latents re-solved with the reliability-weighted prior'
            call flush(logfhandle)
            if( present(rho_out) ) rho_out = rho
            deallocate(rho, Gcache, bcache, ccache, zhalf, Ghf, bhf, chf)
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
        type(image) :: img_o
        real,                allocatable :: rho_e(:,:,:,:), rho_o(:,:,:,:), filt(:), corrs(:)
        real(dp),            allocatable :: zbatch(:,:), dens(:,:,:), z(:,:), prior(:)
        real(dp),            allocatable :: Gth(:,:,:), Ath(:,:,:), bth(:,:), cth(:,:), zth(:,:)
        !> per-thread posterior covariance A^-1 and the untouched copy of A it is taken from
        real(dp),            allocatable :: Ainvth(:,:,:), Acpth(:,:,:)
        !> per-thread accumulator for the EM Gamma update, and its reduction
        real(dp),            allocatable :: gam_thr(:,:), gam_acc(:)
        integer,             allocatable :: nval_thr(:)
        logical,             allocatable :: valid(:), valid_e(:), valid_o(:)
        integer,             allocatable :: eo(:)
        real, pointer :: rmatp(:,:,:)
        real     :: fc
        real(dp) :: sig2, a, aa, e_mm, myv, mu_q, sd_q
        integer  :: it, q, r, i, ithr, nthr, batchlims(2), batchsz, ibatch, row, d_new, es(3), filtsz, sh
        integer  :: npairs, nval, pstride, npp, ihalf, nkept, nprobe_max
        logical  :: l_probe_distr
        real(dp),            allocatable :: gam_sum(:)
        complex,             allocatable :: cme(:,:,:,:), cmo(:,:,:,:)
        real,                allocatable :: rhe(:,:,:,:), rhoo(:,:,:,:)
        integer,             allocatable :: ppinds(:)
        type(image),         allocatable :: prev_real(:)      ! previous iteration's orthonormal basis
        real(dp),            allocatable :: Mconv(:,:), sconv(:)
        real(dp) :: cos_mean
        logical  :: l_converged
        type(string) :: fname
        integer(timer_int_kind) :: t_it
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
        ! ---- ABSOLUTE CAP, not a fixed ratio ----
        ! The probe refines ncomp band-limited, FSC-regularised volumes. That parameter count is FIXED
        ! -- it does not grow when the dataset does -- so the number of particles needed to determine
        ! it is fixed too, and a constant stride is the wrong control: at stride 4 a million-particle
        ! project would still refine on 250k particles and the probe would still scale linearly,
        ! staying the dominant stage (projected 5228 s of a 4.1 h run). Capping the COUNT instead
        ! makes the probe O(1) in dataset size.
        !
        ! COV_PROBE_MAX_PTCLS is the total across all processes, so each one takes its share; a worker
        ! sees only its own partition and would otherwise take the whole budget nparts times over.
        ! 25000 is what the stride-4 benchmark actually validated at 100k particles: basis principal
        ! angles 0.988557 vs 0.989457 against every-particle refinement, and 71 vs 70 delivered states.
        ! SIMPLE_COV_PROBE_STRIDE still wins if set, so the old knob keeps working.
        nprobe_max = max(1, COV_PROBE_MAX_PTCLS / max(1, params%nparts))
        pstride    = 1
        if( nptcls > nprobe_max ) pstride = (nptcls + nprobe_max - 1) / nprobe_max
        call cov_env_int('SIMPLE_COV_PROBE_STRIDE', pstride)
        pstride = max(1, pstride)
        allocate(ppinds(nptcls))
        npp = 0
        do ihalf = 0, 1
            nkept = 0
            do i = 1, nptcls
                if( build%spproj_field%get_eo(pinds(i)) /= ihalf ) cycle
                if( mod(nkept, pstride) == 0 )then
                    npp = npp + 1
                    ppinds(npp) = pinds(i)
                endif
                nkept = nkept + 1
            end do
        end do
        if( npp < 2 ) THROW_HARD('probe stride left too few particles; lower SIMPLE_COV_PROBE_STRIDE')
        call hpsort(ppinds(:npp))   ! restore project order so batched image reads stay sequential
        if( pstride > 1 )then
            write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA PROBE stride ',pstride, &
                &': refining the basis on ',npp,' of ',nptcls,' particles'
            call flush(logfhandle)
        endif
        allocate(z(npp,ncomp))
        do it = 1, niters
            t_it = tic()
            write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA PROBE SUBSPACE ITERATION ',it,' / ',niters, &
                &'  basis dim=',ncomp
            call flush(logfhandle)
            allocate(prior(ncomp))
            do q = 1, ncomp
                prior(q) = 1.d0 / max(eigvals(q), DTINY)
            end do
            ! even/odd Y_q accumulators (half-set FSC regularization) + the COUPLED latent normal matrix.
            ! rho carries one entry per (q,r) pair, not one shared density: the M-step below solves the
            ! components together at every grid point.
            allocate(Yeven(ncomp), Yodd(ncomp))
            do q = 1, ncomp
                call init_column_reconstructor(params, build, Yeven(q)); call Yeven(q)%reset; call Yeven(q)%reset_exp
                call init_column_reconstructor(params, build, Yodd(q));  call Yodd(q)%reset;  call Yodd(q)%reset_exp
            end do
            es     = shape(Yeven(1)%cmat_exp)
            npairs = (ncomp*(ncomp+1))/2
            allocate(rho_e(npairs,es(1),es(2),es(3)), rho_o(npairs,es(1),es(2),es(3)), source=0.)
            write(logfhandle,'(A,I0,A,F8.2,A)') '>>> FLEX_PCA PROBE coupled normal matrix pairs=',npairs, &
                &'  rho even+odd ', 8.d0*real(npairs,dp)*real(es(1),dp)*real(es(2),dp)*real(es(3),dp)/1.d9,' GB'
            call flush(logfhandle)
            allocate(Gth(ncomp,ncomp,nthr), Ath(ncomp,ncomp,nthr), bth(ncomp,nthr), cth(ncomp,nthr), zth(ncomp,nthr))
            allocate(Ainvth(ncomp,ncomp,nthr), Acpth(ncomp,ncomp,nthr))
            allocate(basis_fpls(ncomp,nthr), mean_fpl(nthr), orientations(MAXIMGBATCHSZ))
            allocate(zbatch(ncomp,MAXIMGBATCHSZ), dens(ncomp,ncomp,MAXIMGBATCHSZ))
            allocate(valid(MAXIMGBATCHSZ), valid_e(MAXIMGBATCHSZ), valid_o(MAXIMGBATCHSZ), eo(MAXIMGBATCHSZ))
            allocate(gam_thr(ncomp,nthr), source=0.d0)
            allocate(gam_acc(ncomp), source=0.d0)
            allocate(nval_thr(nthr), source=0)
            dens = 0.d0
            call mean_rec%expand_exp
            do q = 1, ncomp
                call basis_recs(q)%expand_exp
            end do
            call init_rec(params, build, MAXIMGBATCHSZ, fpls)
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
                call reduce_probe_parts(params, flex_pca_nparts(), cme, rhe, cmo, rhoo, &
                    &rho_e, rho_o, gam_sum, nval, ncomp)
                do q = 1, ncomp
                    Yeven(q)%cmat_exp = cme(:,:,:,q); Yeven(q)%rho_exp = rhe(:,:,:,q)
                    Yodd(q)%cmat_exp  = cmo(:,:,:,q); Yodd(q)%rho_exp  = rhoo(:,:,:,q)
                end do
                deallocate(cme, cmo, rhe, rhoo)
            else
            call prepimgbatch(params, build, MAXIMGBATCHSZ)
            z = 0.d0
            do ibatch = 1, npp, MAXIMGBATCHSZ
                batchlims = [ibatch, min(npp, ibatch + MAXIMGBATCHSZ - 1)]
                batchsz   = batchlims(2) - batchlims(1) + 1
                call discrete_read_imgbatch(params, build, npp, ppinds, batchlims)
                call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                    &ppinds(batchlims(1):batchlims(2)), fpls(:batchsz), mskrad=cov_image_mask_radius(params))
                do i = 1, batchsz
                    call build%spproj_field%get_ori(ppinds(batchlims(1)+i-1), orientations(i))
                    eo(i) = build%spproj_field%get_eo(ppinds(batchlims(1)+i-1))
                end do
                valid(:batchsz) = .false.
                zbatch(:,:batchsz) = 0.d0
                dens(:,:,:batchsz) = 0.d0
                !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
                !$omp& private(i,ithr,q,r,a,aa,e_mm,myv,row)
                do i = 1, batchsz
                    if( orientations(i)%isstatezero() ) cycle
                    ithr = omp_get_thread_num() + 1
                    row  = batchlims(1) + i - 1
                    call project_fplanes_mean_basis(mean_rec, basis_recs, orientations(i), fpls(i), &
                        &mean_fpl(ithr), basis_fpls(:,ithr), apply_ctf_amp=.true.)
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
                    ! Posterior precision A = (a^2/sig2) G + Gamma^-1. The whole normal system is already
                    ! scaled by 1/sig2, so Cov[z|y] = A^-1 exactly -- no further sig2 factor.
                    Ath(:,:,ithr) = (aa/sig2)*Gth(:,:,ithr)
                    do q = 1, ncomp
                        Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                        zth(q,ithr)   = (a*bth(q,ithr) - aa*cth(q,ithr))/sig2
                    end do
                    ! A^-1 comes off a COPY: spd_solve_dp rescales and ridges A in place.
                    Acpth(:,:,ithr) = Ath(:,:,ithr)
                    call spd_inv_dp(Acpth(:,:,ithr), Ainvth(:,:,ithr), ncomp)
                    call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), ncomp)
                    z(row,:)          = zth(:,ithr)
                    zbatch(:,i)       = zth(:,ithr)
                    ! EM sufficient statistic E[z z'|y] = z z' + Cov[z|y]. BOTH the coupled M-step normal
                    ! matrix and the Gamma update below need it. Dropping Cov underestimates Gamma, which
                    ! tightens the prior, which shrinks z further: the bias compounds across iterations.
                    do r = 1, ncomp
                        do q = 1, ncomp
                            dens(q,r,i) = zth(q,ithr)*zth(r,ithr) + Ainvth(q,r,ithr)
                        end do
                    end do
                    do q = 1, ncomp
                        gam_thr(q,ithr) = gam_thr(q,ithr) + dens(q,q,i)
                    end do
                    nval_thr(ithr)    = nval_thr(ithr) + 1
                    valid(i)          = .true.
                    ! residual observation r_i = y - a*(T mu) in place (transfer/ctfsq intact for backprojection)
                    fpls(i)%cmplx_plane = fpls(i)%cmplx_plane - real(a)*mean_fpl(ithr)%cmplx_plane
                end do
                !$omp end parallel do
                ! M-step by halfset: Y_q += sum_i z_iq * backproject(r_i), and the coupled normal matrix
                ! rho(q,r) += sum_i |CTF|^2 E[z_iq z_ir]   (batched KB)
                do i = 1, batchsz
                    valid_e(i) = valid(i) .and. eo(i)==0
                    valid_o(i) = valid(i) .and. eo(i)==1
                end do
                call insert_planes_oversamp_coupled_batch_scaled(Yeven, rho_e, build%pgrpsyms, &
                    &orientations(:batchsz), fpls(:batchsz), zbatch(:,:batchsz), dens(:,:,:batchsz), &
                    &valid_e(:batchsz), batchsz)
                call insert_planes_oversamp_coupled_batch_scaled(Yodd, rho_o, build%pgrpsyms, &
                    &orientations(:batchsz), fpls(:batchsz), zbatch(:,:batchsz), dens(:,:,:batchsz), &
                    &valid_o(:batchsz), batchsz)
                if( batchlims(2)==nptcls .or. mod(batchlims(2), 5*MAXIMGBATCHSZ)==0 )then
                    write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA PROBE PASS PARTICLES: ',batchlims(2),' / ',npp
                    call flush(logfhandle)
                endif
            end do
            ! reduce the EM Gamma accumulator BEFORE ncomp is replaced by d_new below.
            ! Gamma travels between parts as a SUM and is divided by the REDUCED nval: dividing per
            ! part would weight a small part equally with a large one.
            nval = sum(nval_thr)
            do q = 1, ncomp
                gam_sum(q) = sum(gam_thr(q,:))
            end do
            if( flex_pca_is_worker() )then
                allocate(cme(es(1),es(2),es(3),ncomp), cmo(es(1),es(2),es(3),ncomp))
                allocate(rhe(es(1),es(2),es(3),ncomp), rhoo(es(1),es(2),es(3),ncomp))
                do q = 1, ncomp
                    cme(:,:,:,q) = Yeven(q)%cmat_exp; rhe(:,:,:,q)  = Yeven(q)%rho_exp
                    cmo(:,:,:,q) = Yodd(q)%cmat_exp;  rhoo(:,:,:,q) = Yodd(q)%rho_exp
                end do
                call write_probe_part(flex_pca_part_fname('probe', params%part, params%numlen), &
                    &cme, rhe, cmo, rhoo, rho_e, rho_o, gam_sum, nval, ncomp)
                deallocate(cme, cmo, rhe, rhoo, gam_sum)
                call cleanup_rec_buffers(build, fpls)
                do q = 1, size(Yeven)
                    call Yeven(q)%dealloc_rho; call Yeven(q)%kill
                    call Yodd(q)%dealloc_rho;  call Yodd(q)%kill
                end do
                deallocate(Yeven, Yodd, rho_e, rho_o, prior)
                deallocate(Gth, Ath, bth, cth, zth, basis_fpls, mean_fpl, orientations, zbatch, dens)
                deallocate(valid, valid_e, valid_o, eo, Ainvth, Acpth, gam_thr, gam_acc, nval_thr)
                deallocate(z, ppinds)
                return
            endif
            endif
            do q = 1, ncomp
                gam_acc(q) = gam_sum(q) / real(max(1,nval),dp)
            end do
            deallocate(gam_sum)
            ! COUPLED M-step solve. At every grid point the components share ONE k x k normal matrix
            ! sum_i |CTF|^2 E[z_i z_i'], so the basis volumes have to be solved together. Dividing each
            ! Y_q by the scalar density sum_i |CTF|^2 instead drops the cross-component coupling and
            ! degrades the update to plain subspace iteration. Solve per halfset, then hand a unit
            ! density to compress_exp and skip sampl_dens_correct -- the divide has already happened.
            call solve_coupled_basis_exp(Yeven, rho_e, ncomp)
            call solve_coupled_basis_exp(Yodd,  rho_o, ncomp)
            ! finalize even/odd Y_q, half-set FSC Wiener-merge -> clean band-limited masked basis volume.
            allocate(realvols(ncomp))
            filtsz = max(1, fdim(params%box_crop) - 1)
            allocate(filt(filtsz), corrs(filtsz))
            do q = 1, ncomp
                ! even half -> band-limited real image (UNmasked, for an unbiased FSC)
                Yeven(q)%rho_exp = 1.0
                call Yeven(q)%compress_exp; call Yeven(q)%ifft
                call Yeven(q)%div(real(params%box))
                call realvols(q)%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
                call Yeven(q)%get_rmat_ptr(rmatp); call realvols(q)%set_rmat(rmatp, .false.)
                ! odd half
                Yodd(q)%rho_exp = 1.0
                call Yodd(q)%compress_exp; call Yodd(q)%ifft
                call Yodd(q)%div(real(params%box))
                call img_o%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
                call Yodd(q)%get_rmat_ptr(rmatp); call img_o%set_rmat(rmatp, .false.)
                ! half-set FSC -> Wiener filter 2F/(1+F)
                call realvols(q)%fft; call img_o%fft
                call realvols(q)%fsc(img_o, corrs)
                do sh = 1, filtsz
                    fc = max(0., min(0.999, corrs(sh)))
                    filt(sh) = 2.*fc/(1.+fc)
                end do
                ! merged, FSC-Wiener + band-limit filtered, back to real, masked
                call realvols(q)%add(img_o); call realvols(q)%mul(0.5)
                call realvols(q)%apply_filter(filt)
                if( params%lp > 2.0*params%smpd_crop + TINY ) call realvols(q)%bp(0., params%lp)
                call realvols(q)%ifft
                if( params%msk_crop > TINY ) call realvols(q)%mask3D_soft(params%msk_crop, backgr=0.)
                call img_o%kill
            end do
            deallocate(filt, corrs)
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
                if( cos_mean >= COV_PROBE_CONV ) l_converged = .true.
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
            end do
            call cleanup_rec_buffers(build, fpls)
            do q = 1, size(Yeven); call Yeven(q)%dealloc_rho; call Yeven(q)%kill; call Yodd(q)%dealloc_rho; call Yodd(q)%kill; end do
            do q = 1, size(utilde); call utilde(q)%dealloc_rho; call utilde(q)%kill; call utilde_real(q)%kill; end do
            do q = 1, size(realvols); call realvols(q)%kill; end do
            deallocate(Yeven, Yodd, utilde, utilde_real, realvols, rho_e, rho_o, prior)
            deallocate(Gth, Ath, bth, cth, zth, basis_fpls, mean_fpl, orientations, zbatch, dens, valid, valid_e, valid_o, eo)
            deallocate(Ainvth, Acpth, gam_thr, gam_acc, nval_thr)
            if( l_converged )then
                write(logfhandle,'(A,I0,A,F6.4,A)') '>>> FLEX_PCA PROBE converged after ',it, &
                    &' iterations (principal-angle cosine >= ',COV_PROBE_CONV,'); stopping early'
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
        deallocate(z, ppinds)
    end subroutine probe_subspace_iteration

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
        integer  :: i, j, k, attempt
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
        if( COV_UNIT_CONTRAST )then
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
        integer  :: pf, h, k, hmin, hmax, kmin, kmax, hp, kp, sh_lo, sh
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

    subroutine write_column_diagnostics( col_hkl, col_fsc, ncol )
        integer, intent(in) :: col_hkl(:,:), ncol
        real,    intent(in) :: col_fsc(ncol)
        integer :: u, s
        call del_file('flex_pca_selected_frequencies.txt')
        open(newunit=u,file='flex_pca_selected_frequencies.txt',status='replace',action='write')
        write(u,'(A)') '# column h k l halfset_column_fsc'
        do s = 1, ncol
            write(u,'(I5,3(1X,I5),1X,F10.6)') s, col_hkl(1,s), col_hkl(2,s), col_hkl(3,s), col_fsc(s)
        end do
        close(u)
        call del_file('flex_pca_column_fsc.txt')
        open(newunit=u,file='flex_pca_column_fsc.txt',status='replace',action='write')
        write(u,'(A)') '# column halfset_column_fsc'
        do s = 1, ncol
            write(u,'(I5,1X,F10.6)') s, col_fsc(s)
        end do
        close(u)
    end subroutine write_column_diagnostics

    subroutine cleanup_plane( fpl )
        type(fplane_type), intent(inout) :: fpl
        if( allocated(fpl%cmplx_plane)    ) deallocate(fpl%cmplx_plane)
        if( allocated(fpl%ctfsq_plane)    ) deallocate(fpl%ctfsq_plane)
        if( allocated(fpl%transfer_plane) ) deallocate(fpl%transfer_plane)
    end subroutine cleanup_plane

    ! ============================ SELF-CONTAINED TESTS ============================ Synthetic inputs only.

    !> svec packing must be an ISOMETRY: <X,Y>_Frobenius == <svec(X),svec(Y)>, and the index map a
    !! bijection onto 1..NP. The packed solve returns a wrong Sigma if either fails.
    subroutine test_flex_pca_svec_isometry()
        integer,  parameter :: D = 5, NP = D*(D+1)/2
        real(dp) :: X(D,D), Y(D,D), xs(NP), ys(NP), fro, pack_ip
        integer  :: i, j, k
        write(logfhandle,'(A)') '>>> TEST flex_pca svec packing isometry'
        do i = 1, D
            do j = 1, D
                X(i,j) = sin(real(i*D+j,dp))        ! deterministic
                Y(i,j) = cos(real(i+2*j,dp))
            end do
        end do
        X = 0.5d0*(X + transpose(X))                 ! symmetrize
        Y = 0.5d0*(Y + transpose(Y))
        do j = 1, D
            do i = 1, j
                k = svec_idx(i,j,D)
                xs(k) = merge(1.d0, sqrt(2.d0), i == j)*X(i,j)
                ys(k) = merge(1.d0, sqrt(2.d0), i == j)*Y(i,j)
            end do
        end do
        fro = sum(X*Y)
        pack_ip = sum(xs*ys)
        if( abs(fro - pack_ip) > 1.d-10*max(1.d0,abs(fro)) ) &
            &THROW_HARD('svec packing is not an isometry; the packed solve would return a wrong Sigma')
        ! index map must be a bijection onto 1..NP
        k = 0
        do j = 1, D
            do i = 1, j
                k = k + 1
                if( svec_idx(i,j,D) /= svec_idx(j,i,D) ) THROW_HARD('svec_idx is not symmetric in its arguments')
                if( svec_idx(i,j,D) < 1 .or. svec_idx(i,j,D) > NP ) THROW_HARD('svec_idx out of range')
            end do
        end do
        if( k /= NP ) THROW_HARD('svec index count does not match d(d+1)/2')
        write(logfhandle,'(A)') '>>>   PASSED (isometry to 1e-10, index map bijective)'
    end subroutine test_flex_pca_svec_isometry

    !> End-to-end recovery through the packed operator and CG.
    !! D >= 4 matters: with fewer columns no orbit has four distinct indices, so the generic
    !! three-cycle of rearrange_packed_selfsum would never be exercised.
    subroutine test_flex_pca_packed_solve()
        integer,  parameter :: D = 5, NPK = D*(D+1)/2, NG = 3
        real(dp) :: G(D,D,NG), Sig(D,D), Dmat(D,D), Ms(NPK,NPK), rhs(NPK), x(NPK), Sig_est(D,D)
        real(dp) :: sv(NPK,NG), tmp(D,D), err, relres
        real(dp) :: Aref(NPK,NPK), ei(NPK), col(NPK), err_rea
        integer  :: i, j, q, r, ig, k, l, iters
        write(logfhandle,'(A)') '>>> TEST flex_pca packed CG solve against a known Sigma'
        ! symmetric positive definite G_g = I*(g+1) + small symmetric perturbation
        G = 0.d0
        do ig = 1, NG
            do i = 1, D
                do j = 1, D
                    G(i,j,ig) = 0.1d0*sin(real(i*j + 3*ig,dp))
                end do
            end do
            G(:,:,ig) = 0.5d0*(G(:,:,ig) + transpose(G(:,:,ig)))
            do i = 1, D
                G(i,i,ig) = G(i,i,ig) + real(ig + 1,dp)    ! diagonally dominant -> PD
            end do
        end do
        ! known symmetric Sigma
        do i = 1, D
            do j = 1, D
                Sig(i,j) = 0.3d0*cos(real(2*i + j,dp))
            end do
        end do
        Sig = 0.5d0*(Sig + transpose(Sig))
        do i = 1, D
            Sig(i,i) = Sig(i,i) + 2.d0
        end do
        ! exact RHS D = sum_g G_g Sigma G_g
        Dmat = 0.d0
        do ig = 1, NG
            tmp  = matmul(G(:,:,ig), matmul(Sig, G(:,:,ig)))
            Dmat = Dmat + tmp
        end do
        ! packed accumulator, exactly as reduced_covariance_solve builds it
        do ig = 1, NG
            do j = 1, D
                do i = 1, j
                    sv(svec_idx(i,j,D),ig) = merge(1.d0, sqrt(2.d0), i == j)*G(i,j,ig)
                end do
            end do
        end do
        Ms = 0.d0
        do ig = 1, NG
            do l = 1, NPK
                do k = 1, NPK
                    Ms(k,l) = Ms(k,l) + sv(k,ig)*sv(l,ig)
                end do
            end do
        end do
        do r = 1, D
            do q = 1, r
                rhs(svec_idx(q,r,D)) = merge(1.d0, sqrt(2.d0), q == r)*Dmat(q,r)
            end do
        end do
        ! the rearrangement must reproduce apply_packed_A entry for entry -- build the reference
        ! operator column by column from the gathered routine, then rearrange in place and compare
        do l = 1, NPK
            ei = 0.d0; ei(l) = 1.d0
            call apply_packed_A(Ms, NPK, D, ei, col)
            Aref(:,l) = col
        end do
        call rearrange_packed_selfsum(Ms, NPK, D)
        err_rea = maxval(abs(Ms - Aref))/max(maxval(abs(Aref)), DTINY)
        write(logfhandle,'(A,ES10.3)') '>>>   rearranged vs gathered operator, max relative error=',err_rea
        if( err_rea > 1.d-12 ) THROW_HARD('packed rearrangement does not reproduce the packed operator')
        call cg_solve_packed(Ms, NPK, rhs, 0.d0, COV_CG_MAXIT, COV_CG_TOL, x, iters, relres)
        if( relres > 1.d-8 ) THROW_HARD('packed CG did not converge on a well-conditioned test system')
        do r = 1, D
            do q = 1, r
                Sig_est(q,r) = x(svec_idx(q,r,D)) / merge(1.d0, sqrt(2.d0), q == r)
                Sig_est(r,q) = Sig_est(q,r)
            end do
        end do
        err = maxval(abs(Sig_est - Sig)) / maxval(abs(Sig))
        write(logfhandle,'(A,I0,A,ES10.3,A,ES10.3)') '>>>   CG iterations=',iters, &
            &'  relative residual=',relres,'  max relative Sigma error=',err
        if( err > 1.d-6 ) THROW_HARD('packed CG solve did not recover the known Sigma')
        write(logfhandle,'(A)') '>>>   PASSED (Sigma recovered to 1e-6)'
    end subroutine test_flex_pca_packed_solve
end module simple_flex_pca_columns
