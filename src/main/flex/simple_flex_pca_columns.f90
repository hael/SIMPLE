! @descr:
module simple_flex_pca_columns
use simple_core_module_api
use simple_builder,          only: builder
use simple_image,            only: image
use simple_parameters,       only: parameters
use simple_reconstructor,    only: reconstructor
use simple_kbinterpol,       only: kbinterpol
use simple_gridding,         only: prep3D_inv_instrfun4mul
use simple_linalg,           only: jacobi, eigsrt, svd_solve
use simple_math,             only: ceil_div, floor_div
use simple_srch_sort_loc,    only: hpsort
use simple_matcher_3Drec,    only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io,  only: discrete_read_imgbatch, prepimgbatch
use simple_flex_reconstructor_latent_ops, only: project_fplane_mean, project_fplanes_mean_basis, &
    &insert_planes_oversamp_coupled_batch_scaled
use simple_flex_projected_latent_model,   only: prep_imgs4projected_model
implicit none
private
#include "simple_local_flags.inc"

public :: build_covariance_eigenbasis, embed_latents_with_contrast, estimate_covariance_mean
public :: svec_idx, apply_packed_A, cg_solve_packed
public :: probe_subspace_iteration, align_basis_to_reference, probe_external_basis
public :: cov_env_int_pub
public :: test_flex_pca_svec_isometry, test_flex_pca_packed_solve

! Density observability floor, matching simple_image_arith::div_cmat_at_1 and the projected-latent coupled
! solve.
real(dp), parameter :: COV_DENSITY_FLOOR = 1.0d-6
! Relative ridge used ONLY for the covariance diagonal in the S.B SNR proxy, which runs before the S.C
! weights exist (Algorithm 1 precedes Algorithm 2).
real,     parameter :: COV_RIDGE_REL     = 5.0e-2
! Relative eigenvalue floor for retaining direct-column PCA components.
real(dp), parameter :: COV_EIG_REL_FLOOR = 1.0d-6
! Cap on the column-subspace dimension. The accumulation is now a batched dsyrk on the Van Loan-Pitsianis
! rearrangement (see unrearrange_kron_selfsum), which needs ONE shared d^4 array regardless of thread count,
! so the budget below buys d ~ nthr^(1/4) times more dimension -- 53 -> 126 at nthr=30 -- and the
! accumulation runs at BLAS-3 rather than scalar-loop speed.
integer,  parameter :: COV_MAX_DTILDE    = 320
! Memory budget for the shared A accumulator, in bytes.
real(dp), parameter :: COV_ATHR_BUDGET   = 8.0d9
! When .true. the per-particle contrast subtraction was helpful only while the generalized-FSC ridge was
! inverted (fix 2.0). With the ridge corrected the 2.8 re-test shows a==1 is strictly better on IgG-100k --
! subspace capture of the GT hinge 0.395 (a=a_i) -> 0.420 (a==1) -- because subtracting a_i*T*mu also
! deletes the component of the conformational signal parallel to T*mu (proposal 2.8, point 1).
logical,  parameter :: COV_UNIT_CONTRAST  = .true.
! contrast mean 2.000, sd 0.000 over all 100k particles, i.e. Earlier IgG-20k measurement also showed no
! z-tau lift (a==1 0.045 vs grid 0.047).
logical,  parameter :: COV_EMBED_CONTRAST_GRID = .false.
integer,  parameter :: COV_CG_MAXIT = 2000     ! CG iteration cap; convergence is reported, not assumed
real(dp), parameter :: COV_CG_TOL   = 1.d-10   ! relative residual target
integer,  parameter :: GRAM_DIAG_STRIDE = 200   ! subsample for the projected-Gram spectrum
integer,  parameter :: NCONTRAST_GRID = 50
real(dp), parameter :: A_GRID_HI = 2.0d0
real(dp), parameter :: COV_PINV_RCOND = 1.0d-6

! Source of the covariance mean mu.
logical, parameter :: COV_MEAN_FROM_DATA = .false.

! At box_crop=64 with mskdiam=200 the disc keeps ~21 % of the frame before the margin below, i.e. ~42 %
! after it.
logical, parameter :: COV_MASK_IMAGES = .false.
! Radial margin on that disc, as a multiple of the model radius. delocalisation is lambda*defocus/d, which
! on this benchmark reaches ~70 A at 5.5 um defocus and 15 A resolution, against a model radius of 100 A.
! 1.4 keeps the molecule plus that margin and still discards ~60 % of the frame.
real, parameter :: COV_MASK_MARGIN = 1.4

! Subtract the analytic per-sample noise bias K_R(.,q_s)|T|^2 from the column numerator. The per-shell
! diagnostic shows the half-set column FSC collapsing there (0.00 / 0.00 / 0.01 / 0.13 / 0.35 at shells 0-4,
! against 0.99), which drives the Wiener shrinkage <H>/(<H>+R) down to 0.001-0.13 and deletes the band.
logical, parameter :: COV_COLUMN_NOISE_DEBIAS = .true.

! Width of the RIGHT kernel -- the one that reads each image's value AT the column frequency
! (gather_column_values). SIMPLE historically used ONE shared 3-tap KB stencil for both, whose support is
! |d| <= 1.5 per axis against |d| < 2, so it gathers roughly half as many image samples into each column.
real :: COV_RIGHT_KERNEL_W = 0.0
logical :: cov_rkw_read = .false.

contains

    !> HeteroPCA (Zhang, Cai & Wu, Ann. the top latent component carries 53% of the latent variance but its
    !! between-conformation SNR is 0.007 -- it is a pure nuisance mode, and its eigenvalue is 4.7x the next.
    !! the analytic debias Dmat = Sbb - 0.5*sig2_eff*SG subtracts a modelled noise term, and whatever it
    !! fails to remove piles onto the diagonal, inflating the leading eigenvalue and tilting its eigenvector
    !! toward the lowest-frequency mode.
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

    !> a typical particle carried weight in 5.6 of 15 states.
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
        real(dp), allocatable :: vred(:,:)
        integer :: ncol, nreal, s, lb(3), ub(3), nyq_rec, d_tilde, q, directsvd
        real(dp), allocatable :: svals(:)
        ! one work reconstructor defines the expanded lattice / Nyquist / grid correction
        call init_column_reconstructor(params, build, work)
        lb      = lbound(work%cmat_exp)
        ub      = ubound(work%cmat_exp)
        nyq_rec = work%get_lfny(1)
        ! column selection:
        if( trim(params%column_sampling) == 'lowfreq' )then
            call select_covariance_columns_lowfreq(params, ncols_req, col_sep, col_hkl, ncol)
        else
            call select_covariance_columns_snr(params, build, mean_rec, pinds, nptcls, &
                &lb, ub, nyq_rec, ncols_req, col_sep, col_hkl, ncol)
        endif
        write(logfhandle,'(A,I0,A,A,A,I0)') '>>> FLEX_PCA selected covariance columns=',ncol, &
            &' sampling=',trim(params%column_sampling),' separation=',col_sep
        call flush(logfhandle)
        call build_column_lookup(col_hkl, ncol, lb, ub, col_lookup)
        allocate(Bcol_e(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol), source=cmplx(0.,0.))
        allocate(Bcol_o(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol), source=cmplx(0.,0.))
        allocate(Hcol_e(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol), source=0.)
        allocate(Hcol_o(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol), source=0.)
        call accumulate_covariance_columns(params, build, mean_rec, pinds, nptcls, &
            &col_hkl, col_lookup, ncol, lb, ub, nyq_rec, Bcol_e, Hcol_e, Bcol_o, Hcol_o)
        ! regularized merge + half-column FSC diagnostics
        allocate(colvol(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),ncol), col_fsc(ncol))
        call regularize_and_merge_columns(Bcol_e, Hcol_e, Bcol_o, Hcol_o, ncol, lb, colvol, col_fsc)
        call write_column_diagnostics(col_hkl, col_fsc, ncol)
        deallocate(Bcol_e, Bcol_o, Hcol_e, Hcol_o, col_lookup, col_hkl)
        ! Friedel-symmetrize each complex column into two real (cos/sin) volumes
        call columns_to_real_representatives(params, work, colvol, ncol, lb, ub, realvols, nreal)
        deallocate(colvol)
        call work%dealloc_rho; call work%kill
        write(logfhandle,'(A,I0)') '>>> FLEX_PCA real column representatives=',nreal
        call flush(logfhandle)
        ! NOTE:
        call orthonormalize_representatives(params, build, realvols, nreal, utilde, utilde_real, d_tilde, svals)
        do s = 1, nreal
            call realvols(s)%kill
        end do
        deallocate(realvols, col_fsc)
        write(logfhandle,'(A,I0)') '>>> FLEX_PCA column subspace dimension d_tilde=',d_tilde
        call flush(logfhandle)
        ! S.B reduced projected-covariance solve (eqs S.6-S.9):
        call reduced_covariance_solve(params, build, mean_rec, utilde, d_tilde, pinds, nptcls, &
            &neigs_req, vred, eigvals, ncomp_out, sig2_out)
        ! DIRECT-SVD BASIS (SIMPLE_COV_DIRECTSVD=1). SIMPLE instead uses the columns only to define a SPAN
        ! and re-estimates the covariance inside it, which discards that weighting -- and measured, its
        ! projected-Gram spectrum is more top-heavy (lam1/lam5 2.24 vs 1.70, effective rank 9.35 vs 10.52)
        ! and it resolves 1 component per particle against 7.
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
        call form_eigenbasis_from_reduced(params, build, mean_rec, utilde_real, d_tilde, vred, eigvals, &
            &ncomp_out, basis_recs, basis_imgs=basis_imgs, fprefix=fprefix)
        do s = 1, d_tilde
            call utilde(s)%dealloc_rho; call utilde(s)%kill
            call utilde_real(s)%kill
        end do
        deallocate(utilde, utilde_real, vred)
    end subroutine build_covariance_eigenbasis

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
        ! Report the radial placement of the selected columns. greedy choice on this benchmark lands
        ! entirely inside |xi| <= 8 of 64 shells ("Largest frequency computed:
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

    !> Per-voxel signal-variance (SNR proxy) from the whitened adjoint residuals:
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
        integer(timer_int_kind) :: t_phase
        pf = OSMPL_PAD_FAC; inv_pf2 = 1.0/real(pf*pf)
        allocate(var_acc(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), source=0.)
        allocate(dens_acc(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), source=0.)
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
                var_acc  = var_acc  + (real(scratch%cmat_exp)*inv_pf2)**2 + (aimag(scratch%cmat_exp)*inv_pf2)**2
                dens_acc = dens_acc + scratch%rho_exp
            end do
        end do
        write(logfhandle,'(A,F8.1)') '>>> FLEX_PCA SNR VARIANCE SECONDS: ', toc(t_phase)
        call flush(logfhandle)
        call orientation%kill
        call cleanup_plane(mean_fpl); call cleanup_plane(adj_fpl)
        call scratch%dealloc_rho; call scratch%kill
        call cleanup_rec_buffers(build, fpls)
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

    !> Deterministic low-frequency column selection:
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

    !> EXTERNAL-BASIS PROBE. passing the known eigenvolumes through it gave rho = 0.30 where this machinery
    !! gives 0.801 on the same particles.
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
        ! nprobe = 0 is legitimate and important: With a probe volume appended, its projected norm can
        ! dominate the Gram spectrum and the reported conditioning then describes the probe rather than the
        ! basis -- measured, an appended volume contributed a leading eigenvalue 60x the next.
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
            call estimate_mean_scale(params, build, mean_rec, pinds, nptcls)
        endif
    end subroutine estimate_covariance_mean

    !> Kernel-regression mean of eq. on IgG-RL/100k the self-estimated mean moves the SNR-greedy column
    !! choice from |xi| mean 6.2 / max 8.2 (which is where greedy choice lands, "Largest frequency computed:
    !! 8") out to |xi| mean 9.4 / max 14.5, and the resulting eigenbasis captures 0.363 of the band-limited
    !! ground-truth motion subspace against 0.596 for the vol1 path.
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

    !> Self-estimate the global amplitude scale of the consensus mean map relative to the whitened data
    !! (proposal task:
    subroutine estimate_mean_scale( params, build, mean_rec, pinds, nptcls )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: pinds(:), nptcls
        integer, parameter :: NSAMPLE = 4000
        type(fplane_type), allocatable :: fpls(:)
        type(fplane_type) :: mean_fpl
        type(ori) :: orientation
        integer  :: batchlims(2), batchsz, ibatch, i, iptcl, stride, used, nyq, sh
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
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz), mskrad=cov_image_mask_radius(params))
            do i = 1, batchsz
                if( mod(batchlims(1)+i-2, stride) /= 0 ) cycle
                iptcl = pinds(batchlims(1)+i-1)
                call build%spproj_field%get_ori(iptcl, orientation)
                if( orientation%isstatezero() ) cycle
                call project_fplane_mean(mean_rec, orientation, fpls(i), mean_fpl, apply_ctf_amp=.true.)
                s_my = s_my + real(cov_herm_inner(mean_fpl, fpls(i)), dp)
                s_mm = s_mm + real(cov_herm_inner(mean_fpl, mean_fpl), dp)
                call plane_shell_cross_accum(mean_fpl, fpls(i), nyq, smy_sh, smm_sh)
                used = used + 1
            end do
        end do
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
        ! Per-shell mean/data amplitude scale (fix 2.2, diagnostic D5).
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
        deallocate(smy_sh, smm_sh, sprof, filt)
    end subroutine estimate_mean_scale
    !>  Accumulate per-shell mean/data cross power Re<T mu, y> and mean auto power |T mu|^2 over the
    !!  native k<=0 half. The per-shell ratio s(sh)=sum my_sh/sum mm_sh is the ML mean amplitude scale
    !!  at each shell (fix 2.2); the k=0 double-count cancels in the ratio so no weighting is needed.
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
    !> Reconstruction-mode plane from a whitened observation-model plane:
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
        integer :: h_sq, k_max_h, k_lo, k_hi, win(2,3), plb(2), pub(2)
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
                call cov_kb_weights(kbwin, loc, win, wx, wy, wz)
                if( any(win(1,:) < lb) .or. any(win(2,:) > ub) ) cycle
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
            ! reference-style independent right kernel:
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
        integer :: s, j, ix, iy, iz, hx, ky, mz, qh, qk, ql, iqx, iqy, iqz
        real    :: pf2, w, wqs, wyz, hs
        complex :: cg, cnum
        logical :: in_win
        pf2 = real(OSMPL_PAD_FAC*OSMPL_PAD_FAC)
        !$omp parallel do default(shared) schedule(static) proc_bind(close) &
        !$omp& private(s,j,ix,iy,iz,hx,ky,mz,qh,qk,ql,iqx,iqy,iqz,w,wqs,wyz,hs,cg,cnum,in_win)
        do s = 1, ncol
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

    !> KB separable stencil, replicating simple_flex_reconstructor_latent_ops::
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
        ! ---- Regularizer initialization. Taking w as the mean map's per-shell Fourier power and sigma^2 =
        ! 1 (the supplement's whitened convention) was tried and MEASURED to break the stage -- v^2 lands
        ! below the 1e-6 density floor, so P pins at the floor, R pins at ~1e6, the per-shell |Sigma|^2 sums
        ! underflow the DTINY guard, and every column FSC comes out identically 0.0000 instead of a profile.
        ! SIMPLE's non-unitary FFT/gridding convention means sigma^2 is ~3.9e-6 here, not 1, and the mean
        ! map carries yet another amplitude convention, so v v^T cannot be formed without definition of w.
        do s = 1, ncol
            hbar_sum = 0.d0; num = 0.d0; nvox_sh = 0.d0
            do i3 = 1, n3; do i2 = 1, n2; do i1 = 1, n1
                sh = shell(i1,i2,i3)
                hbar_sum(sh) = hbar_sum(sh) + 0.5d0*(real(Hcol_e(i1,i2,i3,s),dp)+real(Hcol_o(i1,i2,i3,s),dp))
                num(sh) = num(sh) + 1.d0
                nvox_sh(sh) = nvox_sh(sh) + 1.d0
            end do; end do; end do
            ! Initial ridge. the per-shell Wiener factor was 0.001-0.13 at shells 0-3 while the merged
            ! column magnitude at those shells was 1e-6 of its shell-13 value.
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
                    ! `den > DTINY` (1e-10) on a sum of |B/(H+R)|^2, a quantity carrying the pipeline's
                    ! non-unitary Fourier convention (sig2_eff ~ 4e-6 here, not 1) and evaluated where H is
                    ! LARGEST -- the lowest shells. Because the ridge is R = 1/P with P proportional to
                    ! FSC/(1-FSC), that zero drove the Wiener factor <H>/(<H>+R) to 0.001-0.13 and DELETED
                    ! the conformational band from every column.
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
            ! is REMOVED:
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
                ! Rsh, not 2*Rsh.
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
        ! operation. at box_crop=64 with lp=15 (kstop=25) the out-of-band shells 26-27 come back at FSC
        ! 0.997, i.e.
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
        integer :: i, q, nrot, keep, d_budget, cg_budget
        real(dp) :: lam_max, nrm
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
        ! Largest d the reduced solve's accumulator can afford. one shared d^2 x d^2 array -> 8*d^4 bytes
        ! packed + CG : one shared npk x npk array, npk = d(d+1)/2, and the operator is never formed ->
        ! 8*[d(d+1)/2]^2 ~ 2*d^4 bytes The packed accumulator is 4x smaller, which is sqrt(2) in d -- at the
        ! 8 GB budget, d 177 -> 250, past rank cap of 200.
        cg_budget = 0
        call cov_env_int('SIMPLE_COV_CGSOLVE', cg_budget)
        if( cg_budget /= 0 )then
            ! d(d+1)/2 = sqrt(BUDGET/8)  =>  d = (-1 + sqrt(1 + 8*sqrt(BUDGET/8)))/2
            d_budget = max(1, int((-1.d0 + sqrt(1.d0 + 8.d0*sqrt(COV_ATHR_BUDGET/8.d0)))/2.d0))
        else
            d_budget = max(1, int((COV_ATHR_BUDGET / 8.d0)**0.25d0))
        endif
        ! SIMPLE_COV_DTILDE pins the subspace dimension, so a change to the accumulation or the
        ! solver can be A/B'd at a FIXED d against an earlier run instead of confounding the two.
        call cov_env_int('SIMPLE_COV_DTILDE', d_budget)
        d_tilde  = max(1, min(keep, COV_MAX_DTILDE, d_budget))
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA d_tilde=',d_tilde, &
            &'  (above energy floor=',keep,', memory cap=',d_budget,', rank cap=',COV_MAX_DTILDE
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
            ! set_rmat(...,.false.) -- NOT add() -- see form_eigenbasis_from_reduced:
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
    !! FORMING A_s.
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
    subroutine cg_solve_packed( Ms, np, d, rhs, ridge, maxit, tol, x, iters, relres )
        integer,  intent(in)    :: np, d, maxit
        real(dp), intent(in)    :: Ms(np,np), rhs(np), ridge, tol
        real(dp), intent(out)   :: x(np)
        integer,  intent(out)   :: iters
        real(dp), intent(out)   :: relres
        real(dp), allocatable :: r(:), p(:), Ap(:), z(:), diag(:)
        real(dp) :: rz, rz_new, alpha, beta_cg, pAp, nrm0, dg
        integer  :: it, k, al, be
        allocate(r(np), p(np), Ap(np), z(np), diag(np))
        ! Jacobi preconditioner: the diagonal of A_s, read straight out of Ms
        do k = 1, np
            be = int((sqrt(8.d0*real(k,dp)-7.d0)+1.d0)/2.d0)
            if( svec_idx(1,be,d) > k ) be = be - 1
            al = k - (be-1)*be/2
            if( al == be )then
                dg = packed_T(Ms, np, al, al, al, al, d)
            else
                dg = packed_T(Ms, np, al, al, be, be, d) + packed_T(Ms, np, al, be, be, al, d)
            endif
            diag(k) = max(dg + ridge, DTINY)
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
            call apply_packed_A(Ms, np, d, p, Ap)
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
        real(dp), allocatable :: hfpw_thr(:), hfcnt_thr(:)
        real(dp), allocatable :: cvpw_thr(:), cvcnt_thr(:)   ! DIAGNOSTIC (see cov_herm_selfpower)
        real(dp), allocatable :: babs_thr(:), bre_thr(:), gdi_thr(:)   ! DIAGNOSTIC: |b|^2 vs Re(b)^2 vs G_qq
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
        ! PACKED accumulation (SIMPLE_COV_CGSOLVE=1).
        cgsolve = 0
        call cov_env_int('SIMPLE_COV_CGSOLVE', cgsolve)
        npk = d_tilde*(d_tilde+1)/2
        if( cgsolve /= 0 )then
            allocate(Mspk(npk,npk), source=0.d0)
            write(logfhandle,'(A,I0,A,F8.3,A,F8.3,A)') '>>> FLEX_PCA PACKED+CG solve: npk=',npk, &
                &'  accumulator ', 8.d0*real(npk,dp)**2/1.d9, ' GB (dense would be ', &
                &8.d0*real(dd,dp)**2/1.d9, ' GB)'
            call flush(logfhandle)
        else
            allocate(A(dd,dd), source=0.d0)
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
        allocate(hfpw_thr(nthr), hfcnt_thr(nthr), source=0.d0)
        allocate(cvpw_thr(nthr), cvcnt_thr(nthr), source=0.d0)   ! DIAGNOSTIC: noise in cov_herm_inner's index convention
        allocate(babs_thr(nthr), bre_thr(nthr), gdi_thr(nthr), source=0.d0)   ! DIAGNOSTIC
        nyq_rec = utilde(1)%get_lfny(1)
        allocate(pwsh_thr(0:nyq_rec,nthr), cntsh_thr(0:nyq_rec,nthr), source=0.d0)
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
            !$omp& private(i,ithr,q,r,alpha,beta,gam1,delta,a1,a2,pw,cnt)
            do i = 1, batchsz
                if( orientations(i)%isstatezero() ) cycle
                ithr = omp_get_thread_num() + 1
                call project_fplanes_mean_basis(mean_rec, utilde, orientations(i), fpls(i), &
                    &mean_fpl(ithr), basis_fpls(:,ithr), apply_ctf_amp=.true.)
                call form_residual_plane(fpls(i), mean_fpl(ithr), resid_fpl(ithr), &
                    &particle_contrast(mean_fpl(ithr), fpls(i)))
                ! signal-free noise variance from the outer residual shells
                call plane_hf_power(resid_fpl(ithr), nyq_rec, 0.7, pw, cnt)
                hfpw_thr(ithr)  = hfpw_thr(ithr)  + pw
                hfcnt_thr(ithr) = hfcnt_thr(ithr) + cnt
                ! DIAGNOSTIC: same residual, measured in cov_herm_inner's index convention
                call cov_herm_selfpower(resid_fpl(ithr), pw, cnt)
                cvpw_thr(ithr)  = cvpw_thr(ithr)  + pw
                cvcnt_thr(ithr) = cvcnt_thr(ithr) + cnt
                call plane_shell_power_accum(resid_fpl(ithr), nyq_rec, pwsh_thr(:,ithr), cntsh_thr(:,ithr))
                do q = 1, d_tilde
                    bc_thr(q,ithr) = cov_herm_inner(basis_fpls(q,ithr), resid_fpl(ithr))
                    ! DIAGNOSTIC: is Re(b) anomalously small vs |b|, and how big is G_qq?
                    babs_thr(ithr) = babs_thr(ithr) + real(bc_thr(q,ithr)*conjg(bc_thr(q,ithr)),dp)
                    bre_thr(ithr)  = bre_thr(ithr)  + real(bc_thr(q,ithr))**2
                    do r = q, d_tilde
                        Gi_thr(q,r,ithr) = real(cov_herm_inner(basis_fpls(q,ithr), basis_fpls(r,ithr)), dp)
                        Gi_thr(r,q,ithr) = Gi_thr(q,r,ithr)
                    end do
                    gdi_thr(ithr) = gdi_thr(ithr) + Gi_thr(q,q,ithr)
                end do
                do q = 1, d_tilde
                    do r = 1, d_tilde
                        ! properness fix (2.1): E[Re(b_q)Re(b_r)]_noise = 0.5*sig2*G (see Dmat).
                        Sbb_thr(q,r,ithr) = Sbb_thr(q,r,ithr) + real(bc_thr(q,ithr))*real(bc_thr(r,ithr))
                        SG_thr(q,r,ithr)  = SG_thr(q,r,ithr) + Gi_thr(q,r,ithr)
                    end do
                end do
                ! Stage vec(G_i) for the batched rank-k update below.
                if( cgsolve /= 0 )then
                    ! svec:
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
            if( cgsolve /= 0 )then
                call dsyrk('U', 'N', npk, batchsz, 1.d0, Vg, npk, 1.d0, Mspk, npk)
            else
                call dsyrk('U', 'N', dd, batchsz, 1.d0, Vg, dd, 1.d0, A, dd)
            endif
            Vg(:,1:batchsz) = 0.d0
        end do
        ! dsyrk touched only the upper triangle; mirror it before the un-rearrangement.
        if( cgsolve /= 0 )then
            !$omp parallel do default(shared) private(a1,a2) schedule(static)
            do a2 = 1, npk
                do a1 = a2+1, npk
                    Mspk(a1,a2) = Mspk(a2,a1)
                end do
            end do
            !$omp end parallel do
            ! NO un-rearrangement: apply_packed_A reads the rearranged packed accumulator directly
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
        if( nvalid < 1 ) THROW_HARD('flex_pca reduced covariance solve found no valid particles')
        ! Noise-bias scale sig2_eff = signal-free per-coefficient whitened-noise variance E|n|^2,
        ! measured from the outer residual shells. With the properness-consistent RHS the debias is
        ! D = Sbb - 0.5*sig2_eff*SG (see below), which removes the noise bias WITHOUT over-subtracting
        ! the signal (unlike trace matching, which folds the signal into sig2_eff and collapses the rank).
        sig2_eff = sum(hfpw_thr) / max(sum(hfcnt_thr), 1.d0)
        ! per-shell whitened residual variance profile (DIAGNOSTIC):
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
        ! than a rounding residue. with a dead (identically zero) projected basis tr_g underflowed the
        ! DTINY=1e-10 floor (measured tr_g=6.05e-18), so the "1.66e-13" it reported was simply tr_bb/1e-10
        ! -- a pure floor artefact that read as a plausible noise scale and cost hours.
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
        ! DIAGNOSTIC: The debias identity E[Re(b)Re(b)] = 0.5*sig2*G requires sig2 in THIS convention
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
        ! Eq. SIMPLE's FFT/gridding convention is not unitary, so after whitening the noise covariance is
        ! Lambda_i = sig2_conv * I with sig2_conv a pure convention constant (box, padding, gen_fplane4rec
        ! scaling) -- measured at ~3.9e-6 at box 128, not 1. The 0.5 IS exact:
        Dmat = Sbb - 0.5d0 * sig2_eff * SG
        ! Same identity gives the embedding's noise scale: b = G z + n with Cov(n) = 0.5*sig2*G, so eq.
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
        ! squares of symmetric matrices is PSD
        if( cgsolve /= 0 )then
            ! pack the RHS and solve matrix-free
            allocate(rhspk(npk), spk(npk))
            do r = 1, d_tilde
                do q = 1, r
                    rhspk(svec_idx(q,r,d_tilde)) = merge(1.d0, sqrt(2.d0), q == r)*Dmat(q,r)
                end do
            end do
            call cg_solve_packed(Mspk, npk, d_tilde, rhspk, ridge, COV_CG_MAXIT, COV_CG_TOL, &
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
        ! Higham PSD projection:
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
            ! reliable direct write of the real-space eigenvolume
            if( present(fprefix) )then
                fname = trim(fprefix)//int2str_pad(k,3)//MRC_EXT
            else
                fname = 'flex_pca_pc'//int2str_pad(k,3)//MRC_EXT
            endif
            call eigimg%write(fname, del_if_exists=.true.)
            call fname%kill
            if( present(basis_imgs) ) call basis_imgs(k)%copy(eigimg)
            ! reconstructor for embedding projection (mean_rec idiom). add() then deposits into the shared
            ! rmat/cmat buffer while the subsequent fft() is a no-op, so expand_exp propagates an
            ! untransformed grid and the projected basis comes out IDENTICALLY ZERO for ~99% of particles
            ! (measured:
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
        &pinds, nptcls, z, contrast, precision, resid_energy, resid_mean_energy )
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
        real(dp), parameter :: A_LO = 0.1d0, A_HI = 5.0d0
        type(fplane_type), allocatable :: fpls(:)
        type(fplane_type), allocatable :: basis_fpls(:,:), mean_fpl(:), data_fpl(:)
        type(ori), allocatable :: orientations(:)
        real(dp), allocatable :: prior(:), rho(:), Gcache(:,:,:), bcache(:,:), ccache(:,:)
        real(dp), allocatable :: zhalf(:,:,:), Ghf(:,:,:,:), bhf(:,:,:), chf(:,:,:)
        real(dp), parameter   :: RHO_FLOOR = 1.d-3
        real(dp) :: rho_max, rrel
        logical :: l_relprior
        integer :: ihf
        integer :: batchlims(2), batchsz, ibatch, i, iptcl, q, r, ithr, nthr, ia, row
        integer, allocatable :: nzeroG_thr(:), nzeroR_thr(:), nzeroZ_thr(:)   ! DIAGNOSTIC D2
        real(dp) :: a, a_best, e_yy, e_mm, best_res, res, aa, sig2
        integer(timer_int_kind) :: t_phase
        real(dp), allocatable :: Gth(:,:,:), Ath(:,:,:), zth(:,:), zbest(:,:), cth(:,:), bth(:,:), myth(:)
        real(dp), allocatable :: Gtilth(:,:,:)   ! per-thread noise-whitened projected Gram
        real(dp), allocatable :: gwork(:,:,:), gvec(:,:,:), gev(:,:), gspec_thr(:,:), gspec(:)
        integer,  allocatable :: gcnt_thr(:)
        integer :: nrot_t, gcnt
        real(dp) :: gsum
        nthr = nthr_glob
        sig2 = max(sig2_eff, DTINY)      ! whitened-noise variance for the MAP shrinkage (fix 2.4)
        allocate(nzeroG_thr(nthr), nzeroR_thr(nthr), nzeroZ_thr(nthr), source=0)   ! DIAGNOSTIC D2
        l_relprior = .not. cov_env_int_off('SIMPLE_COV_RELPRIOR')
        allocate(prior(ncomp))
        if( l_relprior )then
            allocate(Gcache(ncomp,ncomp,nptcls), bcache(ncomp,nptcls), ccache(ncomp,nptcls), source=0.d0)
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
                ! SPLIT-HALF:
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
                ! Contrast: For each a on linspace(0,2,51)[1:] solve the fixed-a MAP (a^2 G/sig2 + Gamma^-1)
                ! z = (a b - a^2 c)/sig2 (fix 2.4, S.E
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
                ! PROJECTED-GRAM CONDITIONING.
                if( mod(row, GRAM_DIAG_STRIDE) == 0 )then
                    gwork(:,:,ithr) = Gth(:,:,ithr)
                    call jacobi(gwork(:,:,ithr), ncomp, ncomp, gev(:,ithr), gvec(:,:,ithr), nrot_t)
                    call eigsrt(gev(:,ithr), gvec(:,:,ithr), ncomp, ncomp)
                    do q = 1, ncomp
                        gspec_thr(q,ithr) = gspec_thr(q,ithr) + max(gev(q,ithr), 0.d0)
                    end do
                    gcnt_thr(ithr) = gcnt_thr(ithr) + 1
                endif
                ! DIAGNOSTIC D2: is the projected basis / rhs numerically dead for this particle?
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
            ! NORMALISE FIRST.
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
        ! ---- RELIABILITY-WEIGHTED PRIOR ---- prior_precision = sig2/Gamma_q hands the LARGEST eigenvalue
        ! the WEAKEST prior. component 1 has Gamma 4.66x the next and the smallest prior_precision of all
        ! 20, so the MAP barely constrains it and it becomes near-unregularized least squares along a
        ! direction the data hardly measures -- 53% of the latent variance at between-conformation SNR
        ! 0.007, correlating with nothing (defocus 0.001, per-particle scale 0.037).
        if( l_relprior )then
            allocate(rho(ncomp))
            do q = 1, ncomp
                rho(q) = corr_dp(zhalf(:,q,1), zhalf(:,q,2), nptcls)
                rho(q) = max(0.d0, rho(q))
                rho(q) = 2.d0*rho(q) / (1.d0 + rho(q))            ! Spearman-Brown to full length
            end do
            ! RELATIVE, not absolute. it fixed the nuisance mode (z1 left the top-8 by variance,
            ! conformational share 0.352 -> 0.426) but over-shrank the informative components too (z2's
            ! prior +34%), compressing the latent spread the trajectory depends on and costing swing
            ! coverage 74.5% -> 64.0% and trajectory amplitude 0.726 -> 0.376.
            rho_max = maxval(rho)
            if( rho_max <= DTINY ) rho_max = 1.d0
            do q = 1, ncomp
                rrel = (rho(q)*rho(q)) / (rho_max*rho_max)
                prior(q) = 1.d0 / max(max(rrel, RHO_FLOOR) * eigvals(q), DTINY)
            end do
            write(logfhandle,'(A)') '>>> FLEX_PCA split-half reliability per component (rho, corrected):'
            do q = 1, min(ncomp,10)
                write(logfhandle,'(A,I3,A,F7.4,A,ES11.3,A,ES11.3)') '>>>   z',q,'  rho=',rho(q), &
                    &'  eigval=',eigvals(q),'  prior_precision=',prior(q)
            end do
            call flush(logfhandle)
            ! re-solve every particle in closed form from the cached sufficient statistics
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
            write(logfhandle,'(A)') '>>> FLEX_PCA latents re-solved with the reliability-weighted prior'
            call flush(logfhandle)
            deallocate(rho, Gcache, bcache, ccache, zhalf, Ghf, bhf, chf)
        endif
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
    end subroutine embed_latents_with_contrast

    !> Native probe-based subspace iteration (proposal 4.1).
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
        real(dp),            allocatable :: zbatch(:,:), dens_dummy(:,:,:), z(:,:), prior(:)
        real(dp),            allocatable :: Gth(:,:,:), Ath(:,:,:), bth(:,:), cth(:,:), zth(:,:)
        logical,             allocatable :: valid(:), valid_e(:), valid_o(:)
        integer,             allocatable :: eo(:)
        real, pointer :: rmatp(:,:,:)
        real     :: fc
        real(dp) :: sig2, a, aa, e_mm, myv, mu_q, sd_q
        integer  :: it, q, r, i, ithr, nthr, batchlims(2), batchsz, ibatch, row, d_new, es(3), filtsz, sh
        type(string) :: fname
        integer(timer_int_kind) :: t_it
        nthr = nthr_glob
        sig2 = max(sig2_eff, DTINY)
        allocate(z(nptcls,ncomp))
        do it = 1, niters
            t_it = tic()
            write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA PROBE SUBSPACE ITERATION ',it,' / ',niters, &
                &'  basis dim=',ncomp
            call flush(logfhandle)
            allocate(prior(ncomp))
            do q = 1, ncomp
                prior(q) = 1.d0 / max(eigvals(q), DTINY)
            end do
            ! even/odd Y_q accumulators (half-set FSC regularization) + shared sampling density
            allocate(Yeven(ncomp), Yodd(ncomp))
            do q = 1, ncomp
                call init_column_reconstructor(params, build, Yeven(q)); call Yeven(q)%reset; call Yeven(q)%reset_exp
                call init_column_reconstructor(params, build, Yodd(q));  call Yodd(q)%reset;  call Yodd(q)%reset_exp
            end do
            es = shape(Yeven(1)%cmat_exp)
            allocate(rho_e(1,es(1),es(2),es(3)), rho_o(1,es(1),es(2),es(3)), source=0.)
            allocate(Gth(ncomp,ncomp,nthr), Ath(ncomp,ncomp,nthr), bth(ncomp,nthr), cth(ncomp,nthr), zth(ncomp,nthr))
            allocate(basis_fpls(ncomp,nthr), mean_fpl(nthr), orientations(MAXIMGBATCHSZ))
            allocate(zbatch(ncomp,MAXIMGBATCHSZ), dens_dummy(ncomp,ncomp,MAXIMGBATCHSZ))
            allocate(valid(MAXIMGBATCHSZ), valid_e(MAXIMGBATCHSZ), valid_o(MAXIMGBATCHSZ), eo(MAXIMGBATCHSZ))
            dens_dummy = 0.d0
            call mean_rec%expand_exp
            do q = 1, ncomp
                call basis_recs(q)%expand_exp
            end do
            call init_rec(params, build, MAXIMGBATCHSZ, fpls)
            call prepimgbatch(params, build, MAXIMGBATCHSZ)
            z = 0.d0
            do ibatch = 1, nptcls, MAXIMGBATCHSZ
                batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
                batchsz   = batchlims(2) - batchlims(1) + 1
                call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
                call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                    &pinds(batchlims(1):batchlims(2)), fpls(:batchsz), mskrad=cov_image_mask_radius(params))
                do i = 1, batchsz
                    call build%spproj_field%get_ori(pinds(batchlims(1)+i-1), orientations(i))
                    eo(i) = build%spproj_field%get_eo(pinds(batchlims(1)+i-1))
                end do
                valid(:batchsz) = .false.
                zbatch(:,:batchsz) = 0.d0
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
                    Ath(:,:,ithr) = (aa/sig2)*Gth(:,:,ithr)
                    do q = 1, ncomp
                        Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                        zth(q,ithr)   = (a*bth(q,ithr) - aa*cth(q,ithr))/sig2
                    end do
                    call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), ncomp)
                    z(row,:)          = zth(:,ithr)
                    zbatch(:,i)       = zth(:,ithr)
                    valid(i)          = .true.
                    ! residual observation r_i = y - a*(T mu) in place (transfer/ctfsq intact for backprojection)
                    fpls(i)%cmplx_plane = fpls(i)%cmplx_plane - real(a)*mean_fpl(ithr)%cmplx_plane
                end do
                !$omp end parallel do
                ! M-step by halfset: Y_q += sum_i z_iq * backproject(r_i)   (shared density, batched KB)
                do i = 1, batchsz
                    valid_e(i) = valid(i) .and. eo(i)==0
                    valid_o(i) = valid(i) .and. eo(i)==1
                end do
                call insert_planes_oversamp_coupled_batch_scaled(Yeven, rho_e, build%pgrpsyms, &
                    &orientations(:batchsz), fpls(:batchsz), zbatch(:,:batchsz), dens_dummy(:,:,:batchsz), &
                    &valid_e(:batchsz), batchsz)
                call insert_planes_oversamp_coupled_batch_scaled(Yodd, rho_o, build%pgrpsyms, &
                    &orientations(:batchsz), fpls(:batchsz), zbatch(:,:batchsz), dens_dummy(:,:,:batchsz), &
                    &valid_o(:batchsz), batchsz)
                if( batchlims(2)==nptcls .or. mod(batchlims(2), 5*MAXIMGBATCHSZ)==0 )then
                    write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA PROBE PASS PARTICLES: ',batchlims(2),' / ',nptcls
                    call flush(logfhandle)
                endif
            end do
            ! finalize even/odd Y_q, half-set FSC Wiener-merge -> clean band-limited masked basis volume.
            allocate(realvols(ncomp))
            filtsz = max(1, fdim(params%box_crop) - 1)
            allocate(filt(filtsz), corrs(filtsz))
            do q = 1, ncomp
                ! even half -> band-limited real image (UNmasked, for an unbiased FSC)
                Yeven(q)%rho_exp = rho_e(1,:,:,:)
                call Yeven(q)%compress_exp; call Yeven(q)%sampl_dens_correct; call Yeven(q)%ifft
                call Yeven(q)%div(real(params%box))
                call realvols(q)%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
                call Yeven(q)%get_rmat_ptr(rmatp); call realvols(q)%set_rmat(rmatp, .false.)
                ! odd half
                Yodd(q)%rho_exp = rho_o(1,:,:,:)
                call Yodd(q)%compress_exp; call Yodd(q)%sampl_dens_correct; call Yodd(q)%ifft
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
                mu_q = sum(z(:,min(q,ncomp))) / real(nptcls,dp)
                sd_q = sum((z(:,min(q,ncomp))-mu_q)**2) / real(max(1,nptcls-1),dp)
                eigvals(q) = max(sd_q, DTINY)
                ! overwrite the eigenvolume MRC with the refined basis vector
                fname = 'flex_pca_pc'//int2str_pad(q,3)//MRC_EXT
                call utilde_real(q)%write(fname, del_if_exists=.true.); call fname%kill
            end do
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
            deallocate(Gth, Ath, bth, cth, zth, basis_fpls, mean_fpl, orientations, zbatch, dens_dummy, valid, valid_e, valid_o, eo)
        end do
        deallocate(z)
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

    !> In-place symmetric positive-definite solve A x = b (b overwritten by x) via Cholesky the
    !! factorization then fails on rounding, the fixed 1e-8 ridge swamps the matrix by ~1e8, and the solve
    !! returns the b=0 fallback for essentially every particle -- which collapsed 99% of particles onto one
    !! latent and put 99.6% of them into a single state.
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

    !> Complex Hermitian-half plane inner product over the native k<=0 half (the planes are stored
    !! half-plane, k in [kmin,0]). The properness fix (2.1) is applied at the Sbb accumulation (use
    !! Re(b_q)Re(b_r), not |b|^2) together with the 0.5*sig2 noise scaling in reduced_covariance_solve,
    !! since E[Re(b_q)Re(b_r)]_noise = 0.5*sig2*G.
    function cov_herm_inner( lhs, rhs, half ) result( val )
        type(fplane_type), intent(in) :: lhs, rhs
        integer, optional, intent(in) :: half
        complex(dp) :: val, acc
        integer :: h, k, hmin, hmax, kmin, kmax, nyq_eff, pf, h_hi, hlf, par
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
        do k = kmin, kmax, pf
            ! 2.6:
            h_hi = hmax
            if( k == 0 ) h_hi = 0
            do h = hmin, h_hi, pf
                if( nint(sqrt(real(h*h + k*k))) > nyq_eff ) cycle
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

    !> DIAGNOSTIC: The reduced-solve debias assumes E[Re(b_q)Re(b_r)]_noise = 0.5*sig2*G with b,G from
    !! cov_herm_inner, so the sig2 fed to it must be the noise variance in THIS index convention.
    subroutine cov_herm_selfpower( fpl, pw, cnt )
        type(fplane_type), intent(in)  :: fpl
        real(dp),          intent(out) :: pw, cnt
        integer     :: h, k, hmin, hmax, kmin, kmax, nyq_eff, pf, h_hi
        complex(dp) :: c
        pf  = OSMPL_PAD_FAC
        pw  = 0.d0; cnt = 0.d0
        nyq_eff = fpl%nyq
        if( nyq_eff <= 0 ) nyq_eff = ubound(fpl%cmplx_plane,1)
        hmin = max(pf*ceil_div(lbound(fpl%cmplx_plane,1),pf), pf*ceil_div(-nyq_eff,pf))
        hmax = min(pf*floor_div(ubound(fpl%cmplx_plane,1),pf), pf*floor_div(nyq_eff,pf))
        kmin = max(pf*ceil_div(lbound(fpl%cmplx_plane,2),pf), pf*ceil_div(-nyq_eff,pf))
        kmax = min(0, pf*floor_div(nyq_eff,pf))
        do k = kmin, kmax, pf
            h_hi = hmax
            if( k == 0 ) h_hi = 0
            do h = hmin, h_hi, pf
                if( nint(sqrt(real(h*h + k*k))) > nyq_eff ) cycle
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

    !> svec packing is an ISOMETRY:
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
    subroutine test_flex_pca_packed_solve()
        integer,  parameter :: D = 3, NPK = D*(D+1)/2, NG = 3
        real(dp) :: G(D,D,NG), Sig(D,D), Dmat(D,D), Ms(NPK,NPK), rhs(NPK), x(NPK), Sig_est(D,D)
        real(dp) :: sv(NPK,NG), tmp(D,D), err, relres
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
        call cg_solve_packed(Ms, NPK, D, rhs, 0.d0, COV_CG_MAXIT, COV_CG_TOL, x, iters, relres)
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
