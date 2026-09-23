!@descr: unit test routines for SVD, probabilistic and kernel PCA (simple_pca_svd, simple_ppca, simple_kpca_svd)
! The analysers behind ppca_denoise, cluster2D and the class split: SVD PCA on both the D >= N and the
! transposed D < N branches, probabilistic PCA (EM against the Tipping & Bishop maximum-likelihood solution,
! BIC rank suggestion, external reconstruction) and kernel PCA (exact and Nystroem backends, cosine and RBF
! kernels: kernel eigenvalues, feature vectors and pre-images pinned on a double-precision emulation of the
! same pipeline).
module simple_pca_tester
use simple_test_utils      ! assertions etc.
use simple_defs            ! sp, dp, logfhandle
use simple_string,         only: string
use simple_string_utils,   only: int2str
use simple_syslib,         only: del_file
use simple_pca_svd,        only: pca_svd
use simple_ppca,           only: ppca
use simple_kpca_svd,       only: kpca_svd, suggest_kpca_nystrom_neigs
implicit none
private
public :: run_all_pca_tests

!---------------- SVD PCA fixtures (references: pca_svd_ref.py, stats batch 2026-09-23) ----------------

! A: D = 5 >= N = 4 (master_ori); singular values of the centred data 10.71309, 5.69645, 2.45564, 0
real, parameter :: PCA_A(5,4) = reshape([ 1., 3., -1., 2., 0.,   2., 1., 0., -2., 5., &
                                         &3., 5.,  4., 1., -3.,  4., 8., 10., 3., 1.], [5,4])
real, parameter :: A_ENERGIES(3) = [114.77033, 32.4495163, 6.0301537]   ! w^2 per component
real, parameter :: A_FEAT1(4)    = [3.749356, 6.004916, -1.954746, -7.799526]
real, parameter :: A_RANK1(5,4)  = reshape([ 1.8040614, 2.4559543,  0.3241060, 0.0243782, 1.6676305, &
                                            &1.3853944, 1.3766814, -1.4360711,-0.5625423, 2.2196643, &
                                            &2.8628312, 5.1853349,  4.7754299, 1.5086454, 0.2715886, &
                                            &3.9477131, 7.9820293,  9.3365352, 3.0295186,-1.1588834], [5,4])
real, parameter :: A_RANK2(5,4)  = reshape([ 1.485989, 2.628556, -0.660858, 0.726860,-0.629858, &
                                            &1.805257, 1.148843, -0.135900,-1.489832, 5.252394, &
                                            &2.499673, 5.382402,  3.650852, 2.310701,-2.351559, &
                                            &4.209082, 7.840198, 10.145906, 2.452271, 0.729023], [5,4])
! B: D = 3 < N = 6 (master_T); singular values 6.95192, 3.25915, 2.17149
real, parameter :: PCA_B(3,6) = reshape([ 2., 1., -3.,  1., 3., 2.,  -1., 2., 1., &
                                         &4., -2., 0.,  0., 5., 2.,   3., 1., -1.], [3,6])
real, parameter :: B_ENERGIES(3) = [48.3292353, 10.6220518, 4.7153795]
real, parameter :: B_FEAT1(6)    = [-2.215050, 2.063195, 1.956683, -3.997955, 4.007651, -1.814524]
real, parameter :: B_RANK1(3,6)  = reshape([ 2.6805986, 0.1034326,-0.8672252,  0.4003385, 3.1227316, 1.1296790, &
                                            &0.4571082, 3.0475628, 1.0799639,  3.6308685,-1.1548226,-1.6994101, &
                                            &-0.6360362,4.4949985, 2.0372690,  2.4671224, 0.3860971,-0.6802766], [3,6])
real, parameter :: B_RANK2(3,6)  = reshape([ 2.057324, 1.042113,-2.998216,  0.651721, 2.744139, 1.989161, &
                                            &0.446891, 3.062951, 1.045031,  4.129089,-1.905165, 0.004018, &
                                            &-0.652880,4.520365, 1.979681,  2.367855, 0.535598,-1.019674], [3,6])

!---------------- probabilistic PCA fixture (references: ppca_ref.py) ----------------

! D = 5, N = 16: a rank-2 signal plus noise; the sample covariance of the centred data has eigenvalues
! 7.175061, 2.488612, 0.103408, 0.065567, 0.033428. The maximum-likelihood PPCA with Q = 2 (Tipping &
! Bishop 1999) has sigma^2 = mean of the three discarded eigenvalues, retained eigenvalues equal to the
! top two, and reconstructs x as sum_k (1 - sigma^2/lambda_k) u_k u_k^T x
integer, parameter :: PPCA_D = 5, PPCA_N = 16, PPCA_Q = 2
real, parameter :: PPCA_X(5,16) = reshape([ &
        &  0.481,  -2.178,  -1.396,   0.843,  -1.019, &
        &  0.870,   0.471,  -0.737,   0.406,   0.539, &
        &  0.730,  -2.900,  -1.131,   1.366,  -1.957, &
        & -2.101,  -3.598,   0.038,   0.388,  -3.082, &
        & -0.218,  -3.318,  -1.698,   1.024,  -1.818, &
        & -2.448,  -1.591,   0.384,  -0.797,  -2.445, &
        &  0.745,  -1.822,  -0.951,   0.967,  -0.356, &
        &  3.716,   2.594,  -1.633,   0.429,   3.645, &
        & -1.497,  -0.515,   0.409,   0.005,  -1.233, &
        & -1.679,  -0.946,   0.345,   0.004,  -1.568, &
        &  2.540,  -2.442,  -2.222,   2.136,   0.009, &
        &  1.399,  -0.498,  -0.550,   0.840,   0.675, &
        &  0.814,   0.155,  -0.250,  -0.021,   0.134, &
        & -3.332,  -1.358,   0.788,  -0.260,  -1.929, &
        &  1.149,  -2.224,  -1.419,   1.110,  -0.462, &
        &  2.367,  -0.094,  -0.501,   0.964,   1.403 &
        &], [5,16])
real, parameter :: PPCA_REC(5,16) = reshape([ &
        & 0.305628, -0.936582, -0.431144,  0.421683, -0.388691, &
        & 0.704368,  1.521119,  0.135618, -0.253968,  1.269092, &
        & 0.212138, -1.795984, -0.650930,  0.674428, -0.943978, &
        &-2.395187, -2.161589,  0.462878, -0.139431, -2.544042, &
        &-0.223748, -2.035748, -0.519551,  0.597881, -1.316577, &
        &-2.795938, -0.335200,  1.211828, -0.891673, -1.682362, &
        & 0.566764, -0.412258, -0.393022,  0.340723,  0.058504, &
        & 3.630366,  3.765692, -0.551389,  0.048302,  4.143914, &
        &-1.686446,  0.603411,  0.978181, -0.806207, -0.540793, &
        &-1.877018,  0.199706,  0.943896, -0.740108, -0.879539, &
        & 2.400217, -1.201332, -1.497308,  1.261536,  0.568149, &
        & 1.122366,  0.896144, -0.252734,  0.104232,  1.123423, &
        & 0.313115,  1.364620,  0.271563, -0.342237,  0.969191, &
        &-3.169644, -0.022956,  1.483378, -1.129798, -1.697160, &
        & 0.950087, -0.834375, -0.702813,  0.618902,  0.013766, &
        & 1.942934,  1.385331, -0.488450,  0.235734,  1.847104 &
        &], [5,16])
real, parameter :: PPCA_EIGVALS(2) = [7.175061, 2.488612]
real, parameter :: PPCA_SIGMA2     = 0.0674676
real, parameter :: PPCA_BIC_Q2     = 171.46     ! -2 ln L + p ln N from the PPCA marginal likelihood at the ML solution

!---------------- kernel PCA fixture (references: kpca_ref.py / kpca_nys.py) ----------------

! D = 4, N = 12: two tight clusters of six points; Q = 2; the automatic RBF gamma (inverse mean squared
! distance over all pairs) is 0.1086098
integer, parameter :: KPCA_D = 4, KPCA_N = 12, KPCA_Q = 2
real, parameter :: KPCA_X(4,12) = reshape([ &
        &  2.10,   0.30,  -0.95,   0.15, &
        &  1.85,   0.60,  -0.80,  -0.05, &
        &  2.05,   0.75,  -1.10,  -0.20, &
        &  2.20,   0.55,  -0.85,   0.10, &
        &  1.90,   0.35,  -1.05,   0.25, &
        &  2.15,   0.40,  -1.20,  -0.10, &
        & -1.10,   1.70,   0.45,   1.85, &
        & -0.85,   1.40,   0.30,   2.05, &
        & -1.05,   1.25,   0.60,   2.20, &
        & -1.20,   1.45,   0.35,   1.90, &
        & -0.90,   1.65,   0.55,   1.75, &
        & -1.15,   1.60,   0.70,   2.10 &
        &], [4,12])
real, parameter :: CLUSTER1(4) = [2.0, 0.5, -1.0, 0.0]
real, parameter :: CLUSTER2(4) = [-1.0, 1.5, 0.5, 2.0]
! exact backend, cosine kernel: eigenvalues of the centred kernel, features sqrt(lambda_k) v_k, and the
! converged pre-images of the fixed-point rule with non-negative weights (three iterations each)
real, parameter :: COS_EIGVALS(2) = [7.62999, 0.078635]
real, parameter :: COS_FEAT1(12) = [ &
        &-0.806915, -0.785007, -0.797745, -0.779530, -0.770534, -0.842254,  0.795400,  0.750786,  0.815735,  0.823342,  0.775491,  0.821231 &
        &]
real, parameter :: COS_FEAT2(12) = [ &
        &-0.097485,  0.081444,  0.134318, -0.005442, -0.106983, -0.007005,  0.081774, -0.063999, -0.117686,  0.004946,  0.090818,  0.005300 &
        &]
real, parameter :: COS_PREIMG(4,12) = reshape([ &
        & 2.052100,  0.370292, -1.010764,  0.137840, &
        & 2.023993,  0.625338, -0.988393, -0.098190, &
        & 2.011934,  0.641415, -0.987370, -0.108628, &
        & 2.044542,  0.479554, -0.994690,  0.033029, &
        & 2.049198,  0.368041, -1.010397,  0.141088, &
        & 2.044724,  0.477437, -0.994803,  0.034992, &
        &-1.054667,  1.617982,  0.504905,  1.872520, &
        &-1.037696,  1.393185,  0.490939,  2.083412, &
        &-1.024739,  1.367367,  0.491804,  2.102494, &
        &-1.046295,  1.516103,  0.494222,  1.968812, &
        &-1.055665,  1.623874,  0.508269,  1.866533, &
        &-1.046336,  1.516703,  0.494240,  1.968253 &
        &], [4,12])
! exact backend, RBF kernel with the automatic gamma: converged pre-images (4 to 7 iterations)
real, parameter :: RBF_EIGVALS(2) = [4.901204, 0.061867]
real, parameter :: RBF_FEAT1(12) = [ &
        & 0.643752,  0.613861,  0.646519,  0.644045,  0.620036,  0.664932, -0.643752, -0.613861, -0.646519, -0.644045, -0.620036, -0.664932 &
        &]
real, parameter :: RBF_FEAT2(12) = [ &
        &-0.080966,  0.055207,  0.117325, -0.006518, -0.086732,  0.000532,  0.080966, -0.055207, -0.117325,  0.006518,  0.086732, -0.000532 &
        &]
real, parameter :: RBF_PREIMG(4,12) = reshape([ &
        & 1.991156,  0.393848, -0.971766,  0.176139, &
        & 2.038978,  0.610322, -0.996873, -0.090008, &
        & 1.885713,  0.685395, -0.931291, -0.023948, &
        & 2.045466,  0.477712, -0.993568,  0.036366, &
        & 1.977836,  0.392973, -0.966942,  0.190585, &
        & 2.044784,  0.492173, -0.993903,  0.022567, &
        &-0.991156,  1.606152,  0.471766,  1.823861, &
        &-1.038978,  1.389678,  0.496873,  2.090008, &
        &-0.885713,  1.314605,  0.431291,  2.023948, &
        &-1.045466,  1.522288,  0.493568,  1.963634, &
        &-0.977836,  1.607027,  0.466942,  1.809415, &
        &-1.044784,  1.507827,  0.493903,  1.977433 &
        &], [4,12])
! Nystroem backend with every point a landmark: the same eigenvalues, features and (RBF) pre-images as the
! exact backend; the cosine pre-image is the power-weighted cosine average of the landmarks, independent of
! the components
real, parameter :: NYS_COS_PREIMG(4,12) = reshape([ &
        & 2.041797,  0.491465, -0.991465,  0.025103, &
        & 2.041793,  0.491505, -0.991462,  0.025066, &
        & 2.041792,  0.491516, -0.991467,  0.025052, &
        & 2.041795,  0.491485, -0.991460,  0.025085, &
        & 2.041796,  0.491462, -0.991468,  0.025104, &
        & 2.041795,  0.491485, -0.991470,  0.025080, &
        &-1.041879,  1.508660,  0.491704,  1.974868, &
        &-1.041877,  1.508641,  0.491703,  1.974886, &
        &-1.041879,  1.508632,  0.491706,  1.974896, &
        &-1.041880,  1.508649,  0.491703,  1.974879, &
        &-1.041878,  1.508662,  0.491705,  1.974867, &
        &-1.041879,  1.508649,  0.491707,  1.974880 &
        &], [4,12])
character(len=*), parameter :: SCRATCH_LOG = 'pca_tester_scratch.log'
integer :: logfhandle_saved = 0

contains

    subroutine run_all_pca_tests()
        write(*,'(A)') '**** running all PCA tests ****'
        call test_pca_svd_tall()
        call test_pca_svd_wide()
        call test_ppca_ml_solution()
        call test_ppca_rank_suggestion()
        call test_kpca_exact_cosine()
        call test_kpca_exact_rbf()
        call test_kpca_nystrom_all_landmarks()
        call test_kpca_nystrom_landmark_subset()
        call test_kpca_suggest_neigs()
    end subroutine run_all_pca_tests

    !---------------- helpers ----------------

    subroutine centre( x, avg, xc )
        real, intent(in)  :: x(:,:)
        real, intent(out) :: avg(size(x,1)), xc(size(x,1),size(x,2))
        integer :: j
        avg = sum(x, dim=2) / real(size(x,2))
        do j = 1,size(x,2)
            xc(:,j) = x(:,j) - avg
        end do
    end subroutine centre

    ! a feature row is defined up to its sign: align on the first sizeable entry, then compare
    subroutine assert_row_upto_sign( expected, actual, tol, label )
        real,             intent(in) :: expected(:), actual(:), tol
        character(len=*), intent(in) :: label
        real    :: s
        integer :: i, ipiv
        logical :: l_ok
        ipiv = maxloc(abs(expected), dim=1)
        s    = sign(1., expected(ipiv)) * sign(1., actual(ipiv))
        l_ok = .true.
        do i = 1,size(expected)
            if( abs(s*actual(i) - expected(i)) > tol ) l_ok = .false.
        end do
        call assert_true(l_ok, label)
        if( .not. l_ok ) write(*,'(A,12F10.5)') '   got (sign-aligned): ', s*actual
    end subroutine assert_row_upto_sign

    subroutine assert_columns_equal( expected, actual, tol, label )
        real,             intent(in) :: expected(:,:), actual(:,:), tol
        character(len=*), intent(in) :: label
        integer :: j
        real    :: maxdev
        maxdev = 0.
        do j = 1,size(expected,2)
            maxdev = max(maxdev, maxval(abs(expected(:,j) - actual(:,j))))
        end do
        call assert_true(maxdev <= tol, label)
        if( maxdev > tol ) write(*,'(A,ES10.3)') '   largest deviation: ', maxdev
    end subroutine assert_columns_equal

    ! how many pre-images lie closer to their own cluster centre than the input point did
    integer function count_closer( gen )
        real, intent(in) :: gen(KPCA_D,KPCA_N)
        integer :: j
        count_closer = 0
        do j = 1,KPCA_N
            if( j <= 6 )then
                if( sqrt(sum((gen(:,j)-CLUSTER1)**2)) < sqrt(sum((KPCA_X(:,j)-CLUSTER1)**2)) ) count_closer = count_closer + 1
            else
                if( sqrt(sum((gen(:,j)-CLUSTER2)**2)) < sqrt(sum((KPCA_X(:,j)-CLUSTER2)**2)) ) count_closer = count_closer + 1
            endif
        end do
    end function count_closer

    ! the Nystroem backend reports its pre-image progress on logfhandle; keep that out of the gate log
    subroutine divert_log( unit )
        integer, intent(out) :: unit
        logfhandle_saved = logfhandle
        open(newunit=unit, file=SCRATCH_LOG, status='replace', action='write')
        logfhandle = unit
    end subroutine divert_log

    subroutine restore_log( unit )
        integer, intent(in) :: unit
        close(unit)
        logfhandle = logfhandle_saved
        call del_file(string(SCRATCH_LOG))
    end subroutine restore_log

    !---------------- SVD PCA ----------------

    ! D >= N: the left singular vectors come straight from the SVD of the centred data; with all N
    ! components the data are reproduced exactly (the centred N columns span at most N-1 dimensions),
    ! with Q = 1 or 2 the best rank-Q approximation; the feature energies are the squared singular values
    subroutine test_pca_svd_tall()
        type(pca_svd) :: pca
        real    :: avg(5), xc(5,4), gen(5,4), feats(4,4), energies(4), tmp(5)
        integer :: j, k
        write(*,'(A)') 'test_pca_svd_tall'
        call centre(PCA_A, avg, xc)
        call pca%new(4, 5, 4)
        call pca%master(xc)
        do j = 1,4
            call pca%generate(j, avg, tmp)
            gen(:,j) = tmp
        end do
        call assert_columns_equal(PCA_A, gen, 1.e-4, 'Q = N: generate reproduces the data (exact reconstruction)')
        do j = 1,4
            feats(:,j) = pca%get_feat(j)
        end do
        do k = 1,4
            energies(k) = sum(feats(k,:)**2)
        end do
        call assert_real(A_ENERGIES(1), energies(1), 1.e-2, 'feature energy 1 is the largest squared singular value')
        call assert_real(A_ENERGIES(2), energies(2), 1.e-2, 'feature energy 2 is the second squared singular value')
        call assert_real(A_ENERGIES(3), energies(3), 1.e-2, 'feature energy 3 is the third squared singular value')
        call assert_real(0., energies(4), 1.e-3, 'the fourth component of four centred columns carries no energy')
        call assert_row_upto_sign(A_FEAT1, feats(1,:), 1.e-3, 'first feature vector is u_1^T x (up to sign)')
        call pca%new(4, 5, 1)
        call pca%master(xc)
        do j = 1,4
            call pca%generate(j, avg, tmp)
            gen(:,j) = tmp
        end do
        call assert_columns_equal(A_RANK1, gen, 1.e-3, 'Q = 1: generate is the best rank-1 approximation plus the mean')
        call pca%new(4, 5, 2)
        call pca%master(xc)
        do j = 1,4
            call pca%generate(j, avg, tmp)
            gen(:,j) = tmp
        end do
        call assert_columns_equal(A_RANK2, gen, 1.e-3, 'Q = 2: generate is the best rank-2 approximation plus the mean')
        call pca%kill
    end subroutine test_pca_svd_tall

    ! D < N takes the transposed route (SVD of X^T); the results must be the same PCA
    subroutine test_pca_svd_wide()
        type(pca_svd) :: pca
        real    :: avg(3), xc(3,6), gen(3,6), feats(3,6), energies(3), tmp(3)
        integer :: j, k
        write(*,'(A)') 'test_pca_svd_wide'
        call centre(PCA_B, avg, xc)
        call pca%new(6, 3, 3)
        call pca%master(xc)
        do j = 1,6
            call pca%generate(j, avg, tmp)
            gen(:,j) = tmp
        end do
        call assert_columns_equal(PCA_B, gen, 1.e-4, 'Q = D: generate reproduces the data (exact reconstruction)')
        do j = 1,6
            feats(:,j) = pca%get_feat(j)
        end do
        do k = 1,3
            energies(k) = sum(feats(k,:)**2)
        end do
        call assert_real(B_ENERGIES(1), energies(1), 1.e-2, 'transposed route: feature energy 1')
        call assert_real(B_ENERGIES(2), energies(2), 1.e-2, 'transposed route: feature energy 2')
        call assert_real(B_ENERGIES(3), energies(3), 1.e-2, 'transposed route: feature energy 3')
        call assert_row_upto_sign(B_FEAT1, feats(1,:), 1.e-3, 'transposed route: first feature vector (up to sign)')
        call pca%new(6, 3, 1)
        call pca%master(xc)
        do j = 1,6
            call pca%generate(j, avg, tmp)
            gen(:,j) = tmp
        end do
        call assert_columns_equal(B_RANK1, gen, 1.e-3, 'transposed route, Q = 1: best rank-1 approximation')
        call pca%new(6, 3, 2)
        call pca%master(xc)
        do j = 1,6
            call pca%generate(j, avg, tmp)
            gen(:,j) = tmp
        end do
        call assert_columns_equal(B_RANK2, gen, 1.e-3, 'transposed route, Q = 2: best rank-2 approximation')
        call pca%kill
    end subroutine test_pca_svd_wide

    !---------------- probabilistic PCA ----------------

    ! EM from a random start converges to the unique maximum of the PPCA likelihood; the tolerances
    ! reflect the linear convergence along the slow (largest-eigenvalue) direction at the built-in
    ! stopping thresholds
    subroutine test_ppca_ml_solution()
        type(ppca) :: prob
        real     :: avg(PPCA_D), xc(PPCA_D,PPCA_N), gen(PPCA_D,PPCA_N), tmp(PPCA_D), ext(PPCA_D)
        real     :: feats(PPCA_Q,PPCA_N)
        real, allocatable :: eigvals(:), signal(:)
        real(dp) :: sigma2
        integer  :: j
        write(*,'(A)') 'test_ppca_ml_solution'
        call centre(PPCA_X, avg, xc)
        call prob%new(PPCA_N, PPCA_D, PPCA_Q)
        call prob%set_verbose(.false.)
        call prob%master(xc, 2000)
        eigvals = prob%get_eigvals()
        signal  = prob%get_signal_eigvals()
        sigma2  = prob%get_sigma2()
        call assert_int(PPCA_Q, size(eigvals), 'one retained eigenvalue per component')
        call assert_real(PPCA_EIGVALS(1), eigvals(1), 0.05,  'largest retained eigenvalue is the top covariance eigenvalue')
        call assert_real(PPCA_EIGVALS(2), eigvals(2), 0.01,  'second retained eigenvalue is the second covariance eigenvalue')
        call assert_real(PPCA_SIGMA2, real(sigma2), 1.e-3,   'sigma^2 is the mean of the discarded covariance eigenvalues')
        call assert_real(real(sigma2), eigvals(1) - signal(1), 1.e-5, 'eigenvalue = signal eigenvalue + sigma^2 (component 1)')
        call assert_real(real(sigma2), eigvals(2) - signal(2), 1.e-5, 'eigenvalue = signal eigenvalue + sigma^2 (component 2)')
        call assert_true(eigvals(1) > eigvals(2), 'retained eigenvalues are sorted in descending order')
        do j = 1,PPCA_N
            call prob%generate(j, avg, tmp)
            gen(:,j) = tmp - avg
            feats(:,j) = prob%get_feat(j)
        end do
        call assert_columns_equal(PPCA_REC, gen, 2.e-3, 'generate reconstructs sum_k (1 - sigma^2/lambda_k) u_k u_k^T x')
        call assert_real(0.526628, gen(1,1) + avg(1), 2.e-3, 'generate adds the mean back (column 1)')
        call assert_real(2.163934, gen(1,16) + avg(1), 2.e-3, 'generate adds the mean back (column 16)')
        ! posterior latent means: |z_k| = sqrt(lambda_k - sigma^2)/lambda_k |u_k^T x|
        call assert_real(1.528147, abs(feats(1,4)),  0.02, 'latent mean of component 1, column 4')
        call assert_real(2.438430, abs(feats(1,8)),  0.02, 'latent mean of component 1, column 8')
        call assert_real(1.814753, abs(feats(2,11)), 0.02, 'latent mean of component 2, column 11')
        ! reconstruct_external on a training column reproduces generate (the model is the same)
        call prob%reconstruct_external(xc(:,3), ext)
        call assert_columns_equal(reshape(gen(:,3),[PPCA_D,1]), reshape(ext,[PPCA_D,1]), 2.e-3, &
            &'reconstruct_external of a training column equals its generate output')
        call assert_real(PPCA_BIC_Q2, real(prob%calc_bic(xc)), 0.5, 'BIC of the converged Q = 2 model (marginal likelihood, 10 parameters)')
        call prob%kill
    end subroutine test_ppca_ml_solution

    ! suggest_rank fits each candidate (EM to its own tolerances within the caller's cap) and takes the
    ! smallest rank whose BIC is within 2 of the best. With the marginal-likelihood BIC the converged values
    ! on this fixture are 249.8 (rank 1), 171.5 (rank 2), 176.7 (rank 3): the scan stops where the spectrum
    ! flattens (ten iterations, the former hard cap, left the margin within the tolerance on Linux)
    subroutine test_ppca_rank_suggestion()
        type(ppca) :: prob
        real     :: avg(PPCA_D), xc(PPCA_D,PPCA_N)
        integer  :: best_q
        integer,  allocatable :: qs(:)
        real(dp), allocatable :: bics(:), sigma2s(:)
        write(*,'(A)') 'test_ppca_rank_suggestion'
        call centre(PPCA_X, avg, xc)
        best_q = prob%suggest_rank(xc, [1, 2, 3], 500, qs, bics, sigma2s)
        call assert_int(2, best_q, 'the rank-2-plus-noise fixture is recognised as rank 2')
        call assert_int(3, size(qs), 'one entry per candidate')
        call assert_int(1, qs(1), 'candidate ranks are reported')
        call assert_int(2, qs(2), 'candidate ranks are reported (middle)')
        call assert_int(3, qs(3), 'candidate ranks are reported (last)')
        call assert_true(bics(2) < bics(3) .and. bics(3) < bics(1), 'BIC ranks the candidates 2 < 3 < 1')
        call assert_true(bics(3) - bics(2) > 2._dp, 'rank 3 loses to rank 2 by more than the tolerance')
        call assert_real(171.46, real(bics(2)), 0.5, 'rank 2 BIC at convergence')
        call assert_real(176.70, real(bics(3)), 0.5, 'rank 3 BIC at convergence')
        call assert_true(bics(1) - bics(2) > 50._dp, 'rank 1 loses by far (the second signal component unexplained)')
        call assert_true(sigma2s(1) > sigma2s(2) .and. sigma2s(2) > sigma2s(3), 'sigma^2 shrinks with the rank (converged fits)')
        call assert_real(0.6728, real(sigma2s(1)), 5.e-3, 'rank 1: sigma^2 is the mean of the four discarded eigenvalues')
        call assert_real(0.0675, real(sigma2s(2)), 2.e-3, 'rank 2: sigma^2 is the mean of the three discarded eigenvalues')
        best_q = prob%suggest_rank(xc, [2, 2, 2, 4], 500, qs, bics, sigma2s)
        call assert_int(2, best_q, 'repeated candidates are fitted once; rank 4 does not beat rank 2')
        call assert_int(0, qs(2), 'a repeated candidate is reported as rank 0')
        call assert_int(0, qs(3), 'a twice-repeated candidate is reported as rank 0 too (compared with the previous fit, not the previous slot)')
        call assert_int(4, qs(4), 'the next distinct candidate is fitted')
        call prob%kill
    end subroutine test_ppca_rank_suggestion

    !---------------- kernel PCA ----------------

    ! exact backend, cosine kernel: the centred cosine-similarity kernel's leading eigenpairs, the
    ! features sqrt(lambda_k) v_k(i), and the converged pre-images of the fixed-point rule with
    ! non-negative weights (max(0, projected column) x max(0, cosine); with the sign-mixed weights and
    ! the L1 denominator of before 2026-09-23 one cluster never converged and landed in the other)
    subroutine test_kpca_exact_cosine()
        type(kpca_svd) :: kpca
        real, allocatable :: eigvals(:)
        real    :: feats(KPCA_Q,KPCA_N), gen(KPCA_D,KPCA_N), tmp(KPCA_D), zero(KPCA_D)
        integer :: j
        write(*,'(A)') 'test_kpca_exact_cosine'
        zero = 0.
        call kpca%new(KPCA_N, KPCA_D, KPCA_Q)
        call kpca%set_params(nthr=1, kpca_ker='cosine', kpca_backend='exact')
        call kpca%master(KPCA_X)
        eigvals = kpca%get_eigvals()
        call assert_int(KPCA_Q, size(eigvals), 'one kernel eigenvalue per component')
        call assert_real(COS_EIGVALS(1), eigvals(1), 2.e-4, 'cosine kernel: leading eigenvalue')
        call assert_real(COS_EIGVALS(2), eigvals(2), 1.e-4, 'cosine kernel: second eigenvalue')
        do j = 1,KPCA_N
            feats(:,j) = kpca%get_feat(j)
            call kpca%generate(j, zero, tmp)
            gen(:,j) = tmp
        end do
        call assert_row_upto_sign(COS_FEAT1, feats(1,:), 2.e-3, 'cosine kernel: first features are sqrt(lambda_1) v_1')
        call assert_row_upto_sign(COS_FEAT2, feats(2,:), 2.e-3, 'cosine kernel: second features are sqrt(lambda_2) v_2')
        call assert_real(COS_EIGVALS(1), sum(feats(1,:)**2), 2.e-3, 'cosine kernel: feature energy equals the eigenvalue')
        call assert_columns_equal(COS_PREIMG, gen, 3.e-3, 'cosine kernel: converged pre-images')
        call assert_int(KPCA_N, count_closer(gen), 'cosine kernel: every pre-image is closer to its cluster centre than the input')
        call kpca%kill
    end subroutine test_kpca_exact_cosine

    ! exact backend, RBF kernel with the automatic gamma: the pre-image iteration converges in a few steps
    ! and pulls every point towards its cluster
    subroutine test_kpca_exact_rbf()
        type(kpca_svd) :: kpca
        real, allocatable :: eigvals(:)
        real    :: feats(KPCA_Q,KPCA_N), gen(KPCA_D,KPCA_N), tmp(KPCA_D), zero(KPCA_D), avg(KPCA_D)
        integer :: j
        write(*,'(A)') 'test_kpca_exact_rbf'
        zero = 0.
        call kpca%new(KPCA_N, KPCA_D, KPCA_Q)
        call kpca%set_params(nthr=1, kpca_ker='rbf', kpca_backend='exact', kpca_rbf_gamma=0.)
        call kpca%master(KPCA_X)
        eigvals = kpca%get_eigvals()
        call assert_real(RBF_EIGVALS(1), eigvals(1), 2.e-4, 'RBF kernel: leading eigenvalue (automatic gamma)')
        call assert_real(RBF_EIGVALS(2), eigvals(2), 1.e-4, 'RBF kernel: second eigenvalue')
        do j = 1,KPCA_N
            feats(:,j) = kpca%get_feat(j)
            call kpca%generate(j, zero, tmp)
            gen(:,j) = tmp
        end do
        call assert_row_upto_sign(RBF_FEAT1, feats(1,:), 2.e-3, 'RBF kernel: first features separate the clusters')
        call assert_row_upto_sign(RBF_FEAT2, feats(2,:), 2.e-3, 'RBF kernel: second features')
        call assert_columns_equal(RBF_PREIMG, gen, 3.e-3, 'RBF kernel: converged pre-images')
        call assert_int(KPCA_N, count_closer(gen), 'RBF kernel: every pre-image is closer to its cluster centre than the input')
        ! generate adds the supplied mean
        avg = 1.
        call kpca%generate(1, avg, tmp)
        call assert_columns_equal(reshape(gen(:,1)+1.,[KPCA_D,1]), reshape(tmp,[KPCA_D,1]), 1.e-6, 'generate adds the mean to the pre-image')
        call kpca%kill
    end subroutine test_kpca_exact_rbf

    ! Nystroem backend with every point a landmark (npts = N, no local support) reproduces the exact
    ! backend: eigenvalues, features (sqrt(lambda_k) v_k, up to sign) and the RBF pre-images (the same
    ! fixed point, with early stopping); the cosine pre-image follows the backend's own power-weighted
    ! average of the landmarks
    subroutine test_kpca_nystrom_all_landmarks()
        type(kpca_svd) :: kpca
        real, allocatable :: eigvals(:)
        real    :: feats(KPCA_Q,KPCA_N), gen(KPCA_D,KPCA_N), tmp(KPCA_D), zero(KPCA_D)
        integer :: j, unit
        write(*,'(A)') 'test_kpca_nystrom_all_landmarks'
        zero = 0.
        call divert_log(unit)
        ! cosine
        call kpca%new(KPCA_N, KPCA_D, KPCA_Q)
        call kpca%set_params(nthr=1, kpca_ker='cosine', kpca_backend='nystrom', kpca_nystrom_npts=KPCA_N)
        call kpca%master(KPCA_X)
        eigvals = kpca%get_eigvals()
        call assert_real(COS_EIGVALS(1), eigvals(1), 2.e-4, 'Nystroem/cosine: leading eigenvalue equals the exact one')
        call assert_real(COS_EIGVALS(2), eigvals(2), 1.e-4, 'Nystroem/cosine: second eigenvalue equals the exact one')
        do j = 1,KPCA_N
            feats(:,j) = kpca%get_feat(j)
            call kpca%generate(j, zero, tmp)
            gen(:,j) = tmp
        end do
        call assert_row_upto_sign(COS_FEAT1, feats(1,:), 2.e-3, 'Nystroem/cosine: first features equal the exact ones')
        call assert_row_upto_sign(COS_FEAT2, feats(2,:), 2.e-3, 'Nystroem/cosine: second features equal the exact ones')
        call assert_real(COS_EIGVALS(1), sum(feats(1,:)**2), 2.e-3, 'Nystroem/cosine: feature energy equals the eigenvalue')
        call assert_columns_equal(NYS_COS_PREIMG, gen, 2.e-3, 'Nystroem/cosine: power-weighted pre-images (cluster means)')
        ! rbf
        call kpca%new(KPCA_N, KPCA_D, KPCA_Q)
        call kpca%set_params(nthr=1, kpca_ker='rbf', kpca_backend='nystrom', kpca_nystrom_npts=KPCA_N, kpca_rbf_gamma=0.)
        call kpca%master(KPCA_X)
        eigvals = kpca%get_eigvals()
        call assert_real(RBF_EIGVALS(1), eigvals(1), 2.e-4, 'Nystroem/RBF: leading eigenvalue equals the exact one')
        call assert_real(RBF_EIGVALS(2), eigvals(2), 1.e-4, 'Nystroem/RBF: second eigenvalue equals the exact one')
        do j = 1,KPCA_N
            feats(:,j) = kpca%get_feat(j)
            call kpca%generate(j, zero, tmp)
            gen(:,j) = tmp
        end do
        call assert_row_upto_sign(RBF_FEAT1, feats(1,:), 2.e-3, 'Nystroem/RBF: first features equal the exact ones')
        call assert_row_upto_sign(RBF_FEAT2, feats(2,:), 2.e-3, 'Nystroem/RBF: second features equal the exact ones')
        call assert_columns_equal(RBF_PREIMG, gen, 3.e-3, 'Nystroem/RBF: pre-images equal the exact ones (early stopping allowed)')
        call kpca%kill
        call restore_log(unit)
    end subroutine test_kpca_nystrom_all_landmarks

    ! Nystroem backend with a landmark subset (6 of 12) and local support: the eigenvalues stay ordered and
    ! positive, the pre-images are convex combinations of data points (inside the data's bounding box)
    ! and land on the side of their own cluster
    subroutine test_kpca_nystrom_landmark_subset()
        type(kpca_svd) :: kpca
        real, allocatable :: eigvals(:)
        real    :: gen(KPCA_D,KPCA_N), tmp(KPCA_D), zero(KPCA_D), lo(KPCA_D), hi(KPCA_D), d_own, d_other
        integer :: j, unit, nown, ninside
        write(*,'(A)') 'test_kpca_nystrom_landmark_subset'
        zero = 0.
        lo   = minval(KPCA_X, dim=2)
        hi   = maxval(KPCA_X, dim=2)
        call divert_log(unit)
        call kpca%new(KPCA_N, KPCA_D, KPCA_Q)
        call kpca%set_params(nthr=1, kpca_ker='rbf', kpca_backend='nystrom', kpca_nystrom_npts=6, kpca_rbf_gamma=0.)
        call kpca%master(KPCA_X)
        call restore_log(unit)
        eigvals = kpca%get_eigvals()
        call assert_true(eigvals(1) > eigvals(2) .and. eigvals(2) > 0., 'landmark subset: eigenvalues ordered and positive')
        call assert_true(eigvals(1) > 1., 'landmark subset: the cluster split dominates the spectrum')
        nown    = 0
        ninside = 0
        do j = 1,KPCA_N
            call kpca%generate(j, zero, tmp)
            gen(:,j) = tmp
            if( all(tmp >= lo - 1.e-5) .and. all(tmp <= hi + 1.e-5) ) ninside = ninside + 1
            if( j <= 6 )then
                d_own   = sqrt(sum((tmp-CLUSTER1)**2))
                d_other = sqrt(sum((tmp-CLUSTER2)**2))
            else
                d_own   = sqrt(sum((tmp-CLUSTER2)**2))
                d_other = sqrt(sum((tmp-CLUSTER1)**2))
            endif
            if( d_own < d_other ) nown = nown + 1
        end do
        call assert_int(KPCA_N, ninside, 'landmark subset: pre-images stay inside the bounding box of the data')
        call assert_int(KPCA_N, nown,    'landmark subset: every pre-image lands on the side of its own cluster')
        call assert_true(all(abs(gen) < 10.), 'landmark subset: pre-images are finite and of data scale')
        call kpca%kill
    end subroutine test_kpca_nystrom_landmark_subset

    ! suggest_kpca_nystrom_neigs: the number of centred-feature gram eigenvalues carrying 99 % of the energy,
    ! but at least min(8, rank)
    subroutine test_kpca_suggest_neigs()
        write(*,'(A)') 'test_kpca_suggest_neigs'
        call assert_int(4, suggest_kpca_nystrom_neigs(KPCA_X, 'cosine'), 'cosine: rank 4 kernel (D = 4) caps the suggestion at 4')
        call assert_int(8, suggest_kpca_nystrom_neigs(KPCA_X, 'rbf'),    'rbf: 99 % energy needs 5, the floor of min(8, rank) wins')
        call assert_int(8, suggest_kpca_nystrom_neigs(KPCA_X, 'rbf', kpca_nystrom_npts=KPCA_N, kpca_rbf_gamma=0.), &
            &'rbf: explicit landmark count and automatic gamma give the same suggestion')
        call assert_int(1, suggest_kpca_nystrom_neigs(KPCA_X(:,1:1), 'rbf'), 'a single sample suggests one component')
    end subroutine test_kpca_suggest_neigs

end module simple_pca_tester
