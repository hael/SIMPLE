!@descr: unit tests of the PCG atomicity prior (simple_pcg_atomic_prior) and its hook in simple_reconstructor_pcg
! The projector onto the atom kernels on a synthetic box, the prior term, the width fit, and in the reconstructor:
! strength zero bit-identical, the term identical in both operator modes, amplitudes fitted through the data.
! Tolerances follow from single-precision maps (~1e-7) and the relative Gram ridge (1e-7), not from a run.
module simple_pcg_atomic_prior_tester
use simple_defs,              only: dp, PI
use simple_oris,              only: oris
use simple_image,             only: image
use simple_pcg_atomic_prior,  only: pcg_atomic_prior, pcg_atomic_prior_stats
use simple_reconstructor_pcg, only: reconstructor_pcg, PCG_OP_KERNEL, PCG_OP_MATRIXFREE, PCG_LAMBDA
use simple_test_utils,        only: assert_true
use ieee_arithmetic,          only: ieee_is_finite
implicit none
private
public :: run_all_pcg_atomic_prior_tests

integer,  parameter :: BOX          = 32
real,     parameter :: SMPD         = 0.4
real,     parameter :: ATOM_RAD     = 9.0    !< px: atoms lie within this radius of the centre
real,     parameter :: ATOM_DMIN    = 3.5    !< px: minimum atom spacing (kernels of ~1 px overlap weakly)
integer,  parameter :: NATOMS       = 60
real,     parameter :: SUPPORT_RAD  = 15.0   !< px: the solve support
real(dp), parameter :: PROJ_TOL     = 1.0e-5_dp  !< projector identities, relative
real(dp), parameter :: SYM_TOL      = 1.0e-6_dp  !< symmetry of the prior term, relative
real(dp), parameter :: GRAM_TOL     = 1.0e-6_dp  !< Gram solve against a known solution, relative
real,     parameter :: WIDTH_TOL    = 0.05       !< width fit on a noise-free map, relative
! reconstructor fixture
integer,  parameter :: RBOX         = 32
real,     parameter :: RSMPD        = 0.4
real,     parameter :: R_ATOM_RAD   = 7.0    !< atoms and their footprints stay inside the support
real,     parameter :: R_SUPPORT    = 12.0   !< four voxels from the box edge, the interior margin of test=pcg_recon
integer,  parameter :: R_NATOMS     = 25
real,     parameter :: R_SIGMA_VOX  = 1.2    !< phantom atom width = kernel width of the amplitude test
real(dp), parameter :: MODE_TOL     = 1.0e-5_dp  !< the prior term is the same in both operator modes, relative
real(dp), parameter :: INTERIOR_TOL = 5.0e-2_dp  !< kernel vs matrix-free with the prior on (EPS_INTERIOR of test=pcg_recon)
real(dp), parameter :: AMP_CORR_MIN = 0.99_dp    !< fitted against true amplitudes, dense noise-free views

contains

    subroutine run_all_pcg_atomic_prior_tests()
        write(*,'(A)') '**** running all pcg_atomic_prior tests ****'
        call test_projector()
        call test_apply_symmetric()
        call test_column_scaling()
        call test_fitted_model_projector()
        call test_gram_solve()
        call test_width_fit()
        call test_reconstructor_zero_strength()
        call test_reconstructor_operator_modes()
        call test_reconstructor_amplitudes()
    end subroutine run_all_pcg_atomic_prior_tests

    ! ---- fixtures --------------------------------------------------------------

    !> deterministic atom positions (3-D additive recurrence, no RNG) within rad of the box centre, dmin apart
    subroutine build_atoms( n, rad, dmin, nbox, pos )
        integer,           intent(in)  :: n, nbox
        real,              intent(in)  :: rad, dmin
        real, allocatable, intent(out) :: pos(:,:)
        real(dp), parameter :: G = 1.2207440846057595_dp   ! the root of x**4 = x + 1
        real(dp) :: a(3), u(3)
        real     :: p(3), acc(3,n)
        integer  :: it, na
        a  = [1._dp/G, 1._dp/G**2, 1._dp/G**3]
        na = 0
        it = 0
        do while( na < n .and. it < 100000 )
            it = it + 1
            u = modulo(0.5_dp + real(it,dp)*a, 1._dp)
            p = real((2._dp*u - 1._dp) * real(rad,dp))
            if( norm2(p) > rad ) cycle
            if( na > 0 )then
                if( any(norm2(acc(:,1:na) - spread(p,2,na), dim=1) < dmin) ) cycle
            endif
            na = na + 1
            acc(:,na) = p
        end do
        allocate(pos(3,na), source=acc(:,1:na) + real(nbox/2 + 1))
    end subroutine build_atoms

    !> 1 inside a sphere about the box centre
    function sphere_support( nbox, rad ) result( supp )
        integer, intent(in) :: nbox
        real,    intent(in) :: rad
        real, allocatable :: supp(:,:,:)
        real    :: c
        integer :: i, j, k
        c = real(nbox/2 + 1)
        allocate(supp(nbox,nbox,nbox), source=0.)
        do k = 1, nbox
            do j = 1, nbox
                do i = 1, nbox
                    if( (real(i)-c)**2 + (real(j)-c)**2 + (real(k)-c)**2 <= rad**2 ) supp(i,j,k) = 1.
                end do
            end do
        end do
    end function sphere_support

    !> a deterministic map with content everywhere, seed selects the pattern
    function probe_map( nbox, seed ) result( x )
        integer, intent(in) :: nbox, seed
        real, allocatable :: x(:,:,:)
        integer :: i, j, k
        allocate(x(nbox,nbox,nbox))
        do k = 1, nbox
            do j = 1, nbox
                do i = 1, nbox
                    x(i,j,k) = sin(0.37*real(seed)*real(i) + 0.11*real(j*k)) + cos(0.23*real(i*j) - 0.31*real(seed*k))
                end do
            end do
        end do
    end function probe_map

    !> B of a Gaussian of width sigma voxels
    real function bfac_of_sigma( sigma_vox, smpd_here )
        real, intent(in) :: sigma_vox, smpd_here
        bfac_of_sigma = 8. * PI * PI * (sigma_vox * smpd_here)**2
    end function bfac_of_sigma

    real(dp) function dotv( a, b )
        real, intent(in) :: a(:,:,:), b(:,:,:)
        dotv = sum(real(a,dp) * real(b,dp))
    end function dotv

    real(dp) function normv( a )
        real, intent(in) :: a(:,:,:)
        normv = sqrt(dotv(a, a))
    end function normv

    !> the prior on the BOX fixture with kernels of width sigma_vox, optional column scales
    subroutine new_fixture_prior( prior, sigma_vox, amps )
        class(pcg_atomic_prior), intent(inout) :: prior
        real,                    intent(in)    :: sigma_vox
        real,          optional, intent(in)    :: amps(:)
        real, allocatable :: pos(:,:)
        call build_atoms(NATOMS, ATOM_RAD, ATOM_DMIN, BOX, pos)
        if( present(amps) )then
            call prior%new(BOX, SMPD, pos, sphere_support(BOX, SUPPORT_RAD), amps)
        else
            call prior%new(BOX, SMPD, pos, sphere_support(BOX, SUPPORT_RAD))
        endif
        call prior%build(bfac_of_sigma(sigma_vox, SMPD))
    end subroutine new_fixture_prior

    !> deterministic amplitudes in [lo, hi]
    function varied( n, lo, hi ) result( c )
        integer, intent(in) :: n
        real,    intent(in) :: lo, hi
        real(dp), allocatable :: c(:)
        integer :: i
        allocate(c(n))
        do i = 1, n
            c(i) = real(lo,dp) + real(hi-lo,dp) * modulo(real(i,dp) * 0.6180339887498949_dp, 1._dp)
        end do
    end function varied

    ! ---- the prior on its own ----------------------------------------------------

    !> P_G is a projector: P_G G = G, P_G^2 = P_G, symmetric
    subroutine test_projector()
        class(pcg_atomic_prior), allocatable :: prior
        real, allocatable :: gc(:,:,:), pgc(:,:,:), a(:,:,:), b(:,:,:), pa(:,:,:), pb(:,:,:), ppa(:,:,:)
        real(dp) :: err_span, err_idem, err_sym
        write(*,'(A)') 'test_projector'
        allocate(prior)
        call new_fixture_prior(prior, 1.2)
        call assert_true(prior%get_natoms() == NATOMS, 'every fixture atom is kept')
        gc  = prior%scatter(varied(prior%get_natoms(), 0.5, 1.5))
        pgc = prior%project(gc)
        err_span = normv(pgc - gc) / normv(gc)
        a   = probe_map(BOX, 1)
        b   = probe_map(BOX, 2)
        pa  = prior%project(a)
        pb  = prior%project(b)
        ppa = prior%project(pa)
        err_idem = normv(ppa - pa) / normv(pa)
        err_sym  = abs(dotv(a, pb) - dotv(pa, b)) / (normv(a) * normv(pb))
        write(*,'(A,3ES11.3)') '  |P G c - G c|, |P P a - P a|, <a,Pb>-<Pa,b> (relative): ', err_span, err_idem, err_sym
        call assert_true(err_span <= PROJ_TOL, 'P_G reproduces the span of the kernels')
        call assert_true(err_idem <= PROJ_TOL, 'P_G is idempotent')
        call assert_true(err_sym  <= SYM_TOL,  'P_G is symmetric')
        call prior%kill
    end subroutine test_projector

    !> lambda (I - P_G) is symmetric, positive semi-definite and zero on span(G)
    subroutine test_apply_symmetric()
        real, parameter :: LAMBDA = 2.5
        class(pcg_atomic_prior), allocatable :: prior
        real, allocatable :: a(:,:,:), b(:,:,:), aa(:,:,:), ab(:,:,:), gc(:,:,:), agc(:,:,:)
        real(dp) :: err_sym, quad, err_null
        write(*,'(A)') 'test_apply_symmetric'
        allocate(prior)
        call new_fixture_prior(prior, 1.2)
        a = probe_map(BOX, 3)
        b = probe_map(BOX, 4)
        allocate(aa(BOX,BOX,BOX), ab(BOX,BOX,BOX), source=0.)
        call prior%apply(LAMBDA, a, aa)
        call prior%apply(LAMBDA, b, ab)
        err_sym = abs(dotv(a, ab) - dotv(aa, b)) / (normv(a) * normv(ab))
        quad    = dotv(a, aa) / (normv(a) * normv(aa))
        gc = prior%scatter(varied(prior%get_natoms(), 0.5, 1.5))
        allocate(agc(BOX,BOX,BOX), source=0.)
        call prior%apply(LAMBDA, gc, agc)
        err_null = normv(agc) / (LAMBDA * normv(gc))
        write(*,'(A,3ES11.3)') '  asymmetry, <a,Aa>/(|a||Aa|), |A G c|/(lambda |G c|): ', err_sym, quad, err_null
        call assert_true(err_sym  <= SYM_TOL,   'the prior term is symmetric')
        call assert_true(quad     >= -SYM_TOL,  'the prior term is positive semi-definite')
        call assert_true(err_null <= PROJ_TOL,  'the prior term vanishes on the atom kernels')
        call prior%kill
    end subroutine test_apply_symmetric

    !> P_G does not change when the columns of G are rescaled
    subroutine test_column_scaling()
        class(pcg_atomic_prior), allocatable :: p1, p2
        real(dp), allocatable :: amps(:)
        real, allocatable :: x(:,:,:), px1(:,:,:), px2(:,:,:)
        real(dp) :: err
        write(*,'(A)') 'test_column_scaling'
        allocate(p1, p2)
        call new_fixture_prior(p1, 1.2)
        amps = varied(NATOMS, 0.1, 10.)
        call new_fixture_prior(p2, 1.2, real(amps))
        x   = probe_map(BOX, 5)
        px1 = p1%project(x)
        px2 = p2%project(x)
        err = normv(px2 - px1) / normv(px1)
        write(*,'(A,ES11.3)') '  |P_G x - P_G'' x| / |P_G x| with columns scaled 0.1 to 10: ', err
        call assert_true(err <= PROJ_TOL, 'P_G is invariant to column scaling')
        call p1%kill
        call p2%kill
    end subroutine test_column_scaling

    !> fitted B shapes the kernels, while occupancy times model intensity fixes their relative weights up to one scale
    subroutine test_fitted_model_projector()
        class(pcg_atomic_prior), allocatable :: p1, p2, p3, p4
        type(pcg_atomic_prior_stats) :: stats
        real, allocatable :: pos(:,:), supp(:,:,:), q(:), q2(:), inten(:), inten2(:), b(:), b2(:)
        real, allocatable :: model(:,:,:), pmodel(:,:,:), probe(:,:,:), pp1(:,:,:), pp2(:,:,:), pp3(:,:,:), pp4(:,:,:),&
            &amodel(:,:,:)
        real(dp), allocatable :: coeff(:)
        real(dp) :: err_model, err_q, err_i, err_b, err_null
        integer :: i
        write(*,'(A)') 'test_fitted_model_projector'
        allocate(p1, p2, p3, p4)
        call build_atoms(NATOMS, ATOM_RAD, ATOM_DMIN, BOX, pos)
        supp = sphere_support(BOX, SUPPORT_RAD)
        allocate(q(NATOMS), q2(NATOMS), inten(NATOMS), inten2(NATOMS), b(NATOMS), b2(NATOMS))
        do i = 1, NATOMS
            q(i)     = 0.35 + 0.65 * modulo(real(i) * 0.381966, 1.)
            inten(i) = 0.6 + 1.8 * modulo(real(i) * 0.618034, 1.)
            b(i)     = bfac_of_sigma(1.0 + 0.7 * modulo(real(i) * 0.414214, 1.), SMPD)
            b2(i)    = bfac_of_sigma(1.0 + 0.7 * modulo(real(i) * 0.732051, 1.), SMPD)
        enddo
        q2     = q(NATOMS:1:-1)
        inten2 = inten(NATOMS:1:-1)
        call p1%new(BOX, SMPD, pos, supp, bfacs=b,  occupancies=q,  intensities=inten)
        call p2%new(BOX, SMPD, pos, supp, bfacs=b,  occupancies=q2, intensities=inten)
        call p3%new(BOX, SMPD, pos, supp, bfacs=b,  occupancies=q,  intensities=inten2)
        call p4%new(BOX, SMPD, pos, supp, bfacs=b2, occupancies=q,  intensities=inten)
        coeff = real(q * inten * (4. * PI / b)**1.5,dp)
        model = p1%scatter(coeff)
        pmodel = p1%project(model)
        err_model = normv(pmodel - model) / normv(model)
        allocate(amodel(BOX,BOX,BOX), source=0.)
        call p1%apply(2.5, model, amodel)
        err_null = normv(amodel) / (2.5_dp * normv(model))
        probe = probe_map(BOX, 11)
        pp1 = p1%project(probe)
        pp2 = p2%project(probe)
        pp3 = p3%project(probe)
        pp4 = p4%project(probe)
        err_q = normv(pp2 - pp1) / normv(pp1)
        err_i = normv(pp3 - pp1) / normv(pp1)
        err_b = normv(pp4 - pp1) / normv(pp1)
        stats = p1%get_stats()
        write(*,'(A,5ES11.3)') '  model error, null error, q effect, I effect, B effect: ', &
            &err_model, err_null, err_q, err_i, err_b
        call assert_true(p1%uses_fitted_model(), 'the fitted-model mode is active')
        call assert_true(stats%bfac_model_max > stats%bfac_model_min, 'the fitted-model B range is retained')
        call assert_true(err_model <= PROJ_TOL, 'the fitted model is in the active prior subspace')
        call assert_true(err_null <= PROJ_TOL, 'the prior vanishes on a scaled fitted model')
        call assert_true(err_q  > 1.e-3_dp, 'occupancy changes the active prior subspace')
        call assert_true(err_i  > 1.e-3_dp, 'model intensity changes the active prior subspace')
        call assert_true(err_b  > 1.e-3_dp, 'per-atom B factors change the active prior subspace')
        call p1%kill
        call p2%kill
        call p3%kill
        call p4%kill
    end subroutine test_fitted_model_projector

    !> the Gram solve returns a known solution of G^T G c = y
    subroutine test_gram_solve()
        class(pcg_atomic_prior), allocatable :: prior
        real(dp), allocatable :: c_true(:), c(:), y(:)
        real(dp) :: err, rres
        integer  :: nit
        write(*,'(A)') 'test_gram_solve'
        allocate(prior)
        call new_fixture_prior(prior, 1.6)
        c_true = varied(prior%get_natoms(), -1., 2.)
        y = prior%gather(prior%scatter(c_true))
        call prior%solve_gram(y, c, nit, rres)
        err = sqrt(sum((c - c_true)**2) / sum(c_true**2))
        write(*,'(A,I0,A,2ES11.3)') '  iterations ', nit, ', residual and solution error: ', rres, err
        call assert_true(err <= GRAM_TOL, 'the Gram solve is accurate')
        call prior%kill
    end subroutine test_gram_solve

    !> a noise-free map of wider, untruncated atoms: the fit returns their B less the narrowest kernel's
    subroutine test_width_fit()
        real, parameter :: SIGMA_TRUE = 1.6
        class(pcg_atomic_prior), allocatable :: prior
        real, allocatable :: pos(:,:), x(:,:,:), supp(:,:,:)
        real(dp), allocatable :: c(:)
        real    :: bdata, bexpect
        integer :: nsh, i, j, k, ia
        write(*,'(A)') 'test_width_fit'
        allocate(prior)
        call build_atoms(NATOMS, ATOM_RAD, ATOM_DMIN, BOX, pos)
        c    = varied(size(pos,2), 0.5, 1.5)
        supp = sphere_support(BOX, SUPPORT_RAD)
        allocate(x(BOX,BOX,BOX), source=0.)
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    if( supp(i,j,k) <= 0. ) cycle
                    do ia = 1, size(pos,2)
                        x(i,j,k) = x(i,j,k) + real(c(ia)) * exp(-sum(([real(i),real(j),real(k)] - pos(:,ia))**2) / &
                            &(2. * SIGMA_TRUE**2))
                    end do
                end do
            end do
        end do
        call prior%new(BOX, SMPD, pos, supp)
        call prior%fit_width(x, bdata, nsh)
        bexpect = bfac_of_sigma(SIGMA_TRUE, SMPD) - prior%get_bfac_min()
        write(*,'(A,I0,A,2F9.3)') '  shells ', nsh, ', fitted and true B_data (A^2): ', bdata, bexpect
        call assert_true(nsh >= 3, 'the width fit has shells to fit')
        call assert_true(abs(bdata - bexpect) <= WIDTH_TOL * bexpect, 'the width fit recovers B_data')
        call prior%kill
    end subroutine test_width_fit

    ! ---- the hook in the reconstructor -------------------------------------------

    !> exact line integrals of Gaussian atoms of width R_SIGMA_VOX (voxel coordinates, 1-based), amplitudes c
    subroutine project_atoms( pos, c, rot, img )
        real,     intent(in)  :: pos(:,:), rot(3,3)
        real(dp), intent(in)  :: c(:)
        real,     intent(out) :: img(RBOX,RBOX)
        real, allocatable :: xy(:,:), p0(:,:)
        real    :: ctr, amp
        integer :: i, j, ia
        ctr = real(RBOX/2 + 1)
        amp = sqrt(2. * PI) * R_SIGMA_VOX
        p0  = pos - ctr
        allocate(xy(2,size(pos,2)))
        xy(1,:) = matmul(rot(1,:), p0)
        xy(2,:) = matmul(rot(2,:), p0)
        img = 0.
        do j = 1, RBOX
            do i = 1, RBOX
                do ia = 1, size(pos,2)
                    img(i,j) = img(i,j) + real(c(ia)) * amp * exp(-((real(i) - ctr - xy(1,ia))**2 + &
                        &(real(j) - ctr - xy(2,ia))**2) / (2. * R_SIGMA_VOX**2))
                end do
            end do
        end do
    end subroutine project_atoms

    !> an accumulated reconstructor on the RBOX fixture (kernel built, sphere support) from views of kind
    !! tilt (single-axis, -60 to +60) or dense (spiral), and the fixture's atom positions
    subroutine new_fixture_reconstructor( kind, pcgop, pos, c_true )
        character(len=*),         intent(in)    :: kind
        class(reconstructor_pcg), intent(inout) :: pcgop
        real, allocatable,        intent(out)   :: pos(:,:)
        real(dp), allocatable,    intent(out)   :: c_true(:)
        type(oris)  :: os
        type(image) :: proj
        complex, allocatable :: planes(:,:,:)
        real    :: img(RBOX,RBOX,1)
        integer :: i, n, lims2(2,2)
        call build_atoms(R_NATOMS, R_ATOM_RAD, ATOM_DMIN, RBOX, pos)
        c_true = varied(size(pos,2), 0.5, 1.5)
        select case( kind )
            case( 'tilt' )
                n = 31
                call os%new(n, .false.)
                do i = 1, n
                    call os%set_euler(i, [0., -60. + real(i-1) * 4., 0.])
                end do
            case default
                n = 120
                call os%new(n, .false.)
                call os%spiral()
        end select
        call pcgop%new(RBOX, RSMPD, PCG_LAMBDA)
        lims2 = pcgop%get_lims2()
        allocate(planes(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),n))
        call proj%new([RBOX,RBOX,1], RSMPD, wthreads=.false.)
        do i = 1, n
            call project_atoms(pos, c_true, os%get_mat(i), img(:,:,1))
            call proj%set_rmat(img, .false.)
            call proj%fft()
            planes(:,:,i) = pcgop%extract_native_plane(proj)
        end do
        call proj%kill()
        call pcgop%set_window_sphere(R_SUPPORT)
        call pcgop%prep_particles(os, use_ctf=.false.)
        call pcgop%begin_accum
        call pcgop%accumulate_batch(planes, n, 1)
        call pcgop%end_accum(.true.)
        call pcgop%set_op_mode(PCG_OP_KERNEL)
        call os%kill()
    end subroutine new_fixture_reconstructor

    subroutine new_fixture_reconstructor_prior( pcgop, pos, prior )
        class(reconstructor_pcg), intent(in)    :: pcgop
        real,                     intent(in)    :: pos(:,:)
        class(pcg_atomic_prior),  intent(inout) :: prior
        real, allocatable :: supp(:,:,:)
        call pcgop%get_solve_support(supp)
        call prior%new(RBOX, RSMPD, pos, supp)
        call prior%build(bfac_of_sigma(R_SIGMA_VOX, RSMPD))
    end subroutine new_fixture_reconstructor_prior

    !> strength zero leaves the solve bit for bit as without the prior
    subroutine test_reconstructor_zero_strength()
        class(reconstructor_pcg), allocatable :: pcgop
        class(pcg_atomic_prior),  allocatable, target :: prior
        real,     allocatable :: pos(:,:), x0(:,:,:), x1(:,:,:)
        real(dp), allocatable :: c_true(:)
        write(*,'(A)') 'test_reconstructor_zero_strength'
        allocate(pcgop, prior)
        call new_fixture_reconstructor('tilt', pcgop, pos, c_true)
        allocate(x0(RBOX,RBOX,RBOX), x1(RBOX,RBOX,RBOX), source=0.)
        call pcgop%solve_accum(x0, maxits=6, rtol=0.)
        call new_fixture_reconstructor_prior(pcgop, pos, prior)
        call pcgop%set_atomic_prior(prior, 0.)
        call pcgop%solve_accum(x1, maxits=6, rtol=0.)
        call assert_true(pcgop%get_atomic_lambda() == 0., 'strength zero attaches no prior')
        call assert_true(all(x1 == x0), 'strength zero reproduces the prior-free solve bit for bit')
        call pcgop%kill
        call prior%kill
    end subroutine test_reconstructor_zero_strength

    !> with the prior on, the kernel and matrix-free operators agree as without it, and differ by the same amount
    subroutine test_reconstructor_operator_modes()
        class(reconstructor_pcg), allocatable :: pcgop
        class(pcg_atomic_prior),  allocatable, target :: prior
        real,     allocatable :: pos(:,:), probe(:,:,:), hm(:,:,:), hk(:,:,:), hm0(:,:,:), hk0(:,:,:)
        real(dp), allocatable :: c_true(:)
        real(dp) :: err_mode, err_kern
        write(*,'(A)') 'test_reconstructor_operator_modes'
        allocate(pcgop, prior)
        call new_fixture_reconstructor('tilt', pcgop, pos, c_true)
        call new_fixture_reconstructor_prior(pcgop, pos, prior)
        probe = probe_map(RBOX, 6)
        call pcgop%set_op_mode(PCG_OP_MATRIXFREE)
        hm0 = pcgop%apply_normal(probe)
        call pcgop%set_op_mode(PCG_OP_KERNEL)
        hk0 = pcgop%apply_normal(probe)
        call pcgop%set_atomic_prior(prior, 1.)
        call assert_true(pcgop%get_atomic_lambda() > 0., 'a positive strength attaches the prior')
        call pcgop%set_op_mode(PCG_OP_MATRIXFREE)
        hm = pcgop%apply_normal(probe)
        call pcgop%set_op_mode(PCG_OP_KERNEL)
        hk = pcgop%apply_normal(probe)
        err_mode = normv((hk - hm) - (hk0 - hm0)) / normv(hm)
        err_kern = normv(hk - hm) / normv(hm)
        write(*,'(A,2ES11.3)') '  mode difference change, kernel vs matrix-free (relative): ', err_mode, err_kern
        call assert_true(err_mode <= MODE_TOL,     'the prior term is the same in both operator modes')
        call assert_true(err_kern <= INTERIOR_TOL, 'kernel and matrix-free operators agree with the prior on')
        call pcgop%kill
        call prior%kill
    end subroutine test_reconstructor_operator_modes

    !> amplitudes fitted through the normal operator track the true amplitudes on dense noise-free views
    subroutine test_reconstructor_amplitudes()
        class(reconstructor_pcg), allocatable :: pcgop
        class(pcg_atomic_prior),  allocatable :: prior
        real,     allocatable :: pos(:,:)
        real(dp), allocatable :: c_true(:), c(:)
        real(dp) :: corr
        real     :: rres
        integer  :: nit
        write(*,'(A)') 'test_reconstructor_amplitudes'
        allocate(pcgop, prior)
        call new_fixture_reconstructor('dense', pcgop, pos, c_true)
        call new_fixture_reconstructor_prior(pcgop, pos, prior)
        call pcgop%fit_atom_amplitudes(prior, c, 50, 1.e-4, nit, rres)
        corr = sum((c - sum(c)/size(c)) * (c_true - sum(c_true)/size(c_true))) / &
            &sqrt(sum((c - sum(c)/size(c))**2) * sum((c_true - sum(c_true)/size(c_true))**2))
        write(*,'(A,I0,A,ES11.3,A,F8.5,A,F8.4)') '  iterations ', nit, ', residual ', rres, &
            &', correlation with the truth ', corr, ', mean amplitude ', sum(c)/size(c)
        call assert_true(all(ieee_is_finite(c)), 'the fitted amplitudes are finite')
        call assert_true(corr >= AMP_CORR_MIN, 'amplitudes fitted to the data track the true amplitudes')
        call pcgop%kill
        call prior%kill
    end subroutine test_reconstructor_amplitudes

end module simple_pcg_atomic_prior_tester
