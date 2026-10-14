!@descr: white-box tests of the cls_expansion kernel route (submodule of simple_flex_cls_kernel): coefficient split, ML latents, representatives, Stein correction, the class route end to end
! Synthetic classes on the half-plane of a 32-pixel box: a mean plus one latent direction, a CTF-like factor per
! member and coefficient, complex Gaussian noise of unit power (w = 1, wq = 2), so the model's likelihood is
! calibrated. Expected values are closed forms or properties with wide margins.
submodule (simple_flex_cls_kernel) simple_flex_cls_kernel_tester
use simple_rnd,                only: gasdev, ran3
use simple_flex_cls_expansion, only: flex_cls_ctf_mean, flex_cls_weighted_mean
use simple_test_utils,         only: assert_true, assert_false, assert_int, assert_real
implicit none

integer, parameter :: TBOX  = 32
integer, parameter :: KFIT  = 8         !< fit band (shell)
integer, parameter :: SEED0 = 20261013

contains

    module subroutine run_all_flex_cls_kernel_tests()
        write(*,'(A)') '**** running all flex cls_expansion kernel tests ****'
        call test_coeff_sets()
        call test_ml_latents_exact()
        call test_representatives()
        call test_stein_correction()
        call test_kernel_null()
        call test_kernel_two_states()
        call test_kernel_budget()
    end subroutine run_all_flex_cls_kernel_tests

    !> the half-plane lattice of the box: shells 1..nyq, the h = 0 column for k >= 0 only
    subroutine lattice( hidx, kidx, shell, wq )
        integer,  allocatable, intent(out) :: hidx(:), kidx(:), shell(:)
        real(sp), allocatable, intent(out) :: wq(:)
        integer :: h, k, sh, n, pass
        do pass = 1, 2
            n = 0
            do h = 0, TBOX / 2
                do k = -TBOX / 2, TBOX / 2 - 1
                    if( h == 0 .and. k < 0 ) cycle
                    sh = nint(sqrt(real(h*h + k*k)))
                    if( sh < 1 .or. sh > TBOX / 2 - 1 ) cycle
                    n = n + 1
                    if( pass == 2 )then
                        hidx(n) = h; kidx(n) = k; shell(n) = sh; wq(n) = 2.
                    endif
                end do
            end do
            if( pass == 1 ) allocate(hidx(n), kidx(n), shell(n), wq(n))
        end do
    end subroutine lattice

    !> complex Gaussian of unit power
    complex(sp) function cgauss()
        cgauss = cmplx(gasdev(), gasdev(), kind=sp) / sqrt(2.)
    end function cgauss

    !> a class of n members, y_i = c_i (mu + z_i v) + noise; ainfo the latent information of v over the fit band at c = 1
    subroutine make_class( n, ztrue, ainfo, seed, y, c, w, mu, v, fitmask )
        integer,                  intent(in)  :: n, seed
        real(dp),                 intent(in)  :: ztrue(:), ainfo
        complex(sp), allocatable, intent(out) :: y(:,:), mu(:), v(:)
        real(sp),    allocatable, intent(out) :: c(:,:), w(:,:)
        logical,     allocatable, intent(out) :: fitmask(:)
        integer,  allocatable :: hidx(:), kidx(:), shell(:)
        real(sp), allocatable :: wq(:)
        real(sp) :: defoc
        integer  :: nc, i, j
        call lattice(hidx, kidx, shell, wq)
        nc = size(hidx)
        call seed_rnd_fixed(seed)
        allocate(mu(nc), v(nc), y(nc,n), c(nc,n), w(nc,n), fitmask(nc))
        fitmask = shell <= KFIT
        do j = 1, nc
            mu(j) = 2. * exp(-real(shell(j)) / 5.) * cgauss()
            v(j)  = exp(-real(shell(j)) / 5.) * cgauss()
        end do
        v = v * real(sqrt(ainfo / real(sum(wq * abs(v)**2, mask=fitmask),dp)), sp)
        w = 1.
        do i = 1, n
            defoc = ran3()
            do j = 1, nc
                c(j,i) = 0.4 + 0.6 * abs(cos(0.2 * real(shell(j)) * (1. + defoc)))
                y(j,i) = c(j,i) * (mu(j) + real(ztrue(i), sp) * v(j)) + cgauss()
            end do
        end do
    end subroutine make_class

    !> a model by hand: mean and basis columns on the fit band, no fit
    subroutine hand_model( model, mu, u, fitmask, nptcls )
        type(flex_cls_model), intent(inout) :: model
        complex(sp),          intent(in)    :: mu(:), u(:,:)
        logical,              intent(in)    :: fitmask(:)
        integer,              intent(in)    :: nptcls
        integer :: j
        call model%new(count(fitmask), size(u,2), nptcls)
        model%fitidx = pack([(j, j=1,size(fitmask))], fitmask)
        model%mu     = cmplx(mu(model%fitidx), kind=dp)
        model%U      = cmplx(u(model%fitidx,:), kind=dp)
    end subroutine hand_model

    !> the class route as cls_expansion runs it: a rank-2 model with the class mean as nuisance, new and select
    subroutine split_class( y, c, w, nsub, kern, model )
        complex(sp),           intent(in)    :: y(:,:)
        real(sp),              intent(in)    :: c(:,:), w(:,:)
        integer,               intent(in)    :: nsub
        type(flex_cls_kernel), intent(inout) :: kern
        type(flex_cls_model),  intent(inout) :: model
        integer,     allocatable :: hidx(:), kidx(:), shell(:)
        real(sp),    allocatable :: wq(:)
        complex(dp), allocatable :: mu(:)
        complex(sp), allocatable :: nuis(:,:)
        call lattice(hidx, kidx, shell, wq)
        allocate(mu(size(y,1)), nuis(size(y,1),1))
        call flex_cls_weighted_mean(y, c, w, mu)
        nuis(:,1) = cmplx(mu, kind=sp)
        call seed_rnd_fixed(SEED0 + 10)
        call flex_cls_fit(model, y, c, w, wq, shell <= KFIT, 2, verbose=.false., nuis=nuis)
        call kern%new(model, nuis, y, c, w, wq, hidx, kidx, SEED0 + 11)
        call kern%select(nsub, 20)
    end subroutine split_class

    !> the checkerboard: disjoint, covering the fit band, balanced, constant on a block
    subroutine test_coeff_sets()
        integer,  allocatable :: hidx(:), kidx(:), shell(:)
        real(sp), allocatable :: wq(:)
        logical,  allocatable :: fitmask(:), seta(:), setb(:)
        integer :: j, l, na, nb
        logical :: l_const
        write(*,'(A)') 'test_coeff_sets'
        call lattice(hidx, kidx, shell, wq)
        allocate(fitmask(size(hidx)), seta(size(hidx)), setb(size(hidx)))
        fitmask = shell <= KFIT
        call coeff_sets(hidx, kidx, fitmask, seta, setb)
        na = count(seta); nb = count(setb)
        call assert_false(any(seta .and. setb), 'the two sets are disjoint')
        call assert_true(all((seta .or. setb) .eqv. fitmask), 'the two sets cover the fit band and nothing else')
        call assert_true(abs(na - nb) <= (na + nb) / 5, 'the two sets are balanced within 20%')
        l_const = .true.
        do j = 1, size(hidx)
            do l = j + 1, size(hidx)
                if( .not. (fitmask(j) .and. fitmask(l)) ) cycle
                if( floor(real(hidx(j)) / BLOCK) /= floor(real(hidx(l)) / BLOCK) ) cycle
                if( floor(real(kidx(j)) / BLOCK) /= floor(real(kidx(l)) / BLOCK) ) cycle
                if( seta(j) .neqv. seta(l) ) l_const = .false.
            end do
        end do
        call assert_true(l_const, 'a Fourier block lies in one set')
    end subroutine test_coeff_sets

    !> noise-free members: the ML latent is the generating latent, its covariance (sum wq w c^2 Re(U^H U))^-1, and a
    !! nuisance column is profiled out without biasing it
    subroutine test_ml_latents_exact()
        integer, parameter :: N = 5
        type(flex_cls_model) :: model
        complex(sp), allocatable :: y(:,:), mu(:), v(:), u(:,:)
        real(sp),    allocatable :: c(:,:), w(:,:), wq(:)
        integer,     allocatable :: hidx(:), kidx(:), shell(:)
        logical,     allocatable :: fitmask(:), info(:)
        real(dp) :: ztrue(N), z(N,2), cov(2,2,N), a(2,2), ainv(2,2), nu(N), z1(N,1), cov1(1,1,N)
        integer  :: i, j, jj
        write(*,'(A)') 'test_ml_latents_exact'
        ztrue = [-1.5d0, -0.3d0, 0.d0, 0.7d0, 2.d0]
        call make_class(N, ztrue, 4.d0, SEED0 + 1, y, c, w, mu, v, fitmask)
        call lattice(hidx, kidx, shell, wq)
        allocate(u(size(mu),2))
        u(:,1) = v
        do j = 1, size(mu)
            u(j,2) = exp(-real(shell(j)) / 4.) * cgauss()
        end do
        do i = 1, N
            y(:,i) = c(:,i) * (mu + real(ztrue(i),sp) * u(:,1) - 0.5 * real(ztrue(i),sp) * u(:,2))
        end do
        call hand_model(model, mu, u, fitmask, N)
        call ml_latents(model, y, c, w, wq, z, cov, info)
        call assert_true(all(info), 'every noise-free member has a latent')
        do i = 1, N
            call assert_real(real(ztrue(i)), real(z(i,1)), 2.e-4, 'ML latent recovers the generating latent (1)')
            call assert_real(real(-0.5d0*ztrue(i)), real(z(i,2)), 2.e-4, 'ML latent recovers the generating latent (2)')
            a = 0.d0
            do jj = 1, model%nfit
                j = model%fitidx(jj)
                a(1,1) = a(1,1) + wq(j) * w(j,i) * c(j,i)**2 * abs(cmplx(u(j,1),kind=dp))**2
                a(2,2) = a(2,2) + wq(j) * w(j,i) * c(j,i)**2 * abs(cmplx(u(j,2),kind=dp))**2
                a(1,2) = a(1,2) + wq(j) * w(j,i) * c(j,i)**2 * real(conjg(cmplx(u(j,1),kind=dp)) * u(j,2), dp)
            end do
            ainv = reshape([a(2,2), -a(1,2), -a(1,2), a(1,1)], [2,2]) / (a(1,1) * a(2,2) - a(1,2)**2)
            call assert_real(1., real(cov(1,1,i) / ainv(1,1)), 1.e-4, 'ML covariance is the inverse information (11)')
            call assert_real(1., real(cov(2,2,i) / ainv(2,2)), 1.e-4, 'ML covariance is the inverse information (22)')
            call assert_real(real(ainv(1,2)), real(cov(1,2,i)), 1.e-4 * real(sqrt(ainv(1,1)*ainv(2,2))), &
                &'ML covariance is the inverse information (12)')
        end do
        call model%kill
        ! one latent column, the second a nuisance with its own coefficients
        nu = [0.4d0, -1.d0, 2.5d0, 0.d0, -0.2d0]
        do i = 1, N
            y(:,i) = c(:,i) * (mu + real(ztrue(i),sp) * u(:,1) + real(nu(i),sp) * u(:,2))
        end do
        call hand_model(model, mu, u(:,1:1), fitmask, N)
        model%nnuis = 1
        model%N     = cmplx(u(model%fitidx,2:2), kind=dp)
        allocate(model%nu(N,1), source=0.d0)
        call ml_latents(model, y, c, w, wq, z1, cov1, info)
        do i = 1, N
            call assert_real(real(ztrue(i)), real(z1(i,1)), 2.e-4, 'the nuisance is profiled out of the ML latent')
        end do
        call model%kill
    end subroutine test_ml_latents_exact

    !> two clusters of latents, a budget of four: four pre-images of members, the first two at the cluster centres
    !! (coverage before redundancy), every member labelled with a representative of its cluster, every one holding
    !! minpop; the budget filled when no pre-image reaches minpop
    subroutine test_representatives()
        integer,  parameter :: NA = 200, NB = 150
        real(dp), parameter :: CA(2) = [-1.d0, 0.d0], CB(2) = [1.2d0, 0.3d0], S2 = 0.02d0
        type(flex_cls_kernel) :: k
        integer, allocatable  :: lab(:)
        real(dp) :: cen(2,4)
        integer  :: i, r, ia, ib
        logical  :: l_ina(4)
        write(*,'(A)') 'test_representatives'
        call seed_rnd_fixed(SEED0 + 2)
        k%n = NA + NB; k%nk = 2; k%nsig = 2; k%tau = 0.1d0; k%exists = .true.
        allocate(k%info(NA+NB), source=.true.)
        allocate(k%zr(NA+NB,2), k%covr(2,2,NA+NB), source=0.d0)
        allocate(k%covc(2,2), source=reshape([S2, 0.d0, 0.d0, S2], [2,2]))
        do i = 1, NA + NB
            k%covr(:,:,i) = k%covc
            if( i <= NA )then
                k%zr(i,:) = CA + sqrt(S2) * [real(gasdev(),dp), real(gasdev(),dp)]
            else
                k%zr(i,:) = CB + sqrt(S2) * [real(gasdev(),dp), real(gasdev(),dp)]
            endif
        end do
        call k%select(4, 20)
        call assert_int(4, k%get_nrep(), 'the budget of four is filled')
        if( k%get_nrep() == 4 )then
            call assert_true(all(k%imem >= 1 .and. k%imem <= NA + NB), 'every representative is a member''s pre-image')
            do r = 1, 4
                cen(:,r) = k%zr(k%get_member(r),:)
                l_ina(r) = cen(1,r) < 0.d0
            end do
            ia = merge(1, 2, l_ina(1))
            ib = 3 - ia
            call assert_true(l_ina(1) .neqv. l_ina(2), 'the first two representatives are in different clusters')
            call assert_true(sqrt(sum((cen(:,ia) - CA)**2)) < 0.15d0 .and. sqrt(sum((cen(:,ib) - CB)**2)) < 0.15d0, &
                &'the first two representatives are the cluster centres')
            call assert_true(l_ina(1), 'the denser cluster comes first')
            call assert_true(all([(k%neff(r) >= 20.d0, r=1,4)]), 'every representative holds minpop effective members')
            lab = k%labels()
            call assert_int(NA + NB, count(l_ina(lab(1:NA))) + count(.not. l_ina(lab(NA+1:))), &
                &'every member is labelled with a representative of its cluster')
        endif
        call k%select(4, 100000)
        call assert_int(4, k%get_nrep(), 'the budget is filled when no pre-image reaches minpop')
        call k%select(1, 20)
        call assert_false(k%l_pre, 'a budget of one gives the class mean')
        call k%kill
    end subroutine test_representatives

    !> A homogeneous class weighted around an off-centre location: the plain kernel mean deviates along the latent
    !! direction by cen cov / (2 cov + covc + tau^2) (tau = 0; 0.5 cen without covc, cen / 3 with covc = cov);
    !! Stein's identity with the kernel's own precision removes it
    subroutine test_stein_correction()
        integer,  parameter :: N = 6000
        real(dp), parameter :: CEN0 = 1.d0
        type(flex_cls_model)  :: model
        type(flex_cls_kernel) :: k
        complex(sp), allocatable :: y(:,:), mu(:), v(:), avg(:)
        real(sp),    allocatable :: c(:,:), w(:,:), wq(:)
        integer,     allocatable :: hidx(:), kidx(:), shell(:)
        logical,     allocatable :: fitmask(:), members(:), okl(:)
        real(dp),    allocatable :: z(:,:), lw(:,:), wts(:)
        real(dp) :: ztrue(N), d_plain, d_stein, expected, unorm, covmean
        integer  :: icase
        write(*,'(A)') 'test_stein_correction'
        ztrue = 0.d0
        call make_class(N, ztrue, 1.d0, SEED0 + 3, y, c, w, mu, v, fitmask)
        call lattice(hidx, kidx, shell, wq)
        call hand_model(model, mu, reshape(v, [size(v),1]), fitmask, N)
        ! the kernel by hand: one latent direction, not whitened, so the location is in latent units
        k%n = N; k%nk = 1; k%nn = 0; k%nsig = 1; k%tau = 0.d0; k%l_pre = .true.; k%exists = .true.
        allocate(z(N,1), k%cov(1,1,N))
        call ml_latents(model, y, c, w, wq, z, k%cov, k%info)
        allocate(k%snz(1,1,N), source=0.d0)
        allocate(k%wsig(1,1), source=1.d0)
        k%zr = z; k%covr = k%cov; k%fitidx = model%fitidx; k%basis = model%U
        covmean = sum(k%cov(1,1,:)) / real(N,dp)
        allocate(members(N), source=.true.)
        allocate(avg(size(mu)), lw(N,1), wts(N))
        unorm = sum(real(wq,dp) * abs(cmplx(v,kind=dp))**2, mask=fitmask)
        do icase = 1, 2
            k%covc = reshape([merge(0.d0, covmean, icase == 1)], [1,1])
            call k%logweights(reshape([CEN0], [1,1]), lw, okl)
            wts = exp(lw(:,1) - maxval(lw(:,1)))
            call k%kernel_mean([CEN0], wts, y, c, w, members, avg, .false.)
            d_plain = latent_dev(avg)
            call k%kernel_mean([CEN0], wts, y, c, w, members, avg, .true.)
            d_stein = latent_dev(avg)
            expected = CEN0 * covmean / (2.d0 * covmean + k%covc(1,1))
            write(*,'(A,I0,A,F8.4,A,F8.4,A,F8.4)') '    case ', icase, ': selection bias (latent units): plain ', d_plain, &
                &' Stein ', d_stein, ' expected plain about ', expected
            call assert_true(d_plain > 0.6d0 * expected, 'the plain kernel mean is biased along the latent direction')
            call assert_true(abs(d_stein) < 0.12d0, 'Stein''s identity removes the selection bias')
        end do
        call model%kill
        call k%kill

      contains

        !> the deviation of an average from the true mean along v, in latent units
        real(dp) function latent_dev( a )
            complex(sp), intent(in) :: a(:)
            latent_dev = sum(real(wq,dp) * real(conjg(cmplx(v,kind=dp)) * (a - mu), dp), mask=fitmask) / unorm
        end function latent_dev

    end subroutine test_stein_correction

    !> a homogeneous class: no significant CV gain over its mean, still the budget of four pre-images at a finite
    !! bandwidth; with a budget of one, the class mean, also as every member's pre-image
    subroutine test_kernel_null()
        integer, parameter :: N = 400
        type(flex_cls_kernel) :: k
        type(flex_cls_model)  :: model
        complex(sp), allocatable :: y(:,:), mu(:), v(:), avg(:), ref(:)
        real(sp),    allocatable :: c(:,:), w(:,:)
        logical,     allocatable :: fitmask(:), members(:)
        real(dp) :: ztrue(N)
        write(*,'(A)') 'test_kernel_null'
        ztrue = 0.d0
        call make_class(N, ztrue, 6.d0, SEED0 + 4, y, c, w, mu, v, fitmask)
        call split_class(y, c, w, 4, k, model)
        write(*,'(A,F8.3,A,F6.2,A,F6.2,A,I0)') '    CV z-score=', k%zscore, ' tau=', k%tau_cv, ' pre-image tau=', k%tau, &
            &' representatives=', k%get_nrep()
        call assert_true(k%zscore < 3.d0, 'a homogeneous class shows no significant CV gain over its mean')
        call assert_int(4, k%get_nrep(), 'the budget of four is filled')
        call assert_true(k%l_pre .and. k%tau >= 0.d0, 'the representatives are pre-images at a finite bandwidth')
        allocate(members(N), source=.true.)
        allocate(avg(size(mu)), ref(size(mu)))
        call flex_cls_ctf_mean(y, c, w, ref)
        call k%select(1, 20)
        call assert_int(1, k%get_nrep(), 'a budget of one gives one representative')
        call k%average(1, y, c, w, members, avg)
        call assert_true(maxval(abs(avg - ref)) <= 1.e-5 * maxval(abs(ref)), 'it is the CTF-corrected class mean')
        call k%preimage(7, y, c, w, members, avg)
        call assert_true(maxval(abs(avg - ref)) <= 1.e-5 * maxval(abs(ref)), 'so is every member''s pre-image')
        call model%kill
        call k%kill
    end subroutine test_kernel_null

    !> two states (60/40 at z = +-1, latent information 6 at c = 1) and a budget of two: a significant CV gain at
    !! rank 1, two pre-images at the CV optimum, hard labels at least 90% pure (the best is about 96%), the averages
    !! on either side of the class mean along the state direction
    subroutine test_kernel_two_states()
        integer, parameter :: N = 400
        type(flex_cls_kernel) :: k
        type(flex_cls_model)  :: model
        complex(sp), allocatable :: y(:,:), mu(:), v(:), avg(:)
        real(sp),    allocatable :: c(:,:), w(:,:), wq(:)
        integer,     allocatable :: hidx(:), kidx(:), shell(:), lab(:)
        logical,     allocatable :: fitmask(:), members(:)
        real(dp) :: ztrue(N), zmean, proj, pmax, pmin, vnorm
        integer  :: r, nplus, ncorrect
        write(*,'(A)') 'test_kernel_two_states'
        call two_states(ztrue)
        zmean = sum(ztrue) / real(N,dp)
        call make_class(N, ztrue, 6.d0, SEED0 + 6, y, c, w, mu, v, fitmask)
        call split_class(y, c, w, 2, k, model)
        write(*,'(A,F8.3,A,I0,A,F6.2,A,I0,A,F6.3)') '    CV z-score=', k%zscore, ' rank=', k%nsig, ' tau=', k%tau_cv, &
            &' representatives=', k%get_nrep(), ' overlap=', k%overlap
        call assert_true(k%zscore >= 3.d0, 'the CV gain over the class mean is significant')
        call assert_true(k%tau_cv >= 0.d0 .and. k%tau == k%tau_cv, 'the pre-images take the finite CV optimum')
        call assert_int(1, k%nsig, 'the CV keeps the one latent direction the states differ in')
        call assert_int(2, k%get_nrep(), 'two representatives')
        if( k%get_nrep() /= 2 ) return
        lab = k%labels()
        ncorrect = 0
        do r = 1, 2
            nplus    = count(lab == r .and. ztrue > 0.d0)
            ncorrect = ncorrect + max(nplus, count(lab == r) - nplus)
        end do
        write(*,'(A,F6.3)') '    label purity ', real(ncorrect) / real(N)
        call assert_true(real(ncorrect) / real(N) >= 0.9, 'the hard labels are at least 90% pure')
        call lattice(hidx, kidx, shell, wq)
        allocate(members(N), source=.true.)
        allocate(avg(size(mu)))
        vnorm = sum(real(wq,dp) * abs(cmplx(v,kind=dp))**2, mask=fitmask)
        pmax = -huge(1.d0); pmin = huge(1.d0)
        do r = 1, 2
            call k%average(r, y, c, w, members, avg)
            proj = sum(real(wq,dp) * real(conjg(cmplx(v,kind=dp)) * (avg - (mu + real(zmean,sp) * v)), dp), &
                &mask=fitmask) / vnorm
            pmax = max(pmax, proj); pmin = min(pmin, proj)
        end do
        write(*,'(A,F7.3,A,F7.3,A,F7.3,A,F7.3,A)') '    state deviations of the averages: ', pmin, ' and ', pmax, &
            &' (true ', -1.d0 - zmean, ' and ', 1.d0 - zmean, ')'
        call assert_true(pmax > 0.5d0 * (1.d0 - zmean), 'one average lies on the side of the + state')
        call assert_true(pmin < 0.5d0 * (-1.d0 - zmean), 'one average lies on the side of the - state')
        call model%kill
        call k%kill
    end subroutine test_kernel_two_states

    !> the same two states and a budget of six: six pre-images, each holding minpop effective members, both states
    !! represented, every representative's average the pre-image of its member
    subroutine test_kernel_budget()
        integer, parameter :: N = 400, NSUB = 6
        type(flex_cls_kernel) :: k
        type(flex_cls_model)  :: model
        complex(sp), allocatable :: y(:,:), mu(:), v(:), avg(:), pre(:)
        real(sp),    allocatable :: c(:,:), w(:,:)
        integer,     allocatable :: lab(:)
        logical,     allocatable :: fitmask(:), members(:)
        real(dp) :: ztrue(N)
        integer  :: r, nplus
        logical  :: l_plus, l_minus, l_same
        write(*,'(A)') 'test_kernel_budget'
        call two_states(ztrue)
        call make_class(N, ztrue, 6.d0, SEED0 + 6, y, c, w, mu, v, fitmask)
        call split_class(y, c, w, NSUB, k, model)
        call assert_int(NSUB, k%get_nrep(), 'the budget of six is filled')
        if( k%get_nrep() /= NSUB ) return
        write(*,'(A,F8.1,A,F8.1,A,F6.3)') '    neff min ', minval([(k%neff(r), r=1,NSUB)]), ' max ', &
            &maxval([(k%neff(r), r=1,NSUB)]), ' overlap ', k%overlap
        call assert_true(all([(k%neff(r) >= 20.d0, r=1,NSUB)]), 'every representative holds minpop effective members')
        lab = k%labels()
        l_plus = .false.; l_minus = .false.
        do r = 1, NSUB
            nplus = count(lab == r .and. ztrue > 0.d0)
            if( count(lab == r) == 0 ) cycle
            if( 2 * nplus > count(lab == r) )then
                l_plus = .true.
            else
                l_minus = .true.
            endif
        end do
        call assert_true(l_plus .and. l_minus, 'both states are represented')
        allocate(members(N), source=.true.)
        allocate(avg(size(mu)), pre(size(mu)))
        l_same = .true.
        do r = 1, NSUB
            call k%average(r, y, c, w, members, avg)
            call k%preimage(k%get_member(r), y, c, w, members, pre)
            if( maxval(abs(avg - pre)) > 1.e-5 * maxval(abs(avg)) ) l_same = .false.
        end do
        call assert_true(l_same, 'a representative''s average is its member''s pre-image')
        call model%kill
        call k%kill
    end subroutine test_kernel_budget

    !> 60/40 at z = +1/-1
    subroutine two_states( ztrue )
        real(dp), intent(out) :: ztrue(:)
        integer :: i
        call seed_rnd_fixed(SEED0 + 5)
        do i = 1, size(ztrue)
            ztrue(i) = merge(1.d0, -1.d0, ran3() < 0.6)
        end do
    end subroutine two_states

end submodule simple_flex_cls_kernel_tester
