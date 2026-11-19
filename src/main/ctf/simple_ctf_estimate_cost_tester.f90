!@descr: unit tests for the continuous CTF fit's cost (ctf_estimate_cost4Dcont): curvature oracle and both solvers
! A synthetic |CTF| spectrum on the fit's layout (h, k from -box/2, the resolution annulus without the central
! cross): at a zero residual (data proportional to the model) the Levenberg-Marquardt curvature is the Hessian;
! from an isotropic start the defocus improves and the angle stays; L-BFGS-B and 'lm' fit the same starts.
submodule (simple_ctf_estimate_cost) simple_ctf_estimate_cost_tester
use simple_test_utils
use simple_memoize_ft_maps, only: memoize_ft_maps
use simple_linalg,          only: cholesky
implicit none

integer,  parameter :: CBOX = 128
real,     parameter :: CSMPD = 1.5, CKV = 300., CCS = 2.7, CFRACA = 0.1, CHP = 30., CLP = 5.
real,     parameter :: TRUE_DFX = 2.0, TRUE_DFY = 2.3, TRUE_ANGAST = 27., CASTIGTOL = 0.05
real(dp), parameter :: ORACLE_RTOL = 1.d-3   !< guessed: the ftexp gradient check's relative tolerance
real,     parameter :: ANGLE_KEEP_TOL = 1.e-4 !< degrees, guessed: single-precision radians-degrees round trip
! starts inside the cost's capture basin, as the grid search leaves them: at CLP the phase moves about 25 rad per um
real,     parameter :: ISO_OFFSET     = 0.02  !< um

contains

    module subroutine run_all_ctf_estimate_cost_tests()
        write(*,'(A)') '**** running all continuous CTF cost tests ****'
        call test_cont_curvature_oracle()
        call test_cont_isotropic_start()
        call test_cont_solver_comparison()
    end subroutine run_all_ctf_estimate_cost_tests

    !> the fit's pixel set and a synthetic spectrum of the model at (dfx, dfy, angast, phshift): zero_residual gives
    !! the model scaled to the norm sqrt(n), else the fit's normalisation (zero mean, unit variance over the pixels)
    subroutine build_spectrum( pspec, inds, dfx, dfy, angast, phshift, zero_residual )
        class(image),         intent(inout) :: pspec
        integer, allocatable, intent(inout) :: inds(:,:)
        real,                 intent(in)    :: dfx, dfy, angast, phshift
        logical,              intent(in)    :: zero_residual
        real, pointer :: prmat(:,:,:)
        real(dp), allocatable :: vals(:)
        real(dp) :: wl, cs, amp, df_eff, phi, s2, ang, mean, sdev
        integer  :: h, k, l, n, sh, kfrom, kto
        call memoize_ft_maps([CBOX,CBOX,1], CSMPD)
        kfrom = int(real(CBOX)*CSMPD/CHP)
        kto   = min(int(real(CBOX)*CSMPD/CLP), CBOX/2 - 1)
        n = 0
        do h = -CBOX/2, CBOX/2
            do k = -CBOX/2, CBOX/2 - 1
                if( in_fit(h, k) ) n = n + 1
            enddo
        enddo
        if( allocated(inds) ) deallocate(inds)
        allocate(inds(2,n), vals(n))
        wl  = 12.26d0 / dsqrt(1.d3*real(CKV,dp) + 0.9784d0*real(CKV,dp)**2) / real(CSMPD,dp)
        cs  = real(CCS,dp)*1.d7 / real(CSMPD,dp)
        amp = datan(real(CFRACA,dp) / dsqrt(1.d0 - real(CFRACA,dp)**2))
        l = 0
        do h = -CBOX/2, CBOX/2
            do k = -CBOX/2, CBOX/2 - 1
                if( .not. in_fit(h, k) ) cycle
                l = l + 1
                inds(:,l) = [h, k]
                s2     = real(ft_map_spafreqsq(h, k), dp)
                ang    = real(ft_map_astigang(h, k), dp)
                df_eff = 0.5d0*(real(dfx+dfy,dp) + dcos(2.d0*(ang - real(angast,dp)*DPI/180.d0))*real(dfx-dfy,dp)) &
                    &* 1.d4/real(CSMPD,dp)
                phi    = DPI*wl*s2*(df_eff - 0.5d0*wl*wl*s2*cs) + real(phshift,dp) + amp
                vals(l) = dabs(dsin(phi))
            enddo
        enddo
        if( zero_residual )then
            vals = vals * dsqrt(real(n,dp)) / dsqrt(sum(vals**2))
        else
            mean = sum(vals) / real(n,dp)
            sdev = dsqrt(sum((vals - mean)**2) / real(n,dp))
            vals = (vals - mean) / sdev
        endif
        call pspec%new([CBOX,CBOX,1], CSMPD)
        call pspec%get_rmat_ptr(prmat)
        prmat = 0.
        do l = 1, n
            prmat(inds(1,l) + CBOX/2 + 1, inds(2,l) + CBOX/2 + 1, 1) = real(vals(l))
        enddo

    contains

        logical function in_fit( h, k )
            integer, intent(in) :: h, k
            real :: r
            r      = sqrt(real(h*h + k*k))
            sh     = nint(r)
            in_fit = h /= 0 .and. k /= 0 .and. r >= real(kfrom) .and. r <= real(int(real(CBOX)*CSMPD/CLP)) .and. &
                &sh >= kfrom .and. sh <= kto
        end function in_fit

    end subroutine build_spectrum

    subroutine set_parms( parms, dfx, dfy, angast, phshift, fit_phshift )
        type(ctfparams), intent(out) :: parms
        real,            intent(in)  :: dfx, dfy, angast, phshift
        logical,         intent(in)  :: fit_phshift
        parms%smpd    = CSMPD
        parms%kv      = CKV
        parms%cs      = CCS
        parms%fraca   = CFRACA
        parms%dfx     = dfx
        parms%dfy     = dfy
        parms%angast  = angast
        parms%phshift = phshift
        parms%l_fit_phshift = fit_phshift
    end subroutine set_parms

    !> the fit's box: defocus 0.5-5 um, the angle within 30 degrees of its start, the phase shift within pi/6
    subroutine fit_limits( ndim, angast, phshift, limits )
        integer, intent(in)  :: ndim
        real,    intent(in)  :: angast, phshift
        real,    intent(out) :: limits(2,ndim)
        limits(:,1) = [0.5, 5.]
        limits(:,2) = [0.5, 5.]
        if( ndim >= 3 ) limits(:,3) = [angast - 30., angast + 30.]
        if( ndim == 4 ) limits(:,4) = [phshift - real(PI)/6., phshift + real(PI)/6.]
    end subroutine fit_limits

    !> at a zero residual (3 and 4 parameters) the curvature equals the central-difference Hessian of the
    !! analytic gradient within ORACLE_RTOL, and it is symmetric positive semi-definite
    subroutine test_cont_curvature_oracle()
        class(image), allocatable, target :: pspec
        type(ctf_estimate_cost4Dcont) :: cost
        type(ctfparams), target       :: parms
        integer, allocatable, target :: inds(:,:)
        real(dp), allocatable :: vec(:), g(:), gp(:), gm(:), curv(:,:), hfd(:,:), chol(:,:), ridge(:,:)
        real(dp) :: f, err, step
        real     :: limits(2,4), phshift
        integer  :: ndim, d, j
        logical  :: ok
        write(*,'(A)') 'test_cont_curvature_oracle'
        allocate(pspec)
        do ndim = 3, 4
            phshift = merge(0.6, 0., ndim == 4)
            call build_spectrum(pspec, inds, TRUE_DFX, TRUE_DFY, TRUE_ANGAST, phshift, .true.)
            call set_parms(parms, TRUE_DFX, TRUE_DFY, TRUE_ANGAST, phshift, ndim == 4)
            call fit_limits(ndim, TRUE_ANGAST, phshift, limits(:,1:ndim))
            call cost%init(pspec, parms, inds, ndim, limits(:,1:ndim), CASTIGTOL, solver='lm')
            allocate(vec(ndim), g(ndim), gp(ndim), gm(ndim), curv(ndim,ndim), hfd(ndim,ndim), chol(ndim,ndim), &
                &ridge(ndim,ndim))
            vec(1:3) = [real(TRUE_DFX,dp), real(TRUE_DFY,dp), real(TRUE_ANGAST,dp)*DPI/180.d0]
            if( ndim == 4 ) vec(4) = real(phshift,dp)
            step = 1.d-5
            do j = 1, ndim
                vec(j) = vec(j) + step
                call cost%fdf(ndim, vec, f, gp)
                vec(j) = vec(j) - 2.d0*step
                call cost%fdf(ndim, vec, f, gm)
                vec(j) = vec(j) + step
                hfd(:,j) = (gp - gm) / (2.d0*step)
            enddo
            call cost%fdf(ndim, vec, f, g, curv)
            err = maxval(abs(curv - hfd)) / maxval(abs(hfd))
            write(*,'(A,I0,A,ES12.4,A,ES10.3)') '  CTF oracle ndim ', ndim, ': cost ', f, &
                &', max relative difference to the Hessian ', err
            call assert_true(err <= ORACLE_RTOL, 'the CTF curvature is the Hessian at zero residual')
            ridge = 0.d0
            do d = 1, ndim
                ridge(d,d) = 1.d-12*maxval(abs(curv))
            enddo
            call cholesky(curv + ridge, chol, ndim, ok)
            call assert_true(ok .and. all(curv == transpose(curv)), 'the CTF curvature is symmetric positive semi-definite')
            call cost%kill
            call pspec%kill
            deallocate(vec, g, gp, gm, curv, hfd, chol, ridge)
        enddo
    end subroutine test_cont_curvature_oracle

    !> isotropic truth and start (dfx = dfy): the angle has no gradient or curvature; under 'lm' the defocus improves
    !! and the angle does not move (ANGLE_KEEP_TOL)
    subroutine test_cont_isotropic_start()
        class(image), allocatable, target :: pspec
        type(ctf_estimate_cost4Dcont) :: cost
        type(ctfparams), target       :: parms
        integer, allocatable, target :: inds(:,:)
        real :: limits(2,3), cc, err_before, err_after
        write(*,'(A)') 'test_cont_isotropic_start'
        allocate(pspec)
        call build_spectrum(pspec, inds, TRUE_DFX, TRUE_DFX, TRUE_ANGAST, 0., .false.)
        call set_parms(parms, TRUE_DFX + ISO_OFFSET, TRUE_DFX + ISO_OFFSET, TRUE_ANGAST, 0., .false.)
        call fit_limits(3, TRUE_ANGAST, 0., limits)
        call cost%init(pspec, parms, inds, 3, limits, CASTIGTOL, solver='lm')
        err_before = max(abs(parms%dfx - TRUE_DFX), abs(parms%dfy - TRUE_DFX))
        call cost%minimize(parms, cc)
        err_after  = max(abs(parms%dfx - TRUE_DFX), abs(parms%dfy - TRUE_DFX))
        write(*,'(A,2F10.6,A,F12.6,A,I0)') '  CTF isotropic start: dfx dfy ', parms%dfx, parms%dfy, ', angast ', &
            &parms%angast, ', evaluations ', cost%ospec%nevals
        call assert_true(err_after < err_before, 'from an isotropic start the defocus improves')
        call assert_true(abs(parms%angast - TRUE_ANGAST) <= ANGLE_KEEP_TOL, 'from an isotropic start the angle does not move')
        call cost%kill
        call pspec%kill
    end subroutine test_cont_isotropic_start

    !> both solvers from the same eight starts on an astigmatic spectrum normalised as the fit normalises it, each start
    !! timed NTIME times; one LM_COMPARE row per solver (evaluations, iterations, median seconds, defocus and angle error;
    !! reported, not gated)
    subroutine test_cont_solver_comparison()
        use simple_stat, only: median
        integer,          parameter :: NTIME = 5
        character(len=6), parameter :: SOLVERS(2) = ['lbfgsb', 'lm    ']
        real, parameter :: DSTART(2,8) = reshape([0.02,0.01, -0.015,0.02, 0.01,-0.02, -0.02,-0.01, 0.015,0.015, &
            &-0.01,0.015, 0.015,-0.005, -0.02,0.005], [2,8])
        real, parameter :: ASTART(8)   = [4., -5., 2., -3., 5., -4., 1., 3.]
        class(image), allocatable, target :: pspec
        type(ctf_estimate_cost4Dcont) :: cost
        type(ctfparams), target       :: parms
        integer, allocatable, target :: inds(:,:)
        real    :: limits(2,3), cc, err_df, err_ang, sum_df, max_df, sum_ang
        integer :: isolver, istart, nevals, niter, irep, ntimed
        integer(kind=8) :: c0, c1, crate, csum
        real(kind=8)    :: times(NTIME*size(ASTART))
        write(*,'(A)') 'test_cont_solver_comparison'
        allocate(pspec)
        call build_spectrum(pspec, inds, TRUE_DFX, TRUE_DFY, TRUE_ANGAST, 0., .false.)
        do isolver = 1, 2
            nevals = 0
            niter  = 0
            csum   = 0
            sum_df = 0.
            max_df = 0.
            sum_ang = 0.
            ntimed  = 0
            do irep = 1, NTIME
            do istart = 1, size(ASTART)
                call set_parms(parms, TRUE_DFX + DSTART(1,istart), TRUE_DFY + DSTART(2,istart), &
                    &TRUE_ANGAST + ASTART(istart), 0., .false.)
                call fit_limits(3, parms%angast, 0., limits)
                call cost%init(pspec, parms, inds, 3, limits, CASTIGTOL, solver=trim(SOLVERS(isolver)))
                call system_clock(c0, crate)
                call cost%minimize(parms, cc)
                call system_clock(c1)
                csum    = csum + (c1 - c0)
                ntimed  = ntimed + 1
                times(ntimed) = real(c1 - c0, 8) / real(crate, 8)
                if( irep > 1 )then
                    call cost%kill
                    cycle
                endif
                nevals  = nevals + cost%ospec%nevals
                niter   = niter  + cost%ospec%niter
                err_df  = max(abs(parms%dfx - TRUE_DFX), abs(parms%dfy - TRUE_DFY))
                err_ang = abs(modulo(parms%angast - TRUE_ANGAST + 90., 180.) - 90.)
                sum_df  = sum_df + err_df
                max_df  = max(max_df, err_df)
                sum_ang = sum_ang + err_ang
                call cost%kill
            enddo
            enddo
            write(*,'(A,A6,A,F7.2,A,F6.2,A,ES10.3,A,ES10.3,A,F8.5,A,F8.5,A,F8.4)') '  LM_COMPARE ctf_cost4Dcont ', &
                &SOLVERS(isolver), ' evals/solve ', real(nevals)/8., ' iters/solve ', real(niter)/8., &
                &' s/solve median ', median(times(1:ntimed)), ' mean ', real(csum,8)/real(crate,8)/real(ntimed,8), &
                &' mean df err ', sum_df/8., ' max df err ', max_df, ' mean angle err ', sum_ang/8.
        enddo
        call pspec%kill
    end subroutine test_cont_solver_comparison

end submodule simple_ctf_estimate_cost_tester
