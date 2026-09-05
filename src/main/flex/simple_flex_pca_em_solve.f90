!@descr: flex_pca EM: dense SPD/ECM/MCFA solvers and subspace-angle diagnostics
submodule (simple_flex_pca_em) simple_flex_pca_em_solve
implicit none
#include "simple_local_flags.inc"

contains

    !> log(det A) for symmetric positive-definite A, by Cholesky on a private copy.
    !!
    !! This is the term the probe's `resid_energy` has always been missing. resid_energy is the JOINT
    !! MAP objective ||.||^2/sig2 + z'Gamma^-1 z evaluated at zhat, which decreases by construction and
    !! therefore says nothing about convergence. The MARGINAL likelihood needs log det A_i and
    !! log det Gamma as well, and without them the EM has no objective to watch -- which is why the
    !! iteration count was a tuned constant standing in for a stopping rule.
    pure module subroutine spd_logdet_dp( A, n, logdet, ok )
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
    module subroutine deflate_against_basis( imgs, n, basis, nb )
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
    module subroutine cross_half_subspace_angles( eimgs, oimgs, n, svals )
        integer,     intent(in)    :: n
        type(image), intent(inout) :: eimgs(n), oimgs(n)
        real(dp), allocatable, intent(out) :: svals(:)
        real, pointer :: ri(:,:,:), rj(:,:,:)
        real(dp), allocatable :: Gee(:,:), Goo(:,:), Geo(:,:), We(:,:), Wo(:,:)
        real(dp), allocatable :: M(:,:), MtM(:,:), V2(:,:), ev2(:)
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
    pure module function quad_form( M, z, n ) result( val )
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
    module subroutine spd_solve_dp( A, b, n )
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
    module subroutine probe_solve_ecm( n, G, b, c, myv, e_mm, prior_, sig2, nml, a, z_, Ainv_, ldA, lok, quad )
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
    module subroutine probe_solve_mix( n, kmix, G, b, c, myv, e_mm, nml, a, sig2, Ominv, Omxi, lpi, xiOx, &
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
    module subroutine mcfa_init( z, nptcls, ncomp, kmix, gam_sum, nval, xi, ppi, Om )
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
    module subroutine mcfa_condition( ncomp, diag_only, Om, Ominv, ldOm )
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
    module subroutine mcfa_mstep( ncomp, kmix, nval, sr, sm, smm, sai, &
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


    module subroutine spd_inv_dp( A, Ainv, n )
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

end submodule simple_flex_pca_em_solve
