!@descr: polar-Fourier shared-direction basis bank for flex_pca
! Mean + basis sections are projected once per shared direction (cov_polar_ndir) on polar rings; particles are
! sampled at their continuous relative in-plane angle and ring-mean |T|^2 factorises G_qr = sum_k w_i(k) C_qr(k).
! Approximations: direction snap and radial |T|^2 (tazim). Quadrature measure, KB weights and CTF adjoint must
! match the Cartesian cov_herm_inner path.
module simple_flex_pca_polar
use simple_core_module_api, only: cmplx_zero, dp, kbalpha, kbinterpol, kbwinsz, osmpl_pad_fac, pi, &
    &simple_exception, tiny
use simple_reconstructor,                 only: reconstructor
use simple_kbinterpol,                    only: kbinterpol
use simple_math,                          only: ceil_div, floor_div
use simple_flex_reconstructor_latent_ops, only: latent_projection_weights, weighted_expanded_cmat, &
    &LATENT_WDIM
implicit none
private
#include "simple_local_flags.inc"

public :: polar_grid_t, flex_polar_bank, polar_grid_build, polar_grid_kill
public :: polar_project_recs, polar_relative_inplane
public :: polar_assign_directions, polar_sample_particle_fused

integer, parameter :: POLAR_NANG_MIN = 6

!> Polar sampling grid over the half-plane k<=0, ring by ring.
!! `nsamp` samples carry the model at one-lattice-unit angular spacing.
type :: polar_grid_t
    integer :: kfrom = 0, kto = 0, nk = 0, nsamp = 0
    integer :: hlo = 0, hhi = 0, klo = 0                  !< particle-plane bounds, unpadded units
    integer :: ph0 = 0, pk0 = 0                           !< padded array lower bounds of cmplx_plane
    integer,  allocatable :: rbeg(:), rend(:)             !< (nk) sample range of each ring
    real,     allocatable :: rad(:), cs(:), sn(:), wq(:)  !< (nsamp)
    real,     allocatable :: sqwq(:)                      !< (nsamp) sqrt(wq), hoisted out of the inner loop
    real,     allocatable :: wq_ring(:)                   !< (nk) per-sample weight of that ring
    logical :: exists = .false.
  contains
    procedure :: build => polar_grid_build
    procedure :: kill  => polar_grid_kill
end type polar_grid_t

!> One fit's polar E-step service. The fit owns this bank explicitly; no ring geometry,
!! projected basis, direction assignment or thread scratch is shared through module state.
type :: flex_polar_bank
    logical :: l_pol_grid = .false.
    logical :: l_pol_bank_it = .false., l_pol_hyb = .false.
    integer :: ndir_es = 0, nsamp_es = 0, nsamp2_es = 0, nk_es = 0
    integer :: ph0_es = 0, pk0_es = 0, hlo_es = 0, hhi_es = 0, klo_es = 0
    integer :: nyqr_es = 0, nyqb_es = 0, rhyb_es = 0, npos_es = 0
    integer,  allocatable :: hex_es(:), kex_es(:)
    type(polar_grid_t) :: pg_es
    real,     allocatable :: rmatb_es(:,:,:), nrmb_es(:,:)
    real,     allocatable :: cae(:), sae(:)
    integer,  allocatable :: dir_es(:)
    logical,  allocatable :: dused_es(:)
    real,     allocatable :: UsallE(:,:,:)
    real(dp), allocatable :: CfE(:,:,:), Cm0E(:,:,:), c00E(:,:)
    complex,  allocatable :: UbankE(:,:,:)
    real,     allocatable :: CspE(:,:,:)
    real,     allocatable :: xws_es(:,:), wr_es(:,:), Reb_es(:,:)
    real(dp), allocatable :: wrd_es(:,:)
    real :: sec_bank = 0.
  contains
    procedure :: kill => flex_polar_bank_kill
end type flex_polar_bank

contains

    !> Ring-wise polar grid over kfrom..kto. gate_lo makes the measure exclude h^2+k^2 <= gate_lo
    !! (nint shells: above shell r starts at r*(r+1)+1) instead of h^2+k^2 < kfrom^2.
    subroutine polar_grid_build( g, kfrom, kto, hlo, hhi, klo, ph0, pk0, gate_lo )
        class(polar_grid_t), intent(inout) :: g
        integer,            intent(in)    :: kfrom, kto, hlo, hhi, klo, ph0, pk0
        integer, optional,  intent(in)    :: gate_lo
        integer :: r, t, nang, j, ncart, h, k, nyq_disk, hk2
        real    :: phi, dphi, wtot, scal
        logical :: l_gate_lo
        call polar_grid_kill(g)
        if( kfrom < 1 .or. kto < kfrom ) THROW_HARD('invalid band; polar_grid_build')
        l_gate_lo = present(gate_lo)
        g%kfrom = kfrom; g%kto = kto; g%nk = kto - kfrom + 1
        g%hlo = hlo; g%hhi = hhi; g%klo = klo
        g%ph0 = ph0; g%pk0 = pk0
        ! --- band rings
        g%nsamp = 0
        do r = kfrom, kto
            g%nsamp = g%nsamp + polar_nang(r)
        end do
        allocate(g%rbeg(g%nk), g%rend(g%nk), g%wq_ring(g%nk))
        allocate(g%rad(g%nsamp), g%cs(g%nsamp), g%sn(g%nsamp), g%wq(g%nsamp), g%sqwq(g%nsamp))
        j = 0
        do r = kfrom, kto
            nang = polar_nang(r)
            dphi = PI / real(nang)
            g%rbeg(r-kfrom+1) = j + 1
            ! Retain the established even-then-odd ring layout used by every bank consumer.
            do t = 1, nang, 2
                phi = real(t-1) * dphi
                j   = j + 1
                g%rad(j) = real(r); g%cs(j) = cos(phi); g%sn(j) = sin(phi)
                g%wq(j)  = PI * real(r) / real(nang)      ! half-annulus area / #samples
            end do
            do t = 2, nang, 2
                phi = real(t-1) * dphi
                j   = j + 1
                g%rad(j) = real(r); g%cs(j) = cos(phi); g%sn(j) = sin(phi)
                g%wq(j)  = PI * real(r) / real(nang)
            end do
            g%rend(r-kfrom+1) = j
            g%wq_ring(r-kfrom+1) = PI * real(r) / real(nang)
        end do
        ! Renormalise so the polar quadrature carries EXACTLY the total measure of the Cartesian
        ! half-plane lattice inside the same disc. cov_herm_inner counts one unit per lattice point
        ! with the integer disc gate h^2+k^2 <= nyq*(nyq+1); reproducing that total is what keeps
        ! sig2, the eigenvalues and the priors on the same scale between the two paths.
        ncart    = 0
        nyq_disk = kto * (kto + 1)
        do k = -kto, 0
            do h = -kto, merge(0, kto, k == 0)
                hk2 = h*h + k*k
                if( hk2 > nyq_disk ) cycle
                if( l_gate_lo )then
                    if( hk2 <= gate_lo ) cycle
                else
                    if( hk2 < kfrom*kfrom ) cycle
                endif
                ncart = ncart + 1
            end do
        end do
        wtot = sum(g%wq)
        if( wtot > TINY .and. ncart > 0 )then
            scal       = real(ncart) / wtot
            g%wq       = g%wq * scal
            g%wq_ring  = g%wq_ring * scal
        endif
        g%sqwq = sqrt(g%wq)
        g%exists = .true.
    end subroutine polar_grid_build

    pure integer function polar_nang( r )
        integer, intent(in) :: r
        polar_nang = max(POLAR_NANG_MIN, nint(PI * real(r)))
    end function polar_nang

    subroutine polar_grid_kill( g )
        class(polar_grid_t), intent(inout) :: g
        if( allocated(g%rbeg)    ) deallocate(g%rbeg)
        if( allocated(g%rend)    ) deallocate(g%rend)
        if( allocated(g%wq_ring) ) deallocate(g%wq_ring)
        if( allocated(g%rad)     ) deallocate(g%rad)
        if( allocated(g%cs)      ) deallocate(g%cs)
        if( allocated(g%sn)      ) deallocate(g%sn)
        if( allocated(g%wq)      ) deallocate(g%wq)
        if( allocated(g%sqwq)    ) deallocate(g%sqwq)
        g%exists = .false.
    end subroutine polar_grid_kill

    subroutine flex_polar_bank_kill( self )
        class(flex_polar_bank), intent(inout) :: self
        call self%pg_es%kill
        if( allocated(self%UsallE) ) deallocate(self%UsallE)
        if( allocated(self%CfE)    ) deallocate(self%CfE)
        if( allocated(self%Cm0E)   ) deallocate(self%Cm0E)
        if( allocated(self%c00E)   ) deallocate(self%c00E)
        if( allocated(self%UbankE) ) deallocate(self%UbankE)
        if( allocated(self%CspE)   ) deallocate(self%CspE)
        if( allocated(self%xws_es) ) deallocate(self%xws_es)
        if( allocated(self%wr_es)  ) deallocate(self%wr_es)
        if( allocated(self%wrd_es) ) deallocate(self%wrd_es)
        if( allocated(self%Reb_es) ) deallocate(self%Reb_es)
        if( allocated(self%rmatb_es) ) deallocate(self%rmatb_es)
        if( allocated(self%nrmb_es)  ) deallocate(self%nrmb_es)
        if( allocated(self%dir_es) ) deallocate(self%dir_es)
        if( allocated(self%cae)    ) deallocate(self%cae)
        if( allocated(self%sae)    ) deallocate(self%sae)
        if( allocated(self%dused_es) ) deallocate(self%dused_es)
        if( allocated(self%hex_es) ) deallocate(self%hex_es)
        if( allocated(self%kex_es) ) deallocate(self%kex_es)
        self%l_pol_grid = .false.
        self%l_pol_bank_it = .false.
        self%l_pol_hyb = .false.
        self%ndir_es = 0; self%nsamp_es = 0; self%nsamp2_es = 0; self%nk_es = 0
        self%ph0_es = 0; self%pk0_es = 0; self%hlo_es = 0; self%hhi_es = 0; self%klo_es = 0
        self%nyqr_es = 0; self%nyqb_es = 0; self%rhyb_es = 0; self%npos_es = 0
        self%sec_bank = 0.
    end subroutine flex_polar_bank_kill

    !> Polar central sections of `nrec` reconstructors at ONE direction, all at once. The sample
    !! geometry -- 3D location, KB window, weights, in/out-of-lattice test -- depends on the sample
    !! and the orientation alone, so it is built once and the volume loop is hoisted outside it,
    !! exactly as project_fplanes_mean_basis does for the Cartesian sweep.
    subroutine polar_project_recs( rec0, recs, nrec, rotmat, g, out )
        type(reconstructor), intent(in)  :: rec0        !< slot 0 (the mean)
        type(reconstructor), intent(in)  :: recs(nrec)  !< slots 1..nrec (the basis)
        integer,             intent(in)  :: nrec
        real,                intent(in)  :: rotmat(3,3)
        type(polar_grid_t),  intent(in)  :: g
        complex,             intent(out) :: out(:,0:)   !< (nsamp, 0:nrec)
        type(kbinterpol) :: kbwin
        integer, allocatable :: swin(:,:,:), jok(:)
        real,    allocatable :: swx(:,:), swy(:,:), swz(:,:)
        logical, allocatable :: scj(:)
        integer :: ns_ok
        integer :: exp_lb(3), exp_ub(3), j, jj, q, win(2,3)
        real    :: loc(3), hb, kb, wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM)
        complex :: val
        kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        exp_lb = lbound(rec0%cmat_exp)
        exp_ub = ubound(rec0%cmat_exp)
        allocate(swin(2,3,g%nsamp), jok(g%nsamp), scj(g%nsamp))
        allocate(swx(LATENT_WDIM,g%nsamp), swy(LATENT_WDIM,g%nsamp), swz(LATENT_WDIM,g%nsamp))
        ns_ok  = 0
        do j = 1, g%nsamp
            hb     =  g%rad(j) * g%cs(j)
            kb     = -g%rad(j) * g%sn(j)
            loc(1) = hb*rotmat(1,1) + kb*rotmat(2,1)
            loc(2) = hb*rotmat(1,2) + kb*rotmat(2,2)
            loc(3) = hb*rotmat(1,3) + kb*rotmat(2,3)
            scj(j) = loc(1) < 0.
            if( scj(j) ) loc = -loc
            call latent_projection_weights(kbwin, loc, win, wx, wy, wz)
            if( any(win(1,:) < exp_lb) .or. any(win(2,:) > exp_ub) )then
                out(j,:) = CMPLX_ZERO
                cycle
            endif
            ns_ok         = ns_ok + 1
            jok(ns_ok)    = j
            swin(:,:,j)   = win
            swx(:,j)      = wx
            swy(:,j)      = wy
            swz(:,j)      = wz
        end do
        do jj = 1, ns_ok
            j   = jok(jj)
            val = weighted_expanded_cmat(rec0, swin(:,:,j), swx(:,j), swy(:,j), swz(:,j))
            if( scj(j) ) val = conjg(val)
            out(j,0) = val
        end do
        do q = 1, nrec
            do jj = 1, ns_ok
                j   = jok(jj)
                val = weighted_expanded_cmat(recs(q), swin(:,:,j), swx(:,j), swy(:,j), swz(:,j))
                if( scj(j) ) val = conjg(val)
                out(j,q) = val
            end do
        end do
        deallocate(swin, jok, scj, swx, swy, swz)
    end subroutine polar_project_recs

    !> Relative in-plane angle alpha between a particle orientation and a bank direction, defined by
    !!     R_ptcl(1,:) = cos(alpha) R_bank(1,:) + sin(alpha) R_bank(2,:).
    !! A bank sample at plane angle phi then corresponds to particle-plane angle phi + alpha. Taken
    !! from the matrices rather than from e3 so it is independent of the Euler convention and stays
    !! correct when the bank direction is not exactly the particle direction.
    pure subroutine polar_relative_inplane( rot_p, rot_b, ca, sa )
        real, intent(in)  :: rot_p(3,3), rot_b(3,3)
        real, intent(out) :: ca, sa
        real :: c, s, n
        c = rot_p(1,1)*rot_b(1,1) + rot_p(1,2)*rot_b(1,2) + rot_p(1,3)*rot_b(1,3)
        s = rot_p(1,1)*rot_b(2,1) + rot_p(1,2)*rot_b(2,2) + rot_p(1,3)*rot_b(2,3)
        n = sqrt(c*c + s*s)
        if( n > TINY )then
            ca = c / n; sa = s / n
        else
            ca = 1.0;   sa = 0.0
        endif
    end subroutine polar_relative_inplane

    !> Allocation-free polar sampler for the shared-direction E-step: one KB window per ring sample
    !! serves the data and transfer gathers, and xws is written sqrt(wq)-packed.
    subroutine polar_sample_particle_fused( cplane, tplane, g, ca, sa, xws, wr, tazim )
        complex,            intent(in)  :: cplane(:,:)
        complex,            intent(in)  :: tplane(:,:)
        type(polar_grid_t), intent(in)  :: g
        real,               intent(in)  :: ca, sa
        real,               intent(out) :: xws(:)
        real,               intent(out) :: wr(:)
        real,               intent(out) :: tazim
        type(kbinterpol) :: kbwin
        integer :: j, ir, nang, i, iwinsz, wlox, wloy, ix, iy, hx, ky, pf
        real    :: hu, ku, c1, s1, t2, tm, tv, bx, by, sx, sy, w, wyy, inv_wdim, eps_norm
        real    :: wx(LATENT_WDIM), wy(LATENT_WDIM)
        complex :: yv, tv_c, xwj
        kbwin    = kbinterpol(KBWINSZ, KBALPHA)
        pf       = OSMPL_PAD_FAC
        iwinsz   = ceiling(KBWINSZ - 0.5)
        inv_wdim = 1.0 / real(LATENT_WDIM)
        eps_norm = epsilon(1.0)
        tazim    = 0.
        do ir = 1, g%nk
            tm = 0.; tv = 0.
            do j = g%rbeg(ir), g%rend(ir)
                c1 = g%cs(j)*ca - g%sn(j)*sa            ! cos(phi + alpha)
                s1 = g%sn(j)*ca + g%cs(j)*sa            ! sin(phi + alpha)
                hu =  g%rad(j) * c1
                ku = -g%rad(j) * s1
                ! window geometry once per sample (latent_projection_weights' x/y axes verbatim)
                wlox = nint(hu) - iwinsz
                wloy = nint(ku) - iwinsz
                bx   = real(wlox) - hu
                by   = real(wloy) - ku
                do i = 1, LATENT_WDIM
                    wx(i) = kbwin%apod(bx + real(i-1))
                    wy(i) = kbwin%apod(by + real(i-1))
                end do
                sx = sum(wx)
                sy = sum(wy)
                if( abs(sx) > eps_norm )then
                    wx = wx * (1.0 / sx)
                else
                    wx = inv_wdim
                endif
                if( abs(sy) > eps_norm )then
                    wy = wy * (1.0 / sy)
                else
                    wy = inv_wdim
                endif
                ! fused gather: both planes through the one tap set, per-tap Friedel as before
                yv   = CMPLX_ZERO
                tv_c = CMPLX_ZERO
                do iy = 1, LATENT_WDIM
                    ky  = wloy + iy - 1
                    wyy = wy(iy)
                    do ix = 1, LATENT_WDIM
                        hx = wlox + ix - 1
                        w  = wx(ix) * wyy
                        if( ky > 0 )then
                            if( -hx < g%hlo .or. -hx > g%hhi .or. -ky < g%klo ) cycle
                            yv   = yv   + w * conjg(cplane(pf*(-hx) - g%ph0 + 1, pf*(-ky) - g%pk0 + 1))
                            tv_c = tv_c + w * conjg(tplane(pf*(-hx) - g%ph0 + 1, pf*(-ky) - g%pk0 + 1))
                        else
                            if( hx < g%hlo .or. hx > g%hhi .or. ky < g%klo ) cycle
                            yv   = yv   + w * cplane(pf*hx - g%ph0 + 1, pf*ky - g%pk0 + 1)
                            tv_c = tv_c + w * tplane(pf*hx - g%ph0 + 1, pf*ky - g%pk0 + 1)
                        endif
                    end do
                end do
                xwj        = conjg(tv_c) * yv
                xws(2*j-1) = g%sqwq(j)*real (xwj)
                xws(2*j)   = g%sqwq(j)*aimag(xwj)
                t2 = real(tv_c*conjg(tv_c))
                tm = tm + t2
                tv = tv + t2*t2
            end do
            nang   = g%rend(ir) - g%rbeg(ir) + 1
            tm     = tm / real(nang)
            wr(ir) = tm
            tv     = max(0., tv/real(nang) - tm*tm)
            if( tm > TINY ) tazim = tazim + sqrt(tv)/tm
        end do
        tazim = tazim / real(max(1,g%nk))
    end subroutine polar_sample_particle_fused

    !> Nearest bank direction for every particle, by the plane normal (row 3 of the rotation
    !! matrix). Done as a BLAS-3 sweep because a per-particle scan over the direction grid is
    !! O(N*nspace) scalar work and nspace is deliberately large here.
    subroutine polar_assign_directions( nrm_p, nptcls, nrm_b, ndir, dir_of )
        real,    intent(in)  :: nrm_p(3,nptcls)
        integer, intent(in)  :: nptcls, ndir
        real,    intent(in)  :: nrm_b(3,ndir)
        integer, intent(out) :: dir_of(nptcls)
        integer, parameter :: CHUNK = 2048
        real,    allocatable :: dots(:,:)
        integer :: i0, i1, nc, i, j, jbest
        real    :: best
        allocate(dots(ndir,CHUNK))
        do i0 = 1, nptcls, CHUNK
            i1 = min(nptcls, i0 + CHUNK - 1)
            nc = i1 - i0 + 1
            call sgemm('T','N', ndir, nc, 3, 1.0, nrm_b, 3, nrm_p(1,i0), 3, 0.0, dots, ndir)
            !$omp parallel do default(shared) private(i,j,jbest,best) schedule(static) proc_bind(close)
            do i = 1, nc
                jbest = 1
                best  = dots(1,i)
                do j = 2, ndir
                    if( dots(j,i) > best )then
                        best  = dots(j,i)
                        jbest = j
                    endif
                end do
                dir_of(i0+i-1) = jbest
            end do
            !$omp end parallel do
        end do
        deallocate(dots)
    end subroutine polar_assign_directions



end module simple_flex_pca_polar
