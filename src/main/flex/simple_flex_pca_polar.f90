!@descr: polar-Fourier shared-direction basis bank for flex_pca
!
! WHY THIS EXISTS
! ---------------
! Every per-particle stage of flex_pca (reduced solve, embedding, probe E-step) extracts
! `d+1` Cartesian central sections at the particle's own continuous orientation and then forms
! `G_qr = <B_q,B_r>` and `b_q = <B_q,y>` over that plane. Because every particle sits at its own
! orientation, nothing is shared: the projection cost is `N * (d+1) * plane_samples * kb_taps`
! and the Gram cost is `N * d^2 * plane_samples`. Measured on Ribosembly (335k particles,
! d_tilde=177) those two terms are 3786 and 3763 thread-seconds against 87 s for everything else
! in the solve -- 98.9 % of the stage.
!
! Both terms collapse if the model is evaluated on a DISCRETE grid of projection DIRECTIONS
! instead of per particle:
!
!   1. the `d+1` sections at one direction serve every particle assigned to it, so the projection
!      cost becomes `nspace * (d+1) * ...` and stops scaling with N;
!   2. the noise/CTF weight `w_i` is RADIAL, so the theta-sum in `G` factorises:
!         G_qr(i) = sum_k w_i(k) C_qr(k),   C_qr(k) = sum_{theta in ring k} U_q conj(U_r)
!      and `C` is shared by every particle at that direction. The per-particle Gram drops from
!      `d^2 * nsamp` to `d^2 * nk` -- a factor of ~40 at box_crop 64.
!
! The in-plane angle stays CONTINUOUS: the bank is built at `e3 = 0` and the particle is polar-
! sampled at the relative in-plane angle recovered from the two rotation matrices, so no in-plane
! discretisation is introduced. The ONLY approximation is snapping the out-of-plane direction to
! the grid, and `nspace` is a free knob that controls it (the bank is streamed direction by
! direction, so nspace costs almost nothing in memory).
!
! CONVENTIONS -- these must match simple_flex_pca_em exactly or the two paths are not
! comparable:
!   * the Cartesian inner product `cov_herm_inner` sums over the HALF-plane k<=0 of the UNPADDED
!     lattice (the plane is stored on a stride-OSMPL_PAD_FAC grid) with unit weight per lattice
!     point. The polar grid therefore carries a quadrature weight `wq` normalised so that
!     `sum(wq)` equals the number of Cartesian lattice points in the same disc.
!   * the basis is interpolated out of `cmat_exp` with `latent_projection_weights` /
!     `weighted_expanded_cmat`, i.e. the SAME normalised KB kernel `project_fplanes_mean_basis`
!     uses -- NOT `interp_cmat_exp`, whose weights are unnormalised.
!   * the CTF is NOT baked into the bank (it cannot be, the bank is shared). The Cartesian path
!     computes `G = <T U_q, T U_r>` and `b = <T U_q, y>`; here the same numbers come from
!     `|T|^2` folded into the radial weight and `conj(T) y` folded into the particle.
!
module simple_flex_pca_polar
use simple_core_module_api
use simple_reconstructor, only: reconstructor
use simple_kbinterpol,    only: kbinterpol
use simple_math,          only: ceil_div, floor_div
use simple_flex_reconstructor_latent_ops, only: latent_projection_weights, weighted_expanded_cmat, &
    &LATENT_WDIM
implicit none
private
#include "simple_local_flags.inc"

public :: polar_grid_t, polar_grid_build, polar_grid_kill
public :: polar_project_recs, polar_sample_particle, polar_relative_inplane
public :: polar_assign_directions, polar_sample_at_pose, polar_apply_shift
public :: polar_dir_neighbours, polar_sample_particle_fused

!> Angular sampling density relative to the Cartesian lattice. 1.0 puts `nint(PI*r)` angles on
!! ring `r`, i.e. an arc spacing of one lattice unit -- the same sample density as the Cartesian
!! half-plane. Raising it oversamples the ring and correlates neighbouring samples without adding
!! information; the noise identity E[Re(b)Re(b)] = 0.5*sig2*G assumes uncorrelated samples, so
!! 1.0 is the value that keeps the polar path on the same footing as the Cartesian one.
real,    parameter, public :: POLAR_ANG_OVERSAMP = 1.0
integer, parameter         :: POLAR_NANG_MIN     = 6

!> Polar sampling grid over the half-plane k<=0, ring by ring.
!! `nsamp` band samples carry the model; `nsamp_n` noise-ring samples exist only so the whitened
!! noise variance can be measured in the SAME representation the Gram is measured in (see
!! polar_sample_particle).
type :: polar_grid_t
    integer :: kfrom = 0, kto = 0, nk = 0, nsamp = 0
    integer :: knfrom = 0, knto = 0, nsamp_n = 0
    integer :: hlo = 0, hhi = 0, klo = 0                  !< particle-plane bounds, unpadded units
    integer :: ph0 = 0, pk0 = 0                           !< padded array lower bounds of cmplx_plane
    integer,  allocatable :: rbeg(:), rend(:)             !< (nk) sample range of each ring
    !> First sample of the SECOND half within each ring. Angles are stored even-indices-first, so
    !! [rbeg,rmid-1] and [rmid,rend] are two interleaved half-sets of the ring -- the polar analogue
    !! of the Cartesian checkerboard split cov_herm_inner takes with `half`. Both are contiguous, so
    !! the half-set Grams and b-vectors are still one BLAS call each, and the ring block as a whole
    !! stays contiguous so nothing in the solve changes shape.
    integer,  allocatable :: rmid(:)                      !< (nk)
    !> (2*nsamp) half-set label of each INTERLEAVED row of the packed real representation, so a
    !! split-half score can be accumulated in one pass without touching the ring bookkeeping.
    integer,  allocatable :: hrow(:)
    real,     allocatable :: rad(:), cs(:), sn(:), wq(:)  !< (nsamp)
    real,     allocatable :: sqwq(:)                      !< (nsamp) sqrt(wq), hoisted out of the inner loop
    real,     allocatable :: wq_ring(:)                   !< (nk) per-sample weight of that ring
    real,     allocatable :: nrad(:), ncs(:), nsn(:), nwq(:)  !< (nsamp_n)
    logical :: exists = .false.
end type polar_grid_t

contains

    !> Build the ring-wise polar grid. `kto` is the band edge (the shell beyond which every model
    !! plane is identically zero) and [knfrom,knto] are the pure-noise rings used for sig2.
    !! `ang_osamp` multiplies the per-ring angular sample count (accuracy knob; default 1 is the
    !! historical grid). `gate_lo` replaces the inner measure gate: when present, the Cartesian
    !! point count the quadrature is normalized to excludes h^2+k^2 <= gate_lo instead of
    !! h^2+k^2 < kfrom^2 -- the SHELL-complement convention the hybrid exact/ring split needs
    !! (shells are nint(sqrt())-assigned, so the annulus above shell r starts at r*(r+1)+1, not
    !! at (r+1)^2; using the default gate there would double-count the boundary points' measure).
    subroutine polar_grid_build( g, kfrom, kto, knfrom, knto, hlo, hhi, klo, ph0, pk0, &
            &ang_osamp, gate_lo )
        type(polar_grid_t), intent(inout) :: g
        integer,            intent(in)    :: kfrom, kto, knfrom, knto, hlo, hhi, klo, ph0, pk0
        integer, optional,  intent(in)    :: ang_osamp, gate_lo
        integer :: r, t, nang, j, ncart, h, k, nyq_disk, osamp, hk2
        real    :: phi, dphi, wtot, scal
        logical :: l_gate_lo
        call polar_grid_kill(g)
        if( kfrom < 1 .or. kto < kfrom ) THROW_HARD('invalid band; polar_grid_build')
        osamp = 1
        if( present(ang_osamp) ) osamp = max(1, ang_osamp)
        l_gate_lo = present(gate_lo)
        g%kfrom = kfrom; g%kto = kto; g%nk = kto - kfrom + 1
        g%knfrom = knfrom; g%knto = knto
        g%hlo = hlo; g%hhi = hhi; g%klo = klo
        g%ph0 = ph0; g%pk0 = pk0
        ! --- band rings
        g%nsamp = 0
        do r = kfrom, kto
            g%nsamp = g%nsamp + osamp*polar_nang(r)
        end do
        allocate(g%rbeg(g%nk), g%rend(g%nk), g%rmid(g%nk), g%wq_ring(g%nk))
        allocate(g%rad(g%nsamp), g%cs(g%nsamp), g%sn(g%nsamp), g%wq(g%nsamp), g%sqwq(g%nsamp))
        j = 0
        do r = kfrom, kto
            nang = osamp*polar_nang(r)
            dphi = PI / real(nang)
            g%rbeg(r-kfrom+1) = j + 1
            ! even angle indices first, then odd: an interleaved half-set split that is contiguous in
            ! storage. Splitting by a contiguous arc instead would correlate the halves with
            ! orientation, which is exactly what a reliability estimate must not do.
            do t = 1, nang, 2
                phi = real(t-1) * dphi
                j   = j + 1
                g%rad(j) = real(r); g%cs(j) = cos(phi); g%sn(j) = sin(phi)
                g%wq(j)  = PI * real(r) / real(nang)      ! half-annulus area / #samples
            end do
            g%rmid(r-kfrom+1) = j + 1
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
        allocate(g%hrow(2*g%nsamp))
        do r = 1, g%nk
            do j = g%rbeg(r), g%rend(r)
                t = merge(1, 2, j < g%rmid(r))
                g%hrow(2*j-1) = t
                g%hrow(2*j)   = t
            end do
        end do
        ! --- noise rings (particle only; every model plane is zero out here)
        g%nsamp_n = 0
        if( knto >= knfrom .and. knfrom >= 1 )then
            do r = knfrom, knto
                g%nsamp_n = g%nsamp_n + polar_nang(r)
            end do
            allocate(g%nrad(g%nsamp_n), g%ncs(g%nsamp_n), g%nsn(g%nsamp_n), g%nwq(g%nsamp_n))
            j = 0
            do r = knfrom, knto
                nang = polar_nang(r)
                dphi = PI / real(nang)
                do t = 1, nang
                    phi = real(t-1) * dphi
                    j   = j + 1
                    g%nrad(j) = real(r)
                    g%ncs(j)  = cos(phi)
                    g%nsn(j)  = sin(phi)
                    g%nwq(j)  = PI * real(r) / real(nang)
                end do
            end do
        else
            allocate(g%nrad(0), g%ncs(0), g%nsn(0), g%nwq(0))
        endif
        g%exists = .true.
    end subroutine polar_grid_build

    pure integer function polar_nang( r )
        integer, intent(in) :: r
        polar_nang = max(POLAR_NANG_MIN, nint(PI * real(r) * POLAR_ANG_OVERSAMP))
    end function polar_nang

    subroutine polar_grid_kill( g )
        type(polar_grid_t), intent(inout) :: g
        if( allocated(g%rbeg)    ) deallocate(g%rbeg)
        if( allocated(g%rend)    ) deallocate(g%rend)
        if( allocated(g%rmid)    ) deallocate(g%rmid)
        if( allocated(g%hrow)    ) deallocate(g%hrow)
        if( allocated(g%wq_ring) ) deallocate(g%wq_ring)
        if( allocated(g%rad)     ) deallocate(g%rad)
        if( allocated(g%cs)      ) deallocate(g%cs)
        if( allocated(g%sn)      ) deallocate(g%sn)
        if( allocated(g%wq)      ) deallocate(g%wq)
        if( allocated(g%sqwq)    ) deallocate(g%sqwq)
        if( allocated(g%nrad)    ) deallocate(g%nrad)
        if( allocated(g%ncs)     ) deallocate(g%ncs)
        if( allocated(g%nsn)     ) deallocate(g%nsn)
        if( allocated(g%nwq)     ) deallocate(g%nwq)
        g%exists = .false.
    end subroutine polar_grid_kill

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

    !> Polar-sample one particle plane at the in-plane angle (ca,sa) relative to its bank direction.
    !!
    !! Returns, on the BAND rings,
    !!    xw(j) = conj(T_i(j)) * y_i(j)      (the whitened, CTF-adjoint observation)
    !!    wr(r) = ring mean of |T_i|^2       (the radial weight the Gram factorisation needs)
    !! and on the NOISE rings the weighted power of the observation, which above the band equals the
    !! residual because every model plane is zero there. sig2 measured this way is in the SAME
    !! representation as G and b -- including whatever variance the polar resampling removes -- so
    !! the debias identity E[Re(b_q)Re(b_r)] = 0.5*sig2*G stays self-consistent.
    !!
    !! `tazim` reports the mean relative azimuthal spread of |T|^2 within a ring. It is 0 for a
    !! non-astigmatic CTF and is the quantity that invalidates the radial factorisation if it is not.
    subroutine polar_sample_particle( cplane, tplane, g, ca, sa, xw, wr, hfpw, hfcnt, tazim, xw1, xw2 )
        complex,            intent(in)  :: cplane(:,:)
        complex,            intent(in)  :: tplane(:,:)
        type(polar_grid_t), intent(in)  :: g
        real,               intent(in)  :: ca, sa
        complex,            intent(out) :: xw(:)
        real,               intent(out) :: wr(:)
        real(dp),           intent(out) :: hfpw, hfcnt
        real,               intent(out) :: tazim
        !> the two INDEPENDENT half-fields, by Cartesian lattice parity (see below)
        complex, optional,  intent(out) :: xw1(:), xw2(:)
        type(kbinterpol) :: kbwin
        integer :: j, r, ir, nang
        real    :: hu, ku, c1, s1, t2, tm, tv
        complex :: yv, tv_c
        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        tazim = 0.
        do ir = 1, g%nk
            tm = 0.; tv = 0.
            do j = g%rbeg(ir), g%rend(ir)
                c1 = g%cs(j)*ca - g%sn(j)*sa            ! cos(phi + alpha)
                s1 = g%sn(j)*ca + g%cs(j)*sa            ! sin(phi + alpha)
                hu =  g%rad(j) * c1
                ku = -g%rad(j) * s1
                yv   = polar_interp_plane(cplane, g, kbwin, hu, ku)
                tv_c = polar_interp_plane(tplane, g, kbwin, hu, ku)
                xw(j) = conjg(tv_c) * yv
                ! THE HALF-SPLIT. Alternating polar SAMPLES does not give independent halves:
                ! consecutive samples on a ring are ~1 lattice unit apart by construction, while the
                ! 3-tap KB kernel spans +-1, so neighbours interpolate from almost the same Cartesian
                ! points and their noise is shared. Measured on EMPIAR-10076, that drove every
                ! split-half rho to 0.95-1.00 and flattened the reliability prior to nothing.
                ! Splitting the underlying LATTICE by parity -- the same convention cov_herm_inner
                ! uses -- gives halves that are functions of disjoint samples.
                if( present(xw1) ) xw1(j) = conjg(tv_c) * polar_interp_plane(cplane, g, kbwin, hu, ku, 1)
                if( present(xw2) ) xw2(j) = conjg(tv_c) * polar_interp_plane(cplane, g, kbwin, hu, ku, 2)
                t2    = real(tv_c*conjg(tv_c))
                tm    = tm + t2
                tv    = tv + t2*t2
            end do
            nang   = g%rend(ir) - g%rbeg(ir) + 1
            tm     = tm / real(nang)
            wr(ir) = tm
            tv     = max(0., tv/real(nang) - tm*tm)
            if( tm > TINY ) tazim = tazim + sqrt(tv)/tm
        end do
        tazim = tazim / real(max(1,g%nk))
        hfpw = 0.d0; hfcnt = 0.d0
        do j = 1, g%nsamp_n
            c1 = g%ncs(j)*ca - g%nsn(j)*sa
            s1 = g%nsn(j)*ca + g%ncs(j)*sa
            hu =  g%nrad(j) * c1
            ku = -g%nrad(j) * s1
            yv = polar_interp_plane(cplane, g, kbwin, hu, ku)
            hfpw  = hfpw  + real(g%nwq(j),dp) * real(yv*conjg(yv), dp)
            hfcnt = hfcnt + real(g%nwq(j),dp)
        end do
    end subroutine polar_sample_particle

    !> Zero-allocation fused variant of polar_sample_particle for the shared-direction E-step:
    !! the KB window geometry of each ring sample is computed ONCE and shared by the data-plane
    !! and transfer-plane gathers (the general sampler interpolates each plane separately, which
    !! recomputes the identical geometry twice per sample), the unused z-axis weights are never
    !! evaluated (polar_interp_plane computes them and drops them), and the sqrt(wq)-packed real
    !! representation is written directly so no complex temporary -- and no per-call allocation --
    !! exists at all. Numerics: each accumulator (yv, tv_c, tm, tv, tazim) sees exactly the taps,
    !! weights and operation order of the original path, so xws/wr/tazim are bit-identical to
    !! polar_sample_particle_packed without halves. No noise rings (the E-step grid has none) and
    !! no parity half-fields (its callers never request them) -- the general sampler stays the
    !! entry point for the embed path, which needs both.
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

    !> Polar-sample a particle at a TRIAL pose: in-plane angle (ca,sa) plus an extra Fourier-space
    !! shift phase. The plane already carries the project's own shift, so (px,py) is the increment.
    !!
    !! This is the operation that makes pose refinement cheap in polar coordinates. The in-plane
    !! angle is just a different set of sampling angles -- no interpolation of the model, no
    !! re-projection of any volume -- and the shift is a phase multiply that needs no resampling at
    !! all, so a whole shift grid can be scored from ONE angular resampling.
    subroutine polar_sample_at_pose( cplane, tplane, g, ca, sa, px, py, xw )
        complex,            intent(in)  :: cplane(:,:), tplane(:,:)
        type(polar_grid_t), intent(in)  :: g
        real,               intent(in)  :: ca, sa, px, py
        complex,            intent(out) :: xw(:)
        type(kbinterpol) :: kbwin
        integer :: j
        real    :: c1, s1, hu, ku, ph
        complex :: yv, tv_c
        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        do j = 1, g%nsamp
            c1 = g%cs(j)*ca - g%sn(j)*sa
            s1 = g%sn(j)*ca + g%cs(j)*sa
            hu =  g%rad(j) * c1
            ku = -g%rad(j) * s1
            yv    = polar_interp_plane(cplane, g, kbwin, hu, ku)
            tv_c  = polar_interp_plane(tplane, g, kbwin, hu, ku)
            ph    = hu*px + ku*py
            xw(j) = conjg(tv_c) * yv * cmplx(cos(ph), sin(ph))
        end do
    end subroutine polar_sample_at_pose

    !> Apply only the shift phase increment to samples already taken at some angle. No resampling,
    !! so an entire shift grid costs one complex multiply per sample per trial.
    pure subroutine polar_apply_shift( g, ca, sa, px, py, xw_in, xw_out )
        type(polar_grid_t), intent(in)  :: g
        real,               intent(in)  :: ca, sa, px, py
        complex,            intent(in)  :: xw_in(:)
        complex,            intent(out) :: xw_out(:)
        integer :: j
        real    :: c1, s1, hu, ku, ph
        do j = 1, g%nsamp
            c1 = g%cs(j)*ca - g%sn(j)*sa
            s1 = g%sn(j)*ca + g%cs(j)*sa
            hu =  g%rad(j) * c1
            ku = -g%rad(j) * s1
            ph = hu*px + ku*py
            xw_out(j) = xw_in(j) * cmplx(cos(ph), sin(ph))
        end do
    end subroutine polar_apply_shift

    !> KB interpolation of a stored Fourier half-plane at continuous UNPADDED lattice coordinates.
    !! gen_fplane4rec fills only indices that are multiples of OSMPL_PAD_FAC, so the data lattice is
    !! the unpadded one; the Friedel mate is taken per tap, because a window centred near k=0
    !! straddles the stored boundary.
    !> `parity`: 0 uses every tap; 1 and 2 use only the taps whose UNPADDED lattice index has
    !! (hx+ky) even / odd, renormalised. This is how an INDEPENDENT half-split is obtained in the
    !! polar representation -- see polar_sample_particle for why alternating polar samples is not.
    pure complex function polar_interp_plane( plane, g, kbwin, hu, ku, parity ) result( val )
        complex,            intent(in) :: plane(:,:)
        type(polar_grid_t), intent(in) :: g
        type(kbinterpol),   intent(in) :: kbwin
        real,               intent(in) :: hu, ku
        integer, optional,  intent(in) :: parity
        real    :: loc(3), wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM), w
        integer :: win(2,3), ix, iy, hx, ky, pf, ipar
        real    :: wsum
        complex :: cv
        pf     = OSMPL_PAD_FAC
        ipar   = 0
        if( present(parity) ) ipar = parity
        loc    = [hu, ku, 0.]
        call latent_projection_weights(kbwin, loc, win, wx, wy, wz)
        val  = CMPLX_ZERO
        wsum = 0.
        do iy = 1, LATENT_WDIM
            ky = win(1,2) + iy - 1
            do ix = 1, LATENT_WDIM
                hx = win(1,1) + ix - 1
                w  = wx(ix) * wy(iy)
                if( ipar /= 0 )then
                    if( modulo(hx + ky, 2) + 1 /= ipar ) cycle
                endif
                if( ky > 0 )then
                    ! Friedel mate, per tap: a window centred near k=0 straddles the stored half
                    if( -hx < g%hlo .or. -hx > g%hhi .or. -ky < g%klo ) cycle
                    cv = conjg(plane(pf*(-hx) - g%ph0 + 1, pf*(-ky) - g%pk0 + 1))
                else
                    if( hx < g%hlo .or. hx > g%hhi .or. ky < g%klo ) cycle
                    cv = plane(pf*hx - g%ph0 + 1, pf*ky - g%pk0 + 1)
                endif
                val  = val + w * cv
                wsum = wsum + w
            end do
        end do
        ! renormalise so each parity sub-kernel is still a partition of unity: the two halves then
        ! carry the same signal and independent noise, which is what a reliability split requires
        if( ipar /= 0 .and. abs(wsum) > 1.e-12 ) val = val / wsum
    end function polar_interp_plane

    !> The `k` nearest bank directions to each bank direction, by plane-normal dot product. This is
    !! the candidate set for the P2 direction search, and restricting to it IS the cone limit: this
    !! is refinement, not global search -- a particle that wants to move across the sphere is telling
    !! you about the consensus, not about conformation.
    subroutine polar_dir_neighbours( ndir, k, nrm_b, nnmat )
        integer, intent(in)  :: ndir, k
        real,    intent(in)  :: nrm_b(3,ndir)
        integer, intent(out) :: nnmat(k,ndir)
        real,    allocatable :: d(:,:)
        integer :: i, j, m, jbest, ithr
        real    :: best
        allocate(d(ndir, omp_get_max_threads()))
        !$omp parallel do default(shared) private(i,j,m,jbest,best,ithr) schedule(static) proc_bind(close)
        do i = 1, ndir
            ithr = omp_get_thread_num() + 1
            do j = 1, ndir
                d(j,ithr) = nrm_b(1,i)*nrm_b(1,j) + nrm_b(2,i)*nrm_b(2,j) + nrm_b(3,i)*nrm_b(3,j)
            end do
            do m = 1, k
                best = -2.0; jbest = 1
                do j = 1, ndir
                    if( d(j,ithr) > best )then
                        best = d(j,ithr); jbest = j
                    endif
                end do
                nnmat(m,i)   = jbest
                d(jbest,ithr) = -3.0
            end do
        end do
        !$omp end parallel do
        deallocate(d)
    end subroutine polar_dir_neighbours


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
