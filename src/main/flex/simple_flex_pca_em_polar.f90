!@descr: flex_pca EM: the polar E-step accumulation and its ring/band helpers
submodule (simple_flex_pca_em) simple_flex_pca_em_polar
use simple_matcher_3Drec,   only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch
use simple_flex_reconstructor_latent_ops, only: latent_projection_weights, weighted_expanded_cmat,&
    &LATENT_WDIM
use simple_flex_reconstructor_latent_ops, only: prep_imgs4projected_model
use simple_flex_gpu,        only: flex_gpu_available, flex_gpu_prep_begin_f, flex_gpu_prep_free_f,&
    &flex_gpu_psample_begin_f, flex_gpu_psample_batch_f, flex_gpu_psample_free_f,&
    &flex_gpu_psample_batch_res_f
use simple_flex_pca_polar,  only: polar_grid_build, polar_grid_kill, polar_project_recs,&
    &polar_sample_particle, polar_relative_inplane, polar_assign_directions, polar_sample_at_pose,&
    &polar_apply_shift, polar_dir_neighbours
implicit none
#include "simple_local_flags.inc"

contains

    !> Is the polar (shared-direction) former requested for the reduced solve?
    logical module function cov_polar_enabled()
        cov_polar_enabled = .true.
    end function cov_polar_enabled

    !> Number of bank directions.
    !!
    !! MEASURED (IgG-RL 10k, box_crop 64, lp 15, d_tilde 64): the largest reduced eigenvalue is
    !! 1194.2 / 1195.5 / 1196.2 / 1198.2 at ndir = 1000 / 2000 / 8000 / 32000 and 1198.3 with the
    !! grid removed entirely (SIMPLE_COV_POLAR_EXACT=1), and ground-truth basis capture is 0.5828
    !! at ndir=2000 against 0.5819 exact and 0.5848 Cartesian. A 6 degree direction grid is
    !! indistinguishable from no discretisation at all here, because this stage never uses a
    !! per-particle b on its own -- it accumulates Sbb and sum_i G_i (x) G_i over 10^4-10^5
    !! particles, and the Gram is additionally a sum over ~10^3 plane samples, so direction error
    !! enters suppressed by 1/sqrt(nsamp) rather than as a per-sample decorrelation.
    !!
    !! So the default targets AMORTISATION (~40 particles per direction) rather than resolution,
    !! with a floor so small datasets still get a reasonable grid. Raise it with
    !! SIMPLE_COV_POLAR_NDIR if a dataset ever shows direction sensitivity -- the bank is streamed
    !! direction by direction, so ndir costs no memory, only bank-build time.
    integer module function cov_polar_ndir( nptcls )
        integer, intent(in) :: nptcls
        integer :: v
        v = min(4000, max(1000, nptcls/40))
        v = 2*((v+1)/2)                                  ! build_refspiral needs an even count
        cov_polar_ndir = v
    end function cov_polar_ndir

    !> One ring's contribution to the Gram of the basis (columns 1..ncomp of Us) and to the mean
    !! cross term (column 0). `Cout` is the full ncomp x ncomp block flattened column-major, `Mout`
    !! the ncomp-vector <U_q, T mu>.
    module subroutine polar_ring_gram( Us, ldu, ncomp, row0, nrow, Csp, Cout, Mout )
        integer,  intent(in)    :: ldu, ncomp, row0, nrow
        real,     intent(in)    :: Us(ldu,0:ncomp)
        real,     intent(inout) :: Csp(0:ncomp,0:ncomp)      !< caller-owned scratch
        real(dp), intent(out)   :: Cout(ncomp*ncomp), Mout(ncomp)
        integer :: q, r, i0, n2
        i0 = 2*row0 - 1
        n2 = 2*nrow
        if( n2 <= 0 )then
            Cout = 0.d0; Mout = 0.d0
            return
        endif
        call ssyrk('U','T', ncomp+1, n2, 1.0, Us(i0,0), ldu, 0.0, Csp, ncomp+1)
        do r = 1, ncomp
            do q = 1, r
                Cout((r-1)*ncomp+q) = real(Csp(q,r), dp)
                Cout((q-1)*ncomp+r) = real(Csp(q,r), dp)
            end do
            Mout(r) = real(Csp(0,r), dp)
        end do
    end subroutine polar_ring_gram

    real(dp) module function polar_ring_selfpower( Us, ldu, row0, nrow )
        integer, intent(in) :: ldu, row0, nrow
        real,    intent(in) :: Us(ldu,0:*)
        integer :: j
        polar_ring_selfpower = 0.d0
        do j = 2*row0-1, 2*(row0+nrow-1)
            polar_ring_selfpower = polar_ring_selfpower + real(Us(j,0),dp)*real(Us(j,0),dp)
        end do
    end function polar_ring_selfpower

    !> <y,y> in the polar measure. xws already carries sqrt(wq) and the CTF adjoint, so the CTF has
    !! to be divided back out ring by ring to recover the observation's own energy.
    real(dp) module function polar_self_energy( xws, wr, pg ) result( e )
        real,               intent(in) :: xws(:), wr(:)
        type(polar_grid_t), intent(in) :: pg
        integer  :: ir, j
        real(dp) :: acc
        e = 0.d0
        do ir = 1, pg%nk
            if( real(wr(ir),dp) <= DTINY ) cycle
            acc = 0.d0
            do j = 2*pg%rbeg(ir)-1, 2*pg%rend(ir)
                acc = acc + real(xws(j),dp)*real(xws(j),dp)
            end do
            e = e + acc / real(wr(ir),dp)
        end do
    end function polar_self_energy

    pure real(dp) module function sum_dp_safe( acc, n ) result( v )
        real(dp), intent(in) :: acc
        integer,  intent(in) :: n
        v = acc / real(max(1,n), dp)
    end function sum_dp_safe

    !> thin wrapper so the OpenMP body stays readable; folds sqrt(wq) into the stored samples
    !> resample a stored half-plane by an in-plane rotation into the bank frame: unit-tap 2D KB
    !! at pf-multiples with per-tap Friedel, per-axis normalized weights, OOB taps dropped --
    !! the polar former's interpolation scheme (polar_interp_plane) on the Cartesian lattice.
    !! Positions outside the nyq disk come back ZERO (nothing downstream reads them).
    module subroutine align_halfplane_inplane( frlims, nyq_eff, src, ca, sa, dst )
        integer, intent(in)  :: frlims(3,2), nyq_eff
        complex, intent(in)  :: src(frlims(1,1):frlims(1,2), frlims(2,1):0)
        real,    intent(in)  :: ca, sa
        complex, intent(out) :: dst(frlims(1,1):frlims(1,2), frlims(2,1):0)
        type(kbinterpol) :: kbwin
        real    :: hu, ku, w, wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM)
        integer :: win(2,3), h, k, hlo2, hhi2, klo2, hx, ky, ix, iy, nd, pf
        complex :: acc, cv
        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        pf    = OSMPL_PAD_FAC
        hlo2  = ceil_div (frlims(1,1), pf); hhi2 = floor_div(frlims(1,2), pf)
        klo2  = ceil_div (frlims(2,1), pf)
        nd    = nyq_eff*(nyq_eff+1)
        dst   = CMPLX_ZERO
        do k = klo2, 0
            do h = hlo2, hhi2
                if( h*h + k*k > nd ) cycle
                hu =  h*ca + k*sa
                ku = -h*sa + k*ca
                call latent_projection_weights(kbwin, [hu, ku, 0.], win, wx, wy, wz)
                acc = CMPLX_ZERO
                do iy = 1, LATENT_WDIM
                    ky = win(1,2) + iy - 1
                    do ix = 1, LATENT_WDIM
                        hx = win(1,1) + ix - 1
                        w  = wx(ix)*wy(iy)
                        if( pf*ky <= 0 )then
                            if( pf*hx < frlims(1,1) .or. pf*hx > frlims(1,2) .or. pf*ky < frlims(2,1) ) cycle
                            cv = src(pf*hx, pf*ky)
                        else
                            if( -pf*hx < frlims(1,1) .or. -pf*hx > frlims(1,2) .or. -pf*ky < frlims(2,1) ) cycle
                            cv = conjg(src(-pf*hx, -pf*ky))
                        endif
                        acc = acc + w*cv
                    end do
                end do
                dst(pf*h, pf*k) = acc
            end do
        end do
    end subroutine align_halfplane_inplane

    module subroutine polar_sample_particle_packed( fpl, pg, ca, sa, xws, wr, hfpw, hfcnt, tazim, xws1, xws2 )
        type(fplane_type),  intent(in)    :: fpl
        type(polar_grid_t), intent(in)    :: pg
        real,               intent(in)    :: ca, sa
        real,               intent(out)   :: xws(:)
        real,               intent(out)   :: wr(:)
        real(dp),           intent(inout) :: hfpw, hfcnt
        real,               intent(out)   :: tazim
        !> packed forms of the two lattice-parity half-fields, for the reliability prior
        real, optional,     intent(out)   :: xws1(:), xws2(:)
        complex, allocatable :: xw(:), xw1(:), xw2(:)
        real(dp) :: pw, cnt
        integer  :: j
        allocate(xw(pg%nsamp), xw1(pg%nsamp), xw2(pg%nsamp))
        if( present(xws1) )then
            call polar_sample_particle(fpl%cmplx_plane, fpl%transfer_plane, pg, ca, sa, xw, wr, &
                &pw, cnt, tazim, xw1, xw2)
        else
            call polar_sample_particle(fpl%cmplx_plane, fpl%transfer_plane, pg, ca, sa, xw, wr, &
                &pw, cnt, tazim)
        endif
        do j = 1, pg%nsamp
            xws(2*j-1) = pg%sqwq(j)*real(xw(j))
            xws(2*j)   = pg%sqwq(j)*aimag(xw(j))
        end do
        if( present(xws1) )then
            do j = 1, pg%nsamp
                xws1(2*j-1) = pg%sqwq(j)*real(xw1(j))
                xws1(2*j)   = pg%sqwq(j)*aimag(xw1(j))
                xws2(2*j-1) = pg%sqwq(j)*real(xw2(j))
                xws2(2*j)   = pg%sqwq(j)*aimag(xw2(j))
            end do
        endif
        hfpw  = hfpw  + pw
        hfcnt = hfcnt + cnt
        deallocate(xw, xw1, xw2)
    end subroutine polar_sample_particle_packed

    !> Banded mean projection for the polar E-step: reconstructor%project_fplane's numerics --
    !! the SAME banded (h,k) sweep, apod_mat_3d interpolation weights (including their final
    !! global renormalization), per-sample Friedel conjugation and transfer multiply -- with the
    !! per-call full-plane work removed. project_fplane zero-fills the whole PADDED plane and
    !! copies the reference ctfsq and transfer planes into the output on EVERY call; at the
    !! native padded lattice that is several MB of memory traffic per particle, which measured
    !! as ~80% of the polar E-step's project bucket, all spent on values the polar branch never
    !! reads (only mean_fpl%cmplx_plane is consumed, by the residual subtraction). Here the
    !! plane is zero-filled once at (re)allocation; every call rewrites exactly the in-band disc
    !! samples, and out-of-disc positions stay zero -- the same invariant the Cartesian former's
    !! ensure_latent_projection_plane relies on. The interpolated values are bit-identical to
    !! project_fplane's (same expressions, same kbwin), so the residual planes the M-step
    !! consumes are unchanged wherever the mean is nonzero and unchanged-because-zero elsewhere.
    module subroutine project_fplane_mean_banded( rec, o, fpl_ref, fpl_out )
        type(reconstructor), intent(in)    :: rec
        class(ori),          intent(inout) :: o
        type(fplane_type),   intent(in)    :: fpl_ref
        type(fplane_type),   intent(inout) :: fpl_out
        type(kbinterpol) :: kbwin
        real    :: rotmat(3,3), loc(3), loc_friedel(3), hrow(3)
        real    :: w3(LATENT_WDIM,LATENT_WDIM,LATENT_WDIM)
        integer :: fpllims_pd(3,2), fpllims(3,2), h, k, hp, kp, pf, iwinsz, win(2,3)
        integer :: h_sq, k_max_h, k_lo, k_hi, nyq_disk, nyq_eff
        logical :: l_conjg, l_realloc
        complex :: comp
        kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        iwinsz = ceiling(KBWINSZ - 0.5)
        fpl_out%frlims  = fpl_ref%frlims
        fpl_out%shconst = fpl_ref%shconst
        fpl_out%nyq     = fpl_ref%nyq
        l_realloc = .not. allocated(fpl_out%cmplx_plane)
        if( .not. l_realloc )then
            l_realloc = any(lbound(fpl_out%cmplx_plane) /= lbound(fpl_ref%cmplx_plane)) .or. &
                &any(ubound(fpl_out%cmplx_plane) /= ubound(fpl_ref%cmplx_plane))
        endif
        if( l_realloc )then
            if( allocated(fpl_out%cmplx_plane) ) deallocate(fpl_out%cmplx_plane)
            allocate(fpl_out%cmplx_plane(lbound(fpl_ref%cmplx_plane,1):ubound(fpl_ref%cmplx_plane,1), &
                &lbound(fpl_ref%cmplx_plane,2):ubound(fpl_ref%cmplx_plane,2)))
            fpl_out%cmplx_plane = CMPLX_ZERO
        endif
        rotmat      = o%get_mat()
        pf          = OSMPL_PAD_FAC
        fpllims_pd  = fpl_ref%frlims
        fpllims     = fpllims_pd
        fpllims(1,1)= ceil_div (fpllims_pd(1,1), pf)
        fpllims(1,2)= floor_div(fpllims_pd(1,2), pf)
        fpllims(2,1)= ceil_div (fpllims_pd(2,1), pf)
        fpllims(2,2)= floor_div(fpllims_pd(2,2), pf)
        nyq_eff = rec%get_lfny(1)
        if( fpl_ref%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpl_ref%nyq / pf))
        nyq_disk = nyq_eff * (nyq_eff + 1)
        do h = fpllims(1,1), fpllims(1,2)
            h_sq = h*h
            if( h_sq > nyq_disk ) cycle
            k_max_h = int(sqrt(real(nyq_disk - h_sq)))
            k_lo    = max(fpllims(2,1), -k_max_h)
            k_hi    = min(0, min(fpllims(2,2), k_max_h))
            hp      = h * pf
            hrow(1) = real(h) * rotmat(1,1)
            hrow(2) = real(h) * rotmat(1,2)
            hrow(3) = real(h) * rotmat(1,3)
            do k = k_lo, k_hi
                kp     = k * pf
                loc(1) = hrow(1) + real(k) * rotmat(2,1)
                loc(2) = hrow(2) + real(k) * rotmat(2,2)
                loc(3) = hrow(3) + real(k) * rotmat(2,3)
                ! interp_cmat_exp, verbatim (it is private to simple_reconstructor)
                l_conjg     = loc(1) < 0.
                loc_friedel = loc
                if( l_conjg ) loc_friedel = -loc_friedel
                win(1,:) = nint(loc_friedel)
                win(2,:) = win(1,:) + iwinsz
                win(1,:) = win(1,:) - iwinsz
                call kbwin%apod_mat_3d(loc_friedel, iwinsz, LATENT_WDIM, w3)
                comp = sum(w3 * rec%cmat_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)))
                if( l_conjg ) comp = conjg(comp)
                ! apply_ctf_amp=.true. semantics of project_fplane
                if( allocated(fpl_ref%transfer_plane) )then
                    fpl_out%cmplx_plane(hp,kp) = fpl_ref%transfer_plane(hp,kp) * comp
                else
                    fpl_out%cmplx_plane(hp,kp) = sqrt(max(0., fpl_ref%ctfsq_plane(hp,kp))) * comp
                endif
            end do
        end do
    end subroutine project_fplane_mean_banded

    !> Exact Cartesian statistics of the low-k shells for the HYBRID polar E-step, added on top
    !! of the ring statistics. Per lattice position the KB window geometry is computed once and
    !! all ncomp+1 volumes are gathered through it (the Cartesian former's hoist); the data value,
    !! CTF/whitening transfer and quadrature weight (1 per lattice point) are exactly the
    !! Cartesian former's, so the shells this covers contribute to G/b/c/e_mm/myv precisely what
    !! project_fplanes_mean_basis + cov_herm_inner would contribute for them -- including the DC
    !! sample. This is what removes the ring quadrature's multiplicative posterior-variance bias:
    !! after whitening the low-k shells still anchor the latent scale, and rings sample them
    !! worst (few samples, steep integrand).
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
        type(kbinterpol) :: kbwin
        real        :: rotmat(3,3), loc(3), wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM)
        integer     :: j, q, r, win(2,3), hp, kp, exp_lb(3), exp_ub(3), pf
        logical     :: l_conjg, l_tf
        complex     :: tf, yv, u0, val
        complex     :: uq(ncomp)
        complex(dp) :: u0d, yd
        kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        rotmat = o%get_mat()
        pf     = OSMPL_PAD_FAC
        exp_lb = lbound(rec0%cmat_exp)
        exp_ub = ubound(rec0%cmat_exp)
        l_tf   = allocated(fpl%transfer_plane)
        do j = 1, npos
            loc(1) = real(hex(j))*rotmat(1,1) + real(kex(j))*rotmat(2,1)
            loc(2) = real(hex(j))*rotmat(1,2) + real(kex(j))*rotmat(2,2)
            loc(3) = real(hex(j))*rotmat(1,3) + real(kex(j))*rotmat(2,3)
            l_conjg = loc(1) < 0.
            if( l_conjg ) loc = -loc
            call latent_projection_weights(kbwin, loc, win, wx, wy, wz)
            if( any(win(1,:) < exp_lb) .or. any(win(2,:) > exp_ub) ) cycle
            hp = pf*hex(j)
            kp = pf*kex(j)
            if( l_tf )then
                tf = fpl%transfer_plane(hp,kp)
            else
                tf = cmplx(sqrt(max(0., fpl%ctfsq_plane(hp,kp))), 0.)
            endif
            yv  = fpl%cmplx_plane(hp,kp)
            val = weighted_expanded_cmat(rec0, win, wx, wy, wz)
            if( l_conjg ) val = conjg(val)
            u0 = tf * val
            do q = 1, ncomp
                val = weighted_expanded_cmat(recs(q), win, wx, wy, wz)
                if( l_conjg ) val = conjg(val)
                uq(q) = tf * val
            end do
            u0d  = cmplx(u0, kind=dp)
            yd   = cmplx(yv, kind=dp)
            e_mm = e_mm + real(conjg(u0d)*u0d, dp)
            myv  = myv  + real(conjg(u0d)*yd,  dp)
            do q = 1, ncomp
                bd(q) = bd(q) + real(conjg(cmplx(uq(q),kind=dp))*yd,  dp)
                cd(q) = cd(q) + real(conjg(cmplx(uq(q),kind=dp))*u0d, dp)
                do r = q, ncomp
                    Gd(q,r) = Gd(q,r) + real(conjg(cmplx(uq(q),kind=dp))*cmplx(uq(r),kind=dp), dp)
                end do
            end do
        end do
        ! mirror the accumulated upper triangle (the ring dgemv filled both triangles already;
        ! the exact increments above touched q<=r only)
        do r = 1, ncomp
            do q = r+1, ncomp
                Gd(q,r) = Gd(r,q)
            end do
        end do
    end subroutine polar_hybrid_exact_accum

    !> Banded residual subtraction, fpl = fpl - a*mean over EXACTLY the disc the banded (or any
    !! full-plane) mean projection wrote. Everywhere outside that disc the mean plane is
    !! identically zero, so the full-array statement this replaces only rewrote unchanged values
    !! there -- another few MB of per-particle traffic at the native padded lattice for no effect.
    !! The loop bounds are the same expressions as project_fplane_mean_banded's, so written and
    !! subtracted sample sets coincide by construction.
    module subroutine subtract_mean_banded( fpl, mean_fpl, a, rec_nyq )
        type(fplane_type), intent(inout) :: fpl
        type(fplane_type), intent(in)    :: mean_fpl
        real,              intent(in)    :: a
        integer,           intent(in)    :: rec_nyq
        integer :: fpllims_pd(3,2), fpllims(3,2), h, k, hp, kp, pf
        integer :: h_sq, k_max_h, k_lo, k_hi, nyq_disk, nyq_eff
        pf          = OSMPL_PAD_FAC
        fpllims_pd  = fpl%frlims
        fpllims     = fpllims_pd
        fpllims(1,1)= ceil_div (fpllims_pd(1,1), pf)
        fpllims(1,2)= floor_div(fpllims_pd(1,2), pf)
        fpllims(2,1)= ceil_div (fpllims_pd(2,1), pf)
        fpllims(2,2)= floor_div(fpllims_pd(2,2), pf)
        nyq_eff = rec_nyq
        if( fpl%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpl%nyq / pf))
        nyq_disk = nyq_eff * (nyq_eff + 1)
        do h = fpllims(1,1), fpllims(1,2)
            h_sq = h*h
            if( h_sq > nyq_disk ) cycle
            k_max_h = int(sqrt(real(nyq_disk - h_sq)))
            k_lo    = max(fpllims(2,1), -k_max_h)
            k_hi    = min(0, min(fpllims(2,2), k_max_h))
            hp      = h * pf
            do k = k_lo, k_hi
                kp = k * pf
                fpl%cmplx_plane(hp,kp) = fpl%cmplx_plane(hp,kp) - a*mean_fpl%cmplx_plane(hp,kp)
            end do
        end do
    end subroutine subtract_mean_banded

end submodule simple_flex_pca_em_polar
