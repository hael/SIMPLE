!@descr: the lattice family of the nano class-average bank (Fourier spot measurements of one class average) and the model gate
! Per class average: mean removed; band presence on the unmasked image; then the image soft-masked (mask2D_soft geometry
! at box/2 - 2), its full-plane power P = |F|^2, smoothed by a periodic 3 x 3 box (Ps), ring means of Ps on every shell
! (mu), the smoothed ratio R = Ps / mu over the band (3 A to Nyquist minus the graphene shells). lattice_peak, the log10
! of the largest R, is a normalized measurement and enters the bank raw; LAT_PEAK_MIN is the lattice grade
! (doc/policies/sieving_and_rejection/model_cavgs_rejection.md, nano context). The model gate, trajectory_solve3D's
! second rejection, cuts gate 1's survivors on their correlation to the simulated start model (K_MODEL).
module simple_cavg_quality_lattice
!$ use omp_lib
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs,               only: dp, GRAPHENE_BAND1, GRAPHENE_BAND2, GRAPHENE_BAND3, COSMSKHALFWIDTH, DPI
use simple_error,              only: simple_exception
use simple_image,              only: image
use simple_math_ft,            only: calc_graphene_mask
use simple_stat,               only: median, mad_gau
implicit none
private
#include "simple_local_flags.inc"

public :: cavg_quality_lattice
public :: LAT_NFEATS, LAT_LATTICE_PEAK, LAT_SPOT_COUNT, LAT_SPOT_ARC, LAT_SPOT_DIAM, LAT_SPOT_SNR, LAT_BAND_PRESENCE
public :: LAT_PEAK_MIN, LAT_PEAK_MIN_RAW, K_MODEL
public :: model_gate_states

! the six lattice features, in the order of the nano bank
integer, parameter :: LAT_NFEATS = 6
integer, parameter :: LAT_LATTICE_PEAK = 1, LAT_SPOT_COUNT = 2, LAT_SPOT_ARC = 3, LAT_SPOT_DIAM = 4, LAT_SPOT_SNR = 5
integer, parameter :: LAT_BAND_PRESENCE = 6
integer,  parameter :: NLOC          = 64        !< local maxima kept per class, strongest first (spot_count)
integer,  parameter :: SPOT_HALFWIN  = 7         !< the 15 x 15 window about the strongest spot
integer,  parameter :: NWIN          = (2 * SPOT_HALFWIN + 1)**2
real,     parameter :: BAND_LOW_RES  = 3.        !< the band runs from this resolution (A) to Nyquist
!> the lattice grade: a smoothed peak at least this far (log10) above its ring mean, 4.5 times the ring. The nano
!! preset's threshold on lattice_peak and the bar a local maximum must clear to count in spot_count. Set for recall:
!! on three labeled platinum trajectories 0.65 keeps 98-99% of the maintainer's selected classes and lets 27-56% of
!! his rejected ones through, which is the trade wanted for these time series (losing useful frames costs more than
!! carrying some junk into 3D); the balanced-accuracy optimum, 0.90-0.95, rejected faint-fringe classes the
!! maintainer keeps (nano_cavgs_rejection_report.md)
real,     parameter :: LAT_PEAK_MIN  = 0.65
!> the lattice grade of class averages of raw (not denoised) graphene-subtracted frames, which trajectory_solve3D
!! gates before denoising: the same rule as LAT_PEAK_MIN (the strictest threshold that keeps at least 98% of the
!! maintainer's selected classes on each of the three labeled trajectories, here rebuilt on the raw frames) gives
!! 0.52 (recall 0.995, 0.995, 0.985; 0.53 gave 0.977 on 71-6/traj04), where 0.65 keeps 0.984, 0.962, 0.939
!! (trajectory_solve3D.md, phase 1, ruling 6). The program passes it as thres
real,     parameter :: LAT_PEAK_MIN_RAW = 0.52
!> the model gate (gate 2 of trajectory_solve3D, model_gate_states): a class gate 1 kept is rejected when its
!! correlation to its best reprojection of the simulated start model lies more than this many robust spreads
!! (1.4826 x MAD) below the median of the kept classes. Set for recall on the three labeled platinum trajectories:
!! the smallest of 1.5, 2, 2.5, 3, 4 that keeps at least 98% of the maintainer's selected classes on each with the
!! lattice gate in front gave 3 on denoised averages (phase 0). On the raw averages the program judges, behind
!! LAT_PEAK_MIN_RAW, every K keeps the same recall (0.995, 0.995, 0.985); a gate set for recall breaks that tie
!! toward the permissive end, so 3 stays (ruling 9; 1.5 removed 2 of 13 classes of a clean simulated trajectory)
real,     parameter :: K_MODEL       = 3.0
integer,  parameter :: ND   = 105                !< sphere diameters 8 to 60 A, step 0.5
integer,  parameter :: NOFF = 9                  !< sub-pixel centre offsets -1 to 1 px, step 0.25
real(dp), parameter :: D0 = 8._dp, DSTEP = 0.5_dp, OFF0 = -1._dp, OFFSTEP = 0.25_dp

type :: cavg_quality_lattice
    private
    integer :: box = 0, c = 0, shmax = 0
    real    :: smpd = 0.
    integer,  allocatable :: sh(:,:)          !< (-c:c-1,-c:c-1) shell of (h,k)
    logical,  allocatable :: band(:,:), bandhalf(:,:)
    integer,  allocatable :: shell_count(:)   !< (0:shmax)
    real(dp), allocatable :: mask(:,:)        !< (box,box) the soft mask
    logical,  allocatable :: disc(:,:), annulus(:,:)
    real(dp), allocatable :: mc(:,:)          !< (NWIN, ND*NOFF*NOFF) centred sphere models
    real(dp), allocatable :: mbar(:), smm(:)
    logical :: exists = .false.
contains
    procedure :: new
    procedure :: extract
    procedure, private :: measure
    procedure :: kill
end type cavg_quality_lattice

contains

    !> the band, shells, mask and sphere-model bank for class averages of size box at smpd
    subroutine new( self, box, smpd )
        class(cavg_quality_lattice), intent(inout) :: self
        integer,                     intent(in)    :: box
        real,                        intent(in)    :: smpd
        logical, allocatable :: gkeep(:)
        logical, allocatable :: half(:,:)
        real(dp) :: rmsk, q, dq, x, dx, dy, mdl(NWIN)
        integer  :: h, k, i, j, s, id, jy, lx, im, e, iy, ix
        call self%kill
        if( mod(box, 2) /= 0 ) THROW_HARD('even box required; cavg_quality_lattice%new')
        self%box  = box
        self%smpd = smpd
        self%c    = box / 2
        allocate(self%sh(-self%c:self%c-1,-self%c:self%c-1))
        do k = -self%c, self%c - 1
            do h = -self%c, self%c - 1
                self%sh(h,k) = nint(hypot(real(h, dp), real(k, dp)))
            end do
        end do
        self%shmax = maxval(self%sh)
        allocate(self%shell_count(0:self%shmax), source=0)
        do k = -self%c, self%c - 1
            do h = -self%c, self%c - 1
                self%shell_count(self%sh(h,k)) = self%shell_count(self%sh(h,k)) + 1
            end do
        end do
        ! the band: 3 A to Nyquist, minus the graphene shells of calc_graphene_mask (shells 1..box/2)
        gkeep = calc_graphene_mask(box, smpd, [GRAPHENE_BAND1, GRAPHENE_BAND2, GRAPHENE_BAND3])
        allocate(self%band(-self%c:self%c-1,-self%c:self%c-1), self%bandhalf(-self%c:self%c-1,-self%c:self%c-1), &
            &half(-self%c:self%c-1,-self%c:self%c-1), source=.false.)
        do k = -self%c, self%c - 1
            do h = -self%c, self%c - 1
                s = self%sh(h,k)
                if( s >= ceiling(real(box) * smpd / BAND_LOW_RES - 1.e-6) .and. s <= self%c )then
                    if( s <= size(gkeep) ) self%band(h,k) = gkeep(s)
                endif
                half(h,k) = k > 0 .or. (k == 0 .and. h > 0)
            end do
        end do
        self%bandhalf = self%band .and. half
        ! soft mask (mask2D_soft / cosedge_r2_2d geometry) and the real-space disc and annulus of band presence
        rmsk = real(box, dp) / 2._dp - 2._dp
        allocate(self%mask(box,box), self%disc(box,box), self%annulus(box,box))
        do j = 1, box
            do i = 1, box
                self%mask(i,j)    = soft_edge(real((i-1-self%c)**2 + (j-1-self%c)**2, dp), box, rmsk)
                x = hypot(real(i-1-self%c, dp), real(j-1-self%c, dp))
                self%disc(i,j)    = x < real(box, dp) / 4._dp
                self%annulus(i,j) = x >= real(box, dp) / 4._dp .and. x < real(box, dp) / 2._dp
            end do
        end do
        ! projected-sphere models over the window, D outermost, then the row offset, then the column offset
        dq = 1._dp / (real(box, dp) * real(smpd, dp))
        allocate(self%mc(NWIN, ND*NOFF*NOFF), self%mbar(ND*NOFF*NOFF), self%smm(ND*NOFF*NOFF))
        im = 0
        do id = 1, ND
            do jy = 1, NOFF
                dy = OFF0 + real(jy - 1, dp) * OFFSTEP
                do lx = 1, NOFF
                    dx = OFF0 + real(lx - 1, dp) * OFFSTEP
                    im = im + 1
                    e  = 0
                    do iy = -SPOT_HALFWIN, SPOT_HALFWIN
                        do ix = -SPOT_HALFWIN, SPOT_HALFWIN
                            e = e + 1
                            q = hypot(real(ix, dp) - dx, real(iy, dp) - dy) * dq
                            mdl(e) = sphere2(DPI * q * (D0 + real(id - 1, dp) * DSTEP))
                        end do
                    end do
                    self%mbar(im)  = sum(mdl) / real(NWIN, dp)
                    self%mc(:,im)  = mdl - self%mbar(im)
                    self%smm(im)   = sum(self%mc(:,im)**2)
                end do
            end do
        end do
        self%exists = .true.
    end subroutine new

    !> raw lattice features of the class averages (columns in the LAT_* order); a valid class without a measurable
    !! spectrum (blank) is returned in blank(:). Parallel over classes.
    subroutine extract( self, imgs, valid, raw, blank )
        class(cavg_quality_lattice), intent(in)    :: self
        class(image),                intent(inout) :: imgs(:)
        logical,                     intent(in)    :: valid(:)
        real,                        intent(inout) :: raw(:,:)
        logical,                     intent(out)   :: blank(:)
        type(image), allocatable :: fimgs(:)
        real(dp),    allocatable :: p(:,:), ps(:,:), r(:,:), mu(:), x(:,:)
        real,        allocatable :: rbuf(:,:,:), fbuf(:,:,:), locmax(:,:)
        real    :: vals(5), locs(NLOC)
        integer :: n, nthr, ithr, i, ldim(3)
        logical :: invalid
        n = size(imgs)
        if( .not. self%exists ) THROW_HARD('bank not made; cavg_quality_lattice%extract')
        if( size(valid) /= n .or. size(blank) /= n ) THROW_HARD('mask sizes; cavg_quality_lattice%extract')
        if( size(raw,1) /= n .or. size(raw,2) /= LAT_NFEATS ) THROW_HARD('raw table shape; cavg_quality_lattice%extract')
        ldim = imgs(1)%get_ldim()
        if( ldim(1) /= self%box .or. ldim(2) /= self%box ) THROW_HARD('image size; cavg_quality_lattice%extract')
        allocate(locmax(n,NLOC), source=0.)
        raw   = 0.
        blank = .false.
        ! one image (an FFTW plan) per thread of the team, made outside the region
        nthr = 1
        !$ nthr = omp_get_max_threads()
        allocate(fimgs(nthr))
        do ithr = 1, nthr
            call fimgs(ithr)%new([self%box,self%box,1], self%smpd, wthreads=.false.)
        end do
        !$omp parallel num_threads(nthr) default(shared) private(i,ithr,p,ps,r,mu,x,rbuf,fbuf,vals,locs,invalid) &
        !$omp& proc_bind(close)
        allocate(p(-self%c:self%c-1,-self%c:self%c-1), ps(-self%c:self%c-1,-self%c:self%c-1), &
            &r(-self%c:self%c-1,-self%c:self%c-1), mu(0:self%shmax), x(self%box,self%box), &
            &rbuf(self%box,self%box,1), fbuf(self%box,self%box,1))
        !$omp do schedule(static)
        do i = 1, n
            if( .not. valid(i) ) cycle
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            call imgs(i)%get_rmat_sub(fbuf)
            call self%measure(fbuf(:,:,1), fimgs(ithr), p, ps, r, mu, x, rbuf, invalid, vals, locs)
            blank(i) = invalid
            raw(i,LAT_LATTICE_PEAK)  = vals(1)
            raw(i,LAT_BAND_PRESENCE) = vals(2)
            raw(i,LAT_SPOT_ARC)      = vals(3)
            raw(i,LAT_SPOT_DIAM)     = vals(4)
            raw(i,LAT_SPOT_SNR)      = vals(5)
            locmax(i,:)              = locs
        end do
        !$omp end do
        deallocate(p, ps, r, mu, x, rbuf, fbuf)
        !$omp end parallel
        do ithr = 1, nthr
            call fimgs(ithr)%kill
        end do
        deallocate(fimgs)
        ! spot_count: the local maxima of lattice grade
        do i = 1, n
            if( valid(i) .and. .not. blank(i) ) raw(i,LAT_SPOT_COUNT) = real(count(locmax(i,:) > 10.**LAT_PEAK_MIN))
        end do
    end subroutine extract

    !> one class average; the work arrays are the calling thread's. vals: lattice_peak, band_presence, spot_arc,
    !! spot_diam, spot_snr; locs the strongest local maxima of R; invalid for non-finite pixels or a blank image
    subroutine measure( self, frame, img, p, ps, r, mu, x, rbuf, invalid, vals, locs )
        class(cavg_quality_lattice), intent(in)    :: self
        real,                        intent(in)    :: frame(:,:)
        type(image),                 intent(inout) :: img
        real(dp),                    intent(inout) :: p(-self%c:,-self%c:), ps(-self%c:,-self%c:), r(-self%c:,-self%c:)
        real(dp),                    intent(inout) :: mu(0:), x(:,:)
        real,                        intent(inout) :: rbuf(:,:,:)
        logical,                     intent(out)   :: invalid
        real,                        intent(out)   :: vals(5), locs(NLOC)
        real,    pointer :: rmat(:,:,:)
        complex  :: cval
        real(dp) :: rmax, mdisc, mann, v, ymean, a, bgr, best, res, yy
        real(dp) :: y(NWIN), w(NWIN), yc(NWIN), u(2), cx, cy, sw, vr, vt, along, across, offx, offy
        integer  :: h, k, dh, dk, nloc, g1h, g1k, e, iy, ix, ibest, im, nd_disc, nd_ann
        logical  :: found
        vals    = 0.
        locs    = 0.
        invalid = .not. all(ieee_is_finite(frame))
        if( invalid ) return
        x = real(frame, dp)
        x = x - sum(x) / real(size(x), dp)
        ! band presence: the unmasked image band-passed, disc over annulus mean square
        rbuf(:,:,1) = real(x)
        call img%set_rmat(rbuf, .false.)
        call img%fft
        do k = -self%c, self%c - 1
            do h = 0, self%c
                if( .not. self%band(wrap(h, self%c, self%box), k) ) &
                    &call img%set_cmat_at(h + 1, k + 1 + merge(self%box, 0, k < 0), 1, cmplx(0., 0.))
            end do
        end do
        call img%ifft
        call img%get_rmat_ptr(rmat)
        mdisc   = 0._dp
        mann    = 0._dp
        nd_disc = 0
        nd_ann  = 0
        do k = 1, self%box
            do h = 1, self%box
                v = real(rmat(h,k,1), dp)**2
                if( self%disc(h,k) )then
                    mdisc   = mdisc + v
                    nd_disc = nd_disc + 1
                else if( self%annulus(h,k) )then
                    mann   = mann + v
                    nd_ann = nd_ann + 1
                endif
            end do
        end do
        mdisc = mdisc / real(nd_disc, dp)
        mann  = mann  / real(nd_ann,  dp)
        ! the masked spectrum, smoothed, ring means, ratio
        rbuf(:,:,1) = real(x * self%mask)
        call img%set_rmat(rbuf, .false.)
        call img%fft
        do k = -self%c, self%c - 1
            do h = -self%c, self%c - 1
                cval   = img%get_fcomp2D(h, k)
                p(h,k) = real(real(cval), dp)**2 + real(aimag(cval), dp)**2
            end do
        end do
        do k = -self%c, self%c - 1
            do h = -self%c, self%c - 1
                v = 0._dp
                do dk = -1, 1
                    do dh = -1, 1
                        v = v + p(wrap(h + dh, self%c, self%box), wrap(k + dk, self%c, self%box))
                    end do
                end do
                ps(h,k) = v / 9._dp
            end do
        end do
        mu = 0._dp
        do k = -self%c, self%c - 1
            do h = -self%c, self%c - 1
                mu(self%sh(h,k)) = mu(self%sh(h,k)) + ps(h,k)
            end do
        end do
        mu = mu / real(max(self%shell_count, 1), dp)
        rmax = 0._dp
        do k = -self%c, self%c - 1
            do h = -self%c, self%c - 1
                r(h,k) = 0._dp
                if( self%band(h,k) )then
                    r(h,k) = ps(h,k) / max(mu(self%sh(h,k)), 1.e-300_dp)
                    rmax   = max(rmax, r(h,k))
                endif
            end do
        end do
        if( rmax <= 0._dp .or. mdisc <= 0._dp .or. mann <= 0._dp )then
            invalid = .true.          ! a blank image: nothing to measure
            return
        endif
        ! local maxima of R on the half plane (at least their 8 periodic neighbours), the strongest first; g1 is the
        ! first of equal values in raster order (k outer, h inner)
        nloc  = 0
        found = .false.
        g1h   = 0
        g1k   = 0
        best  = -1._dp
        do k = -self%c, self%c - 1
            do h = -self%c, self%c - 1
                if( .not. self%bandhalf(h,k) ) cycle
                if( r(h,k) <= 0._dp ) cycle
                if( .not. is_local_max(r, h, k, self%c, self%box) ) cycle
                if( r(h,k) > best )then
                    best  = r(h,k)
                    g1h   = h
                    g1k   = k
                    found = .true.
                endif
                call insert_desc(locs, nloc, real(r(h,k)))
            end do
        end do
        vals(1) = real(log10(rmax))
        vals(2) = real(log10(mdisc / mann))
        if( .not. found ) return
        ! the 15 x 15 window about g1: raw power for the fit, Ps minus the ring mean for the arc
        e = 0
        do iy = -SPOT_HALFWIN, SPOT_HALFWIN
            do ix = -SPOT_HALFWIN, SPOT_HALFWIN
                e = e + 1
                h = wrap(g1h + ix, self%c, self%box)
                k = wrap(g1k + iy, self%c, self%box)
                y(e) = p(h,k)
                w(e) = max(ps(h,k) - mu(self%sh(h,k)), 0._dp)
            end do
        end do
        ! spot_arc: second moments across and along the radius of g1 about the weighted centroid
        sw = sum(w)
        if( sw > 0._dp )then
            cx = 0._dp
            cy = 0._dp
            e  = 0
            do iy = -SPOT_HALFWIN, SPOT_HALFWIN
                do ix = -SPOT_HALFWIN, SPOT_HALFWIN
                    e  = e + 1
                    cx = cx + w(e) * real(ix, dp)
                    cy = cy + w(e) * real(iy, dp)
                end do
            end do
            cx = cx / sw
            cy = cy / sw
            u  = [real(g1h, dp), real(g1k, dp)] / hypot(real(g1h, dp), real(g1k, dp))
            vr = 0._dp
            vt = 0._dp
            e  = 0
            do iy = -SPOT_HALFWIN, SPOT_HALFWIN
                do ix = -SPOT_HALFWIN, SPOT_HALFWIN
                    e      = e + 1
                    offx   = real(ix, dp) - cx
                    offy   = real(iy, dp) - cy
                    along  = offx * u(1) + offy * u(2)
                    across = -offx * u(2) + offy * u(1)
                    vr = vr + w(e) * along**2
                    vt = vt + w(e) * across**2
                end do
            end do
            if( vr > 0._dp .and. vt > 0._dp ) vals(3) = real(0.5_dp * log10(vt / vr))
        endif
        ! spot_diam, spot_snr: least squares a * model + b, a >= 0, over the model bank
        ymean = sum(y) / real(NWIN, dp)
        yc    = y - ymean
        yy    = sum(yc**2)
        ibest = 1
        best  = huge(1._dp)
        a     = 0._dp
        do im = 1, size(self%smm)
            v   = max(dot_product(self%mc(:,im), yc) / self%smm(im), 0._dp)
            res = yy - v * v * self%smm(im)
            if( res < best )then
                best  = res
                ibest = im
                a     = v
            endif
        end do
        bgr     = ymean - a * self%mbar(ibest)
        vals(4) = real(D0 + real((ibest - 1) / (NOFF * NOFF), dp) * DSTEP)
        vals(5) = real(log10(max(a / max(bgr, 1.e-300_dp), 1.e-6_dp)))
    end subroutine measure

    subroutine kill( self )
        class(cavg_quality_lattice), intent(inout) :: self
        if( allocated(self%sh)          ) deallocate(self%sh)
        if( allocated(self%band)        ) deallocate(self%band)
        if( allocated(self%bandhalf)    ) deallocate(self%bandhalf)
        if( allocated(self%shell_count) ) deallocate(self%shell_count)
        if( allocated(self%mask)        ) deallocate(self%mask)
        if( allocated(self%disc)        ) deallocate(self%disc)
        if( allocated(self%annulus)     ) deallocate(self%annulus)
        if( allocated(self%mc)          ) deallocate(self%mc)
        if( allocated(self%mbar)        ) deallocate(self%mbar)
        if( allocated(self%smm)         ) deallocate(self%smm)
        self%exists = .false.
    end subroutine kill

    !> the model gate: among the classes with states_in > 0 (gate 1's survivors), the ones whose correlation to the
    !! start model lies below med - nsd * sd (med and sd = 1.4826 x MAD over the survivors' correlations) get state 0;
    !! every other state is returned as given. With no survivor nothing changes and med, sd and cut are 0
    subroutine model_gate_states( corrs, states_in, nsd, states_out, med, sd, cut )
        real,    intent(in)  :: corrs(:)
        integer, intent(in)  :: states_in(:)
        real,    intent(in)  :: nsd
        integer, intent(out) :: states_out(size(states_in))
        real,    intent(out) :: med, sd, cut
        real, allocatable :: surv(:)
        if( size(corrs) /= size(states_in) ) THROW_HARD('correlations and states differ in size; model_gate_states')
        states_out = states_in
        med = 0.
        sd  = 0.
        cut = 0.
        surv = pack(corrs, states_in > 0)
        if( size(surv) == 0 ) return
        med = median(surv)
        sd  = mad_gau(surv, med)
        cut = med - nsd * sd
        where( states_in > 0 .and. corrs < cut ) states_out = 0
    end subroutine model_gate_states

    !> R(h,k) at least as large as its 8 periodic neighbours
    pure logical function is_local_max( r, h, k, c, box )
        integer,  intent(in) :: h, k, c, box
        real(dp), intent(in) :: r(-c:,-c:)
        integer :: dh, dk
        is_local_max = .true.
        do dk = -1, 1
            do dh = -1, 1
                if( dh == 0 .and. dk == 0 ) cycle
                if( r(wrap(h + dh, c, box), wrap(k + dk, c, box)) > r(h,k) ) is_local_max = .false.
            end do
        end do
    end function is_local_max

    !> the descending list of the strongest local maxima; an equal value goes after the earlier ones
    pure subroutine insert_desc( locs, nloc, val )
        real,    intent(inout) :: locs(:)
        integer, intent(inout) :: nloc
        real,    intent(in)    :: val
        integer :: pos
        if( nloc == size(locs) )then
            if( val <= locs(nloc) ) return
        else
            nloc = nloc + 1
        endif
        pos = nloc
        do while( pos > 1 )
            if( locs(pos-1) >= val ) exit
            locs(pos) = locs(pos-1)
            pos = pos - 1
        end do
        locs(pos) = val
    end subroutine insert_desc

    !> frequency index i on the periodic grid -c..c-1
    pure elemental integer function wrap( i, c, box )
        integer, intent(in) :: i, c, box
        wrap = modulo(i + c, box) - c
    end function wrap

    !> cosedge_r2_2d with box and mskrad, as mask2D_soft applies it (unit where above 0.9999)
    pure real(dp) function soft_edge( r2, box, mskrad )
        real(dp), intent(in) :: r2, mskrad
        integer,  intent(in) :: box
        real(dp) :: maxrad, width, r0
        integer  :: minlen
        minlen = min(nint(2._dp * (mskrad + real(COSMSKHALFWIDTH, dp))), box)
        maxrad = 0.5_dp * real(minlen, dp)
        width  = 2._dp * (maxrad - mskrad)
        r0     = maxrad - width
        if( r2 >= maxrad * maxrad )then
            soft_edge = 0._dp
        else if( r2 <= r0 * r0 )then
            soft_edge = 1._dp
        else
            soft_edge = 0.5_dp * (cos((sqrt(r2) - r0) / width * DPI) + 1._dp)
        endif
        if( soft_edge > 0.9999_dp ) soft_edge = 1._dp
    end function soft_edge

    !> the projected-sphere shape transform [3 (sin x - x cos x) / x^3]^2
    pure real(dp) function sphere2( x )
        real(dp), intent(in) :: x
        if( x > 1.e-6_dp )then
            sphere2 = (3._dp * (sin(x) - x * cos(x)) / x**3)**2
        else
            sphere2 = 1._dp
        endif
    end function sphere2

end module simple_cavg_quality_lattice
