!@descr: the per-frame Fourier measurements of the frame-quality bank: lattice peak, band maximum, spot list, spot shape and size, band presence
! Per frame: mean removed; band presence on the unmasked frame; then the frame soft-masked (mask2D_soft geometry),
! its full-plane power P = |F|^2 (F scaled 1/box^2), smoothed by a periodic 3 x 3 box (Ps), ring means of Ps on every
! shell (mu), the smoothed ratio R = Ps / mu over the band (3 A to Nyquist minus the graphene shells). The reference
! implementation of these definitions is the run's numpy frame_feats.py.
module single_frame_quality_feats
!$ use omp_lib
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs,               only: dp, GRAPHENE_BAND1, GRAPHENE_BAND2, GRAPHENE_BAND3, COSMSKHALFWIDTH, DPI
use simple_error,              only: simple_exception
use simple_image,              only: image
use simple_math_ft,            only: calc_graphene_mask
use single_frame_quality_types, only: frame_quality_table, FQ_NLOC, FQ_SPOT_HALFWIN, FQ_BAND_LOW_RES, &
    I_LATTICE_PEAK, I_BAND_MAX, I_SPOT_ARC, I_SPOT_DIAM, I_SPOT_SNR, I_BAND_PRESENCE
implicit none
private
#include "simple_local_flags.inc"

public :: frame_quality_bank

integer,  parameter :: NWIN    = (2 * FQ_SPOT_HALFWIN + 1)**2
integer,  parameter :: ND      = 105        !< sphere diameters 8 to 60 A, step 0.5
integer,  parameter :: NOFF    = 9          !< sub-pixel centre offsets -1 to 1 px, step 0.25
real(dp), parameter :: D0 = 8._dp, DSTEP = 0.5_dp, OFF0 = -1._dp, OFFSTEP = 0.25_dp

type :: frame_quality_bank
    private
    integer :: box = 0, c = 0, shmax = 0
    real    :: smpd = 0.
    integer, public :: sh_lo = 0, sh_hi = 0
    integer,  allocatable :: sh(:,:)          !< (-c:c-1,-c:c-1) shell of (h,k)
    logical,  allocatable :: band(:,:), half(:,:), bandhalf(:,:)
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
    procedure :: white_noise
    procedure :: kill
end type frame_quality_bank

contains

    !> the band, shells, mask and model bank for frames of size box at smpd; mskdiam (A) > 0 sets the mask radius,
    !! else the box's inscribed circle (box/2 - 2) as trajectory_denoise uses it
    subroutine new( self, box, smpd, mskdiam )
        class(frame_quality_bank), intent(inout) :: self
        integer,                   intent(in)    :: box
        real,                      intent(in)    :: smpd
        real, optional,            intent(in)    :: mskdiam
        logical, allocatable :: gkeep(:)
        real(dp) :: rmsk, q, dq, x, dx, dy, mdl(NWIN)
        integer  :: h, k, i, j, s, id, jy, lx, im, e, iy, ix
        call self%kill
        if( mod(box, 2) /= 0 ) THROW_HARD('even box required; frame_quality_bank%new')
        self%box  = box
        self%smpd = smpd
        self%c    = box / 2
        allocate(self%sh(-self%c:self%c-1,-self%c:self%c-1))
        do k = -self%c, self%c - 1
            do h = -self%c, self%c - 1
                self%sh(h,k) = shell_of(h, k)
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
        self%sh_lo = ceiling(real(box) * smpd / FQ_BAND_LOW_RES - 1.e-6)
        self%sh_hi = self%c
        gkeep = calc_graphene_mask(box, smpd, [GRAPHENE_BAND1, GRAPHENE_BAND2, GRAPHENE_BAND3])
        allocate(self%band(-self%c:self%c-1,-self%c:self%c-1), self%half(-self%c:self%c-1,-self%c:self%c-1), &
            &self%bandhalf(-self%c:self%c-1,-self%c:self%c-1), source=.false.)
        do k = -self%c, self%c - 1
            do h = -self%c, self%c - 1
                s = self%sh(h,k)
                if( s >= self%sh_lo .and. s <= self%sh_hi )then
                    if( s <= size(gkeep) ) self%band(h,k) = gkeep(s)
                endif
                self%half(h,k) = k > 0 .or. (k == 0 .and. h > 0)
            end do
        end do
        self%bandhalf = self%band .and. self%half
        ! soft mask (mask2D_soft / cosedge_r2_2d geometry) and the real-space disc and annulus of band presence
        rmsk = real(box, dp) / 2._dp - 2._dp
        if( present(mskdiam) )then
            if( mskdiam > 0. ) rmsk = min(rmsk, 0.5_dp * real(mskdiam, dp) / real(smpd, dp))
        endif
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
                    do iy = -FQ_SPOT_HALFWIN, FQ_SPOT_HALFWIN
                        do ix = -FQ_SPOT_HALFWIN, FQ_SPOT_HALFWIN
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

    !> raw features of frames(:,:,i) into rows i0, i0+1, ... of tab (lattice_peak, band_max, spot_arc, spot_diam,
    !! spot_snr, band_presence, the local maxima and the invalid flag); parallel over frames
    subroutine extract( self, frames, tab, i0 )
        class(frame_quality_bank), intent(in)    :: self
        real,                      intent(in)    :: frames(:,:,:)
        type(frame_quality_table), intent(inout) :: tab
        integer,                   intent(in)    :: i0
        type(image), allocatable :: imgs(:)
        real(dp),    allocatable :: p(:,:), ps(:,:), r(:,:), mu(:), x(:,:)
        real,        allocatable :: rbuf(:,:,:)
        real    :: vals(6), locs(FQ_NLOC)
        integer :: n, nthr, ithr, i, irow, shell
        logical :: invalid
        n = size(frames, 3)
        if( .not. self%exists ) THROW_HARD('bank not made; frame_quality_bank%extract')
        if( size(frames,1) /= self%box .or. size(frames,2) /= self%box ) THROW_HARD('frame size; frame_quality_bank%extract')
        if( i0 < 1 .or. i0 + n - 1 > tab%n ) THROW_HARD('rows out of range; frame_quality_bank%extract')
        ! one image (an FFTW plan) per thread of the team, made outside the region
        nthr = 1
        !$ nthr = omp_get_max_threads()
        allocate(imgs(nthr))
        do ithr = 1, nthr
            call imgs(ithr)%new([self%box,self%box,1], self%smpd, wthreads=.false.)
        end do
        !$omp parallel num_threads(nthr) default(shared) private(i,ithr,irow,p,ps,r,mu,x,rbuf,vals,locs,invalid,shell) &
        !$omp& proc_bind(close)
        allocate(p(-self%c:self%c-1,-self%c:self%c-1), ps(-self%c:self%c-1,-self%c:self%c-1), &
            &r(-self%c:self%c-1,-self%c:self%c-1), mu(0:self%shmax), x(self%box,self%box), &
            &rbuf(self%box,self%box,1))
        !$omp do schedule(static)
        do i = 1, n
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            call self%measure(frames(:,:,i), imgs(ithr), p, ps, r, mu, x, rbuf, invalid, vals, locs, shell)
            irow = i0 + i - 1
            tab%invalid(irow)             = invalid
            tab%raw(irow,I_LATTICE_PEAK)  = vals(1)
            tab%raw(irow,I_BAND_MAX)      = vals(2)
            tab%raw(irow,I_BAND_PRESENCE) = vals(3)
            tab%raw(irow,I_SPOT_ARC)      = vals(4)
            tab%raw(irow,I_SPOT_DIAM)     = vals(5)
            tab%raw(irow,I_SPOT_SNR)      = vals(6)
            tab%locmax(irow,:)            = locs
            tab%spot_shell(irow)          = shell
        end do
        !$omp end do
        deallocate(p, ps, r, mu, x, rbuf)
        !$omp end parallel
        do ithr = 1, nthr
            call imgs(ithr)%kill
        end do
        deallocate(imgs)
    end subroutine extract

    !> one frame; the work arrays are the calling thread's. vals: lattice_peak, band_max, band_presence, spot_arc,
    !! spot_diam, spot_snr; locs the strongest local maxima of R; shell that of the strongest spot; invalid for
    !! non-finite pixels or a blank frame
    subroutine measure( self, frame, img, p, ps, r, mu, x, rbuf, invalid, vals, locs, shell )
        class(frame_quality_bank), intent(in)    :: self
        real,                      intent(in)    :: frame(:,:)
        type(image),               intent(inout) :: img
        real(dp),                  intent(inout) :: p(-self%c:,-self%c:), ps(-self%c:,-self%c:), r(-self%c:,-self%c:)
        real(dp),                  intent(inout) :: mu(0:), x(:,:)
        real,                      intent(inout) :: rbuf(:,:,:)
        logical,                   intent(out)   :: invalid
        real,                      intent(out)   :: vals(6), locs(FQ_NLOC)
        integer,                   intent(out)   :: shell
        real,    pointer :: rmat(:,:,:)
        complex  :: cval
        real(dp) :: rmax, bmax, mdisc, mann, v, ymean, a, bgr, best, res, yy
        real(dp) :: y(NWIN), w(NWIN), yc(NWIN), u(2), cx, cy, sw, vr, vt, along, across, offx, offy
        integer  :: h, k, dh, dk, nloc, g1h, g1k, e, iy, ix, ibest, im, nd_disc, nd_ann
        logical  :: found
        vals    = 0.
        locs    = 0.
        shell   = 0
        invalid = .not. all(ieee_is_finite(frame))
        if( invalid ) return
        x = real(frame, dp)
        x = x - sum(x) / real(size(x), dp)
        ! band presence: the unmasked frame band-passed, disc over annulus mean square
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
        bmax = 0._dp
        do k = -self%c, self%c - 1
            do h = -self%c, self%c - 1
                r(h,k) = 0._dp
                if( self%band(h,k) )then
                    r(h,k) = ps(h,k) / max(mu(self%sh(h,k)), 1.e-300_dp)
                    rmax   = max(rmax, r(h,k))
                    bmax   = max(bmax, ps(h,k))
                endif
            end do
        end do
        if( rmax <= 0._dp .or. mdisc <= 0._dp .or. mann <= 0._dp )then
            invalid = .true.          ! a blank frame: nothing to measure
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
        vals(2) = real(log10(bmax))
        vals(3) = real(log10(mdisc / mann))
        if( .not. found ) return
        shell = self%sh(g1h,g1k)
        ! the 15 x 15 window about g1: raw power for the fit, Ps minus the ring mean for the arc
        e = 0
        do iy = -FQ_SPOT_HALFWIN, FQ_SPOT_HALFWIN
            do ix = -FQ_SPOT_HALFWIN, FQ_SPOT_HALFWIN
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
            do iy = -FQ_SPOT_HALFWIN, FQ_SPOT_HALFWIN
                do ix = -FQ_SPOT_HALFWIN, FQ_SPOT_HALFWIN
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
            do iy = -FQ_SPOT_HALFWIN, FQ_SPOT_HALFWIN
                do ix = -FQ_SPOT_HALFWIN, FQ_SPOT_HALFWIN
                    e      = e + 1
                    offx   = real(ix, dp) - cx
                    offy   = real(iy, dp) - cy
                    along  = offx * u(1) + offy * u(2)
                    across = -offx * u(2) + offy * u(1)
                    vr = vr + w(e) * along**2
                    vt = vt + w(e) * across**2
                end do
            end do
            if( vr > 0._dp .and. vt > 0._dp ) vals(4) = real(0.5_dp * log10(vt / vr))
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
        vals(5) = real(D0 + real((ibest - 1) / (NOFF * NOFF), dp) * DSTEP)
        vals(6) = real(log10(max(a / max(bgr, 1.e-300_dp), 1.e-6_dp)))
    end subroutine measure

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

    !> median lattice_peak of nframes white Gaussian frames through the same pipeline (deterministic stream)
    real function white_noise( self, nframes )
        class(frame_quality_bank), intent(in) :: self
        integer,                   intent(in) :: nframes
        type(frame_quality_table) :: tab
        real,    allocatable :: frames(:,:,:), lp(:)
        integer(kind=8) :: state
        integer :: i, j, l, m
        real(dp) :: u1, u2
        allocate(frames(self%box,self%box,nframes))
        state = 20261027_8
        do l = 1, nframes
            do j = 1, self%box
                do i = 1, self%box
                    u1 = next_uniform(state)
                    u2 = next_uniform(state)
                    frames(i,j,l) = real(sqrt(-2._dp * log(u1)) * cos(2._dp * DPI * u2))
                end do
            end do
        end do
        call tab%new(nframes, self%box, self%smpd)
        call self%extract(frames, tab, 1)
        lp = pack(tab%raw(:,I_LATTICE_PEAK), .not. tab%invalid)
        m  = size(lp)
        call sort_real(lp)
        if( mod(m, 2) == 1 )then
            white_noise = lp(m/2 + 1)
        else
            white_noise = 0.5 * (lp(m/2) + lp(m/2 + 1))
        endif
        call tab%kill
    end function white_noise

    subroutine kill( self )
        class(frame_quality_bank), intent(inout) :: self
        if( allocated(self%sh)          ) deallocate(self%sh)
        if( allocated(self%band)        ) deallocate(self%band)
        if( allocated(self%half)        ) deallocate(self%half)
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

    pure integer function shell_of( h, k )
        integer, intent(in) :: h, k
        shell_of = nint(hypot(real(h, dp), real(k, dp)))
    end function shell_of

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

    real(dp) function next_uniform( state )
        integer(kind=8), intent(inout) :: state
        state = modulo(state * 48271_8, 2147483647_8)
        next_uniform = (real(state, dp) + 0.5_dp) / 2147483647._dp
    end function next_uniform

    subroutine sort_real( v )
        real, intent(inout) :: v(:)
        integer :: i, j
        real    :: t
        do i = 2, size(v)
            t = v(i)
            j = i - 1
            do while( j >= 1 )
                if( v(j) <= t ) exit
                v(j+1) = v(j)
                j = j - 1
            end do
            v(j+1) = t
        end do
    end subroutine sort_real

end module single_frame_quality_feats
