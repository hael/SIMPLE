!@descr: the pose student's input: a particle's CTF-multiplied polar half ring in the band packed by the shell weights, its conjugate wrap, the shift and rotation augmentation, and the feature cache
! The tensor is what the search's objective sees (doc/implementation_notes/completed/learned_pose_state_assignment.md,
! section 3): the half ring of the calculator's rotation grid, the CTF left in the data, the real and imaginary
! parts scaled by sqrt(w_k) so that a Euclidean product of two packed tensors is the objective's weighted
! Hermitian product, the whole scaled by the particle's weighted power, the ring mean of the CTF squared beside it.
module simple_polar_feats
use, intrinsic :: iso_fortran_env, only: int64, int32
use simple_defs,         only: sp, dp, DPI, DTWOPI, TINY, logfhandle
use simple_error,        only: simple_exception
use simple_string,       only: string
use simple_fileio,       only: fopen, fclose, fileiochk, del_file
use simple_syslib,       only: file_exists, simple_atomic_replace, fs_avail_bytes, get_process_id, simple_abspath
use simple_string_utils, only: int2str
use simple_polarft_calc, only: polarft_calc
implicit none

public :: polar_feats, polar_feat, polar_feats_cache
public :: polar_feats_contract, polar_feats_hash_ints, polar_feats_net_pftsz, POLAR_FEATS_NHEAD
private
#include "simple_local_flags.inc"

!> the cache's contract: one integer(int64) per slot, record 1 of the cache file
integer,        parameter :: POLAR_FEATS_NHEAD   = 16
integer(int64), parameter :: FEATS_CACHE_MAGIC   = 6180339887_int64
integer(int64), parameter :: FEATS_CACHE_VERSION = 2_int64   !< 2: the tensors on the network's capped half ring
integer, parameter :: IH_MAGIC = 1, IH_VERSION = 2, IH_PATH = 3, IH_MTIME = 4, IH_ROWS = 5, IH_BOX = 6, IH_BOXCROP = 7, &
    &IH_SMPD = 8, IH_NSEL = 9, IH_SELHASH = 10, IH_HIGHROW = 11, IH_PFTSZ = 12, IH_KFROM = 13, IH_KTO = 14, IH_WHASH = 15, &
    &IH_CTF = 16
integer(int64), parameter :: FREE_SPACE_DENOM = 4_int64  !< the cache takes at most a quarter of the free space
real(dp),       parameter :: NET_ANGLE_OVERSAMPLING = 2._dp !< the network's half ring holds at most this many times pi k_max angles, the outer ring's sampling limit
! the overflow-safe rolling hashes of the selection and the weights: two Lehmer generators on primes below 2^31
integer(int64), parameter :: HASH_M1 = 2147483647_int64, HASH_A1 = 48271_int64
integer(int64), parameter :: HASH_M2 = 2147483629_int64, HASH_A2 = 16807_int64

!> one particle's input: the packed half ring, the ring-averaged CTF squared beside it, the side scalar
type :: polar_feat
    real(sp), allocatable :: x(:,:,:)       !< (pftsz, kfromto(1):kfromto(2), 2): Re and Im of the CTF-multiplied ring, each scaled by sqrt(w_k) and by the band's scale
    real(sp), allocatable :: ctf2(:)        !< (kfromto(1):kfromto(2)): the ring mean of the CTF squared, the first layer's normaliser
    real(sp)              :: tazim  = 0.    !< the azimuthal variation of the CTF squared: the mean over rings of its coefficient of variation
    real(sp)              :: logpow = 0.    !< the side scalar: log of the weighted power
    real(dp)              :: wsq    = 0.d0  !< the weighted power sum_k w_k sum_p |X_kp|^2 over the half ring (the calculator's wsqsum)
    logical               :: iseven = .true.
  contains
    procedure :: kill => feat_kill
end type polar_feat

!> the geometry and the shell weights the tensors are built with: the calculator's half ring (pftsz_in) and band,
!! the network's half ring (pftsz, at most the band's cap), the weights k/sigma2(k) per half frozen in the model, and
!! the calculator's shift arguments for the augmentation
type :: polar_feats
    private
    integer               :: pftsz = 0, pftsz_in = 0, nrots = 0, box = 0, nk = 0
    integer               :: kfromto(2) = 0
    real(sp), allocatable :: wk(:,:)        !< (kfromto(1):kfromto(2), 2) the shell weights k/sigma2(k) per half
    real(sp), allocatable :: sqw(:,:)       !< their square roots
    real(dp), allocatable :: ax(:,:), ay(:,:)       !< (pftsz, kfromto(1):kfromto(2)) the phase ramp's shift arguments on the network's half ring
    real(dp), allocatable :: ax_in(:,:), ay_in(:,:) !< (pftsz_in, kfromto(1):kfromto(2)) the same on the calculator's half ring
    real(sp), allocatable :: rs_re(:,:), rs_im(:,:) !< (pftsz, pftsz_in) the band-limited resampling of the real and imaginary channels
    logical               :: l_resample = .false.
    logical               :: exists = .false.
  contains
    procedure :: new
    procedure :: get_pftsz
    procedure :: get_pftsz_in
    procedure :: get_nrots
    procedure :: get_kfromto
    procedure :: alloc_feat
    procedure :: build
    procedure :: build_from_pft
    procedure :: shift_pft
    procedure :: shift_feat
    procedure :: rotate
    procedure :: pad
    procedure :: kill
end type polar_feats

!> The feature cache: one direct-access record per particle, record 1 the contract (polar_feats_contract). Written
!! under a process token and renamed into place after the contract, so a reader finds no cache or a complete one;
!! refused when the project, the selection, the geometry or the weights differ; at most a quarter of the free space
type :: polar_feats_cache
    private
    type(string)   :: fname, tmpname
    integer        :: funit   = -1
    integer        :: pftsz   = 0, nk = 0, kfrom = 0, reclen = 0
    integer(int64) :: highrow = 0_int64
    logical        :: l_open    = .false.
    logical        :: l_writing = .false.
  contains
    procedure :: new       => cache_new
    procedure :: create    => cache_create
    procedure :: put       => cache_put
    procedure :: seal      => cache_seal
    procedure :: probe     => cache_probe
    procedure :: open_read => cache_open_read
    procedure :: get       => cache_get
    procedure :: close     => cache_close
    procedure :: delete    => cache_delete
    procedure :: kill      => cache_kill
end type polar_feats_cache

contains

    ! ===== THE FEATURES =====

    !> the geometry of the calculator's half ring (pftsz_in angles over the half turn, the band kfromto, the box the
    !! shifts are measured in) and the shell weights wk(kfromto(1):kfromto(2), 2) of the two halves; the network's
    !! half ring is the calculator's capped by polar_feats_net_pftsz
    subroutine new( self, pftsz_in, kfromto, box, wk )
        class(polar_feats), intent(inout) :: self
        integer,            intent(in)    :: pftsz_in, kfromto(2), box
        real(sp),           intent(in)    :: wk(kfromto(1):kfromto(2),2)
        call self%kill
        if( pftsz_in < 1 ) THROW_HARD('pftsz must be positive; polar_feats%new')
        if( kfromto(1) < 1 .or. kfromto(2) < kfromto(1) ) THROW_HARD('invalid band; polar_feats%new')
        if( box < 2 ) THROW_HARD('box must be at least 2; polar_feats%new')
        if( any(wk < 0.) ) THROW_HARD('negative shell weights; polar_feats%new')
        self%pftsz_in   = pftsz_in
        self%pftsz      = polar_feats_net_pftsz(pftsz_in, kfromto(2))
        self%l_resample = self%pftsz < pftsz_in
        self%nrots      = 2 * self%pftsz
        self%box        = box
        self%kfromto    = kfromto
        self%nk         = kfromto(2) - kfromto(1) + 1
        allocate(self%wk(kfromto(1):kfromto(2),2), source=wk)
        allocate(self%sqw(kfromto(1):kfromto(2),2), source=sqrt(wk))
        call ramp_args(self%pftsz, self%ax, self%ay)
        call ramp_args(pftsz_in, self%ax_in, self%ay_in)
        if( self%l_resample ) call resampling_matrices(pftsz_in, self%pftsz, self%rs_re, self%rs_im)
        self%exists = .true.

      contains

        !> the calculator's shift arguments (argtransf) on a half ring of n angles: 2 pi / box times the polar
        !! coordinate (h, k) of each sample, h = k sin(theta_p), kc = -k cos(theta_p), theta_p = (p-1) pi / n
        subroutine ramp_args( n, ax, ay )
            integer,               intent(in)  :: n
            real(dp), allocatable, intent(out) :: ax(:,:), ay(:,:)
            real(dp) :: theta, a
            integer  :: p, k
            allocate(ax(n,kfromto(1):kfromto(2)), ay(n,kfromto(1):kfromto(2)))
            a = DTWOPI / real(box,dp)
            do p = 1,n
                theta = real(p-1,dp) * DPI / real(n,dp)
                do k = kfromto(1),kfromto(2)
                    ax(p,k) =  real(k,dp) * sin(theta) * a
                    ay(p,k) = -real(k,dp) * cos(theta) * a
                enddo
            enddo
        end subroutine ramp_args

    end subroutine new

    !> The network's half ring for a calculator's of pftsz_in angles and the band's last ring kmax: at most
    !! NET_ANGLE_OVERSAMPLING pi kmax angles (even), the outer ring's angular sampling limit being pi kmax
    pure integer function polar_feats_net_pftsz( pftsz_in, kmax ) result( pftsz )
        integer, intent(in) :: pftsz_in, kmax
        integer :: cap
        cap   = 2 * ceiling(NET_ANGLE_OVERSAMPLING * DPI * real(kmax,dp) / 2._dp)
        pftsz = min(pftsz_in, cap)
    end function polar_feats_net_pftsz

    !> The band-limited resampling of a Friedel-extended ring from n to m half-ring angles (m < n): the full ring's
    !! angular harmonics |j| < m kept, as the Dirichlet kernel D. With f(theta + pi) = conj(f(theta)), entry q is
    !! sum_p (D(dq_p) + D(dq_p - pi)) Re f_p / 2n for the real channel and (D(dq_p) - D(dq_p - pi)) for the imaginary
    subroutine resampling_matrices( n, m, rs_re, rs_im )
        integer,               intent(in)  :: n, m
        real(sp), allocatable, intent(out) :: rs_re(:,:), rs_im(:,:)
        real(dp) :: d0, d1, phi
        integer  :: p, q
        allocate(rs_re(m,n), rs_im(m,n))
        do p = 1,n
            do q = 1,m
                phi = real(q-1,dp) * DPI / real(m,dp) - real(p-1,dp) * DPI / real(n,dp)
                d0  = dirichlet(phi)
                d1  = dirichlet(phi - DPI)
                rs_re(q,p) = real((d0 + d1) / real(2*n,dp))
                rs_im(q,p) = real((d0 - d1) / real(2*n,dp))
            enddo
        enddo

      contains

        !> sum over |j| <= m-1 of exp(i j phi)
        real(dp) function dirichlet( x )
            real(dp), intent(in) :: x
            real(dp) :: sh
            sh = sin(0.5_dp * x)
            if( abs(sh) < 1.e-12_dp )then
                dirichlet = real(2*m-1,dp)      ! x a multiple of 2 pi: every term is one
            else
                dirichlet = sin((real(m,dp) - 0.5_dp) * x) / sh
            endif
        end function dirichlet

    end subroutine resampling_matrices

    pure integer function get_pftsz( self )
        class(polar_feats), intent(in) :: self
        get_pftsz = self%pftsz
    end function get_pftsz

    pure integer function get_pftsz_in( self )
        class(polar_feats), intent(in) :: self
        get_pftsz_in = self%pftsz_in
    end function get_pftsz_in

    pure integer function get_nrots( self )
        class(polar_feats), intent(in) :: self
        get_nrots = self%nrots
    end function get_nrots

    pure function get_kfromto( self ) result( kfromto )
        class(polar_feats), intent(in) :: self
        integer :: kfromto(2)
        kfromto = self%kfromto
    end function get_kfromto

    !> a feature of this geometry, zeroed
    subroutine alloc_feat( self, feat )
        class(polar_feats), intent(in)    :: self
        type(polar_feat),   intent(inout) :: feat
        call feat%kill
        allocate(feat%x(self%pftsz,self%kfromto(1):self%kfromto(2),2), source=0.)
        allocate(feat%ctf2(self%kfromto(1):self%kfromto(2)), source=0.)
    end subroutine alloc_feat

    !> the particle's tensor from the calculator: its half ring over the band, its CTF matrix and its half, read
    !! through the calculator's getters; the one routine training and prediction build their input with
    subroutine build( self, pftc, iptcl, feat, normalise )
        class(polar_feats),  intent(in)    :: self
        class(polarft_calc), intent(in)    :: pftc
        integer,             intent(in)    :: iptcl
        type(polar_feat),    intent(inout) :: feat
        logical, optional,   intent(in)    :: normalise
        complex(sp), allocatable :: pft(:,:)
        real(sp),    allocatable :: ctfmat(:,:)
        integer :: kfromto(2)
        if( .not. self%exists ) THROW_HARD('polar_feats not constructed; build')
        kfromto = pftc%get_kfromto()
        if( pftc%get_pftsz() /= self%pftsz_in .or. any(kfromto /= self%kfromto) )then
            write(logfhandle,*) 'calculator pftsz, kfromto: ', pftc%get_pftsz(), kfromto
            write(logfhandle,*) 'features   pftsz, kfromto: ', self%pftsz_in, self%kfromto
            THROW_HARD('the calculator does not match the features'' geometry; polar_feats%build')
        endif
        allocate(pft(self%pftsz_in,self%kfromto(1):self%kfromto(2)), ctfmat(self%pftsz_in,self%kfromto(1):self%kfromto(2)))
        call pftc%get_ptcl_pft(iptcl, pft)
        call pftc%get_ptcl_ctfmat(iptcl, ctfmat)
        call self%build_from_pft(pft, ctfmat, pftc%ptcl_iseven(iptcl), feat, normalise)
        deallocate(pft, ctfmat)
    end subroutine build

    !> The tensor from a calculator's half ring and its CTF matrix: the CTF-multiplied ring's real and imaginary
    !! parts, each scaled by sqrt(w_k) of the particle's half, resampled onto the network's half ring, the whole
    !! divided by the square root of the weighted power (normalise, the default) so that the band has unit weighted
    !! norm; the ring mean of the CTF squared and its azimuthal variation beside it; the log of the weighted power
    !! as the side scalar
    subroutine build_from_pft( self, pft, ctfmat, iseven, feat, normalise )
        class(polar_feats), intent(in)    :: self
        complex(sp),        intent(in)    :: pft(self%pftsz_in,self%kfromto(1):self%kfromto(2))
        real(sp),           intent(in)    :: ctfmat(self%pftsz_in,self%kfromto(1):self%kfromto(2))
        logical,            intent(in)    :: iseven
        type(polar_feat),   intent(inout) :: feat
        logical, optional,  intent(in)    :: normalise
        complex(sp) :: c(self%pftsz_in)
        real(sp)    :: xin(self%pftsz_in,self%kfromto(1):self%kfromto(2),2)
        real(dp)    :: wsq, scale, sk, t2(self%pftsz_in), tm, tv, taz, s
        integer     :: k, ih, nrings
        logical     :: l_norm
        if( .not. self%exists ) THROW_HARD('polar_feats not constructed; build_from_pft')
        l_norm = .true.
        if( present(normalise) ) l_norm = normalise
        call self%alloc_feat(feat)
        ih  = merge(1, 2, iseven)
        wsq = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            sk  = sum(real(pft(:,k) * conjg(pft(:,k)), dp))
            wsq = wsq + real(self%wk(k,ih),dp) * sk
        enddo
        scale = 1.d0
        ! the resampled ring holds pftsz / pftsz_in of the samples' power (the band's harmonics are kept)
        if( l_norm .and. wsq > 0.d0 ) scale = 1.d0 / sqrt(wsq * real(self%pftsz,dp) / real(self%pftsz_in,dp))
        taz    = 0.d0
        nrings = 0
        do k = self%kfromto(1),self%kfromto(2)
            c  = pft(:,k) * ctfmat(:,k)
            s  = real(self%sqw(k,ih),dp) * scale
            xin(:,k,1) = real(real(real(c),dp)  * s)
            xin(:,k,2) = real(real(aimag(c),dp) * s)
            t2 = real(ctfmat(:,k),dp)**2
            tm = sum(t2) / real(self%pftsz_in,dp)
            tv = max(0.d0, sum(t2*t2) / real(self%pftsz_in,dp) - tm*tm)
            feat%ctf2(k) = real(tm)
            if( tm > real(TINY,dp) )then
                taz    = taz + sqrt(tv) / tm
                nrings = nrings + 1
            endif
        enddo
        if( self%l_resample )then
            feat%x(:,:,1) = matmul(self%rs_re, xin(:,:,1))
            feat%x(:,:,2) = matmul(self%rs_im, xin(:,:,2))
        else
            feat%x = xin
        endif
        feat%iseven = iseven
        feat%wsq    = wsq
        feat%tazim  = 0.
        if( nrings > 0 ) feat%tazim = real(taz / real(nrings,dp))
        feat%logpow = real(log(max(wsq, real(TINY,dp))))
    end subroutine build_from_pft

    !> the calculator's half ring shifted by shift (pixels of the box) through its phase ramp exp(i (ax sx + ay sy)),
    !! as shift_ptcl applies it
    subroutine shift_pft( self, pft, shift, pft_sh )
        class(polar_feats), intent(in)  :: self
        complex(sp),        intent(in)  :: pft(self%pftsz_in,self%kfromto(1):self%kfromto(2))
        real(sp),           intent(in)  :: shift(2)
        complex(sp),        intent(out) :: pft_sh(self%pftsz_in,self%kfromto(1):self%kfromto(2))
        real(dp) :: arg
        integer  :: p, k
        if( .not. self%exists ) THROW_HARD('polar_feats not constructed; shift_pft')
        do k = self%kfromto(1),self%kfromto(2)
            do p = 1,self%pftsz_in
                arg = self%ax_in(p,k) * real(shift(1),dp) + self%ay_in(p,k) * real(shift(2),dp)
                pft_sh(p,k) = pft(p,k) * cmplx(cmplx(cos(arg), sin(arg), kind=dp), kind=sp)
            enddo
        enddo
    end subroutine shift_pft

    !> the packed tensor shifted by shift (pixels of the box): the phase ramp turns each (Re, Im) pair, the real
    !! scalings of the packing commuting with it; the augmentation's shift on a cached tensor
    subroutine shift_feat( self, x, shift, x_sh )
        class(polar_feats), intent(in)  :: self
        real(sp),           intent(in)  :: x(self%pftsz,self%kfromto(1):self%kfromto(2),2)
        real(sp),           intent(in)  :: shift(2)
        real(sp),           intent(out) :: x_sh(self%pftsz,self%kfromto(1):self%kfromto(2),2)
        real(dp) :: arg, c, s
        integer  :: p, k
        if( .not. self%exists ) THROW_HARD('polar_feats not constructed; shift_feat')
        do k = self%kfromto(1),self%kfromto(2)
            do p = 1,self%pftsz
                arg = self%ax(p,k) * real(shift(1),dp) + self%ay(p,k) * real(shift(2),dp)
                c   = cos(arg)
                s   = sin(arg)
                x_sh(p,k,1) = real(real(x(p,k,1),dp) * c - real(x(p,k,2),dp) * s)
                x_sh(p,k,2) = real(real(x(p,k,1),dp) * s + real(x(p,k,2),dp) * c)
            enddo
        enddo
    end subroutine shift_feat

    !> the packed tensor rotated by grid rotation irot (1: unchanged): entry p of the result is entry p - irot + 1
    !! of the Friedel-extended ring, whose second half is the conjugate of the first, so what passes the half ring's
    !! end returns with its imaginary channel negated (the conjugate wrap)
    subroutine rotate( self, x, irot, xrot )
        class(polar_feats), intent(in)  :: self
        real(sp),           intent(in)  :: x(self%pftsz,self%kfromto(1):self%kfromto(2),2)
        integer,            intent(in)  :: irot
        real(sp),           intent(out) :: xrot(self%pftsz,self%kfromto(1):self%kfromto(2),2)
        integer :: p, q
        do p = 1,self%pftsz
            q = modulo(p - irot, self%nrots) + 1
            if( q <= self%pftsz )then
                xrot(p,:,:) = x(q,:,:)
            else
                xrot(p,:,1) =  x(q-self%pftsz,:,1)
                xrot(p,:,2) = -x(q-self%pftsz,:,2)
            endif
        enddo
    end subroutine rotate

    !> the packed tensor padded along the angle by npad entries at both ends with the conjugate of the opposite end:
    !! the convolution's padding, which makes the half ring compute what the full ring would
    subroutine pad( self, x, npad, xpad )
        class(polar_feats), intent(in)  :: self
        real(sp),           intent(in)  :: x(self%pftsz,self%kfromto(1):self%kfromto(2),2)
        integer,            intent(in)  :: npad
        real(sp),           intent(out) :: xpad(1-npad:self%pftsz+npad,self%kfromto(1):self%kfromto(2),2)
        integer :: j
        if( npad < 0 .or. npad > self%pftsz ) THROW_HARD('npad must lie in [0, pftsz]; polar_feats%pad')
        xpad(1:self%pftsz,:,:) = x
        do j = 1,npad
            xpad(1-j,:,1)          =  x(self%pftsz+1-j,:,1)
            xpad(1-j,:,2)          = -x(self%pftsz+1-j,:,2)
            xpad(self%pftsz+j,:,1) =  x(j,:,1)
            xpad(self%pftsz+j,:,2) = -x(j,:,2)
        enddo
    end subroutine pad

    subroutine kill( self )
        class(polar_feats), intent(inout) :: self
        if( allocated(self%wk)  ) deallocate(self%wk)
        if( allocated(self%sqw) ) deallocate(self%sqw)
        if( allocated(self%ax)  ) deallocate(self%ax)
        if( allocated(self%ay)  ) deallocate(self%ay)
        if( allocated(self%ax_in) ) deallocate(self%ax_in)
        if( allocated(self%ay_in) ) deallocate(self%ay_in)
        if( allocated(self%rs_re) ) deallocate(self%rs_re)
        if( allocated(self%rs_im) ) deallocate(self%rs_im)
        self%l_resample = .false.
        self%pftsz_in = 0
        self%pftsz   = 0
        self%nrots   = 0
        self%box     = 0
        self%nk      = 0
        self%kfromto = 0
        self%exists  = .false.
    end subroutine kill

    subroutine feat_kill( self )
        class(polar_feat), intent(inout) :: self
        if( allocated(self%x)    ) deallocate(self%x)
        if( allocated(self%ctf2) ) deallocate(self%ctf2)
        self%tazim  = 0.
        self%logpow = 0.
        self%wsq    = 0.d0
        self%iseven = .true.
    end subroutine feat_kill

    ! ===== THE CONTRACT AND ITS HASHES =====

    !> one step of the two Lehmer hashes on a value: h1 < 2^31, A1 < 2^16, so no product leaves 2^47
    pure subroutine hash_step( h1, h2, v )
        integer(int64), intent(inout) :: h1, h2
        integer(int64), intent(in)    :: v
        h1 = modulo(h1 * HASH_A1 + modulo(v, HASH_M1), HASH_M1)
        h2 = modulo(h2 * HASH_A2 + modulo(v, HASH_M2), HASH_M2)
    end subroutine hash_step

    !> the hash of an integer list (order matters); the two generators' states packed into 62 bits
    pure function polar_feats_hash_ints( vals ) result( h )
        integer,        intent(in) :: vals(:)
        integer(int64) :: h, h1, h2
        integer :: i
        h1 = 1_int64
        h2 = 1_int64
        call hash_step(h1, h2, int(size(vals),int64))
        do i = 1,size(vals)
            call hash_step(h1, h2, int(vals(i),int64))
        enddo
        h = h1 + ishft(h2, 31)
    end function polar_feats_hash_ints

    !> the hash of a real list by its single-precision bit patterns (order matters)
    pure function polar_feats_hash_reals( vals ) result( h )
        real(sp),       intent(in) :: vals(:)
        integer(int64) :: h, h1, h2
        integer(int32) :: bits
        integer :: i
        h1 = 1_int64
        h2 = 1_int64
        call hash_step(h1, h2, int(size(vals),int64))
        do i = 1,size(vals)
            bits = transfer(vals(i), bits)
            call hash_step(h1, h2, int(bits,int64))
        enddo
        h = h1 + ishft(h2, 31)
    end function polar_feats_hash_reals

    !> the low 60 bits of the FNV-1a hash of a path
    function path_hash( path ) result( h )
        character(len=*), intent(in) :: path
        type(string)      :: str, hex
        character(len=16) :: digits
        integer(int64)    :: h
        str    = path
        hex    = str%to_fnv1a_hash64()
        digits = hex%to_char()
        read(digits(2:16), '(Z15)') h
        call str%kill
        call hex%kill
    end function path_hash

    !> The cache's contract: the project (its absolute path's hash and its modification time), the rows and the
    !! selection (count, hash, highest row), the geometry (box, cropped box, sampling, pftsz, the band), the shell
    !! weights' hash and whether the CTF is applied. A changed project, selection, geometry or spectrum changes
    !! the key and the cache is refused
    function polar_feats_contract( project_path, project_mtime, project_rows, box, box_crop, smpd, pinds, pftsz, kfromto, &
        &wk, with_ctf ) result( key )
        character(len=*), intent(in) :: project_path
        integer(int64),   intent(in) :: project_mtime
        integer,          intent(in) :: project_rows, box, box_crop, pinds(:), pftsz, kfromto(2)
        real,             intent(in) :: smpd
        real(sp),         intent(in) :: wk(:,:)
        logical,          intent(in) :: with_ctf
        integer(int64) :: key(POLAR_FEATS_NHEAD)
        key(IH_MAGIC)   = FEATS_CACHE_MAGIC
        key(IH_VERSION) = FEATS_CACHE_VERSION
        key(IH_PATH)    = path_hash(project_path)
        key(IH_MTIME)   = project_mtime
        key(IH_ROWS)    = int(project_rows,int64)
        key(IH_BOX)     = int(box,int64)
        key(IH_BOXCROP) = int(box_crop,int64)
        key(IH_SMPD)    = int(nint(smpd * 1.e6),int64)
        key(IH_NSEL)    = int(size(pinds),int64)
        key(IH_SELHASH) = polar_feats_hash_ints(pinds)
        key(IH_HIGHROW) = 0_int64
        if( size(pinds) > 0 ) key(IH_HIGHROW) = int(maxval(pinds),int64)
        key(IH_PFTSZ)   = int(pftsz,int64)
        key(IH_KFROM)   = int(kfromto(1),int64)
        key(IH_KTO)     = int(kfromto(2),int64)
        key(IH_WHASH)   = polar_feats_hash_reals(reshape(wk, [size(wk)]))
        key(IH_CTF)     = merge(1_int64, 0_int64, with_ctf)
    end function polar_feats_contract

    ! ===== THE CACHE =====

    subroutine cache_new( self, fname )
        class(polar_feats_cache), intent(inout) :: self
        class(string),            intent(in)    :: fname
        call self%kill
        self%fname   = fname
        self%tmpname = fname%to_char()//'_'//int2str(int(get_process_id()))//'_part'
    end subroutine cache_new

    !> the record length of one particle of the geometry, by inquiry
    subroutine cache_reclen( self, feats, feat )
        class(polar_feats_cache), intent(inout) :: self
        class(polar_feats),       intent(in)    :: feats
        type(polar_feat),         intent(inout) :: feat
        integer :: kfromto(2), reclen, ieven
        kfromto    = feats%get_kfromto()
        self%pftsz = feats%get_pftsz()
        self%kfrom = kfromto(1)
        self%nk    = kfromto(2) - kfromto(1) + 1
        call feats%alloc_feat(feat)
        ieven = 1
        inquire(iolength=reclen) feat%x, feat%ctf2, feat%tazim, feat%logpow, feat%wsq, ieven
        self%reclen = reclen
    end subroutine cache_reclen

    !> the in-progress file opened under the process token, after the space check, with a zero header: an
    !! unfinished cache never matches a contract; status 0, or -1 when the space is not there (nothing opened)
    subroutine cache_create( self, feats, highrow, status )
        class(polar_feats_cache), intent(inout) :: self
        class(polar_feats),       intent(in)    :: feats
        integer,                  intent(in)    :: highrow
        integer,                  intent(out)   :: status
        type(polar_feat) :: feat
        integer(int64)   :: key(POLAR_FEATS_NHEAD), want, avail
        integer          :: io_stat, head_len
        type(string)     :: dir
        status = -1
        if( self%l_open ) THROW_HARD('the cache is open; polar_feats_cache%create')
        if( highrow < 1 ) THROW_HARD('highrow must be positive; polar_feats_cache%create')
        call cache_reclen(self, feats, feat)
        inquire(iolength=head_len) key
        if( head_len > self%reclen ) THROW_HARD('the contract does not fit one record; polar_feats_cache%create')
        want  = int(highrow + 1,int64) * int(self%reclen,int64)
        dir   = simple_abspath('.')
        avail = fs_avail_bytes(dir)
        if( avail >= 0_int64 .and. want > avail / FREE_SPACE_DENOM )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> POLAR FEATS CACHE: ', want, ' bytes wanted, ', avail, &
                &' available; more than a quarter: no cache'
            call feat%kill
            return
        endif
        call del_file(self%tmpname)
        call fopen(self%funit, self%tmpname, 'replace', 'readwrite', iostat=io_stat, access='direct', &
            &form='unformatted', recl=self%reclen)
        call fileiochk('polar_feats_cache%create; open '//self%tmpname%to_char(), io_stat)
        key = 0_int64
        write(self%funit, rec=1) key
        self%highrow   = int(highrow,int64)
        self%l_open    = .true.
        self%l_writing = .true.
        status = 0
        call feat%kill
    end subroutine cache_create

    !> particle iptcl's feature into its record
    subroutine cache_put( self, iptcl, feat )
        class(polar_feats_cache), intent(inout) :: self
        integer,                  intent(in)    :: iptcl
        type(polar_feat),         intent(in)    :: feat
        integer :: ieven
        if( .not. (self%l_open .and. self%l_writing) ) THROW_HARD('the cache is not open for writing; polar_feats_cache%put')
        if( iptcl < 1 .or. int(iptcl,int64) > self%highrow ) THROW_HARD('iptcl outside the cache; polar_feats_cache%put')
        if( .not. allocated(feat%x) ) THROW_HARD('empty feature; polar_feats_cache%put')
        ieven = merge(1, 0, feat%iseven)
        write(self%funit, rec=iptcl+1) feat%x, feat%ctf2, feat%tazim, feat%logpow, feat%wsq, ieven
    end subroutine cache_put

    !> the publication: the contract last, then the rename into place
    subroutine cache_seal( self, key )
        class(polar_feats_cache), intent(inout) :: self
        integer(int64),           intent(in)    :: key(POLAR_FEATS_NHEAD)
        if( .not. (self%l_open .and. self%l_writing) ) THROW_HARD('the cache is not open for writing; polar_feats_cache%seal')
        if( key(IH_HIGHROW) > self%highrow ) THROW_HARD('the contract''s highest row exceeds the cache; polar_feats_cache%seal')
        write(self%funit, rec=1) key
        call fclose(self%funit)
        self%funit     = -1
        self%l_open    = .false.
        self%l_writing = .false.
        call simple_atomic_replace(self%tmpname, self%fname)
    end subroutine cache_seal

    !> whether a complete cache of this contract is in place: the file exists, its header equals the key and it
    !! holds every record up to the contract's highest row
    logical function cache_probe( self, feats, key ) result( l_avail )
        class(polar_feats_cache), intent(inout) :: self
        class(polar_feats),       intent(in)    :: feats
        integer(int64),           intent(in)    :: key(POLAR_FEATS_NHEAD)
        type(polar_feat) :: feat
        integer(int64)   :: stored(POLAR_FEATS_NHEAD), fsize
        integer          :: funit, io_stat
        l_avail = .false.
        if( self%l_open ) THROW_HARD('the cache is open; polar_feats_cache%probe')
        if( .not. file_exists(self%fname) ) return
        call cache_reclen(self, feats, feat)
        call feat%kill
        call fopen(funit, self%fname, 'old', 'read', iostat=io_stat, access='direct', form='unformatted', recl=self%reclen)
        if( io_stat /= 0 ) return
        read(funit, rec=1, iostat=io_stat) stored
        if( io_stat == 0 )then
            inquire(unit=funit, size=fsize)
            l_avail = all(stored == key) .and. fsize >= (key(IH_HIGHROW) + 1_int64) * int(self%reclen,int64)
        endif
        call fclose(funit)
    end function cache_probe

    !> the cache opened for reading after a successful probe
    subroutine cache_open_read( self, feats, key )
        class(polar_feats_cache), intent(inout) :: self
        class(polar_feats),       intent(in)    :: feats
        integer(int64),           intent(in)    :: key(POLAR_FEATS_NHEAD)
        integer :: io_stat
        if( .not. self%probe(feats, key) ) THROW_HARD('no cache of this contract: '//self%fname%to_char()//'; polar_feats_cache%open_read')
        call fopen(self%funit, self%fname, 'old', 'read', iostat=io_stat, access='direct', form='unformatted', recl=self%reclen)
        call fileiochk('polar_feats_cache%open_read; open '//self%fname%to_char(), io_stat)
        self%highrow   = key(IH_HIGHROW)
        self%l_open    = .true.
        self%l_writing = .false.
    end subroutine cache_open_read

    !> particle iptcl's feature from its record (the feature's arrays must match the cache's geometry)
    subroutine cache_get( self, iptcl, feat )
        class(polar_feats_cache), intent(in)    :: self
        integer,                  intent(in)    :: iptcl
        type(polar_feat),         intent(inout) :: feat
        integer :: ieven
        if( .not. (self%l_open .and. .not. self%l_writing) ) THROW_HARD('the cache is not open for reading; polar_feats_cache%get')
        if( iptcl < 1 .or. int(iptcl,int64) > self%highrow ) THROW_HARD('iptcl outside the cache; polar_feats_cache%get')
        if( .not. allocated(feat%x) ) THROW_HARD('allocate the feature first (alloc_feat); polar_feats_cache%get')
        if( size(feat%x,1) /= self%pftsz .or. size(feat%x,2) /= self%nk ) THROW_HARD('the feature does not match the cache''s geometry; polar_feats_cache%get')
        read(self%funit, rec=iptcl+1) feat%x, feat%ctf2, feat%tazim, feat%logpow, feat%wsq, ieven
        feat%iseven = ieven == 1
    end subroutine cache_get

    subroutine cache_close( self )
        class(polar_feats_cache), intent(inout) :: self
        if( .not. self%l_open ) return
        call fclose(self%funit)
        if( self%l_writing ) call del_file(self%tmpname)
        self%funit     = -1
        self%l_open    = .false.
        self%l_writing = .false.
    end subroutine cache_close

    !> the cache file removed (after training read its last epoch)
    subroutine cache_delete( self )
        class(polar_feats_cache), intent(inout) :: self
        call self%close
        if( file_exists(self%fname) ) call del_file(self%fname)
    end subroutine cache_delete

    subroutine cache_kill( self )
        class(polar_feats_cache), intent(inout) :: self
        call self%close
        call self%fname%kill
        call self%tmpname%kill
        self%pftsz   = 0
        self%nk      = 0
        self%kfrom   = 0
        self%reclen  = 0
        self%highrow = 0_int64
    end subroutine cache_kill

end module simple_polar_feats
