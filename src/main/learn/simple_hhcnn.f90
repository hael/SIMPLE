!@descr: the pose student's Hermitian network on the polar half ring: circular convolutions with real kernels and the conjugate wrap, the modulus gate, ring pooling, the rotation-searching template readout, its trainer and its weights file
! A model of the Fourier transform of a real image (doc/implementation_notes/completed/learned_pose_state_assignment.md,
! section 6): every layer commutes with conjugation, so the half ring with the conjugate wrap computes what the
! full ring would. The readout (ruling R3, 2026-11-07) keeps the last block's maps unpooled: a pointwise
! complex-linear map with real weights to CNN_NREAD channels, an average over ring bands, the band-limited angular
! spectrum of the full ring (even orders from the real part, odd from the imaginary: the Hermitian extension), and
! per direction class a learned template in that domain; the class score is the real Hermitian inner product
! maximised over the cyclic shifts of the ring (the rotation grid), times a learned scale; the softmax over the
! classes is the posterior and the gradient flows through the maximising shift. Forward and backward as im2col plus
! sgemm; a mini-batch is cut into fixed chunks computed in parallel into per-chunk gradient buffers and reduced
! serially in chunk order, so the same sum comes out at any thread count. The network is invariant to the mirror:
! its kernels are symmetric along the angle, so a reflected input gives reflected maps, and a class scores the map
! and its reflection, the larger; a class is a pair of mirror directions, which the search tells apart.
module simple_hhcnn
use, intrinsic :: iso_fortran_env, only: int64, int32
use simple_defs,   only: sp, dp, logfhandle, DTWOPI
use simple_error,  only: simple_exception
use simple_string, only: string
use simple_fileio, only: fopen, fclose, fileiochk
use simple_syslib, only: file_exists
!$ use omp_lib,    only: omp_get_max_threads, omp_get_thread_num, omp_in_parallel
implicit none

! the BLAS product (the library the build links; simple_linalg keeps its interface private)
interface
    subroutine sgemm( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
        character(len=1), intent(in)    :: transa, transb
        integer(kind=4),  intent(in)    :: m, n, k, lda, ldb, ldc
        real(kind=4),     intent(in)    :: alpha, beta
        real(kind=4),     intent(in)    :: a(lda,*), b(ldb,*)
        real(kind=4),     intent(inout) :: c(ldc,*)
    end subroutine sgemm
end interface

public :: hhcnn, run_all_hhcnn_tests
public :: CNN_NBLK, CNN_NREAD, CNN_KREAD, CNN_MREAD, CNN_CHUNK, CNN_WEIGHTS_VERSION
private
#include "simple_local_flags.inc"

integer,        parameter :: CNN_NBLK            = 3       !< convolution blocks
integer,        parameter :: CNN_NREAD           = 8       !< the readout map's channels
integer,        parameter :: CNN_KREAD           = 8       !< its ring bands (fewer when the last block has fewer rings)
integer,        parameter :: CNN_MREAD           = 16      !< the templates' band limit: full-ring angular orders 0..CNN_MREAD
integer,        parameter :: CNN_CHUNK           = 8       !< samples per chunk of a mini-batch (the reduction's unit)
integer(int64), parameter :: CNN_WEIGHTS_MAGIC   = 5417389302_int64
integer(int64), parameter :: CNN_WEIGHTS_VERSION = 2_int64
real(sp),       parameter :: SCALE_INIT_LOGIT    = 8.      !< the learned scale starts so that a class's own reprojection scores this logit on average
real(sp),       parameter :: GATE_DELTA2  = 1.e-12  !< the modulus is sqrt(re^2 + im^2 + delta^2): the gate is smooth at zero
real(sp),       parameter :: NORM_EPS     = 1.e-3   !< the first block's normaliser: the filter's energy under the CTF squared plus this share of its band mean
real(sp),       parameter :: ADAM_BETA1   = 0.9
real(sp),       parameter :: ADAM_BETA2   = 0.999
real(sp),       parameter :: ADAM_EPS     = 1.e-8
#ifdef SIMPLE_BLAS_TEAM_SAFE
logical, parameter :: L_BLAS_TEAM_SAFE = .true.
#else
logical, parameter :: L_BLAS_TEAM_SAFE = .false.
#endif

!> one convolution block: real kernels w(cout, 9 cin) over the 3 by 3 taps (ci fastest, then the angle offset,
!! then the ring offset), the modulus gate's bias b(cout), and their Adam moments
type :: conv_layer
    integer :: cin = 0, cout = 0, nkin = 0, nkout = 0, npos = 0
    real(sp), allocatable :: w(:,:), b(:)
    real(sp), allocatable :: mw(:,:), vw(:,:), mb(:), vb(:)
end type conv_layer

!> the template readout: the pointwise map wread(CNN_NREAD, cout) with real weights, the templates tmpl in the
!! angular harmonic domain as one array indexed by tix (part, channel, band, order 0..CNN_MREAD, class; the order-0
!! imaginary part is never used), the scale, and their Adam moments
type :: readout_layer
    integer :: cout = 0, kread = 0
    integer,  allocatable :: kb0(:), kb1(:)        !< (kread) the ring bands' first and last ring of the last block
    real(sp), allocatable :: wread(:,:), mwread(:,:), vwread(:,:)
    real(sp), allocatable :: tmpl(:), mtmpl(:), vtmpl(:)
    real(sp)              :: scale = 1., mscale = 0., vscale = 0.
end type readout_layer

!> one block's activations of one sample, forward and backward
type :: blk_buf
    real(sp), allocatable :: col(:,:)       !< (9 cin, 2 npos) the im2col matrix
    real(sp), allocatable :: raw(:,:,:,:)   !< (cout, pftsz, nkin, 2) the first block's convolution before its normalisation
    real(sp), allocatable :: dnorm(:,:)     !< (cout, nkin) the first block's normaliser per channel and ring
    real(sp), allocatable :: pre(:,:,:,:)   !< (cout, pftsz, nkin, 2) the convolution's output (normalised in block 1)
    real(sp), allocatable :: gated(:,:,:,:) !< after the modulus gate
    real(sp), allocatable :: pooled(:,:,:,:) !< (cout, pftsz, nkout, 2) after the ring pooling, the next block's input
    real(sp), allocatable :: dpre(:,:,:,:), dgated(:,:,:,:), dpooled(:,:,:,:) !< their gradients
    real(sp), allocatable :: dcol(:,:)      !< (9 cin, 2 npos)
end type blk_buf

!> one thread's workspace
type :: cnn_work
    type(blk_buf) :: blk(CNN_NBLK)
    real(sp), allocatable :: fmap(:,:,:,:)  !< (CNN_NREAD, pftsz, kread, 2) the readout map on the half ring
    real(sp), allocatable :: fhat(:,:,:,:)  !< (2, CNN_NREAD, kread, 0:CNN_MREAD) its full-ring angular spectrum (Re, Im)
    real(sp), allocatable :: dfmap(:,:,:,:), dfhat(:,:,:,:) !< their gradients
    real(sp), allocatable :: smax(:), logits(:), prob(:), dlogits(:) !< (nclasses) the maximised scores, the logits, the posterior
    integer,  allocatable :: rbest(:)       !< (nclasses) the maximising rotation of each class, 1..nrots
    logical,  allocatable :: refl(:)        !< (nclasses) whether the reflected map scored it
    real(sp), allocatable :: g(:)           !< (nk) the ring-averaged CTF squared of the sample (the shell weights are in the tensor)
    real(sp), allocatable :: wnorm(:)       !< (cout1) the first block's filter norms
    real(sp), allocatable :: cosm(:,:), sinm(:,:) !< (pftsz, 0:CNN_MREAD) the half-ring angles' full-ring harmonics
    real(sp), allocatable :: cosr(:,:), sinr(:,:) !< (nrots, 0:CNN_MREAD) the rotation grid's
end type cnn_work

!> one block's gradient sums
type :: blk_grad
    real(sp), allocatable :: dw(:,:), db(:)
end type blk_grad

!> one block's gradient sums in double, where a step reduces the chunks
type :: dw_set
    real(dp), allocatable :: dw(:,:), db(:)
end type dw_set

!> one chunk's gradient sums and loss
type :: cnn_grad
    type(blk_grad) :: blk(CNN_NBLK)
    real(sp), allocatable :: dwread(:,:), dtmpl(:)
    real(sp)              :: dscale = 0.
    real(dp)              :: loss  = 0.d0
    integer               :: nmode = 0      !< samples whose predicted mode is the target's
end type cnn_grad

type :: hhcnn
    private
    integer :: pftsz = 0, nrots = 0, nk = 0, nclasses = 0, width = 0, nfeat = 0, ntmpl = 0
    integer :: adam_step = 0
    integer(int64) :: seed = 0_int64
    type(conv_layer)    :: conv(CNN_NBLK)
    type(readout_layer) :: rd
    type(cnn_work), allocatable :: work(:)
    type(cnn_grad), allocatable :: gchunk(:)
    logical :: exists = .false.
  contains
    procedure :: new
    procedure :: get_nweights
    procedure :: get_shape
    procedure :: get_seed
    procedure :: prepare_work
    procedure :: predict
    procedure :: get_features
    procedure :: get_scores
    procedure :: init_templates
    procedure :: train_batch
    procedure :: eval_batch
    procedure :: write_weights
    procedure :: read_weights
    procedure :: copy_weights_from
    procedure :: weights_identical
    procedure :: kill
    procedure, private :: forward_sample
    procedure, private :: backward_sample
    procedure, private :: grad_chunk
    procedure, private :: zero_grad
    procedure, private :: reduce_and_step
    procedure, private :: alloc_work
    procedure, private :: alloc_grad
    procedure, private :: ensure_chunks
    procedure, private :: readout_forward
    procedure, private :: readout_backward
    procedure, private :: impose_symmetry
    procedure, private :: tie_kernel_grads
    procedure, private, nopass :: tix
    ! the layer kernels the tester drives on their own (bound so that the submodule links them)
    procedure, private, nopass :: im2col
    procedure, private, nopass :: sample_loss
    procedure, private, nopass :: grad_kill
end type hhcnn

interface
    module subroutine run_all_hhcnn_tests()
    end subroutine run_all_hhcnn_tests
end interface

contains

    ! ===== CONSTRUCTION =====

    !> the network of width channels in the first block (doubling per block) over pftsz angles and nk rings, with
    !! nclasses outputs (one template each), its weights drawn from seed (recorded in the model); the templates start
    !! at zero until init_templates sets them from the classes' reprojections
    subroutine new( self, pftsz, nk, nclasses, width, seed )
        class(hhcnn),   intent(inout) :: self
        integer,        intent(in)    :: pftsz, nk, nclasses, width
        integer(int64), intent(in)    :: seed
        integer(int64) :: state
        integer :: b, nkin, cin, kb, nper
        real(sp) :: scale
        call self%kill
        if( pftsz < 2 .or. nk < 1 .or. nclasses < 2 .or. width < 1 ) THROW_HARD('invalid network shape; hhcnn%new')
        self%pftsz    = pftsz
        self%nrots    = 2 * pftsz
        self%nk       = nk
        self%nclasses = nclasses
        self%width    = width
        self%seed     = seed
        state         = modulo(seed, 2147483646_int64) + 1_int64
        nkin = nk
        cin  = 1
        do b = 1,CNN_NBLK
            self%conv(b)%cin   = cin
            self%conv(b)%cout  = width * 2**(b-1)
            self%conv(b)%nkin  = nkin
            self%conv(b)%nkout = nkin / 2
            self%conv(b)%npos  = pftsz * nkin
            if( self%conv(b)%nkout < 1 ) THROW_HARD('too few rings for the ring pooling; hhcnn%new')
            allocate(self%conv(b)%w(self%conv(b)%cout, 9*cin), self%conv(b)%b(self%conv(b)%cout))
            allocate(self%conv(b)%mw(self%conv(b)%cout, 9*cin), self%conv(b)%vw(self%conv(b)%cout, 9*cin), &
                &self%conv(b)%mb(self%conv(b)%cout), self%conv(b)%vb(self%conv(b)%cout))
            ! He-type scale over the taps and input channels; the gate's bias starts at zero (a plain modulus)
            scale = sqrt(2. / real(9*cin))
            call fill_normal(state, self%conv(b)%w, scale)
            if( b == CNN_NBLK ) call self%impose_symmetry
            self%conv(b)%b  = 0.
            self%conv(b)%mw = 0.
            self%conv(b)%vw = 0.
            self%conv(b)%mb = 0.
            self%conv(b)%vb = 0.
            nkin = self%conv(b)%nkout
            cin  = self%conv(b)%cout
        enddo
        ! the readout: the pointwise map (a He-type scale over the last block's channels), the ring bands (contiguous
        ! groups of the last block's rings, as many as CNN_KREAD or the rings), the templates at zero, the scale at one
        self%rd%cout  = self%conv(CNN_NBLK)%cout
        self%rd%kread = min(CNN_KREAD, self%conv(CNN_NBLK)%nkout)
        allocate(self%rd%kb0(self%rd%kread), self%rd%kb1(self%rd%kread))
        nper = self%conv(CNN_NBLK)%nkout / self%rd%kread
        do kb = 1,self%rd%kread
            self%rd%kb0(kb) = (kb - 1) * nper + 1
            self%rd%kb1(kb) = kb * nper
        enddo
        self%rd%kb1(self%rd%kread) = self%conv(CNN_NBLK)%nkout
        allocate(self%rd%wread(CNN_NREAD, self%rd%cout), self%rd%mwread(CNN_NREAD, self%rd%cout), &
            &self%rd%vwread(CNN_NREAD, self%rd%cout))
        call fill_normal(state, self%rd%wread, sqrt(1. / real(self%rd%cout)))
        self%rd%mwread = 0.
        self%rd%vwread = 0.
        self%ntmpl = 2 * CNN_NREAD * CNN_KREAD * (CNN_MREAD + 1) * nclasses
        allocate(self%rd%tmpl(self%ntmpl), self%rd%mtmpl(self%ntmpl), self%rd%vtmpl(self%ntmpl), source=0.)
        self%rd%scale  = 1.
        self%rd%mscale = 0.
        self%rd%vscale = 0.
        self%nfeat     = 2 * CNN_NREAD * pftsz * self%rd%kread
        self%adam_step = 0
        self%exists    = .true.
    end subroutine new

    !> the templates' index: part (1 Re, 2 Im), channel n, band kb, order m (0..CNN_MREAD), class c
    pure integer function tix( n, kb, m, c, part )
        integer, intent(in) :: n, kb, m, c, part
        tix = part + 2 * ((n - 1) + CNN_NREAD * ((kb - 1) + CNN_KREAD * (m + (CNN_MREAD + 1) * (c - 1))))
    end function tix

    !> the next uniform variate in (0,1) of a Lehmer generator (modulus 2^31 - 1, no overflow in int64)
    real(dp) function lehmer_uniform( state )
        integer(int64), intent(inout) :: state
        integer(int64), parameter :: M = 2147483647_int64, A = 48271_int64
        state = modulo(state * A, M)
        if( state == 0_int64 ) state = 1_int64
        lehmer_uniform = real(state,dp) / real(M,dp)
    end function lehmer_uniform

    !> normal variates of standard deviation scale (Box-Muller on the generator)
    subroutine fill_normal( state, arr, scale )
        integer(int64), intent(inout) :: state
        real(sp),       intent(out)   :: arr(:,:)
        real(sp),       intent(in)    :: scale
        real(dp) :: u1, u2, r
        integer  :: i, j
        do j = 1,size(arr,2)
            do i = 1,size(arr,1)
                u1 = lehmer_uniform(state)
                u2 = lehmer_uniform(state)
                r  = sqrt(-2.d0 * log(u1))
                arr(i,j) = real(scale,sp) * real(r * cos(6.283185307179586d0 * u2),sp)
            enddo
        enddo
    end subroutine fill_normal

    pure integer function get_nweights( self )
        class(hhcnn), intent(in) :: self
        integer :: b
        get_nweights = 0
        do b = 1,CNN_NBLK
            get_nweights = get_nweights + size(self%conv(b)%w) + size(self%conv(b)%b)
        enddo
        get_nweights = get_nweights + size(self%rd%wread) + 2 * CNN_NREAD * self%rd%kread * (CNN_MREAD + 1) * self%nclasses + 1
    end function get_nweights

    !> the shape; nfeat the readout map's size (the features get_features returns)
    pure subroutine get_shape( self, pftsz, nk, nclasses, width, nfeat )
        class(hhcnn), intent(in)  :: self
        integer,      intent(out) :: pftsz, nk, nclasses, width, nfeat
        pftsz    = self%pftsz
        nk       = self%nk
        nclasses = self%nclasses
        width    = self%width
        nfeat    = self%nfeat
    end subroutine get_shape

    pure integer(int64) function get_seed( self )
        class(hhcnn), intent(in) :: self
        get_seed = self%seed
    end function get_seed

    ! ===== WORKSPACES =====

    !> one workspace per thread of the team that will run the regions (read omp_get_max_threads before them)
    subroutine prepare_work( self, nthr )
        class(hhcnn), intent(inout) :: self
        integer,      intent(in)    :: nthr
        integer :: ithr
        if( .not. self%exists ) THROW_HARD('network not constructed; hhcnn%prepare_work')
        if( allocated(self%work) )then
            if( size(self%work) >= nthr ) return
            do ithr = 1,size(self%work)
                call work_kill(self%work(ithr))
            enddo
            deallocate(self%work)
        endif
        allocate(self%work(nthr))
        do ithr = 1,nthr
            call self%alloc_work(self%work(ithr))
        enddo
    end subroutine prepare_work

    subroutine alloc_work( self, wk )
        class(hhcnn),   intent(in)    :: self
        type(cnn_work), intent(inout) :: wk
        integer :: b, cin, cout, nkin, nkout, npos, m, p
        do b = 1,CNN_NBLK
            cin   = self%conv(b)%cin
            cout  = self%conv(b)%cout
            nkin  = self%conv(b)%nkin
            nkout = self%conv(b)%nkout
            npos  = self%conv(b)%npos
            allocate(wk%blk(b)%col(9*cin, 2*npos), wk%blk(b)%dcol(9*cin, 2*npos), source=0.)
            allocate(wk%blk(b)%pre(cout, self%pftsz, nkin, 2), wk%blk(b)%gated(cout, self%pftsz, nkin, 2), source=0.)
            allocate(wk%blk(b)%pooled(cout, self%pftsz, nkout, 2), source=0.)
            allocate(wk%blk(b)%dpre(cout, self%pftsz, nkin, 2), wk%blk(b)%dgated(cout, self%pftsz, nkin, 2), source=0.)
            allocate(wk%blk(b)%dpooled(cout, self%pftsz, nkout, 2), source=0.)
            if( b == 1 )then
                allocate(wk%blk(b)%raw(cout, self%pftsz, nkin, 2), source=0.)
                allocate(wk%blk(b)%dnorm(cout, nkin), source=1.)
            endif
        enddo
        allocate(wk%fmap(CNN_NREAD, self%pftsz, self%rd%kread, 2), wk%dfmap(CNN_NREAD, self%pftsz, self%rd%kread, 2), source=0.)
        allocate(wk%fhat(2, CNN_NREAD, self%rd%kread, 0:CNN_MREAD), wk%dfhat(2, CNN_NREAD, self%rd%kread, 0:CNN_MREAD), source=0.)
        allocate(wk%smax(self%nclasses), wk%logits(self%nclasses), wk%prob(self%nclasses), wk%dlogits(self%nclasses), source=0.)
        allocate(wk%rbest(self%nclasses), source=1)
        allocate(wk%refl(self%nclasses),  source=.false.)
        allocate(wk%g(self%nk), source=0.)
        allocate(wk%wnorm(self%conv(1)%cout), source=1.)
        allocate(wk%cosm(self%pftsz, 0:CNN_MREAD), wk%sinm(self%pftsz, 0:CNN_MREAD), source=0.)
        allocate(wk%cosr(self%nrots, 0:CNN_MREAD), wk%sinr(self%nrots, 0:CNN_MREAD), source=0.)
        ! the full ring's harmonics (order m over nrots = 2 pftsz angles) at the half ring's angles and at the
        ! rotation grid, computed in double
        do m = 0,CNN_MREAD
            do p = 1,self%pftsz
                wk%cosm(p,m) = real(cos(DTWOPI * real(m,dp) * real(p-1,dp) / real(self%nrots,dp)))
                wk%sinm(p,m) = real(sin(DTWOPI * real(m,dp) * real(p-1,dp) / real(self%nrots,dp)))
            enddo
            do p = 1,self%nrots
                wk%cosr(p,m) = real(cos(DTWOPI * real(m,dp) * real(p-1,dp) / real(self%nrots,dp)))
                wk%sinr(p,m) = real(sin(DTWOPI * real(m,dp) * real(p-1,dp) / real(self%nrots,dp)))
            enddo
        enddo
    end subroutine alloc_work

    subroutine alloc_grad( self, gr )
        class(hhcnn),   intent(in)    :: self
        type(cnn_grad), intent(inout) :: gr
        integer :: b
        do b = 1,CNN_NBLK
            allocate(gr%blk(b)%dw(self%conv(b)%cout, 9*self%conv(b)%cin), gr%blk(b)%db(self%conv(b)%cout), source=0.)
        enddo
        allocate(gr%dwread(CNN_NREAD, self%rd%cout), source=0.)
        allocate(gr%dtmpl(self%ntmpl), source=0.)
        gr%dscale = 0.
        gr%loss   = 0.d0
        gr%nmode  = 0
    end subroutine alloc_grad

    subroutine zero_grad( self, gr )
        class(hhcnn),   intent(in)    :: self
        type(cnn_grad), intent(inout) :: gr
        integer :: b
        do b = 1,CNN_NBLK
            gr%blk(b)%dw = 0.
            gr%blk(b)%db = 0.
        enddo
        gr%dwread = 0.
        gr%dtmpl  = 0.
        gr%dscale = 0.
        gr%loss   = 0.d0
        gr%nmode  = 0
    end subroutine zero_grad

    !> the chunk buffers of a batch of n samples, and the team's workspaces
    subroutine ensure_chunks( self, n, nchunks )
        class(hhcnn), intent(inout) :: self
        integer,      intent(in)    :: n
        integer,      intent(out)   :: nchunks
        integer :: ic, nthr
        if( .not. self%exists ) THROW_HARD('network not constructed; hhcnn%ensure_chunks')
        if( n < 1 ) THROW_HARD('empty batch; hhcnn%ensure_chunks')
        nthr = 1
        !$ nthr = omp_get_max_threads()
        call self%prepare_work(nthr)
        nchunks = (n + CNN_CHUNK - 1) / CNN_CHUNK
        if( allocated(self%gchunk) )then
            if( size(self%gchunk) < nchunks )then
                do ic = 1,size(self%gchunk)
                    call grad_kill(self%gchunk(ic))
                enddo
                deallocate(self%gchunk)
            endif
        endif
        if( .not. allocated(self%gchunk) )then
            allocate(self%gchunk(nchunks))
            do ic = 1,nchunks
                call self%alloc_grad(self%gchunk(ic))
            enddo
        endif
    end subroutine ensure_chunks

    ! ===== THE LAYERS =====

    !> the BLAS may be called: outside a team, or inside one when the build knows it stays on the calling thread
    logical function blas_ok()
        logical :: l_team
        l_team = .false.
        !$ l_team = omp_in_parallel()
        blas_ok = L_BLAS_TEAM_SAFE .or. .not. l_team
    end function blas_ok

    !> c = beta c + a b (a: m x k, b: k x n), on the calling thread
    subroutine gemm_nn( m, n, k, a, lda, b, ldb, c, ldc, beta )
        integer,  intent(in)    :: m, n, k, lda, ldb, ldc
        real(sp), intent(in)    :: a(lda,*), b(ldb,*), beta
        real(sp), intent(inout) :: c(ldc,*)
        integer :: i, j, l
        real(sp) :: s
        if( blas_ok() )then
            call sgemm('N', 'N', m, n, k, 1.0, a, lda, b, ldb, beta, c, ldc)
        else
            do j = 1,n
                do i = 1,m
                    s = 0.
                    do l = 1,k
                        s = s + a(i,l) * b(l,j)
                    enddo
                    c(i,j) = beta * c(i,j) + s
                enddo
            enddo
        endif
    end subroutine gemm_nn

    !> c = c + a b^T (a: m x k, b: n x k)
    subroutine gemm_nt_acc( m, n, k, a, lda, b, ldb, c, ldc )
        integer,  intent(in)    :: m, n, k, lda, ldb, ldc
        real(sp), intent(in)    :: a(lda,*), b(ldb,*)
        real(sp), intent(inout) :: c(ldc,*)
        integer :: i, j, l
        real(sp) :: s
        if( blas_ok() )then
            call sgemm('N', 'T', m, n, k, 1.0, a, lda, b, ldb, 1.0, c, ldc)
        else
            do j = 1,n
                do i = 1,m
                    s = 0.
                    do l = 1,k
                        s = s + a(i,l) * b(j,l)
                    enddo
                    c(i,j) = c(i,j) + s
                enddo
            enddo
        endif
    end subroutine gemm_nt_acc

    !> c = a^T b (a: k x m, b: k x n)
    subroutine gemm_tn( m, n, k, a, lda, b, ldb, c, ldc )
        integer,  intent(in)    :: m, n, k, lda, ldb, ldc
        real(sp), intent(in)    :: a(lda,*), b(ldb,*)
        real(sp), intent(inout) :: c(ldc,*)
        integer :: i, j, l
        real(sp) :: s
        if( blas_ok() )then
            call sgemm('T', 'N', m, n, k, 1.0, a, lda, b, ldb, 0.0, c, ldc)
        else
            do j = 1,n
                do i = 1,m
                    s = 0.
                    do l = 1,k
                        s = s + a(l,i) * b(l,j)
                    enddo
                    c(i,j) = s
                enddo
            enddo
        endif
    end subroutine gemm_tn

    !> The im2col matrix of a(cin, nang, nring, 2) over the 3 by 3 taps: column q = p + nang (k-1) + npos (part-1)
    !! holds the taps of output position (p, k) of the real (part 1) or imaginary (part 2) channel; row
    !! r = ci + cin (ia+1) + 3 cin (ik+1) for the angle offset ia and the ring offset ik. Rings outside 1..nring
    !! are zero. Along the angle the ring wraps: with l_conj the half ring's conjugate wrap (what passes an end
    !! returns from the other end with its imaginary part negated), otherwise the plain cyclic wrap of a full ring
    subroutine im2col( cin, nang, nring, a, l_conj, col )
        integer,  intent(in)  :: cin, nang, nring
        real(sp), intent(in)  :: a(cin, nang, nring, *)
        logical,  intent(in)  :: l_conj
        real(sp), intent(out) :: col(9*cin, *)
        integer  :: part, k, p, ik, ia, kk, pp, q, r0, npos
        real(sp) :: sgn
        npos = nang * nring
        do part = 1,2
            do k = 1,nring
                do p = 1,nang
                    q = p + nang*(k-1) + npos*(part-1)
                    do ik = -1,1
                        kk = k + ik
                        if( kk < 1 .or. kk > nring )then
                            r0 = 3*cin*(ik+1)
                            col(r0+1:r0+3*cin, q) = 0.
                            cycle
                        endif
                        do ia = -1,1
                            pp  = p + ia
                            sgn = 1.
                            if( pp == 0 )then
                                pp = nang
                                if( l_conj .and. part == 2 ) sgn = -1.
                            else if( pp == nang + 1 )then
                                pp = 1
                                if( l_conj .and. part == 2 ) sgn = -1.
                            endif
                            r0 = cin*(ia+1) + 3*cin*(ik+1)
                            col(r0+1:r0+cin, q) = sgn * a(1:cin, pp, kk, part)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine im2col

    !> the transpose of im2col: the column gradients scattered back onto the input, the wrapped taps with their sign
    subroutine col2im( cin, nang, nring, dcol, l_conj, da )
        integer,  intent(in)  :: cin, nang, nring
        real(sp), intent(in)  :: dcol(9*cin, *)
        logical,  intent(in)  :: l_conj
        real(sp), intent(out) :: da(cin, nang, nring, 2)
        integer  :: part, k, p, ik, ia, kk, pp, q, r0, npos
        real(sp) :: sgn
        npos = nang * nring
        da   = 0.
        do part = 1,2
            do k = 1,nring
                do p = 1,nang
                    q = p + nang*(k-1) + npos*(part-1)
                    do ik = -1,1
                        kk = k + ik
                        if( kk < 1 .or. kk > nring ) cycle
                        do ia = -1,1
                            pp  = p + ia
                            sgn = 1.
                            if( pp == 0 )then
                                pp = nang
                                if( l_conj .and. part == 2 ) sgn = -1.
                            else if( pp == nang + 1 )then
                                pp = 1
                                if( l_conj .and. part == 2 ) sgn = -1.
                            endif
                            r0 = cin*(ia+1) + 3*cin*(ik+1)
                            da(1:cin, pp, kk, part) = da(1:cin, pp, kk, part) + sgn * dcol(r0+1:r0+cin, q)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine col2im

    !> the modulus gate z' = z relu(|z| + b) / |z|, in place; |z| is smoothed at zero
    subroutine gate_forward( cout, nang, nring, b, z )
        integer,  intent(in)    :: cout, nang, nring
        real(sp), intent(in)    :: b(cout)
        real(sp), intent(inout) :: z(cout, nang, nring, 2)
        real(sp) :: m, f
        integer  :: c, p, k
        do k = 1,nring
            do p = 1,nang
                do c = 1,cout
                    m = sqrt(z(c,p,k,1)**2 + z(c,p,k,2)**2 + GATE_DELTA2)
                    if( m + b(c) > 0. )then
                        f = (m + b(c)) / m
                        z(c,p,k,1) = z(c,p,k,1) * f
                        z(c,p,k,2) = z(c,p,k,2) * f
                    else
                        z(c,p,k,1) = 0.
                        z(c,p,k,2) = 0.
                    endif
                enddo
            enddo
        enddo
    end subroutine gate_forward

    !> the gate's backward from its input z and its output gradient dout: dz, and db accumulated
    subroutine gate_backward( cout, nang, nring, b, z, dout, dz, db )
        integer,  intent(in)    :: cout, nang, nring
        real(sp), intent(in)    :: b(cout), z(cout, nang, nring, 2), dout(cout, nang, nring, 2)
        real(sp), intent(out)   :: dz(cout, nang, nring, 2)
        real(sp), intent(inout) :: db(cout)
        real(sp) :: m, f, zg
        integer  :: c, p, k
        do k = 1,nring
            do p = 1,nang
                do c = 1,cout
                    m = sqrt(z(c,p,k,1)**2 + z(c,p,k,2)**2 + GATE_DELTA2)
                    if( m + b(c) > 0. )then
                        f  = (m + b(c)) / m
                        zg = z(c,p,k,1) * dout(c,p,k,1) + z(c,p,k,2) * dout(c,p,k,2)
                        dz(c,p,k,1) = f * dout(c,p,k,1) - b(c) * zg * z(c,p,k,1) / m**3
                        dz(c,p,k,2) = f * dout(c,p,k,2) - b(c) * zg * z(c,p,k,2) / m**3
                        db(c) = db(c) + zg / m
                    else
                        dz(c,p,k,1) = 0.
                        dz(c,p,k,2) = 0.
                    endif
                enddo
            enddo
        enddo
    end subroutine gate_backward

    !> average pooling of ring pairs (an odd last ring is dropped)
    subroutine pool_forward( cout, nang, nkin, nkout, a, pooled )
        integer,  intent(in)  :: cout, nang, nkin, nkout
        real(sp), intent(in)  :: a(cout, nang, nkin, 2)
        real(sp), intent(out) :: pooled(cout, nang, nkout, 2)
        integer :: kk
        do kk = 1,nkout
            pooled(:,:,kk,:) = 0.5 * (a(:,:,2*kk-1,:) + a(:,:,2*kk,:))
        enddo
    end subroutine pool_forward

    subroutine pool_backward( cout, nang, nkin, nkout, dpooled, da )
        integer,  intent(in)  :: cout, nang, nkin, nkout
        real(sp), intent(in)  :: dpooled(cout, nang, nkout, 2)
        real(sp), intent(out) :: da(cout, nang, nkin, 2)
        integer :: kk
        da = 0.
        do kk = 1,nkout
            da(:,:,2*kk-1,:) = 0.5 * dpooled(:,:,kk,:)
            da(:,:,2*kk,:)   = 0.5 * dpooled(:,:,kk,:)
        enddo
    end subroutine pool_backward

    !> the first block's normaliser: the filter's energy under the ring-averaged CTF squared (g over the rings, zero
    !! outside the band; the Gram's diagonal of the packed domain, whose tensor carries the shell weights) plus
    !! NORM_EPS times its band mean
    subroutine norm_factors( cout, nk, w, g, d )
        integer,  intent(in)  :: cout, nk
        real(sp), intent(in)  :: w(cout, 9), g(nk)
        real(sp), intent(out) :: d(cout, nk)
        real(sp) :: wsq(cout, -1:1), gbar(-1:1), e
        integer  :: c, k, ik, ia, kk
        do ik = -1,1
            do c = 1,cout
                wsq(c,ik) = 0.
                do ia = -1,1
                    wsq(c,ik) = wsq(c,ik) + w(c, 1 + (ia+1) + 3*(ik+1))**2
                enddo
            enddo
            gbar(ik) = 0.
            do k = 1,nk
                kk = k + ik
                if( kk >= 1 .and. kk <= nk ) gbar(ik) = gbar(ik) + g(kk)
            enddo
            gbar(ik) = gbar(ik) / real(nk)
        enddo
        do k = 1,nk
            do c = 1,cout
                e = 0.
                do ik = -1,1
                    kk = k + ik
                    if( kk >= 1 .and. kk <= nk ) e = e + wsq(c,ik) * g(kk)
                    e = e + NORM_EPS * wsq(c,ik) * gbar(ik)
                enddo
                if( e > 0. )then
                    d(c,k) = e
                else
                    d(c,k) = 1.
                endif
            enddo
        enddo
    end subroutine norm_factors

    !> The forward pass of one sample in workspace wk: x(pftsz, nk, 2) the packed half ring, g(nk) the ring-averaged
    !! CTF squared, side the side scalar (not read by the template readout; kept in the interface); the class
    !! probabilities land in wk%prob
    subroutine forward_sample( self, wk, x, g, side )
        class(hhcnn),   intent(in)    :: self
        type(cnn_work), intent(inout) :: wk
        real(sp),       intent(in)    :: x(self%pftsz, self%nk, 2), g(self%nk), side
        integer  :: b, cin, cout, nkin, nkout, npos, c, k, p
        wk%g = g
        do b = 1,CNN_NBLK
            cin   = self%conv(b)%cin
            cout  = self%conv(b)%cout
            nkin  = self%conv(b)%nkin
            nkout = self%conv(b)%nkout
            npos  = self%conv(b)%npos
            if( b == 1 )then
                call im2col(1, self%pftsz, nkin, x, .true., wk%blk(1)%col)
            else
                call im2col(cin, self%pftsz, nkin, wk%blk(b-1)%pooled, .true., wk%blk(b)%col)
            endif
            call gemm_nn(cout, 2*npos, 9*cin, self%conv(b)%w, cout, wk%blk(b)%col, 9*cin, wk%blk(b)%pre, cout, 0.0)
            if( b == 1 )then
                ! the normalisation of the unit-norm filter: the response times the filter's norm, divided by its
                ! energy under the weights and the CTF squared; the raw response kept for the backward pass
                wk%blk(1)%raw = wk%blk(1)%pre
                call norm_factors(cout, self%nk, self%conv(1)%w, wk%g, wk%blk(1)%dnorm)
                do c = 1,cout
                    wk%wnorm(c) = sqrt(sum(self%conv(1)%w(c,:)**2))
                enddo
                do k = 1,nkin
                    do p = 1,self%pftsz
                        do c = 1,cout
                            wk%blk(1)%pre(c,p,k,1) = wk%blk(1)%raw(c,p,k,1) * wk%wnorm(c) / wk%blk(1)%dnorm(c,k)
                            wk%blk(1)%pre(c,p,k,2) = wk%blk(1)%raw(c,p,k,2) * wk%wnorm(c) / wk%blk(1)%dnorm(c,k)
                        enddo
                    enddo
                enddo
            endif
            wk%blk(b)%gated = wk%blk(b)%pre
            call gate_forward(cout, self%pftsz, nkin, self%conv(b)%b, wk%blk(b)%gated)
            call pool_forward(cout, self%pftsz, nkin, nkout, wk%blk(b)%gated, wk%blk(b)%pooled)
        enddo
        call self%readout_forward(wk)
    end subroutine forward_sample

    !> The template readout of the last block's pooled maps: the readout map (the pointwise map with real weights,
    !! the mean over each ring band), its full-ring angular spectrum to order CNN_MREAD (even orders from the real
    !! part, odd orders from the imaginary part: the Hermitian extension's), and per class the real inner product
    !! with the class's template at every rotation of the grid, S(r) = (1/nrots) [Re P(0) + 2 sum_m Re(P(m) w^{m r})]
    !! with P(m) the spectra's product summed over channels and bands, maximised over r and over the map and its
    !! reflection (the conjugate spectrum, conj F); the logits are the scale times the maxima, the posterior their softmax
    subroutine readout_forward( self, wk )
        class(hhcnn),   intent(in)    :: self
        type(cnn_work), intent(inout) :: wk
        integer  :: cout, kb, part, p, n, kk, c, m, ic, r, nb, ir
        real(sp) :: pr(0:CNN_MREAD), pi(0:CNN_MREAD), sc, smx, zmax, inv_n, tr, ti, fr, fi, sgn
        real(dp) :: acc, sr, si, ssum
        cout  = self%rd%cout
        inv_n = 1. / real(self%nrots)
        ! the readout map
        do part = 1,2
            do kb = 1,self%rd%kread
                nb = self%rd%kb1(kb) - self%rd%kb0(kb) + 1
                do p = 1,self%pftsz
                    do n = 1,CNN_NREAD
                        acc = 0.d0
                        do kk = self%rd%kb0(kb),self%rd%kb1(kb)
                            do c = 1,cout
                                acc = acc + real(self%rd%wread(n,c) * wk%blk(CNN_NBLK)%pooled(c,p,kk,part),dp)
                            enddo
                        enddo
                        wk%fmap(n,p,kb,part) = real(acc / real(nb,dp))
                    enddo
                enddo
            enddo
        enddo
        ! its spectrum over the full ring: F(m) = sum_p [z_p + (-1)^m conj z_p] e^{-i theta_{m,p}}
        do m = 0,CNN_MREAD
            do kb = 1,self%rd%kread
                do n = 1,CNN_NREAD
                    sr = 0.d0
                    si = 0.d0
                    if( mod(m,2) == 0 )then
                        do p = 1,self%pftsz
                            sr = sr + real(wk%fmap(n,p,kb,1) * wk%cosm(p,m),dp)
                            si = si - real(wk%fmap(n,p,kb,1) * wk%sinm(p,m),dp)
                        enddo
                    else
                        do p = 1,self%pftsz
                            sr = sr + real(wk%fmap(n,p,kb,2) * wk%sinm(p,m),dp)
                            si = si + real(wk%fmap(n,p,kb,2) * wk%cosm(p,m),dp)
                        enddo
                    endif
                    wk%fhat(1,n,kb,m) = real(2.d0 * sr)
                    wk%fhat(2,n,kb,m) = real(2.d0 * si)
                enddo
            enddo
        enddo
        ! the classes: the spectra's product with the map and with its reflection (the conjugate spectrum), the rotation
        ! search of both, the maximum
        do ic = 1,self%nclasses
            smx = -huge(smx)
            wk%rbest(ic) = 1
            wk%refl(ic)  = .false.
            do ir = 0,1
                sgn = merge(1., -1., ir == 0)
                do m = 0,CNN_MREAD
                    sr = 0.d0
                    si = 0.d0
                    do kb = 1,self%rd%kread
                        do n = 1,CNN_NREAD
                            tr = self%rd%tmpl(tix(n,kb,m,ic,1))
                            ti = self%rd%tmpl(tix(n,kb,m,ic,2))
                            fr = wk%fhat(1,n,kb,m)
                            fi = sgn * wk%fhat(2,n,kb,m)
                            sr = sr + real(fr * tr + fi * ti,dp)
                            si = si + real(fi * tr - fr * ti,dp)
                        enddo
                    enddo
                    pr(m) = real(sr)
                    pi(m) = real(si)
                enddo
                do r = 1,self%nrots
                    ssum = real(pr(0),dp)
                    do m = 1,CNN_MREAD
                        ssum = ssum + 2.d0 * real(pr(m) * wk%cosr(r,m) - pi(m) * wk%sinr(r,m),dp)
                    enddo
                    sc = real(ssum) * inv_n
                    if( sc > smx )then
                        smx = sc
                        wk%rbest(ic) = r
                        wk%refl(ic)  = ir == 1
                    endif
                enddo
            enddo
            wk%smax(ic)   = smx
            wk%logits(ic) = self%rd%scale * smx
        enddo
        zmax    = maxval(wk%logits)
        wk%prob = exp(wk%logits - zmax)
        ssum    = sum(real(wk%prob,dp))
        wk%prob = real(real(wk%prob,dp) / ssum)
    end subroutine readout_forward

    !> the readout's backward from the logits' gradient wk%dlogits: the scale's, the templates' and the map's
    !! gradients into gr, the last block's dpooled; the gradient flows through each class's maximising rotation and,
    !! when the reflection scored it, through the conjugate spectrum
    subroutine readout_backward( self, wk, gr )
        class(hhcnn),   intent(in)    :: self
        type(cnn_work), intent(inout) :: wk
        type(cnn_grad), intent(inout) :: gr
        integer  :: cout, kb, part, p, n, kk, c, m, ic, r, nb
        real(sp) :: ds, dpr, dpi, inv_n, tr, ti, fr, fi, cr, sr, sgn
        real(dp) :: acc, dsc
        cout  = self%rd%cout
        inv_n = 1. / real(self%nrots)
        dsc   = 0.d0
        wk%dfhat = 0.
        do ic = 1,self%nclasses
            dsc = dsc + real(wk%dlogits(ic) * wk%smax(ic),dp)
            ds  = wk%dlogits(ic) * self%rd%scale
            r   = wk%rbest(ic)
            sgn = merge(-1., 1., wk%refl(ic))
            do m = 0,CNN_MREAD
                if( m == 0 )then
                    dpr = ds * inv_n
                    dpi = 0.
                else
                    cr  = wk%cosr(r,m)
                    sr  = wk%sinr(r,m)
                    dpr =  ds * 2. * inv_n * cr
                    dpi = -ds * 2. * inv_n * sr
                endif
                do kb = 1,self%rd%kread
                    do n = 1,CNN_NREAD
                        tr = self%rd%tmpl(tix(n,kb,m,ic,1))
                        ti = self%rd%tmpl(tix(n,kb,m,ic,2))
                        fr = wk%fhat(1,n,kb,m)
                        fi = sgn * wk%fhat(2,n,kb,m)
                        gr%dtmpl(tix(n,kb,m,ic,1)) = gr%dtmpl(tix(n,kb,m,ic,1)) + dpr * fr + dpi * fi
                        if( m > 0 ) gr%dtmpl(tix(n,kb,m,ic,2)) = gr%dtmpl(tix(n,kb,m,ic,2)) + dpr * fi - dpi * fr
                        wk%dfhat(1,n,kb,m) = wk%dfhat(1,n,kb,m) + dpr * tr - dpi * ti
                        wk%dfhat(2,n,kb,m) = wk%dfhat(2,n,kb,m) + sgn * (dpr * ti + dpi * tr)
                    enddo
                enddo
            enddo
        enddo
        gr%dscale = gr%dscale + real(dsc)
        ! the spectrum's backward onto the half-ring map
        wk%dfmap = 0.
        do m = 0,CNN_MREAD
            do kb = 1,self%rd%kread
                do n = 1,CNN_NREAD
                    fr = 2. * wk%dfhat(1,n,kb,m)
                    fi = 2. * wk%dfhat(2,n,kb,m)
                    if( mod(m,2) == 0 )then
                        do p = 1,self%pftsz
                            wk%dfmap(n,p,kb,1) = wk%dfmap(n,p,kb,1) + fr * wk%cosm(p,m) - fi * wk%sinm(p,m)
                        enddo
                    else
                        do p = 1,self%pftsz
                            wk%dfmap(n,p,kb,2) = wk%dfmap(n,p,kb,2) + fr * wk%sinm(p,m) + fi * wk%cosm(p,m)
                        enddo
                    endif
                enddo
            enddo
        enddo
        ! the map's backward: the pointwise weights and the pooled maps
        wk%blk(CNN_NBLK)%dpooled = 0.
        do part = 1,2
            do kb = 1,self%rd%kread
                nb = self%rd%kb1(kb) - self%rd%kb0(kb) + 1
                do p = 1,self%pftsz
                    do kk = self%rd%kb0(kb),self%rd%kb1(kb)
                        do c = 1,cout
                            acc = 0.d0
                            do n = 1,CNN_NREAD
                                acc = acc + real(self%rd%wread(n,c) * wk%dfmap(n,p,kb,part),dp)
                            enddo
                            wk%blk(CNN_NBLK)%dpooled(c,p,kk,part) = wk%blk(CNN_NBLK)%dpooled(c,p,kk,part) + real(acc / real(nb,dp))
                        enddo
                    enddo
                    do c = 1,cout
                        acc = 0.d0
                        do kk = self%rd%kb0(kb),self%rd%kb1(kb)
                            acc = acc + real(wk%blk(CNN_NBLK)%pooled(c,p,kk,part),dp)
                        enddo
                        acc = acc / real(nb,dp)
                        do n = 1,CNN_NREAD
                            gr%dwread(n,c) = gr%dwread(n,c) + real(acc * real(wk%dfmap(n,p,kb,part),dp))
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine readout_backward

    !> The templates from geometry (ruling R3): each class's template is the readout spectrum of that class's
    !! reprojection xrefs(:,:,:,class) under g through the network as it stands, and the scale is set so that a
    !! class's own reprojection scores SCALE_INIT_LOGIT on average; on thread 1's workspace (prepare_work first)
    subroutine init_templates( self, xrefs, g )
        class(hhcnn), intent(inout) :: self
        real(sp),     intent(in)    :: xrefs(self%pftsz, self%nk, 2, self%nclasses), g(self%nk)
        integer  :: ic, m, kb, n
        real(dp) :: sself, ssum
        if( .not. allocated(self%work) ) THROW_HARD('no workspace: prepare_work first; hhcnn%init_templates')
        self%rd%tmpl  = 0.
        self%rd%scale = 1.
        ssum = 0.d0
        do ic = 1,self%nclasses
            call self%forward_sample(self%work(1), xrefs(:,:,:,ic), g, 0.)
            sself = 0.d0
            do m = 0,CNN_MREAD
                do kb = 1,self%rd%kread
                    do n = 1,CNN_NREAD
                        self%rd%tmpl(tix(n,kb,m,ic,1)) = self%work(1)%fhat(1,n,kb,m)
                        if( m > 0 ) self%rd%tmpl(tix(n,kb,m,ic,2)) = self%work(1)%fhat(2,n,kb,m)
                        if( m == 0 )then
                            sself = sself + real(self%work(1)%fhat(1,n,kb,m),dp)**2
                        else
                            sself = sself + 2.d0 * (real(self%work(1)%fhat(1,n,kb,m),dp)**2 + real(self%work(1)%fhat(2,n,kb,m),dp)**2)
                        endif
                    enddo
                enddo
            enddo
            ssum = ssum + sself / real(self%nrots,dp)
        enddo
        ssum = ssum / real(self%nclasses,dp)
        if( ssum > 0.d0 ) self%rd%scale = real(real(SCALE_INIT_LOGIT,dp) / ssum)
        self%rd%mtmpl  = 0.
        self%rd%vtmpl  = 0.
        self%rd%mscale = 0.
        self%rd%vscale = 0.
    end subroutine init_templates

    !> the kernels symmetric under the reflection of the angle, their two angle offsets averaged (idempotent)
    subroutine impose_symmetry( self )
        class(hhcnn), intent(inout) :: self
        integer :: b, cin, ci, ik, rlo, rhi
        do b = 1,CNN_NBLK
            cin = self%conv(b)%cin
            do ik = -1,1
                do ci = 1,cin
                    rlo = ci + 3*cin*(ik+1)
                    rhi = rlo + 2*cin
                    self%conv(b)%w(:,rlo) = 0.5 * (self%conv(b)%w(:,rlo) + self%conv(b)%w(:,rhi))
                    self%conv(b)%w(:,rhi) = self%conv(b)%w(:,rlo)
                enddo
            enddo
        enddo
    end subroutine impose_symmetry

    !> the gradient of a symmetric kernel's shared weight onto both offsets (the chain rule), so that Adam, acting per
    !! element, keeps the symmetry exact
    subroutine tie_kernel_grads( self, dw )
        class(hhcnn), intent(in)    :: self
        type(dw_set), intent(inout) :: dw(CNN_NBLK)
        integer :: b, cin, ci, ik, rlo, rhi
        do b = 1,CNN_NBLK
            cin = self%conv(b)%cin
            do ik = -1,1
                do ci = 1,cin
                    rlo = ci + 3*cin*(ik+1)
                    rhi = rlo + 2*cin
                    dw(b)%dw(:,rlo) = dw(b)%dw(:,rlo) + dw(b)%dw(:,rhi)
                    dw(b)%dw(:,rhi) = dw(b)%dw(:,rlo)
                enddo
            enddo
        enddo
    end subroutine tie_kernel_grads

    !> the backward pass of one sample after forward_sample, for the soft targets t: the gradients accumulated
    !! into gr; the loss (the cross-entropy against t) returned
    subroutine backward_sample( self, wk, t, gr, loss )
        class(hhcnn),   intent(in)    :: self
        type(cnn_work), intent(inout) :: wk
        real(sp),       intent(in)    :: t(self%nclasses)
        type(cnn_grad), intent(inout) :: gr
        real(dp),       intent(out)   :: loss
        real(sp) :: s_norm(self%conv(1)%cout, self%nk), ysum(self%conv(1)%cout), wsum, gterm, gbar(-1:1)
        integer  :: b, cin, cout, nkin, nkout, npos, c, k, p, kk, ik, ia, r
        ! the loss and the logits' gradient, then the readout's backward
        call sample_loss(self%nclasses, wk%prob, t, loss)
        wk%dlogits = wk%prob - t
        if( maxloc(wk%prob, dim=1) == maxloc(t, dim=1) ) gr%nmode = gr%nmode + 1
        call self%readout_backward(wk, gr)
        ! the blocks, last first
        do b = CNN_NBLK,1,-1
            cin   = self%conv(b)%cin
            cout  = self%conv(b)%cout
            nkin  = self%conv(b)%nkin
            nkout = self%conv(b)%nkout
            npos  = self%conv(b)%npos
            call pool_backward(cout, self%pftsz, nkin, nkout, wk%blk(b)%dpooled, wk%blk(b)%dgated)
            call gate_backward(cout, self%pftsz, nkin, self%conv(b)%b, wk%blk(b)%pre, wk%blk(b)%dgated, wk%blk(b)%dpre, gr%blk(b)%db)
            if( b == 1 )then
                ! through the normalisation: the response's gradient, the normaliser's and the filter norm's
                ! dependence on the weights (y = raw |w| / d: dy/dd = -y/d, dy/d|w| = y/|w|)
                ysum = 0.
                do k = 1,nkin
                    do c = 1,cout
                        s_norm(c,k) = 0.
                        do p = 1,self%pftsz
                            s_norm(c,k) = s_norm(c,k) - (wk%blk(1)%dpre(c,p,k,1) * wk%blk(1)%pre(c,p,k,1) + &
                                &wk%blk(1)%dpre(c,p,k,2) * wk%blk(1)%pre(c,p,k,2))
                            wk%blk(1)%dpre(c,p,k,1) = wk%blk(1)%dpre(c,p,k,1) * wk%wnorm(c) / wk%blk(1)%dnorm(c,k)
                            wk%blk(1)%dpre(c,p,k,2) = wk%blk(1)%dpre(c,p,k,2) * wk%wnorm(c) / wk%blk(1)%dnorm(c,k)
                        enddo
                        ysum(c)     = ysum(c) - s_norm(c,k)
                        s_norm(c,k) = s_norm(c,k) / wk%blk(1)%dnorm(c,k)
                    enddo
                enddo
                do r = 1,9
                    do c = 1,cout
                        gr%blk(1)%dw(c,r) = gr%blk(1)%dw(c,r) + self%conv(1)%w(c,r) * ysum(c) / wk%wnorm(c)**2
                    enddo
                enddo
                do ik = -1,1
                    gbar(ik) = 0.
                    do k = 1,self%nk
                        kk = k + ik
                        if( kk >= 1 .and. kk <= self%nk ) gbar(ik) = gbar(ik) + wk%g(kk)
                    enddo
                    gbar(ik) = gbar(ik) / real(self%nk)
                enddo
                do ik = -1,1
                    do ia = -1,1
                        r = 1 + (ia+1) + 3*(ik+1)
                        do c = 1,cout
                            wsum = 0.
                            do k = 1,self%nk
                                kk    = k + ik
                                gterm = NORM_EPS * gbar(ik)
                                if( kk >= 1 .and. kk <= self%nk ) gterm = gterm + wk%g(kk)
                                wsum = wsum + s_norm(c,k) * gterm
                            enddo
                            gr%blk(1)%dw(c,r) = gr%blk(1)%dw(c,r) + 2. * self%conv(1)%w(c,r) * wsum
                        enddo
                    enddo
                enddo
            endif
            ! the convolution's weights and, except for the first block, its input
            call gemm_nt_acc(cout, 9*cin, 2*npos, wk%blk(b)%dpre, cout, wk%blk(b)%col, 9*cin, gr%blk(b)%dw, cout)
            if( b > 1 )then
                call gemm_tn(9*cin, 2*npos, cout, self%conv(b)%w, cout, wk%blk(b)%dpre, cout, wk%blk(b)%dcol, 9*cin)
                call col2im(cin, self%pftsz, nkin, wk%blk(b)%dcol, .true., wk%blk(b-1)%dpooled)
            endif
        enddo
    end subroutine backward_sample

    pure subroutine sample_loss( nclasses, prob, t, loss )
        integer,  intent(in)  :: nclasses
        real(sp), intent(in)  :: prob(nclasses), t(nclasses)
        real(dp), intent(out) :: loss
        integer :: i
        loss = 0.d0
        do i = 1,nclasses
            if( t(i) > 0. ) loss = loss - real(t(i),dp) * log(max(real(prob(i),dp), 1.d-30))
        enddo
    end subroutine sample_loss

    !> the gradients of samples i0..i1 of the batch into chunk buffer ic, on thread ithr's workspace
    subroutine grad_chunk( self, ic, ithr, i0, i1, x, g, side, t )
        class(hhcnn), intent(inout) :: self
        integer,      intent(in)    :: ic, ithr, i0, i1
        real(sp),     intent(in)    :: x(self%pftsz, self%nk, 2, *), g(self%nk, *), side(*), t(self%nclasses, *)
        real(dp) :: loss
        integer  :: i
        call self%zero_grad(self%gchunk(ic))
        do i = i0,i1
            call self%forward_sample(self%work(ithr), x(:,:,:,i), g(:,i), side(i))
            call self%backward_sample(self%work(ithr), t(:,i), self%gchunk(ic), loss)
            self%gchunk(ic)%loss = self%gchunk(ic)%loss + loss
        enddo
    end subroutine grad_chunk

    !> One training step on a mini-batch of n samples: x(pftsz, nk, 2, n) the packed rings, g(nk, n) the ring-averaged
    !! CTF squared, side(n), t(nclasses, n) the soft targets. The batch is cut into
    !! chunks of CNN_CHUNK computed in parallel into their own buffers and reduced serially in chunk order, then
    !! one Adam step at rate lr with decoupled weight decay wd; the mean loss and the mode agreement returned
    subroutine train_batch( self, n, x, g, side, t, lr, wd, loss_mean, frac_mode )
        class(hhcnn), intent(inout) :: self
        integer,      intent(in)    :: n
        real(sp),     intent(in)    :: x(self%pftsz, self%nk, 2, n), g(self%nk, n), side(n), t(self%nclasses, n)
        real(sp),     intent(in)    :: lr, wd
        real(sp),     intent(out)   :: loss_mean, frac_mode
        integer :: nchunks, ic, ithr, i0, i1
        call self%ensure_chunks(n, nchunks)
        !$omp parallel do default(shared) private(ic,ithr,i0,i1) schedule(dynamic) proc_bind(close)
        do ic = 1,nchunks
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            i0 = (ic-1) * CNN_CHUNK + 1
            i1 = min(n, ic * CNN_CHUNK)
            call self%grad_chunk(ic, ithr, i0, i1, x, g, side, t)
        enddo
        !$omp end parallel do
        call self%reduce_and_step(nchunks, n, lr, wd, loss_mean, frac_mode)
    end subroutine train_batch

    !> the chunks' sums in chunk order, in double, then the Adam step
    subroutine reduce_and_step( self, nchunks, n, lr, wd, loss_mean, frac_mode )
        class(hhcnn), intent(inout) :: self
        integer,      intent(in)    :: nchunks, n
        real(sp),     intent(in)    :: lr, wd
        real(sp),     intent(out)   :: loss_mean, frac_mode
        type(dw_set) :: acc(CNN_NBLK)
        real(dp), allocatable :: dwread(:,:), dtmpl(:)
        real(dp) :: loss, scl, bc1, bc2, dscale
        real(sp) :: scale1(1), mscale1(1), vscale1(1), gscale1(1)
        integer  :: ic, b, nmode
        do b = 1,CNN_NBLK
            allocate(acc(b)%dw(self%conv(b)%cout, 9*self%conv(b)%cin), acc(b)%db(self%conv(b)%cout), source=0.d0)
        enddo
        allocate(dwread(CNN_NREAD, self%rd%cout), dtmpl(self%ntmpl), source=0.d0)
        dscale = 0.d0
        loss   = 0.d0
        nmode  = 0
        do ic = 1,nchunks
            do b = 1,CNN_NBLK
                acc(b)%dw = acc(b)%dw + real(self%gchunk(ic)%blk(b)%dw,dp)
                acc(b)%db = acc(b)%db + real(self%gchunk(ic)%blk(b)%db,dp)
            enddo
            dwread = dwread + real(self%gchunk(ic)%dwread,dp)
            dtmpl  = dtmpl  + real(self%gchunk(ic)%dtmpl,dp)
            dscale = dscale + real(self%gchunk(ic)%dscale,dp)
            loss   = loss + self%gchunk(ic)%loss
            nmode  = nmode + self%gchunk(ic)%nmode
        enddo
        call self%tie_kernel_grads(acc)
        scl       = 1.d0 / real(n,dp)
        loss_mean = real(loss * scl)
        frac_mode = real(nmode) / real(n)
        self%adam_step = self%adam_step + 1
        bc1 = 1.d0 - real(ADAM_BETA1,dp)**self%adam_step
        bc2 = 1.d0 - real(ADAM_BETA2,dp)**self%adam_step
        do b = 1,CNN_NBLK
            call adam_update(self%conv(b)%w, self%conv(b)%mw, self%conv(b)%vw, real(acc(b)%dw * scl), lr, wd, bc1, bc2)
            call adam_update_vec(self%conv(b)%b, self%conv(b)%mb, self%conv(b)%vb, real(acc(b)%db * scl), lr, bc1, bc2)
            deallocate(acc(b)%dw, acc(b)%db)
        enddo
        call adam_update(self%rd%wread, self%rd%mwread, self%rd%vwread, real(dwread * scl), lr, wd, bc1, bc2)
        call adam_update_vec(self%rd%tmpl, self%rd%mtmpl, self%rd%vtmpl, real(dtmpl * scl), lr, bc1, bc2)
        scale1  = self%rd%scale
        mscale1 = self%rd%mscale
        vscale1 = self%rd%vscale
        gscale1 = real(dscale * scl)
        call adam_update_vec(scale1, mscale1, vscale1, gscale1, lr, bc1, bc2)
        self%rd%scale  = scale1(1)
        self%rd%mscale = mscale1(1)
        self%rd%vscale = vscale1(1)
        deallocate(dwread, dtmpl)
    end subroutine reduce_and_step

    !> one Adam step with decoupled weight decay on a weight matrix
    subroutine adam_update( w, m, v, grad, lr, wd, bc1, bc2 )
        real(sp), intent(inout) :: w(:,:), m(:,:), v(:,:)
        real(sp), intent(in)    :: grad(:,:), lr, wd
        real(dp), intent(in)    :: bc1, bc2
        m = ADAM_BETA1 * m + (1. - ADAM_BETA1) * grad
        v = ADAM_BETA2 * v + (1. - ADAM_BETA2) * grad**2
        w = w - lr * (real(real(m,dp) / bc1) / (sqrt(real(real(v,dp) / bc2)) + ADAM_EPS) + wd * w)
    end subroutine adam_update

    !> one Adam step on a bias vector (no decay)
    subroutine adam_update_vec( w, m, v, grad, lr, bc1, bc2 )
        real(sp), intent(inout) :: w(:), m(:), v(:)
        real(sp), intent(in)    :: grad(:), lr
        real(dp), intent(in)    :: bc1, bc2
        m = ADAM_BETA1 * m + (1. - ADAM_BETA1) * grad
        v = ADAM_BETA2 * v + (1. - ADAM_BETA2) * grad**2
        w = w - lr * real(real(m,dp) / bc1) / (sqrt(real(real(v,dp) / bc2)) + ADAM_EPS)
    end subroutine adam_update_vec

    !> the mean loss and the mode agreement of n samples without a step (the held-out evaluation), the same
    !! chunked reduction as training
    subroutine eval_batch( self, n, x, g, side, t, loss_mean, frac_mode )
        class(hhcnn), intent(inout) :: self
        integer,      intent(in)    :: n
        real(sp),     intent(in)    :: x(self%pftsz, self%nk, 2, n), g(self%nk, n), side(n), t(self%nclasses, n)
        real(sp),     intent(out)   :: loss_mean, frac_mode
        real(dp) :: loss
        integer  :: nchunks, ic, ithr, i0, i1, i, nmode
        call self%ensure_chunks(n, nchunks)
        !$omp parallel do default(shared) private(ic,ithr,i0,i1,i,loss) schedule(dynamic) proc_bind(close)
        do ic = 1,nchunks
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            i0 = (ic-1) * CNN_CHUNK + 1
            i1 = min(n, ic * CNN_CHUNK)
            self%gchunk(ic)%loss  = 0.d0
            self%gchunk(ic)%nmode = 0
            do i = i0,i1
                call self%forward_sample(self%work(ithr), x(:,:,:,i), g(:,i), side(i))
                call sample_loss(self%nclasses, self%work(ithr)%prob, t(:,i), loss)
                self%gchunk(ic)%loss = self%gchunk(ic)%loss + loss
                if( maxloc(self%work(ithr)%prob, dim=1) == maxloc(t(:,i), dim=1) ) self%gchunk(ic)%nmode = self%gchunk(ic)%nmode + 1
            enddo
        enddo
        !$omp end parallel do
        loss  = 0.d0
        nmode = 0
        do ic = 1,nchunks
            loss  = loss + self%gchunk(ic)%loss
            nmode = nmode + self%gchunk(ic)%nmode
        enddo
        loss_mean = real(loss / real(n,dp))
        frac_mode = real(nmode) / real(n)
    end subroutine eval_batch

    !> the class probabilities of one sample on thread ithr's workspace (prepare_work first)
    subroutine predict( self, ithr, x, g, side, prob )
        class(hhcnn), intent(inout) :: self
        integer,      intent(in)    :: ithr
        real(sp),     intent(in)    :: x(self%pftsz, self%nk, 2), g(self%nk), side
        real(sp),     intent(out)   :: prob(self%nclasses)
        if( .not. allocated(self%work) ) THROW_HARD('no workspace: prepare_work first; hhcnn%predict')
        if( ithr < 1 .or. ithr > size(self%work) ) THROW_HARD('thread index outside the workspaces; hhcnn%predict')
        call self%forward_sample(self%work(ithr), x, g, side)
        prob = self%work(ithr)%prob
    end subroutine predict

    !> the readout map of thread ithr's last forward pass (predict first), flattened: (CNN_NREAD, pftsz, kread, 2)
    subroutine get_features( self, ithr, f )
        class(hhcnn), intent(in)  :: self
        integer,      intent(in)  :: ithr
        real(sp),     intent(out) :: f(self%nfeat)
        if( .not. allocated(self%work) ) THROW_HARD('no workspace: prepare_work first; hhcnn%get_features')
        if( ithr < 1 .or. ithr > size(self%work) ) THROW_HARD('thread index outside the workspaces; hhcnn%get_features')
        f = reshape(self%work(ithr)%fmap, [self%nfeat])
    end subroutine get_features

    !> the classes' maximised scores and maximising rotations (1..2 pftsz) of thread ithr's last forward pass
    subroutine get_scores( self, ithr, smax, rbest )
        class(hhcnn), intent(in)  :: self
        integer,      intent(in)  :: ithr
        real(sp),     intent(out) :: smax(self%nclasses)
        integer,      intent(out) :: rbest(self%nclasses)
        if( .not. allocated(self%work) ) THROW_HARD('no workspace: prepare_work first; hhcnn%get_scores')
        if( ithr < 1 .or. ithr > size(self%work) ) THROW_HARD('thread index outside the workspaces; hhcnn%get_scores')
        smax  = self%work(ithr)%smax
        rbest = self%work(ithr)%rbest
    end subroutine get_scores

    ! ===== THE WEIGHTS FILE =====

    !> the weights as a binary stream: the magic, the version, the shape, the seed, then every weight array
    subroutine write_weights( self, fname )
        class(hhcnn),  intent(in) :: self
        class(string), intent(in) :: fname
        integer :: funit, io_stat, b
        if( .not. self%exists ) THROW_HARD('network not constructed; hhcnn%write_weights')
        call fopen(funit, fname, 'replace', 'write', iostat=io_stat, access='stream', form='unformatted')
        call fileiochk('hhcnn%write_weights; open '//fname%to_char(), io_stat)
        write(funit) CNN_WEIGHTS_MAGIC, CNN_WEIGHTS_VERSION
        write(funit) self%pftsz, self%nk, self%nclasses, self%width, CNN_NBLK, CNN_NREAD, CNN_KREAD, CNN_MREAD
        write(funit) self%seed
        do b = 1,CNN_NBLK
            write(funit) self%conv(b)%w, self%conv(b)%b
        enddo
        write(funit) self%rd%wread, self%rd%tmpl, self%rd%scale
        call fclose(funit)
    end subroutine write_weights

    !> the network rebuilt from its weights file (the Adam moments start afresh)
    subroutine read_weights( self, fname )
        class(hhcnn),  intent(inout) :: self
        class(string), intent(in)    :: fname
        integer(int64) :: magic, version, seed
        integer :: funit, io_stat, b, pftsz, nk, nclasses, width, nblk, nread, kread, mread
        if( .not. file_exists(fname) ) THROW_HARD('no weights file '//fname%to_char()//'; hhcnn%read_weights')
        call fopen(funit, fname, 'old', 'read', iostat=io_stat, access='stream', form='unformatted')
        call fileiochk('hhcnn%read_weights; open '//fname%to_char(), io_stat)
        read(funit) magic, version
        if( magic /= CNN_WEIGHTS_MAGIC )     THROW_HARD(fname%to_char()//' is no weights file; hhcnn%read_weights')
        if( version /= CNN_WEIGHTS_VERSION ) THROW_HARD(fname%to_char()//' has another version (the readout changed); hhcnn%read_weights')
        read(funit) pftsz, nk, nclasses, width, nblk, nread, kread, mread
        if( nblk /= CNN_NBLK ) THROW_HARD(fname%to_char()//' has another block count; hhcnn%read_weights')
        if( nread /= CNN_NREAD .or. kread /= CNN_KREAD .or. mread /= CNN_MREAD ) THROW_HARD(fname%to_char()//' has another readout shape; hhcnn%read_weights')
        read(funit) seed
        call self%new(pftsz, nk, nclasses, width, seed)
        do b = 1,CNN_NBLK
            read(funit) self%conv(b)%w, self%conv(b)%b
        enddo
        read(funit) self%rd%wread, self%rd%tmpl, self%rd%scale
        call fclose(funit)
    end subroutine read_weights

    !> the weights of another network of the same shape (the early stopping's best epoch)
    subroutine copy_weights_from( self, other )
        class(hhcnn), intent(inout) :: self
        class(hhcnn), intent(in)    :: other
        integer :: b
        if( .not. (self%exists .and. other%exists) ) THROW_HARD('both networks must exist; hhcnn%copy_weights_from')
        do b = 1,CNN_NBLK
            if( any(shape(self%conv(b)%w) /= shape(other%conv(b)%w)) ) THROW_HARD('shapes differ; hhcnn%copy_weights_from')
            self%conv(b)%w = other%conv(b)%w
            self%conv(b)%b = other%conv(b)%b
        enddo
        if( any(shape(self%rd%wread) /= shape(other%rd%wread)) .or. size(self%rd%tmpl) /= size(other%rd%tmpl) )&
            &THROW_HARD('shapes differ; hhcnn%copy_weights_from')
        self%rd%wread = other%rd%wread
        self%rd%tmpl  = other%rd%tmpl
        self%rd%scale = other%rd%scale
    end subroutine copy_weights_from

    !> whether two networks hold the same weights bit for bit
    logical function weights_identical( self, other )
        class(hhcnn), intent(in) :: self, other
        integer :: b
        weights_identical = .false.
        if( .not. (self%exists .and. other%exists) ) return
        do b = 1,CNN_NBLK
            if( any(shape(self%conv(b)%w) /= shape(other%conv(b)%w)) ) return
            if( .not. same_bits2(self%conv(b)%w, other%conv(b)%w) ) return
            if( .not. same_bits1(self%conv(b)%b, other%conv(b)%b) ) return
        enddo
        if( any(shape(self%rd%wread) /= shape(other%rd%wread)) .or. size(self%rd%tmpl) /= size(other%rd%tmpl) ) return
        if( .not. same_bits2(self%rd%wread, other%rd%wread) ) return
        if( .not. same_bits1(self%rd%tmpl, other%rd%tmpl) ) return
        if( transfer(self%rd%scale, 0_int32) /= transfer(other%rd%scale, 0_int32) ) return
        weights_identical = .true.
    end function weights_identical

    logical function same_bits2( a, b )
        real(sp), intent(in) :: a(:,:), b(:,:)
        same_bits2 = all(transfer(a, [0_int32], size(a)) == transfer(b, [0_int32], size(b)))
    end function same_bits2

    logical function same_bits1( a, b )
        real(sp), intent(in) :: a(:), b(:)
        same_bits1 = all(transfer(a, [0_int32], size(a)) == transfer(b, [0_int32], size(b)))
    end function same_bits1

    subroutine kill( self )
        class(hhcnn), intent(inout) :: self
        integer :: b, i
        do b = 1,CNN_NBLK
            if( allocated(self%conv(b)%w)  ) deallocate(self%conv(b)%w)
            if( allocated(self%conv(b)%b)  ) deallocate(self%conv(b)%b)
            if( allocated(self%conv(b)%mw) ) deallocate(self%conv(b)%mw)
            if( allocated(self%conv(b)%vw) ) deallocate(self%conv(b)%vw)
            if( allocated(self%conv(b)%mb) ) deallocate(self%conv(b)%mb)
            if( allocated(self%conv(b)%vb) ) deallocate(self%conv(b)%vb)
            self%conv(b)%cin   = 0
            self%conv(b)%cout  = 0
            self%conv(b)%nkin  = 0
            self%conv(b)%nkout = 0
            self%conv(b)%npos  = 0
        enddo
        if( allocated(self%rd%kb0)    ) deallocate(self%rd%kb0, self%rd%kb1)
        if( allocated(self%rd%wread)  ) deallocate(self%rd%wread, self%rd%mwread, self%rd%vwread)
        if( allocated(self%rd%tmpl)   ) deallocate(self%rd%tmpl, self%rd%mtmpl, self%rd%vtmpl)
        self%rd%cout   = 0
        self%rd%kread  = 0
        self%rd%scale  = 1.
        self%rd%mscale = 0.
        self%rd%vscale = 0.
        if( allocated(self%work) )then
            do i = 1,size(self%work)
                call work_kill(self%work(i))
            enddo
            deallocate(self%work)
        endif
        if( allocated(self%gchunk) )then
            do i = 1,size(self%gchunk)
                call grad_kill(self%gchunk(i))
            enddo
            deallocate(self%gchunk)
        endif
        self%pftsz     = 0
        self%nrots     = 0
        self%nk        = 0
        self%nclasses  = 0
        self%width     = 0
        self%nfeat     = 0
        self%ntmpl     = 0
        self%adam_step = 0
        self%seed      = 0_int64
        self%exists    = .false.
    end subroutine kill


    subroutine work_kill( wk )
        type(cnn_work), intent(inout) :: wk
        integer :: b
        do b = 1,CNN_NBLK
            if( allocated(wk%blk(b)%col)     ) deallocate(wk%blk(b)%col)
            if( allocated(wk%blk(b)%raw)     ) deallocate(wk%blk(b)%raw)
            if( allocated(wk%blk(b)%dnorm)   ) deallocate(wk%blk(b)%dnorm)
            if( allocated(wk%blk(b)%pre)     ) deallocate(wk%blk(b)%pre)
            if( allocated(wk%blk(b)%gated)   ) deallocate(wk%blk(b)%gated)
            if( allocated(wk%blk(b)%pooled)  ) deallocate(wk%blk(b)%pooled)
            if( allocated(wk%blk(b)%dpre)    ) deallocate(wk%blk(b)%dpre)
            if( allocated(wk%blk(b)%dgated)  ) deallocate(wk%blk(b)%dgated)
            if( allocated(wk%blk(b)%dpooled) ) deallocate(wk%blk(b)%dpooled)
            if( allocated(wk%blk(b)%dcol)    ) deallocate(wk%blk(b)%dcol)
        enddo
        if( allocated(wk%fmap)    ) deallocate(wk%fmap, wk%dfmap, wk%fhat, wk%dfhat)
        if( allocated(wk%smax)    ) deallocate(wk%smax, wk%logits, wk%prob, wk%dlogits, wk%rbest, wk%refl)
        if( allocated(wk%g)       ) deallocate(wk%g)
        if( allocated(wk%wnorm)   ) deallocate(wk%wnorm)
        if( allocated(wk%cosm)    ) deallocate(wk%cosm, wk%sinm, wk%cosr, wk%sinr)
    end subroutine work_kill

    subroutine grad_kill( gr )
        type(cnn_grad), intent(inout) :: gr
        integer :: b
        do b = 1,CNN_NBLK
            if( allocated(gr%blk(b)%dw) ) deallocate(gr%blk(b)%dw)
            if( allocated(gr%blk(b)%db) ) deallocate(gr%blk(b)%db)
        enddo
        if( allocated(gr%dwread) ) deallocate(gr%dwread)
        if( allocated(gr%dtmpl)  ) deallocate(gr%dtmpl)
        gr%dscale = 0.
        gr%loss   = 0.d0
        gr%nmode  = 0
    end subroutine grad_kill

end module simple_hhcnn
