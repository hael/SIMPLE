!@descr: the class-average CNN of scripts/sscnn (junk, secondary structure present, likely protein): its weights file and its forward pass, averaged over the 8 symmetries of the square
! Inference only: the network is trained by scripts/sscnn/cnn_ss.py and exported by export_fortran.py, whose header
! gives the weights file. Arrays are PyTorch's read column-major: an activation (w, h, c) is PyTorch's (c, h, w), and
! SIMPLE's rmat(x, y) is (w, h). A convolution is im2col and one product: sgemm where the build's BLAS stays on the
! calling thread inside teams, explicit loops otherwise; one forward pass runs on one thread.
module simple_cavg_cnn
use, intrinsic :: iso_fortran_env, only: int64, int32
use simple_defs,   only: sp, dp
use simple_error,  only: simple_exception
use simple_string, only: string
use simple_fileio, only: fopen, fclose
use simple_syslib, only: file_exists
implicit none

! the BLAS product (the library the build links)
interface
    subroutine sgemm( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
        character(len=1), intent(in)    :: transa, transb
        integer(kind=4),  intent(in)    :: m, n, k, lda, ldb, ldc
        real(kind=4),     intent(in)    :: alpha, beta
        real(kind=4),     intent(in)    :: a(lda,*), b(ldb,*)
        real(kind=4),     intent(inout) :: c(ldc,*)
    end subroutine sgemm
end interface

public :: cavg_cnn, run_all_cavg_cnn_tests
private
#include "simple_local_flags.inc"

integer(int64), parameter :: CAVG_CNN_MAGIC   = 7741227103_int64
integer(int64), parameter :: CAVG_CNN_VERSION = 1_int64
integer,        parameter :: CAVG_CNN_NBLK    = 5       !< residual blocks, each halving the map
integer,        parameter :: CAVG_CNN_GROUPS  = 8       !< GroupNorm's groups
integer,        parameter :: CAVG_CNN_NSYM    = 8       !< the symmetries of the square the logits are averaged over
real(dp),       parameter :: CAVG_CNN_GN_EPS  = 1.d-5   !< GroupNorm's epsilon (PyTorch's default; the exporter checks it)
#ifdef SIMPLE_BLAS_TEAM_SAFE
logical, parameter :: L_BLAS_TEAM_SAFE = .true.
#else
logical, parameter :: L_BLAS_TEAM_SAFE = .false.
#endif

!> a GroupNorm's per-channel weight and bias
type :: cnn_gn
    real(sp), allocatable :: g(:), b(:)
end type cnn_gn

!> a residual block: 3x3 convolutions w1 (9 cin, cout) at stride 2 and w2 (9 cout, cout), the skip's 1x1 ws (cin, cout)
type :: cnn_block
    integer :: cin = 0, cout = 0
    real(sp), allocatable :: w1(:,:), w2(:,:), ws(:,:)
    type(cnn_gn) :: n1, n2, ns
end type cnn_block

type :: cavg_cnn
    private
    integer  :: width = 0, cin = 0, box = 0
    real(sp) :: smpd  = 0.
    logical  :: l_blas = L_BLAS_TEAM_SAFE
    real(sp), allocatable :: w0(:,:)            !< the stem's 3x3 convolution (9 cin, width)
    type(cnn_gn)          :: n0
    type(cnn_block)       :: blk(CAVG_CNN_NBLK)
    real(sp), allocatable :: wh(:,:), bh(:)     !< the head (channels of the last block, 2) and its bias
    logical  :: exists = .false.
  contains
    procedure :: new
    procedure :: get_box
    procedure :: get_smpd
    procedure :: get_nchannels
    procedure :: predict
    procedure :: kill
    procedure, private :: read_weights
    procedure, private :: logits_sym
    procedure, private :: forward
    procedure, private, nopass :: conv3x3
    procedure, private, nopass :: conv1x1_s2
    procedure, private, nopass :: group_norm
    procedure, private, nopass :: symmetry
    procedure, private, nopass :: gemm
end type cavg_cnn

interface
    module subroutine run_all_cavg_cnn_tests()
    end subroutine run_all_cavg_cnn_tests
end interface

contains

    !> the network of a weights file written by scripts/sscnn/export_fortran.py
    subroutine new( self, fname )
        class(cavg_cnn), intent(inout) :: self
        class(string),   intent(in)    :: fname
        character(len=:), allocatable  :: msg
        logical :: ok
        call self%read_weights(fname, ok, msg)
        if( .not. ok ) THROW_HARD(msg)
    end subroutine new

    integer function get_box( self )
        class(cavg_cnn), intent(in) :: self
        get_box = self%box
    end function get_box

    real function get_smpd( self )
        class(cavg_cnn), intent(in) :: self
        get_smpd = self%smpd
    end function get_smpd

    !> 3: the class average and its even and odd halves; 1: the average alone
    integer function get_nchannels( self )
        class(cavg_cnn), intent(in) :: self
        get_nchannels = self%cin
    end function get_nchannels

    !> probabilities (junk, secondary structure, likely protein) of prepared inputs x(box, box, channels, n): the
    !! logits averaged over the 8 symmetries; P(protein) = sigmoid(l1), P(ss | protein) = sigmoid(l2)
    subroutine predict( self, x, probs )
        class(cavg_cnn), intent(in)  :: self
        real(sp),        intent(in)  :: x(:,:,:,:)
        real(dp),        intent(out) :: probs(:,:)
        real(sp), allocatable :: lsym(:,:,:)
        real(dp) :: l(2), pp, ps
        integer  :: n, itask, i, s
        if( .not. self%exists ) THROW_HARD('network not constructed; cavg_cnn%predict')
        n = size(x,4)
        if( size(x,1) /= self%box .or. size(x,2) /= self%box ) THROW_HARD('inputs are not of the network box; cavg_cnn%predict')
        if( size(x,3) /= self%cin ) THROW_HARD('inputs have another channel count; cavg_cnn%predict')
        if( size(probs,1) /= 3 .or. size(probs,2) /= n ) THROW_HARD('probs must be (3, inputs); cavg_cnn%predict')
        allocate(lsym(2, CAVG_CNN_NSYM, n))
        !$omp parallel do default(shared) private(itask, i, s) schedule(dynamic,1)
        do itask = 1,CAVG_CNN_NSYM * n
            i = (itask - 1) / CAVG_CNN_NSYM + 1
            s = mod(itask - 1, CAVG_CNN_NSYM)
            call one_view(x(:,:,:,i), s, lsym(:,s+1,i))
        enddo
        !$omp end parallel do
        do i = 1,n
            l = 0.d0
            do s = 1,CAVG_CNN_NSYM
                l = l + real(lsym(:,s,i), dp)
            enddo
            l  = l / real(CAVG_CNN_NSYM, dp)
            pp = 1.d0 / (1.d0 + exp(-l(1)))
            ps = 1.d0 / (1.d0 + exp(-l(2)))
            probs(:,i) = [1.d0 - pp, pp * ps, pp * (1.d0 - ps)]
        enddo
        deallocate(lsym)

      contains

        subroutine one_view( xi, s, lv )
            real(sp), intent(in)  :: xi(:,:,:)
            integer,  intent(in)  :: s
            real(sp), intent(out) :: lv(2)
            real(sp), allocatable :: xs(:,:,:)
            allocate(xs, mold=xi)
            call self%symmetry(xi, s, xs)
            call self%forward(xs, lv)
            deallocate(xs)
        end subroutine one_view

    end subroutine predict

    subroutine kill( self )
        class(cavg_cnn), intent(inout) :: self
        integer :: b
        if( allocated(self%w0) )    deallocate(self%w0)
        if( allocated(self%n0%g) )  deallocate(self%n0%g, self%n0%b)
        do b = 1,CAVG_CNN_NBLK
            if( allocated(self%blk(b)%w1) ) deallocate(self%blk(b)%w1, self%blk(b)%w2, self%blk(b)%ws)
            if( allocated(self%blk(b)%n1%g) ) deallocate(self%blk(b)%n1%g, self%blk(b)%n1%b)
            if( allocated(self%blk(b)%n2%g) ) deallocate(self%blk(b)%n2%g, self%blk(b)%n2%b)
            if( allocated(self%blk(b)%ns%g) ) deallocate(self%blk(b)%ns%g, self%blk(b)%ns%b)
            self%blk(b)%cin  = 0
            self%blk(b)%cout = 0
        enddo
        if( allocated(self%wh) ) deallocate(self%wh, self%bh)
        self%width  = 0
        self%cin    = 0
        self%box    = 0
        self%smpd   = 0.
        self%l_blas = L_BLAS_TEAM_SAFE
        self%exists = .false.
    end subroutine kill

    ! ===== THE WEIGHTS FILE =====

    !> the weights file read whole; ok is false with the reason in msg on any mismatch, a short or a long file
    subroutine read_weights( self, fname, ok, msg )
        class(cavg_cnn),               intent(inout) :: self
        class(string),                 intent(in)    :: fname
        logical,                       intent(out)   :: ok
        character(len=:), allocatable, intent(out)   :: msg
        integer :: funit, ios
        call self%kill
        ok = .false.
        if( .not. file_exists(fname) )then
            msg = 'no weights file '//fname%to_char()
            return
        endif
        call fopen(funit, fname, 'old', 'read', iostat=ios, access='stream', form='unformatted')
        if( ios /= 0 )then
            msg = 'cannot open '//fname%to_char()
            return
        endif
        call read_all
        call fclose(funit)
        if( ok )then
            self%exists = .true.
        else
            msg = fname%to_char()//': '//msg
            call self%kill
        endif

      contains

        subroutine read_all
            integer(int64) :: magic, version
            integer  :: width, cin, box, groups, nblk, b, ch(0:CAVG_CNN_NBLK), extra
            real(sp) :: smpd
            read(funit, iostat=ios) magic, version
            if( ios /= 0 )then
                msg = 'too short for a weights file'
                return
            endif
            if( magic /= CAVG_CNN_MAGIC )then
                msg = 'not a class-average CNN weights file'
                return
            endif
            if( version /= CAVG_CNN_VERSION )then
                msg = 'another weights file version than this build reads'
                return
            endif
            read(funit, iostat=ios) width, cin, box, groups, nblk, smpd
            if( ios /= 0 )then
                msg = 'truncated header'
                return
            endif
            if( groups /= CAVG_CNN_GROUPS .or. nblk /= CAVG_CNN_NBLK )then
                msg = 'another architecture (GroupNorm groups or block count)'
                return
            endif
            if( cin /= 1 .and. cin /= 3 )then
                msg = 'input channels must be 1 or 3'
                return
            endif
            if( width < CAVG_CNN_GROUPS .or. mod(width, CAVG_CNN_GROUPS) /= 0 .or. box < 1 .or. smpd <= 0. )then
                msg = 'invalid width, box or pixel size'
                return
            endif
            self%width = width
            self%cin   = cin
            self%box   = box
            self%smpd  = smpd
            ch = [width, 2*width, 4*width, 4*width, 8*width, 8*width]
            allocate(self%w0(9*cin, width), self%n0%g(width), self%n0%b(width))
            if( .not. tensor([3,3,cin,width], self%w0, size(self%w0)) ) return
            if( .not. tensor([width], self%n0%g, width) ) return
            if( .not. tensor([width], self%n0%b, width) ) return
            do b = 1,CAVG_CNN_NBLK
                associate( bk => self%blk(b) )
                bk%cin  = ch(b-1)
                bk%cout = ch(b)
                allocate(bk%w1(9*bk%cin, bk%cout), bk%w2(9*bk%cout, bk%cout), bk%ws(bk%cin, bk%cout))
                allocate(bk%n1%g(bk%cout), bk%n1%b(bk%cout), bk%n2%g(bk%cout), bk%n2%b(bk%cout))
                allocate(bk%ns%g(bk%cout), bk%ns%b(bk%cout))
                if( .not. tensor([3,3,bk%cin,bk%cout],  bk%w1,   size(bk%w1)) ) return
                if( .not. tensor([bk%cout],             bk%n1%g, bk%cout) )     return
                if( .not. tensor([bk%cout],             bk%n1%b, bk%cout) )     return
                if( .not. tensor([3,3,bk%cout,bk%cout], bk%w2,   size(bk%w2)) ) return
                if( .not. tensor([bk%cout],             bk%n2%g, bk%cout) )     return
                if( .not. tensor([bk%cout],             bk%n2%b, bk%cout) )     return
                if( .not. tensor([1,1,bk%cin,bk%cout],  bk%ws,   size(bk%ws)) ) return
                if( .not. tensor([bk%cout],             bk%ns%g, bk%cout) )     return
                if( .not. tensor([bk%cout],             bk%ns%b, bk%cout) )     return
                end associate
            enddo
            allocate(self%wh(ch(CAVG_CNN_NBLK), 2), self%bh(2))
            if( .not. tensor([ch(CAVG_CNN_NBLK),2], self%wh, size(self%wh)) ) return
            if( .not. tensor([2], self%bh, 2) ) return
            read(funit, iostat=ios) extra
            if( ios == 0 )then
                msg = 'data after the last tensor'
                return
            endif
            ok = .true.
        end subroutine read_all

        !> the next tensor: its rank and extents must be ext; values read into a(1:n)
        logical function tensor( ext, a, n )
            integer,  intent(in)  :: ext(:), n
            real(sp), intent(out) :: a(n)
            integer :: rank
            integer, allocatable :: ext_file(:)
            tensor = .false.
            read(funit, iostat=ios) rank
            if( ios /= 0 )then
                msg = 'truncated before a tensor'
                return
            endif
            if( rank /= size(ext) )then
                msg = 'a tensor of another rank than this architecture'
                return
            endif
            allocate(ext_file(rank))
            read(funit, iostat=ios) ext_file
            if( ios /= 0 )then
                msg = 'truncated in a tensor''s extents'
                return
            endif
            if( any(ext_file /= ext) )then
                msg = 'a tensor of other extents than this architecture'
                return
            endif
            read(funit, iostat=ios) a
            if( ios /= 0 )then
                msg = 'truncated in a tensor''s values'
                return
            endif
            tensor = .true.
        end function tensor

    end subroutine read_weights

    ! ===== THE FORWARD PASS =====

    !> the logits of one input x(n, n, channels) averaged over the 8 symmetries, in their order
    subroutine logits_sym( self, x, l )
        class(cavg_cnn), intent(in)  :: self
        real(sp),        intent(in)  :: x(:,:,:)
        real(sp),        intent(out) :: l(2)
        real(sp), allocatable :: xs(:,:,:)
        real(sp) :: lv(2)
        real(dp) :: acc(2)
        integer  :: s
        allocate(xs, mold=x)
        acc = 0.d0
        do s = 0,CAVG_CNN_NSYM-1
            call self%symmetry(x, s, xs)
            call self%forward(xs, lv)
            acc = acc + real(lv, dp)
        enddo
        l = real(acc / real(CAVG_CNN_NSYM, dp), sp)
    end subroutine logits_sym

    !> the logits of one input x(n, n, channels), any n the stride-2 stages accept
    subroutine forward( self, x, l )
        class(cavg_cnn), intent(in)  :: self
        real(sp),        intent(in)  :: x(:,:,:)
        real(sp),        intent(out) :: l(2)
        real(sp), allocatable :: a(:,:), y(:,:), y2(:,:), sk(:,:)
        real(dp), allocatable :: feat(:)
        integer :: n, na, ny, ny2, b, c
        n = size(x,1)
        call self%conv3x3(reshape(x, [n*n, self%cin]), n, self%cin, self%w0, self%width, 1, self%l_blas, a, na)
        call self%group_norm(a, self%n0%g, self%n0%b)
        a = max(a, 0.)
        do b = 1,CAVG_CNN_NBLK
            associate( bk => self%blk(b) )
            call self%conv3x3(a, na, bk%cin, bk%w1, bk%cout, 2, self%l_blas, y, ny)
            call self%group_norm(y, bk%n1%g, bk%n1%b)
            y = max(y, 0.)
            call self%conv3x3(y, ny, bk%cout, bk%w2, bk%cout, 1, self%l_blas, y2, ny2)
            call self%group_norm(y2, bk%n2%g, bk%n2%b)
            call self%conv1x1_s2(a, na, bk%cin, bk%ws, bk%cout, self%l_blas, sk, ny2)
            call self%group_norm(sk, bk%ns%g, bk%ns%b)
            a  = max(y2 + sk, 0.)
            na = ny
            end associate
        enddo
        allocate(feat(size(a,2)))
        do c = 1,size(a,2)
            feat(c) = sum(real(a(:,c), dp)) / real(na*na, dp)
        enddo
        l(1) = real(sum(feat * real(self%wh(:,1), dp)), sp) + self%bh(1)
        l(2) = real(sum(feat * real(self%wh(:,2), dp)), sp) + self%bh(2)
    end subroutine forward

    !> PyTorch's 3x3 convolution with zero padding 1 at stride 1 or 2 (a cross-correlation, no bias): x(n*n, cin) to
    !! y(ny*ny, cout), ny = (n - 1) / stride + 1; w(9 cin, cout) holds the taps kw fastest, then kh, then cin
    subroutine conv3x3( x, n, cin, w, cout, stride, l_blas, y, ny )
        integer,               intent(in)    :: n, cin, cout, stride
        real(sp),              intent(in)    :: x(n*n, cin), w(9*cin, cout)
        logical,               intent(in)    :: l_blas
        real(sp), allocatable, intent(inout) :: y(:,:)
        integer,               intent(out)   :: ny
        real(sp), allocatable :: col(:,:)
        integer :: ci, kh, kw, k, oh, ow, ih, iw, p0
        ny = (n - 1) / stride + 1
        allocate(col(ny*ny, 9*cin))
        do ci = 1,cin
            do kh = 0,2
                do kw = 0,2
                    k = kw + 3*kh + 9*(ci-1) + 1
                    do oh = 1,ny
                        ih = (oh-1) * stride + kh
                        p0 = (oh-1) * ny
                        if( ih < 1 .or. ih > n )then
                            col(p0+1:p0+ny,k) = 0.
                            cycle
                        endif
                        do ow = 1,ny
                            iw = (ow-1) * stride + kw
                            if( iw < 1 .or. iw > n )then
                                col(p0+ow,k) = 0.
                            else
                                col(p0+ow,k) = x((ih-1)*n + iw, ci)
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if( allocated(y) ) deallocate(y)
        allocate(y(ny*ny, cout))
        call gemm(ny*ny, cout, 9*cin, col, w, y, l_blas)
        deallocate(col)
    end subroutine conv3x3

    !> PyTorch's 1x1 convolution at stride 2 (no padding, no bias): x(n*n, cin) to y(ny*ny, cout), w(cin, cout)
    subroutine conv1x1_s2( x, n, cin, w, cout, l_blas, y, ny )
        integer,               intent(in)    :: n, cin, cout
        real(sp),              intent(in)    :: x(n*n, cin), w(cin, cout)
        logical,               intent(in)    :: l_blas
        real(sp), allocatable, intent(inout) :: y(:,:)
        integer,               intent(out)   :: ny
        real(sp), allocatable :: sub(:,:)
        integer :: ci, oh, ow
        ny = (n - 1) / 2 + 1
        allocate(sub(ny*ny, cin))
        do ci = 1,cin
            do oh = 1,ny
                do ow = 1,ny
                    sub((oh-1)*ny + ow, ci) = x(2*(oh-1)*n + 2*(ow-1) + 1, ci)
                enddo
            enddo
        enddo
        if( allocated(y) ) deallocate(y)
        allocate(y(ny*ny, cout))
        call gemm(ny*ny, cout, cin, sub, w, y, l_blas)
        deallocate(sub)
    end subroutine conv1x1_s2

    !> GroupNorm over CAVG_CNN_GROUPS groups of channels of a(npix, c), in place: statistics in double precision,
    !! the variance biased as PyTorch's
    subroutine group_norm( a, g, b )
        real(sp), intent(inout) :: a(:,:)
        real(sp), intent(in)    :: g(:), b(:)
        real(dp) :: mean, var, inv
        integer  :: cg, grp, c0, c1, c
        cg = size(a,2) / CAVG_CNN_GROUPS
        do grp = 1,CAVG_CNN_GROUPS
            c0   = (grp-1) * cg + 1
            c1   = grp * cg
            mean = sum(real(a(:,c0:c1), dp)) / real(size(a,1) * cg, dp)
            var  = sum((real(a(:,c0:c1), dp) - mean)**2) / real(size(a,1) * cg, dp)
            inv  = 1.d0 / sqrt(var + CAVG_CNN_GN_EPS)
            do c = c0,c1
                a(:,c) = real((real(a(:,c), dp) - mean) * inv, sp) * g(c) + b(c)
            enddo
        enddo
    end subroutine group_norm

    !> symmetry s (0..7) of the square applied to x(n, n, channels): s mod 4 quarter turns of the map, after a mirror
    !! of its first axis when s >= 4
    subroutine symmetry( x, s, y )
        real(sp), intent(in)  :: x(:,:,:)
        integer,  intent(in)  :: s
        real(sp), intent(out) :: y(:,:,:)
        integer :: n, c, r
        n = size(x,1)
        do c = 1,size(x,3)
            if( s >= 4 )then
                y(:,:,c) = x(n:1:-1,:,c)
            else
                y(:,:,c) = x(:,:,c)
            endif
            do r = 1,mod(s,4)
                y(:,:,c) = transpose(y(n:1:-1,:,c))
            enddo
        enddo
    end subroutine symmetry

    !> c = a b (a: m x k, b: k x n) on the calling thread: sgemm or the loops down the columns of c
    subroutine gemm( m, n, k, a, b, c, l_blas )
        integer,  intent(in)  :: m, n, k
        real(sp), intent(in)  :: a(m,k), b(k,n)
        real(sp), intent(out) :: c(m,n)
        logical,  intent(in)  :: l_blas
        integer :: j, l
        if( l_blas )then
            call sgemm('N', 'N', m, n, k, 1.0, a, m, b, k, 0.0, c, m)
        else
            do j = 1,n
                c(:,j) = 0.
                do l = 1,k
                    c(:,j) = c(:,j) + a(:,l) * b(l,j)
                enddo
            enddo
        endif
    end subroutine gemm

end module simple_cavg_cnn
