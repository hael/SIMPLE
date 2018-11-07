! expanded Fourier transform class for improved cache utilisation
module simple_ft_expanded
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,  only: image
implicit none
private
public :: ft_expanded, ft_exp_reset_tmp_pointers
#include "simple_local_flags.inc"

complex(dp), parameter   :: J = CMPLX(0.0_dp, 1.0_dp, kind=dp)
real(dp),    parameter   :: denom = 0.00075_dp ! denominator for rescaling of cost function
real(dp),    allocatable :: ft_exp_tmpmat_re_2d(:,:,:)
real(dp),    allocatable :: ft_exp_tmpmat_im_2d(:,:,:)
complex(dp), allocatable :: ft_exp_tmp_cmat12(:,:,:)

type :: ft_expanded
    private
    integer              :: lims(3,2)          !< physical limits for the Fourier transform
    integer              :: flims(3,2)         !< shifted limits (2 make transfer 2 GPU painless)
    integer              :: ldim(3)=[1,1,1]    !< logical dimension of originating image
    real                 :: shconst(3)         !< shift constant
    real                 :: hp                 !< high-pass limit
    real                 :: lp                 !< low-pass limit
    real                 :: smpd=0.            !< sampling distance of originating image
    real,    allocatable :: transfmat(:,:,:,:) !< shift transfer matrix
    complex, allocatable :: cmat(:,:,:)        !< Fourier components
    logical              :: existence=.false.  !< existence
  contains
    ! constructors
    procedure, private :: new_1
    procedure, private :: new_2
    procedure, private :: new_3
    generic            :: new => new_1, new_2, new_3
    ! getters
    procedure          :: get_flims
    procedure          :: get_cmat
    ! setters
    procedure          :: set_cmat
    procedure          :: zero
    ! arithmetics
    procedure          :: add
    procedure          :: subtr
    ! modifiers
    procedure          :: shift
    ! calculators
    procedure          :: corr
    procedure          :: corr_shifted_8
    procedure          :: corr_gshifted_8
    procedure          :: corr_fdfshifted_8
    procedure          :: corr_normalize
    ! destructor
    procedure          :: kill
end type ft_expanded

type ftexp_ptr
    type(ft_expanded), pointer :: p => null()
end type ftexp_ptr

type(ftexp_ptr), allocatable :: ft_exp_tmp_cmat12_self1(:)
type(ftexp_ptr), allocatable :: ft_exp_tmp_cmat12_self2(:)

contains

    ! CONSTRUCTORS

    subroutine new_1( self, img, hp, lp )
        class(ft_expanded), intent(inout) :: self
        class(image),       intent(inout) :: img
        real,               intent(in)    :: hp, lp
        integer :: h,k,l,i,hcnt,kcnt,lcnt
        integer :: lplim,hplim,hh,kk,ll,sqarg,phys(3)
        logical :: didft
        ! kill pre-existing object
        call self%kill
        ! set constants
        self%ldim = img%get_ldim()
        if( self%ldim(3) > 1 ) THROW_HARD('only 4 2D images; new_1')
        self%smpd = img%get_smpd()
        self%hp   = hp
        self%lp   = lp
        self%lims = img%loop_lims(1,lp)
        ! shift the limits 2 make transfer 2 GPU painless
        self%flims = 1
        do i=1,3
            self%flims(i,2) = self%lims(i,2) - self%lims(i,1) + 1
        end do
        ! set the squared filter limits
        hplim = img%get_find(hp)
        hplim = hplim*hplim
        lplim = img%get_find(lp)
        lplim = lplim*lplim
        ! set shift constant (shconst)
        do i=1,3
            if( self%ldim(i) == 1 )then
                self%shconst(i) = 0.
                cycle
            endif
            if( is_even(self%ldim(i)) )then
                self%shconst(i) = PI/real(self%ldim(i)/2.)
            else
                self%shconst(i) = PI/real((self%ldim(i)-1)/2.)
            endif
        end do
        ! prepare image
        didft = .false.
        if( .not. img%is_ft() )then
            call img%fft()
            didft = .true.
        endif
        ! allocate instance variables
        allocate(    self%cmat(  self%flims(1,1):self%flims(1,2),&
                                 self%flims(2,1):self%flims(2,2),&
                                 self%flims(3,1):self%flims(3,2)),&
                  self%transfmat(self%flims(1,1):self%flims(1,2),&
                                 self%flims(2,1):self%flims(2,2),&
                                 self%flims(3,1):self%flims(3,2), 3), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In: new_1; simple_ft_expanded, 2",alloc_stat)
        self%cmat      = cmplx(0.,0.)
        self%transfmat = 0.
        hcnt = 0
        do h=self%lims(1,1),self%lims(1,2)
            hh   = h * h
            hcnt = hcnt + 1
            kcnt = 0
            do k=self%lims(2,1),self%lims(2,2)
                kk   = k * k
                kcnt = kcnt + 1
                lcnt = 0
                do l=self%lims(3,1),self%lims(3,2)
                    ll = l * l
                    lcnt = lcnt+1
                    sqarg = hh + kk + ll
                    if( sqarg <= lplim .and. sqarg >= hplim  )then
                        phys = img%comp_addr_phys([h,k,l])
                        self%transfmat(hcnt,kcnt,lcnt,:) = real([h,k,l])*self%shconst
                        self%cmat(hcnt,kcnt,lcnt) = img%get_fcomp([h,k,l],phys)
                    endif
                end do
            end do
        end do
        if( didft ) call img%ifft()
        ! allocate class variables
        if( .not. allocated(ft_exp_tmpmat_re_2d) )then
            allocate( ft_exp_tmpmat_re_2d( 1:self%ldim(1), 1:self%ldim(2), nthr_glob ),&
                      ft_exp_tmpmat_im_2d( 1:self%ldim(1), 1:self%ldim(2), nthr_glob ),&
                      ft_exp_tmp_cmat12(   1:self%ldim(1), 1:self%ldim(2), nthr_glob ),&
                      stat=alloc_stat                               )
            if(alloc_stat/=0)call allocchk("In: new_1; simple_ft_expanded, 3")
        end if
        ! pointer arrays for bookkeeping
        if( .not. allocated(ft_exp_tmp_cmat12_self1) )then
            allocate(ft_exp_tmp_cmat12_self1(nthr_glob), ft_exp_tmp_cmat12_self2(nthr_glob))
        endif
        self%existence = .true.
    end subroutine new_1

    subroutine new_2( self, ldim, smpd, hp, lp )
        class(ft_expanded), intent(inout) :: self
        integer,            intent(in)    :: ldim(3)
        real,               intent(in)    :: smpd
        real,               intent(in)    :: hp
        real,               intent(in)    :: lp
        type(image) :: img
        call img%new(ldim,smpd)
        img = cmplx(0.,0.)
        call self%new_1(img, hp, lp)
        call img%kill
    end subroutine new_2

    subroutine new_3( self, self_in )
        class(ft_expanded), intent(inout) :: self
        class(ft_expanded), intent(in)    :: self_in
        if( self_in%existence )then
            call self%new_2(self_in%ldim, self_in%smpd, self_in%hp, self_in%lp)
            self%cmat = cmplx(0.,0.)
        else
            THROW_HARD('self_in does not exists; new_3')
        endif
    end subroutine new_3

    ! GETTERS

    pure function get_flims( self ) result( flims)
        class(ft_expanded), intent(in) :: self
        integer :: flims(3,2)
        flims = self%flims
    end function get_flims

    pure subroutine get_cmat( self, cmat )
        class(ft_expanded), intent(in) :: self
        complex,            intent(out) :: cmat(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2),self%flims(3,1):self%flims(3,2))
        cmat = self%cmat
    end subroutine get_cmat

    ! SETTERS

    pure subroutine set_cmat( self, cmat )
        class(ft_expanded), intent(inout) :: self
        complex,            intent(in) :: cmat(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2),self%flims(3,1):self%flims(3,2))
        self%cmat = cmat
    end subroutine set_cmat

    pure subroutine zero( self )
        class(ft_expanded), intent(inout) :: self
        self%cmat = cmplx(0.,0.)
    end subroutine zero

    ! ARITHMETICS

    subroutine add( self, self2add, w )
        class(ft_expanded), intent(inout) :: self
        class(ft_expanded), intent(in)    :: self2add
        real, optional,     intent(in)    :: w
        real :: ww
        ww =1.0
        if( present(w) ) ww = w
        if( ww > 0. ) self%cmat = self%cmat + self2add%cmat*ww
    end subroutine add

    subroutine subtr( self, self2subtr, w )
        class(ft_expanded), intent(inout) :: self
        class(ft_expanded), intent(in)    :: self2subtr
        real, optional,     intent(in)    :: w
        real :: ww
        ww = 1.0
        if( present(w) ) ww = w
        if( ww > 0. ) self%cmat = self%cmat-ww*self2subtr%cmat
    end subroutine subtr

    ! MODIFIERS

    subroutine shift( self, shvec, self_out )
        class(ft_expanded), intent(in)    :: self
        real,               intent(in)    :: shvec(3)
        class(ft_expanded), intent(inout) :: self_out
        integer :: hind,kind,lind
        real    :: shvec_here(3), arg
        shvec_here    = shvec
        shvec_here(3) = 0. ! only for 2D images, see constructor (new_1)
        do hind=self%flims(1,1),self%flims(1,2)
            do kind=self%flims(2,1),self%flims(2,2)
                do lind=self%flims(3,1),self%flims(3,2)
                    arg = sum(shvec_here*self%transfmat(hind,kind,lind,:))
                    self_out%cmat(hind,kind,lind) = self%cmat(hind,kind,lind) * cmplx(cos(arg),sin(arg))
                end do
            end do
        end do
    end subroutine shift

    ! CALCULATORS

    function corr( self1, self2 ) result( r )
        class(ft_expanded), intent(in) :: self1, self2
        real :: r,sumasq,sumbsq
        ! corr is real part of the complex mult btw 1 and 2*
        r = sum(real(self1%cmat*conjg(self2%cmat)))
        ! normalisation terms
        sumasq = sum(csq(self1%cmat))
        sumbsq = sum(csq(self2%cmat))
        ! finalise the correlation coefficient
        if( sumasq > 0. .and. sumbsq > 0. )then
            r = r / sqrt(sumasq * sumbsq)
        else
            r = 0.
        endif
    end function corr

    function corr_shifted_8( self1, self2, shvec ) result( r )
        class(ft_expanded), intent(inout) :: self1, self2
        real(dp),           intent(in)    :: shvec(2)
        real(dp) :: r
        integer  :: ithr
        call calc_tmpmat_re(self1, self2, shvec)
        ithr = omp_get_thread_num() + 1
        r = sum(ft_exp_tmpmat_re_2d(self1%flims(1,1):self1%flims(1,2),self1%flims(2,1):self1%flims(2,2),ithr)) / denom
    end function corr_shifted_8

    subroutine corr_gshifted_8( self1, self2, shvec, grad )
        class(ft_expanded), target, intent(inout) :: self1, self2
        real(dp),                   intent(in)    :: shvec(2)
        real(dp),                   intent(out)   :: grad(2)
        real(dp) :: grad1, grad2
        integer  :: ithr
        call calc_tmpmat_im(self1, self2, shvec)
        ithr = omp_get_thread_num() + 1
        grad(1) = sum(ft_exp_tmpmat_im_2d(self1%flims(1,1):self1%flims(1,2),self1%flims(2,1):self1%flims(2,2),ithr) *&
                   &self1%transfmat(self1%flims(1,1):self1%flims(1,2),self1%flims(2,1):self1%flims(2,2),1,1))
        grad(2) = sum(ft_exp_tmpmat_im_2d(self1%flims(1,1):self1%flims(1,2),self1%flims(2,1):self1%flims(2,2),ithr) *&
                   &self1%transfmat(self1%flims(1,1):self1%flims(1,2),self1%flims(2,1):self1%flims(2,2),1,2))
        grad = grad / denom
    end subroutine corr_gshifted_8

    subroutine corr_fdfshifted_8( self1, self2, shvec, f, grad )
        class(ft_expanded), intent(inout) :: self1, self2
        real(dp),           intent(in)    :: shvec(2)
        real(dp),           intent(out)   :: grad(2), f
        real(dp) :: grad1, grad2
        integer  :: hind,kind,ithr
        call calc_tmpmat_re_im(self1, self2, shvec)
        ithr    = omp_get_thread_num() + 1
        f       = sum(ft_exp_tmpmat_re_2d(self1%flims(1,1):self1%flims(1,2),self1%flims(2,1):self1%flims(2,2),ithr)) / denom
        grad(1) = sum(ft_exp_tmpmat_im_2d(self1%flims(1,1):self1%flims(1,2),self1%flims(2,1):self1%flims(2,2),ithr) *&
                   &self1%transfmat(self1%flims(1,1):self1%flims(1,2),self1%flims(2,1):self1%flims(2,2),1,1))
        grad(2) = sum(ft_exp_tmpmat_im_2d(self1%flims(1,1):self1%flims(1,2),self1%flims(2,1):self1%flims(2,2),ithr) *&
                   &self1%transfmat(self1%flims(1,1):self1%flims(1,2),self1%flims(2,1):self1%flims(2,2),1,2))
        grad    = grad / denom
    end subroutine corr_fdfshifted_8

    subroutine corr_normalize( self1, self2, corr )
        class(ft_expanded), intent(inout) :: self1, self2
        real(sp),           intent(inout) :: corr
        real(dp) :: sumasq,sumbsq
        sumasq = sum(csq(self1%cmat))
        sumbsq = sum(csq(self2%cmat))
        corr   = real(real(corr,dp) * denom / sqrt(sumasq * sumbsq), sp)
    end subroutine corr_normalize

    subroutine calc_tmpmat_re(self1, self2, shvec)
        class(ft_expanded), target, intent(inout) :: self1, self2
        real(dp),                   intent(in)    :: shvec(2)
        real(dp) :: arg
        integer  :: hind,kind,ithr
        ithr = omp_get_thread_num() + 1
        if (associated(ft_exp_tmp_cmat12_self1(ithr)%p, self1) .and. associated(ft_exp_tmp_cmat12_self2(ithr)%p, self2)) then
            do hind=self1%flims(1,1),self1%flims(1,2)
                do kind=self1%flims(2,1),self1%flims(2,2)
                    arg  = dot_product(shvec(:), self1%transfmat(hind,kind,1,1:2))
                    ft_exp_tmpmat_re_2d(hind,kind,ithr) = real(ft_exp_tmp_cmat12(hind,kind,ithr) * exp(-J * arg),kind=dp)
                end do
            end do
        else
            do hind=self1%flims(1,1),self1%flims(1,2)
                do kind=self1%flims(2,1),self1%flims(2,2)
                    arg  = dot_product(shvec(:), self1%transfmat(hind,kind,1,1:2))
                    ft_exp_tmp_cmat12(hind,kind,ithr)   = self1%cmat(hind,kind,1) * conjg(self2%cmat(hind,kind,1))
                    ft_exp_tmpmat_re_2d(hind,kind,ithr) = real(ft_exp_tmp_cmat12(hind,kind,ithr) * exp(-J * arg),kind=dp)
                end do
            end do
            ithr = omp_get_thread_num() + 1
            ft_exp_tmp_cmat12_self1(ithr)%p => self1
            ft_exp_tmp_cmat12_self2(ithr)%p => self2
        end if
    end subroutine calc_tmpmat_re

    subroutine calc_tmpmat_im(self1, self2, shvec)
        class(ft_expanded), target, intent(inout) :: self1, self2
        real(dp),                   intent(in)    :: shvec(2)
        real(dp) :: arg
        integer  :: hind,kind,ithr
        ithr = omp_get_thread_num() + 1
        if (associated(ft_exp_tmp_cmat12_self1(ithr)%p, self1) .and. associated(ft_exp_tmp_cmat12_self2(ithr)%p, self2)) then
            do hind=self1%flims(1,1),self1%flims(1,2)
                do kind=self1%flims(2,1),self1%flims(2,2)
                    arg  = dot_product(shvec(:), self1%transfmat(hind,kind,1,1:2))
                    ft_exp_tmpmat_im_2d(hind,kind,ithr) = aimag(ft_exp_tmp_cmat12(hind,kind,ithr) * exp(-J * arg))
                end do
            end do
        else
            do hind=self1%flims(1,1),self1%flims(1,2)
                do kind=self1%flims(2,1),self1%flims(2,2)
                    arg  = dot_product(shvec(:), self1%transfmat(hind,kind,1,1:2))
                    ft_exp_tmp_cmat12(hind,kind,ithr)   = self1%cmat(hind,kind,1) * conjg(self2%cmat(hind,kind,1))
                    ft_exp_tmpmat_im_2d(hind,kind,ithr) = aimag(ft_exp_tmp_cmat12(hind,kind,ithr) * exp(-J * arg))
                end do
            end do
            ithr = omp_get_thread_num() + 1
            ft_exp_tmp_cmat12_self1(ithr)%p => self1
            ft_exp_tmp_cmat12_self2(ithr)%p => self2
        end if
    end subroutine calc_tmpmat_im

    subroutine calc_tmpmat_re_im(self1, self2, shvec)
        class(ft_expanded), target, intent(inout) :: self1, self2
        real(dp),                   intent(in)    :: shvec(2)
        real(dp)    :: arg
        complex(dp) :: tmp
        integer     :: hind,kind,ithr
        ithr = omp_get_thread_num() + 1
        if (associated(ft_exp_tmp_cmat12_self1(ithr)%p, self1) .and. associated(ft_exp_tmp_cmat12_self2(ithr)%p, self2)) then
            do hind=self1%flims(1,1),self1%flims(1,2)
                do kind=self1%flims(2,1),self1%flims(2,2)
                    arg  = dot_product(shvec(:), self1%transfmat(hind,kind,1,1:2))
                    tmp  = ft_exp_tmp_cmat12(hind,kind,ithr) * exp(-J * arg)
                    ft_exp_tmpmat_re_2d(hind,kind,ithr) = real(tmp,kind=dp)
                    ft_exp_tmpmat_im_2d(hind,kind,ithr) = aimag(tmp)
                end do
            end do
        else
            do hind=self1%flims(1,1),self1%flims(1,2)
                do kind=self1%flims(2,1),self1%flims(2,2)
                    arg  = dot_product(shvec(:), self1%transfmat(hind,kind,1,1:2))
                    ft_exp_tmp_cmat12(hind,kind,ithr) = self1%cmat(hind,kind,1) * conjg(self2%cmat(hind,kind,1))
                    tmp  = ft_exp_tmp_cmat12(hind,kind,ithr) * exp(-J * arg)
                    ft_exp_tmpmat_re_2d(hind,kind,ithr) = real(tmp,kind=dp)
                    ft_exp_tmpmat_im_2d(hind,kind,ithr) = aimag(tmp)
                end do
            end do
            ithr = omp_get_thread_num() + 1
            ft_exp_tmp_cmat12_self1(ithr)%p => self1
            ft_exp_tmp_cmat12_self2(ithr)%p => self2
        end if
    end subroutine calc_tmpmat_re_im

    subroutine ft_exp_reset_tmp_pointers
        integer :: ithr
        do ithr=1,nthr_glob
            ft_exp_tmp_cmat12_self1(ithr)%p => null()
            ft_exp_tmp_cmat12_self2(ithr)%p => null()
        end do
    end subroutine ft_exp_reset_tmp_pointers

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(ft_expanded), intent(inout) :: self
        if( self%existence )then
            deallocate(self%cmat, self%transfmat)
            call ft_exp_reset_tmp_pointers
            self%existence = .false.
        endif
    end subroutine kill

end module simple_ft_expanded
