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

complex(dp), parameter     :: J = CMPLX(0.0_dp, 1.0_dp, kind=dp)
real(dp),    parameter     :: denom = 0.00075_dp ! denominator for rescaling of cost function
real(dp),    allocatable   :: ft_exp_tmpmat_re_2d(:,:,:)
real(dp),    allocatable   :: ft_exp_tmpmat_im_2d(:,:,:)
complex(dp), allocatable   :: ft_exp_tmp_cmat12(:,:,:)

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
    procedure          :: new_1
    procedure          :: new_2
    procedure          :: new_3
    generic            :: new => new_1, new_2, new_3
   procedure           :: copy
    ! checkers
    procedure          :: exists
    procedure, private :: same_dims
    generic            :: operator(.eqdims.) => same_dims
    ! getters
    procedure          :: get_ldim
    procedure          :: get_flims
    procedure          :: get_lims
    ! arithmetics
    procedure, private :: assign
    generic            :: assignment(=) => assign
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

    !>  \brief  is a constructor
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

    !>  \brief  is a constructor
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

    !>  \brief  is a constructor
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

    !>  \brief  is a constructor that copies the input object
    subroutine copy( self, self_in )
        class(ft_expanded), intent(inout) :: self
        class(ft_expanded), intent(in)    :: self_in
        if( self_in%existence )then
            call self%new_2(self_in%ldim, self_in%smpd, self_in%hp, self_in%lp)
            self%cmat = self_in%cmat
        else
            THROW_HARD('self_in does not exists; copy')
        endif
    end subroutine copy

    ! CHECKERS

    !>  \brief  checks if an instance exists
    pure function exists( self ) result( yep )
        class(ft_expanded), intent(in) :: self
        logical :: yep
        yep = self%existence
    end function exists

    !>  \brief  checks for same dimensions, overloaded as (.eqdims.)
    pure function same_dims( self1, self2 ) result( yep )
        class(ft_expanded), intent(in) :: self1, self2
        logical :: yep
        yep = all(self1%lims == self2%lims)
    end function same_dims

    ! GETTERS

    !>  \brief  is a getter
    pure function get_ldim( self ) result( ldim )
        class(ft_expanded), intent(in) :: self
        integer :: ldim(3)
        ldim = self%ldim
    end function get_ldim

    !>  \brief  is a getter
    pure function get_flims( self ) result( flims)
        class(ft_expanded), intent(in) :: self
        integer :: flims(3,2)
        flims = self%flims
    end function get_flims

    !>  \brief  is a getter
    pure function get_lims( self ) result( lims)
        class(ft_expanded), intent(in) :: self
        integer :: lims(3,2)
        lims = self%lims
    end function get_lims

    ! ARITHMETICS

    !>  \brief  polymorphic assignment (=)
    subroutine assign( selfout, selfin )
        class(ft_expanded), intent(inout) :: selfout
        class(ft_expanded), intent(in)    :: selfin
        call selfout%copy(selfin)
    end subroutine assign

    !>  \brief  is for ft_expanded summation
    subroutine add( self, self2add, w )
        class(ft_expanded), intent(inout) :: self
        class(ft_expanded), intent(in)    :: self2add
        real, optional,     intent(in)    :: w
        real :: ww
        ww =1.0
        if( present(w) ) ww = w
        if( self%existence )then
            if( self.eqdims.self2add )then
                !$omp parallel workshare proc_bind(close)
                self%cmat = self%cmat + self2add%cmat*ww
                !$omp end parallel workshare
            else
                THROW_HARD('cannot sum ft_expanded objects of different dims; add')
            endif
        else
            self = self2add
            self%cmat = self%cmat*ww
        endif
    end subroutine add

    !>  \brief is for image subtraction,  not overloaded
    subroutine subtr( self, self2subtr, w )
        class(ft_expanded), intent(inout) :: self
        class(ft_expanded), intent(in)    :: self2subtr
        real, optional,     intent(in)    :: w
        real :: ww
        ww = 1.0
        if( present(w) ) ww = w
        if( self%existence )then
            if( self.eqdims.self2subtr )then
                !$omp parallel workshare proc_bind(close)
                self%cmat = self%cmat-ww*self2subtr%cmat
                !$omp end parallel workshare
            else
                THROW_HARD('cannot subtract ft_expanded objects of different dims; subtr')
            endif
        else
            THROW_HARD('the object to subtract from does not exist; subtr')
        endif
    end subroutine subtr

    ! MODIFIERS

    !>  \brief  is 4 shifting an ft_expanded instance
    subroutine shift( self, shvec, self_out )
        class(ft_expanded), intent(in)    :: self
        real,               intent(in)    :: shvec(3)
        class(ft_expanded), intent(inout) :: self_out
        integer :: hind,kind,lind
        real    :: shvec_here(3), arg
        if( self%existence )then
            if( self_out%existence )then
                if( self.eqdims.self_out )then
                    shvec_here = shvec
                    if( self%ldim(3) == 1 ) shvec_here(3) = 0.
                    !$omp parallel do collapse(3) schedule(static) default(shared) &
                    !$omp private(hind,kind,lind,arg) proc_bind(close)
                    do hind=self%flims(1,1),self%flims(1,2)
                        do kind=self%flims(2,1),self%flims(2,2)
                            do lind=self%flims(3,1),self%flims(3,2)
                                arg = sum(shvec_here*self%transfmat(hind,kind,lind,:))
                                self_out%cmat(hind,kind,lind) = self%cmat(hind,kind,lind)*cmplx(cos(arg),sin(arg))
                            end do
                        end do
                    end do
                    !$omp end parallel do
                else
                    write(*,*) 'self     lims: ', self%lims
                    write(*,*) 'self_out lims: ', self_out%lims
                    THROW_HARD('input/output objects have nonconforming dims; shift')
                endif
            else
                THROW_HARD('output object does not exist; shift')
            endif
        else
            THROW_HARD('cannot shift non-existent object; shift')
        endif
    end subroutine shift

    ! CALCULATORS

    !>  \brief  is a correlation calculator
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

    !>  \brief  is a correlation calculator with origin shift of self2, double precision
    function corr_shifted_8( self1, self2, shvec ) result( r )
        class(ft_expanded), intent(inout) :: self1, self2
        real(dp),           intent(in)    :: shvec(2)
        real(dp) :: r
        integer  :: hind,kind, ithr
        if ( self1.eqdims.self2 ) then
            call calc_tmpmat_re(self1, self2, shvec)
            r = 0.0_dp
            !$omp parallel do collapse(2) schedule(static) reduction(+:r) private(hind,kind,ithr) default(shared)
            do hind=self1%flims(1,1),self1%flims(1,2)
                do kind=self1%flims(2,1),self1%flims(2,2)
                    ithr = omp_get_thread_num() + 1
                    r = r + ft_exp_tmpmat_re_2d(hind,kind,ithr)
                end do
            end do
            !$omp end parallel do
            ! finalise the correlation coefficient
            r = r / denom
        else
            write(*,*) 'self1 flims: ', self1%flims(1,1), self1%flims(1,2), self1%flims(2,1),&
                self1%flims(2,2), self1%flims(3,1), self1%flims(3,2)
            write(*,*) 'self2 flims: ', self2%flims(1,1), self2%flims(1,2), self2%flims(2,1),&
                self2%flims(2,2), self2%flims(3,1), self2%flims(3,2)
            THROW_HARD('cannot correlate expanded_ft:s with different dims; corr_shifted_8')
        endif ! end of if( self1.eqdims.self2 ) statement
    end function corr_shifted_8

    !>  \brief  is a correlation calculator with origin shift of self2, double precision
    subroutine corr_gshifted_8( self1, self2, shvec, grad )
        class(ft_expanded), target, intent(inout) :: self1, self2
        real(dp),                   intent(in)    :: shvec(2)
        real(dp),                   intent(out)   :: grad(2)
        real(dp) :: grad1, grad2
        integer  :: hind,kind,ithr
        if( self1.eqdims.self2 ) then
            call calc_tmpmat_im(self1, self2, shvec)
            grad1 = 0.0_dp
            grad2 = 0.0_dp
            !$omp parallel do collapse(2) schedule(static) reduction(+:grad1,grad2) private(hind,kind,ithr) default(shared)
            do hind=self1%flims(1,1),self1%flims(1,2)
                do kind=self1%flims(2,1),self1%flims(2,2)
                    ithr = omp_get_thread_num() + 1
                    grad1 = grad1 + ft_exp_tmpmat_im_2d(hind,kind,ithr) * self1%transfmat(hind,kind,1,1)
                    grad2 = grad2 + ft_exp_tmpmat_im_2d(hind,kind,ithr) * self1%transfmat(hind,kind,1,2)
                end do
            end do
            !$omp end parallel do
            grad(1) = grad1
            grad(2) = grad2
            ! finalise the correlation coefficient
            grad = grad / denom
        else
            write(*,*) 'self1 flims: ', self1%flims(1,1), self1%flims(1,2), self1%flims(2,1),&
                self1%flims(2,2), self1%flims(3,1), self1%flims(3,2)
            write(*,*) 'self2 flims: ', self2%flims(1,1), self2%flims(1,2), self2%flims(2,1),&
                self2%flims(2,2), self2%flims(3,1), self2%flims(3,2)
            THROW_HARD('cannot correlate expanded_ft:s with different dims; :corr_gshifted_8')
        endif ! end of if( self1.eqdims.self2 ) statement
    end subroutine corr_gshifted_8

    !>  \brief  is a correlation calculator with origin shift of self2, double precision
    subroutine corr_fdfshifted_8( self1, self2, shvec, f, grad )
        class(ft_expanded), intent(inout) :: self1, self2
        real(dp),           intent(in)    :: shvec(2)
        real(dp),           intent(out)   :: grad(2), f
        real(dp) :: grad1, grad2
        integer  :: hind,kind,ithr
        if ( self1.eqdims.self2 ) then
            call calc_tmpmat_re_im(self1, self2, shvec)
            ! corr is real part of the complex mult btw 1 and 2*
            f     = 0.0_dp
            grad1 = 0.0_dp
            grad2 = 0.0_dp
            !$omp parallel do collapse(2) schedule(static) reduction(+:f,grad1,grad2) private(hind,kind,ithr) default(shared)
            do hind=self1%flims(1,1),self1%flims(1,2)
                do kind=self1%flims(2,1),self1%flims(2,2)
                    ithr  = omp_get_thread_num() + 1
                    f     = f     + ft_exp_tmpmat_re_2d(hind,kind,ithr)
                    grad1 = grad1 + ft_exp_tmpmat_im_2d(hind,kind,ithr) * self1%transfmat(hind,kind,1,1)
                    grad2 = grad2 + ft_exp_tmpmat_im_2d(hind,kind,ithr) * self1%transfmat(hind,kind,1,2)
                end do
            end do
            !$omp end parallel do
            grad(1) = grad1
            grad(2) = grad2
            ! finalise the correlation coefficient
            f    = f    / denom
            grad = grad / denom
        else
            write(*,*) 'self1 flims: ', self1%flims(1,1), self1%flims(1,2), self1%flims(2,1),&
                self1%flims(2,2), self1%flims(3,1), self1%flims(3,2)
            write(*,*) 'self2 flims: ', self2%flims(1,1), self2%flims(1,2), self2%flims(2,1),&
                self2%flims(2,2), self2%flims(3,1), self2%flims(3,2)
            THROW_HARD('cannot correlate expanded_ft:s with different dims; corr_fdfshifted_8')
        endif ! end of if( self1.eqdims.self2 ) statement
    end subroutine corr_fdfshifted_8

    !> \brief  correctly normalize correlations after minimization
    subroutine corr_normalize( self1, self2, corr )
        class(ft_expanded), intent(inout) :: self1, self2
        real(sp),           intent(inout) :: corr
        real(dp)                          :: sumasq,sumbsq
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
            !$omp parallel do collapse(2) schedule(static) default(shared) &
            !$omp private(hind,kind,arg,ithr) proc_bind(close)
            do hind=self1%flims(1,1),self1%flims(1,2)
                do kind=self1%flims(2,1),self1%flims(2,2)
                    ithr = omp_get_thread_num() + 1
                    arg  = dot_product(shvec(:), self1%transfmat(hind,kind,1,1:2))
                    ft_exp_tmpmat_re_2d(hind,kind,ithr) = real(ft_exp_tmp_cmat12(hind,kind,ithr) * exp(-J * arg),kind=dp)
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do collapse(2) schedule(static) default(shared) &
            !$omp private(hind,kind,arg,ithr) proc_bind(close)
            do hind=self1%flims(1,1),self1%flims(1,2)
                do kind=self1%flims(2,1),self1%flims(2,2)
                    ithr = omp_get_thread_num() + 1
                    arg  = dot_product(shvec(:), self1%transfmat(hind,kind,1,1:2))
                    ft_exp_tmp_cmat12(hind,kind,ithr)   = self1%cmat(hind,kind,1) * conjg(self2%cmat(hind,kind,1))
                    ft_exp_tmpmat_re_2d(hind,kind,ithr) = real(ft_exp_tmp_cmat12(hind,kind,ithr) * exp(-J * arg),kind=dp)
                end do
            end do
            !$omp end parallel do
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
            !$omp parallel do collapse(2) schedule(static) default(shared) &
            !$omp private(hind,kind,arg,ithr) proc_bind(close)
            do hind=self1%flims(1,1),self1%flims(1,2)
                do kind=self1%flims(2,1),self1%flims(2,2)
                    ithr = omp_get_thread_num() + 1
                    arg  = dot_product(shvec(:), self1%transfmat(hind,kind,1,1:2))
                    ft_exp_tmpmat_im_2d(hind,kind,ithr) = aimag(ft_exp_tmp_cmat12(hind,kind,ithr) * exp(-J * arg))
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do collapse(2) schedule(static) default(shared) &
            !$omp private(hind,kind,arg,ithr) proc_bind(close)
            do hind=self1%flims(1,1),self1%flims(1,2)
                do kind=self1%flims(2,1),self1%flims(2,2)
                    ithr = omp_get_thread_num() + 1
                    arg  = dot_product(shvec(:), self1%transfmat(hind,kind,1,1:2))
                    ft_exp_tmp_cmat12(hind,kind,ithr)   = self1%cmat(hind,kind,1) * conjg(self2%cmat(hind,kind,1))
                    ft_exp_tmpmat_im_2d(hind,kind,ithr) = aimag(ft_exp_tmp_cmat12(hind,kind,ithr) * exp(-J * arg))
                end do
            end do
            !$omp end parallel do
            ithr = omp_get_thread_num() + 1
            ft_exp_tmp_cmat12_self1(ithr)%p => self1
            ft_exp_tmp_cmat12_self2(ithr)%p => self2
        end if
    end subroutine calc_tmpmat_im

    subroutine calc_tmpmat_re_im(self1, self2, shvec)
        class(ft_expanded), target, intent(inout) :: self1, self2 !< instances
        real(dp),                   intent(in)    :: shvec(2)
        real(dp)    :: arg
        complex(dp) :: tmp
        integer     :: hind,kind,ithr
        ithr = omp_get_thread_num() + 1
        if (associated(ft_exp_tmp_cmat12_self1(ithr)%p, self1) .and. associated(ft_exp_tmp_cmat12_self2(ithr)%p, self2)) then
            !$omp parallel do collapse(2) schedule(static) default(shared) &
            !$omp private(hind,kind,arg,tmp,ithr) proc_bind(close)
            do hind=self1%flims(1,1),self1%flims(1,2)
                do kind=self1%flims(2,1),self1%flims(2,2)
                    ithr = omp_get_thread_num() + 1
                    arg  = dot_product(shvec(:), self1%transfmat(hind,kind,1,1:2))
                    tmp  = ft_exp_tmp_cmat12(hind,kind,ithr) * exp(-J * arg)
                    ft_exp_tmpmat_re_2d(hind,kind,ithr) = real(tmp,kind=dp)
                    ft_exp_tmpmat_im_2d(hind,kind,ithr) = aimag(tmp)
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do collapse(2) schedule(static) default(shared) &
            !$omp private(hind,kind,arg,tmp,ithr) proc_bind(close)
            do hind=self1%flims(1,1),self1%flims(1,2)
                do kind=self1%flims(2,1),self1%flims(2,2)
                    ithr = omp_get_thread_num() + 1
                    arg  = dot_product(shvec(:), self1%transfmat(hind,kind,1,1:2))
                    ft_exp_tmp_cmat12(hind,kind,ithr) = self1%cmat(hind,kind,1) * conjg(self2%cmat(hind,kind,1))
                    tmp  = ft_exp_tmp_cmat12(hind,kind,ithr) * exp(-J * arg)
                    ft_exp_tmpmat_re_2d(hind,kind,ithr) = real(tmp,kind=dp)
                    ft_exp_tmpmat_im_2d(hind,kind,ithr) = aimag(tmp)
                end do
            end do
            !$omp end parallel do
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
