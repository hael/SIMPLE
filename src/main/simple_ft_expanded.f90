! expanded Fourier transform class for improved cache utilisation
module simple_ft_expanded
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,  only: image
implicit none
private
public :: ft_expanded
#include "simple_local_flags.inc"

real(dp),         parameter   :: denom = 0.00075_dp  ! denominator for rescaling of cost function
real(dp),         parameter   :: num   = 1.0d8       ! numerator for rescaling of cost function

type :: ft_expanded
    private
    complex, allocatable :: cmat(:,:,:)        !< Fourier components
    real,    allocatable :: transfmat(:,:,:,:) !< shift transfer matrix (TODO: make this global to save memory)
    real                 :: shconst(3)         !< shift constant
    real                 :: hp                 !< high-pass limit
    real                 :: lp                 !< low-pass limit
    real                 :: smpd  = 0.         !< sampling distance of originating image
    real                 :: sumsq = 0.         !< sum of squares
    integer              :: lims(3,2)          !< physical limits for the Fourier transform
    integer              :: flims(3,2)         !< shifted limits (2 make transfer 2 GPU painless)
    integer              :: ldim(3)=[1,1,1]    !< logical dimension of originating image
    integer              :: flims_nyq(3,2)     !< nyquist limits
    logical              :: existence=.false.  !< existence
  contains
    ! constructors
    procedure          :: new
    ! getters
    procedure          :: get_flims
    procedure          :: get_lims
    procedure          :: get_ldim
    procedure          :: get_flims_nyq
    procedure          :: get_cmat
    procedure          :: get_cmat_ptr
    procedure          :: get_transfmat_ptr
    procedure          :: get_sumsq
    ! setters
    procedure          :: set_cmat
    procedure          :: zero
    ! arithmetics
    procedure          :: add
    procedure          :: subtr
    ! modifiers
    procedure          :: shift
    ! calculators
    procedure, private :: calc_sumsq
    procedure          :: corr
    procedure, private :: corr_normalize_sp
    procedure, private :: corr_normalize_dp
    generic            :: corr_normalize => corr_normalize_sp, corr_normalize_dp
    ! destructor
    procedure          :: kill
end type ft_expanded

contains

    ! CONSTRUCTORS

    subroutine new( self, img, hp, lp, fetch_comps )
        class(ft_expanded), intent(inout) :: self
        class(image),       intent(inout) :: img
        real,               intent(in)    :: hp, lp
        logical,            intent(in)    :: fetch_comps
        real    :: lp_nyq
        integer :: h,k,l,i,hcnt,kcnt,lcnt
        integer :: lplim,hplim,hh,kk,ll,sqarg,phys(3)
        logical :: didft
        ! kill pre-existing object
        call self%kill
        ! set constants
        self%ldim      = img%get_ldim()
        if( self%ldim(3) > 1 ) THROW_HARD('only for 2D images; new_1')
        self%smpd      = img%get_smpd()
        self%hp        = hp
        lp_nyq         = 2.*self%smpd
        self%lp        = max(lp, lp_nyq)
        self%lims      = img%loop_lims(1,self%lp)
        self%flims_nyq = img%loop_lims(1,lp_nyq)
        ! shift the limits 2 make transfer 2 GPU painless
        self%flims     = 1
        do i=1,3
            self%flims(i,2)     = self%lims(i,2)      - self%lims(i,1)      + 1
            self%flims_nyq(i,2) = self%flims_nyq(i,2) - self%flims_nyq(i,1) + 1
        end do
        self%flims_nyq(:,1) = 1
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
        if( .not. img%is_ft() .and. fetch_comps )then
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
                    if( sqarg < 0.5 )then
                        cycle ! excludes zero
                    elseif( sqarg <= lplim .and. sqarg >= hplim  )then
                        phys = img%comp_addr_phys([h,k,l])
                        self%transfmat(hcnt,kcnt,lcnt,:) = real([h,k,l])*self%shconst
                        if( fetch_comps ) self%cmat(hcnt,kcnt,lcnt) = img%get_fcomp([h,k,l],phys)
                    end if
                end do
            end do
        end do
        if( fetch_comps ) call self%calc_sumsq
        if( didft ) call img%ifft()
        ! allocate class variables
        self%existence = .true.
    end subroutine new

    ! GETTERS

    pure function get_flims( self ) result( flims )
        class(ft_expanded), intent(in) :: self
        integer :: flims(3,2)
        flims = self%flims
    end function get_flims

    pure function get_lims( self ) result( lims )
        class(ft_expanded), intent(in) :: self
        integer :: lims(3,2)
        lims = self%lims
    end function get_lims

    pure function get_ldim( self ) result( ldim )
        class(ft_expanded), intent(in) :: self
        integer :: ldim(3)
        ldim = self%ldim
    end function get_ldim

    pure function get_flims_nyq( self ) result( flims_nyq )
        class(ft_expanded), intent(in) :: self
        integer :: flims_nyq(3,2)
        flims_nyq = self%flims_nyq
    end function get_flims_nyq

    pure subroutine get_cmat( self, cmat )
        class(ft_expanded), intent(in) :: self
        complex,            intent(out) :: cmat(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2),self%flims(3,1):self%flims(3,2))
        cmat = self%cmat
    end subroutine get_cmat

    pure subroutine get_cmat_ptr( self, cmat_ptr )
        class(ft_expanded), target,  intent(inout) :: self
        complex,            pointer, intent(out)   :: cmat_ptr(:,:,:)
        cmat_ptr => self%cmat
    end subroutine get_cmat_ptr

    pure subroutine get_transfmat_ptr( self, transfmat_ptr )
        class(ft_expanded), target,  intent(inout) :: self
        real,               pointer, intent(out)   :: transfmat_ptr(:,:,:,:)
        transfmat_ptr => self%transfmat
    end subroutine get_transfmat_ptr

    pure real function get_sumsq( self )
        class(ft_expanded), intent(in) :: self
        get_sumsq = self%sumsq
    end function get_sumsq

    ! SETTERS

    pure subroutine set_cmat( self, cmat )
        class(ft_expanded), intent(inout) :: self
        complex,            intent(in) :: cmat(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2),self%flims(3,1):self%flims(3,2))
        self%cmat = cmat
        call self%calc_sumsq
    end subroutine set_cmat

    pure subroutine zero( self )
        class(ft_expanded), intent(inout) :: self
        self%cmat  = cmplx(0.,0.)
        self%sumsq = 0.
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
        call self%calc_sumsq
    end subroutine add

    subroutine subtr( self, self2subtr, w )
        class(ft_expanded), intent(inout) :: self
        class(ft_expanded), intent(in)    :: self2subtr
        real, optional,     intent(in)    :: w
        real :: ww
        ww = 1.0
        if( present(w) ) ww = w
        if( ww > 0. ) self%cmat = self%cmat-ww*self2subtr%cmat
        call self%calc_sumsq
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
        call self_out%calc_sumsq
    end subroutine shift

    ! CALCULATORS

    pure subroutine calc_sumsq( self )
        class(ft_expanded), intent(inout) :: self
        self%sumsq =                 sum(csq(self%cmat(                  1,1:self%flims(2,2)-1,1)))
        self%sumsq = self%sumsq +    sum(csq(self%cmat(  self%flims(1,2)  ,1:self%flims(2,2)-1,1)))
        self%sumsq = self%sumsq + 2.*sum(csq(self%cmat(2:self%flims(1,2)-1,1:self%flims(2,2)-1,1)))
    end subroutine calc_sumsq

    pure function corr( self1, self2 ) result( r )
        class(ft_expanded), intent(in) :: self1, self2
        real :: r,sumasq,sumbsq
        ! corr is real part of the complex mult btw 1 and 2*
        r =        sum(real(self1%cmat(                   1,1:self1%flims(2,2)-1,1)*&
                      conjg(self2%cmat(                   1,1:self1%flims(2,2)-1,1))))
        r = r +    sum(real(self1%cmat(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)*&
                      conjg(self2%cmat(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1))))
        r = r + 2.*sum(real(self1%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)*&
                      conjg(self2%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1))))
        ! normalise the correlation coefficient
        if( self1%sumsq > 0. .and. self2%sumsq > 0. )then
            r = r / sqrt(self1%sumsq * self2%sumsq)
        else
            r = 0.
        endif
    end function corr

    elemental subroutine corr_normalize_sp( self1, self2, corr )
        class(ft_expanded), intent(in)    :: self1, self2
        real(sp),           intent(inout) :: corr
        corr   = corr / real(num,sp)
    end subroutine corr_normalize_sp

    elemental subroutine corr_normalize_dp( self1, self2, corr )
        class(ft_expanded), intent(in)    :: self1, self2
        real(dp),           intent(inout) :: corr
        corr   = corr / num
    end subroutine corr_normalize_dp

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(ft_expanded), intent(inout) :: self
        if( self%existence )then
            deallocate(self%cmat, self%transfmat)
            self%existence = .false.
        endif
    end subroutine kill

end module simple_ft_expanded
