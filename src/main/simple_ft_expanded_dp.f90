! expanded Fourier transform class for improved cache utilisation
module simple_ft_expanded_dp
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,  only: image
implicit none
private
public :: ft_expanded_dp
#include "simple_local_flags.inc"

complex(dp), parameter :: JJ = dcmplx(0.d0, 1.d0)

type :: ft_expanded_dp
    private
    complex(dp), allocatable :: cmat(:,:)      !< Fourier components
    logical,     allocatable :: bandmsk(:,:)   !< for band-passed correlation
    real(dp)             :: shconst(2)
    real                 :: hp                 !< high-pass limit
    real                 :: lp                 !< low-pass limit
    real                 :: smpd  = 0.         !< sampling distance of originating image
    integer              :: flims(3,2)         !< shifted limits
    integer              :: ldim(3)=[1,1,1]    !< logical dimension of originating image
    logical              :: existence=.false.  !< existence
  contains
    ! constructors
    procedure          :: new
    ! getters
    procedure          :: get_flims
    procedure          :: get_ldim
    procedure          :: get_cmat
    ! setters
    procedure          :: set_cmat
    procedure          :: zero
    procedure          :: copy
    ! arithmetics
    procedure          :: add
    procedure          :: subtr
    procedure          :: div
    ! modifiers
    procedure          :: shift
    procedure          :: normalize_mat
    ! calculators
    procedure          :: corr
    procedure          :: corr_unnorm
    procedure          :: gen_grad_noshift
    procedure          :: f
    ! procedure          :: df3
    procedure          :: fdf
    ! destructor
    procedure          :: kill
end type ft_expanded_dp

contains

    ! CONSTRUCTORS

    subroutine new( self, img, ldim, hp, lp, fetch_comps, bfac )
        class(ft_expanded_dp), intent(inout) :: self
        class(image),       intent(inout) :: img
        integer,            intent(in)    :: ldim(2)
        real,               intent(in)    :: hp, lp
        logical,            intent(in)    :: fetch_comps
        real, optional,     intent(in)    :: bfac
        real(dp), allocatable :: wh(:)
        real(dp) :: spafreqh,spafreqk,wk
        real     :: lp_nyq
        integer  :: h,k,i,lplim,hplim,sqarg
        logical  :: didft, l_bfac
        ! kill pre-existing object
        call self%kill
        ! set constants
        self%ldim = img%get_ldim()
        if( self%ldim(3) > 1 ) THROW_HARD('only for 2D images; new_1')
        self%smpd  = img%get_smpd()
        self%hp    = hp
        lp_nyq     = 2.*self%smpd
        self%lp    = max(lp, lp_nyq)
        self%flims = img%loop_lims(1,self%lp)
        ! set the squared filter limits
        hplim = img%get_find(hp)
        hplim = hplim*hplim
        lplim = img%get_find(self%lp)
        lplim = lplim*lplim
        ! set shift constant
        do i=1,2
            if( self%ldim(i) > 1 )then
                if( is_even(self%ldim(i)) )then
                    self%shconst(i) = DPI/real(ldim(i)/2,dp)
                else
                    self%shconst(i) = DPI/real((ldim(i)-1)/2,dp)
                endif
            endif
        end do
        ! bfactor
        l_bfac = .false.
        if( present(bfac) )then
            if( bfac > 1.e-6) l_bfac = .true.
        endif
        ! prepare image
        didft = .false.
        if( .not. img%is_ft() .and. fetch_comps )then
            call img%fft()
            didft = .true.
        endif
        ! allocate instance variables
        allocate(self%cmat(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2)),&
                 &self%bandmsk(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2)),&
                 &stat=alloc_stat)
        self%cmat    = cmplx(0.,0.)
        self%bandmsk = .false.
        if(alloc_stat.ne.0)call allocchk("In: new_1; simple_ft_expanded",alloc_stat)
        ! init matrices
        if( fetch_comps .and. l_bfac )then
            allocate(wh(self%flims(1,1):self%flims(1,2)))
            do h=self%flims(1,1),self%flims(1,2)
                spafreqh = real(h,dp)/real(self%ldim(1),dp)/real(self%smpd,dp)
                wh(h)    = exp(-spafreqh*spafreqh*bfac/4.d0)
            enddo
        endif
        do k=self%flims(2,1),self%flims(2,2)
            if( fetch_comps .and. l_bfac )then
                spafreqk = real(k,dp)/real(self%ldim(2),dp)/real(self%smpd,dp)
                wk = exp(-spafreqk*spafreqk*bfac/4.d0)
            endif
            do h=self%flims(1,1),self%flims(1,2)
                sqarg = h*h + k*k
                if( sqarg < 1 )then
                    cycle ! excludes zero
                elseif( sqarg <= lplim .and. sqarg >= hplim  )then
                    if( fetch_comps )then
                        if( l_bfac )then
                            self%cmat(h,k) = wh(h) * wk * img%get_fcomp2D(h,k)
                        else
                            self%cmat(h,k) = img%get_fcomp2D(h,k)
                        endif
                    endif
                    self%bandmsk(h,k) = .true.
                end if
            end do
        end do
        if( didft ) call img%ifft()
        self%existence = .true.
    end subroutine new

    ! GETTERS

    pure function get_flims( self ) result( flims )
        class(ft_expanded_dp), intent(in) :: self
        integer :: flims(3,2)
        flims = self%flims
    end function get_flims

    pure function get_ldim( self ) result( ldim )
        class(ft_expanded_dp), intent(in) :: self
        integer :: ldim(3)
        ldim = self%ldim
    end function get_ldim

    pure subroutine get_cmat( self, cmat )
        class(ft_expanded_dp), intent(in) :: self
        complex(dp),        intent(out) :: cmat(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2))
        cmat = self%cmat
    end subroutine get_cmat

    ! SETTERS

    subroutine set_cmat( self, cmat )
        class(ft_expanded_dp), intent(inout) :: self
        complex(dp),        intent(in)    :: cmat(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2))
        self%cmat = cmat
    end subroutine set_cmat

    pure subroutine zero( self )
        class(ft_expanded_dp), intent(inout) :: self
        self%cmat  = dcmplx(0.d0,0.d0)
    end subroutine zero

    ! not a constructor
    pure subroutine copy( self, self2copy )
        class(ft_expanded_dp), intent(inout) :: self, self2copy
        self%hp    = self2copy%hp
        self%lp    = self2copy%lp
        self%smpd  = self2copy%smpd
        self%ldim  = self2copy%ldim
        self%flims = self2copy%flims
        self%cmat  = self2copy%cmat
        self%shconst = self2copy%shconst
        self%bandmsk = self2copy%bandmsk
        self%existence = self%existence
    end subroutine copy

    ! ARITHMETICS

    subroutine add( self, self2add, w )
        class(ft_expanded_dp), intent(inout) :: self
        class(ft_expanded_dp), intent(in)    :: self2add
        real(dp), optional, intent(in)    :: w
        if( present(w) )then
            self%cmat  = self%cmat + self2add%cmat*w
        else
            self%cmat  = self%cmat + self2add%cmat
        endif
    end subroutine add

    subroutine subtr( self, self2subtr, w )
        class(ft_expanded_dp), intent(inout) :: self
        class(ft_expanded_dp), intent(in)    :: self2subtr
        real(dp), optional,  intent(in)    :: w
        if( present(w) )then
             self%cmat = self%cmat-w*self2subtr%cmat
        else
            self%cmat  = self%cmat-self2subtr%cmat
        endif
    end subroutine subtr

    subroutine div( self, val )
        class(ft_expanded_dp), intent(inout) :: self
        real(dp),              intent(in)    :: val
        if( abs(val) < TINY ) return
        self%cmat  = self%cmat / val
    end subroutine div

    ! MODIFIERS

    subroutine shift( self, shvec, self_out )
        class(ft_expanded_dp),           intent(inout) :: self
        real(dp),                        intent(in)    :: shvec(2)
        class(ft_expanded_dp), optional, intent(inout) :: self_out
        integer  :: h,k
        real(dp) :: ch(self%flims(1,1):self%flims(1,2)),sh(self%flims(1,1):self%flims(1,2))
        real(dp) :: argh,argk,ck,sk
        do h = self%flims(1,1),self%flims(1,2)
            argh  = real(h,dp) * shvec(1) * self%shconst(1)
            ch(h) = cos(argh)
            sh(h) = sin(argh)
        enddo
        if( present(self_out) )then
            do k = self%flims(2,1),self%flims(2,2)
                argk  = real(k,dp) * shvec(2) * self%shconst(2)
                ck    = cos(argk)
                sk    = sin(argk)
                do h = self%flims(1,1),self%flims(1,2)
                    if( self%bandmsk(h,k) )then
                        self_out%cmat(h,k) = self%cmat(h,k) * dcmplx(ck*ch(h)-sk*sh(h), (sk*ch(h)+ck*sh(h)))
                    endif
                end do
            end do
        else
            do k = self%flims(2,1),self%flims(2,2)
                argk  = real(k,dp) * shvec(2) * self%shconst(2)
                ck    = cos(argk)
                sk    = sin(argk)
                do h = self%flims(1,1),self%flims(1,2)
                    if( self%bandmsk(h,k) )then
                        self%cmat(h,k) = self%cmat(h,k) * dcmplx(ck*ch(h)-sk*sh(h), (sk*ch(h)+ck*sh(h)))
                    endif
                end do
            end do
        endif
    end subroutine shift

    pure subroutine normalize_mat( self )
        class(ft_expanded_dp), intent(inout) :: self
        real(dp) :: normsq
        normsq = self%corr_unnorm(self)
        if( normsq > DTINY )then
            self%cmat = merge(self%cmat/dsqrt(normsq), dcmplx(0.d0,0.d0), self%bandmsk)
        else
            self%cmat = dcmplx(0.d0,0.d0)
        endif
    end subroutine normalize_mat

    ! CALCULATORS

    pure real(dp) function corr( self1, self2 )
        class(ft_expanded_dp), intent(in) :: self1, self2
        real(dp) :: sumsq1, sumsq2
        sumsq1 = self1%corr_unnorm(self1)
        sumsq2 = self2%corr_unnorm(self2)
        if( sumsq1 * sumsq2 > DTINY )then
            corr = self1%corr_unnorm(self2) / dsqrt(sumsq1 * sumsq2)
        else
            corr = 0.d0
        endif
    end function corr

    ! unnormalized correlations, i.e. only the numerator term
    pure function corr_unnorm( self1, self2 ) result( r )
        class(ft_expanded_dp), intent(in) :: self1, self2
        real(dp) :: r
        r = sum(real(self1%cmat(0,:)*conjg(self2%cmat(0,:))), mask=self1%bandmsk(0,:))
        r = r + 2.d0*sum(real(self1%cmat(1:self1%flims(1,2),:)*conjg(self2%cmat(1:self1%flims(1,2),:))),&
            &mask=self1%bandmsk(1:self1%flims(1,2),:) )
    end function corr_unnorm

    subroutine gen_grad_noshift( self, self_gradx, self_grady )
        class(ft_expanded_dp), intent(in)    :: self
        class(ft_expanded_dp), intent(inout) :: self_gradx, self_grady
        real(dp) :: transfh, transfk
        integer  :: h,k
        do k = self%flims(2,1),self%flims(2,2)
            transfk = real(k,dp) * self%shconst(2)
            do h = self%flims(1,1),self%flims(1,2)
                if( self%bandmsk(h,k) )then
                    transfh = real(h,dp) * self%shconst(1)
                    self_gradx%cmat(h,k) = self%cmat(h,k) * JJ * transfh
                    self_grady%cmat(h,k) = self%cmat(h,k) * JJ * transfk
                end if
            end do
        end do
    end subroutine gen_grad_noshift

    ! gradients of correlations with respect to shifts, self is assumed pre-normalized to 1.
    ! selfR is destroyed on exit
    real(dp) function f( selfR, self, shift)
        class(ft_expanded_dp), intent(inout) :: selfR
        class(ft_expanded_dp), intent(in)    :: self
        real(dp),            intent(in)   :: shift(2)
        complex(dp) :: fcomp, fcompI
        real(dp)    ::  mag, denom, magR0, magR, f0
        integer     :: h,k
        real(dp) :: ch(self%flims(1,1):self%flims(1,2)),sh(self%flims(1,1):self%flims(1,2))
        real(dp) :: argh,argk,ck,sk
        magR0 = 0.d0
        magR  = 0.d0
        f0    = 0.d0
        f     = 0.d0
        do h = self%flims(1,1),self%flims(1,2)
            argh  = real(h,dp) * shift(1) * self%shconst(1)
            ch(h) = cos(argh)
            sh(h) = sin(argh)
        enddo
        do k = self%flims(2,1),self%flims(2,2)
            argk  = real(k,dp) * shift(2) * self%shconst(2)
            ck    = cos(argk)
            sk    = sin(argk)
            do h = self%flims(1,1),self%flims(1,2)
                if( self%bandmsk(h,k) )then
                    fcompI = self%cmat(h,k)
                    fcomp  = selfR%cmat(h,k) - fcompI
                    mag    = csq_fast(fcomp)
                    fcomp  = (fcomp * conjg(fcompI)) * dcmplx(ck*ch(h)-sk*sh(h), (sk*ch(h)+ck*sh(h)))
                    if( h == 0 )then
                        f0    = f0    + real(fcomp)
                        magR0 = magR0 + mag
                    else
                        f       = f    + real(fcomp)
                        magR    = magR + mag
                    endif
                endif
            end do
        enddo
        denom = (magR0 + 2.d0*magR)
        if( denom > 0.d0 )then
            denom = dsqrt(denom)
            f = (f0 + 2.d0*f) / denom
        else
            f = 0.d0
        endif
    end function  f

    ! gradients of correlations with respect to shifts, self is assumed pre-normalized to 1.
    ! selfR is destroyed on exit
    subroutine fdf( selfR, self, shift, f, g )
        class(ft_expanded_dp), intent(inout) :: selfR
        class(ft_expanded_dp), intent(in)    :: self
        real(dp),            intent(in)   :: shift(2)
        real(dp),            intent(out)   :: f, g(2)
        complex(dp) :: fcomp, fcompI
        real(dp)    :: g0(2), mag, denom, magR0, magR, f0, transfh, transfk
        integer     :: h,k
        real(dp) :: ch(self%flims(1,1):self%flims(1,2)),sh(self%flims(1,1):self%flims(1,2))
        real(dp) :: argh,argk,ck,sk
        g0    = 0.d0
        g     = 0.d0
        magR0 = 0.d0
        magR  = 0.d0
        f0    = 0.d0
        f     = 0.d0
        do h = self%flims(1,1),self%flims(1,2)
            argh  = real(h,dp) * shift(1) * self%shconst(1)
            ch(h) = cos(argh)
            sh(h) = sin(argh)
        enddo
        do k = self%flims(2,1),self%flims(2,2)
            transfk = real(k,dp) * self%shconst(2)
            argk  = real(k,dp) * shift(2) * self%shconst(2)
            ck    = cos(argk)
            sk    = sin(argk)
            do h = self%flims(1,1),self%flims(1,2)
                if( self%bandmsk(h,k) )then
                    fcompI = self%cmat(h,k)
                    fcomp  = selfR%cmat(h,k) - fcompI
                    mag    = csq_fast(fcomp)
                    transfh = real(h,dp) * self%shconst(1)
                    fcomp  = (fcomp * conjg(fcompI)) * dcmplx(ck*ch(h)-sk*sh(h), (sk*ch(h)+ck*sh(h)))
                    if( h == 0 )then
                        f0    = f0    + real(fcomp)
                        g0    = g0    + dimag(fcomp) * [transfh, transfk]
                        magR0 = magR0 + mag
                    else
                        f       = f    + real(fcomp)
                        g       = g    + dimag(fcomp) * [transfh, transfk]
                        magR    = magR + mag
                    endif
                endif
            end do
        enddo
        denom = (magR0 + 2.d0*magR)
        if( denom > 0.d0 )then
            denom = dsqrt(denom)
            f = (f0 + 2.d0*f) / denom
            g = (g0 + 2.d0*g) / denom
        else
            f = 0.d0
            g = 0.d0
        endif
    end subroutine fdf

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(ft_expanded_dp), intent(inout) :: self
        if( self%existence )then
            deallocate(self%cmat,self%bandmsk)
            self%existence = .false.
        endif
    end subroutine kill

end module simple_ft_expanded_dp
