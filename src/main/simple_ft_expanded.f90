! expanded Fourier transform class for improved cache utilisation
module simple_ft_expanded
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,  only: image
implicit none
private
public :: ft_expanded, ftexp_gen_bfacweights, ftexp_get_bfacweights_ptr, ftexp_has_bfacweights
#include "simple_local_flags.inc"

real(sp), target, allocatable :: bfac_weights(:,:,:) !< 2D b-factor weights matrix
real(dp),         parameter   :: denom = 0.00075_dp  ! denominator for rescaling of cost function
logical                       :: l_bfac = .false.

type :: ft_expanded
    private
    integer              :: lims(3,2)          !< physical limits for the Fourier transform
    integer              :: flims(3,2)         !< shifted limits (2 make transfer 2 GPU painless)
    integer              :: ldim(3)=[1,1,1]    !< logical dimension of originating image
    integer              :: flims_nyq(3,2)     !< nyquist limits
    real                 :: shconst(3)         !< shift constant
    real                 :: hp                 !< high-pass limit
    real                 :: lp                 !< low-pass limit
    real                 :: smpd=0.            !< sampling distance of originating image
    real,    allocatable :: transfmat(:,:,:,:) !< shift transfer matrix (TODO: make this global to save memory)
    complex, allocatable :: cmat(:,:,:)        !< Fourier components
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
    ! setters
    procedure          :: set_cmat
    procedure          :: extract_img
    procedure          :: zero
    ! arithmetics
    procedure          :: add
    procedure          :: subtr
    ! modifiers
    procedure          :: shift
    ! calculators
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
                    if( sqarg <= lplim .and. sqarg >= hplim  )then
                        phys = img%comp_addr_phys([h,k,l])
                        self%transfmat(hcnt,kcnt,lcnt,:) = real([h,k,l])*self%shconst
                        if( fetch_comps ) self%cmat(hcnt,kcnt,lcnt) = img%get_fcomp([h,k,l],phys)
                    end if
                end do
            end do
        end do
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

    ! SETTERS

    pure subroutine set_cmat( self, cmat )
        class(ft_expanded), intent(inout) :: self
        complex,            intent(in) :: cmat(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2),self%flims(3,1):self%flims(3,2))
        self%cmat = cmat
    end subroutine set_cmat

    subroutine extract_img( self, img, hp, lp )
        class(ft_expanded), intent(inout) :: self
        class(image),       intent(inout) :: img
        real,               intent(in)    :: hp, lp
        logical :: didft
        integer :: hcnt, h, hh, kcnt, k, kk, lcnt, l, ll, sqarg
        integer :: hplim, lplim
        integer :: phys(3)
        if (do_call_new()) then
            call self%new(img, hp, lp, .true.)
            return
        end if
        hplim = img%get_find(hp)
        hplim = hplim*hplim
        lplim = img%get_find(lp)
        lplim = lplim*lplim
        ! prepare image
        didft = .false.
        if( .not. img%is_ft() )then
            call img%fft()
            didft = .true.
        endif
        self%cmat = cmplx(0.,0.)
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
                        self%cmat(hcnt,kcnt,lcnt) = img%get_fcomp([h,k,l],phys)
                    end if
                end do
            end do
        end do
        if( didft ) call img%ifft()

    contains

        function do_call_new() result( res )
            logical :: res
            integer :: ldim_here(3)
            real    :: smpd_here
            real    :: lp_nyq
            real    :: lp_here
            integer :: lims_here(3,2)
            integer :: flims_nyq_here(3,2)
            integer :: flims_here(3,2)
            integer :: i
            res = .true.
            if (.not. self%existence)            return
            if (.not. allocated(self%cmat     )) return
            if (.not. allocated(self%transfmat)) return
            ldim_here      = img%get_ldim()
            if( ldim_here(3) > 1 ) THROW_HARD('only for 2D images; extract_img')
            if (any(ldim_here(1:2) /= self%ldim(1:2))) return
            smpd_here      = img%get_smpd()
            if (smpd_here /= self%smpd) return
            if (self%hp /= hp) return
            lp_nyq         = 2.*self%smpd
            lp_here        = max(lp, lp_nyq)
            if (lp_here /= self%lp) return
            lims_here      = img%loop_lims(1,self%lp)
            if (any(lims_here(:,:) /= self%lims(:,:))) return
            flims_nyq_here = img%loop_lims(1,lp_nyq)
            flims_here     = 1
            do i=1,3
                flims_here(i,2)     = lims_here(i,2)      - lims_here(i,1)      + 1
                flims_nyq_here(i,2) = flims_nyq_here(i,2) - flims_nyq_here(i,1) + 1
            end do
            if (any(flims_nyq_here(:,:) /= self%flims_nyq(:,:))) return
            if (any(flims_here(:,:)     /= self%flims(:,:)    )) return
            res = .false.
        end function do_call_new

    end subroutine extract_img

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
        if( l_bfac )then
            ! b-factor weighted correlation
            r = sum(bfac_weights(1,1:self1%flims(2,2)-1,1)*&
            &    real(self1%cmat(1,1:self1%flims(2,2)-1,1)*&
            &   conjg(self2%cmat(1,1:self1%flims(2,2)-1,1))) )
            r = r + sum( bfac_weights(self1%flims(1,2),1:self1%flims(2,2)-1,1)*&
            &         real(self1%cmat(self1%flims(1,2),1:self1%flims(2,2)-1,1)*&
            &        conjg(self2%cmat(self1%flims(1,2),1:self1%flims(2,2)-1,1))) )
            r = r + 2.*sum( bfac_weights(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)*&
            &            real(self1%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)*&
            &           conjg(self2%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1))) )
            sumasq =             sum(bfac_weights(                   1,1:self1%flims(2,2)-1,1)*csq(self1%cmat(                   1,1:self1%flims(2,2)-1,1)))
            sumasq = sumasq +    sum(bfac_weights(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)*csq(self1%cmat(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)))
            sumasq = sumasq + 2.*sum(bfac_weights(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)*csq(self1%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)))
            sumbsq =             sum(bfac_weights(                   1,1:self1%flims(2,2)-1,1)*csq(self2%cmat(                   1,1:self1%flims(2,2)-1,1)))
            sumbsq = sumbsq +    sum(bfac_weights(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)*csq(self2%cmat(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)))
            sumbsq = sumbsq + 2.*sum(bfac_weights(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)*csq(self2%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)))
        else
            ! corr is real part of the complex mult btw 1 and 2*
            r =        sum(real(self1%cmat(                   1,1:self1%flims(2,2)-1,1)*&
                conjg(self2%cmat(                   1,1:self1%flims(2,2)-1,1))))
            r = r +    sum(real(self1%cmat(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)*&
                conjg(self2%cmat(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1))))
            r = r + 2.*sum(real(self1%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)*&
                conjg(self2%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1))))
            ! normalisation terms
            sumasq =             sum(csq(self1%cmat(                   1,1:self1%flims(2,2)-1,1)))
            sumasq = sumasq +    sum(csq(self1%cmat(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)))
            sumasq = sumasq + 2.*sum(csq(self1%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)))
            sumbsq =             sum(csq(self2%cmat(                   1,1:self1%flims(2,2)-1,1)))
            sumbsq = sumbsq +    sum(csq(self2%cmat(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)))
            sumbsq = sumbsq + 2.*sum(csq(self2%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)))
        endif
        ! finalise the correlation coefficient
        if( sumasq > 0. .and. sumbsq > 0. )then
            r = r / sqrt(sumasq * sumbsq)
        else
            r = 0.
        endif
    end function corr

    subroutine corr_normalize_sp( self1, self2, corr )
        class(ft_expanded), intent(inout) :: self1, self2
        real(sp),           intent(inout) :: corr
        real(dp) :: sumasq,sumbsq
        if( l_bfac )then
            sumasq =             sum(bfac_weights(                   1,1:self1%flims(2,2)-1,1)*csq(self1%cmat(                   1,1:self1%flims(2,2)-1,1)))
            sumasq = sumasq +    sum(bfac_weights(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)*csq(self1%cmat(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)))
            sumasq = sumasq + 2.*sum(bfac_weights(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)*csq(self1%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)))
            sumbsq =             sum(bfac_weights(                   1,1:self1%flims(2,2)-1,1)*csq(self2%cmat(                   1,1:self1%flims(2,2)-1,1)))
            sumbsq = sumbsq +    sum(bfac_weights(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)*csq(self2%cmat(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)))
            sumbsq = sumbsq + 2.*sum(bfac_weights(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)*csq(self2%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)))
        else
            sumasq =             sum(csq(self1%cmat(                   1,1:self1%flims(2,2)-1,1)))
            sumasq = sumasq +    sum(csq(self1%cmat(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)))
            sumasq = sumasq + 2.*sum(csq(self1%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)))
            sumbsq =             sum(csq(self2%cmat(                   1,1:self1%flims(2,2)-1,1)))
            sumbsq = sumbsq +    sum(csq(self2%cmat(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)))
            sumbsq = sumbsq + 2.*sum(csq(self2%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)))
        endif
        corr   = real(real(corr,dp) * denom / sqrt(sumasq * sumbsq), sp)
    end subroutine corr_normalize_sp

    subroutine corr_normalize_dp( self1, self2, corr )
        class(ft_expanded), intent(inout) :: self1, self2
        real(dp),           intent(inout) :: corr
        real(dp) :: sumasq,sumbsq
        if( l_bfac )then
            sumasq =             sum(bfac_weights(                   1,1:self1%flims(2,2)-1,1)*csq(self1%cmat(                   1,1:self1%flims(2,2)-1,1)))
            sumasq = sumasq +    sum(bfac_weights(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)*csq(self1%cmat(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)))
            sumasq = sumasq + 2.*sum(bfac_weights(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)*csq(self1%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)))
            sumbsq =             sum(bfac_weights(                   1,1:self1%flims(2,2)-1,1)*csq(self2%cmat(                   1,1:self1%flims(2,2)-1,1)))
            sumbsq = sumbsq +    sum(bfac_weights(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)*csq(self2%cmat(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)))
            sumbsq = sumbsq + 2.*sum(bfac_weights(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)*csq(self2%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)))
        else
            sumasq =             sum(csq(self1%cmat(                   1,1:self1%flims(2,2)-1,1)))
            sumasq = sumasq +    sum(csq(self1%cmat(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)))
            sumasq = sumasq + 2.*sum(csq(self1%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)))
            sumbsq =             sum(csq(self2%cmat(                   1,1:self1%flims(2,2)-1,1)))
            sumbsq = sumbsq +    sum(csq(self2%cmat(  self1%flims(1,2)  ,1:self1%flims(2,2)-1,1)))
            sumbsq = sumbsq + 2.*sum(csq(self2%cmat(2:self1%flims(1,2)-1,1:self1%flims(2,2)-1,1)))
        endif
        corr   = corr * denom / sqrt(sumasq * sumbsq)
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

    ! GLOBAL MODULE ROUTINES

    logical function ftexp_has_bfacweights()
        ftexp_has_bfacweights = l_bfac
    end function ftexp_has_bfacweights

    !>   \brief  Initializes B-factor weights
    subroutine ftexp_gen_bfacweights( img, bfac )
        class(image),       intent(inout) :: img
        real,               intent(in)    :: bfac
        real    :: g2, Ah, Ak, Al, B, smpd
        integer :: i,h,k,l, alloc_stat, ldim(3), flims_nyq(3,2)
        smpd = img%get_smpd()
        ldim = img%get_ldim()
        flims_nyq = img%loop_lims(1,2.*smpd)
        do i=1,3
            flims_nyq(i,2) = flims_nyq(i,2) - flims_nyq(i,1) + 1
        end do
        flims_nyq(:,1) = 1
        Ah = (real(ldim(1))*smpd)**2.
        Ak = (real(ldim(2))*smpd)**2.
        Al = (real(ldim(3))*smpd)**2.
        if( allocated(bfac_weights) ) deallocate(bfac_weights)
        allocate(bfac_weights(1:flims_nyq(1,2),1:flims_nyq(2,2),1:flims_nyq(3,2)),stat=alloc_stat)
        B = - bfac / 4.
        if( ldim(3) == 1 )then
            ! 2D
            do k=flims_nyq(2,1),flims_nyq(2,2)
                do h=flims_nyq(1,1),flims_nyq(1,2)
                    g2 = real(h*h)/Ah + real(k*k)/Ak
                    bfac_weights(h,k,1) = max(0., exp(B*g2))
                end do
            end do
        else
            ! 3D
            do l=flims_nyq(3,1),flims_nyq(3,2)
                do k=flims_nyq(2,1),flims_nyq(2,2)
                    do h=flims_nyq(1,1),flims_nyq(1,2)
                        g2 = real(h*h)/Ah + real(k*k)/Ak + real(l*l)/Al
                        bfac_weights(h,k,l) = max(0., exp(B*g2))
                    enddo
                end do
            end do
        endif
        l_bfac = .true.
    end subroutine ftexp_gen_bfacweights

    subroutine ftexp_get_bfacweights_ptr( bfacweights_ptr )
        real, pointer, intent(out) :: bfacweights_ptr(:,:,:)
        bfacweights_ptr => bfac_weights
    end subroutine ftexp_get_bfacweights_ptr

    !>  \brief  is a destructor
    subroutine ftexp_kill_bfac_weights
        if(allocated(bfac_weights))deallocate(bfac_weights)
        l_bfac    = .false.
    end subroutine ftexp_kill_bfac_weights

end module simple_ft_expanded
