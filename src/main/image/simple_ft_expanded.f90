!@descr: expanded Fourier transform class for improved cache utilisation
module simple_ft_expanded
use simple_core_module_api
use simple_image, only: image
implicit none
private
public :: ft_expanded
public :: ftexp_transfmat, ftexp_transfmat_init, ftexp_transfmat_kill
#include "simple_local_flags.inc"

real(dp), parameter :: num   = 1.0d8            ! numerator for rescaling of cost function
complex,  parameter :: JJ    = complex(0., 1.)

real,   allocatable :: ftexp_transfmat(:,:,:)   ! transfer matrix
integer             :: ftexp_transf_kzero       ! index shift for transfer & b-factor matrix

type :: ft_expanded
    private
    complex,     allocatable :: cmat(:,:)       !< Fourier components
    logical,     allocatable :: bandmsk(:,:)    !< for band-passed correlation
    complex(dp), allocatable :: tmp_cmat12(:,:) !< temporary matrix for shift search
    real                 :: hp                  !< high-pass limit
    real                 :: lp                  !< low-pass limit
    real                 :: smpd  = 0.          !< sampling distance of originating image
    real                 :: sumsq = 0.          !< sum of squares
    integer              :: lims(3,2)           !< physical limits for the Fourier transform
    integer              :: flims(3,2)          !< shifted limits
    integer              :: ldim(3)=[1,1,1]     !< logical dimension of originating image
    integer              :: kzero               !< to determine shift index to be applied to access the transfer matrix
    logical              :: existence=.false.   !< existence
  contains
    ! constructors
    procedure          :: new
    ! getters
    procedure          :: get_flims
    procedure          :: get_sumsq
    ! setters
    procedure          :: set_cmat
    procedure          :: zero
    procedure          :: copy
    ! arithmetics
    procedure          :: add
    procedure          :: add_uncond
    procedure          :: add2cmat
    procedure          :: subtr
    procedure          :: div
    ! modifiers
    procedure          :: shift
    procedure          :: shift_and_add
    procedure          :: normalize_mat
    ! calculators
    procedure          :: corr
    procedure          :: corr_unnorm
    procedure          :: corr_unnorm_serial
    procedure, private :: corr_normalize_sp
    procedure, private :: corr_normalize_dp
    generic            :: corr_normalize => corr_normalize_sp, corr_normalize_dp
    procedure          :: get_hp_lp
    ! shift optimization
    procedure          :: alloc_and_calc_tmp_cmat12, dealloc_tmp_cmat12
    procedure          :: corr_shifted_cost_8, corr_gshifted_cost_8, corr_fdfshifted_cost_8
    ! destructor
    procedure          :: kill
end type ft_expanded

contains

    ! CONSTRUCTORS

    subroutine new( self, img, hp, lp, fetch_comps, bfac )
        class(ft_expanded), intent(inout) :: self
        class(image),       intent(inout) :: img
        real,               intent(in)    :: hp, lp
        logical,            intent(in)    :: fetch_comps
        real, optional,     intent(in)    :: bfac
        real, allocatable :: wh(:)
        real    :: lp_nyq, spafreqh,spafreqk,wk
        integer :: h,k,i,hcnt,kcnt,lplim,hplim,sqarg
        logical :: didft, l_bfac
        ! kill pre-existing object
        call self%kill
        ! set constants
        self%ldim = img%get_ldim()
        if( self%ldim(3) > 1 ) THROW_HARD('only for 2D images; new_1')
        self%smpd = img%get_smpd()
        self%hp   = hp
        lp_nyq    = 2.*self%smpd
        self%lp   = max(lp, lp_nyq)
        self%lims = img%loop_lims(1,self%lp)
        ! shift the limits 2 make transfer 2 GPU painless
        self%flims = 1
        do i=1,3
            self%flims(i,2)= self%lims(i,2) - self%lims(i,1) + 1
        end do
        ! set the squared filter limits
        hplim = img%get_find(hp)
        hplim = hplim*hplim
        lplim = img%get_find(self%lp)
        lplim = lplim*lplim
        ! indexing to transfer matrix
        kcnt = 0
        do k=self%lims(2,1),self%lims(2,2)
            kcnt = kcnt + 1
            if( k==0 ) self%kzero = kcnt
        enddo
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
                 &self%bandmsk(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2)))
        self%cmat    = cmplx(0.,0.)
        self%bandmsk = .false.
        ! init matrices
        if( fetch_comps .and. l_bfac )then
            allocate(wh(self%lims(1,1):self%lims(1,2)))
            do h=self%lims(1,1),self%lims(1,2)
                spafreqh = real(h)/real(self%ldim(1))/self%smpd
                wh(h)    =  exp(-spafreqh*spafreqh*bfac/4.)
            enddo
        endif
        kcnt = 0
        do k=self%lims(2,1),self%lims(2,2)
            kcnt = kcnt + 1
            hcnt = 0
            if( fetch_comps .and. l_bfac )then
                spafreqk = real(k)/real(self%ldim(2))/self%smpd
                wk = exp(-spafreqk*spafreqk*bfac/4.)
            endif
            do h=self%lims(1,1),self%lims(1,2)
                sqarg = h*h + k*k
                hcnt  = hcnt + 1
                if( sqarg < 1 )then
                    cycle ! excludes zero
                elseif( sqarg <= lplim .and. sqarg >= hplim  )then
                    if( fetch_comps )then
                        if( l_bfac )then
                            self%cmat(hcnt,kcnt) = wh(h) * wk * img%get_fcomp2D(h,k)
                        else
                            self%cmat(hcnt,kcnt) = img%get_fcomp2D(h,k)
                        endif
                    endif
                    self%bandmsk(hcnt,kcnt) = .true.
                end if
            end do
        end do
        if( fetch_comps ) self%sumsq = real(self%corr_unnorm_serial(self))
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

    pure real function get_sumsq( self )
        class(ft_expanded), intent(in) :: self
        get_sumsq = self%sumsq
    end function get_sumsq

    ! SETTERS

    subroutine set_cmat( self, cmat )
        class(ft_expanded), intent(inout) :: self
        complex,            intent(in) :: cmat(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2))
        self%cmat = cmat
        self%sumsq = real(self%corr_unnorm_serial(self))
    end subroutine set_cmat

    pure subroutine zero( self )
        class(ft_expanded), intent(inout) :: self
        self%cmat  = cmplx(0.,0.)
        self%sumsq = 0.
    end subroutine zero

    ! not a constructor
    pure subroutine copy( self, self2copy )
        class(ft_expanded), intent(inout) :: self, self2copy
        self%hp    = self2copy%hp
        self%lp    = self2copy%lp
        self%smpd  = self2copy%smpd
        self%ldim  = self2copy%ldim
        self%lims  = self2copy%lims
        self%flims = self2copy%flims
        self%kzero = self2copy%kzero
        self%sumsq = self2copy%sumsq
        self%cmat  = self2copy%cmat
        self%bandmsk = self2copy%bandmsk
        self%existence = self%existence
    end subroutine copy

    ! ARITHMETICS

    subroutine add( self, self2add, w )
        class(ft_expanded), intent(inout) :: self
        class(ft_expanded), intent(in)    :: self2add
        real, optional,     intent(in)    :: w
        if( present(w) )then
             if( w > 0. )then
                 self%cmat  = self%cmat + self2add%cmat*w
                 self%sumsq = real(self%corr_unnorm_serial(self))
             endif
        else
            self%cmat  = self%cmat + self2add%cmat
            self%sumsq = real(self%corr_unnorm_serial(self))
        endif
    end subroutine add

    subroutine add_uncond( self, self2add, w )
        class(ft_expanded), intent(inout) :: self
        class(ft_expanded), intent(in)    :: self2add
        real, optional,     intent(in)    :: w
        if( present(w) )then
             self%cmat = self%cmat + self2add%cmat*w
        else
            self%cmat = self%cmat + self2add%cmat
        endif
    end subroutine add_uncond

    subroutine add2cmat( self, cmat, w )
        class(ft_expanded), intent(in)    :: self
        complex,            intent(inout) :: cmat(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2))
        real,               intent(in)    :: w
        !$omp parallel workshare proc_bind(close)
        cmat = cmat + w * self%cmat
        !$omp end parallel workshare
    end subroutine add2cmat

    subroutine subtr( self, self2subtr, w )
        class(ft_expanded), intent(inout) :: self
        class(ft_expanded), intent(in)    :: self2subtr
        real, optional,     intent(in)    :: w
        real :: ww
        ww = 1.0
        if( present(w) ) ww = w
        if( ww > 0. ) self%cmat = self%cmat-ww*self2subtr%cmat
        self%sumsq = real(self%corr_unnorm_serial(self))
    end subroutine subtr

    subroutine div( self, val )
        class(ft_expanded), intent(inout) :: self
        real,               intent(in)    :: val
        if( abs(val) < TINY ) return
        self%cmat  = self%cmat / val
        self%sumsq = real(self%corr_unnorm_serial(self))
    end subroutine div

    ! MODIFIERS

    subroutine shift( self, shvec, self_out, calc_sumsq )
        class(ft_expanded),           intent(inout) :: self
        real,                         intent(in)    :: shvec(2)
        class(ft_expanded), optional, intent(inout) :: self_out
        logical,            optional, intent(in)    :: calc_sumsq
        integer  :: hind,kind,kind_shift,kkind
        real(dp) :: ch(self%flims(1,1):self%flims(1,2)),sh(self%flims(1,1):self%flims(1,2))
        real(dp) :: argh,argk,ck,sk
        kind_shift = ftexp_transf_kzero - self%kzero
        do hind=self%flims(1,1),self%flims(1,2)
            argh     = real(ftexp_transfmat(hind,1,1) * shvec(1),dp)
            ch(hind) = cos(argh)
            sh(hind) = sin(argh)
        enddo
        if( present(self_out) )then
            do kind=self%flims(2,1),self%flims(2,2)
                kkind = kind+kind_shift
                argk  = real(ftexp_transfmat(1,kkind,2) * shvec(2),dp)
                ck    = cos(argk)
                sk    = sin(argk)
                do hind=self%flims(1,1),self%flims(1,2)
                    if( self%bandmsk(hind,kind) )then
                        self_out%cmat(hind,kind) = self%cmat(hind,kind) * cmplx(ck*ch(hind)-sk*sh(hind), (sk*ch(hind)+ck*sh(hind)))
                    endif
                end do
            end do
        else
            do kind=self%flims(2,1),self%flims(2,2)
                kkind = kind+kind_shift
                argk  = ftexp_transfmat(1,kkind,2) * shvec(2)
                ck    = cos(argk)
                sk    = sin(argk)
                do hind=self%flims(1,1),self%flims(1,2)
                    if( self%bandmsk(hind,kind) )then
                        self%cmat(hind,kind) = self%cmat(hind,kind) * cmplx(ck*ch(hind)-sk*sh(hind), (sk*ch(hind)+ck*sh(hind)))
                    endif
                end do
            end do
        endif
        if( present(calc_sumsq) )then
            if( calc_sumsq ) self_out%sumsq = real(self_out%corr_unnorm_serial(self))
        else
            self_out%sumsq = real(self_out%corr_unnorm_serial(self_out))
        endif
    end subroutine shift

    subroutine shift_and_add( self, shvec, w, self_out )
        class(ft_expanded), intent(in)    :: self
        real,               intent(in)    :: shvec(2), w
        class(ft_expanded), intent(inout) :: self_out
        integer :: hind,kind,kind_shift,kkind
        real    :: ch(self%flims(1,1):self%flims(1,2)),sh(self%flims(1,1):self%flims(1,2))
        real    :: argh,argk,ck,sk
        kind_shift = ftexp_transf_kzero - self%kzero
        do hind=self%flims(1,1),self%flims(1,2)
            argh     = ftexp_transfmat(hind,1,1) * shvec(1)
            ch(hind) = cos(argh)
            sh(hind) = sin(argh)
        enddo
        do kind=self%flims(2,1),self%flims(2,2)
            kkind = kind+kind_shift
            argk  = ftexp_transfmat(1,kkind,2) * shvec(2)
            ck    = cos(argk)
            sk    = sin(argk)
            do hind=self%flims(1,1),self%flims(1,2)
                if( self%bandmsk(hind,kind) )then
                    self_out%cmat(hind,kind) = self_out%cmat(hind,kind) + w * self%cmat(hind,kind) * cmplx(ck*ch(hind)-sk*sh(hind), (sk*ch(hind)+ck*sh(hind)))
                endif
            end do
        end do
        self_out%sumsq = real(self_out%corr_unnorm_serial(self_out))
    end subroutine shift_and_add

    subroutine normalize_mat( self )
        class(ft_expanded), intent(inout) :: self
        real    :: anorm
        anorm = real(sqrt(self%corr_unnorm_serial(self)))
        where(self%bandmsk)
            self%cmat = self%cmat / anorm
        end where
    end subroutine normalize_mat

    ! CALCULATORS

    function corr( self1, self2 ) result( r )
        class(ft_expanded), intent(inout) :: self1, self2
        real :: r
        ! corr is real part of the complex mult btw 1 and 2*
        r =        sum(real(self1%cmat(                 1,:)*conjg(self2%cmat(                 1,:))), mask=self1%bandmsk(1,:))
        r = r + 2.*sum(real(self1%cmat(2:self1%flims(1,2),:)*conjg(self2%cmat(2:self1%flims(1,2),:))), mask=self1%bandmsk(2:self1%flims(1,2),:) )
        ! normalise the correlation coefficient
        if( self1%sumsq > 0. .and. self2%sumsq > 0. )then
            r = r / sqrt(self1%sumsq * self2%sumsq)
        else
            r = 0.
        endif
    end function corr

    ! unnormalized correlations, i.e. only the numerator term
    function corr_unnorm( self1, self2 ) result( r )
        class(ft_expanded), intent(in) :: self1, self2
        real(dp) :: r, tmp ! we need double precision here to avoid round-off errors in openmp loop
        integer  :: hind, kind
        ! corr is real part of the complex mult btw 1 and 2*
        r = sum(real(self1%cmat(1,:)*conjg(self2%cmat(1,:)),dp), mask=self1%bandmsk(1,:))
        tmp = 0._dp
        !$omp parallel do collapse(2) default(shared) private(kind,hind) reduction(+:tmp) proc_bind(close) schedule(static)
        do kind = self1%flims(2,1), self1%flims(2,2)
            do hind = 2, self1%flims(1,2)
                if (self1%bandmsk(hind,kind)) then
                    tmp = tmp + real(self1%cmat(hind,kind)*conjg(self2%cmat(hind,kind)),dp)
                end if
            end do
        end do
        !$omp end parallel do
        r = r + 2._dp * tmp
    end function corr_unnorm

    ! unnormalized correlations, i.e. only the numerator term
    pure function corr_unnorm_serial( self1, self2 ) result( r )
        class(ft_expanded), intent(in) :: self1, self2
        real(dp) :: r
        r = sum(real(self1%cmat(1,:)*conjg(self2%cmat(1,:))), mask=self1%bandmsk(1,:))
        r = r + 2.d0*sum(real(self1%cmat(2:self1%flims(1,2),:)*conjg(self2%cmat(2:self1%flims(1,2),:))), mask=self1%bandmsk(2:self1%flims(1,2),:) )
    end function corr_unnorm_serial

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

    subroutine get_hp_lp( self, hp, lp )
        class(ft_expanded), intent(in) :: self
        real, intent(out) :: hp, lp
        hp = self%hp
        lp = self%lp
    end subroutine get_hp_lp

    ! FOR SHIFT OPTIMIZATION

    !< allocate & calculate tmp matrix for cost function
    subroutine alloc_and_calc_tmp_cmat12( self, self_reference, denominator )
        class(ft_expanded), intent(inout) :: self           !< is the particle instance
        class(ft_expanded), intent(in)    :: self_reference
        real(dp),           intent(out)   :: denominator
        logical :: do_alloc
        do_alloc = .true.
        if ( allocated( self%tmp_cmat12 ) ) then
            if ( (ubound(self%tmp_cmat12, 1) == self%flims(1,2)) .and. &
                 (ubound(self%tmp_cmat12, 2) == self%flims(2,2)) ) then
                do_alloc = .false.
            else
                deallocate( self%tmp_cmat12 )
            end if
        end if
        if( do_alloc ) allocate(self%tmp_cmat12(1:self%flims(1,2),1:self%flims(2,2)))
        ! fill matrix
        denominator = dsqrt(real(self_reference%get_sumsq(),dp) * real(self%get_sumsq(),dp))
        self%tmp_cmat12 = merge(self_reference%cmat * conjg(self%cmat(:,:)), cmplx(0.,0.), self%bandmsk)
    end subroutine alloc_and_calc_tmp_cmat12

    subroutine dealloc_tmp_cmat12( self )
        class(ft_expanded), intent(inout) :: self
        if(allocated(self%tmp_cmat12)) deallocate(self%tmp_cmat12)
    end subroutine dealloc_tmp_cmat12

    !< cost function for minimizer, cost only
    pure function corr_shifted_cost_8( self, shvec, denominator )result( r )
        class(ft_expanded), intent(in) :: self
        real(dp),           intent(in) :: shvec(2), denominator
        real(dp) :: ch(self%flims(1,1):self%flims(1,2)),sh(self%flims(1,1):self%flims(1,2))
        real(dp) :: r, r1, r2, ck,sk, argh,argk
        integer  :: hind,kind,kkind,kind_shift
        do hind=self%flims(1,1),self%flims(1,2)
            argh     = real(ftexp_transfmat(hind,1,1),dp) * shvec(1)
            ch(hind) = cos(argh)
            sh(hind) = sin(argh)
        enddo
        r1 = 0.d0
        r2 = 0.d0
        kind_shift = ftexp_transf_kzero - self%kzero
        do kind=self%flims(2,1),self%flims(2,2)
            kkind = kind+kind_shift
            argk  = real(ftexp_transfmat(1,kkind,2),dp) * shvec(2)
            ck    = cos(argk)
            sk    = sin(argk)
            do hind=self%flims(1,1),self%flims(1,2)
                if( self%bandmsk(hind,kind) )then
                    if( hind == 1 )then
                        ! h = 0
                        r1  = r1 + real(self%tmp_cmat12(1,kind) * dcmplx(ck*ch(hind)-sk*sh(hind), -(sk*ch(hind)+ck*sh(hind))),kind=dp)
                    else
                        ! h > 0
                        r2  = r2 + real(self%tmp_cmat12(hind,kind) * dcmplx(ck*ch(hind)-sk*sh(hind), -(sk*ch(hind)+ck*sh(hind))),kind=dp)
                    endif
                endif
            end do
        enddo
        r = (r1 + 2.d0*r2) * num / denominator
    end function corr_shifted_cost_8

    !< cost function for minimizer, gradient only
    subroutine corr_gshifted_cost_8( self, shvec, denominator, grad )
        class(ft_expanded), intent(inout) :: self
        real(dp),           intent(in)    :: shvec(2), denominator
        real(dp),           intent(out)   :: grad(2)
        real(dp)    :: ch(self%flims(1,1):self%flims(1,2)),sh(self%flims(1,1):self%flims(1,2))
        real(dp)    :: g1(2),g2(2),transf_vec(2), ck,sk, argh,argk
        integer     :: hind,kind,kkind,kind_shift
        do hind = self%flims(1,1),self%flims(1,2)
            argh     = real(ftexp_transfmat(hind,1,1),dp) * shvec(1)
            ch(hind) = cos(argh)
            sh(hind) = sin(argh)
        enddo
        g1 = 0.d0
        g2 = 0.d0
        kind_shift = ftexp_transf_kzero - self%kzero
        do kind = self%flims(2,1),self%flims(2,2)
            kkind = kind + kind_shift
            argk  = real(ftexp_transfmat(1,kkind,2),dp) * shvec(2)
            ck    = cos(argk)
            sk    = sin(argk)
            do hind = self%flims(1,1),self%flims(1,2)
                if( self%bandmsk(hind,kind) )then
                    transf_vec = real(ftexp_transfmat(hind,kkind,:),dp)
                    if( hind == 1 )then ! h = 0
                        g1(:) = g1(:) + dimag(self%tmp_cmat12(hind,kind) * dcmplx(ck*ch(hind)-sk*sh(hind), -(sk*ch(hind)+ck*sh(hind))))*transf_vec
                    else ! h > 0
                        g2(:) = g2(:) + dimag(self%tmp_cmat12(hind,kind) * dcmplx(ck*ch(hind)-sk*sh(hind), -(sk*ch(hind)+ck*sh(hind))))*transf_vec
                    endif
                endif
            end do
        enddo
        grad(1) = (g1(1)+ 2.d0*g2(1)) * num / denominator
        grad(2) = (g1(2)+ 2.d0*g2(2)) * num / denominator
    end subroutine corr_gshifted_cost_8

    !< cost function for minimizer, f and gradient
    subroutine corr_fdfshifted_cost_8( self, shvec, denominator, f, grad )
        class(ft_expanded), intent(inout) :: self
        real(dp),            intent(in)    :: shvec(2), denominator
        real(dp),            intent(out)   :: grad(2), f
        complex(dp) :: tmp
        real(dp)    :: f1,f2,g1(2),g2(2), transf_vec(2), argh,argk, ck,sk
        real(dp)    :: ch(self%flims(1,1):self%flims(1,2)),sh(self%flims(1,1):self%flims(1,2))
        integer     :: hind,kind,kkind,kind_shift
        f1 = 0.d0
        f2 = 0.d0
        g1 = 0.d0
        g2 = 0.d0
        kind_shift = ftexp_transf_kzero - self%kzero
        do hind = self%flims(1,1),self%flims(1,2)
            argh     = real(ftexp_transfmat(hind,1,1),dp) * shvec(1)
            ch(hind) = cos(argh)
            sh(hind) = sin(argh)
        enddo
        do kind = self%flims(2,1),self%flims(2,2)
            kkind = kind + kind_shift
            argk  = real(ftexp_transfmat(1,kkind,2),dp) * shvec(2)
            ck    = cos(argk)
            sk    = sin(argk)
            do hind=self%flims(1,1),self%flims(1,2)
                if( self%bandmsk(hind,kind) )then
                    transf_vec = real(ftexp_transfmat(hind,kkind,:),dp)
                    tmp        = self%tmp_cmat12(hind,kind) * dcmplx(ck*ch(hind)-sk*sh(hind), -(sk*ch(hind)+ck*sh(hind)))
                    if( hind == 1 )then ! h = 0
                        f1    = f1    + real(tmp,dp)
                        g1(:) = g1(:) + dimag(tmp) * transf_vec
                    else ! h > 0
                        f2    = f2    + real(tmp,dp)
                        g2(:) = g2(:) + dimag(tmp) * transf_vec
                    endif
                endif
            end do
        enddo
        f       = (f1   + 2.d0*f2)    * num / denominator
        grad(1) = (g1(1)+ 2.d0*g2(1)) * num / denominator
        grad(2) = (g1(2)+ 2.d0*g2(2)) * num / denominator
    end subroutine corr_fdfshifted_cost_8

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(ft_expanded), intent(inout) :: self
        if( self%existence )then
            deallocate(self%cmat,self%bandmsk)
            call self%dealloc_tmp_cmat12
            self%existence = .false.
        endif
    end subroutine kill

    ! TRANSFER MATRIX ROUTINES

    subroutine ftexp_transfmat_init( img, lp )
        class(image),   intent(in) :: img
        real, optional, intent(in) :: lp
        real    :: lp_here, shconst(2)
        integer :: h,k,i,hcnt,kcnt
        integer :: ldim(3),flims(3,2),ftexp_transf_flims(3,2)
        call ftexp_transfmat_kill
        ! dimensions
        ldim = img%get_ldim()
        if( ldim(3) > 1 ) THROW_HARD("In: ftexp_transfmat_init; simple_ft_expanded, 1")
        lp_here = 2.*img%get_smpd()
        if( present(lp) ) lp_here = lp
        flims = img%loop_lims(1,lp_here)
        ftexp_transf_flims = flims
        do i=1,3
            ftexp_transf_flims(i,2) = ftexp_transf_flims(i,2) - ftexp_transf_flims(i,1) + 1
        end do
        ftexp_transf_flims(:,1) = 1
        ! set shift constant
        shconst = 0.
        do i=1,2
            if( ldim(i) > 1 )then
                if( is_even(ldim(i)) )then
                    shconst(i) = PI/real(ldim(i)/2)
                else
                    shconst(i) = PI/real((ldim(i)-1)/2)
                endif
            endif
        end do
        allocate(ftexp_transfmat(ftexp_transf_flims(1,1):ftexp_transf_flims(1,2),&
                                &ftexp_transf_flims(2,1):ftexp_transf_flims(2,2), 2), source=0.)
        !$omp parallel do collapse(2) default(shared) private(h,k,hcnt,kcnt) proc_bind(close) schedule(static)
        do k=flims(2,1),flims(2,2)
            do h=flims(1,1),flims(1,2)
                hcnt = h-flims(1,1)+1
                kcnt = k-flims(2,1)+1
                if( k==0 ) ftexp_transf_kzero = kcnt ! for transfer matrix indexing
                ftexp_transfmat(hcnt,kcnt,:)  = real([h,k])*shconst
            end do
        end do
        !$omp end parallel do
    end subroutine ftexp_transfmat_init

    subroutine ftexp_transfmat_kill
        if( allocated(ftexp_transfmat) ) deallocate(ftexp_transfmat)
        ftexp_transf_kzero = 0
    end subroutine ftexp_transfmat_kill

end module simple_ft_expanded
