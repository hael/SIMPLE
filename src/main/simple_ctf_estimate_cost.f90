module simple_ctf_estimate_cost
include 'simple_lib.f08'
use simple_image,     only: image
use simple_ctf,       only: ctf
use simple_opt_spec,  only: opt_spec
use simple_opt_de,    only: opt_de
use simple_optimizer, only: optimizer
implicit none

public :: ctf_estimate_cost1D, ctf_estimate_cost2D, ctf_estimate_cost4Dcont
private
#include "simple_local_flags.inc"

type ctf_estimate_cost1D
    private
    class(image),     pointer :: pspec              !< reference SQRT( power spectrum )
    class(ctfparams), pointer :: parms              !< For microscope characteristics
    real,             pointer :: spec1d(:)
    type(ctf)                 :: tfun               !< Transfer function
    integer                   :: ldim(3)      = 0
    integer                   :: reslims1d(2) = 0
    contains
        procedure, private :: init1D
        procedure, private :: cost1D
        procedure, private :: kill1D
        procedure          :: init => init1D
        procedure          :: cost => cost1D
        procedure          :: kill => kill1D
end type ctf_estimate_cost1D

type ctf_estimate_cost2D
    private
    class(image),     pointer :: pspec              !< reference SQRT( power spectrum )
    class(ctfparams), pointer :: parms              !< For microscope characteristics
    integer,          pointer :: inds(:,:)          !< Mask indices
    type(opt_spec)            :: ospec              !< optimizer specification
    type(opt_de)              :: diffevol           !< DE optimizer
    type(ctf)                 :: tfun               !< Transfer function
    real                      :: astigtol     = 0.05  !< Astigmatic tolerance
    integer                   :: ndim         = 0     !< # of dimensions to optimize
    integer                   :: ldim(3)      = 0
    integer                   :: flims(3,2)   = 0
    integer                   :: ninds        = 0
    contains
        procedure, private :: init2D
        procedure          :: minimize
        procedure, private :: cost2D
        procedure, private :: calc_cost2D
        procedure, private :: kill2D
        procedure          :: init => init2D
        procedure          :: kill => kill2D
end type ctf_estimate_cost2D

type ctf_estimate_cost4Dcont
    private
    class(image),     pointer :: pspec => null()    !< reference SQRT( power spectrum ) normalized
    class(ctfparams), pointer :: parms => null()    !< For microscope characteristics
    class(optimizer), pointer :: nlopt => null()    !< optimizer object
    integer,          pointer :: inds(:,:)          !< Mask indices
    type(opt_spec)            :: ospec              !< optimizer specification object
    real(dp)                  :: kv, cs, fraca, angast, df, phshift, smpd, amp_contr_const, wl, wlsq ! CTF
    real(dp)                  :: conversionfactor = 0.d0    ! unit conversion
    real(dp)                  :: astigtol         = 0.05d0  !< Astigmatic tolerance (microns)
    integer                   :: ndim             = 0       !< # of dimensions to optimize
    integer                   :: ldim(3)          = 0       !< dimensions of spectrum
    integer                   :: flims(3,2)       = 0       !< fourier limits
    integer                   :: ninds            = 0       !< number of pixels to correlate
    contains
        procedure, private :: init4Dcont
        procedure, private :: minimize4Dcont
        procedure, private :: fdf
        procedure, private :: kill4Dcont
        procedure          :: init     => init4Dcont
        procedure          :: minimize => minimize4Dcont
        procedure          :: kill     => kill4Dcont
end type ctf_estimate_cost4Dcont

contains

    !>  Constructors
    subroutine init1D( self, spec_img, lims, reslims, spec,  parms )
        class(ctf_estimate_cost1D), intent(inout) :: self
        class(image),             intent(in)    :: spec_img
        integer,                  intent(in)    :: lims(2), reslims(2)
        real,             target, intent(inout) :: spec(lims(1):lims(2))
        class(ctfparams), target, intent(inout) :: parms
        call self%kill
        self%spec1d => spec
        self%parms  => parms
        self%reslims1d = reslims
        self%ldim      = spec_img%get_ldim()
        self%tfun      = ctf(self%parms%smpd, self%parms%kV, self%parms%Cs, self%parms%fraca)
    end subroutine init1D

    subroutine init2D( self, pspec, parms, inds, ndim, limits, astigtol, tol )
        class(ctf_estimate_cost2D), intent(inout) :: self
        class(image),     target, intent(inout) :: pspec
        class(ctfparams), target, intent(inout) :: parms
        integer,          target, intent(inout) :: inds(:,:)
        integer,                  intent(in)    :: ndim
        real,                     intent(in)    :: limits(ndim,2)
        real,                     intent(in)    :: astigtol, tol
        call self%kill2D
        self%pspec => pspec
        self%parms => parms
        self%inds  => inds
        self%ldim  = self%pspec%get_ldim()
        self%flims = self%pspec%loop_lims(3)
        self%ninds = size(self%inds,dim=2)
        self%ndim  = ndim
        self%astigtol = astigtol
        self%tfun = ctf(self%parms%smpd, self%parms%kV, self%parms%Cs, self%parms%fraca)
        call self%ospec%specify('de',self%ndim,limits=limits, maxits=400, ftol=tol)
        call self%ospec%set_costfun(cost2D_wrapper)
    end subroutine init2D

    subroutine init4Dcont( self, pspec, parms, inds, ndim, limits, astigtol)
        use simple_opt_factory, only: opt_factory
        class(ctf_estimate_cost4Dcont), intent(inout) :: self
        class(image),     target, intent(inout) :: pspec
        class(ctfparams), target, intent(inout) :: parms
        integer,          target, intent(inout) :: inds(:,:)
        integer,                  intent(in)    :: ndim
        real,                     intent(in)    :: limits(ndim,2)
        real,                     intent(in)    :: astigtol
        type(opt_factory) :: opt_fact
        real(dp)          :: amp_contr, phaseq
        real              :: optlims(ndim,2)
        call self%kill4Dcont
        self%pspec    => pspec
        self%parms    => parms
        self%inds     => inds
        self%ldim     = self%pspec%get_ldim()
        self%flims    = self%pspec%loop_lims(3)
        self%ninds    = size(self%inds,dim=2)
        self%ndim     = ndim
        if( self%ndim < 2 ) THROW_HARD('Multidimensional optimization only!')
        self%astigtol = real(astigtol,dp)
        ! set constants, cf. ctf type
        self%smpd             = real(self%parms%smpd,dp)
        self%conversionfactor = 1.d4/self%smpd ! microns to pixels
        self%kV               = real(self%parms%kV,dp)
        self%wl               = 12.26d0 / dsqrt(1.d3*self%kV + 0.9784d0*self%kV*self%kv) / self%smpd
        self%wlsq             = self%wl*self%wl
        self%cs               = real(self%parms%Cs,dp)*1.d7 / self%smpd
        amp_contr             = real(self%parms%fraca,dp)
        phaseq                = dsqrt(1.d0-amp_contr*amp_contr)
        self%amp_contr_const  = datan(amp_contr / phaseq)
        self%angast           = real(self%parms%angast,dp)*DPI/180.d0 ! astigmatism in radians
        self%phshift          = 0.d0
        if( self%parms%l_phaseplate ) self%phshift = real(self%parms%phshift,dp)
        ! optimizer
        optlims = limits
        if( self%ndim >= 3 ) optlims(3,:) = deg2rad(optlims(3,:)) ! astigmatism in radians
        call self%ospec%specify('lbfgsb', self%ndim, ftol=1e-1, gtol=1e-3,&
            &limits=optlims, max_step=0.01, maxits=100)
        call opt_fact%new(self%ospec, self%nlopt)
        self%ospec%costfun_8    => f_costfun
        self%ospec%gcostfun_8   => df_costfun
        self%ospec%fdfcostfun_8 => fdf_costfun
    end subroutine init4Dcont

    ! 1D ROUTINE

    real function cost1D( self, df )
        class(ctf_estimate_cost1D), intent(inout) :: self
        real,                       intent(in)    :: df         !< average defocus
        real          :: ang, spaFreqSq, hinv, inv_ldim, phshift
        real(dp)      :: ctf_sqsum,dotproduct,tvalsq,tval,corr,ctf_sum
        integer       :: h, n
        ! assumes that the 1d spectrum has zero mean and
        ! unit variance over the resolution range
        phshift = 0
        if( self%parms%l_phaseplate )phshift = self%parms%phshift
        call self%tfun%init(df, df, 0.)
        n          = self%reslims1d(2)-self%reslims1d(1)+1
        inv_ldim   = 1./real(self%ldim(1))
        ctf_sqsum  = 0.d0
        dotproduct = 0.d0
        ctf_sum    = 0.d0
        do h = self%reslims1d(1),self%reslims1d(2)
            ! |CTF|
            hinv      = real(h) * inv_ldim
            spaFreqSq = hinv * hinv
            ang       = atan2(0.,real(h))
            tval      = real(self%tfun%eval(spaFreqSq, 0., phshift),dp)
            tvalsq    = min(1.d0,max(tval*tval,0.000001d0))
            tval      = dsqrt(tvalsq)
            ! correlation sums
            ctf_sum    = ctf_sum    + tval
            ctf_sqsum  = ctf_sqsum  + tvalsq
            dotproduct = dotproduct + tval*real(self%spec1d(h),dp)
        end do
        ! one pass correlation
        corr = dotproduct / dsqrt(real(n,dp)*ctf_sqsum-ctf_sum*ctf_sum)
        ! cost
        cost1D = -real(corr,sp)
    end function cost1D

    ! 2D DISCRETE ROUTINES

    !>  Optimization routine
    subroutine minimize( self, parms, cc )
        class(ctf_estimate_cost2D), intent(inout) :: self
        class(ctfparams),           intent(inout) :: parms
        real,                       intent(out)   :: cc
        real :: cost_here
        self%parms%dfx     = parms%dfx
        self%parms%dfy     = parms%dfy
        self%parms%angast  = parms%angast
        self%parms%phshift = parms%phshift
        self%parms%l_phaseplate = parms%l_phaseplate
        select case(self%ndim)
            case(2)
                ! defocus
                self%ospec%x = [parms%dfx,parms%dfy]
            case(3)
                ! defocus & astigmatism
                self%ospec%x = [parms%dfx,parms%dfy,parms%angast]
            case(4)
                ! defocus, astigmatism & phase shift
                if( .not.parms%l_phaseplate ) THROW_HARD('Inconsistency; minimize')
                self%ospec%x = [parms%dfx,parms%dfy,parms%angast,parms%phshift]
        end select
        call self%diffevol%new(self%ospec)
        call self%diffevol%minimize(self%ospec, self, cost_here)
        ! report
        cc = -cost_here
        parms%dfx = self%ospec%x(1)
        parms%dfy = self%ospec%x(2)
        select case(self%ndim)
            case(3)
                parms%angast = self%ospec%x(3)
            case(4)
                parms%phshift = self%ospec%x(4)
        end select
    end subroutine minimize

    !>  Wrapper for optimizer
    real function cost2D_wrapper( self, vec, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        cost2D_wrapper = 1.
        select type(self)
        class is(ctf_estimate_cost2D)
                cost2D_wrapper = self%cost2D(D, vec)
            class DEFAULT
                THROW_HARD('Forbidden type; cost_wrapper')
        end select
    end function cost2D_wrapper

    !>  Cost evaluation fork
    real function cost2D( self, D, vec )
        class(ctf_estimate_cost2D), intent(inout) :: self
        integer,                    intent(in)    :: D
        real,                       intent(in)    :: vec(D)
        cost2D = 1.
        if( abs(vec(1)-vec(2)) > 2.*self%astigtol )then
            return
        else
            select case(D)
                case(1)
                    ! 1D defocus
                    cost2D = self%calc_cost2D(vec(1),vec(1),0.,0.)
                case(2)
                    ! 2D defocus
                    if( self%parms%l_phaseplate )then
                        cost2D = self%calc_cost2D(vec(1),vec(2),self%parms%angast,add_phshift=self%parms%phshift)
                    else
                        cost2D = self%calc_cost2D(vec(1),vec(2),self%parms%angast)
                    endif
                case(3)
                    ! defocus & astigmatism
                    if( self%parms%l_phaseplate )then
                        cost2D = self%calc_cost2D(vec(1),vec(2),vec(3),add_phshift=self%parms%phshift)
                    else
                        cost2D = self%calc_cost2D(vec(1),vec(2),vec(3))
                    endif
                case(4)
                    ! defocus,astigmatism & phase-shift
                    cost2D = self%calc_cost2D(vec(1),vec(2),vec(3),vec(4))
            end select
        endif
    end function cost2D

    ! cost function is correlation within resolution mask between the CTF
    ! spectrum (the model) and the pre-processed micrograph spectrum (the data)
    real function calc_cost2D( self, dfx, dfy, angast, add_phshift )
        class(ctf_estimate_cost2D), intent(inout) :: self
        real,                       intent(in)    :: dfx         !< defocus x-axis
        real,                       intent(in)    :: dfy         !< defocus y-axis
        real,                       intent(in)    :: angast      !< angle of astigmatism
        real,             optional, intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
        real, pointer :: prmat(:,:,:)
        real          :: ang, spaFreqSq, hinv, aadd_phshift, kinv, inv_ldim(3)
        real(dp)      :: ctf_sqsum,dotproduct,tvalsq,tval,corr
        integer       :: h,mh,k,mk, i,j, l
        ! initialize, it assumes that the reference(img) ins zero mean and
        ! unit variance over the resolution range
        aadd_phshift = 0.
        if( present(add_phshift) ) aadd_phshift = add_phshift
        call self%tfun%init(dfx, dfy, angast)
        call self%pspec%get_rmat_ptr(prmat)
        mh = abs(self%flims(1,1))
        mk = abs(self%flims(2,1))
        inv_ldim   = 1./real(self%ldim)
        ctf_sqsum  = 0.d0
        dotproduct = 0.d0
        !$omp parallel do default(shared) private(h,hinv,k,kinv,i,j,spaFreqSq,ang,tval,tvalsq,l) &
        !$omp schedule(static) proc_bind(close) reduction(+:ctf_sqsum,dotproduct)
        do l = 1,self%ninds
            h = self%inds(1,l)
            k = self%inds(2,l)
            i = min(max(1,h+mh+1),self%ldim(1))
            j = min(max(1,k+mk+1),self%ldim(2))
            ! |CTF|
            hinv      = real(h) * inv_ldim(1)
            kinv      = real(k) * inv_ldim(2)
            spaFreqSq = hinv * hinv + kinv * kinv
            ang       = atan2(real(k),real(h))
            tval      = real(self%tfun%eval(spaFreqSq, ang, aadd_phshift),dp)
            tvalsq    = min(1.d0,max(tval*tval,0.000001d0))
            tval      = dsqrt(tvalsq)
            ! correlation sums
            ctf_sqsum  = ctf_sqsum  + tvalsq
            dotproduct = dotproduct + tval*real(prmat(i,j,1),dp)
        end do
        !$omp end parallel do
        ! because reference is normalized (zero-mean, unit variance)
        corr = dotproduct / dsqrt(ctf_sqsum*real(self%ninds,dp))
        ! cost
        calc_cost2D = -real(corr,sp)
    end function calc_cost2D

    ! 4D CONTINUOUS ROUTINES

    !>  Wrappers for optimizer
    subroutine fdf_costfun( self, vec, f, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: f, grad(D)
        select type(self)
        class is(ctf_estimate_cost4Dcont)
                call self%fdf(D, vec, f, grad)
            class DEFAULT
                THROW_HARD('Forbidden type; fdf_costfun')
        end select
    end subroutine fdf_costfun

    subroutine df_costfun( self, vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        real(dp) :: f
        select type(self)
        class is(ctf_estimate_cost4Dcont)
                call self%fdf(D, vec, f, grad )
            class DEFAULT
                THROW_HARD('Forbidden type; fdf_costfun')
        end select
    end subroutine df_costfun

    real(dp) function f_costfun( self, vec, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(in)    :: vec(D)
        real(dp) :: grad(D)
        f_costfun = 1.d0
        select type(self)
        class is(ctf_estimate_cost4Dcont)
                call self%fdf(D, vec, f_costfun, grad )
            class DEFAULT
                THROW_HARD('Forbidden type; f_costfun')
        end select
    end function f_costfun

    !>  Cost & gradients evaluation
    subroutine fdf( self, D, vec, f, grad )
        class(ctf_estimate_cost4Dcont), intent(inout) :: self
        integer,                        intent(in)    :: D
        real(dp),                       intent(in)    :: vec(D)
        real(dp),                       intent(inout) :: f, grad(D)
        real, pointer :: prmat(:,:,:)
        real(dp)      :: dcc_dangast(2),dcc_ddfx(2),dcc_ddfy(2),dcc_dphshift(2),gpen(2), inv_ldim(3)
        real(dp)      :: dfx,dfy,angast,phshift, df_eff,ang,spaFreqSq,costerm,signctf,  denom,denomsq
        real(dp)      :: ctfsq_sum,dotproduct,rn,abstval,tval,tvalsh,hinv,kinv,ddf_dangast,corr,fpen
        real(dp)      :: dabsctf_ddf,ddf_ddfx,ddf_ddfy,phshift_eff,absctf_obs,A,dabsctf_dphshift
        integer       :: h,mh,k,mk, i,j, l
        ! variables
        dfx = vec(1) * self%conversionfactor    ! pixels
        dfy = vec(2) * self%conversionfactor    ! pixels
        angast = self%angast
        if( D >= 3 ) angast = vec(3)            ! radians
        phshift = self%phshift
        if( D == 4 ) phshift = vec(4)           ! radians
        ! init
        call self%pspec%get_rmat_ptr(prmat)
        mh           = abs(self%flims(1,1))
        mk           = abs(self%flims(2,1))
        inv_ldim     = 1.d0/real(self%ldim,dp)
        rn           = real(self%ninds,dp)
        ctfsq_sum    = 0.d0
        dotproduct   = 0.d0
        dcc_ddfx     = 0.d0
        dcc_ddfy     = 0.d0
        dcc_dangast  = 0.d0
        dcc_dphshift = 0.d0
        !$omp parallel do default(shared) private(h,hinv,k,kinv,i,j,l,spaFreqSq,ang,costerm,df_eff,phshift_eff) &
        !$omp private(tval,abstval,absctf_obs,signctf,tvalsh,dabsctf_ddf,ddf_ddfx,ddf_ddfy,ddf_dangast,dabsctf_dphshift)&
        !$omp schedule(static) proc_bind(close) reduction(+:ctfsq_sum,dotproduct,dcc_ddfx,dcc_ddfy,dcc_dangast,dcc_dphshift)
        do l = 1,self%ninds
            h = self%inds(1,l)
            k = self%inds(2,l)
            i = min(max(1,h+mh+1),self%ldim(1))
            j = min(max(1,k+mk+1),self%ldim(2))
            ! |CTF|
            hinv        = real(h,dp) * inv_ldim(1)
            kinv        = real(k,dp) * inv_ldim(2)
            spaFreqSq   = hinv * hinv + kinv * kinv
            ang         = datan2(real(k,dp),real(h,dp))
            costerm     = dcos(2.d0*(ang-angast))
            df_eff      = 0.5d0*(dfx+dfy + costerm*(dfx-dfy))
            phshift_eff = DPI*self%wl*spaFreqSq*(df_eff-0.5d0*self%wlsq*spaFreqSq*self%Cs) + phshift + self%amp_contr_const
            tval        = dsin(phshift_eff)
            abstval     = dabs(tval)
            absctf_obs  = real(prmat(i,j,1),dp)
            ! sums
            dotproduct = dotproduct + abstval*absctf_obs
            ctfsq_sum  = ctfsq_sum  + abstval*abstval
            ! gradients
            signctf     = dsign(1.d0, tval)
            tvalsh      = dcos(phshift_eff)
            dabsctf_ddf = signctf * DPI*self%wl*spaFreqSq * tvalsh
            ! gradients defocus
            ddf_ddfx    = 0.5d0*(1.d0+costerm)
            ddf_ddfy    = 0.5d0*(1.d0-costerm)
            dcc_ddfx    = dcc_ddfx + dabsctf_ddf*ddf_ddfx*[absctf_obs, abstval]
            dcc_ddfy    = dcc_ddfy + dabsctf_ddf*ddf_ddfy*[absctf_obs, abstval]
            if( D >= 3 )then
                ! gradients astigmatism
                ddf_dangast = (dfx-dfy) * dsin(2.d0*(ang-angast))
                dcc_dangast = dcc_dangast + dabsctf_ddf*ddf_dangast*[absctf_obs, abstval]
                if( D == 4 )then
                    ! gradients phase shift
                    dabsctf_dphshift = signctf * tvalsh
                    dcc_dphshift     = dcc_dphshift + dabsctf_dphshift*[absctf_obs, abstval]
                endif
            endif
        end do
        !$omp end parallel do
        denomsq = rn * ctfsq_sum
        denom   = dsqrt(denomsq)
        A       = dotproduct*rn/denomsq**1.5d0
        ! correlation & gradients
        corr    = dotproduct / denom
        grad(1) = self%conversionfactor * (dcc_ddfx(1)/denom - dcc_ddfx(2)*A)
        grad(2) = self%conversionfactor * (dcc_ddfy(1)/denom - dcc_ddfy(2)*A)
        if( D >= 3 )then
            grad(3) = dcc_dangast(1)/denom - dcc_dangast(2)*A
            if( D == 4 ) grad(4) = dcc_dphshift(1)/denom - dcc_dphshift(2)*A
        endif
        ! atigmatism penalty term
        fpen    = ((vec(1)-vec(2))/self%astigtol)**2.d0 / (2.d0*rn)
        gpen(1) = (vec(1)-vec(2)) / (rn*self%astigtol**2.d0)
        gpen(2) = -gpen(1)
        ! the end
        f         = fpen-corr
        grad(1:2) = grad(1:2)-gpen
        grad      = -grad
    end subroutine fdf

    !>  Optimization routine
    subroutine minimize4Dcont( self, parms, cc )
        class(ctf_estimate_cost4Dcont), intent(inout) :: self
        class(ctfparams),               intent(inout) :: parms
        real,                           intent(out)   :: cc
        real :: cost_here
        self%parms%dfx     = parms%dfx
        self%parms%dfy     = parms%dfy
        self%parms%angast  = parms%angast
        self%parms%phshift = parms%phshift
        self%parms%l_phaseplate = parms%l_phaseplate
        select case(self%ndim)
            case(2)
                ! defocus
                self%ospec%x = [parms%dfx,parms%dfy]
            case(3)
                ! defocus + astimgatism
                self%ospec%x = [parms%dfx,parms%dfy,deg2rad(parms%angast)]
            case(4)
                ! defocus + astimgatism + phase shift
                self%ospec%x = [parms%dfx,parms%dfy,deg2rad(parms%angast),parms%phshift]
        end select
        self%ospec%x_8 = real(self%ospec%x,dp)
        call self%nlopt%minimize(self%ospec, self, cost_here)
        ! report
        cc = -real(cost_here)
        parms%dfx = self%ospec%x(1)
        parms%dfy = self%ospec%x(2)
        select case(self%ndim)
            case(3)
                parms%angast  = rad2deg(self%ospec%x(3)) ! radians to degrees
            case(4)
                parms%angast  = rad2deg(self%ospec%x(3)) ! radians to degrees
                parms%phshift = self%ospec%x(4)          ! radians
        end select
    end subroutine minimize4Dcont

    !>  DESTRUCTORS

    subroutine kill1D( self )
        class(ctf_estimate_cost1D), intent(inout) :: self
        self%pspec => null()
        self%parms => null()
    end subroutine kill1D

    subroutine kill2D( self )
        class(ctf_estimate_cost2D), intent(inout) :: self
        self%pspec => null()
        self%parms => null()
        self%inds  => null()
        call self%ospec%kill
        call self%diffevol%kill
    end subroutine kill2D

    subroutine kill4Dcont( self )
        class(ctf_estimate_cost4Dcont), intent(inout) :: self
        self%pspec => null()
        self%parms => null()
        self%inds  => null()
        if(associated(self%nlopt)) call self%nlopt%kill
        self%nlopt => null()
    end subroutine kill4Dcont

end module simple_ctf_estimate_cost
