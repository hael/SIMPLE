module simple_ctf_estimate_patch
include 'simple_lib.f08'
use simple_image,    only: image
use simple_ctf,      only: ctf
use simple_opt_spec, only: opt_spec
use simple_opt_de,   only: opt_de
implicit none

public :: ctf_estimate_patch
private
#include "simple_local_flags.inc"

type ctf_estimate_patch
    private
    class(image), pointer :: pspec                   ! all micrograph powerspec
    type(image)           :: pspec_ctf               ! CTF powerspec
    type(image)           :: pspec_ctf_roavg         ! rotationally averaged CTF powerspec
    type(image)           :: pspec_roavg             ! rotationally averaged all micrograph powerspec
    type(image)           :: imgmsk                  ! mask image
    type(ctf)             :: tfun                    ! transfer function object
    type(opt_spec)        :: ospec_de                ! optimiser specification differential evolution (DE)
    type(opt_de)          :: diffevol                ! DE search object
    type(ctfparams)       :: parms                   ! for storing ctf parameters
    integer, allocatable  :: inds_msk(:,:)
    logical, allocatable  :: cc_msk(:,:,:)           ! corr mask
    real                  :: df_lims(2)   = [0.3,5.0]
    real                  :: df_step      = 0.05     ! defocus step for grid search
    real                  :: astigtol     = 0.05     ! tolerated astigmatism
    real                  :: hp           = 0.       ! high-pass limit
    real                  :: lp           = 0.       ! low-pass limit
    real                  :: limits(4,2)  = 0.       ! barrier limits
    real                  :: cc_fit       = -1.
    real                  :: cc90         = -1.
    real                  :: ctfscore     = -1.
    integer               :: flims(3,2)   = 0
    integer               :: ldim(3)      = 0        ! logical dimensions
    integer               :: ndim         = 0        ! # optimisation dims
    integer               :: npix_msk     = 0
    logical               :: l_patch      = .false.
    logical               :: exists       = .false.
contains
    ! constructor
    procedure          :: new
    ! getters
    procedure          :: get_cc90
    procedure          :: get_ccfit
    procedure          :: get_ctfscore
    ! doers
    procedure          :: fit
    procedure, private :: norm_pspec
    procedure, private :: gen_resmsk
    procedure, private :: calc_ctfscore
    procedure, private :: write_diagnostic
    procedure, private :: grid_srch
    procedure, private :: refine
    procedure, private :: refine_patch
    procedure, private :: ctf2pspecimg
    ! cost function
    procedure, private :: calc_cost
    procedure, private :: cost
    procedure, private :: cost_patch
    ! destructor
    procedure          :: kill
end type ctf_estimate_patch

real,    parameter :: FTOL_REFINE = 1.e-4
integer, parameter :: IARES = 5, NSTEPS = 100

contains

    subroutine new(self, pspec, parms, dfrange, resrange, astigtol_in, l_patch)
        class(ctf_estimate_patch), intent(inout) :: self
        class(image), target, intent(inout) :: pspec       !< all micrograph powerspec
        class(ctfparams),     intent(in)    :: parms
        real,                 intent(in)    :: dfrange(2)  !< defocus range, [30.0,5.0] default
        real,                 intent(in)    :: resrange(2) !< resolution range, [30.0,5.0] default
        real,                 intent(in)    :: astigtol_in !< tolerated astigmatism, 0.05 microns default
        logical,              intent(in)    :: l_patch
        type(image) :: pspec_rot90 ! 90 deg rotated pspec
        call self%kill
        ! set constants
        self%parms%smpd         = parms%smpd
        self%parms%cs           = parms%Cs
        self%parms%kv           = parms%kV
        self%parms%fraca        = parms%fraca
        self%parms%l_phaseplate = parms%l_phaseplate
        self%pspec => pspec
        self%flims = self%pspec%loop_lims(3)
        self%ldim  = self%pspec%get_ldim()
        if( dfrange(1) < dfrange(2) )then
            self%df_lims = dfrange
            self%df_step = (self%df_lims(2) - self%df_lims(1)) / real(NSTEPS)
        else
            THROW_HARD('invalid defocus range; ctf_estimate_init')
        endif
        if( resrange(1) > resrange(2) )then
            self%hp = resrange(1)
            self%lp = resrange(2)
        else
            THROW_HARD('invalid resolution range; new')
        endif
        self%l_patch  = l_patch
        self%astigtol = astigtol_in
        if( self%l_patch )then
            self%ndim = 2
        else
            if(self%parms%l_phaseplate)then
                self%ndim = 4
            else
                self%ndim = 3
            endif
        endif
        ! construct CTF objects
        self%tfun = ctf(self%parms%smpd, self%parms%kV, self%parms%Cs, self%parms%fraca)
        ! generate correlation mask
        call self%gen_resmsk
        ! normalize power spectrum with respect to resolution range
        call self%norm_pspec(self%pspec)
        ! prepare CTF power spectra
        call self%pspec_ctf%new(self%ldim, self%parms%smpd)
        ! prepare rotationally averaged power spectra & CTF power spectrum
        if( self%l_patch )then
            ! done
        else
            call self%pspec%roavg(IARES, self%pspec_roavg, 180)
            call self%norm_pspec(self%pspec_roavg)
        endif
        call self%pspec_ctf_roavg%new(self%ldim, self%parms%smpd)
        ! calculate CTF quality score based on corr with 90 deg rotated
        call self%pspec%rtsq(90., 0., 0., pspec_rot90)
        self%cc90 = pspec%real_corr(pspec_rot90, self%cc_msk)
        call pspec_rot90%kill
        ! limits
        self%limits = 0.
        ! random seed
        call seed_rnd
        self%exists = .true.
    end subroutine new

    ! GETTERS

    real function get_cc90(self)
        class(ctf_estimate_patch), intent(inout) :: self
        get_cc90 = self%cc90
    end function get_cc90

    real function get_ccfit(self)
        class(ctf_estimate_patch), intent(inout) :: self
        get_ccfit = self%cc_fit
    end function get_ccfit

    real function get_ctfscore(self)
        class(ctf_estimate_patch), intent(inout) :: self
        get_ctfscore = self%ctfscore
    end function get_ctfscore

    ! DOERS

    subroutine fit( self, parms, diagfname )
        class(ctf_estimate_patch), intent(inout) :: self
        type(ctfparams),           intent(inout) :: parms
        character(len=*),          intent(in)    :: diagfname
        if( self%l_patch )then
            ! init
            self%parms%dfx     = parms%dfx
            self%parms%dfy     = parms%dfy
            self%parms%angast  = parms%angast
            self%parms%phshift = parms%phshift
            self%parms%l_phaseplate = parms%l_phaseplate
            ! refine
            call self%refine_patch
            ! make a half-n-half diagnostic
            call self%write_diagnostic(diagfname)
        else
            ! 1D grid search with rotational average
            call self%grid_srch
            ! 3/4D refinement of grid solution
            call self%refine
            ! make a half-n-half diagnostic
            call self%write_diagnostic(diagfname)
            ! calculate CTF score diagnostic
            call self%calc_ctfscore
        endif
        parms%dfx          = self%parms%dfx
        parms%dfy          = self%parms%dfy
        parms%angast       = self%parms%angast
        parms%phshift      = self%parms%phshift
        parms%l_phaseplate = self%parms%l_phaseplate
    end subroutine fit

    !>  \brief  Normalize to zero mean and unit variance the reference power spectrum
    !>  within the relevent resolution range
    subroutine norm_pspec( self, img )
        class(ctf_estimate_patch), intent(inout) :: self
        class(image),              intent(inout) :: img
        real, pointer :: prmat(:,:,:)
        real(dp)      :: avg, sdev, val
        integer       :: i,j,k, mh,mk
        call img%get_rmat_ptr(prmat)
        mh   = maxval(self%flims(1,:))
        mk   = maxval(self%flims(2,:))
        sdev = 0.d0
        avg  = 0.d0
        do i=1,self%npix_msk
            j = min(max(1,self%inds_msk(1,i)+mh+1),self%ldim(1))
            k = min(max(1,self%inds_msk(2,i)+mk+1),self%ldim(2))
            avg = avg + real(prmat(j,k,1),dp)
        enddo
        avg = avg / real(self%npix_msk,dp)
        do i=1,self%npix_msk
            j = min(max(1,self%inds_msk(1,i)+mh+1),self%ldim(1))
            k = min(max(1,self%inds_msk(2,i)+mk+1),self%ldim(2))
            val  = real(prmat(j,k,1),dp) - avg
            sdev = sdev + val*val
        enddo
        sdev = dsqrt(sdev/real(self%npix_msk,dp))
        if(sdev <= TINY) sdev = 1.d0
        do i=1,self%npix_msk
            j = min(max(1,self%inds_msk(1,i)+mh+1),self%ldim(1))
            k = min(max(1,self%inds_msk(2,i)+mk+1),self%ldim(2))
            prmat(j,k,1) = real((real(prmat(j,k,1),dp) - avg)/sdev)
        enddo
        nullify(prmat)
    end subroutine norm_pspec

    ! builds resolution dependent mask and indices for correlation calculation
    subroutine gen_resmsk( self )
        class(ctf_estimate_patch), intent(inout) :: self
        logical, allocatable :: nr_msk(:,:,:)
        integer :: h,k, i,j, cnt, mh,mk
        call self%imgmsk%new(self%ldim, self%parms%smpd)
        call self%imgmsk%resmsk(self%hp, self%lp)
        self%cc_msk = self%imgmsk%bin2logical()
        mh   = maxval(self%flims(1,:))
        mk   = maxval(self%flims(2,:))
        nr_msk = self%cc_msk
        nr_msk(1:self%ldim(1)/2-1,:,1)           = .false.
        nr_msk(self%ldim(1)/2,self%ldim(2)/2:,1) = .false.
        self%npix_msk = count(nr_msk)
        allocate(self%inds_msk(2,self%npix_msk))
        cnt = 0
        do h=self%flims(1,1),self%flims(1,2)
            do k=self%flims(2,1),self%flims(2,2)
                i = min(max(1,h+mh+1),self%ldim(1))
                j = min(max(1,k+mk+1),self%ldim(2))
                if( nr_msk(i,j,1) )then
                    cnt = cnt + 1
                    self%inds_msk(:,cnt) = [h,k]
                endif
            enddo
        enddo
    end subroutine gen_resmsk

    ! calculate CTF score diagnostic
    subroutine calc_ctfscore( self )
        class(ctf_estimate_patch), intent(inout) :: self
        real, allocatable :: corrs(:)
        real              :: df_avg
        integer           :: filtsz, hpfind, lpfind
        df_avg = (self%parms%dfx + self%parms%dfy) / 2.0
        call self%ctf2pspecimg(self%pspec_ctf_roavg, df_avg, df_avg, 0.)
        hpfind = self%pspec_roavg%get_find(self%hp)
        lpfind = self%pspec_roavg%get_find(self%lp)
        filtsz = self%pspec_roavg%get_filtsz()
        call self%pspec_roavg%mask(real(lpfind), 'soft', inner=real(hpfind))
        call self%pspec_ctf_roavg%mask(real(lpfind), 'soft', inner=real(hpfind))
        call self%pspec_roavg%norm_bin
        call self%pspec_ctf_roavg%norm_bin
        allocate(corrs(filtsz))
        call self%pspec_roavg%frc_pspec(self%pspec_ctf_roavg, corrs)
        self%ctfscore = real(count(corrs(hpfind:lpfind) > 0.)) / real(lpfind - hpfind + 1)
    end subroutine calc_ctfscore

    ! make & write half-n-half diagnostic
    subroutine write_diagnostic( self, diagfname )
        class(ctf_estimate_patch), intent(inout) :: self
        character(len=*),          intent(in)    :: diagfname
        type(image) :: pspec_half_n_half
        if( self%parms%l_phaseplate )then
            call self%ctf2pspecimg(self%pspec_ctf, self%parms%dfx, self%parms%dfy, self%parms%angast, add_phshift=self%parms%phshift)
        else
            call self%ctf2pspecimg(self%pspec_ctf, self%parms%dfx, self%parms%dfy, self%parms%angast)
        endif
        call self%pspec_ctf%norm()
        call self%pspec%norm()
        call self%pspec%before_after(self%pspec_ctf, pspec_half_n_half, self%cc_msk)
        call pspec_half_n_half%scale_pspec4viz
        call pspec_half_n_half%write_jpg(trim(diagfname), norm=.true.)
        call pspec_half_n_half%kill
    end subroutine write_diagnostic

    subroutine grid_srch( self )
        class(ctf_estimate_patch), intent(inout) :: self
        real :: cost, df, cost_lowest, df_best
        if( self%l_patch ) THROW_HARD('wrong routine; grid_srch')
        df          = self%df_lims(1)
        df_best     = df
        cost_lowest = huge(cost_lowest)
        ! do a first grid search assuming no astigmatism
        do while( df <= self%df_lims(2) )
            if( self%parms%l_phaseplate )then
                cost = self%calc_cost(self%pspec_roavg, df, df, 0., add_phshift=PIO2)
            else
                cost = self%calc_cost(self%pspec_roavg, df, df, 0.)
            endif
            if( cost < cost_lowest )then
                cost_lowest = cost
                df_best     = df
            endif
            df = df + self%df_step
        end do
        self%parms%dfx     = df_best
        self%parms%dfy     = df_best
        self%parms%angast  = 0.
        self%parms%phshift = 0.
        if( self%parms%l_phaseplate ) self%parms%phshift = PIO2
        self%cc_fit = -cost_lowest
    end subroutine grid_srch

    subroutine refine_patch( self )
        class(ctf_estimate_patch), intent(inout) :: self
        real :: cost
        self%limits(1,1) = max(self%df_lims(1),self%parms%dfx-0.1)
        self%limits(1,2) = min(self%df_lims(2),self%parms%dfx+0.1)
        self%limits(2,1) = max(self%df_lims(1),self%parms%dfy-0.1)
        self%limits(2,2) = min(self%df_lims(2),self%parms%dfy+0.1)
        ! init minimizer
        call self%ospec_de%specify('de', self%ndim, limits=self%limits(1:self%ndim,:), maxits=400, ftol=FTOL_REFINE)
        call self%ospec_de%set_costfun(cost_patch_wrapper)
        self%ospec_de%x = [self%parms%dfx,self%parms%dfy]
        call self%diffevol%new(self%ospec_de)
        call self%diffevol%minimize(self%ospec_de, self, cost)
        self%parms%dfx = self%ospec_de%x(1)
        self%parms%dfy = self%ospec_de%x(2)
        self%cc_fit    = -cost
    end subroutine refine_patch

    subroutine refine( self )
        class(ctf_estimate_patch), intent(inout) :: self
        real :: cost, half_range
        if( self%l_patch ) THROW_HARD('wrong routine; refine')
        ! re-init limits for local search
        half_range       = max(self%astigtol, self%df_step)
        self%limits      = 0.
        self%limits(1,1) = max(self%df_lims(1),self%parms%dfx - half_range)
        self%limits(2,1) = max(self%df_lims(1),self%parms%dfy - half_range)
        self%limits(1,2) = min(self%df_lims(2),self%parms%dfx + half_range)
        self%limits(2,2) = min(self%df_lims(2),self%parms%dfy + half_range)
        self%limits(3,1) = 0.
        self%limits(3,2) = 180.
        if( self%parms%l_phaseplate )then
            self%limits(4,1) = 0.
            self%limits(4,2) = 3.142
        endif
        ! init minimizer
        call self%ospec_de%specify('de',self%ndim,limits=self%limits(1:self%ndim,:),maxits=400,ftol=FTOL_REFINE)
        call self%ospec_de%set_costfun(cost_wrapper)
        if( self%parms%l_phaseplate )then
            self%ospec_de%x = [self%parms%dfx,self%parms%dfy,self%parms%angast,self%parms%phshift]
        else
            self%ospec_de%x = [self%parms%dfx,self%parms%dfy,self%parms%angast]
        endif
        call self%diffevol%new(self%ospec_de)
        ! search
        call self%diffevol%minimize(self%ospec_de, self, cost)
        self%cc_fit       = -cost
        self%parms%dfx    = self%ospec_de%x(1)
        self%parms%dfy    = self%ospec_de%x(2)
        self%parms%angast = self%ospec_de%x(3)
        self%parms%phshift = 0.
        if( self%parms%l_phaseplate ) self%parms%phshift = self%ospec_de%x(4)
    end subroutine refine

    !>  \brief  is for making a CTF power-spec image
    subroutine ctf2pspecimg( self, img, dfx, dfy, angast, add_phshift )
        class(ctf_estimate_patch), intent(inout) :: self
        class(image),   intent(inout) :: img         !< image (output)
        real,           intent(in)    :: dfx         !< defocus x-axis
        real,           intent(in)    :: dfy         !< defocus y-axis
        real,           intent(in)    :: angast      !< angle of astigmatism
        real, optional, intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
        real, pointer :: prmat(:,:,:)
        real    :: ang, tval, spaFreqSq, hinv, aadd_phshift, kinv, inv_ldim(3)
        integer :: lims(3,2),h,mh,k,mk,ldim(3), i,j
        ! initialize
        aadd_phshift = 0.
        if( present(add_phshift) ) aadd_phshift = add_phshift
        call self%tfun%init(dfx, dfy, angast)
        call img%get_rmat_ptr(prmat)
        prmat    = 0.
        lims     = img%loop_lims(3)
        mh       = maxval(lims(1,:))
        mk       = maxval(lims(2,:))
        ldim     = img%get_ldim()
        inv_ldim = 1./real(ldim)
        !$omp parallel do collapse(2) default(shared) private(h,hinv,k,kinv,i,j,spaFreqSq,ang,tval) &
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                i         = min(max(1,h+mh+1),ldim(1))
                j         = min(max(1,k+mk+1),ldim(2))
                hinv      = real(h) * inv_ldim(1)
                kinv      = real(k) * inv_ldim(2)
                spaFreqSq = hinv * hinv + kinv * kinv
                ang       = atan2(real(k),real(h))
                tval      = self%tfun%eval(spaFreqSq, ang, aadd_phshift)
                tval      = min(1.,max(tval * tval,0.0001))
                prmat(i,j,1) = sqrt(tval)
            end do
        end do
        !$omp end parallel do
    end subroutine ctf2pspecimg

    ! COST FUNCTION RELATED

    real function cost_wrapper( self, vec, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        cost_wrapper = 1.
        select type(self)
            class is (ctf_estimate_patch)
                cost_wrapper = self%cost(D, vec)
            class DEFAULT
                THROW_HARD('Forbidden type; cost_wrapper')
        end select
    end function cost_wrapper

    real function cost( self, D, vec )
        class(ctf_estimate_patch), intent(inout) :: self
        integer,                   intent(in)    :: D
        real,                      intent(in)    :: vec(D)
        cost = 1.0
        if( abs(vec(1) - vec(2)) > self%astigtol )then
            return
        else
            if( self%parms%l_phaseplate )then
                if( D /= 4 ) THROW_HARD('wrong routine; cost 1')
                cost = self%calc_cost(self%pspec,vec(1),vec(2),vec(3),add_phshift=vec(4))
            else
                if( D /= 3 ) THROW_HARD('wrong routine; cost 2')
                cost = self%calc_cost(self%pspec,vec(1),vec(2),vec(3))
            endif
        endif
    end function cost

    real function cost_patch_wrapper( self, vec, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        cost_patch_wrapper = 1.
        select type(self)
            class is (ctf_estimate_patch)
                cost_patch_wrapper = self%cost_patch(D, vec)
            class DEFAULT
                THROW_HARD('Forbidden type; cost_patch_wrapper')
        end select
    end function cost_patch_wrapper

    real function cost_patch( self, D, vec )
        class(ctf_estimate_patch), intent(inout) :: self
        integer,                   intent(in)    :: D
        real,                      intent(in)    :: vec(D)
        cost_patch = 1.0
        if( abs(vec(1) - vec(2)) > self%astigtol )then
            return
        else
            if( D /= 2 ) THROW_HARD('wrong routine; cost_patch')
            cost_patch = self%calc_cost(self%pspec,vec(1),vec(2),self%parms%angast,add_phshift=self%parms%phshift)
        endif
    end function cost_patch

    ! cost function is real-space correlation within resolution mask between the CTF
    ! powerspectrum (the model) and the pre-processed micrograph powerspectrum (the data)
    real function calc_cost( self, img, dfx, dfy, angast, add_phshift )
        class(ctf_estimate_patch), intent(inout) :: self
        class(image),   intent(inout) :: img         !< image (output)
        real,           intent(in)    :: dfx         !< defocus x-axis
        real,           intent(in)    :: dfy         !< defocus y-axis
        real,           intent(in)    :: angast      !< angle of astigmatism
        real, optional, intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
        real, pointer :: prmat(:,:,:)
        real          :: ang, spaFreqSq, hinv, aadd_phshift, kinv, inv_ldim(3)
        real(dp)      :: ctf_sqsum,dotproduct,tvalsq,tval,corr,ctf_sum
        integer       :: h,mh,k,mk,ldim(3), i,j, l
        ! initialize, it assumes that the reference(img) ins zero mean and
        ! unit variance over the resolution range
        aadd_phshift = 0.
        if( present(add_phshift) ) aadd_phshift = add_phshift
        call self%tfun%init(dfx, dfy, angast)
        call img%get_rmat_ptr(prmat)
        ldim       = img%get_ldim()
        mh         = maxval(self%flims(1,:))
        mk         = maxval(self%flims(2,:))
        inv_ldim   = 1./real(ldim)
        ctf_sqsum  = 0.d0
        dotproduct = 0.d0
        ctf_sum    = 0.d0
        !$omp parallel do default(shared) private(h,hinv,k,kinv,i,j,spaFreqSq,ang,tval,tvalsq,l) &
        !$omp schedule(static) proc_bind(close) reduction(+:ctf_sqsum,dotproduct,ctf_sum)
        do l = 1,self%npix_msk
            h = self%inds_msk(1,l)
            k = self%inds_msk(2,l)
            i = min(max(1,h+mh+1),ldim(1))
            j = min(max(1,k+mk+1),ldim(2))
            ! |CTF|
            hinv      = real(h) * inv_ldim(1)
            kinv      = real(k) * inv_ldim(2)
            spaFreqSq = hinv * hinv + kinv * kinv
            ang       = atan2(real(k),real(h))
            tval      = real(self%tfun%eval(spaFreqSq, ang, aadd_phshift),dp)
            tvalsq    = min(1.d0,max(tval*tval,0.0001d0))
            tval      = dsqrt(tvalsq)
            ! correlation sums
            ctf_sum    = ctf_sum    + tval
            ctf_sqsum  = ctf_sqsum  + tvalsq
            dotproduct = dotproduct + tval*real(prmat(i,j,1),dp)
        end do
        !$omp end parallel do
        ! one pass correlation
        ! because reference is normalized (zero-mean, unit variance)
        corr = dotproduct / dsqrt(real(self%npix_msk,dp)*ctf_sqsum-ctf_sum*ctf_sum)
        ! cost
        calc_cost = -real(corr,sp)
    end function calc_cost

    ! DESTRUCTOR

    subroutine kill( self )
        class(ctf_estimate_patch), intent(inout) :: self
        self%cc_fit       = -1.
        self%cc90         = -1.
        self%ctfscore     = -1.
        nullify(self%pspec)
        call self%pspec_ctf%kill
        call self%pspec_ctf_roavg%kill
        call self%pspec_roavg%kill
        call self%imgmsk%kill
        call self%ospec_de%kill
        call self%diffevol%kill
        if( allocated(self%cc_msk) ) deallocate(self%cc_msk)
        if( allocated(self%inds_msk) ) deallocate(self%inds_msk)
        self%exists = .false.
    end subroutine kill

end module simple_ctf_estimate_patch
