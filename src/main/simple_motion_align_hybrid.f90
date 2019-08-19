module simple_motion_align_hybrid
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_error
use simple_image,        only: image
use simple_parameters,   only: params_glob
use simple_ft_expanded,  only: ft_expanded
use simple_ftexp_shsrch, only: ftexp_shsrch
implicit none
public :: motion_align_hybrid
private
#include "simple_local_flags.inc"

real,    parameter :: SMALLSHIFT    = 1.
real,    parameter :: SRCH_TOL      = 1.e-6
real,    parameter :: NIMPROVED_TOL = 1.e-7
integer, parameter :: MAXITS_PHCORR = 30
integer, parameter :: MAXITS_CORR   = 10
integer, parameter :: POLYDIM       = 4

type :: motion_align_hybrid
    private
    type(image),           pointer :: frames(:)                         !< pointer to stack of frames
    type(image),       allocatable :: frames_sh(:)
    type(image)                    :: reference
    type(ft_expanded), allocatable :: frames_ftexp(:)             !< movie frames
    type(ft_expanded), allocatable :: frames_ftexp_sh(:)          !< shifted movie frames
    type(ft_expanded), allocatable :: references_ftexp(:) !< array of global movie sums for parallel refinement
    real,              allocatable :: shifts_toplot(:,:)                !< for plotting
    real,              allocatable :: opt_shifts(:,:)                   !< shifts identified
    real,              allocatable :: frameweights(:)                   !< array of frameweights
    real,              allocatable :: corrs(:)                          !< per-frame correlations
    real(dp)                       :: polyx(POLYDIM), polyy(POLYDIM)    !< polynomial coefficients
    real                           :: ftol           = 1.e-6            !< tolerance parameter for minimizer
    real                           :: gtol           = 1.e-6            !< tolerance parameter for minimizer
    real                           :: hp             = -1.              !< high pass value
    real                           :: lp             = -1.              !< low pass value
    real                           :: bfactor        = -1.              !< b-factor for alignment weights
    real                           :: corr           = -1.              !< correlation
    real                           :: shsrch_tol     = SRCH_TOL         !< tolerance parameter for continuous srch update
    real                           :: smpd           = 0.
    real                           :: resstep        = 0.
    integer                        :: trs            = 30               !< half phase-correlation search bound
    integer                        :: ldim(3)
    integer                        :: maxits_phcorr  = MAXITS_PHCORR    !< maximum number of iterations for phase-correlation search
    integer                        :: maxits_corr    = MAXITS_CORR      !< maximum number of iterations for continuous search
    integer                        :: nframes        = 0                !< number of frames
    integer                        :: iter           = 0                !< iteration number
    integer                        :: fixed_frame    = 1                !< fixed (non-shifted) frame
    integer                        :: px=0 , py=0                    !< patch x/y id
    integer                        :: lp_updates       = 0              !< # of resolution updates performed [0;3]
    logical                        :: l_bfac           = .false.        !< whether to use b-factor weights
    logical                        :: group_frames     = .false.        !< whether to group frames
    logical                        :: fitshifts        = .false.        ! whether to perform iterative incremental shifts fitting
    logical                        :: rand_init_shifts = .false.        !< randomize initial condition?
    logical, public                :: existence        = .false.

contains
    ! Constructor
    procedure          :: new
    ! Frames & memory management
    procedure, private :: init_images
    procedure, private :: dealloc_images
    procedure, private :: init_ftexps
    procedure, private :: dealloc_ftexps
    ! Doers
    procedure, private :: calc_corr2ref
    procedure, private :: calc_shifts
    procedure, private :: shift_frames_gen_ref
    procedure          :: align
    procedure, private :: align_phcorr
    procedure, private :: align_corr
    procedure, private :: shift_wsum_and_calc_corrs
    procedure, private :: recenter_shifts                   !< put shifts rel. to fixed frame
    procedure, private :: calc_rmsd
    ! Trajectory fitting related
    procedure, private :: fit_polynomial
    procedure, private :: polynomial2shift
    ! Getters & setters
    procedure          :: set_shsrch_tol
    procedure          :: set_hp_lp
    procedure          :: set_trs
    procedure          :: get_iter
    procedure          :: set_weights
    procedure          :: get_weights
    procedure          :: set_rand_init_shifts
    procedure          :: get_corr
    procedure          :: get_corrs
    procedure          :: get_opt_shifts
    procedure          :: get_shifts_toplot
    procedure          :: set_fixed_frame
    procedure          :: set_fitshifts
    procedure          :: set_maxits
    procedure          :: set_coords
    procedure          :: get_coords
    procedure          :: is_fitshifts
    procedure          :: set_bfactor
    procedure          :: set_group_frames
    ! Destructor
    procedure          :: kill
end type motion_align_hybrid

contains

    subroutine new( self, frames_ptr )
        class(motion_align_hybrid),       intent(inout) :: self
        type(image), allocatable, target, intent(in)    :: frames_ptr(:)
        call self%kill
        self%trs                   = nint(params_glob%trs)
        self%ftol                  = params_glob%motion_correctftol
        self%gtol                  = params_glob%motion_correctgtol
        self%fixed_frame           = 1
        self%lp_updates            = 0
        self%maxits_phcorr  = MAXITS_PHCORR
        self%maxits_corr    = MAXITS_CORR
        self%hp                    = -1.
        self%lp                    = -1.
        self%shsrch_tol            = SRCH_TOL
        self%bfactor               = -1.
        self%fitshifts             = .false.
        self%group_frames          = .false.
        self%rand_init_shifts      = .false.
        self%frames                => frames_ptr
        self%nframes               =  size(frames_ptr, 1)
        self%resstep               = (params_glob%lpstart - params_glob%lpstop) / 3.
        if ( self%nframes < 2 ) then
            THROW_HARD('nframes < 2; simple_motion_align_hybrid: align')
        end if
        self%smpd = self%frames(1)%get_smpd()
        self%ldim = self%frames(1)%get_ldim()
        allocate(self%frames_sh(self%nframes),&
                &self%shifts_toplot(self%nframes,2), self%opt_shifts(self%nframes,2),&
                &self%corrs(self%nframes), self%frameweights(self%nframes), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('new; simple_motion_align_hybrid')
        self%shifts_toplot = 0.
        self%opt_shifts    = 0.
        self%corrs         = -1.
        self%frameweights  = 1./real(self%nframes)
        self%l_bfac        = .false.
        self%existence     = .true.
    end subroutine new

    subroutine init_images( self )
        class(motion_align_hybrid), intent(inout) :: self
        integer :: iframe
        call self%dealloc_images
        call self%reference%new(self%ldim,self%smpd)
        call self%reference%zero_and_flag_ft
        allocate(self%frames_sh(self%nframes))
        !$omp parallel do default(shared) private(iframe) proc_bind(close)
        do iframe = 1,self%nframes
            call self%frames(iframe)%fft
            call self%frames_sh(iframe)%new(self%ldim,self%smpd,wthreads=.false.)
            call self%frames_sh(iframe)%zero_and_flag_ft
            call self%frames_sh(iframe)%set_cmat(self%frames(iframe))
        enddo
        !$omp end parallel do
    end subroutine init_images

    subroutine dealloc_images( self )
        class(motion_align_hybrid), intent(inout) :: self
        integer :: iframe
        if(allocated(self%frames_sh))then
            do iframe=1,self%nframes
                call self%frames_sh(iframe)%kill
            enddo
            deallocate(self%frames_sh)
        endif
        call self%reference%kill
    end subroutine dealloc_images

    subroutine init_ftexps( self )
        class(motion_align_hybrid), intent(inout) :: self
        integer :: iframe
        call self%dealloc_ftexps
        allocate(self%frames_ftexp(self%nframes),self%frames_ftexp_sh(self%nframes),self%references_ftexp(self%nframes))
        !$omp parallel do default(shared) private(iframe) schedule(static) proc_bind(close)
        do iframe=1,self%nframes
            call self%frames_ftexp(iframe)%new(self%frames(iframe), self%hp, self%lp, .true., bfac=self%bfactor)
            call self%frames_ftexp_sh(iframe)%new(self%frames(iframe), self%hp, self%lp, .false.)
            call self%references_ftexp(iframe)%new(self%frames(iframe), self%hp, self%lp, .false.)
        end do
        !$omp end parallel do
    end subroutine init_ftexps

    subroutine dealloc_ftexps( self )
        class(motion_align_hybrid), intent(inout) :: self
        integer :: iframe
        if( allocated(self%frames_ftexp) )then
            do iframe=1,self%nframes
                call self%frames_ftexp(iframe)%kill
                call self%frames_ftexp_sh(iframe)%kill
                call self%references_ftexp(iframe)%kill
            end do
        endif
    end subroutine dealloc_ftexps


    ! Alignment routines

    subroutine align( self, ini_shifts, frameweights )
        class(motion_align_hybrid), intent(inout) :: self
        real,             optional, intent(in)    :: ini_shifts(self%nframes,2), frameweights(self%nframes)
        if ( .not. self%existence ) then
            THROW_HARD('not instantiated; simple_motion_align_hybrid: align')
        end if
        if (( self%hp < 0. ) .or. ( self%lp < 0.)) then
            THROW_HARD('hp or lp < 0; simple_motion_align_hybrid: align')
        end if
        ! prep for phase-correlation search
        call self%init_images
        ! phase-correlation search
        call self%align_phcorr( ini_shifts, frameweights )
        call self%dealloc_images
        ! prep ftexps for continuous search
        ! call self%init_ftexps
        ! ! continuous search
        ! write(logfhandle,'(a)') '>>> PERFORMING CORRELATION-BASED CONTINUOUS OPTIMIZATION'
        ! call self%align_corr( frameweights )
        ! call self%dealloc_ftexps
    end subroutine align

    ! Phase-correlation based alignment
    subroutine align_phcorr( self, ini_shifts, frameweights )
        class(motion_align_hybrid), intent(inout) :: self
        real,             optional, intent(in)    :: ini_shifts(self%nframes,2), frameweights(self%nframes)
        real    :: opt_shifts_prev(self%nframes, 2), rmsd
        integer :: iter, iframe
        logical :: l_calc_frameweights
        ! frameweights
        l_calc_frameweights = .not.present(frameweights)
        self%frameweights   = 1./real(self%nframes)
        if( .not.l_calc_frameweights ) self%frameweights = frameweights
        self%lp_updates = 0
        self%opt_shifts = 0.
        if( present(ini_shifts) ) self%opt_shifts = ini_shifts
        ! random initialization
        if ( self%rand_init_shifts ) then
            do iframe = 1, self%nframes
                if( iframe == self%fixed_frame ) cycle
                self%opt_shifts(iframe,1) = self%opt_shifts(iframe,1) + (ran3()-.5)*SMALLSHIFT
                self%opt_shifts(iframe,2) = self%opt_shifts(iframe,2) + (ran3()-.5)*SMALLSHIFT
            end do
        end if
        ! shift frames, generate reference & calculates correlation
        call self%shift_frames_gen_ref
        ! main loop
        do iter=1,self%maxits_phcorr
            self%iter       = iter
            opt_shifts_prev = self%opt_shifts
            ! individual optimizations
            !$omp parallel do default(shared) private(iframe) proc_bind(close)
            do iframe = 1,self%nframes
                call self%calc_shifts(self%frames_sh(iframe), self%frameweights(iframe), self%opt_shifts(iframe,:))
            end do
            !$omp end parallel do
            ! recenter shifts
            call self%recenter_shifts(self%opt_shifts)
            ! shift frames, generate reference & calculates correlation
            call self%shift_frames_gen_ref
            ! updates weights
            if( l_calc_frameweights ) self%frameweights = corrs2weights(self%corrs)
            ! convergence
            rmsd = self%calc_rmsd(opt_shifts_prev, self%opt_shifts)
            if( iter > 1 .and. rmsd < 1. )then
                self%lp_updates = self%lp_updates+1
                if( self%lp_updates > 3 )then
                    self%lp_updates = 3
                    exit
                endif
                if( self%fitshifts .and. self%lp_updates <= 2 )then
                    ! optional shifts fitting
                    call self%fit_polynomial(self%opt_shifts)
                    do iframe = 1,self%nframes
                        call self%polynomial2shift(iframe, self%opt_shifts(iframe,:))
                    enddo
                    call self%recenter_shifts(self%opt_shifts)
                endif
                ! resolution update
                self%lp = self%lp - self%resstep
                write(logfhandle,*)'>>> Updating low-pass to', self%lp,' at patch:', self%px, ' - ', self%py
            endif
        enddo
        self%shifts_toplot = self%opt_shifts
        ! cleanup
        do iframe=1,self%nframes
            call self%frames_sh(iframe)%kill
        enddo
        deallocate(self%frames_sh)
        call self%reference%kill
    end subroutine align_phcorr

    subroutine align_corr( self, frameweights )
        class(motion_align_hybrid), intent(inout) :: self
        real, optional,             intent(in)    :: frameweights(self%nframes)
        type(ftexp_shsrch) :: ftexp_srch(self%nframes)
        real    :: opt_shifts_saved(self%nframes,2), opt_shifts_prev(self%nframes, 2), corrfrac
        real    :: frameweights_saved(self%nframes), cxy(3), rmsd, corr_prev, corr_saved, trs, frac_improved
        integer :: iter, iframe, nimproved
        logical :: l_calc_frameweights
        ! frameweights
        l_calc_frameweights = .not.present(frameweights)
        self%frameweights   = 1./real(self%nframes)
        if( .not.l_calc_frameweights ) self%frameweights = frameweights
        ! shift boundaries
        trs = 3.
        ! allocation & extraction
        do iframe=1,self%nframes
            call ftexp_srch(iframe)%new(self%references_ftexp(iframe),&
                self%frames_ftexp_sh(iframe),trs, motion_correct_ftol=self%ftol, motion_correct_gtol=self%gtol)
            call ftexp_srch(iframe)%set_shsrch_tol(self%shsrch_tol)
            ftexp_srch(iframe)%ospec%maxits = 100
        end do
        self%iter = 0
        ! generate movie sum for refinement
        call self%shift_wsum_and_calc_corrs
        ! main loop
        corr_saved = self%corr
        do iter=1,self%maxits_corr
            self%iter       = iter
            nimproved       = 0
            opt_shifts_prev = self%opt_shifts
            corr_prev       = self%corr
            ! individual optimizations
            !$omp parallel do default(shared) private(iframe,cxy) proc_bind(close) reduction(+:nimproved)
            do iframe = 1,self%nframes
                call self%frames_ftexp(iframe)%shift([-self%opt_shifts(iframe,1), -self%opt_shifts(iframe,2), 0.],&
                    &self%frames_ftexp_sh(iframe))
                cxy = ftexp_srch(iframe)%minimize(self%corrs(iframe))
                if( cxy(1) - self%corrs(iframe) > NIMPROVED_TOL ) nimproved = nimproved + 1
                self%opt_shifts(iframe,:) = self%opt_shifts(iframe,:) + cxy(2:3)
                self%corrs(iframe) = cxy(1)
            end do
            !$omp end parallel do
            frac_improved = real(nimproved) / real(self%nframes) * 100.
            ! recenter shifts
            call self%recenter_shifts(self%opt_shifts)
            ! updates weights
            if( l_calc_frameweights ) self%frameweights = corrs2weights(self%corrs)
            ! build new reference
            call self%shift_wsum_and_calc_corrs
            ! convergence
            rmsd = self%calc_rmsd(opt_shifts_prev, self%opt_shifts)
            if( self%corr >= corr_saved ) then
                ! save the local optimum
                corr_saved         = self%corr
                frameweights_saved = self%frameweights
                opt_shifts_saved   = self%opt_shifts
            endif
            corrfrac = corr_prev / self%corr
            print *,iter, 'rmsd',rmsd
            if( corrfrac > 0.999 .and. rmsd < 0.2 ) exit
        end do
        ! best local optimum
        self%corr          = corr_saved
        self%opt_shifts    = opt_shifts_saved
        self%frameweights  = frameweights_saved
        self%shifts_toplot = self%opt_shifts
        ! cleanup
        do iframe = 1, self%nframes
            call ftexp_srch(iframe)%kill
        end do
    end subroutine align_corr

    ! shifts frames, generate reference and calculates correlations
    subroutine shift_frames_gen_ref( self )
        class(motion_align_hybrid), intent(inout) :: self
        complex, allocatable :: cmat_sum(:,:,:)
        complex,     pointer :: pcmat(:,:,:)
        integer :: iframe
        cmat_sum = self%reference%get_cmat()
        cmat_sum = cmplx(0.,0.)
        !$omp parallel default(shared) private(iframe,pcmat) proc_bind(close)
        !$omp do schedule(static) reduction(+:cmat_sum)
        do iframe=1,self%nframes
            call self%frames_sh(iframe)%set_cmat(self%frames(iframe))
            call self%frames_sh(iframe)%shift([-self%opt_shifts(iframe,1),-self%opt_shifts(iframe,2),0.])
            call self%frames_sh(iframe)%get_cmat_ptr(pcmat)
            cmat_sum = cmat_sum + pcmat * self%frameweights(iframe)
        enddo
        !$omp end do
        !$omp single
        call self%reference%set_cmat(cmat_sum)
        !$omp end single
        !$omp do schedule(static)
        do iframe = 1,self%nframes
            self%corrs(iframe) = self%calc_corr2ref(self%frames_sh(iframe), self%frameweights(iframe))
        end do
        !$omp end do
        !$omp end parallel
        self%corr = sum(self%corrs) / real(self%nframes)
    end subroutine shift_frames_gen_ref

    ! band-passed correlation to frame-subtracted reference
    real function calc_corr2ref( self, frame, weight )
        class(motion_align_hybrid), intent(inout) :: self
        class(image),               intent(inout) :: frame
        real,                       intent(in)    :: weight
        complex(dp) :: cref, cframe
        real(dp)    :: corr,sumrefsq,sumframesq,w,bfacw,spafreq_denom_sq,spafreq_sq,bfac_sc
        integer     :: nr_lims(3,2),hh,h,k,shsq,hplimsq,lplimsq
        bfac_sc    = real(self%bfactor,dp) / 4.d0
        spafreq_denom_sq = (real(self%ldim(1),dp)*real(self%smpd,dp))**2.d0
        corr       = 0.d0
        sumrefsq   = 0.d0
        sumframesq = 0.d0
        w = real(weight,dp)
        hplimsq = calc_fourier_index(self%hp,self%ldim(1),self%smpd)**2
        lplimsq = calc_fourier_index(self%lp,self%ldim(1),self%smpd)**2
        nr_lims = self%reference%loop_lims(2)
        do h = nr_lims(1,1),nr_lims(1,2)
            hh = h*h
            if( hh > lplimsq ) cycle
            do k = nr_lims(2,1),nr_lims(2,2)
                shsq = hh+k*k
                if( shsq <= lplimsq .and. shsq >= hplimsq  )then
                    bfacw = 1.d0
                    if( self%l_bfac )then
                        spafreq_sq = real(shsq,dp) / spafreq_denom_sq
                        bfacw      = max(1.d-8,exp(-bfac_sc*spafreq_sq))
                    endif
                    cref       = dcmplx(self%reference%get_fcomp2D(h,k))
                    cframe     = dcmplx(frame%get_fcomp2D(h,k))
                    cref       = cref - w*cframe
                    corr       = corr + bfacw * real(cref*dconjg(cframe),dp)
                    sumrefsq   = sumrefsq   + bfacw * csq(cref)
                    sumframesq = sumframesq + bfacw * csq(cframe)
                endif
            end do
        end do
        if( sumrefsq < DTINY .or. sumframesq < DTINY )then
            corr = 0.d0
        else
            corr = corr / dsqrt(sumrefsq*sumframesq)
        endif
        calc_corr2ref = real(corr)
    end function calc_corr2ref

    ! identifies discretized shifts within search range, frame destroyed on exit
    subroutine calc_shifts( self, frame, weight, shift )
        class(motion_align_hybrid), intent(inout) :: self
        class(image),               intent(inout) :: frame
        real,                       intent(in)    :: weight
        real,                       intent(inout) :: shift(2)
        integer, parameter :: BPWIDTH = 8
        ! complex, pointer   :: pcref(:,:,:), pcframe(:,:,:)
        real,    pointer   :: prframe(:,:,:)
        complex :: cref, cframe, comp
        real    :: prev_shift(2), corr, w, bfacw, spafreq_denom_sq, spafreq_sq
        real    :: rsh, rhplim, rlplim, width, alpha, beta, gamma
        integer :: center(2),pos(2),nr_lims(3,2),phys(3)
        integer :: bphplimsq,hplimsq,bplplimsq,lplimsq, shsq, hh,h,k, i,j, hplim,lplim
        prev_shift = shift
        ! phase correlations
        ! call frame%get_cmat_ptr(pcframe)
        ! call self%reference%get_cmat_ptr(pcref)
        ! pcframe(:,:,1) = (pcref(:,:,1)-weight*pcframe(:,:,1)) * conjg(pcframe(:,:,1))
        ! call frame%ifft
        width     = real(BPWIDTH)
        spafreq_denom_sq = (real(self%ldim(1))*self%smpd)**2.
        nr_lims   = self%reference%loop_lims(2)
        hplim     = calc_fourier_index(self%hp,self%ldim(1),self%smpd)
        lplim     = calc_fourier_index(self%lp,self%ldim(1),self%smpd)
        rhplim    = real(hplim)
        rlplim    = real(lplim)
        hplimsq   = hplim*hplim
        lplimsq   = lplim*lplim
        bphplimsq = max(0,hplim-BPWIDTH)**2
        bplplimsq = min(minval(nr_lims(1:2,2)),lplimsq+BPWIDTH)**2
        do h = nr_lims(1,1),nr_lims(1,2)
            hh = h*h
            do k = nr_lims(2,1),nr_lims(2,2)
                shsq = hh+k*k
                phys = frame%comp_addr_phys([h,k,0])
                if( shsq < bphplimsq .or. shsq > bplplimsq )then
                    call frame%set_cmat_at(phys, cmplx(0.,0.))
                else
                    rsh  = sqrt(real(shsq))
                    ! B-factor weight
                    bfacw = 1.
                    if( self%l_bfac )then
                        spafreq_sq = real(shsq) / spafreq_denom_sq
                        bfacw      = max(0.,min(1.,exp(-0.25*self%bfactor*spafreq_sq)))
                    endif
                    ! filter weight
                    w = 1.
                    if( shsq < hplimsq )then
                        w = (1. - cos(((rsh-rhplim)/width)*PI))/2.          ! high-pass
                    elseif( shsq > lplimsq )then
                        w = (cos(((rsh-(rlplim-width))/width)*pi) + 1.)/2. ! low-pass
                    endif
                    ! product
                    cref   = self%reference%get_cmat_at(phys)
                    cframe = frame%get_cmat_at(phys)
                    comp   = bfacw * w*w * (cref - weight*cframe) * conjg(cframe)
                    call frame%set_cmat_at(phys, comp)
                endif
            end do
        end do
        call frame%ifft
        ! find peak
        call frame%get_rmat_ptr(prframe)
        center = [self%ldim(1)/2+1, self%ldim(2)/2+1]
        corr   = -huge(corr)
        pos    = center
        do i = center(1)-self%trs,center(1)+self%trs
            do j = center(2)-self%trs,center(2)+self%trs
                if( prframe(i,j,1) > corr )then
                    pos  = [i,j]-center
                    corr = prframe(i,j,1)
                endif
            enddo
        enddo
        shift = prev_shift + real(pos)
        ! interpolate
        alpha = prframe(pos(1)+center(1)-1,pos(2)+center(2),1)
        beta  = prframe(pos(1)+center(1)  ,pos(2)+center(2),1)
        gamma = prframe(pos(1)+center(1)+1,pos(2)+center(2),1)
        shift(1) = shift(1) + interp_peak()
        alpha = prframe(pos(1)+center(1),pos(2)+center(2)-1,1)
        beta  = prframe(pos(1)+center(1),pos(2)+center(2)  ,1)
        gamma = prframe(pos(1)+center(1),pos(2)+center(2)+1,1)
        shift(2) = shift(2) + interp_peak()
        ! cleanup
        call frame%zero_and_flag_ft
        contains

            real function interp_peak()
                real :: denom
                interp_peak = 0.
                denom = alpha+gamma-2.*beta
                if( abs(denom) < TINY )return
                interp_peak = 0.5 * (alpha-gamma) / denom
            end function interp_peak

    end subroutine calc_shifts

    ! center shifts with respect to fixed_frame
    subroutine recenter_shifts( self, shifts )
        class(motion_align_hybrid), intent(inout) :: self
        real,     intent(inout) :: shifts(self%nframes,2)
        integer :: iframe
        do iframe=1,self%nframes
            shifts(iframe,:) = shifts(iframe,:) - shifts(self%fixed_frame,:)
            if( abs(shifts(iframe,1)) < 1.e-6 ) shifts(iframe,1) = 0.
            if( abs(shifts(iframe,2)) < 1.e-6 ) shifts(iframe,2) = 0.
        end do
    end subroutine recenter_shifts

    real function calc_rmsd( self, prev_shifts, shifts )
        class(motion_align_hybrid), intent(in) :: self
        real,                       intent(in) :: prev_shifts(self%nframes,2), shifts(self%nframes,2)
        integer :: iframe
        calc_rmsd = 0.
        do iframe = 1,self%nframes
            calc_rmsd = calc_rmsd + sum((shifts(iframe,:)-prev_shifts(iframe,:))**2.)
        enddo
        calc_rmsd = sqrt(calc_rmsd/real(self%nframes))
    end function calc_rmsd

    ! Continuous search routines

    ! shifts frames, generate references, substracts self from references and calculates correlations
    subroutine shift_wsum_and_calc_corrs( self )
        class(motion_align_hybrid), intent(inout) :: self
        complex, allocatable :: cmat_sum(:,:,:)
        complex,     pointer :: pcmat(:,:,:)
        integer :: iframe, flims(3,2)
        flims = self%references_ftexp(1)%get_flims()
        allocate(cmat_sum(flims(1,1):flims(1,2),flims(2,1):flims(2,2),1),source=cmplx(0.,0.))
        !$omp parallel default(shared) private(iframe,pcmat) proc_bind(close)
        ! accumulate shifted sum
        !$omp do schedule(static) reduction(+:cmat_sum)
        do iframe=1,self%nframes
            call self%frames_ftexp(iframe)%shift([-self%opt_shifts(iframe,1),-self%opt_shifts(iframe,2),0.],&
                &self%frames_ftexp_sh(iframe))
            call self%frames_ftexp_sh(iframe)%get_cmat_ptr(pcmat)
            cmat_sum = cmat_sum + pcmat * self%frameweights(iframe)
        end do
        !$omp end do
        !$omp do schedule(static)
        do iframe=1,self%nframes
            ! set references
            call self%references_ftexp(iframe)%set_cmat(cmat_sum)
            ! subtract frame
            call self%references_ftexp(iframe)%subtr(self%frames_ftexp_sh(iframe),w=self%frameweights(iframe))
            ! calc corr
            self%corrs(iframe) = self%references_ftexp(iframe)%corr(self%frames_ftexp_sh(iframe))
        end do
        !$omp end do
        !$omp end parallel
        self%corr = sum(self%corrs)/real(self%nframes)
    end subroutine shift_wsum_and_calc_corrs

    ! Getters/setters

    subroutine set_hp_lp( self, hp, lp )
        class(motion_align_hybrid), intent(inout) :: self
        real,                       intent(in)    :: hp, lp
        if (.not. self%existence) then
            THROW_HARD('not instantiated; simple_motion_align_hybrid: set_hp_lp')
        end if
        self%hp = hp
        self%lp = lp
    end subroutine set_hp_lp

    subroutine set_trs( self, trs )
        class(motion_align_hybrid), intent(inout) :: self
        real,                    intent(in)    :: trs
        self%trs = nint(trs)
    end subroutine set_trs

    function get_iter( self ) result( iter )
        class(motion_align_hybrid), intent(inout) :: self
        integer :: iter
        iter = self%iter
    end function get_iter

    subroutine set_weights( self, frameweights )
        class(motion_align_hybrid), intent(inout) :: self
        real, allocatable,       intent(in)    :: frameweights(:)
        if (size(frameweights) /= self%nframes) then
            THROW_HARD('inconsistency; simple_motion_align_hybrid: set_weights')
        end if
        self%frameweights(:) = frameweights(:)
    end subroutine set_weights

    subroutine get_weights( self, frameweights )
        class(motion_align_hybrid), intent(inout) :: self
        real, allocatable,          intent(out)   :: frameweights(:)
        allocate(frameweights(self%nframes), source=self%frameweights)
    end subroutine get_weights

    subroutine set_rand_init_shifts( self, rand_init_shifts )
        class(motion_align_hybrid), intent(inout) :: self
        logical,                    intent(in)    :: rand_init_shifts
        self%rand_init_shifts = rand_init_shifts
    end subroutine set_rand_init_shifts

    real function get_corr( self )
        class(motion_align_hybrid), intent(in) :: self
        get_corr = self%corr
    end function get_corr

    subroutine get_corrs( self, corrs )
        class(motion_align_hybrid), intent(inout) :: self
        real, allocatable,          intent(out)   :: corrs(:)
        allocate( corrs(self%nframes), source=self%corrs )
    end subroutine get_corrs

    subroutine get_opt_shifts( self, opt_shifts )
        class(motion_align_hybrid), intent(inout) :: self
        real, allocatable,          intent(out)   :: opt_shifts(:,:)
        allocate( opt_shifts(self%nframes, 2), source=self%opt_shifts )
    end subroutine get_opt_shifts

    subroutine get_shifts_toplot( self, shifts_toplot )
        class(motion_align_hybrid), intent(inout) :: self
        real, allocatable,          intent(out)   :: shifts_toplot(:,:)
        allocate( shifts_toplot(self%nframes, 2), source=self%shifts_toplot )
    end subroutine get_shifts_toplot

    subroutine set_fixed_frame( self, fixed_frame )
        class(motion_align_hybrid), intent(inout) :: self
        integer,                    intent(in)    :: fixed_frame
        self%fixed_frame = fixed_frame
    end subroutine set_fixed_frame

    subroutine set_fitshifts( self, fitshifts )
        class(motion_align_hybrid), intent(inout) :: self
        logical,                 intent(in)    :: fitshifts
        self%fitshifts = fitshifts
    end subroutine set_fitshifts

    subroutine set_shsrch_tol( self, shsrch_tol )
        class(motion_align_hybrid), intent(inout) :: self
        real, intent(in) :: shsrch_tol
        self%shsrch_tol = shsrch_tol
    end subroutine set_shsrch_tol

    subroutine set_maxits( self, maxits )
        class(motion_align_hybrid), intent(inout) :: self
        integer, intent(in) :: maxits
        !self%maxits = maxits
    end subroutine set_maxits

    subroutine set_coords( self, x, y )
        class(motion_align_hybrid), intent(inout) :: self
        integer, intent(in) :: x, y
        self%px = x
        self%py = y
    end subroutine set_coords

    subroutine get_coords(self, x, y )
        class(motion_align_hybrid), intent(inout) :: self
        integer, intent(out) :: x, y
        x = self%px
        y = self%py
    end subroutine get_coords

    logical function is_fitshifts( self )
        class(motion_align_hybrid), intent(in) :: self
        is_fitshifts = self%fitshifts
    end function is_fitshifts

    subroutine set_bfactor( self, bfac )
        class(motion_align_hybrid), intent(inout) :: self
        real,                       intent(in)    :: bfac
        self%l_bfac = bfac > 1.e-6
        if( self%l_bfac ) self%bfactor = bfac
    end subroutine set_bfactor

    subroutine set_group_frames( self, group_frames )
        class(motion_align_hybrid), intent(inout) :: self
        logical,                    intent(in)    :: group_frames
        self%group_frames = group_frames
    end subroutine set_group_frames

    ! FITTING RELATED

    ! fit shifts to 2 polynomials
    subroutine fit_polynomial( self, shifts )
        class(motion_align_hybrid), intent(inout) :: self
        real     :: shifts(self%nframes,2)
        real(dp) :: x(self%nframes), y(self%nframes,2), sig(self%nframes)
        real(dp) :: v(POLYDIM,POLYDIM), w(POLYDIM), chisq
        integer  :: iframe
        x      = (/(real(iframe-self%fixed_frame,dp),iframe=1,self%nframes)/)
        y(:,1) = real(shifts(:,1),dp)
        y(:,2) = real(shifts(:,2),dp)
        sig    = 1.d0
        call svdfit(x, y(:,1), sig, self%polyx, v, w, chisq, poly)
        call svdfit(x, y(:,2), sig, self%polyy, v, w, chisq, poly)
    end subroutine fit_polynomial

    ! evaluate fitted shift
    subroutine polynomial2shift( self, t, shvec )
        class(motion_align_hybrid), intent(inout) :: self
        integer,                    intent(in)    :: t
        real,                       intent(out)   :: shvec(2)
        real(dp) :: rt
        rt = real(t-self%fixed_frame,dp)
        shvec(1) = poly2val(self%polyx,rt)
        shvec(2) = poly2val(self%polyy,rt)
        contains
            real function poly2val(p,x)
                real(dp), intent(in) :: p(POLYDIM), x
                poly2val = real(dot_product(p, [1.d0, x, x*x, x*x*x]))
            end function poly2val
    end subroutine polynomial2shift

    function poly(p,n) result( res )
        real(dp), intent(in) :: p
        integer,  intent(in) :: n
        real(dp) :: res(n)
        res = [1.d0, p, p*p, p*p*p]
    end function poly

    ! Destructor
    subroutine kill( self )
        class(motion_align_hybrid), intent(inout) :: self
        nullify(self%frames)
        self%bfactor = -1.
        self%smpd    = 0.
        self%ldim    = [0,0,0]
        call self%dealloc_images
        call self%dealloc_ftexps
        if( allocated(self%opt_shifts) )    deallocate(self%opt_shifts)
        if( allocated(self%shifts_toplot) ) deallocate(self%shifts_toplot)
        if( allocated(self%corrs) )         deallocate(self%corrs)
        if( allocated(self%frameweights) )  deallocate(self%frameweights)
        self%l_bfac    = .false.
        self%fitshifts = .false.
        self%existence = .false.
    end subroutine kill

end module simple_motion_align_hybrid
