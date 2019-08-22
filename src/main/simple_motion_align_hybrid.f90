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
integer, parameter :: MAXITS_PHCORR = 15
integer, parameter :: MAXITS_CORR   = 5
integer, parameter :: POLYDIM       = 4
integer, parameter :: NRESUPDATES   = 3

type :: motion_align_hybrid
    private
    type(image),           pointer :: frames_orig(:)                    !< pointer to stack of frames
    type(image),       allocatable :: frames(:)                         !< cropped frames
    type(image),       allocatable :: frames_sh(:)                      !< shifted cropped frames
    type(image)                    :: reference                         !< reference image
    real,              allocatable :: weights(:,:)                      !< weight matrix (b-factor*band-pass)
    logical,           allocatable :: resmask(:,:)                      !< band-pass mask for correlation calculation
    type(ft_expanded), allocatable :: frames_ftexp(:)                   !< ft expanded of frames
    type(ft_expanded), allocatable :: frames_ftexp_sh(:)                !< ft expanded of shifted frames
    type(ft_expanded), allocatable :: references_ftexp(:)               !< ft expanded of references
    real,              allocatable :: shifts_toplot(:,:)                !< for plotting
    real,              allocatable :: opt_shifts(:,:)                   !< shifts identified
    real,              allocatable :: frameweights(:)                   !< array of frameweights
    real,              allocatable :: corrs(:)                          !< per-frame correlations
    real(dp)                       :: polyx(POLYDIM), polyy(POLYDIM)    !< polynomial coefficients
    real                           :: ftol=1.e-6,  gtol=1.e-6           !< tolerance parameters for minimizer
    real                           :: hp=-1.,      lp=-1.               !< high/low pass value
    real                           :: lpstart=-1., lpstop=-1.           !< resolutions limits
    real                           :: bfactor = -1., bfactor_sc = 1.    !< b-factor for alignment weights
    real                           :: resstep        = 0.               !< resolution step
    real                           :: corr           = -1.              !< correlation
    real                           :: shsrch_tol     = SRCH_TOL         !< tolerance parameter for continuous srch update
    real                           :: smpd = 0., smpd_sc =0.            !< sampling distance
    real                           :: scale_factor   = 1.
    real                           :: trs            = 30.              !< half phase-correlation search bound
    integer                        :: ldim(3) = 0, ldim_sc(3) = 0       !< frame dimensions
    integer                        :: maxits_phcorr  = MAXITS_PHCORR    !< maximum number of iterations for phase-correlation search
    integer                        :: maxits_corr    = MAXITS_CORR      !< maximum number of iterations for continuous search
    integer                        :: nframes        = 0                !< number of frames
    integer                        :: fixed_frame    = 1                !< fixed (non-shifted) frame
    integer                        :: px=0 , py=0                       !< patch x/y id
    integer                        :: lp_updates       = 1              !< # of resolution updates performed [0;3]
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
    procedure, private :: gen_weights
    procedure, private :: recenter_shifts
    procedure, private :: calc_rmsd
    ! Trajectory fitting related
    procedure, private :: fit_polynomial
    procedure, private :: polynomial2shift
    ! Getters & setters
    procedure          :: set_shsrch_tol
    procedure          :: set_reslims
    procedure          :: set_trs
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
        self%trs            = params_glob%trs
        self%ftol           = params_glob%motion_correctftol
        self%gtol           = params_glob%motion_correctgtol
        self%fixed_frame    = 1
        self%lp_updates     = 1
        self%maxits_phcorr  = MAXITS_PHCORR
        self%maxits_corr    = MAXITS_CORR
        self%shsrch_tol     = SRCH_TOL
        self%bfactor        = -1.
        self%nframes        =  size(frames_ptr, 1)
        if ( self%nframes < 2 ) then
            THROW_HARD('nframes < 2; simple_motion_align_hybrid: align')
        end if
        self%frames_orig => frames_ptr
        self%smpd    = self%frames_orig(1)%get_smpd()
        self%ldim    = self%frames_orig(1)%get_ldim()
        self%hp      = min((real(minval(self%ldim(1:2))) * self%smpd)/4.,2000.)
        self%lp      = params_glob%lpstart
        self%lpstart = params_glob%lpstart
        self%lpstop  = params_glob%lpstop
        self%resstep = (self%lpstart-self%lpstop) / real(NRESUPDATES-1)
        allocate(self%frames_sh(self%nframes),self%frames(self%nframes),&
                &self%shifts_toplot(self%nframes,2), self%opt_shifts(self%nframes,2),&
                &self%corrs(self%nframes), self%frameweights(self%nframes), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('new; simple_motion_align_hybrid')
        self%shifts_toplot    = 0.
        self%opt_shifts       = 0.
        self%corrs            = -1.
        self%frameweights     = 1./real(self%nframes)
        self%fitshifts        = .false.
        self%group_frames     = .false.
        self%rand_init_shifts = .false.
        self%l_bfac           = .false.
        self%existence        = .true.
    end subroutine new

    subroutine init_images( self )
        class(motion_align_hybrid), intent(inout) :: self
        integer :: cdim(3), iframe,box
        call self%dealloc_images
        ! works out dimensions for Fourier cropping
        self%bfactor_sc   = max(30.,self%bfactor*params_glob%scale**2.)
        self%scale_factor = sqrt(-log(1.e-9) / (2.*self%bfactor_sc))
        box = round2even(self%scale_factor*real(maxval(self%ldim(1:2))))
        self%scale_factor = min(1.,real(box)/real(maxval(self%ldim(1:2))))
        self%smpd_sc = self%smpd / self%scale_factor
        if( 2.*self%smpd_sc > self%lpstop )then
            self%smpd_sc = self%lpstop/2.
            self%scale_factor = min(1.,self%smpd/self%smpd_sc)
            box = round2even(self%scale_factor*real(maxval(self%ldim(1:2))))
        endif
        if( self%ldim(1) > self%ldim(2) )then
            self%ldim_sc = [box, round2even(self%scale_factor*real(self%ldim(2))), 1]
        else
            self%ldim_sc = [round2even(self%scale_factor*real(self%ldim(1))), box, 1]
        endif
        self%smpd_sc = self%smpd / self%scale_factor
        self%trs     = self%scale_factor*real(self%trs)
        ! allocate & set
        call self%reference%new(self%ldim_sc,self%smpd_sc,wthreads=.false.)
        call self%reference%zero_and_flag_ft
        cdim = self%reference%get_array_shape()
        allocate(self%weights(cdim(1),cdim(2)),self%frames(self%nframes),&
            &self%frames_sh(self%nframes),self%resmask(cdim(1),cdim(2)))
        self%weights = 0.
        self%resmask = .false.
        !$omp parallel do default(shared) private(iframe) proc_bind(close)
        do iframe = 1,self%nframes
            call self%frames_orig(iframe)%fft
            call self%frames(iframe)%new(self%ldim_sc,self%smpd,wthreads=.false.)
            call self%frames_sh(iframe)%new(self%ldim_sc,self%smpd_sc,wthreads=.false.)
            call self%frames(iframe)%zero_and_flag_ft
            call self%frames_sh(iframe)%zero_and_flag_ft
            call self%frames_orig(iframe)%clip(self%frames(iframe))
        enddo
        !$omp end parallel do
    end subroutine init_images

    subroutine dealloc_images( self )
        class(motion_align_hybrid), intent(inout) :: self
        integer :: iframe
        if(allocated(self%frames_sh))then
            do iframe=1,self%nframes
                call self%frames(iframe)%kill
                call self%frames_sh(iframe)%kill
            enddo
            deallocate(self%frames_sh,self%frames)
        endif
        call self%reference%kill
        if( allocated(self%weights) )deallocate(self%weights)
        if( allocated(self%resmask) )deallocate(self%resmask)
    end subroutine dealloc_images

    subroutine init_ftexps( self )
        class(motion_align_hybrid), intent(inout) :: self
        integer :: iframe
        call self%dealloc_ftexps
        allocate(self%frames_ftexp(self%nframes),self%frames_ftexp_sh(self%nframes),self%references_ftexp(self%nframes))
        !$omp parallel do default(shared) private(iframe) schedule(static) proc_bind(close)
        do iframe=1,self%nframes
            call self%frames_ftexp(iframe)%new(self%frames_orig(iframe), self%hp, self%lp, .true., bfac=self%bfactor)
            call self%frames_ftexp_sh(iframe)%new(self%frames_orig(iframe), self%hp, self%lp, .false.)
            call self%references_ftexp(iframe)%new(self%frames_orig(iframe), self%hp, self%lp, .false.)
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
        ! phase-correlation search
        call self%init_images

        write(logfhandle,'(A,2I3)') '>>> PERFORMING PHASE-CORRELATION OPTIMIZATION FOR PATCH',self%px,self%py
        call self%align_phcorr( ini_shifts, frameweights )
        self%opt_shifts = self%opt_shifts / self%scale_factor
        call self%dealloc_images
        ! correlation continuous search
        call self%init_ftexps
        write(logfhandle,'(A,2I3)') '>>> PERFORMING CONTINUOUS CORRELATION OPTIMIZATION FOR PATCH',self%px,self%py
        call self%align_corr( frameweights )
        call self%dealloc_ftexps
        ! the end
        self%shifts_toplot = self%opt_shifts
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
        ! resolution related
        self%lp         = self%lpstart
        self%lp_updates = 1
        ! init shifts
        self%opt_shifts = 0.
        if( present(ini_shifts) ) self%opt_shifts = ini_shifts
        if ( self%rand_init_shifts ) then
            ! random initialization
            do iframe = 1, self%nframes
                if( iframe == self%fixed_frame ) cycle
                self%opt_shifts(iframe,1) = self%opt_shifts(iframe,1) + (ran3()-.5)*SMALLSHIFT
                self%opt_shifts(iframe,2) = self%opt_shifts(iframe,2) + (ran3()-.5)*SMALLSHIFT
            end do
        end if
        ! init weights matrix
        call self%gen_weights
        ! shift frames, generate reference & calculates correlation
        call self%shift_frames_gen_ref
        ! main loop
        do iter=1,self%maxits_phcorr
            opt_shifts_prev = self%opt_shifts
            ! individual optimizations
            !$omp parallel do schedule(static) default(shared) private(iframe) proc_bind(close)
            do iframe = 1,self%nframes
                call self%calc_shifts(iframe)
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
                if( self%lp_updates > NRESUPDATES )then
                    self%lp_updates = NRESUPDATES
                    exit
                endif
                if( self%fitshifts .and. self%lp_updates <= NRESUPDATES-1 )then
                    ! optional shifts fitting
                    call self%fit_polynomial(self%opt_shifts)
                    do iframe = 1,self%nframes
                        call self%polynomial2shift(iframe, self%opt_shifts(iframe,:))
                    enddo
                    call self%recenter_shifts(self%opt_shifts)
                endif
                ! resolution & weights update
                self%lp = max(self%lp-self%resstep, params_glob%lpstop)
                call self%gen_weights
            endif
        enddo
        ! cleanup
        do iframe=1,self%nframes
            call self%frames_sh(iframe)%kill
        enddo
        deallocate(self%frames_sh)
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
        frameweights_saved = self%frameweights
        ! shift boundaries
        trs = 3.
        ! search object allocation
        do iframe=1,self%nframes
            call ftexp_srch(iframe)%new(self%references_ftexp(iframe),&
                self%frames_ftexp_sh(iframe),trs, motion_correct_ftol=self%ftol, motion_correct_gtol=self%gtol)
            call ftexp_srch(iframe)%set_shsrch_tol(self%shsrch_tol)
            ftexp_srch(iframe)%ospec%maxits = 100
        end do
        ! generate movie sum for refinement
        opt_shifts_saved = self%opt_shifts
        call self%shift_wsum_and_calc_corrs
        corr_saved = self%corr
        ! main loop
        do iter=1,self%maxits_corr
            nimproved       = 0
            opt_shifts_prev = self%opt_shifts
            corr_prev       = self%corr
            ! individual optimizations
            !$omp parallel do schedule(static) default(shared) private(iframe,cxy) proc_bind(close)&
            !$omp reduction(+:nimproved)
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
            if( iter > 1 .and. corrfrac > 0.999 .and. rmsd < 0.1 ) exit
        end do
        ! best local optimum
        self%corr          = corr_saved
        self%opt_shifts    = opt_shifts_saved
        self%frameweights  = frameweights_saved
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

    ! shifts frames, generate reference and calculates correlations
    subroutine gen_frames_group( self )
        class(motion_align_hybrid), intent(inout) :: self
        ! complex, allocatable :: cmat_sum(:,:,:)
        ! complex,     pointer :: pcmat(:,:,:)
        ! real    :: shconst(2),ds(2), w
        ! integer :: nrlims(3,2),cdim(3), iframe
        ! ! cdim   = self%reference%get_array_shape()
        ! nrlims = self%reference%loops_lims(2)
        ! shconst(1) = merge(PI/real(self%ldim(1))/2., PI/real(self%ldim(1)-1)/2., is_even(self%ldim(1)))
        ! shconst(2) = merge(PI/real(self%ldim(2))/2., PI/real(self%ldim(2)-1)/2., is_even(self%ldim(2)))
        ! do iframe = 1,self%nframes
        !     call self%frames_sh(iframe)%get_cmt_ptr(psum)
        !     call self%frames_sh(iframe-1)%get_cmt_ptr(p)
        !     dsh = self%opt_shifts(iframe) - self%opt_shifts(iframe-1)
        !
        !
        !
        ! enddo
        !
        ! contains
        !
        !     subroutine add_shifted_weighted_frame_into(pin,pout,s,w)
        !         real, intent(in) :: s(2), w
        !         complex :: comp
        !         real    :: sh_here(2), arg
        !         sh_here = sh*shconst
        !         do h = nrlims(1,1),nrlims(1,2)
        !             do k = nrlims(2,1),nrlims(2,2)
        !                 phys = self%reference%comp_addr_phys([h,k,0])
        !                 arg  = sum(real([h,k]) * sh_here)
        !                 comp = cmplx(cos(arg),sin(arg)) * pin(phys(1),phys(2),1)
        !                 pout(phys(1),phys(2),1) = pout(phys(1),phys(2),1) + w*comp
        !             enddo
        !         enddo
        !     end subroutine add_shifted_weighted_frame_into

    end subroutine gen_frames_group

    ! band-passed correlation to frame-subtracted reference
    real function calc_corr2ref( self, frame, weight )
        class(motion_align_hybrid), intent(inout) :: self
        class(image),               intent(inout) :: frame
        real,                       intent(in)    :: weight
        complex,     pointer :: pref(:,:,:), B(:,:,:)
        complex, allocatable :: A(:,:)
        real    :: num,sumrefsq,sumframesq
        integer :: cdim(3)
        cdim = self%reference%get_array_shape()
        allocate(A(cdim(1),cdim(2)))
        call self%reference%get_cmat_ptr(pref)
        call frame%get_cmat_ptr(B)
        A = pref(:,:,1) - weight*B(:,:,1)
        num        = sum(self%weights * real(A*conjg(B(:,:,1))),mask=self%resmask)
        sumrefsq   = sum(self%weights * csq(A),                 mask=self%resmask)
        sumframesq = sum(self%weights * csq(B(:,:,1)),          mask=self%resmask)
        calc_corr2ref = 0.
        if( sumrefsq > TINY .and. sumframesq > TINY ) calc_corr2ref = num / sqrt(sumrefsq*sumframesq)
    end function calc_corr2ref

    ! identifies interpolated shifts within search range, frame destroyed on exit
    subroutine calc_shifts( self, iframe )
        class(motion_align_hybrid), intent(inout) :: self
        integer,                    intent(inout) :: iframe
        complex, pointer :: pcref(:,:,:), pcmat(:,:,:)
        real    :: prev_shift(2), shift(2), corr, alpha, beta, gamma, weight
        integer :: center(2),pos(2),i,j,trs
        ! init
        prev_shift = self%opt_shifts(iframe,:)
        weight     = self%frameweights(iframe)
        trs        = min(ceiling(self%trs),minval(self%ldim_sc(1:2)/2))
        ! phase correlation
        call self%reference%get_cmat_ptr(pcref)
        call self%frames_sh(iframe)%get_cmat_ptr(pcmat)
        pcmat(:,:,1) = self%weights(:,:) * (pcref(:,:,1) - weight*pcmat(:,:,1)) * conjg(pcmat(:,:,1))
        call self%frames_sh(iframe)%ifft
        ! find peak
        center = [self%ldim_sc(1)/2+1, self%ldim_sc(2)/2+1]
        corr   = -huge(corr)
        pos    = center
        do i = center(1)-trs,center(1)+trs
            do j = center(2)-trs,center(2)+trs
                if( self%frames_sh(iframe)%get([i,j,1]) > corr )then
                    pos  = [i,j]-center
                    corr = self%frames_sh(iframe)%get([i,j,1])
                endif
            enddo
        enddo
        shift = prev_shift + real(pos)
        ! interpolate
        alpha = self%frames_sh(iframe)%get([pos(1)+center(1)-1,pos(2)+center(2),1])
        beta  = self%frames_sh(iframe)%get([pos(1)+center(1)  ,pos(2)+center(2),1])
        gamma = self%frames_sh(iframe)%get([pos(1)+center(1)+1,pos(2)+center(2),1])
        shift(1) = shift(1) + interp_peak()
        alpha = self%frames_sh(iframe)%get([pos(1)+center(1),pos(2)+center(2)-1,1])
        beta  = self%frames_sh(iframe)%get([pos(1)+center(1),pos(2)+center(2)  ,1])
        gamma = self%frames_sh(iframe)%get([pos(1)+center(1),pos(2)+center(2)+1,1])
        shift(2) = shift(2) + interp_peak()
        self%opt_shifts(iframe,:) = shift
        ! cleanup
        call self%frames_sh(iframe)%zero_and_flag_ft
        contains

            real function interp_peak()
                real :: denom
                interp_peak = 0.
                denom = alpha+gamma-2.*beta
                if( abs(denom) < TINY )return
                interp_peak = 0.5 * (alpha-gamma) / denom
            end function interp_peak

    end subroutine calc_shifts

    ! generates weights and mask matrix
    subroutine gen_weights( self )
        class(motion_align_hybrid), intent(inout) :: self
        integer, parameter :: BPWIDTH = 3
        real    :: w, bfacw, rsh, rhplim, rlplim, width, spafreqsq,spadenomsq
        integer :: phys(3),nr_lims(3,2), bphplimsq,hplimsq,bplplimsq,lplimsq, shsq, h,k, hplim,lplim
        self%weights = 0.
        self%resmask = .false.
        width   = real(BPWIDTH)
        nr_lims = self%reference%loop_lims(2)
        hplim   = max(1,calc_fourier_index(self%hp,self%ldim(1),self%smpd))
        lplim   = calc_fourier_index(self%lp,self%ldim(1),self%smpd)
        rhplim  = real(hplim)
        rlplim  = real(lplim)
        hplimsq = hplim*hplim
        lplimsq = lplim*lplim
        bphplimsq  = max(0,hplim-BPWIDTH)**2
        bplplimsq  = min(minval(nr_lims(1:2,2)),lplim+BPWIDTH)**2
        spadenomsq = real((maxval(self%ldim(1:2)))**2)
        !$omp parallel do collapse(2) schedule(static) default(shared)&
        !$omp private(h,k,shsq,phys,rsh,bfacw,w,spafreqsq) proc_bind(close)
        do h = nr_lims(1,1),nr_lims(1,2)
            do k = nr_lims(2,1),nr_lims(2,2)
                shsq = h*h+k*k
                if( shsq < bphplimsq ) cycle
                if( shsq > bplplimsq ) cycle
                phys = self%reference%comp_addr_phys([h,k,0])
                rsh  = sqrt(real(shsq))
                spafreqsq = real(shsq) / spadenomsq
                ! B-factor weight
                bfacw = min(1.,exp(-spafreqsq*self%bfactor_sc))
                ! filter weight
                w = 1.
                if( shsq < hplimsq )then
                    ! high-pass
                    w = (1. - cos(((rsh-rhplim)/width)*PI))/2.
                else if( shsq > lplimsq )then
                    ! low_pass
                    w = (cos(((rsh-(rlplim-width))/width)*PI) + 1.)/2.
                else
                    self%resmask(phys(1),phys(2)) = .true.
                endif
                self%weights(phys(1),phys(2)) = bfacw*bfacw * w*w
            end do
        end do
        !$omp end parallel do
    end subroutine gen_weights

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

    subroutine set_reslims( self, hp, lpstart, lpstop )
        class(motion_align_hybrid), intent(inout) :: self
        real,                       intent(in)    :: hp, lpstart, lpstop
        if (.not. self%existence) then
            THROW_HARD('not instantiated; simple_motion_align_hybrid: set_reslims')
        end if
        if( hp<0. .or. lpstart<0. .or. lpstop<0. .or. lpstop>lpstart .or. hp<lpstart)then
            THROW_HARD('inconsistent resolutions limits; simple_motion_align_hybrid: set_reslims')
        endif
        self%hp      = hp
        self%lpstart = max(2.*self%smpd, lpstart)
        self%lpstop  = max(2.*self%smpd, lpstop)
        self%resstep = (self%lpstart-self%lpstop) / real(NRESUPDATES-1)
        self%lp      = self%lpstart
    end subroutine set_reslims

    subroutine set_trs( self, trs )
        class(motion_align_hybrid), intent(inout) :: self
        real,                    intent(in)    :: trs
        self%trs = trs
    end subroutine set_trs

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
        nullify(self%frames_orig)
        self%ftol    =1.e-6
        self%gtol    =1.e-6
        self%hp      =-1.
        self%lp      =-1.
        self%lpstart =-1.
        self%lpstop  =-1.
        self%resstep = 0.
        self%bfactor = -1.
        self%smpd    = 0.
        self%smpd_sc = 0.
        self%ldim    = 0
        self%ldim_sc = 0
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
