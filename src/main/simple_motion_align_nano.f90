module simple_motion_align_nano
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_error
use simple_image,        only: image
use simple_parameters,   only: params_glob
implicit none
public :: motion_align_nano
private
#include "simple_local_flags.inc"

integer, parameter :: MINITS_DCORR  = 3, MAXITS_DCORR  = 10
integer, parameter :: NRESUPDATES   = 2

type :: motion_align_nano
    private
    type(image), pointer     :: frames_orig(:)              !< pointer to stack of frames
    type(image), allocatable :: frames(:)                   !< cropped frames
    type(image), allocatable :: frames_sh(:)                !< shifted cropped frames
    type(image)              :: reference                   !< reference image
    real,        allocatable :: weights(:,:)                !< weight matrix (b-factor*band-pass)
    real,        allocatable :: opt_shifts(:,:)             !< shifts identified
    real,        allocatable :: frameweights(:)             !< array of frameweights
    real,        allocatable :: corrs(:)                    !< per-frame correlations
    real                     :: hp=-1., lp=-1.              !< high/low pass value
    real                     :: bfactor      = -1.          !< b-factor for alignment weights
    real                     :: resstep      = 0.           !< resolution step
    real                     :: corr         = -1.          !< correlation
    real                     :: smpd         = 0.           !< sampling distance
    real                     :: trs          = 10.          !< half correlation discrete search bound
    integer                  :: ldim(3)      = 0            !< frame dimensions
    integer                  :: maxits_dcorr = MAXITS_DCORR !< maximum number of iterations for discrete search
    integer                  :: nframes      = 0            !< number of frames
    logical                  :: l_bfac       = .false.      !< whether to use b-factor weights
    logical                  :: existence    = .false.

contains
    ! Constructor
    procedure          :: new
    ! Frames & memory management
    procedure, private :: init_images
    procedure, private :: dealloc_images
    ! Doers
    procedure, private :: calc_shifts
    procedure, private :: shift_frames_gen_ref
    procedure          :: align
    procedure, private :: align_dcorr
    procedure, private :: gen_weights
    procedure, private :: calc_rmsd
    ! Getters & setters
    procedure          :: set_reslims
    procedure          :: set_trs
    procedure          :: get_reference
    procedure          :: get_corr
    procedure          :: get_corrs
    procedure          :: get_opt_shifts
    procedure          :: set_bfactor
    ! Destructor
    procedure          :: kill
end type motion_align_nano

contains

    subroutine new( self, frames_ptr )
        class(motion_align_nano),       intent(inout) :: self
        type(image), allocatable, target, intent(in)    :: frames_ptr(:)
        call self%kill
        self%trs            = params_glob%trs
        self%maxits_dcorr   = MAXITS_DCORR
        self%bfactor        = -1.
        self%nframes        =  size(frames_ptr, 1)
        if ( self%nframes < 2 ) then
            THROW_HARD('nframes < 2; simple_motion_align_nano: align')
        end if
        self%frames_orig => frames_ptr
        self%smpd        =  self%frames_orig(1)%get_smpd()
        self%ldim        =  self%frames_orig(1)%get_ldim()
        self%hp          =  min((real(minval(self%ldim(1:2))) * self%smpd)/4.,2000.)
        self%hp          =  min(params_glob%hp, self%hp)
        self%lp          =  params_glob%lp
        allocate(self%frames_sh(self%nframes),self%frames(self%nframes),&
                &self%opt_shifts(self%nframes,2),&
                &self%corrs(self%nframes), self%frameweights(self%nframes), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('new; simple_motion_align_nano')
        self%opt_shifts       = 0.
        self%corrs        = -1.
        self%frameweights = 1./real(self%nframes)
        self%l_bfac       = .false.
        self%existence        = .true.
    end subroutine new

    subroutine init_images( self )
        class(motion_align_nano), intent(inout) :: self
        integer :: cdim(3), iframe
        call self%dealloc_images
        ! allocate & set
        call self%reference%new(self%ldim,self%smpd,wthreads=.false.)
        call self%reference%zero_and_flag_ft
        cdim = self%reference%get_array_shape()
        allocate(self%weights(cdim(1),cdim(2)),self%frames(self%nframes),&
            &self%frames_sh(self%nframes))
        self%weights = 0.
        !$omp parallel do default(shared) private(iframe) proc_bind(close)
        do iframe = 1,self%nframes
            call self%frames_orig(iframe)%fft
            call self%frames(iframe)%new(self%ldim,self%smpd,wthreads=.false.)
            call self%frames_sh(iframe)%new(self%ldim,self%smpd,wthreads=.false.)
            call self%frames(iframe)%zero_and_flag_ft
            call self%frames_sh(iframe)%zero_and_flag_ft
            call self%frames_orig(iframe)%clip(self%frames(iframe))
        enddo
        !$omp end parallel do
    end subroutine init_images

    subroutine dealloc_images( self )
        class(motion_align_nano), intent(inout) :: self
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
    end subroutine dealloc_images

    ! Alignment routines

    subroutine align( self, reference, ini_shifts )
        class(motion_align_nano), intent(inout)   :: self
        class(image),     optional, intent(inout) :: reference
        real,             optional, intent(in)    :: ini_shifts(self%nframes,2)
        if ( .not. self%existence ) then
            THROW_HARD('not instantiated; simple_motion_align_nano: align')
        end if
        if (( self%hp < 0. ) .or. ( self%lp < 0.)) then
            THROW_HARD('hp or lp < 0; simple_motion_align_nano: align')
        end if
        ! discrete correlation search
        call self%init_images
        call self%align_dcorr( reference, ini_shifts )
    end subroutine align

    ! semi-discrete correlation based alignment
    subroutine align_dcorr( self, reference, ini_shifts )
        class(motion_align_nano), intent(inout) :: self
        class(image),   optional, intent(inout) :: reference
        real,           optional, intent(in)    :: ini_shifts(self%nframes,2)
        real    :: frameweights_saved(self%nframes), opt_shifts_prev(self%nframes,2), rmsd
        integer :: iter, iframe
        logical :: l_refset
        l_refset = present(reference)
        ! init shifts & generate groups
        self%opt_shifts = 0.
        if( present(ini_shifts) ) self%opt_shifts = ini_shifts
        ! init weights matrix
        call self%gen_weights
        ! shift frames, generate reference & calculates correlation
        call self%shift_frames_gen_ref
        if( l_refset ) call self%reference%copy(reference)
        ! main loop
        do iter=1,self%maxits_dcorr
            opt_shifts_prev = self%opt_shifts
            ! individual optimizations
            if( l_refset .and. iter==1 )then
                ! setting weight to zero if reference provided on the first iteration only
                ! to circumvent frame subtraction from reference
                frameweights_saved = self%frameweights
                self%frameweights = 0.
            endif
            !$omp parallel do schedule(static) default(shared) private(iframe) proc_bind(close)
            do iframe = 1,self%nframes
                call self%calc_shifts(iframe)
            end do
            !$omp end parallel do
            if( l_refset .and. iter==1) self%frameweights = frameweights_saved
            ! shift frames, generate reference & calculates correlations
            call self%shift_frames_gen_ref
            ! convergence
            rmsd = self%calc_rmsd(opt_shifts_prev, self%opt_shifts)
            if( iter > MINITS_DCORR .and. rmsd < 0.1 ) exit
        enddo
    end subroutine align_dcorr

    ! shifts frames, generate reference and calculates correlations
    subroutine shift_frames_gen_ref( self )
        class(motion_align_nano), intent(inout) :: self
        complex, allocatable :: cmat_sum(:,:,:)
        complex,     pointer :: pcmat(:,:,:)
        integer :: iframe
        cmat_sum = self%reference%get_cmat()
        cmat_sum = cmplx(0.,0.)
        !$omp parallel do default(shared) private(iframe,pcmat) proc_bind(close) schedule(static) reduction(+:cmat_sum)
        do iframe=1,self%nframes
            call self%frames_sh(iframe)%set_cmat(self%frames(iframe))
            call self%frames_sh(iframe)%shift2Dserial(-self%opt_shifts(iframe,:))
            call self%frames_sh(iframe)%get_cmat_ptr(pcmat)
            cmat_sum = cmat_sum + pcmat * self%frameweights(iframe)
        enddo
        !$omp end parallel do
        call self%reference%set_cmat(cmat_sum)
    end subroutine shift_frames_gen_ref

    subroutine calc_shifts( self, iframe )
        class(motion_align_nano), intent(inout) :: self
        integer,                  intent(inout) :: iframe
        real,    pointer :: pcorrs(:,:,:)
        complex, pointer :: pcmat(:,:,:), pcref(:,:,:)
        real     :: dshift(2),alpha,beta,gamma,weight
        integer  :: pos(2),center(2),trs
        weight = self%frameweights(iframe)
        trs    = max(1, min(floor(self%trs),minval(self%ldim(1:2)/2)))
        ! correlations
        call self%frames_sh(iframe)%get_cmat_ptr(pcmat)
        call self%reference%get_cmat_ptr(pcref)
        pcmat(:,:,1) = self%weights(:,:) * (pcref(:,:,1)- weight*pcmat(:,:,1)) * conjg(pcmat(:,:,1))
        call self%frames_sh(iframe)%ifft
        ! find peak
        call self%frames_sh(iframe)%get_rmat_ptr(pcorrs)
        center = self%ldim(1:2)/2+1
        pos    = maxloc(pcorrs(center(1)-trs:center(1)+trs, center(2)-trs:center(2)+trs, 1))-trs-1
        dshift = real(pos)
        ! interpolate
        beta  = pcorrs(pos(1)+center(1), pos(2)+center(2), 1)
        alpha = pcorrs(pos(1)+center(1)-1,pos(2)+center(2),1)
        gamma = pcorrs(pos(1)+center(1)+1,pos(2)+center(2),1)
        if( alpha<beta .and. gamma<beta ) dshift(1) = dshift(1) + interp_peak()
        alpha = pcorrs(pos(1)+center(1),pos(2)+center(2)-1,1)
        gamma = pcorrs(pos(1)+center(1),pos(2)+center(2)+1,1)
        if( alpha<beta .and. gamma<beta ) dshift(2) = dshift(2) + interp_peak()
        ! update shift
        self%opt_shifts(iframe,:) = self%opt_shifts(iframe,:) + dshift
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
        class(motion_align_nano), intent(inout) :: self
        integer, parameter :: BPWIDTH = 3
        real    :: w, bfacw, rsh, rhplim, rlplim, width, spafreqsq,spafreqh,spafreqk
        integer :: phys(3),nr_lims(3,2), bphplimsq,hplimsq,bplplimsq,lplimsq, shsq, h,k, hplim,lplim
        self%weights = 0.
        width   = real(BPWIDTH)
        nr_lims = self%reference%loop_lims(2)
        hplim     = max(1,calc_fourier_index(self%hp,minval(self%ldim(1:2)),self%smpd))
        rhplim    = real(hplim)
        hplimsq   = hplim*hplim
        bphplimsq = max(0,hplim+BPWIDTH)**2
        lplim     = calc_fourier_index(self%lp,minval(self%ldim(1:2)),self%smpd)
        rlplim    = real(lplim)
        lplimsq   = lplim*lplim
        bplplimsq = min(minval(nr_lims(1:2,2)),lplim-BPWIDTH)**2
        !$omp parallel do collapse(2) schedule(static) default(shared) proc_bind(close)&
        !$omp private(h,k,shsq,phys,rsh,bfacw,w,spafreqsq,spafreqh,spafreqk)
        do h = nr_lims(1,1),nr_lims(1,2)
            do k = nr_lims(2,1),nr_lims(2,2)
                shsq = h*h+k*k
                if( shsq < hplimsq ) cycle
                if( shsq > lplimsq ) cycle
                if( shsq == 0 )      cycle
                phys = self%reference%comp_addr_phys([h,k,0])
                ! B-factor weight
                spafreqh  = real(h) / real(self%ldim(1)) / self%smpd
                spafreqk  = real(k) / real(self%ldim(2)) / self%smpd
                spafreqsq = spafreqh*spafreqh + spafreqk*spafreqk
                bfacw     = max(0.,exp(-spafreqsq*self%bfactor/4.))
                ! filter weight
                w = 1.
                if( shsq < bphplimsq )then
                    ! high-pass
                    rsh = sqrt(real(shsq))
                    w   = 0.5 * (1.-cos(PI*(rsh-rhplim)/width))
                else if( shsq > bplplimsq )then
                    ! low_pass
                    rsh = sqrt(real(shsq))
                    w   = 0.5*(1.+cos(PI*(rsh-(rlplim-width))/width))
                endif
                self%weights(phys(1),phys(2)) = bfacw*bfacw * w*w
            end do
        end do
        !$omp end parallel do
    end subroutine gen_weights

    real function calc_rmsd( self, prev_shifts, shifts )
        class(motion_align_nano), intent(in) :: self
        real,                       intent(in) :: prev_shifts(self%nframes,2), shifts(self%nframes,2)
        integer :: iframe
        calc_rmsd = 0.
        do iframe = 1,self%nframes
            calc_rmsd = calc_rmsd + sum((shifts(iframe,:)-prev_shifts(iframe,:))**2.)
        enddo
        calc_rmsd = sqrt(calc_rmsd/real(self%nframes))
    end function calc_rmsd

    ! Getters/setters

    subroutine set_reslims( self, hp, lp )
        class(motion_align_nano), intent(inout) :: self
        real,                       intent(in)    :: hp, lp
        if (.not. self%existence) then
            THROW_HARD('not instantiated; simple_motion_align_nano: set_reslims')
        end if
        if( hp<0. .or. lp<0. .or.  hp<lp)then
            write(logfhandle,*)'HP: ',hp
            write(logfhandle,*)'LPSTART: ',lp
            THROW_HARD('inconsistent resolutions limits; simple_motion_align_nano: set_reslims')
        endif
        self%hp = hp
        self%lp = max(2.*self%smpd, lp)
    end subroutine set_reslims

    subroutine set_trs( self, trs )
        class(motion_align_nano), intent(inout) :: self
        real,                     intent(in)    :: trs
        self%trs = trs
    end subroutine set_trs

    subroutine get_reference( self, img )
        class(motion_align_nano), intent(inout) :: self
        class(image),             intent(inout) :: img
        call img%copy(self%reference)
    end subroutine get_reference

    real function get_corr( self )
        class(motion_align_nano), intent(in) :: self
        get_corr = self%corr
    end function get_corr

    subroutine get_corrs( self, corrs )
        class(motion_align_nano), intent(inout) :: self
        real, allocatable,          intent(out)   :: corrs(:)
        allocate( corrs(self%nframes), source=self%corrs )
    end subroutine get_corrs

    subroutine get_opt_shifts( self, opt_shifts )
        class(motion_align_nano), intent(inout) :: self
        real, allocatable,          intent(out)   :: opt_shifts(:,:)
        allocate( opt_shifts(self%nframes, 2), source=self%opt_shifts )
    end subroutine get_opt_shifts

    subroutine set_maxits( self, maxits )
        class(motion_align_nano), intent(inout) :: self
        integer,                    intent(in)    :: maxits
        self%maxits_dcorr = maxits
    end subroutine set_maxits

    subroutine set_bfactor( self, bfac )
        class(motion_align_nano), intent(inout) :: self
        real,                       intent(in)    :: bfac
        self%l_bfac = bfac > 1.e-6
        if( self%l_bfac ) self%bfactor = bfac
    end subroutine set_bfactor

    ! Destructor
    subroutine kill( self )
        class(motion_align_nano), intent(inout) :: self
        nullify(self%frames_orig)
        self%hp      =-1.
        self%lp      =-1.
        self%bfactor = -1.
        self%smpd    = 0.
        self%ldim    = 0
        call self%dealloc_images
        if( allocated(self%opt_shifts) )    deallocate(self%opt_shifts)
        if( allocated(self%corrs) )         deallocate(self%corrs)
        if( allocated(self%frameweights) )  deallocate(self%frameweights)
        self%l_bfac    = .false.
        self%existence = .false.
    end subroutine kill

end module simple_motion_align_nano
