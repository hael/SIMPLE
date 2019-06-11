module simple_ctf_estimate_fit
include 'simple_lib.f08'
use simple_oris,              only: oris
use simple_image,             only: image
use simple_ctf,               only: ctf
use simple_ctf_estimate_cost, only: ctf_estimate_cost
use simple_starfile_wrappers
use CPlot2D_wrapper_module

implicit none

public :: ctf_estimate_fit
private
#include "simple_local_flags.inc"

character(len=STDLEN), parameter :: SPECKIND = 'sqrt'
integer,               parameter :: NPATCH   = 5, IARES = 5, NSTEPS = 100, POLYDIM = 6

type ctf_estimate_fit
    private
    class(image),    pointer :: micrograph
    type(image), allocatable :: boxes(:,:)              ! for storing all boxes used to build power spectra
    type(image)              :: pspec_patch(NPATCH,NPATCH)   ! patches micrograph powerspec
    type(image)              :: pspec                   ! all micrograph powerspec
    type(image)              :: pspec_ctf               ! CTF powerspec
    type(image)              :: pspec_ctf_roavg         ! rotationally averaged CTF powerspec
    type(image)              :: pspec_roavg             ! rotationally averaged all micrograph powerspec
    type(image)              :: imgmsk                  ! mask image
    type(ctf)                :: tfun                    ! transfer function object
    type(ctfparams)          :: parms                   ! for storing ctf parameters
    type(ctfparams)          :: parms_patch(NPATCH,NPATCH) ! for storing patch ctf parameters
    type(ctf_estimate_cost)  :: ctf_cost_patch(NPATCH,NPATCH) ! patch optimization objects
    type(ctf_estimate_cost)  :: ctf_cost                ! optimization object for whole micrograph
    integer, allocatable     :: inds_msk(:,:)           ! indices of pixels within resolution mask
    logical, allocatable     :: cc_msk(:,:,:)           ! redundant (including Friedel symmetry) resolution mask
    real                     :: polyx(POLYDIM), polyy(POLYDIM)
    integer                  :: centers(NPATCH,NPATCH,2)
    real                     :: smpd         = 0.
    real                     :: df_lims(2)   = [0.3,5.0]! defocus range
    real                     :: df_step      = 0.05     ! defocus step for grid search
    real                     :: astigtol     = 0.05     ! tolerated astigmatism
    real                     :: hp           = 0.       ! high-pass limit
    real                     :: lp           = 0.       ! low-pass limit
    real                     :: cc_fit       = -1.
    real                     :: cc90         = -1.
    real                     :: ctfscore     = -1.
    integer                  :: box          = 0        ! box size
    integer                  :: nbox(2)      = 0        ! # boxes along x/y
    integer                  :: npatch       = NPATCH*NPATCH
    integer                  :: flims(3,2)   = 0        ! fourier dimensions
    integer                  :: ldim_box(3)  = 0        ! box logical dimensions
    integer                  :: ldim_mic(3)  = 0        ! logical dimensions
    integer                  :: npix_msk     = 0        ! # pixels in non-redudant resolution mask
    logical                  :: exists       = .false.
contains
    ! constructor
    procedure          :: new
    procedure          :: read_doc
    procedure, private :: gen_resmsk
    procedure, private :: gen_boxes
    ! getters
    procedure          :: get_ccfit
    procedure          :: get_cc90
    procedure          :: get_ctfscore
    procedure          :: get_pspec
    ! CTF fitting
    procedure          :: fit
    procedure, private :: mic2spec
    procedure, private :: grid_srch
    procedure, private :: refine
    procedure          :: fit_patches
    procedure, private :: mic2spec_patch
    procedure, private :: gen_centers
    procedure, private :: norm_pspec
    ! scoring, display & output
    procedure          :: plot_parms
    procedure          :: write_doc
    procedure          :: write_star
    procedure, private :: calc_ctfscore
    procedure          :: write_diagnostic
    procedure          :: write_diagnostic_patch
    procedure, private :: ctf2pspecimg
    ! polynomial fitting
    procedure, private :: fit_polynomial
    procedure, private :: pix2poly
    procedure          :: pix2polyvals

    ! destructor
    procedure          :: kill
end type ctf_estimate_fit

contains

    subroutine new( self, micrograph, box, parms, dfrange, resrange, astigtol_in)
        class(ctf_estimate_fit), intent(inout) :: self
        class(image), target, intent(inout) :: micrograph       !< all micrograph powerspec
        integer,              intent(in)    :: box
        class(ctfparams),     intent(in)    :: parms
        real,                 intent(in)    :: dfrange(2)  !< defocus range, [30.0,5.0] default
        real,                 intent(in)    :: resrange(2) !< resolution range, [30.0,5.0] default
        real,                 intent(in)    :: astigtol_in !< tolerated astigmatism, 0.05 microns default
        integer :: i,j
        call self%kill
        ! set constants
        self%parms%smpd         = parms%smpd
        self%parms%cs           = parms%Cs
        self%parms%kv           = parms%kV
        self%parms%fraca        = parms%fraca
        self%parms%l_phaseplate = parms%l_phaseplate
        self%micrograph => micrograph
        call self%micrograph%ifft
        self%smpd     = self%micrograph%get_smpd()
        self%ldim_mic = self%micrograph%get_ldim()
        ! power spectrum
        self%box      = box
        self%ldim_box = [self%box,self%box,1]
        call self%pspec%new(self%ldim_box, self%smpd)
        call self%pspec_roavg%new(self%ldim_box, self%smpd)
        call self%pspec_ctf%new(self%ldim_box, self%smpd)
        call self%pspec_ctf_roavg%new(self%ldim_box, self%smpd)
        self%flims = self%pspec%loop_lims(3)
        ! init patches power spectra images
        call self%gen_centers
        do i=1,NPATCH
            do j=1,NPATCH
                call self%pspec_patch(i,j)%new(self%ldim_box, self%smpd)
            enddo
        enddo
        ! generate windows
        call self%gen_boxes
        ! search related
        if( dfrange(1) < dfrange(2) )then
            self%df_lims = dfrange
            self%df_step = (self%df_lims(2) - self%df_lims(1)) / real(NSTEPS)
        else
            THROW_HARD('invalid defocus range; ctf_estimate_init')
        endif
        self%astigtol = astigtol_in
        if( resrange(1) > resrange(2) )then
            self%hp = resrange(1)
            self%lp = resrange(2)
        else
            THROW_HARD('invalid resolution range; new')
        endif
        ! construct CTF objects
        self%tfun = ctf(self%parms%smpd, self%parms%kV, self%parms%Cs, self%parms%fraca)
        ! generate correlation mask
        call self%gen_resmsk
        ! random seed
        call seed_rnd
        self%exists = .true.
    end subroutine new

    ! constructor for reading and evaluating the polynomials only
    ! with routine pix2polyvals
    subroutine read_doc( self, fname )
        class(ctf_estimate_fit), intent(inout) :: self
        character(len=*),        intent(in)    :: fname
        type(oris) :: os
        character(len=STDLEN) :: phaseplate
        integer    :: i,xdim,ydim
        if( nlines(fname) /= 3 ) THROW_HARD('Invalid document; read_doc')
        call os%new(3)
        call os%read(fname)
        self%parms%smpd    = os%get(1,'smpd')
        self%parms%cs      = os%get(1,'cs')
        self%parms%kv      = os%get(1,'kv')
        self%parms%fraca   = os%get(1,'fraca')
        self%parms%dfx     = os%get(1,'dfx')
        self%parms%dfy     = os%get(1,'dfy')
        self%parms%angast  = os%get(1,'angast')
        self%parms%fraca   = os%get(1,'fraca')
        self%parms%phshift = os%get(1,'phshift')
        phaseplate         = os%get_static(1,'phaseplate')
        self%parms%l_phaseplate = trim(phaseplate).eq.'yes'
        ! micrograph dimensions
        self%ldim_mic(1) = nint(os%get(1,'xdim'))
        self%ldim_mic(2) = nint(os%get(1,'ydim'))
        self%ldim_mic(3) = 1
        ! polynomes
        do i = 1,POLYDIM
            self%polyx(i) = os%get(2,'px'//int2str(i))
            self%polyy(i) = os%get(3,'py'//int2str(i))
        enddo
        ! clean
        call os%kill
    end subroutine read_doc

    ! stores tiled windows
    subroutine gen_boxes( self )
        class(ctf_estimate_fit), intent(inout) :: self
        type(image) :: tmp
        integer     :: xind,yind, i,j, firstx,lastx, firsty,lasty
        logical     :: outside
        call tmp%new(self%ldim_box, self%smpd)
        self%nbox(1) = floor(real(self%ldim_mic(1))/real(self%box/2))
        self%nbox(2) = floor(real(self%ldim_mic(2))/real(self%box/2))
        allocate(self%boxes(self%nbox(1),self%nbox(2)))
        firstx = 1
        lastx  = self%ldim_mic(1)-self%box+1
        firsty = 1
        lasty  = self%ldim_mic(2)-self%box+1
        do i = 1,self%nbox(1)
            do j = 1,self%nbox(2)
                xind = firstx + floor(real((i-1)*(lastx-firstx))/real(self%nbox(1)-1)) - 1
                yind = firsty + floor(real((j-1)*(lasty-firsty))/real(self%nbox(2)-1)) - 1
                call self%boxes(i,j)%new(self%ldim_box, self%smpd)
                call self%micrograph%window_slim([xind,yind],self%box,tmp,outside)
                call tmp%norm
                call tmp%zero_edgeavg
                call tmp%fft
                call tmp%ft2img(SPECKIND, self%boxes(i,j))
                call tmp%zero_and_unflag_ft
                call self%boxes(i,j)%dampen_pspec_central_cross
            enddo
        enddo
        call tmp%kill
    end subroutine gen_boxes

    ! GETTERS

    real function get_cc90(self)
        class(ctf_estimate_fit), intent(inout) :: self
        get_cc90 = self%cc90
    end function get_cc90

    real function get_ccfit(self)
        class(ctf_estimate_fit), intent(inout) :: self
        get_ccfit = self%cc_fit
    end function get_ccfit

    real function get_ctfscore(self)
        class(ctf_estimate_fit), intent(inout) :: self
        get_ctfscore = self%ctfscore
    end function get_ctfscore

    subroutine get_pspec(self, pspec_out)
        class(ctf_estimate_fit), intent(inout) :: self
        class(image),            intent(inout) :: pspec_out
        ! self%pspec may have gone through arbitrary normalization
        call pspec_out%copy(self%pspec)
    end subroutine get_pspec

    ! DOERS

    !>  Performs initial grid search & 2D refinement, calculate stats
    subroutine fit( self, parms, diagfname )
        class(ctf_estimate_fit), intent(inout) :: self
        type(ctfparams),           intent(inout) :: parms
        character(len=*),          intent(in)    :: diagfname
        type(image) :: pspec_rot90
        ! generate power spectrum from boxes
        call self%mic2spec
        ! calculate CTF quality score based on corr with 90 deg rotated
        call self%pspec%rtsq(90., 0., 0., pspec_rot90)
        self%cc90 = self%pspec%real_corr(pspec_rot90, self%cc_msk)
        call pspec_rot90%kill
        ! normalize power spectrum with respect to resolution range
        call self%norm_pspec(self%pspec)
        ! prepare rotationally averaged power spectra & CTF power spectrum
        call self%pspec%roavg(IARES, self%pspec_roavg, 180)
        call self%norm_pspec(self%pspec_roavg)
        ! 1D grid search with rotational average
        call self%grid_srch
        ! 3/4D refinement of grid solution
        call self%refine
        ! calculate CTF score diagnostic
        call self%calc_ctfscore
        parms%dfx          = self%parms%dfx
        parms%dfy          = self%parms%dfy
        parms%angast       = self%parms%angast
        parms%phshift      = self%parms%phshift
        parms%l_phaseplate = self%parms%l_phaseplate
    end subroutine fit

    !>  Performs patch based refinement
    subroutine fit_patches( self )
        class(ctf_estimate_fit), intent(inout) :: self
        real    :: limits(2,2), cc, cc_sum
        integer :: pi,pj
        limits(1,1) = max(self%df_lims(1),self%parms%dfx-0.5)
        limits(1,2) = min(self%df_lims(2),self%parms%dfx+0.5)
        limits(2,1) = max(self%df_lims(1),self%parms%dfy-0.5)
        limits(2,2) = min(self%df_lims(2),self%parms%dfy+0.5)
        ! init
        do pi=1,NPATCH
            do pj=1,NPATCH
                self%parms_patch(pi,pj)%kv      = self%parms%kv
                self%parms_patch(pi,pj)%cs      = self%parms%cs
                self%parms_patch(pi,pj)%fraca   = self%parms%fraca
                self%parms_patch(pi,pj)%smpd    = self%parms%smpd
                self%parms_patch(pi,pj)%dfx     = self%parms%dfx
                self%parms_patch(pi,pj)%dfy     = self%parms%dfy
                self%parms_patch(pi,pj)%angast  = self%parms%angast
                self%parms_patch(pi,pj)%phshift = self%parms%phshift
                self%parms_patch(pi,pj)%l_phaseplate = self%parms%l_phaseplate
            enddo
        enddo
        ! generate sqrt(power spectra)
        call self%mic2spec_patch
        ! normalize, minimize, clean
        cc_sum = 0.
        !$omp parallel do collapse(2) default(shared) private(pi,pj,cc) &
        !$omp schedule(static) proc_bind(close) reduction(+:cc_sum)
        do pi=1,NPATCH
            do pj=1,NPATCH
                call self%norm_pspec(self%pspec_patch(pi,pj))
                call self%ctf_cost_patch(pi,pj)%init(self%pspec_patch(pi,pj), self%parms_patch(pi,pj),&
                    &self%inds_msk, 2, limits, self%astigtol)
                call self%ctf_cost_patch(pi,pj)%minimize(self%parms_patch(pi,pj), cc)
                cc_sum = cc_sum + cc
                call self%ctf_cost_patch(pi,pj)%kill
            enddo
        enddo
        !$omp end parallel do
        self%cc_fit = cc_sum/real(self%npatch)
        ! polynomial fit
        call self%fit_polynomial
    end subroutine fit_patches

    !> mic2spec calculates the average powerspectrum over a micrograph
    !!          the resulting spectrum has dampened central cross and subtracted background
    subroutine mic2spec( self )
        class(ctf_estimate_fit), intent(inout) :: self
        integer     :: i,j,n
        self%pspec = 0.
        do i = 1,self%nbox(1)
            do j = 1,self%nbox(2)
                call self%pspec%add(self%boxes(i,j))
            end do
        end do
        n = product(self%nbox)
        call self%pspec%div(real(n))
        call self%pspec%dampen_pspec_central_cross
        call self%pspec%subtr_backgr(self%hp)
    end subroutine mic2spec

    subroutine gen_centers( self )
        class(ctf_estimate_fit), intent(inout) :: self
        integer :: lims_patches(NPATCH,NPATCH,2,2)
        integer :: i,j, ldim_patch(2)
        real    :: cen, dist
        ldim_patch(1) = round2even(real(self%ldim_mic(1))/real(NPATCH))
        ldim_patch(2) = round2even(real(self%ldim_mic(2))/real(NPATCH))
        ! along X
        ! limits & center first patches
        lims_patches(1,:,1,1) = 1
        lims_patches(1,:,1,2) = ldim_patch(1)
        self%centers(1,:,1)  = sum(lims_patches(1,:,1,1:2),dim=2) / 2
        ! limits & center last patches
        lims_patches(NPATCH,:,1,1) = self%ldim_mic(1)-ldim_patch(1)+1
        lims_patches(NPATCH,:,1,2) = self%ldim_mic(1)
        self%centers(NPATCH,:,1)   = sum(lims_patches(NPATCH,:,1,1:2),dim=2) / 2
        ! adjust other patch centers to be evenly spread
        dist = real(self%centers(NPATCH,1,1)-self%centers(1,1,1)+1) / real(NPATCH-1)
        do i=2,NPATCH-1
            cen = self%centers(1,1,1) + real(i-1)*dist
            lims_patches(i,:,1,1) = ceiling(cen) - ldim_patch(1)/2
            lims_patches(i,:,1,2) = lims_patches(i,:,1,1) + ldim_patch(1) - 1
            self%centers(i,:,1)   = sum(lims_patches(i,:,1,1:2),dim=2) / 2
        enddo
        ! along Y
        lims_patches(:,1,2,1) = 1
        lims_patches(:,1,2,2) = ldim_patch(2)
        self%centers(:,1,2)  = sum(lims_patches(:,1,2,1:2),dim=2) / 2
        lims_patches(:,NPATCH,2,1) = self%ldim_mic(2)-ldim_patch(2)+1
        lims_patches(:,NPATCH,2,2) = self%ldim_mic(2)
        self%centers(:,NPATCH,2)  = sum(lims_patches(:,NPATCH,2,1:2),dim=2) / 2
        dist = real(self%centers(1,NPATCH,2)-self%centers(1,1,2)+1) / real(NPATCH-1)
        do j=2,NPATCH-1
            cen = self%centers(1,1,2) + real(j-1)*dist
            lims_patches(:,j,2,1) = ceiling(cen) - ldim_patch(2)/2
            lims_patches(:,j,2,2) = lims_patches(:,j,2,1) + ldim_patch(2) - 1
            self%centers(:,j,2)  = sum(lims_patches(:,j,2,1:2),dim=2) /2
        enddo
    end subroutine gen_centers

    !> mic2spec calculates the average powerspectrum over a micrograph
    !!          the resulting spectrum has dampened central cross and subtracted background
    subroutine mic2spec_patch( self )
        class(ctf_estimate_fit), intent(inout) :: self
        real        :: dists(product(self%nbox))
        real        :: dist, dist_thresh, w, sumw
        integer     :: pi,pj, xind,yind, cnt, center_win(2), i,j, nbox, nthresh
        nbox    = product(self%nbox)
        nthresh = max(nint(real(nbox)*0.30),min(nbox,15))
        !$omp parallel do collapse(2) default(shared) schedule(static) proc_bind(close) &
        !$omp private(i,j,pi,pj,cnt,w,sumw,xind,yind,dists,dist_thresh,center_win,dist)
        do pi = 1,NPATCH
            do pj = 1,NPATCH
                self%pspec_patch(pi,pj) = 0.
                cnt = 0
                do xind=0,self%ldim_mic(1)-self%box,self%box/2
                    do yind=0,self%ldim_mic(2)-self%box,self%box/2
                        cnt = cnt+1
                        center_win = [xind, yind] + self%box/2
                        dists(cnt) = sqrt(real(sum((self%centers(pi,pj,:)-center_win)**2)))
                    end do
                end do
                call hpsort(dists)
                dist_thresh = dists(nthresh)
                sumw = 0.
                i = 0
                do xind=0,self%ldim_mic(1)-self%box,self%box/2
                    i = i+1
                    j = 0
                    do yind=0,self%ldim_mic(2)-self%box,self%box/2
                        j = j+1
                        center_win = [xind, yind] + self%box/2
                        dist = sqrt(real(sum((self%centers(pi,pj,:)-center_win)**2)))
                        if( dist > dist_thresh ) cycle
                        w    = exp(-0.25*(dist/self%box)**2.)
                        sumw = sumw + w
                        call self%pspec_patch(pi,pj)%add(self%boxes(i,j),w)
                    end do
                end do
                call self%pspec_patch(pi,pj)%div(sumw)
                call self%pspec_patch(pi,pj)%dampen_pspec_central_cross
            enddo
        enddo
        !$omp end parallel do
        do pi = 1,NPATCH
            do pj = 1,NPATCH
                call self%pspec_patch(pi,pj)%subtr_backgr(self%hp)
            enddo
        enddo
    end subroutine mic2spec_patch

    !>  \brief  Normalize to zero mean and unit variance the reference power spectrum
    !>  within the relevent resolution range
    subroutine norm_pspec( self, img )
        class(ctf_estimate_fit), intent(inout) :: self
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
            j = min(max(1,self%inds_msk(1,i)+mh+1),self%ldim_box(1))
            k = min(max(1,self%inds_msk(2,i)+mk+1),self%ldim_box(2))
            avg = avg + real(prmat(j,k,1),dp)
        enddo
        avg = avg / real(self%npix_msk,dp)
        do i=1,self%npix_msk
            j = min(max(1,self%inds_msk(1,i)+mh+1),self%ldim_box(1))
            k = min(max(1,self%inds_msk(2,i)+mk+1),self%ldim_box(2))
            val  = real(prmat(j,k,1),dp) - avg
            sdev = sdev + val*val
        enddo
        sdev = dsqrt(sdev/real(self%npix_msk,dp))
        if(sdev <= TINY) sdev = 1.d0
        prmat(1:self%box,1:self%box,1) = real((real(prmat(1:self%box,1:self%box,1),dp) - avg)/sdev)
        nullify(prmat)
    end subroutine norm_pspec

    ! builds resolution dependent mask and indices for correlation calculation
    subroutine gen_resmsk( self )
        class(ctf_estimate_fit), intent(inout) :: self
        logical, allocatable :: nr_msk(:,:,:)
        integer :: h,k, i,j, cnt, mh,mk
        ! resolution mask
        call self%imgmsk%new(self%ldim_box, self%smpd)
        call self%imgmsk%resmsk(self%hp, self%lp)
        self%cc_msk = self%imgmsk%bin2logical()
        mh   = maxval(self%flims(1,:))
        mk   = maxval(self%flims(2,:))
        nr_msk = self%cc_msk
        ! builds mask indices
        self%npix_msk = count(nr_msk)
        allocate(self%inds_msk(2,self%npix_msk))
        cnt = 0
        do h=self%flims(1,1),self%flims(1,2)
            do k=self%flims(2,1),self%flims(2,2)
                i = min(max(1,h+mh+1),self%ldim_box(1))
                j = min(max(1,k+mk+1),self%ldim_box(2))
                if( nr_msk(i,j,1) )then
                    cnt = cnt + 1
                    self%inds_msk(:,cnt) = [h,k]
                endif
            enddo
        enddo
    end subroutine gen_resmsk

    ! calculate CTF score diagnostic
    subroutine calc_ctfscore( self )
        class(ctf_estimate_fit), intent(inout) :: self
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
        class(ctf_estimate_fit), intent(inout) :: self
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

    ! make & write half-n-half diagnostic
    subroutine write_diagnostic_patch( self, diagfname )
        class(ctf_estimate_fit), intent(inout) :: self
        character(len=*),          intent(in)    :: diagfname
        type(image) :: pspec_half_n_half, tmp
        integer     :: pi,pj
        call tmp%new(self%ldim_box,self%parms%smpd)
        self%pspec_ctf = 0.
        do pi = 1,NPATCH
            do pj = 1,NPATCH
                tmp = 0.
                if( self%parms%l_phaseplate )then
                    call self%ctf2pspecimg(tmp, self%parms_patch(pi,pj)%dfx, self%parms_patch(pi,pj)%dfy, self%parms_patch(pi,pj)%angast, add_phshift=self%parms_patch(pi,pj)%phshift)
                else
                    call self%ctf2pspecimg(tmp, self%parms_patch(pi,pj)%dfx, self%parms_patch(pi,pj)%dfy, self%parms_patch(pi,pj)%angast)
                endif
                call self%pspec_ctf%add(tmp)
            enddo
        enddo
        call self%pspec_ctf%norm()
        call self%pspec%norm()
        call self%pspec%before_after(self%pspec_ctf, pspec_half_n_half, self%cc_msk)
        call pspec_half_n_half%scale_pspec4viz
        call pspec_half_n_half%write_jpg(trim(diagfname), norm=.true.)
        call pspec_half_n_half%kill
        call tmp%kill
    end subroutine write_diagnostic_patch

    subroutine grid_srch( self )
        class(ctf_estimate_fit), intent(inout) :: self
        real :: cost, df, cost_lowest, df_best
        df          = self%df_lims(1)
        df_best     = df
        cost_lowest = huge(cost_lowest)
        call self%ctf_cost%init(self%pspec_roavg, self%parms, self%inds_msk, 1, [0.,0.], self%astigtol)
        ! do a first grid search assuming no astigmatism
        do while( df <= self%df_lims(2) )
            if( self%parms%l_phaseplate )then
                cost = self%ctf_cost%eval(df, df, 0., PIO2)
            else
                cost = self%ctf_cost%eval(df, df, 0.)
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
        call self%ctf_cost%kill
    end subroutine grid_srch

    subroutine refine( self )
        class(ctf_estimate_fit), intent(inout) :: self
        real :: limits(4,2), half_range
        ! re-init limits for local search
        half_range       = max(self%astigtol, self%df_step)
        limits      = 0.
        limits(1,1) = max(self%df_lims(1),self%parms%dfx - half_range)
        limits(2,1) = max(self%df_lims(1),self%parms%dfy - half_range)
        limits(1,2) = min(self%df_lims(2),self%parms%dfx + half_range)
        limits(2,2) = min(self%df_lims(2),self%parms%dfy + half_range)
        limits(3,1) = 0.
        limits(3,2) = 180.
        if( self%parms%l_phaseplate )then
            limits(4,1) = 0.
            limits(4,2) = 3.142
            call self%ctf_cost%init(self%pspec, self%parms, self%inds_msk, 4, limits, self%astigtol)
        else
            call self%ctf_cost%init(self%pspec, self%parms, self%inds_msk, 3, limits(1:3,:), self%astigtol)
        endif
        call self%ctf_cost%minimize(self%parms, self%cc_fit)
        call self%ctf_cost%kill
    end subroutine refine

    !>  \brief  is for making a |CTF| power-spec image
    subroutine ctf2pspecimg( self, img, dfx, dfy, angast, add_phshift )
        class(ctf_estimate_fit), intent(inout) :: self
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

    ! fit dfx/y to 2 polynomials
    subroutine fit_polynomial( self )
        class(ctf_estimate_fit), intent(inout) :: self
        real :: x(2,self%npatch), yx(self%npatch), yy(self%npatch), sig(self%npatch)
        real :: a(POLYDIM), v(POLYDIM,POLYDIM), w(POLYDIM)
        real :: dfxavg, dfyavg, chi_sq
        integer :: pi,pj, cnt
        dfxavg = 0.
        dfyavg = 0.
        cnt    = 0
        do pi=1,NPATCH
            do pj = 1,NPATCH
                cnt = cnt + 1
                dfxavg = dfxavg + self%parms_patch(pi,pj)%dfx
                dfyavg = dfyavg + self%parms_patch(pi,pj)%dfy
                call self%pix2poly(real(self%centers(pi,pj,1)), real(self%centers(pi,pj,2)), x(1,cnt), x(2,cnt))
                yx(cnt) = self%parms_patch(pi,pj)%dfx
                yy(cnt) = self%parms_patch(pi,pj)%dfy
            enddo
        enddo
        dfxavg = dfxavg / real(self%npatch)
        dfyavg = dfyavg / real(self%npatch)
        yx     = yx - dfxavg
        yy     = yy - dfyavg
        sig    = 1.
        call svd_multifit(x,yx,sig,a,v,w,chi_sq,poly)
        self%polyx(1)  = dfxavg+a(1)
        self%polyx(2:) = a(2:)
        call svd_multifit(x,yy,sig,a,v,w,chi_sq,poly)
        self%polyy(1)  = dfyavg+a(1)
        self%polyy(2:) = a(2:)
    end subroutine fit_polynomial

    function poly(p,n) result( res )
        real,    intent(in)  :: p(:)
        integer, intent(in)  :: n
        real :: res(n), x,y
        x = p(1)
        y = p(2)
        res = [1., x, x*x, y, y*y, x*y]
    end function poly

    ! real space coordinates to polynomial coordinates
    subroutine pix2poly( self, xin, yin, xout, yout )
        class(ctf_estimate_fit), intent(inout) :: self
        real,                      intent(in)  :: xin,yin
        real,                      intent(out) :: xout,yout
        xout = (xin-1.) / (self%ldim_mic(1)-1.) - 0.5
        yout = (yin-1.) / (self%ldim_mic(2)-1.) - 0.5
    end subroutine pix2poly

    ! evaluate fitted defocus
    subroutine pix2polyvals( self, xin, yin, dfx, dfy )
        class(ctf_estimate_fit), intent(inout) :: self
        real,                    intent(in)    :: xin,yin
        real,                    intent(out)   :: dfx,dfy
        real :: xp,yp
        call self%pix2poly(real(xin),real(yin), xp,yp)
        dfx = poly2val(self%polyx,xp,yp)
        dfy = poly2val(self%polyy,xp,yp)
        contains

            real function poly2val(p,x,y)
                real, intent(in) :: p(POLYDIM),x,y
                poly2val = dot_product(p, [1., x, x*x, y, y*y, x*y])
            end function poly2val

    end subroutine pix2polyvals

    subroutine plot_parms( self, fname )
        class(ctf_estimate_fit), intent(inout) :: self
        character(len=*),        intent(in)    :: fname
        real, parameter       :: SCALE = 200.
        type(str4arr)         :: title
        type(CPlot2D_type)    :: plot2D
        type(CDataSet_type)   :: calc, fit, center
        type(CDataPoint_type) :: p1
        real                  :: cx,cy,dfx,dfy,avgdfx,avgdfy,dfxmin,dfxmax,dfymin,dfymax
        integer               :: pi,pj
        avgdfx = 0.
        avgdfy = 0.
        dfxmin = huge(dfxmin)
        dfxmax = -huge(dfxmax)
        dfymin = huge(dfymin)
        dfymax = -huge(dfymax)
        do pi = 1, NPATCH
            do pj = 1, NPATCH
                avgdfx = avgdfx + self%parms_patch(pi,pj)%dfx
                avgdfy = avgdfy + self%parms_patch(pi,pj)%dfy
                dfxmin = min(dfxmin, self%parms_patch(pi,pj)%dfx)
                dfxmax = max(dfxmax, self%parms_patch(pi,pj)%dfx)
                dfymin = min(dfymin, self%parms_patch(pi,pj)%dfy)
                dfymax = max(dfymax, self%parms_patch(pi,pj)%dfy)
            enddo
        enddo
        avgdfx = avgdfx / real(self%npatch)
        avgdfy = avgdfy / real(self%npatch)
        call CPlot2D__new(plot2D, fname)
        call CPlot2D__SetDrawXAxisGridLines(plot2D, C_FALSE)
        call CPlot2D__SetDrawYAxisGridLines(plot2D, C_FALSE)
        call CPlot2D__SetXAxisSize(plot2D, 600._c_double)
        call CPlot2D__SetYAxisSize(plot2D, 600._c_double)
        call CPlot2D__SetDrawLegend(plot2D, C_FALSE)
        call CPlot2D__SetFlipY(plot2D, C_TRUE)
        do pi = 1, NPATCH
            do pj = 1, NPATCH
                ! center
                cx = real(self%centers(pi,pj,1))
                cy = real(self%centers(pi,pj,2))
                call CDataSet__new(center)
                call CDataSet__SetDrawMarker(center, C_TRUE)
                call CDataSet__SetMarkerSize(center, real(5., c_double))
                call CDataSet__SetDatasetColor(center, 1.0_c_double,0.0_c_double,0.0_c_double)
                call CDataPoint__new2(real(cx, c_double), real(cy, c_double), p1)
                call CDataSet__AddDataPoint(center, p1)
                call CDataPoint__delete(p1)
                call CPlot2D__AddDataSet(plot2D, center)
                call CDataSet__delete(center)
                ! calculated
                call CDataSet__new(calc)
                call CDataSet__SetDrawMarker(calc, C_FALSE)
                call CDataSet__SetDatasetColor(calc, 0.0_c_double,0.0_c_double,1.0_c_double)
                call CDataPoint__new2(real(cx, c_double), real(cy, c_double), p1)
                call CDataSet__AddDataPoint(calc, p1)
                call CDataPoint__delete(p1)
                dfx = SCALE * (self%parms_patch(pi,pj)%dfx-dfxmin)/(dfxmax-dfxmin)
                dfy = SCALE * (self%parms_patch(pi,pj)%dfy-dfymin)/(dfymax-dfymin)
                call CDataPoint__new2(real(cx+dfx, c_double), real(cy+dfy, c_double), p1)
                call CDataSet__AddDataPoint(calc, p1)
                call CDataPoint__delete(p1)
                call CPlot2D__AddDataSet(plot2D, calc)
                call CDataSet__delete(calc)
                ! fit
                call CDataSet__new(fit)
                call CDataSet__SetDrawMarker(fit, C_FALSE)
                call CDataSet__SetDatasetColor(fit, 0.0_c_double,0.0_c_double,0.0_c_double)
                call CDataPoint__new2(real(cx, c_double), real(cy, c_double), p1)
                call CDataSet__AddDataPoint(fit, p1)
                call CDataPoint__delete(p1)
                call self%pix2polyvals(cx,cy,dfx,dfy)
                dfx = SCALE * (dfx-dfxmin)/(dfxmax-dfxmin)
                dfy = SCALE * (dfy-dfymin)/(dfymax-dfymin)
                call CDataPoint__new2(real(cx+dfx, c_double), real(cy+dfy, c_double), p1)
                call CDataSet__AddDataPoint(fit, p1)
                call CDataPoint__delete(p1)
                call CPlot2D__AddDataSet(plot2D, fit)
                call CDataSet__delete(fit)
            end do
        end do
        title%str = 'Need to simplify. Blue: calculated; black: interpolated'//C_NULL_CHAR
        call CPlot2D__SetXAxisTitle(plot2D, title%str)
        call CPlot2D__OutputPostScriptPlot(plot2D, fname)
        call CPlot2D__delete(plot2D)
    end subroutine plot_parms

    subroutine write_doc( self, moviename, fname )
        class(ctf_estimate_fit), intent(inout) :: self
        character(len=*),        intent(in)    :: moviename, fname
        type(oris) :: os
        integer    :: i
        call os%new(3)
        call os%set(1,'smpd',   self%parms%smpd)
        call os%set(1,'cs',     self%parms%cs)
        call os%set(1,'kv',     self%parms%kv)
        call os%set(1,'fraca',  self%parms%fraca)
        call os%set(1,'dfx',    self%parms%dfx)
        call os%set(1,'dfy',    self%parms%dfy)
        call os%set(1,'angast', self%parms%angast)
        call os%set(1,'phshift',self%parms%phshift)
        call os%set(1,'forctf',    moviename)
        call os%set(1,'xdim',   real(self%ldim_mic(1)))
        call os%set(1,'ydim',   real(self%ldim_mic(2)))
        if( self%parms%l_phaseplate )then
            call os%set(1,'phaseplate','yes')
        else
            call os%set(1,'phaseplate','no')
        endif
        do i = 1, POLYDIM
            call os%set(2,'px'//int2str(i),self%polyx(i))
            call os%set(3,'py'//int2str(i),self%polyy(i))
        enddo
        call os%write(fname)
        call os%kill
    end subroutine write_doc

    subroutine write_star( self, fname )
        class(ctf_estimate_fit), intent(inout) :: self
        character(len=*),        intent(in)    :: fname
        type(starfile_table_type) :: mc_starfile
        integer :: i
        call starfile_table__new( mc_starfile )
        call starfile_table__open_ofile(mc_starfile, fname)
        call starfile_table__addObject(mc_starfile)
        call starfile_table__setIsList(mc_starfile, .true.)
        call starfile_table__setname(mc_starfile, "general")
        call starfile_table__setValue_int(mc_starfile, EMDL_IMAGE_SIZE_X, self%ldim_mic(1))
        call starfile_table__setValue_int(mc_starfile, EMDL_IMAGE_SIZE_Y, self%ldim_mic(2))
        call starfile_table__setValue_int(mc_starfile, EMDL_IMAGE_SIZE_Z, self%ldim_mic(3))
        call starfile_table__setValue_string(mc_starfile, EMDL_MICROGRAPH_NAME, simple_abspath(fname))
        call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, real(self%parms%smpd, dp))
        ! whole micrograph model
        call starfile_table__setValue_double(mc_starfile, EMDL_CTF_CS,            real(self%parms%cs, dp))
        call starfile_table__setValue_double(mc_starfile, EMDL_CTF_VOLTAGE,       real(self%parms%kv, dp))
        call starfile_table__setValue_double(mc_starfile, EMDL_CTF_Q0,            real(self%parms%fraca, dp))
        call starfile_table__setValue_double(mc_starfile, EMDL_CTF_DEFOCUSU,      real(self%parms%dfx*10000., dp))
        call starfile_table__setValue_double(mc_starfile, EMDL_CTF_DEFOCUSV,      real(self%parms%dfy*10000., dp))
        call starfile_table__setValue_double(mc_starfile, EMDL_CTF_DEFOCUS_ANGLE, real(self%parms%angast,dp))
        call starfile_table__setValue_double(mc_starfile, EMDL_CTF_PHASESHIFT,    real(self%parms%phshift,dp))
        call starfile_table__write_ofile(mc_starfile)
        call starfile_table__clear(mc_starfile)
        ! local model
        call starfile_table__setIsList(mc_starfile, .false.)
        call starfile_table__setName(mc_starfile, "local_ctf_model")
        do i = 1, POLYDIM
            call starfile_table__addObject(mc_starfile)
            call starfile_table__setValue_double(mc_starfile, EMDL_CTF_MODEL_DEFOCUSU, real(self%polyx(i)*10000.,dp))
            call starfile_table__setValue_double(mc_starfile, EMDL_CTF_MODEL_DEFOCUSV, real(self%polyy(i)*10000.,dp))
        end do
        call starfile_table__write_ofile(mc_starfile)
        call starfile_table__clear(mc_starfile)
        ! close & clean
        call starfile_table__close_ofile(mc_starfile)
        call starfile_table__delete(mc_starfile)
    end subroutine write_star

    ! DESTRUCTOR

    subroutine kill( self )
        class(ctf_estimate_fit), intent(inout) :: self
        integer :: i,j
        self%cc_fit       = -1.
        self%cc90         = -1.
        self%ctfscore     = -1.
        nullify(self%micrograph)
        call self%pspec%kill
        call self%pspec_ctf%kill
        call self%pspec_ctf_roavg%kill
        call self%pspec_roavg%kill
        call self%imgmsk%kill
        call self%ctf_cost%kill
        if( allocated(self%cc_msk) ) deallocate(self%cc_msk)
        if( allocated(self%inds_msk) ) deallocate(self%inds_msk)
        if( allocated(self%boxes) )then
            do i=1,self%nbox(1)
                do j=1,self%nbox(2)
                    call self%boxes(i,j)%kill
                enddo
            enddo
            deallocate(self%boxes)
        endif
        do i = 1,NPATCH
            do j = 1,NPATCH
                call self%pspec_patch(i,j)%kill
                call self%ctf_cost_patch(i,j)%kill
            enddo
        enddo
        self%exists = .false.
    end subroutine kill

end module simple_ctf_estimate_fit
