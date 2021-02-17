module simple_motion_align_poly
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_error
use simple_image,           only: image
use simple_ft_expanded_dp,  only: ft_expanded_dp
use simple_parameters,      only: params_glob
use CPlot2D_wrapper_module
use simple_timer
implicit none
public :: motion_align_poly
private
#include "simple_local_flags.inc"

integer,  parameter :: POLYDIM              = 18         ! One-dimensional size of polynomial
integer,  parameter :: MAXITS_ISO           = 5          ! Maximum # of iterations for isotropic alignement
integer,  parameter :: MAXITS_ANISO         = 10         ! Maximum # of iterations for anisotropic alignement
real,     parameter :: RMSD_THRESHOLD_ISO   = 0.2        ! Convergence RMSD for isotropic alignement
real,     parameter :: RMSD_THRESHOLD_ANISO = 0.2        ! Convergence RMSD for anisotropic alignement
real,     parameter :: RMSD_POLY_THRESHOLD  = 4.         ! RMSD acceptance threshold for anisotropic alignement
real,     parameter :: MIN_ANISO_WEIGHT     = 0.5        ! Minimum distance-based weight for anisotropic alignement
real,     parameter :: RMSD_THRESHOLD_POLY  = 0.1        ! Convergence RMSD for polynomial refinement
real,     parameter :: PBOUND               = 2.0        ! Half bound for polynomial coefficient search
integer,  parameter :: LBFGSB_MAXITS = 60, MAXITS = 3    ! Iterations parameters for polynomial refinement

type :: motion_align_poly
    private
    type(image),            pointer :: frames_orig(:)                    !< pointer to frames
    type(image),        allocatable :: tiles(:,:,:), tiles_sh(:,:,:)
    type(image),        allocatable :: tmp_imgs(:)
    type(ft_expanded_dp), allocatable :: ftexp_tiles(:,:,:), ftexp_tiles_sh(:,:)
    type(ft_expanded_dp), allocatable :: ftexp_R(:), ftexp_dR(:,:)
    type(ft_expanded_dp), allocatable :: ftexp_Rhat(:), ftexp_dRhat(:,:), ftexp_Rhat2(:)
    integer,           allocatable :: lims_tiles(:,:,:,:)
    real,              allocatable :: tile_centers(:,:,:)
    real,              allocatable :: patch_pos(:,:,:)
    real(dp),          allocatable :: patch_coords(:,:,:)
    real,              allocatable :: shifts(:,:,:,:)                   !< global shifts
    real,              allocatable :: iso_shifts(:,:,:,:)               !< isotropic shifts
    real,              allocatable :: aniso_shifts(:,:,:,:)             !< anisotropic shifts
    real,              allocatable :: patch_shifts(:,:,:,:)             !< anisotropic shifts for display only
    real,              allocatable :: frameweights(:)                   !< array of frameweights
    real,              allocatable :: corrs(:)                          !< per-frame correlations
    logical,           allocatable :: resolution_mask(:,:,:)
    real(dp)                       :: fit_poly_coeffs(POLYDIM,2)        !< fitted polynomial coefficients
    real(dp)                       :: poly_coeffs(POLYDIM,2)            !< optimized polynomial coefficients
    real                           :: hp=-1.,      lp=-1.               !< high/low pass value
    real                           :: bfactor        = 0.              !< b-factor for alignment weights
    real                           :: corr           = -1.              !< correlation
    real                           :: smpd           = 0.               !< sampling distance
    real                           :: smpd_tile      = 0.               !< sampling distance
    real                           :: scale_factor   = 1.               !< local frame scaling
    real                           :: trs            = 20.              !< half correlation disrete search bound
    integer                        :: tilesz         = 768
    integer                        :: nthr           = 1
    integer                        :: nxpatch = 0, nypatch = 0
    integer                        :: ldim(3) = 0, ldim_sc(3) = 0       !< frame dimensions
    integer                        :: ldim_tile(3)   = 0
    integer                        :: nframes        = 0                !< number of frames
    integer                        :: align_frame    = 1                !< reference frame for alignement
    integer                        :: fixed_frame    = 1                !< reference frame for interpolation
    logical                        :: l_aniso_success  = .false.
    logical                        :: existence        = .false.

contains
    ! Constructor
    procedure          :: new
    ! Frames & memory management
    procedure, private :: alloc_tmp_objs
    procedure, private :: dealloc_tmp_objs
    procedure          :: gen_patches_dimensions
    procedure          :: gen_tiles
    procedure, private :: calc_correlations
    procedure, private :: interp_peak
    procedure          :: align
    procedure          :: align_iso
    procedure          :: align_aniso
    procedure          :: refine_poly
    procedure, private :: poly_refine_f, poly_refine_fdf
    procedure          :: plot_shifts
    procedure, private :: calc_shifts
    procedure, private :: center_shifts
    procedure, private :: center_shifts_to_frame
    ! Trajectory fitting related
    procedure, private :: pix2polycoords
    procedure, private :: fit_polynomial
    ! Getters & setters
    procedure          :: get_weights
    procedure          :: get_corr
    procedure          :: get_iso_shifts
    ! Destructor
    procedure          :: kill
end type motion_align_poly

type image_ptr
    real,    public, pointer :: rmat(:,:,:)
    complex, public, pointer :: cmat(:,:,:)
end type image_ptr

contains

    subroutine new( self, frames_ptr, fixed_frame )
        class(motion_align_poly),          intent(inout) :: self
        type(image), allocatable, target, intent(in)    :: frames_ptr(:)
        integer,                          intent(in)    :: fixed_frame
        call self%kill
        self%nframes     =  size(frames_ptr, 1)
        if ( self%nframes < 2 ) then
            THROW_HARD('nframes < 2; simple_motion_align_poly: align')
        end if
        self%frames_orig => frames_ptr
        self%fixed_frame = fixed_frame
        self%align_frame = nint(real(self%nframes)/2.)
        self%smpd        = self%frames_orig(1)%get_smpd()
        self%ldim        = self%frames_orig(1)%get_ldim()
        self%hp          = min((real(minval(self%ldim(1:2))) * self%smpd)/4.,800.)
        self%nthr        = params_glob%nthr
        self%bfactor     = max(0.,params_glob%bfac)
        self%lp          = params_glob%lpstop
        self%trs         = params_glob%trs
        self%nxpatch     = params_glob%nxpatch
        self%nypatch     = params_glob%nypatch
        self%l_aniso_success = .false.
        allocate(self%corrs(self%nframes), source=-1., stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('new; simple_motion_align_poly')
        self%existence        = .true.
    end subroutine new

    subroutine alloc_tmp_objs( self )
        class(motion_align_poly), intent(inout) :: self
        integer :: ithr
        allocate(self%tmp_imgs(self%nthr))
        ! $omp parallel do private(ithr) default(shared) proc_bind(close) schedule(static)
        do ithr = 1,self%nthr
            call self%tmp_imgs(ithr)%new(self%ldim_tile, self%smpd, wthreads=.false.)
        enddo
        ! $omp end do
    end subroutine alloc_tmp_objs

    subroutine gen_patches_dimensions( self )
        class(motion_align_poly), intent(inout) :: self
        integer :: i,j, tilesz
        real    :: cen, dist
        tilesz = ceiling(maxval(real(self%ldim(1:2)) / real([self%nxpatch,self%nypatch])))
        self%tilesz = find_larger_magic_box(tilesz)
        allocate(self%lims_tiles(self%nxpatch,self%nypatch,2,2),source=0)
        allocate(self%tile_centers(self%nxpatch,self%nypatch,2),&
            &self%patch_pos(self%nxpatch,self%nypatch,2),&
            &self%shifts(self%nxpatch,self%nypatch,self%nframes,2),&
            &self%iso_shifts(self%nxpatch,self%nypatch,self%nframes,2),&
            &self%patch_shifts(self%nxpatch,self%nypatch,self%nframes,2),&
            &self%aniso_shifts(self%nxpatch,self%nypatch,self%nframes,2),source=0.)
        ! along X
        ! limits & center first patches
        self%lims_tiles(1,:,1,1) = 1
        self%lims_tiles(1,:,1,2) = self%tilesz
        self%tile_centers(1,:,1) = sum(self%lims_tiles(1,:,1,1:2),dim=2) / 2.
        ! limits & center last patches
        self%lims_tiles(self%nxpatch,:,1,1) = self%ldim(1)-self%tilesz+1
        self%lims_tiles(self%nxpatch,:,1,2) = self%ldim(1)
        self%tile_centers(self%nxpatch,:,1)  = sum(self%lims_tiles(self%nxpatch,:,1,1:2),dim=2) / 2.
        ! adjust other patch centers to be evenly spread
        dist = real(self%tile_centers(self%nxpatch,1,1)-self%tile_centers(1,1,1)+1) / real(self%nxpatch-1)
        do i=2,self%nxpatch-1
            cen = self%tile_centers(1,1,1) + real(i-1)*dist
            self%lims_tiles(i,:,1,1) = ceiling(cen) - self%tilesz/2
            self%lims_tiles(i,:,1,2) = self%lims_tiles(i,:,1,1) + self%tilesz - 1
            self%tile_centers(i,:,1)  = sum(self%lims_tiles(i,:,1,1:2),dim=2) / 2.
        enddo
        ! along Y
        self%lims_tiles(:,1,2,1) = 1
        self%lims_tiles(:,1,2,2) = self%tilesz
        self%tile_centers(:,1,2) = sum(self%lims_tiles(:,1,2,1:2),dim=2) / 2.
        self%lims_tiles(:,self%nypatch,2,1) = self%ldim(2)-self%tilesz+1
        self%lims_tiles(:,self%nypatch,2,2) = self%ldim(2)
        self%tile_centers(:,self%nypatch,2)  = sum(self%lims_tiles(:,self%nypatch,2,1:2),dim=2) / 2.
        dist = real(self%tile_centers(1,self%nypatch,2)-self%tile_centers(1,1,2)+1) / real(self%nypatch-1)
        do j=2,self%nypatch-1
            cen = self%tile_centers(1,1,2) + real(j-1)*dist
            self%lims_tiles(:,j,2,1) = ceiling(cen) - self%tilesz/2
            self%lims_tiles(:,j,2,2) = self%lims_tiles(:,j,2,1) + self%tilesz - 1
            self%tile_centers(:,j,2)  = sum(self%lims_tiles(:,j,2,1:2),dim=2) /2.
        enddo
        ! patch positions
        self%patch_pos = self%tile_centers
        allocate(self%patch_coords(self%nxpatch,self%nypatch,2))
        do i=1,self%nxpatch
            do j = 1,self%nypatch
                call self%pix2polycoords(real(self%patch_pos(i,j,1),dp),real(self%patch_pos(i,j,2),dp),&
                    &self%patch_coords(i,j,1), self%patch_coords(i,j,2))
            enddo
        enddo
    end subroutine gen_patches_dimensions

    subroutine gen_tiles( self )
        class(motion_align_poly), intent(inout) :: self
        type(image)            :: tmp_imgs(self%nthr)
        type(image_ptr)        :: ptmp_imgs(self%nthr), pframes(self%nframes), ptiles(self%nthr)
        real, allocatable      :: bfac_weights(:,:)
        real          :: spafreqsq, spafreqh, spafreqk, shsq, sumsq, w, sumw, d, sigma
        integer       :: nr_lims(3,2),phys(3),cdim(3),find,magic_ldim,hplim,lplim,hplimsq,lplimsq,ithr
        integer       :: pi,pj,i,j,t,h,k
        allocate(self%tiles(self%nframes,self%nxpatch,self%nypatch),&
                &self%tiles_sh(self%nframes,self%nxpatch,self%nypatch))
        ! generates FFT images for alignement
        !$omp parallel private(i,j,t,ithr,sumsq,w,sumw,d,pi,pj,sigma) default(shared) proc_bind(close)
        !$omp do schedule(static)
        do ithr = 1,self%nthr
            call tmp_imgs(ithr)%new([self%tilesz,self%tilesz,1], self%smpd, wthreads=.false.)
        enddo
        !$omp end do nowait
        ! works out dimensions of FFTW-friendly cropped fourier images
        !$omp single
        find            = min(tmp_imgs(1)%get_nyq(), calc_fourier_index(self%lp, self%tilesz, self%smpd))
        magic_ldim      = min(self%tilesz, find_larger_magic_box(2*(find+1)) )
        self%ldim_tile = [magic_ldim, magic_ldim, 1]
        !$omp end single
        !$omp do schedule(static)
        do t = 1,self%nframes
            ithr = omp_get_thread_num() + 1
            call self%frames_orig(t)%get_rmat_ptr(pframes(t)%rmat)
            call tmp_imgs(ithr)%get_rmat_ptr(ptmp_imgs(ithr)%rmat)
            do i = 1,self%nxpatch
                do j = 1,self%nypatch
                    ! initialize to correct dimensions
                    call self%tiles(t,i,j)%new(self%ldim_tile, self%smpd, wthreads=.false.)
                    call self%tiles(t,i,j)%zero_and_flag_ft
                    ! copy real space patches
                    call tmp_imgs(ithr)%zero_and_unflag_ft
                    ptmp_imgs(ithr)%rmat(1:self%tilesz, 1:self%tilesz, 1) = pframes(t)%rmat(self%lims_tiles(i,j,1,1):self%lims_tiles(i,j,1,2),&
                        &self%lims_tiles(i,j,2,1):self%lims_tiles(i,j,2,2), 1)
                    ! tile forward transform
                    call tmp_imgs(ithr)%fft
                    ! crop to FFTW-friendly dimensions & magnitude correction
                    call tmp_imgs(ithr)%clip(self%tiles(t,i,j))
                    if( magic_ldim /= self%tilesz ) call self%tiles(t,i,j)%mul((real(self%tilesz)/real(magic_ldim))**2.)
                enddo
            enddo
        enddo
        !$omp end do nowait
        !$omp single
        self%smpd_tile = self%tiles(1,1,1)%get_smpd()
        !$omp end single
        !$omp do schedule(static)
        do ithr = 1,self%nthr
            call tmp_imgs(ithr)%kill
        enddo
        !$omp end do nowait
        !$omp single
        self%smpd_tile    = self%tiles(1,1,1)%get_smpd()
        self%scale_factor = self%smpd / self%smpd_tile
        cdim              = self%tiles(1,1,1)%get_array_shape()
        ! shift boundary
        self%trs = self%scale_factor*self%trs
        ! weight and normalize once and for all
        nr_lims = self%tiles(1,1,1)%loop_lims(2)
        hplim   = max(1,                           calc_fourier_index(self%hp, self%tilesz, self%smpd))
        lplim   = min(self%tiles(1,1,1)%get_nyq(), calc_fourier_index(self%lp, self%tilesz, self%smpd))
        hplimsq = hplim*hplim
        lplimsq = lplim*lplim
        ! generate b-factor matrix & resolution mask
        allocate(self%resolution_mask(cdim(1),cdim(2),1),source=.false.)
        allocate(bfac_weights(cdim(1),cdim(2)),source=0.0)
        do k = nr_lims(2,1),nr_lims(2,2)
            spafreqk = real(k) / real(self%tilesz) / self%smpd
            do h = nr_lims(1,1),nr_lims(1,2)
                shsq = h*h+k*k
                if( shsq >= hplimsq .and. shsq <= lplimsq )then
                    spafreqh  = real(h) / real(self%tilesz) / self%smpd
                    spafreqsq = spafreqh*spafreqh + spafreqk*spafreqk
                    phys      = self%tiles(1,1,1)%comp_addr_phys([h,k,0])
                    bfac_weights(phys(1),phys(2)) = max(0.,exp(-spafreqsq*self%bfactor/4.))
                    self%resolution_mask(phys(1),phys(2),1) = .true.
                endif
            enddo
        enddo
        !$omp end single
        ! b-factor weights & normalization
        !$omp do collapse(3) schedule(static)
        do t = 1,self%nframes
            do i = 1,self%nxpatch
                do j = 1,self%nypatch
                    ithr = omp_get_thread_num() + 1
                    ! apply b-factor weights
                    call self%tiles(t,i,j)%get_cmat_ptr(ptiles(ithr)%cmat)
                    ptiles(ithr)%cmat(:,:,1) = ptiles(ithr)%cmat(:,:,1) * bfac_weights
                    ! normalize within resolution range
                    sumsq = sum(csq_fast(ptiles(ithr)%cmat(1,:,1)), mask=self%resolution_mask(1,:,1))
                    sumsq = sumsq + 2.*sum(csq_fast(ptiles(ithr)%cmat(2:,:,1)), mask=self%resolution_mask(2:,:,1))
                    if( sumsq > TINY ) ptiles(ithr)%cmat = ptiles(ithr)%cmat / sqrt(sumsq)
                    ! shifted tiles & correlations images
                    call self%tiles_sh(t,i,j)%new(self%ldim_tile,self%smpd,wthreads=.false.)
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine gen_tiles

    ! Alignment routines

    subroutine align( self, algorithm, poly_rmsd, aniso_success, poly_coeffs, poly_coeffs_star )
        class(motion_align_poly), intent(inout) :: self
        character(len=*),        intent(in)    :: algorithm
        real,                    intent(inout) :: poly_rmsd
        logical,                 intent(inout) :: aniso_success
        real(dp),   allocatable, intent(inout) :: poly_coeffs(:)
        real(dp),                intent(out)   :: poly_coeffs_star(2*POLYDIM)
        ! real(dp) :: global_shifts(self%nxpatch,self%nypatch,self%nframes,2)
        real(dp) :: dt
        integer  :: i,j!,t
        select case(trim(algorithm))
        case('iso')
            call self%align_iso
            aniso_success    = .false.
            self%poly_coeffs = 0.d0
        case('wpatch')
            call self%align_iso
            call self%align_aniso(poly_rmsd, aniso_success)
        case('poly')
            call self%align_iso
            call self%align_aniso(poly_rmsd, aniso_success)
            if( aniso_success ) call self%refine_poly
        case DEFAULT
            THROW_HARD('Unsupported algorithm: '//trim(algorithm))
        end select
        if(allocated(poly_coeffs)) deallocate(poly_coeffs)
        allocate(poly_coeffs(2*POLYDIM), source=0.d0)
        poly_coeffs_star = 0.d0
        if( .not.aniso_success ) return !!
        ! adjusting polynomial coefficients for interpolation, after which poly_coeffs is with respect to fixed_frame
        poly_coeffs(:POLYDIM)   = self%poly_coeffs(:,1)
        poly_coeffs(POLYDIM+1:) = self%poly_coeffs(:,2)
        if( self%fixed_frame /= self%align_frame )then
            dt = real(self%fixed_frame-self%align_frame,dp)
            do i = 1,POLYDIM,3
                self%poly_coeffs(i,1)   = dot_product(poly_coeffs(  i:i+2), [1.d0, 2.d0*dt, 3.d0*dt*dt])
                self%poly_coeffs(i+1,1) = dot_product(poly_coeffs(i+1:i+2), [1.d0, 3.d0*dt])
                self%poly_coeffs(i+2,1) = poly_coeffs(i+2)
                j = i + POLYDIM
                self%poly_coeffs(i,2)   = dot_product(poly_coeffs(  j:j+2), [1.d0, 2.d0*dt, 3.d0*dt*dt])
                self%poly_coeffs(i+1,2) = dot_product(poly_coeffs(j+1:j+2), [1.d0, 3.d0*dt])
                self%poly_coeffs(i+2,2) = poly_coeffs(j+2)
            enddo
            poly_coeffs(:POLYDIM)   = self%poly_coeffs(:,1)
            poly_coeffs(POLYDIM+1:) = self%poly_coeffs(:,2)
        endif
        ! adjusting polynomial coefficients for star output
        if( self%fixed_frame /= 1 )then
            dt = real(1-self%fixed_frame,dp)
            do i = 1,2*POLYDIM,3
                poly_coeffs_star(i)   = dot_product(poly_coeffs(  i:i+2), [1.d0, 2.d0*dt, 3.d0*dt*dt])
                poly_coeffs_star(i+1) = dot_product(poly_coeffs(i+1:i+2), [1.d0, 3.d0*dt])
                poly_coeffs_star(i+2) = poly_coeffs(i+2)
            enddo
        else
            poly_coeffs_star = poly_coeffs
        endif
        !!!!!!!!!!!! DEBUG
        ! self%align_frame = self%fixed_frame
        ! call self%calc_shifts(self%align_frame, global_shifts)
        ! ! global_shifts = self%iso_shifts + global_shifts
        ! print *,'----------- ISO X'
        ! do t=1,self%nframes
        !     write(*,'(I6)',advance='no')t
        !     do i = 1,1
        !         do j = 1,1
        !             write(*,'(2F12.6)',advance='no')self%iso_shifts(i,j,t,:)
        !         enddo
        !     enddo
        !     write(*,*)
        ! enddo
        ! print *,'----------- POLY X'
        ! do t=1,self%nframes
        !     write(*,'(I6)',advance='no')t
        !     do i = 1,1
        !         do j = 1,1
        !             write(*,'(2F12.6)',advance='no')self%aniso_shifts(i,j,t,:)
        !         enddo
        !     enddo
        !     write(*,*)
        ! enddo
        ! print *,'----------- GLOBAL X'
        ! do t=1,self%nframes
        !     write(*,'(I6)',advance='no')t
        !     do i = 1,1
        !         do j = 1,1
        !             write(*,'(2F12.6)',advance='no')global_shifts(i,j,t,:)
        !         enddo
        !     enddo
        !     write(*,*)
        ! enddo
    end subroutine align

    ! Calculate all correlations images
    subroutine calc_correlations( self, shifts )
        class(motion_align_poly), intent(inout) :: self
        real,                    intent(in)    :: shifts(self%nxpatch,self%nypatch,self%nframes,2)
        type(image_ptr) :: ptiles(self%nthr), prefs(self%nthr)
        real            :: sumsq
        integer         :: i,j,t,ithr
        !$omp parallel do private(i,j,t,ithr,sumsq) default(shared) proc_bind(close) schedule(static) collapse(2)
        do i = 1,self%nxpatch
            do j = 1,self%nypatch
                ithr = omp_get_thread_num() + 1
                ! reference generation, tmp_imgs are the references
                call self%tmp_imgs(ithr)%zero_and_flag_ft
                do t = 1,self%nframes
                    call self%tiles_sh(t,i,j)%zero_and_flag_ft
                    call self%tiles(t,i,j)%shift2Dserial( -shifts(i,j,t,:), self%tiles_sh(t,i,j))
                    call self%tmp_imgs(ithr)%add(self%tiles_sh(t,i,j))
                enddo
                ! correlations matrices
                call self%tmp_imgs(ithr)%get_cmat_ptr(prefs(ithr)%cmat)
                do t = 1,self%nframes
                    call self%tiles_sh(t,i,j)%get_cmat_ptr(ptiles(ithr)%cmat)
                    sumsq =             sum(csq_fast(prefs(ithr)%cmat(1,:,1)  - ptiles(ithr)%cmat(1,:,1)),  mask=self%resolution_mask(1,:,1) )
                    sumsq =  sumsq + 2.*sum(csq_fast(prefs(ithr)%cmat(2:,:,1) - ptiles(ithr)%cmat(2:,:,1)), mask=self%resolution_mask(2:,:,1))
                    where( self%resolution_mask )
                        ptiles(ithr)%cmat = (prefs(ithr)%cmat - ptiles(ithr)%cmat) * conjg(ptiles(ithr)%cmat)
                    else where
                        ptiles(ithr)%cmat = cmplx(0.,0.)
                    end where
                    call self%tiles_sh(t,i,j)%ifft ! now tiles_sh holds the correlations
                    if( sumsq > TINY ) call self%tiles_sh(t,i,j)%mul(1./sqrt(sumsq))
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine calc_correlations

    ! Quadratic peak interpolation
    subroutine interp_peak( self, img, shift, cc )
        class(motion_align_poly), intent(in)    :: self
        class(image),            intent(inout) :: img
        real,                    intent(out)   :: shift(2), cc
        type(image_ptr) :: pimg
        real            :: alpha, beta, gamma, denom, p
        integer         :: pos(2), center, itrs
        call img%get_rmat_ptr(pimg%rmat)
        itrs   = max(1, min(floor(self%trs), self%ldim_tile(1)/2))
        center = self%ldim_tile(1)/2+1
        pos    = maxloc(pimg%rmat(center-itrs:center+itrs, center-itrs:center+itrs, 1))-itrs-1
        cc     = pimg%rmat(pos(1)+center, pos(2)+center, 1)
        shift  = real(pos)
        beta   = cc
        ! along x
        alpha = pimg%rmat(pos(1)+center-1, pos(2)+center, 1)
        gamma = pimg%rmat(pos(1)+center+1, pos(2)+center, 1)
        if( alpha<beta .and. gamma<beta )then
            denom = alpha + gamma - 2.*beta
            if( abs(denom) > TINY )then
                p        = 0.5 * (alpha-gamma) / denom
                shift(1) = shift(1) + p
                cc       = min(1., max(cc, beta-0.25*(alpha-gamma)*p))
            endif
        endif
        ! along y
        alpha = pimg%rmat(pos(1)+center, pos(2)+center-1, 1)
        gamma = pimg%rmat(pos(1)+center, pos(2)+center+1, 1)
        if( alpha<beta .and. gamma<beta )then
            denom = alpha + gamma - 2.*beta
            if( abs(denom) > TINY )then
                p        = 0.5 * (alpha-gamma) / denom
                shift(2) = shift(2) + p
                cc       = min(1., max(cc, beta-0.25*(alpha-gamma)*p))
            endif
        endif
    end subroutine interp_peak

    ! Isotropic alignement
    subroutine align_iso( self )
        class(motion_align_poly), intent(inout) :: self
        real    :: shifts(self%nxpatch,self%nypatch,self%nframes,2)
        real    :: prev_shifts(self%nxpatch,self%nypatch,self%nframes,2)
        real    :: dshift(2), sumsq, corr, rmsd, cumul_rmsd
        integer :: i,j,t, ithr, iter
        if ( .not. self%existence ) then
            THROW_HARD('not instantiated; simple_motion_align_poly: align')
        end if
        if (( self%hp < 0. ) .or. ( self%lp < 0.)) then
            THROW_HARD('hp or lp < 0; simple_motion_align_poly: align')
        end if
        write(logfhandle,'(A)') '>>> ISOTROPIC ALIGNMENT'
        shifts = 0.
        ! MAIN LOOP
        call self%alloc_tmp_objs
        do iter = 1,MAXITS_ISO
            ! book-keeping
            prev_shifts = shifts
            ! generate all 2D correlations matrices
            call self%calc_correlations( shifts )
            !$omp parallel do private(i,j,t,ithr,dshift,sumsq,corr) default(shared) proc_bind(close) schedule(static)
            do t = 1,self%nframes
                ithr = omp_get_thread_num() + 1
                ! sum over patches at time t, tmp_imgs holds the patch sum of 2D correlations
                call self%tmp_imgs(ithr)%zero_and_unflag_ft
                do i = 1,self%nxpatch
                    do j = 1,self%nypatch
                        call self%tmp_imgs(ithr)%add(self%tiles_sh(t,i,j))
                    enddo
                enddo
                call self%tmp_imgs(ithr)%div(real(self%nxpatch*self%nypatch))
                ! quadratic interpolation
                call self%interp_peak(self%tmp_imgs(ithr), dshift, self%corrs(t))
                ! shift update
                shifts(:,:,t,1) = shifts(:,:,t,1) + dshift(1)
                shifts(:,:,t,2) = shifts(:,:,t,2) + dshift(2)
            enddo
            !$omp end parallel do
            call self%center_shifts(shifts)
            ! RMSD
            rmsd       = sum((shifts-prev_shifts)**2.)
            rmsd       = sqrt(rmsd/real(self%nxpatch*self%nypatch*self%nframes)) / self%scale_factor
            cumul_rmsd = sum(shifts**2.)
            cumul_rmsd = sqrt(cumul_rmsd/real(self%nxpatch*self%nypatch*self%nframes)) / self%scale_factor
            self%corr  = sum(self%corrs)/real(self%nframes)
            ! print *,'iso iter rmsd cumul_rmsd corr: ',iter, rmsd, cumul_rmsd, self%corr
            ! convergence
            if( iter > 1 )then
                if( rmsd < RMSD_THRESHOLD_ISO )exit
            endif
        enddo
        ! frame of reference & scale
        call self%center_shifts_to_frame(shifts, self%fixed_frame)
        self%shifts       = shifts / self%scale_factor
        self%iso_shifts   = self%shifts
        self%frameweights = corrs2weights(real(self%corrs), params_glob%wcrit_enum)
        ! cleanup
        call self%dealloc_tmp_objs
    end subroutine align_iso

    ! shifts are assumed determined from iso
    subroutine align_aniso( self, poly_rmsd, success )
        class(motion_align_poly), intent(inout) :: self
        real,                    intent(out)   :: poly_rmsd
        logical,                 intent(out)   :: success
        real    :: corrs(self%nxpatch,self%nypatch,self%nframes)
        real    :: prev_corrs(self%nxpatch,self%nypatch,self%nframes)
        real    :: iso_shifts(self%nxpatch,self%nypatch,self%nframes,2)
        real    :: shifts(self%nxpatch,self%nypatch,self%nframes,2)
        real    :: aniso_shifts(self%nxpatch,self%nypatch,self%nframes,2)
        real    :: prev_shifts(self%nxpatch,self%nypatch,self%nframes,2)
        real    :: dshift(2), sumsq, rmsd, sigma, w, sumw, d, prmsd
        integer :: i,j,t,pi,pj, ithr, iter
        if( .not. self%existence ) THROW_HARD('not instantiated; simple_motion_align_poly: align')
        write(logfhandle,'(A)') '>>> ANISOTROPIC ALIGNMENT'
        call self%alloc_tmp_objs
        shifts     = self%shifts     * self%scale_factor
        iso_shifts = self%iso_shifts * self%scale_factor
        call self%center_shifts( shifts )
        ! MAIN LOOP
        do iter = 1,MAXITS_ANISO
            ! sigma update
            sigma = max(MIN_ANISO_WEIGHT, 3. / real((iter-1)**2+1))
            ! book-keeping
            prev_shifts = shifts
            prev_corrs  = corrs
            ! generate all 2D correlations matrices
            call self%calc_correlations( shifts )
            !$omp parallel private(i,j,pi,pj,t,ithr,dshift,sumsq,w,sumw,d) default(shared) proc_bind(close)
            !$omp do schedule(static) collapse(3)
            do t = 1,self%nframes                   ! time loop
                do pi = 1,self%nxpatch              ! patches loops
                    do pj = 1,self%nypatch
                        ithr = omp_get_thread_num() + 1
                        ! Sum correlations matrices with distance weight into tmp_imgs
                        sumw = 0.
                        call self%tmp_imgs(ithr)%zero_and_unflag_ft
                        do i = 1,self%nxpatch       ! tiles  loops
                            do j = 1,self%nypatch
                                d = sqrt(sum((self%patch_pos(pi,pj,:)-self%tile_centers(i,j,:))**2.))
                                d = d / real(self%tilesz) ! in 1/TILESZ unit
                                w = exp(-0.5*(d/sigma)**2.)
                                if( w < 0.001 )cycle
                                call self%tmp_imgs(ithr)%add(self%tiles_sh(t,i,j), w=w)
                                sumw = sumw + w
                            enddo
                        enddo
                        call self%tmp_imgs(ithr)%div(sumw)
                        ! quadratic interpolation
                        call self%interp_peak(self%tmp_imgs(ithr), dshift, corrs(pi,pj,t))
                        ! shift update
                        shifts(pi,pj,t,:) = shifts(pi,pj,t,:) + dshift(:)
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp single
            call self%center_shifts(shifts)
            !$omp end single
            !$omp do schedule(static)
            do t = 1,self%nframes
                self%corrs(t)       = sum(corrs(:,:,t)) / real(self%nxpatch*self%nypatch)
                iso_shifts(:,:,t,1) = sum(shifts(:,:,t,1)) / real(self%nxpatch*self%nypatch)
                iso_shifts(:,:,t,2) = sum(shifts(:,:,t,2)) / real(self%nxpatch*self%nypatch)
                aniso_shifts(:,:,t,1) = shifts(:,:,t,1) - iso_shifts(:,:,t,1)
                aniso_shifts(:,:,t,2) = shifts(:,:,t,2) - iso_shifts(:,:,t,2)
            enddo
            !$omp end do
            !$omp end parallel
            ! shifts RMSD
            rmsd = sum((shifts-prev_shifts)**2.)
            rmsd = sqrt(rmsd/real(self%nxpatch*self%nypatch*self%nframes)) / self%scale_factor
            ! deformation model RMSD
            call self%fit_polynomial(self%align_frame, aniso_shifts, poly_rmsd)
            poly_rmsd = poly_rmsd / self%scale_factor
            ! convergence
            self%corr = sum(self%corrs) / real(self%nframes)
            ! print *,'aniso iter sigma rmsd polyrmsd corr: ',iter, sigma, rmsd, poly_rmsd, self%corr
            if( poly_rmsd > RMSD_POLY_THRESHOLD  )then
                ! revert to previous solution when fitting is poor
                corrs  = prev_corrs
                shifts = prev_shifts
                exit
            endif
            if( iter > 3 )then
                if( rmsd < RMSD_THRESHOLD_ANISO ) exit
            endif
        enddo
        ! updates solution
        self%corr = sum(corrs) / real(self%nxpatch*self%nypatch*self%nframes)
        self%shifts = shifts
        if( iter == 1 )then
            ! alignemment failed on first iteration
            success = .false.
            call self%center_shifts_to_frame(self%shifts, self%fixed_frame)
            self%shifts = self%shifts / self%scale_factor
            do t = 1,self%nframes
                self%corrs(t) = sum(corrs(:,:,t)) / real(self%nxpatch*self%nypatch)
            enddo
            self%aniso_shifts    = 0.
            self%patch_shifts    = 0.
            self%fit_poly_coeffs = 0.d0
            self%poly_coeffs     = 0.d0
        else
            success = .true.
            do t = 1,self%nframes
                self%corrs(t) = sum(corrs(:,:,t)) / real(self%nxpatch*self%nypatch)
                self%iso_shifts(:,:,t,1)   = sum(shifts(:,:,t,1)) / real(self%nxpatch*self%nypatch)
                self%iso_shifts(:,:,t,2)   = sum(shifts(:,:,t,2)) / real(self%nxpatch*self%nypatch)
                self%aniso_shifts(:,:,t,1) = shifts(:,:,t,1) - self%iso_shifts(:,:,t,1)
                self%aniso_shifts(:,:,t,2) = shifts(:,:,t,2) - self%iso_shifts(:,:,t,2)
            enddo
            ! for display only
            self%patch_shifts = self%aniso_shifts / self%scale_factor
            call self%center_shifts_to_frame(self%patch_shifts, 1)
            call self%fit_polynomial(1, self%patch_shifts, prmsd)
            self%fit_poly_coeffs = self%poly_coeffs
            ! adjust reference & recover scale
            call self%center_shifts_to_frame(self%iso_shifts,   self%fixed_frame)
            call self%center_shifts_to_frame(self%shifts,       self%fixed_frame)
            call self%center_shifts_to_frame(self%aniso_shifts, self%fixed_frame)
            self%iso_shifts   = self%iso_shifts   / self%scale_factor
            self%shifts       = self%shifts       / self%scale_factor
            self%aniso_shifts = self%aniso_shifts / self%scale_factor
            call self%fit_polynomial(self%align_frame, self%aniso_shifts, prmsd)
        endif
        self%l_aniso_success = success
        ! cleanup
        call self%dealloc_tmp_objs
    end subroutine align_aniso

    ! Refinement of polynomial coefficients
    subroutine refine_poly( self )
        use simple_opt_factory, only: opt_factory
        use simple_opt_spec,    only: opt_spec
        ! use simple_opt_lbfgsb,  only: PRINT_NEVALS
        use simple_optimizer,   only: optimizer
        class(motion_align_poly), intent(inout) :: self
        type(opt_factory)         :: ofac
        type(opt_spec)            :: ospec
        class(optimizer), pointer :: nlopt
        ! DEBUG VARIABLES
        ! real(dp) :: tmp_shifts(self%nxpatch,self%nypatch,self%nframes,2)
        ! real(dp) :: poly_coeffs(POLYDIM,2), dt
        ! real     :: minrmsd
        ! integer  :: align_frame
        real(dp) :: ini_shifts_dp(self%nxpatch,self%nypatch,self%nframes,2), shifts(self%nxpatch,self%nypatch,self%nframes,2)
        real(dp) :: prev_shifts(self%nxpatch,self%nypatch,self%nframes,2)
        real     :: iso_shifts(self%nxpatch,self%nypatch,self%nframes,2), ini_shifts(self%nxpatch,self%nypatch,self%nframes,2)
        real     :: opt_lims(POLYDIM*2,2), lowest_cost, rmsd_cumul, rmsd
        integer  :: t, i, j, iter, ithr
        iso_shifts = self%iso_shifts
        call self%center_shifts_to_frame(iso_shifts, self%align_frame)
        allocate(self%ftexp_tiles(self%nframes,self%nxpatch,self%nypatch),&
                &self%ftexp_tiles_sh(self%nframes,self%nthr), self%ftexp_R(self%nthr),&
                &self%ftexp_dR(self%nthr,2), self%ftexp_Rhat(self%nthr), self%ftexp_dRhat(self%nthr,2),&
                &self%ftexp_Rhat2(self%nthr))
        !$omp parallel default(shared) private(i,j,t,ithr) proc_bind(close)
        !$omp do collapse(3) schedule(static)
        do t = 1,self%nframes
            do i = 1, self%nxpatch
                do j = 1, self%nypatch
                    call self%tiles(t,i,j)%set_smpd(self%smpd_tile) !!
                    call self%ftexp_tiles(t,i,j)%new(self%tiles(t,i,j), [self%tilesz,self%tilesz], self%hp, self%lp, .true., bfac=0.)
                    call self%ftexp_tiles(t,i,j)%shift(-real(iso_shifts(i,j,t,:),dp))
                    call self%ftexp_tiles(t,i,j)%normalize_mat
                    call self%tiles(t,i,j)%kill
                    call self%tiles_sh(t,i,j)%kill
                enddo
            enddo
        enddo
        !$omp end do nowait
        !$omp do schedule(static)
        do ithr = 1,self%nthr
            call self%ftexp_R(ithr)%copy( self%ftexp_tiles(1,1,1) )
            call self%ftexp_R(ithr)%zero
            call self%ftexp_dR(ithr,1)%copy( self%ftexp_R(ithr) )
            call self%ftexp_dR(ithr,2)%copy( self%ftexp_R(ithr) )
            call self%ftexp_Rhat(ithr)%copy( self%ftexp_R(ithr) )
            call self%ftexp_Rhat2(ithr)%copy( self%ftexp_R(ithr) )
            call self%ftexp_dRhat(ithr,1)%copy( self%ftexp_R(ithr) )
            call self%ftexp_dRhat(ithr,2)%copy( self%ftexp_R(ithr) )
            do t = 1,self%nframes
                call self%ftexp_tiles_sh(t,ithr)%copy( self%ftexp_R(ithr) )
            enddo
        enddo
        !$omp end do
        !$omp end parallel
        ! initial polynomial coefficients
        ini_shifts = self%aniso_shifts
        call self%center_shifts_to_frame(ini_shifts, self%align_frame)
        call self%fit_polynomial(self%align_frame, ini_shifts, rmsd)
        call self%calc_shifts(self%align_frame, ini_shifts_dp)
        shifts = ini_shifts_dp
        !!!!!!!!!!!!!!!!!! gradients check
        ! vec(:polydim)   = self%poly_coeffs(:,1)
        ! vec(polydim+1:) = self%poly_coeffs(:,2)
        ! call self%poly_refine_fdf(vec, cc, grad)
        ! do i=1,2*POLYDIM
        !     vec(:polydim)   = self%poly_coeffs(:,1)
        !     vec(polydim+1:) = self%poly_coeffs(:,2)
        !     vec(i) = vec(i) - 1.d-6
        !     ccm    = self%poly_refine_f(vec)
        !     vec(:polydim)   = self%poly_coeffs(:,1)
        !     vec(polydim+1:) = self%poly_coeffs(:,2)
        !     vec(i) = vec(i) + 1.d-6
        !     ccp    = self%poly_refine_f(vec)
        !     vec(:polydim)   = self%poly_coeffs(:,1)
        !     vec(polydim+1:) = self%poly_coeffs(:,2)
        !     print *,i, vec(i), grad(i), (ccp-ccm)/2.d-6, ccm, cc, ccp, (ccp-ccm)/2.d-6/grad(i)
        ! enddo
        ! stop 'gradients check'
        !!!!!!!!!!!!!!!!!!!
        ! Optimization
        ini_shifts = real(ini_shifts_dp)
        do iter = 1,MAXITS
            prev_shifts = shifts
            call minimize
            call self%calc_shifts(self%align_frame, shifts)
            ! convergence
            rmsd       = calc_rmsd(real(prev_shifts), real(shifts))
            rmsd_cumul = calc_rmsd(ini_shifts,        real(shifts))
            self%corr  = -lowest_cost / real(self%nxpatch*self%nypatch*self%nframes)
            ! print *,'poly iter rmsd rmsdcumul corr: ', iter,rmsd,rmsd_cumul,self%corr
            if( iter>=2 .and. rmsd<RMSD_THRESHOLD_POLY ) exit
        enddo
        self%aniso_shifts = real(shifts)
        call self%center_shifts_to_frame(self%aniso_shifts, self%fixed_frame)
        ! Average correlation & frame weights
        self%corr         = sum(self%corrs) / real(self%nframes)
        self%frameweights = corrs2weights(real(self%corrs), params_glob%wcrit_enum)
        ! cleanup
        !$omp parallel default(shared) private(i,j,t,ithr) proc_bind(close)
        !$omp do schedule(static)
        do t = 1,self%nframes
            do i = 1, self%nxpatch
                do j = 1, self%nypatch
                    call self%ftexp_tiles(t,i,j)%kill
                enddo
            enddo
            do ithr = 1,self%nthr
                call self%ftexp_tiles_sh(t,ithr)%kill
            enddo
        enddo
        !$omp end do nowait
        !$omp do schedule(static)
        do ithr = 1,self%nthr
            call self%ftexp_R(ithr)%kill
            call self%ftexp_Rhat(ithr)%kill
            call self%ftexp_Rhat2(ithr)%kill
            call self%ftexp_dR(ithr,1)%kill
            call self%ftexp_dR(ithr,2)%kill
            call self%ftexp_dRhat(ithr,1)%kill
            call self%ftexp_dRhat(ithr,2)%kill
        enddo
        !$omp end do
        !$omp end parallel
        deallocate(self%ftexp_tiles, self%ftexp_tiles_sh, self%ftexp_R, self%ftexp_Rhat,&
            &self%ftexp_Rhat2, self%ftexp_dR, self%ftexp_dRhat)
        contains

            subroutine minimize
                ! search limits
                opt_lims(1:POLYDIM,1) = real(self%poly_coeffs(:,1)) - PBOUND
                opt_lims(1:POLYDIM,2) = real(self%poly_coeffs(:,1)) + PBOUND
                opt_lims(POLYDIM+1:2*POLYDIM,1) = real(self%poly_coeffs(:,2)) - PBOUND
                opt_lims(POLYDIM+1:2*POLYDIM,2) = real(self%poly_coeffs(:,2)) + PBOUND
                ! init
                call ospec%specify('lbfgsb', 2*POLYDIM, ftol=1e-3, gtol=1e-3, limits=opt_lims, maxits=LBFGSB_MAXITS)
                call ospec%set_costfun_8(poly_refine_f_wrapper)
                call ospec%set_gcostfun_8(poly_refine_df_wrapper)
                call ospec%set_fdfcostfun_8(poly_refine_fdf_wrapper)
                call ofac%new(ospec, nlopt)
                ospec%x_8(:POLYDIM)   = self%poly_coeffs(:,1)
                ospec%x_8(POLYDIM+1:) = self%poly_coeffs(:,2)
                ospec%x = real(ospec%x_8)
                ! minimization
                call nlopt%minimize(ospec, self, lowest_cost)
                ! report solution
                self%poly_coeffs(:,1) = ospec%x_8(:POLYDIM  )
                self%poly_coeffs(:,2) = ospec%x_8(POLYDIM+1:)
                call nlopt%kill
                nullify(nlopt)
            end subroutine minimize

            real function calc_rmsd( sh1, sh2 )
                real, intent(in) :: sh1(self%nxpatch,self%nypatch,self%nframes,2)
                real, intent(in) :: sh2(self%nxpatch,self%nypatch,self%nframes,2)
                calc_rmsd = sum((sh1-sh2)**2.)
                calc_rmsd = sqrt(calc_rmsd/real(self%nxpatch*self%nypatch*self%nframes))
            end function calc_rmsd

    end subroutine refine_poly

    subroutine calc_shifts( self, iframe, sh )
        class(motion_align_poly), intent(in)  :: self
        integer,                 intent(in)  :: iframe
        real(dp),                intent(out) :: sh(self%nxpatch,self%nypatch,self%nframes,2)
        real(dp) :: x,y,rt
        integer  :: i,j,t
        do i = 1, self%nxpatch
            do j = 1, self%nypatch
                x = self%patch_coords(i,j,1)
                y = self%patch_coords(i,j,2)
                do t = 1, self%nframes
                    rt = real(t-iframe, dp)
                    sh(i,j,t,1) = apply_patch_poly_dp(self%poly_coeffs(:,1), x, y, rt)
                    sh(i,j,t,2) = apply_patch_poly_dp(self%poly_coeffs(:,2), x, y, rt)
                end do
            end do
        end do
    end subroutine calc_shifts

    ! center shifts with respect to average shift such that variance is minimal
    subroutine center_shifts( self, shifts )
        class(motion_align_poly), intent(in   ) :: self
        real,                    intent(inout) :: shifts(self%nxpatch, self%nypatch, self%nframes, 2)
        real    :: avg_shift(2)
        avg_shift(1) = sum(shifts(:,:,:,1)) / real(self%nxpatch*self%nypatch*self%nframes)
        avg_shift(2) = sum(shifts(:,:,:,2)) / real(self%nxpatch*self%nypatch*self%nframes)
        shifts(:,:,:,1) = shifts(:,:,:,1) - avg_shift(1)
        shifts(:,:,:,2) = shifts(:,:,:,2) - avg_shift(2)
    end subroutine center_shifts

    ! center shifts with respect to fixed_frame
    subroutine center_shifts_to_frame( self, shifts, iframe )
        class(motion_align_poly), intent(in)    :: self
        real,                    intent(inout) :: shifts(self%nxpatch,self%nypatch,self%nframes,2)
        integer,                 intent(in)    :: iframe
        integer :: t
        do t = 1,self%nframes
            if( t == iframe )cycle
            shifts(:,:,t,1) = shifts(:,:,t,1) - shifts(:,:,iframe,1)
            shifts(:,:,t,2) = shifts(:,:,t,2) - shifts(:,:,iframe,2)
        enddo
        shifts(:,:,iframe,:) = 0.
    end subroutine center_shifts_to_frame

    ! Getters/setters

    subroutine get_weights( self, frameweights )
        class(motion_align_poly), intent(inout) :: self
        real, allocatable,       intent(out)   :: frameweights(:)
        if( allocated(frameweights) )deallocate(frameweights)
        allocate(frameweights(self%nframes), source=self%frameweights)
    end subroutine get_weights

    real function get_corr( self )
        class(motion_align_poly), intent(in) :: self
        get_corr = self%corr
    end function get_corr

    subroutine get_iso_shifts( self, shifts )
        class(motion_align_poly), intent(inout) :: self
        real, allocatable,       intent(out)   :: shifts(:,:)
        integer :: t
        allocate(shifts(self%nframes,2), source=0.)
        do t=1,self%nframes
            shifts(t,1) = sum(self%iso_shifts(:,:,t,1)) / real(self%nxpatch*self%nypatch)
            shifts(t,2) = sum(self%iso_shifts(:,:,t,2)) / real(self%nxpatch*self%nypatch)
        enddo
    end subroutine get_iso_shifts

    subroutine plot_shifts( self, ofname )
        class(motion_align_poly), intent(inout) :: self
        character(len=*),        intent(inout) :: ofname
        real, parameter           :: SCALE = 40.
        type(str4arr)             :: title
        type(CPlot2D_type)        :: plot2D
        type(CDataSet_type)       :: dataSet
        character(len=LONGSTRLEN) :: ps2pdf_cmd, fname_pdf, ps2jpeg_cmd, fname_jpeg
        real(dp) :: tt
        real     :: shifts(self%nxpatch,self%nypatch,self%nframes,2), ref_shift(2)
        real     :: x,y,cx,cy
        integer  :: i, t, j, iostat
        call CPlot2D__new(plot2D, trim(ofname)//C_NULL_CHAR)
        call CPlot2D__SetXAxisSize(plot2D, 600.d0)
        call CPlot2D__SetYAxisSize(plot2D, 600.d0)
        call CPlot2D__SetDrawLegend(plot2D, C_FALSE)
        call CPlot2D__SetFlipY(plot2D, C_TRUE)
        ! center of micrograph
        cx = real(self%ldim(1)/2+1)
        cy = real(self%ldim(2)/2+1)
        call CDataSet__new(dataSet)
        call CDataSet__SetDrawMarker(dataSet, C_TRUE)
        call CDataSet__SetMarkerSize(dataSet,5.d0)
        call CDataSet__SetDatasetColor(dataSet, 1.d0,0.d0,0.d0)
        call CDataSet_addpoint(dataSet, cx, cy)
        call CPlot2D__AddDataSet(plot2D, dataSet)
        call CDataSet__delete(dataSet)
        ! "iso" shifts
        shifts = self%iso_shifts
        call self%center_shifts_to_frame(shifts, 1)
        call CDataSet__new(dataSet)
        call CDataSet__SetDrawMarker(dataSet, C_FALSE)
        call CDataSet__SetDatasetColor(dataSet, 0.d0,0.d0,1.d0)
        do t = 1,self%nframes
            x = cx + SCALE * shifts(1,1,t,1)
            y = cy + SCALE * shifts(1,1,t,2)
            call CDataSet_addpoint(dataSet, x, y)
        end do
        call CPlot2D__AddDataSet(plot2D, dataset)
        call CDataSet__delete(dataset)
        if( self%l_aniso_success )then
            ! anisotropic shifts prior to polynomial optimization
            shifts = self%patch_shifts
            call self%center_shifts_to_frame(shifts, 1)
            do i = 1, self%nxpatch
                do j = 1, self%nypatch
                    cx = self%patch_pos(i,j,1)
                    cy = self%patch_pos(i,j,2)
                    ! patch position
                    call CDataSet__new(dataSet)
                    call CDataSet__SetDrawMarker(dataSet,C_TRUE)
                    call CDataSet__SetMarkerSize(dataSet,5.0_c_double)
                    call CDataSet__SetDatasetColor(dataSet,1.d0,0.d0,0.d0)
                    call CDataSet_addpoint(dataSet, cx, cy)
                    call CPlot2D__AddDataSet(plot2D, dataSet)
                    call CDataSet__delete(dataSet)
                    ! trajectory
                    call CDataSet__new(dataset)
                    call CDataSet__SetDrawMarker(dataset, C_FALSE)
                    call CDataSet__SetDatasetColor(dataset, 0.5d0,0.5d0,0.5d0)
                    do t = 1,self%nframes
                        x = cx + SCALE * shifts(i,j,t,1)
                        y = cy + SCALE * shifts(i,j,t,2)
                        call CDataSet_addpoint(dataset, x, y)
                    enddo
                    call CPlot2D__AddDataSet(plot2D, dataset)
                    call CDataSet__delete(dataset)
                enddo
            enddo
            ! fitted polynomial shifts
            do i = 1, self%nxpatch
                do j = 1, self%nypatch
                    cx = self%patch_pos(i,j,1)
                    cy = self%patch_pos(i,j,2)
                    ! trajectory
                    call CDataSet__new(dataset)
                    call CDataSet__SetDrawMarker(dataset, C_FALSE)
                    call CDataSet__SetDatasetColor(dataset, 0.d0,1.d0,0.d0)
                    ref_shift(1) = real(apply_patch_poly_dp(self%fit_poly_coeffs(:,1),self%patch_coords(i,j,1),self%patch_coords(i,j,2),0.d0))
                    ref_shift(2) = real(apply_patch_poly_dp(self%fit_poly_coeffs(:,1),self%patch_coords(i,j,1),self%patch_coords(i,j,2),0.d0))
                    do t = 1,self%nframes
                        tt = real(t-1,dp)
                        x = real(apply_patch_poly_dp(self%fit_poly_coeffs(:,1),self%patch_coords(i,j,1),self%patch_coords(i,j,2),tt))
                        y = real(apply_patch_poly_dp(self%fit_poly_coeffs(:,2),self%patch_coords(i,j,1),self%patch_coords(i,j,2),tt))
                        x  = cx + SCALE * (x -ref_shift(1))
                        y  = cy + SCALE * (y -ref_shift(2))
                        call CDataSet_addpoint(dataset, x, y)
                    enddo
                    call CPlot2D__AddDataSet(plot2D, dataset)
                    call CDataSet__delete(dataset)
                enddo
            enddo
            ! optimized polynomial shifts
            shifts = self%aniso_shifts
            call self%center_shifts_to_frame(shifts, 1)
            do i = 1, self%nxpatch
                do j = 1, self%nypatch
                    cx = self%patch_pos(i,j,1)
                    cy = self%patch_pos(i,j,2)
                    ! trajectory
                    call CDataSet__new(dataset)
                    call CDataSet__SetDrawMarker(dataset, C_FALSE)
                    call CDataSet__SetDatasetColor(dataset, 0.d0,0.d0,0.d0)
                    do t = 1,self%nframes
                        x = cx + SCALE * shifts(i,j,t,1)
                        y = cy + SCALE * shifts(i,j,t,2)
                        call CDataSet_addpoint(dataset, x, y)
                    enddo
                    call CPlot2D__AddDataSet(plot2D, dataset)
                    call CDataSet__delete(dataset)
                enddo
            enddo
        endif
        title%str = 'X (in pixels; trajectory scaled by '//trim(int2str(nint(SCALE)))//')'//C_NULL_CHAR
        call CPlot2D__SetXAxisTitle(plot2D, title%str)
        title%str(1:1) = 'Y'
        call CPlot2D__SetYAxisTitle(plot2D, title%str)
        call CPlot2D__OutputPostScriptPlot(plot2D, trim(ofname)//C_NULL_CHAR)
        call CPlot2D__delete(plot2D)
        ! conversion to JPEG
        fname_jpeg  = trim(get_fbody(ofname,'eps'))//'.jpeg'
        ps2jpeg_cmd = 'gs -q -sDEVICE=jpeg -dJPEGQ=98 -dNOPAUSE -dBATCH -dSAFER -dDEVICEWIDTHPOINTS=760 -dDEVICEHEIGHTPOINTS=760 -sOutputFile='&
            //trim(fname_jpeg)//' '//trim(ofname)
        call exec_cmdline(trim(adjustl(ps2jpeg_cmd)), suppress_errors=.true., exitstat=iostat)
        ! conversion to PDF
        fname_pdf  = trim(get_fbody(ofname,'eps'))//'.pdf'
        ps2pdf_cmd = 'gs -q -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dDEVICEWIDTHPOINTS=760 -dDEVICEHEIGHTPOINTS=760 -sOutputFile='&
            //trim(fname_pdf)//' '//trim(ofname)
        call exec_cmdline(trim(adjustl(ps2pdf_cmd)), suppress_errors=.true., exitstat=iostat)
        ! update name
        if( iostat == 0 )then
            call del_file(ofname)
            ofname = trim(fname_pdf)
        endif
    end subroutine plot_shifts

    ! FITTING RELATED

    subroutine pix2polycoords( self, xin, yin, x, y )
        class(motion_align_poly), intent(inout) :: self
        real(dp),                intent(in)    :: xin, yin
        real(dp),                intent(inout) :: x, y
        x = (xin-1.d0) / real(self%ldim(1)-1,dp) - 0.5d0
        y = (yin-1.d0) / real(self%ldim(2)-1,dp) - 0.5d0
    end subroutine pix2polycoords

    ! Polynomial for patch motion
    function patch_poly(p, n) result(res)
        real(dp), intent(in) :: p(:)
        integer,  intent(in) :: n
        real(dp) :: res(n)
        real(dp) :: x, y, t
        x = p(1)
        y = p(2)
        t = p(3)
        res(    1) = t
        res(    2) = t*t
        res(    3) = t*t*t
        res( 4: 6) = x * res( 1: 3)  ! x   * {t,t^2,t^3}
        res( 7: 9) = x * res( 4: 6)  ! x^2 * {t,t^2,t^3}
        res(10:12) = y * res( 1: 3)  ! y   * {t,t^2,t^3}
        res(13:15) = y * res(10:12)  ! y^2 * {t,t^2,t^3}
        res(16:18) = y * res( 4: 6)  ! x*y * {t,t^2,t^3}
    end function patch_poly

    subroutine fit_polynomial( self, iframe, shifts, poly_rmsd )
        class(motion_align_poly), intent(inout) :: self
        integer,                 intent(in)    :: iframe
        real,                    intent(in)    :: shifts(self%nxpatch,self%nypatch,self%nframes,2)
        real,                    intent(out)   :: poly_rmsd
        real(dp) :: yx(self%nframes*self%nxpatch*self%nypatch)      ! along x
        real(dp) :: yy(self%nframes*self%nxpatch*self%nypatch)      ! along y
        real(dp) :: x(3,self%nframes*self%nxpatch*self%nypatch)     ! x,y,t
        real(dp) :: sig(self%nframes*self%nxpatch*self%nypatch)
        real(dp) :: tt, v(POLYDIM,POLYDIM), w(POLYDIM), chisq, shift(2)
        integer  :: idx, t, i, j
        sig = 1.d0
        ! fitting
        idx = 0
        do t = 1,self%nframes
            do i = 1,self%nxpatch
                do j = 1,self%nypatch
                    idx = idx+1
                    yx(idx)  = real(shifts(i,j,t,1) - shifts(i,j,iframe,1),dp)
                    yy(idx)  = real(shifts(i,j,t,2) - shifts(i,j,iframe,2),dp)
                    x(1,idx) = self%patch_coords(i,j,1)
                    x(2,idx) = self%patch_coords(i,j,2)
                    x(3,idx) = real(t-iframe,dp)
                end do
            end do
        end do
        call svd_multifit(x,yx,sig,self%poly_coeffs(:,1),v,w,chisq,patch_poly)
        call svd_multifit(x,yy,sig,self%poly_coeffs(:,2),v,w,chisq,patch_poly)
        ! goodness of fit
        poly_rmsd = 0.
        idx       = 0
        do t = 1,self%nframes
            tt = real(t-iframe, dp)
            do i = 1,self%nxpatch
                do j = 1,self%nypatch
                    idx = idx+1
                    shift(1) = apply_patch_poly_dp(self%poly_coeffs(:,1), self%patch_coords(i,j,1), self%patch_coords(i,j,2), tt)
                    shift(2) = apply_patch_poly_dp(self%poly_coeffs(:,2), self%patch_coords(i,j,1), self%patch_coords(i,j,2), tt)
                    poly_rmsd = poly_rmsd + real(sum(shift-[yx(idx),yy(idx)])**2.)
                end do
            end do
        end do
        poly_rmsd = sqrt(poly_rmsd/real(self%nframes*self%nxpatch*self%nypatch))
    end subroutine fit_polynomial

    pure real(dp) function apply_patch_poly_dp(c, x, y, t)
        real(dp), intent(in) :: c(POLYDIM), x, y, t
        real(dp) :: res, t2, t3
        t2 =  t * t
        t3 = t2 * t
        res =              c( 1) * t + c( 2) * t2 + c( 3) * t3
        res = res + x   * (c( 4) * t + c( 5) * t2 + c( 6) * t3)
        res = res + x*x * (c( 7) * t + c( 8) * t2 + c( 9) * t3)
        res = res + y   * (c(10) * t + c(11) * t2 + c(12) * t3)
        res = res + y*y * (c(13) * t + c(14) * t2 + c(15) * t3)
        res = res + x*y * (c(16) * t + c(17) * t2 + c(18) * t3)
        apply_patch_poly_dp = res
    end function apply_patch_poly_dp

    ! Destructor

    subroutine dealloc_tmp_objs( self )
        class(motion_align_poly), intent(inout) :: self
        integer :: ithr
        if( allocated(self%tmp_imgs ))then
            ! $omp parallel do private(ithr) default(shared) proc_bind(close) schedule(static)
            do ithr = 1,self%nthr
                call self%tmp_imgs(ithr)%kill
            enddo
            ! $omp end do
            deallocate(self%tmp_imgs)
        endif
    end subroutine dealloc_tmp_objs

    subroutine kill( self )
        class(motion_align_poly), intent(inout) :: self
        integer :: i,j,t
        nullify(self%frames_orig)
        call self%dealloc_tmp_objs
        self%hp      = -1.
        self%lp      = -1.
        self%bfactor = -1.
        self%scale_factor = 1.
        self%corr        = -1.
        self%trs         = 20.
        self%smpd        = 0.
        self%smpd_tile   = 0.
        self%nthr        = 1
        self%ldim        = 0
        self%ldim_sc     = 0
        self%fixed_frame = 1
        self%align_frame = 1
        self%fit_poly_coeffs = 0.d0
        self%poly_coeffs     = 0.d0
        if( allocated(self%tiles) )then
            !$omp parallel do private(i,j,t) default(shared) proc_bind(close) schedule(static)
            do i = 1,self%nxpatch
                do j = 1,self%nypatch
                    do t = 1,self%nframes
                        call self%tiles(t,i,j)%kill
                        call self%tiles_sh(t,i,j)%kill
                    enddo
                enddo
            enddo
            !$omp end parallel do
            deallocate(self%tiles,self%tiles_sh)
        endif
        if(allocated(self%lims_tiles))      deallocate(self%lims_tiles)
        if(allocated(self%tile_centers))    deallocate(self%tile_centers)
        if(allocated(self%patch_pos))       deallocate(self%patch_pos)
        if(allocated(self%patch_coords))    deallocate(self%patch_coords)
        if(allocated(self%resolution_mask)) deallocate(self%resolution_mask)
        if(allocated(self%shifts))          deallocate(self%shifts)
        if(allocated(self%iso_shifts))      deallocate(self%iso_shifts)
        if(allocated(self%patch_shifts))    deallocate(self%patch_shifts)
        if(allocated(self%aniso_shifts))    deallocate(self%aniso_shifts)
        if(allocated(self%corrs))           deallocate(self%corrs)
        if(allocated(self%frameweights))    deallocate(self%frameweights)
        self%existence = .false.
    end subroutine kill

    ! DIRECT CONTINUOUS POLYNOMIAL REFINEMENT

    real(dp) function poly_refine_f_wrapper( self, vec, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(in)    :: vec(D)
        poly_refine_f_wrapper = 0.d0
        select type(self)
        class is(motion_align_poly)
            poly_refine_f_wrapper = self%poly_refine_f( vec )
        class DEFAULT
            THROW_HARD('unknown type; patched_refine_f_wrapper')
        end select
    end function poly_refine_f_wrapper

    subroutine poly_refine_fdf_wrapper( self, vec, f, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: f, grad(D)
        f    = 0.d0
        grad = 0.d0
        select type(self)
        class is(motion_align_poly)
                call self%poly_refine_fdf( vec, f, grad )
        class DEFAULT
            THROW_HARD('unknown type; patched_direct_fdf_wrapper')
        end select
    end subroutine poly_refine_fdf_wrapper

    subroutine poly_refine_df_wrapper( self, vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        real(dp) :: f
        grad = 0.d0
        select type(self)
        class is(motion_align_poly)
            call self%poly_refine_fdf( vec, f, grad )
        class DEFAULT
            THROW_HARD('unknown type; patched_direct_fdf_wrapper')
        end select
    end subroutine poly_refine_df_wrapper

    real(dp) function poly_refine_f( self, vec )
        class(motion_align_poly), intent(inout) :: self
        real(dp),                intent(in)    :: vec(2*POLYDIM)
        real(dp) :: Es(self%nframes), Ds(self%nframes), ccs(self%nframes)
        real(dp) :: RR,x,y,rt,sx,sy
        integer  :: i,j,t,ithr
        ccs = 0.d0
        !$omp parallel do collapse(2) default(shared) private(i,j,t,rt,x,y,ithr,RR,Es,Ds,sx,sy)&
        !$omp proc_bind(close) schedule(static) reduction(+:ccs)
        do i = 1, self%nxpatch
            do j = 1, self%nypatch
                ithr = omp_get_thread_num() + 1
                x    = self%patch_coords(i,j,1)
                y    = self%patch_coords(i,j,2)
                call self%ftexp_R(ithr)%zero
                do t = 1,self%nframes
                    rt = real(t-self%align_frame, dp)
                    sx = apply_patch_poly_dp(vec(1:POLYDIM),  x,y, rt)
                    sy = apply_patch_poly_dp(vec(POLYDIM+1:), x,y, rt)
                    call self%ftexp_tiles(t,i,j)%shift(-[sx,sy], self%ftexp_tiles_sh(t,ithr))
                    call self%ftexp_R(ithr)%add( self%ftexp_tiles_sh(t,ithr) )
                end do
                do t = 1, self%nframes
                    Es(t) = self%ftexp_R(ithr)%corr_unnorm( self%ftexp_tiles_sh(t,ithr) )
                end do
                RR  = sum(Es)
                Ds  = sqrt(RR - 2.d0 * Es + 1.d0)
                ccs = ccs + (Es-1.d0)/Ds
            end do
        end do
        !$omp end parallel do
        poly_refine_f = -sum(ccs)
        self%corrs = real(ccs) / real(self%nxpatch*self%nypatch)
    end function poly_refine_f

    subroutine poly_refine_fdf( self, vec, f, grads )
        class(motion_align_poly), intent(inout) :: self
        real(dp),                intent(in)    :: vec(2*POLYDIM)
        real(dp),                intent(out)   :: grads(2*POLYDIM)
        real(dp),                intent(out)   :: f
        real(dp) :: ccs(self%nframes),Es(self%nframes), Ds(self%nframes), Fs(self%nframes)
        real(dp) :: g(6), x,y,rt,sx,sy, sumF, RR, w
        integer  :: i,j,t,ithr,p
        grads = 0.d0
        ccs   = 0.d0
        !$omp parallel do collapse(2) default(shared) private(i,j,t,rt,x,y,ithr,p,sx,sy,w,g,RR,Es,Ds,Fs,sumF)&
        !$omp proc_bind(close) schedule(static) reduction(+:ccs,grads)
        do i = 1,self%nxpatch
            do j = 1,self%nypatch
                ithr = omp_get_thread_num() + 1
                x    = self%patch_coords(i,j,1)
                y    = self%patch_coords(i,j,2)
                call self%ftexp_R(ithr)%zero
                do t = 1,self%nframes
                    rt = real(t-self%align_frame, dp)
                    sx = apply_patch_poly_dp(vec(1:POLYDIM),  x,y, rt)
                    sy = apply_patch_poly_dp(vec(POLYDIM+1:), x,y, rt)
                    call self%ftexp_tiles(t,i,j)%shift( -[sx,sy], self%ftexp_tiles_sh(t,ithr) )
                    call self%ftexp_R(ithr)%add( self%ftexp_tiles_sh(t,ithr) )
                end do
                do t = 1,self%nframes
                    Es(t) = self%ftexp_R( ithr )%corr_unnorm( self%ftexp_tiles_sh(t,ithr) )
                end do
                RR   = sum(Es)
                Ds   = sqrt(RR - 2.d0 * Es + 1.d0)
                Fs   = (Es - 1.d0) / Ds**3.d0
                sumF = sum(Fs)
                ccs  = ccs + (Es-1.d0) / Ds
                call self%ftexp_Rhat(ithr)%zero
                do t = 1,self%nframes
                    call self%ftexp_Rhat(ithr)%add( self%ftexp_tiles_sh(t,ithr), w=1.d0/Ds(t) )
                end do
                ! calc gradient of Rhat, R
                call self%ftexp_Rhat(ithr)%gen_grad_noshift(self%ftexp_dRhat(ithr,1), self%ftexp_dRhat(ithr,2))
                call self%ftexp_R(ithr)%gen_grad_noshift(   self%ftexp_dR(ithr,1),    self%ftexp_dR(ithr,2))
                do p = 1,3
                    call self%ftexp_Rhat2(ithr)%zero
                    do t = 1,self%nframes
                        if( t == self%align_frame ) cycle  ! w = 0.
                        w = real(t-self%align_frame,dp)**p
                        call self%ftexp_Rhat2(ithr)%add( self%ftexp_tiles_sh(t,ithr), w=w)
                    end do
                    g(p)   = self%ftexp_dRhat(ithr,1)%corr_unnorm( self%ftexp_Rhat2(ithr) )
                    g(3+p) = self%ftexp_dRhat(ithr,2)%corr_unnorm( self%ftexp_Rhat2(ithr) )
                    call self%ftexp_Rhat2(ithr)%zero
                    do t = 1,self%nframes
                        if( t == self%align_frame ) cycle  ! w = 0.
                        w = real(t-self%align_frame,dp)**p * (sumF - Fs(t) - 1.d0 / Ds(t))
                        call self%ftexp_Rhat2(ithr)%add( self%ftexp_tiles_sh(t,ithr), w=w)
                    end do
                    g(p)   = g(p)   - self%ftexp_dR(ithr,1)%corr_unnorm( self%ftexp_Rhat2(ithr) )
                    g(3+p) = g(3+p) - self%ftexp_dR(ithr,2)%corr_unnorm( self%ftexp_Rhat2(ithr) )
                end do
                grads( 1:18) = grads( 1:18) - [g(1:3), x*g(1:3), x*x*g(1:3), y*g(1:3), y*y*g(1:3), x*y*g(1:3)]
                grads(19:36) = grads(19:36) - [g(4:6), x*g(4:6), x*x*g(4:6), y*g(4:6), y*y*g(4:6), x*y*g(4:6)]
            end do
        end do
        !$omp end parallel do
        self%corrs = real(ccs) / real(self%nxpatch*self%nypatch)
        f = -sum(ccs)
    end subroutine poly_refine_fdf

end module simple_motion_align_poly
