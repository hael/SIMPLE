! patched-based anisotropic motion correction
module simple_motion_patched
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,                    only: params_glob
use simple_opt_factory,                   only: opt_factory
use simple_opt_spec,                      only: opt_spec
use simple_opt_lbfgsb,                    only: PRINT_NEVALS
use simple_optimizer,                     only: optimizer
use simple_image,                         only: image, imstack_type
use simple_ft_expanded,                   only: ftexp_transfmat_init, ftexp_transfmat_kill
use simple_motion_align_hybrid,           only: motion_align_hybrid
use CPlot2D_wrapper_module
implicit none
private
public :: motion_patched
#include "simple_local_flags.inc"

! module global constants
real,    parameter :: TOL            = 1e-6 !< tolerance parameter
real,    parameter :: TRS_DEFAULT    = 5.
integer, parameter :: PATCH_PDIM     = 18   ! dimension of fitted polynomial

type :: rmat_ptr_type
   real, pointer :: rmat_ptr(:,:,:)
end type rmat_ptr_type

type :: motion_patched
    private
    logical                             :: existence
    type(imstack_type),                  allocatable :: frame_patches(:,:)
    real,                   allocatable :: shifts_patches(:,:,:,:)
    real,                   allocatable :: shifts_patches_for_fit(:,:,:,:)
    real,                   allocatable :: lp(:,:)
    real,                   allocatable :: global_shifts(:,:)       ! isotropic solution
    real,                   allocatable :: frameweights(:)
    real,                   allocatable :: patch_centers(:,:,:)     ! patches centers
    integer,                allocatable :: updateres(:,:)
    integer,                allocatable :: lims_patches(:,:,:,:)    ! patches corners
    character(len=:),       allocatable :: shift_fname
    integer                             :: nframes
    integer                             :: ldim(3)                  ! size of entire frame, reference
    integer                             :: ldim_patch(3)            ! size of one patch
    integer                             :: fixed_frame        = 1   ! frame of reference for alignment
    integer                             :: interp_fixed_frame = 1   ! frame of reference for interpolation
    real                                :: smpd               = 0.
    real                                :: bfactor            = -1.
    real(dp)                            :: poly_coeffs(PATCH_PDIM,2)  ! coefficients of fitted polynomial
    real                                :: polyfit_rmsd(2)            ! polynomial fitting goodness of fit
    real                                :: trs
    real, public                        :: hp
    real                                :: resstep
    logical                             :: has_global_shifts
    logical                             :: has_frameweights  = .false.
    logical                             :: fitshifts         = .false.
contains
    ! Constructors
    procedure                           :: new
    procedure, private                  :: allocate_fields
    procedure, private                  :: deallocate_fields
    ! Drivers & algorithms
    procedure                           :: correct
    procedure                           :: correct_poly
    ! Specific methods
    procedure, private                  :: det_shifts
    ! Generic routine
    procedure, private                  :: set_size_frames_ref
    procedure, private                  :: gen_patch
    procedure, private                  :: fit_polynomial
    procedure, private                  :: get_local_shift
    procedure, private                  :: plot_shifts
    procedure, private                  :: pix2polycoords
    procedure                           :: set_frameweights
    procedure                           :: set_fitshifts
    procedure                           :: set_poly_coeffs
    procedure                           :: set_fixed_frame
    procedure                           :: set_interp_fixed_frame
    procedure                           :: set_nframes
    procedure                           :: set_bfactor
    procedure                           :: get_poly4star
    procedure                           :: get_polyfit_rmsd
    procedure                           :: polytransfo
    ! Destructor
    procedure                           :: kill
end type motion_patched

contains

    ! CONSTRUCTORS

    ! subroutine new( self, motion_correct_ftol, motion_correct_gtol, trs )
    subroutine new( self, trs )
        class(motion_patched), intent(inout) :: self
        real, optional,        intent(in)    :: trs
        call self%kill()
        self%trs = TRS_DEFAULT
        if (present(trs)) self%trs = trs
        allocate(self%lp(params_glob%nxpatch,params_glob%nypatch),&
                &self%updateres(params_glob%nxpatch,params_glob%nypatch),&
                &self%patch_centers(params_glob%nxpatch,params_glob%nypatch,2),&
                &self%lims_patches(params_glob%nxpatch,params_glob%nypatch,2,2))
        self%updateres   = 0
        self%lp          = -1.
        self%hp          = -1.
        self%resstep     = -1.
        self%fixed_frame = 1
        self%bfactor     = -1.
        self%existence = .true.
    end subroutine new

    subroutine allocate_fields( self )
        class(motion_patched), intent(inout) :: self
        logical :: do_allocate
        integer :: i, j
        do_allocate = .true.
        if (allocated(self%shifts_patches)) then
            if (size(self%shifts_patches, dim=1) < self%nframes) then
                do_allocate = .false.
            else
                call self%deallocate_fields()
            end if
        end if
        if (do_allocate) then
            allocate(self%shifts_patches   (2, self%nframes, params_glob%nxpatch, params_glob%nypatch),&
                self%shifts_patches_for_fit(2, self%nframes, params_glob%nxpatch, params_glob%nypatch),&
                self%frame_patches(params_glob%nxpatch, params_glob%nypatch) )
            self%shifts_patches         = 0.
            self%shifts_patches_for_fit = 0.
            do i = 1, params_glob%nxpatch
                do j = 1, params_glob%nypatch
                    allocate( self%frame_patches(i, j)%stack(self%nframes) )
                end do
            end do
        end if
    end subroutine allocate_fields

    subroutine deallocate_fields( self )
        class(motion_patched), intent(inout) :: self
        integer :: iframe, i, j
        if (allocated(self%shifts_patches)) deallocate(self%shifts_patches)
        if (allocated(self%shifts_patches_for_fit)) deallocate(self%shifts_patches_for_fit)
        if (allocated(self%frame_patches)) then
            do j = 1, params_glob%nypatch
                do i = 1, params_glob%nxpatch
                    do iframe = 1, self%nframes
                        call self%frame_patches(i,j)%stack(iframe)%kill()
                    end do
                    deallocate(self%frame_patches(i,j)%stack)
                end do
            end do
            deallocate(self%frame_patches)
        end if
    end subroutine deallocate_fields

    ! DRIVERS

    subroutine correct( self, hp, resstep, frames, shift_fname, global_shifts )
        class(motion_patched),           intent(inout) :: self
        real,                            intent(in)    :: hp, resstep
        type(image),        allocatable, intent(inout) :: frames(:)
        character(len=:),   allocatable, intent(inout) :: shift_fname
        real,     optional, allocatable, intent(in)    :: global_shifts(:,:)
        integer :: ldim_frames(3)
        integer :: iframe
        ! prep
        self%hp          = hp
        self%lp          = params_glob%lpstart
        self%resstep     = resstep
        self%updateres   = 0
        self%shift_fname = shift_fname // C_NULL_CHAR
        if (allocated(self%global_shifts)) deallocate(self%global_shifts)
        self%has_global_shifts = .false.
        if (present(global_shifts)) then
            allocate(self%global_shifts(size(global_shifts, 1), size(global_shifts, 2)))
            self%global_shifts     = global_shifts
            self%has_global_shifts = .true.
        end if
        self%nframes = size(frames,dim=1)
        self%ldim    = frames(1)%get_ldim()
        self%smpd    = frames(1)%get_smpd()
        do iframe = 1,self%nframes
            ldim_frames = frames(iframe)%get_ldim()
            if (any(ldim_frames(1:2) /= self%ldim(1:2))) then
                THROW_HARD('error in motion_patched_correct: frame dimensions do not match reference dimension; simple_motion_patched')
            end if
        end do
        call self%allocate_fields()
        ! determines patch geometry
        call self%set_size_frames_ref()
        ! determine shifts for patches
        call self%det_shifts(frames)
        ! deals with frame of reference convention
        select case(trim(params_glob%mcconvention))
        case('first','relion')
            self%fixed_frame        = 1
            self%interp_fixed_frame = 1
            do iframe = 2, self%nframes
                self%shifts_patches(1,iframe,:,:) = self%shifts_patches(1,iframe,:,:) - self%shifts_patches(1,1,:,:)
                self%shifts_patches(2,iframe,:,:) = self%shifts_patches(2,iframe,:,:) - self%shifts_patches(2,1,:,:)
            enddo
            self%shifts_patches(:,1,:,:) = 0.
            self%shifts_patches_for_fit = self%shifts_patches
        case DEFAULT
            ! all good
        end select
        ! fit the polynomial model against determined shifts
        call self%fit_polynomial()
        ! report visual results
        call self%plot_shifts()
        shift_fname = trim(self%shift_fname) // C_NULL_CHAR
    end subroutine correct

    subroutine correct_poly( self, hp, resstep, rmsd_threshold, frames, shift_fname, patched_polyn, global_shifts )
        use simple_motion_align_poly2, only: motion_align_poly2
        class(motion_patched),           intent(inout) :: self
        real,                            intent(in)    :: hp, resstep, rmsd_threshold
        type(image),        allocatable, intent(inout) :: frames(:)
        character(len=:),   allocatable, intent(inout) :: shift_fname
        real(dp),           allocatable, intent(inout) :: patched_polyn(:)
        real,     optional, allocatable, intent(in)    :: global_shifts(:,:)
        type(motion_align_hybrid), allocatable :: align_hybrid(:,:)
        type(motion_align_poly2) :: align_poly
        real,        allocatable :: opt_shifts(:,:), res(:)
        real(dp)                 :: poly_coeffs(2*PATCH_PDIM)
        real                     :: corr_avg, rmsd(2)
        integer                  :: ldim_frames(3), iframe, i, j, fixed_frame_bak
        self%hp          = hp
        self%lp          = params_glob%lpstart
        self%resstep     = resstep
        self%updateres   = 0
        self%shift_fname = shift_fname // C_NULL_CHAR
        if (allocated(self%global_shifts)) deallocate(self%global_shifts)
        self%has_global_shifts = .false.
        if (present(global_shifts)) then
            allocate(self%global_shifts(size(global_shifts, 1), size(global_shifts, 2)))
            self%global_shifts     = global_shifts
            self%has_global_shifts = .true.
        end if
        self%nframes = size(frames,dim=1)
        self%ldim    = frames(1)%get_ldim()
        self%smpd    = frames(1)%get_smpd()
        do iframe = 1,self%nframes
            ldim_frames = frames(iframe)%get_ldim()
            if (any(ldim_frames(1:2) /= self%ldim(1:2))) then
                THROW_HARD('error in motion_patched_correct: frame dimensions do not match reference dimension; simple_motion_patched')
            end if
        end do
        call self%allocate_fields()
        ! determines patch geometry
        call self%set_size_frames_ref()
        ! determine shifts for patches
        self%shifts_patches = 0.
        allocate(align_hybrid(params_glob%nxpatch, params_glob%nypatch))
        ! initialize transfer matrix to correct dimensions
        call self%frame_patches(1,1)%stack(1)%new(self%ldim_patch, self%smpd, wthreads=.false.)
        call ftexp_transfmat_init(self%frame_patches(1,1)%stack(1), params_glob%lpstop)
        res               = self%frame_patches(1,1)%stack(1)%get_res()
        self%hp           = min(self%hp,res(1))
        corr_avg          = 0.0
        self%frameweights = 1.0 ! need to be one here
        write(logfhandle,'(A,F6.1)')'>>> PATCH HIGH-PASS: ',self%hp
        !$omp parallel do collapse(2) default(shared) private(i,j,iframe,opt_shifts)&
        !$omp proc_bind(close) schedule(dynamic) reduction(+:corr_avg)
        do i = 1,params_glob%nxpatch
            do j = 1,params_glob%nypatch
                ! init
                self%lp(i,j) = params_glob%lpstart
                call self%gen_patch(frames,i,j)
                call align_hybrid(i,j)%new(self%frame_patches(i,j)%stack)
                call align_hybrid(i,j)%set_group_frames(.false.)
                call align_hybrid(i,j)%set_rand_init_shifts(.true.)
                call align_hybrid(i,j)%set_reslims(self%hp, self%lp(i,j), params_glob%lpstop)
                call align_hybrid(i,j)%set_trs(params_glob%scale*params_glob%trs)
                call align_hybrid(i,j)%set_coords(i,j)
                call align_hybrid(i,j)%set_fitshifts(self%fitshifts)
                call align_hybrid(i,j)%set_fixed_frame(self%fixed_frame)
                call align_hybrid(i,j)%set_bfactor(self%bfactor)
                ! align
                call align_hybrid(i,j)%align(frameweights=self%frameweights)
                ! fetch info
                corr_avg = corr_avg + align_hybrid(i,j)%get_corr()
                call align_hybrid(i,j)%get_opt_shifts(opt_shifts)
                ! making sure the shifts are in reference to fixed_frame
                do iframe = 1, self%nframes
                    self%shifts_patches(:,iframe,i,j) = opt_shifts(iframe,:) - opt_shifts(self%fixed_frame,:)
                end do
                ! cleanup
                call align_hybrid(i,j)%kill
                ! do iframe=1,self%nframes
                !     call self%frame_patches(i,j)%stack(iframe)%kill
                ! end do
            end do
        end do
        !$omp end parallel do
        call ftexp_transfmat_kill
        self%shifts_patches_for_fit = self%shifts_patches
        corr_avg = corr_avg / real(params_glob%nxpatch*params_glob%nypatch)
        write(logfhandle,'(A,F6.3)')'>>> AVERAGE PATCH & FRAMES CORRELATION: ', corr_avg
        deallocate(align_hybrid,res)
        call self%fit_polynomial() ! therefore whether polynomial fit was successful is determined with central frame!
        rmsd = self%polyfit_rmsd
        if(allocated(patched_polyn)) deallocate(patched_polyn)
        allocate(patched_polyn(2*PATCH_PDIM),source=0.d0)
        if( all(self%polyfit_rmsd < rmsd_threshold) )then
            ! for polynomial refinement we use fixed_frame=1
            fixed_frame_bak  = self%fixed_frame
            self%fixed_frame = 1 !!
            do iframe = 2, self%nframes
                self%shifts_patches_for_fit(1,iframe,:,:) = self%shifts_patches_for_fit(1,iframe,:,:) - self%shifts_patches_for_fit(1,1,:,:)
                self%shifts_patches_for_fit(2,iframe,:,:) = self%shifts_patches_for_fit(2,iframe,:,:) - self%shifts_patches_for_fit(2,1,:,:)
            enddo
            self%shifts_patches_for_fit(:,1,:,:) = 0.
            call self%fit_polynomial()
            ! refine polynomial
            call align_poly%new(self%frame_patches, self%poly_coeffs, self%patch_centers, self%ldim, self%hp, self%fixed_frame)
            do i = 1,params_glob%nxpatch
            do j = 1,params_glob%nypatch
            do iframe=1,self%nframes
                call self%frame_patches(i,j)%stack(iframe)%kill
            end do
            end do
            end do
            call align_poly%refine(self%poly_coeffs, self%frameweights)
            patched_polyn(:PATCH_PDIM)   = self%poly_coeffs(:,1)
            patched_polyn(PATCH_PDIM+1:) = self%poly_coeffs(:,2)
            write(logfhandle,'(A,F6.3)')'>>> AVERAGE POLYNOMIAL CORRELATION: ', align_poly%get_corr()
            call align_poly%kill
            select case(trim(params_glob%mcconvention))
            case('first','relion')
                self%interp_fixed_frame = 1
            case DEFAULT
                ! need to update frame of reference & polynomial
                do i = 1,params_glob%nxpatch
                do j = 1,params_glob%nypatch
                do iframe=1,self%nframes
                    call self%get_local_shift(iframe, self%patch_centers(i,j,1),self%patch_centers(i,j,2), self%shifts_patches_for_fit(:,iframe,i,j))
                enddo
                enddo
                enddo
                self%fixed_frame = fixed_frame_bak
                do iframe = 1, self%nframes
                    if( iframe == self%fixed_frame ) cycle
                    self%shifts_patches_for_fit(1,iframe,:,:) = self%shifts_patches_for_fit(1,iframe,:,:) - self%shifts_patches_for_fit(1,self%fixed_frame,:,:)
                    self%shifts_patches_for_fit(2,iframe,:,:) = self%shifts_patches_for_fit(2,iframe,:,:) - self%shifts_patches_for_fit(2,self%fixed_frame,:,:)
                enddo
                self%shifts_patches_for_fit(:,self%fixed_frame,:,:) = 0.
                call self%fit_polynomial()
                self%shifts_patches_for_fit = self%shifts_patches
            end select
            self%polyfit_rmsd = rmsd
        endif
        do i = 1,params_glob%nxpatch
        do j = 1,params_glob%nypatch
        do iframe=1,self%nframes
            call self%frame_patches(i,j)%stack(iframe)%kill
        end do
        end do
        end do
        ! report visual results
        call self%plot_shifts()
        shift_fname = trim(self%shift_fname) // C_NULL_CHAR
    end subroutine correct_poly

    ! SPECIFIC METHODS

    subroutine det_shifts( self, frames )
        class(motion_patched), target, intent(inout) :: self
        type(image),      allocatable, intent(inout) :: frames(:)
        type(motion_align_hybrid), allocatable :: align_hybrid(:,:)
        real, allocatable :: opt_shifts(:,:), res(:)
        real              :: corr_avg
        integer           :: iframe, i, j
        logical           :: l_groupframes
        l_groupframes = trim(params_glob%groupframes).eq.'yes'
        self%shifts_patches = 0.
        allocate(align_hybrid(params_glob%nxpatch, params_glob%nypatch))
        corr_avg = 0.
        ! initialize transfer matrix to correct dimensions
        call self%frame_patches(1,1)%stack(1)%new(self%ldim_patch, self%smpd, wthreads=.false.)
        call ftexp_transfmat_init(self%frame_patches(1,1)%stack(1), params_glob%lpstop)
        res      = self%frame_patches(1,1)%stack(1)%get_res()
        self%hp  = min(self%hp,res(1))
        corr_avg = 0.
        write(logfhandle,'(A,F6.1)')'>>> PATCH HIGH-PASS: ',self%hp
        !$omp parallel do collapse(2) default(shared) private(i,j,iframe,opt_shifts)&
        !$omp proc_bind(close) schedule(dynamic) reduction(+:corr_avg)
        do i = 1,params_glob%nxpatch
            do j = 1,params_glob%nypatch
                ! init
                self%lp(i,j) = params_glob%lpstart
                call self%gen_patch(frames,i,j)
                call align_hybrid(i,j)%new(self%frame_patches(i,j)%stack)
                call align_hybrid(i,j)%set_group_frames(l_groupframes)
                call align_hybrid(i,j)%set_rand_init_shifts(.true.)
                call align_hybrid(i,j)%set_reslims(self%hp, self%lp(i,j), params_glob%lpstop)
                call align_hybrid(i,j)%set_trs(params_glob%scale*params_glob%trs)
                call align_hybrid(i,j)%set_coords(i,j)
                call align_hybrid(i,j)%set_fitshifts(self%fitshifts)
                call align_hybrid(i,j)%set_fixed_frame(self%fixed_frame)
                call align_hybrid(i,j)%set_bfactor(self%bfactor)
                ! align
                call align_hybrid(i,j)%align(frameweights=self%frameweights)
                ! fetch info
                corr_avg = corr_avg + align_hybrid(i,j)%get_corr()
                call align_hybrid(i,j)%get_opt_shifts(opt_shifts)
                ! making sure the shifts are in reference to fixed_frame
                do iframe = 1, self%nframes
                    self%shifts_patches(:,iframe,i,j) = opt_shifts(iframe,:) - opt_shifts(self%fixed_frame,:)
                end do
                ! cleanup
                call align_hybrid(i,j)%kill
                do iframe=1,self%nframes
                    call self%frame_patches(i,j)%stack(iframe)%kill
                end do
            end do
        end do
        !$omp end parallel do
        self%shifts_patches_for_fit = self%shifts_patches
        corr_avg = corr_avg / real(params_glob%nxpatch*params_glob%nypatch)
        write(logfhandle,'(A,F6.3)')'>>> AVERAGE PATCH & FRAMES CORRELATION: ', corr_avg
        deallocate(align_hybrid,res)
        call ftexp_transfmat_kill
    end subroutine det_shifts

    ! OTHER

    subroutine fit_polynomial( self )
        class(motion_patched), intent(inout) :: self
        real(dp) :: yx(self%nframes*params_glob%nxpatch*params_glob%nypatch)      ! along x
        real(dp) :: yy(self%nframes*params_glob%nxpatch*params_glob%nypatch)      ! along y
        real(dp) :: x(3,self%nframes*params_glob%nxpatch*params_glob%nypatch)     ! x,y,t
        real(dp) :: sig(self%nframes*params_glob%nxpatch*params_glob%nypatch)
        real(dp) :: v(PATCH_PDIM,PATCH_PDIM), w(PATCH_PDIM), chisq
        real     :: fitted_shift(2)
        integer  :: idx, iframe, i, j
        ! fitting
        sig = 1.d0
        idx = 0
        do iframe = 1, self%nframes
            do i = 1, params_glob%nxpatch
                do j = 1, params_glob%nypatch
                    idx     = idx+1
                    yx(idx) = real(self%shifts_patches_for_fit(1,iframe,i,j),dp)
                    yy(idx) = real(self%shifts_patches_for_fit(2,iframe,i,j),dp)
                    call self%pix2polycoords(real(self%patch_centers(i,j,1),dp),real(self%patch_centers(i,j,2),dp),x(1,idx),x(2,idx))
                    x(3,idx) = real(iframe-self%fixed_frame,dp)
                end do
            end do
        end do
        call svd_multifit(x,yx,sig,self%poly_coeffs(:,1),v,w,chisq,patch_poly)
        call svd_multifit(x,yy,sig,self%poly_coeffs(:,2),v,w,chisq,patch_poly)
        ! goodness of fit
        idx = 0
        self%polyfit_rmsd = 0.
        do iframe = 1,self%nframes
            do i = 1,params_glob%nxpatch
                do j = 1,params_glob%nypatch
                    idx = idx+1
                    call self%get_local_shift(iframe, self%patch_centers(i,j,1),self%patch_centers(i,j,2),fitted_shift)
                    self%polyfit_rmsd(1) = self%polyfit_rmsd(1) + (fitted_shift(1)-yx(idx))**2.
                    self%polyfit_rmsd(2) = self%polyfit_rmsd(2) + (fitted_shift(2)-yy(idx))**2.
                end do
            end do
        end do
        self%polyfit_rmsd = sqrt(self%polyfit_rmsd/real(self%nframes*params_glob%nxpatch*params_glob%nypatch))
    end subroutine fit_polynomial

    subroutine plot_shifts( self )
        class(motion_patched), intent(inout) :: self
        real, parameter           :: SCALE = 40.
        type(str4arr)             :: title
        type(CPlot2D_type)        :: plot2D
        type(CDataSet_type)       :: dataSetStart, dataSet, fit, obs, patch_start
        type(CDataPoint_type)     :: point2, p_obs, p_fit, point
        character(len=LONGSTRLEN) :: ps2pdf_cmd, fname_pdf, ps2jpeg_cmd, fname_jpeg
        integer :: l, ipx,ipy, iframe, j, iostat
        real    :: shifts(self%nframes,2), loc_shift(2), ref_shift(2), xcenter,ycenter, cx,cy
        call CPlot2D__new(plot2D, self%shift_fname)
        call CPlot2D__SetXAxisSize(plot2D, 600._c_double)
        call CPlot2D__SetYAxisSize(plot2D, 600._c_double)
        call CPlot2D__SetDrawLegend(plot2D, C_FALSE)
        call CPlot2D__SetFlipY(plot2D, C_TRUE)
        if (self%has_global_shifts) then
            ! centering to first frame for display only
            shifts      = self%global_shifts
            shifts(:,1) = shifts(:,1) - shifts(1,1)
            shifts(:,2) = shifts(:,2) - shifts(1,2)
            ! plot
            call CDataSet__new(dataSet)
            call CDataSet__SetDrawMarker(dataSet, C_FALSE)
            call CDataSet__SetDatasetColor(dataSet, 0.0_c_double, 0.0_c_double, 1.0_c_double)
            xcenter = real(self%ldim(1))/2.
            ycenter = real(self%ldim(2))/2.
            do j = 1, self%nframes
                call CDataPoint__new2( real(xcenter + SCALE * shifts(j, 1), c_double), &
                                      &real(ycenter + SCALE * shifts(j, 2), c_double), point)
                call CDataSet__AddDataPoint(dataSet, point)
                call CDataPoint__delete(point)
            end do
            call CPlot2D__AddDataSet(plot2D, dataset)
            call CDataSet__delete(dataset)
        end if
        call CDataSet__new(dataSetStart)
        call CDataSet__SetDrawMarker(dataSetStart, C_TRUE)
        call CDataSet__SetMarkerSize(dataSetStart,5._c_double)
        call CDataSet__SetDatasetColor(dataSetStart, 1.0_c_double,0.0_c_double,0.0_c_double)
        call CDataPoint__new2(real(xcenter,c_double), real(ycenter,c_double), point2)
        call CDataSet__AddDataPoint(dataSetStart, point2)
        call CPlot2D__AddDataSet(plot2D, dataSetStart)
        do ipx = 1, params_glob%nxpatch
            do ipy = 1, params_glob%nypatch
                cx = self%patch_centers(ipx, ipy, 1)
                cy = self%patch_centers(ipx, ipy, 2)
                ! centering to first frame for display only
                shifts(:,1) = self%shifts_patches_for_fit(1,:,ipx,ipy)
                shifts(:,2) = self%shifts_patches_for_fit(2,:,ipx,ipy)
                shifts(:,1) = shifts(:,1) - shifts(1,1)
                shifts(:,2) = shifts(:,2) - shifts(1,2)
                call self%get_local_shift(1,cx,cy, ref_shift)
                ref_shift = ref_shift - shifts(1,:)
                ! plot
                call CDataSet__new(patch_start)
                call CDataSet__SetDrawMarker(patch_start,C_TRUE)
                call CDataSet__SetMarkerSize(patch_start,5.0_c_double)
                call CDataSet__SetDatasetColor(patch_start,1.0_c_double,0.0_c_double,0.0_c_double)
                call CDataSet__new(fit)
                call CDataSet__new(obs)
                call CDataSet__SetDrawMarker(fit, C_FALSE)
                call CDataSet__SetDatasetColor(fit, 0.0_c_double,0.0_c_double,0.0_c_double)
                call CDataSet__SetDrawMarker(obs, C_FALSE)
                call CDataSet__SetDatasetColor(obs, 0.5_c_double,0.5_c_double,0.5_c_double)
                call CDataPoint__new2(real(cx,c_double), real(cy, c_double), p_fit)
                call CDataSet__AddDataPoint(patch_start, p_fit)
                call CDataPoint__delete(p_fit)
                do iframe = 1, self%nframes
                    call CDataPoint__new2(real(cx + SCALE*shifts(iframe,1), c_double),&
                                         &real(cy + SCALE*shifts(iframe,2), c_double), p_obs)
                    call CDataSet__AddDataPoint(obs, p_obs)
                    call self%get_local_shift(iframe,cx,cy, loc_shift)
                    loc_shift = loc_shift - ref_shift
                    call CDataPoint__new2(real(cx + SCALE*loc_shift(1), c_double),&
                                         &real(cy + SCALE*loc_shift(2), c_double), p_fit)
                    call CDataSet__AddDataPoint(fit, p_fit)
                    call CDataPoint__delete(p_fit)
                    call CDataPoint__delete(p_obs)
                end do
                call CPlot2D__AddDataSet(plot2D, obs)
                call CPlot2D__AddDataSet(plot2D, fit)
                call CPlot2D__AddDataSet(plot2D, patch_start)
                call CDataSet__delete(patch_start)
                call CDataSet__delete(fit)
                call CDataSet__delete(obs)
            end do
        end do
        title%str = 'X (in pixels; trajectory scaled by '//trim(int2str(nint(SCALE)))//')'//C_NULL_CHAR
        call CPlot2D__SetXAxisTitle(plot2D, title%str)
        title%str(1:1) = 'Y'
        call CPlot2D__SetYAxisTitle(plot2D, title%str)
        call CPlot2D__OutputPostScriptPlot(plot2D, self%shift_fname)
        call CPlot2D__delete(plot2D)
        ! Format conversion
        l = len_trim(self%shift_fname)
        self%shift_fname = self%shift_fname(:l-1) ! removing trailing C NULL character
        ! conversion to JPEG
        fname_jpeg  = trim(get_fbody(self%shift_fname,'eps'))//'.jpeg'
        ps2jpeg_cmd = 'gs -q -sDEVICE=jpeg -dJPEGQ=92 -dNOPAUSE -dBATCH -dSAFER -dDEVICEWIDTHPOINTS=760 -dDEVICEHEIGHTPOINTS=760 -sOutputFile='&
            //trim(fname_jpeg)//' '//trim(self%shift_fname)
        call exec_cmdline(trim(adjustl(ps2jpeg_cmd)), suppress_errors=.true., exitstat=iostat)
        ! conversion to PDF
        fname_pdf  = trim(get_fbody(self%shift_fname,'eps'))//'.pdf'
        ps2pdf_cmd = 'gs -q -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dDEVICEWIDTHPOINTS=760 -dDEVICEHEIGHTPOINTS=760 -sOutputFile='&
            //trim(fname_pdf)//' '//trim(self%shift_fname)
        call exec_cmdline(trim(adjustl(ps2pdf_cmd)), suppress_errors=.true., exitstat=iostat)
        ! update name
        if( iostat == 0 )then
            call del_file(self%shift_fname)
            self%shift_fname = trim(fname_pdf)
        endif
    end subroutine plot_shifts

    ! produces shifts for 'polishing' close to relion 3 convention
    subroutine get_poly4star( self, polycoeffs, patch_shifts, patch_centers )
        class(motion_patched), intent(inout) :: self
        real(dp), allocatable, intent(inout) :: polycoeffs(:), patch_shifts(:,:,:,:), patch_centers(:,:,:)
        real(dp) :: yx(self%nframes*params_glob%nxpatch*params_glob%nypatch)      ! along x
        real(dp) :: yy(self%nframes*params_glob%nxpatch*params_glob%nypatch)      ! along y
        real(dp) :: x(3,self%nframes*params_glob%nxpatch*params_glob%nypatch)     ! x,y,t
        real(dp) :: sig(self%nframes*params_glob%nxpatch*params_glob%nypatch)
        real(dp) :: v(PATCH_PDIM,PATCH_PDIM), w(PATCH_PDIM), chisq
        integer  :: idx, iframe, i, j
        if( allocated(polycoeffs) )    deallocate(polycoeffs)
        if( allocated(patch_shifts) )  deallocate(patch_shifts)
        if( allocated(patch_centers) ) deallocate(patch_centers)
        allocate(polycoeffs(2*PATCH_PDIM),patch_shifts(2,self%nframes,params_glob%nxpatch,params_glob%nypatch),&
            &patch_centers(params_glob%nxpatch,params_glob%nypatch,2),source=0.d0)
        idx = 0
        do iframe = 1, self%nframes
            do i = 1, params_glob%nxpatch
                do j = 1, params_glob%nypatch
                    idx     = idx+1
                    yx(idx) = real(self%shifts_patches_for_fit(1,iframe,i,j)-self%shifts_patches_for_fit(1,1,i,j),dp)
                    yy(idx) = real(self%shifts_patches_for_fit(2,iframe,i,j)-self%shifts_patches_for_fit(2,1,i,j),dp)
                    call self%pix2polycoords(real(self%patch_centers(i,j,1),dp),&
                                            &real(self%patch_centers(i,j,2),dp),x(1,idx),x(2,idx))
                    x(3,idx) = real(iframe-1,dp)
                    patch_shifts(1,iframe,i,j) = yx(idx)
                    patch_shifts(2,iframe,i,j) = yy(idx)
                    patch_centers(i,j,:)       = real(self%patch_centers(i,j,:),dp)
                end do
            end do
        end do
        sig = 1.d0
        call svd_multifit(x,yx,sig,polycoeffs(           1:  PATCH_PDIM),v,w,chisq,patch_poly)
        call svd_multifit(x,yy,sig,polycoeffs(PATCH_PDIM+1:2*PATCH_PDIM),v,w,chisq,patch_poly)
    end subroutine get_poly4star

    function get_polyfit_rmsd( self )result( rmsd )
        class(motion_patched), intent(in) :: self
        real :: rmsd(2)
        rmsd = self%polyfit_rmsd
    end function get_polyfit_rmsd

    elemental subroutine pix2polycoords( self, xin, yin, x, y )
        class(motion_patched), intent(in)  :: self
        real(dp),              intent(in)  :: xin, yin
        real(dp),              intent(out) :: x, y
        x = (xin-1.d0) / real(self%ldim(1)-1,dp) - 0.5d0
        y = (yin-1.d0) / real(self%ldim(2)-1,dp) - 0.5d0
    end subroutine pix2polycoords

    pure subroutine get_local_shift( self, iframe, x, y, shift )
        class(motion_patched), intent(inout) :: self
        integer,               intent(in)  :: iframe
        real,                  intent(in)  :: x, y
        real,                  intent(out) :: shift(2)
        real(dp) :: t, xx, yy
        t  = real(iframe-self%fixed_frame, dp)
        call self%pix2polycoords(real(x,dp),real(y,dp), xx,yy)
        shift(1) = apply_patch_poly(self%poly_coeffs(:,1), xx,yy,t)
        shift(2) = apply_patch_poly(self%poly_coeffs(:,2), xx,yy,t)
    end subroutine get_local_shift

    !>  Per frame real space polynomial interpolation
    subroutine polytransfo( self, frames, weights, frame_output )
        class(motion_patched), intent(inout) :: self
        type(image),           intent(inout) :: frames(self%nframes)
        real,                  intent(in)    :: weights(self%nframes)
        type(image),           intent(inout) :: frame_output
        real, pointer :: rmatin(:,:,:), rmatout(:,:,:)
        real(dp)      :: t,ti, dt,dt2,dt3, x,x2,y,y2,xy, A1,A2, B1x,B1x2,B1xy,B2x,B2x2,B2xy
        integer       :: ldim(3), i, j, iframe
        real          :: w, pixx,pixy
        ldim    = frames(1)%get_ldim()
        call frame_output%zero_and_unflag_ft
        call frame_output%get_rmat_ptr(rmatout)
        ti = real(self%interp_fixed_frame-self%fixed_frame, dp)
        do iframe = 1,self%nframes
            call frames(iframe)%get_rmat_ptr(rmatin)
            w = weights(iframe)
            t = real(iframe-self%fixed_frame, dp)
            dt  = ti-t
            dt2 = ti*ti - t*t
            dt3 = ti*ti*ti - t*t*t
            B1x  = sum(self%poly_coeffs(4:6,1)   * [dt,dt2,dt3])
            B1x2 = sum(self%poly_coeffs(7:9,1)   * [dt,dt2,dt3])
            B1xy = sum(self%poly_coeffs(16:18,1) * [dt,dt2,dt3])
            B2x  = sum(self%poly_coeffs(4:6,2)   * [dt,dt2,dt3])
            B2x2 = sum(self%poly_coeffs(7:9,2)   * [dt,dt2,dt3])
            B2xy = sum(self%poly_coeffs(16:18,2) * [dt,dt2,dt3])
            !$omp parallel do default(shared) private(i,j,x,x2,y,y2,xy,A1,A2,pixx,pixy)&
            !$omp proc_bind(close) schedule(static)
            do j = 1, ldim(2)
                y  = real(j-1,dp) / real(ldim(2)-1,dp) - 0.5d0
                y2 = y*y
                A1 =           sum(self%poly_coeffs(1:3,1)   * [dt,dt2,dt3])
                A1 = A1 + y  * sum(self%poly_coeffs(10:12,1) * [dt,dt2,dt3])
                A1 = A1 + y2 * sum(self%poly_coeffs(13:15,1) * [dt,dt2,dt3])
                A2 =           sum(self%poly_coeffs(1:3,2)   * [dt,dt2,dt3])
                A2 = A2 + y  * sum(self%poly_coeffs(10:12,2) * [dt,dt2,dt3])
                A2 = A2 + y2 * sum(self%poly_coeffs(13:15,2) * [dt,dt2,dt3])
                do i = 1, ldim(1)
                    x  = real(i-1,dp) / real(ldim(1)-1,dp) - 0.5d0
                    x2 = x*x
                    xy = x*y
                    pixx = real(i) + real(A1 + B1x*x + B1x2*x2 + B1xy*xy)
                    pixy = real(j) + real(A2 + B2x*x + B2x2*x2 + B2xy*xy)
                    rmatout(i,j,1) = rmatout(i,j,1) + w*interp_bilin(pixx,pixy)
                end do
            end do
            !$omp end parallel do
        enddo
    contains

        pure real function interp_bilin( xval, yval )
            real, intent(in) :: xval, yval
            integer  :: x1_h,  x2_h,  y1_h,  y2_h
            real     :: t, u
            logical  :: outside
            outside = .false.
            x1_h = floor(xval)
            x2_h = x1_h + 1
            if( x1_h<1 .or. x2_h<1 )then
                x1_h    = 1
                outside = .true.
            endif
            if( x1_h>ldim(1) .or. x2_h>ldim(1) )then
                x1_h    = ldim(1)
                outside = .true.
            endif
            y1_h = floor(yval)
            y2_h = y1_h + 1
            if( y1_h<1 .or. y2_h<1 )then
                y1_h    = 1
                outside = .true.
            endif
            if( y1_h>ldim(2) .or. y2_h>ldim(2) )then
                y1_h    = ldim(2)
                outside = .true.
            endif
            if( outside )then
                interp_bilin = rmatin(x1_h, y1_h, 1)
                return
            endif
            t  = xval - real(x1_h)
            u  = yval - real(y1_h)
            interp_bilin =  (1. - t) * (1. - u) * rmatin(x1_h, y1_h, 1) + &
                                 &t  * (1. - u) * rmatin(x2_h, y1_h, 1) + &
                                 &t  *       u  * rmatin(x2_h, y2_h, 1) + &
                           &(1. - t) *       u  * rmatin(x1_h, y2_h, 1)
        end function interp_bilin

    end subroutine polytransfo

    subroutine set_size_frames_ref( self )
        class(motion_patched), intent(inout) :: self
        integer :: i,j
        real    :: cen, dist
        self%ldim_patch(1) = round2even(real(self%ldim(1)) / real(params_glob%nxpatch))
        self%ldim_patch(2) = round2even(real(self%ldim(2)) / real(params_glob%nypatch))
        self%ldim_patch(3) = 1
        ! fftw friendly size
        self%ldim_patch(1) = find_larger_magic_box(self%ldim_patch(1))
        self%ldim_patch(2) = find_larger_magic_box(self%ldim_patch(2))
        ! along X
        ! limits & center first patches
        self%lims_patches(1,:,1,1) = 1
        self%lims_patches(1,:,1,2) = self%ldim_patch(1)
        self%patch_centers(1,:,1)  = sum(self%lims_patches(1,:,1,1:2),dim=2) / 2.
        ! limits & center last patches
        self%lims_patches(params_glob%nxpatch,:,1,1) = self%ldim(1)-self%ldim_patch(1)+1
        self%lims_patches(params_glob%nxpatch,:,1,2) = self%ldim(1)
        self%patch_centers(params_glob%nxpatch,:,1)  = sum(self%lims_patches(params_glob%nxpatch,:,1,1:2),dim=2) / 2.
        ! adjust other patch centers to be evenly spread
        dist = real(self%patch_centers(params_glob%nxpatch,1,1)-self%patch_centers(1,1,1)+1) / real(params_glob%nxpatch-1)
        do i=2,params_glob%nxpatch-1
            cen = self%patch_centers(1,1,1) + real(i-1)*dist
            self%lims_patches(i,:,1,1) = ceiling(cen) - self%ldim_patch(1)/2
            self%lims_patches(i,:,1,2) = self%lims_patches(i,:,1,1) + self%ldim_patch(1) - 1
            self%patch_centers(i,:,1)  = sum(self%lims_patches(i,:,1,1:2),dim=2) / 2.
        enddo
        ! along Y
        self%lims_patches(:,1,2,1) = 1
        self%lims_patches(:,1,2,2) = self%ldim_patch(2)
        self%patch_centers(:,1,2)  = sum(self%lims_patches(:,1,2,1:2),dim=2) / 2.
        self%lims_patches(:,params_glob%nypatch,2,1) = self%ldim(2)-self%ldim_patch(2)+1
        self%lims_patches(:,params_glob%nypatch,2,2) = self%ldim(2)
        self%patch_centers(:,params_glob%nypatch,2)  = sum(self%lims_patches(:,params_glob%nypatch,2,1:2),dim=2) / 2.
        dist = real(self%patch_centers(1,params_glob%nypatch,2)-self%patch_centers(1,1,2)+1) / real(params_glob%nypatch-1)
        do j=2,params_glob%nypatch-1
            cen = self%patch_centers(1,1,2) + real(j-1)*dist
            self%lims_patches(:,j,2,1) = ceiling(cen) - self%ldim_patch(2)/2
            self%lims_patches(:,j,2,2) = self%lims_patches(:,j,2,1) + self%ldim_patch(2) - 1
            self%patch_centers(:,j,2)  = sum(self%lims_patches(:,j,2,1:2),dim=2) /2.
        enddo
    end subroutine set_size_frames_ref

    subroutine gen_patch( self, stack, pi, pj )
        class(motion_patched),          intent(inout) :: self
        type(image),       allocatable, intent(inout) :: stack(:)
        integer,                        intent(in)    :: pi, pj
        integer             :: iframe, k, l, ip, jp           ! ip, jp: i_patch, j_patch
        real, pointer       :: prmat_patch(:,:,:), prmat_frame(:,:,:)
        do iframe=1,self%nframes
            ! init
            call self%frame_patches(pi,pj)%stack(iframe)%new(self%ldim_patch, self%smpd, wthreads=.false.)
            call self%frame_patches(pi,pj)%stack(iframe)%get_rmat_ptr(prmat_patch)
            call stack(iframe)%get_rmat_ptr(prmat_frame)
            ! copy
            do l = self%lims_patches(pi,pj,2,1), self%lims_patches(pi,pj,2,2)
                jp = l - self%lims_patches(pi,pj,2,1) + 1
                do k = self%lims_patches(pi,pj,1,1), self%lims_patches(pi,pj,1,2)
                    ip = k - self%lims_patches(pi,pj,1,1) + 1
                    prmat_patch(ip,jp,1) = prmat_frame(k,l,1)
                end do
            end do
        end do
    end subroutine gen_patch

    subroutine set_frameweights( self, frameweights )
        class(motion_patched), intent(inout) :: self
        real, allocatable,     intent(in) :: frameweights(:)
        integer :: nlen
        nlen = size(frameweights)
        if (allocated(self%frameweights)) deallocate(self%frameweights)
        allocate(self%frameweights(nlen), source=frameweights)
        self%has_frameweights = .true.
    end subroutine set_frameweights

    subroutine set_fitshifts( self, fitshifts )
        class(motion_patched), intent(inout) :: self
        logical,               intent(in)    :: fitshifts
        self%fitshifts = fitshifts
    end subroutine set_fitshifts

    subroutine set_poly_coeffs( self, p )
        class(motion_patched), intent(inout) :: self
        real(dp),              intent(in)    :: p(PATCH_PDIM,2)
        self%poly_coeffs = p
    end subroutine set_poly_coeffs

    subroutine set_fixed_frame( self, fixed_frame )
        class(motion_patched), intent(inout) :: self
        integer,               intent(in)    :: fixed_frame
        self%fixed_frame = fixed_frame
    end subroutine set_fixed_frame

    subroutine set_interp_fixed_frame( self, fixed_frame )
        class(motion_patched), intent(inout) :: self
        integer,               intent(in)    :: fixed_frame
        self%interp_fixed_frame = fixed_frame
    end subroutine set_interp_fixed_frame

    subroutine set_nframes( self, n )
        class(motion_patched), intent(inout) :: self
        integer,               intent(in)    :: n
        self%nframes = n
    end subroutine set_nframes

    subroutine set_bfactor( self, bfac )
        class(motion_patched), intent(inout) :: self
        real,                  intent(in)    :: bfac
        self%bfactor = bfac
    end subroutine set_bfactor

    ! DESTRUCTOR

    subroutine kill( self )
        class(motion_patched), intent(inout) :: self
        call self%deallocate_fields()
        if (allocated(self%frameweights)) deallocate(self%frameweights)
        if (allocated(self%updateres)) deallocate(self%updateres)
        if (allocated(self%lp)) deallocate(self%lp)
        if (allocated(self%patch_centers)) deallocate(self%patch_centers)
        if (allocated(self%lims_patches)) deallocate(self%lims_patches)
        self%has_frameweights = .false.
        self%existence        = .false.
        call ftexp_transfmat_kill
    end subroutine kill

    ! POLYNOMIAL DEFORMATION MODEL UTILITIES

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
        res(    2) = t**2
        res(    3) = t**3
        res( 4: 6) = x * res( 1: 3)  ! x   * {t,t^2,t^3}
        res( 7: 9) = x * res( 4: 6)  ! x^2 * {t,t^2,t^3}
        res(10:12) = y * res( 1: 3)  ! y   * {t,t^2,t^3}
        res(13:15) = y * res(10:12)  ! y^2 * {t,t^2,t^3}
        res(16:18) = y * res( 4: 6)  ! x*y * {t,t^2,t^3}
    end function patch_poly

    pure function apply_patch_poly(c, x, y, t) result(res_sp)
        real(dp), intent(in) :: c(PATCH_PDIM), x, y, t
        real(sp) :: res_sp
        real(dp) :: res
        real(dp) :: x2, y2, xy, t2, t3
        x2 = x * x
        y2 = y * y
        xy = x * y
        t2 = t * t
        t3 = t2 * t
        res =       c( 1) * t      + c( 2) * t2      + c( 3) * t3
        res = res + c( 4) * t * x  + c( 5) * t2 * x  + c( 6) * t3 * x
        res = res + c( 7) * t * x2 + c( 8) * t2 * x2 + c( 9) * t3 * x2
        res = res + c(10) * t * y  + c(11) * t2 * y  + c(12) * t3 * y
        res = res + c(13) * t * y2 + c(14) * t2 * y2 + c(15) * t3 * y2
        res = res + c(16) * t * xy + c(17) * t2 * xy + c(18) * t3 * xy
        res_sp = real(res)
    end function apply_patch_poly

end module simple_motion_patched
