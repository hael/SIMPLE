! patched-based anisotropic motion correction
module simple_motion_patched
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,      only: params_glob
use simple_opt_factory,     only: opt_factory
use simple_opt_spec,        only: opt_spec
use simple_optimizer,       only: optimizer
use simple_image,           only: image
use simple_ft_expanded,     only: ft_expanded
use simple_ftexp_shsrch,    only: ftexp_shsrch
implicit none
private
public :: motion_patched
#include "simple_local_flags.inc"

! module global constants
integer, parameter :: NX_PATCHED     = 5   ! number of patches in x-direction
integer, parameter :: NY_PATCHED     = 5   !       "      "       y-direction
integer, parameter :: X_OVERLAP      = 0   ! number of overlapping pixels per patch in x-direction
integer, parameter :: Y_OVERLAP      = 0   !       "      "        "         "         y-direction
real,    parameter :: TOL            = 1e-6 !< tolerance parameter
real,    parameter :: TRS_DEFAULT    = 7.
integer, parameter :: PATCH_PDIM     = 18  ! dimension of fitted polynomial

type :: rmat_ptr_type
    real, pointer :: rmat_ptr(:,:,:)
end type rmat_ptr_type

type :: motion_patched
    private
    logical                          :: existence
    type(ft_expanded),  allocatable  :: frame_patches(:,:,:)
    type(ft_expanded),  allocatable  :: ref_patches(:,:,:)
    type(image),        allocatable  :: patches_imgs(:,:,:)
    type(ftexp_shsrch), allocatable  :: shsearch_patches(:,:,:)
    real,               allocatable  :: shifts_patches(:,:,:,:)
    integer                          :: nframes
    integer                          :: ldim(3)       ! size of entire frame, reference
    integer                          :: ldim_patch(3) ! size of one patch
    integer                          :: lims_patches(NX_PATCHED,NY_PATCHED,2,2) ! corners of the patches
    real                             :: motion_correct_ftol
    real                             :: motion_correct_gtol
    real                             :: poly_coeffs(PATCH_PDIM,2)  ! coefficients of fitted polynomial
    real                             :: trs
    real                             :: hp
    real                             :: lp
contains
    procedure, private               :: allocate_fields
    procedure, private               :: deallocate_fields
    procedure, private               :: set_size_frames_ref
    procedure, private               :: set_patches
    procedure, private               :: det_shifts
    procedure, private               :: fit_polynomial
    procedure, private               :: apply_polytransfo
    procedure, private               :: write_shifts
    procedure                        :: new             => motion_patched_new
    procedure                        :: correct         => motion_patched_correct
    procedure                        :: kill            => motion_patched_kill
end type motion_patched

contains

    ! Polynomial for patch motion
    function patch_poly(p, n) result(res)
        real, dimension(:), intent(in) :: p
        integer,            intent(in) :: n
        real :: res(n)
        real :: x, y, t
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

    function apply_patch_poly(c, x, y, t) result(res)
        real, intent(in) :: c(PATCH_PDIM), x, y, t
        real :: res
        real :: x2, y2, xy, t2, t3
        x2 = x * x
        y2 = y * y
        xy = x * y
        t2 = t * t
        t3 = t2 * t
        res = 0.
        res = res + c( 1) * t      + c( 2) * t2      + c( 3) * t3
        res = res + c( 4) * t * x  + c( 5) * t2 * x  + c( 6) * t3 * x
        res = res + c( 7) * t * x2 + c( 8) * t2 * x2 + c( 9) * t3 * x2
        res = res + c(10) * t * y  + c(11) * t2 * y  + c(12) * t3 * y
        res = res + c(13) * t * y2 + c(14) * t2 * y2 + c(15) * t3 * y2
        res = res + c(16) * t * xy + c(17) * t2 * xy + c(18) * t3 * xy
    end function apply_patch_poly

    subroutine fit_polynomial( self )
        class(motion_patched), intent(inout) :: self
        real :: y(  self%nframes*NX_PATCHED*NY_PATCHED)
        real :: x(3,self%nframes*NX_PATCHED*NY_PATCHED)    ! x,y,t
        real :: sig(self%nframes*NX_PATCHED*NY_PATCHED)
        real :: a(PATCH_PDIM), v(PATCH_PDIM,PATCH_PDIM), w(PATCH_PDIM), chisq
        integer :: alloc_stat
        integer :: iframe, i, j
        integer :: idx
        real :: patch_cntr_x, patch_cntr_y   ! patch centers, x-/y-coordinate
        do iframe = 1, self%nframes
            do i = 1, NX_PATCHED
                do j = 1, NY_PATCHED
                    patch_cntr_x = real(self%lims_patches(i,j,1,1) + self%lims_patches(i,j,1,2)) / 2.
                    patch_cntr_y = real(self%lims_patches(i,j,2,1) + self%lims_patches(i,j,2,2)) / 2.
                    idx = (iframe-1) * (NX_PATCHED * NY_PATCHED) + (i-1) * NY_PATCHED + j
                    y(idx) = self%shifts_patches(2,iframe,i,j)   ! shift in x-direction first
                    x(1,idx) = patch_cntr_x
                    x(2,idx) = patch_cntr_y
                    x(3,idx) = real(iframe)
                end do
            end do
        end do
        ! fit polynomial for shifts in x-direction
        call svd_multifit(x,y,sig,a,v,w,chisq,patch_poly)
        ! store polynomial coefficients
        self%poly_coeffs(:,1) = a
        do iframe = 1, self%nframes
            do i = 1, NX_PATCHED
                do j = 1, NY_PATCHED
                    patch_cntr_x = real(self%lims_patches(i,j,1,1) + self%lims_patches(i,j,1,2)) / 2.
                    patch_cntr_y = real(self%lims_patches(i,j,2,1) + self%lims_patches(i,j,2,2)) / 2.
                    idx = (iframe-1) * (NX_PATCHED * NY_PATCHED) + (i-1) * NY_PATCHED + j
                    y(idx) = self%shifts_patches(3,iframe,i,j)   ! shift in y-direction first
                end do
            end do
        end do
        ! fit polynomial for shifts in y-direction
        call svd_multifit(x,y,sig,a,v,w,chisq,patch_poly)
        self%poly_coeffs(:,2) = a
    end subroutine fit_polynomial

    subroutine apply_polytransfo( self, frames, frames_output )
        class(motion_patched),    intent(inout) :: self
        type(image), allocatable, intent(inout) :: frames(:)
        type(image), allocatable, intent(inout) :: frames_output(:)
        integer :: i, j, iframe
        real    :: x, y, t
        real    :: x_trafo, y_trafo
        type(rmat_ptr_type) :: rmat_ins(self%nframes), rmat_outs(self%nframes)
        do iframe = 1, self%nframes
            call frames_output(iframe)%new(self%ldim, params_glob%smpd)
            if (frames(iframe)%is_ft()) call frames(iframe)%ifft()
            call frames(iframe)%get_rmat_ptr(rmat_ins(iframe)%rmat_ptr)
            call frames_output(iframe)%get_rmat_ptr(rmat_outs(iframe)%rmat_ptr)
        end do
        !$omp parallel do collapse(3) default(shared) private(iframe,j,i,x,y,t,x_trafo,y_trafo) proc_bind(close) schedule(static)
        do iframe = 1, self%nframes
            do i = 1, self%ldim(1)
                do j = 1, self%ldim(2)
                    t = real(iframe)
                    x = real(i)
                    y = real(j)
                    x_trafo = apply_patch_poly(self%poly_coeffs(:,1),x,y,t)
                    y_trafo = apply_patch_poly(self%poly_coeffs(:,2),x,y,t)
                    rmat_outs(iframe)%rmat_ptr(i,j,1) = interp_bilin(x_trafo, y_trafo)
                end do
            end do
        end do
        !$omp end parallel do
    contains

        function interp_bilin( xval, yval ) result(val)
            real,  intent(in)    :: xval, yval
            real     :: val
            logical  :: x1_valid, x2_valid, y1_valid, y2_valid
            integer  :: x1_h,  x2_h,  y1_h,  y2_h
            integer  :: x1_hh, x2_hh, y1_hh, y2_hh
            real     :: y1, y2, y3, y4, t, u
            ! if outside of image
!!$        if ((x(1) < 1._dp) .or. (x(1) >= self%ldim_out(1)) .or. (x(2) < 1._dp) .or. (x(2) >= self%ldim_out(2))) then
!!$            val  = 0._dp
!!$            return
!!$        end if
            x1_h = floor(xval)
            x1_hh = x1_h
            if (x1_h < 1) then
                x1_hh = x1_h + self%ldim(1)
            else if (x1_h > self%ldim(1)) then
                x1_hh = x1_h - self%ldim(1)
            end if
            x2_h = x1_h + 1
            x2_hh = x2_h
            if (x2_h < 1) then
                x2_hh = x2_h + self%ldim(1)
            else if (x2_h > self%ldim(1)) then
                x2_hh = x2_h - self%ldim(1)
            end if
            y1_h = floor(yval)
            y1_hh = y1_h
            if (y1_h < 1) then
                y1_hh = y1_h + self%ldim(2)
            else if (y1_h > self%ldim(2)) then
                y1_hh = y1_h - self%ldim(2)
            end if
            y2_h = y1_h + 1
            y2_hh = y2_h
            if (y2_h < 1) then
                y2_hh = y2_h + self%ldim(2)
            else if (y2_h > self%ldim(2)) then
                y2_hh = y2_h - self%ldim(2)
            end if
!!$            y1 = rmat_in(x1_hh, y1_hh, 1)
!!$            y2 = rmat_in(x2_hh, y1_hh, 1)
!!$            y3 = rmat_in(x2_hh, y2_hh, 1)
!!$            y4 = rmat_in(x1_hh, y2_hh, 1)
            y1 = rmat_ins(iframe)%rmat_ptr(x1_hh, y1_hh, 1)
            y2 = rmat_ins(iframe)%rmat_ptr(x2_hh, y1_hh, 1)
            y3 = rmat_ins(iframe)%rmat_ptr(x2_hh, y2_hh, 1)
            y4 = rmat_ins(iframe)%rmat_ptr(x1_hh, y2_hh, 1)

            t    = xval - x1_h
            u    = yval - y1_h
            val  =  (1._dp - t) * (1._dp - u) * y1 + &
               t  * (1._dp - u) * y2 + &
               t  *          u  * y3 + &
                (1._dp - t) *          u  * y4
        end function interp_bilin

    end subroutine apply_polytransfo

    subroutine write_shifts( self )
        class(motion_patched), intent(inout) :: self
        integer :: i,j
        open(123, file='shifts.txt')
        write (123,*) 'shifts_x=[...'
        do i = 1, NX_PATCHED
            do j = 1, NY_PATCHED
                write (123,'(A)',advance='no') real2str(self%shifts_patches(2,1,i,j))
                if (j < NY_PATCHED) write (123,'(A)',advance='no') ', '
            end do
            if (i < NX_PATCHED) write (123,*) '; ...'
        end do
        write (123,*) '];'
        write (123,*) 'shifts_y=[...'
        do i = 1, NX_PATCHED
            do j = 1, NY_PATCHED
                write (123,'(A)',advance='no') real2str(self%shifts_patches(3,1,i,j))
                if (j < NY_PATCHED) write (123,'(A)',advance='no') ', '
            end do
            if (i < NX_PATCHED) write (123,*) '; ...'
        end do
        write (123,*) '];'
        close(123)
    end subroutine write_shifts

    subroutine allocate_fields( self )
        class(motion_patched), intent(inout) :: self
        integer :: alloc_stat
        logical :: do_allocate
        integer :: i,j,iframe
        do_allocate = .true.
        if (allocated(self%shifts_patches)) then
            if (size(self%shifts_patches, dim=1) < self%nframes) then
                do_allocate = .false.
            else
                call self%deallocate_fields()
            end if
        end if
        if (do_allocate) then
            allocate(self%shifts_patches  (3, self%nframes, NX_PATCHED, NY_PATCHED),&
                self%frame_patches   (   self%nframes, NX_PATCHED, NY_PATCHED),&
                self%ref_patches     (   self%nframes, NX_PATCHED, NY_PATCHED),&
                self%patches_imgs    (   self%nframes, NX_PATCHED, NY_PATCHED),&
                self%shsearch_patches(   self%nframes, NX_PATCHED, NY_PATCHED),&
                stat=alloc_stat )
            if (alloc_stat /= 0) call allocchk('allocate_fields 1; simple_motion_patched')
        end if
    end subroutine allocate_fields

    subroutine deallocate_fields( self )
        class(motion_patched), intent(inout) :: self
        integer :: iframe, i, j
        write (*,*) 'deallocate_fields ' ; call flush(6)
        if (allocated(self%shifts_patches)) deallocate(self%shifts_patches)
        if (allocated(self%frame_patches)) then
            do iframe = 1, self%nframes
                do j = 1, NY_PATCHED
                    do i = 1, NX_PATCHED
                        call self%frame_patches(iframe,i,j)%kill()
                    end do
                end do
            end do
            deallocate(self%frame_patches)
        end if
        if (allocated(self%ref_patches)) then
            do iframe = 1, self%nframes
                do j = 1, NY_PATCHED
                    do i = 1, NX_PATCHED
                        call self%ref_patches(iframe,i,j)%kill()
                    end do
                end do
            end do
            deallocate(self%ref_patches)
        end if
        if (allocated(self%patches_imgs)) then
            do iframe = 1, self%nframes
                do j = 1, NY_PATCHED
                    do i = 1, NX_PATCHED
                        call self%patches_imgs(iframe,i,j)%kill()
                    end do
                end do
            end do
            deallocate(self%patches_imgs)
        end if
        if (allocated(self%shsearch_patches)) then
            do iframe = 1, self%nframes
                do j = 1, NY_PATCHED
                    do i = 1, NX_PATCHED
                        call self%shsearch_patches(iframe,i,j)%kill()
                    end do
                end do
            end do
            deallocate(self%shsearch_patches)
        end if
    end subroutine deallocate_fields

    subroutine set_size_frames_ref( self )
        class(motion_patched), intent(inout) :: self
        integer :: ldim_nooverlap(2)
        integer :: i, j
        self%ldim_patch(1) = int(self%ldim(1) / NX_PATCHED) + 2 * X_OVERLAP
        self%ldim_patch(2) = int(self%ldim(2) / NY_PATCHED) + 2 * Y_OVERLAP
        self%ldim_patch(3) = 1
        ldim_nooverlap(1) = int(self%ldim(1) / NX_PATCHED)
        ldim_nooverlap(2) = int(self%ldim(2) / NY_PATCHED)
        do j = 1, NY_PATCHED
            do i = 1, NX_PATCHED
                self%lims_patches(i,j,1,1) = (i-1) * ldim_nooverlap(1) - X_OVERLAP + 1
                self%lims_patches(i,j,1,2) =  i    * ldim_nooverlap(1) + X_OVERLAP
                self%lims_patches(i,j,2,1) = (j-1) * ldim_nooverlap(2) - Y_OVERLAP + 1
                self%lims_patches(i,j,2,2) =  j    * ldim_nooverlap(2) + Y_OVERLAP
            end do
        end do
    end subroutine set_size_frames_ref

    subroutine set_patches( self, stack, patches_ftexp, abc )
        class(motion_patched),          intent(inout) :: self
        type(image),       allocatable, intent(inout) :: stack(:)
        type(ft_expanded), allocatable, intent(inout) :: patches_ftexp(:,:,:)
        integer, intent(in) :: abc
        integer :: i,j,iframe, k, l, kk, ll
        integer :: ip, jp           ! ip, jp: i_patch, j_patch
        integer :: lims_patch(2,2)
        integer :: frame_cd(2,2)    ! coordinates in the system of the frame
        integer :: patch_cd(2,2)    ! coordinates in the system of the patch
        type(rmat_ptr_type) :: rmat_ptrs(self%nframes)
        real, pointer :: rmat_patch(:,:,:)
        do iframe=1,self%nframes
            do j = 1, NY_PATCHED
                do i = 1, NX_PATCHED
                    call self%patches_imgs(iframe,i,j)%new(self%ldim_patch, params_glob%smpd)
                end do
            end do
        end do
        do iframe=1,self%nframes
            call stack(iframe)%get_rmat_ptr(rmat_ptrs(iframe)%rmat_ptr)
        end do
        !$omp parallel do collapse(3) default(shared) private(iframe,j,i,lims_patch,rmat_patch,k,l,kk,ll,ip,jp) proc_bind(close) schedule(static)
        do iframe=1,self%nframes
            do j = 1, NY_PATCHED
                do i = 1, NX_PATCHED
                    call self%patches_imgs(iframe,i,j)%get_rmat_ptr(rmat_patch)
                    lims_patch(:,:) = self%lims_patches(i,j,:,:)
                    do k = lims_patch(1,1), lims_patch(1,2)
                        kk = k
                        if (kk < 1) then
                            kk = kk + self%ldim(1)
                        else if (kk > self%ldim(1)) then
                            kk = kk - self%ldim(1)
                        end if
                        ip = k - lims_patch(1,1) + 1
                        do l = lims_patch(2,1), lims_patch(2,2)
                            ll = l
                            if (ll < 1) then
                                ll = ll + self%ldim(2)
                            else if (ll > self%ldim(2)) then
                                ll = ll - self%ldim(2)
                            end if
                            jp = l - lims_patch(2,1) + 1
                            ! now copy the value
                            rmat_patch(ip,jp,1) = rmat_ptrs(iframe)%rmat_ptr(kk,ll,1)
                        end do
                    end do
                    call patches_ftexp(iframe,i,j)%extract_img(self%patches_imgs(iframe,i,j), self%hp, self%lp)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine set_patches

    subroutine det_shifts( self )
        class(motion_patched), intent(inout) :: self
        integer :: iframe, i, j
        self%shifts_patches = 0.
        !$omp parallel do collapse(3) default(shared) private(iframe,j,i) proc_bind(close) schedule(static)
        do iframe = 1, self%nframes
            do i = 1, NX_PATCHED
                do j = 1, NY_PATCHED
                    !if ((i /= 1).or.(j /= 5)) cycle
                    call self%shsearch_patches(iframe,i,j)%new(self%ref_patches(iframe,i,j),self%frame_patches(iframe,i,j),&
                        self%trs, self%motion_correct_ftol, self%motion_correct_gtol)
                    self%shifts_patches(:,iframe,i,j) = self%shsearch_patches(iframe,i,j)%minimize()
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine det_shifts

    subroutine motion_patched_new( self, motion_correct_ftol, motion_correct_gtol, trs )
        class(motion_patched), intent(inout) :: self
        real, optional,        intent(in)    :: motion_correct_ftol, motion_correct_gtol
        real, optional,        intent(in)    :: trs
        integer :: alloc_stat
        call self%kill()
        if (present(motion_correct_ftol)) then
            self%motion_correct_ftol = motion_correct_ftol
        else
            self%motion_correct_ftol = TOL
        end if
        if (present(motion_correct_gtol)) then
            self%motion_correct_gtol = motion_correct_gtol
        else
            self%motion_correct_gtol = TOL
        end if
        if (present(trs)) then
            self%trs = trs
        else
            self%trs = TRS_DEFAULT
        end if
        self%existence = .true.
    end subroutine motion_patched_new

    subroutine motion_patched_correct( self, hp, lp, references, frames, frames_output )
        class(motion_patched),    intent(inout) :: self
        real,                     intent(in)    :: hp, lp
        type(image), allocatable, intent(inout) :: references(:)
        type(image), allocatable, intent(inout) :: frames(:)
        type(image), allocatable, intent(inout) :: frames_output(:)
        integer                  :: ldim_frames(3)
        integer                  :: i
        self%hp = hp
        self%lp = lp
        self%nframes = size(frames,dim=1)
        self%ldim   = references(1)%get_ldim()
        do i = 1,self%nframes
            ldim_frames = frames(i)%get_ldim()
            if (any(ldim_frames(1:2) /= self%ldim(1:2))) then
                THROW_HARD('error in motion_patched_correct: frame dimensions do not match reference dimension; simple_motion_patched')
            end if
        end do
        call self%allocate_fields()
        call self%set_size_frames_ref()
        ! divide the reference into patches
        call self%set_patches(references, self%ref_patches  ,1)
        call self%set_patches(frames,     self%frame_patches,2)
        ! determine shifts for patches
        call self%det_shifts()
        ! fit the polynomial model against determined shifts
        call self%fit_polynomial()
        ! apply transformation
        call self%apply_polytransfo(frames, frames_output)
        ! write shifts to file
        call self%write_shifts()
    end subroutine motion_patched_correct

    subroutine motion_patched_kill( self )
        class(motion_patched), intent(inout) :: self
        integer :: i, j, iframe
        call self%deallocate_fields()
        self%existence = .false.
    end subroutine motion_patched_kill

end module simple_motion_patched
