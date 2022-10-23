module simple_motion_align_poly2
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,           only: image, imstack_type
use simple_ft_expanded_dp,  only: ft_expanded_dp
use simple_parameters,      only: params_glob
use simple_timer
implicit none

public :: motion_align_poly2
private
#include "simple_local_flags.inc"

integer,  parameter :: POLYDIM              = 18         ! One-dimensional size of polynomial
real,     parameter :: RMSD_THRESHOLD_POLY  = 0.1        ! Convergence RMSD for polynomial refinement
real,     parameter :: PBOUND               = 2.0        ! Half bound for polynomial coefficient search
integer,  parameter :: LBFGSB_MAXITS = 60, MAXITS = 4    ! Iterations parameters for polynomial refinement

type :: motion_align_poly2
    private
    type(ft_expanded_dp), allocatable :: ftexp_tiles(:,:,:), ftexp_tiles_sh(:,:)
    type(ft_expanded_dp), allocatable :: ftexp_R(:), ftexp_dR(:,:)
    type(ft_expanded_dp), allocatable :: ftexp_Rhat(:), ftexp_dRhat(:,:), ftexp_Rhat2(:)
    real,              allocatable :: patch_pos(:,:,:)
    real(dp),          allocatable :: patch_coords(:,:,:)
    real,              allocatable :: frameweights(:)                   !< array of frameweights
    real,              allocatable :: corrs(:)                          !< per-frame correlations
    real(dp)                       :: poly_coeffs(POLYDIM,2)            !< optimized polynomial coefficients
    real                           :: hp=-1.,      lp=-1.               !< high/low pass value
    real                           :: bfactor        = 0.              !< b-factor for alignment weights
    real                           :: corr           = -1.              !< correlation
    real                           :: smpd           = 0.               !< sampling distance
    real                           :: scale_factor   = 1.               !< local frame scaling
    real                           :: trs            = 20.              !< half correlation disrete search bound
    integer                        :: nthr           = 1
    integer                        :: nxpatch = 0, nypatch = 0
    integer                        :: ldim(3) = 0
    integer                        :: ldim_tile(3)   = 0
    integer                        :: nframes        = 0                !< number of frames
    integer                        :: fixed_frame    = 1                !< reference frame for interpolation
    logical                        :: l_aniso_success  = .false.
    logical                        :: existence        = .false.

contains
    ! Constructor
    procedure          :: new
    ! Frames & memory management
    procedure          :: refine
    procedure, private :: refine_poly
    procedure, private :: poly_refine_f, poly_refine_fdf
    procedure, private :: calc_shifts
    ! Trajectory fitting related
    procedure, private :: pix2polycoords
    ! Getters & setters
    procedure          :: get_corr
    ! Destructor
    procedure          :: kill
end type motion_align_poly2

contains

    subroutine new( self, frames_patches_ptr, poly_coeffs, patch_centers, ldim, hp, fixed_frame )
        class(motion_align_poly2),         intent(inout)    :: self
        type(imstack_type), allocatable, target, intent(inout) :: frames_patches_ptr(:,:)
        real(dp),                                intent(in) :: poly_coeffs(:,:)
        real,                                    intent(in)    :: patch_centers(:,:,:)
        integer,                                 intent(in) :: ldim(3)
        real,                             intent(in)    :: hp
        integer,                          intent(in)    :: fixed_frame
        integer :: i,j,t,ithr, itmp(2), ldim_patch(3)
        call self%kill
        self%nframes     =  size(frames_patches_ptr(1,1)%stack, 1)
        if ( self%nframes < 2 ) then
            THROW_HARD('nframes < 2; simple_motion_align_poly: align')
        end if
        self%fixed_frame = fixed_frame
        self%smpd        = frames_patches_ptr(1,1)%stack(1)%get_smpd()
        ldim_patch       = frames_patches_ptr(1,1)%stack(1)%get_ldim()
        self%ldim        = ldim
        self%hp          = hp
        self%nthr        = params_glob%nthr
        self%bfactor     = max(0.,params_glob%bfac)
        self%lp          = params_glob%lpstop
        self%nxpatch     = params_glob%nxpatch
        self%nypatch     = params_glob%nypatch
        self%poly_coeffs = poly_coeffs
        allocate(self%corrs(self%nframes), source=-1.)
        ! patch positions
        self%patch_pos = patch_centers
        allocate(self%patch_coords(self%nxpatch,self%nypatch,2))
        do i=1,self%nxpatch
            do j = 1,self%nypatch
                call self%pix2polycoords(real(self%patch_pos(i,j,1),dp),real(self%patch_pos(i,j,2),dp),&
                    &self%patch_coords(i,j,1), self%patch_coords(i,j,2))
            enddo
        enddo
        ! allocations
        allocate(self%ftexp_tiles(self%nframes,self%nxpatch,self%nypatch),&
                &self%ftexp_tiles_sh(self%nframes,self%nthr), self%ftexp_R(self%nthr),&
                &self%ftexp_dR(self%nthr,2), self%ftexp_Rhat(self%nthr), self%ftexp_dRhat(self%nthr,2),&
                &self%ftexp_Rhat2(self%nthr))
        itmp = 0
        !$omp parallel default(shared) private(i,j,t,ithr) proc_bind(close)
        !$omp do collapse(3) schedule(static)
        do t = 1,self%nframes
            do i = 1, self%nxpatch
                do j = 1, self%nypatch
                    call self%ftexp_tiles(t,i,j)%new(frames_patches_ptr(i,j)%stack(t), ldim_patch, self%hp, self%lp, .true., bfac=self%bfactor)
                    call self%ftexp_tiles(t,i,j)%normalize_mat
                    if( itmp(1) == 0 )itmp = [i,j]
                enddo
            enddo
        enddo
        !$omp end do
        !$omp do schedule(static)
        do ithr = 1,self%nthr
            call self%ftexp_R(ithr)%copy( self%ftexp_tiles(1, itmp(1),itmp(2)))
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
        self%existence   = .true.
    end subroutine new


    ! Alignment routines

    subroutine refine( self, poly_coeffs, frameweights )
        class(motion_align_poly2), intent(inout) :: self
        real(dp),                  intent(inout) :: poly_coeffs(POLYDIM,2)
        real,                      intent(inout) :: frameweights(self%nframes)
        call self%refine_poly
        self%corr         = sum(self%corrs) / real(self%nframes)
        self%frameweights = corrs2weights(real(self%corrs), params_glob%wcrit_enum)
        poly_coeffs  = self%poly_coeffs
        frameweights = self%frameweights
    end subroutine refine

    ! Refinement of polynomial coefficients
    subroutine refine_poly( self )
        use simple_opt_factory, only: opt_factory
        use simple_opt_spec,    only: opt_spec
        ! use simple_opt_lbfgsb,  only: PRINT_NEVALS
        use simple_optimizer,   only: optimizer
        class(motion_align_poly2), intent(inout) :: self
        type(opt_factory)         :: ofac
        type(opt_spec)            :: ospec
        class(optimizer), pointer :: nlopt
        real(dp) :: ini_shifts_dp(self%nxpatch,self%nypatch,self%nframes,2), shifts(self%nxpatch,self%nypatch,self%nframes,2)
        real(dp) :: prev_shifts(self%nxpatch,self%nypatch,self%nframes,2)
        real     :: ini_shifts(self%nxpatch,self%nypatch,self%nframes,2)
        real     :: opt_lims(POLYDIM*2,2), lowest_cost, rmsd_cumul, rmsd
        integer  :: t, i, j, iter, ithr, ntot
        ! real(dp) :: vec(POLYDIM*2), grads(POLYDIM*2), f,fl,fr,g,dp
        ntot = self%nframes*self%nxpatch*self%nypatch
        ! initial shifts
        call self%calc_shifts(self%fixed_frame, ini_shifts_dp)
        shifts = ini_shifts_dp
        ! Gradients check
        ! vec(1:POLYDIM)  = self%poly_coeffs(:,1)
        ! vec(POLYDIM+1:) = self%poly_coeffs(:,2)
        ! call self%poly_refine_fdf(vec, f, grads)
        ! dp = 1.0d-7
        ! do i = 1,2*POLYDIM
        !     vec(1:POLYDIM)  = self%poly_coeffs(:,1)
        !     vec(POLYDIM+1:) = self%poly_coeffs(:,2)
        !     vec(i) = vec(i) - dp
        !     fl     = self%poly_refine_f(vec)
        !     vec(1:POLYDIM)  = self%poly_coeffs(:,1)
        !     vec(POLYDIM+1:) = self%poly_coeffs(:,2)
        !     vec(i) = vec(i) + dp
        !     fr     = self%poly_refine_f(vec)
        !     g = (fr-fl) / (2.0d0*dp)
        !     print *,i, f, fl, fr, grads(i), g, grads(i)/g
        ! enddo
        ! stop
        ! Optimization
        ini_shifts = real(ini_shifts_dp)
        do iter = 1,MAXITS
            prev_shifts = shifts
            call minimize
            call self%calc_shifts(self%fixed_frame, shifts)
            ! convergence
            rmsd       = calc_rmsd(real(prev_shifts), real(shifts))
            rmsd_cumul = calc_rmsd(ini_shifts,        real(shifts))
            self%corr  = -lowest_cost / real(ntot)
            if( iter>=2 .and. rmsd<RMSD_THRESHOLD_POLY ) exit
        enddo
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
        !$omp end do
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
                integer :: i,j
                calc_rmsd = 0.
                do i = 1,self%nxpatch
                    do j = 1,self%nypatch
                        calc_rmsd = calc_rmsd + sum((sh1(i,j,:,:)-sh2(i,j,:,:))**2.)
                    enddo
                enddo
                calc_rmsd = sqrt(calc_rmsd/real(ntot))
            end function calc_rmsd

    end subroutine refine_poly

    ! polynomial model -> tiles shifts
    subroutine calc_shifts( self, iframe, sh )
        class(motion_align_poly2), intent(in)  :: self
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

    ! Getters/setters

    real function get_corr( self )
        class(motion_align_poly2), intent(in) :: self
        get_corr = self%corr
    end function get_corr

    ! FITTING RELATED

    subroutine pix2polycoords( self, xin, yin, x, y )
        class(motion_align_poly2), intent(inout) :: self
        real(dp),                intent(in)    :: xin, yin
        real(dp),                intent(inout) :: x, y
        x = (xin-1.d0) / real(self%ldim(1)-1,dp) - 0.5d0
        y = (yin-1.d0) / real(self%ldim(2)-1,dp) - 0.5d0
    end subroutine pix2polycoords

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

    subroutine kill( self )
        class(motion_align_poly2), intent(inout) :: self
        self%hp      = -1.
        self%lp      = -1.
        self%bfactor = -1.
        self%corr        = -1.
        self%smpd        = 0.
        self%nthr        = 1
        self%ldim        = 0
        self%fixed_frame = 1
        self%poly_coeffs     = 0.d0
        if(allocated(self%patch_pos))       deallocate(self%patch_pos)
        if(allocated(self%patch_coords))    deallocate(self%patch_coords)
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
        class is(motion_align_poly2)
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
        class is(motion_align_poly2)
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
        class is(motion_align_poly2)
            call self%poly_refine_fdf( vec, f, grad )
        class DEFAULT
            THROW_HARD('unknown type; patched_direct_fdf_wrapper')
        end select
    end subroutine poly_refine_df_wrapper

    real(dp) function poly_refine_f( self, vec )
        class(motion_align_poly2), intent(inout) :: self
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
                    rt = real(t-self%fixed_frame, dp)
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
                !!!!!!!!!!!!
                ! if( i == 1 .or. i == self%nxpatch .or. j == 1 .or. j == self%nypatch )then
                !     ccs = ccs + 2.d0*(Es-1.d0)/Ds
                ! else
                !     ccs = ccs + (Es-1.d0)/Ds
                ! endif
                !!!!!!!!!!!!!!!!
            end do
        end do
        !$omp end parallel do
        poly_refine_f = -sum(ccs)
        self%corrs = real(ccs) / real(self%nxpatch*self%nypatch)
    end function poly_refine_f

    subroutine poly_refine_fdf( self, vec, f, grads )
        class(motion_align_poly2), intent(inout) :: self
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
                    rt = real(t-self%fixed_frame, dp)
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
                !!!!!!!!!!
                ! if( i == 1 .or. i == self%nxpatch .or. j == 1 .or. j == self%nypatch )then
                !     ccs = ccs + 2.d0*(Es-1.d0)/Ds
                ! else
                !     ccs = ccs + (Es-1.d0)/Ds
                ! endif
                !!!!!!!!!!!!!!!
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
                        if( t == self%fixed_frame ) cycle  ! w = 0.
                        w = real(t-self%fixed_frame,dp)**p
                        call self%ftexp_Rhat2(ithr)%add( self%ftexp_tiles_sh(t,ithr), w=w)
                    end do
                    g(p)   = self%ftexp_dRhat(ithr,1)%corr_unnorm( self%ftexp_Rhat2(ithr) )
                    g(3+p) = self%ftexp_dRhat(ithr,2)%corr_unnorm( self%ftexp_Rhat2(ithr) )
                    call self%ftexp_Rhat2(ithr)%zero
                    do t = 1,self%nframes
                        if( t == self%fixed_frame ) cycle  ! w = 0.
                        w = real(t-self%fixed_frame,dp)**p * (sumF - Fs(t) - 1.d0 / Ds(t))
                        call self%ftexp_Rhat2(ithr)%add( self%ftexp_tiles_sh(t,ithr), w=w)
                    end do
                    g(p)   = g(p)   - self%ftexp_dR(ithr,1)%corr_unnorm( self%ftexp_Rhat2(ithr) )
                    g(3+p) = g(3+p) - self%ftexp_dR(ithr,2)%corr_unnorm( self%ftexp_Rhat2(ithr) )
                end do
                !!!!!!!!!
                ! if( i == 1 .or. i == self%nxpatch .or. j == 1 .or. j == self%nypatch ) g = g*2.d0
                !!!!!!!!!!!!!
                grads( 1:18) = grads( 1:18) - [g(1:3), x*g(1:3), x*x*g(1:3), y*g(1:3), y*y*g(1:3), x*y*g(1:3)]
                grads(19:36) = grads(19:36) - [g(4:6), x*g(4:6), x*x*g(4:6), y*g(4:6), y*y*g(4:6), x*y*g(4:6)]
            end do
        end do
        !$omp end parallel do
        self%corrs = real(ccs) / real(self%nxpatch*self%nypatch)
        f = -sum(ccs)
    end subroutine poly_refine_fdf

end module simple_motion_align_poly2
