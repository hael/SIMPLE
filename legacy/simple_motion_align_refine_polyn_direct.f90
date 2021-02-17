! align frames in a maximum-likelihood way, with direct optimization of polynomial coefficients
! no frameweights are being applied here
module simple_motion_align_refine_polyn_direct
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_error
use simple_ft_expanded,           only: ft_expanded
use simple_image,                 only: image
use simple_opt_lbfgsb,            only: PRINT_NEVALS
implicit none
public :: align_refine_polyn_direct, POLYDIM, POLYDIM2
private
#include "simple_local_flags.inc"

integer,  parameter :: POLYDIM  = 3
integer,  parameter :: POLYDIM2 = 2*POLYDIM

type :: align_refine_polyn_direct
    ! private
    type(image), pointer                              :: frames(:)                         !< pointer to stack of frames
    type(ft_expanded), allocatable                    :: movie_frames_Ij(:)                !< normalised movie frames
    type(ft_expanded), allocatable                    :: movie_frames_Ij_sh(:)             !< shifted movie frames
    type(ft_expanded), allocatable                    :: movie_frames_R                    !< reference
    type(ft_expanded), allocatable                    :: movie_frames_dR12(:)              !< gradient of movie_frames_R
    type(ft_expanded), allocatable                    :: movie_frames_Rhat                 !< reference movie frame
    type(ft_expanded), allocatable                    :: movie_frames_dRhat12(:)           !< gradient of reference movie frame
    type(ft_expanded), allocatable                    :: movie_frames_Rhat2                !< reference movie frame
    type(ft_expanded), allocatable                    :: movie_frames_Rhat3                !< reference movie frame
    real(dp),          allocatable, public            :: shifts(:,:)                       !< shifts, to be updated externally
    real(dp),          allocatable                    :: corrs(:)                          !< per-frame correlations
    real                                              :: hp = -1., lp = -1.                !< resolution limits
    real                                              :: bfactor   = -1.                   !< b-factor for alignment weights
    real                                              :: corr      = -1.                   !< average correlation
    integer                                           :: nframes   = 0                     !< number of frames
    integer                                           :: coord_x   = 0, coord_y = 0        !< x,y coordinates for patch-based alignment callback
    integer                                           :: fixed_frame = 1
    logical                                           :: existence = .false.
contains
    ! Constructor
    procedure,                                private :: new_1, new_2
    generic                                           :: new => new_1, new_2
    procedure                                         :: create_ftexp_objs
    procedure                                         :: transfer_ftexp_obj
    procedure                                         :: create_tmp_ftexp_objs
    ! Calculators
    procedure                                         :: contribs_cost
    procedure                                         :: contribs_gcost
    procedure                                         :: contribs_fdf
    ! Getters/Setters
    procedure                                         :: set_hp_lp
    procedure                                         :: set_fixed_frame
    procedure                                         :: set_coords
    procedure                                         :: set_bfactor
    procedure                                         :: get_corrs
    ! Destructor
    procedure                                         :: kill
end type align_refine_polyn_direct

contains

    subroutine new_1( self, frames_ptr, nframes )
        class(align_refine_polyn_direct), intent(inout) :: self
        type(image), allocatable, target, intent(in)    :: frames_ptr(:)
        integer,                          intent(in)    :: nframes
        integer :: alloc_stat
        call self%kill
        self%frames  => frames_ptr
        if (size(frames_ptr, 1) /= nframes) then
            THROW_HARD('nframes > #frames provided; simple_align_refine_polyn_direct: new')
        end if
        self%nframes = nframes
        allocate(self%movie_frames_Ij(self%nframes),  self%movie_frames_Ij_sh(self%nframes),                &
            &self%movie_frames_R, self%movie_frames_Rhat, self%movie_frames_Rhat2, self%movie_frames_Rhat3, &
            &self%movie_frames_dR12(2), self%movie_frames_dRhat12(2), self%shifts(self%nframes,2),          &
            &self%corrs(self%nframes), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('new; simple_align_refine_polyn_direct')
        self%existence = .true.
    end subroutine new_1

    subroutine new_2( self, nframes, frame )
        class(align_refine_polyn_direct), intent(inout) :: self
        integer,                          intent(in)    :: nframes
        class(image),                     intent(inout) :: frame
        integer :: alloc_stat
        call self%kill
        self%nframes = nframes
        allocate(self%movie_frames_Ij(self%nframes),  self%movie_frames_Ij_sh(self%nframes),                &
            &self%movie_frames_R, self%movie_frames_Rhat, self%movie_frames_Rhat2, self%movie_frames_Rhat3, &
            &self%movie_frames_dR12(2), self%movie_frames_dRhat12(2), self%shifts(self%nframes,2),          &
            &self%corrs(self%nframes), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('new; simple_align_refine_polyn_direct')
        self%existence = .true.
    end subroutine new_2

    subroutine transfer_ftexp_obj( self, frame, t )
        class(align_refine_polyn_direct), intent(inout) :: self
        class(image),                     intent(inout) :: frame
        integer,                          intent(in)    :: t
        call frame%fft
        call self%movie_frames_Ij(t)%new(frame, self%hp, self%lp, .true., bfac=self%bfactor)
        call self%movie_frames_Ij(t)%normalize_mat
        call self%movie_frames_Ij_sh(t)%copy(self%movie_frames_Ij(t))
        call self%movie_frames_Ij_sh(t)%zero
    end subroutine transfer_ftexp_obj

    subroutine create_tmp_ftexp_objs( self )
        class(align_refine_polyn_direct), intent(inout) :: self
        call self%movie_frames_R      %copy(self%movie_frames_Ij_sh(1))
        call self%movie_frames_Rhat   %copy(self%movie_frames_R)
        call self%movie_frames_Rhat2  %copy(self%movie_frames_R)
        call self%movie_frames_Rhat3  %copy(self%movie_frames_R)
        call self%movie_frames_dR12(1)%copy(self%movie_frames_R)
        call self%movie_frames_dR12(2)%copy(self%movie_frames_R)
        call self%movie_frames_dRhat12(1)%copy(self%movie_frames_R)
        call self%movie_frames_dRhat12(2)%copy(self%movie_frames_R)
    end subroutine create_tmp_ftexp_objs

    subroutine create_ftexp_objs( self )
        class(align_refine_polyn_direct), intent(inout) :: self
        integer :: j
        call self%movie_frames_R%new(self%frames(1), self%hp, self%lp, .false.)
        do j=1,self%nframes
            call self%frames(j)%fft
            call self%movie_frames_Ij(j)%new(self%frames(j), self%hp, self%lp, .true., bfac=self%bfactor)
            call self%movie_frames_Ij(j)%normalize_mat
            call self%movie_frames_Ij_sh(j)%copy(self%movie_frames_R)
        end do
        call self%movie_frames_Rhat   %copy(self%movie_frames_R)
        call self%movie_frames_Rhat2  %copy(self%movie_frames_R)
        call self%movie_frames_Rhat3  %copy(self%movie_frames_R)
        call self%movie_frames_dR12(1)%copy(self%movie_frames_R)
        call self%movie_frames_dR12(2)%copy(self%movie_frames_R)
        call self%movie_frames_dRhat12(1)%copy(self%movie_frames_R)
        call self%movie_frames_dRhat12(2)%copy(self%movie_frames_R)
    end subroutine create_ftexp_objs

    ! CALCULATORS

    real(dp) function contribs_cost( self )
        class(align_refine_polyn_direct), intent(inout) :: self
        integer  :: j
        real(dp) :: val_e (self%nframes), val_D (self%nframes), shvec(2), val_RR
        ! shifts (need to be updated externally) & calculate reference
        call self%movie_frames_R%zero
        do j = 1, self%nframes
            shvec = self%shifts(j,:)
            call self%movie_frames_Ij(j)%shift( real(shvec), self%movie_frames_Ij_sh(j), calc_sumsq=.false. )
            call self%movie_frames_R%add_uncond( self%movie_frames_Ij_sh(j) )
        end do
        ! calculate val_e
        do j = 1, self%nframes
            val_e(j) = self%movie_frames_R%corr_unnorm_serial( self%movie_frames_Ij_sh(j) )
        end do
        val_RR        = sum(val_e)
        val_D         = sqrt(val_RR - 2._dp * val_e + 1.0_dp)
        self%corrs    = (val_e-1.d0) / val_D
        contribs_cost = -sum( self%corrs )
        self%corr     = -real(contribs_cost)/real(self%nframes)
    end function contribs_cost

    subroutine contribs_gcost( self, grad )
        class(align_refine_polyn_direct), intent(inout) :: self
        real(dp),                         intent(out)   :: grad(POLYDIM2)
        real(dp) :: f
        call self%contribs_fdf( f, grad )
    end subroutine contribs_gcost

    subroutine contribs_fdf( self, f, grad )
        class(align_refine_polyn_direct), intent(inout) :: self
        real(dp),                         intent(out)   :: f
        real(dp),                         intent(out)   :: grad(POLYDIM2)
        integer  :: j, p
        real(dp) :: val_e (self%nframes), val_F (self%nframes), val_D (self%nframes), val_EE(self%nframes)
        real(dp) :: shvec(2), val_RR, tmp1, tmp2, wp, wj
        call self%movie_frames_R%zero
        ! shifts (need to be updated externally) & calculate reference
        do j = 1, self%nframes
            shvec = self%shifts(j,:)
            call self%movie_frames_Ij(j)%shift(real(shvec), self%movie_frames_Ij_sh(j), calc_sumsq=.false.)
            call self%movie_frames_R%add_uncond( self%movie_frames_Ij_sh(j) )
        end do
        ! calculate val_e
        do j = 1, self%nframes
            val_e(j) = self%movie_frames_R%corr_unnorm_serial( self%movie_frames_Ij_sh(j) )
        end do
        val_RR = sum(val_e)
        ! Denominator terms
        val_D  = sqrt(val_RR - 2._dp * val_e + 1.0_dp)
        val_EE = val_D**3
        val_F  = (val_e - 1._dp) / val_EE
        call self%movie_frames_Rhat%zero
        do j = 1, self%nframes
            call self%movie_frames_Rhat%add_uncond(self%movie_frames_Ij_sh(j), real( 1._dp / val_D(j)  ) )
        end do
        ! Cost
        f = sum( (val_e-1.d0) / val_D )
        ! calc gradient of Rhat, R
        call self%movie_frames_Rhat%gen_grad_noshift( self%movie_frames_dRhat12(1), self%movie_frames_dRhat12(2) )
        call self%movie_frames_R   %gen_grad_noshift( self%movie_frames_dR12(1),    self%movie_frames_dR12(2)    )
        do p = 1, POLYDIM
            call self%movie_frames_Rhat2%zero
            call self%movie_frames_Rhat3%zero
            do j = 2, self%nframes
                wp = real(j-1,dp)**p
                wj = sum(val_F) - val_F(j) - 1._dp / val_D(j)
                call self%movie_frames_Rhat2%add_uncond(self%movie_frames_Ij_sh(j), real(wp) )
                call self%movie_frames_Rhat3%add_uncond(self%movie_frames_Ij_sh(j), real(wp * wj))
            end do
            tmp1 = self%movie_frames_dRhat12(1)%corr_unnorm_serial( self%movie_frames_Rhat2 )
            tmp2 = self%movie_frames_dR12(1)   %corr_unnorm_serial( self%movie_frames_Rhat3 )
            grad(p) = tmp1 - tmp2
            tmp1 = self%movie_frames_dRhat12(2)%corr_unnorm_serial( self%movie_frames_Rhat2 )
            tmp2 = self%movie_frames_dR12(2)   %corr_unnorm_serial( self%movie_frames_Rhat3 )
            grad(POLYDIM+p) = tmp1 - tmp2
        end do
        f    = - f
        grad = - grad
    end subroutine contribs_fdf

    ! GETTERS/SETTERS

    subroutine set_hp_lp( self, hp, lp )
        class(align_refine_polyn_direct), intent(inout) :: self
        real,                             intent(in)    :: hp, lp
        if (.not. self%existence) then
            THROW_HARD('not instantiated; simple_align_refine_polyn_direct: set_hp_lp')
        end if
        self%hp = hp
        self%lp = lp
    end subroutine set_hp_lp

    subroutine set_fixed_frame( self, iframe )
        class(align_refine_polyn_direct), intent(inout) :: self
        integer,                             intent(in) :: iframe
        self%fixed_frame = iframe
    end subroutine set_fixed_frame

    subroutine set_coords( self, x, y )
        class(align_refine_polyn_direct), intent(inout) :: self
        integer,                             intent(in) :: x, y
        self%coord_x = x
        self%coord_y = y
    end subroutine set_coords

    subroutine set_bfactor( self, bfac )
        class(align_refine_polyn_direct), intent(inout) :: self
        real,                             intent(in)    :: bfac
        self%bfactor = bfac
    end subroutine set_bfactor

    subroutine get_corrs( self, corrs )
        class(align_refine_polyn_direct), intent(in)  :: self
        real(dp),                         intent(out) :: corrs(self%nframes)
        corrs = self%corrs
    end subroutine get_corrs

    ! DESTRUCTOR

    subroutine kill( self )
        class(align_refine_polyn_direct), intent(inout) :: self
        integer :: iframe
        self%hp      = -1.
        self%lp      = -1.
        self%bfactor = -1.
        self%corr    = -1.
        self%coord_x = 0
        self%coord_y = 0
        self%fixed_frame = 1
        if( allocated(self%movie_frames_Ij) )then
            do iframe = 1,self%nframes
                call self%movie_frames_Ij(iframe)%kill
                call self%movie_frames_Ij_sh(iframe)%kill
            enddo
            deallocate(self%movie_frames_Ij,  self%movie_frames_Ij_sh, self%movie_frames_R, &
                &self%movie_frames_Rhat, self%movie_frames_Rhat2, self%movie_frames_Rhat3,  &
                &self%movie_frames_dR12, self%movie_frames_dRhat12, self%shifts, self%corrs)
            self%nframes   = 0
        endif
        self%existence = .false.
    end subroutine kill

end module simple_motion_align_refine_polyn_direct
