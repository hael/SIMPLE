module simple_tvlam_opt
!$ use omp_lib
!$ use omp_lib_kinds
use simple_opt_spec,    only: opt_spec
use simple_opt_simplex, only: opt_simplex
use simple_image,       only: image
use simple_tvfilter,    only: tvfilter
implicit none

public :: tvlam_opt
private
#include "simple_local_flags.inc"

type tvlam_opt
    private
    type(opt_spec)        :: ospec          !< optimizer specification object
    type(opt_simplex)     :: nlopt          !< optimizer object
    class(image), pointer :: img1 => null() !< img1 to compare with img2
    class(image), pointer :: img2 => null() !< img2 to compare with img1
    type(image)           :: img1tv         !< tv regularized version of img1
    type(image)           :: img2tv         !< tv regularized version of img2
    type(tvfilter)        :: tvfilt         !< TV filter instance
    logical, allocatable  :: lmsk(:,:,:)    !< logical mask for distance calc
    logical               :: is3D = .false. !< indicates whether 3D/2D
contains
    procedure :: new
    procedure :: set_img_ptrs
    procedure :: minimize
    procedure :: kill
end type tvlam_opt

real, parameter :: lam_bounds(2) = [0.5,5.0]

contains

    subroutine new( self, ldim, smpd, msk )
        class(tvlam_opt), intent(inout) :: self
        integer, intent(in) :: ldim(3)
        real,    intent(in) :: smpd, msk
        type(image) :: mskimg
        real        :: lims(1,2)
        ! logical mask for distance calculation
        call mskimg%disc(ldim, smpd, msk, self%lmsk)
        ! tv regularized versions of img1/img2
        call img1tv%new(ldim, smpd)
        call img2tv%new(ldim, smpd)
        ! optimizer specification
        lims(1,1) = lam_bounds(1)
        lims(1,2) = lam_bounds(2)
        call self%ospec%specify('simplex', 1, ftol=1e-4, gtol=1e-4, limits=lims, nrestarts=3, maxits=100)
        call self%ospec%set_costfun(tvlam_costfun)
        call self%nlopt%new(self%ospec)
        ! make TV filter
        call self%tvfilt%new
        ! set 2D/3D flag
        self%is3D = ldim(3) /= 1
        call mskimg%kill
    end subroutine new

    subroutine set_img_ptrs( self, img1, img2 )
        class(tvlam_opt),     intent(inout) :: self
        class(image), target, intent(inout) :: img1, img2
        self%img1 => img1
        self%img2 => img2
    end subroutine set_img_ptrs

    subroutine minimize( self, lam )
        class(tvlam_opt), intent(inout) :: self
        real,             intent(inout) :: lam
        real, parameter :: stepsz = 0.2
        real :: lam_trial(1), cost, cost_best, cost_init
        cost_best    = huge(lam_trial(1))
        lam_trial(1) = lam_bounds(1)
        lam          = lam_trial(1)
        do while( lam_trial(1) < lam_bounds(2) )
            cost = self%tvlam_costfun(lam_trial, 1)

            print *, 'lambda / cost: ', lam_trial(1), cost

            if( cost < cost_best )then
                lam = lam_trial(1)
                cost_best = cost
            endif
            lam_trial(1) = lam_trial(1) + stepsz
        end do
        lam_trial(1)    = lam
        self%ospec%x(1) = lam
        cost_init = self%tvlam_costfun(lam_trial, 1)
        call self%nlopt%minimize(self%ospec, self, cost)
        if( cost < cost_init )then
            lam = self%ospec%x(1) 
            cost_best = cost

            print *, 'simplex minimization identified a better solution'
            print *, 'lambda / cost: ', lam, cost_best

        endif
    end subroutine minimize

    subroutine kill( self )
        class(tvlam_opt), intent(inout) :: self
        call self%ospec%kill
        call nlopt%kill
        self%img1 => null()
        self%img2 => null()
        call self%img1tv%kill
        call self%img2tv%kill
        call self%tvfilt%kill
        if( allocated(self%lmsk) ) deallocate(self%lmsk)
    end subroutine kill

    ! accessory functions

    function tvlam_costfun( self, vec, D ) result( cost )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real :: cost
        call img1tv%copy(self%img1)
        call img2tv%copy(self%img2) 
        if( self%is3D )then
            call tvfilt%apply_filter_3d(img1tv, vec(1)) ! vec(1) is lambda
            call tvfilt%apply_filter_3d(img2tv, vec(1))
        else
            call tvfilt%apply_filter(img1tv, vec(1))    ! vec(1) is lambda
            call tvfilt%apply_filter(img2tv, vec(1))
        endif
        select type( self )
            class is (tvlam_opt)
                cost = img1%sqeuclid(img2, mask=self%lmsk)
            class default
                THROW_HARD('error in tvlam_costfun: unknown type')
        end select
    end function tvlam_costfun

end module simple_tvlam_opt
