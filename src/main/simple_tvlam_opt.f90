module simple_tvlam_opt
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_image,    only: image
use simple_tvfilter, only: tvfilter
implicit none

public :: tvlam_opt
private
#include "simple_local_flags.inc"

type tvlam_opt
    private
    class(image), pointer :: img_e => null() !< img_e to compare with img_o
    class(image), pointer :: img_o => null() !< img_o to compare with img_e
    type(image)           :: img_e_tv        !< tv regularized version of img_e
    type(image)           :: img_o_tv        !< tv regularized version of img_o
    type(tvfilter)        :: tvfilt          !< TV filter instance
    logical, allocatable  :: lmsk(:,:,:)     !< logical mask for distance calc
    logical               :: is3D = .false.  !< indicates whether 3D/2D
contains
    procedure :: new
    procedure :: set_img_ptrs
    procedure :: get_tvfiltered
    procedure :: minimize
    procedure :: kill
end type tvlam_opt

real, parameter :: lam_bounds(2) = [0.5,5.0]
real, parameter :: CORR_THRES = 0.96

contains

    subroutine new( self, ldim, smpd, msk )
        class(tvlam_opt), intent(inout) :: self
        integer, intent(in) :: ldim(3)
        real,    intent(in) :: smpd, msk
        type(image) :: mskimg
        ! logical mask for distance calculation
        call mskimg%disc(ldim, smpd, msk, self%lmsk)
        ! tv regularized versions of img_e/img_o
        call self%img_e_tv%new(ldim, smpd)
        call self%img_o_tv%new(ldim, smpd)
        ! make TV filter
        call self%tvfilt%new
        ! set 2D/3D flag
        self%is3D = ldim(3) /= 1
        call mskimg%kill
    end subroutine new

    subroutine set_img_ptrs( self, img_e, img_o )
        class(tvlam_opt),     intent(inout) :: self
        class(image), target, intent(inout) :: img_e, img_o
        self%img_e => img_e
        self%img_o => img_o
    end subroutine set_img_ptrs

    subroutine get_tvfiltered( self, img, is_even )
        class(tvlam_opt), intent(inout) :: self
        class(image),     intent(inout) :: img
        logical, optional,   intent(in) :: is_even
        logical :: l_is_even
        l_is_even = .true.
        if( present(is_even) ) l_is_even = is_even
        if( l_is_even )then
            call img%copy(self%img_e_tv)
        else
            call img%copy(self%img_o_tv)
        endif
    end subroutine get_tvfiltered

    subroutine minimize( self, lam )
        class(tvlam_opt), intent(inout) :: self
        real,             intent(inout) :: lam
        real, parameter :: stepsz = 0.1
        real            :: lam_trial(1), corr, corr_best, corr_init
        logical         :: found_opt
        corr_best    = huge(lam_trial(1))
        lam_trial(1) = lam_bounds(1)
        lam          = lam_trial(1)
        found_opt    = .false.
        do while( lam_trial(1) < lam_bounds(2) )
            corr = tvlam_corr(self, lam_trial, 1)
            print *, corr
            if( corr > CORR_THRES )then
                lam = lam_trial(1)
                corr_best = corr
                found_opt = .true.
                exit
            endif
            lam_trial(1) = lam_trial(1) + stepsz
        end do
        if( .not. found_opt ) lam = lam_bounds(2)
        print *, 'lambda: ', lam
    end subroutine minimize

    subroutine kill( self )
        class(tvlam_opt), intent(inout) :: self
        self%img_e => null()
        self%img_o => null()
        call self%img_e_tv%kill
        call self%img_o_tv%kill
        call self%tvfilt%kill
        if( allocated(self%lmsk) ) deallocate(self%lmsk)
    end subroutine kill

    function tvlam_corr( self, vec, D ) result( corr )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real :: corr
        select type( self )
            class is (tvlam_opt)
                call self%img_e_tv%copy(self%img_e)
                call self%img_o_tv%copy(self%img_o)
                if( self%is3D )then
                    call self%tvfilt%apply_filter_3d(self%img_e_tv, vec(1)) ! vec(1) is lambda
                    call self%tvfilt%apply_filter_3d(self%img_o_tv, vec(1))
                else
                    call self%tvfilt%apply_filter(self%img_e_tv, vec(1))    ! vec(1) is lambda
                    call self%tvfilt%apply_filter(self%img_o_tv, vec(1))
                endif
                corr = self%img_e_tv%real_corr(self%img_o_tv, mask=self%lmsk)
            class default
                THROW_HARD('error in tvlam_corrfun: unknown type')
        end select
    end function tvlam_corr

end module simple_tvlam_opt
