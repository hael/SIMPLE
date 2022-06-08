module simple_eval_cartftcc
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,           only: image
use simple_projector,       only: projector
use simple_parameters,      only: params_glob
use simple_oris,            only: oris
use simple_cartft_corrcalc, only: cartftcc_glob
use simple_parameters,      only: params_glob
implicit none

public :: eval_cartftcc
private
#include "simple_local_flags.inc"

type :: eval_cartftcc
    private
    integer                  :: nspace  = 0
    integer                  :: ldim(3) = 0
    type(projector), pointer :: vol_even => null(), vol_odd => null() ! prepared e/o vols
    type(image), allocatable :: projs(:) ! heap var for thread safety
    type(oris)               :: orispace
    logical                  :: exists = .false.
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! SETTERS
    procedure          :: set_eo_ptrs
    procedure          :: set_ori
    ! CALCULATORS      
    procedure          :: project_and_correlate
    ! DESTRUCTOR
    procedure          :: kill
end type eval_cartftcc

contains

    ! CONSTRUCTOR

    subroutine new( self, vol_even, vol_odd, nspace )
        class(eval_cartftcc),    intent(inout) :: self
        type(projector), target, intent(in)    :: vol_even, vol_odd
        integer,                 intent(in)    :: nspace
        integer :: ithr
        call self%set_eo_ptrs(vol_even, vol_odd)
        self%nspace = nspace
        self%ldim   = vol_even%get_ldim()
        call self%orispace%new(self%nspace, is_ptcl=.false.)
        allocate(self%projs(nthr_glob))
        do ithr = 1,nthr_glob
            call self%projs(ithr)%new([self%ldim(1),self%ldim(2),1], params_glob%smpd)
        end do
        self%exists = .true.
    end subroutine new

    ! SETTERS

    subroutine set_eo_ptrs( self, vol_even, vol_odd )
        class(eval_cartftcc),    intent(inout) :: self
        type(projector), target, intent(in)    :: vol_even, vol_odd
        if( .not. vol_even%is_expanded() ) THROW_HARD('input vol_even expected to be prepared for interpolation') 
        if( .not. vol_odd%is_expanded()  ) THROW_HARD('input vol_odd expected to be prepared for interpolation')
        self%vol_even => vol_even
        self%vol_odd  => vol_odd
    end subroutine set_eo_ptrs

    subroutine set_ori( self, i, euls, shvec )
        class(eval_cartftcc), intent(inout) :: self
        integer,              intent(in)    :: i
        real,                 intent(in)    :: euls(3), shvec(2)
        if( i < 1 .or. i > self%nspace ) THROW_HARD('index i out of range')
        call self%orispace%set_euler(i, euls)
        call self%orispace%set_shift(i, shvec)
    end subroutine set_ori

    ! CALCULATORS

    subroutine project_and_correlate( self, iptcl, corrs )
        class(eval_cartftcc), intent(inout) :: self
        integer,              intent(in)    :: iptcl
        real,                 intent(inout) :: corrs(self%nspace)
        type(projector), pointer :: vol_ptr => null()
        integer :: iref, ithr
        logical :: iseven
        iseven = cartftcc_glob%ptcl_iseven(iptcl)
        if( iseven )then
            vol_ptr => self%vol_even
        else
            vol_ptr => self%vol_odd
        endif
        ithr = omp_get_thread_num() + 1 ! needs to be moved into the loop if we want to parallelize here
        do iref = 1,self%nspace
            call vol_ptr%fproject_serial(self%orispace, iref, self%projs(ithr), params_glob%kstop)
            call cartftcc_glob%set_ref(iref, self%projs(ithr), iseven)
            corrs(iref) = cartftcc_glob%calc_corr(iref, iptcl, self%orispace%get_2Dshift(iptcl))
        end do
    end subroutine project_and_correlate

    ! DESTRUCTOR

    subroutine kill( self )
        class(eval_cartftcc), intent(inout) :: self
        integer :: iref
        if( self%exists )then
            self%nspace = 0
            self%ldim   = 0
            self%vol_even => null()
            self%vol_odd  => null()
            call self%orispace%kill
            do iref = 1,self%nspace
                call self%projs(iref)%kill
            end do
            self%exists = .false.
        endif
    end subroutine kill

end module simple_eval_cartftcc
