!@descr: the class student's alignment after prediction: one particle at its picked position against the teacher's class averages of a list of classes, every in-plane rotation and shifts over the disk, the best alignment kept
! The search is the classic searches' polar shift map (shiftmap_srch): the classes scored at the stored shift, a coarse
! map over the disk of three shift_sigma under its prior, the refinement from each class's best cell; a class keeps its
! refined pose when it beats its stored-shift score, as the classic 2D search commits (doc/algorithms/class_student.md).
module simple_class_student_align
use simple_core_module_api, only: sp, dp
use simple_error,           only: simple_exception
use simple_parameters,      only: parameters
use simple_builder,         only: builder
use simple_shiftmap_srch,   only: shiftmap_srch
implicit none

public :: class_student_align
private
#include "simple_local_flags.inc"

!> one per thread; the shift maps on the calculator are set up once before the parallel region (setup_shiftmap_srch)
type :: class_student_align
    private
    class(builder), pointer :: b_ptr => null()
    type(shiftmap_srch)  :: smap
    real,    allocatable :: vals(:)          !< a class's objective over the rotations at the stored shift
    real,    allocatable :: corrs0(:)        !< per class of the list: the stored shift's best score
    integer, allocatable :: irots0(:)        !< and its rotation
    integer              :: nmax   = 0       !< the longest list align takes
    logical              :: exists = .false.
contains
    procedure :: new
    procedure :: align
    procedure :: kill
end type class_student_align

contains

    !> the search objects for lists of at most nmax classes
    subroutine new( self, params, build, nmax )
        class(class_student_align), intent(inout) :: self
        class(parameters), target,  intent(in)    :: params
        class(builder),    target,  intent(in)    :: build
        integer,                    intent(in)    :: nmax
        call self%kill
        if( nmax < 1 ) THROW_HARD('a list holds one class or more; class_student_align%new')
        self%b_ptr => build
        self%nmax  =  nmax
        call self%smap%new(params, build)
        allocate(self%vals(build%pftc%get_nrots()), self%corrs0(nmax), source=0.)
        allocate(self%irots0(nmax), source=0)
        self%exists = .true.
    end subroutine new

    !> particle iptcl against the classes of refs: the class of the best alignment, its rotation index, score (larger is
    !! better), shift increment from the stored shift in the particle's axes and its map's posterior width sdev (pixels)
    subroutine align( self, iptcl, refs, cls, irot, corr, sh, sdev )
        class(class_student_align), intent(inout) :: self
        integer,                    intent(in)    :: iptcl, refs(:)
        integer,                    intent(out)   :: cls, irot
        real,                       intent(out)   :: corr, sh(2), sdev
        real    :: corr_k, sh_k(2), sdev_k
        integer :: n, k, irot_k
        logical :: found, valid
        n = size(refs)
        if( n < 1 .or. n > self%nmax ) THROW_HARD('a list of one to nmax classes; class_student_align%align')
        call self%smap%prep(iptcl)
        do k = 1, n
            call self%b_ptr%pftc%gen_objfun_vals(refs(k), iptcl, [0.,0.], self%vals)
            self%irots0(k) = maxloc(self%vals, dim=1)
            self%corrs0(k) = self%vals(self%irots0(k))
        enddo
        call self%smap%srch(refs, self%corrs0(1:n), n)
        cls  = 0
        irot = 0
        corr = -huge(corr)
        sh   = 0.
        sdev = -1.
        do k = 1, n
            call self%smap%get_result(refs(k), found, valid, irot_k, corr_k, sh_k, sdev_k)
            if( .not. (found .and. valid) )then
                ! the stored shift at its best rotation beat the map's pose
                irot_k = self%irots0(k)
                corr_k = self%corrs0(k)
                sh_k   = 0.
            endif
            if( corr_k > corr )then
                cls  = refs(k)
                irot = irot_k
                corr = corr_k
                sh   = sh_k
                sdev = sdev_k
            endif
        enddo
    end subroutine align

    subroutine kill( self )
        class(class_student_align), intent(inout) :: self
        if( .not. self%exists ) return
        call self%smap%kill
        deallocate(self%vals, self%corrs0, self%irots0)
        nullify(self%b_ptr)
        self%nmax   = 0
        self%exists = .false.
    end subroutine kill

end module simple_class_student_align
