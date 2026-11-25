!@descr: the pose student's alignment after prediction: one particle at its stored shift against the map's reprojections of a list of grid directions and their mirror partners, every in-plane rotation and shifts over the disk, the best alignment kept
! The search is the classic searches' polar shift map (shiftmap_srch), as the class student's alignment
! (simple_class_student_align): the references scored at the stored shift, a coarse map over the disk of three
! shift_sigma under its prior, the refinement from each reference's best cell; a reference keeps its refined pose when
! it beats its stored-shift score. The pose student's cone names one direction of each mirror pair: both are searched
! (doc/algorithms/pose_student.md).
module simple_pose_student_align
use simple_core_module_api, only: sp
use simple_error,           only: simple_exception
use simple_parameters,      only: parameters
use simple_builder,         only: builder
use simple_shiftmap_srch,   only: shiftmap_srch
implicit none

public :: pose_student_align, pose_align_result
private
#include "simple_local_flags.inc"

!> an alignment's result: the reference (the grid direction), its rotation index, score (larger is better), the shift
!! increment from the stored shift in the particle's axes and its map's posterior width (pixels)
type :: pose_align_result
    integer  :: iref  = 0
    integer  :: irot  = 0
    real(sp) :: corr  = 0.
    real(sp) :: sh(2) = 0.
    real(sp) :: sdev  = -1.
end type pose_align_result

!> one per thread; the shift maps on the calculator are set up once before the parallel region (setup_shiftmap_srch)
type :: pose_student_align
    private
    class(builder), pointer     :: b_ptr  => null()
    type(shiftmap_srch)         :: smap
    real,           allocatable :: vals(:)    !< a reference's objective over the rotations at the stored shift
    real,           allocatable :: corrs0(:)  !< per reference of the list: the stored shift's best score
    integer,        allocatable :: irots0(:)  !< and its rotation
    integer,        allocatable :: refs(:)    !< the list: the directions and their partners, each once
    integer                     :: nmax   = 0 !< the longest list align takes
    logical                     :: exists = .false.
contains
    procedure :: new
    procedure :: align
    procedure :: kill
end type pose_student_align

contains

    !> the search objects for lists of at most nmax references (directions and partners together)
    subroutine new( self, params, build, nmax )
        class(pose_student_align), intent(inout)      :: self
        class(parameters),         target, intent(in) :: params
        class(builder),            target, intent(in) :: build
        integer,                   intent(in)         :: nmax
        call self%kill
        if( nmax < 1 ) THROW_HARD('a list holds one reference or more; pose_student_align%new')
        self%b_ptr => build
        self%nmax  =  nmax
        call self%smap%new(params, build)
        allocate(self%vals(build%pftc%get_nrots()), self%corrs0(nmax), source=0.)
        allocate(self%irots0(nmax), self%refs(nmax), source=0)
        self%exists = .true.
    end subroutine new

    !> particle iptcl against the directions dirs and their partners mirr(dirs) (each reference once, in that order):
    !! the best alignment over every rotation and the disk; res%iref 0 when the list is empty or exceeds nmax
    subroutine align( self, iptcl, dirs, mirr, res )
        class(pose_student_align), intent(inout) :: self
        integer,                   intent(in)    :: iptcl, dirs(:), mirr(:)
        type(pose_align_result),   intent(out)   :: res
        real    :: corr_k, sh_k(2), sdev_k
        integer :: n, j, k, irot_k
        logical :: found, valid, l_over
        n      = 0
        l_over = .false.
        do j = 1, size(dirs)
            call add(dirs(j))
            call add(mirr(dirs(j)))
        enddo
        res%iref = 0
        if( n < 1 .or. l_over ) return
        call self%smap%prep(iptcl)
        do k = 1, n
            call self%b_ptr%pftc%gen_objfun_vals(self%refs(k), iptcl, [0.,0.], self%vals)
            self%irots0(k) = maxloc(self%vals, dim=1)
            self%corrs0(k) = self%vals(self%irots0(k))
        enddo
        call self%smap%srch(self%refs(1:n), self%corrs0(1:n), n)
        res%corr = -huge(res%corr)
        do k = 1, n
            call self%smap%get_result(self%refs(k), found, valid, irot_k, corr_k, sh_k, sdev_k)
            if( .not. (found .and. valid) )then
                ! the stored shift at its best rotation beat the map's pose
                irot_k = self%irots0(k)
                corr_k = self%corrs0(k)
                sh_k   = 0.
            endif
            if( corr_k > res%corr )then
                res%iref = self%refs(k)
                res%irot = irot_k
                res%corr = corr_k
                res%sh   = sh_k
                res%sdev = sdev_k
            endif
        enddo

    contains

        subroutine add( iref )
            integer, intent(in) :: iref
            if( any(self%refs(1:n) == iref) ) return
            if( n == self%nmax )then
                l_over = .true.
                return
            endif
            n = n + 1
            self%refs(n) = iref
        end subroutine add

    end subroutine align

    subroutine kill( self )
        class(pose_student_align), intent(inout) :: self
        if( .not. self%exists ) return
        call self%smap%kill
        deallocate(self%vals, self%corrs0, self%irots0, self%refs)
        nullify(self%b_ptr)
        self%nmax   = 0
        self%exists = .false.
    end subroutine kill

end module simple_pose_student_align
