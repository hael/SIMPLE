!@descr: immutable slab of decoded real-space images, the storage of the particle buffer
! One contiguous allocation of n images of box x box reals, filled slot by slot, sealed, then read by copy:
! an image read out of a slot equals, bit for bit and with its padding zeroed, the image that was put in.
! The slab is allocated only within a byte budget and never stops the run: a refused or failed allocation
! leaves it empty for the caller's fallback. Policy: section 3.1 of particle_io_few_particle_stacks.md.
module simple_image_slab
use, intrinsic :: iso_c_binding,   only: c_float
use, intrinsic :: iso_fortran_env, only: int64
use simple_core_module_api
use simple_image, only: image
implicit none

public :: image_slab, SLAB_OK, SLAB_OVER_BUDGET, SLAB_NO_MEMORY
private
#include "simple_local_flags.inc"

integer, parameter :: SLAB_OK = 0, SLAB_OVER_BUDGET = 1, SLAB_NO_MEMORY = 2

type :: image_slab
    private
    real(c_float), allocatable :: pix(:,:,:)   !< box x box x n decoded pixels
    integer :: box = 0, n = 0
    logical :: l_sealed = .false.
contains
    procedure :: new
    procedure :: put
    procedure :: seal
    procedure :: get
    procedure :: bytes
    procedure :: exists
    procedure :: kill
end type image_slab

contains

    !> n slots of box x box: status SLAB_OK, or the slab stays empty with SLAB_OVER_BUDGET (more than budget
    !! bytes) or SLAB_NO_MEMORY (the allocation failed)
    subroutine new( self, box, n, budget, status )
        class(image_slab), intent(inout) :: self
        integer,           intent(in)    :: box, n
        integer(int64),    intent(in)    :: budget
        integer,           intent(out)   :: status
        integer :: ierr
        call self%kill
        if( int(box,int64)**2 * int(n,int64) * 4_int64 > budget )then
            status = SLAB_OVER_BUDGET
            return
        endif
        allocate(self%pix(box,box,n), stat=ierr)
        if( ierr /= 0 )then
            status = SLAB_NO_MEMORY
            return
        endif
        self%box = box
        self%n   = n
        status   = SLAB_OK
    end subroutine new

    !> slot k receives the pixels of img, a real-space image of the slab's box; only before the seal
    subroutine put( self, k, img )
        class(image_slab), intent(inout) :: self
        integer,           intent(in)    :: k
        class(image),      intent(inout) :: img
        real(c_float), pointer :: rmat(:,:,:)
        integer :: ldim(3)
        if( self%l_sealed ) THROW_HARD('the slab is sealed; put')
        if( k < 1 .or. k > self%n ) THROW_HARD('slot out of range; put')
        ldim = img%get_ldim()
        if( any(ldim /= [self%box,self%box,1]) .or. img%is_ft() ) THROW_HARD('not a real-space image of the slab box; put')
        call img%get_rmat_ptr(rmat)
        self%pix(:,:,k) = rmat(1:self%box,1:self%box,1)
    end subroutine put

    !> no slot changes after this
    subroutine seal( self )
        class(image_slab), intent(inout) :: self
        self%l_sealed = .true.
    end subroutine seal

    !> img receives slot k: its pixels, padding zeroed, flagged real space
    subroutine get( self, k, img )
        class(image_slab), intent(in)    :: self
        integer,           intent(in)    :: k
        class(image),      intent(inout) :: img
        real(c_float), pointer :: rmat(:,:,:)
        integer :: ldim(3)
        if( .not. self%l_sealed ) THROW_HARD('the slab is not sealed; get')
        if( k < 1 .or. k > self%n ) THROW_HARD('slot out of range; get')
        ldim = img%get_ldim()
        if( any(ldim /= [self%box,self%box,1]) ) THROW_HARD('image does not match the slab box; get')
        call img%get_rmat_ptr_padded(rmat)
        rmat(1:self%box,1:self%box,1) = self%pix(:,:,k)
        rmat(self%box+1:,:,1) = 0.
        call img%set_ft(.false.)
    end subroutine get

    integer(int64) function bytes( self )
        class(image_slab), intent(in) :: self
        bytes = int(self%box,int64)**2 * int(self%n,int64) * 4_int64
    end function bytes

    logical function exists( self )
        class(image_slab), intent(in) :: self
        exists = allocated(self%pix)
    end function exists

    subroutine kill( self )
        class(image_slab), intent(inout) :: self
        if( allocated(self%pix) ) deallocate(self%pix)
        self%box      = 0
        self%n        = 0
        self%l_sealed = .false.
    end subroutine kill

end module simple_image_slab
