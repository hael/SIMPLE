module simple_radial_medians
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image, only: image
implicit none

public :: radial_medians
private
#include "simple_local_flags.inc"

type radial_medians
    private
    integer              :: ldim(3) = 0, i_rad_max = 0, max_npix_per_ring = 0
    real,    allocatable :: cis(:), cjs(:)
    integer, allocatable :: npix_per_ring(:)
    logical              :: exists = .false.
  contains
    procedure :: new
    procedure :: get_rad_max
    procedure :: calc_radial_medians
    procedure :: kill
end type radial_medians

contains

    subroutine new( self, ldim )
        class(radial_medians), intent(inout) :: self
        integer,               intent(in)    :: ldim(3)
        integer :: i, j, i_rad
        if( ldim(3) /= 1 ) THROW_HARD('only for 2D images')
        call self%kill
        self%ldim = ldim
        ! init center as origin
        allocate(self%cis(self%ldim(1)), self%cjs(self%ldim(2)), source=0.)
        forall(i=1:self%ldim(1)) self%cis(i) = -real(self%ldim(1))/2. + real(i-1)
        forall(i=1:self%ldim(2)) self%cjs(i) = -real(self%ldim(2))/2. + real(i-1)
        ! find maximum integer radius
        self%i_rad_max = 0
        do j = 1,self%ldim(2)
            do i = 1,self%ldim(1)
                ! find integer radius
                i_rad = nint(hyp(self%cis(i),self%cjs(j)))
                ! update maximum integer radius
                if( i_rad > self%i_rad_max ) self%i_rad_max = i_rad
            end do
        end do
        ! count number of pixels per ring
        allocate(self%npix_per_ring(self%i_rad_max), source=0)
        do j = 1,self%ldim(2)
            do i = 1,self%ldim(1)
                ! find integer radius
                i_rad = nint(hyp(self%cis(i),self%cjs(j)))
                if( i_rad < 1 ) cycle
                ! update number of pixels per ring
                self%npix_per_ring(i_rad) = self%npix_per_ring(i_rad) + 1
            end do
        end do
        self%max_npix_per_ring = maxval(self%npix_per_ring)
        self%exists = .true.
    end subroutine new

    pure function get_rad_max( self ) result( i_rad_max )
        class(radial_medians), intent(in) :: self
        integer :: i_rad_max
        i_rad_max = self%i_rad_max
    end function get_rad_max

    subroutine calc_radial_medians( self, img, stats, medians )
        class(radial_medians), intent(inout) :: self
        class(image),          intent(in)    :: img
        type(stats_struct),    intent(inout) :: stats
        real,                  intent(inout) :: medians(self%i_rad_max)
        real(kind=c_float), pointer :: rmat(:,:,:)=>null() !< image pixels in img
        integer :: npix_per_ring(self%i_rad_max), i, j, i_rad
        real    :: ring_vals(self%i_rad_max,self%max_npix_per_ring), var
        logical :: mask(self%i_rad_max,self%max_npix_per_ring), err
        call img%get_rmat_ptr(rmat)
        npix_per_ring = 0
        mask          = .false.
        do j = 1,self%ldim(2)
            do i = 1,self%ldim(1)
                ! find integer radius
                i_rad = nint(hyp(self%cis(i),self%cjs(j)))
                if( i_rad < 1 ) cycle
                ! update number of pixels per ring
                npix_per_ring(i_rad) = npix_per_ring(i_rad) + 1
                ! extract pixel value
                ring_vals(i_rad,npix_per_ring(i_rad)) = rmat(i,j,1)
                ! update mask
                mask(i_rad,npix_per_ring(i_rad)) = .true.
            end do
        end do
        ! overall stats
        call avg_sdev(ring_vals, stats%avg, stats%sdev, mask)
        ! ring medians
        do i_rad = 1,self%i_rad_max
            medians(i_rad) = median_nocopy(ring_vals(i_rad,:npix_per_ring(i_rad)))
        end do
    end subroutine calc_radial_medians

    subroutine kill( self )
        class(radial_medians), intent(inout) :: self
        if( self%exists )then
            self%ldim              = 0
            self%i_rad_max         = 0
            self%max_npix_per_ring = 0
            deallocate(self%cis, self%cjs, self%npix_per_ring)
        endif
    end subroutine kill

end module simple_radial_medians
