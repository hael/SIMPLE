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
    integer              :: ldim(3) = 0, i_rad_max = 0, maxpix = 0
    real,    allocatable :: cis(:), cjs(:), ring_vals(:,:)
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
                ! update number of pixels per ring
                self%npix_per_ring(i_rad) = self%npix_per_ring(i_rad) + 1
            end do
        end do
        self%maxpix = maxval(self%npix_per_ring)
        allocate(self%ring_vals(self%i_rad_max,self%maxpix), source=0.) ! heap array for threading?
        self%exists = .true.
    end subroutine new

    pure function get_rad_max( self ) result( i_rad_max )
        class(radial_medians), intent(in) :: self
        integer :: i_rad_max
        i_rad_max = self%i_rad_max
    end function get_rad_max

    subroutine calc_radial_medians( self, img, medians )
        class(radial_medians), intent(inout) :: self
        class(image),          intent(in)    :: img
        real,                  intent(inout) :: medians(self%i_rad_max)
        real(kind=c_float), pointer :: rmat(:,:,:)=>null() !< image pixels in img
        integer :: i, j, i_rad
        call img%get_rmat_ptr(rmat)
        self%npix_per_ring = 0
        do j = 1,self%ldim(2)
            do i = 1,self%ldim(1)
                ! find integer radius
                i_rad = nint(hyp(self%cis(i),self%cjs(j)))
                ! update number of pixels per ring
                self%npix_per_ring(i_rad) = self%npix_per_ring(i_rad) + 1
                ! extract pixel value
                self%ring_vals(i_rad,self%npix_per_ring(i_rad)) = rmat(i,j,1)
            end do
        end do
        do i_rad = 1,self%i_rad_max
            medians(i_rad) = median_nocopy(self%ring_vals(i_rad,:self%npix_per_ring(i_rad)))
        end do
    end subroutine calc_radial_medians

    subroutine kill( self )
        class(radial_medians), intent(inout) :: self
        if( self%exists )then
            self%ldim      = 0
            self%i_rad_max = 0
            self%maxpix    = 0
            deallocate(self%cis, self%cjs, self%ring_vals, self%npix_per_ring)
        endif
    end subroutine kill

end module simple_radial_medians
