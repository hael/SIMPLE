module simple_kbinterpol_mem
#include "simple_lib.f08"
implicit none

public :: kbinterpol_mem
private

type kbinterpol_mem
	private
	integer           :: nrots, wdim, lims(2,2)
	real, allocatable :: inpl_rots(:), w(:,:,:,:,:)
	logical           :: exists = .false.
contains
	procedure :: new
	procedure :: get_wdim
	procedure :: memoize_kb
	procedure :: fetch
	procedure :: kill
end type kbinterpol_mem

contains 

	subroutine new( self, inpl_rots, lims )
		class(kbinterpol_mem), intent(inout) :: self
		real,                  intent(in)    :: inpl_rots(:)
		real,                  intent(in)    :: lims(2,2)
		! set constants
		self%nrots = size(inpl_rots)
        self%wdim  = ceiling(KBALPHA * KBWINSZ) + 1
        self%lims  = lims
        ! allocate
		allocate( self%w(self%nrots,self%lims(1,1):self%lims(1,2),&
			&self%lims(2,1):self%lims(2,2),self%wdim,self%wdim), source=1.0 )
		allocate( self%inpl_rots(self%nrots), source=inpl_rots )
		! flag existence
		self%exists = .true.
	end subroutine new

	pure integer function get_wdim( self )
		class(kbinterpol_mem), intent(in) :: self
		get_wdim = self%wdim
	end function get_wdim

	subroutine memoize_kb( self, kbwin )
		use simple_kbinterpol, only: kbinterpol
		class(kbinterpol_mem), intent(inout) :: self
		class(kbinterpol),     intent(inout) :: kbwin
		integer :: irot, h, k, l, incr, win(2,2)
		real    :: mat(2,2), vec(2), loc(2)
		! Rotations loop
		!$omp parallel do default(shared) private(irot,mat,h,k,vec,loc,win,l,incr)&
		!$omp schedule(static) proc_bind(close)
		do irot=1,self%nrots
			mat = rotmat2d( -self%inpl_rots(irot) )
			! Fourier components loop
	        do h=self%lims(1,1),self%lims(1,2)
	            do k=self%lims(2,1),self%lims(2,2)
	                ! calculate non-uniform sampling location
	                vec   = [real(h), real(k)]
	                loc   = matmul(vec,mat)
	                ! window
	                win = sqwin_2d(loc(1),loc(2), KBWINSZ)
	                ! kernel values
	                self%w(irot,h,k,:,:) = 1.
	                do l=1,self%wdim
	                    incr = l - 1
	                    ! interpolation kernel matrix
	                    self%w(irot,h,k,l,:) = self%w(irot,h,k,l,:) * &
	                    	&kbwin%apod( real(win(1,1) + incr) - loc(1))
	                    self%w(irot,h,k,:,l) = self%w(irot,h,k,:,l) * &
	                    	&kbwin%apod( real(win(2,1) + incr) - loc(2))
	                end do
	            end do
	        end do
	    end do
	    !$omp end parallel do
	end subroutine memoize_kb

	function fetch( self, irot, h, k ) result( w )
		class(kbinterpol_mem), intent(in) :: self
		integer,               intent(in) :: irot, h, k
		real :: w(self%wdim,self%wdim)
		w = self%w(irot,h,k,:,:)
	end function fetch

	subroutine kill( self )
		class(kbinterpol_mem), intent(inout) :: self
		if( self%exists )then
			deallocate(self%w, self%inpl_rots)
			self%exists = .false.
		endif
	end subroutine kill

	subroutine test_kbinterpol_mem
		use simple_ftiter, only: ftiter
		use simple_timer   ! use all in there
		type(kbinterpol_mem) :: kbmem
		type(ftiter)         :: fit
		integer, parameter   :: NROTS=359, NTST=100
		real    :: lims(3,2)
		real    :: inpl_rots(NROTS)
		integer :: jrot, itst
		integer(timer_int_kind) :: t_init, t_prep_pftcc
		real(timer_int_kind)    :: rt_init, rt_prep_pftcc
		do jrot=1,NROTS
			inpl_rots(jrot) = real(jrot)
		end do
		call fit%new([240,240,1], 1.06)
		lims = fit%loop_lims(3)
		call kbmem%new(inpl_rots, lims(:2,:))
		do itst=1,NTST

		end do



	end subroutine test_kbinterpol_mem

end module simple_kbinterpol_mem