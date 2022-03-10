program simple_test_openmp
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
implicit none
#include "simple_local_flags.inc"

integer :: vals(100), i,j, counts(10), correct_counts(10)

! reference, no openmp
vals    = 0
counts  = 0
do i=1,100
    vals(i) = nint(real(i)/real(10)) + 1
enddo
do j=1,10
	do i = 1,100
		if(vals(i)==j)then
			counts(j) = counts(j)+1
		endif
	enddo
enddo
correct_counts = counts

! safe
vals    = 0
counts  = 0
!$omp parallel default(shared) proc_bind(close) private(i,j)
!$omp do schedule(static)
do i=1,100
    vals(i) = nint(real(i)/real(10)) + 1
enddo
!$omp end do
!$omp do schedule(static)
do j=1,10
	do i = 1,100
		if(vals(i)==j)then
			counts(j) = counts(j)+1
		endif
	enddo
enddo
!$omp end do
!$omp end parallel

if(all(counts==correct_counts))then
	print *,'passed scenario one'
else
	print *,'failed scenario one'
endif

! unsafe, nowait
vals    = 0
counts  = 0
!$omp parallel default(shared) proc_bind(close) private(i,j)
!$omp do schedule(static)
do i=1,100
    vals(i) = nint(real(i)/real(10)) + 1
enddo
!$omp end do nowait
!$omp do schedule(static)
do j=1,10
	do i = 1,100
		if(vals(i)==j)then
			counts(j) = counts(j)+1
		endif
	enddo
enddo
!$omp end do
!$omp end parallel

if(all(counts==correct_counts))then
	print *,'passed scenario two'
else
	print *,'failed scenario two'
endif
end program simple_test_openmp
