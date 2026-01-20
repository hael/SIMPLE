program simple_test_msk_routines
! testing msk routines 
use simple_core_module_api
use simple_image
type(image), allocatable    :: stk(:), stk_copy(:)
type(image) :: img
integer :: nimgs  = 10, i, ldim(3)
real    :: smpd 

ldim = [1024,1024,1]
smpd = 3.0

! testing serial 
call img%new(ldim, smpd)
call img%ran()
call img%mask(ldim(1)/2., 'soft')
call img%ran()
call img%mask(ldim(1)/2., 'softavg')
call img%ran()
call img%mask(ldim(1)/2., 'hard')

allocate(stk(nimgs))
allocate(stk_copy(nimgs))
do i = 1, nimgs
    call stk(i)%new(ldim, smpd)
    call stk(i)%ran()
end do 
! testing parallel 

!$omp parallel do default(shared)  private(i) proc_bind(close) schedule(static)
do i = 1, nimgs
    call stk(i)%mask(ldim(1)/2., 'soft')
end do
!$omp end parallel do

do i = 1, nimgs
    call stk(i)%ran()
end do 
!$omp parallel do default(shared)  private(i) proc_bind(close) schedule(static)
do i = 1, nimgs
    call stk(i)%mask(ldim(1)/2., 'softavg')
end do
!$omp end parallel do

do i = 1, nimgs
    call stk(i)%ran()
end do 

!$omp parallel do default(shared)  private(i) proc_bind(close) schedule(static)
do i = 1, nimgs
    call stk(i)%mask(ldim(1)/2., 'hard')
end do
!$omp end parallel do

end program simple_test_msk_routines 





