program test_indexing
integer, parameter :: nx=3, ny=3, nxr=2, nyr=4
integer :: i,j,i_new,j_new
integer :: mat(3,3),cnt=0,matt(3,3)
integer :: rectmat(2,4), rectmatt(4,2)

cnt=0
do i=1,nx
    do j=1,ny
        cnt = cnt+1
        mat(i,j) = cnt
!         print *, i, j, mat(i,j)
    end do
end do
! print *, '******************'
do i=1,nx
    do j=1,ny
!         matt(i,j) = mat(nx-(i-1),j)
        matt(i,j) = mat(i,ny-(j-1))
!         print *, i, j, matt(i,j)
!         matt(i,j) = mat(ny-(j-1),i)
!         matt(i,j) = mat(nx-(i-1),ny-(j-1))
!         matt(i,j) = mat(ny-(j-1),nx-(i-1))
!         print *, i, j, matt(i,j)
    end do
end do

cnt=0
do i=1,nxr
    do j=1,nyr
        cnt = cnt+1
        rectmat(i,j) = cnt
        print *, i, j, rectmat(i,j)
    end do
end do
print *, '******************'
do i=1,nxr
    do j=1,nyr
        rectmatt(j,i) = rectmat(i,nyr-(j-1))
        print *, j, i, rectmatt(j,i)
    end do
end do




stop
print *, '******************'
print *, mat(1,:)
print *, mat(2,:)  
print *, mat(3,:)
print *, '******** modified **********'
! matt = transpose(mat)
print *, matt(1,:)
print *, matt(2,:)  
print *, matt(3,:)
print *, '******************'
stop
do i=1,nx
    i_new = nx-(i-1)
    do j=1,ny
        j_new = i
        print *, i, j, 'maps to: ', i_new, j_new
    end do
end do
end program
