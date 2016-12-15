!     FortranCode.f

integer function sumsquaredf(n,a)
  integer a(*)
  !write(*,'(a)')"-- We are now in the FORTRAN program FortranCode --"
  !write(*,'(a)')"Print contents of array A() copied from a[] in C "

  do i=1,n
     !write(*,'(2i5)')i,a(i)
  enddo
  !write(*,*)!'("calculate then print out squares of elements with accumulated sums ")')
  isum=0
  do i=1,n
     a(i)=exp( 1.d0*a(i) )
     isum=isum+a(i)
     !write(*,'(3i5)')i,a(i),isum
  enddo
  sumsquaredf=isum
end function sumsquaredf
