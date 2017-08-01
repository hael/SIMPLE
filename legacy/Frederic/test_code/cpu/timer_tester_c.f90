! ============================================================================
! Name        : timer_tester.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 6th of August 2013
! Copyright   : Your copyright notice
! Description : tests the fnuctionality of the timing routine in C using a 
!             : matrix multiplication of nxm mxn dense matrix.
! ============================================================================
program timer_tester_c

  use simple_SU3_tester

  implicit none

  !global variables
  double precision,dimension(3)   :: elps_tm
  double precision,dimension(3)   :: t1s,t2s

  !local variables
  !timing variables
  double precision                :: elapse_time_mm
  double precision                :: t1_mm,t2_mm

  double precision                :: elapse_time_u1, elapse_time_u2
  double precision                :: t1_u1,t2_u1
  double precision                :: t1_u2,t2_u2

  !start of the excution commands

  call su3randomlinks_mul_timing(elps_tm,t1s,t2s)

  elapse_time_u1 = elps_tm(1)
  elapse_time_u2 = elps_tm(2)
  elapse_time_mm = elps_tm(3)

  t1_u1 = t1s(1)
  t1_u2 = t1s(2)
  t1_mm = t1s(3)

  t2_u1 = t2s(1)
  t2_u2 = t2s(2)
  t2_mm = t2s(3)

  !printing the results on screen

  write(*,'(a)')' '
  write(*,'(a)')'========================================================================='
  write(*,'(4(a,i2))')'The lattice size: nx= ',nx,', ny= ',ny,', nz= ',nz,', nt= ',nt 
  write(*,'(a)')' '
  write(*,'(a)')'The generation of the ---1st--- SU(3) matrix using su3random'
  write(*,'(a,f16.8,a)')'Using timming_c.cpp method, elapse_time: ',elapse_time_u1,' Secs'
  write(*,'(a,f16.8,a)')'Using cpu_time() method: ',t2_u1 - t1_u1,' Secs'
  write(*,'(a)')'The generation of the ---2nd--- SU(3) matrix using su3random'
  write(*,'(a,f16.8,a)')'Using timming_c.cpp method, elapse_time: ',elapse_time_u2,' Secs'
  write(*,'(a,f16.8,a)')'Using cpu_time() method: ',t2_u2 - t1_u2,' Secs'
  write(*,'(a)')' '
  write(*,'(a)')'The generation of the ---total time--- SU(3) matrix using su3random'
  write(*,'(a,f16.8,a)')'Using timming_c.cpp method, elapse_time: ',elapse_time_u2+elapse_time_u1,' Secs'
  write(*,'(a,f16.8,a)')'Using cpu_time() method: ',(t2_u2 - t1_u2)+(t2_u1 - t1_u1),' Secs'
  write(*,'(a)')' '
  write(*,'(a)')'Multiplication of two SU(3) matrix on random links (U(x)_1 x U(x)_2)'
  write(*,'(a,f16.8,a)')'Using timming_c.cpp method, elapse_time: ',elapse_time_mm,' Secs'
  write(*,'(a,f16.8,a)')'Using cpu_time() method: ',t2_mm - t1_mm,' Secs'
  write(*,'(a)')'========================================================================='
  write(*,'(a)')' '

end program timer_tester_c
