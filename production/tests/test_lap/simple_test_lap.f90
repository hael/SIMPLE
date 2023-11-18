program simple_test_lap
use simple_lap, only: lap
implicit none  ! Ensures all variables must have an explicit type (all modern Fortran programs should include this statement)
integer, parameter :: N_PTCLS = 100, N_REFS = 3
type(lap) :: lap_obj
integer   :: sol(N_PTCLS), i
real      :: mat(N_REFS,N_PTCLS), total, cost_mat(N_PTCLS,N_PTCLS)
! mat(:,1) = [3, 6, 4]
! mat(:,2) = [3, 2, 7]
! mat(:,3) = [1, 5, 8]
! mat(:,4) = [4, 1, 5]
! mat(:,5) = [8, 3, 4]
! mat(:,6) = [2, 6, 1]
! cost_mat(1:3,:) = mat
! cost_mat(4:6,:) = mat
call random_number(cost_mat)
call lap_obj%new(cost_mat)
call lap_obj%solve_lap(sol)
total = 0.
do i = 1, N_PTCLS
    total = total + cost_mat(sol(i), i)
    ! print *, i, ' -> ', mod(sol(i)-1, N_REFS) + 1
enddo
print *, total
end program simple_test_lap