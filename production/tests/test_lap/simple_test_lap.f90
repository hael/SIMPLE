program simple_test_lap
use simple_lap, only: lap
implicit none  ! Ensures all variables must have an explicit type (all modern Fortran programs should include this statement)
type(lap) :: lap_obj
integer   :: sol(3)
call lap_obj%new
call lap_obj%solve_lap(sol)
print *, sol
end program simple_test_lap