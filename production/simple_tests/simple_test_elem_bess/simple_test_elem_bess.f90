program simple_test_elem_bess
use simple_winfuns, only: winfuns
use simple_kbinterpol ! use all in there
use simple_math,    only: euclid
implicit none
real, allocatable :: grid(:), vals_elem(:), vals(:)
real, parameter   :: Whalf=1.5, alpha=2.0, stepsz=0.001
type(winfuns)     :: wfs
integer           :: npoints, ipoint
real              :: point, val
wfs = winfuns('kb', Whalf, alpha)
call init_kbiterpol( Whalf, alpha )
point   = -Whalf
npoints = 0
do while( point < Whalf )
    npoints = npoints + 1
    point   = point + stepsz
end do
print *, 'NUMBER OF MEASUREMENTS: ', npoints
allocate(grid(npoints), vals(npoints), vals_elem(npoints))
point   = -Whalf
npoints = 0
do while( point < Whalf )
    npoints       = npoints + 1 
    grid(npoints) = point
    point         = point + stepsz
end do
vals_elem = kb_apod(grid)
do ipoint=1,npoints
    vals(ipoint) = wfs%eval_apod(grid(ipoint))
end do
print *, 'APOD DIFFERENCE: ', euclid(vals,vals_elem)
vals_elem = kb_instr(grid)
do ipoint=1,npoints
    vals(ipoint) = wfs%eval_instr(grid(ipoint))
end do
print *, 'INSTR DIFFERENCE: ', euclid(vals,vals_elem)

end program simple_test_elem_bess


