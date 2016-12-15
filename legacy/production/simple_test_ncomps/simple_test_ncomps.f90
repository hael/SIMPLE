program simple_test_ncomps
use simple_params,  only: params
use simple_image,   only: image
use simple_cmdline, only: cmdline
use simple_defs
use simple_timing
implicit none
type(image)   :: img
type(params)  :: p
type(cmdline) :: cline
integer       :: find, nradial_lines, cnt, h, k, lims(3,2), i, j
if( command_argument_count() < 3 )then
    write(*,'(a)', advance='no') 'SIMPLE_TEST_NCOMPS box=<box size (in pixels)> smpd=<sampling distance(in A)>'
    write(*,'(a)') ' lp=<low-pass limit(in A)>'
    stop
endif
call cline%parse
call cline%checkvar('box',  1)
call cline%checkvar('smpd', 2)
call cline%checkvar('lp',   3)
call cline%check
call cline%set('prg', 'test_ncomps')
p = params(cline)
find = int((real(p%box-1)*p%smpd)/p%lp)
nradial_lines = int(twopi*(real(p%box)/2.))
cnt = 0
print *, 'nradial:', nradial_lines
do i=1,nradial_lines/2
    do j=2,find
        cnt = cnt+1
    end do 
end do
print *, 'nr of comps in polar representation:', cnt
call img%new([p%box,p%box,1], p%smpd)
cnt = 0
lims = img%loop_lims(2,p%lp)
print *, 'lims1:', lims(1,1), lims(1,2)
print *, 'lims2:', lims(2,1), lims(2,2)
do h=lims(1,1),lims(1,2)
    do k=lims(2,1),lims(2,2)
        cnt = cnt+1
    end do
end do
print *, 'nr of comps in cartesian representation:', cnt
end program
