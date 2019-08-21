program simple_test_graphene_mask
use simple_math
use simple_defs
implicit none
integer, parameter   :: box  = 160
real,    parameter   :: smpd = 0.358
real,    allocatable :: res(:)
logical, allocatable :: graphene_mask(:)
integer              :: i
res           = get_resarr( BOX, SMPD )
graphene_mask = calc_graphene_mask( box, smpd )
write(*,*) 'RES OF GRAPHENE BAND 1: ', GRAPHENE_BAND1
write(*,*) 'RES OF GRAPHENE BAND 2: ', GRAPHENE_BAND2
do i=1,size(res)
    write(*,*) '>>> RESOLUTION:', res(i), graphene_mask(i)
end do
end program simple_test_graphene_mask
