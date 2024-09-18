program simple_test_kpca
include 'simple_lib.f08'
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
use simple_kpca_svd,   only: kpca_svd
implicit none
integer, parameter :: D = 2, Q = 4
type(kpca_svd)     :: kpca_obj
type(cmdline)      :: cline
type(parameters)   :: p
real, allocatable  :: XY(:,:), data_here(:,:)
integer            :: N, funit, i
real               :: avg(D)
character(len=120) :: binfname
call cline%parse_oldschool
call p%new(cline)
binfname = '/home/vanc2/Downloads/output.txt'
open(unit=funit,file=binfname)
read(funit,*) N
allocate(XY(D,N))
read(funit,*) XY
call fclose(funit)
call kpca_obj%new(N, D, Q)
call kpca_obj%master(XY)
allocate( data_here(D,N), source=0.)
avg = 0.
do i = 1, N
    call kpca_obj%generate(i, avg, data_here(:,i))
enddo
binfname = '/home/vanc2/result.txt'
open(unit=funit,file=binfname)
write(funit,*) transpose(data_here)
call fclose(funit)
end program simple_test_kpca
