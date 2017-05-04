program simple_test_ptcl_conv
use simple_defs          ! use all in there
use simple_math          ! use all in there
use simple_filehandling, only: read_filetable, nlines
use simple_oris,         only: oris
use simple_jiffys,       only: progress
implicit none
character(len=STDLEN), parameter   :: DOCLIST='mylistofdocs.txt'
character(len=STDLEN), allocatable :: oritabs(:)
type(oris),            allocatable :: osarr(:)
integer,               allocatable :: classes(:,:), conv_points(:), predicted_conv_points(:)
real,                  allocatable :: tmp(:), mi_joint(:)
integer :: ndocs, noris, idoc, iori, final_class
real    :: med_conv
! read & allocate
call read_filetable(DOCLIST, oritabs)
ndocs = size(oritabs)
noris = nlines(oritabs(1))
allocate( osarr(ndocs), classes(ndocs,noris), conv_points(noris), predicted_conv_points(noris))
print *, 'extracting ori and class info'
do idoc=1,ndocs
    call progress(idoc, ndocs)
    call osarr(idoc)%new(noris)
    call osarr(idoc)%read(oritabs(idoc))
    tmp = osarr(idoc)%get_all('class')
    classes(idoc,:) = nint(tmp)
    deallocate(tmp)
end do
print *, 'generating convergence stats'
! looking at one orientation over all docs (iterations)
do iori=1,noris
    call progress(iori, noris)
    ! walk backwards to identify the point of convergence
    final_class = classes(ndocs-1,iori)
    do idoc=ndocs-1,1,-1
        if( classes(idoc,iori) == final_class )then
            cycle
        else
            conv_points(iori) = idoc - 1
            exit
        endif
    end do
end do
print *, 'calculating iteration median'
med_conv = median(real(conv_points))
print *, 'median point of convergence: ', nint(med_conv)
print *, 'attempt to predict covergence points per particle'
mi_joint = osarr(1)%get_all('mi_joint')
predicted_conv_points = 0
do idoc=2,ndocs
    tmp = osarr(idoc)%get_all('mi_joint')
    mi_joint = 0.5 * mi_joint + 0.5 * tmp
    do iori=1,noris
        if( predicted_conv_points(iori) == 0 .and. mi_joint(iori) >= 0.95 )then
            predicted_conv_points(iori) = idoc
        endif
    end do
end do
where( predicted_conv_points == 0 ) predicted_conv_points = ndocs
print *, 'predicted_conv - conv = ', real(sum(predicted_conv_points - conv_points))/real(noris)
med_conv = median(real(predicted_conv_points))
print *, 'median point of predicted convergence: ', nint(med_conv)
end program simple_test_ptcl_conv
