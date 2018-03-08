!------------------------------------------------------------------------------!
! SIMPLE v2.5         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
program test_binomial_para
include 'simple_lib.f08'

implicit none

integer,               parameter   :: NPARTS=10 
character(len=STDLEN), allocatable :: labels(:)
type(chash),           allocatable :: work_assignments(:)
integer :: i, nlayers, njobs, nworkers

! figure out how many workers are needed
if( mod(NPARTS,2) == 0. )then
    nworkers = NPARTS/2
else
    nworkers = (NPARTS-1)/2 + 1
endif
allocate(work_assignments(nworkers))
! prepare the labels
call prep_labels(NPARTS)
! first pass to count the number of layers & jobs
nlayers = 0
njobs   = 0
call fibonacci_partition(labels, nlayers, njobs)
print *, '# workers: ', nworkers
print *, '# layers:  ', nlayers
print *, '# jobs:    ', njobs

do i=1,nworkers
    ! # layers is max # assignments
    call work_assignments(i)%new(nlayers)
end do

! re-prep labels
deallocate(labels)
call prep_labels(NPARTS)

contains

    subroutine prep_labels( nparts )
        integer, intent(in) :: nparts
        integer :: i
        allocate(labels(nparts))
        do i=1,nparts
            labels(i) = int2str(i)
        end do
    end subroutine prep_labels

    recursive subroutine fibonacci_partition( labels, nlayers, njobs )
        character(len=STDLEN), allocatable :: labels(:)
        character(len=STDLEN), allocatable :: reduced_labels(:)
        integer, intent(inout)             :: nlayers, njobs
        integer :: n, nred, cnt
        n = size(labels)
        if( n == 1 ) return
        if( mod(n,2) == 0. )then
            nred = n/2
        else
            nred = (n-1)/2 + 1
        endif
        allocate(reduced_labels(nred))
        njobs = njobs + nred
        cnt   = 0
        do i=1,n,2
            cnt = cnt+1
            if( i+1 > n )then
                reduced_labels(cnt) = trim(labels(i))
            else
                reduced_labels(cnt) = trim(labels(i))//trim(labels(i+1))
            endif
            write(*,'(a,a)',advance='no') trim(reduced_labels(cnt)), ' '
        end do
        write(*,*) ''
        ! count the number of layers
        nlayers = nlayers + 1
        ! labels now become reduced labels
        deallocate(labels)
        allocate(labels(nred), source=reduced_labels)
        deallocate(reduced_labels)
        ! recursion
        call fibonacci_partition(labels, nlayers, njobs)
    end subroutine fibonacci_partition

end program test_binomial_para
