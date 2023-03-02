! Tests both implementations of gauss2Dfit, which fits each particle
! in a micrograph with a 2D multivariate Gaussian distribution
program simple_test_gauss2Dfit

include 'simple_lib.f08'
use simple_picker_utils, only: picker_utils
use simple_image,        only: image
use simple_gauss2Dfit,   only: gauss2Dfit
implicit none

type(image),    allocatable :: refs(:), fits(:)
real,             parameter :: SMPD=2.0, TOLERANCE=0.01
real,          allocatable  :: corrs_sort(:), corrs(:), centers_sort(:,:)
real,          allocatable  :: centers(:,:), sigmas(:,:,:), sigmas_sort(:,:,:)
real                        :: sdevs(2), eigenvals(2), eigenvecs(2,2)
integer,          parameter :: NTESTS=2
integer,       allocatable  :: irefs(:)
integer                     :: ldim(3), ldim_refs(3), iref, nrefs, ifoo, nrot
integer                     :: i, j, funit, tpassed
integer(timer_int_kind)     :: t
real(timer_int_kind)        :: rt1, rt2
logical                     :: match
character(len=*), parameter :: micname = '20170629_00021_frameImage_intg_after_filters.mrcs'
character(len=*), parameter :: micrefs = 'refs.mrcs'
character(len=*), parameter :: micfits = 'fit.mrcs'
character(len=*), parameter :: out = 'out.csv'

character(len=*), parameter :: csv_head = 'INDEX'//CSV_DELIM//'CORR'//&
    &CSV_DELIM//'X'//CSV_DELIM//'Y'//CSV_DELIM//'SDEV_MAJ'//CSV_DELIM//'SDEV_MIN'

tpassed = 0
write(logfhandle,'(a)') 'Test 1: Testing an unallocated array of references'
call gauss2Dfit(refs, centers, sigmas, corrs)
if (.not. allocated(corrs)) then 
    tpassed = tpassed + 1
    write(logfhandle,'(a)') 'Test 1 passed'
else
    write(logfhandle,'(a)') 'Test 1 failed'
end if

call find_ldim_nptcls(micname, ldim, ifoo)
ldim_refs = [ldim(1), ldim(2), 1]
nrefs = ldim(3)
allocate(refs(nrefs), fits(nrefs))
allocate(irefs(nrefs), source=0)
allocate(centers(2,nrefs), corrs(nrefs), sigmas(2,2,nrefs), source=0.)
allocate(centers_sort(2,nrefs), corrs_sort(nrefs), sigmas_sort(2,2,nrefs),&
    &source=0.)

! Calls to both implementations of gauss2Dfit with timing tests
do iref = 1,nrefs
    call refs(iref)%new(ldim_refs, SMPD)
    call refs(iref)%read(micname, iref)
end do
t = tic()
call gauss2Dfit(refs, centers, sigmas, corrs)
rt1 = toc(t)
call gauss2Dfit(refs, centers_sort, sigmas_sort, corrs_sort, irefs, fits)
rt2 = toc(t) - rt1

! Here we assume that if correlations match then all other stats match
write(logfhandle,'(a)') 'Test 2: Testing that both implementations agree.'
match = .true.
do i=1, nrefs
    if (abs(corrs_sort(i) - corrs(irefs(i))) > TOLERANCE) then
        match = .false.
        write(logfhandle, *) "Test 2 Failure at sorted index: ", i
        write(logfhandle, *) "corrs: ", corrs_sort(i), corrs(irefs(i))
    end if
end do
if (match) then
    write(logfhandle,'(a)') 'Test 2 passed'
    tpassed = tpassed + 1
end if

write(logfhandle,'(a)') 'Printing out sorted stats for user sanity check.'
call fopen(funit, status='REPLACE', action='WRITE', file=trim(out))
write(funit, '(a)') csv_head
601 format(F8.3,A2)
602 format(F8.3)
do i=1, nrefs
    call jacobi(sigmas_sort(:2,:2,i), 2, 2, eigenvals, eigenvecs, nrot)
    call eigsrt(eigenvals, eigenvecs, 2, 2)
    sdevs(:2) = sqrt(eigenvals(:2))
    ! Write stats
    write(funit, 601, advance='no') real(irefs(i)),     CSV_DELIM   ! INDEX
    write(funit, 601, advance='no') corrs_sort(i),      CSV_DELIM   ! CORR
    write(funit, 601, advance='no') centers_sort(1,i),  CSV_DELIM   ! X
    write(funit, 601, advance='no') centers_sort(2,i),  CSV_DELIM   ! Y
    write(funit, 601, advance='no') sdevs(1),           CSV_DELIM   ! SDEV_MAJ
    write(funit, 602)               sdevs(2)                        ! SDEV_MIN
    ! Last particle in output micrographs has greatest correlation
    call refs(irefs(i))%write(micrefs, i)
    call fits(irefs(i))%write(micfits, i)
end do

call fclose(funit)
deallocate(refs, fits, corrs, corrs_sort, sigmas, sigmas_sort, centers, &
    &centers_sort, irefs)

! End gracefully
if (tpassed == 2) then
    write(logfhandle,'(a)') 'SIMPLE_TEST_GAUSS2DFIT COMPLETE W/ SUCCESS'
    write(logfhandle,'(a,i1,a,i1,a)') 'PASSED ', tpassed, '/', NTESTS, ' TESTS'
else
    write(logfhandle,'(a)') 'test_gauss2Dfit COMPLETE W/ FAILURE'
    write(logfhandle,'(a,i1,a,i1,a)') 'PASSED ', tpassed, '/', NTESTS, ' TESTS'
end if
write(logfhandle,'(a,f10.6,a)') 'gauss2Dfit_1 Running Time: ', rt1, 's'
write(logfhandle,'(a,f10.6,a)') 'gauss2Dfit_2 Running Time: ', rt2, 's'
write(logfhandle,'(a)') '**** END SIMPLE_TEST_GAUSS2DFIT ****'
end program simple_test_gauss2Dfit