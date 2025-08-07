program simple_test_euler
include 'simple_lib.f08'
implicit none
type(oris) :: oris_obj
type(sym)  :: pgrpsyms               !< symmetry elements object
type(ori)  :: o1, o2
real       :: R(3,3),eullims(3,2), euls(3), euls2(3), diff, error, threshold
integer    :: isample, n
integer, parameter :: NSAMPLE = 100
call seed_rnd
call pgrpsyms%new('c4')
oris_obj = oris(1, is_ptcl=.false.)
call oris_obj%set_euler(1, [45., 90., 180.])
call oris_obj%get_ori(1, o1)
eullims = pgrpsyms%get_eullims()
print *, 'lim1: ', eullims(1,1), eullims(1,2)
print *, 'lim2: ', eullims(2,1), eullims(2,2)
print *, 'lim3: ', eullims(3,1), eullims(3,2)
print *, 'RND EULS'
call o2%new(.true.)
do isample = 1,NSAMPLE
    call pgrpsyms%rnd_euler(o2)
    euls = o2%get_euler()
    print *, euls(1), euls(2), euls(3)
end do
! consistency between spider and new m2euler routines
n = 0
error = 0.
threshold = 0.01 ! degrees
print *,'--- random'
do isample = 1,NSAMPLE
    call pgrpsyms%rnd_euler(o2)
    euls  = o2%get_euler()
    R     = o2%get_mat()
    euls  = m2euler(R); euls2 = m2euler_fast(R)
    diff  = angular_error(euls, euls2)
    error = error + diff
    if( diff > threshold )then
        n = n+1
        print *,n,isample,diff,euls,euls2
    endif
enddo
print *,'--- phi=0'
do isample = 1,NSAMPLE
    call pgrpsyms%rnd_euler(o2)
    euls  = o2%get_euler()
    euls(1) = 0.
    call o2%set_euler(euls)
    R     = o2%get_mat()
    euls  = m2euler(R); euls2 = m2euler_fast(R)
    diff  = angular_error(euls, euls2)
    error = error + diff
    if( diff > threshold )then
        n = n+1
        print *,n,isample,diff,euls,euls2
    endif
enddo
print *,'--- phi=90'
do isample = 1,NSAMPLE
    call pgrpsyms%rnd_euler(o2)
    euls  = o2%get_euler()
    euls(1) = 90.
    call o2%set_euler(euls)
    R     = o2%get_mat()
    euls  = m2euler(R); euls2 = m2euler_fast(R)
    diff  = angular_error(euls, euls2)
    error = error + diff
    if( diff > threshold )then
        n = n+1
        print *,n,isample,diff,euls,euls2
    endif
enddo
print *,'--- psi=0'
do isample = 1,NSAMPLE
    call pgrpsyms%rnd_euler(o2)
    euls  = o2%get_euler()
    euls(3) = 0.
    call o2%set_euler(euls)
    R     = o2%get_mat()
    euls  = m2euler(R); euls2 = m2euler_fast(R)
    diff  = angular_error(euls, euls2)
    error = error + diff
    if( diff > threshold )then
        n = n+1
        print *,n,isample,diff,euls,euls2
    endif
enddo
print *,'--- psi=270'
do isample = 1,NSAMPLE
    call pgrpsyms%rnd_euler(o2)
    euls  = o2%get_euler()
    euls(3) = 270.
    call o2%set_euler(euls)
    R     = o2%get_mat()
    euls  = m2euler(R); euls2 = m2euler_fast(R)
    diff  = angular_error(euls, euls2)
    error = error + diff
    if( diff > threshold )then
        n = n+1
        print *,n,isample,diff,euls,euls2
    endif
enddo
print *,'--- theta=0'
do isample = 1,NSAMPLE
    call pgrpsyms%rnd_euler(o2)
    euls  = o2%get_euler()
    euls(2) = 0.
    call o2%set_euler(euls)
    R     = o2%get_mat()
    euls  = m2euler(R); euls2 = m2euler_fast(R)
    diff  = min(angular_error(euls, euls2), angular_error(euls, [euls2(3),euls2(2),euls2(1)]))
    error = error + diff
    if( diff > threshold )then
        n = n+1
        print *,n,isample,diff,euls,euls2
    endif
enddo
if( n > 0 )then
    write(*,'(A,I6,A)')'M2EULER TEST NOT PASSED: ',n, ' ERRORS'
else
    write(*,*)'M2EULER TEST PASSED'
endif
write(*,'(A,F9.6,A)')'AVERAGE DIFFERENCE: ',error/real(6*NSAMPLE),' degrees'
contains

    real function angular_error( angs1, angs2 )
        real, intent(in) :: angs1(3), angs2(3)
        angular_error = min((angs1(1)-angs2(1))**2, (360.-angs1(1)-angs2(1))**2)
        angular_error = angular_error + (angs1(2)-angs2(2))**2
        angular_error = angular_error + min((angs1(3)-angs2(3))**2, (360.-angs1(3)-angs2(3))**2)
        angular_error = sqrt(angular_error/3.)
    end function angular_error

end program simple_test_euler
