program simple_test_euler_consensus
use simple_orialgn_ensemble, only: orialgn_ensemble
use simple_orialgn_pair,     only: orialgn_pair
use simple_params,           only: params
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_ran_tabu,         only: ran_tabu
use simple_cmdline           ! singleton
use simple_jiffys            ! singleton
use simple_defs              ! singleton
use simple_rnd               ! singleton
use simple_defs              ! singleton
use simple_math              ! singleton
implicit none
type(params)            :: p
type(orialgn_ensemble)  :: dockem
type(ori)               :: o
integer                 :: i, j, which, irnd
real                    :: angdist
type(oris)              :: consensus, facit
type(oris), allocatable :: members(:)
type(orialgn_pair)      :: op
integer, allocatable    :: pinds(:)
type(ran_tabu)          :: rt
character(len=STDLEN)   :: dig
if( command_argument_count() < 1 )then
    write(*,'(a)') 'simple_test_euler_consensus nmembers=<nr ensemble members> noris=<nr oris>'
    stop
endif
call parse_cmdline
call cmdcheckvar('nmembers', 1)
call cmdcheckvar('noris',    2)
call cmdcheck
p = params()
write(dig,*) p%part
call seed_rnd
allocate(members(p%nmembers), pinds(p%noris))
call dockem%new(p%nmembers)
!goto 99
! make the members for part 1
do i=1,p%nmembers/2
    call members(i)%new(p%noris)
    call members(i)%rnd_oris(p%eullims, p%trs)
end do
do i=p%nmembers/2+1,p%nmembers
    call members(i)%new(p%noris)
    call members(i)%spiral
    call members(i)%introd_alig_err(20., 0.)
end do
! set the orientations
call dockem%construct_ensemble(p%noris)
do i=1,p%nmembers
    do j=1,p%noris
        o = members(i)%get_ori(j)
        call dockem%set_ori(i,j,o)
    end do
end do
! generate the L1-median
consensus = dockem%l1median()
! find vote and report success/failure
call dockem%find_closest(consensus, angdist, which)
if( which <= 5 )then
    write(*,'(A)') '>>> VOTE SUCCESS: .NO.'
else
    write(*,'(A)') '>>> VOTE SUCCESS: .YES.'
endif
! THIS WILL IDENTIFY BAD ALIGNMENTS, THE PER PARTICLE MEASURE COULD ELIMINATE BAD PARTICLES

! IDEA: CLUSTER PARTILCES INTO TWO GROUPS AND SEE IF WE CAN REMOVE THE (OUTLIERS)

do i=1,p%nmembers
    write(*,*) 'DISTANCE 2 MEMBER:', i, 'IS:', rad2deg(dockem%dist2member(consensus, i)/real(p%noris))
end do
stop

! write consensus_1
call consensus%write('consensus_part'//trim(adjustl(dig))//'_1.txt')
! make the members for part 2
call facit%new(p%noris)
call facit%spiral
do i=1,p%nmembers
    call members(i)%new(p%noris)
    call members(i)%rnd_oris(p%eullims, p%trs)
end do
rt = ran_tabu(p%noris)
do i=1,p%nmembers
    call rt%reset
    call rt%ne_ran_iarr(pinds)
    do j=1,nint(0.7*real(p%noris))
        o = facit%get_ori(pinds(j))
        call members(i)%set_ori(pinds(j),o)
    end do
end do
do i=1,p%nmembers
    call members(i)%introd_alig_err(20., 0.)
end do
! set the orientations
call dockem%construct_ensemble(p%noris)
do i=1,p%nmembers
    do j=1,p%noris
        o = members(i)%get_ori(j)
        call dockem%set_ori(i,j,o)
    end do
end do
! select a member at random and print its angdoc
irnd = irnd_uni(p%nmembers)
call members(irnd)%write('member_'//trim(adjustl(dig))//'.txt')
! generate the L1-median
consensus = dockem%l1median()
! report distance
call op%new(facit, consensus)
call op%dock(angdist)
write(*,'(A,f7.4)') '>>> ANGULAR DISTANCE, SECOND:', rad2deg(angdist/real(p%noris))
! write consensus_3
call consensus%write('consensus_part'//trim(adjustl(dig))//'_2.txt')
call simple_end('**** SIMPLE_TEST_EULER_CONSENSUS NORMAL STOP ****')
end program simple_test_euler_consensus