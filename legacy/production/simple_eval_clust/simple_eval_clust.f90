program simple_eval_clust
use simple_oris,   only: oris
use simple_ori,    only: ori
use simple_hac,    only: hac
use simple_params, only: params
use simple_jiffys, only: simple_end, alloc_err, nlines
use simple_cmdline ! singleton
use simple_math,   only: hpsort
implicit none
type(params)         :: p
type(oris)           :: o
integer              :: i, cnt, ncls, alloc_stat
real                 :: avg, sdev, med, sdevmed, avg10best, sdev10best, avg10worst, sdev10worst
real, allocatable    :: avgs(:), sdevs(:)
integer, allocatable :: order(:)
if( command_argument_count() < 1 )then
    write(*,'(a)', advance='no') 'SIMPLE_EVAL_CLUST oritab=<SIMPLE alignment doc>'
    write(*,'(a)') ' clsdoc=<SPIDER clustering doc> [minp=<minimum ptcls in cluster{10}>]'
    stop
endif
call parse_cmdline
call cmdcheckvar('oritab', 1)
call cmdcheckvar('clsdoc', 2)
call cmdcheck ! checks args and prints cmdline to cmdline.dat
p = params() ! constants & derived constants produced
o = oris(nlines(p%clsdoc))
call o%read_clsdoc(p%clsdoc)
call o%read(p%oritab)
call o%zero('e3')
ncls = o%get_ncls(p%minp)
allocate(avgs(ncls), sdevs(ncls), order(ncls), stat=alloc_stat)
call alloc_err('utst_calc_stat; simple_cluster', alloc_stat)
! collect statistics
cnt = 0
do i=1,o%get_ncls()
    if( o%get_clspop(i) >= p%minp )then
        cnt = cnt+1
        order(cnt) = cnt
        call o%class_dist_stat(i, avgs(cnt), sdevs(cnt))
    endif
end do
call hpsort(ncls, avgs, order)
! median
med = avgs(ncls/2)
sdevmed = sdevs(order(ncls/2))
! average
avg  = sum(avgs)/real(ncls)
sdev = sum(sdevs)/real(ncls)
! 10 best
avg10best  = 0.
sdev10best = 0.
do i=1,10
    avg10best  = avg10best+avgs(i)
    sdev10best = sdev10best+sdevs(order(i))
end do
avg10best  = avg10best/10.
sdev10best = sdev10best/10.
! 10 worst
avg10worst  = 0.
sdev10worst = 0.
do i=ncls,ncls-9,-1
    avg10worst  = avg10worst+avgs(i)
    sdev10worst = sdev10worst+sdevs(order(i))
end do
avg10worst  = avg10worst/10.
sdev10worst = sdev10worst/10.
write(*,'(a)') '**** TEST RESULTS, ANGULAR SPREAD WITHIN CLASSES ****'
write(*,'(a)') '|  median   |  average  |  10 best  |  10 worst |'
write(*,'(a)') '|avg    sdev|avg    sdev|avg    sdev|avg    sdev|'
write(*,'(1x,f4.1,3x,f4.1,1x,f4.1,3x,f4.1,1x,f4.1,3x,f4.1,1x,f4.1,3x,f4.1)') &
med, sdevmed, avg, sdev, avg10best, sdev10best, avg10worst, sdev10worst
deallocate( avgs, sdevs, order )
end program simple_eval_clust