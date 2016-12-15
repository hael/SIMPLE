!==Program simple_cluster_smat
!
! <cluster\_smat/begin> is a program for clustering a similarity matrix and use an combined cluster validation
! index to assess the quality of the clustering based on the number of clusters. <cluster\_smat/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2016
!
program simple_cluster_smat
use simple_defs           ! singleton
use simple_jiffys,        ! singleton
use simple_cmdline,       only: cmdline
use simple_jiffys,        only: simple_end
use simple_shc_cluster,   only: shc_cluster
use simple_oris,          only: oris
use simple_build,         only: build
use simple_params,        only: params
use simple_cluster_valid, only: cluster_valid
implicit none
type(params)         :: p
type(build), target  :: b
type(cmdline)        :: cline
type(shc_cluster)    :: shcc
type(cluster_valid)  :: cvalid
real, allocatable    :: smat(:,:)
integer              :: funit, io_stat, ncls_min, loc(1), ncls_stop, icls
integer              :: pop, ncls, alloc_stat, irestart, ntot, cnt=0, numlen
real                 :: avg_ratio, min_ratio, ratio, x, sim
real, allocatable    :: validinds(:)
integer, parameter   :: NRESTARTS=10
logical              :: debug=.false., done=.false.
if( command_argument_count() < 2 )then
    write(*,'(a)', advance='no') 'SIMPLE_CLUSTER_SMAT nptcls=<nr particles> fname=<smat.bin> ncls=<max nr'
    write(*,'(a)', advance='no') ' clusters to test> label=<class|state|subclass{class}>'
    write(*,'(a)') ' [nthr=<nr OpenMP threads{1}>]'
    stop
endif
call cline%parse
call cline%checkvar('nptcls', 1)
call cline%checkvar('fname',  2)
call cline%checkvar('ncls',   3)
call cline%checkvar('label',  4)
call cline%check    ! checks args and prints cmdline to cmdline.dat
call cline%set('prg', 'cluster_smat')
p = params(cline,.false.) ! constants & derived constants produced
call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built (no oritab reading)
! obtain similarity matrix
allocate(smat(p%nptcls,p%nptcls), stat=alloc_stat)
call alloc_err('In: simple_cluster_smat, 1', alloc_stat)
smat = 1.
funit = get_fileunit()
open(unit=funit, status='OLD', action='READ', file=p%fname, access='STREAM')
read(unit=funit,pos=1,iostat=io_stat) smat
if( io_stat .ne. 0 )then
    write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when reading: ', p%fname
    stop 'I/O error; simple_cluster_smat'
endif
close(funit)
allocate(validinds(2:p%ncls), stat=alloc_stat)
call alloc_err("In: simple_cluster_smat", alloc_stat)
validinds = 0
ntot = (p%ncls-1)*NRESTARTS
cnt = 0
numlen = len(int2str(p%ncls))
do ncls=2,p%ncls
    avg_ratio = 0.
    min_ratio = huge(x)
    do irestart=1,NRESTARTS
        cnt = cnt+1
        call progress(cnt,ntot)
        call shcc%new(p%nptcls, ncls, smat, b%a)
        call shcc%shc(.false., p%label, sim)
        call cvalid%new(b%a, ncls, p%label, smat)
        ratio = cvalid%elmlunds_ratio_index()
        avg_ratio = avg_ratio+ratio
        if( ratio < min_ratio )then
            min_ratio = ratio
            call b%a%write('shc_clustering_ncls'//int2str_pad(ncls,numlen)//'.txt')
        endif
    end do
    validinds(ncls) = avg_ratio/real(NRESTARTS)
end do
ncls_stop = 0
done = .false. 
do ncls=2,p%ncls
    write(*,'(a,1x,f9.3,8x,a,1x,i3)') 'COHESION/SEPARATION RATIO INDEX: ', validinds(ncls), ' NCLS: ', ncls
    call b%a%read('shc_clustering_ncls'//int2str_pad(ncls,numlen)//'.txt')
    do icls=1,ncls
        pop = b%a%get_pop(icls, p%label)
        write(*,'(a,3x,i5,1x,a,1x,i3)') '  CLUSTER POPULATION:', pop, 'CLUSTER:', icls
    end do
    write(*,'(a)') '***************************************************'
    if( ncls < p%ncls )then
        if( validinds(ncls+1) >= validinds(ncls) .and. .not. done )then
            ncls_stop = ncls
            done = .true.
        endif
    endif
end do
loc = minloc(validinds)
ncls_min = loc(1)+1
write(*,'(a,i3)') 'NUMBER OF CLUSTERS FOUND BY MINIMIZING ELMLUNDS RATIO INDEX: ', ncls_min
if( ncls_stop /= 0 )then
    write(*,'(a,i3)') 'NUMBER OF CLUSTERS FOUND BY STOPPING CRITERIUM:              ', ncls_stop
endif
call simple_end('**** SIMPLE_CLUSTER_SMAT NORMAL STOP ****')
end program simple_cluster_smat