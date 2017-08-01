!==Program simple_unblur_movies
!
! <unblur_movies/begin> is a program for movie alignment or unblurring.
! Input is a textfile with absolute paths to movie files in addition to a few obvious input
! parameters. Output is (x,y) shift parameters for every frame of the movie. <unblur_movies/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_unblur_movies
use simple_cmdline  ! singleton
use simple_defs     ! singleton
use simple_params,  only: params, find_ldim
use simple_imgfile, only: imgfile
use simple_jiffys,  only: simple_end, nlines, get_fileunit, int2str
use simple_unblur   ! singleton
use simple_oris,    only: oris
implicit none
type(params)                       :: p
integer                            :: nmovies, imovie, iframe, funit_movies, cnt, imovie_start, imovie_stop, ntot
character(len=STDLEN), allocatable :: movienames(:)
real                               :: corr
real, allocatable                  :: shifts(:,:)
type(oris)                         :: outoris
logical                            :: debug=.true.
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'SIMPLE_UNBLUR_MOVIES filetab=<movies.txt> smpd=<sampling distance(in A)>'
    write(*,'(a)',advance='no') ' box=<box size(in pixels)> [opt=<optimizer(bfgs|powell|simplex|oasis|'
    write(*,'(a)',advance='no') 'bforce){simplex}>] [trs=<maximum halfwidth shift(in pixels){5}>] '
    write(*,'(a)') '[nrestarts=<number of restarts{3}>] [nthr=<nr of OpenMP threads{1}>]'
    stop
endif
call parse_cmdline
call cmdcheckvar('filetab', 1)
call cmdcheckvar('smpd',    2)
call cmdcheckvar('box',     3)
if( .not. defined_cmd_arg('trs') )then
    call set_cmdline('trs', 5.)
endif
if( .not. defined_cmd_arg('nrestarts') )then
    call set_cmdline('nrestarts', 3.)
endif
call cmdcheck
p = params(checkpara=.false.) ! constants & derived constants produced
! read movienames
nmovies = nlines(p%filetab)
allocate(movienames(nmovies))
funit_movies = get_fileunit()
open(unit=funit_movies, status='old', file=p%filetab)
do imovie=1,nmovies
    read(funit_movies,*) movienames(imovie)
end do
close(funit_movies)
if( defined_cmd_arg('part') )then
    imovie_start = p%fromp
    imovie_stop  = p%top
else
    imovie_start = 1
    imovie_stop  = nmovies
endif
ntot = imovie_stop-imovie_start+1
! create output orietntations
call outoris%new(ntot)
! align
cnt = 0
do imovie=imovie_start,imovie_stop
    cnt = cnt+1
    call unblur_init( movienames(imovie), p%opt, p%smpd, p%trs, p%box, p%nrestarts )
    corr = unblur_minimize()
    write(*,'(a,1x,f7.4,1x,a,i5)') '>>> OPTIMAL CORRELATION:', corr, 'MOVIE:', imovie
    shifts = get_unblur_shifts()
    do iframe=1,size(shifts,1)
        call outoris%set(cnt, 'x'//int2str(iframe), shifts(iframe,1))
        call outoris%set(cnt, 'y'//int2str(iframe), shifts(iframe,2))
    end do
    write(*,*) 100.*(real(cnt)/real(ntot)), ' percent of the movies processed'
end do
if( defined_cmd_arg('part') )then
    call outoris%write('unblur_movies_params'//'_part'//int2str(p%part)//'.txt')
else
    call outoris%write('unblur_movies_params.txt')
endif
call simple_end('**** SIMPLE_UNBLUR_MOVIES NORMAL STOP ****')
end program simple_unblur_movies
