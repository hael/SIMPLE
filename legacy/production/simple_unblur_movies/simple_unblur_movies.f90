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
use simple_defs     ! singleton
use simple_jiffys   ! singleton
use simple_unblur   ! singleton
use simple_cmdline, only: cmdline
use simple_build,   only: build
use simple_params,  only: params
use simple_imgfile, only: imgfile
use simple_oris,    only: oris
use simple_image,   only: image
use simple_timing
use simple_cuda_defs
use simple_cuda
implicit none
type(params)                       :: p
type(build)                        :: b
type(cmdline)                      :: cline
type(image)                        :: movie_sum, pspec_sum, pspec_corrected, movie_sum_ctf
type(image)                        :: movie_sum_corrected, pspec_half_n_half
integer                            :: nmovies, imovie, imovie_start, imovie_stop, file_stat
integer                            :: funit_movies, cnt, numlen,  ntot, alloc_stat, fnr
character(len=:), allocatable      :: cpcmd, new_name
character(len=STDLEN), allocatable :: movienames(:)
real                               :: corr
logical                            :: debug=.true.
! CUDA err variable for the return function calls
integer                            :: err
call timestamp()
! call start_Alltimers_cpu()
! starting the cuda environment
call simple_cuda_init(err)
if( err .ne. RC_SUCCESS ) write(*,*) 'cublas init failed'
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'SIMPLE_UNBLUR_MOVIES filetab=<movies.txt> smpd=<sampling distance(in A)>'
    write(*,'(a)',advance='no') ' [fbody=<body of output files>] [lpstart=<low-pass limit{15}>] [lpstop=<low'
    write(*,'(a)',advance='no') '-pass limit{8}>] [trs=<maximum halfwidth shift(in pixels){5}>]'
    write(*,'(a)',advance='no') ' [pspecsz=<box size for boxcovolution(pixels){512}>]'
    write(*,'(a)',advance='no') ' [use_gpu=<yes|no{no}>] [numlen=<length of nr str>]'
    write(*,'(a)') ' [nthr=<nr of OpenMP threads{1}>] [startit=<start from here>] [shrink=<shrink factor{1}>]'
    stop
endif
call cline%parse
call cline%checkvar('filetab', 1)
call cline%checkvar('smpd',    2)
if( .not. cline%defined('trs') )then
    call cline%set('trs', 5.)
endif
if( .not. cline%defined('lpstart') )then
    call cline%set('lpstart', 15.)
endif
if( .not. cline%defined('lpstop') )then
    call cline%set('lpstop', 8.)
endif
call cline%check
p = params(cline, checkpara=.false.) ! constants & derived constants produced
call b%build_general_tbox(p,cline,do3d=.false.)
! read movienames
nmovies = nlines(p%filetab)
if( debug ) write(*,*) 'nmovies: ', nmovies
allocate( movienames(nmovies), stat=alloc_stat )
call alloc_err('In: simple_integrate_movies', alloc_stat)
funit_movies = get_fileunit()
open(unit=funit_movies, status='old', file=p%filetab)
do imovie=1,nmovies
    read(funit_movies,'(a256)') movienames(imovie)
end do
close(funit_movies)
if( debug ) write(*,*) 'read the movie filenames'
! DETERMINE LOOP RANGE
if( cline%defined('part') )then
    if( cline%defined('fromp') .and. cline%defined('top') )then
        imovie_start = p%fromp
        imovie_stop  = p%top
    else
        stop 'fromp & top args need to be defined in parallel execution; simple_unblur_movies'
    endif
else
    imovie_start = 1
    if( cline%defined('startit') ) imovie_start = p%startit
    imovie_stop  = nmovies
endif
if( debug ) write(*,*) 'fromto: ', imovie_start, imovie_stop
ntot = imovie_stop-imovie_start+1
if( cline%defined('numlen') )then
    numlen = p%numlen
else
    numlen = len(int2str(nmovies))
endif
if( debug ) write(*,*) 'length of number string: ', numlen
! create output orientations
call b%a%new(ntot)
! align
cnt = 0
do imovie=imovie_start,imovie_stop
    if( .not. file_exists(movienames(imovie)) )then
        write(*,*) 'inputted movie stack does not exist: ', trim(adjustl(movienames(imovie)))
    endif
    write(*,'(a,1x,i5)') '>>> PROCESSING MOVIE:', imovie
    cnt = cnt+1
    call unblur_movie(movienames(imovie), p, cnt, b%a, corr, movie_sum, movie_sum_corrected, movie_sum_ctf)
    if( debug ) print *, 'ldim of output movie_sum:           ', movie_sum%get_ldim()
    if( debug ) print *, 'ldim of output movie_sum_corrected: ', movie_sum_corrected%get_ldim()
    if( debug ) print *, 'ldim of output movie_sum_ctf      : ', movie_sum_ctf%get_ldim()
    if( cline%defined('fbody') )then
        call movie_sum_corrected%write(trim(adjustl(p%fbody))//'_intg'//int2str_pad(imovie, numlen)//p%ext)
        call movie_sum_ctf%write(trim(adjustl(p%fbody))//'_ctf'//int2str_pad(imovie, numlen)//p%ext)
    else
        call movie_sum_corrected%write(int2str_pad(imovie, numlen)//'_intg'//p%ext)
        call movie_sum_ctf%write(int2str_pad(imovie, numlen)//'_ctf'//p%ext)
    endif
    pspec_sum         = movie_sum%mic2spec(p%pspecsz, trim(adjustl(p%speckind)))
    pspec_corrected   = movie_sum_corrected%mic2spec(p%pspecsz, trim(adjustl(p%speckind)))
    pspec_half_n_half = pspec_sum%before_after(pspec_corrected)
    if( debug ) print *, 'ldim of pspec_sum:         ', pspec_sum%get_ldim()
    if( debug ) print *, 'ldim of pspec_corrected:   ', pspec_corrected%get_ldim()
    if( debug ) print *, 'ldim of pspec_half_n_half: ', pspec_half_n_half%get_ldim()
    if( cline%defined('fbody') )then
        call pspec_half_n_half%write(trim(adjustl(p%fbody))//'_pspec'//int2str_pad(imovie, numlen)//p%ext)
    else
        call pspec_half_n_half%write(int2str_pad(imovie, numlen)//'_pspec'//p%ext)
    endif
    call movie_sum%kill
    call movie_sum_corrected%kill
    call movie_sum_ctf%kill
    call pspec_sum%kill
    call pspec_corrected%kill
    call pspec_half_n_half%kill
    ! RENAME MOVIE
    if( cline%defined('fbody') )then
        allocate(new_name, source=trim(adjustl(p%fbody))//int2str_pad(imovie, numlen)//p%ext)
    else
        allocate(new_name, source=int2str_pad(imovie, numlen)//p%ext)
    endif
    if( .not. file_exists(new_name))then
        allocate(cpcmd, source='cp '//trim(adjustl(movienames(imovie)))//' ./'//new_name)
        call system(cpcmd)
        deallocate(cpcmd)
    endif
    deallocate(new_name)
    if( .not. cline%defined('part') ) call b%a%write(imovie, 'unblur_movies_params.txt')
    write(*,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the movies processed'
end do
if( cline%defined('part') )then
    call b%a%write('unblur_movies_params'//'_part'//int2str_pad(p%part,p%numlen)//'.txt')
    fnr = get_fileunit()
    open(unit=fnr, FILE='JOB_FINISHED_'//int2str_pad(p%part,p%numlen), STATUS='REPLACE', action='WRITE', iostat=file_stat)
    call fopen_err( 'In: simple_unblur_movies', file_stat )
    write(fnr,'(A)') '**** SIMPLE_COMLIN_SMAT NORMAL STOP ****'
    close(fnr)
endif
call simple_end('**** SIMPLE_UNBLUR_MOVIES NORMAL STOP ****')
!*******************************************************************************
!    Environment shutdown
!
!*******************************************************************************
!shutting down the environment
call simple_cuda_shutdown()
!shutting down the timers
! call stop_Alltimers_cpu()
end program simple_unblur_movies
