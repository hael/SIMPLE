!==Program simple_phflip
!
! <simple_phflip/begin> is a program for phase flipping of movies/micrographs <simple_phflip/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2011-08-16.
!
program simple_phflip
use simple_cmdline    ! singleton
use simple_defs       ! singleton
use simple_jiffys     ! singleton
use simple_build,     only: build
use simple_params,    only: params
use simple_image,     only: image
implicit none
type(params)                       :: p
type(build)                        :: b
integer                            :: nmovies, funit_movies, nframes, frame, alloc_stat
integer                            :: ldim(3), movie, niter, fromto(2), ntot, ifoo, lfoo(3)
character(len=STDLEN)              :: mode
character(len=:), allocatable      :: flipname
character(len=STDLEN), allocatable :: movienames(:), boxfilenames(:)
type(image)                        :: img_frame
logical, parameter                 :: debug = .false.
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'SIMPLE_PHFLIP filetab=<movies.txt> smpd=<sampling'
    write(*,'(a)',advance='no') ' distance(in A)> oritab=<CTF doc> [neg=<yes|no{yes}>]'
    write(*,'(a)',advance='no') ' [kv=<acceleration voltage(kV){300.}>]'
    write(*,'(a)',advance='no') ' [cs=<spherical aberration constant(mm){2.7}>]'
    write(*,'(a)') ' [fraca=<frac amp contrast{0.07}>]'
    stop
endif
call parse_cmdline
call cmdcheckvar('filetab', 1)
call cmdcheckvar('smpd',    2)
call set_cmdline('ctf', 'flip')
if( .not. defined_cmd_arg('neg') )then
    call set_cmdline('neg', 'yes')
endif
call cmdcheck
p = params(checkpara=.false.) ! constants & derived constants produced
call b%build_general_tbox(p,do3d=.false.)
nmovies = nlines(p%filetab)
if( debug ) write(*,*) 'nmovies: ', nmovies
funit_movies = get_fileunit()
open(unit=funit_movies, status='old', file=p%filetab)
allocate( movienames(nmovies), stat=alloc_stat )
call alloc_err('In: simple_phflip; boxdata etc., 1', alloc_stat)

! read the filenames
do movie=1,nmovies
    read(funit_movies,'(a256)') movienames(movie)
end do
close(funit_movies)
if( debug ) write(*,*) 'read the filenames'

! determine loop range
if( defined_cmd_arg('part') )then
    if( defined_cmd_arg('fromp') .and. defined_cmd_arg('top') )then
        fromto(1) = p%fromp
        fromto(2) = p%top
    else
        stop 'fromp & top args need to be defined in parallel execution; simple_extr_ptcls'
    endif
else
    fromto(1) = 1
    fromto(2) = nmovies
endif
if( debug ) write(*,*) 'fromto: ', fromto(1), fromto(2)

! find ldim
call find_ldim_nptcls(movienames(1), ldim, ifoo)
if( debug ) write(*,*) 'ldim: ', ldim

! create frame
call img_frame%new([ldim(1),ldim(2),1], p%smpd)

! set CTF mode
if( p%neg .eq. 'yes' )then
    mode = 'flipneg'
else
    mode = 'flip'
endif
if( debug ) write(*,*) 'ctfmode: ', trim(mode)

! loop over exposures (movies)
niter = 0
do movie=fromto(1),fromto(2)
    
    ! show progress
    if( .not. defined_cmd_arg('part') .and. niter > 1 )then
        call progress(movie,nmovies)
    endif
    
    ! update iteration counter
    niter = niter+1
    
    ! get number of frames from stack
    call find_ldim_nptcls(movienames(movie), lfoo, nframes)
    if( debug ) write(*,*) 'nframes: ', nframes
    
    ! allocate the name of the new stack
    allocate(flipname, source=add2fbody(trim(adjustl(movienames(movie))), p%ext, 'flip'))
    if( debug ) write(*,*) flipname
    
    ! loop over frames
    do frame=1,nframes
        ! read frame
        call img_frame%read(movienames(movie),frame)
        call img_frame%fwd_ft
        if( b%a%isthere('dfy') )then ! astigmatic CTF
            call b%tfun%apply(img_frame, b%a%get(movie,'dfx'), trim(adjustl(mode)),&
            dfy=b%a%get(movie,'dfy'), angast=b%a%get(movie,'angast'))
        else                         ! non-astigmatic CTF
            call b%tfun%apply(img_frame, b%a%get(movie,'dfx'), trim(adjustl(mode)))
        endif
        call img_frame%bwd_ft
        call img_frame%write(flipname,frame)
    end do
    deallocate(flipname)
end do
call simple_end('**** SIMPLE_PHFLIP NORMAL STOP ****')
end program simple_phflip
