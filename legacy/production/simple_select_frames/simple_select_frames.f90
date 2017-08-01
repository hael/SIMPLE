!<select_frames/begin> is a program for selecting contiguous segments of frames from DDD movies <select_frames/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2016-04-07.
!
program simple_select_frames
use simple_defs       ! singleton
use simple_jiffys     ! singleton
use simple_cmdline,   only: cmdline
use simple_build,     only: build
use simple_params,    only: params
use simple_imgfile,   only: imgfile
use simple_image,     only: image
implicit none
type(params)                       :: p
type(build)                        :: b
type(cmdline)                      :: cline
integer                            :: nmovies, nframes, frame, numlen, ldim(3)
integer                            :: fromto(2), movie, cnt, cnt2, ntot, lfoo(3), ifoo
character(len=STDLEN), allocatable :: movienames(:)
character(len=:), allocatable      :: new_name
type(image)                        :: img_frame
logical, parameter                 :: debug = .false.
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'SIMPLE_SELECT_FRAMES filetab=<movies.txt> fbody=<body of output files>'
    write(*,'(a)',advance='no') ' fromf=<from frame> tof=<to frame> smpd=<sampling distance(in A)>'
    write(*,'(a)') ' [startit=<start from here>]'
    stop
endif
call cline%parse
call cline%checkvar('filetab', 1)
call cline%checkvar('fbody',   2)
call cline%checkvar('fromf',   3)
call cline%checkvar('tof',     4)
call cline%checkvar('smpd',    5)
call cline%check
p = params(cline, checkpara=.false.) ! constants & derived constants produced
call b%build_general_tbox(p,cline,do3d=.false.)
call read_filetable(p%filetab, movienames)
nmovies = size(movienames)
! FIND LDIM AND NUMLEN (LENGTH OF NUMBER STRING)
if( cline%defined('startit') )then
    call find_ldim_nptcls(movienames(p%startit), ldim, ifoo)
endif
if( debug ) write(*,*) 'logical dimension: ', ldim
ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
numlen = len(int2str(nmovies))
if( debug ) write(*,*) 'length of number string: ', numlen
! DETERMINE LOOP RANGE
if( cline%defined('part') )then
    if( cline%defined('fromp') .and. cline%defined('top') )then
        fromto(1) = p%fromp
        fromto(2) = p%top
        ntot = fromto(2)-fromto(1)+1
    else
        stop 'fromp & top args need to be defined in parallel execution; simple_select_frames'
    endif
else
    fromto(1) = 1
    if( cline%defined('startit') ) fromto(1) = p%startit
    fromto(2) = nmovies
    ntot      = nmovies
endif
if( debug ) write(*,*) 'fromto: ', fromto(1), fromto(2)
call img_frame%new([ldim(1),ldim(2),1], p%smpd)
! LOOP OVER EXPOSURES (MOVIES)
cnt2 = 0
do movie=fromto(1),fromto(2)
    if( .not. file_exists(movienames(movie)) )then
        write(*,*) 'inputted movie stack does not exist: ', trim(adjustl(movienames(movie)))
    endif
    cnt2 = cnt2+1
    ! GET NUMBER OF FRAMES FROM STACK
    call find_ldim_nptcls(movienames(movie), lfoo, nframes)
    if( debug ) write(*,*) 'number of frames: ', nframes
    ! CREATE NEW NAME
    allocate(new_name, source=trim(adjustl(p%fbody))//int2str_pad(movie, numlen)//p%ext)
    cnt = 0
    do frame=p%fromf,p%tof
        cnt = cnt+1
        call img_frame%read(movienames(movie),frame)
        call img_frame%write(new_name,cnt)
    end do
    deallocate(new_name)
    write(*,'(f4.0,1x,a)') 100.*(real(cnt2)/real(ntot)), 'percent of the movies processed'
end do
call img_frame%kill
call simple_end('**** SIMPLE_SELECT_FRAMES NORMAL STOP ****')
end program simple_select_frames