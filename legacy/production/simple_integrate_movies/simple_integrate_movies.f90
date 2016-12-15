!==Program simple_integrate_movies
!
!<integrate_movies/begin> is a program for integrating DDD movies <integrate_movies/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2011-08-16.
!
program simple_integrate_movies
use simple_build,     only: build
use simple_params,    only: params
use simple_imgfile,   only: imgfile
use simple_image,     only: image
use simple_cmdline,   only: cmdline
use simple_defs       ! singleton
use simple_jiffys     ! singleton
implicit none
type(params)                       :: p
type(build)                        :: b
type(cmdline)                      :: cline
integer                            :: nmovies, nframes, frame, alloc_stat, lfoo(3)
integer                            :: numlen, ldim(3), fromto(2), movie, ifoo
character(len=STDLEN), allocatable :: movienames(:)
character(len=:), allocatable      :: cpcmd, new_name
real                               :: x, y
type(image), allocatable           :: img_frames(:)
type(image)                        :: img_sum, pspec
logical, parameter                 :: debug = .false.
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'SIMPLE_INTEGRATE_MOVIES filetab=<movies.txt> fbody=<body of output files>'
    write(*,'(a)',advance='no') ' smpd=<sampling distance(in A)> [oritab=<shift parameters 4 frames>]'
    write(*,'(a)') ' [pspecsz=<box size for boxcovolution(pixels){512}>]'
    stop
endif
call cline%parse
call cline%checkvar('filetab', 1)
call cline%checkvar('fbody',   2)
call cline%checkvar('smpd',    3)
call cline%check
if( .not. cline%defined('pspecsz') )then
    call cline%set('pspecsz', 512.)
endif
call cline%set('prg', 'integrate_movies')
p = params(cline,checkdistr=.false.) ! constants & derived constants produced
call b%build_general_tbox(p,cline,do3d=.false.)
call read_filetable(p%filetab, movienames)
nmovies = size(movienames)
if( debug ) write(*,*) 'read the movie filenames'
! FIND LDIM AND NUMLEN (LENGTH OF NUMBER STRING)
call find_ldim_nptcls(movienames(1), ldim, ifoo)
if( debug ) write(*,*) 'logical dimension: ', ldim
ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
numlen = len(int2str(nmovies))
if( debug ) write(*,*) 'length of number string: ', numlen
! DETERMINE LOOP RANGE
if( cline%defined('part') )then
    if( cline%defined('fromp') .and. cline%defined('top') )then
        fromto(1) = p%fromp
        fromto(2) = p%top
    else
        stop 'fromp & top args need to be defined in parallel execution; simple_integrate_movies'
    endif
else
    fromto(1) = 1
    fromto(2) = nmovies
endif
if( debug ) write(*,*) 'fromto: ', fromto(1), fromto(2)
! CREATE SUM
call img_sum%new([ldim(1),ldim(2),1], p%smpd)
! LOOP OVER EXPOSURES (MOVIES)
do movie=fromto(1),fromto(2)
    if( .not. file_exists(movienames(movie)) )then
        write(*,*) 'inputted movie stack does not exist: ', trim(adjustl(movienames(movie)))
    endif
    ! GET NUMBER OF FRAMES FROM STACK
    call find_ldim_nptcls(movienames(movie), lfoo, nframes)
    if( debug ) write(*,*) 'number of frames: ', nframes
    ! CREATE FRAMES & READ
    allocate(img_frames(nframes), stat=alloc_stat)
    call alloc_err('In: simple_integrate_movies', alloc_stat)
    img_sum = 0.
    do frame=1,nframes
        call img_frames(frame)%new([ldim(1),ldim(2),1], p%smpd)
        call img_frames(frame)%read(movienames(movie),frame)
        if( cline%defined('oritab') )then
            ! SHIFT FRAME ACCORDING TO GLOBAL SHIFT (DRIFT CORRECTION)
            if( b%a%isthere(movie, 'x'//int2str(frame)) .and. b%a%isthere(movie, 'y'//int2str(frame)) )then
                call img_frames(frame)%fwd_ft
                x = b%a%get(movie, 'x'//int2str(frame))
                y = b%a%get(movie, 'y'//int2str(frame))
                if( debug ) print *, 'shifting frame: ', x, y
                call img_frames(frame)%shift(-x,-y)
                call img_frames(frame)%bwd_ft
                call img_sum%add(img_frames(frame))
            else
                stop 'no shift parameters available for alignment in oritab'
            endif
        else
            call img_sum%add(img_frames(frame))
        endif
    end do
    ! RENAME MOVIE
    allocate(new_name, source=trim(adjustl(p%fbody))//int2str_pad(movie, numlen)//p%ext)
    if( .not. file_exists(new_name))then
        allocate(cpcmd, source='cp '//trim(adjustl(movienames(movie)))//' ./'//new_name)
        call system(cpcmd)
        deallocate(cpcmd)
    endif
    ! DESTROY OBJECTS AND DEALLOCATE
    do frame=1,nframes
        call img_frames(frame)%kill
    end do
    deallocate(new_name,img_frames)
    call img_sum%write(trim(adjustl(p%fbody))//'_intg'//int2str_pad(movie, numlen)//p%ext)
    pspec = img_sum%mic2spec(p%pspecsz, trim(adjustl(p%speckind)))
    call pspec%write(trim(adjustl(p%fbody))//'_pspec'//int2str_pad(movie, numlen)//p%ext)
end do
call simple_end('**** SIMPLE_INTEGRATE_MOVIES NORMAL STOP ****')
end program simple_integrate_movies
