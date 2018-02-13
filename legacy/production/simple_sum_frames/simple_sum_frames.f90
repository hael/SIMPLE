program simple_sum_frames
use simple_cmdline    ! singleton
use simple_defs       ! singleton
use simple_build,     only: build
use simple_params,    only: params, find_ldim
use simple_jiffys,    only: nlines, get_fileunit, simple_end, same_format
use simple_imgfile,   only: imgfile
use simple_image,     only: image
use simple_oris,      only: oris
implicit none
type(params)          :: p
type(build)           :: b
integer               :: nmovies, funit_movies, funit_sumfiles, movie
integer               :: ldim(3), nframes, frame, nsumfiles, i
character(len=STDLEN) :: moviename, sumfilename
type(image)           :: img_frame, img_sum
type(imgfile)         :: imgstk
logical, parameter    :: debug = .true.
if( command_argument_count() < 2 )then
    write(*,'(a)') 'SIMPLE_SUM_FRAMES filetab=<movies.txt> sumtab=<sumfilenames.txt> [vis=<yes|no{no}>]'
    stop
endif
call parse_cmdline
call cmdcheckvar('filetab', 1)
call cmdcheckvar('sumtab',  2)
call cmdcheck
p = params() ! constants & derived constants produced
nmovies   = nlines(p%filetab)
nsumfiles = nlines(p%sumtab)
if( nmovies /= nsumfiles ) stop 'number of entries in inputted files do not match!'
funit_movies = get_fileunit()
open(unit=funit_movies, status='old', file=p%filetab)
funit_sumfiles = get_fileunit()
open(unit=funit_sumfiles, status='old', file=p%sumtab)
! LOOP OVER EXPOSURES (MOVIES)
do movie=1,nmovies
    ! READ & ALLOCATE
    read(funit_movies,'(a256)')   moviename
    read(funit_sumfiles,'(a256)') sumfilename
    if( .not. same_format(moviename,sumfilename) )then
        stop 'incompatible formats in filetab & sumtab!'
    endif
    ldim = find_ldim(moviename)
    call img_frame%new([ldim(1),ldim(2),1], p%smpd)
    call img_sum%new([ldim(1),ldim(2),1], p%smpd)
    ! GET NUMBER OF FRAMES FROM STACK
    call imgstk%open(moviename)
    nframes = imgstk%getStackSz()
    call imgstk%close
    if( debug ) write(*,*) 'number of frames: ', nframes
    ! SUM THE FRAMES
    img_sum = 0.
    do i=1,nframes
        ! READ FRAME
        call img_frame%read(moviename,i)
        if( p%vis .eq. 'yes' ) call img_frame%vis()
        ! AGGLOMERATE SUM
        call img_sum%add(img_frame)
        if( p%vis .eq. 'yes' ) call img_sum%vis()
    end do
    call img_sum%write(sumfilename, 1)
    if( p%vis .eq. 'yes' ) exit
end do
close(funit_movies)
call simple_end('**** SIMPLE_SUM_FRAMES NORMAL STOP ****')
end program simple_sum_frames
