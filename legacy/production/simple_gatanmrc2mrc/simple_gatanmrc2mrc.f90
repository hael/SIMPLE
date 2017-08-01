program simple_gatanmrc2mrc
use simple_jiffys   ! singleton
use simple_params,  only: params
use simple_imgfile, only: imgfile
use simple_image,   only: image
use simple_cmdline, only: cmdline
implicit none
type(params)                  :: p
type(cmdline)                 :: cline
type(image)                   :: frameimg
integer                       :: ldim(3), nmrcs, funit_mrc, imrc, ldim_read(3)
integer                       :: p_read, nmovies, imovie, iframe, numlen
character(len=STDLEN)         :: mrcfnam
character(len=:), allocatable :: moviename
if( command_argument_count() < 4 )then
    write(*,'(a)',advance='no') 'SIMPLE_GATANMRC2MRC filetab=<mrcfiles.txt> xdim=<nr pixel in x>'
    write(*,'(a)',advance='no') ' ydim=<nr pixels in y> smpd=<sampling distance(in A)>'
    write(*,'(a)',advance='no') ' [endian=<big|little|native{native}>] [nframes=<number of movie frames>]'
    write(*,'(a)') ' [fbody=<body of output stacks>] [numlen=<length of nr str>]'
    stop
endif
call cline%parse
call cline%checkvar('filetab', 1)
call cline%checkvar('xdim',    2)
call cline%checkvar('ydim',    3)
call cline%checkvar('smpd',    4)
call cline%check
p         = params(cline, checkdistr=.false.) ! constants & derived constants produced
ldim      = [p%xdim,p%ydim,1]
nmrcs     = nlines(p%filetab)
funit_mrc = get_fileunit()
open(unit=funit_mrc, status='old', file=p%filetab)
call frameimg%new(ldim,p%smpd)
if( cline%defined('nframes') )then
    if( .not. cline%defined('fbody') ) stop 'need fbody (file body of output stacks) on the command line'
    if( mod(nmrcs,p%nframes) .eq. 0 )then
        ! fine, go ahead
    else
        stop 'Number of mrc files not a multiple of nframes!'
    endif
    nmovies = nmrcs/p%nframes
    if( cline%defined('numlen') )then
        numlen = p%numlen
    else
        numlen = len(int2str(nmovies))
    endif
    imrc = 0
    do imovie=1,nmovies
        allocate(moviename, source=trim(adjustl(p%fbody))//int2str_pad(imovie, numlen)//'.mrcs')
        call del_binfile(moviename)
        do iframe=1,p%nframes
            imrc = imrc+1
            call progress(imrc, nmrcs)
            read(funit_mrc,'(a256)') mrcfnam
            call frameimg%read(mrcfnam,1,readhead=.false.)
            call frameimg%write(moviename,iframe)
        end do
    end do
else
    do imrc=1,nmrcs 
        read(funit_mrc,'(a256)') mrcfnam
        call frameimg%read(mrcfnam,1,readhead=.false.)
        call frameimg%write(mrcfnam,1)
    end do
endif
close(funit_mrc)
call simple_end('**** SIMPLE_GATANMRC2MRC NORMAL STOP ****')
end program simple_gatanmrc2mrc
