program simple_boxconvs
use simple_jiffys  ! singleton
use simple_defs    ! singleton
use simple_cmdline, only: cmdline
use simple_build,   only: build
use simple_params,  only: params
use simple_image,   only: image
implicit none
type(params)                       :: p
type(build)                        :: b
type(cmdline)                      :: cline
type(image)                        :: tmp
character(len=STDLEN), allocatable :: imgnames(:)
integer                            :: iimg, nimgs, ldim(3), iimg_start, iimg_stop, ifoo
logical                            :: debug=.false.
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'SIMPLE_BOXCONVS [stk=<input stack>] [filetab=<imagelist.txt>]'
    write(*,'(a)',advance='no') 'fbody=<body of output stack> [boxconvsz=<box size for boxcovolution'
    write(*,'(a)') '(pixels){256}>] [startit=<start from here>]'
    stop
endif
call cline%parse
call cline%checkvar('fbody', 1)
call cline%check
if( cline%defined('stk') .and. cline%defined('filetab') )then
    stop 'stk and filetab cannot both be defined; input either or!'
endif
if( .not. cline%defined('stk') .and. .not. cline%defined('filetab') )then
    stop 'either stk or filetab need to be defined!'
endif
if( .not. cline%defined('boxconvsz') ) call cline%set('boxconvsz', 512.)
call cline%set('prg', 'boxconvs')
p = params(cline, checkdistr=.false.)               ! constants & derived constants produced
call b%build_general_tbox( p, cline, do3d=.false. ) ! general toolbox built
! do the work
if( cline%defined('stk') )then
    call b%img%new(p%ldim, p%smpd) ! img re-generated (to account for possible non-square)
    tmp = 0.0
    do iimg=1,p%nptcls
        call b%img%read(p%stk, iimg)
        tmp = b%img%boxconv(p%boxconvsz)
        call tmp%write(trim(adjustl(p%fbody))//p%ext, iimg)
        call progress(iimg, p%nptcls)
    end do
else
    call read_filetable(p%filetab, imgnames)
    nimgs = size(imgnames)
    if( debug ) write(*,*) 'read the img filenames'
    ! get logical dimension of micrographs
    call find_ldim_nptcls(imgnames(1), ldim, ifoo)
    ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
    if( debug ) write(*,*) 'logical dimension: ', ldim
    call b%img%new(ldim, p%smpd) ! img re-generated (to account for possible non-square)
    ! determine loop range
    iimg_start = 1
    if( cline%defined('startit') ) iimg_start = p%startit
    iimg_stop  = nimgs
    if( debug ) write(*,*) 'fromto: ', iimg_start, iimg_stop
    ! do it
    tmp = 0.0
    do iimg=iimg_start,iimg_stop
        if( .not. file_exists(imgnames(iimg)) )then
            write(*,*) 'inputted img file does not exist: ', trim(adjustl(imgnames(iimg)))
        endif
        call b%img%read(imgnames(iimg), 1)
        tmp = b%img%boxconv(p%boxconvsz)
        call tmp%write(trim(adjustl(p%fbody))//p%ext, iimg)
        call progress(iimg, nimgs)
    end do
endif
call simple_end('**** SIMPLE_BOXCONVS NORMAL STOP ****')
end program