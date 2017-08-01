program simple_simsubtomo
use simple_jiffys   ! singleton
use simple_cmdline, only: cmdline
use simple_build,   only: build
use simple_params,  only: params
use simple_image,   only: image
use simple_ori,     only: ori
implicit none
type(params)  :: p
type(build)   :: b
type(cmdline) :: cline
type(image)   :: vol_rot
type(ori)     :: o
integer       :: iptcl, numlen
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'SIMPLE_SIMSUBTOMO vol1=<reference.ext> smpd=<sampling distance(in A)>'
    write(*,'(a)') ' nptcls=<nr of subtomos> snr=<signal to noise ratio> [nthr=<nr of openMP threads{1}>]'
    stop
endif
call cline%parse
call cline%checkvar('vol1', 1)
call cline%checkvar('smpd', 3)
call cline%checkvar('msk',  4)
call cline%checkvar('snr',  5)
call cline%check
p = params(cline, checkpara=.false.) ! constants & derived constants produced, mode=2
call b%build_general_tbox(p, cline)  ! general objects built
call vol_rot%new([p%box,p%box,p%box], p%smpd)
call b%vol%new([p%box,p%box,p%box],   p%smpd)
call b%vol%read(p%vols(1))
call b%vol%add_gauran(p%snr)
call o%new
numlen = len(int2str(p%nptcls))
do iptcl=1,p%nptcls
    call o%rnd_ori
    vol_rot = b%proj%rotvol(b%vol, o, p)
    call vol_rot%write('subtomo'//int2str_pad(iptcl,numlen)//p%ext)
end do    
end program