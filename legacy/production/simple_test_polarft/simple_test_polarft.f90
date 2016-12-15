program simple_test_polarft
use simple_build,   only: build
use simple_params,  only: params
use simple_cmdline, only: cmdline
use simple_polarft_tester
use simple_jiffys
type(params)  :: p
type(build)   :: b
type(cmdline) :: cline
integer       :: i
real          :: lp20e,lp15e,lp10e,lp7e,snr02e,snr01e,snr005e,snr001e
real          :: fslp20e,fslp15e,fslp10e,fslp7e,fssnr02e,fssnr01e,fssnr005e,fssnr001e
if( command_argument_count() < 2 )then
    write(*,'(a)', advance='no') 'SIMPLE_TEST_POLARFT stk=stack.spi box=<image size(in pixels)>'
    write(*,'(a)') ' smpd=<sampling distance(in A)> msk=<mask radius> [fraczero=<fraction of zeros>]'
endif
call cline%parse
call cline%checkvar('stk',  1)
call cline%checkvar('box',  2)
call cline%checkvar('smpd', 3)
call cline%checkvar('msk',  4)
call cline%check                   ! checks args and prints cmdline to cmdline.txt
call cline%set('prg', 'test_polarft')
p = params(cline)                  ! constants & derived constants produced
call b%build_general_tbox(p,cline) ! general objects built
do i=1,p%nptcls
    call progress(i,p%nptcls)
    call b%img%read(p%stk, i)
    call init_polarft_tester(b%img, p%msk)
    call run_polarft_tests( lp20e, lp15e, lp10e, lp7e, snr02e, snr01e, snr005e, snr001e )
    fslp20e   = fslp20e+lp20e
    fslp15e   = fslp15e+lp15e
    fslp10e   = fslp10e+lp10e
    fslp7e    = fslp7e+lp7e
    fssnr02e  = fssnr02e+snr02e
    fssnr01e  = fssnr01e+snr01e
    fssnr005e = fssnr005e+snr005e
    fssnr001e = fssnr001e+snr001e
end do
fslp20e   = fslp20e/real(p%nptcls)
fslp15e   = fslp15e/real(p%nptcls)
fslp10e   = fslp10e/real(p%nptcls)
fslp7e    = fslp7e/real(p%nptcls)
fssnr02e  = fssnr02e/real(p%nptcls)
fssnr01e  = fssnr01e/real(p%nptcls)
fssnr005e = fssnr005e/real(p%nptcls)
fssnr001e = fssnr001e/real(p%nptcls)
write(*,'(a)') '>>> AVERAGE ANGLE ERROR'
write(*,*) 'SNR=0.2:',  fssnr02e 
write(*,*) 'SNR=0.1:',  fssnr01e 
write(*,*) 'SNR=0.05:', fssnr005e 
write(*,*) 'SNR=0.01:', fssnr001e 
write(*,*) 'LP=20:',    fslp20e 
write(*,*) 'LP=15:',    fslp15e
write(*,*) 'LP=10:',    fslp10e 
write(*,*) 'LP=7:',     fslp7e
end program 
