program simple_test_wiener
use simple_cmdline  ! singleton
use simple_defs     ! singleton
use simple_image,   only: image
use simple_build,   only: build
use simple_params,  only: params
use simple_jiffys,  only: simple_end, progress, file2rarr, file_exists
use simple_math,    only: fsc2ssnr
implicit none
type(build)       :: b
type(params)      :: p
type(image)       :: filter
integer           :: i 
real, allocatable :: ssnr(:)
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_TEST_WIENER stk=<ptclstk.ext> smpd=<sampling dist(in A)>'
    write(*,'(a)',advance='no') ' oritab=<ori/ctf-params.txt> [fsc=<fsc_state1.bin>]'
    write(*,'(a)') ' [ssnr=<ssnr_state1.bin>] [ctf=<yes|flip|mul|{flip}>] [outstk=<outstk.ext>]'
    stop
endif
call parse_cmdline
call cmdcheckvar('stk',    1)
call cmdcheckvar('smpd',   2)
call cmdcheckvar('oritab', 3)
if( .not. defined_cmd_arg('ctf') )then
    call set_cmdline('ctf', 'flip')
endif
call cmdcheck
p = params() ! parameters generated
call b%build_general_tbox(p,do3d=.false.)
if( file_exists(p%ssnr) )then
    ssnr = file2rarr(p%ssnr)
else if( file_exists(p%fsc) )then
    ssnr = fsc2ssnr(file2rarr(p%fsc))
    stop 'need either SSNR file (ssnr_state1.bin) or FSC file (fsc_state1.bin) to run simple_test_wiener'
endif
do i=1,p%nptcls
    call b%img%read(p%stk, i)
    call b%tfun%ctf2img(b%img_copy, b%a%get(i,'dfx'), 'ctf', b%a%get(i,'dfy'), b%a%get(i,'angast'), 100.)
    call b%img_copy%ft2img('real', b%img_copy)
    call b%img_copy%write('ctfimgs'//p%ext, i)
    filter = b%tfun%wiener(b%img%get_ldim(), p%smpd, ssnr, p%ctf, b%a%get(i,'dfx'), b%a%get(i,'dfy'), b%a%get(i,'angast'))
    call filter%ft2img('real', b%wiener)
    call b%wiener%write('wienerimgs'//p%ext, i)
    call b%img%apply_filter(filter)
    call b%img%write(p%outstk, i)
    call progress(i, p%nptcls)
end do
call simple_end('**** SIMPLE_TEST_WIENER NORMAL STOP ****')
end program