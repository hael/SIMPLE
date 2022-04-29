program simple_test_cartft_corrcalc
include 'simple_lib.f08'
 use simple_cartft_corrcalc, only: cartft_corrcalc
use simple_cmdline,           only: cmdline
use simple_builder,           only: builder
use simple_parameters,        only: parameters
use simple_timer
use simple_oris
use simple_image
implicit none
integer,      parameter :: BOX=256
type(parameters)        :: p
type(cartft_corrcalc)   :: cftcc
type(cmdline)           :: cline
type(builder)           :: b
type(image)             :: img
real                    :: cc

if( command_argument_count() < 3 )then
    write(logfhandle,'(a)',advance='no') 'simple_test_cartft_corrcalc lp=xx smpd=yy nthr=zz>'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('lp',   1)
call cline%checkvar('smpd', 2)
call cline%checkvar('nthr', 3)
call cline%set('ctf', 'no')
call cline%set('match_filt', 'no')

call cline%check
call p%new(cline)
p%box         = BOX
p%kfromto(1)  = 3
p%kfromto(2)  = calc_fourier_index(p%lp, p%box, p%smpd)
p%kstop       = p%kfromto(2)
p%nptcls      = 1
call b%build_general_tbox(p, cline)

call cftcc%new(1, [1, 1], .false.)

call img%new([BOX,BOX,1], p%smpd)
img = 0.0
call img%add_gauran(0.5)
call img%fft
call cftcc%set_ptcl(1,img)
img = img * 5.91
call cftcc%set_ref(1, img, .true.)

print *,'cc=',cftcc%calc_corr( 1, 1, [0.0,0.0] )

end program simple_test_cartft_corrcalc
