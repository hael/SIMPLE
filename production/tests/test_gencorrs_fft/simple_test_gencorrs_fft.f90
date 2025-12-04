program simple_test_gen_corrs_fft
include 'simple_lib.f08'
 use simple_polarft_calc, only: polarft_calc
use simple_cmdline,           only: cmdline
use simple_builder,           only: builder
use simple_parameters,        only: parameters
use simple_timer
implicit none
type(parameters)        :: p
type(polarft_calc)  :: pftc
type(cmdline)           :: cline
type(builder)           :: b
real, allocatable       :: cc(:), cc_fft(:)
integer                 :: iptcl, jptcl
integer(timer_int_kind) :: tfft
if( command_argument_count() < 3 )then
    write(logfhandle,'(a)',advance='no') 'simple_test_srch stk=<particles.mrc> msk=<mask radius(in pixels)>'
    write(logfhandle,'(a)') ' smpd=<sampling distance(in A)>'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('stk',  1)
call cline%checkvar('msk',  2)
call cline%checkvar('smpd', 3)
call cline%check
call p%new(cline)
p%kfromto(1) = 2
p%kfromto(2) = 100
call b%build_general_tbox(p, cline)
call pftc%new(p%nptcls, [1, p%nptcls], p%kfromto)
call b%img_crop_polarizer%init_polarizer(pftc, p%alpha)
do iptcl=1,p%nptcls
    call b%img_crop_polarizer%read(p%stk, iptcl)
    call b%img_crop_polarizer%fft()
    ! transfer to polar coordinates
    call b%img_crop_polarizer%polarize(pftc, iptcl, isptcl=.false., iseven=.true., mask=b%l_resmsk)
    call b%img_crop_polarizer%polarize(pftc, iptcl, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
end do
allocate(cc(pftc%get_nrots()), cc_fft(pftc%get_nrots()))

!### TIMING

tfft = tic()
do iptcl=1,p%nptcls - 1
    do jptcl=iptcl + 1, p%nptcls
        call pftc%gen_objfun_vals(iptcl, jptcl, [0.,0.], cc_fft)
    end do
end do
print *, 'time of fft_mod: ', toc(tfft)
end program simple_test_gen_corrs_fft
