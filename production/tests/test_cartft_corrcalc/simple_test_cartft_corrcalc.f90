program simple_test_cartft_corrcalc
include 'simple_lib.f08'
use simple_cartft_corrcalc, only: cartft_corrcalc
use simple_cmdline,         only: cmdline
use simple_builder,         only: builder
use simple_parameters,      only: parameters
use simple_timer
use simple_oris
use simple_image
implicit none
type(parameters)      :: p
type(cartft_corrcalc) :: cftcc
type(cmdline)         :: cline
type(builder)         :: b
type(image)           :: img
real                  :: cc, shifts(2,2)
integer               :: iptcl, iref, loc(1)
real,    allocatable  :: ccmat(:,:), pshifts(:,:)
integer, allocatable  :: passign(:), ground_truth(:)
if( command_argument_count() < 3 )then
    write(logfhandle,'(a)') 'simple_test_cartft_corrcalc lp=xx smpd=yy nthr=zz stk=reprojs>'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('lp',   1)
call cline%checkvar('smpd', 2)
call cline%checkvar('nthr', 3)
call cline%checkvar('stk',  4)
call cline%set('ctf', 'no')
call cline%set('match_filt', 'no')
call cline%check
call p%new(cline)
p%kfromto(1)  = 3
p%kfromto(2)  = calc_fourier_index(p%lp, p%box, p%smpd)
p%kstop       = p%kfromto(2)
call b%build_general_tbox(p, cline)
call cftcc%new(1, [1, 1], .false.)
call img%new([p%box,p%box,1], p%smpd)
img = 0.0
call img%add_gauran(0.5)
call img%fft
call cftcc%set_ptcl(1,img)
img = img * 5.91
call cftcc%set_ref(1, img, .true.)

print *,'cc=',cftcc%calc_corr( 1, 1, [0.0,0.0] )

call cftcc%new(p%nptcls, [1, p%nptcls], .false.)
do iref = 1, p%nptcls
    call img%read(p%stk, iref)
    call img%fft
    call cftcc%set_ref(iref, img, .true.)
end do
do iptcl = 1,p%nptcls
    call img%read(p%stk, iptcl)
    call img%fft
    call cftcc%set_ptcl(iptcl,img)
end do
allocate(ccmat(p%nptcls,p%nptcls), source=-1.)
do iref = 1, p%nptcls
    do iptcl = 1, p%nptcls
        ccmat(iref,iptcl) = cftcc%calc_corr(iref, iptcl, [0.0,0.0])
    end do
end do
allocate(passign(p%nptcls), ground_truth(p%nptcls), source=0)
do iref = 1, p%nptcls
    loc = maxloc(ccmat(iref,:))
    passign(iref) = loc(1)
end do
ground_truth = (/(iptcl,iptcl=1,p%nptcls)/)
if( all(passign == ground_truth) )then
    print *, '*** assignment test PASSED ***'
else
    print *, '*** assignment test FAILED ***'
    stop
endif

! shift array
shifts(1,1) = -0.5
shifts(1,2) =  0.5
shifts(2,1) = -0.5
shifts(2,2) =  0.5

end program simple_test_cartft_corrcalc
