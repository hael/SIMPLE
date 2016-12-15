program simple_test_pftcc_inplsrch
use simple_image,            only: image
use simple_build,            only: build
use simple_params,           only: params
use simple_cmdline,          only: cmdline
use simple_simulator,        only: simimg
use simple_ori,              only: ori
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_jiffys          ! singleton
use simple_hadamard_common ! singleton
use simple_pftcc_inplsrch  ! singleton
use simple_math            ! singleton
implicit none
type(polarft_corrcalc) :: pftcc
type(ori)              :: orientation, orientation_opt
type(params)           :: p
type(build)            :: b
type(cmdline)          :: cline
integer                :: iptcl
real                   :: shlims(2,2), snr_pink, snr_detector, crxy(4)
real                   :: sherr, roerr, sherr_avg, roerr_avg
real, parameter        :: dfx=2.2, dfy=2.5, angast=30., bfac=50.
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'simple_test_pftcc_inplsrch stk=<stack.ext> msk=<mask radius(in pixels)>'
    write(*,'(a)',advance='no') ' smpd=<sampling distance(in A)> lp=<low-pass limit{20}>'
    write(*,'(a)',advance='no') ' opt=<powell|simplex|oasis|bforce|pso|de> trs=<maximum halfwidth shift(in pixels)>'
    write(*,'(a)') ' snr=<signal2noise ratio> [nrestarts=<number of opt restarts{1}>] [shbarrier=<yes|no{yes}>]'
    stop
endif
call cline%parse
call cline%checkvar('stk',  1)
call cline%checkvar('msk',  2)
call cline%checkvar('smpd', 3)
call cline%checkvar('lp',   4)
call cline%checkvar('opt',  5)
call cline%checkvar('trs',  6)
call cline%checkvar('snr',  7)
call cline%check
p = params(cline)                   ! parameters generated
call b%build_general_tbox(p, cline) ! general objects built
! build corr calculator
p%kfromto(1) = 2
p%kfromto(2) = b%img%get_find(p%lp)
call pftcc%new(1, [1,1], [p%box,p%box,1], p%kfromto, p%ring2, 'yes')
! set shift limits
shlims(:,1) = -p%trs
shlims(:,2) =  p%trs
! build searcher
call pftcc_inplsrch_init(pftcc, shlims)
call pftcc_inplsrch_set_indices(1,1)
! set snrs
snr_pink     = p%snr/0.2
snr_detector = p%snr/0.8
sherr_avg    = 0.
roerr_avg    = 0.
do iptcl=1,p%nptcls
    call progress(iptcl,p%nptcls)
    call b%img%read(p%stk, iptcl)
    call insert_ref
    call pftcc%apply_ctf(b%tfun, dfx, dfy, angast)
    call orientation%rnd_ori(p%trs)
    call orientation%set('dfx',       dfx)
    call orientation%set('dfy',       dfy)
    call orientation%set('angast', angast)
    call b%img%read(p%stk, iptcl)
    call b%img%fwd_ft
    call b%proj%rotimg(b%img, orientation%e3get(), p%msk, b%img_copy)
    b%img = b%img_copy
    call b%img%shift(orientation%get('x'),orientation%get('y'))    
    call simimg(b%img, orientation, b%tfun, 'ctf', p%snr, snr_pink, snr_detector, bfac)
    call insert_ptcl
    crxy = pftcc_inplsrch_minimize()
    call orientation_opt%set('corr', crxy(1))
    call orientation_opt%e3set(-crxy(2))
    call orientation_opt%set('x', crxy(3))
    call orientation_opt%set('y', crxy(4))
    sherr     = euclid(crxy(3:4),[orientation%get('x'),orientation%get('y')])
    sherr_avg = sherr_avg + sherr
    roerr     = rad2deg(orientation.inplrotdist.orientation_opt)
    roerr_avg = roerr_avg + roerr
end do
sherr_avg = sherr_avg/real(p%nptcls)
roerr_avg = roerr_avg/real(p%nptcls)
write(*,'(a,1x,f5.2)') 'SHIFT ERROR (IN PIXELS ): ', sherr_avg
write(*,'(a,1x,f5.2)') 'ROT   ERROR (IN DEGREES): ', roerr_avg

contains

    subroutine insert_ref
        ! normalise
        call b%img%norm
        ! apply mask
        call b%img%mask(p%msk, 'soft')
        ! move to Fourier space
        call b%img%fwd_ft
        ! transfer to polar coordinates
        call b%proj%img2polarft(1, b%img, pftcc, isptcl=.false.)
    end subroutine insert_ref

    subroutine insert_ptcl
        call prepimg4align(b, p, 1)
        call b%proj%img2polarft(1, b%img, pftcc)
    end subroutine insert_ptcl

end program simple_test_pftcc_inplsrch