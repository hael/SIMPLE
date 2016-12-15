program simple_test_shsrch
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_build,            only: build
use simple_params,           only: params
use simple_rnd,              only: ran3
use simple_math,             only: euclid
use simple_cmdline,          only: cmdline
use simple_pftcc_shsrch      ! singleton
use simple_defs              ! singleton
use simple_timing            ! singleton
implicit none
type(polarft_corrcalc) :: pftcc
type(params)           :: p
type(build)            :: b
type(cmdline)          :: cline
integer                :: i, j
real                   :: cxy(3), lims(2,2), x, y, corrav, dev, ep, var, nr, avdist, nevals
integer, parameter     :: nsample = 10
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'SIMPLE_TEST_SHSRCH stk=projections.ext smpd=<saming dist(in A)> '
    write(*,'(a)',advance='no') 'msk=<mask radius(pixels)> trs=<shift halfwidth(pixels)> lp=<low-pass '
    write(*,'(a)') 'limit(in A)> [snr=<signal2noise ratio>] [opt=<optimizer(powell|simplex|oasis|pso|de)>'
    stop
endif
call cline%parse
call cline%checkvar('stk',  1)
call cline%checkvar('smpd', 2)
call cline%checkvar('msk',  3)
call cline%checkvar('trs',  4)
call cline%checkvar('lp',   5)
call cline%check
call cline%set('prg', 'test_shsrch')
p = params(cline) ! constants & derived constants produced
call b%build_general_tbox(p, cline, do3d=.false.)
p%kfromto(1) = max(2,b%img%get_find(p%hp))
p%kfromto(2) = b%img%get_find(p%lp)
! create pftcc object
call pftcc%new(1, [1,1], [p%box,p%box,1], p%kfromto, p%ring2, 'no')
! initialize search objects
lims(1,1) = -p%trs
lims(1,2) = p%trs
lims(2,1) = -p%trs
lims(2,2) = p%trs
call pftcc_shsrch_init(pftcc, lims)
call pftcc_shsrch_set_indices(1, 1, 1)
nr     = real(p%nptcls*nsample)
ep     = 0.
var    = 0.
corrav = 0.
avdist = 0.
nevals = 0
call start_Alltimers_cpu()
call start_timer_cpu("shiftsrch")
do i=1,p%nptcls
    do j=1,nsample
        ! create random shifts
        x = 2.*p%trs*ran3()-p%trs
        y = 2.*p%trs*ran3()-p%trs
        ! create shifted particle and set ref/ptcl in pftcc
        call b%img%read(p%stk,i)
        call b%img%shift(x,y,imgout=b%img_copy)
        if( cline%defined('snr') )then
            call b%img_copy%add_gauran(p%snr)
        endif
        call b%proj%img2polarft(1, b%img,      pftcc, isptcl=.false.)
        call b%proj%img2polarft(1, b%img_copy, pftcc, isptcl=.true.)
        ! minimize
        cxy = pftcc_shsrch_minimize()
        ! calc stats
        dev    = euclid([x,y],[cxy(2),cxy(3)])
        avdist = avdist+dev
        corrav = corrav+cxy(1)
        nevals = nevals+real(pftcc_shsrch_get_nevals())
!         write(*,*) 'correct shift: ', x, y, 'shift identified by ', trim(p%opt)//' ', cxy(2), cxy(3), 'correlation: ', cxy(1)
    end do
end do
call stop_timer_cpu("shiftsrch")
call stop_Alltimers_cpu()
avdist = avdist/nr
corrav = corrav/real(p%nptcls*nsample)
nevals = nevals/nr
write(*,*) 'average correlation: ', corrav
write(*,*) 'average error: ', avdist
write(*,*) 'nevals: ', int(nevals)
end program simple_test_shsrch