!==Program simple_fsc_filt
!
! <fsc_filt/begin> is a program for calculating the Fourier Shell Correlation (FSC) between inputted volumes 
! \texttt{vol2} and \texttt{vol3} and apply an optimal low-pass filter $\frac{2FSC}{FSC+1}$  to the 
! \texttt{vol1} input volume. The \texttt{msk} parameter controls the radius (in pixels) of the soft-edged 
! mask applied to the even/odd volumes. <fsc_filt/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_fsc_filt
use simple_cmdline   ! singleton
use simple_math      ! singleton
use simple_build,    only: build
use simple_params,   only: params
use simple_image,    only: image
use simple_jiffys,   only: simple_end
use simple_filterer, only: apply_optlp, optlp2file
implicit none
type(params)      :: p
type(build)       :: b
type(image)       :: vol2
real, allocatable :: res(:), corrs(:)
if( command_argument_count() < 6 )then
    write(*,'(a)', advance='no') 'SIMPLE_FSC_FILT vol1=<raw.ext> vol2=<even.ext> vol3=<odd.ext> '
    write(*,'(a)', advance='no') ' smpd=<sampling distance(in A)> msk=<mask radius(in pixels)>'
    write(*,'(a)', advance='no') ' outvol=<filtered.ext> [norec=<yes|no{no}>] [mw=<molecular weight'
    write(*,'(a)', advance='no') ' (in kD)>] [nthr=<nr of openMP threads{1}>] [box=<image size'
    write(*,'(a)', advance='no') '(in pixels)>] [inner=<inner mask radius(in pixels)>]'
    write(*,'(a)', advance='no') ' [width=<pixels falloff inner mask{10}>] [lpstop=<stay at this'
    write(*,'(a)') ' low-pass limit(in A)>] [dens=<density(e.g. 9.368 Da/A3 4 gold clusters){0.}>]'
    stop
endif
call parse_cmdline
call cmdcheckvar('vol1',   1)
call cmdcheckvar('vol2',   2)
call cmdcheckvar('vol3',   3)
call cmdcheckvar('smpd',   4)
call cmdcheckvar('msk',    5)
call cmdcheckvar('outvol', 6)
call cmdcheck
p = params()                                         ! constants & derived constants produced, mode=2
call b%build_general_tbox(p)                         ! general objects built 
call b%build_hadamard_prime_tbox(p)                  ! refinement objects built
call b%vol%read(p%vols(2))
call b%vol%mask(p%msk, 'soft')
call b%vol%fwd_ft
call vol2%copy(b%vol)
call vol2%read(p%vols(3))
call vol2%mask(p%msk, 'soft')
call vol2%fwd_ft
call b%vol%fsc(vol2, res, corrs)
call optlp2file(b%vol, vol2, 1, p%fny)          ! OBS: state is 1
call b%vol%read(p%vols(1))
call apply_optlp(b%vol, 1)                          ! OBS: state is 1
call b%vol%write(p%outvol, del_if_exists=.true.)
call simple_end('**** SIMPLE_FSC_FILT NORMAL STOP ****')
end program simple_fsc_filt