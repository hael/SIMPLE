program simple_test_shiftsrch
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_cmdline,          only: cmdline
use simple_build,            only: build
use simple_image,            only: image
use simple_params,           only: params
use simple_polarizer,        only: polarizer
implicit none
type(cmdline)          :: cline
type(build)            :: b
type(params)           :: p
type(polarft_corrcalc) :: pftcc
type(polarizer)        :: img_copy
logical                :: be_verbose=.false.
real, parameter        :: SHMAG=1.0
real, allocatable      :: corrs(:)
real                   :: corrmax, corr
integer                :: xsh, ysh, xbest, ybest, i
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'simple_test_shiftsrch stk=<particles.ext> msk=<mask radius(in pixels)>'
    write(*,'(a)') ' smpd=<sampling distance(in A)> [nthr=<number of threads{1}>] [verbose=<yes|no{no}>]'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('stk',  1)
call cline%checkvar('msk',  2)
call cline%checkvar('smpd', 3)
call cline%check
be_verbose = .false.
if( cline%defined('verbose') )then
    if( trim(cline%get_carg('verbose')) .eq. 'yes' )then
        be_verbose = .true.
    endif
endif
p = params(cline)
p%kfromto(1) = 2
p%kfromto(2) = 40
call b%build_general_tbox(p, cline)
call pftcc%new(8, p)
allocate(corrs(pftcc%get_nrots()))
call img_copy%init_polarizer(pftcc, p%alpha)
call b%img%read(p%stk, 1)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 1, isptcl=.false., iseven=.true.)
call img_copy%shift([SHMAG,0.,0.])   ! left
call img_copy%polarize(pftcc, 1, isptcl=.true., iseven=.true.)
call img_copy%ifft()
call img_copy%write('shifted.mrc', 1)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 2, isptcl=.false., iseven=.true.)
call img_copy%shift([0.,SHMAG,0.])   ! down
call img_copy%polarize(pftcc, 2, isptcl=.true., iseven=.true.)
call img_copy%ifft()
call img_copy%write('shifted.mrc', 2)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 3, isptcl=.false., iseven=.true.)
call img_copy%shift([-SHMAG,0.,0.])   ! right
call img_copy%polarize(pftcc, 3, isptcl=.true., iseven=.true.)
call img_copy%ifft()
call img_copy%write('shifted.mrc', 3)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 4, isptcl=.false., iseven=.true.)
call img_copy%shift([0.,-SHMAG,0.])   ! up
call img_copy%polarize(pftcc, 4, isptcl=.true., iseven=.true.)
call img_copy%ifft()
call img_copy%write('shifted.mrc', 4)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 5, isptcl=.false., iseven=.true.)
call img_copy%shift([SHMAG,SHMAG,0.])   ! left + down
call img_copy%polarize(pftcc, 5, isptcl=.true., iseven=.true.)
call img_copy%ifft()
call img_copy%write('shifted.mrc', 5)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 6, isptcl=.false., iseven=.true.)
call img_copy%shift([-SHMAG,-SHMAG,0.]) ! right + up
call img_copy%polarize(pftcc, 6, isptcl=.true., iseven=.true.)
call img_copy%ifft()
call img_copy%write('shifted.mrc', 6)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 7, isptcl=.false., iseven=.true.)
call img_copy%shift([-SHMAG,SHMAG,0.])  ! right + down
call img_copy%polarize(pftcc, 7, isptcl=.true., iseven=.true.)
call img_copy%ifft()
call img_copy%write('shifted.mrc', 7)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 8, isptcl=.false., iseven=.true.)
call img_copy%shift([SHMAG,-SHMAG,0.])  ! left + up
call img_copy%polarize(pftcc, 8, isptcl=.true., iseven=.true.)
call img_copy%ifft()
call img_copy%write('shifted.mrc', 8)
call pftcc%memoize_ffts
do i=1,8
    corrmax = -1.
    do xsh=-2,2
        do ysh=-2,2
            call pftcc%gencorrs(i, i, real([xsh,ysh]), corrs)
            corr  = maxval(corrs)

            print *, 'corr: ', corr, xsh, ysh

            if( corr > corrmax )then
                corrmax = corr
                xbest   = xsh
                ybest   = ysh
            endif
        end do
    end do
    print *, xbest, ybest
end do
end program simple_test_shiftsrch
