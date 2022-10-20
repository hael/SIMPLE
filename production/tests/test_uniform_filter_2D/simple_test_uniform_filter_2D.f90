program simple_test_uniform_filter_2D
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_commander_volops,   only: reproject_commander
use simple_image,              only: image
use simple_opt_filter,         only: uniform_filter_2D
implicit none
type(parameters)              :: p
type(cmdline)                 :: cline
type(image)                   :: img_noisy, img_filt, weights_img
integer                       :: nptcls, iptcl, m, n
real                          :: best_ind, val
real,    parameter            :: LP_LOWRES_PHASE = 7.
integer, parameter            :: NSEARCH = 100, SMOOTH_EXT = 8
real,    allocatable          :: butterworth_fil(:)
if( command_argument_count() < 5 )then
    write(logfhandle,'(a)') 'Usage: simple_test_uniform_filter_2D smpd=xx nthr=yy stk=stk.mrc mskdiam=zz lp=ll'
    write(logfhandle,'(a)') 'Example: projections of https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
    stop
else
    call cline%parse_oldschool
endif
call cline%checkvar('smpd',    1)
call cline%checkvar('nthr',    2)
call cline%checkvar('stk' ,    3)
call cline%checkvar('mskdiam', 4)
call cline%checkvar('lp',      5)
call cline%check
call p%new(cline)
call find_ldim_nptcls(p%stk, p%ldim, nptcls)
p%ldim(3) = 1 ! because we operate on stacks
!write(*, *) 'Filtering in progress...'
call img_noisy%new(p%ldim, p%smpd)
call img_filt %new(p%ldim, p%smpd)
call weights_img%new(p%ldim, p%smpd, .false.)
call weights_img%zero_and_unflag_ft()
do m = -SMOOTH_EXT, SMOOTH_EXT
    do n = -SMOOTH_EXT, SMOOTH_EXT
        val = -hyp(real(m), real(n))/(SMOOTH_EXT + 1) + 1.
        if( val > 0 ) call weights_img%set_rmat_at(p%ldim(1)/2+m+1, p%ldim(1)/2+n+1, 1, val)
    enddo
enddo
call weights_img%fft()
allocate(butterworth_fil(p%ldim(1)), source=0.)
do iptcl = 1, min(10, p%nptcls)
    !write(*, *) 'Particle # ', iptcl
    ! comparing the nonuniform result with the original data
    call img_noisy%read(p%stk, iptcl)
    best_ind = uniform_filter_2D(img_noisy, img_noisy, img_filt, weights_img, butterworth_fil, p%lp, NSEARCH, SMOOTH_EXT)
    print *, calc_lowpass_lim(nint(best_ind), p%ldim(1), p%smpd)
    call img_filt%write('stk_img_filt.mrc', iptcl)
    ! spherical masking
    call img_filt%mask(p%msk, 'soft')
    call img_noisy%zero_and_unflag_ft
    call img_filt %zero_and_unflag_ft
enddo
call img_noisy%kill()
call img_filt %kill()
call weights_img%kill()
end program simple_test_uniform_filter_2D