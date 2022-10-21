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
type(image)                   :: img_clean, img_noisy, img_filt, noise_img, weights_img
integer                       :: nptcls, iptcl, m, n, find_start, find_stop, noise_n, noise_i
real                          :: best_ind, val, ave, sdev, maxv, minv, noise_lvl
real,    parameter            :: LP_LOWRES_PHASE = 7., NOISE_MIN = 0.1, NOISE_MAX = 1., NOISE_DEL = 0.1
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
call img_clean%new(p%ldim, p%smpd)
call img_filt %new(p%ldim, p%smpd)
call noise_img%new(p%ldim, p%smpd)
call weights_img%new(p%ldim, p%smpd, .false.)
call weights_img%zero_and_unflag_ft()
do m = -SMOOTH_EXT, SMOOTH_EXT
    do n = -SMOOTH_EXT, SMOOTH_EXT
        val = -hyp(real(m), real(n))/(SMOOTH_EXT + 1) + 1.
        if( val > 0 ) call weights_img%set_rmat_at(p%ldim(1)/2+m+1, p%ldim(1)/2+n+1, 1, val)
    enddo
enddo
call weights_img%fft()
find_start = calc_fourier_index(p%lp,     p%ldim(1), p%smpd)
find_stop  = calc_fourier_index(p%smpd*2, p%ldim(1), p%smpd)
noise_n    = int((NOISE_MAX - NOISE_MIN)/NOISE_DEL + 1.)
allocate(butterworth_fil(p%ldim(1)), source=0.)
do iptcl = 1, min(10, p%nptcls)
    write(*, *) 'Particle # ', iptcl
    ! comparing the nonuniform result with the original data
    do noise_i = 1, noise_n
        noise_lvl = NOISE_MIN + (noise_i - 1)*NOISE_DEL
        print *, 'noise level = ', noise_lvl
        call img_clean%read(p%stk, iptcl)
        call img_clean%mask(p%msk, 'soft')
        call img_clean%stats('foreground', ave, sdev, maxv, minv)
        ! adding noise
        call noise_img%gauran(0., noise_lvl * sdev)
        call noise_img%mask(1.5 * p%msk, 'soft')
        call img_noisy%copy(img_clean)
        call img_noisy%add(noise_img)
        call img_noisy%write('stk_img_noisy.mrc', iptcl)
        best_ind = uniform_filter_2D(img_noisy, img_clean, img_filt, weights_img, butterworth_fil, find_start, find_stop, NSEARCH, SMOOTH_EXT)
        print *, 'best resolution cut_off = ', calc_lowpass_lim(nint(best_ind), p%ldim(1), p%smpd)
        call img_noisy%write('stk_img_filt.mrc', iptcl)
        ! spherical masking
        call img_noisy%zero_and_unflag_ft
        call img_clean%zero_and_unflag_ft
        call img_filt %zero_and_unflag_ft
    enddo
enddo
call img_noisy%kill()
call img_clean%kill()
call img_filt %kill()
call noise_img%kill()
call weights_img%kill()
end program simple_test_uniform_filter_2D