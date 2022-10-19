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
type(image)                   :: img_noisy, img_filt
integer                       :: nptcls, iptcl
real                          :: best_ind
real,    parameter            :: LP_LOWRES_PHASE = 7.
integer, parameter            :: NSEARCH = 100
if( command_argument_count() < 6 )then
    write(logfhandle,'(a)') 'Usage: simple_test_uniform_filter_2D smpd=xx nthr=yy stk=stk.mrc mskdiam=zz lp=ll hp=hh'
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
call cline%checkvar('hp',      6)
call cline%check
call p%new(cline)
call find_ldim_nptcls(p%stk, p%ldim, nptcls)
p%ldim(3) = 1 ! because we operate on stacks
!write(*, *) 'Filtering in progress...'
call img_noisy%new(p%ldim, p%smpd)
call img_filt %new(p%ldim, p%smpd)
do iptcl = 1, p%nptcls
    !write(*, *) 'Particle # ', iptcl
    ! comparing the nonuniform result with the original data
    call img_noisy%read(p%stk, iptcl)
    best_ind = uniform_filter_2D(img_noisy, img_noisy, img_filt, p%lp, p%hp, NSEARCH)
    print *, calc_lowpass_lim(nint(best_ind), p%ldim(1), p%smpd)
    call img_filt%write('stk_img_filt.mrc', iptcl)
    ! spherical masking
    call img_filt%mask(p%msk, 'soft')
    call img_noisy%zero_and_unflag_ft
    call img_filt %zero_and_unflag_ft
enddo
call img_noisy%kill()
call img_filt %kill()
end program simple_test_uniform_filter_2D