program test_dock_vols
! include 'simple_lib.f08'
! use simple_dock_vols,       only: dock_vols
! use simple_parameters,      only: parameters
! use simple_cmdline,         only: cmdline
! use simple_image,           only: image
! use simple_simple_volinterp,  only: rotvol
! use simple_ori,             only: ori
! implicit none

! character(len=*), parameter :: VREF = 'vol1.mrc', VTARG = 'vol2.mrc', VTARG_DOCKED = 'vol2_docked.mrc'
! type(dock_vols)  :: dvols
! type(parameters) :: params
! type(cmdline)    :: cline
! type(image)      :: rotvol_shifted, vol_shifted
! type(ori)        :: o
! real    :: euls(3), shift(3)
! integer :: ldim(3),nptcls

! if( command_argument_count() < 4 )then
!     write(logfhandle,'(a)',advance='no') 'ERROR! Usage: simple_test_dov_vols vol1=<volume.ext> mskdiam=<mask radius(in A)>'
!     write(logfhandle,'(a)') ' smpd=<sampling distance(in A)> nthr=<number of threads{1}> [hp=<freq in A>] [lp=<freq in A>]'
! endif
! call cline%parse_oldschool
! call cline%set('prg',      'xxx')
! call cline%checkvar('vol1',     1)
! call cline%checkvar('mskdiam',  2)
! call cline%checkvar('smpd',     3)
! call cline%checkvar('nthr',     4)
! ! call cline%checkvar('hp',       5) ! optional
! ! call cline%checkvar('lp',       6) ! optional
! call cline%check
! call params%new(cline)
! ! random shift & orientation
! euls  = [360.*ran3(), 180.*ran3(), 360.*ran3()]
! shift = 10.*[ran3()-0.5, ran3()-0.5, ran3()-0.5]
! ! some singular values
! ! euls  = [ 52.3584137,      53.7275314,      325.139801]
! ! shift = [ -4.54392052,     -1.00544930,      -3.04297638]
! ! euls  =[   217.004868,       106.567047,       335.972809]
! ! shift =[  -2.45968294,      0.471755266,      -6.28441572E-02]
! ! euls    =[ 29.5008869,       99.0748367,       194.840454]
! ! shift   =[3.88919878,       4.76009178,       3.89923334]
! o = ori(is_ptcl=.true.)
! call o%set_euler(euls)
! ! shifted rotated volume
! call find_ldim_nptcls(params%vols(1), ldim, nptcls)
! call vol_shifted%new(ldim,params%smpd)
! call vol_shifted%read(params%vols(1))
! call vol_shifted%write(VREF)
! call vol_shifted%fft
! call vol_shifted%shift(shift)
! call vol_shifted%ifft
! rotvol_shifted = rotvol(vol_shifted, o)
! call rotvol_shifted%write(VTARG)
! ! search
! call dvols%new(VREF, VTARG, params%smpd, params%hp, params%lp, params%mskdiam, mag=.true.)
! call dvols%srch
! ! target docked to reference (vol1)
! call dvols%rotate_target(VTARG, VTARG_DOCKED)
end program test_dock_vols
