!USAGE : simple_test_chiara_NLmean prg=dummy smpd=1. fname='/home/chiara/Desktop/Chiara/mrcsFePtNanoparticles/BiMetal6_dm4_Hour_00_Minute_00_Second_00_Frame_0000.mrc'
program simple_test_chiara_NLmean
  include 'simple_lib.f08'
  use simple_image, only : image
  use simple_parameters, only: parameters
  use simple_cmdline,    only: cmdline
  implicit none
  type(image)        :: img
  integer            :: ldim(3), nptcls
  type(cmdline)      :: cline
  type(parameters)   :: params
  real               :: smpd
  if( command_argument_count() < 2 )then
      write(*,'(a)',advance='no') 'simple_test_chiara_NLmean smpd=<sampling distance(in A)> [fname = file name]'
      stop
  endif
  call cline%parse_oldschool
  call cline%checkvar('fname', 1)
  call cline%checkvar('smpd', 2)
  !Set defaults
  call params%new(cline)  !<read cline parameters
  call find_ldim_nptcls(params%fname, ldim, nptcls, smpd)
  call img%new(ldim, smpd)
  call img%read(params%fname)
  call img%NLmean()
  call img%write('NLmeanDenoised.mrc')
end program simple_test_chiara_NLmean
