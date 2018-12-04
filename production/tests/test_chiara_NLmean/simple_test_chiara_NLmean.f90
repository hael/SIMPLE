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
  if( command_argument_count() < 2 )then
      write(*,'(a)',advance='no') 'simple_test_chiara_try smpd=<sampling distance(in A)> [fname = file name]'
      stop
  endif
  call cline%parse_oldschool
  call cline%checkvar('fname', 1)
  call cline%checkvar('smpd', 2)
  !Set defaults
  call params%new(cline)  !<read cline parameters
end program simple_test_chiara_NLmean
