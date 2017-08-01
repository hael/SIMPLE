program testing_gen_polar_coords

  use simple_defs
  use simple_cuda_defs
  use simple_cuda
  use simple_math_gpu
  use simple_timing
  use greeting_version

  implicit none
#define devptr_t integer*8

  integer                          :: ring2

  devptr_t                        :: devPtr_kfromto
  devptr_t                        :: devPtr_coords
  devptr_t                        :: devptr_angtab

  !start of the execution commands
  !start of the greeting message
  call hello_gpu_math()
  call timestamp()
  call start_Alltimers_cpu()

  ring2 = 2
  !temporary test values
  devPtr_kfromto = 5.0
  devPtr_coords = 6.0
  devPtr_angtab = 7.0
  
  write(*,*) 
  write(*,*) "in the testing_gen_polar_coords"
  write(*,*) ring2, devPtr_kfromto, devPtr_coords, devPtr_angtab
  write(*,*) 


  call gen_polar_coords_cuda_gpu( devPtr_kfromto, ring2, &
				  devPtr_coords, devPtr_angtab )

  !shutting down the environment
  call simple_cuda_shutdown()

  !shutting down the timers
  call stop_Alltimers_cpu()

  !end of greeting message
  call bye_gpu_math()



end program testing_gen_polar_coords
