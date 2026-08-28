!@descr: A/B the CUDA-C flex insertion kernel against the CPU batch path (P1 gate)
program simple_test_flex_gpu
use simple_flex_gpu, only: test_flex_gpu_insert, test_flex_gpu_coupled, test_flex_gpu_coupled_banked, &
    &test_flex_gpu_psample, test_flex_gpu_estep
implicit none
call test_flex_gpu_insert
call test_flex_gpu_coupled
call test_flex_gpu_coupled_banked
call test_flex_gpu_psample
call test_flex_gpu_estep
write(*,'(A)') '**** SIMPLE_TEST_FLEX_GPU NORMAL STOP ****'
end program simple_test_flex_gpu
