######################## Fortran Code ########################
./f90g95_DCUDA_local.csh ./test_code/cpu/testing_yaml.f90;
./f90g95_DCUDA_local.csh ./test_code/cpu/testing_summing.f90;
./f90g95_DCUDA_local.csh ./test_code/cpu/testing_RowColMajF90.f90;
./f90g95_DCUDA_local.csh ./test_code/cpu/testing_congruence.f90;
./f90g95_DCUDA_local.csh ./test_code/cpu/testing_BFGSNum_2D_cpu.f90;
./f90g95_DCUDA_local.csh ./test_code/cpu/testing_corr_cpu.f90;
./f90g95_DCUDA_local.csh ./test_code/cpu/testing_corr_gpu.f90;
./f90g95_DCUDA_local.csh ./test_code/cpu/testing_cshift.f90;
./f90g95_DCUDA_local.csh ./test_code/cpu/testing_NumDiff_1D_cpu.f90;
./f90g95_DCUDA_local.csh ./test_code/cpu/testing_NumGrad_2D_cpu.f90;
./f90g95_DCUDA_local.csh ./test_code/cpu/timer_tester_c.f90;
./f90g95_DCUDA_local.csh ./test_code/cpu/timer_tester.f90;

############################## C++ Code ########################
./gccg++_local.csh ./test_code/cpu/testing_RowColMajCpp.cpp;

