#!/bin/bash
set TARGET_HOME_SIDE = /home/frederic/SimpleSideBrch/19Jul16/Simple_Restruct.projet

set DIFF_PATH = $TARGET_HOME_SIDE/simple_utils/common

echo "::::::::::::::"
echo "diff_common.csh"
echo "::::::::::::::"
diff diff_common.csh $DIFF_PATH
echo "::::::::::::::"
echo "diff_test_common.csh"
echo "::::::::::::::"
diff diff_test_common.csh $DIFF_PATH
echo "::::::::::::::"
echo "fft123D_cpu.f90"
echo "::::::::::::::"
diff fft123D_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "fft123D_gpu.f90"
echo "::::::::::::::"
diff fft123D_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "get_cpu_time_c.c"
echo "::::::::::::::"
diff get_cpu_time_c.c $DIFF_PATH
echo "::::::::::::::"
echo "get_deviceQuery_gpu_c.c"
echo "::::::::::::::"
diff get_deviceQuery_gpu_c.c $DIFF_PATH
echo "::::::::::::::"
echo "get_deviceQuery_gpu.cpp"
echo "::::::::::::::"
diff get_deviceQuery_gpu.cpp $DIFF_PATH
echo "::::::::::::::"
echo "get_fft123D_cpu_c.c"
echo "::::::::::::::"
diff get_fft123D_cpu_c.c $DIFF_PATH
echo "::::::::::::::"
echo "get_fft123D_cpu.cpp"
echo "::::::::::::::"
diff get_fft123D_cpu.cpp $DIFF_PATH
echo "::::::::::::::"
echo "get_fft123D_gpu_c.c"
echo "::::::::::::::"
diff get_fft123D_gpu_c.c $DIFF_PATH
echo "::::::::::::::"
echo "get_fft123D_gpu.cpp"
echo "::::::::::::::"
diff get_fft123D_gpu.cpp $DIFF_PATH
echo "::::::::::::::"
echo "get_systemQuery_cpu_c.c"
echo "::::::::::::::"
diff get_systemQuery_cpu_c.c $DIFF_PATH
echo "::::::::::::::"
echo "get_systemQuery_cpu.cpp"
echo "::::::::::::::"
diff get_systemQuery_cpu.cpp $DIFF_PATH
echo "::::::::::::::"
echo "greeting_version.f90"
echo "::::::::::::::"
diff greeting_version.f90 $DIFF_PATH
echo "::::::::::::::"
echo "initialu.f90"
echo "::::::::::::::"
diff initialu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "Makefile_target"
echo "::::::::::::::"
diff Makefile_target $DIFF_PATH
echo "::::::::::::::"
echo "matrixGetter.f90"
echo "::::::::::::::"
diff matrixGetter.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_DataParallel_IO.f90"
echo "::::::::::::::"
diff simple_DataParallel_IO.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_deviceQuery_gpu.f90"
echo "::::::::::::::"
diff simple_deviceQuery_gpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_dynamic_memory.f90"
echo "::::::::::::::"
diff simple_dynamic_memory.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_eglossary.f90"
echo "::::::::::::::"
diff simple_eglossary.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_eglossary_lowlev.f90"
echo "::::::::::::::"
diff simple_eglossary_lowlev.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_error_handling.f90"
echo "::::::::::::::"
diff simple_error_handling.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_error_handling-incFl.f90"
echo "::::::::::::::"
diff simple_error_handling-incFl.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_file_highlev.f90"
echo "::::::::::::::"
diff simple_file_highlev.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_file_utils.f90"
echo "::::::::::::::"
diff simple_file_utils.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_memory_profiling.f90"
echo "::::::::::::::"
diff simple_memory_profiling.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_module_file_malloc.f90"
echo "::::::::::::::"
diff simple_module_file_malloc.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_random.f90"
echo "::::::::::::::"
diff simple_random.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_sorting.f90"
echo "::::::::::::::"
diff simple_sorting.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_SU3_tester.f90"
echo "::::::::::::::"
diff simple_SU3_tester.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_systemQuery_cpu.f90"
echo "::::::::::::::"
diff simple_systemQuery_cpu.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_testfuns.f90"
echo "::::::::::::::"
diff simple_testfuns.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_textHandler.f90"
echo "::::::::::::::"
diff simple_textHandler.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_timing.f90"
echo "::::::::::::::"
diff simple_timing.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_yaml_open-incFl.f90"
echo "::::::::::::::"
diff simple_yaml_open-incFl.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_yaml_output.f90"
echo "::::::::::::::"
diff simple_yaml_output.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_yaml_strings.f90"
echo "::::::::::::::"
diff simple_yaml_strings.f90 $DIFF_PATH
echo "::::::::::::::"
echo "strlen.f90"
echo "::::::::::::::"
diff strlen.f90 $DIFF_PATH
echo "::::::::::::::"
echo "su2random.f90"
echo "::::::::::::::"
diff su2random.f90 $DIFF_PATH
echo "::::::::::::::"
echo "su3random.f90"
echo "::::::::::::::"
diff su3random.f90 $DIFF_PATH
echo "::::::::::::::"
echo "timestamp.f90"
echo "::::::::::::::"
diff timestamp.f90 $DIFF_PATH
echo "::::::::::::::"
echo "timming_c.c"
echo "::::::::::::::"
diff timming_c.c $DIFF_PATH
echo "::::::::::::::"
echo "timming.cpp"
echo "::::::::::::::"
diff timming.cpp $DIFF_PATH
echo "::::::::::::::"
echo "simple_yaml_toa-incFl.f90"
echo "::::::::::::::"
diff simple_yaml_toa-incFl.f90 $DIFF_PATH
