################################################################
# Set CUDA
################################################################
find_package(CUDA QUIET REQUIRED)
set(CUDA_USE_STATIC_CUDA_RUNTIME OFF)
if (CUDA_FOUND)
  try_run(RUN_RESULT_VAR COMPILE_RESULT_VAR
    ${CMAKE_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/has_cuda_gpu.c
    CMAKE_FLAGS
    -DINCLUDE_DIRECTORIES:STRING=${CUDA_TOOLKIT_INCLUDE}
    -DLINK_LIBRARIES:STRING=${CUDA_CUDART_LIBRARY}
    COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT_VAR
    RUN_OUTPUT_VARIABLE RUN_OUTPUT_VAR)
  message(STATUS "${RUN_OUTPUT_VAR}") # Display number of GPUs found
  # COMPILE_RESULT_VAR is TRUE when compile succeeds
  # RUN_RESULT_VAR is zero when a GPU is found
  if(COMPILE_RESULT_VAR AND NOT RUN_RESULT_VAR)
    set(CUDA_HAVE_GPU TRUE CACHE BOOL "Whether CUDA-capable GPU is present")
  else()
    set(CUDA_HAVE_GPU FALSE CACHE BOOL "Whether CUDA-capable GPU is present")
  endif()

  find_program(_nvidia_smi "nvidia-smi")
  if (_nvidia_smi)
    set(DETECT_GPU_COUNT_NVIDIA_SMI 0)
    # execute nvidia-smi -L to get a short list of GPUs available
    exec_program(${_nvidia_smi_path} ARGS -L
      OUTPUT_VARIABLE _nvidia_smi_out
      RETURN_VALUE    _nvidia_smi_ret)
    # process the stdout of nvidia-smi
    if (_nvidia_smi_ret EQUAL 0)
      # convert string with newlines to list of strings
      string(REGEX REPLACE "\n" ";" _nvidia_smi_out "${_nvidia_smi_out}")
      foreach(_line ${_nvidia_smi_out})
        if (_line MATCHES "^GPU [0-9]+:")
          math(EXPR DETECT_GPU_COUNT_NVIDIA_SMI "${DETECT_GPU_COUNT_NVIDIA_SMI}+1")
          # the UUID is not very useful for the user, remove it
          string(REGEX REPLACE " \\(UUID:.*\\)" "" _gpu_info "${_line}")
          if (NOT _gpu_info STREQUAL "")
            list(APPEND DETECT_GPU_INFO "${_gpu_info}")
          endif()
        endif()
      endforeach()
      check_num_gpu_info(${DETECT_GPU_COUNT_NVIDIA_SMI} DETECT_GPU_INFO)
      set(DETECT_GPU_COUNT ${DETECT_GPU_COUNT_NVIDIA_SMI})
    endif()
  endif()
  if (APPLE AND NOT CUDA_HOST_COMPILER)
    # https://gitlab.kitware.com/cmake/cmake/issues/13674
    message(WARNING "CUDA_HOST_COMPILER is set to the default /usr/bin/clang; CUDA 8.0 requires 7.3 or older.")
    set(CUDA_HOST_COMPILER "/usr/bin/clang" CACHE FILEPATH "" FORCE)
  endif()
if (NOT DEFINED CUDA_ARCH)
   set(CUDA_ARCH "61")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wall -Werror -DCUDA_ARCH=${CUDA_ARCH} -std=c++11 ${OpenMP_CXX_FLAGS}")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(NVCC_FLAGS "-G -g")
endif()
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch sm_${CUDA_ARCH} -Xptxas=-v -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES")
if (CMAKE_MAJOR_VERSION LESS 4 AND CMAKE_MINOR_VERSION LESS 3)
  # workaround https://github.com/Kitware/CMake/commit/99abebdea01b9ef73e091db5594553f7b1694a1b
  message(STATUS "Applied CUDA C++11 workaround on CMake < 3.3")
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --std c++11")
endif()

endif(CUDA_FOUND)
