# production/CMakeLists.txt
#
# Build production executables and (optionally) tests.

# ------------------------------------------------------------------------------
# Main executables
# ------------------------------------------------------------------------------
set(SIMPLE_EXE          simple_exec)
set(SIMPLE_TEST_EXE     simple_test_exec)
set(SINGLE_EXE          single_exec)
set(STREAM_EXE          simple_stream)
set(SIMPLE_PRIVATE_EXE  simple_private_exec)
set(SIMPLE_WORKER_EXE   simple_persistent_worker)

add_executable(${SIMPLE_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_exec.f90")
add_executable(${SINGLE_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/single_exec.f90")
add_executable(${STREAM_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_stream.f90")
add_executable(${SIMPLE_PRIVATE_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_private_exec.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_private_exec_driver.f90")
add_executable(${SIMPLE_WORKER_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/simple_persistent_worker.f90")

set(SIMPLE_EXECUTABLES
    ${SIMPLE_EXE}
    ${SINGLE_EXE}
    ${STREAM_EXE}
    ${SIMPLE_PRIVATE_EXE}
    ${SIMPLE_WORKER_EXE}
)

# simple_test_exec is the one test executable: built only with BUILD_TESTS=ON
# (its library sources are filtered out in src/ otherwise).
if(BUILD_TESTS)
    add_executable(${SIMPLE_TEST_EXE}
        "${CMAKE_CURRENT_SOURCE_DIR}/simple_test_exec.f90")
    list(APPEND SIMPLE_EXECUTABLES ${SIMPLE_TEST_EXE})
endif()

foreach(exe IN LISTS SIMPLE_EXECUTABLES)
    set_target_properties(${exe} PROPERTIES
        LINKER_LANGUAGE Fortran
    )
    if(USE_COARRAYS AND exe STREQUAL SIMPLE_PRIVATE_EXE)
        set_target_properties(${exe} PROPERTIES
            Fortran_MODULE_DIRECTORY "${SIMPLE_COARRAY_MODULE_DIRECTORY}"
        )
        target_link_libraries(${exe}
            PRIVATE
                ${SIMPLE_COARRAYLIB}
                ${SIMPLE_LIBRARIES}
        )
    else()
        target_link_libraries(${exe}
            PRIVATE
                ${SIMPLELIB}
                ${SIMPLE_LIBRARIES}
        )
    endif()
    if( USE_OPENMP )
      target_link_options(${exe} PUBLIC -fopenmp)
    endif()
endforeach()

install(TARGETS ${SIMPLE_EXECUTABLES}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# ------------------------------------------------------------------------------
# Complete UI JSON validation
# ------------------------------------------------------------------------------
# The UI descriptor is a client contract.  Regenerate and validate the complete
# emitted registry whenever simple_exec is rebuilt, independently of the
# optional CTest suite. Use the serial front end so metadata validation does not
# initialize an optional coarray/MPI runtime.
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(UI_JSON_VALIDATION_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/simple_ui.json")
set(UI_JSON_VALIDATION_STAMP "${CMAKE_CURRENT_BINARY_DIR}/simple_ui.json.valid")
add_custom_command(
  OUTPUT "${UI_JSON_VALIDATION_STAMP}"
  COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/validate_complete_ui_json.py"
    --executable "$<TARGET_FILE:${SIMPLE_EXE}>"
    --output "${UI_JSON_VALIDATION_OUTPUT}"
    --stamp "${UI_JSON_VALIDATION_STAMP}"
  DEPENDS
    ${SIMPLE_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/validate_complete_ui_json.py"
  BYPRODUCTS "${UI_JSON_VALIDATION_OUTPUT}"
  COMMENT "Generating and validating the complete UI JSON registry"
  VERBATIM
)
add_custom_target(validate_complete_ui_json ALL
  DEPENDS "${UI_JSON_VALIDATION_STAMP}"
)

# The stream interface is a reviewed acquisition-facing screen. Its exact
# snapshot is deliberately opt-in: ordinary builds verify valid JSON without
# blocking routine wording or layout work on a separate approval step.
set(STREAM_UI_CONTRACT_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/simple_stream_ui_contract.json")
set(STREAM_UI_CONTRACT_STAMP "${CMAKE_CURRENT_BINARY_DIR}/simple_stream_ui_contract.valid")
set(STREAM_UI_CONTRACT_UPDATE_STAMP "${CMAKE_CURRENT_BINARY_DIR}/simple_stream_ui_contract.updated")
set(STREAM_UI_CONTRACT_REFERENCE "${CMAKE_CURRENT_SOURCE_DIR}/stream_ui_contract.json")
add_custom_command(
  OUTPUT "${STREAM_UI_CONTRACT_STAMP}"
  COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/validate_complete_ui_json.py"
    --executable "$<TARGET_FILE:${SIMPLE_PRIVATE_EXE}>"
    --output "${STREAM_UI_CONTRACT_OUTPUT}"
    --stamp "${STREAM_UI_CONTRACT_STAMP}"
    --check-stream-contract
  DEPENDS
    ${SIMPLE_PRIVATE_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/validate_complete_ui_json.py"
    "${STREAM_UI_CONTRACT_REFERENCE}"
  BYPRODUCTS "${STREAM_UI_CONTRACT_OUTPUT}"
  COMMENT "Reviewing the stream UI snapshot"
  VERBATIM
)
add_custom_target(validate_stream_ui_contract
  DEPENDS "${STREAM_UI_CONTRACT_STAMP}"
)
add_custom_target(update_stream_ui_contract
  COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/validate_complete_ui_json.py"
    --executable "$<TARGET_FILE:${SIMPLE_PRIVATE_EXE}>"
    --output "${STREAM_UI_CONTRACT_OUTPUT}"
    --stamp "${STREAM_UI_CONTRACT_UPDATE_STAMP}"
    --update-stream-contract
  DEPENDS
    ${SIMPLE_PRIVATE_EXE}
    "${CMAKE_CURRENT_SOURCE_DIR}/validate_complete_ui_json.py"
  COMMENT "Updating the reviewed stream UI JSON reference"
  VERBATIM
)

# ------------------------------------------------------------------------------
# UI display-default audit
# ------------------------------------------------------------------------------
# This is intentionally opt-in: it is review evidence, not a runtime artifact,
# and it reports ambiguities for human resolution without failing normal builds.
set(UI_DEFAULT_AUDIT_JSON "${CMAKE_CURRENT_BINARY_DIR}/simple_ui_default_audit.json")
set(UI_DEFAULT_AUDIT_MARKDOWN "${CMAKE_CURRENT_BINARY_DIR}/simple_ui_default_audit.md")
add_custom_command(
  OUTPUT "${UI_DEFAULT_AUDIT_JSON}" "${UI_DEFAULT_AUDIT_MARKDOWN}"
  COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_SOURCE_DIR}/scripts/default_audit.py"
    --source-root "${CMAKE_SOURCE_DIR}"
    --ui-executable "$<TARGET_FILE:${SIMPLE_PRIVATE_EXE}>"
    --output-json "${UI_DEFAULT_AUDIT_JSON}"
    --output-markdown "${UI_DEFAULT_AUDIT_MARKDOWN}"
  DEPENDS
    ${SIMPLE_PRIVATE_EXE}
    "${CMAKE_SOURCE_DIR}/scripts/default_audit.py"
  COMMENT "Generating UI display-default audit report"
  VERBATIM
)
add_custom_target(audit_ui_defaults
  DEPENDS "${UI_DEFAULT_AUDIT_JSON}" "${UI_DEFAULT_AUDIT_MARKDOWN}"
)

# ------------------------------------------------------------------------------
# Tests  (doc/refactoring_notes/completed/uniform_test_environment_refactoring.md)
# ------------------------------------------------------------------------------
# Two gates behind one front door, simple_test_exec:
#   fast         runs on every compile_*.sh build under a 30 s budget
#                (Phase 2 declares it; until then `provisional` stands in)
#   library      nightly, one process per coherent part of the library
#   highlevel    only by explicit CTest command; long simulated workflows
#   platform     coarray/MPI/GPU/socket/child-process cases, own launcher
# Every registration has a label, a timeout, its own working directory under
# build/test_runs/<name>, an explicit OMP_NUM_THREADS and the environment the
# SIMPLE runtime reads (SIMPLE_PATH, SIMPLE_QSYS, SIMPLE_EMAIL, PATH to the
# installed bin/). A CTest entry is an isolation unit, not an assertion:
# checks grow inside existing suites, and a new process needs a reason here
# and a raised SIMPLE_CTEST_BUDGET. The fast label as a whole must stay
# under 30 s of real time (scripts/run_fast_gate.sh); there is no per-entry
# time rule.
if(BUILD_TESTS)
  enable_testing()

  # ---- registration helper -------------------------------------------------
  # simple_add_test(<name> <label> <timeout> <omp_threads> [RUN_SERIAL] COMMAND <args...>)
  # registers `simple_test_exec <args...>` (or the given launcher) with the
  # properties every entry must carry.
  set(SIMPLE_TEST_RUNS_DIR "${CMAKE_BINARY_DIR}/test_runs")
  # SIMPLE_SEED makes seed_rnd (called by parameters%new, so by every commander)
  # draw from a fixed seed instead of /dev/urandom; distributed workers inherit
  # it, so the library and high-level entries are reproducible (plan, section 9.7,
  # stream).
  set(SIMPLE_TEST_ENV
    "SIMPLE_PATH=${CMAKE_INSTALL_PREFIX}"
    "SIMPLE_SEED=20260923"
    "SIMPLE_QSYS=local"
    "SIMPLE_EMAIL=test@localhost"
    "GFORTRAN_UNBUFFERED_ALL=y"
  )
  function(simple_add_test name label timeout nthr)
    cmake_parse_arguments(ARG "RUN_SERIAL" "" "COMMAND" ${ARGN})
    file(MAKE_DIRECTORY "${SIMPLE_TEST_RUNS_DIR}/${name}")
    add_test(NAME ${name}
      COMMAND ${ARG_COMMAND}
      WORKING_DIRECTORY "${SIMPLE_TEST_RUNS_DIR}/${name}"
    )
    set_tests_properties(${name} PROPERTIES
      LABELS "${label}"
      TIMEOUT ${timeout}
      ENVIRONMENT "OMP_NUM_THREADS=${nthr};${SIMPLE_TEST_ENV}"
      ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR};PATH=path_list_prepend:${CMAKE_INSTALL_PREFIX}/scripts"
    )
    if(ARG_RUN_SERIAL)
      set_tests_properties(${name} PROPERTIES RUN_SERIAL TRUE)
    endif()
  endfunction()

  # ---- the fast gate ------------------------------------------------------
  # Thirteen area suites over the unit tests, one process each, label `fast`;
  # every sub-suite asserts through simple_test_utils, is localhost-only and
  # runs on one thread (plan, section 5.1). `units` runs all of them in one
  # process as a developer convenience and is not registered: the entries
  # below are the gate. Measured 2026-09-22 (Debug, 1 thread): the seven
  # together about 4 s; the longest, unit_image and unit_ori, about 1.3 s.
  foreach(_area IN ITEMS core ori image numerics project ui ipc reconstruction pftc_align2D3D cart_align3D heterogeneity parallel single)
    simple_add_test(unit_${_area} fast 120 1
      COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=unit_${_area})
  endforeach()

  # ---- library suites (nightly) --------------------------------------------
  # One process per coherent part of the library, label `library`, the same
  # runner and assertions as the fast gate without its 30 s budget (plan,
  # section 5.2.1). lib_reconstruction: half-set PCG solves against gridding
  # at box 24, minutes on one thread. Review verdict 2026-09-23 (plan,
  # section 9.7, reconstruction).
  simple_add_test(lib_reconstruction library 3600 8 RUN_SERIAL
    COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=lib_reconstruction)
  # lib_cart_align3D: 5 000 simulated 1JYX particles at box 144 refined by
  # the Cartesian pose refiner and reconstructed; the run directory keeps the
  # volumes and TSVs. Review verdict 2026-09-23 (plan, section 9.7, pose).
  simple_add_test(lib_cart_align3D library 7200 8 RUN_SERIAL
    COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=lib_cart_align3D)
  # lib_heterogeneity: flex_pca deconvolution of 20 000 particles at realistic
  # noise, the PCG M-step operator at box 64 with 400 samples against the exact
  # nonuniform-DFT Gram, and the twelve-setting PCG solve sweep at box 32.
  # Review verdict 2026-09-23 (plan, section 9.7, heterogeneity).
  simple_add_test(lib_heterogeneity library 3600 8 RUN_SERIAL
    COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=lib_heterogeneity)
  # lib_single: in-process pdb2mrc coverage of the built-in molecular models.
  # The longer nanoparticle atom-statistics pipeline is an explicit highlevel entry.
  simple_add_test(lib_single library 7200 8 RUN_SERIAL
    COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=lib_single)
  # lib_stream: the stream stages that run in-process, one thread, with the
  # arguments the stream gives them (Ruben's stream tests): optics assignment
  # (a minute: the stream watcher waits LONGTIME = 60 s), picking references,
  # pick and extract. Review verdict 2026-09-24 (plan, section 9.7, stream).
  simple_add_test(lib_stream library 1800 1
    COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=lib_stream)

  # ---- high-level gates (explicit CTest runs only) --------------------------
  # These workflows use generated data and quantitative output gates. The
  # simulated suites and single workflow register the final map to known truth
  # and check dimensions, sampling, whole-volume correlation, and masked FSC.
  # Each owns the machine: RUN_SERIAL, an OpenMP team, and it may start
  # simple_exec / simple_private_exec workers, hence PATH and SIMPLE_PATH above.
  # These entries are deliberately excluded from scripts/run_fast_gate.sh,
  # which selects only the `fast` label during compilation/installation. Run
  # them explicitly with: ctest -L highlevel --output-on-failure
  simple_add_test(mini_stream_6vxx highlevel 7200 8 RUN_SERIAL
    COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=mini_stream suite=6vxx)
  simple_add_test(mini_stream_1jxy highlevel 7200 8 RUN_SERIAL
    COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=mini_stream suite=1jxy)
  simple_add_test(simulated_workflow_6vxx highlevel 7200 8 RUN_SERIAL
    COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=simulated_workflow suite=6vxx)
  simple_add_test(simulated_workflow_1jxy highlevel 7200 8 RUN_SERIAL
    COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=simulated_workflow suite=1jxy)
  simple_add_test(single_workflow highlevel 7200 8 RUN_SERIAL
    COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=single_workflow smpd=0.358 element=Pt)
  simple_add_test(pcg_recon highlevel 3600 8 RUN_SERIAL
    COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=pcg_recon)
  # reproject + simulate_particles on the embedded 6VXX volume; the only
  # explicit high-level run of either commander (simulated_workflow uses simulate_movie).
  # Review verdict 2026-09-23 (plan, section 9.7, highlevel).
  simple_add_test(simulate_particles highlevel 3600 8 RUN_SERIAL
    COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=simulate_particles nthr=8)
  simple_add_test(single_atoms_stats highlevel 3600 8 RUN_SERIAL
    COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=single_atoms_stats smpd=0.358 element=Pt)
  # stream preprocessing: five simulated movies through the stream's
  # preprocessing stage, which submits its jobs to the local queue; the other
  # stream stage tests moved to lib_stream, unit_project and forked_process,
  # abinitio2D_stream was retired (plan, section 9.7, stream).
  simple_add_test(stream_preproc highlevel 3600 8 RUN_SERIAL
    COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=stream_preproc)

  # ---- platform cases -------------------------------------------------------
  # Own process, own launcher; not part of the build gate.
  # forked process: real child processes polled on a clock (12.9 s of it),
  # excluded from the fast gate by decision (plan, section 4.6); its second
  # sub-suite, stream heartbeat, forks the seven stream stages' children.
  simple_add_test(forked_process platform 300 1 RUN_SERIAL
    COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=forked_process)
  # Registered only when CMake has the capability and the launcher.
  if(USE_FLEX_CUDA)
    # CUDA-C flex kernels against the CPU batch path; needs a device at run time
    simple_add_test(flex_gpu platform 1200 1 RUN_SERIAL
      COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=flex_gpu)
  endif()
  if(USE_OPENMP_OFFLOAD)
    # OpenMP target offload, cuFFT against FFTW, cuBLAS, the KB device forms (Cyril's,
    # simple_openmp_offload_tester since the utils review); a failed check stops
    # through THROW_HARD, and the entry also fails on its failure messages
    simple_add_test(openmp_offload platform 1200 8 RUN_SERIAL
      COMMAND $<TARGET_FILE:${SIMPLE_TEST_EXE}> test=openmp_offload nthr=8 device=0)
    set_property(TEST openmp_offload APPEND PROPERTY ENVIRONMENT "OMP_TARGET_OFFLOAD=MANDATORY")
    set_tests_properties(openmp_offload PROPERTIES FAIL_REGULAR_EXPRESSION "Fatal error|FAIL: ")
  endif()
  if(USE_COARRAYS AND Coarray_CAFRUN_EXECUTABLE)
    # Minimal capability check compiled into the coarray-linked private executable:
    # two images synchronize and transfer one known integer.
    simple_add_test(coarrays platform 60 1 RUN_SERIAL
      COMMAND "${Coarray_CAFRUN_EXECUTABLE}" -np 2 $<TARGET_FILE:${SIMPLE_PRIVATE_EXE}> --coarray-smoke)
  endif()

  # ---- process budget -------------------------------------------------------
  # A CTest entry is an isolation unit, not an assertion: checks grow inside
  # existing suites, a new process needs a stated reason, and this number only
  # goes down. Raising it is an owner decision recorded in the plan. Entries
  # that depend on a capability (coarrays, flex_gpu, openmp_offload) are counted only
  # when registered.
  set(SIMPLE_CTEST_BUDGET 28)   # 13 fast + 1 platform + 9 highlevel + 5 library
  set(_simple_expected ${SIMPLE_CTEST_BUDGET})
  if(TEST coarrays)
    math(EXPR _simple_expected "${_simple_expected} + 1")
  endif()
  if(TEST flex_gpu)
    math(EXPR _simple_expected "${_simple_expected} + 1")
  endif()
  if(TEST openmp_offload)
    math(EXPR _simple_expected "${_simple_expected} + 1")
  endif()
  get_property(_simple_tests DIRECTORY PROPERTY TESTS)
  list(LENGTH _simple_tests _simple_ntests)
  if(NOT _simple_ntests EQUAL _simple_expected)
    message(FATAL_ERROR "${_simple_ntests} CTest entries do not match SIMPLE_CTEST_BUDGET (${_simple_expected}): "
      "add checks to an existing suite, justify a new process in production/CMakeLists.txt, "
      "or lower the budget after an intentional consolidation")
  endif()
  message(STATUS "CTest entries:        ${_simple_ntests} (budget ${_simple_expected})")

endif()
