if(CPACK_GENERATOR)
  get_filename_component(cpack_build_dir "${CMAKE_BINARY_DIR}" ABSOLUTE)
  get_filename_component(cpack_source_dir "${CMAKE_SOURCE_DIR}" ABSOLUTE)
  # if("${cpack_build_dir}" STREQUAL "${cpack_source_dir}")
endif()

if(CPACK_GENERATOR)

  #if("${CPACK_GENERATOR}" STREQUAL "ON")
  # default settings
  set(CPACK_GENERATOR_TGZ ON)
  set(CPACK_SOURCE_TBZ2  OFF)
  set(CPACK_SOURCE_TZ    OFF)
  set(CPACK_SOURCE_TXZ   OFF)
  set(CPACK_SOURCE_ZIP   ON)
  #endif()
  # common package information
  set(CPACK_PACKAGE_NAME "SIMPLE")
  set(CPACK_PACKAGE_VENDOR "Elmlund & Elmlund Lab -- Monash University -- 2017")
  set(CPACK_PACKAGE_VERSION "${SIMPLE_VERSION}")
  set(CPACK_PACKAGE_VERSION_MAJOR "${SIMPLE_VERSION_MAJOR}")
  set(CPACK_PACKAGE_VERSION_MINOR "${SIMPLE_VERSION_MINOR}")
  set(CPACK_PACKAGE_VERSION_PATCH "${SIMPLE_GIT_VERSION}")
  set(CPACK_SET_DESTDIR "/usr/local")
  # Useful descriptions for components
  SET(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "SIMPLE cryoEM library")
  SET(CPACK_COMPONENT_DOCUMENTATION_NAME "Doxygen documentation")
  SET(CPACK_COMPONENT_HEADERS_NAME "Developmental headers")
  SET(CPACK_COMPONENT_CMAKE_NAME "CMake support")
  set(CPACK_PACKAGE_DESCRIPTION
    "SIMPLE is a SIngle-particle cryo electron Microscope Image Processing
Engine, focusing on ab initio 3D reconstruction of low-symmetry single-particles.
The SIMPLE back-end consists of an object-oriented numerical library written in
modern Fortran. The SIMPLE front-end consists of a few standalone, interoperable
components developed according to the 'Unix toolkit philosophy'.

HOMEPAGE: https://simplecryoem.com
LICENCE: GNU Public licence (version 3 or newer)
AUTHORS: Hans Elmlund, Dominika Elmlund, Cyril Reboul, Michael Eager, Fred Bonnet

 What is new in SIMPLE release ${SIMPLE_VERSION}

– A new DDD movie pre-processing program unblur that implements motion
  correction based the same principal strategy as Grigorieff’s program (hence
  the name). There are two important differences: automatic weighting of the
  frames using a correlation-based M-estimator and stochastic continuous
  optimisation of the shift parameters. This enables analysis of movies with
  severe pathologies due to radiation damage or extreme drift.

– A new program unblur_tomo for movie processing of tomographic tilt-series.

– Improved simultaneous 2D alignment and clustering with prime2D using a hybrid
  extremal/stochastic hill-climbing search approach, Wiener restoration-based
  CTF correction and acceleration of the search using Hadamard projection
  matching. It is now possible to generate a sub-nanometer resolution ab initio
  3D reconstruction from class averages obtained with prime2D in a about 10
  minutes on a laptop (MacBook Pro mid 2015, 2.8 GHz Intel i7, four physical
  cores).

– Improved ab initio 3D reconstruction from class averages using stochastic
  neighbourhood hill-climbing for initialisation of the 3D orientation search,
  improving the success rate from around 40% to 90-100% for more challenging
  starting model generation problems, executed with program ini3D_from_cavgs

– Serial CPU code optimisation through data re-organisation and pipelining.

– Improved parallel CPU performance through load balancing as well as data and
  algorithm re-organisation. It is now possible to process data sets of
  realistic size on laptops or lightweight workstations that cost less than
  2,000 USD.

– High-level workflows for 2D analysis and initial 3D model generation that
  automate initialisation, update of search parameters and dynamic down-scaling
  of the images for improved performance.

")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "SIMPLE is a SIngle-particle cryo electron Microscope Image Processing Engine")
  set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "${CPACK_PACKAGE_DESCRIPTION_SUMMARY}")
  #set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/README.txt")
  set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
  set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.txt")
  set(CPACK_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
  set(CPACK_PACKAGE_CONTACT "Hans Elmlund <Hans.Elmlund@monash.edu>")
  set(CPACK_OUTPUT_FILE_PREFIX packages)
  set(CPACK_PACKAGE_RELOCATABLE true)
  set(CPACK_MONOLITHIC_INSTALL true)
  SET(CPACK_PACKAGE_INSTALL_DIRECTORY "simple ${SIMPLE_VERSION}.${SIMPLE_GIT_VERSION}")
  SET(CPACK_SOURCE_PACKAGE_FILE_NAME "simple-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
  # Prefix Debug/Nightly release
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
  if(NIGHTLY)
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}-nightly")
    #execute_process(COMMAND "date +%Y.%m.%d" OUTPUT_VARIABLE NIGHTLY_DATE)
  endif()
  if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}-dbg")
  endif()

  # default package generators
  if(APPLE)
    set(PACKAGE_GENERATOR "PackageMaker")
    set(PACKAGE_SOURCE_GENERATOR "TGZ;ZIP")
  elseif(UNIX)
    set(PACKAGE_GENERATOR "DEB;RPM")
    set(PACKAGE_SOURCE_GENERATOR "TGZ;ZIP")
  else()
    set(PACKAGE_GENERATOR "ZIP")
    set(PACKAGE_SOURCE_GENERATOR "ZIP")
  endif()

  # Mac OS X package
  if(CPACK_GENERATOR MATCHES "PackageMaker|DragNDrop")
    set(CPACK_PACKAGE_FILE_NAME
      "${CPACK_PACKAGE_FILE_NAME}-${CPACK_PACKAGE_VERSION}")
    set(CPACK_PACKAGING_INSTALL_PREFIX /usr/local)
    # Debian package
  elseif(CPACK_GENERATOR MATCHES "DEB")
    # https://github.com/hhvm/packaging/tree/master/hhvm/deb
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "gcc, build-essential, binutils, libfftw3-dev, gfortran, cmake, gnuplot")
    set(CPACK_DEBIAN_PACKAGE_OPTIONAL "nvidia-cuda-dev, libopenmpi-dev, libjpeg9-dev, libsqlite3-dev ")
    set(CPACK_DEBIAN_PACKAGE_SECTION "science")
    set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
    set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "https://simplecryoem.org")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}"
      "-${CPACK_PACKAGE_VERSION}"
      "-${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")
    # RPM package
  elseif(CPACK_GENERATOR MATCHES "RPM")
    # https://github.com/pld-linux/hhvm
    # https://github.com/hhvm/packaging/tree/master/hhvm/rpm/fedora20/rpmbuild/
    set(CPACK_RPM_PACKAGE_REQUIRES "binutils-devel, cmake >= 2.8.7, "
      "gcc >= 6:4.6.0, libfftw3-devel, gnuplot   ")
    set(CPACK_RPM_PACKAGE_GROUP "Applications/Engineering")
    set(CPACK_RPM_PACKAGE_LICENSE "GPL3+ ")
    set(CPACK_RPM_PACKAGE_URL "https://simplecryoem.org")
    set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}"
      "-${CPACK_PACKAGE_VERSION}"
      "-${CPACK_RPM_PACKAGE_ARCHITECTURE}")
  endif()
  include(CPack)
  # message( STATUS "CPack setup complete")
endif()
