if(CPACK_GENERATOR)
  # common package information
  set(CPACK_PACKAGE_NAME "${PACKAGE_NAME}")
  set(CPACK_PACKAGE_VENDOR "Elmlund Lab")
  set(CPACK_PACKAGE_VERSION "${PACKAGE_VERSION}")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "cryo Electron Microscope Image Processor")
  set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "${CPACK_PACKAGE_DESCRIPTION_SUMMARY}")
  #set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_CURRENT_LIST_DIR}/README.md")
  #set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_LIST_DIR}/LICENSE")
  set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_LIST_DIR}/README.md")
  set(CPACK_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
  set(CPACK_PACKAGE_CONTACT "Hans Elmlund <Hans.Elmlund@monash.edu>")
  set(CPACK_OUTPUT_FILE_PREFIX packages)
  set(CPACK_PACKAGE_RELOCATABLE true)
  set(CPACK_MONOLITHIC_INSTALL true)

  # Prefix Debug/Nightly release
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}")
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
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "binutils, libfftw3-dev, libfftw3-bin,"
      "gfortran, cmake,"
      " libc6, libgcc1, libpng12-0, libstdc++6, ")
    set(CPACK_DEBIAN_PACKAGE_SECTION "web")
    set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
    set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "http://simplecryoem.org")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}"
      "-${CPACK_PACKAGE_VERSION}"
      "-${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")
    # RPM package
  elseif(CPACK_GENERATOR MATCHES "RPM")
    # https://github.com/pld-linux/hhvm
    # https://github.com/hhvm/packaging/tree/master/hhvm/rpm/fedora20/rpmbuild/
    set(CPACK_RPM_PACKAGE_REQUIRES "binutils-devel, cmake >= 2.8.7, "
      "gcc >= 6:4.6.0, libstdc++-devel >= 6:4.3, ")
    set(CPACK_RPM_PACKAGE_GROUP "Development/Languages")
    set(CPACK_RPM_PACKAGE_LICENSE "GPL3+ and BSD")
    set(CPACK_RPM_PACKAGE_URL "http://simplecryoem.org")
    set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}"
      "-${CPACK_PACKAGE_VERSION}"
      "-${CPACK_RPM_PACKAGE_ARCHITECTURE}")
  endif()
  include(CPack)
  endif()
