set(PATH $ENV{PATH})
set(LD_LIBRARY_PATH $ENV{LD_LIBRARY_PATH})
configure_file(nice_project/wsgi.py.in nice_project/wsgi.py @ONLY)
configure_file(${SRC_PATH}/nice_local.py.in ${SIMPLE_PATH}/nice_local @ONLY FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)