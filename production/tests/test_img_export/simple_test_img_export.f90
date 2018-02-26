 !! Modified by Michael Eager, Feb 2018
program simple_test_img_export
    include 'simple_lib.f08'
#ifdef LIBGD
    use simple_test_libgd_io
#endif
    use simple_test_pnm_io
    implicit none


    character(len=8)      :: datestr
    character(len=STDLEN) :: folder, cmd
    character(len=300)    :: command


    call seed_rnd
    call date_and_time(date=datestr)
    folder = trim('./SIMPLE_TEST_PNG_'//datestr)
    command = 'mkdir ' // trim( folder )//'|| true'
    call exec_cmdline( trim(command) )
    call simple_chdir( trim(folder) )

#ifdef LIBGD
    call test_jpeg_io
    call test_png_io
#endif

    call test_pnm_io

end program simple_test_img_export
