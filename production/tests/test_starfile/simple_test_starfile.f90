! just write out a few fields into a starfile and retrieve one string
program simple_test_starfile
    use simple_starfile_wrappers
    implicit none

    type(starfile_table_type) :: sfile

    logical :: aresult
    character(len=:), allocatable :: retrieved_string
    
    call starfile_table__new(sfile)
    call starfile_table__clear(sfile)
    call starfile_table__setislist(sfile, .true.)
    
    call starfile_table__addObject(sfile)
    call starfile_table__setName(sfile, "this_is_a_a_name")    
    call starfile_table__setComment(sfile, "this_is_a_comment")
    
    call starfile_table__setValue_string(sfile, EMDL_MICROGRAPH_NAME, "this_is_a_string")
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_TOTAL, 99._8)
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_EARLY, 42._8)

    ! we use the getValue routine here
    write (*,*) 'the first one should be retrieved, the second one shouldn''t'
    aresult = starfile_table__getValue_string(sfile, EMDL_MICROGRAPH_NAME, retrieved_string)
    write (*,*) 'result=', aresult, ' retrieved_string=', retrieved_string

    aresult = starfile_table__getValue_string(sfile, EMDL_MLMODEL_REF_IMAGE, retrieved_string)
    if (aresult) then
       write (*,*) 'result=', aresult, ' retrieved_string=', retrieved_string
    else
       write (*,*) 'result=', aresult
    end if
    
    call starfile_table__open_ofile(sfile, "outputfile.txt")
    call starfile_table__write_ofile(sfile)
    call starfile_table__clear(sfile)

    call starfile_table__addObject(sfile)
    call starfile_table__setName(sfile, "name2")
    call starfile_table__setIsList(sfile, .false.)
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, 123._8)
    call starfile_table__write_ofile(sfile)
    call starfile_table__close_ofile(sfile)
    call starfile_table__delete(sfile)
    
  end program simple_test_starfile
