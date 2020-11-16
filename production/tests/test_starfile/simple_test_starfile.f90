! just write out a few fields into a starfile and retrieve one string
program simple_test_starfile
    use, intrinsic :: ISO_C_Binding, only: C_long
    use simple_starfile_wrappers
    implicit none

    type(starfile_table_type) :: sfile

    logical :: aresult
    character(len=:), allocatable :: retrieved_string
    real(kind=8) :: amt, aml
    integer(C_long) :: object_id, num_objects

    ! step 1: write star-file
    ! alloc and open output file
    call starfile_table__new(sfile)
    call starfile_table__open_ofile(sfile, "outputfile.star")
    call starfile_table__clear(sfile)
    ! first segment ("name1")
    call starfile_table__setName(sfile, "name1")
    ! first segment is list
    call starfile_table__setislist(sfile, .true.)
    call starfile_table__addObject(sfile)
    ! add a comment for good measure
    call starfile_table__setComment(sfile, "this_is_a_comment")
    ! add 3 fields, 1 string and 2 doubles
    call starfile_table__setValue_string(sfile, EMDL_MICROGRAPH_NAME, "this_is_a_string")
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_TOTAL, 99._8)
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_EARLY, 42._8)
    ! write first segment
    call starfile_table__write_ofile(sfile)
    call starfile_table__clear(sfile)
    ! create second segment ("name2")
    call starfile_table__setName(sfile, "name2")
    call starfile_table__addObject(sfile)
    ! this one is not a list
    call starfile_table__setIsList(sfile, .false.)
    ! 4 double values
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, 123._8)
    call starfile_table__addObject(sfile)
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, 456._8)
    call starfile_table__addObject(sfile)
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, 789._8)
    call starfile_table__addObject(sfile)
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, 101112._8)
    ! write second segment
    call starfile_table__write_ofile(sfile)
    call starfile_table__clear(sfile)
    ! create third segment ("name3")
    call starfile_table__setName(sfile, "name3")
    ! this one is a list again
    call starfile_table__setIsList(sfile, .true.)
    call starfile_table__addObject(sfile)
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, 523._8)

    ! write and deallocate
    call starfile_table__write_ofile(sfile)
    call starfile_table__close_ofile(sfile)
    call starfile_table__delete(sfile)

    ! reallocate and read first segment
    call starfile_table__new(sfile)
    call starfile_table__read(sfile, "outputfile.star", "name1")
    ! see if we can retrieve string correctly
    aresult = starfile_table__getValue_string(sfile, EMDL_MICROGRAPH_NAME, retrieved_string)
    write (*,*) 'result=', aresult, ' retrieved_string=', retrieved_string
    ! this one should fail
    aresult = starfile_table__getValue_string(sfile, EMDL_MLMODEL_REF_IMAGE, retrieved_string)
    if (aresult) then
       write (*,*) 'result=', aresult, ' retrieved_string=', retrieved_string
    else
       write (*,*) 'result=', aresult
    end if
    ! now read the other two fields; the first one should go through, the other should fail
    aresult = starfile_table__getValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_TOTAL, amt)
    write (*,*) 'aresult = ', aresult, ' ; amt = ', amt
    aresult = starfile_table__getValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, aml)
    write (*,*) 'aresult = ', aresult, ' ; aml = ', aml
    write (*,*) '-----------------'
    ! now read in second segment
    call starfile_table__read(sfile, "outputfile.star","name2")
    ! iterate through list
    object_id = starfile_table__firstobject(sfile)
    num_objects = starfile_table__numberofobjects(sfile)
    do while( (object_id < num_objects) .and. (object_id >= 0) )
            aresult = starfile_table__getValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, aml)
            write (*,*) 'aresult = ', aresult, ' ; aml = ', aml
        object_id = starfile_table__nextobject(sfile)
    end do
    ! now read in third segment
    call starfile_table__read(sfile, "outputfile.star","name3")
    write (*,*) '-----------------'
    aresult = starfile_table__getValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, aml)
    write (*,*) 'aresult = ', aresult, ' ; aml = ', aml
    ! deallocate
    call starfile_table__delete(sfile)

  end program simple_test_starfile
