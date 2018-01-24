program example2

        use,intrinsic :: iso_fortran_env, only: wp => real64
        use json_module

        implicit none

        type(json_core) :: json
        type(json_value),pointer :: p, inp, inp2

        ! initialize the class
        call json%initialize()

        ! initialize the structure:
        call json%create_object(p,'')

        ! add an "inputs" object to the structure:
        call json%create_object(inp,'inputs')
        call json%create_object(inp2,'test')
        call json%add(p, inp) !add it to the root
        call json%add(inp, inp2) !add it to the root

        ! add some data to inputs:
        call json%add(inp, 't0', 0.1_wp)
        call json%add(inp, 'tf', 1.1_wp)
        call json%add(inp, 'x0', 9999.0000d0)
        call json%add(inp, 'integer_scalar', 787)
        call json%add(inp, 'integer_array', [2,4,99])
        call json%add(inp, 'names', ['aaa','bbb','ccc'])
        call json%add(inp, 'logical_scalar', .true.)
        call json%add(inp, 'logical_vector', [.true., .false., .true.])
        call json%add(inp2, 'tkljrhgkue0', 0.1_wp)

        ! write the file:
        call json%print(p,'example2.json')

        !cleanup:
        call json%destroy(p)
        if (json%failed()) stop 1

    end program example2
