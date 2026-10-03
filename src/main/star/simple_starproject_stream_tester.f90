!@descr: unit tests for the stream STAR export (simple_starproject_stream)
! The fixture is a project of two micrographs, each with a stack of two particles, in a fresh
! fixture directory; one particle is deselected. Its particles are exported as particles2D.star,
! and the rows are read back as text: the field holding '@' is the image name, the field ending
! in '_intg.mrc' the micrograph name.
module simple_starproject_stream_tester
use simple_test_utils
use simple_string,             only: string
use simple_string_utils,       only: int2str
use simple_fileio,             only: file_exists, simple_getcwd
use simple_parameters,         only: parameters
use simple_sp_project,         only: sp_project
use simple_starproject_stream, only: starproject_stream
implicit none
private
public :: run_all_starproject_stream_tests

integer, parameter :: LINELEN = 2048
integer, parameter :: MAXROWS = 10

contains

    subroutine run_all_starproject_stream_tests()
        write(*,'(A)') '**** running all stream STAR export tests ****'
        call test_particles2D_names()
    end subroutine run_all_starproject_stream_tests

    !> every exported particle names its image in its stack (index@stack) and its micrograph; a
    !! deselected particle is not exported
    subroutine test_particles2D_names()
        type(starproject_stream) :: star
        type(parameters)         :: params
        type(sp_project)         :: spproj
        type(string)             :: cwd_saved, root, cwd
        character(len=LINELEN)   :: images(MAXROWS), mics(MAXROWS)
        integer                  :: nfail0, n
        write(*,'(A)') 'test_particles2D_names'
        nfail0 = tests_failed
        call enter_fixture('star_stream_particles2D', cwd_saved, root)
        call simple_getcwd(cwd)
        call make_two_stack_project(spproj, cwd)
        call star%stream_export_particles_2D(params, spproj, cwd, optics_set=.true.)
        call assert_true(file_exists(string('particles2D.star')), 'particles2D.star is written')
        call read_particle_rows('particles2D.star', images, mics, n)
        call assert_int(3, n, 'one row per selected particle')
        if( n == 3 )then
            call assert_true(is_image(images(1), 1, '/stk_1.mrcs'), 'particle 1 is image 1 of stack 1')
            call assert_true(ends_with(mics(1), '/mic_1_intg.mrc'), 'on micrograph 1')
            call assert_true(is_image(images(2), 1, '/stk_2.mrcs'), 'particle 3 is image 1 of stack 2')
            call assert_true(ends_with(mics(2), '/mic_2_intg.mrc'), 'on micrograph 2')
            call assert_true(is_image(images(3), 2, '/stk_2.mrcs'), 'particle 4 is image 2 of stack 2')
            call assert_true(ends_with(mics(3), '/mic_2_intg.mrc'), 'on micrograph 2')
        endif
        call spproj%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_particles2D_names

    ! ---- fixtures ------------------------------------------------------------

    ! two micrographs in @p dir, each with a stack of two particles, with the stacks' image counts
    ! (nptcls_stk) and the particles' image indices (indstk) that release 4 requires; particle 2 is
    ! deselected; one optics group
    subroutine make_two_stack_project( spproj, dir )
        type(sp_project), intent(inout) :: spproj
        class(string),    intent(in)    :: dir
        integer, parameter :: STKINDS(4) = [1, 1, 2, 2], INDSTK(4) = [1, 2, 1, 2], STATES(4) = [1, 0, 1, 1]
        integer :: i
        call spproj%os_optics%new(1, is_ptcl=.false.)
        call spproj%os_optics%set(1, 'ogid', 1)
        call spproj%os_optics%set(1, 'smpd', 1.0)
        call spproj%os_optics%set(1, 'kv',   300.)
        call spproj%os_optics%set_state(1, 1)
        call spproj%os_mic%new(2, is_ptcl=.false.)
        call spproj%os_stk%new(2, is_ptcl=.false.)
        do i = 1,2
            call spproj%os_mic%set(i, 'intg',    dir//'/mic_'//int2str(i)//'_intg.mrc')
            call spproj%os_mic%set(i, 'imgkind', 'mic')
            call spproj%os_mic%set_state(i, 1)
            call spproj%os_stk%set(i, 'stk',   dir//'/stk_'//int2str(i)//'.mrcs')
            call spproj%os_stk%set(i, 'box',   64)
            call spproj%os_stk%set(i, 'fromp', 2 * i - 1)
            call spproj%os_stk%set(i, 'top',   2 * i)
            call spproj%os_stk%set(i, 'nptcls_stk', 2)
            call spproj%os_stk%set_state(i, 1)
        enddo
        call spproj%os_ptcl2D%new(4, is_ptcl=.true.)
        do i = 1,4
            call spproj%os_ptcl2D%set(i, 'stkind', STKINDS(i))
            call spproj%os_ptcl2D%set(i, 'indstk', INDSTK(i))
            call spproj%os_ptcl2D%set(i, 'ogid',   1)
            call spproj%os_ptcl2D%set_state(i, STATES(i))
        enddo
    end subroutine make_two_stack_project

    ! the image and micrograph names of the particle rows of @p fname (the lines with an '@'),
    ! in file order; @p n rows
    subroutine read_particle_rows( fname, images, mics, n )
        character(len=*),       intent(in)  :: fname
        character(len=LINELEN), intent(out) :: images(MAXROWS), mics(MAXROWS)
        integer,                intent(out) :: n
        character(len=LINELEN) :: line
        integer :: funit, ios
        n      = 0
        images = ''
        mics   = ''
        open(newunit=funit, file=fname, status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        do
            read(funit, '(A)', iostat=ios) line
            if( ios /= 0 ) exit
            if( index(line, '@') == 0 ) cycle
            if( n == MAXROWS ) exit
            n = n + 1
            images(n) = field_with(line, '@')
            mics(n)   = field_with(line, '_intg.mrc')
        enddo
        close(funit)
    end subroutine read_particle_rows

    ! the blank-delimited field of @p line that contains @p key; '' when none does
    function field_with( line, key ) result( field )
        character(len=*), intent(in) :: line, key
        character(len=LINELEN) :: field
        integer :: pos, first, last
        field = ''
        pos   = index(line, key)
        if( pos == 0 ) return
        first = pos
        do while( first > 1 )
            if( is_blank(line(first-1:first-1)) ) exit
            first = first - 1
        enddo
        last = pos + len(key) - 1
        do while( last < len_trim(line) )
            if( is_blank(line(last+1:last+1)) ) exit
            last = last + 1
        enddo
        field = line(first:last)
    end function field_with

    logical function is_blank( c )
        character(len=1), intent(in) :: c
        is_blank = c == ' ' .or. c == achar(9)
    end function is_blank

    ! @p field names image @p ind of a stack whose path ends in @p stk_suffix
    logical function is_image( field, ind, stk_suffix )
        character(len=*), intent(in) :: field, stk_suffix
        integer,          intent(in) :: ind
        character(len=:), allocatable :: prefix
        prefix   = int2str(ind)//'@'
        is_image = .false.
        if( len_trim(field) < len(prefix) ) return
        if( field(1:len(prefix)) /= prefix ) return
        is_image = ends_with(field, stk_suffix)
    end function is_image

    logical function ends_with( field, suffix )
        character(len=*), intent(in) :: field, suffix
        integer :: lf
        lf        = len_trim(field)
        ends_with = lf >= len(suffix)
        if( ends_with ) ends_with = field(lf-len(suffix)+1:lf) == suffix
    end function ends_with

end module simple_starproject_stream_tester
