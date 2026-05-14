!@descr: unit tests for simple_starfile
!        Covers: init/complete lifecycle, tmp-file staging, atomic rename,
!        write_optics_table, write_mics_table, write_ptcl2D_table,
!        write_ptcl2D_table_parallel, empty-table fast paths, state-filtered
!        rows, parallel write ordering, and relative-path relativisation.
module simple_starfile_tester
    use omp_lib
    use simple_core_module_api
    use simple_starfile
    use simple_starfile_wrappers
    use simple_test_utils
    implicit none

    private
    public :: run_all_starfile_tests

    character(len=*), parameter :: TMPDIR = "test_starfile_tmp"

contains

    !=======================================================================
    !  ENTRY POINT
    !=======================================================================

    subroutine run_all_starfile_tests()
        write(*,'(A)') '**** running all simple_starfile tests ****'
        call setup_tmpdir()
        ! --- lifecycle ---
        call test_init_creates_no_files()
        call test_complete_renames_tmp()
        call test_complete_noop_without_tmp()
        call test_complete_overwrites_existing()
        call test_init_removes_stale_tmp()
        ! --- optics table ---
        call test_write_optics_table_basic()
        call test_write_optics_table_state_filter()
        call test_write_optics_table_empty()
        call test_write_optics_table_all_fields()
        ! --- micrographs table ---
        call test_write_mics_table_basic()
        call test_write_mics_table_state_filter()
        call test_write_mics_table_empty()
        call test_write_mics_table_all_fields()
        ! --- particles table (serial) ---
        call test_write_ptcl2D_table_basic()
        call test_write_ptcl2D_table_state_filter()
        call test_write_ptcl2D_table_bad_stkind()
        call test_write_ptcl2D_table_indstk_explicit()
        call test_write_ptcl2D_table_with_mics()
        ! --- particles table (parallel) ---
        call test_write_ptcl2D_parallel_basic()
        call test_write_ptcl2D_parallel_empty()
        call test_write_ptcl2D_parallel_count()
        call test_write_ptcl2D_parallel_vs_serial()
        ! --- verbose flag does not alter output ---
        call test_verbose_flag_no_content_change()
        call report_summary()
    end subroutine run_all_starfile_tests

    !=======================================================================
    !  HELPERS
    !=======================================================================

    subroutine setup_tmpdir()
        logical :: ex
        inquire(file=TMPDIR, exist=ex)
        if( ex ) call exec_cmdline('rm -rf ' // TMPDIR)
        call exec_cmdline('mkdir -p ' // TMPDIR)
    end subroutine setup_tmpdir

    ! Return a string path inside TMPDIR
    function tmp(fname) result(s)
        character(len=*), intent(in) :: fname
        type(string) :: s
        s = TMPDIR // '/' // fname
    end function tmp

    ! Count objects in a named table inside a STAR file
    integer function count_objects(path, tablename)
        type(string),     intent(in)           :: path
        character(len=*), intent(in), optional :: tablename
        type(starfile_table_type) :: tbl
        integer :: obj
        count_objects = 0
        call starfile_table__new(tbl)
        if( present(tablename) ) then
            call starfile_table__read(tbl, path, tablename)
        else
            call starfile_table__read(tbl, path)
        end if
        obj = starfile_table__firstobject(tbl)
        if( obj < 0 ) then
            call starfile_table__delete(tbl)
            return
        end if
        count_objects = 1
        do
            obj = starfile_table__nextobject(tbl)
            if( obj < 0 ) exit
            count_objects = count_objects + 1
        end do
        call starfile_table__delete(tbl)
    end function count_objects

    ! Read back an int field from the first object of a named table
    integer function readback_int(path, emdl_label, tablename)
        type(string),     intent(in)           :: path
        integer,          intent(in)           :: emdl_label
        character(len=*), intent(in), optional :: tablename
        type(starfile_table_type) :: tbl
        integer :: iv
        logical :: ok
        readback_int = -999999
        call starfile_table__new(tbl)
        if( present(tablename) ) then
            call starfile_table__read(tbl, path, tablename)
        else
            call starfile_table__read(tbl, path)
        end if
        if( starfile_table__firstobject(tbl) < 0 ) then
            call starfile_table__delete(tbl)
            return
        end if
        ok = starfile_table__getValue_int(tbl, emdl_label, iv)
        if( ok ) readback_int = iv
        call starfile_table__delete(tbl)
    end function readback_int

    ! Read back a double field from the first object of a named table
    real(dp) function readback_double(path, emdl_label, tablename)
        type(string),     intent(in)           :: path
        integer,          intent(in)           :: emdl_label
        character(len=*), intent(in), optional :: tablename
        type(starfile_table_type) :: tbl
        real(dp) :: dv
        logical  :: ok
        readback_double = -1.0_dp
        call starfile_table__new(tbl)
        if( present(tablename) ) then
            call starfile_table__read(tbl, path, tablename)
        else
            call starfile_table__read(tbl, path)
        end if
        if( starfile_table__firstobject(tbl) < 0 ) then
            call starfile_table__delete(tbl)
            return
        end if
        ok = starfile_table__getValue_double(tbl, emdl_label, dv)
        if( ok ) readback_double = dv
        call starfile_table__delete(tbl)
    end function readback_double

    ! Build a minimal oris object with n rows, all state=1
    subroutine make_oris_n(o, n, is_ptcl)
        type(oris),  intent(out) :: o
        integer,     intent(in)  :: n
        logical,     intent(in)  :: is_ptcl
        call o%new(n, is_ptcl=is_ptcl)
        call o%set_all2single('state', 1.0)
    end subroutine make_oris_n

    !=======================================================================
    !  SECTION 1: LIFECYCLE (init / complete)
    !=======================================================================

    ! 1.1 — init creates neither the final file nor the .tmp file
    subroutine test_init_creates_no_files()
        type(starfile) :: sf
        type(string)   :: fname
        write(*,'(A)') 'test_init_creates_no_files'
        fname = tmp('init_no_files.star')
        call sf%init(fname)
        call assert_false(file_exists(fname),           'init: final file must not exist')
        call assert_false(file_exists(fname%to_char()//'.tmp'), 'init: .tmp must not exist')
    end subroutine test_init_creates_no_files

    ! 1.2 — complete renames .tmp → final
    subroutine test_complete_renames_tmp()
        type(starfile) :: sf
        type(string)   :: fname
        write(*,'(A)') 'test_complete_renames_tmp'
        fname = tmp('complete_renames.star')
        call sf%init(fname)
        ! Manufacture a .tmp file manually
        call write_dummy_file(fname%to_char()//'.tmp')
        call sf%complete()
        call assert_true (file_exists(fname),                    'complete: final file must exist')
        call assert_false(file_exists(fname%to_char()//'.tmp'),  'complete: .tmp must be gone')
    end subroutine test_complete_renames_tmp

    ! 1.3 — complete is a no-op when no .tmp file exists
    subroutine test_complete_noop_without_tmp()
        type(starfile) :: sf
        type(string)   :: fname
        write(*,'(A)') 'test_complete_noop_without_tmp'
        fname = tmp('complete_noop.star')
        call sf%init(fname)
        call sf%complete()        ! nothing to rename
        call assert_false(file_exists(fname), 'complete noop: final file must not be created')
    end subroutine test_complete_noop_without_tmp

    ! 1.4 — complete overwrites a pre-existing final file
    subroutine test_complete_overwrites_existing()
        type(starfile) :: sf
        type(string)   :: fname
        write(*,'(A)') 'test_complete_overwrites_existing'
        fname = tmp('complete_overwrite.star')
        call write_dummy_file(fname%to_char())           ! pre-existing final
        call sf%init(fname)
        call write_dummy_file(fname%to_char()//'.tmp')   ! new staging content
        call sf%complete()
        call assert_true (file_exists(fname),                   'complete: final exists after overwrite')
        call assert_false(file_exists(fname%to_char()//'.tmp'), 'complete: .tmp gone after overwrite')
    end subroutine test_complete_overwrites_existing

    ! 1.5 — init removes a stale .tmp left from a previous run
    subroutine test_init_removes_stale_tmp()
        type(starfile) :: sf
        type(string)   :: fname
        write(*,'(A)') 'test_init_removes_stale_tmp'
        fname = tmp('init_stale_tmp.star')
        call write_dummy_file(fname%to_char()//'.tmp')   ! stale .tmp
        call sf%init(fname)
        call assert_false(file_exists(fname%to_char()//'.tmp'), 'init: stale .tmp must be removed')
    end subroutine test_init_removes_stale_tmp

    ! Helper: write a one-line dummy text file
    subroutine write_dummy_file(path)
        character(len=*), intent(in) :: path
        integer :: fh
        open(newunit=fh, file=path, status='replace', action='write')
        write(fh,'(A)') '# dummy'
        close(fh)
    end subroutine write_dummy_file

    !=======================================================================
    !  SECTION 2: write_optics_table
    !=======================================================================

    ! 2.1 — basic round-trip: one optics group, verify ogid is written
    subroutine test_write_optics_table_basic()
        type(starfile) :: sf
        type(oris)     :: o
        type(string)   :: fname
        integer        :: iv
        write(*,'(A)') 'test_write_optics_table_basic'
        fname = tmp('optics_basic.star')
        call sf%init(fname)
        call make_oris_n(o, 1, .false.)
        call o%set(1, 'ogid', 7.0)
        call sf%write_optics_table(o)
        call sf%complete()
        call assert_true(file_exists(fname), 'optics basic: output file exists')
        iv = readback_int(fname, EMDL_IMAGE_OPTICS_GROUP, 'optics')
        call assert_int(7, iv, 'optics basic: ogid == 7')
        call o%kill()
    end subroutine test_write_optics_table_basic

    ! 2.2 — rows with state == 0 must be excluded
    subroutine test_write_optics_table_state_filter()
        type(starfile) :: sf
        type(oris)     :: o
        type(string)   :: fname
        integer        :: n
        write(*,'(A)') 'test_write_optics_table_state_filter'
        fname = tmp('optics_state_filter.star')
        call sf%init(fname)
        call make_oris_n(o, 4, .false.)
        call o%set(1, 'ogid', 1.0) ; call o%set_state(1, 0)   ! dead
        call o%set(2, 'ogid', 2.0)                              ! alive
        call o%set(3, 'ogid', 3.0) ; call o%set_state(3, 0)   ! dead
        call o%set(4, 'ogid', 4.0)                              ! alive
        call sf%write_optics_table(o)
        call sf%complete()
        n = count_objects(fname, 'optics')
        call assert_int(2, n, 'optics state filter: 2 active rows written')
        call o%kill()
    end subroutine test_write_optics_table_state_filter

    ! 2.3 — empty oris → empty table (but file is created)
    subroutine test_write_optics_table_empty()
        type(starfile) :: sf
        type(oris)     :: o
        type(string)   :: fname
        write(*,'(A)') 'test_write_optics_table_empty'
        fname = tmp('optics_empty.star')
        call sf%init(fname)
        call make_oris_n(o, 0, .false.)
        call sf%write_optics_table(o)
        call sf%complete()
        call assert_true(file_exists(fname), 'optics empty: staging file promoted')
        call assert_int(1, count_objects(fname, 'optics'), 'optics empty: 1 objects') ! STAR files must have at least one object, even if it's empty
        call o%kill()
    end subroutine test_write_optics_table_empty

    ! 2.4 — verify all numeric fields survive a round-trip
    subroutine test_write_optics_table_all_fields()
        type(starfile) :: sf
        type(oris)     :: o
        type(string)   :: fname
        real(dp)       :: dv
        write(*,'(A)') 'test_write_optics_table_all_fields'
        fname = tmp('optics_all_fields.star')
        call sf%init(fname)
        call make_oris_n(o, 1, .false.)
        call o%set(1, 'ogid',  1.0)
        call o%set(1, 'pop',   42.0)
        call o%set(1, 'kv',    300.0)
        call o%set(1, 'smpd',  1.1)
        call o%set(1, 'cs',    2.7)
        call o%set(1, 'fraca', 0.1)
        call o%set(1, 'opcx',  100.0)
        call o%set(1, 'opcy',  200.0)
        call sf%write_optics_table(o)
        call sf%complete()
        dv = readback_double(fname, EMDL_CTF_VOLTAGE, 'optics')
        call assert_double(300.0_dp, dv, 'optics all_fields: kv==300')
        dv = readback_double(fname, EMDL_IMAGE_PIXEL_SIZE, 'optics')
        call assert_double(1.1_dp,   dv, 'optics all_fields: smpd==1.1', ulp_tol=1e6_dp)
        dv = readback_double(fname, EMDL_CTF_CS, 'optics')
        call assert_double(2.7_dp,   dv, 'optics all_fields: cs==2.7',   ulp_tol=1e6_dp)
        call assert_int(42, readback_int(fname, SMPL_OPTICS_POPULATION, 'optics'), &
                        'optics all_fields: pop==42')
        call o%kill()
    end subroutine test_write_optics_table_all_fields

    !=======================================================================
    !  SECTION 3: write_mics_table
    !=======================================================================

    ! 3.1 — basic: one micrograph, verify an int field survives round-trip
    subroutine test_write_mics_table_basic()
        type(starfile) :: sf
        type(oris)     :: o
        type(string)   :: fname
        integer        :: iv
        write(*,'(A)') 'test_write_mics_table_basic'
        fname = tmp('mics_basic.star')
        call sf%init(fname)
        call make_oris_n(o, 1, .false.)
        call o%set(1, 'ogid', 3.0)
        call o%set(1, 'xdim', 4096.0)
        call sf%write_mics_table(o)
        call sf%complete()
        call assert_true(file_exists(fname), 'mics basic: file exists')
        iv = readback_int(fname, EMDL_IMAGE_SIZE_X, 'micrographs')
        call assert_int(4096, iv, 'mics basic: xdim==4096')
        call o%kill()
    end subroutine test_write_mics_table_basic

    ! 3.2 — state == 0 rows excluded
    subroutine test_write_mics_table_state_filter()
        type(starfile) :: sf
        type(oris)     :: o
        type(string)   :: fname
        write(*,'(A)') 'test_write_mics_table_state_filter'
        fname = tmp('mics_state_filter.star')
        call sf%init(fname)
        call make_oris_n(o, 5, .false.)
        call o%set(1, 'ogid', 1.0) ; call o%set(1, 'xdim', 4096.0)
        call o%set(2, 'ogid', 1.0) ; call o%set(2, 'xdim', 4096.0)
        call o%set(3, 'ogid', 1.0) ; call o%set(3, 'xdim', 4096.0)
        call o%set(4, 'ogid', 1.0) ; call o%set(4, 'xdim', 4096.0)
        call o%set(5, 'ogid', 1.0) ; call o%set(5, 'xdim', 4096.0)
        call o%set_state(2, 0)
        call o%set_state(4, 0)
        call sf%write_mics_table(o)
        call sf%complete()
        call assert_int(3, count_objects(fname, 'micrographs'), &
                        'mics state filter: 3 active rows')
        call o%kill()
    end subroutine test_write_mics_table_state_filter

    ! 3.3 — empty input → empty table
    subroutine test_write_mics_table_empty()
        type(starfile) :: sf
        type(oris)     :: o
        type(string)   :: fname
        write(*,'(A)') 'test_write_mics_table_empty'
        fname = tmp('mics_empty.star')
        call sf%init(fname)
        call make_oris_n(o, 0, .false.)
        call sf%write_mics_table(o)
        call sf%complete()
        call assert_true(file_exists(fname), 'mics empty: file exists')
        call assert_int(1, count_objects(fname, 'micrographs'), 'mics empty: 1 object') ! STAR files must have at least one object, even if it's empty
        call o%kill()
    end subroutine test_write_mics_table_empty

    ! 3.4 — verify key double fields survive round-trip (defocus unit conversion)
    subroutine test_write_mics_table_all_fields()
        type(starfile) :: sf
        type(oris)     :: o
        type(string)   :: fname
        real(dp)       :: dv
        write(*,'(A)') 'test_write_mics_table_all_fields'
        fname = tmp('mics_all_fields.star')
        call sf%init(fname)
        call make_oris_n(o, 1, .false.)
        call o%set(1, 'dfx',    1.0)   ! 1 µm → 10000 Å in STAR
        call o%set(1, 'dfy',    2.0)
        call o%set(1, 'angast', 45.0)
        call o%set(1, 'nframes', 40.0)
        call sf%write_mics_table(o)
        call sf%complete()
        ! dfx stored as dfx/0.0001 → 1.0/0.0001 = 10000
        dv = readback_double(fname, EMDL_CTF_DEFOCUSU, 'micrographs')
        call assert_double(10000.0_dp, dv, 'mics all_fields: dfx unit conversion', ulp_tol=1e6_dp)
        call assert_int(40, readback_int(fname, EMDL_IMAGE_SIZE_Z, 'micrographs'), &
                        'mics all_fields: nframes==40')
        call o%kill()
    end subroutine test_write_mics_table_all_fields

    !=======================================================================
    !  SECTION 4: write_ptcl2D_table (serial)
    !=======================================================================

    ! 4.1 — basic: single particle, verify class written
    subroutine test_write_ptcl2D_table_basic()
        type(starfile) :: sf
        type(oris)     :: ptcl, stk
        type(string)   :: fname
        integer        :: iv
        write(*,'(A)') 'test_write_ptcl2D_table_basic'
        fname = tmp('ptcl2D_basic.star')
        call sf%init(fname)
        call make_oris_n(stk,  1, .false.)
        call stk%set(1, 'box', 128.0)
        call stk%set(1, 'stk', 'dummy.mrcs')
        call make_oris_n(ptcl, 1, .true.)
        call ptcl%set(1, 'stkind', 1.0)
        call ptcl%set(1, 'indstk', 1.0)
        call ptcl%set(1, 'class',  5.0)
        call sf%write_ptcl2D_table(ptcl, stk)
        call sf%complete()
        call assert_true(file_exists(fname), 'ptcl2D basic: file exists')
        iv = readback_int(fname, EMDL_PARTICLE_CLASS, 'particles')
        call assert_int(5, iv, 'ptcl2D basic: class==5')
        call ptcl%kill() ; call stk%kill()
    end subroutine test_write_ptcl2D_table_basic

    ! 4.2 — particles with state == 0 are excluded
    subroutine test_write_ptcl2D_table_state_filter()
        type(starfile) :: sf
        type(oris)     :: ptcl, stk
        type(string)   :: fname
        integer        :: i
        write(*,'(A)') 'test_write_ptcl2D_table_state_filter'
        fname = tmp('ptcl2D_state_filter.star')
        call sf%init(fname)
        call make_oris_n(stk, 1, .false.)
        call stk%set(1, 'box', 64.0)
        call stk%set(1, 'stk', 'a.mrcs')
        call make_oris_n(ptcl, 6, .true.)
        do i = 1, 6
            call ptcl%set(i, 'stkind', 1.0)
            call ptcl%set(i, 'indstk', real(i))
        end do
        call ptcl%set_state(2, 0)
        call ptcl%set_state(5, 0)
        call sf%write_ptcl2D_table(ptcl, stk)
        call sf%complete()
        call assert_int(4, count_objects(fname, 'particles'), &
                        'ptcl2D state filter: 4 active particles')
        call ptcl%kill() ; call stk%kill()
    end subroutine test_write_ptcl2D_table_state_filter

    ! 4.3 — stkind <= 0 rows are skipped
    subroutine test_write_ptcl2D_table_bad_stkind()
        type(starfile) :: sf
        type(oris)     :: ptcl, stk
        type(string)   :: fname
        write(*,'(A)') 'test_write_ptcl2D_table_bad_stkind'
        fname = tmp('ptcl2D_bad_stkind.star')
        call sf%init(fname)
        call make_oris_n(stk, 1, .false.)
        call stk%set(1, 'box', 64.0)
        call stk%set(1, 'stk', 'b.mrcs')
        call make_oris_n(ptcl, 3, .true.)
        call ptcl%set(1, 'stkind',  0.0)  ! invalid
        call ptcl%set(1, 'indstk',  1.0)
        call ptcl%set(2, 'stkind', -1.0)  ! invalid
        call ptcl%set(2, 'indstk',  1.0)
        call ptcl%set(3, 'stkind',  1.0)  ! valid
        call ptcl%set(3, 'indstk',  1.0)
        call sf%write_ptcl2D_table(ptcl, stk)
        call sf%complete()
        call assert_int(1, count_objects(fname, 'particles'), &
                        'ptcl2D bad stkind: only 1 valid particle written')
        call ptcl%kill() ; call stk%kill()
    end subroutine test_write_ptcl2D_table_bad_stkind

    ! 4.4 — indstk taken from ptcl oris when present
    subroutine test_write_ptcl2D_table_indstk_explicit()
        type(starfile) :: sf
        type(oris)     :: ptcl, stk
        type(string)   :: fname
        write(*,'(A)') 'test_write_ptcl2D_table_indstk_explicit'
        fname = tmp('ptcl2D_indstk_explicit.star')
        call sf%init(fname)
        call make_oris_n(stk, 1, .false.)
        call stk%set(1, 'box', 64.0)
        call stk%set(1, 'stk', 'c.mrcs')
        call make_oris_n(ptcl, 2, .true.)
        ! explicit indstk
        call ptcl%set(1, 'stkind', 1.0)
        call ptcl%set(1, 'indstk', 3.0)
        call ptcl%set(2, 'stkind', 1.0)
        call ptcl%set(2, 'indstk', 7.0)
        call sf%write_ptcl2D_table(ptcl, stk)
        call sf%complete()
        call assert_int(2, count_objects(fname, 'particles'), &
                        'ptcl2D indstk explicit: 2 particles written')
        call ptcl%kill() ; call stk%kill()
    end subroutine test_write_ptcl2D_table_indstk_explicit

    ! 4.5 — micrograph name written when mics_oris present and size matches
    subroutine test_write_ptcl2D_table_with_mics()
        type(starfile) :: sf
        type(oris)     :: ptcl, stk, mics
        type(string)   :: fname
        type(starfile_table_type)          :: tbl
        character(len=:), allocatable      :: sv
        logical :: ok
        write(*,'(A)') 'test_write_ptcl2D_table_with_mics'
        fname = tmp('ptcl2D_with_mics.star')
        call sf%init(fname)
        call make_oris_n(stk, 1, .false.)
        call stk%set(1, 'box', 64.0)
        call stk%set(1, 'stk', 'stack.mrcs')
        call make_oris_n(mics, 1, .false.)
        call mics%set(1, 'intg', 'micrograph_1.mrc')
        call make_oris_n(ptcl, 1, .true.)
        call ptcl%set(1, 'stkind', 1.0)
        call ptcl%set(1, 'indstk', 1.0)
        call sf%write_ptcl2D_table(ptcl, stk, mics_oris=mics)
        call sf%complete()
        ! Verify EMDL_MICROGRAPH_NAME present
        call starfile_table__new(tbl)
        call starfile_table__read(tbl, fname, 'particles')
        if( starfile_table__firstobject(tbl) >= 0 ) then
            ok = starfile_table__getValue_string(tbl, EMDL_MICROGRAPH_NAME, sv)
            call assert_true(ok, 'ptcl2D with_mics: micrograph name field present')
        else
            call assert_true(.false., 'ptcl2D with_mics: no first object found')
        end if
        call starfile_table__delete(tbl)
        call ptcl%kill() ; call stk%kill() ; call mics%kill()
    end subroutine test_write_ptcl2D_table_with_mics

    !=======================================================================
    !  SECTION 5: write_ptcl2D_table_parallel
    !=======================================================================

    ! 5.1 — parallel write with small dataset: file must exist and be valid
    subroutine test_write_ptcl2D_parallel_basic()
        type(starfile) :: sf
        type(oris)     :: ptcl, stk
        type(string)   :: fname
        write(*,'(A)') 'test_write_ptcl2D_parallel_basic'
        fname = tmp('ptcl2D_parallel_basic.star')
        call omp_set_num_threads(2)
        call sf%init(fname)
        call make_oris_n(stk, 1, .false.)
        call stk%set(1, 'box', 64.0)
        call stk%set(1, 'stk', 'par.mrcs')
        call make_oris_n(ptcl, 8, .true.)
        call fill_ptcl_stk(ptcl, stk, 8)
        call sf%write_ptcl2D_table_parallel(ptcl, stk)
        call sf%complete()
        call assert_true(file_exists(fname), 'ptcl2D parallel basic: file exists')
        call assert_int(8, count_objects(fname, 'particles'), &
                        'ptcl2D parallel basic: 8 particles')
        call ptcl%kill() ; call stk%kill()
    end subroutine test_write_ptcl2D_parallel_basic

    ! 5.2 — empty input → empty table (fast path)
    subroutine test_write_ptcl2D_parallel_empty()
        type(starfile) :: sf
        type(oris)     :: ptcl, stk
        type(string)   :: fname
        write(*,'(A)') 'test_write_ptcl2D_parallel_empty'
        fname = tmp('ptcl2D_parallel_empty.star')
        call sf%init(fname)
        call make_oris_n(stk, 1, .false.)
        call stk%set(1, 'box', 64.0)
        call make_oris_n(ptcl, 0, .true.)
        call sf%write_ptcl2D_table_parallel(ptcl, stk)
        call sf%complete()
        call assert_true(file_exists(fname), 'ptcl2D parallel empty: file exists')
        call assert_int(1, count_objects(fname, 'particles'), 'ptcl2D parallel empty: 1 object') ! STAR files must have at least one object, even if it's empty
        call ptcl%kill() ; call stk%kill()
    end subroutine test_write_ptcl2D_parallel_empty

    ! 5.3 — verify exact row count for a larger dataset
    subroutine test_write_ptcl2D_parallel_count()
        type(starfile) :: sf
        type(oris)     :: ptcl, stk
        type(string)   :: fname
        integer, parameter :: NPTCLS = 200
        write(*,'(A)') 'test_write_ptcl2D_parallel_count'
        fname = tmp('ptcl2D_parallel_count.star')
        call omp_set_num_threads(4)
        call sf%init(fname)
        call make_oris_n(stk, 1, .false.)
        call stk%set(1, 'box', 64.0)
        call stk%set(1, 'stk', 'big.mrcs')
        call make_oris_n(ptcl, NPTCLS, .true.)
        call fill_ptcl_stk(ptcl, stk, NPTCLS)
        call sf%write_ptcl2D_table_parallel(ptcl, stk)
        call sf%complete()
        call assert_int(NPTCLS, count_objects(fname, 'particles'), &
                        'ptcl2D parallel count: correct row count')
        call ptcl%kill() ; call stk%kill()
    end subroutine test_write_ptcl2D_parallel_count

    ! 5.4 — parallel output must match serial output row-for-row
    subroutine test_write_ptcl2D_parallel_vs_serial()
        type(starfile) :: sf_ser, sf_par
        type(oris)     :: ptcl, stk
        type(string)   :: fname_ser, fname_par
        integer, parameter :: N = 50
        integer :: nser, npar
        write(*,'(A)') 'test_write_ptcl2D_parallel_vs_serial'
        fname_ser = tmp('ptcl2D_pvs_serial.star')
        fname_par = tmp('ptcl2D_pvs_parallel.star')
        call omp_set_num_threads(3)
        call make_oris_n(stk, 1, .false.)
        call stk%set(1, 'box', 64.0)
        call stk%set(1, 'stk', 'pvs.mrcs')
        call make_oris_n(ptcl, N, .true.)
        call fill_ptcl_stk(ptcl, stk, N)
        call ptcl%set_state(10, 0)   ! one dead row — both paths must skip it
        call sf_ser%init(fname_ser)
        call sf_ser%write_ptcl2D_table(ptcl, stk)
        call sf_ser%complete()
        call sf_par%init(fname_par)
        call sf_par%write_ptcl2D_table_parallel(ptcl, stk)
        call sf_par%complete()
        nser = count_objects(fname_ser, 'particles')
        npar = count_objects(fname_par, 'particles')
        call assert_int(nser, npar, 'ptcl2D parallel vs serial: same object count')
        call ptcl%kill() ; call stk%kill()
    end subroutine test_write_ptcl2D_parallel_vs_serial

    ! Helper: populate ptcl/stk oris for N particles assigned to stk index 1
    subroutine fill_ptcl_stk(ptcl, stk, n)
        type(oris), intent(inout) :: ptcl, stk
        integer,    intent(in)    :: n
        integer :: i
        call stk%set(1, 'fromp', 1.0)
        call stk%set(1, 'top',   real(n))
        do i = 1, n
            call ptcl%set(i, 'stkind', 1.0)
            call ptcl%set(i, 'indstk', real(i))
            call ptcl%set(i, 'xpos',   real(mod(i, 64)))
            call ptcl%set(i, 'ypos',   real(mod(i, 64)))
            call ptcl%set(i, 'class',  real(mod(i, 10) + 1))
        end do
    end subroutine fill_ptcl_stk

    !=======================================================================
    !  SECTION 6: verbose flag
    !=======================================================================

    ! 6.1 — verbose=.true. must not change the number of written rows
    subroutine test_verbose_flag_no_content_change()
        type(starfile) :: sf_quiet, sf_verbose
        type(oris)     :: o
        type(string)   :: fq, fv
        integer        :: nq, nv
        write(*,'(A)') 'test_verbose_flag_no_content_change'
        fq = tmp('verbose_quiet.star')
        fv = tmp('verbose_verbose.star')
        call make_oris_n(o, 3, .false.)
        call o%set(1, 'ogid', 1.0)
        call o%set(2, 'ogid', 2.0)
        call o%set(3, 'ogid', 3.0)
        call sf_quiet%init(fq, verbose=.false.)
        call sf_quiet%write_optics_table(o)
        call sf_quiet%complete()
        call sf_verbose%init(fv, verbose=.true.)
        call sf_verbose%write_optics_table(o)
        call sf_verbose%complete()
        nq = count_objects(fq, 'optics')
        nv = count_objects(fv, 'optics')
        call assert_int(nq, nv, 'verbose flag: same row count regardless of verbose')
        call o%kill()
    end subroutine test_verbose_flag_no_content_change

end module simple_starfile_tester
