!@descr: unit test routines for starproject
module simple_starproject_tester
!=======================================================================
!   SIMPLE STAR PROJECT INTEGRATION TEST SUITE
!
!   Covers:
!      - simple_starfile
!      - simple_relion
!      - simple_starproject_utils
!      - simple_starproject
!
!=======================================================================
use omp_lib
use simple_cmdline
use simple_defs
use simple_relion
use simple_sp_project
use simple_starfile
use simple_starfile_wrappers
use simple_starproject
use simple_starproject_utils
use simple_string
use simple_test_utils
use simple_type_defs, only: canonical_phshift
use simple_image,     only: image
use simple_oris,      only: oris
use simple_fileio,    only: simple_getcwd
implicit none

private
public :: run_all_starproject_tests

character(len=*), parameter :: TMPDIR = "test_star_tmp"
! the origins fixture: particle i alone in stack i, optics group i, at ORIGINS_SMPDS(i) A per pixel
real,             parameter :: ORIGINS_SMPDS(2) = [1.3, 0.9]
real,             parameter :: ORIGINS_PX(2,2) = reshape([2., -1., -3., 0.5], [2,2]) ! (x|y, particle)
real,             parameter :: ORIGINS_EUL(3,2) = reshape([10., 20., 30., 40., 50., 60.], [3,2])

contains

    !=======================================================================
    !  MAIN ENTRY
    !=======================================================================
    subroutine run_all_starproject_tests()
        integer :: nthr_orig
        write(*,'(A)') "**** Running SIMPLE STARPROJECT test suite ****"
        call setup_tmpdir()
        ! deterministic threading for the OpenMP part; restored at the end so the
        ! sub-suites that follow in the same process see the thread count they started with
        nthr_orig = omp_get_max_threads()
        call force_openmp_threads(4)
        write(*,'(A)') "**** Running PART 1/4 (starfile core tests) ****"
        call test_starfile_basic_init()
        call test_starfile_write_and_readback()
        call test_starfile_append_mode()
        call test_starfile_multiple_tables()
        write(*,'(A)') "**** Completed PART 1/4 (starfile core tests) ****"
        write(*,'(A)') "**** Running RELION interface tests (Part 2/4) ****"
        call test_relion_find_movienames()
        call test_relion_epu_tiltgroups()
        call test_relion_single_tiltgroup()
        call test_relion_allocate_opticsgroups()
        call test_relion_write_corrected_micrographs_star()
        call test_relion_write_micrographs_star()
        call test_relion_write_particles2D_star()
        call test_relion_phase_shift()
        call test_split_dataline_terminal_field()
        call test_resolve_stack_path()
        call test_relion_export_poses()
        write(*,'(A)') "**** Completed PART 2/4 (simple_relion tests) ****"
        write(*,'(A)') "**** Running OpenMP tests (Part 3/4) ****"
        call test_large_particle_table_export()
        call test_write_omem_parallel_stability()
        write(*,'(A)') "**** Completed OpenMP tests (Part 3/4) ****"
        write(*,'(A)') "**** Running big integration tests (Part 4/4) ****"
        call test_star_micrographs_readback()
        call test_star_export_micrographs()
        call test_optics_clustering_basic()
        call test_xml_tiltinfo()
        call test_relion_writer_micrographs()
        call test_roundtrip_micrographs()
        call test_import_origins_angstrom()
        call test_roundtrip_origins()
        call test_relion30_layout_refused()
        write(*,'(A)') "**** Completed big integration tests (Part 4/4) ****"
        call omp_set_num_threads(nthr_orig)
        ! the caller reports and fails the run (simple_test_utils counters are process-wide)
    end subroutine run_all_starproject_tests

    !> STAR stores rlnPhaseShift in degrees, SIMPLE the CTF phase in radians, with one factor for the
    !! micrograph, particle-2D and particle-3D maps. RELION does not constrain rlnPhaseShift to [0,180):
    !! its aberration fit can move a constant gamma offset into that column, and folding such a value at
    !! pi would negate the transfer function relative to what RELION computes from the same STAR file
    subroutine test_relion_phase_shift()
        real, parameter :: PHASE_DEG = 45., BEYOND_PI_DEG = 250.
        real :: phase_internal, phase_exported
        write(*,'(A)') 'test_relion_phase_shift'
        phase_internal = PHASE_DEG * RELION_PHASE_DEG2RAD
        call assert_real(PI/4., phase_internal, 1.e-6, 'RELION phase-shift import converts degrees to SIMPLE radians')
        phase_exported = phase_internal / RELION_PHASE_DEG2RAD
        call assert_real(PHASE_DEG, phase_exported, 1.e-6, 'RELION phase-shift export converts SIMPLE radians to degrees')
        phase_internal = canonical_phshift(BEYOND_PI_DEG * RELION_PHASE_DEG2RAD)
        call assert_real(BEYOND_PI_DEG * RELION_PHASE_DEG2RAD, phase_internal, 1.e-6, &
            &'a RELION phase shift beyond 180 degrees is not folded on import')
        phase_exported = phase_internal / RELION_PHASE_DEG2RAD
        call assert_real(BEYOND_PI_DEG, phase_exported, 1.e-4, &
            &'a phase shift beyond 180 degrees round-trips back to RELION degrees')
    end subroutine test_relion_phase_shift

    subroutine test_split_dataline_terminal_field()
        type(string) :: line, splitline(3)
        logical      :: fieldmask(3)
        write(*,'(A)') 'test_split_dataline_terminal_field'
        line = 'particle.mrc 1 29'
        call split_dataline(line, splitline)
        call assert_char('particle.mrc', splitline(1)%to_char(), 'split_dataline: first field without trailing whitespace')
        call assert_char('1',            splitline(2)%to_char(), 'split_dataline: middle field without trailing whitespace')
        call assert_char('29',           splitline(3)%to_char(), 'split_dataline: multi-character final field')
        line = 'particle.mrc 2 4'
        call split_dataline(line, splitline)
        call assert_char('4', splitline(3)%to_char(), 'split_dataline: one-character final field')
        line = 'particle.mrc'//achar(9)//'3 22 '
        call split_dataline(line, splitline)
        call assert_char('particle.mrc', splitline(1)%to_char(), 'split_dataline: tab delimiter')
        call assert_char('22',           splitline(3)%to_char(), 'split_dataline: trailing whitespace')
        splitline(2) = 'not copied'
        fieldmask = [.true., .false., .true.]
        line = 'particle.mrc 4 29'
        call split_dataline(line, splitline, fieldmask)
        call assert_char('not copied', splitline(2)%to_char(), 'split_dataline: field mask skips unneeded strings')
        call assert_char('29',         splitline(3)%to_char(), 'split_dataline: field mask preserves final field')
        call line%kill
        call splitline(:)%kill
    end subroutine test_split_dataline_terminal_field

    subroutine test_resolve_stack_path()
        type(string) :: path
        logical      :: exists
        write(*,'(A)') 'test_resolve_stack_path'
        call write_textfile(TMPDIR//'/stack.mrc', 'dummy')
        path = TMPDIR//'/missing/stack.mrc'
        call resolve_stack_path(path, exists)
        call assert_true(exists, 'resolve_stack_path: repairs a path with an extra directory')
        call assert_char(TMPDIR//'/stack.mrc', path%to_char(), 'resolve_stack_path: returns repaired path')
        path = TMPDIR//'/missing/first.mrc'
        call resolve_stack_path(path, exists)
        call assert_false(exists, 'resolve_stack_path: first missing stack remains missing')
        path = TMPDIR//'/missing/second.mrc'
        call resolve_stack_path(path, exists)
        call assert_false(exists, 'resolve_stack_path: consecutive missing stacks remain bounds-safe')
        call path%kill
    end subroutine test_resolve_stack_path

    !=======================================================================
    !  TEST ENVIRONMENT HELPERS
    !=======================================================================
    subroutine setup_tmpdir()
        logical :: ex
        inquire(file=TMPDIR, exist=ex)
        if ( ex ) then
            call delete_directory(TMPDIR)
        end if
        call create_directory(TMPDIR)
    end subroutine setup_tmpdir

    subroutine create_directory(dirname)
        character(len=*), intent(in) :: dirname
        call exec_cmdline("mkdir -p " // dirname)
    end subroutine create_directory

    subroutine delete_directory(dirname)
        character(len=*), intent(in) :: dirname
        call exec_cmdline("rm -rf " // dirname)
    end subroutine delete_directory

    ! Force deterministic OpenMP threading
    subroutine force_openmp_threads(n)
        integer, intent(in) :: n
        call omp_set_num_threads(n)
        write(*,*) "  [OpenMP] Forced thread count =", n
    end subroutine force_openmp_threads

    ! Utility: create path inside tmpdir
    function tmpfile(fname) result(out)
        type(string) :: out
        character(len=*), intent(in) :: fname
        out = TMPDIR // "/" // fname
    end function tmpfile

    ! Utility: write a simple text file
    subroutine write_textfile(path, contents)
        character(len=*), intent(in) :: path, contents
        integer :: fh, ios
        open(newunit=fh, file=path, status='replace', action='write', iostat=ios)
        call assert_int(0, ios, "write_textfile: open")
        write(fh,'(A)') contents
        close(fh)
    end subroutine write_textfile

    ! Helper constructor for tilt_info
    function make_tiltinfo(basename) result(t)
        character(len=*), intent(in) :: basename
        type(tilt_info) :: t
        t%basename = basename
        t%initialtiltgroupid = 1
        t%finaltiltgroupid   = 1
        t%tiltx = 0.0
        t%tilty = 0.0
        t%smpd  = 1.0
        t%kv    = 300.0
        t%fraca = 0.1
        t%cs    = 2.7
        t%box   = 200
    end function make_tiltinfo

     !-----------------------------------------------------------------------
    ! Helper: Verify the particles2D.star file is valid
    !-----------------------------------------------------------------------
    subroutine test_particle_export_readback(fname, nexpected)
        type(string), intent(in) :: fname
        integer,      intent(in) :: nexpected
        type(starfile_table_type) :: tbl
        integer :: obj, count
        write(*,'(A)') "test_particle_export_readback"
        count = 0
        call starfile_table__new(tbl)
        call starfile_table__read(tbl, fname, "particles")
        ! iterate objects by hand
         obj = int(starfile_table__firstobject(tbl))
        if (obj < 0) then
            call assert_true(.false., "particle readback: no first object")
            return
        end if
        do
            count = count + 1
             obj = int(starfile_table__nextobject(tbl))
            if (obj < 0) exit
        end do
        ! sanity check: must equal nexpected
        call assert_int(nexpected, count, "particle export: number of objects == nptcls")
        call starfile_table__delete(tbl)
    end subroutine test_particle_export_readback


    !=======================================================================
    !  SECTION 1: BASIC TESTS FOR simple_starfile
    !=======================================================================

    !-----------------------------------------------------------------------
    ! 1.1: starfile %init + %complete
    !-----------------------------------------------------------------------
    subroutine test_starfile_basic_init()
        type(starfile) :: sf
        type(string)   :: fname
        write(*,'(A)') "test_starfile_basic_init"
        fname = tmpfile("basic_init.star")
        call sf%init(fname)
        call assert_true(.not. file_exists(fname), "starfile%init should not yet create output file")
        call assert_true(file_exists(fname%to_char()//".tmp") .eqv. .false., "no tmp file yet")
        ! Write a tiny optics block
        call write_textfile(fname%to_char()//".tmp", "dummy")
        call sf%complete()
        call assert_true(file_exists(fname), "starfile%complete should move .tmp → final")
    end subroutine test_starfile_basic_init

    !-----------------------------------------------------------------------
    ! 1.2: Write single-table optics block and read it back
    !-----------------------------------------------------------------------
    subroutine test_starfile_write_and_readback()
        type(starfile)           :: sf
        type(starfile_table_type):: tbl
        type(string)             :: fname
        real(dp)                 :: dv
        integer                  :: iv
        logical                  :: ok
        write(*,'(A)') "test_starfile_write_and_readback"
        fname = tmpfile("optics_write_test.star")
        call sf%init(fname, verbose=.false.)
        ! build single-object optics table
        call starfile_table__new(tbl)
        call starfile_table__setIsList(tbl, .false.)
        call starfile_table__setname(tbl, "optics")
        call starfile_table__addObject(tbl)
        call starfile_table__setValue_int(tbl, EMDL_IMAGE_OPTICS_GROUP, 3)
        call starfile_table__setValue_double(tbl, EMDL_CTF_VOLTAGE, 300.0_dp)
        call starfile_table__open_ofile(tbl, sf%ftmp%to_char(), 0)
        call starfile_table__write_ofile(tbl)
        call starfile_table__close_ofile(tbl)
        call starfile_table__delete(tbl)
        call sf%complete()
        ! read it back
        call starfile_table__new(tbl)
        call starfile_table__read(tbl, fname)
        ok = starfile_table__getValue_int(tbl, EMDL_IMAGE_OPTICS_GROUP, iv)
        call assert_true(ok, "readback: ogid present?")
        call assert_int(3, iv, "readback: ogid == 3")
        ok = starfile_table__getValue_double(tbl, EMDL_CTF_VOLTAGE, dv)
        call assert_true(ok, "readback: voltage present?")
        call assert_double(300.0_dp, dv, "readback voltage")
        call starfile_table__delete(tbl)
    end subroutine test_starfile_write_and_readback

    !-----------------------------------------------------------------------
    ! 1.3: Append mode test
    !-----------------------------------------------------------------------
    subroutine test_starfile_append_mode()
        type(starfile)            :: sf
        type(starfile_table_type) :: t1, t2
        type(string)              :: fname
        integer                   :: nobj1, nobj2
        write(*,'(A)') "test_starfile_append_mode"
        fname = tmpfile("append_mode.star")
        call sf%init(fname)
        ! First block
        call starfile_table__new(t1)
        call starfile_table__setname(t1, "optics")
        call starfile_table__addObject(t1)
        call starfile_table__setValue_int(t1, EMDL_IMAGE_OPTICS_GROUP, 1)
        call starfile_table__open_ofile(t1, sf%ftmp%to_char(), 0)
        call starfile_table__write_ofile(t1, tableend=.false.)
        call starfile_table__close_ofile(t1)
        call starfile_table__delete(t1)
        ! Second block appended
        call starfile_table__new(t2)
        call starfile_table__setname(t2, "optics")
        call starfile_table__addObject(t2)
        call starfile_table__setValue_int(t2, EMDL_IMAGE_OPTICS_GROUP, 2)
        call starfile_table__open_ofile(t2, sf%ftmp%to_char(), 1)
        call starfile_table__write_ofile(t2, ignoreheader=.true.)
        call starfile_table__close_ofile(t2)
        call starfile_table__delete(t2)
        call sf%complete()
        ! Readback — now file contains TWO objects
        call starfile_table__new(t1)
        call starfile_table__read(t1, fname)
         nobj1 = int(starfile_table__firstobject(t1))
         nobj2 = int(starfile_table__nextobject(t1))
        call assert_true(nobj1 >= 0, "first object exists")
        call assert_true(nobj2 >= 0, "second object exists")
        call starfile_table__delete(t1)
    end subroutine test_starfile_append_mode

    !-----------------------------------------------------------------------
    ! 1.4: Multi-table STAR file test
    !-----------------------------------------------------------------------
    subroutine test_starfile_multiple_tables()
        type(starfile)            :: sf
        type(starfile_table_type) :: t
        type(string) :: fname
        integer :: ok
        write(*,'(A)') "test_starfile_multiple_tables"
        fname = tmpfile("multi_table.star")
        call sf%init(fname)
        ! Write OPTICS table
        call starfile_table__new(t)
        call starfile_table__setname(t, "optics")
        call starfile_table__addObject(t)
        call starfile_table__setValue_int(t, EMDL_IMAGE_OPTICS_GROUP, 5)
        call starfile_table__open_ofile(t, sf%ftmp%to_char(), 0)
        call starfile_table__write_ofile(t)
        call starfile_table__close_ofile(t)
        call starfile_table__delete(t)
        ! Write MICROGRAPHS table
        call starfile_table__new(t)
        call starfile_table__setname(t, "micrographs")
        call starfile_table__addObject(t)
        call starfile_table__setValue_int(t, EMDL_IMAGE_SIZE_X, 4096)
        call starfile_table__open_ofile(t, sf%ftmp%to_char(), 1)
        call starfile_table__write_ofile(t)
        call starfile_table__close_ofile(t)
        call starfile_table__delete(t)
        call sf%complete()
        ! Readback first (optics)
        call starfile_table__new(t)
        call starfile_table__read(t, fname, "optics")
         ok = int(starfile_table__firstobject(t))
        call assert_true(ok >= 0, "optics block found")
        call starfile_table__delete(t)
        ! Readback second (micrographs)
        call starfile_table__new(t)
        call starfile_table__read(t, fname, "micrographs")
         ok = int(starfile_table__firstobject(t))
        call assert_true(ok >= 0, "micrographs block found")
        call starfile_table__delete(t)
    end subroutine test_starfile_multiple_tables

    !-----------------------------------------------------------------------
    ! 2.1  find_movienames
    !-----------------------------------------------------------------------
    subroutine test_relion_find_movienames()
        type(relion_project) :: rp
        type(sp_project)     :: proj
        type(cmdline)        :: cl
        integer              :: i
        write(*,'(A)') "test_relion_find_movienames"
        ! create 3 fake micrographs with intg keys
        call proj%os_mic%new(3, is_ptcl=.false.)
        do i=1,3
            call proj%os_mic%set(i, "intg", "path/mic"//int2str(i)//"_intg.mrc")
            call proj%os_mic%set_state(i, 1)
        end do
        call rp%find_movienames(cl, proj)
        call assert_true(allocated(rp%movienames), "movienames allocated")
        call assert_int(3, size(rp%movienames), "3 movienames detected")
        call assert_char("mic1", rp%movienames(1)%to_char(), "first basename ok")
    end subroutine test_relion_find_movienames

    !-----------------------------------------------------------------------
    ! 2.2  generate_epu_tiltgroups
    !-----------------------------------------------------------------------
    subroutine test_relion_epu_tiltgroups()
        type(relion_project) :: rp
        type(sp_project)     :: proj
        type(cmdline)        :: cl
        integer              :: i
        write(*,'(A)') "test_relion_epu_tiltgroups"
        call proj%os_mic%new(4, .false.)
        ! Fake filenames containing EPU “Data_XYZ”
        call proj%os_mic%set(1,"intg",".../FoilHole_123_Data_AA_frames.mrc")
        call proj%os_mic%set(2,"intg",".../FoilHole_123_Data_AA_frames.mrc")
        call proj%os_mic%set(3,"intg",".../FoilHole_123_Data_BB_frames.mrc")
        call proj%os_mic%set(4,"intg",".../FoilHole_123_Data_CC_frames.mrc")
        do i=1,4
            call proj%os_mic%set_state(i,1)
        end do
        call rp%find_movienames(cl, proj)
        call rp%generate_epu_tiltgroups(cl, proj)
        call assert_int(3, rp%opticsgroups, "3 distinct EPU tilt groups")
        call assert_int(1, rp%moviegroup(1), "AA group id=1")
        call assert_int(1, rp%moviegroup(2), "AA group id=1 again")
        call assert_int(2, rp%moviegroup(3), "BB group id=2")
        call assert_int(3, rp%moviegroup(4), "CC group id=3")
    end subroutine test_relion_epu_tiltgroups

    !-----------------------------------------------------------------------
    ! 2.3  generate_single_tiltgroup
    !-----------------------------------------------------------------------
    subroutine test_relion_single_tiltgroup()
        type(relion_project) :: rp
        type(sp_project)     :: proj
        type(cmdline)        :: cl
        integer              :: i
        write(*,'(A)') "test_relion_single_tiltgroup"
        call proj%os_mic%new(3, .false.)
        do i=1,3
            call proj%os_mic%set(i,"intg","abc"//int2str(i)//".mrc")
            call proj%os_mic%set_state(i,1)
        end do
        call rp%find_movienames(cl,proj)
        call rp%generate_single_tiltgroup(cl,proj)
        call assert_int(1, rp%opticsgroups, "single group")
        do i=1,3
            call assert_int(1, rp%moviegroup(i), "all assigned to group 1")
        end do
    end subroutine test_relion_single_tiltgroup

    !-----------------------------------------------------------------------
    ! 2.4  allocate_opticsgroups
    !-----------------------------------------------------------------------
    subroutine test_relion_allocate_opticsgroups()
        type(relion_project) :: rp
        type(sp_project)     :: proj
        type(cmdline)        :: cl
        integer              :: i
        write(*,'(A)') "test_relion_allocate_opticsgroups"
        call proj%os_mic%new(3, .false.)
        do i=1,3
            call proj%os_mic%set(i,"intg","mic"//int2str(i)//"_frames.mrc")
            call proj%os_mic%set_state(i,1)
        end do
        call cl%set("tiltgroupmax", 1.0)   ! force splitting after more than 1 mic/group
        call rp%find_movienames(cl,proj)
        call rp%generate_single_tiltgroup(cl,proj)
        call rp%allocate_opticsgroups(cl,proj)
        ! Now 3 mics, tiltgroupmax=1 ⇒ expect splitting:
        call assert_true(rp%opticsgroups >= 2, "splitting opticsgroups triggered")
    end subroutine test_relion_allocate_opticsgroups

    !-----------------------------------------------------------------------
    ! 2.5 write_corrected_micrographs_star
    !-----------------------------------------------------------------------
    subroutine test_relion_write_corrected_micrographs_star()
        type(relion_project)      :: rp
        type(sp_project)          :: proj
        type(cmdline)             :: cl
        type(starfile_table_type) :: checktbl
        type(string)              :: fname
        integer :: i, ok
        write(*,'(A)') "test_relion_write_corrected_micrographs_star"
        call proj%os_mic%new(2, .false.)
        do i=1,2
            call proj%os_mic%set(i,"intg","mic"//int2str(i)//".mrc")
            call proj%os_mic%set(i,"movie","mov"//int2str(i)//".eer")
            call proj%os_mic%set(i,"smpd",1.5)
            call proj%os_mic%set(i,"kv",300)
            call proj%os_mic%set(i,"cs",2.7)
            call proj%os_mic%set(i,"fraca",0.1)
            call proj%os_mic%set(i,"opticsgroup",real(i))
            call proj%os_mic%set_state(i,1)
        end do
        rp%opticsgroups = 2
        call cl%set("optics_offset", 0.0)
        fname = string("micrographs_corrected.star")
        if (file_exists(fname)) call del_file(fname)
        call rp%write_corrected_micrographs_star(cl,proj)
        call assert_true(file_exists(fname), "corrected micrographs star created")
        call starfile_table__new(checktbl)
        call starfile_table__read(checktbl,fname,"optics")
         ok = int(starfile_table__firstobject(checktbl))
        call assert_true(ok >= 0,"optics block exists")

        call starfile_table__delete(checktbl)
    end subroutine test_relion_write_corrected_micrographs_star

    !-----------------------------------------------------------------------
    ! 2.6 write_micrographs_star
    !-----------------------------------------------------------------------
    subroutine test_relion_write_micrographs_star()
        type(relion_project)      :: rp
        type(sp_project)          :: proj
        type(cmdline)             :: cl
        type(starfile_table_type) :: checktbl
        type(string)              :: fname
        integer :: i, ok
        write(*,'(A)') "test_relion_write_micrographs_star"
        call proj%os_mic%new(2, .false.)
        do i=1,2
            call proj%os_mic%set(i,"intg","mic"//int2str(i)//".mrc")
            call proj%os_mic%set(i,"forctf","for"//int2str(i)//".mrc")
            call proj%os_mic%set(i,"pspec","pspec"//int2str(i)//".mrc")
            call proj%os_mic%set(i,"box",256)
            call proj%os_mic%set(i,"smpd",1.2)
            call proj%os_mic%set(i,"kv",300)
            call proj%os_mic%set(i,"cs",2.7)
            call proj%os_mic%set(i,"dfx",1.0)
            call proj%os_mic%set(i,"dfy",2.0)
            call proj%os_mic%set(i,"angast",45.0)
            call proj%os_mic%set(i,"opticsgroup",1.0)
            call proj%os_mic%set_state(i,1)
        end do
        rp%opticsgroups = 1
        call cl%set("optics_offset",0.0)
        fname = string("micrographs.star")
        if(file_exists(fname)) call del_file(fname)
        call rp%write_micrographs_star(cl,proj)
        call assert_true(file_exists(fname),"micrographs.star created")
        call starfile_table__new(checktbl)
        call starfile_table__read(checktbl,fname,"micrographs")
         ok = int(starfile_table__firstobject(checktbl))
        call assert_true(ok >= 0, "micrographs table exists")
        call starfile_table__delete(checktbl)
    end subroutine test_relion_write_micrographs_star

    !-----------------------------------------------------------------------
    ! 2.7 write_particles2D_star
    !-----------------------------------------------------------------------
    subroutine test_relion_write_particles2D_star()
        type(relion_project)      :: rp
        type(sp_project)          :: proj
        type(cmdline)             :: cl
        type(starfile_table_type) :: checktbl
        type(string)              :: fname
        integer :: i, ok
        write(*,'(A)') "test_relion_write_particles2D_star"
        call proj%os_ptcl2D%new(3,.true.)
        call proj%os_stk%new(2,.false.)
        ! stk definitions
        call proj%os_stk%set(1,"stk","stack1.mrcs")
        call proj%os_stk%set(1,"box",256)
        call proj%os_stk%set_state(1,1)
        call proj%os_stk%set(2,"stk","stack2.mrcs")
        call proj%os_stk%set(2,"box",256)
        call proj%os_stk%set_state(2,1)
        ! assign stkind
        call proj%os_ptcl2D%set(1,"stkind",1)
        call proj%os_ptcl2D%set(2,"stkind",1)
        call proj%os_ptcl2D%set(3,"stkind",2)
        do i=1,3
            call proj%os_ptcl2D%set(i,"xpos",10.0)
            call proj%os_ptcl2D%set(i,"ypos",20.0)
            call proj%os_ptcl2D%set(i,"dfx",1.0)
            call proj%os_ptcl2D%set(i,"dfy",2.0)
            call proj%os_ptcl2D%set(i,"angast",30.0)
            call proj%os_ptcl2D%set(i,"opticsgroup",1.0)
            call proj%os_ptcl2D%set_state(i,1)
        end do
        rp%opticsgroups = 1
        call cl%set("reliongroups",0.0)
        call cl%set("optics_offset",0.0)
        fname = string("particles2D.star")
        if(file_exists(fname)) call del_file(fname)
        call rp%write_particles2D_star(cl,proj)
        call assert_true(file_exists(fname),"particles2D.star created")
        call starfile_table__new(checktbl)
        call starfile_table__read(checktbl,fname,"particles")
         ok = int(starfile_table__firstobject(checktbl))
        call assert_true(ok >= 0, "particles table exists")
        call starfile_table__delete(checktbl)
    end subroutine test_relion_write_particles2D_star

    !-----------------------------------------------------------------------
    !  PART 3: Large-table and OpenMP serialization tests
    !-----------------------------------------------------------------------
    subroutine test_large_particle_table_export()
        type(starfile_table_type) :: tbl
        type(string) :: outdir, fname
        integer :: nptcls, i
        logical :: exists
        write(*,'(A)') "test_large_particle_table_export"
        ! Create directory for output
        outdir = tmpfile("large_table_out")
        call exec_cmdline("mkdir -p " // outdir%to_char())
        ! Exercise table growth and serialization with a realistic row count.
        nptcls = 20000
        call starfile_table__new(tbl)
        call starfile_table__setname(tbl, "particles")
        do i = 1, nptcls
            call starfile_table__addObject(tbl)
            call starfile_table__setValue_int(tbl, EMDL_IMAGE_OPTICS_GROUP, 1)
            call starfile_table__setValue_double(tbl, EMDL_IMAGE_COORD_X, real(mod(i,100), dp))
            call starfile_table__setValue_double(tbl, EMDL_IMAGE_COORD_Y, real(mod(i,200), dp))
            call starfile_table__setValue_double(tbl, EMDL_ORIENT_PSI, real(mod(i,360), dp))
            call starfile_table__setValue_string(tbl, EMDL_IMAGE_NAME, int2str(i)//"@dummy_stk.mrcs")
        end do
        fname = outdir // "/particles2D.star"
        call starfile_table__open_ofile(tbl, fname%to_char(), 0)
        call starfile_table__write_ofile(tbl)
        call starfile_table__close_ofile(tbl)
        call starfile_table__delete(tbl)
        inquire(file=fname%to_char(), exist=exists)
        call assert_true(exists, "large particle table: particles2D.star must exist")
        call test_particle_export_readback(fname, nptcls)
    end subroutine test_large_particle_table_export

    !-----------------------------------------------------------------------
    ! Additional test:
    !   verify that write_omem produces stable buffer ordering
    !-----------------------------------------------------------------------
    subroutine test_write_omem_parallel_stability()
        use omp_lib
        type(starfile_table_type) :: tbl
        type(string) :: tmp
        character(len=:), allocatable :: s
        integer :: len1
        integer :: i
        write(*,'(A)') "test_write_omem_parallel_stability"
        tmp = tmpfile("test_write_omem.star")
        ! Build a big table for write_omem testing
        call starfile_table__new(tbl)
        call starfile_table__setname(tbl, "particles")
        do i = 1, 2000
            call starfile_table__addObject(tbl)
            call starfile_table__setValue_int(tbl, EMDL_PARTICLE_CLASS, mod(i,10))
            call starfile_table__setValue_double(tbl, EMDL_ORIENT_PSI, real(mod(i,360),dp))
        end do
        ! Parallel call write_omem (we simulate repeated calls)-
        !$omp parallel default(shared)
        !$omp single
        call starfile_table__write_omem(tbl, s, len1, ignoreheader=.false.)
        !$omp end single
        !$omp end parallel
        call assert_true(len1 > 0, "write_omem produced non-zero buffer size")
        ! optional consistency check: must contain "loop_" once
        call assert_true(index(s, "loop_") > 0, "write_omem buffer contains STAR loop block")
        call starfile_table__delete(tbl)
    end subroutine test_write_omem_parallel_stability

    !-----------------------------------------------------------------------
    !  PART 4.1 — STAR IMPORT TEST
    !-----------------------------------------------------------------------
    subroutine test_star_micrographs_readback()
        type(starfile_table_type)     :: tbl
        type(string)                  :: fname, tmpdir
        character(len=:), allocatable :: micname
        real(dp)                      :: value
        integer                       :: unit, obj
        logical                       :: ok
        write(*,'(A)') "test_star_micrographs_readback"
        ! Create temporary directory
        tmpdir = tmpfile("import_mic_test")
        call exec_cmdline("mkdir -p " // tmpdir%to_char())
        ! Create a small Relion-style micrographs.star
        fname = tmpdir // "/micrographs.star"
        open(newunit=unit, file=fname%to_char(), status="replace")
        write(unit,'(A)') "data_optics"
        write(unit,'(A)') "loop_"
        write(unit,'(A)') "_rlnOpticsGroupName #1"
        write(unit,'(A)') "_rlnVoltage #2"
        write(unit,'(A)') "_rlnSphericalAberration #3"
        write(unit,'(A)') "_rlnAmplitudeContrast #4"
        write(unit,'(A)') "_rlnImagePixelSize #5"
        write(unit,'(A)') "_rlnImageSize #6"
        write(unit,'(A)') "_rlnOpticsGroup #7"
        write(unit,'(A)') "opticsGroup1 300.0 2.7 0.1 1.0 200 1"
        write(unit,'(A)') ""
        write(unit,'(A)') "data_micrographs"
        write(unit,'(A)') "loop_"
        write(unit,'(A)') "_rlnMicrographName #1"
        write(unit,'(A)') "_rlnDefocusU #2"
        write(unit,'(A)') "_rlnDefocusV #3"
        write(unit,'(A)') "_rlnDefocusAngle #4"
        write(unit,'(A)') "_rlnPhaseShift #5"
        write(unit,'(A)') "_rlnOpticsGroup #6"
        write(unit,'(A)') "intg1.mrc 10000 12000 15 45 1"
        write(unit,'(A)') "intg2.mrc 11000 13000 18 90 1"
        close(unit)
        ! Read the micrograph block directly.  Project-level import and optics
        ! remapping have their own lifecycle requirements and are not part of
        ! the STAR parser contract exercised by this suite.
        call starfile_table__new(tbl)
        call starfile_table__read(tbl, fname, "micrographs")
         obj = int(starfile_table__firstobject(tbl))
        call assert_true(obj >= 0, "micrograph readback: first row exists")
        ok = starfile_table__getValue_string(tbl, EMDL_MICROGRAPH_NAME, micname)
        call assert_true(ok, "micrograph readback: name present")
        call assert_true(micname == "intg1.mrc", "micrograph readback: first name")
        ok = starfile_table__getValue_double(tbl, EMDL_CTF_DEFOCUSU, value)
        call assert_true(ok, "micrograph readback: defocus U present")
        call assert_double(10000.0_dp, value, "micrograph readback: defocus U")
        ok = starfile_table__getValue_double(tbl, EMDL_CTF_PHASESHIFT, value)
        call assert_true(ok, "micrograph readback: phase shift present")
        call assert_double(45.0_dp, value, "micrograph readback: phase shift")
         obj = int(starfile_table__nextobject(tbl))
        call assert_true(obj >= 0, "micrograph readback: second row exists")
         obj = int(starfile_table__nextobject(tbl))
        call assert_true(obj < 0, "micrograph readback: exactly two rows")
        call starfile_table__delete(tbl)
    end subroutine test_star_micrographs_readback

    !-----------------------------------------------------------------------
    !  PART 4.2 — STAR EXPORT TEST
    !-----------------------------------------------------------------------
    subroutine test_star_export_micrographs()
        type(starproject) :: sp
        type(sp_project)  :: proj
        type(string)      :: out
        logical :: exists
        write(*,'(A)') "test_star_export_micrographs"
        out = tmpfile("export_mic_test")
        call exec_cmdline("mkdir -p " // out%to_char())
        ! Create dummy micrographs in project
        call proj%os_mic%new(2, .false.)
        call proj%os_mic%set_state(1,1)
        call proj%os_mic%set(1,"intg","intg1.mrc")
        call proj%os_mic%set(1,"kv",300.0)
        call proj%os_mic%set(1,"smpd",1.0)
        call proj%os_mic%set(1,"fraca",0.07)
        call proj%os_mic%set(1,"cs",2.7)
        call proj%os_mic%set(1,"dfx",10000.0)
        call proj%os_mic%set(1,"dfy",12000.0)
        call proj%os_mic%set_state(2,1)
        call proj%os_mic%set(2,"intg","intg2.mrc")
        call proj%os_mic%set(2,"kv",300.0)
        call proj%os_mic%set(2,"smpd",1.0)
        call proj%os_mic%set(2,"fraca",0.07)
        call proj%os_mic%set(2,"cs",2.7)
        call proj%os_mic%set(2,"dfx",11000.0)
        call proj%os_mic%set(2,"dfy",13000.0)
        ! Export STAR
        call sp%export_mics(proj)
        inquire(file="micrographs.star", exist=exists)
        call assert_true(exists, "export_micrographs: file exists after export")
    end subroutine test_star_export_micrographs

    !-----------------------------------------------------------------------
    !  PART 4.3 — OPTICS + TILT CLUSTERING TESTS
    !-----------------------------------------------------------------------
    subroutine test_optics_clustering_basic()
        integer              :: labels(4)
        integer, allocatable :: populations(:)
        real                 :: tilts(4,2)
        real,    allocatable :: centroids(:,:)
        write(*,'(A)') "test_optics_clustering_basic"
        ! Exercise the numerical clustering used by optics assignment without
        ! coupling this unit test to starproject filename/lifecycle behavior.
        tilts(:,1) = [0.0, 0.02, 5.0, 5.1]
        tilts(:,2) = [0.0, 0.02, 5.0, 5.1]
        call h_clust(tilts, 0.2, labels, centroids, populations)
        ! group 1: mic1, mic2    (near zero)
        ! group 2: mic3, mic4    (near 5)
        call assert_int(2, size(populations), "cluster: two optics groups")
        call assert_true(labels(1) == labels(2), &
            "cluster: mic1 and mic2 must share ogid")
        call assert_true(labels(3) == labels(4), &
            "cluster: mic3 and mic4 must share ogid")
        call assert_true(labels(1) /= labels(3), &
            "cluster: group1 != group2")
        deallocate(centroids, populations)
    end subroutine test_optics_clustering_basic

    !-----------------------------------------------------------------------
    !  PART 4.4 — XML TILT FILE HANDLING (FoX)
    !-----------------------------------------------------------------------
    subroutine test_xml_tiltinfo()
        type(starproject) :: sp
        type(string)      :: xmldir, basename
        type(tilt_info)   :: t
        type(string)      :: fn
        integer :: unit, i
        write(*,'(A)') "test_xml_tiltinfo"
        xmldir = tmpfile("xml_test")
        call exec_cmdline("mkdir -p " // xmldir%to_char())
        ! Create minimal XML files:
        ! <BeamShift><a:_x>...</a:_x><a:_y>...</a:_y></BeamShift>
        do i = 1,2
            basename = "mic"//int2str(i)
            fn = xmldir // "/" // basename // ".xml"
            open(newunit=unit, file=fn%to_char(), status="replace")
            write(unit,'(A)') '<BeamShift xmlns:a="urn:simple-test">'
            write(unit,'(A)') "  <a:_x>" // trim(real2str(0.1*i)) // "</a:_x>"
            write(unit,'(A)') "  <a:_y>" // trim(real2str(0.2*i)) // "</a:_y>"
            write(unit,'(A)') "</BeamShift>"
            close(unit)
        end do
        ! Add dummy tiltinfo to project
        allocate(sp%tiltinfo(0))
        do i=1,2
            t = make_tiltinfo("mic"//int2str(i))
            sp%tiltinfo = [ sp%tiltinfo,t ]
        end do
        ! read XML info
        call sp%assign_xml_tiltinfo(xmldir)
        ! Check values
        call assert_real(0.1, sp%tiltinfo(1)%tiltx, 1e-6, "xml tilt x #1")
        call assert_real(0.2, sp%tiltinfo(1)%tilty, 1e-6, "xml tilt y #1")
        call assert_real(0.2, sp%tiltinfo(2)%tiltx, 1e-6, "xml tilt x #2")
        call assert_real(0.4, sp%tiltinfo(2)%tilty, 1e-6, "xml tilt y #2")
    end subroutine test_xml_tiltinfo

    !-----------------------------------------------------------------------
    !  PART 4.5 — RELION COMPATIBLE STAR WRITER (simple_relion)
    !-----------------------------------------------------------------------
    subroutine test_relion_writer_micrographs()
        type(relion_project) :: rp
        type(sp_project)     :: proj
        type(cmdline)        :: cl
        logical :: exists
        integer :: i
        write(*,'(A)') "test_relion_writer_micrographs"
        ! Create micrographs in project
        call proj%os_mic%new(2,.false.)
        do i=1,2
            call proj%os_mic%set_state(i,1)
            call proj%os_mic%set(i,"intg","mic"//int2str(i)//".mrc")
            call proj%os_mic%set(i,"movie","mic"//int2str(i)//".eer")
            call proj%os_mic%set(i,"smpd",1.0)
            call proj%os_mic%set(i,"kv",300.0)
            call proj%os_mic%set(i,"cs",2.7)
            call proj%os_mic%set(i,"fraca",0.08)
            call proj%os_mic%set(i,"opticsgroup", 1)
            call proj%os_mic%set(i,"state",1)
        end do
        ! A cmdline object with optics_offset=0
        call cl%set("optics_offset", 0.0)
        ! run writer
        rp%opticsgroups = 1
        call rp%write_micrographs_star(cl, proj)
        inquire(file="micrographs.star", exist=exists)
        call assert_true(exists, "relion writer: micrographs.star created")
    end subroutine test_relion_writer_micrographs

    !-----------------------------------------------------------------------
    !  PART 4.6 — END-TO-END ROUNDTRIP
    !-----------------------------------------------------------------------
    subroutine test_roundtrip_micrographs()
        type(starproject)     :: sp1, sp2
        type(sp_project)      :: proj1, proj2
        type(cmdline)         :: cl
        type(string)          :: fname, outdir, str_intg
        logical :: exists
        integer :: i
        write(*,'(A)') "test_roundtrip_micrographs"
        outdir = tmpfile("rtrip_mic")
        call exec_cmdline("mkdir -p " // outdir%to_char())
        ! Create data
        call proj1%os_mic%new(3,.false.)
        call proj1%os_optics%new(3,.false.)
        do i=1,3
            call proj1%os_mic%set_state(i,1)
            call proj1%os_mic%set(i,"intg","intg"//int2str(i)//".mrc")
            call proj1%os_mic%set(i,"ogid",i)
            call proj1%os_mic%set(i,"kv",300.0)
            call proj1%os_mic%set(i,"smpd",1.0*i)
            call proj1%os_mic%set(i,"fraca",i*0.1)
            call proj1%os_mic%set(i,"cs",2.0+0.3*i)
            call proj1%os_mic%set(i,"dfx",10000.+100.*i)
            call proj1%os_mic%set(i,"dfy",11000.+100.*i)
            call proj1%os_mic%set(i,"angast",10.*i)
            call proj1%os_mic%set(i,"phshift",deg2rad(20.*i))
            call proj1%os_optics%set_state(i,1)
            call proj1%os_optics%set(i,"ogid",i)
            call proj1%os_optics%set(i,"ogname","opticsgroup"//int2str(i))
            call proj1%os_optics%set(i,"kv",300.0)
            call proj1%os_optics%set(i,"smpd",1.0*i)
            call proj1%os_optics%set(i,"fraca",i*0.1)
            call proj1%os_optics%set(i,"cs",2.0+0.3*i)
            call proj1%os_optics%set(i,"opcx",0.0)
        end do
        ! Export
        call sp1%export_mics(proj1)
        fname = "micrographs.star"
        inquire(file=fname%to_char(), exist=exists)
        call assert_true(exists, "roundtrip: micrographs.star written")
        ! Re-import into proj2
        call cl%set("import_dir", ".")
        call sp2%import_mics(cl, proj2, fname)
        ! Validate all fields match
        do i=1,3
            str_intg = proj2%os_mic%get_str(i,"intg")
            call assert_char("intg"//int2str(i)//".mrc", str_intg%to_char(), "roundtrip intg")
            call assert_real(1.0*i, proj2%os_optics%get(i,"smpd"), 1e-6, "roundtrip optics smpd")
            call assert_real(i*0.1, proj2%os_mic%get(i,"fraca"), 1e-6, "roundtrip fraca")
            call assert_real(deg2rad(20.*i), proj2%os_mic%get(i,"phshift"), 1e-5, &
                &"roundtrip phase shift preserves SIMPLE radians through RELION degrees")
            call str_intg%kill
        end do
    end subroutine test_roundtrip_micrographs

    !-----------------------------------------------------------------------
    !  PARTICLE ORIGINS: Angstrom in STAR, pixels of the particle's stack in SIMPLE
    !-----------------------------------------------------------------------

    ! the absolute fixture directory and its two stacks of one 16x16 image, at ORIGINS_SMPDS
    subroutine make_origins_stacks( dir, stknames )
        type(string), intent(out) :: dir, stknames(2)
        type(image)  :: img
        type(string) :: cwd
        integer      :: i
        call simple_getcwd(cwd)
        dir = cwd//'/'//TMPDIR//'/origins'
        call create_directory(dir%to_char())
        do i = 1, 2
            stknames(i) = dir//'/ptcls_'//int2str(i)//'.mrcs'
            call img%new([16,16,1], ORIGINS_SMPDS(i))
            call img%gauran(0., 1.)
            call img%write(stknames(i), 1)
            call img%kill
        end do
    end subroutine make_origins_stacks

    ! a RELION 3.1 particle STAR file of the fixture: poses, and origins in Angstrom
    subroutine write_origins_star( fname, stknames )
        class(string), intent(in) :: fname, stknames(2)
        integer :: funit, i
        open(newunit=funit, file=fname%to_char(), status='replace', action='write')
        write(funit,'(A)') 'data_optics'
        write(funit,'(A)') 'loop_'
        write(funit,'(A)') '_rlnOpticsGroupName #1'
        write(funit,'(A)') '_rlnOpticsGroup #2'
        write(funit,'(A)') '_rlnVoltage #3'
        write(funit,'(A)') '_rlnSphericalAberration #4'
        write(funit,'(A)') '_rlnAmplitudeContrast #5'
        write(funit,'(A)') '_rlnImagePixelSize #6'
        write(funit,'(A)') '_rlnImageSize #7'
        do i = 1, 2
            write(funit,'(A,I0,A,I0,A,F8.4,A)') 'opticsGroup', i, ' ', i, ' 300.0 2.7 0.1 ', ORIGINS_SMPDS(i), ' 16'
        end do
        write(funit,'(A)') ''
        write(funit,'(A)') 'data_particles'
        write(funit,'(A)') 'loop_'
        write(funit,'(A)') '_rlnImageName #1'
        write(funit,'(A)') '_rlnAngleRot #2'
        write(funit,'(A)') '_rlnAngleTilt #3'
        write(funit,'(A)') '_rlnAnglePsi #4'
        write(funit,'(A)') '_rlnOriginXAngst #5'
        write(funit,'(A)') '_rlnOriginYAngst #6'
        write(funit,'(A)') '_rlnOpticsGroup #7'
        do i = 1, 2
            write(funit,'(2A,5F12.6,A,I0)') '1@', stknames(i)%to_char(), ORIGINS_EUL(:,i), &
                &ORIGINS_SMPDS(i) * ORIGINS_PX(:,i), ' ', i
        end do
        close(funit)
    end subroutine write_origins_star

    ! both particle fields hold the fixture's pixel origins; the 3D field its Euler angles when l3D
    subroutine check_origins_px( proj, l3D, label )
        class(sp_project), intent(inout) :: proj
        logical,           intent(in)    :: l3D
        character(len=*),  intent(in)    :: label
        real    :: err, eerr
        integer :: i
        call assert_int(2, proj%os_ptcl2D%get_noris(), label//': two particles in the 2D field')
        call assert_int(2, proj%os_ptcl3D%get_noris(), label//': two particles in the 3D field')
        if( proj%os_ptcl2D%get_noris() /= 2 .or. proj%os_ptcl3D%get_noris() /= 2 ) return
        err  = 0.
        eerr = 0.
        do i = 1, 2
            err = max(err, abs(proj%os_ptcl2D%get(i, 'x') - ORIGINS_PX(1,i)), abs(proj%os_ptcl2D%get(i, 'y') - ORIGINS_PX(2,i)))
            err = max(err, abs(proj%os_ptcl3D%get(i, 'x') - ORIGINS_PX(1,i)), abs(proj%os_ptcl3D%get(i, 'y') - ORIGINS_PX(2,i)))
            if( l3D ) eerr = max(eerr, maxval(abs(proj%os_ptcl3D%get_euler(i) - ORIGINS_EUL(:,i))))
        end do
        call assert_real(0., err, 1.e-5, label//': Angstrom origins import to the known pixel shifts')
        if( l3D ) call assert_real(0., eerr, 1.e-4, label//': Euler angles import unchanged')
    end subroutine check_origins_px

    ! every particle import path converts the origins once, by the pixel size of the particle's stack:
    ! import_starproject with import_type ptcl2D and ptcl3D (the commander's sequence), import_particles
    subroutine test_import_origins_angstrom()
        type(starproject) :: sp
        type(sp_project)  :: proj
        type(cmdline)     :: cl
        type(string)      :: dir, stknames(2), fname
        integer           :: ipath
        write(*,'(A)') 'test_import_origins_angstrom'
        call make_origins_stacks(dir, stknames)
        fname = dir//'/particles31.star'
        call write_origins_star(fname, stknames)
        call cl%set('import_dir', dir%to_char())
        do ipath = 1, 3
            select case(ipath)
                case(1)
                    call sp%import_ptcls2D(cl, proj, fname)
                    call proj%os_ptcl3D%copy(proj%os_ptcl2D, is_ptcl=.true.)
                    call sp%finalise_ptcl_import(proj)
                    call check_origins_px(proj, .false., 'import_starproject ptcl2D')
                case(2)
                    call sp%import_ptcls3D(cl, proj, fname)
                    call proj%os_ptcl2D%copy(proj%os_ptcl3D, is_ptcl=.true.)
                    call sp%finalise_ptcl_import(proj)
                    call check_origins_px(proj, .true., 'import_starproject ptcl3D')
                case(3)
                    call sp%import_ptcls(fname, proj, string('no'))
                    call check_origins_px(proj, .true., 'import_particles')
            end select
            call sp%kill
            call proj%kill
        end do
    end subroutine test_import_origins_angstrom

    ! pixels -> STAR (the writer of export_starproject) -> pixels is the identity, two stacks at
    ! different pixel sizes
    subroutine test_roundtrip_origins()
        type(starfile)    :: sf
        type(starproject) :: sp
        type(sp_project)  :: proj, proj2
        type(cmdline)     :: cl
        type(string)      :: dir, stknames(2), fname
        integer           :: i
        write(*,'(A)') 'test_roundtrip_origins'
        call make_origins_stacks(dir, stknames)
        fname = dir//'/roundtrip.star'
        call proj%os_optics%new(2, is_ptcl=.false.)
        call proj%os_stk%new(2, is_ptcl=.false.)
        call proj%os_ptcl2D%new(2, is_ptcl=.true.)
        do i = 1, 2
            call proj%os_optics%set(i, 'ogid',   i)
            call proj%os_optics%set(i, 'ogname', 'opticsgroup'//int2str(i))
            call proj%os_optics%set(i, 'kv',     300.)
            call proj%os_optics%set(i, 'cs',     2.7)
            call proj%os_optics%set(i, 'fraca',  0.1)
            call proj%os_optics%set(i, 'smpd',   ORIGINS_SMPDS(i))
            call proj%os_optics%set(i, 'box',    16)
            call proj%os_optics%set_state(i, 1)
            call proj%os_stk%set(i, 'stk',   stknames(i))
            call proj%os_stk%set(i, 'box',   16)
            call proj%os_stk%set(i, 'smpd',  ORIGINS_SMPDS(i))
            call proj%os_stk%set(i, 'ogid',  i)
            call proj%os_stk%set(i, 'fromp', i)
            call proj%os_stk%set(i, 'top',   i)
            call proj%os_stk%set_state(i, 1)
            call proj%os_ptcl2D%set(i, 'stkind', i)
            call proj%os_ptcl2D%set(i, 'indstk', 1)
            call proj%os_ptcl2D%set(i, 'ogid',   i)
            call proj%os_ptcl2D%set(i, 'x',      ORIGINS_PX(1,i))
            call proj%os_ptcl2D%set(i, 'y',      ORIGINS_PX(2,i))
            call proj%os_ptcl2D%set_state(i, 1)
        end do
        call sf%init(fname)
        call sf%write_optics_table(proj%os_optics)
        call sf%write_ptcl2D_table(proj%os_ptcl2D, proj%os_stk)
        call sf%complete()
        call cl%set('import_dir', dir%to_char())
        call sp%import_ptcls2D(cl, proj2, fname)
        call proj2%os_ptcl3D%copy(proj2%os_ptcl2D, is_ptcl=.true.)
        call sp%finalise_ptcl_import(proj2)
        call check_origins_px(proj2, .false., 'pixels to STAR to pixels')
        call sp%kill
        call proj%kill
        call proj2%kill
    end subroutine test_roundtrip_origins

    ! the rows' values of label in the particles table of fname; found when every row has it
    subroutine read_particle_doubles( fname, label, vals, found )
        class(string), intent(in)  :: fname
        integer,       intent(in)  :: label
        real(dp),      intent(out) :: vals(:)
        logical,       intent(out) :: found
        type(starfile_table_type) :: tbl
        real(dp) :: v
        integer  :: obj, n
        vals  = 0._dp
        found = .true.
        n     = 0
        call starfile_table__new(tbl)
        call starfile_table__read(tbl, fname, 'particles')
        obj = int(starfile_table__firstobject(tbl))
        do while( obj >= 0 )
            n = n + 1
            if( n > size(vals) ) exit
            if( starfile_table__getValue_double(tbl, label, v) )then
                vals(n) = v
            else
                found = .false.
            endif
            obj = int(starfile_table__nextobject(tbl))
        end do
        if( n /= size(vals) ) found = .false.
        call starfile_table__delete(tbl)
    end subroutine read_particle_doubles

    ! export_relion writes the 3D poses and origins when the 3D field holds poses, else the
    ! 2D field's in-plane angle and origins; origins in Angstrom
    subroutine test_relion_export_poses()
        real, parameter :: X2D(2) = [1.5, -2.5], Y2D(2) = [0.5, 3.0], PSI2D(2) = [15., 25.]
        type(relion_project) :: rp
        type(sp_project)     :: proj
        type(cmdline)        :: cl
        type(string)         :: fname
        real(dp)             :: vals(2)
        logical              :: found
        integer              :: i, icase
        write(*,'(A)') 'test_relion_export_poses'
        fname = string('particles2D.star')
        call cl%set('reliongroups',  0.0)
        call cl%set('optics_offset', 0.0)
        do icase = 1, 2
            call proj%os_stk%new(1, is_ptcl=.false.)
            call proj%os_stk%set(1, 'stk',  'ptcls_from_poses_intg.mrc')
            call proj%os_stk%set(1, 'box',  16)
            call proj%os_stk%set(1, 'smpd', ORIGINS_SMPDS(1))
            call proj%os_stk%set_state(1, 1)
            call proj%os_ptcl2D%new(2, is_ptcl=.true.)
            call proj%os_ptcl3D%new(2, is_ptcl=.true.)
            do i = 1, 2
                call proj%os_ptcl2D%set(i, 'stkind', 1)
                call proj%os_ptcl2D%set(i, 'e3',     PSI2D(i))
                call proj%os_ptcl2D%set(i, 'x',      X2D(i))
                call proj%os_ptcl2D%set(i, 'y',      Y2D(i))
                call proj%os_ptcl2D%set_state(i, 1)
                call proj%os_ptcl3D%set(i, 'stkind', 1)
                call proj%os_ptcl3D%set_state(i, 1)
                if( icase == 1 )then
                    call proj%os_ptcl3D%set_euler(i, ORIGINS_EUL(:,i))
                    call proj%os_ptcl3D%set(i, 'x', ORIGINS_PX(1,i))
                    call proj%os_ptcl3D%set(i, 'y', ORIGINS_PX(2,i))
                endif
            end do
            rp%opticsgroups = 1
            if( file_exists(fname) ) call del_file(fname)
            call rp%write_particles2D_star(cl, proj)
            if( icase == 1 )then
                call read_particle_doubles(fname, EMDL_ORIENT_ROT, vals, found)
                call assert_true(found, 'export_relion 3D poses: rlnAngleRot written')
                call assert_real(0., real(maxval(abs(vals - ORIGINS_EUL(1,:)))), 1.e-4, 'export_relion 3D poses: rot = e1')
                call read_particle_doubles(fname, EMDL_ORIENT_TILT, vals, found)
                call assert_true(found, 'export_relion 3D poses: rlnAngleTilt written')
                call assert_real(0., real(maxval(abs(vals - ORIGINS_EUL(2,:)))), 1.e-4, 'export_relion 3D poses: tilt = e2')
                call read_particle_doubles(fname, EMDL_ORIENT_PSI, vals, found)
                call assert_real(0., real(maxval(abs(vals - ORIGINS_EUL(3,:)))), 1.e-4, 'export_relion 3D poses: psi = e3')
                call read_particle_doubles(fname, EMDL_ORIENT_ORIGIN_X_ANGSTROM, vals, found)
                call assert_true(found, 'export_relion 3D poses: rlnOriginXAngst written')
                call assert_real(0., real(maxval(abs(vals - ORIGINS_SMPDS(1) * ORIGINS_PX(1,:)))), 1.e-4, &
                    &'export_relion 3D poses: origin x in Angstrom')
                call read_particle_doubles(fname, EMDL_ORIENT_ORIGIN_Y_ANGSTROM, vals, found)
                call assert_real(0., real(maxval(abs(vals - ORIGINS_SMPDS(1) * ORIGINS_PX(2,:)))), 1.e-4, &
                    &'export_relion 3D poses: origin y in Angstrom')
            else
                call read_particle_doubles(fname, EMDL_ORIENT_ROT, vals, found)
                call assert_false(found, 'export_relion 2D poses: no rlnAngleRot')
                call read_particle_doubles(fname, EMDL_ORIENT_PSI, vals, found)
                call assert_true(found, 'export_relion 2D poses: rlnAnglePsi written')
                call assert_real(0., real(maxval(abs(vals - PSI2D))), 1.e-4, 'export_relion 2D poses: psi = 2D e3')
                call read_particle_doubles(fname, EMDL_ORIENT_ORIGIN_X_ANGSTROM, vals, found)
                call assert_real(0., real(maxval(abs(vals - ORIGINS_SMPDS(1) * X2D))), 1.e-4, &
                    &'export_relion 2D poses: origin x in Angstrom')
                call read_particle_doubles(fname, EMDL_ORIENT_ORIGIN_Y_ANGSTROM, vals, found)
                call assert_real(0., real(maxval(abs(vals - ORIGINS_SMPDS(1) * Y2D))), 1.e-4, &
                    &'export_relion 2D poses: origin y in Angstrom')
            endif
            call proj%kill
        end do
    end subroutine test_relion_export_poses

    ! a STAR file without a data_optics block (the RELION 3.0 layout) is refused with a message
    ! that names the missing block; a RELION 3.1 file is not
    subroutine test_relion30_layout_refused()
        type(starproject) :: sp
        type(string)      :: dir, stknames(2), fname, msg
        integer           :: funit, iblock
        write(*,'(A)') 'test_relion30_layout_refused'
        call make_origins_stacks(dir, stknames)
        do iblock = 1, 2
            fname = dir//'/relion30_'//int2str(iblock)//'.star'
            open(newunit=funit, file=fname%to_char(), status='replace', action='write')
            ! RELION 3.0 writes one unnamed block; a particles block alone is refused alike
            if( iblock == 1 )then
                write(funit,'(A)') 'data_'
            else
                write(funit,'(A)') 'data_particles'
            endif
            write(funit,'(A)') 'loop_'
            write(funit,'(A)') '_rlnImageName #1'
            write(funit,'(A)') '_rlnAnglePsi #2'
            write(funit,'(A)') '_rlnOriginXAngst #3'
            write(funit,'(A)') '_rlnImagePixelSize #4'
            write(funit,'(3A)') '1@', stknames(1)%to_char(), ' 10.0 2.6 1.3'
            close(funit)
            sp%starfile%filename = fname
            call sp%read_layout(msg)
            call assert_true(index(msg%to_char(), 'no data_optics block') > 0, 'RELION 3.0 layout: the refusal names data_optics')
            call assert_true(index(msg%to_char(), 'RELION 3.1 or later') > 0, 'RELION 3.0 layout: the refusal names RELION 3.1')
            call sp%kill
        end do
        fname = dir//'/relion31.star'
        call write_origins_star(fname, stknames)
        sp%starfile%filename = fname
        call sp%read_layout(msg)
        call assert_int(0, msg%strlen(), 'a RELION 3.1 file is not refused')
        call sp%kill
    end subroutine test_relion30_layout_refused

end module simple_starproject_tester
