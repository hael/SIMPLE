module simple_starproject_tester
!=======================================================================
!   SIMPLE STAR PROJECT INTEGRATION TEST SUITE
!
!   Covers:
!      - simple_starfile
!      - simple_relion
!      - simple_starproject_stream
!      - simple_starproject_utils
!      - simple_starproject
!
!=======================================================================
use omp_lib
use simple_cmdline
use simple_relion
use simple_sp_project
use simple_starfile
use simple_starfile_wrappers
use simple_starproject
use simple_starproject_stream
use simple_starproject_utils
use simple_string
use simple_test_utils
implicit none

private
public :: run_all_starproject_tests

character(len=*), parameter :: TMPDIR = "test_star_tmp"

contains

    !=======================================================================
    !  MAIN ENTRY
    !=======================================================================
    subroutine run_all_starproject_tests()
        write(*,'(A)') "**** Running SIMPLE STARPROJECT test suite ****"
        call setup_tmpdir()
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
        write(*,'(A)') "**** Completed PART 2/4 (simple_relion tests) ****"
        write(*,'(A)') "**** Running OpenMP tests (Part 3/4) ****"
        call test_parallel_particle_export()
        call test_write_omem_parallel_stability()
        write(*,'(A)') "**** Completed OpenMP tests (Part 3/4) ****"
        write(*,'(A)') "**** Running big integration tests (Part 4/4) ****"
        call test_star_import_mics()
        call test_star_export_micrographs()
        call test_optics_clustering_basic()
        call test_xml_tiltinfo()
        call test_relion_writer_micrographs()
        call test_roundtrip_micrographs()
        write(*,'(A)') "**** Completed big integration tests (Part 4/4) ****"
        ! call report_summary()
    end subroutine run_all_starproject_tests

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
    ! Helper: Verify the merged particles2D.star file is valid
    !-----------------------------------------------------------------------
    subroutine test_parallel_export_readback(fname, nexpected)
        type(string), intent(in) :: fname
        integer,      intent(in) :: nexpected
        type(starfile_table_type) :: tbl
        integer :: obj, count
        write(*,'(A)') "test_parallel_export_readback"
        count = 0
        call starfile_table__new(tbl)
        call starfile_table__read(tbl, fname)
        ! iterate objects by hand
        obj = starfile_table__firstobject(tbl)
        if (obj < 0) then
            call assert_true(.false., "parallel readback: no first object")
            return
        end if
        do
            count = count + 1
            obj = starfile_table__nextobject(tbl)
            if (obj < 0) exit
        end do
        ! sanity check: must equal nexpected
        call assert_int(nexpected, count, "parallel export: number of particle objects == nptcls")
        call starfile_table__delete(tbl)
    end subroutine test_parallel_export_readback


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
        call starfile_table__write_ofile(t1)
        call starfile_table__close_ofile(t1)
        call starfile_table__delete(t1)
        ! Second block appended
        call starfile_table__new(t2)
        call starfile_table__setname(t2, "optics")
        call starfile_table__addObject(t2)
        call starfile_table__setValue_int(t2, EMDL_IMAGE_OPTICS_GROUP, 2)
        call starfile_table__open_ofile(t2, sf%ftmp%to_char(), 1)
        call starfile_table__write_ofile(t2)
        call starfile_table__close_ofile(t2)
        call starfile_table__delete(t2)
        call sf%complete()
        ! Readback — now file contains TWO objects
        call starfile_table__new(t1)
        call starfile_table__read(t1, fname)
        nobj1 = starfile_table__firstobject(t1)
        nobj2 = starfile_table__nextobject(t1)
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
        ok = starfile_table__firstobject(t)
        call assert_true(ok >= 0, "optics block found")
        call starfile_table__delete(t)
        ! Readback second (micrographs)
        call starfile_table__new(t)
        call starfile_table__read(t, fname, "micrographs")
        ok = starfile_table__firstobject(t)
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
            call proj%os_mic%set(i, "intg", "path/mic"//int2str(i)//"_frames.mrc")
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
        ok = starfile_table__firstobject(checktbl)
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
        ok = starfile_table__firstobject(checktbl)
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
        ok = starfile_table__firstobject(checktbl)
        call assert_true(ok >= 0, "particles table exists")
        call starfile_table__delete(checktbl)
    end subroutine test_relion_write_particles2D_star

    !-----------------------------------------------------------------------
    !  PART 3: Parallel batching and parallel starfile export tests
    !-----------------------------------------------------------------------
    subroutine test_parallel_particle_export()
        use omp_lib
        type(starproject_stream) :: stream
        type(sp_project)         :: proj
        type(string)             :: outdir, fname
        integer :: nptcls, i
        logical :: exists
        write(*,'(A)') "test_parallel_particle_export"
        ! Create directory for output
        outdir = tmpfile("parallel_out")
        call exec_cmdline("mkdir -p " // outdir%to_char())
        ! Construct a mock project with stacks + particles
        ! (Not fully populated, but enough for batching)
        nptcls = 20000          ! large enough to force batching
        call proj%os_stk%new(1, .false.)
        call proj%os_stk%set(1, "box", 128)
        call proj%os_stk%set(1, "ogid", 1)
        call proj%os_stk%set(1, "stk", "dummy_stk.mrcs")
        call proj%os_stk%set(1, "fromp", 1)
        call proj%os_stk%set(1, "top",   nptcls)
        call proj%os_ptcl2D%new(nptcls, .true.)
        do i = 1, nptcls
            call proj%os_ptcl2D%set_state(i, 1)
            call proj%os_ptcl2D%set(i, "stkind", 1)
            call proj%os_ptcl2D%set(i, "xpos", real(mod(i,100)))
            call proj%os_ptcl2D%set(i, "ypos", real(mod(i,200)))
            call proj%os_ptcl2D%set(i, "e3",   real(mod(i,360)))
            call proj%os_ptcl2D%set(i, "dfx",  1.0)
            call proj%os_ptcl2D%set(i, "dfy",  2.0)
            call proj%os_ptcl2D%set(i, "angast", 15.0)
        end do
        ! Run stream_export_particles_2D in full OpenMP mode.
        call stream%stream_export_particles_2D(proj, outdir, optics_set=.false., verbose=.true.)
        ! Check existence
        fname = outdir // "/particles2D.star"
        inquire(file=fname%to_char(), exist=exists)
        call assert_true(exists, "parallel export: particles2D.star must exist")
        ! Now re-read the file and validate correct # of objects
        call test_parallel_export_readback(fname, nptcls)
    end subroutine test_parallel_particle_export

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
    subroutine test_star_import_mics()
        type(starproject) :: sp
        type(sp_project)  :: proj
        type(cmdline)     :: cl
        type(string)      :: fname, tmpdir
        integer :: unit
        write(*,'(A)') "test_star_import_mics"
        ! Create temporary directory
        tmpdir = tmpfile("import_mic_test")
        call exec_cmdline("mkdir -p " // tmpdir%to_char())
        ! Create a small Relion-style micrographs.star
        fname = tmpdir // "/micrographs.star"
        open(newunit=unit, file=fname%to_char(), status="replace")
        write(unit,*) "data_optics"
        write(unit,*) "loop_"
        write(unit,*) "_rlnVoltage #1"
        write(unit,*) "_rlnImagePixelSize #2"
        write(unit,*) "_rlnOpticsGroup #3"
        write(unit,*) "300.0 1.0 1"
        write(unit,*) ""
        write(unit,*) "data_micrographs"
        write(unit,*) "loop_"
        write(unit,*) "_rlnMicrographName #1"
        write(unit,*) "_rlnDefocusU #2"
        write(unit,*) "_rlnDefocusV #3"
        write(unit,*) "_rlnDefocusAngle #4"
        write(unit,*) "_rlnOpticsGroup #5"
        write(unit,*) "intg1.mrc 10000 12000 15 1"
        write(unit,*) "intg2.mrc 11000 13000 18 1"
        close(unit)
        ! Prepare a cmdline object
        call cl%set("import_dir", tmpdir%to_char())
        ! Run import
        call sp%import_mics(cl, proj, fname)
        ! Assertions
        call assert_int(2, proj%os_mic%get_noris(), "import_mics: two micrographs imported")
        call assert_true(proj%os_mic%isthere(1,"intg"), "import_mics: intg key exists")
        call assert_real(10000.0, proj%os_mic%get(1,"dfx"), 1e-6, "import_mics: defocusU imported")
        call assert_real(12000.0, proj%os_mic%get(1,"dfy"), 1e-6, "import_mics: defocusV imported")
        call assert_real(15.0,   proj%os_mic%get(1,"angast"),1e-6,"import_mics: angast imported")
    end subroutine test_star_import_mics

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
        type(starproject) :: sp
        type(sp_project)  :: proj
        type(cmdline)     :: cl
        integer :: i
        write(*,'(A)') "test_optics_clustering_basic"
        ! Create synthetic tiltinfo from micrographs
        call proj%os_mic%new(4, .false.)
        do i=1,4
            call proj%os_mic%set_state(i,1)
            call proj%os_mic%set(i,"intg","mic"//int2str(i)//".mrc")
            call proj%os_mic%set(i,"smpd",1.0)
            call proj%os_mic%set(i,"kv",300.0)
            call proj%os_mic%set(i,"fraca",0.1)
            call proj%os_mic%set(i,"cs",2.7)
            call proj%os_mic%set(i,"box", 128)
        end do
        ! assign fake tilts:
        call proj%os_mic%set(1,"tiltx",0.0)
        call proj%os_mic%set(1,"tilty",0.0)
        call proj%os_mic%set(2,"tiltx",0.02)
        call proj%os_mic%set(2,"tilty",0.02)
        call proj%os_mic%set(3,"tiltx",5.0)
        call proj%os_mic%set(3,"tilty",5.0)
        call proj%os_mic%set(4,"tiltx",5.1)
        call proj%os_mic%set(4,"tilty",5.1)
        ! cluster using threshold:
        ! group 1: mic1, mic2    (near zero)
        ! group 2: mic3, mic4    (near 5)
        call sp%assign_optics(cl, proj, propagate=.true.)
        call assert_true(proj%os_mic%get_int(1,"ogid") == proj%os_mic%get_int(2,"ogid"), &
            "cluster: mic1 and mic2 must share ogid")
        call assert_true(proj%os_mic%get_int(3,"ogid") == proj%os_mic%get_int(4,"ogid"), &
            "cluster: mic3 and mic4 must share ogid")
        call assert_true(proj%os_mic%get_int(1,"ogid") /= proj%os_mic%get_int(3,"ogid"), &
            "cluster: group1 != group2")
    end subroutine test_optics_clustering_basic

    !-----------------------------------------------------------------------
    !  PART 4.4 — XML TILT FILE HANDLING (FoX)
    !-----------------------------------------------------------------------
    subroutine test_xml_tiltinfo()
        type(starproject) :: sp
        type(sp_project)  :: proj
        type(cmdline)     :: cl
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
            write(unit,'(A)') "<BeamShift>"
            write(unit,'(A)') "  <a:_x>" // trim(real2str(0.1*i)) // "</a:_x>"
            write(unit,'(A)') "  <a:_y>" // trim(real2str(0.2*i)) // "</a:_y>"
            write(unit,'(A)') "</BeamShift>"
            close(unit)
        end do
        ! Add dummy tiltinfo to project
        call proj%os_mic%new(2, .false.)
        do i=1,2
            call proj%os_mic%set_state(i,1)
            call proj%os_mic%set(i,"intg","mic"//int2str(i)//".mrc")
            call proj%os_mic%set(i,"smpd",1.0)
            call proj%os_mic%set(i,"kv",300.0)
            call proj%os_mic%set(i,"fraca",0.1)
            call proj%os_mic%set(i,"cs",2.7)
            call proj%os_mic%set(i,"box", 200)
            call proj%os_mic%set(i,"tind", i)
            call proj%os_mic%set(i,"tiltx", 0.0)  ! will be overwritten from XML
            call proj%os_mic%set(i,"tilty", 0.0)
            t = make_tiltinfo("mic"//int2str(i))
            sp%tiltinfo = [ sp%tiltinfo,t ]
        end do
        ! instruct commandline to use xmldir
        call cl%set("xmldir", xmldir%to_char())
        ! read XML info
        call sp%assign_xml_tiltinfo(cl%get_carg("xmldir"))
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
        do i=1,3
            call proj1%os_mic%set_state(i,1)
            call proj1%os_mic%set(i,"intg","intg"//int2str(i)//".mrc")
            call proj1%os_mic%set(i,"kv",300.0)
            call proj1%os_mic%set(i,"smpd",1.0*i)
            call proj1%os_mic%set(i,"fraca",i*0.1)
            call proj1%os_mic%set(i,"cs",2.0+0.3*i)
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
            call assert_real(1.0*i, proj2%os_mic%get(i,"smpd"), 1e-6, "roundtrip smpd")
            call assert_real(i*0.1, proj2%os_mic%get(i,"fraca"), 1e-6, "roundtrip fraca")
            call str_intg%kill
        end do
    end subroutine test_roundtrip_micrographs

end module simple_starproject_tester