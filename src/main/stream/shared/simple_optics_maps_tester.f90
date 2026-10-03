!@descr: unit tests for the optics-map publication protocol (simple_optics_maps)
module simple_optics_maps_tester
use simple_test_utils
use simple_defs_fname,   only: OPTICS_MAP_PREFIX, TXT_EXT, METADATA_EXT
use simple_string,       only: string
use simple_string_utils, only: int2str
use simple_fileio,       only: file_exists, simple_getcwd, simple_rmdir
use simple_syslib,       only: simple_mkdir, get_process_id, dir_exists
use simple_sp_project,   only: sp_project
use simple_optics_maps,  only: publish_optics_map, import_latest_optics_map, latest_optics_map_id,&
                               &copy_project_with_optics_map
implicit none
private
public :: run_all_optics_maps_tests

integer, parameter :: NKEEP = 5

contains

    subroutine run_all_optics_maps_tests()
        write(*,'(A)') '**** running all optics map tests ****'
        call test_publish_keeps_newest()
        call test_import_latest()
        call test_empty_directory()
        call test_copy_with_optics_map()
    end subroutine run_all_optics_maps_tests

    !> seven publications with nkeep=5 leave maps 3 to 7
    subroutine test_publish_keeps_newest()
        type(sp_project) :: spproj
        type(string)     :: dir
        integer          :: id
        write(*,'(A)') 'test_publish_keeps_newest'
        dir = new_test_dir('publish')
        call make_grouped_mics(spproj)
        do id = 1,7
            call publish_optics_map(spproj, dir, id, NKEEP)
        enddo
        do id = 1,2
            call assert_false(map_exists(dir, id), 'map '//int2str(id)//' was pruned')
        enddo
        do id = 3,7
            call assert_true(map_exists(dir, id), 'map '//int2str(id)//' is kept')
        enddo
        call assert_int(7, latest_optics_map_id(dir), 'the newest map id')
        call spproj%kill
        call simple_rmdir(dir)
    end subroutine test_publish_keeps_newest

    !> a project with the same import indices and no groups takes the groups of the newest map
    subroutine test_import_latest()
        type(sp_project) :: spproj, reader
        type(string)     :: dir
        integer          :: id, imic
        write(*,'(A)') 'test_import_latest'
        dir = new_test_dir('import')
        call make_grouped_mics(spproj)
        call publish_optics_map(spproj, dir, 1, NKEEP)
        call spproj%os_mic%set(1, 'ogid', 2.0) ! map 2 moves micrograph 1 to group 2
        call publish_optics_map(spproj, dir, 2, NKEEP)
        call reader%os_mic%new(3, is_ptcl=.false.)
        do imic = 1,3
            call reader%os_mic%set(imic, 'importind', real(imic))
        enddo
        id = import_latest_optics_map(reader, dir)
        call assert_int(2, id, 'the newest map is imported')
        call assert_int(2, reader%os_mic%get_int(1, 'ogid'), 'micrograph 1 takes its group from the newest map')
        call assert_int(2, reader%os_mic%get_int(2, 'ogid'), 'micrograph 2 takes its group')
        call assert_int(2, reader%os_mic%get_int(3, 'ogid'), 'micrograph 3 takes its group')
        call spproj%kill
        call reader%kill
        call simple_rmdir(dir)
    end subroutine test_import_latest

    subroutine test_empty_directory()
        type(sp_project) :: reader
        type(string)     :: dir
        write(*,'(A)') 'test_empty_directory'
        dir = new_test_dir('empty')
        call reader%os_mic%new(1, is_ptcl=.false.)
        call reader%os_mic%set(1, 'importind', 1.0)
        call assert_int(0, latest_optics_map_id(dir),             'no map id in an empty directory')
        call assert_int(0, import_latest_optics_map(reader, dir), 'nothing is imported from an empty directory')
        call reader%kill
        call simple_rmdir(dir)
    end subroutine test_empty_directory

    !> a project copied with a map takes its groups on every segment; without a map it is copied
    !! unchanged; both copies are renamed into place from their .tmp
    subroutine test_copy_with_optics_map()
        integer, parameter :: STKINDS(4) = [1, 1, 2, 3]
        type(sp_project) :: spproj, src, copied
        type(string)     :: dir, empty_dir, src_file, dst_file
        integer          :: imic, iptcl
        write(*,'(A)') 'test_copy_with_optics_map'
        dir       = new_test_dir('copy')
        empty_dir = new_test_dir('copy_empty')
        call make_grouped_mics(spproj)
        call publish_optics_map(spproj, dir, 1, NKEEP)
        ! three micrographs with a stack each and four particles, with no groups
        call src%os_mic%new(3, is_ptcl=.false.)
        call src%os_stk%new(3, is_ptcl=.false.)
        do imic = 1,3
            call src%os_mic%set(imic, 'importind', real(imic))
            call src%os_stk%set(imic, 'fromp', 1)
            call src%os_stk%set(imic, 'top',   1)
        enddo
        call src%os_ptcl2D%new(4, is_ptcl=.true.)
        do iptcl = 1,4
            call src%os_ptcl2D%set_stkind(iptcl, STKINDS(iptcl))
        enddo
        src%os_ptcl3D = src%os_ptcl2D
        src_file = dir//'/source'//METADATA_EXT
        call src%update_projinfo(src_file)
        call src%write(src_file)
        dst_file = dir//'/copied'//METADATA_EXT
        call copy_project_with_optics_map(src_file, dst_file, dir)
        call assert_false(file_exists(dir//'/copied.tmp'), 'the temporary copy is renamed into place')
        call copied%read(dst_file)
        call assert_int(2, copied%os_optics%get_noris(),         'the copy holds the map''s optics groups')
        call assert_int(1, copied%os_mic%get_int(1, 'ogid'),     'micrograph 1 takes group 1')
        call assert_int(2, copied%os_mic%get_int(3, 'ogid'),     'micrograph 3 takes group 2')
        call assert_int(2, copied%os_stk%get_int(2, 'ogid'),     'stack 2 takes its micrograph''s group')
        call assert_int(1, copied%os_ptcl2D%get_int(1, 'ogid'),  'a particle of stack 1 takes group 1')
        call assert_int(2, copied%os_ptcl2D%get_int(4, 'ogid'),  'a particle of stack 3 takes group 2')
        call assert_int(2, copied%os_ptcl3D%get_int(4, 'ogid'),  'in the 3D segment as well')
        call copied%kill
        dst_file = dir//'/copied_without_map'//METADATA_EXT
        call copy_project_with_optics_map(src_file, dst_file, empty_dir)
        call assert_false(file_exists(dir//'/copied_without_map.tmp'), 'so is the plain copy')
        call copied%read(dst_file)
        call assert_int(3, copied%os_mic%get_noris(),            'without a map the project is copied')
        call assert_int(0, copied%os_optics%get_noris(),         'unchanged, without optics groups')
        call copied%kill
        call src%kill
        call spproj%kill
        call simple_rmdir(dir)
        call simple_rmdir(empty_dir)
    end subroutine test_copy_with_optics_map

    ! three micrographs with import indices 1-3 in groups 1, 2, 2, and the two optics rows
    subroutine make_grouped_mics( spproj )
        type(sp_project), intent(inout) :: spproj
        integer, parameter :: OGIDS(3) = [1, 2, 2]
        integer :: imic
        call spproj%os_mic%new(3, is_ptcl=.false.)
        do imic = 1,3
            call spproj%os_mic%set(imic, 'importind', real(imic))
            call spproj%os_mic%set(imic, 'ogid',      real(OGIDS(imic)))
        enddo
        call spproj%os_optics%new(2, is_ptcl=.false.)
        call spproj%os_optics%set(1, 'ogid', 1.0)
        call spproj%os_optics%set(2, 'ogid', 2.0)
    end subroutine make_grouped_mics

    ! an empty directory with an absolute path, as the stream stages pass it
    function new_test_dir( tag ) result( dir )
        character(len=*), intent(in) :: tag
        type(string) :: dir, cwd
        call simple_getcwd(cwd)
        dir = cwd//'/optics_maps_test_'//tag//'_'//int2str(get_process_id())
        if( dir_exists(dir) ) call simple_rmdir(dir)
        call simple_mkdir(dir)
    end function new_test_dir

    logical function map_exists( dir, id )
        type(string), intent(in) :: dir
        integer,      intent(in) :: id
        type(string) :: prefix
        prefix     = dir//'/'//OPTICS_MAP_PREFIX//int2str(id)
        map_exists = file_exists(prefix//TXT_EXT) .and. file_exists(prefix//METADATA_EXT)
    end function map_exists

end module simple_optics_maps_tester
