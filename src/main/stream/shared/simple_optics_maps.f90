!@descr: the versioned optics-map files the optics stage publishes and the later stages import
!==============================================================================
! MODULE: simple_optics_maps
!
! PURPOSE:
!   Owns the file protocol between the optics stage (writer) and the stages
!   that apply its optics groups (readers):
!     <dir>/optics_map_<id>.txt     importind -> ogid table
!     <dir>/optics_map_<id>.simple  the optics segment
!   written by sp_project%write_optics_map and read by
!   sp_project%import_optics_map. Ids grow by one per publication; the writer
!   keeps the newest nkeep maps, and readers take the highest id present.
!
!   The writer used to prune by testing every id from 1 up on each
!   publication, and readers found the newest map with three copies of the
!   same directory scan (get_latest_optics_map_id in simple_stream_utils and
!   two internal get_latest_optics_map in simple_stream_refine2D_utils).
!   Readers switch to import_latest_optics_map when their stages move over.
!   copy_project_with_optics_map hands a project to the next stage with the
!   newest groups applied (the particle sieve uses it for finished chunks).
!   latest_optics_map_table reads the newest map's table as it is (the
!   optics stage restores its groups from it on a restart).
!
! HOME:
!   In src/main/stream/shared for now; it belongs beside sp_project's
!   write_optics_map / import_optics_map.
!
! TESTS:
!   simple_optics_maps_tester
!==============================================================================
module simple_optics_maps
use simple_defs_fname,   only: OPTICS_MAP_PREFIX, TXT_EXT, METADATA_EXT
use simple_string,       only: string
use simple_string_utils, only: int2str
use simple_fileio,       only: del_file, file_exists, simple_copy_file, simple_rename, swap_suffix
use simple_nrtxtfile,    only: nrtxtfile
use simple_sp_project,   only: sp_project
use simple_stream_utils, only: get_latest_optics_map_id
implicit none

public :: publish_optics_map, import_latest_optics_map, latest_optics_map_id, copy_project_with_optics_map,&
          &latest_optics_map_table
private

contains

    !> Writes the optics groups of @p spproj as map @p id in @p dir, and removes map
    !! @p id - @p nkeep, so the newest @p nkeep maps remain.
    subroutine publish_optics_map( spproj, dir, id, nkeep )
        class(sp_project), intent(inout) :: spproj
        class(string),     intent(in)    :: dir
        integer,           intent(in)    :: id, nkeep
        type(string) :: prefix
        prefix = map_prefix(dir, id)
        call spproj%write_optics_map(prefix%to_char())
        if( id - nkeep < 1 ) return
        prefix = map_prefix(dir, id - nkeep)
        if( file_exists(prefix//TXT_EXT)      ) call del_file(prefix//TXT_EXT)
        if( file_exists(prefix//METADATA_EXT) ) call del_file(prefix//METADATA_EXT)
    end subroutine publish_optics_map

    !> Applies the newest map in @p dir to @p spproj; returns its id, 0 when there is none.
    function import_latest_optics_map( spproj, dir ) result( id )
        class(sp_project), intent(inout) :: spproj
        class(string),     intent(in)    :: dir
        integer :: id
        id = latest_optics_map_id(dir)
        if( id > 0 ) call spproj%import_optics_map(map_prefix(dir, id))
    end function import_latest_optics_map

    !> Writes the project @p src as @p dst with the groups of the newest map in @p dir applied to
    !! its micrographs, stacks and particles, and the map's optics segment; an exact copy when
    !! there is no map yet. The groups are applied by import index, so a project merged from
    !! several (whose group ids were offset per source) gets the map's ids back. The copy is
    !! written as <dst stem>.tmp and renamed, so a stage watching for @p dst sees it whole.
    subroutine copy_project_with_optics_map( src, dst, dir )
        class(string), intent(in) :: src, dst, dir
        type(sp_project) :: spproj
        type(string)     :: tmp
        integer          :: id
        if( latest_optics_map_id(dir) == 0 )then
            tmp = swap_suffix(dst, '.tmp', METADATA_EXT)
            call simple_copy_file(src, tmp)
            call simple_rename(tmp, dst)
            return
        endif
        call spproj%read(src)
        id = import_latest_optics_map(spproj, dir)
        call spproj%write(dst, tempfile=.true.)
        call spproj%kill
    end subroutine copy_project_with_optics_map

    !> The table of the newest map in @p dir: the import indices it lists and their group ids
    !! (empty when there is no map, or it cannot be read). Returns the map's id, 0 when none.
    function latest_optics_map_table( dir, importinds, ogids ) result( id )
        class(string),        intent(in)  :: dir
        integer, allocatable, intent(out) :: importinds(:), ogids(:)
        integer         :: id
        type(nrtxtfile) :: mapfile
        type(string)    :: fname
        real            :: entry(2)
        integer         :: il, nl
        allocate(importinds(0), ogids(0))
        id = latest_optics_map_id(dir)
        if( id == 0 ) return
        fname = map_prefix(dir, id)//TXT_EXT
        if( .not. file_exists(fname) ) return
        call mapfile%new(fname, 1)
        if( mapfile%get_nrecs_per_line() /= 2 )then
            call mapfile%kill()
            return
        endif
        nl = mapfile%get_ndatalines()
        deallocate(importinds, ogids)
        allocate(importinds(nl), ogids(nl))
        do il = 1,nl
            call mapfile%readNextDataLine(entry)
            importinds(il) = nint(entry(1))
            ogids(il)      = nint(entry(2))
        enddo
        call mapfile%kill()
    end function latest_optics_map_table

    !> The highest map id in @p dir, 0 when there is none.
    function latest_optics_map_id( dir ) result( id )
        class(string), intent(in) :: dir
        integer :: id
        id = get_latest_optics_map_id(dir)
    end function latest_optics_map_id

    function map_prefix( dir, id ) result( prefix )
        class(string), intent(in) :: dir
        integer,       intent(in) :: id
        type(string) :: prefix
        prefix = dir//'/'//OPTICS_MAP_PREFIX//int2str(id)
    end function map_prefix

end module simple_optics_maps
