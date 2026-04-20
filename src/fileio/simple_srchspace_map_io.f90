!@descr: for storing search-space maps on disk
module simple_srchspace_map_io
use simple_error,         only: simple_exception
use simple_fileio,        only: imat2file, file2imat
use simple_srchspace_map, only: srchspace_map
use simple_string,        only: string
implicit none

public :: write_srchspace_map, read_srchspace_map, test_io
private
#include "simple_local_flags.inc"

contains

    subroutine write_srchspace_map(map, fname)
        class(srchspace_map), intent(in) :: map
        class(string),        intent(in) :: fname
        integer, allocatable :: outmat(:,:), full2sub_map(:), sub2full_map(:)
        integer :: nspace, nsub, ncols
        if (.not. map%has_data()) THROW_HARD('write_srchspace_map: empty srchspace_map')
        nspace = map%get_nspace()
        nsub   = map%get_nsub()
        ncols = max(max(nspace, nsub), 2)
        full2sub_map = map%get_full2sub_map()
        sub2full_map = map%get_sub2full_map()
        allocate(outmat(3, ncols), source=0)
        outmat(1,1) = nspace
        outmat(1,2) = nsub
        outmat(2,1:nsub) = sub2full_map
        outmat(3,1:nspace) = full2sub_map
        call imat2file(outmat, fname)
        if (allocated(outmat))       deallocate(outmat)
        if (allocated(full2sub_map)) deallocate(full2sub_map)
        if (allocated(sub2full_map)) deallocate(sub2full_map)
    end subroutine write_srchspace_map

    subroutine read_srchspace_map(map, fname)
        class(srchspace_map), intent(inout) :: map
        class(string),        intent(in)    :: fname
        integer, allocatable :: inmat(:,:), full2sub_map(:), sub2full_map(:)
        integer :: nrows, ncols, nspace, nsub
        call map%kill()
        call file2imat(fname, inmat)
        nrows = size(inmat,1)
        ncols = size(inmat,2)
        if (nrows < 3 .or. ncols < 2) THROW_HARD('read_srchspace_map: empty or invalid file')
        nspace = inmat(1,1)
        nsub   = inmat(1,2)
        if (nspace < 1 .or. nsub < 1) THROW_HARD('read_srchspace_map: header contains invalid sizes')
        if (ncols < max(nspace, nsub)) THROW_HARD('read_srchspace_map: file payload too small for header sizes')
        allocate(sub2full_map(nsub), source=inmat(2,1:nsub))
        allocate(full2sub_map(nspace), source=inmat(3,1:nspace))
        if (minval(sub2full_map) < 1 .or. maxval(sub2full_map) > nspace) then
            THROW_HARD('read_srchspace_map: sub2full_map values outside [1,nspace]')
        endif
        if (minval(full2sub_map) < 1 .or. maxval(full2sub_map) > nsub) then
            THROW_HARD('read_srchspace_map: full2sub_map values outside [1,nsub]')
        endif
        call map%new(full2sub_map, sub2full_map)
        if (allocated(inmat))        deallocate(inmat)
        if (allocated(full2sub_map)) deallocate(full2sub_map)
        if (allocated(sub2full_map)) deallocate(sub2full_map)
    end subroutine read_srchspace_map

    subroutine test_io
        type(srchspace_map) :: map_orig, map_read
        type(string) :: fname
        integer, allocatable :: full2sub_map(:), sub2full_map(:)
        fname = 'srchspace_map_test.bin'
        call map_orig%new([1, 1, 2, 2, 3, 3], [1, 3, 5])
        call write_srchspace_map(map_orig, fname)
        call read_srchspace_map(map_read, fname)
        full2sub_map = map_read%get_full2sub_map()
        sub2full_map = map_read%get_sub2full_map()
        if (any(full2sub_map /= [1, 1, 2, 2, 3, 3])) THROW_HARD('simple_srchspace_map_io::test_io full2sub mismatch')
        if (any(sub2full_map /= [1, 3, 5])) THROW_HARD('simple_srchspace_map_io::test_io sub2full mismatch')
        call map_orig%kill()
        call map_read%kill()
        if (allocated(full2sub_map)) deallocate(full2sub_map)
        if (allocated(sub2full_map)) deallocate(sub2full_map)
    end subroutine test_io

end module simple_srchspace_map_io
