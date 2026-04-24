!@descr: for storing search-space maps on disk
module simple_srchspace_map2D_io
use simple_error,         only: simple_exception
use simple_fileio,        only: imat2file, file2imat
use simple_srch_sort_loc, only: hpsort
use simple_string,        only: string
implicit none

public :: write_srchspace_map2D, read_srchspace_map2D, test_srchspace_map2D_io
private
#include "simple_local_flags.inc"

contains

    ! Build a search-space map from 2D class->cluster labels.
    ! Input labels can be sparse (e.g. parent class ids); they are remapped to
    ! dense [1..nsub] so the stored map is a valid srchspace_map:
    !   nspace = #classes (full space)
    !   nsub   = #clusters (coarse space)
    subroutine write_srchspace_map2D(cluster_of_class, fname)
        integer,       intent(in) :: cluster_of_class(:)
        class(string), intent(in) :: fname
        integer, allocatable :: unique_clusters(:), full2sub_map(:), sub2full_map(:), outmat(:,:)
        integer :: i, j, nspace, nsub, nfound, ncols
        logical :: found
        nspace = size(cluster_of_class)
        if( nspace < 1 ) THROW_HARD('write_srchspace_map2D: empty class->cluster labels')
        if( minval(cluster_of_class) < 1 ) THROW_HARD('write_srchspace_map2D: cluster labels must be >= 1')
        allocate(unique_clusters(nspace), source=0)
        nsub = 0
        do i = 1, nspace
            found = .false.
            do j = 1, nsub
                if( unique_clusters(j) == cluster_of_class(i) )then
                    found = .true.
                    exit
                endif
            end do
            if( .not. found )then
                nsub = nsub + 1
                unique_clusters(nsub) = cluster_of_class(i)
            endif
        end do
        call hpsort(unique_clusters(1:nsub))
        allocate(full2sub_map(nspace), source=0)
        allocate(sub2full_map(nsub), source=0)
        do i = 1, nspace
            found = .false.
            do j = 1, nsub
                if( cluster_of_class(i) == unique_clusters(j) )then
                    full2sub_map(i) = j
                    found = .true.
                    exit
                endif
            end do
            if( .not. found )then
                THROW_HARD('write_srchspace_map2D: failed to remap cluster label to dense index')
            endif
        end do
        do j = 1, nsub
            nfound = 0
            do i = 1, nspace
                if( full2sub_map(i) == j )then
                    sub2full_map(j) = i
                    nfound = nfound + 1
                    exit
                endif
            end do
            if( nfound == 0 ) THROW_HARD('write_srchspace_map2D: dense cluster index without members')
        end do
        ! 2D map format extends generic srchspace_map layout with row 4:
        ! sparse cluster labels for each dense sub-index.
        ncols = max(max(nspace, nsub), 2)
        allocate(outmat(4, ncols), source=0)
        outmat(1,1) = nspace
        outmat(1,2) = nsub
        outmat(2,1:nsub) = sub2full_map
        outmat(3,1:nspace) = full2sub_map
        outmat(4,1:nsub) = unique_clusters(1:nsub)
        call imat2file(outmat, fname)
        if( allocated(unique_clusters) ) deallocate(unique_clusters)
        if( allocated(outmat) )          deallocate(outmat)
        if( allocated(full2sub_map) )    deallocate(full2sub_map)
        if( allocated(sub2full_map) )    deallocate(sub2full_map)
    end subroutine write_srchspace_map2D

    subroutine read_srchspace_map2D(cluster_of_class, fname)
        integer, allocatable, intent(out) :: cluster_of_class(:)
        class(string),        intent(in)  :: fname
        integer, allocatable :: inmat(:,:), full2sub_map(:), dense2cluster(:)
        integer :: nrows, ncols, nspace, nsub, i
        call file2imat(fname, inmat)
        nrows = size(inmat,1)
        ncols = size(inmat,2)
        if (nrows < 3 .or. ncols < 2) THROW_HARD('read_srchspace_map2D: empty or invalid file')
        nspace = inmat(1,1)
        nsub   = inmat(1,2)
        if (nspace < 1 .or. nsub < 1) THROW_HARD('read_srchspace_map2D: header contains invalid sizes')
        if (ncols < max(nspace, nsub)) THROW_HARD('read_srchspace_map2D: file payload too small for header sizes')
        allocate(full2sub_map(nspace), source=inmat(3,1:nspace))
        if (minval(full2sub_map) < 1 .or. maxval(full2sub_map) > nsub) then
            THROW_HARD('read_srchspace_map2D: full2sub_map values outside [1,nsub]')
        endif
        if( nrows >= 4 )then
            allocate(dense2cluster(nsub), source=inmat(4,1:nsub))
            if( minval(dense2cluster) < 1 ) THROW_HARD('read_srchspace_map2D: invalid sparse cluster labels')
        else
            allocate(dense2cluster(nsub))
            dense2cluster = (/(i, i=1,nsub)/)
        endif
        allocate(cluster_of_class(nspace), source=0)
        do i = 1, nspace
            cluster_of_class(i) = dense2cluster(full2sub_map(i))
        end do
        if (allocated(inmat))        deallocate(inmat)
        if (allocated(full2sub_map)) deallocate(full2sub_map)
        if (allocated(dense2cluster)) deallocate(dense2cluster)
    end subroutine read_srchspace_map2D

    ! Sandbox test modelled on the apply_split_project_updates use case:
    ! simulate nptcls particles assigned to nsplit subclasses across ncls parent
    ! classes, write the resulting cluster_of_class map, read it back and verify
    ! the full roundtrip recovers the original sparse parent-class labels.
    subroutine test_srchspace_map2D_io
        type(string) :: fname
        integer, allocatable :: cluster_of_class_in(:), cluster_of_class_out(:)
        ! --- case 1: basic sparse labels (same as original test) ---
        fname = 'srchspace_map_test.bin'
        call write_srchspace_map2D([10, 10, 40, 40, 80, 80], fname)
        call read_srchspace_map2D(cluster_of_class_out, fname)
        if( any(cluster_of_class_out /= [10, 10, 40, 40, 80, 80]) ) &
            THROW_HARD('simple_srchspace_map2D_io::test_srchspace_map2D_io case 1: cluster map roundtrip mismatch')
        deallocate(cluster_of_class_out)
        call fname%kill
        ! --- case 2: simulate apply_split_project_updates output ---
        ! 3 parent classes (ids 5, 12, 99), each split into 2 subclasses => 6 subclasses total.
        ! cluster_of_class(iglob) = parent class id, matching parent_of_subcls(:) in the strategy.
        allocate(cluster_of_class_in(6))
        cluster_of_class_in = [5, 5, 12, 12, 99, 99]
        fname = 'srchspace_map_test2.bin'
        call write_srchspace_map2D(cluster_of_class_in, fname)
        call read_srchspace_map2D(cluster_of_class_out, fname)
        if( size(cluster_of_class_out) /= 6 ) &
            THROW_HARD('simple_srchspace_map2D_io::test_srchspace_map2D_io case 2: wrong nspace on readback')
        if( any(cluster_of_class_out /= cluster_of_class_in) ) &
            THROW_HARD('simple_srchspace_map2D_io::test_srchspace_map2D_io case 2: sparse parent labels not recovered')
        deallocate(cluster_of_class_in)
        deallocate(cluster_of_class_out)
        call fname%kill
        ! --- case 3: single parent class (edge case: nsub=1) ---
        fname = 'srchspace_map_test3.bin'
        call write_srchspace_map2D([7, 7, 7], fname)
        call read_srchspace_map2D(cluster_of_class_out, fname)
        if( any(cluster_of_class_out /= [7, 7, 7]) ) &
            THROW_HARD('simple_srchspace_map2D_io::test_srchspace_map2D_io case 3: single-cluster map roundtrip mismatch')
        deallocate(cluster_of_class_out)
        call fname%kill
        write(*,'(a)') 'SIMPLE_SRCHSPACE_MAP2D_IO: TEST_SRCHSPACE_MAP2D_IO COMPLETED'
    end subroutine test_srchspace_map2D_io

end module simple_srchspace_map2D_io
