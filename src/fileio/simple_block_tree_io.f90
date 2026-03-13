!@descr: for storing block tree information on disk
module simple_block_tree_io
use simple_string_utils
use simple_error,   only: simple_exception
use simple_fileio,  only: imat2file, file2imat
use simple_multi_dendro, only: multi_dendro, serialize_multi_dendro, deserialize_multi_dendro
implicit none
public :: write_block_tree, read_block_tree, test_io
private
#include "simple_local_flags.inc"

contains
! -------------------------------------------------------------------
! Write/read multi_dendro to/from single integer file
! -------------------------------------------------------------------
    subroutine write_block_tree(block_tree, fname)
        class(multi_dendro), intent(in)  :: block_tree
        class(string),       intent(in)  :: fname
        integer, allocatable :: bigmat(:,:), trees_meta(:,:), offsets(:), lengths(:), map_int(:,:), tree_pops(:)
        integer, allocatable :: outmat(:,:)
        integer :: n_trees, n_refs, total_nodes, mat_cols, meta_cols, map_cols
        integer :: ncols, nrows, row, pos, i, r, c
        ! Serialize the object into component integer matrices
        call serialize_multi_dendro(block_tree, bigmat, trees_meta, offsets, lengths, map_int, tree_pops)
        ! Gather sizes (handle empty)
        n_trees = size(trees_meta,1)
        if (n_trees < 0) n_trees = 0
        n_refs  = 0
        if (allocated(map_int)) n_refs = size(map_int,2)
        total_nodes = 0
        if (allocated(bigmat)) total_nodes = size(bigmat,1)
        mat_cols = 0
        if (allocated(bigmat)) mat_cols = size(bigmat,2)
        meta_cols = 0
        if (allocated(trees_meta)) meta_cols = size(trees_meta,2)
        map_cols = n_refs
        if (mat_cols == 0) mat_cols = 6   ! fallback/default
        if (meta_cols == 0) meta_cols = 3 ! fallback/default
        ! choose number of columns for the packed matrix
        ncols = max(mat_cols, max(meta_cols, max(map_cols, 1)))
        ! compute total rows needed:
        ! 1 header row + rows for each block (padded to ncols)
        nrows = 1                                    & ! header
                + total_nodes                        & ! bigmat rows
                + n_trees                           & ! trees_meta rows
                + n_trees                           & ! offsets (col)
                + n_trees                           & ! lengths (col)
                + n_trees * n_refs                  & ! map_int rows
                + n_trees                            ! tree_pops (col)

        if (nrows <= 0) then
            allocate(outmat(0, ncols))
        else
            allocate(outmat(nrows, ncols))
            outmat = 0
        end if
        ! Fill header row (row 1)
        if (nrows > 0) then
            outmat(1,1) = n_trees
            outmat(1,2) = n_refs
            outmat(1,3) = total_nodes
            outmat(1,4) = mat_cols
            outmat(1,5) = meta_cols
            outmat(1,6) = map_cols
        end if
        ! Fill blocks
        pos = 2  ! next row to write
        ! 1) bigmat
        if (total_nodes > 0) then
            if (size(bigmat,2) > ncols) THROW_HARD("write_block_tree_to_file: unexpected bigmat column size")
            do r = 1, total_nodes
                outmat(pos + r - 1, 1:size(bigmat,2)) = bigmat(r, :)
            end do
            pos = pos + total_nodes
        end if
        ! 2) trees_meta (n_trees x meta_cols)
        if (n_trees > 0) then
            if (size(trees_meta,2) > ncols) THROW_HARD("write_block_tree_to_file: unexpected trees_meta column size")
            do r = 1, n_trees
                outmat(pos + r - 1, 1:size(trees_meta,2)) = trees_meta(r, :)
            end do
            pos = pos + n_trees
        end if
        ! 3) offsets (n_trees x 1)
        if (n_trees > 0) then
            do r = 1, n_trees
                outmat(pos + r - 1, 1) = offsets(r)
            end do
            pos = pos + n_trees
        end if
        ! 4) lengths (n_trees x 1)
        if (n_trees > 0) then
            do r = 1, n_trees
                outmat(pos + r - 1, 1) = lengths(r)
            end do
            pos = pos + n_trees
        end if
        ! 5) map_int (n_trees x n_refs)
        if (n_trees > 0 .and. n_refs > 0) then
            if (size(map_int,2) > ncols) THROW_HARD("write_block_tree_to_file: unexpected map_int column size")
            do r = 1, n_trees
                outmat(pos + (r-1), 1:n_refs) = map_int(r, 1:n_refs)
            end do
            pos = pos + n_trees * n_refs
        end if
        ! 6) tree_pops (n_trees x 1)
        if (n_trees > 0) then
            do r = 1, n_trees
                outmat(pos + r - 1, 1) = tree_pops(r)
            end do
            pos = pos + n_trees
        end if
        ! Write single matrix to file
        call imat2file(outmat, fname)
        ! cleanup
        if (allocated(outmat)) deallocate(outmat)
        if (allocated(bigmat))  deallocate(bigmat)
        if (allocated(trees_meta)) deallocate(trees_meta)
        if (allocated(offsets)) deallocate(offsets)
        if (allocated(lengths)) deallocate(lengths)
        if (allocated(map_int)) deallocate(map_int)
        if (allocated(tree_pops)) deallocate(tree_pops)
    end subroutine write_block_tree

    subroutine read_block_tree(block_tree, fname)
        class(multi_dendro), intent(inout) :: block_tree
        class(string),       intent(in)      :: fname
        integer, allocatable :: inmat(:,:)
        integer, allocatable :: bigmat(:,:), trees_meta(:,:), offsets(:), lengths(:), map_int(:,:), tree_pops(:)
        integer :: nrows, ncols
        integer :: n_trees, n_refs, total_nodes, mat_cols, meta_cols, map_cols
        integer :: pos, r, c, required_rows
        ! Read full integer matrix from file
        call file2imat(fname, inmat)
        nrows = size(inmat,1)
        ncols = size(inmat,2)
        if (nrows < 1 .or. ncols < 1) then
            THROW_HARD("read_block_tree_from_file: empty or invalid file")
        end if
        ! header is first row
        n_trees   = inmat(1,1)
        n_refs    = inmat(1,2)
        total_nodes = inmat(1,3)
        mat_cols  = inmat(1,4)
        meta_cols = inmat(1,5)
        map_cols  = inmat(1,6)
        if (n_trees < 0 .or. n_refs < 0 .or. total_nodes < 0) then
            THROW_HARD("read_block_tree_from_file: header contains negative sizes")
        end if
        ! compute expected rows to basic-check (non-exhaustive)
        required_rows = 1 + total_nodes + n_trees + n_trees + n_trees * n_refs + n_trees
        if (nrows < required_rows - n_trees * n_refs) then
            ! allow the case n_refs==0 which reduces required rows
        end if
        pos = 2
        ! 1) bigmat
        if (total_nodes > 0) then
            if (mat_cols <= 0) THROW_HARD("read_block_tree_from_file: invalid mat_cols in header")
            allocate(bigmat(total_nodes, mat_cols))
            bigmat(:, :) = 0
            do r = 1, total_nodes
                bigmat(r, 1:mat_cols) = inmat(pos + r - 1, 1:mat_cols)
            end do
            pos = pos + total_nodes
        else
            allocate(bigmat(0, mat_cols))
        end if
        ! 2) trees_meta (n_trees x meta_cols)
        if (n_trees > 0) then
            if (meta_cols <= 0) THROW_HARD("read_block_tree_from_file: invalid meta_cols in header")
            allocate(trees_meta(n_trees, meta_cols))
            do r = 1, n_trees
                trees_meta(r, 1:meta_cols) = inmat(pos + r - 1, 1:meta_cols)
            end do
            pos = pos + n_trees
        else
            allocate(trees_meta(0,0))
        end if
        ! 3) offsets
        if (n_trees > 0) then
            allocate(offsets(n_trees))
            do r = 1, n_trees
                offsets(r) = inmat(pos + r - 1, 1)
            end do
            pos = pos + n_trees
        else
            allocate(offsets(0))
        end if
        ! 4) lengths
        if (n_trees > 0) then
            allocate(lengths(n_trees))
            do r = 1, n_trees
                lengths(r) = inmat(pos + r - 1, 1)
            end do
            pos = pos + n_trees
        else
            allocate(lengths(0))
        end if
        ! 5) map_int (n_trees x n_refs)
        if (n_trees > 0 .and. n_refs > 0) then
            allocate(map_int(n_trees, n_refs))
            do r = 1, n_trees
                map_int(r, 1:n_refs) = inmat(pos + (r - 1), 1:n_refs)
            end do
            pos = pos + n_trees * n_refs
        else
            allocate(map_int(0,0))
        end if
        ! 6) tree_pops
        if (n_trees > 0) then
            allocate(tree_pops(n_trees))
            do r = 1, n_trees
                tree_pops(r) = inmat(pos + r - 1, 1)
            end do
            pos = pos + n_trees
        else
            allocate(tree_pops(0))
        end if
        ! Finally reconstruct multi_dendro from these integer components
        call deserialize_multi_dendro(block_tree, bigmat, trees_meta, offsets, lengths, map_int, tree_pops)
        ! cleanup
        if (allocated(inmat))      deallocate(inmat)
        if (allocated(bigmat))     deallocate(bigmat)
        if (allocated(trees_meta)) deallocate(trees_meta)
        if (allocated(offsets))    deallocate(offsets)
        if (allocated(lengths))    deallocate(lengths)
        if (allocated(map_int))    deallocate(map_int)
        if (allocated(tree_pops))  deallocate(tree_pops)
    end subroutine read_block_tree

    subroutine test_io 
        type(multi_dendro) :: md_orig, md_read
        integer, allocatable           :: labels(:), refs(:)
        real(kind=real32), allocatable :: subdist(:,:)
        integer :: nrefs, i, itree, j
        type(string) :: fname
        logical :: ok
        ! local filename
        fname = 'multi_dendro_test.bin'
        ! 1) Build a simple multi_dendro:
        !    - n_refs = 4, two trees: labels = [1,1,2,2]
        nrefs = 4
        allocate(labels(nrefs))
        labels = [1, 1, 2, 2]
        call md_orig%kill()
        call md_orig%new(labels)
        ! For each tree, build a simple 2x2 sub-distance matrix:
        do itree = 1, md_orig%get_n_trees()
            refs = md_orig%get_tree_refs(itree)  ! returns GLOBAL ref ids for this tree
            if (size(refs) <= 0) cycle
            ! Construct simple distance: diagonal 0, off-diag 1.0
            allocate(subdist(size(refs), size(refs)))
            subdist = 0.0
            if (size(refs) > 1) then
                do i = 1, size(refs)
                    do j = 1, size(refs)
                        if (i /= j) subdist(i,j) = 1.0
                    end do
                end do
            end if
            ! linkage code: supply 1 (single) — exact value not important for this tiny test
            call md_orig%build_tree_from_subdistmat(itree, refs, subdist, 1)
            if (allocated(subdist)) deallocate(subdist)
            if (allocated(refs)) deallocate(refs)
        end do
        ! 2) Write to file
        call write_block_tree(md_orig, fname)  ! uses module routine you previously added
        ! 3) Read from file into md_read
        call md_read%kill()
        call read_block_tree(md_read, fname)
        ! 4) Compare original and read back via serialization
        ok = compare_multi_dendro(md_orig, md_read)
        if (ok) then
            print *, 'TEST PASS: reconstructed multi_dendro matches original'
        else
            print *, 'TEST FAIL: reconstructed multi_dendro differs from original'
            stop 1
        end if
        ! cleanup
        call md_orig%kill()
        call md_read%kill()
        contains
            !----------------------------------------------------------------
            ! Compare two multi_dendro objects by serializing both into integer
            ! blocks and checking exact equality of all integer blocks.
            !----------------------------------------------------------------
            logical function compare_multi_dendro(a, b)
                class(multi_dendro), intent(in) :: a, b
                integer, allocatable :: mat_a(:,:), tmeta_a(:,:), offs_a(:), lens_a(:), map_a(:,:), pops_a(:)
                integer, allocatable :: mat_b(:,:), tmeta_b(:,:), offs_b(:), lens_b(:), map_b(:,:), pops_b(:)
                logical :: same

                same = .false.

                ! serialize a
                call serialize_multi_dendro(a, mat_a, tmeta_a, offs_a, lens_a, map_a, pops_a)
                ! serialize b
                call serialize_multi_dendro(b, mat_b, tmeta_b, offs_b, lens_b, map_b, pops_b)

                ! Compare sizes first
                if (size(mat_a,1) /= size(mat_b,1) .or. size(mat_a,2) /= size(mat_b,2)) return
                if (size(tmeta_a,1) /= size(tmeta_b,1) .or. size(tmeta_a,2) /= size(tmeta_b,2)) return
                if (size(offs_a) /= size(offs_b)) return
                if (size(lens_a) /= size(lens_b)) return
                if (size(map_a,1) /= size(map_b,1) .or. size(map_a,2) /= size(map_b,2)) return
                if (size(pops_a) /= size(pops_b)) return

                ! Elementwise comparisons
                if (any(mat_a /= mat_b)) return
                if (any(tmeta_a /= tmeta_b)) return
                if (any(offs_a /= offs_b)) return
                if (any(lens_a /= lens_b)) return
                if (any(map_a /= map_b)) return
                if (any(pops_a /= pops_b)) return

                same = .true.

                ! cleanup
                if (allocated(mat_a)) deallocate(mat_a)
                if (allocated(tmeta_a)) deallocate(tmeta_a)
                if (allocated(offs_a)) deallocate(offs_a)
                if (allocated(lens_a)) deallocate(lens_a)
                if (allocated(map_a)) deallocate(map_a)
                if (allocated(pops_a)) deallocate(pops_a)

                if (allocated(mat_b)) deallocate(mat_b)
                if (allocated(tmeta_b)) deallocate(tmeta_b)
                if (allocated(offs_b)) deallocate(offs_b)
                if (allocated(lens_b)) deallocate(lens_b)
                if (allocated(map_b)) deallocate(map_b)
                if (allocated(pops_b)) deallocate(pops_b)

                compare_multi_dendro = same
            end function compare_multi_dendro

    end subroutine test_io

end module simple_block_tree_io 