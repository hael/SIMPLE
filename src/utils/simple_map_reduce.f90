! routines for distributed SIMPLE execution
module simple_map_reduce
use simple_defs
use simple_strings, only: int2str, int2str_pad
use simple_error,   only: allocchk
use simple_fileio,  only: fopen, fileiochk, fclose, file2rarr
use simple_jiffys,  only: progress
use simple_math,    only: hpsort
implicit none
private
#include "simple_local_flags.inc"
public ::  split_nobjs_even,split_pairs_in_parts,merge_similarities_from_parts,&
&merge_rmat_from_parts,merge_nnmat_from_parts, approx_balanced_partitioning
contains

    !>  \brief  for generating balanced partitions of nobjs objects
    function split_nobjs_even( nobjs, nparts, szmax ) result( parts )
        integer,           intent(in)  :: nobjs, nparts
        integer, optional, intent(out) :: szmax
        integer, allocatable :: parts(:,:)
        integer :: nobjs_per_part, leftover, istop, istart, ipart, sszmax
        allocate(parts(nparts,2), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('In: simple_map_reduce :: split_nobjs_even',alloc_stat)
        nobjs_per_part = nobjs/nparts
        leftover = nobjs-nobjs_per_part*nparts
        istop  = 0
        istart = 0
        sszmax = 0
        do ipart=1,nparts
            if( ipart == nparts )then
                istart = istop+1
                istop  = nobjs
            else
                if( leftover == 0 )then
                    istart = istop+1
                    istop  = istart+nobjs_per_part-1
                else
                    istop    = ipart*(nobjs_per_part+1)
                    istart   = istop-(nobjs_per_part+1)+1
                    leftover = leftover-1
                endif
            endif
            parts(ipart,1) = istart
            parts(ipart,2) = istop
            sszmax         = max(sszmax,istop - istart + 1)
        end do
        if( present(szmax) ) szmax = sszmax
    end function split_nobjs_even

    subroutine approx_balanced_partitioning( nrs, n, k, kassgn )
        integer,  intent(in)    :: n, k
        real,     intent(inout) :: nrs(n)
        integer,  intent(inout) :: kassgn(n)
        integer, allocatable :: inds(:)
        real,    allocatable :: sums(:)
        integer :: i, loc(1)
        allocate(inds(n), source=(/(i,i=1,n)/))
        allocate(sums(k), source=0.)
        call hpsort(nrs, inds)
        do i=n,1,-1 ! traverse in descending order
            if( nrs(i) > TINY )then
                loc = minloc(sums)
                ! assignment
                kassgn(inds(i)) = loc(1)
                ! update sum
                sums(loc(1)) = sums(loc(1)) + nrs(i)
            else
                kassgn(inds(i)) = 0
            endif
        end do
    end subroutine approx_balanced_partitioning

    !>  \brief  for generating balanced partitions for pairwise calculations on nobjs ojects
    subroutine split_pairs_in_parts( nobjs, nparts )
        integer, intent(in)  :: nobjs  !< number objects to analyse in pairs
        integer, intent(in)  :: nparts !< number of partitions (nodes) for parallel execution
        integer              :: npairs, cnt, funit, i, j
        integer              :: ipart, io_stat, numlen
        integer, allocatable :: pairs(:,:), parts(:,:)
        character(len=:), allocatable :: fname
        ! generate all pairs
        DebugPrint  ' split_pairs_in_parts nobjs: ', nobjs
        npairs = (nobjs*(nobjs-1))/2
        DebugPrint  'split_pairs_in_parts, npairs: ', npairs
        allocate( pairs(npairs,2), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('mapreduce ;split_pairs_in_parts  1',alloc_stat)
        cnt = 0
        do i=1,nobjs-1
            do j=i+1,nobjs
                cnt = cnt+1
                pairs(cnt,1) = i
                pairs(cnt,2) = j
            end do
        end do
        ! generate balanced partitions of pairs
        parts  = split_nobjs_even(npairs, nparts)     ! realloc lhs
        numlen = len(int2str(nparts))
        ! write the partitions
        do ipart=1,nparts
            call progress(ipart,nparts)
            allocate(fname, source='pairs_part' // int2str_pad(ipart,numlen) // '.bin', stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("mapreduce ;split_pairs_in_parts creating fname ",alloc_stat)
            call fopen(funit, status='REPLACE', action='WRITE', file=fname, access='STREAM',iostat=io_stat)
            call fileiochk('mapreduce ;split_pairs_in_parts '//trim(fname), io_stat)
            DebugPrint   'writing pairs in range: ', parts(ipart,1), parts(ipart,2)
            write(unit=funit,pos=1,iostat=io_stat) pairs(parts(ipart,1):parts(ipart,2),:)
            ! Check if the write was successful
            if( io_stat .ne. 0 )&
                call fileiochk('mapreduce ;split_pairs_in_parts writing to '//trim(fname), io_stat)
            call fclose(funit,errmsg='mapreduce ;split_pairs_in_parts ')
            deallocate(fname)
        end do
        deallocate(pairs, parts)
    end subroutine split_pairs_in_parts

    !>  \brief  for merging partial calculations of similarities
    function merge_similarities_from_parts( nobjs, nparts ) result( smat )
        integer, intent(in) :: nobjs  !< number objects to analyse in pairs
        integer, intent(in) :: nparts !< number of partitions (nodes) for parallel execution
        integer :: npairs, funit, ipart, i, io_stat, numlen
        real,    allocatable :: smat(:,:), sims(:)
        integer, allocatable :: pairs(:,:), parts(:,:)
        character(len=:), allocatable :: fname
        ! allocate pairs and similarities
        npairs = (nobjs*(nobjs-1))/2
        DebugPrint   'analysing this number of objects: ', nobjs
        DebugPrint   'analysing this number of pairs: ', npairs
        allocate( smat(nobjs,nobjs), pairs(npairs,2), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('In: simple_map_reduce::merge_similarities_from_parts, 1',alloc_stat)
        ! initialise similarities
        smat = -1.
        do i=1,nobjs
            smat(i,i) = 1
        end do
        ! repeat the generation of balanced partitions
        parts  = split_nobjs_even(npairs, nparts)   !! realloc lhs
        numlen = len(int2str(nparts))
        ! compress the partial similarities into a single matrix
        do ipart=1,nparts
            ! read the index mapping
            allocate(fname, source='pairs_part'//int2str_pad(ipart,numlen)//'.bin')
            call fopen(funit, status='OLD', action='READ', file=fname, access='STREAM',iostat=io_stat)
            call fileiochk("In map_reduce merge_similarities_from_parts opening "//trim(fname), io_stat)
            read(unit=funit,pos=1,iostat=io_stat) pairs(parts(ipart,1):parts(ipart,2),:)
            ! Check if the read was successful
            if( io_stat .ne. 0 )&
                 call fileiochk('**ERROR(merge_similarities_from_parts): I/O error when reading file: '&
                 //trim(fname),io_stat)
            call fclose(funit, errmsg="In map_reduce merge_similarities_from_parts closing "//trim(fname))
            deallocate(fname)
            ! retrieve the similarities
            DebugPrint   'allocating this number of similarities: ', parts(ipart,2)-parts(ipart,1)+1
            allocate(sims(parts(ipart,1):parts(ipart,2)), stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("In: simple_map_reduce::merge_similarities_from_parts, sims(:) 2",alloc_stat)
            allocate(fname, source='similarities_part'//int2str_pad(ipart,numlen)//'.bin')
            call fopen(funit, status='OLD', action='READ', file=fname, access='STREAM',iostat=io_stat)
            call fileiochk("In map_reduce merge_rmat_from_parts opening "//trim(fname), io_stat)
            read(unit=funit,pos=1,iostat=io_stat) sims(parts(ipart,1):parts(ipart,2))
            ! check if the read was successful
            if( io_stat .ne. 0 )then
               call fileiochk('**ERROR(merge_similarities_from_parts): I/O error when reading file: '&
                    //trim(fname),io_stat)
            endif
            call fclose(funit, errmsg="In map_reduce merge_similarities_from_parts opening "//trim(fname))
            deallocate(fname)
            if(alloc_stat.ne.0)call allocchk('In: simple_map_reduce::merge_similarities_from_parts dealloc',alloc_stat)
            ! set the similarity matrix components
            do i=parts(ipart,1),parts(ipart,2)
                smat(pairs(i,1),pairs(i,2)) = sims(i)
                smat(pairs(i,2),pairs(i,1)) = sims(i)
            end do
            deallocate(sims)
            if(alloc_stat.ne.0)call allocchk('In: simple_map_reduce::merge_similarities_from_parts, 2',alloc_stat)
        end do
        deallocate(parts,pairs)
    end function merge_similarities_from_parts

    !>  \brief  for merging partial calculations of the nearest neighbour matrix
    function merge_nnmat_from_parts( nobjs, nparts, nnn ) result( nnmat )
        integer, intent(in)  :: nobjs  !< number of objects (images)
        integer, intent(in)  :: nparts !< number of partitions, assuming even split (CPU exec)
        integer, intent(in)  :: nnn    !< number of nearest neighbours
        integer, allocatable :: nnmat(:,:), parts(:,:)
        character(len=:), allocatable :: fname
        integer :: numlen, ipart, funit, io_stat
        ! allocate nearest neighbour matrix
        allocate(nnmat(nobjs,nnn), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In: simple_map_reduce :: merge_nnmat_from_parts",alloc_stat)
        ! repeat the generation of balanced partitions
        parts  = split_nobjs_even(nobjs, nparts)                                        !! realloc lhs
        numlen = len(int2str(nparts))
        ! compress the partial nearest neighbour matrices into a single matrix
        do ipart=1,nparts
            allocate(fname, source='nnmat_part'//int2str_pad(ipart,numlen)//'.bin')
            if(alloc_stat.ne.0)call allocchk('In: simple_map_reduce::merge_nnmat_from_parts, 1',alloc_stat)
            call fopen(funit, status='OLD', action='READ', file=fname, access='STREAM',iostat=io_stat)
            call fileiochk("In map_reduce merge_nnmat_from_parts opening "//trim(fname), io_stat)
            read(unit=funit,pos=1,iostat=io_stat) nnmat(parts(ipart,1):parts(ipart,2),:)
            ! check if the read was successful
            if( io_stat .ne. 0 )then
               call fileiochk('**ERROR(merge_nnmat_from_parts): I/O error when reading file: '&
                    //trim(fname),io_stat)
            endif
            call fclose(funit,errmsg="In map_reduce merge_nnmat_from_parts closing "//trim(fname))
            deallocate(fname)
        end do
        deallocate(parts)
    end function merge_nnmat_from_parts

    !>  \brief  for merging partial calculations of the nearest neighbour matrix
    function merge_rmat_from_parts( nobjs, nparts, ydim, fbody ) result( mat_merged )
        integer,          intent(in)  :: nobjs  !< number of objects (images)
        integer,          intent(in)  :: nparts !< number of partitions, assuming even split (CPU exec)
        integer,          intent(in)  :: ydim   !< y dimension
        character(len=*), intent(in)  :: fbody  !< file body of partial files
        integer,          allocatable :: parts(:,:)
        real,             allocatable :: mat_merged(:,:)
        character(len=:), allocatable :: fname
        integer :: numlen, ipart, funit, io_stat
        ! allocate merged matrix
        allocate(mat_merged(nobjs,ydim), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In: simple_map_reduce :: merge_mat_from_parts",alloc_stat)
        ! repeat the generation of balanced partitions
        parts  = split_nobjs_even(nobjs, nparts)
        numlen = len(int2str(nparts))
        ! compress the partial matrices into a single matrix
        do ipart=1,nparts
            allocate(fname, source=trim(fbody)//int2str_pad(ipart,numlen)//'.bin', stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("In: simple_map_reduce :: merge_mat_from_parts ",alloc_stat)
            call fopen(funit, status='OLD', action='READ', file=fname, access='STREAM',iostat=io_stat)
            call fileiochk("In map_reduce merge_rmat_from_parts opening "//trim(fname), io_stat)
            read(unit=funit,pos=1,iostat=io_stat) mat_merged(parts(ipart,1):parts(ipart,2),:)
            ! check if the read was successful
            if( io_stat .ne. 0 )&
                 call fileiochk('**ERROR(merge_rmat_from_parts): I/O error when reading file: '&
                 //trim(fname),io_stat)
            call fclose(funit,errmsg="In map_reduce merge_rmat_from_parts opening "//trim(fname))
            deallocate(fname)
        end do
        deallocate(parts)
    end function merge_rmat_from_parts

    subroutine merge_proj_weights_from_parts( nparts, proj_weights )
        integer,           intent(in)  :: nparts
        real, allocatable, intent(out) :: proj_weights(:)
        character(len=:), allocatable  :: fname
        real, allocatable :: pwpart(:)
        logical :: file_exists
        integer :: ipart, numlen
        if( allocated(proj_weights) ) deallocate(proj_weights)
        numlen = len(int2str(nparts))
        do ipart=1,nparts
            allocate(fname, source=PROJ_WEIGHTS_FBODY//int2str_pad(ipart, numlen)//BIN_EXT)
            if( file_exists )then
                pwpart = file2rarr(fname)
                if( allocated(proj_weights) )then
                    proj_weights = proj_weights + pwpart
                else
                    proj_weights = pwpart
                endif
            else
                write(logfhandle,'(A)') 'WARNING! '//fname//' does not exist; merge_proj_weights_from_parts'
            endif
            deallocate(fname)
        end do
    end subroutine merge_proj_weights_from_parts

end module simple_map_reduce
