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
public ::  split_nobjs_even,split_pairs_in_parts,merge_similarities_from_parts

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

    !>  \brief  for generating balanced partitions for pairwise calculations on nobjs ojects
    subroutine split_pairs_in_parts( nobjs, nparts )
        integer, intent(in)  :: nobjs  !< number objects to analyse in pairs
        integer, intent(in)  :: nparts !< number of partitions (nodes) for parallel execution
        integer              :: npairs, cnt, funit, i, j
        integer              :: ipart, io_stat, numlen
        integer, allocatable :: pairs(:,:), parts(:,:)
        character(len=:), allocatable :: fname
        ! generate all pairs
        npairs = (nobjs*(nobjs-1))/2
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

end module simple_map_reduce
