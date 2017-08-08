!------------------------------------------------------------------------------!
! SIMPLE v3.0         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Simple map reduce module:  distributed resource allocation
module simple_map_reduce
use simple_defs
use simple_strings,      only: int2str, int2str_pad
use simple_filehandling, only: get_fileunit
use simple_jiffys,       only: alloc_err
implicit none

#include "simple_local_flags.inc"
public ::  split_nobjs_even,split_nobjs_in_chunks,split_pairs_in_parts,merge_similarities_from_parts,merge_rmat_from_parts
contains
    
    !>  \brief  for generating balanced partitions of nobjs objects
    function split_nobjs_even( nobjs, nparts ) result( parts )
        integer, intent(in)  :: nobjs, nparts
        integer, allocatable :: parts(:,:)
        integer :: nobjs_per_part, leftover, istop, istart, alloc_stat, ipart
        allocate(parts(nparts,2), stat=alloc_stat)
        call alloc_err('In: simple_map_reduce :: split_nobjs_even', alloc_stat)
        nobjs_per_part = nobjs/nparts
        DebugPrint   'nobjs_per_part: ', nobjs_per_part
        leftover = nobjs-nobjs_per_part*nparts
        DebugPrint   'leftover: ', leftover
        istop  = 0
        istart = 0
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
        end do
    end function split_nobjs_even
    
    !>  \brief  for generating partitions of nobjs objects of pre-determined size
    function split_nobjs_in_chunks( nobjs, chunksz ) result( parts )
        integer, intent(in)  :: nobjs, chunksz
        integer, allocatable :: parts(:,:)
        integer :: alloc_stat, ipart, nparts, cnt
        logical :: leftover
        ! count the number of parts
        cnt   = 0
        nparts = 0
        do 
            cnt = cnt+chunksz
            if( cnt <= nobjs )then
                nparts = nparts+1
                cycle
            else
                exit
            endif
        end do
        leftover = .false.
        if( chunksz*nparts == nobjs )then
            ! the objects were evenly partitioned, i.e. mod(nobjs,chunksz) .eq. 0
        else
            ! there's a leftover
            leftover = .true.
            ! and we need to increment the number of partitions with 1
            nparts = nparts+1
        endif
        ! generate output
        allocate(parts(nparts,2), stat=alloc_stat)
        call alloc_err('In: simple_map_reduce :: split_nobjs_in_chunks', alloc_stat)
        ! count the number of parts
        cnt   = 0
        ipart = 0
        do 
            cnt = cnt+chunksz
            if( cnt <= nobjs )then
                ipart = ipart+1
                parts(ipart,1) = cnt-chunksz+1
                parts(ipart,2) = cnt
                cycle
            else
                exit
            endif
        end do
        if( leftover )then
            parts(ipart+1,1) = parts(ipart,2)+1
            parts(ipart+1,2) = nobjs
        endif
    end function split_nobjs_in_chunks
    
    !>  \brief  for generating balanced partitions for pairwise calculations on nobjs ojects
    subroutine split_pairs_in_parts( nobjs, nparts )
        use simple_jiffys, only: progress
        integer, intent(in)  :: nobjs  !< number objects to analyse in pairs
        integer, intent(in)  :: nparts !< number of partitions (nodes) for parallel execution
        integer              :: npairs, alloc_stat, cnt, funit, i, j
        integer              :: ipart, io_stat, numlen
        integer, allocatable :: pairs(:,:), parts(:,:)
        character(len=:), allocatable :: fname
        
        ! generate all pairs
        DebugPrint  ' split_pairs_in_parts nobjs: ', nobjs
        npairs = (nobjs*(nobjs-1))/2
        DebugPrint  'split_pairs_in_parts, npairs: ', npairs
        allocate( pairs(npairs,2), stat=alloc_stat )
        cnt = 0
        do i=1,nobjs-1
            do j=i+1,nobjs
                cnt = cnt+1
                pairs(cnt,1) = i
                pairs(cnt,2) = j
            end do
        end do
        ! generate balanced partitions of pairs
        parts  = split_nobjs_even(npairs, nparts)
        numlen = len(int2str(nparts))
        ! write the partitions
        do ipart=1,nparts
            call progress(ipart,nparts)
            funit = get_fileunit()
            allocate(fname, source='pairs_part' // int2str_pad(ipart,numlen) // '.bin')
            open(unit=funit, status='REPLACE', action='WRITE', file=fname, access='STREAM')
            DebugPrint   'writing pairs in range: ', parts(ipart,1), parts(ipart,2)
            write(unit=funit,pos=1,iostat=io_stat) pairs(parts(ipart,1):parts(ipart,2),:)
            ! Check if the write was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(split_pairs_in_parts): I/O error ',&
                io_stat, ' when writing to: ', fname
                stop 'I/O error; split_pairs_in_part; simple_map_reduce'
            endif
            close(funit)
            deallocate(fname)
        end do
        deallocate(pairs, parts)
    end subroutine split_pairs_in_parts
    
    !>  \brief  for merging partial calculations of similarities
    function merge_similarities_from_parts( nobjs, nparts ) result( smat )
        integer, intent(in) :: nobjs  !< number objects to analyse in pairs
        integer, intent(in) :: nparts !< number of partitions (nodes) for parallel execution
        integer :: npairs, alloc_stat, funit, ipart, i, io_stat, numlen
        real,    allocatable :: smat(:,:), sims(:)
        integer, allocatable :: pairs(:,:), parts(:,:)
        character(len=:), allocatable :: fname
        ! allocate pairs and similarities
        npairs = (nobjs*(nobjs-1))/2
        DebugPrint   'analysing this number of objects: ', nobjs
        DebugPrint   'analysing this number of pairs: ', npairs
        allocate( smat(nobjs,nobjs), pairs(npairs,2), stat=alloc_stat )
        call alloc_err('In: simple_map_reduce::merge_similarities_from_parts, 1', alloc_stat)
        ! initialise similarities
        smat = -1.
        do i=1,nobjs
            smat(i,i) = 1
        end do
        ! repeat the generation of balanced partitions
        parts  = split_nobjs_even(npairs, nparts)
        numlen = len(int2str(nparts))
        ! compress the partial similarities into a single matrix
        do ipart=1,nparts
            ! read the index mapping
            funit = get_fileunit()
            allocate(fname, source='pairs_part'//int2str_pad(ipart,numlen)//'.bin')
            open(unit=funit, status='OLD', action='READ', file=fname, access='STREAM')
            read(unit=funit,pos=1,iostat=io_stat) pairs(parts(ipart,1):parts(ipart,2),:)
            ! Check if the read was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(merge_similarities_from_parts): I/O error ',&
                io_stat, ' when reading file: ', fname
                stop 'I/O error; merge_similarities_from_parts; simple_map_reduce'
            endif
            close(funit)
            deallocate(fname)
            ! retrieve the similarities
            DebugPrint   'allocating this number of similarities: ', parts(ipart,2)-parts(ipart,1)+1
            allocate(sims(parts(ipart,1):parts(ipart,2)), stat=alloc_stat)
            call alloc_err('In: simple_map_reduce::merge_similarities_from_parts, 2', alloc_stat)
            allocate(fname, source='similarities_part'//int2str_pad(ipart,numlen)//'.bin')
            open(unit=funit, status='OLD', action='READ', file=fname, access='STREAM')
            read(unit=funit,pos=1,iostat=io_stat) sims(parts(ipart,1):parts(ipart,2))
            ! check if the read was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(merge_similarities_from_parts): I/O error ',&
                io_stat, ' when reading: ', fname
                stop 'I/O error; merge_similarities_from_parts; simple_map_reduce'
            endif
            close(funit)
            deallocate(fname)
            ! set the similarity matrix components
            do i=parts(ipart,1),parts(ipart,2)
                smat(pairs(i,1),pairs(i,2)) = sims(i)
                smat(pairs(i,2),pairs(i,1)) = sims(i)
            end do
            deallocate(sims)
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
        integer :: alloc_stat, numlen, ipart, funit, io_stat
        ! allocate nearest neighbour matrix
        allocate(nnmat(nobjs,nnn), stat=alloc_stat)
        call alloc_err("In: simple_map_reduce :: merge_nnmat_from_parts", alloc_stat)
        ! repeat the generation of balanced partitions
        parts  = split_nobjs_even(nobjs, nparts)
        numlen = len(int2str(nparts))
        ! compress the partial nearest neighbour matrices into a single matrix
        funit = get_fileunit()
        do ipart=1,nparts
            allocate(fname, source='nnmat_part'//int2str_pad(ipart,numlen)//'.bin')
            open(unit=funit, status='OLD', action='READ', file=fname, access='STREAM')
            read(unit=funit,pos=1,iostat=io_stat) nnmat(parts(ipart,1):parts(ipart,2),:)
            ! check if the read was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(merge_nnmat_from_parts): I/O error ',&
                io_stat, ' when reading: ', fname
                stop 'I/O error; merge_nnmat_from_parts; simple_map_reduce'
            endif
            close(funit)
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
        integer :: alloc_stat, numlen, ipart, funit, io_stat
        ! allocate merged matrix
        allocate(mat_merged(nobjs,ydim), stat=alloc_stat)
        call alloc_err("In: simple_map_reduce :: merge_mat_from_parts", alloc_stat)
        ! repeat the generation of balanced partitions
        parts  = split_nobjs_even(nobjs, nparts)
        numlen = len(int2str(nparts))
        ! compress the partial matrices into a single matrix
        funit = get_fileunit()
        do ipart=1,nparts
            allocate(fname, source=trim(fbody)//int2str_pad(ipart,numlen)//'.bin')
            open(unit=funit, status='OLD', action='READ', file=fname, access='STREAM')
            read(unit=funit,pos=1,iostat=io_stat) mat_merged(parts(ipart,1):parts(ipart,2),:)
            ! check if the read was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(merge_rmat_from_parts_1): I/O error ',&
                io_stat, ' when reading: ', fname
                stop 'I/O error; merge_rmat_from_parts_1; simple_map_reduce'
            endif
            close(funit)
            deallocate(fname)
        end do
        deallocate(parts)
    end function merge_rmat_from_parts
    
end module simple_map_reduce
