! concrete commander: routines for managing distributed SIMPLE execution
module simple_commander_distr
use simple_defs
use simple_syslib,         only: alloc_errchk
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_strings,        only: int2str, int2str_pad
use simple_fileio,         only: fopen,fclose, fileio_errmsg
use simple_imgfile,        only: find_ldim_nptcls
use simple_jiffys,         only: simple_end, progress
implicit none

public :: merge_algndocs_commander
public :: merge_binalgndocs_commander
public :: merge_nnmat_commander
public :: merge_similarities_commander
public :: split_pairs_commander
public :: split_commander
private

!> generator type
type, extends(commander_base) :: merge_algndocs_commander
  contains
    procedure :: execute      => exec_merge_algndocs
end type merge_algndocs_commander
type, extends(commander_base) :: merge_binalgndocs_commander
  contains
    procedure :: execute      => exec_merge_binalgndocs
end type merge_binalgndocs_commander
type, extends(commander_base) :: merge_nnmat_commander
  contains
    procedure :: execute      => exec_merge_nnmat
end type merge_nnmat_commander
type, extends(commander_base) :: merge_similarities_commander
  contains
    procedure :: execute      => exec_merge_similarities
end type merge_similarities_commander
type, extends(commander_base) :: split_pairs_commander
  contains
    procedure :: execute      => exec_split_pairs
end type split_pairs_commander
type, extends(commander_base) :: split_commander
  contains
    procedure :: execute      => exec_split
end type split_commander

contains

    !> for merging alignment documents from SIMPLE runs in distributed mode
    subroutine exec_merge_algndocs( self, cline )
        use simple_oris, only: oris
        use simple_map_reduce ! use all in there
        class(merge_algndocs_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(params)          :: p
        integer               :: i, j, nj, numlen, funit, funit_merge, io_stat, nentries_all
        character(len=STDLEN) :: fname
        integer, allocatable  :: parts(:,:)
        character(len=1024)   :: line
        p = params(cline) ! parameters generated
        parts = split_nobjs_even(p%nptcls, p%ndocs)
        if(.not.fopen(funit_merge, file=p%outfile, iostat=io_stat, status='replace',&
             &action='write', position='append', access='sequential'))then
            call fileio_errmsg("Error opening file"//trim(adjustl(p%outfile)), io_stat)
        endif
        numlen = len(int2str(p%ndocs))
        
        do i=1,p%ndocs
            fname = trim(adjustl(p%fbody))//int2str_pad(i,numlen)//'.txt'
            nj = nlines(fname)
            nentries_all = parts(i,2) - parts(i,1) + 1
            if( nentries_all /= nj ) then
                write(*,*) 'nr of entries in partition: ', nentries_all
                write(*,*) 'nr of lines in file: ', nj
                write(*,*) 'filename: ', trim(fname)
                stop 'number of lines in file not consistent with the size of the partition'
            endif
            if(.not.fopen(funit, file=fname, iostat=io_stat, status='old', action='read', access='sequential'))then
                call fileio_errmsg("Error opening file "//trim(adjustl(fname)), io_stat)
                stop
            endif
            do j=1,nj
                read(funit,fmt='(A)') line
                write(funit_merge,fmt='(A)') trim(line)
            end do
            if(.not.fclose(funit, iostat=io_stat))&
                 call fileio_errmsg("Error closing file ", io_stat)
        end do
        if(.not.fclose(funit_merge, iostat=io_stat))&
             call fileio_errmsg("Error closing outfile ", io_stat)
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_ALGNDOCS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_merge_algndocs

    !> for merging binary alignment documents from SIMPLE runs in distributed mode
    subroutine exec_merge_binalgndocs( self, cline )
        use simple_oris,    only: oris
        use simple_binoris, only: binoris
        use simple_map_reduce ! use all in there
        class(merge_binalgndocs_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(params)          :: p
        integer               :: i, j, nj, numlen, nentries_all, cnt, fromto(2)
        character(len=STDLEN) :: fname
        integer, allocatable  :: parts(:,:)
        type(binoris)         :: fhandle_doc, fhandle_merged
        p      = params(cline) ! parameters generated
        parts  = split_nobjs_even(p%nptcls, p%ndocs)
        numlen = len(int2str(p%ndocs))
        ! generate a merged filehandling object based on the first file
        fname = trim(adjustl(p%fbody))//int2str_pad(1,numlen)//'.bin'
        call fhandle_merged%open(fname)
        call fhandle_merged%set_fromto([1,p%nptcls])
        call fhandle_merged%close()
        ! make a new header based on the modified template
        call fhandle_merged%open(p%outfile, del_if_exists=.true.)
        call fhandle_merged%write_header()
        ! loop over documents
        cnt = 0
        do i=1,p%ndocs
            fname = trim(adjustl(p%fbody))//int2str_pad(i,numlen)//'.bin'
            call fhandle_doc%open(fname)
            nj = fhandle_doc%get_n_records()
            nentries_all = parts(i,2) - parts(i,1) + 1
            if( nentries_all /= nj ) then
                write(*,*) 'in: commander_distr :: exec_merge_binalgndocs'
                write(*,*) 'nr of entries in partition: ', nentries_all
                write(*,*) 'nr of records in file     : ', nj
                write(*,*) 'filename                  : ', trim(fname)
                stop 'number of records in file not consistent with the size of the partition'
            endif
            fromto = fhandle_doc%get_fromto()
            if( any(fromto /= parts(i,:)) )then
                write(*,*) 'in: commander_distr :: exec_merge_binalgndocs'
                write(*,*) 'range in partition: ', parts(i,:)
                write(*,*) 'range in file     : ', fromto
                stop 'range in file not consistent with the range in the partition'
            endif
            do j=fromto(1),fromto(2)
                cnt = cnt + 1
                call fhandle_doc%read_record(j)
                call fhandle_merged%write_record(cnt, fhandle_doc)
            end do
        end do
        call fhandle_merged%close()
        call fhandle_doc%close()
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_BINALGNDOCS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_merge_binalgndocs

    !> merge_nnmat is a program for merging partial nearest neighbour matrices calculated in distributed mode
    subroutine exec_merge_nnmat( self, cline )
        use simple_map_reduce, only: merge_nnmat_from_parts
        class(merge_nnmat_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)         :: p
        integer, allocatable :: nnmat(:,:)
        integer :: filnum, io_stat
        p      = params(cline) ! parameters generated
        nnmat  = merge_nnmat_from_parts(p%nptcls, p%nparts, p%nnn)
        if(.not.fopen(filnum, status='REPLACE', action='WRITE', file='nnmat.bin', access='STREAM', iostat=io_stat))&
             call fileio_errmsg('simple_merge_nnmat ; fopen error when opening nnmat.bin  ', io_stat)
        write(unit=filnum,pos=1,iostat=io_stat) nnmat
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to nnmat.bin'
            stop 'I/O error; simple_merge_nnmat'
        endif
        if(.not.fclose(filnum, iostat=io_stat))&
             call fileio_errmsg('simple_merge_nnmat ; fopen error when closing nnmat.bin  ', io_stat)
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_NNMAT NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_merge_nnmat

    subroutine exec_merge_similarities( self, cline )
        use simple_map_reduce, only: merge_similarities_from_parts
        class(merge_similarities_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(params)      :: p
        real, allocatable :: simmat(:,:)
        integer           :: filnum, io_stat
        p      = params(cline) ! parameters generated
        simmat = merge_similarities_from_parts(p%nptcls, p%nparts)
       
        if(.not.fopen(filnum, status='REPLACE', action='WRITE', file='smat.bin', access='STREAM', iostat=io_stat))&
             call fileio_errmsg('simple_merge_nnmat ; fopen error when opening smat.bin  ', io_stat)
        write(unit=filnum,pos=1,iostat=io_stat) simmat
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to smat.bin'
            stop 'I/O error; simple_merge_similarities'
        endif
         if(.not.fclose(filnum, iostat=io_stat))&
             call fileio_errmsg('simple_merge_nnmat ; fopen error when opening smat.bin  ', io_stat)
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_SIMILARITIES NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_merge_similarities

    !> split_pairs is a program for splitting calculations between pairs of objects into balanced partitions
    subroutine exec_split_pairs( self, cline )
        use simple_map_reduce, only: split_pairs_in_parts
        class(split_pairs_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params) :: p
        p = params(cline) ! parameters generated
        call split_pairs_in_parts(p%nptcls, p%nparts)
        call simple_end('**** SIMPLE_SPLIT_PAIRS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_split_pairs
    
    !> split is a program for splitting of image stacks into partitions for parallel execution.
    !! This is done to reduce I/O latency
    subroutine exec_split( self, cline )
        use simple_map_reduce ! use all in there
        use simple_image,  only: image
        use simple_jiffys, only: progress
        class(split_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(params)         :: p
        type(image)          :: img
        integer              :: iptcl, ipart, ldim(3), cnt, nimgs
        integer, allocatable :: parts(:,:)
        logical              :: either_defined
        p = params(cline) ! parameters generated
        call find_ldim_nptcls(p%stk, ldim, nimgs)
        ldim(3) = 1
        call img%new(ldim,1.)
        parts = split_nobjs_even(nimgs, p%nparts)
        if( size(parts,1) /= p%nparts ) stop 'ERROR! generated number of parts not same as inputted nparts'
        if( .not. stack_is_split() )then
            call exec_cmdline('mkdir -p '//trim(STKPARTSDIR))
            do ipart=1,p%nparts
                call progress(ipart,p%nparts)
                cnt = 0
                do iptcl=parts(ipart,1),parts(ipart,2)
                    cnt = cnt+1
                    call img%read(p%stk, iptcl)
                    if( p%neg .eq. 'yes' ) call img%neg
                    call img%write(trim(STKPARTFBODY)//int2str_pad(ipart,p%numlen)//p%ext, cnt)
                end do
            end do
        else
            write(*,'(a)') '>>> NO SPLIT PERFORMED, STACK ALREADY PARTITIONED CORRECTLY'
        endif
        deallocate(parts)
        call img%kill
        call simple_end('**** SIMPLE_SPLIT NORMAL STOP ****', print_simple=.false.)

        contains

            logical function stack_is_split()
                character(len=:), allocatable :: stack_part_fname
                logical,          allocatable :: stack_parts_exist(:) 
                integer :: ipart, numlen, sz, sz_correct, ldim(3)
                logical :: is_split, is_correct, is_split_correctly
                allocate( stack_parts_exist(p%nparts) )
                numlen = len(int2str(p%nparts))
                do ipart=1,p%nparts
                    allocate(stack_part_fname, source=trim(STKPARTFBODY)//int2str_pad(ipart,numlen)//p%ext)
                    stack_parts_exist(ipart) = file_exists(stack_part_fname)
                    deallocate(stack_part_fname)
                end do
                is_split   = all(stack_parts_exist)
                is_correct = .true.
                if( is_split )then
                    do ipart=1,p%nparts
                        sz_correct = parts(ipart,2)-parts(ipart,1)+1
                        allocate(stack_part_fname, source=trim(STKPARTFBODY)//int2str_pad(ipart,numlen)//p%ext)
                        call find_ldim_nptcls(stack_part_fname, ldim, sz)
                        if( sz /= sz_correct )then
                            is_correct = .false.
                            exit
                        endif
                        if( ldim(1) == p%box .and. ldim(2) == p%box )then
                            ! dimension ok
                        else
                            is_correct = .false.
                            exit
                        endif
                        deallocate(stack_part_fname)
                    end do
                endif
                stack_is_split = is_split .and. is_correct
                if( .not. stack_is_split ) call del_files(trim(STKPARTFBODY), p%nparts, ext=p%ext)
            end function stack_is_split

    end subroutine exec_split

end module simple_commander_distr
