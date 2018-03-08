! concrete commander: routines for managing distributed SIMPLE execution
module simple_commander_distr
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_oris,           only: oris
use simple_binoris,        only: binoris
implicit none

public :: merge_algndocs_commander
public :: merge_nnmat_commander
public :: merge_similarities_commander
public :: split_pairs_commander
private

type, extends(commander_base) :: merge_algndocs_commander
  contains
    procedure :: execute      => exec_merge_algndocs
end type merge_algndocs_commander
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

contains

    !> for merging alignment documents from SIMPLE runs in distributed mode
    subroutine exec_merge_algndocs( self, cline )
        use simple_map_reduce, only: split_nobjs_even
        use simple_sp_project, only: transfer_sp_project_segment
        class(merge_algndocs_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(params)               :: p
        type(str4arr), allocatable :: os_strings(:)
        integer,       allocatable :: parts(:,:)
        integer                    :: i, j, nj, numlen, funit, funit_merge, n_records
        integer                    :: io_stat, partsz, cnt, fromto(2), isegment, strlen, strlen_max
        type(binoris)              :: bos_doc, bos_merged
        character(len=STDLEN)      :: fname
        character(len=1024)        :: line
        p      = params(cline) ! parameters generated
        !$ allocate(parts(p%nptcls,2))
        parts  = split_nobjs_even(p%nptcls, p%ndocs)
        if( cline%defined('numlen') )then
            numlen = p%numlen
        else
            numlen = len(int2str(p%ndocs))
        endif
        ! if projfile is present, it will be both the input and output
        ! to preserve meta-data information
        if( cline%defined('projfile') )then
            if( cline%defined('outfile') )then
                write(*,*) 'ERROR! outfile cannot be defined simultaneously with projfile'
                stop 'commander_distr :: exec_merge_algndocs'
            endif
            p%outfile = 'tmpfile_from_merge_algndocs.simple'
        endif
        select case(trim(METADATA_EXT))
            case('.simple')
                ! allocate merged string representation
                allocate( os_strings(p%nptcls) )
                ! convert from flag to enumerator to integer
                select case(trim(p%oritype))
                    case('mic')
                        isegment = MIC_SEG
                    case('stk')
                        isegment = STK_SEG
                    case('ptcl2D')
                        isegment = PTCL2D_SEG
                    case('cls3D')
                        isegment = CLS3D_SEG
                    case('ptcl3D')
                        isegment = PTCL3D_SEG
                    case DEFAULT
                        isegment = GENERIC_SEG
                end select
                ! read into string representation
                do i=1,p%ndocs
                    fname     = trim(adjustl(p%fbody))//int2str_pad(i,numlen)//trim(METADATA_EXT)
                    call bos_doc%open(trim(fname))
                    n_records = bos_doc%get_n_records(isegment)
                    fromto    = bos_doc%get_fromto(isegment)
                    partsz    = parts(i,2) - parts(i,1) + 1
                    if( n_records /= partsz .or. .not. all(fromto == parts(i,:)) )then
                        write(*,*) 'ERROR, # records does not match expectation'
                        write(*,*) 'EXTRACTED FROM file: ', trim(fname)
                        write(*,*) 'fromto   : ', fromto(1), fromto(2)
                        write(*,*) 'n_records: ', n_records
                        write(*,*) 'CALCULATED FROM input p%nptcls/p%ndocs'
                        write(*,*) 'fromto: ', parts(i,1), parts(i,2)
                        write(*,*) 'partsz: ', partsz
                        stop
                    endif
                    call bos_doc%read_segment(isegment, os_strings)
                    call bos_doc%close
                end do
                ! find maxium string lenght
                strlen_max = 0
                !$omp parallel do schedule(static) default(shared) proc_bind(close)&
                !$omp private(i,strlen) reduction(max:strlen_max)
                do i=1,p%nptcls
                    strlen = len_trim(os_strings(i)%str)
                    if( strlen > strlen_max ) strlen_max = strlen
                end do
                !$omp end parallel do
                ! write as one (merged) segment
                call bos_merged%open(p%outfile, del_if_exists=.true.)
                call bos_merged%write_segment(isegment, [1,p%nptcls], strlen_max, os_strings)
                call bos_merged%write_header
                call bos_merged%close
                ! better be explicit about deallocating the derived type
                do i=1,p%nptcls
                    deallocate(os_strings(i)%str)
                end do
                deallocate(os_strings)
                if( cline%defined('projfile') )then
                    call transfer_sp_project_segment(p%outfile, p%projfile, trim(p%oritype))
                    call del_file(p%outfile)
                endif
            case('.txt')
                call fopen(funit_merge, file=p%outfile, iostat=io_stat, status='replace',&
                &action='write', position='append', access='sequential')
                call fileiochk("Error opening file"//trim(adjustl(p%outfile)), io_stat)
                do i=1,p%ndocs
                    fname = trim(adjustl(p%fbody))//int2str_pad(i,numlen)//trim(METADATA_EXT)
                    nj = nlines(fname)
                    partsz = parts(i,2) - parts(i,1) + 1
                    if( partsz /= nj ) then
                        write(*,*) 'nr of entries in partition: ', partsz
                        write(*,*) 'nr of lines in file: ', nj
                        write(*,*) 'filename: ', trim(fname)
                        stop 'number of lines in file not consistent with the size of the partition'
                    endif
                    call fopen(funit, file=fname, iostat=io_stat, status='old', action='read', access='sequential')
                    call fileiochk("Error opening file "//trim(adjustl(fname)), io_stat)
                    do j=1,nj
                        read(funit,fmt='(A)') line
                        write(funit_merge,fmt='(A)') trim(line)
                    end do
                    call fclose(funit, iostat=io_stat,errmsg="Error closing file ")
                end do
                call fclose(funit_merge, iostat=io_stat,errmsg="Error closing outfile ")
        end select
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_ALGNDOCS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_merge_algndocs

    !> for merging partial nearest neighbour matrices calculated in distributed mode
    subroutine exec_merge_nnmat( self, cline )
        use simple_map_reduce, only: merge_nnmat_from_parts
        class(merge_nnmat_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)         :: p
        integer, allocatable :: nnmat(:,:)
        integer :: filnum, io_stat
        p      = params(cline) ! parameters generated
        nnmat  = merge_nnmat_from_parts(p%nptcls, p%nparts, p%nnn)            !! intel realloc warning
        call fopen(filnum, status='REPLACE', action='WRITE', file='nnmat.bin', access='STREAM', iostat=io_stat)
        call fileiochk('simple_merge_nnmat ; fopen error when opening nnmat.bin  ', io_stat)
        write(unit=filnum,pos=1,iostat=io_stat) nnmat
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to nnmat.bin'
            stop 'I/O error; simple_merge_nnmat'
        endif
        call fclose(filnum,errmsg='simple_merge_nnmat ;error when closing nnmat.bin  ')
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
        simmat = merge_similarities_from_parts(p%nptcls, p%nparts)           !! intel realloc warning
        call fopen(filnum, status='REPLACE', action='WRITE', file='smat.bin', access='STREAM', iostat=io_stat)
        call fileiochk('simple_merge_nnmat ; fopen error when opening smat.bin  ', io_stat)
        write(unit=filnum,pos=1,iostat=io_stat) simmat
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to smat.bin'
            stop 'I/O error; simple_merge_similarities'
        endif
        call fclose(filnum,errmsg='simple_merge_nnmat ; error when closing smat.bin ')
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_SIMILARITIES NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_merge_similarities

    !> for splitting calculations between pairs of objects into balanced partitions
    subroutine exec_split_pairs( self, cline )
        use simple_map_reduce, only: split_pairs_in_parts
        class(split_pairs_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params) :: p
        p = params(cline) ! parameters generated
        call split_pairs_in_parts(p%nptcls, p%nparts)
        call simple_end('**** SIMPLE_SPLIT_PAIRS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_split_pairs

end module simple_commander_distr
