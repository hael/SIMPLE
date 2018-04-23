! concrete commander: routines for managing distributed SIMPLE execution
module simple_commander_distr
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_commander_base, only: commander_base
implicit none

public :: merge_nnmat_commander
public :: merge_similarities_commander
public :: split_pairs_commander
private

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
