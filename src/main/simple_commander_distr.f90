! concrete commander: routines for managing distributed SIMPLE execution
module simple_commander_distr
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: parameters
implicit none

public :: merge_nnmat_commander
public :: merge_similarities_commander
public :: split_pairs_commander
public :: split_commander
private
#include "simple_local_flags.inc"

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

    !> for merging partial nearest neighbour matrices calculated in distributed mode
    subroutine exec_merge_nnmat( self, cline )
        use simple_map_reduce, only: merge_nnmat_from_parts
        class(merge_nnmat_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters) :: params
        integer, allocatable :: nnmat(:,:)
        integer :: filnum, io_stat
        call params%new(cline)
        nnmat  = merge_nnmat_from_parts(params%nptcls, params%nparts, params%nnn) !! intel realloc warning
        call fopen(filnum, status='REPLACE', action='WRITE', file='nnmat.bin', access='STREAM', iostat=io_stat)
        call fileiochk('simple_merge_nnmat ; fopen error when opening nnmat.bin  ', io_stat)
        write(unit=filnum,pos=1,iostat=io_stat) nnmat
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to nnmat.bin'
            THROW_HARD('I/O; exec_merge_nnmat')
        endif
        call fclose(filnum,errmsg='simple_merge_nnmat ;error when closing nnmat.bin  ')
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_NNMAT NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_merge_nnmat

    subroutine exec_merge_similarities( self, cline )
        use simple_map_reduce, only: merge_similarities_from_parts
        class(merge_similarities_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters)  :: params
        real, allocatable :: simmat(:,:)
        integer           :: filnum, io_stat
        call params%new(cline)
        simmat = merge_similarities_from_parts(params%nptcls, params%nparts) !! intel realloc warning
        call fopen(filnum, status='REPLACE', action='WRITE', file='smat.bin', access='STREAM', iostat=io_stat)
        call fileiochk('simple_merge_nnmat ; fopen error when opening smat.bin  ', io_stat)
        write(unit=filnum,pos=1,iostat=io_stat) simmat
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to smat.bin'
            THROW_HARD('I/O; exec_merge_similarities')
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
        type(parameters) :: params
        call params%new(cline)
        call split_pairs_in_parts(params%nptcls, params%nparts)
        call simple_end('**** SIMPLE_SPLIT_PAIRS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_split_pairs

    !> split is a program for splitting of image stacks into partitions for parallel execution.
    !! This is done to reduce I/O latency
    subroutine exec_split( self, cline )
        use simple_map_reduce ! use all in there
        use simple_image, only: image
        class(split_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(parameters)     :: p
        type(image)          :: img
        integer              :: iptcl, ipart, ldim(3), cnt, nimgs
        integer, allocatable :: parts(:,:)
        call p%new(cline)
        call find_ldim_nptcls(p%stk, ldim, nimgs)
        ldim(3) = 1
        call img%new(ldim,p%smpd)
        parts = split_nobjs_even(nimgs, p%nparts)
        if( size(parts,1) /= p%nparts ) THROW_HARD('generated number of parts not same as inputted nparts')
        do ipart=1,p%nparts
            call progress(ipart,p%nparts)
            cnt = 0
            do iptcl=parts(ipart,1),parts(ipart,2)
                cnt = cnt+1
                call img%read(p%stk, iptcl)
                call img%write('stack_part'//int2str_pad(ipart,p%numlen)//p%ext, cnt)
            end do
        end do
        deallocate(parts)
        call img%kill
        call simple_end('**** SIMPLE_SPLIT NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_split

end module simple_commander_distr
