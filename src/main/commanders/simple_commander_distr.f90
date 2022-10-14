! concrete commander: routines for managing distributed SIMPLE execution
module simple_commander_distr
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: parameters
use simple_stack_io,       only: stack_io
implicit none

public :: split_pairs_commander
public :: split_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: split_pairs_commander
  contains
    procedure :: execute      => exec_split_pairs
end type split_pairs_commander

type, extends(commander_base) :: split_commander
  contains
    procedure :: execute      => exec_split
end type split_commander

contains

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
        type(stack_io)       :: stkio_r, stkio_w
        integer              :: iptcl, ipart, ldim(3), cnt, nimgs
        integer, allocatable :: parts(:,:)
        call p%new(cline)
        call find_ldim_nptcls(p%stk, ldim, nimgs)
        ldim(3) = 1
        call img%new(ldim,p%smpd)
        parts = split_nobjs_even(nimgs, p%nparts)
        if( size(parts,1) /= p%nparts ) THROW_HARD('generated number of parts not same as inputted nparts')
        call stkio_r%open(p%stk, p%smpd, 'read')
        do ipart=1,p%nparts
            call progress(ipart,p%nparts)
            call stkio_w%open('stack_part'//int2str_pad(ipart,p%numlen)//p%ext, p%smpd, 'write', box=p%box, is_ft=.false.)
            cnt = 0
            do iptcl=parts(ipart,1),parts(ipart,2)
                cnt = cnt + 1
                call stkio_r%read(iptcl, img)
                call stkio_w%write(cnt, img)
            end do
            call stkio_w%close
        end do
        call stkio_r%close
        deallocate(parts)
        call img%kill
        call simple_end('**** SIMPLE_SPLIT NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_split

end module simple_commander_distr
