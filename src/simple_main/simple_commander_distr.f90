!==Class simple_commander_distr
!
! This class contains the set of concrete distr commanders of the SIMPLE library used to provide pre/post processing routines
! for SIMPLE when executed in distributed mode. This class provides the glue between the reciver (main reciever is simple_exec 
! program) and the abstract action, which is simply execute (defined by the base class: simple_commander_base). 
! Later we can use the composite pattern to create MacroCommanders (or workflows)
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Cyril Reboul & Hans Elmlund 2016
!
module simple_commander_distr
use simple_defs
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_strings,        only: int2str, int2str_pad
use simple_filehandling    ! use all in there
use simple_jiffys          ! use all in there
implicit none

public :: merge_algndocs_commander
public :: merge_nnmat_commander
public :: merge_shellweights_commander
public :: merge_similarities_commander
public :: split_pairs_commander
public :: split_commander
private

type, extends(commander_base) :: merge_algndocs_commander
  contains
    procedure :: execute      => exec_merge_algndocs
end type merge_algndocs_commander
type, extends(commander_base) :: merge_nnmat_commander
  contains
    procedure :: execute      => exec_merge_nnmat
end type merge_nnmat_commander
type, extends(commander_base) :: merge_shellweights_commander
  contains
    procedure :: execute      => exec_merge_shellweights
end type merge_shellweights_commander
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
    
    subroutine exec_merge_algndocs( self, cline )
        use simple_oris, only: oris
        use simple_map_reduce ! use all in there
        class(merge_algndocs_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(params)          :: p
        type(oris)            :: o, o_read
        integer               :: i, j, nentries, cnt, nentries_all, numlen
        integer, allocatable  :: parts(:,:)
        character(len=STDLEN) :: fname
        logical               :: here
        logical, parameter    :: print_part_info = .false.
        p = params(cline) ! parameters generated
        select case(p%split_mode)
            case('chunk')
                parts = split_nobjs_in_chunks(p%nptcls, p%chunksz)
            case('even')
                parts = split_nobjs_even(p%nptcls, p%ndocs)
            case DEFAULT
                write(*,*) 'split mode: ', trim(p%split_mode)
                stop 'unsupported split_mode; simple_commander_distr :: exec_merge_algndocs'
        end select
        numlen = len(int2str(p%ndocs))
        do i=1,p%ndocs
            fname = trim(adjustl(p%fbody))//int2str_pad(i,numlen)//'.txt'
            ! calculate the number of all entries
            nentries_all = parts(i,2)-parts(i,1)+1 
            ! calculate the actual number of entries
            inquire(FILE=fname, EXIST=here)
            if( here )then
                nentries = nlines(fname)
            else
                nentries = 0
            endif
            ! check if oritab is there to fill-in blanks
            if( nentries < nentries_all )then
                write(*,*) 'nentries: ',     nentries
                write(*,*) 'nentries_all: ', nentries_all
                write(*,*) 'number of entries in the doc to be merged not equal to the expected number'
                stop 'simple_commander_distr :: exec_merge_algndocs'
            endif
            if( print_part_info )then
                write(*,'(a,1x,i3,1x,a,1x,i6,1x,i6)') 'partition:', i, 'from/to:', parts(i,1), parts(i,2)
            endif
            o_read = oris(nentries)
            call o_read%read(fname)
            call o%merge(o_read)
        end do
        call o%write(p%outfile)
        deallocate(parts)
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_ALGNDOCS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_merge_algndocs

    subroutine exec_merge_nnmat( self, cline )
        use simple_map_reduce, only: merge_nnmat_from_parts
        class(merge_nnmat_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)         :: p
        integer, allocatable :: nnmat(:,:)
        integer :: filnum, io_stat
        p      = params(cline) ! parameters generated
        nnmat  = merge_nnmat_from_parts(p%nptcls, p%nparts, p%nnn)
        filnum = get_fileunit()
        open(unit=filnum, status='REPLACE', action='WRITE', file='nnmat.bin', access='STREAM')
        write(unit=filnum,pos=1,iostat=io_stat) nnmat
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to nnmat.bin'
            stop 'I/O error; simple_merge_nnmat'
        endif
        close(filnum)
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_NNMAT NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_merge_nnmat

    subroutine exec_merge_shellweights( self, cline )
        use simple_map_reduce, only: merge_rmat_from_parts
        use simple_filterer,   only: normalise_shellweights
        class(merge_shellweights_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(params)      :: p
        type(build)       :: b
        real, allocatable :: wmat(:,:), wmat_states(:,:,:)
        integer :: filnum, io_stat, filtsz
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built (assumes stk input)
        filtsz = b%img%get_filtsz()         ! nr of resolution elements
        if( p%refine .eq. 'isw' )then
            wmat_states = merge_rmat_from_parts(p%nstates, p%nptcls, p%nparts, filtsz, 'shellweights_part')
            call normalise_shellweights(wmat_states, p%npeaks)
            filnum = get_fileunit()
            open(unit=filnum, status='REPLACE', action='WRITE', file='shellweights.bin', access='STREAM')
            write(unit=filnum,pos=1,iostat=io_stat) wmat_states
            deallocate(wmat_states)
        else
            wmat = merge_rmat_from_parts(p%nptcls, p%nparts, filtsz, 'shellweights_part')
            call normalise_shellweights(wmat)
            filnum = get_fileunit()
            open(unit=filnum, status='REPLACE', action='WRITE', file='shellweights.bin', access='STREAM')
            write(unit=filnum,pos=1,iostat=io_stat) wmat
            deallocate(wmat)
        endif
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to shellweights.bin'
            stop 'I/O error; merge_shellweights'
        endif
        close(filnum)
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_SHELLWEIGHTS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_merge_shellweights
    
    subroutine exec_merge_similarities( self, cline )
        use simple_map_reduce, only: merge_similarities_from_parts
        class(merge_similarities_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(params)      :: p
        real, allocatable :: simmat(:,:)
        integer           :: filnum, io_stat
        p      = params(cline) ! parameters generated
        simmat = merge_similarities_from_parts(p%nptcls, p%nparts)
        filnum = get_fileunit()
        open(unit=filnum, status='REPLACE', action='WRITE', file='smat.bin', access='STREAM')
        write(unit=filnum,pos=1,iostat=io_stat) simmat
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to smat.bin'
            stop 'I/O error; simple_merge_similarities'
        endif
        close(filnum)
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_SIMILARITIES NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_merge_similarities

    subroutine exec_split_pairs( self, cline )
        use simple_map_reduce, only: split_pairs_in_parts
        class(split_pairs_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params) :: p
        p = params(cline) ! parameters generated
        call split_pairs_in_parts(p%nptcls, p%nparts)
        call simple_end('**** SIMPLE_SPLIT_PAIRS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_split_pairs

    subroutine exec_split( self, cline )
        use simple_map_reduce ! use all in there
        use simple_image, only: image
        class(split_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(params)         :: p
        type(image)          :: img
        integer              :: iptcl, ipart, ldim(3), cnt, nimgs
        character(len=4)     :: ext
        integer, allocatable :: parts(:,:)
        logical              :: either_defined
        p = params(cline) ! parameters generated
        ext = '.'//fname2ext(p%stk)
        call find_ldim_nptcls(p%stk, ldim, nimgs)
        ldim(3) = 1
        call img%new(ldim,1.)
        select case(p%split_mode)
            case('chunk')
                parts = split_nobjs_in_chunks(nimgs, p%chunksz)
            case('even')
                parts = split_nobjs_even(nimgs, p%nparts)
            case DEFAULT
                write(*,*) 'split mode: ', trim(p%split_mode)
                stop 'unsupported split_mode; simple_commander_distr :: exec_split'
        end select
        if( size(parts,1) /= p%nparts ) stop 'ERROR! generated number of parts not same as inputted nparts'
        do ipart=1,p%nparts
            call progress(ipart,p%nparts)
            cnt = 0
            do iptcl=parts(ipart,1),parts(ipart,2)
                cnt = cnt+1
                call img%read(p%stk, iptcl)
                if( p%neg .eq. 'yes' ) call img%neg
                call img%write('stack_part'//int2str_pad(ipart,p%numlen)//ext, cnt)
            end do
        end do
        deallocate(parts)
        call img%kill
        call simple_end('**** SIMPLE_SPLIT NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_split

end module simple_commander_distr
