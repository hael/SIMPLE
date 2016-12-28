!==Class simple_commander_prime2D
!
! This class contains the set of concrete prime2D commanders of the SIMPLE library. This class provides the glue between the reciver 
! (main reciever is simple_exec program) and the abstract action, which is simply execute (defined by the base class: simple_commander_base). 
! Later we can use the composite pattern to create MacroCommanders (or workflows)
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Cyril Reboul & Hans Elmlund 2016
!
module simple_commander_prime2D
use simple_defs            ! singleton
use simple_jiffys          ! singleton
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
implicit none

public :: prime2D_init_commander
public :: prime2D_commander
public :: cavgassemble_commander
public :: check2D_conv_commander
public :: rank_cavgs_commander
private

type, extends(commander_base) :: prime2D_init_commander 
  contains
    procedure :: execute      => exec_prime2D_init
end type prime2D_init_commander 
type, extends(commander_base) :: prime2D_commander 
  contains
    procedure :: execute      => exec_prime2D
end type prime2D_commander
type, extends(commander_base) :: cavgassemble_commander
  contains
    procedure :: execute      => exec_cavgassemble
end type cavgassemble_commander
type, extends(commander_base) :: check2D_conv_commander
  contains
    procedure :: execute      => exec_check2D_conv
end type check2D_conv_commander
type, extends(commander_base) :: rank_cavgs_commander
 contains
   procedure :: execute      => exec_rank_cavgs
end type rank_cavgs_commander

contains

    subroutine exec_prime2D_init( self, cline )
        use simple_hadamard2D_matcher, only: prime2D_assemble_sums, prime2D_write_sums, &
        & prime2D_write_partial_sums
        class(prime2D_init_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)  :: p
        type(build)   :: b
        integer       :: ncls_in_oritab, icls, fnr, file_stat
        p = params(cline)                     ! parameters generated
        p%boxmatch = p%box                    !!!!!!!!!!!!!!!!!! 4 NOW
        call b%build_general_tbox(p, cline)   ! general objects built
        call b%build_hadamard_prime2D_tbox(p) ! 2D Hadamard matcher built
        write(*,'(a)') '>>> GENERATING INITIAL CLUSTER CENTERS'
        if( cline%defined('oritab') )then
            call b%a%remap_classes
            ncls_in_oritab = b%a%get_ncls()
            if( p%ncls < ncls_in_oritab ) stop 'Inputted ncls < ncls_in_oritab; not allowed!'
            if( p%ncls > ncls_in_oritab )then
                call b%a%expand_classes(p%ncls)
            endif
        else
            if( p%srch_inpl .eq. 'yes' )then
                call b%a%rnd_cls(p%ncls)
            else
                call b%a%rnd_cls(p%ncls, srch_inpl=.false.)
            endif
        endif
        p%oritab = 'prime2D_startdoc.txt'
        if( p%mul > 1. ) call b%a%mul_shifts(p%mul)
        if( p%l_distr_exec .and. nint(cline%get_rarg('part')).ne.1 )then
        else
            call b%a%write(p%oritab)
        endif
        if( cline%defined('filwidth') )then
            if( p%l_distr_exec)then
                stop 'filwidth mode not implemented for distributed mode; simple_commander_prime2D.f90; exec_prime2D_init'
            endif
            do icls=1,p%ncls
                call b%cavgs(icls)%bin_filament(p%filwidth)
            end do
            call prime2D_write_sums(b, p)
        else
            call prime2D_assemble_sums(b, p)
            if( p%l_distr_exec)then
                call prime2D_write_partial_sums( b, p )
                fnr = get_fileunit()
                open(unit=fnr, FILE='JOB_FINISHED_'//int2str_pad(p%part,p%numlen),&
                STATUS='REPLACE', action='WRITE', iostat=file_stat)
                call fopen_err( 'In: simple_commander_prime2D :: exec_prime2D_init', file_stat )
                close(fnr)
            else
                call prime2D_write_sums(b, p)
            endif
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_PRIME2D_INIT NORMAL STOP ****')
    end subroutine exec_prime2D_init
    
    subroutine exec_prime2D( self, cline )
        use simple_hadamard2D_matcher, only: prime2D_exec
        class(prime2D_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer      :: i, startit, ncls_from_refs, lfoo(3)
        logical      :: converged=.false.
        p = params(cline)                     ! parameters generated
        p%boxmatch = p%box                    !!!!!!!!!!!!!!!!!! 4 NOW
        call b%build_general_tbox(p, cline)   ! general objects built
        call b%build_hadamard_prime2D_tbox(p) ! 2D Hadamard matcher built
        if( p%srch_inpl .eq. 'no' )then
            if( .not. cline%defined('oritab') )then
                stop 'need oritab for this mode (srch_inpl=no) of execution!'
            endif
        endif
        if( cline%defined('refs') )then
            call find_ldim_nptcls(p%refs, lfoo, ncls_from_refs)
            ! consistency check
            if( p%ncls /=  ncls_from_refs ) stop 'nrefs /= inputted ncls'
        endif
        ! execute
        if( cline%defined('part') )then
            if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
            call prime2D_exec(b, p, cline, 0, converged) ! partition or not, depending on 'part'       
        else
            startit = 1
            if( cline%defined('startit') ) startit = p%startit
            do i=startit,p%maxits
                call prime2D_exec(b, p, cline, i, converged)
                if(converged) exit
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_PRIME2D NORMAL STOP ****')
    end subroutine exec_prime2D
    
    subroutine exec_cavgassemble( self, cline )
        use simple_hadamard2D_matcher, only: prime2D_assemble_sums_from_parts, prime2D_write_sums
        class(cavgassemble_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)                                 :: p
        type(build)                                  :: b
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        call b%build_hadamard_prime2D_tbox(p)
        call prime2D_assemble_sums_from_parts(b, p)
        if( cline%defined('which_iter') )then
            call prime2D_write_sums( b, p, p%which_iter)
        else
            call prime2D_write_sums( b, p )
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CAVGASSEMBLE NORMAL STOP ****')
    end subroutine exec_cavgassemble
    
    subroutine exec_check2D_conv( self, cline )
        class(check2D_conv_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        logical      :: converged
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        p%ncls    = b%a%get_ncls()
        converged = b%conv%check_conv2D()   ! convergence check
        call cline%set('frac', b%conv%get('frac'))
        if( p%doshift )then
            ! activates shift serach
            call cline%set('trs', p%trs)
        endif
        if( converged )then
            call cline%set('converged', 'yes')
        else
            call cline%set('converged', 'no')
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK2D_CONV STOP ****')
    end subroutine exec_check2D_conv
    
    subroutine exec_rank_cavgs( self, cline )
        use simple_oris, only: oris
        class(rank_cavgs_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params)         :: p
        type(build)          :: b
        integer              :: iclass
        integer, allocatable :: order(:)
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        p%ncls   = p%nptcls
        p%nptcls = nlines(p%oritab) 
        call b%a%new(p%nptcls)
        call b%a%read(p%oritab)
        order = b%a%order_cls()
        do iclass=1,p%ncls
            write(*,'(a,1x,i5,1x,a,i5)') 'CLASS:', order(iclass), 'POP:', b%a%get_clspop(order(iclass)) 
            call b%img%read(p%stk, order(iclass))
            call b%img%write(p%outstk, iclass)
        end do
        call simple_end('**** SIMPLE_RANK_CAVGS NORMAL STOP ****')
    end subroutine exec_rank_cavgs

end module simple_commander_prime2D
