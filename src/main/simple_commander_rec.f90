! concrete commander: 3D reconstruction routines
module simple_commander_rec
use simple_cmdline,         only: cmdline
use simple_params,          only: params
use simple_build,           only: build
use simple_commander_base,  only: commander_base
use simple_strings,         only: int2str_pad
use simple_syslib,          only: wait_for_closure
use simple_hadamard_common  ! use all in there
use simple_fileio           ! use all in there
use simple_jiffys           ! use all in there
use simple_defs             ! use all in there
implicit none

public :: recvol_commander
public :: eo_volassemble_commander
public :: volassemble_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: recvol_commander
  contains
    procedure :: execute      => exec_recvol
end type recvol_commander
type, extends(commander_base) :: eo_volassemble_commander
  contains
    procedure :: execute      => exec_eo_volassemble
end type eo_volassemble_commander
type, extends(commander_base) :: volassemble_commander
  contains
    procedure :: execute      => exec_volassemble
end type volassemble_commander

contains

    !> for reconstructing volumes from image stacks and their estimated orientations
    subroutine exec_recvol( self, cline )
        use simple_rec_master, only: exec_rec_master
        class(recvol_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer      :: fnr, file_stat
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        select case(p%eo)
            case( 'yes', 'aniso' )
                call b%build_eo_rec_tbox(p) ! eo_reconstruction objs built
            case( 'no' )
                call b%build_rec_tbox(p)    ! reconstruction objects built
            case DEFAULT
                stop 'unknonw eo flag; simple_commander_rec :: exec_recvol'
        end select
        call exec_rec_master(b, p, cline)
        ! end gracefully
        call simple_end('**** SIMPLE_RECVOL NORMAL STOP ****', print_simple=.false.) 
    end subroutine exec_recvol

    !> for assembling even/odd volumes generated with distributed execution
    subroutine exec_eo_volassemble( self, cline )
        use simple_eo_reconstructor, only: eo_reconstructor
        class(eo_volassemble_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(params)                  :: p
        type(build)                   :: b
        type(eo_reconstructor)        :: eorecvol_read
        character(len=:), allocatable :: fname
        real, allocatable             :: res05s(:), res0143s(:)
        real                          :: res
        integer                       :: part, s, alloc_stat, n, ss, state4name, file_stat, fnr
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        call b%build_eo_rec_tbox(p)         ! reconstruction toolbox built
        call b%eorecvol%kill_exp            ! reduced meory usage
        call b%mskvol%kill                  ! reduced memory usage
        allocate(res05s(p%nstates), res0143s(p%nstates), stat=alloc_stat)
        call alloc_errchk("In: simple_eo_volassemble", alloc_stat)
        res0143s = 0.
        res05s   = 0.
        ! rebuild b%vol according to box size (beacuse it is otherwise boxmatch)
        call b%vol%new([p%box,p%box,p%box], p%smpd)
        call eorecvol_read%new(p)
        call eorecvol_read%kill_exp ! reduced memory usage
        n = p%nstates*p%nparts
        do ss=1,p%nstates
            if( cline%defined('state') )then
                s = 1
            else
                s = ss
            endif
            DebugPrint  'processing state: ', s
            if( b%a%get_pop( s, 'state' ) == 0 )cycle ! Empty state
            call b%eorecvol%reset_all
            do part=1,p%nparts
                if( cline%defined('state') )then
                    state4name = p%state
                else
                    state4name = s
                endif
                allocate(fname, source='recvol_state'//int2str_pad(state4name,2)//'_part'//int2str_pad(part,p%numlen))
                DebugPrint  'processing file: ', fname
                call assemble(fname)
                deallocate(fname)
            end do
            call normalize('recvol_state'//int2str_pad(state4name,2))
        end do
        ! set the resolution limit according to the worst resolved model
        res  = maxval(res0143s)
        p%lp = max( p%lpstop,res )
        write(*,'(a,1x,F6.2)') '>>> LOW-PASS LIMIT:', p%lp
        call eorecvol_read%kill
        ! end gracefully
        call simple_end('**** SIMPLE_EO_VOLASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        if(.not.fopen(fnr, FILE='VOLASSEMBLE_FINISHED', STATUS='REPLACE', action='WRITE', iostat=file_stat))&
             call  fileio_errmsg('In: commander_rec :: eo_volassemble', file_stat )
        if(.not.fclose( fnr , iostat=file_stat))&
             call  fileio_errmsg('In: commander_rec :: eo_volassemble', file_stat )
        call wait_for_closure('VOLASSEMBLE_FINISHED')
        
        contains

            subroutine assemble( fbody )
                character(len=*), intent(in) :: fbody
                call eorecvol_read%read_eos(trim(fbody))
                ! sum the Fourier coefficients
                call b%eorecvol%sum(eorecvol_read)
            end subroutine assemble
            
            subroutine normalize( recname )
                use simple_filterer, only: gen_anisotropic_optlp
                character(len=*), intent(in) :: recname
                character(len=STDLEN) :: volname
                character(len=32)     :: eonames(2)
                volname = trim(recname)//trim(p%ext)
                call b%eorecvol%sum_eos
                if( p%eo .eq. 'aniso' )then
                    ! anisotropic resolution model
                    eonames(1) = 'tmpvoleven'//p%ext
                    eonames(2) = 'tmpvolodd'//p%ext
                    call b%eorecvol%sampl_dens_correct_eos(s, eonames)
                    call gen_projection_frcs( p, eonames(1), eonames(2), s, b%projfrcs)
                    call del_file(eonames(1))
                    call del_file(eonames(2))
                    call b%projfrcs%write('frcs_state'//int2str_pad(s,2)//'.bin')
                    ! generate the anisotropic 3D optimal low-pass filter
                    call gen_anisotropic_optlp(b%vol, b%projfrcs, b%e, s)
                    call b%vol%write('aniso_optlp_state'//int2str_pad(s,2)//p%ext)
                else
                    call b%eorecvol%sampl_dens_correct_eos(s)
                endif
                call b%eorecvol%get_res(res05s(s), res0143s(s))
                call b%eorecvol%sampl_dens_correct_sum( b%vol )
                call b%eorecvol%write_eos(recname)
                call b%vol%write( volname, del_if_exists=.true. )
                call wait_for_closure( volname )
            end subroutine normalize

    end subroutine exec_eo_volassemble

    !> for assembling a volume generated with distributed execution
    subroutine exec_volassemble( self, cline )
        use simple_reconstructor, only: reconstructor
        class(volassemble_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)                  :: p
        type(build)                   :: b
        character(len=:), allocatable :: fbody
        character(len=STDLEN)         :: recvolname, rho_name
        integer                       :: part, s, ss, i, state4name, file_stat, fnr ! endit
        type(reconstructor)           :: recvol_read
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        call b%build_rec_tbox(p)            ! reconstruction toolbox built
        ! rebuild b%vol according to box size (because it is otherwise boxmatch)
        call b%vol%new([p%box,p%box,p%box], p%smpd)
        if( cline%defined('find') )then
            p%lp = b%img%get_lp(p%find)
        endif
        call recvol_read%new([p%boxpd,p%boxpd,p%boxpd], p%smpd)
        call recvol_read%alloc_rho(p)
        do ss=1,p%nstates
            if( cline%defined('state') )then
                s = 1
            else
                s = ss
            endif
            DebugPrint  'processing state: ', s
            if( b%a%get_pop( s, 'state' ) == 0 ) cycle ! Empty state
            call b%recvol%reset
            do part=1,p%nparts
                if( cline%defined('state') )then
                    state4name = p%state
                else
                    state4name = s
                endif
                allocate(fbody, source='recvol_state'//int2str_pad(state4name,2)//'_part'//int2str_pad(part,p%numlen))
                DebugPrint  'processing fbody: ', fbody
                p%vols(s) = fbody//p%ext
                rho_name  = 'rho_'//fbody//p%ext
                call assemble(p%vols(s), trim(rho_name))
                deallocate(fbody)
            end do
            if( p%nstates == 1 .and. cline%defined('outvol') )then
                recvolname = trim(p%outvol)
            else
                recvolname = 'recvol_state'//int2str_pad(state4name,2)//p%ext
            endif
            call normalize( trim(recvolname) )
        end do
        call recvol_read%dealloc_rho
        call recvol_read%kill
        ! end gracefully
        call simple_end('**** SIMPLE_VOLASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        if(.not.fopen(fnr, FILE='VOLASSEMBLE_FINISHED', STATUS='REPLACE', action='WRITE', iostat=file_stat))&
             call  fileio_errmsg('In: commander_rec :: volassemble', file_stat )
        if(.not.fclose( fnr , iostat=file_stat))&
             call  fileio_errmsg('In: commander_rec :: volassemble', file_stat )
        call wait_for_closure('VOLASSEMBLE_FINISHED')

        contains

            subroutine assemble( recnam, kernam )
                character(len=*), intent(in) :: recnam
                character(len=*), intent(in) :: kernam
                logical                      :: here(2)
                here(1)=file_exists(recnam)
                here(2)=file_exists(kernam)
                if( all(here) )then     
                    call recvol_read%read(recnam)
                    call recvol_read%read_rho(kernam)
                    call b%recvol%sum(recvol_read)
                else
                    if( .not. here(1) ) write(*,'(A,A,A)') 'WARNING! ', adjustl(trim(recnam)), ' missing'
                    if( .not. here(2) ) write(*,'(A,A,A)') 'WARNING! ', adjustl(trim(kernam)), ' missing'
                    return
                endif
            end subroutine assemble
    
            subroutine normalize( recname )
                character(len=*), intent(in) :: recname
                call b%recvol%gridding_correct
                call b%recvol%bwd_ft
                call b%recvol%clip(b%vol)
                call b%vol%write(recname, del_if_exists=.true.)
                call wait_for_closure(recname)
            end subroutine normalize

    end subroutine exec_volassemble

end module simple_commander_rec
