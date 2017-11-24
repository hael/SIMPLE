! concrete commander: 3D reconstruction routines
module simple_commander_rec
#include "simple_lib.f08"
use simple_defs_fname
use simple_cmdline,         only: cmdline
use simple_params,          only: params
use simple_build,           only: build
use simple_commander_base,  only: commander_base
use simple_projection_frcs, only: projection_frcs
use simple_hadamard_common  ! use all in there
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
        character(len=:), allocatable :: fbody, finished_fname
        real, allocatable             :: res05s(:), res0143s(:)
        real                          :: res
        integer                       :: part, s, n, ss, state, ldim(3), find4eoavg
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        call b%build_eo_rec_tbox(p)         ! reconstruction toolbox built
        allocate(res05s(p%nstates), res0143s(p%nstates), stat=alloc_stat)
        allocchk("In: simple_eo_volassemble res05s res0143s")
        res0143s = 0.
        res05s   = 0.
        ! rebuild b%vol according to box size (beacuse it is otherwise boxmatch)
        call b%vol%new([p%box,p%box,p%box], p%smpd)
        call eorecvol_read%new(p)
        n = p%nstates*p%nparts
        do ss=1,p%nstates
            if( cline%defined('state') )then
                s     = 1        ! index in recvol
                state = p%state  ! actual state
            else
                s     = ss
                state = ss
            endif
            if( b%a%get_pop(state, 'state' ) == 0 )cycle ! Empty state
            call b%eorecvol%reset_all
            do part=1,p%nparts
                allocate(fbody, source='recvol_state'//int2str_pad(state,2)//'_part'//int2str_pad(part,p%numlen))
                DebugPrint  'processing file: ', fbody
                call assemble(fbody)
                deallocate(fbody)
            end do
            call correct_for_sampling_density_and_estimate_res('recvol_state'//int2str_pad(state,2))
            if( cline%defined('state') )exit
        end do
        ! set the resolution limit according to the worst resolved model
        res  = maxval(res0143s)
        p%lp = max( p%lpstop,res )
        write(*,'(a,1x,F6.2)') '>>> LOW-PASS LIMIT:', p%lp
        call eorecvol_read%kill
        ! end gracefully
        call simple_end('**** SIMPLE_EO_VOLASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        if( cline%defined('state') )then
            allocate( finished_fname, source='VOLASSEMBLE_FINISHED_STATE'//int2str_pad(state,2))
        else
            allocate( finished_fname, source='VOLASSEMBLE_FINISHED' )
        endif
        call simple_touch( finished_fname , errmsg='In: commander_rec::eo_volassemble')

        contains

            subroutine assemble( fbody )
                character(len=*), intent(in) :: fbody
                call eorecvol_read%read_eos_exp(trim(fbody))
                ! sum the Fourier coefficients
                call b%eorecvol%sum_exp(eorecvol_read)
            end subroutine assemble

            subroutine correct_for_sampling_density_and_estimate_res( recname )
                use simple_filterer, only: gen_anisotropic_optlp
                character(len=*), intent(in) :: recname
                character(len=STDLEN)        :: volname
                character(len=32)            :: eonames(2), resmskname
                volname    = trim(recname)//trim(p%ext)
                eonames(1) = trim(recname)//'_even'//trim(p%ext)
                eonames(2) = trim(recname)//'_odd'//trim(p%ext)
                resmskname = 'resmask'//p%ext
                call b%eorecvol%compress_exp
                call b%eorecvol%sum_eos
                call b%eorecvol%sampl_dens_correct_eos(state, eonames(1), eonames(2), resmskname, find4eoavg)
                call gen_projection_frcs( b, p, cline, eonames(1), eonames(2), resmskname, s, b%projfrcs)
                call b%projfrcs%write('frcs_state'//int2str_pad(state,2)//'.bin')
                call gen_anisotropic_optlp(b%vol2, b%projfrcs, b%e_bal, s, p%pgrp)
                call b%vol2%write('aniso_optlp_state'//int2str_pad(state,2)//p%ext)
                call b%eorecvol%get_res(res05s(s), res0143s(s))
                call b%eorecvol%sampl_dens_correct_sum( b%vol )
                call b%vol%write( volname, del_if_exists=.true. )
                call wait_for_closure( volname )
                ! need to put the sum back at lowres for the eo pairs 
                call b%vol%fwd_ft
                call b%vol2%zero_and_unflag_ft
                call b%vol2%read(eonames(1))
                call b%vol2%fwd_ft
                call b%vol2%insert_lowres(b%vol, find4eoavg)
                call b%vol2%bwd_ft
                call b%vol2%write(eonames(1), del_if_exists=.true.)
                call b%vol2%zero_and_unflag_ft
                call b%vol2%read(eonames(2))
                call b%vol2%fwd_ft
                call b%vol2%insert_lowres(b%vol, find4eoavg)
                call b%vol2%bwd_ft
                call b%vol2%write(eonames(2), del_if_exists=.true.)
            end subroutine correct_for_sampling_density_and_estimate_res

    end subroutine exec_eo_volassemble

    !> for assembling a volume generated with distributed execution
    subroutine exec_volassemble( self, cline )
        use simple_reconstructor, only: reconstructor
        class(volassemble_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)                  :: p
        type(build)                   :: b
        character(len=:), allocatable :: fbody, finished_fname
        character(len=STDLEN)         :: recvolname, rho_name
        integer                       :: part, s, ss, state, ldim(3)
        type(reconstructor)           :: recvol_read
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        call b%build_rec_tbox(p)            ! reconstruction toolbox built
        ! rebuild b%vol according to box size (because it is otherwise boxmatch)
        call b%vol%new([p%box,p%box,p%box], p%smpd)
        call recvol_read%new([p%boxpd,p%boxpd,p%boxpd], p%smpd)
        call recvol_read%alloc_rho(p)
        do ss=1,p%nstates
            if( cline%defined('state') )then
                s     = 1        ! index in recvol
                state = p%state  ! actual state
            else
                s     = ss
                state = ss
            endif
            if( b%a%get_pop(state, 'state' ) == 0 ) cycle ! Empty state
            call b%recvol%reset
            call b%recvol%reset_exp
            do part=1,p%nparts
                allocate(fbody, source='recvol_state'//int2str_pad(state,2)//'_part'//int2str_pad(part,p%numlen))
                call assemble(fbody)
                deallocate(fbody)
            end do
            if( p%nstates == 1 .and. cline%defined('outvol') )then
                recvolname = trim(p%outvol)
            else
                recvolname = 'recvol_state'//int2str_pad(state,2)//p%ext
            endif
            call correct_for_sampling_density( trim(recvolname) )
            if( cline%defined('state') )exit
        end do
        call recvol_read%dealloc_rho
        call recvol_read%kill
        ! end gracefully
        call simple_end('**** SIMPLE_VOLASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        if( cline%defined('state') )then
            allocate( finished_fname, source='VOLASSEMBLE_FINISHED_STATE'//int2str_pad(state,2) )
        else
            allocate( finished_fname, source='VOLASSEMBLE_FINISHED' )
        endif
        call simple_touch( finished_fname, errmsg='In: commander_rec :: volassemble')

        contains

            subroutine assemble( fbody )
                character(len=*), intent(in)  :: fbody
                character(len=:), allocatable :: recname, rhoname
                logical :: here(2)
                allocate(recname, source=trim(adjustl(fbody))//BIN_EXT)
                allocate(rhoname, source='rho_'//trim(adjustl(fbody))//BIN_EXT)
                here(1)=file_exists(recname)
                here(2)=file_exists(rhoname)
                if( all(here) )then
                    call recvol_read%read_cmat_exp(recname)
                    call recvol_read%read_rho_exp(rhoname)
                    call b%recvol%sum_exp(recvol_read)
                else
                    if( .not. here(1) ) write(*,'(A,A,A)') 'WARNING! ', recname, ' missing'
                    if( .not. here(2) ) write(*,'(A,A,A)') 'WARNING! ', rhoname, ' missing'
                    return
                endif
            end subroutine assemble

            subroutine correct_for_sampling_density( recname )
                character(len=*), intent(in) :: recname
                call b%recvol%compress_exp
                call b%recvol%sampl_dens_correct
                call b%recvol%bwd_ft
                call b%recvol%clip(b%vol)
                call b%vol%write(recname, del_if_exists=.true.)
                call wait_for_closure(recname)
            end subroutine correct_for_sampling_density

    end subroutine exec_volassemble

end module simple_commander_rec
