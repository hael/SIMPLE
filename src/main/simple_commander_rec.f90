! concrete commander: 3D reconstruction routines
module simple_commander_rec
use simple_defs
use simple_cmdline,         only: cmdline
use simple_params,          only: params
use simple_build,           only: build
use simple_commander_base,  only: commander_base
use simple_strings,         only: int2str_pad
use simple_syslib,          only: wait_for_closure
use simple_hadamard_common  ! use all in there
use simple_fileio           ! use all in there
use simple_jiffys           ! use all in there
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

    !> RECVOL is a SIMPLE program to reconstruct volumes from EM stacks and their estimated orientations
    !!
    !! \see http://simplecryoem.com/tutorials.html?#resolution-estimate-from-single-particle-images
    !!
    !! `map2ptcls' generates a document named mapped_params_ptcls.txt that we can directly use to reconstruct a new map from the particle images and calculate the resolution
    !!
    !! ```sh
    !! simple_distr_exec prg=recvol eo=yes stk=../stack/sumstack.mrc \
    !!    oritab=mapped_params_ptcls.txt smpd=1.62 msk=88 ctf=yes \
    !!    pgrp=d2 nparts=2 nthr=4 >& EOREC
    !!```
    !!    This will only takes a few minutes and will print out the FSC values
    !!    in the EOREC file. The eo=yes option specifies that the resolution
    !!    will be calculated from the FSC between the even/odd halves of the
    !!    dataset. The final reconstruction is the recvol_state01.mrc. The end
    !!    of the RESOLUTION file produced gives you the calculated resolution
    !!
    !!    >>> RESOLUTION AT FSC=0.143 DETERMINED TO:    9.87
    !!    >>> RESOLUTION AT FSC=0.500 DETERMINED TO:  13.38
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
            case( 'yes' )
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

    !> EO_VOLASSEMBLE is a SIMPLE program to reconstruct volume with EO enabled
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
        integer                       :: part, s, n, ss, state4name, file_stat, fnr
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
        call eorecvol_read%kill_exp        ! reduced memory usage
        n = p%nstates*p%nparts
        do ss=1,p%nstates
            if( cline%defined('state') )then
                s = 1
            else
                s = ss
            endif
            DebugPrint  'processing state: ', s
            if( b%a%get_pop( s, 'state' ) == 0 )cycle           ! Empty state
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
        
        call fopen(fnr, FILE='VOLASSEMBLE_FINISHED', STATUS='REPLACE', action='WRITE', iostat=file_stat)
        call fileio_errmsg('In: commander_rec :: eo_volassemble', file_stat )
        call fclose( fnr , errmsg='In: commander_rec :: eo_volassemble')
        call wait_for_closure('VOLASSEMBLE_FINISHED')
        
        contains

            subroutine assemble( fbody )
                character(len=*), intent(in) :: fbody
                call eorecvol_read%read_eos(trim(fbody))
                ! sum the Fourier coefficients
                call b%eorecvol%sum(eorecvol_read)
            end subroutine assemble
            
            subroutine normalize( recname )
                character(len=*), intent(in)  :: recname
                character(len=STDLEN) :: volname
                volname = trim(recname)//trim(p%ext)
                call b%eorecvol%sum_eos
                call b%eorecvol%sampl_dens_correct_eos(s)
                call b%eorecvol%get_res(res05s(s), res0143s(s))
                call b%eorecvol%sampl_dens_correct_sum( b%vol )
                call b%eorecvol%write_eos(recname)
                call b%vol%write( volname, del_if_exists=.true. )
                call wait_for_closure( volname )
            end subroutine normalize

    end subroutine exec_eo_volassemble

    !> VOLASSEMBLE is a SIMPLE program to reconstruct volumes
    subroutine exec_volassemble( self, cline )
        use simple_reconstructor, only: reconstructor
        class(volassemble_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)                  :: p
        type(build)                   :: b
        character(len=:), allocatable :: fbody
        character(len=STDLEN)         :: recvolname, rho_name
        integer                       :: part, s, ss, endit, i, state4name, file_stat, fnr
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
        endit = 1
        if( p%eo .eq. 'yes' ) endit = 2
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
                do i=1,endit
                    if( cline%defined('even') .or. cline%defined('odd') )then
                        if( p%even .eq. 'yes' .and. p%odd .eq. 'no' )then
                            p%vols(s) = fbody//'_even'//p%ext
                            rho_name  = 'rho_'//fbody//'_even'//p%ext
                        else if( p%odd .eq. 'yes' .and. p%even .eq. 'no' )then
                            p%vols(s) = fbody//'_odd'//p%ext
                            rho_name  = 'rho_'//fbody//'_odd'//p%ext
                        else if( p%odd .eq. 'yes' .and. p%even .eq. 'yes' )then
                            stop 'ERROR! Cannot have even=yes and odd=yes simultaneously'
                        endif
                    else
                        if( p%eo .eq. 'yes' )then
                            if( i == 1 )then
                                p%vols(s) = fbody//'_odd'//p%ext
                                rho_name  = 'rho_'//fbody//'_odd'//p%ext
                            else
                                p%vols(s) = fbody//'_even'//p%ext
                                rho_name  = 'rho_'//fbody//'_even'//p%ext
                            endif   
                        else
                            p%vols(s) = fbody//p%ext
                            rho_name  = 'rho_'//fbody//p%ext
                        endif
                    endif
                    call assemble(p%vols(s), trim(rho_name))
                end do
                deallocate(fbody)
            end do
            if( p%nstates==1 .and. cline%defined('outvol') )then
                recvolname = trim(p%outvol)
            else
                if( p%even .eq. 'yes' .and. p%odd .eq. 'no' )then
                    recvolname = 'recvol_state'//int2str_pad(state4name,2)//'_even'//p%ext
                else if( p%odd .eq. 'yes' .and. p%even .eq. 'no' )then
                    recvolname = 'recvol_state'//int2str_pad(state4name,2)//'_odd'//p%ext
                else if( p%odd .eq. 'yes' .and. p%even .eq. 'yes' )then
                    stop 'ERROR! Cannot have even=yes and odd=yes simultaneously'
                else
                    recvolname = 'recvol_state'//int2str_pad(state4name,2)//p%ext
                endif
            endif
            call normalize( trim(recvolname) )
        end do
        call recvol_read%dealloc_rho
        call recvol_read%kill
        ! end gracefully
        call simple_end('**** SIMPLE_VOLASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        call fopen(fnr, FILE='VOLASSEMBLE_FINISHED', STATUS='REPLACE', action='WRITE', iostat=file_stat)
        call fileio_errmsg('In: commander_rec :: volassemble', file_stat )
        call fclose( fnr ,errmsg='In: commander_rec :: volassemble')
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
                call b%recvol%sampl_dens_correct
                call b%recvol%bwd_ft
                call b%recvol%clip(b%vol)
                call b%vol%write(recname, del_if_exists=.true.)
                call wait_for_closure(recname)
            end subroutine normalize

    end subroutine exec_volassemble

end module simple_commander_rec
