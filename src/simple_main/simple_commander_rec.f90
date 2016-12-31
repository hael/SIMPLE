!==Class simple_commander_rec
!
! This class contains the set of concrete 3D reconstruction commanders of the SIMPLE library. This class provides the glue between the reciver 
! (main reciever is simple_exec program) and the abstract action, which is simply execute (defined by the base class: simple_commander_base). 
! Later we can use the composite pattern to create MacroCommanders (or workflows)
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Cyril Reboul & Hans Elmlund 2016
!
module simple_commander_rec
use simple_defs             ! singleton
use simple_jiffys           ! singleton
use simple_hadamard_common  ! singleton
use simple_cmdline,         only: cmdline
use simple_params,          only: params
use simple_build,           only: build
use simple_commander_base,  only: commander_base
implicit none

public :: recvol_commander
public :: eo_volassemble_commander
public :: volassemble_commander
private

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

    subroutine exec_recvol( self, cline )
        use simple_rec_master, only: exec_rec_master
        class(recvol_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params)      :: p
        type(build)       :: b
        logical           :: doshellweight
        real, allocatable :: wmat(:,:)
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
        call setup_shellweights(b, p, doshellweight, wmat) ! shell-weights setup
        if( doshellweight )then
            call exec_rec_master(b, p, cline, wmat=wmat)
        else
            call exec_rec_master(b, p, cline)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_RECVOL NORMAL STOP ****')    
    end subroutine exec_recvol

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
        integer                       :: part, s, alloc_stat, cnt, n, ss, state4name
        logical                       :: debug=.false.
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        if( cline%defined('nstates') )then
            if( p%nstates /= b%a%get_nstates() ) stop 'Inconsistent number of states between command-line and oritab'
        endif
        call b%build_eo_rec_tbox(p)   ! reconstruction toolbox built
        allocate(res05s(p%nstates), res0143s(p%nstates), stat=alloc_stat)
        call alloc_err("In: simple_eo_volassemble", alloc_stat)
        ! rebuild b%vol according to box size (beacuse it is otherwise boxmatch)
        call b%vol%new([p%box,p%box,p%box], p%smpd, p%imgkind)
        call eorecvol_read%new(p)
        n = p%nstates*p%nparts
        cnt = 0
        do ss=1,p%nstates
            if( cline%defined('state') )then
                s = 1
            else
                s = ss
            endif
            if( debug ) write(*,*) 'processing state: ', s
            call b%eorecvol%reset_all
            do part=1,p%nparts
                cnt = cnt+1
                call progress(cnt,n)
                if( cline%defined('state') )then
                    state4name = p%state
                else
                    state4name = s
                endif
                allocate(fname, source='recvol'//'_state'//int2str_pad(state4name,2)//'_part'//int2str_pad(part,p%numlen))
                if( debug ) write(*,*) 'processing file: ', fname
                call assemble(fname)
                deallocate(fname)
            end do
            call normalize('recvol_state'//int2str_pad(state4name,2))
        end do
        ! set the resolution limit according to the worst resolved model
        res  = maxval(res0143s)
        p%lp = max( p%lpstop,res )
        write(*,'(a,1x,F6.2)') '>>> LOW-PASS LIMIT:', p%lp
        write(0,'(a)') "GENERATED VOLUMES: recvol*.ext"
        ! end gracefully
        call simple_end('**** SIMPLE_EO_VOLASSEMBLE NORMAL STOP ****')
        
        contains

            subroutine assemble( fbody )
                character(len=*), intent(in) :: fbody
                call eorecvol_read%read_eos(fbody)
                ! sum the Fourier coefficients
                call b%eorecvol%sum(eorecvol_read)
            end subroutine
            
            subroutine normalize( recnam )
                use simple_image, only: image
                character(len=*), intent(in)  :: recnam
                call b%eorecvol%sum_eos
                call b%eorecvol%sampl_dens_correct_eos(s)
                call b%eorecvol%get_res(res05s(s), res0143s(s))
                call b%eorecvol%sampl_dens_correct_sum(b%vol)
                call b%eorecvol%write_eos(recnam)
                call b%vol%write(recnam//p%ext, del_if_exists=.true.)
            end subroutine

    end subroutine exec_eo_volassemble
    
    subroutine exec_volassemble( self, cline )
        use simple_reconstructor, only: reconstructor
        class(volassemble_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)                  :: p
        type(build)                   :: b
        character(len=:), allocatable :: fbody
        character(len=STDLEN)         :: recvolname
        integer                       :: part, s, ss, endit, i, cnt, state4name
        type(reconstructor)           :: recvol_read
        logical                       :: here(2)
        logical, parameter            :: debug=.false.
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        call b%build_rec_tbox(p)            ! reconstruction toolbox built
        ! rebuild b%vol according to box size (because it is otherwise boxmatch)
        call b%vol%new([p%box,p%box,p%box], p%smpd, p%imgkind)
        if( cline%defined('find') )then
            p%lp = b%img%get_lp(p%find)
        endif
        call recvol_read%new([p%boxpd,p%boxpd,p%boxpd], p%smpd, p%imgkind)
        call recvol_read%alloc_rho(p)
        endit = 1
        if( p%eo .eq. 'yes' ) endit = 2
        cnt = 0
        do ss=1,p%nstates
            if( cline%defined('state') )then
                s = 1
            else
                s = ss
            endif
            if( debug ) write(*,*) 'processing state: ', s
            call b%recvol%reset
            do part=1,p%nparts
                cnt = cnt+1
                call progress(cnt,p%nstates*p%nparts)
                if( cline%defined('state') )then
                    state4name = p%state
                else
                    state4name = s
                endif
                allocate(fbody, source='recvol'//'_state'//int2str_pad(state4name,2)//'_part'//int2str_pad(part,p%numlen))
                if( debug ) write(*,*) 'processing fbody: ', fbody
                do i=1,endit
                    if( cline%defined('even') .or. cline%defined('odd') )then
                        if( p%even .eq. 'yes' .and. p%odd .eq. 'no' )then
                            p%vols(s) = fbody//'_even'//p%ext
                            p%masks(s) = 'rho_'//fbody//'_even'//p%ext
                        else if( p%odd .eq. 'yes' .and. p%even .eq. 'no' )then
                            p%vols(s) = fbody//'_odd'//p%ext
                            p%masks(s) = 'rho_'//fbody//'_odd'//p%ext
                        else if( p%odd .eq. 'yes' .and. p%even .eq. 'yes' )then
                            stop 'ERROR! Cannot have even=yes and odd=yes simultaneously'
                        endif
                    else
                        if( p%eo .eq. 'yes' )then
                            if( i == 1 )then
                                p%vols(s) = fbody//'_odd'//p%ext
                                p%masks(s) = 'rho_'//fbody//'_odd'//p%ext
                            else
                                p%vols(s) = fbody//'_even'//p%ext
                                p%masks(s) = 'rho_'//fbody//'_even'//p%ext
                            endif   
                        else
                            p%vols(s)  = fbody//p%ext
                            p%masks(s) = 'rho_'//fbody//p%ext
                        endif
                    endif
                    call assemble(p%vols(s), p%masks(s))
                end do
                deallocate(fbody)
            end do
            if( p%nstates==1 .and. cline%defined('outvol') )then
                recvolname = trim(p%outvol)
            else
                if( p%even .eq. 'yes' .and. p%odd .eq. 'no' )then
                    recvolname = 'recvol_state'//int2str_pad(s,2)//'_even'//p%ext
                else if( p%odd .eq. 'yes' .and. p%even .eq. 'no' )then
                    recvolname = 'recvol_state'//int2str_pad(s,2)//'_odd'//p%ext
                else if( p%odd .eq. 'yes' .and. p%even .eq. 'yes' )then
                    stop 'ERROR! Cannot have even=yes and odd=yes simultaneously'
                else
                    recvolname = 'recvol_state'//int2str_pad(s,2)//p%ext
                endif
            endif
            call normalize( trim(recvolname) )
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_VOLASSEMBLE NORMAL STOP ****')

        contains

            subroutine assemble( recnam, kernam )
                character(len=*), intent(in) :: recnam
                character(len=*), intent(in) :: kernam
                inquire(FILE=recnam, EXIST=here(1))
                inquire(FILE=kernam, EXIST=here(2))
                if( all(here) )then     
                    call recvol_read%read(recnam)
                    if( debug )then
                        if( recvol_read%contains_nans() )then
                            write(*,*) 'WARNING! recvol: ', recnam, 'contains NaN:s'
                        endif
                    endif
                    call recvol_read%read_rho(kernam)
                    if( debug )then
                        if( recvol_read%rho_contains_nans() )then
                            write(*,*) 'WARNING! kernel: ', kernam, 'contains NaN:s'
                        endif
                    endif
                    call b%recvol%sum(recvol_read)
                    if( debug )then
                        if( b%recvol%contains_nans() )then
                            write(*,*) 'WARRNING! summed image part contains NaN:s'
                        endif 
                        if( b%recvol%rho_contains_nans() )then
                            write(*,*) 'WARNING! summed kernel part contains NaN:s'
                        endif
                    endif
                else
                    if( .not. here(1) ) write(*,'(A,A,A)') 'WARNING! ', adjustl(trim(recnam)), ' missing'
                    if( .not. here(2) ) write(*,'(A,A,A)') 'WARNING! ', adjustl(trim(kernam)), ' missing'
                    return
                endif
            end subroutine
    
            subroutine normalize( recnam )
                character(len=*), intent(in) :: recnam
                call b%recvol%sampl_dens_correct
                call b%recvol%bwd_ft
                call b%recvol%clip(b%vol)
                call b%vol%write(recnam, del_if_exists=.true.)
            end subroutine
    end subroutine exec_volassemble

end module simple_commander_rec
