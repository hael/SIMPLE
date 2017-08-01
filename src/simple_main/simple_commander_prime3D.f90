!------------------------------------------------------------------------------!
! SIMPLE v3.0         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> simple_commander_prime3D
!
!! This class contains the set of concrete prime3D commanders of the SIMPLE
!! library. This class provides the glue between the reciver (main reciever is
!! simple_exec program) and the abstract action, which is simply execute
!! (defined by the base class: simple_commander_base). Later we can use the
!! composite pattern to create MacroCommanders (or workflows)
!!
!!
!! @see  http://simplecryoem.com/tutorials.html?#ab-initio-3d-reconstruction-from-class-averages-using-prime3d
!! Ab initio 3D Reconstruction from Class Averages Using PRIME3D
!!
!! A major obstacle to achieving near-atomic resolution with single-particle
!! cryo-EM is the problem of generating an accurate de novo 3D reconstruction.
!! Many cryo-EM structures are therefore solved by alignment of the images to a
!! priori models. The use of prior information, in the form of either a starting
!! model from an independent source or an assumption of a particular point-group
!! symmetry, is associated with the risk of introducing model bias. The model
!! bias phenomenon is often illustrated by alignment of pure noise images to an
!! image of Einstein. The image of Einstein is almost perfectly reproduced when
!! the aligned noise images are averaged. It is often stated that low-pass
!! filtering of X-ray maps, before they are used as starting models, eliminates
!! model bias. This is a misunderstanding, as any model can be convincingly
!! reproduced from noisy images . Most refinement software, such as FREALIGN ,
!! RELION , or projection matching , depends on an accurate starting model for
!! convergence to a high-resolution map. If the starting model is not supported
!! by the images, there is a potent risk of introducing model bias. To what
!! degree a starting model can bias the final 3D structure needs to be better
!! characterised by methodological studies. We introduced the PRIME3D algorithm
!! to remove the requirement for a priori structural knowledge and open the
!! method to the study of particles with novel structure. Robust algorithms for
!! ab initio 3D reconstruction are particularly important for the analysis of
!! small particles with low symmetry.
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Cyril Reboul & Hans Elmlund 2016
!
module simple_commander_prime3D
use simple_defs
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_qsys_funs       ! use all in there
use simple_filehandling    ! use all in there
use simple_jiffys          ! use all in there
implicit none

public :: resrange_commander
public :: npeaks_commander
public :: nspace_commander
public :: prime3D_init_commander
public :: multiptcl_init_commander
public :: prime3D_commander
public :: cont3D_commander
public :: check3D_conv_commander
private

type, extends(commander_base) :: resrange_commander
  contains
    procedure :: execute      => exec_resrange
end type resrange_commander
type, extends(commander_base) :: npeaks_commander
  contains
    procedure :: execute      => exec_npeaks
end type npeaks_commander
type, extends(commander_base) :: nspace_commander
 contains
   procedure :: execute      => exec_nspace
end type nspace_commander
type, extends(commander_base) :: prime3D_init_commander
  contains
    procedure :: execute      => exec_prime3D_init
end type prime3D_init_commander
type, extends(commander_base) :: multiptcl_init_commander
  contains
    procedure :: execute      => exec_multiptcl_init
end type multiptcl_init_commander
type, extends(commander_base) :: prime3D_commander
  contains
    procedure :: execute      => exec_prime3D
end type prime3D_commander
type, extends(commander_base) :: cont3D_commander
  contains
    procedure :: execute      => exec_cont3D
end type cont3D_commander
type, extends(commander_base) :: check3D_conv_commander
  contains
    procedure :: execute      => exec_check3D_conv
end type check3D_conv_commander

contains

    subroutine exec_resrange( self, cline )
        use simple_hadamard3D_matcher, only: prime3D_find_resrange
        class(resrange_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        if( cline%defined('box') .or. cline%defined('moldiam') )then
            p = params(cline)                     ! parameters generated
            call b%build_general_tbox(p, cline)   ! general objects built
            call b%build_hadamard_prime3D_tbox(p) ! prime objects built
            call prime3D_find_resrange( b, p, p%lp, p%lpstop )
            p%lpstart = p%lp
            call cline%set('lpstart',p%lpstart)   ! for reporting
            call cline%set('lpstop',p%lpstop)     ! for reporting
            write(*,'(A,1X,F5.1)') '>>> LP START:', p%lpstart
            write(*,'(A,2X,F5.1)') '>>> LP STOP:',  p%lpstop
            write(*,'(A,2X,F5.1)') '>>> HP:',       p%hp
        else
            stop 'need either box size or moldiam to estimate resrange; simple_resrange'
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_RESRANGE NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_resrange

    subroutine exec_npeaks( self, cline )
        use simple_oris, only: oris
        class(npeaks_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(build)  :: b
        type(params) :: p
        p = params(cline) ! parameters generated
        call b%build_general_tbox(p, cline)
        p%npeaks = min(10,b%e%find_npeaks(p%lp, p%moldiam))
        write(*,'(A,1X,I4)') '>>> NPEAKS:', p%npeaks
        ! end gracefully
        call simple_end('**** SIMPLE_NPEAKS NORMAL STOP ****')
    end subroutine exec_npeaks

    subroutine exec_nspace(self,cline)
        use simple_math, only: resang
        use simple_oris, only: oris
        class(nspace_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(oris)   :: o
        type(params) :: p
        real         :: ares
        integer      :: i
        p = params(cline) ! parameters generated
        do i=500,5000,500
            o = oris(i)
            call o%spiral
            ares = o%find_angres()
            write(*,'(A,1X,I7,1X,A,1X,F5.2)') 'NR OF PROJDIRS:', i, 'RESOLUTION:', resang(ares, p%moldiam)
        end do
        call simple_end('**** SIMPLE_NSPACE NORMAL STOP ****')
    end subroutine exec_nspace

    subroutine exec_prime3D_init( self, cline )
        use simple_hadamard3D_matcher, only: gen_random_model, prime3D_find_resrange
        class(prime3D_init_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)       :: p
        type(build)        :: b
        integer, parameter :: MAXIMGS=1000
        p = params(cline) ! parameters generated
        if( p%l_xfel )then
            if( cline%defined('msk') .or. cline%defined('mw') .or.&
            cline%defined('nvox') .or. cline%defined('mskfile') )then
                stop 'no mask allowed when processing XFEL patterns; simple_prime3D_init'
            endif
        endif
        if( p%ctf .ne. 'no')then
            if( .not. cline%defined('deftab') )&
            &stop 'need texfile with defocus/astigmatism values for ctf .ne. no mode exec'
        endif
        call b%build_general_tbox(p, cline)   ! general objects built
        call b%build_hadamard_prime3D_tbox(p) ! prime3D objects built
        ! determine resolution range
        if( cline%defined('lp') ) call prime3D_find_resrange( b, p, p%lp, p%lpstop )
        ! determine the number of peaks
        if( .not. cline%defined('npeaks') ) p%npeaks = min(10,b%e%find_npeaks(p%lp, p%moldiam))
        ! generate the random model
        if( cline%defined('nran') )then
            call gen_random_model( b, p, p%nran )
        else
            if( p%nptcls > MAXIMGS )then
                call gen_random_model( b, p, MAXIMGS )
            else
                call gen_random_model( b, p )
            endif
        endif
        ! end gracefully
        call qsys_job_finished( p, cline%get_carg('prg') )
        call simple_end('**** SIMPLE_PRIME3D_INIT NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prime3D_init

    subroutine exec_multiptcl_init( self, cline )
        use simple_rec_master, only: exec_rec_master
        class(multiptcl_init_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        p = params(cline) ! constants & derived constants produced
        call b%build_general_tbox(p, cline)
        if( .not. cline%defined('oritab') )then
            call b%a%rnd_oris
            call b%a%zero_shifts
        endif
        if( cline%defined('state2split') )then
            if( cline%defined('oritab') )then
                p%nstates = b%a%get_nstates()
                call b%a%split_state(p%state2split)
                p%nstates = p%nstates + 1
            else
                stop 'Need oritab to be defined when state2split is defined on command line; simple_multiptcl_init'
            endif
        else if( p%tseries .eq. 'yes' )then
            call b%a%ini_tseries(p%nstates, 'state')
        else
            call b%a%rnd_states(p%nstates)
            if( p%nstates < 2 ) stop 'Nonsensical to have nstates < 2; simple_multiptcl_init'
        endif
        if( p%norec .eq. 'no' )then
            if( cline%defined('lp') )then
                call b%build_rec_tbox(p)
                p%eo = 'no'
            else
                call b%build_eo_rec_tbox(p)
                p%eo = 'yes'
            endif
            call exec_rec_master(b, p, cline, 'startvol')
        endif
        if( p%zero .eq. 'yes' ) call b%a%set_all2single('corr', 0.)
        call b%a%write('multiptcl_startdoc.txt')
        ! end gracefully
        call simple_end('**** SIMPLE_MULTIPTCL_INIT NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_multiptcl_init

    !> PRIME3D is a SIMPLE program

    subroutine exec_prime3D( self, cline )
        use simple_math, only: calc_lowpass_lim, calc_fourier_index
        use simple_hadamard3D_matcher, only: prime3D_exec, prime3D_find_resrange
        use simple_strings,            only: str_has_substr
        class(prime3D_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)      :: p
        type(build)       :: b
        integer           :: i, startit
        logical           :: update_res, converged
        real              :: lpstart, lpstop, corr, corr_prev
        update_res = .false.
        converged  = .false.
        p = params(cline) ! parameters generated
        if( p%l_xfel )then
            if( cline%defined('msk') .or. cline%defined('mw') .or.&
            cline%defined('nvox') .or. cline%defined('automsk') )then
                stop 'no mask allowed when processing XFEL patterns; simple_prime3D'
            endif
        endif
        if( str_has_substr(p%refine,'neigh') .or. p%refine .eq. 'shift' )then
            if( .not. cline%defined('oritab') )then
                stop 'need oritab input for execution of prime3D with this refine mode'
            endif
        endif
        call b%build_general_tbox(p, cline)   ! general objects built
        if( .not. cline%defined('eo') ) p%eo = 'no' ! default
        if( p%eo .eq. 'yes' ) p%dynlp = 'no'
        if( cline%defined('lp') .or. cline%defined('find')&
        .or. p%eo .eq. 'yes' .or. p%dynlp .eq. 'yes' )then
            ! alles ok!
        else
           stop 'need a starting low-pass limit (set lp or find)!'
        endif
        call b%build_hadamard_prime3D_tbox(p) ! prime3D objects built
        startit = 1
        if( cline%defined('startit') )startit = p%startit
        if( cline%defined('part') )then
            if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
            if( cline%defined('find') )then
                p%lp = calc_lowpass_lim( p%find, p%boxmatch, p%smpd )
            endif
            call prime3D_exec(b, p, cline, startit, update_res, converged) ! partition or not, depending on 'part'
        else
            if( p%dynlp .eq. 'yes' )then
                call prime3D_find_resrange( b, p, lpstart, lpstop ) ! determine resolution range
                if( cline%defined('lpstart') )then
                    p%lp = p%lpstart
                else
                    p%lp = lpstart
                endif
                if( cline%defined('lpstop') )then
                    ! defined aleady
                else
                    p%lpstop = lpstop
                endif
            endif
            p%find = calc_fourier_index( p%lp, p%boxmatch, p%smpd )
            ! extremal dynamics
            if( cline%defined('extr_thresh') )then
                ! all is well
            else
                ! starts from the top
                p%extr_thresh = EXTRINITHRESH/p%rrate
                if( startit > 1 )then
                    ! need to update the randomization rate
                    do i=1,startit-1
                         p%extr_thresh = p%extr_thresh * p%rrate
                    end do
                endif
            endif
            corr = -1
            do i=startit,p%maxits
                p%extr_thresh = p%extr_thresh * p%rrate
                call prime3D_exec(b, p, cline, i, update_res, converged)
                if( .not. p%l_distr_exec .and. p%refine .eq. 'snhc' .and. .not. cline%defined('szsn') )then
                    ! update stochastic neighborhood size if corr is not improving
                    corr_prev = corr
                    corr      = b%a%get_avg('corr')
                    if( i > 1 .and. corr <= corr_prev )then
                        p%szsn = min(SZSN_MAX,p%szsn + SZSN_STEP)
                    endif
                endif
                if( update_res )then
                    ! dynamic low-pass
                    p%find = p%find+p%fstep
                    p%lp   = max(p%lpstop, calc_lowpass_lim( p%find, p%boxmatch, p%smpd ))
                endif
                if( converged )exit
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_PRIME3D NORMAL STOP ****')
    end subroutine exec_prime3D
    !> CONT3D is a SIMPLE program
    subroutine exec_cont3D( self, cline )
        use simple_cont3D_matcher,         only: cont3D_exec
        class(cont3D_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer      :: i, startit
        logical      :: converged
        converged = .false.
        p = params(cline) ! parameters generated
        select case(p%refine)
            case('yes','de','ada')
                ! alles klar  
            case DEFAULT
                stop 'unknown refinement mode; simple_commander_prime3D%exec_cont3D'
        end select
        if( p%xfel .eq. 'yes' )then
            if( cline%defined('msk') .or. cline%defined('mw') .or.&
            cline%defined('nvox') .or. cline%defined('automsk') )then
                stop 'no mask allowed when processing XFEL patterns; simple_prime3D'
            endif
        endif
        call b%build_general_tbox(p, cline)
        call b%build_cont3D_tbox(p)
        startit = 1
        if( cline%defined('startit') )startit = p%startit
        if( cline%defined('part') )then
            if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
            call cont3D_exec(b, p, cline, startit, converged)   
        else
            startit = 1
            if( cline%defined('startit') )startit = p%startit
            do i=startit,p%maxits
                call cont3D_exec(b, p, cline, i, converged)
                if(converged) exit
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CONT3D NORMAL STOP ****')
    end subroutine exec_cont3D
    !> CHECK3D_CONV is a SIMPLE program
    subroutine exec_check3D_conv( self, cline )
        use simple_math,    only: rad2deg, get_lplim
        class(check3D_conv_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)      :: p
        type(build)       :: b
        real, allocatable :: maplp(:)
        integer           :: istate, loc(1)
        logical           :: here, limset, converged, update_res
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        ! nstates consistency check
        ! Now incompatible with empty states
        !if( cline%defined('nstates') )then
        !    if( p%nstates /= b%a%get_nstates() ) stop 'Inconsistent number of states between command-line and oritab'
        !endif
        limset = .false.
        if( p%eo .eq. 'yes' )then
            allocate( maplp(p%nstates) )
            maplp = 0.
            do istate=1,p%nstates
                if( b%a%get_statepop( istate ) == 0 )cycle ! empty state
                p%fsc = 'fsc_state'//int2str_pad(istate,2)//'.bin'
                inquire(file=p%fsc, exist=here)
                if( here )then
                    b%fsc(istate,:) = file2rarr(p%fsc)
                    maplp(istate)   = max(b%img%get_lp(get_lplim(b%fsc(istate,:))),2.*p%smpd)
                else
                    write(*,*) 'Tried to open the fsc file: ', trim(p%fsc)
                    stop 'but it does not exist!'
                endif
            enddo
            loc     = maxloc( maplp )
            p%state = loc(1)            ! state with worst low-pass
            p%lp    = maplp( p%state )  ! worst lp
            p%fsc   = 'fsc_state'//int2str_pad(p%state,2)//'.bin'
            deallocate(maplp)
            limset = .true.
        endif
        ! Let find override the command line input lp (if given)
        if( .not. limset .and. cline%defined('find') )then
            p%lp = b%img%get_lp(p%find)
            limset = .true.
        endif
        ! Method for setting lp with lowest priority is lp on the command line
        if( cline%defined('lp') ) limset = .true.
        ! If we arrived here and the limit wasn't set: fall over
        if( limset )then
            ! we are happy
        else
            ! we fall over
            stop 'No method available to set low-pass limit! ABORTING...'
        endif
        ! calculate angular threshold
        select case(p%refine)
            case('yes','de')
                p%athres = max(p%lp, ATHRES_LIM)
            case DEFAULT
                p%athres = rad2deg(atan(p%lp/(p%moldiam/2.)))
        end select
        ! check convergence
        if( cline%defined('update_res') )then
            update_res = .false.
            if( cline%get_carg('update_res').eq.'yes' )update_res = .true.
            if( cline%get_carg('update_res').eq.'no' .and. p%refine.eq.'het')then
                converged = b%conv%check_conv_het()
            else
                select case(p%refine)
                    case('yes','de')
                        converged = b%conv%check_conv_cont3D()
                    case DEFAULT
                        converged = b%conv%check_conv3D()
                end select
            endif
        else
            select case(p%refine)
                case('het')
                    converged = b%conv%check_conv_het()
                case('yes','de')
                    converged = b%conv%check_conv_cont3D()
                case DEFAULT
                    converged = b%conv%check_conv3D()
            end select                 
        endif
        ! reports convergence, shift activation, resolution update and
        ! fraction of search space scanned to the distr commander
        if( p%doshift )then
            call cline%set('trs', p%trs)        ! activates shift search
        endif
        if( converged )then
            call cline%set('converged', 'yes')
        else
            call cline%set('converged', 'no')
        endif
        if( update_res )then
            call cline%set('update_res', 'yes') ! fourier index to be updated in distr commander
        else
            call cline%set('update_res', 'no')
        endif
        call cline%set('frac', b%conv%get('frac'))
        ! end gracefully
        call b%kill_general_tbox
        call simple_end('**** SIMPLE_CHECK3D_CONV NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check3D_conv

end module simple_commander_prime3D
