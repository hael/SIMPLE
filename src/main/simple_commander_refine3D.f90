! concrete commander: prime3D for ab initio 3D reconstruction and 3D refinement
module simple_commander_refine3D
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_ori,            only: ori
use simple_oris,           only: oris
use simple_commander_base, only: commander_base

implicit none

public :: npeaks_commander
public :: nspace_commander
public :: refine3D_init_commander
! public :: multiptcl_init_commander
public :: prime3D_commander
public :: check_3Dconv_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: npeaks_commander
  contains
    procedure :: execute      => exec_npeaks
end type npeaks_commander
type, extends(commander_base) :: nspace_commander
 contains
   procedure :: execute      => exec_nspace
end type nspace_commander
type, extends(commander_base) :: refine3D_init_commander
  contains
    procedure :: execute      => exec_refine3D_init
end type refine3D_init_commander
! type, extends(commander_base) :: multiptcl_init_commander
!   contains
!     procedure :: execute      => exec_multiptcl_init
! end type multiptcl_init_commander
type, extends(commander_base) :: prime3D_commander
  contains
    procedure :: execute      => exec_refine3D
end type prime3D_commander
type, extends(commander_base) :: check_3Dconv_commander
  contains
    procedure :: execute      => exec_check_3Dconv
end type check_3Dconv_commander

contains

    subroutine exec_npeaks( self, cline )
        class(npeaks_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(build)  :: b
        type(params) :: p
        integer :: npeaks
        p = params(cline) ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.)
        npeaks = min(10,b%e%find_npeaks(p%lp, p%moldiam))
        write(*,'(A,1X,I4)') '>>> NPEAKS:', npeaks
        ! end gracefully
        call simple_end('**** SIMPLE_NPEAKS NORMAL STOP ****')
    end subroutine exec_npeaks

    subroutine exec_nspace(self,cline)
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

    subroutine exec_refine3D_init( self, cline )
        use simple_qsys_funs,      only: qsys_job_finished
        class(refine3D_init_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)       :: p
        type(build)        :: b
        integer, parameter :: MAXIMGS=1000
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        call b%build_strategy3D_tbox(p)     ! prime3D objects built
        ! generate the random model
        if( cline%defined('nran') )then
            call gen_random_model( p%nran )
        else
            if( p%nptcls > MAXIMGS )then
                call gen_random_model( MAXIMGS )
            else
                call gen_random_model()
            endif
        endif
        ! end gracefully
        call qsys_job_finished( p, cline%get_carg('prg') )
        call simple_end('**** SIMPLE_REFINE3D_INIT NORMAL STOP ****', print_simple=.false.)

        contains

            subroutine gen_random_model( nsamp_in )
                use simple_kbinterpol,         only: kbinterpol
                use simple_prep4cgrid,         only: prep4cgrid
                use simple_strategy2D3D_common ! use all in there
                integer, optional, intent(in) :: nsamp_in  !< num input samples
                type(ran_tabu)       :: rt
                type(ori)            :: orientation
                type(kbinterpol)     :: kbwin
                type(prep4cgrid)     :: gridprep
                type(ctfparams)      :: ctfvars
                integer, allocatable :: sample(:)
                integer              :: i, nsamp

                ! init volumes
                call preprecvols(b, p)
                if( trim(p%refine).eq.'tseries' )then
                    call b%a%spiral
                else
                    call b%a%rnd_oris
                    call b%a%zero_shifts
                endif
                p%vols(1) = 'startvol'//p%ext
                nsamp = p%top - p%fromp + 1
                if( present(nsamp_in) ) nsamp = nsamp_in
                allocate( sample(nsamp) )
                if( present(nsamp_in) )then
                    rt = ran_tabu(p%top - p%fromp + 1)
                    call rt%ne_ran_iarr(sample)
                    call rt%kill
                else
                    forall(i=1:nsamp) sample(i) = i
                endif
                write(*,'(A)') '>>> RECONSTRUCTING RANDOM MODEL'
                ! make the gridding prepper
                kbwin = b%recvols(1)%get_kbwin()
                call gridprep%new(b%img, kbwin, [p%boxpd,p%boxpd,1])
                do i=1,nsamp
                    call progress(i, nsamp)
                    orientation = b%a%get_ori(sample(i) + p%fromp - 1)
                    ctfvars     = b%spproj%get_ctfparams(p%oritype, sample(i) + p%fromp - 1)
                    call read_img_and_norm( b, p, sample(i) + p%fromp - 1 )
                    call gridprep%prep(b%img, b%img_pad)
                    call b%recvols(1)%insert_fplane(b%se, orientation, ctfvars, b%img_pad, pwght=1.0)
                end do
                deallocate(sample)
                call norm_struct_facts(b, p)
                call killrecvols(b, p)
                if( p%part .ne. 1 )then
                    ! so random oris only written once in distributed mode
                else
                    ! update the spproj on disk
                    call b%spproj%write()
                endif
            end subroutine gen_random_model

    end subroutine exec_refine3D_init

    ! subroutine exec_multiptcl_init( self, cline )
    !     use simple_rec_master, only: exec_rec_master
    !     use simple_binoris_io, only: binwrite_oritab
    !     class(multiptcl_init_commander), intent(inout) :: self
    !     class(cmdline),                  intent(inout) :: cline
    !     type(params) :: p
    !     type(build)  :: b
    !     p = params(cline) ! constants & derived constants produced
    !     call b%build_general_tbox(p, cline)
    !     if( .not. cline%defined('oritab') )then
    !         call b%a%rnd_oris
    !         call b%a%zero_shifts
    !     endif
    !     if( cline%defined('state2split') )then
    !         if( cline%defined('oritab') )then
    !             p%nstates = b%a%get_n('state')
    !             call b%a%split_state(p%state2split)
    !             p%nstates = p%nstates + 1
    !         else
    !             stop 'Need oritab to be defined when state2split is defined on command line; simple_multiptcl_init'
    !         endif
    !     else if( p%tseries .eq. 'yes' )then
    !         call b%a%ini_tseries(p%nstates, 'state')
    !     else
    !         call b%a%rnd_states(p%nstates)
    !         if( p%nstates < 2 ) stop 'Nonsensical to have nstates < 2; simple_multiptcl_init'
    !     endif
    !     if( p%norec .ne. 'yes' )then
    !         if( cline%defined('lp') )then
    !             call b%build_rec_tbox(p)
    !             p%eo = 'no'
    !         else
    !             call b%build_rec_eo_tbox(p)
    !             p%eo = 'yes'
    !         endif
    !         call exec_rec_master(b, p, cline, 'startvol')
    !     endif
    !     if( p%zero .eq. 'yes' ) call b%a%set_all2single('corr', 0.)
    !     call binwrite_oritab('multiptcl_startdoc'//trim(METADATA_EXT), b%spproj, b%a, [1,b%a%get_noris()])
    !     ! end gracefully
    !     call simple_end('**** SIMPLE_MULTIPTCL_INIT NORMAL STOP ****', print_simple=.false.)
    ! end subroutine exec_multiptcl_init

    subroutine exec_refine3D( self, cline )
        use simple_strategy3D_matcher, only: refine3D_exec
        class(prime3D_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)               :: p
        type(build)                :: b
        integer                    :: i, startit
        logical                    :: converged
        real                       :: corr, corr_prev
        converged  = .false.
        p = params(cline) ! parameters generated
        if( p%neigh.eq.'yes' .and. .not. cline%defined('oritab') )then
            stop 'need oritab input for execution of prime3D with this refine mode'
        endif
        call b%build_general_tbox(p, cline)   ! general objects built
        if( .not. cline%defined('eo') ) p%eo = 'no' ! default
        if( cline%defined('lp') .or. cline%defined('find').or. p%eo .ne. 'no')then
            ! alles ok!
        else
           stop 'need a starting low-pass limit (set lp or find)!'
        endif
        call b%build_strategy3D_tbox(p) ! prime3D objects built
        startit = 1
        if( cline%defined('startit') )startit = p%startit
        if( startit == 1 )call b%a%clean_updatecnt
        if( p%l_distr_exec )then
            if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
            if( cline%defined('find') )then
                p%lp = calc_lowpass_lim( p%find, p%boxmatch, p%smpd )
            endif
            call refine3D_exec(b, p, cline, startit, converged) ! partition or not, depending on 'part'
        else
            p%find = calc_fourier_index( p%lp, p%boxmatch, p%smpd )
            ! init extremal dynamics
            if( cline%defined('extr_iter') )then
                ! all is well
            else
                p%extr_iter = startit
            endif
            corr = -1
            do i=startit,p%maxits
                call refine3D_exec(b, p, cline, i, converged)
                ! updates extremal iteration
                p%extr_iter = p%extr_iter + 1
                if( .not. p%l_distr_exec .and. p%refine .eq. 'snhc' .and. .not. cline%defined('szsn') )then
                    ! update stochastic neighborhood size if corr is not improving
                    corr_prev = corr
                    corr      = b%a%get_avg('corr')
                    if( i > 1 .and. corr <= corr_prev )then
                        p%szsn = min(SZSN_MAX,p%szsn + SZSN_STEP)
                    endif
                endif
                if( converged )exit
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_PRIME3D NORMAL STOP ****')
    end subroutine exec_refine3D

    subroutine exec_check_3Dconv( self, cline )
        class(check_3Dconv_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)      :: p
        type(build)       :: b
        real, allocatable :: maplp(:)
        integer           :: istate, loc(1)
        logical           :: limset, converged, update_res
        p = params(cline) ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        limset = .false. ;  update_res = .false.
        if( p%eo .ne. 'no' )then
            allocate( maplp(p%nstates), stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("In simple_commander_refine3D:: exec_check3D_conv", alloc_stat)
            maplp = 0.
            do istate=1,p%nstates
                if( b%a%get_pop( istate, 'state' ) == 0 )cycle ! empty state
                p%fsc = 'fsc_state'//int2str_pad(istate,2)//'.bin'
                if( file_exists(p%fsc) )then
                    b%fsc(istate,:) = file2rarr(p%fsc)
                    maplp(istate)   = max(b%img%get_lp(get_lplim_at_corr(b%fsc(istate,:),p%lplim_crit)),2.*p%smpd)
                else
                    write(*,*) 'Tried to check the fsc file: ', trim(p%fsc)
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
        ! check convergence
        if( cline%defined('update_res') )then
            update_res = .false.
            if( cline%get_carg('update_res').eq.'yes' )update_res = .true.
            if( cline%get_carg('update_res').eq.'no' .and. str_has_substr(p%refine,'cluster') )then
                converged = b%conv%check_conv_cluster()
            else
                converged = b%conv%check_conv3D()
            endif
        else
            select case(p%refine)
            case('cluster','clustersym')
                    converged = b%conv%check_conv_cluster()
                case DEFAULT
                    converged = b%conv%check_conv3D()
            end select
        endif
        ! reports convergence, shift activation, resolution update and
        ! fraction of search space scanned to the distr commander
        if( p%l_doshift )then
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
        call simple_end('**** SIMPLE_CHECK_3DCONV NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check_3Dconv

end module simple_commander_refine3D
