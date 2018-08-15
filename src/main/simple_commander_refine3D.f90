! concrete commander: refine3D for ab initio 3D reconstruction and 3D refinement
module simple_commander_refine3D
include 'simple_lib.f08'
use simple_parameters,     only: parameters
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_ori,            only: ori
use simple_oris,           only: oris
use simple_commander_base, only: commander_base
implicit none

public :: nspace_commander
public :: refine3D_init_commander
public :: refine3D_commander
public :: check_3Dconv_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: nspace_commander
 contains
   procedure :: execute      => exec_nspace
end type nspace_commander
type, extends(commander_base) :: refine3D_init_commander
  contains
    procedure :: execute      => exec_refine3D_init
end type refine3D_init_commander
type, extends(commander_base) :: refine3D_commander
  contains
    procedure :: execute      => exec_refine3D
end type refine3D_commander
type, extends(commander_base) :: check_3Dconv_commander
  contains
    procedure :: execute      => exec_check_3Dconv
end type check_3Dconv_commander

contains

    subroutine exec_nspace(self,cline)
        class(nspace_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters) :: params
        type(oris)       :: o
        real             :: ares
        integer          :: i
        call params%new(cline)
        do i=500,5000,500
            o = oris(i)
            call o%spiral
            ares = o%find_angres()
            write(*,'(A,1X,I7,1X,A,1X,F5.2)') 'NR OF PROJDIRS:', i, 'RESOLUTION:', resang(ares, params%moldiam)
        end do
        call simple_end('**** SIMPLE_NSPACE NORMAL STOP ****')
    end subroutine exec_nspace

    subroutine exec_refine3D_init( self, cline )
        use simple_qsys_funs,      only: qsys_job_finished
        class(refine3D_init_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)   :: params
        type(builder)      :: build
        integer, parameter :: MAXIMGS=1000
        call build%init_params_and_build_strategy3D_tbox(cline, params)
        ! generate the random model
        if( cline%defined('nran') )then
            call gen_random_model( params%nran )
        else
            if( params%nptcls > MAXIMGS )then
                call gen_random_model( MAXIMGS )
            else
                call gen_random_model()
            endif
        endif
        ! end gracefully
        call qsys_job_finished( cline%get_carg('prg') )
        call simple_end('**** SIMPLE_REFINE3D_INIT NORMAL STOP ****', print_simple=.false.)

        contains

            subroutine gen_random_model( nsamp_in )
                use simple_strategy2D3D_common ! use all in there
                use simple_image,              only: image
                integer, optional, intent(in) :: nsamp_in  !< num input samples
                type(ran_tabu)       :: rt
                type(ori)            :: orientation
                type(ctfparams)      :: ctfvars
                type(image)          :: mskimg
                integer, allocatable :: sample(:)
                integer              :: i, nsamp, ind
                ! init volumes
                call preprecvols()
                if( trim(params%refine).eq.'tseries' )then
                    call build%spproj_field%spiral
                else
                    call build%spproj_field%rnd_oris
                    call build%spproj_field%zero_shifts
                endif
                params%vols(1) = 'startvol'//params%ext
                nsamp = params%top - params%fromp + 1
                if( present(nsamp_in) ) nsamp = nsamp_in
                allocate( sample(nsamp) )
                if( present(nsamp_in) )then
                    rt = ran_tabu(params%top - params%fromp + 1)
                    call rt%ne_ran_iarr(sample)
                    call rt%kill
                else
                    forall(i=1:nsamp) sample(i) = i
                endif
                write(*,'(A)') '>>> RECONSTRUCTING RANDOM MODEL'
                call mskimg%disc(build%img%get_ldim(), params%smpd, params%msk)
                do i=1,nsamp
                    call progress(i, nsamp)
                    ind         = sample(i) + params%fromp - 1
                    orientation = build%spproj_field%get_ori(ind)
                    ctfvars     = build%spproj%get_ctfparams(params%oritype, ind)
                    call read_img(ind)
                    call build%img%norm_subtr_backgr_pad_fft(mskimg, build%img_pad)
                    call build%recvols(1)%insert_fplane(build%pgrpsyms, orientation, ctfvars, build%img_pad, pwght=1.0)
                end do
                deallocate(sample)
                call norm_struct_facts()
                call killrecvols()
                if( params%part .ne. 1 )then
                    ! so random oris only written once in distributed mode
                else
                    ! update the spproj on disk
                    call build%spproj%write_segment_inside(params%oritype)
                endif
                call mskimg%kill
            end subroutine gen_random_model

    end subroutine exec_refine3D_init

    subroutine exec_refine3D( self, cline )
        use simple_strategy3D_matcher, only: refine3D_exec
        class(refine3D_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer :: startit
        logical :: converged
        call build%init_params_and_build_strategy3D_tbox(cline,params)
        if( cline%defined('lp') .or. cline%defined('find').or. params%eo .ne. 'no')then
            ! alles ok!
        else
           stop 'need a starting low-pass limit (set lp or find)!'
        endif
        startit = 1
        if( cline%defined('startit') )startit = params%startit
        if( startit == 1 ) call build%spproj_field%clean_updatecnt
        if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
        call refine3D_exec( cline, startit, converged) ! partition or not, depending on 'part'
        ! end gracefully
        call simple_end('**** SIMPLE_REFINE3D NORMAL STOP ****')
    end subroutine exec_refine3D

    subroutine exec_check_3Dconv( self, cline )
        use simple_convergence, only: convergence
        use simple_parameters,  only: params_glob
        class(check_3Dconv_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)  :: params
        type(builder)     :: build
        type(convergence) :: conv
        real, allocatable :: maplp(:)
        integer           :: istate, loc(1)
        logical           :: limset, converged, update_res
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        limset = .false. ;  update_res = .false.
        if( params%eo .ne. 'no' )then
            allocate( maplp(params%nstates), stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("In simple_commander_refine3D:: exec_check3D_conv", alloc_stat)
            maplp = 0.
            do istate=1,params%nstates
                if( build%spproj_field%get_pop( istate, 'state' ) == 0 )cycle ! empty state
                params%fsc = 'fsc_state'//int2str_pad(istate,2)//'.bin'
                if( file_exists(params%fsc) )then
                    build%fsc(istate,:) = file2rarr(params%fsc)
                    maplp(istate)   = max(build%img%get_lp(get_lplim_at_corr(build%fsc(istate,:),params%lplim_crit)),2.*params%smpd)
                else
                    write(*,*) 'Tried to check the fsc file: ', trim(params%fsc)
                    stop 'but it does not exist!'
                endif
            enddo
            loc     = maxloc( maplp )
            params%state = loc(1)            ! state with worst low-pass
            params%lp    = maplp( params%state )  ! worst lp
            params%fsc   = 'fsc_state'//int2str_pad(params%state,2)//'.bin'
            deallocate(maplp)
            limset = .true.
        endif
        ! Let find override the command line input lp (if given)
        if( .not. limset .and. cline%defined('find') )then
            params%lp = build%img%get_lp(params%find)
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
            if( cline%get_carg('update_res').eq.'no' .and. str_has_substr(params%refine,'cluster') )then
                converged = conv%check_conv_cluster(cline)
            else
                converged = conv%check_conv3D(cline, params%msk)
            endif
        else
            select case(params%refine)
            case('cluster','clustersym')
                    converged = conv%check_conv_cluster(cline)
                case DEFAULT
                    converged = conv%check_conv3D(cline, params%msk)
            end select
        endif
        ! reports convergence, shift activation, resolution update and
        ! fraction of search space scanned to the distr commander
        if( params_glob%l_doshift )then
            call cline%set('trs', params_glob%trs)        ! activates shift search
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
        call cline%set('frac', conv%get('frac'))
        ! end gracefully
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_CHECK_3DCONV NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check_3Dconv

end module simple_commander_refine3D
