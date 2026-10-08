!@descr: reconstruction service: state maps from explicit particle rows, for reconstruct3D and its callers
! A request names the particle rows, the weight source, the backend (gridding or kernel PCG), the
! output kind (raw half maps or delivered maps), a file-name prefix, the publication policy and the
! dispatch (in this process or the standard queue rounds). The caller owns sampling, sigma2 setup and
! the project; the service owns insertion, assembly, post-processing, naming and registration.
module simple_rec3D_service
use simple_core_module_api
use simple_builder,              only: builder
use simple_parameters,           only: parameters
use simple_cmdline,              only: cmdline
use simple_qsys_env,             only: qsys_env
use simple_matcher_3Drec,        only: calc_3Drec
use simple_commanders_rec_distr, only: commander_volassemble, filter_pcg_nonuniform_maps
use simple_refine3D_fnames,      only: refine3D_fsc_fname, refine3D_state_vol_fname, refine3D_pcg_raw_accum_fname
use simple_rec3D_pcg_strategy,   only: execute_rec3D_pcg_worker, execute_rec3D_pcg_distributed_master
use simple_sigma2_files,         only: load_sigma2_groups
use simple_state_weight_set,     only: state_weight_set, STATE_WEIGHTS_KIND_PARTITION
implicit none

public :: rec3D_request, rec3D_service
public :: rec3D_backend_id, rec3D_backend_is_wired
public :: REC3D_BACKEND_INVALID, REC3D_BACKEND_GRIDDING, REC3D_BACKEND_PCG
public :: REC3D_WEIGHTS_HARD, REC3D_WEIGHTS_SET, REC3D_WEIGHTS_TABLE, REC3D_OUTPUT_RAW, REC3D_OUTPUT_DELIVERED
public :: REC3D_DISPATCH_INPROC, REC3D_DISPATCH_QUEUE
private
#include "simple_local_flags.inc"

integer,          parameter :: REC3D_BACKEND_INVALID  = 0
integer,          parameter :: REC3D_BACKEND_GRIDDING = 1
integer,          parameter :: REC3D_BACKEND_PCG      = 2
!> weight source: the hard state labels of the particle field
integer,          parameter :: REC3D_WEIGHTS_HARD     = 1
!> weight source: the state weight set the project registers (m_estimator=flex; workers open it themselves)
integer,          parameter :: REC3D_WEIGHTS_SET      = 2
!> weight source: a transient table the caller owns (FLEX's trial weights; in-process dispatch only)
integer,          parameter :: REC3D_WEIGHTS_TABLE    = 3
!> raw: assembled half maps and merged map only; delivered: also post-processing
integer,          parameter :: REC3D_OUTPUT_RAW       = 1
integer,          parameter :: REC3D_OUTPUT_DELIVERED = 2
integer,          parameter :: REC3D_DISPATCH_INPROC  = 1
integer,          parameter :: REC3D_DISPATCH_QUEUE   = 2
character(len=*), parameter :: PINDS_FNAME            = 'rec3D_service_pinds.txt'

!> What to reconstruct and how
type :: rec3D_request
    integer, allocatable :: pinds(:)                         !< particle rows of the field, increasing
    real,    allocatable :: weights_table(:,:)               !< REC3D_WEIGHTS_TABLE: (size(pinds), nstates)
    integer              :: weights  = REC3D_WEIGHTS_HARD
    integer              :: backend  = REC3D_BACKEND_GRIDDING
    integer              :: output   = REC3D_OUTPUT_DELIVERED
    type(string)         :: prefix                           !< blank: the canonical recvol_stateNN names
    logical              :: register = .false.               !< state maps and FSCs into the project out segment
    integer              :: dispatch = REC3D_DISPATCH_INPROC
end type rec3D_request

!> The service; the queue environment lives here when the dispatch is the queue
type :: rec3D_service
    private
    type(qsys_env), allocatable :: qenv
    type(chash),    allocatable :: job_descr
    integer                     :: dispatch    = REC3D_DISPATCH_INPROC
    integer                     :: nthr_master = 1
    logical                     :: exists      = .false.
  contains
    procedure :: new
    procedure :: execute
    procedure :: kill
end type rec3D_service

contains

    pure integer function rec3D_backend_id(name) result(backend_id)
        character(len=*), intent(in) :: name
        select case(trim(name))
            case('gridding')
                backend_id = REC3D_BACKEND_GRIDDING
            case('pcg')
                backend_id = REC3D_BACKEND_PCG
            case DEFAULT
                backend_id = REC3D_BACKEND_INVALID
        end select
    end function rec3D_backend_id

    pure logical function rec3D_backend_is_wired(backend_id) result(l_wired)
        integer, intent(in) :: backend_id
        l_wired = backend_id == REC3D_BACKEND_GRIDDING .or. backend_id == REC3D_BACKEND_PCG
    end function rec3D_backend_is_wired

    !> For the queue dispatch, the partitions (balanced over the particles with state > 0) and the
    !! job description of cline; nthr_master is the thread count of master-side assembly
    subroutine new( self, params, build, cline, dispatch, nthr_master )
        class(rec3D_service), intent(inout) :: self
        type(parameters),     intent(inout) :: params
        type(builder),        intent(inout) :: build
        class(cmdline),       intent(inout) :: cline
        integer,              intent(in)    :: dispatch, nthr_master
        call self%kill
        self%dispatch    = dispatch
        self%nthr_master = nthr_master
        select case(dispatch)
            case(REC3D_DISPATCH_INPROC)
            case(REC3D_DISPATCH_QUEUE)
                allocate(self%qenv, self%job_descr)
                call self%qenv%new(params, params%nparts, l_active=build%spproj_field%included())
                call cline%gen_job_descr(self%job_descr)
            case DEFAULT
                THROW_HARD('unknown reconstruction dispatch; rec3D_service%new')
        end select
        self%exists = .true.
    end subroutine new

    !> Reconstruct every state of params from the requested rows, then name, post-process and register
    subroutine execute( self, params, build, cline, request )
        class(rec3D_service), intent(inout) :: self
        type(parameters),     intent(inout) :: params
        type(builder),        intent(inout) :: build
        class(cmdline),       intent(inout) :: cline
        type(rec3D_request),  intent(in)    :: request
        if( .not. self%exists ) THROW_HARD('service not constructed; rec3D_service%execute')
        if( request%dispatch /= self%dispatch ) THROW_HARD('request dispatch differs from the service; rec3D_service%execute')
        select case(request%weights)
            case(REC3D_WEIGHTS_HARD)
            case(REC3D_WEIGHTS_SET)
                if( .not. params%l_m_estimator_flex ) THROW_HARD('a weight set request needs m_estimator=flex; rec3D_service%execute')
            case(REC3D_WEIGHTS_TABLE)
                if( request%dispatch /= REC3D_DISPATCH_INPROC ) THROW_HARD('a weight table is reconstructed in process only; rec3D_service%execute')
                if( .not. allocated(request%weights_table) ) THROW_HARD('weight table request without a table; rec3D_service%execute')
            case DEFAULT
                THROW_HARD('unsupported weight source; rec3D_service%execute')
        end select
        if( .not. rec3D_backend_is_wired(request%backend) ) THROW_HARD('unsupported backend; rec3D_service%execute')
        if( .not. allocated(request%pinds) ) THROW_HARD('request without particle rows; rec3D_service%execute')
        if( request%prefix%strlen_trim() > 0 .and. (request%register .or. request%output /= REC3D_OUTPUT_RAW) )then
            THROW_HARD('a file-name prefix applies to unregistered raw output only; rec3D_service%execute')
        endif
        if( params%l_state_defined )then
            if( params%state < 1 .or. params%state > params%nstates ) THROW_HARD('state= lies outside 1..nstates; rec3D_service%execute')
        endif
        if( request%weights == REC3D_WEIGHTS_SET ) call report_weight_set(params, build)
        select case(self%dispatch)
            case(REC3D_DISPATCH_INPROC)
                call execute_inproc(self, params, build, cline, request)
            case(REC3D_DISPATCH_QUEUE)
                call execute_queue(self, params, build, cline, request)
        end select
        if( request%register ) call register_outputs(params, build)
        if( request%output == REC3D_OUTPUT_DELIVERED ) call postprocess_states(params, build, cline)
        if( request%prefix%strlen_trim() > 0 ) call rename_outputs(params, request%prefix)
    end subroutine execute

    subroutine kill( self )
        class(rec3D_service), intent(inout) :: self
        if( allocated(self%qenv) )then
            call self%qenv%kill
            deallocate(self%qenv)
        endif
        if( allocated(self%job_descr) )then
            call self%job_descr%kill
            deallocate(self%job_descr)
        endif
        self%dispatch    = REC3D_DISPATCH_INPROC
        self%nthr_master = 1
        self%exists      = .false.
    end subroutine kill

    ! ---- private ----

    subroutine execute_inproc( self, params, build, cline, request )
        class(rec3D_service), intent(inout) :: self
        type(parameters),     intent(inout) :: params
        type(builder),        intent(inout) :: build
        class(cmdline),       intent(inout) :: cline
        type(rec3D_request),  intent(in)    :: request
        type(state_weight_set) :: wset
        type(string) :: volname
        integer      :: state
        logical      :: l_sigma_loaded, l_set, l_table
        l_set   = request%weights == REC3D_WEIGHTS_SET
        l_table = request%weights == REC3D_WEIGHTS_TABLE
        if( l_set ) call wset%new(build%spproj, build%spproj_field)
        select case(request%backend)
            case(REC3D_BACKEND_PCG)
                call remove_pcg_raw_files(params)
                if( l_set )then
                    call execute_rec3D_pcg_worker(params, build, cline, request%pinds, wset)
                else if( l_table )then
                    call execute_rec3D_pcg_worker(params, build, cline, request%pinds, wtab=request%weights_table)
                else
                    call execute_rec3D_pcg_worker(params, build, cline, request%pinds)
                endif
                call assemble_pcg(params, build, cline)
            case(REC3D_BACKEND_GRIDDING)
                ! sigma weighting belongs to the Euclidean data objective; ML regularization is a separate
                ! FSC/SSNR prior applied by volassemble
                ! (a table caller works in its own process and has its noise model loaded already)
                if( params%cc_objfun == OBJFUN_EUCLID .and. .not. (l_table .and. allocated(build%esig%sigma2_noise)) )then
                    call load_sigma2_groups(params, build%pftc, build%esig, build%spproj, build%spproj_field, &
                        &l_sigma_loaded)
                    if( .not. l_sigma_loaded ) THROW_HARD('gridding objfun=euclid requires sigma2 files')
                endif
                if( l_set )then
                    call calc_3Drec(params, build, size(request%pinds), request%pinds, wset)
                else if( l_table )then
                    call calc_3Drec(params, build, size(request%pinds), request%pinds, wtab=request%weights_table)
                else
                    call calc_3Drec(params, build, size(request%pinds), request%pinds)
                endif
                call assemble_gridding(params, cline, params%nthr)
                do state = 1, params%nstates
                    volname = refine3D_state_vol_fname(state)
                    params%vols(state) = volname
                    call cline%set('vol'//int2str(state), volname)
                end do
                call volname%kill
        end select
        call wset%kill
    end subroutine execute_inproc

    !> The parts run the private reconstruct3D worker on their range of the requested rows (pindfile)
    subroutine execute_queue( self, params, build, cline, request )
        class(rec3D_service), intent(inout) :: self
        type(parameters),     intent(inout) :: params
        type(builder),        intent(inout) :: build
        class(cmdline),       intent(inout) :: cline
        type(rec3D_request),  intent(in)    :: request
        type(string) :: pindfile
        if( request%backend == REC3D_BACKEND_PCG ) call remove_pcg_raw_files(params)
        call write_pinds(request%pinds)
        pindfile = simple_abspath(PINDS_FNAME)
        call self%job_descr%set('pindfile', pindfile%to_char())
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        call self%job_descr%delete('pindfile')
        call del_file(PINDS_FNAME)
        call pindfile%kill
        if( request%backend == REC3D_BACKEND_PCG )then
            call assemble_pcg(params, build, cline)
        else
            call assemble_gridding(params, cline, self%nthr_master)
        endif
    end subroutine execute_queue

    !> volassemble on the gridding partials; a vol<N> naming this run's own output is dropped when the
    !! file does not exist yet, so the assembly does not take it as a previous reference
    subroutine assemble_gridding( params, cline, nthr )
        type(parameters), intent(in)    :: params
        class(cmdline),   intent(inout) :: cline
        integer,          intent(in)    :: nthr
        type(commander_volassemble) :: xvolassemble
        type(cmdline)               :: cline_volassemble
        type(string)                :: volname, vol_in
        integer                     :: state
        cline_volassemble = cline
        call cline_volassemble%set('prg',  'volassemble')
        call cline_volassemble%set('nthr', nthr)
        do state = 1, params%nstates
            volname = refine3D_state_vol_fname(state)
            if( cline_volassemble%defined('vol'//int2str(state)) )then
                vol_in = cline_volassemble%get_carg('vol'//int2str(state))
                if( trim(vol_in%to_char()) == trim(volname%to_char()) )then
                    if( .not. file_exists(volname) ) call cline_volassemble%delete('vol'//int2str(state))
                endif
            endif
        end do
        call xvolassemble%execute(cline_volassemble)
        call cline_volassemble%kill
        call volname%kill
        call vol_in%kill
    end subroutine assemble_gridding

    !> The PCG master: reduce the parts' raw accumulators and solve; with nonuniform filtering also the
    !! NU competition and the matching low-pass handoff, so the outputs equal a refinement iteration's
    subroutine assemble_pcg( params, build, cline )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        real, allocatable :: nu_align_lps(:)
        if( params%l_nonuniform )then
            allocate(nu_align_lps(params%nstates), source=0.0)
            call execute_rec3D_pcg_distributed_master(params, build, cline, nu_align_lps=nu_align_lps)
            call filter_pcg_nonuniform_maps(params, build, cline%defined('frozen_rec'), nu_align_lps)
            deallocate(nu_align_lps)
        else
            call execute_rec3D_pcg_distributed_master(params, build, cline)
        endif
    end subroutine assemble_pcg

    !> A stale raw accumulator must never pass for a completed worker of this launch; workers publish
    !! through .tmp and an atomic rename
    subroutine remove_pcg_raw_files( params )
        type(parameters), intent(in) :: params
        type(string) :: raw_fname
        integer      :: state, part, eo
        do state = 1, params%nstates
            do eo = 0, 1
                do part = 1, params%nparts
                    raw_fname = refine3D_pcg_raw_accum_fname(state, part, params%numlen, &
                        &merge('odd ', 'even', eo == 1))
                    call del_file(raw_fname)
                    call del_file(raw_fname//'.tmp')
                enddo
            enddo
        enddo
        call raw_fname%kill
    end subroutine remove_pcg_raw_files

    !> The state maps and FSCs into the project out segment. A state this run did not reconstruct (no
    !! particles, or no weight) has no FSC file: assembly dropped it, and its entries are left as they are
    subroutine register_outputs( params, build )
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        type(string) :: fsc_file
        integer      :: state
        do state = 1, params%nstates
            if( params%l_state_defined .and. state /= params%state ) cycle
            fsc_file = refine3D_fsc_fname(state)
            if( .not. file_exists(fsc_file) .or. .not. file_exists(refine3D_state_vol_fname(state)) )then
                write(logfhandle,'(A,I0,A)') '>>> RECONSTRUCT3D: state ', state, ' was not reconstructed; not registered'
                call fsc_file%kill
                cycle
            endif
            call build%spproj%add_fsc2os_out(fsc_file, state, params%box_crop)
            if( trim(params%oritype).eq.'cls3D' )then
                call build%spproj%add_vol2os_out(refine3D_state_vol_fname(state), &
                    &params%smpd_crop, state, 'vol_cavg')
            else
                call build%spproj%add_vol2os_out(refine3D_state_vol_fname(state), &
                    &params%smpd_crop, state, 'vol')
            endif
            call fsc_file%kill
        enddo
        call build%spproj%write_segment_inside('out', params%projfile)
    end subroutine register_outputs

    !> Post-processing of every populated state map (postprocess=yes; not in a worker part)
    subroutine postprocess_states( params, build, cline )
        use simple_commanders_volops, only: postprocess_volume_from_files
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(state_weight_set) :: wset
        type(string) :: fname_vol, fname_fsc, fsc_saved
        real         :: bfac_saved
        integer      :: state, nptcls, ldim(3)
        if( trim(params%postprocess) /= 'yes' )return
        if( cline%defined('part') )return
        if( params%l_m_estimator_flex ) call wset%new(build%spproj, build%spproj_field)
        if( params%l_nonuniform )then
            write(logfhandle,'(A)') &
                &'>>> reconstruct3D postprocess: using classical postprocessing'
        endif
        ! postprocessing sets fsc and bfac in params for each state; the caller's values are kept
        fsc_saved  = params%fsc
        bfac_saved = params%bfac
        do state = 1, params%nstates
            if( params%l_state_defined .and. state /= params%state ) cycle
            ! a state without particles (hard labels) or without weight (weight set) is not post-processed
            if( params%l_m_estimator_flex )then
                if( wset%get_mass(state) <= 0._dp ) cycle
            else if( .not. cline%defined('frozen_rec') )then
                if( build%spproj_field%get_pop(state, 'state') == 0 ) cycle
            endif
            fname_vol = refine3D_state_vol_fname(state)
            fname_fsc = refine3D_fsc_fname(state)
            if( .not. file_exists(fname_vol) )then
                call fname_vol%kill
                call fname_fsc%kill
                cycle
            endif
            call find_ldim_nptcls(fname_vol, ldim, nptcls)
            params%fsc  = fsc_saved
            params%bfac = bfac_saved
            call postprocess_volume_from_files(fname_vol, fname_fsc, ldim(1), params%smpd_crop, params, cline, state)
            call fname_vol%kill
            call fname_fsc%kill
        enddo
        params%fsc  = fsc_saved
        params%bfac = bfac_saved
        call fsc_saved%kill
        call wset%kill
    end subroutine postprocess_states

    !> the weight set a weighted reconstruction uses, recorded in its log (D3: the kind with every result)
    subroutine report_weight_set( params, build )
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        type(state_weight_set) :: wset
        character(len=9) :: kind_name
        integer :: state
        call wset%new(build%spproj, build%spproj_field)
        kind_name = 'kernel'
        if( wset%get_kind() == STATE_WEIGHTS_KIND_PARTITION ) kind_name = 'partition'
        write(logfhandle,'(A,A,A,I0,A,A)') '>>> RECONSTRUCTION WEIGHTS: state weight set of kind ', trim(kind_name), &
            &', generation ', wset%get_generation(), ', producer ', wset%get_producer()
        do state = 1, params%nstates
            if( params%l_state_defined .and. state /= params%state ) cycle
            write(logfhandle,'(A,I0,A,F14.3,A,F14.3,A,I0)') '>>>   STATE ', state, ' APPLIED MASS ', wset%get_mass(state), &
                &' EFFECTIVE SAMPLE SIZE ', wset%get_ess(state), ' HARD POPULATION ', wset%get_pop(state)
        enddo
        call wset%kill
    end subroutine report_weight_set

    !> Caller-owned names of raw output: <prefix>_stateNN{,_even,_odd,_even_unfil,_odd_unfil}.mrc and <prefix>_fsc_stateNN.bin
    subroutine rename_outputs( params, prefix )
        type(parameters), intent(in) :: params
        class(string),    intent(in) :: prefix
        type(string) :: src, dst, stem, fsc_src, fsc_dst
        integer      :: state, i
        character(len=11), parameter :: SUFFIXES(5) = [character(len=11) :: '', '_even', '_odd', '_even_unfil', '_odd_unfil']
        do state = 1, params%nstates
            stem = prefix//'_state'//int2str_pad(state,2)
            do i = 1, size(SUFFIXES)
                src = add2fbody(refine3D_state_vol_fname(state), MRC_EXT, trim(SUFFIXES(i)))
                dst = stem//trim(SUFFIXES(i))//MRC_EXT
                if( file_exists(src) ) call simple_rename(src, dst)
            enddo
            fsc_src = refine3D_fsc_fname(state)
            fsc_dst = prefix//'_fsc_state'//int2str_pad(state,2)//'.bin'
            if( file_exists(fsc_src) ) call simple_rename(fsc_src, fsc_dst)
        enddo
        call src%kill
        call dst%kill
        call stem%kill
        call fsc_src%kill
        call fsc_dst%kill
    end subroutine rename_outputs

    subroutine write_pinds( pinds )
        integer, intent(in) :: pinds(:)
        integer :: funit, i, io_stat
        call fopen(funit, file=string(PINDS_FNAME), status='REPLACE', action='WRITE', iostat=io_stat)
        call fileiochk('rec3D_service: cannot write '//PINDS_FNAME, io_stat)
        do i = 1, size(pinds)
            write(funit,'(I0)') pinds(i)
        enddo
        call fclose(funit)
    end subroutine write_pinds

end module simple_rec3D_service
