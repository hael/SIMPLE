!@descr: reconstruct3D execution strategies (shared memory, kernel PCG, distributed master) around the reconstruction service
! A strategy parses the parameters, partitions even/odd, prepares the canonical sigma2 state, samples
! the particle rows and hands them to simple_rec3D_service, which reconstructs, assembles, post-processes
! and registers the state maps.
module simple_rec3D_strategy
use simple_core_module_api
use simple_builder,           only: builder
use simple_parameters,        only: parameters
use simple_cmdline,           only: cmdline
use simple_rec3D_service,     only: rec3D_service, rec3D_request, rec3D_backend_id, REC3D_BACKEND_GRIDDING, &
    &REC3D_BACKEND_PCG, REC3D_DISPATCH_INPROC, REC3D_DISPATCH_QUEUE, REC3D_WEIGHTS_HARD, REC3D_WEIGHTS_SET
implicit none

public :: rec3D_strategy, rec3D_inmem_strategy, rec3D_distr_strategy, create_rec3D_strategy
public :: rec3D_pcg_inmem_strategy
private
#include "simple_local_flags.inc"

! --------------------------------------------------------------------
! Strategy interface
! --------------------------------------------------------------------

type, abstract :: rec3D_strategy
    type(rec3D_service), allocatable :: service
contains
    procedure(init_interface),    deferred :: initialize
    procedure(exec_interface),    deferred :: execute
    procedure(cleanup_interface), deferred :: cleanup
end type rec3D_strategy

! Shared-memory gridding
type, extends(rec3D_strategy) :: rec3D_inmem_strategy
contains
    procedure :: initialize => inmem_initialize
    procedure :: execute    => inmem_execute
    procedure :: cleanup    => inmem_cleanup
end type rec3D_inmem_strategy

! Shared-memory kernel PCG: the distributed route in one process (this
! process is the only worker, then the master), on the in-memory lifecycle
type, extends(rec3D_inmem_strategy) :: rec3D_pcg_inmem_strategy
contains
    procedure :: execute => pcg_inmem_execute
end type rec3D_pcg_inmem_strategy

! Distributed-memory master (gridding or kernel PCG)
type, extends(rec3D_strategy) :: rec3D_distr_strategy
    integer :: nthr_master = 1
contains
    procedure :: initialize => distr_initialize
    procedure :: execute    => distr_execute
    procedure :: cleanup    => distr_cleanup
end type rec3D_distr_strategy

abstract interface
    subroutine init_interface(self, params, build, cline)
        import :: rec3D_strategy, parameters, builder, cmdline
        class(rec3D_strategy), intent(inout) :: self
        type(parameters),      intent(inout) :: params
        type(builder),         intent(inout) :: build
        class(cmdline),        intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, build, cline)
        import :: rec3D_strategy, parameters, builder, cmdline
        class(rec3D_strategy), intent(inout) :: self
        type(parameters),      intent(inout) :: params
        type(builder),         intent(inout) :: build
        class(cmdline),        intent(inout) :: cline
    end subroutine exec_interface

    subroutine cleanup_interface(self, params, build, cline)
        import :: rec3D_strategy, parameters, builder, cmdline
        class(rec3D_strategy), intent(inout) :: self
        type(parameters),      intent(in)    :: params
        type(builder),         intent(inout) :: build
        class(cmdline),        intent(inout) :: cline
    end subroutine cleanup_interface
end interface

contains

    ! --------------------------------------------------------------------
    ! Strategy selection
    ! --------------------------------------------------------------------

    function create_rec3D_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(rec3D_strategy), allocatable :: strategy
        type(string) :: backend
        backend = string('gridding')
        if( cline%defined('rec_backend') ) backend = cline%get_carg('rec_backend')
        select case(rec3D_backend_id(backend%to_char()))
            case(REC3D_BACKEND_GRIDDING)
                ! Distributed master iff: nparts defined AND part not defined.
                ! Keep this branch identical to the pre-selector strategy choice.
                if( cline%defined('nparts') .and. (.not.cline%defined('part')) )then
                    allocate(rec3D_distr_strategy :: strategy)
                    if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') &
                        &'>>> DISTRIBUTED-MEMORY REC3D EXECUTION (gridding)'
                else
                    allocate(rec3D_inmem_strategy :: strategy)
                    if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') &
                        &'>>> SHARED-MEMORY REC3D EXECUTION (gridding)'
                endif
            case(REC3D_BACKEND_PCG)
                if( cline%defined('nparts') .and. (.not.cline%defined('part')) )then
                    if( cline%get_iarg('nparts') > 1 )then
                        allocate(rec3D_distr_strategy :: strategy)
                        if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') &
                            &'>>> DISTRIBUTED-MEMORY REC3D EXECUTION (kernel PCG)'
                    else
                        allocate(rec3D_pcg_inmem_strategy :: strategy)
                    endif
                else
                    allocate(rec3D_pcg_inmem_strategy :: strategy)
                endif
                if( L_VERBOSE_GLOB .and. .not. cline%defined('nparts') )then
                    write(logfhandle,'(A)') '>>> SHARED-MEMORY REC3D EXECUTION (kernel PCG)'
                endif
            case DEFAULT
                THROW_HARD('rec_backend must be gridding or pcg')
        end select
        call backend%kill
    end function create_rec3D_strategy

    ! =====================================================================
    ! SHARED-MEMORY IMPLEMENTATION
    ! =====================================================================

    subroutine inmem_initialize(self, params, build, cline)
        class(rec3D_inmem_strategy), intent(inout) :: self
        type(parameters),            intent(inout) :: params
        type(builder),               intent(inout) :: build
        class(cmdline),              intent(inout) :: cline
        call build%init_params_and_build_general_tbox(cline, params)
        call sync_resolved_rec_params(params, cline)
        call build%build_strategy3D_tbox(params)
        ! Even/odd partitioning
        if( build%spproj_field%get_nevenodd() == 0 ) call build%spproj_field%partition_eo
        ! Update eo flags in project
        call build%spproj%write_segment_inside(params%oritype)
        call ensure_canonical_sigma_state(params, build, cline)
        if( .not. allocated(self%service) ) allocate(self%service)
        call self%service%new(params, build, cline, REC3D_DISPATCH_INPROC, params%nthr)
    end subroutine inmem_initialize

    !> gridding in this process; the outputs are not registered in the project
    subroutine inmem_execute(self, params, build, cline)
        class(rec3D_inmem_strategy), intent(inout) :: self
        type(parameters),            intent(inout) :: params
        type(builder),               intent(inout) :: build
        class(cmdline),              intent(inout) :: cline
        type(rec3D_request) :: request
        call sample_rows(params, build, [params%fromp,params%top], request%pinds)
        request%backend  = REC3D_BACKEND_GRIDDING
        request%weights  = merge(REC3D_WEIGHTS_SET, REC3D_WEIGHTS_HARD, params%l_m_estimator_flex)
        request%dispatch = REC3D_DISPATCH_INPROC
        request%register = .false.
        call self%service%execute(params, build, cline, request)
    end subroutine inmem_execute

    !> kernel PCG in this process; the outputs are registered when the run has its own directory
    subroutine pcg_inmem_execute(self, params, build, cline)
        class(rec3D_pcg_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        type(rec3D_request) :: request
        call sample_rows(params, build, [params%fromp,params%top], request%pinds)
        request%backend  = REC3D_BACKEND_PCG
        request%weights  = merge(REC3D_WEIGHTS_SET, REC3D_WEIGHTS_HARD, params%l_m_estimator_flex)
        request%dispatch = REC3D_DISPATCH_INPROC
        request%register = params%mkdir.eq.'yes'
        call self%service%execute(params, build, cline, request)
    end subroutine pcg_inmem_execute

    subroutine inmem_cleanup(self, params, build, cline)
        use simple_qsys_funs, only: qsys_declare_part_finished
        class(rec3D_inmem_strategy), intent(inout) :: self
        type(parameters),            intent(in)    :: params
        type(builder),               intent(inout) :: build
        class(cmdline),              intent(inout) :: cline
        if( allocated(self%service) )then
            call self%service%kill
            deallocate(self%service)
        endif
        call build%esig%kill
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
        call qsys_declare_part_finished(params, string('simple_rec3D_strategy :: exec_rec3D'))
    end subroutine inmem_cleanup

    ! =====================================================================
    ! DISTRIBUTED-MEMORY IMPLEMENTATION
    ! =====================================================================

    subroutine distr_initialize(self, params, build, cline)
        use simple_exec_helpers, only: set_master_num_threads
        class(rec3D_distr_strategy), intent(inout) :: self
        type(parameters),            intent(inout) :: params
        type(builder),               intent(inout) :: build
        class(cmdline),              intent(inout) :: cline
        logical :: fall_over
        ! master thread count
        call set_master_num_threads(self%nthr_master, string('rec3D'))
        ! parse parameters and project
        call build%init_params_and_build_spproj(cline, params)
        call sync_resolved_rec_params(params, cline)
        ! sanity check
        fall_over = .false.
        select case(trim(params%oritype))
            case('ptcl3D')
                fall_over = build%spproj%get_nptcls() == 0
            case('cls3D')
                fall_over = build%spproj%os_out%get_noris() == 0
            case DEFAULT
                THROW_HARD('unsupported ORITYPE')
        end select
        if( fall_over ) THROW_HARD('No images found!')
        ! avoid nested directory structure for jobs
        call cline%set('mkdir', 'no')
        ! Even/odd partitioning
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
        endif
        ! Update eo flags in project
        call build%spproj%write_segment_inside(params%oritype)
        call ensure_canonical_sigma_state(params, build, cline)
        ! the queue rounds: partitions balance the particles with state > 0
        if( .not. allocated(self%service) ) allocate(self%service)
        call self%service%new(params, build, cline, REC3D_DISPATCH_QUEUE, self%nthr_master)
    end subroutine distr_initialize

    !> the parts reconstruct their range of the sampled rows; the outputs are registered when the run
    !! has its own directory
    subroutine distr_execute(self, params, build, cline)
        class(rec3D_distr_strategy), intent(inout) :: self
        type(parameters),            intent(inout) :: params
        type(builder),               intent(inout) :: build
        class(cmdline),              intent(inout) :: cline
        type(rec3D_request) :: request
        call sample_rows(params, build, [1,build%spproj_field%get_noris()], request%pinds)
        request%backend  = rec3D_backend_id(params%rec_backend)
        request%weights  = merge(REC3D_WEIGHTS_SET, REC3D_WEIGHTS_HARD, params%l_m_estimator_flex)
        request%dispatch = REC3D_DISPATCH_QUEUE
        request%register = params%mkdir.eq.'yes'
        call self%service%execute(params, build, cline, request)
    end subroutine distr_execute

    subroutine distr_cleanup(self, params, build, cline)
        use simple_qsys_funs, only: qsys_cleanup
        class(rec3D_distr_strategy), intent(inout) :: self
        type(parameters),            intent(in)    :: params
        type(builder),               intent(inout) :: build
        class(cmdline),              intent(inout) :: cline
        call qsys_cleanup(params)
        if( allocated(self%service) )then
            call self%service%kill
            deallocate(self%service)
        endif
        call build%spproj_field%kill
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
    end subroutine distr_cleanup

    !> The rows a reconstruct3D run inserts: the previous sampling (update_frac) or every particle with
    !! state > 0 (and updatecnt > 0 once any was updated); a row filter, so a part's range of the master's
    !! rows equals the part's own sampling
    subroutine sample_rows(params, build, fromto, pinds)
        type(parameters),     intent(in)    :: params
        type(builder),        intent(inout) :: build
        integer,              intent(in)    :: fromto(2)
        integer, allocatable, intent(inout) :: pinds(:)
        integer, allocatable :: labels(:)
        integer :: nptcls2update, i
        if( params%l_update_frac .and. build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod(fromto, nptcls2update, pinds)
        else
            call build%spproj_field%sample4rec(fromto, nptcls2update, pinds)
        endif
        ! state= with hard labels: only the rows labelled with that state (weighted runs select by weight)
        if( params%l_state_defined .and. .not. params%l_m_estimator_flex )then
            allocate(labels(size(pinds)))
            do i = 1, size(pinds)
                labels(i) = build%spproj_field%get_state(pinds(i))
            enddo
            pinds = pack(pinds, labels == params%state)
            deallocate(labels)
        endif
    end subroutine sample_rows

    !> The canonical sigma2 state at reconstruction entry, through the one rebuild entry; the parts
    !! of a distributed run consume what their master ensured
    subroutine ensure_canonical_sigma_state(params, build, cline)
        use simple_sigma2_bootstrap, only: ensure_sigma2_for_iteration
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        logical :: l_rebuilt
        if( params%cc_objfun /= OBJFUN_EUCLID ) return
        if( cline%defined('part') ) return
        if( trim(params%oritype) /= 'ptcl3D' ) &
            &THROW_HARD('canonical sigma2 reconstruction requires oritype=ptcl3D')
        call ensure_sigma2_for_iteration(cline, params%projfile, params%oritype, params%which_iter, params%box, &
            &params%smpd, 'RECONSTRUCT3D', l_rebuilt)
        if( l_rebuilt ) call build%spproj%read_segment('projinfo', params%projfile)
    end subroutine ensure_canonical_sigma_state

    subroutine sync_resolved_rec_params(params, cline)
        type(parameters), intent(in)    :: params
        class(cmdline),   intent(inout) :: cline
        if( params%box > 0 ) call cline%set('box', params%box)
        if( params%smpd > TINY ) call cline%set('smpd', params%smpd)
        if( params%box_crop > 0 ) call cline%set('box_crop', params%box_crop)
        if( params%smpd_crop > TINY ) call cline%set('smpd_crop', params%smpd_crop)
        if( params%mskdiam > 0. )then
            call cline%set('mskdiam', params%mskdiam)
        else
            call cline%delete('mskdiam')
        endif
    end subroutine sync_resolved_rec_params

end module simple_rec3D_strategy
