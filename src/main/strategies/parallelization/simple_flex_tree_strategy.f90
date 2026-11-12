!@descr: FLEX tree lifecycle and sequential node FLEX/refinement/state splitting
module simple_flex_tree_strategy
use simple_flex_tree, only: flex_tree, flex_tree_node, flex_tree_decision, &
    &TREE_UNASSESSED, TREE_SPLIT, TREE_LEAF, TREE_JUNK
use simple_string,    only: string
use simple_error,     only: simple_exception
implicit none
private
public :: flex_tree_strategy
#include "simple_local_flags.inc"

type :: flex_tree_strategy
    private
    type(flex_tree), allocatable :: tree
    type(string) :: workdir
    logical :: initialized = .false.
    logical :: completed   = .false.
  contains
    procedure :: initialize
    procedure :: execute
    procedure :: finalize_run
    procedure :: cleanup
    procedure :: get_tree
    procedure, private, nopass :: validate_input
    procedure, private :: initialize_root
    procedure, private :: run_workflow
    procedure, private :: run_node_flex
    procedure, private :: refine_node
    procedure, private :: split_node
    procedure, private :: apply_node_decision
    procedure, private :: publish_results
    procedure, private :: stage_directory
    procedure, private, nopass :: stage_command
    procedure, private, nopass :: execute_stage
end type flex_tree_strategy

contains

    ! Initialize parameters once and read the project through its builder owner.
    ! The caller supplies normalized defaults and retains ownership of params/build.
    subroutine initialize( self, params, build, cline )
        use simple_parameters, only: parameters
        use simple_builder,    only: builder
        use simple_cmdline,    only: cmdline
        class(flex_tree_strategy), intent(inout) :: self
        type(parameters),         intent(inout) :: params
        type(builder),            intent(inout) :: build
        class(cmdline),           intent(inout) :: cline
        if( self%initialized ) THROW_HARD('FLEX tree strategy is already initialized; call cleanup first')
        if( .not. cline%defined('projfile') ) THROW_HARD('FLEX tree initialization requires projfile')
        if( cline%defined('part') ) THROW_HARD('A partition worker cannot initialize a FLEX tree')
        call build%init_params_and_build_spproj(cline, params)
        call self%validate_input(params, build)
        allocate(self%tree)
        call self%tree%new(nroots=1)
        call self%initialize_root(params)
        self%initialized = .true.
        self%completed   = .false.
    end subroutine initialize

    ! Inspect the loaded project only; stage-specific model/identity checks belong to the gateway.
    subroutine validate_input( params, build )
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use simple_parameters, only: parameters
        use simple_builder,    only: builder
        use simple_type_defs,  only: PTCL3D_SEG
        type(parameters), intent(in) :: params
        type(builder),    intent(in) :: build
        integer :: nptcls, iptcl
        logical :: has_selection
        if( params%projfile%is_blank() ) THROW_HARD('FLEX tree input project is blank')
        if( params%spproj_iseg /= PTCL3D_SEG ) THROW_HARD('FLEX tree requires oritype=ptcl3D')
        if( params%box <= 0 ) THROW_HARD('FLEX tree requires a positive project box size')
        if( .not. ieee_is_finite(params%smpd) ) THROW_HARD('FLEX tree requires finite project sampling')
        if( params%smpd <= 0. ) THROW_HARD('FLEX tree requires positive project sampling')
        nptcls = build%spproj%os_ptcl3D%get_noris()
        if( nptcls == 0 ) THROW_HARD('FLEX tree input has no ptcl3D records')
        if( build%spproj%os_ptcl2D%get_noris() /= nptcls ) THROW_HARD('FLEX tree particle fields have different sizes')
        has_selection = .false.
        do iptcl = 1, nptcls
            if( build%spproj%os_ptcl2D%get_state(iptcl) <= 0 ) cycle
            if( build%spproj%os_ptcl3D%get_state(iptcl) <= 0 ) cycle
            has_selection = .true.
            exit
        enddo
        if( .not. has_selection ) THROW_HARD('FLEX tree input has no selected 3D particles')
    end subroutine validate_input

    ! Root membership stays in the project; no particle remapping or file writes here.
    subroutine initialize_root( self, params )
        use simple_parameters, only: parameters
        use simple_syslib,     only: simple_abspath
        class(flex_tree_strategy), intent(inout) :: self
        type(parameters),         intent(in)    :: params
        type(flex_tree_node) :: root
        root%name         = 'root'
        root%project_file = simple_abspath(params%projfile)
        if( .not. params%vols(1)%is_blank() ) root%map_file = simple_abspath(params%vols(1))
        call self%tree%set_source_project(root%project_file)
        call self%tree%set_node(1, root)
    end subroutine initialize_root

    ! Whole-program command objects are supplied by the caller to avoid a commander/module cycle.
    subroutine execute( self, params, cline, flex_runner, refine_runner )
        use simple_parameters,     only: parameters
        use simple_cmdline,        only: cmdline
        use simple_commander_base, only: commander_base
        use simple_syslib,         only: simple_getcwd, simple_mkdir, dir_exists
        class(flex_tree_strategy), intent(inout) :: self
        type(parameters),         intent(in)    :: params
        class(cmdline),           intent(in)    :: cline
        class(commander_base),    intent(inout) :: flex_runner, refine_runner
        type(string) :: cwd
        if( .not. self%initialized ) THROW_HARD('FLEX tree strategy is not initialized')
        if( self%completed ) THROW_HARD('FLEX tree strategy has already completed')
        call simple_getcwd(cwd)
        self%workdir = cwd%to_char()//'/flex_tree'
        if( dir_exists(self%workdir) ) THROW_HARD('FLEX tree work directory exists; use a fresh run directory')
        call simple_mkdir(self%workdir)
        call self%run_workflow(params, cline, flex_runner, refine_runner)
        self%completed = .true.
    end subroutine execute

    ! Breadth-first execution; each command completes before its dependent stage begins.
    subroutine run_workflow( self, params, cline, flex_runner, refine_runner )
        use simple_parameters,     only: parameters
        use simple_cmdline,        only: cmdline
        use simple_commander_base, only: commander_base
        class(flex_tree_strategy), intent(inout) :: self
        type(parameters),         intent(in)    :: params
        class(cmdline),           intent(in)    :: cline
        class(commander_base),    intent(inout) :: flex_runner, refine_runner
        type(flex_tree_node) :: node
        type(flex_tree_decision) :: decision
        type(flex_tree_node), allocatable :: children(:)
        integer :: inode
        logical :: ready
        inode = 1
        do while( inode <= self%tree%get_nnodes() )
            node = self%tree%get_node(inode)
            if( allocated(node%decision) )then
                if( node%decision%outcome /= TREE_UNASSESSED )then
                    inode = inode + 1
                    cycle
                endif
            endif
            if( size(self%tree%get_children(inode)) /= 0 ) THROW_HARD('An unresolved FLEX tree node already has children')
            call self%run_node_flex(params, cline, flex_runner, inode, node, ready)
            if( ready )then
                call self%refine_node(params, cline, refine_runner, inode, node)
                call self%split_node(node, decision, children)
            else
                ! A small child retains its inherited refined map and original membership.
                decision = flex_tree_decision()
                decision%outcome = TREE_LEAF
                decision%reason = 'Insufficient particles for the paired FLEX fit (100 per mod-4 half)'
                if( allocated(children) ) deallocate(children)
                allocate(children(0))
            endif
            call self%apply_node_decision(inode, node, decision, children)
            inode = inode + 1
        enddo
    end subroutine run_workflow

    ! Fit one population and reconstruct its candidate state maps using the existing FLEX command.
    subroutine run_node_flex( self, params, cline, flex_runner, inode, node, ready )
        use simple_parameters,     only: parameters
        use simple_cmdline,        only: cmdline
        use simple_commander_base, only: commander_base
        use simple_flex_tree_project_gateway, only: prepare_tree_flex, collect_tree_fit
        class(flex_tree_strategy), intent(inout) :: self
        type(parameters),         intent(in)    :: params
        class(cmdline),           intent(in)    :: cline
        class(commander_base),    intent(inout) :: flex_runner
        integer,                  intent(in)    :: inode
        type(flex_tree_node),     intent(inout) :: node
        logical,                  intent(out)   :: ready
        type(cmdline), allocatable :: flex_cline
        type(string) :: directory, project, mean
        character(len=24), parameter :: keys(*) = [character(len=24) :: &
            'neigs', 'npreimages', 'preimage_auto', 'n_probe_iters', 'nkern', 'state_axis', &
            'state_placement', 'min_neff', 'min_state_frac', 'column_separation', 'rec_backend', &
            'rec_states_backend', 'maxits_pcg', 'rtol', 'lp', 'box_crop', 'smpd_target', 'umap']
        integer :: i
        directory = self%stage_directory(inode, 'flex')
        project = directory%to_char()//'/flex.simple'
        call prepare_tree_flex(node, project, mean, ready)
        if( .not. ready ) return
        allocate(flex_cline)
        call self%stage_command(params, project, 'flex_pca', flex_cline)
        do i = 1, size(keys)
            call flex_cline%copy_arg(cline, trim(keys(i)))
        enddo
        call flex_cline%set('vol1', mean)
        call flex_cline%set('rec_states', 'yes')
        if( .not. params%pcg_mskfile%is_blank() ) call flex_cline%set('pcg_mskfile', params%pcg_mskfile)
        call self%execute_stage(flex_runner, flex_cline, directory)
        call flex_cline%kill
        node%project_file = project
        if( allocated(node%project_states) ) deallocate(node%project_states)
        call collect_tree_fit(node, mean, directory)
    end subroutine run_node_flex

    ! Ordinary hard-label multistate refinement, initialized from the FLEX state maps and poses.
    subroutine refine_node( self, params, cline, refine_runner, inode, node )
        use simple_parameters,     only: parameters
        use simple_cmdline,        only: cmdline
        use simple_commander_base, only: commander_base
        use simple_string_utils,   only: int2str
        use simple_flex_tree_project_gateway, only: prepare_tree_refine
        class(flex_tree_strategy), intent(inout) :: self
        type(parameters),         intent(in)    :: params
        class(cmdline),           intent(in)    :: cline
        class(commander_base),    intent(inout) :: refine_runner
        integer,                  intent(in)    :: inode
        type(flex_tree_node),     intent(inout) :: node
        type(cmdline), allocatable :: refine_cline
        type(string), allocatable :: volumes(:)
        type(string) :: directory, project
        character(len=24), parameter :: keys(*) = [character(len=24) :: &
            'maxits', 'nspace', 'nspace_sub', 'trs', 'lp', 'ml_reg', 'automsk', 'filt_mode', &
            'prob_beta', 'prob_beta_end', 'athres_cont']
        integer :: i
        directory = self%stage_directory(inode, 'refine')
        project = directory%to_char()//'/refine.simple'
        call prepare_tree_refine(node, project, volumes)
        if( size(volumes) == 0 ) return
        allocate(refine_cline)
        call self%stage_command(params, project, 'refine3D', refine_cline)
        call refine_cline%set('refine', 'prob_neigh')
        call refine_cline%set('prob_neigh_mode', 'state')
        call refine_cline%set('multivol_mode', 'input_oris_refine')
        call refine_cline%set('center', 'no')
        call refine_cline%set('trail_rec', 'no')
        call refine_cline%set('m_estimator', 'no')
        call refine_cline%set('objfun', 'euclid')
        call refine_cline%set('maxits', 10)
        do i = 1, size(keys)
            call refine_cline%copy_arg(cline, trim(keys(i)))
        enddo
        call refine_cline%set('nstates', size(volumes))
        do i = 1, size(volumes)
            call refine_cline%set('vol'//int2str(i), volumes(i))
        enddo
        call self%execute_stage(refine_runner, refine_cline, directory)
        call refine_cline%kill
        node%project_file = project
        if( allocated(node%project_states) ) deallocate(node%project_states)
    end subroutine refine_node

    ! Branching is determined by occupied states, with no fixed binary-tree assumption.
    subroutine split_node( self, node, decision, children )
        use simple_flex_tree_project_gateway, only: split_tree_states
        class(flex_tree_strategy), intent(in)    :: self
        type(flex_tree_node),     intent(inout) :: node
        type(flex_tree_decision), intent(out)   :: decision
        type(flex_tree_node), allocatable, intent(out) :: children(:)
        call split_tree_states(node, decision, children)
    end subroutine split_node

    function stage_directory( self, inode, stage ) result( directory )
        use simple_string_utils, only: int2str_pad
        use simple_syslib,       only: simple_mkdir, dir_exists
        class(flex_tree_strategy), intent(in) :: self
        integer,                  intent(in) :: inode
        character(len=*),         intent(in) :: stage
        type(string) :: directory, parent
        parent = self%workdir%to_char()//'/node_'//int2str_pad(inode, 6)
        call simple_mkdir(parent)
        directory = parent%to_char()//'/'//stage
        if( dir_exists(directory) ) THROW_HARD('FLEX tree stage directory already exists')
        call simple_mkdir(directory)
    end function stage_directory

    ! Fresh commands avoid inheriting tree mode, resume flags or a parent's partition range.
    subroutine stage_command( params, project, program, cline )
        use simple_parameters, only: parameters
        use simple_cmdline,    only: cmdline
        type(parameters), intent(in)    :: params
        type(string),     intent(in)    :: project
        character(len=*), intent(in)    :: program
        type(cmdline),    intent(inout) :: cline
        call cline%set('prg', program)
        call cline%set('projfile', project)
        call cline%set('mkdir', 'no')
        call cline%set('oritype', 'ptcl3D')
        if( params%nparts > 1 ) call cline%set('nparts', params%nparts)
        call cline%set('nthr', params%nthr)
        call cline%set('qsys_name', trim(params%qsys_name))
        call cline%set('pgrp', trim(params%pgrp))
        call cline%set('mskdiam', params%mskdiam)
    end subroutine stage_command

    ! Each full command owns its builder and may use shared-memory or distributed workers.
    subroutine execute_stage( runner, cline, directory )
        use simple_commander_base, only: commander_base
        use simple_cmdline,        only: cmdline
        use simple_syslib,         only: simple_getcwd, simple_chdir
        use simple_defs,           only: CWD_GLOB, CWD_GLOB_ORIG, nthr_glob
        use omp_lib,               only: omp_get_max_threads, omp_set_num_threads
        class(commander_base), intent(inout) :: runner
        type(cmdline),         intent(inout) :: cline
        type(string),          intent(in)    :: directory
        type(string) :: previous_cwd
        character(len=:), allocatable :: previous_global, previous_origin
        integer :: nthr, previous_nthr_glob
        call simple_getcwd(previous_cwd)
        if( allocated(CWD_GLOB) ) previous_global = CWD_GLOB
        if( allocated(CWD_GLOB_ORIG) ) previous_origin = CWD_GLOB_ORIG
        nthr = omp_get_max_threads()
        previous_nthr_glob = nthr_glob
        call simple_chdir(directory)
        CWD_GLOB = directory%to_char()
        call runner%execute(cline)
        call simple_chdir(previous_cwd)
        if( allocated(CWD_GLOB) ) deallocate(CWD_GLOB)
        if( allocated(CWD_GLOB_ORIG) ) deallocate(CWD_GLOB_ORIG)
        if( allocated(previous_global) ) CWD_GLOB = previous_global
        if( allocated(previous_origin) ) CWD_GLOB_ORIG = previous_origin
        call omp_set_num_threads(nthr)
        nthr_glob = previous_nthr_glob
    end subroutine execute_stage

    function get_tree( self ) result( tree )
        class(flex_tree_strategy), intent(in) :: self
        type(flex_tree) :: tree
        if( .not. self%initialized ) THROW_HARD('FLEX tree strategy is not initialized')
        tree = self%tree
    end function get_tree

    ! Commit only the assessed outcome; child creation invalidates the stored parent decision.
    subroutine apply_node_decision( self, inode, node, decision, children )
        class(flex_tree_strategy), intent(inout) :: self
        integer,                  intent(in)    :: inode
        type(flex_tree_node),     intent(inout) :: node
        type(flex_tree_decision), intent(in)    :: decision
        type(flex_tree_node), allocatable, intent(in) :: children(:)
        integer :: ichild, child_id
        select case( decision%outcome )
            case(TREE_SPLIT)
                if( .not. allocated(children) ) THROW_HARD('An accepted split needs child payloads')
                if( size(children) < 2 ) THROW_HARD('An accepted split needs at least two children')
                do ichild = 1, size(children)
                    if( children(ichild)%project_file%is_blank() ) THROW_HARD('A split child needs a project reference')
                    if( allocated(children(ichild)%decision) )then
                        if( children(ichild)%decision%outcome /= TREE_UNASSESSED )then
                            THROW_HARD('A new split child must be unassessed')
                        endif
                    endif
                enddo
                do ichild = 1, size(children)
                    child_id = self%tree%add_child(inode)
                    call self%tree%set_node(child_id, children(ichild))
                enddo
            case(TREE_LEAF, TREE_JUNK)
                if( allocated(children) )then
                    if( size(children) /= 0 ) THROW_HARD('A terminal node cannot create children')
                endif
            case default
                THROW_HARD('Split assessment must return SPLIT, LEAF or JUNK')
        end select
        node%decision = decision
        call self%tree%set_node(inode, node)
    end subroutine apply_node_decision

    subroutine finalize_run( self )
        class(flex_tree_strategy), intent(inout) :: self
        if( .not. self%initialized ) THROW_HARD('FLEX tree strategy is not initialized')
        if( .not. self%completed ) THROW_HARD('Cannot publish an incomplete FLEX tree')
        call self%publish_results()
    end subroutine finalize_run

    ! Publication must not report success before the delivery and checkpoint contracts exist.
    subroutine publish_results( self )
        class(flex_tree_strategy), intent(inout) :: self
        THROW_HARD('FLEX tree result publication is not implemented')
    end subroutine publish_results

    ! Release only strategy-owned state; the caller releases its project and builder.
    subroutine cleanup( self )
        class(flex_tree_strategy), intent(inout) :: self
        if( allocated(self%tree) )then
            call self%tree%kill
            deallocate(self%tree)
        endif
        call self%workdir%kill
        self%initialized = .false.
        self%completed   = .false.
    end subroutine cleanup

end module simple_flex_tree_strategy
