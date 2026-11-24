!@descr: FLEX tree lifecycle and sequential node FLEX/refinement/state splitting
module simple_flex_tree_strategy
use simple_flex_tree,      only: flex_tree, flex_tree_node, flex_tree_decision, &
    &TREE_UNASSESSED, TREE_SPLIT, TREE_LEAF, TREE_JUNK
use simple_string,         only: string
use simple_error,          only: simple_exception
use simple_parameters,     only: parameters
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
implicit none
private
public :: flex_tree_strategy, flex_tree_inmem_strategy, flex_tree_distr_strategy, create_flex_tree_strategy
#include "simple_local_flags.inc"

logical, parameter :: ROOT_MASK_GROUPS = .true.
integer, parameter :: PTCLS_PER_PART = 20000, NPARTS_MIN = 4, FLEX_MASTER_NTHR = 32, FLEX_WIDE_PTCLS = 100000, POLL_SECS = 10
character(len=*), parameter :: TREE_STATE = 'tree_state.txt', STEP_DONE = 'TREE_STEP_DONE'

type, abstract :: flex_tree_strategy
    private
    type(flex_tree), allocatable :: tree
    type(flex_tree_node), allocatable :: jparent
    type(string) :: workdir
    class(commander_base), pointer :: flex_runner => null(), refine_runner => null()
    integer :: nrows = 0, part_nthr = 0
    logical :: initialized = .false.
    logical :: completed   = .false.
  contains
    procedure(init_interface),     deferred :: initialize
    procedure(exec_interface),     deferred :: execute
    procedure :: finalize_run => tree_finalize_run
    procedure :: cleanup      => tree_cleanup
    procedure :: end_message  => tree_end_message
    procedure :: set_runners
    procedure :: get_tree
    procedure, private :: tree_init
    procedure, private, nopass :: validate_input
    procedure, private :: initialize_root
    procedure, private :: run_node
    procedure, private :: write_state
    procedure, private :: read_state
    procedure, private :: node_dir
    procedure, private :: run_node_flex
    procedure, private :: refine_node
    procedure, private :: split_node
    procedure, private :: apply_node_decision
    procedure, private :: publish_results
    procedure, private :: stage_directory
    procedure, private :: stage_command
    procedure, private, nopass :: execute_stage
end type flex_tree_strategy

! shared memory, or a distributed worker running the node given by part
type, extends(flex_tree_strategy) :: flex_tree_inmem_strategy
  contains
    procedure :: initialize   => inmem_initialize
    procedure :: execute      => inmem_execute
end type flex_tree_inmem_strategy

type, extends(flex_tree_strategy) :: flex_tree_distr_strategy
  contains
    procedure :: initialize   => distr_initialize
    procedure :: execute      => distr_execute
end type flex_tree_distr_strategy

abstract interface
    subroutine init_interface( self, params, build, cline )
        import :: flex_tree_strategy, parameters, builder, cmdline
        class(flex_tree_strategy), intent(inout) :: self
        type(parameters),          intent(inout) :: params
        type(builder),             intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface( self, params, cline )
        import :: flex_tree_strategy, parameters, cmdline
        class(flex_tree_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
        class(cmdline),            intent(in)    :: cline
    end subroutine exec_interface

end interface

contains

    function create_flex_tree_strategy( cline ) result( strategy )
        class(cmdline), intent(in) :: cline
        class(flex_tree_strategy), allocatable :: strategy
        if( cline%defined('nparts') .and. .not. cline%defined('part') )then
            allocate(flex_tree_distr_strategy :: strategy)
        else
            allocate(flex_tree_inmem_strategy :: strategy)
        endif
    end function create_flex_tree_strategy

    subroutine set_runners( self, flex_runner, refine_runner )
        class(flex_tree_strategy),     intent(inout) :: self
        class(commander_base), target, intent(inout) :: flex_runner, refine_runner
        self%flex_runner   => flex_runner
        self%refine_runner => refine_runner
    end subroutine set_runners


    ! Initialize parameters once and read the project through its builder owner.
    ! The caller supplies normalized defaults and retains ownership of params/build.
    subroutine tree_init( self, params, build, cline )
        class(flex_tree_strategy), intent(inout) :: self
        type(parameters),         intent(inout) :: params
        type(builder),            intent(inout) :: build
        class(cmdline),           intent(inout) :: cline
        if( self%initialized ) THROW_HARD('FLEX tree strategy is already initialized; call cleanup first')
        if( .not. cline%defined('projfile') ) THROW_HARD('FLEX tree initialization requires projfile')
        call build%init_params_and_build_spproj(cline, params)
        call self%validate_input(params, build)
        allocate(self%tree)
        call self%tree%new(nroots=1)
        call self%initialize_root(params)
        self%initialized = .true.
        self%completed   = .false.
    end subroutine tree_init

    subroutine inmem_initialize( self, params, build, cline )
        class(flex_tree_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        call self%tree_init(params, build, cline)
    end subroutine inmem_initialize

    subroutine distr_initialize( self, params, build, cline )
        class(flex_tree_distr_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        call self%tree_init(params, build, cline)
    end subroutine distr_initialize

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

    subroutine distr_execute( self, params, cline )
        use simple_syslib,         only: simple_getcwd, simple_mkdir, dir_exists, file_exists
        class(flex_tree_distr_strategy), intent(inout) :: self
        type(parameters),               intent(in)    :: params
        class(cmdline),                 intent(in)    :: cline
        type(string) :: cwd
        if( .not. self%initialized ) THROW_HARD('FLEX tree strategy is not initialized')
        if( self%completed ) THROW_HARD('FLEX tree strategy has already completed')
        call simple_getcwd(cwd)
        self%workdir = cwd%to_char()//'/flex_tree'
        if( dir_exists(self%workdir) )then
            if( .not. file_exists(self%workdir//'/'//TREE_STATE) ) THROW_HARD('FLEX tree work directory exists without a tree state')
            call self%read_state
        else
            call simple_mkdir(self%workdir)
            call self%write_state
        endif
        call run_queue(self, params, cline)
        self%completed = .true.
    end subroutine distr_execute

    subroutine run_queue( self, params, cline )
        use simple_qsys_env,       only: qsys_env
        use simple_qsys_async_job, only: qsys_async_job, ASYNC_JOB_DONE, ASYNC_JOB_FAILED
        use simple_string_utils,   only: int2str
        use simple_defs,           only: logfhandle
        use simple_flex_tree_project_gateway, only: tree_node_rows, write_tree_record, read_tree_record, tree_value, &
            &get_tree_node
        class(flex_tree_distr_strategy), intent(inout) :: self
        type(parameters),               intent(in)    :: params
        class(cmdline),                 intent(in)    :: cline
        type(qsys_env) :: qstd, qwide
        type(qsys_async_job), allocatable :: jobs(:)
        type(flex_tree_node), allocatable :: children(:), recs(:)
        type(string),         allocatable :: lines(:), head(:), tags(:)
        integer,              allocatable :: jnode(:)
        type(flex_tree_node)     :: node
        type(flex_tree_decision) :: decision
        type(cmdline) :: jcline
        type(string)  :: dir
        integer :: inode, slot, nrows, nwide, i, nch, iparent
        nwide = max(params%nthr, FLEX_MASTER_NTHR)
        call qstd%new(params, 1, exec_bin=string('simple_exec'), qsys_nthr=params%nthr)
        call qwide%new(params, 1, exec_bin=string('simple_exec'), qsys_nthr=nwide)
        allocate(jobs(max(1, params%nodepar)), jnode(max(1, params%nodepar)))
        jnode = 0
        do
            do slot = 1, size(jnode)
                if( jnode(slot) == 0 ) cycle
                select case( jobs(slot)%status() )
                    case(ASYNC_JOB_FAILED)
                        dir = jobs(slot)%get_log()
                        THROW_HARD('FLEX tree node job failed; see '//dir%to_char())
                    case(ASYNC_JOB_DONE)
                        inode = jnode(slot)
                        call read_tree_record(self%node_dir(inode)//'/node_out.txt', lines)
                        call get_tree_node(lines, 'node', node)
                        decision = node%decision
                        deallocate(node%decision)
                        nch = str2int(tree_value(lines, 'job', 'nchildren'))
                        if( allocated(children) ) deallocate(children)
                        allocate(children(nch))
                        do i = 1, nch
                            call get_tree_node(lines, 'child'//int2str(i), children(i))
                        enddo
                        call self%apply_node_decision(inode, node, decision, children)
                        jnode(slot) = 0
                        call self%write_state
                        write(logfhandle,'(A,I0,A)') '>>> FLEX TREE node ', inode, ': '//decision%reason%to_char()
                end select
            enddo
            do inode = 1, self%tree%get_nnodes()
                if( count(jnode > 0) == size(jnode) ) exit
                if( any(jnode == inode) ) cycle
                node = self%tree%get_node(inode)
                if( allocated(node%decision) )then
                    if( node%decision%outcome /= TREE_UNASSESSED ) cycle
                endif
                if( params%maxdepth > 0 .and. depth(inode) > params%maxdepth )then
                    decision = flex_tree_decision()
                    decision%outcome = TREE_LEAF
                    decision%reason  = 'Below maxdepth'
                    if( allocated(children) ) deallocate(children)
                    allocate(children(0))
                    call self%apply_node_decision(inode, node, decision, children)
                    call self%write_state
                    cycle
                endif
                slot  = findloc(jnode, 0, dim=1)
                nrows = count(tree_node_rows(node))
                dir   = self%node_dir(inode)
                head  = [string('job id '//int2str(inode)), string('job workdir '//self%workdir%to_char()), &
                    &string('job part_nthr '//int2str(params%nthr))]
                tags  = [string('node')]
                recs  = [node]
                iparent = self%tree%get_parent(inode)
                if( iparent > 0 )then
                    tags = [tags, string('parent')]
                    recs = [recs, self%tree%get_node(iparent)]
                endif
                call write_tree_record(dir//'/node_in.txt', head, tags, recs)
                jcline = cline
                call jcline%set('part', inode)
                call jcline%set('mkdir', 'no')
                if( nrows >= FLEX_WIDE_PTCLS )then
                    call jcline%set('nthr', nwide)
                    call jobs(slot)%start(qwide, jcline, dir, 'node')
                else
                    call jobs(slot)%start(qstd, jcline, dir, 'node')
                endif
                jnode(slot) = inode
                write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX TREE node ', inode, ' submitted (', nrows, ' particles)'
            enddo
            if( all(jnode == 0) ) exit
            call sleep(POLL_SECS)
        enddo
        call qstd%kill
        call qwide%kill

    contains

        integer function depth( i )
            integer, intent(in) :: i
            integer :: j
            depth = 0
            j = self%tree%get_parent(i)
            do while( j > 0 )
                depth = depth + 1
                j = self%tree%get_parent(j)
            enddo
        end function depth

        integer function str2int( v )
            use simple_string_utils, only: s2i => str2int
            type(string), intent(in) :: v
            str2int = 0
            if( .not. v%is_blank() ) str2int = s2i(v%to_char())
        end function str2int

    end subroutine run_queue

    subroutine inmem_execute( self, params, cline )
        use simple_syslib,       only: simple_getcwd, simple_mkdir, dir_exists
        use simple_string_utils, only: int2str, str2int
        use simple_flex_tree_project_gateway, only: tree_node_rows, write_tree_record, read_tree_record, tree_value, &
            &get_tree_node
        class(flex_tree_inmem_strategy), intent(inout) :: self
        type(parameters),               intent(in)    :: params
        class(cmdline),                 intent(in)    :: cline
        type(flex_tree_node), allocatable :: children(:)
        type(string),         allocatable :: lines(:), head(:), tags(:)
        type(flex_tree_node)     :: node
        type(flex_tree_decision) :: decision
        type(string) :: cwd, v
        integer :: inode, i
        if( .not. self%initialized ) THROW_HARD('FLEX tree strategy is not initialized')
        if( self%completed ) THROW_HARD('FLEX tree strategy has already completed')
        call simple_getcwd(cwd)
        if( params%l_distr_worker )then
            call read_tree_record(cwd//'/node_in.txt', lines)
            inode = params%part
            self%workdir = tree_value(lines, 'job', 'workdir')
            v = tree_value(lines, 'job', 'part_nthr')
            self%part_nthr = str2int(v%to_char())
            call get_tree_node(lines, 'node', node)
            v = tree_value(lines, 'parent', 'project_file')
            if( .not. v%is_blank() )then
                allocate(self%jparent)
                call get_tree_node(lines, 'parent', self%jparent)
            endif
            self%nrows = count(tree_node_rows(node))
            call self%run_node(params, cline, inode, node, decision, children)
            allocate(node%decision, source=decision)
            head = [string('job nchildren '//int2str(size(children)))]
            tags = [string('node')]
            do i = 1, size(children)
                tags = [tags, string('child'//int2str(i))]
            enddo
            call write_tree_record(self%node_dir(inode)//'/node_out.txt', head, tags, [node, children])
        else
            self%workdir = cwd%to_char()//'/flex_tree'
            if( dir_exists(self%workdir) ) THROW_HARD('FLEX tree work directory exists; use a fresh run directory')
            call simple_mkdir(self%workdir)
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
                if( self%tree%get_parent(inode) > 0 )then
                    if( allocated(self%jparent) ) deallocate(self%jparent)
                    allocate(self%jparent, source=self%tree%get_node(self%tree%get_parent(inode)))
                endif
                self%nrows = count(tree_node_rows(node))
                call self%run_node(params, cline, inode, node, decision, children)
                call self%apply_node_decision(inode, node, decision, children)
                inode = inode + 1
            enddo
        endif
        self%completed = .true.
    end subroutine inmem_execute

    subroutine run_node( self, params, cline, inode, node, decision, children )
        class(flex_tree_strategy), intent(inout) :: self
        type(parameters),         intent(in)    :: params
        class(cmdline),           intent(in)    :: cline
        integer,                  intent(in)    :: inode
        type(flex_tree_node),     intent(inout) :: node
        type(flex_tree_decision), intent(out)   :: decision
        type(flex_tree_node), allocatable, intent(out) :: children(:)
        logical :: ready
        if( .not. associated(self%flex_runner) ) THROW_HARD('FLEX tree runners are not set')
        call self%run_node_flex(params, cline, self%flex_runner, inode, node, ready)
        if( ready )then
            call self%refine_node(params, cline, self%refine_runner, inode, node)
            call self%split_node(node, decision, children)
        else
            ! A small child retains its inherited refined map and original membership.
            decision = flex_tree_decision()
            decision%outcome = TREE_LEAF
            decision%reason = 'Insufficient particles for the paired FLEX fit (100 per mod-4 half)'
            allocate(children(0))
        endif
    end subroutine run_node

    function node_dir( self, inode ) result( dir )
        use simple_string_utils, only: int2str_pad
        use simple_syslib,       only: simple_mkdir
        class(flex_tree_strategy), intent(in) :: self
        integer,                  intent(in) :: inode
        type(string) :: dir
        dir = self%workdir%to_char()//'/node_'//int2str_pad(inode, 6)
        call simple_mkdir(dir)
    end function node_dir

    subroutine write_state( self )
        use simple_string_utils, only: int2str
        use simple_flex_tree_project_gateway, only: write_tree_record
        class(flex_tree_strategy), intent(inout) :: self
        type(string),         allocatable :: head(:), tags(:)
        type(flex_tree_node), allocatable :: nodes(:)
        integer :: i, n
        n = self%tree%get_nnodes()
        allocate(head(n), tags(n), nodes(n))
        do i = 1, n
            head(i)  = 'tree parent '//int2str(i)//' '//int2str(self%tree%get_parent(i))
            tags(i)  = 'n'//int2str(i)
            nodes(i) = self%tree%get_node(i)
        enddo
        call write_tree_record(self%workdir//'/'//TREE_STATE, head, tags, nodes)
    end subroutine write_state

    subroutine read_state( self )
        use simple_string_utils, only: int2str, str2int
        use simple_flex_tree_project_gateway, only: read_tree_record, get_tree_node
        class(flex_tree_strategy), intent(inout) :: self
        type(string),         allocatable :: lines(:)
        type(flex_tree_node), allocatable :: nodes(:)
        integer, allocatable :: parents(:)
        character(len=:), allocatable :: line
        integer :: i, n, id, ip, k
        call read_tree_record(self%workdir//'/'//TREE_STATE, lines)
        n = count([(index(lines(i)%to_char(), 'tree parent ') == 1, i=1,size(lines))])
        allocate(parents(n), nodes(n))
        do i = 1, size(lines)
            line = lines(i)%to_char()
            if( index(line, 'tree parent ') /= 1 ) cycle
            read(line(13:),*) id, ip
            parents(id) = ip
        enddo
        do i = 1, n
            call get_tree_node(lines, 'n'//int2str(i), nodes(i))
        enddo
        do i = 1, n
            if( i > 1 )then
                k = self%tree%add_child(parents(i))
                if( k /= i ) THROW_HARD('FLEX tree state does not rebuild in node order')
            endif
            if( allocated(nodes(i)%decision) )then
                call self%tree%set_node(i, strip(nodes(i)))
            else
                call self%tree%set_node(i, nodes(i))
            endif
        enddo
        do i = n, 1, -1
            call self%tree%set_node(i, nodes(i))
        enddo

    contains

        function strip( nd ) result( out )
            type(flex_tree_node), intent(in) :: nd
            type(flex_tree_node) :: out
            out = nd
            deallocate(out%decision)
        end function strip

    end subroutine read_state

    ! Fit one population and reconstruct its candidate state maps using the existing FLEX command.
    subroutine run_node_flex( self, params, cline, flex_runner, inode, node, ready )
        use simple_parameters,     only: parameters
        use simple_cmdline,        only: cmdline
        use simple_commander_base, only: commander_base
        use simple_fileio,         only: get_fpath
        use simple_defs,           only: logfhandle
        use simple_commanders_rec, only: commander_rec3D
        use simple_flex_tree_project_gateway, only: prepare_tree_flex, collect_tree_fit, write_tree_mask, write_tree_groups, &
            &tree_node_rows
        class(flex_tree_strategy), intent(inout) :: self
        type(parameters),         intent(in)    :: params
        class(cmdline),           intent(in)    :: cline
        class(commander_base),    intent(inout) :: flex_runner
        integer,                  intent(in)    :: inode
        type(flex_tree_node),     intent(inout) :: node
        logical,                  intent(out)   :: ready
        type(commander_rec3D)      :: xrec3D
        type(cmdline), allocatable :: rec_cline
        type(flex_tree_node)       :: parent
        type(string) :: directory, project, mean, mask, gdir, gproj
        integer      :: ngroups, nmaps
        logical      :: l_done, l_gdone
        character(len=24), parameter :: keys(*) = [character(len=24) :: &
            'neigs', 'npreimages', 'preimage_auto', 'n_probe_iters', 'nkern', 'state_axis', &
            'state_placement', 'min_neff', 'min_state_frac', 'column_separation', 'rec_backend', &
            'rec_states_backend', 'maxits_pcg', 'rtol', 'lp', 'box_crop', 'smpd_target', 'umap']
        directory = self%stage_directory(inode, 'flex', l_done)
        project = directory%to_char()//'/flex.simple'
        call prepare_tree_flex(node, project, mean, ready, l_write=.not. l_done)
        if( .not. ready ) return
        if( .not. params%pcg_mskfile%is_blank() ) mask = params%pcg_mskfile
        if( allocated(self%jparent) )then
            ! a child's mask: the parent's flex states its rows occupy (the parent's stA0 maps once stA0 is ported)
            parent = self%jparent
            mask   = get_fpath(directory)//'variance_mask.mrc'
            call write_tree_mask(parent%models(1)%project_file, params%mskdiam, mask, nmaps, tree_node_rows(node))
            if( nmaps < 2 )then
                mask = parent%mask_file
                write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX TREE node ', inode, ' occupies ', nmaps, &
                    &' parent flex state(s): it inherits its parent''s mask'
            endif
            node%mask_file = mask
        endif
        if( .not. l_done ) call fit(directory, project)
        if( .not. allocated(self%jparent) )then
            mask = get_fpath(directory)//'variance_mask.mrc'
            if( ROOT_MASK_GROUPS )then
                gdir  = self%stage_directory(inode, 'groups', l_gdone)
                gproj = gdir%to_char()//'/groups.simple'
                if( .not. l_gdone )then
                    call write_tree_groups(directory, project, gproj, ngroups)
                    allocate(rec_cline)
                    call self%stage_command(params, gproj, 'reconstruct3D', rec_cline)
                    call rec_cline%set('nstates', ngroups)
                    call rec_cline%set('ml_reg', 'no')
                    call rec_cline%set('postprocess', 'no')
                    call self%execute_stage(xrec3D, rec_cline, gdir)
                    call rec_cline%kill
                endif
                call write_tree_mask(gproj, params%mskdiam, mask, nmaps)
            else
                call write_tree_mask(project, params%mskdiam, mask, nmaps)
            endif
            if( nmaps < 2 ) THROW_HARD('The root variance mask needs at least two populated state maps')
            directory = self%stage_directory(inode, 'fflex', l_done)
            project = directory%to_char()//'/fflex.simple'
            call prepare_tree_flex(node, project, mean, ready, l_write=.not. l_done)
            if( .not. l_done ) call fit(directory, project)
            node%mask_file = mask
        endif
        node%project_file = project
        if( allocated(node%project_states) ) deallocate(node%project_states)
        call collect_tree_fit(node, mean, directory)

    contains

        subroutine fit( stage_dir, stage_project )
            type(string), intent(in) :: stage_dir, stage_project
            type(cmdline), allocatable :: flex_cline
            integer :: i
            allocate(flex_cline)
            call self%stage_command(params, stage_project, 'flex_pca', flex_cline)
            do i = 1, size(keys)
                call flex_cline%copy_arg(cline, trim(keys(i)))
            enddo
            call flex_cline%set('vol1', mean)
            call flex_cline%set('rec_states', 'yes')
            if( self%nrows >= FLEX_WIDE_PTCLS ) call flex_cline%set('nthr', params%nthr)
            if( .not. mask%is_blank() ) call flex_cline%set('pcg_mskfile', mask)
            call self%execute_stage(flex_runner, flex_cline, stage_dir)
            call flex_cline%kill
        end subroutine fit

    end subroutine run_node_flex

    ! Ordinary hard-label multistate refinement, initialized from the FLEX state maps and poses.
    ! the node's whole refinement procedure, swapped by editing this routine: stage A 15 A box 64, then stage B 10 A box 88
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
        character(len=24), parameter :: keys(*) = [character(len=24) :: &
            'nspace_sub', 'trs', 'ml_reg', 'automsk', 'filt_mode', 'athres_cont']
        type(string) :: mask
        logical :: l_done
        mask = node%mask_file
        if( mask%is_blank() .and. .not. params%pcg_mskfile%is_blank() ) mask = params%pcg_mskfile
        call refine_stage('stageA', 15., 64, 2500, 3)
        if( l_done ) return
        call refine_stage('stageB', 10., 88, 5000, 3)

    contains

        subroutine refine_stage( stage, lp, box_crop, nspace, maxits )
            character(len=*), intent(in) :: stage
            real,             intent(in) :: lp
            integer,          intent(in) :: box_crop, nspace, maxits
            type(cmdline), allocatable :: refine_cline
            type(string),  allocatable :: volumes(:)
            type(string) :: directory, project
            integer :: i
            logical :: l_sdone
            directory = self%stage_directory(inode, stage, l_sdone)
            project = directory%to_char()//'/refine.simple'
            call prepare_tree_refine(node, project, volumes, l_write=.not. l_sdone)
            l_done = size(volumes) == 0
            if( l_done ) return
            if( l_sdone )then
                node%project_file = project
                if( allocated(node%project_states) ) deallocate(node%project_states)
                return
            endif
            allocate(refine_cline)
            call self%stage_command(params, project, 'refine3D', refine_cline)
            call refine_cline%set('refine', 'prob_neigh')
            call refine_cline%set('prob_neigh_mode', 'geom')
            call refine_cline%set('multivol_mode', 'input_oris_refine')
            call refine_cline%set('prob_beta', 10.)
            call refine_cline%set('center', 'no')
            call refine_cline%set('trail_rec', 'no')
            call refine_cline%set('m_estimator', 'no')
            call refine_cline%set('objfun', 'euclid')
            call refine_cline%set('rec_backend', 'pcg')
            ! the tree's reference runs also focused the reconstruction (pcg_focus), which master lacks
            if( .not. mask%is_blank() ) call refine_cline%set('pcg_mskfile', mask)
            do i = 1, size(keys)
                call refine_cline%copy_arg(cline, trim(keys(i)))
            enddo
            call refine_cline%set('lp', lp)
            call refine_cline%set('box_crop', box_crop)
            call refine_cline%set('nspace', nspace)
            call refine_cline%set('maxits', maxits)
            call refine_cline%set('nstates', size(volumes))
            do i = 1, size(volumes)
                call refine_cline%set('vol'//int2str(i), volumes(i))
            enddo
            call self%execute_stage(refine_runner, refine_cline, directory)
            call refine_cline%kill
            node%project_file = project
            if( allocated(node%project_states) ) deallocate(node%project_states)
        end subroutine refine_stage

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

    ! a finished stage (its done marker) is reused; an unfinished one is moved aside and started afresh
    function stage_directory( self, inode, stage, l_done ) result( directory )
        use simple_string_utils, only: int2str
        use simple_syslib,       only: simple_mkdir, dir_exists, file_exists, simple_rename
        class(flex_tree_strategy), intent(in)  :: self
        integer,                  intent(in)  :: inode
        character(len=*),         intent(in)  :: stage
        logical,                  intent(out) :: l_done
        type(string) :: directory
        integer :: k
        directory = self%node_dir(inode)//'/'//stage
        l_done = file_exists(directory//'/'//STEP_DONE)
        if( l_done ) return
        if( dir_exists(directory) )then
            k = 1
            do while( dir_exists(directory//'_aborted'//int2str(k)) )
                k = k + 1
            enddo
            call simple_rename(directory, directory//'_aborted'//int2str(k))
        endif
        call simple_mkdir(directory)
    end function stage_directory

    ! Fresh commands avoid inheriting tree mode, resume flags or a parent's partition range.
    subroutine stage_command( self, params, project, program, cline )
        use simple_parameters, only: parameters
        use simple_cmdline,    only: cmdline
        class(flex_tree_strategy), intent(in) :: self
        type(parameters), intent(in)    :: params
        type(string),     intent(in)    :: project
        character(len=*), intent(in)    :: program
        type(cmdline),    intent(inout) :: cline
        call cline%set('prg', program)
        call cline%set('projfile', project)
        call cline%set('mkdir', 'no')
        call cline%set('oritype', 'ptcl3D')
        if( params%nparts > 1 ) call cline%set('nparts', min(params%nparts, max(NPARTS_MIN, self%nrows / PTCLS_PER_PART)))
        call cline%set('nthr', merge(self%part_nthr, params%nthr, self%part_nthr > 0))
        call cline%set('qsys_name', trim(params%qsys_name))
        call cline%set('pgrp', trim(params%pgrp))
        call cline%set('mskdiam', params%mskdiam)
    end subroutine stage_command

    ! Each full command owns its builder and may use shared-memory or distributed workers.
    subroutine execute_stage( runner, cline, directory )
        use simple_commander_base, only: commander_base
        use simple_cmdline,        only: cmdline
        use simple_syslib,         only: simple_getcwd, simple_chdir, simple_touch
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
        call simple_touch(directory//'/'//STEP_DONE)
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
        if( allocated(node%decision) ) deallocate(node%decision)
        allocate(node%decision, source=decision)
        call self%tree%set_node(inode, node)
    end subroutine apply_node_decision

    subroutine tree_finalize_run( self )
        class(flex_tree_strategy), intent(inout) :: self
        if( .not. self%initialized ) THROW_HARD('FLEX tree strategy is not initialized')
        if( .not. self%completed ) THROW_HARD('Cannot publish an incomplete FLEX tree')
        call self%publish_results()
    end subroutine tree_finalize_run

    ! Publication must not report success before the delivery and checkpoint contracts exist.
    subroutine publish_results( self )
        class(flex_tree_strategy), intent(inout) :: self
        THROW_HARD('FLEX tree result publication is not implemented')
    end subroutine publish_results

    ! Release only strategy-owned state; the caller releases its project and builder.
    subroutine tree_cleanup( self )
        class(flex_tree_strategy), intent(inout) :: self
        if( allocated(self%tree) )then
            call self%tree%kill
            deallocate(self%tree)
        endif
        call self%workdir%kill
        self%initialized = .false.
        self%completed   = .false.
    end subroutine tree_cleanup

    function tree_end_message( self ) result( msg )
        class(flex_tree_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_FLEX_PCA TREE NORMAL STOP ****'
    end function tree_end_message

end module simple_flex_tree_strategy
