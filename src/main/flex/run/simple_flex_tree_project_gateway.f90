!@descr: FLEX tree stage projects, fitted-model references and state-derived children
module simple_flex_tree_project_gateway
use simple_flex_tree,   only: flex_tree_node, flex_tree_model, flex_tree_decision, TREE_SPLIT, TREE_LEAF, TREE_JUNK
use simple_sp_project,  only: sp_project
use simple_string,      only: string
use simple_error,       only: simple_exception
use simple_defs,        only: MAXS
implicit none
private
public :: prepare_tree_flex, prepare_tree_refine, collect_tree_fit, split_tree_states
#include "simple_local_flags.inc"

contains

    ! Read a frozen node project; selection retains the original row namespace.
    subroutine load_node( node, project, selected, states )
        type(flex_tree_node), intent(in)    :: node
        type(sp_project),    intent(inout) :: project
        logical, allocatable, intent(out)  :: selected(:)
        integer, allocatable, intent(out)  :: states(:)
        integer :: i, label, n, populations(MAXS)
        call project%read(node%project_file)
        n = project%os_ptcl3D%get_noris()
        if( project%os_ptcl2D%get_noris() /= n ) THROW_HARD('FLEX tree particle fields differ in size')
        allocate(selected(n), source=.false.)
        populations = 0
        do i = 1, n
            if( project%os_ptcl2D%get_state(i) <= 0 ) cycle
            label = project%os_ptcl3D%get_state(i)
            if( label <= 0 ) cycle
            if( label > MAXS ) THROW_HARD('FLEX tree state label exceeds MAXS')
            if( allocated(node%project_states) )then
                if( size(node%project_states) > 0 )then
                    if( .not. any(node%project_states == label) ) cycle
                endif
            endif
            selected(i) = .true.
            populations(label) = populations(label) + 1
        enddo
        states = pack([(i, i=1,MAXS)], populations > 0)
    end subroutine load_node

    ! Write a new full-row project; never compact rows or mutate the source project.
    subroutine write_stage_project( project, destination, selected, states, single_state )
        type(sp_project), intent(inout) :: project
        type(string),     intent(in)    :: destination
        logical,          intent(in)    :: selected(:), single_state
        integer,          intent(in)    :: states(:)
        integer :: i, label, j
        do i = 1, size(selected)
            label = 0
            if( selected(i) )then
                label = 1
                if( .not. single_state )then
                    do j = 1, size(states)
                        if( project%os_ptcl3D%get_state(i) == states(j) )then
                            label = j
                            exit
                        endif
                    enddo
                endif
            endif
            call project%os_ptcl3D%set_state(i, label)
            call project%os_ptcl2D%set_state(i, merge(1, 0, selected(i)))
            call project%os_ptcl3D%set(i, 'pind', i)
            call project%os_ptcl2D%set(i, 'pind', i)
        enddo
        ! New selection/labels invalidate inherited weights and stage artifact registrations.
        call project%os_out%kill
        call project%projinfo%delete_entry('sigma2_state')
        call project%update_projinfo(destination)
        call project%write(destination, tempfile=.true.)
    end subroutine write_stage_project

    ! The existing paired FLEX fit needs at least 100 rows in each mod-4 fit.
    subroutine prepare_tree_flex( node, destination, mean, ready )
        use simple_flex_pca_stages, only: flex_pca_half_of, FLEX_FIT_A, FLEX_FIT_B
        use simple_syslib,          only: file_exists, simple_abspath
        type(flex_tree_node), intent(in)  :: node
        type(string),        intent(in)  :: destination
        type(string),        intent(out) :: mean
        logical,             intent(out) :: ready
        type(sp_project), allocatable :: project
        logical, allocatable :: selected(:)
        integer, allocatable :: states(:)
        integer :: i, counts(2), box
        real :: smpd
        allocate(project)
        call load_node(node, project, selected, states)
        counts = 0
        do i = 1, size(selected)
            if( selected(i) ) counts(flex_pca_half_of(i)) = counts(flex_pca_half_of(i)) + 1
        enddo
        ready = counts(FLEX_FIT_A) >= 100 .and. counts(FLEX_FIT_B) >= 100
        if( .not. ready )then
            call project%kill
            return
        endif
        mean = node%map_file
        if( mean%is_blank() )then
            if( size(states) /= 1 ) THROW_HARD('A multi-state root needs an explicit consensus mean')
            call project%get_vol('vol', states(1), mean, smpd, box)
        endif
        if( .not. file_exists(mean) ) THROW_HARD('FLEX tree node mean volume is missing')
        mean = simple_abspath(mean)
        call write_stage_project(project, destination, selected, states, .true.)
        call project%kill
    end subroutine prepare_tree_flex

    ! Compact occupied labels for refine3D, returning corresponding map references in that order.
    subroutine prepare_tree_refine( node, destination, volumes )
        use simple_syslib, only: file_exists
        type(flex_tree_node), intent(in) :: node
        type(string),        intent(in) :: destination
        type(string), allocatable, intent(out) :: volumes(:)
        type(sp_project), allocatable :: project
        logical, allocatable :: selected(:)
        integer, allocatable :: states(:), boxes(:)
        real, allocatable :: smpds(:)
        integer :: i
        allocate(project)
        call load_node(node, project, selected, states)
        allocate(volumes(size(states)), boxes(size(states)), smpds(size(states)))
        do i = 1, size(states)
            call project%get_vol('vol', states(i), volumes(i), smpds(i), boxes(i))
            if( .not. file_exists(volumes(i)) ) THROW_HARD('A FLEX state has no map for refinement')
        enddo
        if( size(states) > 0 )then
            call write_stage_project(project, destination, selected, states, .false.)
            do i = 1, size(states)
                call project%add_vol2os_out(volumes(i), smpds(i), i, 'vol', box=boxes(i))
            enddo
            call project%write_segment_inside('out', destination)
        endif
        call project%kill
    end subroutine prepare_tree_refine

    ! Keep the fitting project frozen while refinement writes a different project.
    subroutine collect_tree_fit( node, mean, directory )
        use simple_syslib, only: file_exists
        type(flex_tree_node), intent(inout) :: node
        type(string),        intent(in)    :: mean, directory
        type(sp_project), allocatable :: project
        type(flex_tree_model) :: model
        type(string) :: kind, path
        integer :: i
        logical :: found
        allocate(project)
        call project%read_segment('out', node%project_file)
        model%name = 'initial'
        model%project_file = node%project_file
        model%mean_file = mean
        do i = 1, project%os_out%get_noris()
            call project%os_out%getter(i, 'imgkind', kind)
            if( kind%to_char() /= 'flex_embedding' ) cycle
            call project%os_out%getter(i, 'flex_embedding', model%embedding_file)
        enddo
        if( model%embedding_file%is_blank() ) THROW_HARD('FLEX did not register a node embedding')
        if( .not. file_exists(model%embedding_file) ) THROW_HARD('FLEX node embedding file is missing')
        call project%get_state_weights(path, found)
        if( found ) model%weights_file = path
        path = directory%to_char()//'/flex_pca_state_targets.txt'
        if( file_exists(path) ) model%targets_file = path
        path = directory%to_char()//'/flex_pca_manifest.txt'
        if( file_exists(path) ) model%provenance_file = path
        if( allocated(node%models) ) deallocate(node%models)
        allocate(node%models(1))
        node%models(1) = model
        call project%kill
    end subroutine collect_tree_fit

    ! Each occupied positive state becomes one child; empty labels never create children.
    subroutine split_tree_states( node, decision, children )
        use simple_string_utils, only: int2str
        type(flex_tree_node),     intent(inout) :: node
        type(flex_tree_decision), intent(out)   :: decision
        type(flex_tree_node), allocatable, intent(out) :: children(:)
        type(sp_project), allocatable :: project
        logical, allocatable :: selected(:)
        integer, allocatable :: states(:)
        integer :: i, nstates
        allocate(project)
        call load_node(node, project, selected, states)
        nstates = size(states)
        select case(nstates)
            case(0)
                decision%outcome = TREE_JUNK
                decision%reason = 'No populated states after refinement'
                allocate(children(0))
            case(1)
                decision%outcome = TREE_LEAF
                decision%reason = 'One populated state after refinement'
                node%project_states = states
                call state_artifacts(project, states(1), node)
                allocate(children(0))
            case default
                decision%outcome = TREE_SPLIT
                decision%reason = int2str(nstates)//' populated states after refinement'
                allocate(children(nstates))
                do i = 1, nstates
                    children(i)%name = node%name%to_char()//'_'//int2str(states(i))
                    children(i)%project_file = node%project_file
                    children(i)%project_states = [states(i)]
                    call state_artifacts(project, states(i), children(i))
                enddo
        end select
        call project%kill
    end subroutine split_tree_states

    subroutine state_artifacts( project, state, node )
        use simple_syslib,           only: file_exists
        use simple_fileio,           only: get_fpath
        use simple_refine3D_fnames, only: refine3D_state_halfvol_fname
        type(sp_project),     intent(in)    :: project
        integer,              intent(in)    :: state
        type(flex_tree_node), intent(inout) :: node
        type(string) :: directory, name, even, odd
        real :: smpd
        integer :: box
        call project%get_vol('vol', state, node%map_file, smpd, box)
        if( .not. file_exists(node%map_file) ) THROW_HARD('Refined state map is missing')
        if( project%isthere_in_osout('fsc', state) ) call project%get_fsc(state, node%fsc_file, box)
        directory = get_fpath(node%map_file)
        name = refine3D_state_halfvol_fname(state, 'even', unfil=.true.)
        even = directory%to_char()//name%to_char()
        name = refine3D_state_halfvol_fname(state, 'odd', unfil=.true.)
        odd = directory%to_char()//name%to_char()
        call node%even_file%kill
        call node%odd_file%kill
        if( file_exists(even) .and. file_exists(odd) )then
            node%even_file = even
            node%odd_file = odd
        endif
    end subroutine state_artifacts

end module simple_flex_tree_project_gateway
