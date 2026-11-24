!@descr: FLEX tree stage projects, fitted-model references and state-derived children
module simple_flex_tree_project_gateway
use simple_flex_tree,   only: flex_tree_node, flex_tree_model, flex_tree_decision, TREE_SPLIT, TREE_LEAF, TREE_JUNK
use simple_sp_project,  only: sp_project
use simple_string,      only: string
use simple_error,       only: simple_exception
use simple_defs,        only: MAXS
implicit none
private
public :: prepare_tree_flex, prepare_tree_refine, collect_tree_fit, split_tree_states, write_tree_mask, write_tree_groups, tree_node_rows
public :: write_tree_record, read_tree_record, tree_value, get_tree_node
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
        call project%os_ptcl3D%clean_entry('updatecnt', 'sampled')
        call project%os_ptcl2D%clean_entry('updatecnt', 'sampled')
        call project%os_out%kill
        call project%projinfo%delete_entry('sigma2_state')
        call project%update_projinfo(destination)
        call project%write(destination, tempfile=.true.)
    end subroutine write_stage_project

    ! The existing paired FLEX fit needs at least 100 rows in each mod-4 fit.
    subroutine prepare_tree_flex( node, destination, mean, ready, l_write )
        use simple_flex_pca_stages, only: flex_pca_half_of, FLEX_FIT_A, FLEX_FIT_B
        use simple_syslib,          only: file_exists, simple_abspath
        type(flex_tree_node), intent(in)  :: node
        type(string),        intent(in)  :: destination
        type(string),        intent(out) :: mean
        logical,             intent(out) :: ready
        logical, optional,   intent(in)  :: l_write
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
        if( .not. present(l_write) )then
            call write_stage_project(project, destination, selected, states, .true.)
        else if( l_write )then
            call write_stage_project(project, destination, selected, states, .true.)
        endif
        call project%kill
    end subroutine prepare_tree_flex

    ! Compact occupied labels for refine3D, returning corresponding map references in that order.
    subroutine prepare_tree_refine( node, destination, volumes, l_write )
        use simple_syslib, only: file_exists
        type(flex_tree_node), intent(in) :: node
        type(string),        intent(in) :: destination
        type(string), allocatable, intent(out) :: volumes(:)
        logical, optional,   intent(in) :: l_write
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
        if( present(l_write) )then
            if( .not. l_write ) states = [integer ::]
        endif
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

    subroutine write_tree_mask( project_file, mskdiam, fname, nmaps, rows )
        use simple_defs,             only: dp
        use simple_syslib,           only: file_exists
        use simple_fileio,           only: get_fpath
        use simple_refine3D_fnames, only: refine3D_state_halfvol_fname
        use simple_variance_mask,    only: write_variance_mask
        type(string),      intent(in)  :: project_file, fname
        real,              intent(in)  :: mskdiam
        integer,           intent(out) :: nmaps
        logical, optional, intent(in)  :: rows(:)
        type(sp_project), allocatable :: project
        type(string),     allocatable :: evens(:), odds(:)
        integer,          allocatable :: labels(:)
        real(dp),         allocatable :: w(:)
        type(string) :: directory, even, odd
        integer :: i, state
        allocate(project)
        call project%read_segment('ptcl3D', project_file)
        labels = [(project%os_ptcl3D%get_state(i), i=1,project%os_ptcl3D%get_noris())]
        call project%kill
        if( present(rows) ) where( .not. rows ) labels = 0
        directory = get_fpath(project_file)
        allocate(evens(0), odds(0), w(0))
        do state = 1, maxval(labels)
            if( count(labels == state) == 0 ) cycle
            even = directory//refine3D_state_halfvol_fname(state, 'even', unfil=.true.)
            odd  = directory//refine3D_state_halfvol_fname(state, 'odd',  unfil=.true.)
            if( .not. file_exists(even) .or. .not. file_exists(odd) ) THROW_HARD('A populated FLEX state has no half maps')
            evens = [evens, even]
            odds  = [odds,  odd]
            w     = [w, real(count(labels == state), dp)]
        enddo
        nmaps = size(w)
        if( nmaps >= 2 ) call write_variance_mask(evens, odds, w, mskdiam, fname%to_char())
    end subroutine write_tree_mask

    function tree_node_rows( node ) result( rows )
        type(flex_tree_node), intent(in) :: node
        logical, allocatable :: rows(:)
        type(sp_project), allocatable :: project
        integer, allocatable :: states(:)
        allocate(project)
        call load_node(node, project, rows, states)
        call project%kill
    end function tree_node_rows

    subroutine write_tree_groups( fit_dir, fit_project, destination, ngroups )
        use simple_defs,          only: PI
        use simple_fileio,        only: fopen, fclose, fileiochk
        use simple_variance_mask, only: variance_mask_groups
        type(string), intent(in)  :: fit_dir, fit_project, destination
        integer,      intent(out) :: ngroups
        type(sp_project), allocatable :: project
        real,    allocatable :: z(:,:), dirs(:,:)
        integer, allocatable :: pinds(:), labels(:), child(:), groups(:)
        character(len=8192)  :: line
        real    :: e(3), u(2)
        integer :: funit, ios, n, q, k, j, g, i, ieo, ilab
        call fopen(funit, string(fit_dir%to_char()//'/flex_pca_coordinates.txt'), status='old', action='read', iostat=ios)
        call fileiochk('write_tree_groups: cannot read the fit''s flex_pca_coordinates.txt', ios)
        n = 0
        q = 0
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( line(1:1) == '#' .or. len_trim(line) == 0 ) cycle
            n = n + 1
            if( q == 0 ) q = nwords(line) - 5
        end do
        rewind(funit)
        allocate(pinds(n), z(n,q))
        k = 0
        do while( k < n )
            read(funit,'(A)') line
            if( line(1:1) == '#' .or. len_trim(line) == 0 ) cycle
            k = k + 1
            read(line,*) pinds(k), ieo, ilab, u(1), u(2), z(k,:)
        end do
        call fclose(funit)
        allocate(project)
        call project%read(fit_project)
        allocate(dirs(3,n))
        do j = 1, n
            e         = project%os_ptcl3D%get_euler(pinds(j)) * PI / 180.
            dirs(:,j) = [sin(e(2))*cos(e(1)), sin(e(2))*sin(e(1)), cos(e(2))]
        end do
        labels = variance_mask_groups(z, dirs)
        allocate(child(project%os_ptcl3D%get_noris()), source=0)
        child(pinds) = labels
        do i = 1, size(child)
            call project%os_ptcl3D%set_state(i, child(i))
        end do
        groups  = pack([(g, g=1,maxval(child))], [(count(child == g) > 0, g=1,maxval(child))])
        ngroups = size(groups)
        call write_stage_project(project, destination, child > 0, groups, .false.)
        call project%kill
        deallocate(project)

    contains

        pure integer function nwords( str )
            character(len=*), intent(in) :: str
            integer :: c
            nwords = 0
            do c = 1, len_trim(str)
                if( str(c:c) /= ' ' .and. (c == 1 .or. str(max(1,c-1):max(1,c-1)) == ' ') ) nwords = nwords + 1
            end do
        end function nwords

    end subroutine write_tree_groups

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

    subroutine write_tree_record( fname, head, tags, nodes )
        use simple_fileio, only: fopen, fclose, fileiochk
        use simple_syslib, only: simple_rename
        use simple_string_utils, only: int2str
        type(string),         intent(in) :: fname
        type(string),         intent(in) :: head(:), tags(:)
        type(flex_tree_node), intent(in) :: nodes(:)
        type(string) :: tmp
        integer :: funit, ios, i, j
        character(len=:), allocatable :: t, st
        tmp = fname%to_char()//'.tmp'
        call fopen(funit, tmp, status='replace', action='write', iostat=ios)
        call fileiochk('write_tree_record: cannot write '//tmp%to_char(), ios)
        do i = 1, size(head)
            write(funit,'(A)') head(i)%to_char()
        end do
        do i = 1, size(nodes)
            t = tags(i)%to_char()
            call put(t, 'name', nodes(i)%name)
            call put(t, 'project_file', nodes(i)%project_file)
            if( allocated(nodes(i)%project_states) )then
                st = ''
                do j = 1, size(nodes(i)%project_states)
                    st = st//int2str(nodes(i)%project_states(j))//','
                end do
                write(funit,'(A)') t//' project_states '//st
            endif
            call put(t, 'map_file', nodes(i)%map_file)
            call put(t, 'even_file', nodes(i)%even_file)
            call put(t, 'odd_file', nodes(i)%odd_file)
            call put(t, 'mask_file', nodes(i)%mask_file)
            call put(t, 'fsc_file', nodes(i)%fsc_file)
            if( allocated(nodes(i)%models) )then
                call put(t, 'model_name', nodes(i)%models(1)%name)
                call put(t, 'model_project', nodes(i)%models(1)%project_file)
                call put(t, 'model_embedding', nodes(i)%models(1)%embedding_file)
                call put(t, 'model_mean', nodes(i)%models(1)%mean_file)
                call put(t, 'model_weights', nodes(i)%models(1)%weights_file)
                call put(t, 'model_targets', nodes(i)%models(1)%targets_file)
                call put(t, 'model_provenance', nodes(i)%models(1)%provenance_file)
            endif
            if( allocated(nodes(i)%decision) )then
                write(funit,'(A)') t//' outcome '//int2str(nodes(i)%decision%outcome)
                call put(t, 'reason', nodes(i)%decision%reason)
            endif
        end do
        call fclose(funit)
        call simple_rename(tmp, fname, overwrite=.true.)

    contains

        subroutine put( tag, key, val )
            character(len=*), intent(in) :: tag, key
            type(string),     intent(in) :: val
            if( .not. val%is_blank() ) write(funit,'(A)') tag//' '//key//' '//val%to_char()
        end subroutine put

    end subroutine write_tree_record

    subroutine read_tree_record( fname, lines )
        use simple_fileio, only: fopen, fclose, fileiochk, nlines
        type(string),              intent(in)  :: fname
        type(string), allocatable, intent(out) :: lines(:)
        character(len=8192) :: buf
        integer :: funit, ios, i
        allocate(lines(nlines(fname)))
        call fopen(funit, fname, status='old', action='read', iostat=ios)
        call fileiochk('read_tree_record: cannot read '//fname%to_char(), ios)
        do i = 1, size(lines)
            read(funit,'(A)') buf
            lines(i) = trim(buf)
        end do
        call fclose(funit)
    end subroutine read_tree_record

    function tree_value( lines, tag, key ) result( val )
        type(string),     intent(in) :: lines(:)
        character(len=*), intent(in) :: tag, key
        type(string) :: val
        character(len=:), allocatable :: line, prefix
        integer :: i
        val    = ''
        prefix = tag//' '//key//' '
        do i = 1, size(lines)
            line = lines(i)%to_char()
            if( len(line) <= len(prefix) ) cycle
            if( line(1:len(prefix)) == prefix ) val = line(len(prefix)+1:)
        end do
    end function tree_value

    subroutine get_tree_node( lines, tag, node )
        use simple_string_utils, only: list_of_ints2arr, str2int
        type(string),         intent(in)  :: lines(:)
        character(len=*),     intent(in)  :: tag
        type(flex_tree_node), intent(out) :: node
        type(string) :: v
        call setv(node%name, 'name')
        call setv(node%project_file, 'project_file')
        v = tree_value(lines, tag, 'project_states')
        if( .not. v%is_blank() ) node%project_states = list_of_ints2arr(v%to_char())
        call setv(node%map_file, 'map_file')
        call setv(node%even_file, 'even_file')
        call setv(node%odd_file, 'odd_file')
        call setv(node%mask_file, 'mask_file')
        call setv(node%fsc_file, 'fsc_file')
        v = tree_value(lines, tag, 'model_project')
        if( .not. v%is_blank() )then
            allocate(node%models(1))
            node%models(1)%project_file = v
            call setv(node%models(1)%name, 'model_name')
            call setv(node%models(1)%embedding_file, 'model_embedding')
            call setv(node%models(1)%mean_file, 'model_mean')
            call setv(node%models(1)%weights_file, 'model_weights')
            call setv(node%models(1)%targets_file, 'model_targets')
            call setv(node%models(1)%provenance_file, 'model_provenance')
        endif
        v = tree_value(lines, tag, 'outcome')
        if( .not. v%is_blank() )then
            allocate(node%decision)
            node%decision%outcome = str2int(v%to_char())
            call setv(node%decision%reason, 'reason')
        endif

    contains

        subroutine setv( dst, key )
            type(string),     intent(inout) :: dst
            character(len=*), intent(in)    :: key
            type(string) :: val
            val = tree_value(lines, tag, key)
            if( .not. val%is_blank() ) dst = val
        end subroutine setv

    end subroutine get_tree_node

end module simple_flex_tree_project_gateway
