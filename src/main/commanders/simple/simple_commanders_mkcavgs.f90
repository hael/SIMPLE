!@descr: for producing class averages
module simple_commanders_mkcavgs
use simple_commanders_api
use simple_pftc_srch_api
use simple_classaverager
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_make_cavgs_distr
  contains
    procedure :: execute      => exec_make_cavgs_distr
end type commander_make_cavgs_distr

type, extends(commander_base) :: commander_make_cavgs
  contains
    procedure :: execute      => exec_make_cavgs
end type commander_make_cavgs

type, extends(commander_base) :: commander_bootstrap_cavgs
  contains
    procedure :: execute      => exec_bootstrap_cavgs
end type commander_bootstrap_cavgs

type, extends(commander_base) :: commander_unbootstrap_cavgs
    contains
        procedure :: execute      => exec_unbootstrap_cavgs
end type commander_unbootstrap_cavgs

type, extends(commander_base) :: commander_cavgassemble
  contains
    procedure :: execute      => exec_cavgassemble
end type commander_cavgassemble

type, extends(commander_base) :: commander_write_classes
  contains
    procedure :: execute      => exec_write_classes
end type commander_write_classes

contains

    subroutine exec_make_cavgs_distr( self, cline )
        class(commander_make_cavgs_distr), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        call run_make_cavgs_workflow(cline, from_distr_cmd=.true.)
    end subroutine exec_make_cavgs_distr

    subroutine exec_make_cavgs( self, cline )
        class(commander_make_cavgs), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        call run_make_cavgs_workflow(cline, from_distr_cmd=.false.)
    end subroutine exec_make_cavgs

    subroutine exec_bootstrap_cavgs( self, cline )
        class(commander_bootstrap_cavgs), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type bootstrap_parent_sample
            integer :: parent_cls = 0
            integer :: out_cls    = 0
            integer :: pop        = 0
            integer :: nanchor    = 0
            integer, allocatable :: ranked(:)
            integer, allocatable :: anchor(:)
            integer, allocatable :: rest(:)
            integer, allocatable :: parts(:,:)
        end type bootstrap_parent_sample
        type bootstrap_row
            integer :: source_pind    = 0
            integer :: synthetic_pind = 0
            integer :: parent_cls     = 0
            integer :: out_cls        = 0
            integer :: child_id       = 0
            integer :: source_rank    = 0
            integer :: eo             = 0
            integer :: stkind_src     = 0
            integer :: indstk         = 0
            real    :: corr           = 0.
            character(len=16) :: role = ''
        end type bootstrap_row
        type bootstrap_skipped_class
            integer :: parent_cls = 0
            integer :: pop        = 0
            integer :: state      = -1
            character(len=64) :: reason = ''
        end type bootstrap_skipped_class
        type(parameters) :: params
        type(sp_project) :: src_proj, boot_proj
        type(cmdline)    :: cline_make
        type(string)     :: manifest_file, membership_file
        type(class_sample), allocatable :: clssmp(:)
        type(bootstrap_parent_sample), allocatable :: parents(:), parents_work(:)
        type(bootstrap_parent_sample)              :: parent_candidate
        type(bootstrap_row),           allocatable :: rows(:)
        type(bootstrap_skipped_class), allocatable :: skipped(:)
        real,    allocatable :: pop_reals(:)
        integer, allocatable :: clsinds(:), pops(:)
        integer :: ncls_acc, nout, nrows, y_child, anchor_target, n_med
        integer :: ncls_in, ncls_cand, nskipped
        integer :: i, j, k, ipart, row_count, child_cls
        logical :: class_ok
        character(len=64) :: skip_reason
        if( .not. cline%defined('osmpl_fac') ) call cline%set('osmpl_fac', 2)
        if( .not. cline%defined('frac_best') ) call cline%set('frac_best', 0.5)
        if( .not. cline%defined('refs')      ) call cline%set('refs', 'cavgs_bootstrap_001.mrc')
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir', 'yes')
        call cline%set('prg',     'bootstrap_cavgs')
        call cline%set('oritype', 'ptcl2D')
        call params%new(cline)
        call cline%set('mkdir', 'no')
        if( params%osmpl_fac < 1 ) THROW_HARD('osmpl_fac must be >= 1; exec_bootstrap_cavgs')
        if( params%frac_best <= 0. .or. params%frac_best > 1. )then
            THROW_HARD('frac_best must be > 0 and <= 1; exec_bootstrap_cavgs')
        endif
        call src_proj%read(params%projfile)
        if( src_proj%os_ptcl2D%get_noris() == 0 ) THROW_HARD('ptcl2D segment is empty; exec_bootstrap_cavgs')
        if( src_proj%os_cls2D%get_noris()  == 0 ) THROW_HARD('cls2D segment is empty; exec_bootstrap_cavgs')
        if( src_proj%os_stk%get_noris()    == 0 ) THROW_HARD('stk segment is empty; exec_bootstrap_cavgs')
        call collect_bootstrap_classes(src_proj, clsinds, pops, skipped, nskipped, ncls_in)
        ncls_cand = size(clsinds)
        y_child  = params%osmpl_fac - 1
        allocate(pop_reals(ncls_cand), source=real(pops))
        n_med         = max(1, nint(median(pop_reals)))
        anchor_target = max(2, nint(params%frac_best * real(n_med)))
        deallocate(pop_reals)
        call src_proj%os_ptcl2D%get_class_sample_stats(clsinds, clssmp)
        allocate(parents_work(ncls_cand))
        ncls_acc = 0
        do i = 1,ncls_cand
            call clear_parent_sample(parent_candidate)
            parent_candidate%parent_cls = clsinds(i)
            parent_candidate%pop        = clssmp(i)%pop
            if( .not. allocated(clssmp(i)%pinds) .or. parent_candidate%pop == 0 )then
                call append_skipped_class(skipped, nskipped, clsinds(i), 0, &
                    class_state_for_report(src_proj, clsinds(i)), 'missing_particle_membership')
                cycle
            endif
            allocate(parent_candidate%ranked(parent_candidate%pop), source=clssmp(i)%pinds)
            call validate_evenodd_present(src_proj%os_ptcl2D, parent_candidate%ranked, &
                parent_candidate%parent_cls, class_ok, skip_reason)
            if( .not. class_ok )then
                call append_skipped_class(skipped, nskipped, parent_candidate%parent_cls, parent_candidate%pop, &
                    class_state_for_report(src_proj, parent_candidate%parent_cls), skip_reason)
                cycle
            endif
            if( y_child > 0 )then
                parent_candidate%nanchor = min(anchor_target, max(0, parent_candidate%pop - y_child))
                if( parent_candidate%nanchor < 2 )then
                    call append_skipped_class(skipped, nskipped, parent_candidate%parent_cls, parent_candidate%pop, &
                        class_state_for_report(src_proj, parent_candidate%parent_cls), &
                        'too_few_particles_for_osmpl_fac')
                    cycle
                endif
                call select_bootstrap_anchor(src_proj%os_ptcl2D, parent_candidate%ranked, parent_candidate%nanchor, &
                    parent_candidate%anchor, parent_candidate%rest, class_ok, skip_reason)
                if( .not. class_ok )then
                    call append_skipped_class(skipped, nskipped, parent_candidate%parent_cls, parent_candidate%pop, &
                        class_state_for_report(src_proj, parent_candidate%parent_cls), skip_reason)
                    cycle
                endif
                if( size(parent_candidate%rest) < y_child )then
                    call append_skipped_class(skipped, nskipped, parent_candidate%parent_cls, parent_candidate%pop, &
                        class_state_for_report(src_proj, parent_candidate%parent_cls), &
                        'too_few_non_anchor_particles')
                    cycle
                endif
                call shuffle_ints(parent_candidate%rest)
                parent_candidate%parts = split_nobjs_even(size(parent_candidate%rest), y_child)
            endif
            ncls_acc = ncls_acc + 1
            parent_candidate%out_cls = ncls_acc
            parents_work(ncls_acc) = parent_candidate
        enddo
        call clear_parent_sample(parent_candidate)
        if( ncls_acc == 0 ) THROW_HARD('no usable cls2D classes found for bootstrap_cavgs')
        allocate(parents(ncls_acc))
        parents = parents_work(:ncls_acc)
        deallocate(parents_work)
        nout  = ncls_acc * params%osmpl_fac
        nrows = sum_parent_pops(parents)
        if( y_child > 0 )then
            do i = 1,ncls_acc
                nrows = nrows + y_child * parents(i)%nanchor + size(parents(i)%rest)
            enddo
        endif
        allocate(rows(nrows))
        row_count = 0
        do i = 1,ncls_acc
            do j = 1,parents(i)%pop
                call append_bootstrap_row(rows, row_count, src_proj, parents(i)%ranked(j), parents(i)%parent_cls, &
                    parents(i)%out_cls, 0, 'original', j)
            enddo
            if( y_child == 0 ) cycle
            do ipart = 1,y_child
                child_cls = ncls_acc + (i - 1) * y_child + ipart
                do j = 1,parents(i)%nanchor
                    call append_bootstrap_row(rows, row_count, src_proj, parents(i)%anchor(j), parents(i)%parent_cls, &
                        child_cls, ipart, 'anchor', rank_of(parents(i)%ranked, parents(i)%anchor(j)))
                enddo
                do k = parents(i)%parts(ipart,1),parents(i)%parts(ipart,2)
                    call append_bootstrap_row(rows, row_count, src_proj, parents(i)%rest(k), parents(i)%parent_cls, &
                        child_cls, ipart, 'member', rank_of(parents(i)%ranked, parents(i)%rest(k)))
                enddo
            enddo
        enddo
        if( row_count /= nrows ) THROW_HARD('internal row count mismatch; exec_bootstrap_cavgs')
        call build_bootstrap_project(src_proj, parents, rows, nout, y_child, boot_proj)
        call derive_bootstrap_output_names(params%refs, manifest_file, membership_file)
        call boot_proj%update_projinfo(params%projfile)
        call boot_proj%write(params%projfile)
        call write_bootstrap_outputs(manifest_file, membership_file, params, parents, rows, skipped, nskipped, &
            ncls_in, nout, y_child, anchor_target)
        call report_bootstrap_class_accounting(manifest_file, ncls_in, ncls_acc, nskipped, nout)
        cline_make = cline
        call cline_make%set('prg',      'make_cavgs')
        call cline_make%set('projfile', params%projfile%to_char())
        call cline_make%set('refs',     params%refs%to_char())
        call cline_make%set('ncls',     nout)
        call cline_make%set('mkdir',    'no')
        call cline_make%set('oritype',  'ptcl2D')
        call run_make_cavgs_workflow(cline_make, from_distr_cmd=.true.)
        call cline_make%kill
        call deallocate_class_samples(clssmp)
        call src_proj%kill
        call boot_proj%kill

    contains

        subroutine collect_bootstrap_classes( spproj, clsinds, pops, skipped, nskipped, ncls_in )
            type(sp_project),              intent(inout) :: spproj
            integer, allocatable,          intent(inout) :: clsinds(:), pops(:)
            type(bootstrap_skipped_class), allocatable, intent(inout) :: skipped(:)
            integer,                       intent(inout) :: nskipped, ncls_in
            integer, allocatable :: cls_tmp(:), pop_tmp(:)
            integer :: icls, nacc, pop, cls_state
            logical :: have_state
            if( allocated(clsinds) ) deallocate(clsinds)
            if( allocated(pops)    ) deallocate(pops)
            if( allocated(skipped) ) deallocate(skipped)
            ncls_in = spproj%os_cls2D%get_noris()
            allocate(cls_tmp(ncls_in), source=0)
            allocate(pop_tmp(ncls_in), source=0)
            allocate(skipped(ncls_in), source=bootstrap_skipped_class())
            nacc = 0
            nskipped = 0
            do icls = 1,ncls_in
                pop = spproj%os_ptcl2D%get_pop(icls, 'class')
                have_state = spproj%os_cls2D%isthere(icls, 'state')
                cls_state  = -1
                if( have_state ) cls_state = spproj%os_cls2D%get_state(icls)
                if( have_state .and. cls_state == 0 )then
                    call append_skipped_class(skipped, nskipped, icls, pop, cls_state, 'inactive_cls2D_state')
                    cycle
                endif
                if( pop == 0 )then
                    call append_skipped_class(skipped, nskipped, icls, pop, cls_state, 'zero_population')
                    cycle
                endif
                nacc = nacc + 1
                cls_tmp(nacc) = icls
                pop_tmp(nacc) = pop
            enddo
            if( nacc == 0 ) THROW_HARD('no active populated cls2D classes found; collect_bootstrap_classes')
            allocate(clsinds(nacc), pops(nacc))
            clsinds = cls_tmp(:nacc)
            pops    = pop_tmp(:nacc)
            deallocate(cls_tmp, pop_tmp)
        end subroutine collect_bootstrap_classes

        subroutine validate_evenodd_present( os, pinds, parent_cls, ok, reason )
            class(oris), intent(inout) :: os
            integer,     intent(in)    :: pinds(:), parent_cls
            logical,          intent(out) :: ok
            character(len=*), intent(out) :: reason
            integer :: i, eo
            ok     = .true.
            reason = ''
            do i = 1,size(pinds)
                eo = os%get_eo(pinds(i))
                if( eo /= 0 .and. eo /= 1 )then
                    ok     = .false.
                    reason = 'invalid_eo_assignment'
                    write(logfhandle,'(A,I0,A,I0,A,I0)') &
                        '>>> BOOTSTRAP_CAVGS SKIPPING CLASS ', parent_cls, &
                        ': particle ', pinds(i), ' has invalid eo=', eo
                    return
                endif
            enddo
        end subroutine validate_evenodd_present

        subroutine select_bootstrap_anchor( os, ranked, nanchor, anchor, rest, ok, reason )
            class(oris),              intent(inout) :: os
            integer,                  intent(in)    :: ranked(:), nanchor
            integer, allocatable,     intent(inout) :: anchor(:), rest(:)
            logical,                  intent(out)   :: ok
            character(len=*),         intent(out)   :: reason
            integer :: n, i, cnt, missing_eo, swap_pind
            logical :: have_even, have_odd
            ok     = .false.
            reason = ''
            if( allocated(anchor) ) deallocate(anchor)
            if( allocated(rest)   ) deallocate(rest)
            n = size(ranked)
            if( nanchor >= n )then
                reason = 'anchor_leaves_no_child_particles'
                return
            endif
            allocate(anchor(nanchor), source=ranked(:nanchor))
            have_even = count_eo(os, anchor, 0) > 0
            have_odd  = count_eo(os, anchor, 1) > 0
            if( .not.(have_even .and. have_odd) )then
                if( have_even )then
                    missing_eo = 1
                else
                    missing_eo = 0
                endif
                swap_pind = 0
                do i = nanchor + 1,n
                    if( os%get_eo(ranked(i)) == missing_eo )then
                        swap_pind = ranked(i)
                        exit
                    endif
                enddo
                if( swap_pind == 0 )then
                    reason = 'cannot_make_parity_complete_anchor'
                    return
                endif
                anchor(nanchor) = swap_pind
            endif
            if( count_eo(os, anchor, 0) == 0 .or. count_eo(os, anchor, 1) == 0 )then
                reason = 'cannot_make_parity_complete_anchor'
                return
            endif
            allocate(rest(n - nanchor), source=0)
            cnt = 0
            do i = 1,n
                if( any(anchor == ranked(i)) ) cycle
                cnt = cnt + 1
                rest(cnt) = ranked(i)
            enddo
            if( cnt /= n - nanchor ) THROW_HARD('anchor/rest partition mismatch; select_bootstrap_anchor')
            ok = .true.
        end subroutine select_bootstrap_anchor

        integer function count_eo( os, pinds, eo )
            class(oris), intent(inout) :: os
            integer,     intent(in)    :: pinds(:), eo
            integer :: i
            count_eo = 0
            do i = 1,size(pinds)
                if( os%get_eo(pinds(i)) == eo ) count_eo = count_eo + 1
            enddo
        end function count_eo

        subroutine shuffle_ints( vals )
            integer, allocatable, intent(inout) :: vals(:)
            type(ran_tabu) :: rt
            if( size(vals) < 2 ) return
            rt = ran_tabu(size(vals))
            call rt%shuffle(vals)
            call rt%kill
        end subroutine shuffle_ints

        integer function rank_of( ranked, pind )
            integer, intent(in) :: ranked(:), pind
            integer :: i
            rank_of = 0
            do i = 1,size(ranked)
                if( ranked(i) == pind )then
                    rank_of = i
                    return
                endif
            enddo
            THROW_HARD('particle rank lookup failed; rank_of')
        end function rank_of

        subroutine append_bootstrap_row( rows, row_count, src_proj, source_pind, parent_cls, out_cls, child_id, role, source_rank )
            type(bootstrap_row), intent(inout) :: rows(:)
            integer,             intent(inout) :: row_count
            type(sp_project),    intent(inout) :: src_proj
            integer,             intent(in)    :: source_pind, parent_cls, out_cls, child_id, source_rank
            character(len=*),    intent(in)    :: role
            row_count = row_count + 1
            if( row_count > size(rows) ) THROW_HARD('row_count exceeds rows size; append_bootstrap_row')
            rows(row_count)%source_pind = source_pind
            rows(row_count)%parent_cls  = parent_cls
            rows(row_count)%out_cls     = out_cls
            rows(row_count)%child_id    = child_id
            rows(row_count)%role        = role
            rows(row_count)%source_rank = source_rank
            rows(row_count)%eo          = src_proj%os_ptcl2D%get_eo(source_pind)
            rows(row_count)%corr        = src_proj%os_ptcl2D%get(source_pind, 'corr')
            call src_proj%map_ptcl_ind2stk_ind('ptcl2D', source_pind, rows(row_count)%stkind_src, rows(row_count)%indstk)
        end subroutine append_bootstrap_row

        subroutine build_bootstrap_project( src_proj, parents, rows, nout, y_child, boot_proj )
            type(sp_project),                  intent(inout) :: src_proj, boot_proj
            type(bootstrap_parent_sample),     intent(in)    :: parents(:)
            type(bootstrap_row),               intent(inout) :: rows(:)
            integer,                           intent(in)    :: nout, y_child
            integer, allocatable :: stk_counts(:), stk_newinds(:), stk_next(:), row_order(:)
            integer :: nstks_src, nstks_used, istk, i, row_out, stkind_new, fromp, top, row_pos, src_row
            integer :: icls, ipart, child_cls, child_sz
            nstks_src = src_proj%os_stk%get_noris()
            allocate(stk_counts(nstks_src), source=0)
            allocate(stk_newinds(nstks_src), source=0)
            do i = 1,size(rows)
                if( rows(i)%stkind_src < 1 .or. rows(i)%stkind_src > nstks_src )then
                    THROW_HARD('source stack index out of range; build_bootstrap_project')
                endif
                stk_counts(rows(i)%stkind_src) = stk_counts(rows(i)%stkind_src) + 1
            enddo
            nstks_used = count(stk_counts > 0)
            allocate(stk_next(nstks_src),    source=0)
            allocate(row_order(size(rows)),  source=0)
            call boot_proj%copy(src_proj)
            call boot_proj%os_mic%kill
            call boot_proj%os_out%kill
            call boot_proj%os_stk%new(nstks_used, is_ptcl=.false.)
            call boot_proj%os_ptcl2D%new(size(rows), is_ptcl=.true.)
            call boot_proj%os_cls2D%new(nout, is_ptcl=.false.)
            stkind_new = 0
            top        = 0
            do istk = 1,nstks_src
                if( stk_counts(istk) == 0 ) cycle
                stkind_new = stkind_new + 1
                stk_newinds(istk) = stkind_new
                fromp = top + 1
                top   = top + stk_counts(istk)
                call boot_proj%os_stk%transfer_ori(stkind_new, src_proj%os_stk, istk)
                call boot_proj%os_stk%set(stkind_new, 'fromp',  fromp)
                call boot_proj%os_stk%set(stkind_new, 'top',    top)
                call boot_proj%os_stk%set(stkind_new, 'nptcls', stk_counts(istk))
                stk_next(istk) = fromp
            enddo
            do i = 1,size(rows)
                istk = rows(i)%stkind_src
                row_order(stk_next(istk)) = i
                stk_next(istk) = stk_next(istk) + 1
            enddo
            row_out = 0
            do row_pos = 1,size(row_order)
                src_row = row_order(row_pos)
                if( src_row == 0 ) THROW_HARD('stack row ordering failed; build_bootstrap_project')
                istk    = rows(src_row)%stkind_src
                row_out = row_out + 1
                rows(src_row)%synthetic_pind = row_out
                call boot_proj%os_ptcl2D%transfer_ori(row_out, src_proj%os_ptcl2D, rows(src_row)%source_pind)
                call boot_proj%os_ptcl2D%set(row_out, 'stkind',    stk_newinds(istk))
                call boot_proj%os_ptcl2D%set(row_out, 'indstk',    rows(src_row)%indstk)
                call boot_proj%os_ptcl2D%set(row_out, 'class',     rows(src_row)%out_cls)
                call boot_proj%os_ptcl2D%set(row_out, 'state',     1)
                call boot_proj%os_ptcl2D%set(row_out, 'pind',      row_out)
                call boot_proj%os_ptcl2D%set(row_out, 'pind_prev', rows(src_row)%source_pind)
            enddo
            if( row_out /= size(rows) ) THROW_HARD('synthetic row count mismatch; build_bootstrap_project')
            boot_proj%os_ptcl3D = boot_proj%os_ptcl2D
            do icls = 1,size(parents)
                call set_bootstrap_class_row(boot_proj, src_proj, parents(icls)%out_cls, parents(icls)%parent_cls, &
                    0, parents(icls)%pop, parents(icls)%nanchor)
                if( y_child == 0 ) cycle
                do ipart = 1,y_child
                    child_cls = size(parents) + (icls - 1) * y_child + ipart
                    child_sz  = parents(icls)%nanchor + parents(icls)%parts(ipart,2) - parents(icls)%parts(ipart,1) + 1
                    call set_bootstrap_class_row(boot_proj, src_proj, child_cls, parents(icls)%parent_cls, &
                        ipart, child_sz, parents(icls)%nanchor)
                enddo
            enddo
            boot_proj%os_cls3D = boot_proj%os_cls2D
            deallocate(stk_counts, stk_newinds, stk_next, row_order)
        end subroutine build_bootstrap_project

        subroutine set_bootstrap_class_row( boot_proj, src_proj, out_cls, parent_cls, child_id, pop, nanchor )
            type(sp_project), intent(inout) :: boot_proj, src_proj
            integer,          intent(in)    :: out_cls, parent_cls, child_id, pop, nanchor
            call boot_proj%os_cls2D%transfer_ori(out_cls, src_proj%os_cls2D, parent_cls)
            call boot_proj%os_cls2D%set(out_cls, 'class', out_cls)
            call boot_proj%os_cls2D%set(out_cls, 'state', 1)
            call boot_proj%os_cls2D%set(out_cls, 'pop', pop)
            call boot_proj%os_cls2D%set(out_cls, 'bootstrap_parent', parent_cls)
            call boot_proj%os_cls2D%set(out_cls, 'bootstrap_child',  child_id)
            call boot_proj%os_cls2D%set(out_cls, 'bootstrap_anchor', nanchor)
        end subroutine set_bootstrap_class_row

        subroutine derive_bootstrap_output_names( refs, manifest_file, membership_file )
            class(string), intent(in)    :: refs
            type(string),  intent(inout) :: manifest_file, membership_file
            type(string) :: refs_ext, refs_body
            refs_ext = fname2ext(refs)
            if( refs_ext%strlen_trim() > 0 )then
                refs_body = get_fbody(refs, refs_ext)
            else
                refs_body = refs
            endif
            manifest_file  = refs_body%to_char()//'_manifest.txt'
            membership_file = refs_body%to_char()//'_membership.txt'
            call refs_ext%kill
            call refs_body%kill
        end subroutine derive_bootstrap_output_names

        subroutine write_bootstrap_outputs( manifest_file, membership_file, params, parents, rows, skipped, nskipped, &
                ncls_in, nout, y_child, anchor_target )
            type(string),                    intent(in) :: manifest_file, membership_file
            type(parameters),                intent(in) :: params
            type(bootstrap_parent_sample),   intent(in) :: parents(:)
            type(bootstrap_row),             intent(in) :: rows(:)
            type(bootstrap_skipped_class),   intent(in) :: skipped(:)
            integer,                         intent(in) :: nskipped, ncls_in, nout, y_child, anchor_target
            integer :: funit, i, ipart, child_cls, child_sz
            call fopen(funit, FILE=manifest_file, STATUS='REPLACE', action='WRITE')
            write(funit,'(a)') '# SIMPLE bootstrap_cavgs manifest'
            write(funit,'(a,a)') '# refs: ', params%refs%to_char()
            write(funit,'(a,i0)') '# ncls_input: ', ncls_in
            write(funit,'(a,i0)') '# ncls_used: ', size(parents)
            write(funit,'(a,i0)') '# ncls_skipped: ', nskipped
            write(funit,'(a,i0)') '# ncls_expected_if_all_input_valid: ', ncls_in * params%osmpl_fac
            write(funit,'(a,i0)') '# ncls_out: ', nout
            write(funit,'(a,i0)') '# osmpl_fac: ', params%osmpl_fac
            write(funit,'(a,f8.4)') '# frac_best: ', params%frac_best
            write(funit,'(a,i0)') '# anchor_target: ', anchor_target
            if( nskipped > 0 )then
                write(funit,'(a)') '# skipped_classes: parent_cls parent_pop cls_state reason'
                do i = 1,nskipped
                    write(funit,'(a,1x,i0,1x,i0,1x,i0,1x,a)') '# skipped', skipped(i)%parent_cls, &
                        skipped(i)%pop, skipped(i)%state, trim(skipped(i)%reason)
                enddo
            endif
            write(funit,'(a)') 'out_cls role parent_cls child_id parent_pop n_anchor n_members'
            do i = 1,size(parents)
                write(funit,'(i0,1x,a,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0)') parents(i)%out_cls, 'original', &
                    parents(i)%parent_cls, 0, parents(i)%pop, parents(i)%nanchor, parents(i)%pop
                if( y_child == 0 ) cycle
                do ipart = 1,y_child
                    child_cls = size(parents) + (i - 1) * y_child + ipart
                    child_sz  = parents(i)%nanchor + parents(i)%parts(ipart,2) - parents(i)%parts(ipart,1) + 1
                    write(funit,'(i0,1x,a,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0)') child_cls, 'child', &
                        parents(i)%parent_cls, ipart, parents(i)%pop, parents(i)%nanchor, child_sz
                enddo
            enddo
            call fclose(funit)
            call fopen(funit, FILE=membership_file, STATUS='REPLACE', action='WRITE')
            write(funit,'(a)') 'synthetic_pind out_cls role parent_cls child_id source_pind source_rank source_corr eo stkind indstk'
            do i = 1,size(rows)
                write(funit,'(i0,1x,i0,1x,a,1x,i0,1x,i0,1x,i0,1x,i0,1x,f12.6,1x,i0,1x,i0,1x,i0)') &
                    rows(i)%synthetic_pind, rows(i)%out_cls, rows(i)%role, rows(i)%parent_cls, rows(i)%child_id, &
                    rows(i)%source_pind, rows(i)%source_rank, rows(i)%corr, rows(i)%eo, rows(i)%stkind_src, rows(i)%indstk
            enddo
            call fclose(funit)
        end subroutine write_bootstrap_outputs

        subroutine append_skipped_class( skipped, nskipped, parent_cls, pop, state, reason )
            type(bootstrap_skipped_class), intent(inout) :: skipped(:)
            integer,                       intent(inout) :: nskipped
            integer,                       intent(in)    :: parent_cls, pop, state
            character(len=*),              intent(in)    :: reason
            nskipped = nskipped + 1
            if( nskipped > size(skipped) ) THROW_HARD('skipped-class accounting overflow; append_skipped_class')
            skipped(nskipped)%parent_cls = parent_cls
            skipped(nskipped)%pop        = pop
            skipped(nskipped)%state      = state
            skipped(nskipped)%reason     = reason
        end subroutine append_skipped_class

        integer function class_state_for_report( spproj, icls )
            type(sp_project), intent(inout) :: spproj
            integer,          intent(in)    :: icls
            class_state_for_report = -1
            if( spproj%os_cls2D%isthere(icls, 'state') )then
                class_state_for_report = spproj%os_cls2D%get_state(icls)
            endif
        end function class_state_for_report

        integer function sum_parent_pops( parents )
            type(bootstrap_parent_sample), intent(in) :: parents(:)
            integer :: i
            sum_parent_pops = 0
            do i = 1,size(parents)
                sum_parent_pops = sum_parent_pops + parents(i)%pop
            enddo
        end function sum_parent_pops

        subroutine clear_parent_sample( parent )
            type(bootstrap_parent_sample), intent(inout) :: parent
            parent%parent_cls = 0
            parent%out_cls    = 0
            parent%pop        = 0
            parent%nanchor    = 0
            if( allocated(parent%ranked) ) deallocate(parent%ranked)
            if( allocated(parent%anchor) ) deallocate(parent%anchor)
            if( allocated(parent%rest)   ) deallocate(parent%rest)
            if( allocated(parent%parts)  ) deallocate(parent%parts)
        end subroutine clear_parent_sample

        subroutine report_bootstrap_class_accounting( manifest_file, ncls_in, ncls_used, nskipped, nout )
            type(string), intent(in) :: manifest_file
            integer,      intent(in) :: ncls_in, ncls_used, nskipped, nout
            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') &
                '>>> BOOTSTRAP_CAVGS CLASSES INPUT/USED/SKIPPED/OUTPUT: ', ncls_in, '/', &
                ncls_used, '/', nskipped, '/', nout
            if( nskipped > 0 )then
                write(logfhandle,'(A,A)') '>>> BOOTSTRAP_CAVGS SKIPPED CLASS DETAILS WRITTEN TO: ', &
                    manifest_file%to_char()
            endif
            call flush(logfhandle)
        end subroutine report_bootstrap_class_accounting

    end subroutine exec_bootstrap_cavgs

    subroutine exec_unbootstrap_cavgs( self, cline )
        class(commander_unbootstrap_cavgs), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: boot_proj, src_proj
        type(string)     :: src_projfile
        integer, allocatable :: seen_parent(:)
        real, allocatable    :: states(:)
        integer :: ncls_src, ncls_boot, i, parent_cls, child_id, nmap
        real    :: corr, proj, state
        call params%new(cline)
        if( .not. cline%defined('projfile_orig') )then
            THROW_HARD('missing required key projfile_orig; usage: simple_exec prg=unbootstrap_cavgs projfile=<bootstrap.simple> projfile_orig=<original.simple>')
        endif
        src_projfile = cline%get_carg('projfile_orig')
        call boot_proj%read(params%projfile)
        call src_proj%read(src_projfile)
        if( boot_proj%os_cls2D%get_noris() == 0 ) THROW_HARD('empty cls2D in bootstrap project; exec_unbootstrap_cavgs')
        if( boot_proj%os_cls3D%get_noris() == 0 ) THROW_HARD('empty cls3D in bootstrap project; exec_unbootstrap_cavgs')
        if( .not. boot_proj%os_cls2D%isthere('bootstrap_parent') )then
            THROW_HARD('bootstrap_parent missing in cls2D; not a bootstrap_cavgs project')
        endif
        if( .not. boot_proj%os_cls2D%isthere('bootstrap_child') )then
            THROW_HARD('bootstrap_child missing in cls2D; not a bootstrap_cavgs project')
        endif
        ncls_src  = src_proj%os_cls2D%get_noris()
        ncls_boot = boot_proj%os_cls2D%get_noris()
        if( ncls_src == 0 ) THROW_HARD('empty cls2D in source project; exec_unbootstrap_cavgs')
        if( src_proj%os_cls3D%get_noris() /= ncls_src )then
            call src_proj%os_cls3D%new(ncls_src, is_ptcl=.false.)
            states = src_proj%os_cls2D%get_all('state')
            call src_proj%os_cls3D%set_all('state', states)
            deallocate(states)
        endif
        call src_proj%os_cls3D%delete_3Dalignment()
        allocate(states(ncls_src), source=0.)
        call src_proj%os_cls3D%set_all('state', states)
        deallocate(states)
        allocate(seen_parent(ncls_src), source=0)
        nmap = 0
        do i = 1,ncls_boot
            child_id = nint(boot_proj%os_cls2D%get(i, 'bootstrap_child'))
            if( child_id /= 0 ) cycle
            parent_cls = nint(boot_proj%os_cls2D%get(i, 'bootstrap_parent'))
            if( parent_cls < 1 .or. parent_cls > ncls_src )then
                THROW_HARD('bootstrap_parent out of range in cls2D; exec_unbootstrap_cavgs')
            endif
            if( seen_parent(parent_cls) /= 0 )then
                THROW_HARD('duplicate bootstrap_parent with bootstrap_child=0; exec_unbootstrap_cavgs')
            endif
            seen_parent(parent_cls) = 1
            corr  = boot_proj%os_cls3D%get(i, 'corr')
            proj  = boot_proj%os_cls3D%get(i, 'proj')
            state = boot_proj%os_cls3D%get(i, 'state')
            call src_proj%os_cls3D%set(parent_cls, 'corr', corr)
            call src_proj%os_cls3D%set(parent_cls, 'proj', proj)
            call src_proj%os_cls3D%set_euler(parent_cls, boot_proj%os_cls3D%get_euler(i))
            call src_proj%os_cls3D%set_shift(parent_cls, boot_proj%os_cls3D%get_2Dshift(i))
            call src_proj%os_cls3D%set(parent_cls, 'state', state)
            nmap = nmap + 1
        enddo
        if( nmap == 0 ) THROW_HARD('no bootstrap_child=0 classes found to map; exec_unbootstrap_cavgs')
        call src_proj%map2ptcls
        call src_proj%write(src_projfile)
        call boot_proj%kill
        call src_proj%kill
        call simple_end('**** SIMPLE_UNBOOTSTRAP_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_unbootstrap_cavgs

    ! ------------------------------------------------------------------
    ! Unified runtime-polymorphic workflow
    ! ------------------------------------------------------------------

    subroutine run_make_cavgs_workflow( cline, from_distr_cmd )
        use simple_make_cavgs_strategy, only: make_cavgs_strategy, make_cavgs_hooks, create_make_cavgs_strategy
        use simple_cmdline,             only: cmdline
        use simple_parameters,          only: parameters
        class(cmdline), intent(inout) :: cline
        logical,        intent(in)    :: from_distr_cmd
        class(make_cavgs_strategy), allocatable :: strategy
        type(make_cavgs_hooks) :: hooks
        type(parameters) :: params
        ! Ensure distributed scripts see the correct program name.
        call cline%set('prg',    'make_cavgs')
        call cline%set('oritype','ptcl2D')
        ! Provide master-side hook for cavgassemble (avoids strategy importing commander modules).
        hooks%run_cavgassemble => make_cavgs_exec_cavgassemble
        strategy = create_make_cavgs_strategy(cline, hooks, from_distr_cmd=from_distr_cmd)
        call strategy%apply_defaults(cline)
        call strategy%initialize(params, cline)
        if( params%l_nonuniform ) THROW_HARD('2D nonuniform filtering has been removed; run_make_cavgs_workflow')
        call strategy%execute(params, cline)
        call strategy%finalize_run(params, cline)
        call strategy%cleanup(params, cline)
        call simple_end(strategy%end_message(), print_simple=.false.)
        ! exec_make_cavgs_distr behavior: async touch marker
        call strategy%after_end(params, cline)
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine run_make_cavgs_workflow

    ! ----------------------------------------------------------------------
    ! Hook implementation: run cavgassemble using the existing commander
    ! ----------------------------------------------------------------------

    subroutine make_cavgs_exec_cavgassemble( cline, nthr )
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: nthr
        type(commander_cavgassemble) :: xcavgassemble
        type(cmdline)               :: cline_cavgassemble
        cline_cavgassemble = cline
        call cline_cavgassemble%set('prg',  'cavgassemble')
        call cline_cavgassemble%set('nthr', nthr)
        call xcavgassemble%execute(cline_cavgassemble)
        call cline_cavgassemble%kill
    end subroutine make_cavgs_exec_cavgassemble

    subroutine exec_cavgassemble( self, cline )
        use simple_timer, only: timer_int_kind, tic, toc
        class(commander_cavgassemble), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)   :: params
        type(builder)      :: build
        type(starproject)  :: starproj
        type(string)       :: benchfname
        real, allocatable  :: states(:)
        integer            :: iterstr_start, iterstr_end, iter, io_stat, fnr
        integer(timer_int_kind) :: t_tot, t_phase
        real(timer_int_kind)    :: rt_init, rt_reduce_partials, rt_classdoc, rt_write_cavgs
        real(timer_int_kind)    :: rt_export_star, rt_project_update, rt_cleanup, rt_tot
        if( L_BENCH_GLOB )then
            rt_init           = 0.
            rt_reduce_partials= 0.
            rt_classdoc       = 0.
            rt_write_cavgs    = 0.
            rt_export_star    = 0.
            rt_project_update = 0.
            rt_cleanup        = 0.
            t_tot             = tic()
            t_phase           = tic()
        endif
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
        if( L_BENCH_GLOB ) rt_init = toc(t_phase)
        if( cline%defined('which_iter') )then
            params%refs      = CAVGS_ITER_FBODY//int2str_pad(params%which_iter,3)//params%ext%to_char()
            params%refs_even = CAVGS_ITER_FBODY//int2str_pad(params%which_iter,3)//'_even'//params%ext%to_char()
            params%refs_odd  = CAVGS_ITER_FBODY//int2str_pad(params%which_iter,3)//'_odd'//params%ext%to_char()
        else if( .not. cline%defined('refs') )then
            params%refs      = 'start2Drefs'//params%ext%to_char()
            params%refs_even = 'start2Drefs_even'//params%ext%to_char()
            params%refs_odd  = 'start2Drefs_odd'//params%ext%to_char()
        endif
        if( L_BENCH_GLOB ) t_phase = tic()
        call cavger_new(params, build)
        call cavger_assemble_sums_from_parts
        if( L_BENCH_GLOB ) rt_reduce_partials = toc(t_phase)
        if( L_BENCH_GLOB ) t_phase = tic()
        call cavger_gen2Dclassdoc
        if( L_BENCH_GLOB ) rt_classdoc = toc(t_phase)
        call terminate_stream(params, 'SIMPLE_CAVGASSEMBLE HARD STOP')
        if( L_BENCH_GLOB ) t_phase = tic()
        call cavger_write_all(params%refs, params%refs_even, params%refs_odd)
        call cavger_kill
        if( L_BENCH_GLOB ) rt_write_cavgs = toc(t_phase)
        ! get iteration from which_iter else from refs filename and write cavgs starfile
        if( cline%defined('which_iter') ) then
            if( L_BENCH_GLOB ) t_phase = tic()
            call starproj%export_cls2D(build%spproj, params%which_iter)
            if( L_BENCH_GLOB ) rt_export_star = toc(t_phase)
            iter = params%which_iter
        else if( cline%defined('refs') .and. params%refs%substr_ind(CAVGS_ITER_FBODY) > 0 ) then
            iterstr_start = params%refs%substr_ind(CAVGS_ITER_FBODY) + 10
            iterstr_end   = params%refs%substr_ind(params%ext) - 1
            iter = str2int(params%refs%to_char([iterstr_start,iterstr_end]), io_stat)
            if( L_BENCH_GLOB ) t_phase = tic()
            call starproj%export_cls2D(build%spproj, iter)
            if( L_BENCH_GLOB ) rt_export_star = toc(t_phase)
        else
            iter = max(1, params%which_iter)
        end if
        ! updates project
        ! cls2D and state congruent cls3D
        if( L_BENCH_GLOB ) t_phase = tic()
        call build%spproj%os_cls3D%new(params%ncls, is_ptcl=.false.)
        states = build%spproj%os_cls2D%get_all('state')
        call build%spproj%os_cls3D%set_all('state',states)
        deallocate(states)
        call build%spproj%add_frcs2os_out( string(FRCS_FILE), 'frc2D')
        call build%spproj%add_cavgs2os_out(params%refs, build%spproj%get_smpd(), imgkind='cavg')
        ! multiple fields updated, do a full write
        call build%spproj%write(params%projfile)
        if( L_BENCH_GLOB ) rt_project_update = toc(t_phase)
        ! end gracefully
        if( L_BENCH_GLOB ) t_phase = tic()
        call starproj%kill
        call build%spproj%kill
        call build%kill_general_tbox
        call build%kill_strategy2D_tbox
        if( L_BENCH_GLOB )then
            rt_cleanup = toc(t_phase)
            rt_tot     = toc(t_tot)
            benchfname = string('CAVGASSEMBLE_BENCH_ITER')//int2str_pad(iter,3)//'.txt'
            call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
            write(fnr,'(a,a)')  'cavgassemble assembly mode          : class'
            write(fnr,'(a,a)')  'cavgassemble refs                   : ', params%refs%to_char()
            write(fnr,'(a,i0)') 'cavgassemble nclasses               : ', params%ncls
            write(fnr,'(a,i0)') 'cavgassemble nparts                 : ', params%nparts
            write(fnr,'(a,i0)') 'cavgassemble kfrom                  : ', params%kfromto(1)
            write(fnr,'(a,i0)') 'cavgassemble kto                    : ', params%kfromto(2)
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,t52,f9.2)') 'cavgassemble setup/init             : ', rt_init
            write(fnr,'(a,t52,f9.2)') 'cavgassemble reduce partial sums    : ', rt_reduce_partials
            write(fnr,'(a,t52,f9.2)') 'cavgassemble class document         : ', rt_classdoc
            write(fnr,'(a,t52,f9.2)') 'cavgassemble write class averages   : ', rt_write_cavgs
            write(fnr,'(a,t52,f9.2)') 'cavgassemble export cls2D star      : ', rt_export_star
            write(fnr,'(a,t52,f9.2)') 'cavgassemble project update         : ', rt_project_update
            write(fnr,'(a,t52,f9.2)') 'cavgassemble cleanup                : ', rt_cleanup
            write(fnr,'(a,t52,f9.2)') 'cavgassemble total time             : ', rt_tot
            write(fnr,'(a,t52,f9.2)') 'cavgassemble % accounted for        : ', &
                &((rt_init + rt_reduce_partials + rt_classdoc + rt_write_cavgs + &
                &  rt_export_star + rt_project_update + rt_cleanup) / rt_tot) * 100.
            call fclose(fnr)
            call benchfname%kill
        endif
        call simple_end('**** SIMPLE_CAVGASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        call simple_touch('CAVGASSEMBLE_FINISHED')
    end subroutine exec_cavgassemble

    subroutine exec_write_classes( self, cline )
        class(commander_write_classes), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        type(image)      :: img_cavg
        type(string)     :: cavgsstk, classname, stkname
        type(image),        allocatable :: imgs_class(:)
        real,               allocatable :: states(:), inpls(:,:)
        integer,            allocatable :: pops(:), pinds(:)
        real(kind=c_float), allocatable :: rmat_rot(:,:,:)
        integer :: ncls, n, ldim(3), icls, pop_max, ind_in_stk, i, cnt
        real    :: smpd
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call spproj%read(params%projfile)
        ! get class average stack
        call spproj%get_cavgs_stk(cavgsstk, ncls, smpd)
        call find_ldim_nptcls(cavgsstk, ldim, n)
        ldim(3) = 1
        if( n /= ncls ) THROW_HARD('Incosistent # classes in project file vs cavgs stack; exec_write_classes')
        ! get state flag array
        states = spproj%os_cls2D%get_all('state')
        ! get ncls from ptcl2D field
        n = spproj%os_ptcl2D%get_n('class')
        if( n /= ncls ) THROW_HARD('Incosistent # classes in ptcl2D field of spproj vs cavgs stack; exec_write_classes')
        ! find out maximum population and allocate image arrays accordingly
        call spproj%os_ptcl2D%get_pops(pops, 'class')
        pop_max = maxval(pops)
        write(logfhandle,'(A,I5)') '>>> MAXIMUM CLASS POPULATION: ', pop_max
        allocate(imgs_class(pop_max), inpls(pop_max,3), rmat_rot(ldim(1),ldim(2),1))
        rmat_rot = 0.
        inpls    = 0.
        do i=1,pop_max
            call imgs_class(i)%new(ldim, smpd, wthreads=.false.)
        end do
        call img_cavg%new(ldim, smpd)
        ! loop over classes
        do icls=1,ncls
            if( states(icls) < 0.5 ) cycle
            ! get particle indices of class
            call spproj%os_ptcl2D%get_pinds(icls, 'class', pinds)
            if( .not. allocated(pinds) ) cycle
            ! read the class average
            call img_cavg%read(cavgsstk, icls)
            ! read the images and get the in-plane parameters
            do i=1,size(pinds)
                ! read
                call spproj%get_stkname_and_ind('ptcl2D', pinds(i), stkname, ind_in_stk)
                call imgs_class(i)%read(stkname, ind_in_stk)
                ! get params
                inpls(i,1)  = spproj%os_ptcl2D%e3get(pinds(i))
                inpls(i,2:) = spproj%os_ptcl2D%get_2Dshift(pinds(i))
            end do
            ! rotate the images (in parallel)
            !$omp parallel do default(shared) private(i,rmat_rot) schedule(static) proc_bind(close)
            do i=1,size(pinds)
                call imgs_class(i)%fft
                call imgs_class(i)%shift2Dserial([-inpls(i,2),-inpls(i,3)])
                call imgs_class(i)%ifft
                call imgs_class(i)%rtsq_serial(inpls(i,1), 0., 0., rmat_rot)
                call imgs_class(i)%set_rmat(rmat_rot,.false.)
            end do
            !$omp end parallel do
            ! make a filename for the class
            classname = 'class'//int2str_pad(icls,5)//STK_EXT
            ! write the class average first, followed by the rotated and shifted particles
            call img_cavg%write(classname, 1)
            cnt = 1
            do i=1,size(pinds)
                cnt = cnt + 1
                call imgs_class(i)%write(classname, cnt)
            end do
        end do
        ! destruct
        call spproj%kill
        call img_cavg%kill
        do i=1,size(imgs_class)
            call imgs_class(i)%kill
        end do
        deallocate(imgs_class, inpls, rmat_rot)
        if( allocated(states) ) deallocate(states)
        if( allocated(pops)   ) deallocate(pops)
        if( allocated(pinds)  ) deallocate(pinds)
        ! end gracefully
        call simple_end('**** SIMPLE_WRITE_CLASSES NORMAL STOP ****')
    end subroutine exec_write_classes

end module simple_commanders_mkcavgs
