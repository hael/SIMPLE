!@descr: direction-balanced top-score particle selection and congruent project, stack and parameter export
module simple_ptcl_subset
use simple_sp_project, only: sp_project
use simple_string,     only: string
implicit none
private
public :: select_ptcl_subset, write_ptcl_subset

contains

    subroutine select_ptcl_subset(project, nkeep, nspace, pgrp, score_key, state, selected, status, message, rebalance_eo)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use simple_oris,         only: oris
        use simple_ori,          only: ori
        use simple_sym,          only: sym, is_valid_pointgroup
        use simple_type_defs,    only: class_sample
        use simple_eo_partition, only: eo_partition_status, EO_COMPLETE, repartition_selected_halves
        class(sp_project), intent(inout) :: project
        integer, intent(in) :: nkeep, nspace, state
        character(len=*), intent(in) :: pgrp, score_key
        logical, allocatable, intent(out) :: selected(:)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        logical, optional, intent(in) :: rebalance_eo
        class(oris), allocatable :: directions, ranked
        type(sym) :: symmetry
        type(ori) :: orientation
        type(class_sample), allocatable :: groups(:)
        integer, allocatable :: bins(:), flags(:)
        integer :: np, nstks, particle, stack, bin, half, physical_index, nbins
        real :: euler(3), score
        status = 1
        message = ''
        np = project%os_ptcl3D%get_noris()
        allocate(selected(np), source=.false.)
        if( nkeep < 1 .or. state < 1 .or. nspace < 2 .or. mod(nspace,2) /= 0 )then
            message = 'selection requires positive nsample and state, and a positive even nspace'; return
        endif
        if( score_key /= 'corr' .and. score_key /= 'corr_cart' )then
            message = 'subset_score must be corr or corr_cart'; return
        endif
        if( .not.is_valid_pointgroup(pgrp) )then
            message = 'invalid point-group symmetry'; return
        endif
        nstks = project%get_nstks()
        if( np < 1 .or. nstks < 1 .or. project%os_ptcl2D%get_noris() /= np )then
            message = 'selection requires congruent ptcl2D and ptcl3D rows and particle stacks'; return
        endif
        if( eo_partition_status(project,message) /= EO_COMPLETE )then
            if( len_trim(message) == 0 ) message = 'selection requires complete existing particle half-set provenance'
            return
        endif
        allocate(directions,ranked)
        call directions%new(nspace, is_ptcl=.false.)
        call symmetry%new(pgrp)
        call symmetry%build_refspiral(directions)
        nbins = directions%get_noris()
        allocate(bins(nbins), flags(np), source=0)
        bins = [(bin, bin=1,nbins)]
        call ranked%new(np, is_ptcl=.false.)
        call ranked%set_all2single('state',0)
        do particle = 1,np
            if( project%os_ptcl2D%get_state(particle) <= 0 ) cycle
            if( project%os_ptcl3D%get_state(particle) /= state ) cycle
            if( project%os_ptcl3D%get_int(particle,'updatecnt') <= 0 ) cycle
            if( .not.project%os_ptcl3D%isthere(particle,score_key) ) cycle
            score = project%os_ptcl3D%get(particle,score_key)
            euler = project%os_ptcl3D%get_euler(particle)
            if( .not.ieee_is_finite(score) .or. .not.all(ieee_is_finite(euler)) ) cycle
            if( euler(2) < 0. .or. euler(2) > 180. ) cycle
            stack = project%os_ptcl3D%get_int(particle,'stkind')
            if( stack < 1 .or. stack > nstks )then
                message = 'eligible particle has invalid stkind'; exit
            endif
            if( project%os_stk%get_state(stack) <= 0 ) cycle
            physical_index = project%os_ptcl3D%get_int(particle,'indstk')
            if( project%os_ptcl2D%get_int(particle,'stkind') /= stack .or. &
                project%os_ptcl2D%get_int(particle,'indstk') /= physical_index )then
                message = 'ptcl2D and ptcl3D physical stack mappings disagree'; exit
            endif
            if( .not.project%os_stk%isthere(stack,'nptcls_stk') )then
                message = 'eligible stack lacks nptcls_stk'; exit
            endif
            if( physical_index < 1 .or. physical_index > project%os_stk%get_int(stack,'nptcls_stk') .or. &
                particle < project%os_stk%get_fromp(stack) .or. particle > project%os_stk%get_top(stack) )then
                message = 'eligible particle lies outside its physical stack or logical row range'; exit
            endif
            half = project%os_ptcl3D%get_eo(particle)
            if( half < 0 .or. half > 1 )then
                message = 'eligible particle has an invalid half label'; exit
            endif
            if( project%os_ptcl2D%get_eo(particle) /= half )then
                message = 'ptcl2D and ptcl3D half assignments disagree'; exit
            endif
            call project%os_ptcl3D%get_ori(particle,orientation)
            bin = symmetry%find_closest_proj(directions,orientation)
            call ranked%set(particle,'proj',bin)
            call ranked%set(particle,score_key,score)
            call ranked%set_state(particle,1)
        enddo
        if( len_trim(message) == 0 )then
            if( ranked%get_pop(1,'state') < nkeep )then
                message = 'nsample exceeds the number of eligible scored particles in the requested state'
            else
                call ranked%get_class_sample_stats(bins,groups,label='proj',score_key=score_key)
                call ranked%sample_balanced(groups,nkeep,.true.,flags,exact=.true.)
                selected = flags == 1
                if( present(rebalance_eo) )then
                    if( rebalance_eo ) call repartition_selected_halves(project,groups)
                endif
                status = 0
            endif
        endif
        call orientation%kill
        call directions%kill
        call ranked%kill
        call symmetry%kill
    end subroutine select_ptcl_subset

    subroutine write_ptcl_subset(project, selected, stkout, oritab)
        use simple_ori,       only: ori
        use simple_type_defs, only: ctfparams
        use simple_fileio,    only: fopen, fclose
        use simple_syslib,    only: file_exists
        use simple_error,     only: simple_exception
        class(sp_project), intent(inout) :: project
        logical, intent(in) :: selected(:)
        class(string), intent(in) :: stkout, oritab
        type(ori) :: orientation
        type(ctfparams) :: ctfvars
        integer, allocatable :: rows(:)
        integer :: particle, stack, np, box, unit
        real :: smpd
        if( size(selected) /= project%os_ptcl3D%get_noris() .or. .not.any(selected) )then
            call simple_exception('invalid particle selection mask', __FILENAME__, __LINE__)
        endif
        if( file_exists(stkout) .or. file_exists(oritab) )then
            call simple_exception('selected stack or parameter table already exists', __FILENAME__, __LINE__)
        endif
        do particle = 1,size(selected)
            if( selected(particle) ) cycle
            call project%os_ptcl2D%set_state(particle,0)
            call project%os_ptcl3D%set_state(particle,0)
        enddo
        call project%prune_particles
        np = project%os_ptcl3D%get_noris()
        if( np /= count(selected) )then
            call simple_exception('pruned project disagrees with the particle selection', __FILENAME__, __LINE__)
        endif
        box = project%get_box()
        smpd = project%get_smpd()
        do stack = 1,project%get_nstks()
            if( project%os_stk%get_int(stack,'box') /= box .or. &
                abs(project%os_stk%get(stack,'smpd')-smpd) > epsilon(smpd)*max(1.,smpd) )then
                call simple_exception('selected particle stacks must share box and sampling', __FILENAME__, __LINE__)
            endif
        enddo
        rows = [(particle, particle=1,np)]
        call project%write_substk(rows,stkout)
        do stack = 1,project%get_nstks()
            call project%os_stk%set(stack,'stk',stkout)
            call project%os_stk%set(stack,'nptcls_stk',np)
        enddo
        call fopen(unit,oritab,status='new',action='write')
        do particle = 1,np
            call project%os_ptcl2D%set(particle,'indstk',particle)
            call project%os_ptcl3D%set(particle,'indstk',particle)
            ctfvars = project%get_ctfparams('ptcl3D',particle)
            call project%os_ptcl3D%get_ori(particle,orientation)
            call orientation%set_ctfvars(ctfvars)
            call orientation%write(unit)
        enddo
        call fclose(unit)
        call orientation%kill
        call project%os_out%kill
        call project%projinfo%delete_entry('sigma2_state')
    end subroutine write_ptcl_subset

end module simple_ptcl_subset