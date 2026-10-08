!@descr: particle-layout identity: the digest that ties per-particle side files to a project's particle rows
module simple_ptcl_layout
use, intrinsic :: iso_fortran_env, only: int64
use simple_oris,              only: oris
use simple_sp_project,        only: sp_project
use simple_string,            only: string
use simple_sigma2_state_file, only: sigma2_state_digest_begin, sigma2_state_digest_text, sigma2_state_digest_integer
implicit none
private

public :: ptcl_layout_digest

!> Digest of (lineage, stack reference and index of every row); 0 when the layout is undefined
interface ptcl_layout_digest
    module procedure ptcl_layout_digest_project
    module procedure ptcl_layout_digest_arrays
end interface ptcl_layout_digest

contains

    function ptcl_layout_digest_arrays(lineage, stack_refs, stack_ids, stack_indices, nrows) result(digest)
        character(len=*),  intent(in) :: lineage
        type(string),      intent(in) :: stack_refs(:)
        integer,           intent(in) :: stack_ids(:), stack_indices(:)
        integer, optional, intent(in) :: nrows
        integer(int64) :: digest
        integer :: i, n
        n = size(stack_ids)
        if( present(nrows) ) n = nrows
        if( n < 0 .or. n > size(stack_ids) .or. n > size(stack_indices) )then
            digest = 0_int64
            return
        endif
        digest = sigma2_state_digest_begin()
        call sigma2_state_digest_text(digest, trim(lineage))
        do i = 1, n
            if( stack_ids(i) < 1 .or. stack_ids(i) > size(stack_refs) )then
                digest = 0_int64
                return
            endif
            call sigma2_state_digest_text(digest, trim(stack_refs(stack_ids(i))%to_char()))
            call sigma2_state_digest_integer(digest, stack_indices(i))
        enddo
        if( digest == 0_int64 ) digest = 1_int64
    end function ptcl_layout_digest_arrays

    !> Layout of the first nrows (default all) rows of particles; the lineage is the project name
    function ptcl_layout_digest_project(project, particles, nrows) result(digest)
        type(sp_project),  intent(in) :: project
        class(oris),       intent(in) :: particles
        integer, optional, intent(in) :: nrows
        integer(int64) :: digest
        type(string), allocatable :: stack_refs(:)
        integer,      allocatable :: stack_ids(:), stack_indices(:)
        type(string) :: lineage, stack_ref
        integer :: i, nptcls, nstks
        digest = 0_int64
        nptcls = particles%get_noris(consider_state=.false.)
        nstks  = project%os_stk%get_noris(consider_state=.false.)
        if( nptcls < 1 .or. nstks < 1 ) return
        if( project%projinfo%get_noris() /= 1 ) return
        if( project%projinfo%isthere(1, 'projname') )then
            lineage = project%projinfo%get_str(1, 'projname')
        else if( project%projinfo%isthere(1, 'projfile') )then
            lineage = project%projinfo%get_str(1, 'projfile')
        else
            return
        endif
        allocate(stack_refs(nstks), stack_ids(nptcls), stack_indices(nptcls))
        do i = 1, nstks
            stack_ref = project%os_stk%get_str(i, 'stk')
            stack_refs(i) = trim(adjustl(stack_ref%to_char()))
        enddo
        do i = 1, nptcls
            stack_ids(i)     = particles%get_int(i, 'stkind')
            stack_indices(i) = particles%get_int(i, 'indstk')
        enddo
        if( present(nrows) )then
            digest = ptcl_layout_digest_arrays(lineage%to_char(), stack_refs, stack_ids, stack_indices, nrows)
        else
            digest = ptcl_layout_digest_arrays(lineage%to_char(), stack_refs, stack_ids, stack_indices)
        endif
        call lineage%kill
        call stack_ref%kill
        call stack_refs(:)%kill
        deallocate(stack_refs, stack_ids, stack_indices)
    end function ptcl_layout_digest_project

end module simple_ptcl_layout
