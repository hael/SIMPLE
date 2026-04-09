!@descr: particle image batch I/O routines shared by matcher workflows
module simple_matcher_ptcl_io
use simple_pftc_srch_api
use simple_builder,           only: builder
use simple_discrete_stack_io, only: dstack_io
implicit none

public :: prepimgbatch, killimgbatch, read_imgbatch, discrete_read_imgbatch
private

interface read_imgbatch
    module procedure read_imgbatch_1
    module procedure read_imgbatch_2
    module procedure read_imgbatch_3
end interface read_imgbatch

type(stack_io) :: stkio_r

contains

    subroutine prepimgbatch( params, build, batchsz, box )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: batchsz
        integer, optional, intent(in)    :: box
        integer :: currsz, ibatch, box_here
        logical :: doprep
        if( .not. allocated(build%imgbatch) )then
            doprep = .true.
        else
            currsz = size(build%imgbatch)
            if( batchsz > currsz )then
                call killimgbatch(build)
                doprep = .true.
            else
                doprep = .false.
            endif
        endif
        if( doprep )then
            box_here = params%box
            if( present(box) ) box_here = box
            allocate(build%imgbatch(batchsz))
            !$omp parallel do default(shared) private(ibatch) schedule(static) proc_bind(close)
            do ibatch = 1,batchsz
                call build%imgbatch(ibatch)%new([box_here, box_here, 1], params%smpd, wthreads=.false.)
            end do
            !$omp end parallel do
        endif
    end subroutine prepimgbatch

    subroutine killimgbatch( build )
        type(builder), intent(inout) :: build
        integer :: ibatch
        if( allocated(build%imgbatch) )then
            do ibatch=1,size(build%imgbatch)
                call build%imgbatch(ibatch)%kill
            end do
            deallocate(build%imgbatch)
        endif
    end subroutine killimgbatch

    subroutine read_imgbatch_1( params, build, fromptop )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: fromptop(2)
        type(string) :: stkname
        integer :: iptcl, ind_in_batch, ind_in_stk
        do iptcl=fromptop(1),fromptop(2)
            ind_in_batch = iptcl - fromptop(1) + 1
            call build%spproj%get_stkname_and_ind(params%oritype, iptcl, stkname, ind_in_stk)
            if( .not. stkio_r%stk_is_open() )then
                call stkio_r%open(stkname, params%smpd, 'read')
            else if( .not. stkio_r%same_stk(stkname, [params%box,params%box,1]) )then
                call stkio_r%close
                call stkio_r%open(stkname, params%smpd, 'read')
            endif
            call stkio_r%read(ind_in_stk, build%imgbatch(ind_in_batch))
        end do
        call stkio_r%close
    end subroutine read_imgbatch_1

    subroutine read_imgbatch_2( params, build, n, pinds, batchlims )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: n, pinds(n), batchlims(2)
        type(string) :: stkname
        integer :: ind_in_stk, i, ii
        do i=batchlims(1),batchlims(2)
            ii = i - batchlims(1) + 1
            call build%spproj%get_stkname_and_ind(params%oritype, pinds(i), stkname, ind_in_stk)
            if( .not. stkio_r%stk_is_open() )then
                call stkio_r%open(stkname, params%smpd, 'read')
            else if( .not. stkio_r%same_stk(stkname, [params%box,params%box,1]) )then
                call stkio_r%close
                call stkio_r%open(stkname, params%smpd, 'read')
            endif
            call stkio_r%read(ind_in_stk, build%imgbatch(ii))
        end do
        call stkio_r%close
    end subroutine read_imgbatch_2

    subroutine read_imgbatch_3( params, build, iptcl, img )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: iptcl
        class(image),      intent(inout) :: img
        type(string) :: stkname
        integer :: ind_in_stk
        call build%spproj%get_stkname_and_ind(params%oritype, iptcl, stkname, ind_in_stk)
        if( .not. stkio_r%stk_is_open() )then
            call stkio_r%open(stkname, params%smpd, 'read')
        else if( .not. stkio_r%same_stk(stkname, [params%box,params%box,1]) )then
            call stkio_r%close
            call stkio_r%open(stkname, params%smpd, 'read')
        endif
        call stkio_r%read(ind_in_stk, img)
        call stkio_r%close
    end subroutine read_imgbatch_3

    subroutine discrete_read_imgbatch( params, build, n, pinds, batchlims )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: n, pinds(n), batchlims(2)
        type(dstack_io) :: dstkio_r
        type(string) :: stkname
        integer :: ind_in_stk, i, ii
        call dstkio_r%new(params%smpd, params%box)
        do i=batchlims(1),batchlims(2)
            ii = i - batchlims(1) + 1
            call build%spproj%get_stkname_and_ind(params%oritype, pinds(i), stkname, ind_in_stk)
            call dstkio_r%read(stkname, ind_in_stk, build%imgbatch(ii))
        end do
        call dstkio_r%kill
    end subroutine discrete_read_imgbatch

end module simple_matcher_ptcl_io