!@descr: particle image batch I/O routines shared by matcher workflows
module simple_matcher_ptcl_io
use simple_pftc_srch_api
use simple_builder,           only: builder
use simple_discrete_stack_io, only: dstack_io
use simple_imghead,           only: find_ldim_nptcls
implicit none
#include "simple_local_flags.inc"

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
        type(dstack_io), allocatable :: dstkios(:)
        type(string), allocatable :: stknames(:), uniq_stknames(:)
        integer,      allocatable :: inds_in_stk(:), stk_ids(:), uniq_ldims(:,:), uniq_nptcls(:)
        integer :: i, ii, istk, nbatch, nstks, nthr_read
        logical :: l_known_stack, l_verbose
        if( batchlims(1) < 1 .or. batchlims(2) > n .or. batchlims(1) > batchlims(2) )then
            write(logfhandle,*) 'batchlims: ', batchlims
            write(logfhandle,*) 'n        : ', n
            THROW_HARD('invalid batchlims; discrete_read_imgbatch')
        endif
        nbatch = batchlims(2) - batchlims(1) + 1
        l_verbose = params%prg .eq. 'cls_split'
        if( l_verbose )then
            write(logfhandle,'(A,I8,A,I8,A,I8,A,I8)') 'discrete_read_imgbatch: n=', n, &
                &' pinds[min,max]=', minval(pinds(batchlims(1):batchlims(2))), &
                &' ', maxval(pinds(batchlims(1):batchlims(2))), ' batch_to=', batchlims(2)
            call flush(logfhandle)
        endif
        allocate(stknames(nbatch), inds_in_stk(nbatch), stk_ids(nbatch))
        allocate(uniq_stknames(nbatch), uniq_ldims(3,nbatch), uniq_nptcls(nbatch))
        nstks = 0
        do i=batchlims(1),batchlims(2)
            ii = i - batchlims(1) + 1
            call build%spproj%get_stkname_and_ind(params%oritype, pinds(i), stknames(ii), inds_in_stk(ii))
            l_known_stack = .false.
            do istk = 1,nstks
                if( uniq_stknames(istk) .eq. stknames(ii) )then
                    l_known_stack = .true.
                    stk_ids(ii) = istk
                    exit
                endif
            enddo
            if( .not.l_known_stack )then
                nstks = nstks + 1
                uniq_stknames(nstks) = stknames(ii)
                stk_ids(ii) = nstks
            endif
            if( l_verbose .and. n <= 32 )then
                write(logfhandle,'(A,I8,A,I8,A,I8,A,A)') 'discrete_read_imgbatch item: batch_i=', i, &
                    &' iptcl=', pinds(i), ' indstk=', inds_in_stk(ii), ' stk=', trim(stknames(ii)%to_char())
                call flush(logfhandle)
            endif
        end do
        do istk = 1,nstks
            call find_ldim_nptcls(uniq_stknames(istk), uniq_ldims(:,istk), uniq_nptcls(istk))
            if( (uniq_ldims(1,istk) /= params%box) .or. (uniq_ldims(2,istk) /= params%box) )then
                write(logfhandle,*) 'ldim ', uniq_ldims(:,istk)
                write(logfhandle,*) 'box ', params%box
                write(logfhandle,*) 'stkname ', uniq_stknames(istk)%to_char()
                THROW_HARD('Incompatible dimensions! discrete_read_imgbatch')
            endif
        enddo
        allocate(dstkios(nstks))
        do istk = 1,nstks
            call dstkios(istk)%new(params%smpd, params%box)
            call dstkios(istk)%cache_stack_info(uniq_stknames(istk), uniq_ldims(:,istk), uniq_nptcls(istk))
        enddo
        nthr_read = min(max(1,nthr_glob), nstks)
        !$omp parallel do default(shared) private(istk,ii) schedule(static) proc_bind(close) &
        !$omp& num_threads(nthr_read) if(nthr_read > 1)
        do istk = 1,nstks
            do ii = 1,nbatch
                if( stk_ids(ii) == istk ) call dstkios(istk)%read(stknames(ii), inds_in_stk(ii), build%imgbatch(ii))
            enddo
        enddo
        !$omp end parallel do
        do istk = 1,nstks
            call dstkios(istk)%kill
        enddo
        call stknames(:)%kill
        call uniq_stknames(:)%kill
        deallocate(dstkios, stknames, inds_in_stk, stk_ids, uniq_stknames, uniq_ldims, uniq_nptcls)
    end subroutine discrete_read_imgbatch

end module simple_matcher_ptcl_io
