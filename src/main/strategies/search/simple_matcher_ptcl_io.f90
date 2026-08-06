!@descr: particle image batch I/O routines shared by matcher workflows
module simple_matcher_ptcl_io
use simple_pftc_srch_api
use simple_builder,           only: builder
use simple_discrete_stack_io, only: dstack_io
use simple_imghead,           only: find_ldim_nptcls
use simple_syslib,            only: io_read_nstreams
implicit none
#include "simple_local_flags.inc"

public :: prepimgbatch, killimgbatch, read_imgbatch, discrete_read_imgbatch, discrete_read_imgbatch_source
private

interface read_imgbatch
    module procedure read_imgbatch_1
    module procedure read_imgbatch_2
    module procedure read_imgbatch_3
end interface read_imgbatch

type(stack_io) :: stkio_r

! Per-stack dimensions, memoized by stack index for the lifetime of the process.
! The physical MRC header stays the source of truth -- os_stk's 'box' is never
! repaired and 'nptcls_stk' may hold a project range count rather than the real
! image count, so trusting either would risk wrong record offsets after a stack is
! replaced. Memoizing just removes the repeated open+header read per batch.
integer, allocatable :: memo_ldim(:,:), memo_nptcls(:)

contains

    !>  Logical dimensions and image count of a particle stack, read from the file once
    !!  per stack per process. stkind < 1 means the caller has no os_stk row to key on
    !!  (cls3D, alternate particle sources) and the header is read every time.
    subroutine stk_dims( params, build, stkind, stkname, ldim, nptcls )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: stkind
        class(string),     intent(in)    :: stkname
        integer,           intent(out)   :: ldim(3), nptcls
        integer :: nstks
        if( stkind < 1 )then
            call find_ldim_nptcls(stkname, ldim, nptcls)
        else
            nstks = build%spproj%os_stk%get_noris()
            if( allocated(memo_nptcls) )then
                if( size(memo_nptcls) /= nstks ) call forget_stk_dims
            endif
            if( .not. allocated(memo_nptcls) )then
                allocate(memo_ldim(3,nstks), source=0)
                allocate(memo_nptcls(nstks), source=0)
            endif
            if( stkind > nstks )then
                call find_ldim_nptcls(stkname, ldim, nptcls)
            else
                if( memo_nptcls(stkind) < 1 )then
                    call find_ldim_nptcls(stkname, memo_ldim(:,stkind), memo_nptcls(stkind))
                endif
                ldim   = memo_ldim(:,stkind)
                nptcls = memo_nptcls(stkind)
            endif
        endif
        if( (ldim(1) /= params%box) .or. (ldim(2) /= params%box) )then
            write(logfhandle,*) 'ldim ', ldim
            write(logfhandle,*) 'box ', params%box
            write(logfhandle,*) 'stkname ', stkname%to_char()
            THROW_HARD('Incompatible dimensions! stk_dims')
        endif
    end subroutine stk_dims

    subroutine forget_stk_dims
        if( allocated(memo_ldim)   ) deallocate(memo_ldim)
        if( allocated(memo_nptcls) ) deallocate(memo_nptcls)
    end subroutine forget_stk_dims

    subroutine prepimgbatch( params, build, batchsz, box, smpd )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: batchsz
        integer, optional, intent(in)    :: box
        real,    optional, intent(in)    :: smpd
        integer :: currsz, ibatch, box_here, ldim_curr(3)
        real    :: smpd_here
        logical :: doprep
        box_here  = params%box
        smpd_here = params%smpd
        if( present(box)  ) box_here  = box
        if( present(smpd) ) smpd_here = smpd
        if( .not. allocated(build%imgbatch) )then
            doprep = .true.
        else
            currsz    = size(build%imgbatch)
            ldim_curr = build%imgbatch(1)%get_ldim()
            ! a batch read from the downscaled cache needs box_crop buffers at
            ! smpd_crop, so both have to be honoured on reuse, not just the batch size
            if( batchsz > currsz .or. ldim_curr(1) /= box_here .or. &
               &abs(build%imgbatch(1)%get_smpd() - smpd_here) > 1.e-6 )then
                call killimgbatch(build)
                doprep = .true.
            else
                doprep = .false.
            endif
        endif
        if( doprep )then
            allocate(build%imgbatch(batchsz))
            !$omp parallel do default(shared) private(ibatch) schedule(static) proc_bind(close)
            do ibatch = 1,batchsz
                call build%imgbatch(ibatch)%new([box_here, box_here, 1], smpd_here, wthreads=.false.)
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

    subroutine discrete_read_imgbatch_source( params, build, source, n, pinds, batchlims, imgs )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        character(len=*),  intent(in)    :: source
        integer,           intent(in)    :: n, pinds(n), batchlims(2)
        class(image),      intent(inout) :: imgs(:)
        type(dstack_io), allocatable :: dstkios(:)
        type(string), allocatable :: stknames(:), uniq_stknames(:)
        integer,      allocatable :: inds_in_stk(:), stk_ids(:), uniq_ldims(:,:), uniq_nptcls(:)
        integer :: i, ii, istk, nbatch, nstks, nthr_read, stk_from, stk_to, iopen, nopen
        logical :: l_known_stack, l_verbose
        if( batchlims(1) < 1 .or. batchlims(2) > n .or. batchlims(1) > batchlims(2) )then
            write(logfhandle,*) 'batchlims: ', batchlims
            write(logfhandle,*) 'n        : ', n
            THROW_HARD('invalid batchlims; discrete_read_imgbatch_source')
        endif
        nbatch = batchlims(2) - batchlims(1) + 1
        if( size(imgs) < nbatch )then
            write(logfhandle,*) 'size(imgs): ', size(imgs)
            write(logfhandle,*) 'nbatch    : ', nbatch
            THROW_HARD('image buffer too small; discrete_read_imgbatch_source')
        endif
        l_verbose = params%prg .eq. 'cls_split'
        if( l_verbose )then
            write(logfhandle,'(A,A,A,I8,A,I8,A,I8,A,I8)') 'discrete_read_imgbatch_source(', trim(source), '): n=', n, &
                &' pinds[min,max]=', minval(pinds(batchlims(1):batchlims(2))), &
                &' ', maxval(pinds(batchlims(1):batchlims(2))), ' batch_to=', batchlims(2)
            call flush(logfhandle)
        endif
        allocate(stknames(nbatch), inds_in_stk(nbatch), stk_ids(nbatch))
        allocate(uniq_stknames(nbatch), uniq_ldims(3,nbatch), uniq_nptcls(nbatch))
        nstks = 0
        do i=batchlims(1),batchlims(2)
            ii = i - batchlims(1) + 1
            call build%spproj%get_ptcl_source_stkname_and_ind(params%oritype, pinds(i), source, stknames(ii), inds_in_stk(ii))
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
                write(logfhandle,'(A,I8,A,I8,A,I8,A,A)') 'discrete_read_imgbatch_source item: batch_i=', i, &
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
                THROW_HARD('Incompatible dimensions! discrete_read_imgbatch_source')
            endif
        enddo
        do ii = 1,nbatch
            if( inds_in_stk(ii) > uniq_nptcls(stk_ids(ii)) )then
                write(logfhandle,*) 'iptcl/n/image index: ', pinds(batchlims(1)+ii-1), uniq_nptcls(stk_ids(ii)), inds_in_stk(ii)
                write(logfhandle,*) 'stkname           : ', stknames(ii)%to_char()
                THROW_HARD('particle indstk out of source stack range; discrete_read_imgbatch_source')
            endif
        enddo
        nthr_read = io_read_nstreams(uniq_stknames(1), nstks)
        allocate(dstkios(nthr_read))
        ! A batch can touch hundreds of stacks; keep simultaneous file handles bounded.
        do stk_from = 1,nstks,nthr_read
            stk_to = min(stk_from + nthr_read - 1, nstks)
            nopen  = stk_to - stk_from + 1
            do iopen = 1,nopen
                istk = stk_from + iopen - 1
                call dstkios(iopen)%new(params%smpd, params%box)
                call dstkios(iopen)%cache_stack_info(uniq_stknames(istk), uniq_ldims(:,istk), uniq_nptcls(istk))
                call dstkios(iopen)%open(uniq_stknames(istk))
            enddo
            !$omp parallel do default(shared) private(iopen,istk,ii) schedule(static) proc_bind(close) &
            !$omp& num_threads(nopen) if(nopen > 1)
            do iopen = 1,nopen
                istk = stk_from + iopen - 1
                do ii = 1,nbatch
                    if( stk_ids(ii) == istk ) call dstkios(iopen)%read(stknames(ii), inds_in_stk(ii), imgs(ii))
                enddo
            enddo
            !$omp end parallel do
            do iopen = 1,nopen
                call dstkios(iopen)%kill
            enddo
        enddo
        call stknames(:)%kill
        call uniq_stknames(:)%kill
        deallocate(dstkios, stknames, inds_in_stk, stk_ids, uniq_stknames, uniq_ldims, uniq_nptcls)
    end subroutine discrete_read_imgbatch_source

    subroutine discrete_read_imgbatch( params, build, n, pinds, batchlims )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: n, pinds(n), batchlims(2)
        type(dstack_io), allocatable :: dstkios(:)
        type(string), allocatable :: stknames(:), uniq_stknames(:)
        integer,      allocatable :: inds_in_stk(:), stk_ids(:), uniq_ldims(:,:), uniq_nptcls(:), uniq_stkinds(:)
        integer :: i, ii, istk, nbatch, nstks, nthr_read, stk_from, stk_to, iopen, nopen, stkind, ind_in_stk
        logical :: l_known_stack, l_verbose, l_cls3D
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
        if( read_sorted_stack_runs() ) return
        allocate(stknames(nbatch), inds_in_stk(nbatch), stk_ids(nbatch))
        allocate(uniq_stknames(nbatch), uniq_ldims(3,nbatch), uniq_nptcls(nbatch), uniq_stkinds(nbatch))
        l_cls3D = trim(params%oritype) == 'cls3D'
        nstks   = 0
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
                ! os_stk row backing this stack, so its dimensions can be read from
                ! the project rather than from the file header; cls3D has no such row
                uniq_stkinds(nstks) = 0
                if( .not. l_cls3D )then
                    call build%spproj%map_ptcl_ind2stk_ind(params%oritype, pinds(i), stkind, ind_in_stk)
                    uniq_stkinds(nstks) = stkind
                endif
            endif
            if( l_verbose .and. n <= 32 )then
                write(logfhandle,'(A,I8,A,I8,A,I8,A,A)') 'discrete_read_imgbatch item: batch_i=', i, &
                    &' iptcl=', pinds(i), ' indstk=', inds_in_stk(ii), ' stk=', trim(stknames(ii)%to_char())
                call flush(logfhandle)
            endif
        end do
        do istk = 1,nstks
            call stk_dims(params, build, uniq_stkinds(istk), uniq_stknames(istk), &
                &uniq_ldims(:,istk), uniq_nptcls(istk))
        enddo
        nthr_read = io_read_nstreams(uniq_stknames(1), nstks)
        allocate(dstkios(nthr_read))
        ! A batch can touch hundreds of stacks; keep simultaneous file handles bounded.
        do stk_from = 1,nstks,nthr_read
            stk_to = min(stk_from + nthr_read - 1, nstks)
            nopen  = stk_to - stk_from + 1
            do iopen = 1,nopen
                istk = stk_from + iopen - 1
                call dstkios(iopen)%new(params%smpd, params%box)
                call dstkios(iopen)%cache_stack_info(uniq_stknames(istk), uniq_ldims(:,istk), uniq_nptcls(istk))
                call dstkios(iopen)%open(uniq_stknames(istk))
            enddo
            !$omp parallel do default(shared) private(iopen,istk,ii) schedule(static) proc_bind(close) &
            !$omp& num_threads(nopen) if(nopen > 1)
            do iopen = 1,nopen
                istk = stk_from + iopen - 1
                do ii = 1,nbatch
                    if( stk_ids(ii) == istk ) call dstkios(iopen)%read(stknames(ii), inds_in_stk(ii), build%imgbatch(ii))
                enddo
            enddo
            !$omp end parallel do
            do iopen = 1,nopen
                call dstkios(iopen)%kill
            enddo
        enddo
        call stknames(:)%kill
        call uniq_stknames(:)%kill
        deallocate(dstkios, stknames, inds_in_stk, stk_ids, uniq_stknames, uniq_ldims, uniq_nptcls, uniq_stkinds)

    contains

        logical function read_sorted_stack_runs() result(handled)
            type(dstack_io), allocatable :: dstkios_run(:)
            type(string),    allocatable :: run_stknames(:)
            integer,         allocatable :: run_stkinds(:), run_from(:), run_to(:), run_ldims(:,:), run_nptcls(:)
            integer,         allocatable :: inds_in_stk_run(:)
            integer :: i, ii, irun, stkind, prev_stkind, ind_in_stk, nthr_run, run_from_i, run_to_i, iopen, nopen
            handled = .false.
            if( trim(params%oritype) == 'cls3D' ) return
            do i = batchlims(1)+1, batchlims(2)
                if( pinds(i) < pinds(i-1) ) return
            end do
            allocate(inds_in_stk_run(nbatch), run_stkinds(nbatch), run_from(nbatch), run_to(nbatch))
            nstks = 0
            prev_stkind = -1
            do i = batchlims(1), batchlims(2)
                ii = i - batchlims(1) + 1
                call build%spproj%map_ptcl_ind2stk_ind(params%oritype, pinds(i), stkind, ind_in_stk)
                if( stkind < prev_stkind )then
                    deallocate(inds_in_stk_run, run_stkinds, run_from, run_to)
                    return
                endif
                inds_in_stk_run(ii) = ind_in_stk
                if( stkind /= prev_stkind )then
                    nstks = nstks + 1
                    run_stkinds(nstks) = stkind
                    run_from(nstks)    = ii
                    if( nstks > 1 ) run_to(nstks-1) = ii - 1
                    prev_stkind = stkind
                endif
            end do
            if( nstks < 1 )then
                deallocate(inds_in_stk_run, run_stkinds, run_from, run_to)
                handled = .true.
                return
            endif
            run_to(nstks) = nbatch
            allocate(run_stknames(nstks), run_ldims(3,nstks), run_nptcls(nstks))
            do irun = 1, nstks
                run_stknames(irun) = build%spproj%os_stk%get_str(run_stkinds(irun), 'stk')
                call stk_dims(params, build, run_stkinds(irun), run_stknames(irun), &
                    &run_ldims(:,irun), run_nptcls(irun))
            end do
            nthr_run = io_read_nstreams(run_stknames(1), nstks)
            allocate(dstkios_run(nthr_run))
            do run_from_i = 1,nstks,nthr_run
                run_to_i = min(run_from_i + nthr_run - 1, nstks)
                nopen    = run_to_i - run_from_i + 1
                do iopen = 1,nopen
                    irun = run_from_i + iopen - 1
                    call dstkios_run(iopen)%new(params%smpd, params%box)
                    call dstkios_run(iopen)%cache_stack_info(run_stknames(irun), run_ldims(:,irun), run_nptcls(irun))
                    call dstkios_run(iopen)%open(run_stknames(irun))
                end do
                !$omp parallel do default(shared) private(iopen,irun,ii) schedule(static) proc_bind(close) &
                !$omp& num_threads(nopen) if(nopen > 1)
                do iopen = 1,nopen
                    irun = run_from_i + iopen - 1
                    do ii = run_from(irun), run_to(irun)
                        call dstkios_run(iopen)%read(run_stknames(irun), inds_in_stk_run(ii), build%imgbatch(ii))
                    end do
                end do
                !$omp end parallel do
                do iopen = 1,nopen
                    call dstkios_run(iopen)%kill
                end do
            end do
            call run_stknames(:)%kill
            deallocate(dstkios_run, run_stknames, run_stkinds, run_from, run_to, run_ldims, run_nptcls, inds_in_stk_run)
            handled = .true.
        end function read_sorted_stack_runs
    end subroutine discrete_read_imgbatch

end module simple_matcher_ptcl_io
