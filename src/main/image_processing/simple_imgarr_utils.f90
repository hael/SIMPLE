!@descr: utilities for image arrays
module simple_imgarr_utils
use simple_core_module_api
use simple_image,      only: image
use simple_sp_project, only: sp_project
use simple_stack_io,   only: stack_io
implicit none
#include "simple_local_flags.inc"

interface write_imgarr
    module procedure write_imgarr_1
    module procedure write_imgarr_2
    module procedure write_imgarr_3
end interface

contains

    subroutine alloc_imgarr( n, ldim, smpd, imgs, wthreads )
        integer,                  intent(in)    :: n, ldim(3)
        real,                     intent(in)    :: smpd
        type(image), allocatable, intent(inout) :: imgs(:)
        logical,        optional, intent(in)    :: wthreads
        integer :: i
        logical :: with_threads
        with_threads = .false.
        if( present(wthreads) ) with_threads = wthreads
        if( allocated(imgs) ) call dealloc_imgarr(imgs)
        allocate(imgs(n))
        !$omp parallel do schedule(static) proc_bind(close) private(i) default(shared)
        do i = 1, n
            call imgs(i)%new(ldim, smpd, wthreads=with_threads)
        end do
        !$omp end parallel do
    end subroutine alloc_imgarr

    function pack_imgarr( imgs, mask ) result( imgs_packed )
        class(image), intent(in) :: imgs(:)
        logical,      intent(in) :: mask(:)
        type(image), allocatable :: imgs_packed(:)
        integer :: n, cnt, n_pack, i
        n = size(imgs)
        if( n /= size(mask) ) THROW_HARD('Incongruent mask: '//int2str(n)//' vs '//int2str(size(mask)))
        n_pack = count(mask)
        if( n_pack == 0 ) return
        allocate(imgs_packed(n_pack))
        cnt = 0
        do i = 1, n
            if( mask(i) )then
                cnt = cnt + 1
                call imgs_packed(cnt)%copy(imgs(i))
            endif
        end do
    end function pack_imgarr

    function extract_imgarr( imgs, inds ) result( imgs_extr )
        class(image), intent(in) :: imgs(:)
        integer,      intent(in) :: inds(:)
        type(image), allocatable :: imgs_extr(:)
        integer :: n, n_extr, i
        n      = size(imgs)
        n_extr = size(inds)
        allocate(imgs_extr(n_extr))
        do i = 1, n_extr
            call imgs_extr(i)%copy(imgs(inds(i)))
        end do
    end function extract_imgarr

    function copy_imgarr( imgarr_in ) result( imgarr_copy )
        class(image), intent(in) :: imgarr_in(:)
        type(image), allocatable :: imgarr_copy(:)
        integer :: n, i 
        n = size(imgarr_in)
        allocate(imgarr_copy(n))
        !$omp parallel do schedule(static) proc_bind(close) private(i) default(shared)
        do i = 1, n
            call imgarr_copy(i)%copy(imgarr_in(i))
        end do
        !$omp end parallel do
    end function copy_imgarr

    function join_imgarrs( imgs1, imgs2 ) result( imgs_joined )
        class(image), intent(in) :: imgs1(:), imgs2(:)
        type(image), allocatable :: imgs_joined(:)
        integer :: n1, n2, n, cnt, i
        n1 = size(imgs1)
        n2 = size(imgs2)
        n  = n1 + n2
        call alloc_imgarr(n, imgs1(1)%get_ldim(), imgs1(1)%get_smpd(), imgs_joined)
        cnt = 0
        do i = 1, n
            if( i <= n1 )then
                call imgs_joined(i)%copy(imgs1(i))
            else
                cnt = cnt + 1
                call imgs_joined(i)%copy(imgs2(cnt))
            endif
        end do
    end function join_imgarrs

    subroutine dealloc_imgarr( imgs )
        type(image), allocatable, intent(inout) :: imgs(:)
        integer :: n , i
        if( allocated(imgs) )then
            n = size(imgs)
            !$omp parallel do schedule(static) proc_bind(close) private(i) default(shared)
            do i = 1, n
                call imgs(i)%kill
            end do
            !$omp end parallel do
            deallocate(imgs)
        endif
    end subroutine dealloc_imgarr

    function read_cavgs_into_imgarr( spproj, mask ) result( imgs )
        class(sp_project), intent(inout) :: spproj
        logical, optional, intent(in)    :: mask(:)
        type(image),       allocatable   :: imgs(:)
        type(string)   :: cavgsstk
        type(stack_io) :: stkio_r
        integer :: icls, ncls, ldim_read(3), cnt, ncls_sel
        real    :: smpd
        call spproj%get_cavgs_stk(cavgsstk, ncls, smpd, imgkind='cavg')
        if(.not. file_exists(cavgsstk)) THROW_HARD('cavgs stk does not exist')
        call stkio_r%open(cavgsstk, smpd, 'read', bufsz=min(1024,ncls))
        ldim_read    = stkio_r%get_ldim()
        ldim_read(3) = 1
        if( present(mask) )then
            if( size(mask) /= ncls ) THROW_HARD('Nonconforming mask size')
            ncls_sel = count(mask)
            allocate(imgs(ncls_sel))
            cnt = 0
            do icls = 1,ncls
                if( mask(icls) )then
                    cnt = cnt + 1
                    call imgs(cnt)%new(ldim_read,smpd,wthreads=.false.)
                    call stkio_r%read(icls, imgs(cnt))
                endif
            end do
        else
            allocate(imgs(ncls))
            do icls = 1,ncls
                call imgs(icls)%new(ldim_read,smpd,wthreads=.false.)
                call stkio_r%read(icls, imgs(icls))
            end do
        endif
        call stkio_r%close
    end function read_cavgs_into_imgarr

    function read_stk_into_imgarr( stkname, mask ) result( imgs )
        class(string),     intent(in)  :: stkname
        logical, optional, intent(in)  :: mask(:)
        type(image),       allocatable :: imgs(:)
        type(stack_io) :: stkio_r
        integer :: icls, ncls, ldim_read(3), cnt, ncls_sel
        real    :: smpd
        if(.not. file_exists(stkname)) THROW_HARD('stk stk does not exist')
        call find_ldim_nptcls(stkname, ldim_read, ncls)
        smpd = find_img_smpd(stkname)
        ldim_read(3) = 1
        call stkio_r%open(stkname, smpd, 'read', bufsz=min(1024,ncls))
        if( present(mask) )then
            if( size(mask) /= ncls ) THROW_HARD('Nonconforming mask size')
            ncls_sel = count(mask)
            allocate(imgs(ncls_sel))
            cnt = 0
            do icls = 1,ncls
                if( mask(icls) )then
                    cnt = cnt + 1
                    call imgs(cnt)%new(ldim_read,smpd,wthreads=.false.)
                    call stkio_r%read(icls, imgs(cnt))
                endif
            end do
        else
            allocate(imgs(ncls))
            do icls = 1,ncls
                call imgs(icls)%new(ldim_read,smpd,wthreads=.false.)
                call stkio_r%read(icls, imgs(icls))
            end do
        endif
        call stkio_r%close
    end function read_stk_into_imgarr

    subroutine write_imgarr_1( n, imgs, labels, fbody, ext )
        integer,          intent(in)    :: n
        class(image),     intent(inout) :: imgs(n)
        integer,          intent(in)    :: labels(n)
        character(len=*), intent(in)    :: fbody, ext
        integer, allocatable :: cnts(:)
        type(string) :: fname 
        integer      :: i, maxlab, pad_len
        maxlab = maxval(labels)
        allocate(cnts(maxlab), source=0)
        pad_len = 2
        if( maxlab > 99 ) pad_len = 3
        do i = 1, n
            if( labels(i) > 0 )then
                fname = trim(fbody)//int2str_pad(labels(i),pad_len)//'_cavgs'//trim(ext)
                cnts(labels(i)) = cnts(labels(i)) + 1
                call imgs(i)%write(fname, cnts(labels(i)))
            endif
        end do
        deallocate(cnts)
    end subroutine write_imgarr_1

    subroutine write_imgarr_2( imgs, fname )
        class(image),  intent(inout) :: imgs(:)
        class(string), intent(in)    :: fname
        integer :: n, i
        n = size(imgs)
        do i = 1, n
            call imgs(i)%write(fname, i)
        end do
    end subroutine write_imgarr_2

    subroutine write_imgarr_3( imgs, fname, inds )
        class(image),  intent(inout) :: imgs(:)
        class(string), intent(in)    :: fname
        integer,       intent(in)    :: inds(:)
        integer :: n, i, ni, cnt, ind
        n   = size(imgs)
        ni  = size(inds)
        cnt = 0
        do i = 1, ni
            ind = inds(i)
            if( ind < 0 .or. ind > n ) THROW_HARD('fetched index ind out of range')
            cnt = cnt + 1
            call imgs(ind)%write(fname, cnt)
        end do
    end subroutine write_imgarr_3

    subroutine write_junk_cavgs( n, imgs, labels, ext )
        integer,          intent(in)    :: n
        class(image),     intent(inout) :: imgs(n)
        integer,          intent(in)    :: labels(n)
        character(len=*), intent(in)    :: ext
        type(string) :: fname
        integer :: i, cnt
        cnt = 0
        do i = 1, n
            if( labels(i) == 0 )then
                fname = 'junk_cavgs'//trim(ext)
                cnt = cnt + 1
                call imgs(i)%write(fname, cnt)
            endif
        end do
    end subroutine write_junk_cavgs

    subroutine write_selected_cavgs( n, imgs, labels, ext )
        integer,           intent(in)    :: n
        class(image),      intent(inout) :: imgs(n)
        integer,           intent(in)    :: labels(n)
        character(len=*),  intent(in)    :: ext
        integer, allocatable :: cnt(:)
        type(string) :: fname 
        integer      :: i, maxlab
        maxlab = maxval(labels)
        allocate(cnt(0:maxlab), source=0)
        cnt = 0
        do i = 1, n
            if( labels(i) == 0 )then
                fname = 'unselected_cavgs'//trim(ext)
                cnt(0) = cnt(0) + 1
                call imgs(i)%write(fname, cnt(0))
            else
                fname  = 'rank'//int2str_pad(labels(i),2)//'_cavgs'//trim(ext)
                cnt(labels(i)) = cnt(labels(i)) + 1
                call imgs(i)%write(fname, cnt(labels(i)))
            endif
        end do
    end subroutine write_selected_cavgs

    !> The class averages of @p stk (sampling @p smpd) written to @p outstk in the order of @p flag
    !! ('res' and 'pop' ascending, 'corr' descending), with the ranks of the selected classes
    !! (state > 0) in @p classdoc_fname. Nothing is done unless @p os_cls holds one row per class
    !! average. The core of commander_rank_cavgs, which the stream pool calls directly.
    subroutine rank_cavgs_stk( os_cls, stk, smpd, flag, outstk, classdoc_fname )
        class(oris),      intent(in) :: os_cls
        class(string),    intent(in) :: stk, outstk, classdoc_fname
        real,             intent(in) :: smpd
        character(len=*), intent(in) :: flag
        type(oris)           :: clsdoc_ranked
        type(stack_io)       :: stkio_r, stkio_w
        type(image)          :: img
        integer, allocatable :: order(:)
        real,    allocatable :: vals(:), rstates(:)
        logical, allocatable :: mask(:)
        integer              :: ldim(3), ncls, icls, i, rank
        select case(trim(flag))
            case('res','corr','pop')
                ! supported
            case DEFAULT
                THROW_HARD('Unsupported cavg flag: '//trim(flag))
        end select
        call find_ldim_nptcls(stk, ldim, ncls)
        if( os_cls%get_noris() /= ncls ) return
        call img%new([ldim(1),ldim(2),1], smpd)
        call clsdoc_ranked%new(ncls, is_ptcl=.false.)
        vals    = os_cls%get_all(flag)
        order   = (/(icls,icls=1,ncls)/)
        rstates = os_cls%get_all('state')
        mask    = rstates > 0.5
        call hpsort(vals, order)
        if( trim(flag) == 'corr' ) call reverse(order)
        call stkio_r%open(stk, smpd, 'read', bufsz=ncls)
        call stkio_r%read_whole ! because need asynchronous access
        call stkio_w%open(outstk, smpd, 'write', box=ldim(1), bufsz=ncls)
        rank = 0
        do icls=1,ncls
            i = order(icls)
            if( mask(i) )then
                rank = rank + 1
                call clsdoc_ranked%set(rank, 'class', i)
                call clsdoc_ranked%set(rank, 'rank',  rank)
                select case(trim(flag))
                case('corr')
                    call clsdoc_ranked%set(rank, 'corr', os_cls%get(i, 'corr'))
                    write(logfhandle,'(a,1x,i5,1x,a,1x,i5,1x,a,1x,f6.2)') 'CLASS:', i,&
                    &'RANK:', rank ,'CORR:', os_cls%get(i, 'corr')
                case DEFAULT
                    call clsdoc_ranked%set(rank, 'pop',  os_cls%get(i, 'pop'))
                    call clsdoc_ranked%set(rank, 'res',  os_cls%get(i, 'res'))
                    call clsdoc_ranked%set(rank, 'corr', os_cls%get(i, 'corr'))
                    write(logfhandle,'(a,1x,i5,1x,a,1x,i5,1x,a,i5,1x,a,1x,f6.2)') 'CLASS:', i,&
                    &'RANK:', rank ,'POP:', os_cls%get_int(i, 'pop'),&
                    &'RES:', os_cls%get(i, 'res')
                end select
                call flush(logfhandle)
            endif
            call stkio_r%get_image(i, img)
            call stkio_w%write(icls, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call clsdoc_ranked%write(classdoc_fname)
        call clsdoc_ranked%kill
        call img%kill
    end subroutine rank_cavgs_stk

end module simple_imgarr_utils
