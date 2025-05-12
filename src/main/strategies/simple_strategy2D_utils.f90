module simple_strategy2D_utils
include 'simple_lib.f08'
use simple_image,      only: image
use simple_masker,     only: density_outside_mask
use simple_stack_io,   only: stack_io
use simple_sp_project, only: sp_project
implicit none

public :: read_cavgs_into_imgarr, flag_non_junk_cavgs, write_cavgs, write_junk_cavgs, write_selected_cavgs, align_imgs2ref
private
#include "simple_local_flags.inc"

interface read_cavgs_into_imgarr
    module procedure read_cavgs_into_imgarr_1
    module procedure read_cavgs_into_imgarr_2
end interface

contains

    function read_cavgs_into_imgarr_1( spproj, mask ) result( imgs )
        class(sp_project), intent(inout) :: spproj
        logical, optional, intent(in)    :: mask(:)
        type(image),       allocatable   :: imgs(:)
        character(len=:),  allocatable   :: cavgsstk, stkpath
        type(stack_io) :: stkio_r
        integer :: icls, ncls, n, ldim_read(3), cnt, ncls_sel
        real    :: smpd
        call spproj%get_cavgs_stk(cavgsstk, ncls, smpd, imgkind='cavg', stkpath=stkpath)
        if(.not. file_exists(cavgsstk)) cavgsstk = trim(stkpath) // '/' // trim(cavgsstk)
        if(.not. file_exists(cavgsstk)) THROW_HARD('cavgs stk does not exist')
        call stkio_r%open(trim(cavgsstk), smpd, 'read', bufsz=min(1024,ncls))
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
    end function read_cavgs_into_imgarr_1

    function read_cavgs_into_imgarr_2( cavgsstk, mask ) result( imgs )
        character(len=*),  intent(in)  :: cavgsstk
        logical, optional, intent(in)  :: mask(:)
        type(image),       allocatable :: imgs(:)
        type(stack_io) :: stkio_r
        integer :: icls, ncls, n, ldim_read(3), cnt, ncls_sel
        real    :: smpd
        if(.not. file_exists(cavgsstk)) THROW_HARD('cavgs stk does not exist')
        call find_ldim_nptcls(cavgsstk, ldim_read, ncls, smpd)
        ldim_read(3) = 1
        call stkio_r%open(trim(cavgsstk), smpd, 'read', bufsz=min(1024,ncls))
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
    end function read_cavgs_into_imgarr_2

    subroutine flag_non_junk_cavgs( cavgs, lp_bin, msk, l_non_junk, os_cls2D )
        class(image),          intent(inout) :: cavgs(:)
        real,                  intent(in)    :: lp_bin, msk
        logical, allocatable,  intent(inout) :: l_non_junk(:)
        class(oris), optional, intent(in)    :: os_cls2D
        real,        parameter   :: DYNRANGE_THRES = 1e-6
        real,        parameter   :: HP_SPEC        = 20.
        real,        parameter   :: LP_SPEC        = 6.
        integer,     parameter   :: MINPOP         = 20
        type(image), allocatable :: cavg_threads(:)
        real,        allocatable :: pspec(:)
        integer :: ncls, icls, ldim(3), kfromto(2), ithr
        real    :: dynrange, smpd
        logical :: l_dens_outside, l_os2D_present
        ncls = size(cavgs)
        l_os2D_present = present(os_cls2D)
        if( l_os2D_present )then
            if( os_cls2D%get_noris() /= ncls ) THROW_HARD('# cavgs /= # entries in os_cls2D')
        endif
        ldim = cavgs(1)%get_ldim()
        smpd = cavgs(1)%get_smpd()
        if( allocated(l_non_junk) ) deallocate(l_non_junk)
        allocate(l_non_junk(ncls), source=.false.)
        kfromto(1) = calc_fourier_index(HP_SPEC, ldim(1), smpd)
        kfromto(2) = calc_fourier_index(LP_SPEC, ldim(1), smpd)
        allocate(cavg_threads(nthr_glob))
        do ithr = 1, nthr_glob
            call cavg_threads(ithr)%new(ldim, smpd)
        end do
        !$omp parallel do default(shared) private(icls,ithr,pspec,dynrange,l_dens_outside) proc_bind(close) schedule(static)
        do icls = 1, ncls
            ithr = omp_get_thread_num() + 1
            call cavg_threads(ithr)%copy(cavgs(icls))
            call cavg_threads(ithr)%norm
            l_dens_outside = density_outside_mask(cavg_threads(ithr), lp_bin, msk)
            call cavg_threads(ithr)%mask(msk, 'soft')
            call cavg_threads(ithr)%spectrum('sqrt', pspec)
            dynrange = pspec(kfromto(1)) - pspec(kfromto(2))
            if( l_os2D_present )then
                if( dynrange > DYNRANGE_THRES .and. os_cls2D%get_int(icls, 'pop') >= MINPOP )then
                    if( .not. l_dens_outside ) l_non_junk(icls) = .true.
                endif
            else
                if( dynrange > DYNRANGE_THRES .and. .not. l_dens_outside ) l_non_junk(icls) = .true.
            endif
        enddo
        !$omp end parallel do
        do ithr = 1, nthr_glob
            call cavg_threads(ithr)%kill
        end do
        deallocate(cavg_threads)
    end subroutine flag_non_junk_cavgs

    subroutine write_cavgs( n, imgs, labels, fbody, ext )
        integer,          intent(in)    :: n
        class(image),     intent(inout) :: imgs(n)
        integer,          intent(in)    :: labels(n)
        character(len=*), intent(in)    :: fbody, ext
        character(len=:), allocatable   :: fname
        integer,          allocatable   :: cnts(:)
        integer :: i, maxlab, pad_len
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
    end subroutine write_cavgs

    subroutine write_junk_cavgs( n, imgs, labels, ext )
        integer,          intent(in)    :: n
        class(image),     intent(inout) :: imgs(n)
        integer,          intent(in)    :: labels(n)
        character(len=*), intent(in)    :: ext
        character(len=:), allocatable   :: fname
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
        integer,          intent(in)    :: n
        class(image),     intent(inout) :: imgs(n)
        integer,          intent(in)    :: labels(n)
        character(len=*), intent(in)    :: ext
        character(len=:), allocatable   :: fname
        integer :: i, cnt(0:1)
        cnt = 0
        do i = 1, n
            if( labels(i) == 0 )then
                fname = 'unselected_cavgs'//trim(ext)
                cnt(0) = cnt(0) + 1
                call imgs(i)%write(fname, cnt(0))
            else
                fname  = 'selected_cavgs'//trim(ext)
                cnt(1) = cnt(1) + 1
                call imgs(i)%write(fname, cnt(1))
            endif
        end do
    end subroutine write_selected_cavgs

    subroutine align_imgs2ref( n, hp, lp, trs, imgs, img_ref, imgs_aligned )
        use simple_polarizer,         only: polarizer
        use simple_polarft_corrcalc,  only: polarft_corrcalc
        use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
        integer,      intent(in)    :: n
        real,         intent(in)    :: hp, lp, trs
        class(image), intent(inout) :: imgs(n), img_ref, imgs_aligned(n)
        integer,      parameter :: MAXITS_SH = 60
        real,           allocatable :: inpl_corrs(:)
        type(pftcc_shsrch_grad)     :: grad_shsrch_obj(nthr_glob)
        type(polarizer)             :: polartransform
        type(polarft_corrcalc)      :: pftcc
        type(image)                 :: imgs_mirr(n)
        type(inpl_struct)           :: algninfo_mirr(n), algninfo(n)
        real(kind=c_float), pointer :: rmat_ptr(:,:,:) => null()
        integer :: ldim(3), ldim_ref(3), box, kfromto(2), ithr, i, loc(1), nrots, irot
        real    :: smpd, lims(2,2), lims_init(2,2), cxy(3)
        logical :: l_mirr(n)
        ldim       = imgs(1)%get_ldim()
        ldim_ref   = img_ref%get_ldim()
        if( .not. all(ldim == ldim_ref) ) THROW_HARD('Incongruent logical image dimensions (imgs & img_ref)')
        box        = ldim(1)
        smpd       = imgs(1)%get_smpd()
        kfromto(1) = max(2, calc_fourier_index(hp, box, smpd))
        kfromto(2) =        calc_fourier_index(lp, box, smpd)
        ! create mirrored versions of the images
        do i = 1, n
            call imgs_mirr(i)%new(ldim, smpd, wthreads=.false.)
        end do
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, n
            call imgs_aligned(i)%new(ldim, smpd, wthreads=.false.)
            call imgs_mirr(i)%copy(imgs(i))
            call imgs_mirr(i)%mirror('x')
        end do
        !$omp end parallel do
        ! initialize pftcc, polarizer
        call pftcc%new(1, [1,2*n], kfromto) ! 2*n because of mirroring
        call polartransform%new([box,box,1], smpd)
        call polartransform%init_polarizer(pftcc, KBALPHA)
        ! ! in-plane search object objects for parallel execution
        lims(:,1)      = -trs
        lims(:,2)      =  trs
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%new(lims, lims_init=lims_init, shbarrier='yes',&
            &maxits=MAXITS_SH, opt_angle=.true.)
        end do
        ! set the reference transform
        call polartransform%polarize(pftcc, img_ref, 1, isptcl=.false., iseven=.true.)
        ! set the particle transforms
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, 2 * n
            if( i <= n )then
                call imgs(i)%fft()
                call polartransform%polarize(pftcc, imgs(i),        i, isptcl=.true., iseven=.true.)
            else
                call imgs_mirr(i-n)%fft()
                call polartransform%polarize(pftcc, imgs_mirr(i-n), i, isptcl=.true., iseven=.true.)
            endif
        end do
        !$omp end parallel do
        ! register imgs to img_ref
        nrots = pftcc%get_nrots()
        allocate(inpl_corrs(nrots), source=0.)
        !$omp parallel do default(shared) private(i,ithr,inpl_corrs,loc,irot,cxy) schedule(static) proc_bind(close)
        do i = 1, 2 * n
            ithr = omp_get_thread_num() + 1
            call pftcc%gencorrs(1, i, inpl_corrs)
            loc  = maxloc(inpl_corrs)
            irot = loc(1) 
            call grad_shsrch_obj(ithr)%set_indices(1, i)
            cxy = grad_shsrch_obj(ithr)%minimize(irot=irot)
            if( irot == 0 )then ! no improved solution found, put back the old one
                cxy(1) = inpl_corrs(loc(1))
                cxy(2) = 0.
                cxy(3) = 0.
                irot   = loc(1)
            endif
            if( i <= n )then
                algninfo(i)%e3            = 360. - pftcc%get_rot(irot)
                algninfo(i)%corr          = cxy(1)
                algninfo(i)%x             = cxy(2)
                algninfo(i)%y             = cxy(3)
                algninfo(i)%l_mirr        = .false.
            else
                algninfo_mirr(i-n)%e3     = 360. - pftcc%get_rot(irot)
                algninfo_mirr(i-n)%corr   = cxy(1)
                algninfo_mirr(i-n)%x      = cxy(2)
                algninfo_mirr(i-n)%y      = cxy(3)
                algninfo_mirr(i-n)%l_mirr = .true.
            endif
        end do
        !$omp end parallel do
        ! set mirror flags
        where( algninfo_mirr(:)%corr > algninfo(:)%corr ) algninfo = algninfo_mirr
        ! shift and rotate the images
        !$omp parallel do default(shared) private(i,rmat_ptr) schedule(static) proc_bind(close)
        do i = 1, n
            call imgs_aligned(i)%get_rmat_ptr(rmat_ptr)
            if( algninfo(i)%l_mirr )then
                call imgs_mirr(i)%shift2Dserial([-algninfo(i)%x,-algninfo(i)%y])
                call imgs_mirr(i)%ifft
                call imgs_mirr(i)%rtsq_serial(algninfo(i)%e3, 0., 0., rmat_ptr)
            else
                call imgs(i)%shift2Dserial([-algninfo(i)%x,-algninfo(i)%y])
                call imgs(i)%ifft
                call imgs(i)%rtsq_serial(algninfo(i)%e3, 0., 0., rmat_ptr)
            endif
        end do
        ! destruct
        do i = 1, n
            call imgs_mirr(i)%kill
        end do
    end subroutine align_imgs2ref

end module simple_strategy2D_utils
