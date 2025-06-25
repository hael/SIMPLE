module simple_strategy2D_utils
include 'simple_lib.f08'
use simple_image,      only: image
use simple_masker,     only: density_inoutside_mask
use simple_stack_io,   only: stack_io
use simple_sp_project, only: sp_project
implicit none

public :: read_cavgs_into_imgarr, flag_non_junk_cavgs, match_imgs2ref, match_imgs2refs, rtsq_imgs, prep_cavgs4clustering
public :: pack_imgarr, alloc_imgarr, dealloc_imgarr, write_cavgs, write_junk_cavgs, write_selected_cavgs
public :: align_clusters2medoids, write_aligned_cavgs
private
#include "simple_local_flags.inc"

interface read_cavgs_into_imgarr
    module procedure read_cavgs_into_imgarr_1
    module procedure read_cavgs_into_imgarr_2
end interface

interface write_cavgs
    module procedure write_cavgs_1
    module procedure write_cavgs_2
end interface

contains

    function read_cavgs_into_imgarr_1( spproj, mask ) result( imgs )
        class(sp_project), intent(inout) :: spproj
        logical, optional, intent(in)    :: mask(:)
        type(image),       allocatable   :: imgs(:)
        character(len=:),  allocatable   :: cavgsstk, stkpath
        type(stack_io) :: stkio_r
        integer :: icls, ncls, ldim_read(3), cnt, ncls_sel
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
        integer :: icls, ncls, ldim_read(3), cnt, ncls_sel
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

    subroutine prep_cavgs4clustering( spproj, cavg_imgs, mskdiam, clspops, clsinds, l_non_junk, mm )
        use simple_class_frcs, only: class_frcs
        class(sp_project),        intent(inout) :: spproj
        type(image), allocatable, intent(inout) :: cavg_imgs(:)
        real,                     intent(in)    :: mskdiam
        integer,     allocatable, intent(inout) :: clspops(:), clsinds(:)
        logical,     allocatable, intent(inout) :: l_non_junk(:)
        real,        allocatable, intent(inout) :: mm(:,:)
        real,              parameter  :: LP_BIN = 20.
        logical,           parameter  :: DEBUG = .true.
        character(len=:), allocatable :: frcs_fname
        real,             allocatable :: frcs(:,:), filter(:)
        logical,          allocatable :: l_msk(:,:,:)
        integer,          allocatable :: states(:)
        type(image)                   :: img_msk
        type(class_frcs)              :: clsfrcs
        integer                       :: ncls, ldim(3), box, filtsz, ncls_sel, cnt, i, j
        real                          :: smpd, mskrad
        call dealloc_imgarr(cavg_imgs)
        ncls      = spproj%os_cls2D%get_noris()
        cavg_imgs = read_cavgs_into_imgarr(spproj)
        if( allocated(clspops) ) deallocate(clspops)
        clspops   = spproj%os_cls2D%get_all_asint('pop')
        smpd      = cavg_imgs(1)%get_smpd()
        ldim      = cavg_imgs(1)%get_ldim()
        box       = ldim(1)
        mskrad    = min(real(box/2) - COSMSKHALFWIDTH - 1., 0.5 * mskdiam/smpd)
        filtsz    = fdim(box) - 1
        ! get FRCs
        call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
        if( file_exists(frcs_fname) )then
            call clsfrcs%read(frcs_fname)
            filtsz = clsfrcs%get_filtsz()
        else
            THROW_HARD('FRC file: '//trim(frcs_fname)//' does not exist!')
        endif
        if( allocated(l_non_junk) ) deallocate(l_non_junk)
        call flag_non_junk_cavgs(cavg_imgs, LP_BIN, mskrad, l_non_junk, spproj%os_cls2D)
        if( DEBUG )then
            cnt = 0
            do i = 1, ncls
                if( .not. l_non_junk(i) )then
                    cnt = cnt + 1
                    call cavg_imgs(i)%write('cavgs_junk.mrc', cnt)
                endif
            enddo
        endif
        ! re-create cavg_imgs
        ncls_sel = count(l_non_junk)
        write(logfhandle,'(A,I5)') '# classes left after junk rejection ', ncls_sel
        call dealloc_imgarr(cavg_imgs)
        cavg_imgs = read_cavgs_into_imgarr(spproj, mask=l_non_junk)
        ! keep track of the original class indices
        if( allocated(clsinds) ) deallocate(clsinds)
        clsinds = pack((/(i,i=1,ncls)/), mask=l_non_junk)
        ! update class populations
        clspops = pack(clspops, mask=l_non_junk)
        ! create the stuff needed in the loop
        allocate(frcs(filtsz,ncls_sel), filter(filtsz), mm(ncls_sel,2), source=0.)
        ! prep mask
        call img_msk%new([box,box,1], smpd)
        img_msk = 1.
        call img_msk%mask(mskrad, 'hard')
        l_msk = img_msk%bin2logical()
        call img_msk%kill
        write(logfhandle,'(A)') '>>> PREPARING CLASS AVERAGES'
        !$omp parallel do default(shared) private(i,j,filter) schedule(static) proc_bind(close)
        do i = 1, ncls_sel
            j = clsinds(i)
            ! FRC-based filter
            call clsfrcs%frc_getter(j, frcs(:,i))
            if( any(frcs(:,i) > 0.143) )then
                call fsc2optlp_sub(clsfrcs%get_filtsz(), frcs(:,i), filter)
                where( filter > TINY ) filter = sqrt(filter) ! because the filter is applied to the average not the even or odd
                call cavg_imgs(i)%fft()
                call cavg_imgs(i)%apply_filter_serial(filter)
                call cavg_imgs(i)%ifft()
            endif
            ! normalization
            call cavg_imgs(i)%norm_within(l_msk)
            ! mask
            call cavg_imgs(i)%mask(mskrad, 'soft', backgr=0.)
            ! stash minmax
            mm(i,:) = cavg_imgs(i)%minmax(mskrad)
        end do
        !$omp end parallel do
        if( DEBUG )then
            do i = 1, ncls_sel
                call cavg_imgs(i)%write('cavgs_prepped.mrc', i)
            enddo
        endif
        call clsfrcs%kill
    end subroutine prep_cavgs4clustering

    function pack_imgarr( imgs, mask ) result( imgs_packed )
        class(image), intent(in) :: imgs(:)
        logical,      intent(in) :: mask(:)
        type(image), allocatable :: imgs_packed(:)
        integer :: n, cnt, n_pack, i
        n = size(imgs)
        if( n /= size(mask) ) THROW_HARD('Incongruent mask')
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

    subroutine flag_non_junk_cavgs( cavgs, lp_bin, msk, l_non_junk, os_cls2D )
        class(image),          intent(inout) :: cavgs(:)
        real,                  intent(in)    :: lp_bin, msk
        logical, allocatable,  intent(inout) :: l_non_junk(:)
        class(oris), optional, intent(in)    :: os_cls2D
        real,        parameter   :: DYNRANGE_THRES  = 1e-6
        real,        parameter   :: HP_SPEC         = 20.
        real,        parameter   :: LP_SPEC         = 6.
        real,        parameter   :: RATIO_THRESHOLD = 0.55
        real,        parameter   :: CENMSKFACTOR    = 5.0
        integer,     parameter   :: MINPOP          = 20
        type(image), allocatable :: cavg_threads(:)
        real,        allocatable :: pspec(:)
        integer,     allocatable :: states(:)
        integer :: ncls, icls, ldim(3), kfromto(2), ithr, nin, nout, nmsk
        real    :: cen(2),dynrange, smpd
        logical :: l_dens_inoutside, l_os2D_present
        ncls = size(cavgs)
        l_os2D_present = present(os_cls2D)
        if( l_os2D_present )then
            if( os_cls2D%get_noris() /= ncls ) THROW_HARD('# cavgs /= # entries in os_cls2D')
            states = os_cls2D%get_all_asint('state')
        else
            allocate(states(ncls), source=1)
        endif
        ldim = cavgs(1)%get_ldim()
        smpd = cavgs(1)%get_smpd()
        if( allocated(l_non_junk) ) deallocate(l_non_junk)
        allocate(l_non_junk(ncls), source=.false.)
        kfromto(1) = calc_fourier_index(HP_SPEC, ldim(1), smpd)
        kfromto(2) = calc_fourier_index(LP_SPEC, ldim(1), smpd)
        call alloc_imgarr(nthr_glob, ldim, smpd, cavg_threads)
        !$omp parallel do default(shared) proc_bind(close) schedule(static)&
        !$omp private(icls,ithr,pspec,dynrange,l_dens_inoutside,nin,nout,nmsk,cen)
        do icls = 1, ncls
            ithr = omp_get_thread_num() + 1
            call cavg_threads(ithr)%copy(cavgs(icls))
            call cavg_threads(ithr)%div_below(0., 10.) ! reduce influence of negative values
            call cavg_threads(ithr)%norm
            call density_inoutside_mask(cavg_threads(ithr), lp_bin, msk, nin, nout, nmsk, cen)
            l_dens_inoutside = (nout > 0) .or.&                           ! object oustide of mask
                              &(real(nin)/real(nmsk) > RATIO_THRESHOLD)   ! or object too big inside mask
            call cavg_threads(ithr)%mask(msk, 'soft', backgr=0.)
            call cavg_threads(ithr)%spectrum('sqrt', pspec)
            dynrange = pspec(kfromto(1)) - pspec(kfromto(2))
            l_non_junk(icls) = .false. ! exclusion by default
            if( states(icls) == 0 )then
                ! do nothing
            else
                if( l_os2D_present )then
                    if( dynrange > DYNRANGE_THRES .and. os_cls2D%get_int(icls, 'pop') >= MINPOP )then
                        if( .not. l_dens_inoutside ) l_non_junk(icls) = .true.
                    endif
                else
                    if( dynrange > DYNRANGE_THRES .and. .not. l_dens_inoutside ) l_non_junk(icls) = .true.
                endif
            endif
            ! center of identified connected-component
            if( .not.l_non_junk(icls) )then
                if( arg(cen) > real(floor(msk/CENMSKFACTOR)) )then
                    ! center too far from image center
                    l_non_junk(icls) = .false.
                endif
            endif
        enddo
        !$omp end parallel do
        call dealloc_imgarr(cavg_threads)
    end subroutine flag_non_junk_cavgs

    subroutine write_cavgs_1( n, imgs, labels, fbody, ext )
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
    end subroutine write_cavgs_1

    subroutine write_cavgs_2( imgs, fname )
        class(image),     intent(inout) :: imgs(:)
        character(len=*), intent(in)    :: fname
        integer :: n, i
        n = size(imgs)
        do i = 1, n
            call imgs(i)%write(fname, i)
        end do
    end subroutine write_cavgs_2

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
        integer,           intent(in)    :: n
        class(image),      intent(inout) :: imgs(n)
        integer,           intent(in)    :: labels(n)
        character(len=*),  intent(in)    :: ext
        character(len=:), allocatable    :: fname
        integer,          allocatable    :: cnt(:)
        integer :: i, maxlab
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

    function align_clusters2medoids( labels, i_medoids, cavg_imgs, hp, lp, trs, l_msk ) result( clust_info_arr )
        integer,          intent(in)    :: labels(:), i_medoids(:)
        class(image),     intent(inout) :: cavg_imgs(:)
        real,             intent(in)    :: hp, lp, trs
        logical,          intent(in)    :: l_msk(:,:,:)
        type(clust_info), allocatable   :: clust_info_arr(:)
        real,             allocatable   :: frc(:)
        type(image),      allocatable   :: cluster_imgs(:), cluster_imgs_aligned(:)
        real,             allocatable   :: resvals(:), resarr(:)
        real,             parameter     :: FRAC_BEST_CAVGS=0.5
        integer :: cnt, i, filtsz, ldim(3), iclust, nclust
        real    :: smpd, rfoo, best_res, worst_res
        write(logfhandle,'(A)') '>>> ALIGNING THE CLUSTERS OF CLASS AVERAGES TO THEIR MEDOIDS'
        filtsz = cavg_imgs(1)%get_filtsz()
        smpd   = cavg_imgs(1)%get_smpd()
        ldim   = cavg_imgs(1)%get_ldim()
        resarr = get_resarr(ldim(1), smpd)
        nclust = maxval(labels)
        allocate(frc(filtsz), clust_info_arr(nclust))
        nclust = maxval(labels)
        do iclust = 1, nclust
            clust_info_arr(iclust)%pop             = count(labels == iclust)
            cluster_imgs                           = pack_imgarr(cavg_imgs, mask=labels == iclust)
            clust_info_arr(iclust)%algninfo%params = match_imgs2ref(clust_info_arr(iclust)%pop, hp, lp, trs, cavg_imgs(i_medoids(iclust)), cluster_imgs)
            cluster_imgs_aligned                   = rtsq_imgs(clust_info_arr(iclust)%pop, clust_info_arr(iclust)%algninfo%params, cluster_imgs)
            ! estimate resolution
            cnt = 0
            call cavg_imgs(i_medoids(iclust))%fft
            allocate(resvals(clust_info_arr(iclust)%pop), source=0.)
            do i = 1, clust_info_arr(iclust)%pop
                call cluster_imgs_aligned(i)%fft
                call cavg_imgs(i_medoids(iclust))%fsc(cluster_imgs_aligned(i), frc)
                if( .not. all(frc > 0.5) )then ! excluding the medoid
                    cnt = cnt + 1
                    call get_resolution(frc, resarr, rfoo, resvals(cnt))
                endif
                call cluster_imgs_aligned(i)%ifft
            end do
            call cavg_imgs(i_medoids(iclust))%ifft
            ! calculate Euclidean distances
            clust_info_arr(iclust)%euclid = 0.
            do i = 1, clust_info_arr(iclust)%pop
                clust_info_arr(iclust)%euclid = clust_info_arr(iclust)%euclid + cavg_imgs(i_medoids(iclust))%sqeuclid(cluster_imgs_aligned(i), l_msk)
            end do
            clust_info_arr(iclust)%euclid = clust_info_arr(iclust)%euclid / real(clust_info_arr(iclust)%pop)
            ! report resolution as the average of the best agreeing 25% within a cluster
            clust_info_arr(iclust)%res = avg_frac_smallest(resvals(:cnt), FRAC_BEST_CAVGS)
            ! destruct
            call dealloc_imgarr(cluster_imgs)
            call dealloc_imgarr(cluster_imgs_aligned)
            deallocate(resvals)
        end do
        best_res  = minval(clust_info_arr(:)%res)
        worst_res = maxval(clust_info_arr(:)%res)
        where( clust_info_arr(:)%pop < 2 ) clust_info_arr(:)%res = worst_res ! nothing else makes sense
    end function align_clusters2medoids

    subroutine write_aligned_cavgs( labels, cavg_imgs, clust_info_arr, fbody, ext )
        integer,          intent(in)    :: labels(:)
        class(image),     intent(inout) :: cavg_imgs(:)
        type(clust_info), intent(in)    :: clust_info_arr(:)
        character(len=*), intent(in)    :: fbody, ext
        type(image), allocatable :: cluster_imgs(:), cluster_imgs_aligned(:)
        integer :: iclust, nclust
        write(logfhandle,'(A)') '>>> ROTATING & SHIFTING CLASS AVERAGES'
        nclust = size(clust_info_arr)
        do iclust = 1, nclust
            cluster_imgs         = pack_imgarr(cavg_imgs, mask=labels == iclust)
            cluster_imgs_aligned = rtsq_imgs(clust_info_arr(iclust)%pop, clust_info_arr(iclust)%algninfo%params, cluster_imgs)
            call write_cavgs(cluster_imgs_aligned, trim(fbody)//int2str_pad(iclust,2)//trim(ext))
            call dealloc_imgarr(cluster_imgs)
            call dealloc_imgarr(cluster_imgs_aligned)
        end do
    end subroutine write_aligned_cavgs

    function match_imgs2ref( n, hp, lp, trs, img_ref, imgs ) result( algninfo )
        use simple_polarizer,         only: polarizer
        use simple_polarft_corrcalc,  only: polarft_corrcalc
        use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
        integer,                  intent(in)    :: n
        real,                     intent(in)    :: hp, lp, trs
        class(image),             intent(inout) :: imgs(n), img_ref
        integer,     parameter          :: MAXITS_SH = 60
        real,        allocatable        :: inpl_corrs(:)
        type(image), allocatable        :: imgs_mirr(:)
        type(pftcc_shsrch_grad)         :: grad_shsrch_obj(nthr_glob)
        type(polarizer)                 :: polartransform
        type(polarft_corrcalc)          :: pftcc
        type(inpl_struct)               :: algninfo_mirr(n)
        type(inpl_struct), allocatable  :: algninfo(:)
        integer :: ldim(3), ldim_ref(3), box, kfromto(2), ithr, i, loc(1), nrots, irot
        real    :: smpd, lims(2,2), lims_init(2,2), cxy(3)
        ldim       = imgs(1)%get_ldim()
        ldim_ref   = img_ref%get_ldim()
        if( .not. all(ldim == ldim_ref) ) THROW_HARD('Incongruent logical image dimensions (imgs & img_ref)')
        box        = ldim(1)
        smpd       = imgs(1)%get_smpd()
        kfromto(1) = max(2, calc_fourier_index(hp, box, smpd))
        kfromto(2) =        calc_fourier_index(lp, box, smpd)
        ! create mirrored versions of the images
        call alloc_imgarr(n, ldim, smpd, imgs_mirr)
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, n
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
                call imgs(i)%ifft
            else
                call imgs_mirr(i-n)%fft()
                call polartransform%polarize(pftcc, imgs_mirr(i-n), i, isptcl=.true., iseven=.true.)
                call imgs_mirr(i-n)%ifft()
            endif
        end do
        !$omp end parallel do
        call pftcc%memoize_refs
        call pftcc%memoize_ptcls
        ! register imgs to img_ref
        nrots = pftcc%get_nrots()
        allocate(inpl_corrs(nrots), algninfo(n))
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
                algninfo(i)%e3            = pftcc%get_rot(irot)
                algninfo(i)%corr          = cxy(1)
                algninfo(i)%x             = cxy(2)
                algninfo(i)%y             = cxy(3)
                algninfo(i)%l_mirr        = .false.
            else
                algninfo_mirr(i-n)%e3     = pftcc%get_rot(irot)
                algninfo_mirr(i-n)%corr   = cxy(1)
                algninfo_mirr(i-n)%x      = cxy(2)
                algninfo_mirr(i-n)%y      = cxy(3)
                algninfo_mirr(i-n)%l_mirr = .true.
            endif
        end do
        !$omp end parallel do
        ! set mirror flags
        where( algninfo_mirr(:)%corr > algninfo(:)%corr ) algninfo = algninfo_mirr
        ! destruct
        call dealloc_imgarr(imgs_mirr)
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        call pftcc%kill
        call polartransform%kill
    end function match_imgs2ref

    function match_imgs2refs( nrefs, n, hp, lp, trs, refimgs, imgs ) result( algninfo )
        use simple_polarizer,         only: polarizer
        use simple_polarft_corrcalc,  only: polarft_corrcalc
        use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
        integer,                  intent(in)    :: nrefs, n
        real,                     intent(in)    :: hp, lp, trs
        class(image),             intent(inout) :: refimgs(nrefs), imgs(n) 
        integer,     parameter          :: MAXITS_SH = 60
        real,        allocatable        :: inpl_corrs(:)
        type(image), allocatable        :: imgs_mirr(:)
        type(pftcc_shsrch_grad)         :: grad_shsrch_obj(nthr_glob)
        type(polarizer)                 :: polartransform
        type(polarft_corrcalc)          :: pftcc
        type(inpl_struct)               :: algninfo_mirr(n,nrefs)
        type(inpl_struct), allocatable  :: algninfo(:,:)
        integer :: ldim(3), ldim_ref(3), box, kfromto(2), ithr, i, loc(1), nrots, irot, iref
        real    :: smpd, lims(2,2), lims_init(2,2), cxy(3)
        ldim       = imgs(1)%get_ldim()
        ldim_ref   = refimgs(1)%get_ldim()
        if( .not. all(ldim == ldim_ref) ) THROW_HARD('Incongruent logical image dimensions (imgs & imgs_ref)')
        box        = ldim(1)
        smpd       = imgs(1)%get_smpd()
        kfromto(1) = max(2, calc_fourier_index(hp, box, smpd))
        kfromto(2) =        calc_fourier_index(lp, box, smpd)
        ! create mirrored versions of the images
        call alloc_imgarr(n, ldim, smpd, imgs_mirr)
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, n
            call imgs_mirr(i)%copy(imgs(i))
            call imgs_mirr(i)%mirror('x')
        end do
        !$omp end parallel do
        ! initialize pftcc, polarizer
        call pftcc%new(nrefs, [1,2*n], kfromto) ! 2*n because of mirroring
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
        ! set the reference transforms
        do iref = 1, nrefs
            call polartransform%polarize(pftcc, refimgs(iref), iref, isptcl=.false., iseven=.true.)
        enddo
        ! set the particle transforms
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, 2 * n
            if( i <= n )then
                call imgs(i)%fft()
                call polartransform%polarize(pftcc, imgs(i),        i, isptcl=.true., iseven=.true.)
                call imgs(i)%ifft
            else
                call imgs_mirr(i-n)%fft()
                call polartransform%polarize(pftcc, imgs_mirr(i-n), i, isptcl=.true., iseven=.true.)
                call imgs_mirr(i-n)%ifft()
            endif
        end do
        !$omp end parallel do
        call pftcc%memoize_refs
        call pftcc%memoize_ptcls
        ! register imgs to references
        nrots = pftcc%get_nrots()
        allocate(inpl_corrs(nrots), algninfo(n,nrefs))
        !$omp parallel do default(shared) private(i,ithr,inpl_corrs,loc,irot,cxy) schedule(static) proc_bind(close)
        do i = 1, 2 * n
            ithr = omp_get_thread_num() + 1
            do iref = 1, nrefs
                call pftcc%gencorrs(iref, i, inpl_corrs)
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
                    algninfo(i,iref)%e3            = pftcc%get_rot(irot)
                    algninfo(i,iref)%corr          = cxy(1)
                    algninfo(i,iref)%x             = cxy(2)
                    algninfo(i,iref)%y             = cxy(3)
                    algninfo(i,iref)%l_mirr        = .false.
                else
                    algninfo_mirr(i-n,iref)%e3     = pftcc%get_rot(irot)
                    algninfo_mirr(i-n,iref)%corr   = cxy(1)
                    algninfo_mirr(i-n,iref)%x      = cxy(2)
                    algninfo_mirr(i-n,iref)%y      = cxy(3)
                    algninfo_mirr(i-n,iref)%l_mirr = .true.
                endif
            end do
        end do
        !$omp end parallel do
        ! set mirror flags
        do iref = 1, nrefs
            where( algninfo_mirr(:,iref)%corr > algninfo(:,iref)%corr ) algninfo(:,iref) = algninfo_mirr(:,iref)
        enddo
        ! destruct
        call dealloc_imgarr(imgs_mirr)
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        call pftcc%kill
        call polartransform%kill
    end function match_imgs2refs

    function rtsq_imgs( n, algninfo, imgs ) result( imgs_aligned )
        integer,            intent(in)    :: n
        type(inpl_struct),  intent(in)    :: algninfo(n)
        class(image),       intent(inout) :: imgs(n)
        type(image),        allocatable   :: imgs_aligned(:)
        real(kind=c_float), allocatable   :: rmat_rot(:,:,:)
        type(image),        allocatable   :: imgs_heap(:)
        integer :: ldim(3), i, ithr
        real    :: smpd
        ldim = imgs(1)%get_ldim()
        smpd = imgs(1)%get_smpd()
        call alloc_imgarr(n,         ldim, smpd, imgs_aligned)
        call alloc_imgarr(nthr_glob, ldim, smpd, imgs_heap)
        allocate(rmat_rot(ldim(1),ldim(2),1), source=0.)
        !$omp parallel do default(shared) private(i,ithr,rmat_rot) schedule(static) proc_bind(close)
        do i = 1, n
            if( algninfo(i)%l_mirr )then
                ithr = omp_get_thread_num() + 1
                call imgs_heap(ithr)%copy(imgs(i))
                call imgs_heap(ithr)%mirror('x')
                call imgs_heap(ithr)%fft
                call imgs_heap(ithr)%shift2Dserial([-algninfo(i)%x,-algninfo(i)%y])
                call imgs_heap(ithr)%ifft
                call imgs_heap(ithr)%rtsq_serial(algninfo(i)%e3, 0., 0., rmat_rot)
            else
                call imgs(i)%fft
                call imgs(i)%shift2Dserial([-algninfo(i)%x,-algninfo(i)%y])
                call imgs(i)%ifft
                call imgs(i)%rtsq_serial(algninfo(i)%e3, 0., 0., rmat_rot)
            endif
            call imgs_aligned(i)%set_rmat(rmat_rot, .false.)
        end do
        !$omp end parallel do
        call dealloc_imgarr(imgs_heap)
    end function rtsq_imgs

end module simple_strategy2D_utils
