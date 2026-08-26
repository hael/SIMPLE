!@descr: utility routines for 2D class average processing
module simple_strategy2D_utils
use simple_pftc_srch_api
use simple_clustering_utils,  only: cluster_dmat
use simple_cmdline,           only: cmdline
use simple_histogram,         only: histogram
use simple_image_bin,         only: image_bin
use simple_image_msk,         only: density_inoutside_mask
use simple_pftc_shsrch_grad,  only: pftc_shsrch_grad
use simple_pspecs,            only: pspecs
use simple_segmentation,      only: otsu_img
use simple_sp_project,        only: sp_project
use simple_builder,           only: builder
implicit none

public :: id_junk_and_prep_cavgs4clust, prep_cavgs4clust, id_junk, flag_non_junk_cavgs, calc_cluster_cavgs_dmat
public :: calc_match_cavgs_dmat, align_and_score_cavg_clusters, write_aligned_cavgs, calc_cavg_offset
public :: test_cavg_registration
public :: calc_cavg_pairwise_algninfo, calc_cavg_sigstats_components
private
#include "simple_local_flags.inc"

! objective function weights
real,    private, parameter :: W_HIST_SIG=0.5, W_POW_SIG=(1. - W_HIST_SIG)
real,    private, parameter :: W_SIG=0.5, W_CC=0.2, W_RES=0.3
real,    private, parameter :: W_SCORE_HOMO=0.2, W_SCORE_RES=0.5, W_SCORE_CLUST=0.3
! common control params
integer, private, parameter :: NHISTBINS = 128
real,    private, parameter :: HP_SPEC=20., LP_SPEC=6.

type, public :: cavg_sigstats_matrices
    integer :: histogram_bins = NHISTBINS
    real    :: power_hp       = HP_SPEC
    real    :: power_lp       = LP_SPEC
    real, allocatable :: power_distance(:,:)
    real, allocatable :: histogram_tvd(:,:)
    real, allocatable :: histogram_jsd(:,:)
    real, allocatable :: histogram_hellinger(:,:)
    real, allocatable :: histogram_distance(:,:)
    real, allocatable :: signal_stats_distance(:,:)
contains
    procedure :: kill => kill_cavg_sigstats_matrices
end type cavg_sigstats_matrices

contains

    subroutine id_junk_and_prep_cavgs4clust( spproj, cavg_imgs, mskdiam, clspops, clsinds, l_non_junk, mm )
        class(sp_project),        intent(inout) :: spproj
        type(image), allocatable, intent(inout) :: cavg_imgs(:)
        real,                     intent(in)    :: mskdiam
        integer,     allocatable, intent(inout) :: clspops(:), clsinds(:)
        logical,     allocatable, intent(inout) :: l_non_junk(:)
        real,        allocatable, intent(inout) :: mm(:,:)
        real,              parameter  :: LP_BIN = 20.
        logical,           parameter  :: DEBUG = .true.
        real,    allocatable :: frcs(:,:), filter(:)
        logical, allocatable :: l_msk(:,:,:)
        type(string)         :: frcs_fname
        type(image)          :: img_msk
        type(class_frcs)     :: clsfrcs
        integer              :: ncls, ldim(3), box, filtsz, ncls_sel, cnt, i, j
        real                 :: smpd, mskrad
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
            THROW_HARD('FRC file: '//frcs_fname%to_char()//' does not exist!')
        endif
        if( allocated(l_non_junk) ) deallocate(l_non_junk)
        call flag_non_junk_cavgs(cavg_imgs, LP_BIN, mskrad, l_non_junk, spproj%os_cls2D)
        if( DEBUG )then
            cnt = 0
            do i = 1, ncls
                if( .not. l_non_junk(i) )then
                    cnt = cnt + 1
                    call cavg_imgs(i)%write(string('cavgs_junk.mrc'), cnt)
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
        call img_msk%memoize_mask_coords
        call img_msk%mask2D_hard(mskrad)
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
            call cavg_imgs(i)%mask2D_soft(mskrad, backgr=0.)
            ! stash minmax
            mm(i,:) = cavg_imgs(i)%minmax(mskrad)
        end do
        !$omp end parallel do
        call clsfrcs%kill
    end subroutine id_junk_and_prep_cavgs4clust

    subroutine id_junk( spproj, mskdiam, l_non_junk )
        class(sp_project),        intent(inout) :: spproj
        real,                     intent(in)    :: mskdiam
        logical,     allocatable, intent(inout) :: l_non_junk(:)
        type(image), allocatable                :: cavg_imgs(:)
        real,        parameter                  :: LP_BIN = 20.
        integer                                 :: ldim(3), box
        real                                    :: smpd, mskrad
        cavg_imgs = read_cavgs_into_imgarr(spproj)
        smpd      = cavg_imgs(1)%get_smpd()
        ldim      = cavg_imgs(1)%get_ldim()
        box       = ldim(1)
        mskrad    = min(real(box/2) - COSMSKHALFWIDTH - 1., 0.5 * mskdiam/smpd)
        if( allocated(l_non_junk) ) deallocate(l_non_junk)
        call flag_non_junk_cavgs(cavg_imgs, LP_BIN, mskrad, l_non_junk)
        call dealloc_imgarr(cavg_imgs)
    end subroutine id_junk

    subroutine prep_cavgs4clust( spproj, cavg_imgs, mskdiam, clspops, clsinds, mask, mm )
        class(sp_project),        intent(inout) :: spproj
        type(image), allocatable, intent(inout) :: cavg_imgs(:)
        real,                     intent(in)    :: mskdiam
        integer,     allocatable, intent(inout) :: clspops(:), clsinds(:)
        logical,                  intent(in)    :: mask(:)
        real,        allocatable, intent(inout) :: mm(:,:)
        real,              parameter  :: LP_BIN = 20.
        real,    allocatable :: frcs(:,:), filter(:)
        logical, allocatable :: l_msk(:,:,:)
        type(string)         :: frcs_fname
        type(image)          :: img_msk
        type(class_frcs)     :: clsfrcs
        integer              :: ncls, ldim(3), box, filtsz, ncls_sel, i, j
        real                 :: smpd, mskrad
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
            THROW_HARD('FRC file: '//frcs_fname%to_char()//' does not exist!')
        endif
        ! re-create cavg_imgs
        ncls_sel = count(mask)
        write(logfhandle,'(A,I5)') '# classes left after mask application ', ncls_sel
        call dealloc_imgarr(cavg_imgs)
        cavg_imgs = read_cavgs_into_imgarr(spproj, mask=mask)
        ! keep track of the original class indices
        if( allocated(clsinds) ) deallocate(clsinds)
        clsinds = pack((/(i,i=1,ncls)/), mask=mask)
        ! update class populations
        clspops = pack(clspops, mask=mask)
        ! create the stuff needed in the loop
        allocate(frcs(filtsz,ncls_sel), filter(filtsz), mm(ncls_sel,2), source=0.)
        ! prep mask
        call img_msk%new([box,box,1], smpd)
        img_msk = 1.
        call img_msk%memoize_mask_coords
        call img_msk%mask2D_hard(mskrad)
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
            call cavg_imgs(i)%mask2D_soft(mskrad, backgr=0.)
            ! stash minmax
            mm(i,:) = cavg_imgs(i)%minmax(mskrad)
        end do
        !$omp end parallel do
        call clsfrcs%kill
    end subroutine prep_cavgs4clust

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
        real    :: cen(2),dynrange, smpd, large_msk, ave, sdev, maxv, minv
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
        large_msk = real(ldim(1)/2-1)
        if( allocated(l_non_junk) ) deallocate(l_non_junk)
        allocate(l_non_junk(ncls), source=.false.)
        kfromto(1) = calc_fourier_index(HP_SPEC, ldim(1), smpd)
        kfromto(2) = calc_fourier_index(LP_SPEC, ldim(1), smpd)
        call alloc_imgarr(nthr_glob, ldim, smpd, cavg_threads)
        call cavg_threads(1)%memoize_mask_coords
        !$omp parallel do default(shared) proc_bind(close) schedule(static)&
        !$omp private(icls,ithr,pspec,dynrange,l_dens_inoutside,nin,nout,nmsk,cen,ave,sdev,maxv,minv)
        do icls = 1, ncls
            ithr = omp_get_thread_num() + 1
            call cavg_threads(ithr)%copy(cavgs(icls))
            ! variance criterion
            call cavg_threads(ithr)%stats('foreground', ave, sdev, maxv, minv, large_msk)
            l_non_junk(icls) = sdev*sdev <=  ABS_VAR_THRESHOLD
            if( .not.l_non_junk(icls) ) cycle
            ! size & position criteria
            call cavg_threads(ithr)%div_below(0., 10.) ! reduce influence of negative values
            call cavg_threads(ithr)%norm
            call density_inoutside_mask(cavg_threads(ithr), lp_bin, msk, nin, nout, nmsk, cen)
            l_dens_inoutside = (nout > 0) .or.&                           ! object oustide of mask
                              &(real(nin)/real(nmsk) > RATIO_THRESHOLD)   ! or object too big inside mask
            call cavg_threads(ithr)%mask2D_soft(msk, backgr=0.)
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

    subroutine calc_sigstats_dmats( params, cavg_imgs, oa_minmax, dmat_sig )
        class(parameters),    intent(in)    :: params
        class(image),         intent(inout) :: cavg_imgs(:)
        real,                 intent(in)    :: oa_minmax(2)
        real, allocatable,    intent(inout) :: dmat_sig(:,:)
        type(cavg_sigstats_matrices) :: components
        call calc_cavg_sigstats_components(params, cavg_imgs, oa_minmax, components)
        if( allocated(dmat_sig) ) deallocate(dmat_sig)
        call move_alloc(components%signal_stats_distance, dmat_sig)
        call components%kill()
    end subroutine calc_sigstats_dmats

    subroutine calc_cavg_sigstats_components( params, cavg_imgs, oa_minmax, components )
        class(parameters),              intent(in)    :: params
        class(image),                   intent(inout) :: cavg_imgs(:)
        real,                           intent(in)    :: oa_minmax(2)
        type(cavg_sigstats_matrices),   intent(inout) :: components
        type(histogram),      allocatable   :: hists(:)
        real,                 allocatable   :: scaled_power(:,:)
        type(pspecs) :: pows
        integer      :: ncls_sel, i, j
        call components%kill()
        ncls_sel = size(cavg_imgs)
        ! generate histograms
        allocate(hists(ncls_sel))
        do i = 1, ncls_sel
            call hists(i)%new(cavg_imgs(i), NHISTBINS, minmax=oa_minmax, radius=params%msk)
        end do
        ! generate power spectra and associated distance/similarity matrix
        ! create pspecs object
        call pows%new(cavg_imgs, params%msk, HP_SPEC, LP_SPEC)
        ! create a joint similarity matrix for clustering based on spectral profile and in-plane invariant correlation
        components%power_distance = pows%calc_distmat()
        scaled_power = components%power_distance
        call normalize_minmax(scaled_power)
        ! calculate histogram-based distance matrices
        allocate(components%histogram_tvd(ncls_sel,ncls_sel), &
                 components%histogram_jsd(ncls_sel,ncls_sel), &
                 components%histogram_hellinger(ncls_sel,ncls_sel), source=0.)
        !$omp parallel do default(shared) private(i,j)&
        !$omp schedule(dynamic) proc_bind(close)
        do i = 1, ncls_sel - 1
            do j = i + 1, ncls_sel
                components%histogram_tvd(i,j) = hists(i)%TVD(hists(j))
                components%histogram_tvd(j,i) = components%histogram_tvd(i,j)
                components%histogram_jsd(i,j) = hists(i)%JSD(hists(j))
                components%histogram_jsd(j,i) = components%histogram_jsd(i,j)
                components%histogram_hellinger(i,j) = hists(i)%HD(hists(j))
                components%histogram_hellinger(j,i) = components%histogram_hellinger(i,j)
            end do
        end do
        !$omp end parallel do
        call hists(:)%kill
        components%histogram_distance = merge_dmats(components%histogram_tvd, &
            components%histogram_jsd, components%histogram_hellinger)
        components%signal_stats_distance = W_HIST_SIG * components%histogram_distance + &
                                           W_POW_SIG  * scaled_power
        ! destruct
        call pows%kill
        if( allocated(hists) ) deallocate(hists)
        if( allocated(scaled_power) ) deallocate(scaled_power)
    end subroutine calc_cavg_sigstats_components

    subroutine kill_cavg_sigstats_matrices( self )
        class(cavg_sigstats_matrices), intent(inout) :: self
        if( allocated(self%power_distance)       ) deallocate(self%power_distance)
        if( allocated(self%histogram_tvd)        ) deallocate(self%histogram_tvd)
        if( allocated(self%histogram_jsd)        ) deallocate(self%histogram_jsd)
        if( allocated(self%histogram_hellinger)  ) deallocate(self%histogram_hellinger)
        if( allocated(self%histogram_distance)   ) deallocate(self%histogram_distance)
        if( allocated(self%signal_stats_distance)) deallocate(self%signal_stats_distance)
    end subroutine kill_cavg_sigstats_matrices

    subroutine calc_sigstats_dmats_ref( params, cavg_imgs_ref, cavg_imgs_match, oa_minmax, dmat_sig )
        use simple_clustering_utils, only: cluster_dmat, labels2smat
        class(parameters),    intent(in)    :: params
        class(image),         intent(inout) :: cavg_imgs_ref(:), cavg_imgs_match(:)
        real,                 intent(in)    :: oa_minmax(2)
        real, allocatable,    intent(inout) :: dmat_sig(:,:)
        type(histogram),      allocatable   :: hists_ref(:), hists_match(:)
        real,                 allocatable   :: dmat_pow(:,:), dmat_tvd(:,:), dmat_jsd(:,:)
        real,                 allocatable   :: dmat_hd(:,:), dmat_hist(:,:)
        type(pspecs) :: pows_ref, pows_match
        integer      :: ncls_ref, ncls_match, i, j
        ncls_ref   = size(cavg_imgs_ref)
        ncls_match = size(cavg_imgs_match)
        if( allocated(dmat_sig) ) deallocate(dmat_sig)
        ! generate histograms
        allocate(hists_ref(ncls_ref), hists_match(ncls_match))
        do i = 1, ncls_ref
            call hists_ref(i)%new(cavg_imgs_ref(i),     NHISTBINS, minmax=oa_minmax, radius=params%msk)
        end do
        do i = 1, ncls_match
            call hists_match(i)%new(cavg_imgs_match(i), NHISTBINS, minmax=oa_minmax, radius=params%msk)
        end do
        ! generate power spectra and associated distance/similarity matrix
        call pows_ref%new(cavg_imgs_ref,     params%msk, HP_SPEC, LP_SPEC)
        call pows_match%new(cavg_imgs_match, params%msk, HP_SPEC, LP_SPEC)
        dmat_pow = pows_ref%calc_distmat(pows_match)
        call normalize_minmax(dmat_pow)
        ! calculate histogram-based distance matrices
        allocate(dmat_tvd(ncls_ref,ncls_match), dmat_jsd(ncls_ref,ncls_match), dmat_hd(ncls_ref,ncls_match), source=0.)
        !$omp parallel do default(shared) private(i,j)&
        !$omp schedule(static) proc_bind(close)
        do i = 1, ncls_ref
            do j = 1, ncls_match 
                dmat_tvd(i,j) = hists_ref(i)%TVD(hists_match(j))
                dmat_jsd(i,j) = hists_ref(i)%JSD(hists_match(j))
                dmat_hd(i,j)  = hists_ref(i)%HD(hists_match(j))
            end do
        end do
        !$omp end parallel do
        dmat_hist = merge_dmats(dmat_tvd, dmat_jsd, dmat_hd)     ! the different histogram distances are given equal weight
        dmat_sig = W_HIST_SIG * dmat_hist + W_POW_SIG * dmat_pow ! this is the joint signal distance matrix
        ! destruct
        call hists_ref(:)%kill
        call hists_match(:)%kill
        call pows_ref%kill
        call pows_match%kill
        if( allocated(dmat_pow)  ) deallocate(dmat_pow)
        if( allocated(dmat_hist) ) deallocate(dmat_hist)
    end subroutine calc_sigstats_dmats_ref

    subroutine calc_cc_and_res_dmats( params, cavg_imgs, hp, lp, trs, dmat_cc, dmat_res )
        class(parameters), intent(in)    :: params
        class(image),      intent(inout) :: cavg_imgs(:)
        real,              intent(in)    :: hp, lp, trs
        real,allocatable,  intent(inout) :: dmat_cc(:,:), dmat_res(:,:)
        type(inpl_struct), allocatable   :: algninfo(:,:)
        real,              allocatable   :: ccmat(:,:)
        integer :: ncavgs, i, j
        ncavgs   = size(cavg_imgs)
        algninfo = calc_cavg_pairwise_algninfo(params, cavg_imgs, hp, lp, trs)
        if( allocated(dmat_res) ) deallocate(dmat_res)
        if( allocated(dmat_cc)  ) deallocate(dmat_cc)
        allocate(dmat_res(ncavgs,ncavgs), ccmat(ncavgs,ncavgs), source=0.)
        !$omp parallel do default(shared) private(i,j) schedule(dynamic) proc_bind(close)
        do i = 1, ncavgs - 1
            do j = i + 1, ncavgs
                dmat_res(i,j) = 1./real(algninfo(i,j)%find_fsc05)
                dmat_res(j,i) = dmat_res(i,j)
                ccmat(i,j)    = algninfo(i,j)%corr
                ccmat(j,i)    = ccmat(i,j)
            enddo 
        enddo
        !$omp end parallel do
        forall( i = 1:ncavgs ) dmat_res(i,i) = 0.
        call normalize_minmax(dmat_res)
        dmat_cc = smat2dmat(ccmat)
        deallocate(algninfo, ccmat)
    end subroutine calc_cc_and_res_dmats

    function calc_cavg_pairwise_algninfo( params, cavg_imgs, hp, lp, trs ) result( algninfo )
        class(parameters), intent(in)    :: params
        class(image),      intent(inout) :: cavg_imgs(:)
        real,              intent(in)    :: hp, lp, trs
        type(inpl_struct), allocatable   :: algninfo(:,:)
        if( trim(params%objfun) /= 'cc' ) &
            THROW_HARD('calc_cavg_pairwise_algninfo requires objfun=cc')
        if( trim(params%ctf) /= 'no' ) &
            THROW_HARD('calc_cavg_pairwise_algninfo requires ctf=no')
        algninfo = match_imgs(params, hp, lp, trs, cavg_imgs, cavg_imgs)
    end function calc_cavg_pairwise_algninfo

    subroutine calc_cc_and_res_dmats_ref( params, cavg_imgs_ref, cavg_imgs_match, hp, lp, trs, dmat_cc, dmat_res )
        class(parameters), intent(in)    :: params
        class(image),      intent(inout) :: cavg_imgs_ref(:), cavg_imgs_match(:)
        real,              intent(in)    :: hp, lp, trs
        real,allocatable,  intent(inout) :: dmat_cc(:,:), dmat_res(:,:)
        type(inpl_struct), allocatable   :: algninfo(:,:)
        real,              allocatable   :: ccmat(:,:)
        integer :: ncls_ref, ncls_match, i, j
        ncls_ref   = size(cavg_imgs_ref)
        ncls_match = size(cavg_imgs_match)
        algninfo   = match_imgs(params, hp, lp, trs, cavg_imgs_ref, cavg_imgs_match)
        if( allocated(dmat_res) ) deallocate(dmat_res)
        if( allocated(dmat_cc)  ) deallocate(dmat_cc)
        allocate(dmat_res(ncls_ref,ncls_match), ccmat(ncls_ref,ncls_match), source=0.)
        !$omp parallel do default(shared) private(i,j) schedule(static) collapse(2) proc_bind(close)
        do i = 1, ncls_ref
            do j = 1, ncls_match
                dmat_res(i,j) = 1./real(algninfo(i,j)%find_fsc05)
                ccmat(i,j)    = algninfo(i,j)%corr
            enddo 
        enddo
        !$omp end parallel do
        call normalize_minmax(dmat_res)
        dmat_cc = smat2dmat(ccmat)
        deallocate(algninfo, ccmat)
    end subroutine calc_cc_and_res_dmats_ref

    function calc_cluster_cavgs_dmat( params, cavg_imgs, oa_minmax, which ) result( dmat )
        class(parameters),    intent(in)    :: params
        class(image),         intent(inout) :: cavg_imgs(:)
        real,                 intent(in)    :: oa_minmax(2)
        character(len=*),     intent(in)    :: which
        real,                 allocatable   :: dmat_sig(:,:), dmat_cc(:,:), dmat_res(:,:), dmat(:,:)
        write(logfhandle,'(A)') '>>> GENERATING DISTANCE MATRICES FOR SIGNAL STATISTICS'
        call calc_sigstats_dmats(params, cavg_imgs, oa_minmax, dmat_sig)
        write(logfhandle,'(A)') '>>> PAIRWISE CORRELATIONS & FRCS THROUGH FULL IN-PLANE SEARCH'
        call calc_cc_and_res_dmats(params, cavg_imgs, params%hp, params%lp, params%trs, dmat_cc, dmat_res)
        select case(trim(which))    
            case('sig')
                dmat = dmat_sig
            case('cc')
                dmat = dmat_cc
            case('res')
                dmat = dmat_res
            case DEFAULT
                dmat = W_SIG * dmat_sig + W_CC * dmat_cc + W_RES * dmat_res
        end select
        call normalize_minmax(dmat)
    end function calc_cluster_cavgs_dmat

    function calc_match_cavgs_dmat( params, cavg_imgs_ref, cavg_imgs_match, oa_minmax, which ) result( dmat )
        class(parameters),    intent(in)    :: params
        class(image),         intent(inout) :: cavg_imgs_ref(:), cavg_imgs_match(:)
        real,                 intent(in)    :: oa_minmax(2)
        character(len=*),     intent(in)    :: which
        real,                 allocatable   :: dmat_sig(:,:), dmat_cc(:,:), dmat_res(:,:), dmat(:,:)
        write(logfhandle,'(A)') '>>> GENERATING DISTANCE MATRICES FOR SIGNAL STATISTICS'
        call calc_sigstats_dmats_ref(params, cavg_imgs_ref, cavg_imgs_match, oa_minmax, dmat_sig)
        write(logfhandle,'(A)') '>>> PAIRWISE CORRELATIONS & FRCS THROUGH FULL IN-PLANE SEARCH'
        call calc_cc_and_res_dmats_ref(params, cavg_imgs_ref, cavg_imgs_match, params%hp, params%lp, params%trs, dmat_cc, dmat_res)
        select case(trim(which))    
            case('sig')
                dmat = dmat_sig
            case('cc')
                dmat = dmat_cc
            case('res')
                dmat = dmat_res
            case DEFAULT
                dmat = W_SIG * dmat_sig + W_CC * dmat_cc + W_RES * dmat_res
        end select
        call normalize_minmax(dmat)
    end function calc_match_cavgs_dmat

    function align_and_score_cavg_clusters( params, dmat, cavg_imgs, clspops, i_medoids, labels) result( clust_info_arr ) ! , clustscores )
        class(parameters),    intent(in)    :: params
        real,                 intent(in)    :: dmat(:,:)
        class(image),         intent(inout) :: cavg_imgs(:)
        integer,              intent(in)    :: clspops(:)
        integer,              intent(inout) :: i_medoids(:), labels(:)
        real,                 parameter     :: SCORE_THRES=60.
        type(clust_info),     allocatable   :: clust_info_arr(:)
        logical,              allocatable   :: l_msk(:,:,:)
        type(image)  :: img_msk 
        integer      :: ncls_sel, nclust, iclust
        ncls_sel = size(cavg_imgs)
        ! prep mask
        call img_msk%new([params%box,params%box,1], params%smpd)
        img_msk = 1.
        call img_msk%memoize_mask_coords
        call img_msk%mask2D_hard(params%msk)
        l_msk = img_msk%bin2logical()
        call img_msk%kill
        ! align clusters to medoids and gather information
        clust_info_arr = align_clusters2medoids(params, i_medoids, labels, cavg_imgs, params%hp, params%lp, params%trs, l_msk )
        nclust         = maxval(labels)
        ! set particle populations
        do iclust = 1, nclust
            clust_info_arr(iclust)%nptcls = sum(clspops, mask=labels == iclust)
        end do
        ! calculate scores
        call calc_scores
        ! rank clusters
        call rank_clusters
        ! good/bad cluster identification through tresholding
        call identify_good_bad_clusters

    contains

        subroutine calc_scores
            integer :: iclust, icls, cnt
            real    :: euclid_max, res_max, clustscore_min
            ! HOMOGENEITY SCORE
            clust_info_arr(:)%homogeneity = clust_info_arr(:)%euclid
            euclid_max = maxval(clust_info_arr(:)%euclid)
            where( clust_info_arr(:)%pop < 2 ) clust_info_arr(:)%homogeneity = euclid_max
            call dists2scores_percen(clust_info_arr(:)%homogeneity)
            ! RESOLUTION SCORE
            clust_info_arr(:)%resscore = clust_info_arr(:)%res
            res_max = maxval(clust_info_arr(:)%res)
            where( clust_info_arr(:)%pop < 2 ) clust_info_arr(:)%resscore = res_max
            call dists2scores_percen(clust_info_arr(:)%resscore)
            ! CLUSTSCORE
            do iclust = 1, nclust
                clust_info_arr(iclust)%clustscore  = 0.
                cnt = 0
                do icls = 1, ncls_sel 
                    if( labels(icls) == iclust )then
                        clust_info_arr(iclust)%clustscore  = clust_info_arr(iclust)%clustscore  +      dmat(icls,i_medoids(iclust))
                        cnt = cnt + 1
                    endif
                enddo
                clust_info_arr(iclust)%clustscore  = clust_info_arr(iclust)%clustscore  / real(cnt)
            end do
            clustscore_min  = minval(clust_info_arr(:)%clustscore)
            where( clust_info_arr(:)%pop < 2 )
                clust_info_arr(:)%clustscore  = clustscore_min
            endwhere
            call dists2scores_percen(clust_info_arr(:)%clustscore)
            ! JOINT SCORE
            clust_info_arr(:)%jointscore = W_SCORE_HOMO * clust_info_arr(:)%homogeneity + W_SCORE_RES * clust_info_arr(:)%resscore + W_SCORE_CLUST * clust_info_arr(:)%clustscore
            call scores2scores_percen(clust_info_arr(:)%jointscore)
        end subroutine calc_scores

        subroutine rank_clusters
            integer              :: i_medoids_ranked(nclust), rank_assign(ncls_sel)
            type(clust_info)     :: clust_info_arr_ranked(nclust)
            integer              :: rank, icls
            integer, allocatable :: clust_order(:)
            clust_order = scores2order(clust_info_arr(:)%jointscore)
            do rank = 1, nclust
                i_medoids_ranked(rank)      = i_medoids(clust_order(rank))
                clust_info_arr_ranked(rank) = clust_info_arr(clust_order(rank))
                do icls = 1, ncls_sel
                    if( labels(icls) == clust_order(rank) ) rank_assign(icls) = rank
                end do
            end do
            i_medoids      = i_medoids_ranked
            clust_info_arr = clust_info_arr_ranked
            labels         = rank_assign
            deallocate(clust_order)
        end subroutine rank_clusters

        subroutine identify_good_bad_clusters
            integer :: ngood
            clust_info_arr(:)%good_bad = 0
            ngood = count(clust_info_arr(:)%jointscore > SCORE_THRES)
            if( ngood < 2 )then
                clust_info_arr(:2)%good_bad     = 1
            else
                clust_info_arr(:ngood)%good_bad = 1
            endif
        end subroutine identify_good_bad_clusters

    end function align_and_score_cavg_clusters

    function align_clusters2medoids( params, i_medoids, labels, cavg_imgs, hp, lp, trs, l_msk ) result( clust_info_arr )
        class(parameters),    intent(in)    :: params
        integer,          intent(in)    :: i_medoids(:), labels(:)
        class(image),     intent(inout) :: cavg_imgs(:)
        real,             intent(in)    :: hp, lp, trs
        logical,          intent(in)    :: l_msk(:,:,:)
        type(clust_info), allocatable   :: clust_info_arr(:)
        real,             allocatable   :: frc(:)
        type(image),      allocatable   :: cluster_imgs(:)
        real,             allocatable   :: resvals(:), resarr(:)
        real,             parameter     :: FRAC_BEST_CAVGS=0.3
        integer :: cnt, i, filtsz, ldim(3), iclust, nclust
        real    :: smpd, rfoo, best_res, worst_res
        write(logfhandle,'(A)') '>>> ALIGNING THE CLUSTERS OF CLASS AVERAGES TO THEIR MEDOIDS'
        filtsz = cavg_imgs(1)%get_filtsz()
        smpd   = cavg_imgs(1)%get_smpd()
        ldim   = cavg_imgs(1)%get_ldim()
        resarr = get_resarr(ldim(1), smpd)
        nclust = maxval(labels)
        allocate(frc(filtsz), clust_info_arr(nclust))
        do iclust = 1, nclust
            clust_info_arr(iclust)%pop = count(labels == iclust)
            if( clust_info_arr(iclust)%pop == 0 ) cycle
            cluster_imgs = pack_imgarr(cavg_imgs, mask=labels == iclust)
            clust_info_arr(iclust)%algninfo%params = match_imgs2ref(params, hp, lp, trs, cavg_imgs(i_medoids(iclust)), cluster_imgs)
            call rtsq_imgs(clust_info_arr(iclust)%pop, clust_info_arr(iclust)%algninfo%params, cluster_imgs)
            ! estimate resolution
            ! FWD FT
            !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static) 
            do i = 1, size(cavg_imgs)
                if( i <= clust_info_arr(iclust)%pop ) call cluster_imgs(i)%fft
                call cavg_imgs(i)%fft
            end do
            !$omp end parallel do
            cnt = 0
            allocate(resvals(clust_info_arr(iclust)%pop), source=0.)
            do i = 1, clust_info_arr(iclust)%pop
                call cavg_imgs(i_medoids(iclust))%fsc(cluster_imgs(i), frc)
                if( .not. all(frc > 0.5) )then ! excluding the medoid
                    cnt = cnt + 1
                    call get_resolution(frc, resarr, rfoo, resvals(cnt))
                endif
            end do
            ! BWD FT
            !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static) 
            do i = 1, size(cavg_imgs)
                if( i <= clust_info_arr(iclust)%pop ) call cluster_imgs(i)%ifft
                call cavg_imgs(i)%ifft
            end do
            !$omp end parallel do
            ! calculate Euclidean distances
            clust_info_arr(iclust)%euclid = 0.
            do i = 1, clust_info_arr(iclust)%pop
                clust_info_arr(iclust)%euclid = clust_info_arr(iclust)%euclid + cavg_imgs(i_medoids(iclust))%sqeuclid(cluster_imgs(i), l_msk)
            end do
            clust_info_arr(iclust)%euclid = clust_info_arr(iclust)%euclid / real(clust_info_arr(iclust)%pop)
            ! report resolution as the average of the best agreeing 25% within a cluster
            clust_info_arr(iclust)%res = avg_frac_smallest(resvals(:cnt), FRAC_BEST_CAVGS)
            ! destruct
            call dealloc_imgarr(cluster_imgs)
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
        type(image), allocatable :: cluster_imgs(:)
        integer :: iclust, nclust
        write(logfhandle,'(A)') '>>> ROTATING & SHIFTING CLASS AVERAGES'
        nclust = size(clust_info_arr)
        do iclust = 1, nclust
            if( clust_info_arr(iclust)%pop == 0 ) cycle
            cluster_imgs = pack_imgarr(cavg_imgs, mask=labels == iclust)
            call rtsq_imgs(clust_info_arr(iclust)%pop, clust_info_arr(iclust)%algninfo%params, cluster_imgs)
            call write_imgarr(cluster_imgs, string(trim(fbody)//int2str_pad(iclust,2)//trim(ext)))
            call dealloc_imgarr(cluster_imgs)
        end do
    end subroutine write_aligned_cavgs

    function match_imgs2ref( params, hp, lp, trs, img_ref, imgs ) result( algninfo )
        class(parameters), intent(in)    :: params
        real,              intent(in)    :: hp, lp, trs
        class(image),      intent(inout) :: imgs(:), img_ref
        integer,     parameter           :: MAXITS_SH = 60
        real,        allocatable         :: inpl_corrs(:)
        type(image), allocatable         :: imgs_mirr(:)
        type(pftc_shsrch_grad)           :: grad_shsrch_obj(nthr_glob)
        type(builder)                    :: build
        type(inpl_struct), allocatable   :: algninfo(:), algninfo_mirr(:)
        integer :: ldim(3), ldim_ref(3), box, kfromto(2), ithr, i, loc(1), nrots, irot, irot_seed, n, pdim_srch(3)
        real    :: smpd, lims(2,2), lims_init(2,2), joint_lims(3,2), joint_lims_loc(3,2), cxy(3), inpl_angle
        real(dp):: rotind_frac
        logical :: didft
        n        = size(imgs)
        ldim     = imgs(1)%get_ldim()
        ldim_ref = img_ref%get_ldim()
        if( .not. all(ldim == ldim_ref) )then
            print *, 'ldim     ', ldim
            print *, 'ldim_ref ', ldim_ref
            THROW_HARD('Incongruent logical image dimensions (imgs & img_ref)')
        endif
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
        ! initialize pftc, polarizer
        call build%pftc%new(params, 1, [1,2*n], kfromto) ! 2*n because of mirroring
        pdim_srch = build%pftc%get_pdim_srch()
        call img_ref%memoize4polarize(pdim_srch)
        ! in-plane search object objects for parallel execution
        lims(:,1)      = -trs
        lims(:,2)      =  trs
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        nrots = build%pftc%get_nrots()
        joint_lims(1:2,:) = lims
        joint_lims(3,:) = [1.-2., real(nrots)+2.]
        do ithr = 1, nthr_glob
            if( trim(params%inpl_cont) == 'yes' )then
                call grad_shsrch_obj(ithr)%new_joint(build, joint_lims, MAXITS_SH)
            else
                call grad_shsrch_obj(ithr)%new_legacy(build, lims, lims_init=lims_init, maxits=MAXITS_SH)
            endif
        end do
        ! set the reference transform
        didft = .false.
        if( .not. img_ref%is_ft() )then
            call img_ref%fft
            didft = .true.
        endif
        call build%pftc%polarize_ref_pft(img_ref, 1, iseven=.true., pdim=pdim_srch, oversamp=.false.)
        if( didft ) call img_ref%ifft
        ! set the particle transforms
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, 2 * n
            if( i <= n )then
                call imgs(i)%fft()
                call build%pftc%polarize_ptcl_pft(imgs(i), i, pdim=pdim_srch, oversamp=.false.)
                call imgs(i)%ifft
            else
                call imgs_mirr(i-n)%fft()
                call build%pftc%polarize_ptcl_pft(imgs_mirr(i-n), i, pdim=pdim_srch, oversamp=.false.)
                call imgs_mirr(i-n)%ifft()
            endif
        end do
        !$omp end parallel do
        call build%pftc%memoize_refs
        call build%pftc%memoize_ptcls
        ! register imgs to img_ref
        allocate(inpl_corrs(nrots), algninfo(n), algninfo_mirr(n))
        !$omp parallel do default(shared) private(i,ithr,inpl_corrs,loc,irot,irot_seed,cxy,rotind_frac,inpl_angle,joint_lims_loc)&
        !$omp schedule(static) proc_bind(close)
        do i = 1, 2 * n
            ithr = omp_get_thread_num() + 1
            call build%pftc%gen_objfun_vals(1, i, [0.,0.], inpl_corrs)
            loc  = maxloc(inpl_corrs)
            irot = loc(1)
            inpl_angle = build%pftc%get_rot(irot)
            call grad_shsrch_obj(ithr)%set_indices(1, i)
            if( trim(params%inpl_cont) == 'yes' )then
                joint_lims_loc = joint_lims
                joint_lims_loc(3,:) = [real(irot)-2., real(irot)+2.]
                call grad_shsrch_obj(ithr)%set_limits(joint_lims_loc)
                irot_seed = irot
                cxy = grad_shsrch_obj(ithr)%minimize_joint(irot, [0.,0.], &
                    &sh_rot=.true., rotind_frac=rotind_frac, irot_in=irot_seed)
                if( irot > 0 ) inpl_angle = real(modulo((rotind_frac-1.d0) * &
                    &real(build%pftc%get_dang(),dp), 360.d0))
            else
                cxy = grad_shsrch_obj(ithr)%minimize(irot=irot)
            endif
            if( irot == 0 )then ! no improved solution found, put back the old one
                cxy(1) = inpl_corrs(loc(1))
                cxy(2) = 0.
                cxy(3) = 0.
                irot   = loc(1)
            endif
            if( i <= n )then
                algninfo(i)%e3            = inpl_angle
                algninfo(i)%corr          = cxy(1)
                algninfo(i)%x             = cxy(2)
                algninfo(i)%y             = cxy(3)
                algninfo(i)%l_mirr        = .false.
            else
                algninfo_mirr(i-n)%e3     = inpl_angle
                algninfo_mirr(i-n)%corr   = cxy(1)
                algninfo_mirr(i-n)%x      = cxy(2)
                algninfo_mirr(i-n)%y      = cxy(3)
                algninfo_mirr(i-n)%l_mirr = .true.
            endif
        end do
        !$omp end parallel do
        ! select for mirroring
        where( algninfo_mirr(:)%corr > algninfo(:)%corr ) algninfo = algninfo_mirr
        ! destruct
        call dealloc_imgarr(imgs_mirr)
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        call build%pftc%kill
    end function match_imgs2ref

    function match_imgs( params, hp, lp, trs, imgs_ref, imgs_targ ) result( algninfo )
        use simple_parameters, only: parameters
        class(parameters), intent(in)    :: params
        real,              intent(in)    :: hp, lp, trs
        class(image),       intent(inout) :: imgs_ref(:), imgs_targ(:)
        integer,      parameter         :: MAXITS_SH = 60
        real,         parameter         :: FRC_CRIT = 0.5
        real,         allocatable       :: inpl_corrs(:), frc(:)
        type(image),  allocatable       :: imgs_targ_mirr(:)
        type(pftc_shsrch_grad)          :: grad_shsrch_obj(nthr_glob)
        type(builder)                   :: build
        type(inpl_struct), allocatable  :: algninfo(:,:), algninfo_mirr(:,:)
        integer :: ldim(3), box, kfromto(2), ithr, i, j, k, m, loc, nrots, irot, irot_seed, nrefs, ntargets, pdim_srch(3)
        real    :: smpd, lims(2,2), lims_init(2,2), joint_lims(3,2), joint_lims_loc(3,2)
        real    :: cxy(3), rotmat(2,2), inpl_angle
        real(dp):: rotind_frac
        logical :: l_rot_shvec
        nrefs      = size(imgs_ref)
        ntargets   = size(imgs_targ)
        ldim       = imgs_ref(1)%get_ldim()
        box        = ldim(1)
        smpd       = imgs_ref(1)%get_smpd()
        kfromto(1) = max(2, calc_fourier_index(hp, box, smpd))
        kfromto(2) =        calc_fourier_index(lp, box, smpd)
        ! initialize mirrores images, pftc, polarizer
        call alloc_imgarr(ntargets, ldim, smpd, imgs_targ_mirr)
        call build%pftc%new(params, nrefs, [1,2*ntargets], kfromto) ! 2*ntargets because of mirroring
        pdim_srch = build%pftc%get_pdim_srch()
        call imgs_ref(1)%memoize4polarize(pdim_srch)
        ! in-plane search object objects for parallel execution
        lims(:,1)      = -trs
        lims(:,2)      =  trs
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        nrots = build%pftc%get_nrots()
        joint_lims(1:2,:) = lims
        joint_lims(3,:) = [1.-2., real(nrots)+2.]
        do ithr = 1, nthr_glob
            if( trim(params%inpl_cont) == 'yes' )then
                call grad_shsrch_obj(ithr)%new_joint(build, joint_lims, MAXITS_SH)
            else
                call grad_shsrch_obj(ithr)%new_legacy(build, lims, lims_init=lims_init, &
                    &maxits=MAXITS_SH)
            endif
        end do
        !$omp parallel default(shared)  private(i,m)  proc_bind(close)
        ! set the reference transforms
        !$omp do schedule(static)
        do i = 1, nrefs
            call imgs_ref(i)%fft()
            call build%pftc%polarize_ref_pft(imgs_ref(i), i, iseven=.true., pdim=pdim_srch, oversamp=.false.)
            call imgs_ref(i)%ifft
        end do
        !$omp end do nowait
        ! set the particle transforms
        !$omp do schedule(static)
        do i = 1,ntargets
            ! target
            call imgs_targ(i)%fft()
            call build%pftc%polarize_ptcl_pft(imgs_targ(i), i, pdim=pdim_srch, oversamp=.false.)
            call imgs_targ(i)%ifft
            ! target mirror
            m = ntargets+i
            call imgs_targ_mirr(i)%copy(imgs_targ(i))
            call imgs_targ_mirr(i)%mirror('x')
            call imgs_targ_mirr(i)%fft()
            call build%pftc%polarize_ptcl_pft(imgs_targ_mirr(i), m, pdim=pdim_srch, oversamp=.false.)
        end do
        !$omp end do
        !$omp end parallel
        call dealloc_imgarr(imgs_targ_mirr) ! no longer needed
        call build%pftc%memoize_refs
        call build%pftc%memoize_ptcls
        ! register imgs
        allocate(inpl_corrs(nrots), algninfo(nrefs,ntargets),&
        &algninfo_mirr(nrefs,ntargets), frc(kfromto(1):kfromto(2)))
        !$omp parallel do private(i,j,k,m,ithr,inpl_corrs,loc,irot,cxy,frc,l_rot_shvec,rotind_frac,inpl_angle,joint_lims_loc)&
        !$omp default(shared) collapse(2) schedule(static) proc_bind(close)
        do i = 1, 2 * ntargets ! ptcls
            do j = 1, nrefs    ! refs
                ithr = omp_get_thread_num() + 1
                call build%pftc%gen_objfun_vals(j, i, [0.,0.], inpl_corrs)
                loc  = maxloc(inpl_corrs,dim=1)
                irot = loc
                inpl_angle = build%pftc%get_rot(irot)
                call grad_shsrch_obj(ithr)%set_indices(j, i)
                if( trim(params%inpl_cont) == 'yes' )then
                    joint_lims_loc = joint_lims
                    joint_lims_loc(3,:) = [real(irot)-2., real(irot)+2.]
                    call grad_shsrch_obj(ithr)%set_limits(joint_lims_loc)
                    irot_seed = irot
                    cxy = grad_shsrch_obj(ithr)%minimize_joint(irot, [0.,0.], &
                        &sh_rot=.false., rotind_frac=rotind_frac, irot_in=irot_seed)
                    if( irot > 0 ) inpl_angle = real(modulo((rotind_frac-1.d0) * &
                        &real(build%pftc%get_dang(),dp), 360.d0))
                else
                    cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.false.) ! FRC needs native shift
                endif
                l_rot_shvec = .true.
                if( irot == 0 )then ! no improved solution found, put back the old one
                    cxy(1)      = inpl_corrs(loc)
                    cxy(2)      = 0.
                    cxy(3)      = 0.
                    irot        = loc
                    l_rot_shvec = .false.
                endif
                call build%pftc%calc_frc(j, i, irot, cxy(2:3), frc)
                if( i <= ntargets )then
                    algninfo(j,i)%e3     = inpl_angle
                    algninfo(j,i)%corr   = cxy(1)
                    algninfo(j,i)%l_mirr = .false.
                    algninfo(j,i)%find_fsc05 = 1
                    do k = kfromto(1),kfromto(2)
                        if( frc(k) >= FRC_CRIT )then
                            algninfo(j,i)%find_fsc05 = k
                        else
                            exit
                        endif
                    end do
                    ! rotate the shift vector to the frame of reference
                    if( l_rot_shvec )then
                        call rotmat2d(inpl_angle, rotmat)
                        cxy(2:) = matmul(cxy(2:), rotmat)
                    endif
                    algninfo(j,i)%x = cxy(2)
                    algninfo(j,i)%y = cxy(3)
                else
                    ! mirror
                    m = i - ntargets
                    algninfo_mirr(j,m)%e3     = inpl_angle
                    algninfo_mirr(j,m)%corr   = cxy(1)
                    algninfo_mirr(j,m)%l_mirr = .true.
                    algninfo_mirr(j,m)%find_fsc05 = 1
                    do k = kfromto(1),kfromto(2)
                        if( frc(k) >= FRC_CRIT )then
                            algninfo_mirr(j,m)%find_fsc05 = k
                        else
                            exit
                        endif
                    end do
                    ! rotate the shift vector to the frame of reference
                    if( l_rot_shvec )then
                        call rotmat2d(inpl_angle, rotmat)
                        cxy(2:) = matmul(cxy(2:), rotmat)
                    endif
                    algninfo_mirr(j,m)%x = cxy(2)
                    algninfo_mirr(j,m)%y = cxy(3)
                endif
            end do
        end do
        !$omp end parallel do
        ! select for mirroring
        where( algninfo_mirr(:,:)%corr > algninfo(:,:)%corr ) algninfo = algninfo_mirr
        ! destruct
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        call build%pftc%kill
        deallocate(algninfo_mirr)
    end function match_imgs

    subroutine rtsq_imgs( n, algninfo, imgs )
        integer,            intent(in)    :: n
        type(inpl_struct),  intent(in)    :: algninfo(n)
        class(image),       intent(inout) :: imgs(n)
        real(kind=c_float), allocatable   :: rmat_rot(:,:,:)
        integer :: ldim(3), i
        real    :: smpd
        ldim = imgs(1)%get_ldim()
        smpd = imgs(1)%get_smpd()
        allocate(rmat_rot(ldim(1),ldim(2),1), source=0.)
        !$omp parallel do default(shared) private(i,rmat_rot) schedule(static) proc_bind(close)
        do i = 1, n
            if( algninfo(i)%l_mirr ) call imgs(i)%mirror('x')                
            call imgs(i)%fft
            call imgs(i)%shift2Dserial([-algninfo(i)%x,-algninfo(i)%y])
            call imgs(i)%ifft
            call imgs(i)%rtsq_serial(algninfo(i)%e3, 0., 0., rmat_rot)
            call imgs(i)%set_rmat(rmat_rot, .false.)
        end do
        !$omp end parallel do
    end subroutine rtsq_imgs

    subroutine calc_cavg_offset( cavg, lp, msk, offset, ind )
        class(image),   intent(in)    :: cavg
        real,           intent(in)    :: lp, msk
        real,           intent(inout) :: offset(2)
        integer, optional, intent(in) :: ind
        type(image_bin)    :: bincavg, bincc
        integer, allocatable :: ccsz(:)
        real    :: hp
        integer :: loc
        if( cavg%is_ft() ) THROW_HARD('Real space only! calc_cavg_offset')
        if( cavg%is_3D() ) THROW_HARD('2D only! calc_cavg_offset')
        offset = 0.
        hp     = max(1.5*lp, 2.*msk*cavg%get_smpd())
        call bincavg%transfer2bimg(cavg)
        ! band-pass
        call bincavg%fft
        call bincavg%bpgau2D(hp, lp)
        call bincavg%ifft
        ! mask
        call bincavg%mask2D_hard(msk)
        ! ignore negative values
        call bincavg%zero_below(0.)
        ! normalize
        call bincavg%norm_minmax
        ! threshold
        call otsu_img(bincavg)
        call bincavg%set_imat
        ! find the largest connected component
        call bincavg%find_ccs(bincc)
        ccsz = bincc%size_ccs()
        loc  = maxloc(ccsz,dim=1)
        ! detemines position of center
        call bincc%masscen_cc(loc, offset)
        ! cleanup
        call bincc%kill_bimg
        call bincavg%kill_bimg
        if( allocated(ccsz) ) deallocate(ccsz)
    end subroutine calc_cavg_offset

    subroutine test_cavg_registration
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        type(inpl_struct), allocatable :: alg_info1(:), alg_info2(:,:)
        type(image),       allocatable :: imgs_ref(:), imgs_targ(:)
        type(parameters)    :: params
        type(cmdline)       :: cline
        real, pointer       :: rmat(:,:,:)
        integer             :: i, j, x, y
        real                :: rx, ry
        integer, parameter  :: NIMGS = 5, BOX = 64
        real,    parameter  :: SMPD = 1.5, CORR_TOL = 0.95
        call cline%set('smpd',  SMPD)
        call cline%set('lp',      6.)
        call cline%set('hp',     20.)
        call cline%set('nthr',     1)
        call cline%set('trs',      8.)
        call cline%set('ctf',    'no')
        call cline%set('objfun', 'cc')
        call cline%set('box',      BOX)
        call cline%set('mskdiam',  48.)
        call params%new(cline)
        allocate(imgs_ref(NIMGS), imgs_targ(NIMGS))
        call imgs_ref(1)%new([BOX,BOX,1], SMPD)
        call imgs_ref(1)%get_rmat_ptr(rmat)
        do y = 1, BOX
            ry = real(y - (BOX / 2 + 1))
            do x = 1, BOX
                rx = real(x - (BOX / 2 + 1))
                rmat(x,y,1) = exp(-((rx - 8.)**2 + (ry + 6.)**2) / 18.) + &
                    &0.7 * exp(-((rx + 7.)**2 + (ry - 4.)**2) / 32.) + &
                    &0.4 * exp(-((rx - 2.)**2 + (ry - 10.)**2) / 8.)
            enddo
        enddo
        call imgs_ref(1)%norm
        call imgs_targ(1)%copy(imgs_ref(1))
        do i = 2, NIMGS
            call imgs_ref(i)%copy(imgs_ref(1))
            call imgs_ref(i)%rtsq(real(i - 1) * 30., 0., 0.)
            call imgs_targ(i)%copy(imgs_ref(i))
        enddo
        alg_info2 = match_imgs(params, params%hp, params%lp, params%trs, imgs_ref, imgs_targ)
        do i = 1, NIMGS
            do j = 1, NIMGS
                if( .not.ieee_is_finite(alg_info2(i,j)%corr) .or. alg_info2(i,j)%corr < CORR_TOL )then
                    THROW_HARD('MATCH_IMGS FAILED')
                endif
            enddo
        enddo
        do i = 2, NIMGS
            call imgs_targ(i)%copy(imgs_ref(1))
            call imgs_targ(i)%rtsq(real(i - 1) * 30., 0.25 * real(i - 1), 0.25 * real(i - 1))
        enddo
        alg_info1 = match_imgs2ref(params, params%hp, params%lp, params%trs, imgs_ref(1), imgs_targ)
        do i = 1, NIMGS
            write(logfhandle,'(A,I0,A,F8.4,A,2F8.3)') 'cavg_registration image ',i, &
                &' corr=',alg_info1(i)%corr,' shift=',alg_info1(i)%x,alg_info1(i)%y
            if( .not.ieee_is_finite(alg_info1(i)%corr) .or. alg_info1(i)%corr < CORR_TOL )then
                THROW_HARD('MATCH_IMGS2REF FAILED')
            endif
        enddo
        call dealloc_imgarr(imgs_ref)
        call dealloc_imgarr(imgs_targ)
        deallocate(alg_info1, alg_info2)
        call cline%kill
    end subroutine test_cavg_registration

end module simple_strategy2D_utils
