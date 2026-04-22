!@descr: simultaneous 2D alignment and clustering of single-particle images
module simple_commanders_cluster2D
use simple_commanders_api
use simple_pftc_srch_api
use simple_classaverager
use simple_commanders_cavgs,   only: commander_rank_cavgs
use simple_commanders_mkcavgs, only: commander_make_cavgs, commander_make_cavgs_distr, commander_cavgassemble
use simple_commanders_euclid,  only: commander_calc_pspec, commander_calc_group_sigmas
use simple_gui_utils,          only: mrc2jpeg_tiled
use simple_procimgstk,         only: selection_from_tseries_imgfile, random_selection_from_imgfile, copy_imgfile, noise_imgfile
use simple_progress,           only: progressfile_init, progressfile_update
use simple_commanders_imgops,  only: commander_scale
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_cluster2D
  contains
    procedure :: execute      => exec_cluster2D
end type commander_cluster2D

type, extends(commander_base) :: commander_cluster2D_distr_worker
  contains
    procedure :: execute      => exec_cluster2D_distr_worker
end type commander_cluster2D_distr_worker

type, extends(commander_base) :: commander_ppca_denoise_classes
  contains
    procedure :: execute      => exec_ppca_denoise_classes
end type commander_ppca_denoise_classes

type, extends(commander_base) :: commander_ppca_cls_split
  contains
    procedure :: execute      => exec_ppca_cls_split
end type commander_ppca_cls_split

type, extends(commander_base) :: commander_diffusion_cls_split
  contains
    procedure :: execute      => exec_diffusion_cls_split
end type commander_diffusion_cls_split

contains

    !> Single entrypoint (shared-memory OR distributed master), driven by a strategy.
    subroutine exec_cluster2D( self, cline )
        use simple_cluster2D_strategy
        class(commander_cluster2D), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        class(cluster2D_strategy), allocatable    :: strategy
        type(parameters)       :: params
        type(builder)          :: build
        type(simple_nice_comm) :: nice_comm
        logical                :: converged
        integer                :: niters
        ! local defaults (kept consistent with previous distributed master)
        call cline%set('prg', 'cluster2D')
        call set_cluster2D_defaults(cline)
        if( cline%defined('polar') )then
            if( cline%get_carg('polar') .ne. 'no' )then
                THROW_HARD('polar=yes is no longer supported for 2D workflows; use cartesian cluster2D')
            endif
        endif
        call cline%set('polar', 'no')
        ! Select execution strategy (shared-memory vs distributed master)
        strategy = create_cluster2D_strategy(cline)
        call strategy%initialize(params, build, cline)
        ! Nice communicator
        call nice_comm%init(params%niceprocid, params%niceserver)
        nice_comm%stat_root%stage = "initialising"
        call nice_comm%cycle()
        if( cline%defined("niceserver") ) call cline%delete('niceserver')
        if( cline%defined("niceprocid") ) call cline%delete('niceprocid')
        if( trim(params%stream2d).eq.'no' ) call progressfile_init()
        ! Main loop counter semantics:
        !   - params%maxits is the *number of iterations to run* in this invocation.
        !   - params%which_iter starts at params%startit.
        niters            = 0
        params%which_iter = params%startit - 1
        if( cline%defined('extr_iter') )then
            params%extr_iter = params%extr_iter - 1
        else
            params%extr_iter = params%startit   - 1
        endif
        do
            niters            = niters + 1
            params%which_iter = params%which_iter + 1
            params%extr_iter  = params%extr_iter  + 1
            nice_comm%stat_root%stage = "iteration " // int2str(params%which_iter)
            call nice_comm%cycle()
            ! Strategy handles everything: alignment + cavgs + convergence
            call strategy%execute_iteration(params, build, cline, converged)
            call strategy%finalize_iteration(params, build)
            if( converged .or. niters >= params%maxits ) exit
        end do
        ! Cleanup
        nice_comm%stat_root%stage = "terminating"
        call nice_comm%cycle()
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params)
        call nice_comm%terminate()
        if( allocated(strategy) ) deallocate(strategy)
        ! Global teardown (strategies may have built different toolboxes)
        call build%kill_general_tbox()
        call build%kill_strategy2D_tbox()
        call simple_touch(CLUSTER2D_FINISHED)
        call simple_end('**** SIMPLE_CLUSTER2D NORMAL STOP ****')
    end subroutine exec_cluster2D

    !> Distributed worker (single-iteration execution). This should be the command
    !> invoked by the scheduler for each partition.
    subroutine exec_cluster2D_distr_worker( self, cline )
        use simple_strategy2D_matcher, only: cluster2D_exec
        class(commander_cluster2D_distr_worker), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        logical          :: converged
        call set_cluster2D_defaults(cline)
        if( cline%defined('polar') )then
            if( cline%get_carg('polar') .ne. 'no' )then
                THROW_HARD('polar=yes is no longer supported for 2D workflows; use cartesian cluster2D')
            endif
        endif
        call cline%set('polar', 'no')
        ! Flags required for worker execution
        if( .not. cline%defined('part')    ) THROW_HARD('PART must be defined for distributed worker execution')
        if( .not. cline%defined('outfile') ) THROW_HARD('OUTFILE must be defined for distributed worker execution')
        ! Worker needs the alignment toolboxes
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
        if( cline%defined('which_iter') )then
            params%which_iter = max(1, params%which_iter)
        else
            params%which_iter = max(1, params%startit)
        endif
        if( .not. cline%defined('extr_iter') ) params%extr_iter = params%which_iter
        call cline%set('which_iter', int2str(params%which_iter))
        call cluster2D_exec(params, build, cline, params%which_iter, converged)
        call build%kill_strategy2D_tbox
        call build%kill_general_tbox
    end subroutine exec_cluster2D_distr_worker

    subroutine exec_ppca_denoise_classes( self, cline )
        use simple_imgarr_utils, only: alloc_imgarr, dealloc_imgarr
        use simple_imgproc,       only: make_pcavecs
        use simple_pca,           only: pca
        use simple_pca_svd,       only: pca_svd
        use simple_kpca_svd,      only: kpca_svd, suggest_kpca_nystrom_neigs
        use simple_ppca,          only: ppca
        use simple_mppca,         only: mppca
        use simple_strategy2D_utils, only: prep_cavgs4clust, calc_cluster_cavgs_dmat
        class(commander_ppca_denoise_classes), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        integer,          parameter   :: MAXPCAITS = 15
        ! Denoising-by-class is intentionally biased toward low-rank PPCA. A broad rank
        ! search tends to select pathological edge models (q=1 or q~n-1), which is both
        ! slow and visually inconsistent. Keep the reusable PPCA selector general, but
        ! hand this workflow a small denoising-oriented candidate grid.
        integer,          parameter   :: PPCA_AUTO_NCAND = 10
        integer,          parameter   :: PPCA_AUTO_CAND(PPCA_AUTO_NCAND) = [1,2,3,4,5,6,8,10,12,16]
        class(pca),       pointer     :: pca_ptr  => null()
        type(parameters)              :: params
        type(builder)                 :: build
        type(ppca)                    :: ppca_rank_selector
        type(ppca),       allocatable :: ppca_local_models(:)
        type(image),      allocatable :: imgs(:), imgs_ori(:)
        type(image),      allocatable :: mix_cavg_imgs(:)
        type(image)                   :: cavg, img, timg
        type(oris)                    :: os
        type(sp_project), target      :: spproj
        type(string)                  :: label, fname, fname_denoised, fname_cavgs, fname_cavgs_denoised
        type(string)                  :: fname_oris, fname_denoised_ori, fname_ori, fname_class_ptcls_den
        integer,          allocatable :: cls_inds(:), pinds(:), cls_pops(:), ori_map(:), cls_inds_mix(:), clspops_mix(:), mix_neigh(:)
        logical,          allocatable :: mix_mask(:)
        real,             allocatable :: avg(:), avg_pix(:), pcavecs(:,:), tmpvec(:)
        integer,          allocatable :: rank_scan_qs(:)
        real(dp),         allocatable :: rank_scan_bics(:), rank_scan_sigma2(:)
        real,             allocatable :: avg_store(:,:), rawvec(:), selfvec(:), neighvec(:), blendvec(:)
        real,             allocatable :: mix_dmat(:,:), mm_mix(:,:)
        real    :: shift(2), loc(2), dist(2), e3, kw, mat(2,2), mat_inv(2,2)
        complex :: fcompl, fcompll
        integer :: npix, i, j, ncls, nptcls, cnt1, cnt2, neigs, h, k, win_corner(2),&
                  &l, ll, m, mm, phys(2), logi_lims(3,2), cyc_lims(3,2), cyc_limsR(2,2), errflg
        logical :: l_phflip, l_transp_pca, l_pre_norm, l_ppca_local_mix ! pixel-wise learning
        integer :: ineigh
        integer :: mix_ldim(3)
        real    :: oa_min, oa_max, mix_w_neigh, mix_w_self
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        if( .not. cline%defined('neigs')   ) call cline%set('neigs',    10)
        if( cline%get_carg('pca_mode') .eq. 'ppca_local_mix' )then
            if( .not. cline%defined('objfun') ) call cline%set('objfun', 'cc')
            if( .not. cline%defined('trs')    ) call cline%set('trs',    10.)
            if( .not. cline%defined('lp')     ) call cline%set('lp',      6.)
        endif
        call build%init_params_and_build_general_tbox(cline, params, do3d=(trim(params%oritype) .eq. 'ptcl3D'))
        call spproj%read(params%projfile)
        select case(trim(params%oritype))
            case('ptcl2D')
                label = 'class'
            case('ptcl3D')
                label = 'proj'
                call build%spproj_field%proj2class
            case DEFAULT
                THROW_HARD('ORITYPE not supported!')
        end select
        l_transp_pca = (trim(params%transp_pca) .eq. 'yes')
        l_pre_norm   = (trim(params%pre_norm)   .eq. 'yes')
        l_ppca_local_mix = trim(params%pca_mode) .eq. 'ppca_local_mix'
        if( l_ppca_local_mix )then
            if( trim(params%pca_img_ori) .eq. 'yes' ) THROW_HARD('ppca_local_mix currently supports pca_img_ori=no only')
            if( l_transp_pca ) THROW_HARD('ppca_local_mix currently supports transp_pca=no only')
        endif
        l_phflip     = .false.
        select case( spproj%get_ctfflag_type(params%oritype) )
            case(CTFFLAG_NO)
                THROW_WARN('No CTF information could be found, phase flipping is deactivated')
            case(CTFFLAG_FLIP)
                THROW_WARN('Images have already been phase-flipped, phase flipping is deactivated')
            case(CTFFLAG_YES)
                l_phflip = .true.
            case DEFAULT
                THROW_HARD('UNSUPPORTED CTF FLAG')
        end select
        cls_inds = build%spproj_field%get_label_inds(label%to_char())
        if( cline%defined('class') .and. cline%defined('ncls') )then
            THROW_HARD('EITHER class OR ncls CAN BE DEFINED')
        endif
        if( cline%defined('class') )then
            cls_inds = pack(cls_inds, mask=(cls_inds == params%class))
        endif
        ncls = size(cls_inds)
        nptcls = 0
        if( cline%defined('ncls') )then
            ncls     = params%ncls
            cls_inds = cls_inds(1:ncls)
        endif
        allocate(cls_pops(ncls), source=0)
        do i = 1, ncls
            call build%spproj_field%get_pinds(cls_inds(i), label%to_char(), pinds)
            if( allocated(pinds) )then
                cls_pops(i) = size(pinds)
                nptcls = nptcls + cls_pops(i)
                deallocate(pinds)
            endif
        end do
        cls_inds = pack(cls_inds, mask=cls_pops > 2)
        nptcls   = sum(cls_pops,  mask=cls_pops > 2)
        ncls     = size(cls_inds)
        if( trim(params%pca_ori_stk) .eq. 'yes' ) allocate(ori_map(nptcls))
        call os%new(nptcls, is_ptcl=.true.)
        fname                = 'ptcls.mrcs'
        fname_denoised       = 'ptcls_denoised.mrcs'
        fname_cavgs          = 'cavgs.mrcs'
        fname_cavgs_denoised = 'cavgs_denoised.mrcs'
        fname_oris           = 'oris_denoised.txt'
        fname_ori            = 'ptcls_ori_order.mrcs'
        fname_denoised_ori   = 'ptcls_denoised_ori_order.mrcs'
        cnt1 = 0
        cnt2 = 0
        if( l_ppca_local_mix )then
            allocate(mix_mask(spproj%os_cls2D%get_noris()), source=.false.)
            mix_mask(cls_inds) = .true.
            call prep_cavgs4clust(spproj, mix_cavg_imgs, params%mskdiam, clspops_mix, cls_inds_mix, mix_mask, mm_mix)
            if( size(cls_inds_mix) /= ncls ) THROW_HARD('ppca_local_mix class subset mismatch after cavg preparation')
            if( any(cls_inds_mix /= cls_inds) ) THROW_HARD('ppca_local_mix expected class ordering to match current class subset')
            mix_ldim = mix_cavg_imgs(1)%get_ldim()
            params%smpd = mix_cavg_imgs(1)%get_smpd()
            params%box  = mix_ldim(1)
            params%msk  = min(real(params%box/2) - COSMSKHALFWIDTH - 1., 0.5 * params%mskdiam / params%smpd)
            oa_min = minval(mm_mix(:,1))
            oa_max = maxval(mm_mix(:,2))
            mix_dmat = calc_cluster_cavgs_dmat(params, mix_cavg_imgs, [oa_min,oa_max], params%clust_crit)
            allocate(mix_neigh(ncls), source=0)
            do i = 1, ncls
                ineigh = 0
                do j = 1, ncls
                    if( i == j ) cycle
                    if( ineigh == 0 )then
                        ineigh = j
                    else if( mix_dmat(i,j) < mix_dmat(i,ineigh) )then
                        ineigh = j
                    endif
                enddo
                mix_neigh(i) = ineigh
                if( ineigh > 0 )then
                    write(logfhandle,'(A,I8,A,I8,A,F7.4)') 'PPCA local mix neighborhood: class=', cls_inds(i), ' neigh=', cls_inds(ineigh), ' d=', mix_dmat(i,ineigh)
                    call flush(logfhandle)
                endif
            enddo
            allocate(ppca_local_models(ncls))
        endif
        ! pca allocation
        select case(trim(params%pca_mode))
            case('ppca','ppca_local_mix')
                allocate(ppca :: pca_ptr)
            case('mppca')
                allocate(mppca :: pca_ptr)
            case('pca_svd')
                allocate(pca_svd    :: pca_ptr)
            case('kpca')
                allocate(kpca_svd   :: pca_ptr)
        end select
        select type(pca_ptr)
            type is(mppca)
                call pca_ptr%set_params(params%mppca_k, params%nthr)
            type is(kpca_svd)
                call pca_ptr%set_params(params%nthr, params%kpca_ker, params%kpca_backend,&
                    &params%kpca_nystrom_npts, params%kpca_rbf_gamma, params%kpca_nystrom_local_nbrs, params%kpca_cosine_weight_power)
        end select
        do i = 1, ncls
            call progress_gfortran(i,ncls)
            write(logfhandle,'(A,I8,A,I8)') 'ppca_denoise_classes entering class loop: iclass=', i, ' cls_id=', cls_inds(i)
            call flush(logfhandle)
            if( trim(params%pca_img_ori) .eq. 'yes' )then
                call transform_ptcls(params, build, spproj, params%oritype, cls_inds(i), imgs, pinds, phflip=l_phflip, cavg=cavg, imgs_ori=imgs_ori)
                do j = 1, size(imgs)
                    call imgs(j)%copy_fast(imgs_ori(j))
                enddo
            else
                call transform_ptcls(params, build, spproj, params%oritype, cls_inds(i), imgs, pinds, phflip=l_phflip, cavg=cavg)
            endif
            write(logfhandle,'(A,I8,A,I8)') 'ppca_denoise_classes transform_ptcls done: iclass=', i, ' nimgs=', size(imgs)
            call flush(logfhandle)
            nptcls = size(imgs)
            if( trim(params%neigs_per).eq.'yes' )then
                if( params%neigs >= 99 )then
                    THROW_WARN('neigs is greater than 99% the number of particles within this class. All eigens are used now!')
                    neigs = nptcls - 1
                else
                    neigs = max(2, nint(real(params%neigs * nptcls) / 100.))
                endif
            else
                neigs = params%neigs
                if( neigs >= nptcls )then
                    THROW_WARN('neigs is greater than the number of particles within this class. All eigens are used now!')
                    neigs = nptcls - 1
                endif
            endif
            if( l_pre_norm )then
                do j = 1, nptcls
                    call imgs(j)%norm
                end do
            endif
            do j = 1, nptcls
                cnt1 = cnt1 + 1
                call imgs(j)%write(fname, cnt1)
            end do
            call cavg%write(fname_cavgs, i)
            ! performs ppca
            if( trim(params%projstats).eq.'yes' )then
                call make_pcavecs(imgs, npix, avg, pcavecs, transp=l_transp_pca, avg_pix=avg_pix)
            else
                call make_pcavecs(imgs, npix, avg, pcavecs, transp=l_transp_pca)
            endif
            if( (trim(params%pca_mode) .eq. 'ppca' .or. l_ppca_local_mix) .and. params%neigs <= 0 )then
                if( allocated(rank_scan_qs) )     deallocate(rank_scan_qs)
                if( allocated(rank_scan_bics) )   deallocate(rank_scan_bics)
                if( allocated(rank_scan_sigma2) ) deallocate(rank_scan_sigma2)
                neigs = ppca_rank_selector%suggest_rank(pcavecs, PPCA_AUTO_CAND, MAXPCAITS, rank_scan_qs, rank_scan_bics, rank_scan_sigma2)
                write(logfhandle,'(A,I8,A,I8,A,I8)') 'PPCA classes auto-selected neigs: class=', cls_inds(i), ' size=', nptcls, ' neigs=', neigs
                call flush(logfhandle)
                call log_ppca_rank_scan(cls_inds(i), nptcls, neigs, rank_scan_qs, rank_scan_bics, rank_scan_sigma2)
            endif
            if( trim(params%pca_mode) .eq. 'kpca' .and. trim(params%kpca_backend) .eq. 'nystrom' .and. neigs <= 0 )then
                neigs = suggest_kpca_nystrom_neigs(pcavecs, params%kpca_ker, params%kpca_nystrom_npts, params%kpca_rbf_gamma)
                write(logfhandle,'(A,I8,A,I8,A)') 'kPCA classes auto-selected neigs: ', neigs, ' for class size ', nptcls, &
                    ' (Nyström spectrum 99% energy)'
                call flush(logfhandle)
            endif
            if( l_transp_pca )then
                neigs = min(max(neigs, 1), max(npix-1, 1))
            else
                neigs = min(max(neigs, 1), max(nptcls-1, 1))
            endif
            if( allocated(tmpvec) ) deallocate(tmpvec)
            if( l_ppca_local_mix )then
                if( allocated(avg_store) )then
                    if( size(avg_store,1) /= npix )then
                        deallocate(avg_store)
                        allocate(avg_store(npix,ncls), source=0.)
                    endif
                else
                    allocate(avg_store(npix,ncls), source=0.)
                endif
                avg_store(:,i) = avg
                call ppca_local_models(i)%new(nptcls, npix, neigs)
                call ppca_local_models(i)%set_verbose(.false.)
                call ppca_local_models(i)%master(pcavecs, MAXPCAITS)
                call ppca_local_models(i)%slim()
                cycle
            else if( l_transp_pca )then
                call pca_ptr%new(npix, nptcls, neigs)
                call pca_ptr%master(pcavecs, MAXPCAITS)
                allocate(tmpvec(nptcls))
                !$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
                do j = 1, npix
                    call pca_ptr%generate(j, avg, tmpvec)
                    pcavecs(:,j) = tmpvec
                end do
                !$omp end parallel do
                pcavecs = transpose(pcavecs)
            else
                call pca_ptr%new(nptcls, npix, neigs)
                call pca_ptr%master(pcavecs, MAXPCAITS)
                allocate(tmpvec(npix))
                !$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
                do j = 1, nptcls
                    call pca_ptr%generate(j, avg, tmpvec)
                    pcavecs(:,j) = tmpvec
                end do
                !$omp end parallel do
            endif
            if( trim(params%projstats).eq.'yes' )then
                call cavg%unserialize(avg_pix)
                call cavg%write(string('cavgs_unserialized.mrcs'), i)
            endif
            ! output
            call cavg%zero_and_unflag_ft
            if( trim(params%pca_img_ori) .eq. 'yes' )then
                do j = 1, nptcls
                    cnt2 = cnt2 + 1
                    call imgs_ori(j)%unserialize(pcavecs(:,j))
                    call os%transfer_ori(cnt2, build%spproj_field, pinds(j))
                    call imgs_ori(j)%write(fname_denoised, cnt2)
                    if( trim(params%pca_ori_stk) .eq. 'yes' ) ori_map(pinds(j)) = cnt2
                end do
                call transform_ptcls(params, build, spproj, params%oritype, cls_inds(i), imgs, pinds, phflip=l_phflip, cavg=cavg, imgs_ori=imgs_ori)
            else
                fname_class_ptcls_den = 'class'//int2str_pad(i,4)//'ptcls.mrcs'
                do j = 1, nptcls
                    cnt2 = cnt2 + 1
                    call imgs(j)%unserialize(pcavecs(:,j))
                    call cavg%add(imgs(j))
                    call os%transfer_ori(cnt2, build%spproj_field, pinds(j))
                    call imgs(j)%write(fname_class_ptcls_den, j)
                    call imgs(j)%write(fname_denoised, cnt2)
                    if( trim(params%pca_ori_stk) .eq. 'yes' ) ori_map(pinds(j)) = cnt2
                    call imgs(j)%kill
                end do
                call cavg%div(real(nptcls))
            endif
            call cavg%write(fname_cavgs_denoised, i)
        end do
        if( l_ppca_local_mix )then
            write(logfhandle,'(A)') 'PPCA local mix reconstruct/write pass'
            call flush(logfhandle)
            do i = 1, ncls
                call progress_gfortran(i,ncls)
                call transform_ptcls(params, build, spproj, params%oritype, cls_inds(i), imgs, pinds, phflip=l_phflip, cavg=cavg)
                nptcls = size(imgs)
                if( l_pre_norm )then
                    do j = 1, nptcls
                        call imgs(j)%norm
                    end do
                endif
                if( trim(params%projstats).eq.'yes' )then
                    call make_pcavecs(imgs, npix, avg, pcavecs, transp=l_transp_pca, avg_pix=avg_pix)
                else
                    call make_pcavecs(imgs, npix, avg, pcavecs, transp=l_transp_pca)
                endif
                if( .not. allocated(rawvec) ) then
                    allocate(rawvec(npix), selfvec(npix), neighvec(npix), blendvec(npix), source=0.)
                else if( size(rawvec) /= npix ) then
                    deallocate(rawvec, selfvec, neighvec, blendvec)
                    allocate(rawvec(npix), selfvec(npix), neighvec(npix), blendvec(npix), source=0.)
                endif
                ineigh = mix_neigh(i)
                call cavg%zero_and_unflag_ft
                fname_class_ptcls_den = 'class'//int2str_pad(i,4)//'ptcls.mrcs'
                do j = 1, nptcls
                    rawvec = avg + pcavecs(:,j)
                    call ppca_local_models(i)%reconstruct_external(pcavecs(:,j), selfvec)
                    blendvec = avg + selfvec
                    if( ineigh > 0 )then
                        neighvec = rawvec - avg_store(:,ineigh)
                        call ppca_local_models(ineigh)%reconstruct_external(neighvec, neighvec)
                        mix_w_neigh = 0.20 * max(0., 1.0 - mix_dmat(i,ineigh))
                        mix_w_self  = 1.0 - mix_w_neigh
                        blendvec = mix_w_self * blendvec + mix_w_neigh * (avg_store(:,ineigh) + neighvec)
                    endif
                    cnt2 = cnt2 + 1
                    call imgs(j)%unserialize(blendvec)
                    call cavg%add(imgs(j))
                    call os%transfer_ori(cnt2, build%spproj_field, pinds(j))
                    call imgs(j)%write(fname_class_ptcls_den, j)
                    call imgs(j)%write(fname_denoised, cnt2)
                    if( trim(params%pca_ori_stk) .eq. 'yes' ) ori_map(pinds(j)) = cnt2
                    call imgs(j)%kill
                enddo
                call cavg%div(real(nptcls))
                call cavg%write(fname_cavgs_denoised, i)
            enddo
        endif
        if( trim(params%pca_ori_stk) .eq. 'yes' )then
            call  img%new([params%boxpd,params%boxpd,1],params%smpd, wthreads=.false.)
            call timg%new([params%boxpd,params%boxpd,1],params%smpd, wthreads=.false.)
            logi_lims      = img%loop_lims(2)
            cyc_lims       = img%loop_lims(3)
            cyc_limsR(:,1) = cyc_lims(1,:)
            cyc_limsR(:,2) = cyc_lims(2,:)
            do i = 1, cnt2
                shift = build%spproj_field%get_2Dshift(i)
                e3    = build%spproj_field%e3get(i)
                do j = 1, 2
                    call  img%zero_and_flag_ft
                    call timg%zero_and_flag_ft
                    call cavg%zero_and_unflag_ft
                    if( j == 1 )then
                        call cavg%read(fname_denoised, ori_map(i))
                    else
                        call cavg%read(fname,          ori_map(i))
                    endif
                    call cavg%pad_fft(img)
                    ! particle rotation
                    call rotmat2d(-e3, mat)
                    call matinv(mat, mat_inv, 2, errflg)
                    !$omp parallel do collapse(2) private(h,k,loc,win_corner,dist,l,ll,m,mm,phys,kw,fcompl,fcompll) default(shared) proc_bind(close) schedule(static)
                    do h = logi_lims(1,1),logi_lims(1,2)
                        do k = logi_lims(2,1),logi_lims(2,2)
                            ! Rotation
                            loc        = matmul(real([h,k]),mat_inv)
                            win_corner = floor(loc) ! bottom left corner
                            dist       = loc - real(win_corner)
                            ! Bi-linear interpolation
                            l      = cyci_1d(cyc_limsR(:,1), win_corner(1))
                            ll     = cyci_1d(cyc_limsR(:,1), win_corner(1)+1)
                            m      = cyci_1d(cyc_limsR(:,2), win_corner(2))
                            mm     = cyci_1d(cyc_limsR(:,2), win_corner(2)+1)
                            ! l, bottom left corner
                            phys   = img%comp_addr_phys(l,m)
                            kw     = (1.-dist(1))*(1.-dist(2))   ! interpolation kernel weight
                            fcompl = kw * img%get_cmat_at(phys(1), phys(2),1)
                            ! l, bottom right corner
                            phys   = img%comp_addr_phys(l,mm)
                            kw     = (1.-dist(1))*dist(2)
                            fcompl = fcompl + kw * img%get_cmat_at(phys(1), phys(2),1)
                            if( l < 0 ) fcompl = conjg(fcompl) ! conjugation when required!
                            ! ll, upper left corner
                            phys    = img%comp_addr_phys(ll,m)
                            kw      = dist(1)*(1.-dist(2))
                            fcompll = kw * img%get_cmat_at(phys(1), phys(2),1)
                            ! ll, upper right corner
                            phys    = img%comp_addr_phys(ll,mm)
                            kw      = dist(1)*dist(2)
                            fcompll = fcompll + kw * img%get_cmat_at(phys(1), phys(2),1)
                            if( ll < 0 ) fcompll = conjg(fcompll) ! conjugation when required!
                            ! update with interpolated values
                            phys = img%comp_addr_phys(h,k)
                            call timg%set_cmat_at(phys(1),phys(2),1, fcompl + fcompll)
                        end do
                    end do
                    !$omp end parallel do
                    ! shift
                    call timg%shift2Dserial(shift)
                    call timg%ifft
                    call timg%clip(cavg)
                    if( j == 1 )then
                        call cavg%write(fname_denoised_ori, i)
                    else
                        call cavg%write(fname_ori, i)
                    endif
                enddo
            enddo
        endif
        call os%zero_inpl
        call os%write(fname_oris)
        if( trim(params%projstats).eq.'yes' ) call build%spproj_field%write(string('ptcl_field.txt'))
        ! cleanup
        if( allocated(ppca_local_models) )then
            do i = 1, size(ppca_local_models)
                call ppca_local_models(i)%kill()
            enddo
            deallocate(ppca_local_models)
        endif
        if( allocated(mix_cavg_imgs) ) call dealloc_imgarr(mix_cavg_imgs)
        deallocate(imgs)
        call build%kill_general_tbox
        call os%kill
        ! end gracefully
        call simple_end('**** SIMPLE_PPCA_DENOISE_CLASSES NORMAL STOP ****')

    contains

        subroutine log_ppca_rank_scan(cls_id, cls_size, selected_q, qs, bics, sigma2s)
            integer,  intent(in) :: cls_id, cls_size, selected_q
            integer,  intent(in) :: qs(:)
            real(dp), intent(in) :: bics(:), sigma2s(:)
            real(dp), parameter :: SIGMA_FLOOR_WARN = 1.1e-8_dp
            logical, allocatable :: used(:)
            integer :: k, i, idx
            real(dp) :: best_bic, delta_bic, selected_sigma2
            if( size(qs) == 0 ) return
            best_bic = huge(1._dp)
            do i = 1, size(qs)
                if( qs(i) <= 0 ) cycle
                if( bics(i) < best_bic ) best_bic = bics(i)
            enddo
            if( best_bic >= huge(1._dp) / 2._dp ) return
            allocate(used(size(qs)), source=.false.)
            do k = 1, min(3, count(qs > 0))
                idx = 0
                do i = 1, size(qs)
                    if( qs(i) <= 0 .or. used(i) ) cycle
                    if( idx == 0 )then
                        idx = i
                    else if( bics(i) < bics(idx) )then
                        idx = i
                    endif
                enddo
                if( idx == 0 ) exit
                used(idx) = .true.
                delta_bic = bics(idx) - best_bic
                write(logfhandle,'(A,I8,A,I8,A,I1,A,I8,A,ES10.3,A,ES10.3)') &
                    'PPCA rank scan top: class=', cls_id, ' size=', cls_size, ' rank=', k, ' q=', qs(idx), ' dBIC=', delta_bic, ' sigma2=', sigma2s(idx)
                call flush(logfhandle)
            enddo
            deallocate(used)
            selected_sigma2 = -1._dp
            do i = 1, size(qs)
                if( qs(i) == selected_q )then
                    selected_sigma2 = sigma2s(i)
                    exit
                endif
            enddo
            if( selected_q >= max(1, cls_size - 1) )then
                write(logfhandle,'(A,I8,A,I8,A)') 'PPCA rank scan warning: class=', cls_id, ' selected q hit class ceiling at ', selected_q, ''
                call flush(logfhandle)
            endif
            if( selected_sigma2 > 0._dp .and. selected_sigma2 <= SIGMA_FLOOR_WARN )then
                write(logfhandle,'(A,I8,A,ES10.3,A)') 'PPCA rank scan warning: class=', cls_id, ' selected sigma2 is near floor: ', selected_sigma2, ''
                call flush(logfhandle)
            endif
        end subroutine log_ppca_rank_scan

    end subroutine exec_ppca_denoise_classes

    subroutine exec_ppca_cls_split( self, cline )
        use simple_ppca_cls_split_strategy
        class(commander_ppca_cls_split), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        class(ppca_cls_split_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        strategy = create_ppca_cls_split_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params)
        if( allocated(strategy) ) deallocate(strategy)
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_PPCA_CLS_SPLIT NORMAL STOP ****')
    end subroutine exec_ppca_cls_split

    subroutine exec_diffusion_cls_split( self, cline )
        use simple_diffusion_cls_split_strategy
        class(commander_diffusion_cls_split), intent(inout) :: self
        class(cmdline),                             intent(inout) :: cline
        class(diffusion_cls_split_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        strategy = create_diffusion_cls_split_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params)
        if( allocated(strategy) ) deallocate(strategy)
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_DIFFUSION_CLS_SPLIT NORMAL STOP ****')
    end subroutine exec_diffusion_cls_split

end module simple_commanders_cluster2D
