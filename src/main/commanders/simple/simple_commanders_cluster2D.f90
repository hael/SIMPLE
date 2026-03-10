!@descr: simultaneous 2D alignment and clustering of single-particle images
module simple_commanders_cluster2D
use simple_commanders_api
use simple_pftc_srch_api
use simple_classaverager
use simple_commanders_cavgs,   only: commander_rank_cavgs
use simple_commanders_mkcavgs, only: commander_make_cavgs, commander_make_cavgs_distr, commander_cavgassemble
use simple_commanders_euclid,  only: commander_calc_pspec_distr, commander_calc_group_sigmas
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

contains

    subroutine exec_cluster2D(self, cline)
        use simple_cluster2D_strategy
        use simple_cluster2D_common
        class(commander_cluster2D), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
         class(cluster2D_strategy), allocatable   :: strategy
        type(parameters)       :: params
        type(builder)          :: build
        type(simple_nice_comm) :: nice_comm
        type(string)           :: finalcavgs
        logical                :: converged
        ! Initialize
        call cline%set('prg','cluster2D')
        call set_cluster2D_defaults(cline)
        call params%new(cline)
        call build%build_spproj(params, cline, wthreads=.true.)
        call build%build_general_tbox(params, cline, do3d=.false.)
        call build%build_strategy2D_tbox(params)
        if( build%spproj%get_nptcls() == 0 ) THROW_HARD('no particles found!')
        call cline%set('mkdir', 'no')
        ! Nice communicator
        call nice_comm%init(params%niceprocid, params%niceserver)
        nice_comm%stat_root%stage = "initialising"
        call nice_comm%cycle()
        if(cline%defined("niceserver")) call cline%delete('niceserver')
        if(cline%defined("niceprocid")) call cline%delete('niceprocid')
        ! This needs to be before init_cluster2D_refs since it uses which_iter
        params%which_iter = max(1, params%startit)
        call init_cluster2D_refs(cline, params, build)
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        ! Create strategy
        strategy = create_cluster2D_strategy(params, cline)
        call strategy%initialize(params, build, cline)
        if( trim(params%stream2d).eq.'no' ) call progressfile_init()
        ! Main loop
        params%which_iter = params%which_iter   - 1
        if( cline%defined('extr_iter') )then
            params%extr_iter = params%extr_iter - 1
        else
            params%extr_iter = params%startit   - 1
        endif
        do
            params%which_iter = params%which_iter + 1
            params%extr_iter  = params%extr_iter  + 1
            write(logfhandle,'(A)')   '>>>'
            write(logfhandle,'(A,I6)')'>>> ITERATION ', params%which_iter
            write(logfhandle,'(A)')   '>>>'
            nice_comm%stat_root%stage = "iteration " // int2str(params%which_iter)
            call nice_comm%cycle()
            ! Strategy handles everything: alignment + cavgs + convergence
            call strategy%execute_iteration(params, build, cline, converged)
            call strategy%finalize_iteration(params, build)
            if( converged .or. params%which_iter >= params%maxits ) exit
        end do
        ! Cleanup
        nice_comm%stat_root%stage = "terminating"
        call nice_comm%cycle()
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params)
        call nice_comm%terminate()
        call simple_end('**** SIMPLE_CLUSTER2D NORMAL STOP ****')
        if (allocated(strategy)) deallocate(strategy)
        call build%kill_general_tbox()
        call build%kill_strategy2D_tbox()
        call simple_touch(CLUSTER2D_FINISHED)
    end subroutine exec_cluster2D

    subroutine exec_cluster2D_distr_worker(self, cline)
        use simple_strategy2D_matcher, only: cluster2D_exec
        class(commander_cluster2D_distr_worker), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        logical          :: converged
        call set_cluster2D_defaults(cline)
        call params%new(cline)
        ! flag subprocess executed through simple_private_exec 
        if( .not. cline%defined('part')    ) THROW_HARD('PART must be defined for distributed worker execution')
        if( .not. cline%defined('outfile') ) THROW_HARD('OUTFILE must be defined for distributed worker execution')
        call build%build_spproj(params, cline, wthreads=.true.)
        call build%build_general_tbox(params, cline, do3d=.false.)
        call build%build_strategy2D_tbox(params)
        params%which_iter = max(1, params%startit)
        params%extr_iter  = params%which_iter
        call cline%set('which_iter', int2str(params%which_iter))
        call cluster2D_exec(params, build, cline, params%which_iter, converged)
        call build%kill_general_tbox()
        call build%kill_strategy2D_tbox()
    end subroutine exec_cluster2D_distr_worker

    subroutine exec_ppca_denoise_classes( self, cline )
        use simple_imgproc,       only: make_pcavecs
        use simple_pca,           only: pca
        use simple_pca_svd,       only: pca_svd
        use simple_kpca_svd,      only: kpca_svd
        use simple_ppca_inmem,    only: ppca_inmem
        class(commander_ppca_denoise_classes), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        integer,          parameter   :: MAXPCAITS = 15
        class(pca),       pointer     :: pca_ptr  => null()
        type(parameters)              :: params
        type(builder)                 :: build
        type(image),      allocatable :: imgs(:), imgs_ori(:)
        type(image)                   :: cavg, img, timg
        type(oris)                    :: os
        type(sp_project), target      :: spproj
        type(string)                  :: label, fname, fname_denoised, fname_cavgs, fname_cavgs_denoised
        type(string)                  :: fname_oris, fname_denoised_ori, fname_ori, fname_class_ptcls_den
        integer,          allocatable :: cls_inds(:), pinds(:), cls_pops(:), ori_map(:)
        real,             allocatable :: avg(:), avg_pix(:), pcavecs(:,:), tmpvec(:)
        real    :: shift(2), loc(2), dist(2), e3, kw, mat(2,2), mat_inv(2,2)
        complex :: fcompl, fcompll
        integer :: npix, i, j, ncls, nptcls, cnt1, cnt2, neigs, h, k, win_corner(2),&
                  &l, ll, m, mm, phys(2), logi_lims(3,2), cyc_lims(3,2), cyc_limsR(2,2), errflg
        logical :: l_phflip, l_transp_pca, l_pre_norm ! pixel-wise learning
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        if( .not. cline%defined('neigs')   ) call cline%set('neigs',    10)
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
        ! pca allocation
        select case(trim(params%pca_mode))
            case('ppca')
                allocate(ppca_inmem :: pca_ptr)
            case('pca_svd')
                allocate(pca_svd    :: pca_ptr)
            case('kpca')
                allocate(kpca_svd   :: pca_ptr)
        end select
        do i = 1, ncls
            call progress_gfortran(i,ncls)
            if( trim(params%pca_img_ori) .eq. 'yes' )then
                call transform_ptcls(params, build, spproj, params%oritype, cls_inds(i), imgs, pinds, phflip=l_phflip, cavg=cavg, imgs_ori=imgs_ori)
                do j = 1, size(imgs)
                    call imgs(j)%copy_fast(imgs_ori(j))
                enddo
            else
                call transform_ptcls(params, build, spproj, params%oritype, cls_inds(i), imgs, pinds, phflip=l_phflip, cavg=cavg)
            endif
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
            if( allocated(tmpvec) ) deallocate(tmpvec)
            if( l_transp_pca )then
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
                    ! call imgs(j)%write(fname_denoised, cnt2)
                    if( trim(params%pca_ori_stk) .eq. 'yes' ) ori_map(pinds(j)) = cnt2
                    call imgs(j)%kill
                end do
                call cavg%div(real(nptcls))
            endif
            call cavg%write(fname_cavgs_denoised, i)
        end do
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
        deallocate(imgs)
        call build%kill_general_tbox
        call os%kill
        ! end gracefully
        call simple_end('**** SIMPLE_PPCA_DENOISE_CLASSES NORMAL STOP ****')
    end subroutine exec_ppca_denoise_classes

end module simple_commanders_cluster2D
