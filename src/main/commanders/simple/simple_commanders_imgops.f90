!@descr: standard image operations: binarize, filter, denoise, normalize, scale etc.
module simple_commanders_imgops
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_binarize
  contains
    procedure :: execute      => exec_binarize
end type commander_binarize

type, extends(commander_base) :: commander_edge_detect
  contains
    procedure :: execute      => exec_edge_detect
end type commander_edge_detect

type, extends(commander_base) :: commander_filter
  contains
    procedure :: execute      => exec_filter
end type commander_filter

type, extends(commander_base) :: commander_ppca_denoise
  contains
    procedure :: execute      => exec_ppca_denoise
end type commander_ppca_denoise

type, extends(commander_base) :: commander_ppca_denoise_polarft_lines
  contains
    procedure :: execute      => exec_ppca_denoise_polarft_lines
end type commander_ppca_denoise_polarft_lines

type, extends(commander_base) :: commander_normalize
  contains
    procedure :: execute      => exec_normalize
end type commander_normalize

type, extends(commander_base) :: commander_scale
  contains
    procedure :: execute      => exec_scale
end type commander_scale

contains

    !> for binarization of stacks and volumes
    subroutine exec_binarize( self, cline )
        use simple_segmentation
        class(commander_binarize), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(parameters) :: params
        type(image_bin)  :: img_or_vol
        type(image)      :: img_sdevs
        integer :: iptcl
        logical :: is2D, fill_holes, l_sauvola
        ! error check
        if( .not. cline%defined('stk') .and. .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! stk or vol1 needs to be present; binarize')
        endif
        if( cline%defined('stk') .and. cline%defined('vol1') )then
            THROW_HARD('either stk or vol1 key can be present, not both; binarize')
        endif
        if( cline%defined('thres') .and. cline%defined('npix') )then
            THROW_HARD('either thres-based or npix-based binarisation; both keys cannot be present; binarize')
        endif
        call params%new(cline)
        fill_holes = .false.
        if( cline%defined('fill_holes') )then
            if( params%fill_holes .eq. 'yes' ) fill_holes = .true.
        endif
        l_sauvola = .false.
        if( cline%defined('winsz') ) l_sauvola = .true.
        if( cline%defined('stk') )then
            is2D = .true.
            call img_or_vol%new_bimg([params%box,params%box,1], params%smpd)
            do iptcl=1,params%nptcls
                call img_or_vol%read(params%stk, iptcl)
                call doit
                if( fill_holes ) call img_or_vol%set_edgecc2background
                if( l_sauvola  ) call img_sdevs%write(string('local_sdevs.mrc'), iptcl)
                call img_or_vol%write(params%outstk, iptcl)
            end do
        else if( cline%defined('vol1') )then
            is2D = .false.
            call img_or_vol%new_bimg([params%box,params%box,params%box], params%smpd)
            call img_or_vol%read(params%vols(1))
            call doit
            call img_or_vol%write(params%outvol)
        endif
        call img_or_vol%kill
        ! end gracefully
        call simple_end('**** SIMPLE_BINARIZE NORMAL STOP ****')

        contains

            subroutine doit
                type(image_bin) :: cos_img
                real :: ave, sdev, maxv, minv
                if( cline%defined('thres') )then
                    call img_or_vol%binarize(params%thres)
                else if( cline%defined('npix') )then
                    call img_or_vol%binarize(params%npix)
                else if( cline%defined('ndev') )then
                    call img_or_vol%stats( ave, sdev, maxv, minv)
                    call img_or_vol%binarize(ave + params%ndev * sdev)
                else
                    if( cline%defined('winsz') )then
                        call sauvola(img_or_vol, nint(params%winsz), img_sdevs, params%nsig)
                    else
                        if( trim(params%omit_neg) .eq. 'yes' )then
                            call otsu_img(img_or_vol)
                        else
                            call otsu_img(img_or_vol)
                        endif
                    endif
                endif
                write(logfhandle,'(a,1x,i9)') '# FOREGROUND PIXELS:', img_or_vol%nforeground()
                write(logfhandle,'(a,1x,i9)') '# BACKGROUND PIXELS:', img_or_vol%nbackground()
                if( cline%defined('grow') ) call img_or_vol%grow_bins(params%grow)
                if( cline%defined('edge') ) then
                    call img_or_vol%cos_edge(params%edge,cos_img)
                    call img_or_vol%copy_bimg(cos_img)
                    call cos_img%kill
                endif
                if( cline%defined('neg')  ) call img_or_vol%inv_bimg
            end subroutine doit

    end subroutine exec_binarize

    !> for edge detection of stacks
    subroutine exec_edge_detect( self, cline )
        class(commander_edge_detect), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer          :: iptcl
        ! error check
        if( .not. cline%defined('stk') )then
            THROW_HARD('ERROR! stk needs to be present; exec_edge_detect')
        endif
        if( .not. cline%defined('detector') )then
            THROW_HARD('ERROR! detector needs to be present; exec_edge_detect')
        endif
        if( .not. cline%defined('automatic') )then
            THROW_HARD('ERROR! automatic needs to be present; exec_edge_detect')
        endif
        if(.not. cline%defined('outstk') )then
            params%outstk = 'outstk.mrc'
        endif
        if( cline%defined('thres') .and. cline%defined('npix') )then
            THROW_HARD('either thres-based or npix-based edge detection; both keys cannot be present; exec_edge_detect')
        endif
        if( cline%defined('thres') .and. params%automatic .eq. 'yes') then
            THROW_HARD('cannot chose thres in automatic mode; exec_edge_detect')
        endif
        if( cline%defined('npix') .and. params%automatic .eq. 'yes') then
            THROW_HARD('cannot chose npix in automatic mode; exec_edge_detect')
        endif
        if( cline%defined('thres_low') .and. params%automatic .eq. 'yes') then
            THROW_HARD('cannot chose thres_low in automatic mode; exec_edge_detect')
        endif
        if( cline%defined('thres_up') .and. params%automatic .eq. 'yes') then
            THROW_HARD('cannot chose thres_up in automatic mode; exec_edge_detect')
        endif
        if( cline%defined('thres') .and. params%detector .eq. 'canny') then
            THROW_HARD('canny needs double thresholding; exec_edge_detect')
        endif
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        do iptcl=1,params%nptcls
            call build%img%new([params%box,params%box,1],params%smpd,wthreads=.false.)
            call build%img%read(params%stk, iptcl)
            call doit(build%img)
            call build%img%write(params%outstk, iptcl)
            call build%img%kill
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_EDGE_DETECT NORMAL STOP ****')

    contains

            subroutine doit( img )
                use simple_segmentation
                class(image), intent(inout) :: img
                type (image)      :: img_grad
                real, allocatable :: grad(:,:,:)
                real    :: thresh(1), ave, sdev, maxv, minv, lp(1)
                real    :: smpd
                integer :: ldim(3)
                thresh = 0. !initialise
                ldim = img%get_ldim()
                smpd = img%get_smpd()
                allocate(grad(ldim(1), ldim(2), ldim(3)), source = 0.)
                if(cline%defined('lp')) lp(1) = params%lp
                select case ( params%detector )
                case ('sobel')
                    if( cline%defined('thres') )then
                        thresh(1) = params%thres
                        call sobel(img,thresh)
                    else if( cline%defined('npix') )then
                        call automatic_thresh_sobel(img,real(params%npix)/(real(ldim(1)*ldim(2))))
                    elseif( params%automatic .eq. 'yes') then
                        call img%scale_pixels([0.,255.])
                        call img%calc_gradient(grad)
                        call img_grad%new(ldim, smpd)
                        call img_grad%set_rmat(grad,.false.)
                        call img_grad%stats( ave, sdev, maxv, minv )
                        call img_grad%kill
                        thresh(1) = ave + 0.7*sdev
                        write(logfhandle,*) 'Selected threshold: ', thresh
                        call sobel(img,thresh)
                        deallocate(grad)
                    else
                        THROW_HARD('If not automatic threshold needed; exec_edge_detect')
                    endif
                case('canny')
                    if(ldim(3) .ne. 1) THROW_HARD('Canny for vol is not implemented; exec_edge_detect')
                    if( params%automatic .eq. 'no' ) then
                        if(.not. cline%defined('thres_low') .or. .not. cline%defined('thres_up') )then
                            THROW_HARD('both upper and lower threshold needed; exec_edge_detect')
                        else
                            if( cline%defined('lp')) then
                                call canny(img,img,thresh=[params%thres_low, params%thres_up],lp=lp(1)) !inout/output image coincide
                            else
                                call canny(img,img,thresh=[params%thres_low, params%thres_up])
                            endif
                        endif
                    elseif( params%automatic .eq. 'yes') then
                        if(cline%defined('thres_low') .or. cline%defined('thres_up')) then
                            THROW_HARD('cannot define thresholds in automatic mode; exec_edge_detect')
                        else
                            if( .not. cline%defined('lp')) then
                              THROW_HARD('Canny automatic requires lp in input; exec_edge_detect')
                            else
                                call canny(img,lp = lp(1))
                            endif
                        endif
                    else
                      THROW_HARD('Wrong input for automatic parameter!; exec_edge_detect')
                    endif
                case DEFAULT
                    THROW_HARD('Unknown detector argument; exec_edge_detect')
               end select
            end subroutine

    end subroutine exec_edge_detect

    subroutine exec_filter( self, cline )
        use simple_procimgstk
        use simple_tvfilter, only: tvfilter
        class(commander_filter), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters)  :: params
        type(builder)     :: build
        type(tvfilter)    :: tvfilt
        real, allocatable :: fsc(:), optlp(:), res(:)
        real, parameter   :: SIGMA_DEFAULT=1.0
        real              :: fsc05, fsc0143
        integer           :: find
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',  'no')
        if( .not. cline%defined('outstk') ) call cline%set('outstk', 'filtered'//STK_EXT)
        if( .not. cline%defined('outvol') ) call cline%set('outvol', 'filtered.mrc')
        if( cline%defined('stk') )then
            ! 2D
            call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
            if( .not.file_exists(params%stk) ) THROW_HARD('cannot find input stack (stk)')
            if( cline%defined('filter') )then
                select case(trim(params%filter))
                    case('tv')
                        call tvfilter_imgfile(params%stk, params%outstk, params%smpd, params%lambda)
                    case('nlmean')
                        if( cline%defined('sigma') )then
                            call nlmean_imgfile(params%stk, params%outstk, params%smpd, params%sigma)
                        else
                            call nlmean_imgfile(params%stk, params%outstk, params%smpd)
                        endif
                    case('icm')
                        call icm_imgfile(params%stk, params%outstk, params%smpd, params%lambda)
                    case DEFAULT
                        THROW_HARD('Unknown filter!')
                end select
            else if( params%phrand .eq. 'no')then
                ! Band pass
                if( cline%defined('lp') .and. cline%defined('hp') )then
                    call bp_imgfile(params%stk, params%outstk, params%smpd, params%hp, params%lp)
                else if( cline%defined('lp') )then
                    call bp_imgfile(params%stk, params%outstk, params%smpd, 0., params%lp)
                else if( cline%defined('hp') )then
                    call bp_imgfile(params%stk, params%outstk, params%smpd, params%hp, 0.)
                ! real-space
                else if( cline%defined('real_filter') )then
                    if( .not. cline%defined('winsz') .and. trim(params%real_filter) .ne. 'NLmean') THROW_HARD('need winsz input for real-space filtering')
                    call real_filter_imgfile(params%stk, params%outstk, params%smpd, trim(params%real_filter), nint(params%winsz))
                else
                    THROW_HARD('Nothing to do!')
                endif
            else if ( params%phrand.eq.'yes' )then
                ! Phase randomization
                if( .not. cline%defined('lp') ) THROW_HARD('low-pass limit needed 4 phase randomization')
                call phase_rand_imgfile(params%stk, params%outstk, params%smpd, params%lp)
            endif
        else
            ! 3D
            call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
            if( .not.file_exists(params%vols(1)) ) THROW_HARD('Cannot find input volume (vol1)')
            call build%vol%read(params%vols(1))
            if( params%phrand.eq.'no')then
                if( cline%defined('bfac') )then
                    ! bfactor
                    call build%vol%apply_bfac(params%bfac)
                    ! Band pass
                else if( cline%defined('hp') .and. cline%defined('lp') )then
                    call build%vol%bp(params%hp, params%lp)
                else if( cline%defined('hp') )then
                    call build%vol%bp(params%hp, 0.)
                else if( cline%defined('lp') )then
                    call build%vol%bp(0., params%lp)
                else if( params%tophat .eq. 'yes' .and. cline%defined('find') )then
                    call build%vol%tophat(params%find)
                    ! real-space
                else if( cline%defined('real_filter') )then
                    if( .not. cline%defined('winsz') ) THROW_HARD('need winsz input for real-space filtering')
                    call build%vol%real_space_filter(nint(params%winsz), params%real_filter)
                else if( cline%defined('fsc') )then
                    if( .not.file_exists(params%fsc)) THROW_HARD('Cannot find FSC filter (vol_filt)')
                    ! resolution & optimal low-pass filter from FSC
                    fsc   = file2rarr(params%fsc)
                    optlp = fsc2optlp(fsc)
                    res   = build%vol%get_res()
                    call get_resolution( fsc, res, fsc05, fsc0143 )
                    where(res < TINY) optlp = 0.
                    call build%vol%apply_filter(optlp)
                else if( cline%defined('filter') )then
                    select case(trim(params%filter))
                        case('tv')
                            call tvfilt%new
                            call tvfilt%apply_filter_3d(build%vol, params%lambda)
                            call tvfilt%kill
                        case('nlmean')
                            call build%vol%NLmean3D()
                        case('icm')
                            call build%vol%ICM3D(params%lambda)
                         case('tent')
                            if( .not. cline%defined('winsz') ) THROW_HARD('need winsz input for tent filter')
                            call build%vol%bartlett_reg_3D(nint(params%winsz))
                        case DEFAULT
                            THROW_HARD('Unknown filter!')
                    end select
                else
                    THROW_HARD('Nothing to do!')
                endif
            else
                if( .not. cline%defined('lp') )THROW_HARD('low-pass limit needed 4 phase randomization')
                call build%vol%phase_rand(params%lp)
            endif
            if( params%outvol .ne. '' ) call build%vol%write(params%outvol, del_if_exists=.true.)
        endif
        ! cleanup
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_FILTER NORMAL STOP ****')
    end subroutine exec_filter

    subroutine exec_ppca_denoise( self, cline )
        use simple_imgproc,     only: make_pcavecs
        use simple_pca,         only: pca
        use simple_ppca,        only: ppca
        use simple_pca_svd,     only: pca_svd
        use simple_kpca_svd,    only: kpca_svd, suggest_kpca_nystrom_neigs
        class(commander_ppca_denoise), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        integer,           parameter   :: MAXPCAITS = 15
        class(pca),        pointer     :: pca_ptr         => null()
        type(ppca),        pointer     :: ppca_ptr_typed  => null()
        type(kpca_svd),    pointer     :: kpca_ptr        => null()
        type(image),       allocatable :: imgs(:)
        real,              allocatable :: avg(:), pcavecs(:, :), tmpvec(:), zavg(:), corrvec(:)
        type(simple_nice_comm)         :: nice_comm
        type(parameters)               :: params
        type(builder)                  :: build
        integer(int64)                 :: t0, t1
        real(real64)                   :: trate
        integer                        :: npix, iptcl, j, neigs
        logical                        :: l_transp_pca, l_profile_pca, l_hybrid_resid
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',  'no')
        if( .not. cline%defined('outstk') ) call cline%set('outstk', 'ppca_denoised'//STK_EXT)
        ! doesn't work if projfile given - may need to mod in future
        if(cline%defined('projfile')) call cline%delete('projfile')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        if( .not.file_exists(params%stk) ) THROW_HARD('cannot find input stack (stk)')
        l_hybrid_resid = trim(params%pca_mode) .eq. 'ppca_kpca_resid'
        l_profile_pca = trim(params%pca_mode) .eq. 'kpca' .or. trim(params%pca_mode) .eq. 'ppca' .or. l_hybrid_resid
        ! nice communicator init
        call nice_comm%init(params%niceprocid, params%niceserver)
        call nice_comm%cycle()
        call system_clock(t0, trate)
        allocate(imgs(params%nptcls))
        do iptcl = 1, params%nptcls
            call imgs(iptcl)%new([params%box,params%box,1], params%smpd)
            call imgs(iptcl)%read(params%stk, iptcl)
        end do
        call system_clock(t1)
        if( l_profile_pca )then
            if( trim(params%pca_mode) .eq. 'kpca' )then
                write(logfhandle,'(A,F8.3,A,I8)') 'kPCA denoise read stack: ', real(t1-t0)/real(trate), ' s; nptcls=', params%nptcls
            else if( l_hybrid_resid )then
                write(logfhandle,'(A,F8.3,A,I8)') 'PPCA+kPCA residual denoise read stack: ', real(t1-t0)/real(trate), ' s; nptcls=', params%nptcls
            else
                write(logfhandle,'(A,F8.3,A,I8)') 'PPCA denoise read stack: ', real(t1-t0)/real(trate), ' s; nptcls=', params%nptcls
            endif
        endif
        l_transp_pca = (trim(params%transp_pca) .eq. 'yes')
        call system_clock(t0)
        call make_pcavecs(imgs, npix, avg, pcavecs, transp=l_transp_pca)
        call system_clock(t1)
        if( l_profile_pca )then
            if( trim(params%pca_mode) .eq. 'kpca' )then
                write(logfhandle,'(A,F8.3,A,I8)') 'kPCA denoise make_pcavecs: ', real(t1-t0)/real(trate), ' s; npix=', npix
            else if( l_hybrid_resid )then
                write(logfhandle,'(A,F8.3,A,I8)') 'PPCA+kPCA residual denoise make_pcavecs: ', real(t1-t0)/real(trate), ' s; npix=', npix
            else
                write(logfhandle,'(A,F8.3,A,I8)') 'PPCA denoise make_pcavecs: ', real(t1-t0)/real(trate), ' s; npix=', npix
            endif
        endif
        neigs = params%neigs
        if( (trim(params%pca_mode) .eq. 'kpca' .or. l_hybrid_resid) .and. trim(params%kpca_backend) .eq. 'nystrom' .and. neigs <= 0 )then
            neigs = max(8, min(64, max(1, params%kpca_nystrom_npts / 2)))
            if( l_hybrid_resid )then
                write(logfhandle,'(A,I8,A)') 'PPCA+kPCA residual auto-selected neigs: ', neigs, ' (Nyström safe heuristic)'
            else
                write(logfhandle,'(A,I8,A)') 'kPCA denoise auto-selected neigs: ', neigs, ' (Nyström safe heuristic)'
            endif
            call flush(logfhandle)
        endif
        if( l_transp_pca )then
            neigs = min(max(neigs, 1), max(npix-1, 1))
        else
            neigs = min(max(neigs, 1), max(params%nptcls-1, 1))
        endif
        if( l_hybrid_resid )then
            ! Residual hybrid denoiser:
            ! 1. Fit PPCA to the centered particle stack and reconstruct the linear denoised estimate.
            ! 2. Form a residual stack r = x - x_ppca that contains what the linear model did not explain.
            ! 3. Run kPCA only on that residual stack, not on the full images.
            ! 4. Write x_hybrid = avg + x_ppca + alpha * r_kpca.
            !
            ! Why this exists as a separate mode:
            ! - A direct PPCA -> kPCA image-domain cascade tended to make images blurrier.
            ! - The intent here is different: PPCA handles the global linear denoising,
            !   while kPCA is restricted to modeling structured nonlinear leftovers.
            ! - Keeping this in its own pca_mode keeps the existing ppca and kpca
            !   behavior untouched and makes the hybrid easy to benchmark independently.
            !
            ! Current scope:
            ! - Implemented only for transp_pca=no in the image-denoise commander.
            ! - Uses the existing PPCA and kPCA solvers as black-box stages.
            ! - alpha controls how strongly the residual kPCA correction is blended back.
            if( l_transp_pca ) THROW_HARD('ppca_kpca_resid currently supports transp_pca=no only')
            allocate(ppca_ptr_typed, kpca_ptr)
            call ppca_ptr_typed%new(params%nptcls, npix, neigs)
            call kpca_ptr%new(params%nptcls, npix, neigs)
            call kpca_ptr%set_params(params%nthr, params%kpca_ker, params%kpca_backend,&
            &params%kpca_nystrom_npts, params%kpca_rbf_gamma, params%kpca_nystrom_local_nbrs, params%kpca_cosine_weight_power)
            allocate(zavg(npix), tmpvec(npix), corrvec(npix), source=0.)
            call system_clock(t0)
            call ppca_ptr_typed%master(pcavecs, MAXPCAITS)
            do iptcl = 1, params%nptcls
                call ppca_ptr_typed%generate(iptcl, zavg, tmpvec)
                pcavecs(:,iptcl) = pcavecs(:,iptcl) - tmpvec
            enddo
            if( l_profile_pca ) write(logfhandle,'(A)') 'PPCA+kPCA residual: residual stack formed from PPCA reconstruction'
            call flush(logfhandle)
            call kpca_ptr%master(pcavecs, MAXPCAITS)
            call system_clock(t1)
            if( l_profile_pca ) write(logfhandle,'(A,F8.3,A)') 'PPCA+kPCA residual denoise master: ', real(t1-t0)/real(trate), ' s'
            call system_clock(t0)
            block
                ! Hybrid residual diagnostics:
                ! - residual RMS stats: size of the residual left after PPCA alone.
                ! - correction RMS stats: size of the damped kPCA correction being added back.
                ! - remaining residual RMS stats: size of (residual - alpha * correction).
                !   If this drops meaningfully, the residual kPCA stage is explaining structure
                !   that PPCA did not capture.
                ! - residual/correction cosine stats: directional alignment between the PPCA
                !   residual and the kPCA correction. High positive values mean the correction
                !   is targeting the residual rather than acting like an unrelated smoother.
                !
                ! Practical interpretation:
                ! - correction RMS << residual RMS  => hybrid is barely changing the PPCA result
                ! - remaining residual RMS much lower than residual RMS => hybrid is active/useful
                ! - cosine near 1 => correction tracks the residual well
                ! - cosine near 0 or negative => correction is weakly aligned or potentially harmful
                real(dp) :: resid_rms, corr_rms, remain_rms, align_cos
                real(dp) :: resid_norm2, corr_norm2, dot_rc
                real(dp) :: resid_sum, resid_sumsq, resid_max
                real(dp) :: corr_sum, corr_sumsq, corr_max
                real(dp) :: remain_sum, remain_sumsq, remain_max
                real(dp) :: align_sum, align_sumsq, align_max
                real(dp) :: inv_np, alpha
                alpha = real(params%ppca_kpca_resid_alpha, dp)
                inv_np = 1._dp / real(npix, dp)
                resid_sum = 0._dp;  resid_sumsq = 0._dp;  resid_max = 0._dp
                corr_sum  = 0._dp;  corr_sumsq  = 0._dp;  corr_max  = 0._dp
                remain_sum = 0._dp; remain_sumsq = 0._dp; remain_max = 0._dp
                align_sum = 0._dp;  align_sumsq = 0._dp;  align_max = -huge(1._dp)
                do iptcl = 1, params%nptcls
                    call ppca_ptr_typed%generate(iptcl, zavg, tmpvec)
                    call kpca_ptr%generate(iptcl, zavg, corrvec)
                    resid_norm2 = sum(real(pcavecs(:,iptcl), dp)**2)
                    corr_norm2  = sum(real(corrvec, dp)**2)
                    dot_rc      = sum(real(pcavecs(:,iptcl), dp) * real(corrvec, dp))
                    resid_rms = sqrt(resid_norm2 * inv_np)
                    corr_rms  = sqrt((alpha * alpha) * corr_norm2 * inv_np)
                    remain_rms = sqrt(sum((real(pcavecs(:,iptcl), dp) - alpha * real(corrvec, dp))**2) * inv_np)
                    if( resid_norm2 > DTINY .and. corr_norm2 > DTINY )then
                        align_cos = dot_rc / sqrt(resid_norm2 * corr_norm2)
                    else
                        align_cos = 0._dp
                    endif
                    resid_sum   = resid_sum + resid_rms
                    resid_sumsq = resid_sumsq + resid_rms * resid_rms
                    resid_max   = max(resid_max, resid_rms)
                    corr_sum    = corr_sum + corr_rms
                    corr_sumsq  = corr_sumsq + corr_rms * corr_rms
                    corr_max    = max(corr_max, corr_rms)
                    remain_sum   = remain_sum + remain_rms
                    remain_sumsq = remain_sumsq + remain_rms * remain_rms
                    remain_max   = max(remain_max, remain_rms)
                    align_sum   = align_sum + align_cos
                    align_sumsq = align_sumsq + align_cos * align_cos
                    align_max   = max(align_max, align_cos)
                    tmpvec = avg + tmpvec + params%ppca_kpca_resid_alpha * corrvec
                    call imgs(iptcl)%unserialize(tmpvec)
                    call imgs(iptcl)%write(params%outstk, iptcl)
                    call imgs(iptcl)%kill
                end do
                if( l_profile_pca )then
                    write(logfhandle,'(A,3(A,F9.4))') 'PPCA+kPCA residual RMS stats:', &
                        ' mean=', resid_sum / real(params%nptcls, dp), &
                        ' std=', sqrt(max(0._dp, resid_sumsq / real(params%nptcls, dp) - (resid_sum / real(params%nptcls, dp))**2)), &
                        ' max=', resid_max
                    write(logfhandle,'(A,3(A,F9.4))') 'PPCA+kPCA correction RMS stats:', &
                        ' mean=', corr_sum / real(params%nptcls, dp), &
                        ' std=', sqrt(max(0._dp, corr_sumsq / real(params%nptcls, dp) - (corr_sum / real(params%nptcls, dp))**2)), &
                        ' max=', corr_max
                    write(logfhandle,'(A,3(A,F9.4))') 'PPCA+kPCA remaining residual RMS stats:', &
                        ' mean=', remain_sum / real(params%nptcls, dp), &
                        ' std=', sqrt(max(0._dp, remain_sumsq / real(params%nptcls, dp) - (remain_sum / real(params%nptcls, dp))**2)), &
                        ' max=', remain_max
                    write(logfhandle,'(A,3(A,F9.4))') 'PPCA+kPCA residual/correction cosine stats:', &
                        ' mean=', align_sum / real(params%nptcls, dp), &
                        ' std=', sqrt(max(0._dp, align_sumsq / real(params%nptcls, dp) - (align_sum / real(params%nptcls, dp))**2)), &
                        ' max=', align_max
                endif
            end block
            call system_clock(t1)
            if( l_profile_pca ) write(logfhandle,'(A,F8.3,A,F6.3)') 'PPCA+kPCA residual denoise reconstruct/write: ', real(t1-t0)/real(trate), ' s; alpha=', params%ppca_kpca_resid_alpha
            deallocate(imgs)
            deallocate(zavg, tmpvec, corrvec)
            deallocate(ppca_ptr_typed, kpca_ptr)
            call build%kill_general_tbox
            call nice_comm%terminate()
            call simple_end('**** SIMPLE_PPCA_DENOISE NORMAL STOP ****')
            return
        endif
        ! pca allocation
        select case(trim(params%pca_mode))
            case('ppca')
                allocate(ppca :: pca_ptr)
            case('pca_svd')
                allocate(pca_svd     :: pca_ptr)
            case('kpca')
                allocate(kpca_svd    :: pca_ptr)
            case DEFAULT
                THROW_HARD('pca_mode must be ppca, ppca_kpca_resid, pca_svd, or kpca')
        end select
        if( l_transp_pca )then
            call pca_ptr%new(npix, params%nptcls, neigs)
            select type(pca_ptr)
                type is(kpca_svd)
                    call pca_ptr%set_params(params%nthr, params%kpca_ker,&
                    &params%kpca_backend, params%kpca_nystrom_npts, params%kpca_rbf_gamma, params%kpca_nystrom_local_nbrs, params%kpca_cosine_weight_power)
            end select
            call system_clock(t0)
            select type(pca_ptr)
                type is(kpca_svd)
                    call pca_ptr%master(pcavecs, MAXPCAITS)
                class default
                    call pca_ptr%master(pcavecs, MAXPCAITS)
            end select
            call system_clock(t1)
            if( l_profile_pca )then
                if( trim(params%pca_mode) .eq. 'kpca' )then
                    write(logfhandle,'(A,F8.3,A)') 'kPCA denoise master: ', real(t1-t0)/real(trate), ' s'
                else
                    write(logfhandle,'(A,F8.3,A)') 'PPCA denoise master: ', real(t1-t0)/real(trate), ' s'
                endif
            endif
            allocate(tmpvec(params%nptcls))
            call system_clock(t0)
            !$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
            do j = 1, npix
                call pca_ptr%generate(j, avg, tmpvec)
                pcavecs(:,j) = tmpvec
            end do
            !$omp end parallel do
            pcavecs = transpose(pcavecs)
            do iptcl = 1, params%nptcls
                call imgs(iptcl)%unserialize(pcavecs(:,iptcl))
                call imgs(iptcl)%write(params%outstk, iptcl)
                call imgs(iptcl)%kill
            end do
            call system_clock(t1)
            if( l_profile_pca )then
                if( trim(params%pca_mode) .eq. 'kpca' )then
                    write(logfhandle,'(A,F8.3,A)') 'kPCA denoise reconstruct/write: ', real(t1-t0)/real(trate), ' s'
                else
                    write(logfhandle,'(A,F8.3,A)') 'PPCA denoise reconstruct/write: ', real(t1-t0)/real(trate), ' s'
                endif
            endif
        else
            call pca_ptr%new(params%nptcls, npix, neigs)
            select type(pca_ptr)
                type is(kpca_svd)
                    call pca_ptr%set_params(params%nthr, params%kpca_ker, params%kpca_backend,&
                    &params%kpca_nystrom_npts, params%kpca_rbf_gamma, params%kpca_nystrom_local_nbrs, params%kpca_cosine_weight_power)
            end select
            call system_clock(t0)
            select type(pca_ptr)
                type is(kpca_svd)
                    call pca_ptr%master(pcavecs, MAXPCAITS)
                class default
                    call pca_ptr%master(pcavecs, MAXPCAITS)
            end select
            call system_clock(t1)
            if( l_profile_pca )then
                if( trim(params%pca_mode) .eq. 'kpca' )then
                    write(logfhandle,'(A,F8.3,A)') 'kPCA denoise master: ', real(t1-t0)/real(trate), ' s'
                else
                    write(logfhandle,'(A,F8.3,A)') 'PPCA denoise master: ', real(t1-t0)/real(trate), ' s'
                endif
            endif
            call system_clock(t0)
            !$omp parallel do private(iptcl) default(shared) proc_bind(close) schedule(static)
            do iptcl = 1, params%nptcls
                call pca_ptr%generate(iptcl, avg, pcavecs(:,iptcl))
            end do
            !$omp end parallel do
            do iptcl = 1, params%nptcls
                call imgs(iptcl)%unserialize(pcavecs(:,iptcl))
                call imgs(iptcl)%write(params%outstk, iptcl)
                call imgs(iptcl)%kill
            end do
            call system_clock(t1)
            if( l_profile_pca )then
                if( trim(params%pca_mode) .eq. 'kpca' )then
                    write(logfhandle,'(A,F8.3,A)') 'kPCA denoise reconstruct/write: ', real(t1-t0)/real(trate), ' s'
                else
                    write(logfhandle,'(A,F8.3,A)') 'PPCA denoise reconstruct/write: ', real(t1-t0)/real(trate), ' s'
                endif
            endif
        endif
        ! cleanup
        deallocate(imgs)
        call build%kill_general_tbox
        ! end gracefully
        call nice_comm%terminate()
        call simple_end('**** SIMPLE_PPCA_DENOISE NORMAL STOP ****')
    end subroutine exec_ppca_denoise

    subroutine exec_ppca_denoise_polarft_lines( self, cline )
        use simple_complex_ppca,               only: complex_ppca
        use simple_matcher_pftc_prep,          only: prep_pftc4align2D
        use simple_matcher_ptcl_batch,         only: prep_batch_particles2D, clean_batch_particles2D
        ! use simple_polarft_lines_ppca_stream,  only: stream_pft_lines_ppca, denoise_write_pft_lines_ppca
        class(commander_ppca_denoise_polarft_lines), intent(inout) :: self
        class(cmdline),                              intent(inout) :: cline
        type(parameters)               :: params
        type(builder)                  :: build
        type(complex_ppca)             :: model
        type(image), allocatable       :: ptcl_imgs(:), ptcl_match_imgs(:), ptcl_match_imgs_pad(:)
        real(dp),    allocatable       :: eigvals(:)
        integer                        :: batchsz_max, qfit, which_iter
        integer                        :: nptcls_total
        integer(kind=8)                :: nlines_fit, nlines_den
        if( .not. cline%defined('outfile') ) call cline%set('outfile', 'ppca_denoised_polarfts'//BIN_EXT)
        if( .not. cline%defined('neigs') )   call cline%set('neigs', 16.0)
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        nptcls_total = params%top - params%fromp + 1
        batchsz_max = min(nptcls_total, params%nthr * BATCHTHRSZ)
        qfit = max(1, params%neigs)
        which_iter = 1
        if( cline%defined('which_iter') ) which_iter = max(1, params%which_iter)
        call prep_batch_particles2D(params, build, batchsz_max, ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad)
        call prep_pftc4align2D(params, build, ptcl_match_imgs_pad, batchsz_max, which_iter, .false.)
        call clean_batch_particles2D(build, ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad)
        ! call stream_pft_lines_ppca(params, build, qfit, model, nlines_fit)
        eigvals = model%get_eigvals()
        write(logfhandle,'(A,I12)')      'PPCA polarft line fit count: ', nlines_fit
        write(logfhandle,'(A,F12.5)')    'PPCA polarft sigma2: ', model%get_sigma2()
        write(logfhandle,'(A,10(1X,ES12.5))') 'PPCA polarft leading eigvals:', eigvals(1:min(10,size(eigvals)))
        deallocate(eigvals)
        ! call denoise_write_pft_lines_ppca(params, build, model, params%outfile, nlines_den)
        write(logfhandle,'(A,I12)')   'PPCA polarft denoised line count: ', nlines_den
        write(logfhandle,'(A,A)')     'PPCA polarft output: ', params%outfile%to_char()
        call model%kill()
        call build%pftc%kill()
        call build%esig%kill()
        call build%kill_general_tbox()
        call simple_end('**** SIMPLE_PPCA_DENOISE_POLARFT_LINES NORMAL STOP ****')
    end subroutine exec_ppca_denoise_polarft_lines

    !> normalize is a program for normalization of MRC or SPIDER stacks and volumes.
    !! If you want to normalize your images inputted with stk, set norm=yes.
    !! hfun (e.g. hfun=sigm) controls the normalization function. If you want to
    !! perform noise normalization of the images set noise_norm=yes given a mask
    !! radius msk (pixels).
    subroutine exec_normalize( self, cline )
        use simple_procimgstk, only: norm_imgfile, noise_norm_imgfile
        class(commander_normalize), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters)  :: params
        type(builder)     :: build
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('noise_norm') ) call cline%set('noise_norm', 'no')
        if( cline%defined('stk')  .and. cline%defined('vol1') )THROW_HARD('Cannot operate on images AND volume at once')
        if( cline%defined('stk') )then
            ! 2D
            call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
            if( params%norm.eq.'yes' )then
                ! Normalization
                call norm_imgfile(params%stk, params%outstk, params%smpd)
            else if( params%noise_norm.eq.'yes' )then
                ! Noise normalization
                if( cline%defined('mskdiam') )then
                    call noise_norm_imgfile(params%stk, params%msk, params%outstk, params%smpd)
                else
                    THROW_HARD('need msk parameter for noise normalization')
                endif
            endif
        else if( cline%defined('vol1') )then
            ! 3D
            call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
            if( .not.file_exists(params%vols(1)) )THROW_HARD('Cannot find input volume')
            call build%vol%read(params%vols(1))
            if( params%norm.eq.'yes' )then
                call build%vol%norm()
                call build%vol%write(params%outvol, del_if_exists=.true.)
            else
                THROW_HARD('Normalization type not implemented yet')
            endif
        else
            THROW_HARD('No input images(s) or volume provided')
        endif
        ! cleanup
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_NORMALIZE NORMAL STOP ****')
    end subroutine exec_normalize

    !> provides re-scaling and clipping routines for MRC or SPIDER stacks and volumes
    subroutine exec_scale( self, cline )
        use simple_procimgstk, only: scale_and_clip_imgfile, scale_imgfile, pad_imgfile, clip_imgfile
        class(commander_scale), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(string), allocatable :: filenames(:)
        type(parameters) :: params
        type(builder)    :: build
        type(image)      :: vol2, img, img2
        type(stack_io)   :: stkio_r, stkio_w
        real             :: ave, sdev, var, med, smpd_new, scale, smpd_sc
        integer          :: ldim(3), ldim_scaled(3), nfiles, nframes, iframe, ifile
        integer          :: istk, nstks, ptcl_fromp, ptcl_top
        type(string)     :: fname, stkin, stkout, ext
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( cline%defined('stk') .and. cline%defined('vol1') ) THROW_HARD('Cannot operate on images AND volume at once')
        if( cline%defined('projfile') .and. cline%defined('fromp') .and. cline%defined('top')&
            &.and..not.cline%defined('stk') .and. .not.cline%defined('vol1') .and. .not.cline%defined('filetab'))then
            ! input comes from distributed execution scale_project
            call cline%set('oritype', 'stk')
            call build%init_params_and_build_spproj(cline, params)
            call build%spproj%read_segment('stk',params%projfile)
            nstks = build%spproj%os_stk%get_noris()
            if( nstks == 0 ) THROW_HARD('NO EXISTING STACK IN PROJECT')
            if( cline%defined('newbox') )then
                ldim_scaled  = [params%newbox,params%newbox,1]
                params%scale = real(params%newbox) / real(params%box)
            else if( cline%defined('scale') )then
                ldim_scaled(1) = round2even(real(params%box)*params%scale)
                ldim_scaled(2) = round2even(real(params%box)*params%scale)
                ldim_scaled(3) = 1
            else
                THROW_HARD('MISSING NEW DIMENSION!')
            endif
            do istk = params%fromp, params%top
                call progress(istk, nstks)
                stkin      = build%spproj%get_stkname(istk)
                ptcl_fromp = build%spproj%os_stk%get_fromp(istk)
                ptcl_top   = build%spproj%os_stk%get_top(istk)
                ext        = string('.')//fname2ext(stkin)
                if( cline%defined('dir_target') )then
                    stkout = filepath(params%dir_target,add2fbody(basename(stkin),ext,SCALE_SUFFIX))
                else
                    stkout = add2fbody(stkin,ext,SCALE_SUFFIX)
                endif
                call scale_imgfile(stkin, stkout, params%smpd, ldim_scaled, smpd_new)
            enddo
        else
            if( cline%defined('stk') )then
                ! 2D
                call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
                if( cline%defined('scale') .or. cline%defined('newbox') )then
                    call del_file(params%outstk)
                    ! Rescaling
                    ldim_scaled = [params%newbox,params%newbox,1] ! dimension of scaled
                    if( cline%defined('clip') )then
                        call scale_and_clip_imgfile(params%stk,params%outstk,params%smpd,ldim_scaled,&
                        [params%clip,params%clip,1],smpd_new)
                    else
                        call scale_imgfile(params%stk,params%outstk,params%smpd,ldim_scaled,smpd_new)
                        write(logfhandle,'(a,1x,i5)') 'BOX SIZE AFTER SCALING:', ldim_scaled(1)
                    endif
                    write(logfhandle,'(a,1x,f9.4)') 'SAMPLING DISTANCE AFTER SCALING:', smpd_new
                else if( cline%defined('clip') )then
                    call del_file(params%outstk)
                    ! Clipping
                    if( params%clip > params%box )then
                        call pad_imgfile(params%stk,params%outstk,[params%clip,params%clip,1],params%smpd)
                    else
                        call clip_imgfile(params%stk,params%outstk,[params%clip,params%clip,1],params%smpd)
                    endif
                endif
            else if( cline%defined('vol1') )then
                ! 3D
                call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
                if( .not.file_exists(params%vols(1)) ) THROW_HARD('Cannot find input volume')
                call build%vol%read(params%vols(1))
                if( cline%defined('scale') .or. cline%defined('newbox') )then
                    ! Rescaling
                    call vol2%new([params%newbox,params%newbox,params%newbox],params%smpd)
                    call build%vol%fft()
                    call vol2%set_ft(.true.)
                    if( params%newbox < params%box )then
                        call build%vol%clip(vol2)
                    else if( params%newbox > params%box )then
                        call build%vol%pad(vol2)
                    else
                        call vol2%copy(build%vol)
                    endif
                    call build%vol%copy(vol2)
                    call build%vol%ifft()
                    scale = real(params%newbox)/real(params%box)
                    params%box = params%newbox
                    smpd_new = params%smpd/scale
                    write(logfhandle,'(a,1x,f9.4)') 'SAMPLING DISTANCE AFTER SCALING:', smpd_new
                endif
                if( cline%defined('clip') )then
                    ! Clipping
                    call vol2%new([params%clip,params%clip,params%clip],params%smpd)
                    if( params%clip < params%box )then
                        call build%vol%clip(vol2)
                    else
                        call build%vol%stats( 'background', ave, sdev, var, med )
                        call build%vol%pad(vol2, backgr=med)
                    endif
                    call build%vol%copy(vol2)
                else
                     write(logfhandle,'(a,1x,i5)') 'BOX SIZE AFTER SCALING:', params%newbox
                endif
                if( params%outvol .ne. '' )call build%vol%write(params%outvol, del_if_exists=.true.)
            else if( cline%defined('filetab') )then
                call params%new(cline)
                call read_filetable(params%filetab, filenames)
                nfiles = size(filenames)
                call find_ldim_nptcls(filenames(1),ldim,nframes)
                ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
                if( cline%defined('newbox') )then
                    ldim_scaled = [params%newbox,params%newbox,1] ! dimension of scaled
                    params%scale = real(params%newbox) / real(ldim(1))
                else if( cline%defined('scale') )then
                    ldim_scaled(1) = round2even(real(ldim(1))*params%scale)
                    ldim_scaled(2) = round2even(real(ldim(2))*params%scale)
                    ldim_scaled(3)   = 1
                else
                    THROW_HARD('filetab key only in combination with scale or newbox!')
                endif
                call img%new(ldim,params%smpd)
                smpd_sc = params%smpd/params%scale
                call img2%new(ldim_scaled, smpd_sc)
                do ifile=1,nfiles
                    call progress(ifile, nfiles)
                    if( cline%defined('dir_target') )then
                        fname = filepath(params%dir_target,add2fbody(basename(filenames(ifile)), params%ext, SCALE_SUFFIX))
                    else
                        fname = add2fbody(filenames(ifile), params%ext, SCALE_SUFFIX)
                    endif
                    call stkio_w%open(fname, smpd_sc, 'write', box=ldim_scaled(1))
                    call find_ldim_nptcls(filenames(ifile),ldim,nframes)
                    ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
                    do iframe= 1, nframes
                        if( .not. stkio_r%stk_is_open() )then
                            call stkio_r%open(filenames(ifile), params%smpd, 'read')
                        else if( .not. stkio_r%same_stk(filenames(ifile), ldim) )then
                            call stkio_r%close
                            call stkio_r%open(filenames(ifile), params%smpd, 'read')
                        endif
                        call stkio_r%read(iframe, img)
                        call img%fft()
                        if( ldim_scaled(1) <= ldim(1) .and. ldim_scaled(2) <= ldim(2) .and. ldim_scaled(3) <= ldim(3) )then
                            call img%clip(img2)
                        else
                            call img%pad(img2)
                        endif
                        call img2%ifft()
                        call stkio_w%write(iframe, img2)
                    end do
                    call fname%kill
                    call stkio_w%close
                    call stkio_r%close
                end do
            else
                THROW_HARD('SIMPLE_SCALE needs input image(s) or volume or filetable!')
            endif
        endif
        ! cleanup
        call vol2%kill
        call img%kill
        call img2%kill
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_SCALE NORMAL STOP ****', print_simple=.false.)
        call qsys_job_finished(params, string('simple_commanders_imgops :: exec_scale'))
    end subroutine exec_scale

end module simple_commanders_imgops
