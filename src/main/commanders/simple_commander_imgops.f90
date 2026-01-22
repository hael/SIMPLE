!@descr: commanders for standard image operations: binarize, filter, denoise, normalize, scale etc.
module simple_commanders_imgops
use simple_commander_module_api
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
        use simple_imgproc,    only: make_pcavecs
        use simple_pca,        only: pca
        use simple_ppca_inmem, only: ppca_inmem
        use simple_pca_svd,    only: pca_svd
        use simple_kpca_svd,   only: kpca_svd
        class(commander_ppca_denoise), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        integer,           parameter   :: MAXPCAITS = 15
        class(pca),        pointer     :: pca_ptr  => null()
        type(image),       allocatable :: imgs(:)
        real,              allocatable :: avg(:), gen(:), pcavecs(:,:), tmpvec(:)
        type(simple_nice_communicator) :: nice_communicator
        type(parameters)               :: params
        type(builder)                  :: build
        integer                        :: npix, iptcl, j
        logical                        :: l_transp_pca
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',  'no')
        if( .not. cline%defined('outstk') ) call cline%set('outstk', 'ppca_denoised'//STK_EXT)
        ! doesn't work if projfile given - may need to mod in future
        if(cline%defined('projfile')) call cline%delete('projfile')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        if( .not.file_exists(params%stk) ) THROW_HARD('cannot find input stack (stk)')
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        call nice_communicator%cycle()
        allocate(imgs(params%nptcls))
        do iptcl = 1, params%nptcls
            call imgs(iptcl)%new([params%box,params%box,1], params%smpd)
            call imgs(iptcl)%read(params%stk, iptcl)
        end do
        l_transp_pca = (trim(params%transp_pca) .eq. 'yes')
        call make_pcavecs(imgs, npix, avg, pcavecs, transp=l_transp_pca)
        ! pca allocation
        select case(trim(params%pca_mode))
            case('ppca')
                allocate(ppca_inmem :: pca_ptr)
            case('pca_svd')
                allocate(pca_svd    :: pca_ptr)
            case('kpca')
                allocate(kpca_svd   :: pca_ptr)
        end select
        if( l_transp_pca )then
            call pca_ptr%new(npix, params%nptcls, params%neigs)
            call pca_ptr%master(pcavecs, MAXPCAITS)
            allocate(tmpvec(params%nptcls))
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
        else
            call pca_ptr%new(params%nptcls, npix, params%neigs)
            call pca_ptr%master(pcavecs, MAXPCAITS)
            allocate(gen(npix))
            do iptcl = 1, params%nptcls
                call pca_ptr%generate(iptcl, avg, gen)
                call imgs(iptcl)%unserialize(gen)
                call imgs(iptcl)%write(params%outstk, iptcl)
                call imgs(iptcl)%kill
            end do
        endif
        ! cleanup
        deallocate(imgs)
        call build%kill_general_tbox
        ! end gracefully
        call nice_communicator%terminate()
        call simple_end('**** SIMPLE_PPCA_DENOISE NORMAL STOP ****')
    end subroutine exec_ppca_denoise

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
        call qsys_job_finished(string('simple_commanders_imgops :: exec_scale'))
    end subroutine exec_scale

end module simple_commanders_imgops
