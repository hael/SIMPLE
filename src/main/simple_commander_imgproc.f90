! concrete commander: general image processing routines
module simple_commander_imgproc
include 'simple_lib.f08'
!use simple_binoris_io      ! use all in there
use simple_procimgfile     ! use all in there
use simple_image,          only: image
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base

implicit none

public :: binarise_commander
public :: convert_commander
public :: corrcompare_commander
public :: ctfops_commander
public :: filter_commander
public :: image_diff_commander
public :: image_smat_commander
public :: normalize_commander
public :: scale_commander
public :: stack_commander
public :: stackops_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: binarise_commander
  contains
    procedure :: execute      => exec_binarise
end type binarise_commander
type, extends(commander_base) :: convert_commander
  contains
    procedure :: execute      => exec_convert
end type convert_commander
type, extends(commander_base) :: corrcompare_commander
  contains
    procedure :: execute      => exec_corrcompare
end type corrcompare_commander
type, extends(commander_base) :: ctfops_commander
  contains
    procedure :: execute      => exec_ctfops
end type ctfops_commander
type, extends(commander_base) :: filter_commander
  contains
    procedure :: execute      => exec_filter
end type filter_commander
type, extends(commander_base) :: image_smat_commander
 contains
   procedure :: execute      => exec_image_smat
end type image_smat_commander
type, extends(commander_base) :: image_diff_commander
 contains
   procedure :: execute      => exec_image_diff
end type image_diff_commander
type, extends(commander_base) :: normalize_commander
  contains
    procedure :: execute      => exec_normalize
end type normalize_commander
type, extends(commander_base) :: scale_commander
  contains
    procedure :: execute      => exec_scale
end type scale_commander
type, extends(commander_base) :: stack_commander
  contains
    procedure :: execute      => exec_stack
end type stack_commander
type, extends(commander_base) :: stackops_commander
  contains
    procedure :: execute      => exec_stackops
end type stackops_commander

contains

    !> for binarisation of stacks and volumes
    subroutine exec_binarise( self, cline )
        class(binarise_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer      :: igrow, iptcl
        ! error check
        if( .not. cline%defined('stk') .and. .not. cline%defined('vol1') )then
            stop 'ERROR! stk or vol1 needs to be present; simple_binarise'
        endif
        if( cline%defined('stk') .and. cline%defined('vol1') )then
            stop 'ERROR! either stk or vol1 key can be present, not both; simple_binarise'
        endif
        if( cline%defined('thres') .and. cline%defined('npix') )then
            stop 'ERROR! either thres-based or npix-based binarisation; both keys cannot be present; simple_binarise'
        endif
        p = params(cline)                                     ! parameters generated
        if( cline%defined('stk') )then
            call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
            do iptcl=1,p%nptcls
                call b%img%read(p%stk, iptcl)
                call doit(b%img)
                call b%img%write(p%outstk, iptcl)
            end do
        else if( cline%defined('vol1') )then
            call b%build_general_tbox(p, cline)               ! general objects built
            call b%vol%read(p%vols(1)%str)
            call doit(b%vol)
            call b%vol%write(p%outvol)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_BINARISE NORMAL STOP ****')

        contains

            subroutine doit( img_or_vol )
                class(image), intent(inout) :: img_or_vol
                if( cline%defined('thres') )then
                    call img_or_vol%bin(p%thres)
                else if( cline%defined('npix') )then
                    call img_or_vol%bin(p%npix)
                else
                   call img_or_vol%bin_kmeans
                endif
                write(*,'(a,1x,i9)') '# FOREGROUND PIXELS:', img_or_vol%nforeground()
                write(*,'(a,1x,i9)') '# BACKGROUND PIXELS:', img_or_vol%nbackground()
                if( cline%defined('grow') )then
                    do igrow=1,p%grow
                        call img_or_vol%grow_bin
                    end do
                endif
                if( cline%defined('edge') ) call img_or_vol%cos_edge(p%edge)
                if( cline%defined('neg')  ) call img_or_vol%bin_inv
            end subroutine

    end subroutine exec_binarise

    !> for converting between SPIDER and MRC formats
    subroutine exec_convert( self, cline )
        class(convert_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params), target :: p
        type(build),  target :: b
        integer              :: iptcl
        p = params(cline, allow_mix=.true.) ! parameters generated
        if( cline%defined('stk') )then
            call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
            do iptcl=1,p%nptcls
                call progress(iptcl, p%nptcls)
                call b%img%read(p%stk, iptcl)
                call b%img%write(p%outstk, iptcl)
            end do
        else if( cline%defined('vol1') )then
            call b%build_general_tbox(p, cline)               ! general objects built
            call b%vol%read(p%vols(1)%str)
            call b%img%write(p%outvol)
        else
            stop 'either vol1 or stk argument required to execute simple_convert'
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CONVERT NORMAL STOP ****')
    end subroutine exec_convert

    subroutine exec_corrcompare( self, cline )
        class(corrcompare_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)      :: p
        type(build)       :: b
        integer           :: iptcl, j
        real              :: corr, ave, sdev, var, fsc05, fsc0143
        real, allocatable :: res(:), corrs(:), corrs_sum(:)
        logical           :: err
        p = params(cline)                                 ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        res = b%img%get_res()
        allocate(corrs(b%img%get_filtsz()))
        if( cline%defined('msk') )then
            if( p%stats .eq. 'yes' )then
                deallocate(corrs)
                allocate(corrs(p%nptcls), stat=alloc_stat)
                if(alloc_stat.ne.0)call allocchk('In: simple_commander_imgproc::corrcompare 1 ',alloc_stat)
            endif
            do iptcl=1,p%nptcls
                call b%img%read(p%stk, iptcl)
                call b%img_copy%read(p%stk2, iptcl)
                if( cline%defined('lp') )then
                    call b%img%mask(p%msk, 'soft')
                    call b%img_copy%mask(p%msk, 'soft')
                    corr = b%img%corr(b%img_copy, lp_dyn=p%lp)
                    if( p%stats .eq. 'yes' )then
                        corrs(iptcl) = corr
                    else
                        write(*,'(A,1X,F7.3)') '>>> FOURIER CORRELATION:', corr
                    endif
                else
                    call b%img%mask(p%msk, 'hard')
                    call b%img_copy%mask(p%msk, 'hard')
                    corr = b%img%real_corr(b%img_copy)
                    if( p%stats .eq. 'yes' )then
                        corrs(iptcl) = corr
                    else
                        write(*,'(A,1X,F7.3)') '>>> REAL-SPACE CORRELATION:', corr
                    endif
                endif
                call progress(iptcl,p%nptcls)
            end do
            if( p%stats .eq. 'yes' )then
                call moment(corrs, ave, sdev, var, err )
                if( cline%defined('lp') )then
                    write(*,'(A,1X,F7.3)') '>>> FOURIER CORRELATION AVERAGE:', ave
                    write(*,'(A,1X,F7.3)') '>>> FOURIER CORRELATION STDEV:  ', sdev
                else
                    write(*,'(A,1X,F7.3)') '>>> REAL-SPACE CORRELATION AVERAGE:', ave
                    write(*,'(A,1X,F7.3)') '>>> REAL-SPACE CORRELATION STDEV:  ', sdev
                endif
                if( allocated(corrs) ) deallocate(corrs)
            endif
        else
            if( .not. cline%defined('smpd') ) stop 'need smpd 4 FRC comarison!'
            do iptcl=1,p%nptcls
                call b%img%read(p%stk, iptcl)
                call b%img_copy%read(p%stk2, iptcl)
                call b%img%fsc(b%img_copy, corrs)
                if( .not. allocated(corrs_sum) )then
                    allocate(corrs_sum(size(corrs)), stat=alloc_stat)
                    if(alloc_stat.ne.0)call allocchk('In: simple_commander_imgproc:: corrcompare , 2',alloc_stat)
                    corrs_sum = 0.
                endif
                corrs_sum = corrs_sum+corrs
            end do
            corrs_sum = corrs_sum/real(p%nptcls)
            do j=1,size(res)
                write(*,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs_sum(j)
            end do
            call get_resolution( corrs_sum, res, fsc05, fsc0143 )
            write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', fsc0143
            write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', fsc05
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CORRCOMPARE NORMAL STOP ****')
    end subroutine exec_corrcompare

    !> for applying CTF to stacked images
    subroutine exec_ctfops( self, cline )
        use simple_ctf,         only: ctf
        class(ctfops_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        type(ctf)    :: tfun
        real         :: dfx, dfy, angast
        p = params(cline)                                 ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        if( cline%defined('oritab') .or. cline%defined('deftab') )then
            call b%raise_hard_ctf_exception(p)
        else
            stop 'oritab/deftab with CTF info needed for phase flipping/multiplication/CTF image generation'
        endif
        if( .not. cline%defined('stk') )then
            dfx = b%a%get(1,'dfx')
            if( b%a%isthere('dfy') )then
                dfy    = b%a%get(1,'dfy')
                angast = b%a%get(1,'angast')
            else
                dfy = dfx
                angast = 0.
            endif
            tfun = ctf(p%smpd, b%a%get(1,'kv'), b%a%get(1,'cs'), b%a%get(1,'fraca'))
            call tfun%ctf2img(b%img, dfx, 'ctf', dfy, angast)
            call b%img%ft2img('real', b%img_copy)
            call b%img_copy%write(p%outstk, 1)
            return
        endif
        if( p%ctf .ne. 'no' )then
            select case( p%ctf )
                case( 'flip' )
                    if( p%neg .eq. 'yes' )then
                        call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, 'flipneg')
                    else
                        call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, 'flip')
                    endif
                case( 'yes' )
                    if( p%neg .eq. 'yes' )then
                        call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, 'neg')
                    else
                        call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, 'ctf')
                    endif
                case DEFAULT
                    stop 'Unknown ctf argument'
            end select
        else
            stop 'Nothing to do!'
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CTFOPS NORMAL STOP ****')
    end subroutine exec_ctfops

    subroutine exec_filter( self, cline )
        class(filter_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        real :: width
        p     = params(cline) ! parameters generated
        width = 10.
        if( cline%defined('width') ) width = p%width
        if( cline%defined('stk') )then
            ! 2D
            call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
            if( .not.file_exists(p%stk) )stop 'Cannot find input stack (stk)'
            if( p%phrand .eq. 'no')then
                ! Band pass
                if( cline%defined('lp') .and. cline%defined('hp') )then
                    call bp_imgfile(p%stk, p%outstk, p%smpd, p%hp, p%lp, width=width)
                else if( cline%defined('lp') )then
                    call bp_imgfile(p%stk, p%outstk, p%smpd, 0., p%lp, width=width)
                else if( cline%defined('hp') )then
                    call bp_imgfile(p%stk, p%outstk, p%smpd, p%hp, 0., width=width)
                ! real-space
                else if( cline%defined('real_filter') )then
                    if( .not. cline%defined('winsz') )&
                    call simple_stop('need winsz input for real-space filtering; commander_imgproc :: exec_filter')
                    call real_filter_imgfile(p%stk, p%outstk, p%smpd, trim(p%real_filter), nint(p%winsz))
                else
                    stop 'Nothing to do!'
                endif
            else if ( p%phrand.eq.'yes' )then
                ! Phase randomization
                if( .not. cline%defined('lp') )stop 'low-pass limit needed 4 phase randomization'
                call phase_rand_imgfile(p%stk, p%outstk, p%smpd, p%lp)
            endif
        else
            ! 3D
            call b%build_general_tbox(p, cline) ! general objects built
            if( .not.file_exists(p%vols(1)%str) )stop 'Cannot find input volume (vol1)'
            call b%vol%read(p%vols(1)%str)
            if( p%phrand.eq.'no')then
                if( cline%defined('bfac') )then
                ! bfactor
                    call b%vol%apply_bfac(p%bfac)
                ! Band pass
                else if( cline%defined('hp') .and. cline%defined('lp') )then
                    call b%vol%bp(p%hp, p%lp, width=width)
                else if( cline%defined('hp') )then
                    call b%vol%bp(p%hp, 0., width=width)
                else if( cline%defined('lp') )then
                    call b%vol%bp(0., p%lp, width=width)
                ! real-space
                else if( cline%defined('real_filter') )then
                    if( .not. cline%defined('winsz') )&
                        call simple_stop('need winsz input for real-space filtering; commander_imgproc :: exec_filter')
                    call b%vol%real_space_filter(nint(p%winsz), p%real_filter)
                else
                    stop 'Nothing to do!'
                endif
            else
                if( .not. cline%defined('lp') )stop 'low-pass limit needed 4 phase randomization'
                call b%vol%phase_rand(p%lp)
            endif
            if( p%outvol .ne. '' )call b%vol%write(p%outvol, del_if_exists=.true.)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_FILTER NORMAL STOP ****')
    end subroutine exec_filter

    !> for creating a similarity matrix based on image2image correlation
    subroutine exec_image_smat(self, cline)
        use simple_corrmat, only: calc_cartesian_corrmat
        class(image_smat_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params)         :: p
        type(build)          :: b
        integer              :: iptcl, funit, io_stat
        real, allocatable    :: corrmat(:,:)
        p = params(cline, .false.) ! constants & derived constants produced
        call b%build_general_tbox(p, cline, do3d=.false., nooritab=.true.) ! general objects built (no oritab reading)
        allocate(b%imgs_sym(p%nptcls), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('In: simple_image_smat, 1',alloc_stat)
        do iptcl=1,p%nptcls
            call b%imgs_sym(iptcl)%new([p%box,p%box,1], p%smpd)
            call b%imgs_sym(iptcl)%read(p%stk, iptcl)
        end do
        write(*,'(a)') '>>> CALCULATING CORRELATIONS'
        if( cline%defined('lp') )then
            if( .not. cline%defined('msk') ) stop 'need mask radius (msk) 4 Fourier corr calc!'
            call calc_cartesian_corrmat(b%imgs_sym, corrmat, p%msk, p%lp)
        else
            if( cline%defined('msk') )then
                call calc_cartesian_corrmat(b%imgs_sym, corrmat, p%msk)
            else
                call calc_cartesian_corrmat(b%imgs_sym, corrmat)
            endif
        endif
        call fopen(funit, status='REPLACE', action='WRITE', file='img_smat.bin', access='STREAM',iostat=io_stat)
        call fileiochk("commander_imgproc; image_smat failed to open image_smat.bin ", io_stat)
        write(unit=funit,pos=1,iostat=io_stat) corrmat
        call fileiochk("commander_imgproc; image_smat failed to write image_smat.bin ", io_stat)
        call fclose(funit,errmsg="commander_imgproc; image_smat failed to close image_smat.bin ")
        ! end gracefully
        call simple_end('**** SIMPLE_IMAGE_SMAT NORMAL STOP ****')
    end subroutine exec_image_smat

    !> volume/image_diff is a program for creating a volume based on volume difference
    !! this is purely for direct comparison of images in debugging
    subroutine exec_image_diff(self, cline)
        class(image_diff_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params)         :: p
        type(build)          :: b
        p = params(cline, .false.) ! constants & derived constants produced
        call b%build_general_tbox(p, cline, do3d=.false., nooritab=.true.) ! general objects built (no oritab reading)
        if( cline%defined('stk') .and. cline%defined('vol1') )stop 'Cannot operate on images AND volume at once'
        call diff_imgfiles(p%stk, p%stk2, p%outstk, p%smpd)
        ! end gracefully
        call simple_end('**** SIMPLE_IMAGE_DIFF NORMAL STOP ****')
    end subroutine exec_image_diff

    !> normalize is a program for normalization of MRC or SPIDER stacks and volumes.
    !! If you want to normalize your images inputted with stk, set norm=yes.
    !! hfun (e.g. hfun=sigm) controls the normalization function. If you want to
    !! perform noise normalization of the images set noise_norm=yes given a mask
    !! radius msk (pixels). If you want to normalize your images or volume
    !! (vol1) with respect to their power spectrum set shell_norm=yes
    subroutine exec_normalize( self, cline )
        class(normalize_commander), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        type(build)       :: b
        type(params)      :: p
        real, allocatable :: spec(:)
        integer           :: k
        p = params(cline) ! parameters generated
        if( cline%defined('stk')  .and. cline%defined('vol1') )stop 'Cannot operate on images AND volume at once'
        if( p%norm.eq.'yes'       .and. p%noise_norm.eq.'yes' )stop 'Invalid normalization type'
        if( p%norm.eq.'yes'       .and. p%shellnorm .eq.'yes' )stop 'Invalid normalization type'
        if( p%noise_norm.eq.'yes' .and. p%shellnorm .eq.'yes' )stop 'Invalid normalization type'
        if( cline%defined('stk') )then
            ! 2D
            call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
            if( p%norm.eq.'yes' )then
                ! Normalization
                call norm_imgfile(p%stk, p%outstk, p%smpd)
            else if( p%noise_norm.eq.'yes' )then
                ! Noise normalization
                if( cline%defined('msk') )then
                    call noise_norm_imgfile(p%stk, p%msk, p%outstk, p%smpd)
                else
                    stop 'need msk parameter for noise normalization'
                endif
            else if( p%shellnorm.eq.'yes' )then
                ! shell normalization
                print *,'in'
                call shellnorm_imgfile( p%stk, p%outstk, p%smpd)
            endif
        else if( cline%defined('vol1') )then
            ! 3D
            call b%build_general_tbox(p, cline)         ! general objects built
            call b%vol%new([p%box,p%box,p%box], p%smpd) ! reallocate vol (boxmatch issue)
            if( .not.file_exists(p%vols(1)%str) )stop 'Cannot find input volume'
            call b%vol%read(p%vols(1)%str)
            if( p%norm.eq.'yes' )then
                call b%vol%norm()
                call b%vol%write(p%outvol, del_if_exists=.true.)
            else if( p%shellnorm.eq.'yes' )then
                ! shell normalization
                call b%vol%shellnorm()
                call b%vol%spectrum('power', spec)
                do k=1,size(spec)
                    print *, k, spec(k)
                end do
                call b%vol%write(p%outvol, del_if_exists=.true.)
            else
                stop 'Normalization type not implemented yet'
            endif
        else
            stop 'No input images(s) or volume provided'
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_NORMALIZE NORMAL STOP ****')
    end subroutine exec_normalize

    !> provides re-scaling and clipping routines for MRC or SPIDER stacks and volumes
    subroutine exec_scale( self, cline )
        use simple_qsys_funs,      only: qsys_job_finished
        class(scale_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        type(image)  :: vol2, img, img2
        real         :: ave, sdev, var, med, smpd_new, smpds_new(2), scale
        integer      :: ldim(3), ldim_scaled(3), nfiles, nframes, iframe, ifile
        integer      :: ldims_scaled(2,3)
        character(len=:), allocatable      :: fname
        character(len=STDLEN), allocatable :: filenames(:)
        p = params(cline)                     ! parameters generated
        call img%new([p%box,p%box,1],p%smpd)  ! image created
        call img2%new([p%box,p%box,1],p%smpd) ! image created
        if( cline%defined('stk') .and. cline%defined('vol1') ) stop 'Cannot operate on images AND volume at once'
        if( cline%defined('stk') )then
            ! 2D
            call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
            if( cline%defined('scale2') )then
                call del_file(p%outstk)
                call del_file(p%outstk2)
                ! Rescaling, double
                if( cline%defined('clip')        ) stop 'clip is not allowed in double scaling'
                if( .not. cline%defined('scale') ) stop 'need scale to be part of command line as well 4 double scaling'
                ldims_scaled(1,:) = [p%newbox,p%newbox,1]   ! dimension of scaled
                ldims_scaled(2,:) = [p%newbox2,p%newbox2,1] ! dimension of scaled
                call resize_imgfile_double(p%stk, p%outstk, p%outstk2, p%smpd, ldims_scaled, smpds_new)
                write(*,'(a,1x,f9.4)') 'SAMPLING DISTANCE AFTER SCALING (OUTSTK) :', smpds_new(1)
                write(*,'(a,1x,f9.4)') 'SAMPLING DISTANCE AFTER SCALING (OUTSTK2):', smpds_new(2)
                write(*,'(a,1x,i5)')   'BOX SIZE AFTER SCALING (OUTSTK) :', ldims_scaled(1,1)
                write(*,'(a,1x,i5)')   'BOX SIZE AFTER SCALING (OUTSTK2):', ldims_scaled(2,1)
            else if( cline%defined('scale') .or. cline%defined('newbox') )then
                call del_file(p%outstk)
                ! Rescaling
                ldim_scaled = [p%newbox,p%newbox,1] ! dimension of scaled
                if( cline%defined('clip') )then
                    call resize_and_clip_imgfile(p%stk,p%outstk,p%smpd,ldim_scaled,&
                    [p%clip,p%clip,1],smpd_new)
                else
                    call resize_imgfile(p%stk,p%outstk,p%smpd,ldim_scaled,smpd_new)
                    write(*,'(a,1x,i5)') 'BOX SIZE AFTER SCALING:', ldim_scaled(1)
                endif
                write(*,'(a,1x,f9.4)') 'SAMPLING DISTANCE AFTER SCALING:', smpd_new
            else if( cline%defined('clip') )then
                call del_file(p%outstk)
                ! Clipping
                if( p%clip > p%box )then
                    call pad_imgfile(p%stk,p%outstk,[p%clip,p%clip,1],p%smpd)
                else
                    call clip_imgfile(p%stk,p%outstk,[p%clip,p%clip,1],p%smpd)
                endif
            endif
        else if( cline%defined('vol1') )then
            ! 3D
            call b%build_general_tbox(p, cline) ! general objects built
            ! reallocate vol (boxmatch issue)
            call b%vol%new([p%box,p%box,p%box], p%smpd)
            if( .not.file_exists(p%vols(1)%str) ) stop 'Cannot find input volume'
            call b%vol%read(p%vols(1)%str)
            if( cline%defined('scale') .or. cline%defined('newbox') )then
                ! Rescaling
                call vol2%new([p%newbox,p%newbox,p%newbox],p%smpd)
                call b%vol%fft()
                call vol2%set_ft(.true.)
                if( p%newbox < p%box )then
                    call b%vol%clip(vol2)
                else if( p%newbox > p%box )then
                    call b%vol%pad(vol2)
                else
                    call vol2%copy(b%vol)
                endif
                call b%vol%copy(vol2)
                call b%vol%ifft()
                scale = real(p%newbox)/real(p%box)
                p%box = p%newbox
                smpd_new = p%smpd/scale
                write(*,'(a,1x,f9.4)') 'SAMPLING DISTANCE AFTER SCALING:', smpd_new
            endif
            if( cline%defined('clip') )then
                ! Clipping
                call vol2%new([p%clip,p%clip,p%clip],p%smpd)
                if( p%clip < p%box )then
                    call b%vol%clip(vol2)
                else
                    call b%vol%stats( 'background', ave, sdev, var, med )
                    call b%vol%pad(vol2, backgr=med)
                endif
                call b%vol%copy(vol2)
            else
                 write(*,'(a,1x,i5)') 'BOX SIZE AFTER SCALING:', p%newbox
            endif
            if( p%outvol .ne. '' )call b%vol%write(p%outvol, del_if_exists=.true.)
        else if( cline%defined('filetab') )then
            call read_filetable(p%filetab, filenames)
            nfiles = size(filenames)
            call find_ldim_nptcls(filenames(1),ldim,nframes)
            ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
            if( cline%defined('newbox') )then
                ldim_scaled = [p%newbox,p%newbox,1] ! dimension of scaled
                p%scale = real(p%newbox) / real(ldim(1))
            else if( cline%defined('scale') )then
                ldim_scaled(1) = round2even(real(ldim(1))*p%scale)
                ldim_scaled(2) = round2even(real(ldim(2))*p%scale)
                ldim_scaled(3)   = 1
            else
                stop 'filetab key only in combination with scale or newbox!'
            endif
            call img%new(ldim,p%smpd)
            call img2%new(ldim_scaled,p%smpd/p%scale)
            do ifile=1,nfiles
                call progress(ifile, nfiles)
                if( cline%defined('dir_target') )then
                    fname = trim(p%dir_target)//'/'//add2fbody(basename(filenames(ifile)), p%ext, SCALE_SUFFIX)
                else
                    fname = add2fbody(trim(filenames(ifile)), p%ext, '_sc')
                endif
                call find_ldim_nptcls(filenames(ifile),ldim,nframes)
                ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
                do iframe= 1, nframes
                    call img%read(filenames(ifile), iframe)
                    call img%fft()
                    if( ldim_scaled(1) <= ldim(1) .and. ldim_scaled(2) <= ldim(2) .and. ldim_scaled(3) <= ldim(3) )then
                        call img%clip(img2)
                    else
                        call img%pad(img2)
                    endif
                    call img2%ifft()
                    call img2%write(fname, iframe)
                end do
                deallocate(fname)
            end do
        else
            stop 'SIMPLE_SCALE needs input image(s) or volume or filetable!'
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_SCALE NORMAL STOP ****', print_simple=.false.)
        call qsys_job_finished( p, 'simple_commander_imgproc :: exec_scale' )
    end subroutine exec_scale

    !>  for stacking individual images or multiple stacks into one
    subroutine exec_stack( self, cline )
        class(stack_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(params)                          :: p
        type(build)                           :: b
        integer                               :: nfiles, ldim(3), ifile, ifoo, cnt
        integer                               :: lfoo(3), nimgs, iimg
        character(len=STDLEN), allocatable    :: filenames(:)
        type(image)                           :: tmp, frameimg
        real                                  :: mm(2)
        if( cline%defined('lp') )then
            if( .not. cline%defined('smpd') ) stop 'smpd (sampling distance) needs to be defined if lp is'
        endif
        p = params(cline) ! constants & derived constants produced
        call b%build_general_tbox(p,cline,do3d=.false.) ! general stuff built
        call read_filetable(p%filetab, filenames)
        nfiles = size(filenames)
        call find_ldim_nptcls(filenames(1),ldim,ifoo)
        ldim(3) = 1 ! to correct for the stupid 3:d dim of mrc stacks
        p%box = ldim(1)
        ! prepare b%img and tmp for reading
        call b%img%new([p%box,p%box,1], p%smpd)
        if( cline%defined('clip') ) call tmp%new([p%clip,p%clip,1], p%smpd)
        ! loop over files
        cnt = 0
        do ifile=1,nfiles
            if( .not. file_exists(filenames(ifile)) )then
                write(*,*) 'inputted spec file does not exist: ', trim(adjustl(filenames(ifile)))
            endif
            call find_ldim_nptcls(filenames(ifile),lfoo,nimgs)
            do iimg=1,nimgs
                cnt = cnt+1
                call b%img%read(filenames(ifile), iimg, readhead=.false.)
                if( cline%defined('clip') )then
                    call b%img%clip(tmp)
                    mm = tmp%minmax()
                    DebugPrint 'min/max: ', mm(1), mm(2)
                    call tmp%write(p%outstk, cnt)
                else
                    call b%img%write(p%outstk, cnt)
                endif
            end do
            call progress(ifile, nfiles)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_STACK NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_stack

    !> provides standard single-particle image
    !> processing routines that are applied to MRC or SPIDER stacks.
    subroutine exec_stackops( self, cline )
        use simple_oris, only: oris
        use simple_stackops
        class(stackops_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(params)                             :: p
        type(build)                              :: b
        type(ran_tabu)                           :: rt
        type(oris)                               :: o_here, os_ran
        integer,          allocatable            :: pinds(:)
        character(len=:), allocatable            :: fname
        integer :: i, s, ipst, cnt, cnt2, cnt3, nincl
        real    :: p_ctf, p_dfx
        p = params(cline)                                ! parameters generated
        p%boxmatch = p%box
        call b%build_general_tbox(p, cline,do3d=.false.) ! general objects built
        ! random selection
        if( cline%defined('nran') )then
            write(*,'(a)') '>>> RANDOMLY SELECTING IMAGES'
            allocate( pinds(p%nran), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk('In: simple_commander; stackops',alloc_stat)
            rt = ran_tabu(p%nptcls)
            call rt%ne_ran_iarr(pinds)
            if( cline%defined('oritab') .or. cline%defined('deftab') )then
                call del_file(p%outfile)
            endif
            call os_ran%new(p%nran)
            do i=1,p%nran
                call progress(i, p%nran)
                call b%img%read(p%stk, pinds(i))
                call b%img%write(p%outstk, i)
                if( cline%defined('oritab') .or. cline%defined('deftab') )then
                    call os_ran%set_ori(i, b%a%get_ori(pinds(i)))
                endif
            end do
            if( cline%defined('oritab') .or. cline%defined('deftab') )then
                call b%a%write(p%outfile, [1,p%nran])
            endif
            goto 999
        endif
        ! fishing expeditions
        ! frac only
        if( cline%defined('frac') )then
            if( p%oritab == '' ) stop 'need input orientation doc for fishing expedition; simple_stackops'
            ! determine how many particles to include
            if( p%frac < 0.99 )then
                nincl = nint(real(p%nptcls)*p%frac)
            else
                nincl = p%nptcls
            endif
            ! order the particles
            pinds = b%a%order()
            ! fish the best ones out
            if( cline%defined('state') )then
                cnt = 0
                do i=1,nincl
                    call progress(i, nincl)
                    s = nint(b%a%get(pinds(i), 'state'))
                    if( s == p%state )then
                        cnt = cnt+1
                        call b%img%read(p%stk, pinds(i))
                        call b%img%write(p%outstk, cnt)
                    endif
                end do
                ! make orientation structure for the best ones
                o_here = oris(cnt)
                cnt = 0
                do i=1,nincl
                    call progress(i, nincl)
                    s = nint(b%a%get(pinds(i), 'state'))
                    if( s == p%state )then
                        cnt = cnt+1
                        call o_here%set_ori(cnt, b%a%get_ori(pinds(i)))
                    endif
                end do
                allocate(fname, source='extracted_oris_state'//int2str_pad(p%state,2)//trim(TXT_EXT))
            else if( cline%defined('class') )then
                cnt = 0
                do i=1,nincl
                    call progress(i, nincl)
                    s = nint(b%a%get(pinds(i), 'class'))
                    if( s == p%class )then
                        cnt = cnt+1
                        call b%img%read(p%stk, pinds(i))
                        call b%img%write(p%outstk, cnt)
                    endif
                end do
                ! make orientation structure for the best ones
                o_here = oris(cnt)
                cnt = 0
                do i=1,nincl
                    call progress(i, nincl)
                    s = nint(b%a%get(pinds(i), 'class'))
                    if( s == p%class )then
                        cnt = cnt+1
                        call o_here%set_ori(cnt, b%a%get_ori(pinds(i)))
                    endif
                end do
                allocate(fname, source='extracted_oris_class'//int2str_pad(p%class,5)//trim(TXT_EXT))
            else
                o_here = oris(nincl)
                do i=1,nincl
                    call progress(i, nincl)
                    call b%img%read(p%stk, pinds(i))
                    call b%img%write(p%outstk, i)
                    call o_here%set_ori(i, b%a%get_ori(pinds(i)))
                end do
                allocate(fname, source='extracted_oris'//trim(TXT_EXT))
            endif
            if( cline%defined('outfile') )then
                call o_here%write(p%outfile, [1,o_here%get_noris()])
            else
                call o_here%write(fname, [1,o_here%get_noris()])
            endif
            goto 999
        endif
        ! state/class + frac
        if( (cline%defined('state').or.cline%defined('class')) .and. .not.cline%defined('frac') )then
            if( p%oritab == '' ) stop 'need input orientation doc for fishing expedition; simple_stackops'
            if( cline%defined('state') )then
                cnt = 0
                do i=1,p%nptcls
                    call progress(i, p%nptcls)
                    s = nint(b%a%get(i, 'state'))
                    if( s == p%state )then
                        cnt = cnt+1
                        call b%img%read(p%stk, i)
                        call b%img%write(p%outstk, cnt)
                    endif
                end do
                ! make orientation structure for the extracted ones
                o_here = oris(cnt)
                cnt = 0
                do i=1,p%nptcls
                    call progress(i, p%nptcls)
                    s = nint(b%a%get(i, 'state'))
                    if( s == p%state )then
                        cnt = cnt+1
                        call o_here%set_ori(cnt, b%a%get_ori(i))
                    endif
                end do
                allocate(fname, source='extracted_oris_state'//int2str_pad(p%state,2)//trim(TXT_EXT))
            else if( cline%defined('class') )then
                cnt = 0
                do i=1,p%nptcls
                    call progress(i, p%nptcls)
                    s = nint(b%a%get(i, 'class'))
                    if( s == p%class )then
                        cnt = cnt+1
                        call b%img%read(p%stk, i)
                        call b%img%write(p%outstk, cnt)
                    endif
                end do
                ! make orientation structure for the best ones
                o_here = oris(cnt)
                cnt = 0
                do i=1,p%nptcls
                    call progress(i, p%nptcls)
                    s = nint(b%a%get(i, 'class'))
                    if( s == p%class )then
                        cnt = cnt+1
                        call o_here%set_ori(cnt, b%a%get_ori(i))
                    endif
                end do
                allocate(fname, source='extracted_oris_class'//int2str_pad(p%class,5)//trim(TXT_EXT))
            endif
            if( cline%defined('outfile') )then
                call o_here%write(p%outfile, [1,o_here%get_noris()])
            else
                call o_here%write(fname, [1,o_here%get_noris()])
            endif
            goto 999
        endif
        ! invert contrast
        if( p%neg .eq. 'yes' )then
            call neg_imgfile(p%stk, p%outstk, p%smpd)
            goto 999
        endif
        ! auto correlation function
        if( p%acf .eq. 'yes' )then
            call acf_stack(p%stk, p%outstk)
            goto 999
        endif
        ! produce image statistics
        if( p%stats .eq. 'yes' )then
            call stats_imgfile(p%stk, b%a)
            call b%a%write('image_statistics.txt')
            goto 999
        endif
        ! create frame averages
        if( p%nframesgrp > 0 )then
            call frameavg_stack(p%stk, p%outstk, p%nframesgrp, p%smpd)
            goto 999
        endif
        ! visualize
        if( p%vis .eq. 'yes' )then
            do i=1,p%nptcls
                call b%img%read(p%stk, i)
                call b%img%vis()
            end do
            goto 999
        endif
        ! average
        if( p%avg .eq. 'yes' )then
            call make_avg_stack(p%stk, p%outstk, p%smpd)
            goto 999
        endif
        ! add noise
        if( cline%defined('snr') )then
            call add_noise_imgfile(p%stk, p%outstk, p%snr, p%smpd)
            goto 999
        endif
        ! copy
        if( cline%defined('top') .and. .not. cline%defined('part') )then
            call copy_imgfile(p%stk, p%outstk, p%smpd, fromto=[p%fromp,p%top])
            goto 999
        endif
        ! default
        write(*,*)'Nothing to do!'
        ! end gracefully
    999 call simple_end('**** SIMPLE_STACKOPS NORMAL STOP ****')
    end subroutine exec_stackops

end module simple_commander_imgproc
