!==Class simple_commander_imgproc
!
! This class contains the set of concrete general image processing commanders of the SIMPLE library. This class provides the glue 
! between the reciver (main reciever is simple_exec program) and the abstract action, which is simply execute (defined by the base 
! class: simple_commander_base). Later we can use the composite pattern to create MacroCommanders (or workflows)
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Cyril Reboul & Hans Elmlund 2016
!
module simple_commander_imgproc
use simple_defs
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_strings,        only: int2str, int2str_pad
use simple_filehandling    ! use all in there
use simple_jiffys          ! use all in there
implicit none

public :: binarise_commander
public :: convert_commander
public :: corrcompare_commander
public :: ctfops_commander
public :: filter_commander
public :: image_smat_commander
public :: norm_commander
public :: respimg_commander
public :: scale_commander
public :: stack_commander
public :: stackops_commander
public :: fixmapheader_commander
private

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
type, extends(commander_base) :: norm_commander
  contains
    procedure :: execute      => exec_norm
end type norm_commander
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
type, extends(commander_base) :: respimg_commander
  contains
    procedure :: execute      => exec_respimg
end type respimg_commander
type, extends(commander_base) :: fixmapheader_commander
  contains
    procedure :: execute      => exec_fixmapheader
end type fixmapheader_commander

contains
    
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
            call b%build_general_tbox(p, cline)              ! general objects built
            call doit(b%vol)
            call b%vol%read(p%vols(1))
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_BINARISE NORMAL STOP ****')
        
        contains
            
            subroutine doit( img_or_vol )
                use simple_image, only: image
                class(image), intent(inout) :: img_or_vol
                if( cline%defined('thres') )then
                    call img_or_vol%bin(p%thres)
                else if( cline%defined('npix') )then
                    call img_or_vol%bin(p%npix)
                else
                    call img_or_vol%bin
                endif
                write(*,'(a,1x,i9)') 'NO FOREGROUND PIXELS:', img_or_vol%nforeground()
                write(*,'(a,1x,i9)') 'NO BACKGROUND PIXELS:', img_or_vol%nbackground()
                if( cline%defined('grow') )then
                    do igrow=1,p%grow
                        call img_or_vol%grow_bin
                    end do
                endif
                if( cline%defined('edge') ) call img_or_vol%cos_edge(p%edge)
                if( cline%defined('neg')  ) call img_or_vol%bin_inv
            end subroutine
                
    end subroutine exec_binarise
    
    subroutine exec_convert( self, cline )
        class(convert_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params), target :: p
        type(build),  target :: b
        integer              :: iptcl
        p = params(cline, allow_mix=.true.) ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        if( cline%defined('stk') )then
            do iptcl=1,p%nptcls
                call progress(iptcl, p%nptcls)
                call b%img%read(p%stk, iptcl)
                call b%img%write(p%outstk, iptcl)
            end do 
        else if( cline%defined('vol1') )then
            call b%vol%read(p%vols(1))
            call b%img%write(p%outvol)
        else
            stop 'either vol1 or stk argument required to execute simple_convert'
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CONVERT NORMAL STOP ****')
    end subroutine exec_convert

    subroutine exec_corrcompare( self, cline )
        use simple_image, only: image
        use simple_stat,  only: moment
        use simple_math,  only: get_resolution
        class(corrcompare_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)      :: p
        type(build)       :: b
        integer           :: alloc_stat, npix, iptcl, j
        real              :: corr, ave, sdev, var, fsc05, fsc0143
        real, allocatable :: res(:), corrs(:), corrs_sum(:)
        logical           :: err
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        if( cline%defined('msk') )then
            if( p%stats .eq. 'yes' )then
                allocate(corrs(p%nptcls), stat=alloc_stat)
                call alloc_err('In: simple_corrcompare', alloc_stat)
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
                        write(*,'(A,1X,F7.3)') '>>> REAL-SPACE CORRELATION:', b%img%real_corr(b%img_copy)
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
                call b%img%fsc(b%img_copy, res, corrs)
                if( .not. allocated(corrs_sum) )then
                    allocate(corrs_sum(size(corrs)))
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

    subroutine exec_ctfops( self, cline )
        use simple_procimgfile, only: apply_ctf_imgfile
        use simple_ctf,         only: ctf
        class(ctfops_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        type(ctf)    :: tfun
        logical      :: debug = .false.
        real         :: dfx, dfy, angast
        p = params(cline)                     ! parameters generated
        call b%build_general_tbox(p, cline)   ! general objects built
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
        use simple_procimgfile, only: bp_imgfile, phase_rand_imgfile
        class(filter_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        if( cline%defined('stk') )then
            ! 2D
            if( .not.file_exists(p%stk) )stop 'Cannot find input stack (stk)'
            if( p%phrand .eq. 'no')then
                ! Band pass
                if( cline%defined('lp') .and. cline%defined('hp') )then
                    call bp_imgfile(p%stk, p%outstk, p%smpd, p%hp, p%lp)
                else if( cline%defined('lp') )then
                    call bp_imgfile(p%stk, p%outstk, p%smpd, 0., p%lp)
                else if( cline%defined('hp') )then
                    call bp_imgfile(p%stk, p%outstk, p%smpd, p%hp, 0.)
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
            if( .not.file_exists(p%vols(1)) )stop 'Cannot find input volume (vol1)'
            call b%vol%read(p%vols(1))
            if( p%phrand.eq.'no')then
                if( cline%defined('bfac') )then
                ! bfactor
                    call b%vol%apply_bfac(p%bfac)
                ! Band pass
                else if( cline%defined('hp') .and. cline%defined('lp') )then
                    call b%vol%bp(p%hp,p%lp)
                else if( cline%defined('hp') )then
                    call b%vol%bp(p%hp,0.)
                else if( cline%defined('lp') )then
                    call b%vol%bp(0.,p%lp)
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
    
    subroutine exec_image_smat(self, cline)
        use simple_corrmat  ! use all in there
        use simple_ori,     only: ori
        use simple_imgfile, only: imgfile
        use simple_image,   only: image
        class(image_smat_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params)         :: p
        type(build)          :: b
        integer              :: iptcl, alloc_stat, funit, io_stat
        real, allocatable    :: corrmat(:,:)
        logical              :: debug=.false.
        p = params(cline, .false.)                           ! constants & derived constants produced
        call b%build_general_tbox(p, cline, .false., .true.) ! general objects built (no oritab reading)
        allocate(b%imgs_sym(p%nptcls), stat=alloc_stat)
        call alloc_err('In: simple_image_smat, 1', alloc_stat)
        do iptcl=1,p%nptcls
            call b%imgs_sym(iptcl)%new([p%box,p%box,1], p%smpd, p%imgkind)
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
        funit = get_fileunit()
        open(unit=funit, status='REPLACE', action='WRITE', file='img_smat.bin', access='STREAM')
        write(unit=funit,pos=1,iostat=io_stat) corrmat
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to image_smat.bin'
            stop 'I/O error; simple_image_smat'
        endif
        close(funit)
        ! end gracefully
        call simple_end('**** SIMPLE_IMAGE_SMAT NORMAL STOP ****')
    end subroutine exec_image_smat

    subroutine exec_norm( self, cline )
        use simple_procimgfile,   only: norm_imgfile, noise_norm_imgfile, shellnorm_imgfile
        class(norm_commander), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        type(build)       :: b
        type(params)      :: p
        real, allocatable :: spec(:)
        integer           :: k
        p = params(cline)                           ! parameters generated
        call b%build_general_tbox(p, cline)         ! general objects built
        if( cline%defined('stk') .and. cline%defined('vol1') )stop 'Cannot operate on images AND volume at once'
        if( p%norm.eq.'yes' .and. p%noise_norm.eq.'yes' )stop 'Invalid normalization type'
        if( p%norm.eq.'yes' .and. p%shellnorm.eq.'yes' )stop 'Invalid normalization type'
        if( p%noise_norm.eq.'yes' .and. p%shellnorm.eq.'yes' )stop 'Invalid normalization type'
        call b%vol%new([p%box,p%box,p%box], p%smpd) ! reallocate vol (boxmatch issue)
        if( cline%defined('stk') )then
            ! 2D
            if( p%norm.eq.'yes' )then
                ! Normalization
                if( cline%defined('hfun') )then
                    call norm_imgfile(p%stk, p%outstk, p%smpd, hfun=p%hfun)
                else
                    call norm_imgfile(p%stk, p%outstk, p%smpd)
                endif
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
            if( .not.file_exists(p%vols(1)) )stop 'Cannot find input volume'
            call b%vol%read(p%vols(1))
            if( p%shellnorm.eq.'yes' )then
                ! shell normalization
                call b%vol%shellnorm
                spec = b%vol%spectrum('power')
                do k=1,size(spec)
                    print *, k, spec(k)
                end do
                if( p%outvol .ne. '' )call b%vol%write(p%outvol, del_if_exists=.true.)
            else
                stop 'Normalization type not implemented yet'
            endif
        else
            stop 'No input images(s) or volume provided'
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_NORM NORMAL STOP ****')
    end subroutine exec_norm
    
    subroutine exec_scale( self, cline )
        use simple_procimgfile  ! use all in there
        use simple_image,       only: image
        use simple_math,        only: round2even
        class(scale_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        type(image)  :: vol2, img, img2
        real         :: ave, sdev, var, med, smpd_new, smpds_new(2), scale
        integer      :: ldim(3), ldim_scaled(3), nfiles, nframes, iframe, ifile
        integer      :: ldims_scaled(2,3), ldims_clip(2,3)
        character(len=:), allocatable      :: fname
        character(len=STDLEN), allocatable :: filenames(:)
        p = params(cline)                               ! parameters generated
        call b%build_general_tbox(p, cline)             ! general objects built
        call b%vol%new([p%box,p%box,p%box], p%smpd)     ! reallocate vol (boxmatch issue)
        call img%new([p%box,p%box,1],p%smpd,p%imgkind)  ! image created
        call img2%new([p%box,p%box,1],p%smpd,p%imgkind) ! image created
        if( cline%defined('stk') .and. cline%defined('vol1') )stop 'Cannot operate on images AND volume at once'
        if( cline%defined('stk') )then
            ! 2D
            if( cline%defined('scale2') )then
                ! Rescaling, double
                if( cline%defined('clip')  ) stop 'clip is not allowed in double scaling'
                if( .not. cline%defined('scale') ) stop 'need scale to be part of command line as well 4 double scaling'
                ldims_scaled(1,:) = [p%newbox,p%newbox,1]   ! dimension of scaled
                ldims_scaled(2,:) = [p%newbox2,p%newbox2,1] ! dimension of scaled
                if( cline%defined('part') )then
                    p%outstk  = 'outstk_part'//int2str_pad(p%part, p%numlen)//p%ext
                    p%outstk2 = 'outstk2_part'//int2str_pad(p%part, p%numlen)//p%ext
                    call resize_imgfile_double(p%stk, p%outstk, p%outstk2,&
                    p%smpd, ldims_scaled, smpds_new, [p%fromp,p%top])
                else
                    call resize_imgfile_double(p%stk, p%outstk, p%outstk2,&
                    p%smpd, ldims_scaled, smpds_new)
                endif
                write(*,'(a,1x,f9.4)') 'SAMPLING DISTANCE AFTER SCALING (OUTSTK) :', smpds_new(1)
                write(*,'(a,1x,f9.4)') 'SAMPLING DISTANCE AFTER SCALING (OUTSTK2):', smpds_new(2)
                write(*,'(a,1x,i5)')   'BOX SIZE AFTER SCALING (OUTSTK) :', ldims_scaled(1,1)
                write(*,'(a,1x,i5)')   'BOX SIZE AFTER SCALING (OUTSTK2):', ldims_scaled(2,1)
            else if( cline%defined('scale') .or. cline%defined('newbox') )then
                ! Rescaling
                ldim_scaled = [p%newbox,p%newbox,1] ! dimension of scaled
                if( cline%defined('clip') )then
                    if( cline%defined('part') )then
                        p%outstk = 'outstk_part'//int2str_pad(p%part, p%numlen)//p%ext
                        call resize_and_clip_imgfile(p%stk,p%outstk,p%smpd,ldim_scaled,&
                        &[p%clip,p%clip,1],smpd_new,[p%fromp,p%top])
                    else
                        call resize_and_clip_imgfile(p%stk,p%outstk,p%smpd,ldim_scaled,&
                        [p%clip,p%clip,1],smpd_new)
                    endif
                else
                    if( cline%defined('part') )then
                        p%outstk = 'outstk_part'//int2str_pad(p%part, p%numlen)//p%ext
                        call resize_imgfile(p%stk,p%outstk,p%smpd,ldim_scaled,smpd_new,[p%fromp,p%top])
                    else
                        call resize_imgfile(p%stk,p%outstk,p%smpd,ldim_scaled,smpd_new)
                    endif
                    write(*,'(a,1x,i5)') 'BOX SIZE AFTER SCALING:', ldim_scaled(1)
                endif
                write(*,'(a,1x,f9.4)') 'SAMPLING DISTANCE AFTER SCALING:', smpd_new
            else if( cline%defined('clip') )then
                ! Clipping
                call clip_imgfile(p%stk,p%outstk,[p%clip,p%clip,1],p%smpd)
            endif
        else if( cline%defined('vol1') )then
            ! 3D
            if( .not.file_exists(p%vols(1)) ) stop 'Cannot find input volume'
            call b%vol%read(p%vols(1))
            if( cline%defined('scale') .or. cline%defined('newbox') )then
                ! Rescaling
                call vol2%new([p%newbox,p%newbox,p%newbox],p%smpd)
                call b%vol%fwd_ft
                call vol2%set_ft(.true.)
                if( p%newbox < p%box )then
                    call b%vol%clip(vol2)
                else if( p%newbox > p%box )then
                    call b%vol%pad(vol2)
                else
                    vol2 = b%vol
                endif
                b%vol = vol2
                call b%vol%bwd_ft
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
                    if( cline%defined('msk') )then
                        call b%vol%stats( 'background', ave, sdev, var, med, p%msk ) 
                    else
                        call b%vol%stats( 'background', ave, sdev, var, med ) 
                    endif
                    call b%vol%pad(vol2, backgr=med)
                endif
                b%vol = vol2
            else
                 write(*,'(a,1x,i5)') 'BOX SIZE AFTER SCALING:', p%newbox
            endif
            if( p%outvol .ne. '' )call b%vol%write(p%outvol, del_if_exists=.true.)
        else if( cline%defined('filetab') )then
            call read_filetable(p%filetab, filenames)
            nfiles = size(filenames)
            if( cline%defined('scale') )then
                call find_ldim_nptcls(filenames(1),ldim,nframes)
                ldim(3)          = 1 ! to correct for the stupide 3:d dim of mrc stacks
                ldim_scaled(1) = round2even(real(ldim(1))*p%scale)
                ldim_scaled(2) = round2even(real(ldim(2))*p%scale)
                ldim_scaled(3)   = 1
            else
                stop 'need scale factor for this mode of execution; simple_commander_imgproc :: exec_scale'
            endif
            call img%new(ldim,p%smpd,p%imgkind)
            call img2%new(ldim_scaled,p%smpd/p%scale,p%imgkind)
            do ifile=1,nfiles
                call progress(ifile, nfiles)
                fname = add2fbody(remove_abspath(filenames(ifile)), p%ext, '_sc')
                do iframe=1,nframes
                    call img%read(filenames(ifile), iframe)
                    call img%fwd_ft
                    if( ldim_scaled(1) <= ldim(1) .and. ldim_scaled(2) <= ldim(2) .and. ldim_scaled(3) <= ldim(3) )then
                        call img%clip(img2)
                    else
                        call img%pad(img2)
                    endif
                    call img2%bwd_ft
                    call img2%write(fname, iframe)
                end do
                deallocate(fname)
            end do
        else
            stop 'SIMPLE_SCALE needs input image(s) or volume or filetable!'          
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_SCALE NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_scale
    
    subroutine exec_stack( self, cline )
        use simple_imgfile,      only: imgfile
        use simple_image,        only: image
        class(stack_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(params)                          :: p
        type(build)                           :: b
        integer                               :: nfiles, ldim(3), ifile, ifoo, nmovies, cnt
        integer                               :: imovie, iframe, numlen, lfoo(3), nimgs, iimg
        character(len=STDLEN), allocatable    :: filenames(:)
        character(len=:), allocatable         :: moviename
        type(image)                           :: mask, tmp, frameimg
        real                                  :: mm(2)
        logical, parameter                    :: debug = .false.
        if( cline%defined('lp') )then
            if( .not. cline%defined('smpd') ) stop 'smpd (sampling distance) needs to be defined if lp is'
        endif
        p = params(cline,checkdistr=.false.)            ! constants & derived constants produced
        call b%build_general_tbox(p,cline,do3d=.false.) ! general stuff built
        call read_filetable(p%filetab, filenames)
        nfiles = size(filenames)
        if( debug ) write(*,*) 'read the filenames'
        if( cline%defined('xdim') .and. cline%defined('ydim') )then
            ldim = [p%xdim,p%ydim,1]
        else
            call find_ldim_nptcls(filenames(1),ldim,ifoo)
            ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
        endif
        if( debug ) write(*,*) 'logical dimension: ', ldim
        if( cline%defined('nframes') )then
            if( .not. cline%defined('fbody') ) stop 'need fbody (file body of output stacks) on the command line'
            if( mod(nfiles,p%nframes) .eq. 0 )then
                ! fine, go ahead
            else
                stop 'Number of mrc files not a multiple of nframes!'
            endif
            call frameimg%new(ldim,p%smpd)
            nmovies = nfiles/p%nframes
            if( cline%defined('numlen') )then
                numlen = p%numlen
            else
                numlen = len(int2str(nmovies))
            endif
            ifile = 0
            do imovie=1,nmovies
                allocate(moviename, source=trim(adjustl(p%fbody))//int2str_pad(imovie, numlen)//'.mrcs')
                call del_file(moviename)
                do iframe=1,p%nframes
                    ifile = ifile+1
                    call progress(ifile, nfiles)
                    call frameimg%read(filenames(ifile),1,readhead=.false.)
                    call frameimg%write(moviename,iframe)
                end do
                deallocate(moviename)
            end do
        else
            p%box = ldim(1)
            ! create mask
            if( cline%defined('lp') )then
                call tmp%new([p%clip,p%clip,1], p%smpd)
                tmp = cmplx(1.,0.)
                call tmp%bp(0.,p%lp,0.)
                call tmp%ft2img('real', mask)
                call mask%write('resolution_mask.mrc', 1)
                tmp = 0.0
            endif
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
                        if( debug ) print *, 'min/max: ', mm(1), mm(2)
                        call tmp%write(p%outstk, cnt)
                    else
                        call b%img%write(p%outstk, cnt)
                    endif
                end do
                call progress(ifile, nfiles)
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_STACK NORMAL STOP ****')
    end subroutine exec_stack

    subroutine exec_stackops( self, cline )
        use simple_ran_tabu,    only: ran_tabu
        use simple_procimgfile, only: neg_imgfile, acf_imgfile, frameavg_imgfile
        use simple_procimgfile  ! use all in there
        use simple_image,       only: image
        use simple_oris,        only: oris
        class(stackops_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(params)                             :: p
        type(build)                              :: b
        type(ran_tabu)                           :: rt
        type(image)                              :: img
        type(oris)                               :: o_here
        integer,          allocatable            :: pinds(:)
        character(len=:), allocatable            :: fname
        integer :: i, s, cnt, nincl, alloc_stat, lfoo(3), np1,np2,ntot
        p = params(cline)                               ! parameters generated
        call b%build_general_tbox(p, cline)             ! general objects built
        call img%new([p%box,p%box,1],p%smpd,p%imgkind)  ! image created
        ! random selection
        if( cline%defined('nran') )then
            write(*,'(a)') '>>> RANDOMLY SELECTING IMAGES'
            allocate( pinds(p%nran), stat=alloc_stat )
            call alloc_err('In: simple_commander; stackops', alloc_stat)
            rt = ran_tabu(p%nptcls)
            call rt%ne_ran_iarr(pinds)
            if( cline%defined('oritab') .or. cline%defined('deftab') )then
                call del_file(p%outfile)
            endif
            do i=1,p%nran
                call progress(i, p%nran)
                call img%read(p%stk, pinds(i))
                call img%write(p%outstk, i)
                if( cline%defined('oritab') .or. cline%defined('deftab') )then
                    call b%a%write(pinds(i), p%outfile)
                endif
            end do
            goto 999
        endif
        ! fishing expeditions
        ! order
        if( cline%defined('order') )then
            if( p%order .eq. 'yes' )then
                ! order the particles
                pinds = b%a%order()
                cnt = 0
                do i=1,p%nptcls
                    cnt = cnt + 1
                    call progress(i, p%nptcls)
                    call img%read(p%stk, pinds(i))
                    call img%write(p%outstk, cnt)
                end do
            endif
            goto 999
        endif
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
                        call img%read(p%stk, pinds(i))
                        call img%write(p%outstk, cnt)
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
                allocate(fname, source='extracted_oris_state'//int2str_pad(p%state,2)//'.txt')
            else if( cline%defined('class') )then
                cnt = 0
                do i=1,nincl
                    call progress(i, nincl)
                    s = nint(b%a%get(pinds(i), 'class'))
                    if( s == p%class )then
                        cnt = cnt+1
                        call img%read(p%stk, pinds(i))
                        call img%write(p%outstk, cnt)
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
                allocate(fname, source='extracted_oris_class'//int2str_pad(p%class,5)//'.txt')
            else
                o_here = oris(nincl)
                do i=1,nincl
                    call progress(i, nincl)
                    call img%read(p%stk, pinds(i))
                    call img%write(p%outstk, i)
                    call o_here%set_ori(i, b%a%get_ori(pinds(i)))
                end do
                allocate(fname, source='extracted_oris.txt')
            endif
            call o_here%write(fname)
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
                        call img%read(p%stk, i)
                        call img%write(p%outstk, cnt)
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
                allocate(fname, source='extracted_oris_state'//int2str_pad(p%state,2)//'.txt')
            else if( cline%defined('class') )then
                cnt = 0
                do i=1,p%nptcls
                    call progress(i, p%nptcls)
                    s = nint(b%a%get(i, 'class'))
                    if( s == p%class )then
                        cnt = cnt+1
                        call img%read(p%stk, i)
                        call img%write(p%outstk, cnt)
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
                allocate(fname, source='extracted_oris_class'//int2str_pad(p%class,5)//'.txt')
            endif
            call o_here%write(fname)
            goto 999
        endif
        ! invert contrast
        if( p%neg .eq. 'yes' )then
            call neg_imgfile(p%stk, p%outstk, p%smpd)
            goto 999
        endif
        ! auto correlation function
        if( p%acf .eq. 'yes' )then
            call acf_imgfile(p%stk, p%outstk)
            goto 999
        endif
        ! create frame averages
        if( p%frameavg > 0 )then
            call frameavg_imgfile(p%stk, p%outstk, p%frameavg, p%smpd)
            goto 999
        endif
        ! visualize
        if( p%vis .eq. 'yes' )then
            do i=1,p%nptcls
                call img%read(p%stk, i)
                call img%vis
            end do
            goto 999
        endif
        ! average
        if( p%avg .eq. 'yes' )then
            call make_avg_imgfile(p%stk, p%outstk, p%smpd)
            goto 999
        endif
        ! add noise
        if( cline%defined('snr') )then
            call add_noise_imgfile(p%stk, p%outstk, p%snr, p%smpd)
            goto 999
        endif
        ! APPEND STK2 TO STK WHILE PRESERVING THE NAME OF STK
        if( p%append .eq. 'yes' )then
            if( cline%defined('stk') .and. cline%defined('stk2') )then
                ! find out image dimension and number of particles
                call find_ldim_nptcls(p%stk,lfoo,np1)
                call find_ldim_nptcls(p%stk2,lfoo,np2)
                ntot = np1+np2
                cnt = 0
                do i=np1+1,ntot
                    cnt = cnt+1
                    call progress(cnt,np2)
                    call img%read(p%stk2,cnt)
                    call img%write(p%stk,i)
                end do
            else
                stop 'need two stacks (stk & stk2) to append; simple_stackops'
            endif
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

    subroutine exec_respimg( self, cline )
        use simple_stat,             only: pearsn, normalize_sigm, corrs2weights
        use simple_aff_prop,         only: aff_prop
        use simple_ori,              only: ori
        use simple_image,            only: image
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_hadamard_common,  only: preprefs4align, reset_prev_defparms
        use simple_ctf,              only: ctf
        use simple_projector_hlev,   only: rotimg
        class(respimg_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)                  :: p
        type(build)                   :: b
        type(ori)                     :: o
        type(image), allocatable      :: imgs(:), ctfsqs(:)
        type(image)                   :: ctfsqsum, new_img, last
        type(aff_prop)                :: affprop
        type(polarft_corrcalc)        :: pftcc
        real,    allocatable          :: simmat(:,:), resp(:,:), weights(:)
        integer, allocatable          :: proj_inds(:), centers(:), labels(:)
        real                          :: x, y, simm, w
        integer                       :: iptcl, i,j, proj, proj_pop, ncls, fnr
        real,               parameter :: SHTHRESH=0.0001
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        ! Init
        call b%img%new([p%box,p%box,1],p%smpd)
        call new_img%new([p%box,p%box,1],p%smpd)
        call last%new([p%box,p%box,1],p%smpd)
        call ctfsqsum%new([p%box,p%box,1],p%smpd)
        ncls = maxval( nint(b%a%get_all('class')) )
        write(*,'(A,I6)') '>>> NUMBER OF PROJECTIONS DIRECTIONS', ncls
        ! Main Loop over projdirs
        do proj = 1,ncls
            call progress( proj, ncls )
            proj_pop = b%a%get_cls_pop( proj )
            if( proj_pop == 0 )cycle
            proj_inds = b%a%get_cls_pinds( proj )
            if( proj_pop==1 )then
                ! no weighing
                iptcl = proj_inds( 1 )
                o = b%a%get_ori( iptcl )
                call b%img%read( p%stk, iptcl )
                call prep4rec( b%img, o )
                call b%img%write( p%outstk, iptcl )
                deallocate( proj_inds )
            else if( proj_pop>1 )then
                ! Affinity propagation weighing
                ! init 
                allocate( simmat( proj_pop, proj_pop ), weights( proj_pop ) )
                allocate( imgs(proj_pop), ctfsqs( proj_pop ) )
                simmat = 1.
                do i=1,proj_pop
                    call imgs( i )%new([p%box,p%box,1],p%smpd)
                    call ctfsqs( i )%new([p%box,p%box,1],p%smpd)
                enddo
                ! prep images
                call pftcc%new(proj_pop, [p%fromp,p%top], [p%box,p%box,1], p%kfromto, 1, p%ctf)
                do i = 1,proj_pop
                    iptcl = proj_inds( i )
                    o     = b%a%get_ori( iptcl )
                    call b%img%read( p%stk, iptcl )
                    imgs(i) = b%img
                    ! prep for correlation
                    call prep_img4align( b%img, o )
                    call b%img%img2polarft(i, pftcc, isptcl=.true. )
                    call b%img%img2polarft(i, pftcc, isptcl=.false. )
                enddo
                ! build similarity matrix
                do i=1,proj_pop-1
                    call reset_prev_defparms
                    call preprefs4align(b, p, i, pftcc )
                    do j=i+1,proj_pop
                        simm = pftcc%corr(i,j,1)
                        simmat( i,j ) = simm
                        simmat( j,i ) = simm
                    enddo
                enddo
                ! simmat = -(1.-simmat)**2.
                ! if( proj_pop>2 )then
                !     ! exec affinity propagation
                !     call affprop%new( proj_pop, simmat, lam=0.6 )
                !     call affprop%propagate( centers, labels, simm )
                !     resp = affprop%get_resp()
                !     call affprop%kill
                !     !deallocate( centers, labels)
                ! else
                !     allocate( resp(2,2) )
                !     resp = simmat
                ! endif
                ! prep image & ctfsq
                do i=1,proj_pop
                    o = b%a%get_ori( proj_inds(i) )
                    call prep4rec( imgs(i), o )
                    call prep_ctfsq( ctfsqs(i), o )
                enddo
                ! ! exponential weighing
                ! resp = exp(resp)
                ! ! linear weighing
                ! !forall( i=1:proj_pop )resp(i,:) = resp(i,:)+abs(minval(resp(i,:)))
                ! ! forall( i=1:proj_pop )resp(i,:) = (1.+resp(i,:)) / (1.+maxval(resp(i,:)))
                ! forall( i=1:proj_pop )resp(i,:) = resp(i,:) / sum( resp(i,:) )
                ! generate substitute particles
                do i=1,proj_pop
                    iptcl = proj_inds( i )
                    new_img  = 0.
                    ctfsqsum = cmplx(0.,0.)
                    !call normalize_sigm(resp(i,:))
                    !resp(i,:) = exp(resp(i,:))
                    !resp(i,:) = resp(i,:)/sum(resp(i,:))
                    !weights = corrs2weights( simmat( i,:) )
                    weights = simmat(i,:)**12.
                    print *,proj,i,simmat(i,:)
                    print *,proj,i,weights
                    !print *,proj,i,simmat(i,:)**2.
                    do j=1,proj_pop
                        ! if( proj_pop>2 )then
                        !     if( labels(j).ne.labels(i) )cycle
                        ! endif
                        ! w = resp(i,j)
                        w = weights( j )
                        if( w<= 0. )cycle
                        ! accumulates imgs
                        b%img = imgs(j)
                        call b%img%mul( w )
                        call new_img%add( b%img )
                        ! accumulates CTFsq
                        b%img = ctfsqs(j)
                        call b%img%mul( w )
                        call ctfsqsum%add( b%img )
                    enddo
                    ! CTFsq sum division
                    call new_img%fwd_ft
                    call new_img%ctf_dens_correct( ctfsqsum )
                    call new_img%bwd_ft
                    call new_img%norm
                    if( iptcl.eq.p%nptcls )then
                        last = new_img
                    else
                        call new_img%write( p%outstk, iptcl )
                    endif
                enddo
                ! generate cavg
                ! new_img  = 0.
                ! ctfsqsum = cmplx(0.,0.)
                ! do i=1,proj_pop
                !     call new_img%add( imgs(i) )
                !     !call ctfsqsum%add( ctfsqs(i) )
                ! enddo
                ! ! call new_img%fwd_ft
                ! ! call new_img%ctf_dens_correct( ctfsqsum )
                ! ! call new_img%bwd_ft
                ! call new_img%div( real(proj_pop) )
                ! !call new_img%norm
                ! call new_img%write( 'cavgs.mrc', proj )
                ! deallocation
                if( allocated(weights) )deallocate(weights)
                if( allocated(simmat) )deallocate(simmat)
                if( allocated(centers) )deallocate(centers)
                if( allocated(labels) )deallocate(labels)
                do i=1,proj_pop
                    call imgs(i)%kill
                    call ctfsqs(i)%kill
                enddo
                if( allocated(resp) )deallocate( resp )
                deallocate( proj_inds, imgs, ctfsqs )
            endif
        enddo
        call last%write( p%outstk, p%nptcls )
        ! orientations output
        do i=1,b%a%get_noris()
            o = b%a%get_ori( i )
            call o%set('x',0.)
            call o%set('y',0.)
            call o%e3set( 0. )
            call o%set('dfx',0.)
            call o%set('dfy',0.)
            call o%set('angast',0.)
            call b%a%set_ori( i, o )
        enddo
        call b%a%write( p%outfile )
        call simple_end('**** SIMPLE_RESPIMG NORMAL STOP ****')

        contains

        subroutine prep_img4align( img, o )
            class(image), intent(inout) :: img
            type(ori),    intent(inout) :: o
            type(image) :: img_rot
            type(ctf)   :: tfun
            real        :: x, y, dfx, dfy, angast
            if( p%ctf .ne. 'no' )then
                ! create ctf object
                tfun = ctf(img%get_smpd(),o%get('kv'),o%get('cs'),o%get('fraca'))
                select case(p%tfplan%mode)
                    case('astig') ! astigmatic CTF
                        dfx    = o%get('dfx')
                        dfy    = o%get('dfy')
                        angast = o%get('angast')
                    case('noastig') ! non-astigmatic CTF
                        dfx    = o%get('dfx')
                        dfy    = dfx
                        angast = 0.
                    case DEFAULT
                        write(*,*) 'Unsupported p%tfplan%mode: ', trim(p%tfplan%mode)
                        stop 'simple_hadamard_common :: prepimg4align'
                end select
            endif
            call img%fwd_ft
            select case(p%ctf)
                case('mul')  ! images have been multiplied with the CTF, no CTF-dependent weighting of the correlations
                    stop 'ctf=mul is not supported; simple_hadamard_common :: prepimg4align' 
                case('no')   ! do nothing
                case('yes')  ! do nothing
                case('flip') ! flip back
                    call tfun%apply( img, dfx, 'flip', dfy, angast )
                case DEFAULT
                    stop 'Unsupported ctf mode; simple_hadamard_common :: prepimg4align'
            end select
            ! shift image to rotational origin
            x = o%get('x')
            y = o%get('y')
            if( abs(x) > SHTHRESH .or. abs(y) > SHTHRESH ) call img%shift(-x, -y)
            ! rotation
            call rotimg( img, -o%e3get(), p%msk, img_rot )
            img = img_rot
            call img_rot%kill
            call img%bwd_ft
            call img%mask(p%msk, 'soft')
            call img%fwd_ft
        end subroutine

        subroutine prep4rec( img, o )
            class(image), intent(inout) :: img
            type(ori),   intent(inout) :: o
            type(image) :: img_rot
            type(ctf)   :: tfun
            real :: x, y, dfx, dfy, angast
            call img_rot%new( img%get_ldim(), img%get_smpd() )
            if( p%ctf .ne. 'no' )then
                ! create ctf object
                tfun = ctf(img%get_smpd(),o%get('kv'),o%get('cs'),o%get('fraca'))
                select case(p%tfplan%mode)
                    case('astig') ! astigmatic CTF
                        dfx    = o%get('dfx')
                        dfy    = o%get('dfy')
                        angast = o%get('angast')
                    case('noastig') ! non-astigmatic CTF
                        dfx    = o%get('dfx')
                        dfy    = dfx
                        angast = 0.
                end select
            endif
            call img%fwd_ft
            select case(p%tfplan%flag)
                case('yes')  ! multiply with CTF
                    call tfun%apply(img, dfx, 'ctf', dfy, angast)
                case('flip') ! multiply with abs(CTF)
                    call tfun%apply(img, dfx, 'abs', dfy, angast)
                case('mul','no')  ! do nothing
                    ! doing nothing
                case DEFAULT
                    write(*,*) 'Unsupported ctfflag: ', p%tfplan%flag
                    stop 'simple_filterer :: wiener_restore2D_online'
            end select
            ! shift image to rotational origin
            x = o%get('x')
            y = o%get('y')
            if( abs(x) > SHTHRESH .or. abs(y) > SHTHRESH ) call img%shift(-x, -y)
            ! rotation
            call rotimg(img, -o%e3get(), p%msk, img_rot )
            img = img_rot
            !call img%bwd_ft
            call img_rot%kill
        end subroutine prep4rec

        subroutine prep_ctfsq( img, o )
            type(image), intent(inout) :: img
            type(ori),   intent(inout) :: o
            type(ctf) :: tfun
            real      :: dfx, dfy, angast
            if( p%tfplan%flag .ne. 'no' )then
                ! create ctf object
                tfun = ctf(img%get_smpd(),o%get('kv'),o%get('cs'),o%get('fraca'))
                select case(p%tfplan%mode)
                    case('astig') ! astigmatic CTF
                        dfx    = o%get('dfx')
                        dfy    = o%get('dfy')
                        angast = o%get('angast')
                    case('noastig') ! non-astigmatic CTF
                        dfx    = o%get('dfx')
                        dfy    = dfx
                        angast = 0.
                end select
                call tfun%ctf2img( img, dfx, 'square', dfy, angast)
            endif
        end subroutine prep_ctfsq

    end subroutine exec_respimg

    subroutine exec_fixmapheader( self, cline )
        use simple_image,   only: image
        use simple_imgfile, only: imgfile
        class(fixmapheader_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params) :: p
        type(image)  :: img 
        real         :: dev
        integer      :: ldim(3), maxim
        p = params(cline) ! parameters generated
        ! read vol
        call find_ldim_nptcls(p%vols(1), ldim, maxim)
        call img%new(ldim, p%smpd)
        call img%read(p%vols(1))
        ! set RMSD
        dev = img%rmsd()
        ! header will be sorted on write
        call img%write(p%outvol, rmsd=dev)
    end subroutine exec_fixmapheader

end module simple_commander_imgproc
