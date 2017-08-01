
! this is the set of concrete commanders that constitute the SIMPLE library. This class provides the glue between the reciver (main reciever is simple_exec program)
! and the abstract action, which is simply execute (defined by the base class: simple_commander_base). Later we can use the composite pattern to create MacroCommanders
! (or workflows)

module simple_commander
use simple_defs            ! singleton
use simple_jiffys          ! singleton
use simple_timing          ! singleton
use simple_cuda            ! singleton
use simple_cuda_defs       ! singleton
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
implicit none

!********* PUBLIC INTERFACES

! SIMULATORS
public :: simimgs_commander
public :: simmovie_commander
public :: simsubtomo_commander

! PRE-PROCESSING METHODS
public :: select_frames_commander
public :: boxconvs_commander
public :: integrate_movies_commander
public :: powerspecs_commander
public :: unblur_movies_commander
public :: stack_powerspecs_commander
public :: select_commander
public :: extr_ptcls_commander

! PRIME2D METHODS
public :: prime2D_init_commander
public :: prime2D_commander
public :: cavgassemble_commander
public :: check2D_conv_commander
public :: rank_cavgs_commander

! PRIME3D METHODS
public :: resrange_commander
public :: npeaks_commander
public :: nspace_commander
public :: prime3D_init_commander
public :: multiptcl_init_commander
public :: prime3D_commander
public :: cont3D_commander
public :: check3D_conv_commander

! COMMON-LINES METHODS
public :: comlin_smat_commander
public :: symsrch_commander

! MASKING METHODS
public :: mask_commander
public :: automask2D_commander
public :: automask3D_commander

! RECONSTRUCTION METHODS
public :: eo_recvol_commander
public :: eo_volassemble_commander
public :: recvol_commander
public :: volassemble_commander

! CHECKER METHODS
public :: check_box_commander
public :: check_nptcls_commander
public :: iminfo_commander

! VOLOPS METHODS
public :: cenvol_commander
public :: postproc_vol_commander
public :: projvol_commander
public :: volaverager_commander
public :: volops_commander
public :: volume_smat_commander

! MISCELLANOUS METHODS
public :: binarise_commander
public :: cluster_smat_commander
public :: converter_commander
public :: ctfops_commander
public :: filter_commander
public :: gatanmrc2mrc_commander
public :: image_smat_commander
public :: map2ptcls_commander
public :: norm_commander
public :: orisops_commander
public :: print_fsc_commander
public :: res_commander
public :: scale_commander
public :: stackops_commander
public :: tseries_split_commander

! PARALLEL PROCESSING METHODS
public :: merge_algndocs_commander
public :: merge_similarities_commander
public :: split_pairs_commander
public :: split_commander
private

!********* CONRETE TYPE DEFINITIONS

! SIMULATORS (commander_sim)
type, extends(commander_base) :: simimgs_commander
  contains
    procedure :: execute      => exec_simimgs
end type simimgs_commander
type, extends(commander_base) :: simmovie_commander
  contains
    procedure :: execute      => exec_simmovie
end type simmovie_commander
type, extends(commander_base) :: simsubtomo_commander
  contains
    procedure :: execute      => exec_simsubtomo
end type simsubtomo_commander

! PRE-PROCESSING METHODS (commander_preproc)
type, extends(commander_base) :: select_frames_commander 
  contains
    procedure :: execute      => exec_select_frames
end type select_frames_commander
type, extends(commander_base) :: boxconvs_commander
  contains
    procedure :: execute      => exec_boxconvs
end type boxconvs_commander
type, extends(commander_base) :: integrate_movies_commander
 contains
   procedure :: execute      => exec_integrate_movies
end type integrate_movies_commander
type, extends(commander_base) :: powerspecs_commander
 contains
   procedure :: execute      => exec_powerspecs
end type powerspecs_commander
type, extends(commander_base) :: unblur_movies_commander
  contains
    procedure :: execute      => exec_unblur_movies
end type unblur_movies_commander
type, extends(commander_base) :: stack_powerspecs_commander
  contains
    procedure :: execute      => exec_stack_powerspecs
end type stack_powerspecs_commander
type, extends(commander_base) :: select_commander
  contains
    procedure :: execute      => exec_select
end type select_commander
type, extends(commander_base) :: extr_ptcls_commander
  contains
    procedure :: execute      => exec_extr_ptcls
end type extr_ptcls_commander
! PRIME2D METHODS (commander_prime2D)
type, extends(commander_base) :: prime2D_init_commander 
  contains
    procedure :: execute      => exec_prime2D_init
end type prime2D_init_commander 
type, extends(commander_base) :: prime2D_commander 
  contains
    procedure :: execute      => exec_prime2D
end type prime2D_commander 
type, extends(commander_base) :: cavgassemble_commander
  contains
    procedure :: execute      => exec_cavgassemble
end type cavgassemble_commander
type, extends(commander_base) :: check2D_conv_commander
  contains
    procedure :: execute      => exec_check2D_conv
end type check2D_conv_commander
type, extends(commander_base) :: rank_cavgs_commander
 contains
   procedure :: execute      => exec_rank_cavgs
end type rank_cavgs_commander

! PRIME3D METHODS (commander_prime3D)
type, extends(commander_base) :: resrange_commander
  contains
    procedure :: execute      => exec_resrange
end type resrange_commander
type, extends(commander_base) :: npeaks_commander
  contains
    procedure :: execute      => exec_npeaks
end type npeaks_commander
type, extends(commander_base) :: nspace_commander
 contains
   procedure :: execute      => exec_nspace
end type nspace_commander
type, extends(commander_base) :: prime3D_init_commander 
  contains
    procedure :: execute      => exec_prime3D_init
end type prime3D_init_commander
type, extends(commander_base) :: multiptcl_init_commander
  contains
    procedure :: execute      => exec_multiptcl_init
end type multiptcl_init_commander
type, extends(commander_base) :: prime3D_commander
  contains
    procedure :: execute      => exec_prime3D
end type prime3D_commander
type, extends(commander_base) :: cont3D_commander
  contains
    procedure :: execute      => exec_cont3D
end type cont3D_commander
type, extends(commander_base) :: check3D_conv_commander
  contains
    procedure :: execute      => exec_check3D_conv
end type check3D_conv_commander

! COMMON-LINES METHODS (commander_comlin)
type, extends(commander_base) :: comlin_smat_commander
  contains
    procedure :: execute      => exec_comlin_smat
end type comlin_smat_commander
type, extends(commander_base) :: symsrch_commander
  contains
    procedure :: execute      => exec_symsrch
end type symsrch_commander

! MASKING METHODS (commander_mask)
type, extends(commander_base) :: mask_commander
 contains
   procedure :: execute      => exec_mask
end type mask_commander
type, extends(commander_base) :: automask2D_commander
  contains
    procedure :: execute      => exec_automask2D
end type automask2D_commander
type, extends(commander_base) :: automask3D_commander
  contains
    procedure :: execute      => exec_automask3D
end type automask3D_commander

! RECONSTRUCTION METHODS (commander_rec)
type, extends(commander_base) :: eo_recvol_commander
  contains
    procedure :: execute      => exec_eo_recvol
end type eo_recvol_commander
type, extends(commander_base) :: eo_volassemble_commander
  contains
    procedure :: execute      => exec_eo_volassemble
end type eo_volassemble_commander
type, extends(commander_base) :: recvol_commander
  contains
    procedure :: execute      => exec_recvol
end type recvol_commander
type, extends(commander_base) :: volassemble_commander
  contains
    procedure :: execute      => exec_volassemble
end type volassemble_commander

! CHECKER METHODS
type, extends(commander_base) :: check_box_commander
  contains
    procedure :: execute      => exec_check_box
end type check_box_commander
type, extends(commander_base) :: check_nptcls_commander
  contains
    procedure :: execute      => exec_check_nptcls
end type check_nptcls_commander
type, extends(commander_base) :: iminfo_commander
 contains
   procedure :: execute      => exec_iminfo
end type iminfo_commander

! VOLOPS METHODS (commander_volops)
type, extends(commander_base) :: cenvol_commander
  contains
    procedure :: execute      => exec_cenvol
end type cenvol_commander
type, extends(commander_base) :: postproc_vol_commander
 contains
   procedure :: execute      => exec_postproc_vol
end type postproc_vol_commander
type, extends(commander_base) :: projvol_commander
 contains
   procedure :: execute      => exec_projvol
end type projvol_commander
type, extends(commander_base) :: volaverager_commander
  contains
    procedure :: execute      => exec_volaverager
end type volaverager_commander
type, extends(commander_base) :: volops_commander
  contains
    procedure :: execute      => exec_volops
end type volops_commander
type, extends(commander_base) :: volume_smat_commander
  contains
    procedure :: execute      => exec_volume_smat
end type volume_smat_commander

! MISCELLANOUS METHODS (commander_misc)
type, extends(commander_base) :: binarise_commander
  contains
    procedure :: execute      => exec_binarise
end type binarise_commander
type, extends(commander_base) :: cluster_smat_commander
  contains
    procedure :: execute      => exec_cluster_smat
end type cluster_smat_commander
type, extends(commander_base) :: converter_commander
  contains
    procedure :: execute      => exec_converter
end type converter_commander
type, extends(commander_base) :: ctfops_commander
  contains
    procedure :: execute      => exec_ctfops
end type ctfops_commander
type, extends(commander_base) :: filter_commander
  contains
    procedure :: execute      => exec_filter
end type filter_commander
type, extends(commander_base) :: gatanmrc2mrc_commander
  contains
    procedure :: execute      => exec_gatanmrc2mrc
end type gatanmrc2mrc_commander
type, extends(commander_base) :: image_smat_commander
 contains
   procedure :: execute      => exec_image_smat
end type image_smat_commander
type, extends(commander_base) :: map2ptcls_commander
 contains
   procedure :: execute      => exec_map2ptcls
end type map2ptcls_commander
type, extends(commander_base) :: norm_commander
  contains
    procedure :: execute      => exec_norm
end type norm_commander
type, extends(commander_base) :: orisops_commander
 contains
   procedure :: execute      => exec_orisops
end type orisops_commander
type, extends(commander_base) :: print_fsc_commander
 contains
   procedure :: execute      => exec_print_fsc
end type print_fsc_commander
type, extends(commander_base) :: res_commander
 contains
   procedure :: execute      => exec_res
end type res_commander
type, extends(commander_base) :: scale_commander
  contains
    procedure :: execute      => exec_scale
end type scale_commander
type, extends(commander_base) :: stackops_commander
  contains
    procedure :: execute      => exec_stackops
end type stackops_commander
type, extends(commander_base) :: tseries_split_commander
  contains
    procedure :: execute      => exec_tseries_split
end type tseries_split_commander

! METHODS FOR PARALLEL PROCESSING (commader_distr)
type, extends(commander_base) :: merge_algndocs_commander
  contains
    procedure :: execute      => exec_merge_algndocs
end type merge_algndocs_commander
type, extends(commander_base) :: merge_similarities_commander
  contains
    procedure :: execute      => exec_merge_similarities
end type merge_similarities_commander
type, extends(commander_base) :: split_pairs_commander
  contains
    procedure :: execute      => exec_split_pairs
end type split_pairs_commander
type, extends(commander_base) :: split_commander
  contains
    procedure :: execute      => exec_split
end type split_commander

!********* METHODS

contains

!     subroutine exec_template( self, cline )
!         class(template_commander), intent(inout) :: self
!         class(cmdline),            intent(inout) :: cline
!         type(params) :: p
!         type(build)  :: b
!         p = params(cline)                     ! parameters generated
!         call b%build_general_tbox(p, cline)   ! general objects built
!
!         ! template specifics
!
!         ! end gracefully
!         call simple_end('**** SIMPLE_TEMPLATE NORMAL STOP ****')
!     end subroutine exec_template
      
    ! SIMULATORS
    
    subroutine exec_simimgs( self, cline )
        use simple_ori,      only: ori
        use simple_rnd,      only: ran3
        use simple_math,     only: deg2rad
        use simple_gridding, only: prep4cgrid
        class(simimgs_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)       :: p
        type(build)        :: b
        type(ori)          :: orientation
        real               :: snr_pink, snr_detector, bfac, bfacerr, dfx, dfy, angast
        integer            :: i, cnt, ntot
        logical, parameter :: debug=.false.
        p = params(cline, .false.)          ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        if( .not. cline%defined('outstk') ) p%outstk = 'simimgs'//p%ext
        if( cline%defined('part') )then
            if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
        else
            if( .not. cline%defined('outfile') ) p%outfile = 'simoris.txt'
        endif
        if( p%box == 0 ) stop 'box=0, something is fishy! Perhaps forgotten to input volume or stack?'
        ! generate orientation/CTF parameters
        if( cline%defined('ndiscrete') )then
            if( p%ndiscrete > 0 )then
                call b%a%rnd_oris_discrete(p%ndiscrete, p%nsym, p%eullims)
            endif
            call b%a%rnd_inpls(p%trs)
        else if( .not. cline%defined('oritab') .and. p%single .eq. 'no' )then
            call b%a%rnd_oris(p%sherr)
        endif
        if( debug )then
            write(*,*) 'CTF parameters used'
            write(*,*) 'kv = ', p%kv
            write(*,*) 'cs = ', p%cs
            write(*,*) 'fraca = ', p%fraca
        endif
        if( .not. b%a%isthere('dfx') )then
            if( p%ctf .ne. 'no' ) call b%a%rnd_ctf(p%defocus, p%deferr, p%astigerr)
        endif
        if( debug ) write(*,'(A)') '>>> DONE GENERATING ORIENTATION/CTF PARAMETERS'
        call b%a%write(p%outfile)
        ! calculate snr:s
        snr_pink = p%snr/0.2
        snr_detector = p%snr/0.8
        if( debug ) write(*,'(A)') '>>> DONE CALCULATING SNR:S'
        ! prepare for image generation
        call b%vol%read(p%vols(1))
        if( debug ) write(*,'(A)') '>>> DID READ VOL'
        call prep4cgrid(b%vol, b%vol_pad, p%msk, wfuns=b%proj%get_wfuns())
        if( debug ) write(*,'(A)') '>>> DONE PREPARING FOR IMAGE GENERATION'
        write(*,'(A)') '>>> GENERATING IMAGES'
        cnt = 0
        ntot = p%top-p%fromp+1
        do i=p%fromp,p%top
            cnt = cnt+1
            call progress(cnt,ntot)
            ! zero images
            b%img_pad = cmplx(0.,0.)
            b%img = 0.
            ! extract ori
            orientation = b%a%get_ori(i)
            ! project vol
            call b%proj%fproject(b%vol_pad, orientation, b%img_pad)
            ! shift
            call b%img_pad%shift(orientation%get('x'),orientation%get('y'))
            ! back FT
            call b%img_pad%bwd_ft
            ! add pink noise
            if( p%snr < 3. ) call b%img_pad%add_gauran(snr_pink)
            call b%img_pad%fwd_ft
            ! apply ctf/bfactor
            bfacerr = ran3()*p%bfacerr
            if( ran3() < 0.5 )then
                bfac = p%bfac-bfacerr
            else
                bfac = p%bfac+bfacerr
            endif
            if( orientation%isthere('dfx') .and. orientation%isthere('dfy') )then
                dfx = orientation%get('dfx')
                dfy = orientation%get('dfy')
                angast = orientation%get('angast')
                if( cline%defined('bfac') )then
                    call b%tfun%apply(b%img_pad,dfx, 'ctf', dfy, angast, bfac=bfac)
                else
                    call b%tfun%apply(b%img_pad,dfx, 'ctf', dfy, angast)
                endif
            else if( orientation%isthere('dfx') )then
                dfx = orientation%get('dfx')
                dfy = dfx
                angast = 0.
                if( cline%defined('bfac') )then
                    call b%tfun%apply(b%img_pad, orientation%get('dfx'), 'ctf', bfac=bfac)
                else
                    call b%tfun%apply(b%img_pad, orientation%get('dfx'), 'ctf')
                endif
            else
                if( cline%defined('bfac') ) call b%img_pad%apply_bfac(bfac)
            endif
            ! add detector noise
            call b%img_pad%bwd_ft
            if( p%snr < 3. ) call b%img_pad%add_gauran(snr_detector)
            if( p%ctf .eq. 'flip' )then
                ! simulate phase-flipped images
                call b%img_pad%fwd_ft
                call b%tfun%apply(b%img_pad, dfx, 'flip', dfy, angast)
                call b%img_pad%bwd_ft
            else if( p%ctf .eq. 'mul' )then
                ! simulate CTF multiplied images
                call b%img_pad%fwd_ft
                call b%tfun%apply(b%img_pad, dfx, 'ctf', dfy, angast)
                call b%img_pad%bwd_ft
            endif
            ! clip
            call b%img_pad%clip(b%img)
            ! write to stack
            if( cline%defined('part') )then
                call b%img%write('simimgs_part'//int2str_pad(p%part,p%numlen)//p%ext, cnt)
            else
                call b%img%write(p%outstk, i)
            endif
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_SIMIMGS NORMAL STOP ****')
    end subroutine exec_simimgs
    
    subroutine exec_simmovie( self, cline )
        use simple_ori,         only: ori
        use simple_math,        only: deg2rad, gen_ptcl_pos
        use simple_image,       only: image
        use simple_rnd,         only: ran3
        use simple_procimgfile, only: stats_imgfile
        class(simmovie_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(params)         :: p
        type(build)          :: b
        type(image)          :: base_image, shifted_base_image
        real                 :: snr_pink, snr_detector, ave, sdev, var, med, fracarea, x, y, sherr, dfx, dfy, deferr, angast
        integer              :: i, ptclarea, mgrapharea, fixed_frame, alloc_stat
        integer, allocatable :: ptcl_positions(:,:)
        real, allocatable    :: shifts(:,:)
        logical              :: here
        logical, parameter   :: debug=.false.
        p = params(cline)                     ! parameters generated
        if( p%box == 0 ) stop 'box=0, something is fishy!'
        call b%build_general_tbox(p, cline)   ! general objects built
        ! set fixed frame
        fixed_frame = nint(real(p%nframes)/2.)
        ! remake the alignment doc
        call b%a%new(1)
        ! check the fractional area occupied by particles & generate particle positions
        ptclarea   = p%box*p%box*p%nptcls
        mgrapharea = p%xdim*p%ydim
        fracarea   = real(ptclarea)/real(mgrapharea)
        write(*,'(a,1x,f7.3)') 'Fraction of area occupied by ptcls:', fracarea
        if( fracarea > 0.55 )then
            write(*,'(A)') 'It is not recommended that more than 55% of the micrograph area is occupied with particles!'
            write(*,'(A)') 'Please, reduce the number of projection images to place!'
            stop
        endif
        write(*,'(a)') '>>> GENERATING PARTICLE POSITIONS'
        ptcl_positions = gen_ptcl_pos(p%nptcls, p%xdim, p%ydim, p%box)
        ! calculate stack background stats
        inquire(file=p%stk, exist=here)
        if( .not. here )then
            write(*,*) 'file does not exis:', p%stk
            stop
        endif
        call stats_imgfile(p%stk, 'background', ave, sdev, var, med, p%msk)
        write(*,'(a)') '>>> STACK BACKGROUND STATISTICS'
        write(*,'(a,1x,f7.4)') 'ave: ', ave
        write(*,'(a,1x,f7.4)') 'sdev:', sdev
        write(*,'(a,1x,f7.4)') 'med: ', med
        ! make a base image by inserting the projections at ptcl_positions and setting the background to med
        call base_image%new([p%xdim,p%ydim,1], p%smpd, backgr=med)
        do i=1,p%nptcls
            call b%img%read(p%stk, i)
            call b%img%insert(ptcl_positions(i,:), base_image)
        end do
        if( p%vis .eq. 'yes' ) call base_image%vis
        if( debug ) write(*,'(a)') 'inserted projections'
        ! calculate snr:s
        snr_pink     = p%snr/0.2
        snr_detector = p%snr/0.8
        if( debug ) write(*,'(a)') 'calculated SNR:s'
        ! set CTF parameters
        deferr = ran3()*0.2
        if( ran3() < 0.5 )then
            dfx = p%defocus-deferr
        else
            dfx = p%defocus+deferr
        endif
        deferr = ran3()*0.2
        if( ran3() < 0.5 )then
            dfy = p%defocus-deferr
        else
            dfy = p%defocus+deferr
        endif
        angast = ran3()*360.
        if( debug ) write(*,'(a)') 'did set CTF parameters'
        if( debug ) write(*,'(a)') 'initialized shifts'
        ! generate shifts
        allocate( shifts(p%nframes,2), stat=alloc_stat )
        call alloc_err('In: simple_simmovie; shifts', alloc_stat)
        x = 0.
        y = 0.
        do i=1,p%nframes
            ! generate random shifts
            sherr = ran3()*p%trs
            if( ran3() > 0.5 )then
                x = x+sherr
            else
                x = x-sherr
            endif
            sherr = ran3()*p%trs
            if( ran3() > 0.5 )then
                y = y+sherr
            else
                y = y-sherr
            endif
            shifts(i,1) = x
            shifts(i,2) = y
        end do
        ! put the central frame in the series at x,y = (0,0) & fill-up b%a
        do i=1,p%nframes
            shifts(i,:) = shifts(i,:)-shifts(fixed_frame,:)
            call b%a%set(1, 'x'//int2str(i), shifts(i,1))
            call b%a%set(1, 'y'//int2str(i), shifts(i,2))
        end do
        ! make and open a stack for the movie frames
        if( debug ) write(*,'(a)') 'made stack for output movie frames'
        write(*,'(a)') '>>> GENERATING MOVIE FRAMES'
        do i=1,p%nframes
            call progress(i,p%nframes)
            ! shift base image
            call base_image%fwd_ft
            if( i > 1 ) then
                call base_image%shift(shifts(i,1),shifts(i,2),imgout=shifted_base_image)
            else
                shifted_base_image = base_image
            endif
            call shifted_base_image%bwd_ft
            ! add pink noise
            call shifted_base_image%add_gauran(snr_pink)
            ! multiply with CTF
            call shifted_base_image%fwd_ft
            if( p%neg .eq. 'yes' )then
                call b%tfun%apply(shifted_base_image, dfx, 'neg', dfy, angast, p%bfac)
            else
                call b%tfun%apply(shifted_base_image, dfx, 'ctf', dfy, angast, p%bfac)
            endif
            call shifted_base_image%bwd_ft
            ! add the detector noise
            call shifted_base_image%add_gauran(snr_detector)
            if( p%vis .eq. 'yes' ) call shifted_base_image%vis
            call shifted_base_image%write('simmovie'//p%ext, i)
            ! set orientation parameters in object
            call b%a%set(1, 'dfx',    dfx)
            call b%a%set(1, 'dfy',    dfy)
            call b%a%set(1, 'angast', angast)
        end do
        if( debug ) write(*,'(a)') 'generated movie'
        ! generate the optimal average
        base_image = 0.
        write(*,'(a)') '>>> GENERATING OPTIMAL AVERAGE'
        do i=1,p%nframes
            call progress(i,p%nframes)
            call shifted_base_image%read('simmovie'//p%ext, i)
            x = b%a%get(1, 'x'//int2str(i))
            y = b%a%get(1, 'y'//int2str(i))
            call shifted_base_image%shift(-x,-y)
            call base_image%add(shifted_base_image)
        end do
        if( debug ) write(*,'(a,1x,f7.4)') 'constant 4 division:', real(p%nframes)
        call base_image%div(real(p%nframes))
        if( debug ) write(*,'(a)') 'generated optimal average'
        call base_image%write('optimal_movie_average'//p%ext, 1)
        if( p%vis .eq. 'yes' ) call base_image%vis
        ! output orientations
        call b%a%write('simmovie_params.txt')
        ! end gracefully
        call simple_end('**** SIMPLE_SIMMOVIE NORMAL STOP ****')
    end subroutine exec_simmovie
    
    subroutine exec_simsubtomo( self, cline )
        use simple_image, only: image
        use simple_ori,   only: ori
        class(simsubtomo_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        type(image)  :: vol_rot
        type(ori)    :: o
        integer      :: iptcl, numlen
        p = params(cline, checkdistr=.false.) ! constants & derived constants produced, mode=2
        call b%build_general_tbox(p, cline)   ! general objects built
        call vol_rot%new([p%box,p%box,p%box], p%smpd)
        call b%vol%new([p%box,p%box,p%box],   p%smpd)
        call b%vol%read(p%vols(1))
        call b%vol%add_gauran(p%snr)
        call o%new
        numlen = len(int2str(p%nptcls))
        do iptcl=1,p%nptcls
            call o%rnd_ori
            vol_rot = b%proj%rotvol(b%vol, o, p)
            call vol_rot%write('subtomo'//int2str_pad(iptcl,numlen)//p%ext)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_SIMSUBTOMO NORMAL STOP ****')
    end subroutine exec_simsubtomo

    ! PRE-PROCESSING METHODS
    
    subroutine exec_select_frames( self, cline )
        use simple_imgfile, only: imgfile
        use simple_image,   only: image
        class(select_frames_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(params)                       :: p
        type(build)                        :: b
        integer                            :: nmovies, nframes, frame, numlen, ldim(3)
        integer                            :: fromto(2), movie, cnt, cnt2, ntot, lfoo(3), ifoo
        character(len=STDLEN), allocatable :: movienames(:)
        character(len=:), allocatable      :: new_name
        type(image)                        :: img_frame
        logical, parameter                 :: debug = .false.
        p = params(cline, checkdistr=.false.)            ! constants & derived constants produced
        call b%build_general_tbox(p,cline,do3d=.false.) ! general objects built
        call read_filetable(p%filetab, movienames)
        nmovies = size(movienames)
        ! find ldim and numlen (length of number string)
        if( cline%defined('startit') )then
            call find_ldim_nptcls(movienames(p%startit), ldim, ifoo)
        endif
        if( debug ) write(*,*) 'logical dimension: ', ldim
        ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
        numlen = len(int2str(nmovies))
        if( debug ) write(*,*) 'length of number string: ', numlen
        ! determine loop range
        if( cline%defined('part') )then
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = p%fromp
                fromto(2) = p%top
                ntot = fromto(2)-fromto(1)+1
            else
                stop 'fromp & top args need to be defined in parallel execution; simple_select_frames'
            endif
        else
            fromto(1) = 1
            if( cline%defined('startit') ) fromto(1) = p%startit
            fromto(2) = nmovies
            ntot      = nmovies
        endif
        if( debug ) write(*,*) 'fromto: ', fromto(1), fromto(2)
        call img_frame%new([ldim(1),ldim(2),1], p%smpd)
        ! loop over exposures (movies)
        cnt2 = 0
        do movie=fromto(1),fromto(2)
            if( .not. file_exists(movienames(movie)) )then
                write(*,*) 'inputted movie stack does not exist: ', trim(adjustl(movienames(movie)))
            endif
            cnt2 = cnt2+1
            ! get number of frames from stack
            call find_ldim_nptcls(movienames(movie), lfoo, nframes)
            if( debug ) write(*,*) 'number of frames: ', nframes
            ! create new name
            allocate(new_name, source=trim(adjustl(p%fbody))//int2str_pad(movie, numlen)//p%ext)
            cnt = 0
            do frame=p%fromp,p%top
                cnt = cnt+1
                call img_frame%read(movienames(movie),frame)
                call img_frame%write(new_name,cnt)
            end do
            deallocate(new_name)
            write(*,'(f4.0,1x,a)') 100.*(real(cnt2)/real(ntot)), 'percent of the movies processed'
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_SELECT_FRAMES NORMAL STOP ****')
    end subroutine exec_select_frames
    
    subroutine exec_boxconvs( self, cline )
        use simple_image, only: image
        class(boxconvs_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(params)                       :: p
        type(build)                        :: b
        type(image)                        :: tmp
        character(len=STDLEN), allocatable :: imgnames(:)
        integer                            :: iimg, nimgs, ldim(3), iimg_start, iimg_stop, ifoo
        logical                            :: debug=.false.
        if( cline%defined('stk') .and. cline%defined('filetab') )then
            stop 'stk and filetab cannot both be defined; input either or!'
        endif
        if( .not. cline%defined('stk') .and. .not. cline%defined('filetab') )then
            stop 'either stk or filetab need to be defined!'
        endif
        p = params(cline, checkdistr=.false.)               ! parameters generated
        call b%build_general_tbox( p, cline, do3d=.false. ) ! general objects built
        ! do the work
        if( cline%defined('stk') )then
            call b%img%new(p%ldim, p%smpd) ! img re-generated (to account for possible non-square)
            tmp = 0.0
            do iimg=1,p%nptcls
                call b%img%read(p%stk, iimg)
                tmp = b%img%boxconv(p%boxconvsz)
                call tmp%write(trim(adjustl(p%fbody))//p%ext, iimg)
                call progress(iimg, p%nptcls)
            end do
        else
            call read_filetable(p%filetab, imgnames)
            nimgs = size(imgnames)
            if( debug ) write(*,*) 'read the img filenames'
            ! get logical dimension of micrographs
            call find_ldim_nptcls(imgnames(1), ldim, ifoo)
            ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
            if( debug ) write(*,*) 'logical dimension: ', ldim
            call b%img%new(ldim, p%smpd) ! img re-generated (to account for possible non-square)
            ! determine loop range
            iimg_start = 1
            if( cline%defined('startit') ) iimg_start = p%startit
            iimg_stop  = nimgs
            if( debug ) write(*,*) 'fromto: ', iimg_start, iimg_stop
            ! do it
            tmp = 0.0
            do iimg=iimg_start,iimg_stop
                if( .not. file_exists(imgnames(iimg)) )then
                    write(*,*) 'inputted img file does not exist: ', trim(adjustl(imgnames(iimg)))
                endif
                call b%img%read(imgnames(iimg), 1)
                tmp = b%img%boxconv(p%boxconvsz)
                call tmp%write(trim(adjustl(p%fbody))//p%ext, iimg)
                call progress(iimg, nimgs)
            end do
            call tmp%kill
            deallocate(imgnames)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_BOXCONVS NORMAL STOP ****')
    end subroutine exec_boxconvs
    
    subroutine exec_integrate_movies(self,cline)
        use simple_imgfile, only: imgfile
        use simple_image,   only: image
        class(integrate_movies_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(params)                       :: p
        type(build)                        :: b
        integer                            :: nmovies, nframes, frame, alloc_stat, lfoo(3)
        integer                            :: numlen, ldim(3), fromto(2), movie, ifoo
        character(len=STDLEN), allocatable :: movienames(:)
        character(len=:), allocatable      :: cpcmd, new_name
        real                               :: x, y
        type(image), allocatable           :: img_frames(:)
        type(image)                        :: img_sum, pspec
        logical, parameter                 :: debug = .false.
        p = params(cline,checkdistr=.false.) ! constants & derived constants produced
        call b%build_general_tbox(p,cline,do3d=.false.)
        call read_filetable(p%filetab, movienames)
        nmovies = size(movienames)
        if( debug ) write(*,*) 'read the movie filenames'
        ! find ldim and numlen (length of number string)
        call find_ldim_nptcls(movienames(1), ldim, ifoo)
        if( debug ) write(*,*) 'logical dimension: ', ldim
        ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
        numlen  = len(int2str(nmovies))
        if( debug ) write(*,*) 'length of number string: ', numlen
        ! determine loop range
        if( cline%defined('part') )then
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = p%fromp
                fromto(2) = p%top
            else
                stop 'fromp & top args need to be defined in parallel execution; simple_integrate_movies'
            endif
        else
            fromto(1) = 1
            fromto(2) = nmovies
        endif
        if( debug ) write(*,*) 'fromto: ', fromto(1), fromto(2)
        ! create sum
        call img_sum%new([ldim(1),ldim(2),1], p%smpd)
        ! loop over exposures (movies)
        do movie=fromto(1),fromto(2)
            if( .not. file_exists(movienames(movie)) )&
            & write(*,*) 'inputted movie stack does not exist: ', trim(adjustl(movienames(movie)))
            ! get number of frames from stack
            call find_ldim_nptcls(movienames(movie), lfoo, nframes)
            if( debug ) write(*,*) 'number of frames: ', nframes
            ! create frames & read
            allocate(img_frames(nframes), stat=alloc_stat)
            call alloc_err('In: simple_integrate_movies', alloc_stat)
            img_sum = 0.
            do frame=1,nframes
                call img_frames(frame)%new([ldim(1),ldim(2),1], p%smpd)
                call img_frames(frame)%read(movienames(movie),frame)
                if( cline%defined('oritab') )then
                    ! shift frame according to global shift (drift correction)
                    if( b%a%isthere(movie, 'x'//int2str(frame)) .and. b%a%isthere(movie, 'y'//int2str(frame)) )then
                        call img_frames(frame)%fwd_ft
                        x = b%a%get(movie, 'x'//int2str(frame))
                        y = b%a%get(movie, 'y'//int2str(frame))
                        if( debug ) print *, 'shifting frame: ', x, y
                        call img_frames(frame)%shift(-x,-y)
                        call img_frames(frame)%bwd_ft
                        call img_sum%add(img_frames(frame))
                    else
                        stop 'no shift parameters available for alignment in oritab'
                    endif
                else
                    call img_sum%add(img_frames(frame))
                endif
            end do
            ! rename movie
            allocate(new_name, source=trim(adjustl(p%fbody))//int2str_pad(movie, numlen)//p%ext)
            if( .not. file_exists(new_name))then
                allocate(cpcmd, source='cp '//trim(adjustl(movienames(movie)))//' ./'//new_name)
                call system(cpcmd)
                deallocate(cpcmd)
            endif
            ! destroy objects and deallocate
            do frame=1,nframes
                call img_frames(frame)%kill
            end do
            deallocate(new_name,img_frames)
            call img_sum%write(trim(adjustl(p%fbody))//'_intg'//int2str_pad(movie, numlen)//p%ext)
            pspec = img_sum%mic2spec(p%pspecsz, trim(adjustl(p%speckind)))
            call pspec%write(trim(adjustl(p%fbody))//'_pspec'//int2str_pad(movie, numlen)//p%ext)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_INTEGRATE_MOVIES NORMAL STOP ****')
    end subroutine exec_integrate_movies
    
    subroutine exec_powerspecs( self, cline )
        use simple_imgfile, only: imgfile
        use simple_image,   only: image
        class(powerspecs_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params)                       :: p
        type(build)                        :: b
        type(image)                        :: powspec, tmp, mask
        character(len=STDLEN), allocatable :: imgnames(:)
        integer                            :: iimg, nimgs, ldim(3), iimg_start, iimg_stop, ifoo
        logical                            :: debug=.false.
        if( cline%defined('stk') .and. cline%defined('filetab') )then
            stop 'stk and filetab cannot both be defined; input either or!'
        endif
        if( .not. cline%defined('stk') .and. .not. cline%defined('filetab') )then
            stop 'either stk or filetab need to be defined!'
        endif
        p = params(cline, checkdistr=.false.)           ! constants & derived constants produced
        call b%build_general_tbox(p,cline,do3d=.false.) ! general toolbox built
        ! create mask
        call tmp%new([p%clip,p%clip,1], p%smpd)
        tmp = cmplx(1.,0.)
        call tmp%bp(0.,p%lp,0.)
        call tmp%ft2img('real', mask)
        call mask%write('resolution_mask.mrc', 1)
        ! do the work
        if( cline%defined('stk') )then
            call b%img%new(p%ldim, p%smpd) ! img re-generated (to account for possible non-square)
            tmp = 0.0
            do iimg=1,p%nptcls
                call b%img%read(p%stk, iimg)
                powspec = b%img%mic2spec(p%pspecsz, trim(adjustl(p%speckind)))
                call powspec%clip(tmp)
                call tmp%write(trim(adjustl(p%fbody))//p%ext, iimg)
                call progress(iimg, p%nptcls)
            end do
        else
            call read_filetable(p%filetab, imgnames)
            nimgs = size(imgnames)
            if( debug ) write(*,*) 'read the img filenames'
            ! get logical dimension of micrographs
            call find_ldim_nptcls(imgnames(1), ldim, ifoo)
            ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
            if( debug ) write(*,*) 'logical dimension: ', ldim
            call b%img%new(ldim, p%smpd) ! img re-generated (to account for possible non-square)
            ! determine loop range
            iimg_start = 1
            if( cline%defined('startit') ) iimg_start = p%startit
            iimg_stop  = nimgs
            if( debug ) write(*,*) 'fromto: ', iimg_start, iimg_stop
            ! do it
            tmp = 0.0
            do iimg=iimg_start,iimg_stop
                if( .not. file_exists(imgnames(iimg)) )then
                    write(*,*) 'inputted img file does not exist: ', trim(adjustl(imgnames(iimg)))
                endif
                call b%img%read(imgnames(iimg), 1)
                powspec = b%img%mic2spec(p%pspecsz, trim(adjustl(p%speckind)))
                call powspec%clip(tmp)
                call tmp%write(trim(adjustl(p%fbody))//p%ext, iimg)
                call progress(iimg, nimgs)
            end do
        endif
        call simple_end('**** SIMPLE_POWERSPECS NORMAL STOP ****')
        ! end gracefully
        call simple_end('**** SIMPLE_INTEGRATE_MOVIES NORMAL STOP ****')
    end subroutine exec_powerspecs
    
    subroutine exec_unblur_movies( self, cline )
        use simple_unblur   ! singleton
        use simple_imgfile, only: imgfile
        use simple_oris,    only: oris
        use simple_image,   only: image
        class(unblur_movies_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(params)                       :: p
        type(build)                        :: b
        type(image)                        :: movie_sum, pspec_sum, pspec_corrected, movie_sum_ctf
        type(image)                        :: movie_sum_corrected, pspec_half_n_half
        integer                            :: nmovies, imovie, imovie_start, imovie_stop, file_stat
        integer                            :: funit_movies, cnt, numlen,  ntot, alloc_stat, fnr
        character(len=:), allocatable      :: cpcmd, new_name
        character(len=STDLEN), allocatable :: movienames(:)
        real                               :: corr
        logical                            :: debug=.true.
        ! CUDA err variable for the return function calls
        integer                            :: err
        call timestamp()
        ! call start_Alltimers_cpu()
        ! starting the cuda environment
        call simple_cuda_init(err)
        if( err .ne. RC_SUCCESS ) write(*,*) 'cublas init failed'
        p = params(cline, checkdistr=.false.)            ! constants & derived constants produced
        call b%build_general_tbox(p,cline,do3d=.false.) ! general objects built
        ! read movienames
        nmovies = nlines(p%filetab)
        if( debug ) write(*,*) 'nmovies: ', nmovies
        allocate( movienames(nmovies), stat=alloc_stat )
        call alloc_err('In: simple_integrate_movies', alloc_stat)
        funit_movies = get_fileunit()
        open(unit=funit_movies, status='old', file=p%filetab)
        do imovie=1,nmovies
            read(funit_movies,'(a256)') movienames(imovie)
        end do
        close(funit_movies)
        if( debug ) write(*,*) 'read the movie filenames'
        ! determine loop range
        if( cline%defined('part') )then
            if( cline%defined('fromp') .and. cline%defined('top') )then
                imovie_start = p%fromp
                imovie_stop  = p%top
            else
                stop 'fromp & top args need to be defined in parallel execution; simple_unblur_movies'
            endif
        else
            imovie_start = 1
            if( cline%defined('startit') ) imovie_start = p%startit
            imovie_stop  = nmovies
        endif
        if( debug ) write(*,*) 'fromto: ', imovie_start, imovie_stop
        ntot = imovie_stop-imovie_start+1
        if( cline%defined('numlen') )then
            numlen = p%numlen
        else
            numlen = len(int2str(nmovies))
        endif
        if( debug ) write(*,*) 'length of number string: ', numlen
        ! create output orientations
        call b%a%new(ntot)
        ! align
        cnt = 0
        do imovie=imovie_start,imovie_stop
            if( .not. file_exists(movienames(imovie)) )then
                write(*,*) 'inputted movie stack does not exist: ', trim(adjustl(movienames(imovie)))
            endif
            write(*,'(a,1x,i5)') '>>> PROCESSING MOVIE:', imovie
            cnt = cnt+1
            call unblur_movie(movienames(imovie), p, cnt, b%a, corr, movie_sum, movie_sum_corrected, movie_sum_ctf)
            if( debug ) print *, 'ldim of output movie_sum:           ', movie_sum%get_ldim()
            if( debug ) print *, 'ldim of output movie_sum_corrected: ', movie_sum_corrected%get_ldim()
            if( debug ) print *, 'ldim of output movie_sum_ctf      : ', movie_sum_ctf%get_ldim()
            if( cline%defined('fbody') )then
                call movie_sum_corrected%write(trim(adjustl(p%fbody))//'_intg'//int2str_pad(imovie, numlen)//p%ext)
                call movie_sum_ctf%write(trim(adjustl(p%fbody))//'_ctf'//int2str_pad(imovie, numlen)//p%ext)
            else
                call movie_sum_corrected%write(int2str_pad(imovie, numlen)//'_intg'//p%ext)
                call movie_sum_ctf%write(int2str_pad(imovie, numlen)//'_ctf'//p%ext)
            endif
            pspec_sum         = movie_sum%mic2spec(p%pspecsz, trim(adjustl(p%speckind)))
            pspec_corrected   = movie_sum_corrected%mic2spec(p%pspecsz, trim(adjustl(p%speckind)))
            pspec_half_n_half = pspec_sum%before_after(pspec_corrected)
            if( debug ) print *, 'ldim of pspec_sum:         ', pspec_sum%get_ldim()
            if( debug ) print *, 'ldim of pspec_corrected:   ', pspec_corrected%get_ldim()
            if( debug ) print *, 'ldim of pspec_half_n_half: ', pspec_half_n_half%get_ldim()
            if( cline%defined('fbody') )then
                call pspec_half_n_half%write(trim(adjustl(p%fbody))//'_pspec'//int2str_pad(imovie, numlen)//p%ext)
            else
                call pspec_half_n_half%write(int2str_pad(imovie, numlen)//'_pspec'//p%ext)
            endif
            call movie_sum%kill
            call movie_sum_corrected%kill
            call movie_sum_ctf%kill
            call pspec_sum%kill
            call pspec_corrected%kill
            call pspec_half_n_half%kill
            ! RENAME MOVIE
            if( cline%defined('fbody') )then
                allocate(new_name, source=trim(adjustl(p%fbody))//int2str_pad(imovie, numlen)//p%ext)
            else
                allocate(new_name, source=int2str_pad(imovie, numlen)//p%ext)
            endif
            if( .not. file_exists(new_name))then
                allocate(cpcmd, source='cp '//trim(adjustl(movienames(imovie)))//' ./'//new_name)
                call system(cpcmd)
                deallocate(cpcmd)
            endif
            deallocate(new_name)
            if( .not. cline%defined('part') ) call b%a%write(imovie, 'unblur_movies_params.txt')
            write(*,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the movies processed'
        end do
        if( cline%defined('part') )then
            call b%a%write('unblur_movies_params'//'_part'//int2str_pad(p%part,p%numlen)//'.txt')
            fnr = get_fileunit()
            open(unit=fnr, FILE='JOB_FINISHED_'//int2str_pad(p%part,p%numlen), STATUS='REPLACE', action='WRITE', iostat=file_stat)
            call fopen_err( 'In: simple_unblur_movies', file_stat )
            write(fnr,'(A)') '**** SIMPLE_COMLIN_SMAT NORMAL STOP ****'
            close(fnr)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_UNBLUR_MOVIES NORMAL STOP ****')
        !*******************************************************************************
        !    Environment shutdown
        !
        !*******************************************************************************
        !shutting down the environment
        call simple_cuda_shutdown()
        !shutting down the timers
        ! call stop_Alltimers_cpu()
    end subroutine exec_unblur_movies
    
    subroutine exec_stack_powerspecs( self, cline )
        use simple_imgfile, only: imgfile
        use simple_image,   only: image
        class(stack_powerspecs_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(params)                       :: p
        type(build)                        :: b
        integer                            :: nspecs, ldim(3), ispec, ifoo
        character(len=STDLEN), allocatable :: specnames(:)
        type(image)                        :: mask, tmp
        real                               :: mm(2)
        logical, parameter                 :: debug = .false.    
        p = params(cline,checkdistr=.false.)             ! constants & derived constants produced
        call b%build_general_tbox(p,cline,do3d=.false.) ! general stuff built
        call read_filetable(p%filetab, specnames)
        nspecs = size(specnames)
        if( debug ) write(*,*) 'read the spec filenames'
        call find_ldim_nptcls(specnames(1),ldim,ifoo)
        ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
        if( debug ) write(*,*) 'logical dimension: ', ldim
        p%box   = ldim(1)
        ! create mask
        call tmp%new([p%clip,p%clip,1], p%smpd)
        tmp = cmplx(1.,0.)
        call tmp%bp(0.,p%lp,0.)
        call tmp%ft2img('real', mask)
        call mask%write('resolution_mask.mrc', 1)
        ! prepare b%img and tmp for reading
        call b%img%new([p%box,p%box,1], p%smpd)
        tmp = 0.0
        ! loop over spectra
        do ispec=1,nspecs
            if( .not. file_exists(specnames(ispec)) )then
                write(*,*) 'inputted spec file does not exist: ', trim(adjustl(specnames(ispec)))
            endif
            call b%img%read(specnames(ispec))
            call b%img%clip(tmp)  
            mm = tmp%minmax()
            if( debug ) print *, 'min/max: ', mm(1), mm(2)
            call tmp%write(p%outstk, ispec)
            call progress(ispec, nspecs)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_STACK_POWERSPECS NORMAL STOP ****')
    end subroutine exec_stack_powerspecs

    subroutine exec_select( self, cline )
        use simple_image,  only: image
        use simple_corrmat ! singleton
        class(select_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params)                       :: p
        type(build)                        :: b
        type(image), allocatable           :: imgs_sel(:), imgs_all(:)
        type(image)                        :: stk3_img
        character(len=STDLEN), allocatable :: imgnames(:)
        integer                            :: iimg, isel, nall, nsel, loc(1), ios, funit, ldim(3), ifoo, lfoo(3)
        integer, allocatable               :: selected(:)
        real, allocatable                  :: correlations(:,:)
        logical, parameter                 :: debug=.false.
        ! error check
        if( cline%defined('stk3') .or. cline%defined('filetab') )then
            ! all good
        else
            stop 'Need either stk3 or filetab are part of the command line!'
        endif
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        ! find number of selected cavgs
        call find_ldim_nptcls(p%stk2, lfoo, nsel)
        ! find number of original cavgs
        call find_ldim_nptcls(p%stk, lfoo, nall)
        ! read images
        allocate(imgs_sel(nsel), imgs_all(nall))
        do isel=1,nsel
            call imgs_sel(isel)%new([p%box,p%box,1], p%smpd)
            call imgs_sel(isel)%read(p%stk2, isel)
        end do
        do iimg=1,nall
            call imgs_all(iimg)%new([p%box,p%box,1], p%smpd)
            call imgs_all(iimg)%read(p%stk, iimg)
        end do
        write(*,'(a)') '>>> CALCULATING CORRELATIONS'
        call calc_cartesian_corrmat(imgs_sel, imgs_all, correlations)
        ! find selected
        allocate(selected(nsel))
        do isel=1,nsel
            loc = maxloc(correlations(isel,:))
            selected(isel) = loc(1)
            if( debug ) print *, 'selected: ', loc(1), ' with corr: ', correlations(isel,loc(1))
        end do
        if( cline%defined('filetab') )then
            ! read filetable
            call read_filetable(p%filetab, imgnames)
            if( size(imgnames) /= nall ) stop 'nr of entries in filetab and stk not consistent'
            funit = get_fileunit()
            open(unit=funit, file=p%outfile, iostat=ios, status="replace", action="write", access="sequential")
            if( ios /= 0 )then
                write(*,*) "Error opening file name", trim(adjustl(p%outfile))
                stop
            endif
            call execute_command_line('mkdir '//trim(adjustl(p%dir)))
            ! write outoput & move files
            do isel=1,nsel
                write(funit,'(a)') trim(adjustl(imgnames(selected(isel))))
                call execute_command_line('mv '//trim(adjustl(imgnames(selected(isel))))//' '//trim(adjustl(p%dir)))
            end do
            close(funit)
            deallocate(imgnames)
        endif
        if( cline%defined('stk3') )then
            call find_ldim_nptcls(p%stk3, ldim, ifoo)
            ldim(3) = 1
            call stk3_img%new(ldim,1.)
            do isel=1,nsel
                call stk3_img%read(p%stk3, selected(isel))
                call stk3_img%write(p%outstk, isel)
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_SELECT NORMAL STOP ****')
    end subroutine exec_select
    
    subroutine exec_extr_ptcls( self, cline )
        use simple_nrtxtfile, only: nrtxtfile
        use simple_imgfile,   only: imgfile
        use simple_image,     only: image
        use simple_math,      only: euclid, hpsort, median
        use simple_oris,      only: oris
        use simple_stat,      only: moment
        class(extr_ptcls_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params), target               :: p
        type(build),  target               :: b
        type(nrtxtfile)                    :: boxfile
        type(image)                        :: img_frame, mskimg
        type(oris)                         :: outoris
        integer                            :: nmovies, nboxfiles, funit_movies, funit_box, nframes, frame, pind
        integer                            :: i, j, k, alloc_stat, ldim(3), box_current, movie, ndatlines, nptcls
        integer                            :: npix, npix_backgr, npix_tot, cnt, niter, nnn, fromto(2), orig_box
        integer                            :: movie_ind, numlen, ntot, lfoo(3), ifoo
        character(len=STDLEN)              :: framestack, sumstack
        character(len=STDLEN), allocatable :: movienames(:), boxfilenames(:)
        real, allocatable                  :: boxdata(:,:), nndists(:), noise_pixels(:), pixels(:)
        integer, allocatable               :: nninds(:), pinds(:)
        real                               :: x, y, dfx, dfy, angast, ctfres, med, ave,sdev, var
        logical                            :: err
        logical, parameter                 :: debug = .false.
        p = params(cline, checkdistr=.false.) ! parameters generated
        ! check file inout existence
        if( .not. file_exists(p%filetab) ) stop 'inputted filetab does not exist in cwd'
        if( .not. file_exists(p%boxtab) )  stop 'inputted boxtab does not exist in cwd'
        nmovies = nlines(p%filetab)
        if( debug ) write(*,*) 'nmovies: ', nmovies
        nboxfiles = nlines(p%boxtab)
        if( debug ) write(*,*) 'nboxfiles: ', nboxfiles
        if( nmovies /= nboxfiles ) stop 'number of entries in inputted files do not match!'
        funit_movies = get_fileunit()
        open(unit=funit_movies, status='old', file=p%filetab)
        funit_box = get_fileunit()
        open(unit=funit_box, status='old', file=p%boxtab)
        allocate( movienames(nmovies), boxfilenames(nmovies), stat=alloc_stat )
        call alloc_err('In: simple_extr_ptcls; boxdata etc., 1', alloc_stat)

        ! remove output file
        call del_txtfile('extr_ptcls_params.txt')

        ! read the filenames
        do movie=1,nmovies
            read(funit_movies,'(a256)') movienames(movie)
            read(funit_box,'(a256)') boxfilenames(movie)
        end do
        close(funit_movies)
        close(funit_box)
        if( debug ) write(*,*) 'read the filenames'

        ! determine loop range
        fromto(1) = 1
        fromto(2) = nmovies
        ntot = fromto(2)-fromto(1)+1
        if( debug ) write(*,*) 'fromto: ', fromto(1), fromto(2)

        ! count the number of particles & find ldim
        nptcls = 0
        do movie=fromto(1),fromto(2)
            if( file_exists(boxfilenames(movie)) )then 
                call boxfile%new(boxfilenames(movie), 1)
                ndatlines = boxfile%get_ndatalines()
                nptcls = nptcls+ndatlines
                call boxfile%kill
            else
                write(*,*) 'WARNING! The inputted boxfile (below) does not exist'
                write(*,*) trim(boxfilenames(movie))
            endif
        end do
        call find_ldim_nptcls(movienames(1), ldim, ifoo)
        if( debug ) write(*,*) 'number of particles: ', nptcls

        ! create frame
        call img_frame%new([ldim(1),ldim(2),1], p%smpd)

        ! initialize
        call outoris%new(nptcls)
        pind = 0
        if( debug ) write(*,*) 'made outoris'

        ! loop over exposures (movies)
        niter = 0
        do movie=fromto(1),fromto(2)
            
            ! show progress
            if( niter > 1 )then
                call progress(niter,ntot)
            endif
            
            ! get movie index (parsing the number string from the filename)
            call fname2ind(trim(adjustl(movienames(movie))), movie_ind)

            ! process boxfile
            ndatlines = 0
            if( file_exists(boxfilenames(movie)) )then
                call boxfile%new(boxfilenames(movie), 1)
                ndatlines = boxfile%get_ndatalines()        
            endif
            if( ndatlines == 0 ) cycle
            
            ! update iteration counter (has to be after the cycle statements or suffer bug!!!)
            niter = niter+1
            
            ! read box data
            allocate( boxdata(ndatlines,boxfile%get_nrecs_per_line()),&
            nninds(ndatlines), nndists(ndatlines), pinds(ndatlines), stat=alloc_stat )
            call alloc_err('In: simple_extr_ptcls; boxdata etc., 2', alloc_stat)
            do j=1,ndatlines
                call boxfile%readNextDataLine(boxdata(j,:))
                orig_box = nint(boxdata(j,3))
                if( nint(boxdata(j,3)) /= nint(boxdata(j,4)) )then
                    stop 'Only square windows are currently allowed!'
                endif
            end do
            
            ! create particle index list and set movie index
            if( .not. cline%defined('box') ) p%box = nint(boxdata(1,3)) ! use the box size from the box file
            do j=1,ndatlines
                if( box_inside(ldim, nint(boxdata(j,1:2)), p%box) )then
                    pind     = pind+1
                    pinds(j) = pind
                    call outoris%set(pinds(j), 'movie', real(movie_ind))
                else
                    pinds(j) = 0
                endif
            end do
            
            ! create nearest neighbor structure
            if( p%mapnn .eq. 'yes' )then
                if( ndatlines > 1 )then
                    do i=1,ndatlines
                        do j=1,ndatlines
                            if( i == j )then
                                nndists(j) = huge(x)
                            else
                                nndists(j) = euclid(boxdata(i,1:2), boxdata(j,1:2))
                            endif
                        end do
                        nninds = pinds
                        call hpsort(ndatlines, nndists, nninds)
                        nnn = min(20,ndatlines-1)
                        do j=1,nnn
                            if( pinds(i) > 0 .and. nndists(j) > 0 )then
                                call outoris%set(pinds(i), 'nnn',             real(nnn))        ! number of nearest neighbors
                                call outoris%set(pinds(i), 'nn'//int2str(j),  real(nninds(j)))  ! nearest neighbor particle index
                                call outoris%set(pinds(i), 'dnn'//int2str(j), real(nndists(j))) ! nearest neighbor distance
                            endif
                        end do
                    end do
                endif
            endif
            
            ! check box parsing
            if( .not. cline%defined('box') )then
                if( niter == 1 )then
                    p%box = nint(boxdata(1,3))   ! use the box size from the box file
                else
                    box_current = nint(boxdata(1,3))
                    if( box_current /= p%box )then
                        write(*,*) 'box_current: ', box_current, 'box in params: ', p%box
                        stop 'inconsistent box sizes in box files'
                    endif
                endif
                if( p%box == 0 )then
                    write(*,*) 'ERROR, box cannot be zero!'
                    stop
                endif
            endif
            if( debug ) write(*,*) 'did check box parsing'
            
            ! get number of frames from stack
            call find_ldim_nptcls(movienames(movie), lfoo, nframes )
            numlen  = len(int2str(nframes))
            if( debug ) write(*,*) 'number of frames: ', nframes
            
            ! build general objects
            if( niter == 1 )then
                call b%build_general_tbox(p,cline,do3d=.false.)
                call mskimg%disc([p%box,p%box,1], p%smpd, p%msk, npix)
                call mskimg%bin_inv
                npix_backgr = p%box*p%box-npix
            endif
            npix_tot = npix_backgr*ndatlines
            if( npix_tot == 0 )then
                write(*,*) "ERROR, total nr of pixels cannot be zero"
                write(*,*) "npix_backgr: ", npix_backgr
                write(*,*) "ndatlines: ", ndatlines
                stop
            endif
            allocate( noise_pixels(npix_tot), stat=alloc_stat )
            call alloc_err('In: simple_extr_ptcls; boxdata, 3', alloc_stat)

            ! extract ctf info
            if( b%a%isthere('dfx') )then
                dfx = b%a%get(movie,'dfx')
                ctfres = b%a%get(movie,'ctfres')
                angast= 0.
                if( b%a%isthere('dfy') )then ! astigmatic CTF
                    if( .not. b%a%isthere('angast') ) stop 'need angle of astigmatism for CTF correction'
                    dfy = b%a%get(movie,'dfy')
                    angast = b%a%get(movie,'angast')
                endif
                ! set CTF info in outoris
                do i=1,ndatlines
                    if( pinds(i) > 0 )then
                        call outoris%set(pinds(i), 'dfx', dfx)
                        call outoris%set(pinds(i), 'ctfres', ctfres)
                        if( b%a%isthere('dfy') )then
                            call outoris%set(pinds(i), 'angast', angast)
                            call outoris%set(pinds(i), 'dfy', dfy)
                        endif
                    endif
                end do
                if( debug ) write(*,*) 'did set CTF parameters dfx/dfy/angast/ctfres: ', dfx, dfy, angast, ctfres
            endif
            
            ! loop over frames
            do frame=1,nframes
                
                ! read frame
                call img_frame%read(movienames(movie),frame)
                
                if( nframes > 1 )then
                    ! shift frame according to global shift (drift correction)
                    if( b%a%isthere(movie, 'x'//int2str(frame)) .and. b%a%isthere(movie, 'y'//int2str(frame))  )then
                        call img_frame%fwd_ft
                        x = b%a%get(movie, 'x'//int2str(frame))
                        y = b%a%get(movie, 'y'//int2str(frame))
                        call img_frame%shift(-x,-y)
                        call img_frame%bwd_ft
                    else
                        write(*,*) 'no shift parameters available for alignment of frames'
                        stop 'use simple_unblur_movies if you want to integrate frames'
                    endif
                endif
                
                ! extract the particle images & normalize
                if( nframes > 1 )then
                    framestack = 'framestack'//int2str_pad(frame,numlen)//p%ext
                else
                    framestack = 'sumstack'//p%ext
                endif
                cnt = 0
                do j=1,ndatlines ! loop over boxes
                    if( pinds(j) > 0 )then
                        ! modify coordinates if change in box (shift by half)
                        if( debug ) print *, 'original coordinate: ', boxdata(j,1:2) 
                        if( orig_box /= p%box ) boxdata(j,1:2) = boxdata(j,1:2)-real(p%box-orig_box)/2.
                        if( debug ) print *, 'shifted coordinate: ', boxdata(j,1:2) 
                        ! extract the window    
                        call img_frame%window(nint(boxdata(j,1:2)), p%box, b%img)
                        if( p%neg .eq. 'yes' ) call b%img%neg
                        pixels = b%img%extr_pixels(mskimg)
                        do k=1,size(pixels)
                            cnt = cnt+1
                            noise_pixels(cnt) = pixels(k)
                        end do
                        deallocate(pixels)
                        call b%img%write(trim(adjustl(framestack)), pinds(j))
                    endif
                end do
                med = median(noise_pixels(:cnt))
                call moment(noise_pixels(:cnt), ave, sdev, var, err)
                ! normalize the framestack according to the median and
                ! standard devation of the noise over the entire frame
                do j=1,ndatlines ! loop over boxes
                    if( pinds(j) > 0 )then
                        call b%img%read(trim(adjustl(framestack)), pinds(j))
                        if( p%noise_norm .eq. 'yes' )then
                            call b%img%norm_ext(med, sdev)
                        else
                            call b%img%norm
                        endif
                        call b%img%write(trim(adjustl(framestack)), pinds(j))
                    endif
                end do
            end do
            
            if( nframes > 1 )then
                ! create the sum and normalize it
                sumstack = 'sumstack'//p%ext
                cnt = 0
                do j=1,ndatlines ! loop over boxes
                    if( pinds(j) > 0 )then
                        b%img_copy = 0.
                        do frame=1,nframes ! loop over frames
                            framestack = 'framestack'//int2str_pad(frame,numlen)//p%ext
                            call b%img%read(trim(adjustl(framestack)), pinds(j))
                            call b%img_copy%add(b%img)
                        end do
                        call b%img_copy%write(trim(adjustl(sumstack)), pinds(j))
                        pixels = b%img_copy%extr_pixels(mskimg)
                        do k=1,size(pixels)
                            cnt = cnt+1
                            noise_pixels(cnt) = pixels(k)
                        end do
                        deallocate(pixels)
                    endif
                end do
                med = median(noise_pixels(:cnt))
                call moment(noise_pixels(:cnt), ave, sdev, var, err)
                ! normalize the sumstack according to the median and
                ! standard devation of the noise over the entire framesum
                do j=1,ndatlines ! loop over boxes
                    if( pinds(j) > 0 )then
                        call b%img%read(trim(adjustl(sumstack)), pinds(j))
                        if( p%noise_norm .eq. 'yes' )then
                            call b%img%norm_ext(med, sdev)
                        else
                            call b%img%norm
                        endif
                        call b%img%write(trim(adjustl(sumstack)), pinds(j))
                    endif
                end do
            endif
            
            ! write output
            do j=1,ndatlines ! loop over boxes
                if( pinds(j) > 0 )then
                    call outoris%write(pinds(j), 'extr_ptcls_params.txt')
                endif
            end do
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_EXTR_PTCLS NORMAL STOP ****')
        contains
            
            function box_inside( ldim, coord, box ) result( inside )
                integer, intent(in) :: ldim(3), coord(2), box
                integer             :: fromc(2), toc(2)
                logical             :: inside
                if( p%outside .eq. 'yes' )then
                    inside = .true.
                    return
                endif
                fromc = coord+1     ! compensate for the c-range that starts at 0
                toc = fromc+(box-1) ! the lower left corner is 1,1
                inside = .true.
                if( toc(1) > ldim(1) .or. toc(2) > ldim(2) ) inside = .false.
            end function
    end subroutine exec_extr_ptcls
    
    ! PRIME2D METHODS
    
    subroutine exec_prime2D_init( self, cline )
        use simple_hadamard2D_matcher, only: prime2D_assemble_sums, prime2D_write_sums
        class(prime2D_init_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)  :: p
        type(build)   :: b
        integer       :: ncls_in_oritab, icls
        p = params(cline)                     ! parameters generated
        p%boxmatch = p%box                    !!!!!!!!!!!!!!!!!! 4 NOW
        call b%build_general_tbox(p, cline)   ! general objects built
        call b%build_hadamard_prime2D_tbox(p) ! 2D Hadamard matcher built
        write(*,'(a)') '>>> GENERATING INITIAL CLUSTER CENTERS'
        if( cline%defined('oritab') )then
            call b%a%remap_classes
            ncls_in_oritab = b%a%get_ncls()
            if( p%ncls < ncls_in_oritab ) stop 'Inputted ncls < ncls_in_oritab; not allowed!'
            if( p%ncls > ncls_in_oritab )then
                call b%a%expand_classes(p%ncls)
            endif
        else
            if( p%srch_inpl .eq. 'yes' )then
                call b%a%rnd_cls(p%ncls)
            else
                call b%a%rnd_cls(p%ncls, srch_inpl=.false.)
            endif
        endif
        p%oritab = 'prime2D_startdoc.txt'
        call b%a%write(p%oritab)
        if( cline%defined('filwidth') )then
            do icls=1,p%ncls
                call b%cavgs(icls)%bin_filament(p%filwidth)
            end do
        else
            call prime2D_assemble_sums(b, p, mul=p%mul)
        endif
        call prime2D_write_sums(b, p)
        ! end gracefully
        call simple_end('**** SIMPLE_PRIME2D_INIT NORMAL STOP ****')
    end subroutine exec_prime2D_init
    
    subroutine exec_prime2D( self, cline )
        use simple_hadamard2D_matcher, only: prime2D_exec
        class(prime2D_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer      :: i, startit, ncls_from_refs, lfoo(3)
        logical      :: converged=.false.
        integer      :: err ! CUDA err variable for the return function calls
        call timestamp()
        ! call start_Alltimers_cpu()
        ! starting the cuda environment
        call simple_cuda_init(err)
        if( err .ne. RC_SUCCESS ) write(*,*) 'cublas init failed'
        p = params(cline)                     ! parameters generated
        p%boxmatch = p%box                    !!!!!!!!!!!!!!!!!! 4 NOW
        call b%build_general_tbox(p, cline)   ! general objects built
        call b%build_hadamard_prime2D_tbox(p) ! 2D Hadamard matcher built
        if( p%srch_inpl .eq. 'no' )then
            if( .not. cline%defined('oritab') )then
                stop 'need oritab for this mode (srch_inpl=no) of execution!'
            endif
        endif
        if( cline%defined('refs') )then
            call find_ldim_nptcls(p%refs, lfoo, ncls_from_refs)
            ! consistency check
            if( p%ncls /=  ncls_from_refs ) stop 'nrefs /= inputted ncls'
        endif
        ! execute
        if( cline%defined('part') )then
            if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
            call prime2D_exec(b, p, cline, 0, converged) ! partition or not, depending on 'part'       
        else
            startit = 1
            if( cline%defined('startit') ) startit = p%startit
            do i=startit,p%maxits
                call prime2D_exec(b, p, cline, i, converged)
                if(converged) exit
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_PRIME2D NORMAL STOP ****')
        !******************************************************
        !    Environment shutdown
        !******************************************************
        ! shutting down CUDA
        call simple_cuda_shutdown()
        ! shutting down timers
        ! call stop_Alltimers_cpu()
    end subroutine exec_prime2D
    
    subroutine exec_cavgassemble( self, cline )
        use simple_hadamard2D_matcher, only: prime2D_assemble_sums_from_parts, prime2D_write_sums
        class(cavgassemble_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)                                 :: p
        type(build)                                  :: b
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        call b%build_hadamard_prime2D_tbox(p)
        call prime2D_assemble_sums_from_parts(b, p)
        call prime2D_write_sums( b, p, p%which_iter)
        ! end gracefully
        call simple_end('**** SIMPLE_CAVGASSEMBLE NORMAL STOP ****')
    end subroutine exec_cavgassemble
    
    subroutine exec_check2D_conv( self, cline )
        class(check2D_conv_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        logical      :: converged
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        p%ncls    = b%a%get_ncls()
        converged = b%conv%check_conv2D()   ! convergence check
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK2D_CONV STOP ****')
    end subroutine exec_check2D_conv
    
    subroutine exec_rank_cavgs( self, cline )
        use simple_oris, only: oris
        class(rank_cavgs_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params)         :: p
        type(build)          :: b
        integer              :: iclass
        integer, allocatable :: order(:)
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        p%ncls   = p%nptcls
        p%nptcls = nlines(p%oritab) 
        call b%a%new(p%nptcls)
        call b%a%read(p%oritab)
        order = b%a%order_cls()
        do iclass=1,p%ncls
            write(*,'(a,1x,i5,1x,a,i5)') 'CLASS:', order(iclass), 'POP:', b%a%get_clspop(order(iclass)) 
            call b%img%read(p%stk, order(iclass))
            call b%img%write(p%outstk, iclass)
        end do
        call simple_end('**** SIMPLE_RANK_CAVGS NORMAL STOP ****')
    end subroutine exec_rank_cavgs
    
    ! PRIME3D METHODS
    
    subroutine exec_resrange( self, cline )
        use simple_hadamard3D_matcher, only: prime3D_find_resrange
        class(resrange_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        call timestamp()
        call start_Alltimers_cpu()
        if( cline%defined('box') .or. cline%defined('moldiam') )then
            p = params(cline)                     ! parameters generated
            call b%build_general_tbox(p, cline)   ! general objects built
            call b%build_hadamard_prime3D_tbox(p) ! prime objects built
            call prime3D_find_resrange( b, p, p%lp, p%lpstop )
            write(*,'(A,1X,F5.1)') '>>> LP START:', p%lp
            write(*,'(A,2X,F5.1)') '>>> LP STOP:',  p%lpstop
            write(*,'(A,2X,F5.1)') '>>> HP:',       p%hp
        else
            stop 'need either box size or moldiam to estimate resrange; simple_resrange'
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_RESRANGE NORMAL STOP ****')
        !*******************************************************************************
        !    Environment shutdown
        !
        !*******************************************************************************
        !shutting down the timers
        call stop_Alltimers_cpu()
    end subroutine exec_resrange
    
    subroutine exec_npeaks( self, cline )
        use simple_oris, only: oris
        class(npeaks_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(build)  :: b
        type(params) :: p
        p = params(cline) ! parameters generated
        call b%build_general_tbox(p, cline)
        p%npeaks = min(10,b%e%find_npeaks(p%lp, p%moldiam))
        write(*,'(A,1X,I4)') '>>> NPEAKS:', p%npeaks
        ! end gracefully
        call simple_end('**** SIMPLE_NPEAKS NORMAL STOP ****')
    end subroutine exec_npeaks
    
    subroutine exec_nspace(self,cline)
        use simple_math, only: resang
        use simple_oris, only: oris
        class(nspace_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(oris)   :: o
        type(params) :: p
        real         :: ares
        integer      :: i
        p = params(cline) ! parameters generated
        do i=500,5000,500
            o = oris(i)
            call o%spiral
            ares = o%find_angres()
            write(*,'(A,1X,I7,1X,A,1X,F5.2)') 'NR OF PROJDIRS:', i, 'RESOLUTION:', resang(ares, p%moldiam)
        end do
        call simple_end('**** SIMPLE_NSPACE NORMAL STOP ****')
    end subroutine exec_nspace
    
    subroutine exec_prime3D_init( self, cline )
        use simple_hadamard3D_matcher, only: gen_random_model, prime3D_find_resrange
        class(prime3D_init_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)       :: p
        type(build)        :: b
        integer, parameter :: MAXIMGS=1000
        call timestamp()
        call start_Alltimers_cpu()
        p = params(cline) ! parameters generated
        if( p%l_xfel )then
            if( cline%defined('msk') .or. cline%defined('mw') .or.&
            cline%defined('nvox') .or. cline%defined('mskfile') )then
                stop 'no mask allowed when processing XFEL patterns; simple_prime3D_init'
            endif
        endif
        if( p%ctf .ne. 'no')then
            if( .not. cline%defined('deftab') )&
            &stop 'need texfile with defocus/astigmatism values for ctf .ne. no mode exec'
        endif
        call b%build_general_tbox(p, cline)   ! general objects built
        call b%build_hadamard_prime3D_tbox(p) ! prime3D objects built
        ! determine resolution range
        if( cline%defined('lp') ) call prime3D_find_resrange( b, p, p%lp, p%lpstop )
        ! determine the number of peaks
        if( .not. cline%defined('npeaks') ) p%npeaks = min(10,b%e%find_npeaks(p%lp, p%moldiam))
        ! generate the random model
        if( cline%defined('nran') )then
            call gen_random_model( b, p, p%nran )
        else
            if( p%nptcls > MAXIMGS )then
                 call gen_random_model( b, p, MAXIMGS )
            else
                call gen_random_model( b, p )
            endif
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_PRIME3D_INIT NORMAL STOP ****')
        !*******************************************************************************
        !    Environment shutdown
        !
        !*******************************************************************************
        ! shutting down the timers
        call stop_Alltimers_cpu()
    end subroutine exec_prime3D_init
    
    subroutine exec_multiptcl_init( self, cline )
        use simple_rec_master, only: exec_rec, exec_eorec
        class(multiptcl_init_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        p = params(cline) ! constants & derived constants produced
        call b%build_general_tbox(p, cline)
        if( cline%defined('state2split') )then
            if( cline%defined('oritab') )then
                p%nstates = b%a%get_nstates()
                call b%a%split_state(p%state2split)
                p%nstates = p%nstates+1
            else
                stop 'Need oritab to be defined when state2split is defined on command line; simple_multiptcl_init'
            endif
        else
            if( p%nstates < 2 ) stop 'Nonsensical to have nstates < 2; simple_multiptcl_init'
            call b%a%rnd_states(p%nstates)
        endif
        if( p%errify.eq.'yes' ) call errify_oris
        if( p%norec .eq. 'no' )then
            if( cline%defined('lp') )then
                call b%build_rec_tbox(p)
                call exec_rec(b, p, cline, 'startvol')
            else
                call b%build_eo_rec_tbox(p)
                call exec_eorec(b, p, cline, 'startvol')
            endif
        endif
        if( p%zero .eq. 'yes' ) call b%a%set_all('corr', 0.)
        call b%a%write('multiptcl_startdoc.txt')
        ! end gracefully
        call simple_end('**** SIMPLE_MULTIPTCL_INIT NORMAL STOP ****')    

          contains

            subroutine errify_oris()
                real :: sherr,angerr
                sherr = 3.
                if( cline%defined('trs') ) sherr = p%trs
                angerr = 15.
                write(*,'(A,F6.2,A)')'>>> IN-PLANE SHIFT   ERROR INTRODUCED: ',sherr,' PIXELS'
                write(*,'(A,F6.2,A)')'>>> IN-PLANE ANGULAR ERROR INTRODUCED: ',angerr,' DEGREES'
                call b%a%introd_alig_err( angerr, sherr )
            end subroutine errify_oris

    end subroutine exec_multiptcl_init
    
    subroutine exec_prime3D( self, cline )
        use simple_hadamard3D_matcher, only: prime3D_exec, prime3D_find_resrange
        use simple_file_utils
        use simple_file_defs
        !temporary
        !use simple_err_defs
        use simple_file_highlev
        !use simple_eglossary
        !use simple_error_handling
        !use simple_dynamic_memory
        !use simple_systemQuery_cpu
        !use simple_deviceQuery_gpu
        class(prime3D_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)      :: p
        type(build)       :: b
        integer           :: i, startit
        logical           :: update_res=.false., converged=.false.
        real              :: lpstart, lpstop
        ! benchmark control logicals data structure
        type(t_bench)     :: s_bench
        ! filename string
        character(len=10) :: char_out
        character(len=80) :: tmr_name
        ! file handlers
        integer           :: unt
        integer           :: length
        ! CUDA err variable for the return function calls
        integer           :: err
        !function calls
        integer           :: get_length_of_string_c
        integer           :: convert_int2char_pos_c
        integer           :: convert_int2char_indexed_c
        integer           :: strlen
        call timestamp()
        call start_Alltimers_cpu()
        call simple_file_lib_initialise()
        ! starting the cuda environment
        call simple_cuda_init(err)
        if( err .ne. RC_SUCCESS ) write(*,*) 'cublas init failed'
        p = params(cline) ! parameters generated
        if( p%l_xfel )then
            if( cline%defined('msk') .or. cline%defined('mw') .or.&
            cline%defined('nvox') .or. cline%defined('automsk') )then
                stop 'no mask allowed when processing XFEL patterns; simple_prime3D'
            endif
        endif
        if( str_has_substr(p%refine,'neigh') .or. str_has_substr(p%refine,'qcont'))then
            if( .not. cline%defined('oritab') )then
                stop 'need oritab input for execution of prime3D with refine mode'
            endif
        endif
        call b%build_general_tbox(p, cline)   ! general objects built
        if( .not. cline%defined('eo') ) p%eo = 'no' ! default
        if( p%eo .eq. 'yes' ) p%dynlp = 'no'    
        if( cline%defined('lp') .or. cline%defined('find')&
        .or. p%eo .eq. 'yes' .or. p%dynlp .eq. 'yes' )then
            ! alles ok!
        else
           stop 'need a starting low-pass limit (set lp or find)!'
        endif
        call b%build_hadamard_prime3D_tbox(p) ! prime3D objects built
        if( cline%defined('part') )then
            if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
            if( cline%defined('find') ) p%lp = b%img%get_lp(p%find)
            call prime3D_exec(b, p, cline, 0, update_res, converged) ! partition or not, depending on 'part'
        else
            if( p%dynlp .eq. 'yes' )then
                call prime3D_find_resrange( b, p, lpstart, lpstop ) ! determine resolution range
                if( cline%defined('lpstart') )then
                    p%lp = p%lpstart
                else
                    p%lp = lpstart
                endif
                if( cline%defined('lpstop') )then
                    ! defined aleady
                else
                    p%lpstop = lpstop
                endif
            endif
            p%find = int((real(p%box-1)*p%smpd)/p%lp)
            startit = 1
            if( cline%defined('startit') ) startit = p%startit
            call start_timer_cpu('z_prime3D_exec')
            length = get_length_of_string_c(p%maxits)
            unt = 1
            do i=startit,p%maxits
                unt = 100 + i
                !err = convert_int2char_pos_c(char_out,i)
                err = convert_int2char_indexed_c(char_out,i,startit,p%maxits)
                tmr_name = 'z_prime3D_iter'
                tmr_name = tmr_name(1:strlen(tmr_name))//char_out(1:strlen(char_out))
                tmr_name = tmr_name(1:strlen(tmr_name))//".asc"
                write(*,*) "tmr_name: ",tmr_name
                if( b%s_bench%bench_write_i==1 .or. ibench_write.eqv..true. )&
                call file_open(tmr_name,unt,'unknown','asis','readwrite')
                if( b%s_bench%bench_i==1 .or. ibench.eqv..true. ) call start_timer_cpu(tmr_name)
                call prime3D_exec(b, p, cline, i, update_res, converged)
                if( b%s_bench%bench_i==1 .or. ibench.eqv..true. ) call stop_timer_cpu(tmr_name)
                if( update_res )then
                    p%find = p%find+p%fstep
                    p%lp = max(p%lpstop,b%img%get_lp(p%find))
                endif
                if( converged ) exit
            end do
            call stop_timer_cpu('z_prime3D_exec')
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_PRIME3D NORMAL STOP ****')
        !******************************************************
        !    Environment shutdown
        !
        !******************************************************
        !shutting down the environment
        call simple_cuda_shutdown()
        !shutting down the timers
        call stop_Alltimers_cpu()
    end subroutine exec_prime3D

    subroutine exec_cont3D( self, cline )
        use simple_cont3D_matcher, only: cont3D_exec
        class(cont3D_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer       :: i, startit
        logical       :: converged=.false.
        p = params(cline)                     ! parameters generated
        if( p%xfel .eq. 'yes' )then
            if( cline%defined('msk') .or. cline%defined('mw') .or.&
            cline%defined('nvox') .or. cline%defined('automsk') )then
                stop 'no mask allowed when processing XFEL patterns; simple_prime3D'
            endif
        endif
        call b%build_general_tbox(p, cline)
        call b%build_cont3D_tbox(p)
        if( cline%defined('part') )then
            if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
            call cont3D_exec(b, p, cline, 0, converged) ! partition or not, depending on 'part'
        else
            startit = 1
            if( cline%defined('startit') ) startit = p%startit
            do i=startit,p%maxits
                call cont3D_exec(b, p, cline, i, converged)
                if(converged) exit
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CONT3D NORMAL STOP ****')
    end subroutine exec_cont3D
    
    subroutine exec_check3D_conv( self, cline )
        use simple_math,    only: rad2deg, get_lplim
        class(check3D_conv_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)      :: p
        type(build)       :: b
        real, allocatable :: maplp(:)
        integer           :: istate, loc(1)
        logical           :: here, limset, converged
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        ! nstates consistency check
        if( cline%defined('nstates') )then
            if( p%nstates /= b%a%get_nstates() ) stop 'Inconsistent number of states between command-line and oritab'
        endif
        limset = .false.
        if( p%eo .eq. 'yes' )then
            allocate( maplp(p%nstates) )
            do istate=1,p%nstates
                p%fsc = 'fsc_state'//int2str_pad(istate,2)//'.bin'
                inquire(file=p%fsc, exist=here)
                if( here )then
                    b%fsc(istate,:) = file2rarr(p%fsc)
                    maplp(istate)   = max(b%img%get_lp(get_lplim(b%fsc(istate,:))),2.*p%smpd)
                else
                    write(*,*) 'Tried to open the fsc file: ', trim(p%fsc)
                    stop 'but it does not exist!'
                endif
            enddo
            loc     = maxloc( maplp )
            p%state = loc(1)            ! state with worst low-pass
            p%lp    = maplp( p%state )  ! worst lp
            p%fsc   =  'fsc_state'//int2str_pad(p%state,2)//'.bin'
            deallocate(maplp)
            limset = .true.
        endif
        ! Let find override the command line input lp (if given)
        if( .not. limset .and. cline%defined('find') )then
            p%lp = b%img%get_lp(p%find)
            limset = .true.
        endif
        ! Method for setting lp with lowest priority is lp on the command line
        if( cline%defined('lp') ) limset = .true.
        ! If we arrived here and the limit wasn't set: fall over
        if( limset )then
            ! we are happy
        else
            ! we fall over
            stop 'No method available to set low-pass limit! ABORTING...'
        endif
        ! calculate angular threshold
        p%athres = rad2deg(atan(p%lp/(p%moldiam/2.)))
        ! check convergence
        converged = b%conv%check_conv3D()
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK3D_CONV STOP ****')
    end subroutine exec_check3D_conv
    
    ! COMMON-LINES METHODS
    
    subroutine exec_comlin_smat( self, cline )
        use simple_comlin_sym  ! singleton
        use simple_comlin_corr ! singleton
        use simple_ori,        only: ori
        use simple_imgfile,    only: imgfile
        use simple_comlin,     only: comlin
        class(comlin_smat_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params), target          :: p
        type(build), target           :: b
        type(ori)                     :: orientation_best
        integer                       :: iptcl, jptcl, alloc_stat, funit, io_stat
        integer                       :: cnt, ntot, npairs, ipair, fnr
        real, allocatable             :: corrmat(:,:), corrs(:)
        integer, allocatable          :: pairs(:,:)
        logical                       :: debug=.false.
        character(len=:), allocatable :: fname
        p = params(cline, .false.)                            ! constants & derived constants produced
        call b%build_general_tbox(p, cline, .false., .true.)  ! general objects built (no oritab reading)
        allocate(b%imgs_sym(p%nptcls), stat=alloc_stat)
        call alloc_err('In: simple_comlin_smat, 1', alloc_stat)
        if( debug ) print *, 'analysing this number of objects: ', p%nptcls
        do iptcl=1,p%nptcls
            call b%imgs_sym(iptcl)%new([p%box,p%box,1], p%smpd, p%imgkind)
            call b%imgs_sym(iptcl)%read(p%stk, iptcl)
            ! apply a soft-edged mask
            call b%imgs_sym(iptcl)%mask(p%msk, 'soft')
            ! Fourier transform
            call b%imgs_sym(iptcl)%fwd_ft
        end do
        b%clins = comlin(b%a, b%imgs_sym)
        if( cline%defined('part') )then
            npairs = p%top-p%fromp+1
            if( debug ) print *, 'allocating this number of similarities: ', npairs
            allocate(corrs(p%fromp:p%top), pairs(p%fromp:p%top,2), stat=alloc_stat)
            call alloc_err('In: simple_comlin_smat, 1', alloc_stat)
            ! read the pairs
            funit = get_fileunit()
            allocate(fname, source='pairs_part'//int2str_pad(p%part,p%numlen)//'.bin')
            if( .not. file_exists(fname) )then
                write(*,*) 'file: ', fname, ' does not exist!'
                write(*,*) 'If all pair_part* are not in cwd, please execute simple_split_pairs to generate the required files'
                stop 'I/O error; simple_comlin_smat'
            endif
            open(unit=funit, status='OLD', action='READ', file=fname, access='STREAM')
            if( debug ) print *, 'reading pairs in range: ', p%fromp, p%top
            read(unit=funit,pos=1,iostat=io_stat) pairs(p%fromp:p%top,:)
            ! Check if the read was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(simple_comlin_smat): I/O error ', io_stat, ' when reading file: ', fname
                stop 'I/O error; simple_comlin_smat'
            endif
            close(funit)
            deallocate(fname)
            ! calculate the similarities
            call comlin_sym_init(b, p)
            cnt = 0
            do ipair=p%fromp,p%top
                cnt = cnt+1
                call progress(cnt, npairs)
                p%iptcl = pairs(ipair,1)
                p%jptcl = pairs(ipair,2)
                call comlin_sym_axis(orientation_best, 'pair', .false.)
                corrs(ipair) = orientation_best%get('corr')
            end do
            if( debug ) print *, 'did set this number of similarities: ', cnt
            ! write the similarities
            funit = get_fileunit()
            allocate(fname, source='similarities_part'//int2str_pad(p%part,p%numlen)//'.bin')
            open(unit=funit, status='REPLACE', action='WRITE', file=fname, access='STREAM')
            write(unit=funit,pos=1,iostat=io_stat) corrs(p%fromp:p%top)
            ! Check if the write was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(simple_comlin_smat): I/O error ', io_stat, ' when writing to: ', fname
                stop 'I/O error; simple_comlin_smat'
            endif
            close(funit)
            deallocate(fname, corrs, pairs)
            fnr = get_fileunit()
            open(unit=fnr, FILE='JOB_FINISHED_'//int2str_pad(p%part,p%numlen), STATUS='REPLACE', action='WRITE', iostat=io_stat)
            call fopen_err( 'In: simple_comlin_smat', io_stat )
            write(fnr,'(A)') '**** SIMPLE_COMLIN_SMAT NORMAL STOP ****'
            close(fnr)
        else
            allocate(corrmat(p%nptcls,p%nptcls), stat=alloc_stat)
            call alloc_err('In: simple_comlin_smat, 3', alloc_stat)
            corrmat = 1.
            ntot = (p%nptcls*(p%nptcls-1))/2
            call comlin_sym_init(b, p)
            cnt = 0
            do iptcl=1,p%nptcls-1
                do jptcl=iptcl+1,p%nptcls
                    cnt = cnt+1
                    call progress(cnt, ntot)
                    p%iptcl = iptcl
                    p%jptcl = jptcl
                    call comlin_sym_axis(orientation_best, 'pair', .false.)
                    corrmat(iptcl,jptcl) = orientation_best%get('corr')
                    corrmat(jptcl,iptcl) = corrmat(iptcl,jptcl)
                end do
            end do
            call progress(ntot, ntot)
            funit = get_fileunit()
            open(unit=funit, status='REPLACE', action='WRITE', file='clin_smat.bin', access='STREAM')
            write(unit=funit,pos=1,iostat=io_stat) corrmat
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to clin_smat.bin'
                stop 'I/O error; simple_comlin_smat'
            endif
            close(funit)
            deallocate(corrmat)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_COMLIN_SMAT NORMAL STOP ****')
    end subroutine exec_comlin_smat
    
    subroutine exec_symsrch( self, cline )
        use simple_oris,      only: oris
        use simple_symsrcher, only: symsrch_master
        class(symsrch_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        type(oris)   :: o
        p = params(cline) ! parameters generated
        if( cline%defined('stk') )then
            p%nptcls = 1
            p%nspace = 1
        endif
        call b%build_general_tbox(p, cline, .true., .true.) ! general objects built (no oritab reading)
        call b%build_comlin_tbox(p)                         ! objects for common lines based alignment built
        call symsrch_master( cline, p, b, o )
        call o%write(p%outfile)
        ! end gracefully
        call simple_end('**** SIMPLE_SYMSRCH NORMAL STOP ****')
    end subroutine exec_symsrch

    ! MASKING METHODS

    subroutine exec_mask( self, cline )
        use simple_procimgfile,   only: mask_imgfile
        class(mask_commander), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        type(build)                :: b
        type(params)               :: p
        type(automask2D_commander) :: automask2D
        type(automask3D_commander) :: automask3D
        logical                    :: here
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        if( cline%defined('stk') .and. cline%defined('vol1') )stop 'Cannot operate on images AND volume at once'
        if( p%automsk.eq.'yes' .and..not.cline%defined('mw')  )stop 'Missing mw argument for automasking'
        if( p%msktype.ne.'soft' .or. p%msktype.ne.'hard' .or. p%msktype.ne.'cavg' )stop 'Invalid mask type'
        call b%vol%new([p%box,p%box,p%box], p%smpd) ! reallocate vol (boxmatch issue)
        if( cline%defined('stk') )then
            ! 2D
            if( p%automsk.eq.'yes' )then
                ! auto
                if( .not. cline%defined('amsklp') )call cline%set('amsklp', 25.)
                if( .not. cline%defined('edge')   )call cline%set('edge', 20.)
                call exec_automask2D( automask2D, cline )
            else if( cline%defined('msk') .or. cline%defined('inner') )then
                ! spherical
                if( cline%defined('inner') )then
                    if( cline%defined('width') )then
                        call mask_imgfile(p%stk, p%outstk, p%msk, inner=p%inner, width=p%width, which=p%msktype)
                    else
                        call mask_imgfile(p%stk, p%outstk, p%msk, inner=p%inner, which=p%msktype)
                    endif
                else
                    call mask_imgfile(p%stk, p%outstk, p%msk, which=p%msktype)
                endif
            else
                stop 'Nothing to do!'
            endif
        else if( cline%defined('vol1') )then
            ! 3D
            here = .false.
            inquire(FILE=p%vols(1), EXIST=here)
            if( .not.here )stop 'Cannot find input volume'
            call b%vol%read(p%vols(1))
            if( p%automsk.eq.'yes' )then
                ! auto
                call exec_automask3D( automask3D, cline )
            else if( cline%defined('msk') )then
                ! spherical
                if( cline%defined('inner') )then
                    if( cline%defined('width') )then
                        call b%vol%mask(p%msk, p%msktype, inner=p%inner, width=p%width)
                    else
                        call b%vol%mask(p%msk, p%msktype, inner=p%inner)
                    endif
                else
                    call b%vol%mask(p%msk, p%msktype)
                endif
                if( p%outvol .ne. '' )call b%vol%write(p%outvol, del_if_exists=.true.)
            else
                stop 'Nothing to do!'
            endif
        else
            stop 'No input images(s) or volume provided'
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_MASK NORMAL STOP ****')
    end subroutine exec_mask

    subroutine exec_automask2D( self, cline )
        use simple_masker, only: automask2D
        class(automask2D_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer      :: iptcl
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        write(*,'(A,F8.2,A)') '>>> AUTOMASK LOW-PASS:',        p%amsklp, ' ANGSTROMS'
        write(*,'(A,I3,A)')   '>>> AUTOMASK SOFT EDGE WIDTH:', p%edge,   ' PIXELS'
        p%outstk = add2fbody(p%stk, p%ext, 'msk') 
        do iptcl=1,p%nptcls
            call b%img%read(p%stk, iptcl)
            call automask2D(b%img, p, b%img_msk)
            call b%img_msk%write('automasks2D'//p%ext, iptcl)
            call b%img%write(p%outstk, iptcl)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_AUTOMASK2D NORMAL STOP ****')
    end subroutine exec_automask2D
    
    subroutine exec_automask3D( self, cline )
        use simple_masker, only: automask
        class(automask3D_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer      :: istate
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        write(*,'(A,F14.1,A)') '>>> AUTOMASK LOW-PASS:',            p%amsklp,  ' ANGSTROMS'
        write(*,'(A,I7,A)')    '>>> AUTOMASK SOFT EDGE WIDTH:',     p%edge,    ' PIXEL(S)'
        write(*,'(A,I3,A)')    '>>> AUTOMASK BINARY LAYERS WIDTH:', p%binwidth,' PIXEL(S)'
        do istate=1,p%nstates
            p%masks(istate)    = 'automask_state'//int2str_pad(istate,2)//p%ext
            p%vols_msk(istate) = add2fbody(p%vols(istate), p%ext, 'msk')
            call b%vol%read(p%vols(istate))
            call automask(b, p, cline, b%vol,  b%mskvol, p%vols_msk(istate), p%masks(istate))
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_AUTOMASK3D NORMAL STOP ****')
    end subroutine exec_automask3D
    
    ! RECONSTRUCTION METHODS

    subroutine exec_eo_recvol( self, cline )
        use simple_rec_master, only: exec_eorec
        class(eo_recvol_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        call b%build_eo_rec_tbox(p)         ! eo_reconstruction objs built
        call exec_eorec(b, p, cline)
        ! end gracefully
        call simple_end('**** SIMPLE_EO_RECVOL NORMAL STOP ****')    
    end subroutine exec_eo_recvol

    subroutine exec_eo_volassemble( self, cline )
        use simple_eo_reconstructor, only: eo_reconstructor
        class(eo_volassemble_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(params)                  :: p
        type(build)                   :: b
        type(eo_reconstructor)        :: eorecvol_read
        character(len=:), allocatable :: fname
        real, allocatable             :: res05s(:), res0143s(:)
        real                          :: res
        integer                       :: part, s, alloc_stat, cnt, n, ss, state4name
        logical                       :: debug=.false.
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        if( cline%defined('nstates') )then
            if( p%nstates /= b%a%get_nstates() ) stop 'Inconsistent number of states between command-line and oritab'
        endif
        call b%build_eo_rec_tbox(p)   ! reconstruction toolbox built
        allocate(res05s(p%nstates), res0143s(p%nstates), stat=alloc_stat)
        call alloc_err("In: simple_eo_volassemble", alloc_stat)
        ! rebuild b%vol according to box size (beacuse it is otherwise boxmatch)
        call b%vol%new([p%box,p%box,p%box], p%smpd)
        call eorecvol_read%new(p)
        n = p%nstates*p%npart
        cnt = 0
        do ss=1,p%nstates
            if( cline%defined('state') )then
                s = 1
            else
                s = ss
            endif
            if( debug ) write(*,*) 'processing state: ', s
            call b%eorecvol%reset_all
            do part=1,p%npart
                cnt = cnt+1
                call progress(cnt,n)
                if( cline%defined('state') )then
                    state4name = p%state
                else
                    state4name = s
                endif
                allocate(fname, source='recvol'//'_state'//int2str_pad(state4name,2)//'_part'//int2str_pad(part,p%numlen))
                if( debug ) write(*,*) 'processing file: ', fname
                call assemble(fname)
                deallocate(fname)
            end do
            call normalize('recvol_state'//int2str_pad(state4name,2))
        end do
        ! set the resolution limit according to the worst resolved model
        res  = maxval(res0143s)
        p%lp = max( p%lpstop,res )
        write(*,'(a,1x,F6.2)') '>>> LOW-PASS LIMIT:', p%lp
        write(0,'(a)') "GENERATED VOLUMES: recvol*.ext"
        ! end gracefully
        call simple_end('**** SIMPLE_EO_VOLASSEMBLE NORMAL STOP ****')
        
        contains

            subroutine assemble( fbody )
                character(len=*), intent(in) :: fbody
                call eorecvol_read%read_eos(fbody)
                ! sum the Fourier coefficients
                call b%eorecvol%sum(eorecvol_read)
            end subroutine
            
            subroutine normalize( recnam )
                use simple_image, only: image
                character(len=*), intent(in)  :: recnam
                call b%eorecvol%sum_eos
                call b%eorecvol%sampl_dens_correct_eos(s)
                call b%eorecvol%get_res(res05s(s), res0143s(s))
                call b%eorecvol%sampl_dens_correct_sum(b%vol)
                call b%eorecvol%write_eos(recnam)
                call b%vol%write(recnam//p%ext, del_if_exists=.true.)
            end subroutine
    end subroutine exec_eo_volassemble
    
    subroutine exec_recvol( self, cline )
        use simple_rec_master, only: exec_rec
        class(recvol_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        call b%build_rec_tbox(p)            ! reconstruction objects built
        call exec_rec(b, p, cline)
        ! end gracefully
        call simple_end('**** SIMPLE_RECVOL NORMAL STOP ****')    
    end subroutine exec_recvol
    
    subroutine exec_volassemble( self, cline )
        use simple_reconstructor, only: reconstructor
        class(volassemble_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)                  :: p
        type(build)                   :: b
        character(len=:), allocatable :: fbody
        integer                       :: part, s, ss, endit, i, cnt, state4name
        type(reconstructor)           :: recvol_read
        logical                       :: here(2)
        logical, parameter            :: debug=.false.
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        call b%build_rec_tbox(p)            ! reconstruction toolbox built
        ! rebuild b%vol according to box size (beacuse it is otherwise boxmatch)
        call b%vol%new([p%box,p%box,p%box], p%smpd)
        if( cline%defined('find') )then
            p%lp = b%img%get_lp(p%find)
        endif
        call recvol_read%new([p%boxpd,p%boxpd,p%boxpd], p%smpd)
        call recvol_read%alloc_rho(p)
        endit = 1
        if( p%eo .eq. 'yes' ) endit = 2
        cnt = 0
        do ss=1,p%nstates
            if( cline%defined('state') )then
                s = 1
            else
                s = ss
            endif
            if( debug ) write(*,*) 'processing state: ', s
            call b%recvol%reset
            do part=1,p%npart
                cnt = cnt+1
                call progress(cnt,p%nstates*p%npart)
                if( cline%defined('state') )then
                    state4name = p%state
                else
                    state4name = s
                endif
                allocate(fbody, source='recvol'//'_state'//int2str_pad(state4name,2)//'_part'//int2str_pad(part,p%numlen))
                if( debug ) write(*,*) 'processing fbody: ', fbody
                do i=1,endit
                    if( cline%defined('even') .or. cline%defined('odd') )then
                        if( p%even .eq. 'yes' .and. p%odd .eq. 'no' )then
                            p%vols(s) = fbody//'_even'//p%ext
                            p%masks(s) = 'rho_'//fbody//'_even'//p%ext
                        else if( p%odd .eq. 'yes' .and. p%even .eq. 'no' )then
                            p%vols(s) = fbody//'_odd'//p%ext
                            p%masks(s) = 'rho_'//fbody//'_odd'//p%ext
                        else if( p%odd .eq. 'yes' .and. p%even .eq. 'yes' )then
                            stop 'ERROR! Cannot have even=yes and odd=yes simultaneously'
                        endif
                    else
                        if( p%eo .eq. 'yes' )then
                            if( i == 1 )then
                                p%vols(s) = fbody//'_odd'//p%ext
                                p%masks(s) = 'rho_'//fbody//'_odd'//p%ext
                            else
                                p%vols(s) = fbody//'_even'//p%ext
                                p%masks(s) = 'rho_'//fbody//'_even'//p%ext
                            endif   
                        else
                            p%vols(s)  = fbody//p%ext
                            p%masks(s) = 'rho_'//fbody//p%ext
                        endif
                    endif
                    call assemble(p%vols(s), p%masks(s))
                end do
                deallocate(fbody)
            end do
            if( p%even .eq. 'yes' .and. p%odd .eq. 'no' )then
                call normalize('recvol_state'//int2str_pad(s,2)//'_even'//p%ext)
            else if( p%odd .eq. 'yes' .and. p%even .eq. 'no' )then
                call normalize('recvol_state'//int2str_pad(s,2)//'_odd'//p%ext)
            else if( p%odd .eq. 'yes' .and. p%even .eq. 'yes' )then
                stop 'ERROR! Cannot have even=yes and odd=yes simultaneously'
            else
                call normalize('recvol_state'//int2str_pad(s,2)//p%ext)
            endif
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_VOLASSEMBLE NORMAL STOP ****')

        contains

            subroutine assemble( recnam, kernam )
                character(len=*), intent(in) :: recnam
                character(len=*), intent(in) :: kernam
                inquire(FILE=recnam, EXIST=here(1))
                inquire(FILE=kernam, EXIST=here(2))
                if( all(here) )then     
                    call recvol_read%read(recnam)
                    if( debug )then
                        if( recvol_read%contains_nans() )then
                            write(*,*) 'WARNING! recvol: ', recnam, 'contains NaN:s'
                        endif
                    endif
                    call recvol_read%read_rho(kernam)
                    if( debug )then
                        if( recvol_read%rho_contains_nans() )then
                            write(*,*) 'WARNING! kernel: ', kernam, 'contains NaN:s'
                        endif
                    endif
                    call b%recvol%sum(recvol_read)
                    if( debug )then
                        if( b%recvol%contains_nans() )then
                            write(*,*) 'WARRNING! summed image part contains NaN:s'
                        endif 
                        if( b%recvol%rho_contains_nans() )then
                            write(*,*) 'WARNING! summed kernel part contains NaN:s'
                        endif
                    endif
                else
                    if( .not. here(1) ) write(*,'(A,A,A)') 'WARNING! ', adjustl(trim(recnam)), ' missing'
                    if( .not. here(2) ) write(*,'(A,A,A)') 'WARNING! ', adjustl(trim(kernam)), ' missing'
                    return
                endif
            end subroutine
    
            subroutine normalize( recnam )
                character(len=*), intent(in) :: recnam
                call b%recvol%sampl_dens_correct
                call b%recvol%bwd_ft
                call b%recvol%clip(b%vol)
                call b%vol%write(recnam, del_if_exists=.true.)
            end subroutine
    end subroutine exec_volassemble
    
    ! CHECKER METHODS
    
    subroutine exec_check_box( self, cline )
        class(check_box_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(params) :: p
        p = params(cline) ! parameters generated
        write(*,'(A,1X,I7)') '>>> BOX:', p%box
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK_BOX NORMAL STOP ****')
    end subroutine exec_check_box

    subroutine exec_check_nptcls( self, cline )
        class(check_nptcls_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params) :: p
        p = params(cline) ! parameters generated
        write(*,'(A,1X,I7)') '>>> NPTCLS:', p%nptcls
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK_NPTCLS NORMAL STOP ****')
    end subroutine exec_check_nptcls
    
    subroutine exec_iminfo( self, cline)
        use simple_image,   only: image
        use simple_imgfile, only: imgfile
        class(iminfo_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params)      :: p
        type(image)       :: img
        integer           :: ldim(3), maxim, i, n_nans
        real              :: smpd, sdev, ave, minv, maxv
        p = params(cline) ! constants & derived constants produced
        if( cline%defined('fname') )then
            call find_ldim_nptcls(p%fname, ldim, maxim, doprint=.true.)
        endif
        p%box  = ldim(1)
        p%smpd = smpd          !! UNDEFINED
        call img%new([ldim(1),ldim(2),1],p%smpd)
        if( p%vis .eq. 'yes' .or. p%stats .ne. 'no' )then
            do i=1,maxim
                call img%read(p%fname, i)
                if( p%stats .ne. 'no' )then
                    call img%cure(maxv, minv, ave, sdev, n_nans)
                    if( p%stats .eq. 'print' .or. n_nans > 0 )then
                        write(*,*) '*********IMAGE********', i, '*******'
                        write(*,*) 'maxv = ',   maxv
                        write(*,*) 'minv = ',   minv
                        write(*,*) 'ave = ',    ave
                        write(*,*) 'sdev = ',   sdev
                        write(*,*) 'n_nans = ', n_nans
                    endif
                endif
                if( p%vis .eq. 'yes' ) call img%vis
            end do
        endif
        call simple_end('**** SIMPLE_IMINFO NORMAL STOP ****')
    end subroutine exec_iminfo
    
    ! VOLOPS METHODS
    
    subroutine exec_cenvol( self, cline )
        class(cenvol_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params)         :: p
        type(build)          :: b
        real, allocatable    :: shvec(:,:)
        integer              :: istate
        logical, parameter   :: debug=.false.
        p = params(cline)                           ! parameters generated
        call b%build_general_tbox(p, cline, .true.) ! general objects built
        ! center volume(s)
        allocate(shvec(p%nstates,3))
        do istate=1,p%nstates
            call b%vol%read(p%vols(istate))
            shvec(istate,:) = b%vol%center(p%amsklp,p%msk)
            if( debug )then
                call b%vol%shift(-shvec(istate,1), -shvec(istate,2), -shvec(istate,3))
                call b%vol%write('shifted_vol_state'//int2str(istate)//p%ext)
            endif
            ! transfer the 3D shifts to 2D
            call b%a%map3dshift22d(-shvec(istate,:), state=istate)
        end do
        call b%a%write(p%outfile)
        ! end gracefully
        call simple_end('**** SIMPLE_CENVOL NORMAL STOP ****')
    end subroutine exec_cenvol
    
    subroutine exec_postproc_vol(self,cline)
        use simple_estimate_ssnr ! singleton
        use simple_masker,       only: automask
        class(postproc_vol_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)      :: p
        type(build)       :: b
        real, allocatable :: fsc(:), spec_count3D(:), tmparr(:), pssnr_ctfsq3D(:), pssnr3D(:), sqrtssnr(:), optlp(:)
        integer           :: k, state=1
        p = params(cline, checkdistr=.false.) ! constants & derived constants produced, mode=2
        call b%build_general_tbox(p, cline)   ! general objects built
        if( file_exists(p%fsc) )then
            fsc = file2rarr(p%fsc)
            ! allocate CTF**2-dependent PSSNR term
            allocate(pssnr_ctfsq3D(size(fsc)))
            if( file_exists(p%ctfsqspec) )then
                ! get the count spectrum
                spec_count3D = b%vol%spectrum('count')
                ! get the ctfsq spectrum
                tmparr = file2rarr(p%ctfsqspec)
                ! calculate the CTF**2-dependent component of the PSSNR
                where( tmparr > 1e-6 )
                    pssnr_ctfsq3D = spec_count3D/tmparr
                else where
                    pssnr_ctfsq3D = 0.
                end where
            else
                pssnr_ctfsq3D = 1.
            endif
            allocate(pssnr3D(size(fsc)), sqrtssnr(size(fsc)))
            pssnr3D = estimate_pssnr3D(p%avr, fsc)*pssnr_ctfsq3D    
            optlp = ssnr2optlp(pssnr3D)
        else
            write(*,*) 'FSC file: ', trim(p%fsc), ' not in cwd'
            stop
        endif
        call b%vol%read(p%vols(state))
        call b%vol%fwd_ft
        if( cline%defined('bfac') )then
            call b%vol%apply_bfac(p%bfac)
        endif
        call b%vol%apply_filter(optlp)
        call b%vol%bwd_ft
        p%vols_msk(state) = add2fbody(trim(p%vols(state)), p%ext, 'msk')
        if( p%automsk .eq. 'yes' )then
            p%masks(state) = 'automask_state'//int2str_pad(state,2)//p%ext
            call automask(b, p, cline, b%vol, b%mskvol, p%vols_msk(state), p%masks(state))
        else
            call b%vol%mask(p%msk, 'soft')
        endif
        p%outvol = add2fbody(trim(p%vols(state)), p%ext, 'pproc')
        call b%vol%write(p%outvol)
        call simple_end('**** SIMPLE_POSTPROC_VOL NORMAL STOP ****')
    end subroutine exec_postproc_vol
    
    subroutine exec_projvol( self, cline )
        use simple_image, only: image
        class(projvol_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)             :: p
        type(build)              :: b
        type(image), allocatable :: imgs(:)
        integer                  :: i, loop_end
        real                     :: x, y, dfx, dfy, angast
        logical, parameter       :: debug=.false.
        if( .not. cline%defined('oritab') )then
            if( .not. cline%defined('nspace') ) stop 'need nspace (for number of projections)!'
        endif
        p = params(cline) ! parameters generated
        if( cline%defined('oritab') )then
            p%nptcls = nlines(p%oritab)
            call b%build_general_tbox(p, cline)
            call b%a%read(p%oritab)
            p%nspace = b%a%get_noris()
        else if( p%rnd .eq. 'yes' )then
            p%nptcls = p%nspace
            call b%build_general_tbox(p, cline)
            call b%a%rnd_oris(p%trs)
        else
            p%nptcls = p%nspace
            call b%build_general_tbox(p, cline)
            call b%a%spiral(p%nsym, p%eullims)
            if( cline%defined('trs') ) call b%a%rnd_inpls(p%trs)
        endif
        ! fix volumes and stacks
        if( p%xfel .eq. 'yes' )then
            call b%vol%read(p%vols(1), isxfel=.true.)
        else
            call b%vol%read(p%vols(1))
        endif
        if( debug ) print *, 'read volume'
        ! generate projections
        if( p%swap .eq. 'yes' ) call b%a%swape1e3
        if( p%mirr .eq. 'yes' ) call b%a%mirror2d
        if( cline%defined('top') )then
            imgs = b%proj%projvol(b%vol, b%a, p, p%top)
            loop_end = p%top
        else
            imgs = b%proj%projvol(b%vol, b%a, p)
            loop_end = p%nspace
        endif
        if( file_exists(p%outstk) ) call del_binfile(p%outstk)
        do i=1,loop_end
            if( cline%defined('oritab') .or. (p%rnd .eq. 'yes' .or. cline%defined('trs')) )then
                x = b%a%get(i, 'x')
                y = b%a%get(i, 'y')
                call imgs(i)%shift(x, y)
            endif
            if( p%ctf .ne. 'no' )then
                if( cline%defined('oritab') )then
                    dfx = b%a%get(i, 'dfx')
                    if( b%a%isthere('dfy') )then
                        dfy    = b%a%get(i, 'dfy')
                        angast = b%a%get(i, 'angast')
                    else
                        dfy    = dfx
                        angast = 0.
                    endif
                else
                    dfx    = p%defocus
                    dfy    = p%defocus
                    angast = 0.
                    call b%a%set(i, 'dfx', dfx)
                endif
                if( cline%defined('bfac') )then
                    if( p%neg .eq. 'yes' )then
                        call b%tfun%apply(imgs(i), dfx, 'neg', dfy, angast, bfac=p%bfac)
                    else
                        call b%tfun%apply(imgs(i), dfx, 'ctf', dfy, angast, bfac=p%bfac)
                    endif
                else
                    if( p%neg .eq. 'yes' )then
                        call b%tfun%apply(imgs(i), dfx, 'neg', dfy, angast)
                    else
                        call b%tfun%apply(imgs(i), dfx, 'ctf', dfy, angast)
                    endif
                endif
            else if( p%neg .eq. 'yes' )then
                call imgs(i)%neg
            endif
            if( p%mirr .ne. 'no' )then
                if( p%mirr .ne. 'yes' ) call imgs(i)%mirror(p%mirr)
            endif
            call imgs(i)%write(p%outstk,i)
        end do
        call b%a%write('projvol_oris.txt')
        call simple_end('**** SIMPLE_PROJVOL NORMAL STOP ****')
    end subroutine exec_projvol
    
    subroutine exec_volaverager( self, cline )
        use simple_image, only: image
        class(volaverager_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer, allocatable               :: ptcls(:)
        character(len=STDLEN), allocatable :: volnames(:)
        type(image)                        :: vol_avg
        integer                            :: istate, ivol, nvols, funit_vols, numlen, ifoo
        character(len=:), allocatable      :: fname
        character(len=1)                   :: fformat
        logical, parameter                 :: debug=.true.
        p = params(cline) ! parameters generated
        ! read the volnames
        nvols = nlines(p%vollist)
        if( debug ) print *, 'number of volumes: ', nvols
        allocate(volnames(nvols))
        funit_vols = get_fileunit()
        open(unit=funit_vols, status='old', file=p%vollist)
        do ivol=1,nvols
            read(funit_vols,'(a256)') volnames(ivol)
            if( debug ) print *, 'read volname: ', volnames(ivol)
        end do
        close(funit_vols)
        ! find logical dimension
        call find_ldim_nptcls(volnames(1), p%ldim, ifoo)
        p%box  = p%ldim(1)
        ! build general toolbox
        call b%build_general_tbox(p, cline) ! general objects built
        ! figure out the file extension
        fformat = fname2format(volnames(1))
        select case(fformat)
            case('M')
                p%ext = '.mrc'
            case('S')
                p%ext = '.spi'
            case('D')
                p%ext = '.mrc'
            case('B')
                p%ext = '.mrc'
            case DEFAULT
                stop 'This file format is not supported by SIMPLE; simple_volaverager'
        end select
        if( debug ) print *, 'file extension: ', p%ext
        ! average the states
        call vol_avg%copy(b%vol)
        p%nstates = b%a%get_nstates()
        if( debug ) print *, 'number of states: ', p%nstates
        numlen = len(int2str(p%nstates))
        do istate=1,p%nstates
            if( debug ) print *, 'processing state: ', istate
            ptcls = nint(b%a%get_ptcls_in_state(istate))
            vol_avg = 0.
            do ivol=1,size(ptcls)
                call b%vol%read(volnames(ptcls(ivol)))
                if( debug ) print *, 'read volume: ', volnames(ptcls(ivol))
                call vol_avg%add(b%vol)
            end do
            call vol_avg%div(real(size(ptcls)))
            allocate(fname, source='sumvol_state'//int2str_pad(istate, numlen)//p%ext)
            if( debug ) print *, 'trying to write volume to file: ', fname
            call vol_avg%write(fname)
            deallocate(ptcls,fname)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_VOLAVERAGER NORMAL STOP ****')
    end subroutine exec_volaverager
    
    subroutine exec_volops( self, cline )
        class(volops_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        logical      :: here 
        p = params(cline,checkdistr=.false.)        ! constants & derived constants produced, mode=2
        call b%build_general_tbox(p, cline)         ! general objects built
        call b%vol%new([p%box,p%box,p%box], p%smpd) ! reallocate vol (boxmatch issue)
        inquire(FILE=p%vols(1), EXIST=here)
        if( here )then
            call b%vol%read(p%vols(1))
        else
            stop 'vol1 does not exists in cwd'
        endif
        if( p%guinier .eq. 'yes' )then
            if( .not. cline%defined('smpd') ) stop 'need smpd (sampling distance) input for Guinier plot'
            if( .not. cline%defined('hp')   ) stop 'need hp (high-pass limit) input for Guinier plot'
            if( .not. cline%defined('lp')   ) stop 'need lp (low-pass limit) input for Guinier plot'
            p%bfac = b%vol%guinier_bfac(p%hp, p%lp)
            write(*,'(A,1X,F8.2)') '>>> B-FACTOR DETERMINED TO:', p%bfac
        else
            if( cline%defined('neg')  ) call b%vol%neg
            if( cline%defined('snr')  ) call b%vol%add_gauran(p%snr)
            if( cline%defined('mirr') ) call b%vol%mirror(p%mirr)
            call b%vol%write(p%outvol, del_if_exists=.true.)
        endif
        call simple_end('**** SIMPLE_VOLOPS NORMAL STOP ****')
        ! end gracefully
        call simple_end('**** SIMPLE_VOLOPS NORMAL STOP ****')
    end subroutine exec_volops
    
    subroutine exec_volume_smat( self, cline )
        use simple_image, only: image
        class(volume_smat_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params), target               :: p
        integer                            :: funit, io_stat, cnt, npairs, npix, nvols
        integer                            :: ivol, jvol, ldim(3), alloc_stat, ipair, ifoo
        real, allocatable                  :: corrmat(:,:), corrs(:)
        integer, allocatable               :: pairs(:,:)
        type(image)                        :: vol1, vol2, mskimg
        character(len=STDLEN), allocatable :: vollist(:)
        character(len=:), allocatable      :: fname
        logical, parameter                 :: debug=.false.
        p      = params(cline, .false.) ! constants & derived constants produced
        nvols  = nlines(p%vollist)
        npairs = (nvols*(nvols-1))/2
        ! read in list of volumes
        allocate(vollist(nvols))
        funit = get_fileunit()
        open(unit=funit, status='old', file=p%vollist)
        do ivol=1,nvols
            read(funit,'(a256)') vollist(ivol)
            if( debug ) write(*,*) 'read volume: ', vollist(ivol)
        end do
        close(unit=funit)
        ! find logical dimension & make volumes for matching
        call find_ldim_nptcls(vollist(1), ldim, ifoo)
        if( debug ) write(*,*) 'found logical dimension: ', ldim
        call vol1%new(ldim,p%smpd)
        call vol2%new(ldim,p%smpd)
        if( debug ) write(*,*) 'allocated volumes'
        if( cline%defined('part') )then
            npairs = p%top-p%fromp+1
            if( debug ) print *, 'allocating this number of similarities: ', npairs
            allocate(corrs(p%fromp:p%top), pairs(p%fromp:p%top,2), stat=alloc_stat)
            call alloc_err('In: simple_comlin_smat, 1', alloc_stat)
            ! read the pairs
            funit = get_fileunit()
            allocate(fname, source='pairs_part'//int2str_pad(p%part,p%numlen)//'.bin')
            if( .not. file_exists(fname) )then
                write(*,*) 'file: ', fname, 'does not exist!'
                write(*,*) 'If all pair_part* are not in cwd, please execute simple_split_pairs to generate the required files'
                stop 'I/O error; simple_comlin_smat'
            endif
            open(unit=funit, status='OLD', action='READ', file=fname, access='STREAM')
            if( debug ) print *, 'reading pairs in range: ', p%fromp, p%top
            read(unit=funit,pos=1,iostat=io_stat) pairs(p%fromp:p%top,:)
            ! Check if the read was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(simple_volume_smat): I/O error ', io_stat, ' when reading file: ', fname
                stop 'I/O error; simple_comlin_smat'
            endif
            close(funit)
            deallocate(fname)
            ! make real-space mask if needed
            if( .not. cline%defined('lp') .and. cline%defined('msk') )then 
                call mskimg%disc(vol1%get_ldim(), p%smpd, p%msk, npix)
            endif
            ! calculate the similarities
            cnt = 0
            do ipair=p%fromp,p%top
                cnt = cnt+1
                call progress(cnt, npairs)
                ivol = pairs(ipair,1)
                jvol = pairs(ipair,2)
                call vol1%read(vollist(ivol))
                call vol2%read(vollist(jvol))
                if( cline%defined('lp') )then
                    if( cline%defined('msk') )then
                        ! apply a soft-edged mask
                        call vol1%mask(p%msk, 'soft')
                        call vol2%mask(p%msk, 'soft')
                    endif
                    corrs(ipair) = vol1%corr(vol2,lp_dyn=p%lp,hp_dyn=p%hp)
                else
                    if( cline%defined('msk') )then
                        corrs(ipair) = vol1%real_corr(vol2, mskimg)
                    else
                        corrs(ipair) = vol1%real_corr(vol2)
                    endif
                endif
            end do
            if( debug ) print *, 'did set this number of similarities: ', cnt
            ! write the similarities
            funit = get_fileunit()
            allocate(fname, source='similarities_part'//int2str_pad(p%part,p%numlen)//'.bin')
            open(unit=funit, status='REPLACE', action='WRITE', file=fname, access='STREAM')
            write(unit=funit,pos=1,iostat=io_stat) corrs(p%fromp:p%top)
            ! Check if the write was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(simple_volume_smat): I/O error ', io_stat, ' when writing to: ', fname
                stop 'I/O error; simple_comlin_smat'
            endif
            close(funit)
            deallocate(fname, corrs, pairs)
        else
            ! generate similarity matrix
            allocate(corrmat(nvols,nvols))
            corrmat = 1.
            cnt = 0
            do ivol=1,nvols-1
                do jvol=ivol+1,nvols
                    cnt = cnt+1
                    call progress(cnt, npairs)
                    call vol1%read(vollist(ivol))
                    if( debug ) write(*,*) 'read vol1: ', vollist(ivol)
                    call vol2%read(vollist(jvol))
                    if( debug ) write(*,*) 'read vol1: ', vollist(ivol)
                    corrmat(ivol,jvol) = vol1%corr(vol2,lp_dyn=p%lp,hp_dyn=p%hp)
                    if( debug ) write(*,*) 'corr ', ivol, jvol, corrmat(ivol,jvol)
                    corrmat(jvol,ivol) = corrmat(ivol,jvol)
                end do
            end do
            funit = get_fileunit()
            open(unit=funit, status='REPLACE', action='WRITE', file='vol_smat.bin', access='STREAM')
            write(unit=funit,pos=1,iostat=io_stat) corrmat
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to clin_smat.bin'
                stop 'I/O error; simple_volume_smat'
            endif
            close(funit)
            deallocate(corrmat)
        endif     
        ! end gracefully
        call simple_end('**** SIMPLE_VOLUME_SMAT NORMAL STOP ****')        
    end subroutine exec_volume_smat
    
    ! MISCELLANOUS METHODS
    
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

    subroutine exec_cluster_smat( self, cline )
        use simple_shc_cluster,   only: shc_cluster
        use simple_oris,          only: oris
        use simple_cluster_valid, only: cluster_valid
        class(cluster_smat_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)         :: p
        type(build)          :: b
        type(shc_cluster)    :: shcc
        type(cluster_valid)  :: cvalid
        real, allocatable    :: smat(:,:)
        integer              :: funit, io_stat, ncls_min, loc(1), ncls_stop, icls
        integer              :: pop, ncls, alloc_stat, irestart, ntot, cnt=0, numlen
        real                 :: avg_ratio, min_ratio, ratio, x, sim
        real, allocatable    :: validinds(:)
        integer, parameter   :: NRESTARTS=10
        logical              :: debug=.false., done=.false.
        p = params(cline,.false.)                        ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.)! general objects built
        ! obtain similarity matrix
        allocate(smat(p%nptcls,p%nptcls), stat=alloc_stat)
        call alloc_err('In: simple_cluster_smat, 1', alloc_stat)
        smat = 1.
        funit = get_fileunit()
        open(unit=funit, status='OLD', action='READ', file=p%fname, access='STREAM')
        read(unit=funit,pos=1,iostat=io_stat) smat
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when reading: ', p%fname
            stop 'I/O error; simple_cluster_smat'
        endif
        close(funit)
        allocate(validinds(2:p%ncls), stat=alloc_stat)
        call alloc_err("In: simple_cluster_smat", alloc_stat)
        validinds = 0
        ntot = (p%ncls-1)*NRESTARTS
        cnt = 0
        numlen = len(int2str(p%ncls))
        do ncls=2,p%ncls
            avg_ratio = 0.
            min_ratio = huge(x)
            do irestart=1,NRESTARTS
                cnt = cnt+1
                call progress(cnt,ntot)
                call shcc%new(p%nptcls, ncls, smat, b%a)
                call shcc%shc(.false., p%label, sim)
                call cvalid%new(b%a, ncls, p%label, smat)
                ratio = cvalid%elmlunds_ratio_index()
                avg_ratio = avg_ratio+ratio
                if( ratio < min_ratio )then
                    min_ratio = ratio
                    call b%a%write('shc_clustering_ncls'//int2str_pad(ncls,numlen)//'.txt')
                endif
            end do
            validinds(ncls) = avg_ratio/real(NRESTARTS)
        end do
        ncls_stop = 0
        done = .false. 
        do ncls=2,p%ncls
            write(*,'(a,1x,f9.3,8x,a,1x,i3)') 'COHESION/SEPARATION RATIO INDEX: ', validinds(ncls), ' NCLS: ', ncls
            call b%a%read('shc_clustering_ncls'//int2str_pad(ncls,numlen)//'.txt')
            do icls=1,ncls
                pop = b%a%get_pop(icls, p%label)
                write(*,'(a,3x,i5,1x,a,1x,i3)') '  CLUSTER POPULATION:', pop, 'CLUSTER:', icls
            end do
            write(*,'(a)') '***************************************************'
            if( ncls < p%ncls )then
                if( validinds(ncls+1) >= validinds(ncls) .and. .not. done )then
                    ncls_stop = ncls
                    done = .true.
                endif
            endif
        end do
        loc = minloc(validinds)
        ncls_min = loc(1)+1
        write(*,'(a,i3)') 'NUMBER OF CLUSTERS FOUND BY MINIMIZING ELMLUNDS RATIO INDEX: ', ncls_min
        if( ncls_stop /= 0 )then
            write(*,'(a,i3)') 'NUMBER OF CLUSTERS FOUND BY STOPPING CRITERIUM:              ', ncls_stop
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER_SMAT NORMAL STOP ****')
    end subroutine exec_cluster_smat
    
    subroutine exec_converter( self, cline )
        class(converter_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
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
            stop 'not enough arguments to execute simple_converter'
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CONVERTER NORMAL STOP ****')
    end subroutine exec_converter

    subroutine exec_ctfops( self, cline )
        use simple_procimgfile, only: apply_ctf_imgfile, apply_wiener_imgfile
        class(ctfops_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        logical      :: debug = .false.
        p = params(cline)                     ! parameters generated
        call b%build_general_tbox(p, cline)   ! general objects built
        if( cline%defined('oritab') .or. cline%defined('deftab') )then
            if( .not. b%a%isthere('dfx') ) stop 'need at least dfx to be set if ctf is going to be used!'
        else
            stop 'oritab/deftab with CTF info needed for phase flipping/multiplication'
        endif
        if( p%ctf .ne. 'no' )then
            if( debug )then
                write(*,*) 'CTF parameters used in simple_ctfops'
                write(*,*) 'kv = ', p%kv
                write(*,*) 'cs = ', p%cs
                write(*,*) 'fraca = ', p%fraca
            endif
            select case( p%ctf )
                case( 'flip' )
                    if( p%neg .eq. 'yes' )then
                        call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, b%tfun, 'flipneg')
                    else
                        call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, b%tfun, 'flip')
                    endif
                case( 'mul' )
                    if( p%neg .eq. 'yes' )then
                        call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, b%tfun, 'neg')
                    else
                        call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, b%tfun, 'ctf')
                    endif
                case( 'abs' )
                    call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, b%tfun, 'abs')
                case( 'wiener' )
                    print *, 'DOING THE WIENER THING'
                    call apply_wiener_imgfile(p%stk, p%outstk, b%a, p%smpd, b%tfun)    
                case DEFAULT
                    stop 'Unknown ctf argument'
            end select
        else if( p%ctfsq .eq. 'yes' )then
            ! APPLY CFTSQ
            if( cline%defined('bfac') )then
                call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, b%tfun, 'square', bfac=p%bfac)
            else
                call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, b%tfun, 'square')
            endif
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
        p = params(cline)                     ! parameters generated
        call b%build_general_tbox(p, cline)   ! general objects built
        if( cline%defined('stk') )then
            ! 2D
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
            if( .not.file_exists(p%vols(1)) )stop 'Cannot find input volume'
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

    subroutine exec_gatanmrc2mrc( self, cline )
        use simple_imgfile, only: imgfile
        use simple_image,   only: image
        class(gatanmrc2mrc_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)                  :: p
        type(image)                   :: frameimg
        integer                       :: ldim(3), nmrcs, funit_mrc, imrc
        integer                       :: nmovies, imovie, iframe, numlen
        character(len=STDLEN)         :: mrcfnam
        character(len=:), allocatable :: moviename
        p         = params(cline, checkdistr=.false.) ! constants & derived constants produced
        ldim      = [p%xdim,p%ydim,1]
        nmrcs     = nlines(p%filetab)
        funit_mrc = get_fileunit()
        open(unit=funit_mrc, status='old', file=p%filetab)
        call frameimg%new(ldim,p%smpd)
        if( cline%defined('nframes') )then
            if( .not. cline%defined('fbody') ) stop 'need fbody (file body of output stacks) on the command line'
            if( mod(nmrcs,p%nframes) .eq. 0 )then
                ! fine, go ahead
            else
                stop 'Number of mrc files not a multiple of nframes!'
            endif
            nmovies = nmrcs/p%nframes
            if( cline%defined('numlen') )then
                numlen = p%numlen
            else
                numlen = len(int2str(nmovies))
            endif
            imrc = 0
            do imovie=1,nmovies
                allocate(moviename, source=trim(adjustl(p%fbody))//int2str_pad(imovie, numlen)//'.mrcs')
                call del_binfile(moviename)
                do iframe=1,p%nframes
                    imrc = imrc+1
                    call progress(imrc, nmrcs)
                    read(funit_mrc,'(a256)') mrcfnam
                    call frameimg%read(mrcfnam,1,readhead=.false.)
                    call frameimg%write(moviename,iframe)
                end do
                deallocate(moviename)
            end do
        else
            do imrc=1,nmrcs 
                read(funit_mrc,'(a256)') mrcfnam
                call frameimg%read(mrcfnam,1,readhead=.false.)
                call frameimg%write(mrcfnam,1)
            end do
        endif
        close(funit_mrc)
        ! end gracefully
        call simple_end('**** SIMPLE_GATANMRC2MRC NORMAL STOP ****')
    end subroutine exec_gatanmrc2mrc
    
    subroutine exec_image_smat(self, cline)
        use simple_corrmat  ! singleton
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
    
    subroutine exec_map2ptcls( self, cline )
        use simple_oris,    only: oris
        use simple_ori,     only: ori
        use simple_image,   only: image
        use simple_corrmat  ! singleton
        class(map2ptcls_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type state_organiser
            integer, allocatable :: particles(:)
            integer              :: cls_orig = 0 
            integer              :: cls_sel  = 0
            integer              :: istate   = 0
            type(ori)            :: ori3d
        end type state_organiser
        type(params)                       :: p
        type(build)                        :: b
        type(state_organiser), allocatable :: labeler(:)
        type(image), allocatable           :: imgs_sel(:), imgs_cls(:)
        type(oris)                         :: o_comlindoc, o_state, a_copy, o_oritab2
        type(ori)                          :: ori2d, ori_comp
        real, allocatable                  :: correlations(:,:)
        integer, allocatable               :: statepops(:), state_particles(:), rejected_particles(:)
        integer                            :: isel, nsel, loc(1), iptcl, pind, icls
        integer                            :: nlines_oritab, nlines_oritab2, nlines_comlindoc, nlines_deftab
        integer                            :: cnt, istate, funit, iline, nls, lfoo(3)
        logical, allocatable               :: statedoc_exists(:), selected(:)
        character(len=STDLEN)              :: statedoc
        real                               :: corr
        logical, parameter                 :: debug=.false.
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        ! find number of selected cavgs
        call find_ldim_nptcls(p%stk2, lfoo, nsel)
        if( debug ) print *, 'nsel: ', nsel
        ! find number of original cavgs
        call find_ldim_nptcls(p%stk3, lfoo, p%ncls)
        if( debug ) print *, 'ncls: ', p%ncls
        if( p%ncls < nsel ) stop 'nr of original clusters cannot be less than the number of selected ones'
        ! find number of lines in input document
        nlines_oritab = nlines(p%oritab)
        if( debug ) print *, 'nlines_oritab: ', nlines_oritab
        if( nlines_oritab /= p%nptcls ) stop 'nr lines in oritab .ne. nr images in particle stack; must be congruent!'
        if( cline%defined('deftab') )then
            nlines_deftab = nlines(p%deftab)
            if( nlines_oritab /= nlines_deftab ) stop 'nr lines in oritab .ne. nr lines in deftab; must be congruent!'
        endif
        if( cline%defined('doclist') )then
            if( .not. cline%defined('comlindoc') )then
                if( nlines(p%doclist) /= 1 ) stop 'need a comlindoc together with statelist'        
            endif
        endif
        if( cline%defined('comlindoc') .and. cline%defined('oritab2') ) stop 'either comlindoc or oritab2 can be inputted, not both'
        if( cline%defined('doclist')   .and. cline%defined('oritab2') ) stop 'either doclist or oritab2 can be inputted, not both'
        allocate(imgs_sel(nsel), imgs_cls(p%ncls))
        ! read images
        do isel=1,nsel
            call imgs_sel(isel)%new([p%box,p%box,1], p%smpd)
            call imgs_sel(isel)%read(p%stk2, isel)
        end do
        do icls=1,p%ncls
            call imgs_cls(icls)%new([p%box,p%box,1], p%smpd)
            call imgs_cls(icls)%read(p%stk3, icls)
        end do
        write(*,'(a)') '>>> CALCULATING CORRELATIONS'
        call calc_cartesian_corrmat(imgs_sel, imgs_cls, correlations)
        ! find selected clusters & map selected to original clusters & extract the particle indices
        allocate(labeler(nsel), selected(p%ncls))
        ! initialise selection array
        selected = .false.
        write(*,'(a)') '>>> MAPPING SELECTED TO ORIGINAL CLUSTERS'
        do isel=1,nsel
            loc                     = maxloc(correlations(isel,:))
            labeler(isel)%cls_orig  = loc(1)
            selected(loc(1))        = .true.
            if( debug ) print *, 'found orig clsind: ', labeler(isel)%cls_orig
            labeler(isel)%cls_sel   = isel
            if( debug ) print *, 'selected class index: ', labeler(isel)%cls_sel
            labeler(isel)%particles = b%a%get_cls(labeler(isel)%cls_orig)
            if( debug ) print *, 'got this number of partices: ', size(labeler(isel)%particles)
        end do
        ! erase deselected (by setting their state to zero)
        do icls=1,p%ncls
            if( selected(icls) ) cycle
            if( b%a%get_clspop(icls) > 0 )then
                rejected_particles = b%a%get_cls(icls)
                do iptcl=1,size(rejected_particles)
                    call b%a%set(rejected_particles(iptcl), 'state', 0.)
                end do
                deallocate(rejected_particles)
            endif
        end do
        ! parse state info
        if( cline%defined('comlindoc') )then
            write(*,'(a)') '>>> PROCESSING COMLIN STATE ASSIGNMENT DOCUMENT'
            nlines_comlindoc = nlines(p%comlindoc)
            if( nsel /= nlines_comlindoc ) stop 'nr lines in comlindoc .ne. nr of selected clusters; must be congruent!'
            ! make a new oris object and read in the comlin clustering (state) info
            o_comlindoc = oris(nsel)
            call o_comlindoc%read(p%comlindoc)
            if( .not. o_comlindoc%isthere('state') )then
                write(*,*) 'no state labeling in comlindoc, perhaps you clustered with label=class'
                stop 'please, re-cluster with label=state'
            endif
            ! set the number of states
            p%nstates = o_comlindoc%get_nstates()
            ! map states to selected cavgs
            do isel=1,nsel
                labeler(isel)%istate = nint(o_comlindoc%get(isel, 'state'))
            end do
            ! extract the state populations
            allocate( statepops(p%nstates) )
            do istate=1,p%nstates
                statepops(istate) = o_comlindoc%get_statepop(istate)
            end do
        else
            ! set default values for nstates, statepop and state
            p%nstates         = 1
            allocate( statepops(1) )
            statepops(1)      = nsel
            labeler(:)%istate = 1
        endif
        if( cline%defined('oritab2') )then
            if( .not. file_exists(p%oritab2) ) stop 'Inputted oritab2 does not exist in the cwd'
            nlines_oritab2 = nlines(p%oritab2)
            if( nlines_oritab2 /= nsel ) stop 'Nr lines in oritab2 /= nr of selected cavgs'
            o_oritab2 = oris(nsel)
            call o_oritab2%read(p%oritab2)
            ! compose orientations and set states
            do isel=1,nsel
                ! get 3d ori
                labeler(isel)%ori3d  = o_oritab2%get_ori(isel)
                labeler(isel)%istate = nint(labeler(isel)%ori3d%get('state'))
                corr                 = labeler(isel)%ori3d%get('corr')
                do iptcl=1,size(labeler(isel)%particles)
                    ! get particle index 
                    pind = labeler(isel)%particles(iptcl)
                    ! get 2d ori
                    ori2d = b%a%get_ori(pind)
                    if( cline%defined('mul') )then
                        call ori2d%set('x', p%mul*ori2d%get('x'))
                        call ori2d%set('y', p%mul*ori2d%get('y'))
                    endif
                    ! transfer original parameters in b%a 
                    ori_comp = b%a%get_ori(pind)
                    ! compose ori3d and ori2d
                    call labeler(isel)%ori3d%compose3d2d(ori2d, ori_comp)
                    ! set parameters in b%a
                    call b%a%set_ori(pind,ori_comp)
                    call b%a%set(pind, 'corr', corr)
                end do
            end do
        endif
        ! map states to particles
        do isel=1,nsel
            do iptcl=1,size(labeler(isel)%particles)
                ! get particle index
                pind = labeler(isel)%particles(iptcl)
                call b%a%set(pind, 'state', real(labeler(isel)%istate))
            end do
        end do
        ! parse ori info
        if( cline%defined('doclist') )then
            write(*,'(a)') '>>> COMBINING 3D ORIS (CAVGS) WITH 2D ALIGNMENT (PARTICLES)'
            if( nlines(p%doclist) /= p%nstates )then
                stop 'the number of lines in doclist does not match the number of states in comlindoc'
            endif
            allocate(statedoc_exists(p%nstates))
            ! read in 3d orientations
            funit = get_fileunit()
            open(unit=funit, status='old', file=p%doclist)
            do istate=1,p%nstates
                ! read the relevant statedoc
                read(funit,'(a256)') statedoc
                statedoc_exists(istate) = file_exists(statedoc)
                if( statedoc_exists(istate) )then
                    nls = nlines(statedoc)
                    if( nls /= statepops(istate) )then
                        write(*,*) 'the nr of lines in statedoc: ', trim(statedoc),&
                        'does not match pop size: ', statepops(istate), 'in comlindoc'
                        stop
                    endif
                    o_state = oris(nls)
                    call o_state%read(statedoc)
                else
                    ! make a fake o_state
                    o_state = oris(statepops(istate))
                    do iline=1,statepops(istate)
                        call o_state%set(iline, 'state', 0.)
                    end do
                    statepops(istate) = 0
                endif
                cnt = 0
                do isel=1,nsel
                    if( labeler(isel)%istate == istate )then
                        cnt = cnt+1
                        labeler(isel)%ori3d = o_state%get_ori(cnt)
                    endif
                end do
            end do
            close(funit)
            ! wipe out the states for which no docs are provided
            do isel=1,nsel
                do iptcl=1,size(labeler(isel)%particles)
                    ! get particle index
                    pind = labeler(isel)%particles(iptcl)
                    ! get state index
                    istate = nint(b%a%get(pind, 'state'))
                    if( .not. statedoc_exists(istate) )then
                        call b%a%set(pind, 'state', 0.)
                    endif
                end do
            end do
            ! compose orientations
            do isel=1,nsel
                do iptcl=1,size(labeler(isel)%particles)
                    ! get particle index
                    pind = labeler(isel)%particles(iptcl)
                    ! get 2d ori
                    ori2d = b%a%get_ori(pind)
                    if( cline%defined('mul') )then
                        call ori2d%set('x', p%mul*ori2d%get('x'))
                        call ori2d%set('y', p%mul*ori2d%get('y'))
                    endif
                    ! transfer original parameters in b%a
                    ori_comp = b%a%get_ori(pind)
                    ! compose ori3d and ori2d
                    call labeler(isel)%ori3d%compose3d2d(ori2d, ori_comp)
                    ! set parameters in b%a
                    call b%a%set_ori(pind,ori_comp)
                    call b%a%set(pind, 'corr',  labeler(isel)%ori3d%get('corr'))
                end do
            end do
            ! relabel states in consequtive order
            if( any(statepops == 0) )then
                a_copy = b%a
                cnt    = 0
                do istate=1,p%nstates
                    if( statepops(istate) > 0 )then
                        cnt = cnt+1
                        write(*,'(a,1x,i3,1x,a,1x,i3)') '>>> THE STATE THAT WAS FORMERLY:', istate, 'IS NOW:', cnt
                        state_particles = nint(a_copy%get_ptcls_in_state(istate))
                        do iptcl=1,size(state_particles)
                            call b%a%set(state_particles(iptcl),'state',real(cnt))
                        end do
                    else
                        write(*,'(a,1x,i3,1x,a)') '>>> THE STATE THAT WAS FORMERLY:', istate, 'HAS BEEN EXCLUDED'
                    endif
                end do
            endif
        endif
        call b%a%write(p%outfile)
        call simple_end('**** SIMPLE_MAP2PTCLS NORMAL STOP ****')
    end subroutine exec_map2ptcls

    subroutine exec_norm( self, cline )
        use simple_procimgfile,   only: norm_imgfile, noise_norm_imgfile, shell_norm_imgfile
        class(norm_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(build)       :: b
        type(params)      :: p
        real, allocatable :: spec(:)
        integer           :: k
        p = params(cline)                           ! parameters generated
        call b%build_general_tbox(p, cline)         ! general objects built
        if( cline%defined('stk') .and. cline%defined('vol1') )stop 'Cannot operate on images AND volume at once'
        if( p%norm.eq.'yes' .and. p%noise_norm.eq.'yes' )stop 'Invalid normalization type'
        if( p%norm.eq.'yes' .and. p%shell_norm.eq.'yes' )stop 'Invalid normalization type'
        if( p%noise_norm.eq.'yes' .and. p%shell_norm.eq.'yes' )stop 'Invalid normalization type'
        call b%vol%new([p%box,p%box,p%box], p%smpd) ! reallocate vol (boxmatch issue)
        if( cline%defined('stk') )then
            ! 2D

            if( p%norm.eq.'yes' )then
                ! Normalization
                if( cline%defined('hfun') )then
                    call norm_imgfile(p%stk, p%outstk, hfun=p%hfun)
                else
                    call norm_imgfile(p%stk, p%outstk)
                endif
            else if( p%noise_norm.eq.'yes' )then
                ! Noise normalization
                if( cline%defined('msk') )then
                    call noise_norm_imgfile(p%stk, p%msk, p%outstk)
                else
                    stop 'need msk parameter for noise normalization'
                endif
            else if( p%shell_norm.eq.'yes' )then
                ! shell normalization
                call shell_norm_imgfile( p%stk, p%outstk )
            endif
        else if( cline%defined('vol1') )then
            ! 3D
            if( .not.file_exists(p%vols(1)) )stop 'Cannot find input volume'
            call b%vol%read(p%vols(1))
            if( p%shell_norm.eq.'yes' )then
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
   
    subroutine exec_orisops(self,cline)
        use simple_ori,  only: ori
        use simple_oris, only: oris
        use simple_math, only: normvec
        class(orisops_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(build)  :: b
        type(ori)    :: orientation
        type(oris)   :: o
        type(params) :: p
        real         :: normal(3), mind, maxd, avgd, sdevd, sumd, vard
        real         :: mind2, maxd2, avgd2, sdevd2, vard2, homo_cnt, homo_avg
        integer      :: s, i, j, cnt
        logical      :: err
        p = params(cline)
        call b%build_general_tbox(p, cline)
        if( cline%defined('oritab2') )then
            if( .not. cline%defined('oritab') ) stop 'need oritab for comparison'
            if( nlines(p%oritab) .ne. nlines(p%oritab2) )then
                stop 'inconsistent number of lines in the two oritabs!'
            endif
            o = oris(p%nptcls)
            call o%read(p%oritab2)
            call b%a%diststat(o, sumd, avgd, sdevd, mind, maxd)
            write(*,'(a,1x,f15.6)') 'SUM OF ANGULAR DISTANCE BETWEEN ORIENTATIONS  :', sumd
            write(*,'(a,1x,f15.6)') 'AVERAGE ANGULAR DISTANCE BETWEEN ORIENTATIONS :', avgd
            write(*,'(a,1x,f15.6)') 'STANDARD DEVIATION OF ANGULAR DISTANCES       :', sdevd
            write(*,'(a,1x,f15.6)') 'MINIMUM ANGULAR DISTANCE                      :', mind
            write(*,'(a,1x,f15.6)') 'MAXIMUM ANGULAR DISTANCE                      :', maxd
        else if( cline%defined('oritab') )then
            if( cline%defined('hist') )then
                call b%a%histogram(p%hist)
                goto 999
            endif
            if( p%ctfstats .eq. 'yes' )then
                call b%a%stats('ctfres', avgd, sdevd, vard, err )
                call b%a%minmax('ctfres', mind, maxd)
                write(*,'(a,1x,f8.2)') 'AVERAGE CTF RESOLUTION               :', avgd
                write(*,'(a,1x,f8.2)') 'STANDARD DEVIATION OF CTF RESOLUTION :', sdevd
                write(*,'(a,1x,f8.2)') 'MINIMUM CTF RESOLUTION (BEST)        :', mind
                write(*,'(a,1x,f8.2)') 'MAXIMUM CTF RESOLUTION (WORST)       :', maxd
                call b%a%stats('dfx', avgd, sdevd, vard, err )
                call b%a%minmax('dfx', mind, maxd)
                call b%a%stats('dfy', avgd2, sdevd2, vard2, err )
                call b%a%minmax('dfy', mind2, maxd2)
                write(*,'(a,1x,f8.2)') 'AVERAGE DF                           :', (avgd+avgd2)/2.
                write(*,'(a,1x,f8.2)') 'STANDARD DEVIATION OF DF             :', (sdevd+sdevd2)/2.
                write(*,'(a,1x,f8.2)') 'MINIMUM DF                           :', (mind+mind2)/2.
                write(*,'(a,1x,f8.2)') 'MAXIMUM DF                           :', (maxd+maxd2)/2.
                goto 999
            endif
            if( p%trsstats .eq. 'yes' )then
                call b%a%stats('x', avgd, sdevd, vard, err )
                call b%a%minmax('x', mind, maxd)
                call b%a%stats('y', avgd2, sdevd2, vard2, err )
                call b%a%minmax('y', mind2, maxd2)
                write(*,'(a,1x,f8.2)') 'AVERAGE TRS               :', (avgd+avgd2)/2.
                write(*,'(a,1x,f8.2)') 'STANDARD DEVIATION OF TRS :', (sdevd+sdevd2)/2.
                write(*,'(a,1x,f8.2)') 'MINIMUM TRS               :', (mind+mind2)/2.
                write(*,'(a,1x,f8.2)') 'MAXIMUM TRS               :', (maxd+maxd2)/2.
                goto 999
            endif
            if( p%errify .eq. 'yes' )then   ! introduce error in input orientations
                if( cline%defined('angerr') .or. cline%defined('sherr') ) call b%a%introd_alig_err(p%angerr, p%sherr)
                if( p%ctf .eq. 'yes' ) call b%a%introd_ctf_err(p%deferr)
            endif
            if( p%mirr .eq. '2d' )then      ! mirror input Eulers
                call b%a%mirror2d
            endif
            if( p%mirr .eq. '3d' )then      ! mirror input Eulers
                call b%a%mirror3d
            endif
            if( cline%defined('e1') )then ! rotate input Eulers
                call orientation%new
                call orientation%e1set(p%e1)
                call orientation%e2set(p%e2)
                call orientation%e3set(p%e3) 
                if( cline%defined('state') )then
                    do i=1,b%a%get_noris()
                        s = nint(b%a%get(i, 'state'))
                        if( s == p%state )then
                            call b%a%rot(i,orientation)
                        endif
                    end do
                else
                    call b%a%rot(orientation)
                endif
            endif
            if( cline%defined('mul') )then
                call b%a%mul_shifts(p%mul)
            endif
            if( p%zero  .eq. 'yes' ) call b%a%zero_shifts
            if( p%plot  .eq. 'yes' )then ! plot polar vectors                          
                do i=1,b%a%get_noris()
                    normal = b%a%get_normal(i)
                    write(*,'(1x,f7.2,3x,f7.2)') normal(1), normal(2)
                end do
            endif
            if( p%discrete .eq. 'yes' )then
                if( cline%defined('ndiscrete') )then
                    call b%a%discretize(p%ndiscrete)
                else
                    stop 'need ndiscrete to be defined!'
                endif
            endif
            if( cline%defined('xsh') )then
                call b%a%map3dshift22d([p%xsh,p%ysh,p%zsh])
            endif
            if( p%clustvalid .eq. 'yes' )then
                if( cline%defined('ncls') )then
                    write(*,'(a,3x,f5.1)') '>>> COHESION: ',   b%a%cohesion_norm('class',p%ncls)*100.
                    write(*,'(a,1x,f5.1)') '>>> SEPARATION: ', b%a%separation_norm('class',p%ncls)*100.
                else if( cline%defined('nstates') )then
                    write(*,'(a,3x,f5.1)') '>>> COHESION: ',   b%a%cohesion_norm('state',p%nstates)*100.
                    write(*,'(a,1x,f5.1)') '>>> SEPARATION: ', b%a%separation_norm('state',p%nstates)*100.
                else
                    stop 'need ncls/nstates as input for clustvalid'
                endif
            else if( p%clustvalid .eq. 'homo' )then
                if( cline%defined('ncls') )then
                    call b%a%homogeneity('class', p%minp, p%thres, homo_cnt, homo_avg)
                    write(*,'(a,1x,f5.1)') '>>> THIS % OF CLUSTERS CONSIDERED HOMOGENEOUS: ', homo_cnt*100.
                    write(*,'(a,1x,f5.1)') '>>> AVERAGE HOMOGENEITY:                       ', homo_avg*100.
                else if( cline%defined('nstates') )then
                    call b%a%homogeneity('state', p%minp, p%thres, homo_cnt, homo_avg)
                    write(*,'(a,1x,f5.1)') '>>> THIS % OF CLUSTERS CONSIDERED HOMOGENEOUS: ', homo_cnt*100.
                    write(*,'(a,13x,f5.1)') '>>> AVERAGE HOMOGENEITY:                      ', homo_avg*100.
                else
                    stop 'need ncls/nstates as input for clustvalid'
                endif
            endif
            if( cline%defined('nstates') )then
                call b%a%rnd_states(p%nstates)
            endif
        else ! make orientations
            if( cline%defined('ncls') )then
                o = oris(p%ncls)
                call o%spiral(p%nsym, p%eullims)
                call b%a%new(p%ncls*p%minp)
                cnt = 0
                do i=1,p%ncls
                    orientation = o%get_ori(i)
                    do j=1,p%minp
                        cnt = cnt+1
                        call b%a%set_ori(cnt, orientation)
                    end do
                end do
                if( p%zero .ne. 'yes' ) call b%a%rnd_inpls(p%trs)   
            else if( cline%defined('ndiscrete') )then
                if( p%ndiscrete > 0 )then
                    call b%a%rnd_oris_discrete(p%ndiscrete, p%nsym, p%eullims)
                endif
                call b%a%rnd_inpls(p%trs)
            else if( p%even .eq. 'yes' )then
                call b%a%spiral(p%nsym, p%eullims)
                call b%a%rnd_inpls(p%trs)
            else
                call b%a%rnd_oris(p%trs) 
            endif
            if( p%nstates > 1 ) call b%a%rnd_states(p%nstates)
            if( cline%defined('astigerr') )then
                if( p%ctf .eq. 'yes' ) call b%a%rnd_ctf(p%defocus, p%deferr, p%astigerr)
            else
                if( p%ctf .eq. 'yes' ) call b%a%rnd_ctf(p%defocus, p%deferr)
            endif
        endif
        call b%a%write(p%outfile)
        999 call simple_end('**** SIMPLE_ORISOPS NORMAL STOP ****')
    end subroutine exec_orisops
    
    subroutine exec_print_fsc(self,cline)
        use simple_math,    only: get_resolution, get_lplim
        use simple_image,   only: image
        class(print_fsc_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        !local variables
        type(params)      :: p
        type(image)       :: img
        real, allocatable :: res(:), fsc(:)
        integer           :: k
        real              :: res0143, res05
        p = params(cline) ! parameters generated
        call img%new([p%box,p%box,1], p%smpd)
        res = img%get_res() 
        fsc = file2rarr(p%fsc)
        do k=1,size(fsc) 
        write(*,'(A,1X,F6.2,1X,A,1X,F15.3)') '>>> RESOLUTION:', res(k), '>>> FSC:', fsc(k)
        end do
        ! get & print resolution
        call get_resolution(fsc, res, res05, res0143)
        write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res0143
        write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res05
        ! END GRACEFULLY
        call simple_end('**** SIMPLE_PRINT_FSC NORMAL STOP ****')

        return
    end subroutine exec_print_fsc
    
    subroutine exec_res( self, cline )
        class(res_commander), intent(inout) :: self
        class(cmdline),       intent(inout) :: cline
        type(params) :: p
        real         :: lp
        p  = params(cline)
        lp = (real(p%box-1)*p%smpd)/real(p%find)
        write(*,'(A,1X,f7.2)') '>>> LOW-PASS LIMIT:', lp
    end subroutine exec_res
    
    subroutine exec_scale( self, cline )
        use simple_procimgfile, only: resize_and_clip_imgfile, resize_imgfile, clip_imgfile
        use simple_image,       only: image
        class(scale_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        type(image)  :: vol2, img, img2
        real    :: ave, sdev, var, med
        integer :: ldim_scaled(3)
        p = params(cline)                               ! parameters generated
        call b%build_general_tbox(p, cline)             ! general objects built
        call b%vol%new([p%box,p%box,p%box], p%smpd)     ! reallocate vol (boxmatch issue)
        call img%new([p%box,p%box,1],p%smpd,p%imgkind)  ! image created
        call img2%new([p%box,p%box,1],p%smpd,p%imgkind) ! image created
        if( cline%defined('stk') .and. cline%defined('vol1') )stop 'Cannot operate on images AND volume at once'
        if( cline%defined('stk') )then
            ! 2D
            if( cline%defined('clip').and.(.not.cline%defined('newbox').or..not.cline%defined('scale')) )then
                ! Clipping
                call clip_imgfile(p%stk,p%outstk,[p%clip,p%clip,1])
            else if( cline%defined('newbox') .or. cline%defined('scale') )then
                ! Rescaling
                ldim_scaled = [p%newbox,p%newbox,1] ! dimension of scaled
                if( cline%defined('clip') )then
                    if( cline%defined('part') )then
                        p%outstk = 'outstk_part'//int2str_pad(p%part, p%numlen)//p%ext
                        call resize_and_clip_imgfile(p%stk,p%outstk,ldim_scaled,[p%clip,p%clip,1],[p%fromp,p%top])
                    else
                        call resize_and_clip_imgfile(p%stk,p%outstk,ldim_scaled,[p%clip,p%clip,1])
                    endif
                else
                    if( cline%defined('part') )then
                        p%outstk = 'outstk_part'//int2str_pad(p%part, p%numlen)//p%ext
                        call resize_imgfile(p%stk,p%outstk,ldim_scaled,[p%fromp,p%top])
                    else
                        call resize_imgfile(p%stk,p%outstk,ldim_scaled)
                    endif
                endif
            endif
        else if( cline%defined('vol1') )then
            ! 3D
            if( .not.file_exists(p%vols(1)) )stop 'Cannot find input volume'
            call b%vol%read(p%vols(1))
            if( cline%defined('newbox') .or. cline%defined('scale') )then
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
                p%box = p%newbox
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
            endif
            if( p%outvol .ne. '' )call b%vol%write(p%outvol, del_if_exists=.true.)
        else
            stop 'SIMPLE_SCALE needs input image(s) or volume!'          
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_SCALE NORMAL STOP ****')
    end subroutine exec_scale

    subroutine exec_stackops( self, cline )
        use simple_ran_tabu,    only: ran_tabu
        use simple_procimgfile, only: mirror_imgfile, neg_imgfile, acf_imgfile, frameavg_imgfile
        use simple_procimgfile, only: add_noise_imgfile, copy_imgfile,  make_avg_imgfile
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
        p = params(cline)                     ! parameters generated
        call b%build_general_tbox(p, cline)   ! general objects built
        call img%new([p%box,p%box,1],p%smpd,p%imgkind)  ! image created
        ! MIRRORING
        if( cline%defined('mirr') )then
          call mirror_imgfile( p%stk, p%outstk, p%mirr )
          goto 999
        endif
        ! RANDOM SELECTION
        if( cline%defined('nran') )then
            write(*,'(a)') '>>> RANDOMLY SELECTING IMAGES'
            allocate( pinds(p%nran), stat=alloc_stat )
            call alloc_err('In: simple_commander; stackops', alloc_stat)
            rt = ran_tabu(p%nptcls)
            call rt%ne_ran_iarr(pinds)
            do i=1,p%nran
                call progress(i, p%nran)
                call img%read(p%stk, pinds(i))
                call img%write(p%outstk, i)
            end do
            goto 999
        endif
        ! FISHING EXPEDITIONS
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
                allocate(fname, source='extracted_oris_class'//int2str_pad(p%state,2)//'.txt')
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
                allocate(fname, source='extracted_oris_class'//int2str_pad(p%state,2)//'.txt')
            endif
            call o_here%write(fname)
            goto 999
        endif
        ! INVERT CONTRAST
        if( p%neg .eq. 'yes' )then
            call neg_imgfile(p%stk, p%outstk)
            goto 999
        endif
        ! AUTO CORRELATION FUNCTION
        if( p%acf .eq. 'yes' )then
            call acf_imgfile(p%stk, p%outstk)
            goto 999
        endif
        ! CREATE FRAME AVERAGES
        if( p%frameavg > 0 )then
            call frameavg_imgfile(p%stk,p%outstk,p%frameavg)
            goto 999
        endif
        ! VISUALIZE
        if( p%vis .eq. 'yes' )then
            do i=1,p%nptcls
                call img%read(p%stk, i)
                call img%vis
            end do
            goto 999
        endif
        ! AVERAGE
        if( p%avg .eq. 'yes' )then
            call make_avg_imgfile(p%stk, p%outstk)
            goto 999
        endif
        ! ADD NOISE
        if( cline%defined('snr') )then
            call add_noise_imgfile(p%stk, p%outstk, p%snr)
            goto 999
        endif
        ! COPYING
        if( cline%defined('top') .and. .not. cline%defined('part') )then
            call copy_imgfile(p%stk, p%outstk, fromto=[p%fromp,p%top], smpd_in=p%smpd)
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
        ! default
        write(*,*)'Nothing to do!'
        ! end gracefully
     999 call simple_end('**** SIMPLE_STACKOPS NORMAL STOP ****')
    end subroutine exec_stackops

    subroutine exec_tseries_split( self, cline )
        use simple_oris, only: oris
        use simple_ori,  only: ori
        class(tseries_split_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(params)                  :: p
        type(build)                   :: b
        integer                       :: iptcl, istart, istop, cnt, chunkcnt, numlen
        type(oris)                    :: oset
        type(ori)                     :: o
        character(len=:), allocatable :: fname, oname
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        istart   = 1
        istop    = istart+p%chunksz-1
        chunkcnt = 0
        numlen   = len(int2str(p%nptcls/p%chunksz+1))
        do while( istop <= p%nptcls )
            chunkcnt = chunkcnt+1
            allocate(fname, source='stack_chunk'//int2str_pad(chunkcnt,numlen)//p%ext)
            allocate(oname, source='oris_chunk'//int2str_pad(chunkcnt,numlen)//'.txt')
            cnt = 0
            oset = oris(istop-istart+1)
            do iptcl=istart,istop
                cnt = cnt+1
                call b%img%read(p%stk, iptcl)
                call b%img%write(fname, cnt)
                o = b%a%get_ori(iptcl)
                call oset%set_ori(cnt, o)
            end do
            if( cnt == oset%get_noris() )then
                ! all ok
            else
                stop 'wrong number of oris allocated'
            endif
            call oset%write(oname)
            call oset%kill
            istart = istart+p%jumpsz
            istop  = istop+p%jumpsz
            deallocate(fname, oname)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_SPLIT NORMAL STOP ****')
    end subroutine exec_tseries_split

    ! PARALLEL PROCESSING METHODS
    
    subroutine exec_merge_algndocs( self, cline )
        use simple_oris, only: oris
        class(merge_algndocs_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(params)          :: p
        type(oris)            :: o, o_read
        integer               :: i, j, istart, istop, ptcls_per_part, nentries, leftover, cnt, nentries_all, numlen
        character(len=STDLEN) :: fname
        logical               :: here, useoritab
        p = params(cline) ! parameters generated
        useoritab = .false.
        if( cline%defined('oritab') )then
            if( file_exists(p%oritab) ) useoritab = .true.
        endif
        if( useoritab )then
            if( nlines(p%oritab) /= p%nptcls )then
                stop 'the inputted nptcls is not consistent with the nptcls in oritab!'
            endif
            ! create object for orientations
            o = oris(p%nptcls)
            ! read previous orientations
            call o%read(p%oritab)
        endif
        ptcls_per_part = p%nptcls/p%ndocs
        leftover       = p%nptcls-ptcls_per_part*p%ndocs
        istop          = 0
        numlen         = len(int2str(p%ndocs))
        do i=1,p%ndocs
            fname  = trim(adjustl(p%fbody))//int2str_pad(i,numlen)//'.txt'
            if( i == p%ndocs )then
                istart = istop+1
                istop  = p%nptcls
            else
                if( leftover == 0 )then
                    istart = istop+1;
                    istop  = istart+ptcls_per_part-1;
                else
                    istop  = i*(ptcls_per_part+1)
                    istart = istop-(ptcls_per_part+1)+1
                    leftover = leftover-1
                endif
            endif
            ! calculate the number of all entries
            nentries_all = istop-istart+1 
            ! calculate the actual number of entries
            inquire(FILE=fname, EXIST=here)
            if( here )then
                nentries = nlines(fname)
            else
                nentries = 0
            endif
            ! check if oritab is there to fill-in blanks
            if( nentries < nentries_all )then
                if( .not. useoritab )then
                    stop 'need previous oritab to fill-in blanks; simple_merge_algndocs'
                endif
            endif
            ! print partition info
            write(*,'(a,1x,i3,1x,a,1x,i6,1x,i6)') 'partition:', i, 'from/to:', istart, istop
            if( nentries > 0 )then
                o_read = oris(nentries)
                call o_read%read(fname)
            endif
            ! read
            if( useoritab )then ! read and fill-in from oritab
                cnt = 0
                do j=istart,istop
                    cnt = cnt+1
                    if( cnt <= nentries )then
                        call o%set_ori(j,o_read%get_ori(cnt))
                    else
                        exit
                    endif
                end do
            else                                ! just merge (all ptcls is there)
                call o%merge(o_read)
            endif
        end do
        call o%write(p%outfile)
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_ALGNDOCS NORMAL STOP ****')
    end subroutine exec_merge_algndocs
    
    subroutine exec_merge_similarities( self, cline )
        use simple_map_reduce, only: merge_similarities_from_parts
        class(merge_similarities_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(params)      :: p
        real, allocatable :: simmat(:,:)
        integer           :: filnum, io_stat
        p = params(cline) ! parameters generated
        simmat = merge_similarities_from_parts(p%nptcls, p%npart)
        filnum = get_fileunit()
        open(unit=filnum, status='REPLACE', action='WRITE', file='smat.bin', access='STREAM')
        write(unit=filnum,pos=1,iostat=io_stat) simmat
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to smat.bin'
            stop 'I/O error; simple_merge_similarities'
        endif
        close(filnum)
        ! end gracefully
        call simple_end('**** SIMPLE_MERGE_SIMILARITIES NORMAL STOP ****')
    end subroutine exec_merge_similarities

    subroutine exec_split_pairs( self, cline )
        use simple_map_reduce, only: split_pairs_in_parts
        class(split_pairs_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params) :: p
        p = params(cline) ! parameters generated
        call split_pairs_in_parts(p%nptcls, p%npart)
        call simple_end('**** SIMPLE_SPLIT_PAIRS NORMAL STOP ****')
    end subroutine exec_split_pairs

    subroutine exec_split( self, cline )
        use simple_map_reduce, only: splitstk_in_parts
        class(split_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(params) :: p
        p = params(cline) ! parameters generated
        call splitstk_in_parts(p%stk, p%split)
        call simple_end('**** SIMPLE_SPLIT NORMAL STOP ****')
    end subroutine exec_split

end module simple_commander
