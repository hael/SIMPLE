!> simple_commander_sim
!!
!! This class contains the set of concrete simulation commanders of the SIMPLE
!! library. This class provides the glue between the reciver (main reciever is
!! simple_exec program) and the abstract action, which is simply execute
!! (defined by the base class: simple_commander_base). Later we can use the
!! composite pattern to create MacroCommanders (or workflows)
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Cyril Reboul & Hans Elmlund 2016
!
module simple_commander_sim
use simple_defs
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_strings,        only: int2str, int2str_pad
use simple_filehandling    ! use all in there
use simple_jiffys          ! use all in there
implicit none

public :: noiseimgs_commander
public :: simimgs_commander
public :: simmovie_commander
public :: simsubtomo_commander
private
#include "simple_local_flags.inc"
type, extends(commander_base) :: noiseimgs_commander
  contains
    procedure :: execute      => exec_noiseimgs
end type noiseimgs_commander
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


contains

    subroutine exec_noiseimgs( self, cline )
        class(noiseimgs_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer :: i, cnt, ntot
        p = params(cline, .false.)          ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        cnt  = 0
        ntot = p%top-p%fromp+1
        do i=p%fromp,p%top
            cnt = cnt+1
            call progress(cnt,ntot)
            call b%img%ran
            if( cline%defined('part') )then
                call b%img%write('noiseimgs_part'//int2str_pad(p%part,p%numlen)//p%ext, cnt)
            else
                call b%img%write(p%outstk, i)
            endif
        end do
        call simple_end('**** SIMPLE_NOISEIMGS NORMAL STOP ****')
    end subroutine exec_noiseimgs
    
    subroutine exec_simimgs( self, cline )
        use simple_ori,        only: ori
        use simple_rnd,        only: ran3
        use simple_math,       only: deg2rad
        use simple_gridding,   only: prep4cgrid
        use simple_simulator,  only: simimg
        use simple_ctf,        only: ctf
        use simple_kbinterpol, only: kbinterpol
        class(simimgs_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)       :: p
        type(build)        :: b
        type(ori)          :: orientation
        type(ctf)          :: tfun
        type(kbinterpol)   :: kbwin
        real               :: snr_pink, snr_detector, bfac, bfacerr, dfx, dfy, angast
        integer            :: i, cnt, ntot
        p = params(cline, .false.)          ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        tfun  = ctf(p%smpd, p%kv, p%cs, p%fraca)
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
        else if( p%diverse .eq. 'yes' )then
            call b%a%gen_diverse
            call b%a%rnd_trs(p%trs)
        else if( .not. cline%defined('oritab') .and. p%single .eq. 'no' )then
            call b%a%rnd_oris(p%sherr, p%eullims)
        endif
        if( debug )then
            write(*,*) 'CTF parameters used'
            write(*,*) 'kv = ', p%kv
            write(*,*) 'cs = ', p%cs
            write(*,*) 'fraca = ', p%fraca
        endif
        if( .not. b%a%isthere('dfx') )then
            if( p%ctf .ne. 'no' ) call b%a%rnd_ctf(p%kv, p%cs, p%fraca, p%defocus, p%dferr, p%astigerr)
        endif
        DebugPrint  '>>> DONE GENERATING ORIENTATION/CTF PARAMETERS'
        call b%a%write(p%outfile)
        ! calculate snr:s
        snr_pink = p%snr/0.2
        snr_detector = p%snr/0.8
        DebugPrint  '>>> DONE CALCULATING SNR:S'
        ! prepare for image generation
        call b%vol%read(p%vols(1))
        call b%vol%mask(p%msk, 'soft')
        DebugPrint  '>>> DID READ VOL'
        call prep4cgrid(b%vol, b%vol_pad, p%msk, kbwin)
        call b%vol%expand_cmat
        DebugPrint  '>>> DONE PREPARING FOR IMAGE GENERATION'
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
            call b%vol_pad%fproject(orientation, b%img_pad)
            ! shift
            call b%img_pad%shift(orientation%get('x'),orientation%get('y'))
            if( cline%defined('bfac') )then
                ! calculate bfactor
                bfacerr = ran3()*p%bfacerr
                if( ran3() < 0.5 )then
                    bfac = p%bfac-bfacerr
                else
                    bfac = p%bfac+bfacerr
                endif
                call simimg(b%img_pad, orientation, tfun, p%ctf, p%snr, snr_pink, snr_detector)
            else
                call simimg(b%img_pad, orientation, tfun, p%ctf, p%snr, snr_pink, snr_detector, bfac)
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
        call b%vol_pad%kill_expanded
        ! end gracefully
        call simple_end('**** SIMPLE_SIMIMGS NORMAL STOP ****')
    end subroutine exec_simimgs

    subroutine exec_simmovie( self, cline )
        use simple_ori,         only: ori
        use simple_math,        only: deg2rad
        use simple_image,       only: image
        use simple_rnd,         only: ran3
        use simple_procimgfile, only: stats_imgfile
        use simple_ctf,         only: ctf
        class(simmovie_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(params)         :: p
        type(build)          :: b
        type(image)          :: base_image, shifted_base_image
        type(ctf)            :: tfun
        real                 :: snr_pink, snr_detector, ave, sdev, var, med, fracarea, x, y, sherr, dfx, dfy, deferr, angast
        integer              :: i, ptclarea, mgrapharea, fixed_frame, alloc_stat
        integer, allocatable :: ptcl_positions(:,:)
        real, allocatable    :: shifts(:,:)
        logical              :: here
        p = params(cline)                     ! parameters generated
        if( p%box == 0 ) stop 'box=0, something is fishy!'
        call b%build_general_tbox(p, cline)   ! general objects built
        tfun = ctf(p%smpd, p%kv, p%cs, p%fraca)
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
        DebugPrint  'inserted projections'
        ! calculate snr:s
        snr_pink     = p%snr/0.2
        snr_detector = p%snr/0.8
        DebugPrint  'calculated SNR:s'
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
        DebugPrint  'did set CTF parameters'
        DebugPrint  'initialized shifts'
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
        DebugPrint  'made stack for output movie frames'
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
                call tfun%apply(shifted_base_image, dfx, 'neg', dfy, angast, p%bfac)
            else
                call tfun%apply(shifted_base_image, dfx, 'ctf', dfy, angast, p%bfac)
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
        DebugPrint  'generated movie'
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
        DebugPrint  'generated optimal average'
        call base_image%write('optimal_movie_average'//p%ext, 1)
        if( p%vis .eq. 'yes' ) call base_image%vis
        ! output orientations
        call b%a%write('simmovie_params.txt')
        ! end gracefully
        call simple_end('**** SIMPLE_SIMMOVIE NORMAL STOP ****')

        contains

            !> \brief  generate mutually exclusive positions
            function gen_ptcl_pos( npos, xdim, ydim, box ) result( pos )
                use simple_jiffys, only: alloc_err
                use simple_rnd,    only: irnd_uni
                use simple_jiffys, only: progress
                integer, intent(in)           :: npos, xdim, ydim
                integer, intent(in), optional :: box
                integer, allocatable          :: pos(:,:)
                logical                       :: occupied(xdim,ydim)
                integer                       :: alloc_stat, ix, iy, cnt, i, j
                allocate( pos(npos,2), stat=alloc_stat )
                call alloc_err("In: gen_ptcl_pos, simple_math", alloc_stat)
                occupied = .false.
                cnt = 0
                do
                    ix = irnd_uni(xdim)
                    iy = irnd_uni(ydim)
                    if( present(box) )then
                        if( ix < box/2+1 .or. ix > xdim-box/2-1 ) cycle
                        if( iy < box/2+1 .or. iy > ydim-box/2-1 ) cycle
                        do i=ix-box/2,ix+box/2-1
                            do j=iy-box/2,iy+box/2-1
                                if( occupied(i,j) ) cycle
                            end do
                        end do
                    else
                        if(occupied(ix,iy)) cycle
                    endif
                    if( present(box) )then
                        occupied(ix-box/2:ix+box/2-1,iy-box/2:iy+box/2-1) = .true.
                    else
                        occupied(ix,iy) = .true.
                    endif
                    cnt = cnt+1
                    call progress(cnt,npos)
                    pos(cnt,1) = ix
                    pos(cnt,2) = iy
                    if( cnt == 500*npos )then
                        write(*,'(a)') "WARNING! Exiting loop because maximum nr of iterations; gen_ptcl_pos; simple_math"
                        exit 
                    endif
                    if( cnt == npos ) exit
                end do
            end function gen_ptcl_pos

    end subroutine exec_simmovie
    
    subroutine exec_simsubtomo( self, cline )
        use simple_image,          only: image
        use simple_ori,            only: ori
        use simple_projector_hlev, only: rotvol
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
            vol_rot = rotvol(b%vol, o, p)
            call vol_rot%write('subtomo'//int2str_pad(iptcl,numlen)//p%ext)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_SIMSUBTOMO NORMAL STOP ****')
    end subroutine exec_simsubtomo

end module simple_commander_sim
