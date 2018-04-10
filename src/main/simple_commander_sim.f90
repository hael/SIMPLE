! concrete commander: simulation routines

module simple_commander_sim
#include "simple_lib.f08"
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
implicit none

public :: simulate_noise_commander
public :: simulate_particles_commander
public :: simulate_movie_commander
public :: simulate_subtomogram_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: simulate_noise_commander
  contains
    procedure :: execute      => exec_simulate_noise
end type simulate_noise_commander
type, extends(commander_base) :: simulate_particles_commander
  contains
    procedure :: execute      => exec_simulate_particles
end type simulate_particles_commander
type, extends(commander_base) :: simulate_movie_commander
  contains
    procedure :: execute      => exec_simulate_movie
end type simulate_movie_commander
type, extends(commander_base) :: simulate_subtomogram_commander
  contains
    procedure :: execute      => exec_simulate_subtomogram
end type simulate_subtomogram_commander

contains

    subroutine exec_simulate_noise( self, cline )
        class(simulate_noise_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer :: i, cnt, ntot
        p = params(cline, .false.) ! parameters generated
        if( .not. cline%defined('outstk') ) p%outstk = 'simulated_noise'//p%ext
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        cnt  = 0
        ntot = p%top-p%fromp+1
        do i=p%fromp,p%top
            cnt = cnt+1
            call progress(cnt,ntot)
            call b%img%ran
            if( cline%defined('part') )then
                call b%img%write('simulate_noise_part'//int2str_pad(p%part,p%numlen)//p%ext, cnt)
            else
                call b%img%write(p%outstk, i)
            endif
        end do
        call simple_end('**** SIMPLE_SIMULATE_NOISE NORMAL STOP ****')
    end subroutine exec_simulate_noise

    subroutine exec_simulate_particles( self, cline )
        use simple_ori,        only: ori
        use simple_rnd,        only: ran3
        use simple_math,       only: deg2rad
        use simple_gridding,   only: prep4cgrid
        use simple_simulator,  only: simimg
        use simple_ctf,        only: ctf
        use simple_kbinterpol, only: kbinterpol
        use simple_projector,  only: projector
        class(simulate_particles_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)       :: p
        type(build)        :: b
        type(ori)          :: orientation
        type(ctf)          :: tfun
        type(kbinterpol)   :: kbwin
        type(projector)    :: vol_pad
        real               :: snr_pink, snr_detector, bfac, bfacerr
        integer            :: i, cnt, ntot
        debug=.false. ! declared in local flags
        p = params(cline, .false.)          ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        kbwin = kbinterpol(KBWINSZ, p%alpha)
        tfun  = ctf(p%smpd, p%kv, p%cs, p%fraca)
        if( .not. cline%defined('outstk') ) p%outstk = 'simulated_particles'//p%ext
        if( cline%defined('part') )then
            if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
        else
            if( .not. cline%defined('outfile') ) p%outfile = 'simulated_oris'//trim(TXT_EXT)
        endif
        if( p%box == 0 ) stop 'box=0, something is fishy! Perhaps forgotten to input volume or stack?'
        ! generate orientation/CTF parameters
        if( cline%defined('ndiscrete') )then
            if( p%ndiscrete > 0 )then
                call b%a%rnd_oris_discrete(p%ndiscrete, p%nsym, p%eullims)
            endif
            call b%a%rnd_inpls(p%trs)
        else if( .not. cline%defined('oritab') )then
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
        call b%a%write(p%outfile, [1,p%nptcls])
        ! calculate snr:s
        snr_pink = p%snr/0.2
        snr_detector = p%snr/0.8
        DebugPrint  '>>> DONE CALCULATING SNR:S'
        ! prepare for image generation
        call b%vol%read(p%vols(1))
        call b%vol%mask(p%msk, 'soft')
        call vol_pad%new([p%boxpd, p%boxpd, p%boxpd], p%smpd)
        DebugPrint  '>>> DID READ VOL'
        call prep4cgrid(b%vol, vol_pad, kbwin)
        call vol_pad%expand_cmat(p%alpha)
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
            call vol_pad%fproject(orientation, b%img_pad)
            ! shift
            call b%img_pad%shift([orientation%get('x'),orientation%get('y'),0.])
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
                call b%img%write('simulated_particles_part'//int2str_pad(p%part,p%numlen)//p%ext, cnt)
            else
                call b%img%write(p%outstk, i)
            endif
        end do
        call vol_pad%kill_expanded
        ! end gracefully
        call simple_end('**** SIMPLE_SIMULATE_PARTICLES NORMAL STOP ****')
    end subroutine exec_simulate_particles

    subroutine exec_simulate_movie( self, cline )
        use simple_ori,   only: ori
        use simple_math,  only: deg2rad
        use simple_image, only: image
        use simple_rnd,   only: ran3
        use simple_ctf,   only: ctf
        class(simulate_movie_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(params)         :: p
        type(build)          :: b
        type(image)          :: base_image, shifted_base_image
        type(ctf)            :: tfun
        real                 :: snr_pink, snr_detector, ave, sdev, var, med, fracarea, x, y, sherr, dfx, dfy, deferr, angast
        integer              :: i, ptclarea, mgrapharea, fixed_frame
        integer, allocatable :: ptcl_positions(:,:)
        real,    allocatable :: shifts(:,:)
        logical              :: here
        debug=.false.                           ! declared in local flags
        p = params(cline, spproj_a_seg=STK_SEG) ! parameters generated
        if( p%box == 0 ) stop 'box=0, something is fishy!'
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
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
        ! make a base image by inserting the projections at ptcl_positions
        call base_image%new([p%xdim,p%ydim,1], p%smpd)
        do i=1,p%nptcls
            call b%img%read(p%stk, i)
            call b%img%insert(ptcl_positions(i,:), base_image)
        end do
        if( p%vis .eq. 'yes' ) call base_image%vis()
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
        allocchk('In: simple_simulate_movie; shifts')
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
        call base_image%fwd_ft
        do i=1,p%nframes
            call progress(i,p%nframes)
            ! shift base image
            shifted_base_image = base_image
            call shifted_base_image%shift([shifts(i,1),shifts(i,2),0.])
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
            if( p%vis .eq. 'yes' ) call shifted_base_image%vis()
            call shifted_base_image%write('simulate_movie'//p%ext, i)
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
            call shifted_base_image%read('simulate_movie'//p%ext, i)
            x = b%a%get(1, 'x'//int2str(i))
            y = b%a%get(1, 'y'//int2str(i))
            call shifted_base_image%shift([-x,-y,0.])
            call base_image%add(shifted_base_image)
        end do
        if( debug ) write(*,'(a,1x,f7.4)') 'constant 4 division:', real(p%nframes)
        call base_image%div(real(p%nframes))
        DebugPrint  'generated optimal average'
        call base_image%write('optimal_movie_average'//p%ext, 1)
        if( p%vis .eq. 'yes' ) call base_image%vis()
        ! output orientations
        call b%a%write('simulate_movie_params'//trim(TXT_EXT), [1,b%a%get_noris()])
        ! end gracefully
        call simple_end('**** SIMPLE_SIMULATE_MOVIE NORMAL STOP ****')

        contains

            !> \brief  generate mutually exclusive positions
            function gen_ptcl_pos( npos, xdim, ydim, box ) result( pos )
                use simple_rnd,    only: irnd_uni
                integer, intent(in)           :: npos, xdim, ydim
                integer, intent(in), optional :: box
                integer, allocatable          :: pos(:,:)
                logical                       :: occupied(xdim,ydim)
                integer                       :: ix, iy, cnt, i, j
                allocate( pos(npos,2), stat=alloc_stat )
                allocchk("In: gen_ptcl_pos, simple_math")
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

    end subroutine exec_simulate_movie

    subroutine exec_simulate_subtomogram( self, cline )
        use simple_image,          only: image
        use simple_ori,            only: ori
        use simple_projector_hlev, only: rotvol
        class(simulate_subtomogram_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        type(image)  :: vol_rot
        type(ori)    :: o
        integer      :: iptcl, numlen
        p = params(cline) ! constants & derived constants produced, mode=2
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
        call simple_end('**** SIMPLE_SIMULATE_SUBTOMOGRAM NORMAL STOP ****')
    end subroutine exec_simulate_subtomogram

end module simple_commander_sim
