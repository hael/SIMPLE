! concrete commander: simulation routines
module simple_commander_sim
include 'simple_lib.f08'
use simple_parameters,     only: parameters
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_image,          only: image
use simple_ori,            only: ori
use simple_ctf,            only: ctf
use simple_simulator,      only: simimg
use simple_commander_base, only: commander_base
implicit none

public :: simulate_noise_commander
public :: simulate_particles_commander
public :: simulate_movie_commander
public :: simulate_subtomogram_commander
public :: simulate_atoms_commander
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
type, extends(commander_base) :: simulate_atoms_commander
  contains
    procedure :: execute      => exec_simulate_atoms
end type simulate_atoms_commander

contains

    subroutine exec_simulate_noise( self, cline )
        class(simulate_noise_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer :: i, cnt, ntot
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        if( .not. cline%defined('outstk') ) params%outstk = 'simulated_noise'//params%ext
        cnt  = 0
        ntot = params%top-params%fromp+1
        do i=params%fromp,params%top
            cnt = cnt+1
            call progress(cnt,ntot)
            call build%img%ran
            if( cline%defined('part') )then
                call build%img%write('simulate_noise_part'//int2str_pad(params%part,params%numlen)//params%ext, cnt)
            else
                call build%img%write(params%outstk, i)
            endif
        end do
        call simple_end('**** SIMPLE_SIMULATE_NOISE NORMAL STOP ****')
    end subroutine exec_simulate_noise

    subroutine exec_simulate_particles( self, cline )
        use simple_kbinterpol, only: kbinterpol
        use simple_projector,  only: projector
        use simple_oris,       only: oris
        class(simulate_particles_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(ori)        :: orientation
        type(oris)       :: spiral
        type(ctf)        :: tfun
        type(projector)  :: vol_pad
        real             :: snr_pink, snr_detector, bfac, bfacerr
        integer          :: i, cnt, ntot
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir', 'yes')
        call cline%set('nspace', cline%get_rarg('nptcls'))
        if( .not. cline%defined('sherr') .and. .not. cline%defined('oritab') ) call cline%set('sherr', 2.)
        if( .not. cline%defined('ctf')      ) call cline%set('ctf',    'yes')
        if( .not. cline%defined('dferr')    ) call cline%set('dferr',    1.5)
        if( .not. cline%defined('astigerr') ) call cline%set('astigerr', 0.5)
        if( .not. cline%defined('bfacerr')  ) call cline%set('bfacerr',  0.0)
        call cline%set('wfun', 'kb')
        call cline%set('winsz', 1.5)
        call cline%set('alpha',  2.)
        call cline%set('eo', '  no')
        call build%init_params_and_build_general_tbox(cline,params)
        if( .not. cline%defined('outstk') ) params%outstk = 'simulated_particles'//params%ext
        if( cline%defined('part') )then
            if( .not. cline%defined('outfile') ) THROW_HARD('need unique output file for parallel jobs')
        else
            if( .not. cline%defined('outfile') ) params%outfile = 'simulated_oris'//trim(TXT_EXT)
        endif
        if( params%box == 0 ) THROW_HARD('box=0, something is fishy! Perhaps forgotten to input volume or stack?')
        tfun = ctf(params%smpd, params%kv, params%cs, params%fraca)
        ! generate orientation/CTF parameters
        if( cline%defined('ndiscrete') )then
            if( params%ndiscrete > 0 )then
                call spiral%new(params%ndiscrete, is_ptcl=.false.)
                call build%pgrpsyms%build_refspiral(spiral)
                call build%spproj_field%rnd_oris_discrete_from(spiral)
                call spiral%kill
            endif
            call build%spproj_field%rnd_inpls(params%trs)
        else if( .not. cline%defined('oritab') )then
            call build%spproj_field%rnd_oris(params%sherr, params%eullims)
        endif
        if( .not. build%spproj_field%isthere('dfx') )then
            if( params%ctf .ne. 'no' ) call build%spproj_field%rnd_ctf(params%kv, params%cs, params%fraca, params%defocus, params%dferr, params%astigerr)
        endif
        call build%spproj_field%write(params%outfile, [1,params%nptcls])
        ! calculate snr:s
        snr_pink = params%snr/0.2
        snr_detector = params%snr/0.8
        ! prepare for image generation
        call build%vol%read(params%vols(1))
        call build%vol%mask(params%msk, 'soft')
        if( params%gridding.eq.'yes' ) call build%vol%div_w_instrfun(alpha=params%alpha)
        call vol_pad%new([params%boxpd, params%boxpd, params%boxpd], params%smpd)
        call build%vol%pad(vol_pad)
        call vol_pad%fft
        call vol_pad%expand_cmat(params%alpha)
        write(logfhandle,'(A)') '>>> GENERATING IMAGES'
        cnt = 0
        ntot = params%top-params%fromp+1
        do i=params%fromp,params%top
            cnt = cnt+1
            call progress(cnt,ntot)
            ! zero images
            build%img_pad = cmplx(0.,0.)
            build%img = 0.
            ! extract ori
            call build%spproj_field%get_ori(i, orientation)
            ! project vol
            call vol_pad%fproject(orientation, build%img_pad)
            ! shift
            call build%img_pad%shift([orientation%get('x'),orientation%get('y'),0.])
            if( cline%defined('bfac') )then
                ! calculate bfactor
                bfacerr = ran3()*params%bfacerr
                if( ran3() < 0.5 )then
                    bfac = params%bfac-bfacerr
                else
                    bfac = params%bfac+bfacerr
                endif
                call simimg(build%img_pad, orientation, tfun, params%ctf, params%snr, snr_pink, snr_detector)
            else
                call simimg(build%img_pad, orientation, tfun, params%ctf, params%snr, snr_pink, snr_detector, bfac)
            endif
            ! clip
            call build%img_pad%clip(build%img)
            ! write to stack
            if( cline%defined('part') )then
                call build%img%write('simulated_particles_part'//int2str_pad(params%part,params%numlen)//params%ext, cnt)
            else
                call build%img%write(params%outstk, i)
            endif
        end do
        call vol_pad%kill_expanded
        call orientation%kill
        ! end gracefully
        call simple_end('**** SIMLE_SIMULATE_PARTICLES NORMAL STOP ****')
    end subroutine exec_simulate_particles

    subroutine exec_simulate_movie( self, cline )
        class(simulate_movie_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)     :: params
        type(builder)        :: build
        type(image)          :: base_image, shifted_base_image
        type(ctf)            :: tfun
        real                 :: snr_pink, snr_detector, fracarea, x, y, sherr, dfx, dfy, deferr, angast
        integer              :: i, ptclarea, mgrapharea, fixed_frame
        integer, allocatable :: ptcl_positions(:,:)
        real,    allocatable :: shifts(:,:)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('trs')     ) call cline%set('trs',      3.)
        if( .not. cline%defined('ctf')     ) call cline%set('ctf',   'yes')
        if( .not. cline%defined('bfac')    ) call cline%set('bfac',   200.)
        if( .not. cline%defined('nframes') ) call cline%set('nframes', 30.)
        call cline%set('oritype', 'stk')
        call cline%set('wfun', 'kb')
        call cline%set('winsz', 1.5)
        call cline%set('alpha',  2.)
        call cline%set('eo',   'no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        tfun = ctf(params%smpd, params%kv, params%cs, params%fraca)
        ! set fixed frame
        fixed_frame = nint(real(params%nframes)/2.)
        ! remake the alignment doc
        call build%spproj_field%new(1, is_ptcl=.false.)
        ! check the fractional area occupied by particles & generate particle positions
        ptclarea   = params%box*params%box*params%nptcls
        mgrapharea = params%xdim*params%ydim
        fracarea   = real(ptclarea)/real(mgrapharea)
        write(logfhandle,'(a,1x,f7.3)') 'Fraction of area occupied by ptcls:', fracarea
        if( fracarea > 0.55 )then
            write(logfhandle,'(A)') 'It is not recommended that more than 55% of the micrograph area is occupied with particles!'
            write(logfhandle,'(A)') 'Please, reduce the number of projection images to place!'
            stop
        endif
        write(logfhandle,'(a)') '>>> GENERATING PARTICLE POSITIONS'
        ptcl_positions = gen_ptcl_pos(params%nptcls, params%xdim, params%ydim, params%box)
        ! make a base image by inserting the projections at ptcl_positions
        call base_image%new([params%xdim,params%ydim,1], params%smpd)
        do i=1,params%nptcls
            call build%img%read(params%stk, i)
            call build%img%insert(ptcl_positions(i,:), base_image)
        end do
        if( params%vis .eq. 'yes' ) call base_image%vis()
        ! calculate snr:s
        snr_pink     = params%snr/0.2
        snr_detector = params%snr/0.8
        ! set CTF parameters
        deferr = ran3()*0.2
        if( ran3() < 0.5 )then
            dfx = params%defocus-deferr
        else
            dfx = params%defocus+deferr
        endif
        deferr = ran3()*0.2
        if( ran3() < 0.5 )then
            dfy = params%defocus-deferr
        else
            dfy = params%defocus+deferr
        endif
        angast = ran3()*360.
        ! generate shifts
        allocate( shifts(params%nframes,2))
        x = 0.
        y = 0.
        do i=1,params%nframes
            ! generate random shifts
            sherr = ran3()*params%trs
            if( ran3() > 0.5 )then
                x = x+sherr
            else
                x = x-sherr
            endif
            sherr = ran3()*params%trs
            if( ran3() > 0.5 )then
                y = y+sherr
            else
                y = y-sherr
            endif
            shifts(i,1) = x
            shifts(i,2) = y
        end do
        ! put the central frame in the series at x,y = (0,0) & fill-up build%spproj_field
        do i=1,params%nframes
            shifts(i,:) = shifts(i,:)-shifts(fixed_frame,:)
            call build%spproj_field%set(1, 'x'//int2str(i), shifts(i,1))
            call build%spproj_field%set(1, 'y'//int2str(i), shifts(i,2))
        end do
        ! make and open a stack for the movie frames
        write(logfhandle,'(a)') '>>> GENERATING MOVIE FRAMES'
        call base_image%fft()
        do i=1,params%nframes
            call progress(i,params%nframes)
            ! shift base image
            shifted_base_image = base_image
            call shifted_base_image%shift([shifts(i,1),shifts(i,2),0.])
            call shifted_base_image%ifft()
            ! add pink noise
            call shifted_base_image%add_gauran(snr_pink)
            ! multiply with CTF
            call shifted_base_image%fft()
            if( params%neg .eq. 'yes' )then
                call tfun%apply(shifted_base_image, dfx, 'neg', dfy, angast, params%bfac)
            else
                call tfun%apply(shifted_base_image, dfx, 'ctf', dfy, angast, params%bfac)
            endif
            call shifted_base_image%ifft()
            ! add the detector noise
            call shifted_base_image%add_gauran(snr_detector)
            if( params%vis .eq. 'yes' ) call shifted_base_image%vis()
            call shifted_base_image%write('simulate_movie'//params%ext, i)
            ! set orientation parameters in object
            call build%spproj_field%set(1, 'dfx',    dfx)
            call build%spproj_field%set(1, 'dfy',    dfy)
            call build%spproj_field%set(1, 'angast', angast)
        end do
        ! generate the optimal average
        base_image = 0.
        write(logfhandle,'(a)') '>>> GENERATING OPTIMAL AVERAGE'
        do i=1,params%nframes
            call progress(i,params%nframes)
            call shifted_base_image%read('simulate_movie'//params%ext, i)
            x = build%spproj_field%get(1, 'x'//int2str(i))
            y = build%spproj_field%get(1, 'y'//int2str(i))
            call shifted_base_image%shift([-x,-y,0.])
            call base_image%add(shifted_base_image)
        end do
        call base_image%div(real(params%nframes))
        call base_image%write('optimal_movie_average'//params%ext, 1)
        if( params%vis .eq. 'yes' ) call base_image%vis()
        ! output orientations
        call build%spproj_field%write('simulate_movie_params'//trim(TXT_EXT), [1,build%spproj_field%get_noris()])
        ! end gracefully
        call simple_end('**** SIMPLE_SIMULATE_MOVIE NORMAL STOP ****')

        contains

            !> \brief  generate mutually exclusive positions
            function gen_ptcl_pos( npos, xdim, ydim, box ) result( pos )
                integer, intent(in)           :: npos, xdim, ydim
                integer, intent(in), optional :: box
                integer, allocatable          :: pos(:,:)
                logical                       :: occupied(xdim,ydim)
                integer                       :: ix, iy, cnt, i, j
                allocate( pos(npos,2) )
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
                        THROW_WARN('exiting loop because maximum nr of iterations; gen_ptcl_pos')
                        exit
                    endif
                    if( cnt == npos ) exit
                end do
            end function gen_ptcl_pos

    end subroutine exec_simulate_movie

    subroutine exec_simulate_subtomogram( self, cline )
        use simple_projector_hlev, only: rotvol
        class(simulate_subtomogram_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(image)      :: vol_rot
        type(ori)        :: o
        integer          :: iptcl, numlen
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call build%init_params_and_build_general_tbox(cline,params)
        call vol_rot%new([params%box,params%box,params%box], params%smpd)
        call build%vol%new([params%box,params%box,params%box],   params%smpd)
        call build%vol%read(params%vols(1))
        call build%vol%add_gauran(params%snr)
        call o%new(is_ptcl=.false.)
        numlen = len(int2str(params%nptcls))
        do iptcl=1,params%nptcls
            call o%rnd_ori
            vol_rot = rotvol(build%vol, o)
            call vol_rot%write('subtomo'//int2str_pad(iptcl,numlen)//params%ext)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_SIMULATE_SUBTOMOGRAM NORMAL STOP ****')
    end subroutine exec_simulate_subtomogram

    subroutine exec_simulate_atoms( self, cline )
        use simple_atoms, only: atoms
        class(simulate_atoms_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        real, parameter  :: NPIX_CUTOFF = 12.
        type(parameters) :: params
        type(image)      :: vol
        type(atoms)      :: atoms_obj
        character(len=2) :: el1, el2, el_ucase
        character(len=8) :: crystal_system
        real             :: center(3),ha,x,x1,x2,x3,y,y1,y2,y3,z,z1,z2,z3,msksq,cutoff
        real             :: a ! lattice parameter
        integer          :: i, j, k, n, ncubes
        call params%new(cline)
        if( .not.is_even(params%box) ) THROW_HARD('BOX must be even')
        call vol%new([params%box,params%box,params%box], params%smpd)
        if( cline%defined('pdbfile') )then
            if( .not.file_exists(params%pdbfile) ) THROW_HARD('Could not find: '//trim(params%pdbfile))
            call atoms_obj%new(params%pdbfile)
            cutoff = NPIX_CUTOFF * params%smpd
        else if( cline%defined('element') )then
            if( .not. cline%defined('moldiam') ) THROW_HARD('MOLDIAM must be provided for FCC lattice!')
            if( .not. atoms_obj%element_exists(params%element) ) THROW_HARD('Unsupported element for now')
            el1      = params%element(1:2)
            el2      = params%element(3:4)
            msksq    = (params%moldiam / 2.)**2.
            el_ucase = uppercase(trim(adjustl(params%element)))
            call get_lattice_params(el_ucase, crystal_system, a)
            ha       = a / 2.
            ncubes   = floor(real(params%box) * params%smpd / a)
            ! atoms at edges
            n = 0
            center = (real([params%box,params%box,params%box]/2 + 1) - 1.) * params%smpd
            select case( trim(crystal_system) )
                case('fcc')
                    ! count
                    do i = 1, ncubes
                        x  = real(i - 1) * a
                        x1 = x + ha; x2 = x; x3 = x + ha
                        do j = 1, ncubes
                            y  = real(j - 1) * a
                            y1 = y + ha; y2 = y + ha; y3 = y
                            do k = 1, ncubes
                                z  = real(k - 1) * a
                                z1 = z; z2 = z + ha; z3 = z + ha
                                ! edge
                                if( sum(([x ,y ,z ] - center)**2.) <= msksq ) n = n+1
                                ! faces
                                if( sum(([x1,y1,z1] - center)**2.) <= msksq ) n = n+1
                                if( sum(([x2,y2,z2] - center)**2.) <= msksq ) n = n+1
                                if( sum(([x3,y3,z3] - center)**2.) <= msksq ) n = n+1
                            enddo
                        enddo
                    enddo
                    ! generate atoms object
                    call atoms_obj%new(n)
                    do i = 1, n
                        call atoms_obj%set_name(i,params%element//'  ')
                        call atoms_obj%set_element(i,el1)
                    enddo
                    ! generate all coordinates
                    n = 0
                    do i = 1, ncubes
                        x  = real(i - 1) * a ;x1 = x + ha; x2 = x; x3 = x + ha
                        do j = 1, ncubes
                            y  = real(j - 1) * a; y1 = y + ha; y2 = y + ha; y3 = y
                            do k = 1, ncubes
                                z  = real(k - 1) * a; z1 = z; z2 = z + ha; z3 = z + ha
                                ! edge
                                if( sum(([x ,y ,z ] - center)**2.) <= msksq )then
                                     n = n + 1; call atoms_obj%set_coord(n, [x,y,z])
                                endif
                                ! faces
                                if( sum(([x1,y1,z1] - center)**2.) <= msksq )then
                                    n = n + 1; call atoms_obj%set_coord(n, [x1,y1,z1])
                                endif
                                if( sum(([x2,y2,z2] - center)**2.) <= msksq )then
                                    n = n + 1; call atoms_obj%set_coord(n, [x2,y2,z2])
                                endif
                                if( sum(([x3,y3,z3] - center)**2.) <= msksq )then
                                    n = n + 1; call atoms_obj%set_coord(n, [x3,y3,z3])
                                endif
                            enddo
                        enddo
                    enddo
                case('bcc')
                    ! count
                    do i = 1, ncubes
                        x  = real(i - 1) * a; x1 = x + ha
                        do j = 1, ncubes
                            y  = real(j - 1) * a; y1 = y + ha
                            do k= 1, ncubes
                                z  = real(k - 1) * a; z1 = z + ha
                                ! edge
                                if( sum(([x ,y ,z ] - center)**2.) <= msksq ) n = n + 1
                                ! body-centered
                                if( sum(([x1,y1,z1] - center)**2.) <= msksq ) n = n + 1
                            enddo
                        enddo
                    enddo
                    ! generate atoms object
                    call atoms_obj%new(n)
                    do i = 1, n
                        call atoms_obj%set_name(i,params%element//'  ')
                        call atoms_obj%set_element(i,el1)
                    enddo
                    ! generate all coordinates
                    n = 0
                    do i = 1, ncubes
                        x  = real(i - 1) * a; x1 = x + ha
                        do j = 1, ncubes
                            y  = real(j - 1) * a; y1 = y + ha
                            do k = 1, ncubes
                                z  = real(k - 1) * a; z1 = z + ha
                                if( sum(([x ,y ,z ] - center)**2.) <= msksq )then
                                    n = n + 1; call atoms_obj%set_coord(n, [x,y,z])
                                endif
                                if( sum(([x1,y1,z1] - center)**2.) <= msksq )then
                                    n = n + 1; call atoms_obj%set_coord(n, [x1,y1,z1])
                                endif
                            enddo
                        enddo
                    enddo
                case('rocksalt')
                    ! count
                    do i = 1, ncubes
                        x  = real(i - 1) * a; x1 = x + ha
                        do j = 1, ncubes
                            y  = real(j - 1) * a; y1 = y + ha
                            do k=1,ncubes
                                z  = real(k - 1) * a; z1 = z + ha
                                if( sum(([x ,y ,z ] - center)**2.) <= msksq ) n = n + 1
                                if( sum(([x1,y1,z ] - center)**2.) <= msksq ) n = n + 1
                                if( sum(([x ,y1,z1] - center)**2.) <= msksq ) n = n + 1
                                if( sum(([x1,y ,z1] - center)**2.) <= msksq ) n = n + 1
                                if( sum(([x1,y1,z1] - center)**2.) <= msksq ) n = n + 1
                                if( sum(([x1,y ,z ] - center)**2.) <= msksq ) n = n + 1
                                if( sum(([x ,y1,z ] - center)**2.) <= msksq ) n = n + 1
                                if( sum(([x ,y ,z1] - center)**2.) <= msksq ) n = n + 1
                            enddo
                        enddo
                    enddo
                    ! generate all coordinates
                    call atoms_obj%new(n)
                    n = 0
                    do i = 1, ncubes
                        x  = real(i - 1) * a; x1 = x + ha
                        do j = 1, ncubes
                            y  = real(j - 1) * a; y1 = y + ha
                            do k = 1, ncubes
                                z  = real(k - 1) * a; z1 = z + ha
                                ! element 1
                                if( sum(([x ,y ,z ] - center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, [x,y,z])
                                    call atoms_obj%set_name(n,el1//'  ')
                                    call atoms_obj%set_element(n,el1)
                                endif
                                if( sum(([x1 ,y1 ,z ] - center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, [x1,y1,z])
                                    call atoms_obj%set_name(n,el1//'  ')
                                    call atoms_obj%set_element(n,el1)
                                endif
                                if( sum(([x ,y1 ,z1 ] - center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, [x,y1,z1])
                                    call atoms_obj%set_name(n,el1//'  ')
                                    call atoms_obj%set_element(n,el1)
                                endif
                                if( sum(([x1 ,y ,z1 ] - center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, [x1,y,z1])
                                    call atoms_obj%set_name(n,el1//'  ')
                                    call atoms_obj%set_element(n,el1)
                                endif
                                ! element 2
                                if( sum(([x1 ,y1 ,z1] - center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, [x1 ,y1 ,z1])
                                    call atoms_obj%set_name(n,el2//'  ')
                                    call atoms_obj%set_element(n,el2)
                                endif
                                if( sum(([x1 ,y ,z] - center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, [x1 ,y ,z])
                                    call atoms_obj%set_name(n,el2//'  ')
                                    call atoms_obj%set_element(n,el2)
                                endif
                                if( sum(([x ,y1 ,z] - center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, [x ,y1 ,z])
                                    call atoms_obj%set_name(n,el2//'  ')
                                    call atoms_obj%set_element(n,el2)
                                endif
                                if( sum(([x ,y ,z1] - center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, [x ,y ,z1])
                                    call atoms_obj%set_name(n,el2//'  ')
                                    call atoms_obj%set_element(n,el2)
                                endif
                            enddo
                        enddo
                    enddo
            end select
            ! write PDB
            do i=1,n
                call atoms_obj%set_num(i,i)
                call atoms_obj%set_resnum(i,i)
                call atoms_obj%set_chain(i,'A')
                call atoms_obj%set_occupancy(i,1.)
            enddo
            call atoms_obj%writePDB(params%outvol)
            ! for convolution
            cutoff = 3.*a
        else
            THROW_HARD('A PDB file or ELEMENT & MOLDIAM must be defined!')
        endif
        if( cline%defined('lp') )then
            cutoff = max(cutoff,4.*params%lp)
            call atoms_obj%convolve(vol,cutoff,lp=params%lp)
        else
            call atoms_obj%convolve(vol,cutoff)
        endif
        call vol%write(params%outvol)
        ! cleanup
        call atoms_obj%kill
        call vol%kill
        ! end gracefully
        call simple_end('**** SIMPLE_SIMULATE_FCC NORMAL STOP ****')
    end subroutine exec_simulate_atoms

end module simple_commander_sim
