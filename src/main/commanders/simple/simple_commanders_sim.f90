!@descr: for simulating noise, particles, movies, atoms, etc.
module simple_commanders_sim
use simple_commanders_api
use simple_defs_atoms
use simple_ctf,       only: ctf
use simple_simulator, only: simimg, simimg_wave
implicit none
#include "simple_local_flags.inc"

private :: state_fractions, draw_states, oritab_holds_states, simulate_particles_wave

type, extends(commander_base) :: commander_simulate_noise
  contains
    procedure :: execute      => exec_simulate_noise
end type commander_simulate_noise

type, extends(commander_base) :: commander_simulate_particles
  contains
    procedure :: execute      => exec_simulate_particles
end type commander_simulate_particles

type, extends(commander_base) :: commander_simulate_scattering
  contains
    procedure :: execute      => exec_simulate_scattering
end type commander_simulate_scattering

type, extends(commander_base) :: commander_simulate_movie
  contains
    procedure :: execute      => exec_simulate_movie
end type commander_simulate_movie

type, extends(commander_base) :: commander_simulate_subtomogram
  contains
    procedure :: execute      => exec_simulate_subtomogram
end type commander_simulate_subtomogram

type, extends(commander_base) :: commander_simulate_nanoparticle
  contains
    procedure :: execute      => exec_simulate_nanoparticle
end type commander_simulate_nanoparticle

contains

    subroutine exec_simulate_noise( self, cline )
        class(commander_simulate_noise), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer :: i, cnt, ntot
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        if( .not. cline%defined('outstk') ) params%outstk = string('simulated_noise')//params%ext
        cnt  = 0
        ntot = params%top-params%fromp+1
        do i=params%fromp,params%top
            cnt = cnt+1
            call progress(cnt,ntot)
            call build%img%ran
            if( cline%defined('part') )then
                call build%img%write(string('simulate_noise_part')//int2str_pad(params%part,params%numlen)//params%ext%to_char(), cnt)
            else
                call build%img%write(params%outstk, i)
            endif
        end do
        call simple_end('**** SIMPLE_SIMULATE_NOISE NORMAL STOP ****')
    end subroutine exec_simulate_noise

    subroutine exec_simulate_particles( self, cline )
        use simple_projector,           only: projector
        use simple_matcher_ptcl_io,     only: prepimgbatch, killimgbatch
        use simple_memoize_ft_maps,     only: memoize_ft_maps
        class(commander_simulate_particles), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(ori)        :: orientation, o_thres, o_mirr
        type(oris)       :: packed_oris
        type(ctf)        :: tfun
        type(projector), allocatable :: vol_pads(:)
        integer,         allocatable :: states(:)
        real             :: bfac, bfacerr, euls(3), euls_mirr(3)
        integer          :: i, cnt, ntot, nstates, istate, ldim_vol(3), nsec
        logical          :: apply_ctf, l_oritab_states, l_wave
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir',   'no')
        ! wave mode: the keys and defaults of reconstruct3D_wave
        l_wave = cline%defined('wave_model')
        if( l_wave )then
            if( .not. cline%defined('total_dose') ) THROW_HARD('wave_model needs total_dose, the dose of the scale link and the shot noise')
            if( .not. cline%defined('ewald') )then
                if( cline%get_carg('wave_model') .eq. 'projection' )then
                    call cline%set('ewald', 'no')
                else
                    call cline%set('ewald', 'neg')
                endif
            endif
        else
            if( cline%defined('ewald') .or. cline%defined('total_dose') ) THROW_HARD('ewald and total_dose belong to wave_model')
        endif
        if( .not. cline%defined('ctf')      ) call cline%set('ctf',    'yes')
        if( .not. cline%defined('dferr')    ) call cline%set('dferr',    1.5)
        if( .not. cline%defined('astigerr') ) call cline%set('astigerr', 0.5)
        if( .not. cline%defined('bfacerr')  ) call cline%set('bfacerr',  0.0)
        call cline%set('nspace', cline%get_iarg('nptcls'))
        call cline%set('eo', '  no')
        call build%init_params_and_build_general_tbox(cline,params)
        ! the simulator multiplies by the CTF or leaves it out: no other mode is applied
        select case(trim(params%ctf))
            case('yes', 'no')
            case DEFAULT
                THROW_HARD('simulate_particles applies ctf=yes or ctf=no only')
        end select
        apply_ctf     = trim(params%ctf) .eq. 'yes'
        if( l_wave .and. .not. apply_ctf ) THROW_HARD('wave_model needs ctf=yes: a phase object in focus has almost no contrast')
        if( .not. cline%defined('outstk') ) params%outstk = string('simulated_particles')//params%ext
        if( cline%defined('part') )then
            if( .not. cline%defined('outfile') ) THROW_HARD('need unique output file for parallel jobs')
        else
            if( .not. cline%defined('outfile') ) params%outfile = 'simulated_oris'//trim(TXT_EXT)
        endif
        if( params%box == 0 ) THROW_HARD('box=0, something is fishy! Perhaps forgotten to input volume or stack?')
        tfun = ctf(params%smpd, params%kv, params%cs, params%fraca)
        ! generate orientation/CTF parameters
        if( params%even.eq.'yes' )then
            call build%spproj_field%spiral
            call build%spproj_field%rnd_inpls(params%trs)
        else if( .not. cline%defined('oritab') )then
            call build%spproj_field%rnd_oris(params%sherr, params%eullims)
        endif
        if( .not. build%spproj_field%isthere('dfx') )then
            if( apply_ctf ) call build%spproj_field%rnd_ctf(params%kv, params%cs, params%fraca, params%defocus, params%dferr, params%astigerr)
        endif
        ! the state of every particle: from oritab, or drawn with the population fractions
        nstates = max(1, params%nstates)
        do istate = 2, nstates
            if( .not. cline%defined('vol'//int2str(istate)) ) THROW_HARD('nstates=K needs vol1 to volK')
        enddo
        l_oritab_states = .false.
        if( nstates > 1 .and. cline%defined('oritab') )then
            if( fname2format(params%oritab) == 'O' )then
                l_oritab_states = build%spproj_field%isthere('state')
            else
                l_oritab_states = oritab_holds_states(params%oritab)
            endif
        endif
        ! the states come from one source: the populations are not applied to states the oritab holds
        if( l_oritab_states .and. cline%defined('state_pops') )&
            &THROW_HARD('state_pops with an oritab that holds states: the states would come from the oritab')
        if( nstates == 1 )then
            call build%spproj_field%set_all2single('state', 1.0)
        else if( l_oritab_states )then
            states = nint(build%spproj_field%get_all('state'))
            if( any(states < 1 .or. states > nstates) ) THROW_HARD('oritab holds a state outside 1 to nstates')
        else
            states = draw_states(build%spproj_field%get_noris(), state_fractions(params%state_pops, nstates))
            do i = 1, size(states)
                call build%spproj_field%set_state(i, states(i))
            enddo
        endif
        ! every state of the mixture holds a particle
        if( nstates > 1 )then
            do istate = 1, nstates
                if( count(states == istate) == 0 ) THROW_HARD('state '//int2str(istate)//' of the mixture holds no particle')
            enddo
        endif
        if( trim(params%recthres).eq.'yes' )then
            euls = [params%e1, params%e2, params%e3]
            call euler_mirror(euls, euls_mirr)
            call o_thres%new(is_ptcl=.true.)
            call o_thres%set_euler(euls)
            call o_mirr%new(is_ptcl=.true.)
            call o_mirr%set_euler(euls_mirr)
        endif
        if( l_wave )then
            call simulate_particles_wave(params, build, cline, o_thres, o_mirr, cnt)
        else
            ! prepare for image generation: one padded volume per state, held together
            allocate(vol_pads(nstates))
            do istate = 1, nstates
                call find_ldim_nptcls(params%vols(istate), ldim_vol, nsec)
                if( any(ldim_vol /= [params%box, params%box, params%box]) ) THROW_HARD('the volumes of the states differ in box')
                call build%vol%read(params%vols(istate))
                call build%vol%mask3D_soft(params%msk, backgr=0.)
                call vol_pads(istate)%new([params%boxpd, params%boxpd, params%boxpd], params%smpd)
                call build%vol%pad(vol_pads(istate))
                call vol_pads(istate)%fft
                call vol_pads(istate)%expand_cmat()
            enddo
            write(logfhandle,'(A)') '>>> GENERATING IMAGES'
            cnt  = 0
            ntot = params%top-params%fromp+1
            call memoize_ft_maps(build%img_pad%get_ldim(),build%img_pad%get_smpd())
            do i=params%fromp,params%top
                ! extract ori
                call build%spproj_field%get_ori(i, orientation)
                if( trim(params%recthres).eq.'yes' )then
                    if( rad2deg(orientation .euldist. o_thres) < params%rec_athres .or.&
                       &rad2deg(orientation .euldist. o_mirr ) < params%rec_athres ) cycle
                endif
                cnt = cnt+1
                call progress(cnt,ntot)
                ! project vol
                call build%img_pad%zero_and_flag_ft
                call vol_pads(orientation%get_state())%fproject(orientation, build%img_pad)
                ! shift
                call build%img_pad%shift([orientation%get('x'),orientation%get('y'),0.])
                call build%img%zero_and_unflag_ft
                if( cline%defined('bfac') )then
                    ! calculate bfactor
                    bfacerr = ran3()*params%bfacerr
                    if( ran3() < 0.5 )then
                        bfac = params%bfac-bfacerr
                    else
                        bfac = params%bfac+bfacerr
                    endif
                    call simimg(build%img_pad, orientation, tfun, params%ctf, params%snr, bfac=bfac, apply_ctf=apply_ctf)
                else
                    call simimg(build%img_pad, orientation, tfun, params%ctf, params%snr, apply_ctf=apply_ctf)
                endif
                ! clip
                call build%img_pad%clip(build%img)
                ! write to stack
                if( cline%defined('part') )then
                    call build%img%write(string('simulated_particles_part')//int2str_pad(params%part,params%numlen)//params%ext%to_char(), cnt)
                else
                    call build%img%write(params%outstk, cnt)
                endif
            end do
            do istate = 1, nstates
                call vol_pads(istate)%kill_expanded
                call vol_pads(istate)%kill
            enddo
            deallocate(vol_pads)
        endif
        if( trim(params%recthres).eq.'yes' )then
            call packed_oris%new(cnt, is_ptcl=.true.)
            cnt = 0
            do i=params%fromp,params%top
                ! extract ori
                call build%spproj_field%get_ori(i, orientation)
                if( rad2deg(orientation .euldist. o_thres) < params%rec_athres .or.&
                   &rad2deg(orientation .euldist. o_mirr ) < params%rec_athres ) cycle
                cnt = cnt + 1
                call packed_oris%set_ori(cnt, orientation)
            enddo
            call packed_oris%write(params%outfile, [1,cnt])
        else
            call build%spproj_field%write(params%outfile, [1,params%nptcls])
        endif
        call killimgbatch(build)
        call orientation%kill
        ! end gracefully
        call simple_end('**** SIMPLE_SIMULATE_PARTICLES NORMAL STOP ****')
    end subroutine exec_simulate_particles

    subroutine exec_simulate_scattering( self, cline )
        use simple_atoms,                only: atoms
        use simple_wave_imaging,         only: WAVE_PROJECTION, WAVE_FIRST_BORN, WAVE_RYTOV, WAVE_MSLICE, detect_intensity
        use simple_wave_projector,       only: wave_projector
        class(commander_simulate_scattering), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        integer, parameter :: BATCH = 64
        type(parameters) :: params
        type(builder)    :: build
        type(atoms)      :: molecule
        type(ori)        :: orientation
        type(oris)       :: selection
        type(ctfparams)  :: ctfparms
        class(wave_projector),       allocatable :: wproj
        type(string)     :: pdb_arg, pdb_source, potential_file
        complex, allocatable :: planes(:,:,:)
        integer          :: i, ib, nb, ldim(3), nptcls_arg, nsec, vol_dim(3), model, kappa, R
        real             :: smpd
        logical          :: from_pdb, has_angast, has_dfy, has_phshift
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir', 'no')
        if( .not. cline%defined('wave_model') ) call cline%set('wave_model', 'multislice')
        if( .not. cline%defined('ewald')      )then
            ! projection samples the central section; the curved models default to specter's negative sign
            if( cline%get_carg('wave_model') .eq. 'projection' )then
                call cline%set('ewald', 'no')
            else
                call cline%set('ewald', 'neg')
            endif
        endif
        if( .not. cline%defined('bfac')       ) call cline%set('bfac', 64.0)
        if( .not. cline%defined('total_dose') ) call cline%set('total_dose', 20.0)
        if( .not. cline%defined('noise_sdev') ) call cline%set('noise_sdev', 0.0)
        if( .not. cline%defined('dferr')      ) call cline%set('dferr', 0.0)
        if( .not. cline%defined('astigerr')   ) call cline%set('astigerr', 0.0)
        if( .not. cline%defined('outstk')     ) call cline%set('outstk', 'scattering_images.mrcs')
        if( .not. cline%defined('outfile')    ) call cline%set('outfile', 'scattering_oris.txt')
        from_pdb = cline%defined('pdbfile')
        if( from_pdb .eqv. cline%defined('vol1') ) &
            &THROW_HARD('simulate_scattering needs exactly one of pdbfile and vol1')
        if( from_pdb )then
            pdb_arg = cline%get_carg('pdbfile')
            pdb_source = string(pdb_arg%to_char())
            write(logfhandle,'(A)') '>>> CONVERTING PDB POTENTIAL: '//pdb_source%to_char()
            smpd = 1.3
            if( cline%defined('smpd') ) smpd = cline%get_rarg('smpd')
            if( cline%defined('outvol') )then
                potential_file = simple_abspath(cline%get_carg('outvol'), check_exists=.false.)
            else
                potential_file = simple_abspath(string('scattering_potential.mrc'), check_exists=.false.)
            endif
            if( cline%defined('vol_dim') )then
                vol_dim = cline%get_iarg('vol_dim')
                call molecule%pdb2mrc(pdbfile=pdb_source, volfile=potential_file, smpd=smpd, &
                    &center_pdb=.true., vol_dim=vol_dim)
            else
                call molecule%pdb2mrc(pdbfile=pdb_source, volfile=potential_file, smpd=smpd, center_pdb=.true.)
            endif
            call molecule%kill
            call cline%set('vol1', potential_file)
        endif
        nptcls_arg = cline%get_iarg('nptcls')
        call cline%set('nspace', max(2, 2 * ((nptcls_arg + 1) / 2)))
        call cline%set('eo', 'no')
        call build%init_params_and_build_general_tbox(cline, params)
        select case(trim(params%wave_model))
            case('projection')
                model = WAVE_PROJECTION
            case('first_born')
                model = WAVE_FIRST_BORN
            case('rytov')
                model = WAVE_RYTOV
            case('multislice')
                model = WAVE_MSLICE
            case DEFAULT
                THROW_HARD('unknown wave_model; simulate_scattering')
        end select
        select case(trim(params%ewald))
            case('neg')
                kappa = -1
            case('pos')
                kappa = 1
            case DEFAULT
                kappa = 0
        end select
        if( model == WAVE_PROJECTION .and. kappa /= 0 ) THROW_HARD('wave_model=projection samples the central section: ewald=no')
        if( model == WAVE_MSLICE .and. kappa == 0 ) THROW_HARD('wave_model=multislice propagates through the box: ewald=neg or pos')
        if( params%total_dose < 0.0 .or. params%noise_sdev < 0.0 ) &
            &THROW_HARD('total_dose and noise_sdev must be nonnegative')
        call find_ldim_nptcls(params%vols(1), ldim, nsec)
        if( any(ldim /= ldim(1)) ) THROW_HARD('simulate_scattering requires a cubic potential volume')
        call build%vol%read(params%vols(1))
        if( .not. cline%defined('oritab') ) call build%spproj_field%rnd_oris(0.0, params%eullims)
        if( .not. build%spproj_field%isthere('dfx') ) &
            &call build%spproj_field%rnd_ctf(params%kv, params%cs, params%fraca, &
            &params%defocus, params%dferr, params%astigerr)
        has_dfy     = build%spproj_field%isthere('dfy')
        has_angast  = build%spproj_field%isthere('angast')
        has_phshift = build%spproj_field%isthere('phshift')
        ! the particles: orientation, shift and the CTF at the potential's sampling
        call selection%new(params%nptcls, .true.)
        do i = 1, params%nptcls
            call build%spproj_field%get_ori(i, orientation)
            ctfparms%ctfflag = CTFFLAG_YES
            ctfparms%smpd    = params%smpd
            ctfparms%kv      = params%kv
            ctfparms%cs      = params%cs
            ctfparms%fraca   = params%fraca
            ctfparms%dfx     = orientation%get('dfx')
            ctfparms%dfy     = ctfparms%dfx
            ctfparms%angast  = 0.0
            ctfparms%phshift = 0.0
            if( has_dfy     ) ctfparms%dfy     = orientation%get('dfy')
            if( has_angast  ) ctfparms%angast  = orientation%get('angast')
            if( has_phshift ) ctfparms%phshift = orientation%get('phshift')
            call orientation%set_ctfvars(ctfparms)
            call selection%set_ori(i, orientation)
        enddo
        write(logfhandle,'(A)') '>>> GENERATING WAVE-SCATTERING IMAGES: '//trim(params%wave_model)
        ! the shared operator: relative intensity on the native disk
        allocate(wproj)
        call wproj%new(params%box, params%smpd, model, kappa, 1.0, bfac=params%bfac, fft_nthreads=params%nthr)
        call wproj%prep_particles(selection)
        call wproj%set_volume(build%vol%get_rmat())
        R = params%box / 2
        allocate(planes(-R:R,-R:R,BATCH))
        do ib = 1, params%nptcls, BATCH
            nb = min(BATCH, params%nptcls - ib + 1)
            call wproj%predict(nb, ib, planes(:,:,1:nb))
            do i = 1, nb
                call wproj%plane_to_image(planes(:,:,i), build%img)
                call detect_intensity(build%img, params%total_dose, params%noise_sdev)
                call build%img%write(params%outstk, ib + i - 1)
            enddo
        enddo
        call wproj%kill
        deallocate(wproj, planes)
        ! the truth: orientations, shifts and the CTF parameters the images were made with
        call selection%write(params%outfile, [1,params%nptcls])
        call selection%kill
        call orientation%kill
        call simple_end('**** SIMPLE_SIMULATE_SCATTERING NORMAL STOP ****')
    end subroutine exec_simulate_scattering

    subroutine exec_simulate_movie( self, cline )
        use simple_memoize_ft_maps, only: memoize_ft_maps, forget_ft_maps
        class(commander_simulate_movie), intent(inout) :: self
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
        if( .not. cline%defined('nframes') ) call cline%set('nframes',  30)
        call cline%set('oritype', 'stk')
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
        call base_image%mul(-1.) ! invert contrast
        call memoize_ft_maps(base_image%get_ldim(), base_image%get_smpd())
        do i=1,params%nframes
            call progress(i,params%nframes)
            ! shift base image
            shifted_base_image = base_image
            call shifted_base_image%shift([shifts(i,1),shifts(i,2),0.])
            call shifted_base_image%ifft()
            ! add pink noise
            call shifted_base_image%norm
            call shifted_base_image%add_gauran(snr_pink)
            ! multiply with CTF
            call shifted_base_image%fft()
            if( params%neg .eq. 'yes' )then
                call shifted_base_image%apply_ctf_wpad(tfun, dfx, 'neg', 0., dfy, angast, params%bfac)
            else
                call shifted_base_image%apply_ctf_wpad(tfun, dfx, 'ctf', 0., dfy, angast, params%bfac)
            endif
            call shifted_base_image%ifft()
            ! add the detector noise
            call shifted_base_image%norm
            call shifted_base_image%add_gauran(snr_detector)
            if( params%vis .eq. 'yes' ) call shifted_base_image%vis()
            call shifted_base_image%write(string('simulate_movie')//params%ext, i)
            ! set orientation parameters in object
            call build%spproj_field%set_dfx(1,       dfx)
            call build%spproj_field%set_dfy(1,       dfy)
            call build%spproj_field%set(1, 'angast', angast)
            call build%spproj_field%set(1, 'phshift', 0.)
        end do
        ! generate the optimal average
        base_image = 0.
        write(logfhandle,'(a)') '>>> GENERATING OPTIMAL AVERAGE'
        do i=1,params%nframes
            call progress(i,params%nframes)
            call shifted_base_image%read(string('simulate_movie')//params%ext, i)
            x = build%spproj_field%get(1, 'x'//int2str(i))
            y = build%spproj_field%get(1, 'y'//int2str(i))
            call shifted_base_image%shift([-x,-y,0.])
            call base_image%add(shifted_base_image)
        end do
        call base_image%div(real(params%nframes))
        call base_image%write(string('optimal_movie_average')//params%ext, 1)
        if( params%vis .eq. 'yes' ) call base_image%vis()
        ! output orientations
        call build%spproj_field%write(string('simulate_movie_params')//TXT_EXT, [1,build%spproj_field%get_noris()])
        ! Release image/FFTW state when this commander is invoked in-process by a
        ! workflow.  The standalone executable used to hide these missing kills
        ! by letting the operating system reclaim everything at process exit.
        call base_image%kill
        call shifted_base_image%kill
        call build%kill_general_tbox
        call forget_ft_maps
        if( allocated(ptcl_positions) ) deallocate(ptcl_positions)
        if( allocated(shifts)         ) deallocate(shifts)
        ! end gracefully
        call simple_end('**** SIMPLE_SIMULATE_MOVIE NORMAL STOP ****')

        contains

            !> \brief  generate mutually exclusive positions
            function gen_ptcl_pos( npos, xdim, ydim, box ) result( pos )
                integer, intent(in)           :: npos, xdim, ydim
                integer, intent(in), optional :: box
                integer, allocatable          :: pos(:,:)
                logical                       :: occupied(xdim,ydim), skip
                integer                       :: ix, iy, cnt, ntries, max_tries
                integer                       :: x_min, x_max, y_min, y_max, b_half
                allocate( pos(npos,2) )
                occupied = .false.
                cnt = 0
                ntries = 0
                max_tries = max(1000, 500*npos)
                if( present(box) )then
                    if( box < 1 .or. box >= xdim .or. box >= ydim )then
                        THROW_HARD('particle box must fit inside the simulated micrograph; gen_ptcl_pos')
                    endif
                endif
                do while( cnt < npos )
                    ntries = ntries + 1
                    if( ntries >= max_tries )then
                        THROW_HARD('could not place all particles; reduce nptcls or increase xdim/ydim; gen_ptcl_pos')
                    endif
                    skip = .false.
                    ix = irnd_uni(xdim)
                    iy = irnd_uni(ydim)
                    if( present(box) )then
                        b_half = box / 2
                        x_min  = ix - b_half
                        !x_max  = ix + box - b_half - 1
                        x_max  = ix + b_half - 1
                        y_min  = iy - b_half
                        y_max  = iy + b_half - 1
                        if( x_min < 1 .or. x_max > xdim ) cycle
                        if( y_min < 1 .or. y_max > ydim ) cycle
                        skip = any(occupied(x_min:x_max,y_min:y_max))
                    else
                        skip = occupied(ix,iy)
                    endif
                    if( skip ) cycle
                    if( present(box) )then
                        occupied(x_min:x_max,y_min:y_max) = .true.
                    else
                        occupied(ix,iy) = .true.
                    endif
                    cnt = cnt+1
                    call progress(cnt,npos)
                    pos(cnt,1) = ix
                    pos(cnt,2) = iy
                end do
            end function gen_ptcl_pos

    end subroutine exec_simulate_movie

    subroutine exec_simulate_subtomogram( self, cline )
        use simple_simple_volinterp, only: rotvol
        class(commander_simulate_subtomogram), intent(inout) :: self
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
        call build%vol%norm
        call build%vol%add_gauran(params%snr)
        call o%new(is_ptcl=.false.)
        numlen = len(int2str(params%nptcls))
        do iptcl=1,params%nptcls
            call o%rnd_ori
            vol_rot = rotvol(build%vol, o)
            call vol_rot%write(string('subtomo')//int2str_pad(iptcl,numlen)//params%ext%to_char())
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_SIMULATE_SUBTOMOGRAM NORMAL STOP ****')
    end subroutine exec_simulate_subtomogram

    subroutine exec_simulate_nanoparticle( self, cline )
        use simple_atoms, only: atoms
        class(commander_simulate_nanoparticle), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        real, parameter   :: NPIX_CUTOFF = 12.
        type(parameters)  :: params
        type(image)       :: vol
        type(atoms)       :: atoms_obj
        character(len=2)  :: el1, el2
        character(len=5)  :: el_ucase
        character(len=10) :: crystal_system
        real              :: center(3), ha, hc, x, x1, x2, x3, y, y1, y2, y3, z, z1, z2, z3, msksq, cutoff
        real              :: a(3) ! lattice parameters
        real              :: pos(3), u_wz, bw1(3), bw2(3), bw3(3), bw4(3)
        integer           :: i, j, k, n, ncubes, ncubes_c, ncubes_a
        if( .not. cline%defined('pdbout') ) call cline%set('pdbout', 'simatms.pdb')
        call params%new(cline)
        if( .not.is_even(params%box) ) THROW_HARD('BOX must be even')
        call vol%new([params%box,params%box,params%box], params%smpd)
        if( cline%defined('pdbfile') )then
            if( .not.file_exists(params%pdbfile) ) THROW_HARD('Could not find: '//params%pdbfile%to_char())
            call atoms_obj%new(params%pdbfile)
            cutoff = NPIX_CUTOFF * params%smpd
        else if( cline%defined('element') )then
            if( .not. cline%defined('moldiam') ) THROW_HARD('MOLDIAM must be provided for lattice!')
            if( .not. atoms_obj%element_exists(params%element) ) THROW_HARD('Unsupported element for now')
            el1      = params%element(1:2)
            el2      = params%element(3:4)
            msksq    = (params%moldiam / 2.)**2.
            el_ucase = uppercase(trim(adjustl(params%element)))
            call get_lattice_params(el_ucase, crystal_system, a)
            write(logfhandle,'(a)') '>>> ELEMENT: '//el_ucase
            write(logfhandle,'(a)') '>>> CRYSTAL SYSTEM: '//crystal_system
            ha       = a(1) / 2.; hc = a(3) / 2.
            ncubes   = floor(real(params%box) * params%smpd / a(1))
            ncubes_c = floor(real(params%box) * params%smpd / a(3))
            ! atoms at edges
            n = 0
            center = (real([params%box,params%box,params%box]/2 + 1) - 1.) * params%smpd
            select case( trim(crystal_system) )
                case('fcc')
                    ! count
                    do i = 1, ncubes
                        x  = real(i - 1) * a(1)
                        x1 = x + ha; x2 = x; x3 = x + ha
                        do j = 1, ncubes
                            y  = real(j - 1) * a(1)
                            y1 = y + ha; y2 = y + ha; y3 = y
                            do k = 1, ncubes
                                z  = real(k - 1) * a(1)
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
                        call atoms_obj%set_name(i,params%element(1:2)//'  ')
                        call atoms_obj%set_element(i,el1)
                    enddo
                    ! generate all coordinates
                    n = 0
                    do i = 1, ncubes
                        x  = real(i - 1) * a(1); x1 = x + ha; x2 = x; x3 = x + ha
                        do j = 1, ncubes
                            y  = real(j - 1) * a(1); y1 = y + ha; y2 = y + ha; y3 = y
                            do k = 1, ncubes
                                z  = real(k - 1) * a(1); z1 = z; z2 = z + ha; z3 = z + ha
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
                        x  = real(i - 1) * a(1); x1 = x + ha
                        do j = 1, ncubes
                            y  = real(j - 1) * a(1); y1 = y + ha
                            do k= 1, ncubes
                                z  = real(k - 1) * a(1) ; z1 = z + ha
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
                        call atoms_obj%set_name(i,params%element(1:2)//'  ')
                        call atoms_obj%set_element(i,el1)
                    enddo
                    ! generate all coordinates
                    n = 0
                    do i = 1, ncubes
                        x  = real(i - 1) * a(1); x1 = x + ha
                        do j = 1, ncubes
                            y  = real(j - 1) * a(1); y1 = y + ha
                            do k = 1, ncubes
                                z  = real(k - 1) * a(1); z1 = z + ha
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
                    ! Rocksalt structure (cubic, space group Fm-3m, #225):
                    ! element 1 occupies the FCC sites and element 2 occupies
                    ! the octahedral sites displaced by a/2 along a cube axis.
                    ! count
                    do i = 1, ncubes
                        x  = real(i - 1) * a(1); x1 = x + ha
                        do j = 1, ncubes
                            y  = real(j - 1) * a(1); y1 = y + ha
                            do k = 1, ncubes
                                z  = real(k - 1) * a(1); z1 = z + ha
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
                        x  = real(i - 1) * a(1); x1 = x + ha
                        do j = 1, ncubes
                            y  = real(j - 1) * a(1); y1 = y + ha
                            do k = 1, ncubes
                                z  = real(k - 1) * a(1); z1 = z + ha
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
                case('wurtzite')
                    ! Wurtzite structure (hexagonal, space group P6_3mc, #186)
                    ! Lattice vectors (Cartesian):
                    !   a1 = ( a,       0,           0 )
                    !   a2 = (-a/2,     a*sqrt(3)/2, 0 )
                    !   a3 = ( 0,       0,           c )
                    ! Basis atoms in fractional coords (2 of each element per cell):
                    !   el1: (1/3, 2/3, 0  )  and  (2/3, 1/3, 1/2  )
                    !   el2: (1/3, 2/3, u  )  and  (2/3, 1/3, 1/2+u)
                    !   with u = 3/8 (ideal wurtzite internal parameter)
                    ! Cartesian offsets from cell origin:
                    !   bw1 = (0,   a*sqrt(3)/3,  0          )  [el1]
                    !   bw2 = (a/2, a*sqrt(3)/6,  c/2        )  [el1]
                    !   bw3 = (0,   a*sqrt(3)/3,  u*c        )  [el2]
                    !   bw4 = (a/2, a*sqrt(3)/6,  (1/2+u)*c  )  [el2]
                    u_wz = 3. / 8.
                    bw1  = [0.,        a(1)*sqrt(3.)/3.,  0.                    ]
                    bw2  = [a(1)/2.,   a(1)*sqrt(3.)/6.,  a(3)/2.               ]
                    bw3  = [0.,        a(1)*sqrt(3.)/3.,  u_wz*a(3)             ]
                    bw4  = [a(1)/2.,   a(1)*sqrt(3.)/6.,  (0.5 + u_wz)*a(3)     ]
                    ! Loop bounds for hexagonal grid
                    ! j covers y in steps of a*sqrt(3)/2
                    ! i covers x in steps of a, but shifted by -a/2 per j-row (shear)
                    ! k covers z in steps of c
                    ncubes   = ceiling(real(params%box)*params%smpd / (a(1)*sqrt(3.)/2.)) + 1
                    ncubes_a = ceiling(real(params%box)*params%smpd / a(1)) + ncubes/2 + 2
                    ncubes_c = ceiling(real(params%box)*params%smpd / a(3)) + 1
                    ! count atoms inside spherical mask
                    do i = 1, ncubes_a
                        do j = 1, ncubes
                            x = real(i-1)*a(1) + real(j-1)*(-a(1)/2.)
                            y = real(j-1)*(a(1)*sqrt(3.)/2.)
                            do k = 1, ncubes_c
                                z = real(k-1)*a(3)
                                pos = [x,y,z] + bw1; if( sum((pos-center)**2.) <= msksq ) n = n+1
                                pos = [x,y,z] + bw2; if( sum((pos-center)**2.) <= msksq ) n = n+1
                                pos = [x,y,z] + bw3; if( sum((pos-center)**2.) <= msksq ) n = n+1
                                pos = [x,y,z] + bw4; if( sum((pos-center)**2.) <= msksq ) n = n+1
                            enddo
                        enddo
                    enddo
                    ! generate atoms
                    call atoms_obj%new(n)
                    n = 0
                    do i = 1, ncubes_a
                        do j = 1, ncubes
                            x = real(i-1)*a(1) + real(j-1)*(-a(1)/2.)
                            y = real(j-1)*(a(1)*sqrt(3.)/2.)
                            do k = 1, ncubes_c
                                z = real(k-1)*a(3)
                                ! el1 atom 1
                                pos = [x,y,z] + bw1
                                if( sum((pos-center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, pos)
                                    call atoms_obj%set_name(n, el1//'  ')
                                    call atoms_obj%set_element(n, el1)
                                endif
                                ! el1 atom 2
                                pos = [x,y,z] + bw2
                                if( sum((pos-center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, pos)
                                    call atoms_obj%set_name(n, el1//'  ')
                                    call atoms_obj%set_element(n, el1)
                                endif
                                ! el2 atom 1
                                pos = [x,y,z] + bw3
                                if( sum((pos-center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, pos)
                                    call atoms_obj%set_name(n, el2//'  ')
                                    call atoms_obj%set_element(n, el2)
                                endif
                                ! el2 atom 2
                                pos = [x,y,z] + bw4
                                if( sum((pos-center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, pos)
                                    call atoms_obj%set_name(n, el2//'  ')
                                    call atoms_obj%set_element(n, el2)
                                endif
                            enddo
                        enddo
                    enddo
                case('zincblende')
                    ! Zincblende structure (cubic, space group F-43m, #216)
                    ! FCC lattice with two-atom basis:
                    !   el1 at each FCC site: (0,0,0), (a/2,a/2,0), (a/2,0,a/2), (0,a/2,a/2)
                    !   el2 offset by (a/4,a/4,a/4) from each FCC site
                    ! count atoms inside spherical mask
                    do i = 1, ncubes
                        x = real(i-1)*a(1)
                        do j = 1, ncubes
                            y = real(j-1)*a(1)
                            do k = 1, ncubes
                                z = real(k-1)*a(1)
                                ! 4 FCC positions for el1
                                if( sum(([x,      y,      z      ] - center)**2.) <= msksq ) n = n+1
                                if( sum(([x+ha,   y+ha,   z      ] - center)**2.) <= msksq ) n = n+1
                                if( sum(([x+ha,   y,      z+ha   ] - center)**2.) <= msksq ) n = n+1
                                if( sum(([x,      y+ha,   z+ha   ] - center)**2.) <= msksq ) n = n+1
                                ! 4 FCC positions for el2, offset by (a/4,a/4,a/4)
                                if( sum(([x+a(1)/4.,   y+a(1)/4.,   z+a(1)/4.   ] - center)**2.) <= msksq ) n = n+1
                                if( sum(([x+3.*a(1)/4.,y+3.*a(1)/4.,z+a(1)/4.   ] - center)**2.) <= msksq ) n = n+1
                                if( sum(([x+3.*a(1)/4.,y+a(1)/4.,   z+3.*a(1)/4.] - center)**2.) <= msksq ) n = n+1
                                if( sum(([x+a(1)/4.,   y+3.*a(1)/4.,z+3.*a(1)/4.] - center)**2.) <= msksq ) n = n+1
                            enddo
                        enddo
                    enddo
                    ! generate atoms
                    call atoms_obj%new(n)
                    n = 0
                    do i = 1, ncubes
                        x = real(i-1)*a(1)
                        do j = 1, ncubes
                            y = real(j-1)*a(1)
                            do k = 1, ncubes
                                z = real(k-1)*a(1)
                                ! el1 at FCC corner
                                pos = [x, y, z]
                                if( sum((pos - center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, pos)
                                    call atoms_obj%set_name(n, el1//'  ')
                                    call atoms_obj%set_element(n, el1)
                                endif
                                ! el1 at face (xy)
                                pos = [x+ha, y+ha, z]
                                if( sum((pos - center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, pos)
                                    call atoms_obj%set_name(n, el1//'  ')
                                    call atoms_obj%set_element(n, el1)
                                endif
                                ! el1 at face (xz)
                                pos = [x+ha, y, z+ha]
                                if( sum((pos - center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, pos)
                                    call atoms_obj%set_name(n, el1//'  ')
                                    call atoms_obj%set_element(n, el1)
                                endif
                                ! el1 at face (yz)
                                pos = [x, y+ha, z+ha]
                                if( sum((pos - center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, pos)
                                    call atoms_obj%set_name(n, el1//'  ')
                                    call atoms_obj%set_element(n, el1)
                                endif
                                ! el2 at (a/4,a/4,a/4) from corner
                                pos = [x+a(1)/4., y+a(1)/4., z+a(1)/4.]
                                if( sum((pos - center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, pos)
                                    call atoms_obj%set_name(n, el2//'  ')
                                    call atoms_obj%set_element(n, el2)
                                endif
                                ! el2 at (a/4,a/4,a/4) from face (xy)
                                pos = [x+3.*a(1)/4., y+3.*a(1)/4., z+a(1)/4.]
                                if( sum((pos - center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, pos)
                                    call atoms_obj%set_name(n, el2//'  ')
                                    call atoms_obj%set_element(n, el2)
                                endif
                                ! el2 at (a/4,a/4,a/4) from face (xz)
                                pos = [x+3.*a(1)/4., y+a(1)/4., z+3.*a(1)/4.]
                                if( sum((pos - center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, pos)
                                    call atoms_obj%set_name(n, el2//'  ')
                                    call atoms_obj%set_element(n, el2)
                                endif
                                ! el2 at (a/4,a/4,a/4) from face (yz)
                                pos = [x+a(1)/4., y+3.*a(1)/4., z+3.*a(1)/4.]
                                if( sum((pos - center)**2.) <= msksq )then
                                    n = n + 1
                                    call atoms_obj%set_coord(n, pos)
                                    call atoms_obj%set_name(n, el2//'  ')
                                    call atoms_obj%set_element(n, el2)
                                endif
                            enddo
                        enddo
                    enddo
                case default
                    THROW_HARD('Unsupported crystal system: '//trim(crystal_system))
            end select
            ! write PDB
            do i=1,n
                call atoms_obj%set_num(i,i)
                call atoms_obj%set_resnum(i,i)
                call atoms_obj%set_chain(i,'A')
                call atoms_obj%set_occupancy(i,1.)
            enddo
            ! for convolution
            cutoff = 3.*a(3) ! takes the biggest lattice parameter as cutoff for convolution
        else
            THROW_HARD('A PDB file or ELEMENT & MOLDIAM must be defined!')
        endif
        call atoms_obj%writePDB(params%pdbout)
        if( cline%defined('lp') )then
            cutoff = max(cutoff,4.*params%lp)
            call atoms_obj%convolve(vol, cutoff, lp=params%lp, bfac_pdb=params%l_pdb_bfac)
        else
            call atoms_obj%convolve(vol, cutoff, bfac_pdb=params%l_pdb_bfac)
        endif
        call vol%write(params%outvol)
        ! cleanup
        call atoms_obj%kill
        call vol%kill
        ! end gracefully
        call simple_end('**** SIMPLE_SIMULATE_NANOPARTICLE NORMAL STOP ****')
    end subroutine exec_simulate_nanoparticle

    !> whether the orientation document fname gives particles a state (a state= key on a line)
    logical function oritab_holds_states( fname )
        class(string), intent(in) :: fname
        character(len=4096) :: line
        integer :: funit, ios
        oritab_holds_states = .false.
        open(newunit=funit, file=fname%to_char(), status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        do
            read(funit,'(a)',iostat=ios) line
            if( ios /= 0 ) exit
            if( index(' '//line, ' state=') > 0 .or. index(line, achar(9)//'state=') > 0 )then
                oritab_holds_states = .true.
                exit
            endif
        enddo
        close(funit)
    end function oritab_holds_states

    !> The population fractions of nstates states: the comma-separated list pops, or equal when it is
    !! blank; they are non-negative and sum to 1 (within 1e-3)
    function state_fractions( pops, nstates ) result( fracs )
        use simple_string_utils, only: parsestr
        character(len=*), intent(in) :: pops
        integer,          intent(in) :: nstates
        real, allocatable     :: fracs(:)
        character(len=STDLEN) :: list
        character(len=32)     :: items(MAXS)
        integer :: n, i, ios
        allocate(fracs(nstates), source=1. / real(nstates))
        if( len_trim(pops) == 0 ) return
        list = pops
        call parsestr(list, ',', items, n)
        if( n /= nstates ) THROW_HARD('state_pops must hold one fraction per state')
        do i = 1, n
            read(items(i),*,iostat=ios) fracs(i)
            if( ios /= 0 ) THROW_HARD('state_pops is not a list of numbers: '//trim(pops))
        enddo
        if( any(fracs < 0.) ) THROW_HARD('state_pops holds a negative fraction')
        if( abs(sum(fracs) - 1.) > 1.e-3 ) THROW_HARD('the fractions of state_pops must sum to 1')
    end function state_fractions

    !> The states of n particles in exact counts: each fraction times n, rounded, the remainder given
    !! to (or taken from) the states in order of decreasing fraction; the order of the labels is
    !! shuffled with the seeded generator
    function draw_states( n, fracs ) result( states )
        use simple_rnd, only: shuffle
        integer, intent(in) :: n
        real,    intent(in) :: fracs(:)
        integer, allocatable :: states(:)
        integer :: counts(size(fracs)), order(size(fracs)), k, i, j, nleft
        logical :: taken(size(fracs))
        counts = nint(fracs * real(n))
        ! the states by decreasing fraction, ties in index order
        taken = .false.
        do i = 1, size(fracs)
            j = maxloc(fracs, dim=1, mask=.not. taken)
            order(i) = j
            taken(j) = .true.
        enddo
        nleft = n - sum(counts)
        i = 0
        do while( nleft /= 0 )
            k = order(mod(i, size(fracs)) + 1)
            if( nleft > 0 )then
                counts(k) = counts(k) + 1
                nleft     = nleft - 1
            else if( counts(k) > 0 )then
                counts(k) = counts(k) - 1
                nleft     = nleft + 1
            endif
            i = i + 1
        enddo
        allocate(states(n))
        j = 0
        do k = 1, size(fracs)
            states(j+1:j+counts(k)) = k
            j = j + counts(k)
        enddo
        call shuffle(states)
    end function draw_states

    !> simulate_particles in wave mode: the wave operator's image of each state's map, converted to volts with
    !! the scale link, detected with shot and detector noise (simimg_wave); the particles of the part in
    !! output order, cnt of them, recthres applied as in the flat mode
    subroutine simulate_particles_wave( params, build, cline, o_thres, o_mirr, cnt )
        use simple_wave_imaging,   only: WAVE_PROJECTION, WAVE_FIRST_BORN, WAVE_RYTOV, WAVE_MSLICE
        use simple_wave_projector, only: wave_projector
        use simple_math_ctf,       only: electron_scale_link
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        class(ori),        intent(in)    :: o_thres, o_mirr
        integer,           intent(out)   :: cnt
        integer, parameter :: BATCH = 64
        type(wave_projector), allocatable :: wprojs(:)
        type(oris),           allocatable :: sels(:)
        complex,              allocatable :: planes(:,:,:), buf(:,:,:)
        real,                 allocatable :: bfacs(:), vol_volts(:,:,:)
        integer,              allocatable :: inds(:), states(:), isel(:), nsel(:)
        type(ori)       :: orientation
        type(ctfparams) :: ctfparms
        type(string)    :: fname
        real    :: bfacerr, scale_link
        integer :: model, kappa, nstates, n, i, j, j0, j1, nb, istate, a, b, r, ldim_vol(3), nsec
        logical :: has_dfy, has_angast, has_phshift
        select case(trim(params%wave_model))
            case('projection')
                model = WAVE_PROJECTION
            case('first_born')
                model = WAVE_FIRST_BORN
            case('rytov')
                model = WAVE_RYTOV
            case('multislice')
                model = WAVE_MSLICE
            case DEFAULT
                THROW_HARD('unknown wave_model; simulate_particles')
        end select
        select case(trim(params%ewald))
            case('neg')
                kappa = -1
            case('pos')
                kappa = 1
            case DEFAULT
                kappa = 0
        end select
        if( model == WAVE_PROJECTION .and. kappa /= 0 ) THROW_HARD('wave_model=projection samples the central section: ewald=no')
        if( model == WAVE_MSLICE .and. kappa == 0 ) THROW_HARD('wave_model=multislice propagates through the box: ewald=neg or pos')
        if( params%total_dose <= 0. ) THROW_HARD('total_dose must be positive; simulate_particles')
        if( params%snr <= 0. ) THROW_HARD('snr must be positive; simulate_particles')
        nstates = max(1, params%nstates)
        ! the particles of the part in output order, their states and lens B-factors
        allocate(inds(params%top - params%fromp + 1), states(params%top - params%fromp + 1))
        n = 0
        do i = params%fromp, params%top
            call build%spproj_field%get_ori(i, orientation)
            if( trim(params%recthres).eq.'yes' )then
                if( rad2deg(orientation .euldist. o_thres) < params%rec_athres .or.&
                   &rad2deg(orientation .euldist. o_mirr ) < params%rec_athres ) cycle
            endif
            n         = n + 1
            inds(n)   = i
            states(n) = orientation%get_state()
        end do
        allocate(bfacs(n), source=0.)
        if( cline%defined('bfac') )then
            do j = 1, n
                bfacerr = ran3() * params%bfacerr
                if( ran3() < 0.5 )then
                    bfacs(j) = params%bfac - bfacerr
                else
                    bfacs(j) = params%bfac + bfacerr
                endif
            end do
        endif
        ! each state's selection, in output order: particle j is particle isel(j) of its state's
        has_dfy     = build%spproj_field%isthere('dfy')
        has_angast  = build%spproj_field%isthere('angast')
        has_phshift = build%spproj_field%isthere('phshift')
        allocate(sels(nstates), nsel(nstates), isel(n))
        nsel = 0
        do j = 1, n
            nsel(states(j)) = nsel(states(j)) + 1
            isel(j)         = nsel(states(j))
        end do
        do istate = 1, nstates
            call sels(istate)%new(nsel(istate), is_ptcl=.true.)
        end do
        do j = 1, n
            call build%spproj_field%get_ori(inds(j), orientation)
            ctfparms%ctfflag = CTFFLAG_YES
            ctfparms%smpd    = params%smpd
            ctfparms%kv      = params%kv
            ctfparms%cs      = params%cs
            ctfparms%fraca   = params%fraca
            ctfparms%dfx     = orientation%get('dfx')
            ctfparms%dfy     = ctfparms%dfx
            ctfparms%angast  = 0.
            ctfparms%phshift = 0.
            if( has_dfy     ) ctfparms%dfy     = orientation%get('dfy')
            if( has_angast  ) ctfparms%angast  = orientation%get('angast')
            if( has_phshift ) ctfparms%phshift = orientation%get('phshift')
            call orientation%set_ctfvars(ctfparms)
            call sels(states(j))%set_ori(isel(j), orientation)
        end do
        ! one operator per state, its map (masked as in the flat mode) in volts
        scale_link = electron_scale_link(params%total_dose, params%kv, params%smpd, params%box)
        allocate(wprojs(nstates))
        do istate = 1, nstates
            if( nsel(istate) == 0 ) cycle
            call find_ldim_nptcls(params%vols(istate), ldim_vol, nsec)
            if( any(ldim_vol /= [params%box, params%box, params%box]) ) THROW_HARD('the volumes of the states differ in box')
            call build%vol%read(params%vols(istate))
            call build%vol%mask3D_soft(params%msk, backgr=0.)
            vol_volts = build%vol%get_rmat() / scale_link
            call wprojs(istate)%new(params%box, params%smpd, model, kappa, 1.0, fft_nthreads=params%nthr)
            call wprojs(istate)%prep_particles(sels(istate), bfacs=pack(bfacs, states(1:n) == istate))
            call wprojs(istate)%set_volume(vol_volts)
        end do
        write(logfhandle,'(A)') '>>> GENERATING WAVE IMAGES: '//trim(params%wave_model)
        if( cline%defined('part') )then
            fname = string('simulated_particles_part')//int2str_pad(params%part,params%numlen)//params%ext%to_char()
        else
            fname = params%outstk
        endif
        r = params%box / 2
        allocate(planes(-r:r,-r:r,BATCH), buf(-r:r,-r:r,BATCH))
        do j0 = 1, n, BATCH
            j1 = min(n, j0 + BATCH - 1)
            ! a state's particles of the batch are consecutive in its selection
            do istate = 1, nstates
                if( .not. any(states(j0:j1) == istate) ) cycle
                a  = minval(isel(j0:j1), mask=states(j0:j1) == istate)
                b  = maxval(isel(j0:j1), mask=states(j0:j1) == istate)
                nb = b - a + 1
                call wprojs(istate)%predict(nb, a, buf(:,:,1:nb))
                do j = j0, j1
                    if( states(j) == istate ) planes(:,:,j-j0+1) = buf(:,:,isel(j)-a+1)
                end do
            end do
            do j = j0, j1
                call progress(j, n)
                call wprojs(states(j))%plane_to_image(planes(:,:,j-j0+1), build%img)
                call simimg_wave(build%img, params%total_dose, params%snr)
                call build%img%write(fname, j)
            end do
        end do
        cnt = n
        do istate = 1, nstates
            call wprojs(istate)%kill
            call sels(istate)%kill
        end do
        deallocate(wprojs, sels, planes, buf)
        call orientation%kill
    end subroutine simulate_particles_wave

end module simple_commanders_sim
