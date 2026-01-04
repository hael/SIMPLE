module single_commanders_trajectory
use simple_commander_module_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_track_particles_distr
  contains
    procedure :: execute      => exec_track_particles_distr
end type commander_track_particles_distr

type, extends(commander_base) :: commander_track_particles
  contains
    procedure :: execute      => exec_track_particles
end type commander_track_particles

type, extends(commander_base) :: commander_trajectory_backgr_subtr
  contains
    procedure :: execute      => exec_trajectory_backgr_subtr
end type commander_trajectory_backgr_subtr

type, extends(commander_base) :: commander_graphene_subtr
  contains
    procedure :: execute      => exec_graphene_subtr
end type commander_graphene_subtr

type, extends(commander_base) :: commander_import_trajectory
  contains
    procedure :: execute      => exec_import_trajectory
end type commander_import_trajectory

type, extends(commander_base) :: commander_denoise_trajectory
  contains
    procedure :: execute      => exec_denoise_trajectory
end type commander_denoise_trajectory

type, extends(commander_base) :: commander_trajectory_swap_stack
  contains
    procedure :: execute      => exec_trajectory_swap_stack
end type commander_trajectory_swap_stack

type, extends(commander_base) :: commander_extract_substk
  contains
    procedure :: execute      => exec_extract_substk
end type commander_extract_substk

contains

    subroutine exec_track_particles_distr( self, cline )
        class(commander_track_particles_distr), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters)              :: params
        type(qsys_env)                :: qenv
        type(chash)                   :: job_descr
        type(nrtxtfile)               :: boxfile
        real,        allocatable      :: boxdata(:,:)
        type(chash), allocatable      :: part_params(:)
        integer :: ndatlines, numlen, j, orig_box, ipart
        if( .not. cline%defined('neg')       ) call cline%set('neg',      'yes')
        if( .not. cline%defined('lp')        ) call cline%set('lp',         2.3)
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      5.0)
        if( .not. cline%defined('nframesgrp')) call cline%set('nframesgrp',  30)
        if( .not. cline%defined('filter'))     call cline%set('filter',    'tv')
        if( .not. cline%defined('offset'))     call cline%set('offset',     10.)
        call cline%set('oritype','mic')
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        if( .not. file_exists(params%boxfile) ) THROW_HARD('inputted boxfile does not exist in cwd')
        if( nlines(params%boxfile) > 0 )then
            call boxfile%new(params%boxfile, 1)
            ndatlines = boxfile%get_ndatalines()
            numlen    = len(int2str(ndatlines))
            allocate( boxdata(ndatlines,boxfile%get_nrecs_per_line()) )
            do j=1,ndatlines
                call boxfile%readNextDataLine(boxdata(j,:))
                orig_box = nint(boxdata(j,3))
                if( nint(boxdata(j,3)) /= nint(boxdata(j,4)) )then
                    THROW_HARD('Only square windows allowed!')
                endif
            end do
        else
            THROW_HARD('inputted boxfile is empty; exec_track_particles')
        endif
        call boxfile%kill
        params%nptcls  = ndatlines
        params%nparts  = params%nptcls
        params%ncunits = params%nparts
        ! box and numlen need to be part of command line
        if( .not. cline%defined('hp') ) call cline%set('hp', real(orig_box) )
        call cline%set('box',    real(orig_box))
        call cline%set('numlen', real(numlen)  )
        call cline%delete('fbody')
        ! prepare part-dependent parameters
        allocate(part_params(params%nparts))
        do ipart=1,params%nparts
            call part_params(ipart)%new(4)
            call part_params(ipart)%set('xcoord', real2str(boxdata(ipart,1)))
            call part_params(ipart)%set('ycoord', real2str(boxdata(ipart,2)))
            call part_params(ipart)%set('box',    real2str(boxdata(ipart,3)))
            if( params%nparts > 1 )then
                call part_params(ipart)%set('fbody', params%fbody%to_char()//'_'//int2str_pad(ipart,numlen))
            else
                call part_params(ipart)%set('fbody', params%fbody%to_char())
            endif
        end do
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! schedule & clean
        call cline%gen_job_descr(job_descr)
        call qenv%gen_scripts_and_schedule_jobs(job_descr, part_params=part_params, array=L_USE_SLURM_ARR)
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_track_trajectory NORMAL STOP ****')
    end subroutine exec_track_particles_distr

    subroutine exec_track_particles( self, cline )
        use simple_track_trajectory
        class(commander_track_particles), intent(inout) :: self
        class(cmdline),                           intent(inout) :: cline
        type(sp_project)          :: spproj
        type(parameters)          :: params
        type(string)              :: dir, forctf
        type(string), allocatable :: intg_names(:), frame_names(:)
        real,         allocatable :: boxdata(:,:)
        logical,      allocatable :: frames_are_there(:), intgs_are_there(:)
        integer :: i, orig_box, nframes
        call cline%set('oritype','mic')
        call params%new(cline)
        orig_box = params%box
        ! coordinates input
        if( cline%defined('xcoord') .and. cline%defined('ycoord') )then
            if( .not. cline%defined('box') ) THROW_HARD('need box to be part of command line for this mode of execution; exec_track_particles')
            allocate( boxdata(1,2) )
            boxdata(1,1) = real(params%xcoord)
            boxdata(1,2) = real(params%ycoord)
        else
            THROW_HARD('need xcoord/ycoord to be part of command line; exec_track_particles')
        endif
        ! frames input
        call spproj%read(params%projfile)
        nframes = spproj%get_nframes()
        allocate(frames_are_there(nframes), intgs_are_there(nframes), source=.false.)
        do i = 1,nframes
            frames_are_there(i) = spproj%os_mic%isthere(i,'frame')
            intgs_are_there(i)  = spproj%os_mic%isthere(i,'intg')
        enddo
        if( all(frames_are_there) )then
            allocate(frame_names(nframes))
            do i = 1,nframes
                frame_names(i) = spproj%os_mic%get_str(i,'frame')
            enddo
        endif
        if( all(intgs_are_there) )then
            allocate(intg_names(nframes))
            do i = 1,nframes
                intg_names(i) = spproj%os_mic%get_str(i,'intg')
            enddo
        endif
        ! actual tracking
        dir = params%fbody
        call simple_mkdir(dir)
        if( allocated(frame_names) )then
            if( allocated(intg_names) )then
                call init_tracker( nint(boxdata(1,1:2)), intg_names, frame_names, dir, params%fbody)
            else
                call init_tracker( nint(boxdata(1,1:2)), frame_names, frame_names, dir, params%fbody)
            endif
        endif
        if( cline%defined('fromf') )then
            call track_particle( forctf, params%fromf )
        else
            call track_particle( forctf )
        endif
        ! clean tracker
        call kill_tracker
        ! end gracefully
        call qsys_job_finished(string('single_commanders_trajectory :: exec_track_particles'))
        call spproj%kill
        call simple_end('**** SIMPLE_track_trajectory NORMAL STOP ****')
    end subroutine exec_track_particles

    subroutine exec_trajectory_backgr_subtr( self, cline )
        ! for background subtraction in time-series data. The goal is to subtract the two graphene
        ! peaks @ 2.14 A and @ 1.23 A. This is done by band-pass filtering the background image,
        ! recommended (and default settings) are hp=5.0 lp=1.1 and width=5.0.
        use simple_ctf,   only: ctf
        class(commander_trajectory_backgr_subtr), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(image)      :: img_backgr, img_backgr_wctf, ave_img
        type(ctf)        :: tfun
        type(ctfparams)  :: ctfvars
        type(string)     :: ext,imgname
        real             :: ave,sdev,minv,maxv
        integer          :: iptcl, nptcls, ldim(3)
        logical          :: do_flip, err
        call cline%set('oritype','mic')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! dimensions
        call find_ldim_nptcls(params%stk,ldim,nptcls)
        ! CTF logics
        do_flip = .false.
        if( (build%spproj_field%isthere('ctf')) .and. (cline%defined('ctf').and.params%ctf.eq.'flip') )then
            if( build%spproj%get_nintgs() /= nptcls ) THROW_HARD('Incompatible # of images and micrographs!')
            if( .not.build%spproj_field%isthere('dfx') ) THROW_HARD('Missing CTF parameters')
            do_flip = .true.
        endif
        ! get background image
        call img_backgr%new([params%box,params%box,1], params%smpd)
        call img_backgr_wctf%new([params%box,params%box,1], params%smpd)
        call img_backgr%read(params%stk_backgr, 1)
        ! background image is skipped if sdev is small because this neighbour was
        ! outside the micrograph
        call img_backgr%stats(ave, sdev, maxv, minv, errout=err)
        if( sdev>TINY .and. .not.err )then
            call ave_img%new([params%box,params%box,1], params%smpd)
            ave_img = 0.
            ! main loop
            do iptcl=1,nptcls
                call progress(iptcl,nptcls)
                ! read particle image
                call build%img%read(params%stk, iptcl)
                ! copy background image & CTF
                img_backgr_wctf = img_backgr
                ! fwd ft
                call build%img%fft()
                call img_backgr_wctf%fft()
                if( do_flip )then
                    ctfvars = build%spproj_field%get_ctfvars(iptcl)
                    tfun = ctf(ctfvars%smpd, ctfvars%kv, ctfvars%cs, ctfvars%fraca)
                    call tfun%apply_serial(build%img, 'flip', ctfvars)
                    call tfun%apply_serial(img_backgr_wctf, 'flip', ctfvars)
                endif
                ! filter background image
                call img_backgr_wctf%bp(params%hp,params%lp,width=params%width)
                ! subtract background
                call build%img%subtr(img_backgr_wctf)
                ! bwd ft
                call build%img%ifft()
                ! normalise
                call build%img%norm()
                call ave_img%add(build%img)
                ! output corrected image
                call build%img%write(params%outstk, iptcl)
            end do
            ! generates average
            call ave_img%div(real(nptcls))
            ext     = fname2ext(params%outstk)
            imgname = get_fbody(params%outstk,ext,.true.)//'_ave.'//ext%to_char()
            call ave_img%write(imgname)
        endif
        call ave_img%kill
        call img_backgr%kill
        call img_backgr_wctf%kill
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_BACKGR_SUBTR NORMAL STOP ****')
    end subroutine exec_trajectory_backgr_subtr

    subroutine exec_graphene_subtr( self, cline )
        use simple_tseries_graphene_subtr
        class(commander_graphene_subtr), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)   :: params
        type(builder)      :: build
        type(image)        :: ave_pre, ave_post, img_tmp, img_tmp2, img_tmp3
        real,  allocatable :: angles1(:), angles2(:)
        real               :: smpd, ave,var,sdev
        integer            :: iptcl, ldim_ptcl(3), ldim(3), n, nptcls
        logical            :: err
        call cline%set('objfun','cc')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        ! sanity checks & dimensions
        call find_ldim_nptcls(params%stk,ldim_ptcl,nptcls,smpd=smpd)
        if( .not.cline%defined('smpd') )then
            if( smpd < 1.e-4 ) THROW_HARD('Please provide SMPD!')
            params%smpd = smpd
        endif
        ldim_ptcl(3) = 1
        call find_ldim_nptcls(params%stk2,ldim,n)
        ldim(3) = 1
        if( any(ldim-ldim_ptcl/=0) )THROW_HARD('Inconsistent dimensions between stacks')
        if( n /= nptcls )THROW_HARD('Inconsistent number of images between stacks')
        if( .not.cline%defined('outstk') )then
            params%outstk = add2fbody(basename(params%stk),params%ext,'_subtr')
        endif
        ! initialize subtracter
        allocate(angles1(nptcls),angles2(nptcls),source=0.)
        call init_graphene_subtr(params%box,params%smpd)
        ! init images
        call build%img%new(ldim,params%smpd)       ! particle
        call img_tmp%new(ldim,params%smpd)
        call img_tmp2%new(ldim,params%smpd)   ! nn background
        call img_tmp3%new(ldim,params%smpd)
        call ave_pre%new(ldim,params%smpd)
        call ave_post%new(ldim,params%smpd)
        ave_pre  = 0.
        ave_post = 0.
        ! read, subtract & write
        do iptcl = 1,nptcls
            if( mod(iptcl,50)==0 ) call progress(iptcl,nptcls)
            ! particle
            call build%img%read(params%stk,iptcl)
            call build%img%norm()
            ! neighbours background
            call img_tmp2%read(params%stk2,iptcl)
            ! detection
            call calc_peaks(img_tmp2, angles1(iptcl), angles2(iptcl))
            ! pre-subtraction average
            call img_tmp3%copy(build%img)
            call img_tmp3%zero_edgeavg
            call img_tmp3%fft()
            call img_tmp3%ft2img('sqrt', img_tmp)
            call ave_pre%add(img_tmp)
            ! subtraction
            call remove_lattices(build%img, angles1(iptcl), angles2(iptcl))
            call build%img%norm()
            call build%img%write(params%outstk, iptcl)
            ! graphene subtracted average
            call build%img%zero_edgeavg
            call build%img%fft()
            call build%img%ft2img('sqrt', img_tmp)
            call ave_post%add(img_tmp)
        enddo
        call progress(iptcl,nptcls)
        call ave_pre%div(real(nptcls))
        call ave_post%div(real(nptcls))
        call ave_pre%write(string('pre_subtr_ave_pspec.mrc'))
        call ave_post%write(string('subtr_ave_pspec.mrc'))
        ! stats
        call moment(angles1, ave, sdev, var, err)
        write(logfhandle,'(A,F6.2,A2,F6.2,A1)')'>>> POSITION GRAPHENE SHEET 1 (DEGREES): ',ave,' (',sdev,')'
        call moment(angles2, ave, sdev, var, err)
        write(logfhandle,'(A,F6.2,A2,F6.2,A1)')'>>> POSITION GRAPHENE SHEET 2 (DEGREES): ',ave,' (',sdev,')'
        angles1 = angles2 - angles1
        where(angles1>30.) angles1 = -(angles1 - 60.)
        call moment(angles1, ave, sdev, var, err)
        write(logfhandle,'(A,F6.2,A2,F6.2,A1)')'>>> RELATIVE ROTATION (DEGREES): ',ave,' (',sdev,')'
        ! cleanup
        call build%kill_general_tbox
        call kill_graphene_subtr
        call img_tmp%kill
        call img_tmp2%kill
        call img_tmp3%kill
        call ave_pre%kill
        call ave_post%kill
        ! end gracefully
        call simple_end('**** SIMPLE_GRAPHENE_SUBTR NORMAL STOP ****')
    end subroutine exec_graphene_subtr

    subroutine exec_import_trajectory( self, cline )
        use simple_ctf_estimate_fit, only: ctf_estimate_fit
        class(commander_import_trajectory), intent(inout) :: self
        class(cmdline),                            intent(inout) :: cline
        type(parameters)       :: params
        type(sp_project)       :: spproj
        type(ctfparams)        :: ctfvars
        type(ctf_estimate_fit) :: ctffit
        type(oris)             :: os
        integer :: iframe, lfoo(3)
        call cline%set('oritype','mic')
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call spproj%read(params%projfile)
        ! # of particles
        call find_ldim_nptcls(params%stk, lfoo, params%nptcls)
        call os%new(params%nptcls, is_ptcl=.true.)
        ! CTF parameters
        ctfvars%smpd    = params%smpd
        ctfvars%phshift = 0.
        ctfvars%dfx     = 0.
        ctfvars%dfy     = 0.
        ctfvars%angast  = 0.
        ctfvars%l_phaseplate = .false.
        if( cline%defined('deftab') )then
            call ctffit%read_doc(params%deftab)
            call ctffit%get_parms(ctfvars)
            if( .not.is_equal(ctfvars%smpd,params%smpd) )then
                THROW_HARD('Iconsistent sampling distance; exec_import_trajectory')
            endif
            ctfvars%ctfflag = CTFFLAG_YES
            do iframe = 1,spproj%os_mic%get_noris()
                call spproj%os_mic%set(iframe,'ctf','yes')
            enddo
            call os%set_all2single('dfx', ctfvars%dfx)
            call os%set_all2single('dfy', ctfvars%dfy)
            call os%set_all2single('angast', ctfvars%angast)
        else
            ctfvars%ctfflag = CTFFLAG_NO
            do iframe = 1,spproj%os_mic%get_noris()
                call spproj%os_mic%set(iframe,'ctf','no')
            enddo
        endif
        ! import stack
        call os%set_all2single('state', 1.0)
        call spproj%add_single_stk(params%stk, ctfvars, os)
        call spproj%write
        ! end gracefully
        call simple_end('**** TSERIES_IMPORT_PARTICLES NORMAL STOP ****')
    end subroutine exec_import_trajectory

    subroutine exec_denoise_trajectory( self, cline )
        use simple_commanders_imgops, only: commander_ppca_denoise
        class(commander_denoise_trajectory), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(commander_ppca_denoise) :: xkpca_den
        if( .not. cline%defined('neigs')    ) call cline%set('neigs', 500)
        if( .not. cline%defined('pca_mode') ) call cline%set('pca_mode', 'kpca')
        call xkpca_den%execute(cline)
    end subroutine exec_denoise_trajectory

    subroutine exec_trajectory_swap_stack( self, cline )
        class(commander_trajectory_swap_stack), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(sp_project) :: spproj, spproj_tmp
        type(parameters) :: params
        type(ctfparams)  :: ctfparms
        integer :: ldim(3), nimgs, nstks
        call cline%set('oritype','stk')
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call spproj%read(params%projfile)
        nstks = spproj%os_stk%get_noris()
        if( nstks < 1 ) THROW_HARD('No stack could be detected in the project!')
        call find_ldim_nptcls(params%stk, ldim, nimgs)
        ldim(3) = 1
        if( spproj%get_box() /= ldim(1) .or. spproj%get_box() /= ldim(2))then
            THROW_HARD('Incompatible dimensions between stacks')
        endif
        if( nimgs /= spproj%os_ptcl2D%get_noris() ) THROW_HARD('Incompatible number of images and orientation parameters!')
        ctfparms = spproj%get_ctfparams(params%oritype, 1)
        call spproj_tmp%read(params%projfile)
        call spproj%os_stk%kill
        call spproj%os_ptcl2D%kill
        call spproj%os_ptcl3D%kill
        call spproj%add_stk(params%stk, ctfparms)
        spproj%os_ptcl2D = spproj_tmp%os_ptcl2D
        spproj%os_ptcl3D = spproj_tmp%os_ptcl3D
        call spproj_tmp%kill
        if( nstks > 1 )then
            call spproj%os_ptcl2D%set_all2single('stkind',1)
            call spproj%os_ptcl3D%set_all2single('stkind',1)
        endif
        call spproj%write(params%projfile)
        call simple_end('**** SINGLE_trajectory_swap_stack NORMAL STOP ****')
    end subroutine exec_trajectory_swap_stack

    subroutine exec_extract_substk( self, cline )
        class(commander_extract_substk), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        call cline%set('mkdir', 'no')
        ! init params
        call params%new(cline)
        ! read the project file
        call spproj%read(params%projfile)
        call spproj%write_segment_inside('projinfo')
        call spproj%write_substk([params%fromp,params%top], params%outstk)
        ! end gracefully
        call simple_end('**** SINGLE_EXTRACT_SUBSTK NORMAL STOP ****')
    end subroutine exec_extract_substk

end module single_commanders_trajectory
