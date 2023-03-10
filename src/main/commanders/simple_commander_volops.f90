! concrete commander: operations on volumes
module simple_commander_volops
include 'simple_lib.f08'
use simple_binoris_io
use simple_parameters,     only: parameters, params_glob
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_image,          only: image
use simple_projector_hlev, only: reproject, rotvol
use simple_masker,         only: masker
use simple_projector,      only: projector
use simple_volprep,        only: read_and_prep_vol
use simple_volpft_srch
implicit none

public :: centervol_commander
public :: postprocess_commander
public :: reproject_commander
public :: volops_commander
public :: dock_volpair_commander
public :: symaxis_search_commander
public :: symmetrize_map_commander
public :: symmetry_test_commander

private
#include "simple_local_flags.inc"

type, extends(commander_base) :: centervol_commander
  contains
    procedure :: execute      => exec_centervol
end type centervol_commander

type, extends(commander_base) :: postprocess_commander
 contains
   procedure :: execute      => exec_postprocess
end type postprocess_commander

type, extends(commander_base) :: reproject_commander
 contains
   procedure :: execute      => exec_reproject
end type reproject_commander

type, extends(commander_base) :: volops_commander
  contains
    procedure :: execute      => exec_volops
end type volops_commander

type, extends(commander_base) :: dock_volpair_commander
  contains
    procedure :: execute      => exec_dock_volpair
end type dock_volpair_commander

type, extends(commander_base) :: symaxis_search_commander
  contains
    procedure :: execute      => exec_symaxis_search
end type symaxis_search_commander

type, extends(commander_base) :: symmetrize_map_commander
  contains
    procedure :: execute      => exec_symmetrize_map
end type symmetrize_map_commander

type, extends(commander_base) :: symmetry_test_commander
  contains
    procedure :: execute      => exec_symmetry_test
end type symmetry_test_commander

contains

    !> centers a 3D volume and associated particle document
    subroutine exec_centervol( self, cline )
        class(centervol_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters)  :: params
        type(builder)     :: build
        real, allocatable :: shvec(:,:)
        integer           :: istate
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp',   20.)
        call build%init_params_and_build_general_tbox(cline,params)
        ! center volume(s)
        allocate(shvec(params%nstates,3))
        do istate=1,params%nstates
            call build%vol%read(params%vols(istate))
            shvec(istate,:) = build%vol%calc_shiftcen(params%cenlp, params%msk)
            call build%vol%shift([shvec(istate,1),shvec(istate,2),shvec(istate,3)])
            call build%vol%write('shifted_vol_state'//int2str(istate)//params%ext)
            ! transfer the 3D shifts to 2D
            if( cline%defined('oritab') ) call build%spproj_field%map3dshift22d(-shvec(istate,:), state=istate)
        end do
        if( cline%defined('oritab') ) call build%spproj_field%write(params%outfile, [1,build%spproj_field%get_noris()])
        ! end gracefully
        call simple_end('**** SIMPLE_CENTER NORMAL STOP ****')
    end subroutine exec_centervol

    subroutine exec_postprocess( self, cline )
        use simple_sp_project, only: sp_project
        use simple_opt_filter, only: nonuni_filt3D
        class(postprocess_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        character(len=:), allocatable :: fname_vol, fname_fsc, fname_msk, fname_mirr
        character(len=:), allocatable :: fname_even, fname_odd, fname_pproc, fname_lp
        real,             allocatable :: fsc(:), optlp(:), res(:)
        type(parameters) :: params
        type(image)      :: vol, vol_copy, vol_even, vol_odd
        type(masker)     :: mskvol, sphere
        type(sp_project) :: spproj
        real    :: fsc0143, fsc05, smpd, mskfile_smpd, lplim
        integer :: state, box, fsc_box, mskfile_box, ldim(3), flims(3,2)
        logical :: has_fsc, has_mskfile
        ! set defaults
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',     'yes')
        if( .not. cline%defined('nonuniform') ) call cline%set('nonuniform', 'no')
        call cline%set('oritype', 'out')
        ! parse commad-line
        call params%new(cline)
        ! read project segment
        call spproj%read_segment(params%oritype, params%projfile)
        ! state
        if( cline%defined('state') )then
            state = params%state
        else
            state = 1
        endif
        ! check volume, get correct smpd & box
        if( cline%defined('imgkind') )then
            call spproj%get_vol(params%imgkind, state, fname_vol, smpd, box)
        else
            call spproj%get_vol('vol', state, fname_vol, smpd, box)
        endif
        if( .not.file_exists(fname_vol) )then
            THROW_HARD('volume: '//trim(fname_vol)//' does not exist')
        endif
        ! generate file names
        fname_even  = add2fbody(trim(fname_vol), params%ext, '_even')
        fname_odd   = add2fbody(trim(fname_vol), params%ext, '_odd' )
        fname_pproc = basename(add2fbody(trim(fname_vol),   params%ext, PPROC_SUFFIX))
        fname_mirr  = basename(add2fbody(trim(fname_pproc), params%ext, MIRR_SUFFIX))
        fname_lp    = basename(add2fbody(trim(fname_vol),   params%ext, LP_SUFFIX))
        ! read volume(s)
        ldim = [box,box,box]
        call vol%new(ldim, smpd)
        call vol%read(fname_vol)
        if( params%l_nonuniform )then
            if( .not. file_exists(trim(fname_even)) ) THROW_HARD('volume: '//trim(fname_even)//' does not exists')
            if( .not. file_exists(trim(fname_odd))  ) THROW_HARD('volume: '//trim(fname_odd)//' does not exists')
            call vol_even%new(ldim, smpd)
            call vol_even%read(fname_even)
            call vol_odd%new(ldim, smpd)
            call vol_odd%read(fname_odd)
        endif
        ! check fsc filter & determine resolution
        has_fsc = .false.
        if( cline%defined('lp') )then
            lplim = params%lp
        else
            if( .not.cline%defined('fsc') )then
                call spproj%get_fsc(state, fname_fsc, fsc_box)
                params%fsc = trim(fname_fsc)
            endif
            if( .not.file_exists(params%fsc) )then
                THROW_HARD('FSC file: '//trim(params%fsc)//' not found')
            else
                has_fsc = .true.
            endif
        endif
        if( has_fsc )then
            ! resolution & optimal low-pass filter from FSC
            res   = vol%get_res()
            fsc   = file2rarr(params%fsc)
            optlp = fsc2optlp(fsc)
            call get_resolution( fsc, res, fsc05, fsc0143 )
            where( res < TINY ) optlp = 0.
            lplim = fsc0143
        endif
        if( cline%defined('bfac') )then
            ! already in params%bfac
        else
            if( lplim < 5. )then
                params%bfac = vol%guinier_bfac(HPLIM_GUINIER, lplim)
                write(logfhandle,'(A,1X,F8.2)') '>>> B-FACTOR DETERMINED TO:', params%bfac
            else
                params%bfac = 0.
            endif
        endif
        ! check volume mask
        has_mskfile = params%l_filemsk
        if( .not. has_mskfile )then
            call spproj%get_vol('vol_msk', 1, fname_msk, mskfile_smpd, mskfile_box)
            params%mskfile = trim(fname_msk)
            if( file_exists(params%mskfile) ) has_mskfile = .true.
        endif
        if( has_mskfile )then
            call mskvol%new(ldim, smpd)
            call mskvol%read(params%mskfile)
        endif
        ! B-factor
        call vol%fft()
        call vol_copy%copy(vol)
        call vol%apply_bfac(params%bfac)
        ! low-pass filter
        if( cline%defined('lp') )then
            ! low-pass overrides all input
            call vol%bp(0., params%lp)
            call vol_copy%bp(0., params%lp)
        else if( params%l_nonuniform .and. has_fsc )then
            call sphere%new(ldim, smpd)
            sphere = 1.0
            call sphere%mask(params%msk, 'soft', backgr=0.0)
            call sphere%one_at_edge
            call nonuni_filt3D(vol_odd, vol_even, sphere)
            call sphere%kill
            ! merge volumes
            call vol_odd%add(vol_even)
            call vol_odd%mul(0.5)
            vol_copy = vol_odd
            vol      = vol_odd
            call vol%apply_bfac(params%bfac)
            ! final low-pass filtering for smoothness
            call vol%bp(0., fsc0143)
            call vol_copy%bp(0., fsc0143)
        else if( has_fsc )then
            ! optimal low-pass filter of unfiltered volumes from FSC
            if( params%cc_objfun == OBJFUN_CC )then
                call vol%apply_filter(optlp)
                call vol_copy%apply_filter(optlp)
            else
                ! FSC-based regularization is performed upon assembly
            endif
            ! final low-pass filtering for smoothness
            call vol%bp(0., fsc0143)
            call vol_copy%bp(0., fsc0143)
        else
            THROW_HARD('no method for low-pass filtering defined; give fsc|lp on command line; exec_postprocess')
        endif
        ! write low-pass filtered without B-factor or mask & read the original back in
        call vol_copy%ifft
        call vol_copy%write(fname_lp)
        call vol_copy%read(fname_vol)
        ! mask
        call vol%ifft()
        if( params%l_automsk ) has_mskfile = .false. ! turn off masking using file
        if( params%l_automsk )then
            if( .not. cline%defined('thres') .and. cline%defined('mw') )then
                write(logfhandle,*) 'Need a pixel threshold > 0. for the binarisation'
                write(logfhandle,*) 'Procedure for obtaining thresh:'
                write(logfhandle,*) '(1) Use UCSF Chimera to look at the *_lp.mrc volume'
                write(logfhandle,*) '(2) Identify the pixel threshold that excludes any background noise'
                THROW_HARD('postprocess')
            endif
            if( cline%defined('thres') .and. .not. cline%defined('mw') )then
                THROW_HARD('molecular weight must be provided for threshold-based auto-masking; postprocess')
            endif
            if( cline%defined('thres') .and. cline%defined('mw') )then
                call mskvol%automask3D(vol_copy)
            else
                call mskvol%automask3D_otsu(vol_copy)
            endif
            call mskvol%write('automask'//params%ext)
            call vol%zero_background
            call vol%mul(mskvol)
            call mskvol%kill
        else if( has_mskfile )then
            call vol%zero_background
            if( cline%defined('lp_backgr') )then
                call vol%lp_background(mskvol,params%lp_backgr)
            else
                call vol%mul(mskvol)
            endif
            call mskvol%kill
        else
            call vol%mask(params%msk, 'soft')
        endif
        ! output in cwd
        call vol%write(fname_pproc)
        ! also output mirrored by default (unless otherwise stated on command line)
        if( .not. cline%defined('mirr') .or. params%mirr .ne. 'no' )then
            call vol%mirror('x')
            call vol%write(fname_mirr)
        endif
        ! destruct
        call spproj%kill
        call vol%kill
        call vol_copy%kill
        call simple_end('**** SIMPLE_POSTPROCESS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_postprocess

    !> exec_project generate projections from volume
    subroutine exec_reproject( self, cline )
        class(reproject_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters)         :: params
        type(builder)            :: build
        type(image), allocatable :: imgs(:)
        integer                  :: i
        logical                  :: do_zero
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',  'no')
        if( .not. cline%defined('wfun')   ) call cline%set('wfun',   'kb')
        if( .not. cline%defined('winsz')  ) call cline%set('winsz',   1.5)
        if( .not. cline%defined('alpha')  ) call cline%set('alpha',    2.)
        if( .not. cline%defined('outstk') ) call cline%set('outstk', 'reprojs'//trim(STK_EXT))
        if( .not. cline%defined('oritab') )then
            if( .not. cline%defined('nspace') ) THROW_HARD('need nspace (for number of projections)!')
        endif
        call params%new(cline)
        if( cline%defined('oritab') )then
            params%nptcls = binread_nlines(params%oritab)
            call build%build_spproj(params, cline)
            call build%build_general_tbox(params, cline)
            params%nspace = build%spproj_field%get_noris()
        else
            params%nptcls = params%nspace
            call build%build_spproj(params, cline)
            call build%build_general_tbox(params, cline)
            call build%pgrpsyms%build_refspiral(build%spproj_field)
        endif
        ! fix volumes and stacks
        call build%vol%read(params%vols(1))
        ! masking
        if(cline%defined('mskdiam')) call build%vol%mask(params%msk, 'soft')
        ! generate projections
        imgs     = reproject(build%vol, build%spproj_field)
        if( file_exists(params%outstk) ) call del_file(params%outstk)
        do_zero  = build%spproj_field%isthere('state')
        do i=1,params%nspace
            if( params%neg .eq. 'yes' ) call imgs(i)%neg()
            if( do_zero )then
                if( build%spproj_field%get_state(i) < 1 ) call imgs(i)%zero
            endif
            call imgs(i)%write(params%outstk,i)
        end do
        call build%spproj_field%write('reproject_oris'//trim(TXT_EXT), [1,params%nptcls])
        call simple_end('**** SIMPLE_REPROJECT NORMAL STOP ****')
    end subroutine exec_reproject

    !> volume calculations and operations - incl Guinier, snr, mirror or b-factor
    subroutine exec_volops( self, cline )
        class(volops_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        logical          :: here
        type(ori)        :: o
        type(image)      :: vol_rot
        real             :: shvec(3),  ave, sdev, maxv, minv
        call build%init_params_and_build_general_tbox(cline,params)
        inquire(FILE=params%vols(1), EXIST=here)
        if( here )then
            call build%vol%read(params%vols(1))
        else
            THROW_HARD('vol1 does not exists in cwd')
        endif
        if( params%stats .eq. 'yes' )then
            call build%vol%stats('foreground', ave, sdev, maxv, minv)
            write(logfhandle,*) 'ave : ', ave
            write(logfhandle,*) 'sdev: ', sdev
            write(logfhandle,*) 'maxv: ', maxv
            write(logfhandle,*) 'minv: ', minv
            return
        endif
        if( params%guinier .eq. 'yes' )then
            if( .not. cline%defined('smpd') ) THROW_HARD('need smpd (sampling distance) input for Guinier plot')
            if( .not. cline%defined('hp')   ) THROW_HARD('need hp (high-pass limit) input for Guinier plot')
            if( .not. cline%defined('lp')   ) THROW_HARD('need lp (low-pass limit) input for Guinier plot')
            params%bfac = build%vol%guinier_bfac(params%hp, params%lp)
            write(logfhandle,'(A,1X,F8.2)') '>>> B-FACTOR DETERMINED TO:', params%bfac
        else
            if( cline%defined('mul') )then
                call build%vol%mul(params%mul)
            end if
            if( cline%defined('neg')  )then
                call build%vol%neg()
            end if
            if( cline%defined('snr') )then
                call build%vol%add_gauran(params%snr)
            end if
            if( cline%defined('mirr') )then
                call build%vol%mirror(params%mirr)
            end if
            if( cline%defined('bfac') )then
                call build%vol%apply_bfac(params%bfac)
            end if
            if( cline%defined('e1') .or. cline%defined('e2') .or. cline%defined('e3') )then
                if( .not. cline%defined('smpd') ) THROW_HARD('need smpd (sampling distance) input for volume rotation')
                if( .not. cline%defined('mskdiam')  ) THROW_HARD('need mskdiam (mask diameter in A) input for volume rotation')
                call o%new(is_ptcl=.false.)
                call o%set_euler([params%e1,params%e2,params%e3])
                shvec     = [params%xsh,params%ysh,params%zsh]
                vol_rot   = rotvol(build%vol, o,  shvec)
                build%vol = vol_rot
            endif
            if( cline%defined('xsh') .or. cline%defined('ysh') .or. cline%defined('zsh') )then
                call build%vol%shift([params%xsh,params%ysh,params%zsh])
            endif
            if( cline%defined('lp_backgr') .and. params%l_filemsk )then
                call build%vol2%new(build%vol%get_ldim(), build%vol%get_smpd())
                call build%vol2%read(params%mskfile)
                call build%vol%lp_background(build%vol2, params%lp_backgr)
            endif
            call build%vol%write(params%outvol, del_if_exists=.true.)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_VOLOPS NORMAL STOP ****')
    end subroutine exec_volops

    subroutine exec_dock_volpair( self, cline )
        use simple_vol_srch
        class(dock_volpair_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        real,  parameter :: SHSRCH_HWDTH  = 5.0
        type(parameters) :: params
        type(projector)  :: vol1, vol2
        ! type(image)      :: vol_out
        type(projector)      :: vol_out
        type(oris)       :: ori2write
        type(ori)        :: orientation, orientation_best
        real             :: cxyz(4), cxyz2(4), corr_prev
        integer          :: i
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes'       )
        if( .not. cline%defined('trs'  ) ) call cline%set('trs',   SHSRCH_HWDTH)
        call params%new(cline)
        ! prep vols
        call read_and_prep_vol(params%vols(1), vol1)
        call read_and_prep_vol(params%vols(2), vol2)
        select case( trim(params%dockmode) )
        case('shift')
                call vol_srch_init(vol1, vol2, params%hp, params%lpstart, params%trs)
                cxyz = vol_shsrch_minimize()
                write(logfhandle,*) 'corr: ', cxyz(1), 'shift: ', cxyz(2:4)
                call vol2%read(params%vols(2))
                call vol2%shift(cxyz(2:4))
                call vol2%write(params%outvol)
            case('rot')
                ! to be implemented
            case('rotshift')
                ! grid search with volpft (icosahedral sampling geometry) using lpstart low-pass limit
                call volpft_srch_init(vol1, vol2, params%hp, params%lpstart)
                orientation = volpft_srch_minimize()
                ! rotate vol to create reference for shift alignment
                call vol2%ifft
                vol_out = rotvol(vol2, orientation)
                call vol2%fft
                call vol_out%fft
                ! continuous shift alignment using lpstart low-pass limit
                call vol_srch_init(vol1, vol_out, params%hp, params%lpstart, params%trs)
                cxyz = vol_shsrch_minimize()
                ! re-search the angular grid with the shifts in-place
                call volpft_srch_set_shvec(cxyz(2:4))
                orientation = volpft_srch_minimize()
                ! Refinment using lpstop low-pass limit
                call volpft_srch_init(vol1, vol_out, params%hp, params%lpstop)
                ! call volpft_srch_init(vol1, vol2, params%hp, params%lpstop)
                call vol_srch_init(vol1, vol_out, params%hp, params%lpstop, params%trs)
                corr_prev = 0. ! initialise
                cxyz2     = 0.
                do i=1,10
                    ! rotate and shift vol to create reference for shift alignment
                    call vol2%ifft
                    vol_out = rotvol(vol2, orientation, cxyz(2:4))
                    call vol2%fft
                    call vol_out%fft
                    cxyz2 = vol_shsrch_minimize()
                    if( any(abs(cxyz2(2:4)) > 0.) )then ! a better solution was found
                        corr_prev = cxyz2(1)
                        ! obtain joint shifts by vector addition
                        cxyz(2:4) = cxyz(2:4) + cxyz2(2:4)
                        ! update shifts in volpft_srch class and refine angles
                        call volpft_srch_set_shvec(cxyz(2:4))
                        if( i <= 3 )then
                            orientation_best = volpft_srch_refine(orientation)
                        else
                            orientation_best = volpft_srch_refine(orientation, angres=5.)
                        endif
                        orientation = orientation_best
                        call orientation%set('x', cxyz(2))
                        call orientation%set('y', cxyz(3))
                        call orientation%set('z', cxyz(4))
                        call orientation%print_ori
                    endif
                end do
                ! rotate and shift vol for output
                call vol2%ifft
                vol_out = rotvol(vol2, orientation, cxyz(2:4))
                call orientation%print_ori
                if( cline%defined('outfile') )then
                    call ori2write%new(1, is_ptcl=.false.)
                    call ori2write%set_ori(1,orientation)
                    call ori2write%write(params%outfile)
                    call ori2write%kill
                endif
                ! write
                call vol_out%write(params%outvol, del_if_exists=.true.)
                ! destruct
                call volpft_srch_kill
            case('refine')
                ! to be implemented
            case DEFAULT
                write(logfhandle,*) 'dockmode: ', trim(params%dockmode), ' is unsupported'
        end select
        ! cleanup
        call vol1%kill
        call vol2%kill
        call vol_out%kill
        call orientation%kill
        call orientation_best%kill
        ! end gracefully
        call simple_end('**** SIMPLE_DOCK_VOLPAIR NORMAL STOP ****')
    end subroutine exec_dock_volpair

    !> for identification of the principal symmetry axis
    subroutine exec_symaxis_search( self, cline )
        use simple_volpft_symsrch
        class(symaxis_search_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        character(len=32), parameter :: SYMSHTAB   = 'sym_3dshift'//trim(TXT_EXT)
        character(len=32), parameter :: SYMAXISTAB = 'sym_axis'//trim(TXT_EXT)
        type(parameters)      :: params
        type(builder)         :: build
        type(ori)             :: symaxis
        type(oris)            :: symaxis4write
        type(sym)             :: syme
        character(len=STDLEN) :: fbody
        real                  :: shvec(3)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('cenlp')   ) call cline%set('cenlp', 20.)
        if( .not. cline%defined('center')  ) call cline%set('center', 'yes')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'cls3D')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        call build%vol%read(params%vols(1))
        ! center volume
        shvec = 0.
        if( params%center.eq.'yes' )then
            shvec = build%vol%calc_shiftcen(params%cenlp,params%msk)
            call build%vol%shift(shvec)
            fbody = get_fbody(params%vols(1),fname2ext(params%vols(1)))
            call build%vol%write(trim(fbody)//'_centered.mrc')
        endif
        ! mask volume
        call build%vol%mask(params%msk, 'soft')
        ! init search object
        call volpft_symsrch_init(build%vol, params%pgrp, params%hp, params%lp)
        ! search
        call volpft_srch4symaxis(symaxis)
        call symaxis4write%new(1, is_ptcl=.false.)
        call symaxis4write%set_ori(1, symaxis)
        call symaxis4write%write(SYMAXISTAB, [1,1])
        if( cline%defined('projfile') )then
            call build%spproj%read(params%projfile)
            if( .not. build%spproj%is_virgin_field(params%oritype) )then
                ! transfer shift and symmetry to orientations
                call syme%new(params%pgrp)
                ! rotate the orientations & transfer the 3d shifts to 2d
                shvec = -1. * shvec ! the sign is right
                call syme%apply_sym_with_shift(build%spproj_field, symaxis, shvec)
                call build%spproj%write_segment_inside(params%oritype, params%projfile)
            endif
        endif
        ! rotate volume to symaxis for later comparison
        call build%vol%read(params%vols(1))
        build%vol2 = rotvol(build%vol, symaxis)
        call build%vol2%write('vol_aligned2_'//trim(params%pgrp)//'axis'//params%ext)
        ! destruct
        call symaxis%kill
        call symaxis4write%kill
        call syme%kill
        ! end gracefully
        call simple_end('**** SIMPLE_SYMAXIS_SEARCH NORMAL STOP ****')
        call simple_touch('SYMAXIS_SEARCH_FINISHED', errmsg='In: commander_volops :: exec_symsrch')
    end subroutine exec_symaxis_search

    subroutine exec_symmetrize_map( self, cline )
        use simple_symanalyzer
        class(symmetrize_map_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)   :: params
        type(builder)      :: build
        real               :: shvec(3), scale, smpd
        integer            :: ldim(3)
        integer, parameter :: MAXBOX = 128
        character(len=:), allocatable :: fbody
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',  'yes')
        if( .not. cline%defined('cenlp')  ) call cline%set('cenlp',    20.)
        if( .not. cline%defined('center') ) call cline%set('center', 'yes')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        if( .not. cline%defined('outvol') ) params%outvol = 'symmetrized_map'//params%ext
        call build%vol%read(params%vols(1))
        ! possible downscaling of input vol
        ldim = build%vol%get_ldim()
        scale = 1.
        if( ldim(1) > MAXBOX )then
            scale = real(MAXBOX) / real(ldim(1))
            call build%vol%fft
            call build%vol%clip_inplace([MAXBOX,MAXBOX,MAXBOX])
            call build%vol%ifft
            smpd         = build%vol%get_smpd()
            params%msk   = round2even(scale * params%msk)
            params%width = scale * params%width
        endif
        if( params%center.eq.'yes' )then
            shvec = build%vol%calc_shiftcen(params%cenlp,params%msk)
            call build%vol%shift(shvec)
            ! store un-scaled shift parameters
            params%xsh = shvec(1) / scale
            params%ysh = shvec(1) / scale
            params%zsh = shvec(1) / scale
        else
            params%xsh = 0.
            params%ysh = 0.
            params%zsh = 0.
        endif
        ! mask volume
        call build%vol%mask(params%msk, 'soft')
        ! symmetrize
        call symmetrize_map(build%vol, params, build%vol2)
        ! write
        call build%vol2%write(trim(params%outvol))
        ! end gracefully
        call simple_end('**** SIMPLE_SYMMETRIZE_MAP NORMAL STOP ****')
    end subroutine exec_symmetrize_map

    subroutine exec_symmetry_test( self, cline )
        use simple_symanalyzer
        class(symmetry_test_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)      :: params
        type(builder)         :: build
        character(len=STDLEN) :: fbody
        character(len=3)      :: pgrp
        real                  :: shvec(3), scale, smpd
        integer               :: ldim(3)
        integer, parameter    :: MAXBOX = 128
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',  'yes')
        if( .not. cline%defined('cenlp')  ) call cline%set('cenlp',    20.)
        if( .not. cline%defined('center') ) call cline%set('center', 'yes')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        call build%vol%read(params%vols(1))
        ! possible downscaling of input vol
        ldim = build%vol%get_ldim()
        scale = 1.
        if( ldim(1) > MAXBOX )then
            scale = real(MAXBOX) / real(ldim(1))
            call build%vol%fft
            call build%vol%clip_inplace([MAXBOX,MAXBOX,MAXBOX])
            call build%vol%ifft
            smpd         = build%vol%get_smpd()
            params%msk   = round2even(scale * params%msk)
            params%width = scale * params%width
        endif
        ! low-pass limit safety
        params%lp = max(2. * smpd, params%lp)
        ! centering
        shvec = 0.
        if( params%center.eq.'yes' )then
            shvec = build%vol%calc_shiftcen(params%cenlp,params%msk)
            call build%vol%shift(shvec)
        endif
        ! mask volume
        call build%vol%mask(params%msk, 'soft')
        ! run test
        if(cline%defined('fname')) then
          call symmetry_tester(build%vol, params%msk, params%hp,&
          &params%lp, params%cn_stop, params%platonic .eq. 'yes', pgrp, params%fname)
        else
          call symmetry_tester(build%vol, params%msk, params%hp,&
            &params%lp, params%cn_stop, params%platonic .eq. 'yes', pgrp)
          endif
        ! end gracefully
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_SYMMETRY_TEST NORMAL STOP ****')
    end subroutine exec_symmetry_test

end module simple_commander_volops
