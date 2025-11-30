! concrete commander: operations on volumes
module simple_commanders_volops
include 'simple_lib.f08'
use simple_binoris_io
use simple_strategy2D_utils
use simple_parameters,       only: parameters, params_glob
use simple_builder,          only: builder
use simple_cmdline,          only: cmdline
use simple_commander_base,   only: commander_base
use simple_image,            only: image
use simple_simple_volinterp, only: reproject, rotvol
use simple_image_msk,        only: image_msk
use simple_projector,        only: projector
use simple_dock_vols,        only: dock_vols
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_centervol
  contains
    procedure :: execute      => exec_centervol
end type commander_centervol

type, extends(commander_base) :: commander_dock_volpair
  contains
    procedure :: execute      => exec_dock_volpair
end type commander_dock_volpair

type, extends(commander_base) :: commander_noisevol
  contains
    procedure :: execute      => exec_noisevol
end type commander_noisevol

type, extends(commander_base) :: commander_postprocess
 contains
   procedure :: execute      => exec_postprocess
end type commander_postprocess

type, extends(commander_base) :: commander_ppca_volvar
  contains
    procedure :: execute      => exec_ppca_volvar
end type commander_ppca_volvar

type, extends(commander_base) :: commander_reproject
 contains
   procedure :: execute      => exec_reproject
end type commander_reproject

type, extends(commander_base) :: commander_sharpvol
  contains
    procedure :: execute      => exec_sharpvol
end type commander_sharpvol

type, extends(commander_base) :: commander_symaxis_search
  contains
    procedure :: execute      => exec_symaxis_search
end type commander_symaxis_search

type, extends(commander_base) :: commander_symmetrize_map
  contains
    procedure :: execute      => exec_symmetrize_map
end type commander_symmetrize_map

type, extends(commander_base) :: commander_symmetry_test
  contains
    procedure :: execute      => exec_symmetry_test
end type commander_symmetry_test

type, extends(commander_base) :: commander_volanalyze
  contains
    procedure :: execute      => exec_volanalyze
end type commander_volanalyze

type, extends(commander_base) :: commander_volops
  contains
    procedure :: execute      => exec_volops
end type commander_volops

contains

    !> centers a 3D volume and associated particle document
    subroutine exec_centervol( self, cline )
        class(commander_centervol), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters)  :: params
        type(builder)     :: build
        real, allocatable :: shvec(:,:)
        real              :: shift2d(2), xyz(3)
        integer           :: ldim(3), istate, i, nimgs
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp',   20.)
        call build%init_params_and_build_general_tbox(cline,params)
        select case(trim(params%center_type))
        case('params')
            if(.not.cline%defined('oritab') .and. .not.cline%defined('projfile') )then
                THROW_HARD('Alignments parameters are required')
            endif
        end select
        if( cline%defined('stk') )then
            call find_ldim_nptcls(params%stk, ldim, nimgs)
            allocate(shvec(nimgs,3))
            call build%img%construct_thread_safe_tmp_imgs(1)
            do i = 1,nimgs
                call build%img%zero_and_unflag_ft
                call build%img%read(params%stk,i)
                select case(trim(params_glob%center_type))
                case('params')
                    xyz = 0.
                    call build%spproj_field%calc_avg_offset2D(i, xyz(1:2))
                    if( arg(xyz) < CENTHRESH )then
                        shvec(i,:) = 0.
                    else if( arg(xyz) > MAXCENTHRESH2D )then
                        shvec(i,:) = xyz
                    else
                        shvec(i,:) = build%img%calc_shiftcen_serial(params%cenlp, params%msk)
                        if( arg(shvec(i,1:2) - xyz(1:2)) > MAXCENTHRESH2D ) shvec(i,:) = 0.
                    endif
                    call build%spproj_field%add_shift2class(i, -shvec(i,:))
                case('seg')
                    call calc_cavg_offset(build%img, params%cenlp, params%msk, shift2d, i)
                    shvec(i,1:2) = shift2d
                    if( cline%defined('oritab') .or. cline%defined('projfile') )then
                        call build%spproj_field%add_shift2class(i, -shift2d)
                    endif
                case('mass')
                    shvec(i,:) = build%img%calc_shiftcen_serial(params%cenlp, params%msk)
                    if( cline%defined('oritab') .or. cline%defined('projfile') )then
                        call build%spproj_field%add_shift2class(i, -shvec(i,:))
                    endif
                end select
                write(logfhandle,'(A,I4,2F6.1)')'>>> OFFSET: ',i,shvec(i,1:2)
                call build%img%fft
                shift2d = shvec(i,1:2)
                call build%img%shift2Dserial(shift2d)
                call build%img%ifft
                call build%img%write(params%outstk,i)
            enddo
            call build%img%kill_thread_safe_tmp_imgs
        else
            ! center volume(s)
            allocate(shvec(params%nstates,3))
            do istate=1,params%nstates
                call build%vol%read(params%vols(istate))
                select case(trim(params_glob%center_type))
                case('mass')
                    shvec(istate,:) = build%vol%calc_shiftcen(params%cenlp, params%msk)
                    if( cline%defined('oritab') .or. cline%defined('projfile') )then
                        call build%spproj_field%map3dshift22d(-shvec(istate,:), state=istate)
                    endif
                case('params')
                    call build%spproj_field%calc_avg_offset3D(shvec(istate,:), state=istate)
                    call build%spproj_field%map3dshift22d(-shvec(istate,:), state=istate)
                end select
                write(logfhandle,'(A,3F6.1)')'>>> OFFSET: ',shvec(istate,:)
                call build%vol%shift(shvec(istate,1:3))
                call build%vol%write(string('shifted_vol_state'//int2str(istate)//params%ext%to_char()))
            end do
        endif
        if( cline%defined('oritab') )then
            call build%spproj_field%write(params%outfile, [1,build%spproj_field%get_noris()])
        else if( cline%defined('projfile') )then
            params%projfile = basename(params%projfile)
            call build%spproj%write(params%projfile)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CENTER NORMAL STOP ****')
    end subroutine exec_centervol

    ! Sharpening protocol for raw volumes using even/odd volume pairs
    subroutine exec_sharpvol( self, cline )
        !use simple_sp_project, only: sp_project
        use simple_atoms,        only: atoms
        use simple_segmentation, only: otsu_img
        use simple_image_bin,    only: image_bin
        class(commander_sharpvol), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        real,             allocatable :: fsc(:), optlp(:), res(:)
        type(string)          :: fname_vol, fname_even, fname_odd, fname_pdb, fname_pproc, fname_lp, fname_mirr, pdbout_fname
        type(parameters)      :: params
        type(image)           :: vol_bfac, vol_no_bfac, vol_pdb, cos_img
        type(image_bin)       :: vol_pdbin
        type(atoms)           :: pdb
        real    :: fsc0143, fsc05, lplim, lp_min=0., lp_max=20.
        integer :: ldim(3), ifoo
        logical :: has_fsc
        integer, parameter :: EDGE=6
        ! set defaults
        ! parse commad-line
        call params%new(cline)
        ! check volume, get correct smpd & box
        fname_vol = params%vols(1)
        if( .not.file_exists(fname_vol) )then
            THROW_HARD('volume: '//fname_vol%to_char()//' does not exist')
        endif
        ! generate file names
        fname_pdb   = add2fbody(fname_vol, params%ext, '_pdb')
        fname_even  = add2fbody(fname_vol, params%ext, '_even')
        fname_odd   = add2fbody(fname_vol, params%ext, '_odd' )
        fname_pproc = basename(add2fbody(fname_vol,   params%ext, PPROC_SUFFIX))
        fname_mirr  = basename(add2fbody(fname_pproc, params%ext, MIRR_SUFFIX))
        fname_lp    = basename(add2fbody(fname_vol,   params%ext, LP_SUFFIX))
        if( .not.file_exists(fname_even) )then
            THROW_HARD('volume: '//fname_even%to_char()//' does not exist')
        endif
        if( .not.file_exists(fname_odd) )then
            THROW_HARD('volume: '//fname_odd%to_char()//' does not exist')
        endif
        ! read volume(s)
        call find_ldim_nptcls(fname_vol, ldim, ifoo)
        call vol_bfac%new(ldim, params%smpd)
        call vol_bfac%read(fname_vol)
        ! generate & apply pdb mask from pdb
        call pdb%new(params%pdbfile)
        if( pdb%check_center() .or. params%center_pdb .eq. 'yes' )then
            pdbout_fname = get_fbody(params%pdbfile, 'pdb') // '_centered.pdb'
            call pdb%pdb2mrc( params%pdbfile, fname_pdb, params%smpd, center_pdb=.true., pdb_out=pdbout_fname, vol_dim=ldim )
        else 
            call pdb%pdb2mrc( params%pdbfile, fname_pdb, params%smpd, vol_dim=ldim )
        endif
        call vol_pdb%new(ldim, params%smpd)
        call vol_pdb%read(fname_pdb)
        call vol_pdb%bp(lp_min, lp_max)
        call otsu_img(vol_pdb)
        call vol_pdbin%transfer2bimg(vol_pdb)
        call vol_pdbin%grow_bins(1)
        call vol_bfac%mul(vol_pdbin)
        call vol_pdbin%cos_edge(EDGE, cos_img)
        call vol_bfac%mul(cos_img)
        ! filtered function & apply Bfactor
        ! check fsc filter & determine resolution
        has_fsc = .false.        
        if( cline%defined('fsc') )then
            if( .not.file_exists(params%fsc) ) THROW_HARD('FSC file: '//params%fsc%to_char()//' not found')
            has_fsc = .true.
        else
            if( file_exists(params%fsc) )then
                has_fsc = .true.
            else
                THROW_WARN('FSC file: '//params%fsc%to_char()//' not found')
                has_fsc = .false.
                if( .not. cline%defined('lp') ) THROW_HARD('no method for low-pass filtering defined; give fsc|lp on command line; exec_sharpvol')
            endif 
        endif
        if( has_fsc )then
            ! resolution & optimal low-pass filter from FSC
            res   = vol_bfac%get_res()
            fsc   = file2rarr(params%fsc)
            optlp = fsc2optlp(fsc)
            call get_resolution(fsc, res, fsc05, fsc0143)
            where( res < TINY ) optlp = 0.
            lplim = fsc0143
        else
            lplim = params%lp
        endif
        if( cline%defined('bfac') )then
            ! already in params%bfac
        else
            if( lplim < 5. )then
                params%bfac = vol_bfac%guinier_bfac(HPLIM_GUINIER, lplim)
                write(logfhandle,'(A,1X,F8.2)') '>>> B-FACTOR DETERMINED TO:', params%bfac
            else
                params%bfac = 0.
            endif
        endif
        ! B-factor
        call vol_bfac%fft()
        call vol_no_bfac%copy(vol_bfac)
        call vol_bfac%apply_bfac(params%bfac)
        ! low-pass filter    
        if( has_fsc )then
            ! optimal low-pass filter of unfiltered volumes from FSC
            call vol_bfac%apply_filter(optlp)
            call vol_no_bfac%apply_filter(optlp)
            ! final low-pass filtering for smoothness
            call vol_bfac%bp(0., fsc0143)
            call vol_no_bfac%bp(0., fsc0143)
        else
            call vol_bfac%bp(0., lplim)
            call vol_no_bfac%bp(0., lplim)
        endif
        ! write low-pass filtered without B-factor or mask & read the original back in
        call vol_no_bfac%ifft
        call vol_no_bfac%write(fname_lp)
        ! mask
        call vol_bfac%ifft()
        call vol_bfac%mask(params%msk_crop, 'soft')
        ! output in cwd
        call vol_bfac%write(fname_pproc)
        ! also output mirrored by default (unless otherwise stated on command line)
        if( .not. cline%defined('mirr') .or. params%mirr .ne. 'no' )then
            call vol_bfac%mirror('x')
            call vol_bfac%write(fname_mirr)
        endif
        ! destruct
        call pdb%kill
        call vol_bfac%kill
        call vol_pdb%kill
        call vol_pdbin%kill
        call vol_no_bfac%kill
        call simple_end('**** SIMPLE_SHARPVOL NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_sharpvol

    subroutine exec_postprocess( self, cline )
        use simple_sp_project, only: sp_project
        class(commander_postprocess), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        real, allocatable :: fsc(:), optlp(:), res(:)
        type(string)     :: fname_vol, fname_fsc, fname_mirr
        type(string)     :: fname_even, fname_odd, fname_pproc, fname_lp
        type(parameters) :: params
        type(image)      :: vol_bfac, vol_no_bfac
        type(sp_project) :: spproj
        real    :: fsc0143, fsc05, smpd, lplim
        integer :: state, box, fsc_box, ldim(3)
        logical :: has_fsc
        ! set defaults
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
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
            THROW_HARD('volume: '//fname_vol%to_char()//' does not exist')
        endif
        ! using the input volume for postprocessing
        if( cline%defined('vol'//int2str(state)) ) fname_vol = params%vols(state)
        ! generate file names
        fname_even  = add2fbody(fname_vol, params%ext, '_even')
        fname_odd   = add2fbody(fname_vol, params%ext, '_odd' )
        fname_pproc = basename(add2fbody(fname_vol,   params%ext, PPROC_SUFFIX))
        fname_mirr  = basename(add2fbody(fname_pproc, params%ext, MIRR_SUFFIX))
        fname_lp    = basename(add2fbody(fname_vol,   params%ext, LP_SUFFIX))
        ! read volume(s)
        ldim = [box,box,box]
        call vol_bfac%new(ldim, smpd)
        call vol_bfac%read(fname_vol)
        ! check fsc filter & determine resolution
        has_fsc = .false.        
        if( cline%defined('fsc') )then
            if( .not.file_exists(params%fsc) ) THROW_HARD('FSC file: '//params%fsc%to_char()//' not found')
            has_fsc = .true.
        else
            call spproj%get_fsc(state, fname_fsc, fsc_box)
            params%fsc = fname_fsc
            if( file_exists(params%fsc) )then
                has_fsc = .true.
            else
                THROW_WARN('FSC file: '//params%fsc%to_char()//' not found')
                has_fsc = .false.
                if( .not. cline%defined('lp') ) THROW_HARD('no method for low-pass filtering defined; give fsc|lp on command line; exec_postprocess')
            endif 
        endif
        if( has_fsc )then
            ! resolution & optimal low-pass filter from FSC
            res   = vol_bfac%get_res()
            fsc   = file2rarr(params%fsc)
            optlp = fsc2optlp(fsc)
            call get_resolution(fsc, res, fsc05, fsc0143)
            where( fsc < 0.05 ) optlp = 0.
            where( res < TINY ) optlp = 0.
            lplim = fsc0143
        else
            lplim = params%lp
        endif
        if( cline%defined('bfac') )then
            ! already in params%bfac
        else
            if( lplim < 5. )then
                params%bfac = vol_bfac%guinier_bfac(HPLIM_GUINIER, lplim)
                write(logfhandle,'(A,1X,F8.2)') '>>> B-FACTOR DETERMINED TO:', params%bfac
            else
                params%bfac = 0.
            endif
        endif
        ! B-factor
        call vol_bfac%fft()
        call vol_no_bfac%copy(vol_bfac)
        call vol_bfac%apply_bfac(params%bfac)
        ! low-pass filter    
        if( has_fsc )then
            ! optimal low-pass filter of unfiltered volumes from FSC
            call vol_bfac%apply_filter(optlp)
            call vol_no_bfac%apply_filter(optlp)
            ! final low-pass filtering for smoothness
            ! call vol_bfac%bp(0., fsc0143)
            ! call vol_no_bfac%bp(0., fsc0143)
        else
            call vol_bfac%bp(0., lplim)
            call vol_no_bfac%bp(0., lplim)
        endif
        ! write low-pass filtered without B-factor or mask & read the original back in
        call vol_no_bfac%ifft
        call vol_no_bfac%write(fname_lp)
        ! mask
        call vol_bfac%ifft()
        call vol_bfac%mask(params%msk_crop, 'soft')
        ! output in cwd
        call vol_bfac%write(fname_pproc)
        ! also output mirrored by default (unless otherwise stated on command line)
        if( .not. cline%defined('mirr') .or. params%mirr .ne. 'no' )then
            call vol_bfac%mirror('x')
            call vol_bfac%write(fname_mirr)
        endif
        ! destruct
        call spproj%kill
        call vol_bfac%kill
        call vol_no_bfac%kill
        call simple_end('**** SIMPLE_POSTPROCESS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_postprocess

    !> exec_project generate projections from volume
    subroutine exec_reproject( self, cline )
        class(commander_reproject), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters)         :: params
        type(builder)            :: build
        type(image), allocatable :: imgs(:)
        integer,     allocatable :: states(:), tmp(:)
        integer                  :: i, s
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',  'no')
        if( .not. cline%defined('wfun')   ) call cline%set('wfun',   'kb')
        if( .not. cline%defined('winsz')  ) call cline%set('winsz',   1.5)
        if( .not. cline%defined('alpha')  ) call cline%set('alpha',    2.)
        if( .not. cline%defined('outstk') ) call cline%set('outstk', 'reprojs'//STK_EXT)
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
            call build%spproj_field%set_all2single('state',1)
        endif
        ! generate projections
        if( file_exists(params%outstk) ) call del_file(params%outstk)
        states = nint(build%spproj_field%get_all('state'))
        tmp    = states
        do s = 1,params%nstates
            if( cline%defined('state') )then
                if( s /= params%state ) cycle
            endif
            if( .not.cline%defined('vol'//int2str(s)) ) cycle
            if( any(states==s) )then
                ! read and mask
                call build%vol%read(params%vols(s))
                if(cline%defined('mskdiam')) call build%vol%mask(params%msk, 'soft',backgr=0.)
                ! project
                where( states == s )
                    tmp = 1
                elsewhere
                    tmp = 0
                end where
                call build%spproj_field%set_all('state', tmp)
                imgs = reproject(build%vol, build%spproj_field)
                ! write
                do i = 1,params%nspace
                    if( states(i) /= s ) cycle
                    if( params%neg .eq. 'yes' ) call imgs(i)%neg()
                    call imgs(i)%write(params%outstk,i)
                enddo
            endif
        enddo
        call build%spproj_field%set_all('state', states) ! restore
        call build%spproj_field%write(string('reproject_oris'//TXT_EXT), [1,params%nptcls])
        if( cline%defined('state') ) where( states /= params%state ) states = 0
        ! zero state=0
        call imgs(1)%zero
        do i = 1,params%nspace
            if( states(i) == 0 ) call imgs(1)%write(params%outstk,i)
        enddo
        call update_stack_nimgs(params%outstk, params%nspace)
        ! cleanup
        call build%kill_general_tbox
        call dealloc_imgarr(imgs)
        deallocate(tmp,states)
        call simple_end('**** SIMPLE_REPROJECT NORMAL STOP ****')
    end subroutine exec_reproject

    subroutine exec_volanalyze( self, cline )
        use simple_volanalyzer
        class(commander_volanalyze), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters) :: params
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir',    'yes')
        if( .not. cline%defined('gridding') ) call cline%set('gridding', 'yes')
        call params%new(cline)
        call init_volanalyzer(params%filetab)
        call dock_compare_volumes
    end subroutine exec_volanalyze

    !> volume calculations and operations - incl Guinier, snr, mirror or b-factor
    subroutine exec_volops( self, cline )
        class(commander_volops), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        logical          :: here
        type(ori)        :: o
        type(image)      :: vol_rot
        real             :: shvec(3),  ave, sdev, maxv, minv
        call build%init_params_and_build_general_tbox(cline,params)
        here = file_exists(params%vols(1))
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

    subroutine exec_noisevol( self, cline )
        use simple_procimgstk, only: copy_imgfile
        class(commander_noisevol), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(commander_reproject) :: xreproject
        type(cmdline)             :: cline_reproject
        type(parameters)          :: params
        type(image)               :: noisevol
        integer :: s
        call params%new(cline)
        call noisevol%new([params%box,params%box,params%box], params%smpd)
        if( cline%defined('nspace') )then
            call noisevol%ran()
            call noisevol%div(real(params%box))
            call noisevol%write(string('noisevol.mrc'))
            params%refs      = 'start2Drefs.mrc'
            params%refs_even = 'start2Drefs_even.mrc'
            params%refs_odd  = 'start2Drefs_odd.mrc'
            call cline_reproject%set('prg',    'reproject')
            call cline_reproject%set('smpd',   params%smpd)
            call cline_reproject%set('pgrp',   'c1')
            call cline_reproject%set('nthr',   params%nthr)
            call cline_reproject%set('nspace', params%nspace)
            call cline_reproject%set('vol1',   'noisevol.mrc')
            call cline_reproject%set('outstk', params%refs)
            call cline_reproject%set('mskdiam',params%smpd*(real(params%box)-2.0*COSMSKHALFWIDTH))
            call xreproject%execute_safe(cline_reproject)
            call copy_imgfile(params%refs, params%refs_even, params%smpd, [1,params%nspace])
            call copy_imgfile(params%refs, params%refs_odd,  params%smpd, [1,params%nspace])
        else
            do s = 1, params%nstates
                call noisevol%ran()
                call noisevol%write(string('noisevol_state'//int2str_pad(s,2)//'.mrc'))
                call noisevol%ran()
                call noisevol%write(string('noisevol_state'//int2str_pad(s,2)//'_even.mrc'))
                call noisevol%ran()
                call noisevol%write(string('noisevol_state'//int2str_pad(s,2)//'_odd.mrc'))
            end do
        endif
        call noisevol%kill
        call cline_reproject%kill
        ! end gracefully
        call simple_end('**** SIMPLE_NOISEVOL NORMAL STOP ****')
    end subroutine exec_noisevol
    
    subroutine exec_dock_volpair( self, cline )
        class(commander_dock_volpair), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(dock_vols)  :: dvols
        type(parameters) :: params
        type(string)     :: fn_vol_docked
        if( .not. cline%defined('gridding') ) call cline%set('gridding', 'yes')
        call params%new(cline)
        fn_vol_docked = get_fbody(params%vols(2),'mrc')//'_docked.mrc'
        call dvols%new(params%vols(1), params%vols(2), params%smpd, params%hp, params%lp, params%mskdiam)
        call dvols%srch()
        call dvols%rotate_target(params%vols(2), fn_vol_docked)
        ! cleanup
        call dvols%kill()
        ! end gracefully
        call simple_end('**** SIMPLE_DOCK_VOLPAIR NORMAL STOP ****')
    end subroutine exec_dock_volpair

    !> for identification of the principal symmetry axis
    subroutine exec_symaxis_search( self, cline )
        use simple_volpft_symsrch
        class(commander_symaxis_search), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        character(len=*), parameter :: SYMSHTAB   = 'sym_3dshift'//TXT_EXT
        character(len=*), parameter :: SYMAXISTAB = 'sym_axis'//TXT_EXT
        type(parameters)      :: params
        type(builder)         :: build
        type(ori)             :: symaxis
        type(oris)            :: symaxis4write
        type(sym)             :: syme
        type(string)          :: fbody
        real                  :: shvec(3)
        integer               :: box
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('cenlp')   ) call cline%set('cenlp', 20.)
        if( .not. cline%defined('center')  ) call cline%set('center', 'yes')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'cls3D')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        call build%vol%read(params%vols(1))
        ! center volume
        shvec = 0.
        if( params%center.eq.'yes' )then
            shvec = build%vol%calc_shiftcen(params%cenlp,params%msk_crop)
            call build%vol%shift(shvec)
            fbody = get_fbody(params%vols(1),fname2ext(params%vols(1)))
            call build%vol%write(fbody//'_centered.mrc')
        endif
        ! mask volume
        call build%vol%mask(params%msk_crop, 'soft')
        ! init search object
        call volpft_symsrch_init(build%vol, params%pgrp, params%hp, params%lp)
        ! search
        call volpft_srch4symaxis(symaxis)
        call symaxis4write%new(1, is_ptcl=.false.)
        call symaxis4write%set_ori(1, symaxis)
        call symaxis4write%write(string(SYMAXISTAB), [1,1])
        if( cline%defined('projfile') )then
            if( trim(params%mkdir).eq.'yes' )then
                ! updates paths manually as project is not required in this application
                ! this is usually performed in the parameters type
                call simple_copy_file(params%projfile, filepath(PATH_HERE, basename(params%projfile)))
                params%projfile = simple_abspath(filepath(PATH_HERE, basename(params%projfile)))
                params%projname = get_fbody(params%projfile, 'simple')
                call cline%set('projname', params%projname)
                call build%spproj%update_projinfo(cline)
                call build%spproj%write_non_data_segments(params%projfile)
            endif
            call build%spproj%read(params%projfile)
            if( .not. build%spproj%is_virgin_field(params%oritype) )then
                ! transfer shift and symmetry to orientations
                call syme%new(params%pgrp)
                ! rotate the orientations & transfer the 3d shifts to 2d
                shvec = -1. * shvec ! the sign is right
                ! accounting for different volume/particle size
                box   = build%spproj%get_box()
                shvec = shvec * real(box) / real(params%box_crop)
                call syme%apply_sym_with_shift(build%spproj_field, symaxis, shvec)
                if( cline%defined('nspace') )then
                    ! making sure the projection directions assignment
                    ! refers to the input reference space
                    call build%spproj_field%set_projs(build%eulspace)
                endif
                call build%spproj%write_segment_inside(params%oritype, params%projfile)
            endif
        endif
        ! rotate volume to symaxis for later comparison
        call build%vol%read(params%vols(1))
        build%vol2 = rotvol(build%vol, symaxis)
        call build%vol2%write(string('vol_aligned2_')//trim(params%pgrp)//'axis'//params%ext%to_char())
        ! destruct
        call symaxis%kill
        call symaxis4write%kill
        call syme%kill
        ! end gracefully
        call simple_end('**** SIMPLE_SYMAXIS_SEARCH NORMAL STOP ****')
        call simple_touch('SYMAXIS_SEARCH_FINISHED')
    end subroutine exec_symaxis_search

    subroutine exec_symmetrize_map( self, cline )
        use simple_symanalyzer
        class(commander_symmetrize_map), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)   :: params
        type(builder)      :: build
        type(ori)          :: symaxis
        type(sym)          :: syme
        real               :: shvec(3), scale, smpd
        integer            :: ldim(3), box
        integer, parameter :: MAXBOX = 128
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',  'yes')
        if( .not. cline%defined('cenlp')  ) call cline%set('cenlp',    20.)
        if( .not. cline%defined('center') ) call cline%set('center', 'yes')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        if( .not. cline%defined('outvol') ) params%outvol = string('symmetrized_map')//params%ext
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
        shvec = 0.
        if( params%center.eq.'yes' )then
            shvec = build%vol%calc_shiftcen(params%cenlp,params%msk)
            call build%vol%shift(shvec)
            ! offset at scale of inputted volume
            params%xsh = shvec(1) / scale
            params%ysh = shvec(2) / scale
            params%zsh = shvec(3) / scale
        else
            params%xsh = 0.
            params%ysh = 0.
            params%zsh = 0.
        endif
        ! mask volume
        call build%vol%mask(params%msk, 'soft')
        ! symmetrize
        call symmetrize_map(build%vol, params, build%vol2, symaxis)
        if( cline%defined('projfile') )then
            if( trim(params%mkdir).eq.'yes' )then
                ! updates paths manually as project is not required in this application
                ! this is usually performed in the parameters type
                call simple_copy_file(params%projfile, filepath(PATH_HERE, basename(params%projfile)))
                params%projfile = simple_abspath(filepath(PATH_HERE, basename(params%projfile)))
                params%projname = get_fbody(params%projfile, 'simple')
                call cline%set('projname', params%projname)
                call build%spproj%update_projinfo(cline)
                call build%spproj%write_non_data_segments(params%projfile)
            endif
            call build%spproj%read(params%projfile)
            if( .not. build%spproj%is_virgin_field(params%oritype) )then
                ! transfer shift and symmetry to orientations
                call syme%new(params%pgrp)
                ! accounting for different volume/particle size
                box   = build%spproj%get_box()
                shvec = -[params%xsh,params%ysh,params%zsh] * real(box) / real(params%box_crop)
                ! rotate the orientations & transfer the 3d shifts to 2d
                call syme%apply_sym_with_shift(build%spproj_field, symaxis, shvec)
                if( cline%defined('nspace') )then
                    ! making sure the projection directions assignment
                    ! refers to the input reference space
                    call build%spproj_field%set_projs(build%eulspace)
                endif
                call build%spproj%write_segment_inside(params%oritype, params%projfile)
            endif
        endif
        ! write
        call build%vol2%write(params%outvol)
        ! end gracefully
        call simple_end('**** SIMPLE_SYMMETRIZE_MAP NORMAL STOP ****')
    end subroutine exec_symmetrize_map

    subroutine exec_symmetry_test( self, cline )
        use simple_symanalyzer
        class(commander_symmetry_test), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)      :: params
        type(builder)         :: build
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
        else
            smpd         = build%vol%get_smpd()
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

    subroutine exec_ppca_volvar( self, cline )
        use simple_imgproc,    only: make_pcavol
        use simple_pca,        only: pca
        use simple_ppca_inmem, only: ppca_inmem
        use simple_pca_svd,    only: pca_svd
        use simple_kpca_svd,   only: kpca_svd
        use simple_image,      only: image
        class(commander_ppca_volvar), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        integer,     parameter   :: MAXPCAITS = 15
        class(pca),  pointer     :: pca_ptr  => null()
        real,        allocatable :: pcavec(:,:), gen(:), avg(:)
        type(parameters)  :: params
        type(builder)     :: build
        type(image)       :: vol
        integer           :: npix
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',  'no')
        if( .not. cline%defined('outstk') ) call cline%set('outstk', 'ppca_volvar_out'//STK_EXT)
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        if( .not.file_exists(params%vols(1)) ) THROW_HARD('cannot find the inputvolume')
        call build%vol%read(params%vols(1))
        ! masking
        if(cline%defined('mskdiam')) call build%vol%mask(params%msk, 'soft', backgr=0.)
        call make_pcavol(build%vol, npix, avg, pcavec)
        ! pca allocation
        select case(trim(params%pca_mode))
            case('ppca')
                allocate(ppca_inmem :: pca_ptr)
            case('pca_svd')
                allocate(pca_svd    :: pca_ptr)
            case('kpca')
                allocate(kpca_svd   :: pca_ptr)
        end select
        call pca_ptr%new(npix, npix, params%neigs)
        call pca_ptr%master(pcavec, MAXPCAITS)
        allocate(gen(npix))
        call pca_ptr%generate(1, avg, gen)
        call vol%unserialize(gen)
        call vol%write(params%outstk)
        call vol%kill
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_PPCA_VOLVAR NORMAL STOP ****')
    end subroutine exec_ppca_volvar

end module simple_commanders_volops
