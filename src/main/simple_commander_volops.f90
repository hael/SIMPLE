! concrete commander: operations on volumes
module simple_commander_volops
include 'simple_lib.f08'
use simple_parameters,     only: parameters, params_glob
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_image,          only: image
use simple_projector_hlev, only: reproject, rotvol
use simple_ori,            only: ori
use simple_masker,         only: masker
use simple_projector,      only: projector
use simple_volprep,        only: read_and_prep_vol
use simple_volpft_srch
implicit none

public :: centervol_commander
public :: postprocess_commander
public :: automask_commander
public :: reproject_commander
public :: volops_commander
public :: volume_smat_commander
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
type, extends(commander_base) :: automask_commander
 contains
   procedure :: execute      => exec_automask
end type automask_commander
type, extends(commander_base) :: reproject_commander
 contains
   procedure :: execute      => exec_reproject
end type reproject_commander
type, extends(commander_base) :: volops_commander
  contains
    procedure :: execute      => exec_volops
end type volops_commander
type, extends(commander_base) :: volume_smat_commander
  contains
    procedure :: execute      => exec_volume_smat
end type volume_smat_commander
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
        use simple_sp_project,    only: sp_project
        use simple_estimate_ssnr, only: fsc2optlp
        class(postprocess_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters)  :: params
        type(image)       :: vol, vol_copy, vol_filt
        type(masker)      :: mskvol
        type(sp_project)  :: spproj
        character(len=:), allocatable :: vol_fname, fsc_fname, vol_filt_fname, mskfile_fname
        real,             allocatable :: fsc(:), optlp(:), res(:)
        real              :: fsc0143, fsc05, smpd, mskfile_smpd
        integer           :: state, box, fsc_box, mskfile_box, ldim(3)
        logical           :: has_fsc, has_vol_filt, has_mskfile
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
        ! check volume, gets correct smpd & box
        if( cline%defined('imgkind') )then
            call spproj%get_vol(params%imgkind, state, vol_fname, smpd, box)
        else
            call spproj%get_vol('vol', state, vol_fname, smpd, box)
        endif
        if( .not.file_exists(vol_fname) )then
            THROW_HARD('volume: '//trim(vol_fname)//' does not exist; exec_postprocess')
        endif
        ! check fsc filter
        has_fsc = .false.
        if( cline%defined('lp') )then
            ! all good
        else
            call spproj%get_fsc(state, fsc_fname, fsc_box)
            params%fsc = trim(fsc_fname)
            if( .not.file_exists(params%fsc) )then
                THROW_HARD('FSC file: '//trim(params%fsc)//' not found')
            else
                has_fsc = .true.
            endif
        endif
        ! check volume filter
        has_vol_filt = .false.
        if( cline%defined('lp') )then
            ! all good
        else
            call spproj%get_vol('vol_filt', state, vol_filt_fname, smpd, box)
            params%vol_filt = trim(vol_filt_fname)
            if( .not.file_exists(params%vol_filt) )then
                THROW_HARD('volume filter: '//trim(params%vol_filt)//' not found')
            else
                has_vol_filt = .true.
            endif
        endif
        ! check volume mask
        has_mskfile = .false.
        if( cline%defined('mskfile') )then
            if( .not.file_exists(params%mskfile) )then
                THROW_HARD('volume mask file: '//trim(params%vol_filt)//' not found')
            else
                has_mskfile = .true.
            endif
        else
            call spproj%get_vol('vol_msk', 1, mskfile_fname, mskfile_smpd, mskfile_box)
            params%mskfile = trim(mskfile_fname)
            if( file_exists(params%mskfile) ) has_mskfile = .true.
        endif
        ! read volume
        params%outvol = basename(add2fbody(trim(vol_fname), params%ext, PPROC_SUFFIX))
        ldim = [box,box,box]
        call vol%new(ldim, smpd)
        call vol%read(vol_fname)
        call vol%fft()
        if( has_fsc )then
            ! resolution & optimal low-pass filter from FSC
            fsc   = file2rarr(params%fsc)
            optlp = fsc2optlp(fsc)
            res   = vol%get_res()
            call get_resolution( fsc, res, fsc05, fsc0143 )
            where(res < TINY) optlp = 0.
        endif
        if( cline%defined('lp') )then
            ! low-pass overrides all input
            call vol%bp(0., params%lp)
        else if( has_vol_filt )then
            ! optimal low-pass filter from vol_filt
            call vol_filt%new(ldim, smpd)
            call vol_filt%read(params%vol_filt)
            call vol%apply_filter(vol_filt)
            call vol_filt%kill
        else if( has_fsc )then
            ! optimal low-pass filter from FSC
            call vol%apply_filter(optlp)
        else
            THROW_HARD('no method for low-pass filtering defined; give fsc|lp|vol_filt on command line; exec_postprocess')
        endif
        call vol_copy%copy(vol)
        ! B-factor
        if( cline%defined('bfac') ) call vol%apply_bfac(params%bfac)
        ! final low-pass filtering for smoothness
        if( has_fsc ) call vol%bp(0., fsc0143)
        ! masking
        call vol%ifft()
        if( trim(params%automsk) .eq. 'file' .and. has_mskfile )then
            call vol%zero_background
            call mskvol%new(ldim, smpd)
            call mskvol%read(params%mskfile)
            call vol%mul(mskvol)
            call mskvol%kill
        else if( params%automsk .eq. 'yes' )then
            if( .not. cline%defined('thres') )then
                write(logfhandle,*) 'Need a pixel threshold > 0. for the binarisation'
                write(logfhandle,*) 'Procedure for obtaining thresh:'
                write(logfhandle,*) '(1) postproc vol without bfac or automsk'
                write(logfhandle,*) '(2) Use UCSF Chimera to look at the volume'
                write(logfhandle,*) '(3) Identify the pixel threshold that excludes any background noise'
                THROW_HARD('postprocess')
            endif
            if( .not. cline%defined('mw') )then
                THROW_HARD('molecular weight must be provided for auto-masking (MW); postprocess')
            endif
            call vol_copy%ifft
            call mskvol%automask3D(vol_copy)
            call mskvol%write('automask'//params%ext)
            call vol%zero_background
            call vol%mul(mskvol)
            call mskvol%kill
        else
            if( has_mskfile )then
                call mskvol%new(ldim, smpd)
                call mskvol%read(params%mskfile)
                call vol%mul(mskvol)
                call mskvol%kill
            else
                if( params%l_innermsk )then
                    call vol%mask(params%msk, 'soft', inner=params%inner, width=params%width)
                else
                    call vol%mask(params%msk, 'soft')
                endif
            endif
        endif
        ! output in cwd
        call vol%write(params%outvol)
        ! also output mirrored by default (unless otherwise stated on command line)
        if( .not. cline%defined('mirr') .or. params%mirr .ne. 'no' )then
            call vol%mirror('x')
            params%outvol = add2fbody(params%outvol, params%ext, MIRR_SUFFIX)
            call vol%write(params%outvol)
        endif
        ! destruct
        call spproj%kill
        call vol%kill
        call vol_copy%kill
        call simple_end('**** SIMPLE_POSTPROCESS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_postprocess

    subroutine exec_automask( self, cline )
        class(automask_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(builder)    :: build
        type(parameters) :: params
        type(masker)     :: mskvol
        call params%new(cline)
        call build%build_general_tbox(params, cline)
        call build%vol%read(params%vols(1))
        call mskvol%automask3D(build%vol)
        call mskvol%write('automask'//params%ext)
        call simple_end('**** SIMPLE_AUTOMASK NORMAL STOP ****')
    end subroutine exec_automask

    !> exec_project generate projections from volume
    subroutine exec_reproject( self, cline )
        use simple_binoris_io, only: binread_nlines
        class(reproject_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters)         :: params
        type(builder)            :: build
        type(image), allocatable :: imgs(:)
        integer                  :: i
        logical                  :: do_zero
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('wfun')   ) call cline%set('wfun',   'kb')
        if( .not. cline%defined('winsz')  ) call cline%set('winsz',   1.5)
        if( .not. cline%defined('alpha')  ) call cline%set('alpha',    2.)
        if( .not. cline%defined('outstk') ) call cline%set('outstk', 'reprojs.mrcs')
        if( .not. cline%defined('oritab') )then
            if( .not. cline%defined('nspace') ) THROW_HARD('need nspace (for number of projections)!')
        endif
        call params%new(cline)
        if( cline%defined('oritab') )then
            params%nptcls = binread_nlines(params%oritab)
            call build%build_general_tbox(params, cline)
            params%nspace = build%spproj_field%get_noris()
        else
            params%nptcls = params%nspace
            call build%build_general_tbox(params, cline)
            call build%pgrpsyms%build_refspiral(build%spproj_field)
        endif
        ! fix volumes and stacks
        call build%vol%read(params%vols(1))
        ! masking
        if(cline%defined('msk')) call build%vol%mask(params%msk, 'soft')
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
        call build%init_params_and_build_general_tbox(cline,params,boxmatch_off=.true.)
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
                if( .not. cline%defined('msk')  ) THROW_HARD('need msk (mask radius) input for volume rotation')
                call o%new(is_ptcl=.false.)
                call o%set_euler([params%e1,params%e2,params%e3])
                shvec     = [params%xsh,params%ysh,params%zsh]
                vol_rot   = rotvol(build%vol, o,  shvec)
                build%vol = vol_rot
            endif
            if( cline%defined('xsh') .or. cline%defined('ysh') .or. cline%defined('zsh') )then
                call build%vol%shift([params%xsh,params%ysh,params%zsh])
            endif
            if( cline%defined('lp_backgr') .and. cline%defined('mskfile') )then
                call build%vol2%new(build%vol%get_ldim(), build%vol%get_smpd())
                call build%vol2%read(params%mskfile)
                call build%vol%lp_background(build%vol2, params%lp_backgr)
            endif
            call build%vol%write(params%outvol, del_if_exists=.true.)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_VOLOPS NORMAL STOP ****')
    end subroutine exec_volops

    !> calculate similarity matrix between volumes
    subroutine exec_volume_smat( self, cline )
        class(volume_smat_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters) :: params
        integer          :: funit, io_stat, cnt, npairs,  nvols, loc(1)
        integer          :: ivol, jvol, ldim(3), ipair, ifoo, i, spat_med
        integer          :: furthest_from_spat_med
        real             :: corr_max, corr_min, spat_med_corr
        real             :: furthest_from_spat_med_corr
        type(projector)  :: vol1, vol2
        type(ori)        :: o
        real,                      allocatable :: corrmat(:,:), corrs(:), corrs_avg(:)
        integer,                   allocatable :: pairs(:,:)
        character(len=LONGSTRLEN), allocatable :: vollist(:)
        character(len=:),          allocatable :: fname
        call params%new(cline)
        call read_filetable(params%vollist, vollist) ! reads in list of volumes
        nvols  = size(vollist)
        npairs = (nvols*(nvols-1))/2
        ! find logical dimension & make volumes for matching
        call find_ldim_nptcls(vollist(1), ldim, ifoo)
        if( cline%defined('part') )then
            npairs = params%top-params%fromp+1
            allocate(corrs(params%fromp:params%top), pairs(params%fromp:params%top,2), stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('In: simple_volume_smat, 1',alloc_stat)
            ! read the pairs
            allocate(fname, source='pairs_part'//int2str_pad(params%part,params%numlen)//'.bin')
            if( .not. file_exists(fname) ) THROW_HARD('I/O; simple_volume_smat')
            call fopen(funit, status='OLD', action='READ', file=fname, access='STREAM', iostat=io_stat)
            if(io_stat/=0) call fileiochk('volops; vol_smat opening ', io_stat)
            read(unit=funit,pos=1,iostat=io_stat) pairs(params%fromp:params%top,:)
            ! Check if the read was successful
            if(io_stat/=0) then
                call fileiochk('**ERROR(simple_volume_smat): I/O error reading file: '//trim(fname), io_stat)
            endif
            call fclose(funit, errmsg='volops; vol_smat closing '//trim(fname))
            deallocate(fname)
            cnt = 0
            do ipair=params%fromp,params%top
                cnt = cnt + 1
                call progress(cnt, npairs)
                ivol = pairs(ipair,1)
                jvol = pairs(ipair,2)
                call read_and_prep_vol( vollist(ivol), vol1 )
                call read_and_prep_vol( vollist(jvol), vol2 )
                call volpft_srch_init(vol1, vol2, params%hp, params%lp)
                o = volpft_srch_minimize()
                corrs(ipair) = o%get('corr')
            end do
            ! write the similarities
            allocate(fname, source='similarities_part'//int2str_pad(params%part,params%numlen)//'.bin',stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk( 'volops; volume smat ',alloc_stat)
            call fopen(funit, status='REPLACE', action='WRITE', &
                 file=fname, access='STREAM', iostat=io_stat)
            if(io_stat/=0) call fileiochk('volops; volume smat 2  opening ', io_stat)
            write(unit=funit,pos=1,iostat=io_stat) corrs(params%fromp:params%top)
            ! Check if the write was successful
            if(io_stat/=0) &
                call fileiochk('**ERROR(simple_volume_smat): I/O error writing file: '//trim(fname), io_stat)
            call fclose(funit, errmsg='volops; volume smat 2  closing ')
            deallocate(fname, corrs, pairs)
        else
            ! generate similarity matrix
            allocate(corrmat(nvols,nvols), corrs_avg(nvols), stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('In: simple_volume_smat, 2',alloc_stat)
            corrmat = -1.
            forall(i=1:nvols) corrmat(i,i) = 1.0
            cnt = 0
            corr_max = -1.0
            corr_min =  1.0
            do ivol=1,nvols - 1
                do jvol=ivol + 1,nvols
                    cnt = cnt + 1
                    call progress(cnt, npairs)
                    call read_and_prep_vol( vollist(ivol), vol1 )
                    call read_and_prep_vol( vollist(jvol), vol2 )
                    call volpft_srch_init(vol1, vol2, params%hp, params%lp)
                    o = volpft_srch_minimize()
                    corrmat(ivol,jvol) = o%get('corr')
                    corrmat(jvol,ivol) = corrmat(ivol,jvol)
                    if( corrmat(ivol,jvol) > corr_max ) corr_max = corrmat(ivol,jvol)
                    if( corrmat(ivol,jvol) < corr_min ) corr_min = corrmat(ivol,jvol)
                end do
            end do
            do ivol=1,nvols
                corrs_avg(ivol) = (sum(corrmat(ivol,:))-1.0)/real(nvols - 1)
            end do
            loc                         = maxloc(corrs_avg)
            spat_med                    = loc(1)
            spat_med_corr               = corrs_avg(spat_med)
            loc                         = minloc(corrmat(spat_med,:))
            furthest_from_spat_med      = loc(1)
            furthest_from_spat_med_corr = corrmat(spat_med,furthest_from_spat_med)
            write(logfhandle,'(a,1x,f7.4)') 'MAX VOL PAIR CORR          :', corr_max
            write(logfhandle,'(a,1x,f7.4)') 'MIN VOL PAIR CORR          :', corr_min
            write(logfhandle,'(a,1x,i7)'  ) 'SPATIAL MEDIAN             :', spat_med
            write(logfhandle,'(a,1x,f7.4)') 'SPATIAL MEDIAN CORR        :', spat_med_corr
            write(logfhandle,'(a,1x,i7)'  ) 'FURTHEST FROM SPAT MED     :', furthest_from_spat_med
            write(logfhandle,'(a,1x,f7.4)') 'FURTHEST FROM SPAT MED CORR:', furthest_from_spat_med_corr
            call fopen(funit, status='REPLACE', action='WRITE', file='vol_smat.bin', &
                &access='STREAM', iostat=io_stat)
             if(io_stat/=0)call fileiochk('volops; volume smat 3 opening ', io_stat)
            write(unit=funit,pos=1,iostat=io_stat) corrmat
            if(io_stat/=0) &
                call fileiochk('**ERROR(simple_volume_smat): I/O error writing to vol_smat.bin', io_stat)

            call fclose(funit, errmsg='volops; volume smat 3 closing ')
            if(allocated(corrmat))deallocate(corrmat)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_VOLUME_SMAT NORMAL STOP ****')
    end subroutine exec_volume_smat

    subroutine exec_dock_volpair( self, cline )
        use simple_vol_srch
        use simple_oris, only: oris
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
                ! do i=1,10
                !     ! rotate and shift vol to create reference for shift alignment
                !     call vol2%ifft
                !     vol_out = rotvol(vol2, orientation, cxyz(2:4))
                !     call vol2%fft
                !     call vol_out%fft
                !     ! obtain joint shifts by vector addition
                !     cxyz(2:4) = cxyz(2:4) + cxyz2(2:4)
                !     ! update shifts in volpft_srch class and refine angles
                !     call volpft_srch_set_shvec(cxyz(2:4))
                !     if( i <= 3 )then
                !         orientation_best = volpft_srch_refine(orientation)
                !         cxyz2 = vol_shsrch_minimize()
                !         if( cxyz2(1) > corr_prev )then ! a better solution was found
                !             corr_prev = cxyz2(1)
                !         endif
                !     else
                !         orientation_best = volpft_srch_refine(orientation, angres=5.)
                !         if( cxyz2(1) > corr_prev )then ! a better solution was found
                !             corr_prev = cxyz2(1)
                !         endif
                !     endif
                !     orientation = orientation_best
                !     call orientation%set('x', cxyz(2))
                !     call orientation%set('y', cxyz(3))
                !     call orientation%set('z', cxyz(4))
                ! end do
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
        use simple_sym,  only: sym
        use simple_ori,  only: ori
        use simple_oris, only: oris
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
        if( params_glob%l_innermsk )then
            call build%vol%mask(params%msk, 'soft', inner=params%inner, width=params%width)
        else
            call build%vol%mask(params%msk, 'soft')
        endif
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
            params%inner = round2even(scale * params%inner)
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
            params%inner = round2even(scale * params%inner)
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
        if( params_glob%l_innermsk )then
            call build%vol%mask(params%msk, 'soft', inner=params%inner, width=params%width)
        else
            call build%vol%mask(params%msk, 'soft')
        endif
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
