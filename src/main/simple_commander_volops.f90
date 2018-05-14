! concrete commander: operations on volumes
module simple_commander_volops
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_image,          only: image
use simple_projector_hlev, only: reproject, rotvol
use simple_ori,            only: ori
use simple_masker,         only: masker
use simple_projector,      only: projector
use simple_volprep,        only: read_and_prep_vol
use simple_volpft_srch
implicit none

public :: fsc_commander
public :: centervol_commander
public :: postprocess_commander
public :: reproject_commander
public :: volaverager_commander
public :: volops_commander
public :: volume_smat_commander
public :: dock_volpair_commander
public :: make_pickrefs_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: fsc_commander
  contains
    procedure :: execute      => exec_fsc
end type fsc_commander
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
type, extends(commander_base) :: dock_volpair_commander
  contains
    procedure :: execute      => exec_dock_volpair
end type dock_volpair_commander
type, extends(commander_base) :: make_pickrefs_commander
contains
    procedure :: execute      => exec_make_pickrefs
end type make_pickrefs_commander
contains

    !> calculates Fourier shell correlation from Even/Odd Volume pairs
    subroutine exec_fsc( self, cline )
        class(fsc_commander), intent(inout) :: self
        class(cmdline),       intent(inout) :: cline
        type(params)      :: p
        type(image)       :: even, odd
        type(masker)      :: mskvol
        integer           :: j, find_plate
        real              :: res_fsc05, res_fsc0143
        real, allocatable :: res(:), corrs(:)
        p = params(cline)
        ! read even/odd pair
        call even%new([p%box,p%box,p%box], p%smpd)
        call odd%new([p%box,p%box,p%box], p%smpd)
        call odd%read(p%vols(1))
        call even%read(p%vols(2))
        ! always normalise before masking
        call even%norm()
        call odd%norm()
        if( cline%defined('mskfile') )then
            if( file_exists(p%mskfile) )then
                call mskvol%new([p%box,p%box,p%box], p%smpd)
                call mskvol%read(p%mskfile)
                call mskvol%resmask(p)
                call mskvol%write('resmask'//p%ext)
                call even%zero_background
                call odd%zero_background
                call even%mul(mskvol)
                call odd%mul(mskvol)
            else
                write(*,*) 'the inputted mskfile: ', trim(p%mskfile)
                stop 'does not exists in cwd; commander_volops :: exec_fsc'
            endif
        else
            ! spherical masking
            if( p%l_innermsk )then
                call even%mask(p%msk, 'soft', inner=p%inner, width=p%width)
                call odd%mask(p%msk, 'soft', inner=p%inner, width=p%width)
            else
                call even%mask(p%msk, 'soft')
                call odd%mask(p%msk, 'soft')
            endif
        endif
        ! forward FT
        call even%fft()
        call odd%fft()
        ! calculate FSC
        res = even%get_res()
        allocate(corrs(even%get_filtsz()))
        call even%fsc(odd, corrs)
        if( p%l_phaseplate ) call phaseplate_correct_fsc(corrs, find_plate)
        do j=1,size(res)
           write(*,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
        end do
        call get_resolution(corrs, res, res_fsc05, res_fsc0143)
        write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
        write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
        call even%kill
        call odd%kill
        ! end gracefully
        call simple_end('**** SIMPLE_FSC NORMAL STOP ****')
    end subroutine exec_fsc

    !> centers a 3D volume and associated particle document
    subroutine exec_centervol( self, cline )
        class(centervol_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params)         :: p
        type(build)          :: b
        real, allocatable    :: shvec(:,:)
        integer              :: istate
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        ! center volume(s)
        allocate(shvec(p%nstates,3))
        do istate=1,p%nstates
            call b%vol%read(p%vols(istate))
            shvec(istate,:) = b%vol%center(p%cenlp, p%msk)
            call b%vol%shift([shvec(istate,1),shvec(istate,2),shvec(istate,3)])
            call b%vol%write('shifted_vol_state'//int2str(istate)//p%ext)
            ! transfer the 3D shifts to 2D
            if( cline%defined('oritab') ) call b%a%map3dshift22d(-shvec(istate,:), state=istate)
        end do
        if( cline%defined('oritab') ) call b%a%write(p%outfile, [1,b%a%get_noris()])
        ! end gracefully
        call simple_end('**** SIMPLE_CENTER NORMAL STOP ****')
    end subroutine exec_centervol

    subroutine exec_postprocess(self, cline)
        use simple_estimate_ssnr, only: fsc2optlp
        class(postprocess_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)      :: p
        type(build)       :: b
        type(image)       :: vol_copy, vol_filt
        type(masker)      :: mskvol
        real, allocatable :: fsc(:), optlp(:), res(:)
        real              :: fsc0143, fsc05
        integer           :: state, ldim(3)
        ! pre-proc
        p = params(cline) ! constants & derived constants produced, mode=2
        call b%build_general_tbox(p, cline) ! general objects built
        if( cline%defined('projfile') )then

            ! 2do

        else if( cline%defined('vol1') )then
            state = 1
            call b%vol%read(p%vols(state))
            call b%vol%fft()
        else
            stop 'ERROR! either projfile or vol1 needs to be part of command line; commander_volops :: postprocess'
        endif

        if( cline%defined('fsc') )then
            ! optimal low-pass filter from FSC
            if( file_exists(p%fsc) )then
                fsc   = file2rarr(p%fsc)
                optlp = fsc2optlp(fsc)
            else
                write(*,*) 'FSC file: ', trim(p%fsc), ' not in cwd'
                stop
            endif
            res = b%vol%get_res()
            call get_resolution( fsc, res, fsc05, fsc0143 )
            where(res < TINY) optlp = 0.
        endif
        if( cline%defined('vol_filt') )then
            ! optimal low-pass filter from input vol_filt
            call vol_filt%new(b%vol%get_ldim(), p%smpd)
            call vol_filt%read(p%vol_filt)
            call b%vol%apply_filter(vol_filt)
            call vol_filt%kill
        else if( cline%defined('fsc') )then
            ! optimal low-pass filter from FSC
            call b%vol%apply_filter(optlp)
        else if( cline%defined('lp') )then
            ! ad hoc low-pass filter
            call b%vol%bp(0., p%lp)
        else
            write(*,*) 'no method for low-pass filtering defined; give fsc|lp|vol_filt on command line'
            stop 'simple_commander_volops :: exec_postprocess'
        endif
        call vol_copy%copy(b%vol)
        ! B-fact
        if( cline%defined('bfac') ) call b%vol%apply_bfac(p%bfac)
        ! final low-pass filtering for smoothness
        if( cline%defined('fsc')  ) call b%vol%bp(0., fsc0143)
        ! masking
        call b%vol%ifft()
        if( cline%defined('mskfile') )then
            if( file_exists(p%mskfile) )then
                ldim = b%vol%get_ldim()
                call b%vol%zero_background
                call mskvol%new(ldim, p%smpd)
                call mskvol%read(p%mskfile)
                call b%vol%mul(mskvol)
            else
                write(*,*) 'file: ', trim(p%mskfile)
                stop 'maskfile does not exists in cwd'
            endif
        else if( p%automsk .eq. 'yes' )then
            if( .not. cline%defined('thres') )then
                write(*,*) 'Need a pixel threshold > 0. for the binarisation'
                write(*,*) 'Procedure for obtaining thresh:'
                write(*,*) '(1) postproc vol without bfac or automsk'
                write(*,*) '(2) Use UCSF Chimera to look at the volume'
                write(*,*) '(3) Identify the pixel threshold that excludes any background noise'
                stop 'commander_volops :: postprocess'
            endif
            if( .not. cline%defined('mw') )then
                write(*,*) 'Molecular weight must be provided for auto-masking (MW)'
                stop 'commander_volops :: postprocess'
            endif
            call vol_copy%ifft()
            call mskvol%automask3D(p, vol_copy)
            call mskvol%write('automask'//p%ext)
            call b%vol%zero_background
            call b%vol%mul(mskvol)
        else
            if( p%l_innermsk )then
                call b%vol%mask(p%msk, 'soft', inner=p%inner, width=p%width)
            else
                call b%vol%mask(p%msk, 'soft')
            endif
        endif
        ! output
        p%outvol = add2fbody(trim(p%vols(state)), p%ext, PPROC_SUFFIX)
        call b%vol%write(p%outvol)
        ! also output mirrored by default (unless otherwise stated on command line)
        if( .not. cline%defined('mirr') .or. p%mirr .ne. 'no' )then
            call b%vol%mirror('x')
            p%outvol = add2fbody(p%outvol, p%ext, '_mirr')
            call b%vol%write(p%outvol)
        endif
        ! destruct
        call vol_copy%kill
        call mskvol%kill
        call simple_end('**** SIMPLE_POSTPROCESS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_postprocess

    !> exec_project generate projections from volume
    !! \param cline command line
    !!
    subroutine exec_reproject( self, cline )
        use simple_binoris_io, only: binread_nlines
        class(reproject_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)             :: p
        type(build)              :: b
        type(image), allocatable :: imgs(:)
        integer                  :: i, loop_end
        real                     :: x, y
        if( .not. cline%defined('oritab') )then
            if( .not. cline%defined('nspace') ) stop 'need nspace (for number of projections)!'
        endif
        p = params(cline) ! parameters generated
        if( cline%defined('oritab') )then
            p%nptcls = binread_nlines(p, p%oritab)
            call b%build_general_tbox(p, cline)
            p%nspace = b%a%get_noris()
        else
            p%nptcls = p%nspace
            call b%build_general_tbox(p, cline)
            call b%a%spiral(p%nsym, p%eullims)
        endif
        ! fix volumes and stacks
        call b%vol%read(p%vols(1))
        DebugPrint 'read volume'
        ! masking
        if(cline%defined('msk')) call b%vol%mask(p%msk, 'soft')
        ! generate projections
        imgs = reproject(b%vol, b%a, p)
        loop_end = p%nspace
        if( file_exists(p%outstk) ) call del_file(p%outstk)
        do i=1,loop_end
            if( p%neg .eq. 'yes' ) call imgs(i)%neg()
            call imgs(i)%write(p%outstk,i)
        end do
        call b%a%write('project_oris'//trim(TXT_EXT), [1,p%nptcls])
        call simple_end('**** SIMPLE_PROJECT NORMAL STOP ****')
    end subroutine exec_reproject

    !> exec_volaverager create volume average
    !! \param cline
    !!
    subroutine exec_volaverager( self, cline )
        class(volaverager_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer, allocatable               :: ptcls(:)
        character(len=STDLEN), allocatable :: volnames(:)
        type(image)                        :: vol_avg
        integer                            :: istate, ivol, nvols, funit_vols, numlen, ifoo, io_stat
        character(len=:), allocatable      :: fname
        character(len=1)                   :: fformat
        p = params(cline) ! parameters generated
        ! read the volnames
        nvols = nlines(p%vollist)
        DebugPrint   'number of volumes: ', nvols
        allocate(volnames(nvols))
        call fopen(funit_vols, status='old', file=p%vollist, iostat=io_stat)
        if(io_stat.ne.0) call fileiochk('volops; volaverager opening ', io_stat)
        write(*,'(a)') '>>> MAKING PATTERN STACK'
        do ivol=1,nvols
            read(funit_vols,'(a256)') volnames(ivol)
            DebugPrint   'read volname: ', volnames(ivol)
        end do
        call fclose(funit_vols,errmsg='volops; volaverager closing ')
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
        DebugPrint   'file extension: ', p%ext
        ! average the states
        call vol_avg%copy(b%vol)
        p%nstates = b%a%get_n('state')
        DebugPrint   'number of states: ', p%nstates
        numlen = len(int2str(p%nstates))
        do istate=1,p%nstates
            DebugPrint   'processing state: ', istate
            call b%a%get_pinds(istate, 'state', ptcls)
            vol_avg = 0.
            do ivol=1,size(ptcls)
                call b%vol%read(volnames(ptcls(ivol)))
                DebugPrint   'read volume: ', volnames(ptcls(ivol))
                call vol_avg%add(b%vol)
            end do
            call vol_avg%div(real(size(ptcls)))
            allocate(fname, source='sumvol_state'//int2str_pad(istate, numlen)//p%ext)
            DebugPrint   'trying to write volume to file: ', fname
            call vol_avg%write(fname)
            deallocate(ptcls,fname)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_VOLAVERAGER NORMAL STOP ****')
    end subroutine exec_volaverager

    !> volume calculations and operations - incl Guinier, snr, mirror or b-factor
    !! \param cline commandline
    !!
    subroutine exec_volops( self, cline )
        class(volops_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        logical           :: here
        type(ori)         :: o
        type(image)       :: vol_rot
        real              :: shvec(3)
        p = params(cline) ! constants & derived constants produced, mode=2
        call b%build_general_tbox(p, cline)  ! general objects built
        ! reallocate vol (boxmatch issue)
        call b%vol%new([p%box,p%box,p%box], p%smpd)
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
            if( cline%defined('neg')  )then
                call b%vol%neg()
            end if
            if( cline%defined('snr') )then
                call b%vol%add_gauran(p%snr)
            end if
            if( cline%defined('mirr') )then
                call b%vol%mirror(p%mirr)
            end if
            if( cline%defined('bfac') )then
                call b%vol%apply_bfac(p%bfac)
            end if
            if( cline%defined('e1') .or. cline%defined('e2') .or. cline%defined('e3') )then
                if( .not. cline%defined('smpd') ) stop 'need smpd (sampling distance) input for volume rotation'
                if( .not. cline%defined('msk')  ) stop 'need msk (mask radius) input for volume rotation'
                call o%new
                call o%set_euler([p%e1,p%e2,p%e3])
                shvec   = [p%xsh,p%ysh,p%zsh]
                vol_rot = rotvol(b%vol, o, p, shvec)
                b%vol   = vol_rot
            endif
            if( cline%defined('xsh') .or. cline%defined('ysh') .or. cline%defined('zsh') )then
                call b%vol%shift([p%xsh,p%ysh,p%zsh])
            endif
            call b%vol%write(p%outvol, del_if_exists=.true.)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_VOLOPS NORMAL STOP ****')
    end subroutine exec_volops

    !> calculate similarity matrix between volumes
    subroutine exec_volume_smat( self, cline )
        class(volume_smat_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params), target :: p
        integer              :: funit, io_stat, cnt, npairs,  nvols, loc(1)
        integer              :: ivol, jvol, ldim(3), ipair, ifoo, i, spat_med
        integer              :: furthest_from_spat_med
        real                 :: corr_max, corr_min, spat_med_corr
        real                 :: furthest_from_spat_med_corr
        type(projector)      :: vol1, vol2
        type(ori)            :: o
        real,                  allocatable :: corrmat(:,:), corrs(:), corrs_avg(:)
        integer,               allocatable :: pairs(:,:)
        character(len=STDLEN), allocatable :: vollist(:)
        character(len=:),      allocatable :: fname
        p = params(cline, .false.)              ! constants & derived constants produced
        call read_filetable(p%vollist, vollist) ! reads in list of volumes
        nvols  = size(vollist)
        npairs = (nvols*(nvols-1))/2
        ! find logical dimension & make volumes for matching
        call find_ldim_nptcls(vollist(1), ldim, ifoo)
        DebugPrint   'found logical dimension: ', ldim
        if( cline%defined('part') )then
            npairs = p%top-p%fromp+1
            DebugPrint 'allocating this number of similarities: ', npairs
            allocate(corrs(p%fromp:p%top), pairs(p%fromp:p%top,2), stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('In: simple_volume_smat, 1',alloc_stat)
            ! read the pairs
            allocate(fname, source='pairs_part'//int2str_pad(p%part,p%numlen)//'.bin')
            if( .not. file_exists(fname) ) stop 'I/O error; simple_volume_smat'
            call fopen(funit, status='OLD', action='READ', file=fname, access='STREAM', iostat=io_stat)
            if(io_stat/=0) call fileiochk('volops; vol_smat opening ', io_stat)
            DebugPrint   'reading pairs in range: ', p%fromp, p%top
            read(unit=funit,pos=1,iostat=io_stat) pairs(p%fromp:p%top,:)
            ! Check if the read was successful
            if(io_stat/=0) then
                call fileiochk('**ERROR(simple_volume_smat): I/O error reading file: '//trim(fname), io_stat)
            endif
            call fclose(funit, errmsg='volops; vol_smat closing '//trim(fname))
            deallocate(fname)
            cnt = 0
            do ipair=p%fromp,p%top
                cnt = cnt + 1
                call progress(cnt, npairs)
                ivol = pairs(ipair,1)
                jvol = pairs(ipair,2)
                call read_and_prep_vol( p, vollist(ivol), vol1 )
                call read_and_prep_vol( p, vollist(jvol), vol2 )
                call volpft_srch_init(vol1, vol2, p%hp, p%lp, 0.)
                o = volpft_srch_minimize_eul()
                corrs(ipair) = o%get('corr')
            end do
            DebugPrint   'did set this number of similarities: ', cnt
            ! write the similarities
            allocate(fname, source='similarities_part'//int2str_pad(p%part,p%numlen)//'.bin',stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk( 'volops; volume smat ',alloc_stat)
            call fopen(funit, status='REPLACE', action='WRITE', &
                 file=fname, access='STREAM', iostat=io_stat)
            if(io_stat/=0) call fileiochk('volops; volume smat 2  opening ', io_stat)
            write(unit=funit,pos=1,iostat=io_stat) corrs(p%fromp:p%top)
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
                    call read_and_prep_vol( p, vollist(ivol), vol1 )
                    call read_and_prep_vol( p, vollist(jvol), vol2 )
                    call volpft_srch_init(vol1, vol2, p%hp, p%lp, 0.)
                    o = volpft_srch_minimize_eul()
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
            write(*,'(a,1x,f7.4)') 'MAX VOL PAIR CORR          :', corr_max
            write(*,'(a,1x,f7.4)') 'MIN VOL PAIR CORR          :', corr_min
            write(*,'(a,1x,i7)'  ) 'SPATIAL MEDIAN             :', spat_med
            write(*,'(a,1x,f7.4)') 'SPATIAL MEDIAN CORR        :', spat_med_corr
            write(*,'(a,1x,i7)'  ) 'FURTHEST FROM SPAT MED     :', furthest_from_spat_med
            write(*,'(a,1x,f7.4)') 'FURTHEST FROM SPAT MED CORR:', furthest_from_spat_med_corr
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
        class(dock_volpair_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)    :: p
        type(projector) :: vol1, vol2
        type(image)     :: vol_out
        type(ori)       :: orientation
        p = params(cline, .false.) ! constants & derived constants produced
        call read_and_prep_vol( p, p%vols(1), vol1 )
        call read_and_prep_vol( p, p%vols(2), vol2 )
        call volpft_srch_init(vol1, vol2, p%hp, p%lp, 0.)
        select case(p%dockmode)
            case('eul')
                orientation = volpft_srch_minimize_eul()
                vol_out     = rotvol(vol2, orientation, p)
            case('shift')
                orientation = volpft_srch_minimize_shift()
                call vol2%shift(orientation%get_3Dshift())
                vol_out     = vol2
            case('eulshift')
                orientation = volpft_srch_minimize_eul()
                orientation = volpft_srch_minimize_shift()
                vol_out     = rotvol(vol2, orientation, p, orientation%get_3Dshift())
            case('all')
                orientation = volpft_srch_minimize_eul()
                orientation = volpft_srch_minimize_shift()
                orientation = volpft_srch_minimize_all()
                vol_out     = rotvol(vol2, orientation, p, orientation%get_3Dshift())
        end select
        call vol_out%write(p%outvol, del_if_exists=.true.)
        ! end gracefully
        call simple_end('**** SIMPLE_DOCK_VOLPAIR NORMAL STOP ****')
    end subroutine exec_dock_volpair


    !> for making picker references
    subroutine exec_make_pickrefs( self, cline )
        use simple_procimgfile, only: neg_imgfile
        class(make_pickrefs_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline !< command line input
        integer,          parameter :: NREFS=100, NPROJS=20
        character(len=*), parameter :: ORIFILE='pickrefs_oris'//trim(TXT_EXT)
        type(params)                :: p
        type(build)                 :: b
        type(cmdline)               :: cline_project
        type(reproject_commander)   :: xproject
        integer                     :: nrots, cnt, iref, irot, status
        real                        :: ang, rot
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        if( cline%defined('stk') .or. cline%defined('vol1') )then
            p = params(cline) ! constants & derived constants produced
            if( cline%defined('vol1') )then
                p%nptcls = NPROJS
                call b%a%new(NPROJS)
                call b%a%spiral( p%nsym, p%eullims )
                call b%a%write(trim(ORIFILE), [1,NPROJS])
                cline_project = cline
                call cline_project%set('nspace', real(NPROJS))
                p%stk = 'even_projs'//p%ext
                call cline_project%set('outstk', trim(p%stk)  )
                call cline_project%set('oritab', trim(ORIFILE))
                call cline_project%set('smpd',   PICKER_SHRINK)
                call cline_project%set('neg',    'no'         )
                call cline_project%set('msk',    real(p%box/2-5))
                call xproject%execute(cline_project)
            endif
            ! expand in in-plane rotation
            nrots = NREFS/p%nptcls
            if( nrots > 1 )then
                ang = 360./real(nrots)
                rot = 0.
                cnt  = 0
                do iref=1,p%nptcls
                    call b%img%read(p%stk, iref)
                    do irot=1,nrots
                        cnt = cnt + 1
                        call b%img%rtsq(rot, 0., 0., b%img_copy)
                        call b%img_copy%write('rotated_from_make_pickrefs'//p%ext, cnt)
                        rot = rot + ang
                    end do
                end do
                call cline%set('stk', 'rotated_from_make_pickrefs'//p%ext)
            endif
            if( p%pcontrast .eq. 'black' )then
                call neg_imgfile('rotated_from_make_pickrefs'//p%ext, 'pickrefs'//p%ext, p%smpd)
            else
                status= simple_rename('rotated_from_make_pickrefs'//p%ext, 'pickrefs'//p%ext)
            endif
        else
            stop 'need input volume (vol1) or class averages (stk) to generate picking references'
        endif
        call del_file('rotated_from_make_pickrefs'//p%ext)
        ! end gracefully
        call simple_end('**** SIMPLE_MAKE_PICKREFS NORMAL STOP ****')
    end subroutine exec_make_pickrefs


end module simple_commander_volops
