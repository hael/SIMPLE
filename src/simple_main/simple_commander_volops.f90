!==Class simple_commander_volops
!
! This class contains the set of concrete volops (volume operations) commanders of the SIMPLE library. This class provides the glue between the reciver 
! (main reciever is simple_exec program) and the abstract action, which is simply execute (defined by the base class: simple_commander_base). 
! Later we can use the composite pattern to create MacroCommanders (or workflows)
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Cyril Reboul & Hans Elmlund 2016
!
module simple_commander_volops
use simple_defs
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_strings,        only: int2str, int2str_pad
use simple_jiffys          ! use all in there
use simple_filehandling    ! use all in there
implicit none

public :: fsc_commander
public :: cenvol_commander
public :: postproc_vol_commander
public :: projvol_commander
public :: volaverager_commander
public :: volops_commander
public :: volume_smat_commander
private


type, extends(commander_base) :: fsc_commander
  contains
    procedure :: execute      => exec_fsc
end type fsc_commander
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

contains

    subroutine exec_fsc( self, cline )
        use simple_image, only: image
        use simple_math,  only: get_resolution
        class(fsc_commander), intent(inout) :: self
        class(cmdline),       intent(inout) :: cline
        type(params)      :: p
        type(image)       :: even, odd
        integer           :: j
        real              :: res_fsc05, res_fsc0143
        real, allocatable :: res(:), corrs(:)
        p = params(cline)
        ! read even/odd pair
        call even%new([p%box,p%box,p%box], p%smpd)
        call odd%new([p%box,p%box,p%box], p%smpd)
        call odd%read(p%vols(1))
        call even%read(p%vols(2))
        ! forward FT
        call even%fwd_ft
        call odd%fwd_ft
        ! calculate FSC
        call even%fsc(odd, res, corrs)
        do j=1,size(res)
           write(*,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
        end do
        call get_resolution(corrs, res, res_fsc05, res_fsc0143)
        write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
        write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
        call even%kill
        call odd%kill
    end subroutine exec_fsc
    
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
            shvec(istate,:) = b%vol%center(p%cenlp,'no',p%msk)
            if( istate == 1 ) call b%vol%write(p%outvol)
            if( debug )then
                call b%vol%shift(-shvec(istate,1), -shvec(istate,2), -shvec(istate,3))
                call b%vol%write('shifted_vol_state'//int2str(istate)//p%ext)
            endif
            ! transfer the 3D shifts to 2D
            if( cline%defined('oritab') ) call b%a%map3dshift22d(-shvec(istate,:), state=istate)
        end do
        if( cline%defined('oritab') ) call b%a%write(p%outfile)
        ! end gracefully
        call simple_end('**** SIMPLE_CENVOL NORMAL STOP ****')
    end subroutine exec_cenvol
    
    subroutine exec_postproc_vol(self, cline)
        use simple_math, only: get_resolution
        use simple_estimate_ssnr ! use all in there
        class(postproc_vol_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)      :: p
        type(build)       :: b
        real, allocatable :: fsc(:), optlp(:), res(:)
        real              :: fsc0143, fsc05
        integer           :: k, state
        state = 1
        ! pre-proc
        p = params(cline, checkdistr=.false.) ! constants & derived constants produced, mode=2
        call b%build_general_tbox(p, cline)   ! general objects built
        call b%vol%read(p%vols(state))
        call b%vol%fwd_ft
        ! optimal low-pass filt
        if( cline%defined('fsc') )then
            if( file_exists(p%fsc) )then
                fsc   = file2rarr(p%fsc)
                optlp = fsc2optlp(fsc)
            else
                write(*,*) 'FSC file: ', trim(p%fsc), ' not in cwd'
                stop
            endif
            res = b%vol%get_res()
            call get_resolution( fsc, res, fsc05, fsc0143 )
            where(res < fsc0143) optlp = 0.
            call b%vol%apply_filter(optlp)
        else if( cline%defined('lp') )then
            call b%vol%bp(0., p%lp)
        else
            write(*,*) 'no method for low-pass filtering defined; give fsc or lp on command line'
            stop 'comple_commander_volops :: exec_postproc_vol'
        endif
        ! B-fact
        if( cline%defined('bfac') ) call b%vol%apply_bfac(p%bfac)
        ! masking
        call b%vol%bwd_ft
        if( p%automsk .eq. 'yes' )then
            call b%mskvol%automask3D(b%vol, p%msk, p%amsklp, p%mw, p%binwidth, p%edge, p%dens)
            p%masks(state) = 'automask_state'//int2str_pad(state,2)//p%ext
            call b%mskvol%write(p%masks(state))
        else
            call b%vol%mask(p%msk, 'soft')
        endif
        ! output
        p%outvol = add2fbody(trim(p%vols(state)), p%ext, 'pproc')
        call b%vol%write(p%outvol)
        call simple_end('**** SIMPLE_POSTPROC_VOL NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_postproc_vol
    
    subroutine exec_projvol( self, cline )
        use simple_image,          only: image
        use simple_projector_hlev, only: projvol
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
            call b%a%rnd_oris(p%trs, p%eullims)
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
        ! masking
        if(cline%defined('msk'))then
            call b%vol%mask(p%msk, 'soft')
        endif
        ! generate projections
        if( p%swap .eq. 'yes' ) call b%a%swape1e3
        if( cline%defined('top') )then
            imgs = projvol(b%vol, b%a, p, p%top)
            loop_end = p%top
        else
            imgs = projvol(b%vol, b%a, p)
            loop_end = p%nspace
        endif
        if( file_exists(p%outstk) ) call del_file(p%outstk)
        do i=1,loop_end
            if( cline%defined('oritab') .or. (p%rnd .eq. 'yes' .or. cline%defined('trs')) )then
                x = b%a%get(i, 'x')
                y = b%a%get(i, 'y')
                call imgs(i)%shift(x, y)
            endif
            if( p%neg .eq. 'yes' ) call imgs(i)%neg
            if( p%mirr .ne. 'no' ) call imgs(i)%mirror(p%mirr)
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
        logical, parameter                 :: debug=.false.
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
            ptcls = b%a%get_ptcls_in_state(istate)
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
            if( cline%defined('bfac') ) call b%vol%apply_bfac(p%bfac)
            call b%vol%write(p%outvol, del_if_exists=.true.)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_VOLOPS NORMAL STOP ****')
    end subroutine exec_volops
    
    subroutine exec_volume_smat( self, cline )
        use simple_projector, only: projector
        use simple_ori,       only: ori
        use simple_volpft_srch ! singleton
        class(volume_smat_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params), target :: p
        integer              :: funit, io_stat, cnt, npairs, npix, nvols, box_sc, loc(1)
        integer              :: ivol, jvol, ldim(3), alloc_stat, ipair, ifoo, i, spat_med
        integer              :: furthest_from_spat_med
        real                 :: smpd_sc, scale, corr_max, corr_min, spat_med_corr
        real                 :: furthest_from_spat_med_corr
        type(projector)      :: vol1, vol2
        type(ori)            :: o
        logical, parameter   :: debug=.false.
        real,                  allocatable :: corrmat(:,:), corrs(:), corrs_avg(:)
        integer,               allocatable :: pairs(:,:)
        character(len=STDLEN), allocatable :: vollist(:)
        character(len=:),      allocatable :: fname
        complex,               allocatable :: cmat(:,:,:)
        p = params(cline, .false.)              ! constants & derived constants produced
        call read_filetable(p%vollist, vollist) ! reads in list of volumes
        nvols  = size(vollist)
        npairs = (nvols*(nvols-1))/2
        ! find logical dimension & make volumes for matching
        call find_ldim_nptcls(vollist(1), ldim, ifoo)
        if( debug ) write(*,*) 'found logical dimension: ', ldim
        if( cline%defined('part') )then
            npairs = p%top-p%fromp+1
            if( debug ) print *, 'allocating this number of similarities: ', npairs
            allocate(corrs(p%fromp:p%top), pairs(p%fromp:p%top,2), stat=alloc_stat)
            call alloc_err('In: simple_volume_smat, 1', alloc_stat)
            ! read the pairs
            funit = get_fileunit()
            allocate(fname, source='pairs_part'//int2str_pad(p%part,p%numlen)//'.bin')
            if( .not. file_exists(fname) )then
                write(*,*) 'file: ', fname, 'does not exist!'
                write(*,*) 'If all pair_part* are not in cwd, please execute simple_split_pairs to generate the required files'
                stop 'I/O error; simple_volume_smat'
            endif
            open(unit=funit, status='OLD', action='READ', file=fname, access='STREAM')
            if( debug ) print *, 'reading pairs in range: ', p%fromp, p%top
            read(unit=funit,pos=1,iostat=io_stat) pairs(p%fromp:p%top,:)
            ! Check if the read was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(simple_volume_smat): I/O error ', io_stat, ' when reading file: ', fname
                stop 'I/O error; simple_volume_smat'
            endif
            close(funit)
            deallocate(fname)
            cnt = 0
            do ipair=p%fromp,p%top
                cnt = cnt + 1
                call progress(cnt, npairs)
                ivol = pairs(ipair,1)
                jvol = pairs(ipair,2)
                call read_and_prep_vols( ivol, jvol )
                o = volpft_srch_minimize() 
                corrs(ipair) = o%get('corr')
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
                stop 'I/O error; simple_volume_smat'
            endif
            close(funit)
            deallocate(fname, corrs, pairs)
        else
            ! generate similarity matrix
            allocate(corrmat(nvols,nvols), corrs_avg(nvols), stat=alloc_stat)
            call alloc_err('In: simple_volume_smat, 2', alloc_stat)
            corrmat = -1.
            forall(i=1:nvols) corrmat(i,i) = 1.0
            cnt = 0
            corr_max = -1.0
            corr_min =  1.0
            do ivol=1,nvols - 1
                do jvol=ivol + 1,nvols
                    cnt = cnt + 1
                    call progress(cnt, npairs)
                    call read_and_prep_vols( ivol, jvol )
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
            write(*,'(a,1x,f7.4)') 'MAX VOL PAIR CORR          :', corr_max
            write(*,'(a,1x,f7.4)') 'MIN VOL PAIR CORR          :', corr_min
            write(*,'(a,1x,i7)'  ) 'SPATIAL MEDIAN             :', spat_med 
            write(*,'(a,1x,f7.4)') 'SPATIAL MEDIAN CORR        :', spat_med_corr
            write(*,'(a,1x,i7)'  ) 'FURTHEST FROM SPAT MED     :', furthest_from_spat_med
            write(*,'(a,1x,f7.4)') 'FURTHEST FROM SPAT MED CORR:', furthest_from_spat_med_corr
            funit = get_fileunit()
            open(unit=funit, status='REPLACE', action='WRITE', file='vol_smat.bin', access='STREAM')
            write(unit=funit,pos=1,iostat=io_stat) corrmat
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to vol_smat.bin'
                stop 'I/O error; simple_volume_smat'
            endif
            close(funit)
            deallocate(corrmat)
        endif     
        ! end gracefully
        call simple_end('**** SIMPLE_VOLUME_SMAT NORMAL STOP ****')

        contains

            subroutine read_and_prep_vols( ivol, jvol )
                use simple_magic_boxes, only: autoscale
                integer, intent(in) :: ivol, jvol
                real :: smpd_target
                call vol1%new(ldim,p%smpd)
                call vol2%new(ldim,p%smpd)
                call vol1%read(vollist(ivol))
                call vol2%read(vollist(jvol))
                if( p%boxmatch < p%box )then
                    call vol1%clip_inplace([p%boxmatch,p%boxmatch,p%boxmatch])
                    call vol2%clip_inplace([p%boxmatch,p%boxmatch,p%boxmatch]) 
                endif
                call vol1%mask(p%msk,'soft')
                call vol2%mask(p%msk,'soft')
                call vol1%fwd_ft
                call vol2%fwd_ft
                smpd_target = p%lp*LP2SMPDFAC
                call autoscale(p%boxmatch, p%smpd, smpd_target, box_sc, smpd_sc, scale)
                call vol1%clip_inplace([box_sc,box_sc,box_sc])
                call vol2%clip_inplace([box_sc,box_sc,box_sc])
                call vol1%set_smpd(smpd_sc)
                call vol2%set_smpd(smpd_sc)
                call volpft_srch_init(vol1, vol2, p%hp, p%lp)
            end subroutine read_and_prep_vols

    end subroutine exec_volume_smat

end module simple_commander_volops
