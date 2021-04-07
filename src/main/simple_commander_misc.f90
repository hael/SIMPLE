! concrete commander: miscallenaous routines
module simple_commander_misc
include 'simple_lib.f08'
include "starfile/starfile_enum.inc"
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_oris,           only: oris
use simple_ori,            only: ori
use simple_sym,            only: sym
use simple_projector_hlev, only: rotvol
use simple_sp_project,     only: sp_project
use simple_image,          only: image
use simple_binoris_io,     only: binwrite_oritab
use simple_builder,        only: builder
use simple_parameters,     only: parameters
implicit none

public :: cluster_smat_commander
public :: masscen_commander
public :: print_fsc_commander
public :: print_magic_boxes_commander
public :: shift_commander
public :: stk_corr_commander
public :: kstest_commander
public :: mkdir_commander
public :: remoc_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: cluster_smat_commander
  contains
    procedure :: execute      => exec_cluster_smat
end type cluster_smat_commander
type, extends(commander_base) :: masscen_commander
  contains
    procedure :: execute      => exec_masscen
end type masscen_commander
type, extends(commander_base) :: print_fsc_commander
  contains
    procedure :: execute       => exec_print_fsc
end type print_fsc_commander
type, extends(commander_base) :: print_magic_boxes_commander
  contains
    procedure :: execute       => exec_print_magic_boxes
end type print_magic_boxes_commander
type, extends(commander_base) :: shift_commander
  contains
    procedure :: execute       => exec_shift
end type shift_commander
type, extends(commander_base) :: stk_corr_commander
  contains
    procedure :: execute       => exec_stk_corr
end type stk_corr_commander
type, extends(commander_base) :: kstest_commander
  contains
    procedure :: execute       => exec_kstest
end type kstest_commander
type, extends(commander_base) :: mkdir_commander
  contains
    procedure :: execute       => exec_mkdir
end type mkdir_commander
type, extends(commander_base) :: remoc_commander
  contains
    procedure :: execute       => exec_remoc
end type remoc_commander

contains

    subroutine exec_cluster_smat( self, cline )
        use simple_cluster_shc,   only: cluster_shc
        use simple_cluster_valid, only: cluster_valid
        class(cluster_smat_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)    :: params
        type(builder)       :: build
        type(cluster_shc)   :: shcc
        type(cluster_valid) :: cvalid
        real, allocatable   :: smat(:,:)
        integer             :: funit, io_stat, ncls_min, loc(1), ncls_stop, icls
        integer             :: pop, ncls, irestart, ntot, cnt=0, numlen
        real                :: avg_ratio, min_ratio, ratio, x, sim
        real, allocatable   :: validinds(:)
        integer, parameter  :: NRESTARTS=10
        logical             :: done=.false.
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! obtain similarity matrix
        allocate(smat(params%nptcls,params%nptcls), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('In: simple_cluster_smat, 1',alloc_stat)
        smat = 1.
        call fopen(funit, status='OLD', action='READ', file=params%fname, access='STREAM',iostat=io_stat)
        call fileiochk('commander_misc; cluster_smat fopen', io_stat)
        read(unit=funit,pos=1,iostat=io_stat) smat
        if( io_stat .ne. 0 )then
            write(logfhandle,'(a,i0,a)') 'I/O error ', io_stat, ' when reading: ', params%fname
            THROW_HARD('I/O')
        endif
        call fclose(funit,errmsg='commander_misc; cluster_smat fclose ')
        allocate(validinds(2:params%ncls), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In: simple_commander_misc:: cluster_smat",alloc_stat)
        validinds = 0
        ntot = (params%ncls-1)*NRESTARTS
        cnt = 0
        numlen = len(int2str(params%ncls))
        do ncls=2,params%ncls
            avg_ratio = 0.
            min_ratio = huge(x)
            do irestart=1,NRESTARTS
                cnt = cnt+1
                call progress(cnt,ntot)
                call shcc%new(params%nptcls, ncls, smat, build%spproj_field)
                call shcc%shc(.false., params%label, sim)
                call cvalid%new(build%spproj_field, ncls, params%label, smat)
                ratio = cvalid%ratio_index()
                avg_ratio = avg_ratio+ratio
                if( ratio < min_ratio )then
                    min_ratio = ratio
                    call build%spproj_field%write('clustering_shc_ncls'//int2str_pad(ncls,numlen)//trim(TXT_EXT), [1,params%nptcls])
                endif
            end do
            validinds(ncls) = avg_ratio/real(NRESTARTS)
        end do
        ncls_stop = 0
        done = .false.
        do ncls=2,params%ncls
            write(logfhandle,'(a,1x,f9.3,8x,a,1x,i3)') 'COHESION/SEPARATION RATIO INDEX: ', validinds(ncls), ' NCLS: ', ncls
            call build%spproj_field%read('clustering_shc_ncls'//int2str_pad(ncls,numlen)//trim(TXT_EXT), [1,build%spproj_field%get_noris()])
            do icls=1,ncls
                pop = build%spproj_field%get_pop(icls, params%label)
                write(logfhandle,'(a,3x,i5,1x,a,1x,i3)') '  CLUSTER POPULATION:', pop, 'CLUSTER:', icls
            end do
            write(logfhandle,'(a)') '***************************************************'
            if( ncls < params%ncls )then
                if( validinds(ncls+1) >= validinds(ncls) .and. .not. done )then
                    ncls_stop = ncls
                    done = .true.
                endif
            endif
        end do
        loc = minloc(validinds)
        ncls_min = loc(1)+1
        write(logfhandle,'(a,i3)') 'NUMBER OF CLUSTERS FOUND BY MINIMIZING RATIO INDEX: ', ncls_min
        if( ncls_stop /= 0 )then
            write(logfhandle,'(a,i3)') 'NUMBER OF CLUSTERS FOUND BY STOPPING CRITERIUM:     ', ncls_stop
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER_SMAT NORMAL STOP ****')
    end subroutine exec_cluster_smat

    !> centers base on centre of mass
     subroutine exec_masscen( self, cline )
        use simple_procimgfile, only: masscen_imgfile
        class(masscen_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(parameters) :: params
        call params%new(cline)
        params%cenlp = params%lp
        ! center of mass centering
        call masscen_imgfile(params%stk, params%outstk, params%smpd, params%lp, params%msk)
        ! end gracefully
        call simple_end('**** SIMPLE_MASSCEN NORMAL STOP ****')
    end subroutine exec_masscen

    !>  for printing the binary FSC files produced by PRIME3D
    subroutine exec_print_fsc( self, cline )
        class(print_fsc_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters)  :: params
        type(image)       :: img
        real, allocatable :: res(:), fsc(:)
        integer           :: k
        real              :: res0143, res05
        call params%new(cline)
        call img%new([params%box,params%box,1], params%smpd)
        res = img%get_res()
        fsc = file2rarr(params%fsc)
        do k=1,size(fsc)
        write(logfhandle,'(A,1X,F6.2,1X,A,1X,F15.3)') '>>> RESOLUTION:', res(k), '>>> FSC:', fsc(k)
        end do
        ! get & print resolution
        call get_resolution(fsc, res, res05, res0143)
        write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res0143
        write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res05
        ! end gracefully
        call simple_end('**** SIMPLE_PRINT_FSC NORMAL STOP ****')
    end subroutine exec_print_fsc

    !> for printing magic box sizes (fast FFT)
    subroutine exec_print_magic_boxes( self, cline )
        class(print_magic_boxes_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters) :: params
        call params%new(cline)
        call print_magic_box_range(params%smpd, params%moldiam )
        ! end gracefully
        call simple_end('**** SIMPLE_PRINT_MAGIC_BOXES NORMAL STOP ****')
    end subroutine exec_print_magic_boxes

    !> for shifting a stack according to shifts in oritab
    subroutine exec_shift( self, cline )
        use simple_procimgfile, only: shift_imgfile
        class(shift_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        call shift_imgfile(params%stk, params%outstk, build%spproj_field, params%smpd, params%mul)
        call build%spproj_field%zero_shifts
        call binwrite_oritab('shiftdoc'//trim(METADATA_EXT), build%spproj, build%spproj_field, [1,params%nptcls])
        call simple_end('**** SIMPLE_SHIFT NORMAL STOP ****')
    end subroutine exec_shift

    subroutine exec_stk_corr( self, cline )
        class(stk_corr_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(parameters)     :: params
        type(builder)        :: build
        logical, allocatable :: l_mask(:,:,:)
        integer :: i
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        build%img = 1.
        call build%img%mask(params%msk, 'hard')
        l_mask = build%img%bin2logical()
        do i=1,params%nptcls
            call build%img%read(params%stk, i)
            call build%img_copy%read(params%stk2, i)
            if( cline%defined('lp') )then
                call build%img%norm
                call build%img_copy%norm
                call build%img%mask(params%msk, 'soft')
                call build%img_copy%mask(params%msk, 'soft')
                call build%img%fft
                call build%img_copy%fft
                call build%img%bp(0.,params%lp)
                call build%img_copy%bp(0.,params%lp)
                call build%img%ifft
                call build%img_copy%ifft
            endif
            write(logfhandle,'(I6,F8.3)')i,build%img%real_corr(build%img_copy, l_mask)
        enddo
        ! end gracefully
        call simple_end('**** SIMPLE_CONVERT NORMAL STOP ****')
    end subroutine exec_stk_corr

    subroutine exec_kstest( self, cline )
        class(kstest_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters)   :: params
        integer            :: ndat1, ndat2
        real               :: ksstat, prob, ave1, sdev1, var, ave2, sdev2
        real, allocatable  :: dat1(:), dat2(:)
        logical            :: err
        call params%new(cline)
        call read_nrs_dat(params%infile,  dat1, ndat1)
        call read_nrs_dat(params%infile2, dat2, ndat2)
        write(logfhandle,'(a)') '>>> STATISTICS OF THE TWO DISTRIBUTIONS'
        call moment(dat1, ave1, sdev1, var, err)
        call moment(dat2, ave2, sdev2, var, err)
        write(logfhandle,'(a,1x,f4.2,1x,f4.2)') 'mean & sdev for infile : ', ave1, sdev1
        write(logfhandle,'(a,1x,f4.2,1x,f4.2)') 'mean & sdev for infile2: ', ave2, sdev2
        write(logfhandle,'(a)') '>>> KOLMOGOROV-SMIRNOV TEST TO DEDUCE EQUIVALENCE OR NON-EQUIVALENCE BETWEEN TWO DISTRIBUTIONS'
        call kstwo(dat1, ndat1, dat2, ndat2, ksstat, prob)
        write(logfhandle,'(a,1x,f4.2)') 'K-S statistic = ', ksstat
        write(logfhandle,'(a,1x,f4.2)') 'P             = ', prob
        write(logfhandle,'(a)') 'P represents the significance level for the null hypothesis that the two data sets are drawn from the same distribution'
        write(logfhandle,'(a)') 'Small P values show that the cumulative distribution functions of the two data sets differ significantly'
        ! end gracefully
        call simple_end('**** SIMPLE_KSTEST NORMAL STOP ****')

        contains

            subroutine read_nrs_dat( filename, arr, ndat )
                use simple_nrtxtfile, only: nrtxtfile
                character(len=*),  intent(in)  :: filename
                real, allocatable, intent(out) :: arr(:)
                integer,           intent(out) :: ndat
                integer            :: ndatlines, nrecs, i, j, cnt
                real, allocatable  :: line(:)
                type(nrtxtfile)    :: nrsfile
                call nrsfile%new(filename, 1)
                ndatlines = nrsfile%get_ndatalines()
                nrecs     = nrsfile%get_nrecs_per_line()
                ndat = ndatlines * nrecs
                allocate( line(nrecs), arr(ndat) )
                cnt = 0
                do i=1,ndatlines
                    call nrsfile%readNextDataLine(line)
                    do j=1,nrecs
                        cnt = cnt + 1
                        arr(cnt) = line(j)
                    end do
                end do
                deallocate (line)
                call nrsfile%kill
            end subroutine read_nrs_dat

    end subroutine exec_kstest

    !>  dsym_cylinder search intended for symmetry of order D
    subroutine exec_dsym_volinit( dsym_os, cylinder)
        use simple_parameters,   only: params_glob
        use simple_oris,         only: oris
        use simple_image,        only: image
        use simple_sym,          only: sym
        use simple_segmentation, only: otsu_robust_fast
        class(oris),   intent(inout) :: dsym_os
        class(image),  intent(inout) :: cylinder
        type(image)          :: read_img, img_msk, dist_img, roavg_img, topview
        type(sym)            :: se
        real,    allocatable :: corrs(:), forsort(:), e2(:), radii(:)
        integer, allocatable :: labels(:)
        logical, allocatable :: l_msk(:,:,:)
        integer  :: halfnoris, cnt1, cnt2, i, l, noris
        real     :: cen1, cen2, sum1, sum2, sumvals, sdev, sdev_noise
        real     :: minmax(2), width, height, sh1(3), ang, thresh(3)
        if(params_glob%pgrp(1:1).eq.'d' .or. params_glob%pgrp(1:1).eq.'D')then
            ang = 360. / real(se%get_nsym()/2)
        else if(params_glob%pgrp(1:1).eq.'c' .or. params_glob%pgrp(1:1).eq.'C')then
            ang = 360. / real(se%get_nsym())
        else
            THROW_HARD('Not intended for platonic point-groups')
        endif
        ! init
        noris = dsym_os%get_noris()
        allocate(corrs(noris),  source=-1.)
        allocate(radii(noris),e2(noris),  source=0.)
        allocate(labels(noris), source=0)
        call se%new(params_glob%pgrp)
        call cylinder%new( [params_glob%box, params_glob%box, params_glob%box], params_glob%smpd)
        call read_img%new( [params_glob%box, params_glob%box, 1], params_glob%smpd)
        call topview%new( [params_glob%box, params_glob%box, 1], params_glob%smpd)
        call roavg_img%new([params_glob%box, params_glob%box, 1], params_glob%smpd)
        call dist_img%new( [params_glob%box, params_glob%box, 1], params_glob%smpd)
        call dist_img%cendist
        ! prep mask
        call img_msk%new([params_glob%box,params_glob%box,1],params_glob%smpd)
        img_msk = 1.
        call img_msk%mask(params_glob%msk, 'hard')
        l_msk = img_msk%bin2logical()
        call img_msk%kill
        ! centers, calculates self to rotational averages images & radii
        do i = 1, noris
            call read_img%read(params_glob%stk,i)
            call read_img%noise_norm(l_msk, sdev_noise)
            ! center image
            sh1 = read_img%calc_shiftcen(params_glob%cenlp, params_glob%msk)
            call read_img%shift(-sh1)
            call dsym_os%set_shift(i, -sh1)
            ! rotational image
            call read_img%roavg(nint(ang), roavg_img)
            call roavg_img%write('roavg.mrc',i)
            corrs(i) = read_img%real_corr(roavg_img, l_msk)
            ! radii
            call read_img%bp(0., params_glob%cenlp)
            call read_img%mask(params_glob%msk, 'hard')
            call otsu_robust_fast(read_img, is2d=.true., noneg=.false., thresh=thresh)
            call read_img%write('bin.mrc',i)
            call read_img%mul(dist_img)
            minmax = read_img%minmax()
            radii(i) = min(params_glob%msk, minmax(2))
        enddo
        call dsym_os%set_all('corr', corrs)
        ! indentify views along z and x-/y- axes
        forsort = corrs
        call hpsort(forsort)
        halfnoris = nint(real(noris)/2.)
        cen1      = sum(forsort(1:halfnoris)) / real(halfnoris)
        cen2      = sum(forsort(halfnoris+1:noris)) / real(noris-halfnoris)
        sumvals   = sum(forsort)
        ! do 100 iterations of k-means
        do l = 1, 100
            sum1 = 0.
            cnt1 = 0
            do i = 1, noris
                if( (cen1-forsort(i))**2. < (cen2-forsort(i))**2. )then
                    cnt1 = cnt1 + 1
                    sum1 = sum1 + forsort(i)
                endif
            end do
            cnt2 = noris - cnt1
            sum2 = sumvals - sum1
            cen1 = sum1 / real(cnt1)
            cen2 = sum2 / real(cnt2)
        end do
        ! label views
        if(cen1 > cen2)then
            ! cen1: along z-axis
            labels = 2
            where( (cen1-corrs)**2. < (cen2-corrs)**2. )labels = 1
        else
            labels = 1
            where( (cen1-corrs)**2. < (cen2-corrs)**2. )labels = 2
        endif
        e2 = 90.
        where(labels == 1) e2 = 0.
        call dsym_os%set_all('e2', e2)
        call dsym_os%set_all('class', real(labels))
        ! rotates input oris to asymmetric unit
        call se%rotall_to_asym(dsym_os)
        ! estimate height and cylinder radius (0.9 to account for overestimation)
        width  = 0.9 * (sum(radii, mask=(labels==1)) / real(count(labels==1)))
        height = 0.9 * (2. * sum(radii, mask=(labels==2)) / real(count(labels==2)))
        ! dummy top view
        topview = 0.
        do i=1,noris
            if(labels(i)==1)then
                call read_img%read(params_glob%stk,i)
                call read_img%noise_norm(l_msk, sdev_noise)
                call read_img%rtsq(0., -dsym_os%get(i,'x'), -dsym_os%get(i,'y'))
                call topview%add(read_img)
            endif
        enddo
        call topview%div( real(count(labels==1)) )
        call topview%roavg(nint(ang), roavg_img)
        topview = roavg_img
        call topview%norm()
        call topview%mask(params_glob%msk, 'soft')
        cylinder = 0.
        do i=1,params_glob%box
            if( abs( real(i-1)-real(params_glob%box)/2. ) < height/2.)then
                call cylinder%set_slice( i, topview )
            endif
        enddo
        ! cleanup
        deallocate(corrs, e2, radii, labels)
        call dist_img%kill
        call read_img%kill
        call se%kill
        call roavg_img%kill
        call topview%kill
    end subroutine exec_dsym_volinit

    subroutine exec_mkdir( self, cline )
        class(mkdir_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters) :: params
        call cline%set('mkdir', 'yes')
        call params%new(cline)
    end subroutine exec_mkdir

    subroutine exec_remoc( self, cline )
        use simple_starfile_wrappers
        use simple_projector
        use simple_estimate_ssnr, only: fsc2optlp_sub
        class(remoc_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        class(str4arr),   allocatable :: mics(:),mics_from_ptcls(:), star_movies(:)
        character(len=:), allocatable :: fname, movie, vol_even_fname, vol_odd_fname, mask_fname, mic_fname
        real,             allocatable :: fsc(:), optfilter(:)
        type(parameters) :: params
        type(ctfparams)  :: ctf_glob, ctf_mic
        type(oris)       :: os
        type(projector)  :: vol_even, vol_odd
        type(image)      :: vol, vol_mask
        real    :: mov_smpd
        integer :: ldim(3), ldim_mov(3), nmics, nmovies, filtsz
        logical :: l_scale
        call cline%set('mkdir',       'yes')
        if( .not. cline%defined('trs')           ) call cline%set('trs',           20.)
        if( .not. cline%defined('lpstop')        ) call cline%set('lpstop',         5.)
        if( .not. cline%defined('bfac')          ) call cline%set('bfac',          50.)
        if( .not. cline%defined('lplim_crit')    ) call cline%set('lplim_crit',    0.5)
        if( .not. cline%defined('wcrit')         ) call cline%set('wcrit',   'softmax')
        if( .not. cline%defined('eer_fraction')  ) call cline%set('eer_fraction',  20.)
        if( .not. cline%defined('eer_upsampling')) call cline%set('eer_upsampling', 1.)
        call cline%set('oritype', 'mic')
        call params%new(cline)
        ! Particles
        call parse_particles( params%star_ptcl, os, mics, ctf_glob, params%box, params%nptcls )
        params%smpd = ctf_glob%smpd
        ldim = params%box
        ! Micrographs
        call parse_mics( params%star_mic, mics, star_movies, ctf_mic, mov_smpd )
        params%smpd = ctf_glob%smpd
        nmics   = size(mics)
        nmovies = size(star_movies)
        if( abs(params%smpd-ctf_mic%smpd) > 0.01 ) THROW_HARD('Inconsistent pixel sizes between micrographs and particles!')
        params%scale = mov_smpd / ctf_glob%smpd
        l_scale = abs(params%scale-1.) > 0.01
        if( l_scale )then
            params%alpha = 1./params%scale
            ldim_mov = round2even(real(ldim)/params%scale)
        else
            params%alpha = 1.
            ldim_mov     = ldim
        endif
        ! Model
        call parse_model(params%star_model, vol_even_fname, vol_odd_fname, mask_fname, fsc)
        mask_fname = get_fpath(params%star_model)//'/'//trim(mask_fname)
        filtsz = size(fsc)
        allocate(optfilter(filtsz),source=0.)
        call fsc2optlp_sub(filtsz, fsc, optfilter)
        ! volumes prep
        call vol%new(ldim,params%smpd)
        call vol_mask%new(ldim,params%smpd)
        call vol%read(vol_even_fname)
        call vol_mask%read(mask_fname)
        call vol%mul(vol_mask)
        if( l_scale )then
            call vol_even%new(ldim_mov,mov_smpd)
            call vol%pad(vol_even)
        else
            vol_even = vol
        endif
        call vol%read(vol_odd_fname)
        call vol%mul(vol_mask)
        call vol_mask%kill
        if( l_scale )then
            call vol_odd%new(ldim_mov,mov_smpd)
            call vol%pad(vol_odd)
        else
            vol_odd = vol
        endif
        call vol%kill
        call vol_even%fft
        call vol_odd%fft
        call vol_even%apply_filter(optfilter)
        call vol_odd%apply_filter(optfilter)
        call vol_even%expand_cmat(1.,norm4proj=.true.)
        call vol_odd%expand_cmat(1.,norm4proj=.true.)
        ! end gracefully
        call simple_end('**** SIMPLE_REMOC NORMAL STOP ****')
        contains

            integer function parse_int( table, emdl_id )
                class(starfile_table_type) :: table
                integer, intent(in)        :: emdl_id
                logical :: l_ok
                l_ok = starfile_table__getValue_int(table, emdl_id, parse_int)
                if(.not.l_ok) THROW_HARD('Missing value in table!')
            end function parse_int

            real function parse_double( table, emdl_id )
                class(starfile_table_type) :: table
                integer, intent(in)        :: emdl_id
                real(dp) :: v
                logical  :: l_ok
                l_ok = starfile_table__getValue_double(table, emdl_id, v)
                if(.not.l_ok) THROW_HARD('Missing value in table!')
                parse_double = real(v)
            end function parse_double

            subroutine parse_string( table, emdl_id, string )
                class(starfile_table_type)                 :: table
                integer,                       intent(in)  :: emdl_id
                character(len=:), allocatable, intent(out) :: string
                logical :: l_ok
                l_ok = starfile_table__getValue_string(table, emdl_id, string)
                if(.not.l_ok) THROW_HARD('Missing value in table!')
            end subroutine parse_string

            subroutine parse_particles( fname, os, mics, ctf, box, nptcls_out )
                character(len=*),               intent(in)  :: fname
                class(oris),                    intent(out) :: os
                class(str4arr),    allocatable, intent(out) :: mics(:)
                class(ctfparams),               intent(out) :: ctf
                integer,                        intent(out) :: box, nptcls_out
                type(str4arr),    allocatable :: names(:)
                type(str4arr)    ::fname_here
                type(starfile_table_type)     :: star_table
                type(ctfparams)  :: ctfparms
                real             :: euls(3),pos(2),shift(2)
                integer          :: i,n,ind,eo
                integer(C_long)  :: num_objs, object_id
                fname_here%str = trim(fname)//C_NULL_CHAR
                call starfile_table__new(star_table)
                call starfile_table__getnames(star_table, fname_here%str, names)
                n = size(names)
                do i =1,n
                    call starfile_table__read(star_table, fname_here%str, names(i)%str )
                    select case(trim(names(i)%str))
                    case('optics')
                        ! ignoring multiple optics groups for now
                        object_id = starfile_table__firstobject(star_table) ! base 0
                        ! ctf
                        ctf%fraca = parse_double(star_table, EMDL_CTF_Q0)
                        ctf%cs    = parse_double(star_table, EMDL_CTF_CS)
                        ctf%kv    = parse_double(star_table, EMDL_CTF_VOLTAGE)
                        ctf%smpd  = parse_double(star_table, EMDL_IMAGE_PIXEL_SIZE)
                        ! dimensions
                        box = parse_int(star_table, EMDL_IMAGE_SIZE)
                    case('particles')
                        object_id  = starfile_table__firstobject(star_table) ! base 0
                        num_objs   = starfile_table__numberofobjects(star_table)
                        nptcls_out = int(num_objs - object_id)
                        call os%new(nptcls_out)
                        allocate(mics(nptcls_out))
                        ind = 0
                        do while( (object_id < num_objs) .and. (object_id >= 0) )
                            ind = ind+1
                            ! Micrograph name
                            call parse_string(star_table, EMDL_MICROGRAPH_NAME, mics(ind)%str)
                            ! picking coordinates
                            pos(1) = parse_double(star_table, EMDL_IMAGE_COORD_X) - box/2 - 1
                            pos(2) = parse_double(star_table, EMDL_IMAGE_COORD_Y) - box/2 - 1
                            ! CTF
                            ctfparms%dfx    = parse_double(star_table, EMDL_CTF_DEFOCUSU) / 10000.0 ! Microns
                            ctfparms%dfy    = parse_double(star_table, EMDL_CTF_DEFOCUSU) / 10000.0 ! microns
                            ctfparms%angast = parse_double(star_table, EMDL_CTF_DEFOCUS_ANGLE)
                            ! e/o
                            eo = parse_int(star_table, EMDL_PARTICLE_RANDOM_SUBSET) - 1
                            ! Alignement parameters
                            euls(1)  = parse_double(star_table, EMDL_ORIENT_ROT)
                            euls(2)  = parse_double(star_table, EMDL_ORIENT_TILT)
                            euls(3)  = parse_double(star_table, EMDL_ORIENT_PSI)
                            shift(1) = parse_double(star_table, EMDL_ORIENT_ORIGIN_X_ANGSTROM) / ctf_glob%smpd ! in pixels
                            shift(2) = parse_double(star_table, EMDL_ORIENT_ORIGIN_X_ANGSTROM) / ctf_glob%smpd ! in pixels
                            ! transfer to oris object
                            call os%set_euler(ind,   euls)
                            call os%set_shift(ind,   shift)
                            call os%set(ind,'dfx',   ctfparms%dfx)
                            call os%set(ind,'dfy',   ctfparms%dfy)
                            call os%set(ind,'angast',ctfparms%angast)
                            call os%set(ind,'xpos',  pos(1))
                            call os%set(ind,'ypos',  pos(2))
                            call os%set(ind,'eo',    real(eo))
                            ! done
                            object_id = starfile_table__nextobject(star_table)
                        end do
                    case DEFAULT
                        THROW_HARD('Invalid table: '//trim(names(i)%str))
                    end select
                enddo
                call starfile_table__delete(star_table)
            end subroutine parse_particles
    
            subroutine parse_mics( fname, mics, star_movies, ctf, mov_smpd )
                character(len=*),                intent(in) :: fname
                class(str4arr),    allocatable, intent(out) :: mics(:), star_movies(:)
                class(ctfparams),               intent(out) :: ctf
                real,                           intent(out) :: mov_smpd
                type(str4arr),    allocatable :: names(:)
                type(starfile_table_type)     :: star_table
                integer          :: i,n,nmics,ind
                integer(C_long)  :: num_objs, object_id
                call starfile_table__new(star_table)
                call starfile_table__getnames(star_table, trim(fname)//C_NULL_CHAR, names)
                n = size(names)
                do i =1,n
                    call starfile_table__read(star_table, trim(fname)//C_NULL_CHAR, names(i)%str )
                    select case(trim(names(i)%str))
                    case('optics')
                        ! ignoring multiple optics groups for now
                        object_id = starfile_table__firstobject(star_table) ! base 0
                        ! ctf
                        ctf%fraca = parse_double(star_table, EMDL_CTF_Q0)
                        ctf%cs    = parse_double(star_table, EMDL_CTF_CS)
                        ctf%kv    = parse_double(star_table, EMDL_CTF_VOLTAGE)
                        ctf%smpd  = parse_double(star_table, EMDL_MICROGRAPH_PIXEL_SIZE)
                        mov_smpd  = parse_double(star_table, EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE)
                    case('micrographs')
                        object_id = starfile_table__firstobject(star_table) ! base 0
                        num_objs  = starfile_table__numberofobjects(star_table)
                        nmics     = int(num_objs - object_id)
                        allocate(mics(nmics),star_movies(nmics))
                        ind = 0
                        do while( (object_id < num_objs) .and. (object_id >= 0) )
                            ind = ind+1
                            call parse_string(star_table, EMDL_MICROGRAPH_NAME, mics(ind)%str)
                            call parse_string(star_table, EMDL_MICROGRAPH_METADATA_NAME, star_movies(ind)%str)
                            ! done
                            object_id = starfile_table__nextobject(star_table)
                        end do
                    case DEFAULT
                        THROW_HARD('Invalid table: '//trim(names(i)%str))
                    end select
                enddo
                call starfile_table__delete(star_table)
            end subroutine parse_mics

            subroutine parse_model( fname, vol_even, vol_odd, mask, fsc )
                character(len=*),              intent(in)  :: fname
                character(len=:), allocatable, intent(out) :: vol_even, vol_odd, mask
                real,             allocatable, intent(out) :: fsc(:)
                type(str4arr),    allocatable :: names(:)
                type(starfile_table_type)     :: star_table
                integer          :: i,n,ind,nn
                integer(C_long)  :: num_objs, object_id
                call starfile_table__new(star_table)
                call starfile_table__getnames(star_table, trim(fname)//C_NULL_CHAR, names)
                n = size(names)
                do i =1,n
                    call starfile_table__read(star_table, trim(fname)//C_NULL_CHAR, names(i)%str )
                    select case(trim(names(i)%str))
                    case('general')
                        object_id = starfile_table__firstobject(star_table)
                        ! e/o volumes, mask
                        call parse_string(star_table, EMDL_POSTPROCESS_UNFIL_HALFMAP1, vol_even)
                        call parse_string(star_table, EMDL_POSTPROCESS_UNFIL_HALFMAP2, vol_odd)
                        call parse_string(star_table, EMDL_MASK_NAME, mask)
                    case('fsc')
                        object_id = starfile_table__firstobject(star_table)
                        num_objs  = starfile_table__numberofobjects(star_table)
                        nn        = int(num_objs - object_id)
                        if(nn /= params%box/2+1 ) THROW_HARD('Inconsistent image & fsc dimensions!')
                        allocate(fsc(params%box/2),source=0.)
                        do while( (object_id < num_objs) .and. (object_id >= 0) )
                            ind = parse_int(star_table, EMDL_SPECTRAL_IDX)
                            if( ind /= 0 ) fsc(ind) = parse_double(star_table, EMDL_POSTPROCESS_FSC_TRUE)
                            ! done
                            object_id = starfile_table__nextobject(star_table)
                        end do
                    case('guinier')
                        ! Ignored for now
                    case DEFAULT
                        THROW_HARD('Invalid table: '//trim(names(i)%str))
                    end select
                enddo
                call starfile_table__delete(star_table)
            end subroutine parse_model

    end subroutine exec_remoc

end module simple_commander_misc
