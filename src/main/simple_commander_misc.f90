! concrete commander: miscallenaous routines
module simple_commander_misc
include 'simple_lib.f08'
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
public :: print_dose_weights_commander
public :: print_fsc_commander
public :: print_magic_boxes_commander
public :: res_commander
public :: shift_commander
public :: stk_corr_commander
public :: kstest_commander
public :: mkdir_commander
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
type, extends(commander_base) :: print_dose_weights_commander
  contains
    procedure :: execute       => exec_print_dose_weights
end type print_dose_weights_commander
type, extends(commander_base) :: print_fsc_commander
  contains
    procedure :: execute       => exec_print_fsc
end type print_fsc_commander
type, extends(commander_base) :: print_magic_boxes_commander
  contains
    procedure :: execute       => exec_print_magic_boxes
end type print_magic_boxes_commander
type, extends(commander_base) :: res_commander
  contains
    procedure :: execute       => exec_res
end type res_commander
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

    !> for printing the dose weights applied to individual frames
    subroutine exec_print_dose_weights( self, cline )
        use simple_estimate_ssnr, only: acc_dose2filter
        class(print_dose_weights_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters)  :: params
        real, allocatable :: filter(:)
        type(image)       :: dummy_img
        real              :: time_per_frame, current_time, acc_dose
        integer           :: iframe, find, filtsz
        call params%new(cline)
        call dummy_img%new([params%box,params%box,1], params%smpd)
        time_per_frame = params%exp_time/real(params%nframes)
        filtsz = dummy_img%get_filtsz()
        do iframe=1,params%nframes
            current_time = real(iframe)*time_per_frame
            acc_dose     = params%dose_rate*current_time
            filter       = acc_dose2filter(dummy_img, acc_dose, params%kv, filtsz)
            write(logfhandle,'(a)') '>>> PRINTING DOSE WEIGHTS'
            do find=1,size(filter)
                write(logfhandle,'(A,1X,F8.2,1X,A,1X,F6.2)') '>>> RESOLUTION:', dummy_img%get_lp(find), 'WEIGHT: ', filter(find)
            end do
        end do
        call dummy_img%kill
        call simple_end('**** SIMPLE_PRINT_DOSE_WEIGHTS NORMAL STOP ****')
    end subroutine exec_print_dose_weights

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

    !> find the resolution of a Fourier index
    subroutine exec_res( self, cline )
        class(res_commander), intent(inout) :: self
        class(cmdline),       intent(inout) :: cline
        type(parameters) :: params
        real :: lp
        call params%new(cline)
        lp = (real(params%box-1)*params%smpd)/real(params%find)
        write(logfhandle,'(A,1X,f7.2)') '>>> LOW-PASS LIMIT:', lp
        call simple_end('**** SIMPLE_RES NORMAL STOP ****')
    end subroutine exec_res

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
        use simple_parameters, only: params_glob
        use simple_oris,       only: oris
        use simple_image,      only: image
        use simple_sym,        only: sym
        class(oris),   intent(inout) :: dsym_os
        class(image),  intent(inout) :: cylinder
        type(image)          :: read_img, img_msk, dist_img, roavg_img, topview
        type(sym)            :: se
        real,    allocatable :: corrs(:), forsort(:), e2(:), radii(:)
        integer, allocatable :: labels(:)
        logical, allocatable :: l_msk(:,:,:)
        integer  :: halfnoris, cnt1, cnt2, i, l, noris
        real     :: cen1, cen2, sum1, sum2, sumvals
        real     :: minmax(2), width, height, sh1(3), ang
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
            call read_img%noise_norm(l_msk)
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
            call read_img%bin_kmeans
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
                call read_img%noise_norm(l_msk)
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

end module simple_commander_misc
