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
public :: intgpeaks_commander
public :: masscen_commander
public :: print_dose_weights_commander
public :: print_fsc_commander
public :: print_magic_boxes_commander
public :: res_commander
public :: shift_commander
public :: stk_corr_commander
public :: kstest_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: cluster_smat_commander
  contains
    procedure :: execute      => exec_cluster_smat
end type cluster_smat_commander
type, extends(commander_base) :: intgpeaks_commander
  contains
    procedure :: execute       => exec_intgpeaks
end type intgpeaks_commander
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
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when reading: ', params%fname
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
            write(*,'(a,1x,f9.3,8x,a,1x,i3)') 'COHESION/SEPARATION RATIO INDEX: ', validinds(ncls), ' NCLS: ', ncls
            call build%spproj_field%read('clustering_shc_ncls'//int2str_pad(ncls,numlen)//trim(TXT_EXT), [1,build%spproj_field%get_noris()])
            do icls=1,ncls
                pop = build%spproj_field%get_pop(icls, params%label)
                write(*,'(a,3x,i5,1x,a,1x,i3)') '  CLUSTER POPULATION:', pop, 'CLUSTER:', icls
            end do
            write(*,'(a)') '***************************************************'
            if( ncls < params%ncls )then
                if( validinds(ncls+1) >= validinds(ncls) .and. .not. done )then
                    ncls_stop = ncls
                    done = .true.
                endif
            endif
        end do
        loc = minloc(validinds)
        ncls_min = loc(1)+1
        write(*,'(a,i3)') 'NUMBER OF CLUSTERS FOUND BY MINIMIZING RATIO INDEX: ', ncls_min
        if( ncls_stop /= 0 )then
            write(*,'(a,i3)') 'NUMBER OF CLUSTERS FOUND BY STOPPING CRITERIUM:     ', ncls_stop
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER_SMAT NORMAL STOP ****')
    end subroutine exec_cluster_smat

    subroutine exec_intgpeaks( self, cline )
        use simple_intg_atompeak
        use simple_atoms,  only: atoms
        class(intgpeaks_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters)      :: params
        type(builder)         :: build
        type(atoms)           :: mol
        real, allocatable     :: attribute(:)
        character(len=STDLEN) :: fbody, csv_name
        real                  :: xyz(3)
        integer               :: i, natoms, file_stat, fnr
        call build%init_params_and_build_general_tbox(cline,params)
        call build%vol%read(params%vols(1))
        if(cline%defined('msk').and.cline%defined('inner'))then
            ! peak finding
            call set_intgvol(build%vol, params%msk)
            call find_peaks(params%nptcls, params%inner, params%pdbfile)
        else
            ! peak integration
            call mol%new(params%pdbfile)
            natoms = mol%get_n()
            allocate(attribute(natoms), source=0.,stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("In: simple_commander_misc:: exec_intgpeaks attribute")
            call set_intgvol(build%vol)
            do i = 1, natoms
                xyz = mol%get_coord(i)
                if( cline%defined('inner') )then
                    attribute(i) = intg_shell(xyz, params%inner)
                else
                    attribute(i) = intg_nn(xyz)
                endif
            enddo
            ! output
            fbody = trim(get_fbody(trim(params%pdbfile), 'pdb'))
            csv_name = PATH_HERE//trim(adjustl(fbody))//'_intg.csv'
            call fopen(fnr, FILE=csv_name, STATUS='REPLACE', action='WRITE', iostat=file_stat)
            call fileiochk('commander_misc; exec_intgpeaks ', file_stat)
            do i = 1, natoms
                write(fnr,'(I6,A1,I6,A1,F12.6)')i, ',', mol%get_num(i), ',', attribute(i)
            enddo
            call fclose( fnr, errmsg='commander_misc; exec_intgpeaks ')
            deallocate(attribute)
        endif
        ! the end
        call simple_end('**** SIMPLE_INTGPEAKS NORMAL STOP ****')
    end subroutine exec_intgpeaks

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
        integer           :: iframe, find
        call params%new(cline)
        call dummy_img%new([params%box,params%box,1], params%smpd)
        time_per_frame = params%exp_time/real(params%nframes)
        do iframe=1,params%nframes
            current_time = real(iframe)*time_per_frame
            acc_dose     = params%dose_rate*current_time
            filter       = acc_dose2filter(dummy_img, acc_dose, params%kv)
            write(*,'(a)') '>>> PRINTING DOSE WEIGHTS'
            do find=1,size(filter)
                write(*,'(A,1X,F8.2,1X,A,1X,F6.2)') '>>> RESOLUTION:', dummy_img%get_lp(find), 'WEIGHT: ', filter(find)
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
        write(*,'(A,1X,F6.2,1X,A,1X,F15.3)') '>>> RESOLUTION:', res(k), '>>> FSC:', fsc(k)
        end do
        ! get & print resolution
        call get_resolution(fsc, res, res05, res0143)
        write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res0143
        write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res05
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
        write(*,'(A,1X,f7.2)') '>>> LOW-PASS LIMIT:', lp
        call simple_end('**** SIMPLE_RES NORMAL STOP ****')
    end subroutine exec_res

    !> for shifting a stack according to shifts in oritab
    subroutine exec_shift( self, cline )
        use simple_procimgfile, only: shift_imgfile
        class(shift_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
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
            write(*,'(I6,F8.3)')i,build%img%real_corr(build%img_copy, l_mask)
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
        write(*,'(a)') '>>> STATISTICS OF THE TWO DISTRIBUTIONS'
        call moment(dat1, ave1, sdev1, var, err)
        call moment(dat2, ave2, sdev2, var, err)
        write(*,'(a,1x,f4.2,1x,f4.2)') 'mean & sdev for infile : ', ave1, sdev1
        write(*,'(a,1x,f4.2,1x,f4.2)') 'mean & sdev for infile2: ', ave2, sdev2
        write(*,'(a)') '>>> KOLMOGOROV-SMIRNOV TEST TO DEDUCE EQUIVALENCE OR NON-EQUIVALENCE BETWEEN TWO DISTRIBUTIONS'
        call kstwo(dat1, ndat1, dat2, ndat2, ksstat, prob)
        write(*,'(a,1x,f4.2)') 'K-S statistic = ', ksstat
        write(*,'(a,1x,f4.2)') 'P             = ', prob
        write(*,'(a)') 'P represents the significance level for the null hypothesis that the two data sets are drawn from the same distribution'
        write(*,'(a)') 'Small P values show that the cumulative distribution functions of the two data sets differ significantly'
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

end module simple_commander_misc
