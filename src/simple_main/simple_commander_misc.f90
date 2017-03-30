!==Class simple_commander_misc
!
! This class contains the set of concrete miscellanous commanders of the SIMPLE library. This class provides the glue between the reciver 
! (main reciever is simple_exec program) and the abstract action, which is simply execute (defined by the base class: simple_commander_base). 
! Later we can use the composite pattern to create MacroCommanders (or workflows)
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Cyril Reboul & Hans Elmlund 2016
!
module simple_commander_misc
use simple_defs
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_strings,        only: int2str, int2str_pad
use simple_filehandling    ! use all in there
use simple_jiffys          ! use all in there
implicit none

public :: cluster_smat_commander
public :: find_nnimgs_commander
public :: masscen_commander
public :: print_cmd_dict_commander
public :: print_dose_weights_commander
public :: print_fsc_commander
public :: print_magic_boxes_commander
public :: res_commander
public :: shift_commander
private

type, extends(commander_base) :: cluster_smat_commander
  contains
    procedure :: execute      => exec_cluster_smat
end type cluster_smat_commander
type, extends(commander_base) :: find_nnimgs_commander
  contains
    procedure :: execute      => exec_find_nnimgs
end type find_nnimgs_commander
type, extends(commander_base) :: masscen_commander
  contains
    procedure :: execute      => exec_masscen
end type masscen_commander
type, extends(commander_base) :: print_cmd_dict_commander
  contains
    procedure :: execute       => exec_print_cmd_dict
end type print_cmd_dict_commander
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

contains

    subroutine exec_cluster_smat( self, cline )
        use simple_shc_cluster,   only: shc_cluster
        use simple_oris,          only: oris
        use simple_cluster_valid, only: cluster_valid
        class(cluster_smat_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)        :: p
        type(build)         :: b
        type(shc_cluster)   :: shcc
        type(cluster_valid) :: cvalid
        real, allocatable   :: smat(:,:)
        integer             :: funit, io_stat, ncls_min, loc(1), ncls_stop, icls
        integer             :: pop, ncls, alloc_stat, irestart, ntot, cnt=0, numlen
        real                :: avg_ratio, min_ratio, ratio, x, sim
        real, allocatable   :: validinds(:)
        integer, parameter  :: NRESTARTS=10
        logical             :: debug=.false., done=.false.
        p = params(cline,.false.)                        ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.)! general objects built
        ! obtain similarity matrix
        allocate(smat(p%nptcls,p%nptcls), stat=alloc_stat)
        call alloc_err('In: simple_cluster_smat, 1', alloc_stat)
        smat = 1.
        funit = get_fileunit()
        open(unit=funit, status='OLD', action='READ', file=p%fname, access='STREAM')
        read(unit=funit,pos=1,iostat=io_stat) smat
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when reading: ', p%fname
            stop 'I/O error; simple_cluster_smat'
        endif
        close(funit)
        allocate(validinds(2:p%ncls), stat=alloc_stat)
        call alloc_err("In: simple_cluster_smat", alloc_stat)
        validinds = 0
        ntot = (p%ncls-1)*NRESTARTS
        cnt = 0
        numlen = len(int2str(p%ncls))
        do ncls=2,p%ncls
            avg_ratio = 0.
            min_ratio = huge(x)
            do irestart=1,NRESTARTS
                cnt = cnt+1
                call progress(cnt,ntot)
                call shcc%new(p%nptcls, ncls, smat, b%a)
                call shcc%shc(.false., p%label, sim)
                call cvalid%new(b%a, ncls, p%label, smat)
                ratio = cvalid%ratio_index()
                avg_ratio = avg_ratio+ratio
                if( ratio < min_ratio )then
                    min_ratio = ratio
                    call b%a%write('shc_clustering_ncls'//int2str_pad(ncls,numlen)//'.txt')
                endif
            end do
            validinds(ncls) = avg_ratio/real(NRESTARTS)
        end do
        ncls_stop = 0
        done = .false. 
        do ncls=2,p%ncls
            write(*,'(a,1x,f9.3,8x,a,1x,i3)') 'COHESION/SEPARATION RATIO INDEX: ', validinds(ncls), ' NCLS: ', ncls
            call b%a%read('shc_clustering_ncls'//int2str_pad(ncls,numlen)//'.txt')
            do icls=1,ncls
                pop = b%a%get_pop(icls, p%label)
                write(*,'(a,3x,i5,1x,a,1x,i3)') '  CLUSTER POPULATION:', pop, 'CLUSTER:', icls
            end do
            write(*,'(a)') '***************************************************'
            if( ncls < p%ncls )then
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

    subroutine exec_find_nnimgs( self, cline )
        use simple_nnimgs      ! use all in there
        use simple_image,      only: image
        use simple_projector,  only: projector
        use simple_qsys_funs,  only: qsys_job_finished
        class(find_nnimgs_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)                  :: p
        integer,          allocatable :: nnmat(:,:), nnmat_part(:,:)
        type(projector),  allocatable :: imgs(:)
        character(len=:), allocatable :: fname
        integer :: iptcl, ineigh, alloc_stat, funit, io_stat, fnr
        p = params(cline, checkdistr=.false.) ! parameters generated, we don't split stack
        allocate( imgs(p%nptcls) )
        ! read, normalise & mask images
        do iptcl=1,p%nptcls
            call imgs(iptcl)%new([p%box,p%box,1],p%smpd)
            call imgs(iptcl)%read(p%stk, iptcl)
            call imgs(iptcl)%norm
            call imgs(iptcl)%mask(p%msk, 'soft')
        end do
        ! set Fourier index range
        p%kfromto(1) = max(2,imgs(1)%get_find(p%hp))
        p%kfromto(2) = imgs(1)%get_find(p%lp)
        call init_nn_srch( p, imgs )
        if( cline%defined('part') )then
            ! generate the partial nearest neighbor matrix
            allocate(nnmat_part(p%fromp:p%top,p%nnn), stat=alloc_stat)
            call alloc_err('In: simple_commander_misc :: exec_find_nnimgs, 1', alloc_stat)
            ! calculate the nearest neighbors
            call conduct_nn_srch([p%fromp,p%top], nnmat_part)
            ! write the nearest neighbors
            funit = get_fileunit()
            allocate(fname, source='nnmat_part'//int2str_pad(p%part,p%numlen)//'.bin')
            open(unit=funit, status='REPLACE', action='WRITE', file=fname, access='STREAM')
            write(unit=funit,pos=1,iostat=io_stat) nnmat_part(:,:)
            ! check if the write was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(commander_misc :: exec_find_nnimgs): I/O error ',&
                io_stat, ' when writing to: ', fname
                stop
            endif
            close(funit)
            deallocate(fname, nnmat_part)
            call qsys_job_finished( p, 'commander_misc :: exec_find_nnimg' )
        else
            nnmat = img_nearest_neighbors( p, imgs )
            do iptcl=1,p%nptcls
                allocate(fname, source='nnimgs'//int2str(iptcl)//p%ext)
                do ineigh=1,p%nnn
                    call imgs(nnmat(iptcl,ineigh))%bwd_ft
                    call imgs(nnmat(iptcl,ineigh))%write(fname,ineigh)
                end do
                deallocate(fname)
            end do
            do iptcl=1,p%nptcls
                call imgs(iptcl)%kill
            end do
            deallocate(nnmat)
        endif
        deallocate(imgs)
        ! end gracefully
        call simple_end('**** SIMPLE_FIND_NNIMGS NORMAL STOP ****')
    end subroutine exec_find_nnimgs

    subroutine exec_masscen( self, cline )
        use simple_procimgfile, only: masscen_imgfile
        class(masscen_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params) :: p
        p = params(cline) ! parameters generated
        ! center of mass centering
        if( cline%defined('thres') )then
            call masscen_imgfile(p%stk, p%outstk, p%smpd, p%lp, p%neg, p%msk, p%thres)
        else
            call masscen_imgfile(p%stk, p%outstk, p%smpd, p%lp, p%neg, p%msk)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_MASSCEN NORMAL STOP ****')
    end subroutine exec_masscen

    subroutine exec_print_cmd_dict( self, cline )
        use simple_cmd_dict
        class(print_cmd_dict_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(params) :: p
        p = params(cline) ! parameters generated
        call print_cmd_key_descr( p%outfile )
        call simple_end('**** SIMPLE_PRINT_CMD_DICT NORMAL STOP ****')
    end subroutine exec_print_cmd_dict

    subroutine exec_print_dose_weights( self, cline )
        use simple_image,    only: image
        use simple_filterer, only: acc_dose2filter
        class(print_dose_weights_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        real, allocatable :: filter(:)
        type(image)       :: dummy_img
        real              :: time_per_frame, current_time, acc_dose
        integer           :: iframe, alloc_stat, find
        type(params)      :: p
        p = params(cline) ! parameters generated
        call dummy_img%new([p%box,p%box,1], p%smpd)
        time_per_frame = p%exp_time/real(p%nframes)
        do iframe=1,p%nframes
            current_time = real(iframe)*time_per_frame
            acc_dose     = p%dose_rate*current_time
            filter       = acc_dose2filter(dummy_img, acc_dose, p%kv)
            write(*,'(a)') '>>> PRINTING DOSE WEIGHTS'
            do find=1,size(filter)
                write(*,'(A,1X,F8.2,1X,A,1X,F6.2)') '>>> RESOLUTION:', dummy_img%get_lp(find), 'WEIGHT: ', filter(find)
            end do
        end do
        call dummy_img%kill
        call simple_end('**** SIMPLE_PRINT_DOSE_WEIGHTS NORMAL STOP ****')
    end subroutine exec_print_dose_weights

    subroutine exec_print_fsc( self, cline )
        use simple_math,  only: get_resolution, get_lplim
        use simple_image, only: image
        class(print_fsc_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(params)      :: p
        type(image)       :: img
        real, allocatable :: res(:), fsc(:)
        integer           :: k
        real              :: res0143, res05
        p = params(cline) ! parameters generated
        call img%new([p%box,p%box,1], p%smpd)
        res = img%get_res() 
        fsc = file2rarr(p%fsc)
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

    subroutine exec_print_magic_boxes( self, cline )
        use simple_magic_boxes, only: print_magic_box_range
        class(print_magic_boxes_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(params) :: p
        p = params(cline) ! parameters generated
        call print_magic_box_range(p%smpd, p%moldiam )
        ! end gracefully
        call simple_end('**** SIMPLE_PRINT_MAGIC_BOXES NORMAL STOP ****')
    end subroutine exec_print_magic_boxes
    
    subroutine exec_res( self, cline )
        class(res_commander), intent(inout) :: self
        class(cmdline),       intent(inout) :: cline
        type(params) :: p
        real         :: lp
        p  = params(cline)
        lp = (real(p%box-1)*p%smpd)/real(p%find)
        write(*,'(A,1X,f7.2)') '>>> LOW-PASS LIMIT:', lp
        call simple_end('**** SIMPLE_RES NORMAL STOP ****')
    end subroutine exec_res

    subroutine exec_shift( self, cline )
        use simple_procimgfile, only: shift_imgfile
        class(shift_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        p = params(cline)                                 ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        call shift_imgfile(p%stk, p%outstk, b%a, p%smpd, p%mul)
        call b%a%zero_shifts
        call b%a%write('shiftdoc.txt')
        call simple_end('**** SIMPLE_SHIFT NORMAL STOP ****')
    end subroutine exec_shift

end module simple_commander_misc
