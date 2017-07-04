!> simple_commander_misc
!
!! This class contains the set of concrete miscellanous commanders of the SIMPLE
!! library. This class provides the glue between the reciver (main reciever is
!! simple_exec program) and the abstract action, which is simply execute (defined
!! by the base class: simple_commander_base). Later we can use the composite
!! pattern to create MacroCommanders (or workflows)
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
public :: masscen_commander
public :: print_cmd_dict_commander
public :: print_dose_weights_commander
public :: print_fsc_commander
public :: print_magic_boxes_commander
public :: res_commander
public :: shift_commander
public :: sym_aggregate_commander
public :: dsymsrch_commander
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
type, extends(commander_base) :: sym_aggregate_commander
  contains
    procedure :: execute       => exec_sym_aggregate
end type sym_aggregate_commander
type, extends(commander_base) :: dsymsrch_commander
  contains
    procedure :: execute       => exec_dsymsrch
end type dsymsrch_commander

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
        logical             :: done=.false.
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
    !> Calculate centre of mass 
    subroutine exec_masscen( self, cline )
        use simple_procimgfile, only: masscen_imgfile
        class(masscen_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params) :: p
        p = params(cline) ! parameters generated
        p%cenlp = p%lp
        ! center of mass centering
        if( cline%defined('thres') )then
            call masscen_imgfile(p%stk, p%outstk, p%smpd, p%lp, p%neg, p%msk, p%thres)
        else
            call masscen_imgfile(p%stk, p%outstk, p%smpd, p%lp, p%neg, p%msk)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_MASSCEN NORMAL STOP ****')
    end subroutine exec_masscen
    !> for printing the command line key dictonary
    subroutine exec_print_cmd_dict( self, cline )
        use simple_cmd_dict
        class(print_cmd_dict_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(params) :: p
        p = params(cline) ! parameters generated
        call print_cmd_key_descr( p%outfile )
        call simple_end('**** SIMPLE_PRINT_CMD_DICT NORMAL STOP ****')
    end subroutine exec_print_cmd_dict
    !> is a program for printing the dose weights applied to individual frames
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
    !> print_fsc  is a program for printing the binary FSC files produced by PRIME3D
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
    !> print_magic_boxesis a program for printing magic box sizes (fast FFT)
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
    !> shift is a program for shifting a stack according to shifts in oritab
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
    !> sym_aggregateis a program for robust identifiaction of the symmetry axis
    !> of a map using image-to-volume simiarity validation of the axis
    subroutine exec_sym_aggregate( self, cline )
        use simple_oris,  only: oris
        use simple_ori,   only: ori
        use simple_sym,   only: sym
        use simple_image, only: image
        use simple_projector_hlev, only: rotvol
        class(sym_aggregate_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(cmdline)      :: cline_c1
        type(params)       :: p, p_c1
        type(build)        :: b
        type(oris)         :: e, sympeaks, sym_axes
        type(ori)          :: symaxis, o
        type(image)        :: symvol, asym_vol
        type(sym)          :: se_c1
        real               :: cc, rotmat(3,3)
        integer            :: i, fnr, file_stat
        integer, parameter :: MAXLABELS = 10   ! maximum numbers symmetry peaks
        real,    parameter :: ANGTHRESH = 10.  ! maximum half-distance between symmetry peaks
        p = params(cline)                      ! parameters generated
        call b%build_general_tbox(p, cline)    ! general objects built
        ! init
        e = b%a ! b%a contains the orientations of the references projections
        cline_c1 = cline
        call cline_c1%set('pgrp', 'c1')
        p_c1 = params(cline_c1)
        call se_c1%new('c1')
        call asym_vol%new([p%box,p%box,p%box], p%smpd)
        call asym_vol%read(p%vols(1))
        ! spherical mask
        call b%mskvol%new([p%box,p%box,p%box], p%smpd)
        b%mskvol = 1.
        call b%mskvol%mask(p%msk, 'hard')
        ! identify top ranking symmetry peaks
        sym_axes = oris(nlines(p%oritab2))
        call sym_axes%read(p%oritab2)
        call find_sym_peaks(sym_axes, sympeaks)
        ! reconstruct & correlate volumes
        do i = 1, sympeaks%get_noris()
            symaxis = sympeaks%get_ori(i)
            b%a = e
            call b%a%rot(symaxis)
            ! symmetry
            call rec_vol(p, b%se)
            call symvol%copy( b%vol )
            ! call symvol%write('sym_vol'//int2str_pad(i,2)//p%ext)
            ! c1
            rotmat = symaxis%get_mat()
            call o%ori_from_rotmat(transpose(rotmat))
            call b%vol%copy( rotvol(asym_vol, o, p) )
            call b%vol%bp(p%hp, p%lp)
            call b%vol%mask(p%msk, 'soft')
            ! call b%vol%write('asym_vol'//int2str_pad(i,2)//p%ext)
            ! correlation
            cc = symvol%real_corr(b%vol, b%mskvol)
            call sympeaks%set(i, 'corr', cc)
        enddo
        ! output
        call sympeaks%write('sympeaks.txt')
        call sympeaks%write(p%outfile)
        ! the end
        fnr = get_fileunit()
        open(unit=fnr, FILE='SYM_AGGREGATE_FINISHED', STATUS='REPLACE', action='WRITE', iostat=file_stat)
        call fopen_err('In: commander_misc :: sym_aggregate', file_stat )
        close( unit=fnr )
        call simple_end('**** SIMPLE_SYM_AGGREGATE NORMAL STOP ****')

        contains

            subroutine rec_vol( p, se )
                type(params) :: p
                type(sym)    :: se
                call b%build_rec_tbox(p)
                call b%recvol%rec(p%stk, p, b%a, se, 1)
                call b%recvol%clip(b%vol)
                call b%vol%bp(p%hp, p%lp)
                call b%vol%mask(p%msk, 'soft')
            end subroutine rec_vol

            subroutine find_sym_peaks(sym_axes, sym_peaks)
                use simple_math, only: rad2deg
                class(oris), intent(inout) :: sym_axes
                class(oris), intent(out)   :: sym_peaks
                integer, allocatable :: sort_inds(:)
                type(oris) :: axes, tmp_axes
                type(ori)  :: axis, o
                integer    :: iaxis, naxes, loc(1), ilabel, axis_ind, istate, n_sympeaks
                axes       = sym_axes
                naxes      = axes%get_noris()
                tmp_axes   = oris(MAXLABELS)
                n_sympeaks = 0
                ! geometric clustering
                call axes%set_all2single('state', 0.)
                call axes%swape1e3
                sort_inds = axes%order_corr()
                do ilabel = 1, MAXLABELS
                    ! identify next best axis
                    do iaxis = 1, naxes
                        axis_ind = sort_inds(iaxis)
                        if( nint(axes%get(axis_ind,'state') ) == 0)exit
                    enddo
                    if(iaxis > naxes)exit
                    axis = axes%get_ori(axis_ind)
                    n_sympeaks = n_sympeaks + 1
                    call tmp_axes%set_ori(ilabel, axis)
                    ! flags axes within threshold
                    do iaxis = 1,naxes
                        axis_ind = sort_inds(iaxis)
                        o        = axes%get_ori(axis_ind)
                        istate   = nint(o%get('state'))
                        if(istate /= 0)cycle
                        if( rad2deg(axis.euldist.o) <= ANGTHRESH )then
                            axis_ind = sort_inds(iaxis)
                            call axes%set(axis_ind, 'state', real(ilabel))
                        endif
                    enddo
                enddo
                ! prep outcome
                write(*,'(A,I2)') '>>> SYMMETRY AXES PEAKS FOUND: ',n_sympeaks
                sym_peaks = oris(n_sympeaks)
                do iaxis=1,n_sympeaks
                    call sym_peaks%set_ori(iaxis, tmp_axes%get_ori(iaxis))
                enddo
                call sym_peaks%swape1e3
                deallocate(sort_inds)
            end subroutine find_sym_peaks

    end subroutine exec_sym_aggregate
    !> dsymsrch is a program for identifying rotational symmetries in class averages
    !> of D-symmetric molecules and generating a cylinder that matches the shape.
    subroutine exec_dsymsrch( self, cline )
        use simple_symsrcher, only: dsym_cylinder
        class(dsymsrch_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(build)  :: b
        type(params) :: p
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        call dsym_cylinder(p, b%a, b%vol)
        call b%a%write(p%outfile)
        call b%vol%write(p%outvol)
    end subroutine exec_dsymsrch

end module simple_commander_misc
