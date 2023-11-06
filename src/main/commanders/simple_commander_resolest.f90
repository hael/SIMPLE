! concrete commander: resolest for resolution estimation
module simple_commander_resolest
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,     only: parameters, params_glob
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_image,          only: image
use simple_masker,         only: masker
use simple_fsc
implicit none

public :: fsc_commander
public :: whiten_and_filter_commander
public :: nununiform_filter3D_commander
public :: nununiform_filter2D_commander
public :: uniform_filter2D_commander
public :: uniform_filter3D_commander
public :: cavg_filter2D_commander
public :: prune_cavgs_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: fsc_commander
  contains
    procedure :: execute      => exec_fsc
end type fsc_commander

type, extends(commander_base) :: whiten_and_filter_commander
  contains
    procedure :: execute      => exec_whiten_and_filter
end type whiten_and_filter_commander

type, extends(commander_base) :: nununiform_filter2D_commander
  contains
    procedure :: execute      => exec_nununiform_filter2D
end type nununiform_filter2D_commander

type, extends(commander_base) :: uniform_filter2D_commander
  contains
    procedure :: execute      => exec_uniform_filter2D
end type uniform_filter2D_commander

type, extends(commander_base) :: nununiform_filter3D_commander
  contains
    procedure :: execute      => exec_nununiform_filter3D
end type nununiform_filter3D_commander

type, extends(commander_base) :: uniform_filter3D_commander
  contains
    procedure :: execute      => exec_uniform_filter3D
end type uniform_filter3D_commander

type, extends(commander_base) :: cavg_filter2D_commander
  contains
    procedure :: execute      => exec_cavg_filter2D
end type cavg_filter2D_commander

type, extends(commander_base) :: prune_cavgs_commander
  contains
    procedure :: execute      => exec_prune_cavgs
end type prune_cavgs_commander

contains

    !> calculates Fourier shell correlation from Even/Odd Volume pairs
    subroutine exec_fsc( self, cline )
        class(fsc_commander), intent(inout) :: self
        class(cmdline),       intent(inout) :: cline
        type(parameters)               :: params
        type(image)                    :: even, odd
        type(masker)                   :: mskvol
        character(len=LONGSTRLEN) :: fsc_templ
        real,              allocatable :: res(:), fsc(:), fsc_t(:), fsc_n(:)
        integer :: j, find_plate, k_hp, k_lp, nyq
        real    :: res_fsc05, res_fsc0143
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        ! read even/odd pair
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd%new([params%box,params%box,params%box], params%smpd)
        call odd%read(params%vols(1))
        call even%read(params%vols(2))
        res = even%get_res()
        nyq = even%get_filtsz()
        if( params%l_filemsk )then
            call mskvol%new([params%box,params%box,params%box], params%smpd)
            call mskvol%read(params%mskfile)
            call phase_rand_fsc(even, odd, mskvol, params%msk, 1, nyq, fsc, fsc_t, fsc_n)
        else
            ! spherical masking
            call even%mask(params%msk, 'soft')
            call odd%mask(params%msk, 'soft')
            call even%fft()
            call odd%fft()
            allocate(fsc(nyq),source=0.)
            call even%fsc(odd, fsc)
        endif
        if( params%l_phaseplate ) call phaseplate_correct_fsc(fsc, find_plate)
        do j=1,nyq
            write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', fsc(j)
        end do
        call get_resolution(fsc, res, res_fsc05, res_fsc0143)
        k_hp = calc_fourier_index(params%hp, params%box, params%smpd)
        k_lp = calc_fourier_index(params%lp, params%box, params%smpd)
        write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
        write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
        call even%kill
        call odd%kill
        fsc_templ = trim(FSC_FBODY)//int2str_pad(1,2)
        call arr2file(fsc, trim(fsc_templ)//trim(BIN_EXT))
        if( params%l_filemsk )then
            call plot_phrand_fsc(size(fsc), fsc, fsc_t, fsc_n, res, params%smpd, fsc_templ)
        else
            call plot_fsc(size(fsc), fsc, res, params%smpd, fsc_templ)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_FSC NORMAL STOP ****')
    end subroutine exec_fsc

    subroutine exec_whiten_and_filter( self, cline )
        class(whiten_and_filter_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters)  :: params
        type(image)       :: even, odd
        real, allocatable :: res(:), fsc(:), filter(:)
        integer :: j, nyq
        real    :: res_fsc05, res_fsc0143
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        ! read even/odd pair
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd%new([params%box,params%box,params%box], params%smpd)
        call odd%read(params%vols(1))
        call even%read(params%vols(2))
        res = even%get_res()
        nyq = even%get_filtsz()
        ! spherical masking
        call even%mask(params%msk, 'soft')
        call odd%mask(params%msk, 'soft')
        call even%fft()
        call odd%fft()
        allocate(fsc(nyq),filter(nyq),source=0.)
        call even%fsc(odd, fsc)
        do j=1,nyq
            write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', fsc(j)
        end do
        call get_resolution(fsc, res, res_fsc05, res_fsc0143)
        write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
        write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
        ! whiten and apply filter
        call fsc2optlp_sub(nyq, fsc, filter)
        call even%whiten_noise_power(odd, is_ptcl=.false.)
        call even%ifft
        call odd%ifft
        call even%write('even_whitened.mrc')
        call odd%write('odd_whitened.mrc')
        call even%apply_filter(filter)
        call odd%apply_filter(filter)
        call even%write('even_whitened_filtered.mrc')
        call odd%write('odd_whitened_filter.mrc')
        call even%kill
        call odd%kill
        ! end gracefully
        call simple_end('**** SIMPLE_WHITEN_AND_FILTER NORMAL STOP ****')
    end subroutine exec_whiten_and_filter

    subroutine exec_nununiform_filter3D(self, cline)
        use simple_opt_filter, only: nonuni_filt3D
        class(nununiform_filter3D_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters)  :: params
        type(image)       :: even, odd, mskvol
        logical           :: have_mask_file
        character(len=90) :: file_tag
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call odd %new([params%box,params%box,params%box], params%smpd)
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd %read(params%vols(1))
        call even%read(params%vols(2))
        file_tag = 'nonuniform_3D_filter_ext_'//int2str(params%smooth_ext)
        have_mask_file = .false.
        if( params%l_filemsk )then
            call mskvol%new([params%box,params%box,params%box], params%smpd)
            call mskvol%read(params%mskfile)
            call even%zero_background
            call odd%zero_background
            call even%mul(mskvol)
            call odd%mul(mskvol)
            call mskvol%one_at_edge ! to expand before masking of reference internally
            have_mask_file = .true.
        else
            ! spherical masking
            call even%mask(params%msk, 'soft')
            call odd%mask(params%msk, 'soft')
            call mskvol%disc([params%box,params%box,params%box], params%smpd,&
            &real(min(params%box/2, int(params%msk + COSMSKHALFWIDTH))))
        endif
        if( cline%defined('lpstop') )then
            call nonuni_filt3D(odd, even, mskvol, params%lpstop)
        else
            call nonuni_filt3D(odd, even, mskvol)
        endif
        if( have_mask_file )then
            call mskvol%read(params%mskfile) ! restore the soft edge
            call even%mul(mskvol)
            call odd%mul(mskvol)
        else
            call even%mask(params%msk, 'soft')
            call odd%mask(params%msk, 'soft')
        endif
        call odd%write(trim(file_tag)//'_odd.mrc')
        call even%write(trim(file_tag)//'_even.mrc')
        call odd%add(even)
        call odd%mul(0.5)
        call odd%write(trim(file_tag)//'_avg.mrc')
        ! destruct
        call odd%kill
        call even%kill
        call mskvol%kill
        ! end gracefully
        call simple_end('**** SIMPLE_nununiform_filter3D NORMAL STOP ****')
    end subroutine exec_nununiform_filter3D

    subroutine exec_uniform_filter3D(self, cline)
        use simple_opt_filter, only: uni_filt3D
        class(uniform_filter3D_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters)  :: params
        type(image)       :: even, odd, mskvol
        logical           :: have_mask_file
        character(len=90) :: file_tag
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call odd %new([params%box,params%box,params%box], params%smpd)
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd %read(params%vols(1))
        call even%read(params%vols(2))
        file_tag = 'uniform_3D_filter'
        have_mask_file = .false.
        if( params%l_filemsk )then
            call mskvol%new([params%box,params%box,params%box], params%smpd)
            call mskvol%read(params%mskfile)
            call even%zero_background
            call  odd%zero_background
            call even%mul(mskvol)
            call  odd%mul(mskvol)
            call mskvol%one_at_edge ! to expand before masking of reference internally
            have_mask_file = .true.
        else
            ! spherical masking
            call even%mask(params%msk, 'soft')
            call  odd%mask(params%msk, 'soft')
            call mskvol%disc([params%box,params%box,params%box], params%smpd,&
                    &real(min(params%box/2, int(params%msk + COSMSKHALFWIDTH))))
        endif
        call uni_filt3D(odd, even, mskvol)
        if( have_mask_file )then
            call mskvol%read(params%mskfile) ! restore the soft edge
            call even%mul(mskvol)
            call  odd%mul(mskvol)
        else
            call even%mask(params%msk, 'soft')
            call  odd%mask(params%msk, 'soft')
        endif
        call  odd%write(trim(file_tag)//'_odd.mrc')
        call even%write(trim(file_tag)//'_even.mrc')
        call odd%add(even)
        call odd%mul(0.5)
        call odd%write(trim(file_tag)//'_avg.mrc')
        ! destruct
        call    odd%kill
        call   even%kill
        call mskvol%kill
        ! end gracefully
        call simple_end('**** SIMPLE_uniform_filter3D NORMAL STOP ****')
    end subroutine exec_uniform_filter3D

    subroutine exec_nununiform_filter2D( self, cline )
        use simple_opt_filter, only: nonuni_filt2D_sub
        use simple_masker,     only: automask2D
        use simple_default_clines
        class(nununiform_filter2D_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        character(len=:), allocatable :: file_tag
        type(image),      allocatable :: even(:), odd(:), mask(:)
        real,             allocatable :: diams(:)
        type(parameters) :: params
        integer          :: iptcl
        ! init
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call set_automask2D_defaults(cline)
        call params%new(cline)
        call find_ldim_nptcls(params%stk, params%ldim, params%nptcls)
        params%ldim(3) = 1 ! because we operate on stacks
        file_tag = 'nonuniform_filter2D_ext_'//int2str(params%smooth_ext)
        ! allocate
        allocate(odd(params%nptcls), even(params%nptcls), mask(params%nptcls))
        ! construct & read
        do iptcl = 1, params%nptcls
            call odd( iptcl)%new(params%ldim, params%smpd, .false.)
            call even(iptcl)%new(params%ldim, params%smpd, .false.)
            call odd( iptcl)%read(params%stk,  iptcl)
            call even(iptcl)%read(params%stk2, iptcl)
            call mask(iptcl)%copy(odd(iptcl))
            call mask(iptcl)%add(even(iptcl))
            call mask(iptcl)%mul(0.5)
        enddo
        ! filter
        if( params%l_automsk )then
            call automask2D(mask, params%ngrow, nint(params%winsz), params%edge, diams)
            call nonuni_filt2D_sub(even, odd, mask)
        else
            call nonuni_filt2D_sub(even, odd)
        endif
        ! write output and destruct
        do iptcl = 1, params%nptcls
            call odd( iptcl)%write(trim(file_tag)//'_odd.mrc',  iptcl)
            call even(iptcl)%write(trim(file_tag)//'_even.mrc', iptcl)
            call odd(iptcl)%add(even(iptcl))
            call odd(iptcl)%mul(0.5)
            call odd(iptcl)%write(trim(file_tag)//'_avg.mrc', iptcl)
            call odd( iptcl)%kill()
            call even(iptcl)%kill()
        end do
        if( allocated(diams) ) deallocate(diams)
        ! end gracefully
        call simple_end('**** SIMPLE_nununiform_filter2D NORMAL STOP ****')
    end subroutine exec_nununiform_filter2D

    subroutine exec_uniform_filter2D( self, cline )
        use simple_opt_filter, only: uni_filt2D_sub
        use simple_masker,     only: automask2D
        use simple_image,      only: image_ptr 
        use simple_default_clines
        class(uniform_filter2D_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        character(len=:), allocatable :: file_tag
        type(image),      allocatable :: odd(:), even(:), mask(:)
        real,             allocatable :: diams(:)
        type(parameters) :: params
        type(image_ptr)  :: pmask
        integer          :: iptcl
        ! init
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call set_automask2D_defaults(cline)
        call params%new(cline)
        call find_ldim_nptcls(params%stk, params%ldim, params%nptcls)
        params%ldim(3) = 1 ! because we operate on stacks
        file_tag = 'uniform_filter2D'
        ! allocate
        allocate(odd(params%nptcls), even(params%nptcls), mask(params%nptcls))
        ! construct & read
        do iptcl = 1, params%nptcls
            call odd( iptcl)%new( params%ldim, params%smpd, .false.)
            call odd( iptcl)%read(params%stk,  iptcl)
            call even(iptcl)%new( params%ldim, params%smpd, .false.)
            call even(iptcl)%read(params%stk2, iptcl)
            call mask(iptcl)%copy(odd(iptcl))
            call mask(iptcl)%add(even(iptcl))
            call mask(iptcl)%mul(0.5)
        enddo
        ! filter
        if( params%l_automsk )then
            call automask2D(mask, params%ngrow, nint(params%winsz), params%edge, diams)
            do iptcl = 1, params%nptcls
                call mask(iptcl)%get_mat_ptrs(pmask)
            enddo
        else
            do iptcl = 1, params%nptcls
                call mask(iptcl)%get_mat_ptrs(pmask)
                pmask%rmat = 1;
            enddo
        endif
        call uni_filt2D_sub(even, odd, mask)
        ! write output and destruct
        do iptcl = 1, params%nptcls
            call odd( iptcl)%write(trim(file_tag)//'_odd.mrc',  iptcl)
            call even(iptcl)%write(trim(file_tag)//'_even.mrc',  iptcl)
            call odd( iptcl)%kill()
            call even(iptcl)%kill()
        end do
        if( allocated(diams) ) deallocate(diams)
        ! end gracefully
        call simple_end('**** SIMPLE_uniform_filter2D NORMAL STOP ****')
    end subroutine exec_uniform_filter2D

    subroutine exec_cavg_filter2D( self, cline )
        use simple_strategy2D3D_common, only: read_imgbatch, prepimgbatch, prepimg4align
        use simple_polarft_corrcalc,    only: polarft_corrcalc
        use simple_regularizer,         only: regularizer
        class(cavg_filter2D_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        complex,          allocatable :: cmat(:,:)
        integer,          allocatable :: pinds(:)
        character(len=:), allocatable :: cavgsstk
        complex(dp),      allocatable :: cls_avg(:,:), ptcl_rot(:,:)
        real(dp),         allocatable :: ctf_rot(:,:), denom(:,:)
        type(polarft_corrcalc)        :: pftcc
        type(builder)                 :: build
        type(parameters)              :: params
        type(image)                   :: img_cavg, calc_cavg, ptcl_match
        type(regularizer)             :: reg_obj
        integer  :: nptcls, iptcl, nptcls_cls
        logical  :: l_ctf
        integer  :: ncls, j, box, loc
        real     :: smpd
        call cline%set('dir_exec', 'cavg_filter2D')
        call cline%set('mkdir',    'yes')
        call cline%set('oritype',  'ptcl2D')
        call build%init_params_and_build_general_tbox(cline,params)
        call build%spproj%update_projinfo(cline)
        ! reading all from the class 'class'
        call build%spproj%os_ptcl2D%get_pinds(params%class, 'class', pinds)
        nptcls     = build%spproj%get_nptcls()
        nptcls_cls = size(pinds)
        call pftcc%new(nptcls_cls, [1,nptcls_cls], params%kfromto)
        call pftcc%reallocate_ptcls(nptcls_cls, pinds)
        call reg_obj%new(pftcc)
        call build%img_crop_polarizer%init_polarizer(pftcc, params_glob%alpha)
        call ptcl_match%new([params%box_crop, params%box_crop, 1], params%smpd_crop)
        call prepimgbatch(nptcls)
        call read_imgbatch([1, nptcls])
        ! getting the ctfs
        l_ctf = build%spproj%get_ctfflag('ptcl2D',iptcl=pinds(1)).ne.'no'
        ! make CTFs
        if( l_ctf ) call pftcc%create_polar_absctfmats(build%spproj, 'ptcl2D')
        ! computing the class average (mimicking reg's cavg) and comparing to the cluster2D_cavg
        allocate(ctf_rot(pftcc%pftsz, pftcc%kfromto(1):pftcc%kfromto(2)),&
               &ptcl_rot(pftcc%pftsz, pftcc%kfromto(1):pftcc%kfromto(2)),&
                &cls_avg(pftcc%pftsz, pftcc%kfromto(1):pftcc%kfromto(2)),&
                  &denom(pftcc%pftsz, pftcc%kfromto(1):pftcc%kfromto(2)))
        cls_avg = 0.
        denom   = 0.
        do j = 1, nptcls_cls
            iptcl = pinds(j)
            ! prep
            call prepimg4align(iptcl, build%imgbatch(iptcl), ptcl_match)
            ! transfer to polar coordinates
            call build%img_crop_polarizer%polarize(pftcc, ptcl_match, iptcl, .true., .true., mask=build%l_resmsk)
            ! e/o flags
            call pftcc%set_eo(iptcl, .true. )
            ! accumulating the cls_avg
            loc = pftcc%get_roind(build%spproj_field%e3get(iptcl))
            if( loc > pftcc%nrots ) loc = loc - pftcc%nrots
            call reg_obj%rotate_polar(cmplx(pftcc%pfts_ptcls(:,:,j), kind=dp), ptcl_rot, loc)
            call reg_obj%rotate_polar(pftcc%ctfmats(:,:,j), ctf_rot, loc)
            cls_avg = cls_avg + ptcl_rot*ctf_rot
            denom   = denom   +          ctf_rot**2
            ! writing the raw stack
            call pftcc%polar2cartesian(pftcc%pfts_ptcls(:,:,j), cmat, box)
            call calc_cavg%new([box,box,1], params%smpd*real(params%box)/real(box))
            call calc_cavg%zero_and_flag_ft
            call calc_cavg%set_cmat(cmat)
            call calc_cavg%shift_phorig()
            call calc_cavg%ifft
            call calc_cavg%write('ptcls_stk.mrc', j)
            ! writing the ctf stack
            call pftcc%polar2cartesian(cmplx(pftcc%ctfmats(:,:,j), kind=sp), cmat, box)
            call calc_cavg%new([box,box,1], params%smpd*real(params%box)/real(box))
            call calc_cavg%zero_and_flag_ft
            call calc_cavg%set_cmat(cmat)
            call calc_cavg%shift_phorig()
            call calc_cavg%ifft
            call calc_cavg%write('ctfs_stk.mrc', j)
            ! writing the aligned ptcls stack
            call pftcc%polar2cartesian(cmplx(ptcl_rot, kind=sp), cmat, box)
            call calc_cavg%new([box,box,1], params%smpd*real(params%box)/real(box))
            call calc_cavg%zero_and_flag_ft
            call calc_cavg%set_cmat(cmat)
            call calc_cavg%shift_phorig()
            call calc_cavg%ifft
            call calc_cavg%write('aligned_ptcls_stk.mrc', j)
            ! writing the aligned ctf stack
            call pftcc%polar2cartesian(cmplx(ctf_rot, kind=sp), cmat, box)
            call calc_cavg%new([box,box,1], params%smpd*real(params%box)/real(box))
            call calc_cavg%zero_and_flag_ft
            call calc_cavg%set_cmat(cmat)
            call calc_cavg%shift_phorig()
            call calc_cavg%ifft
            call calc_cavg%write('aligned_ctfs_stk.mrc', j)
        enddo
        ! polar class average
        call pftcc%polar2cartesian(cmplx(cls_avg / denom, kind=sp), cmat, box)
        call calc_cavg%new([box,box,1], params%smpd*real(params%box)/real(box))
        call calc_cavg%zero_and_flag_ft
        call calc_cavg%set_cmat(cmat)
        call calc_cavg%shift_phorig()
        call calc_cavg%ifft
        call calc_cavg%write('polar_cavg.mrc')
        ! writing the cluster2D_cavg of the current class
        call build%spproj%get_cavgs_stk(cavgsstk, ncls, smpd)
        call img_cavg%new([params%box,params%box,1], params%smpd)
        call img_cavg%read(cavgsstk, params%class)
        call img_cavg%write('cluster2D_cavg.mrc')
        call img_cavg%kill
        ! end gracefully
        call simple_end('**** SIMPLE_cavg_filter2D NORMAL STOP ****')
    end subroutine exec_cavg_filter2D

    subroutine exec_prune_cavgs( self, cline )
        use simple_strategy2D3D_common, only: read_imgbatch, prepimgbatch, prepimg4align
        use simple_polarft_corrcalc,    only: polarft_corrcalc
        use simple_pftcc_shsrch_grad,   only: pftcc_shsrch_grad
        use simple_ctf,                 only: ctf
        class(prune_cavgs_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        integer, parameter :: NITERS         = 3
        integer, parameter :: NPTCLS_PER_BIN = 50
        type(pftcc_shsrch_grad), allocatable :: grad_shsrch_objs(:)
        type(polarft_corrcalc)        :: pftcc
        type(builder)                 :: build
        type(parameters)              :: params
        type(image),      allocatable :: tmp_imgs(:)
        complex(dp),      allocatable :: cls_avg_bak(:,:),cls_avg(:,:), cls_avg_even(:,:), cls_avg_odd(:,:), num_even(:,:), num_odd(:,:)
        complex(sp),      allocatable :: ptcl(:,:), ptcl_rot(:,:), diff(:,:)
        real(dp),         allocatable :: denom_even(:,:), denom_odd(:,:), R2s(:), RmI2s(:), weights(:)
        real(sp), target, allocatable :: sig2(:,:)
        real(sp),         allocatable :: purity(:),inpl_corrs(:), corrs(:), ctf_rot(:,:), shifts(:,:), dfs(:), bindiff(:)
        real(sp),         allocatable :: frc(:), sig2_even(:,:), sig2_odd(:,:), tmp(:), means(:), sdevs(:), binccs(:)
        integer,          allocatable :: rots(:),pinds(:), states(:), order(:), labels(:), batches(:,:)
        integer,          allocatable :: bin_inds(:,:), bins(:), cls2batch(:), cls_pops(:)
        logical,          allocatable :: selected(:)
        character(len=:), allocatable :: cavgs_stk
        type(ctf)        :: tfun
        type(ctfparams)  :: ctfparms
        real(dp) :: rmi2, r2
        real     :: cxy(3), lims(2,2), lims_init(2,2), threshold, cc, pu, ice_score, score, mean, sdev
        real     :: df, cc_df_corrs, prev_threshold, mdf,mcorr,sdf,scorr, cavgs_smpd, sdev_noise
        integer  :: nstks, nptcls, iptcl, iter, n_lines, icls, nbins, batch_start, batch_end
        integer  :: ibatch, batchsz, ibin, binpop, ithr, nsel, prev_nsel, ini_pop, cavgs_ncls
        integer  :: ncls, i, j, k,fnr, irot, nptcls_sel, pop, nbatches, batchsz_max, pop_sel
        logical  :: l_ctf, l_groundtruth, l_corr_ranking, l_weighted_init, l_write
        call cline%set('oritype', 'ptcl2D')
        call cline%set('mkdir',   'yes')
        if( .not.cline%defined('objfun') )     call cline%set('objfun',  'cc')
        if( .not.cline%defined('reject_cls') ) call cline%set('reject_cls',  'no')
        call build%init_params_and_build_general_tbox(cline, params)
        call build%spproj%get_cavgs_stk(cavgs_stk, cavgs_ncls, cavgs_smpd)
        ! class-based fracture identification
        ! do icls = 1,cavgs_ncls
        !     call build%img%read(cavgs_stk,icls)
        !     call build%img%mask(real(params%box)/2.-COSMSKHALFWIDTH,'soft')
        !     call build%img%fft
        !     call build%img%calc_line_score(ice_score)
        !     print *,icls, ice_score
        ! enddo
        ncls       = build%spproj_field%get_n('class')
        states     = nint( build%spproj_field%get_all('state'))
        nptcls     = size(states)
        nptcls_sel = count(states==1)
        l_groundtruth   = cline%defined('infile')
        l_corr_ranking  = .true.
        l_weighted_init = .true.
        l_write         = .true.
        call build%spproj_field%get_pops(cls_pops, 'class', maxn=ncls)
        ! class based outlier detection
        call detect_class_outliers(cavgs_stk)
        if( l_groundtruth )then
            n_lines = nlines(trim(params%infile))
            allocate(labels(n_lines))
            call fopen(fnr, FILE=trim(params%infile), STATUS='OLD', action='READ')
            do i=1,n_lines
                read(fnr,*) labels(i)
            end do
            call fclose(fnr)
            states = nint(build%spproj_field%get_all('state'))
            call build%spproj_field%set_all('state', real(labels))
            call build%spproj%write('cleaned.simple')
            call build%spproj_field%set_all('state', real(states))
            deallocate(states)
        endif
        nstks = build%spproj%get_nstks()
        params%l_kweight_rot   = .false.
        params%l_kweight_shift = .false.
        allocate(grad_shsrch_objs(params%nthr),tmp_imgs(params%nthr))
        lims(:,1)       = -MINSHIFT
        lims(:,2)       =  MINSHIFT
        lims_init(:,1)  = -MINSHIFT/2.
        lims_init(:,2)  =  MINSHIFT/2.
        do ithr = 1,nthr_glob
            call tmp_imgs(ithr)%new([params%box_crop,params%box_crop,1],params%smpd_crop,wthreads=.false.)
        enddo
        ! Allocations
        call pftcc%new(NITERS, [1,1], params%kfromto)
        allocate(ptcl(pftcc%pftsz,params%kfromto(1):params%kfromto(2)),&
        &ptcl_rot(pftcc%pftsz,params%kfromto(1):params%kfromto(2)),&
        &ctf_rot(pftcc%pftsz,params%kfromto(1):params%kfromto(2)),&
        &cls_avg(pftcc%pftsz,params%kfromto(1):params%kfromto(2)),&
        &cls_avg_even(pftcc%pftsz,params%kfromto(1):params%kfromto(2)),&
        &cls_avg_odd(pftcc%pftsz,params%kfromto(1):params%kfromto(2)),&
        &num_even(pftcc%pftsz,params%kfromto(1):params%kfromto(2)),&
        &num_odd(pftcc%pftsz,params%kfromto(1):params%kfromto(2)),&
        &denom_even(pftcc%pftsz,params%kfromto(1):params%kfromto(2)),&
        &denom_odd(pftcc%pftsz,params%kfromto(1):params%kfromto(2)),&
        &diff(pftcc%pftsz,params%kfromto(1):params%kfromto(2)),&
        &inpl_corrs(pftcc%nrots),frc(params%kfromto(1):params%kfromto(2)),&
        &purity(0:NITERS))
        call build%img_crop_polarizer%init_polarizer(pftcc, params_glob%alpha)
        if( (params%cc_objfun==OBJFUN_EUCLID) )then
            allocate(sig2(params%kfromto(1):params%kfromto(2),nptcls),&
            &sig2_even(params%kfromto(1):params%kfromto(2),nstks),&
            &sig2_odd(params%kfromto(1):params%kfromto(2),nstks))
        endif
        do ithr = 1, params%nthr
            call grad_shsrch_objs(ithr)%new(lims, lims_init=lims_init,&
            &shbarrier=params%shbarrier, maxits=60, opt_angle=.true.)
        enddo
        ! Class loop
        states = nint(build%spproj_field%get_all('state'))
        do icls = 1,ncls
            ! indices
            call build%spproj_field%get_pinds(icls, 'class', pinds)
            pop     = size(pinds)
            ini_pop = pop
            nbins   = ceiling(real(pop)/real(NPTCLS_PER_BIN))
            if( nbins < 1 ) cycle
            ! images
            batchsz_max = pop   ! batchsz_max = min(pop,params_glob%nthr*BATCHTHRSZ)
            nbatches    = ceiling(real(pop)/real(batchsz_max))
            batches     = split_nobjs_even(pop, nbatches)
            batchsz_max = maxval(batches(:,2)-batches(:,1)+1)
            if( allocated(build%imgbatch) )then
                if( batchsz_max > size(build%imgbatch) ) call prepimgbatch(batchsz_max)
            else
                call prepimgbatch(batchsz_max)
            endif
            if( l_groundtruth )then
                ! Initial purity
                purity(0) = real(sum(labels(pinds(:))))
                purity(0) = purity(0) * 100. / real(pop)
                print *, icls,pop,purity(0)
            endif
            allocate(means(ini_pop),sdevs(ini_pop))
            do ibatch=1,nbatches
                ! read
                batch_start = batches(ibatch,1)
                batch_end   = batches(ibatch,2)
                batchsz     = batch_end - batch_start + 1
                call read_imgbatch(batchsz, pinds(batch_start:batch_end), [1,batchsz] )
                cls2batch = (/(i,i=batch_start,batch_end)/)
                ! flags bad ice
                !$omp parallel do private(j,i,iptcl,ithr,ctfparms,tfun,ice_score,sdev_noise) default(shared)&
                !$omp proc_bind(close)
                do j = batch_start,batch_end
                    ithr  = omp_get_thread_num()+1
                    i     = j - batch_start +  1
                    iptcl = pinds(j)
                    call tmp_imgs(ithr)%copy_fast(build%imgbatch(i))
                    call tmp_imgs(ithr)%norm_noise(build%lmsk, sdev_noise)
                    call tmp_imgs(ithr)%mask(real(params%box)/2.-COSMSKHALFWIDTH-1., 'soft', backgr=0.)
                    call tmp_imgs(ithr)%fft
                    ctfparms = build%spproj%get_ctfparams(params%oritype, iptcl)
                    tfun     = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                    call tfun%calc_ice_frac(tmp_imgs(ithr), ctfparms, ice_score)
                    call build%spproj_field%set(iptcl,'ice',ice_score)
                    if( ice_score > 4.0 )then
                        pinds(i)      = 0
                        cls2batch(i)  = 0
                        states(iptcl) = 0
                    endif
                    call tmp_imgs(ithr)%copy_fast(build%imgbatch(i))
                    call tmp_imgs(ithr)%fft
                    call tmp_imgs(ithr)%shift2Dserial(-build%spproj_field%get_2Dshift(iptcl))
                    call tfun%apply_serial(tmp_imgs(ithr), 'flip', ctfparms)
                    call tmp_imgs(ithr)%bp(0.,4.)
                    call tmp_imgs(ithr)%ifft
                    call tmp_imgs(ithr)%norm_noise(build%lmsk, sdev_noise)
                    call tmp_imgs(ithr)%moments_within(build%lmsk, means(i), sdevs(i))
                    call build%spproj_field%set(iptcl,'mean',means(i))
                    call build%spproj_field%set(iptcl,'sdev',sdevs(i))
                enddo
                !$omp end parallel do
            enddo
            deallocate(means,sdevs)
            print *,icls,' Ice rejection: ', count(pinds==0)
            if( l_groundtruth )then
                purity(0) = real(sum(labels(pinds(:)),mask=(pinds>0)))
                purity(0) = purity(0) * 100. / real(count(pinds>0))
                print *, icls,count(pinds>0),purity(0)
            endif
            ! pftcc init
            pop   = count(pinds > 0)
            nbins = ceiling(real(pop)/real(NPTCLS_PER_BIN))
            if( nbins < 3 ) cycle
            call pftcc%new(NITERS, [1,pop], params%kfromto)
            call pftcc%reallocate_ptcls(pop, pack(pinds,mask=(pinds>0)))
            do i = 1,size(pinds)
                if( pinds(i)>0 )then
                    iptcl = pinds(i)
                    exit
                endif
            enddo
            l_ctf = build%spproj%get_ctfflag(params%oritype,iptcl=iptcl).ne.'no'
            if( l_ctf ) call pftcc%create_polar_absctfmats(build%spproj, params%oritype)
            ! polar representation
            !$omp parallel do private(j,i,iptcl,ithr) default(shared) proc_bind(close)
            do i = 1,size(pinds)
                ithr  = omp_get_thread_num()+1
                iptcl = pinds(i)
                if(iptcl == 0) cycle
                call tmp_imgs(ithr)%copy_fast(build%imgbatch(i))
                call prepimg4align(iptcl, build%imgbatch(i), tmp_imgs(ithr))
                call build%img_crop_polarizer%polarize(pftcc, tmp_imgs(ithr), iptcl, .true., .true.)
                call pftcc%set_eo(iptcl, (build%spproj_field%get_eo(iptcl)==0))
            enddo
            !$omp end parallel do
            ! indexing update
            cls2batch = pack(cls2batch,mask=(pinds>0))
            pinds     = pack(pinds,mask=(pinds>0))
            if( l_write )then
                j = 0
                if( cls2batch(1) > 1 )then
                    do i = 1,cls2batch(1)-1
                        j = j+1
                        call build%imgbatch(i)%write('ice_'//int2str_pad(icls,3)//'.mrc',j)
                    enddo
                endif
                do i = 1,pop-1
                    if( cls2batch(i+1) == cls2batch(i)+1 ) cycle
                    do k = cls2batch(i)+1,cls2batch(i+1)-1,1
                        j = j+1
                        call build%imgbatch(k)%write('ice_'//int2str_pad(icls,3)//'.mrc',j)
                    enddo
                enddo
                if( cls2batch(pop) < ini_pop )then
                    do i = cls2batch(pop)+1,ini_pop
                        j = j+1
                        call build%imgbatch(i)%write('ice_'//int2str_pad(icls,3)//'.mrc',j)
                    enddo
                endif
            endif
            ! class allocations
            if( allocated(corrs) ) deallocate(corrs,order,weights,R2s,RmI2s,rots,shifts,dfs,&
                &bin_inds,bins,selected,binccs,bindiff)
            allocate(corrs(pop),order(pop),weights(pop),R2s(nbins),RmI2s(nbins),rots(pop),shifts(2,pop),&
                &dfs(pop),bin_inds(nbins,NPTCLS_PER_BIN),bins(pop),selected(pop),binccs(nbins),bindiff(nbins))
            ! alignement inititialization
            !$omp parallel do private(j,i,iptcl,ithr) default(shared) proc_bind(close)
            do i = 1,pop
                iptcl       = pinds(i)
                corrs(i)    = build%spproj_field%get(iptcl,'corr')
                dfs(i)      = (build%spproj_field%get_dfx(iptcl)+build%spproj_field%get_dfy(iptcl))/2.
                rots(i)     = pftcc%get_roind(360.-build%spproj_field%e3get(iptcl))
                shifts(:,i) = 0.
                weights(i)  = 1.d0
                selected(i) = .true.
            enddo
            !$omp end parallel do
            ! defocus vs. score correlation
            mdf = sum(dfs) / real(pop)
            sdf = sqrt(sum((dfs-mdf)**2)/real(pop))
            dfs = (dfs-mdf)/sdf
            mcorr = sum(corrs) / real(pop)
            scorr = sqrt(sum((corrs-mcorr)**2)/real(pop))
            tmp   = (corrs-mcorr)/scorr
            cc_df_corrs = sum(dfs*tmp) / real(pop)
            if( trim(params%reject_cls).eq.'yes' )then
                ! rejecting entire classes with zero and less correlation
                if( cc_df_corrs < 0.0 )then
                    states(pinds(:)) = 0
                    cycle
                endif
            endif
            ! References & noise power in pftcc
            call restore_cavgs(pop, weights, optfilter=(params%cc_objfun==OBJFUN_EUCLID))
            if( l_weighted_init )then
                r2  = sum(csq_fast(cls_avg))
                !$omp parallel do private(i,irot,ptcl_rot,ctf_rot,diff) default(shared) proc_bind(close)
                do i = 1,pop
                    irot = pftcc%nrots+2-rots(i)
                    if(irot > pftcc%nrots ) irot = irot - pftcc%nrots
                    call pftcc%rotate_ptcl(pftcc%pfts_ptcls(:,:,i), irot, ptcl_rot)
                    call pftcc%rotate_ctf(pinds(i), irot, ctf_rot)
                    diff = ctf_rot * cls_avg - ptcl_rot
                    selected(i) = (sum(csq_fast(diff)) / r2) < 0.9
                enddo
                !$omp end parallel do
                pop_sel = count(selected)
                print *,'Deselected    : ',icls,pop-pop_sel
                if( pop_sel < pop )then
                    nbins = ceiling(real(pop_sel)/real(NPTCLS_PER_BIN))
                    deallocate(R2s,RmI2s,bin_inds)
                    allocate(R2s(nbins),RmI2s(nbins),bin_inds(nbins,NPTCLS_PER_BIN))
                    j = 0
                    do i = 1,pop
                        if(.not.selected(i))then
                            states(pinds(i)) = 0
                            if( l_write )then
                                j = j+1
                                call build%imgbatch(i)%write('bad_'//int2str_pad(icls,3)//'.mrc',j)
                            endif
                        endif
                    enddo
                endif
            else
                pop_sel = pop
            endif
            if( nbins < 3 ) cycle
            where(.not.selected) weights = 0.d0
            call restore_cavgs(pop, weights, optfilter=(params%cc_objfun==OBJFUN_EUCLID))
            call write_cls(cls_avg, 'cls_'//int2str_pad(icls,3)//'_iter.mrc', 1)
            pftcc%pfts_refs_even(:,:,1) = cmplx(cls_avg_even,kind=sp)
            pftcc%pfts_refs_odd(:,:,1)  = cmplx(cls_avg_odd, kind=sp)
            call pftcc%memoize_refs
            call pftcc%memoize_ptcls
            call update_sigmas(1)
            ! first scores
            !$omp parallel do private(i) default(shared) proc_bind(close)
            do i = 1,pop
                if( selected(i) )then
                    corrs(i) = real(pftcc%gencorr_for_rot_8(1, pinds(i), [0.d0,0.d0], rots(i)))
                else
                    corrs(i) = -2.
                endif
            enddo
            !$omp end parallel do
            ! negative correlations
            if( count(selected .and.(corrs<1.e-6)) > 0 )then
                if( l_write )then
                    j = 0
                    do i = 1,pop
                        if( corrs(i) < 0. .and. selected(i))then
                            j = j+1
                            print *,icls,i,cls2batch(i), corrs(i), labels(pinds(i))
                            call build%imgbatch(cls2batch(i))%write('negcorr_'//int2str_pad(icls,3)//'.mrc',j)
                        endif
                    enddo
                endif
                where( selected .and.(corrs<1.e-6) ) selected = .false.
                pop_sel = count(selected)
                print *,'Deselected neg: ',icls,pop-pop_sel
                nbins = ceiling(real(pop_sel)/real(NPTCLS_PER_BIN))
                if( nbins < 3 ) cycle
                deallocate(R2s,RmI2s,bin_inds)
                allocate(R2s(nbins),RmI2s(nbins),bin_inds(nbins,NPTCLS_PER_BIN))
                where(.not. selected)
                    states(pinds(:)) = 0
                    weights(:) = 0.d0
                end where
                call restore_cavgs(pop, weights, optfilter=(params%cc_objfun==OBJFUN_EUCLID))
                call write_cls(cls_avg, 'cls_'//int2str_pad(icls,3)//'_iter.mrc', 1)
                pftcc%pfts_refs_even(:,:,1) = cmplx(cls_avg_even,kind=sp)
                pftcc%pfts_refs_odd(:,:,1)  = cmplx(cls_avg_odd, kind=sp)
                call pftcc%memoize_refs
                call update_sigmas(1)
            endif
            ! Iteration loop
            cls_avg_bak    = cls_avg
            prev_threshold = huge(threshold)
            prev_nsel      = 0
            do iter = 1,NITERS
                ! generate bins
                call partition_cls
                ! bin-based thesholding
                do ibin = 1,nbins
                    binpop = count(bins==ibin)
                    weights = 0.d0
                    where( bins == ibin ) weights = 1.d0
                    call restore_cavgs(pop, weights)
                    call write_cls(cls_avg, 'cls_'//int2str_pad(icls,3)//'.mrc', ibin)
                    R2s(ibin) = sum(csq_fast(cls_avg))
                    rmi2 = 0.d0
                    !$omp parallel do private(i,j,ithr,irot,diff,ptcl,ptcl_rot,ctf_rot)&
                    !$omp reduction(+:rmi2) default(shared) proc_bind(close)
                    do i = 1,binpop
                        j = bin_inds(ibin,i)
                        if( j == 0 ) cycle
                        ithr = omp_get_thread_num()+1
                        call pftcc%gen_shmat(ithr, -real(shifts(:,j)), pftcc%heap_vars(ithr)%shmat)
                        ptcl = pftcc%pfts_ptcls(:,:,j) * pftcc%heap_vars(ithr)%shmat
                        irot = pftcc%nrots+2-rots(j)
                        if(irot > pftcc%nrots ) irot = irot - pftcc%nrots
                        call pftcc%rotate_ptcl(ptcl,    irot, ptcl_rot)
                        call pftcc%rotate_ctf(pinds(j), irot, ctf_rot)
                        diff = ctf_rot * cls_avg - ptcl_rot
                        RmI2 = RmI2 + sum(csq_fast(diff))
                    enddo
                    !$omp end parallel do
                    RmI2s(ibin) = rmi2 / (binpop-1)
                    binccs(ibin) = real(sum(cls_avg_bak*conjg(cls_avg))) / sqrt(sum(csq_fast(cls_avg_bak)) * sum(csq_fast(cls_avg)))
                    bindiff(ibin) = sum(csq_fast(cls_avg_bak-cls_avg))
                enddo
                R2s = R2s / RmI2s
                where( binccs < 0. ) R2s = 0.
                mean = sum(bindiff,mask=(binccs>0.)) / real(count(binccs>0.))
                sdev = sqrt(sum((bindiff-mean)**2,mask=(binccs>0.))/real(count(binccs>0.)))
                where( bindiff > mean+2.*sdev ) R2s = 0.
                ! Threshold
                tmp = real(R2s)
                call hpsort(tmp)
                ! threshold = median(tmp(nbins-4:nbins)) / 2.
                ! threshold = median(tmp(nbins-4:nbins)) / 3.
                ! threshold = median(tmp(nbins-4:nbins)) / 4.
                ! threshold = sum(tmp(nbins-4:nbins)) / 5. / 3.
                threshold = median(tmp(floor(real(nbins)*0.9):nbins)) / 4.
                if( (iter==1) .and. (count(R2s<threshold)==0) )then
                    ! making sure the weaker bin is deactivated on first iteration
                    R2s(minloc(R2s,dim=1)) = threshold - 1.
                endif
                ! Bins rejection
                k = 0
                do ibin = 1,nbins
                    binpop = count(bins==ibin)
                    pu = 0.
                    ! df = 0.0
                    ! cc = 0.0
                    do i = 1,binpop
                        j = bin_inds(ibin,i)
                        if(j == 0) cycle
                        iptcl = pinds(j)
                        if( R2s(ibin) < threshold )then
                            weights(j) = 0.d0
                        else
                            weights(j) = 1.d0
                        endif
                        ! df = df + build%spproj_field%get(iptcl,'dfx')
                        if( l_groundtruth ) pu = pu + real(labels(iptcl))
                        ! k = k + 1
                        ! cc = cc + corrs(k)
                    enddo
                    print *,icls,iter,ibin,100.*pu/real(binpop),R2s(ibin), RmI2s(ibin), binccs(ibin),  bindiff(ibin)
                enddo
                nsel = count(weights > 0.5d0)        
                if( l_groundtruth )then
                    purity(iter) = sum(labels(pinds),mask=(weights>0.5d0))
                    purity(iter) = 100. * purity(iter) / real(nsel)
                endif
                ! Convergence
                if( iter >= 2 )then
                    ! early exit
                    if( (nsel > prev_nsel) .and. (threshold < prev_threshold) ) exit
                endif
                where(weights > 0.5d0)
                    states(pinds(:)) = 1
                elsewhere
                    states(pinds(:)) = 0
                end where
                if( l_groundtruth )then
                    print *,icls,iter,threshold,purity(iter),nsel,pop_sel,cc_df_corrs
                else
                    print *,icls,iter,threshold,nsel,pop_sel,cc_df_corrs
                endif
                if ( iter == NITERS ) exit
                prev_threshold = threshold
                prev_nsel      = nsel
                ! re-scoring
                call restore_cavgs(pop, weights, optfilter=(params%cc_objfun==OBJFUN_EUCLID))
                cls_avg_bak = cls_avg
                call write_cls(cls_avg, 'cls_'//int2str_pad(icls,3)//'_iter.mrc', iter+1)
                pftcc%pfts_refs_even(:,:,iter) = cmplx(cls_avg_even)
                pftcc%pfts_refs_odd(:,:,iter)  = cmplx(cls_avg_odd)
                call pftcc%memoize_refs
                !$omp parallel do private(i,ithr,irot,iptcl,inpl_corrs,cxy) default(shared) proc_bind(close)
                do i = 1,pop
                    if( selected(i) )then
                        ithr     = omp_get_thread_num()+1
                        iptcl    = pinds(i)
                        ! rotation + shift
                        call pftcc%gencorrs(iter, iptcl, inpl_corrs)
                        irot = maxloc(inpl_corrs, dim=1)
                        call grad_shsrch_objs(ithr)%set_indices(iter, iptcl)
                        cxy = grad_shsrch_objs(ithr)%minimize(irot=irot)
                        if( irot > 0 )then
                            corrs(i)    = cxy(1)
                            shifts(:,i) = cxy(2:3)
                        else
                            irot        = maxloc(inpl_corrs, dim=1)
                            corrs(i)    = inpl_corrs(irot)
                            shifts(:,i) = 0.
                        endif
                        rots(i)  = irot
                        ! rotation
                        ! call pftcc%gencorrs(iter, iptcl, inpl_corrs)
                        ! irot     = maxloc(inpl_corrs,dim=1)
                        ! corrs(i) = inpl_corrs(irot)
                        ! rots(i)  = irot
                        ! nothing
                        ! irot = pftcc%get_roind(360.-build%spproj_field%e3get(iptcl))
                        ! corrs(i) = real(pftcc%gencorr_for_rot_8(iter,iptcl, [0.d0,0.d0], irot))
                    else
                        corrs(i) = -2.
                    endif
                enddo
                !$omp end parallel do
                call update_sigmas(iter)
            enddo
            k = 0
            do ibin = 1,nbins
                binpop = count(bins==ibin)
                do i = 1,binpop
                    j = bin_inds(ibin,i)
                    if( j==0 ) cycle
                    k = k+1
                    call build%imgbatch(cls2batch(j))%write('ptcls_'//int2str_pad(icls,3)//'.mrc',k)
                enddo
            enddo
        enddo
        call build%spproj_field%set_all('state', real(states))
        call build%spproj%write_segment_inside(params%oritype, params%projfile)
        call build%spproj_field%write('ptcl2d.txt')
        print *,'NREJECTED     : ', count(states==0), count(states==1)
        if( l_groundtruth )then
            print *,'TRUE REJECTED : ', count(states==0 .and. labels==0)
            print *,'FALSE REJECTED: ', count(states==0 .and. labels==1)
            print *,'TRUE KEPT     : ', count(states==1 .and. labels==1)
            print *,'FALSE KEPT    : ', count(states==1 .and. labels==0)
            print *,'PURITY        : ', 100. * real(count(states==1 .and. labels==1)) / real(count(states==1))
        endif
        where( states == 0 ) states=2
        where( states == 1 ) states=0
        where( states == 2 ) states=1
        call build%spproj_field%set_all('state', real(states))
        call build%spproj%write('inverted.simple')
        call simple_end('**** SIMPLE_PRUNE_CAVGS NORMAL STOP ****')
        contains

            subroutine partition_cls()
                integer, allocatable :: corrbins(:,:)
                integer :: ibin, i, j, k, offset
                bin_inds = 0
                bins     = 0
                order    = (/(i,i=1,pop)/)
                offset = pop - count(selected)
                call hpsort(corrs,order)
                corrbins = split_nobjs_even(pop-offset, nbins)
                do ibin = 1,nbins
                    j = 0
                    do i = corrbins(ibin,1),corrbins(ibin,2)
                        j = j + 1
                        k = order(i+offset)
                        bin_inds(ibin,j) = k
                        bins(k) = ibin
                    enddo
                enddo
            end subroutine partition_cls

            subroutine update_sigmas( iref )
                integer, intent(in) :: iref
                real    :: sig2_contrib(params_glob%kfromto(1):params_glob%kfromto(2))
                integer :: i, iptcl, istk, neven(nstks), nodd(nstks)
                if( params%cc_objfun /= OBJFUN_EUCLID ) return
                call pftcc%assign_sigma2_noise(sig2)
                sig2_even = 0.d0
                sig2_odd  = 0.d0
                neven = 0
                nodd  = 0
                do i =1,pop
                    if(.not.selected(i))cycle
                    iptcl = pinds(i)
                    call pftcc%gencorr_sigma_contrib(iref, iptcl, shifts(:,i), rots(i), sig2_contrib)
                    istk = nint(build%spproj_field%get(iptcl,'stkind'))
                    if( pftcc%iseven(i) )then
                        neven(istk) = neven(istk) + 1
                        sig2_even(:,istk) = sig2_even(:,istk) + sig2_contrib
                    else
                        nodd(istk) = nodd(istk) + 1
                        sig2_odd(:,istk)  = sig2_odd(:,istk)  + sig2_contrib
                    endif
                enddo
                ! do istk = 1,nstks
                !     sig2_even(:,istk) = sig2_even(:,istk) / real(neven(istk))
                !     sig2_odd(:,istk)  = sig2_odd(:,istk)  / real(nodd(istk))
                ! enddo
                sig2_even(:,1) = sum(sig2_even,dim=2) / real(sum(neven))
                sig2_odd(:,1) = sum(sig2_odd,dim=2) / real(sum(nodd))
                do i = 1,pop
                    iptcl = pinds(i)
                    ! istk  = nint(build%spproj_field%get(iptcl,'stkind'))
                    if( pftcc%iseven(i) )then
                        pftcc%sigma2_noise(:,iptcl) = sig2_even(:,1)
                    else
                        pftcc%sigma2_noise(:,iptcl) = sig2_odd(:,1)
                    endif
                    call pftcc%memoize_sqsum_ptcl(iptcl)
                enddo
            end subroutine update_sigmas

            subroutine restore_cavgs( n, weights, optfilter )
                integer,           intent(in) :: n
                real(dp),          intent(in) :: weights(n)
                logical, optional, intent(in) :: optfilter
                real(dp) :: w,cc,g
                integer  :: i,k,ithr
                num_even   = 0.d0
                denom_even = 0.d0
                num_odd    = 0.d0
                denom_odd  = 0.d0
                !$omp parallel do private(i,w,ithr,ptcl,ptcl_rot,ctf_rot,irot) default(shared) proc_bind(close)&
                !$omp reduction(+:num_even,num_odd,denom_even,denom_odd)
                do i = 1,n
                    ithr = omp_get_thread_num()+1
                    w    = weights(i)
                    if( w < 1.d-6 ) cycle
                    call pftcc%gen_shmat(ithr, -real(shifts(:,i)), pftcc%heap_vars(ithr)%shmat)
                    ptcl = pftcc%pfts_ptcls(:,:,i) * pftcc%heap_vars(ithr)%shmat
                    irot = pftcc%nrots+2-rots(i)
                    if(irot > pftcc%nrots ) irot = irot - pftcc%nrots
                    call pftcc%rotate_ptcl(ptcl, irot, ptcl_rot)
                    call pftcc%rotate_ctf(pinds(i), irot, ctf_rot)
                    if(pftcc%iseven(i))then
                        num_even   = num_even   + w * ptcl_rot * ctf_rot
                        denom_even = denom_even + w * ctf_rot**2
                    else
                        num_odd   = num_odd   + w * ptcl_rot * ctf_rot
                        denom_odd = denom_odd + w * ctf_rot**2
                    endif
                enddo
                !$omp end parallel do
                ! e/o drift
                k = min(params%kfromto(1)+2,params%kfromto(2)-1)
                num_even(:,params%kfromto(1):k)   = (num_even(:,params%kfromto(1):k)  +num_odd(:,params%kfromto(1):k))   / 2.d0
                denom_even(:,params%kfromto(1):k) = (denom_even(:,params%kfromto(1):k)+denom_odd(:,params%kfromto(1):k)) / 2.d0
                num_odd(:,params%kfromto(1):k)    = num_even(:,params%kfromto(1):k)
                denom_odd(:,params%kfromto(1):k)  = denom_even(:,params%kfromto(1):k)
                ! restoration & frc
                cls_avg_even = num_even / (denom_even)
                cls_avg_odd  = num_odd / (denom_odd)
                cls_avg      = (num_even+num_odd) / (denom_even+denom_odd)
                do k = params%kfromto(1),params%kfromto(2)
                    frc(k) = sum(real(cls_avg_even(:,k)*conjg(cls_avg_odd(:,k))))
                    frc(k) = frc(k) / sqrt(sum(csq_fast(cls_avg_even(:,k))) * sum(csq_fast(cls_avg_odd(:,k))))
                enddo
                if( present(optfilter) )then
                    ! optimal filtering
                    if( optfilter )then
                        do k = params%kfromto(1),params%kfromto(2)
                            cc = max(frc(k),0.)
                            g  = max(0.d0, min(2.*cc/(1.d0+cc), 1.d0))
                            if( cc < 0.143 ) g=0.d0
                            cls_avg(:,k)      = cls_avg(:,k)      * g
                            cls_avg_even(:,k) = cls_avg_even(:,k) * g
                            cls_avg_odd(:,k)  = cls_avg_odd(:,k)  * g
                        enddo
                    endif
                endif
            end subroutine restore_cavgs

            subroutine write_cls( I, fname, idx )
                complex(dp),      intent(in) :: I(pftcc%pftsz,params%kfromto(1):params%kfromto(2))
                character(len=*), intent(in) :: fname
                integer,          intent(in) :: idx
                complex, allocatable :: cmat(:,:)
                integer :: box
                call pftcc%polar2cartesian(cmplx(I, kind=sp), cmat, box)
                call build%img%new([box,box,1], params%smpd*real(params%box)/real(box))
                call build%img%zero_and_flag_ft
                call build%img%set_cmat(cmat)
                call build%img%shift_phorig()
                call build%img%ifft
                call build%img%write(fname,idx)
            end subroutine write_cls

            subroutine detect_class_outliers(stkname)
                character(len=*), intent(in) :: stkname
                integer, parameter :: NHISTBINS = 256
                real,    parameter :: NSIG = 4.
                type(image)          :: tmpimg
                real,        pointer :: prmat(:,:,:)
                logical, allocatable :: lmsk(:,:,:), cls_mask(:)
                real    :: p(NHISTBINS), ps(NHISTBINS,ncls), vals(ncls), minmax(2)
                real    :: overall_min, overall_max, dr, hrange, mean, std, threshold
                integer :: icls,i,j,b,zerobin,n
                cls_mask = cls_pops > 0
                n = count(cls_mask)
                ! build histograms
                call tmpimg%disc([params%box,params%box,1],params%smpd,real(params%box)/2.-COSMSKHALFWIDTH,lmsk)
                call tmpimg%kill
                overall_min =  999.
                overall_max = -999.
                do icls = 1,ncls
                    if( .not.cls_mask(icls) ) cycle
                    call build%img%read(stkname,icls)
                    call build%img%mask(real(params%box)/2.-COSMSKHALFWIDTH, 'hard')
                    minmax = build%img%minmax()
                    overall_min = min(minmax(1),overall_min)
                    overall_max = max(minmax(2),overall_max)
                enddo
                dr = (overall_max-overall_min) / real(NHISTBINS)
                call build%img%get_rmat_ptr(prmat)
                ps = 0.
                do icls = 1,cavgs_ncls
                    if( .not.cls_mask(icls) ) cycle
                    call build%img%read(stkname,icls)
                    call build%img%mask(real(params%box)/2.-COSMSKHALFWIDTH, 'hard')
                    do i = 1,NHISTBINS
                        minmax(1) = overall_min + real(i-1)*dr
                        minmax(2) = minmax(1) + dr
                        p(i) = real(count(lmsk(:,:,1) .and. (prmat(1:params%box,1:params%box,1) > minmax(1))&
                            &.and. (prmat(1:params%box,1:params%box,1) < minmax(2)) ))
                    enddo
                    p = p / sum(p)
                    ps(:,icls) = p
                enddo
                ! overall histogram
                p = 0.
                do icls = 1,cavgs_ncls
                    if( .not.cls_mask(icls) ) cycle
                    p = p + ps(:,icls)
                enddo
                p = p / real(n)
                zerobin = ceiling(-overall_min/dr)
                ! outliers
                vals = 0.
                do icls = 1,cavgs_ncls
                    if( .not.cls_mask(icls) ) cycle
                    vals(icls) = sum(abs(ps(:,icls)-p(:)))
                enddo
                mean = sum(vals,mask=cls_mask) / real(n)
                std  = sqrt(sum((vals-mean)**2,mask=cls_mask)/real(n))
                print *, 'Brightness mean/std/threshold: ', mean, std, threshold
                threshold = mean + NSIG * std
                do icls = 1,cavgs_ncls
                    if( .not.cls_mask(icls) ) cycle
                    if( abs(vals(icls)) > threshold ) print *,'Bright class to reject: ',icls, vals(icls)
                enddo
            end subroutine detect_class_outliers

    end subroutine exec_prune_cavgs

end module simple_commander_resolest
