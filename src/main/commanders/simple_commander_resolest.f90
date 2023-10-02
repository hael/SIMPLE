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
        type(image)                   :: img_cavg, calc_cavg
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
        call prepimgbatch(nptcls)
        call read_imgbatch([1, nptcls])
        call build%img_match%init_polarizer(pftcc, params%alpha)
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
            call prepimg4align(iptcl, build%imgbatch(iptcl))
            ! transfer to polar coordinates
            call build%img_match%polarize(pftcc, build%imgbatch(iptcl), iptcl, .true., .true., mask=build%l_resmsk)
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
        integer, parameter :: NDFPARTS       = 3
        integer, parameter :: NPTCLS_PER_BIN = 50
        type(pftcc_shsrch_grad), allocatable :: grad_shsrch_objs(:)
        type(polarft_corrcalc)        :: pftcc
        type(builder)                 :: build
        type(parameters)              :: params
        type(image),      allocatable :: tmp_imgs(:)
        complex(dp),      allocatable :: cls_avg(:,:), cls_avg_even(:,:), cls_avg_odd(:,:), num_even(:,:), num_odd(:,:), diff(:,:)
        complex(sp),      allocatable :: ptcl(:,:), ptcl_rot(:,:)
        real(dp),         allocatable :: denom_even(:,:), denom_odd(:,:), R2s(:), RmI2s(:), weights(:)
        real(sp), target, allocatable :: sig2(:,:)
        real(sp),         allocatable :: purity(:),inpl_corrs(:), corrs(:), ctf_rot(:,:), shifts(:,:), scores(:), dfs(:)
        real(sp),         allocatable :: kweights(:,:), frc(:), sig2_even(:,:), sig2_odd(:,:), tmp(:)
        integer,          allocatable :: rots(:),pinds(:), states(:), order(:), labels(:), batches(:,:), df_order(:)
        integer,          allocatable :: bin_inds(:,:), bins(:)
        character(len=:), allocatable :: cavgsstk
        type(ctf)       :: tfun
        type(ctfparams) :: ctfparms
        real(dp) :: rmi2
        real     :: cxy(3), lims(2,2), lims_init(2,2), threshold, cc, pu,var, ccmax, ice_score,df
        integer  :: nstks, nptcls, iptcl, iter, n_lines, icls, nbins, batch_start, batch_end, ibatch, batchsz, ibin
        integer  :: ncls, i, j, k,fnr, istart,iend, irot, nptcls_sel, pop, nbatches, batchsz_max, ithr, ngoodice, nbadice, n, nb, ng
        integer  :: ndfbins, binpop
        logical  :: l_ctf, l_groundtruth, l_corr_ranking
        call cline%set('oritype', 'ptcl2D')
        call cline%set('mkdir',   'yes')
        if( .not.cline%defined('objfun') ) call cline%set('objfun',  'cc')
        call build%init_params_and_build_general_tbox(cline, params)
        ncls       = build%spproj_field%get_n('class')
        states     = nint( build%spproj_field%get_all('state'))
        nptcls     = size(states)
        nptcls_sel = count(states==1)
        l_groundtruth  = cline%defined('infile')
        l_corr_ranking = .true.
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
            call tmp_imgs(ithr)%new([params%box,params%box,1],params%smpd,wthreads=.false.)
        enddo
        n = 0
        nb = 0
        ng = 0
        ! Class loop
        do icls = 1,ncls
            call build%spproj_field%get_pinds(icls, 'class', pinds)
            pop   = size(pinds)
            nbins = ceiling(real(pop)/real(NPTCLS_PER_BIN))
            if( nbins < 5 ) cycle
            ndfbins  = ceiling(real(pop)/real(NDFPARTS))
            ! pftcc init
            call pftcc%new(NITERS, [1,pop], params%kfromto)
            call pftcc%reallocate_ptcls(pop, pinds)
            l_ctf = build%spproj%get_ctfflag(params%oritype,iptcl=pinds(1)).ne.'no'
            if( l_ctf ) call pftcc%create_polar_absctfmats(build%spproj, params%oritype)
            do ithr = 1, params%nthr
                call grad_shsrch_objs(ithr)%new(lims, lims_init=lims_init,&
                &shbarrier=params%shbarrier, maxits=60, opt_angle=.true.)
            enddo
            if( .not.allocated(ptcl) )then
                ! One-time allocations
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
                call build%img_match%init_polarizer(pftcc, params_glob%alpha)
                if( (params%cc_objfun==OBJFUN_EUCLID) )then
                    allocate(sig2(params%kfromto(1):params%kfromto(2),nptcls),&
                    &sig2_even(params%kfromto(1):params%kfromto(2),nstks),&
                    &sig2_odd(params%kfromto(1):params%kfromto(2),nstks))
                endif
            endif
            if( allocated(corrs) ) deallocate(corrs,order,weights,R2s,RmI2s,rots,shifts,df_order,dfs,&
            &bin_inds,bins)
            allocate(corrs(pop),order(pop),weights(pop),R2s(nbins),RmI2s(nbins),rots(pop),shifts(2,pop),df_order(pop),&
                &dfs(pop),bin_inds(nbins,ndfbins),bins(pop))
            ! batchsz_max = min(pop,params_glob%nthr*BATCHTHRSZ)
            batchsz_max = pop
            nbatches    = ceiling(real(pop)/real(batchsz_max))
            batches     = split_nobjs_even(pop, nbatches)
            batchsz_max = maxval(batches(:,2)-batches(:,1)+1)
            if( allocated(build%imgbatch) )then
                if( batchsz_max > size(build%imgbatch) ) call prepimgbatch(batchsz_max)
            else
                call prepimgbatch(batchsz_max)
            endif
            ! images prep
            do ibatch=1,nbatches
                batch_start = batches(ibatch,1)
                batch_end   = batches(ibatch,2)
                batchsz     = batch_end - batch_start + 1
                !print *,icls,batchsz_max,batch_start,batch_end,batchsz,nbatches,pop
                call read_imgbatch(batchsz, pinds(batch_start:batch_end), [1,batchsz] )
                !$omp parallel do private(j,i,iptcl,ithr) default(shared) proc_bind(close)
                do j = batch_start,batch_end
                    ithr  = omp_get_thread_num()+1
                    i     = j - batch_start +  1
                    iptcl = pinds(j)
                    call tmp_imgs(ithr)%copy_fast(build%imgbatch(i))
                    call prepimg4align(iptcl, tmp_imgs(ithr))
                    call build%img_match%polarize(pftcc, tmp_imgs(ithr), iptcl, .true., .true.)
                    call pftcc%set_eo(iptcl, (build%spproj_field%get_eo(iptcl)==0))
                enddo
                !$omp end parallel do
            enddo
            ! !$omp parallel do private(i,iptcl,tfun,ctfparms,ice_score) default(shared) proc_bind(close)
            ! do i = 1,pop
            !     iptcl = pinds(i)
            !     ctfparms = build%spproj%get_ctfparams(params%oritype, iptcl)
            !     tfun     = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
            !     call build%imgbatch(i)%fft
            !     call tfun%calc_line_frac(build%imgbatch(i), ctfparms, ice_score)
            !     call build%spproj_field%set(iptcl,'specscore',ice_score)
            ! enddo
            ! !$omp end parallel do
            ! ngoodice = 0 
            ! nbadice  = 0
            ! pu = 0.
            ! j = 0
            ! do i =1,pop
            !     iptcl = pinds(i)
            !     call build%imgbatch(i)%ifft
            !     if( build%spproj_field%get(iptcl,'specscore') > 0.85 )then
            !         if( labels(iptcl)==1 )then
            !             ng = ng + 1
            !             ngoodice = ngoodice+1
            !             call build%imgbatch(i)%write('gice.mrc',ng)
            !         endif
            !         if( labels(iptcl)==0 )then
            !             nb = nb + 1
            !             nbadice = nbadice+1
            !             call build%imgbatch(i)%write('bice.mrc',nb)
            !             print *,nb,i,build%spproj_field%get(iptcl,'specscore')
            !         endif
            !         n = n + 1
            !     else
            !         j = j+1
            !         call build%imgbatch(i)%write('no_ice_'//int2str_pad(icls,3)//'.mrc',j)
            !     endif
            !     pu = pu + real(labels(iptcl))
            ! enddo
            ! print *,icls, pop, 100.*pu/pop, ngoodice, nbadice, 100.*ngoodice/pop, 100.*nbadice/pop, 100.*(ngoodice+nbadice)/pop
            ! cycle
            ! more init
            !$omp parallel do private(i,iptcl) default(shared) proc_bind(close)
            do i = 1,pop
                iptcl       = pinds(i)
                rots(i)     = pftcc%get_roind(360.-build%spproj_field%e3get(iptcl))
                shifts(:,i) = 0.
                corrs(i)    = build%spproj_field%get(iptcl,'corr')
                dfs(i)      = (build%spproj_field%get_dfx(iptcl)+build%spproj_field%get_dfx(iptcl))/2.
                weights(i)  = 1.d0
                order(i)    = i
                df_order(i) = i
            enddo
            !$omp end parallel do
            if( l_groundtruth )then
                purity(0) = 0.
                do i = 1,pop
                    if( weights(i) > 0.5d0 ) purity(0) = purity(0) + real(labels(pinds(i)))
                enddo
                purity(0) = purity(0) *100./real(pop)
                print *, icls,purity(0)
            endif
            if( .not.l_corr_ranking ) call hpsort(dfs,df_order)
            ! References & noise power in pftcc
            call restore_cavgs(pop, weights, optfilter=(params%cc_objfun==OBJFUN_EUCLID))
            call write_cls(cls_avg,      'cls_'//int2str_pad(icls,3)//'_iter.mrc', 1)
            call write_cls(cls_avg_even, 'cls_even_'//int2str_pad(icls,3)//'_iter.mrc', 1)
            call write_cls(cls_avg_odd,  'cls_odd_'//int2str_pad(icls,3)//'_iter.mrc', 1)
            pftcc%pfts_refs_even(:,:,1) = cmplx(cls_avg_even,kind=sp)
            pftcc%pfts_refs_odd(:,:,1)  = cmplx(cls_avg_odd, kind=sp)
            call pftcc%memoize_refs
            call pftcc%memoize_ptcls
            call update_sigmas(1)
            ! to calculate first scores
            !$omp parallel do private(i,iptcl,inpl_corrs) default(shared) proc_bind(close)
            do i = 1,pop
                iptcl = pinds(i)
                corrs(i) = real(pftcc%gencorr_for_rot_8(1, iptcl, [0.d0,0.d0], rots(i)))
            enddo
            !$omp end parallel do
            do iter = 1,NITERS
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
                    RmI2s(ibin) = rmi2
                    RmI2s(ibin) = RmI2s(ibin) / (pop-1)
                enddo
                R2s = R2s / RmI2s
                ! Threshold
                tmp = real(R2s)
                call hpsort(tmp)
                ! threshold = median(tmp(nbins-4:nbins)) / 2.
                ! threshold = median(tmp(nbins-4:nbins)) / 3.
                ! threshold = median(tmp(nbins-4:nbins)) / 4.
                ! threshold = sum(tmp(nbins-4:nbins)) / 5. / 3.
                threshold = median(tmp(floor(real(nbins)*0.9):nbins)) / 3.
                if( (iter==1) .and. (count(R2s<threshold)==0) )then
                    ! making sure the weaker bin is deactivated on first iteration
                    R2s(minloc(R2s,dim=1)) = threshold - 1.
                endif
                do ibin = 1,nbins
                    binpop = count(bins==ibin)
                    pu = 0.
                    do i = 1,binpop
                        j = bin_inds(ibin,i)
                        if(j == 0) cycle
                        iptcl = pinds(j)
                        if( R2s(ibin) < threshold )then
                            weights(j) = 0.d0
                            call build%spproj_field%set_state(iptcl,0)
                        else
                            weights(j) = 1.d0
                            call build%spproj_field%set_state(iptcl,1)
                        endif
                        pu = pu + real(labels(iptcl))
                    enddo
                    ! print *,icls,iter,ibin,100.*pu/real(binpop),R2s(ibin), threshold
                enddo
                if( l_groundtruth )then
                    purity(iter) = 0.
                    do i = 1,pop
                        if( weights(i) > 0.5d0 ) purity(iter) = purity(iter) + real(labels(pinds(i)))
                    enddo
                    purity(iter) = purity(iter) *100./real(count(weights>0.5d0))
                endif
                ! New class
                call restore_cavgs(pop, weights, optfilter=(params%cc_objfun==OBJFUN_EUCLID))
                call write_cls(cls_avg, 'cls_'//int2str_pad(icls,3)//'_iter.mrc', iter+1)
                print *,icls,iter,threshold,purity(iter),count(weights>0.5d0),pop
                if ( iter == NITERS ) exit
                pftcc%pfts_refs_even(:,:,iter) = cmplx(cls_avg_even)
                pftcc%pfts_refs_odd(:,:,iter)  = cmplx(cls_avg_odd)
                call pftcc%memoize_refs
                !$omp parallel do private(i,ithr,irot,iptcl,inpl_corrs,cxy) default(shared) proc_bind(close)
                do i = 1,pop
                    ithr     = omp_get_thread_num()+1
                    iptcl    = pinds(i)
                    order(i) = i
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

                enddo
                !$omp end parallel do
                call update_sigmas(iter)
            enddo
        enddo
        call build%spproj%write_segment_inside(params%oritype, params%projfile)
        states = nint(build%spproj_field%get_all('state'))
        ! scores = build%spproj_field%get_all('specscore')
        ! itmp = states
        ! where( (scores>1.) ) states = 0
        ! call build%spproj_field%set_all('state', real(states))
        ! call build%spproj%write('all_expunged.simple')
        ! states = itmp
        ! itmp = states
        ! where( (scores>1.) .or. (labels == 0) ) states = 0
        ! call build%spproj_field%set_all('state', real(states))
        ! call build%spproj%write('good_expunged.simple')
        ! states = itmp
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
                real,    allocatable :: tmp(:)
                integer, allocatable :: dfbins(:,:), dfcorrs_bins(:,:), corrbins(:,:)
                integer :: curr_ind(nbins), dfbin, ibin, ndf, i, j, k
                bin_inds = 0
                bins     = 0
                if( l_corr_ranking )then
                    ! is equivalent to NDFPARTS=1
                    call hpsort(corrs,order)
                    corrbins = split_nobjs_even(pop, nbins)
                    do ibin = 1,nbins
                        j = 0
                        do i = corrbins(ibin,1),corrbins(ibin,2)
                            j = j + 1
                            k = order(i)
                            bin_inds(ibin,j) = k
                            bins(k) = ibin
                        enddo
                    enddo
                else
                    curr_ind = 0
                    dfbins   = split_nobjs_even(pop, NDFPARTS)
                    do dfbin = 1,NDFPARTS
                        ndf = dfbins(dfbin,2)-dfbins(dfbin,1)+1
                        if(allocated(tmp)  )deallocate(tmp)
                        if(allocated(order))deallocate(order)
                        allocate(tmp(ndf),order(ndf))
                        k = 0
                        do i = dfbins(dfbin,1),dfbins(dfbin,2)
                            k        = k + 1
                            j        = df_order(i)
                            tmp(k)   = corrs(j)
                            order(k) = j
                        enddo
                        call hpsort(tmp,order)
                        dfcorrs_bins = split_nobjs_even(ndf, nbins)
                        do ibin = 1,nbins
                            do i = dfcorrs_bins(ibin,1),dfcorrs_bins(ibin,2)
                                j = curr_ind(ibin) + 1
                                k = order(i)
                                bin_inds(ibin,j) = k
                                curr_ind(ibin)   = j
                                bins(k) = ibin
                            enddo
                        enddo
                    enddo
                endif
            end subroutine partition_cls

            subroutine update_sigmas( iref )
                integer, intent(in) :: iref
                real    :: sig2_contrib(params_glob%kfromto(1):params_glob%kfromto(2))
                integer :: i, iptcl, k, istk, neven(nstks), nodd(nstks)
                if( params%cc_objfun /= OBJFUN_EUCLID ) return
                call pftcc%assign_sigma2_noise(sig2)
                sig2_even = 0.d0
                sig2_odd  = 0.d0
                neven = 0
                nodd  = 0
                do i =1,pop
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
                do istk = 1,nstks
                    sig2_even(:,istk) = sig2_even(:,istk) / real(neven(istk))
                    sig2_odd(:,istk)  = sig2_odd(:,istk)  / real(nodd(istk))
                enddo
                do i = 1,pop
                    iptcl = pinds(i)
                    istk  = nint(build%spproj_field%get(iptcl,'stkind'))
                    if( pftcc%iseven(i) )then
                        pftcc%sigma2_noise(:,iptcl) = sig2_even(:,istk)
                    else
                        pftcc%sigma2_noise(:,iptcl) = sig2_odd(:,istk)
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

    end subroutine exec_prune_cavgs

end module simple_commander_resolest
