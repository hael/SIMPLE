! concrete commander: resolest for resolution estimation
module simple_commander_resolest
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,     only: parameters
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_image,          only: image
use simple_masker,         only: masker
use simple_fsc
implicit none

public :: fsc_commander
public :: nununiform_filter3D_commander
public :: uniform_filter2D_commander
public :: uniform_filter3D_commander
public :: icm3D_commander
public :: icm2D_commander
public :: denoise_cavgs_commander
public :: cavg_filter2D_commander
public :: score_ptcls_commander
public :: estimate_lpstages_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: fsc_commander
  contains
    procedure :: execute      => exec_fsc
end type fsc_commander

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

type, extends(commander_base) :: icm3D_commander
  contains
    procedure :: execute      => exec_icm3D
end type icm3D_commander

type, extends(commander_base) :: icm2D_commander
  contains
    procedure :: execute      => exec_icm2D
end type icm2D_commander

type, extends(commander_base) :: denoise_cavgs_commander
  contains
    procedure :: execute      => exec_denoise_cavgs
end type denoise_cavgs_commander

type, extends(commander_base) :: cavg_filter2D_commander
  contains
    procedure :: execute      => exec_cavg_filter2D
end type cavg_filter2D_commander

type, extends(commander_base) :: score_ptcls_commander
  contains
    procedure :: execute      => exec_score_ptcls
end type score_ptcls_commander

type, extends(commander_base) :: estimate_lpstages_commander
  contains
    procedure :: execute      => exec_estimate_lpstages
end type estimate_lpstages_commander

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
        integer :: j, k_hp, k_lp, nyq
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
            call even%zero_env_background(mskvol)
            call odd%zero_env_background(mskvol)
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
        call simple_end('**** SIMPLE_NONUNIFORM_FILTER3D NORMAL STOP ****')
    end subroutine exec_nununiform_filter3D

    subroutine exec_uniform_filter3D(self, cline)
        use simple_opt_filter, only: estimate_lplim
        class(uniform_filter3D_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters)  :: params
        type(image)       :: even, odd, mskvol, odd_filt
        real              :: lpopt
        character(len=90) :: file_tag
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call odd%new([params%box,params%box,params%box], params%smpd)
        call odd_filt%new([params%box,params%box,params%box], params%smpd)
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd %read(params%vols(1))
        call even%read(params%vols(2))
        file_tag = 'uniform_3D_filter'
        if( params%l_filemsk )then
            call mskvol%new([params%box,params%box,params%box], params%smpd)
            call mskvol%read(params%mskfile)
            call mskvol%one_at_edge ! to expand before masking of reference internally
        else
            call mskvol%disc([params%box,params%box,params%box], params%smpd,&
                    &real(min(params%box/2, int(params%msk + COSMSKHALFWIDTH))))
        endif
        ! soft masking needed for FT
        call even%mask(params%msk, 'soft')
        call  odd%mask(params%msk, 'soft')
        call estimate_lplim(odd, even, mskvol, [params%lpstart,params%lpstop], lpopt, odd_filt)
        print *, 'found optimal low-pass limit: ', lpopt
        call odd_filt%write('odd_filt.mrc')
        ! destruct
        call odd%kill
        call odd_filt%kill
        call even%kill
        call mskvol%kill
        ! end gracefully
        call simple_end('**** SIMPLE_UNIFORM_FILTER3D NORMAL STOP ****')
    end subroutine exec_uniform_filter3D

    subroutine exec_icm3D( self, cline )
        class(icm3D_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(parameters)     :: params
        type(image)          :: even, odd, even_icm, odd_icm, avg, avg_icm
        type(masker)         :: envmsk
        logical, allocatable :: l_msk(:,:,:)
        real,    allocatable :: fsc(:), res(:), pspec(:), pspec_icm(:)
        character(len=90)    :: file_tag
        real                 :: mskrad
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call odd %new([params%box,params%box,params%box], params%smpd)
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd %read(params%vols(1))
        call even%read(params%vols(2))
        call avg%copy(even)
        call avg%add(odd)
        call avg%mul(0.5)
        call even_icm%copy(even)
        call odd_icm%copy(odd)
        if( params%automsk.ne.'no' )then
            call envmsk%automask3D(even, odd, l_tight=params%automsk.eq.'tight')
            ! apply mask to volumes
            call even_icm%zero_env_background(envmsk)
            call odd_icm%zero_env_background(envmsk)
            call even_icm%mul(envmsk)
            call odd_icm%mul(envmsk)
            call envmsk%one_at_edge ! to expand before masking of reference internally
            l_msk = envmsk%bin2logical()
            call even_icm%ICM3D_eo(odd_icm, params%lambda, l_msk)
            call envmsk%kill
        else
            call even_icm%ICM3D_eo(odd_icm, params%lambda)
        endif
        call avg_icm%copy(even_icm)
        call avg_icm%add(odd_icm)
        call avg_icm%mul(0.5)
        file_tag = 'icm_3D_filter'
        call even_icm%write('vol2_'//trim(file_tag)//'.mrc')
        call odd_icm%write('vol1_'//trim(file_tag)//'.mrc')
        call avg_icm%write('vol_avg_'//trim(file_tag)//'.mrc')
        call even%kill
        call odd%kill
        call even_icm%kill
        call odd_icm%kill
        call avg%kill
        call avg_icm%kill
        call simple_end('**** SIMPLE_ICM3D NORMAL STOP ****')
    end subroutine exec_icm3D

    subroutine exec_uniform_filter2D( self, cline )
        use simple_opt_filter, only: estimate_lplims2D
        class(uniform_filter2D_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(image),      allocatable :: odd(:), even(:), odd_filt(:)
        real,             allocatable :: lpsopt(:)
        type(parameters) :: params
        integer          :: iptcl
        ! init
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call find_ldim_nptcls(params%stk, params%ldim, params%nptcls)
        params%ldim(3) = 1 ! because we operate on stacks
        ! allocate
        allocate(odd(params%nptcls), odd_filt(params%nptcls), even(params%nptcls), lpsopt(params%nptcls))
        ! construct & read
        do iptcl = 1, params%nptcls
            call odd( iptcl)%new( params%ldim, params%smpd, .false.)
            call odd_filt( iptcl)%new( params%ldim, params%smpd, .false.)
            call odd( iptcl)%read(params%stk,  iptcl)
            call even(iptcl)%new( params%ldim, params%smpd, .false.)
            call even(iptcl)%read(params%stk2, iptcl)
        enddo
        ! filter
        call estimate_lplims2D( odd, even, params%msk, [params%lpstart,params%lpstop], lpsopt, odd_filt )
        ! write output and destruct
        do iptcl = 1, params%nptcls
            print *, 'found optimal low-pass limit: ', iptcl, lpsopt(iptcl)
            call odd_filt(iptcl)%write('odd_filt.mrc', iptcl)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_UNIFORM_FILTER2D NORMAL STOP ****')
    end subroutine exec_uniform_filter2D

    subroutine exec_icm2D( self, cline )
        class(icm2D_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        character(len=:), allocatable :: file_tag
        type(image),      allocatable :: odd(:), even(:)
        logical,          allocatable :: mask(:)
        type(parameters) :: params
        real             :: minmax(2)
        integer          :: iptcl
        ! init
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call find_ldim_nptcls(params%stk, params%ldim, params%nptcls)
        params%ldim(3) = 1 ! because we operate on stacks
        file_tag = 'icm_2D_filter'
        ! allocate
        allocate(odd(params%nptcls), even(params%nptcls), mask(params%nptcls))
        ! construct & read
        do iptcl = 1, params%nptcls
            call odd  (iptcl)%new( params%ldim, params%smpd, .false.)
            call odd  (iptcl)%read(params%stk,  iptcl)
            call even (iptcl)%new( params%ldim, params%smpd, .false.)
            call even (iptcl)%read(params%stk2, iptcl)
            minmax      = even(iptcl)%minmax()
            mask(iptcl) = .not.is_equal(minmax(2)-minmax(1),0.) ! empty image
        enddo
        ! filter
        !$omp parallel do schedule(static) default(shared) private(iptcl) proc_bind(close)
        do iptcl = 1, params%nptcls
            if( mask(iptcl) ) call even(iptcl)%ICM2D_eo(odd(iptcl), params%lambda)
        enddo
        !$omp end parallel do
        ! write output and destruct
        do iptcl = 1, params%nptcls
            call even (iptcl)%write(trim(file_tag)//'_even.mrc', iptcl)
            call odd  (iptcl)%write(trim(file_tag)//'_odd.mrc',  iptcl)
            call even (iptcl)%add(odd(iptcl))
            call even (iptcl)%mul(0.5)
            call even (iptcl)%write(trim(file_tag)//'_avg.mrc', iptcl)
            call odd  (iptcl)%kill()
            call even (iptcl)%kill()
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_ICM2D NORMAL STOP ****')
    end subroutine exec_icm2D

    subroutine exec_denoise_cavgs( self, cline )
        use simple_corrmat, only: calc_inplane_invariant_corrmat
        class(denoise_cavgs_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        character(len=:), allocatable :: file_tag
        type(image),      allocatable :: odd(:), even(:), avg(:)
        logical,          allocatable :: mask(:), mask_bincls(:)
        real,             allocatable :: fsc(:), filt(:), corrmat(:,:), res(:), corrs(:)
        real,             allocatable :: resolutions(:), resolutions_w(:), tmp(:)
        integer,          allocatable :: pinds(:), rank(:)
        real,    parameter :: HP_SCORE = 60.
        integer, parameter :: NN = 5
        type(image)        :: empty
        type(parameters)   :: params
        real               :: minmax(2), msk, res_frc0143, lp_score
        integer            :: iptcl, filtsz, n_nonempty, i, i_medoid, cnt_good, cnt_bad
        ! init
        call cline%set('ctf',    'no')
        call cline%set('objfun', 'cc')
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('lambda') ) call cline%set('lambda',  0.1)
        if( .not. cline%defined('trs')    ) call cline%set('trs',     5.0)
        call params%new(cline)
        call find_ldim_nptcls(params%stk, params%ldim, params%nptcls)
        params%box     = params%ldim(1)
        params%ldim(3) = 1 ! because we operate on stacks
        msk            = real(params%box / 2) - COSMSKHALFWIDTH - 1.
        filtsz         = fdim(params%box) - 1
        res            = get_resarr(params%box, params%smpd)
        file_tag       = 'cavgs_denoised'
        ! allocate
        allocate(odd(params%nptcls), even(params%nptcls), mask(params%nptcls), pinds(params%nptcls))
        do iptcl = 1, params%nptcls
            ! construct & read
            call odd (iptcl)%new( params%ldim, params%smpd, .false.)
            call odd (iptcl)%read(params%stk,  iptcl)
            call even(iptcl)%new( params%ldim, params%smpd, .false.)
            call even(iptcl)%read(params%stk2, iptcl)
            minmax       = even(iptcl)%minmax()
            mask(iptcl)  = .not.is_equal(minmax(2)-minmax(1),0.) ! empty image
            if( .not. mask(iptcl) )then
                call odd (iptcl)%kill
                call even(iptcl)%kill
            endif
            pinds(iptcl) = iptcl
        enddo
        n_nonempty = count(mask)
        pinds      = pack(pinds, mask=mask)
        even       = pack(even,  mask=mask)
        odd        = pack(odd,   mask=mask)
        ! filter
        allocate(fsc(filtsz), filt(filtsz), corrs(n_nonempty), resolutions(n_nonempty), tmp(n_nonempty), resolutions_w(n_nonempty), source=0.)
        !$omp parallel do schedule(static) default(shared) private(i,fsc,filt,res_frc0143) proc_bind(close)
        do i = 1, n_nonempty
            ! zero mean of outer pixels
            call odd (i)%zero_edgeavg
            call even(i)%zero_edgeavg
            ! spherical mask
            call odd (i)%mask(msk, 'soft', backgr=0.)
            call even(i)%mask(msk, 'soft', backgr=0.)
            ! fwd FFT
            call odd (i)%fft
            call even(i)%fft
            ! FSC
            call even(i)%fsc(odd(i), fsc)
            call get_resolution(fsc, res, resolutions(i), res_frc0143)
            ! calculate a filter to be applied to the individual e/o pairs
            where( fsc > 0. )
                filt = 2. * fsc / (fsc + 1.)   ! gold standard
                ! filt = fsc                     ! e/o merged
            else where
                filt = 0.
            end where
            where( filt  > 0.99999 ) filt = 0.99999
            ! apply filter
            call even(i)%apply_filter_serial(filt)
            call odd (i)%apply_filter_serial(filt)
            ! bwd FFT
            call odd (i)%ifft
            call even(i)%ifft
            ! ICM
            call even(i)%ICM2D_eo(odd(i), params%lambda)
        end do
        !$omp end parallel do
        ! write output
        do i = 1, n_nonempty
            iptcl = pinds(i)
            call even(i)%write(trim(file_tag)//'_even.mrc', iptcl)
            call odd (i)%write(trim(file_tag)//'_odd.mrc',  iptcl)
            call even(i)%add(odd(i))
            call even(i)%mul(0.5)
            call even(i)%write(trim(file_tag)//'_avg.mrc', iptcl)
        end do
        ! make sure the last image is written
        if( pinds(n_nonempty) /= params%nptcls )then
            call empty%new(params%ldim, params%smpd, .false.)
            call empty%write(trim(file_tag)//'_even.mrc', params%nptcls)
            call empty%write(trim(file_tag)//'_odd.mrc',  params%nptcls)
        endif
        ! rank according to resolution
        allocate(rank(n_nonempty), source=(/(i,i=1,n_nonempty)/))
        tmp = resolutions
        call hpsort(tmp, rank)
        ! write ranked
        do i = 1, n_nonempty
            call even(rank(i))%write(trim(file_tag)//'_ranked.mrc', i)
        end do
        ! calculate in-plane invariant similarity matrix
        lp_score = minval(resolutions)
        call calc_inplane_invariant_corrmat(even, HP_SCORE, lp_score, corrmat)
        ! estimate resolution as weighted average of five nearest neighbors
        do i = 1, n_nonempty
            corrs = corrmat(i,:)
            tmp   = resolutions
            call hpsort(corrs, tmp)
            call reverse(corrs)
            call reverse(tmp)
            ! resolutions_w(i) = sum(tmp(1:NN) * corrs(1:NN)) / sum(corrs(1:NN))
            resolutions_w(i) = sum(tmp(1:NN)) / real(NN)
        end do

        ! rank according to weigted resolution
        rank = (/(i,i=1,n_nonempty)/)
        tmp  = resolutions_w
        call hpsort(tmp, rank)
        ! write ranked
        do i = 1, n_nonempty
            call even(rank(i))%write(trim(file_tag)//'_reranked.mrc', i)

            print *, resolutions_w(rank(i))

        end do


        ! call medoid_ranking_from_smat(corrmat, i_medoid, rank)

        
        ! destruct
        do i = 1, n_nonempty
            call odd (i)%kill()
            call even(i)%kill()
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_DENOISE_CAVGS NORMAL STOP ****')
    end subroutine exec_denoise_cavgs

    subroutine exec_cavg_filter2D( self, cline )
        use simple_strategy2D3D_common, only: read_imgbatch, prepimgbatch, prepimg4align
        use simple_polarft_corrcalc,    only: polarft_corrcalc
        use simple_eul_prob_tab,        only: eul_prob_tab
        class(cavg_filter2D_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        complex,          allocatable :: cmat(:,:), cls_avg(:,:), ptcl_rot(:,:), ctf_rot(:,:)
        integer,          allocatable :: pinds(:)
        character(len=:), allocatable :: cavgsstk
        real,             allocatable :: denom(:,:)
        type(polarft_corrcalc)        :: pftcc
        type(builder)                 :: build
        type(parameters)              :: params
        type(image)                   :: img_cavg, calc_cavg, ptcl_match, img_cavg_pd
        integer  :: nptcls, iptcl, nptcls_cls
        logical  :: l_ctf
        integer  :: ncls, j, box, loc
        real     :: smpd
        call cline%set('dir_exec', 'cavg_filter2D')
        call cline%set('mkdir',    'yes')
        call cline%set('oritype',  'ptcl2D')
        call build%init_params_and_build_general_tbox(cline,params)
        ! reading all from the class 'class'
        call build%spproj%os_ptcl2D%get_pinds(params%class, 'class', pinds)
        nptcls     = build%spproj%get_nptcls()
        nptcls_cls = size(pinds)
        call pftcc%new(nptcls_cls, [1,nptcls_cls], params%kfromto)
        call pftcc%reallocate_ptcls(nptcls_cls, pinds)
        call build%img_crop_polarizer%init_polarizer(pftcc, params%alpha)
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
            call pftcc%rotate_ptcl(      pftcc%pfts_ptcls(:,:,j), loc, ptcl_rot)
            call pftcc%rotate_ptcl(cmplx(pftcc%ctfmats(:,:,j)),   loc, ctf_rot)
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
            call pftcc%polar2cartesian(cmplx(pftcc%ctfmats(:,:,j)), cmat, box)
            call calc_cavg%new([box,box,1], params%smpd*real(params%box)/real(box))
            call calc_cavg%zero_and_flag_ft
            call calc_cavg%set_cmat(cmat)
            call calc_cavg%shift_phorig()
            call calc_cavg%ifft
            call calc_cavg%write('ctfs_stk.mrc', j)
            ! writing the aligned ptcls stack
            call pftcc%polar2cartesian(ptcl_rot, cmat, box)
            call calc_cavg%new([box,box,1], params%smpd*real(params%box)/real(box))
            call calc_cavg%zero_and_flag_ft
            call calc_cavg%set_cmat(cmat)
            call calc_cavg%shift_phorig()
            call calc_cavg%ifft
            call calc_cavg%write('aligned_ptcls_stk.mrc', j)
            ! writing the aligned ctf stack
            call pftcc%polar2cartesian(ctf_rot, cmat, box)
            call calc_cavg%new([box,box,1], params%smpd*real(params%box)/real(box))
            call calc_cavg%zero_and_flag_ft
            call calc_cavg%set_cmat(cmat)
            call calc_cavg%shift_phorig()
            call calc_cavg%ifft
            call calc_cavg%write('aligned_ctfs_stk.mrc', j)
        enddo
        ! polar class average
        call pftcc%polar2cartesian(cls_avg / denom, cmat, box)
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
        ! taking in the corresponding denoised stk
        cls_avg = 0.
        denom   = 0.
        call img_cavg%new([box,box,1], params%smpd*real(params%box)/real(box))
        call img_cavg_pd%new([params%box,params%box,1], params%smpd)
        if( cline%defined('stk2') )then
            print *, nptcls_cls, box
            do j = 1, nptcls_cls
                iptcl = pinds(j)
                call img_cavg%read(trim(params%stk2), j)
                call img_cavg%fft
                call img_cavg_pd%zero_and_flag_ft
                call img_cavg%pad(img_cavg_pd)
                ! transfer to polar coordinates
                call build%img_crop_polarizer%polarize(pftcc, img_cavg_pd, iptcl, .true., .true.)
                ! e/o flags
                call pftcc%set_eo(iptcl, .true. )
                ! accumulating the cls_avg
                loc = pftcc%get_roind(build%spproj_field%e3get(iptcl))
                if( loc > pftcc%nrots ) loc = loc - pftcc%nrots
                call pftcc%rotate_ptcl(cmplx(pftcc%ctfmats(:,:,j)), loc, ctf_rot)
                cls_avg = cls_avg + ctf_rot * pftcc%pfts_ptcls(:,:,j)
                denom   = denom   + ctf_rot**2
            enddo
            ! polar class average
            call pftcc%polar2cartesian(cls_avg / denom, cmat, box)
            call calc_cavg%new([box,box,1], params%smpd*real(params%box)/real(box))
            call calc_cavg%zero_and_flag_ft
            call calc_cavg%set_cmat(cmat)
            call calc_cavg%shift_phorig()
            call calc_cavg%ifft
            call calc_cavg%write('polar_cavg_kPCA.mrc')
        endif
        call img_cavg_pd%kill
        call img_cavg%kill
        call calc_cavg%kill
        ! end gracefully
        call simple_end('**** SIMPLE_CAVG_FILTER2D NORMAL STOP ****')
    end subroutine exec_cavg_filter2D

    subroutine exec_score_ptcls( self, cline )
        use simple_strategy2D3D_common, only: discrete_read_imgbatch, prepimgbatch, prepimg4align, killimgbatch
        use simple_polarft_corrcalc,    only: polarft_corrcalc
        use simple_pftcc_shsrch_grad,   only: pftcc_shsrch_grad
        use simple_class_frcs,          only: class_frcs
        use simple_euclid_sigma2
        use simple_commander_euclid
        class(score_ptcls_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(pftcc_shsrch_grad), allocatable :: grad_shsrch_objs(:)
        type(image),             allocatable :: eimgs(:), oimgs(:), cls_even(:), cls_odd(:)
        type(calc_pspec_commander_distr) :: xcalc_pspec_distr
        type(polarft_corrcalc) :: pftcc
        type(builder)          :: build
        type(parameters)       :: params
        type(cmdline)          :: cline_calc_pspec_distr
        type(euclid_sigma2)    :: eucl_sigma
        type(class_frcs)       :: frcs
        real(sp),  allocatable :: scores(:,:), frc(:), filt(:), corrs(:)
        integer,   allocatable :: pinds(:), batches(:,:), cls(:)
        logical,   allocatable :: cls_mask(:)
        real     :: cxy(3), lims(2,2), lims_init(2,2), minmax(2), msk, best_corr,best_xy(2)
        integer  :: nptcls, iptcl, icls, batch_start, batch_end, filtsz, irot, inpl_ind,best_rot
        integer  :: ibatch, batchsz, ithr, i, j, nbatches, batchsz_max, funit, stat
        integer  :: best_class
        logical  :: l_ctf
        call cline%set('oritype',  'ptcl2D')
        call cline%set('mkdir',    'yes')
        call cline%set('objfun',   'euclid')
        call cline%set('part',     1)
        call cline%set('gridding', 'yes')
        call cline%set('ctf',      'yes')
        if( .not.cline%defined('trs') ) call cline%set('trs', MAXSHIFT)
        ! trs should be 0.1*box?
        call cline%delete('nparts')
        call build%init_params_and_build_general_tbox(cline, params)
        ! some init
        params%which_iter = 1
        call build%spproj%os_ptcl2D%set_all2single('w', 1.0)
        call build%spproj%write_segment_inside(params%oritype)
        ! Particles sampling
        call build%spproj_field%sample4update_all([params%fromp,params%top],nptcls,pinds,.false.)
        ! Number of classes
        call find_ldim_nptcls(params%stk, params%ldim, params%ncls)
        params%ldim(3) = 1
        call frcs%new(params%ncls, params%box, params%smpd, nstates=1)
        ! Batch dimensions
        batchsz_max = min(params%nptcls, 2*params%nthr*BATCHTHRSZ)
        nbatches    = ceiling(real(params%nptcls)/real(batchsz_max))
        batches     = split_nobjs_even(params%nptcls, nbatches)
        batchsz_max = maxval(batches(:,2)-batches(:,1)+1)
        call prepimgbatch(batchsz_max)
        ! Noise sigma2
        cline_calc_pspec_distr  = cline
        call cline_calc_pspec_distr%set('prg',   'calc_pspec' )
        call cline_calc_pspec_distr%set('mkdir', 'no')
        call xcalc_pspec_distr%execute_safe( cline_calc_pspec_distr )
        ! Read classes
        allocate(cls_odd(params%ncls), cls_even(params%ncls), cls_mask(params%nptcls),&
        &eimgs(params%nthr),oimgs(params%nthr))
        do icls = 1, params%ncls
            call cls_even(icls)%new( params%ldim, params%smpd, .false.)
            call cls_odd (icls)%new( params%ldim, params%smpd, .false.)
            call cls_even(icls)%read(params%stk2, icls)
            call cls_odd (icls)%read(params%stk,  icls)
            minmax         = cls_even(icls)%minmax()
            cls_mask(icls) = .not.is_equal(minmax(2)-minmax(1),0.)
        enddo
        do ithr = 1, params%nthr
            call eimgs(ithr)%new(params%ldim, params%smpd, .false.)
            call oimgs(ithr)%new(params%ldim, params%smpd, .false.)
        enddo
        ! Prep classes
        msk    = real(params%box/2-1)-COSMSKHALFWIDTH
        filtsz = fdim(params%box) - 1
        allocate(frc(filtsz), filt(filtsz), source=0.)
        !$omp parallel do schedule(static) default(shared) private(icls,frc,filt,ithr) proc_bind(close)
        do icls = 1, params%ncls
            ithr = omp_get_thread_num() + 1
            if( cls_mask(icls) )then
                call eimgs(ithr)%copy_fast(cls_even(icls))
                call oimgs(ithr)%copy_fast(cls_odd(icls))
                call eimgs(ithr)%mask(msk, 'soft', backgr=0.)
                call oimgs(ithr)%mask(msk, 'soft', backgr=0.)
                call eimgs(ithr)%fft
                call oimgs(ithr)%fft
                call eimgs(ithr)%fsc(oimgs(ithr), frc)
                call frcs%set_frc(icls, frc, 1)
                where( frc > 0. )
                    filt = 2. * frc / (frc + 1.)   ! gold standard
                    ! filt = frc                     ! e/o merged
                else where
                    filt = 0.
                end where
                where( filt  > 0.99999 ) filt = 0.99999
                call cls_even(icls)%fft
                call cls_odd(icls)%fft
                call cls_even(icls)%apply_filter_serial(filt)
                call cls_odd (icls)%apply_filter_serial(filt)
                call cls_even(icls)%ifft
                call cls_odd(icls)%ifft
                call cls_even(icls)%mask(params%msk, 'soft', backgr=0.)
                call cls_odd(icls)%mask( params%msk, 'soft', backgr=0.)
            endif
        enddo
        !$omp end parallel do
        ! Resolution limits
        if( .not. cline%defined('lp') )then
            params%lp = frcs%estimate_lp_for_align(state=1, crit0143=.false.)
        endif
        call frcs%write(FRCS_FILE)
        params%kfromto(1) = max(2,calc_fourier_index(params%hp, params%box, params%smpd))
        params%kfromto(2) = min(fdim(params%box)-1, calc_fourier_index(params%lp, params%box, params%smpd))
        write(logfhandle,'(A,F6.1)')'>>> RESOLUTION LIMIT(ANGS): ', params%lp
        ! PFTCC
        call pftcc%new(params%ncls, [1,batchsz_max], params%kfromto)
        call build%img_crop_polarizer%init_polarizer(pftcc, params%alpha)
        call eucl_sigma%new(SIGMA2_FBODY//int2str_pad(params%part,params%numlen)//'.dat', params%box)
        call eucl_sigma%read_part(  build%spproj_field)
        call eucl_sigma%read_groups(build%spproj_field)
        !$omp parallel do schedule(static) default(shared) private(icls,ithr) proc_bind(close)
        do icls = 1, params%ncls
            ithr = omp_get_thread_num() + 1
            if( cls_mask(icls) )then
                call build%img_crop_polarizer%div_by_instrfun(cls_even(icls))
                call build%img_crop_polarizer%div_by_instrfun(cls_odd(icls))
                call cls_even(icls)%fft
                call cls_odd(icls)%fft
                call build%img_crop_polarizer%polarize(pftcc, cls_even(icls), icls, .false., .true.,  build%l_resmsk)
                call build%img_crop_polarizer%polarize(pftcc, cls_odd(icls),  icls, .false., .false., build%l_resmsk)
            endif
            call cls_even(icls)%kill
            call cls_odd(icls)%kill
        enddo
        !$omp end parallel do
        do ithr = 1, params%nthr
            call oimgs(ithr)%kill
        enddo
        deallocate(cls_even,cls_odd,oimgs)
        call pftcc%memoize_refs
        ! CTF
        do i = 1,size(pinds)
            if( pinds(i)>0 )then
                iptcl = pinds(i)
                exit
            endif
        enddo
        l_ctf = build%spproj%get_ctfflag(params%oritype,iptcl=iptcl).ne.'no'
        ! Optimization allocations
        allocate(scores(params%ncls,nptcls),corrs(pftcc%get_nrots()),source=-1.)
        lims(:,1) = -params%trs
        lims(:,2) =  params%trs
        lims_init = lims / 2.
        allocate(grad_shsrch_objs(params%nthr))
        do ithr = 1, params%nthr
            call grad_shsrch_objs(ithr)%new(lims, lims_init=lims_init,&
                &shbarrier=params%shbarrier, maxits=60, opt_angle=.true., coarse_init=.true.)
        enddo
        ! Particles loop
        do ibatch=1,nbatches
            call progress_gfortran(ibatch, nbatches)
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            ! Refills pftcc
            call discrete_read_imgbatch(batchsz, pinds(batch_start:batch_end), [1,batchsz])
            call pftcc%reallocate_ptcls(batchsz, pinds(batch_start:batch_end))
            !$omp parallel do private(j,i,iptcl,ithr) default(shared) proc_bind(close)
            do j = batch_start,batch_end
                ithr  = omp_get_thread_num()+1
                i     = j - batch_start +  1
                iptcl = pinds(j)
                call eimgs(ithr)%zero_and_flag_ft
                call prepimg4align(iptcl, build%imgbatch(i), eimgs(ithr))
                call build%img_crop_polarizer%polarize(pftcc, eimgs(ithr), iptcl, .true., .true.)
                call pftcc%set_eo(iptcl, (build%spproj_field%get_eo(iptcl)==0))
            enddo
            !$omp end parallel do
            if( l_ctf ) call pftcc%create_polar_absctfmats(build%spproj, params%oritype)
            call pftcc%memoize_ptcls
            ! Scoring
            !$omp parallel do private(j,i,iptcl,icls,ithr,cxy,irot,inpl_ind,corrs,best_class,best_corr,best_xy,best_rot)&
            !$omp schedule(dynamic) proc_bind(close) default(shared)
            do j = batch_start,batch_end
                ithr       = omp_get_thread_num()+1
                i          = j - batch_start +  1
                iptcl      = pinds(j)
                best_corr  = -1.
                best_class = 0
                best_xy    = 0.
                best_rot   = 0
                do icls = 1,params%ncls
                    if( .not. cls_mask(icls) ) cycle
                    call pftcc%gencorrs(icls, iptcl, corrs)
                    irot     = maxloc(corrs, dim=1)
                    inpl_ind = irot
                    call grad_shsrch_objs(ithr)%set_indices(icls, iptcl)
                    cxy = grad_shsrch_objs(ithr)%minimize(irot=inpl_ind)
                    if( inpl_ind == 0 )then
                        inpl_ind = irot
                        cxy      = [real(pftcc%gencorr_for_rot_8(icls, iptcl, irot)), 0.,0.]
                    endif
                    scores(icls,j) = cxy(1)
                    if( cxy(1) > best_corr )then
                        best_corr  = cxy(1)
                        best_class = icls
                        best_xy    = cxy(2:3)
                        best_rot   = inpl_ind
                    endif
                enddo
                call build%spproj_field%e3set(iptcl, 360.-pftcc%get_rot(best_rot))
                call build%spproj_field%set_shift(iptcl,    best_xy)
                call build%spproj_field%set(iptcl, 'inpl',  best_rot)
                call build%spproj_field%set(iptcl, 'class', best_class)
                call build%spproj_field%set(iptcl, 'corr',  best_corr)
            enddo
            !$omp end parallel do
        enddo
        call killimgbatch
        call pftcc%kill
        call build%spproj%write(params%projfile)
        ! Write array
        call fopen(funit, 'scores.mat', form='UNFORMATTED', iostat=stat)
        write(funit) [params%ncls,nptcls]
        write(funit) pinds
        if( l_ctf) write(funit) build%spproj_field%get_all('dfx')
        cls = nint(build%spproj_field%get_all('class'))
        write(funit) cls
        write(funit) scores
        call fclose(funit)
        ! from scipy.io import FortranFile
        ! import numpy as np
        ! f     = FortranFile('scores.mat', 'r')
        ! dims  = f.read_reals(dtype=np.int32)
        ! pinds = f.read_reals(dtype=np.int32)
        ! dfx   = f.read_reals(dtype=np.float32)
        ! pcls  = f.read_reals(dtype=np.int32)
        ! A     = f.read_reals(dtype=np.float32).reshape(dims, order='F')
        ! f.close()
        call simple_end('**** SIMPLE_score_ptcls NORMAL STOP ****')
    end subroutine exec_score_ptcls

    subroutine exec_estimate_lpstages( self, cline )
        class(estimate_lpstages_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(builder)    :: build
        type(parameters) :: params
        real, parameter  :: LPSTART_LB=10., LPSTART_DEFAULT=20.
        real, parameter  :: LPSTOP_BOUNDS(2)  = [4.5,6.0]
        real, parameter  :: LPSTART_BOUNDS(2) = [10.,20.] 
        character(len=:),  allocatable :: frcs_fname
        real,              allocatable :: frcs_avg(:)
        integer,           allocatable :: states(:)
        type(lp_crop_inf), allocatable :: lpinfo(:)
        integer :: filtsz
        real    :: lpfinal
        call build%init_params_and_build_general_tbox(cline, params)
        call build%spproj%get_frcs(frcs_fname, 'frc2D', fail=.true.)
        call build%clsfrcs%read(frcs_fname)
        states  = nint(build%spproj%os_cls2D%get_all('state'))
        filtsz = build%clsfrcs%get_filtsz()
        allocate(frcs_avg(filtsz), source=0.)
        if( params%frc_weight .eq. 'yes' )then
            call build%clsfrcs%avg_frc_getter(frcs_avg, states, cur_oris=build%spproj%os_ptcl2D)
        else
            call build%clsfrcs%avg_frc_getter(frcs_avg, states)
        endif
        allocate(lpinfo(params%nstages))
        lpfinal = max(LPSTOP_BOUNDS(1),calc_lplim_final_stage(3))
        lpfinal = min(LPSTOP_BOUNDS(2),lpfinal)
        call lpstages(params%box, params%nstages, frcs_avg, params%smpd,&
        &LPSTART_BOUNDS(1), LPSTART_BOUNDS(2), lpfinal, lpinfo, l_cavgs=.false.)
        call simple_end('**** SIMPLE_ESTIMATE_LPSTAGES NORMAL STOP ****')

        contains

            function calc_lplim_final_stage( nbest ) result( lplim )
                integer, intent(in)  :: nbest
                real,    allocatable :: res(:), tmp_rarr(:)
                integer, allocatable :: states(:), tmp_iarr(:)
                real :: lplim
                tmp_rarr  = build%spproj%os_cls2D%get_all('res')
                tmp_iarr  = nint(build%spproj%os_cls2D%get_all('state'))
                res       = pack(tmp_rarr, mask=(tmp_iarr>0))
                call hpsort(res)
                lplim = median_nocopy(res(:nbest))
                deallocate(tmp_rarr, tmp_iarr, res)
            end function calc_lplim_final_stage

    end subroutine exec_estimate_lpstages

end module simple_commander_resolest
