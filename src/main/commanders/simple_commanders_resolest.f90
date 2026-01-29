!@descr: for resolution estimation
module simple_commanders_resolest
use simple_commander_module_api
use simple_fsc
use simple_commanders_euclid, only: commander_calc_pspec_distr
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_fsc
  contains
    procedure :: execute      => exec_fsc
end type commander_fsc

type, extends(commander_base) :: commander_clin_fsc
  contains
    procedure :: execute      => exec_clin_fsc
end type commander_clin_fsc

type, extends(commander_base) :: commander_uniform_filter2D
  contains
    procedure :: execute      => exec_uniform_filter2D
end type commander_uniform_filter2D

type, extends(commander_base) :: commander_nununiform_filter3D
  contains
    procedure :: execute      => exec_nununiform_filter3D
end type commander_nununiform_filter3D

type, extends(commander_base) :: commander_uniform_filter3D
  contains
    procedure :: execute      => exec_uniform_filter3D
end type commander_uniform_filter3D

type, extends(commander_base) :: commander_icm3D
  contains
    procedure :: execute      => exec_icm3D
end type commander_icm3D

type, extends(commander_base) :: commander_icm2D
  contains
    procedure :: execute      => exec_icm2D
end type commander_icm2D

type, extends(commander_base) :: commander_estimate_lpstages
  contains
    procedure :: execute      => exec_estimate_lpstages
end type commander_estimate_lpstages

contains

    !> calculates Fourier shell correlation from Even/Odd Volume pairs
    subroutine exec_fsc( self, cline )
        class(commander_fsc), intent(inout) :: self
        class(cmdline),       intent(inout) :: cline
        real, allocatable :: res(:), fsc(:), fsc_t(:), fsc_n(:)
        type(parameters)  :: params
        type(image)       :: even, odd
        type(image_msk)   :: mskvol
        type(string)      :: fsc_templ
        integer           :: j, k_hp, k_lp, nyq
        real              :: res_fsc05, res_fsc0143
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
            call even%mask3D_soft(params%msk)
            call odd%mask3D_soft(params%msk)
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
        fsc_templ = FSC_FBODY//int2str_pad(1,2)
        call arr2file(fsc, fsc_templ//BIN_EXT)
        if( params%l_filemsk )then
            call plot_phrand_fsc(size(fsc), fsc, fsc_t, fsc_n, res, params%smpd, fsc_templ%to_char())
        else
            call plot_fsc(size(fsc), fsc, res, params%smpd, fsc_templ%to_char())
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_FSC NORMAL STOP ****')
    end subroutine exec_fsc

    subroutine exec_clin_fsc( self, cline )
        use simple_strategy2D3D_common
        use simple_polarops
        use simple_polarft_calc, only: polarft_calc
        use simple_imgarr_utils, only: write_imgarr
        class(commander_clin_fsc), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        integer,          allocatable :: pinds(:)
        type(image),      allocatable :: tmp_imgs(:), cavgs(:)
        type(polarft_calc)            :: pftc
        type(builder)                 :: build
        type(parameters)              :: params
        integer :: nptcls, ithr
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
        call set_bp_range( cline )
        call build%spproj_field%sample4update_all([params%fromp,params%top], nptcls, pinds, incr_sampled=.false.)
        ! PREPARATION OF PARTICLES
        call pftc%new(params%nspace, [1,nptcls], params%kfromto)
        call build%img_crop%memoize4polarize(pftc%get_pdim(), params%alpha)
        call prepimgbatch(nptcls)
        allocate(tmp_imgs(nthr_glob))
        !$omp parallel do default(shared) private(ithr) schedule(static) proc_bind(close)
        do ithr = 1,nthr_glob
            call tmp_imgs(ithr)%new([params%box_crop,params%box_crop,1], params%smpd_crop, wthreads=.false.)
        enddo
        ! Build polar particle images
        call pftc%allocate_refs_memoization
        call build_batch_particles(pftc, nptcls, pinds, tmp_imgs)
        ! Dealing with polar cavgs
        call polar_cavger_new(pftc, .true.)
        call polar_cavger_update_sums(nptcls, pinds, build%spproj, pftc, is3D=.true.)
        call polar_cavger_merge_eos_and_norm(reforis=build%eulspace)
        ! write
        allocate(cavgs(params%nspace))
        call polar_cavger_write(string('cavgs_even.bin'), 'even')
        call polar_cavger_write(string('cavgs_odd.bin'),   'odd')
        call polar_cavger_write(string('cavgs.bin'),    'merged')
        call polar_cavger_refs2cartesian(pftc, cavgs,   'even')
        call write_imgarr(cavgs, string('cavgs_even.mrc'))
        call polar_cavger_refs2cartesian(pftc, cavgs,    'odd')
        call write_imgarr(cavgs, string('cavgs_odd.mrc'))
        call polar_cavger_refs2cartesian(pftc, cavgs, 'merged')
        call write_imgarr(cavgs, string('cavgs_merged.mrc'))
        call polar_cavger_kill
        call killimgbatch
        call pftc%kill
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_CLIN_FSC NORMAL STOP ****')
    end subroutine exec_clin_fsc

    subroutine exec_nununiform_filter3D(self, cline)
        use simple_opt_filter, only: nonuni_filt3D
        class(commander_nununiform_filter3D), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters)  :: params
        type(image)       :: even, odd, mskvol
        logical           :: have_mask_file
        type(string)      :: file_tag
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
            call even%mask3D_soft(params%msk)
            call odd%mask3D_soft(params%msk)
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
            call even%mask3D_soft(params%msk)
            call odd%mask3D_soft(params%msk)
        endif
        call odd%write(file_tag//'_odd.mrc')
        call even%write(file_tag//'_even.mrc')
        call odd%add(even)
        call odd%mul(0.5)
        call odd%write(file_tag//'_avg.mrc')
        ! destruct
        call odd%kill
        call even%kill
        call mskvol%kill
        ! end gracefully
        call simple_end('**** SIMPLE_NONUNIFORM_FILTER3D NORMAL STOP ****')
    end subroutine exec_nununiform_filter3D

    subroutine exec_uniform_filter3D(self, cline)
        use simple_opt_filter, only: estimate_lplim
        class(commander_uniform_filter3D), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters)  :: params
        type(image)       :: even, odd, mskvol, odd_filt
        real              :: lpopt
        type(string)      :: file_tag
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
        call even%mask3D_soft(params%msk)
        call  odd%mask3D_soft(params%msk)
        call estimate_lplim(odd, even, mskvol, [params%lpstart,params%lpstop], lpopt, odd_filt)
        print *, 'found optimal low-pass limit: ', lpopt
        call odd_filt%write(string('odd_filt.mrc'))
        ! destruct
        call odd%kill
        call odd_filt%kill
        call even%kill
        call mskvol%kill
        ! end gracefully
        call simple_end('**** SIMPLE_UNIFORM_FILTER3D NORMAL STOP ****')
    end subroutine exec_uniform_filter3D

    subroutine exec_icm3D( self, cline )
        class(commander_icm3D), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(parameters)     :: params
        type(image)          :: even, odd, even_icm, odd_icm, avg, avg_icm
        type(image_msk)         :: envmsk
        logical, allocatable :: l_msk(:,:,:)
        type(string) :: file_tag
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
        call even_icm%write(string('vol2_')//file_tag//'.mrc')
        call odd_icm%write(string('vol1_')//file_tag//'.mrc')
        call avg_icm%write(string('vol_avg_')//file_tag//'.mrc')
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
        class(commander_uniform_filter2D), intent(inout) :: self
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
            call odd_filt(iptcl)%write(string('odd_filt.mrc'), iptcl)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_UNIFORM_FILTER2D NORMAL STOP ****')
    end subroutine exec_uniform_filter2D

    subroutine exec_icm2D( self, cline )
        class(commander_icm2D), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(image), allocatable :: odd(:), even(:)
        logical,     allocatable :: mask(:)
        type(string)     :: file_tag
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
            call even (iptcl)%write(file_tag//'_even.mrc', iptcl)
            call odd  (iptcl)%write(file_tag//'_odd.mrc',  iptcl)
            call even (iptcl)%add(odd(iptcl))
            call even (iptcl)%mul(0.5)
            call even (iptcl)%write(file_tag//'_avg.mrc', iptcl)
            call odd  (iptcl)%kill()
            call even (iptcl)%kill()
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_ICM2D NORMAL STOP ****')
    end subroutine exec_icm2D

    subroutine exec_estimate_lpstages( self, cline )
        class(commander_estimate_lpstages), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(builder)    :: build
        type(parameters) :: params
        real, parameter  :: LPSTART_LB=10., LPSTART_DEFAULT=20.
        real, parameter  :: LPSTOP_BOUNDS(2)  = [4.5,6.0]
        real, parameter  :: LPSTART_BOUNDS(2) = [10.,20.] 
        type(string)     :: frcs_fname
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
        call build%clsfrcs%avg_frc_getter(frcs_avg, states)
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
                integer, allocatable :: tmp_iarr(:)
                real :: lplim
                tmp_rarr  = build%spproj%os_cls2D%get_all('res')
                tmp_iarr  = nint(build%spproj%os_cls2D%get_all('state'))
                res       = pack(tmp_rarr, mask=(tmp_iarr>0))
                call hpsort(res)
                lplim = median_nocopy(res(:nbest))
                deallocate(tmp_rarr, tmp_iarr, res)
            end function calc_lplim_final_stage

    end subroutine exec_estimate_lpstages

end module simple_commanders_resolest
