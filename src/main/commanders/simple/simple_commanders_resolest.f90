!@descr: for resolution estimation
module simple_commanders_resolest
use simple_commanders_api
use simple_pftc_srch_api
use simple_fsc
use simple_commanders_euclid, only: commander_calc_pspec
use simple_refine3D_fnames,  only: refine3D_fsc_fbody, refine3D_fsc_fname
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_fsc
  contains
    procedure :: execute      => exec_fsc
end type commander_fsc

type, extends(commander_base) :: commander_fsc_area_score
  contains
    procedure :: execute      => exec_fsc_area_score
end type commander_fsc_area_score

type, extends(commander_base) :: commander_uniform_filter2D
  contains
    procedure :: execute      => exec_uniform_filter2D
end type commander_uniform_filter2D

type, extends(commander_base) :: commander_uniform_filter3D
  contains
    procedure :: execute      => exec_uniform_filter3D
end type commander_uniform_filter3D

type, extends(commander_base) :: commander_nu_filt3D
    contains
        procedure :: execute      => exec_nu_filt3D
end type commander_nu_filt3D

type, extends(commander_base) :: commander_nu_filt2D
    contains
        procedure :: execute      => exec_nu_filt2D
end type commander_nu_filt2D

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
        if( params%automsk .ne. 'no' )then
            call mskvol%new([params%box,params%box,params%box], params%smpd)
            call mskvol%automask3D(params, even, odd, l_tight=params%automsk.eq.'tight')
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
        fsc_templ = refine3D_fsc_fbody(1)
        call arr2file(fsc, refine3D_fsc_fname(1))
        if( params%automsk .ne. 'no' )then
            call plot_phrand_fsc(size(fsc), fsc, fsc_t, fsc_n, res, params%smpd, fsc_templ%to_char())
        else
            call plot_fsc(size(fsc), fsc, res, params%smpd, fsc_templ%to_char())
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_FSC NORMAL STOP ****')
    end subroutine exec_fsc

    !> calculates a CryoSPARC-like conical FSC area ratio from Even/Odd Volume pairs
    subroutine exec_fsc_area_score( self, cline )
        class(commander_fsc_area_score), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)             :: params
        type(image)                  :: even, odd
        type(image_msk)              :: mskvol
        type(fsc_area_score_result)  :: result
        type(string)                 :: fbody
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('athres')     ) call cline%set('athres', 20.)
        if( .not. cline%defined('nspace')     ) call cline%set('nspace', 256)
        if( .not. cline%defined('lplim_crit') ) call cline%set('lplim_crit', 0.143)
        call params%new(cline)
        call odd%new([params%box,params%box,params%box], params%smpd)
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd%read(params%vols(1))
        call even%read(params%vols(2))
        if( params%automsk .ne. 'no' )then
            call mskvol%new([params%box,params%box,params%box], params%smpd)
            call mskvol%automask3D(params, even, odd, l_tight=params%automsk.eq.'tight')
            call even%zero_env_background(mskvol)
            call odd%zero_env_background(mskvol)
            call even%mul(mskvol)
            call odd%mul(mskvol)
        else
            call even%mask3D_soft(params%msk)
            call odd%mask3D_soft(params%msk)
        endif
        call even%fft()
        call odd%fft()
        call calc_fsc_area_score(even, odd, params%nspace, params%athres, params%lplim_crit, 1, result)
        if( cline%defined('fbody') )then
            fbody = params%fbody
        else
            fbody = 'fsc_area_score'
        endif
        call write_fsc_area_score_outputs(result, fbody%to_char())
        call plot_fsc_area_score(result, fbody%to_char())
        call even%kill
        call odd%kill
        call mskvol%kill
        call simple_end('**** SIMPLE_FSC_AREA_SCORE NORMAL STOP ****')
    end subroutine exec_fsc_area_score

    subroutine exec_uniform_filter3D(self, cline)
        use simple_opt_filter, only: estimate_lplim
        class(commander_uniform_filter3D), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters)  :: params
        type(image)       :: even, odd, odd_filt
        type(image_msk)   :: mskvol
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
        if( params%automsk .ne. 'no' )then
            call mskvol%new([params%box,params%box,params%box], params%smpd)
            call mskvol%automask3D(params, even, odd, l_tight=params%automsk.eq.'tight')
            call mskvol%one_at_edge ! to expand before masking of reference internally
        else
            call mskvol%disc([params%box,params%box,params%box], params%smpd,&
                    &real(min(params%box/2, int(params%msk + COSMSKHALFWIDTH))))
        endif
        ! soft masking needed for FT
        call even%mask3D_soft(params%msk)
        call  odd%mask3D_soft(params%msk)
        call estimate_lplim(odd, even, mskvol, params%lpstart, params%lpstop, lpopt, odd_filt)
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

    subroutine exec_nu_filt3D(self, cline)
        use simple_nu_filter, only: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, cleanup_nu_filter, &
            &print_nu_filtmap_lowpass_stats, analyze_filtmap_neighbor_continuity, write_nu_local_resolution_map
        class(commander_nu_filt3D), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters)     :: params
        type(image)          :: even, odd, even_nu, odd_nu, vol_msk
        type(image_msk)      :: envmsk
        type(string)         :: even_out, odd_out, avg_out, locres_out
        logical, allocatable :: l_mask(:,:,:)
        integer, allocatable :: imat(:,:,:)
        real                 :: mskrad_px
        integer              :: ldim(3)
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call odd%new([params%box,params%box,params%box], params%smpd)
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd%read(params%vols(1))
        call even%read(params%vols(2))
        ldim = even%get_ldim()
        if( params%automsk .ne. 'no' )then
            call envmsk%automask3D(params, even, odd, l_tight=params%automsk.eq.'tight')
            call envmsk%set_imat
            call envmsk%get_imat(imat)
            allocate(l_mask(ldim(1),ldim(2),ldim(3)))
            l_mask = imat > 0
            deallocate(imat)
        else
            mskrad_px = 0.5 * params%mskdiam / params%smpd
            call vol_msk%disc(ldim, params%smpd, mskrad_px, l_mask)
        endif
        call setup_nu_dmats(even, odd, l_mask, [real ::])
        if( allocated(l_mask) ) deallocate(l_mask)
        call optimize_nu_cutoff_finds(histogram_potts=params%l_nu_hist_potts)
        call nu_filter_vols(even_nu, odd_nu)
        call print_nu_filtmap_lowpass_stats()
        call analyze_filtmap_neighbor_continuity()
        odd_out  = add2fbody(params%vols(1), params%ext, NUFILT_SUFFIX)
        even_out = add2fbody(params%vols(2), params%ext, NUFILT_SUFFIX)
        if( params%outvol .ne. '' )then
            avg_out = params%outvol
        else
            avg_out = 'vol_nufilt'//params%ext%to_char()
        endif
        locres_out = add2fbody(avg_out, params%ext, NULOCRES_SUFFIX)
        call odd_nu%write(odd_out, del_if_exists=.true.)
        call even_nu%write(even_out, del_if_exists=.true.)
        call even_nu%add(odd_nu)
        call even_nu%mul(0.5)
        call even_nu%write(avg_out, del_if_exists=.true.)
        call write_nu_local_resolution_map(locres_out)
        call wait_for_closure(avg_out)
        call wait_for_closure(locres_out)
        call cleanup_nu_filter()
        call odd_nu%kill
        call even_nu%kill
        call odd%kill
        call even%kill
        call vol_msk%kill
        call envmsk%kill
        call locres_out%kill
        if( allocated(l_mask) ) deallocate(l_mask)
        call simple_end('**** SIMPLE_nu_filt3D NORMAL STOP ****')
    end subroutine exec_nu_filt3D

    subroutine exec_nu_filt2D(self, cline)
        use simple_nu_filter2D, only: nu_filter2D_state
        use simple_nu_filter2D_stats, only: nu_filter2D_stats, merge_nu_filter2D_stats, print_nu_filter2D_stats, &
            &kill_nu_filter2D_stats
        class(commander_nu_filt2D), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters)     :: params
        type(nu_filter2D_state), allocatable :: nu_states(:)
        type(nu_filter2D_stats) :: nu_stats
        type(string)         :: odd_out, even_out, avg_out, locres_out
        integer              :: ldim(3), ldim2(3), nimgs, nimgs2, iptcl, ithr
        real                 :: align_lp
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('outstk') ) call cline%set('outstk', 'stk_nufilt'//MRC_EXT)
        call params%new(cline)
        call find_ldim_nptcls(params%stk,  ldim,  nimgs)
        call find_ldim_nptcls(params%stk2, ldim2, nimgs2)
        if( any(ldim /= ldim2) ) THROW_HARD('odd/even stacks must have identical dimensions; nu_filt2D')
        if( nimgs /= nimgs2 ) THROW_HARD('odd/even stacks must contain the same number of images; nu_filt2D')
        if( ldim(3) /= 1 ) THROW_HARD('nu_filt2D expects 2D image stacks')
        odd_out    = add2fbody(params%stk,  params%ext, NUFILT_SUFFIX)
        even_out   = add2fbody(params%stk2, params%ext, NUFILT_SUFFIX)
        avg_out    = params%outstk
        locres_out = add2fbody(avg_out, params%ext, NULOCRES_SUFFIX)
        if( file_exists(odd_out)    ) call del_file(odd_out)
        if( file_exists(even_out)   ) call del_file(even_out)
        if( file_exists(avg_out)    ) call del_file(avg_out)
        if( file_exists(locres_out) ) call del_file(locres_out)
        write(logfhandle,'(A)') '>>> 2D nonuniform filter: low-pass bank 30,20,15,12,8,6,5,4 A; whole image'
        write(logfhandle,'(A)') '>>> 2D nonuniform filter: standalone mode uses the full discrete bank'
        write(logfhandle,'(A)') '>>> 2D nonuniform filter: objective smoothing radius = AWF * LP, capped at 30 A'
        write(logfhandle,'(A)') '>>> 2D nonuniform filter: all bank members compete directly'
        write(logfhandle,'(A)') '>>> 2D nonuniform filter: Potts ICM weights one- and two-pixel neighborhoods'
        write(logfhandle,'(A)') '>>> 2D nonuniform filter: Potts penalty includes weak 1-label jumps'
        write(logfhandle,'(A)') '>>> 2D nonuniform filter: output blends bank members over a 10 A tent-smoothed field'
        write(logfhandle,'(A)') '>>> 2D nonuniform filter: standalone mode without auxiliary competitor'
        allocate(nu_states(nthr_glob))
        do ithr = 1, nthr_glob
            call nu_states(ithr)%setup(ldim, params%smpd)
        end do
        !$omp parallel do default(shared) private(iptcl,ithr,align_lp) schedule(static) ordered proc_bind(close)
        do iptcl = 1, nimgs
            ithr = omp_get_thread_num() + 1
            block
                type(image) :: odd, even, avg_raw, odd_nu, even_nu, avg_nu, locres
                call odd%new(ldim, params%smpd, .false.)
                call even%new(ldim, params%smpd, .false.)
                call odd%read(params%stk, iptcl)
                call even%read(params%stk2, iptcl)
                call avg_raw%copy(even)
                call avg_raw%add(odd)
                call avg_raw%mul(0.5)
                call nu_states(ithr)%apply(even, odd, avg_raw, avg_raw, avg_raw, even_nu, odd_nu, avg_nu, &
                    &align_lp, locres_out=locres)
                !$omp ordered
                call odd_nu%write(odd_out, iptcl)
                call even_nu%write(even_out, iptcl)
                call avg_nu%write(avg_out, iptcl)
                call locres%write(locres_out, iptcl)
                !$omp end ordered
                call odd%kill
                call even%kill
                call avg_raw%kill
                call odd_nu%kill
                call even_nu%kill
                call avg_nu%kill
                call locres%kill
            end block
        end do
        !$omp end parallel do
        call update_stack_nimgs(odd_out, nimgs)
        call update_stack_nimgs(even_out, nimgs)
        call update_stack_nimgs(avg_out, nimgs)
        call update_stack_nimgs(locres_out, nimgs)
        do ithr = 1, nthr_glob
            call merge_nu_filter2D_stats(nu_stats, nu_states(ithr)%stats)
        end do
        call print_nu_filter2D_stats(nu_stats)
        write(logfhandle,'(A,1X,A)') '>>> Wrote 2D NU odd stack:',       odd_out%to_char()
        write(logfhandle,'(A,1X,A)') '>>> Wrote 2D NU even stack:',      even_out%to_char()
        write(logfhandle,'(A,1X,A)') '>>> Wrote 2D NU average stack:',   avg_out%to_char()
        write(logfhandle,'(A,1X,A)') '>>> Wrote 2D NU local-res stack:', locres_out%to_char()
        do ithr = 1, nthr_glob
            call nu_states(ithr)%kill
        end do
        call kill_nu_filter2D_stats(nu_stats)
        deallocate(nu_states)
        call odd_out%kill
        call even_out%kill
        call avg_out%kill
        call locres_out%kill
        call simple_end('**** SIMPLE_nu_filt2D NORMAL STOP ****')
    end subroutine exec_nu_filt2D

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
            call envmsk%automask3D(params, even, odd, l_tight=params%automsk.eq.'tight')
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
        call estimate_lplims2D( odd, even, params%msk, params%lpstart, params%lpstop, lpsopt, odd_filt )
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
