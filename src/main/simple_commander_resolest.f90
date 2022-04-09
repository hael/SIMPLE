! concrete commander: resolest for resolution estimation
module simple_commander_resolest
include 'simple_lib.f08'
use simple_parameters,     only: parameters, params_glob
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_image,          only: image
use simple_masker,         only: masker
implicit none

public :: fsc_commander
! public :: local_res_commander
public :: nonuniform_filter_commander
public :: nonuniform_butterworth_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: fsc_commander
  contains
    procedure :: execute      => exec_fsc
end type fsc_commander
! type, extends(commander_base) :: local_res_commander
!   contains
!     procedure :: execute      => exec_local_res
! end type local_res_commander
type, extends(commander_base) :: nonuniform_filter_commander
  contains
    procedure :: execute      => exec_nonuniform_filter
end type nonuniform_filter_commander

type, extends(commander_base) :: nonuniform_butterworth_commander
  contains
    procedure :: execute      => exec_nonuniform_butterworth
end type nonuniform_butterworth_commander

contains

    !> calculates Fourier shell correlation from Even/Odd Volume pairs
    subroutine exec_fsc( self, cline )
        use simple_estimate_ssnr, only: fsc2optlp_sub, fsc2TVfilt
        class(fsc_commander), intent(inout) :: self
        class(cmdline),       intent(inout) :: cline
        type(parameters)  :: params
        type(image)       :: even, odd, e_copy
        type(masker)      :: mskvol
        integer           :: j, find_plate, k_hp, k_lp, nyq, flims(3,2)
        real              :: res_fsc05, res_fsc0143
        real, allocatable :: res(:), corrs(:), filt_optlp(:), filt_tv(:)
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        ! read even/odd pair
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd%new([params%box,params%box,params%box], params%smpd)
        call odd%read(params%vols(1))
        call even%read(params%vols(2))
        if( cline%defined('mskfile') )then
            if( file_exists(params%mskfile) )then
                call mskvol%new([params%box,params%box,params%box], params%smpd)
                call mskvol%read(params%mskfile)
                call even%zero_background
                call odd%zero_background
                call even%mul(mskvol)
                call odd%mul(mskvol)
            else
                THROW_HARD('mskfile: '//trim(params%mskfile)//' does not exist in cwd; exec_fsc')
            endif
        else
            ! spherical masking
            call even%mask(params%msk, 'soft')
            call odd%mask(params%msk, 'soft')
        endif
        ! forward FT
        call even%fft()
        call odd%fft()
        ! calculate FSC
        res = even%get_res()
        nyq = even%get_filtsz()
        allocate(corrs(nyq), filt_optlp(nyq), filt_tv(nyq))
        call even%fsc(odd, corrs)
        if( params%l_phaseplate ) call phaseplate_correct_fsc(corrs, find_plate)
        do j=1,nyq
           write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
        end do
        call get_resolution(corrs, res, res_fsc05, res_fsc0143)
        k_hp = calc_fourier_index(params%hp, params%box, params%smpd)
        k_lp = calc_fourier_index(params%lp, params%box, params%smpd)
        write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
        write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
        write(logfhandle,'(A,1X,F8.4)') '>>> MEDIAN FSC (SPECSCORE):', median_nocopy(corrs(k_hp:k_lp))
        flims = even%loop_lims(2)
        call fsc2optlp_sub(nyq, corrs, filt_optlp)
        call fsc2TVfilt(corrs, flims, filt_tv)
        call e_copy%copy(even)
        call even%apply_filter(filt_optlp)
        call e_copy%apply_filter(filt_tv)
        call even%ifft
        call e_copy%ifft
        call even%write('filtered_optlp.mrc')
        call e_copy%write('filtered_tv.mrc')
        do j=1,nyq
            write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> OPTLP:', filt_optlp(j), '>>> TV:', filt_tv(j)
        end do
        call even%kill
        call odd%kill
        ! end gracefully
        call simple_end('**** SIMPLE_FSC NORMAL STOP ****')
    end subroutine exec_fsc

    !> calculates local resolution from Even/Odd Volume pairs
    ! subroutine exec_local_res( self, cline )
    !     use simple_estimate_ssnr, only: local_res, local_res_lp
    !     class(local_res_commander), intent(inout) :: self
    !     class(cmdline),             intent(inout) :: cline
    !     type(parameters)     :: params
    !     type(image)          :: even, odd, vol2filter
    !     type(masker)         :: mskvol
    !     real                 :: res_fsc05, res_fsc0143
    !     logical              :: have_mask_file
    !     integer, allocatable :: locres_finds(:,:,:)
    !     real,    allocatable :: res(:), corrs(:)
    !     if( .not. cline%defined('mkdir') )      call cline%set('mkdir', 'yes')
    !     call params%new(cline)
    !     ! read even/odd pair
    !     call even%new([params%box,params%box,params%box], params%smpd)
    !     call odd%new([params%box,params%box,params%box], params%smpd)
    !     call odd%read(params%vols(1))
    !     call even%read(params%vols(2))
    !     have_mask_file = .false.
    !     if( cline%defined('mskfile') )then
    !         if( file_exists(params%mskfile) )then
    !             call mskvol%new([params%box,params%box,params%box], params%smpd)
    !             call mskvol%read(params%mskfile)
    !             call even%zero_background
    !             call odd%zero_background
    !             call even%mul(mskvol)
    !             call odd%mul(mskvol)
    !             have_mask_file = .true.
    !         else
    !             THROW_HARD('mskfile: '//trim(params%mskfile)//' does not exist in cwd; exec_local_res')
    !         endif
    !     else
    !         ! spherical masking
    !         call even%mask(params%msk, 'soft')
    !         call odd%mask(params%msk, 'soft')
    !     endif
    !     ! forward FT
    !     call even%fft()
    !     call odd%fft()
    !     if( cline%defined('mskfile') )then
    !         call mskvol%read(params%mskfile)
    !         call mskvol%remove_edge
    !     else
    !         call mskvol%disc([params%box,params%box,params%box], params%smpd, params%msk)
    !     endif
    !     call local_res(even, odd, mskvol, params%lplim_crit, locres_finds)
    !     ! destruct
    !     call even%kill
    !     call odd%kill
    !     ! filter inputted vol3
    !     if( cline%defined('vol3') )then
    !         call vol2filter%new([params%box,params%box,params%box], params%smpd)
    !         call vol2filter%read(params%vols(1))
    !         call local_res_lp(locres_finds, vol2filter)
    !         call vol2filter%write('map_locres_lp.mrc')
    !         call vol2filter%kill
    !     endif
    !     ! end gracefully
    !     call simple_end('**** SIMPLE_LOCAL_RES NORMAL STOP ****')
    ! end subroutine exec_local_res

    subroutine exec_nonuniform_filter( self, cline )
        use simple_estimate_ssnr, only: nonuniform_fsc_filt
        class(nonuniform_filter_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters) :: params
        type(image)      :: even, odd, map2filt
        type(masker)     :: mskvol
        logical          :: have_mask_file, map2filt_present
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        ! read even/odd pair
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd%new([params%box,params%box,params%box],  params%smpd)
        call odd%read(params%vols(1))
        call even%read(params%vols(2))
        map2filt_present = cline%defined('vol3')
        if( map2filt_present )then
            call map2filt%new([params%box,params%box,params%box],  params%smpd)
            call map2filt%read(params%vols(3))
        endif
        have_mask_file = .false.
        if( cline%defined('mskfile') )then
            if( file_exists(params%mskfile) )then
                call mskvol%new([params%box,params%box,params%box], params%smpd)
                call mskvol%read(params%mskfile)
                have_mask_file = .true.
            else
                THROW_HARD('mskfile: '//trim(params%mskfile)//' does not exist in cwd; exec_nonuniform_filter')
            endif
        else
            ! spherical masking
            call even%mask(params%msk, 'soft')
            call odd%mask(params%msk, 'soft')
        endif
        if( have_mask_file )then
            call mskvol%one_at_edge ! to expand before masking of reference internally
        else
            call mskvol%disc([params%box,params%box,params%box], params%smpd, params%msk)
        endif
        if( map2filt_present )then
            call nonuniform_fsc_filt(even, odd, mskvol, .false., map2filt, phran=trim(params%phrand).eq.'yes')
        else
            call nonuniform_fsc_filt(even, odd, mskvol, .false., phran=trim(params%phrand).eq.'yes')
        endif
        if( have_mask_file )then
            call mskvol%read(params%mskfile)
            call even%mul(mskvol)
            call odd%mul(mskvol)
            if( map2filt_present ) call map2filt%mul(mskvol)
            call mskvol%kill
        endif
        if( map2filt_present )then
            call map2filt%write('nonuniformly_filtered.mrc')
        else
            call even%write('nonuniformly_filtered_even.mrc')
            call odd%write('nonuniformly_filtered_odd.mrc')
        endif
        ! destruct
        call even%kill
        call odd%kill
        if( map2filt_present ) call map2filt%kill
        ! end gracefully
        call simple_end('**** SIMPLE_NONUNIFORM_FILTER NORMAL STOP ****')
    end subroutine exec_nonuniform_filter

    subroutine exec_nonuniform_butterworth( self, cline )
        use simple_butterworth
        use simple_optimizer,   only: optimizer
        use simple_opt_factory, only: opt_factory
        use simple_opt_spec,    only: opt_spec
        class(nonuniform_butterworth_commander), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        type(parameters) :: params
        type(image)      :: even, odd, map2opt
        type(masker)     :: mskvol
        logical          :: have_mask_file, map2opt_present

        procedure(fun_butterworth),     pointer   :: costfun_ptr      !< pointer 2 cost function
        class(optimizer),               pointer   :: opt_ptr=>null()  ! the generic optimizer object
        character(len=*),               parameter ::  odd_img_name  = 'NoisyObj1.mrc'
        character(len=*),               parameter :: even_img_name  = 'NoisyObj2.mrc'
        integer,                        parameter :: box    = 202
        real,                           parameter :: smpd   = 1.275
        integer,                        parameter :: ndim   = 1, NRESTARTS = 1
        type(opt_factory)   :: ofac                           ! the optimization factory object
        type(opt_spec)      :: spec                           ! the optimizer specification object
        character(len=8)    :: str_opts                       ! string descriptors for the NOPTS optimizers
        real                :: lims(ndim,2), lowest_cost
    
        call even_img%new([box,box,1], smpd)
        call even_img%read(even_img_name, 1)

        call odd_img%new([box,box,1], smpd)
        call odd_img%read(odd_img_name, 1)

        call ker_odd_img%new([box,box,1], smpd)
        call ker_even_img%new([box,box,1], smpd)
        call ker_der_odd_img%new([box,box,1], smpd)
        call ker_der_even_img%new([box,box,1], smpd)

        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        ! read even/odd pair
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd %new([params%box,params%box,params%box], params%smpd)
        call odd %read(params%vols(1))
        call even%read(params%vols(2))
        map2opt_present = cline%defined('vol3')
        if( map2opt_present )then
            call map2opt%new([params%box,params%box,params%box],  params%smpd)
            call map2opt%read(params%vols(3))
        endif
        have_mask_file = .false.
        if( cline%defined('mskfile') )then
            if( file_exists(params%mskfile) )then
                call mskvol%new([params%box,params%box,params%box], params%smpd)
                call mskvol%read(params%mskfile)
                have_mask_file = .true.
            else
                THROW_HARD('mskfile: '//trim(params%mskfile)//' does not exist in cwd; exec_nonuniform_butterworth')
            endif
        else
            ! spherical masking
            call even%mask(params%msk, 'soft')
            call odd %mask(params%msk, 'soft')
        endif
        if( have_mask_file )then
            call mskvol%one_at_edge ! to expand before masking of reference internally
        else
            call mskvol%disc([params%box,params%box,params%box], params%smpd, params%msk)
        endif
        if( map2opt_present )then
            write(*, *) 'TODO'
        else
            ! do the optimization here
            write(*, *) 'Simple cut-off frequency optimization test:'
            costfun_ptr  => butterworth_cost
            str_opts  = 'lbfgsb'
            lims(1,1) = -50.
            lims(1,2) =  50.
            call spec%specify(str_opts, ndim, limits=lims, nrestarts=NRESTARTS) ! make optimizer spec
            call spec%set_costfun(costfun_ptr)                                  ! set pointer to costfun
            call spec%set_gcostfun(butterworth_gcost)                           ! set pointer to gradient of costfun         
            call ofac%new(spec, opt_ptr)                                        ! generate optimizer object with the factory
            spec%x    = 7.                                                     ! set initial guess
            call opt_ptr%minimize(spec, opt_ptr, lowest_cost)                   ! minimize the test function

            write(*, *) 'cost = ', lowest_cost, '; x = ', spec%x

            call opt_ptr%kill
            deallocate(opt_ptr)
        endif
        if( have_mask_file )then
            call mskvol%read(params%mskfile)
            call even%mul(mskvol)
            call odd %mul(mskvol)
            if( map2opt_present ) call map2opt%mul(mskvol)
            call mskvol%kill
        endif
        if( map2opt_present )then
            call map2opt%write('nonuniformly_butterworth.mrc')
        else
            call even%write('nonuniformly_butterworth_even.mrc')
            call odd %write('nonuniformly_butterworth_odd.mrc')
        endif
        ! destruct
        call even%kill
        call odd %kill
        if( map2opt_present ) call map2opt%kill
        ! end gracefully
        call simple_end('**** SIMPLE_NONUNIFORM_BUTTERWORTH NORMAL STOP ****')
    end subroutine exec_nonuniform_butterworth

end module simple_commander_resolest
