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
    procedure :: execute      => exec_nonuniform_butterworth_discrete
end type nonuniform_butterworth_commander

contains

    !> calculates Fourier shell correlation from Even/Odd Volume pairs
    subroutine exec_fsc( self, cline )
        use simple_estimate_ssnr, only: fsc2optlp_sub, fsc2TVfilt
        class(fsc_commander), intent(inout) :: self
        class(cmdline),       intent(inout) :: cline
        type(parameters)  :: params
        type(image)       :: even, odd
        type(masker)      :: mskvol
        integer           :: j, find_plate, k_hp, k_lp, nyq
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

    subroutine exec_nonuniform_butterworth_continuous( self, cline )
        use simple_butterworth
        use simple_optimizer,   only: optimizer
        use simple_opt_factory, only: opt_factory
        use simple_opt_spec,    only: opt_spec
        class(nonuniform_butterworth_commander), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        type(parameters) :: params
        type(image)      :: map2opt
        type(masker)     :: mskvol
        logical          :: have_mask_file, map2opt_present

        procedure(fun_butterworth),     pointer   :: costfun_ptr      !< pointer 2 cost function
        class(optimizer),               pointer   :: opt_ptr=>null()  ! the generic optimizer object
        integer,                        parameter :: ndim   = 1, NRESTARTS = 1
        type(opt_factory)   :: ofac                           ! the optimization factory object
        type(opt_spec)      :: spec                           ! the optimizer specification object
        character(len=8)    :: str_opts                       ! string descriptors for the NOPTS optimizers
        real                :: lims(ndim,2), lowest_cost

        real                         , pointer :: rmat_odd(:,:,:), rmat_even(:,:,:), rmat_ker(:,:,:),  orig_ker(:,:,:), orig_ker_der(:,:,:)
        complex(kind=c_float_complex), pointer :: cmat_odd(:,:,:), cmat_even(:,:,:), cmat_conv(:,:,:)

        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        ! read even/odd pair
        call even_img%new([params%box,params%box,params%box], params%smpd)
        call odd_img %new([params%box,params%box,params%box], params%smpd)
        call odd_img %read(params%vols(1))
        call even_img%read(params%vols(2))
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
                THROW_HARD('mskfile: '//trim(params%mskfile)//' does not exist in cwd; exec_nonuniform_butterworth_continuous')
            endif
        else
            ! spherical masking [TODO] does it make sense to mask when doing butterworth optimization?
            !call even_img%mask(params%msk, 'soft')
            !call odd_img %mask(params%msk, 'soft')
        endif
        if( have_mask_file )then
            call mskvol%one_at_edge ! to expand before masking of reference internally
        else
            call mskvol%disc([params%box,params%box,params%box], params%smpd, params%msk)
        endif
        if( map2opt_present )then
            write(*, *) 'TODO'
        else
            ! initialize the butterworth kernel and its derivative
            call ker_odd_img     %new([params%box,params%box,params%box], params%smpd)
            call ker_even_img    %new([params%box,params%box,params%box], params%smpd)
            call ker_der_odd_img %new([params%box,params%box,params%box], params%smpd)
            call ker_der_even_img%new([params%box,params%box,params%box], params%smpd)

            ! initialize the temporary arrays
            allocate(rmat_odd(     params%box,params%box,params%box))
            allocate(rmat_even(    params%box,params%box,params%box))
            allocate(rmat_ker(     params%box,params%box,params%box))
            allocate(orig_ker(     params%box,params%box,params%box))
            allocate(orig_ker_der( params%box,params%box,params%box))
            allocate(cmat_odd(     params%box,params%box,params%box))
            allocate(cmat_even(    params%box,params%box,params%box))

            ! do the optimization here to get the optimized cut-off frequency
            write(*, *) 'Cut-off frequency optimization in progress:'
            costfun_ptr  => butterworth_cost
            str_opts  = 'lbfgsb'
            lims(1,1) =  0.00001
            lims(1,2) =  3.
            call spec%specify(str_opts, ndim, limits=lims, nrestarts=NRESTARTS) ! make optimizer spec
            call spec%set_costfun(costfun_ptr)                                  ! set pointer to costfun
            call spec%set_gcostfun(butterworth_gcost)                           ! set pointer to gradient of costfun
            call ofac%new(spec, opt_ptr)                                        ! generate optimizer object with the factory
            spec%x = (lims(1,1) + lims(1,2))/2.                                 ! set initial guess
            call opt_ptr%minimize(spec, opt_ptr, lowest_cost)                   ! minimize the test function

            write(*, *) 'cost = ', lowest_cost, '; x = ', spec%x(1)

            ! apply the optimized cut-off frequency to the odd and even image
            call butterworth_kernel(orig_ker, orig_ker_der, params%box, 8, spec%x(1))  ! WARNING: fix the constants here
            
            ! do the convolution B*even
            call even_img%get_rmat_ptr(rmat_even)
            call even_img%fft()
            call even_img%get_cmat_ptr(cmat_even)

            call ker_odd_img%get_rmat_ptr(rmat_ker)
            call ker_odd_img%set_rmat(orig_ker, .false.)
            call ker_odd_img%fft()
            call ker_odd_img%get_cmat_ptr(cmat_conv)
            cmat_even = cmat_even*cmat_conv
            call even_img%ifft()

            ! do the convolution B*odd
            call odd_img%get_rmat_ptr(rmat_odd)
            call odd_img%fft()
            call odd_img%get_cmat_ptr(cmat_odd)
            cmat_odd = cmat_odd*cmat_conv
            call odd_img%ifft()

            call opt_ptr%kill
            deallocate(opt_ptr)
        endif
        if( have_mask_file )then
            call mskvol%read(params%mskfile)
            call even_img%mul(mskvol)
            call odd_img %mul(mskvol)
            if( map2opt_present ) call map2opt%mul(mskvol)
            call mskvol%kill
        endif
        if( map2opt_present )then
            call map2opt%write('nonuniformly_butterworth.mrc')
        else
            call even_img%write('nonuniformly_butterworth_continuous_even.mrc')
            call odd_img %write('nonuniformly_butterworth_continuous_odd.mrc')
        endif
        ! destruct
        call even_img%kill
        call odd_img %kill
        if( map2opt_present ) call map2opt%kill
        ! end gracefully
        call simple_end('**** SIMPLE_NONUNIFORM_BUTTERWORTH_CONTINUOUS NORMAL STOP ****')
    end subroutine exec_nonuniform_butterworth_continuous

    subroutine exec_nonuniform_butterworth_discrete( self, cline )
        use simple_butterworth, only: butterworth_kernel
        
        class(nonuniform_butterworth_commander), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        type(parameters) :: params
        type(image)      :: even, odd, odd_img, ker_odd_img, map2opt
        type(masker)     :: mskvol
        logical          :: have_mask_file, map2opt_present
        integer          :: k,l,m,n

        ! Fitting constants (constructed in MATLAB) in theta(FT_support) = a*exp(b*x) + c*exp(d*x)
        real   , parameter :: a = 47.27, b = -0.1781, c = 7.69, d = -0.02228, min_sup = 2., max_sup = 100.
        integer, parameter :: N_sup = 50

        ! optimization variables
        real                         , allocatable :: cur_mat(:,:,:), sup, theta
        real                         , pointer     :: rmat_odd(:,:,:), rmat_even(:,:,:), rmat_ker(:,:,:), orig_ker_der(:,:,:),  prev_diff(:,:,:),  cur_diff(:,:,:)
        complex(kind=c_float_complex), pointer     :: cmat_odd(:,:,:), cmat_conv(:,:,:)


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
                THROW_HARD('mskfile: '//trim(params%mskfile)//' does not exist in cwd; exec_nonuniform_butterworth_discrete')
            endif
        else
            ! spherical masking [TODO] comment this for now
            !call even%mask(params%msk, 'soft')
            !call odd %mask(params%msk, 'soft')
        endif
        if( have_mask_file )then
            call mskvol%one_at_edge ! to expand before masking of reference internally
        else
            call mskvol%disc([params%box,params%box,params%box], params%smpd, params%msk)
        endif
        if( map2opt_present )then
            write(*, *) 'TODO'
        else
            call odd_img    %new([params%box,params%box,params%box], params%smpd)
            call ker_odd_img%new([params%box,params%box,params%box], params%smpd)

            allocate(rmat_odd(     params%box,params%box,params%box))
            allocate(rmat_even(    params%box,params%box,params%box))
            allocate(rmat_ker(     params%box,params%box,params%box))
            allocate(orig_ker_der( params%box,params%box,params%box))
            allocate(cmat_odd(     params%box,params%box,params%box))
            allocate(prev_diff(    params%box,params%box,params%box))
            allocate(cur_diff(     params%box,params%box,params%box))
            allocate(cur_mat(      params%box,params%box,params%box))
        
            call odd %get_rmat_sub(rmat_odd)
            call even%get_rmat_sub(rmat_even)
            call odd_img%set_rmat(rmat_odd, .false.)
            call odd_img%fft()
            call odd_img%get_cmat_ptr(cmat_odd)

            cur_mat = 0.
            do n = 1, N_sup
                sup   = min_sup + (n-1.)*(max_sup-min_sup)/(N_sup-1.)
                theta = a*exp(b*sup) + c*exp(d*sup)
                write(*, *) sup, theta
                call butterworth_kernel(rmat_ker, orig_ker_der, params%box, 8, theta)  ! WARNING: fix the constants here
                call ker_odd_img%set_rmat(rmat_ker, .false.)
                call ker_odd_img%fft()
                call ker_odd_img%get_cmat_ptr(cmat_conv)
                cmat_conv = cmat_conv * cmat_odd
                call ker_odd_img%ifft()
                if (sup > min_sup) then
                    cur_diff  = abs(rmat_ker - rmat_even)
                    do k = 1,params%box
                        do l = 1,params%box
                            do m = 1,params%box
                                if (cur_diff(k,l,m) < prev_diff(k,l,m)) then
                                    cur_mat(k,l,m) = rmat_ker(k,l,m)
                                endif
                            enddo
                        enddo
                    enddo
                    prev_diff = cur_diff
                else
                    cur_mat   = rmat_ker
                    prev_diff = abs(rmat_ker - rmat_even)
                endif
            enddo
            call odd_img%kill
            call ker_odd_img%kill
            call odd%set_rmat(cur_mat, .false.)
        endif
        if( have_mask_file )then
            call mskvol%read(params%mskfile)
            call even%mul(mskvol)
            call odd %mul(mskvol)
            if( map2opt_present ) call map2opt%mul(mskvol)
            call mskvol%kill
        endif
        if( map2opt_present )then
            call map2opt%write('nonuniformly_butterworth_discrete.mrc')
        else
            call even%write('nonuniformly_butterworth_discrete_even.mrc')
            call odd %write('nonuniformly_butterworth_discrete_odd.mrc')
        endif
        ! destruct
        call even%kill
        call odd %kill
        if( map2opt_present ) call map2opt%kill
        ! end gracefully
        call simple_end('**** SIMPLE_NONUNIFORM_BUTTERWORTH_DISCRETE NORMAL STOP ****')
    end subroutine exec_nonuniform_butterworth_discrete

end module simple_commander_resolest
