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
implicit none

public :: fsc_commander
! public :: local_res_commander
public :: nonuniform_filter_commander
public :: butterworth_3D_commander
public :: butterworth_2D_commander
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

type, extends(commander_base) :: butterworth_2D_commander
  contains
    procedure :: execute      => exec_butterworth_2D
end type butterworth_2D_commander

type, extends(commander_base) :: butterworth_3D_commander
  contains
    procedure :: execute      => exec_butterworth_3D
end type butterworth_3D_commander


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

    subroutine exec_butterworth(self, cline, dim3)
        use simple_butterworth, only: butterworth_kernel
        
        class(commander_base), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        integer,               intent(in)    :: dim3            ! dim3 is 1 for 2D case
        integer,               parameter     :: CHUNKSZ=20
        type(parameters) :: params
        type(image)      :: even, odd, odd_img, ker_odd_img
        integer          :: k,l,m,n,k1,l1,m1,k_ind,l_ind,m_ind, max_sup
        real             :: rad, cur_min_sum, ref_diff

        real   , parameter   :: A = 47.27, B = -0.1781, C = 7.69, D = -0.02228  ! Fitting constants (constructed in MATLAB) of theta(FT_support) = a*exp(b*x) + c*exp(d*x)
        real   , parameter   :: MIN_SUP = 0.5, RES_LB = 30                      ! lower bound of resolution is 30 Angstrom, upper bound is nyquist, hence .5
        integer, parameter   :: N_SUP = 20                                      ! number of intervals between MIN_SUP and MAX_SUP
        integer, parameter   :: SPA_SUP = 2, MID = 1+SPA_SUP                    ! support of the window function
        real   , allocatable :: weights(:,:,:)                                  ! weights of the neighboring differences

        ! optimization variables
        real                         , allocatable :: cur_mat(:,:,:), sup, theta
        real                         , pointer     :: rmat_odd(:,:,:), rmat_even(:,:,:), rmat_ker(:,:,:), orig_ker(:,:,:), orig_ker_der(:,:,:),  prev_diff(:,:,:), cur_diff(:,:,:)
        complex(kind=c_float_complex), pointer     :: cmat_odd(:,:,:), cmat_conv(:,:,:)


        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        
        call odd %new([params%box,params%box,dim3], params%smpd)
        call even%new([params%box,params%box,dim3], params%smpd)
        call odd %read(params%vols(1))
        call even%read(params%vols(2))

        call odd_img    %new([params%box,params%box,dim3], params%smpd)
        call ker_odd_img%new([params%box,params%box,dim3], params%smpd)

        allocate(rmat_odd(     params%box,params%box,dim3))
        allocate(rmat_even(    params%box,params%box,dim3))
        allocate(rmat_ker(     params%box,params%box,dim3))
        allocate(orig_ker(     params%box,params%box,dim3))
        allocate(orig_ker_der( params%box,params%box,dim3))
        allocate(prev_diff(    params%box,params%box,dim3))
        allocate(cur_diff(     params%box,params%box,dim3))
        allocate(cur_mat(      params%box,params%box,dim3))

        allocate(cmat_conv(    int(params%box/2)+1,params%box,dim3))
        allocate(cmat_odd(     int(params%box/2)+1,params%box,dim3))

        call odd %get_rmat_sub(rmat_odd)
        call even%get_rmat_sub(rmat_even)
        call odd_img%set_rmat(rmat_odd, .false.)
        call odd_img%fft()
        call odd_img%get_cmat_ptr(cmat_odd)

        rmat_odd     = rmat_odd /sum(rmat_odd)      ! Normalize to energy of 1, so B*odd is comparible with even in the cost function
        rmat_even    = rmat_even/sum(rmat_even)
        cur_mat      = 0.
        max_sup      = int(RES_LB/params%smpd)*2    ! multiplication factor depending on the definition of support, set to 2 for now
        prev_diff    = 3.4028235E+38 - 1            ! supposed to be infinity
        cur_min_sum  = 3.4028235E+38 - 1
        
        ! assign the weights of the neighboring voxels
        allocate(weights(SPA_SUP*2+1, SPA_SUP*2+1, SPA_SUP*2+1))
        do k = 1, 2*SPA_SUP+1
            do l = 1, 2*SPA_SUP+1
                do m = 1, 2*SPA_SUP+1
                    rad = hyp(real(k-MID), real(l-MID), real(m-MID))
                    weights(k,l,m) = -rad/(SPA_SUP + 1) + 1.  ! linear function: 1 at rad = 0 and 0 at rad = SPA_SUP + 1
                    if (weights(k,l,m) < 0.) then
                        weights(k,l,m) = 0.
                    endif
                enddo
            enddo
        enddo
        weights = weights/sum(weights) ! weights has energy of 1

        do n = 1, N_SUP
            sup   = MIN_SUP + (n-1.)*(max_sup-MIN_SUP)/(N_SUP-1.)
            theta = A*exp(B*sup) + C*exp(D*sup)
            write(*, *) 'support = ', sup, '; theta = ', theta
            call butterworth_kernel(orig_ker, orig_ker_der, params%box, 8, theta)  ! WARNING: fix the constants here

            ! computing B_kernel 'convolve' odd
            call ker_odd_img%set_rmat(orig_ker, .false.)
            call ker_odd_img%fft()
            call ker_odd_img%get_cmat_ptr(cmat_conv)
            cmat_conv = cmat_conv * cmat_odd
            call ker_odd_img%ifft()
            call ker_odd_img%get_rmat_sub(rmat_ker)

            ! Normalize to energy of 1, so B*odd is comparible with even in the cost function
            rmat_ker = rmat_ker/sum(rmat_ker)

            cur_diff = (rmat_ker - rmat_even)**2

            ! do the non-uniform, i.e. optimizing at each voxel
            if (params%is_uniform == 'no') then
                !$omp parallel do collapse(3) default(shared) private(k,l,m,k1,l1,m1,k_ind,l_ind,m_ind,ref_diff) schedule(dynamic,CHUNKSZ) proc_bind(close)
                do k = 1,params%box
                    do l = 1,params%box
                        do m = 1,dim3
                            ref_diff = 0.
                            ! applying an average window to each diff (eq 7 in the nonuniform paper)
                            do k_ind = 1, 2*SPA_SUP+1
                                k1 = k - SPA_SUP + k_ind - 1
                                do l_ind = 1, 2*SPA_SUP+1
                                    l1 = l - SPA_SUP + l_ind - 1

                                    ! 2D vs 3D cases
                                    if (dim3 == 1) then
                                        if ((k1 >= 1 .and. k1 <= params%box) .and. (l1 >= 1 .and. l1 <= params%box)) then
                                            ref_diff = ref_diff + cur_diff(k1,l1,m)
                                        endif
                                    else
                                        do m_ind = 1, 2*SPA_SUP+1
                                            m1 = m - SPA_SUP + m_ind - 1
                                            if ((k1 >= 1 .and. k1 <= params%box) .and. (l1 >= 1 .and. l1 <= params%box) .and. (m1 >= 1 .and. m1 <= dim3)) then
                                                ref_diff = ref_diff + cur_diff(k1,l1,m1)*weights(k_ind,l_ind,m_ind)
                                            endif
                                        enddo
                                    endif
                                enddo
                            enddo

                            ! prev_diff keeps the lowest cost value at each voxel of the search
                            ! cur_mat   keeps the best voxel of the form B*odd
                            if (ref_diff < prev_diff(k,l,m)) then
                                cur_mat(k,l,m)   = rmat_ker(k,l,m)
                                prev_diff(k,l,m) = ref_diff
                            endif
                        enddo
                    enddo
                enddo
                !$omp end parallel do
            else
                ! keep the theta which gives the lowest cost (over all voxels)
                if (sum(cur_diff) < cur_min_sum) then
                    cur_mat     = rmat_ker
                    cur_min_sum = sum(cur_diff)
                endif
            endif
            
            call odd%set_rmat(cur_mat, .false.)
            call odd %write(params%is_uniform // '_uniform_butterworth_filter_iter' // int2str(n) // '.mrc')
        enddo
        
        ! end gracefully
        call simple_end('**** SIMPLE_BUTTERWORTH_FILTER NORMAL STOP ****')
    end subroutine exec_butterworth


    subroutine exec_butterworth_2D( self, cline )
        class(butterworth_2D_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        call exec_butterworth(self, cline, 1)
    end subroutine exec_butterworth_2D

    subroutine exec_butterworth_3D( self, cline )
        class(butterworth_3D_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters) :: params
        call params%new(cline)
        call exec_butterworth(self, cline, params%box)
    end subroutine exec_butterworth_3D

end module simple_commander_resolest
