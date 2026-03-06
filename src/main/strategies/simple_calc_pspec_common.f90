module simple_calc_pspec_common
use simple_core_module_api
use simple_image,               only: image
use simple_sigma2_binfile,      only: sigma2_binfile
use simple_cmdline,             only: cmdline
use simple_builder,             only: builder
use simple_parameters,          only: parameters
use simple_sp_project,          only: sp_project
use simple_qsys_env,            only: qsys_env
use simple_qsys_funs,           only: qsys_job_finished
use simple_strategy2D3D_common, only: prepimgbatch, discrete_read_imgbatch, killimgbatch
use simple_qsys_funs,           only: qsys_job_finished
implicit none

public :: calc_pspec_exec
private

contains

    !> Core power spectrum calculation logic.
    !> This subroutine is called by the in-memory strategy and by each distributed worker.
    subroutine calc_pspec_exec( params, build )
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        type(image)              :: sum_img
        type(sigma2_binfile)     :: binfile
        type(string)             :: binfname
        complex(dp), allocatable :: cmat_thr_sum(:,:,:)
        complex,     allocatable :: cmat_sum(:,:,:)
        integer,     allocatable :: pinds(:)
        real,        allocatable :: pspec(:), sigma2(:,:)
        integer :: batchlims(2),kfromto(2)
        integer :: i,iptcl,imatch,nyq,nptcls_part_sel,batchsz_max,nbatch
        logical :: l_scale_update_frac
        ! Sampling
        ! Because this is always run prior to reconstruction/search, sampling is not always informed
        ! or may change with workflows. Instead of setting a sampling for the following operations when
        ! l_update_frac, we sample uniformly AND do not write the corresponding field
        l_scale_update_frac = .false.
        if( params%l_update_frac )then
            call build%spproj_field%sample4update_rnd([params%fromp,params%top], params%update_frac, nptcls_part_sel, pinds, .false. )
            l_scale_update_frac = .true.
        else
            call build%spproj_field%sample4update_all([params%fromp,params%top], nptcls_part_sel, pinds, .false.)
        endif
        ! init
        nyq = build%img%get_nyq()
        allocate(sigma2(nyq,params%fromp:params%top),pspec(nyq),source=0.)
        batchsz_max = min(nptcls_part_sel,50 * nthr_glob)
        call prepimgbatch(params, build, batchsz_max)
        call sum_img%new([params%box,params%box,1],params%smpd)
        call sum_img%zero_and_flag_ft
        cmat_sum = sum_img%allocate_cmat()
        allocate(cmat_thr_sum(size(cmat_sum,dim=1),size(cmat_sum,dim=2),1))
        ! mask memoization
        call build%imgbatch(1)%memoize_mask_coords
        do i = 1,nptcls_part_sel,batchsz_max
            batchlims = [i, min(i+batchsz_max-1,nptcls_part_sel)]
            nbatch    = batchlims(2) - batchlims(1) + 1
            call discrete_read_imgbatch(params, build, nbatch, pinds(batchlims(1):batchlims(2)), [1,nbatch])
            cmat_thr_sum = dcmplx(0.d0,0.d0)
            !$omp parallel do default(shared) private(iptcl,imatch,pspec)&
            !$omp schedule(static) proc_bind(close) reduction(+:cmat_thr_sum)
            do imatch = 1,nbatch
                iptcl = pinds(batchlims(1)+imatch-1)
                call build%imgbatch(imatch)%norm_noise_mask_fft_powspec(build%lmsk, params%msk, pspec)
                if( l_scale_update_frac )then
                    ! To account for spectra not included in sampling and yield the correct average
                    sigma2(:,iptcl) = pspec / (2.0 * params%update_frac)
                else
                    sigma2(:,iptcl) = pspec / 2.0
                endif
                ! thread average
                call build%imgbatch(imatch)%add_dble_cmat2mat(cmat_thr_sum(:,:,:))
            end do
            !$omp end parallel do
            ! global average
            cmat_sum(:,:,:) = cmat_sum(:,:,:) + cmplx(cmat_thr_sum(:,:,:),kind=sp)
        end do
        call sum_img%set_cmat(cmat_sum)
        call sum_img%write(string('sum_img_part')//int2str_pad(params%part,params%numlen)//params%ext%to_char())
        ! write to disk
        kfromto  = [1, nyq]
        binfname = 'init_pspec_part'//trim(int2str(params%part))//'.dat'
        call binfile%new(binfname,params%fromp,params%top,kfromto)
        call binfile%write(sigma2)
        ! destruct
        call binfile%kill
        call killimgbatch(build)
        call sum_img%kill
    end subroutine calc_pspec_exec

end module simple_calc_pspec_common
