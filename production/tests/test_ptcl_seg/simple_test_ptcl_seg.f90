program ptcl_seg
include 'simple_lib.f08'
    use simple_image
    use simple_spectral_clustering
    use simple_optimizer
    use simple_segmentation
    use simple_aff_prop
    use simple_afm_image
    use simple_gauss2Dfit
    use simple_math
    use simple_motion_patched
    implicit none 
    type(image)         :: ini_conf, fin_conf, gau2D_test, gau2D_sum
    type(image), allocatable    :: gau2D_stk(:), clus_stk(:), defm_stack(:)
    type(spec_clust)    :: spc_ptcl
    type(aff_prop)      :: aff_prop_ptcl
    character(len = 255)    :: im_directory = '/Users/atifao/Downloads/test_for_ptcl_seg.mrc'
    character(len = 255)    :: stk_directory = '/Users/atifao/Downloads/ordered_picks.mrc'
    real                :: smpd = 5.0, sum, a, cen_gauss(2), cov_gauss(2,2), corr
    integer             :: nptcls, ldim(3), i, j, k, num_clus, iter, counts, label, l, indx, n
    real, allocatable       :: sim_mat(:,:), rmat(:,:,:), bin_mat(:,:), rmat_new(:,:,:), theta(:), r1(:), r2(:), cen_mat(:,:,:), data1(:), means(:)
    integer, allocatable    :: labels(:), centers(:), new_centers(:), center_x(:), center_y(:)
    logical     :: converged
    logical, allocatable    :: mask_labels(:)
    type(motion_patched)    :: defm_model
    CHARACTER(len=:), ALLOCATABLE   :: fname
    type(parameters), target :: params
    ! may need to downsample
    call find_ldim_nptcls(im_directory,ldim,nptcls)
    call ini_conf%new(ldim,smpd)
    call ini_conf%read(im_directory,1)
    allocate(rmat(ldim(1),ldim(2),1))
    allocate(cen_mat(ldim(1),ldim(2),1))
    rmat = ini_conf%get_rmat()
    k = 0
    ! # of nonzero px
    do i = 1, ldim(1)
        do j = 1, ldim(2)
            if(rmat(i,j,1) > 0.1) k = k + 1
        end do 
    end do 
    allocate(data1(k))
    counts = 0
    do i = 1, ldim(1)
        do j = 1, ldim(2)
            if(rmat(i,j,1) > 0.1 .and. counts .le. k) then 
                data1(counts) = sqrt( real(i - 1)**2 + real(j - 1)**2)
                counts = counts + 1
            end if
        end do 
    end do  
    num_clus = 1
    allocate(means(num_clus))
    call sortmeans(data1, 100, means, labels) 
    counts = 0
    do i = 1, ldim(1)
        do j = 1, ldim(2)
            if(rmat(i,j,1) > 0.1) then 
                rmat(i,j,1) = labels(counts)
                counts = counts + 1
            else
                rmat(i,j,1) = 0. 
            end if
        end do 
    end do  
    allocate(clus_stk(num_clus))
    allocate(gau2D_stk(num_clus))
    allocate(rmat_new(ldim(1), ldim(2), 1), source = 0.)
    call gau2D_sum%new(ldim, smpd)
    do l = 1, num_clus
        rmat_new = rmat
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                if( nint(rmat_new(i,j,1)) /= l) rmat_new(i,j,1) = 0.
            end do 
        end do
        call clus_stk(l)%new(ldim, smpd)
        call clus_stk(l)%set_rmat(rmat_new, .false.)
        call gauss2Dfit(clus_stk(l), cen_gauss, cov_gauss, corr, gau2D_stk(l))
        gau2D_stk(l) = gau2D_stk(l)*l
        call gau2D_sum%add( gau2D_stk(l))
    end do 
    ! motion correction to get 
    
    call find_ldim_nptcls(stk_directory,ldim,nptcls)
    allocate(defm_stack(7))
    ldim = [ldim(1),ldim(2),1]
    do i = 871,877
        call defm_stack(i)%new(ldim,smpd)
        call defm_stack(i)%read(stk_directory, i=i)
        call defm_stack(i)%vis()
    end do 
    ! call defm_model%new([50,50])
    ! params_glob => params
    ! params_glob%lpstart = 1.
    ! call defm_model%correct(20., defm_stack, fname)
end program ptcl_seg
