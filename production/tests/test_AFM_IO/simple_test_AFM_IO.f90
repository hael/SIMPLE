program AFM_File_IO
    include 'simple_lib.f08'
    use simple_AFM_image
    use simple_syslib
    use simple_binimage
    use simple_segmentation
    use simple_pickseg
    use simple_corrmat
    use simple_aff_prop
    use simple_spectral_clustering
    use simple_polarft_corrcalc
    type(AFM_image), allocatable :: stack_stack(:)
    type(AFM_image), allocatable :: stack_avg(:)
    real, allocatable            :: mask_array(:,:,:,:), irot_mat(:,:)
    logical, allocatable         :: log_mat(:,:), log_mat_corr(:,:)
    real, allocatable            :: corr_mat(:,:), corr_mat_ex(:,:,:)
    integer, allocatable         ::  max_box(:)
    type(pickseg), allocatable   :: pick_array(:)
    type(image)                  :: img_edg, rot, rot_sh, rot_ref, rot_ref_sh
    type(image)                  :: test_img, test_denoised, test_vec(2), test_1, test_2
    type(image), allocatable     :: pick_mat(:,:)
    type(image), allocatable     :: all_pick_vec(:)
    type(image)                  :: test_vec_corr(2), mrc_test_mat, rot_inv
    type(image), allocatable     :: pick_vec(:)
    type(binimage), allocatable  :: bin_cc_array(:)
    type(parameters), target     :: params
    type(aff_prop)               :: aff_prop_clus   
    type(spec_clust)             :: spec_clust_test    
    integer, allocatable         :: clus_centers(:)
    integer, allocatable         :: clus_labels(:)
    real                         :: clus_simsum, test_shift(2), new_smpd = 15. 
    real                         :: start, finish, count, rmat(150,150,1), hp = 60., lp = 10., comp_val
    real, allocatable            :: corr_mat_test(:,:), rot_test(:,:)
    integer                      :: clip_len, file_num, file_iter, test_dim(3) = [100, 100, 1], i, j, rows, cols, ncls_ini, pick_ldim(3), sums, max, max_i, dims(2), pair(2)
    character(len=LONGSTRLEN), allocatable  :: file_list(:)
    real    :: thresh1(2) = [0.5, 0.2]
    character(len = 255)    :: directory = '/Users/atifao/Downloads/IBW_1/'
    type(cmdline)          :: cline
    call simple_list_files(trim(directory) // '*.ibw', file_list)
    params_glob => params
    params_glob%pcontrast = 'white'
    params_glob%lp  = 10.
    params_glob%nsig  = 1.5 
    call cline%set('objfun','cc')
    call cline%set('ctf',    'no')
    call cline%set('objfun', 'cc')
    call cline%set('mkdir', 'no')
    call cline%set('lambda', 0.)
    call cline%set('trs',     50.0)
    call cline%set('box',     300)
    call cline%set('smpd',    4.89)
    params_glob%cc_objfun = 0
    params_glob%maxits_sh = 200
    ! params_glob%shbarrier = 'yes'
    call params%new(cline)
    call cpu_time(start)
    allocate(stack_stack(size(file_list)))
    
    allocate(stack_avg(size(file_list)))
    allocate(mask_array(1024, 1024, 1, size(file_list)))
    allocate(pick_array(size(file_list)))
    allocate(bin_cc_array(size(file_list)))
    allocate(pick_mat(size(file_list), 500))
    do file_iter = 1, size(file_list)
        call read_ibw(stack_stack(file_iter), file_list(file_iter))
        call align_avg(stack_stack(file_iter), stack_avg(file_iter))
        call stack_avg(file_iter)%img_array(1)%vis()
        if(file_list(file_iter) == trim(directory) // 'Cob_450007.ibw' .or. file_list(file_iter) == trim(directory) // 'Cob_450010.ibw') then 
            do clip_len = 1, size(stack_avg(file_iter)%img_array)
                call stack_avg(file_iter)%img_array(clip_len)%clip_inplace([1024, 900, 1])
            end do
        end if
        call stack_avg(file_iter)%img_array(1)%vis()
        stack_avg(file_iter)%stack_string = get_fbody(basename(trim(file_list(file_iter))), 'ibw')
        call pick_valid(stack_avg(file_iter), stack_avg(file_iter)%stack_string, pick_array(file_iter))
        call zero_padding(stack_avg(file_iter))
        call hough_lines(stack_avg(file_iter)%img_array(9), [PI/2 - PI/180, PI/2 + PI/180], mask_array(:, :, :, file_iter))
        call mask42D(stack_avg(file_iter)%img_array(9), pick_array(file_iter), bin_cc_array(file_iter), mask_array(:, :, :, file_iter), pick_mat(file_iter, :))
        if(file_iter > 0) exit 
    end do
    print *, 'picking worked'
    params_glob%smpd = pick_mat(1,1)%get_smpd()
    allocate(log_mat(size(file_list), 500), source = .false.) 
    allocate(max_box(size(file_list)))
    !need to rewrite this part 
    do rows = 1, size(file_list)
        max_box(rows) = maxval(pick_mat(rows,5)%get_ldim())
    end do 
    params_glob%box = maxval(max_box)
    pick_ldim = [maxval(max_box),maxval(max_box),1]
    ncls_ini = 0 
    do rows = 1, size(file_list)
        do cols = 1, 500
            if(sum(pick_mat(rows,cols)%get_ldim()) > 3) then
                call pick_mat(rows,cols)%pad_inplace(pick_ldim)
                log_mat(rows,cols) = .true.
                ncls_ini = ncls_ini + 1
            end if 
        end do 
    end do 
    pick_vec = pack(pick_mat, log_mat .eqv. .true.)
    do i = 1, size(pick_vec)
        call pick_vec(i)%write('mrc_for_clus.mrc', i)
    end do 
    ! allocate(rot_test(size(pick_vec),size(pick_vec)))
   
    
    allocate(rot_test(2,2))
    call rot_sh%copy(pick_vec(127))
    call rot_sh%rtsq(143.,42.,12.)
    call rot_sh%mirror('y')
    test_vec_corr(1) = pick_vec(127)
    test_vec_corr(2) = rot_sh
    
    ! call calc_inplane_invariant_corrmat(test_vec_corr, hp, lp, rot_test)
    ! print *, rot_test
    ! call calc_inplane_fast(pick_vec, hp, lp, rot_test)
    print *, rot_test

    ! print *, rot_test
    ! call pick_vec(107)%write('/Users/afan/ref.mrc')
    ! call rot_sh%copy(pick_vec(107))
    ! call rot_sh%rtsq(63.,42.,12.)
    ! call rot_sh%write('/Users/afan/orig.mrc')
    
    ! call calc_inplane_mag_corrmat(pick_vec, hp, lp, rot_test)
   ! 10 most similar pairs

    ! comp_val = 0.97
    ! count = 0
    ! do while (count < 20)
    !     pair = maxloc(rot_test, mask = rot_test .lt. comp_val)
    !     comp_val = maxval(rot_test, mask = rot_test .lt. comp_val)
    !     print *, comp_val 
    !     call pick_vec(pair(1))%vis()
    !     call pick_vec(pair(2))%vis()
    !     count = count + 1
    ! end do 

    ! call pick_vec(107)%write('/Users/afan/ref.mrc')
    ! call rot_sh%copy(pick_vec(107))
    ! call rot_sh%rtsq(63.,42.,12.)
    ! call rot_sh%write('/Users/afan/orig.mrc')
    
    ! pick_ldim = [params_glob%box,params_glob%box,1]
    ! call test_1%new(pick_ldim, params_glob%smpd)
    ! call test_2%new(pick_ldim, params_glob%smpd)
    ! call test_1%read('/Users/atifao/Downloads/ref.mrc')
    ! call test_2%read('/Users/atifao/Downloads/orig.mrc')
    ! call test_2%mirror('y')
    ! test_vec_corr(1) = test_1
    ! test_vec_corr(2) = test_2
    ! call test_1%write_jpg('ref.jpg')
    ! call test_2%write_jpg('orig.jpg')
    
    ! allocate(rot_test(2,2))
    ! call calc_inplane_mag_corrmat(test_vec_corr, hp, lp, rot_test)
    ! print *, 'rotations:', rot_test
    
    ! call test_2%rtsq(150.,0.,0.)
    
    ! call test_1%fft()
    ! call test_2%fft()
    ! call test_1%fcorr_shift(test_2, 300., test_shift, peak_interp = .true.)
    ! print *, test_shift 
    ! call rot_sh%rtsq(-maxval(irot_mat), 0., 0.)
    ! call rot_sh%write_jpg('derotated.jpg')
    
    ! test_vec_corr(1) = rot_sh
    
    ! call rot_sh%shift([test_shift(1), test_shift(2), 0.])
    
    ! call rot_sh%write_jpg('derot_deshift.jpg')
    ! print *, 'final sim:', pick_vec(107)%real_corr(rot_sh)
    
    
    
    ! call pick_vec(107)%fft()
    ! call rot_sh%fft()
    ! call pick_vec(107)%fcorr_shift(rot_sh, 200., test_shift)
    ! call pick_vec(107)%ifft()
    ! call rot_sh%ifft()
    ! print *, test_shift 
    
    
    
    
    
    
    ! test corrmat on 107
    ! call calc_inplane_invariant_corrmat(pick_vec, hp, lp, corr_mat_test, .true.)
    ! dims = shape(corr_mat_test)
    ! allocate(log_mat_corr(dims(1), dims(2)), source = .true.)
    
    
    ! call mrc_test_mat%new([dims(1),dims(2),1], 1.)
    ! allocate(corr_mat_ex(dims(1), dims(2), 1))
    ! corr_mat_ex(:,:,1) = corr_mat_test
    ! call mrc_test_mat%set_rmat(corr_mat_ex, .false.)
    ! call mrc_test_mat%write('/Users/atifao/Downloads/python_test/mat.mrc')
    ! call mrc_test_mat%vis()
    ! print *, sum(corr_mat_test)
    ! call pre_proc(pick_vec)
    ! call calc_inplane_invariant_corrmat(pick_vec, hp, lp, corr_mat_test, .true.)

    ! call aff_prop_clus%new(ncls_ini, rot_test)
    ! call aff_prop_clus%propagate(clus_centers,clus_labels,clus_simsum)
    ! print *, clus_labels
    ! do i = 1, size(clus_labels)
    !     if(clus_labels(i) == 4) call pick_vec(i)%vis()
    ! end do 

    ! max = 0 
    ! do i = 1, maxval(clus_labels)
    !     if(max < count(clus_labels == i)) then 
    !         ! max = count(clus_labels == i)
    !         max_i = i 
    !     end if 
    ! end do 
    ! do i = 1, size(clus_labels)
    !     if(clus_labels(i) == max_i) then 
    !         call pick_vec(i)%vis()
    !     end if 
    ! end do 
    ! print *, max, max_i
    
    ! do i = 1, 3
    !     call test_vec(i)%new([params_glob%box, params_glob%box, 1], params_glob%smpd)
    !     call random_number(rmat)
    !     call test_vec(i)%set_rmat(rmat,.false.)
    ! end do 
    ! call calc_inplane_invariant_corrmat(test_vec, hp, lp, corr_mat_test)
    ! print *, corr_mat_test
    ! should remove small particles... 
    call cpu_time(finish)
    print *, finish - start 
end program AFM_File_IO