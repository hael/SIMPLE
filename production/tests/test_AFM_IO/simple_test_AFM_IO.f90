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
type(AFM_image), allocatable :: stack_stack(:)
type(AFM_image), allocatable :: stack_avg(:)
real, allocatable            :: mask_array(:,:,:,:)
logical, allocatable         :: log_mat(:,:), log_mat_corr(:,:)
real, allocatable            :: corr_mat(:,:), corr_mat_ex(:,:,:)
integer, allocatable         ::  max_box(:)
type(pickseg), allocatable   :: pick_array(:)
type(image)                  :: img_edg  
type(image)                  :: test_img, test_denoised     
type(image), allocatable     :: pick_mat(:,:)
type(image), allocatable     :: all_pick_vec(:)
type(image)                  :: test_vec(3), mrc_test_mat 
type(image), allocatable     :: pick_vec(:)
type(binimage), allocatable  :: bin_cc_array(:)
type(parameters), target     :: params
type(aff_prop)               :: aff_prop_clus   
type(spec_clust)             :: spec_clust_test    
integer, allocatable         :: clus_centers(:)
integer, allocatable         :: clus_labels(:)
real                         :: clus_simsum
real                         :: start, finish, count, rmat(150,150,1), hp = 60., lp = 10., comp_val
real, allocatable            :: corr_mat_test(:,:)
integer                      :: clip_len, file_num, file_iter, test_dim(3) = [100, 100, 1], i, j, rows, cols, ncls_ini, pick_ldim(3), sums, max, max_i, dims(2), pair(2)
character(len=LONGSTRLEN), allocatable  :: file_list(:)
real    :: thresh1(2) = [0.5, 0.2]
character(len = 255)    :: directory = '/Users/atifao/Downloads/IBW/'
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
call cline%set('trs',     25.0)
call cline%set('box',     150)
call cline%set('smpd',    4.89)
call params%new(cline)
call cpu_time(start)
allocate(stack_stack(size(file_list)))
allocate(stack_avg(size(file_list)))
allocate(mask_array(1024, 1024, 1, size(file_list)))
allocate(pick_array(size(file_list)))
allocate(bin_cc_array(size(file_list)))
! can adjust max number of particles sampled for clustering.
allocate(pick_mat(size(file_list), 500))
do file_iter = 1, size(file_list)
    call read_ibw(stack_stack(file_iter), file_list(file_iter))
    call align_avg(stack_stack(file_iter), stack_avg(file_iter))
    if(file_list(file_iter) == trim(directory) // 'Cob_450007.ibw' .or. file_list(file_iter) == trim(directory) // 'Cob_450010.ibw') then 
        do clip_len = 1, size(stack_avg(file_iter)%img_array)
            call stack_avg(file_iter)%img_array(clip_len)%clip_inplace([1024, 900, 1])
        end do
    end if
    stack_avg(file_iter)%stack_string = get_fbody(basename(trim(file_list(file_iter))), 'ibw')
    call pick_valid(stack_avg(file_iter), stack_avg(file_iter)%stack_string, pick_array(file_iter))
    call zero_padding(stack_avg(file_iter))
    call hough_lines(stack_avg(file_iter)%img_array(9), [PI/2 - PI/180, PI/2 + PI/180], mask_array(:, :, :, file_iter))
    call mask42D(stack_avg(file_iter)%img_array(9), pick_array(file_iter), bin_cc_array(file_iter), mask_array(:, :, :, file_iter), pick_mat(file_iter, :))
    ! if(file_iter > 1) exit 
end do
params_glob%smpd = pick_mat(1,1)%get_smpd()
allocate(log_mat(size(file_list), 500), source = .false.)
allocate(max_box(size(file_list)))
do rows = 1, size(file_list)
    max_box(rows) = maxval(pick_mat(rows,2)%get_ldim())
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
    call pick_vec(i)%write('/Users/atifao/Downloads/python_test/extracted.mrc', i)
end do 
call calc_inplane_invariant_corrmat(pick_vec, hp, lp, corr_mat_test, .true.)
dims = shape(corr_mat_test)
allocate(log_mat_corr(dims(1), dims(2)), source = .true.)

! 10 most similar pairs
comp_val = 0.9 
count = 0
do while (count < 11)
    pair = maxloc(corr_mat_test, mask = corr_mat_test .lt. comp_val)
    comp_val = maxval(corr_mat_test, mask = corr_mat_test .lt. comp_val)
    print *, comp_val 
    call pick_vec(pair(1))%vis()
    call pick_vec(pair(2))%vis()
    count = count + 1
end do 

call mrc_test_mat%new([dims(1),dims(2),1], 1.)
allocate(corr_mat_ex(dims(1), dims(2), 1))
corr_mat_ex(:,:,1) = corr_mat_test
call mrc_test_mat%set_rmat(corr_mat_ex, .false.)
call mrc_test_mat%write('/Users/atifao/Downloads/python_test/mat.mrc')
! call mrc_test_mat%vis()
! print *, sum(corr_mat_test)
! call pre_proc(pick_vec)
! call calc_inplane_invariant_corrmat(pick_vec, hp, lp, corr_mat_test, .true.)
! call aff_prop_clus%new(ncls_ini, corr_mat_test)
! call aff_prop_clus%propagate(clus_centers,clus_labels,clus_simsum)
! ! print *, clus_labels
! max = 0 
! do i = 1, maxval(clus_labels)
!     if (max < count(clus_labels == i)) then 
!         max = count(clus_labels == i)
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

! testing corr_mat. take particle roate and shift and get corrmat 

call cpu_time(finish)
print *, finish - start 
end program AFM_File_IO