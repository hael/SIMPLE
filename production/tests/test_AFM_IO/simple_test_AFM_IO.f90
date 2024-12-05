program AFM_File_IO
include 'simple_lib.f08'
use simple_AFM_image
use simple_syslib
use simple_binimage
use simple_segmentation
use simple_pickseg

type(AFM_image), allocatable :: stack_stack(:)
type(AFM_image), allocatable :: stack_avg(:)
real, allocatable            :: mask_array(:,:,:,:)
type(pickseg), allocatable   :: pick_array(:)
type(image)                  :: img_edg  
type(image)                  :: test_img, test_denoised     
type(image), allocatable     :: pick_mat(:,:)
type(binimage), allocatable  :: bin_cc_array(:)
type(parameters), target     :: params
real                         :: start, finish, count
integer                      :: clip_len, file_num, file_iter, test_dim(3) = [1000, 1000, 1], i, j
character(len=LONGSTRLEN), allocatable  :: file_list(:)
real    :: thresh1(2) = [0.5, 0.2]
character(len = 255)    :: directory = '/Users/atifao/Downloads/IBW/'

params_glob => params
params_glob%pcontrast = 'white'
params_glob%lp        = 10.
params_glob%nsig      = 1.5 
call cpu_time(start)
call simple_list_files(trim(directory) // '*.ibw', file_list)
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
    if(file_iter > 3) exit 
end do

call corr_clus( pick_mat )
call cpu_time(finish)
print *, finish - start 
end program AFM_File_IO