program AFM_File_IO
include 'simple_lib.f08'
use simple_AFM_image
use simple_syslib
use simple_binimage
use simple_segmentation

type(AFM_image), allocatable :: stack_stack(:)
type(AFM_image), allocatable :: stack_avg(:)
type(image)                  :: img_edg  
type(image)                  :: test_img, test_denoised     
type(binimage)               :: test_bin 
type(parameters), target     :: params
real                         :: start, finish 
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

do file_iter = 4, size(file_list)
    call read_ibw(stack_stack(file_iter), file_list(file_iter))
    call stack_stack(file_iter)%img_array(1)%vis()
    ! call zero_padding(stack_stack(file_iter))
    call hough_lines(stack_stack(file_iter)%img_array(1), img_edg, [PI/2 - PI/180, PI/2 + PI/180])
!     ! call img_edg%vis()
!     ! call stack_stack(file_iter)%img_array(1)%write_jpg('test.jpg', norm = .true. )
!     ! call otsu_img(stack_stack(file_iter)%img_array(1))
!     ! call test_bin%transfer2bimg(stack_stack(file_iter)%img_array(1))
!     ! call test_bin%vis()
!     ! call test_bin%erode()
!     ! call test_bin%dilate()
!     ! call test_bin%vis()
!     ! call stack_stack(file_iter)%img_array(1)%write('Users/atifao/Downloads/' // int2str(file_iter) // '.mrc')
!     ! call stack_stack(file_iter)%img_array(1)%write_jpg('test.jpg')
!     ! call stack_stack(file_iter)%img_array(1)%ICM2D(0.5)
!     ! call stack_stack(file_iter)%img_array(2)%ICM2D(0.5)
!     ! call stack_stack(file_iter)%img_array(1)%vis()
!     ! call stack_stack(file_iter)%img_array(1)%vis()
    
!     ! call align_avg(stack_stack(file_iter), stack_avg(file_iter))
!     ! if(file_list(file_iter) == trim(directory) // 'Cob_450007.ibw' .or. file_list(file_iter) == trim(directory) // 'Cob_450010.ibw') then 
!     !     do clip_len = 1, size(stack_avg(file_iter)%img_array)
!     !         call stack_avg(file_iter)%img_array(clip_len)%clip_inplace([1024, 900, 1])
!     !     end do
!     ! end if
!     ! stack_avg(file_iter)%stack_string = get_fbody(basename(trim(file_list(file_iter))), 'ibw')
!     ! call pick_valid(stack_avg(file_iter), stack_avg(file_iter)%stack_string)
    if(file_iter > 3) exit 
end do

call cpu_time(finish)
print *, finish - start 

! test image with only one line 
! call test_img%new(test_dim, 1.0)
! do i = 1, test_dim(1)
!     do j = 1, test_dim(2)

!         if(j == 300 .and. i > 10 .and. i < 895 ) then 
!             call test_img%set_rmat_at(i, j, 1, 1.0)
!         end if 

!         if( j == 400 .or. j == 600 .or. j == 800) then 
!             call test_img%set_rmat_at(i, j, 1, 1.0)
!         end if 

!         if(j == 500 .and. i > 200 .and. i < 451 ) then 
!             call test_img%set_rmat_at(i, j, 1, 1.0)
!         end if 

!         if(j == 500 .and. i > 600 .and. i < 751 ) then 
!             call test_img%set_rmat_at(i, j, 1, 1.0)
!         end if 

!         if(j == 500 .and. i > 800 .and. i < 895 ) then 
!             call test_img%set_rmat_at(i, j, 1, 1.0)
!         end if 
        

!         if(j == 200 .and. i > 600 .and. i < 625 ) then 
!             call test_img%set_rmat_at(i, j, 1, 1.0)
!         end if 

!         if(i == 500 .and. j > 200 .and. j < 451 ) then 
!             call test_img%set_rmat_at(i, j, 1, 1.0)
!         end if 

!     end do 
! end do 
! call test_img%vis()
! call test_img%write_jpg('test_hough.jpg')
! call hough_lines(test_img, test_denoised)

end program AFM_File_IO