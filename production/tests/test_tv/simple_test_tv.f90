program simple_test_tv
    include 'simple_lib.f08'
    use simple_tvfilter, only: tvfilter
    use simple_image,    only: image
    implicit none
    
    character(len=*), parameter :: test_in_img  = 'test.mrc'
    character(len=*), parameter :: test_out_img = 'test_out.mrc'
    integer,          parameter :: box    = 31
    real,             parameter :: smpd   = 1.275
    real,             parameter :: lambda = 0.01
    type(image)     :: img
    type(tvfilter)  :: tvfilter_obj
    
    call img%new([box,box,1], smpd)
    call img%read(test_in_img, 1)
    call tvfilter_obj%apply_filter(img, lambda)
    call img%write(test_out_img, 1)
    
end program simple_test_tv
    