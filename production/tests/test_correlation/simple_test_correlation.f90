program simple_test_correlation
    include 'simple_lib.f08'
    use simple_nano_detect_atoms
    use simple_stat
    use simple_image

    ! want to compare atom image to essentially blank image
    ! will compare to completely blank image (all zeroes) and outlier position in image
    type(nano_picker) :: np 
    type(image) :: blank, boximg
    real :: smpd, boximg_corr
    real, allocatable :: rmat_blank(:,:,:)
    character(len=100) :: filename_exp, pdbfile_ref
    character(len=2) :: element
    integer :: offset, peak_thres_level, ldim_box(3), center(3), pos(3)

    filename_exp     = 'rec_merged.mrc'
    pdbfile_ref      = 'reference.pdb'
    element          = 'PT'
    smpd             = 0.358
    peak_thres_level = 2
    offset           = 2
    ldim_box         = [8,8,8]
    center           = [4,4,4]

    allocate(rmat_blank(8,8,8), source=0.)
    pos = [106, 106, 70]
    call np%new(smpd, element, filename_exp, peak_thres_level, offset, denoise=.true.)
    call np%simulate_atom()
    call np%extract_img_from_pos(pos, boximg)
    call boximg%write('boximg.mrc')
    call np%one_box_corr(pos, circle=.true., corr=boximg_corr)
    print *, 'boximg_corr = ', boximg_corr

end program simple_test_correlation
