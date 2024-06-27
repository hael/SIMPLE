module AFM_utils
include 'simple_lib.f08'
use simple_image,              only: image
use simple_pca_svd
use simple_linalg
use simple_segmentation

implicit none

contains

    subroutine LRR_2D(self, weight)
        type(image) :: self
        type(pca_svd) :: SVD
        real, allocatable        :: weight, A(:,:,:)
        A = self%get_rmat()
    end subroutine LRR_2D

end module AFM_utils

program simple_test_AFM
    include 'simple_lib.f08'
    use simple_image,              only: image
    use simple_segmentation
    use simple_tvfilter
    use simple_pickseg
    use simple_histogram
    type(image) :: AFM
    type(histogram) :: AFM_h
    real        :: smpd = 5.0
    integer     :: ldim(3), ifoo 
    real, allocatable :: rmat(:,:)
    character(*), parameter     :: fn_in = '/Users/atifao/Downloads/MRC_T/8.mrc'
    character(256)              :: fn_out, hist_out
    fn_out ='/Users/atifao/Downloads/MRC_Inv/8_out.mrc'
    hist_out ='/Users/atifao/Downloads/ABCDE.pdf'
    call find_ldim_nptcls(fn_in, ldim, ifoo, smpd = smpd) 
    call AFM%new(ldim, smpd)
    call AFM%read(fn_in)
    !call AFM%norm_minmax()
    !AFM = AFM*256.
    ! call AFM%ICM2D(1.0)
    call AFM%bin_inv()
    !call otsu_img(AFM, positive=.FALSE.)
    ! call AFM%write(fn_out)
    ! print *, AFM%get_cmat()
    call AFM_h%new(AFM, 100)
    call AFM_h%plot(hist_out)
end program simple_test_AFM
