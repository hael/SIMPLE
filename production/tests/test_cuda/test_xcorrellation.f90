module test_xcorrelation
    include 'simple_lib.f08'
    use simple_cuda_kernels
    use simple_image, only: image
    implicit none

    integer(dp),parameter :: Nmax=10000000
    complex, allocatable :: vec1(:),vec2(:)
    real(dp) :: vcpu, vcuda
    type(image) :: img1, img2

    integer(timer_int_kind) :: t1
    call random_seed(0)
    allocate(vec1(Nmax),vec2(Nmax))
    call random_number(vec1)
    call random_number(vec2)
    t1= tic()
    vcor = cpu_correlation(vec1, vec2, Nmax)
    print *, ' Cross correlation : ', vcor, ' in time: ', toc(t1)

    t1= tic()
    vcor = cuda_correlation(vec1, vec2, Nmax)
    print *, ' Cross correlation : ', vcor, ' in time: ', toc(t1)


contains

    function cuda_correlation(buf1, buf2, N) result(xcorr)
        complex, allocatable :: vec1(:), vec2(:)
        complex
        call kernelrcorr(cmat1,cmat2, r, ldim, sqlp, sqhp, 64,8)
           print *, " corr_cuda r correlation ",r, toc(t2)
            call kernelsumcsq(cmat1, sumasq, ldim, 64, 8 )
            print *, " corr_cuda A sum squared ", sumasq, toc()
            call kernelsumcsq(cmat2, sumbsq, ldim, 64, 8 )
            print *, " corr_cuda A sum squared ",sumbsq,  toc()


            if( sumasq < TINY .or. sumbsq < TINY )then
                r = 0.
            else
                r = r / sqrt(sumasq * sumbsq)
            endif

    end function cuda_correlation

!!  Calculate correlation coefficient between two numeric vectors
    function cpu_correlation(buf1, buf2, N) result(xcorr)
        integer(dp), intent(in) :: N
        complex, intent(in) :: buf1(N),buf2(N)
        complex(dp) :: X, Y, M1 , M2 , D1 , D2, K
        real(dp) :: xcorr
        integer :: i
        M1 = 0.0
        M2 = 0.0
        D1 = 0.0
        D2 = 0.0
        K  = 0.0
        do i = 1, N
            X = real(buf1(i),dp)
            Y = real(buf2(i),dp)
            M1 = M1 + X
            M2 = M2 + Y
            D1 = D1 + X * X
            D2 = D2 + Y * Y
            K  = K + X * Y
        enddo

        M1 = M1/real(N,dp)
        M2 = M2/real(N,dp)
        D1 = D1 - M1 * M1 * N
        D2 = D2 - M2 * M2 * N
        xcorr  = (K - M1 * M2 * N) / sqrt(D1 * D2)

    end function cpu_correlation


end module test_xcorrelation
