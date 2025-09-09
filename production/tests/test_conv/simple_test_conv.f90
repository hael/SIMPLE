program simple_test_conv
include 'simple_lib.f08'
implicit none
integer, parameter :: pftsz=10, NKER=5, HALF=NKER/2, MID=HALF+1
integer :: nums1(pftsz), nums2(pftsz), frc(pftsz), i
do i = 1, pftsz
    nums1(i) = floor(ran3() * pftsz) + 1
    nums2(i) = floor(ran3() * pftsz) + 1
enddo
do i = MID, pftsz-HALF    
    frc(i) = sum(nums1(i-HALF:i+HALF) * nums2(i-HALF:i+HALF))
enddo
do i = 1, MID-1
    frc(i) = sum(nums1(1:i+HALF) * nums2(1:i+HALF)) + sum(nums1(pftsz-HALF+i:pftsz) * nums2(pftsz-HALF+i:pftsz))
enddo
do i = pftsz-HALF+1, pftsz
    frc(i) = sum(nums1(1:HALF-(pftsz-i)) * nums2(1:HALF-(pftsz-i))) + sum(nums1(i-HALF:pftsz) * nums2(i-HALF:pftsz))
enddo
print *, nums1
print *, nums2
print *, nums1 * nums2
print *, frc
end program simple_test_conv
