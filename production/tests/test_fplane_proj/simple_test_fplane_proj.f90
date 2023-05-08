program simple_test_fplane_proj
include 'simple_lib.f08'
implicit none
integer, parameter :: K_MIN = 5, K_MAX = 15
real,    parameter :: THETA = 15.
real    :: sig1(K_MAX), sig2(K_MAX), k_proj
integer :: k, k_low, k_high
sig1 = 0.
sig2 = 0.
sig1(K_MIN:K_MAX) = 1.
sig2(K_MIN:K_MAX) = 1.
print *, sig2
do k = K_MIN, K_MAX
    k_proj = real(k) * cos(THETA * pi / 180.)
    k_low  = floor(  k_proj)
    k_high = ceiling(k_proj)
    if( k_low  >= K_MIN .and. k_low  <= K_MAX ) sig2(k_low)  = sig2(k_low)  + sig1(k)
    if( k_high >= K_MIN .and. k_high <= K_MAX ) sig2(k_high) = sig2(k_high) + sig1(k)
enddo
print *, sig2
end program simple_test_fplane_proj
