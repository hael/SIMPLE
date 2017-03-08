program simple_test_combinatorics
use simple_combinatorics, only: balanced_diverse_labeling
implicit none
integer, allocatable :: configs_diverse(:,:)
configs_diverse = balanced_diverse_labeling(13, 3, 5)

print *, configs_diverse(1,:)
print *, '****************'
print *, configs_diverse(2,:)
print *, '****************'
print *, configs_diverse(3,:)
print *, '****************'
print *, configs_diverse(4,:)
print *, '****************'
print *, configs_diverse(5,:)


end program simple_test_combinatorics