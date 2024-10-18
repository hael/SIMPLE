program simple_test_clustering
include 'simple_lib.f08'
use simple_aff_prop
use simple_spectral_clustering
use simple_linalg
implicit none
call test_aff_prop
call test_spec_clust
end program simple_test_clustering
