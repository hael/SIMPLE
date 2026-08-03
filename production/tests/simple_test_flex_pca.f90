!@descr: validates the flex_pca embedding cache, kernel/state-weight contracts and packed covariance solve
program simple_test_flex_pca
use simple_core_module_api
use simple_flex_pca_columns, only: test_flex_pca_svec_isometry, test_flex_pca_packed_solve
use simple_flex_pca_model,   only: test_flex_pca_embedding_cache_io, test_flex_pca_kernel_bandwidth, &
    &test_flex_pca_state_weights
implicit none

! numerical core of the reduced covariance solve
call test_flex_pca_svec_isometry
call test_flex_pca_packed_solve
! resume path
call test_flex_pca_embedding_cache_io
! state stage
call test_flex_pca_kernel_bandwidth
call test_flex_pca_state_weights

call simple_end('**** SIMPLE_TEST_FLEX_PCA NORMAL STOP ****')

end program simple_test_flex_pca
