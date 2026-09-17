!@descr: validates the flex_pca embedding cache and the kernel/state-weight contracts
program simple_test_flex_pca
use simple_core_module_api
use simple_flex_pca_model,   only: test_flex_pca_embedding_cache_io, test_flex_pca_auto_settings, &
    &test_flex_pca_population_floor
use simple_flex_pca_weights, only: test_flex_pca_kernel_bandwidth, test_flex_pca_state_weights
use simple_flex_pca_deconv,  only: test_flex_pca_deconv
implicit none

! resume path
call test_flex_pca_embedding_cache_io
! state stage
call test_flex_pca_kernel_bandwidth
call test_flex_pca_state_weights
call test_flex_pca_population_floor
! derived settings
call test_flex_pca_auto_settings
! latent deconvolution (calibration, held-out K, posterior means)
call test_flex_pca_deconv

call simple_end('**** SIMPLE_TEST_FLEX_PCA NORMAL STOP ****')

end program simple_test_flex_pca
