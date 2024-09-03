program simple_test_kpca
include 'simple_lib.f08'
use simple_kpca_svd, only: kpca_svd
implicit none
integer, parameter :: NP = 3, NS = 4, NC = 2, MAXPCAITS = 15
type(kpca_svd)     :: kpca_obj
call kpca_obj%new(NS, NP, NC)
call kpca_obj%test_all
end program simple_test_kpca
