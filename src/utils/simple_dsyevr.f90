module simple_dsyevr
implicit none
public :: s_dsyevr
private
interface s_dsyevr
    subroutine dsyevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
                    &abstol, m, w, z, ldz, isuppz, work, lwork,&
                    &iwork, liwork, info )
        character(len=1), intent(in)    :: jobz, range, uplo
        integer,          intent(in)    :: n, il, iu, ldz, lda, lwork, liwork
        integer,          intent(out)   :: m
        integer,          intent(out)   :: isuppz(*)
        real(kind(0d0)),  intent(in)    :: abstol, vl, vu
        integer,          intent(out)   :: iwork(*)
        integer,          intent(out)   :: info
        real(kind(0d0)),  intent(inout) :: a(lda,*)
        real(kind(0d0)),  intent(out)   :: work(*), z(ldz,*)
        real(kind(0d0)),  intent(out)   :: w(*)
    end subroutine dsyevr
end interface s_dsyevr
end module simple_dsyevr