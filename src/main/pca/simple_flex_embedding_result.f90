!@descr: model-neutral particle embedding returned by flex workflows
module simple_flex_embedding_result
use simple_core_module_api
implicit none
private

public :: flex_embedding_result

type :: flex_embedding_result
    character(len=32) :: method = 'none'
    integer :: nptcls = 0
    integer :: ncomp  = 0
    integer, allocatable :: pinds(:)
    real(dp), allocatable :: z(:,:)
    real(dp), allocatable :: eigvals(:)
contains
    procedure :: kill => kill_flex_embedding_result
end type flex_embedding_result

contains

    subroutine kill_flex_embedding_result( self )
        class(flex_embedding_result), intent(inout) :: self
        if( allocated(self%pinds) ) deallocate(self%pinds)
        if( allocated(self%z) ) deallocate(self%z)
        if( allocated(self%eigvals) ) deallocate(self%eigvals)
        self%method='none'; self%nptcls=0; self%ncomp=0
    end subroutine kill_flex_embedding_result

end module simple_flex_embedding_result
