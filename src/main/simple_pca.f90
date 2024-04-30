! abstract strategy3D base class
module simple_pca
implicit none

public :: pca
private

type, abstract :: pca
    integer :: N          !< nr of data vectors
    integer :: D          !< nr of components in each data vec
    integer :: Q          !< nr of components in each latent vec
  contains
    procedure(generic_new),      deferred :: new
    procedure(generic_get_feat), deferred :: get_feat
    procedure(generic_generate), deferred :: generate
    procedure(generic_master),   deferred :: master
    procedure(generic_kill),     deferred :: kill
end type pca

abstract interface

    subroutine generic_new( self, N, D, Q )
        import :: pca
        class(pca), intent(inout) :: self
        integer,    intent(in)    :: N, D, Q
    end subroutine generic_new

    function generic_get_feat( self, i ) result( feat )
        import :: pca
        class(pca), intent(inout) :: self
        integer,    intent(in)    :: i
        real,       allocatable   :: feat(:)
    end function generic_get_feat

    subroutine generic_generate( self, i, avg, dat )
        import :: pca
        class(pca), intent(inout) :: self
        integer,    intent(in)    :: i
        real,       intent(in)    :: avg(self%D)
        real,       intent(inout) :: dat(self%D)
    end subroutine generic_generate

    subroutine generic_master( self, pcavecs, maxpcaits )
        import :: pca
        class(pca),        intent(inout) :: self
        real,              intent(in)    :: pcavecs(self%D,self%N)
        integer, optional, intent(in)    :: maxpcaits
    end subroutine generic_master

    subroutine generic_kill( self )
        import :: pca
        class(pca), intent(inout) :: self
    end subroutine generic_kill

end interface

end module simple_pca
