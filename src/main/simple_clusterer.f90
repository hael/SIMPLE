! high-level clustering module
module simple_clusterer
use simple_defs
implicit none

contains

    ! routine for SHC-based clustering of orientation based on Geodesic distance
    subroutine cluster_shc_oris( os, ncls )
        use simple_cluster_shc, only: cluster_shc
        use simple_oris,        only: oris
        class(oris), intent(inout) :: os
        integer,     intent(in)    :: ncls
        type(oris)        :: labels
        type(cluster_shc) :: shcc
        real, allocatable :: smat(:,:)
        integer :: noris, i
        real    :: sim
        noris = os%get_noris()
        ! calculate the similarities
        smat = os%gen_smat()
        ! cluster
        call labels%new(noris)
        call shcc%new(noris, ncls, smat, labels, minsim=-2.0*sqrt(2.0))
        call shcc%shc(.true., 'class', sim)
        ! transfer the class labels to os
        do i=1,noris
            call os%set(i, 'class', labels%get(i,'class'))
        end do
        deallocate(smat)
    end subroutine cluster_shc_oris

end module simple_clusterer
