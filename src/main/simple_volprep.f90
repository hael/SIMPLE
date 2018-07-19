module simple_volprep
include 'simple_lib.f08'
use simple_parameters, only: params_glob
implicit none

contains

    subroutine read_and_prep_vol( volfname, volprojobj )
        use simple_projector,   only: projector
        class(projector), intent(inout) :: volprojobj
        character(len=*), intent(in)    :: volfname
        real    :: smpd_target, smpd_sc, scale
        integer :: ldim(3), ifoo, box_sc
        ! find logical dimension
        call find_ldim_nptcls(volfname, ldim, ifoo)
        ! create conforming projector object
        call volprojobj%new(ldim,params_glob%smpd)
        ! read
        call volprojobj%read(volfname)
        ! mask
        call volprojobj%mask(params_glob%msk,'soft')
        ! FT
        call volprojobj%fft()
    end subroutine read_and_prep_vol

end module simple_volprep
