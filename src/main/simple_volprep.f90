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
        ! clip
        if( params_glob%boxmatch < params_glob%box ) &
            call volprojobj%clip_inplace([params_glob%boxmatch,params_glob%boxmatch,params_glob%boxmatch])
        ! mask
        call volprojobj%mask(params_glob%msk,'soft')
        ! FT
        call volprojobj%fft()
        ! scale
        smpd_target = params_glob%lp*LP2SMPDFAC
        call autoscale(params_glob%boxmatch, params_glob%smpd, smpd_target, box_sc, smpd_sc, scale)
        call volprojobj%clip_inplace([box_sc,box_sc,box_sc])
        call volprojobj%set_smpd(smpd_sc)
        ! return as FT
    end subroutine read_and_prep_vol

end module simple_volprep
