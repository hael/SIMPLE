module simple_volprep
use simple_defs ! use all in there
implicit none

contains

    subroutine read_and_prep_vol( p, volfname, volprojobj )
        use simple_params,      only: params
        use simple_projector,   only: projector
        use simple_magic_boxes, only: autoscale
        use simple_imgfile,     only: find_ldim_nptcls
        class(params),    intent(in)    :: p
        class(projector), intent(inout) :: volprojobj
        character(len=*), intent(in)    :: volfname
        real    :: smpd_target, smpd_sc, scale
        integer :: ldim(3), ifoo, box_sc
        ! find logical dimension
        call find_ldim_nptcls(volfname, ldim, ifoo)
        ! create conforming projector object
        call volprojobj%new(ldim,p%smpd)
        ! read
        call volprojobj%read(volfname)
        ! clip
        if( p%boxmatch < p%box ) call volprojobj%clip_inplace([p%boxmatch,p%boxmatch,p%boxmatch])
        ! mask
        call volprojobj%mask(p%msk,'soft')
        ! FT
        call volprojobj%fwd_ft
        ! scale
        smpd_target = p%lp*LP2SMPDFAC
        call autoscale(p%boxmatch, p%smpd, smpd_target, box_sc, smpd_sc, scale)
        call volprojobj%clip_inplace([box_sc,box_sc,box_sc])
        call volprojobj%set_smpd(smpd_sc)
        ! return as FT
    end subroutine read_and_prep_vol

end module simple_volprep
