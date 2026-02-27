module simple_cluster2D_common
use simple_core_module_api
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_cmdline,            only: cmdline
use simple_commanders_mkcavgs, only: commander_make_cavgs, commander_make_cavgs_distr
use simple_commanders_imgops,  only: commander_scale
implicit none

public :: init_cluster2D_refs, handle_objfun
private
#include "simple_local_flags.inc"

contains

    !> Initialize references - common to both execution modes
    subroutine init_cluster2D_refs(cline, params, build)
        use simple_procimgstk, only: copy_imgfile
        class(cmdline),   intent(inout) :: cline
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(commander_scale) :: xscale
        type(cmdline)         :: cline_make_cavgs, cline_scalerefs
        type(string)          :: refs_sc
        logical               :: l_scale_inirefs
        integer               :: cnt, iptcl, ptclind
        if( cline%defined('refs') ) return  ! References already provided
        cline_make_cavgs = cline
        params%refs      = 'start2Drefs'     //MRC_EXT
        params%refs_even = 'start2Drefs_even'//MRC_EXT
        params%refs_odd  = 'start2Drefs_odd' //MRC_EXT
        l_scale_inirefs  = .false.
        if( build%spproj%is_virgin_field('ptcl2D') .or. params%which_iter <= 1 )then
            if( params%tseries .eq. 'yes' )then
                call init_tseries_refs(cline, params, build, cline_make_cavgs, l_scale_inirefs)
            else
                call init_standard_refs(cline, params, build, cline_make_cavgs, l_scale_inirefs)
            endif
        else
            ! Use existing classifications
            call cline_make_cavgs%set('refs', params%refs)
            call execute_make_cavgs(cline_make_cavgs, cline, params)
            l_scale_inirefs = .false.
        endif
        ! Scale references if needed
        if( l_scale_inirefs )then
            refs_sc = 'refs'//SCALE_SUFFIX//MRC_EXT
            call cline_scalerefs%set('stk',    params%refs)
            call cline_scalerefs%set('outstk', refs_sc)
            call cline_scalerefs%set('smpd',   params%smpd)
            call cline_scalerefs%set('newbox', params%box_crop)
            call xscale%execute(cline_scalerefs)
            call simple_rename(refs_sc, params%refs)
        endif
        ! Create even/odd copies
        call copy_imgfile(params%refs, params%refs_even, params%smpd_crop, [1,params%ncls])
        call copy_imgfile(params%refs, params%refs_odd,  params%smpd_crop, [1,params%ncls])
        call cline%set('refs', params%refs)
    end subroutine init_cluster2D_refs

    subroutine init_tseries_refs(cline, params, build, cline_make_cavgs, l_scale_inirefs)
        use simple_procimgstk, only: selection_from_tseries_imgfile
        class(cmdline),   intent(inout) :: cline
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(cmdline),    intent(inout) :: cline_make_cavgs
        logical,          intent(out)   :: l_scale_inirefs
        integer :: cnt, iptcl, ptclind
        if( cline%defined('nptcls_per_cls') )then
            if( build%spproj%os_ptcl2D%any_state_zero() )then
                THROW_HARD('cluster2D_nano does not allow state=0 particles, prune project before execution')
            endif
            cnt = 0
            do iptcl=1,params%nptcls,params%nptcls_per_cls
                cnt = cnt + 1
                params%ncls = cnt
                do ptclind=iptcl,min(params%nptcls, iptcl + params%nptcls_per_cls - 1)
                    call build%spproj%os_ptcl2D%set(ptclind, 'class', cnt)
                end do
            end do
            call cline%set('ncls', params%ncls)
            call cline_make_cavgs%set('ncls', params%ncls)
            call cline_make_cavgs%set('refs', params%refs)
            call execute_make_cavgs(cline_make_cavgs, cline, params)
            l_scale_inirefs = .false.
        else
            if( trim(params%refine).eq.'inpl' )then
                params%ncls = build%spproj%os_ptcl2D%get_n('class')
                call cline%set('ncls', params%ncls)
                call cline_make_cavgs%set('ncls', params%ncls)
                call cline_make_cavgs%delete('tseries')
                call cline_make_cavgs%set('refs', params%refs)
                call execute_make_cavgs(cline_make_cavgs, cline, params)
                l_scale_inirefs = .false.
            else
                call selection_from_tseries_imgfile(build%spproj, params%refs, params%box, params%ncls)
                l_scale_inirefs = .true.
            endif
        endif
    end subroutine init_tseries_refs

    subroutine init_standard_refs(cline, params, build, cline_make_cavgs, l_scale_inirefs)
        use simple_procimgstk, only: random_selection_from_imgfile, noise_imgfile
        class(cmdline),   intent(inout) :: cline
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(cmdline),    intent(inout) :: cline_make_cavgs
        logical,          intent(out)   :: l_scale_inirefs
        integer :: iptcl
        select case(trim(params%cls_init))
            case('ptcl')
                ! Initialization from raw images
                call random_selection_from_imgfile(build%spproj, params%refs, params%box, params%ncls)
                l_scale_inirefs = .true.
            case('rand')
                ! From noise
                call noise_imgfile(params%refs, params%ncls, params%box_crop, params%smpd_crop)
                l_scale_inirefs = .false.
            case('randcls')
                if(.not.cline%defined('ncls')) THROW_HARD('NCLS must be provide with CLS_INIT=RANDCLS')
                do iptcl=1,params%nptcls
                    if( build%spproj_field%get_state(iptcl) == 0 ) cycle
                    call build%spproj_field%set(iptcl, 'class', irnd_uni(params%ncls))
                    call build%spproj_field%set(iptcl, 'w',     1.0)
                    call build%spproj_field%e3set(iptcl,ran3()*360.0)
                end do
                call build%spproj%write_segment_inside(params%oritype, params%projfile)
                call cline_make_cavgs%set('refs', params%refs)
                call execute_make_cavgs(cline_make_cavgs, cline, params)
                l_scale_inirefs = .false.
            case DEFAULT
                THROW_HARD('Unsupported mode of initial class generation CLS_INIT='//trim(params%cls_init))
        end select
    end subroutine init_standard_refs

    subroutine execute_make_cavgs(cline_make_cavgs, cline, params)
        type(cmdline),    intent(inout) :: cline_make_cavgs
        class(cmdline),   intent(in)    :: cline
        type(parameters), intent(in)    :: params
        type(commander_make_cavgs)       :: xmake_cavgs
        type(commander_make_cavgs_distr) :: xmake_cavgs_distr
        if( (params%nparts > 1) .and. (.not.cline%defined('part')) )then
            call xmake_cavgs_distr%execute(cline_make_cavgs)
        else
            call xmake_cavgs%execute(cline_make_cavgs)
        endif
    end subroutine execute_make_cavgs

    subroutine handle_objfun(params, cline)
        type(parameters), intent(inout) :: params
        class(cmdline),   intent(inout) :: cline
        select case( params%cc_objfun )
            case( OBJFUN_CC )
                ! Making sure euclid options are turned off
                params%l_needs_sigma = .false.
                params%needs_sigma   = 'no'
                call cline%set('ml_reg','no')
                params%ml_reg        = 'no'
                params%l_ml_reg      = .false.
            case( OBJFUN_EUCLID )
                ! Euclidean distance objective - sigma needed
                ! All set by parameter initialization
            case DEFAULT
                THROW_HARD('Unknown objective function')
        end select
    end subroutine handle_objfun

end module simple_cluster2D_common
