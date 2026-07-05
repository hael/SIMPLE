!@descr: denoising and class split commanders
module simple_commanders_denoise
use simple_commanders_api
implicit none

#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_cls_split
  contains
    procedure :: execute => exec_cls_split
end type commander_cls_split

type, extends(commander_base) :: commander_denoise_project
  contains
        procedure :: execute => exec_denoise_project
end type commander_denoise_project

type, extends(commander_base) :: commander_map_params_from_den
  contains
        procedure :: execute => exec_map_params_from_den
end type commander_map_params_from_den

contains

    subroutine exec_cls_split( self, cline )
        use simple_cls_split_strategy
        use simple_commanders_mkcavgs, only: run_make_cavgs_workflow
        class(commander_cls_split), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        class(cls_split_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        type(cmdline)    :: cline_make
        strategy = create_cls_split_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params)
        if( allocated(strategy) ) deallocate(strategy)
        if( .not. cline%defined('part') )then
            if( trim(params%oritype) == 'ptcl2D' )then
                cline_make = cline
                call cline_make%set('prg',     'make_cavgs')
                call cline_make%set('projfile', params%projfile%to_char())
                if( .not. cline_make%defined('refs') ) call cline_make%set('refs', 'cls_split_cavgs.mrc')
                call cline_make%set('mkdir',   'no')
                call cline_make%set('oritype', 'ptcl2D')
                call cline_make%delete('class')
                call cline_make%delete('ncls')
                call cline_make%delete('gen_model')
                call run_make_cavgs_workflow(cline_make, from_distr_cmd=.true.)
                call cline_make%kill
            else
                THROW_WARN('cls_split class-average regeneration is only available for oritype=ptcl2D')
            endif
        endif
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_CLS_SPLIT NORMAL STOP ****')
    end subroutine exec_cls_split

    subroutine exec_denoise_project( self, cline )
        use simple_denoise_project_strategy
        class(commander_denoise_project), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        class(denoise_project_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        strategy = create_denoise_project_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params)
        if( allocated(strategy) ) deallocate(strategy)
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_DENOISE_PROJECT NORMAL STOP ****')
    end subroutine exec_denoise_project

    subroutine exec_map_params_from_den( self, cline )
        use simple_ori, only: ori
        class(commander_map_params_from_den), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: raw_proj, den_proj, outproj
        type(ori)        :: raw2d, den_ori, mapped_ori
        integer, allocatable :: den2raw(:)
        integer :: iden, iptcl, nptcls, nden2d, nden3d, nraw_active, nmap2d, nmap3d
        logical :: l_use_pind
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('projfile') ) call cline%set('projfile', 'map_params_from_den.simple')
        if( .not. cline%defined('projfile_raw') ) THROW_HARD('map_params_from_den requires projfile_raw')
        if( .not. cline%defined('projfile_den') ) THROW_HARD('map_params_from_den requires projfile_den')
        call params%new(cline)
        if( params%projfile_raw .eq. '' ) THROW_HARD('map_params_from_den requires projfile_raw')
        if( params%projfile_den .eq. '' ) THROW_HARD('map_params_from_den requires projfile_den')
        call raw_proj%read(params%projfile_raw)
        call den_proj%read(params%projfile_den)
        nptcls = raw_proj%os_ptcl2D%get_noris()
        if( nptcls < 1 ) THROW_HARD('projfile_raw has no ptcl2D records; map_params_from_den')
        nden2d = den_proj%os_ptcl2D%get_noris()
        nden3d = den_proj%os_ptcl3D%get_noris()
        if( nden2d < 1 ) THROW_HARD('projfile_den has no ptcl2D records; map_params_from_den')
        if( nden3d > 0 .and. nden3d /= nden2d )then
            THROW_HARD('projfile_den ptcl2D/ptcl3D particle counts differ; map_params_from_den')
        endif
        call build_den2raw_map(raw_proj, den_proj, den2raw, nraw_active, l_use_pind)
        call outproj%copy(raw_proj)
        if( nden3d == nden2d .and. outproj%os_ptcl3D%get_noris() /= nptcls )then
            outproj%os_ptcl3D = outproj%os_ptcl2D
            call outproj%os_ptcl3D%delete_2Dclustering(keepcls=.true.)
        endif
        nmap2d = 0
        nmap3d = 0
        do iden = 1,nden2d
            iptcl = den2raw(iden)
            call raw_proj%os_ptcl2D%get_ori(iptcl, raw2d)
            call den_proj%os_ptcl2D%get_ori(iden, den_ori)
            call den_ori%compose3d2d(raw2d, mapped_ori)
            call outproj%os_ptcl2D%e3set(iptcl, mapped_ori%e3get())
            call outproj%os_ptcl2D%set_shift(iptcl, mapped_ori%get_2Dshift())
            call copy_assignment_keys(den_proj%os_ptcl2D, outproj%os_ptcl2D, iden, iptcl, .false.)
            nmap2d = nmap2d + 1
            if( nden3d == nden2d )then
                call den_proj%os_ptcl3D%get_ori(iden, den_ori)
                call den_ori%compose3d2d(raw2d, mapped_ori)
                call outproj%os_ptcl3D%e1set(iptcl, mapped_ori%e1get())
                call outproj%os_ptcl3D%e2set(iptcl, mapped_ori%e2get())
                call outproj%os_ptcl3D%e3set(iptcl, mapped_ori%e3get())
                call outproj%os_ptcl3D%set_shift(iptcl, mapped_ori%get_2Dshift())
                call copy_assignment_keys(den_proj%os_ptcl3D, outproj%os_ptcl3D, iden, iptcl, .true.)
                nmap3d = nmap3d + 1
            endif
        end do
        call outproj%update_projinfo(params%projfile)
        call outproj%write(params%projfile)
        write(logfhandle,'(A,A)') 'Mapped denoised-project assignments written: ', params%projfile%to_char()
        write(logfhandle,'(A,I10,A,I10)') 'Mapped ptcl2D records: ', nmap2d, ' ptcl3D records: ', nmap3d
        write(logfhandle,'(A,I10,A,A)') 'Raw active particles: ', nraw_active, ' mapping: ', merge('pind','row ',l_use_pind)
        call raw2d%kill
        call den_ori%kill
        call mapped_ori%kill
        call raw_proj%kill
        call den_proj%kill
        call outproj%kill
        if( allocated(den2raw) ) deallocate(den2raw)
        call simple_end('**** SIMPLE_MAP_PARAMS_FROM_DEN NORMAL STOP ****')

        contains

            subroutine build_den2raw_map( raw_proj, den_proj, den2raw, nraw_active, l_use_pind )
                type(sp_project), intent(in)  :: raw_proj, den_proj
                integer, allocatable, intent(out) :: den2raw(:)
                integer, intent(out) :: nraw_active
                logical, intent(out) :: l_use_pind
                integer, allocatable :: pind2raw(:)
                logical, allocatable :: raw_seen(:)
                integer :: nraw, nden, i, raw_pind, den_pind, raw_ind, max_pind, nmiss
                logical :: raw_has_pind, den_has_pind
                nraw = raw_proj%os_ptcl2D%get_noris()
                nden = den_proj%os_ptcl2D%get_noris()
                allocate(den2raw(nden), source=0)
                allocate(raw_seen(nraw), source=.false.)
                raw_seen = .false.
                raw_has_pind = raw_proj%os_ptcl2D%isthere('pind')
                den_has_pind = den_proj%os_ptcl2D%isthere('pind')
                nraw_active = 0
                do i = 1,nraw
                    if( raw_proj%os_ptcl2D%get_state(i) > 0 ) nraw_active = nraw_active + 1
                enddo
                if( den_has_pind )then
                    l_use_pind = .true.
                    max_pind = 0
                    do i = 1,nraw
                        raw_pind = native_row_pind(raw_proj%os_ptcl2D, i, raw_has_pind, 'projfile_raw')
                        max_pind = max(max_pind, raw_pind)
                    enddo
                    do i = 1,nden
                        den_pind = native_row_pind(den_proj%os_ptcl2D, i, .true., 'projfile_den')
                        max_pind = max(max_pind, den_pind)
                    enddo
                    if( max_pind < 1 ) THROW_HARD('non-positive native particle index; map_params_from_den')
                    allocate(pind2raw(max_pind), source=0)
                    do i = 1,nraw
                        raw_pind = native_row_pind(raw_proj%os_ptcl2D, i, raw_has_pind, 'projfile_raw')
                        if( pind2raw(raw_pind) /= 0 ) THROW_HARD('duplicate native particle index in projfile_raw; map_params_from_den')
                        pind2raw(raw_pind) = i
                    enddo
                    do i = 1,nden
                        den_pind = native_row_pind(den_proj%os_ptcl2D, i, .true., 'projfile_den')
                        if( den_pind > size(pind2raw) .or. pind2raw(den_pind) == 0 )then
                            THROW_HARD('projfile_den native particle index absent from projfile_raw; map_params_from_den')
                        endif
                        raw_ind = pind2raw(den_pind)
                        if( raw_proj%os_ptcl2D%get_state(raw_ind) <= 0 )then
                            THROW_HARD('projfile_den maps to a state=0 raw particle; map_params_from_den')
                        endif
                        if( raw_seen(raw_ind) ) THROW_HARD('duplicate projfile_den mapping to raw particle; map_params_from_den')
                        raw_seen(raw_ind) = .true.
                        den2raw(i) = raw_ind
                    enddo
                    nmiss = 0
                    do i = 1,nraw
                        if( raw_proj%os_ptcl2D%get_state(i) > 0 .and. .not.raw_seen(i) ) nmiss = nmiss + 1
                    enddo
                    if( nmiss > 0 )then
                        write(logfhandle,'(A,I10,A,I10)') 'Raw active particles missing denoised assignments: ', nmiss, &
                            ' raw_active=', nraw_active
                        THROW_HARD('projfile_den does not cover every active projfile_raw particle; map_params_from_den')
                    endif
                    deallocate(pind2raw)
                else if( nden == nraw )then
                    l_use_pind = .false.
                    do i = 1,nden
                        den2raw(i) = i
                    enddo
                else
                    write(logfhandle,'(A,I10,A,I10)') 'projfile_raw nptcls=', nraw, ' projfile_den nptcls=', nden
                    THROW_HARD('compressed projfile_den requires ptcl2D:pind native particle indices; map_params_from_den')
                endif
                deallocate(raw_seen)
            end subroutine build_den2raw_map

            integer function native_row_pind( os, i, require_pind, context ) result(pind)
                use simple_oris, only: oris
                class(oris),      intent(in) :: os
                integer,          intent(in) :: i
                logical,          intent(in) :: require_pind
                character(len=*), intent(in) :: context
                if( os%isthere(i, 'pind') )then
                    pind = os%get_int(i, 'pind')
                else if( require_pind )then
                    THROW_HARD(trim(context)//' row missing pind; map_params_from_den')
                else
                    pind = i
                endif
                if( pind < 1 ) THROW_HARD(trim(context)//' has non-positive pind; map_params_from_den')
            end function native_row_pind

            subroutine copy_assignment_keys( src, dst, isrc, idst, include_proj )
                use simple_oris, only: oris
                class(oris), intent(in)    :: src
                class(oris), intent(inout) :: dst
                integer,     intent(in)    :: isrc, idst
                logical,     intent(in)    :: include_proj
                if( src%isthere(isrc, 'class')     ) call dst%set(idst, 'class',     src%get_int(isrc, 'class'))
                if( src%isthere(isrc, 'state')     ) call dst%set(idst, 'state',     src%get_int(isrc, 'state'))
                if( src%isthere(isrc, 'corr')      ) call dst%set(idst, 'corr',      src%get(isrc, 'corr'))
                if( src%isthere(isrc, 'updatecnt') ) call dst%set(idst, 'updatecnt', src%get_int(isrc, 'updatecnt'))
                if( src%isthere(isrc, 'sampled')   ) call dst%set(idst, 'sampled',   src%get_int(isrc, 'sampled'))
                if( include_proj )then
                    if( src%isthere(isrc, 'proj') ) call dst%set(idst, 'proj', src%get_int(isrc, 'proj'))
                endif
            end subroutine copy_assignment_keys

    end subroutine exec_map_params_from_den

end module simple_commanders_denoise
