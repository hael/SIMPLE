!@descr: 3D reconstruction and associated things
module simple_commanders_rec
use simple_commanders_api
use simple_matcher_2Dprep
use simple_matcher_3Drec, only: calc_3Drec
use simple_refine3D_fnames, only: refine3D_fsc_fname, refine3D_state_halfvol_fname, refine3D_state_vol_fname
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_rec3D
  contains
    procedure :: execute => exec_rec3D
end type commander_rec3D

type, extends(commander_base) :: commander_bootstrap_rec3D
  contains
    procedure :: execute => exec_bootstrap_rec3D
end type commander_bootstrap_rec3D

type, extends(commander_base) :: commander_rec3D_worker
  contains
    procedure :: execute      => exec_rec3D_distr_worker
end type commander_rec3D_worker

type, extends(commander_base) :: random_rec_commander
  contains
    procedure :: execute      => exec_random_rec
end type random_rec_commander

contains

    subroutine exec_rec3D( self, cline )
        use simple_rec3D_strategy, only: rec3D_strategy, create_rec3D_strategy
        use simple_parameters,     only: parameters
        use simple_builder,        only: builder
        class(commander_rec3D), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        class(rec3D_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        ! Commander-level defaults (apply to both modes)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('trs')     ) call cline%set('trs', 5.)     ! to assure that shifts are being used
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call cline%delete('refine')
        ! Select and run strategy
        strategy = create_rec3D_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params, build, cline)
        ! End gracefully (single unified termination)
        call simple_end('**** SIMPLE_RECONSTRUCT3D NORMAL STOP ****', print_simple=.false.)
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_rec3D

    subroutine exec_bootstrap_rec3D( self, cline )
        class(commander_bootstrap_rec3D), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(commander_rec3D) :: xrec3D
        type(cmdline)         :: cline_unreg, cline_reg
        type(parameters)      :: params
        type(string)          :: sigma_star
        integer               :: state, which_iter
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',       'yes')
        if( .not. cline%defined('oritype')     ) call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('nstates')     ) call cline%set('nstates',          1)
        call cline%set('sigma_est', 'global')
        if( .not. cline%defined('which_iter')  ) call cline%set('which_iter',       1)
        if( .not. cline%defined('postprocess') ) call cline%set('postprocess',  'yes')
        if( .not. cline%defined('combine_eo')  ) call cline%set('combine_eo',    'no')
        call cline%delete('envfsc')
        call cline%delete('objfun')
        call cline%delete('ml_reg')
        call params%new(cline)
        which_iter = max(1, params%which_iter)
        call cline%set('which_iter', which_iter)
        call cline%set('mkdir', 'no') ! child reconstruct3D calls must not create nested run directories
        cline_unreg = cline
        call prepare_bootstrap_rec_cline(cline_unreg, l_regularized=.false.)
        write(logfhandle,'(A)') '>>> BOOTSTRAP_REC3D PASS 1: UNREGULARIZED EVEN/ODD RECONSTRUCTION'
        call xrec3D%execute(cline_unreg)
        call write_bootstrap_sigma2_from_halfmaps(sigma_star)
        write(logfhandle,'(A,1X,A)') '>>> BOOTSTRAP_REC3D WROTE SIGMA STAR:', sigma_star%to_char()
        cline_reg = cline
        call prepare_bootstrap_rec_cline(cline_reg, l_regularized=.true.)
        write(logfhandle,'(A)') '>>> BOOTSTRAP_REC3D PASS 2: EUCLID ML-REGULARIZED RECONSTRUCTION'
        call xrec3D%execute(cline_reg)
        call register_bootstrap_rec_outputs()
        do state = 1, params%nstates
            call cline%set('vol'//int2str(state), refine3D_state_vol_fname(state))
        enddo
        call sigma_star%kill
        call cline_unreg%kill
        call cline_reg%kill
        call simple_end('**** SIMPLE_BOOTSTRAP_REC3D NORMAL STOP ****', print_simple=.false.)

    contains

        subroutine prepare_bootstrap_rec_cline( cline_rec, l_regularized )
            class(cmdline), intent(inout) :: cline_rec
            logical,        intent(in)    :: l_regularized
            integer :: state
            call cline_rec%set('prg',       'reconstruct3D')
            call cline_rec%set('mkdir',              'no')
            call cline_rec%set('oritype', params%oritype)
            call cline_rec%set('nstates', params%nstates)
            call cline_rec%set('sigma_est',     'global')
            call cline_rec%set('which_iter', which_iter)
            call cline_rec%set('trail_rec',        'no')
            call cline_rec%set('combine_eo',       'no')
            call cline_rec%delete('refine')
            call cline_rec%delete('update_frac')
            call cline_rec%delete('fillin')
            call cline_rec%delete('ufrac_trec')
            call cline_rec%delete('endit')
            call cline_rec%delete('vol_even')
            call cline_rec%delete('vol_odd')
            call cline_rec%delete('refs')
            call cline_rec%delete('refs_even')
            call cline_rec%delete('refs_odd')
            do state = 1, params%nstates
                call cline_rec%delete('vol'//int2str(state))
            enddo
            if( l_regularized )then
                call cline_rec%set('objfun', 'euclid')
                call cline_rec%set('ml_reg',    'yes')
            else
                call cline_rec%set('objfun',       'cc')
                call cline_rec%set('ml_reg',       'no')
                call cline_rec%set('postprocess',  'no')
            endif
        end subroutine prepare_bootstrap_rec_cline

        subroutine register_bootstrap_rec_outputs()
            type(sp_project) :: spproj
            type(string)     :: volname, fscname
            integer          :: state, pop
            character(len=8) :: imgkind
            call spproj%read_segment('out', params%projfile)
            call spproj%read_segment(params%oritype, params%projfile)
            select case(trim(params%oritype))
                case('cls3D')
                    imgkind = 'vol_cavg'
                case DEFAULT
                    imgkind = 'vol'
            end select
            do state = 1, params%nstates
                select case(trim(params%oritype))
                    case('cls3D')
                        pop = spproj%os_cls3D%get_pop(state, 'state')
                    case DEFAULT
                        pop = spproj%os_ptcl3D%get_pop(state, 'state')
                end select
                if( pop == 0 )cycle
                volname = refine3D_state_vol_fname(state)
                if( .not. file_exists(volname) )then
                    call volname%kill
                    cycle
                endif
                fscname = refine3D_fsc_fname(state)
                call spproj%add_vol2os_out(volname, params%smpd_crop, state, trim(imgkind), pop=pop)
                if( file_exists(fscname) ) call spproj%add_fsc2os_out(fscname, state, params%box_crop)
                call volname%kill
                call fscname%kill
            enddo
            call spproj%write_segment_inside('out', params%projfile)
            call spproj%kill
        end subroutine register_bootstrap_rec_outputs

        subroutine write_bootstrap_sigma2_from_halfmaps( sigma_star )
            type(string), intent(out) :: sigma_star
            type(image) :: vol_even, vol_odd, vol_noise
            type(string) :: fname_even, fname_odd
            real, allocatable :: pspec(:), sigma2(:), group_pspecs(:,:,:), state_scales(:), state_weights(:)
            real :: total_weight
            integer :: ldim(3), nptcls, state, nstates_used, last_shell, nspec, filtsz
            filtsz = fdim(params%box) - 1
            if( filtsz < 1 ) THROW_HARD('Invalid box for bootstrap_rec3D sigma estimation')
            allocate(sigma2(filtsz), source=0.)
            call calc_bootstrap_state_scales(state_scales, state_weights)
            nstates_used = 0
            last_shell   = 0
            total_weight = 0.
            do state = 1, params%nstates
                fname_even = refine3D_state_halfvol_fname(state, 'even')
                fname_odd  = refine3D_state_halfvol_fname(state, 'odd')
                if( (.not. file_exists(fname_even)) .or. (.not. file_exists(fname_odd)) )then
                    call fname_even%kill
                    call fname_odd%kill
                    cycle
                endif
                call find_ldim_nptcls(fname_even, ldim, nptcls)
                call vol_even%new(ldim, params%smpd_crop)
                call vol_odd%new(ldim,  params%smpd_crop)
                call vol_even%read(fname_even)
                call vol_odd%read(fname_odd)
                call vol_noise%copy(vol_even)
                call vol_noise%subtr(vol_odd)
                call vol_noise%spectrum('power', pspec, norm=.true.)
                nspec = min(size(pspec), filtsz)
                if( nspec > 0 )then
                    sigma2(1:nspec) = sigma2(1:nspec) + state_weights(state) * state_scales(state) * pspec(1:nspec)
                    last_shell      = max(last_shell, nspec)
                    nstates_used    = nstates_used + 1
                    total_weight    = total_weight + state_weights(state)
                endif
                if( allocated(pspec) ) deallocate(pspec)
                call vol_even%kill
                call vol_odd%kill
                call vol_noise%kill
                call fname_even%kill
                call fname_odd%kill
            enddo
            if( nstates_used == 0 ) THROW_HARD('No even/odd half maps found after bootstrap_rec3D unregularized reconstruction')
            if( total_weight > TINY )then
                sigma2 = sigma2 / total_weight
            else
                sigma2 = sigma2 / real(nstates_used)
            endif
            if( last_shell < filtsz .and. last_shell > 0 ) sigma2(last_shell+1:filtsz) = sigma2(last_shell)
            call condition_bootstrap_sigma2(sigma2, last_shell)
            allocate(group_pspecs(2,1,1:filtsz))
            group_pspecs(1,1,:) = sigma2
            group_pspecs(2,1,:) = sigma2
            sigma_star = sigma2_star_from_iter(which_iter)
            call write_groups_starfile(sigma_star, group_pspecs, 1)
            deallocate(group_pspecs, sigma2, state_scales, state_weights)
        end subroutine write_bootstrap_sigma2_from_halfmaps

        subroutine calc_bootstrap_state_scales( state_scales, state_weights )
            real, allocatable, intent(out) :: state_scales(:)
            real, allocatable, intent(out) :: state_weights(:)
            type(sp_project) :: spproj
            real, allocatable :: eo_weights(:,:)
            real :: denom
            integer :: state
            allocate(state_scales(params%nstates), state_weights(params%nstates), eo_weights(2,params%nstates), source=0.)
            select case(trim(params%oritype))
                case('ptcl3D')
                    call spproj%read_segment('ptcl3D', params%projfile)
                    call accumulate_bootstrap_eo_weights(spproj%os_ptcl3D, eo_weights)
                case('cls3D')
                    call spproj%read_segment('cls3D', params%projfile)
                    call accumulate_bootstrap_eo_weights(spproj%os_cls3D, eo_weights)
                case DEFAULT
                    THROW_HARD('bootstrap_rec3D supports ptcl3D and cls3D orientations only')
            end select
            do state = 1, params%nstates
                denom = eo_weights(1,state) + eo_weights(2,state)
                if( eo_weights(1,state) > TINY .and. eo_weights(2,state) > TINY )then
                    state_scales(state) = (eo_weights(1,state) * eo_weights(2,state)) / denom
                    state_weights(state) = denom
                else
                    state_scales(state) = 0.5
                    state_weights(state) = 1.
                    write(logfhandle,'(A,I0,A)') '>>> WARNING: BOOTSTRAP_REC3D COULD NOT DETERMINE EVEN/ODD WEIGHTS FOR STATE ',&
                        &state, '; USING UNSCALED HALF-MAP DIFFERENCE'
                endif
                write(logfhandle,'(A,I0,A,F12.3)') '>>> BOOTSTRAP_REC3D SIGMA SCALE STATE ', state, ':', state_scales(state)
            enddo
            call spproj%kill
            deallocate(eo_weights)
        end subroutine calc_bootstrap_state_scales

        subroutine accumulate_bootstrap_eo_weights( os, eo_weights )
            class(oris), intent(in)    :: os
            real,        intent(inout) :: eo_weights(:,:)
            real    :: w
            integer :: eo, iptcl, state
            logical :: l_have_weights
            l_have_weights = os%isthere('w')
            do iptcl = 1, os%get_noris()
                state = os%get_state(iptcl)
                if( state < 1 .or. state > size(eo_weights,2) )cycle
                eo = os%get_eo(iptcl) + 1
                if( eo < 1 .or. eo > 2 )cycle
                if( l_have_weights )then
                    w = max(0., os%get(iptcl, 'w'))
                else
                    w = 1.
                endif
                if( w < TINY )cycle
                eo_weights(eo,state) = eo_weights(eo,state) + w
            enddo
        end subroutine accumulate_bootstrap_eo_weights

        subroutine condition_bootstrap_sigma2( sigma2, last_shell )
            real,    intent(inout) :: sigma2(:)
            integer, intent(in)    :: last_shell
            real, allocatable :: logsigma(:), smoothed(:)
            real    :: floor_val, pos_min
            integer :: i, pass, n, anchor
            n = size(sigma2)
            if( n == 0 )return
            pos_min = huge(pos_min)
            do i = 1, n
                if( sigma2(i) > TINY ) pos_min = min(pos_min, sigma2(i))
            enddo
            if( pos_min == huge(pos_min) ) THROW_HARD('Unable to estimate positive sigma2 values from even/odd half-map difference')
            floor_val = max(pos_min * 1.e-3, TINY)
            do i = 1, n
                if( sigma2(i) <= floor_val ) sigma2(i) = floor_val
            enddo
            if( last_shell > 0 .and. last_shell < n ) sigma2(last_shell+1:n) = sigma2(last_shell)
            anchor = min(n, 6)
            if( anchor > 1 ) sigma2(1:anchor-1) = max(sigma2(1:anchor-1), sigma2(anchor))
            allocate(logsigma(n), smoothed(n))
            logsigma = log(max(sigma2, floor_val))
            do pass = 1, 2
                smoothed = logsigma
                do i = 2, n - 1
                    smoothed(i) = 0.25 * logsigma(i-1) + 0.5 * logsigma(i) + 0.25 * logsigma(i+1)
                enddo
                logsigma = smoothed
            enddo
            sigma2 = exp(logsigma)
            deallocate(logsigma, smoothed)
        end subroutine condition_bootstrap_sigma2
    end subroutine exec_bootstrap_rec3D

    subroutine exec_rec3D_distr_worker( self, cline )
        class(commander_rec3D_worker), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)     :: params
        type(builder)        :: build
        type(string)         :: fname
        integer, allocatable :: pinds(:)
        integer              :: nptcls2update
        call build%init_params_and_build_general_tbox(cline, params)
        call build%build_strategy3D_tbox(params)
        if( params%l_update_frac .and. build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls2update, pinds)
        else
            ! we sample all state > 0 and updatecnt > 0
            call build%spproj_field%sample4rec([params%fromp,params%top], nptcls2update, pinds)
        endif
        if( params%l_ml_reg )then
            fname = SIGMA2_FBODY//int2str_pad(params%part,params%numlen)//'.dat'
            call build%esig%new(params, build%pftc, fname, params%box)
            call build%esig%read_groups(build%spproj_field)
        end if
        call calc_3Drec( params, build, cline, nptcls2update, pinds )
        ! cleanup
        call build%esig%kill
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_rec :: exec_rec3D'))
    end subroutine exec_rec3D_distr_worker

    subroutine exec_random_rec( self, cline )
        class(random_rec_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(commander_rec3D) :: xrec3D
        type(parameters)      :: params
        type(builder)         :: build
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes'   )
        call build%init_params_and_build_spproj(cline, params)
        call build%spproj%os_ptcl3D%rnd_oris
        call build%spproj%write_segment_inside('ptcl3D', params%projfile)
        call cline%set('mkdir', 'no') ! to avoid nested dirs
        call cline%set('prg',   'rec3D')
        call xrec3D%execute(cline)
        call build%spproj_field%kill
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_RANDOM_REC NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_random_rec

end module simple_commanders_rec
