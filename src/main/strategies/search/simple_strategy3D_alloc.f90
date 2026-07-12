!@descr: array allocation for concrete strategy3D extensions to improve caching and reduce alloc overheads
module simple_strategy3D_alloc
use simple_pftc_srch_api
use simple_builder, only: builder
implicit none

public :: s3D, clean_strategy3D, prep_strategy3D, prep_strategy3D_thread
public :: ref_state_from_index, ref_proj_from_index
private

type strategy3D_alloc
    ! global parameters
    real,           allocatable :: smpl_refs_athres(:)      !< refine=smpl; angular threshold of projections directions
    real,           allocatable :: smpl_inpl_athres(:)      !< refine=smpl; angular threshold of in-plane rotations
    logical,        allocatable :: state_exists(:)          !< indicates state existence
    ! per thread allocation
    type(ran_tabu), allocatable :: rts(:)                   !< stochastic search order generators
    type(ran_tabu), allocatable :: rts_inpl(:)              !< stochastic search order generators for inpls
    type(ran_tabu), allocatable :: rts_sub(:)               !< stochastic search order generators, subspace
    real,           allocatable :: proj_space_shift(:,:,:)  !< shift vectors
    real,           allocatable :: proj_space_corrs(:,:)    !< reference vs. particle correlations
    integer,        allocatable :: proj_space_inplinds(:,:) !< in-plane indices
    integer,        allocatable :: srch_order(:,:)          !< stochastic search index
    integer,        allocatable :: inpl_order(:,:)          !< stochastic inpl search index
    integer,        allocatable :: srch_order_sub(:,:)      !< stochastic search index, subspace
end type strategy3D_alloc

type(strategy3D_alloc) :: s3D                         ! singleton
logical                :: srch_order_allocated = .false.

contains

    subroutine prep_strategy3D( params, build )
        use simple_eul_prob_tab_utils, only: calc_athres
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer :: istate, ithr, nrefs, nrefs_sub
        real    :: areal
        ! clean all class arrays & types
        call clean_strategy3D()
        ! parameters
        nrefs     = params%nspace     * params%nstates
        nrefs_sub = params%nspace_sub * params%nstates
        ! shared-memory arrays
        allocate(s3D%proj_space_shift(2,nrefs,nthr_glob),&
            &s3D%proj_space_corrs(nrefs,nthr_glob),&
            &s3D%proj_space_inplinds(nrefs,nthr_glob))
        ! states existence
        if( .not.build%spproj%is_virgin_field(params%oritype) )then
            if( str_has_substr(params%refine,'greedy') )then
                allocate(s3D%state_exists(params%nstates), source=.true.)
            else
                s3D%state_exists = build%spproj_field%states_exist(params%nstates)
            endif
        else
            allocate(s3D%state_exists(params%nstates), source=.true.)
        endif
        s3D%proj_space_shift = 0.
        s3D%proj_space_corrs = -HUGE(areal)
        ! search orders allocation
        select case( trim(params%refine) )
            case( 'cluster','clustersym','clustersoft','prob','prob_state','prob_neigh')
                srch_order_allocated = .false.
            case DEFAULT
                allocate(s3D%srch_order(nrefs,nthr_glob), s3D%rts(nthr_glob),&
                    &s3D%rts_inpl(nthr_glob), s3D%inpl_order(build%pftc%get_nrots(),nthr_glob))
                if( params%l_neigh )then
                    allocate(s3D%srch_order_sub(nrefs_sub,nthr_glob), s3D%rts_sub(nthr_glob))
                endif
                do ithr=1,nthr_glob
                    s3D%rts(ithr)      = ran_tabu(nrefs)
                    s3D%rts_inpl(ithr) = ran_tabu(build%pftc%get_nrots())
                    if( params%l_neigh ) s3D%rts_sub(ithr) = ran_tabu(nrefs_sub)
                end do
                srch_order_allocated = .true.
                s3D%srch_order       = 0
                s3D%inpl_order       = 1
                if( allocated(s3D%srch_order_sub) ) s3D%srch_order_sub = 0
        end select
        ! calculate peak thresholds for probabilistic searches
        if( allocated(s3D%smpl_refs_athres) ) deallocate(s3D%smpl_refs_athres)
        if( allocated(s3D%smpl_inpl_athres) ) deallocate(s3D%smpl_inpl_athres)
        allocate(s3D%smpl_refs_athres(params%nstates), s3D%smpl_inpl_athres(params%nstates))
        do istate = 1, params%nstates
            s3D%smpl_refs_athres(istate) = calc_athres(os=build%spproj_field, field_str='dist',      prob_athres=params%prob_athres, state=istate)
            s3D%smpl_inpl_athres(istate) = calc_athres(os=build%spproj_field, field_str='dist_inpl', prob_athres=params%prob_athres, state=istate)
        enddo
    end subroutine prep_strategy3D

    ! init thread specific search arrays
    subroutine prep_strategy3D_thread( ithr )
        integer, intent(in) :: ithr
        real(sp)            :: areal
        s3D%proj_space_corrs(   :,ithr) = -HUGE(areal)
        s3D%proj_space_shift( :,:,ithr) = 0.
        s3D%proj_space_inplinds(:,ithr) = 0
        if(srch_order_allocated)then
            s3D%srch_order(    :,ithr) = 0
            if( allocated(s3D%srch_order_sub) ) s3D%srch_order_sub(:,ithr) = 0
        endif
    end subroutine prep_strategy3D_thread

    subroutine clean_strategy3D
        integer :: ithr
        if( allocated(s3D%state_exists)        ) deallocate(s3D%state_exists)
        if( allocated(s3D%srch_order)          ) deallocate(s3D%srch_order)
        if( allocated(s3D%inpl_order)          ) deallocate(s3D%inpl_order)
        if( allocated(s3D%srch_order_sub)      ) deallocate(s3D%srch_order_sub)
        if( allocated(s3D%proj_space_shift)    ) deallocate(s3D%proj_space_shift)
        if( allocated(s3D%proj_space_corrs)    ) deallocate(s3D%proj_space_corrs)
        if( allocated(s3D%proj_space_inplinds) ) deallocate(s3D%proj_space_inplinds)
        if( allocated(s3D%smpl_refs_athres)    ) deallocate(s3D%smpl_refs_athres)
        if( allocated(s3D%smpl_inpl_athres)    ) deallocate(s3D%smpl_inpl_athres)
        if( allocated(s3D%rts) )then
            do ithr=1,nthr_glob
                call s3D%rts(ithr)%kill
                call s3D%rts_inpl(ithr)%kill
            enddo
            deallocate(s3D%rts, s3D%rts_inpl)
        endif
        if( allocated(s3D%rts_sub) )then
            do ithr=1,nthr_glob
                call s3D%rts_sub(ithr)%kill
            enddo
            deallocate(s3D%rts_sub)
        endif
        srch_order_allocated = .false.
    end subroutine clean_strategy3D

    pure integer function ref_state_from_index( iref, nspace ) result( state )
        integer, intent(in) :: iref, nspace
        state = (iref - 1) / nspace + 1
    end function ref_state_from_index

    pure integer function ref_proj_from_index( iref, nspace ) result( proj )
        integer, intent(in) :: iref, nspace
        proj = modulo(iref - 1, nspace) + 1
    end function ref_proj_from_index

end module simple_strategy3D_alloc
