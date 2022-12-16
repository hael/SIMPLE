! array allocation for concrete strategy3D extensions to improve caching and reduce alloc overheads
module simple_strategy3D_alloc
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_builder,    only: build_glob
implicit none

public :: s3D, clean_strategy3D, prep_strategy3D, prep_strategy3D_thread
private

type strategy3D_alloc
    ! global parameters
    integer,        allocatable :: proj_space_nnmat(:,:)    !< 3 nearest neighbours per reference + self
    ! per-ptcl/ref allocation
    integer,        allocatable :: proj_space_state(:)      !< states
    integer,        allocatable :: proj_space_proj(:)       !< projection directions (1 state assumed)
    logical,        allocatable :: state_exists(:)          !< indicates state existence
    ! per thread allocation
    type(ran_tabu), allocatable :: rts(:)                   !< stochastic search order generators
    type(ran_tabu), allocatable :: rts_sub(:)               !< stochastic search order generators, subspace
    real,           allocatable :: proj_space_euls(:,:,:)   !< euler angles
    real,           allocatable :: proj_space_shift(:,:,:)  !< shift vectors
    real,           allocatable :: proj_space_corrs(:,:)    !< reference vs. particle correlations
    integer,        allocatable :: proj_space_inplinds(:,:) !< in-plane indices
    integer,        allocatable :: srch_order(:,:)          !< stochastic search index
    integer,        allocatable :: srch_order_sub(:,:)      !< stochastic search index, subspace
    logical,        allocatable :: proj_space_mask(:,:)     !< reference vs. particle correlations mask
end type strategy3D_alloc

type(strategy3D_alloc) :: s3D                         ! singleton
real,      allocatable :: master_proj_space_euls(:,:) !< references euler angles
logical                :: srch_order_allocated = .false.

contains

    subroutine prep_strategy3D( ptcl_mask )
        logical, target, intent(in) :: ptcl_mask(params_glob%fromp:params_glob%top)
        integer :: istate, iproj, ithr, cnt, nrefs, nrefs_sub
        real    :: areal

        print *, 'inside prep_strategy3D'

        ! clean all class arrays & types

        print *, 'cleaning strategy3D'

        call clean_strategy3D()
        ! parameters
        nrefs     = params_glob%nspace     * params_glob%nstates
        nrefs_sub = params_glob%nspace_sub * params_glob%nstates

        print *, 'cleaning strategy3D, DONE'

        ! shared-memory arrays

        print *, 'allocating shared memory arrays'

        allocate(master_proj_space_euls(3,nrefs), s3D%proj_space_euls(3,nrefs,nthr_glob),&
            &s3D%proj_space_shift(2,nrefs,nthr_glob), s3D%proj_space_state(nrefs),&
            &s3D%proj_space_corrs(nthr_glob,nrefs),s3D%proj_space_mask(nrefs,nthr_glob),&
            &s3D%proj_space_inplinds(nthr_glob,nrefs),s3D%proj_space_nnmat(4,params_glob%nspace),&
            &s3D%proj_space_proj(nrefs))

        print *, 'allocating shared memory arrays, DONE'

        ! states existence
        if( .not.build_glob%spproj%is_virgin_field(params_glob%oritype) )then
            if( str_has_substr(params_glob%refine,'greedy') )then
                allocate(s3D%state_exists(params_glob%nstates), source=.true.)
            else
                s3D%state_exists = build_glob%spproj_field%states_exist(params_glob%nstates)
            endif
        else
            allocate(s3D%state_exists(params_glob%nstates), source=.true.)
        endif
        ! reference projection directions & state

        print *, 'setting reference projection directions & state '

        cnt = 0
        do istate=1,params_glob%nstates
            do iproj=1,params_glob%nspace
                cnt = cnt + 1
                s3D%proj_space_state(cnt)     = istate
                s3D%proj_space_proj(cnt)      = iproj
                master_proj_space_euls(:,cnt) = build_glob%eulspace%get_euler(iproj)
            enddo
        enddo

        print *, 'setting reference projection directions & state, DONE'

        s3D%proj_space_shift = 0.
        s3D%proj_space_corrs = -HUGE(areal)
        s3D%proj_space_mask  = .false.
        s3D%proj_space_euls  = 0.
        ! search orders allocation
        select case( trim(params_glob%refine) )
            case( 'cluster','clustersym','clustersoft')
                srch_order_allocated = .false.
            case DEFAULT

                print *, 'allocating srch_orders'

                allocate(s3D%srch_order(nthr_glob,nrefs), s3D%srch_order_sub(nthr_glob,nrefs_sub),&
                &s3D%rts(nthr_glob), s3D%rts_sub(nthr_glob))
                do ithr=1,nthr_glob
                    s3D%rts(ithr)     = ran_tabu(nrefs)
                    s3D%rts_sub(ithr) = ran_tabu(nrefs_sub)
                end do
                srch_order_allocated = .true.
                s3D%srch_order       = 0
                s3D%srch_order_sub   = 0

                print *, 'allocating srch_orders, DONE'

        end select
        ! precalculate neaarest neighbour matrix

        print *, 'precalculating proj_space_nnmat'

        call build_glob%eulspace%nearest_proj_neighbors(4, s3D%proj_space_nnmat) ! 4 because self is included

        print *, 'precalculating proj_space_nnmat, DONE'

    end subroutine prep_strategy3D

    ! init thread specific search arrays
    subroutine prep_strategy3D_thread( ithr )
        integer, intent(in)    :: ithr
        real(sp)               :: areal
        s3D%proj_space_euls(:,:,ithr)   = master_proj_space_euls
        s3D%proj_space_corrs(ithr,:)    = -HUGE(areal)
        s3D%proj_space_mask(:,ithr)     = .false.
        s3D%proj_space_shift(:,:,ithr)  = 0.
        s3D%proj_space_inplinds(ithr,:) = 0
        if(srch_order_allocated)then
            s3D%srch_order(ithr,:)     = 0
            s3D%srch_order_sub(ithr,:) = 0
        endif
    end subroutine prep_strategy3D_thread

    subroutine clean_strategy3D
        integer :: ithr
        if( allocated(master_proj_space_euls)  ) deallocate(master_proj_space_euls)
        if( allocated(s3D%state_exists)        ) deallocate(s3D%state_exists)
        if( allocated(s3D%proj_space_state)    ) deallocate(s3D%proj_space_state)
        if( allocated(s3D%proj_space_proj)     ) deallocate(s3D%proj_space_proj)
        if( allocated(s3D%srch_order)          ) deallocate(s3D%srch_order)
        if( allocated(s3D%srch_order_sub)      ) deallocate(s3D%srch_order_sub)
        if( allocated(s3D%proj_space_euls)     ) deallocate(s3D%proj_space_euls)
        if( allocated(s3D%proj_space_shift)    ) deallocate(s3D%proj_space_shift)
        if( allocated(s3D%proj_space_corrs)    ) deallocate(s3D%proj_space_corrs)
        if( allocated(s3D%proj_space_mask)     ) deallocate(s3D%proj_space_mask)
        if( allocated(s3D%proj_space_inplinds) ) deallocate(s3D%proj_space_inplinds)
        if( allocated(s3D%proj_space_nnmat)    ) deallocate(s3D%proj_space_nnmat)
        if( allocated(s3D%rts) )then
            do ithr=1,nthr_glob
                call s3D%rts(ithr)%kill
                call s3D%rts_sub(ithr)%kill
            enddo
            deallocate(s3D%rts, s3D%rts_sub)
        endif
        srch_order_allocated = .false.
    end subroutine clean_strategy3D

end module simple_strategy3D_alloc
