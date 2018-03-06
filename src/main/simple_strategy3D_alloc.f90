! array allocation for concrete strategy3D extensions to improve caching and reduce alloc overheads
module simple_strategy3D_alloc
#include "simple_lib.f08"
use simple_oris, only: oris
implicit none

type(oris), allocatable :: o_peaks(:)              !< solution objects
real,       allocatable :: proj_space_euls(:,:,:)  !< euler angles
real,       allocatable :: proj_space_shift(:,:,:) !< shift vector
real,       allocatable :: proj_space_corrs(:,:)   !< correlations vs. reference orientations
integer,    allocatable :: proj_space_inds(:,:)    !< stochastic index of reference orientations
integer,    allocatable :: proj_space_state(:,:)   !< reference orientations state
integer,    allocatable :: proj_space_proj(:,:)    !< reference orientations projection direction (1 state assumed)
integer,    allocatable :: prev_proj(:)            !< particle previous reference projection direction
integer,    allocatable :: srch_order(:,:)         !< stochastic search index
logical,    allocatable :: state_exists(:)         !< indicates state existence

logical, private, parameter :: DEBUG = .false.

contains

    subroutine clean_strategy3D
        if( allocated(proj_space_euls)  ) deallocate(proj_space_euls)
        if( allocated(proj_space_shift) ) deallocate(proj_space_shift)
        if( allocated(proj_space_corrs) ) deallocate(proj_space_corrs)
        if( allocated(proj_space_inds)  ) deallocate(proj_space_inds)
        if( allocated(proj_space_state) ) deallocate(proj_space_state)
        if( allocated(proj_space_proj)  ) deallocate(proj_space_proj)
        if( allocated(prev_proj)        ) deallocate(prev_proj)
        if( allocated(srch_order)       ) deallocate(srch_order)
        if( allocated(state_exists)     ) deallocate(state_exists)
    end subroutine clean_strategy3D

    subroutine prep_strategy3D( b, p, ptcl_mask )
        use simple_ran_tabu, only: ran_tabu
        use simple_params,   only: params
        use simple_build,    only: build
        class(build),    intent(inout) :: b
        class(params),   intent(inout) :: p
        logical, target, intent(in)    :: ptcl_mask(p%fromp:p%top)
        type(ran_tabu), allocatable :: rts(:)
        type(ran_tabu) :: rt
        integer        :: pinds(p%fromp:p%top)
        integer        :: i, istate, iproj, iptcl, prev_state
        integer        :: nnnrefs, cnt, prev_ref, nrefs, nptcls
        ! clean all class arrays & types
        call clean_strategy3D()
        if( allocated(o_peaks) )then
            do i=p%fromp,p%top
                call o_peaks(i)%kill
            enddo
            deallocate(o_peaks)
        endif
        ! parameters
        nrefs  = p%nspace * p%nstates
        nptcls = count(ptcl_mask)
        ! particle index mapping
        pinds = 0
        cnt   = 0
        do i=p%fromp,p%top
            if( ptcl_mask(i) )then
                cnt = cnt + 1
                pinds(i) = cnt
            endif
        end do
        ! shared-memory arrays
        allocate(proj_space_euls(nptcls,nrefs,3), proj_space_shift(nptcls,nrefs,2),&
            &proj_space_corrs(nptcls,nrefs), proj_space_inds(nptcls,nrefs),&
            &proj_space_state(nptcls,nrefs), proj_space_proj(nptcls,nrefs),&
            &prev_proj(nptcls) )
        ! states existence
        if( p%oritab.ne.'' )then
            state_exists = b%a%states_exist(p%nstates)
        else
            allocate(state_exists(p%nstates), source=.true.)
        endif
        ! The shared memory used in a parallel section should be initialised
        ! with a (redundant) parallel section, because of how pages are organised.
        ! Memory otherwise becomes associated with the single thread used for
        ! allocation, causing load imbalance. This will reduce cache misses.
        !$omp parallel default(shared) private(i,iptcl,cnt,istate,iproj) proc_bind(close)
        !$omp do schedule(static)
        do i=1,nptcls
            proj_space_shift(i,:,:) = 0.
            proj_space_corrs(i,:)   = -1.
            proj_space_inds( i,:)   = 0
            proj_space_state(i,:)   = 0
            proj_space_proj( i,:)   = 0
            ! reference projection directions
            cnt = 0
            do istate=1,p%nstates
                do iproj=1,p%nspace
                    cnt = cnt + 1
                    proj_space_state(i,cnt)   = istate
                    proj_space_proj(i, cnt)   = iproj
                    proj_space_euls(i, cnt,:) = b%e%get_euler(iproj)
                enddo
            enddo
        end do
        !$omp end do nowait
        !$omp do schedule(static)
        do iptcl = p%fromp, p%top
            if( pinds(iptcl) > 0 ) prev_proj(pinds(iptcl)) = b%e%find_closest_proj(b%a%get_ori(iptcl))
        enddo
        !$omp end do nowait
        !$omp end parallel
        ! projection direction peaks, eo & CTF transfer
        allocate(o_peaks(p%fromp:p%top))
        do iptcl = p%fromp, p%top
            if( ptcl_mask(iptcl) )then
                ! orientation peaks
                call o_peaks(iptcl)%new_clean(p%npeaks)
                ! transfer CTF params
                if( p%ctf.ne.'no' )then
                    call o_peaks(iptcl)%set_all2single('kv',    b%a%get(iptcl,'kv')   )
                    call o_peaks(iptcl)%set_all2single('cs',    b%a%get(iptcl,'cs')   )
                    call o_peaks(iptcl)%set_all2single('fraca', b%a%get(iptcl,'fraca'))
                    call o_peaks(iptcl)%set_all2single('dfx',   b%a%get(iptcl,'dfx')  )
                    if( p%tfplan%mode .eq. 'astig' )then
                        call o_peaks(iptcl)%set_all2single('dfy',    b%a%get(iptcl,'dfy')   )
                        call o_peaks(iptcl)%set_all2single('angast', b%a%get(iptcl,'angast'))
                    else
                        call o_peaks(iptcl)%set_all2single('dfy',    b%a%get(iptcl,'dfx'))
                        call o_peaks(iptcl)%set_all2single('angast', 0.)
                    endif
                    if( p%tfplan%l_phaseplate )then
                        call o_peaks(iptcl)%set_all2single('phshift', b%a%get(iptcl,'phshift'))
                    endif
                endif
                ! transfer eo flag
                call o_peaks(iptcl)%set_all2single('eo', b%a%get(iptcl,'eo'))
            endif
        enddo
        ! refine mode specific allocations and initialisations
        select case( trim(p%refine) )
            case( 'cluster','clustersym','clusterdev' )
                ! nothing to do
            case DEFAULT
                if( p%neigh.eq.'yes' )then
                    nnnrefs =  p%nnn * p%nstates
                    allocate(srch_order(nptcls, nnnrefs), rts(nptcls))
                    do i=1,nptcls
                        rts(i) = ran_tabu(nnnrefs)
                    end do
                    !$omp parallel do default(shared) private(iptcl,prev_state,prev_ref,istate,i) schedule(static) proc_bind(close)
                    do iptcl = p%fromp, p%top
                        if( pinds(iptcl) > 0 )then
                            prev_state = b%a%get_state(iptcl)
                            prev_ref   = (prev_state - 1) * p%nspace + prev_proj(pinds(iptcl))
                            do istate = 0, p%nstates - 1
                                i = istate * p%nnn + 1
                                srch_order(pinds(iptcl),i:i + p%nnn - 1) = b%nnmat(prev_proj(pinds(iptcl)),:) + istate * p%nspace
                            enddo
                            call rts(pinds(iptcl))%shuffle(srch_order(pinds(iptcl),:))
                            call put_last(prev_ref, srch_order(pinds(iptcl),:))
                        endif
                    enddo
                    !$omp end parallel do
                    do i=1,nptcls
                        call rts(i)%kill
                    end do
                else
                    allocate(srch_order(nptcls,nrefs), rts(nptcls))
                    do i=1,nptcls
                        rts(i) = ran_tabu(nrefs)
                    end do
                    !$omp parallel do default(shared) private(iptcl,prev_state,prev_ref,istate,i) schedule(static) proc_bind(close)
                    do iptcl = p%fromp, p%top
                        if( pinds(iptcl) > 0 )then
                            call rts(pinds(iptcl))%ne_ran_iarr(srch_order(pinds(iptcl),:))
                            prev_state = b%a%get_state(iptcl)
                            prev_ref   = (prev_state - 1) * p%nspace + prev_proj(pinds(iptcl))
                            call put_last(prev_ref, srch_order(pinds(iptcl),:))
                        endif
                    enddo
                    !$omp end parallel do
                    do i=1,nptcls
                        call rts(i)%kill
                    end do
                endif
        end select
        call rt%kill
        if( allocated(srch_order) )then
            if( any(srch_order == 0) ) stop 'Invalid index in srch_order; simple_strategy3D_alloc :: prep_strategy3D'
        endif
    end subroutine prep_strategy3D

end module simple_strategy3D_alloc
