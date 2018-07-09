! array allocation for concrete strategy3D extensions to improve caching and reduce alloc overheads
module simple_strategy3D_alloc
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_builder,    only: build_glob
use simple_oris,       only: oris
implicit none

public :: s3D, clean_strategy3D, prep_strategy3D
private

type strategy3D_alloc
    type(oris), allocatable :: o_peaks(:)                !< solution objects
    real,       allocatable :: proj_space_euls(:,:,:,:)  !< euler angles
    real,       allocatable :: proj_space_shift(:,:,:,:) !< shift vectors
    real,       allocatable :: proj_space_corrs(:,:,:)   !< reference vs. particle correlations
    integer,    allocatable :: proj_space_refinds(:,:)   !< reference indices
    integer,    allocatable :: proj_space_inplinds(:,:)  !< in-plane indices
    integer,    allocatable :: proj_space_state(:,:)     !< states
    integer,    allocatable :: proj_space_proj(:,:)      !< projection directions (1 state assumed)
    integer,    allocatable :: prev_proj(:)              !< previous projection direction index
    integer,    allocatable :: srch_order(:,:)           !< stochastic search index
    logical,    allocatable :: state_exists(:)           !< indicates state existence
end type strategy3D_alloc

type(strategy3D_alloc) :: s3D ! singleton

contains

    subroutine prep_strategy3D(  ptcl_mask, npeaks )
        logical, target, intent(in) :: ptcl_mask(params_glob%fromp:params_glob%top)
        integer,         intent(in) :: npeaks
        type(ran_tabu), allocatable :: rts(:)
        type(ran_tabu) :: rt
        real    :: eul(3)
        integer :: pinds(params_glob%fromp:params_glob%top)
        integer :: i, istate, iproj, iptcl, prev_state, inpl
        integer :: nnnrefs, cnt, prev_ref, nrefs, nptcls
        ! clean all class arrays & types
        call clean_strategy3D()
        if( allocated(s3D%o_peaks) )then
            do i=params_glob%fromp,params_glob%top
                call s3D%o_peaks(i)%kill
            enddo
            deallocate(s3D%o_peaks)
        endif
        ! parameters
        nrefs  = params_glob%nspace * params_glob%nstates
        nptcls = count(ptcl_mask)
        ! particle index mapping
        pinds = 0
        cnt   = 0
        do i=params_glob%fromp,params_glob%top
            if( ptcl_mask(i) )then
                cnt = cnt + 1
                pinds(i) = cnt
            endif
        end do
        ! shared-memory arrays
        allocate(s3D%proj_space_euls(nptcls,nrefs,MAXNINPLPEAKS,3), s3D%proj_space_shift(nptcls,nrefs,MAXNINPLPEAKS,2),&
            &s3D%proj_space_corrs(nptcls,nrefs,MAXNINPLPEAKS), s3D%proj_space_refinds(nptcls,nrefs),&
            &s3D% proj_space_inplinds(nptcls,MAXNINPLPEAKS), s3D%proj_space_state(nptcls,nrefs),&
            &s3D%proj_space_proj(nptcls,nrefs), s3D%prev_proj(nptcls), stat=alloc_stat )
        if(alloc_stat/=0)call allocchk("strategy3D_alloc failed")
        ! states existence
        if( params_glob%oritab.ne.'' )then
            s3D%state_exists = build_glob%spproj_field%states_exist(params_glob%nstates)
        else
            allocate(s3D%state_exists(params_glob%nstates), source=.true.)
        endif
        ! The shared memory used in a parallel section should be initialised
        ! with a (redundant) parallel section, because of how pages are organised.
        ! Memory otherwise becomes associated with the single thread used for
        ! allocation, causing load imbalance. This will reduce cache misses.
        !$omp parallel default(shared) private(i,iptcl,cnt,istate,iproj,inpl,eul) proc_bind(close)
        !$omp do schedule(static)
        do i=1,nptcls
            s3D%proj_space_shift(i,:,:,:) = 0.
            s3D%proj_space_corrs(i,:,:)   = -1.
            s3D%proj_space_refinds( i,:)     = 0
            s3D%proj_space_state(i,:)     = 0
            s3D%proj_space_proj( i,:)     = 0
            ! reference projection directions
            cnt = 0
            do istate=1,params_glob%nstates
                do iproj=1,params_glob%nspace
                    cnt = cnt + 1
                    s3D%proj_space_state(i,cnt) = istate
                    s3D%proj_space_proj(i,cnt)  = iproj
                    eul = build_glob%eulspace%get_euler(iproj)
                    do inpl=1,MAXNINPLPEAKS
                        s3D%proj_space_euls(i,cnt,inpl,:) = eul
                    end do
                enddo
            enddo
        end do
        !$omp end do nowait
        !$omp do schedule(static)
        do iptcl = params_glob%fromp, params_glob%top
            if( pinds(iptcl) > 0 )s3D%prev_proj(pinds(iptcl)) = build_glob%eulspace%find_closest_proj(build_glob%spproj_field%get_ori(iptcl))
        enddo
        !$omp end do nowait
        !$omp end parallel
        ! projection direction peaks, eo & CTF transfer
        allocate(s3D%o_peaks(params_glob%fromp:params_glob%top))
        do iptcl = params_glob%fromp, params_glob%top
            if( ptcl_mask(iptcl) )then
                ! orientation peaks
                call s3D%o_peaks(iptcl)%new(npeaks * MAXNINPLPEAKS)
                ! transfer CTF params
                if( params_glob%ctf.ne.'no' )then
                    call s3D%o_peaks(iptcl)%set_all2single('kv',     build_glob%spproj_field%get(iptcl,'kv')    )
                    call s3D%o_peaks(iptcl)%set_all2single('cs',     build_glob%spproj_field%get(iptcl,'cs')    )
                    call s3D%o_peaks(iptcl)%set_all2single('fraca',  build_glob%spproj_field%get(iptcl,'fraca') )
                    call s3D%o_peaks(iptcl)%set_all2single('dfx',    build_glob%spproj_field%get(iptcl,'dfx')   )
                    call s3D%o_peaks(iptcl)%set_all2single('dfy',    build_glob%spproj_field%get(iptcl,'dfy')   )
                    call s3D%o_peaks(iptcl)%set_all2single('angast', build_glob%spproj_field%get(iptcl,'angast'))
                    if( params_glob%l_phaseplate )then
                        call s3D%o_peaks(iptcl)%set_all2single('phshift', build_glob%spproj_field%get(iptcl,'phshift'))
                    endif
                endif
                ! transfer eo flag
                call s3D%o_peaks(iptcl)%set_all2single('eo', build_glob%spproj_field%get(iptcl,'eo'))
            endif
        enddo
        ! refine mode specific allocations and initialisations
        select case( trim(params_glob%refine) )
            case( 'cluster','clustersym')
                ! nothing to do
            case DEFAULT
                if( params_glob%neigh.eq.'yes' )then
                    nnnrefs =  params_glob%nnn * params_glob%nstates
                    allocate(s3D%srch_order(nptcls, nnnrefs), rts(nptcls))
                    do i=1,nptcls
                        rts(i) = ran_tabu(nnnrefs)
                    end do
                    !$omp parallel do default(shared) private(iptcl,prev_state,prev_ref,istate,i) schedule(static) proc_bind(close)
                    do iptcl = params_glob%fromp, params_glob%top
                        if( pinds(iptcl) > 0 )then
                            prev_state = build_glob%spproj_field%get_state(iptcl)
                            prev_ref   = (prev_state - 1) * params_glob%nspace + s3D%prev_proj(pinds(iptcl))
                            do istate = 0, params_glob%nstates - 1
                                i = istate * params_glob%nnn + 1
                                s3D%srch_order(pinds(iptcl),i:i + params_glob%nnn - 1) = build_glob%nnmat(s3D%prev_proj(pinds(iptcl)),:) + istate * params_glob%nspace
                            enddo
                            call rts(pinds(iptcl))%shuffle(s3D%srch_order(pinds(iptcl),:))
                            call put_last(prev_ref, s3D%srch_order(pinds(iptcl),:))
                        endif
                    enddo
                    !$omp end parallel do
                    do i=1,nptcls
                        call rts(i)%kill
                    end do
                else
                    allocate(s3D%srch_order(nptcls,nrefs), rts(nptcls))
                    do i=1,nptcls
                        rts(i) = ran_tabu(nrefs)
                    end do
                    !$omp parallel do default(shared) private(iptcl,prev_state,prev_ref,istate,i) schedule(static) proc_bind(close)
                    do iptcl = params_glob%fromp, params_glob%top
                        if( pinds(iptcl) > 0 )then
                            call rts(pinds(iptcl))%ne_ran_iarr(s3D%srch_order(pinds(iptcl),:))
                            prev_state = build_glob%spproj_field%get_state(iptcl)
                            prev_ref   = (prev_state - 1) * params_glob%nspace + s3D%prev_proj(pinds(iptcl))
                            call put_last(prev_ref, s3D%srch_order(pinds(iptcl),:))
                        endif
                    enddo
                    !$omp end parallel do
                    do i=1,nptcls
                        call rts(i)%kill
                    end do
                endif
        end select
        call rt%kill
        if( allocated(s3D%srch_order) )then
            if( any(s3D%srch_order == 0) ) stop 'Invalid index in srch_order; simple_strategy3D_alloc :: prep_strategy3D'
        endif
    end subroutine prep_strategy3D

    subroutine clean_strategy3D
        if( allocated(s3D%proj_space_euls)     ) deallocate(s3D%proj_space_euls)
        if( allocated(s3D%proj_space_shift)    ) deallocate(s3D%proj_space_shift)
        if( allocated(s3D%proj_space_corrs)    ) deallocate(s3D%proj_space_corrs)
        if( allocated(s3D%proj_space_refinds)  ) deallocate(s3D%proj_space_refinds)
        if( allocated(s3D%proj_space_inplinds) ) deallocate(s3D%proj_space_inplinds)
        if( allocated(s3D%proj_space_state)    ) deallocate(s3D%proj_space_state)
        if( allocated(s3D%proj_space_proj)     ) deallocate(s3D%proj_space_proj)
        if( allocated(s3D%prev_proj)           ) deallocate(s3D%prev_proj)
        if( allocated(s3D%srch_order)          ) deallocate(s3D%srch_order)
        if( allocated(s3D%state_exists)        ) deallocate(s3D%state_exists)
    end subroutine clean_strategy3D

end module simple_strategy3D_alloc
