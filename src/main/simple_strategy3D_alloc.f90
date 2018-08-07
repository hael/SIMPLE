! array allocation for concrete strategy3D extensions to improve caching and reduce alloc overheads
module simple_strategy3D_alloc
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_builder,    only: build_glob
use simple_oris,       only: oris
implicit none

public :: s3D, clean_strategy3D, prep_strategy3D, prep_strategy3D_thread
private

type strategy3D_alloc
    ! per-ptcl/ref allocation
    type(oris),       allocatable :: o_peaks(:)                      !< solution objects
    integer,          allocatable :: proj_space_state(:)             !< states
    integer,          allocatable :: proj_space_proj(:)              !< projection directions (1 state assumed)
    logical,          allocatable :: state_exists(:)                 !< indicates state existence
    ! per thread allocation
    type(ran_tabu), allocatable :: rts(:)                            !< stochastic serach order generators
    real,           allocatable :: proj_space_euls(:,:,:,:)          !< euler angles
    real,           allocatable :: proj_space_shift(:,:,:,:)         !< shift vectors
    real,           allocatable :: proj_space_corrs(:,:,:)           !< reference vs. particle correlations
    logical,        allocatable :: proj_space_corrs_srchd(:,:)    !< reference vs. particle correlations
    logical,        allocatable :: proj_space_corrs_calcd(:,:)  !< reference vs. particle correlations    
    integer,        allocatable :: proj_space_refinds(:,:)           !< reference indices
    integer,        allocatable :: proj_space_inplinds(:,:,:)        !< in-plane indices
    integer,        allocatable :: srch_order(:,:)                   !< stochastic search index
end type strategy3D_alloc

type(strategy3D_alloc) :: s3D ! singleton
real,      allocatable :: master_proj_space_euls(:,:,:)     !< references euler angles
logical                :: srch_order_allocated = .false.

contains

    subroutine prep_strategy3D(  ptcl_mask, npeaks )
        logical, target, intent(in) :: ptcl_mask(params_glob%fromp:params_glob%top)
        integer,         intent(in) :: npeaks
        real    :: eul(3)
        integer :: i, istate, iproj, iptcl, inpl, ithr
        integer :: nnnrefs, cnt, nrefs
        real    :: areal
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
        ! shared-memory arrays
        allocate(master_proj_space_euls(nrefs,MAXNINPLPEAKS,3), s3D%proj_space_euls(nthr_glob,nrefs,MAXNINPLPEAKS,3),&
            &s3D%proj_space_shift(nthr_glob,nrefs,MAXNINPLPEAKS,2),s3D%proj_space_state(nrefs),&
            &s3D%proj_space_corrs(nthr_glob,nrefs,MAXNINPLPEAKS), s3D%proj_space_refinds(nthr_glob,nrefs),&
            &s3D%proj_space_corrs_srchd(nthr_glob,nrefs), s3D%proj_space_corrs_calcd(nthr_glob,nrefs),&
            &s3D%proj_space_inplinds(nthr_glob,nrefs,MAXNINPLPEAKS),&
            &s3D%proj_space_proj(nrefs), stat=alloc_stat )
        if(alloc_stat/=0)call allocchk("strategy3D_alloc failed")
        ! states existence
        if( .not.build_glob%spproj%is_virgin_field(params_glob%oritype) )then
            s3D%state_exists = build_glob%spproj_field%states_exist(params_glob%nstates)
        else
            allocate(s3D%state_exists(params_glob%nstates), source=.true.)
        endif
        ! reference projection directions state
        cnt = 0
        do istate=1,params_glob%nstates
            do iproj=1,params_glob%nspace
                cnt = cnt + 1
                s3D%proj_space_state(cnt) = istate
                s3D%proj_space_proj(cnt)  = iproj
                eul = build_glob%eulspace%get_euler(iproj)
                do inpl=1,MAXNINPLPEAKS
                    master_proj_space_euls(cnt,inpl,:) = eul
                end do
            enddo
        enddo
        s3D%proj_space_shift       = 0.
        s3D%proj_space_corrs       = -HUGE(areal)
        s3D%proj_space_refinds     = 0
        s3D%proj_space_euls        = 0.
        s3D%proj_space_corrs_srchd = .false.
        s3D%proj_space_corrs_calcd = .false.
        ! search orders allocation
        select case( trim(params_glob%refine) )
            case( 'cluster','clustersym')
                srch_order_allocated = .false.
            case DEFAULT
                if( params_glob%neigh.eq.'yes' )then
                    nnnrefs =  params_glob%nnn * params_glob%nstates
                    allocate(s3D%srch_order(nthr_glob, nnnrefs), s3D%rts(nthr_glob))
                    do ithr=1,nthr_glob
                        s3D%rts(ithr) = ran_tabu(nnnrefs)
                    end do
                else
                    allocate(s3D%srch_order(nthr_glob,nrefs), s3D%rts(nthr_glob))
                    do ithr=1,nthr_glob
                        s3D%rts(ithr) = ran_tabu(nrefs)
                    end do
                endif
                srch_order_allocated = .true.
                s3D%srch_order = 0
        end select
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
    end subroutine prep_strategy3D

    ! init thread specific search arrays
    subroutine prep_strategy3D_thread( ithr )
        integer, intent(in)    :: ithr
        real(sp)               :: areal
        s3D%proj_space_euls(ithr,:,:,:)                 = master_proj_space_euls
        s3D%proj_space_corrs(ithr,:,:)                  = -HUGE(areal)
        s3D%proj_space_corrs_srchd(ithr,:)              = .false.
        s3D%proj_space_corrs_calcd(ithr,:)              = .false.
        s3D%proj_space_shift(ithr,:,:,:)                = 0.
        s3D%proj_space_refinds(ithr,:)                  = 0
        s3D%proj_space_inplinds(ithr,:,:)               = 0
        if(srch_order_allocated) s3D%srch_order(ithr,:) = 0
    end subroutine prep_strategy3D_thread

    subroutine clean_strategy3D
        integer :: ithr
        if( allocated(master_proj_space_euls)      ) deallocate(master_proj_space_euls)
        if( allocated(s3D%state_exists)            ) deallocate(s3D%state_exists)
        if( allocated(s3D%proj_space_state)        ) deallocate(s3D%proj_space_state)
        if( allocated(s3D%proj_space_proj)         ) deallocate(s3D%proj_space_proj)
        if( allocated(s3D%srch_order)              ) deallocate(s3D%srch_order)
        if( allocated(s3D%proj_space_euls)         ) deallocate(s3D%proj_space_euls)
        if( allocated(s3D%proj_space_shift)        ) deallocate(s3D%proj_space_shift)
        if( allocated(s3D%proj_space_corrs)        ) deallocate(s3D%proj_space_corrs)
        if( allocated(s3D%proj_space_corrs_srchd)  ) deallocate(s3D%proj_space_corrs_srchd)
        if( allocated(s3D%proj_space_corrs_calcd)  ) deallocate(s3D%proj_space_corrs_calcd)
        if( allocated(s3D%proj_space_refinds)      ) deallocate(s3D%proj_space_refinds)
        if( allocated(s3D%proj_space_inplinds)     ) deallocate(s3D%proj_space_inplinds)
        if( allocated(s3D%rts) )then
            do ithr=1,nthr_glob
                call s3D%rts(ithr)%kill
            enddo
            deallocate(s3D%rts)
        endif
        srch_order_allocated = .false.
    end subroutine clean_strategy3D

end module simple_strategy3D_alloc
