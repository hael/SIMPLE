! array allocation for concrete strategy3D extensions to improve caching and reduce alloc overheads
module simple_strategy3D_alloc
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_builder,    only: build_glob
use simple_oris,       only: oris
implicit none

public :: s3D, clean_strategy3D, prep_strategy3D, prep_strategy3D_thread
private

integer, parameter :: NTRS_PEAKS = 81

type inpl_peaks
    integer              :: n    = 0        ! # of peaks
    integer              :: ntrs = 0        ! # of shift peaks per psi
    real,    allocatable :: ccs(:,:)        ! score
    real,    allocatable :: shifts(:,:,:)   ! shifts
    integer, allocatable :: eulinds(:,:,:)  ! indices to euler set: projection direction & in-plane
    contains
    procedure :: allocate
    procedure :: deallocate
end type inpl_peaks

type strategy3D_alloc
    ! per-ptcl/ref allocation
    type(oris),       allocatable :: o_peaks(:)                           !< solution objects
    integer,          allocatable :: proj_space_state(:)                  !< states
    integer,          allocatable :: proj_space_proj(:)                   !< projection directions (1 state assumed)
    logical,          allocatable :: state_exists(:)                      !< indicates state existence
    ! per thread allocation
    type(inpl_peaks),    public :: inplpeaks
    type(ran_tabu), allocatable :: rts(:)                                 !< stochastic search order generators
    real,           allocatable :: proj_space_euls(:,:,:,:)               !< euler angles
    real,           allocatable :: proj_space_shift(:,:,:,:)              !< shift vectors
    real,           allocatable :: proj_space_corrs(:,:,:)                !< reference vs. particle correlations
    logical,        allocatable :: proj_space_corrs_srchd(:,:)            !< has correlation already been searched
    logical,        allocatable :: proj_space_corrs_calcd(:,:)            !< has correlation already been calculated
    integer,        allocatable :: proj_space_inplinds(:,:,:)             !< in-plane indices
    integer,        allocatable :: proj_space_refinds_sorted(:,:)         !< reference indices for shift search
    integer,        allocatable :: proj_space_inplinds_sorted(:,:)        !< in-plane indices for shift search
    integer,        allocatable :: proj_space_refinds_sorted_highest(:,:) !< reference indices for shift search (considering only highest scoring inpl)
    integer,        allocatable :: srch_order(:,:)                        !< stochastic search index
end type strategy3D_alloc

type(strategy3D_alloc) :: s3D                           ! singleton
real,      allocatable :: master_proj_space_euls(:,:,:) !< references euler angles
logical                :: srch_order_allocated = .false.

contains

    subroutine prep_strategy3D( ptcl_mask, npeaks )
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
        allocate(master_proj_space_euls(nrefs,params_glob%ninplpeaks,3), s3D%proj_space_euls(nthr_glob,nrefs,params_glob%ninplpeaks,3),&
            &s3D%proj_space_shift(nthr_glob,nrefs,params_glob%ninplpeaks,2),s3D%proj_space_state(nrefs),&
            &s3D%proj_space_corrs(nthr_glob,nrefs,params_glob%ninplpeaks),&
            &s3D%proj_space_refinds_sorted(nthr_glob,nrefs*params_glob%ninplpeaks),&
            &s3D%proj_space_refinds_sorted_highest(nthr_glob,nrefs),&
            &s3D%proj_space_inplinds_sorted(nthr_glob,nrefs*params_glob%ninplpeaks),&
            &s3D%proj_space_corrs_srchd(nthr_glob,nrefs), s3D%proj_space_corrs_calcd(nthr_glob,nrefs),&
            &s3D%proj_space_inplinds(nthr_glob,nrefs,params_glob%ninplpeaks),&
            &s3D%proj_space_proj(nrefs), stat=alloc_stat )
        if( trim(params_glob%trspeaks).eq.'yes' )then
            call s3D%inplpeaks%allocate(params_glob%npeaks, params_glob%ninplpeaks, NTRS_PEAKS, nthr_glob)
        endif
        if(alloc_stat/=0)call allocchk("strategy3D_alloc failed")
        ! states existence
        if( .not.build_glob%spproj%is_virgin_field(params_glob%oritype) )then
            if( str_has_substr(params_glob%refine,'greedy') )then
                allocate(s3D%state_exists(params_glob%nstates), source=.true.)
                ! allocate(s3D%state_exists(params_glob%nstates), source=.true.)
                ! do istate = 1,params_glob%nstates
                !     s3D%state_exists(istate) = file_exists(params_glob%vols(istate))
                ! enddo
            else
                s3D%state_exists = build_glob%spproj_field%states_exist(params_glob%nstates)
            endif
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
                do inpl=1,params_glob%ninplpeaks
                    master_proj_space_euls(cnt,inpl,:) = eul
                end do
            enddo
        enddo
        s3D%proj_space_shift                  = 0.
        s3D%proj_space_corrs                  = -HUGE(areal)
        s3D%proj_space_refinds_sorted         = 0
        s3D%proj_space_inplinds_sorted        = 0
        s3D%proj_space_refinds_sorted_highest = 0
        s3D%proj_space_euls                   = 0.
        s3D%proj_space_corrs_srchd            = .false.
        s3D%proj_space_corrs_calcd            = .false.
        ! search orders allocation
        select case( trim(params_glob%refine) )
            case( 'cluster','clustersym','clustersoft')
                srch_order_allocated = .false.
            case DEFAULT
                if( str_has_substr(params_glob%refine, 'neigh') )then
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
                call s3D%o_peaks(iptcl)%new(npeaks)
                ! transfer CTF params
                if( params_glob%ctf.ne.'no' )then
                    call s3D%o_peaks(iptcl)%set_all2single('iptcl',  real(iptcl))
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
        s3D%proj_space_inplinds(ithr,:,:)               = 0
        s3D%proj_space_refinds_sorted(ithr,:)           = 0
        s3D%proj_space_inplinds_sorted(ithr,:)          = 0
        s3D%proj_space_refinds_sorted_highest(ithr,:)   = 0
        if(srch_order_allocated) s3D%srch_order(ithr,:) = 0
    end subroutine prep_strategy3D_thread

    subroutine clean_strategy3D
        integer :: ithr
        if( allocated(master_proj_space_euls)               ) deallocate(master_proj_space_euls)
        if( allocated(s3D%state_exists)                     ) deallocate(s3D%state_exists)
        if( allocated(s3D%proj_space_state)                 ) deallocate(s3D%proj_space_state)
        if( allocated(s3D%proj_space_proj)                  ) deallocate(s3D%proj_space_proj)
        if( allocated(s3D%srch_order)                       ) deallocate(s3D%srch_order)
        if( allocated(s3D%proj_space_euls)                  ) deallocate(s3D%proj_space_euls)
        if( allocated(s3D%proj_space_shift)                 ) deallocate(s3D%proj_space_shift)
        if( allocated(s3D%proj_space_corrs)                 ) deallocate(s3D%proj_space_corrs)
        if( allocated(s3D%proj_space_corrs_srchd)           ) deallocate(s3D%proj_space_corrs_srchd)
        if( allocated(s3D%proj_space_corrs_calcd)           ) deallocate(s3D%proj_space_corrs_calcd)
        if( allocated(s3D%proj_space_inplinds)              ) deallocate(s3D%proj_space_inplinds)
        if( allocated(s3D%proj_space_refinds_sorted)        ) deallocate(s3D%proj_space_refinds_sorted)
        if( allocated(s3D%proj_space_inplinds_sorted)       ) deallocate(s3D%proj_space_inplinds_sorted)
        if( allocated(s3D%proj_space_refinds_sorted_highest)) deallocate(s3D%proj_space_refinds_sorted_highest)
        call s3D%inplpeaks%deallocate
        if( allocated(s3D%rts) )then
            do ithr=1,nthr_glob
                call s3D%rts(ithr)%kill
            enddo
            deallocate(s3D%rts)
        endif
        srch_order_allocated = .false.
    end subroutine clean_strategy3D

    subroutine allocate( self, nprojdirs, npsis, nshifts, nthr )
        class(inpl_peaks), intent(inout) :: self
        integer,           intent(in)    :: nprojdirs, npsis, nshifts, nthr
        integer :: alloc_stat
        call self%deallocate
        self%ntrs = nshifts
        self%n    = nprojdirs * npsis * self%ntrs
        allocate(self%ccs(self%n,nthr),self%shifts(2,self%n,nthr),self%eulinds(2,self%n,nthr),stat=alloc_stat)
        if(alloc_stat/=0)call allocchk("inpl_peaks%allocate failed")
        self%ccs     = -1.
        self%shifts  = 0.
        self%eulinds = 0
    end subroutine allocate

    subroutine deallocate( self )
        class(inpl_peaks), intent(inout) :: self
        if( allocated(self%ccs) )then
            deallocate(self%ccs,self%eulinds,self%shifts)
            self%n    = 0
            self%ntrs = 0
        endif
    end subroutine deallocate

end module simple_strategy3D_alloc
