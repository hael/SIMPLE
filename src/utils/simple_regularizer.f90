! orientation regularizer, used in refine3D
module simple_regularizer
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,   only: params_glob
use simple_dist_binfile, only: dist_binfile
use simple_ori,          only: ori
use simple_builder,      only: build_glob
implicit none

public :: regularizer
public :: calc_num2sample, calc_numinpl2sample2D, calc_numcls2sample2D, reg_dist_switch
private
#include "simple_local_flags.inc"

type :: regularizer
    real,    allocatable :: dist_loc_tab(:,:,:) !< nspace, iptcl, 2 (1->dist, 2->inpl_ind)
    integer, allocatable :: ptcl_ref_map(:)     !< ptcl -> ref assignment map
    logical, allocatable :: ptcl_avail(:)       !< for communicating restraints and state selections
    integer, allocatable :: subspace_inds(:)    !< indices for multi-neighborhood search
    logical, allocatable :: l_neigh(:,:)        !< logical representation of multi-neighborhood
    logical              :: do_neigh = .false.  !< indicates neighborhood search or not
  contains
    ! CONSTRUCTOR
    procedure :: new
    ! PROCEDURES
    procedure :: fill_tab
    procedure :: tab_normalize
    procedure :: tab_align
    procedure :: write_tab
    procedure :: read_tab
    procedure :: read_tab_from_glob
    procedure :: read_tab_to_glob
    procedure :: write_assignment
    procedure :: read_assignment
    ! DESTRUCTOR
    procedure :: kill
end type regularizer

contains

    ! CONSTRUCTORS

    subroutine new( self, l_neigh )
        use simple_oris, only: oris
        class(regularizer), intent(inout) :: self
        logical,            intent(in)    :: l_neigh
        type(oris) :: eulspace_sub !< discrete projection direction search space, reduced
        type(ori)  :: o
        integer    :: i, iptcl
        call self%kill
        allocate(self%dist_loc_tab(params_glob%nspace,params_glob%fromp:params_glob%top,2), source=0.)
        allocate(self%ptcl_ref_map(params_glob%fromp:params_glob%top))
        allocate(self%ptcl_avail(params_glob%fromp:params_glob%top), source=.true.)
        do iptcl = params_glob%fromp,params_glob%top
            self%ptcl_avail(iptcl) = (build_glob%spproj_field%get_state(iptcl) > 0) ! state selection, 0 means deselected
        enddo
        if( l_neigh )then
            ! construct projection direction subspace
            call eulspace_sub%new(params_glob%nspace_sub, is_ptcl=.false.)
            call build_glob%pgrpsyms%build_refspiral(eulspace_sub)
            ! get indexing in original projection direction space
            allocate(self%subspace_inds(params_glob%nspace_sub), source=0)
            !$omp parallel do default(shared) proc_bind(close) private(i,o)
            do i = 1, params_glob%nspace_sub
                call eulspace_sub%get_ori(i, o)
                self%subspace_inds(i) = build_glob%pgrpsyms%find_closest_proj(build_glob%eulspace, o)
            end do
            !$omp end parallel do
            call eulspace_sub%kill
            allocate(self%l_neigh(params_glob%nspace,params_glob%fromp:params_glob%top), source=.false.)
            self%do_neigh = .true.
        else
            allocate(self%l_neigh(params_glob%nspace,params_glob%fromp:params_glob%top), source=.true.)
            self%do_neigh = .false.
        endif
    end subroutine new

    subroutine fill_tab( self, pftcc, glob_pinds )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(regularizer),      intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: glob_pinds(pftcc%nptcls)
        integer   :: i, j, iref, iptcl, locn(params_glob%npeaks), refs_ns, ipeak, inpl_ns, ithr
        real      :: dists_inpl(pftcc%nrots,nthr_glob), dists_sub(params_glob%nspace_sub,nthr_glob), dists_inpl_sorted(pftcc%nrots,nthr_glob), x
        integer   :: inds_sorted(pftcc%nrots,nthr_glob)
        type(ori) :: o
        if( self%do_neigh )then
            call calc_num2sample(params_glob%nspace, 'dist', refs_ns)
            self%l_neigh = .false.
            !$omp parallel do default(shared) private(i,iptcl,ithr,j,iref,o,locn,ipeak) proc_bind(close) schedule(static)
            do i = 1, pftcc%nptcls
                iptcl = glob_pinds(i)
                if( self%ptcl_avail(iptcl) )then
                    ithr = omp_get_thread_num() + 1
                    do j = 1, params_glob%nspace_sub
                        iref = self%subspace_inds(j)
                        call pftcc%gencorrs(iref, iptcl, dists_inpl(:,ithr))
                        dists_inpl(:,ithr) = reg_dist_switch(dists_inpl(:,ithr))
                        dists_sub(j,ithr)  = minval(dists_inpl(:,ithr)) ! GREEDY HERE ???
                    enddo
                    locn = minnloc(dists_sub(:,ithr), params_glob%npeaks)
                    do ipeak = 1,params_glob%npeaks
                        iref = self%subspace_inds(locn(ipeak))
                        call build_glob%eulspace%get_ori(iref, o)
                        call build_glob%pgrpsyms%nearest_proj_neighbors(build_glob%eulspace, o, refs_ns, self%l_neigh(:,iptcl))
                    end do
                endif
            end do
            !$omp end parallel do
            call seed_rnd
            call calc_num2sample(pftcc%nrots, 'dist_inpl', inpl_ns)
            !$omp parallel do collapse(2) default(shared) private(i,iref,iptcl,ithr) proc_bind(close) schedule(static)
            do iref = 1, params_glob%nspace
                do i = 1, pftcc%nptcls
                    iptcl = glob_pinds(i)
                    if( self%ptcl_avail(iptcl) )then
                        ithr = omp_get_thread_num() + 1
                        if( self%l_neigh(iref,iptcl) )then
                            call pftcc%gencorrs(iref, iptcl, dists_inpl(:,ithr))
                            dists_inpl(:,ithr) = reg_dist_switch(dists_inpl(:,ithr))
                            self%dist_loc_tab(iref,iptcl,2) = inpl_smpl(ithr) ! contained function, below
                            self%dist_loc_tab(iref,iptcl,1) = dists_inpl(int(self%dist_loc_tab(iref,iptcl,2)),ithr)
                        else
                            self%dist_loc_tab(iref,iptcl,2) = 0
                            self%dist_loc_tab(iref,iptcl,1) = huge(x)
                        endif
                    endif
                enddo
            enddo
            !$omp end parallel do
        else
            self%l_neigh = .true.
            call seed_rnd
            call calc_num2sample(pftcc%nrots, 'dist_inpl', inpl_ns)
            !$omp parallel do collapse(2) default(shared) private(i,ithr,iref,iptcl) proc_bind(close) schedule(static)
            do iref = 1, params_glob%nspace
                do i = 1, pftcc%nptcls
                    iptcl = glob_pinds(i)
                    if( self%ptcl_avail(iptcl) )then
                        ithr = omp_get_thread_num() + 1
                        call pftcc%gencorrs(iref, iptcl, dists_inpl(:,ithr))
                        dists_inpl(:,ithr) = reg_dist_switch(dists_inpl(:,ithr))
                        self%dist_loc_tab(iref,iptcl,2) = inpl_smpl(ithr) ! contained function, below
                        self%dist_loc_tab(iref,iptcl,1) = dists_inpl(int(self%dist_loc_tab(iref,iptcl,2)),ithr)
                    endif
                enddo
            enddo
            !$omp end parallel do
        endif

    contains

        ! inpl greedy sampling based on unnormalized corr values of all inpl rotations
        function inpl_smpl( ithr ) result( which )
            integer, intent(in) :: ithr
            integer :: j, which
            real    :: rnd, bound, sum_dist
            dists_inpl_sorted(:,ithr) = dists_inpl(:,ithr)
            inds_sorted(:,ithr)  = (/(j,j=1,pftcc%nrots)/)
            call hpsort(dists_inpl_sorted(:,ithr), inds_sorted(:,ithr))
            rnd      = ran3()
            sum_dist = sum(dists_inpl_sorted(1:inpl_ns,ithr))
            if( sum_dist < TINY )then
                ! uniform sampling
                which = 1 + floor(real(inpl_ns) * rnd)
            else
                ! normalizing within the hard-limit
                dists_inpl_sorted(1:inpl_ns,ithr) = dists_inpl_sorted(1:inpl_ns,ithr) / sum_dist
                bound = 0.
                do which = 1,inpl_ns
                    bound = bound + dists_inpl_sorted(which,ithr)
                    if( rnd >= bound )exit
                enddo
                which = min(which,inpl_ns)
            endif
            which = inds_sorted(which,ithr)
        end function inpl_smpl

    end subroutine fill_tab

    ! reference normalization (same energy) of the global dist value table
    ! [0,1] normalization of the whole table
    subroutine tab_normalize( self )
        class(regularizer), intent(inout) :: self
        integer :: iptcl
        real    :: sum_dist_all, min_dist, max_dist
        ! normalize so prob of each ptcl is between [0,1] for all refs
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl,sum_dist_all)
        do iptcl = params_glob%fromp, params_glob%top
            if( self%ptcl_avail(iptcl) )then
                self%l_neigh(:,iptcl) = (self%dist_loc_tab(:,iptcl,2) > TINY)    ! communicating partition l_neigh to global l_neigh
                sum_dist_all = sum(self%dist_loc_tab(:,iptcl,1), mask=self%l_neigh(:,iptcl))
                if( sum_dist_all < TINY )then
                    self%dist_loc_tab(:,iptcl,1) = 0.
                else
                    self%dist_loc_tab(:,iptcl,1) = self%dist_loc_tab(:,iptcl,1) / sum_dist_all
                endif
            endif
        enddo
        !$omp end parallel do
        ! min/max normalization to obtain values between 0 and 1
        max_dist = 0.
        min_dist = huge(min_dist)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl)&
        !$omp reduction(min:min_dist) reduction(max:max_dist)
        do iptcl = params_glob%fromp,params_glob%top
            if( self%ptcl_avail(iptcl) )then
                max_dist = max(max_dist, maxval(self%dist_loc_tab(:,iptcl,1), dim=1, mask=self%l_neigh(:,iptcl)))
                min_dist = min(min_dist, minval(self%dist_loc_tab(:,iptcl,1), dim=1, mask=self%l_neigh(:,iptcl)))
            endif
        enddo
        !$omp end parallel do
        if( (max_dist - min_dist) < TINY )then
            self%dist_loc_tab(:,:,1) = 0.
        else
            self%dist_loc_tab(:,:,1) = (self%dist_loc_tab(:,:,1) - min_dist) / (max_dist - min_dist)
        endif
        ! make sure that particles outside the multi-neighborhood are deselected (put in the bottom)
        where( .not. self%l_neigh ) self%dist_loc_tab(:,:,1) = 1.
        ! to place unselected particles in the bottom (see below)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl)
        do iptcl = params_glob%fromp,params_glob%top
            if( .not. self%ptcl_avail(iptcl) )then
                self%dist_loc_tab(:,iptcl,1) = 1.
            endif
        enddo
        !$omp end parallel do
    end subroutine tab_normalize

    ! ptcl -> ref assignment using the global normalized dist value table
    subroutine tab_align( self )
        class(regularizer), intent(inout) :: self
        integer :: iref, iptcl, assigned_iref, assigned_ptcl, refs_ns, ref_dist_inds(params_glob%nspace),&
                   &stab_inds(params_glob%fromp:params_glob%top, params_glob%nspace), inds_sorted(params_glob%nspace)
        real    :: sorted_tab(params_glob%fromp:params_glob%top, params_glob%nspace),&
                   &ref_dist(params_glob%nspace), dists_sorted(params_glob%nspace)
        logical :: ptcl_avail(params_glob%fromp:params_glob%top)
        self%ptcl_ref_map = 1
        ! sorting each columns
        call calc_num2sample(params_glob%nspace, 'dist', refs_ns)
        sorted_tab = transpose(self%dist_loc_tab(:,:,1))
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,iptcl)
        do iref = 1, params_glob%nspace
            stab_inds(:,iref) = (/(iptcl,iptcl=params_glob%fromp,params_glob%top)/)
            call hpsort(sorted_tab(:,iref), stab_inds(:,iref))
        enddo
        !$omp end parallel do
        ! first row is the current best ref distribution
        ref_dist_inds = params_glob%fromp
        ref_dist      = sorted_tab(params_glob%fromp,:)
        ptcl_avail    = self%ptcl_avail
        do while( any(ptcl_avail) )
            ! sampling the ref distribution to choose next iref to assign corresponding iptcl to
            assigned_iref = ref_smpl() ! contained function, below
            assigned_ptcl = stab_inds(ref_dist_inds(assigned_iref), assigned_iref)
            if( .not. self%l_neigh(assigned_iref, assigned_ptcl) ) continue       ! not in the neighborhood
            ptcl_avail(assigned_ptcl)        = .false.
            self%ptcl_ref_map(assigned_ptcl) = assigned_iref
            ! update the ref_dist and ref_dist_inds
            do iref = 1, params_glob%nspace
                do while( ref_dist_inds(iref) < params_glob%top .and. .not.(ptcl_avail(stab_inds(ref_dist_inds(iref), iref))) )
                    ref_dist_inds(iref) = ref_dist_inds(iref) + 1
                    ref_dist(iref)      = sorted_tab(ref_dist_inds(iref), iref)
                enddo
            enddo
        enddo

    contains

        ! ref greedy sampling based on the current reference distribution
        function ref_smpl( ) result( which )
            integer :: i, which
            real    :: rnd, bound, sum_dists_sorted
            rnd          = ran3()
            dists_sorted = ref_dist
            inds_sorted  = (/(i,i=1,params_glob%nspace)/)
            call hpsort(dists_sorted, inds_sorted)
            sum_dists_sorted = sum(dists_sorted(1:refs_ns))
            if( sum_dists_sorted < TINY )then
                ! uniform sampling
                which = 1 + floor(real(refs_ns) * rnd)
            else
                ! normalizing within the hard-limit
                dists_sorted(1:refs_ns) = dists_sorted(1:refs_ns) / sum_dists_sorted
                bound = 0.
                do which = 1,refs_ns
                    bound = bound + dists_sorted(which)
                    if( rnd >= bound )exit
                enddo
                which = min(which, refs_ns)
            endif
            which = inds_sorted(which)
        end function ref_smpl
        
    end subroutine tab_align

    ! FILE IO

    ! write the partition-wise (or global) dist value table to a binary file
    subroutine write_tab( self, binfname )
        class(regularizer), intent(in) :: self
        character(len=*),   intent(in) :: binfname
        type(dist_binfile) :: binfile
        call binfile%new(binfname, params_glob%fromp, params_glob%top, params_glob%nspace)
        call binfile%write(self%dist_loc_tab)
        call binfile%kill
    end subroutine write_tab

    ! read the partition-wise (or global) dist value binary file to partition-wise (or global) reg object's dist value table
    subroutine read_tab( self, binfname )
        class(regularizer), intent(inout) :: self
        character(len=*),   intent(in)    :: binfname
        type(dist_binfile) :: binfile
        if( file_exists(binfname) )then
            call binfile%new_from_file(binfname)
        else
            call binfile%new(binfname, params_glob%fromp, params_glob%top, params_glob%nspace)
        endif
        call binfile%read(self%dist_loc_tab)
        call binfile%kill
    end subroutine read_tab

    ! read the global dist value binary file to partition-wise reg object's dist value table
    ! [fromp, top]: partition particle index range
    subroutine read_tab_from_glob( self, binfname, fromp, top )
        class(regularizer), intent(inout) :: self
        character(len=*),   intent(in)    :: binfname
        integer,            intent(in)    :: fromp, top
        type(dist_binfile) :: binfile
        integer            :: iptcl, iref
        if( file_exists(binfname) )then
            call binfile%new_from_file(binfname)
        else
            call binfile%new(binfname, params_glob%fromp, params_glob%top, params_glob%nspace)
        endif
        call binfile%read_from_glob(fromp, top, self%dist_loc_tab)
        call binfile%kill
    end subroutine read_tab_from_glob

    ! read the partition-wise dist value binary file to global reg object's dist value table
    ! [fromp, top]: global partition particle index range
    subroutine read_tab_to_glob( self, binfname, fromp, top )
        class(regularizer), intent(inout) :: self
        character(len=*),   intent(in)    :: binfname
        integer,            intent(in)    :: fromp, top
        type(dist_binfile) :: binfile
        if( file_exists(binfname) )then
            call binfile%new_from_file(binfname)
        else
            THROW_HARD( 'corr/rot files of partitions should be ready! ' )
        endif
        call binfile%read_to_glob(fromp, top, self%dist_loc_tab)
        call binfile%kill
    end subroutine read_tab_to_glob

    ! write a global assignment map to binary file
    subroutine write_assignment( self, binfname )
        class(regularizer), intent(in) :: self
        character(len=*),   intent(in) :: binfname
        integer(kind=8) :: file_header(3)
        integer :: funit, io_stat, addr, iptcl, datasz
        datasz      = sizeof(iptcl)
        file_header = [params_glob%fromp, params_glob%top,params_glob%nspace]
        call fopen(funit,trim(binfname),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        ! write header
        write(unit=funit,pos=1) file_header
        ! write assignment
        addr = sizeof(file_header) + 1
        do iptcl = params_glob%fromp, params_glob%top
            write(funit, pos=addr) iptcl
            addr = addr + datasz
            write(funit, pos=addr) self%ptcl_ref_map(iptcl)
            addr = addr + datasz
        end do
        call fclose(funit)
    end subroutine write_assignment

    ! read from the global assignment map to local partition for shift search and further refinement
    subroutine read_assignment( self, binfname )
        class(regularizer), intent(inout) :: self
        character(len=*),   intent(in)    :: binfname
        integer(kind=8) :: file_header(3)
        integer :: funit, io_stat, addr, iptcl, datasz, fromp, top, iglob
        datasz = sizeof(iptcl)
        if( .not. file_exists(trim(binfname)) )then
            THROW_HARD('file '//trim(binfname)//' does not exists!')
        else
            call fopen(funit,trim(binfname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
        end if
        call fileiochk('dist_binfile; read_header; file: '//trim(binfname), io_stat)
        read(unit=funit,pos=1) file_header
        fromp = file_header(1)
        top   = file_header(2)
        ! write assignment
        addr = sizeof(file_header) + 1
        do iglob = fromp, top
            read(unit=funit,pos=addr) iptcl
            addr = addr + datasz
            if( iptcl >= params_glob%fromp .and. iptcl <= params_glob%top ) read(unit=funit,pos=addr) self%ptcl_ref_map(iptcl)
            addr = addr + datasz
        end do
        call fclose(funit)
    end subroutine read_assignment

    ! DESTRUCTOR

    subroutine kill( self )
        class(regularizer), intent(inout) :: self
        if( allocated(self%dist_loc_tab)  ) deallocate(self%dist_loc_tab)
        if( allocated(self%ptcl_ref_map)  ) deallocate(self%ptcl_ref_map)
        if( allocated(self%ptcl_avail)    ) deallocate(self%ptcl_avail)
        if( allocated(self%subspace_inds) ) deallocate(self%subspace_inds)
        if( allocated(self%l_neigh)       ) deallocate(self%l_neigh)
    end subroutine kill

    ! PUBLIC UTILITITES

    subroutine calc_num2sample( num_all, field_str, num_smpl)
        use simple_builder, only: build_glob
        integer,          intent(in)  :: num_all
        character(len=*), intent(in)  :: field_str
        integer,          intent(out) :: num_smpl
        real,    allocatable :: vals(:)
        logical, allocatable :: ptcl_mask(:)
        real    :: athres, dist_thres
        integer :: n
        ptcl_mask  = nint(build_glob%spproj_field%get_all('state')) == 1
        n          = count(ptcl_mask)
        vals       = build_glob%spproj_field%get_all(trim(field_str))
        dist_thres = sum(vals, mask=ptcl_mask) / real(n)
        athres     = params_glob%reg_athres
        if( dist_thres > TINY ) athres = min(athres, dist_thres)
        num_smpl   = min(num_all,max(1,int(athres * real(num_all) / 180.)))
    end subroutine calc_num2sample

    subroutine calc_numcls2sample2D( neigh_frac, ncls, ncls2smpl )
        real,    intent(in)  :: neigh_frac
        integer, intent(in)  :: ncls
        integer, intent(out) :: ncls2smpl
        real    :: athres, dist_thres
        dist_thres = neigh_frac * 180.
        athres     = params_glob%reg_athres
        if( dist_thres > TINY ) athres = min(athres, dist_thres)
        ncls2smpl     = min(ncls, max(1,nint(athres * real(ncls) / 180.)))
    end subroutine calc_numcls2sample2D

    subroutine calc_numinpl2sample2D( num_all, num_smpl )
        use simple_builder, only: build_glob
        integer, intent(in)  :: num_all
        integer, intent(out) :: num_smpl
        real,    allocatable :: vals(:)
        logical, allocatable :: ptcl_mask(:)
        real    :: athres, dist_thres
        integer :: n
        ptcl_mask  = nint(build_glob%spproj_field%get_all('state')) == 1
        ptcl_mask  = ptcl_mask .and. (nint(build_glob%spproj_field%get_all('mi_class')) == 1)
        n          = count(ptcl_mask)
        athres     = params_glob%reg_athres
        if( n > 0 )then
            vals       = build_glob%spproj_field%get_all(trim('dist_inpl'))
            dist_thres = sum(vals, mask=ptcl_mask) / real(n)
            if( dist_thres > TINY ) athres = min(athres, dist_thres)
        endif
        num_smpl = min(num_all,max(1,int(athres * real(num_all) / 180.)))
    end subroutine calc_numinpl2sample2D

    ! switch corr in [0,1] to [0, infinity) to do greedy_sampling
    function reg_dist_switch( corr ) result(dist)
        real, intent(in) :: corr(:)
        real :: dist(size(corr))
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                dist = corr
                where(dist < 0.) dist = 0.
                dist = 1. - dist
            case(OBJFUN_EUCLID)
                where( corr < TINY )
                    dist = huge(dist)
                elsewhere
                    dist = - log(corr)
                endwhere
        end select
    end function reg_dist_switch

end module simple_regularizer
