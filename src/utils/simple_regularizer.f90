! orientation regularizer, used in refine3D
module simple_regularizer
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,   only: params_glob
use simple_corr_binfile, only: corr_binfile
implicit none

public :: regularizer, calc_num2sample
private
#include "simple_local_flags.inc"

type :: regularizer
    real,    allocatable :: corr_loc_tab(:,:,:) !< nspace, iptcl, 2 (1->corr, 2->inpl_ind)
    integer, allocatable :: ptcl_ref_map(:)     !< ptcl -> ref assignment map
    logical, allocatable :: ptcl_avail(:)       !< for communicating restraints and state selections
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

    subroutine new( self )
        use simple_builder, only: build_glob
        class(regularizer), target, intent(inout) :: self
        integer :: iptcl
        allocate(self%corr_loc_tab(params_glob%nspace,params_glob%fromp:params_glob%top,2), source=0.)
        allocate(self%ptcl_ref_map(params_glob%fromp:params_glob%top))
        allocate(self%ptcl_avail(params_glob%fromp:params_glob%top), source=.true.)
        do iptcl = params_glob%fromp,params_glob%top
            self%ptcl_avail(iptcl) = (build_glob%spproj_field%get_state(iptcl) > 0) ! state selection, 0 means deselected
        enddo
    end subroutine new

    ! establish the 2D ref/ptcl table of cost values from gencorrs
    ! (partition-wise table if nparts > 1)
    subroutine fill_tab( self, pftcc, glob_pinds )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(regularizer),      intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: glob_pinds(pftcc%nptcls)
        integer :: i, iref, iptcl, inpl_ns, ithr
        real    :: corrs(pftcc%nrots,params_glob%nthr), corrs_sorted(pftcc%nrots,params_glob%nthr) !< inpl sampling corr
        integer :: inds_sorted(pftcc%nrots,params_glob%nthr)                                       !< inpl sampling indices
        call seed_rnd
        call calc_num2sample(pftcc%nrots, 'dist_inpl', inpl_ns)
        !$omp parallel do collapse(2) default(shared) private(i,iref,iptcl,ithr) proc_bind(close) schedule(static)
        do iref = 1, params_glob%nspace
            do i = 1, pftcc%nptcls
                iptcl = glob_pinds(i)
                if( self%ptcl_avail(iptcl) )then
                    ithr = omp_get_thread_num() + 1
                    call pftcc%gencorrs_prob( iptcl, iref, corrs(:,ithr) )
                    self%corr_loc_tab(iref,iptcl,2) = inpl_smpl(ithr) ! contained function, below
                    self%corr_loc_tab(iref,iptcl,1) = corrs(int(self%corr_loc_tab(iref,iptcl,2)),ithr)
                endif
            enddo
        enddo
        !$omp end parallel do

    contains

        ! inpl greedy sampling based on unnormalized corr values of all inpl rotations
        function inpl_smpl( thread ) result( which )
            integer, intent(in) :: thread
            integer :: j, which
            real    :: rnd, bound, sum_corr
            corrs_sorted(:,thread) = corrs(:,thread)
            inds_sorted(:,thread)  = (/(j,j=1,pftcc%nrots)/)
            call hpsort(corrs_sorted(:,thread), inds_sorted(:,thread) )
            rnd      = ran3()
            sum_corr = sum(corrs_sorted(1:inpl_ns,thread))
            if( sum_corr < TINY )then
                ! uniform sampling
                which = 1 + floor(real(inpl_ns) * rnd)
            else
                ! normalizing within the hard-limit
                corrs_sorted(1:inpl_ns,thread) = corrs_sorted(1:inpl_ns,thread) / sum_corr
                bound = 0.
                do which = 1,inpl_ns
                    bound = bound + corrs_sorted(which, thread)
                    if( rnd >= bound )exit
                enddo
                which = min(which,inpl_ns)
            endif
            which = inds_sorted(which, thread)
        end function inpl_smpl
        
    end subroutine fill_tab

    ! reference normalization (same energy) of the global cost value table
    ! [0,1] normalization of the whole table
    subroutine tab_normalize( self )
        class(regularizer), intent(inout) :: self
        integer :: iptcl
        real    :: sum_corr_all, min_corr, max_corr
        ! normalize so prob of each ptcl is between [0,1] for all refs
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl,sum_corr_all)
        do iptcl = params_glob%fromp, params_glob%top
            if( self%ptcl_avail(iptcl) )then
                sum_corr_all = sum(self%corr_loc_tab(:,iptcl,1))
                if( sum_corr_all < TINY )then
                    self%corr_loc_tab(:,iptcl,1) = 0.
                else
                    self%corr_loc_tab(:,iptcl,1) = self%corr_loc_tab(:,iptcl,1) / sum_corr_all
                endif
            endif
        enddo
        !$omp end parallel do
        ! min/max normalization to obtain values between 0 and 1
        max_corr = 0.
        min_corr = huge(min_corr)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl)&
        !$omp reduction(min:min_corr) reduction(max:max_corr)
        do iptcl = params_glob%fromp,params_glob%top
            if( self%ptcl_avail(iptcl) )then
                max_corr = max(max_corr, maxval(self%corr_loc_tab(:,iptcl,1), dim=1))
                min_corr = min(min_corr, minval(self%corr_loc_tab(:,iptcl,1), dim=1))
            endif
        enddo
        !$omp end parallel do
        if( (max_corr - min_corr) < TINY )then
            self%corr_loc_tab(:,:,1) = 0.
        else
            self%corr_loc_tab(:,:,1) = (self%corr_loc_tab(:,:,1) - min_corr) / (max_corr - min_corr)
        endif
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl)
        do iptcl = params_glob%fromp,params_glob%top
            if( .not. self%ptcl_avail(iptcl) )then
                self%corr_loc_tab(:,iptcl,1) = 1. ! to place unselected particles in the bottom (see below)
            endif
        enddo
        !$omp end parallel do
    end subroutine tab_normalize

    ! ptcl -> ref assignment using the global normalized cost value table
    subroutine tab_align( self )
        class(regularizer), intent(inout) :: self
        integer :: iref, iptcl, assigned_iref, assigned_ptcl, refs_ns, ref_dist_inds(params_glob%nspace),&
                   &stab_inds(params_glob%fromp:params_glob%top, params_glob%nspace), inds_sorted(params_glob%nspace)
        real    :: sorted_tab(params_glob%fromp:params_glob%top, params_glob%nspace),&
                   &ref_dist(params_glob%nspace), corrs_sorted(params_glob%nspace)
        logical :: ptcl_avail(params_glob%fromp:params_glob%top)
        self%ptcl_ref_map = 1
        ! sorting each columns
        call calc_num2sample(params_glob%nspace, 'dist', refs_ns)
        sorted_tab = transpose(self%corr_loc_tab(:,:,1))
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
            real    :: rnd, bound, sum_corrs_sorted
            rnd          = ran3()
            corrs_sorted = ref_dist
            inds_sorted  = (/(i,i=1,params_glob%nspace)/)
            call hpsort(corrs_sorted, inds_sorted)
            sum_corrs_sorted = sum(corrs_sorted(1:refs_ns))
            if( sum_corrs_sorted < TINY )then
                ! uniform sampling
                which = 1 + floor(real(refs_ns) * rnd)
            else
                ! normalizing within the hard-limit
                corrs_sorted(1:refs_ns) = corrs_sorted(1:refs_ns) / sum_corrs_sorted
                bound = 0.
                do which = 1,refs_ns
                    bound = bound + corrs_sorted(which)
                    if( rnd >= bound )exit
                enddo
                which = min(which, refs_ns)
            endif
            which = inds_sorted(which)
        end function ref_smpl
        
    end subroutine tab_align

    ! FILE IO

    ! write the partition-wise (or global) cost value table to a binary file
    subroutine write_tab( self, binfname )
        class(regularizer), intent(in) :: self
        character(len=*),   intent(in) :: binfname
        type(corr_binfile) :: binfile
        if( file_exists(binfname) )then
            call binfile%new_from_file(binfname)
        else
            call binfile%new(binfname, params_glob%fromp, params_glob%top, params_glob%nspace)
        endif
        call binfile%write(self%corr_loc_tab)
        call binfile%kill
    end subroutine write_tab

    ! read the partition-wise (or global) cost value binary file to partition-wise (or global) reg object's cost value table
    subroutine read_tab( self, binfname )
        class(regularizer), intent(inout) :: self
        character(len=*),   intent(in)    :: binfname
        type(corr_binfile) :: binfile
        if( file_exists(binfname) )then
            call binfile%new_from_file(binfname)
        else
            call binfile%new(binfname, params_glob%fromp, params_glob%top, params_glob%nspace)
        endif
        call binfile%read(self%corr_loc_tab)
        call binfile%kill
    end subroutine read_tab

    ! read the global cost value binary file to partition-wise reg object's cost value table
    ! [fromp, top]: partition particle index range
    subroutine read_tab_from_glob( self, binfname, fromp, top )
        class(regularizer), intent(inout) :: self
        character(len=*),   intent(in)    :: binfname
        integer,            intent(in)    :: fromp, top
        type(corr_binfile) :: binfile
        integer            :: iptcl, iref
        if( file_exists(binfname) )then
            call binfile%new_from_file(binfname)
        else
            call binfile%new(binfname, params_glob%fromp, params_glob%top, params_glob%nspace)
        endif
        call binfile%read_from_glob(fromp, top, self%corr_loc_tab)
        call binfile%kill
    end subroutine read_tab_from_glob

    ! read the partition-wise cost value binary file to global reg object's cost value table
    ! [fromp, top]: global partition particle index range
    subroutine read_tab_to_glob( self, binfname, fromp, top )
        class(regularizer), intent(inout) :: self
        character(len=*),   intent(in)    :: binfname
        integer,            intent(in)    :: fromp, top
        type(corr_binfile) :: binfile
        if( file_exists(binfname) )then
            call binfile%new_from_file(binfname)
        else
            THROW_HARD( 'corr/rot files of partitions should be ready! ' )
        endif
        call binfile%read_to_glob(fromp, top, self%corr_loc_tab)
        call binfile%kill
    end subroutine read_tab_to_glob

    ! write a global assignment map to binary file
    subroutine write_assignment( self, binfname )
        class(regularizer), intent(in) :: self
        character(len=*),   intent(in) :: binfname
        integer :: funit, io_stat, addr, iptcl, datasz, pfromto(2)
        datasz  = sizeof(iptcl)
        pfromto = [params_glob%fromp, params_glob%top]
        call fopen(funit,trim(binfname),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        ! write header
        write(unit=funit,pos=1) pfromto
        ! write assignment
        addr = sizeof(pfromto) + 1
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
        integer :: funit, io_stat, addr, iptcl, datasz, file_header(2), fromp, top, iglob
        datasz = sizeof(iptcl)
        if( .not. file_exists(trim(binfname)) )then
            THROW_HARD('file '//trim(binfname)//' does not exists!')
        else
            call fopen(funit,trim(binfname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
        end if
        call fileiochk('corr_binfile; read_header; file: '//trim(binfname), io_stat)
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
        deallocate(self%corr_loc_tab,self%ptcl_ref_map,self%ptcl_avail)
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

end module simple_regularizer
