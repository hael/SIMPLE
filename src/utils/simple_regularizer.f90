! regularizer of the cluster2D and refine3D
module simple_regularizer
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,   only: params_glob
use simple_corr_binfile, only: corr_binfile
use simple_image
implicit none

public :: regularizer, calc_num2sample
private
#include "simple_local_flags.inc"

type reg_params
    integer :: iptcl        !< iptcl index
    integer :: iref         !< iref index
    integer :: loc          !< inpl index
    real    :: prob         !< probability/corr
    real    :: sh(2)        !< shift
    real    :: w            !< weight
end type reg_params

type :: regularizer
    integer                       :: nrefs
    real,             allocatable :: corr_loc_tab(:,:,:)         !< 2D corr/loc table
    integer,          allocatable :: ptcl_ref_map(:)             !< ptcl -> ref assignment map
    real,             allocatable :: refs_corr(:,:)
    integer,          allocatable :: refs_inds(:,:)
    logical,          allocatable :: ptcl_avail(:)
    type(reg_params), allocatable :: ref_ptcl_tab(:,:)
    contains
    ! CONSTRUCTOR
    procedure :: new
    ! PROCEDURES
    procedure :: fill_tab
    procedure :: tab_normalize
    procedure :: tab_align
    procedure :: normalize_weight
    procedure :: shift_search
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
        integer :: iptcl, iref
        self%nrefs = params_glob%nspace
        allocate(self%corr_loc_tab(self%nrefs,params_glob%fromp:params_glob%top,2), source=0.)
        allocate(self%ref_ptcl_tab(self%nrefs,params_glob%fromp:params_glob%top))
        allocate(self%ptcl_ref_map(params_glob%fromp:params_glob%top))
        allocate(self%ptcl_avail(params_glob%fromp:params_glob%top), source=.true.)
        do iptcl = params_glob%fromp,params_glob%top
            self%ptcl_avail(iptcl) = (build_glob%spproj_field%get_state(iptcl) > 0)
            if( self%ptcl_avail(iptcl) )then
                do iref = 1,self%nrefs
                    self%ref_ptcl_tab(iref,iptcl)%iptcl = iptcl
                    self%ref_ptcl_tab(iref,iptcl)%iref  = iref
                    self%ref_ptcl_tab(iref,iptcl)%loc   = 0
                    self%ref_ptcl_tab(iref,iptcl)%prob  = 0.
                    self%ref_ptcl_tab(iref,iptcl)%sh    = 0.
                    self%ref_ptcl_tab(iref,iptcl)%w     = 1.
                enddo
            endif
        enddo
    end subroutine new

    subroutine fill_tab( self, pftcc, glob_pinds )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(regularizer),      intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: glob_pinds(pftcc%nptcls)
        integer :: i, iref, iptcl, inpl_ns
        real    :: corrs(pftcc%nrots), inpl_corr(pftcc%nrots,params_glob%nthr)   !< inpl sampling corr
        integer :: inpl_inds(pftcc%nrots,params_glob%nthr)                       !< inpl sampling indices
        call seed_rnd
        call calc_num2sample(pftcc%nrots, 'dist_inpl', inpl_ns)
        !$omp parallel do collapse(2) default(shared) private(i,iref,iptcl,corrs) proc_bind(close) schedule(static)
        do iref = 1, self%nrefs
            do i = 1, pftcc%nptcls
                iptcl = glob_pinds(i)
                if( self%ptcl_avail(iptcl) )then
                    call pftcc%gencorrs_prob( iptcl, iref, corrs )
                    self%corr_loc_tab(iref,iptcl,2) = inpl_smpl()
                    self%corr_loc_tab(iref,iptcl,1) = corrs(int(self%corr_loc_tab(iref,iptcl,2)))
                endif
            enddo
        enddo
        !$omp end parallel do
    contains
        ! inpl greedy sampling based on unnormalized pvec
        function inpl_smpl( ) result( which )
            integer :: j, which, ithr
            real    :: rnd, bound, sum_corr
            ithr = omp_get_thread_num() + 1
            inpl_corr(:,ithr) = corrs
            inpl_inds(:,ithr) = (/(j,j=1,pftcc%nrots)/)
            call hpsort(inpl_corr(:,ithr), inpl_inds(:,ithr) )
            rnd      = ran3()
            sum_corr = sum(inpl_corr(1:inpl_ns,ithr))
            if( sum_corr < TINY )then
                ! uniform sampling
                which = 1 + floor(real(inpl_ns) * rnd)
            else
                ! normalizing within the hard-limit
                inpl_corr(1:inpl_ns,ithr) = inpl_corr(1:inpl_ns,ithr) / sum_corr
                bound = 0.
                do which=1,inpl_ns
                    bound = bound + inpl_corr(which, ithr)
                    if( rnd >= bound )exit
                enddo
                which = min(which,inpl_ns)
            endif
            which = inpl_inds(which, ithr)
        end function inpl_smpl
    end subroutine fill_tab

    subroutine tab_normalize( self )
        class(regularizer), intent(inout) :: self
        integer :: iref, iptcl
        real    :: sum_corr_all, min_corr, max_corr
        ! normalize so prob of each ptcl is between [0,1] for all refs
        if( params_glob%l_reg_norm )then
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
        endif
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
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,iptcl)
        do iptcl = params_glob%fromp,params_glob%top
            if( self%ptcl_avail(iptcl) )then
                do iref = 1, self%nrefs
                    self%ref_ptcl_tab(iref,iptcl)%prob =     self%corr_loc_tab(iref,iptcl,1)
                    self%ref_ptcl_tab(iref,iptcl)%loc  = int(self%corr_loc_tab(iref,iptcl,2))
                enddo
            else
                self%corr_loc_tab(:,iptcl,1) = 1.     ! unselected particles are down on the sorted corrs
            endif
        enddo
        !$omp end parallel do
    end subroutine tab_normalize

    subroutine shift_search( self, glob_pinds )
        use simple_pftcc_shsrch_reg, only: pftcc_shsrch_reg
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: glob_pinds(:)
        type(pftcc_shsrch_reg) :: grad_shsrch_obj(params_glob%nthr)
        integer :: iref, iptcl, ithr, irot, i, nptcls
        real    :: lims(2,2), cxy(3)
        nptcls    = size(glob_pinds)
        lims(1,1) = -params_glob%trs
        lims(1,2) =  params_glob%trs
        lims(2,1) = -params_glob%trs
        lims(2,2) =  params_glob%trs
        do ithr = 1, params_glob%nthr
            call grad_shsrch_obj(ithr)%new(lims, opt_angle=.true.)
        enddo
        !$omp parallel do default(shared) private(i,iref,iptcl,irot,ithr,cxy) proc_bind(close) schedule(static)
        do i = 1, nptcls
            iptcl = glob_pinds(i)
            if( self%ptcl_avail(iptcl) )then
                iref  = self%ptcl_ref_map(iptcl)
                ithr  = omp_get_thread_num() + 1
                call grad_shsrch_obj(ithr)%set_indices(iref, iptcl)
                irot = self%ref_ptcl_tab(iref,iptcl)%loc
                cxy  = grad_shsrch_obj(ithr)%minimize(irot)
                if( irot > 0 )then
                    self%ref_ptcl_tab(iref,iptcl)%sh  = cxy(2:3)
                    self%ref_ptcl_tab(iref,iptcl)%loc = irot
                endif
            endif
        enddo
        !$omp end parallel do
    end subroutine shift_search

    subroutine tab_align( self )
        class(regularizer), intent(inout) :: self
        integer :: iref, iptcl, assigned_iref, assigned_ptcl, refs_ns,&
                  &ref_dist_inds(self%nrefs), stab_inds(params_glob%fromp:params_glob%top, self%nrefs)
        real    :: sorted_tab(params_glob%fromp:params_glob%top, self%nrefs), ref_dist(self%nrefs),&
                  &refs_corr(self%nrefs,params_glob%nthr)
        logical :: ptcl_avail(params_glob%fromp:params_glob%top)
        integer :: refs_inds(self%nrefs,params_glob%nthr)
        self%ptcl_ref_map = 1
        ! sorting each columns
        call calc_num2sample(self%nrefs, 'dist', refs_ns)
        sorted_tab = transpose(self%corr_loc_tab(:,:,1))
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,iptcl)
        do iref = 1, self%nrefs
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
            assigned_iref = ref_smpl()
            assigned_ptcl = stab_inds(ref_dist_inds(assigned_iref), assigned_iref)
            ptcl_avail(assigned_ptcl)        = .false.
            self%ptcl_ref_map(assigned_ptcl) = assigned_iref
            ! update the ref_dist and ref_dist_inds
            do iref = 1, self%nrefs
                do while( ref_dist_inds(iref) < params_glob%top .and. .not.(ptcl_avail(stab_inds(ref_dist_inds(iref), iref))) )
                    ref_dist_inds(iref) = ref_dist_inds(iref) + 1
                    ref_dist(iref)      = sorted_tab(ref_dist_inds(iref), iref)
                enddo
            enddo
        enddo
    contains
        function ref_smpl( ) result( which )
            integer :: i, which, ithr
            real    :: rnd, bound, sum_refs_corr
            ithr = omp_get_thread_num() + 1
            rnd  = ran3()
            refs_corr(:,ithr) = ref_dist
            refs_inds(:,ithr) = (/(i,i=1,self%nrefs)/)
            call hpsort(refs_corr(:,ithr), refs_inds(:,ithr) )
            sum_refs_corr = sum(refs_corr(1:refs_ns,ithr))
            if( sum_refs_corr < TINY )then
                ! uniform sampling
                which = 1 + floor(real(refs_ns) * rnd)
            else
                ! normalizing within the hard-limit
                refs_corr(1:refs_ns,ithr) = refs_corr(1:refs_ns,ithr) / sum_refs_corr
                bound = 0.
                do which=1,refs_ns
                    bound = bound + refs_corr(which, ithr)
                    if( rnd >= bound )exit
                enddo
                which = min(which, refs_ns)
            endif
            which = refs_inds(which, ithr)
        end function ref_smpl
    end subroutine tab_align

    subroutine normalize_weight( self )
        class(regularizer), intent(inout) :: self
        integer :: iptcl, iref
        real    :: min_w, max_w
        min_w = huge(min_w)
        max_w = 0.
        do iptcl = params_glob%fromp, params_glob%top
            if( self%ptcl_avail(iptcl) )then
                iref = self%ptcl_ref_map(iptcl)
                self%ref_ptcl_tab(iref,iptcl)%w = 1. - self%ref_ptcl_tab(iref,iptcl)%prob
                if( self%ref_ptcl_tab(iref,iptcl)%w < min_w ) min_w = self%ref_ptcl_tab(iref,iptcl)%w
                if( self%ref_ptcl_tab(iref,iptcl)%w > max_w ) max_w = self%ref_ptcl_tab(iref,iptcl)%w
            endif
        enddo
        do iptcl = params_glob%fromp, params_glob%top
            if( self%ptcl_avail(iptcl) )then
                iref = self%ptcl_ref_map(iptcl)
                self%ref_ptcl_tab(iref,iptcl)%w = (self%ref_ptcl_tab(iref,iptcl)%w - min_w) / (max_w - min_w)
            endif
        enddo
    end subroutine normalize_weight

    ! FILE IO
    subroutine write_tab( self, binfname )
        class(regularizer), intent(in) :: self
        character(len=*),   intent(in) :: binfname
        type(corr_binfile) :: binfile
        if( file_exists(binfname) )then
            call binfile%new_from_file(binfname)
        else
            call binfile%new(binfname, params_glob%fromp, params_glob%top, self%nrefs)
        endif
        call binfile%write(self%corr_loc_tab)
        call binfile%kill
    end subroutine write_tab

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

    subroutine read_tab( self, binfname )
        class(regularizer), intent(inout) :: self
        character(len=*),   intent(in)    :: binfname
        type(corr_binfile) :: binfile
        if( file_exists(binfname) )then
            call binfile%new_from_file(binfname)
        else
            call binfile%new(binfname, params_glob%fromp, params_glob%top, self%nrefs)
        endif
        call binfile%read(self%corr_loc_tab)
        call binfile%kill
    end subroutine read_tab

    subroutine read_tab_from_glob( self, binfname, fromp, top )
        class(regularizer), intent(inout) :: self
        character(len=*),   intent(in)    :: binfname
        integer,            intent(in)    :: fromp, top
        type(corr_binfile) :: binfile
        integer            :: iptcl, iref
        if( file_exists(binfname) )then
            call binfile%new_from_file(binfname)
        else
            call binfile%new(binfname, params_glob%fromp, params_glob%top, self%nrefs)
        endif
        call binfile%read_from_glob(fromp, top, self%corr_loc_tab)
        call binfile%kill
        ! transfer corr_loc_tab into ref_ptcl_tab
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,iptcl)
        do iptcl = params_glob%fromp,params_glob%top
            if( self%ptcl_avail(iptcl) )then
                do iref = 1, self%nrefs
                    self%ref_ptcl_tab(iref,iptcl)%prob =     self%corr_loc_tab(iref,iptcl,1)
                    self%ref_ptcl_tab(iref,iptcl)%loc  = int(self%corr_loc_tab(iref,iptcl,2))
                enddo
            endif
        enddo
        !$omp end parallel do
    end subroutine read_tab_from_glob

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
        deallocate(self%corr_loc_tab,self%ref_ptcl_tab,self%ptcl_ref_map,self%ptcl_avail)
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
