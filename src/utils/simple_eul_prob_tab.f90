! orientation eul_prob_tab, used in refine3D
module simple_eul_prob_tab
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,   only: params_glob
use simple_dist_binfile, only: dist_binfile
use simple_builder,      only: build_glob
implicit none

public :: eul_prob_tab
public :: calc_num2sample, calc_numinpl2sample2D, eulprob_dist_switch, eulprob_corr_switch, shift_sampling
private
#include "simple_local_flags.inc"

type :: eul_prob_tab
    real,    allocatable :: dist_loc_tab(:,:,:) !< nspace, iptcl, 2 (1->dist, 2->inpl_ind)
    real,    allocatable :: ptcl_ref_map(:,:)   !< nptcls, 3 (1->ref,2->dist, 3->inpl_ind)
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
end type eul_prob_tab

contains

    ! CONSTRUCTORS

    subroutine new( self )
        class(eul_prob_tab), intent(inout) :: self
        integer :: iptcl
        call self%kill
        allocate(self%dist_loc_tab(params_glob%nspace,params_glob%fromp:params_glob%top,2), source=0.)
        allocate(self%ptcl_ref_map(params_glob%fromp:params_glob%top,3))
        allocate(self%ptcl_avail(params_glob%fromp:params_glob%top), source=.true.)
        do iptcl = params_glob%fromp,params_glob%top
            self%ptcl_avail(iptcl) = (build_glob%spproj_field%get_state(iptcl) > 0) ! state selection, 0 means deselected
        enddo
    end subroutine new

    ! partition-wise table filling, used only in shared-memory commander 'exec_prob_tab'
    subroutine fill_tab( self, pftcc, glob_pinds )
        use simple_polarft_corrcalc,  only: polarft_corrcalc
        use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
        class(eul_prob_tab),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: glob_pinds(pftcc%nptcls)
        integer,      parameter :: MAXITS = 60
        integer,    allocatable :: locn(:)
        type(pftcc_shsrch_grad) :: grad_shsrch_obj(nthr_glob) !< origin shift search object, L-BFGS with gradient
        integer :: i, j, iref, iptcl, refs_ns, inpl_ns, ithr, irot, inds_sorted(pftcc%nrots,nthr_glob)
        real    :: dists_inpl(pftcc%nrots,nthr_glob), dists_inpl_sorted(pftcc%nrots,nthr_glob)
        real    :: dists_refs(pftcc%nrefs,nthr_glob), lims(2,2), lims_init(2,2), cxy(3)
        call seed_rnd
        call calc_num2sample(pftcc%nrots,        'dist_inpl', inpl_ns)
        call calc_num2sample(params_glob%nspace, 'dist',      refs_ns)
        if( params_glob%l_prob_sh )then
            allocate(locn(refs_ns), source=0)
            ! make shift search objects
            lims(:,1)      = -params_glob%trs
            lims(:,2)      =  params_glob%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new(lims, lims_init=lims_init, shbarrier=params_glob%shbarrier, maxits=MAXITS, opt_angle=.false.)
            end do
            ! fill the table
            !$omp parallel do default(shared) private(i,j,iptcl,ithr,iref,irot,cxy,locn) proc_bind(close) schedule(static)
            do i = 1, pftcc%nptcls
                iptcl = glob_pinds(i)
                if( .not. self%ptcl_avail(iptcl) ) cycle
                ithr = omp_get_thread_num() + 1
                do iref = 1, params_glob%nspace
                    call pftcc%gencorrs(iref, iptcl, dists_inpl(:,ithr))
                    dists_inpl(:,ithr) = eulprob_dist_switch(dists_inpl(:,ithr))
                    irot = inpl_smpl(ithr) ! contained function, below
                    self%dist_loc_tab(iref,iptcl,1) = dists_inpl(irot,ithr)
                    self%dist_loc_tab(iref,iptcl,2) = irot
                    dists_refs(iref,ithr) = dists_inpl(irot,ithr)
                enddo
                locn = minnloc(dists_refs(:,ithr), refs_ns)
                if( params_glob%l_doshift )then
                    do j = 1,refs_ns
                        iref = locn(j)
                        ! BFGS over shifts
                        call grad_shsrch_obj(ithr)%set_indices(iref, iptcl)
                        irot = nint(self%dist_loc_tab(iref,iptcl,2))
                        cxy  = grad_shsrch_obj(ithr)%minimize(irot=irot)
                        if( irot > 0 )then
                            ! no storing of shifts for now, re-search with in-plane jiggle in strategy3D_prob
                            self%dist_loc_tab(iref,iptcl,1) = eulprob_dist_switch(cxy(1))
                        endif
                    end do
                endif
            enddo
            !$omp end parallel do
        else
            ! fill the table
            !$omp parallel do collapse(2) default(shared) private(i,ithr,iref,iptcl,irot) proc_bind(close) schedule(static)
            do iref = 1, params_glob%nspace
                do i = 1, pftcc%nptcls
                    iptcl = glob_pinds(i)
                    if( self%ptcl_avail(iptcl) )then
                        ithr = omp_get_thread_num() + 1
                        call pftcc%gencorrs(iref, iptcl, dists_inpl(:,ithr))
                        dists_inpl(:,ithr) = eulprob_dist_switch(dists_inpl(:,ithr))
                        irot = inpl_smpl(ithr) ! contained function, below
                        self%dist_loc_tab(iref,iptcl,1) = dists_inpl(irot,ithr)
                        self%dist_loc_tab(iref,iptcl,2) = irot
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
            inds_sorted(:,ithr)       = (/(j,j=1,pftcc%nrots)/)
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
    ! used in the global prob_align commander, in 'exec_prob_align'
    subroutine tab_normalize( self )
        class(eul_prob_tab), intent(inout) :: self
        integer :: iptcl
        real    :: sum_dist_all, min_dist, max_dist
        ! normalize so prob of each ptcl is between [0,1] for all refs
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl,sum_dist_all)
        do iptcl = params_glob%fromp, params_glob%top
            if( self%ptcl_avail(iptcl) )then
                sum_dist_all = sum(self%dist_loc_tab(:,iptcl,1))
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
                max_dist = max(max_dist, maxval(self%dist_loc_tab(:,iptcl,1), dim=1))
                min_dist = min(min_dist, minval(self%dist_loc_tab(:,iptcl,1), dim=1))
            endif
        enddo
        !$omp end parallel do
        if( (max_dist - min_dist) < TINY )then
            self%dist_loc_tab(:,:,1) = 0.
        else
            self%dist_loc_tab(:,:,1) = (self%dist_loc_tab(:,:,1) - min_dist) / (max_dist - min_dist)
        endif
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
    ! used in the global prob_align commander, in 'exec_prob_align'
    subroutine tab_align( self )
        class(eul_prob_tab), intent(inout) :: self
        integer :: iref, iptcl, assigned_iref, assigned_ptcl, refs_ns, ref_dist_inds(params_glob%nspace),&
                   &stab_inds(params_glob%fromp:params_glob%top, params_glob%nspace), inds_sorted(params_glob%nspace)
        real    :: sorted_tab(params_glob%fromp:params_glob%top, params_glob%nspace),&
                   &ref_dist(params_glob%nspace), dists_sorted(params_glob%nspace)
        logical :: ptcl_avail(params_glob%fromp:params_glob%top)
        self%ptcl_ref_map = 1.
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
            ptcl_avail(assigned_ptcl)          = .false.
            self%ptcl_ref_map(assigned_ptcl,1) = assigned_iref
            self%ptcl_ref_map(assigned_ptcl,2) = self%dist_loc_tab(assigned_iref,assigned_ptcl,1)
            self%ptcl_ref_map(assigned_ptcl,3) = self%dist_loc_tab(assigned_iref,assigned_ptcl,2)
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
        class(eul_prob_tab), intent(in) :: self
        character(len=*),    intent(in) :: binfname
        type(dist_binfile) :: binfile
        call binfile%new(binfname, params_glob%fromp, params_glob%top, params_glob%nspace)
        call binfile%write(self%dist_loc_tab)
        call binfile%kill
    end subroutine write_tab

    ! read the partition-wise (or global) dist value binary file to partition-wise (or global) reg object's dist value table
    subroutine read_tab( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        character(len=*),    intent(in)    :: binfname
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
        class(eul_prob_tab), intent(inout) :: self
        character(len=*),    intent(in)    :: binfname
        integer,             intent(in)    :: fromp, top
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
        class(eul_prob_tab), intent(inout) :: self
        character(len=*),    intent(in)    :: binfname
        integer,             intent(in)    :: fromp, top
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
        class(eul_prob_tab), intent(in) :: self
        character(len=*),    intent(in) :: binfname
        integer(kind=8) :: file_header(3), addr
        integer :: funit, io_stat, iptcl, datasz_int, datasz_real
        real    :: x
        datasz_int  = sizeof(iptcl)
        datasz_real = sizeof(x)
        file_header = [params_glob%fromp, params_glob%top,params_glob%nspace]
        call fopen(funit,trim(binfname),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        ! write header
        write(unit=funit,pos=1) file_header
        ! write assignment
        addr = sizeof(file_header) + 1
        do iptcl = params_glob%fromp, params_glob%top
            write(funit, pos=addr) iptcl
            addr = addr + datasz_int
            write(funit, pos=addr) self%ptcl_ref_map(iptcl, 1)
            addr = addr + datasz_real
            write(funit, pos=addr) self%ptcl_ref_map(iptcl, 2)
            addr = addr + datasz_real
            write(funit, pos=addr) self%ptcl_ref_map(iptcl, 3)
            addr = addr + datasz_real
        end do
        call fclose(funit)
    end subroutine write_assignment

    ! read from the global assignment map to local partition for shift search and further refinement
    subroutine read_assignment( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        character(len=*),    intent(in)    :: binfname
        integer(kind=8) :: file_header(3), addr
        integer :: funit, io_stat, iptcl, datasz_int, datasz_real, fromp, top, iglob
        real    :: x
        datasz_int  = sizeof(iptcl)
        datasz_real = sizeof(x)
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
            addr = addr + datasz_int
            if( iptcl >= params_glob%fromp .and. iptcl <= params_glob%top ) read(unit=funit,pos=addr) self%ptcl_ref_map(iptcl,1)
            addr = addr + datasz_real
            if( iptcl >= params_glob%fromp .and. iptcl <= params_glob%top ) read(unit=funit,pos=addr) self%ptcl_ref_map(iptcl,2)
            addr = addr + datasz_real
            if( iptcl >= params_glob%fromp .and. iptcl <= params_glob%top ) read(unit=funit,pos=addr) self%ptcl_ref_map(iptcl,3)
            addr = addr + datasz_real
        end do
        call fclose(funit)
    end subroutine read_assignment

    ! DESTRUCTOR

    subroutine kill( self )
        class(eul_prob_tab), intent(inout) :: self
        if( allocated(self%dist_loc_tab)  ) deallocate(self%dist_loc_tab)
        if( allocated(self%ptcl_ref_map)  ) deallocate(self%ptcl_ref_map)
        if( allocated(self%ptcl_avail)    ) deallocate(self%ptcl_avail)
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
        athres     = params_glob%prob_athres
        if( dist_thres > TINY ) athres = min(athres, dist_thres)
        num_smpl   = min(num_all,max(1,int(athres * real(num_all) / 180.)))
    end subroutine calc_num2sample

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
        athres     = params_glob%prob_athres
        if( n > 0 )then
            vals       = build_glob%spproj_field%get_all(trim('dist_inpl'))
            dist_thres = sum(vals, mask=ptcl_mask) / real(n)
            if( dist_thres > TINY ) athres = min(athres, dist_thres)
        endif
        num_smpl = min(num_all,max(1,int(athres * real(num_all) / 180.)))
    end subroutine calc_numinpl2sample2D

    ! switch corr in [0,1] to [0, infinity) to do greedy_sampling
    elemental function eulprob_dist_switch( corr ) result(dist)
        real, intent(in) :: corr
        real :: dist
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                if( corr < 0. )then
                    dist = 0.
                else
                    dist = corr
                endif
                dist = 1. - dist
            case(OBJFUN_EUCLID)
                if( corr < TINY )then
                    dist = huge(dist)
                else
                    dist = - log(corr)
                endif
        end select
    end function eulprob_dist_switch

    ! switch corr in [0,1] to [0, infinity) to do greedy_sampling
    elemental function eulprob_corr_switch( dist ) result(corr)
        real, intent(in) :: dist
        real :: corr
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                corr = 1 - dist
            case(OBJFUN_EUCLID)
                corr = exp(-dist)
        end select
    end function eulprob_corr_switch

    ! shift multinomal sampling within a threshold, units are in Ang
    function shift_sampling( cur_sh, thres_in ) result(sh)
        real,              intent(in) :: cur_sh(2)
        real,    optional, intent(in) :: thres_in
        integer, parameter   :: N_SMPL = 100
        real,    allocatable :: vals(:)
        logical, allocatable :: ptcl_mask(:)
        integer :: i, sh_signs(2), which, dim
        real    :: sh(2)
        real    :: d_sh, gauss_sh(N_SMPL), sh_vals(N_SMPL), sig2, d_thres, thres
        if( present(thres_in) )then
            thres = thres_in
        else
            ptcl_mask = nint(build_glob%spproj_field%get_all('state')) == 1
            vals      = build_glob%spproj_field%get_all(trim('shincarg'))
            thres     = sum(vals, mask=ptcl_mask) / real(count(ptcl_mask))
            thres     = thres / 2.  ! be more aggressive for convergence
        endif
        sh = cur_sh
        if( thres < TINY ) return
        ! randomly pick the plus/minus for each x,y dimensions
        sh_signs = 1
        if( ran3() < 0.5 ) sh_signs(1) = -1
        if( ran3() < 0.5 ) sh_signs(2) = -1
        ! sampling for x
        d_sh = thres / real(N_SMPL-1)
        sig2 = thres**2.
        do dim = 1, 2
            do i = 1, N_SMPL
                d_thres     = d_sh * (i - 1)
                gauss_sh(i) = gaussian1D(d_thres, avg=0., sigma_sq=sig2)
                sh_vals(i)  = d_thres
            enddo
            gauss_sh = gauss_sh / sum(gauss_sh)
            which    = multinomal( gauss_sh )
            sh(dim)  = cur_sh(dim) + real(sh_signs(dim)) * sh_vals(which)
        enddo
    end function shift_sampling

end module simple_eul_prob_tab
