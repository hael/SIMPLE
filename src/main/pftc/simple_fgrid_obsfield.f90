!@descr: spherical-shell Fourier observation field
module simple_fgrid_obsfield
use simple_core_module_api
use simple_shell_field_geom, only: shell_field_geom
implicit none

public :: fgrid_obs_field, fgrid_obsfield_eo, unsampled_floor
private
#include "simple_local_flags.inc"

integer,  parameter :: OBSFIELD_HEADER_SIZE       = 21
integer,  parameter :: OBSFIELD_FORMAT_VERSION    = 1
integer,  parameter :: OBSFIELD_SHELL_NNODES      = 8
real(dp), parameter :: OBSFIELD_UNSAMPLED_FLOOR   = 1.d3

type :: fgrid_obs_field
    private
    type(shell_field_geom) :: shell_geom
    integer                :: pf            = OSMPL_PAD_FAC
    integer                :: nyq           = 0
    integer                :: lims(3,2)     = 0
    integer                :: nobs          = 0
    integer                :: ncells        = 0
    integer                :: shell_kfromto(2) = 0
    integer                :: shell_nodes      = 0
    logical                :: initialized   = .false.
    logical                :: shell_initialized = .false.
    logical                :: restored      = .false.
    complex(dp), allocatable :: shell_num(:)
    real(dp),    allocatable :: shell_den(:)
    logical,     allocatable :: shell_obs(:)
    integer,     allocatable :: shell_ncells(:)
    integer,     allocatable :: shell_id(:)
    logical,     allocatable :: shell_any_obs(:)
    integer,     allocatable :: shell_obs_cnt(:)
  contains
    procedure, public  :: new                   => obsfield_new
    procedure, public  :: init_shell_cache      => obsfield_init_shell_cache
    procedure, public  :: copy_layout_from      => obsfield_copy_layout_from
    procedure, public  :: reset                 => obsfield_reset
    procedure, public  :: kill                  => obsfield_kill
    procedure, public  :: insert_plane_oversamp
    procedure, public  :: append_field          => obsfield_append_field
    procedure, public  :: append_field_array    => obsfield_append_field_array
    procedure, public  :: restore_field         => obsfield_restore_field
    procedure, public  :: extract_restored_shell_cache_polar => obsfield_extract_restored_shell_cache_polar
    procedure, public  :: count_shell_cells     => obsfield_count_shell_cells
    procedure, public  :: count_shell_cell_counts => obsfield_count_shell_cell_counts
    procedure, public  :: get_nobs              => obsfield_get_nobs
    procedure, public  :: get_ncells            => obsfield_get_ncells
    procedure, private :: compatible_with       => obsfield_compatible_with
    procedure, private :: clear_shell_cache     => obsfield_clear_shell_cache
    procedure, private :: zero_shell_cache      => obsfield_zero_shell_cache
    procedure, private :: refresh_shell_obs_counts => obsfield_refresh_shell_obs_counts
end type fgrid_obs_field

type :: fgrid_obsfield_eo
    type(fgrid_obs_field) :: even
    type(fgrid_obs_field) :: odd
  contains
    procedure, public :: new                   => obsfield_eo_new
    procedure, public :: init_shell_cache      => obsfield_eo_init_shell_cache
    procedure, public :: copy_layout_from      => obsfield_eo_copy_layout_from
    procedure, public :: reset                 => obsfield_eo_reset
    procedure, public :: kill                  => obsfield_eo_kill
    procedure, public :: insert_plane          => obsfield_eo_insert_plane
    procedure, public :: append_field          => obsfield_eo_append_field
    procedure, public :: append_field_array    => obsfield_eo_append_field_array
    procedure, public :: restore_field         => obsfield_eo_restore_field
    procedure, public :: extract_restored_shell_cache_polar => obsfield_eo_extract_restored_shell_cache_polar
    procedure, public :: calc_shell_fsc        => obsfield_eo_calc_shell_fsc
    procedure, public :: write                 => obsfield_eo_write
    procedure, public :: read                  => obsfield_eo_read
end type fgrid_obsfield_eo

contains

    pure real(dp) function unsampled_floor( den )
        real(dp), intent(in) :: den
        unsampled_floor = min(OBSFIELD_UNSAMPLED_FLOOR, OBSFIELD_UNSAMPLED_FLOOR * den)
    end function unsampled_floor

    subroutine obsfield_new( self, lims, nyq )
        class(fgrid_obs_field), intent(inout) :: self
        integer,                intent(in)    :: lims(3,2), nyq
        call self%kill
        if( nyq < 1 ) THROW_HARD('invalid nyq; obsfield_new')
        if( any(lims(:,2) < lims(:,1)) ) THROW_HARD('invalid limits; obsfield_new')
        self%pf          = OSMPL_PAD_FAC
        self%nyq         = nyq
        self%lims        = lims
        self%initialized = .true.
        self%restored    = .false.
    end subroutine obsfield_new

    subroutine obsfield_init_shell_cache( self, geom )
        class(fgrid_obs_field), intent(inout) :: self
        type(shell_field_geom), intent(in)    :: geom
        integer :: shell, ik, nk, node_first, node_last, geom_kfromto(2), total_nodes
        if( .not. self%initialized ) THROW_HARD('obsfield not initialized; obsfield_init_shell_cache')
        if( .not. geom%is_initialized() ) THROW_HARD('shell geometry not initialized; obsfield_init_shell_cache')
        call geom%get_kfromto(geom_kfromto)
        nk = geom_kfromto(2) - geom_kfromto(1) + 1
        if( nk < 1 ) THROW_HARD('invalid shell range; obsfield_init_shell_cache')
        total_nodes = geom%get_total_nodes()
        if( total_nodes < 1 ) THROW_HARD('empty shell geometry; obsfield_init_shell_cache')
        call self%clear_shell_cache
        self%shell_geom    = geom
        self%shell_kfromto = geom_kfromto
        self%shell_nodes   = total_nodes
        self%ncells        = total_nodes
        allocate(self%shell_num(total_nodes),    source=DCMPLX_ZERO)
        allocate(self%shell_den(total_nodes),    source=0.d0)
        allocate(self%shell_obs(total_nodes),    source=.false.)
        allocate(self%shell_ncells(total_nodes), source=0)
        allocate(self%shell_id(total_nodes),     source=0)
        allocate(self%shell_any_obs(nk),         source=.false.)
        allocate(self%shell_obs_cnt(nk),         source=0)
        do ik = 1,nk
            shell = geom_kfromto(1) + ik - 1
            call geom%get_shell_node_range(shell, node_first, node_last)
            if( node_first <= node_last ) self%shell_id(node_first:node_last) = shell
        enddo
        self%shell_initialized = .true.
        self%restored = .false.
        self%nobs = 0
    end subroutine obsfield_init_shell_cache

    subroutine obsfield_copy_layout_from( self, src )
        class(fgrid_obs_field), intent(inout) :: self
        class(fgrid_obs_field), intent(in)    :: src
        if( .not. src%initialized ) THROW_HARD('source obsfield not initialized; obsfield_copy_layout_from')
        if( .not. src%shell_initialized ) THROW_HARD('source shell obsfield not initialized; obsfield_copy_layout_from')
        call self%new(src%lims, src%nyq)
        self%shell_geom    = src%shell_geom
        self%shell_kfromto = src%shell_kfromto
        self%shell_nodes   = src%shell_nodes
        self%ncells        = src%ncells
        allocate(self%shell_num(self%shell_nodes),    source=DCMPLX_ZERO)
        allocate(self%shell_den(self%shell_nodes),    source=0.d0)
        allocate(self%shell_obs(self%shell_nodes),    source=.false.)
        allocate(self%shell_ncells(self%shell_nodes), source=0)
        allocate(self%shell_id(self%shell_nodes),     source=src%shell_id)
        allocate(self%shell_any_obs(size(src%shell_any_obs)), source=.false.)
        allocate(self%shell_obs_cnt(size(src%shell_obs_cnt)), source=0)
        self%shell_initialized = .true.
        self%restored = .false.
        self%nobs = 0
    end subroutine obsfield_copy_layout_from

    subroutine obsfield_reset( self )
        class(fgrid_obs_field), intent(inout) :: self
        self%nobs = 0
        self%restored = .false.
        call self%zero_shell_cache
    end subroutine obsfield_reset

    subroutine obsfield_kill( self )
        class(fgrid_obs_field), intent(inout) :: self
        call self%clear_shell_cache
        self%pf          = OSMPL_PAD_FAC
        self%nyq         = 0
        self%lims        = 0
        self%nobs        = 0
        self%ncells      = 0
        self%initialized = .false.
        self%restored    = .false.
    end subroutine obsfield_kill

    subroutine insert_plane_oversamp( self, se, o, fpl, pwght )
        use simple_math,    only: ceil_div, floor_div
        use simple_math_ft, only: fplane_get_cmplx, fplane_get_ctfsq
        class(fgrid_obs_field), intent(inout) :: self
        class(sym),             intent(inout) :: se
        class(ori),             intent(inout) :: o
        class(fplane_type),     intent(in)    :: fpl
        real,                   intent(in)    :: pwght
        type(ori)   :: o_sym
        complex(sp) :: cmplx_raw
        complex(dp) :: comp
        real(dp)    :: ctfval, pwght_dp, pwght_pf2_dp, node_w(OBSFIELD_SHELL_NNODES)
        real(sp)    :: ctfsq_raw
        real(sp)    :: loc(3), q(3), R(3,3)
        real        :: rotmats(se%get_nsym(),3,3)
        integer     :: node_ids(OBSFIELD_SHELL_NNODES)
        integer     :: fpllims_pd(3,2), fpllims(3,2)
        integer     :: nsym, isym, h, k, hp, kp, pf_local, shell, inode, nfound, j, ik
        integer     :: nyq_disk, h_sq, k_max_h, k_lo, k_hi, nobs_add
        logical     :: l_conj
        if( .not. self%initialized ) THROW_HARD('obsfield not initialized; insert_plane_oversamp')
        if( .not. self%shell_initialized ) THROW_HARD('shell obsfield not initialized; insert_plane_oversamp')
        if( .not. self%shell_geom%is_initialized() ) THROW_HARD('shell geometry not initialized; insert_plane_oversamp')
        if( pwght < TINY ) return
        self%restored = .false.
        nobs_add = 0
        nsym = se%get_nsym()
        rotmats(1,:,:) = o%get_mat()
        if( nsym > 1 )then
            do isym = 2, nsym
                call se%apply(o, isym, o_sym)
                rotmats(isym,:,:) = o_sym%get_mat()
            enddo
        endif
        pf_local     = self%pf
        fpllims_pd   = fpl%frlims
        fpllims      = fpllims_pd
        fpllims(1,1) = ceil_div (fpllims_pd(1,1), pf_local)
        fpllims(1,2) = floor_div(fpllims_pd(1,2), pf_local)
        fpllims(2,1) = ceil_div (fpllims_pd(2,1), pf_local)
        fpllims(2,2) = floor_div(fpllims_pd(2,2), pf_local)
        pwght_dp     = real(pwght, dp)
        pwght_pf2_dp = pwght_dp * real(pf_local*pf_local, dp)
        nyq_disk     = self%nyq * (self%nyq + 1)
        do isym = 1, nsym
            R = rotmats(isym,:,:)
            do h = fpllims(1,1), fpllims(1,2)
                h_sq = h*h
                if( h_sq > nyq_disk ) cycle
                k_max_h = int(sqrt(real(nyq_disk - h_sq, sp)))
                k_lo    = max(fpllims(2,1), -k_max_h)
                k_hi    = min(fpllims(2,2),  k_max_h)
                hp = h * pf_local
                do k = k_lo, k_hi
                    shell = nint(sqrt(real(h*h + k*k, dp)))
                    if( shell < self%shell_kfromto(1) .or. shell > self%shell_kfromto(2) ) cycle
                    kp = k * pf_local
                    cmplx_raw = fplane_get_cmplx(fpl, hp, kp)
                    ctfsq_raw = fplane_get_ctfsq(fpl, hp, kp)
                    if( abs(real(cmplx_raw)) + abs(aimag(cmplx_raw)) <= TINY .and. ctfsq_raw <= TINY ) cycle
                    loc(1) = real(h,sp)*R(1,1) + real(k,sp)*R(2,1)
                    loc(2) = real(h,sp)*R(1,2) + real(k,sp)*R(2,2)
                    loc(3) = real(h,sp)*R(1,3) + real(k,sp)*R(2,3)
                    l_conj = loc(1) < real(self%lims(1,1),sp)
                    if( l_conj )then
                        q    = -loc
                        comp = pwght_pf2_dp * conjg(cmplx(cmplx_raw, kind=dp))
                    else
                        q    = loc
                        comp = pwght_pf2_dp * cmplx(cmplx_raw, kind=dp)
                    endif
                    ctfval = pwght_dp * real(ctfsq_raw, dp)
                    call self%shell_geom%nearest_nodes(shell, q, OBSFIELD_SHELL_NNODES, node_ids, node_w, nfound)
                    ik = shell - self%shell_kfromto(1) + 1
                    do j = 1,nfound
                        inode = node_ids(j)
                        if( inode < 1 .or. inode > self%shell_nodes ) cycle
                        self%shell_num(inode) = self%shell_num(inode) + comp * node_w(j)
                        self%shell_den(inode) = self%shell_den(inode) + ctfval * node_w(j)
                        self%shell_ncells(inode) = self%shell_ncells(inode) + 1
                        if( .not. self%shell_obs(inode) ) self%shell_obs_cnt(ik) = self%shell_obs_cnt(ik) + 1
                        self%shell_obs(inode) = .true.
                        self%shell_any_obs(ik) = .true.
                        nobs_add = nobs_add + 1
                    enddo
                enddo
            enddo
        enddo
        self%nobs = self%nobs + nobs_add
        if( nsym > 1 ) call o_sym%kill
    end subroutine insert_plane_oversamp

    subroutine obsfield_append_field( self, src )
        class(fgrid_obs_field), intent(inout) :: self
        class(fgrid_obs_field), intent(in)    :: src
        if( .not. src%initialized ) return
        if( .not. self%initialized ) THROW_HARD('destination not initialized; obsfield_append_field')
        if( .not. self%compatible_with(src) ) THROW_HARD('incompatible shell observation fields; obsfield_append_field')
        if( self%restored .or. src%restored ) THROW_HARD('cannot append restored observation fields; obsfield_append_field')
        if( src%nobs < 1 ) return
        self%shell_num    = self%shell_num + src%shell_num
        self%shell_den    = self%shell_den + src%shell_den
        self%shell_ncells = self%shell_ncells + src%shell_ncells
        self%shell_obs    = self%shell_obs .or. src%shell_obs
        self%nobs = self%nobs + src%nobs
        call self%refresh_shell_obs_counts
        self%restored = .false.
    end subroutine obsfield_append_field

    subroutine obsfield_append_field_array( self, src )
        class(fgrid_obs_field), intent(inout) :: self
        class(fgrid_obs_field), intent(in)    :: src(:)
        complex(dp) :: num
        real(dp)    :: den
        logical     :: obs
        integer     :: isrc, node, ncells, nobs_add
        if( .not. self%initialized ) THROW_HARD('destination not initialized; obsfield_append_field_array')
        nobs_add = 0
        do isrc = 1,size(src)
            if( .not. src(isrc)%initialized ) cycle
            if( .not. self%compatible_with(src(isrc)) ) THROW_HARD('incompatible shell observation fields; obsfield_append_field_array')
            if( self%restored .or. src(isrc)%restored ) THROW_HARD('cannot append restored observation fields; obsfield_append_field_array')
            nobs_add = nobs_add + src(isrc)%nobs
        enddo
        if( nobs_add < 1 ) return
        !$omp parallel do default(shared) schedule(static) private(node,isrc,num,den,ncells,obs) proc_bind(close)
        do node = 1,self%shell_nodes
            num    = self%shell_num(node)
            den    = self%shell_den(node)
            ncells = self%shell_ncells(node)
            obs    = self%shell_obs(node)
            do isrc = 1,size(src)
                if( src(isrc)%nobs < 1 ) cycle
                num    = num + src(isrc)%shell_num(node)
                den    = den + src(isrc)%shell_den(node)
                ncells = ncells + src(isrc)%shell_ncells(node)
                obs    = obs .or. src(isrc)%shell_obs(node)
            enddo
            self%shell_num(node)    = num
            self%shell_den(node)    = den
            self%shell_ncells(node) = ncells
            self%shell_obs(node)    = obs
        enddo
        !$omp end parallel do
        self%nobs = self%nobs + nobs_add
        call self%refresh_shell_obs_counts
        self%restored = .false.
    end subroutine obsfield_append_field_array

    subroutine obsfield_restore_field( self, kfromto, invtau2, prior_start )
        class(fgrid_obs_field), intent(inout) :: self
        integer,                intent(in)    :: kfromto(2), prior_start
        real(dp),               intent(in)    :: invtau2(kfromto(1):kfromto(2))
        real(dp) :: denom, prior
        integer  :: node, shell, nobs_new
        if( .not. self%initialized ) return
        if( .not. self%shell_initialized ) THROW_HARD('shell obsfield not initialized; obsfield_restore_field')
        if( any(self%shell_kfromto /= kfromto) ) THROW_HARD('shell k-range mismatch; obsfield_restore_field')
        nobs_new = 0
        !$omp parallel do default(shared) schedule(static) private(node,shell,prior,denom) reduction(+:nobs_new) proc_bind(close)
        do node = 1,self%shell_nodes
            if( .not. self%shell_obs(node) ) cycle
            shell = self%shell_id(node)
            prior = 0.d0
            if( shell >= kfromto(1) .and. shell <= kfromto(2) ) prior = invtau2(shell)
            denom = self%shell_den(node) + prior
            if( shell >= kfromto(1) .and. shell <= kfromto(2) .and. shell >= prior_start .and. prior <= DTINY )then
                denom = denom + unsampled_floor(self%shell_den(node))
            endif
            if( denom <= DTINY )then
                self%shell_num(node) = DCMPLX_ZERO
                self%shell_den(node) = 0.d0
                self%shell_obs(node) = .false.
                self%shell_ncells(node) = 0
            else
                self%shell_num(node) = self%shell_num(node) / denom
                self%shell_den(node) = denom
                nobs_new = nobs_new + self%shell_ncells(node)
            endif
        enddo
        !$omp end parallel do
        self%nobs = nobs_new
        call self%refresh_shell_obs_counts
        self%restored = .true.
    end subroutine obsfield_restore_field

    subroutine obsfield_extract_restored_shell_cache_polar( self, eulspace, nrefs, kfromto, polar_x, polar_y, pfts, denw )
        class(fgrid_obs_field), intent(in)    :: self
        class(oris),            intent(in)    :: eulspace
        integer,                intent(in)    :: nrefs, kfromto(2)
        real(sp),               intent(in)    :: polar_x(:,:), polar_y(:,:)
        complex(dp),            intent(inout) :: pfts(:,:,:)
        real(dp), optional,     intent(inout) :: denw(:,:,:)
        real(dp) :: node_w(OBSFIELD_SHELL_NNODES), denacc
        real(sp) :: R(3,3), loc(3), q(3), px, py
        integer  :: node_ids(OBSFIELD_SHELL_NNODES)
        integer  :: iproj, irot, k, kloc, pftsz, inode, nfound, j
        logical  :: l_conj
        if( .not. self%restored ) THROW_HARD('shell obsfield extraction requires restore_field')
        if( .not. self%shell_initialized ) THROW_HARD('shell obsfield not initialized; extraction')
        if( any(self%shell_kfromto /= kfromto) ) THROW_HARD('shell k-range mismatch; extraction')
        if( .not. self%shell_geom%is_initialized() ) THROW_HARD('shell geometry not initialized; extraction')
        pftsz = size(polar_x,1)
        pfts(:,:,1:nrefs) = DCMPLX_ZERO
        if( present(denw) ) denw(:,:,1:nrefs) = 0.d0
        !$omp parallel do default(shared) schedule(static) private(iproj,R,k,kloc,irot,px,py,loc,q,l_conj,node_ids,node_w,nfound,j,inode,denacc) &
        !$omp& proc_bind(close)
        do iproj = 1,nrefs
            R = eulspace%get_mat(iproj)
            do k = kfromto(1),kfromto(2)
                kloc = k - kfromto(1) + 1
                do irot = 1,pftsz
                    px     = polar_x(irot,kloc)
                    py     = polar_y(irot,kloc)
                    loc(1) = px*R(1,1) + py*R(2,1)
                    loc(2) = px*R(1,2) + py*R(2,2)
                    loc(3) = px*R(1,3) + py*R(2,3)
                    l_conj = loc(1) < real(self%lims(1,1),sp)
                    if( l_conj )then
                        q = -loc
                    else
                        q = loc
                    endif
                    call self%shell_geom%nearest_nodes(k, q, OBSFIELD_SHELL_NNODES, node_ids, node_w, nfound)
                    denacc = 0.d0
                    do j = 1,nfound
                        inode = node_ids(j)
                        if( inode < 1 .or. inode > self%shell_nodes ) cycle
                        if( .not. self%shell_obs(inode) ) cycle
                        pfts(irot,kloc,iproj) = pfts(irot,kloc,iproj) + node_w(j) * self%shell_num(inode)
                        denacc = denacc + node_w(j) * self%shell_den(inode)
                    enddo
                    if( l_conj ) pfts(irot,kloc,iproj) = conjg(pfts(irot,kloc,iproj))
                    if( present(denw) ) denw(irot,kloc,iproj) = denacc
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine obsfield_extract_restored_shell_cache_polar

    integer function obsfield_count_shell_cells( self, kfromto )
        class(fgrid_obs_field), intent(in) :: self
        integer,                intent(in) :: kfromto(2)
        integer, allocatable :: shell_cells(:)
        obsfield_count_shell_cells = 0
        call self%count_shell_cell_counts(kfromto, shell_cells)
        if( allocated(shell_cells) )then
            obsfield_count_shell_cells = sum(shell_cells)
            deallocate(shell_cells)
        endif
    end function obsfield_count_shell_cells

    subroutine obsfield_count_shell_cell_counts( self, kfromto, shell_cells )
        class(fgrid_obs_field), intent(in)  :: self
        integer,                intent(in)  :: kfromto(2)
        integer, allocatable,   intent(out) :: shell_cells(:)
        integer :: h, k, l, shell, ik, nk
        nk = kfromto(2) - kfromto(1) + 1
        if( nk < 1 )then
            allocate(shell_cells(0))
            return
        endif
        allocate(shell_cells(nk), source=0)
        if( .not. self%initialized ) return
        !$omp parallel do collapse(3) default(shared) private(h,k,l,shell,ik) reduction(+:shell_cells) schedule(static) proc_bind(close)
        do l = self%lims(3,1), self%lims(3,2)
            do k = self%lims(2,1), self%lims(2,2)
                do h = self%lims(1,1), self%lims(1,2)
                    shell = nint(sqrt(real(h*h + k*k + l*l, dp)))
                    if( shell < kfromto(1) .or. shell > kfromto(2) ) cycle
                    ik = shell - kfromto(1) + 1
                    shell_cells(ik) = shell_cells(ik) + 1
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine obsfield_count_shell_cell_counts

    integer function obsfield_get_nobs( self )
        class(fgrid_obs_field), intent(in) :: self
        obsfield_get_nobs = self%nobs
    end function obsfield_get_nobs

    integer function obsfield_get_ncells( self )
        class(fgrid_obs_field), intent(in) :: self
        obsfield_get_ncells = self%ncells
    end function obsfield_get_ncells

    subroutine obsfield_clear_shell_cache( self )
        class(fgrid_obs_field), intent(inout) :: self
        if( allocated(self%shell_num) ) deallocate(self%shell_num)
        if( allocated(self%shell_den) ) deallocate(self%shell_den)
        if( allocated(self%shell_obs) ) deallocate(self%shell_obs)
        if( allocated(self%shell_ncells) ) deallocate(self%shell_ncells)
        if( allocated(self%shell_id) ) deallocate(self%shell_id)
        if( allocated(self%shell_any_obs) ) deallocate(self%shell_any_obs)
        if( allocated(self%shell_obs_cnt) ) deallocate(self%shell_obs_cnt)
        call self%shell_geom%kill
        self%shell_kfromto = 0
        self%shell_nodes = 0
        self%shell_initialized = .false.
    end subroutine obsfield_clear_shell_cache

    subroutine obsfield_zero_shell_cache( self )
        class(fgrid_obs_field), intent(inout) :: self
        if( allocated(self%shell_num) ) self%shell_num = DCMPLX_ZERO
        if( allocated(self%shell_den) ) self%shell_den = 0.d0
        if( allocated(self%shell_obs) ) self%shell_obs = .false.
        if( allocated(self%shell_ncells) ) self%shell_ncells = 0
        if( allocated(self%shell_any_obs) ) self%shell_any_obs = .false.
        if( allocated(self%shell_obs_cnt) ) self%shell_obs_cnt = 0
    end subroutine obsfield_zero_shell_cache

    subroutine obsfield_refresh_shell_obs_counts( self )
        class(fgrid_obs_field), intent(inout) :: self
        integer :: ik, node
        if( .not. allocated(self%shell_any_obs) ) return
        self%shell_any_obs = .false.
        self%shell_obs_cnt = 0
        do node = 1,self%shell_nodes
            if( .not. self%shell_obs(node) ) cycle
            ik = self%shell_id(node) - self%shell_kfromto(1) + 1
            if( ik < 1 .or. ik > size(self%shell_any_obs) ) cycle
            self%shell_obs_cnt(ik) = self%shell_obs_cnt(ik) + 1
        enddo
        self%shell_any_obs = self%shell_obs_cnt > 0
    end subroutine obsfield_refresh_shell_obs_counts

    logical function obsfield_compatible_with( self, other )
        class(fgrid_obs_field), intent(in) :: self, other
        obsfield_compatible_with = self%initialized .and. other%initialized .and. &
            &self%shell_initialized .and. other%shell_initialized .and. &
            &all(self%shell_kfromto == other%shell_kfromto) .and. self%shell_nodes == other%shell_nodes
        if( obsfield_compatible_with ) obsfield_compatible_with = all(self%shell_id == other%shell_id)
    end function obsfield_compatible_with

    subroutine obsfield_write_local( self, funit, label )
        class(fgrid_obs_field), intent(in) :: self
        integer,                intent(in) :: funit
        character(len=*),       intent(in) :: label
        integer :: header(OBSFIELD_HEADER_SIZE), io_stat
        if( .not. self%initialized ) THROW_HARD('obsfield not initialized; '//trim(label))
        if( .not. self%shell_initialized ) THROW_HARD('shell obsfield not initialized; '//trim(label))
        header = 0
        header(1:2)   = [self%pf, self%nyq]
        header(5:10)  = reshape(self%lims, [6])
        header(11:12) = self%shell_kfromto
        header(17)    = self%shell_nodes
        header(19)    = OBSFIELD_FORMAT_VERSION
        header(20:21) = [sum(self%shell_ncells), self%shell_nodes]
        write(funit, iostat=io_stat) header
        call fileiochk('obsfield_write_local shell header; '//trim(label), io_stat)
        write(funit, iostat=io_stat) self%shell_num
        call fileiochk('obsfield_write_local shell numerator; '//trim(label), io_stat)
        write(funit, iostat=io_stat) self%shell_den
        call fileiochk('obsfield_write_local shell density; '//trim(label), io_stat)
        write(funit, iostat=io_stat) self%shell_obs
        call fileiochk('obsfield_write_local shell observation gate; '//trim(label), io_stat)
        write(funit, iostat=io_stat) self%shell_ncells
        call fileiochk('obsfield_write_local shell cell counts; '//trim(label), io_stat)
        write(funit, iostat=io_stat) self%shell_id
        call fileiochk('obsfield_write_local shell ids; '//trim(label), io_stat)
    end subroutine obsfield_write_local

    subroutine obsfield_read_local( self, funit, label )
        class(fgrid_obs_field), intent(inout) :: self
        integer,                intent(in)    :: funit
        character(len=*),       intent(in)    :: label
        integer :: header(OBSFIELD_HEADER_SIZE), io_stat, lims(3,2), nnodes
        read(funit, iostat=io_stat) header
        call fileiochk('obsfield_read_local shell header; '//trim(label), io_stat)
        if( header(19) /= OBSFIELD_FORMAT_VERSION )then
            THROW_HARD('invalid obsfield shell part format: '//trim(label))
        endif
        lims   = reshape(header(5:10), [3,2])
        nnodes = header(21)
        call self%new(lims, header(2))
        self%shell_kfromto = header(11:12)
        self%shell_nodes   = nnodes
        self%ncells        = nnodes
        allocate(self%shell_num(nnodes))
        allocate(self%shell_den(nnodes))
        allocate(self%shell_obs(nnodes))
        allocate(self%shell_ncells(nnodes))
        allocate(self%shell_id(nnodes))
        allocate(self%shell_any_obs(self%shell_kfromto(2)-self%shell_kfromto(1)+1), source=.false.)
        allocate(self%shell_obs_cnt(self%shell_kfromto(2)-self%shell_kfromto(1)+1), source=0)
        read(funit, iostat=io_stat) self%shell_num
        call fileiochk('obsfield_read_local shell numerator; '//trim(label), io_stat)
        read(funit, iostat=io_stat) self%shell_den
        call fileiochk('obsfield_read_local shell density; '//trim(label), io_stat)
        read(funit, iostat=io_stat) self%shell_obs
        call fileiochk('obsfield_read_local shell observation gate; '//trim(label), io_stat)
        read(funit, iostat=io_stat) self%shell_ncells
        call fileiochk('obsfield_read_local shell cell counts; '//trim(label), io_stat)
        read(funit, iostat=io_stat) self%shell_id
        call fileiochk('obsfield_read_local shell ids; '//trim(label), io_stat)
        self%shell_initialized = .true.
        self%nobs = sum(self%shell_ncells)
        call self%refresh_shell_obs_counts
        self%restored = .false.
    end subroutine obsfield_read_local

    subroutine obsfield_eo_new( self, lims, nyq )
        class(fgrid_obsfield_eo), intent(inout) :: self
        integer,                   intent(in)    :: lims(3,2), nyq
        call self%even%new(lims, nyq)
        call self%odd%new( lims, nyq)
    end subroutine obsfield_eo_new

    subroutine obsfield_eo_init_shell_cache( self, geom )
        class(fgrid_obsfield_eo), intent(inout) :: self
        type(shell_field_geom),   intent(in)    :: geom
        call self%even%init_shell_cache(geom)
        call self%odd%init_shell_cache(geom)
    end subroutine obsfield_eo_init_shell_cache

    subroutine obsfield_eo_copy_layout_from( self, src )
        class(fgrid_obsfield_eo), intent(inout) :: self
        class(fgrid_obsfield_eo), intent(in)    :: src
        call self%even%copy_layout_from(src%even)
        call self%odd%copy_layout_from( src%odd )
    end subroutine obsfield_eo_copy_layout_from

    subroutine obsfield_eo_reset( self )
        class(fgrid_obsfield_eo), intent(inout) :: self
        call self%even%reset
        call self%odd%reset
    end subroutine obsfield_eo_reset

    subroutine obsfield_eo_kill( self )
        class(fgrid_obsfield_eo), intent(inout) :: self
        call self%even%kill
        call self%odd%kill
    end subroutine obsfield_eo_kill

    subroutine obsfield_eo_insert_plane( self, se, o, fpl, eo, pwght )
        class(fgrid_obsfield_eo), intent(inout) :: self
        class(sym),                intent(inout) :: se
        class(ori),                intent(inout) :: o
        class(fplane_type),        intent(in)    :: fpl
        integer,                   intent(in)    :: eo
        real,                      intent(in)    :: pwght
        select case(eo)
            case(-1,0)
                call self%even%insert_plane_oversamp(se, o, fpl, pwght)
            case(1)
                call self%odd%insert_plane_oversamp(se, o, fpl, pwght)
            case DEFAULT
                THROW_HARD('unsupported eo flag; obsfield_eo_insert_plane')
        end select
    end subroutine obsfield_eo_insert_plane

    subroutine obsfield_eo_append_field( self, src )
        class(fgrid_obsfield_eo), intent(inout) :: self
        class(fgrid_obsfield_eo), intent(in)    :: src
        call self%even%append_field(src%even)
        call self%odd%append_field(src%odd)
    end subroutine obsfield_eo_append_field

    subroutine obsfield_eo_append_field_array( self, src )
        class(fgrid_obsfield_eo), intent(inout) :: self
        class(fgrid_obsfield_eo), intent(in)    :: src(:)
        call obsfield_append_eo_half_array(self%even, src, .false.)
        call obsfield_append_eo_half_array(self%odd,  src, .true.)
    end subroutine obsfield_eo_append_field_array

    subroutine obsfield_append_eo_half_array( dest, src, l_odd )
        class(fgrid_obs_field),    intent(inout) :: dest
        class(fgrid_obsfield_eo),  intent(in)    :: src(:)
        logical,                   intent(in)    :: l_odd
        complex(dp) :: num
        real(dp)    :: den
        logical     :: obs
        integer     :: isrc, node, ncells, nobs_add
        if( .not. dest%initialized ) THROW_HARD('destination not initialized; obsfield_append_eo_half_array')
        nobs_add = 0
        if( l_odd )then
            do isrc = 1,size(src)
                if( .not. src(isrc)%odd%initialized ) cycle
                if( .not. dest%compatible_with(src(isrc)%odd) ) THROW_HARD('incompatible odd shell observation fields; obsfield_append_eo_half_array')
                if( dest%restored .or. src(isrc)%odd%restored ) THROW_HARD('cannot append restored odd observation fields; obsfield_append_eo_half_array')
                nobs_add = nobs_add + src(isrc)%odd%nobs
            enddo
            if( nobs_add < 1 ) return
            !$omp parallel do default(shared) schedule(static) private(node,isrc,num,den,ncells,obs) proc_bind(close)
            do node = 1,dest%shell_nodes
                num    = dest%shell_num(node)
                den    = dest%shell_den(node)
                ncells = dest%shell_ncells(node)
                obs    = dest%shell_obs(node)
                do isrc = 1,size(src)
                    if( src(isrc)%odd%nobs < 1 ) cycle
                    num    = num + src(isrc)%odd%shell_num(node)
                    den    = den + src(isrc)%odd%shell_den(node)
                    ncells = ncells + src(isrc)%odd%shell_ncells(node)
                    obs    = obs .or. src(isrc)%odd%shell_obs(node)
                enddo
                dest%shell_num(node)    = num
                dest%shell_den(node)    = den
                dest%shell_ncells(node) = ncells
                dest%shell_obs(node)    = obs
            enddo
            !$omp end parallel do
        else
            do isrc = 1,size(src)
                if( .not. src(isrc)%even%initialized ) cycle
                if( .not. dest%compatible_with(src(isrc)%even) ) THROW_HARD('incompatible even shell observation fields; obsfield_append_eo_half_array')
                if( dest%restored .or. src(isrc)%even%restored ) THROW_HARD('cannot append restored even observation fields; obsfield_append_eo_half_array')
                nobs_add = nobs_add + src(isrc)%even%nobs
            enddo
            if( nobs_add < 1 ) return
            !$omp parallel do default(shared) schedule(static) private(node,isrc,num,den,ncells,obs) proc_bind(close)
            do node = 1,dest%shell_nodes
                num    = dest%shell_num(node)
                den    = dest%shell_den(node)
                ncells = dest%shell_ncells(node)
                obs    = dest%shell_obs(node)
                do isrc = 1,size(src)
                    if( src(isrc)%even%nobs < 1 ) cycle
                    num    = num + src(isrc)%even%shell_num(node)
                    den    = den + src(isrc)%even%shell_den(node)
                    ncells = ncells + src(isrc)%even%shell_ncells(node)
                    obs    = obs .or. src(isrc)%even%shell_obs(node)
                enddo
                dest%shell_num(node)    = num
                dest%shell_den(node)    = den
                dest%shell_ncells(node) = ncells
                dest%shell_obs(node)    = obs
            enddo
            !$omp end parallel do
        endif
        dest%nobs = dest%nobs + nobs_add
        call dest%refresh_shell_obs_counts
        dest%restored = .false.
    end subroutine obsfield_append_eo_half_array

    subroutine obsfield_eo_restore_field( self, kfromto, invtau2_even, invtau2_odd, prior_start )
        class(fgrid_obsfield_eo), intent(inout) :: self
        integer,                  intent(in)    :: kfromto(2), prior_start
        real(dp),                 intent(in)    :: invtau2_even(kfromto(1):kfromto(2))
        real(dp),                 intent(in)    :: invtau2_odd( kfromto(1):kfromto(2))
        call self%even%restore_field(kfromto, invtau2_even, prior_start)
        call self%odd%restore_field( kfromto, invtau2_odd,  prior_start)
    end subroutine obsfield_eo_restore_field

    subroutine obsfield_eo_extract_restored_shell_cache_polar( self, eulspace, nrefs, kfromto, polar_x, polar_y, &
            &pfts_even, pfts_odd, pfts_merg )
        class(fgrid_obsfield_eo), intent(in)    :: self
        class(oris),              intent(in)    :: eulspace
        integer,                  intent(in)    :: nrefs, kfromto(2)
        real(sp),                 intent(in)    :: polar_x(:,:), polar_y(:,:)
        complex(dp),              intent(inout) :: pfts_even(:,:,:), pfts_odd(:,:,:), pfts_merg(:,:,:)
        real(dp)    :: node_w(OBSFIELD_SHELL_NNODES), den_even, den_odd, denom_weight
        complex(dp) :: acc_even, acc_odd
        real(sp)    :: R(3,3), loc(3), q(3), px, py
        integer     :: node_ids(OBSFIELD_SHELL_NNODES)
        integer     :: iproj, irot, k, ik, pftsz, inode, nfound, j
        logical     :: l_conj
        if( .not. self%even%restored .or. .not. self%odd%restored ) THROW_HARD('shell obsfield extraction requires restored fields')
        if( .not. self%even%compatible_with(self%odd) ) THROW_HARD('incompatible even/odd shell fields; obsfield extraction')
        if( any(self%even%shell_kfromto /= kfromto) ) THROW_HARD('even shell k-range mismatch; obsfield extraction')
        if( any(self%odd%shell_kfromto  /= kfromto) ) THROW_HARD('odd shell k-range mismatch; obsfield extraction')
        if( .not. self%even%shell_geom%is_initialized() ) THROW_HARD('shell geometry not initialized; obsfield extraction')
        pftsz = size(polar_x,1)
        pfts_even(:,:,1:nrefs) = DCMPLX_ZERO
        pfts_odd( :,:,1:nrefs) = DCMPLX_ZERO
        pfts_merg(:,:,1:nrefs) = DCMPLX_ZERO
        !$omp parallel do default(shared) schedule(static) private(iproj,R,k,ik,irot,px,py,loc,q,l_conj,node_ids,node_w,nfound,j,inode,acc_even,acc_odd,den_even,den_odd,denom_weight) &
        !$omp& proc_bind(close)
        do iproj = 1,nrefs
            R = eulspace%get_mat(iproj)
            do k = kfromto(1),kfromto(2)
                ik = k - kfromto(1) + 1
                do irot = 1,pftsz
                    px     = polar_x(irot,ik)
                    py     = polar_y(irot,ik)
                    loc(1) = px*R(1,1) + py*R(2,1)
                    loc(2) = px*R(1,2) + py*R(2,2)
                    loc(3) = px*R(1,3) + py*R(2,3)
                    l_conj = loc(1) < real(self%even%lims(1,1),sp)
                    if( l_conj )then
                        q = -loc
                    else
                        q = loc
                    endif
                    call self%even%shell_geom%nearest_nodes(k, q, OBSFIELD_SHELL_NNODES, node_ids, node_w, nfound)
                    acc_even = DCMPLX_ZERO
                    acc_odd  = DCMPLX_ZERO
                    den_even = 0.d0
                    den_odd  = 0.d0
                    do j = 1,nfound
                        inode = node_ids(j)
                        if( inode < 1 .or. inode > self%even%shell_nodes ) cycle
                        if( self%even%shell_obs(inode) )then
                            acc_even = acc_even + node_w(j) * self%even%shell_num(inode)
                            den_even = den_even + node_w(j) * self%even%shell_den(inode)
                        endif
                        if( self%odd%shell_obs(inode) )then
                            acc_odd = acc_odd + node_w(j) * self%odd%shell_num(inode)
                            den_odd = den_odd + node_w(j) * self%odd%shell_den(inode)
                        endif
                    enddo
                    if( l_conj )then
                        acc_even = conjg(acc_even)
                        acc_odd  = conjg(acc_odd)
                    endif
                    pfts_even(irot,ik,iproj) = acc_even
                    pfts_odd( irot,ik,iproj) = acc_odd
                    denom_weight = den_even + den_odd
                    if( denom_weight <= DTINY ) cycle
                    pfts_merg(irot,ik,iproj) = (acc_even * den_even + acc_odd * den_odd) / denom_weight
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine obsfield_eo_extract_restored_shell_cache_polar

    subroutine obsfield_eo_calc_shell_fsc( self, kfromto, fsc )
        class(fgrid_obsfield_eo), intent(in)  :: self
        integer,                  intent(in)  :: kfromto(2)
        real(dp),                 intent(out) :: fsc(kfromto(1):kfromto(2))
        real(dp) :: vare(kfromto(1):kfromto(2)), varo(kfromto(1):kfromto(2))
        integer  :: node, shell
        if( .not. self%even%restored .or. .not. self%odd%restored ) THROW_HARD('shell FSC requires restored obsfields')
        if( .not. self%even%compatible_with(self%odd) ) THROW_HARD('incompatible even/odd shell fields; obsfield_eo_calc_shell_fsc')
        fsc  = 0.d0
        vare = 0.d0
        varo = 0.d0
        do node = 1,self%even%shell_nodes
            if( .not. self%even%shell_obs(node) .or. .not. self%odd%shell_obs(node) ) cycle
            shell = self%even%shell_id(node)
            if( shell < kfromto(1) .or. shell > kfromto(2) ) cycle
            fsc(shell)  = fsc(shell)  + real(self%even%shell_num(node) * conjg(self%odd%shell_num(node)), dp)
            vare(shell) = vare(shell) + real(self%even%shell_num(node) * conjg(self%even%shell_num(node)), dp)
            varo(shell) = varo(shell) + real(self%odd%shell_num(node)  * conjg(self%odd%shell_num(node)),  dp)
        enddo
        vare = vare * varo
        where( vare > DTINY )
            fsc = fsc / sqrt(vare)
        elsewhere
            fsc = 0.d0
        end where
    end subroutine obsfield_eo_calc_shell_fsc

    subroutine obsfield_eo_write( self, fname )
        class(fgrid_obsfield_eo), intent(in) :: self
        class(string),             intent(in) :: fname
        integer :: funit, io_stat
        call fopen(funit, fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('obsfield_eo_write open; '//fname%to_char(), io_stat)
        call obsfield_write_local(self%even, funit, 'even '//fname%to_char())
        call obsfield_write_local(self%odd,  funit, 'odd '//fname%to_char())
        call fclose(funit)
    end subroutine obsfield_eo_write

    subroutine obsfield_eo_read( self, fname )
        class(fgrid_obsfield_eo), intent(inout) :: self
        class(string),             intent(in)    :: fname
        integer :: funit, io_stat
        if( .not. file_exists(fname) ) THROW_HARD('obsfield file does not exist: '//fname%to_char())
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('obsfield_eo_read open; '//fname%to_char(), io_stat)
        call obsfield_read_local(self%even, funit, 'even '//fname%to_char())
        call obsfield_read_local(self%odd,  funit, 'odd '//fname%to_char())
        call fclose(funit)
    end subroutine obsfield_eo_read

end module simple_fgrid_obsfield
