!@descr: stage-independent spherical-shell observation geometry
module simple_shell_field_geom
use simple_pftc_api
implicit none

public :: shell_field_geom
private
#include "simple_local_flags.inc"

real(dp), parameter :: SHELL_SUPPORT_MAX_RAD = 10.d0 * DPI / 180.d0

type :: shell_field_geom
    private
    integer :: nk         = 0
    integer :: kfromto(2) = 0
    integer :: nsym       = 1
    integer :: total_nodes = 0
    integer :: target_nodes = 0
    integer :: shell_budget = 0
    integer,  allocatable :: full_nodes(:), shell_nodes(:), shell_cell_budget(:), offsets(:)
    integer,  allocatable :: bin_nz(:), bin_nphi(:), bin_offsets(:), bin_head(:), bin_next(:)
    real(dp), allocatable :: theta_support(:)
    real(sp), allocatable :: node_dir(:,:) ! (3,total_nodes)
    real(sp), allocatable :: sym_rmat(:,:,:) ! (3,3,nsym)
  contains
    procedure :: new         => shell_geom_new
    procedure :: kill        => shell_geom_kill
    procedure :: log_summary => shell_geom_log_summary
    procedure :: log_stats   => shell_geom_log_stats
    procedure :: get_shell_nodes => shell_geom_get_shell_nodes
    procedure :: get_shell_cell_counts => shell_geom_get_shell_cell_counts
    procedure :: is_initialized  => shell_geom_is_initialized
    procedure :: get_kfromto     => shell_geom_get_kfromto
    procedure :: get_nk          => shell_geom_get_nk
    procedure :: get_nsym        => shell_geom_get_nsym
    procedure :: get_total_nodes => shell_geom_get_total_nodes
    procedure :: get_shell_node_count => shell_geom_get_shell_node_count
    procedure :: get_shell_node_range => shell_geom_get_shell_node_range
    procedure :: get_shell_support    => shell_geom_get_shell_support
    procedure :: get_shell_dirs       => shell_geom_get_shell_dirs
    procedure :: get_node_dir         => shell_geom_get_node_dir
    procedure :: nearest_node         => shell_geom_nearest_node
    procedure :: nearest_nodes        => shell_geom_nearest_nodes
end type shell_field_geom

contains

    subroutine shell_geom_new( self, symop, kfromto, shell_budget, shell_cell_counts )
        class(shell_field_geom), intent(inout) :: self
        class(sym),              intent(in)    :: symop
        integer,                 intent(in)    :: kfromto(2)
        integer,                 intent(in)    :: shell_budget
        integer, optional,        intent(in)    :: shell_cell_counts(:)
        real(sp), allocatable :: tmp_dirs(:,:,:)
        real(sp) :: rsym(3,3), q(3)
        integer, allocatable :: target_shell_nodes(:)
        integer :: isym, ik, ifull, nfull, nkeep, pos, max_full
        call self%kill
        if( kfromto(1) > kfromto(2) ) THROW_HARD('invalid k range; shell_geom_new')
        if( shell_budget < 1 ) THROW_HARD('invalid shell budget; shell_geom_new')
        self%kfromto      = kfromto
        self%nk           = kfromto(2) - kfromto(1) + 1
        self%nsym         = symop%get_nsym()
        self%shell_budget = shell_budget
        self%target_nodes = max(shell_budget, self%nk)
        allocate(self%full_nodes(self%nk),   source=0)
        allocate(self%shell_nodes(self%nk),  source=0)
        allocate(self%shell_cell_budget(self%nk), source=0)
        allocate(self%offsets(self%nk+1),    source=1)
        allocate(self%theta_support(self%nk), source=0.d0)
        allocate(self%sym_rmat(3,3,self%nsym), source=0._sp)
        do isym = 1,self%nsym
            call symop%get_sym_rmat(isym, rsym)
            self%sym_rmat(:,:,isym) = rsym
        enddo
        call distribute_shell_budget(kfromto, self%target_nodes, target_shell_nodes)
        if( present(shell_cell_counts) )then
            if( size(shell_cell_counts) /= self%nk ) THROW_HARD('invalid shell cell count size; shell_geom_new')
            self%shell_cell_budget = shell_cell_counts
        else
            self%shell_cell_budget = target_shell_nodes
        endif
        do ik = 1,self%nk
            self%full_nodes(ik) = max(12, target_shell_nodes(ik) * max(1,self%nsym))
            self%theta_support(ik) = min(SHELL_SUPPORT_MAX_RAD, sqrt(4.d0 * DPI / real(self%full_nodes(ik),dp)))
        enddo
        max_full = maxval(self%full_nodes)
        allocate(tmp_dirs(3,max_full,self%nk), source=0._sp)
        do ik = 1,self%nk
            nfull = self%full_nodes(ik)
            nkeep = 0
            do ifull = 1,nfull
                call fibonacci_dir(ifull, nfull, q)
                if( .not. canonical_under_symmetry(self, q) ) cycle
                nkeep = nkeep + 1
                tmp_dirs(:,nkeep,ik) = q
            enddo
            self%shell_nodes(ik) = nkeep
        enddo
        self%offsets(1) = 1
        do ik = 1,self%nk
            self%offsets(ik+1) = self%offsets(ik) + self%shell_nodes(ik)
        enddo
        self%total_nodes = self%offsets(self%nk+1) - 1
        allocate(self%node_dir(3,self%total_nodes), source=0._sp)
        do ik = 1,self%nk
            pos = self%offsets(ik)
            if( self%shell_nodes(ik) > 0 )then
                self%node_dir(:,pos:pos+self%shell_nodes(ik)-1) = tmp_dirs(:,:self%shell_nodes(ik),ik)
            endif
        enddo
        call build_shell_bins(self)
        deallocate(tmp_dirs)
        deallocate(target_shell_nodes)
    end subroutine shell_geom_new

    subroutine shell_geom_kill( self )
        class(shell_field_geom), intent(inout) :: self
        if( allocated(self%full_nodes) ) deallocate(self%full_nodes)
        if( allocated(self%shell_nodes) ) deallocate(self%shell_nodes)
        if( allocated(self%shell_cell_budget) ) deallocate(self%shell_cell_budget)
        if( allocated(self%offsets) ) deallocate(self%offsets)
        if( allocated(self%bin_nz) ) deallocate(self%bin_nz)
        if( allocated(self%bin_nphi) ) deallocate(self%bin_nphi)
        if( allocated(self%bin_offsets) ) deallocate(self%bin_offsets)
        if( allocated(self%bin_head) ) deallocate(self%bin_head)
        if( allocated(self%bin_next) ) deallocate(self%bin_next)
        if( allocated(self%theta_support) ) deallocate(self%theta_support)
        if( allocated(self%node_dir) ) deallocate(self%node_dir)
        if( allocated(self%sym_rmat) ) deallocate(self%sym_rmat)
        self%nk           = 0
        self%kfromto      = 0
        self%nsym         = 1
        self%total_nodes  = 0
        self%target_nodes = 0
        self%shell_budget = 0
    end subroutine shell_geom_kill

    subroutine shell_geom_log_summary( self, label )
        class(shell_field_geom), intent(in) :: self
        character(len=*),        intent(in) :: label
        real(dp) :: mem_mb, min_support, max_support, budget_ratio
        integer :: min_nodes, max_nodes, min_shell_cells, max_shell_cells, total_shell_cells
        if( .not. allocated(self%shell_nodes) ) return
        mem_mb = 0.d0
        if( allocated(self%node_dir) )then
            mem_mb = real(storage_size(self%node_dir),dp) * real(size(self%node_dir),dp) / (8.d0 * 1024.d0 * 1024.d0)
        endif
        min_nodes   = minval(self%shell_nodes)
        max_nodes   = maxval(self%shell_nodes)
        min_support = minval(self%theta_support)
        max_support = maxval(self%theta_support)
        budget_ratio = real(self%total_nodes,dp) / real(max(1,self%shell_budget),dp)
        write(logfhandle,'(A,1X,A,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,F8.3,2X,A,I0,2X,A,I0,2X,A,F8.3,2X,A,F8.3,2X,A,F10.3)') &
            &'obsfield shell-geom', trim(label), 'nk=', self%nk, 'kmin=', self%kfromto(1), 'kmax=', self%kfromto(2), &
            &'nsym=', self%nsym, 'shell_budget=', self%shell_budget, 'target_nodes=', self%target_nodes, &
            &'total_nodes=', self%total_nodes, 'budget_ratio=', budget_ratio, 'nodes_min=', min_nodes, 'nodes_max=', max_nodes, &
            &'support_min_deg=', rad2deg(real(min_support)), &
            &'support_max_deg=', rad2deg(real(max_support)), 'dir_mem_mb=', mem_mb
        if( allocated(self%shell_cell_budget) )then
            total_shell_cells = sum(self%shell_cell_budget)
            min_shell_cells   = minval(self%shell_cell_budget)
            max_shell_cells   = maxval(self%shell_cell_budget)
            write(logfhandle,'(A,1X,A,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0)') &
                &'obsfield shell-geom cache', trim(label), 'nk=', self%nk, 'kmin=', self%kfromto(1), 'kmax=', self%kfromto(2), &
                &'shell_budget=', self%shell_budget, 'shell_cells_sum=', total_shell_cells, &
                &'shell_cells_min=', min_shell_cells, 'shell_cells_max=', max_shell_cells
        endif
    end subroutine shell_geom_log_summary

    subroutine shell_geom_log_stats( self, label )
        class(shell_field_geom), intent(in) :: self
        character(len=*),        intent(in) :: label
        integer :: ik
        call self%log_summary(label)
        if( .not. allocated(self%node_dir) ) return
        do ik = 1,self%nk
            call log_shell_stats(self, ik, self%kfromto(1)+ik-1)
        enddo
    end subroutine shell_geom_log_stats

    subroutine shell_geom_get_shell_nodes( self, shell_nodes )
        class(shell_field_geom),  intent(in)  :: self
        integer, allocatable,     intent(out) :: shell_nodes(:)
        if( .not. allocated(self%shell_nodes) )then
            allocate(shell_nodes(0))
            return
        endif
        allocate(shell_nodes(size(self%shell_nodes)), source=self%shell_nodes)
    end subroutine shell_geom_get_shell_nodes

    subroutine shell_geom_get_shell_cell_counts( self, shell_cells )
        class(shell_field_geom),  intent(in)  :: self
        integer, allocatable,     intent(out) :: shell_cells(:)
        if( .not. allocated(self%shell_cell_budget) )then
            allocate(shell_cells(0))
            return
        endif
        allocate(shell_cells(size(self%shell_cell_budget)), source=self%shell_cell_budget)
    end subroutine shell_geom_get_shell_cell_counts

    logical function shell_geom_is_initialized( self )
        class(shell_field_geom), intent(in) :: self
        shell_geom_is_initialized = allocated(self%node_dir) .and. allocated(self%shell_nodes) .and. &
            &allocated(self%shell_cell_budget) .and. allocated(self%offsets) .and. allocated(self%bin_head)
    end function shell_geom_is_initialized

    subroutine shell_geom_get_kfromto( self, kfromto )
        class(shell_field_geom), intent(in)  :: self
        integer,                 intent(out) :: kfromto(2)
        kfromto = self%kfromto
    end subroutine shell_geom_get_kfromto

    integer function shell_geom_get_nk( self )
        class(shell_field_geom), intent(in) :: self
        shell_geom_get_nk = self%nk
    end function shell_geom_get_nk

    integer function shell_geom_get_nsym( self )
        class(shell_field_geom), intent(in) :: self
        shell_geom_get_nsym = self%nsym
    end function shell_geom_get_nsym

    integer function shell_geom_get_total_nodes( self )
        class(shell_field_geom), intent(in) :: self
        shell_geom_get_total_nodes = self%total_nodes
    end function shell_geom_get_total_nodes

    integer function shell_geom_get_shell_node_count( self, shell )
        class(shell_field_geom), intent(in) :: self
        integer,                 intent(in) :: shell
        integer :: ik
        shell_geom_get_shell_node_count = 0
        if( .not. allocated(self%shell_nodes) ) return
        ik = shell - self%kfromto(1) + 1
        if( ik < 1 .or. ik > self%nk ) return
        shell_geom_get_shell_node_count = self%shell_nodes(ik)
    end function shell_geom_get_shell_node_count

    real(dp) function shell_geom_get_shell_support( self, shell )
        class(shell_field_geom), intent(in) :: self
        integer,                 intent(in) :: shell
        integer :: ik
        shell_geom_get_shell_support = 0.d0
        if( .not. allocated(self%theta_support) ) return
        ik = shell - self%kfromto(1) + 1
        if( ik < 1 .or. ik > self%nk ) return
        shell_geom_get_shell_support = self%theta_support(ik)
    end function shell_geom_get_shell_support

    subroutine shell_geom_get_shell_node_range( self, shell, node_first, node_last )
        class(shell_field_geom), intent(in)  :: self
        integer,                 intent(in)  :: shell
        integer,                 intent(out) :: node_first, node_last
        integer :: ik
        node_first = 0
        node_last  = -1
        if( .not. allocated(self%offsets) ) return
        ik = shell - self%kfromto(1) + 1
        if( ik < 1 .or. ik > self%nk ) return
        node_first = self%offsets(ik)
        node_last  = self%offsets(ik+1) - 1
    end subroutine shell_geom_get_shell_node_range

    subroutine shell_geom_get_shell_dirs( self, shell, dirs )
        class(shell_field_geom), intent(in)  :: self
        integer,                 intent(in)  :: shell
        real(sp), allocatable,   intent(out) :: dirs(:,:)
        integer :: node_first, node_last, nnodes
        call self%get_shell_node_range(shell, node_first, node_last)
        nnodes = max(0, node_last - node_first + 1)
        allocate(dirs(3,nnodes), source=0._sp)
        if( nnodes > 0 ) dirs = self%node_dir(:,node_first:node_last)
    end subroutine shell_geom_get_shell_dirs

    subroutine shell_geom_get_node_dir( self, shell, inode, dir )
        class(shell_field_geom), intent(in)  :: self
        integer,                 intent(in)  :: shell, inode
        real(sp),                intent(out) :: dir(3)
        integer :: node_first, node_last, node
        dir = 0._sp
        call self%get_shell_node_range(shell, node_first, node_last)
        if( node_last < node_first ) return
        node = node_first + inode - 1
        if( node < node_first .or. node > node_last ) return
        dir = self%node_dir(:,node)
    end subroutine shell_geom_get_node_dir

    subroutine shell_geom_nearest_node( self, shell, q, inode, cos_best )
        class(shell_field_geom), intent(in)  :: self
        integer,                 intent(in)  :: shell
        real(sp),                intent(in)  :: q(3)
        integer,                 intent(out) :: inode
        real(dp),                intent(out) :: cos_best
        real(sp) :: qcanon(3)
        real(dp) :: qnorm, dotv
        integer  :: node_first, node_last, node
        inode    = 0
        cos_best = -1.d0
        if( .not. allocated(self%node_dir) ) return
        qnorm = sqrt(sum(real(q,dp) * real(q,dp)))
        if( qnorm <= DTINY ) return
        call canonicalize_dir(self, real(real(q,dp) / qnorm, sp), qcanon)
        call self%get_shell_node_range(shell, node_first, node_last)
        if( node_last < node_first ) return
        do node = node_first,node_last
            dotv = sum(real(qcanon,dp) * real(self%node_dir(:,node),dp))
            if( dotv > cos_best )then
                cos_best = dotv
                inode    = node
            endif
        enddo
        cos_best = min(1.d0, max(-1.d0, cos_best))
    end subroutine shell_geom_nearest_node

    subroutine shell_geom_nearest_nodes( self, shell, q, max_nodes, inodes, weights, nfound )
        class(shell_field_geom), intent(in)  :: self
        integer,                 intent(in)  :: shell, max_nodes
        real(sp),                intent(in)  :: q(3)
        integer,                 intent(out) :: inodes(:)
        real(dp),                intent(out) :: weights(:)
        integer,                 intent(out) :: nfound
        real(sp) :: qcanon(3)
        real(dp), allocatable :: dots(:)
        real(dp) :: qnorm, dotv, dist, sigma, cos_support, wsum
        integer  :: node_first, node_last, node, ik, j, nkeep
        inodes  = 0
        weights = 0.d0
        nfound  = 0
        if( max_nodes < 1 ) return
        nkeep = min(max_nodes, min(size(inodes), size(weights)))
        if( nkeep < 1 ) return
        if( .not. allocated(self%node_dir) ) return
        qnorm = sqrt(sum(real(q,dp) * real(q,dp)))
        if( qnorm <= DTINY ) return
        ik = shell - self%kfromto(1) + 1
        if( ik < 1 .or. ik > self%nk ) return
        call canonicalize_dir(self, real(real(q,dp) / qnorm, sp), qcanon)
        call self%get_shell_node_range(shell, node_first, node_last)
        if( node_last < node_first ) return
        allocate(dots(nkeep), source=-huge(1.d0))
        call scan_shell_bins(self, ik, qcanon, nkeep, dots, inodes)
        if( count(inodes > 0) < nkeep )then
            dots   = -huge(1.d0)
            inodes = 0
            do node = node_first,node_last
                dotv = sum(real(qcanon,dp) * real(self%node_dir(:,node),dp))
                call insert_shell_candidate(dotv, node, nkeep, dots, inodes)
            enddo
        endif
        sigma = max(self%theta_support(ik), 1.d-6)
        cos_support = cos(sigma)
        wsum = 0.d0
        do j = 1,nkeep
            if( inodes(j) < 1 ) cycle
            dotv = min(1.d0, max(-1.d0, dots(j)))
            if( dotv < cos_support .and. nfound > 0 ) cycle
            if( dotv < cos_support .and. j > 1 ) cycle
            dist = acos(dotv)
            nfound = nfound + 1
            inodes(nfound) = inodes(j)
            weights(j) = exp(-0.5d0 * (dist / sigma)**2)
            if( nfound /= j ) weights(nfound) = weights(j)
            wsum = wsum + weights(nfound)
        enddo
        if( wsum > DTINY )then
            weights(1:nfound) = weights(1:nfound) / wsum
        else if( nfound > 0 )then
            weights(1:nfound) = 1.d0 / real(nfound,dp)
        endif
        deallocate(dots)
    end subroutine shell_geom_nearest_nodes

    subroutine build_shell_bins( self )
        class(shell_field_geom), intent(inout) :: self
        real(sp) :: dir(3)
        integer :: ik, node, total_bins, ibin, iz, iphi
        allocate(self%bin_nz(self%nk), source=0)
        allocate(self%bin_nphi(self%nk), source=0)
        allocate(self%bin_offsets(self%nk+1), source=1)
        do ik = 1,self%nk
            self%bin_nz(ik)   = max(4, int(sqrt(real(max(1,self%shell_nodes(ik)),dp) / 8.d0)))
            self%bin_nphi(ik) = 2 * self%bin_nz(ik)
            self%bin_offsets(ik+1) = self%bin_offsets(ik) + self%bin_nz(ik) * self%bin_nphi(ik)
        enddo
        total_bins = self%bin_offsets(self%nk+1) - 1
        allocate(self%bin_head(total_bins), source=0)
        allocate(self%bin_next(self%total_nodes), source=0)
        do ik = 1,self%nk
            do node = self%offsets(ik),self%offsets(ik+1)-1
                dir = self%node_dir(:,node)
                call shell_dir_bin(self, ik, dir, iz, iphi)
                ibin = self%bin_offsets(ik) + (iz - 1) * self%bin_nphi(ik) + iphi - 1
                self%bin_next(node) = self%bin_head(ibin)
                self%bin_head(ibin) = node
            enddo
        enddo
    end subroutine build_shell_bins

    subroutine shell_dir_bin( self, ik, dir, iz, iphi )
        class(shell_field_geom), intent(in)  :: self
        integer,                 intent(in)  :: ik
        real(sp),                intent(in)  :: dir(3)
        integer,                 intent(out) :: iz, iphi
        real(dp) :: zfrac, phi
        zfrac = 0.5d0 * (real(dir(3),dp) + 1.d0)
        iz = min(self%bin_nz(ik), max(1, int(zfrac * real(self%bin_nz(ik),dp)) + 1))
        phi = atan2(real(dir(2),dp), real(dir(1),dp))
        if( phi < 0.d0 ) phi = phi + 2.d0 * DPI
        iphi = min(self%bin_nphi(ik), max(1, int(phi * real(self%bin_nphi(ik),dp) / (2.d0 * DPI)) + 1))
    end subroutine shell_dir_bin

    subroutine scan_shell_bins( self, ik, qcanon, nkeep, dots, inodes )
        class(shell_field_geom), intent(in)    :: self
        integer,                 intent(in)    :: ik, nkeep
        real(sp),                intent(in)    :: qcanon(3)
        real(dp),                intent(inout) :: dots(:)
        integer,                 intent(inout) :: inodes(:)
        real(dp) :: dotv
        integer :: iz0, iphi0, iz, iphi, dz, dphi, ring, ibin, node, nfilled
        if( .not. allocated(self%bin_head) ) return
        call shell_dir_bin(self, ik, qcanon, iz0, iphi0)
        do ring = 0,3
            do dz = -ring,ring
                iz = iz0 + dz
                if( iz < 1 .or. iz > self%bin_nz(ik) ) cycle
                do dphi = -ring,ring
                    if( ring > 0 .and. abs(dz) < ring .and. abs(dphi) < ring ) cycle
                    iphi = modulo(iphi0 + dphi - 1, self%bin_nphi(ik)) + 1
                    ibin = self%bin_offsets(ik) + (iz - 1) * self%bin_nphi(ik) + iphi - 1
                    node = self%bin_head(ibin)
                    do while( node > 0 )
                        dotv = sum(real(qcanon,dp) * real(self%node_dir(:,node),dp))
                        call insert_shell_candidate(dotv, node, nkeep, dots, inodes)
                        node = self%bin_next(node)
                    enddo
                enddo
            enddo
            nfilled = count(inodes > 0)
            if( nfilled >= nkeep .and. dots(nkeep) > cos(2.d0 * self%theta_support(ik)) ) exit
        enddo
    end subroutine scan_shell_bins

    subroutine insert_shell_candidate( dotv, node, nkeep, dots, inodes )
        real(dp), intent(in)    :: dotv
        integer,  intent(in)    :: node, nkeep
        real(dp), intent(inout) :: dots(:)
        integer,  intent(inout) :: inodes(:)
        integer :: slot
        if( dotv <= dots(nkeep) ) return
        slot = nkeep
        do while( slot > 1 )
            if( dotv <= dots(slot-1) ) exit
            dots(slot)   = dots(slot-1)
            inodes(slot) = inodes(slot-1)
            slot = slot - 1
        enddo
        dots(slot)   = dotv
        inodes(slot) = node
    end subroutine insert_shell_candidate

    subroutine log_shell_stats( self, ik, shell )
        class(shell_field_geom), intent(in) :: self
        integer,                 intent(in) :: ik, shell
        real(dp) :: nn_mean, nn_sumsq, nn_sdev, nn_min, nn_max, dist, dotv, area_ratio
        integer :: stride, isample, inode, node0, node1, node2, nsamples
        node1 = self%offsets(ik)
        node2 = self%offsets(ik+1) - 1
        stride = max(1, self%shell_nodes(ik) / 128)
        nsamples  = 0
        nn_mean   = 0.d0
        nn_sumsq  = 0.d0
        nn_min    = huge(nn_min)
        nn_max    = 0.d0
        !$omp parallel do default(shared) private(isample,inode,node0,dotv,dist) &
        !$omp reduction(+:nsamples,nn_mean,nn_sumsq) reduction(min:nn_min) reduction(max:nn_max) &
        !$omp schedule(static) proc_bind(close)
        do isample = 1,self%shell_nodes(ik),stride
            node0 = node1 + isample - 1
            dotv = -1.d0
            do inode = node1,node2
                if( inode == node0 ) cycle
                dotv = max(dotv, sum(real(self%node_dir(:,node0),dp) * real(self%node_dir(:,inode),dp)))
            enddo
            dotv = min(1.d0, max(-1.d0, dotv))
            dist = acos(dotv)
            nsamples = nsamples + 1
            nn_mean  = nn_mean  + dist
            nn_sumsq = nn_sumsq + dist * dist
            nn_min   = min(nn_min, dist)
            nn_max   = max(nn_max, dist)
        enddo
        !$omp end parallel do
        if( nsamples > 0 )then
            nn_mean  = nn_mean / real(nsamples,dp)
            nn_sumsq = nn_sumsq / real(nsamples,dp)
            nn_sdev  = sqrt(max(0.d0, nn_sumsq - nn_mean*nn_mean))
        else
            nn_min   = 0.d0
            nn_sdev  = 0.d0
        endif
        area_ratio = real(self%shell_nodes(ik) * max(1,self%nsym),dp) / real(max(1,self%full_nodes(ik)),dp)
        write(logfhandle,'(A,I4,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,F8.3,2X,A,F8.3,2X,A,F8.3,2X,A,F8.3,2X,A,F8.3,2X,A,F8.3)') &
            &'obsfield shell-geom shell k=', shell, 'full_nodes=', self%full_nodes(ik), 'asym_nodes=', self%shell_nodes(ik), &
            &'sampled=', nsamples, 'sym_area_ratio=', area_ratio, 'support_deg=', rad2deg(real(self%theta_support(ik))), &
            &'nn_mean_deg=', rad2deg(real(nn_mean)), 'nn_sdev_deg=', rad2deg(real(nn_sdev)), &
            &'nn_min_deg=', rad2deg(real(nn_min)), 'nn_max_deg=', rad2deg(real(nn_max))
    end subroutine log_shell_stats

    subroutine fibonacci_dir( i, n, q )
        integer,  intent(in)  :: i, n
        real(sp), intent(out) :: q(3)
        real(dp), parameter :: GOLDEN_ANGLE = DPI * (3.d0 - sqrt(5.d0))
        real(dp) :: z, r, phi
        z = 1.d0 - 2.d0 * (real(i,dp) - 0.5d0) / real(n,dp)
        r = sqrt(max(0.d0, 1.d0 - z*z))
        phi = GOLDEN_ANGLE * real(i - 1,dp)
        q = real([r*cos(phi), r*sin(phi), z], sp)
    end subroutine fibonacci_dir

    subroutine distribute_shell_budget( kfromto, total_budget, shell_budget )
        integer,              intent(in)  :: kfromto(2), total_budget
        integer, allocatable, intent(out) :: shell_budget(:)
        real(dp), allocatable :: weights(:), exact(:), frac(:)
        real(dp) :: wsum
        integer :: ik, imax, nk, shell, remain
        nk = kfromto(2) - kfromto(1) + 1
        allocate(shell_budget(nk), source=1)
        allocate(weights(nk), exact(nk), frac(nk), source=0.d0)
        do ik = 1,nk
            shell = kfromto(1) + ik - 1
            weights(ik) = real(max(1,shell),dp)**2
        enddo
        wsum = sum(weights)
        exact = real(max(total_budget,nk),dp) * weights / wsum
        shell_budget = max(1, int(floor(exact)))
        frac = exact - real(shell_budget,dp)
        remain = max(total_budget,nk) - sum(shell_budget)
        do while( remain > 0 )
            imax = maxloc(frac, dim=1)
            shell_budget(imax) = shell_budget(imax) + 1
            frac(imax) = -1.d0
            remain = remain - 1
        enddo
        deallocate(weights, exact, frac)
    end subroutine distribute_shell_budget

    logical function canonical_under_symmetry( self, q )
        class(shell_field_geom), intent(in) :: self
        real(sp),                intent(in) :: q(3)
        real(sp) :: qs(3)
        integer :: isym
        canonical_under_symmetry = .true.
        do isym = 2,self%nsym
            qs = matmul(self%sym_rmat(:,:,isym), q)
            if( lex_less(qs, q) )then
                canonical_under_symmetry = .false.
                return
            endif
        enddo
    end function canonical_under_symmetry

    subroutine canonicalize_dir( self, q, qcanon )
        class(shell_field_geom), intent(in)  :: self
        real(sp),                intent(in)  :: q(3)
        real(sp),                intent(out) :: qcanon(3)
        real(sp) :: qs(3)
        integer  :: isym
        qcanon = q
        do isym = 2,self%nsym
            qs = matmul(self%sym_rmat(:,:,isym), q)
            if( lex_less(qs, qcanon) ) qcanon = qs
        enddo
    end subroutine canonicalize_dir

    logical function lex_less( a, b )
        real(sp), intent(in) :: a(3), b(3)
        real(sp), parameter :: TOL = 1.e-6_sp
        lex_less = .false.
        if( a(3) < b(3) - TOL )then
            lex_less = .true.
        elseif( abs(a(3)-b(3)) <= TOL )then
            if( a(2) < b(2) - TOL )then
                lex_less = .true.
            elseif( abs(a(2)-b(2)) <= TOL )then
                lex_less = a(1) < b(1) - TOL
            endif
        endif
    end function lex_less

end module simple_shell_field_geom
