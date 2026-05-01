!@descr: development scaffold for stage-independent spherical-shell geometry
module simple_shell_field_geom
use simple_pftc_api
implicit none

public :: shell_field_geom
private
#include "simple_local_flags.inc"

type :: shell_field_geom
    private
    integer :: nk         = 0
    integer :: kfromto(2) = 0
    integer :: nsym       = 1
    integer :: total_nodes = 0
    integer :: target_nodes = 0
    integer :: cart_budget  = 0
    integer,  allocatable :: full_nodes(:), shell_nodes(:), offsets(:)
    real(dp), allocatable :: theta_support(:)
    real(sp), allocatable :: node_dir(:,:) ! (3,total_nodes)
    real(sp), allocatable :: sym_rmat(:,:,:) ! (3,3,nsym)
  contains
    procedure :: new         => shell_geom_new
    procedure :: kill        => shell_geom_kill
    procedure :: log_summary => shell_geom_log_summary
    procedure :: log_stats   => shell_geom_log_stats
    procedure :: get_shell_nodes => shell_geom_get_shell_nodes
end type shell_field_geom

contains

    subroutine shell_geom_new( self, symop, kfromto, cart_budget )
        class(shell_field_geom), intent(inout) :: self
        class(sym),              intent(in)    :: symop
        integer,                 intent(in)    :: kfromto(2)
        integer,                 intent(in)    :: cart_budget
        real(sp), allocatable :: tmp_dirs(:,:,:)
        real(sp) :: rsym(3,3), q(3)
        integer, allocatable :: target_shell_nodes(:)
        integer :: isym, ik, ifull, nfull, nkeep, pos, max_full
        call self%kill
        if( kfromto(1) > kfromto(2) ) THROW_HARD('invalid k range; shell_geom_new')
        if( cart_budget < 1 ) THROW_HARD('invalid Cartesian budget; shell_geom_new')
        self%kfromto      = kfromto
        self%nk           = kfromto(2) - kfromto(1) + 1
        self%nsym         = symop%get_nsym()
        self%cart_budget  = cart_budget
        self%target_nodes = max(cart_budget, self%nk)
        allocate(self%full_nodes(self%nk),   source=0)
        allocate(self%shell_nodes(self%nk),  source=0)
        allocate(self%offsets(self%nk+1),    source=1)
        allocate(self%theta_support(self%nk), source=0.d0)
        allocate(self%sym_rmat(3,3,self%nsym), source=0._sp)
        do isym = 1,self%nsym
            call symop%get_sym_rmat(isym, rsym)
            self%sym_rmat(:,:,isym) = rsym
        enddo
        call distribute_shell_budget(kfromto, self%target_nodes, target_shell_nodes)
        do ik = 1,self%nk
            self%full_nodes(ik) = max(12, target_shell_nodes(ik) * max(1,self%nsym))
            self%theta_support(ik) = sqrt(4.d0 * DPI / real(self%full_nodes(ik),dp))
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
        deallocate(tmp_dirs)
        deallocate(target_shell_nodes)
    end subroutine shell_geom_new

    subroutine shell_geom_kill( self )
        class(shell_field_geom), intent(inout) :: self
        if( allocated(self%full_nodes) ) deallocate(self%full_nodes)
        if( allocated(self%shell_nodes) ) deallocate(self%shell_nodes)
        if( allocated(self%offsets) ) deallocate(self%offsets)
        if( allocated(self%theta_support) ) deallocate(self%theta_support)
        if( allocated(self%node_dir) ) deallocate(self%node_dir)
        if( allocated(self%sym_rmat) ) deallocate(self%sym_rmat)
        self%nk           = 0
        self%kfromto      = 0
        self%nsym         = 1
        self%total_nodes  = 0
        self%target_nodes = 0
        self%cart_budget  = 0
    end subroutine shell_geom_kill

    subroutine shell_geom_log_summary( self, label )
        class(shell_field_geom), intent(in) :: self
        character(len=*),        intent(in) :: label
        real(dp) :: mem_mb, min_support, max_support, budget_ratio
        integer :: min_nodes, max_nodes
        if( .not. allocated(self%shell_nodes) ) return
        mem_mb = 0.d0
        if( allocated(self%node_dir) )then
            mem_mb = real(storage_size(self%node_dir),dp) * real(size(self%node_dir),dp) / (8.d0 * 1024.d0 * 1024.d0)
        endif
        min_nodes   = minval(self%shell_nodes)
        max_nodes   = maxval(self%shell_nodes)
        min_support = minval(self%theta_support)
        max_support = maxval(self%theta_support)
        budget_ratio = real(self%total_nodes,dp) / real(max(1,self%cart_budget),dp)
        write(logfhandle,'(A,1X,A,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,F8.3,2X,A,I0,2X,A,I0,2X,A,F8.3,2X,A,F8.3,2X,A,F10.3)') &
            &'obsfield shell-geom', trim(label), 'nk=', self%nk, 'kmin=', self%kfromto(1), 'kmax=', self%kfromto(2), &
            &'nsym=', self%nsym, 'cart_budget=', self%cart_budget, 'target_nodes=', self%target_nodes, &
            &'total_nodes=', self%total_nodes, 'budget_ratio=', budget_ratio, 'nodes_min=', min_nodes, 'nodes_max=', max_nodes, &
            &'support_min_deg=', rad2deg(real(min_support)), &
            &'support_max_deg=', rad2deg(real(max_support)), 'dir_mem_mb=', mem_mb
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
