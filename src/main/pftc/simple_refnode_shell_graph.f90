!@descr: development scaffold for exact reference-node shell geometry
module simple_refnode_shell_graph
use simple_pftc_api
implicit none

public :: refnode_shell_graph
private
#include "simple_local_flags.inc"

type :: refnode_shell_graph
    private
    integer :: pftsz       = 0
    integer :: nprojs      = 0
    integer :: nk          = 0
    integer :: kfromto(2)  = 0
    integer :: nsym        = 1
    integer :: nodes_shell = 0
    integer :: ntheta      = 0
    integer :: nphi        = 0
    integer :: ncells      = 0
    real    :: support_mult = 1.5
    real(dp):: theta_support = 0.d0
    real(dp):: dtheta      = 0.d0
    real(dp):: dphi        = 0.d0
    real(sp), allocatable :: node_dir(:,:,:) ! (3,node,shell)
    real(sp), allocatable :: sym_rmat(:,:,:) ! (3,3,nsym)
    integer,  allocatable :: cell_head(:), cell_next(:), node_cell(:)
  contains
    procedure :: new       => refnode_graph_new
    procedure :: kill      => refnode_graph_kill
    procedure :: log_summary => refnode_graph_log_summary
    procedure :: log_stats => refnode_graph_log_stats
    procedure :: query_nodes => refnode_graph_query_nodes
    procedure :: node_to_polar => refnode_graph_node_to_polar
end type refnode_shell_graph

contains

    subroutine refnode_graph_new( self, reforis, symop, polar_x, polar_y, kfromto, nprojs, support_mult )
        class(refnode_shell_graph), intent(inout) :: self
        class(oris),                intent(in)    :: reforis
        class(sym),                 intent(in)    :: symop
        real(sp),                   intent(in)    :: polar_x(:,:), polar_y(:,:)
        integer,                    intent(in)    :: kfromto(2), nprojs
        real,                       intent(in)    :: support_mult
        real(sp), allocatable :: ref_rmat(:,:,:)
        real(sp) :: rmat(3,3), q(3), qnorm, rsym(3,3)
        integer  :: ik, iproj, irot, inode, isym
        call self%kill
        if( nprojs < 1 )then
            THROW_HARD('invalid nprojs; refnode_graph_new')
        endif
        if( kfromto(1) > kfromto(2) )then
            THROW_HARD('invalid k range; refnode_graph_new')
        endif
        if( size(polar_x,1) /= size(polar_y,1) .or. size(polar_x,2) /= size(polar_y,2) )then
            THROW_HARD('polar coordinate array mismatch; refnode_graph_new')
        endif
        if( size(polar_x,2) /= kfromto(2)-kfromto(1)+1 )then
            THROW_HARD('polar coordinate k range mismatch; refnode_graph_new')
        endif
        if( reforis%get_noris() < nprojs )then
            THROW_HARD('not enough reference orientations; refnode_graph_new')
        endif
        self%pftsz        = size(polar_x,1)
        self%nprojs       = nprojs
        self%nk           = kfromto(2) - kfromto(1) + 1
        self%kfromto      = kfromto
        self%nsym         = symop%get_nsym()
        self%nodes_shell  = self%pftsz * self%nprojs
        self%support_mult = support_mult
        self%theta_support = real(support_mult,dp) * sqrt(4.d0 * DPI / real(max(1,self%nodes_shell),dp))
        self%theta_support = min(DPI, max(0.d0, self%theta_support))
        allocate(self%node_dir(3,self%nodes_shell,self%nk), source=0._sp)
        allocate(self%sym_rmat(3,3,self%nsym), source=0._sp)
        allocate(ref_rmat(3,3,self%nprojs), source=0._sp)
        do isym = 1,self%nsym
            call symop%get_sym_rmat(isym, rsym)
            self%sym_rmat(:,:,isym) = rsym
        enddo
        do iproj = 1,self%nprojs
            ref_rmat(:,:,iproj) = reforis%get_mat(iproj)
        enddo
        !$omp parallel do collapse(2) default(shared) private(ik,iproj,irot,inode,rmat,q,qnorm) schedule(static) proc_bind(close)
        do ik = 1,self%nk
            do iproj = 1,self%nprojs
                rmat = ref_rmat(:,:,iproj)
                do irot = 1,self%pftsz
                    inode = (iproj - 1) * self%pftsz + irot
                    q(1) = polar_x(irot,ik) * rmat(1,1) + polar_y(irot,ik) * rmat(2,1)
                    q(2) = polar_x(irot,ik) * rmat(1,2) + polar_y(irot,ik) * rmat(2,2)
                    q(3) = polar_x(irot,ik) * rmat(1,3) + polar_y(irot,ik) * rmat(2,3)
                    qnorm = sqrt(sum(q*q))
                    if( qnorm > TINY ) self%node_dir(:,inode,ik) = q / qnorm
                enddo
            enddo
        enddo
        !$omp end parallel do
        deallocate(ref_rmat)
        call build_index(self)
    end subroutine refnode_graph_new

    subroutine refnode_graph_kill( self )
        class(refnode_shell_graph), intent(inout) :: self
        if( allocated(self%node_dir) ) deallocate(self%node_dir)
        if( allocated(self%sym_rmat) ) deallocate(self%sym_rmat)
        if( allocated(self%cell_head) ) deallocate(self%cell_head)
        if( allocated(self%cell_next) ) deallocate(self%cell_next)
        if( allocated(self%node_cell) ) deallocate(self%node_cell)
        self%pftsz         = 0
        self%nprojs        = 0
        self%nk            = 0
        self%kfromto       = 0
        self%nsym          = 1
        self%nodes_shell   = 0
        self%ntheta        = 0
        self%nphi          = 0
        self%ncells        = 0
        self%support_mult  = 1.5
        self%theta_support = 0.d0
        self%dtheta        = 0.d0
        self%dphi          = 0.d0
    end subroutine refnode_graph_kill

    subroutine refnode_graph_log_summary( self, label )
        class(refnode_shell_graph), intent(in) :: self
        character(len=*),           intent(in) :: label
        real(dp) :: mem_mb
        integer(kind=8) :: total_nodes
        if( .not. allocated(self%node_dir) ) return
        total_nodes = int(self%nodes_shell,kind=8) * int(self%nk,kind=8)
        mem_mb = real(storage_size(self%node_dir),dp) * real(size(self%node_dir),dp) / (8.d0 * 1024.d0 * 1024.d0)
        write(logfhandle,'(A,1X,A,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,F8.3,2X,A,F8.3,2X,A,F10.3)') &
            &'obsfield refnodes graph', trim(label), &
            &'nprojs=', self%nprojs, 'pftsz=', self%pftsz, 'nk=', self%nk, 'nsym=', self%nsym, &
            &'nodes_per_shell=', self%nodes_shell, 'ntheta=', self%ntheta, 'nphi=', self%nphi, &
            &'support_mult=', self%support_mult, &
            &'theta_support_deg=', rad2deg(real(self%theta_support)), 'dir_mem_mb=', mem_mb
        write(logfhandle,'(A,1X,A,2X,A,I0)') 'obsfield refnodes graph', trim(label), 'total_nodes=', total_nodes
    end subroutine refnode_graph_log_summary

    subroutine refnode_graph_log_stats( self, label )
        class(refnode_shell_graph), intent(in) :: self
        character(len=*),           intent(in) :: label
        integer :: ik
        call self%log_summary(label)
        if( .not. allocated(self%node_dir) ) return
        call log_lookup_stats(self)
        do ik = 1,self%nk
            call log_shell_stats(self, ik, self%kfromto(1)+ik-1)
        enddo
    end subroutine refnode_graph_log_stats

    subroutine build_index( self )
        class(refnode_shell_graph), intent(inout) :: self
        integer :: inode, icell
        if( .not. allocated(self%node_dir) ) return
        self%ntheta = max(4, ceiling(DPI   / max(self%theta_support, 1.d-6)))
        self%nphi   = max(8, ceiling(DTWOPI / max(self%theta_support, 1.d-6)))
        self%ncells = self%ntheta * self%nphi
        self%dtheta = DPI / real(self%ntheta,dp)
        self%dphi   = DTWOPI / real(self%nphi,dp)
        allocate(self%cell_head(self%ncells), source=0)
        allocate(self%cell_next(self%nodes_shell), source=0)
        allocate(self%node_cell(self%nodes_shell), source=0)
        do inode = 1,self%nodes_shell
            icell = cell_for_vec(self, self%node_dir(:,inode,1))
            self%node_cell(inode) = icell
            self%cell_next(inode) = self%cell_head(icell)
            self%cell_head(icell) = inode
        enddo
    end subroutine build_index

    integer function cell_for_vec( self, vec ) result( icell )
        class(refnode_shell_graph), intent(in) :: self
        real(sp),                    intent(in) :: vec(3)
        real(dp) :: theta, phi
        integer  :: itheta, iphi
        call vec_angles(vec, theta, phi)
        itheta = min(self%ntheta, max(1, int(theta / self%dtheta) + 1))
        iphi   = min(self%nphi,   max(1, int(phi   / self%dphi)   + 1))
        icell  = (itheta - 1) * self%nphi + iphi
    end function cell_for_vec

    subroutine vec_angles( vec, theta, phi )
        real(sp), intent(in)  :: vec(3)
        real(dp), intent(out) :: theta, phi
        real(dp) :: x, y, z, normv
        x = real(vec(1),dp)
        y = real(vec(2),dp)
        z = real(vec(3),dp)
        normv = sqrt(x*x + y*y + z*z)
        if( normv <= DTINY )then
            theta = 0.d0
            phi   = 0.d0
        else
            z     = min(1.d0, max(-1.d0, z / normv))
            theta = acos(z)
            phi   = atan2(y, x)
            if( phi < 0.d0 ) phi = phi + DTWOPI
        endif
    end subroutine vec_angles

    subroutine query_index( self, q, ik, candidates, accepted )
        class(refnode_shell_graph), intent(in)  :: self
        real(sp),                    intent(in)  :: q(3)
        integer,                     intent(in)  :: ik
        integer,                     intent(out) :: candidates, accepted
        real(dp) :: theta, phi, theta_lo, theta_hi, sin_min, dphi_max, cos_support, dotv
        integer  :: itheta, iphi, iphi0, iphi1, it0, it1, inode, icell, iphi_wrap
        logical  :: full_phi
        candidates = 0
        accepted   = 0
        if( .not. allocated(self%cell_head) ) return
        call vec_angles(q, theta, phi)
        cos_support = cos(self%theta_support)
        it0 = max(1, int(max(0.d0, theta - self%theta_support) / self%dtheta) + 1)
        it1 = min(self%ntheta, int(min(DPI, theta + self%theta_support) / self%dtheta) + 1)
        do itheta = it0,it1
            theta_lo = real(itheta - 1,dp) * self%dtheta
            theta_hi = real(itheta,dp) * self%dtheta
            sin_min  = min(sin(theta), min(sin(theta_lo), sin(theta_hi)))
            full_phi = sin_min <= sin(0.5d0 * self%theta_support)
            if( full_phi )then
                iphi0 = 1
                iphi1 = self%nphi
            else
                dphi_max = 2.d0 * asin(min(1.d0, sin(0.5d0 * self%theta_support) / sqrt(max(DTINY, sin(theta) * sin_min))))
                dphi_max = min(DPI, dphi_max + self%dphi)
                if( dphi_max >= DPI )then
                    full_phi = .true.
                    iphi0 = 1
                    iphi1 = self%nphi
                else
                    iphi0 = int((phi - dphi_max) / self%dphi)
                    iphi1 = int((phi + dphi_max) / self%dphi) + 2
                endif
            endif
            do iphi = iphi0,iphi1
                if( full_phi )then
                    iphi_wrap = iphi
                else
                    iphi_wrap = modulo(iphi - 1, self%nphi) + 1
                endif
                icell = (itheta - 1) * self%nphi + iphi_wrap
                inode = self%cell_head(icell)
                do while( inode > 0 )
                    candidates = candidates + 1
                    dotv = sum(real(q,dp) * real(self%node_dir(:,inode,ik),dp))
                    if( dotv >= cos_support ) accepted = accepted + 1
                    inode = self%cell_next(inode)
                enddo
            enddo
        enddo
    end subroutine query_index

    subroutine refnode_graph_query_nodes( self, q, shell, node_ids, weights, nfound, candidates, overflow )
        class(refnode_shell_graph), intent(in)  :: self
        real(sp),                    intent(in)  :: q(3)
        integer,                     intent(in)  :: shell
        integer,                     intent(out) :: node_ids(:)
        real(dp),                    intent(out) :: weights(:)
        integer,                     intent(out) :: nfound
        integer,           optional, intent(out) :: candidates
        logical,           optional, intent(out) :: overflow
        real(dp) :: theta, phi, theta_lo, theta_hi, sin_min, dphi_max, cos_support, dotv, support_den
        real(sp) :: qdir(3)
        real     :: qnorm
        integer  :: itheta, iphi, iphi0, iphi1, it0, it1, inode, icell, iphi_wrap, ncand
        logical  :: full_phi, l_overflow
        nfound     = 0
        ncand      = 0
        l_overflow = .false.
        if( present(candidates) ) candidates = 0
        if( present(overflow) ) overflow = .false.
        if( .not. allocated(self%cell_head) ) return
        if( shell < self%kfromto(1) .or. shell > self%kfromto(2) ) return
        qnorm = sqrt(sum(q*q))
        if( qnorm <= TINY ) return
        qdir = q / qnorm
        call vec_angles(qdir, theta, phi)
        cos_support = cos(self%theta_support)
        support_den = max(DTINY, 1.d0 - cos_support)
        it0 = max(1, int(max(0.d0, theta - self%theta_support) / self%dtheta) + 1)
        it1 = min(self%ntheta, int(min(DPI, theta + self%theta_support) / self%dtheta) + 1)
        do itheta = it0,it1
            theta_lo = real(itheta - 1,dp) * self%dtheta
            theta_hi = real(itheta,dp) * self%dtheta
            sin_min  = min(sin(theta), min(sin(theta_lo), sin(theta_hi)))
            full_phi = sin_min <= sin(0.5d0 * self%theta_support)
            if( full_phi )then
                iphi0 = 1
                iphi1 = self%nphi
            else
                dphi_max = 2.d0 * asin(min(1.d0, sin(0.5d0 * self%theta_support) / sqrt(max(DTINY, sin(theta) * sin_min))))
                dphi_max = min(DPI, dphi_max + self%dphi)
                if( dphi_max >= DPI )then
                    full_phi = .true.
                    iphi0 = 1
                    iphi1 = self%nphi
                else
                    iphi0 = int((phi - dphi_max) / self%dphi)
                    iphi1 = int((phi + dphi_max) / self%dphi) + 2
                endif
            endif
            do iphi = iphi0,iphi1
                if( full_phi )then
                    iphi_wrap = iphi
                else
                    iphi_wrap = modulo(iphi - 1, self%nphi) + 1
                endif
                icell = (itheta - 1) * self%nphi + iphi_wrap
                inode = self%cell_head(icell)
                do while( inode > 0 )
                    ncand = ncand + 1
                    dotv = sum(real(qdir,dp) * real(self%node_dir(:,inode,shell-self%kfromto(1)+1),dp))
                    if( dotv >= cos_support )then
                        if( nfound < size(node_ids) .and. nfound < size(weights) )then
                            nfound = nfound + 1
                            node_ids(nfound) = inode
                            weights(nfound)  = max(0.d0, (dotv - cos_support) / support_den)
                        else
                            l_overflow = .true.
                        endif
                    endif
                    inode = self%cell_next(inode)
                enddo
            enddo
        enddo
        if( present(candidates) ) candidates = ncand
        if( present(overflow) ) overflow = l_overflow
    end subroutine refnode_graph_query_nodes

    pure subroutine refnode_graph_node_to_polar( self, inode, iproj, irot )
        class(refnode_shell_graph), intent(in)  :: self
        integer,                    intent(in)  :: inode
        integer,                    intent(out) :: iproj, irot
        if( inode < 1 .or. inode > self%nodes_shell .or. self%pftsz < 1 )then
            iproj = 0
            irot  = 0
        else
            iproj = (inode - 1) / self%pftsz + 1
            irot  = mod(inode - 1, self%pftsz) + 1
        endif
    end subroutine refnode_graph_node_to_polar

    subroutine log_lookup_stats( self )
        class(refnode_shell_graph), intent(in) :: self
        real(dp) :: cand_mean, acc_mean, brute_mean, weight_mean, weight_sumsq, cand_sumsq, acc_sumsq
        real(dp) :: cand_sdev, acc_sdev, weight_sdev, dotv, cos_support, wsum
        real(sp) :: q(3)
        integer, allocatable :: node_ids(:)
        real(dp), allocatable :: weights(:)
        integer  :: stride, isample, inode, nsamples, candidates, accepted, brute, nfound
        integer  :: cand_max, acc_max, empty_queries, mismatches, overflows
        logical  :: overflow
        if( .not. allocated(self%cell_head) ) return
        allocate(node_ids(256), source=0)
        allocate(weights(256), source=0.d0)
        stride        = max(1, self%nodes_shell / 256)
        nsamples      = 0
        cand_mean     = 0.d0
        acc_mean      = 0.d0
        brute_mean    = 0.d0
        weight_mean   = 0.d0
        cand_sumsq    = 0.d0
        acc_sumsq     = 0.d0
        weight_sumsq  = 0.d0
        cand_max      = 0
        acc_max       = 0
        empty_queries = 0
        mismatches    = 0
        overflows     = 0
        cos_support   = cos(self%theta_support)
        do isample = 1,self%nodes_shell,stride
            q = self%node_dir(:,isample,1)
            call query_index(self, q, 1, candidates, accepted)
            call self%query_nodes(q, self%kfromto(1), node_ids, weights, nfound, overflow=overflow)
            brute = 0
            do inode = 1,self%nodes_shell
                dotv = sum(real(q,dp) * real(self%node_dir(:,inode,1),dp))
                if( dotv >= cos_support ) brute = brute + 1
            enddo
            nsamples   = nsamples + 1
            cand_mean  = cand_mean  + real(candidates,dp)
            acc_mean   = acc_mean   + real(accepted,dp)
            brute_mean = brute_mean + real(brute,dp)
            wsum       = 0.d0
            if( nfound > 0 ) wsum = sum(weights(:nfound))
            weight_mean  = weight_mean  + wsum
            cand_sumsq = cand_sumsq + real(candidates,dp) * real(candidates,dp)
            acc_sumsq  = acc_sumsq  + real(accepted,dp)   * real(accepted,dp)
            weight_sumsq = weight_sumsq + wsum * wsum
            cand_max   = max(cand_max, candidates)
            acc_max    = max(acc_max, accepted)
            if( accepted == 0 ) empty_queries = empty_queries + 1
            if( accepted /= brute ) mismatches = mismatches + 1
            if( overflow ) overflows = overflows + 1
        enddo
        if( nsamples > 0 )then
            cand_mean  = cand_mean  / real(nsamples,dp)
            acc_mean   = acc_mean   / real(nsamples,dp)
            brute_mean = brute_mean / real(nsamples,dp)
            weight_mean = weight_mean / real(nsamples,dp)
            cand_sumsq = cand_sumsq / real(nsamples,dp)
            acc_sumsq  = acc_sumsq  / real(nsamples,dp)
            weight_sumsq = weight_sumsq / real(nsamples,dp)
            cand_sdev  = sqrt(max(0.d0, cand_sumsq - cand_mean*cand_mean))
            acc_sdev   = sqrt(max(0.d0, acc_sumsq  - acc_mean*acc_mean))
            weight_sdev = sqrt(max(0.d0, weight_sumsq - weight_mean*weight_mean))
        else
            cand_sdev = 0.d0
            acc_sdev  = 0.d0
            weight_sdev = 0.d0
        endif
        write(logfhandle,'(A,2X,A,I0,2X,A,F10.3,2X,A,F10.3,2X,A,I0,2X,A,F10.3,2X,A,F10.3,2X,A,I0,2X,A,F10.3,2X,A,F10.3,2X,A,F10.3,2X,A,I0,2X,A,I0,2X,A,I0)') &
            &'obsfield refnodes lookup', 'sampled=', nsamples, 'candidates_mean=', cand_mean, &
            &'candidates_sdev=', cand_sdev, 'candidates_max=', cand_max, 'accepted_mean=', acc_mean, &
            &'accepted_sdev=', acc_sdev, 'accepted_max=', acc_max, 'brute_mean=', brute_mean, &
            &'weight_sum_mean=', weight_mean, 'weight_sum_sdev=', weight_sdev, &
            &'empty_queries=', empty_queries, 'lookup_mismatch=', mismatches, 'query_overflows=', overflows
        deallocate(node_ids, weights)
    end subroutine log_lookup_stats

    subroutine log_shell_stats( self, ik, shell )
        class(refnode_shell_graph), intent(in) :: self
        integer,                    intent(in) :: ik, shell
        real(dp) :: cos_support, dotv, accepted_mean, accepted_sumsq, accepted_sdev
        integer  :: stride, isample, inode, nsamples, accepted, acc_min, acc_max
        cos_support = cos(self%theta_support)
        stride      = max(1, self%nodes_shell / 256)
        nsamples    = 0
        acc_min     = huge(acc_min)
        acc_max     = 0
        accepted_mean  = 0.d0
        accepted_sumsq = 0.d0
        !$omp parallel do default(shared) private(isample,inode,accepted,dotv) &
        !$omp reduction(+:nsamples,accepted_mean,accepted_sumsq) reduction(min:acc_min) reduction(max:acc_max) &
        !$omp schedule(static) proc_bind(close)
        do isample = 1,self%nodes_shell,stride
            accepted = 0
            do inode = 1,self%nodes_shell
                dotv = sum(real(self%node_dir(:,isample,ik),dp) * real(self%node_dir(:,inode,ik),dp))
                if( dotv >= cos_support ) accepted = accepted + 1
            enddo
            nsamples       = nsamples + 1
            acc_min        = min(acc_min, accepted)
            acc_max        = max(acc_max, accepted)
            accepted_mean  = accepted_mean  + real(accepted,dp)
            accepted_sumsq = accepted_sumsq + real(accepted,dp) * real(accepted,dp)
        enddo
        !$omp end parallel do
        if( nsamples > 0 )then
            accepted_mean  = accepted_mean / real(nsamples,dp)
            accepted_sumsq = accepted_sumsq / real(nsamples,dp)
            accepted_sdev  = sqrt(max(0.d0, accepted_sumsq - accepted_mean*accepted_mean))
        else
            acc_min        = 0
            accepted_mean  = 0.d0
            accepted_sdev  = 0.d0
        endif
        write(logfhandle,'(A,I4,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,F10.3,2X,A,F10.3)') &
            &'obsfield refnodes shell k=', shell, 'nodes=', self%nodes_shell, 'sampled=', nsamples, &
            &'accepted_min=', acc_min, 'accepted_max=', acc_max, 'accepted_mean=', accepted_mean, &
            &'accepted_sdev=', accepted_sdev
    end subroutine log_shell_stats

end module simple_refnode_shell_graph
