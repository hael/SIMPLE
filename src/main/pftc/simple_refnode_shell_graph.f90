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
    real    :: support_mult = 1.5
    real(dp):: theta_support = 0.d0
    real(sp), allocatable :: node_dir(:,:,:) ! (3,node,shell)
    real(sp), allocatable :: sym_rmat(:,:,:) ! (3,3,nsym)
  contains
    procedure :: new       => refnode_graph_new
    procedure :: kill      => refnode_graph_kill
    procedure :: log_stats => refnode_graph_log_stats
end type refnode_shell_graph

contains

    subroutine refnode_graph_new( self, reforis, symop, polar_x, polar_y, kfromto, nprojs, support_mult )
        class(refnode_shell_graph), intent(inout) :: self
        class(oris),                intent(in)    :: reforis
        class(sym),                 intent(in)    :: symop
        real(sp),                   intent(in)    :: polar_x(:,:), polar_y(:,:)
        integer,                    intent(in)    :: kfromto(2), nprojs
        real,                       intent(in)    :: support_mult
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
        do isym = 1,self%nsym
            call symop%get_sym_rmat(isym, rsym)
            self%sym_rmat(:,:,isym) = rsym
        enddo
        do ik = 1,self%nk
            do iproj = 1,self%nprojs
                rmat = reforis%get_mat(iproj)
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
    end subroutine refnode_graph_new

    subroutine refnode_graph_kill( self )
        class(refnode_shell_graph), intent(inout) :: self
        if( allocated(self%node_dir) ) deallocate(self%node_dir)
        if( allocated(self%sym_rmat) ) deallocate(self%sym_rmat)
        self%pftsz         = 0
        self%nprojs        = 0
        self%nk            = 0
        self%kfromto       = 0
        self%nsym          = 1
        self%nodes_shell   = 0
        self%support_mult  = 1.5
        self%theta_support = 0.d0
    end subroutine refnode_graph_kill

    subroutine refnode_graph_log_stats( self, label )
        class(refnode_shell_graph), intent(in) :: self
        character(len=*),           intent(in) :: label
        real(dp) :: mem_mb
        integer(kind=8) :: total_nodes
        integer :: ik
        if( .not. allocated(self%node_dir) ) return
        total_nodes = int(self%nodes_shell,kind=8) * int(self%nk,kind=8)
        mem_mb = real(storage_size(self%node_dir),dp) * real(size(self%node_dir),dp) / (8.d0 * 1024.d0 * 1024.d0)
        write(logfhandle,'(A,1X,A,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,F8.3,2X,A,F8.3,2X,A,F10.3)') &
            &'obsfield refnodes graph', trim(label), &
            &'nprojs=', self%nprojs, 'pftsz=', self%pftsz, 'nk=', self%nk, 'nsym=', self%nsym, &
            &'nodes_per_shell=', self%nodes_shell, 'support_mult=', self%support_mult, &
            &'theta_support_deg=', rad2deg(real(self%theta_support)), 'dir_mem_mb=', mem_mb
        do ik = 1,self%nk
            call log_shell_stats(self, ik, self%kfromto(1)+ik-1)
        enddo
        write(logfhandle,'(A,1X,A,2X,A,I0)') 'obsfield refnodes graph', trim(label), 'total_nodes=', total_nodes
    end subroutine refnode_graph_log_stats

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
            accepted_sumsq = accepted_sumsq + real(accepted*accepted,dp)
        enddo
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
