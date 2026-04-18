!@descr: Bishop/Tipping probabilistic PCA with explicit isotropic noise
module simple_ppca
use simple_core_module_api
use simple_pca, only: pca
implicit none

public :: ppca
private

type, extends(pca) :: ppca
    private
    real, allocatable :: W(:,:)      !< loading matrix
    real, allocatable :: E_zn(:,:)   !< posterior latent means
    real, allocatable :: data(:,:)   !< reconstructed centered data
    real, allocatable :: Wt(:,:), M(:,:), Minv(:,:), MinvWt(:,:), S_xz(:,:), S_zz(:,:), Iq(:,:)
    real              :: sigma2 = 1.
    real(dp)          :: xnorm2_total = 0._dp
    logical           :: existence = .false.
  contains
    procedure :: new      => new_ppca
    procedure :: get_feat => get_feat_ppca
    procedure :: generate => generate_ppca
    procedure :: master   => master_ppca
    procedure, private :: init
    procedure, private :: em_opt
    procedure :: kill     => kill_ppca
end type

contains

    subroutine new_ppca( self, N, D, Q )
        class(ppca), intent(inout) :: self
        integer,     intent(in)    :: N, D, Q
        call self%kill
        self%N = N
        self%D = D
        self%Q = Q
        allocate(self%W(self%D,self%Q), self%E_zn(self%Q,self%N), self%data(self%D,self%N), &
                 self%Wt(self%Q,self%D), self%M(self%Q,self%Q), self%Minv(self%Q,self%Q), &
                 self%MinvWt(self%Q,self%D), self%S_xz(self%D,self%Q), self%S_zz(self%Q,self%Q), &
                 self%Iq(self%Q,self%Q), source=0.)
        self%existence = .true.
    end subroutine new_ppca

    function get_feat_ppca( self, i ) result( feat )
        class(ppca), intent(inout) :: self
        integer,     intent(in)    :: i
        real, allocatable          :: feat(:)
        allocate(feat(self%Q), source=self%E_zn(:,i))
    end function get_feat_ppca

    subroutine generate_ppca( self, i, avg, dat )
        class(ppca), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: avg(self%D)
        real,        intent(inout) :: dat(self%D)
        dat = avg + self%data(:,i)
    end subroutine generate_ppca

    subroutine master_ppca( self, pcavecs, maxpcaits )
        class(ppca),          intent(inout) :: self
        real,                 intent(in)    :: pcavecs(self%D,self%N)
        integer, optional,    intent(in)    :: maxpcaits
        integer, parameter :: MAX_ITS_DEFAULT = 50
        real(dp), parameter :: TOL_SIGMA = 1.e-6_dp
        real(dp), parameter :: TOL_WREL  = 1.e-5_dp
        integer :: it, maxits
        integer(int64) :: t0, t1
        real(real64)   :: trate
        real(dp) :: sigma_prev, w_change, sigma_rel
        if( present(maxpcaits) )then
            maxits = max(1, maxpcaits)
        else
            maxits = MAX_ITS_DEFAULT
        endif
        call system_clock(t0, trate)
        write(logfhandle,'(A,I8,A,I8,A,I8,A,I8)') 'PPCA start: N=', self%N, ' D=', self%D, ' Q=', self%Q, ' maxits=', maxits
        call flush(logfhandle)
        call self%init(pcavecs)
        sigma_prev = huge(1._dp)
        do it = 1, maxits
            call self%em_opt(pcavecs, w_change)
            sigma_rel = abs(real(self%sigma2,dp) - sigma_prev) / max(1.e-8_dp, sigma_prev)
            if( it == 1 .or. mod(it,5) == 0 .or. it == maxits )then
                write(logfhandle,'(A,I4,A,ES10.3,A,ES10.3,A,ES10.3)') &
                    'PPCA iter=', it, ' sigma2=', real(self%sigma2,dp), ' dW=', w_change, ' dsigma=', sigma_rel
                call flush(logfhandle)
            endif
            if( sigma_rel <= TOL_SIGMA .and. w_change <= TOL_WREL )then
                write(logfhandle,'(A,I4,A)') 'PPCA converged at iter=', it, ''
                call flush(logfhandle)
                exit
            endif
            sigma_prev = real(self%sigma2,dp)
        enddo
        self%data = matmul(self%W, self%E_zn)
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A,ES10.3)') 'PPCA total: ', real(t1-t0)/real(trate), ' s; sigma2=', real(self%sigma2,dp)
        call flush(logfhandle)
    end subroutine master_ppca

    subroutine init( self, pcavecs )
        use simple_rnd, only: ran3
        class(ppca), intent(inout) :: self
        real,        intent(in)    :: pcavecs(self%D,self%N)
        integer :: i, j
        real(dp) :: mean_sq
        self%Iq = 0.
        do i = 1,self%Q
            self%Iq(i,i) = 1.
        enddo
        do i = 1,self%D
            do j = 1,self%Q
                self%W(i,j) = 0.01 * ran3()
            enddo
        enddo
        mean_sq = sum(real(pcavecs,dp)**2) / real(max(1,self%D*self%N), dp)
        self%xnorm2_total = sum(real(pcavecs,dp)**2)
        self%sigma2 = real(max(mean_sq * 0.1_dp, real(DTINY,dp)))
        self%E_zn   = 0.
        self%data   = 0.
    end subroutine init

    subroutine em_opt( self, pcavecs, w_change )
        use simple_math, only: matinv
        class(ppca), intent(inout) :: self
        real,        intent(in)    :: pcavecs(self%D,self%N)
        real(dp),    intent(out)   :: w_change
        integer, parameter :: BLOCKSZ = 512
        integer :: i, ibeg, iend, nblk, err
        integer(int64) :: t0, t1
        real(real64)   :: trate
        real(dp) :: sigma_floor, term2, sigma_num
        real, allocatable :: w_old(:,:), wt_w(:,:)
        sigma_floor = max(real(DTINY,dp), 1.e-8_dp)
        allocate(w_old(self%D,self%Q), wt_w(self%Q,self%Q), source=0.)
        w_old = self%W
        self%Wt = transpose(self%W)
        self%M  = matmul(self%Wt, self%W) + self%sigma2 * self%Iq
        call matinv(self%M, self%Minv, self%Q, err)
        if( err == -1 )then
            self%Minv = 0.
            do i = 1,self%Q
                self%Minv(i,i) = 1.
            enddo
        endif
        self%MinvWt = matmul(self%Minv, self%Wt)
        self%S_xz = 0.
        self%S_zz = 0.
        nblk = max(1, (self%N + BLOCKSZ - 1) / BLOCKSZ)
        call system_clock(t0, trate)
        write(logfhandle,'(A,I8,A)') 'PPCA blocked E-step: blocks=', nblk, ''
        call flush(logfhandle)
        !$omp parallel default(shared) private(i,ibeg,iend)
        block
            real, allocatable :: ez_blk(:,:), sxz_blk(:,:), szz_blk(:,:)
            allocate(ez_blk(self%Q,BLOCKSZ), sxz_blk(self%D,self%Q), szz_blk(self%Q,self%Q), source=0.)
            !$omp do schedule(dynamic)
            do i = 1, self%N, BLOCKSZ
                ibeg = i
                iend = min(self%N, i + BLOCKSZ - 1)
                ez_blk(:,1:iend-ibeg+1) = matmul(self%MinvWt, pcavecs(:,ibeg:iend))
                self%E_zn(:,ibeg:iend) = ez_blk(:,1:iend-ibeg+1)
                sxz_blk = matmul(pcavecs(:,ibeg:iend), transpose(ez_blk(:,1:iend-ibeg+1)))
                szz_blk = matmul(ez_blk(:,1:iend-ibeg+1), transpose(ez_blk(:,1:iend-ibeg+1)))
                !$omp critical(ppca_accum)
                self%S_xz = self%S_xz + sxz_blk
                self%S_zz = self%S_zz + szz_blk
                !$omp end critical(ppca_accum)
            enddo
            !$omp end do
            deallocate(ez_blk, sxz_blk, szz_blk)
        end block
        !$omp end parallel
        self%S_zz = self%S_zz + real(self%N,kind(self%sigma2)) * self%sigma2 * self%Minv
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A)') 'PPCA blocked E-step total: ', real(t1-t0)/real(trate), ' s'
        call flush(logfhandle)
        call matinv(self%S_zz, self%M, self%Q, err)
        if( err == -1 )then
            self%M = 0.
            do i = 1,self%Q
                self%M(i,i) = 1.
            enddo
        endif
        self%W = matmul(self%S_xz, self%M)
        self%Wt = transpose(self%W)
        wt_w = matmul(self%Wt, self%W)
        call system_clock(t0)
        term2 = real(self%N,dp) * real(self%sigma2,dp) * sum(real(self%Minv,dp) * transpose(real(wt_w,dp))) + &
            sum(real(self%E_zn,dp) * real(matmul(wt_w, self%E_zn),dp))
        sigma_num = self%xnorm2_total - 2._dp * sum(real(self%E_zn,dp) * real(matmul(self%Wt, pcavecs),dp)) + term2
        self%sigma2 = real(max(sigma_num / real(max(1,self%N*self%D),dp), sigma_floor))
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A)') 'PPCA sigma-step total: ', real(t1-t0)/real(trate), ' s'
        call flush(logfhandle)
        w_change = sqrt(sum((real(self%W,dp) - real(w_old,dp))**2) / real(max(1,self%D*self%Q),dp))
        deallocate(w_old, wt_w)
    end subroutine em_opt

    subroutine kill_ppca( self )
        class(ppca), intent(inout) :: self
        if( self%existence )then
            deallocate(self%W, self%E_zn, self%data, self%Wt, self%M, self%Minv, self%MinvWt, self%S_xz, self%S_zz, self%Iq)
            self%existence = .false.
        endif
    end subroutine kill_ppca

end module simple_ppca
