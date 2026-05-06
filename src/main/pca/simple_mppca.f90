!@descr: Mixture of probabilistic PCA analyzers with soft responsibilities
module simple_mppca
use simple_core_module_api
use simple_pca,  only: pca
use simple_ppca, only: ppca
implicit none

public :: mppca
private
#include "simple_local_flags.inc"

integer, parameter :: MPPCA_DEFAULT_K = 4

type, extends(pca) :: mppca
    private
    integer           :: K = MPPCA_DEFAULT_K
    integer           :: nthr = 1
    real, allocatable :: pi(:)             !< mixture weights
    real, allocatable :: sigma2(:)         !< isotropic noise per component
    real, allocatable :: mu(:,:)           !< component means, D x K
    real, allocatable :: W(:,:,:)          !< loading matrices, D x Q x K
    real, allocatable :: E_zn(:,:,:)       !< posterior latent means, Q x N x K
    real, allocatable :: resp(:,:)         !< responsibilities, K x N
    real, allocatable :: data(:,:)         !< reconstructed centered data, D x N
    real, allocatable :: feat_mix(:,:)     !< responsibility-weighted latent means, Q x N
    real, allocatable :: eigvals(:)        !< responsibility-weighted signal spectrum
    real, allocatable :: Iq(:,:)           !< latent identity
    character(len=16) :: recon_mode = 'soft' !< reconstruction mode: soft responsibility blend or hard max-responsibility chart
    logical           :: existence = .false.
contains
    procedure :: new        => new_mppca
    procedure :: set_params => set_params_mppca
    procedure :: get_feat   => get_feat_mppca
    procedure :: get_eigvals => get_eigvals_mppca
    procedure :: generate   => generate_mppca
    procedure :: generate_from_component => generate_from_component_mppca
    procedure :: generate_from_responsibilities => generate_from_responsibilities_mppca
    procedure :: reconstruct_external => reconstruct_external_mppca
    procedure :: master     => master_mppca
    procedure :: kill       => kill_mppca
    procedure, private :: init_mixture
    procedure, private :: init_from_ppca
    procedure, private :: em_step
    procedure, private :: reconstruct
    procedure, private :: calc_eigvals
    procedure, private :: infer_one
end type mppca

contains

    subroutine new_mppca( self, N, D, Q )
        class(mppca), intent(inout) :: self
        integer,      intent(in)    :: N, D, Q
        call self%kill
        self%N = N
        self%D = D
        self%Q = Q
        self%K = MPPCA_DEFAULT_K
        self%nthr = 1
        call self%init_mixture()
    end subroutine new_mppca

    subroutine set_params_mppca( self, K, nthr, recon_mode )
        class(mppca),      intent(inout) :: self
        integer, optional, intent(in)    :: K, nthr
        character(len=*), optional, intent(in) :: recon_mode
        integer :: knew
        knew = self%K
        if( present(K) ) knew = max(1, K)
        if( present(nthr) ) self%nthr = max(1, nthr)
        if( present(recon_mode) )then
            select case(trim(recon_mode))
                case('soft', 'hard')
                    self%recon_mode = trim(recon_mode)
                case DEFAULT
                    THROW_HARD('mppca_recon must be soft or hard')
            end select
        endif
        if( knew /= self%K )then
            self%K = knew
            if( self%existence ) call self%init_mixture()
        endif
    end subroutine set_params_mppca

    function get_feat_mppca( self, i ) result( feat )
        class(mppca), intent(inout) :: self
        integer,      intent(in)    :: i
        real, allocatable :: feat(:)
        allocate(feat(self%Q), source=self%feat_mix(:,i))
    end function get_feat_mppca

    function get_eigvals_mppca( self ) result( eigvals )
        class(mppca), intent(in) :: self
        real, allocatable :: eigvals(:)
        if( allocated(self%eigvals) )then
            allocate(eigvals(size(self%eigvals)), source=self%eigvals)
        else
            allocate(eigvals(0))
        endif
    end function get_eigvals_mppca

    subroutine generate_mppca( self, i, avg, dat )
        class(mppca), intent(inout) :: self
        integer,      intent(in)    :: i
        real,         intent(in)    :: avg(self%D)
        real,         intent(inout) :: dat(self%D)
        dat = avg + self%data(:,i)
    end subroutine generate_mppca

    subroutine generate_from_component_mppca( self, k, z, avg, dat )
        class(mppca), intent(in)    :: self
        integer,      intent(in)    :: k
        real,         intent(in)    :: z(self%Q), avg(self%D)
        real,         intent(inout) :: dat(self%D)
        integer :: kk
        kk = max(1, min(self%K, k))
        dat = avg + self%mu(:,kk) + matmul(self%W(:,:,kk), z)
    end subroutine generate_from_component_mppca

    subroutine generate_from_responsibilities_mppca( self, resp_in, z_by_comp, avg, dat )
        class(mppca), intent(in)    :: self
        real,         intent(in)    :: resp_in(self%K), z_by_comp(self%Q,self%K), avg(self%D)
        real,         intent(inout) :: dat(self%D)
        integer :: k
        dat = avg
        do k = 1,self%K
            dat = dat + resp_in(k) * (self%mu(:,k) + matmul(self%W(:,:,k), z_by_comp(:,k)))
        enddo
    end subroutine generate_from_responsibilities_mppca

    subroutine reconstruct_external_mppca( self, centered_in, centered_out )
        class(mppca), intent(inout) :: self
        real,         intent(in)    :: centered_in(self%D)
        real,         intent(inout) :: centered_out(self%D)
        real, allocatable :: resp_one(:), z_one(:,:), zeroavg(:)
        integer :: kbest
        allocate(resp_one(self%K), z_one(self%Q,self%K), zeroavg(self%D), source=0.)
        call self%infer_one(centered_in, resp_one, z_one)
        if( trim(self%recon_mode) .eq. 'hard' )then
            kbest = maxloc(resp_one, dim=1)
            call self%generate_from_component(kbest, z_one(:,kbest), zeroavg, centered_out)
        else
            call self%generate_from_responsibilities(resp_one, z_one, zeroavg, centered_out)
        endif
        deallocate(resp_one, z_one, zeroavg)
    end subroutine reconstruct_external_mppca

    subroutine master_mppca( self, pcavecs, maxpcaits )
        class(mppca),      intent(inout) :: self
        real,              intent(in)    :: pcavecs(self%D,self%N)
        integer, optional, intent(in)    :: maxpcaits
        integer, parameter :: MAX_ITS_DEFAULT = 20
        real(dp), parameter :: TOL_LL = 1.e-5_dp
        real(dp), parameter :: TOL_SIGMA = 1.e-5_dp
        integer :: maxits, it
        integer(int64) :: t0, t1
        real(real64)   :: trate
        real(dp) :: ll_prev, ll_cur, sigma_rel
        if( present(maxpcaits) )then
            maxits = max(1, maxpcaits)
        else
            maxits = MAX_ITS_DEFAULT
        endif
        if( .not. self%existence ) call self%init_mixture()
        call system_clock(t0, trate)
        write(logfhandle,'(A,I8,A,I8,A,I8,A,I8,A,I8)') 'mPPCA start: N=', self%N, ' D=', self%D, &
            ' Q=', self%Q, ' K=', self%K, ' maxits=', maxits
        call flush(logfhandle)
        call self%init_from_ppca(pcavecs, max(5, min(10, maxits)))
        ll_prev = -huge(1._dp)
        do it = 1, maxits
            call self%em_step(pcavecs, ll_cur, sigma_rel)
            if( it == 1 .or. mod(it,5) == 0 .or. it == maxits )then
                write(logfhandle,'(A,I4,A,ES12.4,A,ES10.3,A,2F8.4)') 'mPPCA iter=', it, ' loglik=', ll_cur, &
                    ' dsigma=', sigma_rel, ' pi[min,max]=', minval(self%pi), maxval(self%pi)
                call flush(logfhandle)
            endif
            if( it > 1 )then
                if( abs(ll_cur - ll_prev) <= TOL_LL * max(1._dp, abs(ll_prev)) .and. sigma_rel <= TOL_SIGMA )then
                    write(logfhandle,'(A,I4,A)') 'mPPCA converged at iter=', it, ''
                    call flush(logfhandle)
                    exit
                endif
            endif
            ll_prev = ll_cur
        enddo
        call self%reconstruct()
        call self%calc_eigvals()
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A,2F8.4)') 'mPPCA total: ', real(t1-t0)/real(trate), &
            ' s; pi[min,max]=', minval(self%pi), maxval(self%pi)
        call flush(logfhandle)
    end subroutine master_mppca

    subroutine init_mixture( self )
        class(mppca), intent(inout) :: self
        integer :: i
        if( self%existence )then
            if( allocated(self%pi)       ) deallocate(self%pi)
            if( allocated(self%sigma2)   ) deallocate(self%sigma2)
            if( allocated(self%mu)       ) deallocate(self%mu)
            if( allocated(self%W)        ) deallocate(self%W)
            if( allocated(self%E_zn)     ) deallocate(self%E_zn)
            if( allocated(self%resp)     ) deallocate(self%resp)
            if( allocated(self%data)     ) deallocate(self%data)
            if( allocated(self%feat_mix) ) deallocate(self%feat_mix)
            if( allocated(self%eigvals)  ) deallocate(self%eigvals)
            if( allocated(self%Iq)       ) deallocate(self%Iq)
        endif
        allocate(self%pi(self%K), self%sigma2(self%K), self%mu(self%D,self%K), self%W(self%D,self%Q,self%K), &
                 self%E_zn(self%Q,self%N,self%K), self%resp(self%K,self%N), self%data(self%D,self%N), &
                 self%feat_mix(self%Q,self%N), self%eigvals(self%Q), self%Iq(self%Q,self%Q), source=0.)
        self%pi = 1. / real(self%K)
        self%sigma2 = 1.
        self%Iq = 0.
        do i = 1,self%Q
            self%Iq(i,i) = 1.
        enddo
        self%existence = .true.
    end subroutine init_mixture

    subroutine init_from_ppca( self, pcavecs, init_its )
        use simple_rnd, only: ran3
        class(mppca), intent(inout) :: self
        real,         intent(in)    :: pcavecs(self%D,self%N)
        integer,      intent(in)    :: init_its
        type(ppca) :: base_ppca
        real, allocatable :: feats(:,:), centers(:,:), feat(:)
        real(dp), allocatable :: dmin(:)
        real(dp) :: mean_sq, best_d, cur_d
        integer :: i, j, k, seed_idx, assign(self%N), cnt(self%K), center_idx(self%K)
        call base_ppca%new(self%N, self%D, self%Q)
        call base_ppca%set_verbose(.false.)
        call base_ppca%master(pcavecs, init_its)
        allocate(feats(self%Q,self%N), centers(self%Q,self%K), source=0.)
        allocate(dmin(self%N), source=0._dp)
        do i = 1,self%N
            feat = base_ppca%get_feat(i)
            feats(:,i) = feat
            deallocate(feat)
        enddo
        center_idx(1) = maxloc(sum(real(feats,dp)**2, dim=1), dim=1)
        centers(:,1)  = feats(:,center_idx(1))
        dmin = huge(1._dp)
        do k = 2,self%K
            do i = 1,self%N
                cur_d = sum((real(feats(:,i),dp) - real(centers(:,k-1),dp))**2)
                dmin(i) = min(dmin(i), cur_d)
            enddo
            seed_idx = maxloc(dmin, dim=1)
            center_idx(k) = seed_idx
            centers(:,k) = feats(:,seed_idx)
        enddo
        cnt = 0
        do i = 1,self%N
            best_d = huge(1._dp)
            assign(i) = 1
            do k = 1,self%K
                cur_d = sum((real(feats(:,i),dp) - real(centers(:,k),dp))**2)
                if( cur_d < best_d )then
                    best_d = cur_d
                    assign(i) = k
                endif
            enddo
            cnt(assign(i)) = cnt(assign(i)) + 1
        enddo
        do k = 1,self%K
            if( cnt(k) == 0 )then
                assign(center_idx(k)) = k
                cnt(k) = 1
            endif
        enddo
        self%mu = 0.
        self%resp = 0.
        do i = 1,self%N
            k = assign(i)
            self%mu(:,k) = self%mu(:,k) + pcavecs(:,i)
            self%resp(k,i) = 1.
        enddo
        do k = 1,self%K
            self%mu(:,k) = self%mu(:,k) / real(cnt(k))
            self%pi(k) = real(cnt(k)) / real(self%N)
        enddo
        mean_sq = sum(real(pcavecs,dp)**2) / real(max(1,self%D*self%N), dp)
        do k = 1,self%K
            self%sigma2(k) = real(max(mean_sq * 0.1_dp, real(DTINY,dp)))
            do j = 1,self%Q
                do i = 1,self%D
                    self%W(i,j,k) = 0.01 * ran3()
                enddo
            enddo
        enddo
        call base_ppca%kill()
        deallocate(feats, centers, dmin)
    end subroutine init_from_ppca

    subroutine em_step( self, pcavecs, ll_out, sigma_rel )
        use simple_math, only: matinv
        class(mppca), intent(inout) :: self
        real,         intent(in)    :: pcavecs(self%D,self%N)
        real(dp),     intent(out)   :: ll_out, sigma_rel
        integer, parameter :: BLOCKSZ = 256
        real(dp), parameter :: PI_CONST = 3.14159265358979323846_dp
        real(dp), parameter :: SIGMA_FLOOR = 1.e-8_dp
        real(dp), parameter :: MIX_WEIGHT_SHRINK = 0.10_dp
        real(dp), parameter :: SIGMA_REG_FRAC = 0.10_dp
        integer :: k, i, ib, ie, nb, err
        real, allocatable :: M(:,:,:), Minv(:,:,:), MinvWt(:,:,:), S_xz(:,:,:), S_zz(:,:,:), W_new(:,:,:), mu_new(:,:), WTW(:,:,:)
        real, allocatable :: xc_blk(:,:), ez_blk(:,:), mez_blk(:,:)
        real(dp), allocatable :: Nk(:), sigma_new(:), sigma_old(:), logdetM(:)
        real(dp), allocatable :: logp_loc(:,:), xnorm2(:), ez_mez(:)
        real(dp) :: logsumexp, quad, sigma_num, sigma_trace_term, weighted_loglik, sigma_model_floor
        allocate(M(self%Q,self%Q,self%K), Minv(self%Q,self%Q,self%K), MinvWt(self%Q,self%D,self%K), &
                 S_xz(self%D,self%Q,self%K), S_zz(self%Q,self%Q,self%K), W_new(self%D,self%Q,self%K), &
                 mu_new(self%D,self%K), WTW(self%Q,self%Q,self%K), source=0.)
        allocate(Nk(self%K), sigma_new(self%K), sigma_old(self%K), logdetM(self%K), source=0._dp)
        sigma_old = real(self%sigma2, dp)
        do k = 1,self%K
            M(:,:,k) = matmul(transpose(self%W(:,:,k)), self%W(:,:,k)) + self%sigma2(k) * self%Iq
            call matinv(M(:,:,k), Minv(:,:,k), self%Q, err)
            if( err == -1 )then
                Minv(:,:,k) = 0.
                do i = 1,self%Q
                    Minv(i,i,k) = 1.
                enddo
            endif
            MinvWt(:,:,k) = matmul(Minv(:,:,k), transpose(self%W(:,:,k)))
            do i = 1,self%Q
                logdetM(k) = logdetM(k) + log(max(real(M(i,i,k),dp), SIGMA_FLOOR))
            enddo
        enddo
        ll_out = 0._dp
        do ib = 1, self%N, BLOCKSZ
            ie = min(self%N, ib + BLOCKSZ - 1)
            nb = ie - ib + 1
            allocate(logp_loc(self%K,nb), xnorm2(nb), ez_mez(nb), source=0._dp)
            do k = 1,self%K
                allocate(xc_blk(self%D,nb), ez_blk(self%Q,nb), mez_blk(self%Q,nb), source=0.)
                xc_blk = pcavecs(:,ib:ie) - spread(self%mu(:,k), dim=2, ncopies=nb)
                ez_blk = matmul(MinvWt(:,:,k), xc_blk)
                self%E_zn(:,ib:ie,k) = ez_blk
                mez_blk = matmul(M(:,:,k), ez_blk)
                xnorm2 = sum(real(xc_blk,dp)**2, dim=1)
                ez_mez = sum(real(ez_blk,dp) * real(mez_blk,dp), dim=1)
                do i = 1, nb
                    quad = (xnorm2(i) - ez_mez(i)) / max(real(self%sigma2(k),dp), SIGMA_FLOOR)
                    logp_loc(k,i) = log(max(real(self%pi(k),dp), SIGMA_FLOOR)) - 0.5_dp * &
                        ( real(self%D-self%Q,dp) * log(max(real(self%sigma2(k),dp), SIGMA_FLOOR)) + &
                          logdetM(k) + quad + real(self%D,dp)*log(2._dp*PI_CONST) )
                enddo
                deallocate(xc_blk, ez_blk, mez_blk)
            enddo
            do i = 1, nb
                logsumexp = maxval(logp_loc(:,i))
                weighted_loglik = 0._dp
                do k = 1,self%K
                    logp_loc(k,i) = exp(logp_loc(k,i) - logsumexp)
                    weighted_loglik = weighted_loglik + logp_loc(k,i)
                enddo
                logsumexp = logsumexp + log(max(weighted_loglik, SIGMA_FLOOR))
                ll_out = ll_out + logsumexp
                do k = 1,self%K
                    self%resp(k,ib+i-1) = real(logp_loc(k,i) / max(weighted_loglik, SIGMA_FLOOR))
                enddo
            enddo
            deallocate(logp_loc, xnorm2, ez_mez)
        enddo
        Nk = sum(real(self%resp,dp), dim=2)
        do k = 1,self%K
            Nk(k) = max(Nk(k), 1.e-6_dp)
            self%pi(k) = real((1._dp - MIX_WEIGHT_SHRINK) * (Nk(k) / real(self%N,dp)) + MIX_WEIGHT_SHRINK / real(self%K,dp))
        enddo
        mu_new = 0.
        do k = 1,self%K
            !$omp parallel default(shared) private(ib,ie,nb) num_threads(self%nthr)
            block
                real, allocatable :: xc_loc(:,:), mu_loc(:), wresp_loc(:)
                allocate(xc_loc(self%D,BLOCKSZ), mu_loc(self%D), wresp_loc(BLOCKSZ), source=0.)
                !$omp do schedule(dynamic)
                do ib = 1, self%N, BLOCKSZ
                    ie = min(self%N, ib + BLOCKSZ - 1)
                    nb = ie - ib + 1
                    xc_loc(:,1:nb) = pcavecs(:,ib:ie) - matmul(self%W(:,:,k), self%E_zn(:,ib:ie,k))
                    wresp_loc(1:nb) = self%resp(k,ib:ie)
                    mu_loc = mu_loc + matmul(xc_loc(:,1:nb), wresp_loc(1:nb))
                enddo
                !$omp end do
                !$omp critical(mppca_mu_reduce)
                mu_new(:,k) = mu_new(:,k) + mu_loc
                !$omp end critical(mppca_mu_reduce)
                deallocate(xc_loc, mu_loc, wresp_loc)
            end block
            !$omp end parallel
            mu_new(:,k) = mu_new(:,k) / real(Nk(k))
        enddo
        S_xz = 0.
        S_zz = 0.
        sigma_new = 0._dp
        sigma_model_floor = max(SIGMA_FLOOR, SIGMA_REG_FRAC * sum(sigma_old) / real(self%K,dp))
        do k = 1,self%K
            !$omp parallel default(shared) private(ib,ie,nb) num_threads(self%nthr)
            block
                real, allocatable :: xc_loc(:,:), ez_loc(:,:), weighted_ez_loc(:,:), ez_s_loc(:,:)
                real, allocatable :: sxz_loc(:,:), szz_loc(:,:)
                real(dp), allocatable :: rw_loc(:)
                allocate(xc_loc(self%D,BLOCKSZ), ez_loc(self%Q,BLOCKSZ), weighted_ez_loc(self%Q,BLOCKSZ), &
                         ez_s_loc(self%Q,BLOCKSZ), sxz_loc(self%D,self%Q), szz_loc(self%Q,self%Q), source=0.)
                allocate(rw_loc(BLOCKSZ), source=0._dp)
                !$omp do schedule(dynamic)
                do ib = 1, self%N, BLOCKSZ
                    ie = min(self%N, ib + BLOCKSZ - 1)
                    nb = ie - ib + 1
                    xc_loc(:,1:nb) = pcavecs(:,ib:ie) - spread(mu_new(:,k), dim=2, ncopies=nb)
                    ez_loc(:,1:nb) = self%E_zn(:,ib:ie,k)
                    rw_loc(1:nb) = real(self%resp(k,ib:ie), dp)
                    weighted_ez_loc(:,1:nb) = ez_loc(:,1:nb) * spread(real(rw_loc(1:nb)), dim=1, ncopies=self%Q)
                    ez_s_loc(:,1:nb) = ez_loc(:,1:nb) * spread(real(sqrt(max(rw_loc(1:nb),0._dp))), dim=1, ncopies=self%Q)
                    sxz_loc = sxz_loc + matmul(xc_loc(:,1:nb), transpose(weighted_ez_loc(:,1:nb)))
                    szz_loc = szz_loc + matmul(ez_s_loc(:,1:nb), transpose(ez_s_loc(:,1:nb)))
                enddo
                !$omp end do
                !$omp critical(mppca_s_reduce)
                S_xz(:,:,k) = S_xz(:,:,k) + sxz_loc
                S_zz(:,:,k) = S_zz(:,:,k) + szz_loc
                !$omp end critical(mppca_s_reduce)
                deallocate(xc_loc, ez_loc, weighted_ez_loc, ez_s_loc, sxz_loc, szz_loc, rw_loc)
            end block
            !$omp end parallel
            S_zz(:,:,k) = S_zz(:,:,k) + real(self%sigma2(k) * Nk(k)) * Minv(:,:,k)
            call matinv(S_zz(:,:,k), M(:,:,k), self%Q, err)
            if( err == -1 )then
                M(:,:,k) = 0.
                do i = 1,self%Q
                    M(i,i,k) = 1.
                enddo
            endif
            W_new(:,:,k) = matmul(S_xz(:,:,k), M(:,:,k))
            WTW(:,:,k) = matmul(transpose(W_new(:,:,k)), W_new(:,:,k))
            sigma_trace_term = self%sigma2(k) * sum(real(WTW(:,:,k),dp) * transpose(real(Minv(:,:,k),dp)))
            sigma_num = 0._dp
            !$omp parallel default(shared) private(ib,ie,nb) reduction(+:sigma_num) num_threads(self%nthr)
            block
                real, allocatable :: xc_loc(:,:), ez_loc(:,:), wt_xc_loc(:,:)
                real(dp), allocatable :: rw_loc(:), term1_loc(:), term2_loc(:), term3_loc(:)
                allocate(xc_loc(self%D,BLOCKSZ), ez_loc(self%Q,BLOCKSZ), wt_xc_loc(self%Q,BLOCKSZ), source=0.)
                allocate(rw_loc(BLOCKSZ), term1_loc(BLOCKSZ), term2_loc(BLOCKSZ), term3_loc(BLOCKSZ), source=0._dp)
                !$omp do schedule(dynamic)
                do ib = 1, self%N, BLOCKSZ
                    ie = min(self%N, ib + BLOCKSZ - 1)
                    nb = ie - ib + 1
                    xc_loc(:,1:nb) = pcavecs(:,ib:ie) - spread(mu_new(:,k), dim=2, ncopies=nb)
                    ez_loc(:,1:nb) = self%E_zn(:,ib:ie,k)
                    wt_xc_loc(:,1:nb) = matmul(transpose(W_new(:,:,k)), xc_loc(:,1:nb))
                    rw_loc(1:nb) = real(self%resp(k,ib:ie), dp)
                    term1_loc(1:nb) = sum(real(xc_loc(:,1:nb),dp)**2, dim=1)
                    term2_loc(1:nb) = 2._dp * sum(real(ez_loc(:,1:nb),dp) * real(wt_xc_loc(:,1:nb),dp), dim=1)
                    term3_loc(1:nb) = sigma_trace_term + sum(real(ez_loc(:,1:nb),dp) * real(matmul(WTW(:,:,k), ez_loc(:,1:nb)),dp), dim=1)
                    sigma_num = sigma_num + sum(rw_loc(1:nb) * (term1_loc(1:nb) - term2_loc(1:nb) + term3_loc(1:nb)))
                enddo
                !$omp end do
                deallocate(xc_loc, ez_loc, wt_xc_loc, rw_loc, term1_loc, term2_loc, term3_loc)
            end block
            !$omp end parallel
            sigma_new(k) = max(sigma_model_floor, sigma_num / (real(self%D,dp) * Nk(k)))
        enddo
        self%mu = mu_new
        self%W = W_new
        self%sigma2 = real(sigma_new)
        sigma_rel = maxval(abs(sigma_new - sigma_old) / max(SIGMA_FLOOR, sigma_old))
        deallocate(M, Minv, MinvWt, S_xz, S_zz, W_new, mu_new, WTW, Nk, sigma_new, sigma_old, logdetM)
    end subroutine em_step

    subroutine reconstruct( self )
        class(mppca), intent(inout) :: self
        integer :: i, k, kbest
        self%data = 0.
        self%feat_mix = 0.
        if( trim(self%recon_mode) .eq. 'hard' )then
            do i = 1,self%N
                kbest = maxloc(self%resp(:,i), dim=1)
                self%data(:,i) = self%mu(:,kbest) + matmul(self%W(:,:,kbest), self%E_zn(:,i,kbest))
                self%feat_mix(:,i) = self%E_zn(:,i,kbest)
            enddo
        else
            do i = 1,self%N
                do k = 1,self%K
                    self%data(:,i) = self%data(:,i) + self%resp(k,i) * (self%mu(:,k) + matmul(self%W(:,:,k), self%E_zn(:,i,k)))
                    self%feat_mix(:,i) = self%feat_mix(:,i) + self%resp(k,i) * self%E_zn(:,i,k)
                enddo
            enddo
        endif
    end subroutine reconstruct

    subroutine infer_one( self, centered_in, resp_one, z_one )
        use simple_math, only: matinv
        class(mppca), intent(inout) :: self
        real,         intent(in)    :: centered_in(self%D)
        real,         intent(out)   :: resp_one(self%K), z_one(self%Q,self%K)
        real(dp), parameter :: PI_CONST = 3.14159265358979323846_dp
        real(dp), parameter :: SIGMA_FLOOR = 1.e-8_dp
        real, allocatable :: M(:,:), Minv(:,:), MinvWt(:,:), xc(:)
        real(dp), allocatable :: logp(:)
        real(dp) :: logdetM, quad, xnorm2, ez_mez, logsumexp, wsum
        integer :: k, i, err
        allocate(M(self%Q,self%Q), Minv(self%Q,self%Q), MinvWt(self%Q,self%D), xc(self%D), source=0.)
        allocate(logp(self%K), source=0._dp)
        do k = 1,self%K
            M = matmul(transpose(self%W(:,:,k)), self%W(:,:,k)) + self%sigma2(k) * self%Iq
            call matinv(M, Minv, self%Q, err)
            if( err == -1 )then
                Minv = 0.
                do i = 1,self%Q
                    Minv(i,i) = 1.
                enddo
            endif
            MinvWt = matmul(Minv, transpose(self%W(:,:,k)))
            xc = centered_in - self%mu(:,k)
            z_one(:,k) = matmul(MinvWt, xc)
            logdetM = 0._dp
            do i = 1,self%Q
                logdetM = logdetM + log(max(real(M(i,i),dp), SIGMA_FLOOR))
            enddo
            xnorm2 = sum(real(xc,dp)**2)
            ez_mez = sum(real(z_one(:,k),dp) * real(matmul(M, z_one(:,k)),dp))
            quad = (xnorm2 - ez_mez) / max(real(self%sigma2(k),dp), SIGMA_FLOOR)
            logp(k) = log(max(real(self%pi(k),dp), SIGMA_FLOOR)) - 0.5_dp * &
                ( real(self%D-self%Q,dp) * log(max(real(self%sigma2(k),dp), SIGMA_FLOOR)) + &
                  logdetM + quad + real(self%D,dp)*log(2._dp*PI_CONST) )
        enddo
        logsumexp = maxval(logp)
        wsum = 0._dp
        do k = 1,self%K
            logp(k) = exp(logp(k) - logsumexp)
            wsum = wsum + logp(k)
        enddo
        do k = 1,self%K
            resp_one(k) = real(logp(k) / max(wsum, SIGMA_FLOOR))
        enddo
        deallocate(M, Minv, MinvWt, xc, logp)
    end subroutine infer_one

    subroutine calc_eigvals( self )
        class(mppca), intent(inout) :: self
        real(dp), allocatable :: gram(:,:), rot(:,:), vals(:)
        integer :: k, nrot
        self%eigvals = 0.
        if( self%Q <= 0 ) return
        allocate(gram(self%Q,self%Q), rot(self%Q,self%Q), vals(self%Q), source=0._dp)
        do k = 1,self%K
            gram = matmul(transpose(real(self%W(:,:,k),dp)), real(self%W(:,:,k),dp))
            gram = 0.5_dp * (gram + transpose(gram))
            nrot = 0
            call jacobi(gram, self%Q, self%Q, vals, rot, nrot)
            call eigsrt(vals, rot, self%Q, self%Q)
            vals = max(vals, 0._dp)
            self%eigvals = self%eigvals + self%pi(k) * real(vals)
        enddo
        deallocate(gram, rot, vals)
    end subroutine calc_eigvals

    subroutine kill_mppca( self )
        class(mppca), intent(inout) :: self
        if( self%existence )then
            if( allocated(self%pi)       ) deallocate(self%pi)
            if( allocated(self%sigma2)   ) deallocate(self%sigma2)
            if( allocated(self%mu)       ) deallocate(self%mu)
            if( allocated(self%W)        ) deallocate(self%W)
            if( allocated(self%E_zn)     ) deallocate(self%E_zn)
            if( allocated(self%resp)     ) deallocate(self%resp)
            if( allocated(self%data)     ) deallocate(self%data)
            if( allocated(self%feat_mix) ) deallocate(self%feat_mix)
            if( allocated(self%eigvals)  ) deallocate(self%eigvals)
            if( allocated(self%Iq)       ) deallocate(self%Iq)
            self%existence = .false.
        endif
    end subroutine kill_mppca

end module simple_mppca
