!@descr: bounded BLAS-backed Pearson correlation for reference picking
module simple_pickref_corr_batch
use simple_defs,  only: dp, TINY, logfhandle
use simple_image, only: image
use simple_linalg, only: gemm_tn
use simple_timer_omp, only: tic_omp, toc_omp
implicit none
private

public :: pickref_corr_batch

integer, parameter :: DEFAULT_BATCH_SIZE = 256

type :: pickref_corr_batch
    private
    integer :: nx = 0, ny = 0, box = 0, npix = 0, nrefs = 0, batch_size = 0
    real, allocatable :: mic(:,:), refs(:,:), ref_sums(:)
    real, allocatable :: windows(:,:), dots(:,:)
    real(dp), allocatable :: sat(:,:), sat2(:,:)
    real(dp) :: mic_scale = 1.0_dp
    real(dp) :: prepare_seconds = 0.0_dp, gather_seconds = 0.0_dp
    real(dp) :: gemm_seconds = 0.0_dp, reduce_seconds = 0.0_dp
    integer(kind=8) :: candidates_scored = 0_8
    logical, allocatable :: ref_valid(:)
    logical :: exists = .false.
contains
    procedure :: new
    procedure :: score_positions
    procedure :: report_stats
    procedure :: kill
end type pickref_corr_batch

contains

    subroutine new( self, mic, refs, ref_valid, batch_size )
        class(pickref_corr_batch), intent(inout) :: self
        class(image),              intent(in)    :: mic
        class(image),              intent(in)    :: refs(:)
        logical,                   intent(in)    :: ref_valid(:)
        integer, optional,         intent(in)    :: batch_size
        real, allocatable :: rmat(:,:,:)
        real(dp) :: mean_mic, amp, val, t0
        integer :: ldim(3), i, j, iref
        if( self%exists ) call self%kill
        t0 = tic_omp()
        ldim = mic%get_ldim()
        if( ldim(3) /= 1 ) error stop 'pickref_corr_batch: micrograph must be 2D'
        if( size(refs) /= size(ref_valid) ) error stop 'pickref_corr_batch: reference validity size mismatch'
        if( size(refs) == 0 ) error stop 'pickref_corr_batch: empty reference bank'
        self%nx = ldim(1)
        self%ny = ldim(2)
        ldim = refs(1)%get_ldim()
        if( ldim(3) /= 1 .or. ldim(1) /= ldim(2) ) error stop 'pickref_corr_batch: references must be square and 2D'
        self%box = ldim(1)
        self%npix = self%box * self%box
        self%nrefs = size(refs)
        self%batch_size = DEFAULT_BATCH_SIZE
        if( present(batch_size) ) self%batch_size = max(1, batch_size)
        allocate(self%mic(self%nx,self%ny))
        rmat = mic%get_rmat()
        self%mic = rmat(:,:,1)
        deallocate(rmat)
        mean_mic = sum(real(self%mic,dp)) / real(self%nx*self%ny,dp)
        self%mic = self%mic - real(mean_mic)
        amp = maxval(abs(real(self%mic,dp)))
        self%mic_scale = 1.0_dp
        if( amp > real(TINY,dp) )then
            self%mic_scale = amp
            self%mic = self%mic / real(amp)
        endif
        allocate(self%sat(0:self%nx,0:self%ny), source=0.0_dp)
        allocate(self%sat2(0:self%nx,0:self%ny), source=0.0_dp)
        do j = 1,self%ny
            do i = 1,self%nx
                val = real(self%mic(i,j),dp)
                self%sat(i,j)  = val     + self%sat(i-1,j)  + self%sat(i,j-1)  - self%sat(i-1,j-1)
                self%sat2(i,j) = val*val + self%sat2(i-1,j) + self%sat2(i,j-1) - self%sat2(i-1,j-1)
            end do
        end do
        allocate(self%refs(self%npix,self%nrefs), self%ref_sums(self%nrefs))
        allocate(self%ref_valid(self%nrefs), source=ref_valid)
        do iref = 1,self%nrefs
            ldim = refs(iref)%get_ldim()
            if( any(ldim /= [self%box,self%box,1]) ) error stop 'pickref_corr_batch: reference dimensions differ'
            rmat = refs(iref)%get_rmat()
            self%refs(:,iref) = reshape(rmat(:,:,1), [self%npix])
            self%ref_sums(iref) = sum(self%refs(:,iref))
            deallocate(rmat)
        end do
        allocate(self%windows(self%npix,self%batch_size), source=0.0)
        allocate(self%dots(self%nrefs,self%batch_size), source=0.0)
        self%prepare_seconds = toc_omp(t0)
        self%exists = .true.
    end subroutine new

    subroutine score_positions( self, positions, scores )
        class(pickref_corr_batch), intent(inout) :: self
        integer,                   intent(in)    :: positions(:,:)
        real,                      intent(out)   :: scores(size(positions,1))
        real(dp), allocatable :: means(:), sigmas(:)
        logical, allocatable :: valid(:)
        logical :: have_score
        real(dp) :: sum_win, sumsq_win, variance, rnpix, t0
        real :: score, best_score
        integer :: ibeg, iend, nbatch, icol, ipos, iref, x, y
        integer :: i0, i1, j0, j1, jwin, ivec
        if( .not. self%exists ) error stop 'pickref_corr_batch: score_positions before new'
        if( size(positions,2) /= 2 ) error stop 'pickref_corr_batch: positions must have two columns'
        scores = -1.0
        self%candidates_scored = self%candidates_scored + int(size(positions,1),kind=8)
        allocate(means(self%batch_size), sigmas(self%batch_size), valid(self%batch_size))
        rnpix = real(self%npix,dp)
        do ibeg = 1,size(positions,1),self%batch_size
            iend = min(size(positions,1), ibeg+self%batch_size-1)
            nbatch = iend-ibeg+1
            self%windows(:,:nbatch) = 0.0
            means(:nbatch) = 0.0_dp
            sigmas(:nbatch) = 0.0_dp
            valid(:nbatch) = .false.
            t0 = tic_omp()
            do icol = 1,nbatch
                ipos = ibeg+icol-1
                x = positions(ipos,1)
                y = positions(ipos,2)
                if( x < 0 .or. y < 0 ) cycle
                if( x+self%box > self%nx .or. y+self%box > self%ny ) cycle
                i0 = x+1
                i1 = x+self%box
                j0 = y+1
                j1 = y+self%box
                sum_win = self%sat(i1,j1) - self%sat(i0-1,j1) - self%sat(i1,j0-1) + self%sat(i0-1,j0-1)
                sumsq_win = self%sat2(i1,j1) - self%sat2(i0-1,j1) - self%sat2(i1,j0-1) + self%sat2(i0-1,j0-1)
                means(icol) = sum_win / rnpix
                variance = sumsq_win / rnpix - means(icol)*means(icol)
                if( .not.(variance*self%mic_scale*self%mic_scale > real(TINY,dp)) ) cycle
                sigmas(icol) = sqrt(max(0.0_dp,variance))
                ivec = 0
                do jwin = j0,j1
                    self%windows(ivec+1:ivec+self%box,icol) = self%mic(i0:i1,jwin)
                    ivec = ivec+self%box
                end do
                valid(icol) = .true.
            end do
            self%gather_seconds = self%gather_seconds + toc_omp(t0)
            t0 = tic_omp()
            call gemm_tn(self%refs, self%windows(:,:nbatch), self%dots(:,:nbatch))
            self%gemm_seconds = self%gemm_seconds + toc_omp(t0)
            t0 = tic_omp()
            do icol = 1,nbatch
                if( .not. valid(icol) ) cycle
                best_score = -huge(best_score)
                have_score = .false.
                do iref = 1,self%nrefs
                    if( .not. self%ref_valid(iref) ) cycle
                    score = (self%dots(iref,icol) - real(means(icol))*self%ref_sums(iref)) &
                        / (real(self%npix)*real(sigmas(icol)))
                    if( .not. have_score .or. score > best_score )then
                        best_score = score
                        have_score = .true.
                    endif
                end do
                if( have_score ) scores(ibeg+icol-1) = best_score
            end do
            self%reduce_seconds = self%reduce_seconds + toc_omp(t0)
        end do
        deallocate(means, sigmas, valid)
    end subroutine score_positions

    subroutine report_stats( self, label )
        class(pickref_corr_batch), intent(in) :: self
        character(len=*),          intent(in) :: label
        real(dp) :: memory_bytes
        if( .not. self%exists ) return
        memory_bytes = real(size(self%mic)+size(self%refs)+size(self%ref_sums)+size(self%windows)+size(self%dots),dp) * 4.0_dp &
            + real(size(self%sat)+size(self%sat2),dp) * 8.0_dp + real(size(self%ref_valid),dp)
        write(logfhandle,'(a,1x,a,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,f8.2)') 'refpick correlation', trim(label), &
            'candidates=', self%candidates_scored, 'refs=', self%nrefs, 'batch=', self%batch_size, &
            'workspace_mb=', memory_bytes/(1024.0_dp*1024.0_dp)
        write(logfhandle,'(a,1x,a,4(1x,a,es10.3))') 'refpick correlation', trim(label), &
            'prepare=', self%prepare_seconds, 'gather=', self%gather_seconds, &
            'gemm=', self%gemm_seconds, 'reduce=', self%reduce_seconds
    end subroutine report_stats

    subroutine kill( self )
        class(pickref_corr_batch), intent(inout) :: self
        if( allocated(self%mic) ) deallocate(self%mic)
        if( allocated(self%refs) ) deallocate(self%refs)
        if( allocated(self%ref_sums) ) deallocate(self%ref_sums)
        if( allocated(self%windows) ) deallocate(self%windows)
        if( allocated(self%dots) ) deallocate(self%dots)
        if( allocated(self%sat) ) deallocate(self%sat)
        if( allocated(self%sat2) ) deallocate(self%sat2)
        if( allocated(self%ref_valid) ) deallocate(self%ref_valid)
        self%nx = 0
        self%ny = 0
        self%box = 0
        self%npix = 0
        self%nrefs = 0
        self%batch_size = 0
        self%mic_scale = 1.0_dp
        self%prepare_seconds = 0.0_dp
        self%gather_seconds = 0.0_dp
        self%gemm_seconds = 0.0_dp
        self%reduce_seconds = 0.0_dp
        self%candidates_scored = 0_8
        self%exists = .false.
    end subroutine kill

end module simple_pickref_corr_batch
