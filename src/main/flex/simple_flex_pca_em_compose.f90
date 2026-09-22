!@descr: flex_pca EM: multi-band basis composition from finished runs (SIMPLE_COV_COMPOSE)
!!
!! PCA at one band returns the leading eigenvectors of the covariance PROJECTED onto that band, and
!! the eigenvectors of a projection are not the projection of the eigenvectors: a mode that ranks
!! third at 30 A can fall below the rank cut at 16 A, and a mode that needs 20 A detail has no
!! power at 30 A. Marching one basis through the bands (refit or extension) re-orders it at every
!! stage and was measured to erode the coarse-band structure. Composition keeps every band's basis
!! as delivered: the columns of each finished run are loaded on their own grid, Fourier-padded to
!! the composing box (a coarse column stays exactly zero beyond
!! its own band -- it is never noise-fitted in shells it never saw), orthonormalised coarse-first
!! (a fine column keeps only what is new relative to the coarse subspace), and embedded ONCE with
!! the union basis. The latent prior variance of each column follows its rescaling, so the
!! embedding's MAP shrinkage is unchanged in the coarse directions.
!!
!! SIMPLE_COV_COMPOSE=<dir>[,<dir>...]: finished flex_pca run directories; any order (sorted by
!! their crop box, coarse first). Per run the polished namespace is preferred, then merged, then
!! the plain flex_pca_pc*.mrc, each with its own probe meta (rank, sig2, prior variances).
!! Measured and rejected (2026-09-10): per-shell cross-half weighting of each column gutted the
!! leading axes (two 50k-particle half fits disagree beyond 80 A); fine-first ordering was neutral;
!! a per-component half-fit FSC >= 0.5 gate dropped union structure without buying any.
submodule (simple_flex_pca_em) simple_flex_pca_em_compose
use simple_imghead,         only: find_ldim_nptcls
use simple_flex_pca_util,   only: cov_env_flag_on, cov_env_flag_off, cov_env_dp
use simple_flex_pca_deconv, only: calibrate_noise_scale, noise_and_projection, signal_subspace
implicit none
#include "simple_local_flags.inc"

!> a finished run's delivered basis: namespace, rank, grid, prior variances
type compose_src_t
    character(len=LONGSTRLEN) :: dir  = ''
    character(len=32)         :: pfx  = ''
    character(len=64)         :: meta = ''
    integer  :: box   = 0
    integer  :: ncomp = 0
    real(dp) :: sig2  = 0.d0
    real(dp), allocatable :: eigvals(:)
end type compose_src_t

!> a composed column: a fine column whose residual after the coarse projection is below this
!! fraction of its norm is already spanned and is dropped
real(dp), parameter :: COMPOSE_R2_FLOOR = 0.05d0

contains

    module subroutine compose_basis_from_runs( params, build, basis_recs, eigvals, ncomp, sig2_eff )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        real(dp),            allocatable, intent(out) :: eigvals(:)
        integer,             intent(out)   :: ncomp
        real(dp),            intent(out)   :: sig2_eff
        type(compose_src_t), allocatable :: srcs(:)
        type(image),         allocatable :: imgs(:)
        type(image)  :: raw, fine
        type(string) :: fname
        real(dp), allocatable :: VV(:,:), qcol(:), ev(:), s2(:), r2(:)
        integer,  allocatable :: c_src(:), c_comp(:), c_box(:)
        real,     pointer     :: rmat(:,:,:) => null()
        character(len=XLONGSTRLEN) :: env
        real(dp) :: pj
        real     :: smpd_s
        integer  :: envlen, envstat, nsrc, isrc, ntot, bc, bs, nvox, j, k, q, kept
        integer  :: u
        bc   = params%box_crop
        nvox = bc*bc*bc
        env  = ''
        call get_environment_variable('SIMPLE_COV_COMPOSE', env, envlen, envstat)
        if( envstat /= 0 .or. envlen < 1 ) THROW_HARD('compose_basis_from_runs: SIMPLE_COV_COMPOSE is not set')
        call parse_sources(env(1:envlen), srcs, nsrc)
        do isrc = 1, nsrc
            call resolve_source(params, srcs(isrc))
        end do
        call sort_sources_by_box(srcs, nsrc)
        ntot = sum(srcs(1:nsrc)%ncomp)
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA COMPOSE: ', nsrc, ' finished run(s), ', ntot, &
            &' delivered components, composing box=', bc
        do isrc = 1, nsrc
            write(logfhandle,'(A,I0,A,I0,A,I0,A,ES10.3,A,A,A,A)') '>>>   run ', isrc, ': box ', &
                &srcs(isrc)%box, '  ncomp ', srcs(isrc)%ncomp, '  sig2 ', srcs(isrc)%sig2, '  namespace ', &
                &trim(srcs(isrc)%pfx), '  ', trim(srcs(isrc)%dir)
        end do
        call flush(logfhandle)
        allocate(VV(nvox,ntot), ev(ntot), s2(ntot), r2(ntot), c_src(ntot), c_comp(ntot), c_box(ntot))
        ev = 0.d0; s2 = 0.d0; r2 = 0.d0; c_src = 0; c_comp = 0; c_box = 0
        ! ---- load every delivered column on its own grid, weight, pad to the composing grid ----
        j = 0
        do isrc = 1, nsrc
            bs       = srcs(isrc)%box
            smpd_s   = params%smpd_crop*real(bc)/real(bs)
            call raw%new([bs,bs,bs], smpd_s)
            do q = 1, srcs(isrc)%ncomp
                j = j + 1
                c_src(j) = isrc; c_comp(j) = q; c_box(j) = bs
                ev(j)    = srcs(isrc)%eigvals(q)
                fname = trim(srcs(isrc)%dir)//trim(srcs(isrc)%pfx)//int2str_pad(q,3)//MRC_EXT
                call raw%read(fname)
                call fname%kill
                ! Fourier-pad to the composing grid: coefficients copied verbatim, zero beyond the
                ! source band (no antialiasing window -- the band edge is the source's own low-pass)
                if( bs == bc )then
                    call fine%copy(raw)
                else
                    call raw%fft
                    call fine%new([bc,bc,bc], params%smpd_crop)
                    call raw%pad(fine, antialiasing=.false.)
                    call fine%ifft
                    call raw%ifft
                endif
                call fine%set_smpd(params%smpd_crop)
                call fine%get_rmat_ptr(rmat)
                VV(:,j) = reshape(real(rmat(1:bc,1:bc,1:bc),dp), [nvox])
                call fine%kill
            end do
            call raw%kill
        end do
        call flush(logfhandle)
        ! ---- coarse-first Gram-Schmidt on the composing grid; the prior variance follows the
        ! column's rescaling (z_j u_j = (z_j s_j) u_j/s_j) and its retained fraction ----
        allocate(qcol(nvox))
        kept = 0
        do j = 1, ntot
            s2(j) = sum(VV(:,j)**2)
            if( s2(j) <= 1.d-30 )then
                write(logfhandle,'(A,I0,A,I0,A)') '>>>   run ', c_src(j), ' comp ', c_comp(j), ': empty column; dropped'
                r2(j) = 0.d0
                cycle
            endif
            qcol = VV(:,j)/sqrt(s2(j))
            do k = 1, kept
                pj   = dot_product(qcol, VV(:,k))
                qcol = qcol - pj*VV(:,k)
            end do
            r2(j) = sum(qcol**2)
            if( r2(j) < COMPOSE_R2_FLOOR )then
                write(logfhandle,'(A,I0,A,I0,A,F6.3,A)') '>>>   run ', c_src(j), ' comp ', c_comp(j), &
                    &': residual fraction ', real(r2(j)), ' already spanned by the coarser columns; dropped'
                cycle
            endif
            kept = kept + 1
            VV(:,kept)   = qcol/sqrt(r2(j))
            ev(kept)     = ev(j)*s2(j)*r2(j)
            c_src(kept)  = c_src(j); c_comp(kept) = c_comp(j); c_box(kept) = c_box(j)
            s2(kept) = s2(j); r2(kept) = r2(j)
        end do
        deallocate(qcol)
        if( kept < 2 ) THROW_HARD('compose_basis_from_runs: fewer than 2 independent columns survived')
        ncomp    = kept
        sig2_eff = srcs(maxloc(srcs(1:nsrc)%box, dim=1))%sig2   ! the finest run's noise level: its band is the composing band
        allocate(eigvals(ncomp))
        eigvals = ev(1:ncomp)
        write(logfhandle,'(A,I0,A,I0,A,ES10.3)') '>>> FLEX_PCA COMPOSE: ', ncomp, ' of ', ntot, &
            &' columns kept (orthonormal, coarse first); sig2_eff=', sig2_eff
        write(logfhandle,'(A)') '>>>   column  run  comp  box   s2 (fine-grid norm2)   residual   prior var'
        do j = 1, ncomp
            write(logfhandle,'(A,I5,1X,I4,1X,I5,1X,I5,3X,ES12.4,7X,F8.4,3X,ES12.4)') '>>>   ', j, &
                &c_src(j), c_comp(j), c_box(j), s2(j), r2(j), eigvals(j)
        end do
        call flush(logfhandle)
        ! ---- realise: real-space images -> column reconstructors; the composed basis goes to disk
        ! under the plain namespace (the distributed embed workers load flex_pca_pc*.mrc) ----
        allocate(imgs(ncomp))
        do j = 1, ncomp
            call imgs(j)%new([bc,bc,bc], params%smpd_crop)
            call imgs(j)%get_rmat_ptr(rmat)
            rmat(1:bc,1:bc,1:bc) = real(reshape(VV(:,j), [bc,bc,bc]))
            fname = 'flex_pca_pc'//int2str_pad(j,3)//MRC_EXT
            call imgs(j)%write(fname, del_if_exists=.true.)
            call fname%kill
        end do
        deallocate(VV)
        call basis_recs_from_images(params, build, imgs, ncomp, basis_recs)
        do j = 1, ncomp
            call imgs(j)%kill
        end do
        deallocate(imgs)
        call save_probe_state(ncomp, eigvals, sig2_eff)
        ! provenance table
        open(newunit=u, file='flex_pca_compose.txt', status='replace', action='write')
        write(u,'(A)') '# composed basis: column  run  source_comp  source_box  norm2_fine  residual_frac  prior_var'
        do j = 1, ncomp
            write(u,'(I5,1X,I4,1X,I5,1X,I5,1X,ES14.6,1X,F9.5,1X,ES14.6)') j, c_src(j), c_comp(j), &
                &c_box(j), s2(j), r2(j), eigvals(j)
        end do
        do isrc = 1, nsrc
            write(u,'(A,I0,A,I0,A,A,A,A)') '# run ', isrc, ' box ', srcs(isrc)%box, ' namespace ', &
                &trim(srcs(isrc)%pfx), ' dir ', trim(srcs(isrc)%dir)
        end do
        close(u)
        do isrc = 1, nsrc
            if( allocated(srcs(isrc)%eigvals) ) deallocate(srcs(isrc)%eigvals)
        end do
        deallocate(srcs, ev, s2, r2, c_src, c_comp, c_box)
    end subroutine compose_basis_from_runs

    !> comma-separated run directories -> sources (trailing slash normalised)
    subroutine parse_sources( str, srcs, nsrc )
        character(len=*), intent(in) :: str
        type(compose_src_t), allocatable, intent(out) :: srcs(:)
        integer, intent(out) :: nsrc
        character(len=LONGSTRLEN) :: tok
        integer :: i, i0, n, l
        n = 1
        do i = 1, len(str)
            if( str(i:i) == ',' ) n = n + 1
        end do
        allocate(srcs(n))
        nsrc = 0
        i0   = 1
        do i = 1, len(str) + 1
            if( i > len(str) )then
                tok = adjustl(str(i0:len(str)))
            else if( str(i:i) == ',' )then
                tok = adjustl(str(i0:i-1))
            else
                cycle
            endif
            i0 = i + 1
            l  = len_trim(tok)
            if( l < 1 ) cycle
            if( tok(l:l) /= '/' ) tok = tok(1:l)//'/'
            nsrc = nsrc + 1
            srcs(nsrc)%dir = tok
        end do
        if( nsrc < 1 ) THROW_HARD('SIMPLE_COV_COMPOSE names no run directory')
    end subroutine parse_sources

    !> pick the namespace the run embedded with (polished > merged > plain), read its meta,
    !! count its eigenvolumes, read the grid
    subroutine resolve_source( params, src )
        class(parameters),   intent(in)    :: params
        type(compose_src_t), intent(inout) :: src
        type(string) :: fn
        integer :: ldim(3), nvol, q, nfound
        if( file_exists(trim(src%dir)//'flex_pca_polished_pc001.mrc') .and. &
           &file_exists(trim(src%dir)//'flex_pca_probe_polished.txt') )then
            src%pfx = 'flex_pca_polished_pc'; src%meta = 'flex_pca_probe_polished.txt'
        else if( file_exists(trim(src%dir)//'flex_pca_merged_pc001.mrc') .and. &
               &file_exists(trim(src%dir)//'flex_pca_probe_merged.txt') )then
            src%pfx = 'flex_pca_merged_pc'; src%meta = 'flex_pca_probe_merged.txt'
        else if( file_exists(trim(src%dir)//'flex_pca_pc001.mrc') .and. &
               &file_exists(trim(src%dir)//'flex_pca_probe.txt') )then
            src%pfx = 'flex_pca_pc'; src%meta = 'flex_pca_probe.txt'
        else
            THROW_HARD('SIMPLE_COV_COMPOSE: no delivered basis + probe meta in '//trim(src%dir))
        endif
        call load_probe_state(src%ncomp, src%eigvals, src%sig2, fname=trim(src%dir)//trim(src%meta))
        nfound = 0
        do q = 1, src%ncomp
            if( .not. file_exists(trim(src%dir)//trim(src%pfx)//int2str_pad(q,3)//MRC_EXT) ) exit
            nfound = q
        end do
        if( nfound < 1 ) THROW_HARD('SIMPLE_COV_COMPOSE: no eigenvolumes under '//trim(src%dir)//trim(src%pfx))
        if( nfound < src%ncomp )then
            write(logfhandle,'(A,I0,A,I0,A,A)') '>>> FLEX_PCA COMPOSE: meta rank ', src%ncomp, &
                &' but only ', nfound, ' eigenvolumes on disk under ', trim(src%dir)//trim(src%pfx)
            src%ncomp = nfound
        endif
        fn = trim(src%dir)//trim(src%pfx)//'001'//MRC_EXT
        call find_ldim_nptcls(fn, ldim, nvol)
        call fn%kill
        src%box = ldim(1)
        if( src%box > params%box_crop ) THROW_HARD('SIMPLE_COV_COMPOSE: a source box exceeds box_crop; compose at the finest box')
    end subroutine resolve_source

    !> stable insertion sort, coarse box first
    subroutine sort_sources_by_box( srcs, nsrc )
        type(compose_src_t), intent(inout) :: srcs(:)
        integer,             intent(in)    :: nsrc
        type(compose_src_t) :: tmp
        integer :: i, j
        do i = 2, nsrc
            tmp = srcs(i)
            j   = i - 1
            do while( j >= 1 )
                if( srcs(j)%box <= tmp%box ) exit
                srcs(j+1) = srcs(j)
                j = j - 1
            end do
            srcs(j+1) = tmp
        end do
    end subroutine sort_sources_by_box

    !> Cut the composed union to its signal subspace BEFORE the deconvolution, by re-embedding:
    !! one all-N embedding with the union, the noise calibration and the noise-whitened
    !! population-variance criterion of the stage cut (SIMPLE_COV_CUT_SNR), then the basis is
    !! rotated onto the kept directions (U W+), re-orthonormalised and written back under the plain
    !! namespace so the run's embedding proper (and its distributed workers) see k columns. The
    !! deconvolution's cost grows with the square of the latent dimension, so a 58-column union
    !! without this step spends over an hour in the K ladder and its per-particle precision
    !! matrices are mostly noise dimensions.
    module subroutine compose_cut_reembed( params, build, mean_rec, basis_recs, eigvals, ncomp, sig2_eff, &
        &pinds, nptcls, rounds )
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), allocatable, intent(inout) :: basis_recs(:)
        real(dp),            allocatable, intent(inout) :: eigvals(:)
        integer,             intent(inout) :: ncomp
        real(dp),            intent(in)    :: sig2_eff
        integer,             intent(in)    :: pinds(:), nptcls
        type(image), allocatable :: imgs(:)
        type(image)  :: vol
        type(string) :: fname
        real(dp), allocatable :: z(:,:), con(:), prec(:,:,:), re(:), rme(:), zhalf(:,:,:), prior(:), a_comp(:)
        real(dp), allocatable :: R(:,:,:), Nz(:,:,:), mu(:,:), Sig(:,:,:), W(:,:), WWt(:,:), WWti(:,:), Wp(:,:)
        real(dp), allocatable :: VV(:,:), U2(:,:), zc(:,:), evc(:), qcol(:), s2(:), r2(:), nzc(:,:)
        real(dp) :: pik(1), a, snr_thr, pj
        real,     pointer :: rmat(:,:,:) => null()
        integer  :: d, k, i, j, q, n, bc, nvox, errflg, kept
        d    = ncomp
        n    = nptcls
        bc   = params%box_crop
        nvox = bc*bc*bc
        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA COMPOSE CUT: embedding the ', d, &
            &'-column union once to find its signal subspace'
        call flush(logfhandle)
        allocate(z(n,d), con(n), prec(d,d,n), re(n), rme(n), zhalf(n,d,2), prior(d), a_comp(d))
        zhalf = 0.d0
        do q = 1, d
            prior(q) = 1.d0/max(eigvals(q), DTINY)
        end do
        if( rounds%distributed() )then
            call save_probe_state(d, eigvals, sig2_eff)
            call rounds%run_stage(params, PCA_STAGE_EMBED, 'embedding')
            call embed_latents_with_contrast(params, build, mean_rec, basis_recs, d, eigvals, sig2_eff, &
                &pinds, n, z, con, prec, re, rme, from_parts=.true., rounds=rounds, &
                &zhalf_out=zhalf)
        else
            call embed_latents_with_contrast(params, build, mean_rec, basis_recs, d, eigvals, sig2_eff, &
                &pinds, n, z, con, prec, re, rme, rounds=rounds, zhalf_out=zhalf)
        endif
        deallocate(con, re, rme)
        call calibrate_noise_scale(z, zhalf, prec, prior, n, d, a, a_comp)
        deallocate(zhalf, a_comp)
        allocate(R(d,d,n), Nz(d,d,n))
        call noise_and_projection(prec, prior, n, d, a, R, Nz)
        deallocate(R, prec)
        ! method-of-moments population covariance: observed scatter minus the mean measurement noise
        allocate(mu(d,1), Sig(d,d,1))
        mu(:,1) = sum(z, dim=1)/real(n,dp)
        Sig = 0.d0
        do i = 1, n
            do q = 1, d
                Sig(:,q,1) = Sig(:,q,1) + (z(i,:) - mu(:,1))*(z(i,q) - mu(q,1))
            end do
        end do
        Sig(:,:,1) = Sig(:,:,1)/real(max(n-1,1),dp)
        do i = 1, n
            Sig(:,:,1) = Sig(:,:,1) - Nz(:,:,i)/real(n,dp)
        end do
        Sig(:,:,1) = 0.5d0*(Sig(:,:,1) + transpose(Sig(:,:,1)))
        pik(1)  = 1.d0
        snr_thr = 1.d0
        call cov_env_dp('SIMPLE_COV_CUT_SNR', snr_thr)
        call signal_subspace(mu, Sig, pik, 1, Nz, n, d, snr_thr, W, k)
        write(logfhandle,'(A,F7.3,A,I0,A,I0)') '>>> FLEX_PCA COMPOSE CUT: noise scale a=', a, &
            &'  signal subspace ', k, ' of ', d
        call flush(logfhandle)
        if( k >= d )then
            write(logfhandle,'(A)') '>>> FLEX_PCA COMPOSE CUT: every direction carries population variance; no cut'
            deallocate(z, prior, Nz, mu, Sig, W)
            return
        endif
        ! cut coordinates and their population variance (observed minus noise along each direction)
        allocate(zc(n,k), evc(k), nzc(k,k))
        zc  = matmul(z, transpose(W))
        nzc = 0.d0
        do i = 1, n
            nzc = nzc + matmul(W, matmul(Nz(:,:,i), transpose(W)))/real(n,dp)
        end do
        do q = 1, k
            evc(q) = sum((zc(:,q) - sum(zc(:,q))/real(n,dp))**2)/real(max(n-1,1),dp)
            evc(q) = max(evc(q) - nzc(q,q), 0.1d0*evc(q), DTINY)
        end do
        deallocate(z, Nz, zc, nzc, mu, Sig, prior)
        ! U W+ : the composed columns combined onto the kept directions (z = W+ zc on the signal subspace)
        allocate(WWt(k,k), WWti(k,k), Wp(d,k))
        WWt = matmul(W, transpose(W))
        call matinv(WWt, WWti, k, errflg)
        if( errflg /= 0 ) THROW_HARD('compose_cut_reembed: singular W W^T')
        Wp = matmul(transpose(W), WWti)
        allocate(VV(nvox,d), U2(nvox,k))
        call vol%new([bc,bc,bc], params%smpd_crop)
        do j = 1, d
            fname = 'flex_pca_pc'//int2str_pad(j,3)//MRC_EXT
            call vol%read(fname)
            call fname%kill
            call vol%get_rmat_ptr(rmat)
            VV(:,j) = reshape(real(rmat(1:bc,1:bc,1:bc),dp), [nvox])
        end do
        U2 = matmul(VV, Wp)
        deallocate(VV, W, WWt, WWti, Wp)
        ! re-orthonormalise (the kept directions are orthogonal in the whitened frame, not in real space);
        ! the prior variance follows the rescaling as in the composition
        allocate(qcol(nvox), s2(k), r2(k))
        kept = 0
        do j = 1, k
            s2(j) = sum(U2(:,j)**2)
            if( s2(j) <= 1.d-30 ) cycle
            qcol = U2(:,j)/sqrt(s2(j))
            do q = 1, kept
                pj   = dot_product(qcol, U2(:,q))
                qcol = qcol - pj*U2(:,q)
            end do
            r2(j) = sum(qcol**2)
            if( r2(j) < 1.d-6 ) cycle
            kept = kept + 1
            U2(:,kept) = qcol/sqrt(r2(j))
            evc(kept)  = evc(j)*s2(j)*r2(j)
        end do
        deallocate(qcol)
        if( kept < 2 ) THROW_HARD('compose_cut_reembed: fewer than 2 columns survived the cut')
        ! realise: images -> reconstructors, plain namespace rewritten (k files; the surplus removed)
        allocate(imgs(kept))
        do j = 1, kept
            call imgs(j)%new([bc,bc,bc], params%smpd_crop)
            call imgs(j)%get_rmat_ptr(rmat)
            rmat(1:bc,1:bc,1:bc) = real(reshape(U2(:,j), [bc,bc,bc]))
            fname = 'flex_pca_pc'//int2str_pad(j,3)//MRC_EXT
            call imgs(j)%write(fname, del_if_exists=.true.)
            call fname%kill
        end do
        do j = kept+1, d
            fname = 'flex_pca_pc'//int2str_pad(j,3)//MRC_EXT
            call del_file(fname)
            call fname%kill
        end do
        deallocate(U2)
        do q = 1, size(basis_recs)
            call basis_recs(q)%dealloc_rho; call basis_recs(q)%kill
        end do
        deallocate(basis_recs, eigvals)
        call basis_recs_from_images(params, build, imgs, kept, basis_recs)
        do j = 1, kept
            call imgs(j)%kill
        end do
        deallocate(imgs)
        allocate(eigvals(kept))
        eigvals = evc(1:kept)
        ncomp   = kept
        call save_probe_state(ncomp, eigvals, sig2_eff)
        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA COMPOSE CUT: ', ncomp, ' columns re-embedded; prior variances:'
        write(logfhandle,'(A,12(1X,ES9.2))') '>>>   ', (eigvals(q), q=1,min(12,ncomp))
        call flush(logfhandle)
        deallocate(evc, s2, r2)
        call vol%kill
    end subroutine compose_cut_reembed

end submodule simple_flex_pca_em_compose
