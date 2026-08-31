!@descr: flex_pca EM: eigenvolume realization, orthonormalization, alignment and bagging
submodule (simple_flex_pca_em) simple_flex_pca_em_basis
implicit none
#include "simple_local_flags.inc"

contains









    !> Convert each merged complex column C_q into its two real spatial representatives Re(ifft
    !! C_q)=Sigma*cos_q and Im(ifft C_q)=Sigma*sin_q.
    module subroutine basis_to_real_representatives( params, work, colvol, ncol, lb, ub, realvols, nreal )
        class(parameters),   intent(inout) :: params
        type(reconstructor), intent(inout) :: work
        complex,             intent(in)    :: colvol(:,:,:,:)
        integer,             intent(in)    :: ncol, lb(3), ub(3)
        type(image), allocatable, intent(out) :: realvols(:)
        integer,                  intent(out) :: nreal
        type(image)  :: gridcorr_img
        complex, allocatable :: vr(:,:,:), vi(:,:,:)
        integer :: s, i1, i2, i3, n1, n2, n3, hn, kn, ln, h, k, l
        real    :: energy
        n1 = ub(1)-lb(1)+1; n2 = ub(2)-lb(2)+1; n3 = ub(3)-lb(3)+1
        allocate(vr(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
        allocate(vi(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
        gridcorr_img = prep3D_inv_kbenvelope4mul([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        allocate(realvols(2*ncol))
        nreal = 0
        do s = 1, ncol
            do i3 = 1, n3
                l = lb(3)+i3-1; ln = -l
                do i2 = 1, n2
                    k = lb(2)+i2-1; kn = -k
                    do i1 = 1, n1
                        h = lb(1)+i1-1; hn = -h
                        if( hn < lb(1) .or. hn > ub(1) .or. kn < lb(2) .or. kn > ub(2) &
                            &.or. ln < lb(3) .or. ln > ub(3) )then
                            vr(h,k,l) = 0.5*colvol(i1,i2,i3,s)
                            vi(h,k,l) = cmplx(0.,-0.5)*colvol(i1,i2,i3,s)
                        else
                            vr(h,k,l) = 0.5*(colvol(i1,i2,i3,s) + conjg(colvol(n1-i1+1,n2-i2+1,n3-i3+1,s)))
                            vi(h,k,l) = cmplx(0.,-0.5)*(colvol(i1,i2,i3,s) - conjg(colvol(n1-i1+1,n2-i2+1,n3-i3+1,s)))
                        endif
                    end do
                end do
            end do
            call realize_hermitian_volume(params, work, vr, gridcorr_img, energy)
            if( energy > 0. )then
                nreal = nreal + 1
                call realvols(nreal)%copy(work)
            endif
            call realize_hermitian_volume(params, work, vi, gridcorr_img, energy)
            if( energy > 0. )then
                nreal = nreal + 1
                call realvols(nreal)%copy(work)
            endif
        end do
        call gridcorr_img%kill
        deallocate(vr, vi)
    end subroutine basis_to_real_representatives

    !>  Load a Hermitian expanded Fourier volume into the work reconstructor, fold to
    !!  compressed storage, inverse-FFT to a real volume, deapodize, low-pass and mask.
    module subroutine realize_hermitian_volume( params, work, vherm, gridcorr_img, energy )
        class(parameters),   intent(in)    :: params
        type(reconstructor), intent(inout) :: work
        complex,             intent(in)    :: vherm(:,:,:)
        type(image),         intent(inout) :: gridcorr_img
        real,                intent(out)   :: energy
        real, pointer :: rmat(:,:,:)
        integer       :: ldim_work(3)
        call work%reset
        call work%reset_exp
        work%cmat_exp = vherm
        call work%compress_exp
        ! Band-limit the covariance column to the signal band FIRST, in Fourier space, before any real-space
        ! operation, so out-of-band shells never enter the masking or the Gram products downstream.
        if( params%lp > 2.0*params%smpd_crop + TINY ) call work%bp(0., params%lp)
        call work%ifft
        ! deapodize on the native lattice (same correction as production gridding)
        call work%mul(gridcorr_img)
        if( params%msk_crop > TINY ) call work%mask3D_soft(params%msk_crop, backgr=0.)
        if( work%is_ft() ) call work%ifft
        call work%get_rmat_ptr(rmat)
        ldim_work = work%get_ldim()
        energy = sum(rmat(1:ldim_work(1),1:ldim_work(2),1:ldim_work(3))**2)
    end subroutine realize_hermitian_volume

    !> Orthonormalize the real column representatives into the column subspace Utilde by Gram
    !! eigendecomposition, keeping every direction above a relative energy floor.
    module subroutine orthonormalize_representatives( params, build, realvols, nreal, utilde, utilde_real, d_tilde, svals, &
        &nptcls_basis )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(image),         intent(inout) :: realvols(:)
        integer,             intent(in)    :: nreal
        type(reconstructor), allocatable, intent(out) :: utilde(:)
        type(image),         allocatable, intent(out) :: utilde_real(:)
        integer,             intent(out)   :: d_tilde
        !> squared singular values of the representative set.
        real(dp), allocatable, optional, intent(out) :: svals(:)
        !> particles that actually reached the basis stages. Only used to REPORT the
        !! samples-per-parameter rank bound; it selects nothing.
        integer, optional, intent(in) :: nptcls_basis
        real(dp), allocatable :: gram(:,:), evec(:,:), eval(:)
        real, pointer :: rmat_i(:,:,:), rmat_j(:,:,:)
        integer :: i, q, nrot, keep, d_budget, d_cap, d_signal, d_samples, nbasis
        real(dp) :: lam_max, nrm
        logical  :: l_packed
        character(len=9) :: accum_model
        if( nreal < 1 ) THROW_HARD('flex_pca produced no covariance column representatives')
        allocate(gram(nreal,nreal), evec(nreal,nreal), eval(nreal))
        do i = 1, nreal
            call realvols(i)%get_rmat_ptr(rmat_i)
            do q = i, nreal
                call realvols(q)%get_rmat_ptr(rmat_j)
                gram(i,q) = sum(real(rmat_i,dp)*real(rmat_j,dp))
                gram(q,i) = gram(i,q)
            end do
        end do
        call jacobi(gram, nreal, nreal, eval, evec, nrot)
        call eigsrt(eval, evec, nreal, nreal)              ! descending
        lam_max = max(eval(1), DTINY)
        keep = 0
        do q = 1, nreal
            if( eval(q) > COV_EIG_REL_FLOOR*lam_max ) keep = keep + 1
        end do
        ! Largest d the reduced solve's accumulator can afford: one shared d^2 x d^2 array -> 8*d^4 bytes,
        ! against 8*[d(d+1)/2]^2 ~ 2*d^4 bytes for the packed CG path, which never forms the operator.
        ! Size it against the model the solve will ACTUALLY use -- sizing d for the dense accumulator and
        ! then solving packed spends a quarter of the budget and caps the column subspace 41 % below what
        ! the data supports (at 8 GB: d 177 instead of 250).
        l_packed = cov_packed_cgsolve()
        d_budget = cov_dim_budget(l_packed)
        ! data-driven rank, REPORT ONLY: the energy floor and the memory budget never ask how many
        ! directions are real, so log what the data would support and let the discrepancy show
        d_signal = cov_signal_rank(eval, nreal)
        nbasis   = 0
        if( present(nptcls_basis) ) nbasis = nptcls_basis
        d_samples = 0
        if( nbasis > 0 ) d_samples = &
            &max(1, int((-1.d0 + sqrt(1.d0 + 8.d0*real(nbasis,dp)/COV_SAMPLES_PER_PARAM))/2.d0))
        ! memory budget is a GUARD, not the chooser. SIMPLE_COV_DTILDE replaces both, so a fixed-d A/B
        ! is never silently clamped by the box's RAM.
        d_cap = min(d_budget, COV_DEFAULT_DTILDE)
        call cov_env_int('SIMPLE_COV_DTILDE', d_cap)
        d_tilde  = max(1, min(keep, COV_MAX_DTILDE, d_cap))
        if( cov_env_int_on('SIMPLE_COV_DSIGNAL') ) d_tilde = max(1, min(d_tilde, d_signal))
        if( d_samples > 0 )then
            write(logfhandle,'(A,I0,A,I0,A,F4.1,A)') '>>> FLEX_PCA d_samples=',d_samples, &
                &'  (samples-per-parameter bound from N=',nbasis,' at R=',COV_SAMPLES_PER_PARAM, &
                &') -- REPORT ONLY'
        endif
        write(logfhandle,'(A,I0,A,A,A)') '>>> FLEX_PCA d_signal=',d_signal, &
            &' (spectrum noise-bulk estimate; ', &
            &trim(merge('BINDING ','reported',cov_env_int_on('SIMPLE_COV_DSIGNAL'))), &
            &' -- set SIMPLE_COV_DSIGNAL=1 to bind it)'
        if( l_packed )then
            accum_model = 'packed+CG'
        else
            accum_model = 'dense'
        endif
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA d_tilde=',d_tilde, &
            &'  (above energy floor=',keep,', memory cap=',d_budget,', rank cap=',COV_MAX_DTILDE, &
            &', default=',COV_DEFAULT_DTILDE,')'
        write(logfhandle,'(A,A,A,F8.3,A,F6.3,A)') '>>> FLEX_PCA reduced-solve accumulator model: ', &
            &trim(accum_model),', ',cov_accum_bytes(d_tilde, l_packed)/1.d9, &
            &' GB at this d_tilde (budget ',COV_ATHR_BUDGET/1.d9,' GB)'
        if( d_tilde == d_budget .and. keep > d_budget )then
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA NOTE: the column subspace is limited by the &
                &reduced-solve memory budget, not by the data; ',keep,' directions cleared the energy floor.'
        else if( d_tilde == COV_DEFAULT_DTILDE .and. d_budget > COV_DEFAULT_DTILDE )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA NOTE: d_tilde is the measured default; the &
                &memory budget would have allowed ',d_budget,'. Override with SIMPLE_COV_DTILDE=',d_budget, &
                &' to restore the budget-chosen rank.'
        endif
        call flush(logfhandle)
        allocate(utilde(d_tilde), utilde_real(d_tilde))
        do q = 1, d_tilde
            nrm = sqrt(max(eval(q), DTINY))
            ! Build the unit-norm orthonormal basis vector as a CLEAN plain image via image arithmetic on
            ! the (verified band-limited) representatives, rather than raw rmat-pointer math on padded
            ! reconstructor buffers.
            call utilde_real(q)%copy(realvols(1))
            call utilde_real(q)%zero_and_unflag_ft
            do i = 1, nreal
                call utilde_real(q)%add(realvols(i), real(evec(i,q)/nrm))
            end do
            ! set_rmat(...,.false.) then fft then expand_exp -- NEVER add(), which leaves the reconstructor
            ! flagged as Fourier and silently propagates an untransformed grid (see form_eigenbasis_from_reduced)
            call init_basis_reconstructor(params, build, utilde(q))
            call utilde(q)%set_rmat(utilde_real(q)%get_rmat(), .false.)
            call utilde(q)%fft
            call utilde(q)%expand_exp
        end do
        if( present(svals) )then
            allocate(svals(d_tilde))
            do q = 1, d_tilde
                svals(q) = max(eval(q), DTINY)
            end do
        endif
        deallocate(gram, evec, eval)
    end subroutine orthonormalize_representatives



    !>  Cross-halfset basis alignment for the held-out embedding. Both eigenbases are
    !!  first normalized to unit real-space norm, then M(i,j) = <U_ref_i, U_tgt_j>.
    !!  A latent expressed in the TARGET basis is mapped into the REFERENCE frame by
    !!  z_ref = M z_tgt (since x ~ U_tgt z_tgt and z_ref = U_ref^T x). The singular
    !!  values of M are the cosines of the principal angles between the two subspaces,
    !!  i.e. a gold-standard measure of how many latent dimensions actually reproduce
    !!  across independent halves -- unlike per-component FSC, this cannot be fooled by
    !!  a shared basis, because the two bases are estimated from disjoint particles.
    module subroutine align_basis_to_reference( ref_imgs, nref_c, tgt_imgs, ntgt_c, M, svals )
        integer,     intent(in)    :: nref_c, ntgt_c
        type(image), intent(inout) :: ref_imgs(nref_c), tgt_imgs(ntgt_c)
        real(dp), allocatable, intent(out) :: M(:,:), svals(:)
        real, pointer :: rmat_i(:,:,:), rmat_j(:,:,:)
        real(dp), allocatable :: nrm_r(:), nrm_t(:), Mwork(:,:), V(:,:), ev(:)
        integer  :: i, j, nrot, nsv
        allocate(M(nref_c,ntgt_c), source=0.d0)
        allocate(nrm_r(nref_c), nrm_t(ntgt_c), source=0.d0)
        do i = 1, nref_c
            call ref_imgs(i)%get_rmat_ptr(rmat_i)
            nrm_r(i) = sqrt(max(sum(real(rmat_i,dp)*real(rmat_i,dp)), DTINY))
        end do
        do j = 1, ntgt_c
            call tgt_imgs(j)%get_rmat_ptr(rmat_j)
            nrm_t(j) = sqrt(max(sum(real(rmat_j,dp)*real(rmat_j,dp)), DTINY))
        end do
        do i = 1, nref_c
            call ref_imgs(i)%get_rmat_ptr(rmat_i)
            do j = 1, ntgt_c
                call tgt_imgs(j)%get_rmat_ptr(rmat_j)
                M(i,j) = sum(real(rmat_i,dp)*real(rmat_j,dp)) / (nrm_r(i)*nrm_t(j))
            end do
        end do
        ! principal-angle cosines = singular values of M, via the eigenvalues of M^T M
        nsv = min(nref_c, ntgt_c)
        allocate(Mwork(ntgt_c,ntgt_c), V(ntgt_c,ntgt_c), ev(ntgt_c), svals(nsv))
        Mwork = matmul(transpose(M), M)
        call jacobi(Mwork, ntgt_c, ntgt_c, ev, V, nrot)
        call eigsrt(ev, V, ntgt_c, ntgt_c)
        do i = 1, nsv
            svals(i) = sqrt(max(ev(i), 0.d0))
        end do
        deallocate(nrm_r, nrm_t, Mwork, V, ev)
    end subroutine align_basis_to_reference



    !> Bagged basis: pool two independently fitted eigenbases, keep the top-k singular directions.
    !!
    !!  Nominally equivalent fits of the same data vary 5x in explanatory power, so the variance is
    !!  in the fit, not the sample. Directions the two fits AGREE on add coherently in the pooled
    !!  Gram and survive truncation; directions either invented on its own do not. Ceiling is the
    !!  basis-error share of latent error, ~38 % -- the rest is per-particle estimation noise, which
    !!  the deconvolution in simple_flex_pca_model::gmm_state_weights targets instead.
    !!
    !!  Inputs are compared and combined over the WHOLE rmat, matching align_basis_to_reference.
    module subroutine bag_basis_pool( imgs_a, na_c, eig_a, imgs_b, nb_c, eig_b, ncomp_out, pooled, eig_pooled )
        integer,     intent(in)    :: na_c, nb_c, ncomp_out
        type(image), intent(inout) :: imgs_a(na_c), imgs_b(nb_c)
        real(dp),    intent(in)    :: eig_a(na_c), eig_b(nb_c)
        type(image), allocatable, intent(out) :: pooled(:)
        real(dp),    allocatable, intent(out) :: eig_pooled(:)
        real, pointer :: rmat_i(:,:,:), rmat_j(:,:,:), rmat_o(:,:,:)
        real(dp), allocatable :: G(:,:), V(:,:), ev(:), nrm(:), eigin(:)
        integer  :: m, k, i, j, q, nrot, ldim(3)
        real(dp) :: onrm
        real     :: smpd
        m = na_c + nb_c
        k = min(ncomp_out, m)
        if( k < 1 ) THROW_HARD('flex_pca bagging asked for a rank below 1')
        allocate(nrm(m), eigin(m), G(m,m), V(m,m), ev(m))
        do j = 1, m
            call pick(j, rmat_j)
            nrm(j)   = sqrt(max(sum(real(rmat_j,dp)*real(rmat_j,dp)), DTINY))
            eigin(j) = merge(eig_a(min(j,na_c)), eig_b(max(1,j-na_c)), j <= na_c)
        end do
        ! Gram of the unit-normalised union: equal weight per direction is what makes this a bag,
        ! since eigenvalue weighting would reinstate a ranking measured not to track signal.
        do i = 1, m
            call pick(i, rmat_i)
            do j = i, m
                call pick(j, rmat_j)
                G(i,j) = sum(real(rmat_i,dp)*real(rmat_j,dp)) / (nrm(i)*nrm(j))
                G(j,i) = G(i,j)
            end do
        end do
        call jacobi(G, m, m, ev, V, nrot)
        call eigsrt(ev, V, m, m)
        allocate(pooled(k), eig_pooled(k))
        ldim = imgs_a(1)%get_ldim()
        smpd = imgs_a(1)%get_smpd()
        do q = 1, k
            call pooled(q)%new(ldim, smpd)
            call pooled(q)%get_rmat_ptr(rmat_o)
            rmat_o = 0.
            do j = 1, m
                call pick(j, rmat_j)
                rmat_o = rmat_o + real(V(j,q)/nrm(j)) * rmat_j
            end do
            onrm = sqrt(max(sum(real(rmat_o,dp)*real(rmat_o,dp)), DTINY))
            rmat_o = rmat_o / real(onrm)
            ! variance along the pooled direction, propagated from the input spectra: the MAP prior
            ! needs the inputs' units, and the Gram eigenvalue is an overlap, not a variance.
            eig_pooled(q) = sum(V(:,q)*V(:,q)*eigin(:))
        end do
        deallocate(nrm, eigin, G, V, ev)
      contains
        !> leading na_c indices address basis A, the remainder basis B
        subroutine pick( idx, p )
            integer,       intent(in)  :: idx
            real, pointer, intent(out) :: p(:,:,:)
            if( idx <= na_c )then
                call imgs_a(idx)%get_rmat_ptr(p)
            else
                call imgs_b(idx-na_c)%get_rmat_ptr(p)
            endif
        end subroutine pick
    end subroutine bag_basis_pool

    !> Turn a set of real-space basis volumes into embedding-ready column reconstructors.
    module subroutine basis_recs_from_images( params, build, imgs, ncomp, basis_recs )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        integer,             intent(in)    :: ncomp
        type(image),         intent(inout) :: imgs(ncomp)
        type(reconstructor), allocatable, intent(out) :: basis_recs(:)
        integer :: q
        allocate(basis_recs(ncomp))
        do q = 1, ncomp
            call init_basis_reconstructor(params, build, basis_recs(q))
            call basis_recs(q)%set_rmat(imgs(q)%get_rmat(), .false.)
            call basis_recs(q)%fft
            call basis_recs(q)%expand_exp
        end do
    end subroutine basis_recs_from_images

end submodule simple_flex_pca_em_basis
