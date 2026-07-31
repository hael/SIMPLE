!@descr: projection-aware residual Nyström pre-images for flex_analysis
module simple_flex_diffmap_rec3D
use simple_core_module_api
use simple_builder,                only: builder
use simple_sp_project,             only: sp_project
use simple_gridding,               only: prep3D_inv_instrfun4mul
use simple_image,                  only: image
use simple_matcher_3Drec,          only: init_rec, prep_imgs4rec, cleanup_rec_buffers
use simple_matcher_ptcl_io,        only: discrete_read_imgbatch, discrete_read_imgbatch_source, prepimgbatch
use simple_parameters,             only: parameters
use simple_flex_projected_latent_model, only: update_basis_from_latents, write_mstep_stats_part_file, &
    &update_basis_from_mstep_stats_part_files, cleanup_planes, &
    &prep_imgs4projected_model, solve_coupled_basis_exp
use simple_flex_reconstructor_latent_ops, only: project_fplane_mean
use simple_flex_reconstructor_latent_ops, only: insert_planes_oversamp_coupled_batch_scaled
use simple_flex_reconstructor_latent_ops, only: insert_plane_oversamp_multi_scaled
use simple_reconstructor,          only: reconstructor
implicit none
private
#include "simple_local_flags.inc"

public :: reconstruct_flex_diffmap_states, write_flex_diffmap_rec_parts, reduce_flex_diffmap_rec_parts
public :: reconstruct_flex_diffmap_weighted_states
public :: reconstruct_flex_diffmap_local_linear_states
public :: cleanup_flex_diffmap_rec_parts
public :: test_flex_local_linear_preimage

contains

    !> Self-contained numerical validation of the local-linear pre-image
    !! construction, exercising the production per-voxel solver
    !! (solve_coupled_basis_exp) on synthetic normal equations assembled exactly
    !! as reconstruct_flex_diffmap_local_linear_states does (A = sum w_i g_i g_i^T,
    !! b = sum w_i g_i D_i, g_i = [1, u_i]).  No project or particle data required.
    !!  T1: degenerate design (all u_i=0) recovers the Nadaraya-Watson average.
    !!  T2: local-linear intercept removes the O(h) sampling-design bias that
    !!      biases the local-constant (NW) estimator; recovers a linear response
    !!      exactly and beats NW on a curved response.
    !!  T5: a rank-deficient (collinear) design stays finite through the ridge/
    !!      fallback and reproduces the constant data value; never NaN.
    subroutine test_flex_local_linear_preimage()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        integer, parameter :: NP=40
        real(dp) :: u(NP), w(NP), aa, bb, cc, ubar, wsum, h
        complex(dp) :: dvals(NP), sol_nw, sol_ll, sol_ll2, pred
        real(dp) :: nw_err, ll_err, nw_err_q, ll_err_q
        integer  :: i
        h  = 0.5d0
        aa = 2.0d0
        bb = 1.5d0
        cc = 0.8d0
        ! One-sided (boundary) sampling on u in [0,1]; the target pre-image is at
        ! u=0.  Gaussian kernel weights make the weighted mean of u strictly
        ! positive, which is exactly the regime where the local-constant estimator
        ! is biased and the local-linear estimator is not.
        do i=1,NP
            u(i) = real(i-1,dp)/real(NP-1,dp)
            w(i) = exp(-u(i)*u(i)/(h*h))
        end do
        wsum = sum(w)
        ubar = sum(w*u)/wsum
        ! ---- T2a: pure linear response D = a + b*u ----
        do i=1,NP
            dvals(i) = cmplx(aa + bb*u(i), 0.d0, kind=dp)
        end do
        sol_nw = ll_solve_intercept(0, u, w, dvals)   ! d=0 -> Nadaraya-Watson
        sol_ll = ll_solve_intercept(1, u, w, dvals)   ! d=1 -> local-linear
        nw_err = abs(sol_nw - cmplx(aa,0.d0,kind=dp))
        ll_err = abs(sol_ll - cmplx(aa,0.d0,kind=dp))
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,F7.4)') '>>> FLEX LOCAL-LINEAR T2a nw_err=',nw_err, &
            &' ll_err=',ll_err,' weighted_mean_u=',ubar
        if( .not.ieee_is_finite(real(sol_ll)) ) THROW_HARD('local-linear intercept is not finite (T2a)')
        if( ll_err > 1.d-4 ) THROW_HARD('local-linear intercept fails to recover a linear response (T2a)')
        if( nw_err < 20.d0*max(ll_err,1.d-12) ) &
            &THROW_HARD('local-constant estimator is not measurably biased vs local-linear on one-sided sampling (T2a)')
        ! ---- T2b: curved response D = a + b*u + c*u^2 ----
        do i=1,NP
            dvals(i) = cmplx(aa + bb*u(i) + cc*u(i)*u(i), 0.d0, kind=dp)
        end do
        sol_nw = ll_solve_intercept(0, u, w, dvals)
        sol_ll = ll_solve_intercept(1, u, w, dvals)
        nw_err_q = abs(sol_nw - cmplx(aa,0.d0,kind=dp))
        ll_err_q = abs(sol_ll - cmplx(aa,0.d0,kind=dp))
        write(logfhandle,'(A,ES12.4,A,ES12.4)') '>>> FLEX LOCAL-LINEAR T2b nw_err=',nw_err_q,' ll_err=',ll_err_q
        if( ll_err_q > 0.3d0*nw_err_q ) &
            &THROW_HARD('local-linear intercept does not improve on local-constant for a curved response (T2b)')
        ! ---- T1: degenerate design, all u_i = 0 -> intercept == NW average ----
        do i=1,NP
            u(i)     = 0.d0
            w(i)     = 0.5d0 + 0.5d0*real(i,dp)/real(NP,dp)
            dvals(i) = cmplx(aa + 0.1d0*real(i,dp), -0.05d0*real(i,dp), kind=dp)
        end do
        sol_nw = ll_solve_intercept(0, u, w, dvals)
        sol_ll = ll_solve_intercept(1, u, w, dvals)
        write(logfhandle,'(A,2ES12.4,A,2ES12.4)') '>>> FLEX LOCAL-LINEAR T1 nw=',sol_nw,' ll=',sol_ll
        if( .not.ieee_is_finite(real(sol_ll)) .or. .not.ieee_is_finite(aimag(sol_ll)) ) &
            &THROW_HARD('degenerate local-linear intercept is not finite (T1)')
        if( abs(sol_ll - sol_nw) > 1.d-4 ) &
            &THROW_HARD('degenerate local-linear intercept differs from Nadaraya-Watson average (T1)')
        ! ---- T5: rank-deficient collinear design (all u_i equal) stays finite ----
        do i=1,NP
            u(i)     = 0.5d0
            w(i)     = 0.4d0 + 0.6d0*real(i,dp)/real(NP,dp)
            dvals(i) = cmplx(aa + bb*0.5d0, 0.d0, kind=dp)
        end do
        sol_ll  = ll_solve_intercept(1, u, w, dvals)
        sol_ll2 = ll_solve_slope(1, u, w, dvals)
        pred    = sol_ll + sol_ll2*cmplx(0.5d0,0.d0,kind=dp)
        write(logfhandle,'(A,2ES12.4,A,2ES12.4)') '>>> FLEX LOCAL-LINEAR T5 intercept=',sol_ll,' predicted@u0=',pred
        if( .not.ieee_is_finite(real(sol_ll)) .or. .not.ieee_is_finite(aimag(sol_ll)) .or. &
            &.not.ieee_is_finite(real(sol_ll2)) .or. .not.ieee_is_finite(aimag(sol_ll2)) ) &
            &THROW_HARD('rank-deficient local-linear solve produced a non-finite value (T5)')
        if( abs(pred - cmplx(aa + bb*0.5d0, 0.d0, kind=dp)) > 1.d-2 ) &
            &THROW_HARD('rank-deficient local-linear solve does not reproduce the constant data value (T5)')
        write(logfhandle,'(A)') '>>> FLEX LOCAL-LINEAR PRE-IMAGE TEST PASSED'

    contains

        !> Assemble the local-linear normal equations for a single Fourier voxel
        !! and return the intercept (component 1) via the production solver.
        !! ndim=0 yields the local-constant (Nadaraya-Watson) weighted average.
        function ll_solve_intercept( ndim, u, w, dvals ) result( sol1 )
            integer,     intent(in) :: ndim
            real(dp),    intent(in) :: u(:), w(:)
            complex(dp), intent(in) :: dvals(:)
            complex(dp) :: sol1, sol(ndim+1)
            call ll_solve( ndim, u, w, dvals, sol )
            sol1 = sol(1)
        end function ll_solve_intercept

        function ll_solve_slope( ndim, u, w, dvals ) result( slope )
            integer,     intent(in) :: ndim
            real(dp),    intent(in) :: u(:), w(:)
            complex(dp), intent(in) :: dvals(:)
            complex(dp) :: slope, sol(ndim+1)
            call ll_solve( ndim, u, w, dvals, sol )
            slope = sol(min(2,ndim+1))
        end function ll_solve_slope

        subroutine ll_solve( ndim, u, w, dvals, sol )
            integer,     intent(in)  :: ndim
            real(dp),    intent(in)  :: u(:), w(:)
            complex(dp), intent(in)  :: dvals(:)
            complex(dp), intent(out) :: sol(ndim+1)
            integer, parameter :: BOX=8
            type(reconstructor), allocatable :: recs(:)
            real, allocatable :: rho_cross(:,:,:,:)
            real(dp) :: gvec(ndim+1), amat(ndim+1,ndim+1)
            complex(dp) :: bvec(ndim+1)
            integer :: lb(3), ub(3), n, np, q, r, i, ipair, ih, ik, im
            n  = ndim + 1
            np = (n*(n+1))/2
            amat = 0.d0
            bvec = DCMPLX_ZERO
            do i=1,size(u)
                gvec(1) = 1.d0
                do q=1,ndim
                    gvec(q+1) = u(i)
                end do
                do q=1,n
                    bvec(q) = bvec(q) + w(i)*gvec(q)*dvals(i)
                    do r=1,n
                        amat(q,r) = amat(q,r) + w(i)*gvec(q)*gvec(r)
                    end do
                end do
            end do
            lb = [-BOX/2-2,-BOX/2-2,-BOX/2-2]
            ub = [ BOX/2+2, BOX/2+2, BOX/2+2]
            allocate(recs(n))
            do q=1,n
                call recs(q)%new([BOX,BOX,BOX],1.)
                allocate(recs(q)%cmat_exp(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)),source=CMPLX_ZERO)
            end do
            allocate(rho_cross(np,ub(1)-lb(1)+1,ub(2)-lb(2)+1,ub(3)-lb(3)+1),source=0.)
            ih = 1 - lb(1) + 1
            ik = 0 - lb(2) + 1
            im = 0 - lb(3) + 1
            do q=1,n
                recs(q)%cmat_exp(1,0,0) = cmplx(real(bvec(q),sp), real(aimag(bvec(q)),sp))
            end do
            do r=1,n
                do q=1,r
                    ipair = (r*(r-1))/2 + q
                    rho_cross(ipair,ih,ik,im) = real(amat(q,r))
                end do
            end do
            call solve_coupled_basis_exp(recs, rho_cross, n)
            do q=1,n
                sol(q) = cmplx(recs(q)%cmat_exp(1,0,0), kind=dp)
            end do
            do q=1,n
                call recs(q)%kill
            end do
            deallocate(recs, rho_cross)
        end subroutine ll_solve

    end subroutine test_flex_local_linear_preimage

    subroutine reconstruct_flex_diffmap_states( params, build, pinds, z, target_coeffs, nstates )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        integer, intent(in) :: pinds(:), nstates
        real, intent(in) :: z(:,:), target_coeffs(:,:)
        type(parameters) :: model_params
        type(reconstructor) :: mean_rec
        type(reconstructor), allocatable :: basis_recs(:)
        type(fplane_type), allocatable :: fpls(:)
        real(dp), allocatable :: z_dp(:,:)
        integer :: ncomp, q
        call validate_model_tables(pinds,z,target_coeffs,nstates,ncomp)
        call prepare_unfiltered_model_params(params,model_params)
        allocate(z_dp(size(z,1),ncomp),source=real(z,dp))
        allocate(basis_recs(ncomp))
        call init_mean_reconstructor(model_params,build,mean_rec)
        do q=1,ncomp
            call init_basis_reconstructor(model_params,build,basis_recs(q))
        end do
        call update_basis_from_latents(model_params,build,mean_rec,basis_recs,z_dp,pinds,size(pinds),ncomp,fpls, &
            &log_label='FLEX_DIFFMAP')
        call write_preimage_states(params,basis_recs,target_coeffs,nstates,ncomp)
        call cleanup_planes(fpls)
        call cleanup_reconstructors(mean_rec,basis_recs)
        deallocate(z_dp)
    end subroutine reconstruct_flex_diffmap_states

    !> Direct manifold pre-image reconstruction with kernel weights.
    !! Medoid-selected manifold descriptors are converted to soft particle
    !! weights upstream; this routine reconstructs each state from weighted
    !! contributions of all particles rather than from a global residual basis.
    subroutine reconstruct_flex_diffmap_weighted_states( params, build, pinds, state_weights, nstates, fsc_projfile )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:), nstates
        real,              intent(in)    :: state_weights(:,:)
        type(string), optional, intent(in) :: fsc_projfile
        type(reconstructor), allocatable :: state_recs(:)
        type(fplane_type), allocatable :: fpls(:)
        type(image) :: gridcorr_img, state_img
        type(ori) :: orientation
        real(dp), allocatable :: scales(:)
        real, allocatable :: lowpass_filters(:,:)
        logical, allocatable :: has_lowpass_filter(:)
        integer, allocatable :: lowpass_source_state(:)
        integer :: batchlims(2), batchsz, ibatch, i, iptcl, state
        if( size(pinds)<1 .or. nstates<1 ) THROW_HARD('invalid flex weighted state reconstruction dimensions')
        if( any(shape(state_weights)/=[size(pinds),nstates]) ) THROW_HARD('flex weighted state table mismatch')
        allocate(state_recs(nstates),scales(nstates))
        call prepare_project_fsc_lowpass_filters(params,build,nstates,lowpass_filters,has_lowpass_filter,lowpass_source_state, &
            &fsc_projfile)
        do state=1,nstates
            call init_basis_reconstructor(params,build,state_recs(state))
        end do
        call init_rec(params,build,MAXIMGBATCHSZ,fpls)
        call prepimgbatch(params,build,MAXIMGBATCHSZ)
        do ibatch=1,size(pinds),MAXIMGBATCHSZ
            batchlims=[ibatch,min(size(pinds),ibatch+MAXIMGBATCHSZ-1)]
            batchsz=batchlims(2)-batchlims(1)+1
            if( params%l_ptcl_src_den )then
                call discrete_read_imgbatch_source(params,build,'den',batchsz, &
                    &pinds(batchlims(1):batchlims(2)),[1,batchsz],build%imgbatch(:batchsz))
            else
                call discrete_read_imgbatch(params,build,size(pinds),pinds,batchlims)
            endif
            call prep_imgs4rec(params,build,batchsz,build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)),fpls(:batchsz))
            do i=1,batchsz
                iptcl=pinds(batchlims(1)+i-1)
                call build%spproj_field%get_ori(iptcl,orientation)
                if( orientation%isstatezero() ) cycle
                scales=real(state_weights(batchlims(1)+i-1,:),dp)
                call insert_plane_oversamp_multi_scaled(state_recs,build%pgrpsyms,orientation,fpls(i),scales,scales)
            end do
        end do
        call orientation%kill
        call cleanup_rec_buffers(build,fpls)
        gridcorr_img=prep3D_inv_instrfun4mul([params%box_crop,params%box_crop,params%box_crop], &
            &OSMPL_PAD_FAC*[params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
        do state=1,nstates
            call state_recs(state)%compress_exp
            call state_recs(state)%sampl_dens_correct
            call state_recs(state)%ifft
            call state_recs(state)%div(real(params%box))
            call state_recs(state)%mul(gridcorr_img)
            call state_img%copy(state_recs(state))
            if( has_lowpass_filter(state) )then
                call state_img%apply_filter(lowpass_filters(:,state))
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX PRE-IMAGE applied project-FSC low-pass filter to state=',state, &
                    &' using_source_state=',lowpass_source_state(state)
            endif
            call write_state(params,state_img,state)
            call state_img%kill
            call state_recs(state)%dealloc_rho
            call state_recs(state)%kill
        end do
        call gridcorr_img%kill
        deallocate(state_recs,scales,lowpass_filters,has_lowpass_filter,lowpass_source_state)
    end subroutine reconstruct_flex_diffmap_weighted_states

    !> Local-linear diffusion-map pre-image reconstruction.  For each medoid
    !! target z_c, a per-Fourier-voxel weighted least-squares fit of the voxel
    !! value on the centered latent coordinate u_i = z_i - z_c is solved with
    !! design row g_i = [1, u_i(1..d)].  The intercept component (index 1) is
    !! the pre-image; it removes the O(h^2) curvature bias of the local-constant
    !! (Nadaraya-Watson) estimator without imposing a single global linear
    !! subspace.  Each state is solved sequentially, reusing the coupled
    !! accumulation/solve machinery: data_scales = w_i*g_i (RHS b) and
    !! density_scales = w_i*(g_i g_i^T) (normal matrix A), with CTF^2 folded in
    !! by the accumulator and a relative ridge applied by solve_coupled_basis_exp.
    subroutine reconstruct_flex_diffmap_local_linear_states( params, build, pinds, raw_coords, medoids, &
        &state_weights, ndim, nstates, fsc_projfile )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:), medoids(:), ndim, nstates
        real,              intent(in)    :: raw_coords(:,:), state_weights(:,:)
        type(string), optional, intent(in) :: fsc_projfile
        type(reconstructor), allocatable :: basis_recs(:)
        type(fplane_type), allocatable :: model_fpls(:)
        type(ori), allocatable :: orientations(:)
        real(dp), allocatable :: data_scales(:,:), density_scales(:,:,:), zc(:), gvec(:)
        real, allocatable :: rho_cross(:,:,:,:)
        real, allocatable :: lowpass_filters(:,:)
        logical, allocatable :: has_lowpass_filter(:), valid(:)
        integer, allocatable :: lowpass_source_state(:)
        type(image) :: gridcorr_img, state_img
        integer :: nmodes, d, ncomp, exp_shape(3)
        integer :: batchlims(2), batchsz, ibatch, i, iptcl, idx, state, q, r
        real(dp) :: wi
        if( size(pinds)<1 .or. nstates<1 ) THROW_HARD('invalid flex local-linear reconstruction dimensions')
        if( any(shape(state_weights)/=[size(pinds),nstates]) ) THROW_HARD('flex local-linear weight table mismatch')
        if( size(raw_coords,1)/=size(pinds) ) THROW_HARD('flex local-linear coordinate/particle mismatch')
        if( size(medoids)<nstates ) THROW_HARD('flex local-linear medoid table smaller than nstates')
        nmodes = size(raw_coords,2)
        d      = max(1, min(ndim, nmodes))
        ncomp  = d + 1
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX PRE-IMAGE local-linear design dimension d=',d,' ncomp=',ncomp
        allocate(zc(d), gvec(ncomp))
        allocate(orientations(MAXIMGBATCHSZ), valid(MAXIMGBATCHSZ))
        allocate(data_scales(ncomp,MAXIMGBATCHSZ), density_scales(ncomp,ncomp,MAXIMGBATCHSZ), source=0.d0)
        call prepare_project_fsc_lowpass_filters(params,build,nstates,lowpass_filters,has_lowpass_filter, &
            &lowpass_source_state,fsc_projfile)
        call init_rec(params,build,MAXIMGBATCHSZ,model_fpls)
        call prepimgbatch(params,build,MAXIMGBATCHSZ)
        gridcorr_img=prep3D_inv_instrfun4mul([params%box_crop,params%box_crop,params%box_crop], &
            &OSMPL_PAD_FAC*[params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
        do state=1,nstates
            if( medoids(state)<1 .or. medoids(state)>size(pinds) ) &
                &THROW_HARD('flex local-linear medoid outside particle table')
            zc = real(raw_coords(medoids(state),1:d),dp)
            ! Fresh coupled reconstructors per state: component 1 is finalized to
            ! real space below, so reconstructors are not reused across states.
            allocate(basis_recs(ncomp))
            do q=1,ncomp
                call init_basis_reconstructor(params,build,basis_recs(q))
            end do
            if( .not.allocated(rho_cross) )then
                exp_shape = shape(basis_recs(1)%cmat_exp)
                allocate(rho_cross((ncomp*(ncomp+1))/2,exp_shape(1),exp_shape(2),exp_shape(3)),source=0.)
            else
                rho_cross = 0.
            endif
            do ibatch=1,size(pinds),MAXIMGBATCHSZ
                batchlims=[ibatch,min(size(pinds),ibatch+MAXIMGBATCHSZ-1)]
                batchsz=batchlims(2)-batchlims(1)+1
                if( params%l_ptcl_src_den )then
                    call discrete_read_imgbatch_source(params,build,'den',batchsz, &
                        &pinds(batchlims(1):batchlims(2)),[1,batchsz],build%imgbatch(:batchsz))
                else
                    call discrete_read_imgbatch(params,build,size(pinds),pinds,batchlims)
                endif
                call prep_imgs4projected_model(params,build,batchsz,build%imgbatch(:batchsz), &
                    &pinds(batchlims(1):batchlims(2)),model_fpls(:batchsz))
                valid(:batchsz)=.false.
                do i=1,batchsz
                    idx   = batchlims(1)+i-1
                    iptcl = pinds(idx)
                    call orientations(i)%kill
                    call build%spproj_field%get_ori(iptcl,orientations(i))
                    if( orientations(i)%isstatezero() ) cycle
                    if( .not.allocated(model_fpls(i)%transfer_plane) ) &
                        &THROW_HARD('flex local-linear model lacks a forward transfer plane')
                    wi       = real(state_weights(idx,state),dp)
                    gvec(1)  = 1.d0
                    do q=1,d
                        gvec(q+1) = real(raw_coords(idx,q),dp) - zc(q)
                    end do
                    do q=1,ncomp
                        data_scales(q,i) = wi*gvec(q)
                        do r=1,ncomp
                            density_scales(q,r,i) = wi*gvec(q)*gvec(r)
                        end do
                    end do
                    valid(i)=.true.
                end do
                call insert_planes_oversamp_coupled_batch_scaled(basis_recs,rho_cross,build%pgrpsyms, &
                    &orientations(:batchsz),model_fpls(:batchsz),data_scales(:,:batchsz), &
                    &density_scales(:,:,:batchsz),valid(:batchsz),batchsz)
            end do
            ! Per-voxel weighted least-squares solve; component 1 is the intercept.
            call solve_coupled_basis_exp(basis_recs,rho_cross,ncomp)
            call basis_recs(1)%compress_exp
            ! solve_coupled_basis_exp already density-corrects, so no sampl_dens_correct here.
            call basis_recs(1)%ifft
            call basis_recs(1)%div(real(params%box))
            call basis_recs(1)%mul(gridcorr_img)
            call state_img%copy(basis_recs(1))
            if( has_lowpass_filter(state) )then
                call state_img%apply_filter(lowpass_filters(:,state))
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX PRE-IMAGE applied project-FSC low-pass filter to state=',state, &
                    &' using_source_state=',lowpass_source_state(state)
            endif
            call write_state(params,state_img,state)
            call state_img%kill
            do q=1,ncomp
                call basis_recs(q)%dealloc_rho
                call basis_recs(q)%kill
            end do
            deallocate(basis_recs)
        end do
        do i=1,size(orientations)
            call orientations(i)%kill
        end do
        call gridcorr_img%kill
        call cleanup_rec_buffers(build,model_fpls)
        if( allocated(rho_cross) ) deallocate(rho_cross)
        deallocate(zc,gvec,orientations,valid,data_scales,density_scales)
        deallocate(lowpass_filters,has_lowpass_filter,lowpass_source_state)
    end subroutine reconstruct_flex_diffmap_local_linear_states

    subroutine prepare_project_fsc_lowpass_filters( params, build, nstates, lowpass_filters, has_filter, source_state, &
        &fsc_projfile )
        class(parameters), intent(in) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in) :: nstates
        real, allocatable, intent(out) :: lowpass_filters(:,:)
        logical, allocatable, intent(out) :: has_filter(:)
        integer, allocatable, intent(out) :: source_state(:)
        type(string), optional, intent(in) :: fsc_projfile
        type(sp_project) :: spproj
        type(string) :: fsc_fname, imgkind_here, proj_for_fsc
        real, allocatable :: fsc(:)
        integer :: filtsz, state, fsc_box, i, state1_fsc_count
        logical :: out_loaded
        filtsz=fdim(params%box_crop)-1
        allocate(lowpass_filters(filtsz,nstates),has_filter(nstates),source_state(nstates))
        lowpass_filters=0.
        has_filter=.false.
        source_state=0
        if( filtsz<1 ) return
        proj_for_fsc=params%projfile
        if( present(fsc_projfile) )then
            if( len_trim(fsc_projfile%to_char())>0 ) proj_for_fsc=fsc_projfile
        endif
        if( .not.file_exists(proj_for_fsc) )then
            write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE project-FSC low-pass filtering skipped: projfile not found'
            call proj_for_fsc%kill
            return
        endif
        call spproj%read_segment('out',proj_for_fsc)
        out_loaded=spproj%os_out%get_noris()>0
        if( .not.out_loaded )then
            write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE project-FSC low-pass filtering skipped: empty out segment'
            call proj_for_fsc%kill
            call spproj%kill
            return
        endif
        state1_fsc_count=0
        do i=1,spproj%os_out%get_noris()
            if( .not.spproj%os_out%isthere(i,'imgkind') ) cycle
            call spproj%os_out%getter(i,'imgkind',imgkind_here)
            if( imgkind_here%to_char()/='fsc' ) cycle
            if( spproj%os_out%get_state(i)==1 ) state1_fsc_count=state1_fsc_count+1
        end do
        if( state1_fsc_count/=1 )then
            write(logfhandle,'(A,I0)') '>>> FLEX PRE-IMAGE project-FSC low-pass filtering skipped: expected exactly one state=1 FSC in out segment, found=', &
                &state1_fsc_count
            call imgkind_here%kill
            call proj_for_fsc%kill
            call spproj%kill
            return
        endif
        call spproj%get_fsc(1,fsc_fname,fsc_box)
        if( .not.file_exists(fsc_fname) )then
            write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE project-FSC low-pass filtering skipped: state=1 FSC file missing'
            call imgkind_here%kill
            call proj_for_fsc%kill
            call spproj%kill
            return
        endif
        fsc=file2rarr(fsc_fname)
        if( size(fsc)/=filtsz )then
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX PRE-IMAGE project-FSC low-pass filtering skipped: state=1 FSC size mismatch; fsc_nyq=', &
                &size(fsc),' model_nyq=',filtsz
            deallocate(fsc)
            call fsc_fname%kill
            call imgkind_here%kill
            call proj_for_fsc%kill
            call spproj%kill
            return
        endif
        do state=1,nstates
            call fsc2optlp_sub(filtsz,fsc,lowpass_filters(:,state),merged=.false.)
            has_filter(state)=any(lowpass_filters(:,state)>0.)
            source_state(state)=1
        end do
        deallocate(fsc)
        call fsc_fname%kill
        call imgkind_here%kill
        call proj_for_fsc%kill
        call spproj%kill
    end subroutine prepare_project_fsc_lowpass_filters

    subroutine write_flex_diffmap_rec_parts( params, build, pinds, z, ncomp, part )
        class(parameters), intent(inout) :: params
        class(builder), intent(inout) :: build
        integer, intent(in) :: pinds(:), ncomp, part
        real, intent(in) :: z(:,:)
        type(parameters) :: model_params
        type(reconstructor) :: mean_rec
        type(fplane_type), allocatable :: fpls(:)
        real(dp), allocatable :: z_dp(:,:)
        type(string) :: fname
        if( size(pinds)<1 .or. part<1 .or. ncomp<1 ) THROW_HARD('invalid flex residual-model worker assignment')
        if( size(z,1)/=size(pinds) .or. size(z,2)/=ncomp ) THROW_HARD('flex residual-model worker coordinate mismatch')
        call prepare_unfiltered_model_params(params,model_params)
        allocate(z_dp(size(z,1),ncomp),source=real(z,dp))
        call init_mean_reconstructor(model_params,build,mean_rec)
        fname=flex_part_fname(part,params%numlen)
        call write_mstep_stats_part_file(model_params,build,mean_rec,z_dp,pinds,size(pinds),ncomp,[1,size(pinds)], &
            &fname,fpls,log_label='FLEX_DIFFMAP_WORKER')
        call cleanup_planes(fpls)
        call mean_rec%dealloc_rho
        call mean_rec%kill
        call fname%kill
        deallocate(z_dp)
    end subroutine write_flex_diffmap_rec_parts

    subroutine reduce_flex_diffmap_rec_parts( params, build, nparts, ncomp, target_coeffs, nstates )
        class(parameters), intent(inout) :: params
        class(builder), intent(inout) :: build
        integer, intent(in) :: nparts, ncomp, nstates
        real, intent(in) :: target_coeffs(:,:)
        type(parameters) :: model_params
        type(reconstructor), allocatable :: basis_recs(:)
        type(string), allocatable :: part_fnames(:)
        integer :: ipart, q
        if( nparts<1 .or. ncomp<1 .or. nstates<1 ) THROW_HARD('invalid distributed flex residual-model reduction')
        if( any(shape(target_coeffs)/=[nstates,ncomp]) ) THROW_HARD('flex target coefficient shape mismatch')
        call prepare_unfiltered_model_params(params,model_params)
        allocate(basis_recs(ncomp),part_fnames(nparts))
        do q=1,ncomp
            call init_basis_reconstructor(model_params,build,basis_recs(q))
        end do
        do ipart=1,nparts
            part_fnames(ipart)=flex_part_fname(ipart,params%numlen)
        end do
        call update_basis_from_mstep_stats_part_files(model_params,basis_recs,ncomp,part_fnames,nparts, &
            &log_label='FLEX_DIFFMAP')
        call write_preimage_states(params,basis_recs,target_coeffs,nstates,ncomp)
        do q=1,ncomp
            call basis_recs(q)%dealloc_rho
            call basis_recs(q)%kill
        end do
        do ipart=1,nparts
            call part_fnames(ipart)%kill
        end do
        deallocate(basis_recs,part_fnames)
    end subroutine reduce_flex_diffmap_rec_parts

    subroutine cleanup_flex_diffmap_rec_parts( nparts, numlen )
        integer, intent(in) :: nparts, numlen
        type(string) :: fname
        integer :: ipart
        do ipart=1,nparts
            fname=flex_part_fname(ipart,numlen)
            call del_file(fname)
            call fname%kill
        end do
    end subroutine cleanup_flex_diffmap_rec_parts

    subroutine validate_model_tables( pinds, z, target_coeffs, nstates, ncomp )
        integer, intent(in) :: pinds(:), nstates
        real, intent(in) :: z(:,:), target_coeffs(:,:)
        integer, intent(out) :: ncomp
        ncomp=size(z,2)
        if( size(pinds)<1 .or. ncomp<1 .or. nstates<1 ) THROW_HARD('invalid flex residual-model dimensions')
        if( size(z,1)/=size(pinds) ) THROW_HARD('flex residual-model particle coordinate mismatch')
        if( any(shape(target_coeffs)/=[nstates,ncomp]) ) THROW_HARD('flex residual-model target mismatch')
    end subroutine validate_model_tables

    subroutine prepare_unfiltered_model_params( params, model_params )
        class(parameters), intent(in) :: params
        type(parameters), intent(out) :: model_params
        model_params=params
        model_params%lp=0.
    end subroutine prepare_unfiltered_model_params

    subroutine init_mean_reconstructor( params, build, mean_rec )
        class(parameters), intent(inout) :: params
        class(builder), intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        call mean_rec%read_and_crop(params%vols(1),params%smpd,params%box_crop,params%smpd_crop)
        call mean_rec%fft
        call mean_rec%alloc_rho(params,build%spproj,expand=.true.)
        call mean_rec%expand_exp
    end subroutine init_mean_reconstructor

    subroutine init_basis_reconstructor( params, build, basis_rec )
        class(parameters), intent(inout) :: params
        class(builder), intent(inout) :: build
        type(reconstructor), intent(inout) :: basis_rec
        call basis_rec%new([params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
        call basis_rec%alloc_rho(params,build%spproj,expand=.true.)
        call basis_rec%reset
        call basis_rec%reset_exp
    end subroutine init_basis_reconstructor

    subroutine write_preimage_states( params, basis_recs, target_coeffs, nstates, ncomp )
        class(parameters), intent(in) :: params
        integer, intent(in) :: nstates, ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real, intent(in) :: target_coeffs(nstates,ncomp)
        type(image) :: mean_img, state_img
        integer :: state, q
        call mean_img%read_and_crop(params%vols(1),params%smpd,params%box_crop,params%smpd_crop)
        call mean_img%fft
        do state=1,nstates
            call state_img%copy(mean_img)
            do q=1,ncomp
                if( abs(target_coeffs(state,q))>real(DTINY) ) &
                    &call state_img%add(basis_recs(q),target_coeffs(state,q))
            end do
            call state_img%ifft
            call write_state(params,state_img,state)
            call state_img%kill
        end do
        call mean_img%kill
    end subroutine write_preimage_states

    subroutine write_state( params, img, state )
        class(parameters), intent(in) :: params
        type(image), intent(inout) :: img
        integer, intent(in) :: state
        type(string) :: fname, prefix, ext
        character(len=:), allocatable :: stem
        character(len=3) :: tag
        if( state==1 )then
            fname=params%outvol
        else
            ext=fname2ext(params%outvol)
            prefix=get_fbody(params%outvol,ext)
            stem=prefix%to_char()
            if( len_trim(stem)>4 )then
                if( stem(len_trim(stem)-3:len_trim(stem))=='_001' ) stem=stem(:len_trim(stem)-4)
            endif
            prefix=string(stem)
            write(tag,'(I3.3)') state
            fname=prefix//'_'//tag//MRC_EXT
        endif
        call img%write(fname,del_if_exists=.true.)
        write(logfhandle,'(A,I0,A,A)') '>>> FLEX DIFFMAP NYSTROM PRE-IMAGE ',state,': ',fname%to_char()
        call fname%kill
        call prefix%kill
        call ext%kill
    end subroutine write_state

    subroutine cleanup_reconstructors( mean_rec, basis_recs )
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), allocatable, intent(inout) :: basis_recs(:)
        integer :: q
        call mean_rec%dealloc_rho
        call mean_rec%kill
        if( allocated(basis_recs) )then
            do q=1,size(basis_recs)
                call basis_recs(q)%dealloc_rho
                call basis_recs(q)%kill
            end do
            deallocate(basis_recs)
        endif
    end subroutine cleanup_reconstructors

    function flex_part_fname( part, numlen ) result(fname)
        integer, intent(in) :: part, numlen
        type(string) :: fname
        fname=string('flex_diffmap_mstep_part')//int2str_pad(part,numlen)//'.bin'
    end function flex_part_fname

end module simple_flex_diffmap_rec3D
