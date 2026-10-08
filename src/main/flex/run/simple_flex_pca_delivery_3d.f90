!@descr: flex_pca 3-D delivery of latent readouts, figures, covariance tables, eigenvolumes and manifest
!!
!! Everything the run publishes besides the state maps themselves: the UMAP readouts of the
!! embedding (kept so the state stage can colour the same figure by delivered state), the
!! coordinate/weight/target/prior tables, the eigenvalue table and the manifest.
module simple_flex_pca_delivery_3d
use simple_core_module_api,  only: del_file, dp, logfhandle, simple_exception
use simple_flex_pca_records, only: flex_latent_readout, flex_selection, flex_fit_model, flex_latent, flex_state_set
use simple_builder,          only: builder
use simple_parameters,       only: parameters
use simple_umap,             only: umap_embed, umap_subsample
use simple_flex_pca_plot,    only: flex_plot_latent_jpg
implicit none
private
#include "simple_local_flags.inc"

public :: deliver_latent_readouts, write_states_umap_figure
public :: write_covariance_tables, write_covariance_eigenvolumes, write_covariance_manifest

!> UMAP coordinates of the last delivery readout (and their particle indices), kept so the state
!! stage can draw the same embedding coloured by the delivered states (flex_pca_umap_states.jpg)
contains

    !> Delivery readout on the all-N latents: UMAP plot coordinates of a bounded subsample
    !! (flex_pca_umap.txt), opt-in through umap=yes. Dead axes (near-zero latent variance) are
    !! excluded, since whitening explodes them into pure noise.
    subroutine deliver_latent_readouts( readout, sel, model, latent, l_umap, tag, labels )
        type(flex_selection),       intent(in)    :: sel
        type(flex_fit_model),       intent(in)    :: model
        type(flex_latent),          intent(in)    :: latent
        type(flex_latent_readout),  intent(inout) :: readout
        logical,                    intent(in)    :: l_umap
        character(len=*), optional, intent(in)    :: tag   !< suffix on the output names (e.g. '_deconv')
        integer,          optional, intent(in)    :: labels(:) !< per-row labels (e.g. the deconvolution's mixture component) for the figure
        real,    allocatable :: pz1(:), pz2(:)
        integer, allocatable :: plab(:)
        character(len=:), allocatable :: sfx
        integer,  allocatable :: hcols(:), uinds(:)
        real,     allocatable :: Xu(:,:), Yu(:,:)
        real(dp) :: colvar(model%ncomp), vmax, mu_c
        integer  :: vumap, nhc, q, i, u, nsub
        sfx = ''
        if( present(tag) ) sfx = tag
        vumap = merge(1, 0, l_umap)
        if( vumap < 1 ) return
        ! healthy columns: latent variance above 1e-4 of the max
        vmax = 0.d0
        do q = 1, model%ncomp
            mu_c      = sum(latent%z(1:sel%nptcls,q))/real(sel%nptcls,dp)
            colvar(q) = sum((latent%z(1:sel%nptcls,q) - mu_c)**2)/real(max(sel%nptcls-1,1),dp)
            vmax      = max(vmax, colvar(q))
        end do
        allocate(hcols(model%ncomp))
        nhc = 0
        do q = 1, model%ncomp
            if( colvar(q) > 1.d-4*vmax )then
                nhc = nhc + 1
                hcols(nhc) = q
            endif
        end do
        if( nhc < model%ncomp ) write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA delivery readouts: ', &
            &model%ncomp - nhc, ' dead/retired axis(es) excluded (', nhc, ' healthy)'
        if( vumap > 0 .and. nhc >= 3 .and. sel%nptcls > 100 )then
            nsub = min(1000000, sel%nptcls)   ! kd-tree kNN + 200 SGD epochs: every particle, not a subsample
            call umap_subsample(sel%nptcls, nsub, 1234, uinds)
            allocate(Xu(nsub,nhc))
            do q = 1, nhc
                do i = 1, nsub
                    Xu(i,q) = real(latent%z(uinds(i),hcols(q)))
                end do
            end do
            call umap_embed(Xu, 1234, Yu)
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA UMAP delivery plot: ', nsub, &
                &' latents -> flex_pca_umap.txt'
            call flush(logfhandle)
            call del_file('flex_pca_umap'//sfx//'.txt')
            open(newunit=u, file='flex_pca_umap'//sfx//'.txt', status='replace', action='write')
            write(u,'(A)') '# particle  umap1  umap2'
            do i = 1, nsub
                write(u,'(I10,2(1X,ES14.6))') sel%pinds(uinds(i)), Yu(i,1), Yu(i,2)
            end do
            close(u)
            ! the figure: log-density UMAP core, UMAP by label, z1 vs z2 by label (flex_pca_umap<sfx>.jpg);
            ! the coordinates are kept for the state stage's own colouring
            if( allocated(readout%umap_xy) )   deallocate(readout%umap_xy)
            if( allocated(readout%umap_pind) ) deallocate(readout%umap_pind)
            allocate(readout%umap_xy(2,nsub), readout%umap_pind(nsub), pz1(nsub), pz2(nsub), plab(nsub))
            do i = 1, nsub
                readout%umap_xy(1,i) = Yu(i,1); readout%umap_xy(2,i) = Yu(i,2)
                readout%umap_pind(i) = sel%pinds(uinds(i))
                pz1(i) = real(latent%z(uinds(i),hcols(1))); pz2(i) = real(latent%z(uinds(i),hcols(2)))
                plab(i) = 0
                if( present(labels) ) plab(i) = max(0, labels(uinds(i)))
            end do
            call flex_plot_latent_jpg('flex_pca_umap'//sfx//'.jpg', readout%umap_xy(1,:), readout%umap_xy(2,:), &
                &pz1, pz2, nsub, labels=plab)
            write(logfhandle,'(A)') '>>> FLEX_PCA UMAP figure: flex_pca_umap'//sfx//'.jpg'
            deallocate(Xu, Yu, uinds, pz1, pz2, plab)
        endif
        deallocate(hcols)
    end subroutine deliver_latent_readouts

    !> The delivery UMAP coloured by the delivered states (argmax state weight; 0 = no weight), written
    !! as flex_pca_umap_states.jpg next to the state tables. Needs the coordinates kept by the last
    !! deliver_latent_readouts of this process; a states-only worker that never ran it writes nothing.
    subroutine write_states_umap_figure( readout, pinds, z, weights )
        type(flex_latent_readout), intent(in) :: readout
        integer,                   intent(in) :: pinds(:)
        real(dp),                  intent(in) :: z(:,:)
        real,                      intent(in) :: weights(:,:)
        integer, allocatable :: lut(:), lab(:)
        real,    allocatable :: pz1(:), pz2(:)
        integer :: i, row, nsub, pmax, s, sbest
        real    :: wbest
        if( .not. allocated(readout%umap_pind) ) return
        if( size(z,2) < 2 ) return
        nsub = size(readout%umap_pind); pmax = max(maxval(pinds), maxval(readout%umap_pind))
        allocate(lut(pmax), source=0)
        do i = 1, size(pinds); lut(pinds(i)) = i; end do
        allocate(lab(nsub), pz1(nsub), pz2(nsub))
        do i = 1, nsub
            row = lut(readout%umap_pind(i))
            if( row < 1 )then
                lab(i) = 0; pz1(i) = 0.0; pz2(i) = 0.0; cycle
            endif
            pz1(i) = real(z(row,1)); pz2(i) = real(z(row,2))
            sbest = 0; wbest = 0.0
            do s = 1, size(weights,2)
                if( weights(row,s) > wbest )then; wbest = weights(row,s); sbest = s; endif
            end do
            lab(i) = sbest
        end do
        call flex_plot_latent_jpg('flex_pca_umap_states.jpg', readout%umap_xy(1,:), readout%umap_xy(2,:), &
            &pz1, pz2, nsub, labels=lab, nlab=size(weights,2))
        write(logfhandle,'(A)') '>>> FLEX_PCA UMAP figure by delivered state: flex_pca_umap_states.jpg'
        deallocate(lut, lab, pz1, pz2)
    end subroutine write_states_umap_figure

    subroutine write_covariance_tables( readout, build, sel, model, latent, states )
        type(flex_selection),      intent(in)    :: sel
        type(flex_fit_model),      intent(in)    :: model
        type(flex_latent),         intent(in)    :: latent
        type(flex_state_set),      intent(in)    :: states
        type(flex_latent_readout), intent(in)    :: readout
        type(builder),             intent(inout) :: build
        integer :: u, i, q, state
        call del_file('flex_pca_coordinates.txt')
        open(newunit=u,file='flex_pca_coordinates.txt',status='replace',action='write')
        write(u,'(A)',advance='no') '# particle eo label residual mean_residual'
        do q=1,size(latent%z,2); write(u,'(A,I0)',advance='no') ' z',q; end do
        write(u,*)
        do i=1,size(sel%pinds)
            write(u,'(I10,1X,I1,1X,I4,2(1X,ES16.8))',advance='no') sel%pinds(i), &
                &build%spproj_field%get_eo(sel%pinds(i)),states%labels(i),latent%resid_energy(i),latent%resid_mean_energy(i)
            do q=1,size(latent%z,2); write(u,'(1X,ES16.8)',advance='no') latent%z(i,q); end do
            write(u,*)
        end do
        close(u)
        call del_file('flex_pca_contrast.txt')
        open(newunit=u,file='flex_pca_contrast.txt',status='replace',action='write')
        write(u,'(A)') '# particle contrast'
        do i=1,min(size(sel%pinds),size(latent%contrast))
            write(u,'(I10,1X,ES16.8)') sel%pinds(i), latent%contrast(i)
        end do
        close(u)
        call del_file('flex_pca_state_weights.txt')
        open(newunit=u,file='flex_pca_state_weights.txt',status='replace',action='write')
        write(u,'(A)',advance='no') '# particle'
        do state=1,size(states%weights,2); write(u,'(A,I3.3)',advance='no') ' w',state; end do
        write(u,*)
        do i=1,size(sel%pinds)
            write(u,'(I10)',advance='no') sel%pinds(i)
            do state=1,size(states%weights,2); write(u,'(1X,ES16.8)',advance='no') states%weights(i,state); end do
            write(u,*)
        end do
        close(u)
        call write_states_umap_figure(readout, sel%pinds, latent%z, states%weights)
        call del_file('flex_pca_state_targets.txt')
        open(newunit=u,file='flex_pca_state_targets.txt',status='replace',action='write')
        ! targets are full latent-space points, so every coordinate is written
        write(u,'(A)',advance='no') '# state bandwidth effective_particles'
        do q=1,size(states%targets,1); write(u,'(A,I0)',advance='no') ' t',q; end do
        write(u,*)
        do state=1,size(states%targets,2)
            write(u,'(I5,2(1X,ES16.8))',advance='no') state,states%bandwidths(state),states%neff(state)
            do q=1,size(states%targets,1); write(u,'(1X,ES16.8)',advance='no') states%targets(q,state); end do
            write(u,*)
        end do
        close(u)
        call del_file('flex_pca_map_prior.txt')
        open(newunit=u,file='flex_pca_map_prior.txt',status='replace',action='write')
        write(u,'(A)') '# component covariance_eigenvalue prior_precision'
        do q=1,size(model%eigvals)
            write(u,'(I5,2(1X,ES20.10))') q,model%eigvals(q),latent%prior_precision(q)
        end do
        close(u)
    end subroutine write_covariance_tables

    subroutine write_covariance_eigenvolumes( eigvals, ncomp )
        real(dp), intent(in) :: eigvals(ncomp)
        integer,  intent(in) :: ncomp
        character(len=:), allocatable :: fn
        integer :: q, u
        fn = 'flex_pca_eigenvalues.txt'
        ! only the eigenvalue table is written here; the eigenvolume MRCs are not
        call del_file(fn)
        open(newunit=u,file=fn,status='replace',action='write')
        write(u,'(A)') '# component eigenvalue'
        do q = 1, ncomp
            write(u,'(I6,1X,ES20.10)') q,eigvals(q)
        end do
        close(u)
    end subroutine write_covariance_eigenvolumes

    subroutine write_covariance_manifest( params, nptcls, ncomp, nstates, axis, min_neff, sigma_loaded )
        type(parameters), intent(in) :: params
        integer,          intent(in) :: nptcls, ncomp, nstates, axis, min_neff
        logical,          intent(in) :: sigma_loaded
        integer :: u
        call del_file('flex_pca_manifest.txt')
        open(newunit=u,file='flex_pca_manifest.txt',status='replace',action='write')
        write(u,'(A)') 'method=matched_kb_selected_column_covariance'
        write(u,'(A)') 'diffusion_map_dependency=no'
        write(u,'(A)') 'interpolation=matched_simple_kaiser_bessel'
        write(u,'(A,L1)') 'sigma_whitened=',sigma_loaded
        write(u,'(A,I0)') 'particles=',nptcls
        write(u,'(A,I0)') 'box_crop=',params%box_crop
        write(u,'(A,F10.4)') 'smpd_crop=',params%smpd_crop
        write(u,'(A,F10.4)') 'lowpass_angstrom=',params%lp
        write(u,'(A,I0)') 'column_separation=',params%column_separation
        write(u,'(A,I0)') 'probe_iters=',params%n_probe_iters
        ! the covariance band is capped at fdim(box_crop)-1
        write(u,'(A,L1)') 'lowpass_active=',(params%lp > 2.0*params%smpd_crop)
        write(u,'(A,I0)') 'components=',ncomp
        write(u,'(A,I0)') 'states=',nstates
        if( axis < 0 )then
            write(u,'(A)')    'state_placement=density_spread_path'
        else if( axis == 0 .and. trim(params%state_placement) == 'equal_occ' )then
            write(u,'(A)')    'state_placement=equal_occupancy_path'
        else if( axis == 0 )then
            write(u,'(A)')    'state_placement=diffusion_kcenter'
        else
            write(u,'(A)')    'state_placement=single_axis_quantiles'
        endif
        write(u,'(A,I0)') 'state_axis=',axis
        write(u,'(A,I0)') 'minimum_state_neff=',min_neff
        write(u,'(A)') 'half_maps=combined_even_odd'
        write(u,'(A)') 'validation=inspect_half_map_agreement_nuisance_correlations_and_heldout_residuals'
        ! cross-fit-FSC provenance, always written (the ridge is hard-wired on): COUPLED, because the ridge
        ! deflates the cross-fit statistic's independence as gold-standard FSC regularization does
        write(u,'(A,I0)') 'crossfsc_reg_mode=',1
        write(u,'(A,L1)') 'crossfsc_coupled=',.true.
        close(u)
    end subroutine write_covariance_manifest

end module simple_flex_pca_delivery_3d
