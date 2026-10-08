!@descr: flex_pca state inference and bandwidth-CV reconstruction orchestration
!!
!! Everything between the delivered embedding and the weight table the reconstruction consumes:
!! the measurement-error deconvolution of the latents (and its resume adoption), placement of
!! the state kernels (population-floor variant here; the kernel/GMM placement in
!! simple_flex_pca_weights), the occupancy floor, derived settings, trial-map bandwidth CV, and the
!! raw state half maps FLEX infers from, through the reconstruction service.
module simple_flex_pca_state_service
use simple_core_module_api, only: del_file, dp, dtiny, file_exists, hpsort, int2str_pad, irnd_uni, logfhandle, &
    &mrc_ext, simple_exception, string
use simple_builder,               only: builder
use simple_cmdline,               only: cmdline
use simple_image,                 only: image
use simple_flex_pca_records,      only: flex_selection, flex_fit_model, flex_latent, flex_state_set
use simple_parameters,            only: parameters
use simple_srch_sort_loc,         only: hpsort
use simple_flex_pca_deconv,       only: calibrate_noise_scale, deconvolve_latent
use simple_flex_pca_embedding_io, only: read_deconv_block, read_noise_scale, write_deconv_block, write_noise_scale
use simple_flex_pca_rounds,       only: flex_pca_rounds
use simple_flex_pca_weights,      only: build_covariance_state_weights, kernel_weights_at_bandwidth, &
    &bandwidth_cv_grid, bandwidth_cv_error, bandwidth_cv_adopt
implicit none
private
#include "simple_local_flags.inc"

public :: apply_latent_deconvolution, prune_underpopulated_states, place_states_with_population_floor
public :: cv_select_bandwidths, reconstruct_state_halves, state_half_fname, delete_state_halves
public :: auto_box_crop, auto_min_neff, auto_state_count
public :: infile_path
public :: FLEX_AUTO_K_START, FLEX_AUTO_K_MIN, AUTO_NSTATES

!> a weighted PARTITION row sums to one across the states within this (as in simple_state_weight_set)
real(dp), parameter :: PARTITION_ROW_TOL = 1.0e-3_dp

!> Cap and floor of auto_state_count (tester only; preimage_auto's ceiling is AUTO_NSTATES). The cap
!! is bounded by cost: gate 2 of the state merge compares K(K-1)/2 map pairs.
integer, parameter :: FLEX_AUTO_K_START        = 32
integer, parameter :: FLEX_AUTO_K_MIN          = 8
!> Nyquist margin for a derived box_crop; columns are selected inside that band.
real,    parameter :: FLEX_AUTO_BOX_SAFETY     = 1.25
!> Occupancy share of an equal split. Paired with an SNR term because neither alone reproduces both
!! validation datasets.
real,    parameter :: FLEX_AUTO_NEFF_OCCUPANCY = 0.10
! npreimages is a PROVISION CEILING, not a target: state placement lays down that many kernels and
! the two-gate merge collapses the indistinct ones, so the recovered K is only ever <= it.
! preimage_auto=yes raises that ceiling to AUTO_NSTATES and turns the merge on, since over-provisioning
! is the only regime in which the merge can recover K at all.
!> provision cap of the population floor (min_state_frac > 0, refine3D_states flex initialization); independent of AUTO_NSTATES
integer, parameter :: AUTO_NSTATES             = 8
!> provision cap of the population floor (min_state_frac > 0, refine3D_states flex initialization); independent of AUTO_NSTATES
integer, parameter :: POP_FLOOR_MAX_NSTATES    = 32

contains

    !> The raw state half maps FLEX infers from (D11), one per weight column: the reconstruction service
    !! in this process on the covariance box (D9), no delivery filter, the shell density floor of
    !! fractional weights; written as <prefix>_stateNN{,_even,_odd}.mrc (state_half_fname)
    subroutine reconstruct_state_halves( params, build, cline, pinds, weights, prefix )
        use simple_rec3D_service, only: rec3D_service, rec3D_request, rec3D_backend_id, REC3D_WEIGHTS_TABLE, &
            &REC3D_OUTPUT_RAW, REC3D_DISPATCH_INPROC
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(in)    :: cline
        integer,           intent(in)    :: pinds(:)
        real,              intent(in)    :: weights(:,:)   !< (size(pinds), nstates)
        character(len=*),  intent(in)    :: prefix
        type(rec3D_service), allocatable :: service
        type(rec3D_request) :: request
        type(cmdline)       :: cline_rec
        type(string), allocatable :: vols_bak(:)
        character(len=len(params%automsk)) :: automsk_bak
        integer :: nstates_bak, nparts_bak, part_bak, numlen_bak, nstates
        logical :: l_nonuniform_bak, l_trail_rec_bak, l_state_defined_bak
        nstates = size(weights,2)
        if( size(weights,1) /= size(pinds) .or. nstates < 1 ) THROW_HARD('state weight table shape; reconstruct_state_halves')
        ! an all-zero column would leave a state without half maps
        if( any(maxval(weights, dim=1) <= 0.) ) THROW_HARD('a state of the weight table has no weight; reconstruct_state_halves')
        ! The service reconstructs every column as a state of a single in-process part. The PCG path reads
        ! these parameters directly (volassemble rebuilds its own from cline_rec): raw maps carry no
        ! nonuniform filter, envelope support, trailing chain or state selection, whatever the caller's
        ! command line holds.
        nstates_bak            = params%nstates
        nparts_bak             = params%nparts
        part_bak               = params%part
        numlen_bak             = params%numlen
        l_nonuniform_bak       = params%l_nonuniform
        l_trail_rec_bak        = params%l_trail_rec
        l_state_defined_bak    = params%l_state_defined
        automsk_bak            = params%automsk
        vols_bak               = params%vols(1:nstates)
        params%nstates         = nstates
        params%nparts          = 1
        params%part            = 1
        params%numlen          = 1
        params%l_nonuniform    = .false.
        params%l_trail_rec     = .false.
        params%l_state_defined = .false.
        params%automsk         = 'no'
        cline_rec = cline
        call cline_rec%set('nstates',     size(weights,2))
        call cline_rec%set('box_crop',    params%box_crop)
        call cline_rec%set('mkdir',       'no')
        call cline_rec%set('ml_reg',      'no')
        call cline_rec%set('postprocess', 'no')
        call cline_rec%set('rec_backend', trim(params%rec_states_backend))
        call cline_rec%set('rho_floor',   'yes')
        call cline_rec%set('filt_mode',   'none')
        call cline_rec%set('automsk',     'no')
        call cline_rec%delete('nparts')
        call cline_rec%delete('part')
        call cline_rec%delete('vol1')
        call cline_rec%delete('trail_rec')
        call cline_rec%delete('trail_seed')
        call cline_rec%delete('ufrac_trec')
        call cline_rec%delete('frozen_rec')
        call cline_rec%delete('state')
        request%pinds         = pinds
        request%weights_table = weights
        request%weights       = REC3D_WEIGHTS_TABLE
        request%backend       = rec3D_backend_id(trim(params%rec_states_backend))
        request%output        = REC3D_OUTPUT_RAW
        request%prefix        = prefix
        request%register      = .false.
        request%dispatch      = REC3D_DISPATCH_INPROC
        allocate(service)
        call service%new(params, build, cline_rec, REC3D_DISPATCH_INPROC, params%nthr)
        call service%execute(params, build, cline_rec, request)
        call service%kill
        deallocate(service)
        call cline_rec%kill
        params%nstates          = nstates_bak
        params%nparts           = nparts_bak
        params%part             = part_bak
        params%numlen           = numlen_bak
        params%l_nonuniform     = l_nonuniform_bak
        params%l_trail_rec      = l_trail_rec_bak
        params%l_state_defined  = l_state_defined_bak
        params%automsk          = automsk_bak
        params%vols(1:nstates)  = vols_bak
        deallocate(vols_bak)
    end subroutine reconstruct_state_halves

    !> <prefix>_stateNN<suffix>.mrc, suffix '', '_even', '_odd', '_even_unfil' or '_odd_unfil'
    function state_half_fname( prefix, state, suffix ) result( fname )
        character(len=*), intent(in) :: prefix, suffix
        integer,          intent(in) :: state
        type(string) :: fname
        fname = prefix//'_state'//int2str_pad(state,2)//suffix//MRC_EXT
    end function state_half_fname

    !> remove the maps and FSC of reconstruct_state_halves
    subroutine delete_state_halves( prefix, nstates )
        character(len=*), intent(in) :: prefix
        integer,          intent(in) :: nstates
        type(string) :: fname
        integer      :: state
        do state = 1, nstates
            fname = state_half_fname(prefix, state, '')
            call del_file(fname)
            fname = state_half_fname(prefix, state, '_even')
            call del_file(fname)
            fname = state_half_fname(prefix, state, '_odd')
            call del_file(fname)
            fname = state_half_fname(prefix, state, '_even_unfil')
            call del_file(fname)
            fname = state_half_fname(prefix, state, '_odd_unfil')
            call del_file(fname)
            fname = prefix//'_fsc_state'//int2str_pad(state,2)//'.bin'
            call del_file(fname)
        end do
        call fname%kill
    end subroutine delete_state_halves

    !> Trial half maps per bin from the service; weights owns the CV grid, score and final adoption.
    !! The statistic sees the maps under the soft spherical mask (FLEX's own mask, D11).
    subroutine cv_select_bandwidths( params, build, cline, sel, nbins, min_neff, states )
        class(parameters),    intent(inout) :: params
        class(builder),       intent(inout) :: build
        class(cmdline),       intent(in)    :: cline
        type(flex_selection), intent(in)    :: sel
        integer,              intent(in)    :: nbins, min_neff
        type(flex_state_set), intent(inout) :: states
        real,     allocatable :: wbin(:,:), tgt_ev(:,:), tgt_od(:,:)
        real,     allocatable :: rmat(:,:,:), ev_flat(:), od_flat(:)
        real(dp), allocatable :: bins(:,:), hbin(:,:), err(:,:)
        type(image)  :: ev, od
        type(string) :: prefix
        real(dp) :: h_used
        real     :: neff_used
        integer  :: state, ib, nvox, ldim(3), i
        character(len=3) :: bstr
        allocate(wbin(sel%nptcls,states%nstates))
        allocate(bins(nbins,states%nstates), hbin(nbins,states%nstates), err(nbins,states%nstates), source=0.d0)
        call bandwidth_cv_grid(states, nbins, bins)
        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA cross-validated bandwidth selection over ',nbins, &
            &' bins per state'
        call flush(logfhandle)
        ldim = [params%box_crop, params%box_crop, params%box_crop]
        nvox = product(ldim)
        allocate(tgt_ev(nvox,states%nstates), tgt_od(nvox,states%nstates), source=0.)
        allocate(ev_flat(nvox), od_flat(nvox))
        do ib = 1, nbins
            do state = 1, states%nstates
                call kernel_weights_at_bandwidth(states%kdist(:,state), sel%nptcls, sqrt(2.d0*bins(ib,state)), &
                    &min_neff, wbin(:,state), h_used, neff_used)
                hbin(ib,state) = h_used
            end do
            write(bstr,'(I3.3)') ib
            prefix = 'flex_pca_cv'//bstr
            ! the even and odd maps of a bin come from one service call
            call reconstruct_state_halves(params, build, cline, sel%pinds, wbin, prefix%to_char())
            do state = 1, states%nstates
                call ev%new(ldim, params%smpd_crop)
                call od%new(ldim, params%smpd_crop)
                call ev%read(state_half_fname(prefix%to_char(), state, '_even'))
                call od%read(state_half_fname(prefix%to_char(), state, '_odd'))
                call ev%mask3D_soft(params%msk_crop, backgr=0.)
                call od%mask3D_soft(params%msk_crop, backgr=0.)
                rmat    = ev%get_rmat()
                ev_flat = reshape(rmat, [nvox])
                rmat    = od%get_rmat()
                od_flat = reshape(rmat, [nvox])
                if( ib == 1 )then
                    tgt_ev(:,state) = ev_flat
                    tgt_od(:,state) = od_flat
                endif
                err(ib,state) = bandwidth_cv_error(ev_flat, od_flat, tgt_ev(:,state), tgt_od(:,state))
                call ev%kill; call od%kill
            end do
            call delete_state_halves(prefix%to_char(), states%nstates)
            write(logfhandle,'(A,I3,A,ES11.3,A,ES12.4,A,ES12.4)') '>>>   cv bin=',ib, &
                &' h(state1)=',hbin(ib,1),'  cross-halfset error: min=',minval(err(ib,:)),' max=',maxval(err(ib,:))
            call flush(logfhandle)
        end do
        call bandwidth_cv_adopt(states, bins, err, min_neff)
        ! the adopted bandwidths rewrite every weight column: the labels follow the pruning ruling (D10),
        ! a row is active when its weights sum above zero and is labelled by its largest weight
        do i = 1, sel%nptcls
            if( sum(states%weights(i,:)) > 0. )then
                states%labels(i) = maxloc(states%weights(i,:), dim=1)
            else
                states%labels(i) = 0
            endif
        end do
        call prefix%kill
        deallocate(wbin, bins, hbin, err, tgt_ev, tgt_od, ev_flat, od_flat)
    end subroutine cv_select_bandwidths

    !> Calibrate the per-particle noise (from the even/odd half solutions when the run has them,
    !! else from the embedding artifact the original run wrote) and replace z / precision by the
    !! posterior means / precisions under the deconvolved mixture prior.
    subroutine apply_latent_deconvolution( latent, model, sel, applied, labels, resume, adopted, srcfile )
        type(flex_latent),    intent(inout) :: latent  !< z and precision deconvolved in place; zhalf consumed
        type(flex_fit_model), intent(in)    :: model
        type(flex_selection), intent(in)    :: sel
        integer, allocatable, intent(inout) :: labels(:)   !< mixture component per particle
        logical,              intent(in)    :: resume
        logical,              intent(out)   :: adopted
        !> the embedding artifact a resume was given (infile): its deconvolution is adopted when present,
        !! so a states-only resume, also in a fresh directory, does not re-run the K ladder; else its
        !! noise scale is used
        character(len=*),     intent(in)    :: srcfile
        logical,              intent(out)   :: applied
        real(dp) :: prior(model%ncomp), a_comp(model%ncomp), noise_scale
        integer  :: q, k_deconv
        logical  :: l_resume, l_found
        applied = .false.
        adopted = .false.
        l_resume = resume
        ! ---- resume: adopt the deconvolved block of the infile cache (one file carries raw + deconvolved) ----
        if( l_resume )then
            if( len_trim(srcfile) > 0 )then
                if( file_exists(srcfile) )then
                    block
                        real(dp), allocatable :: zb(:,:), pb(:,:,:)
                        integer,  allocatable :: lb(:)
                        real(dp) :: nsb
                        call read_deconv_block(srcfile, sel%nptcls, model%ncomp, zb, pb, lb, nsb, l_found)
                        if( l_found )then
                            latent%z(1:sel%nptcls,1:model%ncomp) = zb
                            latent%precision(1:model%ncomp,1:model%ncomp,1:sel%nptcls) = pb
                            if( allocated(labels) ) deallocate(labels)
                            allocate(labels(sel%nptcls), source=lb)
                            if( any(labels < 1) ) deallocate(labels)
                            write(logfhandle,'(A)') '>>> FLEX_PCA resumed embedding: deconvolved block adopted from '//&
                                &trim(srcfile)//' (no re-deconvolution)'
                            call flush(logfhandle)
                            applied = .true.
                            adopted = .true.
                            deallocate(zb, pb, lb)
                            return
                        endif
                    end block
                endif
            endif
        endif
        do q = 1, model%ncomp
            prior(q) = 1.d0 / max(model%eigvals(q), DTINY)
        end do
        if( allocated(latent%zhalf) )then
            call calibrate_noise_scale(latent%zhalf, latent%precision, prior, sel%nptcls, model%ncomp, noise_scale, a_comp)
            ! the scale joins the run's embedding artifact at once: a resume after a crash in the K ladder
            ! deconvolves at the scale of the fit
            if( file_exists('flex_pca_embedding.bin') ) &
                &call write_noise_scale('flex_pca_embedding.bin', sel%nptcls, model%ncomp, noise_scale)
            deallocate(latent%zhalf)
        else
            ! the calibrated noise scale of the fit, from the resume's embedding artifact
            call read_noise_scale(srcfile, sel%nptcls, model%ncomp, noise_scale, l_found)
            write(logfhandle,'(A,A,A,F8.3)') '>>> FLEX_PCA resumed embedding: noise scale from ', trim(srcfile), &
                &' (1.0 when absent) =', noise_scale
        endif
        call deconvolve_latent(latent%z, latent%precision, prior, sel%nptcls, model%ncomp, noise_scale, 16, k_deconv, &
            &prior_fname='flex_pca_deconv_prior.txt', labels_fname='flex_pca_deconv_labels.txt', pinds=sel%pinds, &
            &labels_out=labels)
        applied = .true.
        ! the deconvolved coordinates join the run's embedding artifact as its trailing block (one file)
        if( file_exists('flex_pca_embedding.bin') )then
            if( allocated(labels) )then
                call write_deconv_block('flex_pca_embedding.bin', sel%nptcls, model%ncomp, latent%z, latent%precision, labels, noise_scale)
            else
                call write_deconv_block('flex_pca_embedding.bin', sel%nptcls, model%ncomp, latent%z, latent%precision, noise_scale=noise_scale)
            endif
        else
            write(logfhandle,'(A)') '>>> FLEX_PCA deconvolved coordinates not cached: no flex_pca_embedding.bin in this directory'
        endif
    end subroutine apply_latent_deconvolution

    !> Drop every state whose effective sample size is below min_neff and compact the per-state arrays;
    !! the weights then decide (pruning ruling): a row's label is its argmax over the kept states, PARTITION
    !! rows are renormalized over them, a row without kept weight is unassigned. At least two states survive.
    subroutine prune_underpopulated_states( min_neff, states )
        type(flex_state_set), intent(inout) :: states
        integer :: nptcls
        integer,              intent(in)    :: min_neff
        logical,  allocatable :: keep(:)
        integer,  allocatable :: map(:), occ(:), ord(:)
        real,     allocatable :: w2(:,:), t2(:,:), b2(:), n2(:)
        real(dp), allocatable :: d2(:,:), f2(:)
        real,     allocatable :: key(:)
        real(dp) :: rowsum
        integer  :: s, i, nkeep, ndrop, nlost, nrelabelled, snew
        logical  :: l_partition
        nptcls = size(states%weights,1)
        if( states%nstates < 2 ) return
        allocate(keep(states%nstates), source=.true.)
        allocate(map(states%nstates), occ(states%nstates), source=0)
        do i = 1, nptcls
            if( states%labels(i) >= 1 .and. states%labels(i) <= states%nstates ) occ(states%labels(i)) = occ(states%labels(i)) + 1
        end do
        do s = 1, states%nstates
            keep(s) = states%neff(s) >= real(min_neff)
        end do
        nkeep = count(keep)
        if( nkeep >= states%nstates ) return
        if( nkeep < 2 )then
            ! nothing clears the floor: keep the two best-supported seats rather than delivering none
            allocate(key(states%nstates), ord(states%nstates))
            key = states%neff
            do s = 1, states%nstates
                ord(s) = s
            end do
            call hpsort(key, ord)
            keep = .false.
            keep(ord(states%nstates))   = .true.
            keep(ord(states%nstates-1)) = .true.
            nkeep = 2
            deallocate(key, ord)
        endif
        ! PARTITION when every weighted row sums to one across the states (the weight set's rule)
        l_partition = .true.
        do i = 1, nptcls
            rowsum = sum(real(states%weights(i,:),dp))
            if( rowsum <= 0._dp ) cycle
            if( abs(rowsum - 1._dp) > PARTITION_ROW_TOL )then
                l_partition = .false.
                exit
            endif
        end do
        ndrop = states%nstates - nkeep
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA OCCUPANCY FLOOR: ', ndrop, ' of ', states%nstates, &
            &' states have fewer than ', min_neff, ' effective particles and are dropped before reconstruction'
        do s = 1, states%nstates
            if( keep(s) ) cycle
            write(logfhandle,'(A,I3,A,F10.1,A,I0)') '>>>   dropped state=', s, '  neff=', states%neff(s), &
                &'  particles=', occ(s)
        end do
        snew = 0
        do s = 1, states%nstates
            if( .not. keep(s) ) cycle
            snew   = snew + 1
            map(s) = snew
        end do
        allocate(w2(nptcls,nkeep), t2(size(states%targets,1),nkeep), b2(nkeep), n2(nkeep))
        snew = 0
        do s = 1, states%nstates
            if( .not. keep(s) ) cycle
            snew = snew + 1
            w2(:,snew) = states%weights(:,s)
            t2(:,snew) = states%targets(:,s)
            b2(snew)   = states%bandwidths(s)
            n2(snew)   = states%neff(s)
        end do
        call move_alloc(w2, states%weights)
        call move_alloc(t2, states%targets)
        call move_alloc(b2, states%bandwidths)
        call move_alloc(n2, states%neff)
        if( allocated(states%kdist) )then
            allocate(d2(nptcls,nkeep))
            snew = 0
            do s = 1, states%nstates
                if( .not. keep(s) ) cycle
                snew = snew + 1
                d2(:,snew) = states%kdist(:,s)
            end do
            call move_alloc(d2, states%kdist)
        endif
        if( allocated(states%kfloor) )then
            allocate(f2(nkeep))
            snew = 0
            do s = 1, states%nstates
                if( .not. keep(s) ) cycle
                snew = snew + 1
                f2(snew) = states%kfloor(s)
            end do
            call move_alloc(f2, states%kfloor)
        endif
        ! the weights are authoritative (pruning ruling): a row with kept weight is relabelled by its
        ! argmax over the kept states (PARTITION rows renormalized over them), a row without is deselected
        nlost       = 0
        nrelabelled = 0
        do i = 1, nptcls
            rowsum = sum(real(states%weights(i,:),dp))
            if( rowsum <= 0._dp )then
                if( states%labels(i) > 0 ) nlost = nlost + 1
                states%labels(i) = 0
                cycle
            endif
            if( l_partition ) states%weights(i,:) = real(real(states%weights(i,:),dp) / rowsum)
            if( states%labels(i) >= 1 .and. states%labels(i) <= states%nstates )then
                if( .not. keep(states%labels(i)) ) nrelabelled = nrelabelled + 1
            endif
            states%labels(i) = maxloc(states%weights(i,:), dim=1)
        end do
        write(logfhandle,'(A,I0,A,F6.2,A,I0,A,I0,A)') '>>> FLEX_PCA OCCUPANCY FLOOR: ', nlost, ' particles (', &
            &100.0*real(nlost)/real(max(nptcls,1)), '%) have no weight in a kept state and feed no map; ', &
            &nrelabelled, ' relabelled by their largest kept weight; ', nkeep, ' states remain'
        call flush(logfhandle)
        states%nstates = nkeep
        deallocate(keep, map, occ)
    end subroutine prune_underpopulated_states

    !> Smallest even crop that still resolves lp with margin: smpd_crop = smpd*box/box_crop and the
    !! crop's Nyquist is 2*smpd_crop, so lp needs box_crop > 2*box*smpd/lp.
    pure integer function auto_box_crop( box, smpd, lp ) result( bc )
        integer, intent(in) :: box
        real,    intent(in) :: smpd, lp
        if( lp <= 0. .or. smpd <= 0. .or. box <= 0 )then
            bc = box
            return
        endif
        bc = 2*nint(0.5*FLEX_AUTO_BOX_SAFETY*2.0*real(box)*smpd/lp)   ! nearest even
        bc = max(32, min(box, bc))
    end function auto_box_crop

    !> Minimum effective particles per state: the larger of the SNR requirement (~1/s particles for
    !! unit conformational SNR) and an occupancy floor. IgG is limited by the first, Ribosembly the
    !! second, so neither term alone suffices.
    pure integer function auto_min_neff( nptcls, nstates, snr_best ) result( mn )
        integer,  intent(in) :: nptcls, nstates
        real(dp), intent(in) :: snr_best          !< best per-component conformational SNR, 0 if unknown
        integer :: n_snr, n_occ
        n_snr = 20
        if( snr_best > 0.d0 ) n_snr = max(20, nint(1.d0/snr_best))
        n_occ = 20
        if( nstates > 0 ) n_occ = nint(FLEX_AUTO_NEFF_OCCUPANCY*real(nptcls)/real(nstates))
        mn = max(20, min(nptcls, max(n_snr, n_occ)))
    end function auto_min_neff

    !> Over-provisioned state count: FLEX_AUTO_K_START, capped by nptcls/(4*min_neff) and floored at FLEX_AUTO_K_MIN.
    pure integer function auto_state_count( nptcls, min_neff ) result( k )
        integer, intent(in) :: nptcls, min_neff
        k = FLEX_AUTO_K_START
        if( min_neff > 0 ) k = min(k, nptcls/(4*min_neff))
        k = max(FLEX_AUTO_K_MIN, k)
    end function auto_state_count

    subroutine place_states_with_population_floor( latent, model, nkern, axis, min_neff, min_state_frac, states, equal_occ )
        use simple_rnd, only: irnd_uni
        type(flex_latent),    intent(in)    :: latent
        type(flex_fit_model), intent(in)    :: model
        type(flex_state_set), intent(inout) :: states  !< nstates in: the requested count; the delivered set out
        logical, optional,    intent(in)    :: equal_occ
        type(flex_latent)    :: lat_r
        type(flex_state_set) :: st_r
        integer :: nptcls, ncomp, nstates_req
        integer,              intent(in)    :: nkern, axis, min_neff
        real,                 intent(in)    :: min_state_frac
        integer, parameter :: ROUND_CAP = 8
        real(dp), allocatable :: sdv(:)
        integer,  allocatable :: idx(:), occ(:), order(:), kept(:), deliver(:)
        logical,  allocatable :: retained(:), qualifies(:)
        integer  :: nmin, K, round, nret, nk, i, q, s, t, nqual, nrand, nsurplus, kbest, itmp
        real(dp) :: d2, dbest, zbar
        logical  :: l_success
        nptcls      = size(latent%z,1)
        ncomp       = size(latent%z,2)
        nstates_req = states%nstates
        if( allocated(states%weights) )    deallocate(states%weights)
        if( allocated(states%targets) )    deallocate(states%targets)
        if( allocated(states%bandwidths) ) deallocate(states%bandwidths)
        if( allocated(states%neff) )       deallocate(states%neff)
        if( allocated(states%labels) )     deallocate(states%labels)
        if( allocated(states%kdist) )      deallocate(states%kdist)
        if( allocated(states%kfloor) )     deallocate(states%kfloor)
        nk   = max(1, min(ncomp, nkern))
        nmin = max(1, nint(min_state_frac * real(nptcls)))
        if( nstates_req * nmin > nptcls )then
            THROW_HARD('min_state_frac is too large for the requested state count: the floors exceed the particle count')
        endif
        allocate(retained(nptcls), source=.true.)
        K         = nstates_req
        nret      = nptcls
        nqual     = 0
        l_success = .false.
        do round = 1, ROUND_CAP
            nret = count(retained)
            if( allocated(idx) ) deallocate(idx)
            allocate(idx(nret))
            t = 0
            do i = 1, nptcls
                if( retained(i) )then
                    t = t + 1
                    idx(t) = i
                endif
            end do
            allocate(lat_r%z(nret,ncomp), lat_r%precision(ncomp,ncomp,nret))
            do i = 1, nret
                lat_r%z(i,:)   = latent%z(idx(i),:)
                lat_r%precision(:,:,i) = latent%precision(:,:,idx(i))
            end do
            st_r%nstates = K
            if( allocated(latent%comp_rho) ) lat_r%comp_rho = latent%comp_rho
            call build_covariance_state_weights(lat_r, nkern, axis, max(20, min(min_neff, nret/2)), st_r, equal_occ=equal_occ)
            deallocate(lat_r%z, lat_r%precision)
            if( allocated(occ) ) deallocate(occ, qualifies)
            allocate(occ(K), source=0)
            allocate(qualifies(K), source=.false.)
            do i = 1, nret
                if( st_r%labels(i) >= 1 ) occ(st_r%labels(i)) = occ(st_r%labels(i)) + 1
            end do
            qualifies = occ >= nmin
            nqual     = count(qualifies)
            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA POPULATION FLOOR round=', round, &
                &' provisioned=', K, ' retained=', nret, ' qualifying=', nqual, ' floor=', nmin
            if( nqual >= nstates_req )then
                l_success = .true.
                exit
            endif
            ! peel: members of under-populated clusters and particles outside every kernel support
            ! leave the placement mass, so the next placement spends its centers on the retained mass
            do i = 1, nret
                if( st_r%labels(i) < 1 )then
                    retained(idx(i)) = .false.
                else if( .not. qualifies(st_r%labels(i)) )then
                    retained(idx(i)) = .false.
                endif
            end do
            if( count(retained) < nstates_req * nmin )then
                write(logfhandle,'(A)') '>>> FLEX_PCA POPULATION FLOOR: the retained mass can no longer hold the floors'
                exit
            endif
            if( K >= POP_FLOOR_MAX_NSTATES )then
                write(logfhandle,'(A)') '>>> FLEX_PCA POPULATION FLOOR: provision cap reached'
                exit
            endif
            K = min(POP_FLOOR_MAX_NSTATES, K + (nstates_req - nqual))
        end do
        if( .not. l_success )then
            THROW_WARN('flex_pca population floor not reached for every requested state; keeping the most populated clusters')
        endif
        ! clusters ordered by population, descending (K <= POP_FLOOR_MAX_NSTATES, selection sort)
        allocate(order(K))
        order = [(s, s=1,K)]
        do s = 1, K-1
            kbest = s
            do t = s+1, K
                if( occ(order(t)) > occ(order(kbest)) ) kbest = t
            end do
            if( kbest /= s )then
                itmp         = order(s)
                order(s)     = order(kbest)
                order(kbest) = itmp
            endif
        end do
        allocate(kept(nstates_req), deliver(K))
        kept    = order(1:nstates_req)
        deliver = 0
        do s = 1, nstates_req
            deliver(kept(s)) = s
        end do
        ! standardized latent metric on the placement components, for attaching surplus clusters
        allocate(sdv(nk))
        do q = 1, nk
            zbar   = sum(latent%z(:,q)) / real(nptcls,dp)
            sdv(q) = max(sqrt(sum((latent%z(:,q) - zbar)**2) / real(nptcls,dp)), 1.d-12)
        end do
        allocate(states%labels(nptcls), source=0)
        nsurplus = 0
        do i = 1, nret
            s = st_r%labels(i)
            if( s < 1 ) cycle
            if( deliver(s) > 0 )then
                states%labels(idx(i)) = deliver(s)
            else if( qualifies(s) )then
                ! surplus qualifying cluster: real mass, attached to the nearest delivered target
                dbest = huge(1.d0)
                kbest = 1
                do t = 1, nstates_req
                    d2 = 0.d0
                    do q = 1, nk
                        d2 = d2 + ((latent%z(idx(i),q) - real(st_r%targets(q,kept(t)),dp)) / sdv(q))**2
                    end do
                    if( d2 < dbest )then
                        dbest = d2
                        kbest = t
                    endif
                end do
                states%labels(idx(i)) = kbest
                nsurplus       = nsurplus + 1
            endif
        end do
        ! members of dropped clusters, particles outside every kernel support and particles peeled in
        ! earlier rounds receive a uniformly random delivered label
        nrand = 0
        do i = 1, nptcls
            if( states%labels(i) < 1 )then
                states%labels(i) = irnd_uni(nstates_req)
                nrand     = nrand + 1
            endif
        end do
        ! delivered tables: hard-label indicator weights, so the state maps are ordinary
        ! reconstructions of the labelled particles
        allocate(states%weights(nptcls,nstates_req), source=0.)
        do i = 1, nptcls
            states%weights(i,states%labels(i)) = 1.
        end do
        allocate(states%targets(ncomp,nstates_req), states%bandwidths(nstates_req), states%neff(nstates_req))
        do s = 1, nstates_req
            states%targets(:,s)  = st_r%targets(:,kept(s))
            states%bandwidths(s) = st_r%bandwidths(kept(s))
            states%neff(s)       = real(count(states%labels == s))
        end do
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA POPULATION FLOOR delivered states=', nstates_req, &
            &' surplus-attached=', nsurplus, ' randomized=', nrand
        do s = 1, nstates_req
            write(logfhandle,'(A,I3,A,I9,A,I9)') '>>>   state=', s, '  particles=', nint(states%neff(s)), '  floor=', nmin
            if( nint(states%neff(s)) < nmin ) THROW_WARN('flex_pca delivered a state below the population floor')
        end do
        call flush(logfhandle)
        deallocate(retained, idx, occ, qualifies, order, kept, deliver, sdv)
        call st_r%kill; call lat_r%kill
    end subroutine place_states_with_population_floor

    !> The resume embedding path itself ('' when not resuming)
    function infile_path( params, l_resume ) result( f )
        type(parameters), intent(in) :: params
        logical,          intent(in) :: l_resume
        character(len=:), allocatable :: f
        f = ''
        if( l_resume ) f = params%infile%to_char()
    end function infile_path

end module simple_flex_pca_state_service
