!@descr: flex_pca state inference: latent deconvolution, state placement, occupancy pruning and the derived state settings
!!
!! Everything between the delivered embedding and the weight table the reconstruction consumes:
!! the measurement-error deconvolution of the latents (and its resume adoption), placement of
!! the state kernels (population-floor variant here; the kernel/GMM placement in
!! simple_flex_pca_weights), the occupancy floor, and the derived box/neff/state-count rules.
!! Pure functions of arrays and files: no project, no builder, no rounds.
module simple_flex_pca_state_service
use simple_core_module_api
use simple_flex_pca_records, only: flex_selection, flex_fit_model, flex_latent, flex_state_set
use simple_parameters,           only: parameters
use simple_srch_sort_loc,        only: hpsort
use simple_flex_pca_deconv,      only: calibrate_noise_scale, deconvolve_latent
use simple_flex_pca_embedding_io, only: read_deconv_block, append_deconv_block
use simple_flex_pca_weights,     only: build_covariance_state_weights
implicit none
private
#include "simple_local_flags.inc"

public :: apply_latent_deconvolution, prune_underpopulated_states, place_states_with_population_floor
public :: auto_box_crop, auto_min_neff, auto_state_count
public :: infile_path, infile_dir
public :: FLEX_AUTO_K_START, FLEX_AUTO_K_MIN, AUTO_NSTATES

!> Cap and floor of auto_state_count (tester only; preimage_auto's ceiling is AUTO_NSTATES). The cap
!! is bounded by cost: gate 2 of the state merge compares K(K-1)/2 map pairs.
integer, parameter :: FLEX_AUTO_K_START = 32
integer, parameter :: FLEX_AUTO_K_MIN   = 8
!> Nyquist margin for a derived box_crop; columns are selected inside that band.
real,    parameter :: FLEX_AUTO_BOX_SAFETY = 1.25
!> Occupancy share of an equal split. Paired with an SNR term because neither alone reproduces both
!! validation datasets.
real,    parameter :: FLEX_AUTO_NEFF_OCCUPANCY = 0.10
! npreimages is a PROVISION CEILING, not a target: state placement lays down that many kernels and
! the two-gate merge collapses the indistinct ones, so the recovered K is only ever <= it.
! preimage_auto=yes raises that ceiling to AUTO_NSTATES and turns the merge on, since over-provisioning
! is the only regime in which the merge can recover K at all.
!> provision cap of the population floor (min_state_frac > 0, refine3D_states flex=yes); independent of AUTO_NSTATES
integer, parameter :: AUTO_NSTATES = 8
!> provision cap of the population floor (min_state_frac > 0, refine3D_states flex=yes); independent of AUTO_NSTATES
integer, parameter :: POP_FLOOR_MAX_NSTATES = 32

contains

    !> Calibrate the per-particle noise (from the even/odd half solutions when the run has them,
    !! else from the scale file the original run wrote) and replace z / precision by the posterior
    !! means / precisions under the deconvolved mixture prior.
    subroutine apply_latent_deconvolution( latent, model, sel, applied, labels, resume, adopted, srcdir, srcfile )
        type(flex_latent),    intent(inout) :: latent  !< z and precision deconvolved in place; zhalf consumed
        type(flex_fit_model), intent(in)    :: model
        type(flex_selection), intent(in)    :: sel
        integer, allocatable, intent(inout) :: labels(:)   !< mixture component per particle
        logical,  intent(in)    :: resume
        logical,  intent(out)   :: adopted
        !> directory of the embedding a resume was given (infile): the deconvolved cache and the labels are
        !! looked up there when the run directory has none, so a states-only resume in a fresh directory
        !! adopts the original run's deconvolution instead of re-running the K ladder
        character(len=*), intent(in) :: srcdir
        character(len=*), intent(in) :: srcfile   !< the infile itself: its trailing deconvolved block is adopted first
        logical,  intent(out)   :: applied
        real(dp) :: prior(model%ncomp), a_comp(model%ncomp), noise_scale
        integer :: q, k_deconv, u_ns, io_ns
        logical  :: l_resume
        character(len=:), allocatable :: ns_fname
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
                        logical  :: l_found
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
            open(newunit=u_ns, file='flex_pca_noise_scale.txt', status='replace', action='write')
            write(u_ns,'(ES16.8)') noise_scale
            close(u_ns)
            deallocate(latent%zhalf)
        else
            noise_scale = 1.d0
            ! the calibrated noise scale of the fit: here, or next to infile (a resume in a fresh directory
            ! used to silently fall back to 1.0 and deconvolve differently from the fit)
            ns_fname = 'flex_pca_noise_scale.txt'
            if( .not. file_exists(ns_fname) )then
                if( len_trim(srcdir) > 0 )then
                    if( file_exists(trim(srcdir)//'/flex_pca_noise_scale.txt') ) ns_fname = trim(srcdir)//'/flex_pca_noise_scale.txt'
                endif
            endif
            open(newunit=u_ns, file=ns_fname, status='old', action='read', iostat=io_ns)
            if( io_ns == 0 )then
                read(u_ns,*,iostat=io_ns) noise_scale
                close(u_ns)
            endif
            if( io_ns /= 0 ) noise_scale = 1.d0
            write(logfhandle,'(A,A,A,F8.3)') '>>> FLEX_PCA resumed embedding: noise scale from ', ns_fname, &
                &' (1.0 when absent) =', noise_scale
        endif
        call deconvolve_latent(latent%z, latent%precision, prior, sel%nptcls, model%ncomp, noise_scale, 16, k_deconv, &
            &prior_fname='flex_pca_deconv_prior.txt', labels_fname='flex_pca_deconv_labels.txt', pinds=sel%pinds, &
            &labels_out=labels)
        applied = .true.
        ! the deconvolved coordinates join the raw cache of this run as a trailing block (one file)
        if( file_exists('flex_pca_embedding.bin') )then
            if( allocated(labels) )then
                call append_deconv_block('flex_pca_embedding.bin', sel%nptcls, model%ncomp, latent%z, latent%precision, labels, noise_scale)
            else
                call append_deconv_block('flex_pca_embedding.bin', sel%nptcls, model%ncomp, latent%z, latent%precision, noise_scale=noise_scale)
            endif
        else
            write(logfhandle,'(A)') '>>> FLEX_PCA deconvolved coordinates not cached: no flex_pca_embedding.bin in this directory'
        endif
    end subroutine apply_latent_deconvolution

    !> Drop every state whose effective sample size is below min_neff and compact the per-state
    !! arrays, so the reconstruction, the bandwidth CV and the merge all see only states that can
    !! support a map. Particles whose argmax state is dropped become unassigned (label 0) and feed no
    !! map -- they are, by construction, the particles the placement could not commit anywhere. At
    !! least two states always survive: with fewer the run has no heterogeneity to deliver.
    subroutine prune_underpopulated_states( min_neff, states )
        type(flex_state_set), intent(inout) :: states
        integer :: nptcls
        integer,              intent(in) :: min_neff
        logical,  allocatable :: keep(:)
        integer,  allocatable :: map(:), occ(:), ord(:)
        real,     allocatable :: w2(:,:), t2(:,:), b2(:), n2(:)
        real(dp), allocatable :: d2(:,:), f2(:)
        real,     allocatable :: key(:)
        integer :: s, i, nkeep, ndrop, nlost, snew
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
        nlost = 0
        do i = 1, nptcls
            if( states%labels(i) < 1 .or. states%labels(i) > states%nstates ) cycle
            if( keep(states%labels(i)) )then
                states%labels(i) = map(states%labels(i))
            else
                states%labels(i) = 0
                nlost     = nlost + 1
            endif
        end do
        write(logfhandle,'(A,I0,A,F6.2,A,I0,A)') '>>> FLEX_PCA OCCUPANCY FLOOR: ', nlost, ' particles (', &
            &100.0*real(nlost)/real(max(nptcls,1)), '%) lost their state and feed no map; ', nkeep, ' states remain'
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
        logical,  optional,   intent(in)    :: equal_occ
        type(flex_latent)    :: lat_r
        type(flex_state_set) :: st_r
        integer :: nptcls, ncomp, nstates_req
        integer,  intent(in) :: nkern, axis, min_neff
        real,     intent(in) :: min_state_frac
        integer,  parameter :: ROUND_CAP = 8
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

    !> Directory of the resume embedding (params%infile), '' when not resuming or when the path has no
    !! directory part; used to adopt the original run's deconvolved cache from a fresh run directory.
    function infile_dir( params, l_resume ) result( d )
        type(parameters), intent(in) :: params
        logical,          intent(in) :: l_resume
        character(len=:), allocatable :: d, f
        integer :: k
        d = ''
        if( .not. l_resume ) return
        f = params%infile%to_char()
        k = index(f, '/', back=.true.)
        if( k > 1 ) d = f(1:k-1)
    end function infile_dir

end module simple_flex_pca_state_service
