!@descr: the species call of detect_atomic_species: the class fit on the atoms of detect_atoms' density model, posteriors with their errors, the _species files
! The design is doc/implementation_notes/planned/atomic_density_model.md (sections 2 and 12): the classes are fitted on
! ln I_i = ln S_k + beta (c_i - 12) with one beta for every class, every atom carrying the variance tau_I(c_i)**2 + v_i
! (tau_I the coordination spread by coordination bin, measured within the classes, v_i the measurement variance of
! ln I_i);
! the products are those of the former discover_species of detect_atoms (doc/implementation_notes/completed/
! species_discovery_completed.md, section 6) with two columns more per atom, the posterior's standard error and the
! coordination spread it carries. The atoms come from _ATMS.pdb, their fields from the _atoms.csv and _model.txt beside it.
module simple_atomic_species
use simple_core_module_api, only: string, logfhandle, PI, simple_exception, file_exists, get_fbody, basename, fopen,&
    &fclose, fileiochk, int2str, hpsort, pearsn_serial, find_ldim_nptcls
use simple_defs,         only: dp
use simple_image,        only: image
use simple_atoms,        only: atoms
use simple_nano_species, only: fit_species_mixture, class_posterior, class_separation
implicit none

public :: call_atomic_species
private
#include "simple_local_flags.inc"

real,    parameter :: CN_CLOSE_PACKED = 12.
real,    parameter :: V_UNDEFINED     = 1.     ! variance of ln I of an atom without a measurement spread (or I <= 0)
real,    parameter :: INT_FLOOR_FRAC  = 1.e-3  ! ln I of an atom with I <= 0 is taken at this fraction of the median I
real,    parameter :: APER_NN_FRAC    = 0.5    ! aperture radius of the half maps' intensities, in d_NN
integer, parameter :: NSHELL_MAX      = 5
integer, parameter :: NSHELL_ATOMS    = 20     ! a radial shell holds at least this many atoms
real,    parameter :: LOW_SNR         = 6.5    ! shells below this intensity signal-to-noise are flagged
real,    parameter :: INNER_FRAC      = 0.3    ! interior atoms for the coordination deficit
real,    parameter :: SURF_CONFINED   = 0.8    ! a class with more of its atoms in the outer zone is flagged
integer, parameter :: NTAU_ITER       = 10     ! refits of the coordination spread within the classes, at most
real,    parameter :: TAU_TOL         = 1.e-3  ! the coordination spread has settled when no bin changes by more
integer, parameter :: NTAU_MIN        = 5      ! a coordination bin with fewer atoms takes the spread of all atoms

contains

    ! the species call on the atoms of pdbfile (the _ATMS.pdb of detect_atoms) in the map vol1 of sampling smpd;
    ! names: the species brightest first (their number is the class count), else the classes are X1, X2, ...;
    ! nspecies: the class count, 0 from the data; with vol_even and vol_odd the measurement variance of ln I comes from
    ! the half maps' aperture intensities, and their label agreement is reported
    subroutine call_atomic_species( vol1, pdbfile, smpd, nspecies_in, names, vol_even, vol_odd )
        class(string),              intent(in) :: vol1, pdbfile
        real,                       intent(in) :: smpd
        integer,                    intent(in) :: nspecies_in
        character(len=2), optional, intent(in) :: names(:)
        class(string),    optional, intent(in) :: vol_even, vol_odd
        type(atoms)  :: atms
        type(string) :: stem, csvfile, modelfile, fbody, errmsg
        real,    allocatable :: xyz(:,:), ints(:), bfacs(:), sd_lnint(:), valid(:), aper(:), x(:), s2(:), v(:), tau(:)
        real,    allocatable :: cov(:), post(:,:), mu(:), var(:), w(:), bic(:), sd_post(:), rad(:), tau_bins(:), tau_new(:)
        real,    allocatable :: ae(:), ao(:)
        integer, allocatable :: coord(:), labels(:), bin_lo(:)
        logical, allocatable :: strong(:), adm(:)
        real    :: beta, d_nn, prune_rad, cen(3), cen_strong(3), med_int, s_i, half_agree, half_corr, s_half
        integer :: n, K, j, nspecies, nstrong, iter, ntau_refits
        logical :: l_halves
        write(logfhandle,'(A)') '>>> ATOMIC SPECIES'
        l_halves = present(vol_even) .and. present(vol_odd)
        ! the atoms and their table, found beside the PDB as the atomicity prior finds it
        call atms%new(pdbfile)
        n = atms%get_n()
        if( n < 2 ) THROW_HARD('the species call needs two atoms or more; call_atomic_species')
        stem = get_fbody(pdbfile, 'pdb')
        if( .not. stem%ends_with_substr('_ATMS') ) THROW_HARD('pdbfile must be the _ATMS.pdb of detect_atoms')
        stem      = stem%to_char([1, stem%strlen_trim() - 5])
        csvfile   = stem//'_atoms.csv'
        modelfile = stem//'_model.txt'
        if( .not. file_exists(csvfile) )then
            errmsg = 'no atom table '//csvfile%to_char()//' beside the PDB; run detect_atoms first'
            THROW_HARD(errmsg%to_char())
        endif
        if( .not. file_exists(modelfile) )then
            errmsg = 'no model parameters '//modelfile%to_char()//' beside the PDB; run detect_atoms first'
            THROW_HARD(errmsg%to_char())
        endif
        fbody = basename(stem)
        call read_atoms_table(csvfile, n, xyz, ints, bfacs, coord, sd_lnint, valid, aper, strong)
        call read_model(modelfile, d_nn, prune_rad, bin_lo, tau_bins)
        nstrong = count(strong)
        ! ln I with its noise: the coordination spread of its bin and the measurement variance of ln I
        med_int = median_of(pack(ints, ints > 0.))
        allocate(x(n), v(n), s2(n))
        do j = 1,n
            if( ints(j) > 0. )then
                x(j) = log(ints(j))
            else
                x(j) = log(INT_FLOOR_FRAC * med_int)
            endif
            v(j) = V_UNDEFINED
            if( ints(j) > 0. .and. sd_lnint(j) > 0. ) v(j) = sd_lnint(j)**2
        enddo
        ! the half maps: aperture intensities scaled to the fitted intensities by the map's own apertures; the noise of the
        ! average, half their difference, gives the measurement variance of ln I in place of the table's
        if( l_halves )then
            call half_apertures(ae, ao)
            s_half = sdev(0.5 * (ae - ao))
            do j = 1,n
                if( ints(j) > 0. ) v(j) = (s_half / ints(j))**2
            enddo
        endif
        cov = real(coord) - CN_CLOSE_PACKED
        nspecies = nspecies_in
        if( present(names) ) nspecies = size(names)
        allocate(labels(n))
        ! the tau_I of _model.txt is measured by detect_atoms under one species: on a particle of several it holds their
        ! separation. The class call measures it within its classes: a first fit on the measurement variances alone,
        ! then tau_I per coordination bin from the residuals about the classes, refitted until it settles; tau is the
        ! spread of the last fit
        s2 = v
        call fit_species_mixture(x, s2, nspecies, K, labels, post, mu, var, bic, adm, cov=cov, slope=beta, w=w)
        tau = within_class_tau()
        do iter = 1,NTAU_ITER
            s2 = tau**2 + v
            call fit_species_mixture(x, s2, nspecies, K, labels, post, mu, var, bic, adm, cov=cov, slope=beta, w=w)
            tau_new = within_class_tau()
            if( maxval(abs(tau_new - tau)) <= TAU_TOL .or. iter == NTAU_ITER ) exit
            tau = tau_new
        enddo
        ntau_refits = min(iter, NTAU_ITER)
        write(logfhandle,'(A,I0,A)') 'coordination spread within the classes: ', ntau_refits, ' refits'
        ! the posterior of every atom, with the standard error of its class's posterior from its measurement noise
        allocate(sd_post(n))
        do j = 1,n
            call class_posterior(x(j), s2(j), v(j), mu + beta * cov(j), var, w, post(j,:), labels(j), sd_post(j))
        enddo
        cen = sum(xyz, dim=2) / real(n)
        cen_strong = cen
        if( nstrong > 0 ) cen_strong = sum(xyz, dim=2, mask=spread(strong, 1, 3)) / real(nstrong)
        allocate(rad(n))
        do j = 1,n
            rad(j) = sqrt(sum((xyz(:,j) - cen)**2))
        enddo
        s_i = median_of(ints * sqrt(v))
        if( l_halves ) call halfmap_agreement
        call write_species_report
        call write_radial_profiles
        call write_species_table
        call write_species_pdb
        write(logfhandle,'(A,I0,A,I0,A,I0)') 'strong atoms: ', nstrong, ', weak atoms: ', n - nstrong, ', classes: ', K
        write(logfhandle,'(A)') '>>> ATOMIC SPECIES, COMPLETED'
        call atms%kill

    contains

        real function sdev( vv )
            real, intent(in) :: vv(:)
            real :: m
            m    = sum(vv) / real(size(vv))
            sdev = sqrt(sum((vv - m)**2) / real(max(size(vv) - 1, 1)))
        end function sdev

        real function median_of( vv )
            real, intent(in) :: vv(:)
            real, allocatable :: vs(:)
            median_of = 0.
            if( size(vv) < 1 ) return
            vs = vv
            call hpsort(vs)
            if( mod(size(vs), 2) == 0 )then
                median_of = 0.5 * (vs(size(vs)/2) + vs(size(vs)/2 + 1))
            else
                median_of = vs(size(vs)/2 + 1)
            endif
        end function median_of

        ! the coordination bin of the density model at coordination c (bin 1 the highest, each from bin_lo up; the last
        ! takes the rest)
        elemental integer function bin_of( c )
            integer, intent(in) :: c
            integer :: kb
            bin_of = max(size(bin_lo), 1)
            do kb = 1,size(bin_lo)
                if( c >= bin_lo(kb) )then
                    bin_of = kb
                    return
                endif
            enddo
        end function bin_of

        ! the coordination spread of bin kb in the final fit, from its first atom
        real function tau_in_bin( kb )
            integer, intent(in) :: kb
            integer :: jj
            tau_in_bin = 0.
            do jj = 1,n
                if( bin_of(coord(jj)) == kb )then
                    tau_in_bin = tau(jj)
                    return
                endif
            enddo
        end function tau_in_bin

        ! the coordination spread per atom from the current fit: in each coordination bin, the posterior-weighted mean
        ! squared residual of ln I about the class models less the mean measurement variance, over the atoms with a
        ! measurement; a bin with fewer than NTAU_MIN of them takes the value of all atoms
        function within_class_tau() result( tau_ )
            real     :: tau_(n)
            real(dp) :: r2(n), tb(max(size(bin_lo), 1)), all_t2
            logical  :: meas(n)
            integer  :: jj, kc, kb, ib(n)
            meas = v < V_UNDEFINED
            do jj = 1,n
                r2(jj) = 0._dp
                do kc = 1,K
                    r2(jj) = r2(jj) + real(post(jj,kc), dp) * (real(x(jj), dp) - real(mu(kc) + beta * cov(jj), dp))**2
                enddo
                ib(jj) = bin_of(coord(jj))
            enddo
            all_t2 = 0._dp
            if( count(meas) > 0 ) all_t2 = max(0._dp, sum(r2, mask=meas) / real(count(meas), dp) -&
                &sum(real(v, dp), mask=meas) / real(count(meas), dp))
            do kb = 1,size(tb)
                tb(kb) = all_t2
                if( count(meas .and. ib == kb) < NTAU_MIN ) cycle
                tb(kb) = max(0._dp, sum(r2, mask=meas .and. ib == kb) / real(count(meas .and. ib == kb), dp) -&
                    &sum(real(v, dp), mask=meas .and. ib == kb) / real(count(meas .and. ib == kb), dp))
            enddo
            do jj = 1,n
                tau_(jj) = real(sqrt(tb(ib(jj))))
            enddo
        end function within_class_tau

        ! ln I the class fit expects of atom j in class kc, and that intensity
        real function expected_int( jj, kc )
            integer, intent(in) :: jj, kc
            expected_int = exp(mu(kc) + beta * cov(jj))
        end function expected_int

        ! the class/environment intensity averaged over the atom's posterior
        real function posterior_model_intensity( jj )
            integer, intent(in) :: jj
            integer :: kc
            posterior_model_intensity = 0.
            do kc = 1,K
                posterior_model_intensity = posterior_model_intensity + post(jj,kc) * expected_int(jj, kc)
            enddo
        end function posterior_model_intensity

        ! the fitted intensity relative to the fully occupied class/environment model, capped at one
        real function species_occupancy( jj )
            integer, intent(in) :: jj
            species_occupancy = min(1., max(0., ints(jj) / expected_int(jj, labels(jj))))
        end function species_occupancy

        ! aperture sums within APER_NN_FRAC d_NN of every atom in the map and the half maps; the halves scaled by the median
        ! ratio of the fitted intensities to the map's apertures
        subroutine half_apertures( ae_, ao_ )
            real, allocatable, intent(out) :: ae_(:), ao_(:)
            real, allocatable :: am(:), ratio(:)
            real :: scale
            am  = apertures(vol1)
            ae_ = apertures(vol_even)
            ao_ = apertures(vol_odd)
            ratio = pack(ints / max(am, tiny(1.)), ints > 0. .and. am > 0.)
            if( size(ratio) < 1 ) THROW_HARD('no atom with a positive intensity and aperture; half_apertures')
            scale = median_of(ratio)
            ae_ = scale * ae_
            ao_ = scale * ao_
        end subroutine half_apertures

        function apertures( fname ) result( ap )
            class(string), intent(in) :: fname
            real, allocatable :: ap(:), rmat(:,:,:)
            type(image) :: img
            integer :: ldim(3), nptcls, jj, lo(3), hi(3), ii, kk, ll
            real    :: c(3), rv
            call find_ldim_nptcls(fname, ldim, nptcls)
            call img%new(ldim, smpd)
            call img%read(fname)
            rmat = img%get_rmat()
            call img%kill
            rv = APER_NN_FRAC * d_nn / smpd
            allocate(ap(n), source=0.)
            do jj = 1,n
                c  = xyz(:,jj) / smpd + 1.
                lo = max(1, floor(c - rv))
                hi = min(ldim, ceiling(c + rv))
                do ll = lo(3),hi(3)
                    do kk = lo(2),hi(2)
                        do ii = lo(1),hi(1)
                            if( sum((real([ii,kk,ll]) - c)**2) <= rv * rv ) ap(jj) = ap(jj) + rmat(ii,kk,ll)
                        enddo
                    enddo
                enddo
            enddo
        end function apertures

        ! the half maps' intensities classified with the final mixture, and their correlation
        subroutine halfmap_agreement
            real    :: pe(K), po(K), sde, sdo
            integer :: jj, nagree, le, lo_
            nagree = 0
            do jj = 1,n
                call class_posterior(log(max(ae(jj), INT_FLOOR_FRAC * med_int)), s2(jj), v(jj), mu + beta * cov(jj), var, w,&
                    &pe, le, sde)
                call class_posterior(log(max(ao(jj), INT_FLOOR_FRAC * med_int)), s2(jj), v(jj), mu + beta * cov(jj), var, w,&
                    &po, lo_, sdo)
                if( le == lo_ ) nagree = nagree + 1
            enddo
            half_agree = real(nagree) / real(n)
            half_corr  = pearsn_serial(ae, ao)
        end subroutine halfmap_agreement

        ! equal-count radial shells of class kc: first and last index into the radius order
        subroutine class_shells( kc, order, nsh, bnd )
            integer,              intent(in)  :: kc
            integer, allocatable, intent(out) :: order(:), bnd(:)
            integer,              intent(out) :: nsh
            real,    allocatable :: rk(:)
            integer :: jj, nk
            order = pack([(jj, jj=1,n)], labels == kc)
            nk    = size(order)
            nsh   = max(1, min(NSHELL_MAX, nk / NSHELL_ATOMS))
            if( nk > 0 )then
                rk = rad(order)
                call hpsort(rk, order)
            endif
            allocate(bnd(nsh+1))
            do jj = 1,nsh+1
                bnd(jj) = ((jj-1) * nk) / nsh
            enddo
        end subroutine class_shells

        ! per class and radial shell: the width and the intensity of the density model
        subroutine write_radial_profiles
            integer, allocatable :: order(:), bnd(:), sel(:)
            real,    allocatable :: sig2(:), ii_(:)
            integer :: funit, kc, ish, nsh, m
            real    :: rms2, snr
            call fopen(funit, file=fbody//'_species_radial.csv', status='REPLACE', action='WRITE')
            write(funit,'(A)') 'CLASS,SHELL,NATOMS,MEAN_RADIUS,RMS_SIGMA,SE_SIGMA,B,MEAN_INT,SE_INT,MEAN_COORD,INT_SNR,LOW_SNR'
            do kc = 1,K
                call class_shells(kc, order, nsh, bnd)
                do ish = 1,nsh
                    m = bnd(ish+1) - bnd(ish)
                    if( m < 1 ) cycle
                    sel  = order(bnd(ish)+1:bnd(ish+1))
                    sig2 = sqrt(bfacs(sel) / (8. * PI**2))
                    ii_  = ints(sel)
                    rms2 = sqrt(sum(sig2**2) / real(m))
                    snr  = 0.
                    if( s_i > 0. ) snr = exp(mu(kc)) / s_i
                    write(funit,'(I0,",",I0,",",I0,8(",",ES16.8),",",I0)') kc, ish, m, sum(rad(sel)) / real(m), rms2,&
                        &sdev(sig2) / sqrt(real(m)), 8. * PI**2 * rms2**2, sum(ii_) / real(m), sdev(ii_) / sqrt(real(m)),&
                        &real(sum(coord(sel))) / real(m), snr, merge(1, 0, snr < LOW_SNR)
                enddo
            enddo
            call fclose(funit)
        end subroutine write_radial_profiles

        ! one row per atom in the order of _ATMS.pdb and _atoms.csv; one posterior column per class, POST1 to POSTK; then
        ! the three columns of the atomicity prior (simple_rec3D_pcg_strategy): the fully occupied intensity of the atom's
        ! class at its coordination, the same averaged over the posterior, and the fitted intensity over the former, capped
        ! at one; then the standard error of the atom's class posterior and the coordination spread its variance carries
        subroutine write_species_table
            type(string) :: header, fmt
            integer      :: funit, jj, kc
            header = 'INDEX,X,Y,Z,RADIUS,COORD,INT,BFAC,APER_INT,CLASS'
            do kc = 1,K
                header = header//',POST'//int2str(kc)
            enddo
            header = header//',VALID_CORR,MODEL_INTENSITY,POSTERIOR_INTENSITY,OCCUPANCY,SD_POSTERIOR,TAU_COORD'
            fmt    = '(I0,8(",",ES16.8),",",I0,'//int2str(K+6)//'(",",ES16.8))'
            call fopen(funit, file=fbody//'_species.csv', status='REPLACE', action='WRITE')
            write(funit,'(A)') header%to_char()
            do jj = 1,n
                write(funit,fmt%to_char()) jj, xyz(:,jj), rad(jj), real(coord(jj)), ints(jj), bfacs(jj), aper(jj), labels(jj),&
                    &post(jj,:), valid(jj), expected_int(jj, labels(jj)), posterior_model_intensity(jj), species_occupancy(jj),&
                    &sd_post(jj), tau(jj)
            enddo
            call fclose(funit)
            write(logfhandle,'(A)') 'output, atoms with their class:   '//fbody%to_char()//'_species.csv'
        end subroutine write_species_table

        ! every atom with its class: element column and chain by class (the named species, or X1, X2, ..., chain A the
        ! brightest), occupancy the fitted intensity over the class/environment model, B column the model's width
        subroutine write_species_pdb
            type(atoms)      :: out
            character(len=2) :: el
            integer          :: jj, kc
            if( K > 36 ) THROW_HARD('more classes than PDB chains A-Z, 0-9; write_species_pdb')
            call out%new(n, dummy=.true.)
            do jj = 1,n
                kc = labels(jj)
                if( present(names) )then
                    el = names(kc)
                else
                    el = 'X'//int2str(kc)
                endif
                call out%set_name(     jj, el//'  ')
                call out%set_element(  jj, el)
                if( kc <= 26 )then
                    call out%set_chain(jj, achar(iachar('A') + kc - 1))
                else
                    call out%set_chain(jj, achar(iachar('0') + kc - 27))
                endif
                call out%set_coord(    jj, xyz(:,jj))
                call out%set_num(      jj, jj)
                call out%set_resnum(   jj, jj)
                call out%set_occupancy(jj, species_occupancy(jj))
                call out%set_beta(     jj, bfacs(jj))
            enddo
            call out%writepdb(fbody//'_species.pdb')
            call out%kill
            write(logfhandle,'(A)') 'output, atoms with their class:   '//fbody%to_char()//'_species.pdb'
        end subroutine write_species_pdb

        subroutine write_species_report
            integer, allocatable :: order(:), bnd(:)
            real,    allocatable :: r(:)
            integer :: funit, kc, l, ninner, nsh, ish, nlow, nin
            real    :: cn_inner, misclass, cn_out, frac_out
            character(len=9) :: source
            call fopen(funit, file=fbody//'_species.txt', status='REPLACE', action='WRITE')
            write(funit,'(A)') '# detect_atomic_species: the class call on the density model of detect_atoms, ln I = ln S_k + beta (c - 12)'
            call write_kv(funit, 'K', int2str(K))
            call write_kv(funit, 'nspecies_fixed', int2str(nspecies))
            do kc = 1,K
                call write_kv(funit, 'class_intensity_'//int2str(kc), rstr(exp(mu(kc))))
                call write_kv(funit, 'class_ratio_'//int2str(kc),     rstr(exp(mu(kc) - mu(1))))
                call write_kv(funit, 'class_fraction_'//int2str(kc),  rstr(real(count(labels == kc)) / real(n)))
                call write_kv(funit, 'class_atoms_'//int2str(kc),     int2str(count(labels == kc)))
                call write_kv(funit, 'class_sd_ln_'//int2str(kc),     rstr(sqrt(var(kc))))
            enddo
            do kc = 1,K-1
                call write_kv(funit, 'separation_D_'//int2str(kc)//'_'//int2str(kc+1),&
                    &rstr(class_separation(mu(kc), var(kc) + median_of(s2), mu(kc+1), var(kc+1) + median_of(s2))))
            enddo
            ! the class intensities are those of an atom with CN_CLOSE_PACKED contacts; beta is the change of ln I per contact
            call write_kv(funit, 'beta', rstr(beta))
            ! the coordination spread per bin, measured within the classes, beside the one-species value of detect_atoms
            call write_kv(funit, 'tau_I_refits', int2str(ntau_refits))
            do l = 1,size(bin_lo)
                call write_kv(funit, 'tau_I_bin_'//int2str(l)//'_cn_min', int2str(bin_lo(l)))
                if( any(bin_of(coord) == l) ) call write_kv(funit, 'tau_I_bin_'//int2str(l), rstr(tau_in_bin(l)))
                call write_kv(funit, 'tau_I_model_bin_'//int2str(l), rstr(tau_bins(l)))
            enddo
            do l = 1,size(bic)
                call write_kv(funit, 'bic_K'//int2str(l),        rstr(bic(l)))
                call write_kv(funit, 'admissible_K'//int2str(l), merge('yes', 'no ', adm(l)))
            enddo
            misclass = sum(1. - maxval(post, dim=2)) / real(n)
            call write_kv(funit, 'expected_misclassification', rstr(misclass))
            call write_kv(funit, 'median_sd_posterior', rstr(median_of(sd_post)))
            call write_kv(funit, 'max_sd_posterior', rstr(maxval(sd_post)))
            call write_kv(funit, 'd_NN', rstr(d_nn))
            source = 'model'
            if( l_halves ) source = 'half_maps'
            call write_kv(funit, 'noise_source', trim(source))
            call write_kv(funit, 's_I', rstr(s_i))
            call write_kv(funit, 'strong_atoms', int2str(nstrong))
            call write_kv(funit, 'weak_atoms', int2str(n - nstrong))
            call write_kv(funit, 'total_atoms', int2str(n))
            ! coordination of the inner atoms: a deficit flags undetected sites
            r = rad
            order = [(l, l=1,n)]
            call hpsort(r, order)
            ninner   = max(1, nint(INNER_FRAC * real(n)))
            cn_inner = real(sum(coord(order(:ninner)))) / real(ninner)
            call write_kv(funit, 'interior_coordination', rstr(cn_inner))
            call write_kv(funit, 'interior_coordination_deficit', rstr(CN_CLOSE_PACKED - cn_inner))
            nlow = 0
            if( s_i > 0. )then
                do kc = 1,K
                    call class_shells(kc, order, nsh, bnd)
                    do ish = 1,nsh
                        if( bnd(ish+1) <= bnd(ish) ) cycle
                        if( exp(mu(kc)) / s_i < LOW_SNR ) nlow = nlow + 1
                    enddo
                enddo
            endif
            call write_kv(funit, 'low_snr_shells', int2str(nlow))
            ! a class confined to the outer zone of the pruning (about the strong atoms' centre) may be partial occupancy
            do kc = 1,K
                nin    = 0
                cn_out = 0.
                do l = 1,n
                    if( labels(l) /= kc ) cycle
                    if( sqrt(sum((xyz(:,l) - cen_strong)**2)) > prune_rad )then
                        nin    = nin + 1
                        cn_out = cn_out + real(coord(l))
                    endif
                enddo
                frac_out = real(nin) / real(max(count(labels == kc), 1))
                call write_kv(funit, 'class_outer_zone_fraction_'//int2str(kc), rstr(frac_out))
                call write_kv(funit, 'class_outer_zone_mean_coord_'//int2str(kc), rstr(cn_out / real(max(nin, 1))))
                call write_kv(funit, 'class_surface_confined_'//int2str(kc), merge('yes', 'no ', frac_out > SURF_CONFINED))
            enddo
            if( l_halves )then
                call write_kv(funit, 'halfmap_label_agreement', rstr(half_agree))
                call write_kv(funit, 'halfmap_intensity_corr', rstr(half_corr))
            endif
            call fclose(funit)
            write(logfhandle,'(A)') 'output, class report:             '//fbody%to_char()//'_species.txt'
        end subroutine write_species_report

        subroutine write_kv( funit, key, val )
            integer,          intent(in) :: funit
            character(len=*), intent(in) :: key, val
            write(funit,'(A,T34,A)') key, '= '//trim(val)
        end subroutine write_kv

        function rstr( val ) result( str )
            real, intent(in) :: val
            character(len=16) :: str
            write(str,'(ES16.8)') val
            str = adjustl(str)
        end function rstr

    end subroutine call_atomic_species

    ! the columns of _atoms.csv the species call reads, found by name: positions (A), intensity, width, coordination,
    ! the measurement spread of ln I, valid_corr, the aperture intensity, and whether the atom is strong
    subroutine read_atoms_table( fname, n, xyz, ints, bfacs, coord, sd_lnint, valid, aper, strong )
        class(string),        intent(in)  :: fname
        integer,              intent(in)  :: n
        real,    allocatable, intent(out) :: xyz(:,:), ints(:), bfacs(:), sd_lnint(:), valid(:), aper(:)
        integer, allocatable, intent(out) :: coord(:)
        logical, allocatable, intent(out) :: strong(:)
        character(len=1024) :: line
        character(len=64)   :: fields(64)
        integer :: funit, ios, nf, j, cx, cy, cz, ci, cb, cc, cs, cv, ca, cl
        call fopen(funit, file=fname, status='old', action='read', iostat=ios)
        call fileiochk('read_atoms_table opening '//fname%to_char(), ios)
        read(funit,'(A)', iostat=ios) line
        call split_fields(line, fields, nf)
        cx = col('X');   cy = col('Y');   cz = col('Z');  ci = col('INT'); cb = col('BFAC'); cc = col('CN')
        cs = col('SD_LNINT'); cv = col('VALID_CORR'); ca = col('APER_INT'); cl = col('LABEL')
        allocate(xyz(3,n), ints(n), bfacs(n), sd_lnint(n), valid(n), aper(n), coord(n), strong(n))
        do j = 1,n
            read(funit,'(A)', iostat=ios) line
            if( ios /= 0 ) THROW_HARD('the atom table has fewer rows than the PDB has atoms; read_atoms_table')
            call split_fields(line, fields, nf)
            read(fields(cx),*) xyz(1,j)
            read(fields(cy),*) xyz(2,j)
            read(fields(cz),*) xyz(3,j)
            read(fields(ci),*) ints(j)
            read(fields(cb),*) bfacs(j)
            read(fields(cc),*) coord(j)
            read(fields(cs),*) sd_lnint(j)
            read(fields(cv),*) valid(j)
            read(fields(ca),*) aper(j)
            strong(j) = trim(adjustl(fields(cl))) == 'strong'
        enddo
        read(funit,'(A)', iostat=ios) line
        if( ios == 0 .and. len_trim(line) > 0 ) THROW_HARD('the atom table has more rows than the PDB has atoms; read_atoms_table')
        call fclose(funit)

    contains

        integer function col( name )
            character(len=*), intent(in) :: name
            type(string) :: msg
            integer :: k
            col = 0
            do k = 1,nf
                if( trim(adjustl(fields(k))) == name ) col = k
            enddo
            if( col == 0 )then
                msg = 'the atom table has no column '//name//'; read_atoms_table'
                THROW_HARD(msg%to_char())
            endif
        end function col

    end subroutine read_atoms_table

    ! d_NN, the pruning zone radius and the coordination spread by bin (lowest coordination of each bin and its spread;
    ! none without profiles) of a _model.txt
    subroutine read_model( fname, d_nn, prune_rad, bin_lo, tau_bins )
        class(string),        intent(in)  :: fname
        real,                 intent(out) :: d_nn, prune_rad
        integer, allocatable, intent(out) :: bin_lo(:)
        real,    allocatable, intent(out) :: tau_bins(:)
        character(len=256) :: line, key
        integer :: funit, ios, ieq, nb, kb
        real    :: val
        d_nn      = 0.
        prune_rad = huge(1.)
        nb        = 0
        allocate(bin_lo(16), source=0)
        allocate(tau_bins(16), source=0.)
        call fopen(funit, file=fname, status='old', action='read', iostat=ios)
        call fileiochk('read_model opening '//fname%to_char(), ios)
        do
            read(funit,'(A)', iostat=ios) line
            if( ios /= 0 ) exit
            ieq = index(line, '=')
            if( ieq < 2 ) cycle
            key = adjustl(line(:ieq-1))
            read(line(ieq+1:),*,iostat=ios) val
            if( ios /= 0 ) cycle
            select case(trim(key))
                case('d_NN')
                    d_nn = val
                case('prune_zone_radius')
                    prune_rad = val
                case('tau_I_bins')
                    nb = min(nint(val), 16)
                case default
                    if( index(key, 'tau_I_bin_') == 1 )then
                        if( index(key, '_cn_min') > 0 )then
                            read(key(len('tau_I_bin_')+1:index(key, '_cn_min')-1),*,iostat=ios) kb
                            if( ios == 0 .and. kb >= 1 .and. kb <= 16 ) bin_lo(kb) = nint(val)
                        else
                            read(key(len('tau_I_bin_')+1:),*,iostat=ios) kb
                            if( ios == 0 .and. kb >= 1 .and. kb <= 16 ) tau_bins(kb) = val
                        endif
                    endif
            end select
        enddo
        call fclose(funit)
        if( .not. d_nn > 0. ) THROW_HARD('the model parameters carry no d_NN; read_model')
        bin_lo   = bin_lo(:nb)
        tau_bins = tau_bins(:nb)
    end subroutine read_model

    ! the comma-separated fields of a line
    subroutine split_fields( line, fields, nf )
        character(len=*), intent(in)  :: line
        character(len=*), intent(out) :: fields(:)
        integer,          intent(out) :: nf
        integer :: j, start
        nf    = 0
        start = 1
        fields = ''
        do j = 1,len_trim(line) + 1
            if( j > len_trim(line) .or. line(j:j) == ',' )then
                nf = nf + 1
                if( nf > size(fields) ) THROW_HARD('too many fields in a table row; split_fields')
                fields(nf) = line(start:j-1)
                start = j + 1
            endif
        enddo
    end subroutine split_fields

end module simple_atomic_species
