module simple_picker
!$ use omp_lib
!$ use omp_lib_kinds
use simple_defs
use simple_jiffys
use simple_math,             only: hpsort
use simple_image,            only: image
use simple_math,             only: euclid
use simple_denspeak_cluster, only: denspeak_cluster
implicit none

public :: init_picker, pick_particles, kill_picker
private

integer                       :: ldim(3), ifoo, ldim_shrink(3), alloc_stat, ntargets, nx, ny
integer                       :: xind, yind, cnt, ldim_refs(3), nrefs, ntop_targets, iref, funit
integer                       :: io_stat, nbest_targets, nworst_targets, pop_peak, pop_fuzz, kmit
integer                       :: npeaks, ipeak, jpeak, ncls, icls, posmax(2), ipos(2), jpos(2)
integer                       :: trs, xrange(2), yrange(2), noutside, xsh, ysh, loc(1), i, j, sz
integer                       :: offset, orig_box
real                          :: corr, smpd_shrunken, corr_peak, corr_fuzz, corrmax, corrmin
real                          :: dist_peak, dist_fuzz, corr_peak_new, corr_fuzz_new, dcrit, dist
real                          :: icorr, jcorr, corr_best, smpd, msk, shrink, lp, distthr, dcrit_rel
type(denspeak_cluster)        :: densclust
type(image)                   :: micrograph, mic_shrunken, ptcl_target
type(image),      allocatable :: refs(:)
real,             allocatable :: corrmat(:,:), corrs(:), corrs_packed(:), dmat(:,:), shcorrs(:,:)
real,             allocatable :: dmat_part(:,:), sxx(:)
real,             allocatable :: clust_diams(:)
logical,          allocatable :: corrmask(:,:)
integer,          allocatable :: peak_positions(:,:), labels(:), centers(:), refmat(:,:), pinds(:)
character(len=:), allocatable :: micname, refsname
character(len=STDLEN)         :: boxname
integer,          parameter   :: MAXKMIT=20
real,             parameter   :: BOXFRAC=0.67
logical,          parameter   :: debug=.false.
logical,          parameter   :: hard_dist_filt=.false.

contains

    subroutine init_picker( micfname, refsfname, smpd_in, msk_in, offset_in,&
        shrink_in, lp_in, dcrit_rel_in, distthr_in )
        character(len=*),  intent(in) :: micfname, refsfname
        real,              intent(in) :: smpd_in, msk_in
        integer, optional, intent(in) :: offset_in
        real,    optional, intent(in) :: shrink_in, lp_in, dcrit_rel_in, distthr_in
        allocate(micname,  source=trim(micfname))
        allocate(refsname, source=trim(refsfname))
        boxname   = remove_abspath(fname_new_ext(micname,'box'))
        smpd      = smpd_in
        msk       = msk_in
        offset    = 3
        if( present(offset_in)    ) offset    = offset_in
        shrink    = 4.0
        if( present(shrink_in)    ) shrink    = shrink_in
        lp        = 20.0
        if( present(lp_in)        ) lp        = lp_in
        dcrit_rel = 0.5
        if( present(dcrit_rel_in) ) dcrit_rel = dcrit_rel_in
        ! read micrograph
        call find_ldim_nptcls(micname, ldim, ifoo)
        call micrograph%new(ldim, smpd)
        call micrograph%read(micname)
        ! find out reference dimensions
        call find_ldim_nptcls(refsname, ldim_refs, nrefs)
        orig_box     = nint(shrink)*ldim_refs(1)
        ldim_refs(3) = 1 ! correct 4 stupid mrc convention
        ! set constants
        ldim_shrink(1) = nint(real(ldim(1))/shrink)
        ldim_shrink(2) = nint(real(ldim(2))/shrink)
        ldim_shrink(3) = 1
        nx             = ldim_shrink(1)-ldim_refs(1)
        ny             = ldim_shrink(2)-ldim_refs(2)
        smpd_shrunken  = shrink*smpd
        msk            = msk/shrink                           ! mask parameter need to be modulated by shrink
        distthr = BOXFRAC*real(ldim_refs(1))
        if( present(distthr_in) ) distthr = distthr_in/shrink ! inputted dist thresh need to be modulated by shrink
        ! read references
        allocate( refs(nrefs), sxx(nrefs), stat=alloc_stat )
        call alloc_err( "In: simple_picker :: init_picker, 1", alloc_stat)
        do iref=1,nrefs
            call refs(iref)%new(ldim_refs, smpd_shrunken)
            call refs(iref)%read(refsname, iref)
            call refs(iref)%mask(msk, 'hard')
            call refs(iref)%prenorm4real_corr(sxx(iref))
        end do
        ! pre-process micrograph
        call micrograph%fwd_ft
        call micrograph%bp(0., lp)
        call mic_shrunken%new(ldim_shrink, smpd_shrunken)
        call micrograph%clip(mic_shrunken)
        call mic_shrunken%bwd_ft
        if( debug ) call mic_shrunken%write('shrunken.mrc')
    end subroutine init_picker

    subroutine pick_particles
        call extr_imgs_and_calc_corrs
        call identify_peaks
        call denspeak_clust_and_distfilt
        call shift_refinement
    end subroutine pick_particles

    subroutine extr_imgs_and_calc_corrs
        ! extract the images & calculate correlations
        allocate( corrmat(0:nx,0:ny), refmat(0:nx,0:ny), corrmask(0:nx,0:ny), corrs(nrefs), stat=alloc_stat )
        call alloc_err( 'In: simple_picker, 2', alloc_stat)
        corrmat  = 0.
        corrmask = .false.
        corrs    = 0.
        ntargets = 0
        ! if( debug .and. file_exists('picker_corrmat.bin') )then
        !     do xind=0,nx,offset
        !         do yind=0,ny,offset
        !             ntargets = ntargets + 1
        !             corrmask(xind,yind) = .true.
        !         end do
        !     end do
        !     ! read in the previously calculated correlations (to save time when debugging)
        !     funit = get_fileunit()
        !     open(unit=funit, status='OLD', action='READ', file='picker_corrmat.bin', access='STREAM')
        !     read(unit=funit, pos=1, iostat=io_stat) corrmat
        !     ! check if the read was successful
        !     if( io_stat .ne. 0 )then
        !         write(*,'(a,i0,2a)') '**ERROR(simple_picker): I/O error ',&
        !         io_stat, ' when reading file picker_corrmat.bin'
        !         stop 'I/O error'
        !     endif
        !     close(funit)
        !     ! read in the previously calculated refmat (to save time when debugging)
        !     open(unit=funit, status='OLD', action='READ', file='picker_refmat.bin', access='STREAM')
        !     read(unit=funit, pos=1, iostat=io_stat) refmat
        !     ! check if the read was successful
        !     if( io_stat .ne. 0 )then
        !         write(*,'(a,i0,2a)') '**ERROR(simple_picker): I/O error ',&
        !         io_stat, ' when reading file picker_refmat.bin'
        !         stop 'I/O error'
        !     endif
        !     close(funit)
        ! else
            corrmax = -1.
            corrmin = 1.
            do xind=0,nx,offset
                do yind=0,ny,offset
                    ntargets = ntargets + 1
                    call mic_shrunken%window([xind,yind], ldim_refs(1), ptcl_target)
                    !$omp parallel do schedule(auto) default(shared) private(iref)
                    do iref=1,nrefs
                        corrs(iref) = refs(iref)%real_corr_prenorm(ptcl_target, sxx(iref))
                    end do
                    !$omp end parallel do
                    loc = maxloc(corrs)
                    refmat(xind,yind)   = loc(1)
                    if( refmat(xind,yind) == 0 ) stop 'refmat cannot have zero elements'
                    corrmat(xind,yind)  = corrs(loc(1))
                    corrmask(xind,yind) = .true.
                    if( corrmat(xind,yind) > corrmax ) corrmax = corrmat(xind,yind)
                    if( corrmat(xind,yind) < corrmin ) corrmin = corrmat(xind,yind)
                end do
            end do
            ! if( debug )then
            !     funit = get_fileunit()
            !     open(unit=funit, status='REPLACE', action='WRITE', file='picker_corrmat.bin', access='STREAM')
            !     write(unit=funit,pos=1,iostat=io_stat) corrmat
            !     ! check if the write was successful
            !     if( io_stat .ne. 0 )then
            !         write(*,'(a,i0,2a)') '**ERROR(simple_picker): I/O error ',&
            !         io_stat, ' when writing to picker_corrmat.bin'
            !         stop 'I/O error'
            !     endif
            !     close(funit)
            !     open(unit=funit, status='REPLACE', action='WRITE', file='picker_refmat.bin', access='STREAM')
            !     write(unit=funit,pos=1,iostat=io_stat) refmat
            !     ! check if the write was successful
            !     if( io_stat .ne. 0 )then
            !         write(*,'(a,i0,2a)') '**ERROR(simple_picker): I/O error ',&
            !         io_stat, ' when writing to picker_refmat.bin'
            !         stop 'I/O error'
            !     endif
            !     close(funit)
            ! endif
        ! endif
    end subroutine extr_imgs_and_calc_corrs

    subroutine identify_peaks
        ! pack and sort the corrmat 
        corrs_packed = pack(corrmat, mask=corrmask)
        call hpsort(size(corrs_packed), corrs_packed)
        ! quantize it into two quantas (peak/fuzz)
        nbest_targets  = nint(real(ntargets)/1000.)
        nworst_targets = ntargets-nbest_targets
        corr_peak      = sum(corrs_packed(nworst_targets+1:))/real(nbest_targets)
        corr_fuzz      = sum(corrs_packed(:nworst_targets))/real(nworst_targets)
        do kmit=1,MAXKMIT
            pop_peak      = 0
            pop_fuzz      = 0
            corr_peak_new = 0.
            corr_fuzz_new = 0.
            do xind=0,nx,offset
                do yind=0,ny,offset
                    dist_peak = sqrt((corr_peak-corrmat(xind,yind))**2.0)
                    dist_fuzz = sqrt((corr_fuzz-corrmat(xind,yind))**2.0)
                    if( dist_peak <= dist_fuzz )then
                        corr_peak_new = corr_peak_new+corrmat(xind,yind)
                        pop_peak      = pop_peak + 1
                    else
                        corr_fuzz_new = corr_fuzz_new+corrmat(xind,yind)
                        pop_fuzz      = pop_fuzz + 1
                    endif
                end do
            end do
            corr_peak = corr_peak_new/real(pop_peak)
            corr_fuzz = corr_fuzz_new/real(pop_fuzz)
        end do
        if( debug ) print *, 'corr_peak: ', corr_peak
        if( debug ) print *, 'corr_fuzz: ', corr_fuzz
        ! count the number of peaks
        npeaks = 0
        do xind=0,nx,offset
            do yind=0,ny,offset
                if( corrmat(xind,yind) >= corr_peak )then
                    npeaks = npeaks + 1
                endif
            end do
        end do
        if( debug ) print *, 'npeaks: ', npeaks
        ! store the peak positions
        allocate( peak_positions(npeaks,2), dmat(npeaks,npeaks), stat=alloc_stat )
        call alloc_err( "In: simple_picker, 3", alloc_stat )
        peak_positions = 0
        dmat           = 0.
        npeaks         = 0
        do xind=0,nx,offset
            do yind=0,ny,offset
                if( corrmat(xind,yind) >= corr_peak )then
                    npeaks = npeaks + 1
                    peak_positions(npeaks,:) = [xind,yind]
                endif
            end do
        end do
        if( debug ) print *, 'stored peak positions'        
    end subroutine identify_peaks

    subroutine denspeak_clust_and_distfilt
        ! create a distance matrix
        !$omp parallel do schedule(auto) default(shared) private(ipeak,jpeak)
        do ipeak=1,npeaks-1
            do jpeak=ipeak+1,npeaks
                dmat(ipeak,jpeak) = euclid(real(peak_positions(ipeak,:)),real(peak_positions(jpeak,:)))
                dmat(jpeak,ipeak) = dmat(ipeak,jpeak)
            end do
        end do
        !$omp end parallel do
        if( debug ) print *, 'calculated distance matrix'
        ! perform density peak clustering
        dcrit = real(ldim_refs(1))*dcrit_rel
        call densclust%new(npeaks, dmat, dcrit)
        call densclust%cluster
        if( debug ) print *, 'clustered correlation density peaks'
        centers = densclust%get_centers()
        labels  = densclust%get_labels()
        ncls    = size(centers)
        write(*,'(i5,1x,a)') ncls, 'clusters of particle coordinates found'
        if( hard_dist_filt )then
            ! measure cluster diameters
            allocate( clust_diams(ncls), stat=alloc_stat )
            call alloc_err( "In: simple_picker, 4", alloc_stat )
            clust_diams = 0.
            do icls=1,ncls
                pinds = from_labels_get_pinds(labels, icls)
                sz = size(pinds)
                if( sz > 1  )then
                    allocate(dmat_part(sz,sz))
                    dmat_part = 0.
                    do i=1,sz-1
                        do j=i+1,sz
                            dmat_part(i,j) = dmat(pinds(i),pinds(j))
                        end do
                    end do
                    clust_diams(icls) = maxval(dmat_part)
                    ! apply cluster diameter distance filter
                    if( clust_diams(icls) > real(ldim_refs(1)) ) centers(icls) = 0
                    deallocate(dmat_part)
                endif
                deallocate(pinds)
            end do
            if( debug ) print *, 'excluded', count(centers == 0), 'clusters based on diameter thresholding'
            if( debug ) print *, ncls-count(centers == 0), 'clusters left'
        endif
        ! re-define the centers to be the corr peak maximum in the cluster
        do icls=1,ncls
            sz = count(labels == icls)
            if( sz > 1 .and. centers(icls) /= 0 )then
                corr_best = -1.
                do ipeak=1,npeaks
                    if( labels(ipeak) == icls )then
                        ipos  = peak_positions(ipeak,:)
                        corr =  corrmat(ipos(1),ipos(2))
                        if( corr > corr_best )then
                            corr_best = corr
                            centers(icls) = ipeak
                        endif
                    endif
                end do
            endif
        end do
        ! apply direct distance filter
        do ipeak=1,ncls-1
            do jpeak=ipeak+1,ncls
                if( centers(ipeak) /= 0 .and. centers(jpeak) /= 0 )then
                    ipos  = peak_positions(centers(ipeak),:)
                    jpos  = peak_positions(centers(jpeak),:)
                    dist  = euclid(real(ipos),real(jpos))
                    if( dist < distthr )then
                        icorr = corrmat(ipos(1),ipos(2))
                        jcorr = corrmat(jpos(1),jpos(2))
                        if( icorr > jcorr )then
                            centers(jpeak) = 0
                        else
                            centers(ipeak) = 0
                        endif
                    endif
                endif
            end do
        end do
        ncls = count(centers /= 0)
        write(*,'(i5,1x,a)') ncls, 'boxes left after direct distance filter'

        contains

            function from_labels_get_pinds( labels, icls ) result( pinds )
                integer, intent(in)  :: labels(:), icls
                integer, allocatable :: pinds(:)
                integer :: pop, ilab, nlabs, cnt
                pop = count(labels == icls)
                if( pop > 0 )then
                    nlabs = size(labels)
                    allocate(pinds(pop))
                    cnt = 0
                    do ilab=1,nlabs
                        if( labels(ilab) == icls )then
                            cnt = cnt + 1
                            pinds(cnt) = ilab
                        endif
                    end do
                endif
            end function from_labels_get_pinds

    end subroutine denspeak_clust_and_distfilt

    subroutine shift_refinement
        if( debug ) print *, 'refining peak positions...'
        ! refine peak positions
        trs = ldim_refs(1)/4
        do icls=1,size(centers)
            if( centers(icls) /= 0 )then
                posmax    = peak_positions(centers(icls),:)
                iref      = refmat(posmax(1),posmax(2))
                xrange(1) = posmax(1)-trs
                xrange(2) = posmax(1)+trs
                yrange(1) = posmax(2)-trs 
                yrange(2) = posmax(2)+trs
                corr_best = -1.
                do xsh=xrange(1),xrange(2)
                    do ysh=yrange(1),yrange(2)
                        noutside = 0
                        call mic_shrunken%window([xsh,ysh], ldim_refs(1), ptcl_target, noutside)
                        if( noutside > 0 )then
                            corr = -1.
                        else
                            corr = refs(iref)%real_corr_prenorm(ptcl_target, sxx(iref))
                        endif
                        if( corr > corr_best )then
                            corr_best = corr
                            posmax    = [xsh,ysh]
                        endif
                    end do
                end do
                peak_positions(centers(icls),:) = posmax
            endif
        end do
        if( debug ) print *, '...done refining peak positions ;-)'
        ! write refined boxfile
        funit = get_fileunit()
        open(unit=funit, status='REPLACE', action='WRITE', file=boxname)
        do icls=1,size(centers)
            if( centers(icls) /= 0 )then
                posmax = peak_positions(centers(icls),:)
                posmax = posmax*nint(shrink)
                write(funit,'(I7,I7,I7,I7,I7)') posmax(1), posmax(2), orig_box, orig_box, -3
            endif
        end do
        close(funit)
    end subroutine shift_refinement

    subroutine kill_picker
        integer :: sz, i
        call densclust%kill
        call micrograph%kill
        call mic_shrunken%kill
        call ptcl_target%kill
        if( allocated(refs) )then
            sz = size(refs)
            do i=1,sz
                call refs(i)%kill
            end do
            deallocate(refs)
        endif
        if( allocated(corrmat)        ) deallocate(corrmat)
        if( allocated(corrs)          ) deallocate(corrs)
        if( allocated(corrs_packed)   ) deallocate(corrs_packed)
        if( allocated(dmat)           ) deallocate(dmat)
        if( allocated(shcorrs)        ) deallocate(shcorrs)
        if( allocated(dmat_part)      ) deallocate(dmat_part)
        if( allocated(sxx)            ) deallocate(sxx)
        if( allocated(clust_diams)    ) deallocate(clust_diams)
        if( allocated(corrmask)       ) deallocate(corrmask)
        if( allocated(peak_positions) ) deallocate(peak_positions)
        if( allocated(labels)         ) deallocate(labels)
        if( allocated(centers)        ) deallocate(centers)
        if( allocated(refmat)         ) deallocate(refmat)
        if( allocated(pinds)          ) deallocate(pinds)
        if( allocated(micname)        ) deallocate(micname)
        if( allocated(refsname)       ) deallocate(refsname)
    end subroutine kill_picker

end module simple_picker