! particle picker
module simple_phasecorr_picker
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,        only: image
implicit none

public :: init_phasecorr_picker, exec_phasecorr_picker, kill_phasecorr_picker
private

! PEAK STATS INDICES
integer,          parameter   :: SDEV     = 2
integer,          parameter   :: DYNRANGE = 3
integer,          parameter   :: SSCORE   = 4
real,             parameter   :: PICKER_SHRINK_HERE   = 2.
! OTHER PARAMS
integer,          parameter   :: NSTAT   = 4
integer,          parameter   :: MAXKMIT = 20
real,             parameter   :: BOXFRAC = 0.5
logical,          parameter   :: DOWRITEIMGS = .true.

! VARS
type(image)                   :: micrograph, mic_shrunken
type(image),      allocatable :: refs(:)
logical,          allocatable :: selected_peak_positions(:)
real,             allocatable :: corrmat(:,:), peak_stats(:,:)
integer,          allocatable :: peak_positions(:,:)
character(len=:), allocatable :: micname, refsname
character(len=LONGSTRLEN)     :: boxname
integer                       :: ldim(3), ldim_refs(3), ldim_shrink(3)
integer                       :: ntargets
integer                       :: nrefs, npeaks, npeaks_sel, orig_box
real                          :: smpd, smpd_shrunken
real                          :: msk, hp,lp, distthr, ndev

contains

    subroutine init_phasecorr_picker( micfname, refsfname, smpd_in, lp_in, distthr_in, ndev_in, dir_out )
        use simple_procimgfile, only :  clip_imgfile
        character(len=*),           intent(in) :: micfname, refsfname
        real,                       intent(in) :: smpd_in
        real,             optional, intent(in) :: lp_in, distthr_in, ndev_in
        character(len=*), optional, intent(in) :: dir_out
        logical, allocatable :: lmsk(:,:,:)
        type(image)          :: refimg, mskimg
        integer              :: ifoo, iref
        real                 :: sdev_noise
        allocate(micname,  source=trim(micfname), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('picker;init, 1',alloc_stat)
        allocate(refsname, source=trim(refsfname), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('picker;init, 2')
        boxname = basename( fname_new_ext(micname,'box') )
        if( present(dir_out) )boxname = trim(dir_out)//trim(boxname)
        smpd = smpd_in
        lp   = 20.0
        if( present(lp_in) ) lp = lp_in
        ndev = 2.0
        if( present(ndev_in)) ndev = ndev_in
        ! read micrograph
        call find_ldim_nptcls(micname, ldim, ifoo)
        call micrograph%new(ldim, smpd)
        call micrograph%read(micname)
        ! find reference dimensions
        call find_ldim_nptcls(refsname, ldim_refs, nrefs)
        orig_box              = ldim_refs(1)
        ! modify according to PICKER_SHRINK_HERE
        ldim_refs(1)          = round2even(real(ldim_refs(1))/PICKER_SHRINK_HERE)
        ldim_refs(2)          = round2even(real(ldim_refs(2))/PICKER_SHRINK_HERE)
        ldim_refs(3)          = 1
        ldim_shrink(1)        = round2even(real(ldim(1))/PICKER_SHRINK_HERE)
        ldim_shrink(2)        = round2even(real(ldim(2))/PICKER_SHRINK_HERE)
        ldim_shrink(3)        = 1
        smpd_shrunken         = PICKER_SHRINK_HERE*smpd
        msk                   = real(ldim_refs(1)/2-5)
        msk                   = max(msk, real(ldim_refs(1)/2-2)) ! for tiny particles
        distthr               = BOXFRAC*real(ldim_refs(1))
        if( present(distthr_in) ) distthr = distthr_in/PICKER_SHRINK_HERE
        ! read and shrink references
        allocate( refs(nrefs), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk( "In: simple_picker :: init_picker, refs etc. ",alloc_stat)
        call mskimg%disc([orig_box,orig_box,1], smpd, msk*real(orig_box)/real(ldim_refs(1)), lmsk)
        call mskimg%kill
        do iref=1,nrefs
            call refs(iref)%new(ldim_refs, smpd_shrunken)
            call refimg%new([orig_box,orig_box,1], smpd)
            call refimg%read(refsname, iref)
            call refimg%noise_norm(lmsk, sdev_noise)
            call refimg%mask(msk*real(orig_box)/real(ldim_refs(1)), 'soft')
            call refimg%fft()
            call refimg%clip(refs(iref))
            call refs(iref)%ifft()
            call refs(iref)%mask(msk, 'soft')
        end do
        call refimg%kill
        ! pre-process micrograph
        call micrograph%fft()
        call mic_shrunken%new(ldim_shrink, smpd_shrunken)
        call mic_shrunken%set_ft(.true.)
        call micrograph%clip(mic_shrunken)
        hp = real(ldim_shrink(1) / 2) * smpd_shrunken
        call mic_shrunken%fft
        call mic_shrunken%bp(hp, lp)
    end subroutine init_phasecorr_picker

    subroutine exec_phasecorr_picker( boxname_out, nptcls_out )
        character(len=LONGSTRLEN), intent(out) :: boxname_out
        integer,                   intent(out) :: nptcls_out
        call extract_peaks
        call distance_filter
        call gather_stats
        call one_cluster_clustering
        nptcls_out = count(selected_peak_positions)
        ! bring back coordinates to original sampling
        peak_positions = nint(PICKER_SHRINK_HERE*(real(peak_positions)))-orig_box/2
        call write_boxfile
        ! returns absolute path
        call make_relativepath(CWD_GLOB, boxname, boxname_out)
    end subroutine exec_phasecorr_picker

    subroutine extract_peaks
        real    :: means(2)
        integer :: xind, yind, alloc_stat, i
        integer, allocatable :: labels(:), target_positions(:,:)
        real,    allocatable :: target_corrs(:)
        real,    pointer     :: rmat_phasecorr(:,:,:)
        type(image)          :: mask_img
        logical, allocatable :: mask(:,:)
        real    :: ave, sdev, maxv, minv
        integer :: border
        write(logfhandle,'(a)') '>>> EXTRACTING PEAKS'
        call gen_phase_correlation(mic_shrunken,mask_img)
        call mic_shrunken%stats( ave=ave, sdev=sdev, maxv=maxv, minv=minv,mskimg=mask_img)
        call mic_shrunken%get_rmat_ptr(rmat_phasecorr)
        allocate(corrmat(1:ldim_shrink(1),1:ldim_shrink(2)))
        corrmat(1:ldim_shrink(1),1:ldim_shrink(2)) = rmat_phasecorr(1:ldim_shrink(1),1:ldim_shrink(2),1)
        call mic_shrunken%binarize(ave+.8*sdev)
        border = max(ldim_refs(1)/2,ldim_refs(2)/2)
        rmat_phasecorr(1:border,:,1) = 0. !set to zero the borders
        rmat_phasecorr(ldim_shrink(1)-border:ldim_shrink(1),:,1) = 0. !set to zero the borders
        rmat_phasecorr(:,1:border,1) = 0. !set to zero the borders
        rmat_phasecorr(:,ldim_shrink(2)-border:ldim_shrink(2),1) = 0. !set to zero the borders
        if(DOWRITEIMGS) call mic_shrunken%write(PATH_HERE//basename(trim(micname))//'_shrunken_bin.mrc')
        allocate(mask(1:ldim_shrink(1), 1:ldim_shrink(2)), source = .false.)
        ntargets = 0
        !$omp parallel do collapse(2) default(shared) private(xind,yind) proc_bind(close) schedule(static) reduction(+:ntargets)
        do xind=1,ldim_shrink(1),PICKER_OFFSET
            do yind=1,ldim_shrink(2),PICKER_OFFSET
                if(rmat_phasecorr(xind,yind,1) > 0.5) then
                    ntargets = ntargets + 1
                    mask(xind,yind) = .true.
                endif
            enddo
        enddo
        !$omp end parallel do
        allocate( target_corrs(ntargets),target_positions(ntargets,2))
        ntargets = 0
        do xind=1,ldim_shrink(1),PICKER_OFFSET
            do yind=1,ldim_shrink(2),PICKER_OFFSET
                if(mask(xind,yind)) then
                    ntargets = ntargets + 1
                    target_positions(ntargets,:) = [xind,yind]
                endif
            enddo
        enddo
        target_corrs = pack(corrmat,mask)
        call sortmeans(target_corrs, MAXKMIT, means, labels)  !TO IMPROVE
        npeaks = count(labels == 2)
        allocate( peak_positions(npeaks,2),  stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk( 'In: simple_picker :: gen_corr_peaks, 2',alloc_stat)
        peak_positions = 0
        npeaks = 0
        do i=1,ntargets
            if( labels(i) == 2 )then
                npeaks = npeaks + 1
                peak_positions(npeaks,:) = target_positions(i,:)
            endif
        end do
        ! cleanup
        call mask_img%kill
    contains
        ! Reference generation and Phase Correlation calculation
        ! FORMULA: phasecorr = ifft2(fft2(field).*conj(fft2(reference)));
        subroutine gen_phase_correlation(field,mask)
            type(image),         intent(inout) :: field
            type(image),optional,intent(inout) :: mask
            type(image)   :: phasecorr, aux, ref_ext
            real, pointer :: mask_rmat(:,:,:)
            integer :: iref
            integer :: border
            border =  max(ldim_refs(1),ldim_refs(2))
            call phasecorr%new(ldim_shrink, smpd_shrunken)
            call aux%new(ldim_shrink, smpd_shrunken)
            call aux%set_ft(.true.)
            call ref_ext%new(ldim_shrink, smpd_shrunken)
            call field%fft
            do iref = 1, nrefs
                call refs(iref)%pad(ref_ext, 0.) ! zero padding
                ! call refs(iref)%fft ! WASSSAAP
                ! call refs(iref)%bp(hp,lp)
                call ref_ext%fft
                call field%phase_corr(ref_ext,aux,lp,border=max(ldim_refs(1)/2,ldim_refs(2)/2)) !phase correlation
                if(iref > 1) then
                    call max_image(phasecorr,phasecorr,aux) !save in phasecorr the maximum value between previous phasecorr and new phasecorr
                else
                    phasecorr   = aux
                endif
                call aux%fft
            enddo
            call field%copy(phasecorr)
            if(DOWRITEIMGS) call field%write(PATH_HERE//basename(trim(micname))//'MaxValPhaseCorr.mrc')
            if(present(mask)) then
                call mask%new(ldim_shrink, smpd_shrunken)
                call mask%get_rmat_ptr(mask_rmat)
                mask_rmat(border+1:ldim_shrink(1)-border,border+1:ldim_shrink(2)-border,1)=1.
            endif
            call phasecorr%kill
            call aux%kill
            call ref_ext%kill
        end subroutine gen_phase_correlation

        ! This subroutine creates an image in which each pixel value
        ! is the max value between the gray level in img1 and img2
        subroutine max_image(img,img1,img2)
            type(image), intent(inout) :: img !output image
            type(image), intent(in)    :: img1, img2
            real, pointer :: rmat(:,:,:)
            real, pointer :: rmat1(:,:,:), rmat2(:,:,:) !matrices for img1 and img2
            integer :: i, j
            call img%get_rmat_ptr(rmat)
            call img1%get_rmat_ptr(rmat1)
            call img2%get_rmat_ptr(rmat2)
            !$omp parallel do collapse(2) schedule(static) default(shared) private(i,j) proc_bind(close)
            do i = 1, ldim_shrink(1)
                do j = 1, ldim_shrink(2)
                    rmat(i,j,1) = max(rmat1(i,j,1), rmat2(i,j,1))
                 enddo
            enddo
            !$omp end parallel do
        end subroutine max_image
    end subroutine extract_peaks

    subroutine distance_filter
        integer :: ipeak, jpeak, ipos(2), jpos(2), loc(1)
        real    :: dist
        logical, allocatable :: mask(:)
        real,    allocatable :: corrs(:)
        write(logfhandle,'(a)') '>>> DISTANCE FILTERING'
        allocate( mask(npeaks), corrs(npeaks), selected_peak_positions(npeaks), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk( 'In: simple_picker :: distance_filter',alloc_stat)
        selected_peak_positions = .true.
        do ipeak=1,npeaks
            ipos = peak_positions(ipeak,:)
            mask = .false.
            !$omp parallel do schedule(static) default(shared) private(jpeak,jpos,dist) proc_bind(close)
            do jpeak=1,npeaks
                jpos = peak_positions(jpeak,:)
                dist = euclid(real(ipos),real(jpos))
                if( dist < distthr ) mask(jpeak) = .true.
                corrs(jpeak) = corrmat(jpos(1),jpos(2))
            end do
            !$omp end parallel do
            ! find best match in the neigh
            loc = maxloc(corrs, mask=mask)
            ! eliminate all but the best
            mask(loc(1)) = .false.
            where( mask )
                selected_peak_positions = .false.
            end where
        end do
        npeaks_sel = count(selected_peak_positions)
        write(logfhandle,'(a,1x,I5)') 'peak positions left after distance filtering: ', npeaks_sel
    end subroutine distance_filter

    subroutine gather_stats
        integer           :: ipeak, cnt, istat
        logical           :: outside
        real              :: ave, maxv, minv
        real, allocatable :: spec(:)
        type(image)       :: ptcl_target
        write(logfhandle,'(a)') '>>> GATHERING REMAINING STATS'
        call ptcl_target%new(ldim_refs, smpd_shrunken)
        cnt = 0
        allocate( peak_stats(npeaks_sel,NSTAT) )
        do ipeak=1,npeaks
            if( selected_peak_positions(ipeak) )then
                cnt = cnt + 1
                call mic_shrunken%window_slim(peak_positions(ipeak,:)-ldim_refs(1),&
                    &ldim_refs(1), ptcl_target, outside)
                call ptcl_target%spectrum('power', spec)
                peak_stats(cnt,SSCORE) = sum(spec)/real(size(spec))
                call ptcl_target%stats('background', ave, peak_stats(cnt,SDEV), maxv, minv)
                peak_stats(cnt,DYNRANGE) = maxv - minv
            endif
        end do
        ! min/max normalise to get all vars on equal footing
        do istat=1,NSTAT
            call normalize_minmax(peak_stats(:,istat))
        end do
        call ptcl_target%kill()
    end subroutine gather_stats

    subroutine one_cluster_clustering
        real, allocatable :: dmat(:,:)
        integer           :: i_median, i, j, cnt, ipeak
        real              :: ddev
        allocate(dmat(npeaks_sel,npeaks_sel), source=0., stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('picker::one_cluster_clustering dmat ',alloc_stat)
        do i=1,npeaks_sel - 1
            do j=i + 1,npeaks_sel
                dmat(i,j) = euclid(peak_stats(i,:), peak_stats(j,:))
                dmat(j,i) = dmat(i,j)
            end do
        end do
        call dev_from_dmat( dmat, i_median, ddev )
        cnt = 0
        do ipeak=1,npeaks
            if( selected_peak_positions(ipeak) )then
                cnt = cnt + 1
                if( dmat(i_median,cnt) <=  ndev * ddev )then
                    ! we are keeping this one
                else
                    ! we are removing this one
                    selected_peak_positions(ipeak) = .false.
                endif
            endif
        end do
        npeaks_sel = count(selected_peak_positions)
        write(logfhandle,'(a,1x,I5)') 'peak positions left after one cluster clustering: ', npeaks_sel
    end subroutine one_cluster_clustering

    subroutine write_boxfile
        integer :: funit, ipeak,iostat
        call fopen(funit, status='REPLACE', action='WRITE', file=trim(adjustl(boxname)),iostat=iostat)
        call fileiochk('phasecorr_picker; write_boxfile ', iostat)
        do ipeak=1,npeaks
            if( selected_peak_positions(ipeak) )then
                write(funit,'(I7,I7,I7,I7,I7)') peak_positions(ipeak,1),&
                peak_positions(ipeak,2), orig_box, orig_box, -3
            endif
        end do
        call fclose(funit,errmsg='picker; write_boxfile end')
    end subroutine write_boxfile

    subroutine kill_phasecorr_picker
        integer :: iref
        if( allocated(micname) )then
            call micrograph%kill
            call mic_shrunken%kill
            deallocate(selected_peak_positions,corrmat, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('phasecorr_picker kill, 1',alloc_stat)
            deallocate(peak_positions,micname,refsname,peak_stats, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('phasecorr_picker kill, 2',alloc_stat)
            do iref=1,nrefs
                call refs(iref)%kill
            end do
            deallocate(refs, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('phasecorr_picker; kill 3',alloc_stat)
        endif
    end subroutine kill_phasecorr_picker
end module simple_phasecorr_picker
