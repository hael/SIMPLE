! particle picker
module simple_phasecorr_picker
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image, only: image
implicit none

public :: init_phasecorr_picker, exec_phasecorr_picker, kill_phasecorr_picker
private

! PEAK STATS INDICES
integer,          parameter   :: SDEV       = 1
integer,          parameter   :: DYNRANGE   = 2
integer,          parameter   :: SSCORE     = 3
! OTHER PARAMS
integer,          parameter   :: NSTAT      = 3
integer,          parameter   :: MAXKMIT    = 20
real,             parameter   :: BOXFRAC    = 0.5
logical,          parameter   :: DEBUG_HERE = .false.
integer,          parameter   :: PICKER_OFFSET_HERE = 3, OFFSET_HWIN = 1
! VARS
type(image)                   :: micrograph, mic_shrunken, mic_saved
type(image),      allocatable :: refs(:)
logical,          allocatable :: selected_peak_positions(:)
real,             allocatable :: corrmat(:,:), peak_stats(:,:)
integer,          allocatable :: peak_positions(:,:)
character(len=:), allocatable :: micname, refsname
character(len=LONGSTRLEN)     :: boxname
integer                       :: ldim(3), ldim_refs(3), ldim_shrink(3)
integer                       :: ntargets
integer                       :: nrefs, nmax, nmax_sel, orig_box
real                          :: smpd, smpd_shrunken
real                          :: msk, hp,lp, distthr, ndev
real                          :: min_rad, max_rad, step_sz

contains

    subroutine init_phasecorr_picker( micfname, refsfname, smpd_in, lp_in, distthr_in, ndev_in, dir_out )
        use simple_procimgstk, only :  clip_imgfile
        character(len=*),           intent(in) :: micfname, refsfname
        real,                       intent(in) :: smpd_in
        real,             optional, intent(in) :: lp_in, distthr_in, ndev_in
        character(len=*), optional, intent(in) :: dir_out
        logical, allocatable :: lmsk(:,:,:)
        type(image)          :: refimg, mskimg
        integer              :: ifoo, iref
        real                 :: sdev_noise
        allocate(micname,  source=trim(micfname))
        allocate(refsname, source=trim(refsfname))
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
        ! modify according to PICKER_SHRINK
        ldim_refs(1)          = round2even(real(ldim_refs(1))/PICKER_SHRINK)
        ldim_refs(2)          = round2even(real(ldim_refs(2))/PICKER_SHRINK)
        ldim_refs(3)          = 1
        ldim_shrink(1)        = round2even(real(ldim(1))/PICKER_SHRINK)
        ldim_shrink(2)        = round2even(real(ldim(2))/PICKER_SHRINK)
        ldim_shrink(3)        = 1
        smpd_shrunken         = PICKER_SHRINK*smpd
        msk                   = real(ldim_refs(1)/2-5)
        msk                   = max(msk, real(ldim_refs(1)/2-2)) ! for tiny particles
        distthr               = BOXFRAC*real(ldim_refs(1))
        if( present(distthr_in) ) distthr = distthr_in / smpd_shrunken ! pixels
        ! read and shrink references
        allocate( refs(nrefs) )
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
        mic_saved = mic_shrunken
        call mic_shrunken%bp(hp, lp)
        call mic_saved%ifft
      end subroutine init_phasecorr_picker

    subroutine exec_phasecorr_picker( boxname_out, nptcls_out )
        character(len=LONGSTRLEN), intent(out) :: boxname_out
        integer,                   intent(out) :: nptcls_out
        real, pointer :: prmat(:,:,:)
        integer       :: i
        real          :: maxv
        call extract_peaks
        if( ntargets==0 )then
            nptcls_out = 0
            return
        endif
        call distance_filter
        call gather_stats
        call one_cluster_clustering
        nptcls_out = count(selected_peak_positions)
        ! bring back coordinates to original sampling
        peak_positions = nint(PICKER_SHRINK*(real(peak_positions)))-orig_box/2
        call write_boxfile
        ! returns absolute path
        call make_relativepath(CWD_GLOB, boxname, boxname_out)
    end subroutine exec_phasecorr_picker

    subroutine extract_peaks
        type(image)              :: mask_img, circ_mask, tmp_mic
        type(image), allocatable :: ptcls(:)
        real,    pointer         :: rmat_phasecorr(:,:,:)
        integer, allocatable     :: labels(:), target_positions(:,:),ref_inds(:,:)
        real,    allocatable     :: target_corrs(:)
        logical, allocatable     :: mask(:,:), l_mask(:,:,:)
        real                     :: means(2), ave, sdev, maxv, minv, msk_shrunken
        integer                  :: xind, yind, i, border, j, l, r, u, d, iref, ithr
        logical :: outside, l_err
        write(logfhandle,'(a)') '>>> EXTRACTING PEAKS'
        write(logfhandle,'(a)') '>>> FOURIER CORRELATIONS'
        allocate(ref_inds(1:ldim_shrink(1), 1:ldim_shrink(2)), source=0)
        call gen_phase_correlation(mic_shrunken,mask_img)
        call mic_shrunken%stats( ave=ave, sdev=sdev, maxv=maxv, minv=minv,mskimg=mask_img)
        call mic_shrunken%get_rmat_ptr(rmat_phasecorr)
        allocate(corrmat(1:ldim_shrink(1),1:ldim_shrink(2)))
        corrmat(1:ldim_shrink(1),1:ldim_shrink(2)) = rmat_phasecorr(1:ldim_shrink(1),1:ldim_shrink(2),1)
        !$omp parallel do default(shared) private(xind,yind,l,r,u,d) proc_bind(close) schedule(static)
        do xind=1,ldim_shrink(1),PICKER_OFFSET_HERE
            l = max(1,xind-OFFSET_HWIN)
            r = min(xind+ OFFSET_HWIN,ldim_shrink(1))
            do yind=1,ldim_shrink(2),PICKER_OFFSET_HERE
                u = max(1,yind- OFFSET_HWIN)
                d = min(yind+OFFSET_HWIN,ldim_shrink(2))
                corrmat(xind,yind) = maxval(rmat_phasecorr(l:r,u:d,1))
            enddo
        enddo
        !$omp end parallel do
        write(logfhandle,'(a)') '>>> BINARIZATION'
        call mic_shrunken%binarize(ave+.8*sdev)
        if( l_err )then
            ! binarization failed because of uniform values
            ntargets = 0
            return
        endif
        border = max(ldim_refs(1)/2,ldim_refs(2)/2)
        rmat_phasecorr(1:border,:,1) = 0. !set to zero the borders
        rmat_phasecorr(ldim_shrink(1)-border:ldim_shrink(1),:,1) = 0. !set to zero the borders
        rmat_phasecorr(:,1:border,1) = 0. !set to zero the borders
        rmat_phasecorr(:,ldim_shrink(2)-border:ldim_shrink(2),1) = 0. !set to zero the borders
        allocate(mask(1:ldim_shrink(1), 1:ldim_shrink(2)), source = .false.)
        ! Select initial peaks
        ntargets = 0
        !$omp parallel do collapse(2) default(shared) private(xind,yind) proc_bind(close) schedule(static) reduction(+:ntargets)
        do xind=1,ldim_shrink(1),PICKER_OFFSET_HERE
            do yind=1,ldim_shrink(2),PICKER_OFFSET_HERE
                if(rmat_phasecorr(xind,yind,1) > 0.5) then
                    ntargets = ntargets + 1
                    mask(xind,yind) = .true.
                endif
            enddo
        enddo
        !$omp end parallel do
        if( ntargets==0 ) return
        allocate( target_corrs(ntargets),target_positions(ntargets,2))
        ntargets = 0
        do xind=1,ldim_shrink(1),PICKER_OFFSET_HERE
            do yind=1,ldim_shrink(2),PICKER_OFFSET_HERE
                if(mask(xind,yind)) then
                    ntargets = ntargets + 1
                    target_positions(ntargets,:) = [xind,yind]
                endif
            enddo
        enddo
        write(logfhandle,'(a)') '>>> REAL SPACE CORRELATIONS'
        allocate(ptcls(max(1,nthr_glob)))
        do i = 1, size(ptcls)
            call ptcls(i)%new(ldim_refs,smpd_shrunken)
        enddo
        do iref = 1,NREFS
            call refs(iref)%bp(hp,lp)
        enddo
        msk_shrunken = msk*real(orig_box)/real(ldim_refs(1))
        call circ_mask%disc(ldim_refs, smpd_shrunken, msk_shrunken, l_mask)
        call circ_mask%kill
        tmp_mic = mic_saved
        call tmp_mic%bp(hp,lp)
        !$omp parallel do default(shared) private(i,outside,xind,yind,ithr) proc_bind(close) schedule(static)
        do i = 1, ntargets
            ithr = omp_get_thread_num()+1
            xind = target_positions(i,1)
            yind = target_positions(i,2)
            call tmp_mic%window_slim([xind,yind]-ldim_refs(1)/2-1, ldim_refs(1), ptcls(ithr), outside)
            target_corrs(i) = ptcls(ithr)%real_corr(refs(ref_inds(xind,yind)),l_mask)
        enddo
        !$omp end parallel do
        write(logfhandle,'(a)') '>>> PEAKS SELECTION'
        call sortmeans(target_corrs, MAXKMIT, means, labels)
        nmax = count(labels == 2)
        if(DEBUG_HERE) then
            write(logfhandle,'(a,1x,I7)') 'peak positions initially identified:         ', ntargets
            write(logfhandle,'(a,1x,I7)') 'peak positions after sortmeans:              ',nmax
         endif
        ! get peak positions
        allocate( peak_positions(nmax,2) )
        peak_positions = 0
        nmax = 0
        do i=1,ntargets
            if( labels(i) == 2 )then
                nmax = nmax + 1
                peak_positions(nmax,:) = target_positions(i,:)
            endif
        end do
        ! cleanup
        call mask_img%kill
        call tmp_mic%kill
        do i = 1, size(ptcls)
            call ptcls(i)%kill
        enddo
        deallocate(ptcls)

    contains

        ! Reference generation and Phase Correlation calculation
        ! FORMULA: phasecorr = ifft2(fft2(field).*conj(fft2(reference)));
        subroutine gen_phase_correlation(field,mask)
            type(image),         intent(inout) :: field
            type(image),optional,intent(inout) :: mask
            type(image)   :: phasecorr, aux, ref_ext
            real, pointer :: mask_rmat(:,:,:)
            real    :: v
            integer :: iref, border
            border =  max(ldim_refs(1),ldim_refs(2))
            call phasecorr%new(ldim_shrink, smpd_shrunken)
            call aux%new(ldim_shrink, smpd_shrunken)
            call aux%set_ft(.true.)
            call ref_ext%new(ldim_shrink, smpd_shrunken)
            call field%fft
            do iref = 1, nrefs
                call refs(iref)%pad(ref_ext, 0.) ! zero padding
                call ref_ext%fft
                call field%phase_corr(ref_ext,aux,lp,border=max(ldim_refs(1)/2,ldim_refs(2)/2)) !phase correlation
                if(iref > 1) then
                    !$omp parallel do collapse(2) schedule(static) default(shared) private(i,j,v) proc_bind(close)
                    do i = 1, ldim_shrink(1)
                        do j = 1, ldim_shrink(2)
                            v = aux%get([i,j,1])
                            if( v > phasecorr%get([i,j,1]))then
                                call phasecorr%set([i,j,1], v)
                                ref_inds(i,j) = iref
                            endif
                         enddo
                    enddo
                    !$omp end parallel do
                else
                    phasecorr = aux
                    ref_inds  = 1
                endif
                call aux%fft
            enddo
            call field%copy(phasecorr)
            if(present(mask)) then
                call mask%new(ldim_shrink, smpd_shrunken)
                call mask%get_rmat_ptr(mask_rmat)
                mask_rmat(border+1:ldim_shrink(1)-border,border+1:ldim_shrink(2)-border,1)=1.
            endif
            call phasecorr%kill
            call aux%kill
            call ref_ext%kill
        end subroutine gen_phase_correlation
    end subroutine extract_peaks

    subroutine distance_filter
        integer :: ipeak, jpeak, ipos(2), jpos(2), loc(1)
        real    :: dist, corr, dist_sq, distthr_sq
        logical, allocatable :: mask(:)
        real,    allocatable :: corrs(:)
        write(logfhandle,'(a)') '>>> DISTANCE FILTERING'
        allocate( mask(nmax), corrs(nmax), selected_peak_positions(nmax) )
        selected_peak_positions = .true.
        distthr_sq = nint(distthr*distthr)
        !$omp parallel do schedule(static) default(shared) private(ipeak,jpeak,ipos,jpos,dist,mask,loc,corr) proc_bind(close)
        do ipeak=1,nmax
            ipos = peak_positions(ipeak,:)
            corr = -HUGE(corr)
            mask = .false.
            do jpeak=1,nmax
                jpos = peak_positions(jpeak,:)
                ! dist = euclid(real(ipos),real(jpos))
                dist_sq = sum((ipos-jpos)**2)
                if( dist_sq < distthr_sq ) then
                  mask(jpeak) = .true.
                  ! find best match in the neigh
                  if(corrmat(jpos(1),jpos(2)) > corr) then
                    corr   = corrmat(jpos(1),jpos(2))
                    loc(1) = jpeak
                  endif
                endif
              end do
            ! eliminate all but the best
            mask(loc(1)) = .false.
            where( mask )
                selected_peak_positions = .false.
            end where
        end do
        !$omp end parallel do
        nmax_sel = count(selected_peak_positions)
        if(DEBUG_HERE) write(logfhandle,'(a,1x,I7)') 'peak positions left after distance filtering: ', nmax_sel
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
        allocate( peak_stats(nmax_sel,NSTAT) )
        do ipeak=1,nmax
            if( selected_peak_positions(ipeak) )then
                cnt = cnt + 1
                call mic_saved%window_slim(peak_positions(ipeak,:)-ldim_refs(1)/2,&
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

    ! In the old school picker there are 4 stats. In picker there are 3.
    ! The disregarded stat is the correlation (CC2REF). The effect it has, if put
    ! back, is that the picking becomes more permissive.
    subroutine one_cluster_clustering
        real, allocatable :: dmat(:,:)
        integer           :: i_median, i, j, cnt, ipeak
        real              :: ddev
        allocate(dmat(nmax_sel,nmax_sel), source=0.)
        do i=1,nmax_sel - 1
            do j=i + 1,nmax_sel
                dmat(i,j) = euclid(peak_stats(i,:), peak_stats(j,:))
                dmat(j,i) = dmat(i,j)
            end do
        end do
        call dev_from_dmat( dmat, i_median, ddev )
        cnt = 0
        do ipeak=1,nmax
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
        nmax_sel = count(selected_peak_positions)
        write(logfhandle,'(a,1x,I5)') 'peak positions after one cluster clustering: ', nmax_sel
    end subroutine one_cluster_clustering

    subroutine write_boxfile
        integer :: funit, ipeak,iostat
        call fopen(funit, status='REPLACE', action='WRITE', file=trim(adjustl(boxname)),iostat=iostat)
        call fileiochk('phasecorr_picker; write_boxfile ', iostat)
        do ipeak=1,nmax
            if( selected_peak_positions(ipeak) )then
                write(funit,'(I7,I7,I7,I7,I7)') peak_positions(ipeak,1),&
                peak_positions(ipeak,2), orig_box, orig_box, -3
            endif
        end do
        call fclose(funit)
    end subroutine write_boxfile

    subroutine kill_phasecorr_picker
        integer :: iref
        if( allocated(micname) )then
            call micrograph%kill
            call mic_shrunken%kill
            call mic_saved%kill
            if( allocated(selected_peak_positions) ) deallocate(selected_peak_positions)
            if( allocated(corrmat)                 ) deallocate(corrmat)
            if( allocated(peak_positions)          ) deallocate(peak_positions)
            if( allocated(peak_stats)              ) deallocate(peak_stats)
            if( allocated(micname)                 ) deallocate(micname)
            if(allocated(refsname)                 ) deallocate(refsname)
            do iref=1,nrefs
                call refs(iref)%kill
            end do
            deallocate(refs)
        endif
    end subroutine kill_phasecorr_picker
end module simple_phasecorr_picker
