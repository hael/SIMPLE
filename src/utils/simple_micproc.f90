! operations on micrographs
module simple_micproc
include 'simple_lib.f08'
use simple_image,     only: image
use simple_image_bin, only: image_bin
implicit none
#include "simple_local_flags.inc"

contains

    function sample_filetab( filetab, nsample ) result(filetab_smpl)
        class(string), intent(in) :: filetab(:)
        integer,       intent(in) :: nsample
        type(string), allocatable :: filetab_smpl(:)
        integer, allocatable :: inds(:)
        type(ran_tabu) :: rt
        integer :: n
        n   = size(filetab)
        rt = ran_tabu(n)
        allocate(inds(nsample), source=0)
        call rt%ne_ran_iarr(inds)
        call rt%kill
        allocate(filetab_smpl(nsample))
        filetab_smpl = filetab(inds)
    end function sample_filetab

    subroutine read_mic_subtr_backgr_shrink( micname, smpd, scale, pcontrast, mic_raw, mic_shrink, l_empty, mic_mask )
        class(string),                  intent(in)    :: micname !< micrograph file name
        real,                           intent(in)    :: smpd    !< sampling distance in A
        real,                           intent(in)    :: scale   !< scale factor
        character(len=*),               intent(in)    :: pcontrast
        class(image),                   intent(inout) :: mic_raw, mic_shrink
        logical,                        intent(out)   :: l_empty
        logical, allocatable, optional, intent(inout) :: mic_mask(:,:)
        integer :: nframes, ldim(3), ldim_shrink(3)
        ! set micrograph info
        call find_ldim_nptcls(micname, ldim, nframes)
        if( ldim(3) /= 1 .or. nframes /= 1 ) THROW_HARD('Only for 2D images')
        ! set shrunked dims
        ldim_shrink(1) = round2even(real(ldim(1)) * scale)
        ldim_shrink(2) = round2even(real(ldim(2)) * scale)
        ldim_shrink(3) = 1
        ! make shrunken micrograph
        call mic_shrink%new(ldim_shrink, smpd/scale)
        ! read micrograph
        call mic_raw%new(ldim, smpd)
        call mic_raw%read(micname)
        l_empty = mic_raw%is_empty()
        if( l_empty ) return
        call mic_raw%subtract_background(HP_BACKGR_SUBTR)
        call mic_raw%fft
        select case(trim(pcontrast))
            case('black')
                ! flip contrast (assuming black particle contrast on input)
                call mic_raw%mul(-1.)
            case('white')
                ! nothing to do
            case DEFAULT
                THROW_HARD('uknown pcontrast parameter, use (black|white)')
        end select
        ! attempt at detecting ice patches
        if( present(mic_mask) )then
            if( allocated(mic_mask) )then
                if( .not.all(shape(mic_mask)==ldim_shrink(1:2)) )then
                    deallocate(mic_mask)
                    allocate(mic_mask(ldim_shrink(1),ldim_shrink(2)))
                endif
            else
                allocate(mic_mask(ldim_shrink(1),ldim_shrink(2)))
            endif
            mic_mask = .true.
            call flag_ice(mic_raw, mic_mask)
        endif
        call mic_raw%mul(real(product(ldim))) ! to prevent numerical underflow when performing FFT
        ! shrink micrograph
        call mic_shrink%set_ft(.true.)
        call mic_raw%clip(mic_shrink)
        call mic_raw%ifft
        call mic_raw%div(real(product(ldim)))
        call mic_shrink%ifft
    end subroutine read_mic_subtr_backgr_shrink

    subroutine read_mic_subtr_backgr( micname, smpd, pcontrast, mic_raw, l_empty )
        class(string),    intent(in)    :: micname !< micrograph file name
        real,             intent(in)    :: smpd    !< sampling distance in A
        character(len=*), intent(in)    :: pcontrast
        class(image),     intent(inout) :: mic_raw
        logical,          intent(out)   :: l_empty
        integer :: nframes, ldim(3)
        ! set micrograph info
        call find_ldim_nptcls(micname, ldim, nframes)
        if( ldim(3) /= 1 .or. nframes /= 1 ) THROW_HARD('Only for 2D images')
        ! read micrograph
        call mic_raw%new(ldim, smpd)
        call mic_raw%read(micname)
        l_empty = mic_raw%is_empty()
        if( l_empty ) return
        call mic_raw%subtract_background(HP_BACKGR_SUBTR)
        call mic_raw%fft
        select case(trim(pcontrast))
            case('black')
                ! flip contrast (assuming black particle contrast on input)
                call mic_raw%mul(-1.)
            case('white')
                ! nothing to do
            case DEFAULT
                THROW_HARD('uknown pcontrast parameter, use (black|white)')
        end select
        call mic_raw%mul(real(product(ldim))) ! to prevent numerical underflow when performing FFT
        call mic_raw%ifft
        call mic_raw%div(real(product(ldim)))
    end subroutine read_mic_subtr_backgr

    subroutine read_mic( micname, mic_out )
        class(string), intent(in)    :: micname !< micrograph file name
        type(image),   intent(inout) :: mic_out
        integer :: nframes, ldim(3)
        real    :: smpd
        call find_ldim_nptcls(micname, ldim, nframes, smpd=smpd)
        if( ldim(3) /= 1 .or. nframes /= 1 ) THROW_HARD('Only for 2D images')
        call mic_out%new(ldim, smpd)
        call mic_out%read(micname)
    end subroutine read_mic

    subroutine cascade_filter_biomol( mic_shrink, mic4viz )
        class(image),           intent(inout) :: mic_shrink
        class(image), optional, intent(inout) :: mic4viz
        real,    parameter :: LP_UB = 15., LAM_ICM = 100., DAMP = 10., LAM_TV = 7.
        integer, parameter :: WINSZ_MED = 3
        call cascade_filter(mic_shrink, DAMP, LP_UB, LAM_TV, LAM_ICM, WINSZ_MED, mic4viz)
    end subroutine cascade_filter_biomol

    subroutine cascade_filter( mic_shrink, damp_below_zero, lp, lam_tv, lam_icm, winsz_med, mic4viz )
        use simple_tvfilter, only: tvfilter
        class(image),           intent(inout) :: mic_shrink
        real,                   intent(in)    :: damp_below_zero, lp, lam_tv, lam_icm
        integer,                intent(in)    :: winsz_med
        class(image), optional, intent(inout) :: mic4viz
        type(tvfilter) :: tvf
        call mic_shrink%zero_edgeavg
        ! dampens below zero
        call mic_shrink%div_below(0.,damp_below_zero)
        ! low-pass filter
        call mic_shrink%bp(0.,lp)
        ! TV denoising
        call tvf%new()
        call tvf%apply_filter(mic_shrink, lam_tv)
        call tvf%kill
        call mic4viz%copy(mic_shrink)
        ! Non-local-means denoising
        call mic_shrink%NLmean2D
        ! Iterated conditional modes denoising
        call mic_shrink%ICM2D(lam_icm)
        ! Median filter
        call mic_shrink%real_space_filter(winsz_med, 'median')
    end subroutine cascade_filter

    subroutine tv_filter_biomol( mic_shrink )
        class(image), intent(inout) :: mic_shrink
        real,    parameter :: LP_UB = 15., DAMP = 10., LAM_TV = 7.
        integer, parameter :: WINSZ_MED = 3
        call tv_filter(mic_shrink, DAMP, LP_UB, LAM_TV)
    end subroutine tv_filter_biomol

    subroutine tv_filter( mic_shrink, damp_below_zero, lp, lam_tv )
        use simple_tvfilter, only: tvfilter
        class(image), intent(inout) :: mic_shrink
        real,         intent(in)    :: damp_below_zero, lp, lam_tv
        type(tvfilter) :: tvf
        call mic_shrink%zero_edgeavg
        ! dampens below zero
        call mic_shrink%div_below(0.,damp_below_zero)
        ! low-pass filter
        call mic_shrink%bp(0.,lp) 
        ! TV denoising
        call tvf%new()
        call tvf%apply_filter(mic_shrink, lam_tv)
        call tvf%kill
    end subroutine tv_filter

    subroutine binarize_mic_den( mic_den, frac_fg, mic_bin )
        class(image),    intent(in)  :: mic_den
        real,            intent(in)  :: frac_fg
        class(image_bin), intent(out) :: mic_bin
        real :: bin_t
        call mic_bin%transfer2bimg(mic_den)
        call mic_bin%calc_bin_thres(frac_fg, bin_t)
        call mic_bin%binarize(bin_t)
        call mic_bin%set_imat            
        call mic_bin%erode() ! -4 A
        call mic_bin%erode() ! -4 A
        call mic_bin%set_largestcc2background
        call mic_bin%inv_bimg()
    end subroutine binarize_mic_den

    subroutine flag_amorphous_carbon( micrograph, picking_mask )
        use simple_histogram, only: histogram
        class(image),         intent(in)    :: micrograph
        logical, allocatable, intent(inout) :: picking_mask(:,:)
        integer, parameter :: K             = 5
        integer, parameter :: BOX           = 32 ! multiple of 4
        integer, parameter :: NBINS         = 64
        real,    parameter :: TVD_THRESHOLD = 0.2
        real,    parameter :: MIN_TVD_DIFF  = 0.05
        real,    parameter :: MIC_LP        = 15.
        logical, parameter :: DEBUG_HERE    = .false.
        type(histogram), allocatable :: hists(:,:)
        logical,         allocatable :: final_mask(:,:)
        type(ran_tabu)  :: rt
        type(histogram) :: khists(K)
        type(image)     :: mic, patches(nthr_glob)
        real    :: smpd
        integer :: dims(3), i
        logical :: found, empty
        found  = .false.
        dims   = micrograph%get_ldim()
        smpd   = micrograph%get_smpd()
        if( allocated(picking_mask) ) deallocate(picking_mask)
        allocate(picking_mask(dims(1),dims(2)),final_mask(dims(1),dims(2)),source=.true.)
        call mic%copy(micrograph)
        call clustering_rejection(found)
        if( found ) picking_mask = picking_mask .and. final_mask
        if( DEBUG_HERE ) print *,' % amorphous carbon: ',100.*real(count(.not.picking_mask))/real(product(dims))
        call mic%kill
        call rt%kill
        do i = 1,nthr_glob
            call patches(i)%kill
        enddo
        call khists(:)%kill
        call hists(:,:)%kill
        deallocate(hists)
        contains

            subroutine clustering_rejection( found )
                logical, intent(out) :: found
                integer, allocatable :: kmeans_labels(:), tmp_labels(:), inds(:)
                type(image)          :: tmpimg, tmpimg2
                type(histogram)      :: hist1, hist2
                real(dp) :: tmp
                real     :: bin_vars(K), tmpvec(K), tvd_inter(K-1), kmeans_score, mean, tvd, min_tvd
                integer  :: b,twob,bon2,bon4,i,j,m,n,ilb,jlb,iub,jub,nx,ny,nxy
                integer  :: bin,repeat,ii,jj,carbon,ithr
                logical  :: mask(dims(1),dims(2)), bin_mask(K), outside, cluster_mask(K-1)
                found = .false.
                ! gradients magnitude image
                call mic%fft
                call mic%bp(0.,MIC_LP)
                call mic%ifft
                call mic%norm
                call mic%gradients_magnitude(tmpimg)
                b     = BOX
                twob  = 2*b
                bon2  = b/2
                bon4  = b/4
                nx    = ceiling(real(dims(1))/real(b))
                ny    = ceiling(real(dims(2))/real(b))
                nxy   = nx*ny
                ! histograms
                allocate(kmeans_labels(nxy), tmp_labels(nxy),hists(nx,ny))
                !$omp parallel do collapse(2) proc_bind(close) private(i,j,ithr,ii,ilb,jlb,iub,jub)&
                !$omp default(shared)
                do i = 1,nx
                    do j =1,ny
                        ithr = omp_get_thread_num() + 1
                        ii  = (i-1)*nx+j
                        ilb = max(1, (i-1)*b-bon2+1)
                        iub = min(ilb+twob-1, dims(1))
                        ilb = iub-twob+1
                        jlb = max(1, (j-1)*b-bon2+1)
                        jub = min(jlb+twob-1, dims(2))
                        jlb = jub-twob+1
                        if(.not.patches(ithr)%exists() ) call patches(ithr)%new([twob,twob,1], smpd, wthreads=.false.)
                        call tmpimg%window_slim([ilb-1,jlb-1],twob, patches(ithr), outside)
                        call hists(i,j)%new(patches(ithr),  NBINS, minmax=[0.,4.])
                    enddo
                enddo
                !$omp end parallel do
                ! Cluster patches into K clusters
                kmeans_score = huge(0.)
                do i = 1,K
                    call khists(i)%new(hists(1,1))
                enddo
                call seed_rnd
                do repeat = 1,300
                    rt = ran_tabu(nxy)
                    call cluster(nxy, nx, ny, K, tmp, tmp_labels)
                    if( real(tmp) < kmeans_score )then
                        kmeans_score  = real(tmp)
                        kmeans_labels = tmp_labels
                    endif
                enddo
                do bin = 1,K
                    bin_mask(bin) = count(kmeans_labels==bin) > 0
                enddo
                if( count(bin_mask) == 1 ) return ! clustering fail, do nothing
                ! histograms for all clusters
                do i = 1,K
                    call khists(i)%zero
                enddo
                bin = 0
                do i = 1,nx
                    do j =1,ny
                        bin = bin+1
                        call khists(kmeans_labels(bin))%add(hists(i,j))
                    enddo
                enddo
                ! deviation from mode
                bin_vars = -1.
                do bin = 1,K
                    if( bin_mask(bin) )then
                        mean          = khists(bin)%hmode()
                        bin_vars(bin) = khists(bin)%variance(mean=mean)
                    endif
                enddo
                ! sorting by cluster variance
                inds   = (/(i,i=1,K)/)
                tmpvec = bin_vars
                call hpsort(tmpvec, inds)
                ! debug
                if( DEBUG_HERE )then
                    do i = 1,K
                        call tmpimg2%copy(tmpimg)
                        j = 0
                        do ii = 1,nx
                            ilb = max(1, (ii-1)*b+1)
                            iub = min(ilb+b-1, dims(1))
                            ilb = iub-b+1
                            do jj =1,ny
                                j = j + 1
                                if( kmeans_labels(j) == inds(i) ) cycle
                                jlb = max(1, (jj-1)*b+1)
                                jub = min(jlb+b-1, dims(2))
                                jlb = jub-b+1
                                do m = ilb,iub
                                    do n = jlb,jub
                                        call tmpimg2%set([m,n,1],0.)
                                    enddo
                                enddo
                            enddo
                        enddo
                        call tmpimg2%write(string('clusters.mrc'),i)
                    enddo
                    call mic%write(string('clusters.mrc'),K+1)
                    call tmpimg2%kill
                endif
                call tmpimg%kill
                ! Segmentation by agglomerative inter-cluster distance
                tvd_inter = -.1
                do i = 1,K-1
                    empty = .true.
                    call hist1%new(khists(1))
                    do j = 1,i
                        jj = inds(j)
                        if( bin_mask(jj) )then
                            call hist1%add(khists(jj))
                            empty = .false.
                        endif
                    enddo
                    if( empty ) cycle
                    empty = .true.
                    call hist2%new(khists(1))
                    do j = i+1,K
                        jj = inds(j)
                        if( bin_mask(jj) )then
                            call hist2%add(khists(jj))
                            empty= .false.
                        endif
                    enddo
                    if( empty ) cycle
                    tvd_inter(i) = hist1%tvd(hist2)
                    if( DEBUG_HERE )then
                        write(*,'(I6,5F9.4)') i, hist1%variance(), hist2%variance(), hist1%tvd(hist2),&
                            &khists(inds(i))%variance(), khists(inds(i))%mean()
                    endif
                enddo
                cluster_mask = tvd_inter > 0.
                carbon       = maxloc(tvd_inter,dim=1)
                tvd          = tvd_inter(carbon)
                min_tvd      = minval(tvd_inter,mask=cluster_mask)
                bin          = findloc(cluster_mask(1:K-1),.true.,dim=1) ! first non-empty cluster
                ! Decision
                found = .false.
                if( bin == 0 )then
                    tvd = 0.
                else
                    if( carbon == bin ) tvd = 0.
                    if( abs(tvd-min_tvd) < MIN_TVD_DIFF ) tvd = 0.
                endif
                found = tvd > TVD_THRESHOLD
                ! Dejection
                if( found )then
                    mask = final_mask
                    do bin = carbon+1,K
                        n = inds(bin)
                        if( .not.bin_mask(n) ) cycle
                        m     = 0
                        do i = 1,nx
                            ii = min((i-1)*b+1, dims(1)-b+1)
                            do j =1,ny
                                jj = min((j-1)*b+1, dims(2)-b+1)
                                m = m+1
                                if( kmeans_labels(m) == n ) mask(ii:ii+b-1,jj:jj+b-1) = .false.
                            enddo
                        enddo
                    enddo
                    ! erosion with half window size
                    final_mask = mask
                    do i = 1,2*nx
                        ilb = max(1, (i-1)*bon2-bon4+1)
                        iub = min(ilb+b-1, dims(1))
                        ilb = iub - b + 1 + bon4
                        iub = iub - bon4
                        do j = 1,2*ny
                            jlb = max(1, (j-1)*bon2-bon4+1)
                            jub = min(jlb+b-1, dims(2))
                            jlb = jub - b + 1 + bon4
                            jub = jub - bon4
                            n = count(mask(ilb-bon4:iub+bon4,jlb-bon4:jub+bon4))
                            if( n == 0   )cycle
                            if( n == b*b )cycle
                            final_mask(ilb:iub,jlb:jub) = .true.
                        enddo
                    enddo
                    ! dilation with half window size x 2
                    mask = final_mask
                    do i = 1,2*nx
                        ilb = max(1, (i-1)*bon2-bon4+1)
                        iub = min(ilb+b-1, dims(1))
                        ilb = iub - b + 1 + bon4
                        iub = iub - bon4
                        do j = 1,2*ny
                            jlb = max(1, (j-1)*bon2-bon4+1)
                            jub = min(jlb+b-1, dims(2))
                            jlb = jub - b + 1 + bon4
                            jub = jub - bon4
                            if( all(mask(ilb:iub,jlb:jub)) )then
                                n = count(.not.mask(ilb-bon4:iub+bon4,jlb-bon4:jub+bon4))
                                if( n > 0 )then
                                    final_mask(ilb:iub,jlb:jub) = .false.
                                endif
                            endif
                        enddo
                    enddo
                    mask = final_mask
                    do i = 1,2*nx
                        ilb = max(1, (i-1)*bon2-bon4+1)
                        iub = min(ilb+b-1, dims(1))
                        ilb = iub - b + 1 + bon4
                        iub = iub - bon4
                        do j = 1,2*ny
                            jlb = max(1, (j-1)*bon2-bon4+1)
                            jub = min(jlb+b-1, dims(2))
                            jlb = jub - b + 1 + bon4
                            jub = jub - bon4
                            if( all(mask(ilb:iub,jlb:jub)) )then
                                n = count(.not.mask(ilb-bon4:iub+bon4,jlb-bon4:jub+bon4))
                                if( n > 0 )then
                                    final_mask(ilb:iub,jlb:jub) = .false.
                                endif
                            endif
                        enddo
                    enddo
                endif
                ! cleanup
                call hist1%kill
                call hist2%kill
            end subroutine clustering_rejection

            subroutine cluster( n, nx, ny, K, overall_score, labels )
                integer, parameter :: nits = 100
                integer,         intent(in)  :: n, nx, ny, K
                real(dp),        intent(out) :: overall_score
                integer,         intent(out) :: labels(n)
                integer  :: pops(K), new_labels(n), i,j,c,cl,it,nk
                real(dp) :: scores(K), kscores(K), prev_score, score
                ! initial labels
                c = 0
                do i = 1,n
                    c = c+1
                    if( c > K ) c = 1
                    labels(i) = c
                enddo
                call rt%shuffle(labels)
                ! main loop
                prev_score = huge(0.0)
                do it = 1,nits
                    scores = 0.d0
                    !$omp parallel private(i,j,c,cl,kscores) default(shared) proc_bind(close)
                    !$omp do
                    ! centers
                    do cl = 1,K
                        call khists(cl)%zero
                        pops(cl) = count(labels==cl)
                        if( pops(cl) == 0 ) cycle
                        do i = 1,nx
                            do j = 1,ny
                                c = (i-1)*ny+j
                                if( labels(c) /= cl )cycle
                                call khists(cl)%add(hists(i,j))
                            enddo
                        enddo
                    enddo
                    !$omp end do
                    !$omp do collapse(2) reduction(+:scores)
                    do i = 1,nx
                        do j = 1,ny
                            c = (i-1)*ny+j
                            do cl = 1,K
                                if( pops(cl) > 0 )then
                                    kscores(cl) = real(khists(cl)%tvd(hists(i,j)),dp)
                                else
                                    kscores(cl) = huge(0.)
                                endif
                            enddo
                            ! current parttioning
                            cl = labels(c)
                            scores(cl) = scores(cl) + kscores(cl)
                            ! new labels
                            new_labels(c) = minloc(kscores,dim=1)
                        enddo
                    enddo
                    !$omp end do
                    !$omp end parallel
                    nk     = count(pops>0)
                    score  = sum(scores/real(pops,dp),mask=pops>0) / real(nk,dp)
                    ! convergence
                    if( score < prev_score )then
                        if( abs(score - prev_score) < 1.d-6 ) exit
                    endif
                    prev_score = score
                    labels     = new_labels
                enddo
                overall_score = score
            end subroutine cluster

    end subroutine flag_amorphous_carbon

    subroutine flag_ice( micrograph, mask )
        class(image), intent(inout) :: micrograph   ! raw micrograph
        logical,      intent(inout) :: mask(:,:)    ! picking mask at the size the first pass of picking is performed
        real,    parameter   :: THRESHOLD = 5.
        real,    parameter   :: SMPD_ICE  = ICE_BAND1/2. - 0.15
        integer, parameter   :: BOX       = 128
        type(image)          :: boximgs_heap(nthr_glob), img
        real,    allocatable :: scores(:,:)
        integer, allocatable :: counts(:,:)
        real        :: score, scale, radius, smpd_raw, smpd_shrink
        integer     :: ldim_raw(3), ldim_shrink(3), ldim(3),i,j,ithr,ii,jj
        logical     :: outside
        smpd_raw = micrograph%get_smpd()
        if( smpd_raw > SMPD_ICE ) return
        ldim_raw         = micrograph%get_ldim()
        ldim_shrink(1:2) = shape(mask)
        ldim_shrink(3)   = 1
        smpd_shrink      = real(ldim_raw(1)) / real(ldim_shrink(1))
        scale            = smpd_shrink / SMPD_ICE
        ldim(1:2)        = round2even(real(ldim_shrink(1:2)) * scale)
        ldim(3)          = 1
        radius  = real(BOX)/2.- COSMSKHALFWIDTH
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new([BOX,BOX,1], SMPD_ICE, wthreads=.false.)
        end do
        call img%new(ldim,SMPD_ICE)
        call img%set_ft(.true.)
        if( micrograph%is_ft() )then
            call micrograph%clip(img)
        else
            call micrograph%fft
            call micrograph%clip(img)
            call micrograph%ifft
        endif
        call img%ifft
        allocate(scores(ldim(1),ldim(2)),source=0.)
        allocate(counts(ldim(1),ldim(2)),source=0)
        !$omp parallel do collapse(2) schedule(static) default(shared) private(i,j,ithr,outside,score) proc_bind(close)
        do i = 1,ldim(1)-BOX+1,BOX/2
            do j = 1,ldim(2)-BOX+1,BOX/2
                ithr = omp_get_thread_num() + 1
                call img%window_slim([i-1, j-1], BOX, boximgs_heap(ithr), outside)
                call boximgs_heap(ithr)%mask(radius,'soft')
                call boximgs_heap(ithr)%fft
                call boximgs_heap(ithr)%calc_ice_score(score)
                !$omp critical
                scores(i:i+BOX-1,j:j+BOX-1) = scores(i:i+BOX-1,j:j+BOX-1) + score
                counts(i:i+BOX-1,j:j+BOX-1) = counts(i:i+BOX-1,j:j+BOX-1) + 1
                !$omp end critical
            end do
        end do
        !$omp end parallel do
        i = ldim(1)-BOX+1
        !$omp parallel do schedule(static) default(shared) private(j,ithr,outside,score) proc_bind(close)
        do j = 1,ldim(2)-BOX+1,BOX/2
            ithr = omp_get_thread_num() + 1
            call img%window_slim([i-1, j-1], BOX, boximgs_heap(ithr), outside)
            call boximgs_heap(ithr)%mask(radius,'soft')
            call boximgs_heap(ithr)%fft
            call boximgs_heap(ithr)%calc_ice_score(score)
            !$omp critical
            scores(i:i+BOX-1,j:j+BOX-1) = scores(i:i+BOX-1,j:j+BOX-1) + score
            counts(i:i+BOX-1,j:j+BOX-1) = counts(i:i+BOX-1,j:j+BOX-1) + 1
            !$omp end critical
        enddo
        !$omp end parallel do
        j = ldim(2)-BOX+1
        !$omp parallel do schedule(static) default(shared) private(i,ithr,outside,score) proc_bind(close)
        do i = 1,ldim(1)-BOX+1,BOX/2
            ithr = omp_get_thread_num() + 1
            call img%window_slim([i-1, j-1], BOX, boximgs_heap(ithr), outside)
            call boximgs_heap(ithr)%mask(radius,'soft')
            call boximgs_heap(ithr)%fft
            call boximgs_heap(ithr)%calc_ice_score(score)
            !$omp critical
            scores(i:i+BOX-1,j:j+BOX-1) = scores(i:i+BOX-1,j:j+BOX-1) + score
            counts(i:i+BOX-1,j:j+BOX-1) = counts(i:i+BOX-1,j:j+BOX-1) + 1
            !$omp end critical
        enddo
        !$omp end parallel do
        where( counts > 0 ) scores = scores / real(counts)
        scale = real(ldim(1)) / real(ldim_shrink(1))
        do i = 1,ldim_shrink(1)
            ii = min(ldim(1), max(1, nint(scale*real(i))))
            do j = 1,ldim_shrink(2)
                jj = min(ldim(2), max(1, nint(scale*real(j))))
                mask(i,j) = mask(i,j) .and. (scores(ii,jj) < THRESHOLD)
            enddo
        enddo
        if( .not.all(mask) )then
            print *,' % ice: ', 100.-100.*count(mask)/real(product(ldim_shrink))
        endif
        ! cleanup
        call img%kill
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%kill
        end do
        deallocate(scores,counts)
    end subroutine flag_ice

end module simple_micproc
