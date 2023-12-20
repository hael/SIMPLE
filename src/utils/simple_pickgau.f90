module simple_pickgau
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
implicit none

public :: read_mic_raw, pickgau
private
#include "simple_local_flags.inc"

! class constants
real,    parameter :: GAUSIG   = 5., BOX_EXP_FAC = 0.111, NDEV_DEFAULT = 2.5
integer, parameter :: OFFSET   = 3 , OFFSET_UB   = 2 * OFFSET, MAXNREFS = 100
logical, parameter :: L_WRITE  = .false.

! class variables
integer     :: ldim_raw(3), ldim_raw_box(3)
real        :: smpd_raw
type(image) :: mic_raw

! instance
type pickgau
    private
    real                     :: smpd_shrink = 0., maxdiam = 0., sig = 0., sxx = 0., t = 0., ndev = 0.
    integer                  :: ldim(3), ldim_box(3), nboxes = 0, nboxes_ub = 0, nx = 0, ny = 0
    integer                  :: nx_offset = 0, ny_offset = 0, npeaks = 0, nrefs = 0
    type(image)              :: mic_shrink, gauref
    type(image), allocatable :: boximgs(:), boxrefs(:)
    logical,     allocatable :: l_mic_mask(:,:), l_err_refs(:)
    integer,     allocatable :: positions(:,:), inds_offset(:,:)
    real,        allocatable :: box_scores(:,:)
    logical                  :: refpick = .false.
    logical                  :: exists  = .false.
contains
    procedure          :: new_gaupicker
    procedure          :: new_refpicker
    procedure, private :: new
    procedure, private :: set_refs
    procedure, private :: setup_iterators
    ! procedure          :: gaupick
    ! procedure          :: refpick
    procedure, private :: extract_boximgs
    procedure, private :: destruct_boximgs
    procedure, private :: gaumatch_boximgs
    procedure, private :: refmatch_boximgs
    procedure, private :: detect_peaks
    procedure, private :: center_filter
    procedure, private :: distance_filter
    procedure, private :: remove_outliers
    procedure          :: get_positions
    procedure          :: refine_upscaled
    procedure          :: kill
end type

contains

    subroutine read_mic_raw( micname, smpd )
        character(len=*), intent(in) :: micname !< micrograph file name
        real,             intent(in) :: smpd    !< sampling distance in A
        integer :: nframes
        call find_ldim_nptcls(micname, ldim_raw, nframes)
        if( ldim_raw(3) /= 1 .or. nframes /= 1 ) THROW_HARD('Only for 2D images')
        smpd_raw = smpd
        call mic_raw%new(ldim_raw, smpd_raw)
        call mic_raw%read(micname)
    end subroutine read_mic_raw

    subroutine new_gaupicker( self, pcontrast, smpd_shrink, moldiam, moldiam_max, ndev )
        class(pickgau),   intent(inout) :: self
        character(len=*), intent(in)    :: pcontrast
        real,             intent(in)    :: smpd_shrink, moldiam
        real, optional,   intent(in)    :: moldiam_max
        real, optional,   intent(in)    :: ndev    !< # std devs for outlier detection
        call self%new( pcontrast, smpd_shrink, moldiam, moldiam_max, ndev )
    end subroutine new_gaupicker

    subroutine new_refpicker( self, pcontrast, smpd_shrink, mskdiam, imgs, ndev )
        class(pickgau),   intent(inout) :: self
        character(len=*), intent(in)    :: pcontrast
        real,             intent(in)    :: smpd_shrink, mskdiam
        class(image),     intent(inout) :: imgs(:)
        real, optional,   intent(in)    :: ndev    !< # std devs for outlier detection
        call self%new( pcontrast, smpd_shrink, mskdiam, mskdiam, mskdiam, ndev )
        call self%set_refs( imgs, mskdiam )
        call self%setup_iterators
    end subroutine new_refpicker

    subroutine new( self, pcontrast, smpd_shrink, moldiam, moldiam_max, mskdiam, ndev )
        class(pickgau),   intent(inout) :: self
        character(len=*), intent(in)    :: pcontrast
        real,             intent(in)    :: smpd_shrink, moldiam, moldiam_max
        real, optional,   intent(in)    :: mskdiam !< reference-based picking if present
        real, optional,   intent(in)    :: ndev    !< # std devs for outlier detection
        character(len=:), allocatable   :: numstr
        integer     :: ldim_pd(3), box_max
        type(image) :: mic_pad, gauimg
        real        :: hpfreq, scale, maxdiam_max
        if( self%exists ) call self%kill
        self%smpd_shrink = smpd_shrink
        if( trim(params_glob%pick_roi).eq.'yes' )then
            ! high-pass micrograph
            ldim_pd(1:2) = find_larger_magic_box(ldim_raw(1:2)+1)
            if( minval(ldim_pd(1:2) - ldim_raw(1:2)) < 16 )then
                ldim_pd(1:2) = find_larger_magic_box(ldim_pd(1:2)+1)
            endif
            ldim_pd(3) = 1
            call mic_pad%new(ldim_pd, smpd_raw)
            call mic_raw%pad_mirr(mic_pad)
            call mic_pad%fft
            hpfreq = real(minval(ldim_pd(1:2)))*self%smpd_shrink/16.
            call mic_pad%bp(hpfreq,0.)
            call mic_pad%ifft
            call mic_pad%clip(mic_raw)
            call mic_pad%kill
        endif
        ! set shrunken logical dimensions
        scale            = smpd_raw / self%smpd_shrink
        self%ldim(1)     = round2even(real(ldim_raw(1)) * scale)
        self%ldim(2)     = round2even(real(ldim_raw(2)) * scale)
        self%ldim(3)     = 1
        ! set logical dimensions of boxes
        if( present(mskdiam) )then ! reference-based picking
            self%maxdiam = mskdiam
        else
            self%maxdiam = moldiam + moldiam * BOX_EXP_FAC
        endif
        ldim_raw_box(1)  = round2even(self%maxdiam / smpd_raw)
        ldim_raw_box(2)  = ldim_raw_box(1)
        ldim_raw_box(3)  = 1
        self%ldim_box(1) = round2even(self%maxdiam / self%smpd_shrink)
        self%ldim_box(2) = self%ldim_box(1)
        self%ldim_box(3) = 1
        ! set # pixels in x/y based on maximum expected box size
        if( present(mskdiam) )then ! reference-based picking
            maxdiam_max  = mskdiam
        else
            maxdiam_max  = moldiam_max + moldiam_max * BOX_EXP_FAC
        endif
        box_max          = round2even(maxdiam_max / self%smpd_shrink)
        self%nx          = self%ldim(1) - box_max
        self%ny          = self%ldim(2) - box_max
        ! set Gaussian
        self%sig = ((self%maxdiam / 2.) / self%smpd_shrink) / GAUSIG
        call self%gauref%new(self%ldim_box, self%smpd_shrink)
        call self%gauref%gauimg2D(self%sig, self%sig)
        call self%gauref%prenorm4real_corr(self%sxx)
        ! set ndev
        self%ndev = NDEV_DEFAULT
        if( present(ndev) ) self%ndev = ndev
        ! shrink micrograph
        call self%mic_shrink%new(self%ldim, self%smpd_shrink)
        call self%mic_shrink%set_ft(.true.)
        call mic_raw%fft
        call mic_raw%clip(self%mic_shrink)
        call mic_raw%clip(self%mic_shrink)
        select case(trim(pcontrast))
            case('black')
                ! flip contrast (assuming black particle contrast on input)
                call self%mic_shrink%mul(-1.)
            case('white')
                ! nothing to do
            case DEFAULT
                THROW_HARD('uknown pcontrast parameter, use (black|white)')
        end select
        ! back to real-space
        call mic_raw%ifft
        call self%mic_shrink%ifft
        if( L_WRITE )then
            if( present(mskdiam) )then ! reference-based picking
                numstr = int2str(nint(mskdiam))
                call self%mic_shrink%write('mic_shrink_mskdiam'//numstr//'.mrc')
                call self%gauref%write('gauref_mskdiam'//numstr//'.mrc')
            else
                numstr = int2str(nint(moldiam))
                call self%mic_shrink%write('mic_shrink_moldiam'//numstr//'.mrc')
                call self%gauref%write('gauref_moldiam'//numstr//'.mrc')
            endif
        endif
        ! denoise mic_shrink
        call gauimg%new(self%ldim, self%smpd_shrink)
        call gauimg%gauimg2D(self%sig,self%sig)
        call gauimg%fft
        call self%mic_shrink%fft
        call self%mic_shrink%mul(gauimg)
        call self%mic_shrink%ifft
        call gauimg%kill
        if( L_WRITE )then
            if( present(mskdiam) )then ! reference-based picking
                call self%mic_shrink%write('gauconv_mic_shrink_mskdiam'//numstr//'.mrc')
            else
                call self%mic_shrink%write('gauconv_mic_shrink_moldiam'//numstr//'.mrc')
            endif
        endif
        if( present(mskdiam) )then ! reference-based picking
            ! iterators are established after calling set_refs because ldim_box, nx & ny are updated
        else
            ! establish iterators
            call self%setup_iterators
        endif
        ! flag existence
        self%exists = .true.
    end subroutine new

    subroutine set_refs( self, imgs, mskdiam )
        class(pickgau), intent(inout) :: self
        class(image),   intent(inout) :: imgs(:)
        real,           intent(in)    :: mskdiam
        character(len=:), allocatable :: numstr
        type(image) :: img_rot
        integer     :: ldim(3), iimg, nimgs, irot, nrots, cnt
        real        :: scale, mskrad, pixrad_shrink1, pixrad_shrink2, smpd, ang, rot
        smpd       = imgs(1)%get_smpd()
        if( abs(smpd - smpd_raw) > 0.01 ) THROW_HARD('raw micrograph and picking references of different scale')
        ldim       = imgs(1)%get_ldim()
        if( ldim(3) /= 1 ) THROW_HARD('box references must be 2D')
        nimgs      = size(imgs)
        nrots      = nint(real(MAXNREFS) / real(nimgs))
        self%nrefs = nimgs * nrots
        ! set shrunken logical dimensions of boxes
        scale            = smpd / self%smpd_shrink
        self%ldim_box(1) = round2even(real(ldim(1)) * scale)
        self%ldim_box(2) = round2even(real(ldim(2)) * scale)
        self%ldim_box(3) = 1
        ! set # pixels in x/y for both box sizes
        self%nx          = self%ldim(1) - self%ldim_box(1)
        self%ny          = self%ldim(2) - self%ldim_box(2)
        mskrad           = (mskdiam / self%smpd_shrink) / 2.
        ! set Gaussian
        call self%gauref%new(self%ldim_box, self%smpd_shrink)
        call self%gauref%gauimg2D(self%sig, self%sig)
        call self%gauref%prenorm4real_corr(self%sxx)
        if( L_WRITE )then
            numstr = int2str(nint(mskdiam))
            call self%gauref%write('gauref_mskdiam'//numstr//'.mrc')
        endif
        if( allocated(self%boxrefs) )then
            do iimg = 1,size(self%boxrefs)
                call self%boxrefs(iimg)%kill
            end do
            deallocate(self%boxrefs)
        endif
        if( allocated(self%l_err_refs) ) deallocate(self%l_err_refs)
        allocate(self%l_err_refs(self%nrefs))
        call img_rot%new(ldim, smpd)
        ang = 360./real(nrots)
        cnt = 0
        do iimg = 1,nimgs
            rot = 0.
            do irot = 1,nrots
                cnt = cnt + 1
                call imgs(iimg)%rtsq(rot, 0., 0., img_rot)
                call img_rot%fft
                call self%boxrefs(cnt)%new(self%ldim_box, self%smpd_shrink, wthreads=.false.)
                call self%boxrefs(cnt)%set_ft(.true.)
                call img_rot%clip(self%boxrefs(cnt))
                rot = rot + ang
            end do
        end do
        !$omp parallel do schedule(static) default(shared) private(cnt) proc_bind(close)
        do cnt = 1,self%nrefs
            ! convolve with Gaussian
            call self%boxrefs(cnt)%mul(self%gauref)
            ! back to real-space
            call self%boxrefs(cnt)%ifft
            call self%boxrefs(cnt)%mask(mskrad, 'hard')
            call self%boxrefs(cnt)%prenorm4real_corr(self%l_err_refs(cnt))
        end do
        !$omp end parallel do
        if( L_WRITE )then
            do cnt = 1,self%nrefs
                call self%boxrefs(cnt)%write('boxrefs_shrink.mrc', cnt)
            enddo
        endif
        call self%gauref%ifft
        call img_rot%kill
        self%refpick = .true.
    end subroutine set_refs

    subroutine setup_iterators( self )
        class(pickgau), intent(inout) :: self
        integer :: xind, yind
        ! allocate and set l_mic_mask
        if( allocated(self%l_mic_mask) ) deallocate(self%l_mic_mask)
        allocate(self%l_mic_mask(self%ldim(1),self%ldim(2)), source=.true.)
        self%nboxes    = 0
        self%nx_offset = 0
        do xind = 0,self%nx,OFFSET
            self%nx_offset = self%nx_offset + 1
            self%ny_offset = 0
            do yind = 0,self%ny,OFFSET
                self%nboxes    = self%nboxes    + 1
                self%ny_offset = self%ny_offset + 1
                if( self%l_mic_mask(xind+1,yind+1) ) self%nboxes = self%nboxes + 1
            end do
        end do
        ! count # boxes, upper bound
        self%nboxes_ub = 0
        do xind = 0,self%nx,OFFSET_UB
            do yind = 0,self%ny,OFFSET_UB
                if( self%l_mic_mask(xind+1,yind+1) ) self%nboxes_ub = self%nboxes_ub + 1
            end do
        end do
        ! allocate and set positions and inds_offset
        if( allocated(self%positions)   ) deallocate(self%positions)
        if( allocated(self%inds_offset) ) deallocate(self%inds_offset)
        allocate(self%positions(self%nboxes,2), self%inds_offset(self%nx_offset,self%ny_offset), source=0)
        ! calculate total # boxes
        self%nboxes    = 0
        self%nx_offset = 0
        do xind = 0,self%nx,OFFSET
            self%nx_offset = self%nx_offset + 1
            self%ny_offset = 0
            do yind = 0,self%ny,OFFSET
                self%ny_offset = self%ny_offset + 1
                if( self%l_mic_mask(xind+1,yind+1) )then
                    self%nboxes = self%nboxes + 1
                    self%positions(self%nboxes,:) = [xind,yind]
                    self%inds_offset(self%nx_offset,self%ny_offset) = self%nboxes
                endif
            end do
        end do
        ! allocate box_scores
        if( allocated(self%box_scores) ) deallocate(self%box_scores)
        allocate(self%box_scores(self%nx_offset,self%ny_offset), source = -1.)
    end subroutine setup_iterators

    subroutine extract_boximgs( self )
        class(pickgau), intent(inout) :: self
        integer :: pos(2), ibox
        logical :: outside
        if( allocated(self%boximgs) )then
            do ibox = 1,self%nboxes
                call self%boximgs(ibox)%kill
            end do
            deallocate(self%boximgs)
        endif
        allocate(self%boximgs(self%nboxes))
        !$omp parallel do schedule(static) default(shared) private(ibox,outside,pos) proc_bind(close)
        do ibox = 1,self%nboxes
            call self%boximgs(ibox)%new(self%ldim_box, self%smpd_shrink)
            pos = self%positions(ibox,:)
            call self%mic_shrink%window_slim(pos, self%ldim_box(1), self%boximgs(ibox), outside)
        end do
        !$omp end parallel do
    end subroutine extract_boximgs

    subroutine destruct_boximgs( self )
        class(pickgau), intent(inout) :: self
        integer :: ibox
        if( allocated(self%boximgs) )then
            do ibox = 1,self%nboxes
                call self%boximgs(ibox)%kill
            end do
            deallocate(self%boximgs)
        endif
    end subroutine destruct_boximgs

    subroutine gaumatch_boximgs( self )
        class(pickgau), intent(inout) :: self
        logical     :: outside
        integer     :: pos(2), ioff, joff, ithr
        type(image) :: boximgs_heap(nthr_glob)
        ! construct heap
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box, self%smpd_shrink)
        end do
        ! calculate box_scores
        !$omp parallel do schedule(static) collapse(2) default(shared) private(ioff,joff,ithr,outside,pos) proc_bind(close)
        do ioff = 1,self%nx_offset
            do joff = 1,self%ny_offset
                if( self%inds_offset(ioff,joff) == 0 )then
                    self%box_scores(ioff,joff) = -1.
                else
                    pos = self%positions(self%inds_offset(ioff,joff),:)
                    if( self%l_mic_mask(pos(1)+1,pos(2)+1) )then
                        ithr = omp_get_thread_num() + 1
                        call self%mic_shrink%window_slim(pos, self%ldim_box(1), boximgs_heap(ithr), outside)
                        self%box_scores(ioff,joff) = self%gauref%real_corr_prenorm(boximgs_heap(ithr), self%sxx)
                    else
                        self%box_scores(ioff,joff) = -1.
                    endif
                endif
            end do
        end do
        !$omp end parallel do
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%kill
        end do
    end subroutine gaumatch_boximgs

    subroutine refmatch_boximgs( self )
        class(pickgau), intent(inout) :: self
        logical     :: outside, l_err
        integer     :: pos(2), ioff, joff, ithr, iref
        real        :: scores(self%nrefs)
        type(image) :: boximgs_heap(nthr_glob)
        ! construct heap
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box, self%smpd_shrink)
        end do
        ! calculate box_scores
        if( self%refpick )then
            !$omp parallel do schedule(static) collapse(2) default(shared) private(ioff,joff,ithr,outside,iref,scores,pos,l_err) proc_bind(close)
            do ioff = 1,self%nx_offset
                do joff = 1,self%ny_offset
                    if( self%inds_offset(ioff,joff) == 0 )then
                        self%box_scores(ioff,joff) = -1.
                    else
                        pos = self%positions(self%inds_offset(ioff,joff),:)
                        if( self%l_mic_mask(pos(1)+1,pos(2)+1) )then
                            ithr = omp_get_thread_num() + 1
                            pos  = self%positions(self%inds_offset(ioff,joff),:)
                            call self%mic_shrink%window_slim(pos, self%ldim_box(1), boximgs_heap(ithr), outside)
                            call boximgs_heap(ithr)%prenorm4real_corr(l_err)
                            if( l_err )then
                                self%box_scores(ioff,joff) = 0.
                            else
                                do iref = 1,self%nrefs
                                    if( self%l_err_refs(iref) )then
                                        scores(iref) = 0.
                                    else
                                        scores(iref) = self%boxrefs(iref)%real_corr_prenorm(boximgs_heap(ithr))
                                    endif
                                end do
                                self%box_scores(ioff,joff) = maxval(scores)
                            endif
                        else
                            self%box_scores(ioff,joff) = -1.
                        endif
                    endif
                end do
            end do
            !$omp end parallel do
        else
            THROW_HARD('Instance not setup for reference-based picking')
        endif
    end subroutine refmatch_boximgs

    subroutine detect_peaks( self )
        class(pickgau), intent(inout) :: self
        real, allocatable :: tmp(:)
        integer :: ioff, joff
        tmp = pack(self%box_scores, mask=(self%box_scores > -1. + TINY))
        call detect_peak_thres(size(tmp), self%nboxes_ub, tmp, self%t)
        deallocate(tmp)
        self%t = max(0.,self%t)
        where( self%box_scores >= self%t )
            ! there's a peak
        elsewhere
            self%box_scores = -1.
        end where
        self%npeaks = count(self%box_scores >= self%t)
        write(logfhandle,'(a,1x,I5)') '# peaks detected:                      ', self%npeaks
    end subroutine detect_peaks

     subroutine center_filter( self )
        class(pickgau), intent(inout) :: self
        type(image) :: boximgs_heap(nthr_glob)
        integer     :: ioff, joff, ithr, npeaks
        real        :: scores_cen(self%nx_offset,self%ny_offset)
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box, self%smpd_shrink)
        end do
        !$omp parallel do schedule(static) collapse(2) default(shared) private(ioff,joff,ithr) proc_bind(close)
        do ioff = 1,self%nx_offset
            do joff = 1,self%ny_offset
                ithr = omp_get_thread_num() + 1
                if( self%box_scores(ioff,joff) >= self%t )then
                    scores_cen(ioff,joff) = self%boximgs(self%inds_offset(ioff,joff))%box_cen_arg(boximgs_heap(ithr))
                else
                    scores_cen(ioff,joff) = real(OFFSET) + 1.
                endif
            end do
        end do
        !$omp end parallel do
        npeaks = count(scores_cen <= real(OFFSET))
        write(logfhandle,'(a,1x,I5)') '# positions before   center filtering: ', count(self%box_scores >= self%t)
        write(logfhandle,'(a,1x,I5)') '# positions after    center filtering: ', npeaks
        ! modify box_scores and update npeaks
        where( scores_cen <= real(OFFSET) )
            ! there's a peak
        elsewhere
            self%box_scores = -1.
        endwhere
        self%npeaks = count(self%box_scores >= self%t)
        write(logfhandle,'(a,1x,I5)') '# npeaks    after    center filtering: ', self%npeaks 
    end subroutine center_filter

    subroutine distance_filter( self, dist_thres )
        class(pickgau), intent(inout) :: self
        real,           intent(in)    :: dist_thres
        integer, allocatable :: pos_inds(:)
        real,    allocatable :: pos_scores(:)
        logical, allocatable :: mask(:), selected_pos(:)
        integer :: nbox, npeaks, ibox, jbox, loc, ioff, joff, ithr, ipeak
        real    :: dist
        logical :: is_peak
        pos_inds   = pack(self%inds_offset(:,:), mask=self%box_scores(:,:) >= self%t)
        pos_scores = pack(self%box_scores(:,:),  mask=self%box_scores(:,:) >= self%t)
        nbox       = size(pos_inds)
        allocate(mask(nbox),         source=.false.)
        allocate(selected_pos(nbox), source=.true. )
        do ibox = 1,nbox
            mask = .false.
            !$omp parallel do schedule(static) default(shared) private(jbox,dist) proc_bind(close)
            do jbox = 1,nbox
                dist = euclid(real(self%positions(pos_inds(ibox),:)),real(self%positions(pos_inds(jbox),:)))
                if( dist <= dist_thres ) mask(jbox) = .true.
            end do
            !$omp end parallel do
            ! find best match in the neigh
            loc = maxloc(pos_scores, mask=mask, dim=1)
            ! eliminate all but the best
            mask(loc) = .false.
            where( mask ) selected_pos = .false.
        end do
        npeaks = count(selected_pos)
        write(logfhandle,'(a,1x,I5)') '# positions before distance filtering: ', nbox
        write(logfhandle,'(a,1x,I5)') '# positions after  distance filtering: ', npeaks
        ! update packed arrays
        pos_inds   = pack(pos_inds,   mask=selected_pos)
        pos_scores = pack(pos_scores, mask=selected_pos)
        ! update box_scores
        !$omp parallel do schedule(static) collapse(2) default(shared) private(ioff,joff,ithr,is_peak,ipeak) proc_bind(close)
        do ioff = 1,self%nx_offset
            do joff = 1,self%ny_offset
                ithr = omp_get_thread_num() + 1
                is_peak = .false.
                do ipeak = 1,npeaks
                    if( pos_inds(ipeak) == self%inds_offset(ioff,joff) )then
                        is_peak = .true.
                        exit
                    endif
                end do
                if( .not. is_peak ) self%box_scores(ioff,joff) = -1.
            end do
        end do
        !$omp end parallel do
        npeaks = count(self%box_scores >= self%t)
        write(logfhandle,'(a,1x,I5)') '# positions after updating box_scores: ', npeaks
    end subroutine distance_filter

    subroutine remove_outliers( self )
        class(pickgau), intent(inout) :: self
        real, allocatable :: tmp(:)
        type(image) :: boximgs_heap(nthr_glob)
        integer     :: ibox, ithr, npeaks, pos(2), ioff, joff
        logical     :: outside
        real        :: loc_sdevs(self%nx_offset,self%ny_offset), avg, sdev, t
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box, self%smpd_shrink)
        end do
        !$omp parallel do schedule(static) collapse(2) default(shared) private(ioff,joff,ithr,outside,pos) proc_bind(close)
        do ioff = 1,self%nx_offset
            do joff = 1,self%ny_offset
                ithr = omp_get_thread_num() + 1
                pos  = self%positions(self%inds_offset(ioff,joff),:)
                if( self%box_scores(ioff,joff) >= self%t )then
                    call self%mic_shrink%window_slim(pos, self%ldim_box(1), boximgs_heap(ithr), outside)
                    loc_sdevs(ioff,joff) = boximgs_heap(ithr)%avg_loc_sdev(OFFSET)
                else
                    loc_sdevs(ioff,joff) = -1.
                endif
            end do
        end do
        !$omp end parallel do
        tmp = pack(loc_sdevs, mask=loc_sdevs > 0.)
        call avg_sdev(tmp, avg, sdev)
        t = avg + self%ndev * sdev
        npeaks = count(tmp < t)
        write(logfhandle,'(a,1x,I5)') '# positions after  outlier    removal: ', npeaks
        ! update box_scores
        !$omp parallel do schedule(static) collapse(2) default(shared) private(ioff,joff,ithr) proc_bind(close)
        do ioff = 1,self%nx_offset
            do joff = 1,self%ny_offset
                ithr = omp_get_thread_num() + 1
                if( loc_sdevs(ioff,joff) < t )then
                    ! it is a peak
                else
                    ! it is not a peak, update box_scores
                    self%box_scores(ioff,joff) = -1.
                endif
            end do
        end do
        !$omp end parallel do
        npeaks = count(self%box_scores >= self%t)
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%kill
        end do
    end subroutine remove_outliers

    subroutine get_positions( self, pos, smpd_new )
        class(pickgau),       intent(in)    :: self
        integer, allocatable, intent(inout) :: pos(:,:)
        real, optional,       intent(in)    :: smpd_new
        integer, allocatable :: pos_inds(:)
        real    :: scale
        integer :: nbox, ibox
        pos_inds = pack(self%inds_offset(:,:), mask=self%box_scores(:,:) >= self%t)
        nbox     = size(pos_inds)
        if( allocated(pos) ) deallocate(pos)
        allocate( pos(nbox,2),  source=0 )
        if( present(smpd_new) )then
            scale = self%smpd_shrink / smpd_new
        else
            scale = 1.
        endif
        do ibox = 1,nbox
            pos(ibox,:) = nint(scale * real(self%positions(pos_inds(ibox),:)))
        end do
    end subroutine get_positions

    subroutine refine_upscaled( self, pos, smpd_old )
        class(pickgau), intent(inout) :: self
        integer,        intent(in)    :: pos(:,:)
        real,           intent(in)    :: smpd_old
        integer, allocatable :: pos_refined(:,:)
        real,    allocatable :: scores_refined(:)
        type(image) :: boximgs_heap(nthr_glob)
        integer     :: nbox, ibox, jbox, ithr, xrange(2), yrange(2), xind, yind, iref, loc(1)
        real        :: factor, rpos(2), box_score, box_score_trial, scores(self%nrefs), dists(self%nboxes)
        logical     :: outside, l_err
        nbox = size(pos, dim=1)
        allocate(pos_refined(nbox,2),  source= 0 )
        allocate(scores_refined(nbox), source=-1.)
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box, self%smpd_shrink)
        end do
        factor = real(OFFSET) * (smpd_old / self%smpd_shrink)
        if( self%refpick )then
            !$omp parallel do schedule(static) default(shared) proc_bind(close)&
            !$omp private(ibox,rpos,xrange,yrange,box_score,xind,yind,ithr,outside,iref,scores,box_score_trial,l_err)
            do ibox = 1,nbox
                rpos      = real(pos(ibox,:))
                xrange(1) = max(0, nint(rpos(1) - factor))
                xrange(2) = min(self%nx, nint(rpos(1) + factor))
                yrange(1) = max(0, nint(rpos(2) - factor))
                yrange(2) = min(self%ny, nint(rpos(2) + factor))
                box_score = -1.
                do xind = xrange(1),xrange(2)
                    do yind = yrange(1),yrange(2)
                        ithr = omp_get_thread_num() + 1
                        call self%mic_shrink%window_slim([xind,yind], self%ldim_box(1), boximgs_heap(ithr), outside)
                        call boximgs_heap(ithr)%prenorm4real_corr(l_err)
                        if( l_err )then
                            box_score_trial = -1.
                        else
                            do iref = 1,self%nrefs
                                if( self%l_err_refs(iref) )then
                                    scores(iref) = -1.
                                else
                                    scores(iref) = self%boxrefs(iref)%real_corr_prenorm(boximgs_heap(ithr))
                                endif
                            end do
                            box_score_trial = maxval(scores)
                        endif
                        if( box_score_trial > box_score )then
                            pos_refined(ibox,:) = [xind,yind]
                            box_score = box_score_trial
                        endif
                    end do
                end do
                scores_refined(ibox) = box_score
            end do
            !$omp end parallel do
        else
            !$omp parallel do schedule(static) default(shared) proc_bind(close)&
            !$omp private(ibox,rpos,xrange,yrange,box_score,xind,yind,ithr,outside,box_score_trial)
            do ibox = 1,nbox
                rpos      = real(pos(ibox,:))
                xrange(1) = max(0, nint(rpos(1) - factor))
                xrange(2) = min(self%nx, nint(rpos(1) + factor))
                yrange(1) = max(0, nint(rpos(2) - factor))
                yrange(2) = min(self%ny, nint(rpos(2) + factor))
                box_score = -1.
                do xind = xrange(1),xrange(2)
                    do yind = yrange(1),yrange(2)
                        ithr = omp_get_thread_num() + 1
                        call self%mic_shrink%window_slim([xind,yind], self%ldim_box(1), boximgs_heap(ithr), outside)
                        box_score_trial = self%gauref%real_corr_prenorm(boximgs_heap(ithr), self%sxx)
                        if( box_score_trial > box_score )then
                            pos_refined(ibox,:) = [xind,yind]
                            box_score = box_score_trial
                        endif
                    end do
                end do
                scores_refined(ibox) = box_score
            end do
            !$omp end parallel do
        endif
        ! translate refined positions into instance
        self%t = minval(scores_refined)
        self%box_scores(:,:) = -1.
        do ibox = 1,nbox
            !$omp parallel do schedule(static) default(shared) private(jbox) proc_bind(close)
            do jbox = 1,self%nboxes
                dists(jbox) = euclid(real(pos_refined(ibox,:)),real(self%positions(jbox,:)))
            end do
            !$omp end parallel do
            loc = minloc(dists)
            where(self%inds_offset(:,:) == loc(1)) self%box_scores(:,:) = scores_refined(ibox) 
        end do
        ! destruct
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%kill
        end do
    end subroutine refine_upscaled

    subroutine kill( self )
        class(pickgau), intent(inout) :: self
        integer :: iimg
        if( self%exists )then
            call self%mic_shrink%kill
            call self%gauref%kill
            call self%destruct_boximgs
            if( allocated(self%l_mic_mask)  ) deallocate(self%l_mic_mask)
            if( allocated(self%l_err_refs)  ) deallocate(self%l_err_refs)
            if( allocated(self%positions)   ) deallocate(self%positions)
            if( allocated(self%inds_offset) ) deallocate(self%inds_offset)
            if( allocated(self%box_scores)  ) deallocate(self%box_scores)
            if( allocated(self%boxrefs) )then
                do iimg = 1,size(self%boxrefs)
                    call self%boxrefs(iimg)%kill
                end do
                deallocate(self%boxrefs)
            endif
            self%refpick = .false.
            self%exists  = .false.
        endif
    end subroutine kill

end module simple_pickgau
