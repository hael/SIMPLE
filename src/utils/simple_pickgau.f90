module simple_pickgau
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
implicit none

public :: pickgau
private
#include "simple_local_flags.inc"

! class constants
real,    parameter :: GAUSIG = 5., BOX_EXP_FAC = 0.111
integer, parameter :: OFFSET      = 3, OFFSET_UB = 2 * OFFSET
real,    parameter :: DIST_THRES1 = real(OFFSET), DIST_THRES2 = real(4*OFFSET)
logical, parameter :: L_WRITE     = .false.

! class variables
integer     :: ldim_raw(3), ldim_raw_box(3)
real        :: smpd_raw
type(image) :: mic_raw

! instance
type pickgau
    private
    real                     :: smpd_shrink = 0., maxdiam = 0., sig = 0., sxx = 0., t = 0.
    integer                  :: ldim(3), ldim_box(3), nboxes = 0, nboxes_ub = 0
    integer                  :: nx = 0, ny = 0, nx_offset = 0, ny_offset = 0, npeaks = 0
    type(image)              :: mic_shrink, gauref
    type(image), allocatable :: boximgs(:)
    logical,     allocatable :: l_mic_mask(:,:)
    integer,     allocatable :: positions(:,:), inds_offset(:,:)
    real,        allocatable :: box_scores(:,:)
    logical                  :: exists = .false.
contains
    procedure :: new
    procedure :: extract_boximgs
    procedure :: destruct_boximgs
    procedure :: gaumatch_boximgs
    procedure :: detect_peaks
    procedure :: center_filter
    procedure :: kill
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

    subroutine new( self, pcontrast, smpd_shrink, moldiam, moldiam_max )
        class(pickgau),   intent(inout) :: self
        character(len=*), intent(in)    :: pcontrast
        real,             intent(in)    :: smpd_shrink, moldiam, moldiam_max
        character(len=:), allocatable   :: numstr
        integer     :: ldim_pd(3), box_max, xind, yind
        type(image) :: mic_pad, gauimg
        real        :: hpfreq, scale, maxdiam_max
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
        self%maxdiam     = moldiam + moldiam * BOX_EXP_FAC
        ldim_raw_box(1)  = round2even(self%maxdiam / smpd_raw)
        ldim_raw_box(2)  = ldim_raw_box(1)
        ldim_raw_box(3)  = 1
        self%ldim_box(1) = round2even(self%maxdiam / self%smpd_shrink)
        self%ldim_box(2) = self%ldim_box(1)
        self%ldim_box(3) = 1
        ! set # pixels in x/y based on maximum expected box size
        maxdiam_max      = moldiam_max + moldiam_max * BOX_EXP_FAC
        box_max          = round2even(maxdiam_max / self%smpd_shrink)
        self%nx          = self%ldim(1) - box_max
        self%ny          = self%ldim(2) - box_max
        ! set Gaussian
        self%sig = ((self%maxdiam / 2.) / self%smpd_shrink) / GAUSIG
        call self%gauref%new(self%ldim_box, self%smpd_shrink)
        call self%gauref%gauimg2D(self%sig, self%sig)
        call self%gauref%prenorm4real_corr(self%sxx)
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
            numstr = int2str(nint(moldiam))
            call self%mic_shrink%write('mic_shrink_moldiam'//numstr//'.mrc')
            call self%gauref%write('gauref_moldiam'//numstr//'.mrc')
        endif
        ! allocate and set l_mic_mask
        allocate(self%l_mic_mask(self%ldim(1),self%ldim(2)), source = .true.)
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
        allocate(self%box_scores(self%nx_offset,self%ny_offset), source = -1.)
        ! denoise mic_shrink
        call gauimg%new(self%ldim, self%smpd_shrink)
        call gauimg%gauimg2D(self%sig,self%sig)
        call gauimg%fft
        call self%mic_shrink%fft
        call self%mic_shrink%mul(gauimg)
        call self%mic_shrink%ifft
        if( L_WRITE )then
            call self%mic_shrink%write('gauconv_mic_shrink_moldiam'//numstr//'.mrc')
        endif
        call gauimg%kill
        self%exists = .true.
    end subroutine new

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
        write(logfhandle,'(a,1x,I5)') '# peaks detected:                    ', self%npeaks
    end subroutine detect_peaks

     subroutine center_filter( self )
        class(pickgau), intent(inout) :: self
        type(image) :: boximgs_heap(nthr_glob)
        integer     :: ioff, joff, ithr, npeaks
        real        :: scores_cen(self%nx_offset,self%ny_offset)
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box, self%smpd_shrink)
        end do
        !$omp parallel do schedule(static) collapse(2) default(shared) private(ioff,joff) proc_bind(close)
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
        write(logfhandle,'(a,1x,I5)') '# positions before center filtering: ', count(self%box_scores >= self%t)
        write(logfhandle,'(a,1x,I5)') '# positions after  center filtering: ', npeaks
        ! modify box_scores and update npeaks
        where( scores_cen <= real(OFFSET) )
            ! there's a peak
        elsewhere
            self%box_scores = -1.
        endwhere
        self%npeaks = count(self%box_scores >= self%t)
        write(logfhandle,'(a,1x,I5)') '# npeaks    after  center filtering: ', self%npeaks
    end subroutine center_filter

    subroutine kill( self )
        class(pickgau), intent(inout) :: self
        if( self%exists )then
            call self%mic_shrink%kill
            call self%gauref%kill
            call self%destruct_boximgs
            if( allocated(self%l_mic_mask)  ) deallocate(self%l_mic_mask)
            if( allocated(self%positions)   ) deallocate(self%positions)
            if( allocated(self%inds_offset) ) deallocate(self%inds_offset)
            if( allocated(self%box_scores)  ) deallocate(self%box_scores)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_pickgau
