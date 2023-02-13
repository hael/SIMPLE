module simple_picker_utils
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,          only: image
use simple_radial_medians, only: radial_medians
use simple_segmentation
implicit none

public :: picker_utils
private
#include "simple_local_flags.inc"

real,    parameter :: SMPD_SHRINK1 = 4.0, SMPD_SHRINK2 = 2.0, GAUSIG = 5.
integer, parameter :: OFFSET       = 3, OFFSET_UB = 2 * OFFSET
real,    parameter :: DIST_THRES   = real(OFFSET_UB)
logical, parameter :: L_WRITE = .true.

type picker_utils
    private
    integer                  :: ldim_raw(3) = 0 , ldim_shrink1(3) = 0 , ldim_shrink2(3) = 0
    integer                  :: ldim_box(3) = 0 , ldim_box1(3)    = 0 , ldim_box2(3)    = 0
    real                     :: smpd_raw    = 0., smpd_shrink1    = 0., smpd_shrink2    = 0.
    real                     :: maxdiam = 0., sig_shrink1 = 0., sig_shrink2 = 0., sxx_shrink1 = 0.
    real                     :: sxx_shrink2 = 0.
    integer(dp)              :: nboxes1 = 0, nboxes2 = 0, nboxes_ub = 0
    integer                  :: nx1 = 0, ny1 = 0, nx2 = 0, ny2 = 0, nx_offset  = 0, ny_offset = 0
    type(image)              :: mic_shrink1, mic_shrink2, imgau_shrink1, imgau_shrink2
    type(image), pointer     :: mic_raw => null()
    type(image), allocatable :: boximgs1(:), boximgs2(:)
    integer,     allocatable :: positions1(:,:), positions2(:,:), inds_offset(:,:)
    real,        allocatable :: box_scores1(:)
    logical                  :: exists = .false.
  contains
    procedure          :: new
    procedure, private :: set_positions_1
    procedure, private :: set_positions_2
    generic            :: set_positions => set_positions_1, set_positions_2
    procedure, private :: set_pos_priv
    procedure          :: gauconv_mics
    procedure, private :: extract_boximgs1
    procedure, private :: extract_boximgs2
    procedure          :: analyze_boximgs1
    procedure          :: distance_filter
    procedure          :: refine_positions
end type picker_utils

contains

    subroutine new( self, mic, smpd, maxdiam, bp_lp )
        class(picker_utils),  intent(inout) :: self
        class(image), target, intent(in)    :: mic
        real,                 intent(in)    :: smpd     !< sampling distance in A
        real,                 intent(in)    :: maxdiam  !< maximum diameter in A
        real, optional,       intent(in)    :: bp_lp(2) !< high- and low-pass limits in A
        real :: scale1, scale2, pixrad_shrink1, pixrad_shrink2
        ! set raw micrograph info
        self%ldim_raw = mic%get_ldim()
        if( self%ldim_raw(3) /= 1 ) THROW_HARD('Only for 2D images')
        self%smpd_raw = smpd
        self%mic_raw  => mic
        ! set shrunken logical dimensions
        scale1 = self%smpd_raw / SMPD_SHRINK1
        scale2 = self%smpd_raw / SMPD_SHRINK2
        self%ldim_shrink1(1) = round2even(real(self%ldim_raw(1)) * scale1)
        self%ldim_shrink1(2) = round2even(real(self%ldim_raw(2)) * scale1)
        self%ldim_shrink1(3) = 1
        self%ldim_shrink2(1) = round2even(real(self%ldim_raw(1)) * scale2)
        self%ldim_shrink2(2) = round2even(real(self%ldim_raw(2)) * scale2)
        self%ldim_shrink2(3) = 1
        ! set logical dimensions of boxes
        self%maxdiam         = maxdiam
        self%ldim_box(1)     = round2even(self%maxdiam / self%smpd_raw)
        self%ldim_box(2)     = self%ldim_box(1)
        self%ldim_box(3)     = 1
        self%ldim_box1(1)    = round2even(self%maxdiam / SMPD_SHRINK1)
        self%ldim_box1(2)    = self%ldim_box1(1)
        self%ldim_box1(3)    = 1
        self%ldim_box2(1)    = round2even(self%maxdiam / SMPD_SHRINK2)
        self%ldim_box2(2)    = self%ldim_box2(1)
        self%ldim_box2(3)    = 1
        ! set # pixels in x/y for both box sizes
        self%nx1 = self%ldim_shrink1(1) - self%ldim_box1(1)
        self%ny1 = self%ldim_shrink1(2) - self%ldim_box1(2)
        self%nx2 = self%ldim_shrink2(1) - self%ldim_box2(1)
        self%ny2 = self%ldim_shrink2(2) - self%ldim_box2(2)
        ! set Gaussians
        pixrad_shrink1       = (self%maxdiam / 2.) / SMPD_SHRINK1
        pixrad_shrink2       = (self%maxdiam / 2.) / SMPD_SHRINK2
        self%sig_shrink1     = pixrad_shrink1 / GAUSIG
        self%sig_shrink2     = pixrad_shrink2 / GAUSIG
        call self%imgau_shrink1%new(self%ldim_box1, SMPD_SHRINK1)
        call self%imgau_shrink1%gauimg2D(self%sig_shrink1, self%sig_shrink1)
        call self%imgau_shrink2%new(self%ldim_box2, SMPD_SHRINK2)
        call self%imgau_shrink2%gauimg2D(self%sig_shrink2, self%sig_shrink2)
        call self%imgau_shrink1%prenorm4real_corr(self%sxx_shrink1)
        call self%imgau_shrink2%prenorm4real_corr(self%sxx_shrink2)
        ! shrink micrograph
        call self%mic_shrink1%new(self%ldim_shrink1, SMPD_SHRINK1)
        call self%mic_shrink2%new(self%ldim_shrink2, SMPD_SHRINK2)
        call self%mic_shrink1%set_ft(.true.)
        call self%mic_shrink2%set_ft(.true.)
        call self%mic_raw%fft
        call self%mic_raw%clip(self%mic_shrink1)
        call self%mic_raw%clip(self%mic_shrink2)
        if( present(bp_lp) )then
            call self%mic_shrink1%bp(bp_lp(1),bp_lp(2))
            call self%mic_shrink2%bp(bp_lp(1),bp_lp(2))
        endif
        ! flip contrast (assuming black particle contrast on input)
        call self%mic_shrink1%mul(-1.)
        call self%mic_shrink2%mul(-1.)
        ! back to real-space
        call self%mic_raw%ifft
        call self%mic_shrink1%ifft
        call self%mic_shrink2%ifft
        if( L_WRITE )then
            call self%mic_shrink1%write('mic_shrink1.mrc')
            call self%mic_shrink2%write('mic_shrink2.mrc')
        endif
        self%exists = .true.
    end subroutine new

    subroutine set_positions_1( self )
        class(picker_utils), intent(inout) :: self
        call self%set_pos_priv
    end subroutine set_positions_1

    subroutine set_positions_2( self, box_raw, box12 )
        class(picker_utils), intent(inout) :: self
        real,                intent(in)    :: box_raw, box12(2)
        ! set logical dimensions of boxes
        self%ldim_box(1)  = box_raw
        self%ldim_box(2)  = self%ldim_box(1)
        self%ldim_box(3)  = 1
        self%ldim_box1(1) = box12(1)
        self%ldim_box1(2) = self%ldim_box1(1)
        self%ldim_box1(3) = 1
        self%ldim_box2(1) = box12(2)
        self%ldim_box2(2) = self%ldim_box2(1)
        self%ldim_box2(3) = 1
        call self%set_pos_priv
    end subroutine set_positions_2

    subroutine set_pos_priv( self )
        class(picker_utils), intent(inout) :: self
        integer :: xind, yind
        ! set # pixels in x/y for both box sizes
        self%nx1 = self%ldim_shrink1(1) - self%ldim_box1(1)
        self%ny1 = self%ldim_shrink1(2) - self%ldim_box1(2)
        self%nx2 = self%ldim_shrink2(1) - self%ldim_box2(1)
        self%ny2 = self%ldim_shrink2(2) - self%ldim_box2(2)
        ! count # boxes
        self%nboxes1   = 0
        self%nx_offset = 0
        do xind = 0,self%nx1,OFFSET
            self%nx_offset = self%nx_offset + 1
            self%ny_offset = 0
            do yind = 0,self%ny1,OFFSET
                self%nboxes1   = self%nboxes1   + 1
                self%ny_offset = self%ny_offset + 1
            end do
        end do
        ! count # boxes, upper bound
        self%nboxes_ub = 0
        do xind = 0,self%nx1,OFFSET_UB
            do yind = 0,self%ny1,OFFSET_UB
                self%nboxes_ub   = self%nboxes_ub + 1
            end do
        end do
        ! allocate and set positions1 
        if( allocated(self%positions1) ) deallocate(self%positions1)
        allocate(self%positions1(self%nboxes1,2), source=0)
        ! allocate and set inds_offset
        if( allocated(self%inds_offset) ) deallocate(self%inds_offset)
        allocate(self%inds_offset(self%nx_offset,self%ny_offset), source=0)
        ! calculate total # boxes
        self%nboxes1   = 0
        self%nx_offset = 0
        do xind = 0,self%nx1,OFFSET
            self%nx_offset = self%nx_offset + 1
            self%ny_offset = 0
            do yind = 0,self%ny1,OFFSET
                self%nboxes1   = self%nboxes1 + 1
                self%ny_offset = self%ny_offset + 1
                self%positions1(self%nboxes1,:) = [xind,yind]
                self%inds_offset(self%nx_offset,self%ny_offset) = self%nboxes1
            end do
        end do
    end subroutine set_pos_priv

    subroutine gauconv_mics( self )
        class(picker_utils), intent(inout) :: self
        type(image) :: img
        ! denoise mic_shrink1
        call img%new(self%ldim_shrink1, SMPD_SHRINK1)
        call img%gauimg2D(self%sig_shrink1,self%sig_shrink1)
        call img%fft
        call self%mic_shrink1%fft
        call self%mic_shrink1%mul(img)
        call self%mic_shrink1%ifft
        if( L_WRITE )then
            call self%mic_shrink1%write('gauconv_mic_shrink1.mrc')
            call img%zero_and_unflag_ft
            call sauvola(self%mic_shrink1, OFFSET, img)
            call img%write('sauvola.mrc')
            call img%real_space_filter(OFFSET * 3, 'average')
            call img%write('sauvola_rfilt.mrc')
            call otsu_img(img)
            call img%write('sauvola_otsu.mrc')
        endif
        ! denoise mic_shrink2
        call img%new(self%ldim_shrink2, SMPD_SHRINK2)
        call img%gauimg2D(self%sig_shrink2,self%sig_shrink2)
        call img%fft
        call self%mic_shrink2%fft
        call self%mic_shrink2%mul(img)
        call self%mic_shrink2%ifft
        if( L_WRITE )then
            call self%mic_shrink2%write('gauconv_mic_shrink2.mrc')
        endif
        call img%kill
    end subroutine gauconv_mics

    subroutine extract_boximgs1( self )
        class(picker_utils), intent(inout) :: self
        integer :: ibox
        logical :: outside
        if( .not. allocated(self%positions1) ) THROW_HARD('positions need to be set before constructing boximgs1')
        if( allocated(self%boximgs1) )then
            do ibox = 1,self%nboxes1
                call self%boximgs1(ibox)%kill
            end do
            deallocate(self%boximgs1)
        endif
        allocate(self%boximgs1(self%nboxes1))
        !$omp parallel do schedule(static) default(shared) private(ibox,outside) proc_bind(close)
        do ibox = 1,self%nboxes1
            call self%boximgs1(ibox)%new(self%ldim_box1, SMPD_SHRINK1)
            call self%mic_shrink1%window_slim(self%positions1(ibox,:), self%ldim_box1(1), self%boximgs1(ibox), outside)
        end do
        !$omp end parallel do
    end subroutine extract_boximgs1

    subroutine extract_boximgs2( self )
        class(picker_utils), intent(inout) :: self
        integer :: ibox
        logical :: outside
        if( .not. allocated(self%positions2) ) THROW_HARD('positions need to be set before constructing boximgs2')
        if( allocated(self%boximgs2) )then
            do ibox = 1,self%nboxes2
                call self%boximgs2(ibox)%kill
            end do
            deallocate(self%boximgs2)
        endif
        allocate(self%boximgs2(self%nboxes2))
        !$omp parallel do schedule(static) default(shared) private(ibox,outside) proc_bind(close)
        do ibox = 1,self%nboxes2
            call self%boximgs2(ibox)%new(self%ldim_box2, SMPD_SHRINK2)
            call self%mic_shrink2%window_slim(self%positions2(ibox,:), self%ldim_box2(1), self%boximgs2(ibox), outside)
        end do
        !$omp end parallel do
    end subroutine extract_boximgs2

    subroutine analyze_boximgs1( self )
        class(picker_utils), intent(inout) :: self
        integer, allocatable :: positions_tmp(:,:)
        real,    allocatable :: tmp(:)
        real        :: box_scores(self%nx_offset,self%ny_offset,1), t
        logical     :: is_peak(self%nx_offset,self%ny_offset,1), outside
        integer     :: ioff, joff, npeaks, cnt, ithr
        type(image) :: boximgs_heap(nthr_glob)
        ! calculate box_scores
        if( .not. allocated(self%positions1) ) THROW_HARD('positions1 need to be set')
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box1, SMPD_SHRINK1)
        end do
        !$omp parallel do schedule(static) collapse(2) default(shared) private(ioff,joff,ithr,outside) proc_bind(close)
        do ioff = 1,self%nx_offset
            do joff = 1,self%ny_offset
                ithr = omp_get_thread_num() + 1
                call self%mic_shrink1%window_slim(self%positions1(self%inds_offset(ioff,joff),:),&
                &self%ldim_box1(1), boximgs_heap(ithr), outside)
                box_scores(ioff,joff,1) = self%imgau_shrink1%real_corr_prenorm(boximgs_heap(ithr), self%sxx_shrink1)
            end do
        end do
        !$omp end parallel do
        tmp = pack(box_scores, mask=.true.)
        call detect_peak_thres(self%nx_offset * self%ny_offset, int(self%nboxes_ub), tmp, t)
        deallocate(tmp)
        is_peak = .false.
        do ioff = 1,self%nx_offset
            do joff = 1,self%ny_offset
                if( box_scores(ioff,joff,1) >= t ) is_peak(ioff,joff,1) = .true.
            end do
        end do
        npeaks = count(is_peak)
        write(logfhandle,'(a,1x,I5)') '# positions considered               : ', self%nx_offset * self%ny_offset
        ! update positions1 and box_scores1
        if( allocated(self%box_scores1) ) deallocate(self%box_scores1)
        allocate(positions_tmp(npeaks,2), self%box_scores1(npeaks))
        positions_tmp    = 0
        self%box_scores1 = 0.
        cnt              = 0
        do ioff = 1,self%nx_offset
            do joff = 1,self%ny_offset
                if( is_peak(ioff,joff,1) )then
                    cnt = cnt + 1
                    positions_tmp(cnt,:)  = self%positions1(self%inds_offset(ioff,joff),:)
                    self%box_scores1(cnt) = box_scores(ioff,joff,1)
                endif
            end do
        end do
        deallocate(self%positions1)
        self%nboxes1 = npeaks
        allocate(self%positions1(self%nboxes1,2), source=positions_tmp)
        deallocate(positions_tmp)
        if( L_WRITE )then
            call self%extract_boximgs1
            call write_boximgs(int(self%nboxes1), self%boximgs1, 'ptcls_before_distfilt.mrcs' )
        endif
    end subroutine analyze_boximgs1

    subroutine distance_filter( self )
        class(picker_utils), intent(inout) :: self
        integer, allocatable :: positions_tmp(:,:)
        real,    allocatable :: box_scores1_tmp(:)
        integer :: ibox, jbox, loc(1), npeaks, cnt
        real    :: dist
        logical :: mask(self%nboxes1), selected_pos1(self%nboxes1)
        selected_pos1 = .true.
        do ibox = 1,self%nboxes1
            mask = .false.
            !$omp parallel do schedule(static) default(shared) private(jbox,dist) proc_bind(close)
            do jbox = 1,self%nboxes1
                dist = euclid(real(self%positions1(ibox,:)),real(self%positions1(jbox,:)))
                if( dist <= DIST_THRES ) mask(jbox) = .true.
            end do
            !$omp end parallel do
            ! find best match in the neigh
            loc = maxloc(self%box_scores1, mask=mask)
            ! eliminate all but the best
            mask(loc(1)) = .false.
            where( mask ) selected_pos1 = .false.
        end do
        npeaks = count(selected_pos1)
        write(logfhandle,'(a,1x,I5)') '# positions before distance filtering: ', self%nboxes1
        write(logfhandle,'(a,1x,I5)') '# positions after  distance filtering: ', npeaks
        ! update positions1 and box_scores1
        allocate(positions_tmp(npeaks,2), box_scores1_tmp(npeaks))
        positions_tmp   = 0
        box_scores1_tmp = 0.
        cnt             = 0
        do ibox = 1,self%nboxes1
            if( selected_pos1(ibox) )then
                cnt = cnt + 1
                positions_tmp(cnt,:) = self%positions1(ibox,:)
                box_scores1_tmp(cnt) = self%box_scores1(ibox)
            endif
        end do
        deallocate(self%positions1, self%box_scores1)
        self%nboxes1 = npeaks
        allocate(self%positions1(self%nboxes1,2), source=positions_tmp)
        allocate(self%box_scores1(self%nboxes1),  source=box_scores1_tmp)
        deallocate(positions_tmp, box_scores1_tmp)
        if( L_WRITE )then
            call self%extract_boximgs1
            call write_boximgs(int(self%nboxes1), self%boximgs1, 'ptcls_after_distfilt.mrcs' )
        endif
    end subroutine distance_filter

    subroutine refine_positions( self )
        class(picker_utils), intent(inout) :: self
        integer     :: ibox, xrange(2), yrange(2), xind, yind, ithr, old_pos(2)
        real        :: box_score, box_score_trial, factor, rpos(2)
        logical     :: outside
        type(image) :: boximgs_heap(nthr_glob)
        if( .not. allocated(self%positions1) ) THROW_HARD('positions1 need to be set')
        self%nboxes2 = self%nboxes1
        allocate(self%positions2(self%nboxes2,2), source=nint((SMPD_SHRINK1/SMPD_SHRINK2) * real(self%positions1)))
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box2, SMPD_SHRINK2)
        end do
        factor = real(OFFSET) * (SMPD_SHRINK1 / SMPD_SHRINK2)
        !$omp parallel do schedule(static) default(shared)&
        !$omp private(ibox,rpos,xrange,yrange,box_score,xind,yind,ithr,outside,box_score_trial) proc_bind(close)
        do ibox= 1,self%nboxes2
            rpos      = real(self%positions2(ibox,:))
            xrange(1) = max(0,        nint(rpos(1) - factor))
            xrange(2) = min(self%nx2, nint(rpos(1) + factor))
            yrange(1) = max(0,        nint(rpos(2) - factor))
            yrange(2) = min(self%ny2, nint(rpos(2) + factor))
            box_score = -1
            do xind = xrange(1),xrange(2)
                do yind = yrange(1),yrange(2)
                    ithr = omp_get_thread_num() + 1
                    call self%mic_shrink2%window_slim([xind,yind], self%ldim_box2(1), boximgs_heap(ithr), outside)
                    box_score_trial = self%imgau_shrink2%real_corr_prenorm(boximgs_heap(ithr), self%sxx_shrink2)
                    if( box_score_trial > box_score )then
                        self%positions2(ibox,:) = [xind,yind]
                        box_score = box_score_trial
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        if( L_WRITE )then
            call self%extract_boximgs2
            call write_boximgs(int(self%nboxes2), self%boximgs2, 'ptcls_refined.mrcs' )
        endif
    end subroutine refine_positions

    ! utilities

    subroutine write_boximgs( n, boximgs, fname, mask )
        integer,           intent(in)    :: n
        class(image),      intent(inout) :: boximgs(n)
        character(len=*),  intent(in)    :: fname
        logical, optional, intent(in)    :: mask(n)
        integer :: ibox, cnt
        if( present(mask) )then
            cnt = 0
            do ibox = 1,n
                if( mask(ibox) )then
                    cnt = cnt + 1
                    call boximgs(ibox)%write(fname, cnt)
                endif
            end do
        else
            do ibox = 1,n
                call boximgs(ibox)%write(fname, ibox)
            end do
        endif
    end subroutine write_boximgs

    subroutine write_boxfile( n, coordinates, box, mask, fname )
        integer,          intent(in) :: n, coordinates(n,2), box
        logical,          intent(in) :: mask(n)
        character(len=*), intent(in) :: fname
        integer :: funit, ibox, iostat
        call fopen(funit, status='REPLACE', action='WRITE', file=trim(adjustl(fname)), iostat=iostat)
        call fileiochk('simple_picker_utils; write_boxfile ', iostat)
        do ibox = 1,n
            if( mask(ibox) )then
                write(funit,'(I7,I7,I7,I7,I7)') coordinates(1,ibox), coordinates(2,ibox), box, box, -3
            endif
        end do
        call fclose(funit)
    end subroutine write_boxfile

end module simple_picker_utils
