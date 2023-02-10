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
integer, parameter :: OFFSET  = 3, MAXKMIT  = 20
logical, parameter :: L_WRITE = .true.

type picker_utils
    private
    integer                  :: ldim_raw(3) = 0 , ldim_shrink1(3) = 0 , ldim_shrink2(3) = 0
    integer                  :: ldim_box(3) = 0 , ldim_box1(3)    = 0 , ldim_box2(3)    = 0
    real                     :: smpd_raw    = 0., smpd_shrink1    = 0., smpd_shrink2    = 0.
    integer(dp)              :: nboxes1 = 0, nboxes2 = 0, nboxes_ub = 0
    integer                  :: nx1 = 0, ny1 = 0, nx2 = 0, ny2 = 0, nx_offset  = 0, ny_offset = 0
    type(image)              :: mic_shrink1, mic_shrink2
    type(image), pointer     :: mic_raw => null()
    type(image), allocatable :: boximgs1(:), boximgs2(:)
    integer,     allocatable :: positions1(:,:), inds_offset(:,:)
    logical,     allocatable :: is_peak(:,:,:)
  contains
    procedure          :: set_mics
    procedure, private :: set_box_ldims
    procedure, private :: set_positions_1
    procedure, private :: set_positions_2
    generic            :: set_positions => set_positions_1, set_positions_2
    procedure, private :: set_pos_priv
    procedure          :: gauconv_mic_shrink1
    procedure, private :: extract_boximgs1
    procedure          :: analyze_boximgs1
end type picker_utils

contains

    subroutine set_mics( self, mic, smpd, bp_lp )
        class(picker_utils),  intent(inout) :: self
        class(image), target, intent(in)    :: mic
        real,                 intent(in)    :: smpd     !< sampling distance in A
        real, optional,       intent(in)    :: bp_lp(2) !< high- and low-pass limits in A
        real :: scale1, scale2
        ! set raw micrograph info
        self%ldim_raw = mic%get_ldim()
        if( self%ldim_raw(3) /= 1 ) THROW_HARD('Only for 2D images')
        self%smpd_raw = smpd
        self%mic_raw  => mic
        ! shrink micrograph
        scale1 = self%smpd_raw / SMPD_SHRINK1
        scale2 = self%smpd_raw / SMPD_SHRINK2
        self%ldim_shrink1(1) = round2even(real(self%ldim_raw(1)) * scale1)
        self%ldim_shrink1(2) = round2even(real(self%ldim_raw(2)) * scale1)
        self%ldim_shrink1(3) = 1
        self%ldim_shrink2(1) = round2even(real(self%ldim_raw(1)) * scale2)
        self%ldim_shrink2(2) = round2even(real(self%ldim_raw(2)) * scale2)
        self%ldim_shrink2(3) = 1
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
    end subroutine set_mics

    subroutine set_box_ldims( self, maxdiam )
        class(picker_utils), intent(inout) :: self
        real,                intent(in)    :: maxdiam !< maximum diameter in A
        ! set logical dimensions of boxes
        self%ldim_box(1)  = round2even(maxdiam / self%smpd_raw)
        self%ldim_box(2)  = self%ldim_box(1)
        self%ldim_box(3)  = 1
        self%ldim_box1(1) = round2even(maxdiam / SMPD_SHRINK1)
        self%ldim_box1(2) = self%ldim_box1(1)
        self%ldim_box1(3) = 1
        self%ldim_box2(1) = round2even(maxdiam / SMPD_SHRINK2)
        self%ldim_box2(2) = self%ldim_box2(1)
        self%ldim_box2(3) = 1
        ! set # pixels in x/y for both box sizes
        self%nx1 = self%ldim_shrink1(1) - self%ldim_box1(1)
        self%ny1 = self%ldim_shrink1(2) - self%ldim_box1(2)
        self%nx2 = self%ldim_shrink2(1) - self%ldim_box2(1)
        self%ny2 = self%ldim_shrink2(2) - self%ldim_box2(2)
    end subroutine set_box_ldims

    subroutine set_positions_1( self, maxdiam )
        class(picker_utils), intent(inout) :: self
        real,                intent(in)    :: maxdiam !< maximum diameter in A
        call self%set_box_ldims(maxdiam)
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
        do xind = 0,self%nx1,3 * OFFSET
            do yind = 0,self%ny1,3 * OFFSET
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

    subroutine gauconv_mic_shrink1( self, maxdiam )
        class(picker_utils), intent(inout) :: self
        real,                intent(in)    :: maxdiam
        type(image) :: img
        real        :: pix_rad, sig
        pix_rad = (maxdiam / 2.) / SMPD_SHRINK1 
        sig     = pix_rad / GAUSIG
        call img%new(self%ldim_shrink1, SMPD_SHRINK1)
        call img%gauimg2D(sig, sig)
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
        call img%kill
    end subroutine gauconv_mic_shrink1

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

    subroutine analyze_boximgs1( self, maxdiam )
        class(picker_utils), intent(inout) :: self
        real,                intent(in)    :: maxdiam
        integer, allocatable :: positions_tmp(:,:)
        real,    allocatable :: tmp(:)
        real        :: box_scores(self%nx_offset,self%ny_offset,1), t, sxx
        logical     :: is_peak(self%nx_offset,self%ny_offset,1), outside
        integer     :: ioff, joff, npeaks, cnt, ithr
        type(image) :: img, boximgs_heap(nthr_glob)
        real        :: pix_rad, sig
        pix_rad = (maxdiam / 2.) / SMPD_SHRINK1 
        sig     = pix_rad / GAUSIG
        call img%new(self%ldim_box1, SMPD_SHRINK1)
        call img%gauimg2D(sig, sig)
        call img%prenorm4real_corr(sxx)
        ! call img%vis
        ! calculate box_scores
        if( .not. allocated(self%positions1) ) THROW_HARD('positions need to be set before constructing boximgs1')
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box1, SMPD_SHRINK1)
        end do
        !$omp parallel do schedule(static) collapse(2) default(shared) private(ioff,joff,ithr,outside) proc_bind(close)
        do ioff = 1,self%nx_offset
            do joff = 1,self%ny_offset
                ithr = omp_get_thread_num() + 1
                call self%mic_shrink1%window_slim(self%positions1(self%inds_offset(ioff,joff),:),&
                &self%ldim_box1(1), boximgs_heap(ithr), outside)
                box_scores(ioff,joff,1) = img%real_corr_prenorm(boximgs_heap(ithr), sxx)
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

        print *, 'npoints ', self%nx_offset * self%ny_offset
        print *, 'npeaks  ', npeaks

        ! update positions1
        allocate(positions_tmp(npeaks,2), source=0)
        cnt = 0
        do ioff = 1,self%nx_offset
            do joff = 1,self%ny_offset
                if( is_peak(ioff,joff,1) )then
                    cnt = cnt + 1
                    positions_tmp(cnt,:) = self%positions1(self%inds_offset(ioff,joff),:)
                endif
            end do
        end do
        deallocate(self%positions1)
        self%nboxes1 = npeaks
        allocate(self%positions1(self%nboxes1,2), source=positions_tmp)
        deallocate(positions_tmp)
        if( L_WRITE )then
            call self%extract_boximgs1
            call write_boximgs(int(self%nboxes1), self%boximgs1, 'peak_particles.mrcs' )
        endif
    end subroutine analyze_boximgs1

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
