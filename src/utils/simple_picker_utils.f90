module simple_picker_utils
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image, only: image
implicit none

public :: picker_utils
private
#include "simple_local_flags.inc"

real,    parameter :: SMPD_SHRINK1 = 4.0, SMPD_SHRINK2 = 2.0, GAUSIG = 5.
integer, parameter :: OFFSET       = 3,   OFFSET_UB    = 2 * OFFSET
integer, parameter :: MAXNREFS     = 100
real,    parameter :: DIST_THRES   = real(OFFSET), NDEV_DEFAULT = 2.5
logical, parameter :: L_WRITE = .true.

type picker_utils
    private
    integer                       :: ldim_raw(3) = 0 , ldim_shrink1(3) = 0 , ldim_shrink2(3) = 0
    integer                       :: ldim_box(3) = 0 , ldim_box1(3)    = 0 , ldim_box2(3)    = 0
    real                          :: smpd_raw    = 0., smpd_shrink1    = 0., smpd_shrink2    = 0.
    real                          :: maxdiam = 0., sig_shrink1 = 0., sig_shrink2 = 0., sxx_shrink1 = 0.
    real                          :: sxx_shrink2 = 0., ndev = 0.
    integer                       :: nboxes1 = 0, nboxes2 = 0, nboxes_ub = 0, nrefs = 1
    integer                       :: nx1 = 0, ny1 = 0, nx2 = 0, ny2 = 0, nx_offset  = 0, ny_offset = 0
    type(image)                   :: mic_raw, mic_shrink1, mic_shrink2, imgau_shrink1, imgau_shrink2
    type(image),      allocatable :: boximgs1(:), boximgs2(:), boxrefs(:,:)
    integer,          allocatable :: positions1(:,:), positions2(:,:), inds_offset(:,:)
    real,             allocatable :: box_scores1(:), sxx_refs(:,:)
    character(len=LONGSTRLEN)     :: boxname
    character(len=:), allocatable :: fbody
    logical                       :: refpick = .false.
    logical                       :: hybrid  = .false.
    logical                       :: exists  = .false.
  contains
    procedure          :: new
    procedure          :: set_refs
    procedure          :: exec_picker
    procedure, private :: set_positions_1
    procedure, private :: set_positions_2
    generic            :: set_positions => set_positions_1, set_positions_2
    procedure, private :: set_pos_priv
    procedure          :: gauconv_mics
    procedure, private :: extract_boximgs1
    procedure, private :: extract_boximgs2
    procedure          :: analyze_boximgs1
    procedure          :: center_filter
    procedure          :: distance_filter
    procedure          :: refine_positions
    procedure          :: remove_outliers
    procedure          :: peak_vs_nonpeak_stats
    procedure          :: kill
end type picker_utils

contains

    subroutine new( self, micname, pcontrast, smpd, moldiam, ndev )
        class(picker_utils),  intent(inout) :: self
        character(len=*),     intent(in)    :: micname, pcontrast
        real,                 intent(in)    :: smpd    !< sampling distance in A
        real,                 intent(in)    :: moldiam !< maximum diameter in A
        real, optional,       intent(in)    :: ndev    !< # std devs for outlier detection
        character(len=:), allocatable :: ext
        real    :: scale1, scale2, pixrad_shrink1, pixrad_shrink2
        integer :: nframes
        if( self%exists ) call self%kill
        ! set raw micrograph info
        call find_ldim_nptcls(micname, self%ldim_raw, nframes)
        if( self%ldim_raw(3) /= 1 .or. nframes /= 1 ) THROW_HARD('Only for 2D images')
        self%smpd_raw = smpd
        call self%mic_raw%new(self%ldim_raw, self%smpd_raw)
        call self%mic_raw%read(micname)
        ! set fbody & boxname
        ext          = fname2ext(trim(micname))
        self%fbody   = trim(get_fbody(basename(trim(micname)), ext))
        self%boxname = basename(fname_new_ext(trim(micname),'box'))
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
        self%maxdiam         = moldiam + moldiam * 0.111
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
        self%nx1             = self%ldim_shrink1(1) - self%ldim_box1(1)
        self%ny1             = self%ldim_shrink1(2) - self%ldim_box1(2)
        self%nx2             = self%ldim_shrink2(1) - self%ldim_box2(1)
        self%ny2             = self%ldim_shrink2(2) - self%ldim_box2(2)
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
        ! set ndev
        self%ndev = NDEV_DEFAULT
        if( present(ndev) ) self%ndev = ndev
        ! shrink micrograph
        call self%mic_shrink1%new(self%ldim_shrink1, SMPD_SHRINK1)
        call self%mic_shrink2%new(self%ldim_shrink2, SMPD_SHRINK2)
        call self%mic_shrink1%set_ft(.true.)
        call self%mic_shrink2%set_ft(.true.)
        call self%mic_raw%fft
        call self%mic_raw%clip(self%mic_shrink1)
        call self%mic_raw%clip(self%mic_shrink2)
        select case(trim(pcontrast))
            case('black')
                ! flip contrast (assuming black particle contrast on input)
                call self%mic_shrink1%mul(-1.)
                call self%mic_shrink2%mul(-1.)
            case('white')
                ! nothing to do
            case DEFAULT
                THROW_HARD('uknown pcontrast parameter, use (black|white)')
        end select
        ! back to real-space
        call self%mic_raw%ifft
        call self%mic_shrink1%ifft
        call self%mic_shrink2%ifft
        if( L_WRITE )then
            call self%mic_shrink1%write('mic_shrink1.mrc')
            call self%mic_shrink2%write('mic_shrink2.mrc')
            call self%imgau_shrink1%write('gau_shrink1.mrc')
            call self%imgau_shrink2%write('gau_shrink2.mrc')
        endif
        self%exists = .true.
    end subroutine new

    subroutine set_refs( self, imgs, mskdiam, hybrid )
        class(picker_utils), intent(inout) :: self
        class(image),        intent(inout) :: imgs(:)
        real,                intent(in)    :: mskdiam
        logical, optional,   intent(in)    :: hybrid
        type(image) :: img_rot
        integer     :: ldim(3), iimg, nimgs, irot, nrots, cnt
        real        :: scale1, scale2, mskrad1, mskrad2, pixrad_shrink1, pixrad_shrink2, smpd, ang, rot
        self%hybrid = .false.
        if( present(hybrid) ) self%hybrid = hybrid
        smpd        = imgs(1)%get_smpd()
        ldim        = imgs(1)%get_ldim()
        if( ldim(3) /= 1 ) THROW_HARD('box references must be 2D')
        nimgs       = size(imgs)
        nrots       = nint(real(MAXNREFS) / real(nimgs))
        self%nrefs  = nimgs * nrots
        ! set shrunken logical dimensions of boxes
        scale1               = smpd / SMPD_SHRINK1
        scale2               = smpd / SMPD_SHRINK2
        self%ldim_box1(1)    = round2even(real(ldim(1)) * scale1)
        self%ldim_box1(2)    = round2even(real(ldim(2)) * scale1)
        self%ldim_box1(3)    = 1
        self%ldim_box2(1)    = round2even(real(ldim(1)) * scale2)
        self%ldim_box2(2)    = round2even(real(ldim(2)) * scale2)
        self%ldim_box2(3)    = 1
        ! set # pixels in x/y for both box sizes
        self%nx1             = self%ldim_shrink1(1) - self%ldim_box1(1)
        self%ny1             = self%ldim_shrink1(2) - self%ldim_box1(2)
        self%nx2             = self%ldim_shrink2(1) - self%ldim_box2(1)
        self%ny2             = self%ldim_shrink2(2) - self%ldim_box2(2)
        ! set Gaussians
        self%maxdiam         = mskdiam
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
        call self%imgau_shrink1%fft
        call self%imgau_shrink2%fft
        mskrad1 = (mskdiam / SMPD_SHRINK1) / 2.
        mskrad2 = (mskdiam / SMPD_SHRINK2) / 2.
        if( allocated(self%boxrefs) )then
            do iimg = 1,size(self%boxrefs,dim=1)
                call self%boxrefs(iimg,1)%kill
                call self%boxrefs(iimg,2)%kill
            end do
            deallocate(self%boxrefs)
        endif
        if( allocated(self%sxx_refs) ) deallocate(self%sxx_refs)
        allocate( self%boxrefs(self%nrefs,2), self%sxx_refs(self%nrefs,2) )
        call img_rot%new(ldim, smpd)
        ang = 360./real(nrots)
        cnt = 0
        do iimg = 1,nimgs
            rot = 0.
            do irot = 1,nrots
                cnt = cnt + 1
                call imgs(iimg)%rtsq(rot, 0., 0., img_rot)
                call img_rot%fft
                call self%boxrefs(cnt,1)%new(self%ldim_box1, SMPD_SHRINK1)
                call self%boxrefs(cnt,2)%new(self%ldim_box2, SMPD_SHRINK2)
                call self%boxrefs(cnt,1)%set_ft(.true.)
                call self%boxrefs(cnt,2)%set_ft(.true.)
                call img_rot%clip(self%boxrefs(cnt,1))
                call img_rot%clip(self%boxrefs(cnt,2))
                ! convolve with Gaussians
                call self%boxrefs(cnt,1)%mul(self%imgau_shrink1)
                call self%boxrefs(cnt,2)%mul(self%imgau_shrink2)
                ! back to real-space
                call self%boxrefs(cnt,1)%ifft
                call self%boxrefs(cnt,2)%ifft
                call self%boxrefs(cnt,1)%mask(mskrad1, 'hard')
                call self%boxrefs(cnt,2)%mask(mskrad2, 'hard')
                call self%boxrefs(cnt,1)%prenorm4real_corr(self%sxx_refs(cnt,1))
                call self%boxrefs(cnt,2)%prenorm4real_corr(self%sxx_refs(cnt,2))
                if( L_WRITE )then
                    call self%boxrefs(cnt,1)%write('boxrefs_shrink1.mrc', cnt)
                    call self%boxrefs(cnt,2)%write('boxrefs_shrink2.mrc', cnt)
                endif
                rot = rot + ang
            end do
        end do
        call self%imgau_shrink1%ifft
        call self%imgau_shrink2%ifft
        call img_rot%kill
        self%refpick = .true.
    end subroutine set_refs

    subroutine exec_picker( self, boxname_out, nptcls, dir_out )
        class(picker_utils),        intent(inout) :: self
        character(len=LONGSTRLEN),  intent(out)   :: boxname_out
        integer,                    intent(out)   :: nptcls
        character(len=*), optional, intent(in)    :: dir_out
        real    :: shrink2
        integer :: box
        if( present(dir_out) )then
            self%fbody   = trim(dir_out)//trim(self%fbody)
            self%boxname = trim(dir_out)//trim(self%boxname)
        endif
        call self%gauconv_mics
        call self%set_positions
        call self%analyze_boximgs1
        call self%center_filter
        call self%distance_filter
        call self%refine_positions
        call self%remove_outliers
        call self%peak_vs_nonpeak_stats
        ! # ptcls
        nptcls          = self%nboxes2
        ! bring back coordinates to original sampling
        shrink2         = SMPD_SHRINK2 / self%smpd_raw
        self%positions2 = nint(shrink2 * real(self%positions2))
        ! write coordinates
        box = find_larger_magic_box(nint(shrink2 * self%ldim_box2(1)))
        call self%extract_boximgs2(box)
        call write_boximgs(self%nboxes2, self%boximgs2, trim(self%fbody)//'_raw.mrcs')
        call write_boxfile(self%nboxes2, self%positions2, box, self%boxname)
        call make_relativepath(CWD_GLOB, self%boxname, boxname_out) ! returns absolute path
    end subroutine exec_picker

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
        ! if( L_WRITE )then
            ! call self%mic_shrink1%write('gauconv_mic_shrink1.mrc')
            ! call img%zero_and_unflag_ft
            ! call sauvola(self%mic_shrink1, OFFSET, img)
            ! call img%write('sauvola.mrc')
            ! call img%real_space_filter(OFFSET * 3, 'average')
            ! call img%write('sauvola_rfilt.mrc')
            ! call otsu_img(img)
            ! call img%write('sauvola_otsu.mrc')
        ! endif
        ! denoise mic_shrink2
        call img%new(self%ldim_shrink2, SMPD_SHRINK2)
        call img%gauimg2D(self%sig_shrink2,self%sig_shrink2)
        call img%fft
        call self%mic_shrink2%fft
        call self%mic_shrink2%mul(img)
        call self%mic_shrink2%ifft
        ! if( L_WRITE )then
        !     call self%mic_shrink2%write('gauconv_mic_shrink2.mrc')
        ! endif
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

    subroutine extract_boximgs2( self, box_in )
        class(picker_utils), target, intent(inout) :: self
        integer, optional,           intent(in)    :: box_in
        type(image), pointer :: mic_ptr => null()
        integer :: ibox, ldim(3), noutside
        logical :: outside
        real    :: smpd
        if( .not. allocated(self%positions2) ) THROW_HARD('positions need to be set before constructing boximgs2')
        if( allocated(self%boximgs2) )then
            do ibox = 1,self%nboxes2
                call self%boximgs2(ibox)%kill
            end do
            deallocate(self%boximgs2)
        endif
        allocate(self%boximgs2(self%nboxes2))
        if( present(box_in) )then
            ldim    =  [box_in,box_in,1]
            smpd    =  self%smpd_raw
            mic_ptr => self%mic_raw
            !$omp parallel do schedule(static) default(shared) private(ibox,noutside) proc_bind(close)
            do ibox = 1,self%nboxes2
                call self%boximgs2(ibox)%new(ldim, smpd)
                call mic_ptr%window(self%positions2(ibox,:), ldim(1), self%boximgs2(ibox), noutside)
            end do
            !$omp end parallel do
        else
            ldim    =  self%ldim_box2
            smpd    =  SMPD_SHRINK2
            mic_ptr => self%mic_shrink2
            !$omp parallel do schedule(static) default(shared) private(ibox,outside) proc_bind(close)
            do ibox = 1,self%nboxes2
                call self%boximgs2(ibox)%new(ldim, smpd)
                call mic_ptr%window_slim(self%positions2(ibox,:), ldim(1), self%boximgs2(ibox), outside)
            end do
            !$omp end parallel do
        endif
    end subroutine extract_boximgs2

    subroutine analyze_boximgs1( self )
        class(picker_utils), intent(inout) :: self
        integer, allocatable :: positions_tmp(:,:)
        real,    allocatable :: tmp(:)
        real        :: box_scores(self%nx_offset,self%ny_offset,1), t, scores(self%nrefs)
        logical     :: is_peak(self%nx_offset,self%ny_offset,1), outside
        integer     :: ioff, joff, npeaks, cnt, ithr, iref
        type(image) :: boximgs_heap(nthr_glob)
        ! calculate box_scores
        if( .not. allocated(self%positions1) ) THROW_HARD('positions1 need to be set')
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box1, SMPD_SHRINK1)
        end do
        if( self%refpick .and. .not. self%hybrid )then
            !$omp parallel do schedule(static) collapse(2) default(shared) private(ioff,joff,ithr,outside,iref,scores) proc_bind(close)
            do ioff = 1,self%nx_offset
                do joff = 1,self%ny_offset
                    ithr = omp_get_thread_num() + 1
                    call self%mic_shrink1%window_slim(self%positions1(self%inds_offset(ioff,joff),:),&
                    &self%ldim_box1(1), boximgs_heap(ithr), outside)
                    do iref = 1,self%nrefs
                        scores(iref) = self%boxrefs(iref,1)%real_corr_prenorm(boximgs_heap(ithr), self%sxx_refs(iref,1))
                    end do
                    box_scores(ioff,joff,1) = maxval(scores)
                end do
            end do
            !$omp end parallel do
        else
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
        endif
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
        call self%extract_boximgs1
        if( L_WRITE )then
            call write_boximgs(int(self%nboxes1), self%boximgs1, trim(self%fbody)//'_before_filters.mrcs')
        endif
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%kill
        end do
    end subroutine analyze_boximgs1

    subroutine center_filter( self )
        class(picker_utils), intent(inout) :: self
        integer, allocatable :: positions_tmp(:,:)
        real,    allocatable :: box_scores1_tmp(:)
        type(image) :: boximgs_heap(nthr_glob)
        integer     :: ibox, ithr, npeaks, cnt
        real        :: scores(self%nboxes1)
        logical     :: mask(self%nboxes1)
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box1, SMPD_SHRINK1)
        end do
        !$omp parallel do schedule(static) default(shared) private(ibox) proc_bind(close)
        do ibox = 1,self%nboxes1
            ithr  = omp_get_thread_num() + 1
            scores(ibox) = self%boximgs1(ibox)%box_cen_arg(boximgs_heap(ithr))
        end do
        !$omp end parallel do
        mask   = scores(:) <= real(OFFSET)
        npeaks = count(mask)
        write(logfhandle,'(a,1x,I5)') '# positions before center   filtering: ', self%nboxes1
        write(logfhandle,'(a,1x,I5)') '# positions after  center   filtering: ', npeaks
       ! update positions1 and box_scores1
        allocate(positions_tmp(npeaks,2), box_scores1_tmp(npeaks))
        positions_tmp   = 0
        box_scores1_tmp = 0.
        cnt             = 0
        do ibox = 1,self%nboxes1
            if( mask(ibox) )then
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
    end subroutine center_filter

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
    end subroutine distance_filter

    subroutine refine_positions( self )
        class(picker_utils), intent(inout) :: self
        integer     :: ibox, xrange(2), yrange(2), xind, yind, ithr, old_pos(2), iref
        real        :: box_score, box_score_trial, factor, rpos(2), scores(self%nrefs)
        logical     :: outside
        type(image) :: boximgs_heap(nthr_glob)
        if( .not. allocated(self%positions1) ) THROW_HARD('positions1 need to be set')
        self%nboxes2 = self%nboxes1
        allocate(self%positions2(self%nboxes2,2), source=nint((SMPD_SHRINK1/SMPD_SHRINK2) * real(self%positions1)))
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box2, SMPD_SHRINK2)
        end do
        factor = real(OFFSET) * (SMPD_SHRINK1 / SMPD_SHRINK2)
        if( self%refpick )then
            !$omp parallel do schedule(static) default(shared) proc_bind(close)&
            !$omp private(ibox,rpos,xrange,yrange,box_score,xind,yind,ithr,outside,iref,scores,box_score_trial)
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
                        do iref = 1,self%nrefs
                            scores(iref) = self%boxrefs(iref,2)%real_corr_prenorm(boximgs_heap(ithr), self%sxx_refs(iref,2))
                        end do
                        box_score_trial = maxval(scores)
                        if( box_score_trial > box_score )then
                            self%positions2(ibox,:) = [xind,yind]
                            box_score = box_score_trial
                        endif
                    end do
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do schedule(static) default(shared) proc_bind(close)&
            !$omp private(ibox,rpos,xrange,yrange,box_score,xind,yind,ithr,outside,box_score_trial)
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
        endif
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%kill
        end do
    end subroutine refine_positions

    subroutine remove_outliers( self )
        class(picker_utils), intent(inout) :: self
        integer, allocatable :: positions_tmp(:,:)
        type(image) :: boximgs_heap(nthr_glob)
        integer     :: ibox, ithr, npeaks, cnt
        logical     :: outside
        real        :: factor, loc_sdevs(self%nboxes2), avg, sdev, t
        if( .not. allocated(self%positions2) ) THROW_HARD('positions2 need to be set')
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box2, SMPD_SHRINK2)
        end do
        factor = real(OFFSET) * (SMPD_SHRINK1 / SMPD_SHRINK2)
        !$omp parallel do schedule(static) default(shared) proc_bind(close) private(ibox,ithr,outside)
        do ibox = 1,self%nboxes2
            ithr            = omp_get_thread_num() + 1
            call self%mic_shrink2%window_slim([self%positions2(ibox,1),self%positions2(ibox,2)], self%ldim_box2(1), boximgs_heap(ithr), outside)
            loc_sdevs(ibox) = boximgs_heap(ithr)%avg_loc_sdev(nint(factor))
        end do
        !$omp end parallel do
        call avg_sdev(loc_sdevs, avg, sdev)
        ! write(logfhandle,'(a,1x,I5)') '# positions after 1.0 sigma outlier removal: ', count(loc_sdevs < avg + 1.0 * sdev)
        ! write(logfhandle,'(a,1x,I5)') '# positions after 1.5 sigma outlier removal: ', count(loc_sdevs < avg + 1.5 * sdev)
        ! write(logfhandle,'(a,1x,I5)') '# positions after 2.0 sigma outlier removal: ', count(loc_sdevs < avg + 2.0 * sdev)
        ! write(logfhandle,'(a,1x,I5)') '# positions after 2.5 sigma outlier removal: ', count(loc_sdevs < avg + 2.5 * sdev)
        ! write(logfhandle,'(a,1x,I5)') '# positions after 3.0 sigma outlier removal: ', count(loc_sdevs < avg + 3.0 * sdev)
        t = avg + self%ndev * sdev
        npeaks = count(loc_sdevs < t)
        write(logfhandle,'(a,1x,I5)') '# positions after  outlier    removal: ', npeaks
        ! update positions2
        allocate(positions_tmp(npeaks,2), source=0)
        cnt = 0
        do ibox = 1,self%nboxes2
            if( loc_sdevs(ibox) < t )then
                cnt = cnt + 1
                positions_tmp(cnt,:) = self%positions2(ibox,:)
            endif
        end do
        deallocate(self%positions2)
        self%nboxes2 = npeaks
        allocate(self%positions2(self%nboxes2,2), source=positions_tmp)
        deallocate(positions_tmp)
        if( L_WRITE )then
            call self%extract_boximgs2
            call write_boximgs(int(self%nboxes2), self%boximgs2, trim(self%fbody)//'_after_filters.mrcs')
        endif
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%kill
        end do
    end subroutine remove_outliers

    subroutine peak_vs_nonpeak_stats( self )
        class(picker_utils), intent(inout) :: self
        integer, allocatable :: positions_backgr(:,:)
        real,    allocatable :: scores_nonpeak(:)
        type(image) :: boximgs_heap(nthr_glob)
        real        :: scores_peak(self%nboxes2), factor, smd, ksstat, prob
        real        :: scores(self%nrefs), pixrad_shrink2, rpos(2), a_peak, a_nonpeak, s_peak, s_nonpeak
        logical     :: mask_backgr(0:self%nx2,0:self%ny2), outside
        integer     :: ibox, jbox, off_here, nbackgr, ithr, iref, xrange(2), yrange(2), xind, yind
        if( .not. allocated(self%positions2) ) THROW_HARD('positions2 need to be set')
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box2, SMPD_SHRINK2)
        end do
        ! prepare background mask
        mask_backgr    = .true. 
        pixrad_shrink2 = (self%maxdiam / 2.) / SMPD_SHRINK2
        do ibox = 1,self%nboxes2
            rpos      = real(self%positions2(ibox,:))
            xrange(1) = max(0,        nint(rpos(1) - pixrad_shrink2))
            xrange(2) = min(self%nx2, nint(rpos(1) + pixrad_shrink2))
            yrange(1) = max(0,        nint(rpos(2) - pixrad_shrink2))
            yrange(2) = min(self%ny2, nint(rpos(2) + pixrad_shrink2))
            do xind = xrange(1),xrange(2)
                do yind = yrange(1),yrange(2)
                    mask_backgr(xind,yind) = .false.
                end do
            end do
        end do
        off_here = nint(real(OFFSET) * (SMPD_SHRINK1 / SMPD_SHRINK2))
        nbackgr  = 0
        do ibox = 0,self%nx2,off_here
            do jbox = 0,self%ny2,off_here
                if( mask_backgr(ibox,jbox) )then
                    nbackgr = nbackgr + 1
                endif
            end do
        end do
        allocate(positions_backgr(nbackgr,2), source=0)
        nbackgr = 0
        do ibox = 0,self%nx2,off_here
            do jbox = 0,self%ny2,off_here
                if( mask_backgr(ibox,jbox) )then
                    nbackgr = nbackgr + 1
                    positions_backgr(nbackgr,:) = [ibox,jbox]
                endif
            end do
        end do
        ! extract background info
        allocate(scores_nonpeak(nbackgr), source=0.)
        factor = real(OFFSET) * (SMPD_SHRINK1 / SMPD_SHRINK2)
        if( self%refpick )then
            !$omp parallel do schedule(static) default(shared) proc_bind(close) private(ibox,ithr,outside,iref,scores)
            do ibox = 1,nbackgr
                ithr = omp_get_thread_num() + 1
                call self%mic_shrink2%window_slim([positions_backgr(ibox,1),positions_backgr(ibox,2)], self%ldim_box2(1), boximgs_heap(ithr), outside)
                do iref = 1,self%nrefs
                    scores(iref) = self%boxrefs(iref,2)%real_corr_prenorm(boximgs_heap(ithr), self%sxx_refs(iref,2))
                end do
                scores_nonpeak(ibox) = maxval(scores)
            end do
            !$omp end parallel do
        else
            !$omp parallel do schedule(static) default(shared) proc_bind(close) private(ibox,ithr,outside)
            do ibox = 1,nbackgr
                ithr = omp_get_thread_num() + 1
                call self%mic_shrink2%window_slim([positions_backgr(ibox,1),positions_backgr(ibox,2)], self%ldim_box2(1), boximgs_heap(ithr), outside)
                scores_nonpeak(ibox)    = self%imgau_shrink2%real_corr_prenorm(boximgs_heap(ithr), self%sxx_shrink2)
            end do
            !$omp end parallel do
        endif
        if( self%refpick )then
            !$omp parallel do schedule(static) default(shared) proc_bind(close) private(ibox,ithr,outside,iref,scores)
            do ibox = 1,self%nboxes2
                ithr = omp_get_thread_num() + 1
                call self%mic_shrink2%window_slim([self%positions2(ibox,1),self%positions2(ibox,2)], self%ldim_box2(1), boximgs_heap(ithr), outside)
                do iref = 1,self%nrefs
                    scores(iref) = self%boxrefs(iref,2)%real_corr_prenorm(boximgs_heap(ithr), self%sxx_refs(iref,2))
                end do
                scores_peak(ibox) = maxval(scores)
            end do
            !$omp end parallel do
        else
            !$omp parallel do schedule(static) default(shared) proc_bind(close) private(ibox,ithr,outside)
            do ibox = 1,self%nboxes2
                ithr = omp_get_thread_num() + 1
                call self%mic_shrink2%window_slim([self%positions2(ibox,1),self%positions2(ibox,2)], self%ldim_box2(1), boximgs_heap(ithr), outside)
                scores_peak(ibox) = self%imgau_shrink2%real_corr_prenorm(boximgs_heap(ithr), self%sxx_shrink2)
            end do
            !$omp end parallel do
        endif
        call avg_sdev(scores_peak,    a_peak,    s_peak)
        call avg_sdev(scores_nonpeak, a_nonpeak, s_nonpeak)
        smd = std_mean_diff(a_peak, a_nonpeak, s_peak, s_nonpeak)
        call kstwo(scores_peak, self%nboxes2, scores_nonpeak, nbackgr, ksstat, prob)
        write(logfhandle,'(a,1x,f4.2)') 'SMD           = ', smd
        write(logfhandle,'(a,1x,f4.2)') 'K-S statistic = ', ksstat
        write(logfhandle,'(a,1x,f4.2)') 'P             = ', prob
        if( smd < 0.2 .and. prob > 0.5 ) write(logfhandle,'(a)') 'peak and non-peak distributions of fom:s are similar' 
    end subroutine peak_vs_nonpeak_stats

    subroutine kill( self )
        class(picker_utils), intent(inout) :: self
        integer :: i
        if( self%exists )then
            call self%mic_raw%kill
            call self%mic_shrink1%kill
            call self%mic_shrink2%kill
            call self%imgau_shrink1%kill
            call self%imgau_shrink2%kill
            call destruct_imarr(self%boximgs1)
            call destruct_imarr(self%boximgs2)
            if( allocated(self%boxrefs) )then
                do i = 1,size(self%boxrefs, dim=1)
                    call self%boxrefs(i,1)%kill
                    call self%boxrefs(i,2)%kill
                end do
                deallocate(self%boxrefs)
            endif
            if( allocated(self%positions1)  ) deallocate(self%positions1)
            if( allocated(self%positions2)  ) deallocate(self%positions2)
            if( allocated(self%inds_offset) ) deallocate(self%inds_offset)
            if( allocated(self%box_scores1) ) deallocate(self%box_scores1)
            if( allocated(self%sxx_refs)    ) deallocate(self%sxx_refs)
            if( allocated(self%fbody)       ) deallocate(self%fbody)
        endif

        contains

            subroutine destruct_imarr( imarr )
                type(image), allocatable, intent(inout) :: imarr(:)
                integer :: i
                if( allocated(imarr) )then
                    do i = 1,size(imarr)
                        call imarr(i)%kill
                    end do
                    deallocate(imarr)
                endif
            end subroutine destruct_imarr

    end subroutine kill

    ! utilities

    subroutine write_boximgs( n, boximgs, fname, mask )
        integer,           intent(in)    :: n
        class(image),      intent(inout) :: boximgs(n)
        character(len=*),  intent(in)    :: fname
        logical, optional, intent(in)    :: mask(n)
        integer :: ibox, cnt

        print *, 'boximgs file: ', trim(fname)

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

    subroutine write_boxfile( n, coordinates, box, fname )
        integer,          intent(in) :: n, coordinates(n,2), box
        character(len=*), intent(in) :: fname
        integer :: funit, ibox, iostat
        call fopen(funit, status='REPLACE', action='WRITE', file=trim(adjustl(fname)), iostat=iostat)
        call fileiochk('simple_picker_utils; write_boxfile ', iostat)
        do ibox = 1,n
            write(funit,'(I7,I7,I7,I7,I7)') coordinates(ibox,1), coordinates(ibox,2), box, box, -3
        end do
        call fclose(funit)
    end subroutine write_boxfile

end module simple_picker_utils
