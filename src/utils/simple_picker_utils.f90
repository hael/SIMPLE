module simple_picker_utils
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
implicit none

public :: picker_utils
private
#include "simple_local_flags.inc"

real,    parameter :: SMPD_SHRINK1 = 4.0, SMPD_SHRINK2 = 2.0, GAUSIG = 5., BOX_EXP_FAC = 0.111
integer, parameter :: OFFSET       = 3,   OFFSET_UB    = 2 * OFFSET
integer, parameter :: MAXNREFS     = 100
real,    parameter :: DIST_THRES1  = real(OFFSET), DIST_THRES2  = real(4*OFFSET)
real,    parameter :: NDEV_DEFAULT = 2.5
logical, parameter :: L_WRITE      = .false.
logical, parameter :: L_DEBUG      = .false.

type picker_utils
    private
    integer                       :: ldim_raw(3) = 0 , ldim_shrink1(3) = 0 , ldim_shrink2(3) = 0
    integer                       :: ldim_box(3) = 0 , ldim_box1(3)    = 0 , ldim_box2(3)    = 0
    real                          :: smpd_raw    = 0.
    real                          :: maxdiam = 0., sig_shrink1 = 0., sig_shrink2 = 0., sxx_shrink1 = 0.
    real                          :: sxx_shrink2 = 0., ndev = 0.
    integer                       :: nboxes1 = 0, nboxes2 = 0, nboxes_ub = 0, nrefs = 1
    integer                       :: nx1 = 0, ny1 = 0, nx2 = 0, ny2 = 0, nx_offset  = 0, ny_offset = 0
    type(image)                   :: mic_raw, mic_shrink1, mic_shrink2, imgau_shrink1, imgau_shrink2
    type(image),      allocatable :: boximgs1(:), boximgs2(:), boxrefs(:,:)
    integer,          allocatable :: positions1(:,:), positions2(:,:), inds_offset(:,:)
    real,             allocatable :: box_scores1(:), box_scores2(:)
    logical,          allocatable :: l_err_refs(:,:), l_mic_mask1(:,:)
    character(len=LONGSTRLEN)     :: boxname
    character(len=:), allocatable :: fbody
    logical                       :: refpick = .false., l_roi = .false.
    logical                       :: exists  = .false.
  contains
    procedure          :: new
    procedure          :: set_refs
    procedure          :: exec_picker
    procedure, private :: set_positions_1, set_positions_2
    generic            :: set_positions => set_positions_1, set_positions_2
    procedure, private :: set_pos_priv
    procedure, private :: gauconv_mics
    procedure, private :: extract_boximgs1
    procedure, private :: extract_boximgs2
    procedure, private :: analyze_boximgs1
    procedure, private :: center_filter
    procedure, private :: distance_filter
    procedure          :: refine_positions
    procedure          :: remove_outliers
    procedure          :: peak_vs_nonpeak_stats
    procedure, private :: flag_ice, flag_amorphous_carbon
    procedure          :: kill
end type picker_utils

contains

    subroutine new( self, micname, pcontrast, smpd, moldiam, ndev )
        class(picker_utils),  intent(inout) :: self
        character(len=*),     intent(in)    :: micname, pcontrast
        real,                 intent(in)    :: smpd    !< sampling distance in A
        real,                 intent(in)    :: moldiam !< maximum diameter in A
        real, optional,       intent(in)    :: ndev    !< # std devs for outlier detection
        type(image)                   :: mic_pad
        character(len=:), allocatable :: ext
        real    :: scale1, scale2, pixrad_shrink1, pixrad_shrink2, hpfreq
        integer :: ldim_pd(3), nframes
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
        self%maxdiam         = moldiam + moldiam * BOX_EXP_FAC
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
        ! init mask
        self%l_roi = trim(params_glob%pick_roi).eq.'yes'
        ! high-pass micrograph
        if( self%l_roi )then
            ldim_pd(1:2) = find_larger_magic_box(self%ldim_raw(1:2)+1)
            if( minval(ldim_pd(1:2)-self%ldim_raw(1:2)) < 16 )then
                ldim_pd(1:2) = find_larger_magic_box(ldim_pd(1:2)+1)
            endif
            ldim_pd(3) = 1
            call mic_pad%new(ldim_pd, self%smpd_raw)
            call self%mic_raw%pad_mirr(mic_pad)
            call mic_pad%fft
            hpfreq = real(minval(ldim_pd(1:2)))*SMPD_SHRINK1/16.
            call mic_pad%bp(hpfreq,0.)
            call mic_pad%ifft
            call mic_pad%clip(self%mic_raw)
            call mic_pad%kill
        endif
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
        allocate(self%l_mic_mask1(self%ldim_shrink1(1),self%ldim_shrink1(2)),source=.true.)
        if( L_WRITE )then
            call self%mic_shrink1%write('mic_shrink1.mrc')
            call self%mic_shrink2%write('mic_shrink2.mrc')
            call self%imgau_shrink1%write('gau_shrink1.mrc')
            call self%imgau_shrink2%write('gau_shrink2.mrc')
        endif
        self%exists = .true.
    end subroutine new

    subroutine set_refs( self, imgs, mskdiam )
        class(picker_utils), intent(inout) :: self
        class(image),        intent(inout) :: imgs(:)
        real,                intent(in)    :: mskdiam
        type(image) :: img_rot
        integer     :: ldim(3), iimg, nimgs, irot, nrots, cnt
        real        :: scale1, scale2, mskrad1, mskrad2, pixrad_shrink1, pixrad_shrink2, smpd, ang, rot
        smpd        = imgs(1)%get_smpd()
        ldim        = imgs(1)%get_ldim()
        if( ldim(3) /= 1 ) THROW_HARD('box references must be 2D')
        nimgs       = size(imgs)
        nrots       = nint(real(MAXNREFS) / real(nimgs))
        self%nrefs  = nimgs * nrots
        ! set shrunken logical dimensions of boxes
        scale1            = smpd / SMPD_SHRINK1
        scale2            = smpd / SMPD_SHRINK2
        self%ldim_box1(1) = round2even(real(ldim(1)) * scale1)
        self%ldim_box1(2) = round2even(real(ldim(2)) * scale1)
        self%ldim_box1(3) = 1
        self%ldim_box2(1) = round2even(real(ldim(1)) * scale2)
        self%ldim_box2(2) = round2even(real(ldim(2)) * scale2)
        self%ldim_box2(3) = 1
        ! set # pixels in x/y for both box sizes
        self%nx1          = self%ldim_shrink1(1) - self%ldim_box1(1)
        self%ny1          = self%ldim_shrink1(2) - self%ldim_box1(2)
        self%nx2          = self%ldim_shrink2(1) - self%ldim_box2(1)
        self%ny2          = self%ldim_shrink2(2) - self%ldim_box2(2)
        ! set Gaussians
        self%maxdiam      = mskdiam
        pixrad_shrink1    = (self%maxdiam / 2.) / SMPD_SHRINK1
        pixrad_shrink2    = (self%maxdiam / 2.) / SMPD_SHRINK2
        self%sig_shrink1  = pixrad_shrink1 / GAUSIG
        self%sig_shrink2  = pixrad_shrink2 / GAUSIG
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
        if( allocated(self%l_err_refs) ) deallocate(self%l_err_refs)
        allocate(self%l_err_refs(self%nrefs,2),self%boxrefs(self%nrefs,2))
        call img_rot%new(ldim, smpd)
        ang = 360./real(nrots)
        cnt = 0
        do iimg = 1,nimgs
            rot = 0.
            do irot = 1,nrots
                cnt = cnt + 1
                call imgs(iimg)%rtsq(rot, 0., 0., img_rot)
                call img_rot%fft
                call self%boxrefs(cnt,1)%new(self%ldim_box1, SMPD_SHRINK1, wthreads=.false.)
                call self%boxrefs(cnt,2)%new(self%ldim_box2, SMPD_SHRINK2, wthreads=.false.)
                call self%boxrefs(cnt,1)%set_ft(.true.)
                call self%boxrefs(cnt,2)%set_ft(.true.)
                call img_rot%clip(self%boxrefs(cnt,1))
                call img_rot%clip(self%boxrefs(cnt,2))
                rot = rot + ang
            end do
        end do
        !$omp parallel do schedule(static) default(shared) private(cnt) proc_bind(close)
        do cnt = 1,self%nrefs
            ! convolve with Gaussians
            call self%boxrefs(cnt,1)%mul(self%imgau_shrink1)
            call self%boxrefs(cnt,2)%mul(self%imgau_shrink2)
            ! back to real-space
            call self%boxrefs(cnt,1)%ifft
            call self%boxrefs(cnt,2)%ifft
            call self%boxrefs(cnt,1)%mask(mskrad1, 'hard')
            call self%boxrefs(cnt,2)%mask(mskrad2, 'hard')
            call self%boxrefs(cnt,1)%prenorm4real_corr(self%l_err_refs(cnt,1))
            call self%boxrefs(cnt,2)%prenorm4real_corr(self%l_err_refs(cnt,2))
        end do
        !$omp end parallel do
        if( L_WRITE )then
            do cnt = 1,self%nrefs
                call self%boxrefs(cnt,1)%write('boxrefs_shrink1.mrc', cnt)
                call self%boxrefs(cnt,2)%write('boxrefs_shrink2.mrc', cnt)
            enddo
        endif
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
        type(image) :: mic,mic_pick,collage
        real    :: shrink2
        integer :: box,i,j,k,ii,jj,pos(2)
        if( present(dir_out) )then
            self%fbody   = trim(dir_out)//trim(self%fbody)
            self%boxname = trim(dir_out)//trim(self%boxname)
        endif
        if( self%l_roi )then
            call self%flag_ice
            call self%flag_amorphous_carbon
            if( count(self%l_mic_mask1) < nint(0.02*product(self%ldim_shrink1)) )then
                nptcls      = 0
                boxname_out = ''
                return
            endif
        endif
        call self%gauconv_mics
        call self%set_positions
        call self%analyze_boximgs1
        if( self%nboxes1 > 0 )then
            call self%center_filter
            call self%distance_filter(self%nboxes1, self%box_scores1, self%positions1, DIST_THRES1)
            call self%refine_positions
            if( self%l_roi )then
                call self%distance_filter(self%nboxes2, self%box_scores2, self%positions2, DIST_THRES2)
            endif
            call self%remove_outliers
            call self%peak_vs_nonpeak_stats
            ! # ptcls
            nptcls          = self%nboxes2
            ! bring back coordinates to original sampling
            shrink2         = SMPD_SHRINK2 / self%smpd_raw
            self%positions2 = nint(shrink2 * real(self%positions2))
            ! debug
            if( self%l_roi .and. L_DEBUG )then
                call mic%copy(self%mic_shrink1)
                call mic%norm
                call mic_pick%copy(mic)
                do i = 1,nptcls
                    pos = nint(real(self%positions2(i,:))*self%smpd_raw/SMPD_SHRINK1)
                    pos = pos + self%ldim_box1(1)/2 + 1
                    do j = pos(1)-3,pos(1)+3
                        ii = min(self%ldim_shrink1(1),max(1,j)) 
                        do k = pos(2)-3,pos(2)+3
                            jj = min(self%ldim_shrink1(2),max(1,k)) 
                            call mic_pick%set([ii,jj,1], -1.)
                        enddo
                    enddo
                enddo
                do i = 1,self%ldim_shrink1(1)
                    do j = 1,self%ldim_shrink1(2)
                        if( self%l_mic_mask1(i,j) )cycle
                        call mic%set([i,j,1], 0.)
                    enddo
                enddo
                call mic%collage(mic_pick,collage)
                call collage%write_jpg(trim(self%fbody)//'_collage.jpg')
                call mic%kill
                call mic_pick%kill
                call collage%kill
            endif
            ! write coordinates
            box = find_larger_magic_box(nint(shrink2 * self%ldim_box2(1)))
            if( L_WRITE )then
                call self%extract_boximgs2(box)
                call write_boximgs(self%nboxes2, self%boximgs2, trim(self%fbody)//'_raw.mrcs')
            endif
            call write_boxfile(self%nboxes2, self%positions2, box, self%boxname)
            call make_relativepath(CWD_GLOB, self%boxname, boxname_out) ! returns absolute path
        else
            ! no particles found
            nptcls      = 0
            boxname_out = ''
        endif
    end subroutine exec_picker

    subroutine flag_ice( self )
        class(picker_utils), intent(inout) :: self
        real, parameter      :: THRESHOLD = 5.
        real, parameter      :: SMPD_ICE  = ICE_BAND1/2. - 0.15
        integer, parameter   :: BOX = 128
        type(image)          :: boximgs_heap(nthr_glob), img
        real,    allocatable :: scores(:,:)
        integer, allocatable :: counts(:,:)
        real        :: score, scale, radius
        integer     :: ldim(3), i,j,ithr,ii,jj
        logical     :: outside
        if( self%smpd_raw > SMPD_ICE ) return
        scale   = self%smpd_raw / SMPD_ICE
        ldim(1) = round2even(real(self%ldim_raw(1)) * scale)
        ldim(2) = round2even(real(self%ldim_raw(2)) * scale)
        ldim(3) = 1
        radius  = real(BOX)/2.- COSMSKHALFWIDTH
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new([BOX,BOX,1], SMPD_ICE, wthreads=.false.)
        end do
        call img%new(ldim,SMPD_ICE)
        call img%set_ft(.true.)
        call self%mic_raw%fft
        call self%mic_raw%clip(img)
        call self%mic_raw%ifft
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
                scores(i:i+BOX-1,j:j+BOX-1) = scores(i:i+BOX-1,j:j+BOX-1) + score
                counts(i:i+BOX-1,j:j+BOX-1) = counts(i:i+BOX-1,j:j+BOX-1) + 1
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
            scores(i:i+BOX-1,j:j+BOX-1) = scores(i:i+BOX-1,j:j+BOX-1) + score
            counts(i:i+BOX-1,j:j+BOX-1) = counts(i:i+BOX-1,j:j+BOX-1) + 1
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
            scores(i:i+BOX-1,j:j+BOX-1) = scores(i:i+BOX-1,j:j+BOX-1) + score
            counts(i:i+BOX-1,j:j+BOX-1) = counts(i:i+BOX-1,j:j+BOX-1) + 1
        enddo
        !$omp end parallel do
        where( counts > 0 ) scores = scores / real(counts)
        scale = real(ldim(1)) / real(self%ldim_shrink1(1))
        do i = 1,self%ldim_shrink1(1)
            ii = min(ldim(1), max(1, nint(scale*real(i))))
            do j = 1,self%ldim_shrink1(2)
                jj = min(ldim(2), max(1, nint(scale*real(j))))
                self%l_mic_mask1(i,j) = self%l_mic_mask1(i,j) .and. (scores(ii,jj) < THRESHOLD)
            enddo
        enddo
        if( .not.all(self%l_mic_mask1(:,:)) )then
            ! do i = 1,ldim(1)
            !     do j = 1,ldim(2)
            !         call img%set([i,j,1], scores(i,j))
            !     end do
            ! end do
            ! call img%write(trim(self%fbody)//'_ice.mrc')
            print *,trim(self%fbody)//' % ice water: ',&
                &100.-100.*count(self%l_mic_mask1(:,:))/real(product(self%ldim_shrink1))
        endif
        call img%kill
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%kill
        end do
    end subroutine flag_ice

    subroutine flag_amorphous_carbon( self )
        use simple_histogram, only: histogram
        class(picker_utils), intent(inout) :: self
        integer, parameter :: K = 5
        integer, parameter :: BOX           = 32 ! multiple of 4
        integer, parameter :: NBINS         = 64
        real,    parameter :: TVD_THRESHOLD = 0.2
        real,    parameter :: MIN_TVD_DIFF  = 0.05
        real,    parameter :: MIC_LP        = 15.
        type(histogram),  allocatable :: hists(:,:)
        type(ran_tabu)                :: rt
        type(histogram)               :: khists(K)
        type(image)                   :: mic, patches(nthr_glob)
        logical,          allocatable :: final_mask(:,:)
        character(len=:), allocatable :: string
        integer :: dims(3), pad, i
        logical :: found, empty
        found = .false.
        dims   = self%mic_shrink1%get_ldim()
        string = trim(self%fbody)
        allocate(final_mask(dims(1),dims(2)),source=.true.)
        call mic%copy(self%mic_shrink1)
        ! call gradients_variance_rejection( found )
        call clustering_rejection( found )
        if( found )then
            self%l_mic_mask1 = self%l_mic_mask1 .and. final_mask
            ! accounting for center of picking boxes
            pad = nint(real(self%ldim_box1(1))/2.)
            final_mask = self%l_mic_mask1
            self%l_mic_mask1 = .false.
            self%l_mic_mask1(1:dims(1)-pad,1:dims(2)-pad) = final_mask(pad+1:,pad+1:)
        endif
        print *,string//' % amorphous carbon: ',100.-100.*count(final_mask(:,:))/real(product(dims))
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
                real     :: bin_vars(K), tmpvec(K), tvd_inter(K-1)
                real     :: kmeans_score, mean, tvd, min_tvd
                logical  :: mask(dims(1),dims(2)), bin_mask(K), outside, cluster_mask(K-1)
                integer  :: b,twob,bon2,bon4,i,j,m,n,ilb,jlb,iub,jub,nx,ny,nxy
                integer  :: bin,repeat,ii,jj,carbon,ithr
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
                        if(.not.patches(ithr)%exists() ) call patches(ithr)%new([twob,twob,1], SMPD_SHRINK1, wthreads=.false.)
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
                if( L_DEBUG )then
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
                        call tmpimg2%write(trim(string)//'_clusters.mrc',i)
                    enddo
                    call mic%write(trim(string)//'_clusters.mrc',K+1)
                    call tmpimg2%kill
                endif
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
                    if( L_DEBUG )then
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
                integer,         intent(in)  :: n, nx, ny, K
                real(dp),        intent(out) :: overall_score
                integer,         intent(out) :: labels(n)
                integer, parameter :: nits = 100
                integer  :: pops(K), new_labels(n)
                real(dp) :: scores(K), kscores(K), prev_score, score
                integer  :: i,j,c,cl,it,nk
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
                                c = (i-1)*nx+j
                                if( labels(c) /= cl )cycle
                                call khists(cl)%add(hists(i,j))
                            enddo
                        enddo
                    enddo
                    !$omp end do
                    !$omp do collapse(2) reduction(+:scores)
                    do i = 1,nx
                        do j = 1,ny
                            c = (i-1)*nx+j
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

    subroutine set_positions_1( self )
        class(picker_utils), intent(inout) :: self
        call self%set_pos_priv
    end subroutine set_positions_1

    subroutine set_positions_2( self, box_raw, box12 )
        class(picker_utils), intent(inout) :: self
        integer,             intent(in)    :: box_raw, box12(2)
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
                if( self%l_mic_mask1(xind+1,yind+1) ) self%nboxes1 = self%nboxes1 + 1
            end do
        end do
        ! count # boxes, upper bound
        self%nboxes_ub = 0
        do xind = 0,self%nx1,OFFSET_UB
            do yind = 0,self%ny1,OFFSET_UB
                if( self%l_mic_mask1(xind+1,yind+1) ) self%nboxes_ub = self%nboxes_ub + 1
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
                self%ny_offset = self%ny_offset + 1
                if( self%l_mic_mask1(xind+1,yind+1) )then
                    self%nboxes1 = self%nboxes1 + 1
                    self%positions1(self%nboxes1,:) = [xind,yind]
                    self%inds_offset(self%nx_offset,self%ny_offset) = self%nboxes1
                endif
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
        integer :: pos(2), ibox
        logical :: outside
        if( .not. allocated(self%positions1) ) THROW_HARD('positions need to be set before constructing boximgs1')
        if( allocated(self%boximgs1) )then
            do ibox = 1,self%nboxes1
                call self%boximgs1(ibox)%kill
            end do
            deallocate(self%boximgs1)
        endif
        allocate(self%boximgs1(self%nboxes1))
        !$omp parallel do schedule(static) default(shared) private(ibox,outside,pos) proc_bind(close)
        do ibox = 1,self%nboxes1
            call self%boximgs1(ibox)%new(self%ldim_box1, SMPD_SHRINK1)
            pos = self%positions1(ibox,:)
            call self%mic_shrink1%window_slim(pos, self%ldim_box1(1), self%boximgs1(ibox), outside)
        end do
        !$omp end parallel do
    end subroutine extract_boximgs1

    subroutine extract_boximgs2( self, box_in )
        class(picker_utils), target, intent(inout) :: self
        integer, optional,           intent(in)    :: box_in
        type(image), pointer :: mic_ptr => null()
        integer :: pos(2), ldim(3), ibox, noutside
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
            !$omp parallel do schedule(static) default(shared) private(pos,ibox,noutside) proc_bind(close)
            do ibox = 1,self%nboxes2
                call self%boximgs2(ibox)%new(ldim, smpd)
                pos = self%positions2(ibox,:)
                call mic_ptr%window(pos, ldim(1), self%boximgs2(ibox), noutside)
            end do
            !$omp end parallel do
        else
            ldim    =  self%ldim_box2
            smpd    =  SMPD_SHRINK2
            mic_ptr => self%mic_shrink2
            !$omp parallel do schedule(static) default(shared) private(pos,ibox,outside) proc_bind(close)
            do ibox = 1,self%nboxes2
                call self%boximgs2(ibox)%new(ldim, smpd)
                pos = self%positions2(ibox,:)
                call mic_ptr%window_slim(pos, ldim(1), self%boximgs2(ibox), outside)
            end do
            !$omp end parallel do
        endif
    end subroutine extract_boximgs2

    subroutine analyze_boximgs1( self )
        class(picker_utils), intent(inout) :: self
        integer, allocatable :: positions_tmp(:,:)
        real,    allocatable :: tmp(:)
        real        :: box_scores(self%nx_offset,self%ny_offset,1), t, scores(self%nrefs)
        logical     :: is_peak(self%nx_offset,self%ny_offset,1), outside, l_err
        integer     :: pos(2), ioff, joff, npeaks, cnt, ithr, iref
        type(image) :: boximgs_heap(nthr_glob)
        ! calculate box_scores
        if( .not. allocated(self%positions1) ) THROW_HARD('positions1 need to be set')
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box1, SMPD_SHRINK1)
        end do
        if( self%refpick )then
            !$omp parallel do schedule(static) collapse(2) default(shared) private(ioff,joff,ithr,outside,iref,scores,pos,l_err) proc_bind(close)
            do ioff = 1,self%nx_offset
                do joff = 1,self%ny_offset
                    if( self%inds_offset(ioff,joff) == 0 )then
                        box_scores(ioff,joff,1) = -1.
                    else
                        pos = self%positions1(self%inds_offset(ioff,joff),:)
                        if( self%l_mic_mask1(pos(1)+1,pos(2)+1) )then
                            ithr = omp_get_thread_num() + 1
                            pos  = self%positions1(self%inds_offset(ioff,joff),:)
                            call self%mic_shrink1%window_slim(pos, self%ldim_box1(1), boximgs_heap(ithr), outside)
                            call boximgs_heap(ithr)%prenorm4real_corr(l_err)
                            if( l_err )then
                                box_scores(ioff,joff,1) = 0.
                            else
                                do iref = 1,self%nrefs
                                    if( self%l_err_refs(iref,1) )then
                                        scores(iref) = 0.
                                    else
                                        scores(iref) = self%boxrefs(iref,1)%real_corr_prenorm(boximgs_heap(ithr))
                                    endif
                                end do
                                box_scores(ioff,joff,1) = maxval(scores)
                            endif
                        else
                            box_scores(ioff,joff,1) = -1.
                        endif
                    endif
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do schedule(static) collapse(2) default(shared) private(ioff,joff,ithr,outside,pos) proc_bind(close)
            do ioff = 1,self%nx_offset
                do joff = 1,self%ny_offset
                    if( self%inds_offset(ioff,joff) == 0 )then
                        box_scores(ioff,joff,1) = -1.
                    else
                        pos = self%positions1(self%inds_offset(ioff,joff),:)
                        if( self%l_mic_mask1(pos(1)+1,pos(2)+1) )then
                            ithr = omp_get_thread_num() + 1
                            call self%mic_shrink1%window_slim(pos, self%ldim_box1(1), boximgs_heap(ithr), outside)
                            box_scores(ioff,joff,1) = self%imgau_shrink1%real_corr_prenorm(boximgs_heap(ithr), self%sxx_shrink1)
                        else
                            box_scores(ioff,joff,1) = -1.
                        endif
                    endif
                end do
            end do
            !$omp end parallel do
        endif
        tmp = pack(box_scores, mask=(box_scores>-1.+TINY))
        call detect_peak_thres(size(tmp), int(self%nboxes_ub), tmp, t)
        deallocate(tmp)
        t = max(0.,t)
        is_peak = .false.
        do ioff = 1,self%nx_offset
            do joff = 1,self%ny_offset
                if( box_scores(ioff,joff,1) >= t ) is_peak(ioff,joff,1) = .true.
            end do
        end do
        npeaks = count(is_peak)
        write(logfhandle,'(a,1x,I5)') '# positions considered               : ', self%nboxes1
        self%nboxes1 = npeaks
        if( self%nboxes1 > 0 )then
            ! update positions1 and box_scores1
            if( allocated(self%box_scores1) ) deallocate(self%box_scores1)
            allocate(positions_tmp(self%nboxes1,2), self%box_scores1(self%nboxes1))
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
            allocate(self%positions1(self%nboxes1,2), source=positions_tmp)
            deallocate(positions_tmp)
            call self%extract_boximgs1
            if( L_WRITE )then
                call write_boximgs(int(self%nboxes1), self%boximgs1, trim(self%fbody)//'_before_filters.mrcs')
            endif
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

    subroutine refine_positions( self )
        class(picker_utils), intent(inout) :: self
        integer     :: ibox, xrange(2), yrange(2), xind, yind, ithr, iref
        real        :: box_score, box_score_trial, factor, rpos(2), scores(self%nrefs)
        logical     :: outside, l_err
        type(image) :: boximgs_heap(nthr_glob)
        if( .not. allocated(self%positions1) ) THROW_HARD('positions1 need to be set')
        self%nboxes2 = self%nboxes1
        allocate(self%positions2(self%nboxes2,2), source=nint((SMPD_SHRINK1/SMPD_SHRINK2) * real(self%positions1)))
        allocate(self%box_scores2(self%nboxes2),  source=-1.)
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box2, SMPD_SHRINK2)
        end do
        factor = real(OFFSET) * (SMPD_SHRINK1 / SMPD_SHRINK2)
        if( self%refpick )then
            !$omp parallel do schedule(static) default(shared) proc_bind(close)&
            !$omp private(ibox,rpos,xrange,yrange,box_score,xind,yind,ithr,outside,iref,scores,box_score_trial,l_err)
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
                        call boximgs_heap(ithr)%prenorm4real_corr(l_err)
                        if( l_err )then
                            box_score_trial = -1.
                        else
                            do iref = 1,self%nrefs
                                if( self%l_err_refs(iref,2) )then
                                    scores(iref) = -1.
                                else
                                    scores(iref) = self%boxrefs(iref,2)%real_corr_prenorm(boximgs_heap(ithr))
                                endif
                            end do
                            box_score_trial = maxval(scores)
                        endif
                        if( box_score_trial > box_score )then
                            self%positions2(ibox,:) = [xind,yind]
                            box_score = box_score_trial
                        endif
                    end do
                end do
                self%box_scores2(ibox) = box_score
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
                self%box_scores2(ibox) = box_score
            end do
            !$omp end parallel do
        endif
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%kill
        end do
    end subroutine refine_positions

    subroutine distance_filter( self, nbox, box_scores, positions, threshold )
        class(picker_utils), intent(inout) :: self
        integer,             intent(inout) :: nbox
        real,                intent(in)    :: threshold
        real, allocatable,               intent(inout) :: box_scores(:)
        integer, allocatable,            intent(inout) :: positions(:,:)
        integer, allocatable :: positions_tmp(:,:)
        real,    allocatable :: box_scores_tmp(:)
        integer :: ibox, jbox, loc, npeaks
        real    :: dist
        logical :: mask(nbox), selected_pos(nbox)
        selected_pos = .true.
        do ibox = 1,nbox
            mask = .false.
            !$omp parallel do schedule(static) default(shared) private(jbox,dist) proc_bind(close)
            do jbox = 1,nbox
                dist = euclid(real(positions(ibox,:)),real(positions(jbox,:)))
                if( dist <= threshold ) mask(jbox) = .true.
            end do
            !$omp end parallel do
            ! find best match in the neigh
            loc = maxloc(box_scores, mask=mask, dim=1)
            ! eliminate all but the best
            mask(loc) = .false.
            where( mask ) selected_pos = .false.
        end do
        npeaks = count(selected_pos)
        write(logfhandle,'(a,1x,I5)') '# positions before distance filtering: ', nbox
        write(logfhandle,'(a,1x,I5)') '# positions after  distance filtering: ', npeaks
        ! update positions2 and box_scores2
        allocate(positions_tmp(npeaks,2), box_scores_tmp(npeaks))
        positions_tmp  = 0
        box_scores_tmp = 0.
        npeaks         = 0
        do ibox = 1,nbox
            if( selected_pos(ibox) )then
                npeaks = npeaks + 1
                positions_tmp(npeaks,:) = positions(ibox,:)
                box_scores_tmp(npeaks)  = box_scores(ibox)
            endif
        end do
        nbox = npeaks
        deallocate(positions, box_scores)
        call move_alloc(positions_tmp, positions)
        call move_alloc(box_scores_tmp, box_scores)
    end subroutine distance_filter

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
        self%nboxes2 = npeaks
        call move_alloc(positions_tmp, self%positions2)
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
        logical     :: mask_backgr(0:self%nx2,0:self%ny2), outside, l_err
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
            !$omp parallel do schedule(static) default(shared) proc_bind(close) private(ibox,ithr,outside,iref,scores,l_err)
            do ibox = 1,nbackgr
                ithr = omp_get_thread_num() + 1
                call self%mic_shrink2%window_slim([positions_backgr(ibox,1),positions_backgr(ibox,2)], self%ldim_box2(1), boximgs_heap(ithr), outside)
                call boximgs_heap(ithr)%prenorm4real_corr(l_err)
                if( l_err )then
                    scores_nonpeak(ibox) = 0.0
                else
                    do iref = 1,self%nrefs
                        if( self%l_err_refs(iref,2))then
                            scores(iref) = 0.0
                        else
                            scores(iref) = self%boxrefs(iref,2)%real_corr_prenorm(boximgs_heap(ithr))
                        endif
                    end do
                    scores_nonpeak(ibox) = maxval(scores)
                endif
            end do
            !$omp end parallel do
        else
            !$omp parallel do schedule(static) default(shared) proc_bind(close) private(ibox,ithr,outside)
            do ibox = 1,nbackgr
                ithr = omp_get_thread_num() + 1
                call self%mic_shrink2%window_slim([positions_backgr(ibox,1),positions_backgr(ibox,2)], self%ldim_box2(1), boximgs_heap(ithr), outside)
                scores_nonpeak(ibox) = self%imgau_shrink2%real_corr_prenorm(boximgs_heap(ithr), self%sxx_shrink2)
            end do
            !$omp end parallel do
        endif
        if( self%refpick )then
            !$omp parallel do schedule(static) default(shared) proc_bind(close) private(ibox,ithr,outside,iref,scores,l_err)
            do ibox = 1,self%nboxes2
                ithr = omp_get_thread_num() + 1
                call self%mic_shrink2%window_slim([self%positions2(ibox,1),self%positions2(ibox,2)], self%ldim_box2(1), boximgs_heap(ithr), outside)
                call boximgs_heap(ithr)%prenorm4real_corr(l_err)
                if( l_err )then
                    scores_nonpeak(ibox) = 0.0
                else
                    do iref = 1,self%nrefs
                        if( self%l_err_refs(iref,2))then
                            scores(iref) = 0.0
                        else
                            scores(iref) = self%boxrefs(iref,2)%real_corr_prenorm(boximgs_heap(ithr))
                        endif
                    end do
                    scores_peak(ibox) = maxval(scores)
                endif
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
        ! cleanup
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%kill
        end do
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
            if( allocated(self%box_scores2) ) deallocate(self%box_scores2)
            if( allocated(self%fbody)       ) deallocate(self%fbody)
            if( allocated(self%l_mic_mask1) ) deallocate(self%l_mic_mask1)
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
