module simple_picker_utils
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,          only: image
use simple_radial_medians, only: radial_medians
implicit none

public :: picker_utils
private
#include "simple_local_flags.inc"

real,    parameter :: SMPD_TARGET1 = 4.0, SMPD_TARGET2 = 2.0, LAMBDA1 = 5., LAMBDA2 = 3.
integer, parameter :: OFFSET  = 3
logical, parameter :: L_WRITE = .true.

type picker_utils
    private
    integer                  :: ldim_raw(3) = 0 , ldim_shrink1(3) = 0 , ldim_shrink2(3) = 0
    integer                  :: ldim_box(3) = 0 , ldim_box1(3)    = 0 , ldim_box2(3)    = 0
    real                     :: smpd_raw    = 0., smpd_shrink1    = 0., smpd_shrink2    = 0.
    integer                  :: nboxes1 = 0, nboxes2 = 0, nx1 = 0, ny1 = 0, nx2 = 0, ny2 = 0
    type(image)              :: mic_shrink1, mic_shrink2
    type(stats_struct)       :: stats_ptcl, stats_bg
    type(radial_medians)     :: radmeds_obj1, radmeds_obj2
    type(image), pointer     :: mic_raw => null()
    type(image), allocatable :: boximgs1(:), boximgs2(:)
    real,        allocatable :: avg_sdev1(:,:), avg_sdev2(:,:), radmeds1(:,:), radmeds2(:,:)
    integer,     allocatable :: positions1(:,:)
    logical,     allocatable :: is_ptcl1(:), is_ptcl2(:)
  contains
    procedure          :: set_mics
    procedure          :: calc_fgbgstats_mic_shrink1
    procedure, private :: set_positions_1
    procedure, private :: set_positions_2
    generic            :: set_positions => set_positions_1, set_positions_2
    procedure, private :: set_pos_priv
    procedure          :: extract_boximgs1
    procedure          :: calc_radmeds1
    procedure          :: classify_ptcls1
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
        scale1 = self%smpd_raw / SMPD_TARGET1
        scale2 = self%smpd_raw / SMPD_TARGET2
        self%ldim_shrink1(1) = round2even(real(self%ldim_raw(1)) * scale1)
        self%ldim_shrink1(2) = round2even(real(self%ldim_raw(2)) * scale1)
        self%ldim_shrink1(3) = 1
        self%ldim_shrink2(1) = round2even(real(self%ldim_raw(1)) * scale2)
        self%ldim_shrink2(2) = round2even(real(self%ldim_raw(2)) * scale2)
        self%ldim_shrink2(3) = 1
        call self%mic_shrink1%new(self%ldim_shrink1, SMPD_TARGET1)
        call self%mic_shrink2%new(self%ldim_shrink2, SMPD_TARGET2)
        call self%mic_shrink1%set_ft(.true.)
        call self%mic_shrink2%set_ft(.true.)
        call self%mic_raw%fft
        call self%mic_raw%clip(self%mic_shrink1)
        call self%mic_raw%clip(self%mic_shrink2)
        if( present(bp_lp) )then
            call self%mic_shrink1%bp(bp_lp(1),bp_lp(2))
            call self%mic_shrink2%bp(bp_lp(1),bp_lp(2))
        endif
        call self%mic_raw%ifft
        call self%mic_shrink1%ifft
        call self%mic_shrink2%ifft
        if( L_WRITE )then
            call self%mic_shrink1%write('mic_shrink1.mrc')
            call self%mic_shrink2%write('mic_shrink2.mrc')
        endif
    end subroutine set_mics

    ! type stats_struct
    !     real :: avg  = 0.
    !     real :: med  = 0.
    !     real :: sdev = 0.
    !     real :: maxv = 0
    !     real :: minv = 0.
    ! end type stats_struct

    subroutine calc_fgbgstats_mic_shrink1( self, bp_lp )
        use simple_segmentation
        use simple_tvfilter, only: tvfilter
        class(picker_utils), intent(inout) :: self
        real,                intent(in)    :: bp_lp(2) !< high- and low-pass limits in A
        real, allocatable  :: vec(:)
        type(tvfilter)     :: tvf
        type(image)        :: img_copy
        real :: t
        call self%mic_shrink1%fft
        call self%mic_shrink1%bp(bp_lp(1),bp_lp(2))
        call tvf%new()
        call tvf%apply_filter(self%mic_shrink1, LAMBDA1)
        call tvf%kill
        call self%mic_shrink1%ifft
        if( L_WRITE ) call self%mic_shrink1%write('tv_mic.mrc')
        vec = self%mic_shrink1%serialize()
        call otsu(size(vec), vec, t)
        if( L_WRITE )then
            call img_copy%copy(self%mic_shrink1)
            call img_copy%binarize(t)
            call img_copy%write('otsu_mic.mrc')
        endif
        call calc_stats(vec, self%stats_bg,   mask=vec >= t) ! background stats, assuming black particle density
        call calc_stats(vec, self%stats_ptcl, mask=vec <  t) ! particle   stats, assuming black particle density

        print *, 'background #pix/avg/sdev ', count(vec >= t), self%stats_bg%avg,   self%stats_bg%sdev
        print *, 'particle   #pix/avg/sdev ', count(vec <  t), self%stats_ptcl%avg, self%stats_ptcl%sdev
        print *, 'std_mean_diff ', std_mean_diff(self%stats_bg%avg, self%stats_ptcl%avg, self%stats_bg%sdev, self%stats_ptcl%sdev)

        ! results from filtered
        ! background #pix/avg/sdev       435326  0.964909315       7.91966170E-03
        ! particle   #pix/avg/sdev       261674  0.940358341       1.32814981E-02
        ! std_mean_diff    2.24531102

        ! results from raw
        ! background #pix/avg/sdev       435326  0.971406519       4.82148230E-02
        ! particle   #pix/avg/sdev       261674  0.929549515       5.25608510E-02
        ! std_mean_diff   0.829925179

        call img_copy%kill
    end subroutine calc_fgbgstats_mic_shrink1

    subroutine set_positions_1( self, maxdiam )
        class(picker_utils), intent(inout) :: self
        real,                intent(in)    :: maxdiam !< maximum diameter in A
        ! set logical dimensions of boxes
        self%ldim_box(1)  = round2even(maxdiam / self%smpd_raw)
        self%ldim_box(2)  = self%ldim_box(1)
        self%ldim_box(3)  = 1
        self%ldim_box1(1) = round2even(maxdiam / SMPD_TARGET1)
        self%ldim_box1(2) = self%ldim_box1(1)
        self%ldim_box1(3) = 1
        self%ldim_box2(1) = round2even(maxdiam / SMPD_TARGET2)
        self%ldim_box2(2) = self%ldim_box2(1)
        self%ldim_box2(3) = 1
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
        ! make radial medians objects
        call self%radmeds_obj1%new(self%ldim_box1)
        call self%radmeds_obj2%new(self%ldim_box2)
        ! set # pixels in x/y for both box sizes
        self%nx1 = self%ldim_shrink1(1) - self%ldim_box1(1)
        self%ny1 = self%ldim_shrink1(2) - self%ldim_box1(2)
        self%nx2 = self%ldim_shrink2(1) - self%ldim_box2(1)
        self%ny2 = self%ldim_shrink2(2) - self%ldim_box2(2)
        ! count # boxes
        self%nboxes1 = 0
        do xind = 0,self%nx1,OFFSET
            do yind = 0,self%ny1,OFFSET
                self%nboxes1 = self%nboxes1 + 1
            end do
        end do
        ! allocate and set positions1 
        if( allocated(self%positions1) ) deallocate(self%positions1)
        allocate(self%positions1(self%nboxes1,2), source=0)
        self%nboxes1 = 0
        do xind = 0,self%nx1,OFFSET
            do yind = 0,self%ny1,OFFSET
                self%nboxes1 = self%nboxes1 + 1
                self%positions1(self%nboxes1,:) = [xind,yind]
            end do
        end do
    end subroutine set_pos_priv

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
            call self%boximgs1(ibox)%new(self%ldim_box1, SMPD_TARGET1)
            call self%mic_shrink1%window_slim(self%positions1(ibox,:), self%ldim_box1(1), self%boximgs1(ibox), outside)
        end do
        !$omp end parallel do
    end subroutine extract_boximgs1

    subroutine calc_radmeds1( self )
        class(picker_utils), intent(inout) :: self
        type(stats_struct) :: stats
        integer :: ibox
        if( .not. allocated(self%positions1) ) THROW_HARD('positions need to be set before caluclating radial medians')
        if( .not. allocated(self%boximgs1)   ) THROW_HARD('boximgs1 need to be extracted before caluclating radial medians')
        if( allocated(self%avg_sdev1) ) deallocate(self%avg_sdev1)
        if( allocated(self%radmeds1)  ) deallocate(self%radmeds1)
        allocate(self%avg_sdev1(self%nboxes1,2), self%radmeds1(self%nboxes1,self%radmeds_obj1%get_rad_max()), source=0.)
        !$omp parallel do schedule(static) default(shared) private(ibox,stats) proc_bind(close)
        do ibox = 1,self%nboxes1
            call self%radmeds_obj1%calc_radial_medians(self%boximgs1(ibox), stats, self%radmeds1(ibox,:))
            self%avg_sdev1(ibox,1) = stats%avg
            self%avg_sdev1(ibox,2) = stats%sdev
        end do
        !$omp end parallel do
    end subroutine calc_radmeds1

    subroutine classify_ptcls1( self )
        class(picker_utils), intent(inout) :: self
        integer :: ibox, nptcls
        real    :: smd_ptcl, smd_bg
        if( .not. allocated(self%avg_sdev1) ) THROW_HARD('avg_sdev1 needs to be set (through calc_radmeds1) before call to classify_ptcls1')
        if( allocated(self%is_ptcl1) ) deallocate(self%is_ptcl1)
        allocate( self%is_ptcl1(self%nboxes1), source=.false.)
        do ibox = 1,self%nboxes1
            smd_ptcl = std_mean_diff(self%avg_sdev1(ibox,1), self%stats_ptcl%avg, self%avg_sdev1(ibox,2), self%stats_ptcl%sdev)
            smd_bg   = std_mean_diff(self%avg_sdev1(ibox,1), self%stats_bg%avg,   self%avg_sdev1(ibox,2), self%stats_bg%sdev)
            if( smd_ptcl < smd_bg ) self%is_ptcl1(ibox) = .true.
        end do
        nptcls = count(self%is_ptcl1)
        print *, 'found ', nptcls, ' particles among ', self%nboxes1, ' positions, % ', 100. * (real(nptcls) / real(self%nboxes1))
        if( L_WRITE )then
            call write_boximgs(self%nboxes1,       self%is_ptcl1, self%boximgs1, 'classified_as_ptcls.mrcs' )
            call write_boximgs(self%nboxes1, .not. self%is_ptcl1, self%boximgs1, 'classified_as_bg.mrcs' )
        endif
    end subroutine classify_ptcls1

    ! utilities

    subroutine write_boximgs( n, mask, boximgs, fname )
        integer,          intent(in)    :: n
        logical,          intent(in)    :: mask(n)
        class(image),     intent(inout) :: boximgs(n)
        character(len=*), intent(in)    :: fname
        integer :: ibox, cnt
        cnt = 0
        do ibox = 1,n
            if( mask(ibox) )then
                cnt = cnt + 1
                call boximgs(ibox)%write(fname, cnt)
            endif
        end do
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
