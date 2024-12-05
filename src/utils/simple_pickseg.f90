module simple_pickseg
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,   only: params_glob
use simple_image,        only: image
use simple_tvfilter,     only: tvfilter
use simple_segmentation, only: otsu_img
use simple_binimage,     only: binimage
use simple_syslib
implicit none

public :: pickseg
private
#include "simple_local_flags.inc"

! class constants
real,    parameter :: SHRINK   = 1.
real,    parameter :: LAMBDA   = 3.
logical, parameter :: L_WRITE  = .true.
logical, parameter :: L_DEBUG  = .false.

! class variables
integer                       :: ldim_raw(3)
real                          :: smpd_raw
type(image)                   :: mic_raw
character(len=:), allocatable :: fbody

! instance
type pickseg
    real               :: smpd_shrink = 0.
    integer            :: ldim(3), ldim_box(3), nboxes = 0, box_raw = 0
    real, allocatable  :: masscens(:,:)
    type(binimage)     :: mic_shrink, img_cc
    type(stats_struct) :: sz_stats, diam_stats
    logical            :: exists = .false.
contains
    procedure :: pick
    procedure :: get_positions
    procedure :: get_nboxes
    procedure :: report_boxfile
end type pickseg

contains

    subroutine pick( self, micname, is_AFM )
        class(pickseg), intent(inout) :: self
        character(len=*), intent(in)  :: micname !< micrograph file name
        real,             allocatable :: diams(:)
        integer,          allocatable :: sz(:)
        character(len=:), allocatable :: ext
        type(tvfilter) :: tvf
        type(image)    :: img_win
        real    :: px(3), otsu_t
        integer :: i, boxcoord(2), sz_max, sz_min, nframes
        logical :: outside, is_AFM_l
        logical, optional, intent(in) :: is_AFM 
        is_AFM_l = .false. 
        if( present(is_AFM) ) is_AFM_l = is_AFM 
        ! set micrograph info
        call find_ldim_nptcls(micname, ldim_raw, nframes, smpd_raw)
        if( ldim_raw(3) /= 1 .or. nframes /= 1 ) THROW_HARD('Only for 2D images')
        ! read micrograph
        call mic_raw%new(ldim_raw, smpd_raw)
        call mic_raw%read(micname)
        ! set fbody
        ext   = fname2ext(trim(micname))
        fbody = trim(get_fbody(basename(trim(micname)), ext))
        ! shrink micrograph
        self%ldim(1)     = round2even(real(ldim_raw(1))/SHRINK)
        self%ldim(2)     = round2even(real(ldim_raw(2))/SHRINK)
        self%ldim(3)     = 1
        self%smpd_shrink = smpd_raw * SHRINK
        call mic_raw%mul(real(product(ldim_raw))) ! to prevent numerical underflow when performing FFT
        call mic_raw%fft
        call self%mic_shrink%new_bimg(self%ldim, self%smpd_shrink)
        call self%mic_shrink%set_ft(.true.)
        call mic_raw%clip(self%mic_shrink)
        select case(trim(params_glob%pcontrast))
            case('black')
                ! flip contrast (assuming black particle contrast on input)
                call self%mic_shrink%mul(-1.)
            case('white')
                ! nothing to do
            case DEFAULT
                THROW_HARD('uknown pcontrast parameter, use (black|white)')
        end select
        ! low-pass filter micrograph
        call self%mic_shrink%bp(0., params_glob%lp)
        call self%mic_shrink%ifft
        call mic_raw%ifft
        if( L_WRITE ) call self%mic_shrink%write('mic_shrink_lp.mrc')
        ! TV denoising
        call tvf%new()
        call tvf%apply_filter(self%mic_shrink, LAMBDA)
        call tvf%kill
        if( L_WRITE ) call self%mic_shrink%write('mic_shrink_lp_tv.mrc')
        call otsu_img(self%mic_shrink, otsu_t)
        call self%mic_shrink%set_imat
        if( L_WRITE ) call self%mic_shrink%write_bimg('mic_shrink_lp_tv_bin.mrc')
        if(is_AFM_l) then
            call self%mic_shrink%erode 
            call self%mic_shrink%dilate
        else 
            call self%mic_shrink%erode
            call self%mic_shrink%erode
        end if 
        if( L_WRITE ) call self%mic_shrink%write_bimg('mic_shrink_lp_tv_bin_erode.mrc')
        ! identify connected components
        call self%mic_shrink%find_ccs(self%img_cc)
        if( L_WRITE ) call self%img_cc%write_bimg('mic_shrink_lp_tv_bin_erode_cc.mrc')
        call self%img_cc%get_nccs(self%nboxes)
        ! eliminate connected components that are too large or too small
        sz = self%img_cc%size_ccs()
       
        call calc_stats(real(sz), self%sz_stats)
        if( L_DEBUG )then
            print *, 'nboxes before elimination: ', self%nboxes
            print *, 'avg size: ', self%sz_stats%avg
            print *, 'med size: ', self%sz_stats%med
            print *, 'sde size: ', self%sz_stats%sdev
            print *, 'min size: ', self%sz_stats%minv
            print *, 'max size: ', self%sz_stats%maxv
        endif
        sz_min = nint(self%sz_stats%avg - params_glob%ndev * self%sz_stats%sdev)
        sz_max = nint(self%sz_stats%avg + params_glob%ndev * self%sz_stats%sdev)
        call self%img_cc%elim_ccs([sz_min,sz_max])
        call self%img_cc%get_nccs(self%nboxes)
        sz = self%img_cc%size_ccs()
        
        call calc_stats(real(sz), self%sz_stats)
        if( L_DEBUG )then
            print *, 'nboxes after  elimination: ', self%nboxes
            print *, 'avg size: ', self%sz_stats%avg
            print *, 'med size: ', self%sz_stats%med
            print *, 'sde size: ', self%sz_stats%sdev
            print *, 'min size: ', self%sz_stats%minv
            print *, 'max size: ', self%sz_stats%maxv
        endif
        allocate(diams(self%nboxes), source=0.)
        call calc_stats(diams, self%diam_stats)
        do i = 1, self%nboxes
            call self%img_cc%diameter_cc(i, diams(i))
        end do
        call calc_stats(diams, self%diam_stats)
        if( L_DEBUG)then 
            print *, 'avg diam: ', self%diam_stats%avg
            print *, 'med diam: ', self%diam_stats%med
            print *, 'sde diam: ', self%diam_stats%sdev
            print *, 'min diam: ', self%diam_stats%minv
            print *, 'max diam: ', self%diam_stats%maxv
        end if 
        self%box_raw = find_magic_box(2 * nint(self%diam_stats%med/smpd_raw))
        if(is_AFM_l) self%box_raw = nint(1.5*find_magic_box(2 * nint(self%diam_stats%med/smpd_raw)))
        call img_win%new([self%box_raw,self%box_raw,1], smpd_raw)
        if( allocated(self%masscens) ) deallocate(self%masscens)
        allocate(self%masscens(self%nboxes,2), source=0.)
        do i = 1, self%nboxes
            px = center_mass_cc(i)
            self%masscens(i,:2) = px(:2)
            if( L_DEBUG )then
                boxcoord = nint((real(SHRINK)*self%masscens(i,:2))-real(self%box_raw)/2.)
                call mic_raw%window_slim(boxcoord, self%box_raw, img_win, outside)
                call img_win%write('extracted.mrc', i)
            endif
        end do
        contains

            function center_mass_cc( i_cc ) result( px )
                integer, intent(in) :: i_cc
                real :: px(3)
                integer, allocatable :: pos(:,:)
                integer, allocatable :: imat_cc(:,:,:)
                imat_cc = int(self%img_cc%get_rmat())
                where(imat_cc .ne. i_cc) imat_cc = 0
                call get_pixel_pos(imat_cc,pos)
                px(1) = sum(pos(1,:))/real(size(pos,dim = 2))
                px(2) = sum(pos(2,:))/real(size(pos,dim = 2))
                px(3) = 1.
                if(allocated(imat_cc)) deallocate(imat_cc)
            end function center_mass_cc

            ! let's center the mass based on the area of particles within the box for crowded micrographs. weight center of mass with area of ccs within box.  
    end subroutine pick

    subroutine get_positions( self, pos, box )
        class(pickseg),       intent(in)    :: self
        integer, allocatable, intent(inout) :: pos(:,:)
        integer, optional,    intent(in)    :: box
        integer :: ibox
        if( allocated(pos) ) deallocate(pos)
        allocate( pos(self%nboxes,2), source=0 )
        if( present(box) )then
            do ibox = 1,self%nboxes
                pos(ibox,:) = nint((real(SHRINK)*self%masscens(ibox,:2))-real(box)/2.)
            end do
        else
            do ibox = 1,self%nboxes
                pos(ibox,:) = nint((real(SHRINK)*self%masscens(ibox,:2))-real(self%box_raw)/2.)
            end do
        endif
    end subroutine get_positions

    pure function get_nboxes( self ) result( nboxes )
        class(pickseg), intent(in) :: self
        integer :: nboxes
        nboxes = self%nboxes
    end function get_nboxes

    ! for writing boxes with arbitrary box size
    subroutine report_boxfile( self, fname, nptcls, box )
        class(pickseg),    intent(in) :: self
        character(len=*),  intent(in) :: fname
        integer,           intent(out):: nptcls
        integer, optional, intent(in) :: box
        integer, allocatable :: pos(:,:)
        integer :: funit, ibox, iostat
        nptcls  = self%nboxes
        if( nptcls == 0 ) return
        call self%get_positions(pos, box)
        call fopen(funit, status='REPLACE', action='WRITE', file=trim(adjustl(fname)), iostat=iostat)
        call fileiochk('simple_pickgau; write_boxfile ', iostat)
        do ibox = 1,size(pos,dim=1)
            write(funit,'(I7,I7,I7,I7,I7)') pos(ibox,1), pos(ibox,2), self%box_raw, self%box_raw, -3
        end do
        call fclose(funit)
    end subroutine report_boxfile

    ! pickseg_multi diams or area. 

end module simple_pickseg
