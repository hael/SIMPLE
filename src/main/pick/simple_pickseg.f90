module simple_pickseg
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,   only: params_glob
use simple_image,        only: image
use simple_tvfilter,     only: tvfilter
use simple_segmentation, only: otsu_img, sauvola
use simple_image_bin,    only: image_bin
use simple_syslib
implicit none

public :: pickseg
private
#include "simple_local_flags.inc"

! class constants
real,    parameter :: SHRINK   = 1.
real,    parameter :: LAMBDA   = 3.
logical, parameter :: L_WRITE  = .true.
logical, parameter :: L_DEBUG  = .true.

! class variables
integer      :: ldim_raw(3)
real         :: smpd_raw
type(image)  :: mic_raw
type(string) :: fbody

! instance
type pickseg
    private
    real               :: smpd_shrink = 0.
    integer            :: ldim(3), nboxes = 0, box_raw = 0
    real, allocatable  :: masscens(:,:)
    type(image_bin)    :: mic_shrink, img_cc
    type(stats_struct) :: sz_stats, diam_stats
    logical            :: exists = .false.
contains
    procedure :: pick
    procedure :: get_boxsize
    procedure :: get_ldim
    procedure :: get_nboxes
    procedure :: get_positions
    procedure :: get_smpd_shrink
    procedure :: report_boxfile
end type pickseg

contains

    subroutine pick( self, micname, is_AFM, moldiam, winsz )
        class(pickseg), intent(inout) :: self
        class(string),  intent(in)    :: micname !< micrograph file name
        real,             allocatable :: diams(:)
        integer,          allocatable :: sz(:)
        type(string)    :: ext
        type(tvfilter)  :: tvf
        type(image)     :: img_win
        type(image_bin) :: img_sdevs
        real    :: px(3), otsu_t
        integer :: i, boxcoord(2), sz_max, sz_min, nframes, box
        logical :: outside, is_AFM_l
        logical, optional, intent(in) :: is_AFM
        real,    optional, intent(in) :: moldiam
        integer, optional, intent(in) :: winsz
        logical :: l_moldiam, l_winsz
        is_AFM_l = .false.; l_moldiam= .false.; l_winsz= .false.
        if( present(moldiam) ) l_moldiam = .true.
        if( present(is_AFM)  ) is_AFM_l  = is_AFM 
        if( present(winsz)   ) l_winsz   = .true.
        ! set micrograph info
        call find_ldim_nptcls(micname, ldim_raw, nframes, smpd_raw)
        if( ldim_raw(3) /= 1 .or. nframes /= 1 ) THROW_HARD('Only for 2D images; pick')
        ! read micrograph
        call mic_raw%new(ldim_raw, smpd_raw)
        call mic_raw%read(micname)
        ! set fbody
        ext   = fname2ext(micname)
        fbody = get_fbody(basename(micname), ext)
        ! shrink micrograph
        self%ldim(1)     = round2even(real(ldim_raw(1))/SHRINK)
        self%ldim(2)     = round2even(real(ldim_raw(2))/SHRINK)
        self%ldim(3)     = 1
        self%smpd_shrink = smpd_raw * SHRINK
        call mic_raw%mul(real(product(ldim_raw))) ! to prevent numerical underflow when performing FFT
        call mic_raw%fft()
        call self%mic_shrink%new_bimg(self%ldim, self%smpd_shrink)
        call self%mic_shrink%set_ft(.true.)
        call mic_raw%clip(self%mic_shrink)
        select case( trim(params_glob%pcontrast) )
            case('black')
                ! flip contrast (assuming black particle contrast on input)
                call self%mic_shrink%mul(-1.)
            case('white')
                ! nothing to do
            case DEFAULT
                THROW_HARD('unknown pcontrast parameter, use (black|white)')
        end select
        ! low-pass filter micrograph
        call self%mic_shrink%bp(0., params_glob%lp)
        call self%mic_shrink%ifft
        call mic_raw%ifft
        if( L_WRITE ) call self%mic_shrink%write(fbody//'_lp.mrc')
        ! TV denoising
        call tvf%new()
        call tvf%apply_filter(self%mic_shrink, LAMBDA)
        call tvf%kill
        if( L_WRITE ) call self%mic_shrink%write(fbody//'_lp_tv.mrc')
        if( l_winsz )then
            call sauvola(self%mic_shrink, winsz, img_sdevs)
            call otsu_img(img_sdevs, otsu_t)
            call self%mic_shrink%transfer2bimg(img_sdevs)
        else
            call otsu_img(self%mic_shrink, otsu_t, tight=.true.)
        endif        
        call self%mic_shrink%set_imat
        if( L_WRITE ) call self%mic_shrink%write_bimg(fbody//'_lp_tv_bin.mrc')
        if( is_AFM_l )then
            call self%mic_shrink%erode 
            call self%mic_shrink%dilate
        else 
            call self%mic_shrink%erode
            call self%mic_shrink%erode
        endif 
        if( l_winsz )then 
            call self%mic_shrink%set_largestcc2background
            call self%mic_shrink%inv_bimg()
        endif
        if( L_WRITE ) call self%mic_shrink%write_bimg(fbody//'_lp_tv_bin_erode.mrc')
        ! identify connected components
        call self%mic_shrink%find_ccs(self%img_cc)
        if( L_WRITE ) call self%img_cc%write_bimg(fbody//'_lp_tv_bin_erode_cc.mrc')
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
        if( l_winsz )then
            ! elimination of connected components that are too small or too big with sauvola is not needed
        else
            sz_min = nint(self%sz_stats%avg - params_glob%ndev * self%sz_stats%sdev)
            sz_max = nint(self%sz_stats%avg + params_glob%ndev * self%sz_stats%sdev)
            ! optimal molecular diameter is inputted
            if( l_moldiam )then
                box      = nint(PI*(moldiam/2.)**2 / smpd_raw**2)
                sz_min   = box / 4
                sz_max   = box * 4
                !write(logfhandle,'("selecting cc [",i0," : ",i0,"]" )') sz_min, sz_max
            endif
            call self%img_cc%elim_ccs([sz_min,sz_max])
        endif
        call self%img_cc%get_nccs(self%nboxes)
        sz = self%img_cc%size_ccs()
        call calc_stats(real(sz), self%sz_stats)
        if( L_DEBUG )then
            print *, 'nboxes after elimination: ', self%nboxes
            print *, 'avg size: ', self%sz_stats%avg
            print *, 'med size: ', self%sz_stats%med
            print *, 'sde size: ', self%sz_stats%sdev
            print *, 'min size: ', self%sz_stats%minv
            print *, 'max size: ', self%sz_stats%maxv
        endif
        allocate(diams(self%nboxes), source=0.)
        do i = 1, self%nboxes
            call self%img_cc%diameter_cc(i, diams(i))
        enddo
        call calc_stats(diams, self%diam_stats)
        if( L_DEBUG )then 
            print *, 'CC diameter (in Angs) statistics'
            print *, 'avg diam: ', self%diam_stats%avg
            print *, 'med diam: ', self%diam_stats%med
            print *, 'sde diam: ', self%diam_stats%sdev
            print *, 'min diam: ', self%diam_stats%minv
            print *, 'max diam: ', self%diam_stats%maxv
        endif 
        self%box_raw = find_magic_box(2 * nint(self%diam_stats%med/smpd_raw))
        if( is_AFM_l  ) self%box_raw = find_magic_box(6 * nint(self%diam_stats%med/smpd_raw))
        if( l_moldiam ) self%box_raw = find_magic_box(3 * nint(moldiam/smpd_raw))
        call img_win%new([self%box_raw,self%box_raw,1], smpd_raw)
        if( allocated(self%masscens) ) deallocate(self%masscens)
        allocate(self%masscens(self%nboxes,2), source=0.)
        do i = 1, self%nboxes
            px = center_mass_cc(i)
            self%masscens(i,:2) = px(:2)
            if( L_DEBUG )then
                boxcoord = nint((real(SHRINK)*self%masscens(i,:2))-real(self%box_raw)/2.)
                call mic_raw%window_slim(boxcoord, self%box_raw, img_win, outside)
                call img_win%write(fbody//'_extracted.mrc', i)
            endif
        enddo
        call mic_raw%kill()
        call img_win%kill()
        call self%mic_shrink%kill_bimg()
        call self%img_cc%kill_bimg()
        if( allocated(sz) ) deallocate(sz)

        contains

            function center_mass_cc( i_cc ) result( px )
                integer, intent(in) :: i_cc
                real :: px(3)
                integer, allocatable :: pos(:,:)
                integer, allocatable :: imat_cc(:,:,:)
                imat_cc = int(self%img_cc%get_rmat())
                where( imat_cc .ne. i_cc ) imat_cc = 0
                call get_pixel_pos(imat_cc,pos)
                px(1) = sum(pos(1,:))/real(size(pos,dim = 2))
                px(2) = sum(pos(2,:))/real(size(pos,dim = 2))
                px(3) = 1.
                if( allocated(imat_cc) ) deallocate(imat_cc)
            end function center_mass_cc

            ! let's center the mass based on the area of particles within the box for crowded micrographs. weight center of mass with area of ccs within box.  
    end subroutine pick

    ! getters/setters

    pure function get_boxsize( self ) result( boxsize )
        class(pickseg), intent(in) :: self
        integer :: boxsize
        boxsize = self%box_raw
    end function get_boxsize

    pure function get_ldim( self ) result( ldim )
        class(pickseg), intent(in) :: self
        integer :: ldim(3)
        ldim = self%ldim
    end function get_ldim

    pure function get_nboxes( self ) result( nboxes )
        class(pickseg), intent(in) :: self
        integer :: nboxes
        nboxes = self%nboxes
    end function get_nboxes

    pure function get_smpd_shrink( self ) result( smpd_shrink )
        class(pickseg), intent(in) :: self
        real :: smpd_shrink
        smpd_shrink = self%smpd_shrink
    end function get_smpd_shrink

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
            enddo
        else
            do ibox = 1,self%nboxes
                pos(ibox,:) = nint((real(SHRINK)*self%masscens(ibox,:2))-real(self%box_raw)/2.)
            enddo
        endif
    end subroutine get_positions

    ! for writing boxes with arbitrary box size
    subroutine report_boxfile( self, fname, nptcls, box )
        class(pickseg),    intent(in) :: self
        class(string),     intent(in) :: fname
        integer,           intent(out):: nptcls
        integer, optional, intent(in) :: box
        integer, allocatable :: pos(:,:)
        integer :: funit, ibox, iostat
        nptcls = self%nboxes
        if( nptcls == 0 ) return
        call self%get_positions(pos, box)
        call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=iostat)
        call fileiochk('simple_pickseg; write_boxfile ', iostat)
        do ibox = 1,size(pos,dim=1)
            write(funit,'(I7,I7,I7,I7,I7)') pos(ibox,1), pos(ibox,2), self%box_raw, self%box_raw, -3
        enddo
        call fclose(funit)
    end subroutine report_boxfile

end module simple_pickseg
