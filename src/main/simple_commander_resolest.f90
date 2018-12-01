! concrete commander: resolest for resolution estimation
module simple_commander_resolest
include 'simple_lib.f08'
use simple_parameters,     only: parameters, params_glob
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_image,          only: image
use simple_masker,         only: masker
implicit none

public :: fsc_commander
public :: local_res_commander
public :: local_res2D_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: fsc_commander
  contains
    procedure :: execute      => exec_fsc
end type fsc_commander
type, extends(commander_base) :: local_res_commander
  contains
    procedure :: execute      => exec_local_res
end type local_res_commander
type, extends(commander_base) :: local_res2D_commander
  contains
    procedure :: execute      => exec_local_res2D
end type local_res2D_commander

contains

    !> calculates Fourier shell correlation from Even/Odd Volume pairs
    subroutine exec_fsc( self, cline )
        class(fsc_commander), intent(inout) :: self
        class(cmdline),       intent(inout) :: cline
        type(parameters)  :: params
        type(image)       :: even, odd
        type(masker)      :: mskvol
        integer           :: j, find_plate
        real              :: res_fsc05, res_fsc0143
        real, allocatable :: res(:), corrs(:)
        call params%new(cline)
        ! read even/odd pair
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd%new([params%box,params%box,params%box], params%smpd)
        call odd%read(params%vols(1))
        call even%read(params%vols(2))
        if( cline%defined('mskfile') )then
            if( file_exists(params%mskfile) )then
                call mskvol%new([params%box,params%box,params%box], params%smpd)
                call mskvol%read(params%mskfile)
                call mskvol%resmask()
                call mskvol%write('resmask'//params%ext)
                call even%zero_background
                call odd%zero_background
                call even%mul(mskvol)
                call odd%mul(mskvol)
            else
                THROW_HARD('mskfile: '//trim(params%mskfile)//' does not exist in cwd; exec_fsc')
            endif
        else
            ! spherical masking
            if( params%l_innermsk )then
                call even%mask(params%msk, 'soft', inner=params%inner, width=params%width)
                call odd%mask(params%msk, 'soft', inner=params%inner, width=params%width)
            else
                call even%mask(params%msk, 'soft')
                call odd%mask(params%msk, 'soft')
            endif
        endif
        ! forward FT
        call even%fft()
        call odd%fft()
        ! calculate FSC
        res = even%get_res()
        allocate(corrs(even%get_filtsz()))
        call even%fsc(odd, corrs)
        if( params%l_phaseplate ) call phaseplate_correct_fsc(corrs, find_plate)
        do j=1,size(res)
           write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
        end do
        call get_resolution(corrs, res, res_fsc05, res_fsc0143)
        write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
        write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
        write(logfhandle,'(A,1X,F6.2)') '>>> MEDIAN FSC (SPECSCORE):', median_nocopy(corrs)
        call even%kill
        call odd%kill
        ! end gracefully
        call simple_end('**** SIMPLE_FSC NORMAL STOP ****')
    end subroutine exec_fsc

    !> calculates local resolution from Even/Odd Volume pairs
    subroutine exec_local_res( self, cline )
        use simple_estimate_ssnr, only: local_res, local_res_lp
        class(local_res_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters)     :: params
        type(image)          :: even, odd, vol2filter
        type(masker)         :: mskvol
        real                 :: res_fsc05, res_fsc0143
        logical              :: have_mask_file
        integer, allocatable :: locres_finds(:,:,:)
        real,    allocatable :: res(:), corrs(:)
        if( .not. cline%defined('lplim_crit') ) call cline%set('lplim_crit', 0.5)
        call params%new(cline)
        ! read even/odd pair
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd%new([params%box,params%box,params%box], params%smpd)
        call odd%read(params%vols(1))
        call even%read(params%vols(2))
        have_mask_file = .false.
        if( cline%defined('mskfile') )then
            if( file_exists(params%mskfile) )then
                call mskvol%new([params%box,params%box,params%box], params%smpd)
                call mskvol%read(params%mskfile)
                call mskvol%resmask()
                call mskvol%write('resmask'//params%ext)
                call even%zero_background
                call odd%zero_background
                call even%mul(mskvol)
                call odd%mul(mskvol)
                have_mask_file = .true.
            else
                THROW_HARD('mskfile: '//trim(params%mskfile)//' does not exist in cwd; exec_local_res')
            endif
        else
            ! spherical masking
            if( params%l_innermsk )then
                call even%mask(params%msk, 'soft', inner=params%inner, width=params%width)
                call odd%mask(params%msk, 'soft', inner=params%inner, width=params%width)
            else
                call even%mask(params%msk, 'soft')
                call odd%mask(params%msk, 'soft')
            endif
        endif
        ! forward FT
        call even%fft()
        call odd%fft()
        if( cline%defined('mskfile') )then
            call mskvol%read(params%mskfile)
            call mskvol%remove_edge
        else
            call mskvol%disc([params%box,params%box,params%box], params%smpd, params%msk)
        endif
        call local_res(even, odd, mskvol, params%lplim_crit, locres_finds)
        ! destruct
        call even%kill
        call odd%kill
        ! filter inputted vol3
        if( cline%defined('vol3') )then
            call vol2filter%new([params%box,params%box,params%box], params%smpd)
            call vol2filter%read(params%vols(1))
            call local_res_lp(locres_finds, vol2filter)
            call vol2filter%write('map_locres_lp.mrc')
            call vol2filter%kill
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_LOCAL_RES NORMAL STOP ****')
    end subroutine exec_local_res

    !> calculates local resolution from Even/Odd class averages
    subroutine exec_local_res2D( self, cline )
        use simple_estimate_ssnr, only: local_res2D, local_res2D_lp
        class(local_res2D_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters)         :: params
        type(image), allocatable :: even_avgs(:), odd_avgs(:), avgs2filter(:)
        integer              :: iptcl
        real                 :: res_fsc05, res_fsc0143
        logical              :: have_mask_file
        integer, allocatable :: locres_finds(:,:,:)
        real,    allocatable :: res(:), corrs(:)
        if( .not. cline%defined('lplim_crit') ) call cline%set('lplim_crit', 0.5)
        call params%new(cline)
        ! read even/odd pairs, mask, FFT
        allocate(even_avgs(params%nptcls), odd_avgs(params%nptcls))
        do iptcl=1,params%nptcls
            ! read
            call even_avgs(iptcl)%new([params%box,params%box,1], params%smpd)
            call odd_avgs(iptcl)%new([params%box,params%box,1], params%smpd)
            call even_avgs(iptcl)%read(params%stk, iptcl)
            call odd_avgs(iptcl)%read(params%stk2, iptcl)
            ! forward FT
            call even_avgs(iptcl)%norm
            call odd_avgs(iptcl)%norm
            call even_avgs(iptcl)%fft()
            call odd_avgs(iptcl)%fft()
        end do
        call local_res2D(even_avgs, odd_avgs, params%lplim_crit, locres_finds)
        ! destruct
        do iptcl=1,params%nptcls
            call even_avgs(iptcl)%kill
            call odd_avgs(iptcl)%kill
        enddo
        ! filter inputted stk3
        if( cline%defined('stk3') )then
            ! read the averages to be filtered
            allocate(avgs2filter(params%nptcls))
            do iptcl=1,params%nptcls
                call avgs2filter(iptcl)%new([params%box,params%box,1], params%smpd, wthreads=.false.)
                call avgs2filter(iptcl)%read(params%stk3, iptcl)
            end do
            call local_res2D_lp(locres_finds, avgs2filter)
            ! write and destruct
            do iptcl=1,params%nptcls
                call avgs2filter(iptcl)%write('cavgs_locres2D_lp.mrc', iptcl)
                call avgs2filter(iptcl)%kill
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_LOCAL_RES2D NORMAL STOP ****')
    end subroutine exec_local_res2D

end module simple_commander_resolest
