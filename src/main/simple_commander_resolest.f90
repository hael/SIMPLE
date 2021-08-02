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
public :: nonuniform_lp_commander
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
type, extends(commander_base) :: nonuniform_lp_commander
  contains
    procedure :: execute      => exec_nonuniform_lp
end type nonuniform_lp_commander

contains

    !> calculates Fourier shell correlation from Even/Odd Volume pairs
    subroutine exec_fsc( self, cline )
        class(fsc_commander), intent(inout) :: self
        class(cmdline),       intent(inout) :: cline
        type(parameters)  :: params
        type(image)       :: even, odd
        type(masker)      :: mskvol
        integer           :: j, find_plate, k_hp, k_lp
        real              :: res_fsc05, res_fsc0143
        real, allocatable :: res(:), corrs(:)
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
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
        k_hp = calc_fourier_index(params%hp, params%box, params%smpd)
        k_lp = calc_fourier_index(params%lp, params%box, params%smpd)
        write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
        write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
        write(logfhandle,'(A,1X,F8.4)') '>>> MEDIAN FSC (SPECSCORE):', median_nocopy(corrs(k_hp:k_lp))
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
        if( .not. cline%defined('mkdir') )      call cline%set('mkdir', 'yes')
        if( .not. cline%defined('lplim_crit') ) call cline%set('lplim_crit', 0.143)
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

    subroutine exec_nonuniform_lp( self, cline )
        use simple_estimate_ssnr, only: nonuniform_lp
        class(nonuniform_lp_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(parameters) :: params
        type(image)      :: even, odd
        type(masker)     :: mskvol
        logical          :: have_mask_file
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        ! read even/odd pair
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd%new([params%box,params%box,params%box],  params%smpd)
        call odd%read(params%vols(1))
        call even%read(params%vols(2))
        have_mask_file = .false.
        if( cline%defined('mskfile') )then
            if( file_exists(params%mskfile) )then
                call mskvol%new([params%box,params%box,params%box], params%smpd)
                call mskvol%read(params%mskfile)
                have_mask_file = .true.
            else
                THROW_HARD('mskfile: '//trim(params%mskfile)//' does not exist in cwd; exec_nonuniform_lp')
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
        if( cline%defined('mskfile') )then
            call mskvol%read(params%mskfile)
            call mskvol%remove_edge
        else
            call mskvol%disc([params%box,params%box,params%box], params%smpd, params%msk)
        endif
        call nonuniform_lp(even, odd, mskvol)
        ! destruct
        call even%kill
        call odd%kill
        ! end gracefully
        call simple_end('**** SIMPLE_nonuniform_lp NORMAL STOP ****')
    end subroutine exec_nonuniform_lp

end module simple_commander_resolest
