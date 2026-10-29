!@descr: standalone commander for model-free local Wiener sharpening of independent half-maps
module simple_commanders_halfmap_sharpen
use simple_defs,             only: logfhandle
use simple_commander_base,   only: commander_base
use simple_cmdline,          only: cmdline
use simple_parameters,       only: parameters
use simple_image,            only: image
use simple_halfmap_enhancer, only: halfmap_enhancer
use simple_imghead,          only: find_ldim_nptcls, find_img_smpd
use simple_fileio,           only: add2fbody, fname2ext, wait_for_closure
use simple_jiffys,           only: simple_end
use simple_string,           only: string
use simple_error,            only: simple_exception
implicit none
private
#include "simple_local_flags.inc"

public :: commander_sharpen

type, extends(commander_base) :: commander_sharpen
  contains
    procedure :: execute => exec_sharpen
end type commander_sharpen

contains

    subroutine exec_sharpen( self, cline )
        class(commander_sharpen), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(parameters) :: params
        type(image) :: even, odd, enhanced, confidence, local_resolution, envelope
        type(halfmap_enhancer) :: enhancer
        type(string) :: ext, confidence_file, locres_file, envelope_file
        integer :: ldim_even(3), ldim_odd(3), nsections
        real :: smpd, smpd_even, smpd_odd, noise_mean, noise_sdev
        logical :: smpd_override
        smpd_override = cline%defined('smpd')
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call find_ldim_nptcls(params%vol_even, ldim_even, nsections)
        call find_ldim_nptcls(params%vol_odd,  ldim_odd,  nsections)
        if( any(ldim_even /= ldim_odd) ) THROW_HARD('half-map dimensions differ; sharpen')
        if( ldim_even(3) < 2 ) THROW_HARD('sharpen requires 3D volumes')
        smpd_even = find_img_smpd(params%vol_even)
        smpd_odd  = find_img_smpd(params%vol_odd)
        if( smpd_override )then
            smpd = params%smpd
        else
            if( smpd_even <= 0. .or. smpd_odd <= 0. ) &
                &THROW_HARD('half-map headers lack sampling; provide smpd; sharpen')
            if( abs(smpd_even - smpd_odd) > 1.e-3 ) &
                &THROW_HARD('half-map sampling distances differ; sharpen')
            smpd = smpd_even
        endif
        call even%new(ldim_even, smpd)
        call odd%new(ldim_even, smpd)
        call even%read(params%vol_even)
        call odd%read(params%vol_odd)
        call enhancer%new(ldim_even, smpd)
        call enhancer%enhance(even, odd, enhanced, confidence, local_resolution, envelope)
        call enhancer%get_noise_stats(noise_mean, noise_sdev)
        ext = fname2ext(params%outvol)
        if( ext .eq. '' ) THROW_HARD('outvol needs an image-file extension; sharpen')
        ext = '.'//ext%to_char()
        confidence_file = add2fbody(params%outvol, ext, '_confidence')
        locres_file     = add2fbody(params%outvol, ext, '_locres')
        envelope_file   = add2fbody(params%outvol, ext, '_envelope')
        call enhanced%write(params%outvol, del_if_exists=.true.)
        call confidence%write(confidence_file, del_if_exists=.true.)
        call local_resolution%write(locres_file, del_if_exists=.true.)
        call envelope%write(envelope_file, del_if_exists=.true.)
        call wait_for_closure(params%outvol)
        call wait_for_closure(confidence_file)
        call wait_for_closure(locres_file)
        call wait_for_closure(envelope_file)
        write(logfhandle,'(A,I0,A,I0,A,F7.3)') '>>> HALF-MAP SHARPENING: box ', ldim_even(1), &
            &', local window ', enhancer%get_window(), ', smpd ', smpd
        write(logfhandle,'(A,ES12.4,A,ES12.4)') '>>> OUTER-SHELL NOISE: mean ', noise_mean, ', merged-map sdev ', noise_sdev
        write(logfhandle,'(A)') '>>> ENHANCED MAP: '//params%outvol%to_char()
        write(logfhandle,'(A)') '>>> CONFIDENCE MAP: '//confidence_file%to_char()
        write(logfhandle,'(A)') '>>> LOCAL RESOLUTION MAP: '//locres_file%to_char()
        write(logfhandle,'(A)') '>>> ENVELOPE MAP: '//envelope_file%to_char()
        write(logfhandle,'(A)') '>>> Interpretation products only; retain the input half-maps for FSC and resolution.'
        call enhancer%kill()
        call even%kill()
        call odd%kill()
        call enhanced%kill()
        call confidence%kill()
        call local_resolution%kill()
        call envelope%kill()
        call ext%kill()
        call confidence_file%kill()
        call locres_file%kill()
        call envelope_file%kill()
        call simple_end('**** SIMPLE_SHARPEN NORMAL STOP ****')
    end subroutine exec_sharpen

end module simple_commanders_halfmap_sharpen
