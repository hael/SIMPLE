!@descr: postprocessing of a state's map (sharpening, low-pass, masking): the postprocess program's body and the
!! refine3D master's per-iteration step on its own command line
module simple_postprocess_strategy
use simple_core_module_api
use simple_parameters,   only: parameters
use simple_cmdline,      only: cmdline
use simple_sp_project,   only: sp_project
use simple_image,        only: image
use simple_string_utils, only: log_item
implicit none

public :: run_postprocess, postprocess_volume_from_files
private
#include "simple_local_flags.inc"

contains

    !> the postprocess program on cline: the state's volume and FSC from the project (or the command line), then
    !! postprocess_volume_from_files
    subroutine run_postprocess( cline )
        class(cmdline),  intent(inout) :: cline
        type(string)     :: fname_vol, fname_fsc, fname_even_unfil, fname_odd_unfil, pair_stem
        type(parameters) :: params
        type(sp_project) :: spproj
        type(image)      :: vol_unfil, vol_unfil_odd
        real    :: smpd
        integer :: state, box, fsc_box
        ! set defaults
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call cline%set('oritype', 'out')
        ! parse commad-line
        call params%new(cline)
        ! read project segment
        call spproj%read_segment(params%oritype, params%projfile)
        ! state
        if( cline%defined('state') )then
            state = params%state
        else
            state = 1
        endif
        ! check volume, get correct smpd & box
        if( cline%defined('imgkind') .and. trim(params%imgkind) /= 'unfil' )then
            call spproj%get_vol(params%imgkind, state, fname_vol, smpd, box)
        else
            call spproj%get_vol('vol', state, fname_vol, smpd, box)
        endif
        if( .not.file_exists(fname_vol) )then
            THROW_HARD('volume: '//fname_vol%to_char()//' does not exist')
        endif
        ! using the input volume for postprocessing
        if( cline%defined('vol'//int2str(state)) ) fname_vol = params%vols(state)
        ! imgkind=unfil: postprocess the average of the unfiltered (prior-free base) half pair beside the project
        ! volume, written as <vol>_unfil.mrc, for comparison with the shipped regularized map. The cutoff still
        ! comes from the _unfil pair beside the project volume (pair_stem).
        pair_stem = fname_vol
        if( cline%defined('imgkind') )then
            if( trim(params%imgkind) == 'unfil' )then
                fname_even_unfil = add2fbody(fname_vol, params%ext, '_even_'//trim(params%imgkind))
                fname_odd_unfil  = add2fbody(fname_vol, params%ext, '_odd_'//trim(params%imgkind))
                if( .not. file_exists(fname_even_unfil) .or. .not. file_exists(fname_odd_unfil) )then
                    write(logfhandle,'(A)') log_item('POSTPROCESS', 'missing '//fname_even_unfil%to_char()//' or '//&
                        &fname_odd_unfil%to_char())
                    THROW_HARD('imgkind='//trim(params%imgkind)//' requires that half pair beside the project volume')
                endif
                call vol_unfil%new([box,box,box], smpd)
                call vol_unfil_odd%new([box,box,box], smpd)
                call vol_unfil%read(fname_even_unfil)
                call vol_unfil_odd%read(fname_odd_unfil)
                call vol_unfil%add(vol_unfil_odd)
                call vol_unfil%mul(0.5)
                fname_vol = add2fbody(fname_vol, params%ext, '_'//trim(params%imgkind))
                call vol_unfil%write(fname_vol, del_if_exists=.true.)
                call vol_unfil%kill
                call vol_unfil_odd%kill
                write(logfhandle,'(A)') log_item('POSTPROCESS', 'source is the '//trim(params%imgkind)//&
                    &' pair average '//fname_vol%to_char())
            endif
        endif
        if( cline%defined('fsc') )then
            if( .not.file_exists(params%fsc) ) THROW_HARD('FSC file: '//params%fsc%to_char()//' not found')
            fname_fsc = params%fsc
        else
            call spproj%get_fsc(state, fname_fsc, fsc_box)
        endif
        call postprocess_volume_from_files(fname_vol, fname_fsc, box, smpd, params, cline, state, pair_stem=pair_stem)
        ! destruct
        call spproj%kill
    end subroutine run_postprocess

    !> Isotropic postprocess: cutoff = FSC=0.143 (FSC file, else the _unfil pair); Guinier B (lp < 5 A) from
    !! the _unfil pair average, capped at the cutoff shell; sqrt(2FSC/(1+FSC)) unless provenance says ML-shrunk;
    !! Butterworth at the cutoff; no mask if the support is already in the map. pair_stem: _unfil pair stem.
    subroutine postprocess_volume_from_files( fname_vol, fname_fsc, box, smpd, params, cline, state, pair_stem )
        use simple_butterworth, only: butterworth_filter
        use simple_vol_pproc_policy,    only: state_mask_is_compatible
        use simple_halfmap_diagnostics, only: read_support_provenance
        use simple_image_msk,           only: image_msk
        class(string),   intent(in)    :: fname_vol, fname_fsc
        integer,         intent(in)    :: box
        real,            intent(in)    :: smpd
        type(parameters),intent(inout) :: params
        class(cmdline),  intent(inout) :: cline
        integer,         intent(in)    :: state
        class(string), optional, intent(in) :: pair_stem
        real, allocatable :: fsc(:), res(:), bwfilter(:), optlp(:)
        type(string)     :: fname_mirr, fname_pproc, fname_lp, fname_envmsk
        type(string)     :: fname_even_unfil, fname_odd_unfil
        type(image)      :: vol_bfac, vol_no_bfac, vol_envmsk, vol_unfil, vol_unfil_odd
        type(image_msk)  :: bfac_envmsk
        real    :: fsc0143, fsc05, lplim, cap_gain
        integer :: ldim(3), ldim_unfil(3), nptcls_unfil, lp_find
        logical :: has_fsc, do_envfsc, msk_exists, msk_compatible, l_unfil_pair
        logical :: l_support_at_source, l_prov_constrained, l_prov_found, l_prov_kind, l_ml_shrunk
        character(len=32) :: prov_solve_kind, prov_support_kind
        character(len=32) :: cutoff_text, bfac_text, gain_text
        if( .not.file_exists(fname_vol) )then
            THROW_HARD('volume: '//fname_vol%to_char()//' does not exist')
        endif
        ! generate file names
        fname_pproc = basename(add2fbody(fname_vol,   params%ext, PPROC_SUFFIX))
        fname_mirr  = basename(add2fbody(fname_pproc, params%ext, MIRR_SUFFIX))
        fname_lp    = basename(add2fbody(fname_vol,   params%ext, LP_SUFFIX))
        ! read volume(s)
        ldim = [box,box,box]
        call vol_bfac%new(ldim, smpd)
        call vol_bfac%read(fname_vol)
        ! protocol and its measurements: doc/policies/3D/refine3D_policy.md section 11
        has_fsc   = .false.
        do_envfsc = .false.
        res = vol_bfac%get_res()
        ! the unfiltered (unregularized, unpriored) pair beside the map: the
        ! B-factor's source, and the cutoff's when no FSC file is given
        if( present(pair_stem) )then
            fname_even_unfil = add2fbody(pair_stem, params%ext, '_even_unfil')
            fname_odd_unfil  = add2fbody(pair_stem, params%ext, '_odd_unfil')
        else
            fname_even_unfil = add2fbody(fname_vol, params%ext, '_even_unfil')
            fname_odd_unfil  = add2fbody(fname_vol, params%ext, '_odd_unfil')
        endif
        l_unfil_pair = file_exists(fname_even_unfil) .and. file_exists(fname_odd_unfil)
        if( l_unfil_pair )then
            ! the pair must come from the assembly that produced this map
            call find_ldim_nptcls(fname_even_unfil, ldim_unfil, nptcls_unfil)
            l_unfil_pair = all(ldim_unfil == ldim)
            if( l_unfil_pair )then
                call find_ldim_nptcls(fname_odd_unfil, ldim_unfil, nptcls_unfil)
                l_unfil_pair = all(ldim_unfil == ldim)
            endif
            if( .not. l_unfil_pair ) write(logfhandle,'(A)') &
                &log_item('POSTPROCESS', 'ignoring the unfiltered pair at a different box than the map')
        endif
        if( l_unfil_pair )then
            call vol_unfil%new(ldim, smpd)
            call vol_unfil_odd%new(ldim, smpd)
            call vol_unfil%read(fname_even_unfil)
            call vol_unfil_odd%read(fname_odd_unfil)
        endif
        params%fsc = fname_fsc
        if( trim(params%fsc%to_char()) /= '' .and. file_exists(params%fsc) )then
            fsc = file2rarr(params%fsc)
            call get_resolution(fsc, res, fsc05, fsc0143)
            has_fsc = .true.
        else if( l_unfil_pair )then
            call vol_unfil%fft()
            call vol_unfil_odd%fft()
            call vol_unfil%fsc(vol_unfil_odd, fsc)
            call vol_unfil%ifft()
            call vol_unfil_odd%ifft()
            call get_resolution(fsc, res, fsc05, fsc0143)
            has_fsc = .true.
            write(logfhandle,'(A)') log_item('POSTPROCESS', 'no FSC file, the FSC of the unfiltered pair')
        else
            THROW_WARN('FSC file: '//params%fsc%to_char()//' not found and no unfiltered pair beside the volume')
            if( .not. cline%defined('lp') )then
                THROW_HARD('no method for low-pass filtering defined; give fsc|lp on command line; postprocess_volume_from_files')
            endif
        endif
        if( has_fsc )then
            lplim     = fsc0143
            do_envfsc = params%l_envfsc
            ! no shell reaches FSC=0.143 (a state without signal): cut at the first shell
            if( lplim < TINY ) lplim = res(1)
        else
            lplim = params%lp
        endif
        ! B from the unfiltered pair: a regularized map's prior steepens the Guinier slope (map = fallback only)
        if( cline%defined('bfac') )then
            ! already in params%bfac
        else if( lplim < 5. )then
            if( l_unfil_pair )then
                call vol_unfil%add(vol_unfil_odd)
                call vol_unfil%mul(0.5)
                ! Guinier over the particle, not the support: solvent noise is
                ! flat in resolution, so the slope of a support-wide fit is set
                ! by where the signal drops under that floor, not by the
                ! protein's falloff. On PfCRT the support is about half solvent
                ! and the un-windowed fit gave -25 against -105 windowed, i.e.
                ! 1.5x instead of 5x sharpening at the cutoff. Masks only this
                ! temporary spectrum-estimation copy; no output is modified.
                call bfac_envmsk%automask3D(params, vol_unfil, .false., &
                    &lp_override=params%envmsklp, l_report=.false.)
                call bfac_envmsk%apply_3Dmask(vol_unfil)
                params%bfac = vol_unfil%guinier_bfac(HPLIM_GUINIER, lplim)
                write(logfhandle,'(A,1X,F8.2)') '>>> B-FACTOR (DENSITY-WINDOWED UNFILTERED PAIR) DETERMINED TO:', &
                    &params%bfac
            else
                params%bfac = vol_bfac%guinier_bfac(HPLIM_GUINIER, lplim)
                write(logfhandle,'(A,1X,F8.2)') '>>> B-FACTOR (THE MAP ITSELF; a regularized map over-estimates) DETERMINED TO:', &
                    &params%bfac
            endif
        else
            params%bfac = 0.
        endif
        if( l_unfil_pair )then
            call vol_unfil%kill
            call vol_unfil_odd%kill
        endif
        call bfac_envmsk%kill_bimg
        ! an ML-regularized map (provenance solve_kind) already carries ~FSC shrinkage per shell from its prior
        prov_solve_kind = ''
        call read_support_provenance(fname_vol, l_prov_constrained, l_prov_found, &
            &solve_kind=prov_solve_kind, l_kind_found=l_prov_kind, support_kind=prov_support_kind)
        l_ml_shrunk = l_prov_kind .and. (trim(prov_solve_kind) == 'regularized' .or. &
            &trim(prov_solve_kind) == 'gridding_regularized')
        call vol_bfac%fft()
        call vol_no_bfac%copy(vol_bfac)
        call vol_bfac%apply_bfac(params%bfac)
        ! sharpening held at the cutoff shell's gain beyond it: exp(-B s^2/4) would outgrow the Butterworth
        lp_find  = max(1, min(box/2, calc_fourier_index(lplim, box, smpd)))
        cap_gain = 1.
        if( params%bfac < 0. )then
            call vol_bfac%apply_filter(bfac_cap_filter(vol_bfac%get_filtsz(), box, smpd, params%bfac, lp_find))
            cap_gain = exp(-(params%bfac / 4.) * (real(lp_find) / (real(box) * smpd))**2)
        endif
        ! FSC weighting (RELION's sqrt(2FSC/(1+FSC))) exactly once: not on a map its prior already shrank
        if( has_fsc .and. .not. l_ml_shrunk )then
            optlp = fsc2cref(fsc)
            where( res < TINY ) optlp = 0.
            call vol_bfac%apply_filter(optlp)
            call vol_no_bfac%apply_filter(optlp)
        endif
        ! Nyquist-sized: apply_filter zeroes shells beyond size(filter), the cube corners where the B gain peaks
        allocate(bwfilter(vol_bfac%get_filtsz()), source=0.)
        call butterworth_filter(lp_find, bwfilter)
        call vol_bfac%apply_filter(bwfilter)
        call vol_no_bfac%apply_filter(bwfilter)
        deallocate(bwfilter)
        write(bfac_text,'(F32.2)') params%bfac
        write(cutoff_text,'(F32.2)') lplim
        if( params%bfac < 0. )then
            write(gain_text,'(F32.2)') cap_gain
            write(logfhandle,'(A)') log_item('POSTPROCESS', 'sharpening gain capped at '//trim(adjustl(gain_text))//&
                &' from shell '//int2str(lp_find)//' to Nyquist')
        endif
        ! one line: the B-factor, the weighting (an ML-regularized map carries it in its prior) and the low-pass
        if( has_fsc .and. l_ml_shrunk )then
            write(logfhandle,'(A)') log_item('POSTPROCESS', 'B-factor '//trim(adjustl(bfac_text))//&
                &' A^2, ML-prior weighting, low-pass '//trim(adjustl(cutoff_text))//' A')
        else if( has_fsc )then
            write(logfhandle,'(A)') log_item('POSTPROCESS', 'B-factor '//trim(adjustl(bfac_text))//&
                &' A^2, FSC weighting, low-pass '//trim(adjustl(cutoff_text))//' A')
        else
            write(logfhandle,'(A)') log_item('POSTPROCESS', 'B-factor '//trim(adjustl(bfac_text))//&
                &' A^2, low-pass '//trim(adjustl(cutoff_text))//' A')
        endif
        ! write low-pass filtered without B-factor or mask & read the original back in
        call vol_no_bfac%ifft
        call vol_no_bfac%write(fname_lp)
        ! An explicit unmasked record overrides the backend's usual support convention.
        l_support_at_source = trim(params%rec_backend) == 'pcg'
        if( l_prov_found ) l_support_at_source = trim(prov_support_kind) /= 'none'
        call vol_bfac%ifft()
        if( l_support_at_source )then
            ! no post-hoc mask: the support is already in the map
        else if( do_envfsc )then
            ! enveloppe mask from FSC
            fname_envmsk = AUTOMASK_FBODY//trim(adjustl(int2str_pad(state,2)))//MRC_EXT
            call state_mask_is_compatible(fname_envmsk, box, smpd, msk_exists, msk_compatible)
            if( msk_exists .and. msk_compatible )then
                call vol_envmsk%new(ldim, smpd)
                call vol_envmsk%read(fname_envmsk)
                call vol_bfac%zero_env_background(vol_envmsk)
                call vol_bfac%mul(vol_envmsk)
            else
                THROW_WARN('enveloppe mask: '//fname_envmsk%to_char()//' has different dimensions&
                & or sampling than volume; using spherical mask instead')
                call vol_bfac%mask3D_soft(params%msk_crop)
            endif
        else
            ! spherical mask
            call vol_bfac%mask3D_soft(params%msk_crop)
        endif
        ! output in cwd
        call vol_bfac%write(fname_pproc)
        ! also output mirrored by default (unless otherwise stated on command line)
        if( .not. cline%defined('mirr') .or. params%mirr .ne. 'no' )then
            call vol_bfac%mirror('x')
            call vol_bfac%write(fname_mirr)
        endif
        ! destruct
        call vol_bfac%kill
        call vol_no_bfac%kill
        call vol_envmsk%kill
        call fname_mirr%kill
        call fname_pproc%kill
        call fname_lp%kill
        call fname_envmsk%kill
        call fname_even_unfil%kill
        call fname_odd_unfil%kill
    end subroutine postprocess_volume_from_files

end module simple_postprocess_strategy
