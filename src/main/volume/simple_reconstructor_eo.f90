!@descr: explicit gridding half-pair I/O, FSC, and restoration operations
module simple_reconstructor_eo
use simple_core_module_api
use simple_reconstructor,   only: reconstructor, gridding_half_restore
use simple_parameters,      only: parameters
use simple_image,           only: image
use simple_refine3D_fnames, only: refine3D_fsc_fname
use simple_fsc
implicit none

public :: gridding_pair_diagnostics, calculate_gridding_pair_fsc, &
    &read_gridding_pair_accumulators, restore_gridding_pair, &
    &write_gridding_pair_accumulators, write_gridding_pair_diagnostics
private
#include "simple_local_flags.inc"

type :: gridding_pair_diagnostics
    real, allocatable :: fsc(:)
    real              :: res_fsc05   = 0.
    real              :: res_fsc0143 = 0.
    real              :: cfar        = 0.
  contains
    procedure :: kill => kill_gridding_pair_diagnostics
end type gridding_pair_diagnostics

contains

    !> Read an explicit even/odd gridding-accumulator artifact. Previous
    !! smaller grids retain the legacy update-fraction zero-padding behavior.
    !! A missing artifact resets the pair (legacy partial semantics) unless the
    !! caller marks the artifact required, in which case it is a hard error.
    subroutine read_gridding_pair_accumulators( params, even_rec, odd_rec, fbody, required )
        use simple_imgfile, only: imgfile
        class(parameters),    intent(in)    :: params
        class(reconstructor), intent(inout) :: even_rec, odd_rec
        class(string),        intent(in)    :: fbody
        logical, optional,    intent(in)    :: required
        type(string)      :: even_vol, even_rho, odd_vol, odd_rho
        type(image)       :: prev_vol_e, prev_vol_o
        type(imgfile)     :: ioimg_e, ioimg_o
        real, allocatable :: rho_e(:,:,:), rho_o(:,:,:)
        integer :: current_ldim(3), cshape(3), prev_ldim(3)
        integer :: fhandle_rho_e, fhandle_rho_o, i, ierr, dummy
        real    :: current_smpd, prev_smpd
        logical :: here(4), l_pad_with_zeros
        current_ldim = even_rec%get_ldim()
        current_smpd = even_rec%get_smpd()
        if( any(odd_rec%get_ldim() /= current_ldim) )then
            THROW_HARD('gridding pair accumulator dimensions do not match')
        endif
        if( abs(odd_rec%get_smpd() - current_smpd) > TINY )then
            THROW_HARD('gridding pair accumulator sampling does not match')
        endif
        even_vol = fbody//'_even'//MRC_EXT
        even_rho = string('rho_')//fbody//'_even'//MRC_EXT
        odd_vol  = fbody//'_odd'//MRC_EXT
        odd_rho  = string('rho_')//fbody//'_odd'//MRC_EXT
        here(1)  = file_exists(even_vol)
        here(2)  = file_exists(even_rho)
        here(3)  = file_exists(odd_vol)
        here(4)  = file_exists(odd_rho)
        if( all(here) )then
            l_pad_with_zeros = .false.
            if( params%l_update_frac )then
                call find_ldim_nptcls(even_vol, prev_ldim, dummy)
                prev_smpd = current_smpd
                if( prev_ldim(1) == current_ldim(1) )then
                    ! matching grids require no compatibility transform
                elseif( prev_ldim(1) > current_ldim(1) )then
                    THROW_HARD('previous gridding accumulator is larger than the current grid')
                else
                    l_pad_with_zeros = .true.
                endif
            endif
            call fopen(fhandle_rho_e, file=even_rho, status='OLD', action='READ', access='STREAM', iostat=ierr)
            call fileiochk('read_gridding_pair_accumulators opening '//even_rho%to_char(), ierr)
            call fopen(fhandle_rho_o, file=odd_rho, status='OLD', action='READ', access='STREAM', iostat=ierr)
            call fileiochk('read_gridding_pair_accumulators opening '//odd_rho%to_char(), ierr)
            if( l_pad_with_zeros )then
                call ioimg_e%open(even_vol, prev_ldim, prev_smpd, formatchar='M', readhead=.false., rwaction='read')
                call ioimg_o%open(odd_vol, prev_ldim, prev_smpd, formatchar='M', readhead=.false., rwaction='read')
                call prev_vol_e%new(prev_ldim, prev_smpd)
                call prev_vol_o%new(prev_ldim, prev_smpd)
                cshape = [fdim(prev_ldim(1)), prev_ldim(2), prev_ldim(3)]
                allocate(rho_e(1:cshape(1),1:cshape(2),1:cshape(3)), &
                    &rho_o(1:cshape(1),1:cshape(2),1:cshape(3)))
                call even_rec%reset
                call odd_rec%reset
                !$omp parallel do default(shared) private(i,ierr) schedule(static) num_threads(4)
                do i = 1, 4
                    select case(i)
                        case(1)
                            call prev_vol_e%read_raw_mrc(ioimg_e)
                        case(2)
                            call prev_vol_o%read_raw_mrc(ioimg_o)
                        case(3)
                            read(fhandle_rho_e, pos=1, iostat=ierr) rho_e
                            if( ierr /= 0 ) call fileiochk('read_gridding_pair_accumulators reading even rho', ierr)
                        case(4)
                            read(fhandle_rho_o, pos=1, iostat=ierr) rho_o
                            if( ierr /= 0 ) call fileiochk('read_gridding_pair_accumulators reading odd rho', ierr)
                    end select
                enddo
                !$omp end parallel do
                call even_rec%pad_with_zeros(prev_vol_e, rho_e)
                call odd_rec%pad_with_zeros(prev_vol_o, rho_o)
                call prev_vol_e%kill
                call prev_vol_o%kill
                deallocate(rho_e, rho_o)
            else
                call ioimg_e%open(even_vol, current_ldim, current_smpd, formatchar='M', &
                    &readhead=.false., rwaction='read')
                call ioimg_o%open(odd_vol, current_ldim, current_smpd, formatchar='M', &
                    &readhead=.false., rwaction='read')
                !$omp parallel do default(shared) private(i) schedule(static) num_threads(4)
                do i = 1, 4
                    select case(i)
                        case(1)
                            call even_rec%read_raw_mrc(ioimg_e)
                        case(2)
                            call odd_rec%read_raw_mrc(ioimg_o)
                        case(3)
                            call even_rec%read_raw_rho(fhandle_rho_e)
                        case(4)
                            call odd_rec%read_raw_rho(fhandle_rho_o)
                    end select
                enddo
                !$omp end parallel do
            endif
            call ioimg_e%close
            call ioimg_o%close
            call fclose(fhandle_rho_e)
            call fclose(fhandle_rho_o)
        else
            if( present(required) )then
                if( required )then
                    THROW_HARD('required gridding pair accumulator artifact is incomplete: '//fbody%to_char())
                endif
            endif
            call even_rec%reset
            call odd_rec%reset
        endif
    end subroutine read_gridding_pair_accumulators

    !> Write an explicit even/odd gridding-accumulator artifact using the
    !! established numerator/rho filenames.
    subroutine write_gridding_pair_accumulators( even_rec, odd_rec, fbody )
        class(reconstructor), intent(inout) :: even_rec, odd_rec
        class(string),        intent(in)    :: fbody
        call even_rec%write_raw_accum(fbody//'_even'//MRC_EXT, &
            &string('rho_')//fbody//'_even'//MRC_EXT)
        call odd_rec%write_raw_accum(fbody//'_odd'//MRC_EXT, &
            &string('rho_')//fbody//'_odd'//MRC_EXT)
    end subroutine write_gridding_pair_accumulators

    !> Restore one explicit gridding half pair. This owns backend-specific
    !! base/replay mechanics but no composite lifetime, partial reduction, or
    !! trailing-chain policy.
    subroutine restore_gridding_pair( params, even_rec, odd_rec, state, fname_even, fname_odd, &
        &diagnostics, fsc_in, cfar_in, cones_in )
        use simple_fsc, only: fsc_area_score_result
        class(parameters),                      intent(in)    :: params
        class(reconstructor),                   intent(inout) :: even_rec, odd_rec
        integer,                                intent(in)    :: state
        class(string),                          intent(in)    :: fname_even, fname_odd
        type(gridding_pair_diagnostics),        intent(out)   :: diagnostics
        real, optional,                         intent(in)    :: fsc_in(:)
        real, optional,                         intent(in)    :: cfar_in
        class(fsc_area_score_result), optional, intent(inout) :: cones_in
        type(gridding_half_restore) :: even_restore, odd_restore
        type(fsc_area_score_result) :: cones_fsc
        real,     allocatable :: res(:)
        real                  :: smpd, fny
        integer               :: box, filtsz
        logical               :: l_have_fsc
        box    = params%box_crop
        smpd   = params%smpd_crop
        filtsz = fdim(box) - 1
        fny    = 2. * smpd
        res    = get_resarr(box, smpd)
        if( present(fsc_in) )then
            if( .not. present(cfar_in) ) THROW_HARD('cfar_in must accompany an input gridding FSC')
            if( size(fsc_in) /= filtsz ) THROW_HARD('input FSC size does not match gridding reconstruction')
            allocate(diagnostics%fsc(filtsz), source=fsc_in)
            diagnostics%cfar = cfar_in
            l_have_fsc = .true.
            if( params%l_ml_reg .and. (params%conical_fsc == 'yes') )then
                if( .not. present(cones_in) )then
                    THROW_HARD('cones_in must be provided if conical regularization is enabled')
                endif
            endif
        else
            allocate(diagnostics%fsc(filtsz), source=0.)
            l_have_fsc = .false.
        endif
        call even_restore%new(even_rec)
        call odd_restore%new(odd_rec)
        ! ML-regularization
        if( params%l_ml_reg )then
            ! preprocessing for FSC calculation
            ! even
            call even_rec%restore_base(even_restore%base, preserve_numerator=.true.)
            ! write a deapodized copy; the in-memory half stays undeapodized so the
            ! FSC below matches the legacy estimate (the inverse envelope's edge
            ! gain up-weights rim noise and shifts the FSC crossing)
            call even_restore%finalize_from_base(even_rec)
            call even_restore%final%write(add2fbody(fname_even,MRC_EXT,'_unfil'), del_if_exists=.true.)
            call even_restore%final%kill
            ! odd
            call odd_rec%restore_base(odd_restore%base, preserve_numerator=.true.)
            call odd_restore%finalize_from_base(odd_rec)
            call odd_restore%final%write(add2fbody(fname_odd,MRC_EXT,'_unfil'), del_if_exists=.true.)
            call odd_restore%final%kill
            ! Regularization
            if( l_have_fsc )then
                if( params%conical_fsc == 'yes' )then
                    call even_rec%add_conical_invtausq2rho(cones_in)
                    call odd_rec%add_conical_invtausq2rho(cones_in)
                else
                    call even_rec%add_invtausq2rho(diagnostics%fsc)
                    call odd_rec%add_invtausq2rho(diagnostics%fsc)
                endif
            else
                call calculate_gridding_pair_fsc(params, even_restore%base, odd_restore%base, &
                    &state, diagnostics%fsc, diagnostics%cfar, cones=cones_fsc)
                ! Regularization
                if( params%conical_fsc == 'yes' )then
                    call even_rec%add_conical_invtausq2rho(cones_fsc)
                    call odd_rec%add_conical_invtausq2rho(cones_fsc)
                else
                    call even_rec%add_invtausq2rho(diagnostics%fsc)
                    call odd_rec%add_invtausq2rho(diagnostics%fsc)
                endif
            endif
            ! Even: uneven sampling density correction, clip, & write
            call even_restore%base%kill
            call even_restore%prepare_final(even_rec)
            call even_rec%restore_final(even_restore%final, preserve_numerator=.true.)
            call even_restore%final%write(fname_even, del_if_exists=.true.)
            call even_restore%final%kill
            ! Odd: uneven sampling density correction, clip, & write
            call odd_restore%base%kill
            call odd_restore%prepare_final(odd_rec)
            call odd_rec%restore_final(odd_restore%final, preserve_numerator=.true.)
            call odd_restore%final%write(fname_odd, del_if_exists=.true.)
            call odd_restore%final%kill
        else
            ! correct for the uneven sampling density
            call even_rec%restore_base(even_restore%base)
            call odd_rec%restore_base(odd_restore%base)
            ! write un-normalised unmasked DEAPODIZED even/odd volumes; the in-memory
            ! halves stay undeapodized so the FSC below matches the legacy estimate
            call even_restore%finalize_from_base(even_rec)
            call even_restore%final%write(fname_even, del_if_exists=.true.)
            call even_restore%final%kill
            call odd_restore%finalize_from_base(odd_rec)
            call odd_restore%final%write(fname_odd, del_if_exists=.true.)
            call odd_restore%final%kill
            if( .not. l_have_fsc )then
                call calculate_gridding_pair_fsc(params, even_restore%base, odd_restore%base, &
                    &state, diagnostics%fsc, diagnostics%cfar, cones=cones_fsc)
            endif
        endif
        ! save, get & print resolution
        call arr2file(diagnostics%fsc, refine3D_fsc_fname(state))
        call get_resolution(diagnostics%fsc, res, diagnostics%res_fsc05, diagnostics%res_fsc0143)
        diagnostics%res_fsc05   = max(diagnostics%res_fsc05, fny)
        diagnostics%res_fsc0143 = max(diagnostics%res_fsc0143, fny)
        deallocate(res)
        call even_restore%kill
        call odd_restore%kill
        call cones_fsc%kill

    end subroutine restore_gridding_pair

    subroutine calc_env_fsc_optlp( params, even, odd, state, fsc_corrected, envmsk )
        class(parameters),     intent(in)    :: params
        class(image),          intent(inout) :: even, odd
        integer,               intent(in)    :: state
        real,     allocatable, intent(inout) :: fsc_corrected(:)
        class(image), optional,intent(inout) :: envmsk
        real, allocatable :: fsc_t(:), fsc_n(:)
        type(image)       :: mskvol
        if( allocated(fsc_corrected) ) deallocate(fsc_corrected)
        if( .not. params%l_envfsc ) return
        call calc_density_envmask(params, even, odd, mskvol)
        call mskvol%write(string(AUTOMASK_FBODY//trim(adjustl(int2str_pad(state,2)))//MRC_EXT))
        ! FSCs: phase-randomized & corrected
        call phase_rand_fsc(even, odd, mskvol, state, even%get_filtsz(), fsc_corrected, fsc_t, fsc_n)
        if( present(envmsk) ) call envmsk%copy(mskvol)
        ! cleanup
        deallocate(fsc_t, fsc_n)
        call mskvol%kill
    end subroutine calc_env_fsc_optlp

    subroutine calc_density_envmask( params, even, odd, envmsk )
        use simple_image_msk, only: image_msk
        class(parameters), intent(in)    :: params
        class(image),      intent(in)    :: even, odd
        class(image),      intent(inout) :: envmsk
        type(image)     :: merged
        type(image_msk) :: mskvol
        call merged%copy(even)
        call merged%add(odd)
        call merged%mul(0.5)
        call mskvol%automask3D(params, merged, .false., lp_override=params%envmsklp)
        call envmsk%copy(mskvol)
        call merged%kill
        call mskvol%kill_bimg
    end subroutine calc_density_envmask

    subroutine calc_masked_cfar( msk, even, odd, state, cones, cfar, envmsk )
        real,                          intent(in)    :: msk
        class(image),                  intent(in)    :: even, odd
        integer,                       intent(in)    :: state
        class(fsc_area_score_result),  intent(inout) :: cones
        real,                          intent(out)   :: cfar
        class(image), optional,        intent(in)    :: envmsk
        type(image) :: even_tmp, odd_tmp
        call even_tmp%copy(even)
        call odd_tmp%copy(odd)
        call even_tmp%ifft
        call odd_tmp%ifft
        if( present(envmsk) )then
            call even_tmp%zero_env_background(envmsk)
            call odd_tmp%zero_env_background(envmsk)
            call even_tmp%mul(envmsk)
            call odd_tmp%mul(envmsk)
        else
            call apply_spherical_fsc_mask(msk, even_tmp, odd_tmp)
        endif
        call cones%new(even_tmp, 256, 20., 0.143, 1)
        call cones%calc_fsc_area_score(even_tmp, odd_tmp, state=state)
        cfar = cones%cfar
        call even_tmp%kill
        call odd_tmp%kill
    end subroutine calc_masked_cfar

    !> Apply the broad spherical mask used when density-envelope FSC is disabled.
    subroutine apply_spherical_fsc_mask( msk, even, odd )
        real,         intent(in)    :: msk
        class(image), intent(inout) :: even, odd
        call even%mask3D_soft(msk, backgr=0.)
        call odd%mask3D_soft(msk, backgr=0.)
    end subroutine apply_spherical_fsc_mask

    !> Calculate the FSC and conical FSC area score for an explicit image pair.
    !! The input images retain the legacy representation expected by the
    !! gridding restoration path; masking and Fourier conversion use copies.
    subroutine calculate_gridding_pair_fsc( params, even, odd, state, fsc, cfar, cones )
        class(parameters),                      intent(in)    :: params
        class(image),                           intent(inout) :: even, odd
        integer,                                intent(in)    :: state
        real, allocatable,                      intent(inout) :: fsc(:)
        real,                                   intent(out)   :: cfar
        class(fsc_area_score_result), optional, intent(inout) :: cones
        type(image)                 :: even_tmp, odd_tmp, envmsk
        type(fsc_area_score_result) :: cones_local
        real    :: msk
        integer :: filtsz
        msk     = real(params%box_crop / 2) - COSMSKHALFWIDTH - 1.
        filtsz  = fdim(params%box_crop) - 1
        if( allocated(fsc) ) deallocate(fsc)
        if( params%l_envfsc )then
            call calc_env_fsc_optlp(params, even, odd, state, fsc, envmsk=envmsk)
            if( present(cones) )then
                call calc_masked_cfar(msk, even, odd, state, cones, cfar, envmsk=envmsk)
            else
                call calc_masked_cfar(msk, even, odd, state, cones_local, cfar, envmsk=envmsk)
                call cones_local%kill
            endif
            call envmsk%kill
        else
            allocate(fsc(filtsz), source=0.)
            if( present(cones) )then
                call calc_masked_cfar(msk, even, odd, state, cones, cfar)
            else
                call calc_masked_cfar(msk, even, odd, state, cones_local, cfar)
                call cones_local%kill
            endif
            ! the temporaries protect the caller's pair from the destructive
            ! spherical masking; the envfsc branch never needs them
            call even_tmp%copy(even)
            call odd_tmp%copy(odd)
            call even_tmp%ifft
            call odd_tmp%ifft
            call apply_spherical_fsc_mask(msk, even_tmp, odd_tmp)
            call even_tmp%fft
            call odd_tmp%fft
            call even_tmp%fsc(odd_tmp, fsc)
        endif
        call even_tmp%kill
        call odd_tmp%kill
    end subroutine calculate_gridding_pair_fsc

    subroutine write_gridding_pair_diagnostics( diagnostics, box, smpd, fname )
        type(gridding_pair_diagnostics), intent(in) :: diagnostics
        integer,                         intent(in) :: box
        real,                            intent(in) :: smpd
        class(string),                   intent(in) :: fname
        real, allocatable :: res(:)
        integer :: k, fnr
        if( .not. allocated(diagnostics%fsc) ) THROW_HARD('No gridding pair FSC available to write')
        res = get_resarr(box, smpd)
        call fopen(fnr, FILE=fname, STATUS='REPLACE', action='WRITE')
        do k=1,size(res)
            write(fnr,'(A,1X,F6.2,1X,A,1X,F7.3)') &
                &'>>> RESOLUTION:', res(k), '>>> CORRELATION:', diagnostics%fsc(k)
        end do
        write(fnr,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', diagnostics%res_fsc05
        write(fnr,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', diagnostics%res_fsc0143
        write(fnr,'(A,1X,F6.2)') '>>> CONICAL FSC AREA RATIO (cFAR) SCORE  :', diagnostics%cfar
        call fclose(fnr)
        deallocate(res)
    end subroutine write_gridding_pair_diagnostics

    ! DIAGNOSTIC LIFECYCLE

    subroutine kill_gridding_pair_diagnostics( self )
        class(gridding_pair_diagnostics), intent(inout) :: self
        if( allocated(self%fsc) ) deallocate(self%fsc)
        self%res_fsc05   = 0.
        self%res_fsc0143 = 0.
        self%cfar        = 0.
    end subroutine kill_gridding_pair_diagnostics

end module simple_reconstructor_eo
