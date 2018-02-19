! common PRIME2D/PRIME3D routines used primarily by the Hadamard matchers
module simple_hadamard_common
#include "simple_lib.f08"
use simple_image,    only: image
use simple_cmdline,  only: cmdline
use simple_build,    only: build
use simple_params,   only: params
use simple_ori,      only: ori
use simple_oris,     only: oris
implicit none

public :: read_img, read_img_and_norm, read_imgbatch, set_bp_range, set_bp_range2D, grid_ptcl,&
&prepimg4align, eonorm_struct_facts, norm_struct_facts, cenrefvol_and_mapshifts2ptcls, preprefvol,&
&prep2Dref, gen2Dclassdoc, preprecvols, killrecvols, gen_projection_frcs, prepimgbatch, grid_ptcl_tst,&
&build_pftcc_particles
private
#include "simple_local_flags.inc"

interface read_imgbatch
    module procedure read_imgbatch_1
    module procedure read_imgbatch_2
end interface read_imgbatch

interface grid_ptcl
    module procedure grid_ptcl_1
    module procedure grid_ptcl_2
end interface grid_ptcl

real, parameter :: SHTHRESH  = 0.001
real, parameter :: CENTHRESH = 0.5    ! threshold for performing volume/cavg centering in pixels

contains

    subroutine read_img( b, p, iptcl )
        class(build),  intent(inout)  :: b
        class(params), intent(inout)  :: p
        integer,       intent(in)     :: iptcl
        character(len=:), allocatable :: stkname
        integer :: ind
        if( p%l_stktab_input )then
            call p%stkhandle%get_stkname_and_ind(iptcl, stkname, ind)
            call b%img%read(stkname, ind)
        else
            if( p%l_distr_exec )then
                call b%img%read(p%stk_part, iptcl - p%fromp + 1)
            else
                call b%img%read(p%stk, iptcl)
            endif
        endif
    end subroutine read_img

    subroutine read_img_and_norm( b, p, iptcl )
        class(build),  intent(inout)  :: b
        class(params), intent(inout)  :: p
        integer,       intent(in)     :: iptcl
        call read_img( b, p, iptcl )
        call b%img%norm()
    end subroutine read_img_and_norm

    subroutine read_imgbatch_1( b, p, fromptop, ptcl_mask )
        class(build),      intent(inout)  :: b
        class(params),     intent(inout)  :: p
        integer,           intent(in)     :: fromptop(2)
        logical, optional, intent(in)     :: ptcl_mask(p%fromp:p%top)
        character(len=:), allocatable :: stkname
        integer :: iptcl, ind_in_batch, ind_in_stk
        if( present(ptcl_mask) )then
            if( p%l_stktab_input )then
                do iptcl=fromptop(1),fromptop(2)
                    if( ptcl_mask(iptcl) )then
                        ind_in_batch = iptcl - fromptop(1) + 1
                        call p%stkhandle%get_stkname_and_ind(iptcl, stkname, ind_in_stk)
                        call b%imgbatch(ind_in_batch)%read(stkname, ind_in_stk)
                    endif
                end do
            else
                if( p%l_distr_exec )then
                    do iptcl=fromptop(1),fromptop(2)
                        if( ptcl_mask(iptcl) )then
                            ind_in_batch = iptcl - fromptop(1) + 1
                            ind_in_stk   = iptcl - p%fromp + 1
                            call b%imgbatch(ind_in_batch)%read(p%stk_part, ind_in_stk)
                        endif
                    end do
                else
                    do iptcl=fromptop(1),fromptop(2)
                        if( ptcl_mask(iptcl) )then
                            ind_in_batch = iptcl - fromptop(1) + 1
                            ind_in_stk   = iptcl
                            call b%imgbatch(ind_in_batch)%read(p%stk, ind_in_stk)
                        endif
                    end do
                endif
            endif
        else
            if( p%l_stktab_input )then
                do iptcl=fromptop(1),fromptop(2)
                    ind_in_batch = iptcl - fromptop(1) + 1
                    call p%stkhandle%get_stkname_and_ind(iptcl, stkname, ind_in_stk)
                    call b%imgbatch(ind_in_batch)%read(stkname, ind_in_stk)
                end do
            else
                if( p%l_distr_exec )then
                    do iptcl=fromptop(1),fromptop(2)
                        ind_in_batch = iptcl - fromptop(1) + 1
                        ind_in_stk   = iptcl - p%fromp + 1
                        call b%imgbatch(ind_in_batch)%read(p%stk_part, ind_in_stk)
                    end do
                else
                    do iptcl=fromptop(1),fromptop(2)
                        ind_in_batch = iptcl - fromptop(1) + 1
                        ind_in_stk   = iptcl
                        call b%imgbatch(ind_in_batch)%read(p%stk, ind_in_stk)
                    end do
                endif
            endif
        endif
    end subroutine read_imgbatch_1

    subroutine read_imgbatch_2( b, p, n, pinds, batchlims )
        class(build),  intent(inout)  :: b
        class(params), intent(inout)  :: p
        integer,       intent(in)     :: n, pinds(n), batchlims(2)
        character(len=:), allocatable :: stkname
        integer :: ind_in_stk, i, ii
        if( p%l_stktab_input )then
            do i=batchlims(1),batchlims(2)
                ii = i - batchlims(1) + 1
                call p%stkhandle%get_stkname_and_ind(pinds(i), stkname, ind_in_stk)
                call b%imgbatch(ii)%read(stkname, ind_in_stk)
            end do
        else
            if( p%l_distr_exec )then
                do i=batchlims(1),batchlims(2)
                    ii         = i - batchlims(1) + 1
                    ind_in_stk = pinds(i) - p%fromp + 1
                    call b%imgbatch(ii)%read(p%stk_part, ind_in_stk)
                end do
            else
                do i=batchlims(1),batchlims(2)
                    ii         = i - batchlims(1) + 1
                    ind_in_stk = pinds(i)
                    call b%imgbatch(ii)%read(p%stk, ind_in_stk)
                end do
            endif
        endif
    end subroutine read_imgbatch_2

    subroutine set_bp_range( b, p, cline )
        use simple_math, only: calc_fourier_index
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        real, allocatable     :: resarr(:), fsc_arr(:)
        real                  :: fsc0143, fsc05, mapres(p%nstates)
        integer               :: s, loc(1), lp_ind
        character(len=STDLEN) :: fsc_fname
        logical               :: fsc_bin_exists(p%nstates), all_fsc_bin_exist, kstop_grid_set
        kstop_grid_set = .false.
        select case(p%eo)
            case('yes','aniso')
                ! check all fsc_state*.bin exist
                all_fsc_bin_exist = .true.
                fsc_bin_exists    = .false.
                do s=1,p%nstates
                    if( str_has_substr(p%refine,'het') )then
                        fsc_fname = CLUSTER3D_FSC
                    else
                        fsc_fname = FSC_FBODY//int2str_pad(s,2)//BIN_EXT
                    endif
                    fsc_bin_exists( s ) = file_exists(trim(adjustl(fsc_fname)))
                    if( b%a%get_pop(s, 'state') > 0 .and. .not.fsc_bin_exists(s))&
                        & all_fsc_bin_exist = .false.
                enddo
                if( p%oritab .eq. '' )all_fsc_bin_exist = (count(fsc_bin_exists)==p%nstates)
                ! set low-pass Fourier index limit
                if( all_fsc_bin_exist )then
                    ! we need the worst resolved fsc
                    resarr = b%img%get_res()
                    do s=1,p%nstates
                        if( fsc_bin_exists(s) )then
                            ! these are the 'classical' resolution measures
                            if( str_has_substr(p%refine,'het') )then
                                fsc_fname = CLUSTER3D_FSC
                            else
                                fsc_fname = FSC_FBODY//int2str_pad(s,2)//BIN_EXT
                            endif
                            fsc_arr    = file2rarr(trim(adjustl(fsc_fname)))
                            b%fsc(s,:) = fsc_arr(:)
                            deallocate(fsc_arr)
                            call get_resolution(b%fsc(s,:), resarr, fsc05, fsc0143)
                            mapres(s)   = fsc0143
                        else
                            ! empty state
                            mapres(s)   = 0.
                            b%fsc(s,:)  = 0.
                        endif
                    end do
                    loc    = maxloc(mapres) ! worst resolved
                    lp_ind = get_lplim_at_corr(b%fsc(loc(1),:), p%lplim_crit)
                    ! low pass limit
                    p%kfromto(2) = calc_fourier_index(resarr(lp_ind), p%boxmatch, p%smpd)
                    if( p%kfromto(2) == 1 )then
                        stop 'simple_hadamard_common, simple_math::get_lplim gives nonsensical result (==1)'
                    endif
                    ! set highest Fourier index for coarse grid search
                    p%kstop_grid   = get_lplim_at_corr(b%fsc(loc(1),:), 0.5)
                    kstop_grid_set = .true.
                    DebugPrint ' extracted FSC info'
                else if( cline%defined('lp') )then
                    p%kfromto(2) = calc_fourier_index(p%lp, p%boxmatch, p%smpd)
                else if( cline%defined('find') )then
                    p%kfromto(2) = min(p%find,p%tofny)
                else if( b%a%isthere(p%fromp,'lp') )then
                    p%kfromto(2) = calc_fourier_index(b%a%get(p%fromp,'lp'), p%boxmatch, p%smpd)
                else
                    write(*,*) 'no method available for setting the low-pass limit'
                    stop 'need fsc file, lp, or find; set_bp_range; simple_hadamard_common'
                endif
                ! lpstop overrides any other method for setting the low-pass limit
                if( cline%defined('lpstop') )then
                    p%kfromto(2) = min(p%kfromto(2), calc_fourier_index(p%lpstop, p%boxmatch, p%smpd))
                endif
                ! set high-pass Fourier index limit
                p%kfromto(1) = max(2,calc_fourier_index( p%hp, p%boxmatch, p%smpd))
                if( p%kfromto(2)-p%kfromto(1) <= 2 )then
                    if( p%hpind_fsc > 0 )then
                        p%kfromto(2) = p%kfromto(1) + 3
                    else
                        write(*,*) 'fromto:', p%kfromto(1), p%kfromto(2)
                        stop 'resolution range too narrow; set_bp_range; simple_hadamard_common'
                    endif
                endif
                ! re-set the low-pass limit
                p%lp = calc_lowpass_lim(p%kfromto(2), p%boxmatch, p%smpd)
                p%lp_dyn = p%lp
                call b%a%set_all2single('lp',p%lp)
            case('no')
                ! set Fourier index range
                p%kfromto(1) = max(2, calc_fourier_index( p%hp, p%boxmatch, p%smpd))
                if( cline%defined('lpstop') )then
                    p%kfromto(2) = min(calc_fourier_index(p%lp, p%boxmatch, p%smpd),&
                    &calc_fourier_index(p%lpstop, p%boxmatch, p%smpd))
                else
                    p%kfromto(2) = calc_fourier_index(p%lp, p%boxmatch, p%smpd)
                endif
                p%lp_dyn = p%lp
                call b%a%set_all2single('lp',p%lp)
            case DEFAULT
                stop 'Unsupported eo flag; simple_hadamard_common'
        end select
        ! set highest Fourier index for coarse grid search
        if( .not. kstop_grid_set )        p%kstop_grid = p%kfromto(2)
        if( p%kstop_grid > p%kfromto(2) ) p%kstop_grid = p%kfromto(2)
        DebugPrint '*** simple_hadamard_common ***: did set Fourier index range'
    end subroutine set_bp_range

    subroutine set_bp_range2D( b, p, cline, which_iter, frac_srch_space )
        use simple_estimate_ssnr, only: fsc2ssnr
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        real,           intent(in)    :: frac_srch_space
        real    :: lplim
        integer :: lpstart_find
        p%kfromto(1) = max(2, calc_fourier_index(p%hp, p%boxmatch, p%smpd))
        if( cline%defined('lp') )then
            p%kfromto(2) = calc_fourier_index(p%lp, p%boxmatch, p%smpd)
            p%lp_dyn     = p%lp
            call b%a%set_all2single('lp',p%lp)
        else
            if( file_exists(p%frcs) )then
                lplim = b%projfrcs%estimate_lp_for_align()
            else
                if( which_iter <= LPLIM1ITERBOUND )then
                    lplim = p%lplims2D(1)
                else if( frac_srch_space >= FRAC_SH_LIM .and. which_iter > LPLIM3ITERBOUND )then
                    lplim = p%lplims2D(3)
                else
                    lplim = p%lplims2D(2)
                endif
            endif
            p%kfromto(2) = calc_fourier_index(lplim, p%boxmatch, p%smpd)
            ! to avoid pathological cases, fall-back on lpstart
            lpstart_find = calc_fourier_index(p%lpstart, p%boxmatch, p%smpd)
            if( lpstart_find > p%kfromto(2) ) p%kfromto(2) = lpstart_find
            p%lp_dyn = lplim
            call b%a%set_all2single('lp',lplim)
        endif
        DebugPrint  '*** simple_hadamard_common ***: did set Fourier index range'
    end subroutine set_bp_range2D

    !>  \brief  grids one particle image to the volume
    subroutine grid_ptcl_1( b, p, img, se, o )
        use simple_sym, only: sym
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        class(image),  intent(inout) :: img
        class(sym),    intent(inout) :: se
        class(ori),    intent(inout) :: o
        real      :: pw
        integer   :: s, eo
        pw = 1.0
        if( o%isthere('w') ) pw = o%get('w')
        eo = 0
        if( p%eo .ne. 'no' ) eo = nint(o%get('eo'))
        if( pw > TINY )then
            ! fwd ft
            call img%fwd_ft
            ! state flag
            s = o%get_state()
            ! gridding
            if( p%eo .ne. 'no' )then
                call b%eorecvols(s)%grid_fplane(se, o, img, eo, pwght=pw)
            else
                call b%recvols(s)%insert_fplane(se, o, img, pwght=pw)
            endif
        endif
    end subroutine grid_ptcl_1

    !>  \brief  grids one particle image to the volume (distribution of weigted oris)
    subroutine grid_ptcl_2( b, p, img, se, o, os )
        use simple_sym, only: sym
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        class(image),  intent(inout) :: img
        class(sym),    intent(inout) :: se
        class(ori),    intent(inout) :: o
        class(oris),   intent(inout) :: os
        real, allocatable :: states(:)
        real    :: pw
        integer :: s, eo
        pw = 1.0
        if( o%isthere('w') ) pw = o%get('w')
        eo = 0
        if( p%eo .ne. 'no' ) eo = nint(o%get('eo'))
        if( pw > TINY )then
            ! fwd ft
            call img%fwd_ft
            ! gridding
            if( p%nstates == 1 )then
                if( p%eo .ne. 'no' )then
                    call b%eorecvols(1)%grid_fplane(se, os, img, eo, pwght=pw)
                else
                    call b%recvols(1)%insert_fplane(se, os, img, pwght=pw)
                endif
            else
                states = os%get_all('state')
                do s=1,p%nstates
                    if( count(nint(states) == s) > 0 )then
                        if( p%eo .ne. 'no' )then
                            call b%eorecvols(s)%grid_fplane(se, os, img, eo, pwght=pw, state=s)
                        else
                            call b%recvols(s)%insert_fplane(se, os, img, pwght=pw, state=s)
                        endif
                    endif
                end do
            endif
        endif
    end subroutine grid_ptcl_2

    !>  \brief  grids one particle image to the volume
    subroutine grid_ptcl_tst( b, p, img, orientation )
        use simple_kbinterpol, only: kbinterpol
        class(build),          intent(inout) :: b
        class(params),         intent(inout) :: p
        class(image),          intent(inout) :: img
        class(ori),            intent(inout) :: orientation
        real       :: pw
        integer    :: s, npeaks, eo, nstates
        npeaks  = 1
        nstates = 1
        ! particle weight
        pw = 1.0
        if( orientation%isthere('w') ) pw = orientation%get('w')
        if( pw <= TINY ) return
        ! e/o flag
        eo = 0
        if( p%eo .ne. 'no' ) eo = nint(orientation%get('eo'))
        s  = 1
        ! fwd ft
        call img%fwd_ft
        ! one peak & one state
         call b%eorecvols(s)%grid_fplane(b%se, orientation, img, eo, pw)
    end subroutine grid_ptcl_tst

    !>  \brief  prepares all particle images for alignment
    subroutine build_pftcc_particles( b, p, pftcc, batchsz_max, match_imgs, is3D,  ptcl_mask )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_polarizer,        only: polarizer
        class(build),            intent(inout) :: b
        class(params),           intent(inout) :: p
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: batchsz_max
        class(polarizer),        intent(inout) :: match_imgs(batchsz_max)
        logical,                 intent(in)    :: is3D
        logical, optional,       intent(in)    :: ptcl_mask(p%fromp:p%top)
        logical :: mask_here(p%fromp:p%top)
        integer :: iptcl_batch, batchlims(2), imatch, iptcl
        if( present(ptcl_mask) )then
            mask_here = ptcl_mask
        else
            mask_here = .true.
        endif
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING PARTICLES'
        call prepimgbatch(b, p, batchsz_max)
        do iptcl_batch=p%fromp,p%top,batchsz_max
            batchlims = [iptcl_batch,min(p%top,iptcl_batch + batchsz_max - 1)]
            call read_imgbatch( b, p, batchlims, mask_here )
            !$omp parallel do default(shared) private(iptcl,imatch)&
            !$omp schedule(static) proc_bind(close)
            do iptcl=batchlims(1),batchlims(2)
                if( .not. mask_here(iptcl) ) cycle
                imatch = iptcl - batchlims(1) + 1
                call prepimg4align(b, p, iptcl, b%imgbatch(imatch), match_imgs(imatch), is3D=is3D)
                ! transfer to polar coordinates
                call match_imgs(imatch)%polarize(pftcc, iptcl, .true., .true.)
            end do
            !$omp end parallel do
        end do
    end subroutine build_pftcc_particles

    !>  \brief  prepares one particle image for alignment
    subroutine prepimg4align( b, p, iptcl, img_in, img_out, is3D )
        use simple_polarizer,     only: polarizer
        use simple_estimate_ssnr, only: fsc2optlp_sub
        use simple_ctf,           only: ctf
        class(build),     intent(inout) :: b
        class(params),    intent(inout) :: p
        integer,          intent(in)    :: iptcl
        class(image),     intent(inout) :: img_in
        class(polarizer), intent(inout) :: img_out
        logical,          intent(in)    :: is3D
        type(ctf) :: tfun
        real      :: frc(b%projfrcs%get_filtsz()), filter(b%projfrcs%get_filtsz())
        real      :: x, y, dfx, dfy, angast, phshift
        integer   :: ifrc
        x = b%a%get(iptcl, 'x')
        y = b%a%get(iptcl, 'y')
        ! normalise
        call img_in%norm()
        ! move to Fourier space
        call img_in%fwd_ft
        ! deal with CTF
        select case(p%ctf)
            case('mul')  ! images have been multiplied with the CTF, no CTF-dependent weighting of the correlations
                stop 'ctf=mul is not supported; simple_hadamard_common :: prepimg4align'
            case('no')   ! do nothing
            case('yes')  ! do nothing
            case('flip') ! flip back
                ! we here need to re-create the CTF object as kV/cs/fraca are now per-particle params
                ! that these parameters are part of the doc is checked in the params class
                tfun = ctf(p%smpd, b%a%get(iptcl,'kv'), b%a%get(iptcl,'cs'), b%a%get(iptcl,'fraca'))
                select case(p%tfplan%mode)
                    case('astig') ! astigmatic CTF
                        dfx    = b%a%get(iptcl,'dfx')
                        dfy    = b%a%get(iptcl,'dfy')
                        angast = b%a%get(iptcl,'angast')
                    case('noastig') ! non-astigmatic CTF
                        dfx    = b%a%get(iptcl,'dfx')
                        dfy    = dfx
                        angast = 0.
                    case DEFAULT
                        write(*,*) 'Unsupported p%tfplan%mode: ', trim(p%tfplan%mode)
                        stop 'simple_hadamard_common :: prepimg4align'
                end select
                phshift = 0.
                if( p%tfplan%l_phaseplate ) phshift = b%a%get(iptcl,'phshift')
                call tfun%apply_serial(img_in, dfx, 'flip', dfy, angast, phshift)
            case DEFAULT
                stop 'Unsupported ctf mode; simple_hadamard_common :: prepimg4align'
        end select
        ! shell normalization & filter
        if( is3D )then
            if( trim(p%eo) .ne. 'no' )then
                ifrc = b%e_bal%find_closest_proj( b%a%get_ori(iptcl) )
            else
                ifrc = 0 ! turns off the matched filter
            endif
        else
            if( trim(p%match_filt) .eq. 'yes' )then
                ifrc = nint(b%a%get(iptcl,'class'))
            else
                ifrc = 0 ! turns off the matched filter
            endif
        endif
        if( ifrc > 0 )then
            call b%projfrcs%frc_getter(ifrc, p%hpind_fsc, p%tfplan%l_phaseplate, frc)
            if( any(frc > 0.143) )then
                call fsc2optlp_sub(b%projfrcs%get_filtsz(), frc, filter)
                call img_in%shellnorm_and_apply_filter_serial(filter)
            endif
        endif
        ! shift image to rotational origin
        if(abs(x) > SHTHRESH .or. abs(y) > SHTHRESH) call img_in%shift2Dserial([-x,-y])
        ! back to real-space
        call img_in%bwd_ft
        ! clip image if needed
        call img_in%clip(img_out)
        ! soft-edged mask
        if( p%l_innermsk )then
            call img_out%mask(p%msk, 'soft', inner=p%inner, width=p%width)
        else
            if( p%l_focusmsk )then
                call img_out%mask(p%focusmsk, 'soft')
            else
                call img_out%mask(p%msk, 'soft')
            endif
        endif
        ! return in Fourier space
        call img_out%fwd_ft
        DebugPrint  '*** simple_hadamard_common ***: finished prepimg4align'
    end subroutine prepimg4align

    !>  \brief  prepares one cluster centre image for alignment
    subroutine prep2Dref( b, p, img_in, img_out, icls, center, xyz_in, xyz_out )
        use simple_estimate_ssnr, only: fsc2optlp_sub
        use simple_polarizer,     only: polarizer
        class(build),      intent(inout) :: b
        class(params),     intent(in)    :: p
        class(image),      intent(inout) :: img_in
        class(polarizer),  intent(inout) :: img_out
        integer,           intent(in)    :: icls
        logical, optional, intent(in)    :: center
        real,    optional, intent(in)    :: xyz_in(3)
        real,    optional, intent(out)   :: xyz_out(3)
        real    :: frc(b%projfrcs%get_filtsz()), filter(b%projfrcs%get_filtsz())
        real    :: xyz(3), sharg
        logical :: do_center
        ! normalise
        call img_in%norm()
        do_center = (p%center .eq. 'yes')
        ! centering only performed if p%center.eq.'yes'
        if( present(center) ) do_center = do_center .and. center
        if( do_center )then
            if( present(xyz_in) )then
                sharg = arg(xyz_in)
                if( sharg > CENTHRESH )then
                    ! apply shift and do NOT update the corresponding class parameters
                    call img_in%fwd_ft
                    call img_in%shift2Dserial(xyz_in(1:2))
                endif
            else
                xyz = img_in%center_serial(p%cenlp, p%msk)
                sharg = arg(xyz)
                if( sharg > CENTHRESH )then
                    ! apply shift and update the corresponding class parameters
                    call img_in%fwd_ft
                    call img_in%shift2Dserial(xyz(1:2))
                    call b%a%add_shift2class(icls, -xyz(1:2))
                endif
                if( present(xyz_out) ) xyz_out = xyz
            endif
        endif
        if( p%l_match_filt )then
            ! anisotropic matched filter
            call b%projfrcs%frc_getter(icls, p%hpind_fsc, p%tfplan%l_phaseplate, frc)
            if( any(frc > 0.143) )then
                call img_in%fwd_ft ! needs to be here in case the shift was never applied (above)
                call fsc2optlp_sub(b%projfrcs%get_filtsz(), frc, filter)
                call img_in%shellnorm_and_apply_filter_serial(filter)
            endif
        endif
        ! ensure we are in real-space before clipping
        call img_in%bwd_ft
        ! clip image if needed
        call img_in%clip(img_out)
        ! apply mask
        if( p%l_innermsk )then
            call img_out%mask(p%msk, 'soft', inner=p%inner, width=p%width)
        else
            call img_out%mask(p%msk, 'soft')
        endif
        ! move to Fourier space
        call img_out%fwd_ft
    end subroutine prep2Dref

    !>  \brief prepares a 2D class document with class index, resolution,
    !!         poulation, average correlation and weight
    subroutine gen2Dclassdoc( b, p, fname )
        class(build),     intent(inout) :: b
        class(params),    intent(inout) :: p
        character(len=*), intent(in)    :: fname
        integer, allocatable :: pops(:)
        integer    :: icls, pop
        real       :: frc05, frc0143
        type(oris) :: classdoc
        call classdoc%new_clean(p%ncls)
        call b%a%get_pops(pops, 'class', maxn=p%ncls)
        do icls=1,p%ncls
            call b%projfrcs%estimate_res(icls, frc05, frc0143)
            pop = pops(icls)
            call classdoc%set(icls, 'class', real(icls))
            call classdoc%set(icls, 'pop',   real(pop))
            call classdoc%set(icls, 'res',   frc0143)
            if( pop > 1 )then
                call classdoc%set(icls, 'corr',  b%a%get_avg('corr', class=icls))
                call classdoc%set(icls, 'w',     b%a%get_avg('w',    class=icls))
            else
                call classdoc%set(icls, 'corr', -1.0)
                call classdoc%set(icls, 'w',     0.0)
            endif
        end do
        call classdoc%write(fname)
        call classdoc%kill
        deallocate(pops)
    end subroutine gen2Dclassdoc

    !>  \brief  initializes all volumes for reconstruction
    subroutine preprecvols( b, p )
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        character(len=:), allocatable :: recname, rhoname, part_str
        integer :: istate
        allocate(part_str, source=int2str_pad(p%part,p%numlen))
        select case(p%eo)
            case('yes','aniso')
                do istate = 1, p%nstates
                    if( b%a%get_pop(istate, 'state') > 0)then
                        call b%eorecvols(istate)%new(p)
                        call b%eorecvols(istate)%reset_all
                        if( p%l_frac_update )then
                            call b%eorecvols(istate)%read_eos(VOL_FBODY//int2str_pad(istate,2)//'_part'//part_str)
                            call b%eorecvols(istate)%expand_exp
                            call b%eorecvols(istate)%apply_weight(1.0 - p%update_frac)
                        endif
                    endif
                end do
            case DEFAULT
                do istate = 1, p%nstates
                    if( b%a%get_pop(istate, 'state') > 0)then
                        call b%recvols(istate)%new([p%boxpd, p%boxpd, p%boxpd], p%smpd)
                        call b%recvols(istate)%alloc_rho(p)
                        call b%recvols(istate)%reset
                        call b%recvols(istate)%reset_exp
                        if( p%l_frac_update )then
                            allocate(recname, source=VOL_FBODY//int2str_pad(istate,2)//'_part'//part_str//p%ext)
                            allocate(rhoname, source='rho_'//VOL_FBODY//int2str_pad(istate,2)//'_part'//part_str//p%ext)
                            if( file_exists(recname) .and. file_exists(rhoname) )then
                                call b%recvols(istate)%read(recname)
                                call b%recvols(istate)%read_rho(rhoname)
                                call b%recvols(istate)%expand_exp
                                call b%recvols(istate)%apply_weight(1.0 - p%update_frac)
                            endif
                            deallocate(recname, rhoname)
                        endif
                    endif
                end do
        end select
    end subroutine preprecvols

    !>  \brief  destructs all volumes for reconstruction
    subroutine killrecvols( b, p )
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        integer :: istate
        if( p%eo .ne. 'no' )then
            do istate = 1, p%nstates
                call b%eorecvols(istate)%kill
            end do
        else
            do istate = 1, p%nstates
                call b%recvols(istate)%dealloc_rho
                call b%recvols(istate)%kill
            end do
        endif
    end subroutine killrecvols

    !>  \brief  prepares a batch of images
    subroutine prepimgbatch( b, p, batchsz )
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        integer,        intent(in)    :: batchsz
        integer :: currsz, ibatch
        logical :: doprep
        if( .not. allocated(b%imgbatch) )then
            doprep = .true.
        else
            currsz = size(b%imgbatch)
            if( batchsz > currsz )then
                do ibatch=1,currsz
                    call b%imgbatch(ibatch)%kill
                end do
                deallocate(b%imgbatch)
                doprep = .true.
            else
                doprep = .false.
            endif
        endif
        if( doprep )then
            allocate(b%imgbatch(batchsz))
            do ibatch=1,batchsz
                call b%imgbatch(ibatch)%new([p%box,p%box,1], p%smpd)
            end do
        endif
    end subroutine prepimgbatch

    !>  \brief  center the reference volume and map shifts back to particles
    subroutine cenrefvol_and_mapshifts2ptcls( b, p, cline, s, volfname, do_center, xyz )
        class(build),     intent(inout) :: b
        class(params),    intent(inout) :: p
        class(cmdline),   intent(inout) :: cline
        integer,          intent(in)    :: s
        character(len=*), intent(in)    :: volfname
        logical,          intent(out)   :: do_center
        real,             intent(out)   :: xyz(3)
        ! ensure correct b%vol dim
        call b%vol%new([p%box,p%box,p%box],p%smpd)
        ! centering
        do_center = .true.
        if( p%center .eq. 'no' .or. p%nstates > 1 .or. .not. p%l_doshift .or.&
        &p%pgrp(:1) .ne. 'c' .or. cline%defined('mskfile') .or. p%l_frac_update )then
            do_center = .false.
            xyz       = 0.
            return
        endif
        call b%vol%read(volfname)
        call b%vol%norm() ! because auto-normalisation on read is taken out
        xyz = b%vol%center(p%cenlp,p%msk) ! find center of mass shift
        if( arg(xyz) <= CENTHRESH )then
            do_center = .false.
            xyz = 0.
            return
        endif
        call b%vol%fwd_ft
        if( p%pgrp .ne. 'c1' ) xyz(1:2) = 0.     ! shifts only along z-axis for C2 and above
        call b%vol%shift([xyz(1),xyz(2),xyz(3)]) ! performs shift
        ! map back to particle oritentations
        if( cline%defined('oritab') ) call b%a%map3dshift22d(-xyz(:), state=s)
    end subroutine cenrefvol_and_mapshifts2ptcls

    !>  \brief  prepares one volume for references extraction
    subroutine preprefvol( b, p, cline, s, volfname, do_center, xyz )
        use simple_estimate_ssnr, only: fsc2optlp
        class(build),     intent(inout) :: b
        class(params),    intent(inout) :: p
        class(cmdline),   intent(inout) :: cline
        integer,          intent(in)    :: s
        character(len=*), intent(in)    :: volfname
        logical,          intent(in)    :: do_center
        real,             intent(in)    :: xyz(3)
        type(image)                     :: mskvol
        real,             allocatable   :: filter(:)
        character(len=:), allocatable   :: fname_vol_filter
        ! ensure correct b%vol dim
        call b%vol%new([p%box,p%box,p%box],p%smpd)
        call b%vol%read(volfname)
        call b%vol%norm() ! because auto-normalisation on read is taken out
        if( do_center )then
            call b%vol%fwd_ft
            call b%vol%shift([xyz(1),xyz(2),xyz(3)])
        endif
        ! Volume filtering
        if( p%eo.ne.'no' )then
            ! anisotropic matched filter
            if( p%nstates.eq.1 )then
                allocate(fname_vol_filter, source=ANISOLP_FBODY//int2str_pad(s,2)//trim(p%ext))
            else
                allocate(fname_vol_filter, source=CLUSTER3D_ANISOLP//trim(p%ext))
            endif
            if( file_exists(fname_vol_filter) )then
                call b%vol2%read(fname_vol_filter)
                call b%vol%fwd_ft ! needs to be here in case the shift was never applied (above)
                call b%vol%shellnorm_and_apply_filter(b%vol2)
            else
                ! matched filter based on Rosenthal & Henderson, 2003
                if( any(b%fsc(s,:) > 0.143) )then
                    call b%vol%fwd_ft ! needs to be here in case the shift was never applied (above)
                    filter = fsc2optlp(b%fsc(s,:))
                    call b%vol%shellnorm_and_apply_filter(filter)
                endif
            endif
            deallocate(fname_vol_filter)
        endif
        ! back to real space
        call b%vol%bwd_ft
        ! clip
        if( p%boxmatch < p%box ) call b%vol%clip_inplace([p%boxmatch,p%boxmatch,p%boxmatch])
        ! masking
        if( cline%defined('mskfile') )then
            ! mask provided
            call mskvol%new([p%box, p%box, p%box], p%smpd)
            call mskvol%read(p%mskfile)
            if( p%boxmatch < p%box ) call mskvol%clip_inplace([p%boxmatch,p%boxmatch,p%boxmatch])
            call b%vol%zero_background
            call b%vol%mul(mskvol)
            call mskvol%kill
        else
            ! circular masking
            if( p%l_innermsk )then
                call b%vol%mask(p%msk, 'soft', inner=p%inner, width=p%width)
            else
                call b%vol%mask(p%msk, 'soft')
            endif
        endif
        ! FT volume
        call b%vol%fwd_ft
        ! expand for fast interpolation
        call b%vol%expand_cmat(p%alpha)
    end subroutine preprefvol

    subroutine norm_struct_facts( b, p, which_iter )
        class(build),      intent(inout) :: b
        class(params),     intent(inout) :: p
        integer, optional, intent(in)    :: which_iter
        integer :: s
        character(len=:), allocatable :: fbody
        character(len=STDLEN) :: pprocvol
        do s=1,p%nstates
            if( b%a%get_pop(s, 'state') == 0 )then
                ! empty space
                cycle
            endif
            if( p%l_distr_exec )then
                allocate(fbody, source='recvol_state'//int2str_pad(s,2)//'_part'//int2str_pad(p%part,p%numlen))
                p%vols(s)  = trim(adjustl(fbody))//p%ext
                call b%recvols(s)%compress_exp
                call b%recvols(s)%write(p%vols(s), del_if_exists=.true.)
                call b%recvols(s)%write_rho('rho_'//trim(adjustl(fbody))//p%ext)
                deallocate(fbody)
            else
                if( p%refine .eq. 'snhc' )then
                     p%vols(s) = trim(SNHCVOL)//int2str_pad(s,2)//p%ext
                else
                    if( present(which_iter) )then
                        p%vols(s) = 'recvol_state'//int2str_pad(s,2)//'_iter'//int2str_pad(which_iter,3)//p%ext
                    else
                        p%vols(s) = 'startvol_state'//int2str_pad(s,2)//p%ext
                    endif
                endif
                call b%recvols(s)%compress_exp
                call b%recvols(s)%sampl_dens_correct
                call b%recvols(s)%bwd_ft
                call b%recvols(s)%clip(b%vol)
                call b%vol%write(p%vols(s), del_if_exists=.true.)
                if( present(which_iter) )then
                    ! post-process volume
                    pprocvol = add2fbody(trim(p%vols(s)), p%ext, '_pproc')
                    call b%vol%fwd_ft
                    ! low-pass filter
                    call b%vol%bp(0., p%lp)
                    call b%vol%bwd_ft
                    ! mask
                    call b%vol%mask(p%msk, 'soft')
                    call b%vol%write(pprocvol)
                endif
            endif
        end do
    end subroutine norm_struct_facts

    subroutine eonorm_struct_facts( b, p, cline, res, which_iter )
        use simple_filterer, only: gen_anisotropic_optlp
        class(build),      intent(inout) :: b
        class(params),     intent(inout) :: p
        class(cmdline),    intent(inout) :: cline
        real,              intent(inout) :: res
        integer, optional, intent(in)    :: which_iter
        integer               :: s, find4eoavg
        real                  :: res05s(p%nstates), res0143s(p%nstates)
        character(len=STDLEN) :: pprocvol
        character(len=32)     :: resmskname
        ! init
        res0143s = 0.
        res05s   = 0.
        ! cycle through states
        do s=1,p%nstates
            if( b%a%get_pop(s, 'state') == 0 )then
                ! empty state
                if( present(which_iter) ) b%fsc(s,:) = 0.
                cycle
            endif
            call b%eorecvols(s)%compress_exp
            if( p%l_distr_exec )then
                call b%eorecvols(s)%write_eos('recvol_state'//int2str_pad(s,2)//'_part'//int2str_pad(p%part,p%numlen))
            else
                if( present(which_iter) )then
                    p%vols(s) = 'recvol_state'//int2str_pad(s,2)//'_iter'//int2str_pad(which_iter,3)//p%ext
                else
                    p%vols(s) = 'startvol_state'//int2str_pad(s,2)//p%ext
                endif
                p%vols_even(s) = add2fbody(p%vols(s), p%ext, '_even')
                p%vols_odd(s)  = add2fbody(p%vols(s), p%ext, '_odd')
                resmskname     = 'resmask'//p%ext
                call b%eorecvols(s)%sum_eos
                call b%eorecvols(s)%sampl_dens_correct_eos(s, p%vols_even(s), p%vols_odd(s), resmskname, find4eoavg)
                call gen_projection_frcs( b, p, cline, p%vols_even(s), p%vols_odd(s), resmskname, s, b%projfrcs)
                call b%projfrcs%write('frcs_state'//int2str_pad(s,2)//'.bin')
                call gen_anisotropic_optlp(b%vol2, b%projfrcs, b%e_bal, s, p%pgrp, p%hpind_fsc, p%tfplan%l_phaseplate)
                call b%vol2%write('aniso_optlp_state'//int2str_pad(s,2)//p%ext)
                call b%eorecvols(s)%get_res(res05s(s), res0143s(s))
                call b%eorecvols(s)%sampl_dens_correct_sum(b%vol)
                call b%vol%write(p%vols(s), del_if_exists=.true.)
                 ! need to put the sum back at lowres for the eo pairs
                call b%vol%fwd_ft
                call b%vol2%zero_and_unflag_ft
                call b%vol2%read(p%vols_even(s))
                call b%vol2%fwd_ft
                call b%vol2%insert_lowres(b%vol, find4eoavg)
                call b%vol2%bwd_ft
                call b%vol2%write(p%vols_even(s), del_if_exists=.true.)
                call b%vol2%zero_and_unflag_ft
                call b%vol2%read(p%vols_odd(s))
                call b%vol2%fwd_ft
                call b%vol2%insert_lowres(b%vol, find4eoavg)
                call b%vol2%bwd_ft
                call b%vol2%write(p%vols_odd(s), del_if_exists=.true.)
                if( present(which_iter) )then
                    ! post-process volume
                    pprocvol   = add2fbody(trim(p%vols(s)), p%ext, '_pproc')
                    b%fsc(s,:) = file2rarr('fsc_state'//int2str_pad(s,2)//'.bin')
                    ! low-pass filter
                    call b%vol%bp(0., p%lp)
                    call b%vol%bwd_ft
                    ! mask
                    call b%vol%mask(p%msk, 'soft')
                    call b%vol%write(pprocvol)
                else
                    call b%vol%zero_and_unflag_ft
                endif
            endif
        end do
        if( .not. p%l_distr_exec )then
            ! set the resolution limit according to the worst resolved model
            res  = maxval(res0143s)
            p%lp = min(p%lp,max(p%lpstop,res))
        endif
    end subroutine eonorm_struct_facts

    !>  \brief generate projection FRCs from even/odd pairs
    subroutine gen_projection_frcs( b, p, cline, ename, oname, resmskname, state, projfrcs )
        use simple_oris,            only: oris
        use simple_projector_hlev,  only: project
        use simple_projection_frcs, only: projection_frcs
        class(build),           intent(inout) :: b
        class(params),          intent(inout) :: p
        class(cmdline),         intent(inout) :: cline
        character(len=*),       intent(in)    :: ename, oname, resmskname
        integer,                intent(in)    :: state
        class(projection_frcs), intent(inout) :: projfrcs
        type(oris)               :: e_space
        type(image)              :: mskvol
        type(image), allocatable :: even_imgs(:), odd_imgs(:)
        real,        allocatable :: frc(:)
        integer :: iproj, find_plate
        ! ensure correct b%vol dim
        call b%vol%new([p%box,p%box,p%box],p%smpd)
        ! read & prep even/odd pair
        call b%vol%read(ename)
        call b%vol2%read(oname)
        call mskvol%new([p%box, p%box, p%box], p%smpd)
        call prepeovol(b%vol)
        call prepeovol(b%vol2)
        ! create e_space
        call e_space%new(NSPACE_BALANCE)
        call e_space%spiral(p%nsym, p%eullims)
        ! generate even/odd projections
        even_imgs = project(b%vol,  e_space, p)
        odd_imgs  = project(b%vol2, e_space, p)
        ! calculate FRCs and fill-in projfrcs object
        allocate(frc(even_imgs(1)%get_filtsz()))
        !$omp parallel do default(shared) private(iproj,frc) schedule(static) proc_bind(close)
        do iproj=1,NSPACE_BALANCE
            call even_imgs(iproj)%fwd_ft
            call odd_imgs(iproj)%fwd_ft
            call even_imgs(iproj)%fsc(odd_imgs(iproj), frc)
            if( p%tfplan%l_phaseplate ) call phaseplate_correct_fsc(frc, find_plate)
            call projfrcs%set_frc(iproj, frc, state)
            call even_imgs(iproj)%kill
            call odd_imgs(iproj)%kill
        end do
        !$omp end parallel do
        deallocate(even_imgs, odd_imgs)
        call e_space%kill
        call mskvol%kill

        contains

            !>  \brief  prepares even/odd volume for FSC/FRC calcualtion
            subroutine prepeovol( vol )
                class(image), intent(inout) :: vol
                call vol%norm() ! because auto-normalisation on read is taken out
                ! masking
                if( cline%defined('mskfile') )then
                    ! mask provided
                    call mskvol%read(resmskname)
                    call vol%zero_background
                    call vol%mul(mskvol)
                else
                    ! circular masking
                    if( p%l_innermsk )then
                        call vol%mask(p%msk, 'soft', inner=p%inner, width=p%width)
                    else
                        call vol%mask(p%msk, 'soft')
                    endif
                endif
        end subroutine prepeovol

    end subroutine gen_projection_frcs

end module simple_hadamard_common
