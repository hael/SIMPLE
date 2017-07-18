module simple_hadamard_common
use simple_defs
use simple_cmdline,      only: cmdline
use simple_build,        only: build
use simple_params,       only: params
use simple_ori,          only: ori
use simple_oris,         only: oris
use simple_rnd,          only: ran3
use simple_gridding,     only: prep4cgrid
use simple_strings       ! use all in there
use simple_math          ! use all in there
use simple_masker        ! use all in there
implicit none

public :: read_img_from_stk, set_bp_range, set_bp_range2D, grid_ptcl, prepimg4align,&
&eonorm_struct_facts, norm_struct_facts, preprefvol, prep2Dref, preprecvols, killrecvols
private

logical, parameter :: DEBUG     = .false.
real,    parameter :: SHTHRESH  = 0.0001
real,    parameter :: CENTHRESH = 0.01   ! threshold for performing volume/cavg centering in pixels
    
contains

    subroutine read_img_from_stk( b, p, iptcl )
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        integer,        intent(in)    :: iptcl
        integer :: stack_index
        if( p%l_distr_exec )then
            stack_index = iptcl - p%fromp + 1
            call b%img%read(p%stk_part, stack_index)
        else
            call b%img%read(p%stk, iptcl)
        endif
    end subroutine read_img_from_stk
    
    subroutine set_bp_range( b, p, cline )
        use simple_math,          only: calc_fourier_index
        use simple_estimate_ssnr, only: fsc2ssnr, fsc2optlp
        use simple_filterer,      only: resample_filter
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        real, allocatable     :: resarr(:), tmp_arr(:), fom_tmp(:), res_match(:)
        real                  :: fsc0143, fsc05, mapres(p%nstates)
        integer               :: s, loc(1), lp_ind
        character(len=STDLEN) :: fsc_fname
        logical               :: fsc_bin_exists(p%nstates), all_fsc_bin_exist
        select case(p%eo)
            case('yes')
                ! check all fsc_state*.bin exist
                all_fsc_bin_exist = .true.
                fsc_bin_exists    = .false.
                do s=1,p%nstates
                    fsc_fname = adjustl('fsc_state'//int2str_pad(s,2)//'.bin')
                    if( file_exists(trim(fsc_fname)) )fsc_bin_exists( s ) = .true.
                    if( b%a%get_statepop( s )>0 .and. .not.fsc_bin_exists(s))&
                        & all_fsc_bin_exist = .false.
                enddo
                if( p%oritab.eq.'')all_fsc_bin_exist = (count(fsc_bin_exists)==p%nstates)
                ! set low-pass Fourier index limit
                if( all_fsc_bin_exist )then
                    ! we need the worst resolved fsc
                    resarr = b%img%get_res()
                    do s=1,p%nstates
                        if( fsc_bin_exists(s) )then
                            ! these are the 'classical' resolution measures
                            fsc_fname   = adjustl('fsc_state'//int2str_pad(s,2)//'.bin')
                            tmp_arr     = file2rarr(trim(fsc_fname))
                            b%fsc(s,:)  = tmp_arr(:)
                            deallocate(tmp_arr)
                            b%ssnr(s,:) = fsc2ssnr(b%fsc(s,:))
                            call get_resolution(b%fsc(s,:), resarr, fsc05, fsc0143)
                            mapres(s)   = fsc0143
                            fom_tmp = fsc2optlp(b%fsc(s,:))
                            if( p%boxmatch .eq. p%box )then
                                b%fom(s,:) = fom_tmp
                            else
                                res_match  = b%img_match%get_res()
                                b%fom(s,:) = resample_filter(fom_tmp, resarr, res_match)
                                deallocate(res_match)
                            endif
                            deallocate(fom_tmp)
                        else
                            ! empty state
                            mapres(s)   = 0.
                            b%fsc(s,:)  = 0.
                            b%fom(s,:)  = 1.
                            b%ssnr(s,:) = 0.
                        endif
                    end do
                    loc    = maxloc(mapres)
                    lp_ind = get_lplim(b%fsc(loc(1),:))
                    p%kfromto(2) = calc_fourier_index( resarr(lp_ind), p%boxmatch, p%smpd )
                    if( p%kfromto(2) == 1 )then
                        stop 'simple_math::get_lplim gives nonsensical result (==1)'
                    endif
                    if( DEBUG ) write(*,*) '*** simple_hadamard_common ***: extracted FSC info'
                else if( cline%defined('lp') )then
                    p%kfromto(2) = calc_fourier_index( p%lp, p%boxmatch, p%smpd )
                else if( cline%defined('find') )then
                    p%kfromto(2) = min(p%find,p%tofny)
                else
                    write(*,*) 'no method available for setting the low-pass limit'
                    stop 'need fsc file, lp, or find; set_bp_range; simple_hadamard_common'
                endif
                if( p%kfromto(2)-p%kfromto(1) <= 2 )then
                    write(*,*) 'fromto:', p%kfromto(1), p%kfromto(2)
                    stop 'resolution range too narrow; set_bp_range; simple_hadamard_common'
                endif
                ! lpstop overrides any other method for setting the low-pass limit
                if( cline%defined('lpstop') )then
                    p%kfromto(2) = min(p%kfromto(2), calc_fourier_index( p%lpstop, p%boxmatch, p%smpd ))
                endif
                ! set high-pass Fourier index limit
                p%kfromto(1) = max(2,calc_fourier_index( p%hp, p%boxmatch, p%smpd ))
                ! re-set the low-pass limit
                p%lp = calc_lowpass_lim( p%kfromto(2), p%boxmatch, p%smpd )
                p%lp_dyn = p%lp
                call b%a%set_all2single('lp',p%lp)
            case('no')
                ! set Fourier index range
                p%kfromto(1) = max(2, calc_fourier_index( p%hp, p%boxmatch, p%smpd ))
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
        p%kstop_grid = calc_fourier_index(p%lp_grid, p%boxmatch, p%smpd)
        if( p%kstop_grid > p%kfromto(2) ) p%kstop_grid = p%kfromto(2)
        if( DEBUG ) write(*,*) '*** simple_hadamard_common ***: did set Fourier index range'
    end subroutine set_bp_range

    subroutine set_bp_range2D( b, p, cline, which_iter, frac_srch_space )
        use simple_estimate_ssnr, only: fsc2ssnr
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        real,           intent(in)    :: frac_srch_space
        real :: lplim
        p%kfromto(1) = max(2, calc_fourier_index(p%hp, p%boxmatch, p%smpd))
        if( cline%defined('lp') )then        
            ! set Fourier index range
            p%kfromto(2) = calc_fourier_index(p%lp, p%boxmatch, p%smpd)
            p%lp_dyn     = p%lp
            call b%a%set_all2single('lp',p%lp)
        else
            ! set Fourier index range
            if( which_iter <= LPLIM1ITERBOUND )then
                lplim = p%lplims2D(1)
            else if( frac_srch_space >= FRAC_SH_LIM .and. which_iter > LPLIM3ITERBOUND )then
                lplim = p%lplims2D(3)
            else
                lplim = p%lplims2D(2)
            endif
            p%kfromto(2) = calc_fourier_index(lplim, p%boxmatch, p%smpd)
            p%lp_dyn = lplim
            call b%a%set_all2single('lp',lplim)
        endif
        if( DEBUG ) write(*,*) '*** simple_hadamard_common ***: did set Fourier index range'
    end subroutine set_bp_range2D

    !>  \brief  grids one particle image to the volume
    subroutine grid_ptcl( b, p, orientation, os, shellweights, ran_eo )
        use simple_kbinterpol, only: kbinterpol
        class(build),              intent(inout) :: b
        class(params),             intent(inout) :: p
        class(ori),                intent(inout) :: orientation
        class(oris),     optional, intent(inout) :: os
        real,            optional, intent(in)    :: shellweights(:)
        real,            optional, intent(in)    :: ran_eo
        real             :: pw, ran, w
        integer          :: jpeak, s, k, npeaks
        type(ori)        :: orisoft, o_sym
        type(kbinterpol) :: kbwin
        logical          :: softrec
        if( p%eo .eq. 'yes' )then
            kbwin = b%eorecvols(1)%get_kbwin()
        else
            kbwin = b%recvols(1)%get_kbwin()
        endif
        softrec = .false.
        npeaks  = 1
        if( present(os) )then
            softrec = .true.
            npeaks  = os%get_noris()
        endif
        pw = 1.0
        if( orientation%isthere('w') ) pw = orientation%get('w')
        if( pw > 0. )then
            ! prepare image for gridding
            ! using the uncorrected/unmodified image as input
            if( p%l_xfel )then
                call b%img%pad(b%img_pad)
            else
                call prep4cgrid(b%img, b%img_pad, p%msk, kbwin)
            endif
            if( DEBUG ) write(*,*) '*** simple_hadamard_common ***: prepared image for gridding'
            ran = ran3()
            if( present(ran_eo) )ran = ran_eo 
            orisoft = orientation
            do jpeak = 1, npeaks
                if( DEBUG ) write(*,*) '*** simple_hadamard_common ***: gridding, iteration:', jpeak
                ! get ori info
                if( softrec )then
                    orisoft =  os%get_ori(jpeak)
                    w = orisoft%get('ow')
                else
                    w = 1.
                endif
                s = nint(orisoft%get('state'))
                if( DEBUG ) write(*,*) '*** simple_hadamard_common ***: got orientation'
                if( p%frac < 0.99 ) w = w*pw
                if( w > 0. )then
                    if( p%pgrp == 'c1' )then
                        if( p%eo .eq. 'yes' )then
                            call b%eorecvols(s)%grid_fplane(orisoft, b%img_pad, pwght=w, ran=ran, shellweights=shellweights)
                        else
                            call b%recvols(s)%inout_fplane(orisoft, .true., b%img_pad, pwght=w, shellweights=shellweights)
                        endif
                    else
                        do k=1,b%se%get_nsym()
                            o_sym = b%se%apply(orisoft, k)
                            if( p%eo .eq. 'yes' )then
                                call b%eorecvols(s)%grid_fplane(o_sym, b%img_pad, pwght=w, ran=ran, shellweights=shellweights)
                            else
                                call b%recvols(s)%inout_fplane(o_sym, .true., b%img_pad, pwght=w, shellweights=shellweights)
                            endif
                        end do
                    endif
                endif
                if( DEBUG ) write(*,*) '*** simple_hadamard_common ***: gridded ptcl'
            end do
        endif
    end subroutine grid_ptcl

    !>  \brief  prepares one particle image for alignment 
    subroutine prepimg4align( b, p, o )
        use simple_ctf, only: ctf
        class(build),      intent(inout) :: b
        class(params),     intent(inout) :: p
        type(ori),         intent(inout) :: o
        type(ctf) :: tfun
        real      :: x, y, dfx, dfy, angast
        integer   :: state, cls
        if( p%l_xfel )then
            ! nothing to do 4 now
            return
        else
            x     = o%get('x')
            y     = o%get('y')
            cls   = nint(o%get('class'))
            state = nint(o%get('state'))
            ! move to Fourier space
            call b%img%fwd_ft
            ! set CTF parameters
            if( p%ctf .ne. 'no' )then
                ! we here need to re-create the CTF object as kV/cs/fraca are now per-particle params
                ! that these parameters are part of the doc is checked in the params class
                tfun = ctf(p%smpd, o%get('kv'), o%get('cs'), o%get('fraca'))
                select case(p%tfplan%mode)
                    case('astig') ! astigmatic CTF
                        dfx    = o%get('dfx')
                        dfy    = o%get('dfy')
                        angast = o%get('angast')
                    case('noastig') ! non-astigmatic CTF
                        dfx    = o%get('dfx')
                        dfy    = dfx
                        angast = 0.
                    case DEFAULT
                        write(*,*) 'Unsupported p%tfplan%mode: ', trim(p%tfplan%mode)
                        stop 'simple_hadamard_common :: prepimg4align'
                end select
            endif
            ! deal with CTF
            select case(p%ctf)
                case('mul')  ! images have been multiplied with the CTF, no CTF-dependent weighting of the correlations
                    stop 'ctf=mul is not supported; simple_hadamard_common :: prepimg4align' 
                case('no')   ! do nothing
                case('yes')  ! do nothing
                case('flip') ! flip back
                    call tfun%apply(b%img, dfx, 'flip', dfy, angast)
                case DEFAULT
                    stop 'Unsupported ctf mode; simple_hadamard_common :: prepimg4align'
            end select
            ! shift image to rotational origin
            if(abs(x) > SHTHRESH .or. abs(y) > SHTHRESH) call b%img%shift(-x, -y)
            ! back to real-space
            call b%img%bwd_ft
            ! clip image if needed
            call b%img%clip(b%img_match) ! SQUARE DIMS ASSUMED
            ! MASKING
            if( p%l_envmsk .and. p%automsk .eq. 'cavg' )then
                ! 2D adaptive cos-edge mask
                call b%mskimg%apply_adamask2ptcl_2D(b%img_match, cls)
            else if( p%l_envmsk )then
                ! 3D adaptive cos-edge mask
                call b%mskvol%apply_adamask2ptcl_3D(o, b%img_match)
            else              
                ! soft-edged mask
                if( p%l_innermsk )then
                    call b%img_match%mask(p%msk, 'soft', inner=p%inner, width=p%width)
                else
                    call b%img_match%mask(p%msk, 'soft')
                endif
            endif
            ! return in Fourier space
            call b%img_match%fwd_ft
        endif
        if( DEBUG ) write(*,*) '*** simple_hadamard_common ***: finished prepimg4align'
    end subroutine prepimg4align

    !>  \brief  prepares one cluster centre image for alignment
    subroutine prep2Dref( b, p, icls )
        use simple_image, only: image
        class(build),      intent(inout) :: b
        class(params),     intent(in)    :: p
        integer,           intent(in)    :: icls
        real :: xyz(3), sharg
        ! p%center:  user gets opportunity to turn on/off centering
        ! p%doshift: search space fraction controlled centering (or not)
        if( p%center .eq. 'yes' .and. p%doshift )then
            ! center the reference
            xyz   = b%img%center(p%cenlp, 'no', p%msk, doshift=.false.)
            sharg = arg(xyz)
            if( sharg > CENTHRESH )then
                ! apply shift and update the corresponding class parameters
                call b%img%shift(xyz(1), xyz(2))
                call b%a%add_shift2class(icls, -xyz(1:2))
            endif
        endif
        ! clip image if needed
        call b%img%clip(b%img_match)
        ! normalise
        call b%img_match%norm
        ! apply mask
        if( p%l_envmsk .and. p%automsk .eq. 'cavg' )then
            ! automasking
            call b%mskimg%apply_2Denvmask22Dref(b%img_match, icls)
            if( p%l_chunk_distr )then
                call b%img_match%write(trim(p%chunktag)//'automasked_refs'//p%ext, icls)
            else if( (p%l_distr_exec .and. p%part.eq.1) .or. (.not. p%l_distr_exec) )then
                call b%img_match%write('automasked_refs'//p%ext, icls)
            endif
        else
            ! soft masking
            if( p%l_innermsk )then
                call b%img_match%mask(p%msk, 'soft', inner=p%inner, width=p%width)
            else 
                call b%img_match%mask(p%msk, 'soft')
            endif
        endif
        ! move to Fourier space
        call b%img_match%fwd_ft
    end subroutine prep2Dref

    !>  \brief  initializes all volumes for reconstruction
    subroutine preprecvols( b, p )
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        integer :: istate
        if( p%eo .eq. 'yes' )then
            do istate = 1, p%nstates
                if( b%a%get_statepop(istate) > 0)then
                    call b%eorecvols(istate)%new(p)
                    call b%eorecvols(istate)%reset_all
                endif
            end do
        else
            do istate = 1, p%nstates
                if( b%a%get_statepop(istate) > 0)then
                    call b%recvols(istate)%new([p%boxpd, p%boxpd, p%boxpd], p%smpd, p%imgkind)
                    call b%recvols(istate)%alloc_rho(p)
                    call b%recvols(istate)%reset
                endif
            end do
        endif
    end subroutine preprecvols

    !>  \brief  destructs all volumes for reconstruction
    subroutine killrecvols( b, p )
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        integer :: istate
        if( p%eo .eq. 'yes' )then
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

    !>  \brief  prepares one volume for references extraction
    subroutine preprefvol( b, p, cline, s, doexpand )
        class(build),      intent(inout) :: b
        class(params),     intent(inout) :: p
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: s
        logical, optional, intent(in)    :: doexpand
        logical :: l_doexpand, do_center
        real    :: shvec(3)
        l_doexpand = .true.
        if( present(doexpand) ) l_doexpand = doexpand
        if( p%boxmatch < p%box )call b%vol%new([p%box,p%box,p%box],p%smpd) ! ensure correct dim
        call b%vol%read(p%vols(s), isxfel=p%l_xfel)
        call b%vol%norm ! because auto-normalisation on read is taken out
        if( p%l_xfel )then
            ! no centering
        else
            do_center = .true.
            if( p%center .eq. 'no' .or. p%nstates > 1 .or. .not. p%doshift .or.&
                p%pgrp(:1) .ne. 'c' .or. cline%defined('mskfile') ) do_center = .false.
            if( do_center )then
                shvec = b%vol%center(p%cenlp,'no',p%msk,doshift=.false.) ! find center of mass shift
                if( arg(shvec) > CENTHRESH )then
                    if( p%pgrp .ne. 'c1' ) shvec(1:2) = 0.       ! shifts only along z-axis for C2 and above
                    call b%vol%shift(shvec(1),shvec(2),shvec(3)) ! performs shift
                    ! map back to particle oritentations
                    if( cline%defined('oritab') )call b%a%map3dshift22d(-shvec(:), state=s)
                endif
            endif
        endif
        ! clip
        if( p%boxmatch < p%box )then
            call b%vol%clip_inplace([p%boxmatch,p%boxmatch,p%boxmatch]) ! SQUARE DIMS ASSUMED
        endif
        ! masking
        if( p%l_xfel )then
            ! no centering or masking
        else
            if( cline%defined('mskfile') )then
                ! mask provided
                call b%mskvol%new([p%box, p%box, p%box], p%smpd, p%imgkind)
                call b%mskvol%read(p%mskfile)
                call b%mskvol%clip_inplace([p%boxmatch,p%boxmatch,p%boxmatch])
                ! no need for the below line anymore as I (HE) removed auto-normalisation on read
                ! call b%mskvol%norm_bin ! ensures [0;1] range
                call b%vol%mul(b%mskvol)
                ! re-initialise the object for 2D envelope masking
                call b%mskvol%init_envmask2D(p%msk)
                ! don't use this for any 3D work from now on, because the soft edge is removed
                ! on initialisation
            else
                ! circular masking
                if( p%l_innermsk )then
                    call b%vol%mask(p%msk, 'soft', inner=p%inner, width=p%width)
                else
                    call b%vol%mask(p%msk, 'soft')
                endif
            endif
        endif
        ! FT volume
        call b%vol%fwd_ft
        ! expand for fast interpolation
        if( l_doexpand )call b%vol%expand_cmat     
    end subroutine preprefvol
    
    subroutine eonorm_struct_facts( b, p, res, which_iter )
        use simple_image, only: image
        class(build),      intent(inout) :: b
        class(params),     intent(inout) :: p
        real,              intent(inout) :: res
        integer, optional, intent(in)    :: which_iter
        integer               :: s
        real                  :: res05s(p%nstates), res0143s(p%nstates)
        character(len=STDLEN) :: pprocvol
        ! init
        res0143s = 0.
        res05s   = 0.
        ! cycle through states
        do s=1,p%nstates
            if( b%a%get_statepop(s) == 0 )then
                ! empty state
                if( present(which_iter) )b%fsc(s,:) = 0.
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
                call b%eorecvols(s)%sum_eos
                call b%eorecvols(s)%sampl_dens_correct_eos(s)
                call b%eorecvols(s)%sampl_dens_correct_sum(b%vol)
                call b%vol%write(p%vols(s), del_if_exists=.true.)
                ! update resolutions for local execution mode
                call b%eorecvols(s)%get_res(res05s(s), res0143s(s))
                if( present(which_iter) )then
                    ! post-process volume
                    pprocvol = add2fbody(trim(p%vols(s)), p%ext, 'pproc')
                    b%fsc(s,:) = file2rarr('fsc_state'//int2str_pad(s,2)//'.bin')
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
        if( .not. p%l_distr_exec )then
            ! set the resolution limit according to the worst resolved model
            res  = maxval(res0143s)
            p%lp = min(p%lp,max(p%lpstop,res))
        endif
    end subroutine eonorm_struct_facts
    
    subroutine norm_struct_facts( b, p, which_iter )
        class(build),      intent(inout) :: b
        class(params),     intent(inout) :: p
        integer, optional, intent(in)    :: which_iter
        integer :: s
        character(len=:), allocatable :: fbody
        character(len=STDLEN) :: pprocvol
        do s=1,p%nstates
            if( b%a%get_statepop(s) == 0 )then
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
                call b%recvols(s)%sampl_dens_correct(self_out=b%vol_pad) ! this preserves the recvol for online update
                if( p%l_xfel )then
                    ! no back transformation of the volume
                else
                    call b%vol_pad%bwd_ft
                endif
                call b%vol_pad%clip(b%vol)
                call b%vol%write(p%vols(s), del_if_exists=.true.)
                if( present(which_iter) )then
                    ! post-process volume
                    pprocvol = add2fbody(trim(p%vols(s)), p%ext, 'pproc')
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

end module simple_hadamard_common
