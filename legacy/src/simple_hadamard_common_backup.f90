module simple_hadamard_common
use simple_cmdline       ! singleton
use simple_jiffys        ! singleton
use simple_math          ! singleton
use simple_defs          ! singleton
use simple_masker        ! singleton
use simple_build,        only: build
use simple_params,       only: params
use simple_ori,          only: ori
use simple_rnd,          only: ran3
use simple_prime3D_srch, only: prime3D_srch
use simple_gridding,     only: prep4cgrid
implicit none

public :: set_bp_range, grid_ptcl, prepimg4align, eonorm_struct_facts, norm_struct_facts, preprefs4align
private

logical, parameter :: debug=.false.
real, parameter    :: SHTHRESH=0.0001
real               :: res, dfx_prev=0., dfy_prev=0., angast_prev=0.
    
contains
    
    subroutine set_bp_range( b, p, cline )
        use simple_estimate_ssnr, only: fsc2ssnr, estimate_pssnr2D, estimate_pssnr3D
        use simple_cmdline,       only: cmdline
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        real, allocatable :: resarr(:), tmparr(:)
        real              :: fsc0143, fsc05, mapres(p%nstates)
        integer :: s, loc(1)
        select case(p%eo)
            case('yes')
                ! set low-pass Fourier index limit
                if( file_exists('fsc_state01.bin') )then
                    ! we need the worst resolved fsc
                    resarr = b%img%get_res()
                    do s=1,p%nstates
                        ! these are the 'classical' resolution measures
                        b%fsc(s,:)  = file2rarr('fsc_state'//int2str_pad(s,2)//'.bin')
                        b%ssnr(s,:) = fsc2ssnr(b%fsc(s,:))
                        call get_resolution(b%fsc(s,:), resarr, fsc05, fsc0143)
                        mapres(s)   = fsc0143
                        if( p%ctf .ne. 'no' )then
                            ! get the ctfsq spectrum from the previous 3D reconstruction
                            tmparr = file2rarr('ctfsqspec_state'//int2str_pad(s,2)//'.bin')
                            ! calculate the CTF**2-dependent component of the PSSNR
                            where( tmparr > 1e-6 )
                                b%pssnr_ctfsq3D(s,:) = b%spec_count3D/tmparr
                            else where
                                b%pssnr_ctfsq3D(s,:) = 0.
                            end where
                            deallocate(tmparr)
                        else
                            b%pssnr_ctfsq3D(s,:) = 1.
                        endif
                        select case(p%filter)
                            case('pssnr')
                                ! calculate the PSSNRS used to design filters in 2D (particles) and 3D (volumes)
                                b%pssnr2D(s,:) = estimate_pssnr2D(p%avr, b%fsc(s,:))*b%pssnr_ctfsq3D(s,:)
                                b%pssnr3D(s,:) = estimate_pssnr3D(p%avr, b%fsc(s,:))*b%pssnr_ctfsq3D(s,:)
                            case DEFAULT
                                b%pssnr2D(s,:) = 1.
                                b%pssnr3D(s,:) = fsc2ssnr(b%fsc(s,:))     
                        end select
                        call arr2file(b%pssnr2D(s,:), 'pssnr2D_state'//int2str_pad(s,2)//'.bin')
                        call arr2file(b%pssnr3D(s,:), 'pssnr3D_state'//int2str_pad(s,2)//'.bin')
                    end do
                    loc = maxloc(mapres)
                    p%kfromto(2) = get_lplim(b%fsc(loc(1),:))
                    if( p%kfromto(2) == 1 )then
                        stop 'simple_math::get_lplim gives nonsensical result (==1)'
                    endif
                    if( debug ) write(*,*) '*** simple_hadamard_common ***: extracted FSC info'
                else if( cline%defined('lp') )then
                    p%kfromto(2) = b%img%get_find(p%lp)
                else if( cline%defined('find') )then
                    p%kfromto(2) = p%find
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
                    p%kfromto(2) = min(p%kfromto(2),b%img%get_find(p%lpstop))
                endif
                ! set high-pass Fourier index limit
                p%kfromto(1) = max(2,b%img%get_find(p%hp))
                ! re-set the low-pass limit
                p%lp = b%img%get_lp(p%kfromto(2))
                p%lp_dyn = p%lp
                call b%a%set_all('lp',p%lp)
            case('no')
                ! set Fourier index range
                p%kfromto(1) = max(2,b%img%get_find(p%hp))
                if( cline%defined('lpstop') )then
                    p%kfromto(2) = min(b%img%get_find(p%lp),b%img%get_find(p%lpstop))
                else
                    p%kfromto(2) = b%img%get_find(p%lp)
                endif
                p%lp_dyn = p%lp
                call b%a%set_all('lp',p%lp)
            case DEFAULT
                stop 'Unsupported eo flag; simple_hadamard_common'
        end select
        if( debug ) write(*,*) '*** simple_hadamard_common ***: did set Fourier index range'
    end subroutine

    !>  \brief  grids one particle image to the volume
    subroutine grid_ptcl( b, p, iptcl, cnt_glob, orientation, primesrch3D )
        class(build), intent(inout)                  :: b
        class(params), intent(inout)                 :: p
        integer, intent(in)                          :: iptcl, cnt_glob
        class(ori), intent(inout)                    :: orientation
        class(prime3D_srch), intent(inout), optional :: primesrch3D
        real      :: pw, ran, w
        integer   :: jpeak, s, k, nbetter
        type(ori) :: orisoft, o_sym
        pw = orientation%get('w')
        if( pw > 0. )then
            if( p%l_distr_exec )then
                call b%img_copy%read(p%stk_part, cnt_glob, p%l_xfel)
            else
                call b%img_copy%read(p%stk, iptcl, p%l_xfel)
            endif
            ! prepare image for gridding
            ! using the uncorrected/unmodified image as input
            if( p%l_xfel )then
                call b%img_copy%pad(b%img_pad)
            else
                if( p%eo .eq. 'yes' )then
                    call prep4cgrid(b%img_copy, b%img_pad, p%msk, b%eorecvols(1)%get_wfuns())
                else
                    call prep4cgrid(b%img_copy, b%img_pad, p%msk, b%recvols(1)%get_wfuns())
                endif
            endif
            if( debug ) write(*,*) '*** simple_hadamard_common ***: prepared image for gridding'
            ran = ran3()
            orisoft = orientation
            do jpeak=1,p%npeaks
                if( debug ) write(*,*) '*** simple_hadamard_common ***: gridding, iteration:', jpeak
                ! get ori info
                select case(p%refine)
                    case('no', 'qcont', 'qcontneigh')
                        if( present(primesrch3D) )then
                            call primesrch3D%get_ori(jpeak, orisoft)
                            w = orisoft%get('ow')
                            s = nint(orisoft%get('state'))
                        else
                            stop 'need primesrch3D object input for this mode of execution; simple_hadamard_common :: grid_ptcl'
                        endif
                    case DEFAULT
                        orisoft = b%a%get_ori(iptcl)
                        w = 1.
                        s = nint(orisoft%get('state'))
                end select
                if( debug ) write(*,*) '*** simple_hadamard_common ***: got orientation'
                if( p%frac < 0.99 ) w = w*pw
                if( w > 0. )then
                    ! grid
                    if( p%pgrp == 'c1' )then
                        if( p%eo .eq. 'yes' )then
                            call b%eorecvols(s)%grid_fplane(orisoft, b%img_pad, pwght=w, ran=ran)
                        else
                            call b%recvols(s)%inout_fplane(orisoft, .true., b%img_pad, pwght=w)
                        endif
                    else
                        do k=1,b%se%get_nsym()
                            o_sym = b%se%apply(orisoft, k)
                            if( p%eo .eq. 'yes' )then
                                call b%eorecvols(s)%grid_fplane(o_sym, b%img_pad, pwght=w, ran=ran)
                            else
                                call b%recvols(s)%inout_fplane(o_sym, .true., b%img_pad, pwght=w)
                            endif
                        end do
                    endif
                endif
                if( debug ) write(*,*) '*** simple_hadamard_common ***: gridded ptcl'
            end do
            call orisoft%kill
        endif
    end subroutine
    
    subroutine prepimg4align( b, p, iptcl, pssnr )
        class(build), intent(inout)  :: b
        class(params), intent(inout) :: p
        integer, intent(in)          :: iptcl
        real, intent(in), optional   :: pssnr(:)
        real                         :: x, y, dfx, dfy, angast
        integer                      :: icls, state
        type(ori)                    :: o
        if( p%l_xfel )then
            ! nothing to do 4 now
            return
        else
            ! parse ori
            x     = b%a%get(iptcl, 'x')
            y     = b%a%get(iptcl, 'y')
            state = nint(b%a%get(iptcl, 'state'))
            icls  = nint(b%a%get(iptcl, 'class'))
            if( p%doautomsk )then
                ! make sure that img_msk is of the correct dimension
                call b%img_msk%new([p%boxmatch,p%boxmatch,1],p%smpd)
                ! create 2D envelope
                o = b%a%get_ori(iptcl)
                call b%proj%fproject(b%vol_pad, o, b%img_pad)
                call b%img_pad%bwd_ft
                call b%img_pad%clip(b%img_msk)
                call b%img_msk%norm('sigm')
            endif
            ! move to Fourier space
            call b%img%fwd_ft
            ! set CTF parameters
            if( p%ctf .ne. 'no' )then
                select case(p%ctfmode)
                    case('astig') ! astigmatic CTF
                        dfx    = b%a%get(iptcl,'dfx')
                        dfy    = b%a%get(iptcl,'dfy')
                        angast = b%a%get(iptcl,'angast')
                    case('noastig') ! non-astigmatic CTF
                        dfx    = b%a%get(iptcl,'dfx')
                        dfy    = dfx
                        angast = 0.
                    case DEFAULT
                        write(*,*) 'Unsupported p%ctfmode: ', trim(p%ctfmode)
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
                    call b%tfun%apply(b%img, dfx, 'flip', dfy, angast)
                case DEFAULT
                    stop 'Unsupported ctf mode; simple_hadamard_common :: prepimg4align'
            end select
            ! SSNR scaling of amplitudes
!             if( present(pssnr) )then
!                 call ssnr_scale_amps(pssnr)
!             else if( p%eo .ne. 'no' )then
!                 call ssnr_scale_amps(b%ssnr(state,:))
!             endif
            ! shift image to rotational origin
            if( abs(x) > SHTHRESH .or. abs(y) > SHTHRESH ) call b%img%shift(-x, -y)
            ! back to real-space
            call b%img%bwd_ft
            ! clip image if needed
            if( p%boxmatch < p%box ) call b%img%clip_inplace([p%boxmatch,p%boxmatch,1]) ! SQUARE DIMS ASSUMED
            ! apply a soft-edged mask
            if( p%l_innermsk )then
                call b%img%mask(p%msk, 'soft', inner=p%inner, width=p%width)
            else 
                call b%img%mask(p%msk, 'soft')
            endif
            if( p%doautomsk )then
                ! multiply with the projected envelope
                call b%img%mul(b%img_msk)
            else if( p%automsk .eq. 'cavg' )then
                ! ab initio mask
                call automask2D(b%img, p)
            endif  
            ! return in Fourier space
            call b%img%fwd_ft
        endif
        if( debug ) write(*,*) '*** simple_hadamard_common ***: finished prepimg4align'
        
      contains
          
!         subroutine ssnr_scale_amps( ssnr )
!             real, intent(in) :: ssnr(:)
!             call b%tfun%ctf2img(b%img_filt, dfx, 'square', dfy, angast)
!             call b%img_filt%apply_filter(ssnr)
!             call b%img_filt%add(1.0)
!             call b%img_filt%square_root
!             call b%img%shellnorm
!             call b%img%apply_filter(b%img_filt)
!         end subroutine

    end subroutine
    
    subroutine eonorm_struct_facts( b, p, res, which_iter )
        use simple_image, only: image
        class(build), intent(inout)    :: b
        class(params), intent(inout)   :: p
        real, intent(inout)            :: res
        integer, intent(in), optional  :: which_iter
        type(image)                    :: vol_tmp
        integer                        :: s, s_loc
        real                           :: res05s(p%nstates), res0143s(p%nstates)
        do s=1,p%nstates
            if( p%l_distr_exec )then
                call b%eorecvols(s)%write_eos('recvol'//'_state'//int2str_pad(s,2)//'_part'//int2str_pad(p%part,p%numlen))
            else
                if( present(which_iter) )then
                    if( which_iter <= 0 )then
                        p%vols(s) = 'recvol'//'_state'//int2str_pad(s,2)//p%ext
                    else
                        p%vols(s) = 'recvol'//'_state'//int2str_pad(s,2)//'_iter'//int2str_pad(which_iter,3)//p%ext
                    endif
                else
                     p%vols(s) = 'startvol'//'_state'//int2str_pad(s,2)//p%ext
                endif
                call b%eorecvols(s)%sum_eos
                call b%eorecvols(s)%sampl_dens_correct_eos(s)
                call b%eorecvols(s)%sampl_dens_correct_sum(b%vol)
                call b%vol%write(p%vols(s), del_if_exists=.true.)
            endif
        end do
        if( .not. p%l_distr_exec )then
            ! set the resolution limit according to the worst resolved model
            res  = maxval(res0143s)
            p%lp = min(p%lp,max(p%lpstop,res))
        endif
    end subroutine
    
    subroutine norm_struct_facts( b, p, which_iter )
        class(build), intent(inout)   :: b
        class(params), intent(inout)  :: p
        integer, intent(in), optional :: which_iter
        integer                       :: s
        character(len=:), allocatable :: fbody
        do s=1,p%nstates
            if( p%l_distr_exec )then
                allocate(fbody, source='recvol_state'//int2str_pad(s,2)//'_part'//int2str_pad(p%part,p%numlen))
                p%vols(s)  = trim(adjustl(fbody))//p%ext
                p%masks(s) = 'rho_'//trim(adjustl(fbody))//p%ext
                call b%recvols(s)%write(p%vols(s), del_if_exists=.true.)
                call b%recvols(s)%write_rho(p%masks(s))
                deallocate(fbody)
            else
                if( present(which_iter) )then
                    if( which_iter <= 0 )then
                        p%vols(s) = 'recvol_state'//int2str_pad(s,2)//p%ext
                    else
                        p%vols(s) = 'recvol_state'//int2str_pad(s,2)//'_iter'//int2str_pad(which_iter,3)//p%ext
                    endif
                else
                     p%vols(s) = 'startvol_state'//int2str_pad(s,2)//p%ext
                endif
                call b%recvols(s)%sampl_dens_correct(self_out=b%vol_pad) ! this preserves the recvol for online update
                if( p%l_xfel )then
                    ! no back transformation of the volume
                else
                    call b%vol_pad%bwd_ft
                endif
                call b%vol_pad%clip(b%vol)
                call b%vol%write(p%vols(s), del_if_exists=.true.)
            endif
        end do
    end subroutine
    
    subroutine preprefs4align( b, p, iptcl, pftcc )
        use simple_math,             only: euclid
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(build),            intent(inout) :: b
        class(params),           intent(inout) :: p
        integer,                 intent(in)    :: iptcl
        class(polarft_corrcalc), intent(inout) :: pftcc
        real :: dfx, dfy, angast, dist
        if( p%ctf .ne. 'no' )then
            if( p%ctfmode .eq. 'astig' )then ! astigmatic CTF
                dfx = b%a%get(iptcl,'dfx')
                dfy = b%a%get(iptcl,'dfy')
                angast = b%a%get(iptcl,'angast')
            else if( p%ctfmode .eq. 'noastig' )then
                dfx = b%a%get(iptcl,'dfx')
                dfy = dfx
                angast = 0.
            else
                stop 'Unsupported ctf mode; preprefs4align; simple_hadamard_common'
            endif
            dist = euclid([dfx,dfy,angast],[dfx_prev,dfy_prev,angast_prev])
            if( dist < 0.001 )then
                ! CTF parameters are the same as for the previous particle & no update is needed
            else
                ! CTF parameters have changed and the reference central sections need to be updated
                call pftcc%apply_ctf(b%tfun, dfx, dfy, angast)
                if( ctfcheck ) print *, 'ctfcheck :: updating reference ctf'
            endif
            dfx_prev    = dfx
            dfy_prev    = dfy
            angast_prev = angast
        endif
    end subroutine

end module simple_hadamard_common
