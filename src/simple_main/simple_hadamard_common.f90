module simple_hadamard_common
use simple_defs
use simple_cmdline,      only: cmdline
use simple_build,        only: build
use simple_params,       only: params
use simple_ori,          only: ori
use simple_rnd,          only: ran3
use simple_prime3D_srch, only: prime3D_srch
use simple_gridding,     only: prep4cgrid
use simple_strings,      ! use all in there
use simple_math          ! use all in there
use simple_masker        ! use all in there
implicit none

public :: read_imgs_from_stk, set_bp_range, setup_shellweights, grid_ptcl, prepimg4align,&
eonorm_struct_facts, norm_struct_facts, preprefs4align, preprefvol, reset_prev_defparms, prep2Dref
private

interface prep2Dref
    module procedure prep2Dref_1
    module procedure prep2Dref_2
end interface

interface setup_shellweights
    module procedure setup_shellweights_1
    module procedure setup_shellweights_2
end interface

logical, parameter :: DEBUG        = .false.
real,    parameter :: SHTHRESH     = 0.0001
real,    parameter :: CENTHRESH    = 0.01   ! threshold for performing volume/cavg centering in pixels
real,    parameter :: MAXCENTHRESH = 0.025  ! max centering shift applied: 2.5% of box size
real               :: dfx_prev     = 0.
real               :: dfy_prev     = 0.
real               :: angast_prev  = 0.
real               :: kV_prev      = 0.
real               :: cs_prev      = 0.
real               :: fraca_prev   = 0.
    
contains

    subroutine read_imgs_from_stk( b, p )
        use simple_imgfile, only: imgfile
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        type(imgfile) :: ioimg
        integer       :: cnt, iptcl, istart, istop, i
        if( allocated(b%imgs) )then
            istart = lbound(b%imgs, dim=1)
            istop  = ubound(b%imgs, dim=1)
            do i=istart,istop
                call b%imgs(i)%kill
            end do
            deallocate(b%imgs)
        endif
        allocate(b%imgs(p%fromp:p%top))
        do iptcl=p%fromp,p%top
            call b%imgs(iptcl)%new([p%box,p%box,1],p%smpd,p%imgkind)
        end do
        if( p%l_distr_exec )then
            call b%imgs(p%fromp)%open(p%stk_part, ioimg)
        else
            call b%imgs(p%fromp)%open(p%stk, ioimg)
        endif
        cnt = 0
        do iptcl=p%fromp,p%top
            cnt = cnt + 1
            if( p%l_distr_exec )then
                call b%imgs(iptcl)%read(p%stk_part, cnt, ioimg=ioimg)
            else
                call b%imgs(iptcl)%read(p%stk, iptcl, ioimg=ioimg)
            endif
        end do
        call ioimg%close
    end subroutine read_imgs_from_stk
    
    subroutine set_bp_range( b, p, cline )
        use simple_estimate_ssnr, only: fsc2ssnr
        use simple_cmdline,       only: cmdline
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        real, allocatable     :: resarr(:)
        real                  :: fsc0143, fsc05, mapres(p%nstates)
        integer               :: s, loc(1)
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
                if( p%oritab.eq.'')all_fsc_bin_exist = ( count(fsc_bin_exists)==p%nstates )
                ! set low-pass Fourier index limit
                if( all_fsc_bin_exist )then
                    ! we need the worst resolved fsc
                    resarr = b%img%get_res()
                    do s=1,p%nstates
                        if( fsc_bin_exists(s) )then
                            ! these are the 'classical' resolution measures
                            fsc_fname   = adjustl('fsc_state'//int2str_pad(s,2)//'.bin')
                            b%fsc(s,:)  = file2rarr( trim(fsc_fname) )
                            b%ssnr(s,:) = fsc2ssnr(b%fsc(s,:))
                            call get_resolution(b%fsc(s,:), resarr, fsc05, fsc0143)
                            mapres(s)   = fsc0143
                        else
                            ! empty state
                            mapres(s)   = 0.
                            b%fsc(s,:)  = 0.
                            b%ssnr(s,:) = 0.
                        endif
                    end do
                    loc = maxloc(mapres)
                    p%kfromto(2) = get_lplim(b%fsc(loc(1),:))
                    if( p%kfromto(2) == 1 )then
                        stop 'simple_math::get_lplim gives nonsensical result (==1)'
                    endif
                    if( DEBUG ) write(*,*) '*** simple_hadamard_common ***: extracted FSC info'
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
                call b%a%set_all2single('lp',p%lp)
            case('no')
                ! set Fourier index range
                p%kfromto(1) = max(2,b%img%get_find(p%hp))
                if( cline%defined('lpstop') )then
                    p%kfromto(2) = min(b%img%get_find(p%lp),b%img%get_find(p%lpstop))
                else
                    p%kfromto(2) = b%img%get_find(p%lp)
                endif
                p%lp_dyn = p%lp
                call b%a%set_all2single('lp',p%lp)
            case DEFAULT
                stop 'Unsupported eo flag; simple_hadamard_common'
        end select
        if( DEBUG ) write(*,*) '*** simple_hadamard_common ***: did set Fourier index range'
    end subroutine set_bp_range

    !>  \brief  constructs the shellweight matrix for 3D search
    subroutine setup_shellweights_1( b, p, doshellweight, wmat, res, res_pad )
        use simple_map_reduce, only: merge_rmat_from_parts
        use simple_filterer,   only: normalise_shellweights
        class(build),                intent(inout) :: b
        class(params),               intent(inout) :: p
        logical,                     intent(out)   :: doshellweight
        real,           allocatable, intent(out)   :: wmat(:,:)
        real, optional, allocatable, intent(out)   :: res(:), res_pad(:)
        logical, allocatable :: files_exist(:)
        integer :: filtsz, filtsz_pad, alloc_stat, filnum, io_stat, ipart
        doshellweight = .false.
        if( .not. p%l_shellw ) return
        filtsz     = b%img%get_filtsz() ! nr of resolution elements
        filtsz_pad = b%img_pad%get_filtsz()
        if( allocated(wmat) ) deallocate(wmat)
        if( present(res) )then
            if( allocated(res) ) deallocate(res)
            res = b%img%get_res()
        endif
        if( present(res_pad) )then
            if( allocated(res_pad) ) deallocate(res_pad)
            res_pad = b%img_pad%get_res()
        endif
        if( p%l_distr_exec )then
            call wmat_from_single_file
            if( doshellweight )then
                ! we are done
                return
            else
                ! we may need to merge partial shellweight files
                allocate( files_exist(p%nparts) )
                do ipart=1,p%nparts
                    files_exist(ipart) = file_exists('shellweights_part'//int2str_pad(ipart,p%numlen)//'.bin')
                end do
                if( all(files_exist) )then
                    wmat = merge_rmat_from_parts(p%nptcls, p%nparts, filtsz, 'shellweights_part')
                    call normalise_shellweights(wmat)
                    doshellweight = .true.
                endif
                deallocate(files_exist)
            endif
        else
            call wmat_from_single_file
        endif

      contains

        subroutine wmat_from_single_file
            if( file_exists(p%shellwfile) )then    
                allocate( wmat(p%nptcls,filtsz), stat=alloc_stat)
                filnum = get_fileunit()
                open(unit=filnum, status='OLD', action='READ', file=p%shellwfile, access='STREAM')
                read(unit=filnum,pos=1,iostat=io_stat) wmat
                ! check if the read was successful
                if( io_stat .ne. 0 )then
                    doshellweight = .false.
                    return  
                endif
                close(filnum)
                call normalise_shellweights(wmat)
                doshellweight = .true.
            endif
        end subroutine wmat_from_single_file

    end subroutine setup_shellweights_1

    !>  \brief  constructs the shellweight matrix for 3D search
    subroutine setup_shellweights_2( b, p, doshellweight, wmat, npeaks, res, res_pad )
        use simple_map_reduce, only: merge_rmat_from_parts
        use simple_filterer,   only: normalise_shellweights
        class(build),                intent(inout) :: b
        class(params),               intent(inout) :: p
        logical,                     intent(out)   :: doshellweight
        real,           allocatable, intent(out)   :: wmat(:,:,:)
        integer,                     intent(in)    :: npeaks
        real, optional, allocatable, intent(out)   :: res(:), res_pad(:)
        logical, allocatable :: files_exist(:)
        integer :: filtsz, filtsz_pad, alloc_stat, filnum, io_stat, ipart
        doshellweight = .false.
        if( .not. p%l_shellw ) return
        filtsz     = b%img%get_filtsz() ! nr of resolution elements
        filtsz_pad = b%img_pad%get_filtsz()
        if( allocated(wmat) ) deallocate(wmat)
        if( present(res) )then
            if( allocated(res) ) deallocate(res)
            res = b%img%get_res()
        endif
        if( present(res_pad) )then
            if( allocated(res_pad) ) deallocate(res_pad)
            res_pad = b%img_pad%get_res()
        endif
        if( p%l_distr_exec )then
            call wmat_from_single_file
            if( doshellweight )then
                ! we are done
                return
            else
                ! we may need to merge partial shellweight files
                allocate( files_exist(p%nparts) )
                do ipart=1,p%nparts
                    files_exist(ipart) = file_exists('shellweights_part'//int2str_pad(ipart,p%numlen)//'.bin')
                end do
                if( all(files_exist) )then
                    wmat = merge_rmat_from_parts(p%nstates, p%nptcls, p%nparts, filtsz, 'shellweights_part')
                    call normalise_shellweights(wmat, npeaks)
                    doshellweight = .true.
                endif
                deallocate(files_exist)
            endif
        else
            call wmat_from_single_file
        endif

      contains

        subroutine wmat_from_single_file
            if( file_exists(p%shellwfile) )then    
                allocate( wmat(p%nstates,p%nptcls,filtsz), stat=alloc_stat)
                filnum = get_fileunit()
                open(unit=filnum, status='OLD', action='READ', file=p%shellwfile, access='STREAM')
                read(unit=filnum,pos=1,iostat=io_stat) wmat
                ! check if the read was successful
                if( io_stat .ne. 0 )then
                    write(*,'(a,i0,2a)') '**ERROR(setup_shellweights_2): I/O error ',&
                    io_stat, ' when reading'//trim(p%shellwfile)
                    stop 'I/O error; setup_shellweights_2; simple_hadamard_common'
                endif
                close(filnum)
                call normalise_shellweights(wmat, npeaks)
                doshellweight = .true.
            endif
        end subroutine wmat_from_single_file

    end subroutine setup_shellweights_2

    !>  \brief  grids one particle image to the volume
    subroutine grid_ptcl( b, p, iptcl, cnt_glob, orientation, os, shellweights )
        use simple_oris, only: oris
        class(build),              intent(inout) :: b
        class(params),             intent(inout) :: p
        integer,                   intent(in)    :: iptcl, cnt_glob
        class(ori),                intent(inout) :: orientation
        class(oris),     optional, intent(inout) :: os
        real,            optional, intent(in)    :: shellweights(:)
        real      :: pw, ran, w
        integer   :: jpeak, s, k
        type(ori) :: orisoft, o_sym
        logical   :: softrec
        softrec = .false.
        if( present(os) )then
            softrec = .true.
            if( p%npeaks /= os%get_noris() )&
                &stop 'non-congruent number of orientations; simple_hadamard_common :: grid_ptcl'
        endif
        if( p%npeaks>1 .and. .not.softrec )&
            stop 'need optional primesrch3D input when npeaks > 1; simple_hadamard_common :: grid_ptcl'
        pw = orientation%get('w')
        if( pw > 0. )then
            if( p%l_distr_exec )then
                call b%img_copy%read(p%stk_part, cnt_glob, isxfel=p%l_xfel)
            else
                call b%img_copy%read(p%stk, iptcl, isxfel=p%l_xfel)
            endif
            ! prepare image for gridding
            ! using the uncorrected/unmodified image as input
            if( p%l_xfel )then
                call b%img_copy%pad(b%img_pad)
            else
                call prep4cgrid(b%img_copy, b%img_pad, p%msk)
            endif
            if( DEBUG ) write(*,*) '*** simple_hadamard_common ***: prepared image for gridding'
            ran = ran3()
            orisoft = orientation
            do jpeak=1,p%npeaks
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
                    if( p%pgrp == 'c1' .or. str_has_substr(p%refine,'adasym') )then
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

    subroutine prepimg4align( b, p, o )
        use simple_ctf, only: ctf
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        type(ori),      intent(inout) :: o
        type(ctf) :: tfun
        real      :: x, y, dfx, dfy, angast
        integer   :: i!, icls
        if( p%l_xfel )then
            ! nothing to do 4 now
            return
        else
            x     = o%get('x')
            y     = o%get('y')
            !icls  = nint(o%get('class'))
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
            if( abs(x) > SHTHRESH .or. abs(y) > SHTHRESH ) call b%img%shift(-x, -y)
            ! back to real-space
            call b%img%bwd_ft
            ! clip image if needed
            if( p%boxmatch < p%box ) call b%img%clip_inplace([p%boxmatch,p%boxmatch,1]) ! SQUARE DIMS ASSUMED
            ! MASKING
            if( p%doautomsk )then
                ! PARTICLE ENVELOPPE MASKING TURNED OFF FOR NOW
                ! ! make sure that img_msk is of the correct dimension
                ! call b%img_msk%new([p%boxmatch,p%boxmatch,1],p%smpd)
                ! ! create 2D envelope
                ! call b%mskvol%env_rproject(o, b%img_msk, p%msk)
                ! do i=1, p%binwidth
                !     call b%img_msk%grow_bin
                ! enddo
                ! call b%img_msk%cos_edge(20)
                ! ! call b%img_msk%norm('sigm')
                ! ! multiply with the projected envelope
                ! call b%img%mul(b%img_msk)
            else if( p%automsk .eq. 'cavg' )then
                ! ab initio mask
                call automask2D(b%img, p)
            else              
                ! apply a soft-edged mask
                if( p%l_innermsk )then
                    call b%img%mask(p%msk, 'soft', inner=p%inner, width=p%width)
                else
                    call b%img%mask(p%msk, 'soft')
                endif
            endif
            ! return in Fourier space
            call b%img%fwd_ft
        endif
        if( DEBUG ) write(*,*) '*** simple_hadamard_common ***: finished prepimg4align'
    end subroutine prepimg4align

    subroutine prep2Dref_1( p, ref )
        use simple_image, only: image
        class(params),  intent(in)    :: p
        class(image),   intent(inout) :: ref
        ! normalise
        call ref%norm
        ! apply mask
        if( p%l_innermsk )then
            call ref%mask(p%msk, 'soft', inner=p%inner, width=p%width)
        else 
            call ref%mask(p%msk, 'soft')
        endif
        if( p%l_automsk ) call automask2D(ref, p)
        ! move to Fourier space
        call ref%fwd_ft
    end subroutine prep2Dref_1

    subroutine prep2Dref_2( p, ref, os, icls )
        use simple_image, only: image
        use simple_oris,  only: oris
        class(params),  intent(in)    :: p
        class(image),   intent(inout) :: ref
        class(oris),    intent(inout) :: os
        integer,        intent(in)    :: icls
        real :: xyz(3), sharg
        if( p%center.eq.'yes' .or. p%doshift )then

            ! TOOK OUT: PARTICLES JUMPING ALL OVER THE SHOP
            ! ! center the reference
            ! xyz   = ref%center(p%cenlp, 'no', p%msk, doshift=.false.)
            ! sharg = arg(xyz)
            ! if(sharg > CENTHRESH)then
            !     ! apply shift  and update the corresponding class parameters
            !     if(sharg > real(p%box)*MAXCENTHRESH) xyz = xyz / sharg
            !     call ref%shift(-xyz(1), -xyz(2))
            !     call os%add_shift2class(icls, xyz(1:2))
            ! endif

            ! center the reference and update the corresponding class parameters
            xyz = ref%center(p%cenlp, 'no', p%msk)
            call os%add_shift2class(icls, -xyz(1:2))
            
        endif
        ! normalise
        call ref%norm
        ! apply mask
        if( p%l_innermsk )then
            call ref%mask(p%msk, 'soft', inner=p%inner, width=p%width)
        else 
            call ref%mask(p%msk, 'soft')
        endif
        if( p%l_automsk ) call automask2D(ref, p)
        ! move to Fourier space
        call ref%fwd_ft
    end subroutine prep2Dref_2

    subroutine preprefvol( b, p, cline, s, doexpand )
        class(build),      intent(inout) :: b
        class(params),     intent(inout) :: p
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: s
        logical, optional, intent(in)    :: doexpand
        logical :: ddoexpand
        real    :: shvec(3)
        ddoexpand = .true.
        if( present(doexpand) ) ddoexpand = doexpand
        if( p%boxmatch < p%box )call b%vol%new([p%box,p%box,p%box],p%smpd) ! ensure correct dim
        call b%vol%read(p%vols(s), isxfel=p%l_xfel)
        if( p%l_xfel )then
            ! no centering
        else
            if(p%doshift) call centervol
        endif
        ! clip
        if( p%boxmatch < p%box )then
            call b%vol%clip_inplace([p%boxmatch,p%boxmatch,p%boxmatch]) ! SQUARE DIMS ASSUMED
        endif
        ! masking
        if( p%l_xfel )then
            ! no centering or masking
        else
            ! mask volume using a spherical soft-edged mask
            p%vols_msk(s) = add2fbody(p%vols(s), p%ext, 'msk')
            if( p%l_innermsk )then
                call b%vol%mask(p%msk, 'soft', inner=p%inner, width=p%width)
            else
                call b%vol%mask(p%msk, 'soft')
            endif
            ! mask using a molecular envelope
            if( p%doautomsk )then
                p%masks(s) = 'automask_state'//int2str_pad(s,2)//p%ext
                if( p%l_distr_exec )then
                    if( p%part == 1 )then
                        ! automask & write files
                        call automask(b, p, cline, b%vol, b%mskvol, p%vols_msk(s), p%masks(s))
                    else
                        ! automask & DO NOT write files
                        call automask(b, p, cline, b%vol, b%mskvol)
                    endif
                else
                    ! automask & write files
                    call automask(b, p, cline, b%vol, b%mskvol, p%vols_msk(s), p%masks(s))
                endif
            endif
        endif
        ! FT volume
        call b%vol%fwd_ft
        ! expand for fast interpolation
        if( ddoexpand ) call b%vol%expand_cmat

        contains

            subroutine centervol
                ! centering only for asymmetric and circular symmetric cases
                if( p%refine.ne.'het' )then
                    if( p%pgrp(:1).eq.'c' )then
                        shvec = b%vol%center(p%cenlp,'no',p%msk,doshift=.false.) ! find center of mass shift
                        if( arg(shvec) > CENTHRESH )then
                            if( p%pgrp.ne.'c1' ) shvec(1:2) = 0.         ! shifts only along z-axis for C2 and above
                            call b%vol%shift(shvec(1),shvec(2),shvec(3)) ! performs shift
                            ! map back to particle oritentations
                            if( cline%defined('oritab') ) call b%a%map3dshift22d(-shvec(:), state=s)
                        endif
                    endif
                endif
            end subroutine centervol

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
            if( p%l_distr_exec )then
                call b%eorecvols(s)%write_eos('recvol_state'//int2str_pad(s,2)//'_part'//int2str_pad(p%part,p%numlen))
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
                call b%eorecvols(s)%sum_eos
                call b%eorecvols(s)%sampl_dens_correct_eos(s)
                call b%eorecvols(s)%sampl_dens_correct_sum(b%vol)
                call b%vol%write(p%vols(s), del_if_exists=.true.)
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

    subroutine preprefs4align( b, p, iptcl, pftcc, ref )
        use simple_math,             only: euclid
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_ctf,              only: ctf
        class(build),            intent(inout) :: b
        class(params),           intent(inout) :: p
        integer,                 intent(in)    :: iptcl
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, optional,       intent(inout) :: ref
        type(ctf) :: tfun
        real      :: kV, cs, fraca, dfx, dfy, angast, dist
        if( p%ctf .ne. 'no' )then
            if( present(ref) )then
                if( ref<1 .or. ref>pftcc%get_nrefs() )&
                    &stop 'reference index out of bounds; simple_hadamard_common%preprefs4align'
            endif
            if( p%tfplan%mode .eq. 'astig' )then ! astigmatic CTF
                dfx = b%a%get(iptcl,'dfx')
                dfy = b%a%get(iptcl,'dfy')
                angast = b%a%get(iptcl,'angast')
            else if( p%tfplan%mode .eq. 'noastig' )then
                dfx = b%a%get(iptcl,'dfx')
                dfy = dfx
                angast = 0.
            else
                stop 'Unsupported ctf mode; preprefs4align; simple_hadamard_common'
            endif
            kV    = b%a%get(iptcl,'kv')
            cs    = b%a%get(iptcl,'cs')
            fraca = b%a%get(iptcl,'fraca')
            dist  = euclid([kV,cs,fraca,dfx,dfy,angast],[kV_prev,cs_prev,fraca_prev,dfx_prev,dfy_prev,angast_prev])
            if( dist < 0.001 )then
                ! CTF parameters are the same as for the previous particle & no update is needed
            else
                ! CTF parameters have changed and ctf object and the reference central sections need to be updated
                tfun = ctf(p%smpd, kV, cs, fraca)
                if( present(ref) )then
                    call pftcc%apply_ctf(iptcl, refvec=[ref,ref])
                else
                    call pftcc%apply_ctf(iptcl)
                endif
            endif
            kV_prev     = kV
            cs_prev     = cs
            fraca_prev  = fraca
            dfx_prev    = dfx
            dfy_prev    = dfy
            angast_prev = angast
        endif
    end subroutine preprefs4align

    subroutine reset_prev_defparms
        kV_prev     = 0.
        cs_prev     = 0.
        fraca_prev  = 0.
        dfx_prev    = 0.
        dfy_prev    = 0.
        angast_prev = 0.
    end subroutine reset_prev_defparms

end module simple_hadamard_common
