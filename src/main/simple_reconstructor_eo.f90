! 3D reconstruction of even-odd pairs for FSC estimation
module simple_reconstructor_eo
include 'simple_lib.f08'
use simple_reconstructor, only: reconstructor
use simple_masker,        only: masker
use simple_parameters,    only: params_glob
use simple_image,         only: image
use simple_sp_project,    only: sp_project
use simple_fsc
implicit none

public :: reconstructor_eo
private
#include "simple_local_flags.inc"

type :: reconstructor_eo
    private
    type(reconstructor) :: even
    type(reconstructor) :: odd
    type(reconstructor) :: eosum
    type(masker)        :: envmask
    character(len=4)    :: ext
    real                :: res_fsc05          !< target resolution at FSC=0.5
    real                :: res_fsc0143        !< target resolution at FSC=0.143
    real                :: smpd, msk, fny, pad_correction=1.
    integer             :: box=0, nstates=1, numlen=2, hpind_fsc=0
    logical             :: phaseplate = .false.
    logical             :: automsk    = .false.
    logical             :: exists     = .false.
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! SETTERS
    procedure          :: set_automsk
    procedure          :: reset_all
    procedure          :: reset_eos
    procedure, private :: reset_eoexp
    procedure, private :: reset_even
    procedure, private :: reset_odd
    procedure          :: reset_sum
    procedure          :: apply_weight
    procedure          :: set_sh_lim
    ! GETTERS
    procedure          :: get_kbwin
    procedure          :: get_res
    ! I/O
    ! writers
    procedure          :: write_eos
    procedure, private :: write_even
    procedure, private :: write_odd
    ! readers
    procedure          :: read_eos
    procedure, private :: read_even
    procedure, private :: read_odd
    ! INTERPOLATION
    procedure          :: grid_plane
    procedure          :: compress_exp
    procedure          :: expand_exp
    procedure          :: sum_eos    !< for merging even and odd into sum
    procedure          :: sum_reduce !< for summing eo_recs obtained by parallel exec
    procedure          :: sampl_dens_correct_eos
    procedure          :: sampl_dens_correct_sum
    ! RECONSTRUCTION
    procedure          :: eorec
    ! DESTRUCTORS
    procedure          :: kill_exp
    procedure          :: kill
end type reconstructor_eo

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine new( self,  spproj )
        class(reconstructor_eo), intent(inout) :: self   !< instance
        class(sp_project),       intent(inout) :: spproj !< project description
        logical     :: neg
        call self%kill
        ! set constants
        neg = .false.
        if( params_glob%neg .eq. 'yes' ) neg = .true.
        self%box        = params_glob%box
        self%smpd       = params_glob%smpd
        self%nstates    = params_glob%nstates
        self%fny        = params_glob%fny
        self%ext        = params_glob%ext
        self%numlen     = params_glob%numlen
        self%msk        = params_glob%msk
        self%automsk    = params_glob%l_filemsk .and. params_glob%l_envfsc
        self%phaseplate = params_glob%l_phaseplate
        self%hpind_fsc  = params_glob%hpind_fsc
        self%pad_correction = (real(params_glob%boxpd)/real(self%box))**3. * real(self%box)
        ! create composites
        if( self%automsk )then
            call self%envmask%new([params_glob%box,params_glob%box,params_glob%box], params_glob%smpd)
            call self%envmask%read(params_glob%mskfile)
        endif
        call self%even%new([params_glob%boxpd,params_glob%boxpd,params_glob%boxpd], params_glob%smpd)
        call self%even%alloc_rho( spproj)
        call self%even%set_ft(.true.)
        call self%odd%new([params_glob%boxpd,params_glob%boxpd,params_glob%boxpd], params_glob%smpd)
        call self%odd%alloc_rho(spproj)
        call self%odd%set_ft(.true.)
        call self%eosum%new([params_glob%boxpd,params_glob%boxpd,params_glob%boxpd], params_glob%smpd)
        call self%eosum%alloc_rho( spproj, expand=.false.)
        ! set existence
        self%exists = .true.
    end subroutine new

    ! SETTERS

    subroutine set_automsk( self, l_which )
        class(reconstructor_eo), intent(inout) :: self
        logical,                 intent(in)    :: l_which
        self%automsk = l_which
    end subroutine set_automsk

    !>  \brief  resets all
    subroutine reset_all( self )
        class(reconstructor_eo), intent(inout) :: self
        call self%reset_eos
        call self%reset_eoexp
        call self%reset_sum
    end subroutine reset_all

    !>  \brief  resets the even odd pairs
    subroutine reset_eos( self )
        class(reconstructor_eo), intent(inout) :: self
        call self%even%reset
        call self%odd%reset
    end subroutine reset_eos

    !>  \brief  resets the even odd pairs expanded matrices
    subroutine reset_eoexp( self )
        class(reconstructor_eo), intent(inout) :: self
        call self%even%reset_exp
        call self%odd%reset_exp
    end subroutine reset_eoexp

    !>  \brief  resets the even
    subroutine reset_even( self )
        class(reconstructor_eo), intent(inout) :: self
        call self%even%reset
    end subroutine reset_even

    !>  \brief  resets the odd
    subroutine reset_odd( self )
        class(reconstructor_eo), intent(inout) :: self
        call self%odd%reset
    end subroutine reset_odd

    !>  \brief  resets the sum
    subroutine reset_sum( self )
        class(reconstructor_eo), intent(inout) :: self
        call self%eosum%reset
    end subroutine reset_sum

    subroutine apply_weight( self, w )
        class(reconstructor_eo), intent(inout) :: self
        real,                    intent(in)    :: w
        call self%even%apply_weight(w)
        call self%odd%apply_weight(w)
    end subroutine apply_weight

    subroutine set_sh_lim(self, sh_lim)
        class(reconstructor_eo), intent(inout) :: self
        integer,                 intent(in)    :: sh_lim
        call self%even%set_sh_lim(sh_lim)
        call self%odd%set_sh_lim(sh_lim)
    end subroutine set_sh_lim

    ! GETTERS

    !>  \brief  return the window functions used by reconstructor_eo
    function get_kbwin( self ) result( wf )
        class(reconstructor_eo), intent(inout) :: self
        type(kbinterpol) :: wf
        wf = self%even%get_kbwin()
    end function get_kbwin

    !> \brief  for getting the resolution
    !> \param res_fsc05  target resolution a FSC=0.5
    !> \param res_fsc0143  target resolution a FSC=0.143
    subroutine get_res( self, res_fsc05, res_fsc0143 )
        class(reconstructor_eo), intent(in)  :: self !< instance
        real,                    intent(out) :: res_fsc05, res_fsc0143
        res_fsc0143 = self%res_fsc0143
        res_fsc05   = self%res_fsc05
    end subroutine get_res

    ! I/O

    !>  \brief  write the even and odd reconstructions
    subroutine write_eos( self, fbody )
        class(reconstructor_eo), intent(inout) :: self
        character(len=*),        intent(in)    :: fbody !< filename
        call self%write_even(fbody)
        call self%write_odd(fbody)
    end subroutine write_eos

    !>  \brief  write the even reconstruction
    subroutine write_even( self, fbody )
        class(reconstructor_eo), intent(inout) :: self
        character(len=*),        intent(in)    :: fbody
        call self%even%write(trim(adjustl(fbody))//'_even'//self%ext, del_if_exists=.true.)
        call self%even%write_rho(trim('rho_'//trim(adjustl(fbody))//'_even'//self%ext))
    end subroutine write_even

    !>  \brief  write the odd reconstruction
    subroutine write_odd( self, fbody )
        class(reconstructor_eo), intent(inout) :: self
        character(len=*),        intent(in)    :: fbody
        call self%odd%write(trim(adjustl(fbody))//'_odd'//self%ext, del_if_exists=.true.)
        call self%odd%write_rho('rho_'//trim(adjustl(fbody))//'_odd'//self%ext)
    end subroutine write_odd

    !>  \brief read the even and odd reconstructions
    subroutine read_eos( self, fbody )
        class(reconstructor_eo), intent(inout) :: self
        character(len=*),        intent(in)    :: fbody
        call self%read_even(fbody)
        call self%read_odd(fbody)
    end subroutine read_eos

    !>  \brief  read the even reconstruction
    subroutine read_even( self, fbody )
        class(reconstructor_eo), intent(inout) :: self
        character(len=*),        intent(in)    :: fbody
        character(len=STDLEN)                  :: even_vol, even_rho
        logical                                :: here(2)
        even_vol = trim(adjustl(fbody))//'_even'//self%ext
        even_rho = 'rho_'//trim(adjustl(fbody))//'_even'//self%ext
        here(1)= file_exists(even_vol)
        here(2)= file_exists(even_rho)
        if( all(here) )then
            call self%even%read(even_vol)
            call self%even%read_rho(even_rho)
        else
            call self%reset_even
        endif
    end subroutine read_even

    !>  \brief  read the odd reconstruction
    subroutine read_odd( self, fbody )
        class(reconstructor_eo), intent(inout) :: self
        character(len=*),        intent(in)    :: fbody
        character(len=STDLEN)                  :: odd_vol, odd_rho
        logical                                :: here(2)
        odd_vol = trim(adjustl(fbody))//'_odd'//self%ext
        odd_rho = 'rho_'//trim(adjustl(fbody))//'_odd'//self%ext
        here(1)= file_exists(odd_vol)
        here(2)= file_exists(odd_rho)
        if( all(here) )then
            call self%odd%read(odd_vol)
            call self%odd%read_rho(odd_rho)
        else
            call self%reset_odd
        endif
    end subroutine read_odd

    ! INTERPOLATION

    !> \brief  for gridding a Fourier plane
    subroutine grid_plane( self, se, o, fpl, eo, pwght )
        use simple_fplane, only: fplane
        class(reconstructor_eo), intent(inout) :: self    !< instance
        class(sym),              intent(inout) :: se      !< symmetry elements
        class(ori),              intent(inout) :: o       !< orientation
        class(fplane),           intent(in)    :: fpl     !< Forurier & ctf planes
        integer,                 intent(in)    :: eo      !< eo flag
        real,                    intent(in)    :: pwght   !< external particle weight (affects both fplane and rho)
        select case(eo)
            case(-1,0)
                call self%even%insert_plane(se, o, fpl, pwght)
            case(1)
                call self%odd%insert_plane(se, o, fpl, pwght)
            case DEFAULT
                THROW_HARD('unsupported eo flag; grid_plane')
        end select
    end subroutine grid_plane

    !> \brief  for summing the even odd pairs, resulting sum in self%even
    subroutine sum_eos( self )
        class(reconstructor_eo), intent(inout) :: self !< instance
        call self%eosum%reset
        call self%eosum%sum_reduce(self%even)
        call self%eosum%sum_reduce(self%odd)
    end subroutine sum_eos

    !> \brief  for summing reconstructors generated by parallel execution
    subroutine sum_reduce( self, self_in )
         class(reconstructor_eo), intent(inout) :: self
         class(reconstructor_eo), intent(in)    :: self_in
         call self%even%sum_reduce(self_in%even)
         call self%odd%sum_reduce(self_in%odd)
    end subroutine sum_reduce

    !>  \brief compress e/o
    subroutine compress_exp( self )
        class(reconstructor_eo), intent(inout) :: self
        call self%even%compress_exp
        call self%odd%compress_exp
    end subroutine compress_exp

    !>  \brief expand e/o
    subroutine expand_exp( self )
        class(reconstructor_eo), intent(inout) :: self
        call self%even%expand_exp
        call self%odd%expand_exp
    end subroutine expand_exp

    !> \brief  for sampling density correction of the eo pairs
    subroutine sampl_dens_correct_eos( self, state, fname_even, fname_odd, find4eoavg )
        class(reconstructor_eo), intent(inout) :: self                  !< instance
        integer,                 intent(in)    :: state                 !< state
        character(len=*),        intent(in)    :: fname_even, fname_odd !< even/odd filenames
        integer,                 intent(out)   :: find4eoavg            !< Fourier index for eo averaging
        type(image)           :: even, odd
        complex,  allocatable :: cmat(:,:,:)
        real,     allocatable :: res(:), fsc(:)
        real                  :: lp_rand, msk
        integer               :: k,k_rand, find_plate, filtsz
        logical               :: l_combined, l_euclid_regularization
        msk    = real(self%box / 2) - COSMSKHALFWIDTH - 1.
        ! msk  = self%msk ! for a tighter spherical mask
        filtsz = fdim(self%box) - 1
        res    = get_resarr(self%box, self%smpd)
        allocate(fsc(filtsz),source=0.)
        ! if e=o then SSNR will be adjusted
        l_combined = trim(params_glob%combine_eo).eq.'yes'
        ! ML-regularization
        l_euclid_regularization = params_glob%l_ml_reg
        if( l_euclid_regularization )then
            ! preprocessing for FSC calculation
            ! even
            cmat = self%even%get_cmat()
            call self%even%sampl_dens_correct(do_gridcorr=.false.)
            even = self%even
            call self%even%set_cmat(cmat)
            deallocate(cmat)
            call even%ifft()
            call even%clip_inplace([self%box,self%box,self%box])
            call even%div(self%pad_correction)
            if( l_combined ) call even%write(add2fbody(fname_even,params_glob%ext,'_unfil'))
            ! odd
            cmat = self%odd%get_cmat()
            call self%odd%sampl_dens_correct(do_gridcorr=.false.)
            odd = self%odd
            call self%odd%set_cmat(cmat)
            deallocate(cmat)
            call odd%ifft()
            call odd%clip_inplace([self%box,self%box,self%box])
            call odd%div(self%pad_correction)
            if( l_combined ) call odd%write(add2fbody(fname_odd,params_glob%ext,'_unfil'))
            ! masking
            if( self%automsk )then
                call even%mul(self%envmask)
                call odd%mul(self%envmask)
            else
                call even%mask(msk, 'soft', backgr=0.)
                call odd%mask(msk, 'soft', backgr=0.)
            endif
            ! calculate FSC
            call even%fft()
            call odd%fft()
            call even%fsc(odd, fsc)
            ! regularization
            call self%even%add_invtausq2rho(fsc)
            call self%odd%add_invtausq2rho(fsc)
            ! Even: uneven sampling density correction, clip, & write
            cmat = self%even%get_cmat()
            call self%even%sampl_dens_correct(do_gridcorr=.false.)
            call self%even%ifft
            call even%zero_and_unflag_ft
            call self%even%clip(even)
            call even%div(self%pad_correction)
            call even%write(trim(fname_even), del_if_exists=.true.)
            call self%even%set_cmat(cmat)
            call even%kill
            deallocate(cmat)
            ! Odd: uneven sampling density correction, clip, & write
            cmat = self%odd%get_cmat()
            call self%odd%sampl_dens_correct(do_gridcorr=.false.)
            call self%odd%ifft
            call odd%zero_and_unflag_ft
            call self%odd%clip(odd)
            call odd%div(self%pad_correction)
            call odd%write(trim(fname_odd), del_if_exists=.true.)
            call self%odd%set_cmat(cmat)
            call odd%kill
            deallocate(cmat)
        else
            ! make clipped volumes
            call even%new([self%box,self%box,self%box],self%smpd)
            call odd%new([self%box,self%box,self%box],self%smpd)
            ! correct for the uneven sampling density
            call self%even%sampl_dens_correct(do_gridcorr=L_DO_GRIDCORR_GLOB)
            call self%odd%sampl_dens_correct(do_gridcorr=L_DO_GRIDCORR_GLOB)
            ! reverse FT
            call self%even%ifft()
            call self%odd%ifft()
            ! clip
            call self%even%clip(even)
            call self%odd%clip(odd)
            ! FFTW padding correction
            call even%div(self%pad_correction)
            call odd%div(self%pad_correction)
            ! write un-normalised unmasked even/odd volumes
            call even%write(trim(fname_even), del_if_exists=.true.)
            call odd%write(trim(fname_odd),   del_if_exists=.true.)
            ! masking
            if( self%automsk )then
                call even%mul(self%envmask)
                call odd%mul(self%envmask)
            else
                call even%mask(msk, 'soft', backgr=0.)
                call odd%mask(msk, 'soft', backgr=0.)
            endif
            ! calculate FSC
            call even%fft()
            call odd%fft()
            call even%fsc(odd, fsc)
        endif
        find_plate = 0
        if( self%phaseplate ) call phaseplate_correct_fsc(fsc, find_plate)
        if( self%hpind_fsc > 0 ) fsc(:self%hpind_fsc) = fsc(self%hpind_fsc + 1)
        ! save, get & print resolution
        call arr2file(fsc, 'fsc_state'//int2str_pad(state,2)//BIN_EXT)
        call get_resolution(fsc, res, self%res_fsc05, self%res_fsc0143)
        self%res_fsc05   = max(self%res_fsc05,self%fny)
        self%res_fsc0143 = max(self%res_fsc0143,self%fny)
        if( trim(params_glob%silence_fsc) .eq. 'no' )then
            do k=1,size(res)
               write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(k), '>>> CORRELATION:', fsc(k)
            end do
            write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', self%res_fsc05
            write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', self%res_fsc0143
        endif
        ! Fourier index for eo averaging
        if( self%hpind_fsc > 0 )then
            find4eoavg = self%hpind_fsc
        else
            find4eoavg = max(K4EOAVGLB,  calc_fourier_index(FREQ4EOAVG3D,self%box,self%smpd))
            find4eoavg = min(find4eoavg, get_lplim_at_corr(fsc, FSC4EOAVG3D))
            find4eoavg = max(find4eoavg, find_plate)
        endif
        deallocate(fsc, res)
        call even%kill
        call odd%kill
    end subroutine sampl_dens_correct_eos

    !> \brief  for sampling density correction, antialiasing, ifft & normalization of the sum
    subroutine sampl_dens_correct_sum( self, reference )
        class(reconstructor_eo), intent(inout) :: self      !< instance
        class(image),            intent(inout) :: reference !< reference volume
        if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> SAMPLING DENSITY (RHO) CORRECTION & WIENER NORMALIZATION'
        call reference%set_ft(.false.)
        call self%eosum%sampl_dens_correct(do_gridcorr=L_DO_GRIDCORR_GLOB)
        call self%eosum%ifft()
        call self%eosum%div(self%pad_correction)
        call self%eosum%clip(reference)
    end subroutine sampl_dens_correct_sum

    ! RECONSTRUCTION

    !> \brief  for distributed reconstruction of even/odd maps
    subroutine eorec( self, spproj, o, se, state, fbody )
        use simple_fplane,   only: fplane
        class(reconstructor_eo),    intent(inout) :: self   !< object
        class(sp_project),          intent(inout) :: spproj !< project description
        class(oris),                intent(inout) :: o      !< orientations
        class(sym),                 intent(inout) :: se     !< symmetry element
        integer,                    intent(in)    :: state  !< state to reconstruct
        character(len=*), optional, intent(in)    :: fbody  !< body of output file
        character(len=:), allocatable :: recname, volname
        character(len=LONGSTRLEN)     :: eonames(2)
        type(fplane)         :: fpl
        type(image)          :: img, mskimg, vol
        logical, allocatable :: lmsk(:,:,:)
        type(ctfparams)      :: ctfvars
        real                 :: sdev_noise
        integer              :: statecnt(params_glob%nstates), i, cnt, state_here, state_glob, find4eoavg
        ! stash global state index
        state_glob = state
        ! make the images
        call img%new([params_glob%box,params_glob%box,1],params_glob%smpd)
        call mskimg%disc([params_glob%box,params_glob%box,1], params_glob%smpd, params_glob%msk, lmsk)
        call vol%new([params_glob%box,params_glob%box,params_glob%box], params_glob%smpd)
        call fpl%new(img)
        ! zero the Fourier volumes and rhos
        call self%reset_all
        call self%reset_eoexp
        if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> KAISER-BESSEL INTERPOLATION'
        statecnt = 0
        cnt      = 0
        do i=1,params_glob%nptcls
            call progress(i, params_glob%nptcls)
            if( i <= params_glob%top .and. i >= params_glob%fromp )then
                cnt = cnt + 1
                state_here = nint(o%get(i,'state'))
                if( state_here > 0 .and. (state_here == state ) )then
                    statecnt(state) = statecnt(state) + 1
                    call rec_dens
                endif
            endif
        end do
        ! undo fourier components expansion
        call self%compress_exp
        ! density correction & output
        if( params_glob%l_distr_exec )then
            if( present(fbody) )then
                call self%write_eos(fbody//int2str_pad(state,2)//'_part'//int2str_pad(params_glob%part,self%numlen))
            else
                call self%write_eos('recvol_state'//int2str_pad(state,2)//'_part'//int2str_pad(params_glob%part,self%numlen))
            endif
        else
            allocate(recname, source=trim(VOL_FBODY)//int2str_pad(state,2))
            allocate(volname, source=recname//params_glob%ext)
            eonames(1) = trim(recname)//'_even'//params_glob%ext
            eonames(2) = trim(recname)//'_odd'//params_glob%ext
            if( params_glob%l_ml_reg )then
                call self%sampl_dens_correct_eos(state, eonames(1), eonames(2), find4eoavg)
                call self%sum_eos
            else
                call self%sum_eos
                call self%sampl_dens_correct_eos(state, eonames(1), eonames(2), find4eoavg)
            endif
            call self%sampl_dens_correct_sum(vol)
            call vol%write(volname, del_if_exists=.true.)
        endif
        call img%kill
        call mskimg%kill
        call vol%kill
        call fpl%kill
        ! report how many particles were used to reconstruct each state
        if( params_glob%nstates > 1 )then
            write(logfhandle,'(a,1x,i3,1x,a,1x,i6)') '>>> NR OF PARTICLES INCLUDED IN STATE:', state, 'WAS:', statecnt(state)
        endif

        contains

            !> \brief  the density reconstruction functionality
            subroutine rec_dens
                character(len=:), allocatable :: stkname
                type(ori) :: orientation
                integer   :: state, ind_in_stk, eo
                real      :: pw
                state = nint(o%get(i, 'state'))
                if( state == 0 ) return
                call o%get_ori(i, orientation)
                ! eo-flag
                eo = nint(orientation%get('eo'))
                ! particle-weight
                pw = 1.
                if( orientation%isthere('w') ) pw = orientation%get('w')
                if( pw > TINY )then
                    call spproj%get_stkname_and_ind(params_glob%oritype, i, stkname, ind_in_stk)
                    call img%read(stkname, ind_in_stk)
                    ctfvars = spproj%get_ctfparams(params_glob%oritype, i)
                    call img%noise_norm(lmsk, sdev_noise)
                    call img%fft
                    call fpl%gen_planes(img, ctfvars)
                    call self%grid_plane(se, orientation, fpl, eo, pw)
                    deallocate(stkname)
                endif
                call orientation%kill
            end subroutine rec_dens

    end subroutine eorec

    ! DESTRUCTORS

    !>  \brief  is the expanded destructor
    subroutine kill_exp( self )
        class(reconstructor_eo), intent(inout) :: self !< instance
        if( self%exists )then
            call self%even%dealloc_exp
            call self%odd%dealloc_exp
            call self%eosum%dealloc_exp
        endif
    end subroutine kill_exp

    !>  \brief  is a destructor
    subroutine kill( self )
        class(reconstructor_eo), intent(inout) :: self !< instance
        if( self%exists )then
            ! kill composites
            call self%envmask%kill
            call self%even%dealloc_rho
            call self%even%kill
            call self%odd%dealloc_rho
            call self%odd%kill
            call self%eosum%dealloc_rho
            call self%eosum%kill
            ! set existence
            self%exists = .false.
        endif
    end subroutine kill

end module simple_reconstructor_eo
