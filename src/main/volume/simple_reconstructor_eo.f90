! 3D reconstruction of even-odd pairs for FSC estimation
module simple_reconstructor_eo
include 'simple_lib.f08'
use simple_reconstructor, only: reconstructor
use simple_image_msk,     only: image_msk
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
    type(image_msk)     :: envmask
    type(string)        :: ext
    real, allocatable   :: fsc(:)
    real                :: res_fsc05          !< target resolution at FSC=0.5
    real                :: res_fsc0143        !< target resolution at FSC=0.143
    real                :: smpd, msk, fny
    real                :: mag_correction=1.  !< scaling factor to correct for slice insertion, cropping & padding
    integer             :: ldim(3), box=0, boxpd=0
    integer             :: nstates=1, numlen=2, filtsz=0
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
    procedure          :: get_rhoexp_ptr
    ! I/O
    ! writers
    procedure          :: write_eos
    procedure, private :: write_even
    procedure, private :: write_odd
    ! readers
    procedure          :: read_eos
    procedure, private :: read_eos_parallel_io
    procedure, private :: read_even
    procedure, private :: read_odd
    ! INTERPOLATION
    procedure          :: grid_plane, test_grid_plane
    procedure          :: compress_exp
    procedure          :: expand_exp
    procedure          :: sum_eos    !< for merging even and odd into sum
    procedure          :: sum_reduce !< for summing eo_recs obtained by parallel exec
    procedure          :: sampl_dens_correct_eos
    procedure          :: calc_fsc4sampl_dens_correct
    procedure          :: write_fsc2txt
    procedure          :: sampl_dens_correct_sum
    ! DESTRUCTORS
    procedure          :: kill_exp
    procedure          :: kill
end type reconstructor_eo

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine new( self,  spproj, expand )
        class(reconstructor_eo), intent(inout) :: self   !< instance
        class(sp_project),       intent(inout) :: spproj !< project description
        logical,       optional, intent(in)    :: expand
        call self%kill
        ! set constants
        self%box     = params_glob%box_crop
        self%boxpd   = params_glob%box_croppd
        self%smpd    = params_glob%smpd_crop
        self%fny     = 2.*self%smpd
        self%nstates = params_glob%nstates
        self%ext     = params_glob%ext
        self%numlen  = params_glob%numlen
        self%filtsz  = fdim(self%box) - 1
        self%msk     = real(self%box / 2) - COSMSKHALFWIDTH - 1.
        self%automsk = params_glob%l_filemsk .and. params_glob%l_envfsc
        ! overall magnitude correction (interpolation padding + downscaling)
        self%mag_correction = real(params_glob%box)                                        ! insertion of 2D slice into 3D
        self%mag_correction = self%mag_correction * (real(self%boxpd) / real(self%box))**3 ! 3D padding factor
        self%ldim           = [self%boxpd,self%boxpd,self%boxpd]
        ! create composites
        if( self%automsk )then
            call self%envmask%new([self%box,self%box,self%box],self%smpd)
            call self%envmask%read(params_glob%mskfile)
        endif
        call self%even%new(self%ldim, params_glob%smpd)
        call self%even%alloc_rho(spproj, expand=expand)
        call self%even%set_ft(.true.)
        call self%odd%new(self%ldim, params_glob%smpd)
        call self%odd%alloc_rho(spproj, expand=expand)
        call self%odd%set_ft(.true.)
        call self%eosum%new(self%ldim, params_glob%smpd)
        call self%eosum%alloc_rho(spproj, expand=.false.)
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

    subroutine get_rhoexp_ptr( self, str_eo, rho_ptr )
        class(reconstructor_eo), target, intent(in)  :: self
        character(len=*),                intent(in)  :: str_eo
        real(kind=c_float), pointer,     intent(out) :: rho_ptr(:,:,:)
        select case(trim(str_eo))
            case('even')
                call self%even%get_rhoexp_ptr(rho_ptr)
            case('odd')
                call self%odd%get_rhoexp_ptr(rho_ptr)
        end select
    end subroutine get_rhoexp_ptr

    ! I/O

    !>  \brief  write the even and odd reconstructions
    subroutine write_eos( self, fbody )
        class(reconstructor_eo), intent(inout) :: self
        class(string),           intent(in)    :: fbody !< filename
        call self%write_even(fbody)
        call self%write_odd(fbody)
    end subroutine write_eos

    !>  \brief  write the even reconstruction
    subroutine write_even( self, fbody )
        class(reconstructor_eo), intent(inout) :: self
        class(string),           intent(in)    :: fbody
        call self%even%write(fbody//'_even'//self%ext%to_char(), del_if_exists=.true.)
        call self%even%write_rho(string('rho_')//fbody//'_even'//self%ext%to_char())
    end subroutine write_even

    !>  \brief  write the odd reconstruction
    subroutine write_odd( self, fbody )
        class(reconstructor_eo), intent(inout) :: self
        class(string),           intent(in)    :: fbody
        call self%odd%write(fbody//'_odd'//self%ext%to_char(), del_if_exists=.true.)
        call self%odd%write_rho(string('rho_')//fbody//'_odd'//self%ext%to_char())
    end subroutine write_odd

    !>  \brief read the even and odd reconstructions
    subroutine read_eos( self, fbody )
        class(reconstructor_eo), intent(inout) :: self
        class(string),           intent(in)    :: fbody
        logical, parameter :: SERIAL_READ = .false.
        if( SERIAL_READ )then
            call self%read_even(fbody)
            call self%read_odd(fbody)
        else
            call self%read_eos_parallel_io(fbody)
        endif
    end subroutine read_eos

    !>  \brief read the even and odd reconstructions
    subroutine read_eos_parallel_io( self, fbody )
        use simple_imgfile, only: imgfile
        class(reconstructor_eo), intent(inout) :: self
        class(string),           intent(in)    :: fbody
        real(kind=c_float),            pointer :: rmat_ptr(:,:,:) => null() !< image pixels/voxels (in data)
        real(kind=c_float),            pointer :: rho_ptr(:,:,:)  => null() !< sampling+CTF**2 density
        complex(kind=c_float_complex), pointer :: pcmate(:,:,:),pcmato(:,:,:), pprevcmate(:,:,:),pprevcmato(:,:,:)
        real(kind=c_float),            pointer :: prhoe(:,:,:), prhoo(:,:,:)
        real,                      allocatable :: rho_e(:,:,:), rho_o(:,:,:)
        type(string)  :: even_vol, even_rho, odd_vol, odd_rho
        type(image)   :: prev_vol_e, prev_vol_o
        type(imgfile) :: ioimg_e, ioimg_o
        integer       :: lims(3,2), cshape(3), prev_ldim(3), phys_out(3),phys_in(3)
        integer       :: h,k,l, fhandle_rho_e, fhandle_rho_o, i, ierr, dummy
        real          :: prev_smpd
        logical       :: here(4), l_pad_with_zeros
        even_vol = fbody//'_even'//self%ext%to_char()
        even_rho = string('rho_')//fbody//'_even'//self%ext%to_char()
        odd_vol  = fbody//'_odd'//self%ext
        odd_rho  = string('rho_')//fbody//'_odd'//self%ext%to_char()
        here(1)  = file_exists(even_vol)
        here(2)  = file_exists(even_rho)
        here(3)  = file_exists(odd_vol)
        here(4)  = file_exists(odd_rho)
        if( all(here) )then
            l_pad_with_zeros = .false.
            if( params_glob%l_update_frac )then
                ! check dimensions
                call find_ldim_nptcls(even_vol, prev_ldim, dummy, smpd=prev_smpd)
                if( prev_ldim(1) == self%ldim(1) )then
                    ! all good
                elseif( prev_ldim(1) > self%ldim(1) )then
                    THROW_HARD('Incorrect dimensions')
                else
                    l_pad_with_zeros = .true.
                endif
            endif
            call fopen(fhandle_rho_e, file=even_rho, status='OLD', action='READ', access='STREAM', iostat=ierr)
            call fileiochk('simple_reconstructor_eo::read_eos_parallel_io, opening '//even_rho%to_char(), ierr)
            call fopen(fhandle_rho_o, file=odd_rho,  status='OLD', action='READ', access='STREAM', iostat=ierr)
            call fileiochk('simple_reconstructor_eo::read_eos_parallel_io, opening '//odd_rho%to_char(),  ierr)
            if( l_pad_with_zeros )then
                call ioimg_e%open(even_vol, prev_ldim, prev_smpd, formatchar='M', readhead=.false., rwaction='read')
                call ioimg_o%open(odd_vol,  prev_ldim, prev_smpd, formatchar='M', readhead=.false., rwaction='read')
                call prev_vol_e%new(prev_ldim, prev_smpd)
                call prev_vol_o%new(prev_ldim, prev_smpd)
                cshape = [fdim(prev_ldim(1)), prev_ldim(2), prev_ldim(3)]
                allocate(rho_e(1:cshape(1),1:cshape(2),1:cshape(3)), rho_o(1:cshape(1),1:cshape(2),1:cshape(3)))
                call self%reset_even
                call self%reset_odd
                ! read
                !$omp parallel do default(shared) private(i,rmat_ptr,rho_ptr,ierr) schedule(static) num_threads(4)
                do i = 1, 4
                    select case(i)
                        case(1)
                            call prev_vol_e%get_rmat_ptr(rmat_ptr)
                            call ioimg_e%rSlices(1,prev_ldim(1),rmat_ptr,is_mrc=.true.)
                        case(2)
                            call prev_vol_o%get_rmat_ptr(rmat_ptr)
                            call ioimg_o%rSlices(1,prev_ldim(1),rmat_ptr,is_mrc=.true.)
                        case(3)
                            read(fhandle_rho_e, pos=1, iostat=ierr) rho_e
                            if( ierr .ne. 0 )&
                                &call fileiochk('simple_reconstructor_eo::read_eos_parallel_io, reading '// even_rho%to_char(), ierr)
                        case(4)
                            read(fhandle_rho_o, pos=1, iostat=ierr) rho_o
                            if( ierr .ne. 0 )&
                                &call fileiochk('simple_reconstructor_eo::read_eos_parallel_io, reading '// odd_rho%to_char(), ierr)
                        end select
                end do
                !$omp end parallel do
                ! pad
                lims = prev_vol_e%loop_lims(2)
                call self%even%get_cmat_ptr(pcmate)
                call self%odd%get_cmat_ptr(pcmato)
                call self%even%get_rho_ptr(prhoe)
                call self%odd%get_rho_ptr(prhoo)
                call prev_vol_e%get_cmat_ptr(pprevcmate)
                call prev_vol_o%get_cmat_ptr(pprevcmato)
                !$omp parallel do collapse(3) schedule(static) default(shared)&
                !$omp private(h,k,l,phys_out,phys_in) proc_bind(close)
                do h=lims(1,1),lims(1,2)
                    do k=lims(2,1),lims(2,2)
                        do l=lims(3,1),lims(3,2)
                            phys_out = self%even%comp_addr_phys(h,k,l)
                            phys_in  = prev_vol_e%comp_addr_phys(h,k,l)
                            pcmate(phys_out(1),phys_out(2),phys_out(3))= pprevcmate(phys_in(1),phys_in(2),phys_in(3))
                            pcmato(phys_out(1),phys_out(2),phys_out(3))= pprevcmato(phys_in(1),phys_in(2),phys_in(3))
                            prhoe(phys_out(1),phys_out(2),phys_out(3)) = rho_e(phys_in(1),phys_in(2),phys_in(3))
                            prhoo(phys_out(1),phys_out(2),phys_out(3)) = rho_o(phys_in(1),phys_in(2),phys_in(3))
                        end do
                    end do
                end do
                !$omp end parallel do
                call prev_vol_e%kill
                call prev_vol_o%kill
                deallocate(rho_e,rho_o)
            else
                call ioimg_e%open(even_vol, self%ldim, self%even%get_smpd(), formatchar='M', readhead=.false., rwaction='read')
                call ioimg_o%open(odd_vol,  self%ldim, self%odd%get_smpd(),  formatchar='M', readhead=.false., rwaction='read')
                !$omp parallel do default(shared) private(i,rmat_ptr,rho_ptr,ierr) schedule(static) num_threads(4)
                do i = 1, 4
                    select case(i)
                        case(1)
                            call self%even%get_rmat_ptr(rmat_ptr)
                            call ioimg_e%rSlices(1,self%ldim(1),rmat_ptr,is_mrc=.true.)
                        case(2)
                            call self%odd%get_rmat_ptr(rmat_ptr)
                            call ioimg_o%rSlices(1,self%ldim(1),rmat_ptr,is_mrc=.true.)
                        case(3)
                            call self%even%get_rho_ptr(rho_ptr)
                            read(fhandle_rho_e, pos=1, iostat=ierr) rho_ptr
                            if( ierr .ne. 0 )&
                                &call fileiochk('simple_reconstructor_eo::read_eos_parallel_io, reading '// even_rho%to_char(), ierr)
                        case(4)
                            call self%odd%get_rho_ptr(rho_ptr)
                            read(fhandle_rho_o, pos=1, iostat=ierr) rho_ptr
                            if( ierr .ne. 0 )&
                                &call fileiochk('simple_reconstructor_eo::read_eos_parallel_io, reading '// odd_rho%to_char(), ierr)
                    end select
                end do
                !$omp end parallel do
            endif
            call ioimg_e%close
            call ioimg_o%close
            call fclose(fhandle_rho_e)
            call fclose(fhandle_rho_o)
        else
            call self%reset_even
            call self%reset_odd
        endif
    end subroutine read_eos_parallel_io

    !>  \brief  read the even reconstruction
    subroutine read_even( self, fbody )
        class(reconstructor_eo), intent(inout) :: self
        class(string),           intent(in)    :: fbody
        type(string) :: even_vol, even_rho
        logical      :: here(2)
        even_vol = fbody//'_even'//self%ext%to_char()
        even_rho = string('rho_')//fbody//'_even'//self%ext%to_char()
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
        class(string),           intent(in)    :: fbody
        type(string) :: odd_vol, odd_rho
        logical      :: here(2)
        odd_vol = fbody//'_odd'//self%ext%to_char()
        odd_rho = string('rho_')//fbody//'_odd'//self%ext%to_char()
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

    !> \brief  for testing gridding a Fourier plane
    subroutine test_grid_plane( self, se, o, fpl, eo, pwght, stride )
        use simple_fplane, only: fplane
        class(reconstructor_eo), intent(inout) :: self
        class(sym),              intent(inout) :: se
        class(ori),              intent(inout) :: o
        class(fplane),           intent(in)    :: fpl
        integer,                 intent(in)    :: eo, stride
        real,                    intent(in)    :: pwght
        select case(eo)
            case(-1,0)
                call self%even%test_insert_plane(se, o, fpl, pwght, stride)
            case(1)
                call self%odd%test_insert_plane(se, o, fpl, pwght, stride)
            case DEFAULT
                THROW_HARD('unsupported eo flag; test_grid_plane')
        end select
    end subroutine test_grid_plane

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
    subroutine sampl_dens_correct_eos( self, state, fname_even, fname_odd, find4eoavg, fsc_in )
        class(reconstructor_eo), intent(inout) :: self                   !< instance
        integer,                 intent(in)    :: state                  !< state
        class(string),           intent(in)    :: fname_even, fname_odd  !< even/odd filenames
        integer,                 intent(out)   :: find4eoavg             !< Fourier index for eo averaging
        real, optional,          intent(in)    :: fsc_in(self%filtsz)    !< inputted fsc
        type(image)           :: even, odd
        complex,  allocatable :: cmat(:,:,:)
        real,     allocatable :: res(:)
        logical               :: l_have_fsc
        res = get_resarr(self%box, self%smpd)
        if( allocated(self%fsc) ) deallocate(self%fsc)
        if( present(fsc_in) )then
            allocate(self%fsc(self%filtsz),source=fsc_in)
            l_have_fsc = .true.
        else
            allocate(self%fsc(self%filtsz),source=0.)
            l_have_fsc = .false.
        endif
        ! ML-regularization
        if( params_glob%l_ml_reg )then
            ! preprocessing for FSC calculation
            ! even
            cmat = self%even%get_cmat()
            call self%even%sampl_dens_correct(do_gridcorr=.false.)
            even = self%even
            call self%even%set_cmat(cmat)
            deallocate(cmat)
            call even%ifft()
            call even%clip_inplace([self%box,self%box,self%box])
            call even%div(self%mag_correction)
            call even%write(add2fbody(fname_even,params_glob%ext,'_unfil'))
            ! odd
            cmat = self%odd%get_cmat()
            call self%odd%sampl_dens_correct(do_gridcorr=.false.)
            odd = self%odd
            call self%odd%set_cmat(cmat)
            deallocate(cmat)
            call odd%ifft()
            call odd%clip_inplace([self%box,self%box,self%box])
            call odd%div(self%mag_correction)
            call odd%write(add2fbody(fname_odd,params_glob%ext,'_unfil'))
            if( .not. l_have_fsc )then
                ! masking
                if( self%automsk )then
                    call even%mul(self%envmask)
                    call odd%mul(self%envmask)
                else if( trim(params_glob%automsk).eq.'yes' )then
                    call self%envmask%automask3D(even, odd, trim(params_glob%automsk).eq.'tight')
                    call even%zero_env_background(self%envmask)
                    call odd%zero_env_background(self%envmask)
                    call even%mul(self%envmask)
                    call odd%mul(self%envmask)
                    call self%envmask%write(string(MSKVOL_FILE))
                else
                    call even%mask(self%msk, 'soft', backgr=0.)
                    call odd%mask(self%msk, 'soft', backgr=0.)
                endif
                ! calculate FSC
                call even%fft()
                call odd%fft()
                call even%fsc(odd, self%fsc)
            endif
            ! regularization
            call self%even%add_invtausq2rho(self%fsc)
            call self%odd%add_invtausq2rho(self%fsc)
            ! Even: uneven sampling density correction, clip, & write
            cmat = self%even%get_cmat()
            call self%even%sampl_dens_correct(do_gridcorr=.false.)
            call self%even%ifft
            call even%zero_and_unflag_ft
            call self%even%clip(even)
            call even%div(self%mag_correction)
            call even%write(fname_even, del_if_exists=.true.)
            call self%even%set_cmat(cmat)
            call even%kill
            deallocate(cmat)
            ! Odd: uneven sampling density correction, clip, & write
            cmat = self%odd%get_cmat()
            call self%odd%sampl_dens_correct(do_gridcorr=.false.)
            call self%odd%ifft
            call odd%zero_and_unflag_ft
            call self%odd%clip(odd)
            call odd%div(self%mag_correction)
            call odd%write(fname_odd, del_if_exists=.true.)
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
            call even%div(self%mag_correction)
            call odd%div(self%mag_correction)
            ! write un-normalised unmasked even/odd volumes
            call even%write(fname_even, del_if_exists=.true.)
            call odd%write(fname_odd,   del_if_exists=.true.)
            if( .not. l_have_fsc )then
                ! masking
                if( self%automsk )then
                    call even%mul(self%envmask)
                    call odd%mul(self%envmask)
                else if( trim(params_glob%automsk).eq.'yes' )then
                    call self%envmask%automask3D(even, odd, trim(params_glob%automsk).eq.'tight')
                    call even%zero_env_background(self%envmask)
                    call odd%zero_env_background(self%envmask)
                    call even%mul(self%envmask)
                    call odd%mul(self%envmask)
                    call self%envmask%write(string(MSKVOL_FILE))
                else
                    call even%mask(self%msk, 'soft', backgr=0.)
                    call odd%mask(self%msk, 'soft', backgr=0.)
                endif
                ! calculate FSC
                call even%fft()
                call odd%fft()
                call even%fsc(odd, self%fsc)
            endif
        endif
        ! save, get & print resolution
        call arr2file(self%fsc, string(FSC_FBODY//int2str_pad(state,2)//BIN_EXT))
        call get_resolution(self%fsc, res, self%res_fsc05, self%res_fsc0143)
        self%res_fsc05   = max(self%res_fsc05,self%fny)
        self%res_fsc0143 = max(self%res_fsc0143,self%fny)
        ! Fourier index for eo averaging
        find4eoavg = max(K4EOAVGLB,  calc_fourier_index(FREQ4EOAVG3D,self%box,self%smpd))
        find4eoavg = min(find4eoavg, get_find_at_crit(self%fsc, FSC4EOAVG3D))
        deallocate(res)
        call even%kill
        call odd%kill
    end subroutine sampl_dens_correct_eos

    subroutine calc_fsc4sampl_dens_correct( self, even, odd, fsc )
        class(reconstructor_eo), intent(inout) :: self
        class(image),            intent(in)    :: even, odd
        real, allocatable,       intent(inout) :: fsc(:)
        type(image) :: even_tmp, odd_tmp
        if( allocated(fsc)      ) deallocate(fsc)
        if( allocated(self%fsc) ) deallocate(self%fsc)
        allocate(fsc(self%filtsz), source=0.)
        ! create temporary e/o:s
        call even_tmp%copy(even)
        call odd_tmp%copy(odd)
        ! masking
        if( self%automsk )then
            call even_tmp%mul(self%envmask)
            call odd_tmp%mul(self%envmask)
        else
            call even_tmp%mask(self%msk, 'soft', backgr=0.)
            call odd_tmp%mask(self%msk, 'soft', backgr=0.)
        endif
        ! calculate FSC
        call even_tmp%fft()
        call odd_tmp%fft()
        call even_tmp%fsc(odd_tmp, fsc)
        allocate(self%fsc(self%filtsz), source=fsc)
        call even_tmp%kill
        call odd_tmp%kill
    end subroutine calc_fsc4sampl_dens_correct

    subroutine write_fsc2txt( self, fname )
        class(reconstructor_eo), intent(in) :: self
        class(string),           intent(in) :: fname
        real, allocatable :: res(:)
        integer :: k, fnr
        if( .not. allocated(self%fsc) ) THROW_HARD('No FSC available to write to text file!')
        res = get_resarr(self%box, self%smpd)
        call fopen(fnr, FILE=fname, STATUS='REPLACE', action='WRITE')
        do k=1,size(res)
            write(fnr,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(k), '>>> CORRELATION:', self%fsc(k)
        end do
        write(fnr,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', self%res_fsc05
        write(fnr,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', self%res_fsc0143
        call fclose(fnr)
    end subroutine write_fsc2txt

    !> \brief  for sampling density correction, antialiasing, ifft & normalization of the sum
    subroutine sampl_dens_correct_sum( self, reference )
        class(reconstructor_eo), intent(inout) :: self      !< instance
        class(image),            intent(inout) :: reference !< reference volume
        if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> SAMPLING DENSITY (RHO) CORRECTION & WIENER NORMALIZATION'
        call reference%set_ft(.false.)
        call self%eosum%sampl_dens_correct(do_gridcorr=L_DO_GRIDCORR_GLOB)
        call self%eosum%ifft()
        call self%eosum%div(self%mag_correction)
        call self%eosum%clip(reference)
    end subroutine sampl_dens_correct_sum

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
            call self%envmask%kill_bimg
            call self%even%dealloc_rho
            call self%even%kill
            call self%odd%dealloc_rho
            call self%odd%kill
            call self%eosum%dealloc_rho
            call self%eosum%kill
            ! set existence
            self%exists = .false.
        endif
        if( allocated(self%fsc) ) deallocate(self%fsc)
    end subroutine kill

end module simple_reconstructor_eo
