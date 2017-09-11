! polar 2D Fourier transform generation by convolution interpolation (gridding)
#include "simple_lib.f08"
module simple_polarizer
!$ use omp_lib
!$ use omp_lib_kinds
use simple_defs        ! use all in there
use simple_syslib,     only: alloc_errchk
use simple_kbinterpol, only: kbinterpol
use simple_image,      only: image
use simple_params,     only: params

implicit none

public :: polarizer
private

type, extends(image) :: polarizer
    private
    type(kbinterpol)      :: kbwin                               !< window function object
    real,    allocatable  :: polweights_mat(:,:,:)               !< polar weights matrix for the image to polar transformer
    integer, allocatable  :: polcyc1_mat(:,:,:)                  !< image cyclic adresses for the image to polar transformer
    integer, allocatable  :: polcyc2_mat(:,:,:)                  !< image cyclic adresses for the image to polar transformer
    real                  :: winsz      = KBWINSZ                !< window half-width
    real                  :: alpha      = KBALPHA                !< oversampling ratio
    real                  :: harwin     = real(ceiling(KBWINSZ)) !< rounded window half-width
    real                  :: harwin_exp = 1.0                    !< rounded window half-width in expanded routines
    integer               :: wdim       = 2*ceiling(1.0) + 1     !< harwin_exp is argument to ceiling!< win dim
  contains
    procedure :: init_polarizer
    procedure :: polarize
    procedure :: kill_polarizer
end type polarizer

contains

    ! IMAGE TO POLAR FT TRANSFORMER

    !> \brief  initialises the image polarizer
    subroutine init_polarizer( self, pftcc )
        use simple_math,             only: sqwin_2d, cyci_1d
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(polarizer),        intent(inout) :: self   !< projector instance
        class(polarft_corrcalc), intent(inout) :: pftcc  !< polarft_corrcalc object to be filled
        real, allocatable :: w(:,:)
        real              :: loc(2)
        integer           :: pdim(3), win(2,2), lims(3,2)
        integer           :: i, k, l, wlen, cnt
        if( .not. pftcc%exists() ) stop 'polarft_corrcalc object needs to be created; init_imgpolarizer; simple_projector'
        call self%kill_polarizer
        self%kbwin = kbinterpol(KBWINSZ, KBALPHA)
        wlen       = self%wdim**2
        pdim       = pftcc%get_pdim(.true.)
        lims       = self%loop_lims(3)
        allocate( self%polcyc1_mat(1:pdim(1), pdim(2):pdim(3), 1:self%wdim),&
                  &self%polcyc2_mat(1:pdim(1), pdim(2):pdim(3), 1:self%wdim),&
                  &self%polweights_mat(1:pdim(1), pdim(2):pdim(3), 1:wlen),&
                  &w(1:self%wdim,1:self%wdim), stat=alloc_stat)
        if(alloc_stat /= 0) allocchk('in simple_projector :: init_imgpolarizer')
        !$omp parallel do collapse(2) schedule(static) default(shared)&
        !$omp private(i,k,l,w,loc,cnt,win) proc_bind(close)
        do i=1,pdim(1)
            do k=pdim(2),pdim(3)
                ! polar coordinates
                loc = pftcc%get_coord(i,k)
                win = sqwin_2d(loc(1), loc(2), self%harwin_exp)
                w   = 1.
                cnt = 0
                do l=1,self%wdim
                    cnt = cnt + 1
                    ! interpolation weights
                    w(l,:) = w(l,:) * self%kbwin%apod( real(win(1,1)+l-1)-loc(1) )
                    w(:,l) = w(:,l) * self%kbwin%apod( real(win(2,1)+l-1)-loc(2) )
                    ! cyclic addresses
                    self%polcyc1_mat(i, k, cnt) = cyci_1d(lims(1,:), win(1,1)+l-1)  !! last system error is temp created here
                    self%polcyc2_mat(i, k, cnt) = cyci_1d(lims(2,:), win(2,1)+l-1)
                end do
                self%polweights_mat(i,k,:) = reshape(w,(/wlen/))
            enddo
        enddo
        !$omp end parallel do
        deallocate(w)
    end subroutine init_polarizer

    !> \brief  creates the polar Fourier transform
    subroutine polarize( self, pftcc, img_ind, isptcl )
        use simple_math, only: sqwin_2d, cyci_1d
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(polarizer),        intent(inout) :: self   !< projector instance
        class(polarft_corrcalc), intent(inout) :: pftcc  !< polarft_corrcalc object to be filled
        integer,                 intent(in)    :: img_ind !< image index
        logical, optional,       intent(in)    :: isptcl !< is the input in polarised coords
        complex, allocatable :: pft(:,:), comps(:,:)
        integer :: i, k, l, m, windim, vecdim, addr_l
        integer :: lims(3,2), ldim_img(3), ldim_pft(3), pdim(3), logi(3), phys(3)
        logical :: iisptcl
        if( .not. allocated(self%polweights_mat) )&
        &stop 'the imgpolarizer has not been initialized!; simple_projector :: imgpolarizer'
        iisptcl = .true.
        if( present(isptcl) ) iisptcl = isptcl
        ldim_img = self%get_ldim()
        if( ldim_img(3) > 1 )      stop 'only for interpolation from 2D images; imgpolarizer; simple_projector'
        if( .not. pftcc%exists() ) stop 'polarft_corrcalc object needs to be created; imgpolarizer; simple_projector'
        ldim_pft = pftcc%get_ldim()
        if( .not. all(ldim_img == ldim_pft) )then
            print *, 'ldim_img: ', ldim_img
            print *, 'ldim_pft: ', ldim_pft
            stop 'logical dimensions do not match; imgpolarizer; simple_projector'
        endif
        if( .not.self%is_ft() ) stop 'image needs to FTed before this operation; simple_projector :: imgpolarizer'
        pdim   = pftcc%get_pdim(iisptcl)
        windim = 2*ceiling(self%harwin_exp) + 1
        vecdim = windim**2
        allocate( pft(pdim(1),pdim(2):pdim(3)), comps(1:windim,1:windim), stat=alloc_stat )
        if(alloc_stat /= 0) allocchk("In: imgpolarizer; simple_projector")
        lims = self%loop_lims(3)
        !$omp parallel do collapse(2) schedule(static) default(shared)&
        !$omp private(i,k,l,m,logi,phys,comps,addr_l) proc_bind(close)
        do i=1,pdim(1)
            do k=pdim(2),pdim(3)
                do l=1,windim
                    addr_l = self%polcyc1_mat(i,k,l)
                    do m=1,windim
                        logi = [addr_l,self%polcyc2_mat(i,k,m),0]
                        phys = self%comp_addr_phys(logi)
                        comps(l,m) = self%get_fcomp(logi,phys)
                    enddo
                enddo
                pft(i,k) = dot_product(self%polweights_mat(i,k,:), reshape(comps,(/vecdim/)))
            end do
        end do
        !$omp end parallel do
        if( iisptcl )then
            call pftcc%set_ptcl_pft(img_ind, pft)
        else
            call pftcc%set_ref_pft(img_ind, pft)
        endif
        ! kill the remains
        deallocate(pft, comps)
    end subroutine polarize

    ! DESTRUCTOR

    !>  \brief  is a destructor of impolarizer
    subroutine kill_polarizer( self )
        class(polarizer), intent(inout) :: self !< projector instance
        if( allocated(self%polweights_mat) ) deallocate(self%polweights_mat)
        if( allocated(self%polcyc1_mat)    ) deallocate(self%polcyc1_mat)
        if( allocated(self%polcyc2_mat)    ) deallocate(self%polcyc2_mat)
    end subroutine kill_polarizer

end module simple_polarizer
