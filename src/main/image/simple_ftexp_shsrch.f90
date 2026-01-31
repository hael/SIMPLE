!@descr: shift search with L-BFGS-B using expanded Fourier transforms (used in motion_correct)
module simple_ftexp_shsrch
use simple_core_module_api
use simple_opt_spec,    only: opt_spec
use simple_optimizer,   only: optimizer
use simple_ft_expanded, only: ft_expanded
implicit none

public :: ftexp_shsrch, test_ftexp_shsrch, test_ftexp_shsrch2
private
#include "simple_local_flags.inc"

real,     parameter :: TOL    = 1e-4    !< tolerance parameter
integer,  parameter :: MAXITS = 30      !< maximum number of iterations

type :: ftexp_shsrch
    private
    type(opt_spec), public      :: ospec                     !< optimizer specification object
    class(optimizer),   pointer :: opt_obj      => null()    !< pointer to nonlinear optimizer
    class(ft_expanded), pointer :: reference    => null()    !< reference ft_exp
    class(ft_expanded), pointer :: particle     => null()    !< particle ft_exp
    real(dp)                    :: denominator        = 0.d0
    real                        :: maxHWshift         = 0.   !< maximum half-width of shift
    real                        :: motion_correctftol = 1e-4 !< function error tolerance
    real                        :: motion_correctgtol = 1e-4 !< gradient error tolerance
    real                        :: shsrch_tol         = TOL
    logical                     :: existence = .false.
  contains
    procedure          :: new            => ftexp_shsrch_new
    procedure          :: minimize       => ftexp_shsrch_minimize
    procedure          :: corr_shifted_8 => ftexp_shsrch_corr_shifted_8
    procedure          :: kill           => ftexp_shsrch_kill
    procedure          :: set_shsrch_tol
    procedure          :: set_factr_pgtol
end type ftexp_shsrch

contains

    !> Initialise  ftexp_shsrch
    subroutine ftexp_shsrch_new( self, ref, ptcl, trs, motion_correct_ftol, motion_correct_gtol )
        use simple_opt_factory, only: opt_factory
        class(ftexp_shsrch),        intent(inout) :: self
        class(ft_expanded), target, intent(in)    :: ref, ptcl
        real,                       intent(in)    :: trs
        real,             optional, intent(in)    :: motion_correct_ftol, motion_correct_gtol
        type(opt_factory) :: ofac
        real              :: opt_lims(2,2)
        call self%kill()
        self%reference  => ref
        self%particle   => ptcl
        self%maxHWshift =  trs
        if( present(motion_correct_ftol) )then
            self%motion_correctftol = motion_correct_ftol
        else
            self%motion_correctftol = TOL
        end if
        if( present(motion_correct_gtol) )then
            self%motion_correctgtol = motion_correct_gtol
        else
            self%motion_correctgtol = TOL
        end if
        opt_lims(1,1) = - self%maxHWshift
        opt_lims(1,2) =   self%maxHWshift
        opt_lims(2,1) = - self%maxHWshift
        opt_lims(2,2) =   self%maxHWshift
        call self%ospec%specify('lbfgsb', 2, ftol=self%motion_correctftol, gtol=self%motion_correctgtol, limits=opt_lims)
        call self%ospec%set_costfun_8(ftexp_shsrch_cost_8)
        call self%ospec%set_gcostfun_8(ftexp_shsrch_gcost_8)
        call self%ospec%set_fdfcostfun_8(ftexp_shsrch_fdfcost_8)
        ! generate optimizer object with the factory
        if( associated(self%opt_obj) )then
            call self%opt_obj%kill
            deallocate(self%opt_obj) ! because this extended type is allocated by opt_factory
            nullify(self%opt_obj)
        end if
        call ofac%new(self%ospec, self%opt_obj)
        self%existence = .true.
    end subroutine ftexp_shsrch_new

    !> Main search routine
    function ftexp_shsrch_minimize( self, prev_corr, prev_shift ) result( cxy )
        class(ftexp_shsrch), intent(inout) :: self
        real, optional,      intent(in)    :: prev_corr, prev_shift(2)
        real :: cxy(3)
        self%ospec%limits(1,1) = - self%maxHWshift
        self%ospec%limits(1,2) =   self%maxHWshift
        self%ospec%limits(2,1) = - self%maxHWshift
        self%ospec%limits(2,2) =   self%maxHWshift
        if( present(prev_shift) )then
            self%ospec%limits(1,:) = self%ospec%limits(1,:) + prev_shift(1)
            self%ospec%limits(2,:) = self%ospec%limits(2,:) + prev_shift(2)
        endif
        if( present(prev_shift) ) then
            self%ospec%x = prev_shift
        else
            self%ospec%x   = 0.
        end if
        self%ospec%x_8 = real(self%ospec%x,dp)
        ! self%kind_shift = self%reference%get_kind_shift()
        ! call self%set_dims_and_alloc()
        ! call self%calc_tmp_cmat12()
        ! the temp matrix is built on the particle only!
        call self%particle%alloc_and_calc_tmp_cmat12(self%reference, self%denominator )
        ! set initial solution to previous shift
        call self%opt_obj%minimize(self%ospec, self, cxy(1))
        call self%reference%corr_normalize(self%particle, cxy(1))
        cxy(1)  = -cxy(1) ! correlation
        cxy(2:) = self%ospec%x ! shift
        if( present(prev_corr) )then
            if( abs(cxy(1)-prev_corr) <= self%shsrch_tol )then
                cxy(1)  = prev_corr
                if( present(prev_shift) ) cxy(2:) = prev_shift
            endif
        endif
    end function ftexp_shsrch_minimize

    subroutine ftexp_shsrch_kill( self )
        class(ftexp_shsrch), intent(inout) :: self
        if ( self%existence ) then
            call self%ospec%kill
            if( associated( self%opt_obj ) )then
                call self%opt_obj%kill
                deallocate(self%opt_obj)
                nullify(self%opt_obj)
            end if
            if( associated(self%reference) ) self%reference => null()
            if( associated(self%particle)  )then
                call self%particle%dealloc_tmp_cmat12
                self%particle  => null()
            endif
            self%existence = .false.
        end if
    end subroutine ftexp_shsrch_kill

    pure subroutine set_shsrch_tol( self, shsrch_tol )
        class(ftexp_shsrch), intent(inout) :: self
        real,                intent(in)    :: shsrch_tol
        self%shsrch_tol = shsrch_tol
    end subroutine set_shsrch_tol

    pure subroutine set_factr_pgtol( self, factr, pgtol )
        class(ftexp_shsrch), intent(inout) :: self
        real(dp),            intent(in)    :: factr, pgtol
        self%ospec%factr = factr
        self%ospec%pgtol = pgtol
    end subroutine set_factr_pgtol

    ! Correlation
    real(dp) function ftexp_shsrch_corr_shifted_8( self, shvec )
        class(ftexp_shsrch), intent(inout) :: self
        real(dp),            intent(in)    :: shvec(2)
        call self%particle%alloc_and_calc_tmp_cmat12(self%reference, self%denominator)
        ftexp_shsrch_corr_shifted_8 = self%particle%corr_shifted_cost_8(shvec, self%denominator)
        call self%particle%corr_normalize(self%reference, ftexp_shsrch_corr_shifted_8)
    end function ftexp_shsrch_corr_shifted_8

    !> Cost function, double precision
    function ftexp_shsrch_cost_8( self, vec, D ) result( cost )
        class(*),     intent(inout) :: self
        integer,      intent(in)    :: D
        real(kind=8), intent(in)    :: vec(D)
        real(kind=8) :: cost
        select type(self)
            class is (ftexp_shsrch)
                cost = -ftexp_shsrch_corr_shifted_8(self, -vec)
            class DEFAULT
                THROW_HARD('unknown type; ftexp_shsrch_cost_8')
        end select
    end function ftexp_shsrch_cost_8

    !> Gradient function, double precision
    subroutine ftexp_shsrch_gcost_8( self, vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        grad = 0.d0
        select type(self)
            class is (ftexp_shsrch)
                call self%particle%corr_gshifted_cost_8( -vec, self%denominator, grad )
            class DEFAULT
                THROW_HARD('unknown type; ftexp_shsrch_gcost_8')
        end select
    end subroutine ftexp_shsrch_gcost_8

    !> Gradient & cost function, double precision
    subroutine ftexp_shsrch_fdfcost_8( self, vec, f, grad, D )
        class(*),     intent(inout) :: self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: f, grad(D)
        f    = 0.d0
        grad = 0.d0
        select type(self)
            class is (ftexp_shsrch)
                call self%particle%corr_fdfshifted_cost_8( -vec, self%denominator, f, grad )
                f = f * (-1.0_dp)
            class DEFAULT
                THROW_HARD('unknown type; ftexp_shsrch_fdfcost_8')
        end select
    end subroutine ftexp_shsrch_fdfcost_8

    ! TEST ROUTINES

    ! Tests correlations/gradients routines & ft_expand interface
    subroutine test_ftexp_shsrch
        use simple_ft_expanded, only: ftexp_transfmat_init
        use simple_image,       only: image
        ! global constants
        integer, parameter :: LDIM(3)=[240,240,1], SQRAD=60, NTST=50, NNOISY=20
        real,    parameter :: SMPD=1.77, TRS=10., HP=100.0, LP=8., SNR=0.2
        ! global variables
        type(ft_expanded)        :: ftexp_img
        type(image)              :: img, img_shifted
        type(image), allocatable :: noisy_imgs(:)
        integer                  :: x, y
        integer :: i
        ! reference images
        call img%new(LDIM, SMPD)
        call img%square(SQRAD)
        call ftexp_img%new(img, HP, LP, .true.)
        call img_shifted%new(LDIM, SMPD)
        ! seed the random number generator
        call seed_rnd
        ! generate noisy images
        if( allocated(noisy_imgs) )then
            do i=1,NNOISY
                call noisy_imgs(i)%kill
            end do
            deallocate(noisy_imgs)
        else
            allocate(noisy_imgs(NNOISY))
            do i=1,NNOISY
                noisy_imgs(i) = img
                call noisy_imgs(i)%add_gauran(SNR)
            end do
        endif
        call test_shifted_correlator
        call profile_corrs
        contains

            subroutine test_shifted_correlator
                type(ft_expanded)  :: ftexp_trial
                type(ftexp_shsrch) :: ftexp_shsrch_trial
                real(dp) :: dbl_corr, grad(2), grad_best(2), vec(2)
                real     :: dist, corr_best, corravg, distavg, diff_corr
                real     :: corr, diff, mag, diff_corr_avg, grad_mag_avg
                integer  :: ysh_best, xsh_best, xsh, ysh, itst
                write(logfhandle,*) 'testing ft_expanded/ftexp_shsrch :: shifted correlator'
                call ftexp_transfmat_init(img, LP)
                corravg       = 0.
                distavg       = 0.
                diff_corr_avg = 0.
                grad_mag_avg  = 0.
                do itst=1,NTST
                    x = nint(ran3()*2*TRS-TRS)
                    y = nint(ran3()*2*TRS-TRS)
                    img_shifted = img
                    call img_shifted%shift([real(x),real(y),0.])
                    ! single search ----
                    call ftexp_trial%new(img_shifted, HP, LP, .true.)
                    call ftexp_shsrch_trial%new(ftexp_img, ftexp_trial, TRS)
                    corr_best = -huge(corr)
                    grad_best = huge(dbl_corr)
                    do xsh=nint(-TRS),nint(TRS)
                        do ysh=nint(-TRS),nint(TRS)
                            vec  = dble([xsh, ysh])
                            corr = -real(ftexp_shsrch_cost_8(ftexp_shsrch_trial, vec, 2))
                            if( corr > corr_best )then
                                corr_best = corr
                                xsh_best  = xsh
                                ysh_best  = ysh
                            endif
                        end do
                    end do
                    dist = euclid(real([xsh_best,ysh_best]), real([-x,-y]))
                    vec  = dble([xsh_best,ysh_best])
                    call ftexp_shsrch_fdfcost_8(ftexp_shsrch_trial, vec, dbl_corr, grad, 2)
                    call ftexp_img%corr_normalize(ftexp_trial, dbl_corr)
                    dbl_corr = -dbl_corr
                    diff = abs(real(dbl_corr)-corr_best)
                    call ftexp_img%corr_normalize(ftexp_trial, grad(1))
                    call ftexp_img%corr_normalize(ftexp_trial, grad(2))
                    print *, x,y, xsh_best,ysh_best, corr_best, dbl_corr, grad
                    mag  = real(sqrt(sum(grad**2)))
                    call ftexp_trial%kill
                    call ftexp_shsrch_trial%kill
                    ! end single search ----
                    corravg       = corravg       + corr_best
                    distavg       = distavg       + dist
                    diff_corr_avg = diff_corr_avg + diff
                    grad_mag_avg  = grad_mag_avg  + mag
                end do
                corravg       = corravg/real(NTST)
                distavg       = distavg/real(NTST)
                diff_corr_avg = diff_corr_avg/real(NTST)
                grad_mag_avg  = grad_mag_avg/real(NTST)
                if( corravg > 0.999 .and. distavg < 0.0001 .and. diff_corr<0.0001 )then
                    write(logfhandle,'(a)')'>>> CORRELATIONS TEST PASSED'
                else
                    print *, 'corravg:   ', corravg
                    print *, 'distavg:   ', distavg
                    print *, 'diff_corr: ', diff_corr

                    THROW_HARD('****ft_expanded_tester FAILURE 1 :: test_shifted_correlator')
                endif
                if( grad_mag_avg < 0.000001 )then
                    write(logfhandle,'(a)')'>>> CORRELATIONS GRADIENTS TEST PASSED'
                else
                    THROW_HARD('****ft_expanded_tester FAILURE 2 :: test_shifted_correlator')
                endif
            end subroutine test_shifted_correlator

            subroutine profile_corrs
                integer, parameter   :: NTSTS=1000, NTHR=8
                integer              :: itst
                type(image)          :: img_ref, img_ptcl
                type(ft_expanded)    :: ftexp_ref, ftexp_ptcl
                type(ftexp_shsrch)   :: ftexp_shsrch1
                real, allocatable    :: shvecs(:,:)
                real(4)    :: corr, actual, delta, tarray(2)
                call img_ref%new([4096,4096,1],SMPD)
                call img_ref%ran
                call img_ref%fft()
                call ftexp_transfmat_init(img_ref, LP)
                call ftexp_ref%new(img_ref, HP, LP, .true.)
                call img_ptcl%new([4096,4096,1],SMPD)
                call img_ptcl%ran
                call img_ptcl%fft()
                call ftexp_ptcl%new(img_ptcl, HP, LP, .true.)
                call ftexp_shsrch1%new(ftexp_ref, ftexp_ptcl, TRS)
                allocate(shvecs(NTSTS,3))
                do itst=1,NTSTS
                    shvecs(itst,1) = ran3()*2*TRS-TRS
                    shvecs(itst,2) = ran3()*2*TRS-TRS
                    shvecs(itst,3) = 0.
                end do
                actual = etime( tarray )
                write(logfhandle,'(A,2X,F9.2)') 'Actual cpu-time:', actual
                delta = dtime( tarray )
                write(logfhandle,'(A,F9.2)') 'Relative cpu-time:', delta

                write(logfhandle,'(a)') '>>> PROFILING STANDARD CORRELATOR'
                do itst=1,NTST
                    corr = img_ref%corr_shifted(img_ptcl, shvecs(itst,:), lp_dyn=LP)
                end do
                actual = etime( tarray )
                write(logfhandle,'(A,2X,F9.2)') 'Actual cpu-time:', actual
                delta = dtime( tarray )
                write(logfhandle,'(A,F9.2)') 'Relative cpu-time:', delta
                write(logfhandle,'(a)') '>>> PROFILING FTEXP CORRELATOR'
                do itst=1,NTST
                    corr = real(ftexp_shsrch1%corr_shifted_8(dble(shvecs(itst,:))))
                end do
                actual = etime( tarray )
                write(logfhandle,'(A,2X,F9.2)') 'Actual cpu-time:', actual
                delta = dtime( tarray )
                write(logfhandle,'(A,F9.2)') 'Relative cpu-time:', delta
            end subroutine profile_corrs

    end subroutine test_ftexp_shsrch

    ! test optmization
    subroutine test_ftexp_shsrch2
        use simple_ft_expanded, only: ftexp_transfmat_init
        use simple_image,       only: image
        type(image)       :: img_ref, img_ptcl
        type(ft_expanded) :: ftexp_ref, ftexp_ptcl
        type(ftexp_shsrch)  :: ftexp_srch
        real, parameter   :: TRS=5.0, lp=6., hp=100.
        real              :: cxy(3), x, y, lims(2,2)
        integer           :: i
        write(logfhandle,*) 'testing ft_expanded :: optimization'
        lims(:,1) = -TRS
        lims(:,2) = TRS
        call img_ref%new([32,32,1], 2.)
        call img_ptcl%new([32,32,1], 2.)
        img_ref = 1.
        call img_ref%memoize_mask_coords
        call img_ref%mask2D_soft(8.,backgr=0.)
        call img_ref%fft()
        call ftexp_ref%new(img_ref, hp, lp, .true.)
        call ftexp_ptcl%new(img_ptcl, hp, lp, .false.)
        call ftexp_srch%new(ftexp_ref, ftexp_ptcl, trs)
        call ftexp_transfmat_init(img_ref,lp)
        do i=1,100
            x = ran3()*2*TRS-TRS
            y = ran3()*2*TRS-TRS
            img_ptcl = img_ref
            call img_ptcl%shift([x,y,0.])
            call ftexp_ptcl%new(img_ptcl, hp, lp, .true.)
            cxy = ftexp_srch%minimize()
            print *,i,x,y,cxy
            if( cxy(1) < 0.999 )then
                THROW_HARD('shift alignment failed; test_ftexp_shsrch')
            endif
        end do
        write(logfhandle,'(a)') 'SIMPLE_ftexp_shsrch_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_ftexp_shsrch2

end module simple_ftexp_shsrch
