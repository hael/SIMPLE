!@descr: for all optimize tests
module simple_commanders_test_optimize
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_lbfgsb
  contains
    procedure :: execute      => exec_test_lbfgsb
end type commander_test_lbfgsb

type, extends(commander_base) :: commander_test_lbfgsb_cosine
  contains
    procedure :: execute      => exec_test_lbfgsb_cosine
end type commander_test_lbfgsb_cosine

type, extends(commander_base) :: commander_test_lplims
  contains
    procedure :: execute      => exec_test_lplims
end type commander_test_lplims

type, extends(commander_base) :: commander_test_lpstages_test
  contains
    procedure :: execute      => exec_test_lpstages_test
end type commander_test_lpstages_test

type, extends(commander_base) :: commander_test_opt_lp
  contains
    procedure :: execute      => exec_test_opt_lp
end type commander_test_opt_lp


contains

subroutine exec_test_lbfgsb( self, cline )
    use simple_optimizer,   only: optimizer
    use simple_opt_factory, only: opt_factory
    use simple_opt_spec,    only: opt_spec
    class(commander_test_lbfgsb),    intent(inout) :: self
    class(cmdline),                  intent(inout) :: cline
    class(optimizer), pointer   :: opt_ptr=>null()      ! the generic optimizer object
    integer,          parameter :: NDIM=1, NRESTARTS=1
    type(opt_factory) :: ofac                           ! the optimization factory object
    type(opt_spec)    :: spec                           ! the optimizer specification object
    character(len=8)  :: str_opts                       ! string descriptors for the NOPTS optimizers
    real              :: lims(NDIM,2), lowest_cost
    str_opts  = 'lbfgsb'
    lims(1,1) = -5.
    lims(1,2) =  5.
    call spec%specify(str_opts, NDIM, limits=lims, nrestarts=NRESTARTS) ! make optimizer spec
    call spec%set_costfun(costfct)                                      ! set pointer to costfun
    call spec%set_gcostfun(gradfct)                                     ! set pointer to gradient of costfun
    call ofac%new(spec, opt_ptr)                                        ! generate optimizer object with the factory
    spec%x = 0.5
    call opt_ptr%minimize(spec, opt_ptr, lowest_cost)                   ! minimize the test function
    write(*, *) lowest_cost, spec%x
    call opt_ptr%kill
    deallocate(opt_ptr)
    call simple_end('**** SIMPLE_TEST_LBFGSB_WORKFLOW NORMAL STOP ****')

    contains

     function costfct( fun_self, x, d ) result( r )
         class(*), intent(inout) :: fun_self
         integer,  intent(in)    :: d
         real,     intent(in)    :: x(d)
         real                    :: r
         r = (x(1) - 1.)**2
     end function
 
     subroutine gradfct( fun_self, x, grad, d )
         class(*), intent(inout) :: fun_self
         integer,  intent(in)    :: d
         real,     intent(inout) :: x(d)
         real,     intent(out)   :: grad(d)
         grad(1) = 2. * (x(1) - 1.)
     end subroutine

end subroutine exec_test_lbfgsb

subroutine exec_test_lbfgsb_cosine( self, cline )
    use simple_optimizer,   only: optimizer
    use simple_opt_factory, only: opt_factory
    use simple_opt_spec,    only: opt_spec
    class(commander_test_lbfgsb_cosine), intent(inout) :: self
    class(cmdline),                      intent(inout) :: cline
    class(optimizer), pointer   :: opt_ptr=>null()      ! the generic optimizer object
    integer,          parameter :: NDIM=2, NRESTARTS=10
    type(opt_factory) :: ofac                           ! the optimization factory object
    type(opt_spec)    :: spec                           ! the optimizer specification object
    character(len=8)  :: str_opts                       ! string descriptors for the NOPTS optimizers
    real              :: lims(NDIM,2), lowest_cost, y_norm(2), y(2)
    str_opts  = 'lbfgsb'
    lims(:,1) = -10.
    lims(:,2) =  10.
    y         = [.75, .25]
    y_norm    = y / sqrt(sum(y**2))
    call spec%specify(str_opts, NDIM, limits=lims, nrestarts=NRESTARTS) ! make optimizer spec
    call spec%set_costfun(costfct)                                      ! set pointer to costfun
    call spec%set_gcostfun(gradfct)                                     ! set pointer to gradient of costfun
    call ofac%new(spec, opt_ptr)                                        ! generate optimizer object with the factory
    spec%x = [-5., -7.5]
    call opt_ptr%minimize(spec, opt_ptr, lowest_cost)                   ! minimize the test function
    write(*, *) lowest_cost, spec%x, y
    call opt_ptr%kill
    deallocate(opt_ptr)
    call simple_end('**** SIMPLE_TEST_LBFGSB_COSINE_WORKFLOW NORMAL STOP ****')

    contains
    
      function costfct( fun_self, x, d ) result( r )
          class(*), intent(inout) :: fun_self
          integer,  intent(in)    :: d
          real,     intent(in)    :: x(d)
          real                    :: r
          real :: x_norm(d)
          x_norm = x / sqrt(sum(x**2))
          r      = acos(sum(x_norm * y_norm))
          print *, 'cost = ', r
      end function
    
      subroutine gradfct( fun_self, x, grad, d )
          class(*), intent(inout) :: fun_self
          integer,  intent(in)    :: d
          real,     intent(inout) :: x(d)
          real,     intent(out)   :: grad(d)
          real :: abs_x, x_norm(d)
          abs_x  = sqrt(sum(x**2))
          x_norm = x / abs_x
          grad   = - (y_norm * abs_x - sum(x_norm * y_norm) * x) / abs_x**2 / sqrt(1. - sum(x_norm * y_norm)**2)
      end subroutine

end subroutine exec_test_lbfgsb_cosine

subroutine exec_test_lplims( self, cline )
    class(commander_test_lplims),    intent(inout) :: self
    class(cmdline),                  intent(inout) :: cline
    real :: mskdiam, lpstart,lpstop, lpcen
    mskdiam = 300.
    do while( mskdiam >= 90.)
        call mskdiam2lplimits( mskdiam, lpstart,lpstop, lpcen )
        print *, 'mskdiam/lpstart/lpstop/lpcen: ', mskdiam, lpstart, lpstop, lpcen
        mskdiam = mskdiam - 10
    end do
    call simple_end('**** SIMPLE_TEST_LPLIMS_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_lplims

subroutine exec_test_lpstages_test( self, cline )
    class(commander_test_lpstages_test), intent(inout) :: self
    class(cmdline),                      intent(inout) :: cline
    integer, parameter :: LDIM(3)=[256,256,1], BOX=LDIM(1), FILTSZ=BOX/2, NSTAGES=10
    real,    parameter :: SMPD=1.3, LPSTART_DEFAULT=20., LPSTART_LB=10., LPFINAL=6.
    real               :: frc(FILTSZ) = 1.
    type(lp_crop_inf)  :: lpinfo(NSTAGES)
    call lpstages(BOX, NSTAGES, frc, SMPD, LPSTART_LB, LPSTART_DEFAULT, LPFINAL, lpinfo, l_cavgs=.false.)
    call simple_end('**** SIMPLE_TEST_LPSTAGES_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_lpstages_test

subroutine exec_test_opt_lp( self, cline )
    use simple_cmdline,           only: cmdline
    use simple_builder,           only: builder
    use simple_parameters,        only: parameters
    use simple_commanders_volops, only: commander_reproject
    use simple_image,             only: image
    use simple_butterworth,       only: butterworth_filter
    use simple_math,              only: create_hist_vector
    use simple_commanders_atoms,  only: commander_pdb2mrc
    class(commander_test_opt_lp),    intent(inout) :: self
    class(cmdline),                  intent(inout) :: cline
    type(parameters)              :: p
    type(cmdline)                 :: cline_projection, cline_pdb2mrc
    type(image)                   :: img, noise, img_noisy, img_filt
    type(commander_reproject)     :: xreproject
    type(commander_pdb2mrc)       :: xpdb2mrc
    integer                       :: nptcls, rc, iptcl, find_stop, find_start, n_bin, n_vec, find_cur
    real                          :: ave, sdev, maxv, minv, noise_mean, noise_std
    type(string)                  :: cmd
    logical                       :: mrc_exists
    real,             allocatable :: cur_fil(:), vec_noise(:), xhist(:), yest(:)
    integer,          allocatable :: yhist(:)
    real,             pointer     :: rmat_img_noisy(:,:,:), rmat_img_filt(:,:,:)
    if( command_argument_count() < 4 )then
        write(logfhandle,'(a)') 'ERROR! Usage: simple_test_opt_lp smpd=xx nthr=yy stk=stk.mrc, mskdiam=zz'
        write(logfhandle,'(a)') 'Example: projections of https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
        write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
        inquire(file="1JYX.mrc", exist=mrc_exists)
        if( .not. mrc_exists )then
            write(*, *) 'Downloading the example dataset...'
            cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
            call execute_command_line(cmd%to_char(), exitstat=rc)
            write(*, *) 'Converting .pdb to .mrc...'
            call cline_pdb2mrc%set('smpd',                            1.)
            call cline_pdb2mrc%set('pdbfile',                 '1JYX.pdb')
            call cline_pdb2mrc%checkvar('smpd',                        1)
            call cline_pdb2mrc%checkvar('pdbfile',                     2)
            call cline_pdb2mrc%check()
            call xpdb2mrc%execute(cline_pdb2mrc)
            call cline_pdb2mrc%kill()
            cmd = 'rm 1JYX.pdb'
            call execute_command_line(cmd%to_char(), exitstat=rc)
            write(*, *) 'Projecting 1JYX.mrc...'
            call cline_projection%set('vol1'      , '1JYX.mrc')
            call cline_projection%set('smpd'      , 1.)
            call cline_projection%set('pgrp'      , 'c1')
            call cline_projection%set('mskdiam'   , 180.)
            call cline_projection%set('nspace'    , 6.)
            call cline_projection%set('nthr'      , 16.)
            call xreproject%execute(cline_projection)
        endif
        call cline%set('smpd'   , 1.)
        call cline%set('nthr'   , 16.)
        call cline%set('stk'    , 'reprojs.mrcs')
        call cline%set('mskdiam', 180.)
    else
        call cline%parse_oldschool
    endif
    call cline%checkvar('smpd',    1)
    call cline%checkvar('nthr',    2)
    call cline%checkvar('stk' ,    3)
    call cline%checkvar('mskdiam', 4)
    call cline%check
    call p%new(cline)
    call find_ldim_nptcls(p%stk, p%ldim, nptcls)
    p%ldim(3) = 1 ! because we operate on stacks
    n_vec     = p%ldim(1)*p%ldim(2)
    n_bin     = int(n_vec/1000.)
    print *, n_vec, n_bin
    call img%new(  p%ldim, p%smpd)
    call noise%new(p%ldim, p%smpd)
    allocate(cur_fil(p%ldim(1)), source=0.)
    allocate(vec_noise(n_vec),   source=0.)
    allocate(xhist(n_bin),       source=0.)
    allocate(yest( n_bin),       source=0.)
    allocate(yhist(n_bin),       source=0)
    find_stop  = calc_fourier_index(p%lpstart,   p%ldim(1), p%smpd)
    find_start = calc_fourier_index(p%lpstart_nonuni, p%ldim(1), p%smpd)
    call img%memoize_mask_coords
    do iptcl = 1, 1
        write(*, *) 'Particle # ', iptcl
        cur_fil = 0.
        call img%read(p%stk, iptcl)
        call img_noisy%copy(img)
        ! spherical masking
        call img%mask2D_soft(p%msk)
        ! img stats
        call img%stats('foreground', ave, sdev, maxv, minv)
        ! add noise in a small center region of the even
        call noise%gauran(0., .2 * sdev)
        call noise%mask2D_soft(1.5 * p%msk)
        call img_noisy%add(noise)
        call img_noisy%write(string('stk_noisy.mrc'), iptcl)
        call img_noisy%get_rmat_ptr(rmat_img_noisy)
        ! find_cur = int((find_start + find_stop)/2.)
        ! find_cur = int(find_start + 50)
        do find_cur = find_start + 10, find_stop - 10, 5
            call img_filt%copy(img_noisy)
            call img_filt%fft
            call butterworth_filter(find_cur - 5, find_cur + 5, cur_fil)
            call img_filt%apply_filter(cur_fil)
            call img_filt%ifft
            call img_filt%add(img_noisy)
            call img_filt%write(string('stk_filt.mrc'), iptcl)
            call img_filt%get_rmat_ptr( rmat_img_filt )
            vec_noise = [rmat_img_noisy] - [rmat_img_filt]
            call avg_sdev(vec_noise, noise_mean, noise_std)
            call create_hist_vector( vec_noise, n_bin, xhist, yhist )
            yest = real(yhist)
            yest = exp( -(yest-noise_mean)**2/2./noise_std**2 )/2./pi/noise_std
            yest = yest * sum(xhist*yhist) / sum(xhist*yest)
            print *, sum( (yest - real(yhist))**2 )
            call img_filt%zero_and_unflag_ft
        enddo
        ! print *, yhist
        call img%zero_and_unflag_ft
        call img_noisy%zero_and_unflag_ft
        call noise%zero_and_unflag_ft
    enddo
    call img%kill()
    call img_noisy%kill()
    call noise%kill()
    call simple_end('**** SIMPLE_TEST_OPT_LP_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_opt_lp


end module simple_commanders_test_optimize
