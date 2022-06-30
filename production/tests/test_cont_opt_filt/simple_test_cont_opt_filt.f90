program simple_test_cont_opt_filt
    include 'simple_lib.f08'
    use simple_cmdline,            only: cmdline
    use simple_parameters,         only: parameters
    use simple_commander_volops,   only: reproject_commander
    use simple_image,              only: image
    use simple_cont_opt_filt
    use simple_optimizer,          only: optimizer
    use simple_opt_factory,        only: opt_factory
    use simple_opt_spec,           only: opt_spec
    implicit none
    type(parameters)              :: p
    type(cmdline)                 :: cline, cline_projection
    type(reproject_commander)     :: xreproject
    type(image)                   :: img, noise
    integer                       :: nptcls, iptcl, rc, ndim
    integer, parameter            :: NRESTARTS = 1
    character(len=:), allocatable :: cmd
    real,             allocatable :: cur_filt(:), lims(:,:), x_mat(:,:)
    logical                       :: mrc_exists
    real                          :: ave, sdev, maxv, minv, lowest_cost
    type(opt_factory)             :: ofac                 ! the optimization factory object
    type(opt_spec)                :: spec                 ! the optimizer specification object
    character(len=8)              :: str_opts             ! string descriptors for the NOPTS optimizers
    procedure(fun_interface), pointer :: costfun_ptr      !< pointer 2 cost function
    class(optimizer),         pointer :: opt_ptr=>null()  ! the generic optimizer object
    if( command_argument_count() < 4 )then
        write(logfhandle,'(a)') 'Usage: simple_test_cont_opt smpd=xx nthr=yy stk=stk.mrc mskdiam=zz'
        write(logfhandle,'(a)') 'Example: projections of https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
        write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
        inquire(file="1JYX.mrc", exist=mrc_exists)
        if( .not. mrc_exists )then
            write(*, *) 'Downloading the example dataset...'
            cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
            call execute_command_line(cmd, exitstat=rc)
            write(*, *) 'Converting .pdb to .mrc...'
            cmd = 'e2pdb2mrc.py 1JYX.pdb 1JYX.mrc'
            call execute_command_line(cmd, exitstat=rc)
            cmd = 'rm 1JYX.pdb'
            call execute_command_line(cmd, exitstat=rc)
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
    ndim      = product(p%ldim)
    call      img%new(p%ldim, p%smpd)
    call    noise%new(p%ldim, p%smpd)
    call  odd_img%new(p%ldim, p%smpd)
    call even_img%new(p%ldim, p%smpd)
    allocate(cur_filt(p%ldim(1)), lims(ndim,2), x_mat(p%ldim(1),p%ldim(2)), source=0.)
    lims(:,1) = calc_fourier_index(30., p%ldim(1), p%smpd)
    lims(:,2) = calc_fourier_index(10., p%ldim(1), p%smpd)
    write(*, *) calc_fourier_index(30., p%ldim(1), p%smpd), calc_fourier_index(10., p%ldim(1), p%smpd)
    do iptcl = 1, p%nptcls
        write(*, *) 'Particle # ', iptcl
        call img%read(p%stk, iptcl)
        call odd_img%copy_fast(img)
        call img%stats('foreground', ave, sdev, maxv, minv)
        ! addding noise
        call noise%gauran(0., 0.5 * sdev)
        call noise%mask(2. * p%msk, 'soft')
        call img%add(noise)
        call img%write('stk_noisy.mrc', iptcl)
        call even_img%copy_fast(img)
        ! do the optimization here to get the optimized cut-off frequency
        write(*, *) 'Cut-off frequency optimization in progress:'
        costfun_ptr  => filt_cost
        str_opts  = 'lbfgsb'
        call spec%specify(str_opts, ndim, limits=lims, nrestarts=NRESTARTS) ! make optimizer spec
        call spec%set_costfun(costfun_ptr)                                  ! set pointer to costfun
        call spec%set_gcostfun(filt_gcost)                                  ! set pointer to gradient of costfun
        call ofac%new(spec, opt_ptr)                                        ! generate optimizer object with the factory
        spec%x = (lims(:,1) + lims(:,2))/2.                                 ! set initial guess
        call opt_ptr%minimize(spec, opt_ptr, lowest_cost)                   ! minimize the test function
        x_mat = reshape(spec%x, [p%ldim(1), p%ldim(2)])
        write(*, *) 'cost = ', lowest_cost, '; x = ', x_mat(p%ldim(1)/2-2:p%ldim(1)/2+2, p%ldim(2)/2-2:p%ldim(2)/2+2)
        call butterworth_filter(cur_filt, 8, x_mat(p%ldim(1)/2, p%ldim(2)/2))
        call even_img%apply_filter(cur_filt)
        call even_img%write('cont_opt_filt_out.mrc', iptcl)
    enddo
end program simple_test_cont_opt_filt