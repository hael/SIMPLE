! testing nonuniform filter using probabilistic optimization (GA, Bayesian)
program simple_test_prob_opt_filt
    include 'simple_lib.f08'
    use simple_cmdline,            only: cmdline
    use simple_parameters,         only: parameters
    use simple_commander_volops,   only: reproject_commander
    use simple_image,              only: image
    use simple_opt_filter,         only: butterworth_filter
    use simple_opt_genetic
    implicit none
    type(parameters)              :: p
    type(cmdline)                 :: cline, cline_projection
    type(reproject_commander)     :: xreproject
    type(image)                   :: img, noise, odd_img, even_img, img_ker, img_pad
    type(image),      allocatable :: filt_img(:)
    integer                       :: npixels, iptcl, rc, ndim
    integer, parameter            :: NRESTARTS = 1, BW_ORDER = 8
    character(len=:), allocatable :: cmd
    real,             allocatable :: cur_filt(:)
    logical                       :: mrc_exists
    real                          :: ave, sdev, maxv, minv, lowest_cost
    integer,          allocatable :: best(:), x_mat(:,:)
    integer                       :: max_iter, pop_size, k, nspace, int_ind, k1, l1, sh
    real                          :: cross_rate, mut_rate, bounds(2), real_ind
    procedure(objective_func), pointer :: obj_func => null()
    if( command_argument_count() < 4 )then
        write(logfhandle,'(a)') 'Usage: simple_test_prob_opt_filt smpd=xx nthr=yy stk=stk.mrc mskdiam=zz'
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
    call find_ldim_nptcls(p%stk, p%ldim, p%nptcls)
    p%ldim(3) = 1 ! because we operate on stacks
    ndim      = product(p%ldim)
    call      img%new(p%ldim, p%smpd)
    call    noise%new(p%ldim, p%smpd)
    call  odd_img%new(p%ldim, p%smpd)
    call even_img%new(p%ldim, p%smpd)
    allocate(cur_filt(p%ldim(1)),        source=0.)
    allocate(x_mat(p%ldim(1),p%ldim(2)), source=0)
    nspace      = 40
    npixels     = p%ldim(1)*p%ldim(2)       ! number of pixels
    max_iter    = 20
    pop_size    = 10                        ! the population size is a magnitude of the number of pixels
    cross_rate  = 0.9
    mut_rate    = 1./npixels
    obj_func    => objective_cont
    allocate(filt_img(nspace))
    allocate(best(npixels), source=0)
    call srand(time())
    bounds      = [calc_fourier_index(30.        , p%ldim(1), p%smpd),&
                  &calc_fourier_index(2. * p%smpd, p%ldim(1), p%smpd)]
    print *, bounds
    print *, p%ldim
    ! caching all the filters (nspace of filters)
    do k = 1, nspace
        real_ind = bounds(1) + (k - 1.)*(bounds(2) - bounds(1))/(nspace - 1.)
        int_ind  = nint(real_ind)
        call filt_img(k)%new([2*int_ind, 2*int_ind, 1], p%smpd)
        call butterworth_filter(cur_filt, BW_ORDER, real_ind)
        call filt_img(k)%set_ft(.true.)
        call filt_img(k)%set_cmat((1., 0.))
        do l1 = -int_ind, int_ind-1
            do k1 = 0, int_ind
                sh = nint(hyp(real(k1),real(l1),0.)*p%ldim(1)/(2*int_ind))
                if( sh == 0 )then 
                    call filt_img(k)%mul([k1,l1,0], maxval(cur_filt))
                elseif( sh <= p%ldim(1) )then
                    call filt_img(k)%mul([k1,l1,0], cur_filt(sh))
                else
                    call filt_img(k)%mul([k1,l1,0], 0.)
                endif
            enddo
        enddo
        call filt_img(k)%ifft()
    enddo
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
        call img_pad%new([p%ldim(1) + 2*int(bounds(2)), p%ldim(2) + 2*int(bounds(2)),1], p%smpd)
        call img_ker%copy(even_img)
        call img_ker%pad(img_pad)
        ! do the optimization here to get the optimized cut-off frequency
        write(*, *) 'Cut-off frequency optimization in progress:'
        call genetic_opt(obj_func, npixels, nspace, max_iter, pop_size, cross_rate, mut_rate, best, lowest_cost)
        write(*, *) 'cost = ', lowest_cost, '; x = ', x_mat(p%ldim(1)/2-2:p%ldim(1)/2+2, p%ldim(2)/2-2:p%ldim(2)/2+2)
        call butterworth_filter(cur_filt, 8, real(to_ind(x_mat(p%ldim(1)/2, p%ldim(2)/2))))
        call even_img%apply_filter(cur_filt)
        call even_img%write('cont_opt_filt_out.mrc', iptcl)
    enddo
contains
    function to_ind(bit) result(val)
        integer, intent(in) :: bit
        integer             :: val
        val = int(bounds(1) + bit*(bounds(2) - bounds(1))/(nspace - 1.))
    end function to_ind

    function objective_cont(bitstring) result(val)
        integer,     intent(in) :: bitstring(:)
        real,        pointer    :: img_rmat(:,:,:), odd_rmat(:,:,:), filt_rmat(:,:,:)
        real                    :: val
        integer                 :: k, l, ldim(3), ext, k_tmp, l_tmp
        ldim  = odd_img%get_ldim()
        x_mat = reshape(bitstring, [ldim(1), ldim(2)])
        val   = 0.
        call odd_img%get_rmat_ptr(odd_rmat)
        do l = 1,ldim(2)
            do k = 1,ldim(1)
                ext   = to_ind(x_mat(k,l))
                k_tmp = k+int(bounds(2) - ext)
                l_tmp = l+int(bounds(2) - ext)
                call img_pad%get_rmat_ptr(img_rmat)
                call filt_img(x_mat(k,l)+1)%get_rmat_ptr(filt_rmat)
                val = val + abs(sum(filt_rmat(1:2*ext, 1:2*ext,1)*img_rmat(k_tmp:k_tmp+2*ext-1,l_tmp:l_tmp+2*ext-1,1)) - odd_rmat(k,l,1))**2
            enddo
        enddo
    end function objective_cont
end program simple_test_prob_opt_filt