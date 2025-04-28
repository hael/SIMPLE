program cavgs_simmat 
include 'simple_lib.f08'
    use simple_image
    use simple_corrmat
    use simple_parameters
    use simple_cmdline
    character(len = 255)    :: dir1 = '/Users/atifao/test_frac_min=0.2ncls_spec=5_2025-04-23_12:37:01/rank2_cavgs.mrc'
    integer     :: ldim(3), n_cavgs, i, N, nthr
    real        :: smpd = 2., hp = 60., lp = 10.
    real, allocatable           :: corrmat(:,:), R(:,:), X(:,:), Y(:,:)
    logical, allocatable        :: M(:,:)
    type(image), allocatable    :: cavgs(:), imgs4sim(:)
    type(parameters)            :: params
    type(cmdline)               :: cline
    ! params
    
    call find_ldim_nptcls(dir1, ldim, n_cavgs)
    ldim = [ldim(1), ldim(2), 1]
    N = ldim(1)

    allocate(cavgs(n_cavgs))
    do i = 1, n_cavgs 
        call cavgs(i)%new(ldim, smpd)
        call cavgs(i)%read(dir1, i)
    end do 
    allocate(imgs4sim(10))
    
    ! Apoferritin View 1 monomer 
    call imgs4sim(1)%copy(cavgs(69))
    call imgs4sim(2)%copy(cavgs(68))
    call imgs4sim(3)%copy(cavgs(64))
    call imgs4sim(4)%copy(cavgs(63))
    call imgs4sim(5)%copy(cavgs(61))
    ! Appoferritin View 2 dimer 
    call imgs4sim(6)%copy(cavgs(67))
    call imgs4sim(7)%copy(cavgs(66))
    call imgs4sim(8)%copy(cavgs(65))
    call imgs4sim(9)%copy(cavgs(59))
    call imgs4sim(10)%copy(cavgs(58))

    ! set params 
    call cline%set('ctf',   'no')
    call cline%set('objfun','cc')
    call cline%set('sh_inv', 'yes')
    call cline%set('box', N)
    call params%new(cline)

    nthr = 11
    call calc_inplane_fm( imgs4sim, hp, lp, nthr, corrmat )

    print *, corrmat
end program 

