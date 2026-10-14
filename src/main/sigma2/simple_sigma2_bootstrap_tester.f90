!@descr: unit tests of the sigma2 bootstrap (simple_sigma2_bootstrap)
! prepare_pspec_cline builds calc_pspec's command line fresh, with the few keys
! it takes from its template; prepare_residual_sigma2_pass_cline keeps the
! template's scoring keys and drops what samples, assembles, trails or seeds.
! The rebuild entry rebuilds a store with the retired stack policy, of another layout or
! grid, or no file, and leaves a valid store alone (a small project on disk); the store resolves
! against the working directory, not the cwd a copied project still records. The per-half sums of the
! power-spectrum bootstrap reproduce the seed spectra of the row-wise assembly (a floored spectrum and a
! repaired shell included) and carry the number of images summed for the average image (the floored one
! left out), and the global spectra summed while the records are written (a retained prefix with an
! invalid record included) equal a full reduction of the store.
module simple_sigma2_bootstrap_tester
use, intrinsic :: iso_fortran_env, only: int64, real32, real64
use simple_defs,              only: STDLEN, DTINY
use simple_string,            only: string
use simple_cmdline,           only: cmdline
use simple_fileio,            only: del_file
use simple_syslib,            only: simple_getcwd, simple_mkdir, simple_chdir, simple_rmdir
use simple_sp_project,        only: sp_project
use simple_ptcl_layout,       only: ptcl_layout_digest
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_init_header, &
    &sigma2_state_create_candidate, sigma2_state_write_particles, sigma2_state_write_global, &
    &sigma2_state_publish, sigma2_state_read_global, SIGMA2_PROV_PSPEC
use simple_sigma2_state,      only: sigma2_state_reduce
use simple_sigma2_pspec_sums, only: pspec_part_sums_fname, write_pspec_part_sums, read_pspec_part_sums, &
    &bootstrap_seed_spectra, sum_valid_records
use simple_sigma2_bootstrap,  only: prepare_pspec_cline, prepare_residual_sigma2_pass_cline, &
    &sigma2_estimate_available, ensure_sigma2_for_iteration
use simple_test_utils
use simple_sigma2_state_tester, only: patch_sigma2_header_word
implicit none
private
public :: run_all_sigma2_bootstrap_tests

character(len=*), parameter :: PROJ_FNAME  = 'tmp_sigma2_bootstrap_tester.simple'
character(len=*), parameter :: SIGMA_FNAME = 'tmp_sigma2_bootstrap_tester_sigma2.bin'
integer,          parameter :: NPTCLS = 12, BOX = 8, KTO = 4 ! shells 1..fdim(BOX)-1
real,             parameter :: SMPD   = 1.3

contains

    subroutine run_all_sigma2_bootstrap_tests()
        write(*,'(A)') '**** running all sigma2 bootstrap tests ****'
        call test_pspec_cline()
        call test_residual_pass_cline()
        call test_rebuild_entry()
        call test_pspec_sums()
    end subroutine run_all_sigma2_bootstrap_tests

    !> the bootstrap in per-half sums against the row-wise assembly it replaces, and the global spectra
    !! summed while writing against a full reduction
    subroutine test_pspec_sums()
        integer, parameter :: NP = 24, NSH = 6, NPART = 2, BLOCKROWS = 7, NPREFIX = 10
        integer, parameter :: SAMPLED(10) = [1, 2, 3, 4, 7, 8, 10, 13, 14, 20]
        character(len=*), parameter :: STORE = 'tmp_sigma2_bootstrap_tester_sums.bin'
        real(real64), allocatable :: part_sums(:,:)
        real(real32), allocatable :: block(:,:), state_globals(:,:)
        type(sigma2_state_header) :: header
        real(real64) :: rows(NSH,NP), row_sums(NSH,2), sums(NSH,2), expected(NSH,2), seeds(NSH,2), gsums(NSH,2)
        real(real32) :: records(NSH,NP), means(NSH,2)
        real    :: pspec_ave(NSH)
        logical :: active(NP)
        integer :: eo(NP), nactive(2), i, k, ipart, fromp, top, first_row, half, status, part_nimgs, nimgs
        integer(int64) :: counts(2), ninvalid
        character(len=STDLEN) :: msg
        write(*,'(A)') 'test_pspec_sums'
        do i = 1, NP
            eo(i)     = mod(i,2)
            active(i) = i /= 5 .and. i /= 11
        enddo
        nactive = [count(active .and. eo == 0), count(active .and. eo == 1)]
        ! the rows the workers wrote before: zero except the sampled particles' scaled spectra; row 7's
        ! spectrum is one that sanitising floored to 1 throughout
        rows = 0.d0
        do i = 1, size(SAMPLED)
            do k = 1, NSH
                rows(k,SAMPLED(i)) = real(real(1. + 0.37*real(k) + 0.011*real(SAMPLED(i))), real64)
            enddo
        enddo
        rows(:,7) = 1.d0
        ! the average image's power over two, large at shell 3: the even seed turns non-positive there
        do k = 1, NSH
            pspec_ave(k) = 0.05*real(k)
        enddo
        pspec_ave(3) = 5.
        ! the row-wise assembly: sums over the active rows by half, over the half's active count, minus the
        ! average power; the non-positive shell takes the previous shell's value
        row_sums = 0.d0
        do i = 1, NP
            if( active(i) ) row_sums(:,eo(i)+1) = row_sums(:,eo(i)+1) + rows(:,i)
        enddo
        do half = 1, 2
            expected(:,half) = row_sums(:,half)/real(nactive(half),real64) - real(pspec_ave,real64)
            do k = 2, NSH
                if( expected(k,half) <= DTINY ) expected(k,half) = expected(k-1,half)
            enddo
        enddo
        call assert_true(row_sums(3,1)/real(nactive(1),real64) - 5.d0 <= DTINY, 'the fixture repairs a shell')
        ! the workers' per-half sums of two partitions, through their files, with the number of images in their
        ! Fourier sums: the sampled ones but row 7, whose spectrum was floored
        sums  = 0.d0
        nimgs = 0
        do ipart = 1, NPART
            fromp = (ipart-1)*NP/NPART + 1
            top   = ipart*NP/NPART
            part_sums = reshape([(0.d0, i = 1, 2*NSH)], [NSH,2])
            do i = 1, size(SAMPLED)
                if( SAMPLED(i) < fromp .or. SAMPLED(i) > top ) cycle
                part_sums(:,eo(SAMPLED(i))+1) = part_sums(:,eo(SAMPLED(i))+1) + rows(:,SAMPLED(i))
            enddo
            part_nimgs = count(SAMPLED >= fromp .and. SAMPLED <= top .and. SAMPLED /= 7)
            call write_pspec_part_sums(pspec_part_sums_fname(ipart), fromp, top, part_nimgs, part_sums)
            deallocate(part_sums)
            part_nimgs = -1
            call read_pspec_part_sums(pspec_part_sums_fname(ipart), NSH, fromp, top, part_nimgs, part_sums, status, msg)
            call assert_int(0, status, 'a partition reads back its sums: '//trim(msg))
            call assert_true(fromp == (ipart-1)*NP/NPART + 1 .and. top == ipart*NP/NPART, 'with its particle range')
            call assert_int(count(SAMPLED >= fromp .and. SAMPLED <= top .and. SAMPLED /= 7), part_nimgs, &
                &'and the number of images its Fourier sum holds')
            if( ipart == 1 ) call check_refused_sums(pspec_part_sums_fname(ipart), NSH)
            sums  = sums + part_sums
            nimgs = nimgs + part_nimgs
            call del_file(pspec_part_sums_fname(ipart))
        enddo
        call assert_int(size(SAMPLED)-1, nimgs, 'the average image divides by the images summed, not the active count')
        call bootstrap_seed_spectra(sums, nactive, pspec_ave, seeds)
        call assert_true(all(abs(seeds - expected) <= 1.d-6*abs(expected)), &
            &'the per-half sums reproduce the row-wise seed spectra')
        call assert_true(all(seeds > DTINY), 'the seed spectra are positive')
        ! the records: a retained prefix (with one invalid active record), then the seeds by half; inactive zero
        records = 0.
        do i = 1, NPREFIX
            do k = 1, NSH
                records(k,i) = 0.5 + 0.1*real(k) + 0.03*real(i)
            enddo
        enddo
        records(2,4) = -1.
        do i = NPREFIX+1, NP
            if( active(i) ) records(:,i) = real(seeds(:,eo(i)+1), real32)
        enddo
        gsums    = 0.d0
        counts   = 0_int64
        ninvalid = 0_int64
        do first_row = 1, NP, BLOCKROWS
            block = records(:,first_row:min(NP,first_row+BLOCKROWS-1))
            call sum_valid_records(block, first_row, active, eo, gsums, counts, ninvalid)
        enddo
        call assert_true(ninvalid == 1_int64, 'the invalid prefix record is skipped')
        do half = 1, 2
            means(:,half) = real(gsums(:,half)/real(counts(half),real64), real32)
        enddo
        ! a full reduction of the same records
        call del_file(STORE)
        call sigma2_state_init_header(header, 1, NSH, NP, 2*NSH, SMPD, 1_int64, 4242_int64, SIGMA2_PROV_PSPEC)
        call sigma2_state_create_candidate(STORE, header, status, msg)
        call assert_int(0, status, 'create the reduction store: '//trim(msg))
        call sigma2_state_write_particles(STORE, 1, records, status, msg)
        call assert_int(0, status, 'write its records: '//trim(msg))
        call sigma2_state_reduce(STORE, active, eo, status, msg)
        call assert_int(0, status, 'reduce it: '//trim(msg))
        call sigma2_state_read_global(STORE, state_globals, status, msg)
        call assert_int(0, status, 'read its global spectra: '//trim(msg))
        if( status == 0 ) call assert_true(all(abs(state_globals - means) <= 1.e-6*abs(means)), &
            &'the global spectra summed while writing equal a full reduction')
        call del_file(STORE)
    end subroutine test_pspec_sums

    !> a partition's sums file is refused, before anything is allocated, when it was summed over another
    !! number of shells, when its image count lies outside its particle range, when its header claims a huge
    !! shell count, and when it is longer than its header says
    subroutine check_refused_sums( fname, nsh )
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: nsh
        integer, parameter :: NSHELL_POS = 25  !< after the magic (16 bytes), fromp and top
        integer, parameter :: NIMGS_POS  = 29  !< after nshell
        real(real64), allocatable :: sums(:,:)
        integer(int64) :: fsize
        integer :: funit, fromp, top, nimgs, nimgs_on_disk, status
        character(len=STDLEN) :: msg
        call read_pspec_part_sums(fname, nsh+1, fromp, top, nimgs, sums, status, msg)
        call assert_true(status /= 0 .and. .not. allocated(sums), 'sums over another number of shells are refused')
        open(newunit=funit, file=fname, access='stream', form='unformatted', status='old', action='readwrite')
        read(funit, pos=NIMGS_POS) nimgs_on_disk
        write(funit, pos=NIMGS_POS) -1
        close(funit)
        call read_pspec_part_sums(fname, nsh, fromp, top, nimgs, sums, status, msg)
        call assert_true(status /= 0 .and. .not. allocated(sums), 'a negative image count is refused')
        open(newunit=funit, file=fname, access='stream', form='unformatted', status='old', action='readwrite')
        write(funit, pos=NIMGS_POS) nimgs_on_disk
        write(funit, pos=NSHELL_POS) huge(nsh)
        close(funit)
        call read_pspec_part_sums(fname, nsh, fromp, top, nimgs, sums, status, msg)
        call assert_true(status /= 0 .and. .not. allocated(sums), 'a header claiming a huge shell count is refused')
        open(newunit=funit, file=fname, access='stream', form='unformatted', status='old', action='readwrite')
        write(funit, pos=NSHELL_POS) nsh
        inquire(unit=funit, size=fsize)
        write(funit, pos=fsize+1_int64) 0.d0
        close(funit)
        call read_pspec_part_sums(fname, nsh, fromp, top, nimgs, sums, status, msg)
        call assert_true(status /= 0 .and. .not. allocated(sums), 'a file longer than its header says is refused')
    end subroutine check_refused_sums

    !> a stage reconstruction command line with every kind of key it carries
    subroutine make_template( cl )
        type(cmdline), intent(inout) :: cl
        call cl%kill
        call cl%set('prg',         'reconstruct3D')
        call cl%set('projfile',    'run.simple')
        call cl%set('oritype',     'ptcl3D')
        call cl%set('mskdiam',     180.)
        call cl%set('nthr',        8)
        call cl%set('nparts',      4)
        call cl%set('qsys_name',   'slurm')
        call cl%set('walltime',    3600)
        call cl%set('box_crop',    128)
        call cl%set('lp',          8.)
        call cl%set('update_frac', 0.5)
        call cl%set('trail_rec',   'yes')
        call cl%set('trail_seed',  'yes')
        call cl%set('frozen_rec',  '/abs/run/frozen_context.txt')
        call cl%set('rec_backend', 'pcg')
        call cl%set('vol1',        'recvol_state01.mrc')
        call cl%set('which_iter',  7)
    end subroutine make_template

    subroutine test_pspec_cline()
        type(cmdline) :: tmpl, cl
        type(string)  :: sval
        write(*,'(A)') 'test_pspec_cline'
        call make_template(tmpl)
        call prepare_pspec_cline(tmpl, string('/abs/run/run.simple'), 'ptcl2D', 0, cl)
        sval = cl%get_carg('prg')
        call assert_string_eq('calc_pspec', sval, 'the calc_pspec program')
        sval = cl%get_carg('projfile')
        call assert_string_eq('/abs/run/run.simple', sval, 'the given project file')
        sval = cl%get_carg('objfun')
        call assert_string_eq('euclid', sval, 'the euclid objective')
        sval = cl%get_carg('mkdir')
        call assert_string_eq('no', sval, 'no run directory')
        call assert_int(1, cl%get_iarg('which_iter'), 'the iteration is at least 1, not the template''s')
        sval = cl%get_carg('oritype')
        call assert_string_eq('ptcl2D', sval, 'the particle segment of the consumer, not the template''s')
        sval = cl%get_carg('qsys_name')
        call assert_string_eq('slurm', sval, 'the queue of the template')
        call assert_real(180., cl%get_rarg('mskdiam'), 0., 'the mask of the template')
        call assert_int(8,    cl%get_iarg('nthr'),     'the threads of the template')
        call assert_int(4,    cl%get_iarg('nparts'),   'the partitions of the template')
        call assert_int(3600, cl%get_iarg('walltime'), 'the walltime of the template')
        call assert_int(11,   cl%get_argcnt(),         'nothing else: six keys set, five copied')
        call assert_false(cl%defined('update_frac') .or. cl%defined('trail_rec') .or. cl%defined('trail_seed') &
            &.or. cl%defined('frozen_rec') .or. cl%defined('rec_backend') .or. cl%defined('vol1') &
            &.or. cl%defined('box_crop') .or. cl%defined('lp'), 'no sampling, trailing, frozen, backend or map key')
        call tmpl%kill
        call cl%kill
        call sval%kill
    end subroutine test_pspec_cline

    subroutine test_residual_pass_cline()
        type(cmdline) :: tmpl, cl
        type(string)  :: sval
        write(*,'(A)') 'test_residual_pass_cline'
        call make_template(tmpl)
        call prepare_residual_sigma2_pass_cline(tmpl, 3, 1, [string('bootstrap_state01.mrc')], cl)
        sval = cl%get_carg('prg')
        call assert_string_eq('refine3D', sval, 'a refine3D pass')
        sval = cl%get_carg('refine')
        call assert_string_eq('sigma', sval, 'in residual sigma2 mode')
        call assert_int(3, cl%get_iarg('which_iter'), 'at the given iteration')
        sval = cl%get_carg('vol1')
        call assert_string_eq('bootstrap_state01.mrc', sval, 'against the given state map')
        call assert_real(8.,   cl%get_rarg('lp'),       0., 'the scoring limit of the template')
        call assert_int(128,   cl%get_iarg('box_crop'),     'the crop of the template')
        call assert_real(180., cl%get_rarg('mskdiam'),  0., 'the mask of the template')
        call assert_false(cl%defined('update_frac') .or. cl%defined('trail_rec') .or. cl%defined('trail_seed') &
            &.or. cl%defined('frozen_rec'), 'no sampling, trailing or frozen key')
        call tmpl%kill
        call cl%kill
        call sval%kill
    end subroutine test_residual_pass_cline

    !> NPTCLS particles in two stacks, written to PROJ_FNAME with SIGMA_FNAME registered
    subroutine make_project( spproj )
        type(sp_project), intent(inout) :: spproj
        integer :: i, istk
        call spproj%kill
        ! as a project read from its file: named after the file, in the file's directory
        call spproj%update_projinfo(string(PROJ_FNAME))
        call spproj%os_stk%new(2, is_ptcl=.false.)
        do istk = 1, 2
            call spproj%os_stk%set(istk, 'stk',        'stack_'//char(48+istk)//'.mrcs')
            call spproj%os_stk%set(istk, 'fromp',      (istk-1)*NPTCLS/2 + 1)
            call spproj%os_stk%set(istk, 'top',        istk*NPTCLS/2)
            call spproj%os_stk%set(istk, 'nptcls_stk', NPTCLS/2)
            call spproj%os_stk%set(istk, 'box',        BOX)
            call spproj%os_stk%set(istk, 'smpd',       SMPD)
        enddo
        call spproj%os_ptcl3D%new(NPTCLS, is_ptcl=.true.)
        do i = 1, NPTCLS
            istk = merge(1, 2, i <= NPTCLS/2)
            call spproj%os_ptcl3D%set(i, 'stkind', istk)
            call spproj%os_ptcl3D%set(i, 'indstk', i - (istk-1)*NPTCLS/2)
            call spproj%os_ptcl3D%set(i, 'eo',     mod(i,2))
            call spproj%os_ptcl3D%set_state(i, 1)
        enddo
        call spproj%projinfo%set(1, 'sigma2_state', SIGMA_FNAME)
        call spproj%write(string(PROJ_FNAME))
    end subroutine make_project

    !> a committed store for spproj; a layout digest off by one stands for a store of another
    !! particle layout
    subroutine write_store( spproj, l_other_layout )
        type(sp_project), intent(inout) :: spproj
        logical,          intent(in)    :: l_other_layout
        type(sigma2_state_header) :: header
        type(string)   :: candidate
        real(real32)   :: spectra(KTO,NPTCLS), globals(KTO,2)
        integer(int64) :: digest
        integer        :: status
        character(len=STDLEN) :: msg
        digest  = ptcl_layout_digest(spproj, spproj%os_ptcl3D)
        if( l_other_layout ) digest = digest + 1_int64
        candidate = SIGMA_FNAME//'.next'
        call del_file(candidate)
        call del_file(SIGMA_FNAME)
        call sigma2_state_init_header(header, 1, KTO, NPTCLS, BOX, SMPD, 1_int64, digest, SIGMA2_PROV_PSPEC)
        call sigma2_state_create_candidate(candidate%to_char(), header, status, msg)
        call assert_int(0, status, 'create the store fixture: '//trim(msg))
        spectra = 1.0_real32
        call sigma2_state_write_particles(candidate%to_char(), 1, spectra, status, msg)
        call assert_int(0, status, 'write the store fixture particles: '//trim(msg))
        globals = 1.0_real32
        call sigma2_state_write_global(candidate%to_char(), globals, status, msg)
        call assert_int(0, status, 'write the store fixture spectra: '//trim(msg))
        call sigma2_state_publish(candidate%to_char(), SIGMA_FNAME, status, msg)
        call assert_int(0, status, 'publish the store fixture: '//trim(msg))
        call candidate%kill
    end subroutine write_store

    subroutine test_rebuild_entry()
        type(sp_project) :: spproj
        type(cmdline)    :: tmpl
        logical          :: l_bootstrapped
        write(*,'(A)') 'test_rebuild_entry'
        call make_project(spproj)
        call write_store(spproj, .false.)
        call assert_true(sigma2_estimate_available(string(PROJ_FNAME), 'ptcl3D', BOX, SMPD), &
            &'a valid global store is consumable')
        ! the rebuild entry leaves it alone: no calc_pspec runs
        call tmpl%set('nthr', 1)
        call ensure_sigma2_for_iteration(tmpl, string(PROJ_FNAME), 'ptcl3D', 1, BOX, SMPD, 'TESTER', l_bootstrapped)
        call assert_false(l_bootstrapped, 'a valid global store is not rebuilt')
        ! a store that is not consumable is exactly what makes ensure_sigma2_for_iteration run calc_pspec,
        ! which needs particle images; the cases below assert that condition
        call assert_false(sigma2_estimate_available(string(PROJ_FNAME), 'ptcl3D', 2*BOX, SMPD), &
            &'a store of another grid triggers the rebuild')
        ! the grouping word (header word 9) of the retired per-stack policy
        call patch_sigma2_header_word(SIGMA_FNAME, 9, 2_int64)
        call assert_false(sigma2_estimate_available(string(PROJ_FNAME), 'ptcl3D', BOX, SMPD), &
            &'a store with the stack policy triggers the rebuild')
        call write_store(spproj, .true.)
        call assert_false(sigma2_estimate_available(string(PROJ_FNAME), 'ptcl3D', BOX, SMPD), &
            &'a store of another particle layout triggers the rebuild')
        call del_file(SIGMA_FNAME)
        call assert_false(sigma2_estimate_available(string(PROJ_FNAME), 'ptcl3D', BOX, SMPD), &
            &'a missing store triggers the rebuild')
        call test_new_exec_dir(spproj)
        call del_file(PROJ_FNAME)
        call tmpl%kill
        call spproj%kill
    end subroutine test_rebuild_entry

    !> a copy of the project in a new execution directory still records the previous run's cwd,
    !! where a valid store lies; the consumers resolve the store against the project file's
    !! directory, so the entry must find none there and rebuild
    subroutine test_new_exec_dir( spproj )
        type(sp_project), intent(inout) :: spproj
        character(len=*), parameter :: EXEC_DIR = 'tmp_sigma2_bootstrap_tester_exec'
        type(string) :: prev_dir
        call make_project(spproj)
        call write_store(spproj, .false.)
        call simple_getcwd(prev_dir)
        call spproj%projinfo%set(1, 'cwd', prev_dir)
        call simple_mkdir(EXEC_DIR)
        call simple_chdir(EXEC_DIR)
        call spproj%write(string(PROJ_FNAME))
        call assert_false(sigma2_estimate_available(string(PROJ_FNAME), 'ptcl3D', BOX, SMPD), &
            &'a new execution directory does not consume the store of the cwd its project records')
        call del_file(PROJ_FNAME)
        call simple_chdir(prev_dir)
        call simple_rmdir(EXEC_DIR)
        call assert_true(sigma2_estimate_available(string(PROJ_FNAME), 'ptcl3D', BOX, SMPD), &
            &'the run that wrote the store still consumes it')
        call del_file(SIGMA_FNAME)
        call prev_dir%kill
    end subroutine test_new_exec_dir

end module simple_sigma2_bootstrap_tester
