!@descr: for all utils tests
module simple_commanders_test_utils
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_ansi_colors
  contains
    procedure :: execute      => exec_test_ansi_colors
end type commander_test_ansi_colors

type, extends(commander_base) :: commander_test_binoris_test
  contains
    procedure :: execute      => exec_test_binoris_test
end type commander_test_binoris_test

type, extends(commander_base) :: commander_test_binoris_io_test
  contains
    procedure :: execute      => exec_test_binoris_io_test
end type commander_test_binoris_io_test

type, extends(commander_base) :: commander_test_cif2mrc
  contains
    procedure :: execute      => exec_test_cif2mrc
end type commander_test_cif2mrc

type, extends(commander_base) :: commander_test_cif2pdb
  contains
    procedure :: execute      => exec_test_cif2pdb
end type commander_test_cif2pdb

type, extends(commander_base) :: commander_test_cmdline
  contains
    procedure :: execute      => exec_test_cmdline
end type commander_test_cmdline

type, extends(commander_base) :: commander_test_install
  contains
    procedure :: execute      => exec_test_install
end type commander_test_install

type, extends(commander_base) :: commander_test_nice
  contains
    procedure :: execute      => exec_test_nice
end type commander_test_nice

type, extends(commander_base) :: commander_test_pdb2mrc 
  contains
    procedure :: execute      => exec_test_pdb2mrc
end type commander_test_pdb2mrc

type, extends(commander_base) :: commander_test_serialize
  contains
    procedure :: execute      => exec_test_serialize
end type commander_test_serialize

type, extends(commander_base) :: commander_test_stringmatch
  contains
    procedure :: execute      => exec_test_stringmatch
end type commander_test_stringmatch

contains

subroutine exec_test_ansi_colors( self, cline )
    use simple_ansi_ctrls
    use simple_defs_fname, only: NEWLINE
    class(commander_test_ansi_colors), intent(inout) :: self
    class(cmdline),                    intent(inout) :: cline
    print '(a)', &
        format_str('Red',     C_RED)     // NEWLINE // &
        format_str('Green',   C_GREEN)   // NEWLINE // &
        format_str('Yellow',  C_YELLOW)  // NEWLINE // &
        format_str('Blue',    C_BLUE)    // NEWLINE // &
        format_str('Magenta', C_MAGENTA) // NEWLINE // &
        format_str('Cyan',    C_CYAN)    // NEWLINE // &
        format_str('White',   C_WHITE)
    call simple_end('**** SIMPLE_TEST_ANSI_COLORS_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_ansi_colors

subroutine exec_test_binoris_test( self, cline )
    class(commander_test_binoris_test),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_BINORIS_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_binoris_test

subroutine exec_test_binoris_io_test( self, cline )
    class(commander_test_binoris_io_test),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_BINORIS_IO_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_binoris_io_test

subroutine exec_test_cif2mrc( self, cline )
    use simple_image,         only : image
    use simple_atoms,         only : atoms
    class(commander_test_cif2mrc), intent(inout) :: self
    class(cmdline),                intent(inout) :: cline
    type(string)    :: cif_file, cmd
    type(atoms)     :: molecule
    real, parameter :: smpd = 1.3
    integer         :: rc
    cif_file = '6VXX.cif'
    print *, 'Downloading 6VXX cif file...'
    cmd = 'curl -L -o '//cif_file%to_char()//' https://files.rcsb.org/download/6VXX.cif'
    call execute_command_line(cmd%to_char(), exitstat=rc)
    call molecule%cif2mrc(cif_file, smpd)
    call simple_end('**** SIMPLE_TEST_CIF2MRC_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_cif2mrc

subroutine exec_test_cif2pdb( self, cline )
    use simple_atoms,         only : atoms
    class(commander_test_cif2pdb), intent(inout) :: self
    class(cmdline),                intent(inout) :: cline
    type(atoms)  :: molecule
    type(string) :: cif_file, pdb_file, cmd
    integer      :: rc
    cif_file = '6VXX.cif'
    pdb_file = '6VXX.pdb'
    print *, 'Downloading 6VXX cif file...'
    cmd = 'curl -L -o '//cif_file%to_char()//' https://files.rcsb.org/download/6VXX.cif'
    call execute_command_line(cmd%to_char(), exitstat=rc)
    call molecule%cif2pdb(cif_file, pdb_file)
    call simple_end('**** SIMPLE_TEST_CIF2PDB_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_cif2pdb

subroutine exec_test_cmdline( self, cline )
    use simple_cmdline
    class(commander_test_cmdline),    intent(inout) :: self
    class(cmdline),                   intent(inout) :: cline
    type(chash)           :: job_descr
    character(len=STDLEN) :: xarg, line
    type(string)          :: fname, str_prg
    logical               :: test_passed
    integer               :: cmdstat, cmdlen, pos
    test_passed  = .true.
    fname        = 'file_command.txt'
    xarg         = "prg=symmetrize_map"
    cmdlen       = len(trim(xarg))
    line         = "projname=system_name smpd=1.3 cs=2.7 kv=300 fraca=0.1 total_dose=53 dir_movies=/usr/local/data/movies &
                  & gainref=gainref.mrc nparts=4 nthr=16 moldiam_max=200"
    cmdstat      = 1
    pos  = index(xarg, '=')
    call cmdline_err(cmdstat, cmdlen, xarg, pos)
    call cline%read(line)
    print *, '>>> CLINE'
    call cline%printline()
    call cline%writeline(fname)
    print *, '>>> CLINE'
    call cline%printline()
    call cline%checkvar('projname',        1)
    call cline%checkvar('smpd',            2)
    call cline%checkvar('cs',              3)
    call cline%checkvar('kv',              4)
    call cline%checkvar('fraca',           5)
    call cline%checkvar('total_dose',      6)
    call cline%checkvar('dir_movies',      7)
    call cline%checkvar('gainref',         8)
    call cline%checkvar('nparts',          9)
    call cline%checkvar('nthr',           10)
    call cline%checkvar('moldiam_max',    11)
    if( .not. cline%defined('projname')   ) test_passed=.false.
    if( .not. cline%defined('smpd')       ) test_passed=.false.
    if( .not. cline%defined('dir_movies') ) test_passed=.false.
    if( .not. cline%defined('gainref')    ) test_passed=.false.
    if( .not. cline%defined('moldiam_max')) test_passed=.false.
    call cline%check()
    call cline%gen_job_descr(job_descr)
    call job_descr%set('prg',      'scale')
    call job_descr%set('autoscale','no')
    print *,'>>> JOB_DESCR'
    call job_descr%print_key_val_pairs(logfhandle) 
    call cline%set('prg', 'list')
    str_prg = cline%get_carg('prg')
    print *,'>>> PROGRAM ', str_prg%to_char()
    if( .not. (cline%get_carg('prg').eq.'list')) test_passed=.false.
    call cline%delete('kv')
    print *,'>>> CS ',cline%get_rarg('cs')
    print *,'>>> NPARTS ',cline%get_iarg('nparts')
    if( test_passed )then
        print *, '>>> TEST PASSED'
    else
        THROW_HARD('>>> TEST FAILED')
    endif
    call cline%kill()
    call job_descr%kill()
    call simple_end('**** SIMPLE_TEST_CMDLINE_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_cmdline

subroutine exec_test_install( self, cline )
    use simple_testfuns ! use all in there
    use simple_image,  only: image
    class(commander_test_install),    intent(inout) :: self
    class(cmdline),                   intent(inout) :: cline
    type( image )         :: cube, img
    real                  :: smpd
    integer               :: box, nspace, msk
    character(len=8)      :: datestr
    character(len=STDLEN) :: folder
    character(len=300)    :: command
    call seed_rnd
    call date_and_time(date=datestr)
    folder = trim('./SIMPLE_TEST_INSTALL_'//datestr)
    call simple_mkdir(trim( folder ))
    call simple_chdir( folder )
    ! dummy data
    box    = 96
    smpd   = 2.
    nspace = 64
    msk    = nint(real(box)/3.)
    ! volume
    call img%new( [box,box,box], smpd )
    call img%square( nint(real(box)/12.) )
    call cube%new( [box,box,box], smpd )
    call cube%square( nint(real(box)/16.) )
    call cube%shift([16.,16.,16.])
    call img%add( cube )
    call cube%new( [box,box,box], smpd )
    call cube%square( nint(real(box)/10.) )
    call cube%shift([4.,-16.,0.])
    call img%add( cube )
    call cube%kill
    call img%write(string('cubes.mrc'))
    call img%kill
    write(logfhandle,*)'>>> WROTE TEST VOLUME cubes.mrc'
    ! test units
    command = 'simple_test_units'
    call exec_cmdline( trim(command) )
    ! test search
    ! command = 'simple_test_srch vol1=cubes.mrc msk='//int2str(msk)//&
    !     & ' smpd='//real2str(smpd)
    ! call exec_cmdline( trim(command) )
    ! end
    call simple_chdir(PATH_PARENT)
    call simple_end('**** SIMPLE_TEST_INSTALL_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_install

subroutine exec_test_nice( self, cline )
    use simple_nice
    class(commander_test_nice),    intent(inout) :: self
    class(cmdline),                intent(inout) :: cline
    type(simple_nice_comm) :: nice_comm
    call nice_comm%init(1, "testserver")
    call nice_comm%cycle()
    call sleep(5)
    nice_comm%view_micrographs%active         = .true.
    nice_comm%view_micrographs%thumbnail%path = "/tmp/cls2D_thumbnail.jpeg"
    nice_comm%view_micrographs%thumbnail%id   = 10 ! should be random
    call nice_comm%cycle()
    call sleep(5)
    nice_comm%view_micrographs%active         = .true.
    nice_comm%view_micrographs%thumbnail%path = "/tmp/cls2D_thumnail.jpeg"
    nice_comm%view_micrographs%thumbnail%id   = 11! should be random
    call nice_comm%cycle()
    call sleep(10)
    call nice_comm%terminate
    call simple_end('**** SIMPLE_TEST_NICE_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_nice

subroutine exec_test_pdb2mrc( self, cline )
    use simple_image,         only : image
    use simple_atoms,         only : atoms
    use simple_molecule_data, only : molecule_data, betagal_1jyx, sars_cov2_spkgp_6vxx
    use simple_imghead,       only : find_ldim_nptcls, find_img_smpd
    class(commander_test_pdb2mrc), intent(inout) :: self
    class(cmdline),                intent(inout) :: cline
    type(string)        :: pdb_file, vol_file, default_vol
    type(atoms)         :: molecule
    type(molecule_data) :: mol
    real, parameter     :: smpd = 1.3
    integer             :: ldim(3), nptcls
    real                :: smpd_out
    logical             :: all_ok
    all_ok = .true.
    ! ---- Test 1: mol-only input (6VXX) ----
    write(logfhandle,'(a)') '>>> Test 1: pdb2mrc with 6VXX mol data (default filenames)'
    mol = sars_cov2_spkgp_6vxx()
    write(logfhandle,'(a,i8)') '    molecule atom count: ', mol%n
    if( mol%n < 1 )then
        write(logfhandle,'(a)') '    FAIL: 6VXX molecule data has no atoms'
        all_ok = .false.
    else
        write(logfhandle,'(a)') '    PASS: 6VXX molecule data loaded'
    endif
    call molecule%pdb2mrc(smpd=smpd, mol=mol)
    default_vol = 'molecule.mrc'
    if( .not. file_exists(default_vol) )then
        write(logfhandle,'(a)') '    FAIL: '//default_vol%to_char()//' not created'
        all_ok = .false.
    else
        call find_ldim_nptcls(default_vol, ldim, nptcls)
        smpd_out = find_img_smpd(default_vol)
        write(logfhandle,'(a,i4,a,i4,a,i4,a,f6.2)') '    volume dims = [', &
            ldim(1),',',ldim(2),',',ldim(3),' ], smpd = ', smpd_out
        if( ldim(1) < 1 .or. ldim(2) < 1 .or. ldim(3) < 1 )then
            write(logfhandle,'(a)') '    FAIL: volume has invalid dimensions'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: volume dimensions valid'
        endif
        if( abs(smpd_out - smpd) > 0.01 )then
            write(logfhandle,'(a,f6.2,a,f6.2)') '    FAIL: smpd mismatch, expected ', smpd, ' got ', smpd_out
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: sampling distance matches'
        endif
    endif
    ! ---- Test 2: mol-only input (1JYX) ----
    write(logfhandle,'(a)') '>>> Test 2: pdb2mrc with 1JYX mol data (default filenames)'
    mol = betagal_1jyx()
    write(logfhandle,'(a,i8)') '    molecule atom count: ', mol%n
    if( mol%n < 1 )then
        write(logfhandle,'(a)') '    FAIL: 1JYX molecule data has no atoms'
        all_ok = .false.
    else
        write(logfhandle,'(a)') '    PASS: 1JYX molecule data loaded'
    endif
    call molecule%pdb2mrc(smpd=smpd, mol=mol)
    if( .not. file_exists(default_vol) )then
        write(logfhandle,'(a)') '    FAIL: '//default_vol%to_char()//' not created'
        all_ok = .false.
    else
        call find_ldim_nptcls(default_vol, ldim, nptcls)
        smpd_out = find_img_smpd(default_vol)
        write(logfhandle,'(a,i4,a,i4,a,i4,a,f6.2)') '    volume dims = [', &
            ldim(1),',',ldim(2),',',ldim(3),' ], smpd = ', smpd_out
        if( ldim(1) < 1 .or. ldim(2) < 1 .or. ldim(3) < 1 )then
            write(logfhandle,'(a)') '    FAIL: volume has invalid dimensions'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: volume dimensions valid'
        endif
        if( abs(smpd_out - smpd) > 0.01 )then
            write(logfhandle,'(a,f6.2,a,f6.2)') '    FAIL: smpd mismatch, expected ', smpd, ' got ', smpd_out
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: sampling distance matches'
        endif
    endif
    ! ---- Test 3: explicit filenames (6VXX) ----
    write(logfhandle,'(a)') '>>> Test 3: pdb2mrc with 6VXX mol data + explicit pdb/vol filenames'
    pdb_file = '6VXX.pdb'
    vol_file = '6VXX.mrc'
    mol = sars_cov2_spkgp_6vxx()
    call molecule%pdb2mrc(pdbfile=pdb_file, volfile=vol_file, smpd=smpd, mol=mol)
    if( .not. file_exists(vol_file) )then
        write(logfhandle,'(a)') '    FAIL: '//vol_file%to_char()//' not created'
        all_ok = .false.
    else
        call find_ldim_nptcls(vol_file, ldim, nptcls)
        smpd_out = find_img_smpd(vol_file)
        write(logfhandle,'(a,i4,a,i4,a,i4,a,f6.2)') '    volume dims = [', &
            ldim(1),',',ldim(2),',',ldim(3),' ], smpd = ', smpd_out
        if( ldim(1) < 1 .or. ldim(2) < 1 .or. ldim(3) < 1 )then
            write(logfhandle,'(a)') '    FAIL: volume has invalid dimensions'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: 6VXX volume dimensions valid'
        endif
        if( abs(smpd_out - smpd) > 0.01 )then
            write(logfhandle,'(a,f6.2,a,f6.2)') '    FAIL: smpd mismatch, expected ', smpd, ' got ', smpd_out
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: 6VXX sampling distance matches'
        endif
    endif
    ! ---- Test 4: explicit filenames (1JYX) ----
    write(logfhandle,'(a)') '>>> Test 4: pdb2mrc with 1JYX mol data + explicit pdb/vol filenames'
    pdb_file = '1JYX.pdb'
    vol_file = '1JYX.mrc'
    mol      = betagal_1jyx()
    call molecule%pdb2mrc(pdbfile=pdb_file, volfile=vol_file, smpd=smpd, mol=mol)
    if( .not. file_exists(vol_file) )then
        write(logfhandle,'(a)') '    FAIL: '//vol_file%to_char()//' not created'
        all_ok = .false.
    else
        call find_ldim_nptcls(vol_file, ldim, nptcls)
        smpd_out = find_img_smpd(vol_file)
        write(logfhandle,'(a,i4,a,i4,a,i4,a,f6.2)') '    volume dims = [', &
            ldim(1),',',ldim(2),',',ldim(3),' ], smpd = ', smpd_out
        if( ldim(1) < 1 .or. ldim(2) < 1 .or. ldim(3) < 1 )then
            write(logfhandle,'(a)') '    FAIL: volume has invalid dimensions'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: 1JYX volume dimensions valid'
        endif
        if( abs(smpd_out - smpd) > 0.01 )then
            write(logfhandle,'(a,f6.2,a,f6.2)') '    FAIL: smpd mismatch, expected ', smpd, ' got ', smpd_out
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: 1JYX sampling distance matches'
        endif
    endif
    ! ---- Final verdict ----
    if( all_ok )then
        call simple_end('**** SIMPLE_TEST_PDB2MRC NORMAL STOP ****')
    else
        THROW_HARD('TEST_PDB2MRC FAILED')
    endif
end subroutine exec_test_pdb2mrc

subroutine exec_test_serialize( self, cline )
    use simple_image
    use simple_stackops
    use simple_ppca
    class(commander_test_serialize),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    type(image)          :: img, img_msk, img_rev
    real,    allocatable :: pcavec(:)
    logical, allocatable :: l_mask(:,:,:)
    integer, parameter   :: box = 256
    call img%new([box,box,1], 1.0)
    call img_rev%new([box,box,1], 1.0)
    call img%square(60)
    call img_msk%disc([box,box,1], 1.0, 80., l_mask)
    call img_msk%vis
    call img%vis
    pcavec = img%serialize(l_mask)
    call img_rev%unserialize(pcavec, l_mask)
    call img_rev%vis
    call simple_end('**** SIMPLE_TEST_SERIALIZE_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_serialize

subroutine exec_test_stringmatch( self, cline )
    class(commander_test_stringmatch),  intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    character(len=*), parameter :: projname = 'state_2'
    character(len=*), parameter :: stkname  = '2_state.mrc'
    character(len=*), parameter :: inds_str = ' 1,3,  5,7,15  '
    integer,        allocatable :: inds(:)
    integer :: i 
    inds = list_of_ints2arr(inds_str)
    print *, 'found # integer numbers: ', size(inds)
    do i = 1, size(inds)
        print *, inds(i)
    end do
    call simple_end('**** SIMPLE_TEST_STRINGMATCH_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_stringmatch

end module simple_commanders_test_utils
