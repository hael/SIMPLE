!@descr: unit tests of the sigma2 bootstrap command lines (simple_sigma2_bootstrap)
! prepare_pspec_cline builds calc_pspec's command line fresh, with the few keys
! it takes from its template; prepare_residual_sigma2_pass_cline keeps the
! template's scoring keys and drops what samples, assembles, trails or seeds.
! In-memory command lines only.
module simple_sigma2_bootstrap_tester
use simple_string,           only: string
use simple_cmdline,          only: cmdline
use simple_sigma2_bootstrap, only: prepare_pspec_cline, prepare_residual_sigma2_pass_cline
use simple_test_utils
implicit none
private
public :: run_all_sigma2_bootstrap_tests

contains

    subroutine run_all_sigma2_bootstrap_tests()
        write(*,'(A)') '**** running all sigma2 bootstrap tests ****'
        call test_pspec_cline()
        call test_residual_pass_cline()
    end subroutine run_all_sigma2_bootstrap_tests

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
        call prepare_pspec_cline(tmpl, string('/abs/run/run.simple'), 0, cl)
        sval = cl%get_carg('prg')
        call assert_string_eq('calc_pspec', sval, 'the calc_pspec program')
        sval = cl%get_carg('projfile')
        call assert_string_eq('/abs/run/run.simple', sval, 'the given project file')
        sval = cl%get_carg('objfun')
        call assert_string_eq('euclid', sval, 'the euclid objective')
        sval = cl%get_carg('sigma_est')
        call assert_string_eq('global', sval, 'global sigma2')
        sval = cl%get_carg('mkdir')
        call assert_string_eq('no', sval, 'no run directory')
        call assert_int(1, cl%get_iarg('which_iter'), 'the iteration is at least 1, not the template''s')
        sval = cl%get_carg('oritype')
        call assert_string_eq('ptcl3D', sval, 'the particle segment of the template')
        sval = cl%get_carg('qsys_name')
        call assert_string_eq('slurm', sval, 'the queue of the template')
        call assert_real(180., cl%get_rarg('mskdiam'), 0., 'the mask of the template')
        call assert_int(8,    cl%get_iarg('nthr'),     'the threads of the template')
        call assert_int(4,    cl%get_iarg('nparts'),   'the partitions of the template')
        call assert_int(3600, cl%get_iarg('walltime'), 'the walltime of the template')
        call assert_int(12,   cl%get_argcnt(),         'nothing else: six keys set, six copied')
        call assert_false(cl%defined('update_frac') .or. cl%defined('trail_rec') .or. cl%defined('trail_seed') &
            &.or. cl%defined('frozen_rec') .or. cl%defined('rec_backend') .or. cl%defined('vol1') &
            &.or. cl%defined('box_crop') .or. cl%defined('lp'), 'no sampling, trailing, frozen, backend or map key')
        ! a consumer that keeps per-stack sigma2 groups asks for them
        call prepare_pspec_cline(tmpl, string('/abs/run/run.simple'), 0, cl, sigma_est='group')
        sval = cl%get_carg('sigma_est')
        call assert_string_eq('group', sval, 'per-stack sigma2 when asked for')
        call assert_int(12,   cl%get_argcnt(),         'the same twelve keys')
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

end module simple_sigma2_bootstrap_tester
