!@descr: reference, sigma and chain provenance identity and mutation checks
module simple_refine3D_provenance_tester
use simple_string,     only: string
use simple_sp_project, only: sp_project
use simple_image,      only: image
use simple_test_utils, only: assert_true, assert_false, assert_int, enter_fixture, leave_fixture, tests_failed
implicit none
private
public :: run_all_refine3D_provenance_tests
integer, parameter :: BOX = 16
real, parameter    :: SMPD = 1.5

contains

    subroutine run_all_refine3D_provenance_tests()
        type(string) :: saved, root
        integer      :: nfail0
        write(*,'(A)') '**** running all refine3D provenance tests ****'
        nfail0 = tests_failed
        call enter_fixture('refine3d_provenance', saved, root)
        call test_provenance()
        call leave_fixture(saved, root, nfail0)
    end subroutine run_all_refine3D_provenance_tests

    ! Small files exercise identity and interrupted publication without a refinement run.
    subroutine test_provenance()
        use simple_parameters,          only: parameters
        use simple_fileio,              only: arr2file
        use simple_syslib,              only: del_file
        use simple_halfmap_diagnostics, only: support_provenance_fname
        use simple_refine3D_provenance, only: publish_reuse_provenance, reuse_provenance_matches, &
            &begin_reference_set, publish_reference_set, refine3D_provenance, &
            &provenance_pending_file, REFSET_OK, REFSET_INVALID, REFSET_LEGACY
        type(refine3D_provenance)      :: reference_record
        class(sp_project), allocatable :: project
        class(parameters), allocatable :: params
        class(image), allocatable      :: volume
        type(string)                   :: sigma, chain, files(3), halves(2), merged, fsc, checked_fsc, record
        character(len=256)             :: message
        integer                        :: i, status
        logical                        :: valid, support_verified
        write(*,'(A)') 'test_provenance'
        allocate(project, params, volume)
        call project%projinfo%new(1, is_ptcl=.false.)
        call project%projinfo%set(1, 'projname', 'provenance')
        call project%os_stk%new(1, is_ptcl=.false.)
        call project%os_stk%set(1, 'stk', 'particles.mrc')
        call project%os_ptcl3D%new(4, is_ptcl=.true.)
        do i = 1,4
            call project%os_ptcl3D%set(i, 'stkind', 1)
            call project%os_ptcl3D%set(i, 'indstk', i)
            call project%os_ptcl3D%set_state(i, 1)
            call project%os_ptcl3D%set(i, 'eo', mod(i,2))
            call project%os_ptcl3D%set(i, 'updatecnt', 1)
            call project%os_ptcl3D%set_euler(i, [real(i), 20., 30.])
            call project%os_ptcl3D%set_shift(i, [0., 0.])
        enddo
        params%oritype = 'ptcl3D'
        params%projfile = 'producer.simple'
        params%box_crop = BOX
        params%smpd_crop = SMPD
        params%mskdiam = 12.
        params%l_ml_reg = .false.
        params%l_nonuniform = .false.
        sigma = 'provenance_sigma.bin'
        call arr2file([1.,2.], sigma)
        call publish_reuse_provenance(params, project, project%os_ptcl3D, 'sigma', sigma, [sigma])
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'sigma', sigma, [sigma], message)
        call assert_true(valid, 'provenance: unchanged sigma producer identity')
        call project%os_ptcl3D%set(1, 'eo', 0)
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'sigma', sigma, [sigma], message)
        call assert_false(valid, 'provenance: changed half assignment')
        call project%os_ptcl3D%set(1, 'eo', 1)
        call project%os_ptcl3D%set_state(1, 0)
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'sigma', sigma, [sigma], message)
        call assert_false(valid, 'provenance: changed active population')
        call project%os_ptcl3D%set_state(1, 1)
        call project%os_ptcl3D%set(1, 'indstk', 2)
        call project%os_ptcl3D%set(2, 'indstk', 1)
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'sigma', sigma, [sigma], message)
        call assert_false(valid, 'provenance: reordered particle rows')
        call project%os_ptcl3D%set(1, 'indstk', 1)
        call project%os_ptcl3D%set(2, 'indstk', 2)
        params%mskdiam = 14.
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'sigma', sigma, [sigma], message)
        call assert_false(valid, 'provenance: changed noise mask setting')
        params%mskdiam = 12.
        call project%os_ptcl3D%set(1, 'dfx', 1.)
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'sigma', sigma, [sigma], message)
        call assert_false(valid, 'provenance: changed observation CTF')
        call project%os_ptcl3D%set(1, 'dfx', 0.)
        call project%os_ptcl3D%set_euler(1, [2.,20.,30.])
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'sigma', sigma, [sigma], message)
        call assert_true(valid, 'provenance: sigma generation may lag poses')
        call project%os_ptcl3D%set_euler(1, [1.,20.,30.])
        chain = 'provenance_chain.txt'
        files(1) = chain
        files(2) = 'provenance_chain_even.bin'
        files(3) = 'provenance_chain_odd.bin'
        do i = 1,3
            call arr2file([real(i),0.], files(i))
        enddo
        call publish_reuse_provenance(params, project, project%os_ptcl3D, 'chain', chain, files)
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'chain', chain, files, message)
        call assert_true(valid, 'provenance: unchanged chain producer identity')
        call project%os_ptcl3D%set(1, 'updatecnt', 0)
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'chain', chain, files, message)
        call assert_false(valid, 'provenance: changed represented population')
        call project%os_ptcl3D%set(1, 'updatecnt', 1)
        call project%os_ptcl3D%set_euler(1, [2.,20.,30.])
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'chain', chain, files, message)
        call assert_false(valid, 'provenance: changed chain producer poses')
        call project%os_ptcl3D%set_euler(1, [1.,20.,30.])
        call arr2file([0.,1.], chain)
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'chain', chain, files, message)
        call assert_false(valid, 'provenance: changed chain generation with equal byte count')
        merged = 'provenance_reference.mrc'
        halves(1) = 'provenance_reference_even.mrc'
        halves(2) = 'provenance_reference_odd.mrc'
        fsc = 'provenance_fsc.bin'
        call volume%new([BOX,BOX,BOX], SMPD)
        call volume%write(merged)
        call volume%write(halves(1))
        call volume%write(halves(2))
        call arr2file([1.,0.5], fsc)
        call begin_reference_set(merged)
        call publish_reference_set(params, project, project%os_ptcl3D, 1, merged, halves, fsc, 0.)
        call reference_record%new(merged, 'reference', message)
        status = reference_record%validate_reference_set(checked_fsc, support_verified, message)
        call assert_int(REFSET_OK, status, 'provenance: complete reference set')
        valid = reference_record%reference_matches_project(project, 1, .false., message)
        call assert_true(valid, 'reference: recorded continuation identity')
        call project%os_ptcl3D%set(1, 'eo', 0)
        call reference_record%new(merged, 'reference', message)
        valid = reference_record%reference_matches_project(project, 1, .true., message)
        call assert_false(valid, 'reference: explicit input cannot bypass same-lineage half changes')
        call project%os_ptcl3D%set(1, 'eo', 1)
        call project%os_ptcl3D%set_state(1, 0)
        call reference_record%new(merged, 'reference', message)
        valid = reference_record%reference_matches_project(project, 1, .false., message)
        call assert_false(valid, 'reference: changed active population')
        call project%os_ptcl3D%set_state(1, 1)
        call reference_record%new(merged, 'reference', message)
        valid = reference_record%reference_matches_project(project, 2, .false., message)
        call assert_false(valid, 'reference: different selected producer state')
        call project%os_ptcl3D%set_euler(1, [2.,20.,30.])
        call reference_record%new(merged, 'reference', message)
        valid = reference_record%reference_matches_project(project, 1, .false., message)
        call assert_true(valid, 'reference: starting poses may differ from lagged reference poses')
        call project%os_ptcl3D%set_euler(1, [1.,20.,30.])
        call project%os_ptcl3D%set(1, 'indstk', 2)
        call project%os_ptcl3D%set(2, 'indstk', 1)
        call reference_record%new(merged, 'reference', message)
        valid = reference_record%reference_matches_project(project, 1, .true., message)
        call assert_false(valid, 'reference: reordered rows are not an independent external lineage')
        call project%os_ptcl3D%set(1, 'indstk', 1)
        call project%os_ptcl3D%set(2, 'indstk', 2)
        call project%projinfo%set(1, 'projname', 'independent')
        call reference_record%new(merged, 'reference', message)
        valid = reference_record%reference_matches_project(project, 1, .false., message)
        call assert_false(valid, 'reference: registered input requires its particle lineage')
        call reference_record%new(merged, 'reference', message)
        valid = reference_record%reference_matches_project(project, 1, .true., message)
        call assert_true(valid, 'reference: independent explicit input keeps its trusted external contract')
        call project%projinfo%set(1, 'projname', 'provenance')
        record = support_provenance_fname(merged)
        call arr2file([1.], record)
        call reference_record%new(merged, 'reference', message)
        status = reference_record%validate_reference_set(checked_fsc, support_verified, message)
        call assert_int(REFSET_OK, status, 'provenance: extra sidecar does not alter the recorded base set')
        call assert_false(support_verified, 'provenance: unrecorded support sidecar is not trusted')
        call arr2file([0.5,1.], fsc)
        call reference_record%new(merged, 'reference', message)
        status = reference_record%validate_reference_set(checked_fsc, support_verified, message)
        call assert_int(REFSET_INVALID, status, 'provenance: stale FSC with equal byte count')
        call begin_reference_set(merged)
        call reference_record%new(merged, 'reference', message)
        status = reference_record%validate_reference_set(checked_fsc, support_verified, message)
        call assert_int(REFSET_INVALID, status, 'provenance: interrupted publication is not legacy')
        record = provenance_pending_file(merged)
        call del_file(record)
        call reference_record%new(merged, 'reference', message)
        status = reference_record%validate_reference_set(checked_fsc, support_verified, message)
        call assert_int(REFSET_LEGACY, status, 'provenance: unrecorded legacy base pair')
        call reference_record%kill
        call reference_record%kill
        status = reference_record%validate_reference_set(checked_fsc, support_verified, message)
        call assert_int(REFSET_INVALID, status, 'provenance: killed record cannot be reused')
        call project%kill()
        call volume%kill()
    end subroutine test_provenance
end module simple_refine3D_provenance_tester
