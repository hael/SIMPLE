!@descr: balanced particle ranking and congruent project, physical image and orientation/CTF export tests
module simple_ptcl_subset_tester
use simple_sp_project,  only: sp_project
use simple_ptcl_subset, only: select_ptcl_subset, write_ptcl_subset
use simple_string,      only: string
use simple_test_utils,  only: assert_true, assert_int, assert_real, assert_string_eq
implicit none
private
public :: run_all_ptcl_subset_tests
character(len=*), parameter :: SOURCE1 = 'ptcl_subset_source1.mrc', SOURCE2 = 'ptcl_subset_source2.mrc'
character(len=*), parameter :: SOURCE_PROJECT = 'ptcl_subset_source.simple', SUBSET_PROJECT = 'ptcl_subset_selected.simple'
character(len=*), parameter :: SUBSET_STACK = 'ptcl_subset_selected.mrc', SUBSET_TABLE = 'ptcl_subset_selected.txt'

contains

    subroutine run_all_ptcl_subset_tests()
        use simple_parameters, only: parameters
        class(sp_project), allocatable :: project
        class(parameters), allocatable :: defaults
        allocate(project,defaults)
        call assert_true(defaults%eo_rebalance=='yes','particle selection rebalances half labels by default')
        call cleanup
        call fixture(project)
        call test_selection(project)
        call test_export(project,.false.)
        call project%kill
        call cleanup
        call fixture(project)
        call test_export(project,.true.)
        call project%kill
        call cleanup
    end subroutine run_all_ptcl_subset_tests

    subroutine cleanup()
        use simple_syslib, only: del_file
        call del_file(string(SOURCE1))
        call del_file(string(SOURCE2))
        call del_file(string(SOURCE_PROJECT))
        call del_file(string(SUBSET_PROJECT))
        call del_file(string(SUBSET_STACK))
        call del_file(string(SUBSET_TABLE))
    end subroutine cleanup

    subroutine fixture(project)
        use simple_image,        only: image
        use simple_oris,         only: oris
        use simple_sym,          only: sym
        use simple_type_defs,    only: ctfparams, CTFFLAG_YES
        use simple_eo_partition, only: eo_transaction
        class(sp_project), intent(inout) :: project
        class(image), allocatable :: img
        class(oris), allocatable :: directions
        type(sym) :: symmetry
        type(ctfparams) :: ctfvars
        type(string) :: filename
        real :: pixels(8,8,1)
        integer :: stack, local_index, particle, bin
        logical :: changed
        allocate(img,directions)
        call img%new([8,8,1],1.5)
        call directions%new(8,is_ptcl=.false.)
        call symmetry%new('c1')
        call symmetry%build_refspiral(directions)
        ctfvars%ctfflag = CTFFLAG_YES
        ctfvars%smpd = 1.5
        ctfvars%fraca = .1
        ctfvars%dfx = 1.
        ctfvars%dfy = 1.2
        ctfvars%angast = 15.
        ctfvars%phshift = .25
        do stack = 1,2
            filename = SOURCE1
            if( stack == 2 ) filename = SOURCE2
            do local_index = 1,4
                pixels = real(100*stack+local_index)
                call img%set_rmat(pixels,ft=.false.)
                call img%write(filename,local_index)
            enddo
            ctfvars%kv = real(100+100*stack)
            ctfvars%cs = real(stack)
            call project%add_stk(filename,ctfvars)
        enddo
        do particle = 1,8
            local_index = mod(particle-1,4)+1
            bin = 1
            if( particle > 6 ) bin = 2
            call project%os_ptcl3D%set_euler(particle,directions%get_euler(bin))
            call project%os_ptcl3D%set_shift(particle,[real(particle)*.25,-real(particle)*.5])
            call project%os_ptcl3D%set_state(particle,2)
            call project%os_ptcl3D%set(particle,'updatecnt',1)
            call project%os_ptcl3D%set(particle,'proj',100+particle)
            call project%os_ptcl3D%set(particle,'corr',real(particle)+merge(10.,-6.,particle<=6))
            call project%os_ptcl3D%set(particle,'corr_cart',real(30-particle))
            call project%os_ptcl3D%set_dfx(particle,1.+real(particle)/10.)
            call project%os_ptcl3D%set_dfy(particle,1.2+real(particle)/10.)
            call project%os_ptcl3D%set(particle,'phshift',real(particle)/10.)
            call project%os_ptcl2D%set(particle,'indstk',5-local_index)
            call project%os_ptcl3D%set(particle,'indstk',5-local_index)
            call project%os_ptcl2D%set(particle,'pind',100+particle)
            call project%os_ptcl3D%set(particle,'pind',100+particle)
        enddo
        call eo_transaction(project,changed)
        call project%update_projinfo(string(SOURCE_PROJECT))
        call project%projinfo%set(1,'sigma2_state','parent_sigma.state')
        call project%os_out%new(1,is_ptcl=.false.)
        call project%os_out%set(1,'vol','parent_map.mrc')
        call project%write(string(SOURCE_PROJECT))
        call img%kill
        call directions%kill
        call symmetry%kill
    end subroutine fixture

    subroutine test_selection(project)
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        class(sp_project), intent(inout) :: project
        logical, allocatable :: selected(:)
        character(len=512) :: message
        integer :: status
        call select_ptcl_subset(project,4,8,'c1','corr',2,selected,status,message)
        call assert_int(0,status,'polar direction-balanced selection succeeds')
        call assert_true(all(selected .eqv. [.false.,.false.,.false.,.false.,.true.,.true.,.true.,.true.]), &
            'two particles per occupied bin, not the four largest global scores')
        call assert_int(2,project%os_ptcl3D%get_state(1),'selection leaves source states unchanged')
        call assert_int(101,project%os_ptcl3D%get_int(1,'proj'),'temporary direction labels do not replace refined projections')
        call select_ptcl_subset(project,4,8,'c1','corr_cart',2,selected,status,message)
        call assert_int(0,status,'Cartesian direction-balanced selection succeeds')
        call assert_true(all(selected .eqv. [.true.,.true.,.false.,.false.,.false.,.false.,.true.,.true.]), &
            'Cartesian score key chooses a different top-ranked subset')
        call project%os_ptcl3D%set(1,'corr_cart',0.)
        call select_ptcl_subset(project,4,8,'c1','corr_cart',2,selected,status,message)
        call assert_true(status==0 .and. .not.selected(1),'unset Cartesian scores are not treated as scored particles')
        call project%os_ptcl3D%set(1,'corr_cart',29.)
        call project%os_ptcl2D%set_state(6,0)
        call project%os_ptcl3D%set_state(5,3)
        call select_ptcl_subset(project,4,8,'c1','corr',2,selected,status,message)
        call assert_int(0,status,'selection excludes rejected particles and other structural states')
        call assert_true(all(selected .eqv. [.false.,.false.,.true.,.true.,.false.,.false.,.true.,.true.]), &
            'inactive and other-state particles are not revived')
        call project%os_ptcl2D%set_state(6,1)
        call project%os_ptcl3D%set_state(5,2)
        call project%os_ptcl3D%set(6,'corr',ieee_value(0.,ieee_quiet_nan))
        call select_ptcl_subset(project,4,8,'c1','corr',2,selected,status,message)
        call assert_true(status==0 .and. .not.selected(6),'nonfinite particle scores are excluded')
        call project%os_ptcl3D%set(6,'corr',16.)
        call project%os_ptcl3D%set(6,'updatecnt',0)
        call select_ptcl_subset(project,4,8,'c1','corr',2,selected,status,message)
        call assert_true(status==0 .and. .not.selected(6),'unrefined particles are excluded')
        call project%os_ptcl3D%set(6,'updatecnt',1)
        call select_ptcl_subset(project,9,8,'c1','corr',2,selected,status,message)
        call assert_true(status/=0 .and. .not.any(selected),'oversized targets are refused without looping')
        call select_ptcl_subset(project,0,8,'c1','corr',2,selected,status,message)
        call assert_true(status/=0,'zero target is refused')
        call select_ptcl_subset(project,4,7,'c1','corr',2,selected,status,message)
        call assert_true(status/=0,'odd direction-grid size is refused')
        call select_ptcl_subset(project,4,8,'c1','none',2,selected,status,message)
        call assert_true(status/=0,'unranked score mode is refused')
        call project%os_ptcl3D%set(1,'stkind',2)
        call select_ptcl_subset(project,4,8,'c1','corr',2,selected,status,message)
        call assert_true(status/=0,'inconsistent physical stack mapping is refused')
        call project%os_ptcl3D%set(1,'stkind',1)
    end subroutine test_selection

    subroutine test_export(project, rebalance_eo)
        use simple_image,   only: image
        use simple_oris,    only: oris
        use simple_imghead, only: find_ldim_nptcls
        use simple_syslib,  only: simple_abspath
        use simple_eo_partition, only: eo_partition_status, EO_COMPLETE
        class(sp_project), intent(inout) :: project
        logical, intent(in) :: rebalance_eo
        class(sp_project), allocatable :: restored, source
        class(oris), allocatable :: table
        class(image), allocatable :: img
        type(string) :: stackfile, tablefile, mappedfile
        logical, allocatable :: selected(:)
        real, allocatable :: pixels(:,:,:)
        real :: eulers(3,4), shifts(2,4)
        integer :: halves(4), old_halves(4), status, row, original, stack, physical_index, dims(3), nimages
        integer, parameter :: ORIGINALS(4) = [1,2,7,8], PIXEL_VALUES(4) = [104,103,202,201]
        character(len=512) :: message
        allocate(restored,source,table,img)
        do row = 1,4
            old_halves(row) = project%os_ptcl3D%get_eo(ORIGINALS(row))
        enddo
        call select_ptcl_subset(project,4,8,'c1','corr_cart',2,selected,status,message,rebalance_eo=rebalance_eo)
        call assert_int(0,status,'export fixture selection succeeds')
        if( status /= 0 ) return
        do row = 1,4
            original = ORIGINALS(row)
            halves(row) = project%os_ptcl3D%get_eo(original)
            eulers(:,row) = project%os_ptcl3D%get_euler(original)
            shifts(:,row) = project%os_ptcl3D%get_2Dshift(original)
        enddo
        if( rebalance_eo )then
            call assert_true(all(halves==[0,1,1,0]),'opt-in pairs adjacent scores into opposite halves per direction')
        else
            call assert_true(all(halves==old_halves),'explicit no preserves the original half labels')
        endif
        stackfile = simple_abspath(string(SUBSET_STACK),check_exists=.false.)
        tablefile = simple_abspath(string(SUBSET_TABLE),check_exists=.false.)
        call write_ptcl_subset(project,selected,stackfile,tablefile)
        call project%update_projinfo(string(SUBSET_PROJECT))
        call project%write(string(SUBSET_PROJECT))
        call restored%read(string(SUBSET_PROJECT))
        call source%read(string(SOURCE_PROJECT))
        call table%new(4,is_ptcl=.false.)
        call table%read(tablefile)
        call find_ldim_nptcls(stackfile,dims,nimages)
        call assert_int(4,nimages,'selected physical stack has exactly four images')
        call assert_int(4,restored%os_ptcl2D%get_noris(),'subset ptcl2D row count matches stack')
        call assert_int(4,restored%os_ptcl3D%get_noris(),'subset ptcl3D row count matches stack')
        call assert_int(4,table%get_noris(),'text table row count matches stack')
        call assert_int(2,restored%get_nstks(),'logical stacks retain their separate microscope metadata')
        call assert_int(0,restored%os_out%get_noris(),'subset does not claim the parent reconstruction')
        call assert_true(.not.restored%projinfo%isthere('sigma2_state'),'subset clears incompatible sigma state')
        call assert_int(8,source%os_ptcl3D%get_noris(),'source project file is unchanged')
        call assert_int(4,source%os_ptcl3D%get_int(5,'indstk'),'source project physical indices are unchanged')
        call assert_true(source%projinfo%isthere('sigma2_state'),'source project retains its sigma registration')
        status = eo_partition_status(restored,message)
        call assert_int(EO_COMPLETE,status,'exported subset has consistent half-set provenance')
        call img%new([8,8,1],1.5)
        do row = 1,4
            original = ORIGINALS(row)
            call restored%get_stkname_and_ind('ptcl3D',row,mappedfile,physical_index)
            call assert_string_eq(stackfile%to_char(),mappedfile,'project resolves to the generated stack')
            call assert_int(row,physical_index,'project physical index equals exported row')
            call restored%map_ptcl_ind2stk_ind('ptcl2D',row,stack,physical_index)
            call assert_int(row,physical_index,'ptcl2D physical index equals exported row')
            call assert_int(row,table%get_int(row,'indstk'),'table physical index equals exported row')
            call assert_int(100+original,table%get_int(row,'pind'),'table retains the original particle identity')
            call assert_int(2,restored%os_ptcl3D%get_state(row),'project retains structural state label')
            call assert_int(2,table%get_state(row),'text table retains structural state label')
            call assert_int(halves(row),table%get_eo(row),'table matches the selected half assignment')
            call assert_int(halves(row),restored%os_ptcl3D%get_eo(row),'project matches the selected half assignment')
            call assert_int(halves(row),restored%os_ptcl2D%get_eo(row),'both project particle segments have matching halves')
            call assert_int(old_halves(row),source%os_ptcl3D%get_eo(original),'source project halves remain unchanged')
            call assert_true(all(abs(table%get_euler(row)-eulers(:,row))<1.e-4),'table preserves Euler angles')
            call assert_true(all(abs(table%get_2Dshift(row)-shifts(:,row))<1.e-4),'table preserves shifts in native pixels')
            call assert_real(real(100+100*((original+3)/4)),table%get(row,'kv'),1.e-4, &
                'table includes the originating stack acceleration voltage')
            call assert_real(real((original+3)/4),table%get(row,'cs'),1.e-4, &
                'table includes the originating stack spherical aberration')
            call assert_real(.1,table%get(row,'fraca'),1.e-4,'table includes amplitude contrast')
            call assert_real(1.5,table%get(row,'smpd'),1.e-4,'table includes native sampling')
            call assert_real(1.+real(original)/10.,table%get_dfx(row),1.e-4,'table includes per-particle defocus')
            call assert_real(1.2+real(original)/10.,table%get_dfy(row),1.e-4,'table includes astigmatic defocus')
            call assert_real(15.,table%get(row,'angast'),1.e-4,'table includes the astigmatism angle')
            call assert_real(real(original)/10.,table%get(row,'phshift'),1.e-4,'table includes per-particle phase')
            call img%read(stackfile,row)
            pixels = img%get_rmat()
            call assert_true(all(abs(pixels-real(PIXEL_VALUES(row)))<1.e-5),'image follows original physical index, not row')
        enddo
        call img%read(string(SOURCE2),1)
        pixels = img%get_rmat()
        call assert_true(all(abs(pixels-201.)<1.e-5),'source particle stack is unchanged')
        call restored%kill
        call source%kill
        call table%kill
        call img%kill
    end subroutine test_export

end module simple_ptcl_subset_tester