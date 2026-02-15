!@descr: the user interface class (it's a beast)
module simple_user_interface
use simple_ui_all
implicit none

public :: make_user_interface, get_prg_ptr, list_simple_prgs_in_ui, list_simple_test_prgs_in_ui
public :: print_ui_json, write_ui_json, list_single_prgs_in_ui, list_stream_prgs_in_ui
public :: print_stream_ui_json, validate_ui_json, test_ui_refactoring_func
private
#include "simple_local_flags.inc"

character(len=26), parameter :: UI_FNAME = 'simple_user_interface.json'
logical,           parameter :: DEBUG    = .false.
type(ui_hash)                :: prgtab
type(string), allocatable    :: prgnames(:)
contains

    subroutine test_ui_refactoring_func
        integer :: i
        do i = 1, size(prgnames)
            call prgnames(i)%print
        end do
    end subroutine test_ui_refactoring_func
    
    ! public class methods

    subroutine make_user_interface
        call set_ui_params

        ! SIMPLE PROGRAMS
        call construct_project_programs(prgtab)
        call construct_preproc_programs(prgtab)
        call construct_cluster2D_programs(prgtab)
        call construct_cavgproc_programs(prgtab)
        call construct_abinitio3D_programs(prgtab)
        call construct_refine3D_programs(prgtab)
        call construct_denoise_programs(prgtab)
        call construct_filter_programs(prgtab)
        call construct_image_programs(prgtab)
        call construct_mask_programs(prgtab)
        call construct_ori_programs(prgtab)
        call construct_print_programs(prgtab)
        call construct_res_programs(prgtab)
        call construct_sim_programs(prgtab)
        call construct_validation_programs(prgtab)
        call construct_symmetry_programs(prgtab)
        call construct_dock_programs(prgtab)
        call construct_volume_programs(prgtab)
        call construct_other_programs(prgtab)

        ! SINGLE PROGRAMS
        call construct_single_atom_programs(prgtab)
        call construct_single_map_programs(prgtab)
        call construct_single_nano2D_programs(prgtab)
        call construct_single_nano3D_programs(prgtab)
        call construct_single_trajectory_programs(prgtab)
        call construct_single_tseries_programs(prgtab)
        call construct_single_validation_programs(prgtab)

        prgnames = prgtab%keys_sorted()
        if( DEBUG ) write(logfhandle,*) '***DEBUG::simple_user_interface; make_user_interface, DONE'
    end subroutine make_user_interface

    subroutine get_prg_ptr( which_program, ptr2prg )
        class(string),             intent(in)    :: which_program
        type(ui_program), pointer, intent(inout) :: ptr2prg
        ptr2prg => null()
        call prgtab%get_ui_program(which_program, ptr2prg)
    end subroutine get_prg_ptr

    subroutine list_simple_prgs_in_ui
        call print_project_programs(logfhandle)
        call print_preproc_programs(logfhandle)
        call print_cluster2D_programs(logfhandle)
        call print_cavgproc_programs(logfhandle)
        call print_abinitio3D_programs(logfhandle)
        call print_refine3D_programs(logfhandle)
        call print_denoise_programs(logfhandle)
        call print_filter_programs(logfhandle)
        call print_image_programs(logfhandle)
        call print_mask_programs(logfhandle)
        call print_ori_programs(logfhandle)
        call print_print_programs(logfhandle)
        call print_res_programs(logfhandle)
        call print_sim_programs(logfhandle)
        call print_validation_programs(logfhandle)
        call print_symmetry_programs(logfhandle)
        call print_dock_programs(logfhandle)
        call print_volume_programs(logfhandle)
        call print_other_programs(logfhandle)
    end subroutine list_simple_prgs_in_ui

    subroutine list_simple_test_prgs_in_ui
        !====================================================================
        ! HIGH-LEVEL
        !====================================================================
        write(logfhandle,'(A)') format_str('HIGH LEVEL:', C_UNDERLINED)
        ! write(logfhandle,'(A)') test_mini_stream%name%to_char()
        write(logfhandle,'(A)') test_sim_workflow%name%to_char()
        ! write(logfhandle,'(A)') test_inside_write%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! INPUT/OUTPUT
        !====================================================================
        write(logfhandle,'(A)') format_str('INPUT/OPTPUT:', C_UNDERLINED)
        ! write(logfhandle,'(A)') test_imgfile%name%to_char()
        ! write(logfhandle,'(A)') test_io%name%to_char()
        ! write(logfhandle,'(A)') test_io_parallel%name%to_char()
        ! write(logfhandle,'(A)') test_stack_io%name%to_char()
        ! write(logfhandle,'(A)') test_mrc_validation%name%to_char()
        ! write(logfhandle,'(A)') test_mrc2jpeg%name%to_char()
        ! write(logfhandle,'(A)') test_starfile%name%to_char()
        ! write(logfhandle,'(A)') test_star_export%name%to_char()
        ! write(logfhandle,'(A)') test_inside_write%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! NETWORK
        !====================================================================
        write(logfhandle,'(A)') format_str('NETWORK:', C_UNDERLINED)
        ! write(logfhandle,'(A)') test_socket_client%name%to_char()
        ! write(logfhandle,'(A)') test_socket_server%name%to_char()
        ! write(logfhandle,'(A)') test_socket_io%name%to_char()
        ! write(logfhandle,'(A)') test_socket_comm_distr%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! PARALLEL
        !====================================================================
        write(logfhandle,'(A)') format_str('PARALLEL:', C_UNDERLINED)
        ! write(logfhandle,'(A)') test_openmp%name%to_char()
        ! write(logfhandle,'(A)') test_openacc%name%to_char()
        ! write(logfhandle,'(A)') test_coarrays%name%to_char()
        ! write(logfhandle,'(A)') test_simd%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! FFT
        !====================================================================
        write(logfhandle,'(A)') format_str('FFT:', C_UNDERLINED)
        ! write(logfhandle,'(A)') test_phasecorr%name%to_char()
        ! write(logfhandle,'(A)') test_order_corr%name%to_char()
        ! write(logfhandle,'(A)') test_gencorrs_fft%name%to_char()
        ! write(logfhandle,'(A)') test_ft_expanded%name%to_char()
        ! write(logfhandle,'(A)') test_eval_polarftcc%name%to_char()
        ! write(logfhandle,'(A)') test_polarops%name%to_char()
        ! write(logfhandle,'(A)') test_corrs2weights%name%to_char()
        ! write(logfhandle,'(A)') test_rank_weights%name%to_char()
        ! write(logfhandle,'(A)') test_rotate_ref%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! GEOMETRY
        !====================================================================
        write(logfhandle,'(A)') format_str('GEOMETRY:', C_UNDERLINED)
        ! write(logfhandle,'(A)') test_angres%name%to_char()
        ! write(logfhandle,'(A)') test_ori%name%to_char()
        ! write(logfhandle,'(A)') test_oris%name%to_char()
        ! write(logfhandle,'(A)') test_uniform_euler%name%to_char()
        ! write(logfhandle,'(A)') test_uniform_rot%name%to_char()
        ! write(logfhandle,'(A)') test_sym%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! MASKS
        !====================================================================
        write(logfhandle,'(A)') format_str('MASKS:', C_UNDERLINED)
        ! write(logfhandle,'(A)') test_mask%name%to_char()
        ! write(logfhandle,'(A)') test_msk_routines%name%to_char()
        ! write(logfhandle,'(A)') test_otsu%name%to_char()
        ! write(logfhandle,'(A)') test_bounds_from_mask3D%name%to_char()
        ! write(logfhandle,'(A)') test_graphene_mask%name%to_char()
        ! write(logfhandle,'(A)') test_nano_mask%name%to_char()
        ! write(logfhandle,'(A)') test_ptcl_center%name%to_char()
        ! write(logfhandle,'(A)') test_image_bin%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! OPTIMIZE
        !====================================================================
        write(logfhandle,'(A)') format_str('OPTIMIZE:', C_UNDERLINED)
        ! write(logfhandle,'(A)') test_lbfgsb%name%to_char()
        ! write(logfhandle,'(A)') test_lbfgsb_cosine%name%to_char()
        ! write(logfhandle,'(A)') test_opt_lp%name%to_char()
        ! write(logfhandle,'(A)') test_lplims%name%to_char()
        ! write(logfhandle,'(A)') test_lpstages%name%to_char()
        ! write(logfhandle,'(A)') test_tree_srch%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! NUMERICS
        !====================================================================
        write(logfhandle,'(A)') format_str('NUMERICS:', C_UNDERLINED)
        ! write(logfhandle,'(A)') test_eigh%name%to_char()
        ! write(logfhandle,'(A)') test_kbinterpol_fast%name%to_char()
        ! write(logfhandle,'(A)') test_neigh%name%to_char()
        ! write(logfhandle,'(A)') test_maxnloc%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! UTILS
        !====================================================================
        write(logfhandle,'(A)') format_str('UTILS:', C_UNDERLINED)
        ! write(logfhandle,'(A)') test_cmdline%name%to_char()
        ! write(logfhandle,'(A)') test_stringmatch%name%to_char()
        ! write(logfhandle,'(A)') test_ansi_colors%name%to_char()
        ! write(logfhandle,'(A)') test_units%name%to_char()
        ! write(logfhandle,'(A)') test_serialize%name%to_char()
        ! write(logfhandle,'(A)') test_install%name%to_char()
        ! write(logfhandle,'(A)') test_nice%name%to_char()
        ! write(logfhandle,'(A)') test_binoris_io%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! STATS
        !====================================================================
        write(logfhandle,'(A)') format_str('STATS:', C_UNDERLINED) 
        ! write(logfhandle,'(A)') test_binoris%name%to_char()
        ! write(logfhandle,'(A)') test_clustering%name%to_char()
        ! write(logfhandle,'(A)') test_pca_all%name%to_char()
        ! write(logfhandle,'(A)') test_pca_imgvar%name%to_char()
        ! write(logfhandle,'(A)') test_class_sample%name%to_char()
        ! write(logfhandle,'(A)') test_multinomal%name%to_char()
        ! write(logfhandle,'(A)') test_extr_frac%name%to_char()
        ! write(logfhandle,'(A)') test_eo_diff%name%to_char()
        ! write(logfhandle,'(A)') test_ctf%name%to_char()
        ! write(logfhandle,'(A)') test_sp_project%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine list_simple_test_prgs_in_ui

    subroutine list_stream_prgs_in_ui
        call print_stream_programs(logfhandle)
    end subroutine list_stream_prgs_in_ui

    subroutine list_single_prgs_in_ui
        write(logfhandle,'(A)') format_str('PROJECT MANAGEMENT PROGRAMS:', C_UNDERLINED)
        write(logfhandle,'(A)') new_project%name%to_char()
        write(logfhandle,'(A)') update_project%name%to_char()
        write(logfhandle,'(A)') print_project_info%name%to_char()
        write(logfhandle,'(A)') print_project_field%name%to_char()
        write(logfhandle,'(A)') tseries_import%name%to_char()
        write(logfhandle,'(A)') import_particles%name%to_char()
        write(logfhandle,'(A)') import_trajectory%name%to_char()
        write(logfhandle,'(A)') prune_project%name%to_char()
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') format_str('TIME-SERIES PRE-PROCESSING PROGRAMS:', C_UNDERLINED)
        write(logfhandle,'(A)') tseries_make_pickavg%name%to_char()
        write(logfhandle,'(A)') tseries_motion_correct%name%to_char()
        write(logfhandle,'(A)') track_particles%name%to_char()
        write(logfhandle,'(A)') graphene_subtr%name%to_char()
        write(logfhandle,'(A)') trajectory_denoise%name%to_char()
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') format_str('PARTICLE 3D RECONSTRUCTION PROGRAMS:', C_UNDERLINED)
        write(logfhandle,'(A)') analysis2D_nano%name%to_char()
        write(logfhandle,'(A)') center2D_nano%name%to_char()
        write(logfhandle,'(A)') cluster2D_nano%name%to_char()
        write(logfhandle,'(A)') map_cavgs_selection%name%to_char()
        write(logfhandle,'(A)') ppca_denoise_classes%name%to_char()
        write(logfhandle,'(A)') ppca_volvar%name%to_char()
        write(logfhandle,'(A)') estimate_diam%name%to_char()
        write(logfhandle,'(A)') simulate_atoms%name%to_char()
        write(logfhandle,'(A)') refine3D_nano%name%to_char()
        write(logfhandle,'(A)') extract_substk%name%to_char()
        write(logfhandle,'(A)') extract_subproj%name%to_char()
        write(logfhandle,'(A)') autorefine3D_nano%name%to_char()
        write(logfhandle,'(A)') trajectory_reconstruct3D%name%to_char()
        write(logfhandle,'(A)') trajectory_swap_stack%name%to_char()
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') format_str('VALIDATION PROGRAMS:', C_UNDERLINED)
        write(logfhandle,'(A)') vizoris%name%to_char()
        write(logfhandle,'(A)') cavgsproc_nano%name%to_char()
        write(logfhandle,'(A)') cavgseoproc_nano%name%to_char()
        write(logfhandle,'(A)') ptclsproc_nano%name%to_char()
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') format_str('MODEL BULDING/ANALYSIS PROGRAMS:', C_UNDERLINED)
        write(logfhandle,'(A)') pdb2mrc%name%to_char()
        write(logfhandle,'(A)') conv_atom_denoise%name%to_char()
        write(logfhandle,'(A)') detect_atoms%name%to_char()
        write(logfhandle,'(A)') atoms_stats%name%to_char()
        write(logfhandle,'(A)') atoms_register%name%to_char()
        write(logfhandle,'(A)') crys_score%name%to_char()
        write(logfhandle,'(A)') atoms_rmsd%name%to_char()
        write(logfhandle,'(A)') core_atoms_analysis%name%to_char()
        write(logfhandle,'(A)') tsegmaps_core_finder%name%to_char()
        write(logfhandle,'(A)') trajectory_make_projavgs%name%to_char()
    end subroutine list_single_prgs_in_ui

    subroutine print_ui_json
        use json_module
        use simple_linked_list, only: linked_list, list_iterator
        type(json_core)           :: json
        type(json_value), pointer :: all_programs
        type(ui_program), pointer :: ptr2prg
        integer                   :: iprg
        ! JSON init
        call json%initialize()
        ! create object of program entries
        call json%create_object(all_programs, 'SIMPLE_UI')
        do iprg = 1, prgtab%count()
            call prgtab%get_ui_program(prgnames(iprg), ptr2prg)
            call create_program_entry(json, ptr2prg, all_programs)
        end do
        ! write & clean
        call json%print(all_programs, logfhandle)
        if ( json%failed() ) then
            THROW_HARD('json input/output error for simple_user_interface')
        endif
        call json%destroy(all_programs)

    contains

        subroutine create_program_entry(json, prg, all_programs)
            type(json_core),           intent(inout) :: json
            class(*),                  intent(in)    :: prg   ! we’ll select type below
            type(json_value), pointer, intent(inout) :: all_programs
            type(json_value), pointer :: program_entry, program
            ! The actual program type is whatever ptr2prg points to.
            ! We only need it to have the components you reference.
            select type (p => prg)
                class is (ui_program)
                    call json%create_object(program_entry, p%name%to_char())
                    call json%create_object(program, 'program')
                    call json%add(program_entry, program)
                    ! program section
                    call json%add(program, 'name',        p%name%to_char())
                    call json%add(program, 'descr_short', p%descr_short%to_char())
                    call json%add(program, 'descr_long',  p%descr_long%to_char())
                    call json%add(program, 'executable',  p%executable%to_char())
                    call json%add(program, 'advanced',    p%advanced)
                    if( p%gui_submenu_list%is_allocated() ) then
                        call json%add(program, 'gui_submenu_list', p%gui_submenu_list%to_char())
                    endif
                    ! all sections (now linked_list)
                    call create_section_list(json, program_entry, 'image input/output',     p%img_ios)
                    call create_section_list(json, program_entry, 'parameter input/output', p%parm_ios)
                    call create_section_list(json, program_entry, 'alternative inputs',     p%alt_ios)
                    call create_section_list(json, program_entry, 'search controls',        p%srch_ctrls)
                    call create_section_list(json, program_entry, 'filter controls',        p%filt_ctrls)
                    call create_section_list(json, program_entry, 'mask controls',          p%mask_ctrls)
                    call create_section_list(json, program_entry, 'computer controls',      p%comp_ctrls)
                    call json%add(all_programs, program_entry)
                class default
                    THROW_HARD('create_program_entry: unsupported prg dynamic type')
            end select
        end subroutine create_program_entry

        subroutine create_section_list(json, program_entry, name, lst)
            type(json_core),           intent(inout) :: json
            type(json_value), pointer, intent(inout) :: program_entry
            character(len=*),          intent(in)    :: name
            type(linked_list),         intent(in)    :: lst
            type(json_value), pointer :: entry, section
            type(list_iterator)       :: it
            class(*), allocatable     :: tmp
            character(len=STDLEN)     :: options_str, before
            character(len=KEYLEN)     :: args(10)
            integer                   :: j, nargs
            logical                   :: found, param_is_multi, param_is_binary, exception
            call json%create_array(section, trim(name))
            it = lst%begin()
            do while ( it%has_value() )
                call it%getter(tmp)
                select type (u => tmp)
                type is (ui_param)
                    call json%create_object(entry, u%key%to_char())
                    call json%add(entry, 'key',               u%key%to_char())
                    call json%add(entry, 'keytype',           u%keytype%to_char())
                    call json%add(entry, 'descr_short',       u%descr_short%to_char())
                    call json%add(entry, 'descr_long',        u%descr_long%to_char())
                    call json%add(entry, 'descr_placeholder', u%descr_placeholder%to_char())
                    call json%add(entry, 'required',          u%required)
                    if ( u%gui_submenu%is_allocated() ) then
                        call json%add(entry, 'gui_submenu', u%gui_submenu%to_char())
                    endif
                    if ( u%exclusive_group%is_allocated() ) then
                        call json%add(entry, 'exclusive_group', u%exclusive_group%to_char())
                    endif
                    if ( u%active_flags%is_allocated() ) then
                        call json%add(entry, 'active_flags', u%active_flags%to_char())
                    endif
                    if ( u%keytype%to_char() == "num" ) then
                        call json%add(entry, 'default', dble(u%rval_default))
                    else if ( u%cval_default%is_allocated() ) then
                        call json%add(entry, 'default', u%cval_default%to_char())
                    else
                        call json%add(entry, 'default', "unknown")
                    endif
                    call json%add(entry, 'advanced', u%advanced)
                    call json%add(entry, 'online',   u%online)
                    param_is_multi  = (u%keytype%to_char() .eq. 'multi')
                    param_is_binary = (u%keytype%to_char() .eq. 'binary')
                    if( param_is_multi .or. param_is_binary )then
                        options_str = u%descr_placeholder%to_char()
                        call split( options_str, '(', before )
                        call split( options_str, ')', before )
                        call parsestr(before, '|', args, nargs)
                        exception = (param_is_binary .and. nargs /= 2) .or. &
                                    (param_is_multi  .and. nargs <  2)
                        if ( exception ) then
                            write(logfhandle,*) 'Poorly formatted options string for entry ', u%key%to_char()
                            THROW_HARD(u%descr_placeholder%to_char())
                        endif
                        call json%add(entry, 'options', args(1:nargs))
                        do j = 1, nargs
                            call json%update(entry, 'options['//int2str(j)//']', trim(args(j)), found)
                        enddo
                    endif
                    call json%add(section, entry)
                class default
                    THROW_HARD('create_section_list: list item is not type(ui_param)')
                end select
                if (allocated(tmp)) deallocate(tmp)
                call it%next()
            end do
            call json%add(program_entry, section)
        end subroutine create_section_list

    end subroutine print_ui_json

    subroutine write_ui_json
        use json_module
        use simple_linked_list, only: linked_list, list_iterator
        implicit none
        type(json_core)           :: json
        type(json_value), pointer :: all_programs
        type(ui_program), pointer :: ptr2prg
        integer :: iprg
        ! JSON init
        call json%initialize()
        ! create array of program entries
        call json%create_array(all_programs, 'SIMPLE User Interface')
        do iprg = 1, prgtab%count()
            call prgtab%get_ui_program(prgnames(iprg), ptr2prg)
            call create_program_entry(json, ptr2prg, all_programs)
        end do
        ! write & clean
        call json%print(all_programs, UI_FNAME)
        if ( json%failed() ) then
            THROW_HARD('json input/output error for simple_user_interface')
        endif
        call json%destroy(all_programs)

    contains

        subroutine create_program_entry(json, prg, all_programs)
            type(json_core),           intent(inout) :: json
            class(*),                  intent(in)    :: prg   ! dynamic program object
            type(json_value), pointer, intent(inout) :: all_programs
            type(json_value), pointer :: program_entry, program
            select type (p => prg)
                class is (ui_program)
                    call json%create_object(program_entry, '')
                    call json%create_object(program, p%name%to_char())
                    call json%add(program_entry, program)
                    ! program section
                    call json%add(program, 'name',        p%name%to_char())
                    call json%add(program, 'descr_short', p%descr_short%to_char())
                    call json%add(program, 'descr_long',  p%descr_long%to_char())
                    call json%add(program, 'executable',  p%executable%to_char())
                    call json%add(program, 'advanced',    p%advanced)
                    if ( p%gui_submenu_list%is_allocated() ) then
                        call json%add(program, 'gui_submenu_list', p%gui_submenu_list%to_char())
                    endif
                    ! all sections (now linked_list)
                    call create_section_list(json, program_entry, 'image input/output',     p%img_ios)
                    call create_section_list(json, program_entry, 'parameter input/output', p%parm_ios)
                    call create_section_list(json, program_entry, 'alternative inputs',     p%alt_ios)
                    call create_section_list(json, program_entry, 'search controls',        p%srch_ctrls)
                    call create_section_list(json, program_entry, 'filter controls',        p%filt_ctrls)
                    call create_section_list(json, program_entry, 'mask controls',          p%mask_ctrls)
                    call create_section_list(json, program_entry, 'computer controls',      p%comp_ctrls)
                    call json%add(all_programs, program_entry)
                class default
                    THROW_HARD('create_program_entry: unsupported prg dynamic type')
            end select
        end subroutine create_program_entry

        subroutine create_section_list(json, program_entry, name, lst)
            type(json_core),           intent(inout) :: json
            type(json_value), pointer, intent(inout) :: program_entry
            character(len=*),          intent(in)    :: name
            type(linked_list),         intent(in)    :: lst
            type(json_value), pointer :: entry, section
            type(list_iterator)       :: it
            class(*), allocatable     :: tmp
            character(len=STDLEN)     :: options_str, before
            character(len=KEYLEN)     :: args(10)
            integer                   :: j, nargs
            logical                   :: found, param_is_multi, param_is_binary, exception
            call json%create_array(section, trim(name))
            it = lst%begin()
            do while ( it%has_value() )
                call it%getter(tmp)
                select type (u => tmp)
                type is (ui_param)
                    call json%create_object(entry, u%key%to_char())
                    call json%add(entry, 'key',               u%key%to_char())
                    call json%add(entry, 'keytype',           u%keytype%to_char())
                    call json%add(entry, 'descr_short',       u%descr_short%to_char())
                    call json%add(entry, 'descr_long',        u%descr_long%to_char())
                    call json%add(entry, 'descr_placeholder', u%descr_placeholder%to_char())
                    call json%add(entry, 'required',          u%required)
                    if ( u%gui_submenu%is_allocated() ) then
                        call json%add(entry, 'gui_submenu', u%gui_submenu%to_char())
                    endif
                    if ( u%exclusive_group%is_allocated() ) then
                        call json%add(entry, 'exclusive_group', u%exclusive_group%to_char())
                    endif
                    if ( u%active_flags%is_allocated() ) then
                        call json%add(entry, 'active_flags', u%active_flags%to_char())
                    endif
                    call json%add(entry, 'advanced', u%advanced)
                    call json%add(entry, 'online',   u%online)
                    if ( u%keytype%to_char() == "num" ) then
                        call json%add(entry, 'default', dble(u%rval_default))
                    else if ( u%cval_default%is_allocated() ) then
                        call json%add(entry, 'default', u%cval_default%to_char())
                    else
                        call json%add(entry, 'default', "unknown")
                    endif
                    param_is_multi  = (u%keytype%to_char() .eq. 'multi')
                    param_is_binary = (u%keytype%to_char() .eq. 'binary')
                    if( param_is_multi .or. param_is_binary ) then
                        options_str = u%descr_placeholder%to_char()
                        call split( options_str, '(', before )
                        call split( options_str, ')', before )
                        call parsestr(before, '|', args, nargs)
                        exception = (param_is_binary .and. nargs /= 2) .or. &
                                    (param_is_multi  .and. nargs <  2)
                        if( exception )then
                            write(logfhandle,*) 'Poorly formatted options string for entry ', u%key%to_char()
                            THROW_HARD(u%descr_placeholder%to_char())
                        endif
                        call json%add(entry, 'options', args(1:nargs))
                        do j = 1, nargs
                            call json%update(entry, 'options['//int2str(j)//']', trim(args(j)), found)
                        enddo
                    endif
                    call json%add(section, entry)
                class default
                    THROW_HARD('create_section_list: list item is not type(ui_param)')
                end select
                if (allocated(tmp)) deallocate(tmp)
                call it%next()
            end do
            call json%add(program_entry, section)
        end subroutine create_section_list

    end subroutine write_ui_json

    subroutine print_stream_ui_json
        ! this is ugly at the moment but paves the way ....
        use json_module
        type(json_core)           :: json
        type(json_value), pointer :: input, user_inputs, ui
        type(json_value), pointer :: processes, process, process_inputs
        ! JSON init
        call json%initialize()
        ! create object of program entries
        call json%create_object(ui, 'STREAM_UI')
        ! master inputs
        call json%create_array(user_inputs, 'user_inputs')
        call json%add(ui, user_inputs)
        !! dir_movies
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'dir_movies')
        call json%add(input, 'keytype',     'dir')
        call json%add(input, 'descr_short', 'Input movies directory')
        call json%add(input, 'descr_long',  'Input movies directory')
        call json%add(input, 'required',    .TRUE.)
        !! dir_meta
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'dir_meta')
        call json%add(input, 'keytype',     'dir')
        call json%add(input, 'descr_short', 'Input metadata directory')
        call json%add(input, 'descr_long',  'Input metadata directory')
        call json%add(input, 'required',    .FALSE.)
        !! gainref
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'gainref')
        call json%add(input, 'keytype',     'file')
        call json%add(input, 'descr_short', 'Gain reference')
        call json%add(input, 'descr_long',  'Gain reference')
        call json%add(input, 'required',    .FALSE.)
        !! cs
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'cs')
        call json%add(input, 'keytype',     'float')
        call json%add(input, 'descr_short', 'Spherical aberration (mm)')
        call json%add(input, 'descr_long',  'Spherical aberration (mm)')
        call json%add(input, 'required',    .TRUE.)
        call json%add(input, 'default',     STREAM_DEFAULT_CS)
        !! fraca
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'fraca')
        call json%add(input, 'keytype',     'float')
        call json%add(input, 'descr_short', 'Amplitude contrast fraction')
        call json%add(input, 'descr_long',  'Amplitude contrast fraction')
        call json%add(input, 'required',    .TRUE.)
        call json%add(input, 'default',     STREAM_DEFAULT_FRACA)
        !! kv
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'kv')
        call json%add(input, 'keytype',     'int')
        call json%add(input, 'descr_short', 'Acceleration voltage (kV)')
        call json%add(input, 'descr_long',  'Acceleration voltage (kV)')
        call json%add(input, 'required',    .TRUE.)
        call json%add(input, 'default',     int2str(STREAM_DEFAULT_KV))
        !! smpd
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'smpd')
        call json%add(input, 'keytype',     'float')
        call json%add(input, 'descr_short', 'Pixel size (A)')
        call json%add(input, 'descr_long',  'Pixel size (A)')
        call json%add(input, 'required',    .TRUE.)
        !! smpd_downscale
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'smpd_downscale')
        call json%add(input, 'keytype',     'hidden')
        call json%add(input, 'descr_short', 'downscale pixel size (A)')
        call json%add(input, 'descr_long',  'downscale pixel size (A)')
        call json%add(input, 'required',    .TRUE.)
        call json%add(input, 'default',     real2str(SMPD4DOWNSCALE))
        !! total_dose
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'total_dose')
        call json%add(input, 'keytype',     'float')
        call json%add(input, 'descr_short', 'Total exposure dose (e/A2)')
        call json%add(input, 'descr_long',  'Total exposure dose (e/A2)')
        call json%add(input, 'required',    .TRUE.)
        !! pickrefs
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'pickrefs')
        call json%add(input, 'keytype',     'file')
        call json%add(input, 'descr_short', '2D averages for use as picking references (optional)')
        call json%add(input, 'descr_long',  '2D averages for use as picking references (optional)')
        call json%add(input, 'required',    .FALSE.)
        !! box size
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'box_extract')
        call json%add(input, 'keytype',     'int')
        call json%add(input, 'descr_short', 'force box size (px, optional)')
        call json%add(input, 'descr_long',  'force a box size (px) eg. to match an existing dataset')
        call json%add(input, 'required',    .FALSE.)
        ! programs
        call json%create_array(processes, 'processes')
        call json%add(ui, processes)
        !! preproc
        call json%create_object(process, 'process')
        call json%add(processes, process)
        call json%add(process, 'name',         PREPROC_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process, 'prg',          'preproc')
        call json%add(process, 'nthr_master',  DEFAULT_NTHR_MASTER)
        call json%create_array(process_inputs, 'user_inputs')
        call json%add(process, process_inputs)
        call json%add(process_inputs, '', 'dir_movies')
        call json%add(process_inputs, '', 'dir_meta')
        call json%add(process_inputs, '', 'gainref')
        call json%add(process_inputs, '', 'cs')
        call json%add(process_inputs, '', 'fraca')
        call json%add(process_inputs, '', 'kv')
        call json%add(process_inputs, '', 'smpd')
        call json%add(process_inputs, '', 'scale')
        call json%add(process_inputs, '', 'total_dose')
        call json%add(process_inputs, '', 'smpd_downscale')
        call json%create_array(process_inputs, 'static_inputs')
        call json%add(process, process_inputs)
        call json%add(process_inputs, '', 'outdir='   // PREPROC_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process_inputs, '', 'nparts='   // int2str(PREPROC_NPARTS))
        call json%add(process_inputs, '', 'nthr='     // int2str(PREPROC_NTHR))
        call json%add(process_inputs, '', 'ninipick=' // int2str(PREPROC_NINIPICK))
        !! assign_optics
        call json%create_object(process, 'process')
        call json%add(processes, process)
        call json%add(process, 'name',         OPTICS_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process, 'prg',          'assign_optics')
        call json%add(process, 'nthr_master',  DEFAULT_NTHR_MASTER)
        call json%create_array(process_inputs, 'static_inputs')
        call json%add(process, process_inputs)
        call json%add(process_inputs, '', 'dir_target=' // PREPROC_JOB_NAME)
        call json%add(process_inputs, '', 'outdir='     // OPTICS_JOB_NAME) !important - directory names and name must match between processes
        !! opening 2D
        call json%create_object(process, 'process')
        call json%add(processes, process)
        call json%add(process, 'name',         OPENING2D_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process, 'prg',          'gen_pickrefs')
        call json%add(process, 'nthr_master',  OPENING2D_NTHR)
        call json%create_array(process_inputs, 'static_inputs')
        call json%add(process, process_inputs)
        call json%add(process_inputs, '', 'dir_target='    // PREPROC_JOB_NAME)
        call json%add(process_inputs, '', 'optics_dir=../' // OPTICS_JOB_NAME)
        call json%add(process_inputs, '', 'outdir='        // OPENING2D_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process_inputs, '', 'nthr='          // int2str(OPENING2D_NTHR))
        !! reference_based_picking
        call json%create_object(process, 'process')
        call json%add(processes, process)
        call json%add(process, 'name',         REFPICK_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process, 'prg',          'pick_extract')
        call json%add(process, 'nthr_master',  DEFAULT_NTHR_MASTER)
        call json%create_array(process_inputs, 'user_inputs')
        call json%add(process, process_inputs)
        call json%add(process_inputs, '', 'pickrefs')
        call json%add(process_inputs, '', 'box_extract')
        call json%create_array(process_inputs, 'static_inputs')
        call json%add(process, process_inputs)
        call json%add(process_inputs, '', 'dir_target='    // PREPROC_JOB_NAME)
        call json%add(process_inputs, '', 'optics_dir=../' // OPTICS_JOB_NAME)
        call json%add(process_inputs, '', 'outdir='        // REFPICK_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process_inputs, '', 'nparts='        // int2str(REFPICK_NPARTS))
        call json%add(process_inputs, '', 'nthr='          // int2str(REFPICK_NTHR))
        !! particle_sieving
        call json%create_object(process, 'process')
        call json%add(processes, process)
        call json%add(process, 'name',         SIEVING_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process, 'prg',          'sieve_cavgs')
        call json%add(process, 'nthr_master',  DEFAULT_NTHR_MASTER)
        call json%create_array(process_inputs, 'static_inputs')
        call json%add(process, process_inputs)
        call json%add(process_inputs, '', 'dir_target='    // REFPICK_JOB_NAME)
        call json%add(process_inputs, '', 'optics_dir=../' // OPTICS_JOB_NAME)
        call json%add(process_inputs, '', 'outdir='        // SIEVING_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process_inputs, '', 'ncls='          // int2str(SIEVING_NCLS))
        call json%add(process_inputs, '', 'nptcls_per_cls='// int2str(SIEVING_NPTCLS_PER_CLASS))
        call json%add(process_inputs, '', 'nchunks='       // int2str(SIEVING_NCHUNKS))
        call json%add(process_inputs, '', 'nparts='        // int2str(SIEVING_NPARTS))
        call json%add(process_inputs, '', 'nthr='          // int2str(SIEVING_NTHR))
        call json%add(process_inputs, '', 'interactive=yes')
        !! 2D classification
        call json%create_object(process, 'process')
        call json%add(processes, process)
        call json%add(process, 'name',         CLASS2D_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process, 'prg',          'abinitio2D_stream')
        call json%add(process, 'nthr_master',  DEFAULT_NTHR_MASTER)
        call json%create_array(process_inputs, 'static_inputs')
        call json%add(process, process_inputs)
        call json%add(process_inputs, '', 'dir_target='    // SIEVING_JOB_NAME)
        call json%add(process_inputs, '', 'optics_dir=../' // OPTICS_JOB_NAME)
        call json%add(process_inputs, '', 'outdir='        // CLASS2D_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process_inputs, '', 'ncls='          // int2str(CLASS2D_NCLS))
        call json%add(process_inputs, '', 'nparts='        // int2str(CLASS2D_NPARTS))
        call json%add(process_inputs, '', 'nthr='          // int2str(CLASS2D_NTHR))
        ! print & clean
        call json%print(ui, logfhandle)
        if( json%failed() )then
            THROW_HARD('json input/output error for simple_user_interface')
        endif
        call json%destroy(ui)
    end subroutine print_stream_ui_json

    subroutine validate_ui_json
        use json_module
        use json_kinds
        use json_file_module
        type(json_core)                       :: json
        type(json_value),             pointer :: p
        character(kind=CK,len=:), allocatable :: fname
        ! Builds & writes UI
        call make_user_interface
        write(*,*)'Constructed UI'
        call write_ui_json
        write(*,*)'Wrote UI'
        ! Parses all values, check for duplicate and invalid types
        fname = UI_FNAME
        call json%parse(fname, p)
        write(*,*)'Completed json_parse_file'
        ! Cleanup
        call json%destroy()
        nullify(p)
    end subroutine validate_ui_json

end module simple_user_interface
