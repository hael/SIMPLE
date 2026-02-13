!@descr: the user interface class (it's a beast)
module simple_user_interface
use simple_ui_api_all
implicit none

public :: make_user_interface, get_prg_ptr, list_simple_prgs_in_ui, list_simple_test_prgs_in_ui
public :: print_ui_json, write_ui_json, list_single_prgs_in_ui, list_stream_prgs_in_ui
public :: print_stream_ui_json, validate_ui_json, test_ui_refactoring_func
private
#include "simple_local_flags.inc"

character(len=26), parameter :: UI_FNAME = 'simple_user_interface.json'
logical,           parameter :: DEBUG    = .false.

! this is for making an array of pointers to all programs
type simple_prg_ptr
    type(ui_program), pointer :: ptr2prg => null()
end type simple_prg_ptr
integer, parameter   :: NMAX_PTRS  = 200
integer              :: n_prg_ptrs = 0
type(simple_prg_ptr) :: prg_ptr_array(NMAX_PTRS)

type(ui_hash) :: prgtab_other
type(ui_hash) :: prgtab_all

contains

    subroutine test_ui_refactoring_func
        type(string), allocatable :: prgnames(:)
        integer :: i
        prgnames = prgtab_all%keys_sorted()
        do i = 1, size(prgnames)
            call prgnames(i)%print
        end do
    end subroutine test_ui_refactoring_func
    
    ! public class methods

    subroutine make_user_interface
        call set_ui_params
        call set_prg_ptr_array
        call new_abinitio2D( prgtab_all )
        call new_abinitio2D_stream( prgtab_all )
        call new_abinitio3D( prgtab_all )
        call new_abinitio3D_cavgs( prgtab_all )
        call new_analysis2D_nano
        call new_assign_optics( prgtab_all )
        call new_assign_optics_groups( prgtab_all ) 
        call new_atoms_register
        call new_atoms_stats
        call new_auto_spher_mask( prgtab_all )
        call new_automask( prgtab_all )
        call new_automask2D( prgtab_all )
        call new_autorefine3D_nano
        call new_binarize( prgtab_all )
        call new_cavgseoproc_nano
        call new_cavgsproc_nano
        call new_center( prgtab_all )
        call new_center2D_nano
        call new_check_refpick( prgtab_all )
        call new_cleanup2D( prgtab_all )
        call new_clin_fsc( prgtab_all )
        call new_cluster2D( prgtab_all )
        call new_cluster2D_nano
        call new_cluster2D_stream( prgtab_all )
        call new_cluster2D_subsets( prgtab_all )
        call new_cluster_cavgs( prgtab_all )
        call new_cluster_cavgs_selection( prgtab_all )
        call new_cluster_stack( prgtab_all )
        call new_conv_atom_denoise
        call new_convert( prgtab_all )
        call new_crys_score
        call new_ctf_estimate( prgtab_all )
        call new_ctf_phaseflip( prgtab_all )
        call new_ctfops( prgtab_all )
        call new_trajectory_denoise
        call new_detect_atoms
        call new_dock_volpair( prgtab_all )
        call new_estimate_diam
        call new_estimate_lpstages( prgtab_all )
        call new_export_relion( prgtab_all )
        call new_export_starproject( prgtab_all )
        call new_extract( prgtab_all )
        call new_extract_subproj( prgtab_all )
        call new_extract_substk
        call new_filter( prgtab_all )
        call new_fsc( prgtab_all )
        call new_gen_pickrefs( prgtab_all )
        call new_gen_pspecs_and_thumbs( prgtab_all )
        call new_graphene_subtr
        call new_icm2D( prgtab_all )
        call new_icm3D( prgtab_all )
        call new_import_boxes( prgtab_all )
        call new_import_cavgs( prgtab_all )
        call new_import_movies( prgtab_all )
        call new_import_particles( prgtab_all )
        call new_import_starproject( prgtab_all )
        call new_info_image( prgtab_all )
        call new_info_stktab( prgtab_all )
        call new_make_cavgs( prgtab_all )
        call new_make_oris( prgtab_all )
        call new_map_cavgs_selection( prgtab_all )
        call new_mask( prgtab_all )
        call new_match_cavgs( prgtab_all )
        call new_match_stacks( prgtab_all )
        call new_merge_projects( prgtab_all )
        call new_mini_stream( prgtab_all )
        call new_mkdir_( prgtab_all )
        call new_model_validation( prgtab_all )
        call new_motion_correct( prgtab_all )
        call new_multivol_assign( prgtab_all )
        call new_new_project( prgtab_all )
        call new_noisevol( prgtab_all )
        call new_normalize( prgtab_all )
        call new_orisops( prgtab_all )
        call new_oristats( prgtab_all )
        call new_pdb2mrc( prgtab_all )
        call new_pick( prgtab_all )
        call new_pick_extract( prgtab_all )
        call new_postprocess( prgtab_all )
        call new_ppca_denoise( prgtab_all )
        call new_ppca_denoise_classes( prgtab_all )
        call new_ppca_volvar( prgtab_all )
        call new_preproc( prgtab_all )
        call new_preprocess( prgtab_all )
        call new_print_dose_weights( prgtab_all )
        call new_print_fsc( prgtab_all )
        call new_print_magic_boxes( prgtab_all )
        call new_print_project_field( prgtab_all )
        call new_print_project_info( prgtab_all )
        call new_prune_project( prgtab_all )
        call new_ptclsproc_nano
        call new_reconstruct3D( prgtab_all )
        call new_reextract( prgtab_all )
        call new_refine3D( prgtab_all )
        call new_refine3D_auto( prgtab_all )
        call new_refine3D_nano
        call new_replace_project_field( prgtab_all )
        call new_reproject( prgtab_all )
        call new_sample_classes( prgtab_all )
        call new_scale( prgtab_all )
        call new_select_( prgtab_all )
        call new_select_clusters( prgtab_all )
        call new_selection( prgtab_all )
        call new_sieve_cavgs( prgtab_all )
        call new_simulate_atoms
        call new_simulate_movie( prgtab_all )
        call new_simulate_noise( prgtab_all )
        call new_simulate_particles( prgtab_all )
        call new_split_( prgtab_all )
        call new_split_stack( prgtab_all )
        call new_stack( prgtab_all )
        call new_stackops( prgtab_all )
        call new_symaxis_search( prgtab_all )
        call new_symmetrize_map( prgtab_all )
        call new_symmetry_test( prgtab_all )
        call new_atoms_rmsd
        call new_core_atoms_analysis
        call new_tsegmaps_core_finder
        call new_tseries_import
        call new_import_trajectory
        call new_tseries_make_pickavg
        call new_trajectory_make_projavgs
        call new_tseries_motion_correct
        call new_trajectory_reconstruct3D
        call new_trajectory_swap_stack
        call new_track_particles
        call new_uniform_filter2D( prgtab_all )
        call new_uniform_filter3D( prgtab_all )
        call new_update_project( prgtab_all )
        call new_vizoris( prgtab_all )
        call new_volanalyze( prgtab_all )
        call new_volops( prgtab_all )
        call new_write_classes( prgtab_all )
        call new_write_mic_filetab( prgtab_all )
        call new_zero_project_shifts( prgtab_all )
        ! test programs
        call new_test_sim_workflow( prgtab_all )
        if( DEBUG ) write(logfhandle,*) '***DEBUG::simple_user_interface; make_user_interface, DONE'
    end subroutine make_user_interface

    subroutine set_prg_ptr_array
        n_prg_ptrs = 0        
        call push2prg_ptr_array(abinitio2D)
        call push2prg_ptr_array(abinitio2D_stream)
        call push2prg_ptr_array(abinitio3D)
        call push2prg_ptr_array(abinitio3D_cavgs)
        call push2prg_ptr_array(analysis2D_nano)
        call push2prg_ptr_array(assign_optics_groups)
        call push2prg_ptr_array(atoms_register)
        call push2prg_ptr_array(atoms_stats)
        call push2prg_ptr_array(auto_spher_mask)
        call push2prg_ptr_array(automask)
        call push2prg_ptr_array(automask2D)
        call push2prg_ptr_array(autorefine3D_nano)
        call push2prg_ptr_array(binarize)
        call push2prg_ptr_array(cavgseoproc_nano)
        call push2prg_ptr_array(cavgsproc_nano)
        call push2prg_ptr_array(center)
        call push2prg_ptr_array(center2D_nano)
        call push2prg_ptr_array(check_refpick)
        call push2prg_ptr_array(clin_fsc)
        call push2prg_ptr_array(cleanup2D)
        call push2prg_ptr_array(cluster2D)
        call push2prg_ptr_array(cluster2D_nano)
        call push2prg_ptr_array(cluster2D_stream)
        call push2prg_ptr_array(cluster2D_subsets)
        call push2prg_ptr_array(cluster_cavgs)
        call push2prg_ptr_array(cluster_cavgs_selection)
        call push2prg_ptr_array(cluster_stack)
        call push2prg_ptr_array(conv_atom_denoise)
        call push2prg_ptr_array(convert)
        call push2prg_ptr_array(crys_score)
        call push2prg_ptr_array(ctf_estimate)
        call push2prg_ptr_array(ctf_phaseflip)
        call push2prg_ptr_array(ctfops)
        call push2prg_ptr_array(trajectory_denoise)
        call push2prg_ptr_array(detect_atoms)
        call push2prg_ptr_array(dock_volpair)
        call push2prg_ptr_array(export_relion)
        call push2prg_ptr_array(export_starproject)
        call push2prg_ptr_array(extract)
        call push2prg_ptr_array(extract_subproj)
        call push2prg_ptr_array(extract_substk)
        call push2prg_ptr_array(filter)
        call push2prg_ptr_array(fsc)
        call push2prg_ptr_array(gen_pickrefs)
        call push2prg_ptr_array(gen_pspecs_and_thumbs)
        call push2prg_ptr_array(graphene_subtr)
        call push2prg_ptr_array(icm2D)
        call push2prg_ptr_array(icm3D)
        call push2prg_ptr_array(import_boxes)
        call push2prg_ptr_array(import_cavgs)
        call push2prg_ptr_array(import_movies)
        call push2prg_ptr_array(import_particles)
        call push2prg_ptr_array(import_starproject)
        call push2prg_ptr_array(info_image)
        call push2prg_ptr_array(info_stktab)
        call push2prg_ptr_array(make_cavgs)
        call push2prg_ptr_array(make_oris)
        call push2prg_ptr_array(map_cavgs_selection)
        call push2prg_ptr_array(mask)
        call push2prg_ptr_array(match_cavgs)
        call push2prg_ptr_array(match_stacks)
        call push2prg_ptr_array(merge_projects)
        call push2prg_ptr_array(mini_stream)
        call push2prg_ptr_array(mkdir_)
        call push2prg_ptr_array(model_validation)
        call push2prg_ptr_array(motion_correct)
        call push2prg_ptr_array(multivol_assign)
        call push2prg_ptr_array(new_project)
        call push2prg_ptr_array(noisevol)
        call push2prg_ptr_array(normalize_)
        call push2prg_ptr_array(orisops)
        call push2prg_ptr_array(oristats)
        call push2prg_ptr_array(pdb2mrc)
        call push2prg_ptr_array(pick)
        call push2prg_ptr_array(pick_extract)
        call push2prg_ptr_array(postprocess)
        call push2prg_ptr_array(ppca_denoise)
        call push2prg_ptr_array(ppca_denoise_classes)
        call push2prg_ptr_array(ppca_volvar)
        call push2prg_ptr_array(preproc)
        call push2prg_ptr_array(preprocess)
        call push2prg_ptr_array(print_dose_weights)
        call push2prg_ptr_array(print_fsc)
        call push2prg_ptr_array(print_magic_boxes)
        call push2prg_ptr_array(print_project_field)
        call push2prg_ptr_array(print_project_info)
        call push2prg_ptr_array(prune_project)
        call push2prg_ptr_array(ptclsproc_nano)
        call push2prg_ptr_array(reconstruct3D)
        call push2prg_ptr_array(reextract)
        call push2prg_ptr_array(refine3D)
        call push2prg_ptr_array(refine3D_auto)
        call push2prg_ptr_array(refine3D_nano)
        call push2prg_ptr_array(replace_project_field)
        call push2prg_ptr_array(reproject)
        call push2prg_ptr_array(sample_classes)
        call push2prg_ptr_array(scale)
        call push2prg_ptr_array(select_)
        call push2prg_ptr_array(select_clusters)
        call push2prg_ptr_array(selection)
        call push2prg_ptr_array(sieve_cavgs)
        call push2prg_ptr_array(simulate_atoms)
        call push2prg_ptr_array(simulate_movie)
        call push2prg_ptr_array(simulate_noise)
        call push2prg_ptr_array(simulate_particles)
        call push2prg_ptr_array(split_)
        call push2prg_ptr_array(split_stack)
        call push2prg_ptr_array(stack)
        call push2prg_ptr_array(stackops)
        call push2prg_ptr_array(symaxis_search)
        call push2prg_ptr_array(symmetrize_map)
        call push2prg_ptr_array(symmetry_test)
        call push2prg_ptr_array(atoms_rmsd)
        call push2prg_ptr_array(core_atoms_analysis)
        call push2prg_ptr_array(tsegmaps_core_finder)
        call push2prg_ptr_array(tseries_import)
        call push2prg_ptr_array(import_trajectory)
        call push2prg_ptr_array(tseries_make_pickavg)
        call push2prg_ptr_array(tseries_motion_correct)
        call push2prg_ptr_array(trajectory_reconstruct3D)
        call push2prg_ptr_array(trajectory_swap_stack)
        call push2prg_ptr_array(track_particles)
        call push2prg_ptr_array(uniform_filter2D)
        call push2prg_ptr_array(uniform_filter3D)
        call push2prg_ptr_array(update_project)
        call push2prg_ptr_array(vizoris)
        call push2prg_ptr_array(volanalyze)
        call push2prg_ptr_array(volops)
        call push2prg_ptr_array(write_classes)
        call push2prg_ptr_array(write_mic_filetab)
        call push2prg_ptr_array(zero_project_shifts)

        call push2prg_ptr_array(test_sim_workflow)
        if( DEBUG ) write(logfhandle,*) '***DEBUG::simple_user_interface; set_prg_ptr_array, DONE'
        
        contains

            subroutine push2prg_ptr_array( prg )
                type(ui_program), target :: prg
                n_prg_ptrs = n_prg_ptrs + 1
                prg_ptr_array(n_prg_ptrs)%ptr2prg => prg
            end subroutine push2prg_ptr_array

    end subroutine set_prg_ptr_array

    subroutine get_prg_ptr( which_program, ptr2prg )
        class(string), intent(in) :: which_program
        type(ui_program), pointer :: ptr2prg
        select case(which_program%to_char())
            case('abinitio2D');                  ptr2prg => abinitio2D
            case('abinitio2D_stream');           ptr2prg => abinitio2D_stream
            case('abinitio3D');                  ptr2prg => abinitio3D
            case('abinitio3D_cavgs');            ptr2prg => abinitio3D_cavgs
            case('analysis2D_nano');             ptr2prg => analysis2D_nano
            case('assign_optics');               ptr2prg => assign_optics   
            case('assign_optics_groups');        ptr2prg => assign_optics_groups
            case('atoms_register');              ptr2prg => atoms_register
            case('atoms_stats');                 ptr2prg => atoms_stats
            case('auto_spher_mask');             ptr2prg => auto_spher_mask
            case('automask');                    ptr2prg => automask
            case('automask2D');                  ptr2prg => automask2D
            case('autorefine3D_nano');           ptr2prg => autorefine3D_nano
            case('binarize');                    ptr2prg => binarize
            case('cavgseoproc_nano');            ptr2prg => cavgseoproc_nano
            case('cavgsproc_nano');              ptr2prg => cavgsproc_nano
            case('center');                      ptr2prg => center
            case('center2D_nano');               ptr2prg => center2D_nano
            case('check_refpick');               ptr2prg => check_refpick
            case('cleanup2D');                   ptr2prg => cleanup2D    
            case('clin_fsc');                    ptr2prg => clin_fsc
            case('cluster2D');                   ptr2prg => cluster2D
            case('cluster2D_nano');              ptr2prg => cluster2D_nano
            case('cluster2D_stream');            ptr2prg => cluster2D_stream
            case('cluster2D_subsets');           ptr2prg => cluster2D_subsets
            case('cluster_cavgs');               ptr2prg => cluster_cavgs
            case('cluster_cavgs_selection');     ptr2prg => cluster_cavgs_selection
            case('cluster_stack');               ptr2prg => cluster_stack
            case('conv_atom_denoise');           ptr2prg => conv_atom_denoise
            case('convert');                     ptr2prg => convert
            case('crys_score');                  ptr2prg => crys_score
            case('ctf_estimate');                ptr2prg => ctf_estimate
            case('ctf_phaseflip');               ptr2prg => ctf_phaseflip
            case('ctfops');                      ptr2prg => ctfops
            case('trajectory_denoise');          ptr2prg => trajectory_denoise
            case('detect_atoms');                ptr2prg => detect_atoms
            case('dock_volpair');                ptr2prg => dock_volpair
            case('estimate_diam');               ptr2prg => estimate_diam
            case('estimate_lpstages');           ptr2prg => estimate_lpstages
            case('export_relion');               ptr2prg => export_relion
            case('export_starproject');          ptr2prg => export_starproject
            case('extract');                     ptr2prg => extract
            case('extract_subproj');             ptr2prg => extract_subproj
            case('extract_substk');              ptr2prg => extract_substk
            case('filter');                      ptr2prg => filter
            case('fsc');                         ptr2prg => fsc
            case('gen_pickrefs');                ptr2prg => gen_pickrefs
            case('gen_pspecs_and_thumbs');       ptr2prg => gen_pspecs_and_thumbs
            case('graphene_subtr');              ptr2prg => graphene_subtr
            case('icm2D');                       ptr2prg => icm2D
            case('icm3D');                       ptr2prg => icm3D
            case('import_boxes');                ptr2prg => import_boxes
            case('import_cavgs');                ptr2prg => import_cavgs
            case('import_movies');               ptr2prg => import_movies
            case('import_particles');            ptr2prg => import_particles
            case('import_starproject');          ptr2prg => import_starproject
            case('info_image');                  ptr2prg => info_image
            case('info_stktab');                 ptr2prg => info_stktab
            case('make_cavgs');                  ptr2prg => make_cavgs
            case('make_oris');                   ptr2prg => make_oris
            case('map_cavgs_selection');         ptr2prg => map_cavgs_selection
            case('mask');                        ptr2prg => mask
            case('match_cavgs');                 ptr2prg => match_cavgs
            case('match_stacks');                ptr2prg => match_stacks
            case('merge_projects');              ptr2prg => merge_projects
            case('mini_stream');                 ptr2prg => mini_stream
            case('mkdir');                       ptr2prg => mkdir_
            case('motion_correct');              ptr2prg => motion_correct
            case('model_validation');            ptr2prg => model_validation
            case('multivol_assign');             ptr2prg => multivol_assign
            case('new_project');                 ptr2prg => new_project
            case('noisevol');                    ptr2prg => noisevol
            case('normalize');                   ptr2prg => normalize_
            case('orisops');                     ptr2prg => orisops
            case('oristats');                    ptr2prg => oristats
            case('pdb2mrc');                     ptr2prg => pdb2mrc
            case('pick');                        ptr2prg => pick
            case('pick_extract');                ptr2prg => pick_extract
            case('postprocess');                 ptr2prg => postprocess
            case('ppca_denoise');                ptr2prg => ppca_denoise
            case('ppca_denoise_classes');        ptr2prg => ppca_denoise_classes
            case('ppca_volvar');                 ptr2prg => ppca_volvar
            case('preproc');                     ptr2prg => preproc
            case('preprocess');                  ptr2prg => preprocess
            case('print_dose_weights');          ptr2prg => print_dose_weights
            case('print_fsc');                   ptr2prg => print_fsc
            case('print_magic_boxes');           ptr2prg => print_magic_boxes
            case('print_project_field');         ptr2prg => print_project_field
            case('print_project_info');          ptr2prg => print_project_info
            case('prune_project');               ptr2prg => prune_project
            case('ptclsproc_nano');              ptr2prg => ptclsproc_nano
            case('reconstruct3D');               ptr2prg => reconstruct3D
            case('reextract');                   ptr2prg => reextract
            case('refine3D');                    ptr2prg => refine3D
            case('refine3D_auto');               ptr2prg => refine3D_auto
            case('refine3D_nano');               ptr2prg => refine3D_nano
            case('replace_project_field');       ptr2prg => replace_project_field
            case('reproject');                   ptr2prg => reproject
            case('sample_classes');              ptr2prg => sample_classes
            case('scale');                       ptr2prg => scale
            case('select');                      ptr2prg => select_
            case('select_clusters');             ptr2prg => select_clusters
            case('selection');                   ptr2prg => selection
            case('sieve_cavgs');                 ptr2prg => sieve_cavgs
            case('simulate_atoms');              ptr2prg => simulate_atoms
            case('simulate_movie');              ptr2prg => simulate_movie
            case('simulate_noise');              ptr2prg => simulate_noise
            case('simulate_particles');          ptr2prg => simulate_particles
            case('split');                       ptr2prg => split_
            case('split_stack');                 ptr2prg => split_stack
            case('stack');                       ptr2prg => stack
            case('stackops');                    ptr2prg => stackops
            case('symaxis_search');              ptr2prg => symaxis_search
            case('symmetrize_map');              ptr2prg => symmetrize_map
            case('symmetry_test');               ptr2prg => symmetry_test
            case('atoms_rmsd');          ptr2prg => atoms_rmsd
            case('core_atoms_analysis'); ptr2prg => core_atoms_analysis
            case('tsegmaps_core_finder');      ptr2prg => tsegmaps_core_finder
            case('tseries_import');              ptr2prg => tseries_import
            case('import_trajectory');    ptr2prg => import_trajectory
            case('tseries_make_pickavg');        ptr2prg => tseries_make_pickavg
            case('trajectory_make_projavgs');    ptr2prg => trajectory_make_projavgs
            case('tseries_motion_correct');      ptr2prg => tseries_motion_correct
            case('trajectory_reconstruct3D');    ptr2prg => trajectory_reconstruct3D
            case('trajectory_swap_stack');       ptr2prg => trajectory_swap_stack
            case('track_particles');             ptr2prg => track_particles
            case('uniform_filter2D');            ptr2prg => uniform_filter2D
            case('uniform_filter3D');            ptr2prg => uniform_filter3D
            case('update_project');              ptr2prg => update_project
            case('vizoris');                     ptr2prg => vizoris
            case('volanalyze');                  ptr2prg => volanalyze
            case('volops');                      ptr2prg => volops
            case('write_classes');               ptr2prg => write_classes
            case('write_mic_filetab');           ptr2prg => write_mic_filetab
            case('zero_project_shifts');         ptr2prg => zero_project_shifts
            ! test programs
            case('test_sim_workflow');           ptr2prg => test_sim_workflow
            case DEFAULT
                ptr2prg => null()
        end select
    end subroutine get_prg_ptr

    subroutine list_simple_prgs_in_ui

        !====================================================================
        ! PROJECT MANAGEMENT
        !====================================================================
        write(logfhandle,'(A)') format_str('PROJECT MANAGEMENT:', C_UNDERLINED)
        write(logfhandle,'(A)') assign_optics_groups%name%to_char()
        write(logfhandle,'(A)') export_relion%name%to_char()
        write(logfhandle,'(A)') export_starproject%name%to_char()
        write(logfhandle,'(A)') extract_subproj%name%to_char()
        write(logfhandle,'(A)') import_boxes%name%to_char()
        write(logfhandle,'(A)') import_cavgs%name%to_char()
        write(logfhandle,'(A)') import_movies%name%to_char()
        write(logfhandle,'(A)') import_particles%name%to_char()
        write(logfhandle,'(A)') import_starproject%name%to_char()
        write(logfhandle,'(A)') merge_projects%name%to_char()
        write(logfhandle,'(A)') new_project%name%to_char()
        write(logfhandle,'(A)') print_project_field%name%to_char()
        write(logfhandle,'(A)') print_project_info%name%to_char()
        write(logfhandle,'(A)') prune_project%name%to_char()
        write(logfhandle,'(A)') replace_project_field%name%to_char()
        write(logfhandle,'(A)') selection%name%to_char()
        write(logfhandle,'(A)') update_project%name%to_char()
        write(logfhandle,'(A)') zero_project_shifts%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! PRE-PROCESSING
        !====================================================================
        write(logfhandle,'(A)') format_str('PRE-PROCESSING:', C_UNDERLINED)
        write(logfhandle,'(A)') preprocess%name%to_char()
        write(logfhandle,'(A)') extract%name%to_char()
        write(logfhandle,'(A)') reextract%name%to_char()
        write(logfhandle,'(A)') motion_correct%name%to_char()
        write(logfhandle,'(A)') gen_pspecs_and_thumbs%name%to_char()
        write(logfhandle,'(A)') ctf_estimate%name%to_char()
        write(logfhandle,'(A)') pick%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! CLUSTER2D WORKFLOWS
        !====================================================================
        write(logfhandle,'(A)') format_str('CLUSTER2D WORKFLOWS:', C_UNDERLINED)
        write(logfhandle,'(A)') abinitio2D%name%to_char()
        write(logfhandle,'(A)') cleanup2D%name%to_char()
        write(logfhandle,'(A)') cluster2D%name%to_char()
        write(logfhandle,'(A)') cluster2D_subsets%name%to_char()
        write(logfhandle,'(A)') cluster_cavgs%name%to_char()
        write(logfhandle,'(A)') cluster_cavgs_selection%name%to_char()
        write(logfhandle,'(A)') cluster_stack%name%to_char()
        write(logfhandle,'(A)') make_cavgs%name%to_char()
        write(logfhandle,'(A)') map_cavgs_selection%name%to_char()
        write(logfhandle,'(A)') match_cavgs%name%to_char()
        write(logfhandle,'(A)') match_stacks%name%to_char()
        write(logfhandle,'(A)') sample_classes%name%to_char()
        write(logfhandle,'(A)') select_clusters%name%to_char()
        write(logfhandle,'(A)') write_classes%name%to_char()
        write(logfhandle,'(A)') write_mic_filetab%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! AB INITIO 3D RECONSTRUCTION
        !====================================================================
        write(logfhandle,'(A)') format_str('AB INITIO 3D RECONSTRUCTION:', C_UNDERLINED)
        write(logfhandle,'(A)') abinitio3D%name%to_char()
        write(logfhandle,'(A)') abinitio3D_cavgs%name%to_char()
        write(logfhandle,'(A)') estimate_lpstages%name%to_char()
        write(logfhandle,'(A)') multivol_assign%name%to_char()
        write(logfhandle,'(A)') noisevol%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! REFINE3D
        !====================================================================
        write(logfhandle,'(A)') format_str('REFINE3D:', C_UNDERLINED)
        write(logfhandle,'(A)') refine3D%name%to_char()
        write(logfhandle,'(A)') refine3D_auto%name%to_char()
        write(logfhandle,'(A)') reconstruct3D%name%to_char()
        write(logfhandle,'(A)') postprocess%name%to_char()
        write(logfhandle,'(A)') automask%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! DENOISING
        !====================================================================
        write(logfhandle,'(A)') format_str('DENOISING:', C_UNDERLINED)
        write(logfhandle,'(A)') icm2D%name%to_char()
        write(logfhandle,'(A)') icm3D%name%to_char()
        write(logfhandle,'(A)') ppca_denoise%name%to_char()
        write(logfhandle,'(A)') ppca_denoise_classes%name%to_char()
        write(logfhandle,'(A)') ppca_volvar%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! FILTERING
        !====================================================================
        write(logfhandle,'(A)') format_str('FILTERING:', C_UNDERLINED)
        write(logfhandle,'(A)') filter%name%to_char()
        write(logfhandle,'(A)') uniform_filter2D%name%to_char()
        write(logfhandle,'(A)') uniform_filter3D%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! GENERAL IMAGE PROCESSING
        !====================================================================
        write(logfhandle,'(A)') format_str('GENERAL IMAGE PROCESSING:', C_UNDERLINED)
        write(logfhandle,'(A)') binarize%name%to_char()
        write(logfhandle,'(A)') convert%name%to_char()
        write(logfhandle,'(A)') ctf_phaseflip%name%to_char()
        write(logfhandle,'(A)') ctfops%name%to_char()
        write(logfhandle,'(A)') normalize_%name%to_char()
        write(logfhandle,'(A)') scale%name%to_char()
        write(logfhandle,'(A)') stack%name%to_char()
        write(logfhandle,'(A)') stackops%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! MASKING
        !====================================================================
        write(logfhandle,'(A)') format_str('MASKING:', C_UNDERLINED)
        write(logfhandle,'(A)') auto_spher_mask%name%to_char()
        write(logfhandle,'(A)') automask2D%name%to_char()
        write(logfhandle,'(A)') mask%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! ORIENTATION PROCESSING
        !====================================================================
        write(logfhandle,'(A)') format_str('ORIENTATION PROCESSING:', C_UNDERLINED)
        write(logfhandle,'(A)') make_oris%name%to_char()
        write(logfhandle,'(A)') orisops%name%to_char()
        write(logfhandle,'(A)') oristats%name%to_char()
        write(logfhandle,'(A)') vizoris%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! PARALLEL UTILITIES
        !====================================================================
        write(logfhandle,'(A)') format_str('PARALLEL UTILITIES:', C_UNDERLINED)
        write(logfhandle,'(A)') split_%name%to_char()
        write(logfhandle,'(A)') split_stack%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! PRINT INFO
        !====================================================================
        write(logfhandle,'(A)') format_str('PRINT INFO:', C_UNDERLINED)
        write(logfhandle,'(A)') info_image%name%to_char()
        write(logfhandle,'(A)') info_stktab%name%to_char()
        write(logfhandle,'(A)') print_dose_weights%name%to_char()
        write(logfhandle,'(A)') print_fsc%name%to_char()
        write(logfhandle,'(A)') print_magic_boxes%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! RESOLUTION ESTIMATION
        !====================================================================
        write(logfhandle,'(A)') format_str('RESOLUTION ESTIMATION:', C_UNDERLINED)
        write(logfhandle,'(A)') fsc%name%to_char()
        write(logfhandle,'(A)') clin_fsc%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! SIMULATION
        !====================================================================
        write(logfhandle,'(A)') format_str('SIMULATION:', C_UNDERLINED)
        write(logfhandle,'(A)') pdb2mrc%name%to_char()
        write(logfhandle,'(A)') simulate_movie%name%to_char()
        write(logfhandle,'(A)') simulate_noise%name%to_char()
        write(logfhandle,'(A)') simulate_particles%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! STREAM VALIDATION
        !====================================================================
        write(logfhandle,'(A)') format_str('STREAM VALIDATION:', C_UNDERLINED)
        write(logfhandle,'(A)') mini_stream%name%to_char()
        write(logfhandle,'(A)') check_refpick%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! SYMMETRY
        !====================================================================
        write(logfhandle,'(A)') format_str('SYMMETRY:', C_UNDERLINED)
        write(logfhandle,'(A)') symaxis_search%name%to_char()
        write(logfhandle,'(A)') symmetrize_map%name%to_char()
        write(logfhandle,'(A)') symmetry_test%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! VOLUME DOCKING
        !====================================================================
        write(logfhandle,'(A)') format_str('VOLUME DOCKING:', C_UNDERLINED)
        write(logfhandle,'(A)') dock_volpair%name%to_char()
        write(logfhandle,'(A)') volanalyze%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! VOLUME PROCESSING
        !====================================================================
        write(logfhandle,'(A)') format_str('VOLUME PROCESSING:', C_UNDERLINED)
        write(logfhandle,'(A)') center%name%to_char()
        write(logfhandle,'(A)') reproject%name%to_char()
        write(logfhandle,'(A)') volops%name%to_char()
        write(logfhandle,'(A)') ''
        !====================================================================
        ! MODEL ANALYSIS
        !====================================================================
        write(logfhandle,'(A)') format_str('MODEL ANALYSIS PROGRAMS:', C_UNDERLINED)
        write(logfhandle,'(A)') model_validation%name%to_char()
        write(logfhandle,'(A)') ''
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
        write(logfhandle,'(A)') abinitio2D_stream%name%to_char()
        write(logfhandle,'(A)') assign_optics%name%to_char()
        write(logfhandle,'(A)') cluster2D_stream%name%to_char()
        write(logfhandle,'(A)') gen_pickrefs%name%to_char()
        write(logfhandle,'(A)') pick_extract%name%to_char()
        write(logfhandle,'(A)') preproc%name%to_char()
        write(logfhandle,'(A)') sieve_cavgs%name%to_char()
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
        integer                   :: iprg
        ! JSON init
        call json%initialize()
        ! create object of program entries
        call json%create_object(all_programs, 'SIMPLE_UI')
        do iprg = 1, n_prg_ptrs
            call create_program_entry(json, prg_ptr_array(iprg)%ptr2prg, all_programs)
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
        integer :: iprg
        ! JSON init
        call json%initialize()
        ! create array of program entries
        call json%create_array(all_programs, 'SIMPLE User Interface')
        do iprg = 1, n_prg_ptrs
            call create_program_entry(json, prg_ptr_array(iprg)%ptr2prg, all_programs)
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
