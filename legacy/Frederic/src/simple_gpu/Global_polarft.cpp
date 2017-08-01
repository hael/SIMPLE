/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 09th of February 2016. 15:13am
 *
 * Name:
 * Global_polarft.cpp - basic definitions used in all modules.
 *
 * Description:
 * file to englodbe all the global variables and
 * class declarion for the header file global_polarft.h
 *******************************************************************************
 */

#if defined (CUDA) /*preprossing for the CUDA environment */

#include "polarft_gpu.h"

// polarft_size global variables
int npart_o, nrot_o, nk_o;
polar_corr_calc *s_polar_o;
polar_corr_calc *s_carte_o;
int n_particles_o;
int n_rotations_o;
int n_rezstions_o;
int vx_o,vy_o,vz_o;

//global objects intanciation for polarft
pfts_Sizes *p_pfts_Sizes = new pfts_Sizes(npart_o,nrot_o,nk_o);
   mesh_3D *p_mesh_3D    = new mesh_3D(s_polar_o,
                                       n_particles_o,
                                       n_rotations_o, n_rezstions_o);
   mesh_1D *p_mesh_1D    = new mesh_1D(s_polar_o,
                                       n_particles_o,
                                       n_rotations_o, n_rezstions_o);
 mesh_1DV *p_mesh_1DV    = new mesh_1DV(s_carte_o, vx_o, vy_o, vz_o);
img_2D_cart_Sizes *p_2D_cart_Sizes = new img_2D_cart_Sizes(vx_o, vy_o, vz_o);
#endif /* CUDA */
