/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 09th of February 2016. 15:13pm
 *
 * Name:
 * global_polarft.h - basic definitions used in all modules.
 *
 * Description:
 * header file for the Global_polarft.cpp file to set the 
 * global variable within the cpp files and routines.
 * class declarion for the header file global_polarft.h
 *******************************************************************************
 */

#if defined (CUDA) /*preprossing for the CUDA environment */

//Inlcude file for the global objects
#include "polarft_gpu.h"

//global variables for the object
extern int npart_o, nrot_o, nk_o;
extern polar_corr_calc *s_polar_o;
extern polar_corr_calc *s_carte_o;
extern int n_particles_o;
extern int n_rotations_o;
extern int n_rezstions_o;
extern int vx_o,vy_o,vz_o;

//global objects declared for polarft calculator
extern pfts_Sizes *p_pfts_Sizes;
extern mesh_3D *p_mesh_3D;
extern mesh_1D *p_mesh_1D;
extern mesh_1DV *p_mesh_1DV;
extern img_2D_cart_Sizes *p_img_2D_cart_Sizes;

#endif /* CUDA */
