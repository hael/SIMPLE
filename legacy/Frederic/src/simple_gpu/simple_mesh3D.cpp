/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 25th of February 2016.13:55pm
 *
 * Name:
 * simple_mesh_3D.cpp - class to get recover the mesh for the 3D data
 *                      Structure for the device
 *
 * Description:
 * classe to set and get the mesh_3D as generic both both CPU/GPU
 *******************************************************************************
 */

#if defined (CUDA) /*preprossing for the CUDA environment */

#include <iostream>
#include "polarft_gpu.h"
#include "global_polarft.h"

/*contructors*/
mesh_3D::mesh_3D(polar_corr_calc *s_polar_in) {
  s_polar_o = s_polar_in;
}
mesh_3D::mesh_3D(polar_corr_calc *s_polar_in,/* overloading the constructors*/
                 int n_particles_in, int n_rotations_in, int n_rezstions_in)
{
  s_polar_o = s_polar_in;
  n_particles_o = n_particles_in;
  n_rotations_o = n_rotations_in;
  n_rezstions_o = n_rezstions_in;
}
//setters
void mesh_3D::set_nx(polar_corr_calc *s_polar) {s_polar_o->nx = s_polar->nx; }
void mesh_3D::set_ny(polar_corr_calc *s_polar) {s_polar_o->ny = s_polar->nz; }
void mesh_3D::set_nz(polar_corr_calc *s_polar) {s_polar_o->nz = s_polar->nz; }
//getters
int mesh_3D::get_mesh3D_nx()    {return s_polar_o->nx;                       }
int mesh_3D::get_mesh3D_ny()    {return s_polar_o->ny;                       }
int mesh_3D::get_mesh3D_nz()    {return s_polar_o->nz;                       }

int mesh_3D::get_mesh3D_gridx() {return  n_particles_o/(float)get_mesh3D_nx()
                                     + ( n_particles_o%get_mesh3D_nx()!=0 ); }
int mesh_3D::get_mesh3D_gridy() {return  n_rotations_o/(float)get_mesh3D_ny()
                                     + ( n_rotations_o%get_mesh3D_ny()!=0 ); }
int mesh_3D::get_mesh3D_gridz() {return  n_rezstions_o/(float)get_mesh3D_nz()
                                     + ( n_rezstions_o%get_mesh3D_nz()!=0 ); }
// the destructor
mesh_3D::~mesh_3D() {
  /* TODO: insert the debug datastrure for object destruction */
  //std::cout << "Object mesh_3D has been destroyed" << std::endl;
}

#endif /* CUDA */
