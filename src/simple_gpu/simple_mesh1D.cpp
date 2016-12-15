/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 25th of February 2016.13:55pm
 *
 * Name:
 * simple_mesh_1D.cpp - class to get recover the mesh for the 3D data
 *                      Structure for the device
 *
 * Description:
 * classe to set and get the mesh_1D as generic both both CPU/GPU
 *******************************************************************************
 */

#if defined (CUDA) /*preprossing for the CUDA environment */

#include <iostream>
#include "polarft_gpu.h"
#include "global_polarft.h"

/*contructors*/
//template <typename S>
mesh_1D::mesh_1D(polar_corr_calc *s_polar_in) {
  s_polar_o = s_polar_in;
}
//template <typename S>
mesh_1D::mesh_1D(polar_corr_calc *s_polar_in,/* overloading the constructors*/
                 int n_particles_in, int n_rotations_in, int n_rezstions_in)
{
  s_polar_o = s_polar_in;
  n_particles_o = n_particles_in;
  n_rotations_o = n_rotations_in;
  n_rezstions_o = n_rezstions_in;
}
//setters
void mesh_1D::set_nx(polar_corr_calc *s_polar) {s_polar_o->nx = s_polar->nx; }
void mesh_1D::set_ny(polar_corr_calc *s_polar) {s_polar_o->ny = s_polar->nz; }
void mesh_1D::set_nz(polar_corr_calc *s_polar) {s_polar_o->nz = s_polar->nz; }
//getters
int mesh_1D::get_mesh1D_nx()    {return s_polar_o->nx;                       }
int mesh_1D::get_mesh1D_ny()    {return s_polar_o->ny;                       }
int mesh_1D::get_mesh1D_nz()    {return s_polar_o->nz;                       }
int mesh_1D::get_mesh1D_chunk() {return n_rotations_o * n_rezstions_o;       }
int mesh_1D::get_mesh1D_threadsPerBlock(){return s_polar_o->threadsPerBlock; }
int mesh_1D::get_mesh1D_blocksPerGrid()  {
  return get_mesh1D_chunk()/(float)get_mesh1D_threadsPerBlock()
    + ( get_mesh1D_chunk()%get_mesh1D_threadsPerBlock()!=0 );
}
int mesh_1D::get_mesh1D_size_pB() {
  return get_mesh1D_blocksPerGrid()*sizeof(float);
}
// the destructor
mesh_1D::~mesh_1D() {
  std::cout << "Object mesh_1D has been destroyed" << std::endl;
}

#endif /* CUDA */
