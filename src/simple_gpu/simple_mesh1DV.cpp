/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 24th of April 2016.14:55pm
 *
 * Name:
 * simple_mesh_1DV.cpp - class to get recover the mesh for the 1D cartesian data
 *                      Structure for the device
 *
 * Description:
 * classe to set and get the mesh_1DV as generic both both CPU/GPU
 *******************************************************************************
 */

#if defined (CUDA) /*preprossing for the CUDA environment */

#include <iostream>
#include "polarft_gpu.h"
#include "global_polarft.h"

/*contructors*/
//template <typename S>
mesh_1DV::mesh_1DV(polar_corr_calc *s_carte_in) {
  s_carte_o = s_carte_in;
}
//template <typename S>
mesh_1DV::mesh_1DV(polar_corr_calc *s_carte_in,/* overloading the constructors*/
                 int vx_in, int vy_in, int vz_in)
{
  s_carte_o = s_carte_in;
  vx_o = vx_in;
  vy_o = vy_in;
  vz_o = vz_in;
}
//setters
void mesh_1DV::set_nx(polar_corr_calc *s_carte) {s_carte_o->nx = s_carte->nx; }
void mesh_1DV::set_ny(polar_corr_calc *s_carte) {s_carte_o->ny = s_carte->nz; }
void mesh_1DV::set_nz(polar_corr_calc *s_carte) {s_carte_o->nz = s_carte->nz; }
//getters
int mesh_1DV::get_mesh1DV_nx()    {return s_carte_o->nx;                       }
int mesh_1DV::get_mesh1DV_ny()    {return s_carte_o->ny;                       }
int mesh_1DV::get_mesh1DV_nz()    {return s_carte_o->nz;                       }
int mesh_1DV::get_mesh1DV_chunk() {return vx_o * vy_o * vz_o;                  }
int mesh_1DV::get_mesh1DV_threadsPerBlock(){return s_carte_o->threadsPerBlock; }
int mesh_1DV::get_mesh1DV_blocksPerGrid()  {
  return get_mesh1DV_chunk()/(float)get_mesh1DV_threadsPerBlock()
     + ( get_mesh1DV_chunk()%get_mesh1DV_threadsPerBlock()!=0 );
}
int mesh_1DV::get_mesh1DV_size_pB() {
  return get_mesh1DV_blocksPerGrid()*sizeof(float);
}
// the destructor
mesh_1DV::~mesh_1DV() {
  std::cout << "Object mesh_1DV has been destroyed" << std::endl;
}

#endif /* CUDA */

