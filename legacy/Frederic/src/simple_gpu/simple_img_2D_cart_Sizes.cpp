/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 09th of February 2016.13:55pm
 *
 * Name:
 * simple_img_2D_cart_Sizes.cpp - class to get recover the sizes of the pfts A and B
 *
 * Description:
 * classe to set and get the 2D_cart_Sizes as generic both both CPU/GPU
 *******************************************************************************
 */

#if defined (CUDA) /*preprossing for the CUDA environment */

#include <iostream>
#include "polarft_gpu.h"
#include "global_polarft.h"

img_2D_cart_Sizes::img_2D_cart_Sizes(int vx_in, int vy_in, int vz_in) {
  vx_o = vx_in;
  vy_o = vy_in;
  vz_o = vz_in;
}
//setters
void img_2D_cart_Sizes::set_2D_vx(int vx) {vx_o = vx; }
void img_2D_cart_Sizes::set_2D_vy(int vy) {vy_o = vy; }
void img_2D_cart_Sizes::set_2D_vz(int vz) {vz_o = vz; }
//getters
int img_2D_cart_Sizes::get_2D_vx() {return vx_o;                               }
int img_2D_cart_Sizes::get_2D_vy() {return vy_o;                               }
int img_2D_cart_Sizes::get_2D_vz() {return vz_o;                               }
int img_2D_cart_Sizes::get_npts()  {return get_2D_vx()*get_2D_vy()*get_2D_vz();}
int img_2D_cart_Sizes::get_npts_A(){return vx_o * (get_2D_vy()) * vz_o;        }
int img_2D_cart_Sizes::get_npts_B(){return vx_o * (get_2D_vy()) * vz_o;        }
// the destructor
img_2D_cart_Sizes::~img_2D_cart_Sizes() {
  std::cout << "Object img_2D_cart_Sizes has been destroyed" << std::endl;
}

#endif /* CUDA */
