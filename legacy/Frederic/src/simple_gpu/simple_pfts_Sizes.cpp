/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 09th of February 2016.13:55pm
 *
 * Name:
 * simple_pfts_Sizes.cpp - class to get recover the sizes of the pfts A and B
 *
 * Description:
 * classe to set and get the pfts_Sizes as generic both both CPU/GPU
 *******************************************************************************
 */

#if defined (CUDA) /*preprossing for the CUDA environment */

#include <iostream>
#include "polarft_gpu.h"
#include "global_polarft.h"

pfts_Sizes::pfts_Sizes(int npart_in, int nrot_in, int nk_in) {
  npart_o = npart_in;
  nrot_o = nrot_in;
  nk_o = nk_in;
}
//setters
void pfts_Sizes::set_npart(int npart) {npart_o = npart;}
void pfts_Sizes::set_nrot(int nrot)   { nrot_o = nrot; }
void pfts_Sizes::set_nk(int nk)       {   nk_o = nk;   }
//getters
int pfts_Sizes::get_npart()    {return npart_o;                              }
int pfts_Sizes::get_nrot()     {return nrot_o;                               }
int pfts_Sizes::get_nk()       {return nk_o;                                 }
int pfts_Sizes::get_nhlfrot()  {return nrot_o / 2;                           }
int pfts_Sizes::get_n2rot()    {return nrot_o * 2;                           }
int pfts_Sizes::get_n2part()   {return npart_o * 2;                          }
int pfts_Sizes::get_npts()     {return npart_o * nrot_o * nk_o;              }
int pfts_Sizes::get_npts_A()   {return npart_o * (get_nhlfrot()) * nk_o;     }
int pfts_Sizes::get_npts_B()   {return (get_n2part()) * (get_n2rot()) * nk_o;}
int pfts_Sizes::get_npts_cr3D(){return (get_npart()) * (get_npart()) *
                                       (get_nrot());                         }
// the destructor
pfts_Sizes::~pfts_Sizes() {
  std::cout << "Object pfts_Sizes has been destroyed" << std::endl;
}

#endif /* CUDA */
