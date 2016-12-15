/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 25th of February 2016.13:55pm
 *
 * Name:
 * wavelet_coef.cpp - class to get the wavelet coefficients
 *
 * Description:
 * classe to set and get the wavelet coefficients
 *******************************************************************************
 */


#include <iostream>
#include "global_wavelets.h"

/*contructors*/
wavelets_coef::wavelets_coef(string name_in,
                             vector<double> lp1_in, vector<double> hp1_in,
                             vector<double> lp2_in, vector<double> hp2_in)
{
  name_o = name_in;
  lp1_o = lp1_in; hp1_o = hp1_in;
  lp2_o = lp2_in; hp2_o = hp2_in;
}
//setters
void wavelets_coef::set_wavelets_name(string name)       {name_o = name;}
void wavelets_coef::set_wavelets_lp1(vector<double> lp1) {lp1_o = lp1;  }
void wavelets_coef::set_wavelets_hp1(vector<double> hp1) {hp1_o = hp1;  }
void wavelets_coef::set_wavelets_lp2(vector<double> lp2) {lp2_o = lp2;  }
void wavelets_coef::set_wavelets_hp2(vector<double> hp2) {hp2_o = hp2;  }
//getters
string wavelets_coef::get_wavelets_name() {name_o;}
vector<double> wavelets_coef::get_wavelets_lp1() {return lp1_o;}
vector<double> wavelets_coef::get_wavelets_hp1() {return hp1_o;}
vector<double> wavelets_coef::get_wavelets_lp2() {return lp2_o;}
vector<double> wavelets_coef::get_wavelets_hp2() {return hp2_o;}
// the destructor
wavelets_coef::~wavelets_coef() {
  /* TODO: insert the debug datastrure for object destruction */
  std::cout << "Object wavelets_coef has been destroyed" << std::endl;
}
