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

#include "waveletsD.h"

// polarft_size global variables
string name_o;
vector<double> lp1_o; vector<double> hp1_o;
vector<double> lp2_o; vector<double> hp2_o;

//global objects intanciation for polarft
wavelets_coef *p_wcoef = new wavelets_coef(name_o,
                                           lp1_o, hp1_o, lp2_o, hp2_o);
