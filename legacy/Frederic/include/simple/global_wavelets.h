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

//Inlcude file for the global objects
#include "waveletsD.h"

//global variables for the object

extern string name_o;
extern vector<double> lp1_o; extern vector<double> hp1_o;
extern vector<double> lp2_o; extern vector<double> hp2_o;

//global objects declared for polarft calculator
extern wavelets_coef *p_wcoef;
