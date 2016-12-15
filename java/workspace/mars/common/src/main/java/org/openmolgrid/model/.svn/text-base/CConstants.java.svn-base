/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.model;

public enum CConstants
{
   R,
   H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, 
   K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, 
   Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, 
   I, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, 
   Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, 
   Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, 
   Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg; 

   public static String getStringByElement(int el) {
      CConstants e = values()[el];
      return e.name();
   }
                              
   public static int getElementByString(String el) {
      for (CConstants e: CConstants.values()) {
         if (e.name().equalsIgnoreCase(el)) {
            return e.ordinal();
         }
      }
      throw new IllegalArgumentException("Unknown element: "+el);
   }
}
