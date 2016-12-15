/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.util;

import java.text.*;
import java.util.*;

public class CFormatter
{
   public static String formatStringLeft(String s, int le)
   {
      StringBuilder sb = new StringBuilder(s.trim());
      if (s.length() > le)
         return s.substring(0, le);
      int l = le - s.length();
      for (int f=0; f < l; f++)
         sb.append(" ");
      return sb.toString();
   }

   public static String formatStringRight(String s, int le)
   {
      StringBuilder fs = new StringBuilder();
      s = s.trim();
      if (s.length() > le)
         return s.substring(s.length()-le, s.length());
      int l = le - s.length();
      for (int f=0; f < l; f++)
         fs.append(" ");
      fs.append(s);
      return fs.toString();
   }

   public static String formatDouble(double fl, int le, String fmt)
   {
      String s = "";
      StringBuilder fs = new StringBuilder();
      NumberFormat nf = NumberFormat.getNumberInstance(Locale.ENGLISH);
      DecimalFormat df = (DecimalFormat)nf;
      df.applyPattern(fmt);
      s = df.format(fl);
      if (le == 0) return s;
      int l = le - s.length();
      for (int f=0; f < l; f++)
         fs.append(" ");
      fs.append(s);
      return fs.toString();
   }
}
