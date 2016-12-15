/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common;

import java.io.File;

import javax.swing.filechooser.FileFilter;

import org.openmolgrid.util.CUtil;

public class DialogExtensionFilter extends FileFilter
{
   private String ext="";
   private String desc="";

   public DialogExtensionFilter(String extension, String description)
   {
      ext=extension;
      desc=description;
   }

   public boolean accept(File f)
   {
      if (f != null)
      {
         if (f.isDirectory())
         {
            return true;
         }
         String extension = CUtil.getExtension(f);
         if (extension != null && ext.equals(extension))
         {
            return true;
         }
      }
      return false;
   }

   public String getDescription()
   {
      return desc;
   }
}
