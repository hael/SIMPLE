/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.util;

import java.io.*;

/**
 * Various utility functions for string and file manipulations.
 *
 * @author Andre Lomaka
 * @author Sulev Sild
 */
public class CUtil
{
   public static String file2String(File f)
   {
      FileReader in;
      try
      {
         in=new FileReader(f);
      } catch (FileNotFoundException e)
      {
         System.err.println("File " + f.toString() + " not found");
         return "";
      }
      return convert(f, in);
   }

   public static void string2File(String s, File f)
   {
      PrintWriter out=null;
      try
      {
         out = new PrintWriter(new BufferedWriter(new FileWriter(f)));
         out.print(s);
      } catch (IOException e)
      {
         throw new RuntimeException("Unable to write: "+f);
      } finally {
    	  if (out != null) out.close();
      }
   }

   private static String convert(File f, Reader in)
   {
      int c;
      ByteArrayOutputStream bo = new ByteArrayOutputStream((int)f.length());
      try
      {
         while ((c = in.read()) != -1)
            bo.write(c);
         in.close();
      } catch (IOException e)
      {
         System.err.println(e.toString());
         return "";
      }
      try
      {
         return bo.toString("UTF8");
      } catch (UnsupportedEncodingException e)
      {
         System.err.println(e.toString());
         return "";
      }
   }

   /**
    * Get the extension of a file
    */
   public static String getExtension(File f)
   {
      String s = f.getName();
      return getExtension(s);
   }

   /**
    * Get the extension of a file name
    */
   public static String getExtension(String fn)
   {
      String ext=null;
      int i = fn.lastIndexOf('.');
      if (i > 0 && i < fn.length()-1)
      {
         ext = fn.substring(i+1).toLowerCase();
      }
      return ext;
   }

   /**
    * Get path of a file name
    */
   public static String getPath(String fn)
   {
      String path = null;
      int i = fn.lastIndexOf(File.separatorChar);
      if (i > 0 && i < fn.length()-1)
      {
         path = fn.substring(0, i);
      }
      return path;
   }

   /**
    * Returns escaped string that is safe for writing to XML format.
    *
    * @param str is a string to be escaped
    * @return a String with escaped XML characters
    */
   static public String escapeXMLString(String str)
   {
		for (int i=0; i<str.length(); ++i) {
			if ("<>&\"'".indexOf(str.charAt(i)) != -1) {
				return "<![CDATA[" + str + "]]>";
			}
		}
		return str;   
	}

   /**
    * Deletes a given file or directory (including its sub-directories).
    * The result will be the same as with "rm -rf target" command under
    * Unix systems. 
    *
    * @param target is the file or directory to be deleted
    * @return true when delete operation was successful, false otherwise
    */
   static public boolean deleteAll(File target)
   {

      // If target directory is read only then symbolic links contained within
      // can't be removed. In this case symbolic links to other directories
      // will be followed and files outside the target might be removed.
      if (!target.canWrite())
      {
         return false;
      }

      // This will work for normal files, empty directories, and symbolic links.
      // Without this symbolic links to other directories will be followed and
      // files outside the target directory might be removed.
      if (target.delete())
      {
         return true;
      }

      // Here we should have a non-empty directory. Empty it by recursion.
      if(target.isDirectory())
      {
         File[] files = target.listFiles();
         for(int i=0; i<files.length; ++i)
         {
            deleteAll(files[i]);
         }
      } else
      {
         return false;
      }

      // Finally there will be left only an empty directory
      return target.delete();
   }
}
