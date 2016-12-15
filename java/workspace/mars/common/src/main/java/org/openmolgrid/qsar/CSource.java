/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import java.io.*;

public class CSource
{
   private String softwareName="N/A";
   private String version="N/A";
   private String architecture=System.getProperty("os.arch");
   private String operatingSystem=System.getProperty("os.name");

   public void setSoftwareName(String arg)
   {
      softwareName=arg;
   }

   public void setVersion(String arg)
   {
      version=arg;
   }

   public void setArchitecture(String arg)
   {
      architecture=arg;
   }

   public void setOperatingSystem(String arg)
   {
      operatingSystem=arg;
   }

   public void writeCDR(PrintWriter pw)
   {
      pw.println("  <source>");
      if (softwareName != null)
         pw.println("    <softwareName>"+softwareName+"</softwareName>");
      if (version != null)
         pw.println("    <version>"+version+"</version>");
      if (architecture != null)
         pw.println("    <architecture>"+architecture+"</architecture>");
      if (operatingSystem != null)
         pw.println("    <operatingSystem>"+operatingSystem+"</operatingSystem>");
      pw.println("  </source>");
   }
}
