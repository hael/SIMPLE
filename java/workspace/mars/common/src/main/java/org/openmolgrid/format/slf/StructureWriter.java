/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.format.slf;

import java.io.*;

import org.openmolgrid.model.CStructure;
import org.openmolgrid.qsar.CDescriptorMetaMap;
import org.openmolgrid.qsar.CPropertyMetaMap;
import org.openmolgrid.qsar.CSource;

/**
 * This class is intended for creating structure list files.<p>
 * 
 * This class is also used as a base class for the FragmentWriter.
 *
 * @author Andre Lomaka
 * @author Sulev Sild
 */
public class StructureWriter
{
   protected PrintWriter outWriter = null;
   private boolean haveMetaData = false;

   public StructureWriter(String outFileName)
   {
      this(new File(outFileName));
   }

   public StructureWriter(File outFile)
   {
      try
      {
         outWriter = new PrintWriter(new FileWriter(outFile));
      } catch (IOException e)
      {
         throw new RuntimeException(e.toString());
      }
      writeHeader();
   }

   protected void writeHeader()
   {
      outWriter.println("<?xml version=\"1.0\"?>");
      outWriter.println("<structureList>");
   }

   protected void writeFooter()
   {
      outWriter.println("</structureList>");
   }

   private void startMetaDataElementIfNeeded()
   {
      if (haveMetaData == false)
      {
         outWriter.println("  <metadata>");
         haveMetaData = true;
      }
   }

   private void endMetaDataElementIfNeeded()
   {
      if (haveMetaData == true)
      {
         outWriter.println("  </metadata>");
         haveMetaData = false;
      }
   }

   /**
    * Writes source information about the structure list. For example, 
    * a software that produced this list.
    *
    * @param source information
    */
   public void writeSource(CSource source)
   {
      source.writeCDR(outWriter);
   }

   /**
    * Writes descriptor metadata.
    *
    * @param metaList is a container with metadata
    */
   public void writeDescriptorMeta(CDescriptorMetaMap metaList)
   {
      startMetaDataElementIfNeeded();
      metaList.writeCDR(outWriter);
   }

   /**
    * Writes property metadata.
    *
    * @param metaList is a container with metadata
    */
   public void writePropertyMeta(CPropertyMetaMap metaList)
   {
      startMetaDataElementIfNeeded();
      metaList.writeCDR(outWriter);
   }

   public void addStructure(CStructure s)
   {
      endMetaDataElementIfNeeded();
      if (outWriter == null) return;
      s.writeCDR(outWriter);
   }

   public void close()
   {
      if (outWriter != null)
      {
         endMetaDataElementIfNeeded();
         writeFooter();
         if (outWriter.checkError())
         {
            throw new RuntimeException("Error while writing to StructureWriter"); 
         }
         outWriter.close();
         outWriter = null;
      }
   }
}
