/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import java.util.*;
import java.io.*;


public class CPropertyList implements Serializable, CAttribute
{
   private Vector<CProperty> properties;

   public CPropertyList()
   {
      properties=new Vector<CProperty>();
   }

   public CPropertyList(CPropertyList p)
   {
      properties=new Vector<CProperty>();
      for (int i=0; i < p.properties.size(); i++)
      {
         properties.add(new CProperty((CProperty)p.properties.get(i)));
      }
   }

   public void addProperty(CProperty p)
   {
      properties.add(p);
   }

   public void addProperty(String id, double value)
   {
      properties.add(new CProperty(id, value));
   }

   public void addProperty(String id, double value, String modelId)
   {
      properties.add(new CProperty(id, value, modelId));
   }

   public CProperty getProperty(int idx)
   {
      return properties.get(idx);
   }
   
   public CProperty getPropertyById(String id)
   {
      for (int i=0; i < properties.size(); i++)
      {
         CProperty p = properties.get(i);
         if (p.getId().equals(id)) return p;
      }
      return null;
   }

   public int getNumberOfProperties()
   {
      return properties.size();
   }

	public CAttribute copy() {
		return new CPropertyList(this);
	}

   public void writeCDR(PrintWriter pw)
   {
	  pw.println("    <propertyList>");
      for (int i=0; i < properties.size(); i++)
      {
         CProperty p=properties.get(i);
         if (p.getModelId().equals(""))
            pw.println("      <property id=\""+p.getId()+"\" value=\""+p.getValue()+"\"/>");
         else
            pw.println("      <property id=\""+p.getId()+"\" modelId=\""+p.getModelId()+"\" value=\""+p.getValue()+"\"/>");
      }
      pw.println("    </propertyList>");
   }
}
