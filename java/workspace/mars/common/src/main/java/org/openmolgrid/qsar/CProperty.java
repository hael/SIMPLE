/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import java.io.Serializable;

public class CProperty implements Serializable
{
   private String id;
   private String modelId="";
   private double value;

   public CProperty(String id, double value)
   {
      this.id=id;
      this.value=value;
   }

   public CProperty(String id, double value, String modelId)
   {
      this.id=id;
      this.value=value;
      this.modelId=modelId;
   }

   public CProperty(CProperty d)
   {
      id=d.id;
      value=d.value;
      modelId=d.modelId;
   }

   public String getId()
   {
      return id;
   }

   public double getValue()
   {
      return value;
   }

   public String getModelId()
   {
      return modelId;
   }

   public void setId(String id)
   {
      this.id=id;
   }

   public void setValue(double value)
   {
      this.value=value;
   }
}
