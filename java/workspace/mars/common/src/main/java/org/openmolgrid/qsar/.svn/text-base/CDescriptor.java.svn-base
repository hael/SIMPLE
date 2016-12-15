/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import java.io.Serializable;

public class CDescriptor implements Serializable
{
   private String id;
   private int pos=0;
   private double value;

   public CDescriptor(String id, double value)
   {
      this.id=id;
      this.value=value;
   }

   public CDescriptor(String id, int pos, double value)
   {
      this.id=id;
      this.pos=pos;
      this.value=value;
   }

   public CDescriptor(CDescriptor d)
   {
      id=d.id;
      pos=d.pos;
      value=d.value;
   }

   public String getId()
   {
      return id;
   }

   public int getPos()
   {
      return pos;
   }

   public double getValue()
   {
      return value;
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
