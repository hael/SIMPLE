/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import org.openmolgrid.model.CStructure;

public interface StructureConsumer
{
   void startHandling();
   void endHandling();
   
   /**
    * Handler that is called for each processed structure. 
    * 
    * @param s - structure to be consumed
    * @throws CancellationException to abort processing
    */
   void consumeStructure(CStructure s);
}
