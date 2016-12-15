/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import org.openmolgrid.format.slf.StructureReader;
import org.openmolgrid.format.slf.StructureWriter;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.ChemLibException;

public class StructureSplitter implements StructureConsumer {

	private StructureReader sr;
	private StructureWriter sw;
	private int fileIndex;
	private String basename;
	private int batchSize;
	private int batchIndex;

	public StructureSplitter(String input, int batchSize, String basename) {
		this.sr = new StructureReader(input, this);
		this.batchSize = batchSize;
		this.batchIndex = 0;
		this.basename = basename;
		this.fileIndex = 0;
		mkWriter();
	}

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String input = args[0];
		int batchSize = Integer.parseInt(args[1]);
		String basename = args[2];

		new StructureSplitter(input, batchSize, basename).split();
	}

	private void split() throws ChemLibException {
		sr.provideStructures();
	}
	
	private void mkWriter() {
		fileIndex++;
		sw = new StructureWriter(String.format("%s_%05d.slf", basename, fileIndex));
	}

	public void consumeStructure(CStructure s) {
		sw.addStructure(s);
		batchIndex++;
		
		if (batchIndex == batchSize) {
			batchIndex = 0;
			sw.close();
			mkWriter();
		}
	}

	public void endHandling() {
		sw.close();
	}

	public void startHandling() {
	}

}
