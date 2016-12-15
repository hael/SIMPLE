/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common;

// TODO: include structure ID & name in SDF file

import java.io.File;

import javax.swing.JFileChooser;

import org.openmolgrid.client.common.table.StructureList;
import org.openmolgrid.client.common.table.StructureListTable;
import org.openmolgrid.format.sdf.SDFWriter;


/**
 * StructureListExporter class is used to export structures from
 * StructureListTable. Currently export to SDF format is implemented.
 *
 * @author Sulev Sild
 */
public class StructureListExporter 
{
	private StructureListTable table;
	
	/**
	 * Creates a new StructureListExporter instance.
	 *
	 * @param sourceTable is a table containing structures to be exported
	 */
	public StructureListExporter(StructureListTable sourceTable)
	{
		table = sourceTable;
	}

	/**
	 * Handle export of structures. Asks the name of an output file and
	 * exports structures.
	 */
	public void handleExport()
	{
		// dialog that asks the output file name
		JFileChooser chooser = new JFileChooser();
		int ret = chooser.showSaveDialog(table);

		if (ret == JFileChooser.APPROVE_OPTION)
		{
			File file = chooser.getSelectedFile();
			doExport(file);
		}
	}

	/**
	 * Exports structures to a given file.
	 *
	 * @param file is a File object with destination file
	 */
	private void doExport(File file) {
		SDFWriter writer = new SDFWriter(file);

		StructureList sl = table.getStructureList();
		for (int i = 0; i < sl.size(); ++i) {
			writer.addStructure(sl.get(i));
		}
		writer.close();
	}
}
