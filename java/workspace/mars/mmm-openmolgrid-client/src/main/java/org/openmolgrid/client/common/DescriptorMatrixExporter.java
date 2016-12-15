/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import org.openmolgrid.client.common.table.StructureList;
import org.openmolgrid.client.common.table.StructureListTable;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.qsar.CDescriptor;
import org.openmolgrid.qsar.CDescriptorList;
import org.openmolgrid.qsar.CDescriptorMetaMap;

/**
 * DescriptorMatrixExporter
 *
 * @author Sulev Sild
 */
public class DescriptorMatrixExporter {
	private StructureListTable table;
	private CDescriptorMetaMap descMeta;

	private BufferedWriter writer;
	private ArrayList<String> currentRow;
	private HashMap<String, Integer> colMap; // maps desc id to corresponding column number
	
	public DescriptorMatrixExporter(StructureListTable sourceTable) {
		table = sourceTable;
		descMeta = table.getStructureList().getDescriptorMetaData();
	}

	/**
	 * Shows file selection dialog and exports calculated descripors to
	 * the selected file.
	 */
	public void handleExport() {
		// dialog for asking for the destination file
		JFileChooser chooser = new JFileChooser();
		int ret = chooser.showSaveDialog(table);
		try {
			if (ret == JFileChooser.APPROVE_OPTION) {
				File file = chooser.getSelectedFile();
				doExport(file);
			}
		} catch (IOException e) {
			JOptionPane.showMessageDialog(table,
				"I/O error occured during file export:\n" + e.getMessage(),
				"Error", JOptionPane.ERROR_MESSAGE);
		}
	}


	/**
	 * Exports descriptors to a given file.
	 *
	 * @param file is a File object with destination file
	 */
	private void doExport(File file) throws IOException
	{
		writer = new BufferedWriter(new FileWriter(file));
		currentRow = new ArrayList<String>();
		createHeaderRow(descMeta);

		StructureList sl = table.getStructureList();
		for (int i=0; i<sl.size(); ++i) {
			exportStructure(sl.get(i));
		}
		writer.close();
	}

	private void createHeaderRow(CDescriptorMetaMap meta) throws IOException {
		colMap = new HashMap<String, Integer>();

		currentRow.add("\"Structure Id\"");

		for (String id: meta) {
			String name = meta.getName(id);
			colMap.put(id, Integer.valueOf(currentRow.size()));
			currentRow.add('"' + name + '"');
		}

		writeRow();
	}

	/**
	 * Writes current row and cleans it.
	 */
	private void writeRow() throws IOException {
		for (int i = 0; i < currentRow.size(); ++i) {
			String cell = (String) currentRow.get(i);
			if (cell != null)
				writer.write(cell);
			writer.write("\t");
		}
		writer.newLine();
	}

	private void exportStructure(CStructure str) throws IOException {
		Collections.fill(currentRow, "");
		currentRow.set(0, str.getId());

		CDescriptorList dl = str.getDescriptorList();

		for (int i = 0; i < dl.getNumberOfDescriptors(); ++i) {
			CDescriptor d = dl.getDescriptor(i);
			Integer col = (Integer) colMap.get(d.getId());
			if (col != null) {
				currentRow.set(col.intValue(), Double.toString(d.getValue()));
			} else {
				System.out.println("Warning: descriptor metadata missing for descripotr " + d.getId());
			}
		}

		writeRow();
	}
	
}
