/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.format.csv;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.ChemLibException;
import org.openmolgrid.qsar.CPropertyMetaMap;
import org.openmolgrid.qsar.DataDictionary;
import org.openmolgrid.qsar.StructureConsumer;
import org.openmolgrid.qsar.StructureProvider;

/**
 * Imports structures from CSV files. Currently tab characters are used as delimiters between columns.
 * 
 * @author Sulev Sild
 */
public class CSVReader implements StructureProvider {

	private File input;
	private StructureConsumer sc;

	public CSVReader(String fileName, StructureConsumer strConsumer) {
		input = new File(fileName);
		sc = strConsumer;
	}

	public CPropertyMetaMap getPropertyMetaData() {
		CPropertyMetaMap map = new CPropertyMetaMap();
		return map;
	}

	public void provideStructures() throws ChemLibException {
		BufferedReader r;
		try {
			r = new BufferedReader(new FileReader(input));
		} catch (FileNotFoundException e) {
			throw new ChemLibException("Can't find file " + input.getAbsolutePath());
		}
		sc.startHandling();
		try {
			int lineno = 1;
			String line = r.readLine();
			if (line == null) {
				throw new IOException("Unexpected end of file");
			}
			line = line.replace("\"", "");
			String[] headers = line.split("\t");
			while ((line = r.readLine()) != null) {
				lineno += 1;
				String[] vals = line.split("\t");
				if (vals.length == 1) {
					break;
				}
				if (headers.length != vals.length) {
					throw new ChemLibException("Incorrect number of columns on the line: " + lineno);
				}

				DataDictionary data = new DataDictionary();
				for (int i = 0; i < headers.length; ++i) {
					data.put(headers[i], vals[i]);
				}

				CStructure s = new CStructure();
				if (data.size() > 0) {
					s.addAttribute(data);
				}
				sc.consumeStructure(s);
			}
		} catch (IOException e) {
			throw new ChemLibException("Parse error in " + input.getAbsolutePath());
		} finally {
			sc.endHandling();
			try {
				r.close();
			} catch (IOException e) {
				throw new ChemLibException("Error while closing");
			}
		}
	}

}
