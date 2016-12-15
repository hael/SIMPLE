/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common;

import java.io.File;

import org.openmolgrid.client.common.table.StructureList;
import org.openmolgrid.format.cml.CMLReader;
import org.openmolgrid.format.csv.CSVReader;
import org.openmolgrid.format.gamess.GamessOutput;
import org.openmolgrid.format.sdf.SDFParser;
import org.openmolgrid.format.slf.StructureReader;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.ChemLibException;
import org.openmolgrid.qsar.CDescriptorList;
import org.openmolgrid.qsar.CDescriptorMetaMap;
import org.openmolgrid.qsar.CPropertyList;
import org.openmolgrid.qsar.CPropertyMetaMap;
import org.openmolgrid.qsar.DataDictionary;
import org.openmolgrid.qsar.StructureConsumer;
import org.openmolgrid.util.CUtil;

//TODO: add file format detection based on the file's content.

/**
 * This class implements importer for chemical file formats that are supported
 * by OpenMolGRID plugins. This class automatically detects imported file
 * formats.  Currently the detection is based on the file name extension.
 * Supported file formats include MDL's MOL and SD files, CML, and custom XML
 * format.
 *
 * @author Sulev Sild
 * @author Andre Lomaka
 */
public class StructureImporter implements StructureConsumer {
	private File file;
	private CPropertyMetaMap propMap = new CPropertyMetaMap();
	private CDescriptorMetaMap descMap = new CDescriptorMetaMap();
	private StructureList sl;
	
	private DataFieldMapper sdfMapper = new DataFieldMapper(null);
	private boolean mergeMetadata = false;
	private StructureReader slfReader;
	
	public StructureImporter(File file, StructureList sl) {
		this.file = file;
		this.sl = sl;
	}

	public CPropertyMetaMap getPropertyMetaMap() {
		return propMap;
	}
	
	public CDescriptorMetaMap getDescriptorMetaMap() {
		return descMap;
	}
	
	/**
	 * Starts parser.
	 * 
	 * @exception ChemLibException on unknown file format or parse error
	 */
	public void provideStructures() throws ChemLibException {
		String ext = CUtil.getExtension(file);
		String fileName = file.getPath();
		if ("slf".equalsIgnoreCase(ext)) {
			try {
				slfReader = new StructureReader(fileName, this);
				mergeMetadata = true;
				slfReader.provideStructures();
			} finally {
				slfReader = null;
			}
		} else if ("sdf".equalsIgnoreCase(ext) || "mol".equalsIgnoreCase(ext) || "mdl".equalsIgnoreCase(ext)) {
			SDFParser sp = new SDFParser(fileName);
			try {
				CStructure str=null;
				while ((str = sp.getNextStructure()) != null) {
					consumeStructure(str);
				}
			} finally {
				sp.close();
			}
		} else if ("cml".equalsIgnoreCase(ext))
		{
			CMLReader cr = new CMLReader(fileName, null, this);
			cr.provideStructures();
		} else if ("tsv".equalsIgnoreCase(ext))
		{
			CSVReader csv = new CSVReader(fileName, this);
			csv.provideStructures();
			propMap = csv.getPropertyMetaData(); // ????
		} else if ("log".equalsIgnoreCase(ext)) {
			CStructure s = new CStructure();
			s.setId(file.getName().split("\\.")[0]);
			GamessOutput data = new GamessOutput(fileName, null);
			s.addAttribute(data);
			consumeStructure(s);
		} else
		{
			throw new ChemLibException("Unrecognized file format!");
		}
	}

	private void handleMdlData(CStructure str) {
		DataDictionary data = str.getAttribute(DataDictionary.class);
		if (data == null) return;
		sdfMapper.apply(data, str);
	}

	private void mergeDescriptors(CStructure dst, CStructure src) {
		CDescriptorList dl_src = src.getAttribute(CDescriptorList.class);
		if (dl_src == null) {
			return;
		}
		CDescriptorList dl_dst = dst.getOrCreateDescriptorList();
		for (int i=0; i<dl_src.getNumberOfDescriptors(); ++i) {
			dl_dst.addDescriptor(dl_src.getDescriptor(i));
		}
	}
	
	private void mergeProperties(CStructure dst, CStructure src) {
		CPropertyList pl_src = src.getAttribute(CPropertyList.class);
		if (pl_src == null) {
			return;
		}
		CPropertyList pl_dst = dst.getOrCreatePropertyList();
		for (int i=0; i<pl_src.getNumberOfProperties(); ++i) {
			pl_dst.addProperty(pl_src.getProperty(i));
		}
	}


	public void startHandling() {
	}
	
	public void endHandling() {
	}
	
	/**
	 * Processes a new record from StructureProvider.
	 */
	public void consumeStructure(CStructure str) {
		// merge meta data on the first structure that is loaded from SLF file
		if (mergeMetadata) {
			mergeMetadata = false;
			sl.getDescriptorMetaData().mergeMetaData(slfReader.getDescriptorMetaData());
			sl.getPropertyMetaData().mergeMetaData(slfReader.getPropertyMetaData());
		}

		handleMdlData(str);
		
		int merge_to = sl.indexById(str.getId());
		if (merge_to != -1) {
			CStructure tmp = sl.get(merge_to);
			mergeDescriptors(tmp, str);
			mergeProperties(tmp, str);
			sl.set(merge_to, tmp);
			return;
		}
		
		if ("".equals(str.getId())) {
			str.setId(Integer.toString(sl.size()+1));
		}
		sl.add(str);
	}

}
