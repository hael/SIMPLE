package org.chemomentum.gridbeans.common;

import java.awt.BorderLayout;
import java.io.InputStream;
import java.util.List;

import javax.xml.namespace.QName;

import org.openmolgrid.client.common.StructureListViewer;

import com.intel.gpe.clients.api.Client;
import com.intel.gpe.clients.api.exceptions.GPEFileException;
import com.intel.gpe.clients.api.transfers.GPEFile;
import com.intel.gpe.gridbeans.IGridBeanModel;
import com.intel.gpe.gridbeans.plugins.DataSetException;
import com.intel.gpe.gridbeans.plugins.TranslationException;
import com.intel.gpe.gridbeans.plugins.swing.panels.GridBeanPanel;

/**
 * Panel for visualising output data in structure list files.
 * 
 * @author Sulev Sild
 */
public class OutputStructurePanel extends GridBeanPanel {
	private StructureListViewer structurePanel;
	private QName structureKey;
	
	/**
	 * Constructs visualisation panel for calculation results.
	 * 
	 * @param client - a GPE client
	 * @param source - the key of the output file parameter
	 */
	public OutputStructurePanel(Client client, QName source) {
		super(client, "Visualise structures");
		structureKey = source;
		setLayout(new BorderLayout());
		structurePanel = new StructureListViewer();
		add(structurePanel);
	}

	/* (non-Javadoc)
	 * @see com.intel.gpe.gridbeans.plugins.GridBeanPanel#load(com.intel.gpe.client2.Client)
	 */
	@Override
	public void load(Client client)
			throws DataSetException, TranslationException {
		super.load(client);
		IGridBeanModel model = getGridBeanModel();
		List<GPEFile> files = model.getOutputFiles(structureKey);
		if (files.size() > 0) {
			GPEFile gf = files.get(0);
			try {
				InputStream is = gf.getInputStream();
				structurePanel.setInput(is);
				structurePanel.updateValues();
			} catch (GPEFileException e) {
				client.showException("Error loading output for: "
						+ structureKey, e);
			}
		}
	}
}