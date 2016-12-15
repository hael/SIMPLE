package org.chemomentum.gridbeans.common;

import java.awt.BorderLayout;
import java.io.InputStream;
import java.util.List;

import com.intel.gpe.clients.api.Client;
import com.intel.gpe.clients.api.exceptions.GPEFileException;
import com.intel.gpe.clients.api.transfers.GPEFile;
import com.intel.gpe.gridbeans.plugins.*;
import com.intel.gpe.gridbeans.plugins.swing.panels.GridBeanPanel;
import com.intel.gpe.gridbeans.*;

import org.openmolgrid.client.common.*;

import javax.xml.namespace.QName;

/**
 * Panel for visualising calculation log files.
 * 
 * @author Sulev Sild
 */
public class OutputSummaryPanel extends GridBeanPanel {

	private CalculationResultViewer crv;
	private QName errorLog;

	/**
	 * Constructs visualisation panel for calculation summary.
	 * 
	 * @param parent - a GPE client
	 * @param source - a key for locating summary file parameter in the GridBean
	 */
	public OutputSummaryPanel(Client parent, QName source) {
		super(parent, "Calculation summary");
		errorLog = source;
		setLayout(new BorderLayout());
		crv = new CalculationResultViewer();
		add(crv);
	}

	@Override
	public void load(Client client)
			throws DataSetException, TranslationException {
		super.load(client);
		IGridBeanModel model = getGridBeanModel();

		List<GPEFile> files = model.getOutputFiles(errorLog);
		if (files.size() > 0) {
			GPEFile gf = files.get(0);
			try {
				InputStream is = gf.getInputStream();
				crv.setInput(is);
				crv.updateValues();
			} catch (GPEFileException e) {
				client.showException("Error loading output for: " + errorLog, e);
			}
		}
	}
}
