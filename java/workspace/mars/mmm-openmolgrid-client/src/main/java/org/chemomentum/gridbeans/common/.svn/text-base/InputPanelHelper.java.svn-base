package org.chemomentum.gridbeans.common;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;

import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;

import org.openmolgrid.client.common.StructureImportPanel;

import com.intel.gpe.clients.api.Client;
import com.intel.gpe.gridbeans.IGridBean;
import com.intel.gpe.gridbeans.plugins.DataSetException;
import com.intel.gpe.gridbeans.plugins.swing.panels.GridBeanPanel;

/**
 * This class is provides a basic GUI for the preparing input
 * panels for GridBeans that process lists of molecular structures.
 *
 * @author Sulev Sild
 */
public abstract class InputPanelHelper extends GridBeanPanel {

	// name of the submitted job
	private JTextField jobTitle;

	// names of the input and output files
//	private JTextField inputFile;
//	private JTextField outputFile;

	private StructureImportPanel uploadPanel;

	private ApplicationFilter appFilter;

	protected JComboBox appCombo;
	protected JTabbedPane tabs;

	public InputPanelHelper(Client client, String title) {
		this(client, title, null);
	}

	public InputPanelHelper(Client client, String title, ApplicationFilter appFilter) {
		super(client, title);
		this.appFilter = appFilter;
		try {
			buildComponents();
		} catch (DataSetException e) {
			throw new RuntimeException("Couldn't create a GUI for "+getTitle());
		}
	}

	
	public abstract String getTitle();
	
//	@Override
//	public void load(Client client) throws DataSetException, TranslationException {
//		// update the list of available application names
//		if (appFilter != null) {
//			List<Application> l = appFilter.filter(getClient().getSelectedTargetSystem());
//			setPossibleValues(GridBeanHelper.APPLICATION, l);
//		}
//		super.load(client);
//
//		IGridBeanModel model = getGridBeanModel();
//		IFileParameterValue tmp;
//		tmp = (IFileParameterValue)model.get(GridBeanHelper.INPUT);
//		inputFile.setText(tmp.getTarget().getRelativePath());
//		tmp = (IFileParameterValue)model.get(GridBeanHelper.OUTPUT);
//		outputFile.setText(tmp.getSource().getRelativePath());
//		
//		
//	}
//	
//	@Override
//	public void store() throws DataSetException, TranslationException {
//		super.store();
//		IGridBeanModel model = getGridBeanModel();
//
//		String infile = inputFile.getText().trim();
//		IFileParameterValue ip;
//		ip = (IFileParameterValue)model.get(GridBeanHelper.INPUT);
//		ip.getTarget().setRelativePath(infile);
//
//		String outfile = outputFile.getText().trim();
//		IFileParameterValue op;
//		op = (IFileParameterValue)model.get(GridBeanHelper.OUTPUT);
//		op.getSource().setRelativePath(outfile);
//		
//		op = (IFileParameterValue)model.get(GridBeanHelper.ERROR);
//		op.getSource().setRelativePath(outfile+".errorlog");
//	}
	
	public void buildComponents() throws DataSetException
	{
		// create a panel on the top of the job preparation area 
		JPanel topPanel = new JPanel(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.insets = new Insets(3, 3, 3, 3);
		c.gridy = 1;
		jobTitle = SwingUtils.addTextEntryField("Job name:", c, topPanel);
		linkJobNameTextField(IGridBean.JOBNAME, jobTitle);
		
//		c.gridy += 1;
//		if (appFilter != null) {
//			appCombo = SwingUtils.addComboBox("Application:", new String[]{}, c, topPanel);
//			linkComboBox(GridBeanHelper.APPLICATION, appCombo);
//			setValueTranslator(GridBeanHelper.APPLICATION, ApplicationTranslator.getInstance());
//			c.gridy += 1;
//		}
//		inputFile = SwingUtils.addTextEntryField("Input file:", c, topPanel);
//		c.gridy += 1;
//		outputFile = SwingUtils.addTextEntryField("Output file:", c, topPanel);
//		topPanel.setBorder(
//			new TitledBorder("Job Preparation for " + getTitle()));

		tabs = new JTabbedPane();
		uploadPanel = new StructureImportPanel();
		StructureImportPanelControl ctrl = new StructureImportPanelControl(getClient(), GridBeanHelper.INPUT, uploadPanel);
		linkDataControl(GridBeanHelper.INPUT, ctrl);
		
		tabs.add("Structures", uploadPanel);

		setLayout(new BorderLayout());
		add(topPanel, BorderLayout.NORTH);
		add(tabs, BorderLayout.CENTER);
	}
	
}
