package org.chemomentum.gridbeans.common;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.IOException;
import java.util.List;

import javax.xml.namespace.QName;

import com.intel.gpe.clients.api.Client;
import com.intel.gpe.clients.api.async.IProgressListener;
import com.intel.gpe.clients.api.transfers.ProtocolConstants;
import com.intel.gpe.clients.common.transfers.GridFileAddress;
import com.intel.gpe.gridbeans.parameters.IFileParameterValue;
import com.intel.gpe.gridbeans.plugins.DataSetException;
import com.intel.gpe.gridbeans.plugins.IDataControlWithExternalDataSource;
import com.intel.gpe.gridbeans.plugins.TranslationException;
import com.intel.gpe.gridbeans.plugins.UnsupportedProtocolException;
import com.intel.gpe.gridbeans.plugins.swing.controls.SwingDataControl;

import org.openmolgrid.client.common.StructureImportPanel;

/**
 * Data control for StructureImportPanel.
 * 
 * @author Sulev Sild
 */
public class StructureImportPanelControl extends SwingDataControl implements IDataControlWithExternalDataSource, PropertyChangeListener {
	
	private IFileParameterValue fileParam;
	private StructureImportPanel importPanel;
	
	private File localFile;
	private boolean modified = false;

	public StructureImportPanelControl(Client client, QName param, StructureImportPanel panel) {
		super(client, param, panel);
		importPanel = panel;
		importPanel.addPropertyChangeListener(this);
	}

	public void saveDataToExternalSource(IProgressListener pl) throws DataSetException {
		importPanel.removePropertyChangeListener(this);
		try {
			if (!isLocalSource(fileParam)) return;
			if (modified) {
				importPanel.saveStructures(localFile);
				modified = false;
			}
		} catch (IOException e) {
			throw new DataSetException("Unable to save input structures: "+localFile);
		} finally {
			importPanel.addPropertyChangeListener(this);
		}
	}

	public Object getValue() throws DataSetException {
		if (isLocalSource(fileParam)) {
			File file = getSourceFile(fileParam);
			if (!file.equals(localFile)) {
				String tempDir = getClient().getFileFactory().getTemporaryDirName();
				
				if (localFile == null) {
					localFile = file;
					modified = true;
					
					// make sure that the content of the old file is not overwritten  
					if (localFile.exists()) {
						try {
							localFile = File.createTempFile("input", ".slf", new File(tempDir));
							localFile.delete();
						} catch (IOException e) {
							throw new DataSetException("Can't create temp file in: "+tempDir);
						}
					}
				}
				
				String fname = localFile.getAbsolutePath();
				fname = fname.startsWith(tempDir) ? fname.replace(tempDir, "") : fname; 
				fileParam.setSource(new GridFileAddress(fname));
			}
		}
		return fileParam;
	}

	@SuppressWarnings("unchecked")
	public void setPossibleValues(List arg0) throws TranslationException, DataSetException {
	}

	public void setValue(Object value) throws TranslationException, UnsupportedProtocolException, DataSetException {
		importPanel.removePropertyChangeListener(this);
		try {
			fileParam = (IFileParameterValue) value;
			
			// clear input panel since, we do not edit remote files
			if (!isLocalSource(fileParam)) {
				importPanel.clearStructures();
				localFile = null;
				modified = false;
				return;
			}

			File file = getSourceFile(fileParam);
			if (!file.canRead()) return;
			if (file.equals(localFile) && file.lastModified() == localFile.lastModified()) 
				return;

			try {
				importPanel.loadStructures(file);
				localFile = file;
				modified = false;
			} catch (IOException e) {
				throw new DataSetException("Could not load: "+file, e);
			}
		}
		finally {
			importPanel.addPropertyChangeListener(this);
		}
	}

	public void propertyChange(PropertyChangeEvent e) {
		if (e.getPropertyName().equals(StructureImportPanel.FILENAME)) {
			// TODO: what happens if file parameter is not local file?
			localFile = (File)e.getNewValue();
			fireValueChange();
		} else if (e.getPropertyName().equals(StructureImportPanel.MODIFIED)) {
			modified = true;
			fireValueChange();
		}
	}
	
	private boolean isLocalSource(IFileParameterValue fp) {
		return ProtocolConstants.LOCAL_PROTOCOL.equals(fp.getSource().getProtocol());
	}
	
	private File getSourceFile(IFileParameterValue fp) {
		File file = new File(fp.getSource().getInternalString());
		if (!file.isAbsolute()) {
			String tmpDir = getClient().getFileFactory().getTemporaryDirName();
			file = new File(tmpDir, file.getPath());
		}
		return file;
	}

}