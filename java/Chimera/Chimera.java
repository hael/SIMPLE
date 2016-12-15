/* vim: set ts=2: */
/**
 * Copyright (c) 2006 The Regents of the University of California.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions, and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above
 *      copyright notice, this list of conditions, and the following
 *      disclaimer in the documentation and/or other materials provided
 *      with the distribution.
 *   3. Redistributions must acknowledge that this software was
 *      originally developed by the UCSF Computer Graphics Laboratory
 *      under support by the NIH National Center for Research Resources,
 *      grant P41-RR01081.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 * OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
package structureViz.actions;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.io.*;
import java.awt.Color;

import cytoscape.Cytoscape;
import cytoscape.CyNode;
import cytoscape.data.CyAttributes;
import cytoscape.logger.CyLogger;
import cytoscape.view.*;

import structureViz.model.ChimeraModel;
import structureViz.model.ChimeraChain;
import structureViz.model.ChimeraResidue;
import structureViz.model.ChimeraStructuralObject;
import structureViz.model.Structure;
import structureViz.actions.CyChimera;
import structureViz.ui.ModelNavigatorDialog;
import structureViz.ui.AlignStructuresDialog;

/**
 * This class provides the main interface to UCSF Chimera
 * 
 * @author scooter
 *
 */
public class Chimera {
  /**
   * Static variables to keep track of the running
   * Chimera instance
   */
  // Chimera process
  static Process chimera;
	static private List<ChimeraModel> models;
	static private Map<Integer,ChimeraModel> modelHash;
	static private ListenerThreads listener;
	static private CyNetworkView networkView;
	static private ModelNavigatorDialog mnDialog = null;
	static private AlignStructuresDialog alDialog = null;
	static private List<ChimeraStructuralObject>selectionList = null;
	static Chimera staticPointer = null;
	static CyLogger logger;

	static int MAX_SUB_MODELS = 1000;

	static public Chimera GetChimeraInstance(CyNetworkView networkView, CyLogger logger) {
		if (staticPointer != null) {
			staticPointer.setNetworkView(networkView);
		} else {
			staticPointer = new Chimera(networkView, logger);
		}
		return staticPointer;
	}
    
	/**
	 * Create the Chimera object
	 *
	 * @param networkView the current CyNetworkView 
	 */
  protected Chimera(CyNetworkView networkView, CyLogger logger) {
		models = new ArrayList<ChimeraModel>();
		modelHash = new HashMap<Integer,ChimeraModel>();
		this.networkView = networkView;
		this.logger = logger;
		selectionList = new ArrayList<ChimeraStructuralObject>();
  }

	/**
	 * Return the list of all open models in this instance of Chimera
	 *
	 * @return list of ChimeraModels
	 */
	public List<ChimeraModel> getChimeraModels () { return models; }

	/**
 	 * Return a specific chimera model based on the model number
 	 *
 	 * @param model model number
 	 * @return the corresponding model
 	 */
	public ChimeraModel getChimeraModel (int model, int subModel) {
		Integer key = makeModelKey(model,subModel);
		if (modelHash.containsKey(key))
			return modelHash.get(key);
		return null;
	}

	/**
 	 * Return a specific chimera model based on the model number
 	 *
 	 * @param model model number
 	 * @return the corresponding model
 	 */
	public ChimeraModel getChimeraModel (int model) {
		return getChimeraModel(model, 0);
	}

	/**
 	 * Return the list of open structures in this instance of Chimera
 	 *
 	 * @return list of Structures
 	 */
	public List<Structure> getOpenStructs() {
		List<Structure>st = new ArrayList<Structure>();
		if (this.models == null)
			return st;

		for (ChimeraModel model: models) {
			Structure structure = model.getStructure();
			if (structure != null)
				st.add(structure);
		}
		return st;
	}

	/**
	 * Return our network view
	 *
	 * @return the network view we were created with
	 */
	public CyNetworkView getNetworkView () { return networkView; }

	/**
	 * Set our network view
	 *
	 * @param networkView the network view to switch to
	 */
	public void setNetworkView (CyNetworkView networkView) { this.networkView = networkView; }

	/**
	 * Provide a handle to the ModelNavigatorDialog that is currently
	 * up.  We need this so that we can update the selections in the
	 * dialog when the selections change in Chimera.
	 *
	 * @param dialog the ModelNavigatorDialog that provides our interface
	 */
	public void setDialog(ModelNavigatorDialog dialog) { mnDialog = dialog; }

	/**
	 * Return the ModelNavigatorDialog that is providing our interface
	 *
	 * @return our ModelNavigatorDialog
	 */
	public ModelNavigatorDialog getDialog() { return mnDialog; }

	public void launchDialog() {
		if (mnDialog == null) {
			mnDialog = ModelNavigatorDialog.LaunchModelNavigator(Cytoscape.getDesktop(), this);
		}
		mnDialog.setVisible(true);
	}

	/**
	 * Provide a handle to the AlignStructuresDialog that is currently
	 * up.  This does not necessarily need to be in the Chimera object,
	 * but it has to be somewhere, and this seemed like a good place.
	 *
	 * @param dialog the AlignStructuresDialog that provides our interface
	 */
	public void setAlignDialog(AlignStructuresDialog dialog) { this.alDialog = dialog; }

	/**
	 * Return the AlignStructuresDialog that provides the interface to our align structures
	 * functionality.
	 *
	 * @return our AlignStructuresDialog
	 */
	public AlignStructuresDialog getAlignDialog() { return this.alDialog; }

	/**
	 * Test to see if we currently have a particular model open.
	 *
	 * @param modelNumber the model number expressed as an integer
	 * @param subModelNumber the subModel number expressed as an integer
	 */
	public boolean containsModel(int modelNumber, int subModelNumber) {
		return modelHash.containsKey(makeModelKey(modelNumber, subModelNumber));
	}

	/**
 	 * Return true if this structure is currently open
 	 *
 	 * @param structure the Structure we're inquiring about
 	 * @return true if open
 	 */
	public boolean containsModel(Structure structure) {
		// Get the model number
		int modelNumber = structure.modelNumber();
		int subModelNumber = structure.subModelNumber();
		return containsModel(modelNumber, subModelNumber);
	}

	/**
	 * Test to see if we currently have a particular model open.
	 *
	 * @param modelNumber the model number expressed as an int
	 */
	public boolean containsModel(int modelNumber) {
		return containsModel(modelNumber, 0);
	}

	/**
	 * Return the ChimeraModel associated with the requested 
	 * model number
	 *
	 * @param modelNumber the model number expressed as an int
	 * @param subModelNumber the subModel number expressed as an int
	 * @return the ChimeraModel with a model number of modelNumber
	 */
	public ChimeraModel getModel(int modelNumber, int subModelNumber) {
		Integer key = makeModelKey(modelNumber, subModelNumber);
		if (modelHash.containsKey(key))
			return modelHash.get(key);
		return null;
	}

	/**
	 * Return the ChimeraModel associated with the requested 
	 * model number
	 *
	 * @param modelNumber the model number expressed as an int
	 * @return the ChimeraModel with a model number of modelNumber
	 */
	public ChimeraModel getModel(int modelNumber) {
		return getModel(modelNumber, 0);
	}

	/**
	 * Return the ChimeraModel associated with the requested 
	 * model name
	 *
	 * @param modelName the model name expressed as a string
	 * @return the ChimeraModel with a model name of modelName
	 */
	public ChimeraModel getModel(String modelName) {
		for (ChimeraModel model: models) {
			if (model.getModelName().equals(modelName))
				return model;
		}
		return null;
	}

	/**
 	 * Return the list of currently selected structural objects
 	 *
 	 * @return selection list
 	 */
	public List<ChimeraStructuralObject> getSelectionList() {
		return selectionList;
	}

	/**
 	 * Add a selection to the selection list.  This is called primarily by
 	 * the Model Navigator Dialog to keep the selections in sync
 	 *
 	 * @param selectionToAdd the selection to add to our list
 	 */
	public void addSelection(ChimeraStructuralObject selectionToAdd) {
		if (selectionToAdd != null && !selectionList.contains(selectionToAdd))
			selectionList.add(selectionToAdd);
	}

	/**
 	 * Remove a selection from the selection list.  This is called primarily by
 	 * the Model Navigator Dialog to keep the selections in sync
 	 *
 	 * @param selectionToRemove the selection to remove from our list
 	 */
	public void removeSelection(ChimeraStructuralObject selectionToRemove) {
		if (selectionToRemove != null && selectionList.contains(selectionToRemove))
			selectionList.remove(selectionToRemove);
	}

	/**
 	 * Clear the list of selected objects
 	 */
	public void clearSelectionList() {
		for (ChimeraStructuralObject cso: selectionList) {
			if (cso != null)
				cso.setSelected(false);
		}
		selectionList.clear();
	}

	/**
	 * Test to see if we currently have a running instance of Chimera
	 *
	 * @return <b>true</b> if Chimera is already running
	 */
	public boolean isLaunched () {
		if (chimera != null) 
			return true;
		return false;
	}
    
  /**
   * Launch (start) an instance of Chimera
   * @return "true" if the launch was successful
   * @throws IOException
 */
  public boolean launch() throws IOException {
 		// See if we already have a chimera instance running
		if (isLaunched() == true)
			return true;
  		
 		// No, get one started
		// Figure out our fallback
		String osPath = "";
		String os = System.getProperty("os.name");
		if (os.startsWith("Linux")) {
			osPath = "/usr/local/chimera/bin/";
		} else if (os.startsWith("Windows")) {
			osPath = "\\Program Files\\Chimera\\bin\\";
		} else if (os.startsWith("Mac")) {
			osPath = "/Applications/Chimera.app/Contents/MacOS/";
		}

		// Get the path
		String path = CyChimera.getProperty("chimeraPath");
		if (path != null) {
			path = path+"chimera";
		} else {
			path = "chimera";
		}
		CyLogger.getLogger(Chimera.class).info("Path: "+path);

		try {
 			List <String> args = new ArrayList<String>();
 			args.add(path);
 			args.add("--start");
 			args.add("ReadStdin");
 			ProcessBuilder pb = new ProcessBuilder(args);
 			chimera = pb.start();
		} catch (Exception e) {
 			List <String> args = new ArrayList<String>();
 			args.add(osPath+"chimera");
 			args.add("--start");
 			args.add("ReadStdin");
 			ProcessBuilder pb = new ProcessBuilder(args);
 			chimera = pb.start();
		}

		// Start up a listener
		listener = new ListenerThreads(chimera, this, logger);
		listener.start();

		// Ask Chimera to give us updates
		chimeraSend("listen start models; listen start selection");
		return true;
  }
  
  /**
   * Open a Chimera model
	 *
   * @param structure the Structure to open
   */
  public void open(Structure structure) {
		structure.setModelNumber(Structure.getNextModel(), 0);
		if (structure.getType() == Structure.StructureType.MODBASE_MODEL)
			chimeraSend("listen stop models; listen stop selection; open "+structure.modelNumber()+" modbase:"+structure.name());
		else if (structure.getType() == Structure.StructureType.SMILES)
			chimeraSend("listen stop models; listen stop selection; open "+structure.modelNumber()+" smiles:"+structure.name());
		else
			chimeraSend("listen stop models; listen stop selection; open "+structure.modelNumber()+" "+structure.name());

		// Now, figure out exactly what model # we got
		List<ChimeraModel> modelList = getModelInfoList(structure);
		if (modelList == null) return;
		for (ChimeraModel newModel: modelList) {

			// Get the color (for our navigator)
			newModel.setModelColor(getModelColor(newModel));

			// Get our properties (default color scheme, etc.)
			// Make the molecule look decent
			// chimeraSend("repr stick "+newModel.toSpec());

			if (structure.getType() != Structure.StructureType.SMILES) {
				// Create the information we need for the navigator
				getResidueInfo(newModel);
			}

			// Add it to our list of models
			models.add(newModel);

			// Add it to the hash table
			modelHash.put(makeModelKey(newModel.getModelNumber(),newModel.getSubModelNumber()),newModel);
		}

		chimeraSend("focus");
		chimeraSend("listen start models; listen start selection");

		if (mnDialog != null)
			mnDialog.modelChanged();

  	return;
  }

	/**
	 * Close a Chimera model
	 *
   * @param structure the Structure to close
	 */
	public void close(Structure structure) {
		int model = structure.modelNumber();
		int subModel = structure.subModelNumber();
		Integer modelKey = makeModelKey(model, subModel);
		if (modelHash.containsKey(modelKey)) {
			ChimeraModel chimeraModel = modelHash.get(modelKey);
			chimeraSend("listen stop models; listen stop select; close "+chimeraModel.toSpec());
			models.remove(chimeraModel);
			modelHash.remove(modelKey);
			selectionList.remove(chimeraModel);
		} else {
			chimeraSend("listen stop models; listen stop select; close #"+model+"."+subModel);
		}
		chimeraSend("listen start models; listen start select");
		Structure.close(structure);
		return;
	}

	/**
	 * Select something in Chimera
	 *
	 * @param command the selection command to pass to Chimera
	 */
	public void select(String command) {
		chimeraSend("listen stop select; "+command+"; listen start select");
	}

	public void selectNoReply(String command) {
		chimeraSendNoReply("listen stop select; "+command+"; listen start select");
	}

	public List<String> commandReply(String text) {
		List<String> r = chimeraSend(text);
		if (r == null) 
			r = new ArrayList<String>();
		return r;
	}

  /**
   * Send a string to the Chimera instance
	 *
   * @param text the text to pass to Chimera
   */
  public List<String> chimeraSend(String command) {
  	if (chimera == null)
  		return null;

		// System.out.println("To Chimera --> "+command);
		listener.clearResponse(command);
		String text = command.concat("\n");

		try {
  		// send the command
  		chimera.getOutputStream().write(text.getBytes());
  		chimera.getOutputStream().flush();
		} catch (IOException e) {
			CyLogger.getLogger(Chimera.class).warning("Unable to execute command: "+text);
			CyLogger.getLogger(Chimera.class).warning("Exiting...");
			chimera = null;
			if (mnDialog != null)
				mnDialog.lostChimera();
		}

		return listener.getResponse(command);
  }

  public void chimeraSendNoReply(String command) {
  	if (chimera == null)
  		return;

		// System.out.println("To Chimera --> "+command);
		listener.clearResponse(command);
		String text = command.concat("\n");
		try {
  		// send the command
  		chimera.getOutputStream().write(text.getBytes());
  		chimera.getOutputStream().flush();
		} catch (IOException e) {
			CyLogger.getLogger(Chimera.class).warning("Unable to execute command: "+text);
			CyLogger.getLogger(Chimera.class).warning("Exiting...");
			chimera = null;
			if (mnDialog != null)
				mnDialog.lostChimera();
		}
		return;
	}
  
  /**
   * Terminate the running Chimera process
   * 
   */
  public void exit() {
  	if (chimera != null) {
  		chimeraSend("stop really");
  		chimera.destroy();
  		chimera = null;
			models = null;
			modelHash = null;
			staticPointer = null;
		}
		if (alDialog != null)
			alDialog.setVisible(false);

		if (mnDialog != null)
			mnDialog.setVisible(false);

		Structure.closeAll();
  }

	/**
	 * Dump and refresh all of our model/chain/residue info
	 */
	public void refresh() {
		// Get a new model list
		HashMap<Integer,ChimeraModel> newHash = new HashMap<Integer,ChimeraModel>();

		// Stop all of our listeners while we try to handle this
		chimeraSend("listen stop select; listen stop models");

		// Get all of the open models
		List<ChimeraModel> newModelList = getModelList();

		// Match them up -- assume that the model #'s haven't changed
		for (ChimeraModel model: newModelList) {
			// Get the color (for our navigator)
			model.setModelColor(getModelColor(model));

			// Get our model info
			int modelNumber = model.getModelNumber();
			int subModelNumber = model.getSubModelNumber();

			// If we already know about this model number, get the Structure,
			// which tells us about the associated CyNode
			if (containsModel(modelNumber, subModelNumber)) {
				ChimeraModel oldModel = getModel(modelNumber, subModelNumber);
				Structure s = oldModel.getStructure();
				if (s.getType() == Structure.StructureType.SMILES)
					model.setModelName(s.name());
				model.setStructure(s);
			} else {
				// This will return a new Structure if we don't know about it
				Structure s = CyChimera.findStructureForModel(networkView, model.getModelName(), true);
				s.setModelNumber(model.getModelNumber(), model.getSubModelNumber());
				model.setStructure(s);
			}

			newHash.put(makeModelKey(model.getModelNumber(), model.getSubModelNumber()),model);

			if (model.getStructure() != null && model.getStructure().getType() != Structure.StructureType.SMILES) {
				// Get the residue information
				getResidueInfo(model);
			}
		}

		// Replace the old model list
		models = newModelList;
		modelHash = newHash;

		// Restart all of our listeners
		chimeraSend("listen start models; listen start select");

		// Done
	}

	/**
	 * Inform our interface that the model has changed
	 */
	public void modelChanged() {
		if (mnDialog != null)
			mnDialog.modelChanged();
	}

	/**
	 * Inform our interface that the selection has changed
	 */
	public void updateSelection(List<ChimeraStructuralObject> selectionList) {
		if (mnDialog != null)
			mnDialog.updateSelection(selectionList);
	}

	/**
	 * This is called by the selectionListener to let us know that
	 * the user has changed their selection in Chimera.  We need
	 * to go back to Chimera to find out what is currently selected
	 * and update our list.
	 */
	public void updateSelection() {
		HashMap<Integer,ChimeraModel>modelSelHash = new HashMap<Integer, ChimeraModel>();
		clearSelectionList();

		// Execute the command to get the list of models with selections
		for (String modelLine: commandReply("lists level molecule")) {
			ChimeraModel chimeraModel = new ChimeraModel(modelLine);
			Integer modelKey = makeModelKey(chimeraModel.getModelNumber(), chimeraModel.getSubModelNumber());
			modelSelHash.put(modelKey, chimeraModel);
		}

		// Now get the residue-level data
		for (String inputLine: commandReply("lists level residue")) {
			ChimeraResidue r = new ChimeraResidue(inputLine);
			Integer modelKey = makeModelKey(r.getModelNumber(), r.getSubModelNumber());
			if (modelSelHash.containsKey(modelKey)) {
				ChimeraModel model = modelSelHash.get(modelKey);
				model.addResidue(r);
			}
		}

		// Get the selected objects
		for (ChimeraModel selectedModel: modelSelHash.values()) {
			int modelNumber = selectedModel.getModelNumber();
			int subModelNumber = selectedModel.getSubModelNumber();
			// Get the corresponding "real" model
			if (containsModel(modelNumber, subModelNumber)) {
				ChimeraModel dataModel = getModel(modelNumber, subModelNumber);
				if (dataModel.getResidueCount() == selectedModel.getResidueCount() ||
				    dataModel.getStructure().getType() == Structure.StructureType.SMILES) {
					// Select the entire model
					selectionList.add(dataModel);
					dataModel.setSelected(true);
				} else {
					for (ChimeraChain selectedChain: selectedModel.getChains()) {
						ChimeraChain dataChain = dataModel.getChain(selectedChain.getChainId());
						if (selectedChain.getResidueCount() == dataChain.getResidueCount()) {
							selectionList.add(dataChain);
							dataChain.setSelected(true);
						} else {
							// Need to select individual residues
							for (ChimeraResidue res: selectedChain.getResidues()) {
								String residueIndex = res.getIndex();
								ChimeraResidue residue = dataChain.getResidue(residueIndex);
								if (residue == null) continue;
								selectionList.add(residue);
								residue.setSelected(true);
							} // resIter.hasNext
						}
					} // chainIter.hasNext()
				}
			}
		} // modelIter.hasNext()

		// Finally, update the navigator panel
		updateSelection(selectionList);
	}

	/**
	 * Return the list of depiction presets available from within Chimera.
	 * Chimera will return the list as a series of lines with the format:
	 *	Preset type number "description"
	 *
	 * @return list of presets
	 */
	public List getPresets() {
		ArrayList<String>presetList = new ArrayList<String>();
		for (String preset: commandReply ("preset list")) {
			preset = preset.substring(7); // Skip over the "Preset"
			preset = preset.replaceFirst("\"", "(");
			preset = preset.replaceFirst("\"", ")");
			// string now looks like: type number (description)
			presetList.add(preset);
		}
		return presetList;
	}

	/**
	 * Return the list of ChimeraModels currently open
	 *
	 * @return List of ChimeraModel's
	 */
	private List getModelList() {
		List<ChimeraModel>modelList = new ArrayList<ChimeraModel>();
		List<String> list = commandReply("listm type molecule");
		if (list != null) {
			for (String modelLine: list) {
				ChimeraModel chimeraModel = new ChimeraModel(modelLine);
				modelList.add(chimeraModel);
			}
		}
		return modelList;
	}

	/**
	 * Return the ChimeraModel associated with this structure.  This
	 * involves asking Chimera to provide us with the information.
	 *
	 * @param structure the Structure we want to get the ChimeraModel for
	 * @return the ChimeraModel for this structure
	 */
	private List<ChimeraModel> getModelInfoList(Structure structure) {
		String name = structure.name();
		int modelNumber = structure.modelNumber();
		int subModelNumber = structure.subModelNumber();
		List<ChimeraModel>infoList = new ArrayList<ChimeraModel>();

		List<String>replyList = null;

		if (subModelNumber == 0)
			replyList = commandReply("listm type molecule spec #"+modelNumber);
		else
			replyList = commandReply("listm type molecule spec #"+modelNumber+"."+subModelNumber);

		for (String modelLine: replyList) {
			if (subModelNumber == 0 && modelLine.contains("id #"+modelNumber)) {
				// got the right model, now get the model number
				ChimeraModel chimeraModel = new ChimeraModel(structure, modelLine);
				// System.out.println("Identified model as "+chimeraModel);
				infoList.add(chimeraModel);
			} else if (modelLine.contains("id #"+modelNumber+"."+subModelNumber)) {
				ChimeraModel chimeraModel = new ChimeraModel(structure, modelLine);
				infoList.add(chimeraModel);
			}
		}
		if (infoList.size() > 0)
			return infoList;
		return null;
	}

	/**
	 * Determine the color that Chimera is using for this model.
	 *
	 * @param model the ChimeraModel we want to get the Color for
	 * @return the default model Color for this model in Chimera
	 */
	private Color getModelColor(ChimeraModel model) {
		List<String> colorLine = chimeraSend ("listm type molecule attr color spec "+model.toSpec());
		if (colorLine == null) return null;
		String inputLine = (String)colorLine.get(0);
		int colorStart = inputLine.indexOf("color ");
		String colorString = inputLine.substring(colorStart+6);
		String[] rgbStrings = colorString.split(",");
		float[] rgbValues = new float[4];
		for (int i = 0; i < rgbStrings.length; i++) {
			Float f = new Float(rgbStrings[i]);
			rgbValues[i] = f.floatValue();
		}
		if (rgbStrings.length == 4) {
			return new Color(rgbValues[0], rgbValues[1], rgbValues[2], rgbValues[3]);
		} else {
			return new Color(rgbValues[0], rgbValues[1], rgbValues[2]);
		}
	}

	/**
	 * Get information about the residues associated with a model.  This
	 * uses the Chimera listr command.  We don't return the resulting
	 * residues, but we add the residues to the model.
	 *
	 * @param model the ChimeraModel to get residue information for
	 * 
	 */
	private void getResidueInfo(ChimeraModel model) {
		int modelNumber = model.getModelNumber();
		int subModelNumber = model.getSubModelNumber();

		// Get the list -- it will be in the reply log
		for (String inputLine: commandReply ("listr spec #"+modelNumber+"."+subModelNumber)) {
			ChimeraResidue r = new ChimeraResidue(inputLine);
			if (r.getModelNumber() == modelNumber || r.getSubModelNumber() == subModelNumber) {
				model.addResidue(r);
			}
		}
	}

	/**
 	 * Create the key to use for forming the model/submodel key into the modelHash
 	 *
 	 * @param model the model number
 	 * @param subModel the submodel number
 	 * @return the model key as an Integer
 	 */
	private Integer makeModelKey(int model, int subModel) {
		return new Integer(model*MAX_SUB_MODELS+subModel);
	}
}
