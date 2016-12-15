package org.openmolgrid.format.yaml;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.util.Map;
import java.util.Set;
import java.util.zip.DataFormatException;

import javax.swing.tree.DefaultMutableTreeNode;

import org.omegahat.Environment.Databases.ObjectNotFoundException;
import org.openmolgrid.model.CSystem;
import org.yaml.snakeyaml.Yaml;

public class YamlParser {

	/** Logger */
	//private static Logger log = Logger.getLogger(YamlParser.class.getName());

	// The Yaml variables
	private Yaml yaml = new Yaml();
	//files variables
	private boolean isFullPath;
	private String yamlfile;
	private String workingDir;
	private String atHome;
	private String atWork;
	private String path;
	private InputStream inputstream;
	//the tree variables
	private DefaultMutableTreeNode root;
	private DefaultMutableTreeNode subTreeNode = null;
	private Iterable<Object> data;
	private Map<?, ?> mapNextLevel0;
	
	/**
	 * The constructor for the class YamlParser
	 * @param filename
	 * @param rootsys
	 */
	public YamlParser(String filename, DefaultMutableTreeNode rootsys) {
		yamlfile = filename;
		root = rootsys;
	}
	public YamlParser(Boolean isFullPath, String filename, DefaultMutableTreeNode rootsys) {
		this.isFullPath = isFullPath;
		yamlfile = filename;
		root = rootsys;
	}
	/**
	 * the methods to construct the tree recursively from the input maps
	 * @return
	 */
	public DefaultMutableTreeNode doContructionTree() {
		MapTreeSystem mapTreeSystem = new MapTreeSystem(root, mapNextLevel0);
		//log.info("mapTreeSystem.getRoot(): "+mapTreeSystem.getRoot());
		//System.out.println("get child" +mapTreeSystem.getChild(root, 0, 0));
		Set<?> set0 = mapNextLevel0.keySet();
        Object[] arrayNextLevel0 = new String[set0.size()];
        arrayNextLevel0 = set0.toArray();
        for ( int elemNextLevel0 = 0 ; elemNextLevel0 < set0.size() ; elemNextLevel0++ ) {
        	//log.info(elemNextLevel0+" arrayNextLevel0: " +arrayNextLevel0[elemNextLevel0]);
        	//log.info(elemNextLevel0+"   getFirstChild: " +mapTreeSystem.getFirstChild(arrayNextLevel0, 0, elemNextLevel0));
        	mapTreeSystem.getFirstChild(arrayNextLevel0, 0, elemNextLevel0);
        }
		return subTreeNode;
	}
	/**
	 * methods to setup environment variables. 
	 * @return
	 */
	public String setWorkingFilePath() {
		if (isFullPath == false) {
			path = workingDir+""+yamlfile;
			//log.info("in class: YamlParser in method setWorkingFilePath path: "+ path);
		} else if (isFullPath == true) {
			path = yamlfile;
			//log.info("in class: YamlParser in method setWorkingFilePath path: "+ path);
		}
		setPath(path);
		return path;
	}
	public InputStream openYamlStreamInput() throws FileNotFoundException {
		inputstream = new FileInputStream(new File(path));
		setInputstream(inputstream);
		return inputstream;
	}
	public Iterable<Object> /**Object*/ objectData() throws ObjectNotFoundException {
		data = yaml.loadAll(getInputstream());
		setData(data);
		return data;
	}
	public Map<?, ?> mapLevelData() throws DataFormatException {
		Object lastYaml = null;
		for (Object currentYaml : data) {
			//System.out.println(currentYaml);
			//log.info(currentYaml);
			lastYaml = currentYaml;
		}
		mapNextLevel0 = (Map<?, ?>) lastYaml;

		setMapNextLevel0(mapNextLevel0);
		return mapNextLevel0;
	}
	/**---------------------------------------------Geometry routines------------------------//
	 * methods to create geometry of the system
	 * @return
	 */
	public CSystem getFinalGeometry() {
		CSystem csystem = new CSystem();
		return csystem;
	}
	/**--------------------------------------------------------------------------------------*/
	/** 
	 * Method to get the toString and check if variables are stored
	 * properly in the class. 
	 */
	public String toString() {
		return "The summary of the variables for the class RecursiveTreeConstruct: \n"+
				"atHome: "+getAtHome()+"\n"+
				"atWork: "+getAtWork()+"\n"+
				"Yaml file considered: "+getYamlfile()+"\n"+
				"path: "+getPath()+"\n"+
				"inputstream: "+getInputstream()+"\n"+
				"data: "+getData()+"\n"+
				"root: "+getRoot()+"\n"+
				"mapNexLvel0: "+getMapNextLevel0()+"\n";
	}
	/**-------------------------------------------*/
	/** The setters and getters for the variables */

	/**
	 * @return the yamlfile
	 */
	public String getYamlfile() {
		return yamlfile;
	}
	/**
	 * @param yamlfile the yamlfile to set
	 */
	public void setYamlfile(String yamlfile) {
		this.yamlfile = yamlfile;
	}
	/**
	 * @return the atHome
	 */
	public String getAtHome() {
		return atHome;
	}
	/**
	 * @param atHome the atHome to set
	 */
	public void setAtHome(String atHome) {
		this.atHome = atHome;
	}
	/**
	 * @return the atWork
	 */
	public String getAtWork() {
		return atWork;
	}
	/**
	 * @param atWork the atWork to set
	 */
	public void setAtWork(String atWork) {
		this.atWork = atWork;
	}
	/**
	 * @return the workingDir
	 */
	public String getWorkingDir() {
		return workingDir;
	}
	/**
	 * @param workingDir the workingDir to set
	 */
	public void setWorkingDir(String workingDir) {
		this.workingDir = workingDir;
	}
	/**
	 * @return the path
	 */
	public String getPath() {
		return path;
	}
	/**
	 * @param path the path to set
	 */
	public void setPath(String path) {
		this.path = path;
	}
	/**
	 * @return the root
	 */
	public DefaultMutableTreeNode getRoot() {
		return root;
	}
	/**
	 * @param root the root to set
	 */
	public void setRoot(DefaultMutableTreeNode root) {
		this.root = root;
	}
	/**
	 * @return the inputstream
	 */
	public InputStream getInputstream() {
		return inputstream;
	}
	/**
	 * @param inputstream the inputstream to set
	 */
	public void setInputstream(InputStream inputstream) {
		this.inputstream = inputstream;
	}
	/**
	 * @return the data
	 */
	public Iterable<Object> getData() {
		return data;
	}
	/**
	 * @param data the data to set
	 */
	public void setData(Iterable<Object> data) {
		this.data = data;
	}
	/**
	 * @return the mapNextLevel0
	 */
	public Map<?, ?> getMapNextLevel0() {
		return mapNextLevel0;
	}
	/**
	 * @param mapNextLevel0 the mapNextLevel0 to set
	 */
	public void setMapNextLevel0(Map<?, ?> mapNextLevel0) {
		this.mapNextLevel0 = mapNextLevel0;
	}

}
