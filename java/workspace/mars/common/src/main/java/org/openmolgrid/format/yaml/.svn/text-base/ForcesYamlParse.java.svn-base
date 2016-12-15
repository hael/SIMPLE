/**
 * 
 */
package org.openmolgrid.format.yaml;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map;

import org.apache.log4j.Logger;

/**
 * @author Frederic Bonnet
 * 
 * Class to parse the forces from the Yaml
 * output
 *
 */
public class ForcesYamlParse {
	/** Logger */
	private static Logger log = Logger.getLogger(ForcesYamlParse.class.getName());

	//local variables
	private int nstressTensor;// = 0;
	private int nints;// = 0;
	private int nmaps;// = 0;
	private int nSymmetries;// = 0;
	private ArrayList<?> stressTensorMatrix;// = null;
	private Collection<?> pressureSet;// = null;
	private ArrayList<Double> pressureArrayList;// = null;
	
	//global variables
	private LinkedHashMap<?, ?> stressTensorsLinkedHashMap;
	private Collection<?> stressTensorSet;
	/**
	 * simple constructor
	 */
	public ForcesYamlParse() {
	}
	/**
	 * constructor of the class.
	 * 
	 * @param stressTensorsLinkedHashMap
	 * @param stressTensorSet
	 */
	public ForcesYamlParse(LinkedHashMap<?, ?> stressTensorsLinkedHashMap, Collection<?> stressTensorSet) {
		this.stressTensorsLinkedHashMap = stressTensorsLinkedHashMap;
		this.stressTensorSet = stressTensorSet;
	}
	/**
	 * method to parse the content of the file.
	 */
	public void ParseForces() {
	
		for (Object currentStressTensorNameString : stressTensorSet ) {
			//System.out.println(currentStressTensorNameString);
			if ( stressTensorsLinkedHashMap.get(currentStressTensorNameString) instanceof ArrayList<?> ) {
				stressTensorMatrix = (ArrayList<?>)stressTensorsLinkedHashMap.get(currentStressTensorNameString);
				//System.out.println(stressTensorMatrix);
				nstressTensor++;
				setNstressTensor(nstressTensor);
				setStressTensorMatrix(stressTensorMatrix);
			} else if (stressTensorsLinkedHashMap.get(currentStressTensorNameString) instanceof Integer) {
				nSymmetries = (Integer)stressTensorsLinkedHashMap.get(currentStressTensorNameString);
				//System.out.println(stressTensorInt);
				nints++;
				setnSymmetries(nSymmetries);
				setNints(nints);
			} else if (stressTensorsLinkedHashMap.get(currentStressTensorNameString) instanceof Map<?, ?>) {
				Map<?, ?> stressTensorMap = (Map<?, ?>)stressTensorsLinkedHashMap.get(currentStressTensorNameString);
				//System.out.println(stressTensorMap);
				pressureSet = ((Map<?,?>)stressTensorMap).keySet();
				pressureArrayList = new ArrayList<Double>();
				int ipressure = 0;
				for ( Object currentPressureStringname : pressureSet) {
					pressureArrayList.add((Double)stressTensorMap.get(currentPressureStringname));
					//System.out.println(pressureArrayList.get(ipressure)+"\n");
					ipressure++;
				}
				log.info("The number of pressure entries: "+ipressure);
				nmaps++;
				setPressureSet(pressureSet);
				setPressureArrayList(pressureArrayList);
			}
		}
		setnSymmetries(nSymmetries);
		setNints(nints);
		setNmaps(nmaps);
		setStressTensorMatrix(stressTensorMatrix);
		setPressureSet(pressureSet);
		setPressureArrayList(pressureArrayList);
	}
	/**
	 * The setters and getters for the:
	 * 	int nstressTensor = 0;
	 *  int nints = 0;
	 *  private int nmaps = 0;
	 *  private int nSymmetries = 0;
	 *  private ArrayList<?> stressTensorMatrix = null;
	 *  private Collection<?> pressureSet = null;
	 *  private ArrayList<Double> pressureArrayList = null;
	 *  //global variables
	 *  private LinkedHashMap<?, ?> stressTensorsLinkedHashMap;
	 *  private Collection<?> stressTensorSet;
	 * 
	 */
	/**
	 * @return the nstressTensor
	 */
	public int getNstressTensor() {
		return nstressTensor;
	}
	/**
	 * @param nstressTensor the nstressTensor to set
	 */
	public void setNstressTensor(int nstressTensor) {
		this.nstressTensor = nstressTensor;
	}
	/**
	 * @return the nints
	 */
	public int getNints() {
		return nints;
	}
	/**
	 * @param nints the nints to set
	 */
	public void setNints(int nints) {
		this.nints = nints;
	}
	/**
	 * @return the nmaps
	 */
	public int getNmaps() {
		return nmaps;
	}
	/**
	 * @param nmaps the nmaps to set
	 */
	public void setNmaps(int nmaps) {
		this.nmaps = nmaps;
	}
	/**
	 * @return the nSymmetries
	 */
	public int getnSymmetries() {
		return nSymmetries;
	}
	/**
	 * @param nSymmetries the nSymmetries to set
	 */
	public void setnSymmetries(int nSymmetries) {
		this.nSymmetries = nSymmetries;
	}
	/**
	 * @return the stressTensorMatrix
	 */
	public ArrayList<?> getStressTensorMatrix() {
		return stressTensorMatrix;
	}
	/**
	 * @param stressTensorMatrix the stressTensorMatrix to set
	 */
	public void setStressTensorMatrix(ArrayList<?> stressTensorMatrix) {
		this.stressTensorMatrix = stressTensorMatrix;
	}
	/**
	 * @return the pressureSet
	 */
	public Collection<?> getPressureSet() {
		return pressureSet;
	}
	/**
	 * @param pressureSet the pressureSet to set
	 */
	public void setPressureSet(Collection<?> pressureSet) {
		this.pressureSet = pressureSet;
	}
	/**
	 * @return the pressureArrayList
	 */
	public ArrayList<Double> getPressureArrayList() {
		return pressureArrayList;
	}
	/**
	 * @param pressureArrayList the pressureArrayList to set
	 */
	public void setPressureArrayList(ArrayList<Double> pressureArrayList) {
		this.pressureArrayList = pressureArrayList;
	}
	/**
	 * @return the stressTensorsLinkedHashMap
	 */
	public LinkedHashMap<?, ?> getStressTensorsLinkedHashMap() {
		return stressTensorsLinkedHashMap;
	}
	/**
	 * @param stressTensorsLinkedHashMap the stressTensorsLinkedHashMap to set
	 */
	public void setStressTensorsLinkedHashMap(
			LinkedHashMap<?, ?> stressTensorsLinkedHashMap) {
		this.stressTensorsLinkedHashMap = stressTensorsLinkedHashMap;
	}
	/**
	 * @return the stressTensorSet
	 */
	public Collection<?> getStressTensorSet() {
		return stressTensorSet;
	}
	/**
	 * @param stressTensorSet the stressTensorSet to set
	 */
	public void setStressTensorSet(Collection<?> stressTensorSet) {
		this.stressTensorSet = stressTensorSet;
	}


}
