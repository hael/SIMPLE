/**
 * 
 */
package org.openmolgrid.model;

import java.util.ArrayList;
import java.util.Collection;

/**
 * @author Frederic Bonnet, Date: 9th of August 2012
 *
 */
public class CForces implements ForcesValues {

	//Global variables
	private String nlocForces;
	private int nSymmetries;
	private ArrayList<?> stressTensorMatrix;
	private Collection<?> pressureSet;
	private ArrayList<Double> pressureArrayList;
	//local variables
	private double pressureHaBohr3;
	private double pressureGPa;
	private double pressurePV_Ha;	
	/**
	 * Simple constructor
	 */
	public CForces() {
	}
	/**
	 * The constructor for the forces
	 * 
	 * @param nlocForces
	 * @param stressTensorMatrix
	 */
	public CForces(String nlocForces, int nSymmetries,
			ArrayList<?> stressTensorMatrix, Collection<?> pressureSet, ArrayList<Double> pressureArrayList) {
		this.nlocForces = nlocForces;
		this.nSymmetries = nSymmetries;
		this.stressTensorMatrix = stressTensorMatrix;
		this.pressureSet = pressureSet;
		this.pressureArrayList = pressureArrayList;
	}

	public void separatePressure() {

		for ( Object currentpressureSet : pressureSet ) {
			if ( currentpressureSet.equals("Ha/Bohr^3") ) {
				setPressureHaBohr3(pressureArrayList.get(0));
			}else if (currentpressureSet.equals("GPa") ) {
				setPressureGPa(pressureArrayList.get(1));
			}else if (currentpressureSet.equals("PV (Ha)") ) {
				setPressurePV_Ha(pressureArrayList.get(2));
			}
		}
	}
	
	/**
	 * The setters and getters for the:
	 * 	Non Local forces calculated   : nlocForces
	 *  Stress Tensor                 : stressTensorMatrix
	 */

	/**
	 * @return the nlocForces
	 */
	public String getNlocForces() {
		return nlocForces;
	}
	/**
	 * @param nlocForces the nlocForces to set
	 */
	public void setNlocForces(String nlocForces) {
		this.nlocForces = nlocForces;
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
	 * @return the pressureHaBohr3
	 */
	public double getPressureHaBohr3() {
		return pressureHaBohr3;
	}
	/**
	 * @param pressureHaBohr3 the pressureHaBohr3 to set
	 */
	public void setPressureHaBohr3(double pressureHaBohr3) {
		this.pressureHaBohr3 = pressureHaBohr3;
	}
	/**
	 * @return the pressureGPa
	 */
	public double getPressureGPa() {
		return pressureGPa;
	}
	/**
	 * @param pressureGPa the pressureGPa to set
	 */
	public void setPressureGPa(double pressureGPa) {
		this.pressureGPa = pressureGPa;
	}
	/**
	 * @return the pressurePV_Ha
	 */
	public double getPressurePV_Ha() {
		return pressurePV_Ha;
	}
	/**
	 * @param pressurePV_Ha the pressurePV_Ha to set
	 */
	public void setPressurePV_Ha(double pressurePV_Ha) {
		this.pressurePV_Ha = pressurePV_Ha;
	}
	

	
}
