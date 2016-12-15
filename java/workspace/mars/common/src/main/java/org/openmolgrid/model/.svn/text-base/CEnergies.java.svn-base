/**
 * 
 */
package org.openmolgrid.model;

import java.util.ArrayList;

import org.openmolgrid.format.cml.CMLWriter;

/**
 * @author Frederic Bonnet
 *
 */
public class CEnergies implements EnergiesValues {

	//Global variables
	private double fermiEnergy;
	private ArrayList<?> energyArrayList;
	private ArrayList<String> energyArrayListName;
	private ArrayList<Double> energyArrayListValues;
	/**
	 * Simple constructor
	 */
	public CEnergies() {
	}	
	public CEnergies(double fermiEnergy) {
		this.fermiEnergy  = fermiEnergy;
	}
	public CEnergies(ArrayList<?> energyArrayList) {
		this.energyArrayList = energyArrayList;
	}
	public CEnergies(double fermiEnergy, ArrayList<?> energyArrayList) {
		this.fermiEnergy  = fermiEnergy;
		this.energyArrayList = energyArrayList;
	}
	public CEnergies(double fermiEnergy, ArrayList<String> energyArrayListName, ArrayList<Double> energyArrayListValues ) {
		this.fermiEnergy  = fermiEnergy;
		this.energyArrayListName = energyArrayListName;
		this.energyArrayListValues = energyArrayListValues;
	}
	public CEnergies(ArrayList<String> energyArrayListName, ArrayList<Double> energyArrayListValues ) {
		this.energyArrayListName = energyArrayListName;
		this.energyArrayListValues = energyArrayListValues;
	}
	/**
	 * Writes a cml representation of this structure as string
	 * 
	 * @return a cml representation of this structure as string
	 */
//	public String writeCML2String() {
//		return CMLWriter.asXmlString(this);
//	}
	/**
	 * The setters and getters for the:
	 * 	Fermi Energy                  : fermiEnergy
	 * 	Arbitrary energy value private ArrayList<?> energyArrayList;
	 *  private ArrayList<String> energyArrayListName;
	 *  private ArrayList<Double> energyArrayListValues;
	 */
	/**
	 * @return the fermiEnergy
	 */
	public double getFermiEnergy() {
		return fermiEnergy;
	}
	/**
	 * @param fermiEnergy the fermiEnergy to set
	 */
	public void setFermiEnergy(double fermiEnergy) {
		this.fermiEnergy = fermiEnergy;
	}
	/**
	 * @return the energyArrayList
	 */
	public ArrayList<?> getEnergyArrayList() {
		return energyArrayList;
	}
	/**
	 * @param energyArrayList the energyArrayList to set
	 */
	public void setEnergyArrayList(ArrayList<?> energyArrayList) {
		this.energyArrayList = energyArrayList;
	}
	/**
	 * @return the energyArrayListName
	 */
	public ArrayList<String> getEnergyArrayListName() {
		return energyArrayListName;
	}
	/**
	 * @param energyArrayListName the energyArrayListName to set
	 */
	public void setEnergyArrayListName(ArrayList<String> energyArrayListName) {
		this.energyArrayListName = energyArrayListName;
	}
	/**
	 * @return the energyArrayListValues
	 */
	public ArrayList<Double> getEnergyArrayListValues() {
		return energyArrayListValues;
	}
	/**
	 * @param energyArrayListValues the energyArrayListValues to set
	 */
	public void setEnergyArrayListValues(ArrayList<Double> energyArrayListValues) {
		this.energyArrayListValues = energyArrayListValues;
	}
}
