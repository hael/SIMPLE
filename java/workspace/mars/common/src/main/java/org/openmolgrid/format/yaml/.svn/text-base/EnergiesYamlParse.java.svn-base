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
 *	Class to parse the energies from the
 *	Yaml file
 */
public class EnergiesYamlParse {

	/** Logger */
	private static Logger log = Logger.getLogger(EnergiesYamlParse.class.getName());

	//local variables
	private int whichFermi;
	private int whichHamil_optimi;
	private ArrayList<String> EnergyArrayListName;
	private ArrayList<Double> EnergyArrayListValues;
	//global variables
	private LinkedHashMap<?, ?> energyMap;
	private String energyString;
	private Boolean energyBoolean;
	private double energyDouble;
	private int energyInt;

	private double fermiEnergy;
	private Map<?, ?> hamil_optimi;
	
	/**
	 * simple constructor
	 */
	public EnergiesYamlParse() {
	}
	/**
	 * constructor for the energy Map
	 * @param energyMap:LinkedHashMap<?, ?>
	 */
	public EnergiesYamlParse(LinkedHashMap<?, ?> energyMap) {
		this.energyMap = energyMap;
	}
	/**
	 * constructor for the energy with input parameter energyString
	 * @param energyString:String
	 */
	public EnergiesYamlParse(String energyString) {
		this.energyString = energyString;
	}
	/**
	 * constructor for the EnergiesYamlParse class
	 * @param energyBoolean:Boolean
	 */
	public EnergiesYamlParse(Boolean energyBoolean) {
		this.energyBoolean = energyBoolean;
	}
	/**
	 * constructor for the EnergiesYamlParse class
	 * @param energyDouble:double
	 */
	public EnergiesYamlParse(double energyDouble) {
		this.energyDouble = energyDouble;
	}
	/**
	 * constructor for the EnergiesYamlParse class
	 * @param energyInt:int
	 */
	public EnergiesYamlParse(int energyInt) {
		this.energyInt = energyInt;
	}
	/**
	 * constructor for the class EnergiesYamlParse
	 * @param energyMap:LinkedHashMap<?, ?>
	 * @param energyString:String
	 * @param energyBoolean:Boolean
	 * @param energyDouble:double
	 * @param energyInt:int
	 */
	public EnergiesYamlParse(LinkedHashMap<?, ?> energyMap, String energyString, Boolean energyBoolean, double energyDouble, int energyInt) {
		this.energyMap = energyMap;
		this.energyString = energyString;
		this.energyBoolean = energyBoolean;
		this.energyDouble = energyDouble;
		this.energyInt = energyInt;
	}
	/**
	 * constructor that takes in the Map<?, ?> hamiltonian_optimization 
	 * @param hamiltonian_optimization
	 */
	public EnergiesYamlParse(Map<?, ?> hamiltonian_optimization) {
		this.hamil_optimi = hamiltonian_optimization;
	}
	/**
	 * method to parse the Fermi Energy from the Yamlmap.
	 */
	public void ParseFermiEnergy() throws NullPointerException, ArrayIndexOutOfBoundsException {

		try {
			Collection<?> hamil_optimiSet = ((Map<?, ?>) hamil_optimi).keySet();

			int i = 0;
			Object[] hamil_optimiStringArray = new String[hamil_optimiSet.size()];
			hamil_optimiStringArray = hamil_optimiSet.toArray();
			i = 0;
			whichFermi = -99;
			whichHamil_optimi = 99;
			for ( Object current : hamil_optimiStringArray ) {
				if ( current.equals("Fermi Energy")) {
					whichFermi = i;
					setWhichFermi(whichFermi);
				}
				if ( current.equals("Hamiltonian Optimization")) {
					whichHamil_optimi = i;
					setWhichHamil_optimi(whichHamil_optimi);
				}
				i++;
			}
			/**The Fermi energy parsing*/
			try {
				fermiEnergy = (Double)hamil_optimi.get(hamil_optimiStringArray[whichFermi]);
				setFermiEnergy(fermiEnergy);
			} catch (Exception e) {
				log.info("The array is out of bound: "+whichFermi);
				log.info("The file does not contain the entry --> Fermi Energy:"); 
				log.info("The exception printout: "+e.toString());
			}
		} catch (Exception e) {
			log.info("Handler for the catch and the exception in testGroundStateOptimization");
			log.info("The exception printout: "+e.toString());
		}
	}
	/**
	 * Method to parse the energy values from the Yamlmap.
	 */
	public void ParseEnergies() throws ArrayIndexOutOfBoundsException {
		
		try {
			Collection<?> hamil_optimiSet = ((Map<?, ?>) hamil_optimi).keySet();

			int i = 0;
			Object[] hamil_optimiStringArray = new String[hamil_optimiSet.size()];
			hamil_optimiStringArray = hamil_optimiSet.toArray();
			i = 0;
			whichFermi = -99;
			whichHamil_optimi = 99;
			for ( Object current : hamil_optimiStringArray ) {
				if ( current.equals("Fermi Energy")) {
					whichFermi = i;
					setWhichFermi(whichFermi);
				}
				if ( current.equals("Hamiltonian Optimization")) {
					whichHamil_optimi = i;
					setWhichHamil_optimi(whichHamil_optimi);
				}
				i++;
			}
			/**The energies parsing*/
			String searchString="Energies";
			GroundStateOptimizationParser groundStateOptimizationParser =
					new GroundStateOptimizationParser(
							searchString,
							hamil_optimi,
							hamil_optimiStringArray,
							whichHamil_optimi);
			
			groundStateOptimizationParser.hamiltonianOptimizationParser();
			
			if ( groundStateOptimizationParser.getEnDouble() != 0.0 ) {
				log.info("This is an instance of (Double)");
				log.info(groundStateOptimizationParser.getEnDouble());
//				System.out.println(energyDouble);
			}else if ( !groundStateOptimizationParser.getEnMap().equals(null) ) {
				log.info("This is an instance of (Map<?,?>)");
				log.info(groundStateOptimizationParser.getEnMap());
//				System.out.println(energyMap);
			}else if ( groundStateOptimizationParser.getEnString() != null ) {
				log.info("This is an instance of (String)");
				log.info(groundStateOptimizationParser.getEnString());
//				System.out.println(energyString);
			}else if ( groundStateOptimizationParser.getEnInt() != -99 ) {
				log.info("This is an instance of (Integer)");
				log.info(groundStateOptimizationParser.getEnInt());
//                System.out.println(energyInt);
			}
			
			/**
			 * Adding the values to the
			 * ArrayLists{<String>EnergyArrayListName,<Double>EnergyArrayListValues}
			 */
			ArrayList<String> EnergyArrayListName = new ArrayList<String>();
			ArrayList<Double> EnergyArrayListValues = new ArrayList<Double>();
			Collection<?> energyCollection = ((LinkedHashMap<?,?>)groundStateOptimizationParser.getEnMap()).keySet();
			int n = 0;
			for (Object currentEnergy : energyCollection ) {
				//System.out.println("for n: "+n+" "+currentEnergy);
				EnergyArrayListName.add((String)currentEnergy);
				EnergyArrayListValues.add((Double)groundStateOptimizationParser.getEnMap().get(currentEnergy));
				//System.out.println(EnergyArrayListName.get(n)+"-->"+EnergyArrayListValues.get(n));
				n++;
			}
			/**
			 * setting up the values in
			 * ArrayLists{<String>EnergyArrayListName,<Double>EnergyArrayListValues}
			 */
			setEnergyArrayListName(EnergyArrayListName);
			setEnergyArrayListValues(EnergyArrayListValues);
			
		} catch (Exception e) {
			log.info("Handler for the catch and the exception in testGroundStateOptimization");
			log.info("The exception printout: "+e.toString());
		}
		
	}
	/**
	 * The setters and getters for the:
	 * 	Fermi Energy                  : fermiEnergy
	 *	private LinkedHashMap<?, ?> energyMap;
	 *	private String energyString;
	 *	private Boolean energyBoolean;
	 *	private double energyDouble;
	 *	private int energyInt;
	 */
	/**
	 * @return the energyMap
	 */
	public LinkedHashMap<?, ?> getEnergyMap() {
		return energyMap;
	}
	/**
	 * @param energyMap the energyMap to set
	 */
	public void setEnergyMap(LinkedHashMap<?, ?> energyMap) {
		this.energyMap = energyMap;
	}
	/**
	 * @return the energyString
	 */
	public String getEnergyString() {
		return energyString;
	}
	/**
	 * @param energyString the energyString to set
	 */
	public void setEnergyString(String energyString) {
		this.energyString = energyString;
	}
	/**
	 * @return the energyBoolean
	 */
	public Boolean getEnergyBoolean() {
		return energyBoolean;
	}
	/**
	 * @param energyBoolean the energyBoolean to set
	 */
	public void setEnergyBoolean(Boolean energyBoolean) {
		this.energyBoolean = energyBoolean;
	}
	/**
	 * @return the energyDouble
	 */
	public double getEnergyDouble() {
		return energyDouble;
	}
	/**
	 * @param energyDouble the energyDouble to set
	 */
	public void setEnergyDouble(double energyDouble) {
		this.energyDouble = energyDouble;
	}
	/**
	 * @return the energyInt
	 */
	public int getEnergyInt() {
		return energyInt;
	}
	/**
	 * @param energyInt the energyInt to set
	 */
	public void setEnergyInt(int energyInt) {
		this.energyInt = energyInt;
	}
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
	 * @return the hamil_optimi
	 */
	public Map<?, ?> getHamil_optimi() {
		return hamil_optimi;
	}
	/**
	 * @param hamil_optimi the hamil_optimi to set
	 */
	public void setHamil_optimi(Map<?, ?> hamil_optimi) {
		this.hamil_optimi = hamil_optimi;
	}
	/**
	 * @return the whichFermi
	 */
	public int getWhichFermi() {
		return whichFermi;
	}
	/**
	 * @param whichFermi the whichFermi to set
	 */
	public void setWhichFermi(int whichFermi) {
		this.whichFermi = whichFermi;
	}
	/**
	 * @return the whichHamil_optimi
	 */
	public int getWhichHamil_optimi() {
		return whichHamil_optimi;
	}
	/**
	 * @param whichHamil_optimi the whichHamil_optimi to set
	 */
	public void setWhichHamil_optimi(int whichHamil_optimi) {
		this.whichHamil_optimi = whichHamil_optimi;
	}
	/**
	 * @return the energyArrayListName
	 */
	public ArrayList<String> getEnergyArrayListName() {
		return EnergyArrayListName;
	}
	/**
	 * @param energyArrayListName the energyArrayListName to set
	 */
	public void setEnergyArrayListName(ArrayList<String> energyArrayListName) {
		EnergyArrayListName = energyArrayListName;
	}
	/**
	 * @return the energyArrayListValues
	 */
	public ArrayList<Double> getEnergyArrayListValues() {
		return EnergyArrayListValues;
	}
	/**
	 * @param energyArrayListValues the energyArrayListValues to set
	 */
	public void setEnergyArrayListValues(ArrayList<Double> energyArrayListValues) {
		EnergyArrayListValues = energyArrayListValues;
	}
}
