/**
 * 
 */
package org.openmolgrid.format.yaml;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map;

import org.yaml.snakeyaml.Yaml;

/**
 * @author Frederic Bonnet
 *
 *	Class to parse the Ground State optimization map
 */
public class GroundStateOptimizationParser {
	
	//local variables
	private int enInt;
	private String enString;
	private Boolean enBoolean;
	private double enDouble;
	private int whichEs;
	private LinkedHashMap<?, ?> enMap;
	
	//global variables
	private String searchString;
	private Map<?, ?> hamil_optimi;
	private Object[] hamil_optimiStringArray;
	private int whichHamil_optimi;
	/**
	 * Simple constructor
	 */
	public GroundStateOptimizationParser() {
	}
	/**
	 * Constructor for the class taking in only the Map<?,?>
	 * 
	 * @param hamiltonian_optimization
	 */
	public GroundStateOptimizationParser(Map<?, ?> hamiltonian_optimization) {
		this.hamil_optimi = hamiltonian_optimization;
	}
	/**
	 * Constructor for the class taking in the Map<?,?> and the
	 * String.
	 * 
	 * @param searchString
	 * @param hamil_optimi
	 */
	public GroundStateOptimizationParser(String searchString,
			Map<?, ?> hamil_optimi,
			Object[] hamil_optimiStringArray,
			int whichHamil_optimi) {
		this.searchString = searchString;
		this.hamil_optimi = hamil_optimi;
		this.hamil_optimiStringArray = hamil_optimiStringArray;
		this.whichHamil_optimi = whichHamil_optimi;
		}
	/**
	 * 
	 */
	public void hamiltonianOptimizationParser() {
	
		Yaml yaml = new Yaml();
		
		enInt = -999;
		enString = null;
		enBoolean = false;
		enDouble = 0.0;
		whichEs = -99;
		enMap = null;
		
		ArrayList<?> hamil_optimiArrayList = (ArrayList<?>)hamil_optimi.get(hamil_optimiStringArray[whichHamil_optimi]);
//		System.out.println(yaml.dump(hamil_optimiArrayList));

		int i = 0;
		for ( Object crtHamil_optimiList : hamil_optimiArrayList ) {
			//System.out.println("for i: "+i+" -->"+crtHamil_optimiList+"<--");
			LinkedHashMap<?, ?> subopt = (LinkedHashMap<?,?>)crtHamil_optimiList;
			//System.out.println("for i: "+i+" ++>"+subopt+"<++");
			Collection<?> suboptCollection = ((LinkedHashMap<?,?>)subopt).keySet();
			//Collection<?> suboptCollection = ((LinkedHashMap<?,?>)((LinkedHashMap<?,?>)crtHamil_optimiList)).keySet();
			int j = 0;
			for ( Object crtsuboptCollection : suboptCollection) {
				//System.out.println("for j: "+j+" "+crtsuboptCollection);
				LinkedHashMap<?, ?> wavefunc = (LinkedHashMap<?,?>)subopt.get(crtsuboptCollection);
				//System.out.println(wavefunc);
				Collection<?> waveCollection = ((LinkedHashMap<?,?>)wavefunc).keySet();
				int k = 0;
				for ( Object crtwaveCollection : waveCollection ) {
					//System.out.println("for k: "+k+" "+crtwaveCollection);
					ArrayList<?> wavefuncArrayList = (ArrayList<?>)wavefunc.get(crtwaveCollection);
					//System.out.println(arraylist);
					int l = 0;
					for ( Object crtWavefuncArrayList : wavefuncArrayList ) {
						//System.out.println("for l: "+l+" "+crtWavefuncArrayList);
						LinkedHashMap<?,?> map = (LinkedHashMap<?, ?>)wavefuncArrayList.get(l);
						//System.out.println(yaml.dump(map));
						Collection<?> collection = ((LinkedHashMap<?,?>)map).keySet();
						int m = 0;
						for (Object crtCollection : collection ) {
							whichEs = -99;
							//System.out.println("for m: "+m+" "+crtCollection);
							if ( crtCollection.equals(searchString ) ){
								whichEs = m;
								//System.out.println("In if search whichEnergies: "+whichEnergies);
								//System.out.println(map.get(cureent));
								//System.out.println("Before the if statements");
								/**The instances of*/
								if ( map.get(crtCollection) instanceof LinkedHashMap<?, ?> ) {
									enMap = (LinkedHashMap<?, ?>)map.get(crtCollection);
//									System.out.println(yaml.dump(enMap));
									setEnMap(enMap); break;
								}else if ( map.get(crtCollection) instanceof Double ) {
									enDouble = (Double)map.get(crtCollection);
									setEnDouble(enDouble); break;
								}else if (map.get(crtCollection) instanceof String ) {
									enString = (String)map.get(crtCollection);
									setEnString(enString); break;
								}else if (map.get(crtCollection) instanceof Boolean ) {
									enBoolean = (Boolean)map.get(crtCollection);
									setEnBoolean(enBoolean); break;
								}else if (map.get(crtCollection) instanceof Integer ) {
									enInt = (Integer)map.get(crtCollection);
									setEnInt(enInt);break;
								}
								/**----------------*/
								//System.out.println("After the if statements");
							}
							m++;
						}
						//System.out.println("location of the Energy entry whichEnergies: "+whichEnergies);
						l++;
					}
					k++;
				}
				j++;
			}
			i++;
		}
	}
	/**
	 * The setters and getters for the class:
	 * private String searchString;
	 * private Map<?, ?> hamiltonian_optimization;
	 */
	/**
	 * @return the enInt
	 */
	public int getEnInt() {
		return enInt;
	}
	/**
	 * @param enInt the enInt to set
	 */
	public void setEnInt(int enInt) {
		this.enInt = enInt;
	}
	/**
	 * @return the enString
	 */
	public String getEnString() {
		return enString;
	}
	/**
	 * @param enString the enString to set
	 */
	public void setEnString(String enString) {
		this.enString = enString;
	}
	/**
	 * @return the enBoolean
	 */
	public Boolean getEnBoolean() {
		return enBoolean;
	}
	/**
	 * @param enBoolean the enBoolean to set
	 */
	public void setEnBoolean(Boolean enBoolean) {
		this.enBoolean = enBoolean;
	}
	/**
	 * @return the enDouble
	 */
	public double getEnDouble() {
		return enDouble;
	}
	/**
	 * @param enDouble the enDouble to set
	 */
	public void setEnDouble(double enDouble) {
		this.enDouble = enDouble;
	}
	/**
	 * @return the whichEs
	 */
	public int getWhichEs() {
		return whichEs;
	}
	/**
	 * @param whichEs the whichEs to set
	 */
	public void setWhichEs(int whichEs) {
		this.whichEs = whichEs;
	}
	/**
	 * @return the enMap
	 */
	public LinkedHashMap<?, ?> getEnMap() {
		return enMap;
	}
	/**
	 * @param enMap the enMap to set
	 */
	public void setEnMap(LinkedHashMap<?, ?> enMap) {
		this.enMap = enMap;
	}
	/**
	 * @return the searchString
	 */
	public String getSearchString() {
		return searchString;
	}
	/**
	 * @param searchString the searchString to set
	 */
	public void setSearchString(String searchString) {
		this.searchString = searchString;
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
	 * @return the hamil_optimiStringArray
	 */
	public Object[] getHamil_optimiStringArray() {
		return hamil_optimiStringArray;
	}
	/**
	 * @param hamil_optimiStringArray the hamil_optimiStringArray to set
	 */
	public void setHamil_optimiStringArray(Object[] hamil_optimiStringArray) {
		this.hamil_optimiStringArray = hamil_optimiStringArray;
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
	
}
