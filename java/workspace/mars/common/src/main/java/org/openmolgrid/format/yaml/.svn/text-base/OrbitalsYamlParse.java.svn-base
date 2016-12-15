/**
 * 
 */
package org.openmolgrid.format.yaml;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * @author Frederic Bonnet
 *
 */
public class OrbitalsYamlParse {

	//local variables
	private ArrayList<Double> arrayList_e = new ArrayList<Double>();
	private ArrayList<Double> arrayList_f = new ArrayList<Double>();
	private ArrayList<Integer> arrayList_k = new ArrayList<Integer>();
	//global variables
	private ArrayList<?> orbitals;
	/**
	 * Simple constructor
	 */
	public OrbitalsYamlParse() {
		
	}
	/**
	 * Constructor for the class
	 * 
	 * @param orbitals
	 */
	public OrbitalsYamlParse(ArrayList<?> orbitals) {
		this.orbitals = orbitals;
	}
	/**
	 * Method to parse the content for the orbitals
	 */
	public void ParseOrbitals() {
		
		System.out.println(orbitals.get(0));
		Map<?,?> items = (Map<?,?>)orbitals.get(0);
		System.out.println(items.get("e"));
		int index = 0;
		for (Object iorbitals : orbitals) {
			Collection<?> curentCollection = ((LinkedHashMap<?,?>)iorbitals).keySet();
			//System.out.println(curentCollection+" "+curentCollection.toArray()[0]+" "+curentCollection.toArray()[1]+" "+curentCollection.toArray()[2]);
			items = (Map<?,?>)orbitals.get(index);
			
//			System.out.println(
//					curentCollection.toArray()[0]+": "+items.get(curentCollection.toArray()[0])+" "
//			+		curentCollection.toArray()[1]+": "+items.get(curentCollection.toArray()[1])+" "
//			+		curentCollection.toArray()[2]+": "+items.get(curentCollection.toArray()[2]) );

			double e = (Double)items.get(curentCollection.toArray()[0]);
			double f = (Double)items.get(curentCollection.toArray()[1]);
			int k = (Integer)items.get(curentCollection.toArray()[2]);

			arrayList_e.add(index, e);
			arrayList_f.add(index, f);
			arrayList_k.add(index, k);

			setArrayList_e(arrayList_e);
			setArrayList_f(arrayList_f);
			setArrayList_k(arrayList_k);
			
			index++;
		}
	}
	/**
	 * The setters and getters for the:
	 * 	Orbitals                      : nlocForces
	 *  arrayList_{e,f,k}             : 
	 */
	/**
	 * @return the orbitals
	 */
	public ArrayList<?> getOrbitals() {
		return orbitals;
	}
	/**
	 * @param orbitals the orbitals to set
	 */
	public void setOrbitals(ArrayList<?> orbitals) {
		this.orbitals = orbitals;
	}
	/**
	 * @return the arrayList_e
	 */
	public ArrayList<Double> getArrayList_e() {
		return arrayList_e;
	}
	/**
	 * @param arrayList_e the arrayList_e to set
	 */
	public void setArrayList_e(ArrayList<Double> arrayList_e) {
		this.arrayList_e = arrayList_e;
	}
	/**
	 * @return the arrayList_f
	 */
	public ArrayList<Double> getArrayList_f() {
		return arrayList_f;
	}
	/**
	 * @param arrayList_f the arrayList_f to set
	 */
	public void setArrayList_f(ArrayList<Double> arrayList_f) {
		this.arrayList_f = arrayList_f;
	}
	/**
	 * @return the arrayList_k
	 */
	public ArrayList<Integer> getArrayList_k() {
		return arrayList_k;
	}
	/**
	 * @param arrayList_k the arrayList_k to set
	 */
	public void setArrayList_k(ArrayList<Integer> arrayList_k) {
		this.arrayList_k = arrayList_k;
	}
	
	
}
