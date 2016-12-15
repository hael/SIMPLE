/**
 * 
 */
package org.openmolgrid.model;

import java.util.ArrayList;

/**
 * This class represents the Khon-Sham Eigenvalue
 * 
 * @author Frederic Bonnet
 */
public class COrbitals implements OrbitalsValues {

	//Global variables
	private ArrayList<Double> e = new ArrayList<Double>();
	private ArrayList<Double> f = new ArrayList<Double>();
	private ArrayList<Integer> k = new ArrayList<Integer>();
	/**
	 * Simple constructor
	 */
	public COrbitals() {
	}
	/**
	 * Constructor for the orbitals for the eigenvalues
	 * @param e
	 */
	public COrbitals(ArrayList<Double> e) {
		this.e = e;
	}
	/**
	 * Constructor for the orbitals for the eigenvalues, occupational Numbers
	 * @param e
	 * @param f
	 */
	public COrbitals(ArrayList<Double> e, ArrayList<Double> f) {
		this.e = e;
		this.f = f;
	}
	/**
	 * Constructor for the orbitals for the eigenvalues, occupational Numbers, k-points
	 * @param e
	 * @param f
	 * @param k
	 */
	public COrbitals(ArrayList<Double> e, ArrayList<Double> f, ArrayList<Integer> k) {
		this.e = e;
		this.f = f;
		this.k = k;		
	}
	public void kohnShamEigenvalues() {
	}
	public void occupationNumbers() {
	}
	public void kptBZCoords() {
	}
	/**
	 * The setters and getters for the:
	 */
	/**
	 * @return the e
	 */
	public ArrayList<Double> getE() {
		return e;
	}
	/**
	 * @param e the e to set
	 */
	public void setE(ArrayList<Double> e) {
		this.e = e;
	}
	/**
	 * @return the f
	 */
	public ArrayList<Double> getF() {
		return f;
	}
	/**
	 * @param f the f to set
	 */
	public void setF(ArrayList<Double> f) {
		this.f = f;
	}
	/**
	 * @return the k
	 */
	public ArrayList<Integer> getK() {
		return k;
	}
	/**
	 * @param k the k to set
	 */
	public void setK(ArrayList<Integer> k) {
		this.k = k;
	}
	@Override
	public String toString() {
		return "The Eigenvalues: "+e+"\n"+
				"The Occupational Numbers: "+f+"\n"+
				"The K-Points: "+k+"\n";
	}
/**End the class COrbitals*/
}
