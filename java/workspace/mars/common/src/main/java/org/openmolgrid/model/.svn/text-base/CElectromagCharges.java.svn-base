/**
 * 
 */
package org.openmolgrid.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * This class represents the values for the
 * 	Total electronic charge
 * 	Electric Dipole moment (AU)
 * 	Electric Dipole moment (Debye)
 * 	Partial charges //TODO: add the partial charges waiting on yaml output
 * 
 * @author Frederic Bonnet
 *
 */
public class CElectromagCharges implements ElectromagChargesValues {

	//Global variables
	private double totalElectronicQ;
	private Map<?, ?> eDipoleMoment_AU, eDipoleMoment_Debye;
	//local variables
	private ArrayList<?> pVector;
	private double px,py,pz;
	private double pnorm;
	private double pnormCalulated;
	/**
	 * Simple constructor
	 */
	public CElectromagCharges() {
	}
	/**
	 * Constructor for the class CElectromagCharges
	 *  
	 * @param totalElectronicQ
	 * @param eDipoleMoment_AU
	 * @param eDipoleMoment_Debye
	 */
	public CElectromagCharges(double totalElectronicQ, Map<?,?> eDipoleMoment_AU,  Map<?,?> eDipoleMoment_Debye ) {
		this.totalElectronicQ = totalElectronicQ;
		this.eDipoleMoment_AU = eDipoleMoment_AU;
		this.eDipoleMoment_Debye = eDipoleMoment_Debye;
	}
	/**
	 * 
	 */
	public void getPVector_AU() {
		getPVector(eDipoleMoment_AU);
		PNormCalculator();
	}
	public void getPNorm_AU() {
		getPNorm(eDipoleMoment_AU);
	}
	public void getPVector_Debye() {
		getPVector(eDipoleMoment_Debye);
		PNormCalculator();
	}
	public void getPNorm_Debye() {
		getPNorm(eDipoleMoment_Debye);
	}
	/**
	 * extracting the vector eDipoleMoment_AU map
	 */
	public void getPVector(Map<?, ?> inMap) {
		Set<?> keyset = (Set<?>)inMap.keySet();
		int nkeys = keyset.size();
		//System.out.println("The number of keysets: "+nkeys);
		Object[] arrayofkeysets = new String[nkeys];
		arrayofkeysets = keyset.toArray();
		for (int ikey = 0 ; ikey < nkeys ; ikey++) {
			if ( inMap.get(arrayofkeysets[ikey]) instanceof ArrayList<?> ) {
				pVector = (ArrayList<?>)inMap.get(arrayofkeysets[ikey]);
				setpVector(pVector);
				setPx((Double)pVector.get(0));
				setPy((Double)pVector.get(1));
				setPz((Double)pVector.get(2));
			}
		}
	}
	/**
	 * extracting the norm of the vector 
	 */
	public void getPNorm(Map<?, ?> inMap) {
		Set<?> keyset = (Set<?>)inMap.keySet();
		int nkeys = keyset.size();
		//System.out.println("The number of keysets: "+nkeys);
		Object[] arrayofkeysets = new String[nkeys];
		arrayofkeysets = keyset.toArray();
		for (int ikey = 0 ; ikey < nkeys ; ikey++) {
			if ( inMap.get(arrayofkeysets[ikey]) instanceof Double ) {
				setPnorm((Double)inMap.get(arrayofkeysets[ikey]));
			}
		}
	}
	/**
	 * norm calculator of the vector P 	
	 */
	public void PNormCalculator() {
		pnormCalulated = Math.sqrt(getPx()*getPx() + getPy()*getPy() + getPz()*getPz() );
		setPnormCalulated(pnormCalulated);
	}
	/**
	 * The setters and getters for the:
	 * 	Total electronic charge:      : totalElectronicQ
	 * 	Electric Dipole Moment (AU)   : eDipoleMoment_AU
	 *  Electric Dipole Moment (Debye): eDipoleMoment_Debye
	 */

	/**
	 * @return the totalElectronicQ
	 */
	public double getTotalElectronicQ() {
		return totalElectronicQ;
	}
	/**
	 * @param totalElectronicQ the totalElectronicQ to set
	 */
	public void setTotalElectronicQ(double totalElectronicQ) {
		this.totalElectronicQ = totalElectronicQ;
	}
	/**
	 * @return the eDipoleMoment_AU
	 */
	public Map<?, ?> geteDipoleMoment_AU() {
		return eDipoleMoment_AU;
	}
	/**
	 * @param eDipoleMoment_AU the eDipoleMoment_AU to set
	 */
	public void seteDipoleMoment_AU(Map<?, ?> eDipoleMoment_AU) {
		this.eDipoleMoment_AU = eDipoleMoment_AU;
	}
	/**
	 * @return the eDipoleMoment_Debye
	 */
	public Map<?, ?> geteDipoleMoment_Debye() {
		return eDipoleMoment_Debye;
	}
	/**
	 * @param eDipoleMoment_Debye the eDipoleMoment_Debye to set
	 */
	public void seteDipoleMoment_Debye(Map<?, ?> eDipoleMoment_Debye) {
		this.eDipoleMoment_Debye = eDipoleMoment_Debye;
	}
	/**
	 * @return the pVector
	 */
	public ArrayList<?> getpVector() {
		return pVector;
	}
	/**
	 * @param pVector the pVector to set
	 */
	public void setpVector(ArrayList<?> pVector) {
		this.pVector = pVector;
	}
	/**
	 * @return the px
	 */
	public double getPx() {
		return px;
	}
	/**
	 * @param px the px to set
	 */
	public void setPx(double px) {
		this.px = px;
	}
	/**
	 * @return the py
	 */
	public double getPy() {
		return py;
	}
	/**
	 * @param py the py to set
	 */
	public void setPy(double py) {
		this.py = py;
	}
	/**
	 * @return the pz
	 */
	public double getPz() {
		return pz;
	}
	/**
	 * @param pz the pz to set
	 */
	public void setPz(double pz) {
		this.pz = pz;
	}
	/**
	 * @return the pnorm
	 */
	public double getPnorm() {
		return pnorm;
	}
	/**
	 * @param pnorm the pnorm to set
	 */
	public void setPnorm(double pnorm) {
		this.pnorm = pnorm;
	}
	/**
	 * @return the pnormCalulated
	 */
	public double getPnormCalulated() {
		return pnormCalulated;
	}
	/**
	 * @param pnormCalulated the pnormCalulated to set
	 */
	public void setPnormCalulated(double pnormCalulated) {
		this.pnormCalulated = pnormCalulated;
	}
	
	
	
}
