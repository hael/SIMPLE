/**
 * 
 */
package org.openmolgrid.model;

import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

/**
 * @author Frederic Bonnet
 *
 */
public class CGeometry implements GeometryValues {

	/** Logger */
	private static Logger log = Logger.getLogger(CGeometry.class.getName());

	//local variables
	private CStructure structure;
	private CStructure shiftedBackStructure;

	private static final int X_INDEX = 0;
	private static final int Y_INDEX = 1;
	private static final int Z_INDEX = 2;
	//global variables
	private ArrayList<Double> XcoordArrayList;
	private ArrayList<Double> YcoordArrayList;
	private ArrayList<Double> ZcoordArrayList;
	private ArrayList<String> AtomNamecoordArrayList;
	private List<?> rigidShiftList;
	/**
	 * simple constructor for CGeometry
	 */
	public CGeometry() {
	}
	/**
	 * constructor for CGeometry
	 *  
	 * @param XcoordArrayList
	 * @param YcoordArrayList
	 * @param ZcoordArrayList
	 * @param AtomNamecoordArrayList
	 */
	public CGeometry(ArrayList<Double> XcoordArrayList, ArrayList<Double> YcoordArrayList, ArrayList<Double> ZcoordArrayList, ArrayList<String> AtomNamecoordArrayList) {
		this.XcoordArrayList = XcoordArrayList;
		this.YcoordArrayList = YcoordArrayList;
		this.ZcoordArrayList = ZcoordArrayList;
		this.AtomNamecoordArrayList = AtomNamecoordArrayList;
	}
	/**
	 * constructor for CGeometry
	 * 
	 * @param XcoordArrayList
	 * @param YcoordArrayList
	 * @param ZcoordArrayList
	 * @param AtomNamecoordArrayList
	 * @param rigidShiftList
	 */
	public CGeometry(ArrayList<Double> XcoordArrayList, ArrayList<Double> YcoordArrayList, ArrayList<Double> ZcoordArrayList, ArrayList<String> AtomNamecoordArrayList, List<?> rigidShiftList) {
		this.XcoordArrayList = XcoordArrayList;
		this.YcoordArrayList = YcoordArrayList;
		this.ZcoordArrayList = ZcoordArrayList;
		this.AtomNamecoordArrayList = AtomNamecoordArrayList;
		this.rigidShiftList = rigidShiftList;
	}
	/**
	 * construct the structure for the given geometry.
	 */
	public void doStructure() {
		int id = 0;
		structure = new CStructure();

		int i = 0;
		double x,y,z;
		for (Object currentAtomNameString : AtomNamecoordArrayList ) {
			x = XcoordArrayList.get(i);
			y = YcoordArrayList.get(i);
			z = ZcoordArrayList.get(i);

			CAtom newAtom = new CAtom((String) currentAtomNameString, x, y, z, id);
			structure.addAtom(newAtom);
			id++;
			i++;
		}
	}
	/**
	 * construct the structure for the given geometry shifted back
	 * from the coordinates back from the Rigid Shift Applied (AU).
	 */
	public void doShiftBackStructure() {
		int id = 0;
		shiftedBackStructure = new CStructure();

		log.info((Double)rigidShiftList.get(X_INDEX)+" "+(Double)rigidShiftList.get(Y_INDEX)+" "+(Double)rigidShiftList.get(Z_INDEX));

		int i = 0;
		double x,y,z;
		for (Object currentAtomNameString : AtomNamecoordArrayList ) {
			x = XcoordArrayList.get(i);
			y = YcoordArrayList.get(i);
			z = ZcoordArrayList.get(i);

			//shifting the coordinates back from the Rigid Shift Applied (AU)
			x = x - (Double)rigidShiftList.get(X_INDEX);
			y = y - (Double)rigidShiftList.get(Y_INDEX);
			z = z - (Double)rigidShiftList.get(Z_INDEX);
			
			log.info(x+" "+y+" "+z+" "+(Double)rigidShiftList.get(X_INDEX)+" "+(Double)rigidShiftList.get(Y_INDEX)+" "+(Double)rigidShiftList.get(Z_INDEX));
			
			CAtom newAtom = new CAtom((String) currentAtomNameString, x, y, z, id);
			shiftedBackStructure.addAtom(newAtom);
			id++;
			i++;
		}
	}
	/**
	 * setters and getters for CGeometry: 
	 * 	//local variables
	 * 	private CStructure structure;
	 * 	//global variables
	 * 	private ArrayList<Double> XcoordArrayList;
	 * 	private ArrayList<Double> YcoordArrayList;
	 * 	private ArrayList<Double> ZcoordArrayList;
	 * 	private ArrayList<String> AtomNamecoordArrayList;
	 * 
	 */
	/**
	 * @return the structure
	 */
	public CStructure getStructure() {
		return structure;
	}
	/**
	 * @param structure the structure to set
	 */
	public void setStructure(CStructure structure) {
		this.structure = structure;
	}
	/**
	 * @return the shiftedBackStructure
	 */
	public CStructure getShiftedBackStructure() {
		return shiftedBackStructure;
	}
	/**
	 * @param shiftedBackStructure the shiftedBackStructure to set
	 */
	public void setShiftedBackStructure(CStructure shiftedBackStructure) {
		this.shiftedBackStructure = shiftedBackStructure;
	}
	/**
	 * @return the rigidShiftList
	 */
	public List<?> getRigidShiftList() {
		return rigidShiftList;
	}
	/**
	 * @param rigidShiftList the rigidShiftList to set
	 */
	public void setRigidShiftList(List<?> rigidShiftList) {
		this.rigidShiftList = rigidShiftList;
	}
	/**
	 * @return the xIndex
	 */
	public static int getxIndex() {
		return X_INDEX;
	}
	/**
	 * @return the yIndex
	 */
	public static int getyIndex() {
		return Y_INDEX;
	}
	/**
	 * @return the zIndex
	 */
	public static int getzIndex() {
		return Z_INDEX;
	}
	/**
	 * @return the xcoordArrayList
	 */
	public ArrayList<Double> getXcoordArrayList() {
		return XcoordArrayList;
	}
	/**
	 * @param xcoordArrayList the xcoordArrayList to set
	 */
	public void setXcoordArrayList(ArrayList<Double> xcoordArrayList) {
		XcoordArrayList = xcoordArrayList;
	}
	/**
	 * @return the ycoordArrayList
	 */
	public ArrayList<Double> getYcoordArrayList() {
		return YcoordArrayList;
	}
	/**
	 * @param ycoordArrayList the ycoordArrayList to set
	 */
	public void setYcoordArrayList(ArrayList<Double> ycoordArrayList) {
		YcoordArrayList = ycoordArrayList;
	}
	/**
	 * @return the zcoordArrayList
	 */
	public ArrayList<Double> getZcoordArrayList() {
		return ZcoordArrayList;
	}
	/**
	 * @param zcoordArrayList the zcoordArrayList to set
	 */
	public void setZcoordArrayList(ArrayList<Double> zcoordArrayList) {
		ZcoordArrayList = zcoordArrayList;
	}
	/**
	 * @return the atomNamecoordArrayList
	 */
	public ArrayList<String> getAtomNamecoordArrayList() {
		return AtomNamecoordArrayList;
	}
	/**
	 * @param atomNamecoordArrayList the atomNamecoordArrayList to set
	 */
	public void setAtomNamecoordArrayList(ArrayList<String> atomNamecoordArrayList) {
		AtomNamecoordArrayList = atomNamecoordArrayList;
	}
}
