package org.openmolgrid.util;

import org.openmolgrid.model.CStructure;

/**
 * The center of mass or mass center is the mean location of all the mass in a system. In the case of a rigid body, the
 * position of the center of mass is fixed in relation to the body. In the case of a loose distribution of masses in
 * free space, such as shot from a shotgun or the planets of the solar system, the position of the center of mass is a
 * point in space among them that may not correspond to the position of any individual mass. The use of the mass center
 * often allows the use of simplified equations of motion, and it is a convenient reference point for many other
 * calculations in physics, such as angular momentum or the moment of inertia. In many applications, such as orbital
 * mechanics, objects can be replaced by point masses located at their mass centers for the purposes of analysis.
 * (Source: http://en.wikipedia.org/wiki/Center_of_mass)
 * 
 * @author Stefan Bozic
 * 
 */
public class CenterOfMass {

	/** The structure to which this center of mass belongs to */
	private CStructure structure;

	/**
	 * the xyz coordinates of the center of mass
	 */
	private double x = 0.0, y = 0.0, z = 0.0;

	/**
	 * Creates a new instance of a Center of Mass
	 * 
	 * @param structure the referending structure
	 * @param x the x coordinate
	 * @param y the y coordinate
	 * @param z the z coordinate
	 */
	public CenterOfMass(CStructure structure, double x, double y, double z) {
		this.structure = structure;
		this.x = x;
		this.y = y;
		this.z = z;
	}

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
	 * @return the x
	 */
	public double getX() {
		return x;
	}

	/**
	 * @param x the x to set
	 */
	public void setX(double x) {
		this.x = x;
	}

	/**
	 * @return the y
	 */
	public double getY() {
		return y;
	}

	/**
	 * @param y the y to set
	 */
	public void setY(double y) {
		this.y = y;
	}

	/**
	 * @return the z
	 */
	public double getZ() {
		return z;
	}

	/**
	 * @param z the z to set
	 */
	public void setZ(double z) {
		this.z = z;
	}

	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "CenterOfMass [structure=" + structure + ", x=" + x + ", y=" + y + ", z=" + z + "]";
	}

}
