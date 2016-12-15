/**
 * 
 */
package org.openmolgrid.model;

/**
 * @author Frederic Bonnet
 * 
 * Handler for the Kohn-Sham Values from
 * a BigDFT calculation.
 * 
 * It includes handler for:
 *  KkonShamEigenvalues()
 * 	OccupationNumbers()
 * 	KptBZCoords()
 */
public interface OrbitalsValues {
	public void kohnShamEigenvalues();
	public void occupationNumbers();
	public void kptBZCoords();
}
