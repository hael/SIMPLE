/**
 * 
 */
package org.openmolgrid.format.bigdft;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openmolgrid.format.yaml.GeometryYamlParse;

/**
 * @author Frederic Bonnet
 *
 */
public interface BigDFTInputWritter {

	public void kptInputWritter(String pth) throws IOException;
	public void dftInputWritter();
	public void xyzInputWritter(String pth) throws IOException;
	public void xyzInputWritter(String path,
			ArrayList<Double> XcoordArrayList,ArrayList<Double> YcoordArrayList, ArrayList<Double> ZcoordArrayList,
			ArrayList<String> AtomNamecoordArrayList, List<?> rigidShiftList,
			GeometryYamlParse poisson_kernelGeoYamlParse,
			GeometryYamlParse simulation_domainGeoYamlParse) throws IOException;
}
