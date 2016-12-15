/**
 * 
 */
package org.openmolgrid.format.yaml;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;

import org.apache.log4j.Logger;

/**
 * Class to parse the Geometry from the yaml file generated from BigDFT Output.
 * @author Frederic Bonnet
 *
 */
public class GeometryYamlParse {
	/** Logger */
	private static Logger log = Logger.getLogger(GeometryYamlParse.class.getName());

	//local variables
	private String boundaryConditions = null;
	private ArrayList<Double> XcoordArrayList = null;
	private ArrayList<Double> YcoordArrayList = null;
	private ArrayList<Double> ZcoordArrayList = null;
	private ArrayList<String> AtomNamecoordArrayList = null;

	private List<?> simulationDomain_AU_List = null;
	private List<?> simulationDomain_Angstroem_List = null;
	private List<?> simulationDomain_GU_List = null;
	private LinkedHashMap<?, ?> poisson_kernel = null;
	private LinkedHashMap<?, ?> simulation_domain = null;
	private Collection<?> simulationDomainNameSet;
	private Collection<?> poissonKernelNameSet;
	//global variables
	private ArrayList<?> atoms;
	private List<?> rigidShiftList;
	private LinkedHashMap<?, ?> linkedHashMap;
	/**
	 * simple constructor
	 */
	public GeometryYamlParse() {
		
	}
	/**
	 * constructor for the GeometryYamlParse class
	 * @param atoms:ArrayList<?>
	 */
	public GeometryYamlParse(ArrayList<?> atoms) {
		this.atoms = atoms;
	}
	/**
	 * constructor for the GeometryYamlParse class
	 * @param rigidShiftList:List<?>
	 */
	public GeometryYamlParse(List<?> rigidShiftList) {
		this.rigidShiftList = rigidShiftList;
	}
	/**
	 * constructor for the GeometryYamlParse class
	 * @param linkedHashMap:LinkedHashMap<?, ?>
	 */
	public GeometryYamlParse(LinkedHashMap<?, ?> linkedHashMap) {
		this.linkedHashMap = linkedHashMap;
	}
	/**
	 * method to parse the geometry from the BigDFT output file.
	 */
	public void ParseGeometry() {

		System.out.println(atoms);
		XcoordArrayList = new ArrayList<Double>();
		YcoordArrayList = new ArrayList<Double>();
		ZcoordArrayList = new ArrayList<Double>();
		AtomNamecoordArrayList = new ArrayList<String>();
		for (Object currentAtom : atoms) {
			Collection<?> currentAtomNameSet = ((LinkedHashMap<?, ?>) currentAtom).keySet();
			for (Object currentAtomNameString : currentAtomNameSet) {

				LinkedHashMap<?, ?> auGuCoordinates = (LinkedHashMap<?, ?>) ((LinkedHashMap<?, ?>) currentAtom)
						.get(currentAtomNameString);
				List<?> coords = (ArrayList<?>) auGuCoordinates.get("AU");

				double x = (Double) coords.get(0);
				double y = (Double) coords.get(1);
				double z = (Double) coords.get(2);
				System.out.println(x+" "+y+" "+" "+z+" "+currentAtomNameString);
				XcoordArrayList.add(x);
				YcoordArrayList.add(y);
				ZcoordArrayList.add(z);
				AtomNamecoordArrayList.add((String)currentAtomNameString);

			}
		}
		setAtomNamecoordArrayList(AtomNamecoordArrayList);
		setXcoordArrayList(XcoordArrayList);
		setYcoordArrayList(YcoordArrayList);
		setZcoordArrayList(ZcoordArrayList);
	}
	/**
	 * method to parse the Poisson Kernel
	 */
	public void ParsePoissonKernel() {
		poisson_kernel = (LinkedHashMap<?, ?>)linkedHashMap;
		setPoisson_kernel(poisson_kernel);
		poissonKernelNameSet = poisson_kernel.keySet();
		setPoissonKernelNameSet(poissonKernelNameSet);
		for (Object currentKernelNameString : poissonKernelNameSet ) {
			log.info(currentKernelNameString);
			if (poisson_kernel.get(currentKernelNameString) instanceof String ) {
				boundaryConditions = (String)poisson_kernel.get(currentKernelNameString);
				setBoundaryConditions(boundaryConditions);
			}
		}
	}
	/**
	 * method to parse the rigid Shift
	 */
	public void ParseRigidShift() {
		setRigidShiftList(rigidShiftList);
	}
	/**
	 * method to parse the simulation domain
	 */
	public void ParseSimulationdomain() {
		simulation_domain = (LinkedHashMap<?, ?>)linkedHashMap;
		setSimulation_domain(simulation_domain);
		log.info(simulation_domain);
		simulationDomainNameSet = simulation_domain.keySet();
		setSimulationDomainNameSet(simulationDomainNameSet);
		for (Object currentSimulationDomainNameString : simulationDomainNameSet ) {
			log.info(currentSimulationDomainNameString);
			if (simulation_domain.get(currentSimulationDomainNameString) instanceof List<?> ) {
				if ( currentSimulationDomainNameString.equals("AU")) {
					simulationDomain_AU_List = (List<?>)simulation_domain.get(currentSimulationDomainNameString);
					setSimulationDomain_AU_List(simulationDomain_AU_List);
				} else if (currentSimulationDomainNameString.equals("Angstroem")) {
					simulationDomain_Angstroem_List = (List<?>)simulation_domain.get(currentSimulationDomainNameString);
					setSimulationDomain_Angstroem_List(simulationDomain_Angstroem_List);
				} else if (currentSimulationDomainNameString.equals("Grid Spacing Units")) {
					simulationDomain_GU_List = (List<?>)simulation_domain.get(currentSimulationDomainNameString);
					setSimulationDomain_GU_List(simulationDomain_GU_List);
				}
			}
		}
	}
	/**
	 * The setters and getters for the:
	 * //local variables
	 * 	private String boundaryConditions = null;
	 * 	private ArrayList<Double> XcoordArrayList = null;
	 * 	private ArrayList<Double> YcoordArrayList = null;
	 * 	private ArrayList<Double> ZcoordArrayList = null;
	 * 	private ArrayList<String> AtomNamecoordArrayList = null;
	 * 	private List<?> simulationDomain_AU_List = null;
	 * 	private List<?> simulationDomain_Angstroem_List = null;
	 * 	private List<?> simulationDomain_GU_List = null;
	 * 	private LinkedHashMap<?, ?> poisson_kernel = null;
	 * 	private LinkedHashMap<?, ?> simulation_domain = null;
	 * 	private Collection<?> simulationDomainNameSet;
	 * 	private Collection<?> poissonKernelNameSet;
	 * 	//global variables
	 * 	private ArrayList<?> atoms;
	 * 	private List<?> rigidShiftList;
	 * 	private LinkedHashMap<?, ?> linkedHashMap;
	 * 
	 */
	/**
	 * @return the atoms
	 */
	public ArrayList<?> getAtoms() {
		return atoms;
	}
	/**
	 * @param atoms the atoms to set
	 */
	public void setAtoms(ArrayList<?> atoms) {
		this.atoms = atoms;
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
	/**
	 * @return the boundaryConditions
	 */
	public String getBoundaryConditions() {
		return boundaryConditions;
	}
	/**
	 * @param boundaryConditions the boundaryConditions to set
	 */
	public void setBoundaryConditions(String boundaryConditions) {
		this.boundaryConditions = boundaryConditions;
	}
	/**
	 * @return the simulationDomain_AU_List
	 */
	public List<?> getSimulationDomain_AU_List() {
		return simulationDomain_AU_List;
	}
	/**
	 * @param simulationDomain_AU_List the simulationDomain_AU_List to set
	 */
	public void setSimulationDomain_AU_List(List<?> simulationDomain_AU_List) {
		this.simulationDomain_AU_List = simulationDomain_AU_List;
	}
	/**
	 * @return the simulationDomain_Angstroem_List
	 */
	public List<?> getSimulationDomain_Angstroem_List() {
		return simulationDomain_Angstroem_List;
	}
	/**
	 * @param simulationDomain_Angstroem_List the simulationDomain_Angstroem_List to set
	 */
	public void setSimulationDomain_Angstroem_List(
			List<?> simulationDomain_Angstroem_List) {
		this.simulationDomain_Angstroem_List = simulationDomain_Angstroem_List;
	}
	/**
	 * @return the simulationDomain_GU_List
	 */
	public List<?> getSimulationDomain_GU_List() {
		return simulationDomain_GU_List;
	}
	/**
	 * @param simulationDomain_GU_List the simulationDomain_GU_List to set
	 */
	public void setSimulationDomain_GU_List(List<?> simulationDomain_GU_List) {
		this.simulationDomain_GU_List = simulationDomain_GU_List;
	}
	/**
	 * @return the poisson_kernel
	 */
	public LinkedHashMap<?, ?> getPoisson_kernel() {
		return poisson_kernel;
	}
	/**
	 * @param poisson_kernel the poisson_kernel to set
	 */
	public void setPoisson_kernel(LinkedHashMap<?, ?> poisson_kernel) {
		this.poisson_kernel = poisson_kernel;
	}
	/**
	 * @return the simulation_domain
	 */
	public LinkedHashMap<?, ?> getSimulation_domain() {
		return simulation_domain;
	}
	/**
	 * @param simulation_domain the simulation_domain to set
	 */
	public void setSimulation_domain(LinkedHashMap<?, ?> simulation_domain) {
		this.simulation_domain = simulation_domain;
	}
	/**
	 * @return the simulationDomainNameSet
	 */
	public Collection<?> getSimulationDomainNameSet() {
		return simulationDomainNameSet;
	}
	/**
	 * @param simulationDomainNameSet the simulationDomainNameSet to set
	 */
	public void setSimulationDomainNameSet(Collection<?> simulationDomainNameSet) {
		this.simulationDomainNameSet = simulationDomainNameSet;
	}
	/**
	 * @return the poissonKernelNameSet
	 */
	public Collection<?> getPoissonKernelNameSet() {
		return poissonKernelNameSet;
	}
	/**
	 * @param poissonKernelNameSet the poissonKernelNameSet to set
	 */
	public void setPoissonKernelNameSet(Collection<?> poissonKernelNameSet) {
		this.poissonKernelNameSet = poissonKernelNameSet;
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
	 * @return the linkedHashMap
	 */
	public LinkedHashMap<?, ?> getLinkedHashMap() {
		return linkedHashMap;
	}
	/**
	 * @param linkedHashMap the linkedHashMap to set
	 */
	public void setLinkedHashMap(LinkedHashMap<?, ?> linkedHashMap) {
		this.linkedHashMap = linkedHashMap;
	}
	
}
