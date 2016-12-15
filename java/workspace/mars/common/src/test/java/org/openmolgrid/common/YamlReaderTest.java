package org.openmolgrid.common;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.JTree;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;

import junit.framework.Assert;

import org.junit.Test;
import org.openmolgrid.format.bigdft.CBigDFTInputWriter;
import org.openmolgrid.format.cml.CMLReader;
import org.openmolgrid.format.yaml.EnergiesYamlParse;
import org.openmolgrid.format.yaml.ForcesYamlParse;
import org.openmolgrid.format.yaml.GeometryYamlParse;
import org.openmolgrid.format.yaml.GroundStateOptimizationParser;
import org.openmolgrid.format.yaml.OrbitalsYamlParse;
import org.openmolgrid.format.yaml.YamlWindowPoper;
import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CElectromagCharges;
import org.openmolgrid.model.CEnergies;
import org.openmolgrid.model.CForces;
import org.openmolgrid.model.CGeometry;
import org.openmolgrid.model.COrbitals;
import org.openmolgrid.model.CStructure;
import org.yaml.snakeyaml.Yaml;

/**
 * @author Frederic Bonnet
 * 
 */
public class YamlReaderTest {
	private Yaml yaml = new Yaml();
	//input test files 
	private static final String INPUT_TIME_YAML_FILE = "time.yaml";
	private static final String INPUT_LOG_AGBULK_YAML_FILE = "log_AgBulk.yaml";
	private static final String INPUT_LOG_MULTIPT_YAML_FILE = "log_multiPt.yaml";
	private static final String INPUT_LOG_SINGLEPT_YAML_FILE = "log_singlePt.yaml";
	//output test files
	private static final String CML_OUTPUT_GEOMETRY_MULTIPT_YAML_FILE = "geometry_multiPt.cml";
	/**
	 * This test parses a bigdft result file and creates an instance of {@link CStructure}
	 * and others classes.
	 */
	@Test
	public void testYamlfile() {
		Reader reader = null;

		try {
			String 
			path = "src" + File.separator
					+ "test" + File.separator
					+ "resources" + File.separator
					+ "yaml" + File.separator
					+ INPUT_TIME_YAML_FILE;
//			path = "/home/unicore/server/unicorex/data/FILESPACE/d750f165-f87a-49e1-a11b-9c682321b0c4/log.yaml";
//			path = "/home/unicore/server/unicorex/data/FILESPACE/97b9675b-18ab-40ca-af4a-44465ffd55c9/log.yaml";
			
			File testFile = new File(path);
	
			reader = new FileReader(testFile);
			Assert.assertTrue(testFile.exists());
			Iterable<Object> allYaml = yaml.loadAll(new FileInputStream(testFile));

			Map<?, ?> yamlMap = null;
			Object lastYaml = null;
			for (Object currentYaml : allYaml) {
				//System.out.println(currentYaml);
				lastYaml = currentYaml;
			}

			yamlMap = (Map<?, ?>) lastYaml;
			//System.out.println(yaml.dump(yamlMap));

			Assert.assertNotNull(yamlMap);
			
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

	}
	@Test
	public void testYamlExist() {
		Reader reader = null;
	
		try {
			File testFile = new File("src" + File.separator + "test" + File.separator + "resources" + File.separator
					+ "yaml" + File.separator + INPUT_TIME_YAML_FILE);
	
			reader = new FileReader(testFile);
			Assert.assertTrue(testFile.exists());
				
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	 }
	/**
	 * This test parses a bigdft result file and creates an instance of {@link CStructure}
	 */
	@Test
	public void testGeometryAtomParsing() throws IOException {
		Reader reader = null;
		int id = 0;
		try {
			System.out.println("/**Start of the testGeometryAtomParsing()*/\n");
			CStructure structure = new CStructure();
			String //logpath = "/home/unicore/server/unicorex/data/FILESPACE/97b9675b-18ab-40ca-af4a-44465ffd55c9/log.yaml";
			logpath = "src" + File.separator + "test" + File.separator + "resources" + File.separator
					+ "yaml" + File.separator + INPUT_LOG_MULTIPT_YAML_FILE;
			File testFile = new File(logpath);

			reader = new FileReader(testFile);
			Assert.assertTrue(testFile.exists());
			Iterable<Object> allYaml = yaml.loadAll(new FileInputStream(testFile));

			Map<?, ?> yamlMap = null;
			Object lastYaml = null;
			for (Object currentYaml : allYaml) {
				//System.out.println(currentYaml);
				lastYaml = currentYaml;
			}

			yamlMap = (Map<?, ?>) lastYaml;
			//System.out.println(yaml.dump(yamlMap));
			
			Assert.assertNotNull(yamlMap);

			ArrayList<?> atoms = (ArrayList<?>) yamlMap.get("Atomic positions within the cell (Atomic and Grid Units)");

			for (Object currentAtom : atoms) {
				Collection<?> currentAtomNameSet = ((LinkedHashMap<?, ?>) currentAtom).keySet();
				for (Object currentAtomNameString : currentAtomNameSet) {

					LinkedHashMap<?, ?> auGuCoordinates = (LinkedHashMap<?, ?>) ((LinkedHashMap<?, ?>) currentAtom)
							.get(currentAtomNameString);
					List<?> coords = (ArrayList<?>) auGuCoordinates.get("AU");

					double x = (Double) coords.get(0);
					double y = (Double) coords.get(1);
					double z = (Double) coords.get(2);
					CAtom newAtom = new CAtom((String) currentAtomNameString, x, y, z, id);
					structure.addAtom(newAtom);
					id++;
				}
			}

			Assert.assertEquals(5, structure.getAtomCount());								
			//Assert.assertEquals(4, structure.getAtomCount());								
			String cmlString = structure.writeCML2String();			
			System.out.println(cmlString);
			
			//Now getting the the boundary conditions of the system
			System.out.println("Now getting the the boundary conditions of the system");

			LinkedHashMap<?, ?> poisson_kernel = (LinkedHashMap<?, ?>) yamlMap.get("Poisson Kernel Creation");
			Collection<?> poissonKernelNameSet = poisson_kernel.keySet();
			String boundaryConditions = null;
			for (Object currentKernelNameString : poissonKernelNameSet ) {
				System.out.println(currentKernelNameString);
				if (poisson_kernel.get(currentKernelNameString) instanceof String ) {
					boundaryConditions = (String)poisson_kernel.get(currentKernelNameString);
				}
			}
			System.out.println(boundaryConditions);
			GeometryYamlParse poisson_kernel_test = new GeometryYamlParse(poisson_kernel);
			poisson_kernel_test.ParsePoissonKernel();
			Assert.assertEquals(poisson_kernel, poisson_kernel_test.getPoisson_kernel());
			Assert.assertEquals(boundaryConditions, poisson_kernel_test.getBoundaryConditions());
			
			//Now getting the shift in the atomic position
			System.out.println("Now getting the shift in the atomic position");
			List<?> rigidShiftList = null;
			rigidShiftList = (List<?>)yamlMap.get("Rigid Shift Applied (AU)");
			System.out.println(rigidShiftList+" entry 0: "+rigidShiftList.get(0));
			GeometryYamlParse rigidShiftList_test = new GeometryYamlParse(rigidShiftList);
			rigidShiftList_test.ParseRigidShift();
			Assert.assertEquals(rigidShiftList, rigidShiftList_test.getRigidShiftList());
			
			//now getting the Sizes of the simulation domain:
			System.out.println("now getting the Sizes of the simulation domain:");
			LinkedHashMap<?, ?> simulation_domain = (LinkedHashMap<?, ?>)yamlMap.get("Sizes of the simulation domain");
			System.out.println(simulation_domain);
			Collection<?> simulationDomainNameSet = simulation_domain.keySet();
			List<?> simulationDomain_AU_List = null;
			List<?> simulationDomain_Angstroem_List = null;
			List<?> simulationDomain_GU_List = null;
			for (Object currentSimulationDomainNameString : simulationDomainNameSet ) {
				System.out.println(currentSimulationDomainNameString);
				if (simulation_domain.get(currentSimulationDomainNameString) instanceof List<?> ) {
					if ( currentSimulationDomainNameString.equals("AU")) {
						simulationDomain_AU_List = (List<?>)simulation_domain.get(currentSimulationDomainNameString);
					} else if (currentSimulationDomainNameString.equals("Angstroem")) {
						simulationDomain_Angstroem_List = (List<?>)simulation_domain.get(currentSimulationDomainNameString);
					} else if (currentSimulationDomainNameString.equals("Grid Spacing Units")) {
						simulationDomain_GU_List = (List<?>)simulation_domain.get(currentSimulationDomainNameString);
					}
				}
			}
			System.out.println(simulationDomain_AU_List);
			System.out.println(simulationDomain_Angstroem_List);
			System.out.println(simulationDomain_GU_List);
			GeometryYamlParse simulation_domain_test = new GeometryYamlParse(simulation_domain);
			simulation_domain_test.ParseSimulationdomain();
			List<?> simulationDomain_Angstroem_List_test = simulation_domain_test.getSimulationDomain_Angstroem_List();
			Assert.assertEquals(simulationDomain_Angstroem_List, simulationDomain_Angstroem_List_test);
						
			//second way via calling GeometryYamlParse.java class
			
			GeometryYamlParse geometryYamlParse = new GeometryYamlParse(atoms);
			geometryYamlParse.ParseGeometry();
			CGeometry cGeometry = new CGeometry(
					geometryYamlParse.getXcoordArrayList(),
					geometryYamlParse.getYcoordArrayList(),
					geometryYamlParse.getZcoordArrayList(),
					geometryYamlParse.getAtomNamecoordArrayList() );
			cGeometry.doStructure();			
			String cmlString2 = cGeometry.getStructure().writeCML2String();
			System.out.println(cmlString2);
			//checking the two methods are identical and returning the same result
			Assert.assertEquals(cmlString, cmlString2);
	
			//writing the cmlString to file as a bufferedWriter
			String path = "src" + File.separator + "test" + File.separator +
					"resources" + File.separator + "yaml" + File.separator +
					CML_OUTPUT_GEOMETRY_MULTIPT_YAML_FILE;

			CBigDFTInputWriter cBigDFTInputWriter = new CBigDFTInputWriter();
			cBigDFTInputWriter.writeCMLString(path, cmlString2);

			//Now testing the reading back in of the cmlString.
			File testIOFile = new File(path);
			CStructure structure2 = new CStructure();
			CMLReader cmlreader = new CMLReader(testIOFile.getAbsolutePath(), null, structure2);
			cmlreader.provideStructures();
			Assert.assertNotNull(structure2);
			
			//printing the entire list
			System.out.println(structure2.getAtomCount()+" "+structure2.getAtoms());

			for(int i = 1 ; i <= structure2.getAtomCount() ; i++) {
				System.out.println(structure2.getAtomCount()+" i:"+i+" "+structure2.getAtom(i));
			}
			
			System.out.println();
			System.out.println("/**End of the testGeometryAtomParsing()*/\n");
			
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	/**
	 * test the parsing of the:
	 * 	charges
	 * 	dipole moment
	 */
	@Test
	public void testElectromagChargesParsing() {
		Reader reader = null;
		try {
			System.out.println("/**Start of the testElectromagChargesParsing()*/\n");
			File testFile = new File("src" + File.separator + "test"
					+ File.separator + "resources"
					+ File.separator + "yaml"
					+ File.separator + INPUT_LOG_AGBULK_YAML_FILE);
			
			reader = new FileReader(testFile);
			Assert.assertTrue(testFile.exists());
			Iterable<Object> allYaml = yaml.loadAll(new FileInputStream(testFile));

			Map<?, ?> yamlMap = null;
			Object lastYaml = null;
			for (Object currentYaml : allYaml) {
				//System.out.println(currentYaml);
				lastYaml = currentYaml;
			}
			yamlMap = (Map<?, ?>) lastYaml;
			Assert.assertNotNull(yamlMap);
			
			double eQTot = (Double) yamlMap.get("Total electronic charge");
			Assert.assertNotNull(eQTot);
			Map<?, ?> eDipole_AU = (Map<?, ?>)yamlMap.get("Electric Dipole Moment (AU)");
			Assert.assertNotNull(eDipole_AU);
			Map<?, ?> eDipole_Debye = (Map<?, ?>)yamlMap.get("Electric Dipole Moment (Debye)");
			Assert.assertNotNull(eDipole_Debye);

			CElectromagCharges electromagCharges = new CElectromagCharges(eQTot, eDipole_AU, eDipole_Debye);
			Assert.assertEquals("Checking the eQtot", eQTot, electromagCharges.getTotalElectronicQ());
			//asserting the AU calculated norm and the parsed norm are the same 
			electromagCharges.getPVector(eDipole_AU);
			Assert.assertEquals(electromagCharges.getPnorm(), electromagCharges.getPnormCalulated(), 0.001);
			//asserting the Debye calculated norm and the parsed norm are the same 
			electromagCharges.getPVector(eDipole_Debye);
			Assert.assertEquals(electromagCharges.getPnorm(), electromagCharges.getPnormCalulated(), 0.001);
			//testing the total eltronic charge is not 0
			Assert.assertEquals(electromagCharges.getTotalElectronicQ(),eQTot);
			
			electromagCharges.setPnormCalulated(0.0);
			electromagCharges.getPVector_AU();			
			electromagCharges.getPNorm_AU();
			System.out.println(electromagCharges.getPnorm());
			System.out.println(electromagCharges.getPnormCalulated());
			
			Assert.assertEquals(electromagCharges.getPnorm(), electromagCharges.getPnormCalulated(), 0.001);
			
			System.out.println();
			System.out.println("/**Start of the testElectromagChargesParsing()*/\n");
			
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		} finally {
			if (reader != null) {
				try {	
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		
	}
	/**
	 * test of the parsing: 
	 *	forces
	 */
	@Test
	public void testForcesParsing() throws NullPointerException {
		Reader reader = null;
		try {
			System.out.println("/**Start of the testForcesParsing()*/\n");
			File testFile = new File("src" + File.separator + "test"
					+ File.separator + "resources"
					+ File.separator + "yaml"
					+ File.separator + INPUT_LOG_AGBULK_YAML_FILE);
			
			reader = new FileReader(testFile);
			Assert.assertTrue(testFile.exists());
			Iterable<Object> allYaml = yaml.loadAll(new FileInputStream(testFile));

			Map<?, ?> yamlMap = null;
			Object lastYaml = null;
			for (Object currentYaml : allYaml) {
				//System.out.println(currentYaml);
				lastYaml = currentYaml;
			}
			yamlMap = (Map<?, ?>) lastYaml;
			Assert.assertNotNull(yamlMap);

			Object nlocForces = (Object)yamlMap.get("Non Local forces calculated");
			Assert.assertEquals(true, nlocForces);

			LinkedHashMap<?, ?> stressTensorsLinkedHashMap = (LinkedHashMap<?, ?>) yamlMap
						.get("Stress Tensor");
			//System.out.println(stressTensorsLinkedHashMap);
			try {
				Collection<?> stressTensorSet = ((LinkedHashMap<?, ?>) ((Map<?, ?>) yamlMap.get("Stress Tensor"))).keySet();

				ForcesYamlParse forcesYamlParse = new ForcesYamlParse(stressTensorsLinkedHashMap, stressTensorSet);
				forcesYamlParse.ParseForces();

				System.out.println("The number of stress tensor: "+forcesYamlParse.getNstressTensor());
				System.out.println("The number of maps: "+forcesYamlParse.getNmaps());
				System.out.println("The number of instances of integer: "+forcesYamlParse.getNints()+" number of symmetries: "+forcesYamlParse.getnSymmetries());

				Assert.assertNotNull(forcesYamlParse.getPressureArrayList()/**pressureArrayList*/);
			
				CForces forces = new CForces("",
						forcesYamlParse.getnSymmetries(),
						forcesYamlParse.getStressTensorMatrix(),
						forcesYamlParse.getPressureSet(),
						forcesYamlParse.getPressureArrayList());

				Assert.assertEquals(forces.getnSymmetries(), forcesYamlParse.getnSymmetries());
				
				Object[] stressTensorTSringArray = new String[stressTensorSet.size()];
				stressTensorTSringArray = stressTensorSet.toArray();
				int whichentry = java.util.Arrays.binarySearch(stressTensorTSringArray, "Total stress tensor matrix (Ha/Bohr^3)");
				Assert.assertEquals(whichentry, stressTensorSet.size()-2);
				Assert.assertEquals("Total stress tensor matrix (Ha/Bohr^3)", stressTensorTSringArray[stressTensorSet.size()-2]);

			} catch (Exception e) {
				System.out.println("Warning There is not an entry \"Stress Tensor\" in this Yaml file");
			}
				
			System.out.println();
			System.out.println("/**End of the testForcesParsing()*/\n");
			
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		} finally {
			if (reader != null) {
				try {	
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		
	}
	/**
	 * test to parse the Kohn-Sham eigenvalues
	 */
	@Test
	public void testOrbitalsParsing() {
		Reader reader = null;
		
		try {
			System.out.println("/**Start of the testOrbitalsParsing()*/\n");
			File testFile = new File("src" + File.separator + "test"
					+ File.separator + "resources"
					+ File.separator + "yaml"
					+ File.separator + INPUT_LOG_AGBULK_YAML_FILE);

			reader = new FileReader(testFile);
			Assert.assertTrue(testFile.exists());
			Iterable<Object> allYaml = yaml.loadAll(new FileInputStream(testFile));
			
			Map<?, ?> yamlMap = null;
			Object lastYaml = null;
			for (Object currentYaml : allYaml) {
				//System.out.println(currentYaml);
				lastYaml = currentYaml;
			}

			yamlMap = (Map<?, ?>) lastYaml;

			Assert.assertNotNull(yamlMap);

			Map<?, ?> occupational_number = (Map<?, ?>)yamlMap.get("Occupation Numbers");
			Assert.assertNotNull(occupational_number);
			System.out.println(occupational_number);
			Map<?, ?> input_Hamiltonian = (Map<?, ?>)yamlMap.get("Input Hamiltonian");
			System.out.println(input_Hamiltonian);

			ArrayList<?> orbitals = (ArrayList<?>)input_Hamiltonian.get("Orbitals");
			//options.setIndent(10);
			//options.setPrettyFlow(false);
			System.out.println(yaml.dump(orbitals));
			//System.out.println(yaml.dump(orbitals, out)(orbitals));
			
			OrbitalsYamlParse orbitalsYamlParse = new OrbitalsYamlParse(orbitals);
			orbitalsYamlParse.ParseOrbitals();
			
			COrbitals corbitals = new COrbitals(
					orbitalsYamlParse.getArrayList_e(),
					orbitalsYamlParse.getArrayList_f(),
					orbitalsYamlParse.getArrayList_k());
			
			corbitals.toString();

			CBigDFTInputWriter cBigDFTInputWriter = new CBigDFTInputWriter(corbitals);
			String cmlWriter = CBigDFTInputWriter.orbitalsAsString(corbitals);
			System.out.println(cmlWriter);			
			/**writing the new kpt_file with the new kpt points*/
			cBigDFTInputWriter.kptInputWritter("");
			
			System.out.println();
			System.out.println("/**End of the testOrbitalsParsing()*/\n");
			
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	/**
	 * Testing the parsing of the ground state optimization
	 */
	@Test
	public void testGroundStateOptimization() throws NullPointerException, ArrayIndexOutOfBoundsException {
		Reader reader = null;
		try {
			System.out.println("/**Start of the testGroundStateOptimization()*/\n");
			File testFile = new File("src" + File.separator + "test"
					+ File.separator + "resources"
					+ File.separator + "yaml"
					+ File.separator + INPUT_LOG_AGBULK_YAML_FILE);
			
			reader = new FileReader(testFile);
			Assert.assertTrue(testFile.exists());
			Iterable<Object> allYaml = yaml.loadAll(new FileInputStream(testFile));

			Map<?, ?> yamlMap = null;
			Object lastYaml = null;
			for (Object currentYaml : allYaml) {
				//System.out.println(currentYaml);
				lastYaml = currentYaml;
			}
			yamlMap = (Map<?, ?>) lastYaml;
			Assert.assertNotNull(yamlMap);

			ArrayList<?> gS_Optimi= (ArrayList<?>)yamlMap.get("Ground State Optimization");
			
			int niter = gS_Optimi.size();
			Map<?, ?> hamil_optimi = (Map<?, ?>)gS_Optimi.get(niter-1);
			System.out.println(hamil_optimi+"\n");
			
			EnergiesYamlParse fermiEnergiesYamlParse = new EnergiesYamlParse(hamil_optimi);
			fermiEnergiesYamlParse.ParseFermiEnergy();
			
			CEnergies cfermienergies = new CEnergies(fermiEnergiesYamlParse.getFermiEnergy());
			
			EnergiesYamlParse energiesYamlParse = new EnergiesYamlParse(hamil_optimi);
			energiesYamlParse.ParseEnergies();
			
			CEnergies cEnergies = new CEnergies(
					energiesYamlParse.getEnergyArrayListName(),
					energiesYamlParse.getEnergyArrayListValues());
			
			try {
				Collection<?> hamil_optimiSet = ((Map<?, ?>) hamil_optimi).keySet();

				int i = 0;
				Object[] hamil_optimiStringArray = new String[hamil_optimiSet.size()];
				hamil_optimiStringArray = hamil_optimiSet.toArray();
				i = 0;
				int whichFermi = -99;
				int whichHamil_optimi = 99;
				for ( Object current : hamil_optimiStringArray ) {
					if ( current.equals("Fermi Energy")) whichFermi = i;
					if ( current.equals("Hamiltonian Optimization")) whichHamil_optimi = i;
					i++;
				}
				/**The Fermi energy parsing*/
				try {
					double fermiEnergy = (Double)hamil_optimi.get(hamil_optimiStringArray[whichFermi]);
					Assert.assertEquals(fermiEnergy, cfermienergies.getFermiEnergy());
				} catch (Exception e) {
					System.out.println("The array is out of bound: "+whichFermi);
					System.out.println("The file does not contain the entry --> Fermi Energy:"); 
					System.out.println("The exception printout: "+e.toString());
					
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
					System.out.println("This is an instance of (Double)");
					System.out.println(groundStateOptimizationParser.getEnDouble());
//					System.out.println(energyDouble);
				}else if ( !groundStateOptimizationParser.getEnMap().equals(null) ) {
					System.out.println("This is an instance of (Map<?,?>)");
					System.out.println(groundStateOptimizationParser.getEnMap());
//					System.out.println(energyMap);
				}else if ( groundStateOptimizationParser.getEnString() != null ) {
					System.out.println("This is an instance of (String)");
					System.out.println(groundStateOptimizationParser.getEnString());
//					System.out.println(energyString);
				}else if ( groundStateOptimizationParser.getEnInt() != -99 ) {
					System.out.println("This is an instance of (Integer)");
					System.out.println(groundStateOptimizationParser.getEnInt());
//	                System.out.println(energyInt);
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
					System.out.println("for n: "+n+" "+currentEnergy);
					EnergyArrayListName.add((String)currentEnergy);
					Assert.assertEquals(energiesYamlParse.getEnergyArrayListName().get(n), EnergyArrayListName.get(n));
					EnergyArrayListValues.add((Double)groundStateOptimizationParser.getEnMap().get(currentEnergy));
					Assert.assertEquals(energiesYamlParse.getEnergyArrayListValues().get(n), EnergyArrayListValues.get(n));
					System.out.println(EnergyArrayListName.get(n)+"-->"+EnergyArrayListValues.get(n));
					n++;
				}
				
			} catch (Exception e) {
				System.out.println("Handler for the catch and the exception in testGroundStateOptimization");
				System.out.println("The exception printout: "+e.toString());
			}
			
			System.out.println();
			System.out.println("/**End of the testGroundStateOptimization()*/\n");
			
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		} finally {
			if (reader != null) {
				try {	
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	/**
	 * Test parsing using loops.
	 */
	@Test
	public void testYamlReader() {
		Reader reader = null;
		try {
			File testFile = new File("src" + File.separator + "test"
					+ File.separator + "resources"
					+ File.separator + "yaml"
					+ File.separator + INPUT_LOG_AGBULK_YAML_FILE);
			
			reader = new FileReader(testFile);
			Assert.assertTrue(testFile.exists());
			Iterable<Object> allYaml = yaml.loadAll(new FileInputStream(testFile));

			Map<?, ?> yamlMap = null;
			Object lastYaml = null;
			for (Object currentYaml : allYaml) {
				//System.out.println(currentYaml);
				lastYaml = currentYaml;
			}
			yamlMap = (Map<?, ?>) lastYaml;
			Assert.assertNotNull(yamlMap);
			
			Set<?> sets = yamlMap.keySet();
	        Object[] arraySets = new String[sets.size()];
	        arraySets = sets.toArray();
			Assert.assertNotNull(arraySets);
			
			boolean actual = false;
			for (int iset = 0 ; iset < sets.size() ; iset++ ) {
				if ( arraySets[iset].equals("Input Hamiltonian")) {
					actual = true;
				}
				//System.out.println("iset= "+iset+" arraySets["+iset+"]:"+arraySets[iset]);
	        }
			System.out.println("\n");
			Assert.assertEquals(true, actual);
			
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		} finally {
			if (reader != null) {
				try {	
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	/**
	 * Testing the GUI interface.
	 */
	//@Test
	public void testMakeGUIPoper() {
		Reader reader = null;
	    //the recursively generated tree
		JTree timeYamlRecursivelyTree;
	    DefaultTreeModel /**TreeModel*/ m_modelRecursivelyTree;  
	    DefaultMutableTreeNode treeNodeTimeYamlParse_RecursivelyTree_root;
		//the full tree
		JTree timeYamlFullTree;
	    DefaultTreeModel m_modelFullTree;  
	    DefaultMutableTreeNode treeNodeTimeYamlParse_FullTree_root;

		try {
			File testFile = new File("src" + File.separator + "test"
					+ File.separator + "resources"
					+ File.separator + "yaml"
					+ File.separator + INPUT_LOG_AGBULK_YAML_FILE);
		
			reader = new FileReader(testFile);
			Assert.assertTrue(testFile.exists());
			Iterable<Object> allYaml = yaml.loadAll(new FileInputStream(testFile));

			Map<?, ?> yamlMap = null;
			Object lastYaml = null;
			for (Object currentYaml : allYaml) {
				//System.out.println(currentYaml);
				lastYaml = currentYaml;
			}
			yamlMap = (Map<?, ?>) lastYaml;
			Assert.assertNotNull(yamlMap);

	        //creating the tree
			treeNodeTimeYamlParse_FullTree_root = new DefaultMutableTreeNode("Time Yaml Parse");
	        m_modelFullTree = new DefaultTreeModel(treeNodeTimeYamlParse_FullTree_root);
			timeYamlFullTree = new JTree(m_modelFullTree);
	        timeYamlFullTree.setEditable(true);  
	        timeYamlFullTree.setSelectionRow(0);
	        timeYamlFullTree.setName("timeYamlFullTree");

		    //creating the tree
			treeNodeTimeYamlParse_RecursivelyTree_root = new DefaultMutableTreeNode("Time Yaml Parse");
		    m_modelRecursivelyTree = new DefaultTreeModel(treeNodeTimeYamlParse_RecursivelyTree_root);
		    timeYamlRecursivelyTree = new JTree(m_modelRecursivelyTree);
		    timeYamlRecursivelyTree.setEditable(true);  
		    timeYamlRecursivelyTree.setSelectionRow(0);
		    timeYamlRecursivelyTree.setName("timeYamlRecursivelyTree");
			
			YamlWindowPoper yamlwindowpoper = new YamlWindowPoper(
					timeYamlFullTree,
					treeNodeTimeYamlParse_RecursivelyTree_root,
					timeYamlRecursivelyTree,
					m_modelRecursivelyTree);
			//yamlwindowpoper.WindowPoper();

			
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		} finally {
			if (reader != null) {
				try {	
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
}
