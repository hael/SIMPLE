/**
 * 
 */
package org.openmolgrid.format.yaml;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;
import org.yaml.snakeyaml.Yaml;
/**
 * @author Frederic Bonnet
 *
 */
public class TimeYamlParse {

	//global variables
	private String path = null;
	private String filename = null;
	//local variables
	/** Logger */
	private static Logger log = Logger.getLogger(TimeYamlParse.class.getName());

	/** The Yaml parsing variables */
	private Yaml yaml = new Yaml();
	private TimeYamlVariables timeyamlparse = new TimeYamlVariables(yaml);
	/**
	 * Simple constructor
	 */
	public TimeYamlParse() {
	}
	public TimeYamlParse(String path, String filename ) {
		this.path = path;
		this.filename = filename;
	}
	/**
	 * method to set the variables from the parsing
	 * @throws Exception
	 */
	public void setTimeYamlVariablesFromStream() throws Exception {
		try {
			InputStream input = new FileInputStream(new File(path));
			//Object dataTimeYaml = yaml.load(input);

			Iterable<Object> allYaml = yaml.loadAll(input);

			Map<?, ?> yamlMap = null;
			Object lastYaml = null;
			for (Object currentYaml : allYaml) {
				//System.out.println(currentYaml);
				lastYaml = currentYaml;
			}

			//yamlMap = (Map<?, ?>) lastYaml;

			
			
			//this is the full data set loaded into a generic Map<?,?>.
			Map<?,?> mapdataTimeYaml = (Map<?,?>)lastYaml;

            //now loading INIT section

			Map<?,?> mapINIT = (Map<?,?>)mapdataTimeYaml.get("INIT");
			timeyamlparse.setMapINIT(mapINIT);
            System.out.println(mapdataTimeYaml.get("INIT"));            
            //now loading INIT section into their respective variables
            {/**start of the loading of INIT Classes:*/
            	Map<?,?> mapINITClasses = (Map<?,?>)mapINIT.get("Classes");
            	timeyamlparse.setMapINITClasses(mapINITClasses);
            	{/**start of the loading of INIT Classes: subclasses*/
            		List<?> listINITClassesCommunications = (List<?>)mapINITClasses.get("Communications");
            		timeyamlparse.setListINITClassesCommunications(listINITClassesCommunications);
            		List<?> listINITClassesConvolutions = (List<?>)mapINITClasses.get("Convolutions");
            		timeyamlparse.setListINITClassesConvolutions(listINITClassesConvolutions);
            		List<?> listINITClassesLinearAlgebra = (List<?>)mapINITClasses.get("Linear Algebra");
            		timeyamlparse.setListINITClassesLinearAlgebra(listINITClassesLinearAlgebra);
            		List<?> listINITClassesOther = (List<?>)mapINITClasses.get("Other");
            		timeyamlparse.setListINITClassesOther(listINITClassesOther);
            		List<?> listINITClassesPotential = (List<?>)mapINITClasses.get("Potential");
            		timeyamlparse.setListINITClassesPotential(listINITClassesPotential);
            		List<?> listINITClassesInitialization = (List<?>)mapINITClasses.get("Initialization");
            		timeyamlparse.setListINITClassesInitialization(listINITClassesInitialization);
            		List<?> listINITClassesFinalization = (List<?>)mapINITClasses.get("Finalization");
            		timeyamlparse.setListINITClassesFinalization(listINITClassesFinalization);
            		List<?> listINITClassesTotal = (List<?>)mapINITClasses.get("Total");
            		timeyamlparse.setListINITClassesTotal(listINITClassesTotal);
            	}/**end of the loading of INIT Classes: subclasses*/
            	Map<?,?> mapINITCategories = (Map<?,?>)mapINIT.get("Categories");
            	timeyamlparse.setMapINITCategories(mapINITCategories);
            	{/**start of the loading of INIT Categories: subclasses*/
            		Map<?,?> mapINITCategoriesWavefunction = (Map<?,?>)mapINITCategories.get("wavefunction");
            		timeyamlparse.setMapINITCategoriesWavefunction(mapINITCategoriesWavefunction);
            		{/**start of the loading of INIT Categories wavefunction : subclasses*/
                		List<?> listINITCategoriesWavefunctionData = (List<?>)mapINITCategoriesWavefunction.get("Data");
                		timeyamlparse.setListINITCategoriesWavefunctionData(listINITCategoriesWavefunctionData);
                		String stringINITCategoriesWavefunctionClass = (String)mapINITCategoriesWavefunction.get("Class");
                		timeyamlparse.setStringINITCategoriesWavefunctionClass(stringINITCategoriesWavefunctionClass);
                		String stringINITCategoriesWavefunctionInfo = (String)mapINITCategoriesWavefunction.get("Info");
                		timeyamlparse.setStringINITCategoriesWavefunctionInfo(stringINITCategoriesWavefunctionInfo);
            		}/**end of the loading of INIT Categories wavefunction : subclasses*/
            		Map<?,?> mapINITCategoriesApplyLocPotKin = (Map<?,?>)mapINITCategories.get("ApplyLocPotKin");
            		timeyamlparse.setMapINITCategoriesApplyLocPotKin(mapINITCategoriesApplyLocPotKin);
            		{/**start of the loading of INIT Categories ApplyLocPotKin : subclasses*/
            			List<?> listINITCategoriesApplyLocPotKinData = (List<?>)mapINITCategoriesApplyLocPotKin.get("Data");
            			timeyamlparse.setListINITCategoriesApplyLocPotKinData(listINITCategoriesApplyLocPotKinData);
            			String stringINITCategoriesApplyLocPotKinClass = (String)mapINITCategoriesApplyLocPotKin.get("Class");
            			timeyamlparse.setStringINITCategoriesApplyLocPotKinClass(stringINITCategoriesApplyLocPotKinClass);
            			String stringINITCategoriesApplyLocPotKinInfo = (String)mapINITCategoriesApplyLocPotKin.get("Info");
            			timeyamlparse.setStringINITCategoriesApplyLocPotKinInfo(stringINITCategoriesApplyLocPotKinInfo);
            		}/**end of the loading of INIT Categories ApplyLocPotKin : subclasses*/
            		Map<?,?> mapINITCategoriesApplyProj = (Map<?,?>)mapINITCategories.get("ApplyProj");
            		timeyamlparse.setMapINITCategoriesApplyProj(mapINITCategoriesApplyProj);
            		{/**start of the loading of INIT Categories  ApplyProj : subclasses*/
            			List<?> listINITCategoriesApplyProjData = (List<?>)mapINITCategoriesApplyProj.get("Data");
            			timeyamlparse.setListINITCategoriesApplyProjData(listINITCategoriesApplyProjData);
            			String stringINITCategoriesApplyProjClass = (String)mapINITCategoriesApplyProj.get("Class");
            			timeyamlparse.setStringINITCategoriesApplyProjClass(stringINITCategoriesApplyProjClass);
            			String stringINITCategoriesApplyProjInfo = (String)mapINITCategoriesApplyProj.get("Info");
            			timeyamlparse.setStringINITCategoriesApplyProjInfo(stringINITCategoriesApplyProjInfo);
            		}/**end of the loading of INIT Categories  ApplyProj : subclasses*/
            		Map<?,?> mapINITCategoriesCrtLocPot = (Map<?,?>)mapINITCategories.get("CrtLocPot");
            		timeyamlparse.setMapINITCategoriesCrtLocPot(mapINITCategoriesCrtLocPot);
            		{/**start of the loading of INIT Categories  CrtLocPot : subclasses*/
            			List<?> listINITCategoriesCrtLocPotData = (List<?>)mapINITCategoriesCrtLocPot.get("Data");
            			timeyamlparse.setListINITCategoriesCrtLocPotData(listINITCategoriesCrtLocPotData);
            			String stringINITCategoriesCrtLocPotClass = (String)mapINITCategoriesCrtLocPot.get("Class");
            			timeyamlparse.setStringINITCategoriesCrtLocPotClass(stringINITCategoriesCrtLocPotClass);
            			String stringINITCategoriesCrtLocPotInfo = (String)mapINITCategoriesCrtLocPot.get("Info");
            			timeyamlparse.setStringINITCategoriesCrtLocPotInfo(stringINITCategoriesCrtLocPotInfo);
            		}/**end of the loading of INIT Categories  CrtLocPot : subclasses*/
            		Map<?,?> mapINITCategoriesRho_comput = (Map<?,?>)mapINITCategories.get("Rho_comput");
            		timeyamlparse.setMapINITCategoriesRho_comput(mapINITCategoriesRho_comput);
            		{/**start of the loading of INIT Categories  Rho_comput : subclasses*/
            			List<?> listINITCategoriesRho_computData = (List<?>)mapINITCategoriesRho_comput.get("Data");
            			timeyamlparse.setListINITCategoriesRho_computData(listINITCategoriesRho_computData);
            			String stringINITCategoriesRho_computClass = (String)mapINITCategoriesRho_comput.get("Class");
            			timeyamlparse.setStringINITCategoriesRho_computClass(stringINITCategoriesRho_computClass);
            			String stringINITCategoriesRho_computInfo = (String)mapINITCategoriesRho_comput.get("Info");
            			timeyamlparse.setStringINITCategoriesRho_computInfo(stringINITCategoriesRho_computInfo);
            		}/**end of the loading of INIT Categories  Rho_comput : subclasses*/
            		Map<?,?> mapINITCategoriesUn_TransSwitch = (Map<?,?>)mapINITCategories.get("Un-TransSwitch");
            		timeyamlparse.setMapINITCategoriesUn_TransSwitch(mapINITCategoriesUn_TransSwitch);
            		{/**start of the loading of INIT Categories  Un-TransSwitch : subclasses*/
            			List<?> listINITCategoriesUn_TransSwitchData = (List<?>)mapINITCategoriesUn_TransSwitch.get("Data");
            			timeyamlparse.setListINITCategoriesUn_TransSwitchData(listINITCategoriesUn_TransSwitchData);
            			String stringINITCategoriesUn_TransSwitchClass = (String)mapINITCategoriesUn_TransSwitch.get("Class");
            			timeyamlparse.setStringINITCategoriesUn_TransSwitchClass(stringINITCategoriesUn_TransSwitchClass);
            			String stringINITCategoriesUn_TransSwitchInfo = (String)mapINITCategoriesUn_TransSwitch.get("Info");
            			timeyamlparse.setStringINITCategoriesUn_TransSwitchInfo(stringINITCategoriesUn_TransSwitchInfo);
            		}/**end of the loading of INIT Categories  Un-TransSwitch : subclasses*/
            		Map<?,?> mapINITCategoriesInput_comput = (Map<?,?>)mapINITCategories.get("Input_comput");
            		timeyamlparse.setMapINITCategoriesInput_comput(mapINITCategoriesInput_comput);
            		{/**start of the loading of INIT Categories  Input_comput : subclasses*/
            			List<?> listINITCategoriesInput_computData = (List<?>)mapINITCategoriesInput_comput.get("Data");
            			timeyamlparse.setListINITCategoriesInput_computData(listINITCategoriesInput_computData);
            			String stringINITCategoriesInput_computClass = (String)mapINITCategoriesInput_comput.get("Class");
            			timeyamlparse.setStringINITCategoriesInput_computClass(stringINITCategoriesInput_computClass);
            			String stringINITCategoriesInput_computInfo = (String)mapINITCategoriesInput_comput.get("Info");
            			timeyamlparse.setStringINITCategoriesInput_computInfo(stringINITCategoriesInput_computInfo);
            		}/**end of the loading of INIT Categories  Input_comput : subclasses*/
            		Map<?,?> mapINITCategoriesPSolv_comput = (Map<?,?>)mapINITCategories.get("PSolv_comput");
            		timeyamlparse.setMapINITCategoriesPSolv_comput(mapINITCategoriesPSolv_comput);
            		{/**start of the loading of INIT Categories  PSolv_comput : subclasses*/
            			List<?> listINITCategoriesPSolv_computData = (List<?>)mapINITCategoriesPSolv_comput.get("Data");
            			timeyamlparse.setListINITCategoriesPSolv_computData(listINITCategoriesPSolv_computData);
            			String stringINITCategoriesPSolv_computClass = (String)mapINITCategoriesPSolv_comput.get("Class");
            			timeyamlparse.setStringINITCategoriesPSolv_computClass(stringINITCategoriesPSolv_computClass);
            			String stringINITCategoriesPSolv_computInfo = (String)mapINITCategoriesPSolv_comput.get("Info");
            			timeyamlparse.setStringINITCategoriesPSolv_computInfo(stringINITCategoriesPSolv_computInfo);
            		}/**end of the loading of INIT Categories  PSolv_comput : subclasses*/
            		Map<?,?> mapINITCategoriesExchangecorr = (Map<?,?>)mapINITCategories.get("Exchangecorr");
            		timeyamlparse.setMapINITCategoriesExchangecorr(mapINITCategoriesExchangecorr);
            		{/**start of the loading of INIT Categories  Exchangecorr : subclasses*/
            			List<?> listINITCategoriesExchangecorrData = (List<?>)mapINITCategoriesExchangecorr.get("Data");
            			timeyamlparse.setListINITCategoriesExchangecorrData(listINITCategoriesExchangecorrData);
            			String stringINITCategoriesExchangecorrClass = (String)mapINITCategoriesExchangecorr.get("Class");
            			timeyamlparse.setStringINITCategoriesExchangecorrClass(stringINITCategoriesExchangecorrClass);
            			String stringINITCategoriesExchangecorrInfo = (String)mapINITCategoriesExchangecorr.get("Info");
            			timeyamlparse.setStringINITCategoriesExchangecorrInfo(stringINITCategoriesExchangecorrInfo);
            		}/**end of the loading of INIT Categories  Exchangecorr : subclasses*/
            		Map<?,?> mapINITCategoriesCrtDescriptors = (Map<?,?>)mapINITCategories.get("CrtDescriptors");
            		timeyamlparse.setMapINITCategoriesCrtDescriptors(mapINITCategoriesCrtDescriptors);
            		{/**start of the loading of INIT Categories  CrtDescriptors : subclasses*/
            			List<?> listINITCategoriesCrtDescriptorsData = (List<?>)mapINITCategoriesCrtDescriptors.get("Data");
            			timeyamlparse.setListINITCategoriesCrtDescriptorsData(listINITCategoriesCrtDescriptorsData);
            			String stringINITCategoriesCrtDescriptorsClass = (String)mapINITCategoriesCrtDescriptors.get("Class");
            			timeyamlparse.setStringINITCategoriesCrtDescriptorsClass(stringINITCategoriesCrtDescriptorsClass);
            			String stringINITCategoriesCrtDescriptorsInfo = (String)mapINITCategoriesCrtDescriptors.get("Info");
            			timeyamlparse.setStringINITCategoriesCrtDescriptorsInfo(stringINITCategoriesCrtDescriptorsInfo);
            		}/**end of the loading of INIT Categories  CrtDescriptors : subclasses*/
            		Map<?,?> mapINITCategoriesPSolvKernel = (Map<?,?>)mapINITCategories.get("PSolvKernel");
            		timeyamlparse.setMapINITCategoriesPSolvKernel(mapINITCategoriesPSolvKernel);
            		{/**start of the loading of INIT Categories  PSolvKernel : subclasses*/
            			List<?> listINITCategoriesPSolvKernelData = (List<?>)mapINITCategoriesPSolvKernel.get("Data");
            			timeyamlparse.setListINITCategoriesPSolvKernelData(listINITCategoriesPSolvKernelData);
            			String stringINITCategoriesPSolvKernelClass = (String)mapINITCategoriesPSolvKernel.get("Class");
            			timeyamlparse.setStringINITCategoriesPSolvKernelClass(stringINITCategoriesPSolvKernelClass);
            			String stringINITCategoriesPSolvKernelInfo = (String)mapINITCategoriesPSolvKernel.get("Info");
            			timeyamlparse.setStringINITCategoriesPSolvKernelInfo(stringINITCategoriesPSolvKernelInfo);
            		}/**end of the loading of INIT Categories  PSolvKernel : subclasses*/
            		Map<?,?> mapINITCategoriesPot_commun = (Map<?,?>)mapINITCategories.get("Pot_commun");
            		timeyamlparse.setMapINITCategoriesPot_commun(mapINITCategoriesPot_commun);
            		{/**start of the loading of INIT Categories  Pot_commun : subclasses*/
            			List<?> listINITCategoriesPot_communData = (List<?>)mapINITCategoriesPot_commun.get("Data");
            			timeyamlparse.setListINITCategoriesPot_communData(listINITCategoriesPot_communData);
            			String stringINITCategoriesPot_communClass = (String)mapINITCategoriesPot_commun.get("Class");
            			timeyamlparse.setStringINITCategoriesPot_communClass(stringINITCategoriesPot_communClass);
            			String stringINITCategoriesPot_communInfo = (String)mapINITCategoriesPot_commun.get("Info");
            			timeyamlparse.setStringINITCategoriesPot_communInfo(stringINITCategoriesPot_communInfo);
            		}/**end of the loading of INIT Categories  Pot_commun : subclasses*/
            	}
            }/**end of the loading of INIT Classes:*/
            
            //now loading the WFN_OPT section

            Map<?,?> mapWFN_OPT = (Map<?,?>)mapdataTimeYaml.get("WFN_OPT");
            timeyamlparse.setMapWFN_OPT(mapWFN_OPT);
            System.out.println(mapdataTimeYaml.get("WFN_OPT"));            
            {/**start of the loading of WFN_OPT Classes:*/
            	Map<?,?> mapWFN_OPTClasses = (Map<?,?>)mapWFN_OPT.get("Classes");
            	timeyamlparse.setMapWFN_OPTClasses(mapWFN_OPTClasses);
            	{/**start of the loading of WFN_OPT Classes: subclasses*/
            		List<?> listWFN_OPTClassesCommunications = (List<?>)mapWFN_OPTClasses.get("Communications");
            		timeyamlparse.setListWFN_OPTClassesCommunications(listWFN_OPTClassesCommunications);
            		List<?> listWFN_OPTClassesConvolutions = (List<?>)mapWFN_OPTClasses.get("Convolutions");
            		timeyamlparse.setListWFN_OPTClassesConvolutions(listWFN_OPTClassesConvolutions);
            		List<?> listWFN_OPTClassesLinearAlgebra = (List<?>)mapWFN_OPTClasses.get("Linear Algebra");
            		timeyamlparse.setListWFN_OPTClassesLinearAlgebra(listWFN_OPTClassesLinearAlgebra);
            		List<?> listWFN_OPTClassesOther = (List<?>)mapWFN_OPTClasses.get("Other");
            		timeyamlparse.setListWFN_OPTClassesOther(listWFN_OPTClassesOther);
            		List<?> listWFN_OPTClassesPotential = (List<?>)mapWFN_OPTClasses.get("Potential");
            		timeyamlparse.setListWFN_OPTClassesPotential(listWFN_OPTClassesPotential);
            		List<?> listWFN_OPTClassesInitialization = (List<?>)mapWFN_OPTClasses.get("Initialization");
            		timeyamlparse.setListWFN_OPTClassesInitialization(listWFN_OPTClassesInitialization);
            		List<?> listWFN_OPTClassesFinalization = (List<?>)mapWFN_OPTClasses.get("Finalization");
            		timeyamlparse.setListWFN_OPTClassesFinalization(listWFN_OPTClassesFinalization);
            		List<?> listWFN_OPTClassesTotal = (List<?>)mapWFN_OPTClasses.get("Total");
            		timeyamlparse.setListWFN_OPTClassesTotal(listWFN_OPTClassesTotal);
            	}/**end of the loading of WFN_OPT Classes: subclasses*/
            	Map<?,?> mapWFN_OPTCategories = (Map<?,?>)mapWFN_OPT.get("Categories");
            	timeyamlparse.setMapWFN_OPTCategories(mapWFN_OPTCategories);
            	{/**start of the loading of WFN_OPT Categories: subclasses*/
            		Map<?,?> mapWFN_OPTCategoriesPrecondition = (Map<?,?>)mapWFN_OPTCategories.get("Precondition");
            		timeyamlparse.setMapWFN_OPTCategoriesPrecondition(mapWFN_OPTCategoriesPrecondition);
            		{/**start of the loading of WFN_OPT Categories Precondition : subclasses*/
            			List<?> listWFN_OPTCategoriesPreconditionData = (List<?>)mapWFN_OPTCategoriesPrecondition.get("Data");
            			timeyamlparse.setListWFN_OPTCategoriesPreconditionData(listWFN_OPTCategoriesPreconditionData);
            			String stringWFN_OPTCategoriesPreconditionClass = (String)mapWFN_OPTCategoriesPrecondition.get("Class");
            			timeyamlparse.setStringWFN_OPTCategoriesPreconditionClass(stringWFN_OPTCategoriesPreconditionClass);
            			String stringWFN_OPTCategoriesPreconditionInfo = (String)mapWFN_OPTCategoriesPrecondition.get("Info");
            			timeyamlparse.setStringWFN_OPTCategoriesPreconditionInfo(stringWFN_OPTCategoriesPreconditionInfo);
            		}/**end of the loading of WFN_OPT Categories Precondition : subclasses*/
            		Map<?,?> mapWFN_OPTCategoriesApplyLocPotKin = (Map<?,?>)mapWFN_OPTCategories.get("ApplyLocPotKin");
            		timeyamlparse.setMapWFN_OPTCategoriesApplyLocPotKin(mapWFN_OPTCategoriesApplyLocPotKin);
            		{/**start of the loading of WFN_OPT Categories ApplyLocPotKin : subclasses*/
            			List<?> listWFN_OPTCategoriesApplyLocPotKinData = (List<?>)mapWFN_OPTCategoriesApplyLocPotKin.get("Data");
            			timeyamlparse.setListWFN_OPTCategoriesApplyLocPotKinData(listWFN_OPTCategoriesApplyLocPotKinData);
            			String stringWFN_OPTCategoriesApplyLocPotKinClass = (String)mapWFN_OPTCategoriesApplyLocPotKin.get("Class");
            			timeyamlparse.setStringWFN_OPTCategoriesApplyLocPotKinClass(stringWFN_OPTCategoriesApplyLocPotKinClass);
            			String stringWFN_OPTCategoriesApplyLocPotKinInfo = (String)mapWFN_OPTCategoriesApplyLocPotKin.get("Info");
            			timeyamlparse.setStringWFN_OPTCategoriesApplyLocPotKinInfo(stringWFN_OPTCategoriesApplyLocPotKinInfo);
            		}/**end of the loading of WFN_OPT Categories ApplyLocPotKin : subclasses*/
            		Map<?,?> mapWFN_OPTCategoriesApplyProj = (Map<?,?>)mapWFN_OPTCategories.get("ApplyProj");
            		timeyamlparse.setMapWFN_OPTCategoriesApplyProj(mapWFN_OPTCategoriesApplyProj);
            		{/**start of the loading of WFN_OPT Categories ApplyProj : subclasses*/
            			List<?> listWFN_OPTCategoriesApplyProjData = (List<?>)mapWFN_OPTCategoriesApplyProj.get("Data");
            			timeyamlparse.setListWFN_OPTCategoriesApplyProjData(listWFN_OPTCategoriesApplyProjData);
            			String stringWFN_OPTCategoriesApplyProjClass = (String)mapWFN_OPTCategoriesApplyProj.get("Class");
            			timeyamlparse.setStringWFN_OPTCategoriesApplyProjClass(stringWFN_OPTCategoriesApplyProjClass);
            			String stringWFN_OPTCategoriesApplyProjInfo = (String)mapWFN_OPTCategoriesApplyProj.get("Info");
            			timeyamlparse.setStringWFN_OPTCategoriesApplyProjInfo(stringWFN_OPTCategoriesApplyProjInfo);
            		}/**end of the loading of WFN_OPT Categories ApplyProj : subclasses*/
            		Map<?,?> mapWFN_OPTCategoriesUn_TransSwitch = (Map<?,?>)mapWFN_OPTCategories.get("Un-TransSwitch");
            		timeyamlparse.setMapWFN_OPTCategoriesUn_TransSwitch(mapWFN_OPTCategoriesUn_TransSwitch);
            		{/**start of the loading of WFN_OPT Categories Un-TransSwitch : subclasses*/
            			List<?> listWFN_OPTCategoriesUn_TransSwitchData = (List<?>)mapWFN_OPTCategoriesUn_TransSwitch.get("Data");
            			timeyamlparse.setListWFN_OPTCategoriesUn_TransSwitchData(listWFN_OPTCategoriesUn_TransSwitchData);
            			String stringWFN_OPTCategoriesUn_TransSwitchClass = (String)mapWFN_OPTCategoriesUn_TransSwitch.get("class");
            			timeyamlparse.setStringWFN_OPTCategoriesUn_TransSwitchClass(stringWFN_OPTCategoriesUn_TransSwitchClass);
            			String stringWFN_OPTCategoriesUn_TransSwitchInfo = (String)mapWFN_OPTCategoriesUn_TransSwitch.get("Info");
            			timeyamlparse.setStringWFN_OPTCategoriesUn_TransSwitchInfo(stringWFN_OPTCategoriesUn_TransSwitchInfo);
            		}/**end of the loading of WFN_OPT Categories Un-TransSwitch : subclasses*/
            		Map<?,?> mapWFN_OPTCategoriesLagrM_comput = (Map<?,?>)mapWFN_OPTCategories.get("LagrM_comput");
            		timeyamlparse.setMapWFN_OPTCategoriesLagrM_comput(mapWFN_OPTCategoriesLagrM_comput);
            		{/**start of the loading of WFN_OPT Categories LagrM_comput : subclasses*/
            			List<?> listWFN_OPTCategoriesLagrM_computData = (List<?>)mapWFN_OPTCategoriesLagrM_comput.get("Data");
            			timeyamlparse.setListWFN_OPTCategoriesLagrM_computData(listWFN_OPTCategoriesLagrM_computData);
            			String stringWFN_OPTCategoriesLagrM_computClass = (String)mapWFN_OPTCategoriesLagrM_comput.get("Class");
            			timeyamlparse.setStringWFN_OPTCategoriesLagrM_computClass(stringWFN_OPTCategoriesLagrM_computClass);
            			String stringWFN_OPTCategoriesLagrM_computInfo = (String)mapWFN_OPTCategoriesLagrM_comput.get("Info");
            			timeyamlparse.setStringWFN_OPTCategoriesLagrM_computInfo(stringWFN_OPTCategoriesLagrM_computInfo);
            		}/**end of the loading of WFN_OPT Categories LagrM_comput : subclasses*/
            		Map<?,?> mapWFN_OPTCategoriesRho_comput = (Map<?,?>)mapWFN_OPTCategories.get("Rho_comput");
            		timeyamlparse.setMapWFN_OPTCategoriesRho_comput(mapWFN_OPTCategoriesRho_comput);
            		{/**start of the loading of WFN_OPT Categories Rho_comput : subclasses*/
            			List<?> listWFN_OPTCategoriesRho_computData = (List<?>)mapWFN_OPTCategoriesRho_comput.get("Data");
            			timeyamlparse.setListWFN_OPTCategoriesRho_computData(listWFN_OPTCategoriesRho_computData);
            			String stringWFN_OPTCategoriesRho_computClass = (String)mapWFN_OPTCategoriesRho_comput.get("Class");
            			timeyamlparse.setStringWFN_OPTCategoriesRho_computClass(stringWFN_OPTCategoriesRho_computClass);
            			String stringWFN_OPTCategoriesRho_computInfo = (String)mapWFN_OPTCategoriesRho_comput.get("Info");
            			timeyamlparse.setStringWFN_OPTCategoriesRho_computInfo(stringWFN_OPTCategoriesRho_computInfo);
            		}/**end of the loading of WFN_OPT Categories Rho_comput : subclasses*/
            		Map<?,?> mapWFN_OPTCategoriesChol_comput = (Map<?,?>)mapWFN_OPTCategories.get("Chol_comput");
            		timeyamlparse.setMapWFN_OPTCategoriesChol_comput(mapWFN_OPTCategoriesChol_comput);
            		{/**start of the loading of WFN_OPT Categories Chol_comput : subclasses*/
            			List<?> listWFN_OPTCategoriesChol_computData = (List<?>)mapWFN_OPTCategoriesChol_comput.get("Data");
            			timeyamlparse.setListWFN_OPTCategoriesChol_computData(listWFN_OPTCategoriesChol_computData);
            			String stringWFN_OPTCategoriesChol_computClass = (String)mapWFN_OPTCategoriesChol_comput.get("Class");
            			timeyamlparse.setStringWFN_OPTCategoriesChol_computClass(stringWFN_OPTCategoriesChol_computClass);
            			String stringWFN_OPTCategoriesChol_computInfo = (String)mapWFN_OPTCategoriesChol_comput.get("Info");
            			timeyamlparse.setStringWFN_OPTCategoriesChol_computInfo(stringWFN_OPTCategoriesChol_computInfo);
            		}/**end of the loading of WFN_OPT Categories Chol_comput : subclasses*/
            		Map<?,?> mapWFN_OPTCategoriesDiis = (Map<?,?>)mapWFN_OPTCategories.get("Diis");
            		timeyamlparse.setMapWFN_OPTCategoriesDiis(mapWFN_OPTCategoriesDiis);
            		{/**start of the loading of WFN_OPT Categories Diis : subclasses*/
            			List<?> listWFN_OPTCategoriesDiisData = (List<?>)mapWFN_OPTCategoriesDiis.get("Data");
            			timeyamlparse.setListWFN_OPTCategoriesDiisData(listWFN_OPTCategoriesDiisData);
            			String stringWFN_OPTCategoriesDiisClass = (String)mapWFN_OPTCategoriesDiis.get("Class");
            			timeyamlparse.setStringWFN_OPTCategoriesDiisClass(stringWFN_OPTCategoriesDiisClass);
            			String stringWFN_OPTCategoriesDiisInfo = (String)mapWFN_OPTCategoriesDiis.get("Info");
            			timeyamlparse.setStringWFN_OPTCategoriesDiisInfo(stringWFN_OPTCategoriesDiisInfo);
            		}/**end of the loading of WFN_OPT Categories Diis : subclasses*/
            		Map<?,?> mapWFN_OPTCategoriesPSolv_comput = (Map<?,?>)mapWFN_OPTCategories.get("PSolv_comput");
            		timeyamlparse.setMapWFN_OPTCategoriesPSolv_comput(mapWFN_OPTCategoriesPSolv_comput);
            		{/**start of the loading of WFN_OPT Categories PSolv_comput : subclasses*/
            			List<?> listWFN_OPTCategoriesPSolv_computData = (List<?>)mapWFN_OPTCategoriesPSolv_comput.get("Data");
            			timeyamlparse.setListWFN_OPTCategoriesPSolv_computData(listWFN_OPTCategoriesPSolv_computData);
            			String stringWFN_OPTCategoriesPSolv_computClass = (String)mapWFN_OPTCategoriesPSolv_comput.get("Class");
            			timeyamlparse.setStringWFN_OPTCategoriesPSolv_computClass(stringWFN_OPTCategoriesPSolv_computClass);
            			String stringWFN_OPTCategoriesPSolv_computInfo = (String)mapWFN_OPTCategoriesPSolv_comput.get("Info");
            			timeyamlparse.setStringWFN_OPTCategoriesPSolv_computInfo(stringWFN_OPTCategoriesPSolv_computInfo);
            		}/**end of the loading of WFN_OPT Categories PSolv_comput : subclasses*/
            		Map<?,?> mapWFN_OPTCategoriesExchangecorr = (Map<?,?>)mapWFN_OPTCategories.get("Exchangecorr");
            		timeyamlparse.setMapWFN_OPTCategoriesExchangecorr(mapWFN_OPTCategoriesExchangecorr);
            		{/**start of the loading of WFN_OPT Categories Exchangecorr : subclasses*/
            			List<?> listWFN_OPTCategoriesExchangecorrData = (List<?>)mapWFN_OPTCategoriesExchangecorr.get("Data");
            			timeyamlparse.setListWFN_OPTCategoriesExchangecorrData(listWFN_OPTCategoriesExchangecorrData);
            			String stringWFN_OPTCategoriesExchangecorrClass = (String)mapWFN_OPTCategoriesExchangecorr.get("Class");
            			timeyamlparse.setStringWFN_OPTCategoriesExchangecorrClass(stringWFN_OPTCategoriesExchangecorrClass);
            			String stringWFN_OPTCategoriesExchangecorrInfo = (String)mapWFN_OPTCategoriesExchangecorr.get("Info");
            			timeyamlparse.setStringWFN_OPTCategoriesExchangecorrInfo(stringWFN_OPTCategoriesExchangecorrInfo);
            		}/**start of the loading of WFN_OPT Categories Exchangecorr : subclasses*/
            		Map<?,?> mapWFN_OPTCategoriesPot_commun = (Map<?,?>)mapWFN_OPTCategories.get("Pot_commun");
            		timeyamlparse.setMapWFN_OPTCategoriesPot_commun(mapWFN_OPTCategoriesPot_commun);
            		{/**start of the loading of WFN_OPT Categories Pot_commun : subclasses*/
            			List<?> listWFN_OPTCategoriesPot_communData = (List<?>)mapWFN_OPTCategoriesPot_commun.get("Data");
            			timeyamlparse.setListWFN_OPTCategoriesPot_communData(listWFN_OPTCategoriesPot_communData);
            			String stringWFN_OPTCategoriesPot_communClass = (String)mapWFN_OPTCategoriesPot_commun.get("Class");
            			timeyamlparse.setStringWFN_OPTCategoriesPot_communClass(stringWFN_OPTCategoriesPot_communClass);
            			String stringWFN_OPTCategoriesPot_communInfo = (String)mapWFN_OPTCategoriesPot_commun.get("Info");
            			timeyamlparse.setStringWFN_OPTCategoriesPot_communInfo(stringWFN_OPTCategoriesPot_communInfo);
            		}/**end of the loading of WFN_OPT Categories Pot_commun : subclasses*/
            	}/**end of the loading of WFN_OPT Categories: subclasses*/
            }/**end of the loading of WFN_OPT Classes:*/

            //now loading the LAST section

            Map<?,?> mapLAST = (Map<?,?>)mapdataTimeYaml.get("LAST");
            timeyamlparse.setMapLAST(mapLAST);
            System.out.println(mapdataTimeYaml.get("LAST"));            
            {/**start of the loading of LAST Classes:*/
            	Map<?,?> mapLASTClasses = (Map<?,?>)mapLAST.get("Classes");
            	timeyamlparse.setMapLASTClasses(mapLASTClasses);
            	{/**start of the loading of LAST Classes: subclasses*/
            		List<?> listLASTClassesCommunications = (List<?>)mapLASTClasses.get("Communications");
            		timeyamlparse.setListLASTClassesCommunications(listLASTClassesCommunications);
            		List<?> listLASTClassesConvolutions = (List<?>)mapLASTClasses.get("Convolutions");
            		timeyamlparse.setListLASTClassesConvolutions(listLASTClassesConvolutions);
            		List<?> listLASTClassesLinearAlgebra = (List<?>)mapLASTClasses.get("Linear Algebra");
            		timeyamlparse.setListLASTClassesLinearAlgebra(listLASTClassesLinearAlgebra);
            		List<?> listLASTClassesOther = (List<?>)mapLASTClasses.get("Other");
            		timeyamlparse.setListLASTClassesOther(listLASTClassesOther);
            		List<?> listLASTClassesPotential = (List<?>)mapLASTClasses.get("Potential");
            		timeyamlparse.setListLASTClassesPotential(listLASTClassesPotential);
            		List<?> listLASTClassesInitialization = (List<?>)mapLASTClasses.get("Initialization");
            		timeyamlparse.setListLASTClassesInitialization(listLASTClassesInitialization);
            		List<?> listLASTClassesFinalization = (List<?>)mapLASTClasses.get("Finalization");
            		timeyamlparse.setListLASTClassesFinalization(listLASTClassesFinalization);
            		List<?> listLASTClassesTotal = (List<?>)mapLASTClasses.get("Total");
            		timeyamlparse.setListLASTClassesTotal(listLASTClassesTotal);
            	}/**end of the loading of LAST Classes: subclasses*/
            	Map<?,?> mapLASTCategories = (Map<?,?>)mapLAST.get("Categories");
            	timeyamlparse.setMapLASTCategories(mapLASTCategories);
            	{/**start of the loading of LAST Categories: subclasses*/
            		Map<?,?> mapLASTCategoriesForces = (Map<?,?>)mapLASTCategories.get("Forces");
            		timeyamlparse.setMapLASTCategoriesForces(mapLASTCategoriesForces);
            		{/**start of the loading of LAST Categories Forces : subclasses*/
            			List<?> listLASTCategoriesForcesData = (List<?>)mapLASTCategoriesForces.get("Data");
            			timeyamlparse.setListLASTCategoriesForcesData(listLASTCategoriesForcesData);
            			String stringLASTCategoriesForcesClass = (String)mapLASTCategoriesForces.get("Class");
            			timeyamlparse.setStringLASTCategoriesForcesClass(stringLASTCategoriesForcesClass);
            			String stringLASTCategoriesForcesInfo = (String)mapLASTCategoriesForces.get("Info");
            			timeyamlparse.setStringLASTCategoriesForcesInfo(stringLASTCategoriesForcesInfo);
            		}/**end of the loading of LAST Categories Forces : subclasses*/
            		Map<?,?> mapLASTCategoriesRho_comput = (Map<?,?>)mapLASTCategories.get("Rho_comput");
            		timeyamlparse.setMapLASTCategoriesRho_comput(mapLASTCategoriesRho_comput);
            		{/**start of the loading of LAST Categories Rho_comput : subclasses*/
            			List<?> listLASTCategoriesRho_computData = (List<?>)mapLASTCategoriesRho_comput.get("Data");
            			timeyamlparse.setListLASTCategoriesRho_computData(listLASTCategoriesRho_computData);
            			String stringLASTCategoriesRho_computClass = (String)mapLASTCategoriesRho_comput.get("Class");
            			timeyamlparse.setStringLASTCategoriesRho_computClass(stringLASTCategoriesRho_computClass);
            			String stringLASTCategoriesRho_computInfo = (String)mapLASTCategoriesRho_comput.get("Info");
            			timeyamlparse.setStringLASTCategoriesRho_computInfo(stringLASTCategoriesRho_computInfo);
            		}/**end of the loading of LAST Categories Rho_comput : subclasses*/
            		Map<?,?> mapLASTCategoriesPSolv_comput = (Map<?,?>)mapLASTCategories.get("PSolv_comput");
            		timeyamlparse.setMapLASTCategoriesPSolv_comput(mapLASTCategoriesPSolv_comput);
            		{/**start of the loading of LAST Categories PSolv_comput : subclasses*/
            			List<?> listLASTCategoriesPSolv_computData = (List<?>)mapLASTCategoriesPSolv_comput.get("Data");
            			timeyamlparse.setListLASTCategoriesPSolv_computData(listLASTCategoriesPSolv_computData);
            			String stringLASTCategoriesPSolv_computClass = (String)mapLASTCategoriesPSolv_comput.get("Class");
            			timeyamlparse.setStringLASTCategoriesPSolv_computClass(stringLASTCategoriesPSolv_computClass);
            			String stringLASTCategoriesPSolv_computInfo = (String)mapLASTCategoriesPSolv_comput.get("Info");
            			timeyamlparse.setStringLASTCategoriesPSolv_computInfo(stringLASTCategoriesPSolv_computInfo);
            		}/**end of the loading of LAST Categories PSolv_comput : subclasses*/
            	}/**end of the loading of LAST Categories: subclasses*/
            }/**end of the loading of LAST Classes:*/

            //now loading the SUMMARY section

            Map<?,?> mapSUMMARY = (Map<?,?>)mapdataTimeYaml.get("SUMMARY");
            timeyamlparse.setMapSUMMARY(mapSUMMARY);
            System.out.println(mapdataTimeYaml.get("SUMMARY"));            
            {/**start of the loading of SUMMARY Classes:*/
            	List<?> listSUMMARYINIT	= (List<?>)mapSUMMARY.get("INIT");
            	timeyamlparse.setListSUMMARYINIT(listSUMMARYINIT);
            	List<?> listSUMMARYWFN_OPT	= (List<?>)mapSUMMARY.get("WFN_OPT");
            	timeyamlparse.setListSUMMARYWFN_OPT(listSUMMARYWFN_OPT);
            	List<?> listSUMMARYLAST	= (List<?>)mapSUMMARY.get("LAST");
            	timeyamlparse.setListSUMMARYLAST(listSUMMARYLAST);
            	List<?> listSUMMARYTotal	= (List<?>)mapSUMMARY.get("Total");
            	timeyamlparse.setListSUMMARYTotal(listSUMMARYTotal);
            }/**end of the loading of SUMMARY Classes:*/
            Map<?,?> mapSUMMARYCPUParallelism = (Map<?,?>)mapSUMMARY.get("CPU Parallelism");
            timeyamlparse.setMapSUMMARYCPUParallelism(mapSUMMARYCPUParallelism);

		} finally {
			log.info(filename+" did yaml.load(input).");
		}
	}
	/**
	 * Setters and getters for the class TimeYamlParse
	 */
	/**
	 * @return the path
	 */
	public String getPath() {
		return path;
	}
	/**
	 * @param path the path to set
	 */
	public void setPath(String path) {
		this.path = path;
	}
	/**
	 * @return the filename
	 */
	public String getFilename() {
		return filename;
	}
	/**
	 * @param filename the filename to set
	 */
	public void setFilename(String filename) {
		this.filename = filename;
	}
	/**
	 * @return the timeyamlparse
	 */
	public TimeYamlVariables getTimeyamlparse() {
		return timeyamlparse;
	}
	/**
	 * @param timeyamlparse the timeyamlparse to set
	 */
	public void setTimeyamlparse(TimeYamlVariables timeyamlparse) {
		this.timeyamlparse = timeyamlparse;
	}
}
