/**
 * 
 */
package org.openmolgrid.format.yaml;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.yaml.snakeyaml.Yaml;

/**
 * class to parse the Time.yaml file.
 * @author Frederic Bonnet
 *
 */
public class TimeYamlVariables {
	private Yaml yaml;
	//main maps ordered
	private Map<?,?> mapINIT;
	private Map<?,?> mapWFN_OPT;
	private Map<?,?> mapLAST;
	private Map<?,?> mapSUMMARY;
	//sub maps ordered 
	private Map<?,?> mapINITClasses;
	//sub sub maps ordered
	private List<?> listINITClassesCommunications;
	private List<?> listINITClassesConvolutions;
	private List<?> listINITClassesLinearAlgebra;
	private List<?> listINITClassesOther;
	private List<?> listINITClassesPotential;
	private List<?> listINITClassesInitialization;
	private List<?> listINITClassesFinalization;
	private List<?> listINITClassesTotal;
	//sub maps ordered 
	private Map<?,?> mapINITCategories;
	//sub sub maps ordered
	private Map<?,?> mapINITCategoriesWavefunction;
			//sub sub sub maps ordered
    		private List<?> listINITCategoriesWavefunctionData;
    		private String stringINITCategoriesWavefunctionClass;
    		private String stringINITCategoriesWavefunctionInfo;
	private Map<?,?> mapINITCategoriesApplyLocPotKin;
			//sub sub sub maps ordered
    		private List<?> listINITCategoriesApplyLocPotKinData;
    		private String stringINITCategoriesApplyLocPotKinClass;
    		private String stringINITCategoriesApplyLocPotKinInfo;
	private Map<?,?> mapINITCategoriesApplyProj;
			//sub sub sub maps ordered
    		private List<?> listINITCategoriesApplyProjData;
    		private String stringINITCategoriesApplyProjClass;
    		private String stringINITCategoriesApplyProjInfo;
	private Map<?,?> mapINITCategoriesCrtLocPot;
			//sub sub sub maps ordered
    		private List<?> listINITCategoriesCrtLocPotData;
    		private String stringINITCategoriesCrtLocPotClass;
    		private String stringINITCategoriesCrtLocPotInfo;
	private Map<?,?> mapINITCategoriesRho_comput;
			//sub sub sub maps ordered
    		private List<?> listINITCategoriesRho_computData;
    		private String stringINITCategoriesRho_computClass;
    		private String stringINITCategoriesRho_computInfo;
	private Map<?,?> mapINITCategoriesUn_TransSwitch;
			//sub sub sub maps ordered
    		private List<?> listINITCategoriesUn_TransSwitchData;
    		private String stringINITCategoriesUn_TransSwitchClass;
    		private String stringINITCategoriesUn_TransSwitchInfo;
	private Map<?,?> mapINITCategoriesInput_comput;
			//sub sub sub maps ordered
    		private List<?> listINITCategoriesInput_computData;
    		private String stringINITCategoriesInput_computClass;
    		private String stringINITCategoriesInput_computInfo;
	private Map<?,?> mapINITCategoriesPSolv_comput;
			//sub sub sub maps ordered
    		private List<?> listINITCategoriesPSolv_computData;
    		private String stringINITCategoriesPSolv_computClass;
    		private String stringINITCategoriesPSolv_computInfo;
	private Map<?,?> mapINITCategoriesExchangecorr;
			//sub sub sub maps ordered
    		private List<?> listINITCategoriesExchangecorrData;
    		private String stringINITCategoriesExchangecorrClass;
    		private String stringINITCategoriesExchangecorrInfo;
	private Map<?,?> mapINITCategoriesCrtDescriptors;
			//sub sub sub maps ordered
    		private List<?> listINITCategoriesCrtDescriptorsData;
    		private String stringINITCategoriesCrtDescriptorsClass;
    		private String stringINITCategoriesCrtDescriptorsInfo;
	private Map<?,?> mapINITCategoriesPSolvKernel;
			//sub sub sub maps ordered
    		private List<?> listINITCategoriesPSolvKernelData;
    		private String stringINITCategoriesPSolvKernelClass;
    		private String stringINITCategoriesPSolvKernelInfo;
	private Map<?,?> mapINITCategoriesPot_commun;
			//sub sub sub maps ordered
    		private List<?> listINITCategoriesPot_communData;
    		private String stringINITCategoriesPot_communClass;
    		private String stringINITCategoriesPot_communInfo;
	//sub maps ordered 
	private Map<?,?> mapWFN_OPTClasses;
	//sub sub maps ordered
	        private List<?> listWFN_OPTClassesCommunications;
	        private List<?> listWFN_OPTClassesConvolutions;
	        private List<?> listWFN_OPTClassesLinearAlgebra;
	        private List<?> listWFN_OPTClassesOther;
	        private List<?> listWFN_OPTClassesPotential;
	        private List<?> listWFN_OPTClassesInitialization;
	        private List<?> listWFN_OPTClassesFinalization;
	        private List<?> listWFN_OPTClassesTotal;
   	//sub maps ordered 
   	private Map<?,?> mapWFN_OPTCategories;
	//sub sub maps ordered
    private Map<?,?> mapWFN_OPTCategoriesPrecondition;
            //sub sub sub maps ordered
            private List<?> listWFN_OPTCategoriesPreconditionData;
            private String stringWFN_OPTCategoriesPreconditionClass;
            private String stringWFN_OPTCategoriesPreconditionInfo;
    private Map<?,?> mapWFN_OPTCategoriesApplyLocPotKin;
            //sub sub sub maps ordered
            private List<?> listWFN_OPTCategoriesApplyLocPotKinData;
            private String stringWFN_OPTCategoriesApplyLocPotKinClass;
            private String stringWFN_OPTCategoriesApplyLocPotKinInfo;
    private Map<?,?> mapWFN_OPTCategoriesApplyProj;
            //sub sub sub maps ordered
            private List<?> listWFN_OPTCategoriesApplyProjData;
            private String stringWFN_OPTCategoriesApplyProjClass;
            private String stringWFN_OPTCategoriesApplyProjInfo;
    private Map<?,?> mapWFN_OPTCategoriesUn_TransSwitch;
            //sub sub sub maps ordered
            private List<?> listWFN_OPTCategoriesUn_TransSwitchData;
            private String stringWFN_OPTCategoriesUn_TransSwitchClass;
            private String stringWFN_OPTCategoriesUn_TransSwitchInfo;
    private Map<?,?> mapWFN_OPTCategoriesLagrM_comput;
            //sub sub sub maps ordered
            private List<?> listWFN_OPTCategoriesLagrM_computData;
            private String stringWFN_OPTCategoriesLagrM_computClass;
            private String stringWFN_OPTCategoriesLagrM_computInfo;
    private Map<?,?> mapWFN_OPTCategoriesRho_comput;
            //sub sub sub maps ordered
            private List<?> listWFN_OPTCategoriesRho_computData;
            private String stringWFN_OPTCategoriesRho_computClass;
            private String stringWFN_OPTCategoriesRho_computInfo;
    private Map<?,?> mapWFN_OPTCategoriesChol_comput;
            //sub sub sub maps ordered
            private List<?> listWFN_OPTCategoriesChol_computData;
            private String stringWFN_OPTCategoriesChol_computClass;
            private String stringWFN_OPTCategoriesChol_computInfo;
    private Map<?,?> mapWFN_OPTCategoriesDiis;
            //sub sub sub maps ordered
            private List<?> listWFN_OPTCategoriesDiisData;
            private String stringWFN_OPTCategoriesDiisClass;
            private String stringWFN_OPTCategoriesDiisInfo;
    private Map<?,?> mapWFN_OPTCategoriesPSolv_comput;
            //sub sub sub maps ordered
            private List<?> listWFN_OPTCategoriesPSolv_computData;
            private String stringWFN_OPTCategoriesPSolv_computClass;
            private String stringWFN_OPTCategoriesPSolv_computInfo;
    private Map<?,?> mapWFN_OPTCategoriesExchangecorr;
            //sub sub sub maps ordered
            private List<?> listWFN_OPTCategoriesExchangecorrData;
            private String stringWFN_OPTCategoriesExchangecorrClass;
            private String stringWFN_OPTCategoriesExchangecorrInfo;
    private Map<?,?> mapWFN_OPTCategoriesPot_commun;
            //sub sub sub maps ordered
            private List<?> listWFN_OPTCategoriesPot_communData;
            private String stringWFN_OPTCategoriesPot_communClass;
            private String stringWFN_OPTCategoriesPot_communInfo;
    //sub maps ordered 
    private Map<?,?> mapLASTClasses;
            //sub sub maps ordered
            private List<?> listLASTClassesCommunications;
            private List<?> listLASTClassesConvolutions;
            private List<?> listLASTClassesLinearAlgebra;
            private List<?> listLASTClassesOther;
            private List<?> listLASTClassesPotential;
            private List<?> listLASTClassesInitialization;
            private List<?> listLASTClassesFinalization;
            private List<?> listLASTClassesTotal;
   	//sub maps ordered 
   	private Map<?,?> mapLASTCategories;
	//sub sub maps ordered
    private Map<?,?> mapLASTCategoriesForces;
            //sub sub sub maps ordered
            private List<?> listLASTCategoriesForcesData;
            private String stringLASTCategoriesForcesClass;
            private String stringLASTCategoriesForcesInfo;
    private Map<?,?> mapLASTCategoriesRho_comput;
            //sub sub sub maps ordered
            private List<?> listLASTCategoriesRho_computData;
            private String stringLASTCategoriesRho_computClass;
            private String stringLASTCategoriesRho_computInfo;
    private Map<?,?> mapLASTCategoriesPSolv_comput;
            //sub sub sub maps ordered
            private List<?> listLASTCategoriesPSolv_computData;
            private String stringLASTCategoriesPSolv_computClass;
            private String stringLASTCategoriesPSolv_computInfo;
    //sub maps ordered 
    private List<?> listSUMMARYINIT;
    private List<?> listSUMMARYWFN_OPT;
    private List<?> listSUMMARYLAST;
    private List<?> listSUMMARYTotal;
    private Map<?,?> mapSUMMARYCPUParallelism;
    		//sub sub sub maps ordered
    		private String stringSUMMARYCPUParallelismMPIprocs;
    		private String stringSUMMARYCPUParallelismOMPthrds;
	
	/** HashMap which stores all the executed modules as key and the success of the module run as value */
	private HashMap<String, Boolean> executedModules;

	/**
	 * default constructor
	 */
	public TimeYamlVariables() {
		executedModules = new HashMap<String, Boolean>();
	}
	/**
	 * constructor with {@link Yaml} instance
	 * @param yaml
	 */
	public TimeYamlVariables(Yaml yaml) {
		this.yaml = yaml;
	}
	/**
	 * constructor 
	 * @param args the application parameters
	 */
	public TimeYamlVariables(String[] args) {
		 yaml = new Yaml();
	}
	
	//The setters and getters for the private variables
	/**
	 * @return the yaml
	 */
	public Yaml getYaml() {
		return yaml;
	}
	/**
	 * @param yaml the yaml to set
	 */
	public void setYaml(Yaml yaml) {
		this.yaml = yaml;
	}
	/**--------------------------------------------------------------*/	
	/**------ The setters and getters for all of the variables ------*/	
	/**--------------------------------------------------------------*/	
	/**
	 * @return the mapINIT
	 */
	public Map<?, ?> getMapINIT() {
		return mapINIT;
	}
	/**
	 * @param mapINIT the mapINIT to set
	 */
	public void setMapINIT(Map<?, ?> mapINIT) {
		this.mapINIT = mapINIT;
	}
	/**
	 * @return the mapWFN_OPT
	 */
	public Map<?, ?> getMapWFN_OPT() {
		return mapWFN_OPT;
	}
	/**
	 * @param mapWFN_OPT the mapWFN_OPT to set
	 */
	public void setMapWFN_OPT(Map<?, ?> mapWFN_OPT) {
		this.mapWFN_OPT = mapWFN_OPT;
	}
	/**
	 * @return the mapLAST
	 */
	public Map<?, ?> getMapLAST() {
		return mapLAST;
	}
	/**
	 * @param mapLAST the mapLAST to set
	 */
	public void setMapLAST(Map<?, ?> mapLAST) {
		this.mapLAST = mapLAST;
	}
	/**
	 * @return the mapSUMMARY
	 */
	public Map<?, ?> getMapSUMMARY() {
		return mapSUMMARY;
	}
	/**
	 * @param mapSUMMARY the mapSUMMARY to set
	 */
	public void setMapSUMMARY(Map<?, ?> mapSUMMARY) {
		this.mapSUMMARY = mapSUMMARY;
	}
	/**
	 * @return the mapINITClasses
	 */
	public Map<?, ?> getMapINITClasses() {
		return mapINITClasses;
	}
	/**
	 * @param mapINITClasses the mapINITClasses to set
	 */
	public void setMapINITClasses(Map<?, ?> mapINITClasses) {
		this.mapINITClasses = mapINITClasses;
	}
	/**
	 * @return the listINITClassesCommunications
	 */
	public List<?> getListINITClassesCommunications() {
		return listINITClassesCommunications;
	}
	/**
	 * @param listINITClassesCommunications the listINITClassesCommunications to set
	 */
	public void setListINITClassesCommunications(
			List<?> listINITClassesCommunications) {
		this.listINITClassesCommunications = listINITClassesCommunications;
	}
	/**
	 * @return the listINITClassesConvolutions
	 */
	public List<?> getListINITClassesConvolutions() {
		return listINITClassesConvolutions;
	}
	/**
	 * @param listINITClassesConvolutions the listINITClassesConvolutions to set
	 */
	public void setListINITClassesConvolutions(List<?> listINITClassesConvolutions) {
		this.listINITClassesConvolutions = listINITClassesConvolutions;
	}
	/**
	 * @return the listINITClassesLinearAlgebra
	 */
	public List<?> getListINITClassesLinearAlgebra() {
		return listINITClassesLinearAlgebra;
	}
	/**
	 * @param listINITClassesLinearAlgebra the listINITClassesLinearAlgebra to set
	 */
	public void setListINITClassesLinearAlgebra(List<?> listINITClassesLinearAlgebra) {
		this.listINITClassesLinearAlgebra = listINITClassesLinearAlgebra;
	}
	/**
	 * @return the listINITClassesOther
	 */
	public List<?> getListINITClassesOther() {
		return listINITClassesOther;
	}
	/**
	 * @param listINITClassesOther the listINITClassesOther to set
	 */
	public void setListINITClassesOther(List<?> listINITClassesOther) {
		this.listINITClassesOther = listINITClassesOther;
	}
	/**
	 * @return the listINITClassesPotential
	 */
	public List<?> getListINITClassesPotential() {
		return listINITClassesPotential;
	}
	/**
	 * @param listINITClassesPotential the listINITClassesPotential to set
	 */
	public void setListINITClassesPotential(List<?> listINITClassesPotential) {
		this.listINITClassesPotential = listINITClassesPotential;
	}
	/**
	 * @return the listINITClassesInitialization
	 */
	public List<?> getListINITClassesInitialization() {
		return listINITClassesInitialization;
	}
	/**
	 * @param listINITClassesInitialization the listINITClassesInitialization to set
	 */
	public void setListINITClassesInitialization(
			List<?> listINITClassesInitialization) {
		this.listINITClassesInitialization = listINITClassesInitialization;
	}
	/**
	 * @return the listINITClassesFinalization
	 */
	public List<?> getListINITClassesFinalization() {
		return listINITClassesFinalization;
	}
	/**
	 * @param listINITClassesFinalization the listINITClassesFinalization to set
	 */
	public void setListINITClassesFinalization(List<?> listINITClassesFinalization) {
		this.listINITClassesFinalization = listINITClassesFinalization;
	}
	/**
	 * @return the listINITClassesTotal
	 */
	public List<?> getListINITClassesTotal() {
		return listINITClassesTotal;
	}
	/**
	 * @param listINITClassesTotal the listINITClassesTotal to set
	 */
	public void setListINITClassesTotal(List<?> listINITClassesTotal) {
		this.listINITClassesTotal = listINITClassesTotal;
	}
	/**
	 * @return the mapINITCategories
	 */
	public Map<?, ?> getMapINITCategories() {
		return mapINITCategories;
	}
	/**
	 * @param mapINITCategories the mapINITCategories to set
	 */
	public void setMapINITCategories(Map<?, ?> mapINITCategories) {
		this.mapINITCategories = mapINITCategories;
	}
	/**
	 * @return the mapINITCategoriesWavefunction
	 */
	public Map<?, ?> getMapINITCategoriesWavefunction() {
		return mapINITCategoriesWavefunction;
	}
	/**
	 * @param mapINITCategoriesWavefunction the mapINITCategoriesWavefunction to set
	 */
	public void setMapINITCategoriesWavefunction(
			Map<?, ?> mapINITCategoriesWavefunction) {
		this.mapINITCategoriesWavefunction = mapINITCategoriesWavefunction;
	}
	/**
	 * @return the listINITCategoriesWavefunctionData
	 */
	public List<?> getListINITCategoriesWavefunctionData() {
		return listINITCategoriesWavefunctionData;
	}
	/**
	 * @param listINITCategoriesWavefunctionData the listINITCategoriesWavefunctionData to set
	 */
	public void setListINITCategoriesWavefunctionData(
			List<?> listINITCategoriesWavefunctionData) {
		this.listINITCategoriesWavefunctionData = listINITCategoriesWavefunctionData;
	}
	/**
	 * @return the stringINITCategoriesWavefunctionClass
	 */
	public String getStringINITCategoriesWavefunctionClass() {
		return stringINITCategoriesWavefunctionClass;
	}
	/**
	 * @param stringINITCategoriesWavefunctionClass the stringINITCategoriesWavefunctionClass to set
	 */
	public void setStringINITCategoriesWavefunctionClass(
			String stringINITCategoriesWavefunctionClass) {
		this.stringINITCategoriesWavefunctionClass = stringINITCategoriesWavefunctionClass;
	}
	/**
	 * @return the stringINITCategoriesWavefunctionInfo
	 */
	public String getStringINITCategoriesWavefunctionInfo() {
		return stringINITCategoriesWavefunctionInfo;
	}
	/**
	 * @param stringINITCategoriesWavefunctionInfo the stringINITCategoriesWavefunctionInfo to set
	 */
	public void setStringINITCategoriesWavefunctionInfo(
			String stringINITCategoriesWavefunctionInfo) {
		this.stringINITCategoriesWavefunctionInfo = stringINITCategoriesWavefunctionInfo;
	}
	/**
	 * @return the mapINITCategoriesApplyLocPotKin
	 */
	public Map<?, ?> getMapINITCategoriesApplyLocPotKin() {
		return mapINITCategoriesApplyLocPotKin;
	}
	/**
	 * @param mapINITCategoriesApplyLocPotKin the mapINITCategoriesApplyLocPotKin to set
	 */
	public void setMapINITCategoriesApplyLocPotKin(
			Map<?, ?> mapINITCategoriesApplyLocPotKin) {
		this.mapINITCategoriesApplyLocPotKin = mapINITCategoriesApplyLocPotKin;
	}
	/**
	 * @return the listINITCategoriesApplyLocPotKinData
	 */
	public List<?> getListINITCategoriesApplyLocPotKinData() {
		return listINITCategoriesApplyLocPotKinData;
	}
	/**
	 * @param listINITCategoriesApplyLocPotKinData the listINITCategoriesApplyLocPotKinData to set
	 */
	public void setListINITCategoriesApplyLocPotKinData(
			List<?> listINITCategoriesApplyLocPotKinData) {
		this.listINITCategoriesApplyLocPotKinData = listINITCategoriesApplyLocPotKinData;
	}
	/**
	 * @return the stringINITCategoriesApplyLocPotKinClass
	 */
	public String getStringINITCategoriesApplyLocPotKinClass() {
		return stringINITCategoriesApplyLocPotKinClass;
	}
	/**
	 * @param stringINITCategoriesApplyLocPotKinClass the stringINITCategoriesApplyLocPotKinClass to set
	 */
	public void setStringINITCategoriesApplyLocPotKinClass(
			String stringINITCategoriesApplyLocPotKinClass) {
		this.stringINITCategoriesApplyLocPotKinClass = stringINITCategoriesApplyLocPotKinClass;
	}
	/**
	 * @return the stringINITCategoriesApplyLocPotKinInfo
	 */
	public String getStringINITCategoriesApplyLocPotKinInfo() {
		return stringINITCategoriesApplyLocPotKinInfo;
	}
	/**
	 * @param stringINITCategoriesApplyLocPotKinInfo the stringINITCategoriesApplyLocPotKinInfo to set
	 */
	public void setStringINITCategoriesApplyLocPotKinInfo(
			String stringINITCategoriesApplyLocPotKinInfo) {
		this.stringINITCategoriesApplyLocPotKinInfo = stringINITCategoriesApplyLocPotKinInfo;
	}
	/**
	 * @return the mapINITCategoriesApplyProj
	 */
	public Map<?, ?> getMapINITCategoriesApplyProj() {
		return mapINITCategoriesApplyProj;
	}
	/**
	 * @param mapINITCategoriesApplyProj the mapINITCategoriesApplyProj to set
	 */
	public void setMapINITCategoriesApplyProj(Map<?, ?> mapINITCategoriesApplyProj) {
		this.mapINITCategoriesApplyProj = mapINITCategoriesApplyProj;
	}
	/**
	 * @return the listINITCategoriesApplyProjData
	 */
	public List<?> getListINITCategoriesApplyProjData() {
		return listINITCategoriesApplyProjData;
	}
	/**
	 * @param listINITCategoriesApplyProjData the listINITCategoriesApplyProjData to set
	 */
	public void setListINITCategoriesApplyProjData(
			List<?> listINITCategoriesApplyProjData) {
		this.listINITCategoriesApplyProjData = listINITCategoriesApplyProjData;
	}
	/**
	 * @return the stringINITCategoriesApplyProjClass
	 */
	public String getStringINITCategoriesApplyProjClass() {
		return stringINITCategoriesApplyProjClass;
	}
	/**
	 * @param stringINITCategoriesApplyProjClass the stringINITCategoriesApplyProjClass to set
	 */
	public void setStringINITCategoriesApplyProjClass(
			String stringINITCategoriesApplyProjClass) {
		this.stringINITCategoriesApplyProjClass = stringINITCategoriesApplyProjClass;
	}
	/**
	 * @return the stringINITCategoriesApplyProjInfo
	 */
	public String getStringINITCategoriesApplyProjInfo() {
		return stringINITCategoriesApplyProjInfo;
	}
	/**
	 * @param stringINITCategoriesApplyProjInfo the stringINITCategoriesApplyProjInfo to set
	 */
	public void setStringINITCategoriesApplyProjInfo(
			String stringINITCategoriesApplyProjInfo) {
		this.stringINITCategoriesApplyProjInfo = stringINITCategoriesApplyProjInfo;
	}
	/**
	 * @return the mapINITCategoriesCrtLocPot
	 */
	public Map<?, ?> getMapINITCategoriesCrtLocPot() {
		return mapINITCategoriesCrtLocPot;
	}
	/**
	 * @param mapINITCategoriesCrtLocPot the mapINITCategoriesCrtLocPot to set
	 */
	public void setMapINITCategoriesCrtLocPot(Map<?, ?> mapINITCategoriesCrtLocPot) {
		this.mapINITCategoriesCrtLocPot = mapINITCategoriesCrtLocPot;
	}
	/**
	 * @return the listINITCategoriesCrtLocPotData
	 */
	public List<?> getListINITCategoriesCrtLocPotData() {
		return listINITCategoriesCrtLocPotData;
	}
	/**
	 * @param listINITCategoriesCrtLocPotData the listINITCategoriesCrtLocPotData to set
	 */
	public void setListINITCategoriesCrtLocPotData(
			List<?> listINITCategoriesCrtLocPotData) {
		this.listINITCategoriesCrtLocPotData = listINITCategoriesCrtLocPotData;
	}
	/**
	 * @return the stringINITCategoriesCrtLocPotClass
	 */
	public String getStringINITCategoriesCrtLocPotClass() {
		return stringINITCategoriesCrtLocPotClass;
	}
	/**
	 * @param stringINITCategoriesCrtLocPotClass the stringINITCategoriesCrtLocPotClass to set
	 */
	public void setStringINITCategoriesCrtLocPotClass(
			String stringINITCategoriesCrtLocPotClass) {
		this.stringINITCategoriesCrtLocPotClass = stringINITCategoriesCrtLocPotClass;
	}
	/**
	 * @return the stringINITCategoriesCrtLocPotInfo
	 */
	public String getStringINITCategoriesCrtLocPotInfo() {
		return stringINITCategoriesCrtLocPotInfo;
	}
	/**
	 * @param stringINITCategoriesCrtLocPotInfo the stringINITCategoriesCrtLocPotInfo to set
	 */
	public void setStringINITCategoriesCrtLocPotInfo(
			String stringINITCategoriesCrtLocPotInfo) {
		this.stringINITCategoriesCrtLocPotInfo = stringINITCategoriesCrtLocPotInfo;
	}
	/**
	 * @return the mapINITCategoriesRho_comput
	 */
	public Map<?, ?> getMapINITCategoriesRho_comput() {
		return mapINITCategoriesRho_comput;
	}
	/**
	 * @param mapINITCategoriesRho_comput the mapINITCategoriesRho_comput to set
	 */
	public void setMapINITCategoriesRho_comput(Map<?, ?> mapINITCategoriesRho_comput) {
		this.mapINITCategoriesRho_comput = mapINITCategoriesRho_comput;
	}
	/**
	 * @return the listINITCategoriesRho_computData
	 */
	public List<?> getListINITCategoriesRho_computData() {
		return listINITCategoriesRho_computData;
	}
	/**
	 * @param listINITCategoriesRho_computData the listINITCategoriesRho_computData to set
	 */
	public void setListINITCategoriesRho_computData(
			List<?> listINITCategoriesRho_computData) {
		this.listINITCategoriesRho_computData = listINITCategoriesRho_computData;
	}
	/**
	 * @return the stringINITCategoriesRho_computClass
	 */
	public String getStringINITCategoriesRho_computClass() {
		return stringINITCategoriesRho_computClass;
	}
	/**
	 * @param stringINITCategoriesRho_computClass the stringINITCategoriesRho_computClass to set
	 */
	public void setStringINITCategoriesRho_computClass(
			String stringINITCategoriesRho_computClass) {
		this.stringINITCategoriesRho_computClass = stringINITCategoriesRho_computClass;
	}
	/**
	 * @return the stringINITCategoriesRho_computInfo
	 */
	public String getStringINITCategoriesRho_computInfo() {
		return stringINITCategoriesRho_computInfo;
	}
	/**
	 * @param stringINITCategoriesRho_computInfo the stringINITCategoriesRho_computInfo to set
	 */
	public void setStringINITCategoriesRho_computInfo(
			String stringINITCategoriesRho_computInfo) {
		this.stringINITCategoriesRho_computInfo = stringINITCategoriesRho_computInfo;
	}
	/**
	 * @return the mapINITCategoriesUn_TransSwitch
	 */
	public Map<?, ?> getMapINITCategoriesUn_TransSwitch() {
		return mapINITCategoriesUn_TransSwitch;
	}
	/**
	 * @param mapINITCategoriesUn_TransSwitch the mapINITCategoriesUn_TransSwitch to set
	 */
	public void setMapINITCategoriesUn_TransSwitch(
			Map<?, ?> mapINITCategoriesUn_TransSwitch) {
		this.mapINITCategoriesUn_TransSwitch = mapINITCategoriesUn_TransSwitch;
	}
	/**
	 * @return the listINITCategoriesUn_TransSwitchData
	 */
	public List<?> getListINITCategoriesUn_TransSwitchData() {
		return listINITCategoriesUn_TransSwitchData;
	}
	/**
	 * @param listINITCategoriesUn_TransSwitchData the listINITCategoriesUn_TransSwitchData to set
	 */
	public void setListINITCategoriesUn_TransSwitchData(
			List<?> listINITCategoriesUn_TransSwitchData) {
		this.listINITCategoriesUn_TransSwitchData = listINITCategoriesUn_TransSwitchData;
	}
	/**
	 * @return the stringINITCategoriesUn_TransSwitchClass
	 */
	public String getStringINITCategoriesUn_TransSwitchClass() {
		return stringINITCategoriesUn_TransSwitchClass;
	}
	/**
	 * @param stringINITCategoriesUn_TransSwitchClass the stringINITCategoriesUn_TransSwitchClass to set
	 */
	public void setStringINITCategoriesUn_TransSwitchClass(
			String stringINITCategoriesUn_TransSwitchClass) {
		this.stringINITCategoriesUn_TransSwitchClass = stringINITCategoriesUn_TransSwitchClass;
	}
	/**
	 * @return the stringINITCategoriesUn_TransSwitchInfo
	 */
	public String getStringINITCategoriesUn_TransSwitchInfo() {
		return stringINITCategoriesUn_TransSwitchInfo;
	}
	/**
	 * @param stringINITCategoriesUn_TransSwitchInfo the stringINITCategoriesUn_TransSwitchInfo to set
	 */
	public void setStringINITCategoriesUn_TransSwitchInfo(
			String stringINITCategoriesUn_TransSwitchInfo) {
		this.stringINITCategoriesUn_TransSwitchInfo = stringINITCategoriesUn_TransSwitchInfo;
	}
	/**
	 * @return the mapINITCategoriesInput_comput
	 */
	public Map<?, ?> getMapINITCategoriesInput_comput() {
		return mapINITCategoriesInput_comput;
	}
	/**
	 * @param mapINITCategoriesInput_comput the mapINITCategoriesInput_comput to set
	 */
	public void setMapINITCategoriesInput_comput(
			Map<?, ?> mapINITCategoriesInput_comput) {
		this.mapINITCategoriesInput_comput = mapINITCategoriesInput_comput;
	}
	/**
	 * @return the listINITCategoriesInput_computData
	 */
	public List<?> getListINITCategoriesInput_computData() {
		return listINITCategoriesInput_computData;
	}
	/**
	 * @param listINITCategoriesInput_computData the listINITCategoriesInput_computData to set
	 */
	public void setListINITCategoriesInput_computData(
			List<?> listINITCategoriesInput_computData) {
		this.listINITCategoriesInput_computData = listINITCategoriesInput_computData;
	}
	/**
	 * @return the stringINITCategoriesInput_computClass
	 */
	public String getStringINITCategoriesInput_computClass() {
		return stringINITCategoriesInput_computClass;
	}
	/**
	 * @param stringINITCategoriesInput_computClass the stringINITCategoriesInput_computClass to set
	 */
	public void setStringINITCategoriesInput_computClass(
			String stringINITCategoriesInput_computClass) {
		this.stringINITCategoriesInput_computClass = stringINITCategoriesInput_computClass;
	}
	/**
	 * @return the stringINITCategoriesInput_computInfo
	 */
	public String getStringINITCategoriesInput_computInfo() {
		return stringINITCategoriesInput_computInfo;
	}
	/**
	 * @param stringINITCategoriesInput_computInfo the stringINITCategoriesInput_computInfo to set
	 */
	public void setStringINITCategoriesInput_computInfo(
			String stringINITCategoriesInput_computInfo) {
		this.stringINITCategoriesInput_computInfo = stringINITCategoriesInput_computInfo;
	}
	/**
	 * @return the mapINITCategoriesPSolv_comput
	 */
	public Map<?, ?> getMapINITCategoriesPSolv_comput() {
		return mapINITCategoriesPSolv_comput;
	}
	/**
	 * @param mapINITCategoriesPSolv_comput the mapINITCategoriesPSolv_comput to set
	 */
	public void setMapINITCategoriesPSolv_comput(
			Map<?, ?> mapINITCategoriesPSolv_comput) {
		this.mapINITCategoriesPSolv_comput = mapINITCategoriesPSolv_comput;
	}
	/**
	 * @return the listINITCategoriesPSolv_computData
	 */
	public List<?> getListINITCategoriesPSolv_computData() {
		return listINITCategoriesPSolv_computData;
	}
	/**
	 * @param listINITCategoriesPSolv_computData the listINITCategoriesPSolv_computData to set
	 */
	public void setListINITCategoriesPSolv_computData(
			List<?> listINITCategoriesPSolv_computData) {
		this.listINITCategoriesPSolv_computData = listINITCategoriesPSolv_computData;
	}
	/**
	 * @return the stringINITCategoriesPSolv_computClass
	 */
	public String getStringINITCategoriesPSolv_computClass() {
		return stringINITCategoriesPSolv_computClass;
	}
	/**
	 * @param stringINITCategoriesPSolv_computClass the stringINITCategoriesPSolv_computClass to set
	 */
	public void setStringINITCategoriesPSolv_computClass(
			String stringINITCategoriesPSolv_computClass) {
		this.stringINITCategoriesPSolv_computClass = stringINITCategoriesPSolv_computClass;
	}
	/**
	 * @return the stringINITCategoriesPSolv_computInfo
	 */
	public String getStringINITCategoriesPSolv_computInfo() {
		return stringINITCategoriesPSolv_computInfo;
	}
	/**
	 * @param stringINITCategoriesPSolv_computInfo the stringINITCategoriesPSolv_computInfo to set
	 */
	public void setStringINITCategoriesPSolv_computInfo(
			String stringINITCategoriesPSolv_computInfo) {
		this.stringINITCategoriesPSolv_computInfo = stringINITCategoriesPSolv_computInfo;
	}
	/**
	 * @return the mapINITCategoriesExchangecorr
	 */
	public Map<?, ?> getMapINITCategoriesExchangecorr() {
		return mapINITCategoriesExchangecorr;
	}
	/**
	 * @param mapINITCategoriesExchangecorr the mapINITCategoriesExchangecorr to set
	 */
	public void setMapINITCategoriesExchangecorr(
			Map<?, ?> mapINITCategoriesExchangecorr) {
		this.mapINITCategoriesExchangecorr = mapINITCategoriesExchangecorr;
	}
	/**
	 * @return the listINITCategoriesExchangecorrData
	 */
	public List<?> getListINITCategoriesExchangecorrData() {
		return listINITCategoriesExchangecorrData;
	}
	/**
	 * @param listINITCategoriesExchangecorrData the listINITCategoriesExchangecorrData to set
	 */
	public void setListINITCategoriesExchangecorrData(
			List<?> listINITCategoriesExchangecorrData) {
		this.listINITCategoriesExchangecorrData = listINITCategoriesExchangecorrData;
	}
	/**
	 * @return the stringINITCategoriesExchangecorrClass
	 */
	public String getStringINITCategoriesExchangecorrClass() {
		return stringINITCategoriesExchangecorrClass;
	}
	/**
	 * @param stringINITCategoriesExchangecorrClass the stringINITCategoriesExchangecorrClass to set
	 */
	public void setStringINITCategoriesExchangecorrClass(
			String stringINITCategoriesExchangecorrClass) {
		this.stringINITCategoriesExchangecorrClass = stringINITCategoriesExchangecorrClass;
	}
	/**
	 * @return the stringINITCategoriesExchangecorrInfo
	 */
	public String getStringINITCategoriesExchangecorrInfo() {
		return stringINITCategoriesExchangecorrInfo;
	}
	/**
	 * @param stringINITCategoriesExchangecorrInfo the stringINITCategoriesExchangecorrInfo to set
	 */
	public void setStringINITCategoriesExchangecorrInfo(
			String stringINITCategoriesExchangecorrInfo) {
		this.stringINITCategoriesExchangecorrInfo = stringINITCategoriesExchangecorrInfo;
	}
	/**
	 * @return the mapINITCategoriesCrtDescriptors
	 */
	public Map<?, ?> getMapINITCategoriesCrtDescriptors() {
		return mapINITCategoriesCrtDescriptors;
	}
	/**
	 * @param mapINITCategoriesCrtDescriptors the mapINITCategoriesCrtDescriptors to set
	 */
	public void setMapINITCategoriesCrtDescriptors(
			Map<?, ?> mapINITCategoriesCrtDescriptors) {
		this.mapINITCategoriesCrtDescriptors = mapINITCategoriesCrtDescriptors;
	}
	/**
	 * @return the listINITCategoriesCrtDescriptorsData
	 */
	public List<?> getListINITCategoriesCrtDescriptorsData() {
		return listINITCategoriesCrtDescriptorsData;
	}
	/**
	 * @param listINITCategoriesCrtDescriptorsData the listINITCategoriesCrtDescriptorsData to set
	 */
	public void setListINITCategoriesCrtDescriptorsData(
			List<?> listINITCategoriesCrtDescriptorsData) {
		this.listINITCategoriesCrtDescriptorsData = listINITCategoriesCrtDescriptorsData;
	}
	/**
	 * @return the stringINITCategoriesCrtDescriptorsClass
	 */
	public String getStringINITCategoriesCrtDescriptorsClass() {
		return stringINITCategoriesCrtDescriptorsClass;
	}
	/**
	 * @param stringINITCategoriesCrtDescriptorsClass the stringINITCategoriesCrtDescriptorsClass to set
	 */
	public void setStringINITCategoriesCrtDescriptorsClass(
			String stringINITCategoriesCrtDescriptorsClass) {
		this.stringINITCategoriesCrtDescriptorsClass = stringINITCategoriesCrtDescriptorsClass;
	}
	/**
	 * @return the stringINITCategoriesCrtDescriptorsInfo
	 */
	public String getStringINITCategoriesCrtDescriptorsInfo() {
		return stringINITCategoriesCrtDescriptorsInfo;
	}
	/**
	 * @param stringINITCategoriesCrtDescriptorsInfo the stringINITCategoriesCrtDescriptorsInfo to set
	 */
	public void setStringINITCategoriesCrtDescriptorsInfo(
			String stringINITCategoriesCrtDescriptorsInfo) {
		this.stringINITCategoriesCrtDescriptorsInfo = stringINITCategoriesCrtDescriptorsInfo;
	}
	/**
	 * @return the mapINITCategoriesPSolvKernel
	 */
	public Map<?, ?> getMapINITCategoriesPSolvKernel() {
		return mapINITCategoriesPSolvKernel;
	}
	/**
	 * @param mapINITCategoriesPSolvKernel the mapINITCategoriesPSolvKernel to set
	 */
	public void setMapINITCategoriesPSolvKernel(
			Map<?, ?> mapINITCategoriesPSolvKernel) {
		this.mapINITCategoriesPSolvKernel = mapINITCategoriesPSolvKernel;
	}
	/**
	 * @return the listINITCategoriesPSolvKernelData
	 */
	public List<?> getListINITCategoriesPSolvKernelData() {
		return listINITCategoriesPSolvKernelData;
	}
	/**
	 * @param listINITCategoriesPSolvKernelData the listINITCategoriesPSolvKernelData to set
	 */
	public void setListINITCategoriesPSolvKernelData(
			List<?> listINITCategoriesPSolvKernelData) {
		this.listINITCategoriesPSolvKernelData = listINITCategoriesPSolvKernelData;
	}
	/**
	 * @return the stringINITCategoriesPSolvKernelClass
	 */
	public String getStringINITCategoriesPSolvKernelClass() {
		return stringINITCategoriesPSolvKernelClass;
	}
	/**
	 * @param stringINITCategoriesPSolvKernelClass the stringINITCategoriesPSolvKernelClass to set
	 */
	public void setStringINITCategoriesPSolvKernelClass(
			String stringINITCategoriesPSolvKernelClass) {
		this.stringINITCategoriesPSolvKernelClass = stringINITCategoriesPSolvKernelClass;
	}
	/**
	 * @return the stringINITCategoriesPSolvKernelInfo
	 */
	public String getStringINITCategoriesPSolvKernelInfo() {
		return stringINITCategoriesPSolvKernelInfo;
	}
	/**
	 * @param stringINITCategoriesPSolvKernelInfo the stringINITCategoriesPSolvKernelInfo to set
	 */
	public void setStringINITCategoriesPSolvKernelInfo(
			String stringINITCategoriesPSolvKernelInfo) {
		this.stringINITCategoriesPSolvKernelInfo = stringINITCategoriesPSolvKernelInfo;
	}
	/**
	 * @return the mapINITCategoriesPot_commun
	 */
	public Map<?, ?> getMapINITCategoriesPot_commun() {
		return mapINITCategoriesPot_commun;
	}
	/**
	 * @param mapINITCategoriesPot_commun the mapINITCategoriesPot_commun to set
	 */
	public void setMapINITCategoriesPot_commun(Map<?, ?> mapINITCategoriesPot_commun) {
		this.mapINITCategoriesPot_commun = mapINITCategoriesPot_commun;
	}
	/**
	 * @return the listINITCategoriesPot_communData
	 */
	public List<?> getListINITCategoriesPot_communData() {
		return listINITCategoriesPot_communData;
	}
	/**
	 * @param listINITCategoriesPot_communData the listINITCategoriesPot_communData to set
	 */
	public void setListINITCategoriesPot_communData(
			List<?> listINITCategoriesPot_communData) {
		this.listINITCategoriesPot_communData = listINITCategoriesPot_communData;
	}
	/**
	 * @return the stringINITCategoriesPot_communClass
	 */
	public String getStringINITCategoriesPot_communClass() {
		return stringINITCategoriesPot_communClass;
	}
	/**
	 * @param stringINITCategoriesPot_communClass the stringINITCategoriesPot_communClass to set
	 */
	public void setStringINITCategoriesPot_communClass(
			String stringINITCategoriesPot_communClass) {
		this.stringINITCategoriesPot_communClass = stringINITCategoriesPot_communClass;
	}
	/**
	 * @return the stringINITCategoriesPot_communInfo
	 */
	public String getStringINITCategoriesPot_communInfo() {
		return stringINITCategoriesPot_communInfo;
	}
	/**
	 * @param stringINITCategoriesPot_communInfo the stringINITCategoriesPot_communInfo to set
	 */
	public void setStringINITCategoriesPot_communInfo(
			String stringINITCategoriesPot_communInfo) {
		this.stringINITCategoriesPot_communInfo = stringINITCategoriesPot_communInfo;
	}
	/**
	 * @return the mapWFN_OPTClasses
	 */
	public Map<?, ?> getMapWFN_OPTClasses() {
		return mapWFN_OPTClasses;
	}
	/**
	 * @param mapWFN_OPTClasses the mapWFN_OPTClasses to set
	 */
	public void setMapWFN_OPTClasses(Map<?, ?> mapWFN_OPTClasses) {
		this.mapWFN_OPTClasses = mapWFN_OPTClasses;
	}
	/**
	 * @return the listWFN_OPTClassesCommunications
	 */
	public List<?> getListWFN_OPTClassesCommunications() {
		return listWFN_OPTClassesCommunications;
	}
	/**
	 * @param listWFN_OPTClassesCommunications the listWFN_OPTClassesCommunications to set
	 */
	public void setListWFN_OPTClassesCommunications(
			List<?> listWFN_OPTClassesCommunications) {
		this.listWFN_OPTClassesCommunications = listWFN_OPTClassesCommunications;
	}
	/**
	 * @return the listWFN_OPTClassesConvolutions
	 */
	public List<?> getListWFN_OPTClassesConvolutions() {
		return listWFN_OPTClassesConvolutions;
	}
	/**
	 * @param listWFN_OPTClassesConvolutions the listWFN_OPTClassesConvolutions to set
	 */
	public void setListWFN_OPTClassesConvolutions(
			List<?> listWFN_OPTClassesConvolutions) {
		this.listWFN_OPTClassesConvolutions = listWFN_OPTClassesConvolutions;
	}
	/**
	 * @return the listWFN_OPTClassesLinearAlgebra
	 */
	public List<?> getListWFN_OPTClassesLinearAlgebra() {
		return listWFN_OPTClassesLinearAlgebra;
	}
	/**
	 * @param listWFN_OPTClassesLinearAlgebra the listWFN_OPTClassesLinearAlgebra to set
	 */
	public void setListWFN_OPTClassesLinearAlgebra(
			List<?> listWFN_OPTClassesLinearAlgebra) {
		this.listWFN_OPTClassesLinearAlgebra = listWFN_OPTClassesLinearAlgebra;
	}
	/**
	 * @return the listWFN_OPTClassesOther
	 */
	public List<?> getListWFN_OPTClassesOther() {
		return listWFN_OPTClassesOther;
	}
	/**
	 * @param listWFN_OPTClassesOther the listWFN_OPTClassesOther to set
	 */
	public void setListWFN_OPTClassesOther(List<?> listWFN_OPTClassesOther) {
		this.listWFN_OPTClassesOther = listWFN_OPTClassesOther;
	}
	/**
	 * @return the listWFN_OPTClassesPotential
	 */
	public List<?> getListWFN_OPTClassesPotential() {
		return listWFN_OPTClassesPotential;
	}
	/**
	 * @param listWFN_OPTClassesPotential the listWFN_OPTClassesPotential to set
	 */
	public void setListWFN_OPTClassesPotential(List<?> listWFN_OPTClassesPotential) {
		this.listWFN_OPTClassesPotential = listWFN_OPTClassesPotential;
	}
	/**
	 * @return the listWFN_OPTClassesInitialization
	 */
	public List<?> getListWFN_OPTClassesInitialization() {
		return listWFN_OPTClassesInitialization;
	}
	/**
	 * @param listWFN_OPTClassesInitialization the listWFN_OPTClassesInitialization to set
	 */
	public void setListWFN_OPTClassesInitialization(
			List<?> listWFN_OPTClassesInitialization) {
		this.listWFN_OPTClassesInitialization = listWFN_OPTClassesInitialization;
	}
	/**
	 * @return the listWFN_OPTClassesFinalization
	 */
	public List<?> getListWFN_OPTClassesFinalization() {
		return listWFN_OPTClassesFinalization;
	}
	/**
	 * @param listWFN_OPTClassesFinalization the listWFN_OPTClassesFinalization to set
	 */
	public void setListWFN_OPTClassesFinalization(
			List<?> listWFN_OPTClassesFinalization) {
		this.listWFN_OPTClassesFinalization = listWFN_OPTClassesFinalization;
	}
	/**
	 * @return the listWFN_OPTClassesTotal
	 */
	public List<?> getListWFN_OPTClassesTotal() {
		return listWFN_OPTClassesTotal;
	}
	/**
	 * @param listWFN_OPTClassesTotal the listWFN_OPTClassesTotal to set
	 */
	public void setListWFN_OPTClassesTotal(List<?> listWFN_OPTClassesTotal) {
		this.listWFN_OPTClassesTotal = listWFN_OPTClassesTotal;
	}
	/**
	 * @return the mapWFN_OPTCategories
	 */
	public Map<?, ?> getMapWFN_OPTCategories() {
		return mapWFN_OPTCategories;
	}
	/**
	 * @param mapWFN_OPTCategories the mapWFN_OPTCategories to set
	 */
	public void setMapWFN_OPTCategories(Map<?, ?> mapWFN_OPTCategories) {
		this.mapWFN_OPTCategories = mapWFN_OPTCategories;
	}
	/**
	 * @return the mapWFN_OPTCategoriesPrecondition
	 */
	public Map<?, ?> getMapWFN_OPTCategoriesPrecondition() {
		return mapWFN_OPTCategoriesPrecondition;
	}
	/**
	 * @param mapWFN_OPTCategoriesPrecondition the mapWFN_OPTCategoriesPrecondition to set
	 */
	public void setMapWFN_OPTCategoriesPrecondition(
			Map<?, ?> mapWFN_OPTCategoriesPrecondition) {
		this.mapWFN_OPTCategoriesPrecondition = mapWFN_OPTCategoriesPrecondition;
	}
	/**
	 * @return the listWFN_OPTCategoriesPreconditionData
	 */
	public List<?> getListWFN_OPTCategoriesPreconditionData() {
		return listWFN_OPTCategoriesPreconditionData;
	}
	/**
	 * @param listWFN_OPTCategoriesPreconditionData the listWFN_OPTCategoriesPreconditionData to set
	 */
	public void setListWFN_OPTCategoriesPreconditionData(
			List<?> listWFN_OPTCategoriesPreconditionData) {
		this.listWFN_OPTCategoriesPreconditionData = listWFN_OPTCategoriesPreconditionData;
	}
	/**
	 * @return the stringWFN_OPTCategoriesPreconditionClass
	 */
	public String getStringWFN_OPTCategoriesPreconditionClass() {
		return stringWFN_OPTCategoriesPreconditionClass;
	}
	/**
	 * @param stringWFN_OPTCategoriesPreconditionClass the stringWFN_OPTCategoriesPreconditionClass to set
	 */
	public void setStringWFN_OPTCategoriesPreconditionClass(
			String stringWFN_OPTCategoriesPreconditionClass) {
		this.stringWFN_OPTCategoriesPreconditionClass = stringWFN_OPTCategoriesPreconditionClass;
	}
	/**
	 * @return the stringWFN_OPTCategoriesPreconditionInfo
	 */
	public String getStringWFN_OPTCategoriesPreconditionInfo() {
		return stringWFN_OPTCategoriesPreconditionInfo;
	}
	/**
	 * @param stringWFN_OPTCategoriesPreconditionInfo the stringWFN_OPTCategoriesPreconditionInfo to set
	 */
	public void setStringWFN_OPTCategoriesPreconditionInfo(
			String stringWFN_OPTCategoriesPreconditionInfo) {
		this.stringWFN_OPTCategoriesPreconditionInfo = stringWFN_OPTCategoriesPreconditionInfo;
	}
	/**
	 * @return the mapWFN_OPTCategoriesApplyLocPotKin
	 */
	public Map<?, ?> getMapWFN_OPTCategoriesApplyLocPotKin() {
		return mapWFN_OPTCategoriesApplyLocPotKin;
	}
	/**
	 * @param mapWFN_OPTCategoriesApplyLocPotKin the mapWFN_OPTCategoriesApplyLocPotKin to set
	 */
	public void setMapWFN_OPTCategoriesApplyLocPotKin(
			Map<?, ?> mapWFN_OPTCategoriesApplyLocPotKin) {
		this.mapWFN_OPTCategoriesApplyLocPotKin = mapWFN_OPTCategoriesApplyLocPotKin;
	}
	/**
	 * @return the listWFN_OPTCategoriesApplyLocPotKinData
	 */
	public List<?> getListWFN_OPTCategoriesApplyLocPotKinData() {
		return listWFN_OPTCategoriesApplyLocPotKinData;
	}
	/**
	 * @param listWFN_OPTCategoriesApplyLocPotKinData the listWFN_OPTCategoriesApplyLocPotKinData to set
	 */
	public void setListWFN_OPTCategoriesApplyLocPotKinData(
			List<?> listWFN_OPTCategoriesApplyLocPotKinData) {
		this.listWFN_OPTCategoriesApplyLocPotKinData = listWFN_OPTCategoriesApplyLocPotKinData;
	}
	/**
	 * @return the stringWFN_OPTCategoriesApplyLocPotKinClass
	 */
	public String getStringWFN_OPTCategoriesApplyLocPotKinClass() {
		return stringWFN_OPTCategoriesApplyLocPotKinClass;
	}
	/**
	 * @param stringWFN_OPTCategoriesApplyLocPotKinClass the stringWFN_OPTCategoriesApplyLocPotKinClass to set
	 */
	public void setStringWFN_OPTCategoriesApplyLocPotKinClass(
			String stringWFN_OPTCategoriesApplyLocPotKinClass) {
		this.stringWFN_OPTCategoriesApplyLocPotKinClass = stringWFN_OPTCategoriesApplyLocPotKinClass;
	}
	/**
	 * @return the stringWFN_OPTCategoriesApplyLocPotKinInfo
	 */
	public String getStringWFN_OPTCategoriesApplyLocPotKinInfo() {
		return stringWFN_OPTCategoriesApplyLocPotKinInfo;
	}
	/**
	 * @param stringWFN_OPTCategoriesApplyLocPotKinInfo the stringWFN_OPTCategoriesApplyLocPotKinInfo to set
	 */
	public void setStringWFN_OPTCategoriesApplyLocPotKinInfo(
			String stringWFN_OPTCategoriesApplyLocPotKinInfo) {
		this.stringWFN_OPTCategoriesApplyLocPotKinInfo = stringWFN_OPTCategoriesApplyLocPotKinInfo;
	}
	/**
	 * @return the mapWFN_OPTCategoriesApplyProj
	 */
	public Map<?, ?> getMapWFN_OPTCategoriesApplyProj() {
		return mapWFN_OPTCategoriesApplyProj;
	}
	/**
	 * @param mapWFN_OPTCategoriesApplyProj the mapWFN_OPTCategoriesApplyProj to set
	 */
	public void setMapWFN_OPTCategoriesApplyProj(
			Map<?, ?> mapWFN_OPTCategoriesApplyProj) {
		this.mapWFN_OPTCategoriesApplyProj = mapWFN_OPTCategoriesApplyProj;
	}
	/**
	 * @return the listWFN_OPTCategoriesApplyProjData
	 */
	public List<?> getListWFN_OPTCategoriesApplyProjData() {
		return listWFN_OPTCategoriesApplyProjData;
	}
	/**
	 * @param listWFN_OPTCategoriesApplyProjData the listWFN_OPTCategoriesApplyProjData to set
	 */
	public void setListWFN_OPTCategoriesApplyProjData(
			List<?> listWFN_OPTCategoriesApplyProjData) {
		this.listWFN_OPTCategoriesApplyProjData = listWFN_OPTCategoriesApplyProjData;
	}
	/**
	 * @return the stringWFN_OPTCategoriesApplyProjClass
	 */
	public String getStringWFN_OPTCategoriesApplyProjClass() {
		return stringWFN_OPTCategoriesApplyProjClass;
	}
	/**
	 * @param stringWFN_OPTCategoriesApplyProjClass the stringWFN_OPTCategoriesApplyProjClass to set
	 */
	public void setStringWFN_OPTCategoriesApplyProjClass(
			String stringWFN_OPTCategoriesApplyProjClass) {
		this.stringWFN_OPTCategoriesApplyProjClass = stringWFN_OPTCategoriesApplyProjClass;
	}
	/**
	 * @return the stringWFN_OPTCategoriesApplyProjInfo
	 */
	public String getStringWFN_OPTCategoriesApplyProjInfo() {
		return stringWFN_OPTCategoriesApplyProjInfo;
	}
	/**
	 * @param stringWFN_OPTCategoriesApplyProjInfo the stringWFN_OPTCategoriesApplyProjInfo to set
	 */
	public void setStringWFN_OPTCategoriesApplyProjInfo(
			String stringWFN_OPTCategoriesApplyProjInfo) {
		this.stringWFN_OPTCategoriesApplyProjInfo = stringWFN_OPTCategoriesApplyProjInfo;
	}
	/**
	 * @return the mapWFN_OPTCategoriesUn_TransSwitch
	 */
	public Map<?, ?> getMapWFN_OPTCategoriesUn_TransSwitch() {
		return mapWFN_OPTCategoriesUn_TransSwitch;
	}
	/**
	 * @param mapWFN_OPTCategoriesUn_TransSwitch the mapWFN_OPTCategoriesUn_TransSwitch to set
	 */
	public void setMapWFN_OPTCategoriesUn_TransSwitch(
			Map<?, ?> mapWFN_OPTCategoriesUn_TransSwitch) {
		this.mapWFN_OPTCategoriesUn_TransSwitch = mapWFN_OPTCategoriesUn_TransSwitch;
	}
	/**
	 * @return the listWFN_OPTCategoriesUn_TransSwitchData
	 */
	public List<?> getListWFN_OPTCategoriesUn_TransSwitchData() {
		return listWFN_OPTCategoriesUn_TransSwitchData;
	}
	/**
	 * @param listWFN_OPTCategoriesUn_TransSwitchData the listWFN_OPTCategoriesUn_TransSwitchData to set
	 */
	public void setListWFN_OPTCategoriesUn_TransSwitchData(
			List<?> listWFN_OPTCategoriesUn_TransSwitchData) {
		this.listWFN_OPTCategoriesUn_TransSwitchData = listWFN_OPTCategoriesUn_TransSwitchData;
	}
	/**
	 * @return the stringWFN_OPTCategoriesUn_TransSwitchClass
	 */
	public String getStringWFN_OPTCategoriesUn_TransSwitchClass() {
		return stringWFN_OPTCategoriesUn_TransSwitchClass;
	}
	/**
	 * @param stringWFN_OPTCategoriesUn_TransSwitchClass the stringWFN_OPTCategoriesUn_TransSwitchClass to set
	 */
	public void setStringWFN_OPTCategoriesUn_TransSwitchClass(
			String stringWFN_OPTCategoriesUn_TransSwitchClass) {
		this.stringWFN_OPTCategoriesUn_TransSwitchClass = stringWFN_OPTCategoriesUn_TransSwitchClass;
	}
	/**
	 * @return the stringWFN_OPTCategoriesUn_TransSwitchInfo
	 */
	public String getStringWFN_OPTCategoriesUn_TransSwitchInfo() {
		return stringWFN_OPTCategoriesUn_TransSwitchInfo;
	}
	/**
	 * @param stringWFN_OPTCategoriesUn_TransSwitchInfo the stringWFN_OPTCategoriesUn_TransSwitchInfo to set
	 */
	public void setStringWFN_OPTCategoriesUn_TransSwitchInfo(
			String stringWFN_OPTCategoriesUn_TransSwitchInfo) {
		this.stringWFN_OPTCategoriesUn_TransSwitchInfo = stringWFN_OPTCategoriesUn_TransSwitchInfo;
	}
	/**
	 * @return the mapWFN_OPTCategoriesLagrM_comput
	 */
	public Map<?, ?> getMapWFN_OPTCategoriesLagrM_comput() {
		return mapWFN_OPTCategoriesLagrM_comput;
	}
	/**
	 * @param mapWFN_OPTCategoriesLagrM_comput the mapWFN_OPTCategoriesLagrM_comput to set
	 */
	public void setMapWFN_OPTCategoriesLagrM_comput(
			Map<?, ?> mapWFN_OPTCategoriesLagrM_comput) {
		this.mapWFN_OPTCategoriesLagrM_comput = mapWFN_OPTCategoriesLagrM_comput;
	}
	/**
	 * @return the listWFN_OPTCategoriesLagrM_computData
	 */
	public List<?> getListWFN_OPTCategoriesLagrM_computData() {
		return listWFN_OPTCategoriesLagrM_computData;
	}
	/**
	 * @param listWFN_OPTCategoriesLagrM_computData the listWFN_OPTCategoriesLagrM_computData to set
	 */
	public void setListWFN_OPTCategoriesLagrM_computData(
			List<?> listWFN_OPTCategoriesLagrM_computData) {
		this.listWFN_OPTCategoriesLagrM_computData = listWFN_OPTCategoriesLagrM_computData;
	}
	/**
	 * @return the stringWFN_OPTCategoriesLagrM_computClass
	 */
	public String getStringWFN_OPTCategoriesLagrM_computClass() {
		return stringWFN_OPTCategoriesLagrM_computClass;
	}
	/**
	 * @param stringWFN_OPTCategoriesLagrM_computClass the stringWFN_OPTCategoriesLagrM_computClass to set
	 */
	public void setStringWFN_OPTCategoriesLagrM_computClass(
			String stringWFN_OPTCategoriesLagrM_computClass) {
		this.stringWFN_OPTCategoriesLagrM_computClass = stringWFN_OPTCategoriesLagrM_computClass;
	}
	/**
	 * @return the stringWFN_OPTCategoriesLagrM_computInfo
	 */
	public String getStringWFN_OPTCategoriesLagrM_computInfo() {
		return stringWFN_OPTCategoriesLagrM_computInfo;
	}
	/**
	 * @param stringWFN_OPTCategoriesLagrM_computInfo the stringWFN_OPTCategoriesLagrM_computInfo to set
	 */
	public void setStringWFN_OPTCategoriesLagrM_computInfo(
			String stringWFN_OPTCategoriesLagrM_computInfo) {
		this.stringWFN_OPTCategoriesLagrM_computInfo = stringWFN_OPTCategoriesLagrM_computInfo;
	}
	/**
	 * @return the mapWFN_OPTCategoriesRho_comput
	 */
	public Map<?, ?> getMapWFN_OPTCategoriesRho_comput() {
		return mapWFN_OPTCategoriesRho_comput;
	}
	/**
	 * @param mapWFN_OPTCategoriesRho_comput the mapWFN_OPTCategoriesRho_comput to set
	 */
	public void setMapWFN_OPTCategoriesRho_comput(
			Map<?, ?> mapWFN_OPTCategoriesRho_comput) {
		this.mapWFN_OPTCategoriesRho_comput = mapWFN_OPTCategoriesRho_comput;
	}
	/**
	 * @return the listWFN_OPTCategoriesRho_computData
	 */
	public List<?> getListWFN_OPTCategoriesRho_computData() {
		return listWFN_OPTCategoriesRho_computData;
	}
	/**
	 * @param listWFN_OPTCategoriesRho_computData the listWFN_OPTCategoriesRho_computData to set
	 */
	public void setListWFN_OPTCategoriesRho_computData(
			List<?> listWFN_OPTCategoriesRho_computData) {
		this.listWFN_OPTCategoriesRho_computData = listWFN_OPTCategoriesRho_computData;
	}
	/**
	 * @return the stringWFN_OPTCategoriesRho_computClass
	 */
	public String getStringWFN_OPTCategoriesRho_computClass() {
		return stringWFN_OPTCategoriesRho_computClass;
	}
	/**
	 * @param stringWFN_OPTCategoriesRho_computClass the stringWFN_OPTCategoriesRho_computClass to set
	 */
	public void setStringWFN_OPTCategoriesRho_computClass(
			String stringWFN_OPTCategoriesRho_computClass) {
		this.stringWFN_OPTCategoriesRho_computClass = stringWFN_OPTCategoriesRho_computClass;
	}
	/**
	 * @return the stringWFN_OPTCategoriesRho_computInfo
	 */
	public String getStringWFN_OPTCategoriesRho_computInfo() {
		return stringWFN_OPTCategoriesRho_computInfo;
	}
	/**
	 * @param stringWFN_OPTCategoriesRho_computInfo the stringWFN_OPTCategoriesRho_computInfo to set
	 */
	public void setStringWFN_OPTCategoriesRho_computInfo(
			String stringWFN_OPTCategoriesRho_computInfo) {
		this.stringWFN_OPTCategoriesRho_computInfo = stringWFN_OPTCategoriesRho_computInfo;
	}
	/**
	 * @return the mapWFN_OPTCategoriesChol_comput
	 */
	public Map<?, ?> getMapWFN_OPTCategoriesChol_comput() {
		return mapWFN_OPTCategoriesChol_comput;
	}
	/**
	 * @param mapWFN_OPTCategoriesChol_comput the mapWFN_OPTCategoriesChol_comput to set
	 */
	public void setMapWFN_OPTCategoriesChol_comput(
			Map<?, ?> mapWFN_OPTCategoriesChol_comput) {
		this.mapWFN_OPTCategoriesChol_comput = mapWFN_OPTCategoriesChol_comput;
	}
	/**
	 * @return the listWFN_OPTCategoriesChol_computData
	 */
	public List<?> getListWFN_OPTCategoriesChol_computData() {
		return listWFN_OPTCategoriesChol_computData;
	}
	/**
	 * @param listWFN_OPTCategoriesChol_computData the listWFN_OPTCategoriesChol_computData to set
	 */
	public void setListWFN_OPTCategoriesChol_computData(
			List<?> listWFN_OPTCategoriesChol_computData) {
		this.listWFN_OPTCategoriesChol_computData = listWFN_OPTCategoriesChol_computData;
	}
	/**
	 * @return the stringWFN_OPTCategoriesChol_computClass
	 */
	public String getStringWFN_OPTCategoriesChol_computClass() {
		return stringWFN_OPTCategoriesChol_computClass;
	}
	/**
	 * @param stringWFN_OPTCategoriesChol_computClass the stringWFN_OPTCategoriesChol_computClass to set
	 */
	public void setStringWFN_OPTCategoriesChol_computClass(
			String stringWFN_OPTCategoriesChol_computClass) {
		this.stringWFN_OPTCategoriesChol_computClass = stringWFN_OPTCategoriesChol_computClass;
	}
	/**
	 * @return the stringWFN_OPTCategoriesChol_computInfo
	 */
	public String getStringWFN_OPTCategoriesChol_computInfo() {
		return stringWFN_OPTCategoriesChol_computInfo;
	}
	/**
	 * @param stringWFN_OPTCategoriesChol_computInfo the stringWFN_OPTCategoriesChol_computInfo to set
	 */
	public void setStringWFN_OPTCategoriesChol_computInfo(
			String stringWFN_OPTCategoriesChol_computInfo) {
		this.stringWFN_OPTCategoriesChol_computInfo = stringWFN_OPTCategoriesChol_computInfo;
	}
	/**
	 * @return the mapWFN_OPTCategoriesDiis
	 */
	public Map<?, ?> getMapWFN_OPTCategoriesDiis() {
		return mapWFN_OPTCategoriesDiis;
	}
	/**
	 * @param mapWFN_OPTCategoriesDiis the mapWFN_OPTCategoriesDiis to set
	 */
	public void setMapWFN_OPTCategoriesDiis(Map<?, ?> mapWFN_OPTCategoriesDiis) {
		this.mapWFN_OPTCategoriesDiis = mapWFN_OPTCategoriesDiis;
	}
	/**
	 * @return the listWFN_OPTCategoriesDiisData
	 */
	public List<?> getListWFN_OPTCategoriesDiisData() {
		return listWFN_OPTCategoriesDiisData;
	}
	/**
	 * @param listWFN_OPTCategoriesDiisData the listWFN_OPTCategoriesDiisData to set
	 */
	public void setListWFN_OPTCategoriesDiisData(
			List<?> listWFN_OPTCategoriesDiisData) {
		this.listWFN_OPTCategoriesDiisData = listWFN_OPTCategoriesDiisData;
	}
	/**
	 * @return the stringWFN_OPTCategoriesDiisClass
	 */
	public String getStringWFN_OPTCategoriesDiisClass() {
		return stringWFN_OPTCategoriesDiisClass;
	}
	/**
	 * @param stringWFN_OPTCategoriesDiisClass the stringWFN_OPTCategoriesDiisClass to set
	 */
	public void setStringWFN_OPTCategoriesDiisClass(
			String stringWFN_OPTCategoriesDiisClass) {
		this.stringWFN_OPTCategoriesDiisClass = stringWFN_OPTCategoriesDiisClass;
	}
	/**
	 * @return the stringWFN_OPTCategoriesDiisInfo
	 */
	public String getStringWFN_OPTCategoriesDiisInfo() {
		return stringWFN_OPTCategoriesDiisInfo;
	}
	/**
	 * @param stringWFN_OPTCategoriesDiisInfo the stringWFN_OPTCategoriesDiisInfo to set
	 */
	public void setStringWFN_OPTCategoriesDiisInfo(
			String stringWFN_OPTCategoriesDiisInfo) {
		this.stringWFN_OPTCategoriesDiisInfo = stringWFN_OPTCategoriesDiisInfo;
	}
	/**
	 * @return the mapWFN_OPTCategoriesPSolv_comput
	 */
	public Map<?, ?> getMapWFN_OPTCategoriesPSolv_comput() {
		return mapWFN_OPTCategoriesPSolv_comput;
	}
	/**
	 * @param mapWFN_OPTCategoriesPSolv_comput the mapWFN_OPTCategoriesPSolv_comput to set
	 */
	public void setMapWFN_OPTCategoriesPSolv_comput(
			Map<?, ?> mapWFN_OPTCategoriesPSolv_comput) {
		this.mapWFN_OPTCategoriesPSolv_comput = mapWFN_OPTCategoriesPSolv_comput;
	}
	/**
	 * @return the listWFN_OPTCategoriesPSolv_computData
	 */
	public List<?> getListWFN_OPTCategoriesPSolv_computData() {
		return listWFN_OPTCategoriesPSolv_computData;
	}
	/**
	 * @param listWFN_OPTCategoriesPSolv_computData the listWFN_OPTCategoriesPSolv_computData to set
	 */
	public void setListWFN_OPTCategoriesPSolv_computData(
			List<?> listWFN_OPTCategoriesPSolv_computData) {
		this.listWFN_OPTCategoriesPSolv_computData = listWFN_OPTCategoriesPSolv_computData;
	}
	/**
	 * @return the stringWFN_OPTCategoriesPSolv_computClass
	 */
	public String getStringWFN_OPTCategoriesPSolv_computClass() {
		return stringWFN_OPTCategoriesPSolv_computClass;
	}
	/**
	 * @param stringWFN_OPTCategoriesPSolv_computClass the stringWFN_OPTCategoriesPSolv_computClass to set
	 */
	public void setStringWFN_OPTCategoriesPSolv_computClass(
			String stringWFN_OPTCategoriesPSolv_computClass) {
		this.stringWFN_OPTCategoriesPSolv_computClass = stringWFN_OPTCategoriesPSolv_computClass;
	}
	/**
	 * @return the stringWFN_OPTCategoriesPSolv_computInfo
	 */
	public String getStringWFN_OPTCategoriesPSolv_computInfo() {
		return stringWFN_OPTCategoriesPSolv_computInfo;
	}
	/**
	 * @param stringWFN_OPTCategoriesPSolv_computInfo the stringWFN_OPTCategoriesPSolv_computInfo to set
	 */
	public void setStringWFN_OPTCategoriesPSolv_computInfo(
			String stringWFN_OPTCategoriesPSolv_computInfo) {
		this.stringWFN_OPTCategoriesPSolv_computInfo = stringWFN_OPTCategoriesPSolv_computInfo;
	}
	/**
	 * @return the mapWFN_OPTCategoriesExchangecorr
	 */
	public Map<?, ?> getMapWFN_OPTCategoriesExchangecorr() {
		return mapWFN_OPTCategoriesExchangecorr;
	}
	/**
	 * @param mapWFN_OPTCategoriesExchangecorr the mapWFN_OPTCategoriesExchangecorr to set
	 */
	public void setMapWFN_OPTCategoriesExchangecorr(
			Map<?, ?> mapWFN_OPTCategoriesExchangecorr) {
		this.mapWFN_OPTCategoriesExchangecorr = mapWFN_OPTCategoriesExchangecorr;
	}
	/**
	 * @return the listWFN_OPTCategoriesExchangecorrData
	 */
	public List<?> getListWFN_OPTCategoriesExchangecorrData() {
		return listWFN_OPTCategoriesExchangecorrData;
	}
	/**
	 * @param listWFN_OPTCategoriesExchangecorrData the listWFN_OPTCategoriesExchangecorrData to set
	 */
	public void setListWFN_OPTCategoriesExchangecorrData(
			List<?> listWFN_OPTCategoriesExchangecorrData) {
		this.listWFN_OPTCategoriesExchangecorrData = listWFN_OPTCategoriesExchangecorrData;
	}
	/**
	 * @return the stringWFN_OPTCategoriesExchangecorrClass
	 */
	public String getStringWFN_OPTCategoriesExchangecorrClass() {
		return stringWFN_OPTCategoriesExchangecorrClass;
	}
	/**
	 * @param stringWFN_OPTCategoriesExchangecorrClass the stringWFN_OPTCategoriesExchangecorrClass to set
	 */
	public void setStringWFN_OPTCategoriesExchangecorrClass(
			String stringWFN_OPTCategoriesExchangecorrClass) {
		this.stringWFN_OPTCategoriesExchangecorrClass = stringWFN_OPTCategoriesExchangecorrClass;
	}
	/**
	 * @return the stringWFN_OPTCategoriesExchangecorrInfo
	 */
	public String getStringWFN_OPTCategoriesExchangecorrInfo() {
		return stringWFN_OPTCategoriesExchangecorrInfo;
	}
	/**
	 * @param stringWFN_OPTCategoriesExchangecorrInfo the stringWFN_OPTCategoriesExchangecorrInfo to set
	 */
	public void setStringWFN_OPTCategoriesExchangecorrInfo(
			String stringWFN_OPTCategoriesExchangecorrInfo) {
		this.stringWFN_OPTCategoriesExchangecorrInfo = stringWFN_OPTCategoriesExchangecorrInfo;
	}
	/**
	 * @return the mapWFN_OPTCategoriesPot_commun
	 */
	public Map<?, ?> getMapWFN_OPTCategoriesPot_commun() {
		return mapWFN_OPTCategoriesPot_commun;
	}
	/**
	 * @param mapWFN_OPTCategoriesPot_commun the mapWFN_OPTCategoriesPot_commun to set
	 */
	public void setMapWFN_OPTCategoriesPot_commun(
			Map<?, ?> mapWFN_OPTCategoriesPot_commun) {
		this.mapWFN_OPTCategoriesPot_commun = mapWFN_OPTCategoriesPot_commun;
	}
	/**
	 * @return the listWFN_OPTCategoriesPot_communData
	 */
	public List<?> getListWFN_OPTCategoriesPot_communData() {
		return listWFN_OPTCategoriesPot_communData;
	}
	/**
	 * @param listWFN_OPTCategoriesPot_communData the listWFN_OPTCategoriesPot_communData to set
	 */
	public void setListWFN_OPTCategoriesPot_communData(
			List<?> listWFN_OPTCategoriesPot_communData) {
		this.listWFN_OPTCategoriesPot_communData = listWFN_OPTCategoriesPot_communData;
	}
	/**
	 * @return the stringWFN_OPTCategoriesPot_communClass
	 */
	public String getStringWFN_OPTCategoriesPot_communClass() {
		return stringWFN_OPTCategoriesPot_communClass;
	}
	/**
	 * @param stringWFN_OPTCategoriesPot_communClass the stringWFN_OPTCategoriesPot_communClass to set
	 */
	public void setStringWFN_OPTCategoriesPot_communClass(
			String stringWFN_OPTCategoriesPot_communClass) {
		this.stringWFN_OPTCategoriesPot_communClass = stringWFN_OPTCategoriesPot_communClass;
	}
	/**
	 * @return the stringWFN_OPTCategoriesPot_communInfo
	 */
	public String getStringWFN_OPTCategoriesPot_communInfo() {
		return stringWFN_OPTCategoriesPot_communInfo;
	}
	/**
	 * @param stringWFN_OPTCategoriesPot_communInfo the stringWFN_OPTCategoriesPot_communInfo to set
	 */
	public void setStringWFN_OPTCategoriesPot_communInfo(
			String stringWFN_OPTCategoriesPot_communInfo) {
		this.stringWFN_OPTCategoriesPot_communInfo = stringWFN_OPTCategoriesPot_communInfo;
	}
	/**
	 * @return the mapLASTClasses
	 */
	public Map<?, ?> getMapLASTClasses() {
		return mapLASTClasses;
	}
	/**
	 * @param mapLASTClasses the mapLASTClasses to set
	 */
	public void setMapLASTClasses(Map<?, ?> mapLASTClasses) {
		this.mapLASTClasses = mapLASTClasses;
	}
	/**
	 * @return the listLASTClassesCommunications
	 */
	public List<?> getListLASTClassesCommunications() {
		return listLASTClassesCommunications;
	}
	/**
	 * @param listLASTClassesCommunications the listLASTClassesCommunications to set
	 */
	public void setListLASTClassesCommunications(
			List<?> listLASTClassesCommunications) {
		this.listLASTClassesCommunications = listLASTClassesCommunications;
	}
	/**
	 * @return the listLASTClassesConvolutions
	 */
	public List<?> getListLASTClassesConvolutions() {
		return listLASTClassesConvolutions;
	}
	/**
	 * @param listLASTClassesConvolutions the listLASTClassesConvolutions to set
	 */
	public void setListLASTClassesConvolutions(List<?> listLASTClassesConvolutions) {
		this.listLASTClassesConvolutions = listLASTClassesConvolutions;
	}
	/**
	 * @return the listLASTClassesLinearAlgebra
	 */
	public List<?> getListLASTClassesLinearAlgebra() {
		return listLASTClassesLinearAlgebra;
	}
	/**
	 * @param listLASTClassesLinearAlgebra the listLASTClassesLinearAlgebra to set
	 */
	public void setListLASTClassesLinearAlgebra(List<?> listLASTClassesLinearAlgebra) {
		this.listLASTClassesLinearAlgebra = listLASTClassesLinearAlgebra;
	}
	/**
	 * @return the listLASTClassesOther
	 */
	public List<?> getListLASTClassesOther() {
		return listLASTClassesOther;
	}
	/**
	 * @param listLASTClassesOther the listLASTClassesOther to set
	 */
	public void setListLASTClassesOther(List<?> listLASTClassesOther) {
		this.listLASTClassesOther = listLASTClassesOther;
	}
	/**
	 * @return the listLASTClassesPotential
	 */
	public List<?> getListLASTClassesPotential() {
		return listLASTClassesPotential;
	}
	/**
	 * @param listLASTClassesPotential the listLASTClassesPotential to set
	 */
	public void setListLASTClassesPotential(List<?> listLASTClassesPotential) {
		this.listLASTClassesPotential = listLASTClassesPotential;
	}
	/**
	 * @return the listLASTClassesInitialization
	 */
	public List<?> getListLASTClassesInitialization() {
		return listLASTClassesInitialization;
	}
	/**
	 * @param listLASTClassesInitialization the listLASTClassesInitialization to set
	 */
	public void setListLASTClassesInitialization(
			List<?> listLASTClassesInitialization) {
		this.listLASTClassesInitialization = listLASTClassesInitialization;
	}
	/**
	 * @return the listLASTClassesFinalization
	 */
	public List<?> getListLASTClassesFinalization() {
		return listLASTClassesFinalization;
	}
	/**
	 * @param listLASTClassesFinalization the listLASTClassesFinalization to set
	 */
	public void setListLASTClassesFinalization(List<?> listLASTClassesFinalization) {
		this.listLASTClassesFinalization = listLASTClassesFinalization;
	}
	/**
	 * @return the listLASTClassesTotal
	 */
	public List<?> getListLASTClassesTotal() {
		return listLASTClassesTotal;
	}
	/**
	 * @param listLASTClassesTotal the listLASTClassesTotal to set
	 */
	public void setListLASTClassesTotal(List<?> listLASTClassesTotal) {
		this.listLASTClassesTotal = listLASTClassesTotal;
	}
	/**
	 * @return the mapLASTCategories
	 */
	public Map<?, ?> getMapLASTCategories() {
		return mapLASTCategories;
	}
	/**
	 * @param mapLASTCategories the mapLASTCategories to set
	 */
	public void setMapLASTCategories(Map<?, ?> mapLASTCategories) {
		this.mapLASTCategories = mapLASTCategories;
	}
	/**
	 * @return the mapLASTCategoriesForces
	 */
	public Map<?, ?> getMapLASTCategoriesForces() {
		return mapLASTCategoriesForces;
	}
	/**
	 * @param mapLASTCategoriesForces the mapLASTCategoriesForces to set
	 */
	public void setMapLASTCategoriesForces(Map<?, ?> mapLASTCategoriesForces) {
		this.mapLASTCategoriesForces = mapLASTCategoriesForces;
	}
	/**
	 * @return the listLASTCategoriesForcesData
	 */
	public List<?> getListLASTCategoriesForcesData() {
		return listLASTCategoriesForcesData;
	}
	/**
	 * @param listLASTCategoriesForcesData the listLASTCategoriesForcesData to set
	 */
	public void setListLASTCategoriesForcesData(List<?> listLASTCategoriesForcesData) {
		this.listLASTCategoriesForcesData = listLASTCategoriesForcesData;
	}
	/**
	 * @return the stringLASTCategoriesForcesClass
	 */
	public String getStringLASTCategoriesForcesClass() {
		return stringLASTCategoriesForcesClass;
	}
	/**
	 * @param stringLASTCategoriesForcesClass the stringLASTCategoriesForcesClass to set
	 */
	public void setStringLASTCategoriesForcesClass(
			String stringLASTCategoriesForcesClass) {
		this.stringLASTCategoriesForcesClass = stringLASTCategoriesForcesClass;
	}
	/**
	 * @return the stringLASTCategoriesForcesInfo
	 */
	public String getStringLASTCategoriesForcesInfo() {
		return stringLASTCategoriesForcesInfo;
	}
	/**
	 * @param stringLASTCategoriesForcesInfo the stringLASTCategoriesForcesInfo to set
	 */
	public void setStringLASTCategoriesForcesInfo(
			String stringLASTCategoriesForcesInfo) {
		this.stringLASTCategoriesForcesInfo = stringLASTCategoriesForcesInfo;
	}
	/**
	 * @return the mapLASTCategoriesRho_comput
	 */
	public Map<?, ?> getMapLASTCategoriesRho_comput() {
		return mapLASTCategoriesRho_comput;
	}
	/**
	 * @param mapLASTCategoriesRho_comput the mapLASTCategoriesRho_comput to set
	 */
	public void setMapLASTCategoriesRho_comput(Map<?, ?> mapLASTCategoriesRho_comput) {
		this.mapLASTCategoriesRho_comput = mapLASTCategoriesRho_comput;
	}
	/**
	 * @return the listLASTCategoriesRho_computData
	 */
	public List<?> getListLASTCategoriesRho_computData() {
		return listLASTCategoriesRho_computData;
	}
	/**
	 * @param listLASTCategoriesRho_computData the listLASTCategoriesRho_computData to set
	 */
	public void setListLASTCategoriesRho_computData(
			List<?> listLASTCategoriesRho_computData) {
		this.listLASTCategoriesRho_computData = listLASTCategoriesRho_computData;
	}
	/**
	 * @return the stringLASTCategoriesRho_computClass
	 */
	public String getStringLASTCategoriesRho_computClass() {
		return stringLASTCategoriesRho_computClass;
	}
	/**
	 * @param stringLASTCategoriesRho_computClass the stringLASTCategoriesRho_computClass to set
	 */
	public void setStringLASTCategoriesRho_computClass(
			String stringLASTCategoriesRho_computClass) {
		this.stringLASTCategoriesRho_computClass = stringLASTCategoriesRho_computClass;
	}
	/**
	 * @return the stringLASTCategoriesRho_computInfo
	 */
	public String getStringLASTCategoriesRho_computInfo() {
		return stringLASTCategoriesRho_computInfo;
	}
	/**
	 * @param stringLASTCategoriesRho_computInfo the stringLASTCategoriesRho_computInfo to set
	 */
	public void setStringLASTCategoriesRho_computInfo(
			String stringLASTCategoriesRho_computInfo) {
		this.stringLASTCategoriesRho_computInfo = stringLASTCategoriesRho_computInfo;
	}
	/**
	 * @return the mapLASTCategoriesPSolv_comput
	 */
	public Map<?, ?> getMapLASTCategoriesPSolv_comput() {
		return mapLASTCategoriesPSolv_comput;
	}
	/**
	 * @param mapLASTCategoriesPSolv_comput the mapLASTCategoriesPSolv_comput to set
	 */
	public void setMapLASTCategoriesPSolv_comput(
			Map<?, ?> mapLASTCategoriesPSolv_comput) {
		this.mapLASTCategoriesPSolv_comput = mapLASTCategoriesPSolv_comput;
	}
	/**
	 * @return the listLASTCategoriesPSolv_computData
	 */
	public List<?> getListLASTCategoriesPSolv_computData() {
		return listLASTCategoriesPSolv_computData;
	}
	/**
	 * @param listLASTCategoriesPSolv_computData the listLASTCategoriesPSolv_computData to set
	 */
	public void setListLASTCategoriesPSolv_computData(
			List<?> listLASTCategoriesPSolv_computData) {
		this.listLASTCategoriesPSolv_computData = listLASTCategoriesPSolv_computData;
	}
	/**
	 * @return the stringLASTCategoriesPSolv_computClass
	 */
	public String getStringLASTCategoriesPSolv_computClass() {
		return stringLASTCategoriesPSolv_computClass;
	}
	/**
	 * @param stringLASTCategoriesPSolv_computClass the stringLASTCategoriesPSolv_computClass to set
	 */
	public void setStringLASTCategoriesPSolv_computClass(
			String stringLASTCategoriesPSolv_computClass) {
		this.stringLASTCategoriesPSolv_computClass = stringLASTCategoriesPSolv_computClass;
	}
	/**
	 * @return the stringLASTCategoriesPSolv_computInfo
	 */
	public String getStringLASTCategoriesPSolv_computInfo() {
		return stringLASTCategoriesPSolv_computInfo;
	}
	/**
	 * @param stringLASTCategoriesPSolv_computInfo the stringLASTCategoriesPSolv_computInfo to set
	 */
	public void setStringLASTCategoriesPSolv_computInfo(
			String stringLASTCategoriesPSolv_computInfo) {
		this.stringLASTCategoriesPSolv_computInfo = stringLASTCategoriesPSolv_computInfo;
	}
	/**
	 * @return the listSUMMARYINIT
	 */
	public List<?> getListSUMMARYINIT() {
		return listSUMMARYINIT;
	}
	/**
	 * @param listSUMMARYINIT the listSUMMARYINIT to set
	 */
	public void setListSUMMARYINIT(List<?> listSUMMARYINIT) {
		this.listSUMMARYINIT = listSUMMARYINIT;
	}
	/**
	 * @return the listSUMMARYWFN_OPT
	 */
	public List<?> getListSUMMARYWFN_OPT() {
		return listSUMMARYWFN_OPT;
	}
	/**
	 * @param listSUMMARYWFN_OPT the listSUMMARYWFN_OPT to set
	 */
	public void setListSUMMARYWFN_OPT(List<?> listSUMMARYWFN_OPT) {
		this.listSUMMARYWFN_OPT = listSUMMARYWFN_OPT;
	}
	/**
	 * @return the listSUMMARYLAST
	 */
	public List<?> getListSUMMARYLAST() {
		return listSUMMARYLAST;
	}
	/**
	 * @param listSUMMARYLAST the listSUMMARYLAST to set
	 */
	public void setListSUMMARYLAST(List<?> listSUMMARYLAST) {
		this.listSUMMARYLAST = listSUMMARYLAST;
	}
	/**
	 * @return the listSUMMARYTotal
	 */
	public List<?> getListSUMMARYTotal() {
		return listSUMMARYTotal;
	}
	/**
	 * @param listSUMMARYTotal the listSUMMARYTotal to set
	 */
	public void setListSUMMARYTotal(List<?> listSUMMARYTotal) {
		this.listSUMMARYTotal = listSUMMARYTotal;
	}
	/**
	 * @return the mapSUMMARYCPUParallelism
	 */
	public Map<?, ?> getMapSUMMARYCPUParallelism() {
		return mapSUMMARYCPUParallelism;
	}
	/**
	 * @param mapSUMMARYCPUParallelism the mapSUMMARYCPUParallelism to set
	 */
	public void setMapSUMMARYCPUParallelism(Map<?, ?> mapSUMMARYCPUParallelism) {
		this.mapSUMMARYCPUParallelism = mapSUMMARYCPUParallelism;
	}
	/**
	 * @return the stringSUMMARYCPUParallelismMPIprocs
	 */
	public String getStringSUMMARYCPUParallelismMPIprocs() {
		return stringSUMMARYCPUParallelismMPIprocs;
	}
	/**
	 * @param stringSUMMARYCPUParallelismMPIprocs the stringSUMMARYCPUParallelismMPIprocs to set
	 */
	public void setStringSUMMARYCPUParallelismMPIprocs(
			String stringSUMMARYCPUParallelismMPIprocs) {
		this.stringSUMMARYCPUParallelismMPIprocs = stringSUMMARYCPUParallelismMPIprocs;
	}
	/**
	 * @return the stringSUMMARYCPUParallelismOMPthrds
	 */
	public String getStringSUMMARYCPUParallelismOMPthrds() {
		return stringSUMMARYCPUParallelismOMPthrds;
	}
	/**
	 * @param stringSUMMARYCPUParallelismOMPthrds the stringSUMMARYCPUParallelismOMPthrds to set
	 */
	public void setStringSUMMARYCPUParallelismOMPthrds(
			String stringSUMMARYCPUParallelismOMPthrds) {
		this.stringSUMMARYCPUParallelismOMPthrds = stringSUMMARYCPUParallelismOMPthrds;
	}
	/**--------------------------------------------------------------*/	
	/**--------------------------------------------------------------*/	
	/**--------------------------------------------------------------*/	
	/**
	 * @return the executedModules
	 */
	public HashMap<String, Boolean> getExecutedModules() {
		return executedModules;
	}
	/**
	 * @param executedModules the executedModules to set
	 */
	public void setExecutedModules(HashMap<String, Boolean> executedModules) {
		this.executedModules = executedModules;
	}

}


