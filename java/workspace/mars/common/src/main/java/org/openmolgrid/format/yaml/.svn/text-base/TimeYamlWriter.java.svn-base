/**
 * 
 */
package org.openmolgrid.format.yaml;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * @author Frederic Bonnet
 *
 */
public class TimeYamlWriter {
	//global variables
	private TimeYamlVariables timeyamlparse;
	private String path;
	//local variables
	/**
	 * simple constructor
	 */
	public TimeYamlWriter() {
	}
	public TimeYamlWriter(String path, TimeYamlVariables timeyamlparse) {
		this.path = path;
		this.timeyamlparse = timeyamlparse;		
	}
	/**
	 * 
	 * @param path
	 * @param timeyamlparse
	 * @throws IOException
	 */
	public void writeTimeYaml2File() /**String path, TimeYamlVariables timeyamlparse)*/ throws IOException {

		FileWriter fw = null;
		BufferedWriter bw = null;
		
		File commandFile = new File(path);
		try {
			fw = new FileWriter(commandFile);
			bw = new BufferedWriter(fw);

			bw.append("---------------"+"-in class wri"+"TimeYaml2File-"+"-----------------"+"--------------------");
			bw.newLine();
			bw.append("---------------"+"--------------"+"--------------"+"-----------------"+"--------------------");
			bw.newLine();
			bw.append("name           "+"percentage (%)"+"      Time (s)"+"          Class  "+"         Info       ");
			bw.newLine();
			bw.append("---------------"+"--------------"+"----------INIT"+"-----------------"+"--------------------");
			bw.newLine();
			bw.append(timeyamlparse.getMapINIT()+"");
			bw.newLine();
			bw.append("---------------"+"--------------"+"--------------"+"-----------------"+"--------------------");
			bw.newLine();
			bw.append("Communications     "+timeyamlparse.getListINITClassesCommunications().get(0)       +"            "+timeyamlparse.getListINITClassesCommunications().get(1));
			bw.newLine();
			bw.append("Convolutions       "+timeyamlparse.getListINITClassesConvolutions().get(0)         +"            "+timeyamlparse.getListINITClassesConvolutions().get(1));
			bw.newLine();
			bw.append("Linear Algebra     "+timeyamlparse.getListINITClassesLinearAlgebra().get(0)        +"            "+timeyamlparse.getListINITClassesLinearAlgebra().get(1));
			bw.newLine();
			bw.append("Other              "+timeyamlparse.getListINITClassesOther().get(0)                +"            "+timeyamlparse.getListINITClassesOther().get(1));
			bw.newLine();
			bw.append("Potential          "+timeyamlparse.getListINITClassesPotential().get(0)            +"            "+timeyamlparse.getListINITClassesPotential().get(1));
			bw.newLine();
			bw.append("Initialization     "+timeyamlparse.getListINITClassesInitialization().get(0)       +"            "+timeyamlparse.getListINITClassesInitialization().get(1));
			bw.newLine();
			bw.append("Finalization       "+timeyamlparse.getListINITClassesFinalization().get(0)         +"            "+timeyamlparse.getListINITClassesFinalization().get(1));
			bw.newLine();
			bw.append("---------------"+"--------------"+"--------------"+"-----------------"+"--------------------");
			bw.newLine();
			bw.append("Total              "+timeyamlparse.getListINITClassesTotal().get(0)                +"            "+timeyamlparse.getListINITClassesTotal().get(1));
			bw.newLine();
			bw.append("---------------"+"--------------"+"--------------"+"-----------------"+"--------------------");
			bw.newLine();
			bw.append("wavefunction       "+timeyamlparse.getListINITCategoriesWavefunctionData().get(0)  +"            "+timeyamlparse.getListINITCategoriesWavefunctionData().get(1)  +"     "+timeyamlparse.getStringINITCategoriesWavefunctionClass()  +"      "+timeyamlparse.getStringINITCategoriesWavefunctionInfo()  );
			bw.newLine();
			bw.append("ApplyLocPotKin     "+timeyamlparse.getListINITCategoriesApplyLocPotKinData().get(0)+"            "+timeyamlparse.getListINITCategoriesApplyLocPotKinData().get(1)+"     "+timeyamlparse.getStringINITCategoriesApplyLocPotKinClass()+"      "+timeyamlparse.getStringINITCategoriesApplyLocPotKinInfo());
			bw.newLine();
			bw.append("ApplyProj          "+timeyamlparse.getListINITCategoriesApplyProjData().get(0)     +"            "+timeyamlparse.getListINITCategoriesApplyProjData().get(1)     +"     "+timeyamlparse.getStringINITCategoriesApplyProjClass()     +"      "+timeyamlparse.getStringINITCategoriesApplyProjInfo()     );
			bw.newLine();
			bw.append("CrtLocPot          "+timeyamlparse.getListINITCategoriesCrtLocPotData().get(0)     +"            "+timeyamlparse.getListINITCategoriesCrtLocPotData().get(1)     +"     "+timeyamlparse.getStringINITCategoriesCrtLocPotClass()     +"      "+timeyamlparse.getStringINITCategoriesCrtLocPotInfo()     );
			bw.newLine();
			bw.append("Rho_comput         "+timeyamlparse.getListINITCategoriesRho_computData().get(0)    +"            "+timeyamlparse.getListINITCategoriesRho_computData().get(1)    +"     "+timeyamlparse.getStringINITCategoriesRho_computClass()    +"      "+timeyamlparse.getStringINITCategoriesRho_computInfo()    );
			bw.newLine();
			bw.append("Un-TransSwitch     "+timeyamlparse.getListINITCategoriesUn_TransSwitchData().get(0)+"            "+timeyamlparse.getListINITCategoriesUn_TransSwitchData().get(1)+"     "+timeyamlparse.getStringINITCategoriesUn_TransSwitchClass()+"      "+timeyamlparse.getStringINITCategoriesUn_TransSwitchInfo());
			bw.newLine();
			bw.append("Input_comput       "+timeyamlparse.getListINITCategoriesInput_computData().get(0)  +"            "+timeyamlparse.getListINITCategoriesInput_computData().get(1)  +"     "+timeyamlparse.getStringINITCategoriesInput_computClass()  +"      "+timeyamlparse.getStringINITCategoriesInput_computInfo()  );
			bw.newLine();
			bw.append("PSolv_comput       "+timeyamlparse.getListINITCategoriesPSolv_computData().get(0)  +"            "+timeyamlparse.getListINITCategoriesPSolv_computData().get(1)  +"     "+timeyamlparse.getStringINITCategoriesPSolv_computClass()  +"      "+timeyamlparse.getStringINITCategoriesPSolv_computInfo()  );
			bw.newLine();
			bw.append("Exchangecorr       "+timeyamlparse.getListINITCategoriesExchangecorrData().get(0)  +"            "+timeyamlparse.getListINITCategoriesExchangecorrData().get(1)  +"     "+timeyamlparse.getStringINITCategoriesExchangecorrClass()  +"      "+timeyamlparse.getStringINITCategoriesExchangecorrInfo()  );
			bw.newLine();
			bw.append("CrtDescriptors     "+timeyamlparse.getListINITCategoriesCrtDescriptorsData().get(0)+"            "+timeyamlparse.getListINITCategoriesCrtDescriptorsData().get(1)+"     "+timeyamlparse.getStringINITCategoriesCrtDescriptorsClass()+"      "+timeyamlparse.getStringINITCategoriesCrtDescriptorsInfo());
			bw.newLine();
			bw.append("PSolvKernel        "+timeyamlparse.getListINITCategoriesPSolvKernelData().get(0)   +"            "+timeyamlparse.getListINITCategoriesPSolvKernelData().get(1)   +"     "+timeyamlparse.getStringINITCategoriesPSolvKernelClass()   +"      "+timeyamlparse.getStringINITCategoriesPSolvKernelInfo()   );
			bw.newLine();
			bw.append("Pot_commun         "+timeyamlparse.getListINITCategoriesPot_communData().get(0)    +"            "+timeyamlparse.getListINITCategoriesPot_communData().get(1)    +"     "+timeyamlparse.getStringINITCategoriesPot_communClass()    +"      "+timeyamlparse.getStringINITCategoriesPot_communInfo()    );
			bw.newLine();
			bw.append("---------------"+"--------------"+"---------WFN_O"+"PT---------------"+"--------------------");
			bw.newLine();
			bw.append(timeyamlparse.getMapWFN_OPT()+"");
			bw.newLine();
			bw.append("---------------"+"--------------"+"--------------"+"-----------------"+"--------------------");
			bw.newLine();
			bw.append("Communications     "+timeyamlparse.getListWFN_OPTClassesCommunications().get(0)       +"            "+timeyamlparse.getListWFN_OPTClassesCommunications().get(1));
			bw.newLine();
			bw.append("Convolutions       "+timeyamlparse.getListWFN_OPTClassesConvolutions().get(0)         +"            "+timeyamlparse.getListWFN_OPTClassesConvolutions().get(1));
			bw.newLine();
			bw.append("Linear Algebra     "+timeyamlparse.getListWFN_OPTClassesLinearAlgebra().get(0)        +"            "+timeyamlparse.getListWFN_OPTClassesLinearAlgebra().get(1));
			bw.newLine();
			bw.append("Other              "+timeyamlparse.getListWFN_OPTClassesOther().get(0)                +"            "+timeyamlparse.getListWFN_OPTClassesOther().get(1));
			bw.newLine();
			bw.append("Potential          "+timeyamlparse.getListWFN_OPTClassesPotential().get(0)            +"            "+timeyamlparse.getListWFN_OPTClassesPotential().get(1));
			bw.newLine();
			bw.append("Initialization     "+timeyamlparse.getListWFN_OPTClassesInitialization().get(0)       +"            "+timeyamlparse.getListWFN_OPTClassesInitialization().get(1));
			bw.newLine();
			bw.append("Finalization       "+timeyamlparse.getListWFN_OPTClassesFinalization().get(0)         +"            "+timeyamlparse.getListWFN_OPTClassesFinalization().get(1));
			bw.newLine();
			bw.append("---------------"+"--------------"+"--------------"+"-----------------"+"--------------------");
			bw.newLine();
			bw.append("Total              "+timeyamlparse.getListWFN_OPTClassesTotal().get(0)                +"            "+timeyamlparse.getListWFN_OPTClassesTotal().get(1));
			bw.newLine();
			bw.append("---------------"+"--------------"+"--------------"+"-----------------"+"--------------------");
			bw.newLine();
			bw.append("Precondition       "+timeyamlparse.getListWFN_OPTCategoriesPreconditionData().get(0)  +"            "+timeyamlparse.getListWFN_OPTCategoriesPreconditionData().get(1)  +"     "+timeyamlparse.getStringWFN_OPTCategoriesPreconditionClass()  +"      "+timeyamlparse.getStringWFN_OPTCategoriesPreconditionInfo()  );
			bw.newLine();
			bw.append("ApplyLocPotKin     "+timeyamlparse.getListWFN_OPTCategoriesApplyLocPotKinData().get(0)+"            "+timeyamlparse.getListWFN_OPTCategoriesApplyLocPotKinData().get(1)+"     "+timeyamlparse.getStringWFN_OPTCategoriesApplyLocPotKinClass()+"      "+timeyamlparse.getStringWFN_OPTCategoriesApplyLocPotKinInfo());
			bw.newLine();
			bw.append("ApplyProj          "+timeyamlparse.getListWFN_OPTCategoriesApplyProjData().get(0)     +"            "+timeyamlparse.getListWFN_OPTCategoriesApplyProjData().get(1)     +"     "+timeyamlparse.getStringWFN_OPTCategoriesApplyProjClass()     +"      "+timeyamlparse.getStringWFN_OPTCategoriesApplyProjInfo()     );
			bw.newLine();
			bw.append("Un-TransSwitch     "+timeyamlparse.getListWFN_OPTCategoriesUn_TransSwitchData().get(0)+"            "+timeyamlparse.getListWFN_OPTCategoriesUn_TransSwitchData().get(1)+"     "+timeyamlparse.getStringWFN_OPTCategoriesUn_TransSwitchClass()+"      "+timeyamlparse.getStringWFN_OPTCategoriesUn_TransSwitchInfo());
			bw.newLine();
			bw.append("LagrM_comput       "+timeyamlparse.getListWFN_OPTCategoriesLagrM_computData().get(0)  +"            "+timeyamlparse.getListWFN_OPTCategoriesLagrM_computData().get(1)  +"     "+timeyamlparse.getStringWFN_OPTCategoriesLagrM_computClass()  +"      "+timeyamlparse.getStringWFN_OPTCategoriesLagrM_computInfo()  );
			bw.newLine();
			bw.append("Rho_comput         "+timeyamlparse.getListWFN_OPTCategoriesRho_computData().get(0)    +"            "+timeyamlparse.getListWFN_OPTCategoriesRho_computData().get(1)    +"     "+timeyamlparse.getStringWFN_OPTCategoriesRho_computClass()    +"      "+timeyamlparse.getStringWFN_OPTCategoriesRho_computInfo()    );
			bw.newLine();
			bw.append("Chol_comput        "+timeyamlparse.getListWFN_OPTCategoriesChol_computData().get(0)   +"            "+timeyamlparse.getListWFN_OPTCategoriesChol_computData().get(1)   +"     "+timeyamlparse.getStringWFN_OPTCategoriesChol_computClass()   +"      "+timeyamlparse.getStringWFN_OPTCategoriesChol_computInfo()   );
			bw.newLine();
			bw.append("Diis               "+timeyamlparse.getListWFN_OPTCategoriesDiisData().get(0)          +"            "+timeyamlparse.getListWFN_OPTCategoriesDiisData().get(1)          +"     "+timeyamlparse.getStringWFN_OPTCategoriesDiisClass()          +"      "+timeyamlparse.getStringWFN_OPTCategoriesDiisInfo()          );
			bw.newLine();
			bw.append("PSolv_comput       "+timeyamlparse.getListWFN_OPTCategoriesPSolv_computData().get(0)  +"            "+timeyamlparse.getListWFN_OPTCategoriesPSolv_computData().get(1)  +"     "+timeyamlparse.getStringWFN_OPTCategoriesPSolv_computClass()  +"      "+timeyamlparse.getStringWFN_OPTCategoriesPSolv_computInfo()  );
			bw.newLine();
			bw.append("Exchangecorr       "+timeyamlparse.getListWFN_OPTCategoriesExchangecorrData().get(0)  +"            "+timeyamlparse.getListWFN_OPTCategoriesExchangecorrData().get(1)  +"     "+timeyamlparse.getStringWFN_OPTCategoriesExchangecorrClass()  +"      "+timeyamlparse.getStringWFN_OPTCategoriesExchangecorrInfo()  );
			bw.newLine();
			bw.append("Pot_commun         "+timeyamlparse.getListWFN_OPTCategoriesPot_communData().get(0)    +"            "+timeyamlparse.getListWFN_OPTCategoriesPot_communData().get(1)    +"     "+timeyamlparse.getStringWFN_OPTCategoriesPot_communClass()    +"      "+timeyamlparse.getStringWFN_OPTCategoriesPot_communInfo()    );
			bw.newLine();
			bw.append("---------------"+"--------------"+"----------LAST"+"-----------------"+"--------------------");
			bw.newLine();
			bw.append(timeyamlparse.getMapLAST()+"");
			bw.newLine();
			bw.append("---------------"+"--------------"+"--------------"+"-----------------"+"--------------------");
			bw.newLine();
			bw.append("Communications     "+timeyamlparse.getListLASTClassesCommunications().get(0)       +"            "+timeyamlparse.getListLASTClassesCommunications().get(1));
			bw.newLine();
			bw.append("Convolutions       "+timeyamlparse.getListLASTClassesConvolutions().get(0)         +"            "+timeyamlparse.getListLASTClassesConvolutions().get(1));
			bw.newLine();
			bw.append("Linear Algebra     "+timeyamlparse.getListLASTClassesLinearAlgebra().get(0)        +"            "+timeyamlparse.getListLASTClassesLinearAlgebra().get(1));
			bw.newLine();
			bw.append("Other              "+timeyamlparse.getListLASTClassesOther().get(0)                +"            "+timeyamlparse.getListLASTClassesOther().get(1));
			bw.newLine();
			bw.append("Potential          "+timeyamlparse.getListLASTClassesPotential().get(0)            +"            "+timeyamlparse.getListLASTClassesPotential().get(1));
			bw.newLine();
			bw.append("Initialization     "+timeyamlparse.getListLASTClassesInitialization().get(0)       +"            "+timeyamlparse.getListLASTClassesInitialization().get(1));
			bw.newLine();
			bw.append("Finalization       "+timeyamlparse.getListLASTClassesFinalization().get(0)         +"            "+timeyamlparse.getListLASTClassesFinalization().get(1));
			bw.newLine();
			bw.append("---------------"+"--------------"+"--------------"+"-----------------"+"--------------------");
			bw.newLine();
			bw.append("Total              "+timeyamlparse.getListLASTClassesTotal().get(0)                +"            "+timeyamlparse.getListLASTClassesTotal().get(1));
			bw.newLine();
			bw.append("---------------"+"--------------"+"--------------"+"-----------------"+"--------------------");
			bw.newLine();
			bw.append("Forces             "+timeyamlparse.getListLASTCategoriesForcesData().get(0)        +"            "+timeyamlparse.getListLASTCategoriesForcesData().get(1)              +"     "+timeyamlparse.getStringLASTCategoriesForcesClass()           +"      "+timeyamlparse.getStringLASTCategoriesForcesInfo()      );
			bw.newLine();
			bw.append("Rho_comput         "+timeyamlparse.getListLASTCategoriesRho_computData().get(0)    +"            "+timeyamlparse.getListLASTCategoriesRho_computData().get(1)          +"     "+timeyamlparse.getStringLASTCategoriesForcesClass()           +"      "+timeyamlparse.getStringLASTCategoriesForcesInfo()      );
			bw.newLine();
			bw.append("PSolv_comput       "+timeyamlparse.getListLASTCategoriesPSolv_computData().get(0)  +"            "+timeyamlparse.getListLASTCategoriesPSolv_computData().get(1)        +"     "+timeyamlparse.getStringLASTCategoriesPSolv_computClass()     +"      "+timeyamlparse.getStringLASTCategoriesPSolv_computInfo());
			bw.newLine();
			bw.append("---------------"+"--------------"+"---------SUMMA"+"RY---------------"+"--------------------");
			bw.newLine();
			bw.append(timeyamlparse.getMapSUMMARY()+"");
			bw.newLine();
			bw.append("---------------"+"--------------"+"--------------"+"-----------------"+"--------------------");
			bw.newLine();
			bw.append("INIT               "+timeyamlparse.getListSUMMARYINIT().get(0)                     +"            "+timeyamlparse.getListSUMMARYINIT().get(1)    );
			bw.newLine();
			bw.append("WFN_OPT            "+timeyamlparse.getListSUMMARYWFN_OPT().get(0)                  +"            "+timeyamlparse.getListSUMMARYWFN_OPT().get(1) );
			bw.newLine();
			bw.append("LAST               "+timeyamlparse.getListSUMMARYLAST().get(0)                     +"            "+timeyamlparse.getListSUMMARYLAST().get(1)    );
			bw.newLine();
			bw.append("---------------"+"--------------"+"--------------"+"-----------------"+"--------------------");
			bw.newLine();
			bw.append("Total              "+timeyamlparse.getListSUMMARYTotal().get(0)                    +"            "+timeyamlparse.getListSUMMARYTotal().get(1)   );
			bw.newLine();
			bw.append("---------------"+"--------------"+"--------------"+"-----------------"+"--------------------");
			bw.newLine();
			bw.append("CPU Parallelism->MPI procs: "+timeyamlparse.getMapSUMMARYCPUParallelism().get("MPI procs")+" OMP thrds: "+timeyamlparse.getMapSUMMARYCPUParallelism().get("OMP thrds") );
			bw.newLine();
			bw.append("---------------"+"--------------"+"--------------"+"-----------------"+"--------------------");
			bw.newLine();
		} finally {
			if (bw != null) {
				bw.close();
			}
			if (fw != null) {
				fw.close();
			}	
		}
	
	}
}
