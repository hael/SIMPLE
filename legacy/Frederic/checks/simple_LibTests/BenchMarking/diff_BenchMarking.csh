#!/bin/bash
set TARGET_HOME_SIDE = /home/frederic/Monash/SourceCode/Simple/Restructured/HansVersion

set DIFF_PATH = $TARGET_HOME_SIDE/Benchmark/BenchMark_tmpData

echo "::::::::::::::"
echo "BenchMrkDtAnlzr-incFile.f90"
echo "::::::::::::::"
diff BenchMrkDtAnlzr-incFile.f90 $DIFF_PATH
echo "::::::::::::::"
echo "bench.py"
echo "::::::::::::::"
diff ../bench.py $DIFF_PATH
echo "::::::::::::::"
echo "create_simulData.csh"
echo "::::::::::::::"
diff create_simulData.csh $DIFF_PATH
echo "::::::::::::::"
echo "diff_BenchMarking.csh"
echo "::::::::::::::"
diff diff_BenchMarking.csh $DIFF_PATH
echo "::::::::::::::"
echo "generatorData_Vol.csh"
echo "::::::::::::::"
diff generatorData_Vol.csh $DIFF_PATH
echo "::::::::::::::"
echo "generatorData_Stk.csh"
echo "::::::::::::::"
diff generatorData_Stk.csh $DIFF_PATH
echo "::::::::::::::"
echo "generator.py"
echo "::::::::::::::"
diff generator.py $DIFF_PATH
echo "::::::::::::::"
echo "header.py"
echo "::::::::::::::"
diff header.py $DIFF_PATH
echo "::::::::::::::"
echo "helper.py"
echo "::::::::::::::"
diff helper.py $DIFF_PATH
echo "::::::::::::::"
echo "parse_simulData.csh"
echo "::::::::::::::"
diff parse_simulData.csh $DIFF_PATH
echo "::::::::::::::"
echo "plotter.py"
echo "::::::::::::::"
diff plotter.py $DIFF_PATH
echo "::::::::::::::"
echo "simple_BenchMrkDtAnlzr.f90"
echo "::::::::::::::"
diff simple_BenchMrkDtAnlzr.f90 $DIFF_PATH
