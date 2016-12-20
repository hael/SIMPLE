#!/bin/bash

set DIFF_PATH = /home/frederic/Monash/SourceCode/Simple/Restructured/HansVersion/Simple_Restruct.projet/simple_utils/common

echo "::::::::::::::"
echo "Makefile_target"
echo "::::::::::::::"
diff Makefile_target                  $DIFF_PATH
echo "::::::::::::::"
echo "greeting_version.f90"
echo "::::::::::::::"
diff greeting_version.f90             $DIFF_PATH
echo "::::::::::::::"
echo "initialu.f90"
echo "::::::::::::::"
diff initialu.f90                     $DIFF_PATH
#diff matrixGetter.f90       $DIFF_PATH
echo "::::::::::::::"
echo "simple_random.f90"
echo "::::::::::::::"
diff simple_random.f90                $DIFF_PATH
echo "::::::::::::::"
echo "simple_sorting.f90"
echo "::::::::::::::"
diff simple_sorting.f90               $DIFF_PATH
echo "::::::::::::::"
echo "simple_SU3_tester.f90"
echo "::::::::::::::"
diff simple_SU3_tester.f90            $DIFF_PATH
echo "::::::::::::::"
echo "simple_testfuns.f90"
echo "::::::::::::::"
diff simple_testfuns.f90              $DIFF_PATH
echo "::::::::::::::"
echo "simple_textHandler.f90"
echo "::::::::::::::"
diff simple_textHandler.f90           $DIFF_PATH
echo "::::::::::::::"
echo "simple_timing.f90"
echo "::::::::::::::"
diff simple_timing.f90                $DIFF_PATH
echo "::::::::::::::"
echo "su2random.f90"
echo "::::::::::::::"
diff su2random.f90                    $DIFF_PATH
echo "::::::::::::::"
echo "su3random.f90"
echo "::::::::::::::"
diff su3random.f90                    $DIFF_PATH
echo "::::::::::::::"
echo "timestamp.f90"
echo "::::::::::::::"
diff timestamp.f90                    $DIFF_PATH
echo "::::::::::::::"
echo "simple_yaml_open-incFl.f90"
echo "::::::::::::::"
diff simple_yaml_open-incFl.f90 $DIFF_PATH
echo "::::::::::::::"
echo "simple_yaml_output.f90"
echo "::::::::::::::"
diff simple_yaml_output.f90           $DIFF_PATH
echo "::::::::::::::"
echo "simple_yaml_strings.f90"
echo "::::::::::::::"
diff simple_yaml_strings.f90          $DIFF_PATH
#diff timming_c.cpp          $DIFF_PATH
#diff timming_c.cpp~         $DIFF_PATH
