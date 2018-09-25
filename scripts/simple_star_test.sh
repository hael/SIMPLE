!#/usr/bin/bash

SIMPLE_PATH=~/el85_scratch/eagerm/SIMPLE3.0/build
module load git cmake fftw  gctf/1.06_cuda8 chimera eman/2.2
export PATH=$SIMPLE_PATH/bin:$PATH
cd ~/el85_scratch/stars_from_matt

[ ! -d simple_convert_star ] && mkdir simple_convert_star

## IMPORTING GCTF PROCESSED MICROGRAPHS
cd ~/el85_scratch/stars_from_matt/simple_convert_star
[ -d micrographs_test1 ] && rm -rf micrographs_test1
simple_exec prg=new_project projname=micrographs_test1

cd micrographs_test1
simple_exec prg=importstar_project starfile=../../micrographs_all_gctf.star
simple_exec prg=print_project_info
simple_exec prg=export_starproject export_type=micrographs starfile=simple_exported.star

cd ..
simple_exec prg=new_project projname=micrographs_test2
cd micrographs_test2
ls -1 CtfFind/job002/Micrographs/00*.mrc | sort -n > filetab.txt
simple_exec prg=import_movies  filetab=filetab.txt cs=2.7 kv=300 fraca=0.1 smpd=3.5 ctf=yes
simple_exec prg=print_project_info



## Cryosparc Movies
cd ~/el85_scratch/stars_from_matt/simple_convert_star
simple_exec prg=new_project projname=cryosparc_movies
cd cryosparc_movies
echo "../../cryosparc_exp000325_002.mrc" > filetab.txt
simple_exec prg=import_movies filetab=filetab.txt cs=2.7 kv=300 fraca=0.1 smpd=1 ctf=yes
simple_exec prg=print_project_info
simple_exec prg=export_starproject export_type=micrographs starfile=simple_exported.star

       

cd ~/el85_scratch/stars_from_matt/simple_convert_star
simple_exec prg=new_project projname=CtfFind
cd CtfFind
simple_exec prg=importstar_project starfile=../../CtfFind/GCTF_Initial/micrographs_ctf.star
simple_exec prg=print_project_info


cd ~/el85_scratch/stars_from_matt/simple_convert_star
simple_exec prg=new_project projname=Class2D
cd Class2D
simple_exec prg=importstar_project starfile=../../Class2D/Class2D_1st
#run_it*_(data|model|optimiser|sampling).star run_it*_classes.mrcs
simple_exec prg=print_project_info

 
cd ~/el85_scratch/stars_from_matt/simple_convert_star
simple_exec prg=new_project projname=InitialModel
cd InitialModel
simple_exec prg=importstar_project starfile=../../InitialModel/ 
simple_exec prg=print_project_info


cd ~/el85_scratch/stars_from_matt/simple_convert_star
simple_exec prg=new_project projname=ManualPick
cd ManualPick
# simple_exec prg=importstar_project starfile=../../ManualPick/ManualPick/coords_suffix_manualpick.star --> this contains one line: ./CtfFind/job002/micrographs_ctf.star
simple_exec prg=importstar_project starfile=../../ManualPick/ManualPick/micrographs_selected.star
simple_exec prg=print_project_info


## Extract -- Relion does not create stacks by default
##
## module load relion/2.1
## mkdir Stack_Extract_Micrographs
## ln -s `pwd`/Micrographs Extract/job006/
## (cd Micrographs; for i in *.mrc; ln -s $i ${i}s; done)
## relion_stack_create --i Extract/job006/particles.star --o Stack_Extract --split_per_micrograph TRUE --apply_transformation TRUE           


cd ~/el85_scratch/stars_from_matt/simple_convert_star
simple_exec prg=new_project projname=Extract
cd Extract
simple_exec prg=importstar_project starfile=../../Extract/364Box_Extract_LocalCTF/particles.star 
simple_exec prg=print_project_info
simple_exec prg=import_particles stktab=../filetab-stardoc.txt oritab=../oritab-stardoc.txt  ctf=yes smpd=3.5

cd ~/el85_scratch/stars_from_matt/simple_convert_star
simple_exec prg=new_project projname=ExtractStk
cd ExtractStk
for i in ../../Stack_Extract_Micrographs/*.mrcs; do
    simple_exec prg=import_cavgs stk=../../Stack_Extract_Micrographs/0001.mrcs  smpd=3.5
done
simple_exec prg=print_project_info



# cd ~/el85_scratch/stars_from_matt/simple_convert_star
# simple_exec prg=new_project projname=MaskCreate
# cd MaskCreate
# simple_exec prg=importstar_project starfile=../../MaskCreate/  

  
cd ~/el85_scratch/stars_from_matt/simple_convert_star
simple_exec prg=new_project projname=Select
cd Select
simple_exec prg=importstar_project starfile=../../Select/1stCut/particles.star
simple_exec prg=importstar_project starfile=../../Select/1stCut/class_averages.star
            
cd ~/el85_scratch/stars_from_matt/simple_convert_star
simple_exec prg=new_project projname=Import
cd Import
simple_exec prg=importstar_project starfile=../../Import/Import/micrographs.star
simple_exec prg=importstar_project starfile=../../Import/GautoMatch_1st/coords_suffix_automatch.star
simple_exec prg=importstar_project starfile=../../Import/LocalCTF_ParticlePositionImport/coords_suffix_local.star


    
cd ~/el85_scratch/stars_from_matt/simple_convert_star
simple_exec prg=new_project projname=Refine3D
cd Refine3D
simple_exec prg=importstar_project starfile=../../Refine3D/Refine3D_1st/run_ct19_data.star
# run_ct19_it019_data.star
# run_it000_data.star
cd ~/el85_scratch/stars_from_matt/simple_convert_star
simple_exec prg=new_project projname=Refine3D_Mask_2nd
cd Refine3D_Mask_2nd
simple_exec prg=importstar_project starfile=../../Refine3D/Refine3D_Mask_2nd/run_it018_data.star





cd ~/el85_scratch/stars_from_matt/simple_convert_star
simple_exec prg=new_project projname=PostProcess
cd PostProcess
simple_exec prg=importstar_project starfile=../../PostProcess/Post/postprocess.star
simple_exec prg=print_project_info

