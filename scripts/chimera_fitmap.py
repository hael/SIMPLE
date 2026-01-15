from chimera import runCommand as rc
import os
rc("open vol_reference.mrc")
rc("open vol_target.mrc")
rc("fitmap #1 #0" )
rc("vop resample #1 onGrid #0 ")
rc("volume #2 save vol_target_docked.mrc")

