This document describes how a single-particle project is defined in the SIMPLE software suite. A project is represented by a binary file (*.simple) divided into segments, where each segment correpsonds to:

Segment 1: StkInfo
contains per-micrograph stack info (many lines) or merged single-stack info (one line)

Segment 2: Ptcl2DInfo
contains per-particle information produced by prime2D

Segment 3: ClsInfo
contains per-cluster information about the 2D clusters produced by prime2D

Segment 4: Cls3DInfo
contains per-cluster information produced by ini3D or prime3D on the class averages produced by prime2D

Segment 5: Ptcl3DInfo
contains per-particle information produced by prime3D

Segment 6: Vol3DInfo
contains per-volume information produced by prime3D: volume filename, anisotropic filter file, fsc_file etc

Segment 7: FSC3DInfo
contains the FSC information produced by prime3D

Segment 8: Het3DInfo
contains per-particle information produced by heterogeneity analysis

Segment 9: HetVol3DInfo
contains per-volume information produced by heterogeneity analysis

Segment 10: HetFSC3DInfo
contains the FSC information produced by heterogeneity analysis

Each segment is represented in the run-time environment by an oris object. The set of oris objects that map to the segments will be implemented in a class (simple_sp_project). The class that deals with turning orientation data into raw binary byte arrays (for writing to disk) and reads raw binary byte arrays from the disk segments, converting them to strings and finally oris objects is simple_binoris.  