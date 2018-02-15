
# SIMPLE API documentation {mainpage}

## Introduction

  This document describes the source code for SIMPLE.
  Visit the [SIMPLE homepage][1] for more information.


## Developer documentation

  You've reached the automically generated documentation of SIMPLE.
  This document contains information about architecture and code structure.
  It is attended for developers wanting to understand or contribute to SIMPLE.

## Modules

  Documentation is generate from comments placed in all parts of the code.
  But you will also find some groups describing specific functional parts:
       -  Definitions
       -  Utilities
       -  Image Operations
       -  Optimisation and search
       -  Command-line and program interfaces

### Definitions

    * singletons: defs, params, args
    * abstract factory
    
### Utilities

   * Math,
   * array,
   * statistics,
   * linked-lists,
   * system functions
   * file I/O
   * image format handling
   

### Image Operations
   * Image class
   * Image operations and transforms
   * Volume operations
   * Projection operations
   * Filtering
   * Fourier transform
   * Image reconstruction


### Optimisation and search

   * Core class
     + Orientation Alignment Optimisation
     + Volume registration
   * Optimisation methods
     + OASIS
     + Simplex
     + Common-line
     + Simulated annealing
     + Particle swarm
     + Brute-force

### Command-line and program interfaces

   * Command-line arguments
     + Files: stk= outstk= deftab=
     + Variables: smpd=
     + Options: verbose='yes' debug='yes'
   * Programs (prg=" ")
     + ctffind
     + unblur
     + convert
     + prime2D
     + prime3D
     + recvol


[1]: https://simplecryoem.com  "Elmlund and Elmlund Lab's --  SIMPLE homepage"


