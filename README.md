oomph-lib-micromagnetics
========================

Code to extend oomph-lib to solve micromagnetics problems.

To use you will need a full install of oomph-lib (see https://github.com/davidshepherd7/oomph-lib and http://oomph-lib.maths.man.ac.uk/doc/html/index.html) then:

* Place this repository in a subdirectory of `user_drivers`
* Build oomph-lib (using `autogen.sh` in the root directory)
* Do `make && make install` in root of this repository (this installs the libraries in oomph-libs lib directory)

The main binaries of interest are `control_scripts/semi_implicit_mm_driver/semi_implicit_mm_driver` and `control_scripts/llg_driver/llg_driver`. Run make in their directories to build them. Run the binaries with --help to see a list of command line flags.

There are also a number of useful scripts in `control_scripts`. In particular `mesh_construction.py` should be run to create the unstructured meshes (requires `tetgen` and `triangle`).
