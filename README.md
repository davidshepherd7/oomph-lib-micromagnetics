oomph-lib-micromagnetics
========================

Code to extend oomph-lib to solve micromagnetics problems.

To use you will need a full install of oomph-lib (see https://github.com/davidshepherd7/oomph-lib and http://oomph-lib.maths.man.ac.uk/doc/html/index.html) then:

* Add ./etc/ directory to your python path variable
* Place this repository in a subdirectory of `user_drivers`
* Build oomph-lib (using `autogen.sh` in the root directory)
* Do `make && make install` in root of this repository (this installs the libraries in oomph-libs lib directory)
* Create the meshes with ./control_scripts/mesh_generation.py

The main binary of interest is `control_scripts/driver/driver`. Run with --help to see a list of command line flags.

There are also a number of useful scripts in `control_scripts`. In particular `mesh_construction.py` can be used to create various unstructured meshes (requires `tetgen` and `triangle`).


Note that the tetgen version for which the self tests are written is:
`Version 1.4.3 (January 19, 2011)`. Newer versions have completely changed
the way they generate meshes and so cannot (currently) be used by oomph-lib.
