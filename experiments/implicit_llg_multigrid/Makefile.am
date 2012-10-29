# Name of executables
noinst_PROGRAMS=multigrid_experiments

# Sources for first executable
multigrid_experiments_SOURCES=multigrid_experiments.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
multigrid_experiments_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -isystem @includedir@
