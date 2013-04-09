# Name of executables
noinst_PROGRAMS=midpoint_llg

# Sources for first executable
midpoint_llg_SOURCES=midpoint_llg.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
midpoint_llg_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -isystem @includedir@
