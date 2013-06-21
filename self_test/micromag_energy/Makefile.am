#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=micromag_energy_test

# Sources for first executable
micromag_energy_test_SOURCES=micromag_energy_test.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
micromag_energy_test_LDADD = -L@libdir@ -lmicromag -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# using -isystem rather than -I suppresses (most) warnings from files in includedir
INCLUDES = -isystem @includedir@
