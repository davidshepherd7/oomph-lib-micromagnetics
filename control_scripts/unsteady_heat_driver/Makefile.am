# Name of executables
noinst_PROGRAMS=unsteady_heat_driver

# Sources for first executable
unsteady_heat_driver_SOURCES=unsteady_heat_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unsteady_heat_driver_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -isystem @includedir@
