# Name of executables
noinst_PROGRAMS=sphere_driver

# Sources for first executable
sphere_driver_SOURCES = sphere_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
sphere_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -I@includedir@
