# Name of executables
noinst_PROGRAMS=llg_driver

# Sources for first executable
llg_driver_SOURCES=llg_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
llg_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -isystem @includedir@
