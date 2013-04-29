# Name of executables
noinst_PROGRAMS=implicit_llg_driver

# Sources for first executable
implicit_llg_driver_SOURCES=implicit_llg_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
implicit_llg_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -isystem @includedir@
