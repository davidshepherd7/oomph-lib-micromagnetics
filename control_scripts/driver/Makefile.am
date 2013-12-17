# Name of executables
noinst_PROGRAMS=driver

# Sources for first executable
driver_SOURCES=driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
driver_LDADD = -L@libdir@ -lgeneric -lmicromag $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -isystem @includedir@
