# Name of executables
noinst_PROGRAMS=mumag4_driver

# Sources for first executable
mumag4_driver_SOURCES = mumag4_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
mumag4_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -I@includedir@
