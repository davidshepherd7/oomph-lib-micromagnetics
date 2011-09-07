# Name of executables 
noinst_PROGRAMS=demag_driver


# Sources for first executable
demag_driver_SOURCES = demag_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
demag_driver_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


INCLUDES = -I@includedir@
