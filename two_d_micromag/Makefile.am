# Name of executables 
noinst_PROGRAMS=two_d_micromag_driver


# Sources for first executable
two_d_micromag_driver_SOURCES = two_d_micromag_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_micromag_driver_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


INCLUDES = -I@includedir@
