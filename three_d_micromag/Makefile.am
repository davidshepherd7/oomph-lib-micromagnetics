# Name of executables
noinst_PROGRAMS=three_d_micromag_driver_driver


# Sources for first executable
three_d_micromag_driver_driver_SOURCES = two_d_micromag_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
three_d_micromag_driver_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


INCLUDES = -I@includedir@
