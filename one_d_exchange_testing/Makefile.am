# Name of executables 
noinst_PROGRAMS=one_d_exchange_testing_driver_driver


# Sources for first executable
one_d_exchange_testing_driver_driver_SOURCES = one_d_micromag_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
one_d_exchange_testing_driver_driver_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


INCLUDES = -I@includedir@
