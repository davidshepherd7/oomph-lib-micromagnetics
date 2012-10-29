# Name of executables
noinst_PROGRAMS=one_d_demag_testing_driver


# Sources for first executable
one_d_demag_testing_driver_SOURCES = one_d_demag_testing_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
one_d_demag_testing_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


INCLUDES = -I@includedir@
