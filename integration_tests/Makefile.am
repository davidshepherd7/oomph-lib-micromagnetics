# Name of executables
noinst_PROGRAMS=integration_tests_driver


# Sources for first executable
integration_tests_driver_SOURCES = integration_tests_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
integration_tests_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


INCLUDES = -I@includedir@