# Name of executables 
noinst_PROGRAMS=discontinuity_tests_driver


# Sources for first executable
discontinuity_tests_driver_SOURCES=discontinuity_tests_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
discontinuity_tests_driver_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


INCLUDES = -I@includedir@
