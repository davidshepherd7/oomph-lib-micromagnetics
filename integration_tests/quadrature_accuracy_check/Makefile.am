# Name of executables
noinst_PROGRAMS=quadrature_accuracy_test_driver

# Sources for first executable
quadrature_accuracy_test_driver_SOURCES = quadrature_accuracy_test_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
quadrature_accuracy_test_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -I@includedir@