# Name of executables
noinst_PROGRAMS=preconditioner_test_driver

# Sources for first executable
preconditioner_test_driver_SOURCES=preconditioner_test_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
preconditioner_test_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -I@includedir@
