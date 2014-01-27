# Name of executables
noinst_PROGRAMS=triangular_integration_driver

# Sources for first executable
triangular_integration_driver_SOURCES = triangular_integration_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
triangular_integration_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

AM_CPPFLAGS += -I@includedir@
