# Name of executables
noinst_PROGRAMS=triangular_elements_driver


# Sources for first executable
triangular_elements_driver_SOURCES = triangular_elements_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
triangular_elements_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


INCLUDES = -I@includedir@
