# Name of executables
noinst_PROGRAMS=hybrid_boundary_element_driver


# Sources for first executable
hybrid_boundary_element_driver_SOURCES = hybrid_boundary_element_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
hybrid_boundary_element_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


INCLUDES = -I@includedir@
