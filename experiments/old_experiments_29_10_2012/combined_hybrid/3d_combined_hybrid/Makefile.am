# Name of executables
noinst_PROGRAMS=3d_combined_hybrid_driver


# Sources for first executable
3d_combined_hybrid_driver_SOURCES = 3d_combined_hybrid_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
3d_combined_hybrid_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


INCLUDES = -I@includedir@
