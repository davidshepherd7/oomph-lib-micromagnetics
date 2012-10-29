# Name of executables
noinst_PROGRAMS=combined_hybrid_driver


# Sources for first executable
combined_hybrid_driver_SOURCES = combined_hybrid_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
combined_hybrid_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


INCLUDES = -I@includedir@
