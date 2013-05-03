# Name of executables
noinst_PROGRAMS=semi_implicit_mm_driver

# Sources for first executable
semi_implicit_mm_driver_SOURCES=semi_implicit_mm_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
semi_implicit_mm_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -isystem @includedir@
