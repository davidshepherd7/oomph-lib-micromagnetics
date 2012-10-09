# Name of executables
noinst_PROGRAMS=generic_poisson_test

# Sources for first executable
generic_poisson_test_SOURCES=generic_poisson_test.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
generic_poisson_test_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# using -isystem rather than -I suppresses (most) warnings from files in includedir
INCLUDES = -isystem @includedir@