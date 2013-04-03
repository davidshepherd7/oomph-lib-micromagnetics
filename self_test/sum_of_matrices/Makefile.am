#Include commands common to every self test Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=sum_of_matrices_test

# Sources for first executable
sum_of_matrices_test_SOURCES=sum_of_matrices_test.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
sum_of_matrices_test_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# using -isystem rather than -I suppresses (most) warnings from files in includedir
INCLUDES = -isystem @includedir@
