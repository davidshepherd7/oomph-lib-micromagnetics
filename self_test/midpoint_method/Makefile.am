#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Names of executables
check_PROGRAMS=midpoint_method_test

# Sources
midpoint_method_test_SOURCES = midpoint_method_test.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
midpoint_method_test_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -isystem @includedir@
