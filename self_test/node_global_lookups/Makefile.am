#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=node_lookup_test

# Sources for first executable
node_lookup_test_SOURCES=node_lookup_test.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
node_lookup_test_LDADD = -L@libdir@ -lpoisson -lgeneric -lmicromag $(EXTERNAL_LIBS) $(FLIBS)

# using -isystem rather than -I suppresses (most) warnings from files in includedir
AM_CPPFLAGS += -isystem @includedir@
