#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=vector_helpers_test

# Sources for first executable
vector_helpers_test_SOURCES=vector_helpers_test.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
vector_helpers_test_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# using -isystem rather than -I suppresses (most) warnings from files in includedir
AM_CPPFLAGS += -isystem @includedir@
