#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=mallinson_test

# Sources for first executable
mallinson_test_SOURCES=mallinson_test.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
mallinson_test_LDADD = -L@libdir@ -lmicromag -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# using -isystem rather than -I suppresses (most) warnings from files in includedir
AM_CPPFLAGS += -isystem @includedir@
