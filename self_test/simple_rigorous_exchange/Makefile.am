#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Names of executables
check_PROGRAMS=simple_rigorous_exchange_test

# Sources
simple_rigorous_exchange_test_SOURCES = simple_rigorous_exchange_test.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
simple_rigorous_exchange_test_LDADD = -L@libdir@ -lgeneric -lmicromag $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -isystem @includedir@
