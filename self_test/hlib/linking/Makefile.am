#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=hlib_test_driver

# Sources for first executable
hlib_test_driver_SOURCES=hlib_test_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
hlib_test_driver_LDADD = -L@libdir@ -lgeneric -lmicromag $(EXTERNAL_LIBS) $(FLIBS)

AM_CPPFLAGS += -isystem @includedir@
