
#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers
# Name of executables
check_PROGRAMS=brick2tet_driver

# Sources for first executable
brick2tet_driver_SOURCES=brick2tet_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
brick2tet_driver_LDADD = -L@libdir@ -lgeneric -lmicromag $(EXTERNAL_LIBS) $(FLIBS)

AM_CPPFLAGS += -isystem @includedir@
