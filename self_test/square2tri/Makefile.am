
#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers
# Name of executables
check_PROGRAMS=square2tri_driver

# Sources for first executable
square2tri_driver_SOURCES=square2tri_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
square2tri_driver_LDADD = -L@libdir@ -lgeneric -lmicromag -lunsteady_heat $(EXTERNAL_LIBS) $(FLIBS)

AM_CPPFLAGS += -isystem @includedir@
