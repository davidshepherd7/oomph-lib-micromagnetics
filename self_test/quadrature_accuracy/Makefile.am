
#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers
# Name of executables
check_PROGRAMS=quadrature_accuracy_driver

# Sources for first executable
quadrature_accuracy_driver_SOURCES=quadrature_accuracy_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
quadrature_accuracy_driver_LDADD = -L@libdir@ -lgeneric -lmicromag $(EXTERNAL_LIBS) $(FLIBS)

AM_CPPFLAGS += -isystem @includedir@
