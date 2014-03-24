#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=hlib_matrix_generate_driver

# Sources for first executable
hlib_matrix_generate_driver_SOURCES=hlib_matrix_generate_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
hlib_matrix_generate_driver_LDADD = -L@libdir@ -lgeneric -lmicromag -lhmatrix -lunsteady_heat $(EXTERNAL_LIBS) $(FLIBS)

AM_CPPFLAGS += -isystem @includedir@ -I /usr/local/include/HLib/
