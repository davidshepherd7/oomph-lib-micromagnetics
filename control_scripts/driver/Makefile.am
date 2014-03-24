# Name of executables
noinst_PROGRAMS=driver

# Sources for driver
driver_SOURCES=driver.cc

# Include hlib if we have it
HMATRIX_LIBS =
if OOMPH_HAS_HLIB
HMATRIX_LIBS += -lhmatrix
AM_CPPFLAGS += -I/usr/local/include/HLib/
endif


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
driver_LDADD = -L@libdir@ -lgeneric -lmicromag -lunsteady_heat $(HMATRIX_LIBS) $(EXTERNAL_LIBS) $(FLIBS)
AM_CPPFLAGS += -isystem @includedir@
