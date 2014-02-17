# Name of executables
noinst_PROGRAMS=singular_quadrature_test_driver

# Sources for first executable
singular_quadrature_test_driver_SOURCES = singular_quadrature_test_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
singular_quadrature_test_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

AM_CPPFLAGS += -I@includedir@
