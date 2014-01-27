# Name of executables
noinst_PROGRAMS=rigorous_degree_test_driver


# Sources for first executable
rigorous_degree_test_driver_SOURCES = rigorous_degree_test_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
rigorous_degree_test_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


AM_CPPFLAGS += -I@includedir@
