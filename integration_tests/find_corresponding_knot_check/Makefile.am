# Name of executables
noinst_PROGRAMS=find_corresponding_knot_check_driver


# Sources for first executable
find_corresponding_knot_check_driver_SOURCES = find_corresponding_knot_check_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
find_corresponding_knot_check_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


INCLUDES = -I@includedir@
