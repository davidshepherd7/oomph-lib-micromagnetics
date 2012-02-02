# Name of executables
noinst_PROGRAMS=generate_gauss_rules_driver


# Sources for first executable
generate_gauss_rules_driver_SOURCES = generate_gauss_rules_driver.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
generate_gauss_rules_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


INCLUDES = -I@includedir@
