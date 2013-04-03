# Name of executables
noinst_PROGRAMS=generate_tetgen_sphere_input

# Sources for first executable
generate_tetgen_sphere_input_SOURCES = generate_tetgen_sphere_input.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
# sphere_driver_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
generate_tetgen_sphere_input_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -I@includedir@
