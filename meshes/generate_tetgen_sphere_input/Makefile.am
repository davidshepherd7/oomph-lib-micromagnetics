# Name of executables
noinst_PROGRAMS=generate_tetgen_sphere_input

generate_tetgen_sphere_input_SOURCES = generate_tetgen_sphere_input.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
generate_tetgen_sphere_input_LDADD = -L@libdir@ $(EXTERNAL_LIBS) $(FLIBS)

AM_CPPFLAGS += -I@includedir@
