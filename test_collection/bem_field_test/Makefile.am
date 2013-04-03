#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Names of executables
noinst_PROGRAMS=bem_field_test generate_tetgen_sphere_input

# Sources
bem_field_test_SOURCES = bem_field_test.cc
generate_tetgen_sphere_input_SOURCES = generate_tetgen_sphere_input.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
bem_field_test_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
generate_tetgen_sphere_input_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -isystem @includedir@
