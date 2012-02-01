# Name of executables
noinst_PROGRAMS=generate_quadrature


# Sources for first executable
generate_quadrature_SOURCES = generate_quadrature.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
generate_quadrature_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


INCLUDES = -I@includedir@
