# Name of executables
noinst_PROGRAMS=circle_bem

# Sources for first executable
circle_bem_SOURCES=circle_bem.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
circle_bem_LDADD = -L@libdir@ -lgeneric -lpoisson $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -I@includedir@
