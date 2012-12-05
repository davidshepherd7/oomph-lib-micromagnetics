# Name of executables
noinst_PROGRAMS=single_element_bem

# Sources for first executable
single_element_bem_SOURCES=single_element_bem.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
single_element_bem_LDADD = -L@libdir@ -lgeneric -lpoisson $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -I@includedir@
