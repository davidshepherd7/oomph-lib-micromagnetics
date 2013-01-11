# Name of executables
noinst_PROGRAMS=spatially_constant_m_length_variations

# Sources for first executable
spatially_constant_m_length_variations_SOURCES=spatially_constant_m_length_variations.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
spatially_constant_m_length_variations_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -isystem @includedir@
