# Names of executables
noinst_PROGRAMS=llg_exchange_test

# Sources
llg_exchange_test_SOURCES = llg_exchange_test.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
llg_exchange_test_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

INCLUDES = -isystem @includedir@
