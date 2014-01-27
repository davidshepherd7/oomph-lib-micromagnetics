# Name of executables
noinst_PROGRAMS=sum_matrix_operation_tests


# Sources for first executable
sum_matrix_operation_tests_SOURCES = sum_matrix_operation_tests.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
sum_matrix_operation_tests_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

AM_CPPFLAGS += -I@includedir@
