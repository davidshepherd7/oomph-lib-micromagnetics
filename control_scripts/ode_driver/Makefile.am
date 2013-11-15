noinst_PROGRAMS= ode_driver

ode_driver_SOURCES = ode_driver.cc
ode_driver_LDADD = -L@libdir@ -lgeneric -lmicromag $(EXTERNAL_LIBS) $(FLIBS)

# Include path for library headers: All library headers live in
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
INCLUDES = -I@includedir@
