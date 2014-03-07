# A few file definitions
#-----------------------

SUBDIRS = self_test

# Define the sources
sources = micromagnetics_element.cc \
	micromagnetics_boundary_element.cc \
	energy_functions.cc \
	template_free_poisson.cc \
	magnetic_materials.cc \
	llg_problem.cc \
	variable_order_quadrature.cc \
	generic_poisson_problem.cc \
	boundary_element_handler.cc \
	template_free_poisson_flux.cc \
	residual_calculator.cc \
	micromag_factories.cc \
	my_generic_problem.cc \
	magnetics_helpers.cc \
	magnetostatics_calculator.cc


# Include files which shouldn't be compiled
incl_cc_files =

# Define the headers
headers = array_interpolator.h \
	micromagnetics_boundary_element.h \
	energy_functions.h \
	micromagnetics_element.h \
	generic_poisson_problem.h \
	micromagnetics_flux_element.h \
	my_general_header.h \
	my_generic_problem.h \
	llg_problem.h \
	magnetic_materials.h \
	template_free_poisson_flux.h \
	magnetics_helpers.h \
	template_free_poisson.h \
	magnetostatic_field_flux_element.h \
	magpar_requirements.h \
	mallinson_solution.h \
	variable_order_quadrature.h \
	vector_helpers.h \
	boundary_element_handler.h \
	prettyprint98.hpp \
	residual_calculator.h \
	poisson_test_problem.h \
	residual_swapping_explicit_timestepper.h \
	multi_mesh.h \
	sum_of_matrices_preconditioner.h \
	single_element_mesh.h \
	ode_problem.h \
	unsteady_heat_problem.h \
	micromag_factories.h \
	micromag_types.h \
	magnetostatics_calculator.h \
       llg_preconditioners.h

# If we have hlib then compile that header too
if OOMPH_HAS_HLIB
headers += hmatrix.h
endif


# Define name of library
libname = micromag

# Combine headers and sources
headers_and_sources = $(headers) $(sources)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = libmicromag.la

# Sources that the library depends on:
#-------------------------------------
libmicromag_la_SOURCES = $(headers_and_sources)

# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/$(libname)

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers) $(incl_cc_files)

# Required libraries
#-------------------
# Include path for library headers
AM_CPPFLAGS += -I/usr/local/include/HLib/

# Combined header file
#---------------------
# Rule for building combined header (dummy target never gets made so
# the combined header is remade every time)
$(libname).h: dummy_$(libname).h

dummy_$(libname).h:  $(headers)
	echo  "./$(libname)" $(headers) > all_$(libname).aux
	$(AWK) -f $(top_srcdir)/bin/headers.awk < \
		all_$(libname).aux > $(libname).h
	sed -e 's/</ "/g' -e 's/>/"/g' -i $(libname).h
	rm all_$(libname).aux


# Extra hook for install: Optionally replace headers by symbolic links
#---------------------------------------------------------------------
if SYMBOLIC_LINKS_FOR_HEADERS
install-data-hook:
	(cd $(library_includedir) && rm -f $(headers) $(incl_cc_files) )
	(echo "$(headers) $(incl_cc_files)" > include_files.list )
	($(top_srcdir)/bin/change_headers_to_links.sh `pwd`)
	($(LN_S) `cat include_files.list.aux` $(library_includedir) )
	(rm -r include_files.list.aux include_files.list )
else
install-data-hook:
endif

# Cleanup: Kill automatically generated combined header file
clean-local:
	rm -f $(libname).h *~
distclean-local:
	rm -f $(libname).h *~
