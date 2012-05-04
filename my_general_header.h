#ifndef OOMPH_MY_GENERAL_HEADER_H
#define OOMPH_MY_GENERAL_HEADER_H

/*
 A header for all my debug and output stuff.
*/


// Include the appropriate version of the pretty print header depending on if we
// are using c++0x or not
#ifdef __GXX_EXPERIMENTAL_CXX0X__
  #include "prettyprint.hpp"
#else
  #include "prettyprint98.hpp"
#endif

// Floating point debugging
#include <fenv.h>

// Quicker to use vector functions
#include "./vector_helpers.h"

// All my micromag headers and .cc
#include "./micromagnetics_element.h"
#include "./micromagnetics_element.cc"
#include "./micromagnetics_boundary_element.h"
#include "./micromagnetics_boundary_element.cc"
#include "./micromagnetics_flux_element.h"

// Basic oomph-lib headers
#include "generic.h"

#endif
