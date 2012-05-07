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

#include <ostream>
#include <utility>


namespace oomph
{

  // standard outputters for doublevector etc Should be "safe" but wrong for
  // distributed version - no idea what order it will output in but shouldn't
  // segfault.
  std::ostream& operator<<(std::ostream& output, const DoubleVector& dv)
  {
#ifdef MPI
    throw OomphLibWarning("I didn't write this with distributed vectors in mind, it might work but I have no idea.",
			  "operator<<(ostream& output, const DoubleVector& dv)",
			  OOMPH_EXCEPTION_LOCATION);
#endif
    output << "[";
    for(unsigned i=0; i<dv.nrow_local()-1; i++)
      output << dv[i] << ", ";
    // Output last entry seperately to get the commas right
    output << dv[dv.nrow_local()-1] << "]";
    return output;
  }



}

#endif
