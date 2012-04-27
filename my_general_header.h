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

#include <fenv.h>

#include "./vector_helpers.h"

#endif
