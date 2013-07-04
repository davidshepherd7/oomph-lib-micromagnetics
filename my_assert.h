#ifndef ASSERT_H
#define ASSERT_H
// // Shamelessly stolen from boost!

#include <string>
#include "../../src/generic/oomph_utilities.h"
#include "../../src/generic/oomph_definitions.h"

// Make sure the macros don't exist already
#if defined(my_assert) || defined(my_assert_msg)
// Intentionally cause a compiler error if they are defined because
// something is likely to go horribly wrong!
###########
#endif


// Only define all this code when PARANOID mode is on
#ifdef PARANOID

// Define functions for what to do when an my_assertion fails. For now just
// pass everything into an OomphLibError so that it all acts the same as it
// used to.
namespace oomph
{
 namespace assert
 {
  using namespace StringConversion;

  /// Error function used by MY_ASSERT macro
  void assertion_failed(char const * expr,
                        char const * function, char const * file, long line)
  {
   throw OomphLibError("Assertion " + std::string(expr) + " failed!",
                       std::string(function),
                       (std::string(file) + ":" + to_string(line)).c_str());
  }

  /// Error function used by ASSERT_MSG macro
  void assertion_failed_msg(char const * expr, const std::string& msg,
                            char const * function, char const * file, long line)
  {
   throw OomphLibError("Assertion " + std::string(expr)
                       + " failed!\n" + std::string(msg),
                       std::string(function),
                       (std::string(file) + ":" + to_string(line)).c_str());
  }

 }
}

/// Call oomph::assertion_failed if (expr) does not evaluate to true
#define my_assert(expr) ((expr) \
  ? ((void)0) \
  : ::oomph::assert::assertion_failed(#expr, OOMPH_CURRENT_FUNCTION, __FILE__, __LINE__))

/// Call oomph::my_assertion_failed_msg if (expr) does not evaluate to true
#define my_assert_msg(expr, msg) ((expr) \
  ? ((void)0) \
  : ::oomph::assert::assertion_failed_msg(#expr, msg, OOMPH_CURRENT_FUNCTION, __FILE__, __L


// Otherwise disable the macros
#else // (else for #ifdef PARANOID)
#define my_assert(expr) ((void)0)
#define my_assert_msg(expr, msg) ((void)0)

#endif // end of #ifdef PARANOID

#endif // end of include guard
