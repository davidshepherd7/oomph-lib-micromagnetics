

// Test chosen integration scheme on a variety of functions with known analytic values.
// C++0x only (used lambda functions to simplify defining the long list of testing functions).

#include "../variable_quadrature.h"
#include <functional>
#include <cmath>

using namespace oomph;

namespace oomph
{

  // Integrate a given function using the given scheme
  // (no check is made that the scheme integrates between -1 and +1, watch out.)
  double integrate_function_1D(const std::function<double(double)> &function,
			       const Integral &integral)
  {
    double integral_value = 0;

    // Loop over integration points
    unsigned n_weights = integral.nweight();
    for(unsigned i=0; i<n_weights; i++)
      {
	// Get knot location
	double x = integral.knot(i,0);

	// Add contribution from this knot
	integral_value += function(x) * integral.weight(i);
      }

    return integral_value;
  }

  // Structure to store the information about each testing function
  struct TestingFn
  {
    std::function<double(double)> Function;
    double Answer;

    TestingFn(const std::function<double(double)> &function, const double &answer)
    {
      Function = function;
      Answer = answer;
    }
  };

}

int main()
{

  // Create a list of testing functions and their exact solutions
  std::vector<TestingFn > function_list =
    {
      // Some very simple functions:
      TestingFn([] (double x) {return 1.0;}, 2.0),
      TestingFn([] (double x) {return x;}, 0.0),

      // Functions from Trefethen2008:
      TestingFn([] (double x) {return pow(x,20);}, 2.0/21),
      TestingFn([] (double x) {return exp(x);}, 2.3504023872876029138),
      TestingFn([] (double x) {return exp(-pow(x,2));}, 1.4936482656248540508),
      TestingFn([] (double x) {return 1/(1 + 16*x*x);}, 0.66290883183401623253),
      TestingFn([] (double x) {return exp(-1/(x*x));}, 0.17814771178156069019),
      TestingFn([] (double x) {return fabs(x*x*x);}, 0.5),
    };

  for(unsigned k=0; k<50; k++)
    {
      VariableClenshawCurtis integral = VariableClenshawCurtis();
      integral.set_dim(1);
      integral.set_order(k);
      std::cout << std::endl << "Order is: " << k << std::endl;

      // Calculate the list of integrals using this method and output the error
      for(unsigned i=0; i<function_list.size(); i++)
	{
	  double integral_value = integrate_function_1D(function_list[i].Function, integral);
	  std:: cout << integral_value - function_list[i].Answer << std::endl;
	}
    }
  return 0;
}
