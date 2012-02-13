

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
  // Exact solutions are taken from Mathematica using the command:
  // N[Integrate[...fn..., {x, -1, +1}], 20]
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

      // Functions related to the hybrid FEM/BEM for micromagnetics
      // =================================================================

      // 1/r
      TestingFn([] (double x) {return 1/fabs(x-2);}, 1.0986122886681096914),
      TestingFn([] (double x) {return 1/fabs(x-5);}, 0.40546510810816438198),
      TestingFn([] (double x) {return 1/fabs(x-20);}, 0.10008345855698253649),

      // 1/(r^2)
      TestingFn([] (double x) {return 1/pow(fabs(x-2),2);}, 0.66666666666666666667),
      TestingFn([] (double x) {return 1/pow(fabs(x-5),2);}, 0.083333333333333333333),
      TestingFn([] (double x) {return 1/pow(fabs(x-20),2);}, 0.0050125313283208020050),
    };

  // Create a list of integration methods to test:
  Vector<VariableQuadrature*> integration_methods(3);
  integration_methods[0] = new VariableGaussLegendre();
  integration_methods[1] = new VariableFejerSecond();
  integration_methods[2] = new VariableClenshawCurtis();

  // Loop over the integration methods
  for(unsigned int_meth=0; int_meth<integration_methods.size(); int_meth++)
    {
      integration_methods[int_meth]->set_dim(1);

      // Set up output file
      std::ofstream outfile;
      outfile.setf(std::ios::fixed,std::ios::floatfield);
      outfile.precision(16);
      char filename[100]; sprintf(filename,"results/integration_method_%u",int_meth);
      outfile.open(filename);

      // Loop over orders of the integration scheme
      for(unsigned k=0; k<50; k++)
	{
	  integration_methods[int_meth]->set_order(k);

	  outfile << k << " ";

	  // Calculate each of th integrals using this method and output the error
	  for(unsigned i=8; i<function_list.size(); i++)
	    {
	      double integral_value =
		integrate_function_1D(function_list[i].Function, *integration_methods[int_meth]);
	      outfile << " " << fabs(integral_value - function_list[i].Answer);
	    }

	  outfile << std::endl;
	}

      outfile.close();
      outfile << std::endl;
    }

  // Clean up the integration methods vector
  for(unsigned i=0; i<integration_methods.size(); i++)
    {
      delete integration_methods[i];
    }

  return 0;
}
