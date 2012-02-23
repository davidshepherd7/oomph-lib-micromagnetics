

// Test chosen integration scheme on a variety of functions with known analytic values.
// C++0x only (used lambda functions to simplify defining the long list of testing functions).

#include "../../variable_quadrature.h"
#include <functional>
#include <cmath>

using namespace oomph;

namespace oomph
{

  // Integrate a given function using the given scheme (no check is made that
  // the scheme integrates between -1 and +1 - watch out).
  double integrate_function_1D(const std::function<double(double)> &function,
			       const VariableOrderQuadrature &integral, const unsigned &order)
  {
    double integral_value = 0;

    // Loop over integration points
    unsigned n_weights = integral.nweight(order);
    for(unsigned i=0; i<n_weights; i++)
      {
	// Get knot location
	double x = integral.knot(i,0,order);

	// Add contribution from this knot
	integral_value += function(x) * integral.weight(i,order);
      }

    return integral_value;
  }

  // Structure to store the information about each testing function
  struct TestingFn
  {
    // Stored values
    std::function<double(double)> Function;
    double Answer;
    std::string Label;

    // Constructor
    TestingFn(const std::function<double(double)> &function,
	      const double &answer, const std::string label)
    {
      Function = function;
      Answer = answer;
      Label = label;
    }
  };

}

int main()
{

  // Create a list of testing functions and their exact solutions Exact
  // solutions are taken from Mathematica using the command:
  // N[Integrate[...fn..., {x, -1, +1}], 20]
  std::vector<TestingFn > function_list =
    {
      // Some very simple functions:
      TestingFn([] (double x) {return 1.0;}, 2.0, "1.0"),
      TestingFn([] (double x) {return x;}, 0.0, "x"),

      // Functions from Trefethen2008:
      TestingFn([] (double x) {return pow(x,20);}, 2.0/21, "x^20"),
      TestingFn([] (double x) {return exp(x);},
		2.3504023872876029138, "exp(x)"),
      TestingFn([] (double x) {return exp(-pow(x,2));},
		1.4936482656248540508, "exp(-x^2)"),
      TestingFn([] (double x) {return 1/(1 + 16*x*x);},
		0.66290883183401623253, "1/(1 + 16 x^2)"),
      TestingFn([] (double x) {return exp(-1/(x*x));},
		0.17814771178156069019, "exp(-1/(x^2))"),
      TestingFn([] (double x) {return fabs(x*x*x);}, 0.5, "|x^3|"),

      // Functions related to the hybrid FEM/BEM for micromagnetics
      // =================================================================

      // 1/r
      TestingFn([] (double x) {return 1/fabs(x-2);},
		1.0986122886681096914, "1/|x-2|"),
      TestingFn([] (double x) {return 1/fabs(x-5);},
		0.40546510810816438198, "1/|x-5|"),
      TestingFn([] (double x) {return 1/fabs(x-20);},
		0.10008345855698253649, "1/|x-20|"),

      // 1/(r^2)
      TestingFn([] (double x) {return 1/pow(fabs(x-2),2);},
		0.6666666666666666666, "1/|x-2|^2"),
      TestingFn([] (double x) {return 1/pow(fabs(x-5),2);},
		0.083333333333333333333, "1/|x-5|^2"),
      TestingFn([] (double x) {return 1/pow(fabs(x-20),2);},
		0.0050125313283208020050, "1/|x-20|^2"),
    };

  // Create a list of integration methods to test:
  Vector<VariableOrderQuadrature*> integration_methods(3);
  integration_methods[0] = new VariableOrderGaussLegendre();
  integration_methods[1] = new VariableOrderFejerSecond();
  integration_methods[2] = new VariableOrderClenshawCurtis();

  // The highest order available:
  unsigned max_order = 50;

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

      // Output list of function labels as first line of file
      outfile << "\"order\" ";
      for(unsigned i=0; i<function_list.size(); i++)
	outfile << "\"" << function_list[i].Label << "\" ";
      outfile << std::endl;

      // Loop over orders of the integration scheme
      for(unsigned k=0; k<=max_order; k++)
	{
	  outfile << k << " ";

	  // Calculate each of th integrals using this method and output the error
	  for(unsigned i=0; i<function_list.size(); i++)
	    {
	      // Do the quadrature
	      double integral_value = integrate_function_1D
		(function_list[i].Function, *integration_methods[int_meth],k);

	      // Output the error
	      outfile << " " << fabs(integral_value - function_list[i].Answer);
	    }

	  outfile << std::endl;
	}

      // Done with this integration method so close the file
      outfile.close();
    }


  // Clean up the integration methods vector
  for(unsigned i=0; i<integration_methods.size(); i++)
    {
      delete integration_methods[i];
    }

  return 0;
}
