

// Test chosen integration scheme on a variety of functions with known analytic values.

#include "../../variable_order_quadrature.h"
#include <functional>
#include <cmath>

using namespace oomph;

namespace oomph
{

  // Function typedef
  typedef double (*fndoubletodouble)(const double&);

  // Integrate a given function using the given scheme (no check is made that
  // the scheme integrates between -1 and +1 - watch out).
  double integrate_function_1D(const fndoubletodouble &function,
			       const QVariableOrderQuadrature<1> &integral,
			       const unsigned &order)
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
    fndoubletodouble Function;
    double Answer;
    std::string Label;

    // Constructor
    TestingFn(const fndoubletodouble function,
	      const double &answer, const std::string label)
    {
      Function = function;
      Answer = answer;
      Label = label;
    }
  };

  // Some analytically known functions for testing:
  //============================================================

  double one(const double &x) {return 1.0;}
  double justx(const double &x) {return x;}

  // Functions from Trefethen2008:
  double xtothe20(const double &x) {return pow(x,20);}
  double expx(const double &x) {return exp(x);}
  double expxmsq(const double &x) {return exp(-pow(x,2));}
  double oneovpoly(const double &x) {return 1/(1 + 16*x*x);}
  double exp1ovxsq(const double &x) {return exp(-1/(x*x));}
  double absxcubed(const double &x) {return fabs(x*x*x);}

  // near 1/r singularity
  double sing2(const double &x) {return 1/fabs(x-2);}
  double sing5(const double &x) {return 1/fabs(x-5);}
  double sing20(const double &x) {return 1/fabs(x-20);}

  // near 1/(r^2) singularity
  double sqsing2(const double &x) {return 1/pow(fabs(x-2),2);}
  double sqsing5(const double &x) {return 1/pow(fabs(x-5),2);}
  double sqsing20(const double &x) {return 1/pow(fabs(x-20),2);}

}

int main()
{

  // Create a list of testing functions, their exact solutions and a
  // label. Exact solutions are taken from Mathematica using the command:
  // N[Integrate[...fn..., {x, -1, +1}], 20]
  std::vector<TestingFn > function_list;

  function_list.push_back(TestingFn(one, 2.0, "1.0"));
  function_list.push_back(TestingFn(justx, 0.0, "x"));

  // Functions from Trefethen2008:
  function_list.push_back(TestingFn(xtothe20, 2.0/21, "x^20"));
  function_list.push_back(TestingFn(expx, 2.3504023872876029138, "exp(x)"));
  function_list.push_back(TestingFn(expxmsq, 1.4936482656248540508, "exp(-x^2)"));
  function_list.push_back(TestingFn(oneovpoly, 0.66290883183401623253, "1/(1 + 16 x^2)"));
  function_list.push_back(TestingFn(exp1ovxsq, 0.17814771178156069019, "exp(-1/(x^2))"));
  function_list.push_back(TestingFn(absxcubed, 0.5, "|x^3|"));

  // 1/r near singularities
  function_list.push_back(TestingFn(sing2, 1.0986122886681096914, "1/|x-2|"));
  function_list.push_back(TestingFn(sing5, 0.40546510810816438198, "1/|x-5|"));
  function_list.push_back(TestingFn(sing20, 0.10008345855698253649, "1/|x-20|"));

  // 1/(r^2) near singularities
  function_list.push_back(TestingFn(sqsing2, 0.6666666666666666666, "1/|x-2|^2"));
  function_list.push_back(TestingFn(sqsing5, 0.083333333333333333333, "1/|x-5|^2"));
  function_list.push_back(TestingFn(sqsing20, 0.0050125313283208020050, "1/|x-20|^2"));



  // Create a list of integration methods to test:
  Vector<QVariableOrderQuadrature<1>*> integration_methods(3);
  integration_methods[0] = new QVariableOrderGaussLegendre<1>();
  integration_methods[1] = new QVariableOrderFejerSecond<1>();
  integration_methods[2] = new QVariableOrderClenshawCurtis<1>();

  // The lsit of orders to test:
  unsigned order_array_length = 53;
  unsigned order_array[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,64,128,256,};

  // Loop over the integration methods
  for(unsigned int_meth=0; int_meth<integration_methods.size(); int_meth++)
    {
      // Set up output file
      std::ofstream outfile;
      outfile.precision(16);
      char filename[100]; sprintf(filename,"results/integration_method_%u",int_meth);
      outfile.open(filename);

      // Output list of function labels as first line of file
      outfile << "\"order\" ";
      for(unsigned i=0; i<function_list.size(); i++)
	outfile << "\"" << function_list[i].Label << "\" ";
      outfile << std::endl;

      // Loop over orders of the integration scheme
      for(unsigned k=0; k<order_array_length; k++)
	{
	  unsigned order = order_array[k];

	  outfile << order << " ";

	  // Calculate each of th integrals using this method and output the error
	  for(unsigned i=0; i<function_list.size(); i++)
	    {
	      // Do the quadrature
	      double integral_value = integrate_function_1D
		(function_list[i].Function, *integration_methods[int_meth],order);

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
