

// Test chosen integration scheme on a variety of functions with known analytic values.

#include "../../variable_order_quadrature.h"
#include <functional>
#include <cmath>

using namespace oomph;

namespace oomph
{

  // Function typedef
  typedef double (*fndoubletodouble)(const double&);

  // Integrate a given function, f(x), using the given scheme (no check is made
  // that the scheme integrates between -1 and +1 - watch out).  If dim = 2
  // integrate f(x,y) = f(x) * f(y), dim = 3 f(x,y,z) = f(x)*f(y)*f(z).
  double integrate_function(const fndoubletodouble &function,
			    const BaseVariableOrderQuadrature &integral,
			    const unsigned &order,
			    const unsigned &dim)
  {
    double integral_value = 0;

    // Loop over integration points
    unsigned n_weights = integral.nweight(order);
    for(unsigned i=0; i<n_weights; i++)
      {

	// Get knot location
	Vector<double> x_kn(dim,0.0);
	for(unsigned j=0; j<dim; j++)
	  {
	    x_kn[j] = integral.knot(i,j,order);
	  }

	// Calculate the value of the product of the function over all
	// dimensions at this knot (i.e. f(x,y,z) = f(x) * f(y) * f(z))
	double function_value = 1.0;
	for(unsigned j=0; j<dim; j++)
	  {
	    function_value *= function(x_kn[j]);
	  }

	integral_value += function_value * integral.weight(i,order);
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
  double xtothe20(const double &x) {return std::pow(x,20);}
  double xtothe30(const double &x) {return std::pow(x,30);}
  double xtothe64(const double &x) {return std::pow(x,64);}
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

  // Mononomials
  function_list.push_back(TestingFn(one, 2.0, "1.0"));
  function_list.push_back(TestingFn(justx, 0.0, "x"));
  function_list.push_back(TestingFn(xtothe20, 2.0/21, "x^20"));
  function_list.push_back(TestingFn(xtothe30, 2.0/31, "x^30"));
  function_list.push_back(TestingFn(xtothe64, 2.0/65, "x^64"));

  // function_list.push_back(TestingFn(expx, 2.3504023872876029138, "exp(x)"));
  // function_list.push_back(TestingFn(expxmsq, 1.4936482656248540508, "exp(-x^2)"));
  // function_list.push_back(TestingFn(oneovpoly, 0.66290883183401623253, "1/(1 + 16 x^2)"));
  // function_list.push_back(TestingFn(exp1ovxsq, 0.17814771178156069019, "exp(-1/(x^2))"));
  // function_list.push_back(TestingFn(absxcubed, 0.5, "|x^3|"));

  // // 1/r near singularities
  // function_list.push_back(TestingFn(sing2, 1.0986122886681096914, "1/|x-2|"));
  // function_list.push_back(TestingFn(sing5, 0.40546510810816438198, "1/|x-5|"));
  // function_list.push_back(TestingFn(sing20, 0.10008345855698253649, "1/|x-20|"));

  // // 1/(r^2) near singularities
  // function_list.push_back(TestingFn(sqsing2, 0.6666666666666666666, "1/|x-2|^2"));
  // function_list.push_back(TestingFn(sqsing5, 0.083333333333333333333, "1/|x-5|^2"));
  // function_list.push_back(TestingFn(sqsing20, 0.0050125313283208020050, "1/|x-20|^2"));



  // Create a list of integration methods to test:
  unsigned n_int_schemes = 3;
  Vector<BaseVariableOrderQuadrature*> integration_methods;
  integration_methods.push_back(new QVariableOrderGaussLegendre<1>());
  integration_methods.push_back(new QVariableOrderFejerSecond<1>());
  integration_methods.push_back(new QVariableOrderClenshawCurtis<1>());

  integration_methods.push_back(new QVariableOrderGaussLegendre<2>());
  integration_methods.push_back(new QVariableOrderFejerSecond<2>());
  integration_methods.push_back(new QVariableOrderClenshawCurtis<2>());

  integration_methods.push_back(new QVariableOrderGaussLegendre<3>());
  integration_methods.push_back(new QVariableOrderFejerSecond<3>());
  integration_methods.push_back(new QVariableOrderClenshawCurtis<3>());


  // The list of orders to test:
  std::vector<unsigned> order_array
    = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,
       28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,
       64,128,256,};

  unsigned order_array_length = order_array.size();


  // Loop over dimensions
  for(unsigned j=0; j<3; j++)
    {
      unsigned dim = j+1;

      // Loop over the integration methods
      for(unsigned int_meth=0; int_meth<n_int_schemes; int_meth++)
	{
	  // Set up output file
	  std::ofstream outfile;
	  outfile.precision(16);
	  char filename[100]; sprintf(filename,"results/integration_method_%u_dim_%u",int_meth,dim);
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
		  double integral_value
		    = integrate_function(function_list[i].Function,
					 *integration_methods[int_meth + n_int_schemes*j],
					 order,dim);

		  // Output the error
		  outfile << " " << fabs(integral_value
					 - std::pow(function_list[i].Answer,dim));
		}

	      outfile << std::endl;
	    }

	  // Done with this integration method so close the file
	  outfile.close();
	}

       }

  // Clean up the integration methods vector
  for(unsigned i=0; i<integration_methods.size(); i++)
    {
      delete integration_methods[i];
    }

  return 0;
}
