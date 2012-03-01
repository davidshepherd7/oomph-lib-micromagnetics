
// Test integration schemes against polynomials of appropriate degree.

#include "../../variable_order_quadrature.h"
#include <functional>

using namespace oomph;

namespace oomph
{

  // Function typedef
  typedef std::function<double(double)> fndoubletodouble;

  struct int_scheme_data
  {
  public:
    BaseVariableOrderQuadrature* Scheme;
    unsigned Degree;
    unsigned Dim;

    int_scheme_data(BaseVariableOrderQuadrature* scheme,
		    const unsigned &degree,
		    const unsigned &dim)
      : Scheme(scheme), Degree(degree), Dim(dim)
    {}
  };

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

  std::pair<fndoubletodouble,double> create_poly(const int &degree)
  {
    fndoubletodouble function = [degree] (const double &x)
      -> double {return std::pow(x,degree);};

    double exact_solution = (1 - std::pow(-1,degree+1) )/(degree + 1);

    // std::cout << degree << " " << exact_solution << std::endl;

    return std::pair<fndoubletodouble,double>(function,exact_solution);
  }

  unsigned scheme_degree(const unsigned &order)
  {
    // // Clenshaw/Fejer
    // return order-1;

    // Gauss-Legendre
    return 2*order+1;
  }

} // End of namespace oomph
int main()
{

  std::vector<unsigned> order_list
    = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,
       28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,
       64,128,256,};

  std::vector<BaseVariableOrderQuadrature*> integration_scheme_list;
  integration_scheme_list.push_back(new QVariableOrderGaussLegendre<1>());
  integration_scheme_list.push_back(new QVariableOrderGaussLegendre<2>());
  integration_scheme_list.push_back(new QVariableOrderGaussLegendre<3>());

  // std::vector<BaseVariableOrderQuadrature*> integration_scheme_list;
  // integration_scheme_list.push_back(new QVariableOrderClenshawCurtis<1>());
  // integration_scheme_list.push_back(new QVariableOrderClenshawCurtis<2>());
  // integration_scheme_list.push_back(new QVariableOrderClenshawCurtis<3>());

  // std::vector<BaseVariableOrderQuadrature*> integration_scheme_list;
  // integration_scheme_list.push_back(new QVariableOrderFejerSecond<1>());
  // integration_scheme_list.push_back(new QVariableOrderFejerSecond<2>());
  // integration_scheme_list.push_back(new QVariableOrderFejerSecond<3>());

  std::cout << "\"order\" \"dim\" \"error: exact poly\" \"error: inexact poly\""
	    << std::endl;

  // For each order in the list
  for(unsigned i_order=0; i_order<order_list.size(); i_order++)
    {
      // Get the order
      unsigned order = order_list[i_order];

      // Construct a polynomial of degree exactly integrable by scheme
      std::pair<fndoubletodouble,double> poly_exact = create_poly(scheme_degree(order));

      // Construct a polynomial of degree not exactly integrable by scheme
      std::pair<fndoubletodouble,double> poly_inexact = create_poly(scheme_degree(order)+1);

      // For each dimension
      for(unsigned dim=1; dim<=3; dim++)
	{
	  // Get the scheme for this dimension
	  BaseVariableOrderQuadrature* integration_scheme = integration_scheme_list[dim-1];

	  // Integrate them using the scheme
	  double solution_exact = integrate_function
	    (poly_exact.first,*integration_scheme,order,dim);
	  double solution_inexact = integrate_function
	    (poly_inexact.first,*integration_scheme,order,dim);

	  // Compare results
	  double error_exact = std::abs(solution_exact - std::pow(poly_exact.second,dim));
	  double error_inexact = std::abs(solution_inexact - std::pow(poly_inexact.second,dim));

	  std::cout << order << " " << dim << " " << error_exact << " "
		    << error_inexact << std::endl;

	  // Error if exactly integrable polynomial is not extremely accurate
	  if (error_exact > 1e-15)
	    {
	      std::ostringstream err_stream;
	      err_stream << "Scheme of order " << order << " and dimension "
	      		 << dim << " gave an error of " << error_exact
	      		 << " for a polynomial that should be exactly integrable.";
	      throw OomphLibWarning(err_stream.str(),
	      			    "rigorous_degree_test::main",
	      			    OOMPH_EXCEPTION_LOCATION);
	    }

	}
    }

  return 1;
}
