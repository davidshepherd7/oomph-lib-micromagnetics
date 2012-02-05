
# include "./quadrule.hpp"

// Details of the function used (from quadrule.cpp):
//
//  Purpose:
//
//    LEGENDRE_DR_COMPUTE: Gauss-Legendre quadrature by Davis-Rabinowitz method.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//   ............
//
//  Parameters:
//
//    Input, int ORDER, the order.
//    ORDER must be greater than 0.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//

int main()
{

  // Set output precision to 16 digits
  std::cout.precision(16);
  std::cout.setf(std::ios::fixed,std::ios::floatfield);


  unsigned max_order = 50;

  // Calculate and output the weights
  for(unsigned order=1; order<=max_order; order++)
    {
      // Calculate
      double weight[order], knot[order];
      legendre_dr_compute (order,knot,weight);

      // Dump weights to stdout
      std::cout << "{";
      for(unsigned i=0; i<order; i++)
	std::cout << weight[i] << ",";
      std::cout << "}," << std::endl;
    }

  std::cout << std::endl;

  // Calcualte and output the knots (redundant calculations but easier to
  // output this way)
  for(unsigned order=1; order<=max_order; order++)
    {
      // Calculate
      double weight[order], knot[order];
      legendre_dr_compute (order,knot,weight);

      // Dump weights to stdout
      std::cout << "{";
      for(unsigned i=0; i<order; i++)
	std::cout << "{" << knot[i] << "},";
      std::cout << "}," << std::endl;
    }

  return 0;
}
