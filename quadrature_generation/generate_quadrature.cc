
# include "generic/quadrule.h"

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


void dump_array(double array[],unsigned length)
{
  std::cout << "{" << array[0];
  for(unsigned i=1; i<length; i++)
    std::cout << "," << array[i];
  std::cout << "}";
}

int main()
{

  unsigned max_order = 50;

  for(unsigned order=2; order<=max_order; order++)
    {
      // Calculate
      double weight[order], knot[order];
      legendre_dr_compute (order,knot,weight);

      // Dump weights to stdout
      std::cout << "Weights[" << order - 2 << "] = {" << weight[0];
      for(unsigned i=1; i<order; i++)
	std::cout << "," << weight[i];
      std::cout << "};" << std::endl;
    }

  for(unsigned order=2; order<=max_order; order++)
    {
      // Calculate
      double weight[order], knot[order];
      legendre_dr_compute (order,knot,weight);

      // Dump weights to stdout
      std::cout << "Knots[" << order - 2 << "] = {{" << knot[0] << "}";
      for(unsigned i=1; i<order; i++)
	std::cout << ",{" << knot[i] << "}";
      std::cout << "};" << std::endl;
    }

  return 0;
}