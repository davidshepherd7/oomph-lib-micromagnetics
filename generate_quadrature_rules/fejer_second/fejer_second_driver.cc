
#include "../quadrule.h"

// From quadrule.cc:
//
// void fejer2_compute ( int n, double x[], double w[] )
//
//  Purpose:
//
//    FEJER2_COMPUTE computes a Fejer type 2 quadrature rule.
//
//  Discussion:
//
//    This method uses a direct approach.  The paper by Waldvogel
//    exhibits a more efficient approach using Fourier transforms.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
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
      fejer2_compute(order,knot,weight);

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
      fejer2_compute(order,knot,weight);

      // Dump weights to stdout
      std::cout << "{";
      for(unsigned i=0; i<order; i++)
	std::cout << "{" << knot[i] << "},";
      std::cout << "}," << std::endl;
    }

  return 0;
}
