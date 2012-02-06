
#include "../quadrule.h"

// From quadrule.cc
//
// void clenshaw_curtis_compute ( int n, double x[], double w[] )
//
//  Purpose:
//
//    CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
//
//  Discussion:
//
//    This method uses a direct approach.  The paper by Waldvogel
//    exhibits a more efficient approach using Fourier transforms.
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
//
//    The abscissas for the rule of order ORDER can be regarded
//    as the cosines of equally spaced angles between 180 and 0 degrees:
//
//      X(I) = cos ( ( I - 1 ) * PI / ( ORDER - 1 ) )
//
//    except for the basic case ORDER = 1, when
//
//      X(1) = 0.
//
//    A Clenshaw-Curtis rule that uses ORDER points will integrate
//    exactly all polynomials of degrees 0 through ORDER-1.  If ORDER
//    is odd, then by symmetry the polynomial of degree ORDER will
//    also be integrated exactly.
//
//    If the value of ORDER is increased in a sensible way, then
//    the new set of abscissas will include the old ones.  One such
//    sequence would be ORDER(K) = 2*K+1 for K = 0, 1, 2, ...

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
      clenshaw_curtis_compute(order,knot,weight);

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
      clenshaw_curtis_compute(order,knot,weight);

      // Dump weights to stdout
      std::cout << "{";
      for(unsigned i=0; i<order; i++)
	std::cout << "{" << knot[i] << "},";
      std::cout << "}," << std::endl;
    }

  return 0;
}
