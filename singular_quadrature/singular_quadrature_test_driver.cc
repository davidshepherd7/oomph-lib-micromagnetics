

// Include infomation on the floating point specifications for our compiler and platform
#include<float.h>

#include<generic.h>
#include<functional>

// Include oomph-lib library for calculation of nodes of the Gauss Legendre polynomials
#include "generic/orthpoly.h"

using namespace oomph;

namespace oomph
{
  // We don't have the ability to use long doubles in the lower level oomph-lib calls
  // so this is useless :(
  // // Check the accuracy of long doubles vs doubles since the implementation is not
  // // standardised.
  // void check_long_double_epsilon()
  // {
  //   // Check what accuracy we have with
  //   std::cout << "Epsilon of a long double on your compiler, machine combination is: "
  // 	      << LDBL_EPSILON << std::endl;
  //   std::cout << "Epsilon of a double on your compiler, machine combination is: "
  // 	      << DBL_EPSILON << std::endl;

  //   if (!(LDBL_EPSILON/DBL_EPSILON < 0.01))
  //     {
  // 	std::cerr << "\nWarning: long doubles may not be sufficiently accurate to get the maximum accuracy.\nTo improve accuracy use a compiler which allows for more accurate long doubles (e.g. gnu, intel)." << std::endl;
  //     }
  // }


  /// Calculate psi_i = P_i * singularity. Singularity order -1 gives P_i, singularity order 0 gives P_i*ln|x-t|, higher orders give 1/((t-x)^sing_order).
    double psi(const unsigned& i, const int& singularity_order,
	     const double& z, const double& x)
    {
    if(singularity_order == -1)
      return Orthpoly::legendre(i,x);
    else if(singularity_order == 0)
      return Orthpoly::legendre(i,x)*log(fabs(z - x)); // log = ln in C++
    else
      return Orthpoly::legendre(i,x)*(1/( pow((z - x), singularity_order)));
  }

  /// Get the finite part integral of a singularity of order one or two at z.
  double singularity_finite_part(const int& singularity_order, const double& z)
  {
    if(singularity_order == 1)
      {
	if(fabs(z) == 1) return z*log(2);
	else return log(fabs( (z+1)/(z-1) ));
      }
    else if(singularity_order == 2)
      {
	if(fabs(z) == 1) return -0.5;
	else return 2/(z*z - 1);
      }
    else
      {
	std::cerr << "Finite part of singularity of order " << singularity_order
		  << " is not handled." << std::endl;
	throw;
      }
  }

  /// Compute the M point Gaussian quadrature of the function fct_pt. One dimensional for now, should probably extend it in the future.
  double gauss_quadrature(const unsigned& M,
			  const std::function<double(const double&)>& integrand,
			  const double& lower_limit, const double& upper_limit)
  {
    // Get the M locations and weights of the abcissas
    Vector<double> locations(M,0.0), weights(M,0.0);
    Orthpoly::gl_nodes(M,locations,weights);

    // Sum the weight*integrand value over the M abcissas
    double approximate_integral = 0.0;
    for(unsigned abcissa=0; abcissa<M; abcissa++)
      approximate_integral += integrand(locations[abcissa]) * weights[abcissa];

    return approximate_integral;
  }

  /// Calculate the integral of psi on [-1,+1] using Brandao's finite part method. The input i is the order of the Legrende polynomial, z is the location of the singularity and M is the order of Gaussian Quadrature to use in the evaluation of the modified integral.
  double psi_brandao_finite_part_integral(unsigned& i, int& singularity_order,
					  double& z, unsigned& M)
  {
    switch(singularity_order)
      {
      case -1:
	if(i==0) return 0.0;
	if(i>0) return 2.0;

      case 0:
	if(fabs(z)!=1.0) return 0; //return 2*(Q(i+1,z) - Q(i-1,z))/(2*i +1);
	else
	  {
	    //??ds do this eventually :(
	    std::cerr << "ln calculation with z = +-1 is horrible so not done yet."
		      << std::endl;
	  }

      case 1:
	{
	  auto integrand = [&i, &z] (const double& x) -> double
	    { return (Orthpoly::legendre(i,x) - Orthpoly::legendre(i,z))/(z-x); };

	  return Orthpoly::legendre(i,z) * singularity_finite_part(1,z)
	    + gauss_quadrature(M,integrand,-1.0,+1.0);
	}

      case 2:
	{
	  auto integrand = [&i, &z] (const double& x) -> double
	    { return (Orthpoly::legendre(i,x) - Orthpoly::legendre(i,z)
		      + (z-x)*Orthpoly::dlegendre(i,z))
	      /(z-x); };

	  return Orthpoly::legendre(i,z) * singularity_finite_part(2,z)
	    - Orthpoly::dlegendre(i,z) * singularity_finite_part(1,z)
	    + gauss_quadrature(M,integrand,-1.0,+1.0);
	}

	// Else give an error
      default:
	std::ostringstream error_stream;
	error_stream << "Singularity of order " << i << " is not yet handled.";
	throw OomphLibError(error_stream.str(),
			    "singular_quadrature::psi_brandao_finite_part_integral",
			    OOMPH_EXCEPTION_LOCATION);
      }

  }

} // End of namespace oomph

  // Compute the quadrature and dump to a file
int main(int argc, char* argv[])
{
  ////////////////////////////////////////////////////////////////
  // Set up
  ////////////////////////////////////////////////////////////////

  // The maximum order of the singularity that this method is known to work with.
  int max_singularity = 2;
  int min_singularity = 1; // singularity order -1 = constant (no singularity)

  // M is a parameter in the calculation of the method. It is the number of abcissas
  // to use when calculating the moments. Higher M should in theory give a more
  // accurate method.
  unsigned M = 4;
  unsigned matrix_height = M*(max_singularity+2);

  // N is a parameter in the method itself. It is the number of abicassas to be
  // used in the finished method.
  unsigned N = 16;

  // The location of the singularity
  double z = -1.0;

  ////////////////////////////////////////////////////////////////
  // Assemble matrix and vectors
  ////////////////////////////////////////////////////////////////

  // Find the abcissas (nodes) of the N-point Gauss-Legendre quadrature using oomph-lib
  Vector<double> abcissa_locations(N,0.0);
  Orthpoly::gl_nodes(N,abcissa_locations);

  // Compute the values of psi_i at each node: psi_i = P_i * 1/(x - t_j),
  // x is the location of the singularity, P_i is the ith Gauss-Legendre polynomial.
  DenseDoubleMatrix psis(matrix_height,N,0.0);
  // Loop over the N point abicassa
  for(unsigned j=0; j<N; j++)
    {
      // Get the location of this abicassa
      double x = abcissa_locations[j];

      // Loop over the list of psi (singularity orders first then Legendre polynomial orders)
      for(int singularity_order=min_singularity;
	  singularity_order<=max_singularity; singularity_order++)
	{
	  for(unsigned i=0; i<M; i++)
	    {
	      // Compute psi_i_so at the abcissa
	      psis(i*(singularity_order+2),j) = psi(i,singularity_order,z,x);
	      std::cout << psis(i*(singularity_order+2),j) << ",";
	    }
	}
      std::cout << std::endl;
    }

  // Compute the values of the moments (= the integral of psi_i over [-1,+1], done using Brandao's method)
  // Loop over singularity orders then Legendre polynomial orders.
  Vector<double> moments(matrix_height,0.0);
  for(int singularity_order=min_singularity;
      singularity_order<=max_singularity; singularity_order++)
    {
      for(unsigned i=0; i<M; i++)
  	{
  	  // Compute the integral of psi_i over [-1,+1]
  	  moments[i] = psi_brandao_finite_part_integral(i,singularity_order,z,M);
  	  //std::cout << moments[i] << std::endl;
  	}
    }

  ////////////////////////////////////////////////////////////////
  // Solve the system
  ////////////////////////////////////////////////////////////////

  // Solve the system psis * weights = moments to get the weights
  Vector<double> weights(N,0.0);
  psis.linear_solver_pt()->solve(&psis,moments,weights);

  return 0;
}

