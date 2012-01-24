

// Include infomation on the floating point specifications for our compiler and platform
#include<float.h>

// Generic oomph-lib header
#include<generic.h>

// Header for the type std::function and related
#include<functional>

// Header for wrappers to the LAPACK dgels functions
#include "./dgels.h"

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
  double psi(const unsigned& legendre_order, const int& singularity_order,
	     const double& z, const double& x)
  {
    if(singularity_order == -1)
      return Orthpoly::legendre(legendre_order,x);
    else if(singularity_order == 0)
      return Orthpoly::legendre(legendre_order,x)*log(fabs(z - x)); // log = ln in C++
    else
      return Orthpoly::legendre(legendre_order,x)*(1/( pow((z - x), singularity_order)));
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
  double psi_brandao_finite_part_integral(const unsigned& legendre_order,
					  const int& singularity_order,
					  const double& z, const unsigned& M)
  {
    switch(singularity_order)
      {
      case -1:
	if(legendre_order==0) return 0.0;
	if(legendre_order>0) return 2.0;

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
	  auto integrand = [&legendre_order, &z] (const double& x) -> double
	    { return (Orthpoly::legendre(legendre_order,x)
		      - Orthpoly::legendre(legendre_order,z))/(z-x); };

	  return Orthpoly::legendre(legendre_order,z) * singularity_finite_part(1,z)
	    + gauss_quadrature(M,integrand,-1.0,+1.0);
	}

      case 2:
	{
	  auto integrand = [&legendre_order, &z] (const double& x) -> double
	    { return (Orthpoly::legendre(legendre_order,x)
		      - Orthpoly::legendre(legendre_order,z)
		      + (z-x)*Orthpoly::dlegendre(legendre_order,z))
	      /(z-x); };

	  return Orthpoly::legendre(legendre_order,z) * singularity_finite_part(2,z)
	    - Orthpoly::dlegendre(legendre_order,z) * singularity_finite_part(1,z)
	    + gauss_quadrature(M,integrand,-1.0,+1.0);
	}

	// Else give an error
      default:
	std::ostringstream error_stream;
	error_stream << "Singularity of order " << singularity_order << " is not yet handled.";
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

  // A list of the orders of singularity that we want to calculate the weights for.
  // The method is only known to work up to order 2. Order -1 = constant (no singularity).
  int singularity_orders[] = {-1,1,2};
  unsigned n_sing_orders = sizeof(singularity_orders)/sizeof(int);

  // M is a parameter in the calculation of the method. It is the number of abcissas
  // to use when calculating the moments. Higher M should in theory give a more
  // accurate method.
  const unsigned M = 4;
  const unsigned matrix_height = M*(n_sing_orders);

  // N is a parameter in the method itself. It is the number of abicassas to be
  // used in the finished method.
  const unsigned N = 6;

  // The location of the singularity
  const double z = -1.0;

  ////////////////////////////////////////////////////////////////
  // Assemble matrix and vectors
  ////////////////////////////////////////////////////////////////

  // Find the abcissas (nodes) of the N-point Gauss-Legendre quadrature using oomph-lib
  Vector<double> abcissa_locations(N,0.0);
  Orthpoly::gl_nodes(N,abcissa_locations);

  // Compute the values of psi_i at each node: psi_i = P_i * sing(z - x),
  // x is the location of the singularity, P_i is the ith Gauss-Legendre polynomial.
  DenseMatrix<double> psis(matrix_height,N,234.2);
  // Loop over the N point abicassa
  for(unsigned j=0; j<N; j++)
    {
      // Get the location of this abicassa
      double x = abcissa_locations[j];

      // Loop over the list of psi(x) (singularity orders first then Legendre polynomial orders)
      for(unsigned i=0; i < n_sing_orders; i++)
	  for(unsigned legendre_order=0; legendre_order<M; legendre_order++)
	      psis(i*n_sing_orders+legendre_order,j)
		= psi(legendre_order,singularity_orders[i],z,x);
    }

  // // dump matrix
  // for(unsigned i=0; i<matrix_height; i++)
  //   {
  //     for(unsigned j=0; j<N; j++)
  // 	{
  // 	  std::cout << psis(i,j) << ",";
  // 	}
  //     std::cout << std::endl;
  //   }
  // std::cout << std:: endl;

  // Compute the values of the moments (= the integral of psi_i over [-1,+1], done using Brandao's method)
  // Loop over singularity orders then Legendre polynomial orders.
  Vector<double> moments(matrix_height,0.0);
  for(unsigned i=0; i < n_sing_orders; i++)
    {
      for(unsigned legendre_order=0; legendre_order<M; legendre_order++)
  	{
  	  // Compute the integral of psi_i over [-1,+1]
  	  moments[i*M+legendre_order]
	    = psi_brandao_finite_part_integral(legendre_order,singularity_orders[i],z,M);
  	}
    }

  ////////////////////////////////////////////////////////////////
  // Solve the system using LAPACK
  ////////////////////////////////////////////////////////////////

  // Solve the system psis * weights = moments to get the weights

  //??ds what should rcond be?

  // Declare some inputs and outputs for DEGELSY
  // (from http://www.netlib.org/lapack/double/dgelsy.f).
  int INFO(0), RANK(0), LWORK(0);
  double TEST_WORK[1], RCOND(-1.0);
  int JVPT[N];
  for(unsigned i=0; i<N; i++) {JVPT[i] = 0;} // worst initialisation method ever...
  int LDB = matrix_height > N ? matrix_height : N; // get max(matrix_height,N)

  // We have to convert to array based formats for the fortran solver
  double X[matrix_height], PSIS[matrix_height*N];
  for(unsigned i=0; i<matrix_height; i++)
    {
      X[i] = moments[i];
      for(unsigned j=0; j<N; j++)
	{
	  PSIS[i*N+j] = psis(i,j);
	}
    }

  for(unsigned i=0; i<matrix_height; i++)
    std::cout << X[i] << std::endl;
  std::cout << std::endl;

  // Dump matrix
  for(unsigned i=0; i<matrix_height; i++)
    {
      for(unsigned j=0; j<N; j++)
	{
	  std::cout << PSIS[i*N+j] << ",";
	}
      std::cout << std::endl;
    }
  std::cout << std::endl;

  // Call with LWORK = -1 to just get the optimal value for LWORK
  // (returned value is in WORK[0]). Then create the workspace array.
  DGELSY(matrix_height,N,1,PSIS,matrix_height,X,LDB,JPVT,RCOND,RANK,TEST_WORK,-1,INFO);
  LWORK = int(TEST_WORK[0]);
  double WORK[LWORK];

  // Solve
  DGELSY(matrix_height,N,1,PSIS,matrix_height,X,LDB,
	 JPVT,RCOND,RANK,WORK,LWORK,INFO);

  // The solution is returned in the double array X
  for(unsigned i=0; i<matrix_height; i++)
    {
      std::cout << X[i] << std::endl;
    }

  return 0;
}

