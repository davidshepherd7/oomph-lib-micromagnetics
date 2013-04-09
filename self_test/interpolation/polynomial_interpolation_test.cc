



#include "generic.h"
#include <fenv.h>

#include "../../poly_interp.h"



using namespace oomph;


namespace input
{
  void vec_function(const double &x, Vector<double>& ys)
  {
    ys.assign(2, 0.0);

    ys[0] = x*x;
    ys[1] = cos(x);
  }

  void d_vec_function(const double &x, Vector<double>& ys)
  {
    ys.assign(2, 0.0);

    ys[0] = 2*x;
    ys[1] = -sin(x);
  }

}


int main()
{
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  // Set up input values (C++98 SUCKS!)
  double a[] = {1.573, 1.58, 1.6, 1.7, 1.8, 1.9, 2.1, 2.14};
  unsigned n = 8;
  Vector<double> xs;
  for(unsigned i=0; i<n; i++)
    {
      xs.push_back(a[i]);
    }
  Vector<Vector<double> > ys;
  for(unsigned i=0, ni=xs.size(); i<ni; i++)
    {
      Vector<double> temp;
      input::vec_function(xs[i], temp);
      ys.push_back(temp);
    }

  // Test interpolation at a point
  BarycentricLagrangeInterpolator b(xs, ys);
  double eval_point = 1.74;

  Vector<double> output, exact;
  b.eval(eval_point, output);
  input::vec_function(eval_point, exact);

  // Check that they are the same
  OOMPH_ASSERT(almost_equal(output, exact, 1e-10));

  Vector<double> output_deriv, exact_deriv;
  b.eval_derivative(eval_point, 1, output_deriv);
  input::d_vec_function(eval_point, exact_deriv);
  OOMPH_ASSERT(almost_equal(output_deriv, exact_deriv, 1e-10));

  std::cout << output << exact << std::endl;

  std::cout << output_deriv << exact_deriv << std::endl;


  return 0;
}
