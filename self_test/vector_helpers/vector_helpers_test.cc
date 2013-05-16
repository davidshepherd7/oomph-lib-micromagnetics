

#include "generic.h"
#include "../../vector_helpers.h"
#include <fenv.h>

using namespace VectorOps;

int main()
{
  // Start with return value of zero, if nothing fails it will stay zero.
  int return_value = 0;

  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);



  // angle diff function
  // ============================================================
  {
  double pi = 4*std::atan2(1,1);

  Vector<Vector<double> > v1s, v2s;
  Vector<double> answers;

  // perp, 2d
  {
    Vector<double> a(2), b(2);
    a[0] = 0; a[1] = 1;
    b[0] = 1; b[1] = 0;
    double c = pi/2;
    v1s.push_back(a); v2s.push_back(b); answers.push_back(c);
  }

  // parallel 2d, with one much longer vector
  {
    Vector<double> a(2), b(2);
    a[0] = 10; a[1] = 0;
    b[0] = 1; b[1] = 0;
    double c = 0;
    v1s.push_back(a); v2s.push_back(b); answers.push_back(c);
  }

  // perp 3d
  {
    Vector<double> a(3), b(3);
    a[0] = 1; a[1] = 0; a[2] = 0;
    b[0] = 0; b[1] = 0; b[2] = 1;
    double c = pi/2;
    v1s.push_back(a); v2s.push_back(b); answers.push_back(c);
  }

  // parallel 3d
  {
    Vector<double> a(3), b(3);
    a[0] = 1; a[1] = 0; a[2] = 0;
    b[0] = 1; b[1] = 0; b[2] = 0;
    double c = 0;
    v1s.push_back(a); v2s.push_back(b); answers.push_back(c);
  }

  {
    Vector<double> a(3), b(3);
    a[0] = 1; a[1] = 0; a[2] = 0;
    b[0] = 1; b[1] = 1; b[2] = 0;
    double c = pi/4;
    v1s.push_back(a); v2s.push_back(b); answers.push_back(c);
  }


  // Run tests
  for(unsigned i=0; i<answers.size(); i++)
    {
      double result = angle_diff(v1s[i], v2s[i]);
      if(!numerical_zero(result - answers[i]))
        {
          std::cerr << "Error:"
                    << "angle_diff("
                    << v1s[i] << ", " << v2s[i] << ") == " << answers[i] <<std::endl
                    << "Got: " << result << std::endl;
          return_value++;
        }
    }

  }

  // Optimised cross product
  // ============================================================
  {
    std::cout << "Checking Optimised cross product." << std::endl;
    Vector<Vector<double> > a, b;

    // Do for some random vectors
    for(unsigned i=0; i<30; i++)
      {
        Vector<double> tempa = random_vector(3);
        Vector<double> tempb = random_vector(3);

        a.push_back(tempa);
        b.push_back(tempb);
      }

    // Run tests
    for(unsigned i=0; i< a.size(); i++)
      {
        Vector<double> v(3);
        v[0] = opt_cross(0, a[i], b[i]);
        v[1] = opt_cross(1, a[i], b[i]);
        v[2] = opt_cross(2, a[i], b[i]);

        // Just compare with the old cross product implementation
        double error = two_norm_diff(v, cross(a[i],b[i]));

        if(!numerical_zero(error))
          {
            std::cerr << "Error:"
                      << "opt_cross("
                      << a[i] << ", " << b[i] << ") == " << cross(a[i], b[i]) <<std::endl
                      << "Got: " << v << std::endl;
            return_value++;
          }
      }
  }


  return return_value;
}
