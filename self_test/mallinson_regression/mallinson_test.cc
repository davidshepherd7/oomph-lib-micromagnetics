

#include "../../mallinson_solution.h"
#include "../../magnetic_parameters.h"
#include "../../magnetics_helpers.h"
#include "../../vector_helpers.h"

#include "../../../../src/generic/Vector.h"
#include <iostream>
#include <math.h>

using namespace oomph;
using namespace CompareSolutions;
using namespace VectorOps;

double toradians(const double& degrees)
  {
    return (degrees/180) * 4.0*std::atan(1.0);
  }

Vector<double> dummy_h(const double& t, const Vector<double>& x)
{
  Vector<double> h(3,0.0);
  h[2] = -1.0;
  return h;
}


int main()
  {
    int return_value = 0;

    // parameters
    MagneticParameters params;
    params.Gilbert_damping = 0.5;
    params.Applied_field_fct_pt = &dummy_h;

    // Don't go quite a full arc because it technically takes infinite
    // time. Use same as mallinson in the paper.
    const double theta_start = toradians(10);
    const double theta_end = toradians(170);
    const Vector<double> m0 = sphpolar_to_cart(rpolarazi(1, theta_start, 0.0));

    const unsigned N = 1000;
    const double dtheta = (theta_end - theta_start)/N;

    // Make a list of magnetisations
    Vector<double> thetas(N);
    thetas[0] = theta_start;
    for(unsigned j=1; j<N; j++)
      {
        thetas[j] = thetas[j-1] + dtheta;
      }

    // Get time and phi for each m, print everything.
    for(unsigned j=0; j<thetas.size(); j++)
      {
        const double theta_now = thetas[j];
        const double phi_now = analytic_phi(params.Gilbert_damping, theta_start,
                                            theta_now);

        const double time = switching_time
          (params.Gilbert_damping,
           1.0,
           two_norm(params.h_app(0.0, Vector<double>())),
           params.normalised_hk(),
           theta_start,
           theta_now);

        std::cout << theta_now << " " << time << " " << phi_now
                  << std::endl;

        // Also compare with "m_exact" function
        const Vector<double> m = m_exact(params, m0, time);
        const Vector<double> old_method_m =
          sphpolar_to_cart(rpolarazi(1, theta_now, phi_now));
        if(two_norm_diff(m, old_method_m) > 1e-13)
          {
            std::cerr << "Different values: " << m << old_method_m
                      << "error norm " << two_norm_diff(m, old_method_m)
                      << std::endl;

            return_value++;
          }
      }

    return return_value;
  }
