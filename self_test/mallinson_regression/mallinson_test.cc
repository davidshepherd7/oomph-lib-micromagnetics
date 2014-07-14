

#include "../../mallinson_solution.h"
#include "../../magnetic_parameters.h"
#include "../../magnetics_helpers.h"
#include "../../../../src/generic/Vector.h"
#include <iostream>
#include <math.h>

using namespace oomph;
using namespace CompareSolutions;

double toradians(const double& degrees)
  {
    return (degrees/180) * 4.0*std::atan(1.0);
  }


int main()
  {

    const double pi = 4.0*std::atan(1.0);

    // parameters
    const double alpha = 0.5;
    const double gamma = 1.0;

    // Don't go quite a full arc because it technically takes infinite
    // time. Use same as mallinson in the paper.
    const double theta_start = toradians(10);
    const double theta_end = toradians(170);

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
        double theta_now = thetas[j];
        std::cout
          << theta_now << " "
          << switching_time(alpha, gamma, 1.0, 0.0,
                            theta_start, theta_now) << " "
          << analytic_phi(alpha, theta_start, theta_now)
          << std::endl;
      }

    return 0;
  }
