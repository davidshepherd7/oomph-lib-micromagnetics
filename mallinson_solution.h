#ifndef OOMPH_MALLINSON_SOLUTION_H
#define OOMPH_MALLINSON_SOLUTION_H

/*
  description of file goes here
*/

#include "../../src/generic/Vector.h"
#include "magnetic_parameters.h"
#include "magnetics_helpers.h"

namespace oomph
{
  using namespace MathematicalConstants;

  // Calculations for the switching time and phi value for switching between
  // two angles in the very simple case when:
  // * H is aligned with the easy axis and constant
  // * There is no spatial variation
  // From Mallinson2000.
  namespace CompareSolutions
  {

    inline double cart2theta(const Vector<double> &m)
    {
      double r = VectorOps::two_norm(m);
      return std::acos(m[2]/r);
    }

    inline double cart2phi(const Vector<double> &m)
    {
      double result_in_mpi_pi = std::atan2(m[1],m[0]);
      if (result_in_mpi_pi > 0)
        return result_in_mpi_pi - 2*Pi;
      else
        return result_in_mpi_pi;
    }

    inline double switching_time(const double &alpha,
                          const double &gamma,
                          const double &H,
                          const double &H_k,
                          const double &theta_start,
                          const double &theta_now)
    {
      using namespace std;
      return ( (pow(alpha,2) + 1) / (alpha * gamma) )
        * ( 1 / ( pow(H,2) - pow(H_k,2)))
        * (H * log((tan(theta_now/2))
                   / (tan(theta_start/2)))

           + H_k * log((H - H_k * cos(theta_start))
                       / (H - H_k * cos(theta_now)))

           + H_k * log( sin(theta_now)
                        / sin(theta_start)));
    }


    inline double switching_time_wrapper(const MagneticParameters* const parameters_pt,
                                         const Vector<double>& initm,
                                         const Vector<double> &m)
    {
      Vector<double> x; //dummy
      Vector<double> H = HApp::minus_z(0, x);

      double theta_start = cart2theta(initm);
      double theta_now = cart2theta(m);

      double analytical_time =
        switching_time(parameters_pt->gilbert_damping(),
                       1.0,
                       std::abs(H[2]),
                       parameters_pt->normalised_hk(),
                       theta_start,
                       theta_now);

      return analytical_time;
    }


    inline double switching_time_wrapper(const MagneticParameters* const parameters_pt,
                                         const Vector<double> &m)
    {
      // Get initial magnetisation for m0 ~= z
      Vector<double> x; //dummy
      Vector<double> initial_values = InitialM::z(0, x);
      Vector<double> initm;
      initm.assign(initial_values.begin()+2, initial_values.end());

      return switching_time_wrapper(parameters_pt, initm, m);
    }


    inline double analytic_phi(const double &alpha,
                               const double &theta_start,
                               const double &theta_now)
    {
      double phi = (-1/alpha) * std::log((std::tan(theta_now/2))
                                         /(std::tan(theta_start/2)));
      return std::fmod(phi,2*Pi); // map into range 0:2pi
    }

    inline double analytic_phi_wrapper(const MagneticParameters* const parameters_pt,
                                       const Vector<double> &m)
    {
      Vector<double> x; // dummy!

      // Get initial magnetisation
      Vector<double> initial_values = InitialM::z(0, x);
      Vector<double> initm;
      initm.assign(initial_values.begin()+2, initial_values.end());
      // (can't use range constructor because they are not implemented for
      // oomph vectors, for some reason...)

      double theta_start = cart2theta(initm);
      double theta_now = cart2theta(m);

      return analytic_phi(parameters_pt->gilbert_damping(),
                          theta_start,
                          theta_now);
    }

    inline double phi_error(const MagneticParameters* const parameters_pt,
                            const Vector<double> &m)
    {
      double phi_now = cart2phi(m);
      double analytical_phi = analytic_phi_wrapper(parameters_pt,m);
      double result = std::abs(analytical_phi - phi_now);
      return result;
    }


    inline Vector<double> m_exact(const MagneticParameters& parameters,
                                  const Vector<double>& m0,
                                  const double& time)
    {
      const Vector<double> dummy;

#ifdef PARANOID
      if(cart2phi(m0) != 0)
        {
          std::string err = "Not sure if this works for phi0 != 0.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      if(parameters.normalised_hk() != 0)
        {
          std::string err = "Won't work with anisotropy.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      if(parameters.h_app(0, dummy) != parameters.h_app(std::sqrt(2.0), dummy))
        {
          std::string err = "Won't work with non-constant field.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      // Get parameters
      const double theta0 = cart2theta(m0);
      const double alpha = parameters.gilbert_damping();
      const double H = VectorOps::two_norm(parameters.h_app(0, dummy));

      // Calculate exact magnetisation
      const double expterm = std::exp(time * H * alpha/(alpha*alpha + 1));
      const double theta = 2*std::atan(std::tan(theta0/2) * expterm);
      const double phi = analytic_phi(alpha, theta0, theta);

      // Return in cartesian coords
      Vector<double> rpolarazi(3, 0.0);
      rpolarazi[0] = 1.0;
      rpolarazi[1] = theta;
      rpolarazi[2] = phi;
      return VectorOps::sphpolar_to_cart(rpolarazi);
    }

  }


} // End of oomph namespace

#endif
