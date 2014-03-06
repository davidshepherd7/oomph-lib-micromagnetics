#ifndef OOMPH_MAGNETICS_HELPERS_H
#define OOMPH_MAGNETICS_HELPERS_H

/*
description of file goes here
*/

#include "vector_helpers.h"
#include "micromag_types.h"

#include <deque>

#include "../../src/generic/Vector.h"
#include "../../src/generic/oomph_utilities.h"
#include "../../src/generic/oomph_definitions.h"


namespace oomph
{

  /// Generic zero initial condition function
  template<unsigned NVAL>
  inline Vector<double> zero_initial_condition(const double& t,
                                               const Vector<double> &x)
  {
    Vector<double> v(NVAL, 0.0);
    return v;
  }

namespace HApp
{
  using namespace oomph;
  using namespace MathematicalConstants;
  using namespace StringConversion;
  using namespace VectorOps;

  inline Vector<double> zero(const double& t, const Vector<double> &x)
  {
    return zero_initial_condition<3>(t, x);
  }

  inline Vector<double> x(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[0] = 1.1;
    return h;
  }

  inline Vector<double> y(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[1] = 1.1;
    return h;
  }

  inline Vector<double> z(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[2] = 1.1;
    return h;
  }

  inline Vector<double> minus_z(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[2] = -1.1;
    return h;
  }

  inline Vector<double> minus_x(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[0] = -1.1;
    return h;
  }

  inline Vector<double> all_directions(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 1.1);
    return h;
  }

  inline Vector<double> z_oscillating_p20(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[0] = 0.1;
    h[2] = 1.1 * std::sin((MathematicalConstants::Pi/20) * t);
    return h;
  }

  inline Vector<double> mumag4_initial(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    double t0 = 5;

    // Strong, decay "slowly" to zero
    h[0] = std::max(2.5 * (t0 - t)/t0, 0.0);
    h[1] = h[0];
    h[2] = h[0];

    return h;
  }

  inline Vector<double> mumag4_field1(const double& t, const Vector<double> &x)
  {
    double mu0 = 4*Pi*1e-7;
    double Ms = 8e5;
    double a = 1e-3/(mu0*Ms);

    Vector<double> h(3, 0.0);
    h[0] = -24.6*a;
    h[1] = 4.3*a;
    h[2] = 0.0;

    return h;
  }

  inline Vector<double> mumag4_field2(const double& t, const Vector<double> &x)
  {
    double mu0 = 4*Pi*1e-7;
    double Ms = 8e5;
    double a = 1e-3/(mu0*Ms);

    Vector<double> h(3, 0.0);
    h[0] = -35.5*a;
    h[1] = -6.3*a;
    h[2] = 0.0;

    return h;
  }

  inline Vector<double> non_uniform_z_helper
  (const double& t, const Vector<double> &x, double l)
  {
    Vector<double> h(3, 0.0);
    h[0] = 0.1;
    h[2] = 1.1 * std::sin(MathematicalConstants::Pi * x[0] / l);
    return h;
  }

  inline Vector<double> non_uniform_z_5(const double& t, const Vector<double> &x)
  {return non_uniform_z_helper(t, x, 5);}
  inline Vector<double> non_uniform_z_50(const double& t, const Vector<double> &x)
  {return non_uniform_z_helper(t, x, 50);}
  inline Vector<double> non_uniform_z_500(const double& t, const Vector<double> &x)
  {return non_uniform_z_helper(t, x, 500);}

  inline Vector<double> tanhx_minus_z(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[2] = -1.1*tanh(5*x[0]);
    return h;
  }

  inline Vector<double> minus_z_above_x0(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    if(x[0] > 0)
      {
        h[2] = -1.1;
      }
    return h;
  }

  inline Vector<double> tanhx_minus_x(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[0] = -1.1*tanh(5*x[0]);
    return h;
  }

  inline Vector<double> minus_x_above_x0(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    if(x[0] > 0)
      {
        h[0] = -1.1;
      }
    return h;
  }


  inline Vector<double> smooth_start_minus_z(const double& t,
                                             const Vector<double> &x)
  {
    double a = 5;
    double smoothing = (1 - std::exp(-a*t));
    Vector<double> h = HApp::z(t, x);
    for(unsigned j=0; j<3; j++)
      {
        h[j] *= smoothing;
      }
    return h;
  }

  inline HAppFctPt h_app_factory(const std::string& field_name)
  {
    if(field_name == "zero")
      {
        return &HApp::zero;
      }
    else if(field_name == "x")
      {
        return &HApp::x;
      }
    else if(field_name == "y")
      {
        return &HApp::y;
      }
    else if(field_name == "z")
      {
        return &HApp::z;
      }
    else if(field_name == "minus_z")
      {
        return &HApp::minus_z;
      }
    else if(field_name == "minus_x")
      {
        return &HApp::minus_x;
      }
    else if(field_name == "all_directions")
      {
        return &HApp::all_directions;
      }
    else if(field_name == "z_oscillating_p20")
      {
        return &HApp::z_oscillating_p20;
      }
    else if(field_name == "non_uniform_z_5")
      {
        return &HApp::non_uniform_z_5;
      }
    else if(field_name == "non_uniform_z_50")
      {
        return &HApp::non_uniform_z_50;
      }
    else if(field_name == "non_uniform_z_500")
      {
        return &HApp::non_uniform_z_500;
      }
    else if(field_name == "tanhx_minus_z")
      {
        return &HApp::tanhx_minus_z;
      }
    else if(field_name == "minus_z_above_x0")
      {
        return &HApp::minus_z_above_x0;
      }
    else if(field_name == "tanhx_minus_x")
      {
        return &HApp::tanhx_minus_x;
      }
    else if(field_name == "minus_x_above_x0")
      {
        return &HApp::minus_x_above_x0;
      }
    else if(field_name == "mumag4_initial")
      {
        return &HApp::mumag4_initial;
      }
    else if(field_name == "mumag4_field1")
      {
        return &HApp::mumag4_field1;
      }
    else if(field_name == "mumag4_field2")
      {
        return &HApp::mumag4_field2;
      }
    else if(field_name == "smooth_start_minus_z")
      {
        return &HApp::smooth_start_minus_z;
      }
    else
      {
        throw OomphLibError("Unrecognised field name " + field_name,
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
  }
}

namespace InitialM
{

  using namespace oomph;
  using namespace MathematicalConstants;
  using namespace StringConversion;
  using namespace VectorOps;

  // Note that dof storage in micromag elements is actually in the order
  // [phi, phi1, mx, my, mz], so we put mx in element 2 of the vector etc.

  inline Vector<double> x(const double& t, const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);
    m[2] = 1.0;
    m[3] = 0.2;
    normalise(m);
    return m;
  }

  inline Vector<double> y(const double& t, const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);
    m[3] = 1.0;
    m[4] = 0.2;
    normalise(m);
    return m;
  }

  inline Vector<double> z(const double& t, const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);
    m[4] = 1.0;
    m[2] = 0.2;
    normalise(m);
    return m;
  }

  inline Vector<double> xyz(const double& t, const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);
    m[2] = 1.0;
    m[3] = 1.0;
    m[4] = 1.0;
    normalise(m);
    return m;
  }

  inline Vector<double> xy(const double& t, const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);
    m[2] = 1.0;
    m[3] = 1.0;
    normalise(m);
    return m;
  }

  inline Vector<double> exactly_z(const double& t, const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);
    m[4] = 1.0;
    return m;
  }

  inline Vector<double> xz(const double& t, const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);
    m[2] = 1.0;
    m[4] = 1.0;
    normalise(m);
    return m;
  }

  inline Vector<double> smoothly_varying_m_helper(double l, const double& t,
                                                  const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);
    m[2] = sin(x[0]*2*Pi/l) + sin(x[1]*2*Pi/l);
    m[3] = cos(x[0]*2*Pi/l) + cos(x[1]*2*Pi/l);
    m[4] = 1.0 - m[2] - m[3];
    normalise(m);

    return m;
  }

  inline Vector<double> smoothly_varying_5(const double& t, const Vector<double> &x)
  {return smoothly_varying_m_helper(5.0, t, x);}
  inline Vector<double> smoothly_varying_50(const double& t, const Vector<double> &x)
  {return smoothly_varying_m_helper(50.0, t, x);}
  inline Vector<double> smoothly_varying_500(const double& t, const Vector<double> &x)
  {return smoothly_varying_m_helper(500.0, t, x);}
    inline Vector<double> smoothly_varying_5000(const double& t, const Vector<double> &x)
  {return smoothly_varying_m_helper(5000.0, t, x);}


  inline InitialMFctPt initial_m_factory(const std::string& m_name)
  {
    if(m_name == "x")
      {
        return &InitialM::x;
      }
    else if(m_name == "y")
      {
        return &InitialM::y;
      }
    else if(m_name == "z")
      {
        return &InitialM::z;
      }
    else if(m_name == "xyz")
      {
        return &InitialM::xyz;
      }
    else if(m_name == "xy")
      {
        return &InitialM::xy;
      }
    else if(m_name == "xz")
      {
        return &InitialM::xz;
      }
    else if(m_name == "exactly_z")
      {
        return &InitialM::exactly_z;
      }
    else if(m_name == "smoothly_varying_5")
      {
        return &InitialM::smoothly_varying_5;
      }
    else if(m_name == "smoothly_varying_50")
      {
        return &InitialM::smoothly_varying_50;
      }
    else if(m_name == "smoothly_varying_500")
      {
        return &InitialM::smoothly_varying_500;
      }
    else if(m_name == "smoothly_varying_5000")
      {
        return &InitialM::smoothly_varying_5000;
      }
    else
      {
        throw OomphLibError("Unrecognised initial m name " + m_name,
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
  }
}


  namespace MagnetostaticFieldFunctions
  {

    inline Vector<double> sphere(const double& time, const Vector<double>& x,
                                 const Vector<double>& m)
    {
      Vector<double> hms(3, 0.0);

      for(unsigned j=0; j<3; j++)
        {
          hms[j] = -m[j]/3;
        }

      return hms;
    }


    inline MagnetostaticFieldFctPt ms_factory(const std::string& ms_name)
    {
      if(ms_name == "sphere")
        {
          return &MagnetostaticFieldFunctions::sphere;
        }
      else
        {
          throw OomphLibError("Unrecognised initial magnetostatic field name "
                              + ms_name,
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }

  }

class MyProblem;

namespace MManipulation
{

  /// \short Compute the effective damping constant (alpha) for the
  /// previous time step (see Albuquerque2001).
  double alt_effective_damping_used(MyProblem* problem_pt,
                                    std::deque<double>& previous_energies);


  /// \short Compute the effective damping constant (alpha) for the
  /// previous time step (see Albuquerque2001).
  double effective_damping_used(MyProblem* problem_pt);


  double exchange_energy(MyProblem* problem_pt);


  double zeeman_energy(MyProblem* problem_pt);

  double crystalline_anisotropy_energy(MyProblem* problem_pt);


  double magnetostatic_energy(MyProblem* problem_pt);

  double integral_of_dmdt_squared(MyProblem* problem_pt);

  double dEnergydt(MyProblem* problem_pt);

  double alt_dEnergydt(MyProblem* problem_pt,
                       std::deque<double>& previous_energies);
}

  class TetMeshBase;

  namespace MeshCreationHelpers
  {
    /// The trivial factory function for an ELEMENT.
    template<class ELEMENT>
    inline  FiniteElement* new_element() { return new ELEMENT;}


    void brick2tet(const Mesh& brick_mesh,
                   ElementFactoryFctPt element_factory_fpt,
                   TetMeshBase& out_mesh);
  }

}


#endif
