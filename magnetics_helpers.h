#ifndef OOMPH_MAGNETICS_HELPERS_H
#define OOMPH_MAGNETICS_HELPERS_H

/*
description of file goes here
*/

#include "./vector_helpers.h"

#include "../../src/generic/Vector.h"
#include "../../src/generic/oomph_utilities.h"
#include "../../src/generic/oomph_definitions.h"

namespace HApp
{
  using namespace oomph;
  using namespace MathematicalConstants;
  using namespace StringConversion;
  using namespace VectorOps;

  typedef Vector<double> (*HAppFctPt)(const double& t, const Vector<double>&x);

  inline Vector<double> zero(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    return h;
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
    else if(field_name == "all_directions")
      {
        return &HApp::all_directions;
      }
    else if(field_name == "z_oscillating_p20")
      {
        return &HApp::z_oscillating_p20;
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

  typedef Vector<double> (*InitialMFctPt)(const double& t, const Vector<double>&x);

  inline Vector<double> x(const double& t, const Vector<double> &x)
  {
    Vector<double> m(3, 0.0);
    m[0] = 1.0;
    m[1] = 0.2;
    normalise(m);
    return m;
  }

  inline Vector<double> y(const double& t, const Vector<double> &x)
  {
    Vector<double> m(3, 0.0);
    m[1] = 1.0;
    m[2] = 0.2;
    normalise(m);
    return m;
  }

  inline Vector<double> z(const double& t, const Vector<double> &x)
  {
    Vector<double> m(3, 0.0);
    m[2] = 1.0;
    m[0] = 0.2;
    normalise(m);
    return m;
  }

  inline Vector<double> exactly_z(const double& t, const Vector<double> &x)
  {
    Vector<double> m(3, 0.0);
    m[2] = 1.0;
    return m;
  }

  inline Vector<double> xz(const double& t, const Vector<double> &x)
  {
    Vector<double> m(3, 0.0);
    m[0] = 1.0;
    m[2] = 1.0;
    normalise(m);
    return m;
  }

  inline Vector<double> smoothly_varying_m_helper(double l, const double& t,
                                                  const Vector<double> &x)
  {
    Vector<double> m(3,0.0);

    m[0] = sin(x[0]*2*Pi/l) + sin(x[1]*2*Pi/l);
    m[1] = cos(x[0]*2*Pi/l) + cos(x[1]*2*Pi/l);
    m[2] = 1.0 - m[0] - m[1];

    normalise(m);
    return m;
  }

  inline Vector<double> smoothly_varying_5(const double& t, const Vector<double> &x)
  {return smoothly_varying_m_helper(5.0, t, x);}
  inline Vector<double> smoothly_varying_50(const double& t, const Vector<double> &x)
  {return smoothly_varying_m_helper(50.0, t, x);}
  inline Vector<double> smoothly_varying_500(const double& t, const Vector<double> &x)
  {return smoothly_varying_m_helper(500.0, t, x);}

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
    else
      {
        throw OomphLibError("Unrecognised initial m name " + m_name,
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
  }
}



#endif
