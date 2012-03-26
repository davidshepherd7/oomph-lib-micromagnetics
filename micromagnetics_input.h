#ifndef OOMPH_MICROMAGNETICS_INPUT_H
#define OOMPH_MICROMAGNETICS_INPUT_H

/*
  description of file goes here
*/

#include "generic.h"

using namespace oomph;

namespace oomph
{

  //======================================================================
  /// A class to contain all the inputs needed
  //======================================================================
  class MicromagInputs
  {
  private:

    const double M_s;

    const double Gilbert_c;
    const double Gymag_c;

    const double Exchange_c;

    const double Cryst_anis_const;
    const Vector<double> Easy_axis;

  public:

    /// Constructor
    MicromagInputs(const double& m_s,
		   const double& gilbert_c,
		   const double& gymag_c,
		   const double& exchange_c,
		   const double& k_1,
		   const Vector<double>& easy_axis) :
      M_s(m_s), Gilbert_c(gilbert_c), Gymag_c(gymag_c),
      Exchange_c(exchange_c), Cryst_anis_const(k_1), Easy_axis(easy_axis)
    {}

    /// Constructor with some default values
    MicromagInputs() :
      M_s(1.0), Gilbert_c(0.05), Gymag_c(1.0),
      Exchange_c(1.0), Cryst_anis_const(1.0), Easy_axis(3)
    {Easy_axis[0] = 0; Easy_axis[1] = 0; Easy_axis[2] = 1;}

    double llg_precession_nn(const double& t, const Vector<double>& x)
    {return Gymag_c/(1 + Gilbert_c*Gilbert_c);}

    double llg_precession(const double& t, const Vector<double>& x)
    {return 1.0;}

    double llg_damping_nn(const double& t, const Vector<double>& x)
    {return llg_precession_nn(t,x) * (Gilbert_c/sat_mag(t,x));}

    double llg_damping(const double& t, const Vector<double>& x)
    {return Gilbert_c;}

    double sat_mag(const double& t, const Vector<double>& x)
    {return m_s;}

    void cryst_anis_field(const double& t, const Vector<double>& x,
			  const Vector<double>& M, Vector<double>& H_cryst_anis)
    {
      H_cryst_anis[0] = 0.0;
      H_cryst_anis[1] = 0.0;
      H_cryst_anis[2] = 0.0;
    }

    void applied_field(const double& t, const Vector<double>& x,
		       Vector<double>& H_applied)
    {
      H_applied[0] = 0.0;
      H_applied[1] = 0.0;
      H_applied[2] = 0.0;
    }

    double exchange_coeff(const double& t, const Vector<double>& x)
    {
      return exchange_c;
    }
  };

} // End of oomph namespace

#endif
