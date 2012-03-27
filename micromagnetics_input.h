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
    Vector<double> Easy_axis; //??ds ideally would be const... I think

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

    /// Check that all values are filled in and valid
    void self_test() const;

    double llg_precession(const double& t, const Vector<double>& x) const
    {return 1.0;}

    double llg_damping(const double& t, const Vector<double>& x) const
    {return Gilbert_c;}

    double sat_mag(const double& t, const Vector<double>& x) const
    {return 1.0;}

    void cryst_anis_field(const double& t, const Vector<double>& x,
			  const Vector<double>& M, Vector<double>& H_cryst_anis) const
    {
      double dot_product = M[0]*Easy_axis[0] + M[1]*Easy_axis[1] + M[2]*Easy_axis[2];
      H_cryst_anis[0] = Easy_axis[0]*dot_product;
      H_cryst_anis[1] = Easy_axis[1]*dot_product;
      H_cryst_anis[2] = Easy_axis[2]*dot_product;
    }

    void applied_field(const double& t, const Vector<double>& x,
		       Vector<double>& H_applied) const
    {
      H_applied[0] = 0.0;
      H_applied[1] = 0.0;
      H_applied[2] = 0.0;
    }

    double exchange_coeff(const double& t, const Vector<double>& x) const
    {
      return Exchange_c;
    }
  };


  //======================================================================
  /// Check that all the inputs are included and valid.
  //======================================================================
  void MicromagInputs::self_test() const
  {
    //??ds do this
  }

} // End of oomph namespace

#endif
