#ifndef OOMPH_MAGNETIC_MATERIALS_H
#define OOMPH_MAGNETIC_MATERIALS_H

#include "generic.h"

namespace mag_parameters
{
  enum enum_crystalline_anisotropy_type
    {
      //??ds check the anisotropy fns...
      CUBIC_CRYSTALLINE_ANISOTROPY
    };

  // Magnetic constant (source: http://physics.nist.gov/)
  double mu0 = 12.566370614e-7; // in N/(A^2) (SI)
}

  // ============================================================
  /// A class to store magnetic material parameters.
  // ============================================================
  class MagneticParameters
  {
  public:
    MagneticParameters()
      : Gamma(1e-15), Gilbert_damping(0.05),
	Crystalline_ansiotropy_type(mag_parameters::CUBIC_CRYSTALLINE_ANISOTROPY)
    {}

    // Standard copy/assign constructors and deconstructor are ok: no pointers and
    // should never be any pointers.

    /// Get functions
    double gamma() const {return Gamma;}
    double gilbert_damping() const {return Gilbert_damping;}
    double saturation_magnetisation() const {return Saturation_magnetisation;}
    double exchange_constant() const {return Exchange_constant;}
    double k1() const {return K1;}
    mag_parameters::enum_crystalline_anisotropy_type crystalline_ansiotropy_type() const
    {return Crystalline_ansiotropy_type;}

    // Get h_ca function and derivative
    void crystalline_ansiotropy_field(const double& t, const Vector<double>& x,
				      const Vector<double>& m, Vector<double>& h_ca)
    {

      switch (crystalline_ansiotropy_type())
	{

	case mag_parameters::CUBIC_CRYSTALLINE_ANISOTROPY:

	  Vector<double> easy_axis(3,0.0);
	  unsigned ax = 0;

	  // Find which direction along the axis we want
	  if(m[ax] < 0)
	    easy_axis[ax] = -1.0;
	  else
	    easy_axis[ax] = +1.0;
	  h_ca = easy_axis;

	  // Multiply by the magnitude
	  double magnitude = VectorOps::dot(easy_axis,m);
	  for(unsigned i=0; i<h_ca.size(); i++)
	    h_ca[i] *= magnitude;

	  break;

	}
    }

    // "shape_fn_l2_at_x" is the shape function of the value we are differentiating
    // with respect to at the point x.
    void crystalline_ansiotropy_field_derivative
    (const double& t, const Vector<double>& x,
     const Vector<double>& m, const double shape_fn_l2_at_x,
     DenseMatrix<double>& dhcadm)
    {
      switch (crystalline_ansiotropy_type())
	{

	case mag_parameters::CUBIC_CRYSTALLINE_ANISOTROPY:

	  Vector<double> easy_axis(3,0.0);
	  unsigned ax = 0;

	  // Find which direction along the axis we want
	  if(m[ax] < 0)
	    easy_axis[ax] = -1.0;
	  else
	    easy_axis[ax] = +1.0;

	  for(unsigned j=0; j<3; j++)
	    for(unsigned i=0; i<3; i++)
	      dhcadm(i,j) = shape_fn_l2_at_x
		* easy_axis[i] * easy_axis[j];
	  break;
	}
    }


    /// Set functions
    double& gamma() {return Gamma;}
    double& gilbert_damping() {return Gilbert_damping;}
    double& saturation_magnetisation() {return Saturation_magnetisation;}
    double& exchange_constant() {return Exchange_constant;}
    double& k1() {return K1;}

    void set_cubic_anisotropy()
    {Crystalline_ansiotropy_type = mag_parameters::CUBIC_CRYSTALLINE_ANISOTROPY;}

    /// set properties for permalloy (Fe_20 Ni_80)
    void set_permalloy()
    {
      saturation_magnetisation() = 1.04/mag_parameters::mu0; // T/mu0 unitsm
      exchange_constant() = 7e-12; // J/m
      k1() = -2e3; // J/(m^3)
    }

    /// Set properties as used in umag standard problem #4
    void set_umag4()
    {
      exchange_constant() = 1.3e-11; //J/m (1.3e-6 erg/cm)
      saturation_magnetisation() = 8.0e5; // A/m (800 emu/cc)
      k1() = 0.0;

      gilbert_damping() = 0.02;
      gamma() = 2.211e5; // m/(As)
    }

    void set_FePt()
    {

    }

    // etc...

  private:
    double Gamma;
    double Gilbert_damping;

    double Saturation_magnetisation;
    double Exchange_constant;
    double K1;

    mag_parameters::enum_crystalline_anisotropy_type Crystalline_ansiotropy_type;
  };


#endif
