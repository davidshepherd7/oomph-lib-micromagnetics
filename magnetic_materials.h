#ifndef OOMPH_MAGNETIC_MATERIALS_H
#define OOMPH_MAGNETIC_MATERIALS_H

#include "generic.h"

namespace mag_parameters
{
  enum enum_crystalline_anisotropy_type
    {
      //??ds check the anisotropy fns...
      UNIAXIAL_CRYSTALLINE_ANISOTROPY
    };

  // Magnetic constant (source: http://physics.nist.gov/)
  double mu0 = 12.566370614e-7; // in N/(A^2) (SI)
}

using namespace mag_parameters;

namespace oomph
{

  // ============================================================
  /// A class to store magnetic material parameters.
  // ============================================================
  // TODO: normalisation, multiple mesh normalisation - divide by static mean values?
  class MagneticParameters
  {

  public:
    MagneticParameters()
      : Gamma(2.211e5), Gilbert_damping(0.05),
        Saturation_magnetisation(1.0),
        Exchange_constant(0.5 * mag_parameters::mu0), // gives exchange strength = 1
        K1(0.0),
        Distance_units(1e-9),
        Magnetostatic_debug_coeff(1.0),
        Exchange_debug_coeff(1.0),
        Boundary_exchange_debug_coeff(1.0),
        Ca_debug_coeff(1.0),
        Crystalline_ansiotropy_type(mag_parameters::UNIAXIAL_CRYSTALLINE_ANISOTROPY),
        Surface_anisotropy_enabled(0)
    {}

    // Standard copy/assign constructors and deconstructor are ok: no pointers and
    // should never be any pointers.
    /// Get functions
    double gamma() const {return Gamma;}
    double gilbert_damping() const {return Gilbert_damping;}
    double saturation_magnetisation() const {return Saturation_magnetisation;}
    double exchange_constant() const {return Exchange_constant;}
    double k1() const {return K1;}
    double magnetostatic_debug_coeff() const {return Magnetostatic_debug_coeff;}
    double exchange_debug_coeff() const {return Exchange_debug_coeff;}
    double boundary_exchange_debug_coeff() const {return Boundary_exchange_debug_coeff;}
    double ca_debug_coeff() const {return Ca_debug_coeff;}
    double distance_units() const {return Distance_units;}
    bool surface_anisotropy_enabled() const {return Surface_anisotropy_enabled;}

    // Coefficients of each of the fields (before normalisation).
    double hk() const
    {return (2 * k1()) / (mag_parameters::mu0 * saturation_magnetisation()); }
    double hex() const
    {return (2 * exchange_constant()) /
        (mag_parameters::mu0 * saturation_magnetisation()
         * distance_units() * distance_units());}
    double hms() const
    {return saturation_magnetisation();}

    mag_parameters::enum_crystalline_anisotropy_type crystalline_ansiotropy_type() const
    {return Crystalline_ansiotropy_type;}

    // Normalised get functions
    double normalised_gamma() const
    {return 1.0;}
    double normalised_gilbert_damping() const
    {return gilbert_damping();}
    double normalised_saturation_magnetisation() const
    {return 1.0;}
    double normalised_hex() const
    {return hex() * field_normalisation_factor() * exchange_debug_coeff();}
    double normalised_hk() const //??ds probably need some distance conversion in here...
    {return hk() * field_normalisation_factor() * ca_debug_coeff();}
    double normalised_hms() const
    {return hms() * field_normalisation_factor() * magnetostatic_debug_coeff();}

    // Normalise by magnetostatic field strength
    double field_normalisation_factor() const
    {
      // return 1/hex_str;
      return 1/hms();
    }

    // ??ds need to make sure all field normalisation factors are the same for all
    // meshes or this falls apart!
    double time_normalisation_factor() const
    {return (1/gamma()) * field_normalisation_factor();}

    // Get h_ca function and derivative
    void crystalline_ansiotropy_field(const double& t, const Vector<double>& x,
                                      const Vector<double>& m, Vector<double>& h_ca)
      const
    {

      switch (crystalline_ansiotropy_type())
        {

        case mag_parameters::UNIAXIAL_CRYSTALLINE_ANISOTROPY:

          Vector<double> easy_axis(3,0.0);
          unsigned ax = 0;

          // Find which direction along the axis we want
          if(m[ax] < 0)
            easy_axis[ax] = -1.0;
          else
            easy_axis[ax] = +1.0;
          h_ca = easy_axis;

          // Multiply by the magnitude
          double magnitude = std::abs(VectorOps::dot(easy_axis,m));
          for(unsigned i=0; i<h_ca.size(); i++)
            h_ca[i] *= magnitude * normalised_hk();

          break;

        }
    }

    // Get the appropriate exchange length for this material according the the
    // nmag user manual.
    double exchange_length() const
    {
      double l1 = std::sqrt( (2* exchange_constant() )
                             / (mag_parameters::mu0 * saturation_magnetisation()
                                * saturation_magnetisation()));
      double l2;
      if (k1() > 0)
        l2 = std::sqrt( exchange_constant() / k1() );
      else
        l2 = 1e30; // infinite (e30 instead of higher in case we ever use floats
      // not doubles).

      return std::min(l1,l2);
    }

    // "shape_fn_l2_at_x" is the shape function of the value we are differentiating
    // with respect to at the point x.
    void crystalline_ansiotropy_field_derivative
    (const double& t, const Vector<double>& x,
     const Vector<double>& m, const double shape_fn_l2_at_x,
     DenseMatrix<double>& dhcadm) const
    {
      switch (crystalline_ansiotropy_type())
        {

        case mag_parameters::UNIAXIAL_CRYSTALLINE_ANISOTROPY:

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
                * easy_axis[i] * easy_axis[j] * normalised_hk() ;
          break;
        }
    }

    void output(std::ostream& stream) const
    {
      stream << "Magnetic parameters are:" << std::endl;
      stream << "Gamma = " << gamma() << std::endl;
      stream << "Damping = " << gilbert_damping() << std::endl;
      stream << "M_s = "  << saturation_magnetisation() << std::endl;
      stream << "Exchange constant = " << exchange_constant() <<std::endl;
      stream << "K_1 = " << k1() << std::endl;
      stream << std::endl;

      stream << "This gives the following normalised coeffs:" << std::endl;
      stream << "Normalised_gamma = " <<    normalised_gamma() << std::endl;
      stream << "Normalised_gilbert_damping = " <<  normalised_gilbert_damping() << std::endl;
      stream << "Normalised_saturation_magnetisation = " <<  normalised_saturation_magnetisation() << std::endl;
      stream << "Normalised_hex = " <<  normalised_hex() << std::endl;
      stream << "Normalised_hk = " <<  normalised_hk() << std::endl;
      stream << "Normalised_hms = " <<  normalised_hms() << std::endl;
      stream << std::endl;

      stream << "Exchange length for this material is: "
             << exchange_length() << std::endl;
      stream << "Your distance units are: " << distance_units()
             << "m" << std::endl;
      stream << "All elements should be smaller than the exchange length (and there should be code to check this...)." << std::endl;
    }

    /// Set functions
    double& gamma() {return Gamma;}
    double& gilbert_damping() {return Gilbert_damping;}
    double& saturation_magnetisation() {return Saturation_magnetisation;}
    double& exchange_constant() {return Exchange_constant;}
    double& k1() {return K1;}
    double& magnetostatic_debug_coeff() {return Magnetostatic_debug_coeff;}
    double& exchange_debug_coeff() {return Exchange_debug_coeff;}
    double& boundary_exchange_debug_coeff() {return Boundary_exchange_debug_coeff;}
    double& ca_debug_coeff() {return Ca_debug_coeff;}
    double& distance_units() {return Distance_units;}

    void set_cubic_anisotropy()
    {Crystalline_ansiotropy_type = mag_parameters::UNIAXIAL_CRYSTALLINE_ANISOTROPY;}

    // /// set properties for permalloy (Fe_20 Ni_80)
    // void set_permalloy()
    // {
    //   saturation_magnetisation() = 1.04/mag_parameters::mu0; // T/mu0 unitsm
    //   exchange_constant() = 7e-12; // J/m
    //   k1() = -2e3; // J/(m^3)
    // }

    /// Set properties as used in umag standard problem #4
    void set_mumag4()
    {
      exchange_constant() = 1.3e-11; //J/m (1.3e-6 erg/cm)
      saturation_magnetisation() = 8.0e5; // A/m (800 emu/cc)
      k1() = 0.0;

      gilbert_damping() = 0.02;
      gamma() = 2.211e5; // m/(As)
    }

    void set_nmag_rectangle()
    {
      saturation_magnetisation() = 0.86e6; // A/m
      exchange_constant() = 13.0e-12; // J/m
      k1() = 0.0;

      gilbert_damping() = 0.5;
      gamma() = 2.210173e5; // m/(As)
    }

    void set_simple_llg_parameters()
    {
      saturation_magnetisation() = 1.0; // normalised units
      exchange_constant() = 0.5 * mag_parameters::mu0; // this gives hex = 1
      k1() = 0.0;
      gamma() = 1;
      distance_units() = 1;

      // The only real parameter left here (can be varied):
      gilbert_damping() = 0.5;
    }


    // etc...

  private:
    double Gamma;
    double Gilbert_damping;

    double Saturation_magnetisation;
    double Exchange_constant;
    double K1;

    double Distance_units;

    /// Debug coefficients
    double Magnetostatic_debug_coeff;
    double Exchange_debug_coeff;
    double Boundary_exchange_debug_coeff;
    double Ca_debug_coeff;

    mag_parameters::enum_crystalline_anisotropy_type Crystalline_ansiotropy_type;

    bool Surface_anisotropy_enabled;
  };

}

#endif
