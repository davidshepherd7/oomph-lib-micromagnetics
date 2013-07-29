#ifndef OOMPH_MAGNETIC_MATERIALS_H
#define OOMPH_MAGNETIC_MATERIALS_H

#include "../../src/generic/Vector.h"
#include "../../src/generic/oomph_utilities.h"
#include "./vector_helpers.h"

namespace mag_parameters
{
  enum enum_crystalline_anisotropy_type
    {
      //??ds check the anisotropy fns...
      UNIAXIAL_CRYSTALLINE_ANISOTROPY
    };
}


namespace oomph
{

  using namespace StringConversion;

  // ============================================================
  /// A class to store magnetic material parameters.
  // ============================================================
  // TODO: multiple mesh normalisation - divide by static mean values?
  class MagneticParameters
  {

  public:
    MagneticParameters()
      : Gamma(2.211e5),
        Mu0(12.566370614e-7),  // in N/(A^2) (SI)
        Gilbert_damping(0.05),
        Saturation_magnetisation(1.0),
        Exchange_constant(0.5 * mu0()), // gives lex = 1
        K1(0.0),
        Easy_axis(3,0),
        Magnetostatic_debug_coeff(1.0),
        Exchange_debug_coeff(1.0),
        Boundary_exchange_debug_coeff(1.0),
        Ca_debug_coeff(1.0),
        Crystalline_ansiotropy_type(mag_parameters::UNIAXIAL_CRYSTALLINE_ANISOTROPY),
        Surface_anisotropy_enabled(0)
    {
      // For now easy axis is always z
      Easy_axis[2] = 1;
    }

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
    double distance_units() const {return magnetostatic_exchange_length();}
    bool surface_anisotropy_enabled() const {return Surface_anisotropy_enabled;}
    double mu0() const {return Mu0;}

    // Coefficients of each of the fields (before normalisation).
    double hk() const
    {return (2 * k1()) / (mu0() * saturation_magnetisation()); }
    double hex() const
    {return (2 * exchange_constant()) /
        (mu0() * saturation_magnetisation()
         * distance_units() * distance_units());}
    double hms() const
    {return saturation_magnetisation();}

    mag_parameters::enum_crystalline_anisotropy_type
    crystalline_ansiotropy_type() const
    {return Crystalline_ansiotropy_type;}

    // Normalised get functions
    double normalised_gamma() const
    {return 1.0;}
    double normalised_gilbert_damping() const
    {return gilbert_damping();}
    double normalised_saturation_magnetisation() const
    {return 1.0;}

    /// \short We are normalising field by Ms: this function allows it to
    /// be used elsewhere (e.g. for applied field).
    double field_normalisation_factor() const
    {
      // Note that all normalised fields will need to be changed (i.e. the
      // normalised_x() functions below) if you want to normalise by
      // something else
      return 1/saturation_magnetisation();
    }

    double normalised_hex() const
    {return 1.0 * exchange_debug_coeff();}

    double normalised_hk() const
    {
      return (2 * k1() / (mu0() *
                      std::pow(saturation_magnetisation(), 2)))
        * ca_debug_coeff();
    }

    double normalised_hms() const
    {return 1.0* magnetostatic_debug_coeff();}

    // ??ds need to make sure all field normalisation factors are the same for all
    // meshes or this falls apart!
    double time_normalisation_factor() const
    {return (1/gamma() * hms());}


    /// \short get easy axis direction unit vector in whatever direction it
    /// was given in.
    Vector<double> easy_axis() const
      {
        return Easy_axis;
      }

    /// \short get easy axis direction unit vector closest to m.
    Vector<double> nearest_easy_axis(const Vector<double> &m) const
      {
        Vector<double> easy_axis(Easy_axis.size(), 0);

        if(VectorOps::dot(m, Easy_axis) < 0)
          {
            for(unsigned j=0; j<Easy_axis.size(); j++)
              {
                easy_axis[j] = Easy_axis[j] * -1;
              }
          }
        else
          {
            easy_axis = Easy_axis;
          }

        return easy_axis;
      }


    // Get h_ca function and derivative
    void crystalline_ansiotropy_field(const double& t, const Vector<double>& x,
                                      const Vector<double>& m, Vector<double>& h_ca)
      const
    {

      if(crystalline_ansiotropy_type() ==
         mag_parameters::UNIAXIAL_CRYSTALLINE_ANISOTROPY)
        {
          Vector<double> e = easy_axis();

          // h_ca is easy direction multiplied by the magnitude (which is
          // negative if closer to -e than e).
          double magnitude = normalised_hk() * VectorOps::dot(e, m);

          h_ca.assign(e.size(), 0.0);
          for(unsigned i=0; i<e.size(); i++)
            {
              h_ca[i] = e[i] * magnitude;
            }
        }
      else
        {
              std::string error_msg
                = "Crystalline anisotropy type not implemented";
              throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
        }
    }


    // Get the appropriate exchange length for this material according to
    // the nmag user manual.
    double exchange_length() const
    {
      double l1 = magnetostatic_exchange_length();
      double l2 = magnetocrystalline_anisotropy_exchange_length();
      return std::min(l1,l2);
    }

    double magnetostatic_exchange_length() const
      {
        return std::sqrt( (2* exchange_constant() )
                          / (mu0() * saturation_magnetisation()
                             * saturation_magnetisation()));
      }

    double magnetocrystalline_anisotropy_exchange_length() const
      {
        if (k1() > 0)
          return std::sqrt( exchange_constant() / k1() );
        else
          return 1e30; // "infinity" (e30 instead of higher in case we ever
                       // use floats not doubles).
        //??ds throw error instead?
      }

    // "shape_fn_l2_at_x" is the shape function of the value we are differentiating
    // with respect to at the point x.
    void crystalline_ansiotropy_field_derivative
    (const double& t, const Vector<double>& x,
     const Vector<double>& m, const double shape_fn_l2_at_x,
     double dhcadm[3][3]) const
    {
      switch (crystalline_ansiotropy_type())
        {

        case mag_parameters::UNIAXIAL_CRYSTALLINE_ANISOTROPY:

          // Get easy axis
          Vector<double> e = easy_axis();

          // Calculate derivatives
          for(unsigned j=0; j<3; j++)
            {
              for(unsigned i=0; i<3; i++)
                {
                  dhcadm[i][j] = shape_fn_l2_at_x
                    * e[i] * e[j] * normalised_hk();
                }
            }
          break;
        }
    }

    void output(std::ostream& stream) const
    {
      stream << std::endl;
      stream << "Magnetic parameters are:" << std::endl;
      stream << "Gamma = " << gamma() << std::endl;
      stream << "Mu0 = " << mu0() << std::endl;
      stream << "Damping = " << gilbert_damping() << std::endl;
      stream << "M_s = "  << saturation_magnetisation() << std::endl;
      stream << "Exchange constant = " << exchange_constant() <<std::endl;
      stream << "K_1 = " << k1() << std::endl;
      stream << "Easy_axis = " << easy_axis() << std::endl;

      stream << "magnetostatic debug = " << magnetostatic_debug_coeff() << std::endl;
      stream << "exchange debug = " << exchange_debug_coeff() << std::endl;
      stream << "boundary exchange debug = " << boundary_exchange_debug_coeff() << std::endl;
      stream << "cryst anisotropy debug = " << ca_debug_coeff() << std::endl;
      stream << std::endl;

      stream << "This gives the following normalised coeffs:" << std::endl;
      stream << "Normalised_gamma = " <<    normalised_gamma() << std::endl;
      stream << "Normalised_gilbert_damping = " <<  normalised_gilbert_damping() << std::endl;
      stream << "Normalised_saturation_magnetisation = "
             <<  normalised_saturation_magnetisation() << std::endl;
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
    double& mu0() {return Mu0;}

    void set_cubic_anisotropy()
    {
      Crystalline_ansiotropy_type =
        mag_parameters::UNIAXIAL_CRYSTALLINE_ANISOTROPY;
    }

  private:
    double Gamma;
    double Mu0;
    double Gilbert_damping;

    double Saturation_magnetisation;
    double Exchange_constant;
    double K1;
    Vector<double> Easy_axis;

    /// Debug coefficients
    double Magnetostatic_debug_coeff;
    double Exchange_debug_coeff;
    double Boundary_exchange_debug_coeff;
    double Ca_debug_coeff;

    mag_parameters::enum_crystalline_anisotropy_type Crystalline_ansiotropy_type;

    bool Surface_anisotropy_enabled;
  };


  namespace Factories
  {

    /// \short Create a MagneticsParameters object.
    inline MagneticParameters* magnetic_parameters_factory(const std::string & parameters_name)
    {
      MagneticParameters* parameters_pt = 0;

      // Set properties as used in umag standard problem #4
      if(to_lower(parameters_name) == "mumag4")
        {
          parameters_pt = new MagneticParameters;
          parameters_pt->exchange_constant() = 1.3e-11; //J/m (1.3e-6 erg/cm)
          parameters_pt->saturation_magnetisation() = 8.0e5; // A/m (800 emu/cc)
          parameters_pt->k1() = 0.0;

          parameters_pt->mu0() = 12.566370614e-7;
          parameters_pt->gilbert_damping() = 0.02;
          parameters_pt->gamma() = 2.211e5; // m/(As)
        }

      // Set parameters as used in nmag's cubeoid example
      else if(to_lower(parameters_name) == "nmag-cubeoid")
        {
          parameters_pt = new MagneticParameters;
          parameters_pt->saturation_magnetisation() = 0.86e6; // A/m
          parameters_pt->exchange_constant() = 13.0e-12; // J/m
          parameters_pt->k1() = 0.0;

          parameters_pt->mu0() = 12.566370614e-7;
          parameters_pt->gilbert_damping() = 0.5;
          parameters_pt->gamma() = 2.210173e5; // m/(As)
        }

      // Set parameters as many coefficients as possible from the LLG. Set
      // damping = 0.5, lex = 1, hk = 0.
      else if(to_lower(parameters_name) == "simple-llg")
        {
          parameters_pt = new MagneticParameters;
          parameters_pt->saturation_magnetisation() = 1.0; // normalised units
          parameters_pt->exchange_constant() = 0.5; // this gives lex = 1
          parameters_pt->k1() = 0.0;
          parameters_pt->mu0() = 1;
          parameters_pt->gamma() = 1;
          parameters_pt->gilbert_damping() = 0.5;
        }

      // Set parameters as many coefficients as possible from the LLG. Set
      // damping = 0.5, lex = 1, hk = 1.
      else if(to_lower(parameters_name) == "simple-llg-anisotropy")
        {
          parameters_pt = magnetic_parameters_factory("simple-llg");
          parameters_pt->k1() = 0.5; // This gives hk = 1
        }

      // Set parameters as many coefficients as possible from the LLG. Set
      // damping = 1, lex = 1, hk = 0.
      else if(to_lower(parameters_name) == "simple-llg-max-damped")
        {
          parameters_pt = magnetic_parameters_factory("simple-llg");
          parameters_pt->gilbert_damping() = 1;
        }

      else
        {
          std::string error_msg = "Unrecognised parameters_name: "
            + to_lower(parameters_name);
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      return parameters_pt;
    }
  }


}

#endif
