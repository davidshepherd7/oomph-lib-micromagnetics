#ifndef OOMPH_MAGNETIC_MATERIALS_H
#define OOMPH_MAGNETIC_MATERIALS_H

#include "../../src/generic/Vector.h"
#include "../../src/generic/oomph_utilities.h"
#include "./vector_helpers.h"
#include "micromag_types.h"


namespace oomph
{

  namespace MagParameterHelpers
  {
    using namespace MathematicalConstants;

    inline double Astar_to_A(const double& Astar, const double& l, const double& Ms)
    {
      double mu0 = 4*Pi*1e-7;
      return (2*Astar)/((l*l) * mu0 * (Ms*Ms));
    }

  }

  /// Base class for magnetic parameter storage.
  class MagneticParametersBase
  {
  public:
    MagneticParametersBase() {}

    virtual ~MagneticParametersBase() {}

    virtual double normalised_hex() const = 0;
    virtual double normalised_hk() const  = 0;
    virtual double normalised_hms() const = 0;
    virtual Vector<double> easy_axis() const = 0;
    virtual double happ_normalisation_factor() const = 0;


    virtual bool surface_anisotropy_enabled() const {return false;}


    /// \short get easy axis direction unit vector closest to m.
    Vector<double> nearest_easy_axis(const Vector<double> &m) const
    {
      Vector<double> easy_axis = this->easy_axis();
      Vector<double> nearest_easy_axis(easy_axis.size(), 0);

      if(VectorOps::dot(m, easy_axis) < 0)
        {
          for(unsigned j=0; j<easy_axis.size(); j++)
            {
              nearest_easy_axis[j] = easy_axis[j] * -1;
            }
        }
      else
        {
          nearest_easy_axis = easy_axis;
        }

      return nearest_easy_axis;
    }

    virtual void output(std::ostream& stream) const
      {
        stream
          << "normalised_hex " << normalised_hex() << std::endl
          << "normalised_hk " << normalised_hk() << std::endl
          << "normalised_hms " << normalised_hms() << std::endl
          << "easy_axis " << easy_axis() << std::endl
          << "happ_normalisation_factor " << happ_normalisation_factor() << std::endl
          << "surface_anisotropy_enabled " << surface_anisotropy_enabled() << std::endl
          ;
      }
  };

  /// A class to store magnetic material parameters for simple cases: when
  /// we don't have things like varying Ms to worry about. More like a
  /// struct really because it's so simple.
  class MagneticParameters : public MagneticParametersBase
  {
  public:

    void build(double A, double K1, double l, double Ms,
               double gilbert_damping)
      {
        double mu0 = 4 * MathematicalConstants::Pi * 1e-7;
        Exchange_coeff = (2*A)/(mu0*Ms*l*l);
        Anisotropy_coeff = K1;
        Gilbert_damping = gilbert_damping;
      }

    MagneticParameters() : Easy_axis(3, 0.0)
    {
      Gilbert_damping = 0.5;
      Exchange_coeff = 1.0;
      Anisotropy_coeff = 0.0;
      Magnetostatic_debug_coeff = 1.0;
      Applied_field_debug_coeff = 1.0;

      // For now easy axis is always z
      Easy_axis[2] = 1;
    }

    double Gilbert_damping;
    double Exchange_coeff;
    double Anisotropy_coeff;
    Vector<double> Easy_axis;
    double Magnetostatic_debug_coeff;
    double Applied_field_debug_coeff;
    HAppFctPt Applied_field_fct_pt;

    double normalised_hex() const {return Exchange_coeff;}
    double normalised_hk() const {return Anisotropy_coeff;}
    double normalised_hms() const {return Magnetostatic_debug_coeff;}
    double gilbert_damping() const {return Gilbert_damping;}
    double damping() const {return Gilbert_damping;}
    Vector<double> easy_axis() const {return Easy_axis;}
    double happ_normalisation_factor() const {return Applied_field_debug_coeff;}

    Vector<double> h_app(const double& time, const Vector<double>& x) const
      {
#ifdef PARANOID
        if(Applied_field_fct_pt == 0)
          {
        std::string err = "Applied_field_fct_pt is null!";
        throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
#endif
        Vector<double> h_app = Applied_field_fct_pt(time, x);

        for(unsigned j=0; j<3; j++)
          {
            h_app[j] *= happ_normalisation_factor();
          }

        return h_app;
      }


    // Get h_ca function and derivative
    void crystalline_ansiotropy_field(const double& t, const Vector<double>& x,
                                      const Vector<double>& m, Vector<double>& h_ca)
      const
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


    // "shape_fn_l2_at_x" is the shape function of the value we are differentiating
    // with respect to at the point x.
    void crystalline_ansiotropy_field_derivative
    (const double& t, const Vector<double>& x,
     const Vector<double>& m, const double shape_fn_l2_at_x,
     double dhcadm[3][3]) const
    {
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
    }


  };


  /// Allowing varying Ms, gamma etc. is much more messy... The class to
  /// handle it will go here:
  class ComplexMagneticParameters : public MagneticParameters
  {
  public:
    ComplexMagneticParameters()
    {
      throw OomphLibError("Class not yet implemented",
                          OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
    }


    // // Get the appropriate exchange length for this material according to
    // // the nmag user manual.
    // double exchange_length() const
    // {
    //   double l1 = magnetostatic_exchange_length();
    //   double l2 = magnetocrystalline_anisotropy_exchange_length();
    //   return std::min(l1,l2);
    // }

    // double magnetostatic_exchange_length() const
    //   {
    //     return std::sqrt( (2* exchange_constant() )
    //                       / (mu0() * saturation_magnetisation()
    //                          * saturation_magnetisation()));
    //   }

    // double magnetocrystalline_anisotropy_exchange_length() const
    //   {
    //     if (k1() > 0)
    //       return std::sqrt( exchange_constant() / k1() );
    //     else
    //       return 1e30; // "infinity" (e30 instead of higher in case we ever
    //                    // use floats not doubles).
    //     //??ds throw error instead?
    //   }


    // void output(std::ostream& stream) const
    // {
    //   stream << std::endl;
    //   stream << "Magnetic parameters are:" << std::endl;
    //   // stream << "Gamma = " << gamma() << std::endl;
    //   // stream << "Mu0 = " << mu0() << std::endl;
    //   stream << "Damping = " << gilbert_damping() << std::endl;
    //   // stream << "M_s = "  << saturation_magnetisation() << std::endl;
    //   // stream << "Exchange constant = " << exchange_constant() <<std::endl;
    //   // stream << "K_1 = " << k1() << std::endl;
    //   stream << "Easy_axis = " << easy_axis() << std::endl;

    //   stream << "magnetostatic debug = " << magnetostatic_debug_coeff() << std::endl;
    //   // stream << "exchange debug = " << exchange_debug_coeff() << std::endl;
    //   // stream << "boundary exchange debug = " << boundary_exchange_debug_coeff() << std::endl;
    //   // stream << "cryst anisotropy debug = " << ca_debug_coeff() << std::endl;
    //   // stream << std::endl;

    //   stream << "This gives the following normalised coeffs:" << std::endl;
    //   // stream << "Normalised_gamma = " <<    normalised_gamma() << std::endl;
    //   // stream << "Normalised_gilbert_damping = " <<  normalised_gilbert_damping() << std::endl;
    //   // stream << "Normalised_saturation_magnetisation = "
    //          // <<  normalised_saturation_magnetisation() << std::endl;
    //   stream << "Normalised_hex = " <<  normalised_hex() << std::endl;
    //   stream << "Normalised_hk = " <<  normalised_hk() << std::endl;
    //   stream << "Normalised_hms = " <<  normalised_hms() << std::endl;
    //   stream << std::endl;

    //   stream << "Exchange length for this material is: "
    //          << exchange_length() << std::endl;
    //   stream << "Your distance units are: " << distance_units()
    //          << "m" << std::endl;
    //   stream << "All elements should be smaller than the exchange length (and there should be code to check this...)." << std::endl;
    // }
  };


}

#endif
