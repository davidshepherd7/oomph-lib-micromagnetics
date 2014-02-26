#ifndef OOMPH_MAGNETOSTATICS_CALCULATOR_H
#define OOMPH_MAGNETOSTATICS_CALCULATOR_H

#include "../../src/generic/Vector.h"
#include "../../src/generic/oomph_utilities.h"

#include "micromag_types.h"

namespace oomph
{

  class MMInterpolator;
  template<unsigned VAL> class MMArrayInterpolator;
  class MagnetostaticFieldEquations;

  /// Class to calculate magnetostatic field
  class MagnetostaticsCalculator
  {
  public:
    virtual ~MagnetostaticsCalculator() {}

    virtual void get_magnetostatic_field(MMArrayInterpolator<5>* intp_pt,
                                         Vector<double> &H_ms) const=0;
    virtual void get_magnetostatic_field_time_derivative
    (MMInterpolator* intp_pt, Vector<double> &H_ms) const=0;
  };



  /// Class to calculate magnetostatic field using phi values in this element
  class ImplicitMagnetostaticsCalculator : public MagnetostaticsCalculator
  {
  public:

    ImplicitMagnetostaticsCalculator() {}

    virtual ~ImplicitMagnetostaticsCalculator() {}

    void get_magnetostatic_field(MMArrayInterpolator<5>* intp_pt,
                                 Vector<double> &H_ms) const;

    void get_magnetostatic_field_time_derivative
    (MMInterpolator* intp_pt, Vector<double> &H_ms) const;

    unsigned phi_index_micromag() const {return 0;}


  };



  /// Class to calculate magnetostatic field using phi values in this element
  class SemiImplicitMagnetostaticsCalculator : public MagnetostaticsCalculator
  {
  public:

    SemiImplicitMagnetostaticsCalculator() {}

    virtual ~SemiImplicitMagnetostaticsCalculator() {}

    void get_magnetostatic_field(MMArrayInterpolator<5>* intp_pt,
                                 Vector<double> &H_ms) const;

    void get_magnetostatic_field_time_derivative
    (MMInterpolator* intp_pt, Vector<double> &H_ms) const;

    /// \short Non-const access function for Magnetostatic_field_element_pt.
    MagnetostaticFieldEquations*& magnetostatic_field_element_pt()
    {return Magnetostatic_field_element_pt;}

    /// \short Const access function for Magnetostatic_field_element_pt.
    MagnetostaticFieldEquations* magnetostatic_field_element_pt() const
    {
#ifdef PARANOID
      if(Magnetostatic_field_element_pt == 0)
        {
          std::ostringstream error_msg;
          error_msg << "Magnetics element pointer not set.";
          throw OomphLibError(error_msg.str(), OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return Magnetostatic_field_element_pt;
    }

  private:

    MagnetostaticFieldEquations* Magnetostatic_field_element_pt;

  };


  /// Class to calculate magnetostatic field from function pt
  class AnalyticalMagnetostatics : public MagnetostaticsCalculator
  {
  public:

    AnalyticalMagnetostatics()
    {
      Magnetostatic_field_fct_pt = 0;
    }

    virtual ~AnalyticalMagnetostatics() {}

    void get_magnetostatic_field(MMArrayInterpolator<5>* intp_pt,
                                 Vector<double> &H_ms) const;

    void get_magnetostatic_field_time_derivative
    (MMInterpolator* intp_pt, Vector<double> &H_ms) const
    {
      //??ds dummy value
      H_ms.assign(3, -1000.0);
    }

    MagnetostaticFieldFctPt Magnetostatic_field_fct_pt;
  };

} // End of oomph namespace

#endif
