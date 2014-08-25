#ifndef OOMPH_BDF2_WITH_PREDICTOR_H
#define OOMPH_BDF2_WITH_PREDICTOR_H

#include "../../src/generic/timesteppers.h"

namespace oomph
{

  /// BDF2 time stepper with predictions by ebdf3 (ie like IMR).
  class BDF2Pred : public BDF<2>
  {
  public:

    BDF2Pred(bool adaptive=false) : BDF<2>(adaptive)
    {
      Predict_by_explicit_step = true;

      // Need one more history value for ebdf3
      Weight.resize(2, Weight.ncol()+1, 0.0);

      // Storing predicted values in slot after other information
      Predictor_storage_index = 2 + 3;
    }

    unsigned nprev_values() const {return 3;}
    unsigned ndt() const {return 3;}

    /// Dummy - just check that the values that
    /// problem::calculate_predicted_values() has been called right.
    void calculate_predicted_values(Data* const &data_pt)
    {
      if(adaptive_flag())
        {
          // Can't do it here, but we can check that the predicted values have
          // been updated.
          check_predicted_values_up_to_date();
        }
    }

    double temporal_error_in_value(Data* const &data_pt, const unsigned &i)
    {
      if(adaptive_flag())
        {
          // predicted value is more accurate so just compare with that
          //??ds is sign important? probably not...
          return data_pt->value(i) - data_pt->value(Predictor_storage_index, i);
        }
      else
        {
          std::string err("Tried to get the temporal error from a non-adaptive");
          err += " time stepper.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        }
    }

    void calculate_predicted_positions(Node* const &node_pt)
    {
      std::string err = "Not implemented";
      throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    double temporal_error_in_position(Node* const &node_pt, const unsigned &i)
    {
      std::string err = "Not implemented";
      throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  };

} // End of oomph namespace

#endif
