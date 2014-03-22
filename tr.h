#ifndef OOMPH_TR_H
#define OOMPH_TR_H


#include "../../src/generic/timesteppers.h"

namespace oomph
{

class TR : public TimeStepper
{
public:
  /// Constructor, storage for two past derivatives, one past value and
  /// present value.
  TR(bool adaptive=false) : TimeStepper(2+2, 1)
  {
    //Set the weight for the zero-th derivative
    Weight(0,0) = 1.0;

    Is_first_step=true;
  }

  /// Virtual destructor
  virtual ~TR() {}

  ///Return the actual order of the scheme
  unsigned order() {return 2;}

  /// Set the weights
  void set_weights()
  {
    double dt = Time_pt->dt(0);
    Weight(1,0) = 2.0/dt;
    Weight(1,1) = -2.0/dt;
    Weight(1,derivative_index(0)) = -1.0; // old derivative
  }

  /// Number of previous values available.
  unsigned nprev_values() {return 1;}

  /// Location in data of derivatives
  unsigned derivative_index(const unsigned& t) const
  {
#ifdef PARANOID
    if(t >= 2)
      {
        std::string err = "Only storing two derivatives.";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif
    return t + 2;
  }

  /// Number of timestep increments that need to be stored by the scheme
  unsigned ndt() {return 2;}

  /// \short Initialise the time-history for the Data values, corresponding
  /// to an impulsive start.
  void assign_initial_values_impulsive(Data* const &data_pt)
  {
    throw OomphLibError("Function not yet implemented",
                        OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
  }

  /// \short Initialise the time-history for the nodal positions
  /// corresponding to an impulsive start.
  void assign_initial_positions_impulsive(Node* const &node_pt)
  {
    throw OomphLibError("Function not yet implemented",
                        OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
  }

  void actions_after_timestep(Problem* problem_pt) {}

  void actions_before_timestep(Problem* problem_pt) {}


  /// \short This function updates the Data's time history so that
  /// we can advance to the next timestep. As for BDF schemes,
  /// we simply push the values backwards...
  void shift_time_values(Data* const &data_pt)
  {
    //Find number of values stored
    unsigned n_value = data_pt->nvalue();

    // Previous step dt, time has already been shifted so it's in slot 1.
    double dtn = time_pt()->dt(1);

    //Loop over the values
    for(unsigned j=0; j<n_value; j++)
      {
        //Set previous values to the previous value, if not a copy
        if(data_pt->is_a_copy(j) == false)
          {
            if(!Is_first_step)
              {

                // Calculate time derivative at step n
                double fnm1 = data_pt->value(derivative_index(0), j);
                double ynm1 = data_pt->value(1, j);
                double yn = data_pt->value(0, j);
                double fn = (2/dtn)*(yn - ynm1) - fnm1;

                data_pt->set_value(derivative_index(0), j, fn);

                // Shift time derivative
                data_pt->set_value(derivative_index(1), j, fnm1);
              }

            // Shift value
            data_pt->set_value(1, j,
                               data_pt->value(0, j));

          }
      }

    Is_first_step = false;
  }


  bool Is_first_step;

  ///\short This function advances the time history of the positions
  ///at a node.
  void shift_time_positions(Node* const &node_pt)
  {
    throw OomphLibError("Function not yet implemented",
                        OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
  }


private:
  /// Broken copy constructor
  TR(const TR& dummy)
  {BrokenCopy::broken_copy("TR");}

  /// Broken assignment operator
  void operator=(const TR& dummy)
  {BrokenCopy::broken_assign("TR");}

};
} // End of oomph namespace

#endif
