#ifndef OOMPH_TR_H
#define OOMPH_TR_H


#include "../../src/generic/timesteppers.h"

#include "llg_problem.h"

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

      Initial_derivative_set = false;
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

    void actions_before_timestep(Problem* problem_pt)
    {
      if(!Initial_derivative_set)
        {
          oomph_info << "Solving for derivative at initial time."
                     << " Warning: if residual is not in the correct form this may fail."
                     << std::endl;

          // Make all the implicit timesteppers steady
          unsigned n_time_steppers = problem_pt->ntime_stepper();
          std::vector<bool> was_steady(n_time_steppers);
          for(unsigned i=0;i<n_time_steppers;i++)
            {
              was_steady[i]=problem_pt->time_stepper_pt(i)->is_steady();
              problem_pt->time_stepper_pt(i)->make_steady();
            }

          // If it's an llg problem then we need some extra hacks so that BEM
          // doesn't interfere with get_jacobian and so that we are using the
          // explicit form. ??ds get rid of this somehow?
          LLGProblem* llg_pt = dynamic_cast<LLGProblem*>(problem_pt);
          bool needs_reset = false;
          if(llg_pt != 0)
            {
              llg_pt->Inside_explicit_timestep = true;
              if(llg_pt->Residual_calculator_pt->use_gilbert_form())
                {
                  llg_pt->Residual_calculator_pt->set_use_ll_form();
                  needs_reset = true;
                }
            }

          // Shift time backwards because we have already shifted it to t_1
          // when this function is called but we actually want the derivative
          // at t_0.
          double backup_time = time_pt()->time();
          time_pt()->time() -= time_pt()->dt();

          // Get the derivative at initial time and store in derivatives slot
          // ready for use in timestepping.
          DoubleVector f0;
          problem_pt->get_inverse_mass_matrix_times_residuals(f0);
          problem_pt->set_dofs(this->derivative_index(0), f0);

          // Revert time value
          time_pt()->time() = backup_time;

          // Revert llg settings
          if(llg_pt != 0)
            {
              llg_pt->Inside_explicit_timestep = false;
              if(needs_reset)
                {
                  llg_pt->Residual_calculator_pt->set_use_gilbert_form();
                  llg_pt->Inside_explicit_timestep = false;
                }
            }

          // Reset the is_steady status of all timesteppers that
          // weren't already steady when we came in here.
          for(unsigned i=0;i<n_time_steppers;i++)
            {
              if (!was_steady[i])
                {
                  problem_pt->time_stepper_pt(i)->undo_make_steady();
                }
            }

          Initial_derivative_set = true;
        }
    }


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


    bool Initial_derivative_set;

    ///\short This function advances the time history of the positions at a
    ///node. ??ds I don't really understand this part so there's a good
    ///chance it will break if you try to do problems with moving nodes...
    void shift_time_positions(Node* const &node_pt)
    {
      unsigned n_dim = node_pt->ndim();
      unsigned n_position_type = node_pt->nposition_type();

      // Previous step dt, time has already been shifted so it's in slot 1.
      double dtn = time_pt()->dt(1);

      //Loop over the position directions
      for(unsigned i=0;i<n_dim;i++)
        {
          //If the position is not a copy
          if(node_pt->position_is_a_copy(i) == false)
            {
              //Loop over the position types
              for(unsigned k=0;k<n_position_type;k++)
                {
                  // Calculate velocity at step n
                  double dxnm1 = node_pt->x_gen(derivative_index(0), k, i);
                  double xnm1 = node_pt->x_gen(1, k, i);
                  double xn = node_pt->x_gen(0, k, i);
                  double dxn = (2/dtn)*(xn - xnm1) - dxnm1;

                  node_pt->x_gen(derivative_index(0), k, i) = dxn;

                  // Shift velocity
                  node_pt->x_gen(derivative_index(1), k, i) = dxnm1;

                  // Shift the stored position
                  node_pt->x_gen(1, k, i) = node_pt->x_gen(0,k,i);
                }
            }

        }
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
