#ifndef OOMPH_ARRAY_INTERPOLATOR_H
#define OOMPH_ARRAY_INTERPOLATOR_H

#include "../../src/generic/Vector.h"
#include "../../src/generic/elements.h"
#include "../../src/generic/shape.h"

namespace oomph
{
  using namespace InterpolatorHelpers;

  /* Implementation notes:

   * Use memoization to ensure we don't duplicate calculations.

   * Lots of code is repeated because c++ function pointers suck

   * Use Vector<Vector<double> > rather than matrix so that we can easily
   (and for "free") return vectors of a single row (i.e. all derivatives
   of one value).

   * No need to template by DIM or Nnode: it has no effect on speed (tested).


   TODO:

   * We really have two classes here: a memoiser and an interpolator.  Split
   them. Memoising interface over an interpolator is probably best, then we
   can swap the interpolator out without breaking any memoising.

   * Pass in parameter for time deriv order to interpolation (instead of
   separate interpolate_dxxxxdt functions).

  */

  // ============================================================
  /// Interpolator that auto detects what to do (possibly slow?).
  // ============================================================
  template<unsigned VAL>
  class GeneralArrayInterpolator
  {

  public:
    /// Default constructor
    GeneralArrayInterpolator(const FiniteElement* const this_element,
                             const Vector<double> &s,
                             bool call_build=true,
                             const unsigned& time_index=0)
      :
      // Cache some loop end conditions
      Nnode(this_element->nnode()),
      Dim(this_element->dim()),
      Nprev_value_derivative(this_element->node_pt(0)->
                             time_stepper_pt()->ntstorage()),
      Nprev_value_pos_derivative(this_element->node_pt(0)->
                                 position_time_stepper_pt()->ntstorage()),

    // Initialise pointers
      This_element(this_element),
      Ts_weights_pt(this_element->node_pt(0)->time_stepper_pt()->weights_pt()),
      Position_ts_weights_pt(this_element->node_pt(0)->position_time_stepper_pt()
                             ->weights_pt()),

    // Initialise storage for shape functions
      Psi(Nnode),
      Test(Nnode),
      Dpsidx(Nnode, Dim),
      Dtestdx(Nnode, Dim),

    // Copy local position
      S(s),

    // Negative time to signify that it has not been calculated yet
      Intp_time(NotYetCalculatedValue),

    // Time index to interpolate values at (0 = present)
      Time_index(time_index)
    {
      Ts_pt = this_element->node_pt(0)->time_stepper_pt();

      // Check for old implementation of IMR: can't handle it in this
      // interpolator.
#ifdef PARANOID
      if(Ts_pt->nprev_values_for_value_at_evaluation_time() != 1)
        {
          std::string err = "Can't handle cases which require history values in value interpolation.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      // Initialise storage. ??ds we could avoid this initialisation and
      // the use of a NotYetCalculatedValue by storing pointers to the
      // values and pointers to the storage separately, then when the value
      // is calculated the pointer to the value is set to point to the
      // storage. In this way storage is statically allocated but not
      // initialised (for speed), and we still get the benefits of
      // memoisation.
      for(unsigned i=0; i<3; i++) X[i] = NotYetCalculatedValue;
      for(unsigned i=0; i<3; i++) Dxdt[i] = NotYetCalculatedValue;
      for(unsigned i=0; i<VAL; i++) Values[i] = NotYetCalculatedValue;
      for(unsigned i=0; i<VAL; i++) Dvaluesdt[i] = NotYetCalculatedValue;
      for(unsigned i=0; i<VAL; i++)
        {
          for(unsigned j=0; j<3; j++) Dvaluesdx[i][j] = NotYetCalculatedValue;
        }

      if(call_build)
        {
          build();
        }
    }

    virtual void build()
      {
        // Set up shape + test functions
        J = This_element->dshape_eulerian(s(), Psi, Dpsidx);
        Test = Psi;
        Dtestdx = Dpsidx;

        // Find out if any nodes in this element are hanging
        Has_hanging_nodes = This_element->has_hanging_nodes();

        // If paranoid check that all nodes have the same nvalues.
        InterpolatorHelpers::check_nvalues_in_element_same(This_element);
      }

    /// Destructor
    ~GeneralArrayInterpolator() {}

    double time()
    {
      if(uninitialised(Intp_time)) {Intp_time = interpolate_time();}
      return Intp_time;
    }

    TimeStepper* ts_pt() const {return Ts_pt;}

    unsigned dim() const {return Dim;}

    double x(const unsigned &i)
    {return x()[i];}

    const double* x()
    {
      if(uninitialised(X)) interpolate_x();
      return X;
    }

    double dxdt(const unsigned &i)
    {return dxdt()[i];}

    const double* dxdt()
    {
      if(uninitialised(Dxdt)) interpolate_dxdt();
      return Dxdt;
    }

    double value(const unsigned &i_val)
    {return value()[i_val];}

    const double* value()
    {
      if(uninitialised(Values)) interpolate_values(0, VAL);
      return Values;
    }

    double dvaluedt(const unsigned &i_val)
    {return dvaluedt()[i_val];}

    const double* dvaluedt()
    {
      if(uninitialised(Dvaluesdt)) interpolate_dvaluesdt(0, VAL);
      return Dvaluesdt;
    }

    double dvaluedx(const unsigned &i_val, const unsigned &direction)
    {return dvaluedx(i_val)[direction];}


    const double* dvaluedx(const unsigned &i_val)
    {
      if(uninitialised(Dvaluesdx[i_val]))
        {
          interpolate_dvaluesdx(i_val);
        }
      return Dvaluesdx[i_val];
    }

    // Access functions for Jacobian and shape/test functions
    double j() const {return J;}
    const Vector<double>& s() const {return S;}
    double psi(const unsigned &i) const {return Psi(i);}
    double test(const unsigned &i) const {return Test(i);}
    double dpsidx(const unsigned &i, const unsigned &i_direc) const
    {return Dpsidx(i, i_direc);}
    double dtestdx(const unsigned &i, const unsigned &i_direc) const
    {return Dtestdx(i, i_direc);}


    /// Interpolate evaluation position
    void interpolate_x()
    {
      if(!Has_hanging_nodes) {interpolate_x_raw();}
      else {interpolate_x_with_hanging_nodes();}
    }

    /// Interpolate derivative of position
    void interpolate_dxdt()
    {
      if(!Has_hanging_nodes) {interpolate_dxdt_raw();}
      else {interpolate_dxdt_with_hanging_nodes();}

    }

    /// Interpolate evaluation time
    double interpolate_time()
    {
      //??ds get rid of?
      double time = 0.0;
      Time* time_pt = This_element->node_pt(0)->time_stepper_pt()->time_pt();
      time += time_pt->time();
      return time;
    }

    // Interpolate values (by default get all, optionally just get some range.
    void interpolate_values(const unsigned &start,
                            const unsigned &end)
    {
      if(!Has_hanging_nodes) {interpolate_values_raw(start, end);}
      else {interpolate_values_with_hanging_nodes(start, end);}
    }

    void interpolate_dvaluesdt(const unsigned &start,
                               const unsigned &end)
    {
      if(!Has_hanging_nodes) {interpolate_dvaluesdt_raw(start, end);}
      else {interpolate_dvaluesdt_with_hanging_nodes(start, end);}
    }

    void interpolate_dvaluesdx(const unsigned &i_val)
    {
      if(!Has_hanging_nodes) {interpolate_dvaluesdx_raw(i_val);}
      else {interpolate_dvaluesdx_with_hanging_nodes(i_val);}
    }

  protected:

    // Loop end conditions etc.
    const unsigned Nnode;
    const unsigned Dim;
    const unsigned Nprev_value_derivative;
    const unsigned Nprev_value_pos_derivative;

    bool Has_hanging_nodes;

    // Cached pointers
    const FiniteElement* This_element;
    const DenseMatrix<double>* Ts_weights_pt;
    const DenseMatrix<double>* Position_ts_weights_pt;
    TimeStepper* Ts_pt;

    // Jacobian + shape/test functions
    double J;
    Shape Psi;
    Shape Test;
    DShape Dpsidx;
    DShape Dtestdx;

    // Local coordinate
    Vector<double> S;

    // Interpolated value storage (note we can't name the time variable
    // "Time" because that is used for the time class).
    double Intp_time;

    /// Time index to interpolate values at, 0=current. If doing history
    /// interpolation then we cannot calculate time derivatives.
    const unsigned Time_index;

    /// Storage for positions
    double X[3];

    /// Storage for position derivatives
    double Dxdt[3];

    /// Storage for values
    double Values[VAL];

    /// Storage for value time derivatives
    double Dvaluesdt[VAL];

    /// Storage for value space derivatives
    double Dvaluesdx[VAL][3];


    // Position interpolation
    // ============================================================

    void interpolate_x_with_hanging_nodes()
    {
      for(unsigned j=0; j<Dim; j++)
        {
          X[j] = 0.0;
          for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
            {
              X[j] += This_element->nodal_position(Time_index, i_nd, j) * Psi(i_nd);
            }
        }
    }

    /// \short Use if no nodes are hanging (raw access functions are faster).
    void interpolate_x_raw()
    {
      for(unsigned j=0; j<Dim; j++)
        {
          X[j] = 0.0;
          for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
            {
              X[j] +=
                This_element->raw_nodal_position(Time_index, i_nd, j) * Psi(i_nd);
            }
        }
    }

    // Position derivatives interpolation
    // ============================================================

    void interpolate_dxdt_raw()
    {
      check_can_calculate_time_derivatives();

      for(unsigned j=0;j<Dim;j++)
        {
          Dxdt[j] = 0;
          for(unsigned l=0;l<Nnode;l++)
            {
              for(unsigned t=0;t<Nprev_value_pos_derivative;t++)
                {
                  Dxdt[j] += (*Position_ts_weights_pt)(1,t)
                    * This_element->raw_nodal_position(t, l, j) * Psi(l);
                }
            }
        }
    }

    void interpolate_dxdt_with_hanging_nodes()
    {
      check_can_calculate_time_derivatives();

      for(unsigned j=0;j<Dim;j++)
        {
          Dxdt[j] = 0;
          for(unsigned l=0;l<Nnode;l++)
            {
              for(unsigned t=0;t<Nprev_value_pos_derivative;t++)
                {
                  Dxdt[j] += (*Position_ts_weights_pt)(1,t)
                    * This_element->nodal_position(t, l, j) * Psi(l);
                }
            }
        }
    }


    // Interpolate values
    // ============================================================

    //??ds copy this to all function's doc

    /// \short Interpolate [something] [with/without] considering hanging nodes (it
    /// is faster to use "raw" access functions which ignore hanging nodes
    /// where possible). Interpolates a range of values [start, end)
    /// (i.e. including start but not end).

    void interpolate_values_raw(const unsigned &start,
                                const unsigned &end)
    {
      for(unsigned j=start; j<end; j++)
        {
          Values[j] = 0;
          for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
            {
              Values[j] +=
                This_element->raw_nodal_value(Time_index, i_nd, j) * Psi(i_nd);
            }
        }
    }

    void interpolate_values_with_hanging_nodes(const unsigned &start,
                                               const unsigned &end)
    {
      for(unsigned j=start; j<end; j++)
        {
          Values[j] = 0;
          for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
            {
              Values[j] +=
                This_element->nodal_value(Time_index, i_nd, j) * Psi(i_nd);
            }
        }
    }


    // Interpolate derivatives of values w.r.t. time
    // ============================================================

    void interpolate_dvaluesdt_raw(const unsigned &start,
                                   const unsigned &end)
    {
      check_can_calculate_time_derivatives();

      for(unsigned j=start; j<end; j++)
        {
          Dvaluesdt[j] = 0;
          for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
            {
              for(unsigned i_tm=0; i_tm<Nprev_value_derivative; i_tm++)
                {
                  Dvaluesdt[j] += This_element->raw_nodal_value(i_tm, i_nd, j)
                    * Psi(i_nd) * (*Ts_weights_pt)(1, i_tm);
                }
            }
        }
    }

    void interpolate_dvaluesdt_with_hanging_nodes(const unsigned &start,
                                                  const unsigned &end)
    {
      check_can_calculate_time_derivatives();

      for(unsigned j=start; j<end; j++)
        {
          Dvaluesdt[j] = 0;
          for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
            {
              for(unsigned i_tm=0; i_tm<Nprev_value_derivative; i_tm++)
                {
                  Dvaluesdt[j] += This_element->nodal_value(i_tm, i_nd, j)
                    * Psi(i_nd) * (*Ts_weights_pt)(1, i_tm);
                }
            }
        }
    }


    // Interpolate derivatives of values w.r.t. position
    // ============================================================


    void interpolate_dvaluesdx_raw(const unsigned &i_value)
    {
      for(unsigned i_direc=0; i_direc<Dim; i_direc++)
        {
          Dvaluesdx[i_value][i_direc] = 0;
          for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
            {
              Dvaluesdx[i_value][i_direc] +=
                This_element->raw_nodal_value(Time_index, i_nd, i_value)
                * Dpsidx(i_nd, i_direc);
            }
        }
    }

    void interpolate_dvaluesdx_with_hanging_nodes(const unsigned &i_value)
    {

      for(unsigned i_direc=0; i_direc<Dim; i_direc++)
        {
          Dvaluesdx[i_value][i_direc] = 0;
          for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
            {
              Dvaluesdx[i_value][i_direc] +=
                This_element->nodal_value(Time_index, i_nd, i_value)
                * Dpsidx(i_nd, i_direc);
            }
        }
    }


    void check_can_calculate_time_derivatives()
    {
#ifdef PARANOID
      if(Time_index != 0)
        {
          std::string err = "Can only calculate time derivatives at present time (because we can't assume that the older history values are what is needed for time derivative calculations).";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
    }


  };



  /// Modified version of the itnerpolator for face elements: 1) can't get x
  /// derivatives very easily here so throw errors if we try. 2) Don't try to
  /// get dpsi etc (because we can't), just get shape and Jacobian separately.
  template<unsigned VAL>
  class FaceElementArrayInterpolator
    : public GeneralArrayInterpolator<VAL>
  {

  public:

    /// Constructor, jsut use underlying interpolator
    FaceElementArrayInterpolator(const FaceElement* this_element,
                                   const Vector<double> &s)
      : GeneralArrayInterpolator<VAL>(this_element, s, false)
    {
      build();
    }


    virtual void build()
    {
      // Set up shape + test functions
      this->J = this->This_element->J_eulerian(this->s());
      this->This_element->shape(this->s(), this->Psi);
      this->Test = this->Psi;

      // Find out if any nodes in this element are hanging
      this->Has_hanging_nodes = this->This_element->has_hanging_nodes();

      // If paranoid check that all nodes have the same nvalues.
      InterpolatorHelpers::check_nvalues_in_element_same(this->This_element);
    }

    double dpsidx(const unsigned &i, const unsigned &i_direc) const
    {
      broken_function_error();
      return 0;
    }

    double dtestdx(const unsigned &i, const unsigned &i_direc) const
    {
      broken_function_error();
      return 0;
    }

  private:

    void broken_function_error() const
    {
      std::string err = "Cannot calculate derivatives w.r.t. x in face elements";
      throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                          OOMPH_CURRENT_FUNCTION);
    }

    // /// Broken constructors
    // // FaceElementArrayInterpolator() {}
    // FaceElementArrayInterpolator(FaceElementArrayInterpolator& d) {}
    // void operator=(FaceElementArrayInterpolator& d) {}

  };

} // End of oomph namespace

#endif
