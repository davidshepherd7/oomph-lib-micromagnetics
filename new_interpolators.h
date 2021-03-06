#ifndef OOMPH_NEW_INTERPOLATORS_H
#define OOMPH_NEW_INTERPOLATORS_H

#include "interpolator.h"

namespace oomph
{

  namespace InterpolatorHelpers
  {
    inline void check_no_history_values_in_current_value(const TimeStepper* ts_pt)
    {
#ifdef PARANOID
      if(ts_pt->nprev_values_for_value_at_evaluation_time() != 1)
        {
          std::string err = "Can't handle cases which require history values";
          err += " in value interpolation.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
    }

    inline void check_position_timestepper_matches_main(const TimeStepper* ts_pt,
                                                        const TimeStepper* pos_ts_pt)
    {
#ifdef PARANOID
      if(ts_pt != pos_ts_pt)
        {
          std::string err = "Not implemented for different position time stepper.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
          // wouldn't be hard to implement, just a bit more messy in places
          // and more code that I don't need.
        }
#endif
    }

  }


  /// A generic cached datatype
  template<class T> class Cached
  {
  public:
    Cached() {Initialised = false;}

    bool is_initialised() const {return Initialised;}

    // Assume that we initialise everything at once
    T& set()
    {
      Initialised = true;
      return Data;
    }

    const T& get() const
    {
#ifdef PARANOID
      if(!is_initialised())
        {
          std::string err = "Got uninitialised value!";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return Data;
    }

    void reset() {Initialised = false;}

  private:
    T Data;
    bool Initialised;
  };


  /// A cached array type (separate from other because we can't use
  /// std::array (c++11) so array isn't really a type and we can't template
  /// by it :( ).
  template<unsigned VAL> class CachedArray
  {
  public:
    CachedArray() {Initialised = false;}

    bool is_initialised() const {return Initialised;}

    // Assume that we initialise everything at once
    double* set()
    {
      Initialised = true;
      return Data_pt;
    }

    const double* get() const
    {
#ifdef PARANOID
      if(!is_initialised())
        {
          std::string err = "Got uninitialised value!";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return Data_pt;
    }

    void reset() {Initialised = false;}

  private:
    double Data_pt[VAL];
    bool Initialised;
  };

  /// Virtual base class for interpolators: class to calculate
  /// values/position/derivatives at a given point in an element.
  class InterpolatorBase
  {
  public:

    InterpolatorBase()
    {
      This_element = 0;
      Ts_weights_pt = 0;
      Ts_pt = 0;
    }

    virtual ~InterpolatorBase() {}

    /// Cache some values from element, set up test functions and check
    /// that the interpolator will work for the given element.
    virtual void build(const Vector<double>& s)
    {
      // Check we have the right pointers
#ifdef PARANOID
      if(This_element == 0)
        {
          std::string err = "This_element pointer is null! Need to set it before building.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      if(Ts_pt == 0)
        {
          std::string err = "Ts_pt is null need to set it before building!";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      // get pointers
      Ts_weights_pt = Ts_pt->weights_pt();

      // get constants
      Dim = This_element->nodal_dimension();
      Nnode = This_element->nnode();
      Nprev_value_derivative = Ts_pt->ntstorage();
      S = s;

      // Allocate shape function memory if needed
      if(Psi.nindex1() != Nnode || Dpsidx.nindex2() != Dim)
        {
          Psi.resize(Nnode);
          Test.resize(Nnode);
          Dpsidx.resize(Nnode, Dim);
          Dtestdx.resize(Nnode, Dim);
        }

      // Calculate shape/test/derivatives as possible according to geometry.
      if(This_element->nodal_dimension() == This_element->dim())
        {
          // Set up shape + test + derivative functions for this position
          J = This_element->dshape_eulerian(s, Psi, Dpsidx);
          Test = Psi;
          Dtestdx = Dpsidx;

          Can_calculate_x_derivatives = true;
        }
      else
        {
          // Set up shape + test functions for this position
          This_element->shape(s, Psi);
          J = This_element->J_eulerian(s);
          Test = Psi;

          Can_calculate_x_derivatives = false;
        }

      // Paranoid checks for this interpolator
      check_interpolator_applicable();
    }

    virtual void check_interpolator_applicable() const=0;

    /// Interpolation functions
    virtual double interpolate_time() const=0;
    virtual double interpolate_x(const unsigned& j) const=0;
    virtual double interpolate_dxdt(const unsigned& j) const=0;
    virtual double interpolate_value(const unsigned& j) const=0;
    virtual double interpolate_dvaluedt(const unsigned& j) const=0;
    virtual double interpolate_dvaluedx(const unsigned& i_val,
                                        const unsigned& i_direc) const=0;

    const DShape& dpsidx() const
    {
#ifdef PARANOID
      if(!Can_calculate_x_derivatives)
        {
          std::string err = "Can't get x derivatives for this geometry.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return Dpsidx;
    }

    const DShape& dtestdx() const
    {
#ifdef PARANOID
      if(!Can_calculate_x_derivatives)
        {
          std::string err = "Can't get x derivatives for this geometry.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return Dtestdx;
    }


    ///
    unsigned Dim;
    unsigned Nnode;
    unsigned Nprev_value_derivative;

    const FiniteElement* This_element;
    const TimeStepper* Ts_pt;
    const DenseMatrix<double>* Ts_weights_pt;

    // Jacobian + shape/test functions
    double J;
    Shape Psi;
    Shape Test;
    Vector<double> S;

    private:
    bool Can_calculate_x_derivatives;
    DShape Dpsidx;
    DShape Dtestdx;

  };


  /// Simplest interpolator implementation: no hanging nodes, no position
  /// time stepping, no history interpolation, no time-interpolation of
  /// current value (i.e midpoint rule).
  class RawInterpolator : public InterpolatorBase
  {
  public:

    virtual void check_interpolator_applicable() const
    {
      using namespace InterpolatorHelpers;

      // If paranoid check that all nodes have the same nvalues.
      check_nvalues_in_element_same(This_element);

      // If paranoid check that all node's time steppers are the same
      check_timesteppers_same(This_element);

      TimeStepper* pos_ts_pt = This_element->node_pt(0)->position_time_stepper_pt();
      check_position_timestepper_matches_main(Ts_pt, pos_ts_pt);

      // If paranoid check that we don't need history values to evaluate
      // current time value.
      check_no_history_values_in_current_value(This_element->node_pt(0)->time_stepper_pt());
    }


    double interpolate_time() const
    {
      return Ts_pt->time();
    }

    double interpolate_x(const unsigned& j) const
    {
      double x = 0;
      for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
        {
          x += This_element->raw_nodal_position(0, i_nd, j) * Psi(i_nd);
        }
      return x;
    }

    double interpolate_dxdt(const unsigned& j) const
    {
      double dxdt = 0;
      for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
        {
          for(unsigned i_tm=0; i_tm<Nprev_value_derivative; i_tm++)
            {
              dxdt += This_element->raw_nodal_position(i_tm, i_nd, j)
                * Psi(i_nd) * (*Ts_weights_pt)(1, i_tm);
            }
        }
      return dxdt;
    }

    double interpolate_value(const unsigned& j) const
    {
      double value = 0;
      for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
        {
          value += This_element->raw_nodal_value(0, i_nd, j)*Psi(i_nd);
        }
      return value;
    }

    double interpolate_dvaluedt(const unsigned& j) const
    {
      double dvaluedt = 0;
      for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
        {
          for(unsigned i_tm=0; i_tm<Nprev_value_derivative; i_tm++)
            {
              dvaluedt += This_element->raw_nodal_value(i_tm, i_nd, j)
                * Psi(i_nd) * (*Ts_weights_pt)(1, i_tm);
            }
        }
      return dvaluedt;
    }

    double interpolate_dvaluedx(const unsigned &i_val,
                                const unsigned& i_direc) const
    {
      double dvaluedx = 0;
      for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
        {
          dvaluedx += This_element->raw_nodal_value(0, i_nd, i_val)
            * dpsidx()(i_nd, i_direc);
        }
      return dvaluedx;
    }
  };


  /// Midpoint rule interpolator: no hanging nodes, no position time
  /// stepping, no history interpolation. But can handle time-interpolation
  /// of current value (i.e midpoint rule).
  class RawTimeInterpolatedValueInterpolator : public InterpolatorBase
  {
    // Maybe this should just inherit from RawInterpolator... mostly the same.
  public:

    unsigned Nprev_value_current_value;

    virtual void build(const Vector<double>& s)
    {
      InterpolatorBase::build(s);

      Nprev_value_current_value = Ts_pt->nprev_values_for_value_at_evaluation_time();
    }

    virtual void check_interpolator_applicable() const
    {
      using namespace InterpolatorHelpers;

      // If paranoid check that all nodes have the same nvalues.
      check_nvalues_in_element_same(This_element);

      // If paranoid check that all node's time steppers are the same
      check_timesteppers_same(This_element);

      TimeStepper* pos_ts_pt = This_element->node_pt(0)->position_time_stepper_pt();
      check_position_timestepper_matches_main(Ts_pt, pos_ts_pt);
    }


    double interpolate_time() const
    {
      double time = 0.0;
      for(unsigned i_tm=0; i_tm<Nprev_value_current_value; i_tm++)
        {
          time += Ts_pt->time_pt()->time(i_tm)
            * (*Ts_weights_pt) (0, i_tm);
        }
      return time;
    }

    double interpolate_x(const unsigned& j) const
    {
      double x = 0;
      for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
        {
          for(unsigned i_tm=0; i_tm<Nprev_value_current_value; i_tm++)
            {
              x += This_element->raw_nodal_position(i_tm, i_nd, j)
                * Psi(i_nd)
                * (*Ts_weights_pt) (0, i_tm);
            }
        }
      return x;
    }

    double interpolate_dxdt(const unsigned& j) const
    {
      double dxdt = 0;
      for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
        {
          for(unsigned i_tm=0; i_tm<Nprev_value_derivative; i_tm++)
            {
              dxdt += This_element->raw_nodal_position(i_tm, i_nd, j)
                * Psi(i_nd)
                * (*Ts_weights_pt)(1, i_tm);
            }
        }
      return dxdt;
    }

    double interpolate_value(const unsigned& j) const
    {
      double value = 0;
      for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
        {
          for(unsigned i_tm=0; i_tm<Nprev_value_current_value; i_tm++)
            {
              value += This_element->raw_nodal_value(i_tm, i_nd, j)
                * Psi(i_nd)
                * (*Ts_weights_pt) (0, i_tm);
            }
        }
      return value;
    }

    double interpolate_dvaluedt(const unsigned& j) const
    {
      double dvaluedt = 0;
      for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
        {
          for(unsigned i_tm=0; i_tm<Nprev_value_derivative; i_tm++)
            {
              dvaluedt += This_element->raw_nodal_value(i_tm, i_nd, j)
                * Psi(i_nd)
                * (*Ts_weights_pt)(1, i_tm);
            }
        }
      return dvaluedt;
    }

    double interpolate_dvaluedx(const unsigned &i_val,
                                const unsigned& i_direc) const
    {
      double dvaluedx = 0;
      for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
        {
          for(unsigned i_tm=0; i_tm<Nprev_value_current_value; i_tm++)
            {
              dvaluedx += This_element->raw_nodal_value(i_tm, i_nd, i_val)
                * dpsidx()(i_nd, i_direc)
                * (*Ts_weights_pt) (0, i_tm);
            }
        }
      return dvaluedx;
    }
  };


  /// Interpolator to interpolate from history values.
  class HistoryInterpolator : public InterpolatorBase
  {
  public:

    unsigned Time_index;

    virtual void check_interpolator_applicable() const
    {
      using namespace InterpolatorHelpers;

      // If paranoid check that all nodes have the same nvalues.
      check_nvalues_in_element_same(This_element);

      // If paranoid check that all node's time steppers are the same
      check_timesteppers_same(This_element);

      TimeStepper* pos_ts_pt = This_element->node_pt(0)->position_time_stepper_pt();
      check_position_timestepper_matches_main(Ts_pt, pos_ts_pt);

      // If paranoid check that we don't need history values to evaluate
      // current time value.
      check_no_history_values_in_current_value(This_element->node_pt(0)->time_stepper_pt());
    }


    double interpolate_time() const
    {
      return Ts_pt->time_pt()->time(Time_index);
    }

    double interpolate_x(const unsigned& j) const
    {
      double x = 0;
      for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
        {
          x += This_element->raw_nodal_position(Time_index, i_nd, j) * Psi(i_nd);
        }
      return x;
    }

    double interpolate_dxdt(const unsigned& j) const
    {
      std::string err = "Can't calculate history time derivatives.";
      throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    double interpolate_value(const unsigned& j) const
    {
      double value = 0;
      for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
        {
          value += This_element->raw_nodal_value(Time_index, i_nd, j)*Psi(i_nd);
        }
      return value;
    }

    double interpolate_dvaluedt(const unsigned& j) const
    {
      std::string err = "Can't calculate history time derivatives.";
      throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    double interpolate_dvaluedx(const unsigned &i_val,
                                const unsigned& i_direc) const
    {
      double dvaluedx = 0;
      for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
        {
          dvaluedx += This_element->raw_nodal_value(Time_index, i_nd, i_val)
            * dpsidx()(i_nd, i_direc);
        }
      return dvaluedx;
    }
  };

  /// Base class for a caching interpolator. Gives interface functions
  /// common to all caching interpolators no matter what array storage they
  /// are using.
  class CachingInterpolatorBase
  {
  public:

    CachingInterpolatorBase() : Intp_pt(0) {}

    virtual ~CachingInterpolatorBase()
    {
      delete Intp_pt; Intp_pt = 0;
    }

    virtual void build(const Vector<double>& s)
    {
      Intp_pt->build(s);
    }

    /// Get time value
    double time()
    {
      if(!Intp_time.is_initialised())
        {
          Intp_time.set() = Intp_pt->interpolate_time();
        }

      return Intp_time.get();
    }


    // Access functions for Jacobian and shape/test functions from
    // interpolator itself.
    double j() const {return Intp_pt->J;}
    const Vector<double>& s() const {return Intp_pt->S;}
    double psi(const unsigned &i) const {return Intp_pt->Psi(i);}
    double test(const unsigned &i) const {return Intp_pt->Test(i);}
    double dpsidx(const unsigned &i, const unsigned &i_direc) const
    {return Intp_pt->dpsidx()(i, i_direc);}
    double dtestdx(const unsigned &i, const unsigned &i_direc) const
    {return Intp_pt->dtestdx()(i, i_direc);}

    unsigned dim() const { return Intp_pt->Dim; }
    const TimeStepper* ts_pt() const { return Intp_pt->Ts_pt; }

    InterpolatorBase* Intp_pt;

  protected:

    Cached<double> Intp_time;

  };

  /// Base class for caching interpolators using c-arrays. Calculates
  /// interpolated position and derivatives using Intp_pt and stores them
  /// in arrays. Provides an interface which allows the user to ignore
  /// whether the value was cached or not. With c++11 we can probably ditch
  /// this and use std::array<double> in CachingArrayInterpolator.
  class CachingArrayInterpolatorBase : public CachingInterpolatorBase
  {
  public:

    virtual void build(const Vector<double>& s)
    {
      CachingInterpolatorBase::build(s);

      // Clear memoised values
      X.reset();
      Dxdt.reset();
      Intp_time.reset();
    }

    /// Get x values
    const double* x()
    {
      if(!X.is_initialised())
        {
          for(unsigned j=0; j<Intp_pt->Dim; j++)
            {
              X.set()[j] = Intp_pt->interpolate_x(j);
            }
        }

      return X.get();
    }

    /// Get dxdt values
    const double* dxdt()
    {
      if(!Dxdt.is_initialised())
        {
          for(unsigned j=0; j<Intp_pt->Dim; j++)
            {
              Dxdt.set()[j] = Intp_pt->interpolate_dxdt(j);
            }
        }

      return Dxdt.get();
    }

  protected:
    CachedArray<3> X;
    CachedArray<3> Dxdt;

  };



  class CachingVectorInterpolatorBase : public CachingInterpolatorBase
  {
  public:

    virtual void build(const Vector<double>& s)
    {
      CachingInterpolatorBase::build(s);

      // Clear memoised values
      X.reset();
      Dxdt.reset();
      Intp_time.reset();
    }

    /// Get x values
    const Vector<double>& x()
    {
      if(!X.is_initialised())
        {
          X.set().assign(Intp_pt->Dim, 0.0);
          for(unsigned j=0; j<Intp_pt->Dim; j++)
            {
              X.set()[j] = Intp_pt->interpolate_x(j);
            }
        }

      return X.get();
    }

    /// Get dxdt values
    const Vector<double>& dxdt()
    {
      if(!Dxdt.is_initialised())
        {
          Dxdt.set().assign(Intp_pt->Dim, 0.0);
          for(unsigned j=0; j<Intp_pt->Dim; j++)
            {
              Dxdt.set()[j] = Intp_pt->interpolate_dxdt(j);
            }
        }

      return Dxdt.get();
    }

  protected:
    Cached<Vector<double> > X;
    Cached<Vector<double> > Dxdt;
  };


  /// A generic caching interpolator which stores values unlabelled in an
  /// array (whose size must be specified by the template parameter).
  template<unsigned VAL>
  class CachingArrayInterpolator : public CachingArrayInterpolatorBase
  {

  public:
    /// Default constructor
    CachingArrayInterpolator() {}

    /// Destructor
    virtual ~CachingArrayInterpolator() {}

    virtual void build(const Vector<double>& s)
    {
      CachingArrayInterpolatorBase::build(s);

      Values.reset();
      Dvaluesdt.reset();
      for(unsigned i=0; i<VAL; i++)
        {
          Dvaluesdx[VAL].reset();
        }
    }

    const double* value()
    {
      if(!Values.is_initialised())
        {
          for(unsigned j=0; j<VAL; j++)
            {
              Values.set()[j] = Intp_pt->interpolate_value(j);
            }
        }

      return Values.get();
    }

    const double* dvaluedt()
    {
      if(!Dvaluesdt.is_initialised())
        {
          for(unsigned j=0; j<VAL; j++)
            {
              Dvaluesdt.set()[j] = Intp_pt->interpolate_dvaluedt(j);
            }
        }

      return Dvaluesdt.get();
    }

    const double* dvaluedx(const unsigned &i_val)
    {
      if(!Dvaluesdx[i_val].is_initialised())
        {
          for(unsigned j=0; j<Intp_pt->Dim; j++)
            {
              Dvaluesdx.set()[j] = Intp_pt->interpolate_dvaluedx(i_val, j);
            }
        }

      return Dvaluesdx[i_val].get();
    }

  protected:

    /// Storage arrays
    CachedArray<VAL> Values;
    CachedArray<VAL> Dvaluesdt;
    CachedArray<VAL> Dvaluesdx[3];
  };



} // End of oomph namespace

#endif
