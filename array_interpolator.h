#ifndef OOMPH_ARRAY_INTERPOLATOR_H
#define OOMPH_ARRAY_INTERPOLATOR_H

/*
  description of file goes here
*/


using namespace oomph;

namespace oomph
{

 // Forward decl.
 class FiniteElement;
 class Time;

 // class InterpolatorBase
 // {
 // public:
 //  virtual double time() = 0;

 //  virtual double x(const unsigned &i) = 0;
 //  virtual const Vector<double>& x() = 0;

 //  virtual double value(const unsigned &i) = 0;
 //  virtual double dvaluedt(const unsigned &i) = 0;
 //  virtual double dvaluedx(const unsigned &direction, const unsigned &i) = 0;
 //  virtual const Vector<double>& dvaluedx(const unsigned &i_val) = 0;

 //  //??ds shape + test + J
 //  virtual double J() = 0;
 // };


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
                           const Vector<double> &s)
   :
   // Cache some loop end conditions
   Nnode(this_element->nnode()),
   Dim(this_element->dim()),
   Nprev_value_zeroth_derivative(this_element->node_pt(0)->
                                 time_stepper_pt()->nprev_values_for_value_at_evaluation_time()),
   Nprev_value_zeroth_pos_derivative
   (this_element->node_pt(0)->position_time_stepper_pt()->nprev_values_for_value_at_evaluation_time()),

   //??ds weird!
   // Nprev_value_derivative(this_element->node_pt(0)->
   // time_stepper_pt()->nprev_values()),
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

  // Negative time to signify that it has not been calculated yet
   Intp_time(NotYetCalculatedValue)
  {
   // Set up shape + test functions
   J = this_element->dshape_eulerian(s, Psi, Dpsidx);
   Test = Psi;
   Dtestdx = Dpsidx;

   // Find out if any nodes in this element are hanging
   Has_hanging_nodes = false;
   for(unsigned nd=0, nnode=This_element->nnode(); nd<nnode; nd++)
    {
     Node* nd_pt = This_element->node_pt(nd);
     if(nd_pt->is_hanging())
      {
       Has_hanging_nodes = true;
       break;
      }
    }

#ifdef PARANOID
   for(unsigned nd=0, nnode=This_element->nnode(); nd<nnode; nd++)
    {
     Node* nd_pt = This_element->node_pt(nd);

     // Check all number of values the same for all nodes
     if(nd_pt->nvalue() != VAL)
      {
       std::ostringstream error_msg;
       error_msg << "Number of values must be the same at all nodes to use this interpolator.";
       throw OomphLibError(error_msg.str(),
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
      }

     // Check ts_pts are all the same
     if(nd_pt->time_stepper_pt() != This_element->node_pt(0)->time_stepper_pt())
      {
       std::ostringstream error_msg;
       error_msg << "Time steppers should all be the same within an element!";
       throw OomphLibError(error_msg.str(),
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
      }

    }
#endif

   // Initialise storage
   for(unsigned i=0; i<3; i++) X[i] = NotYetCalculatedValue;
   for(unsigned i=0; i<3; i++) Dxdt[i] = NotYetCalculatedValue;
   for(unsigned i=0; i<VAL; i++) Values[i] = NotYetCalculatedValue;
   for(unsigned i=0; i<VAL; i++) Dvaluesdt[i] = NotYetCalculatedValue;
   for(unsigned i=0; i<VAL; i++)
     {
       for(unsigned j=0; j<3; j++) Dvaluesdx[i][j] = NotYetCalculatedValue;
     }

  }

  /// Destructor
  ~GeneralArrayInterpolator() {}

  double time()
  {
   if(uninitialised(Intp_time)) {Intp_time = interpolate_time();}
   return Intp_time;
  }

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
   double time = 0.0;
   Time* time_pt = This_element->node_pt(0)->time_stepper_pt()->time_pt();
   for(unsigned i_tm=0; i_tm<Nprev_value_zeroth_derivative; i_tm++)
    {
     // the number of pointers chained together here seems a bit silly... ds
     time += time_pt->time(i_tm) * (*Ts_weights_pt)(0,i_tm);
    }
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
  const unsigned Nprev_value_zeroth_derivative;
  const unsigned Nprev_value_zeroth_pos_derivative;
  const unsigned Nprev_value_derivative;
  const unsigned Nprev_value_pos_derivative;

  bool Has_hanging_nodes;

  // Cached pointers
  const FiniteElement* This_element;
  const DenseMatrix<double>* Ts_weights_pt;
  const DenseMatrix<double>* Position_ts_weights_pt;

  // Jacobian + shape/test functions
  double J;
  Shape Psi;
  Shape Test;
  DShape Dpsidx;
  DShape Dtestdx;

  // Interpolated value storage (note we can't name the time variable
  // "Time" because that is used for the time class).
  double Intp_time;

  double X[3];
  double Dxdt[3];

  double Values[VAL];
  double Dvaluesdt[VAL];
  double Dvaluesdx[VAL][3];


  // Magic number signifying that the time has not been set yet
  static const double NotYetCalculatedValue = -987654.321;

  // Checks to see if private variables are initialised
  bool uninitialised(double var) const {return var == NotYetCalculatedValue;}

  /// \short Check to see if a vector is initialised. Recusive checking
  /// would be faster but takes a fairly large amount of time.
  template<typename T> bool uninitialised(Vector<T> var) const
  {
    bool uninitialised_entry = false;
    for(unsigned j=0; j<var.size(); j++)
     {
      if(uninitialised(var[j]))
       {
        uninitialised_entry = true;
        break;
       }
     }
    return var.empty() || uninitialised_entry;
    // return var.empty();
  }

   bool uninitialised(const double* var) const
     {
       if(var == 0) return true;
       else return uninitialised(var[0]);
     }

  // Position interpolation
  // ============================================================

   void interpolate_x_with_hanging_nodes()
   {
     for(unsigned j=0; j<Dim; j++)
       {
         X[j] = 0.0;
         for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
           {
             for(unsigned i_tm=0; i_tm<Nprev_value_zeroth_pos_derivative; i_tm++)
               {
                 X[j] += This_element->nodal_position(i_tm, i_nd, j)
                   * Psi(i_nd) * (*Position_ts_weights_pt)(0, i_tm);
               }
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
             for(unsigned i_tm=0; i_tm<Nprev_value_zeroth_pos_derivative; i_tm++)
               {
                 X[j] += This_element->raw_nodal_position(i_tm, i_nd, j)
                   * Psi(i_nd) * (*Position_ts_weights_pt)(0, i_tm);
               }
           }
       }
   }

  // Position derivatives interpolation
  // ============================================================

   void interpolate_dxdt_raw()
   {
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
             for(unsigned i_tm=0; i_tm<Nprev_value_zeroth_derivative; i_tm++)
               {
                 Values[j] += This_element->raw_nodal_value(i_tm, i_nd, j)
                   * Psi(i_nd) * (*Ts_weights_pt)(0, i_tm);
               }
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
             for(unsigned i_tm=0; i_tm<Nprev_value_zeroth_derivative; i_tm++)
               {
                 Values[j] += This_element->nodal_value(i_tm, i_nd, j)
                   * Psi(i_nd) * (*Ts_weights_pt)(0, i_tm);
               }
           }
       }
   }


  // Interpolate derivatives of values w.r.t. time
  // ============================================================

  void interpolate_dvaluesdt_raw(const unsigned &start,
                                           const unsigned &end)
  {
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
            for(unsigned i_tm=0; i_tm<Nprev_value_zeroth_derivative; i_tm++)
              {
                Dvaluesdx[i_value][i_direc] += This_element->raw_nodal_value(i_tm, i_nd, i_value)
                  * Dpsidx(i_nd, i_direc) * (*Ts_weights_pt)(0, i_tm);
              }
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
            for(unsigned i_tm=0; i_tm<Nprev_value_zeroth_derivative; i_tm++)
              {
                Dvaluesdx[i_value][i_direc] += This_element->nodal_value(i_tm, i_nd, i_value)
                  * Dpsidx(i_nd, i_direc) * (*Ts_weights_pt)(0, i_tm);
              }
          }
      }
  }




 };



} // End of oomph namespace

#endif
