#ifndef OOMPH_VARIABLE_QUADRATURE_H
#define OOMPH_VARIABLE_QUADRATURE_H

#include "generic.h"
#include <cmath>
/* Weights generated using C++ QUADRULE by John Burkardt.

   You may want to disable line wrapping for this file as the lines
   containing the data are very long. In emacs this is done using M-x
   toggle-truncate-lines.

   Ideally we would use a std::initializer_list to intialise the
   weights_data_structures but it needs c++0x.
*/

using namespace oomph;

namespace oomph
{


  //============================================================
  /// A structure to hold the weights/knots data. We need this so we can have a
  /// constructor to initialise the const static data structures.
  //============================================================
  struct weights_data_structure
  {

  private:

    /// Typedef to keep code lines reasonable and allow easy swapping between
    /// std::vector and oomph::Vector if needed. We need to use std:: because
    /// the oomph:: vector lacks some constructors
    typedef std::vector<double> vdb;

    /// The data structure itself
    std::map<unsigned,vdb > data;

    /// Given the data in arrays, construct the object
    void construction_helper(const double data_array[],
			     const unsigned order_array[],
			     const unsigned &order_array_length)
    {
      // Empty the data (just in case)
      data.clear();

      // Keep track of where we are up to in the (1D) input array
      unsigned k = 0;

      // Keep track of where we are up to in the map
      std::map<unsigned,vdb >::iterator it = data.begin();

      // for each order in the input order_array
      for(unsigned i=0; i<order_array_length; i++)
	{
	  // Get this order from order_array
	  unsigned order = order_array[i];

	  // Put 'order' values from data_array into a vector
	  // Not sure why but the compiler doesn't like [] notation here.
	  vdb temp(data_array +k, data_array + k + order);

	  // Insert data for this order into map just after last insertion
	  it = data.insert(it, std::pair<unsigned,vdb>(order, temp));

	  // Increase k ready for next order
	  k+=order;
	}
    }

    /// Get the data for GaussLegendre


  public:

    weights_data_structure(const unsigned scheme, const bool weight);

    // Return the i-th entry for the given order
    inline double operator()(const unsigned &i, const unsigned&order) const
    {return (*data.find(order)).second[i];}

    // Check whether a given order exists
    inline bool order_existence(const unsigned &order) const
    {return data.count(order) == 1;}
  };


  //==============================================================================
  /// An abstract base class for Gaussian quadrature. Does the error checking and
  /// declares all the virtual functions that will be needed.
  //==============================================================================
  class BaseVariableOrderQuadrature : public Integral
  {

  public:

    /// Default construtor
    BaseVariableOrderQuadrature(){}

    /// Broken copy constructor
    BaseVariableOrderQuadrature(const BaseVariableOrderQuadrature& dummy)
    {BrokenCopy::broken_copy("BaseVariableOrderQuadrature");}

    /// Broken assignment operator
    void operator=(const BaseVariableOrderQuadrature& dummy)
    {BrokenCopy::broken_assign("BaseVariableOrderQuadrature");}

    /// Destructor
    ~BaseVariableOrderQuadrature(){};

    /// Get the current Dim
    virtual unsigned dim() const = 0;

    /// Check that dim, order i and j are all within the appropriate ranges.
    void error_check(const unsigned &i, const unsigned &j,
		     const unsigned &order,
		     const std::string &function_name) const;

    /// Check that the requested order exists
    virtual bool order_existence(const unsigned &order) const = 0;

    /// Get the number of weights for given order (must be implemented in derived class).
    virtual unsigned nweight(const unsigned &order) const = 0;

    /// Get the weights for given order (must be implemented in derived class).
    virtual double weight(const unsigned &i, const unsigned &order) const = 0;

    /// Get the location of knots for given order (must be implemented in derived class).
    virtual double knot(const unsigned &i, const unsigned &j,
			const unsigned &order) const = 0;

    /// Get the weights for the one dimensional case (these are the only ones stored)
    virtual double weight_1d(const unsigned &i, const unsigned &order) const = 0;

    /// Get the knots for the one dimensional case (these are the only ones stored)
    virtual double knot_1d(const unsigned &i, const unsigned &order) const = 0;

    /// Dummy function to override the virtual one from Integral class
    virtual unsigned nweight() const
    {
      throw OomphLibError("Must specify an order for use with variable order integration",
			  "BaseVariableOrderQuadrature::nweight", OOMPH_EXCEPTION_LOCATION);
    }

    /// Dummy function to override the virtual one from Integral class
    virtual double weight(const unsigned &i) const
    {
      throw OomphLibError("Must specify an order for use with variable order integration",
			  "BaseVariableOrderQuadrature::weight", OOMPH_EXCEPTION_LOCATION);
    }

    /// Dummy function to override the virtual one from Integral class
    virtual double knot(const unsigned &i,const unsigned &j) const
    {
      throw OomphLibError("Must specify an order for use with variable order integration",
			  "BaseVariableOrderQuadrature::weight", OOMPH_EXCEPTION_LOCATION);
    }

  };

  //============================================================
  /// The geometry dependant parts of the variable order quadrature for
  /// QElements.
  //============================================================
  template<unsigned DIM>
  class QVariableOrderQuadrature : public BaseVariableOrderQuadrature
  {};


  //============================================================
  /// Specialisation of VariableOrderQuadrature to 1D.
  //============================================================
  template<>
  class QVariableOrderQuadrature<1> : public BaseVariableOrderQuadrature
  {
  public:

    /// Default construtor
    QVariableOrderQuadrature(){}

    /// Broken copy constructor
    QVariableOrderQuadrature(const QVariableOrderQuadrature& dummy)
    {BrokenCopy::broken_copy("QVariableOrderQuadrature");}

    /// Broken assignment operator
    void operator=(const QVariableOrderQuadrature& dummy)
    {BrokenCopy::broken_assign("QVariableOrderQuadrature");}


    inline double weight(const unsigned &i, const unsigned &order) const
    {
#ifdef PARANOID
      unsigned dummy = 0;
      error_check(i,dummy,order,"QVariableOrderQuadrature::weight");
#endif
      return weight_1d(i,order);
    }

    inline double knot(const unsigned &i, const unsigned &j, const unsigned &order) const
    {
#ifdef PARANOID
      error_check(i,j,order,"VariableOrderGaussLegendre::knot");
#endif
      return knot_1d(i,order);
    }

    inline unsigned dim() const {return 1;}

    inline unsigned nweight(const unsigned &order) const
    {return order;}
  };


  //============================================================
  /// Specialisation of QVariableOrderQuadrature to 2D.
  //============================================================
  template<>
  class QVariableOrderQuadrature<2> : public BaseVariableOrderQuadrature
  {
  public:

    /// Default construtor
    QVariableOrderQuadrature(){}

    /// Broken copy constructor
    QVariableOrderQuadrature(const QVariableOrderQuadrature& dummy)
    {BrokenCopy::broken_copy("QVariableOrderQuadrature");}

    /// Broken assignment operator
    void operator=(const QVariableOrderQuadrature& dummy)
    {BrokenCopy::broken_assign("QVariableOrderQuadrature");}

    inline double weight(const unsigned &i, const unsigned &order) const
    {
#ifdef PARANOID
      unsigned dummy = 0;
      error_check(i,dummy,order,"QVariableOrderQuadrature::weight");
#endif
      unsigned i_x = i%order;
      unsigned i_y = i/order;
      return weight_1d(i_y,order)*weight_1d(i_x,order);
    }

    inline double knot(const unsigned &i, const unsigned &j, const unsigned &order) const
    {
#ifdef PARANOID
      error_check(i,j,order,"VariableOrderGaussLegendre::knot");
#endif
      if(j==0)
	{
	  unsigned i_x = i%order;
	  return knot_1d(i_x,order);
	}
      else if(j==1)
	{
	  unsigned i_y = i/order;
	  return knot_1d(i_y,order);
	}
      else
	{
	  std::ostringstream error_stream;
	  error_stream << "Requested knot coordinate in dimension " << j
		       << " which does not exist in a scheme of dimension 2" << std::endl;
	  throw OomphLibError(error_stream.str(),
			      "QVariableOrderQuadrature<2>::knot",
			      OOMPH_EXCEPTION_LOCATION);
	}
    }

    inline unsigned dim() const {return 2;}

    inline unsigned nweight(const unsigned &order) const
    {return order*order;}
  };


  //============================================================
  /// Specialisation of QVariableOrderQuadrature to 3D.
  //============================================================
  template<>
  class QVariableOrderQuadrature<3> : public BaseVariableOrderQuadrature
  {
  public:

    /// Default construtor
    QVariableOrderQuadrature(){}

    /// Broken copy constructor
    QVariableOrderQuadrature(const QVariableOrderQuadrature& dummy)
    {BrokenCopy::broken_copy("QVariableOrderQuadrature");}

    /// Broken assignment operator
    void operator=(const QVariableOrderQuadrature& dummy)
    {BrokenCopy::broken_assign("QVariableOrderQuadrature");}

    inline double weight(const unsigned &i, const unsigned &order) const
    {
#ifdef PARANOID
      unsigned dummy = 0;
      error_check(i,dummy,order,"QVariableOrderQuadrature::weight");
#endif
      unsigned i_x = i%order;
      unsigned i_y = (i/order)%order;
      unsigned i_z = i/(order*order);
      return weight_1d(i_z,order)*weight_1d(i_y,order)*weight_1d(i_x,order);
    }
    inline double knot(const unsigned &i, const unsigned &j, const unsigned &order) const
    {
#ifdef PARANOID
      error_check(i,j,order,"VariableOrderGaussLegendre::knot");
#endif
      if(j==0)
	{
	  unsigned i_x = i%order;
	  return knot_1d(i_x,order);
	}
      else if(j==1)
	{
	  unsigned i_y = (i/order)%order;
	  return knot_1d(i_y,order);
	}
      else if(j==2)
	{
	  unsigned i_z = i/(order*order);
	  return knot_1d(i_z,order);
	}
      else
	{
      std::ostringstream error_stream;
      error_stream << "Requested knot coordinate in dimension " << j
		   << " which does not exist in a scheme of dimension 3" << std::endl;
      throw OomphLibError(error_stream.str(),
			       "QVariableOrderQuadrature<3>::knot",
			       OOMPH_EXCEPTION_LOCATION);
	}
    }

    inline unsigned dim() const {return 3;}

    inline unsigned nweight(const unsigned &order) const
    {return order*order*order;}
  };

  //============================================================
  /// The geometry dependant parts of the variable order quadrature for
  /// TElements.
  //============================================================
  template<unsigned DIM>
  class TVariableOrderQuadrature : public BaseVariableOrderQuadrature
  {};

  //============================================================
  /// Specialisation of VariableOrderQuadrature to 1D triangle/tet
  /// elements. One-dimensional triangles are the same as 1D quads but scaled so
  /// that their local coordinate runs from 0 to 1, rather than -1 to 1.
  // ============================================================
  template<>
  class TVariableOrderQuadrature<1> : public BaseVariableOrderQuadrature
  {
  public:

    /// Default construtor
    TVariableOrderQuadrature(){}

    /// Broken copy constructor
    TVariableOrderQuadrature(const TVariableOrderQuadrature& dummy)
    {BrokenCopy::broken_copy("TVariableOrderQuadrature");}

    /// Broken assignment operator
    void operator=(const TVariableOrderQuadrature& dummy)
    {BrokenCopy::broken_assign("TVariableOrderQuadrature");}


    inline double weight(const unsigned &i, const unsigned &order) const
    {
#ifdef PARANOID
      unsigned dummy = 0;
      error_check(i,dummy,order,"TVariableOrderQuadrature::weight");
#endif
      return weight_1d(i,order);
    }

    inline double knot(const unsigned &i, const unsigned &j, const unsigned &order) const
    {
#ifdef PARANOID
      error_check(i,j,order,"VariableOrderGaussLegendre::knot");
#endif
      return 0.5*(knot_1d(i,order) + 1);
    }

    inline unsigned dim() const {return 1;}

    inline unsigned nweight(const unsigned &order) const
    {return order;}
  };

  //============================================================
  /// Specialisation of TVariableOrderQuadrature to 2D.
  //??ds not implemented yet
  //============================================================
  template<>
  class TVariableOrderQuadrature<2> : public BaseVariableOrderQuadrature
  {
  public:

    /// Default construtor
    TVariableOrderQuadrature(){}

    /// Broken copy constructor
    TVariableOrderQuadrature(const TVariableOrderQuadrature& dummy)
    {BrokenCopy::broken_copy("TVariableOrderQuadrature");}

    /// Broken assignment operator
    void operator=(const TVariableOrderQuadrature& dummy)
    {BrokenCopy::broken_assign("TVariableOrderQuadrature");}

    inline double weight(const unsigned &i, const unsigned &order) const
    {
#ifdef PARANOID
      unsigned dummy = 0;
      error_check(i,dummy,order,"TVariableOrderQuadrature::weight");
#endif
      // Get weight for quadrilateral
      //??ds could replace this by a call to Q quadrature?
      unsigned i_x = i%order;
      unsigned i_y = i/order;
      double qweight = weight_1d(i_y,order)*weight_1d(i_x,order);

      // Get the x-location of this knot
      double knot_x =  this->knot(i,0,order);

      // Transform to triangle
      return qweight * (1.0 - knot_x)/8.0;
    }

    inline double knot(const unsigned &i, const unsigned &j, const unsigned &order) const
    {
#ifdef PARANOID
      error_check(i,j,order,"VariableOrderGaussLegendre::knot");
#endif

      // Get knot for quadrilateral and transform
      if(j==0)
	{
	  unsigned i_x = i%order;
	  return (1.0 + knot_1d(i_x,order))/2.0;
	}
      else if(j==1)
	{
	  unsigned i_x = i%order;
	  unsigned i_y = i/order;

	  return (1.0 - knot_1d(i_x,order))
	    * (1.0 + knot_1d(i_y,order))
	    / 4.0;
	}
      else
	{
	  std::ostringstream error_stream;
	  error_stream << "Requested knot coordinate in dimension " << j
		       << " which does not exist in a scheme of dimension 2" << std::endl;
	  throw OomphLibError(error_stream.str(),
			      "QVariableOrderQuadrature<2>::knot",
			      OOMPH_EXCEPTION_LOCATION);
	}

    }

    inline unsigned dim() const {return 2;}

    inline unsigned nweight(const unsigned &order) const
    {return order*order;}
  };


  //============================================================
  /// Specialisation of TVariableOrderQuadrature to 3D.
  //??ds not implemented yet
  //============================================================
  template<>
  class TVariableOrderQuadrature<3> : public BaseVariableOrderQuadrature
  {
  public:

    /// Default construtor
    TVariableOrderQuadrature(){}

    /// Broken copy constructor
    TVariableOrderQuadrature(const TVariableOrderQuadrature& dummy)
    {BrokenCopy::broken_copy("TVariableOrderQuadrature");}

    /// Broken assignment operator
    void operator=(const TVariableOrderQuadrature& dummy)
    {BrokenCopy::broken_assign("TVariableOrderQuadrature");}

    inline double weight(const unsigned &i, const unsigned &order) const
    {
#ifdef PARANOID
      unsigned dummy = 0;
      error_check(i,dummy,order,"TVariableOrderQuadrature::weight");
#endif
      return 0.0;
    }

    inline double knot(const unsigned &i, const unsigned &j, const unsigned &order) const
    {
#ifdef PARANOID
      error_check(i,j,order,"VariableOrderQuadrature::knot");
#endif
      std::ostringstream error_stream;
      error_stream << "Requested knot coordinate for 3D triangles, not yet implemented for variable order." << std::endl;
      throw OomphLibError(error_stream.str(),
			       "TVariableOrderQuadrature<3>::knot",
			       OOMPH_EXCEPTION_LOCATION);
    }

    inline unsigned dim() const {return 3;}

    inline unsigned nweight(const unsigned &order) const
    {return order*order*order;}
  };

  //============================================================
  /// Gauss-Legendre quadrature (the standard Gaussian quadrature used
  /// in oomph-lib).
  //============================================================
  class VariableOrderGaussLegendre
  {
  private:

    /// Weights for all orders and knots
    static const weights_data_structure Weights;

    /// Positions for all orders and knots
    static const weights_data_structure Knots;

  public:

    /// Return weight of the i-th integration point in 1D.
    inline double weight_1d(const unsigned &i, const unsigned &order) const
    {return Weights(i,order);}

    /// Return local coordinate s of i-th integration point in 1D.
    inline double knot_1d(const unsigned &i, const unsigned &order) const
    {return Knots(i,order);}

    /// Check that this order exists.
    inline bool order_existence(const unsigned &order) const
    {return Weights.order_existence(order);}

    /// Get the first order to use in an adaptive scheme
    inline unsigned adaptive_start_order() const
    {return 2;}

    /// Get the next order to use in an adaptive scheme.
    inline unsigned adaptive_next_order(const unsigned &order) const
    {
      return 2*order;
    }

  };


  //============================================================
  /// Clenshaw-Curtis quadrature
  /// Advantage: higher order methods re-use the same knots.
  //============================================================
  class VariableOrderClenshawCurtis
  {
  private:

    /// Weights for all orders and knots
    static const weights_data_structure Weights;

    /// Positions for all orders and knots
    static const weights_data_structure Knots;

  public:

    /// Return weight of the i-th integration point in 1D.
    inline double weight_1d(const unsigned &i, const unsigned &order) const
    {return Weights(i,order);}

    /// Return local coordinate s of i-th integration point in 1D.
    inline double knot_1d(const unsigned &i, const unsigned &order) const
    {return Knots(i,order);}

    /// Check that this order exists.
    inline bool order_existence(const unsigned &order) const
    {return Knots.order_existence(order);}

    /// \short Get the index of matching knots in higher order schemes. Only
    /// applicable for progressive quadratures (Clenshaw-Curtis, Fejer's
    /// second rule, Patterson) for some orders. Useful in adaptive schemes
    /// when we want to reuse higher order calculations.
    // This works for orders where n_high = 2^a * n_low, where n =
    // (order - 1) and a is some int.  i.e. order_high = 2^a *
    // (order_low - 1) + 1 An example is 4 -> 7 -> 13 -> 25 -> 49 (so
    // max_order = 49 works for all of these orders). Also order 2 works
    // with anything since the only knots are the two endpoints.
    inline unsigned find_corresponding_knot(const unsigned &i,
					    const unsigned &order,
					    const unsigned &high_order) const
    {
      // If order is two then the knots are always the two endpoints
      if(order == 2)
	{
	  if(i==0) return 0;
	  else return (high_order-1);
	}

      // The number of knots from the higher order scheme to skip when using the
      // lower order scheme.
      // (based on rearangement of (high_order - 1) = 2^a * (low_order -1) )
      double a = log2(high_order - 1) - log2(order - 1);

#ifdef PARANOID
      // If a is non-integer then we can't do this
      double dummy;
      if(modf(a,&dummy) > 1e-14)
  	{
  	  std::ostringstream error_stream;
  	  error_stream << "The high order scheme must be such that"
  		       << "n_high = 2^a * n_low, where n = (order - 1)"
  		       << "and a is an integer. Here a = "
  		       << a << ", which is non-integer.";
  	  throw OomphLibError(error_stream.str(),
  			      "VariableOrderClenshawCurtis::find_corresponding_knot",
  			      OOMPH_EXCEPTION_LOCATION);
  	}
#endif
      return i * unsigned(pow(2,unsigned(a)));
    }

    /// Get the first order to use in an adaptive scheme
    inline unsigned adaptive_start_order() const
    {return 2;}

    /// Get the next highest order allowing reuse of all previous knots
    inline unsigned adaptive_next_order(const unsigned &order) const
    {
      if(order == 2) return 4;
      else return (2*order) - 1;
    }
  };


  //============================================================
  /// Fejer's second rule quadrature
  /// Advantage: higher order methods re-use the same knots.
  /// Exactly the same as Clenshaw Curtis except that the endpoints
  /// are not included so it can be used for integrals with endpoint
  /// singularities.
  //============================================================
  class VariableOrderFejerSecond
  {

    /// Weights for all orders and knots
    static const weights_data_structure Weights;

    /// Locations for all orders and all knots
    static const weights_data_structure Knots;

  public:

    /// Constructor with reasonable defaults for min/max adaptive order
    VariableOrderFejerSecond(){}

    /// Return weight of the i-th integration point in 1D.
    inline double weight_1d(const unsigned &i, const unsigned &order) const
    {return Weights(i,order);}

    /// Return local coordinate s of i-th integration point in 1D.
    inline double knot_1d(const unsigned &i, const unsigned &order) const
    {return Knots(i,order);}

    /// Check that this order exists.
    inline bool order_existence(const unsigned &order) const
    {return Weights.order_existence(order);}

    /// \short Get the index of matching knots in higher order
    /// schemes. Only applicable for progressive quadratures
    /// (Clenshaw-Curtis, Fejer's second rule, Patterson) for some
    /// orders. Useful in adaptive schemes when we want to reuse
    /// higher order calculations.
    // Let k = #points = order, x_i = i*pi/(k+1).
    // nested scheme when k_high = 2^a *(k_low + 1) - 1.
    // (Alternatively x_i = i*pi/n, n = #points + 1
    // nested scheme when n_high = 2^a *n+low)
    inline unsigned find_corresponding_knot(const unsigned &i,
  					    const unsigned &order,
  					    const unsigned &high_order) const
    {
      // The number of knots from the higher order scheme to skip when using the
      // lower order scheme.
      // (based on rearangement of (high_order - 1) = 2^a * (low_order -1) )
      double a = log2(high_order + 1) - log2(order + 1);

#ifdef PARANOID
      // If a is non-integer then we can't do this
      double dummy;
      if(modf(a,&dummy) > 1e-14)
  	{
  	  std::ostringstream error_stream;
  	  error_stream << "The high order scheme must be such that"
  		       << "n_high = 2^a * n_low, where n = (order - 1)"
  		       << "and a is an integer. Here a = "
  		       << a << ", which is non-integer.";
  	  throw OomphLibError(error_stream.str(),
  			      "VariableOrderFejerSecond::find_corresponding_knot",
  			      OOMPH_EXCEPTION_LOCATION);
  	}
#endif
      return i * unsigned(pow(2,unsigned(a)));
    }

    /// Get the first order to use in an adaptive scheme
    inline unsigned adaptive_start_order() const
    {return 2;}

    /// Get the next highest order allowing reuse of previous knots
    inline unsigned adaptive_next_order(const unsigned &order) const
    {
      return (2*order) + 1;
    }
  };


  //============================================================
  /// The final class for GaussLegendre on a QElement
  //============================================================
  template <unsigned DIM>
  class QVariableOrderGaussLegendre : public VariableOrderGaussLegendre,
				      public QVariableOrderQuadrature<DIM>
  {
  public:
    // Just make sure we are calling the right functions
    //?? I think this shouldn't be necessary but it seems to be...
    double weight_1d(const unsigned &i, const unsigned &order) const
    {return VariableOrderGaussLegendre::weight_1d(i,order);}

    double knot_1d(const unsigned &i, const unsigned &order) const
    {return VariableOrderGaussLegendre::knot_1d(i,order);}

    bool order_existence(const unsigned &order) const
    {return VariableOrderGaussLegendre::order_existence(order);}
  };


  //============================================================
  /// The final class for ClenshawCurtis on a QElement
  //============================================================
  template <unsigned DIM>
  class QVariableOrderClenshawCurtis : public VariableOrderClenshawCurtis,
				       public QVariableOrderQuadrature<DIM>
  {
  public:
    double weight_1d(const unsigned &i, const unsigned &order) const
    {return VariableOrderClenshawCurtis::weight_1d(i,order);}

    double knot_1d(const unsigned &i, const unsigned &order) const
    {return VariableOrderClenshawCurtis::knot_1d(i,order);}

    bool order_existence(const unsigned &order) const
    {return VariableOrderClenshawCurtis::order_existence(order);}
  };


  //============================================================
  /// The final class for FejerSecond on a QElement
  //============================================================
  template <unsigned DIM>
  class QVariableOrderFejerSecond : public VariableOrderFejerSecond,
				    public QVariableOrderQuadrature<DIM>
  {
  public:
    double weight_1d(const unsigned &i, const unsigned &order) const
    {return VariableOrderFejerSecond::weight_1d(i,order);}

    double knot_1d(const unsigned &i, const unsigned &order) const
    {return VariableOrderFejerSecond::knot_1d(i,order);}

    bool order_existence(const unsigned &order) const
    {return VariableOrderFejerSecond::order_existence(order);}
  };


  //============================================================
  /// The final class for GaussLegendre on a TElement
  //============================================================
  template <unsigned DIM>
  class TVariableOrderGaussLegendre : public VariableOrderGaussLegendre,
				      public TVariableOrderQuadrature<DIM>
  {
  public:
    // Just make sure we are calling the right functions
    //?? I think this shouldn't be necessary but it seems to be...
    double weight_1d(const unsigned &i, const unsigned &order) const
    {return VariableOrderGaussLegendre::weight_1d(i,order);}

    double knot_1d(const unsigned &i, const unsigned &order) const
    {return VariableOrderGaussLegendre::knot_1d(i,order);}

    bool order_existence(const unsigned &order) const
    {return VariableOrderGaussLegendre::order_existence(order);}
  };


  //============================================================
  /// The final class for ClenshawCurtis on a TElement
  //============================================================
  template <unsigned DIM>
  class TVariableOrderClenshawCurtis : public VariableOrderClenshawCurtis,
				       public TVariableOrderQuadrature<DIM>
  {
  public:
    double weight_1d(const unsigned &i, const unsigned &order) const
    {return VariableOrderClenshawCurtis::weight_1d(i,order);}

    double knot_1d(const unsigned &i, const unsigned &order) const
    {return VariableOrderClenshawCurtis::knot_1d(i,order);}

    bool order_existence(const unsigned &order) const
    {return VariableOrderClenshawCurtis::order_existence(order);}
  };


  //============================================================
  /// The final class for FejerSecond on a TElement
  //============================================================
  template <unsigned DIM>
  class TVariableOrderFejerSecond : public VariableOrderFejerSecond,
				    public TVariableOrderQuadrature<DIM>
  {
  public:
    double weight_1d(const unsigned &i, const unsigned &order) const
    {return VariableOrderFejerSecond::weight_1d(i,order);}

    double knot_1d(const unsigned &i, const unsigned &order) const
    {return VariableOrderFejerSecond::knot_1d(i,order);}

    bool order_existence(const unsigned &order) const
    {return VariableOrderFejerSecond::order_existence(order);}
  };


  // this should go into the .cc file eventually
  //////////////////////////////////////////////////////////////////////////////

  //============================================================
  /// General range checking function for use in all weight and knot functions.
  //============================================================
  void BaseVariableOrderQuadrature::
  error_check(const unsigned &i, const unsigned &j,
  	      const unsigned &order, const std::string &function_name) const
  {
    if((dim() > 3) || (dim() < 1))
      {
  	std::ostringstream error_stream;
  	error_stream << "Invalid dimension " << dim();
  	throw OomphLibError(error_stream.str(),function_name,
  			    OOMPH_EXCEPTION_LOCATION);
      }

    if(!(order_existence(order)))
      {
  	std::ostringstream error_stream;
  	error_stream << "Gaussian quadrature of order "
  		     << order << " is not (yet) implemented.";
  	throw OomphLibError(error_stream.str(),function_name,
  			    OOMPH_EXCEPTION_LOCATION);
      }

    if(i >= nweight(order))
      {
  	std::ostringstream error_stream;
  	error_stream << dim() << " dimensional Gaussian quadrature of order "
  		     << order << " does not have a weight number " << i
  		     << ".";
  	throw OomphLibError(error_stream.str(),function_name,
  			    OOMPH_EXCEPTION_LOCATION);
      }

    if(j >= dim())
      {
  	std::ostringstream error_stream;
  	error_stream << "Dimension of quadrature is " << dim()
  		     << " so coordinate " << j << " does not exist";
  	throw OomphLibError(error_stream.str(),function_name,
  			    OOMPH_EXCEPTION_LOCATION);
      }
  }
  // Automatically generated by generate_quadrature_rules_driver, based on QUADRULE.
  //See https://github.com/davidshepherd7/oomph-lib-additions/tree/master/generate_quadrature_rules .
  weights_data_structure::weights_data_structure(const unsigned scheme, const bool weights)
  {
    switch(scheme)
      {
      case 0 :
	{

	  double weights_data_array[] =
	    {
	      2,
	      0.9999999999999998,0.9999999999999998,
	      0.5555555555555556,0.8888888888888888,0.5555555555555556,
	      0.3478548451374538,0.6521451548625461,0.6521451548625461,0.3478548451374538,
	      0.2369268850561891,0.4786286704993664,0.5688888888888889,0.4786286704993664,0.2369268850561891,
	      0.1713244923791703,0.3607615730481386,0.467913934572691,0.467913934572691,0.3607615730481386,0.1713244923791703,
	      0.1294849661688698,0.2797053914892766,0.3818300505051191,0.4179591836734694,0.3818300505051191,0.2797053914892766,0.1294849661688698,
	      0.1012285362903761,0.2223810344533745,0.3137066458778873,0.362683783378362,0.362683783378362,0.3137066458778873,0.2223810344533745,0.1012285362903761,
	      0.08127438836157448,0.1806481606948574,0.2606106964029355,0.3123470770400028,0.3302393550012598,0.3123470770400028,0.2606106964029355,0.1806481606948574,0.08127438836157448,
	      0.06667134430868811,0.1494513491505806,0.2190863625159819,0.2692667193099964,0.2955242247147528,0.2955242247147528,0.2692667193099964,0.2190863625159819,0.1494513491505806,0.06667134430868811,
	      0.05566856711617371,0.1255803694649046,0.1862902109277342,0.2331937645919905,0.2628045445102467,0.2729250867779006,0.2628045445102467,0.2331937645919905,0.1862902109277342,0.1255803694649046,0.05566856711617371,
	      0.04717533638651188,0.1069393259953184,0.1600783285433462,0.203167426723066,0.2334925365383548,0.2491470458134029,0.2491470458134029,0.2334925365383548,0.203167426723066,0.1600783285433462,0.1069393259953184,0.04717533638651188,
	      0.04048400476531593,0.09212149983772849,0.1388735102197872,0.1781459807619457,0.2078160475368886,0.2262831802628972,0.2325515532308739,0.2262831802628972,0.2078160475368886,0.1781459807619457,0.1388735102197872,0.09212149983772849,0.04048400476531593,
	      0.03511946033175192,0.08015808715976019,0.1215185706879032,0.1572031671581935,0.1855383974779378,0.2051984637212955,0.2152638534631577,0.2152638534631577,0.2051984637212955,0.1855383974779378,0.1572031671581935,0.1215185706879032,0.08015808715976019,0.03511946033175192,
	      0.03075324199611711,0.07036604748810804,0.107159220467172,0.1395706779261543,0.1662692058169939,0.1861610000155622,0.1984314853271116,0.2025782419255613,0.1984314853271116,0.1861610000155622,0.1662692058169939,0.1395706779261543,0.107159220467172,0.07036604748810804,0.03075324199611711,
	      0.02715245941175407,0.06225352393864787,0.09515851168249277,0.1246289712555339,0.1495959888165767,0.1691565193950026,0.1826034150449236,0.1894506104550686,0.1894506104550686,0.1826034150449236,0.1691565193950026,0.1495959888165767,0.1246289712555339,0.09515851168249277,0.06225352393864787,0.02715245941175407,
	      0.02414830286854788,0.05545952937398714,0.08503614831717912,0.111883847193404,0.1351363684685255,0.1540457610768103,0.16800410215645,0.1765627053669926,0.1794464703562065,0.1765627053669926,0.16800410215645,0.1540457610768103,0.1351363684685255,0.111883847193404,0.08503614831717912,0.05545952937398714,0.02414830286854788,
	      0.02161601352648345,0.04971454889496989,0.07642573025488909,0.1009420441062871,0.1225552067114785,0.1406429146706507,0.1546846751262654,0.1642764837458327,0.1691423829631436,0.1691423829631436,0.1642764837458327,0.1546846751262654,0.1406429146706507,0.1225552067114785,0.1009420441062871,0.07642573025488909,0.04971454889496989,0.02161601352648345,
	      0.01946178822972657,0.04481422676569961,0.06904454273764124,0.09149002162245001,0.111566645547334,0.1287539625393363,0.1426067021736067,0.1527660420658597,0.1589688433939543,0.1610544498487837,0.1589688433939543,0.1527660420658597,0.1426067021736067,0.1287539625393363,0.111566645547334,0.09149002162245001,0.06904454273764124,0.04481422676569961,0.01946178822972657,
	      0.01761400713915219,0.0406014298003869,0.06267204833410905,0.08327674157670473,0.1019301198172404,0.1181945319615184,0.1316886384491766,0.142096109318382,0.1491729864726038,0.1527533871307259,0.1527533871307259,0.1491729864726038,0.142096109318382,0.1316886384491766,0.1181945319615184,0.1019301198172404,0.08327674157670473,0.06267204833410905,0.0406014298003869,0.01761400713915219,
	      0.01601722825777442,0.03695378977085248,0.05713442542685721,0.07610011362837933,0.09344442345603385,0.1087972991671483,0.1218314160537286,0.1322689386333375,0.1398873947910731,0.14452440398997,0.1460811336496904,0.14452440398997,0.1398873947910731,0.1322689386333375,0.1218314160537286,0.1087972991671483,0.09344442345603385,0.07610011362837933,0.05713442542685721,0.03695378977085248,0.01601722825777442,
	      0.01462799529827204,0.03377490158481423,0.05229333515268329,0.06979646842452042,0.08594160621706774,0.100414144442881,0.1129322960805392,0.1232523768105124,0.1311735047870624,0.1365414983460151,0.139251872855632,0.139251872855632,0.1365414983460151,0.1311735047870624,0.1232523768105124,0.1129322960805392,0.100414144442881,0.08594160621706774,0.06979646842452042,0.05229333515268329,0.03377490158481423,0.01462799529827204,
	      0.01341185948714169,0.03098800585697949,0.04803767173108468,0.06423242140852584,0.07928141177671892,0.09291576606003518,0.1048920914645413,0.1149966402224113,0.1230490843067296,0.1289057221880822,0.1324620394046966,0.1336545721861062,0.1324620394046966,0.1289057221880822,0.1230490843067296,0.1149966402224113,0.1048920914645413,0.09291576606003518,0.07928141177671892,0.06423242140852584,0.04803767173108468,0.03098800585697949,0.01341185948714169,
	      0.01234122979998734,0.02853138862893367,0.04427743881741976,0.05929858491543673,0.07334648141108031,0.08619016153195329,0.09761865210411388,0.1074442701159656,0.1155056680537257,0.1216704729278034,0.1258374563468283,0.1279381953467521,0.1279381953467521,0.1258374563468283,0.1216704729278034,0.1155056680537257,0.1074442701159656,0.09761865210411388,0.08619016153195329,0.07334648141108031,0.05929858491543673,0.04427743881741976,0.02853138862893367,0.01234122979998734,
	      0.01139379850102619,0.0263549866150321,0.04093915670130634,0.05490469597583521,0.06803833381235691,0.08014070033500099,0.09102826198296368,0.1005359490670506,0.1085196244742636,0.1148582591457117,0.1194557635357847,0.1222424429903101,0.1231760537267154,0.1222424429903101,0.1194557635357847,0.1148582591457117,0.1085196244742636,0.1005359490670506,0.09102826198296368,0.08014070033500099,0.06803833381235691,0.05490469597583521,0.04093915670130634,0.0263549866150321,0.01139379850102619,
	      0.01055137261734303,0.02441785109263185,0.03796238329436279,0.05097582529714778,0.0632740463295748,0.07468414976565975,0.08504589431348523,0.09421380035591412,0.1020591610944254,0.1084718405285765,0.1133618165463197,0.1166604434852966,0.1183214152792623,0.1183214152792623,0.1166604434852966,0.1133618165463197,0.1084718405285765,0.1020591610944254,0.09421380035591412,0.08504589431348523,0.07468414976565975,0.0632740463295748,0.05097582529714778,0.03796238329436279,0.02441785109263185,0.01055137261734303,
	      0.009798996051294261,0.02268623159618062,0.03529705375741978,0.04744941252061508,0.05898353685983358,0.0697488237662456,0.07960486777305782,0.08842315854375692,0.09608872737002848,0.1025016378177458,0.1075782857885332,0.1112524883568451,0.1134763461089652,0.114220867378957,0.1134763461089652,0.1112524883568451,0.1075782857885332,0.1025016378177458,0.09608872737002848,0.08842315854375692,0.07960486777305782,0.0697488237662456,0.05898353685983358,0.04744941252061508,0.03529705375741978,0.02268623159618062,0.009798996051294261,
	      0.009124282593094564,0.02113211259277116,0.03290142778230441,0.04427293475900424,0.05510734567571671,0.06527292396699957,0.07464621423456876,0.08311341722890124,0.09057174439303282,0.09693065799792988,0.1021129675780608,0.1060557659228464,0.1087111922582942,0.1100470130164752,0.1100470130164752,0.1087111922582942,0.1060557659228464,0.1021129675780608,0.09693065799792988,0.09057174439303282,0.08311341722890124,0.07464621423456876,0.06527292396699957,0.05510734567571671,0.04427293475900424,0.03290142778230441,0.02113211259277116,0.009124282593094564,
	      0.008516903878746448,0.01973208505612275,0.03074049220209361,0.04140206251868283,0.05159482690249791,0.06120309065707913,0.07011793325505127,0.07823832713576377,0.0854722573661725,0.09173775713925873,0.09696383409440855,0.1010912737599149,0.1040733100777294,0.105876155097321,0.1064793817183142,0.105876155097321,0.1040733100777294,0.1010912737599149,0.09696383409440855,0.09173775713925873,0.0854722573661725,0.07823832713576377,0.07011793325505127,0.06120309065707913,0.05159482690249791,0.04140206251868283,0.03074049220209361,0.01973208505612275,0.008516903878746448,
	      0.007968192496166718,0.01846646831109101,0.02878470788332334,0.03879919256962703,0.04840267283059402,0.0574931562176191,0.06597422988218048,0.07375597473770522,0.08075589522942024,0.08689978720108299,0.09212252223778619,0.0963687371746442,0.09959342058679523,0.1017623897484055,0.1028526528935589,0.1028526528935589,0.1017623897484055,0.09959342058679523,0.0963687371746442,0.09212252223778619,0.08689978720108299,0.08075589522942024,0.07375597473770522,0.06597422988218048,0.0574931562176191,0.04840267283059402,0.03879919256962703,0.02878470788332334,0.01846646831109101,0.007968192496166718,
	      0.007470831579248847,0.01731862079031063,0.0270090191849794,0.03643227391238548,0.0454937075272011,0.05410308242491686,0.06217478656102845,0.06962858323541035,0.0763903865987766,0.08239299176158921,0.08757674060847793,0.0918901138936415,0.09529024291231951,0.09774333538632875,0.09922501122667229,0.09972054479342646,0.09922501122667229,0.09774333538632875,0.09529024291231951,0.0918901138936415,0.08757674060847793,0.08239299176158921,0.0763903865987766,0.06962858323541035,0.06217478656102845,0.05410308242491686,0.0454937075272011,0.03643227391238548,0.0270090191849794,0.01731862079031063,0.007470831579248847,
	      0.007018610009470076,0.01627439473090562,0.02539206530926209,0.03427386291302146,0.04283589802222672,0.0509980592623762,0.05868409347853559,0.06582222277636185,0.07234579410884846,0.07819389578707037,0.08331192422694679,0.08765209300440382,0.09117387869576393,0.09384439908080458,0.09563872007927485,0.09654008851472773,0.09654008851472773,0.09563872007927485,0.09384439908080458,0.09117387869576393,0.08765209300440382,0.08331192422694679,0.07819389578707037,0.07234579410884846,0.06582222277636185,0.05868409347853559,0.0509980592623762,0.04283589802222672,0.03427386291302146,0.02539206530926209,0.01627439473090562,0.007018610009470076,
	      0.00660622784758745,0.01532170151293468,0.02391554810174949,0.03230035863232895,0.04040154133166957,0.04814774281871172,0.05547084663166357,0.06230648253031747,0.0685945728186567,0.07427985484395411,0.07931236479488672,0.08364787606703868,0.08724828761884436,0.0900819586606386,0.09212398664331686,0.09335642606559616,0.09376844616021,0.09335642606559616,0.09212398664331686,0.0900819586606386,0.08724828761884436,0.08364787606703868,0.07931236479488672,0.07427985484395411,0.0685945728186567,0.06230648253031747,0.05547084663166357,0.04814774281871172,0.04040154133166957,0.03230035863232895,0.02391554810174949,0.01532170151293468,0.00660622784758745,
	      0.006229140555908603,0.01445016274859499,0.02256372198549499,0.03049138063844612,0.03816659379638752,0.04552561152335328,0.05250741457267813,0.05905413582752447,0.06511152155407644,0.07062937581425575,0.07556197466003194,0.07986844433977189,0.08351309969984569,0.08646573974703579,0.08870189783569386,0.09020304437064074,0.09095674033025994,0.09095674033025994,0.09020304437064074,0.08870189783569386,0.08646573974703579,0.08351309969984569,0.07986844433977189,0.07556197466003194,0.07062937581425575,0.06511152155407644,0.05905413582752447,0.05250741457267813,0.04552561152335328,0.03816659379638752,0.03049138063844612,0.02256372198549499,0.01445016274859499,0.006229140555908603,
	      0.005883433420443039,0.01365082834836153,0.02132297991148358,0.02882926010889425,0.03611011586346342,0.04310842232617023,0.04976937040135352,0.05604081621237009,0.0618736719660802,0.0672222852690869,0.07204479477256009,0.07630345715544208,0.07996494224232431,0.08300059372885658,0.08538665339209912,0.08710444699718353,0.08814053043027543,0.08848679490710429,0.08814053043027543,0.08710444699718353,0.08538665339209912,0.08300059372885658,0.07996494224232431,0.07630345715544208,0.07204479477256009,0.0672222852690869,0.0618736719660802,0.05604081621237009,0.04976937040135352,0.04310842232617023,0.03611011586346342,0.02882926010889425,0.02132297991148358,0.01365082834836153,0.005883433420443039,
	      0.005565719664245137,0.01291594728406564,0.02018151529773544,0.02729862149856875,0.03421381077030725,0.04087575092364494,0.04723508349026596,0.0532447139777599,0.05886014424532479,0.0640397973550155,0.06874532383573642,0.07294188500565302,0.07659841064587064,0.0796878289120716,0.08218726670433978,0.0840782189796619,0.08534668573933865,0.08598327567039478,0.08598327567039478,0.08534668573933865,0.0840782189796619,0.08218726670433978,0.0796878289120716,0.07659841064587064,0.07294188500565302,0.06874532383573642,0.0640397973550155,0.05886014424532479,0.0532447139777599,0.04723508349026596,0.04087575092364494,0.03421381077030725,0.02729862149856875,0.02018151529773544,0.01291594728406564,0.005565719664245137,
	      0.005273057279498011,0.01223878010030757,0.01912904448908395,0.02588603699055894,0.03246163984752147,0.03880960250193453,0.04488536466243719,0.05064629765482463,0.0560519879982749,0.06106451652322597,0.06564872287275121,0.06977245155570036,0.07340677724848819,0.07652620757052923,0.07910886183752942,0.08113662450846498,0.0825952722364372,0.08347457362586275,0.0837683609931389,0.08347457362586275,0.0825952722364372,0.08113662450846498,0.07910886183752942,0.07652620757052923,0.07340677724848819,0.06977245155570036,0.06564872287275121,0.06106451652322597,0.0560519879982749,0.05064629765482463,0.04488536466243719,0.03880960250193453,0.03246163984752147,0.02588603699055894,0.01912904448908395,0.01223878010030757,0.005273057279498011,
	      0.005002880749639449,0.01161344471646874,0.01815657770961324,0.02457973973823235,0.03083950054517504,0.03689408159402473,0.04270315850467446,0.04822806186075867,0.05343201991033231,0.05828039914699722,0.06274093339213305,0.06678393797914041,0.07038250706689897,0.07351269258474348,0.07615366354844635,0.07828784465821093,0.07990103324352783,0.08098249377059716,0.0815250292803858,0.0815250292803858,0.08098249377059716,0.07990103324352783,0.07828784465821093,0.07615366354844635,0.07351269258474348,0.07038250706689897,0.06678393797914041,0.06274093339213305,0.05828039914699722,0.05343201991033231,0.04822806186075867,0.04270315850467446,0.03689408159402473,0.03083950054517504,0.02457973973823235,0.01815657770961324,0.01161344471646874,0.005002880749639449,
	      0.004752944691635074,0.01103478893916452,0.01725622909372493,0.02336938483217819,0.0293349559839034,0.03511511149813134,0.04067327684793382,0.04597430110891662,0.05098466529212938,0.05567269034091627,0.06000873608859617,0.06396538813868244,0.06751763096623122,0.07064300597060871,0.07332175341426865,0.07553693732283598,0.07727455254468202,0.07852361328737119,0.07927622256836848,0.07952762213944285,0.07927622256836848,0.07852361328737119,0.07727455254468202,0.07553693732283598,0.07332175341426865,0.07064300597060871,0.06751763096623122,0.06396538813868244,0.06000873608859617,0.05567269034091627,0.05098466529212938,0.04597430110891662,0.04067327684793382,0.03511511149813134,0.0293349559839034,0.02336938483217819,0.01725622909372493,0.01103478893916452,0.004752944691635074,
	      0.004521277098533074,0.01049828453115285,0.01642105838190785,0.02224584919416692,0.02793700698002342,0.03346019528254784,0.03878216797447201,0.04387090818567327,0.04869580763507219,0.05322784698393685,0.05743976909939155,0.06130624249292892,0.06480401345660104,0.06791204581523394,0.07061164739128678,0.07288658239580408,0.07472316905796823,0.07611036190062621,0.07703981816424797,0.07750594797842481,0.07750594797842481,0.07703981816424797,0.07611036190062621,0.07472316905796823,0.07288658239580408,0.07061164739128678,0.06791204581523394,0.06480401345660104,0.06130624249292892,0.05743976909939155,0.05322784698393685,0.04869580763507219,0.04387090818567327,0.03878216797447201,0.03346019528254784,0.02793700698002342,0.02224584919416692,0.01642105838190785,0.01049828453115285,0.004521277098533074,
	      0.004306140358164757,0.009999938773905906,0.01564493840781855,0.02120106336877954,0.02663589920711042,0.03191821173169929,0.03701771670350797,0.04190519519590966,0.04655264836901434,0.05093345429461751,0.05502251924257878,0.05879642094987191,0.06223354258096631,0.06531419645352744,0.06802073676087679,0.0703376606208175,0.07225169686102313,0.07375188202722349,0.07482962317622155,0.07547874709271586,0.07569553564729838,0.07547874709271586,0.07482962317622155,0.07375188202722349,0.07225169686102313,0.0703376606208175,0.06802073676087679,0.06531419645352744,0.06223354258096631,0.05879642094987191,0.05502251924257878,0.05093345429461751,0.04655264836901434,0.04190519519590966,0.03701771670350797,0.03191821173169929,0.02663589920711042,0.02120106336877954,0.01564493840781855,0.009999938773905906,0.004306140358164757,
	      0.004105998604648935,0.009536220301748532,0.01492244369735747,0.0202278695690526,0.02542295952611303,0.03047924069960345,0.03536907109759215,0.04006573518069225,0.04454357777196589,0.04877814079280324,0.05274629569917407,0.05642636935801843,0.05979826222758666,0.06284355804500258,0.06554562436490899,0.06788970337652191,0.06986299249259417,0.07145471426517098,0.07265617524380413,0.07346081345346753,0.07386423423217291,0.07386423423217291,0.07346081345346753,0.07265617524380413,0.07145471426517098,0.06986299249259417,0.06788970337652191,0.06554562436490899,0.06284355804500258,0.05979826222758666,0.05642636935801843,0.05274629569917407,0.04877814079280324,0.04454357777196589,0.04006573518069225,0.03536907109759215,0.03047924069960345,0.02542295952611303,0.0202278695690526,0.01492244369735747,0.009536220301748532,0.004105998604648935,
	      0.003919490253844022,0.009103996637401362,0.01424875643157652,0.01931990142368389,0.02429045661383886,0.02913441326149852,0.0338264920868603,0.03834222219413268,0.04265805719798212,0.04675149475434658,0.0506011927843901,0.05418708031888175,0.05749046195691045,0.06049411524999128,0.06318238044939613,0.06554124212632274,0.06755840222936516,0.06922334419365669,0.07052738776508501,0.07146373425251416,0.07202750197142195,0.07221575169379899,0.07202750197142195,0.07146373425251416,0.07052738776508501,0.06922334419365669,0.06755840222936516,0.06554124212632274,0.06318238044939613,0.06049411524999128,0.05749046195691045,0.05418708031888175,0.0506011927843901,0.04675149475434658,0.04265805719798212,0.03834222219413268,0.0338264920868603,0.02913441326149852,0.02429045661383886,0.01931990142368389,0.01424875643157652,0.009103996637401362,0.003919490253844022,
	      0.003745404803112892,0.008700481367524762,0.01361958675557998,0.01847148173681475,0.02323148190201921,0.02787578282128102,0.0323812228120698,0.03672534781380885,0.0408865123103462,0.04484398408197006,0.04857804644835201,0.05207009609170443,0.05530273556372806,0.05825985987759545,0.060926736701562,0.06329007973320386,0.0653381148791814,0.06706063890629371,0.06844907026936667,0.06949649186157261,0.07019768547355819,0.07054915778935406,0.07054915778935406,0.07019768547355819,0.06949649186157261,0.06844907026936667,0.06706063890629371,0.0653381148791814,0.06329007973320386,0.060926736701562,0.05825985987759545,0.05530273556372806,0.05207009609170443,0.04857804644835201,0.04484398408197006,0.0408865123103462,0.03672534781380885,0.0323812228120698,0.02787578282128102,0.02323148190201921,0.01847148173681475,0.01361958675557998,0.008700481367524762,0.003745404803112892,
	      0.003582663155283468,0.008323189296218291,0.01303110499158283,0.01767753525793756,0.02223984755057873,0.02669621396757765,0.03102537493451546,0.03520669220160903,0.03922023672930242,0.04304688070916497,0.04666838771837333,0.05006749923795206,0.05322801673126895,0.0561348787597865,0.05877423271884172,0.06113350083106651,0.06320144007381996,0.06496819575072343,0.06642534844984251,0.06756595416360758,0.06838457737866967,0.06887731697766133,0.06904182482923202,0.06887731697766133,0.06838457737866967,0.06756595416360758,0.06642534844984251,0.06496819575072343,0.06320144007381996,0.06113350083106651,0.05877423271884172,0.0561348787597865,0.05322801673126895,0.05006749923795206,0.04666838771837333,0.04304688070916497,0.03922023672930242,0.03520669220160903,0.03102537493451546,0.02669621396757765,0.02223984755057873,0.01767753525793756,0.01303110499158283,0.008323189296218291,0.003582663155283468,
	      0.003430300868107123,0.007969898229724593,0.01247988377098866,0.01693351400783624,0.02130999875413648,0.02558928639713,0.02975182955220275,0.0337786279991069,0.03765130535738613,0.04135219010967871,0.04486439527731811,0.04817189510171219,0.05125959800714303,0.05411341538585677,0.05672032584399116,0.05906843459554632,0.06114702772465049,0.06294662106439447,0.06445900346713912,0.06567727426778121,0.06659587476845485,0.06721061360067816,0.06751868584903647,0.06751868584903647,0.06721061360067816,0.06659587476845485,0.06567727426778121,0.06445900346713912,0.06294662106439447,0.06114702772465049,0.05906843459554632,0.05672032584399116,0.05411341538585677,0.05125959800714303,0.04817189510171219,0.04486439527731811,0.04135219010967871,0.03765130535738613,0.0337786279991069,0.02975182955220275,0.02558928639713,0.02130999875413648,0.01693351400783624,0.01247988377098866,0.007969898229724593,0.003430300868107123,
	      0.003287453842528111,0.007638616295848896,0.01196284846431233,0.01623533314643302,0.02043693814766842,0.02454921165965883,0.02855415070064337,0.03243423551518473,0.03617249658417498,0.03975258612253103,0.04315884864847955,0.0463763890865059,0.0493911377473612,0.05218991178005719,0.05476047278153028,0.05709158029323159,0.05917304094233885,0.06099575300873965,0.06255174622092172,0.06383421660571703,0.06483755623894574,0.06555737776654977,0.06599053358881048,0.06613512962365548,0.06599053358881048,0.06555737776654977,0.06483755623894574,0.06383421660571703,0.06255174622092172,0.06099575300873965,0.05917304094233885,0.05709158029323159,0.05476047278153028,0.05218991178005719,0.0493911377473612,0.0463763890865059,0.04315884864847955,0.03975258612253103,0.03617249658417498,0.03243423551518473,0.02855415070064337,0.02454921165965883,0.02043693814766842,0.01623533314643302,0.01196284846431233,0.007638616295848896,0.003287453842528111,
	      0.003153346052305979,0.007327553901276245,0.01147723457923454,0.01557931572294384,0.01961616045735553,0.02357076083932438,0.02742650970835695,0.03116722783279812,0.03477722256477047,0.03824135106583069,0.04154508294346471,0.04467456085669427,0.04761665849249046,0.05035903555385442,0.05289018948519365,0.05519950369998413,0.05727729210040319,0.05911483969839568,0.0607044391658939,0.0620394231598926,0.06311419228625405,0.0639242385846482,0.06446616443595007,0.06473769681268392,0.06473769681268392,0.06446616443595007,0.0639242385846482,0.06311419228625405,0.0620394231598926,0.0607044391658939,0.05911483969839568,0.05727729210040319,0.05519950369998413,0.05289018948519365,0.05035903555385442,0.04761665849249046,0.04467456085669427,0.04154508294346471,0.03824135106583069,0.03477722256477047,0.03116722783279812,0.02742650970835695,0.02357076083932438,0.01961616045735553,0.01557931572294384,0.01147723457923454,0.007327553901276245,0.003153346052305979,
	      0.003027278988922775,0.007035099590086438,0.01102055103159355,0.01496214493562463,0.01884359585308948,0.02264920158744667,0.026363618927066,0.02997188462058378,0.03345946679162218,0.03681232096300069,0.04001694576637301,0.04306043698125963,0.04593053935559581,0.04861569588782823,0.05110509433014458,0.05338871070825901,0.05545734967480361,0.05730268153018748,0.05891727576002728,0.060294630953152,0.06142920097919295,0.06231641732005729,0.06295270746519567,0.06333550929649175,0.0634632814047906,0.06333550929649175,0.06295270746519567,0.06231641732005729,0.06142920097919295,0.060294630953152,0.05891727576002728,0.05730268153018748,0.05545734967480361,0.05338871070825901,0.05110509433014458,0.04861569588782823,0.04593053935559581,0.04306043698125963,0.04001694576637301,0.03681232096300069,0.03345946679162218,0.02997188462058378,0.026363618927066,0.02264920158744667,0.01884359585308948,0.01496214493562463,0.01102055103159355,0.007035099590086438,0.003027278988922775,
	      0.002908622553155216,0.006759799195745393,0.01059054838365099,0.01438082276148561,0.01811556071348942,0.02178024317012481,0.0253606735700124,0.02884299358053521,0.03221372822357803,0.03545983561514612,0.03856875661258769,0.04152846309014768,0.04432750433880329,0.04695505130394842,0.0494009384494663,0.05165570306958116,0.05371062188899624,0.05555774480621253,0.05718992564772838,0.05860084981322243,0.05978505870426552,0.06073797084177025,0.06145589959031669,0.06193606742068326,0.06217661665534725,0.06217661665534725,0.06193606742068326,0.06145589959031669,0.06073797084177025,0.05978505870426552,0.05860084981322243,0.05718992564772838,0.05555774480621253,0.05371062188899624,0.05165570306958116,0.0494009384494663,0.04695505130394842,0.04432750433880329,0.04152846309014768,0.03856875661258769,0.03545983561514612,0.03221372822357803,0.02884299358053521,0.0253606735700124,0.02178024317012481,0.01811556071348942,0.01438082276148561,0.01059054838365099,0.006759799195745393,0.002908622553155216,
	      0.001783280721696354,0.004147033260562425,0.006504457968978381,0.00884675982636393,0.01116813946013115,0.01346304789671863,0.01572603047602473,0.01795171577569735,0.02013482315353023,0.02227017380838325,0.02435270256871088,0.02637746971505469,0.02833967261425946,0.03023465707240246,0.03205792835485158,0.0338051618371416,0.03547221325688237,0.03705512854024007,0.03855015317861563,0.03995374113272031,0.04126256324262354,0.04247351512365355,0.04358372452932345,0.04459055816375655,0.04549162792741814,0.04628479658131443,0.04696818281621005,0.0475401657148303,0.04799938859645831,0.048344762234803,0.0485754674415034,0.04869095700913974,0.04869095700913974,0.0485754674415034,0.048344762234803,0.04799938859645831,0.0475401657148303,0.04696818281621005,0.04628479658131443,0.04549162792741814,0.04459055816375655,0.04358372452932345,0.04247351512365355,0.04126256324262354,0.03995374113272031,0.03855015317861563,0.03705512854024007,0.03547221325688237,0.0338051618371416,0.03205792835485158,0.03023465707240246,0.02833967261425946,0.02637746971505469,0.02435270256871088,0.02227017380838325,0.02013482315353023,0.01795171577569735,0.01572603047602473,0.01346304789671863,0.01116813946013115,0.00884675982636393,0.006504457968978381,0.004147033260562425,0.001783280721696354,
	      0.0004493809602922042,0.001045812679340394,0.001642503018669007,0.002238288430962618,0.00283275147145803,0.003425526040910214,0.004016254983738653,0.004604584256702968,0.005190161832676347,0.005772637542865698,0.006351663161707211,0.006926892566898817,0.007497981925634733,0.008064589890486064,0.008626377798616741,0.00918300987166088,0.009734153415006807,0.01027947901583215,0.01081866073950307,0.01135137632408042,0.01187730737274031,0.01239613954395092,0.01290756273926734,0.01341127128861632,0.01390696413295197,0.01439434500416684,0.01487312260214733,0.01534301076886515,0.01580372865939935,0.01625500090978518,0.0166965578015892,0.01712813542311136,0.01754947582711771,0.0179603271850087,0.01836044393733133,0.0187495869405447,0.01912752360995097,0.01949402805870662,0.01984888123283087,0.02019187104213003,0.02052279248696008,0.02084144778075115,0.02114764646822136,0.02144120553920845,0.02172194953805207,0.02198971066846048,0.02224432889379981,0.022485652032745,0.02271353585023647,0.02292784414368688,0.02312844882438702,0.02331522999406275,0.02348807601653592,0.02364688358444763,0.0237915577810034,0.02392201213670344,0.02403816868102403,0.02413995798901929,0.02422731922281526,0.02430020016797182,0.02435855726469062,0.02440235563384958,0.02443156909785007,0.0244461801962625,0.0244461801962625,0.02443156909785007,0.02440235563384958,0.02435855726469062,0.02430020016797182,0.02422731922281526,0.02413995798901929,0.02403816868102403,0.02392201213670344,0.0237915577810034,0.02364688358444763,0.02348807601653592,0.02331522999406275,0.02312844882438702,0.02292784414368688,0.02271353585023647,0.022485652032745,0.02224432889379981,0.02198971066846048,0.02172194953805207,0.02144120553920845,0.02114764646822136,0.02084144778075115,0.02052279248696008,0.02019187104213003,0.01984888123283087,0.01949402805870662,0.01912752360995097,0.0187495869405447,0.01836044393733133,0.0179603271850087,0.01754947582711771,0.01712813542311136,0.0166965578015892,0.01625500090978518,0.01580372865939935,0.01534301076886515,0.01487312260214733,0.01439434500416684,0.01390696413295197,0.01341127128861632,0.01290756273926734,0.01239613954395092,0.01187730737274031,0.01135137632408042,0.01081866073950307,0.01027947901583215,0.009734153415006807,0.00918300987166088,0.008626377798616741,0.008064589890486064,0.007497981925634733,0.006926892566898817,0.006351663161707211,0.005772637542865698,0.005190161832676347,0.004604584256702968,0.004016254983738653,0.003425526040910214,0.00283275147145803,0.002238288430962618,0.001642503018669007,0.001045812679340394,0.0004493809602922042,
	      0.000112789017822257,0.0002625349442964987,0.0004124632544261519,0.0005623489540314241,0.0007121541634733123,0.0008618537014200986,0.001011424393208445,0.001160843557567717,0.001310088681902503,0.00145913733331074,0.001607967130749334,0.001756555736330738,0.001904880853499722,0.00205292022796613,0.002200651649839918,0.002348052956327318,0.002495102034703705,0.002641776825427483,0.002788055325327704,0.002933915590829717,0.003079335741199324,0.003224293961794194,0.003368768507315554,0.003512737705056316,0.003656179958142505,0.003799073748766254,0.003941397641408836,0.004083130286052643,0.004224250421381536,0.004364736877968053,0.004504568581447895,0.00464372455568007,0.004782183925892697,0.004919925921813863,0.005056929880786843,0.005193175250869278,0.005328641593915931,0.005463308588644312,0.005597156033682922,0.005730163850601431,0.00586231208692264,0.005993580919115338,0.006123950655567932,0.006253401739542403,0.006381914752107884,0.006509470415053646,0.006636049593781065,0.00676163330017379,0.006886202695446341,0.007009739092969824,0.007132223961075389,0.00725363892583392,0.00737396577381235,0.007493186454805879,0.007611283084545664,0.007728237947381562,0.007844033498939714,0.007958652368754359,0.008072077362873497,0.008184291466438273,0.008295277846235226,0.008405019853221545,0.008513501025022498,0.00862070508840102,0.008726615961698816,0.008831217757248745,0.008934494783758212,0.00903643154866288,0.009137012760450801,0.009236223330956304,0.009334048377623272,0.009430473225737762,0.009525483410629292,0.009619064679840727,0.00971120299526629,0.00980188453525732,0.00989109569669582,0.009978823097034906,0.01006505357630638,0.01014977419909487,0.01023297225647822,0.01031463526793402,0.01039475098321172,0.0104733073841704,0.01055029268658148,0.01062569534189656,0.01069950403897979,0.01077170770580463,0.0108422955111148,0.01091125686604904,0.01097858142572955,0.0110442590908139,0.01110828000900984,0.01117063457655344,0.01123131343964968,0.01129030749587552,0.01134760789554551,0.01140320604303919,0.01145709359809062,0.01150926247703949,0.01155970485404364,0.01160841316225307,0.01165538009494523,0.01170059860662073,0.01174406191406056,0.01178576349734342,0.01182569710082398,0.01186385673407109,0.01190023667276649,0.01193483145956359,0.01196763590490591,0.01199864508780581,0.01202785435658256,0.01205525932956014,0.01208085589572453,0.01210464021534048,0.01212660872052733,0.01214675811579446,0.01216508537853549,0.01218158775948178,0.01219626278311472,0.01220910824803723,0.01222012222730399,0.01222930306871028,0.01223664939504017,0.0122421601042728,0.0122458343697479,0.01224767164028975,0.01224767164028975,0.0122458343697479,0.0122421601042728,0.01223664939504017,0.01222930306871028,0.01222012222730399,0.01220910824803723,0.01219626278311472,0.01218158775948178,0.01216508537853549,0.01214675811579446,0.01212660872052733,0.01210464021534048,0.01208085589572453,0.01205525932956014,0.01202785435658256,0.01199864508780581,0.01196763590490591,0.01193483145956359,0.01190023667276649,0.01186385673407109,0.01182569710082398,0.01178576349734342,0.01174406191406056,0.01170059860662073,0.01165538009494523,0.01160841316225307,0.01155970485404364,0.01150926247703949,0.01145709359809062,0.01140320604303919,0.01134760789554551,0.01129030749587552,0.01123131343964968,0.01117063457655344,0.01110828000900984,0.0110442590908139,0.01097858142572955,0.01091125686604904,0.0108422955111148,0.01077170770580463,0.01069950403897979,0.01062569534189656,0.01055029268658148,0.0104733073841704,0.01039475098321172,0.01031463526793402,0.01023297225647822,0.01014977419909487,0.01006505357630638,0.009978823097034906,0.00989109569669582,0.00980188453525732,0.00971120299526629,0.009619064679840727,0.009525483410629292,0.009430473225737762,0.009334048377623272,0.009236223330956304,0.009137012760450801,0.00903643154866288,0.008934494783758212,0.008831217757248745,0.008726615961698816,0.00862070508840102,0.008513501025022498,0.008405019853221545,0.008295277846235226,0.008184291466438273,0.008072077362873497,0.007958652368754359,0.007844033498939714,0.007728237947381562,0.007611283084545664,0.007493186454805879,0.00737396577381235,0.00725363892583392,0.007132223961075389,0.007009739092969824,0.006886202695446341,0.00676163330017379,0.006636049593781065,0.006509470415053646,0.006381914752107884,0.006253401739542403,0.006123950655567932,0.005993580919115338,0.00586231208692264,0.005730163850601431,0.005597156033682922,0.005463308588644312,0.005328641593915931,0.005193175250869278,0.005056929880786843,0.004919925921813863,0.004782183925892697,0.00464372455568007,0.004504568581447895,0.004364736877968053,0.004224250421381536,0.004083130286052643,0.003941397641408836,0.003799073748766254,0.003656179958142505,0.003512737705056316,0.003368768507315554,0.003224293961794194,0.003079335741199324,0.002933915590829717,0.002788055325327704,0.002641776825427483,0.002495102034703705,0.002348052956327318,0.002200651649839918,0.00205292022796613,0.001904880853499722,0.001756555736330738,0.001607967130749334,0.00145913733331074,0.001310088681902503,0.001160843557567717,0.001011424393208445,0.0008618537014200986,0.0007121541634733123,0.0005623489540314241,0.0004124632544261519,0.0002625349442964987,0.000112789017822257,
	    };

	  double knots_data_array[] =
	    {
	      0,
	      -0.5773502691896257,0.5773502691896257,
	      -0.7745966692414834,0,0.7745966692414834,
	      -0.8611363115940526,-0.3399810435848563,0.3399810435848563,0.8611363115940526,
	      -0.906179845938664,-0.5384693101056831,0,0.5384693101056831,0.906179845938664,
	      -0.9324695142031521,-0.6612093864662645,-0.2386191860831969,0.2386191860831969,0.6612093864662645,0.9324695142031521,
	      -0.9491079123427585,-0.7415311855993945,-0.4058451513773972,0,0.4058451513773972,0.7415311855993945,0.9491079123427585,
	      -0.9602898564975363,-0.7966664774136267,-0.525532409916329,-0.1834346424956498,0.1834346424956498,0.525532409916329,0.7966664774136267,0.9602898564975363,
	      -0.9681602395076261,-0.8360311073266358,-0.6133714327005904,-0.3242534234038089,0,0.3242534234038089,0.6133714327005904,0.8360311073266358,0.9681602395076261,
	      -0.9739065285171717,-0.8650633666889845,-0.6794095682990244,-0.4333953941292472,-0.1488743389816312,0.1488743389816312,0.4333953941292472,0.6794095682990244,0.8650633666889845,0.9739065285171717,
	      -0.978228658146057,-0.8870625997680953,-0.7301520055740494,-0.5190961292068118,-0.269543155952345,0,0.269543155952345,0.5190961292068118,0.7301520055740494,0.8870625997680953,0.978228658146057,
	      -0.9815606342467192,-0.9041172563704749,-0.7699026741943047,-0.5873179542866175,-0.3678314989981802,-0.1252334085114689,0.1252334085114689,0.3678314989981802,0.5873179542866175,0.7699026741943047,0.9041172563704749,0.9815606342467192,
	      -0.9841830547185881,-0.9175983992229779,-0.8015780907333099,-0.6423493394403402,-0.4484927510364469,-0.2304583159551348,0,0.2304583159551348,0.4484927510364469,0.6423493394403402,0.8015780907333099,0.9175983992229779,0.9841830547185881,
	      -0.9862838086968123,-0.9284348836635735,-0.827201315069765,-0.6872929048116855,-0.5152486363581541,-0.3191123689278897,-0.1080549487073437,0.1080549487073437,0.3191123689278897,0.5152486363581541,0.6872929048116855,0.827201315069765,0.9284348836635735,0.9862838086968123,
	      -0.9879925180204855,-0.937273392400706,-0.8482065834104272,-0.7244177313601701,-0.5709721726085388,-0.3941513470775634,-0.2011940939974345,0,0.2011940939974345,0.3941513470775634,0.5709721726085388,0.7244177313601701,0.8482065834104272,0.937273392400706,0.9879925180204855,
	      -0.9894009349916499,-0.9445750230732326,-0.8656312023878318,-0.755404408355003,-0.6178762444026438,-0.4580167776572274,-0.2816035507792589,-0.09501250983763744,0.09501250983763744,0.2816035507792589,0.4580167776572274,0.6178762444026438,0.755404408355003,0.8656312023878318,0.9445750230732326,0.9894009349916499,
	      -0.9905754753144174,-0.9506755217687678,-0.8802391537269859,-0.7815140038968014,-0.6576711592166907,-0.5126905370864769,-0.3512317634538763,-0.1784841814958479,0,0.1784841814958479,0.3512317634538763,0.5126905370864769,0.6576711592166907,0.7815140038968014,0.8802391537269859,0.9506755217687678,0.9905754753144174,
	      -0.9915651684209309,-0.9558239495713977,-0.8926024664975557,-0.8037049589725231,-0.6916870430603532,-0.5597708310739475,-0.4117511614628426,-0.2518862256915055,-0.08477501304173529,0.08477501304173529,0.2518862256915055,0.4117511614628426,0.5597708310739475,0.6916870430603532,0.8037049589725231,0.8926024664975557,0.9558239495713977,0.9915651684209309,
	      -0.9924068438435844,-0.96020815213483,-0.9031559036148179,-0.8227146565371428,-0.7209661773352294,-0.600545304661681,-0.4645707413759609,-0.3165640999636298,-0.1603586456402254,0,0.1603586456402254,0.3165640999636298,0.4645707413759609,0.600545304661681,0.7209661773352294,0.8227146565371428,0.9031559036148179,0.96020815213483,0.9924068438435844,
	      -0.9931285991850949,-0.9639719272779138,-0.9122344282513259,-0.8391169718222188,-0.7463319064601508,-0.636053680726515,-0.5108670019508271,-0.3737060887154195,-0.2277858511416451,-0.07652652113349732,0.07652652113349732,0.2277858511416451,0.3737060887154195,0.5108670019508271,0.636053680726515,0.7463319064601508,0.8391169718222188,0.9122344282513259,0.9639719272779138,0.9931285991850949,
	      -0.9937521706203895,-0.9672268385663063,-0.9200993341504008,-0.8533633645833173,-0.7684399634756779,-0.6671388041974123,-0.5516188358872198,-0.4243421202074388,-0.2880213168024011,-0.1455618541608951,0,0.1455618541608951,0.2880213168024011,0.4243421202074388,0.5516188358872198,0.6671388041974123,0.7684399634756779,0.8533633645833173,0.9200993341504008,0.9672268385663063,0.9937521706203895,
	      -0.9942945854823994,-0.9700604978354287,-0.926956772187174,-0.8658125777203002,-0.7878168059792081,-0.6944872631866827,-0.5876404035069116,-0.469355837986757,-0.3419358208920842,-0.2078604266882213,-0.06973927331972223,0.06973927331972223,0.2078604266882213,0.3419358208920842,0.469355837986757,0.5876404035069116,0.6944872631866827,0.7878168059792081,0.8658125777203002,0.926956772187174,0.9700604978354287,0.9942945854823994,
	      -0.9947693349975522,-0.9725424712181152,-0.932971086826016,-0.8767523582704416,-0.8048884016188399,-0.7186613631319502,-0.6196098757636461,-0.5095014778460075,-0.3903010380302908,-0.264135680970345,-0.1332568242984661,0,0.1332568242984661,0.264135680970345,0.3903010380302908,0.5095014778460075,0.6196098757636461,0.7186613631319502,0.8048884016188399,0.8767523582704416,0.932971086826016,0.9725424712181152,0.9947693349975522,
	      -0.9951872199970213,-0.9747285559713095,-0.9382745520027328,-0.8864155270044011,-0.8200019859739029,-0.7401241915785544,-0.6480936519369755,-0.5454214713888396,-0.4337935076260451,-0.3150426796961634,-0.1911188674736163,-0.06405689286260563,0.06405689286260563,0.1911188674736163,0.3150426796961634,0.4337935076260451,0.5454214713888396,0.6480936519369755,0.7401241915785544,0.8200019859739029,0.8864155270044011,0.9382745520027328,0.9747285559713095,0.9951872199970213,
	      -0.9955569697904981,-0.9766639214595175,-0.9429745712289743,-0.8949919978782753,-0.833442628760834,-0.7592592630373577,-0.6735663684734684,-0.5776629302412229,-0.473002731445715,-0.3611723058093879,-0.2438668837209884,-0.1228646926107104,0,0.1228646926107104,0.2438668837209884,0.3611723058093879,0.473002731445715,0.5776629302412229,0.6735663684734684,0.7592592630373577,0.833442628760834,0.8949919978782753,0.9429745712289743,0.9766639214595175,0.9955569697904981,
	      -0.9958857011456169,-0.978385445956471,-0.9471590666617142,-0.9026378619843071,-0.8454459427884981,-0.7763859488206789,-0.6964272604199573,-0.6066922930176181,-0.5084407148245057,-0.4030517551234863,-0.2920048394859569,-0.1768588203568902,-0.05923009342931321,0.05923009342931321,0.1768588203568902,0.2920048394859569,0.4030517551234863,0.5084407148245057,0.6066922930176181,0.6964272604199573,0.7763859488206789,0.8454459427884981,0.9026378619843071,0.9471590666617142,0.978385445956471,0.9958857011456169,
	      -0.9961792628889886,-0.9799234759615012,-0.9509005578147049,-0.9094823206774911,-0.8562079080182945,-0.7917716390705082,-0.7170134737394237,-0.6329079719464952,-0.5405515645794569,-0.4411482517500269,-0.3359939036385089,-0.2264593654395368,-0.11397258560953,0,0.11397258560953,0.2264593654395368,0.3359939036385089,0.4411482517500269,0.5405515645794569,0.6329079719464952,0.7170134737394237,0.7917716390705082,0.8562079080182945,0.9094823206774911,0.9509005578147049,0.9799234759615012,0.9961792628889886,
	      -0.9964424975739544,-0.9813031653708728,-0.9542592806289382,-0.9156330263921321,-0.8658925225743951,-0.8056413709171791,-0.7356108780136318,-0.656651094038865,-0.5697204718114017,-0.4758742249551183,-0.3762515160890787,-0.2720616276351781,-0.1645692821333808,-0.05507928988403427,0.05507928988403427,0.1645692821333808,0.2720616276351781,0.3762515160890787,0.4758742249551183,0.5697204718114017,0.656651094038865,0.7356108780136318,0.8056413709171791,0.8658925225743951,0.9156330263921321,0.9542592806289382,0.9813031653708728,0.9964424975739544,
	      -0.9966794422605966,-0.9825455052614132,-0.9572855957780877,-0.9211802329530587,-0.8746378049201028,-0.8181854876152524,-0.7524628517344771,-0.6782145376026865,-0.5962817971382278,-0.5075929551242276,-0.4131528881740086,-0.3140316378676399,-0.2113522861660011,-0.1062782301326792,0,0.1062782301326792,0.2113522861660011,0.3140316378676399,0.4131528881740086,0.5075929551242276,0.5962817971382278,0.6782145376026865,0.7524628517344771,0.8181854876152524,0.8746378049201028,0.9211802329530587,0.9572855957780877,0.9825455052614132,0.9966794422605966,
	      -0.9968934840746495,-0.9836681232797472,-0.9600218649683075,-0.9262000474292743,-0.8825605357920527,-0.8295657623827684,-0.7677774321048262,-0.6978504947933158,-0.6205261829892429,-0.5366241481420199,-0.4470337695380892,-0.3527047255308781,-0.2546369261678899,-0.1538699136085835,-0.05147184255531769,0.05147184255531769,0.1538699136085835,0.2546369261678899,0.3527047255308781,0.4470337695380892,0.5366241481420199,0.6205261829892429,0.6978504947933158,0.7677774321048262,0.8295657623827684,0.8825605357920527,0.9262000474292743,0.9600218649683075,0.9836681232797472,0.9968934840746495,
	      -0.997087481819477,-0.9846859096651525,-0.9625039250929497,-0.9307569978966481,-0.8897600299482711,-0.8399203201462674,-0.7817331484166249,-0.7157767845868532,-0.6427067229242603,-0.5632491614071493,-0.4781937820449025,-0.3883859016082329,-0.2947180699817016,-0.1981211993355706,-0.09955531215234152,0,0.09955531215234152,0.1981211993355706,0.2947180699817016,0.3883859016082329,0.4781937820449025,0.5632491614071493,0.6427067229242603,0.7157767845868532,0.7817331484166249,0.8399203201462674,0.8897600299482711,0.9307569978966481,0.9625039250929497,0.9846859096651525,0.997087481819477,
	      -0.9972638618494816,-0.9856115115452684,-0.9647622555875064,-0.9349060759377397,-0.8963211557660521,-0.84936761373257,-0.7944837959679424,-0.7321821187402897,-0.6630442669302152,-0.5877157572407623,-0.5068999089322294,-0.4213512761306353,-0.3318686022821277,-0.2392873622521371,-0.1444719615827965,-0.04830766568773832,0.04830766568773832,0.1444719615827965,0.2392873622521371,0.3318686022821277,0.4213512761306353,0.5068999089322294,0.5877157572407623,0.6630442669302152,0.7321821187402897,0.7944837959679424,0.84936761373257,0.8963211557660521,0.9349060759377397,0.9647622555875064,0.9856115115452684,0.9972638618494816,
	      -0.9974246942464552,-0.9864557262306425,-0.9668229096899927,-0.9386943726111684,-0.9023167677434336,-0.8580096526765041,-0.8061623562741665,-0.7472304964495622,-0.6817319599697428,-0.610242345836379,-0.5333899047863476,-0.4518500172724507,-0.3663392577480734,-0.277609097152497,-0.1864392988279916,-0.09363106585473338,0,0.09363106585473338,0.1864392988279916,0.277609097152497,0.3663392577480734,0.4518500172724507,0.5333899047863476,0.610242345836379,0.6817319599697428,0.7472304964495622,0.8061623562741665,0.8580096526765041,0.9023167677434336,0.9386943726111684,0.9668229096899927,0.9864557262306425,0.9974246942464552,
	      -0.997571753790842,-0.9872278164063095,-0.9687082625333443,-0.9421623974051071,-0.9078096777183244,-0.8659346383345645,-0.8168842279009336,-0.761064876629873,-0.6989391132162629,-0.6310217270805285,-0.5578755006697467,-0.480106545190327,-0.3983592777586459,-0.3133110813394632,-0.2256666916164495,-0.136152357259183,-0.04550982195310254,0.04550982195310254,0.136152357259183,0.2256666916164495,0.3133110813394632,0.3983592777586459,0.480106545190327,0.5578755006697467,0.6310217270805285,0.6989391132162629,0.761064876629873,0.8168842279009336,0.8659346383345645,0.9078096777183244,0.9421623974051071,0.9687082625333443,0.9872278164063095,0.997571753790842,
	      -0.9977065690996003,-0.9879357644438514,-0.9704376160392298,-0.9453451482078273,-0.9128542613593176,-0.8732191250252224,-0.8267498990922254,-0.7738102522869126,-0.7148145015566287,-0.6502243646658904,-0.5805453447497645,-0.5063227732414887,-0.4281375415178142,-0.346601554430814,-0.262352941209296,-0.1760510611659896,-0.08837134327565926,0,0.08837134327565926,0.1760510611659896,0.262352941209296,0.346601554430814,0.4281375415178142,0.5063227732414887,0.5805453447497645,0.6502243646658904,0.7148145015566287,0.7738102522869126,0.8267498990922254,0.8732191250252224,0.9128542613593176,0.9453451482078273,0.9704376160392298,0.9879357644438514,0.9977065690996003,
	      -0.9978304624840858,-0.9885864789022122,-0.972027691049698,-0.9482729843995076,-0.9174977745156591,-0.8799298008903971,-0.8358471669924753,-0.7855762301322066,-0.7294891715935565,-0.668001236585521,-0.6015676581359806,-0.5306802859262452,-0.4558639444334203,-0.3776725471196892,-0.2966849953440283,-0.2135008923168656,-0.1287361038093848,-0.04301819847370861,0.04301819847370861,0.1287361038093848,0.2135008923168656,0.2966849953440283,0.3776725471196892,0.4558639444334203,0.5306802859262452,0.6015676581359806,0.668001236585521,0.7294891715935565,0.7855762301322066,0.8358471669924753,0.8799298008903971,0.9174977745156591,0.9482729843995076,0.972027691049698,0.9885864789022122,0.9978304624840858,
	      -0.9979445824779136,-0.9891859632143192,-0.9734930300564858,-0.9509723432620948,-0.9217814374124638,-0.8861249621554861,-0.844252987340556,-0.7964592005099023,-0.7430788339819653,-0.6844863091309593,-0.6210926084089244,-0.5533423918615817,-0.4817108778032055,-0.4067005093183261,-0.328837429883707,-0.2486677927913657,-0.166753930239852,-0.0836704089547699,0,0.0836704089547699,0.166753930239852,0.2486677927913657,0.328837429883707,0.4067005093183261,0.4817108778032055,0.5533423918615817,0.6210926084089244,0.6844863091309593,0.7430788339819653,0.7964592005099023,0.844252987340556,0.8861249621554861,0.9217814374124638,0.9509723432620948,0.9734930300564858,0.9891859632143192,0.9979445824779136,
	      -0.9980499305356876,-0.9897394542663855,-0.9748463285901535,-0.9534663309335296,-0.9257413320485844,-0.8918557390046322,-0.8520350219323621,-0.8065441676053168,-0.7556859037539707,-0.6997986803791844,-0.6392544158296817,-0.5744560210478071,-0.5058347179279311,-0.4338471694323765,-0.358972440479435,-0.2817088097901653,-0.2025704538921167,-0.1220840253378674,-0.04078514790457824,0.04078514790457824,0.1220840253378674,0.2025704538921167,0.2817088097901653,0.358972440479435,0.4338471694323765,0.5058347179279311,0.5744560210478071,0.6392544158296817,0.6997986803791844,0.7556859037539707,0.8065441676053168,0.8520350219323621,0.8918557390046322,0.9257413320485844,0.9534663309335296,0.9748463285901535,0.9897394542663855,0.9980499305356876,
	      -0.9981473830664329,-0.990251536854686,-0.9760987093334711,-0.9557752123246522,-0.9294091484867383,-0.8971671192929929,-0.8592529379999062,-0.8159062974301431,-0.7674012429310635,-0.7140444358945347,-0.656173213432011,-0.594153454957278,-0.5283772686604374,-0.4592605123091361,-0.3872401639715615,-0.3127715592481859,-0.2363255124618358,-0.1583853399978378,-0.07944380460875548,0,0.07944380460875548,0.1583853399978378,0.2363255124618358,0.3127715592481859,0.3872401639715615,0.4592605123091361,0.5283772686604374,0.594153454957278,0.656173213432011,0.7140444358945347,0.7674012429310635,0.8159062974301431,0.8592529379999062,0.8971671192929929,0.9294091484867383,0.9557752123246522,0.9760987093334711,0.990251536854686,0.9981473830664329,
	      -0.9982377097105593,-0.990726238699457,-0.9772599499837743,-0.9579168192137917,-0.9328128082786765,-0.9020988069688743,-0.8659595032122596,-0.8246122308333117,-0.7783056514265194,-0.7273182551899271,-0.6719566846141796,-0.6125538896679802,-0.5494671250951282,-0.4830758016861787,-0.413779204371605,-0.3419940908257585,-0.2681521850072537,-0.1926975807013711,-0.1160840706752552,-0.03877241750605082,0.03877241750605082,0.1160840706752552,0.1926975807013711,0.2681521850072537,0.3419940908257585,0.413779204371605,0.4830758016861787,0.5494671250951282,0.6125538896679802,0.6719566846141796,0.7273182551899271,0.7783056514265194,0.8246122308333117,0.8659595032122596,0.9020988069688743,0.9328128082786765,0.9579168192137917,0.9772599499837743,0.990726238699457,0.9982377097105593,
	      -0.9983215885747715,-0.9911671096990163,-0.9783386735610834,-0.9599068917303463,-0.9359769874978539,-0.9066859447581012,-0.8722015116924414,-0.8327212004013613,-0.7884711450474093,-0.7397048030699261,-0.6867015020349513,-0.6297648390721963,-0.5692209416102159,-0.5054165991994061,-0.4387172770514071,-0.3695050226404815,-0.2981762773418249,-0.2251396056334228,-0.1508133548639922,-0.075623258989163,0,0.075623258989163,0.1508133548639922,0.2251396056334228,0.2981762773418249,0.3695050226404815,0.4387172770514071,0.5054165991994061,0.5692209416102159,0.6297648390721963,0.6867015020349513,0.7397048030699261,0.7884711450474093,0.8327212004013613,0.8722015116924414,0.9066859447581012,0.9359769874978539,0.9599068917303463,0.9783386735610834,0.9911671096990163,0.9983215885747715,
	      -0.9983996189900625,-0.9915772883408609,-0.9793425080637482,-0.9617593653382045,-0.9389235573549882,-0.9109597249041275,-0.8780205698121727,-0.8402859832618169,-0.7979620532554874,-0.7512799356894805,-0.7004945905561712,-0.6458833888692478,-0.5877445974851093,-0.5263957499311923,-0.4621719120704219,-0.395423852042975,-0.3265161244654115,-0.2558250793428791,-0.1837368065648546,-0.1106450272085199,-0.03694894316535178,0.03694894316535178,0.1106450272085199,0.1837368065648546,0.2558250793428791,0.3265161244654115,0.395423852042975,0.4621719120704219,0.5263957499311923,0.5877445974851093,0.6458833888692478,0.7004945905561712,0.7512799356894805,0.7979620532554874,0.8402859832618169,0.8780205698121727,0.9109597249041275,0.9389235573549882,0.9617593653382045,0.9793425080637482,0.9915772883408609,0.9983996189900625,
	      -0.9984723322425078,-0.9919595575932442,-0.9802782209802553,-0.96348661301408,-0.9416719568476378,-0.9149479072061387,-0.8834537652186168,-0.8473537162093151,-0.8068359641369386,-0.7621117471949551,-0.7134142352689571,-0.6609973137514982,-0.605134259639601,-0.5461163166600848,-0.4842511767857347,-0.4198613760292693,-0.3532826128643038,-0.2848619980329136,-0.2149562448605182,-0.1439298095107133,-0.07215299087458624,0,0.07215299087458624,0.1439298095107133,0.2149562448605182,0.2848619980329136,0.3532826128643038,0.4198613760292693,0.4842511767857347,0.5461163166600848,0.605134259639601,0.6609973137514982,0.7134142352689571,0.7621117471949551,0.8068359641369386,0.8473537162093151,0.8834537652186168,0.9149479072061387,0.9416719568476378,0.96348661301408,0.9802782209802553,0.9919595575932442,0.9984723322425078,
	      -0.9985402006367742,-0.9923163921385159,-0.981151833077914,-0.9650996504224931,-0.9442395091181941,-0.9186752599841758,-0.8885342382860432,-0.8539665950047104,-0.815144539645135,-0.7722614792487559,-0.725531053660717,-0.6751860706661224,-0.6214773459035758,-0.5646724531854708,-0.5050543913882023,-0.4429201745254115,-0.3785793520147071,-0.3123524665027858,-0.2445694569282013,-0.1755680147755168,-0.1056919017086533,-0.03528923696413536,0.03528923696413536,0.1056919017086533,0.1755680147755168,0.2445694569282013,0.3123524665027858,0.3785793520147071,0.4429201745254115,0.5050543913882023,0.5646724531854708,0.6214773459035758,0.6751860706661224,0.725531053660717,0.7722614792487559,0.815144539645135,0.8539665950047104,0.8885342382860432,0.9186752599841758,0.9442395091181941,0.9650996504224931,0.981151833077914,0.9923163921385159,0.9985402006367742,
	      -0.9986036451819367,-0.9926499984472037,-0.9819687150345405,-0.9666083103968947,-0.9466416909956291,-0.9221639367190004,-0.8932916717532418,-0.8601624759606642,-0.8229342205020863,-0.7817843125939062,-0.7369088489454904,-0.6885216807712006,-0.6368533944532233,-0.5821502125693532,-0.5246728204629161,-0.4646951239196351,-0.4025029438585419,-0.3383926542506022,-0.2726697697523776,-0.2056474897832637,-0.137645205983253,-0.06898698016314417,0,0.06898698016314417,0.137645205983253,0.2056474897832637,0.2726697697523776,0.3383926542506022,0.4025029438585419,0.4646951239196351,0.5246728204629161,0.5821502125693532,0.6368533944532233,0.6885216807712006,0.7369088489454904,0.7817843125939062,0.8229342205020863,0.8601624759606642,0.8932916717532418,0.9221639367190004,0.9466416909956291,0.9666083103968947,0.9819687150345405,0.9926499984472037,0.9986036451819367,
	      -0.998663042133818,-0.9929623489061744,-0.9827336698041669,-0.9680213918539919,-0.9488923634460898,-0.9254337988067539,-0.897752711533942,-0.865975394866858,-0.8302468370660661,-0.7907300570752742,-0.7476053596156661,-0.7010695120204057,-0.6513348462019977,-0.5986282897127152,-0.5431903302618026,-0.4852739183881646,-0.4251433132828284,-0.3630728770209957,-0.29934582270187,-0.2342529222062698,-0.1680911794671035,-0.1011624753055842,-0.03377219001605204,0.03377219001605204,0.1011624753055842,0.1680911794671035,0.2342529222062698,0.29934582270187,0.3630728770209957,0.4251433132828284,0.4852739183881646,0.5431903302618026,0.5986282897127152,0.6513348462019977,0.7010695120204057,0.7476053596156661,0.7907300570752742,0.8302468370660661,0.865975394866858,0.897752711533942,0.9254337988067539,0.9488923634460898,0.9680213918539919,0.9827336698041669,0.9929623489061744,0.998663042133818,
	      -0.9987187285842121,-0.9932552109877686,-0.9834510030716237,-0.9693467873265645,-0.9510039692577085,-0.9285026930123607,-0.9019413294385253,-0.8714360157968963,-0.8371201398999021,-0.799143754167742,-0.7576729184454386,-0.7128889734090643,-0.6649877473903327,-0.6141786999563736,-0.5606840059346642,-0.5047375838635779,-0.4465840731048557,-0.3864777640846672,-0.3246814863377359,-0.2614654592149745,-0.1971061102791118,-0.1318848665545149,-0.06608692391635568,0,0.06608692391635568,0.1318848665545149,0.1971061102791118,0.2614654592149745,0.3246814863377359,0.3864777640846672,0.4465840731048557,0.5047375838635779,0.5606840059346642,0.6141786999563736,0.6649877473903327,0.7128889734090643,0.7576729184454386,0.799143754167742,0.8371201398999021,0.8714360157968963,0.9019413294385253,0.9285026930123607,0.9510039692577085,0.9693467873265645,0.9834510030716237,0.9932552109877686,0.9987187285842121,
	      -0.9987710072524261,-0.9935301722663508,-0.9841245837228269,-0.9705915925462473,-0.9529877031604309,-0.9313866907065543,-0.9058791367155696,-0.8765720202742479,-0.8435882616243935,-0.8070662040294426,-0.7671590325157404,-0.7240341309238146,-0.6778723796326639,-0.6288673967765136,-0.5772247260839727,-0.523160974722233,-0.4669029047509584,-0.4086864819907167,-0.3487558862921608,-0.2873624873554556,-0.2247637903946891,-0.1612223560688917,-0.0970046992094627,-0.03238017096286936,0.03238017096286936,0.0970046992094627,0.1612223560688917,0.2247637903946891,0.2873624873554556,0.3487558862921608,0.4086864819907167,0.4669029047509584,0.523160974722233,0.5772247260839727,0.6288673967765136,0.6778723796326639,0.7240341309238146,0.7671590325157404,0.8070662040294426,0.8435882616243935,0.8765720202742479,0.9058791367155696,0.9313866907065543,0.9529877031604309,0.9705915925462473,0.9841245837228269,0.9935301722663508,0.9987710072524261,
	      -0.9988201506066354,-0.9937886619441678,-0.984757895914213,-0.9717622009015554,-0.9548536586741372,-0.9341002947558101,-0.9095856558280733,-0.8814084455730089,-0.8496821198441658,-0.8145344273598555,-0.7761068943454467,-0.7345542542374027,-0.6900438244251321,-0.6427548324192377,-0.5928776941089007,-0.5406132469917261,-0.486171941452492,-0.4297729933415765,-0.3716435012622849,-0.3120175321197488,-0.2511351786125773,-0.1892415924618136,-0.126585997269672,-0.06342068498268678,0,0.06342068498268678,0.126585997269672,0.1892415924618136,0.2511351786125773,0.3120175321197488,0.3716435012622849,0.4297729933415765,0.486171941452492,0.5406132469917261,0.5928776941089007,0.6427548324192377,0.6900438244251321,0.7345542542374027,0.7761068943454467,0.8145344273598555,0.8496821198441658,0.8814084455730089,0.9095856558280733,0.9341002947558101,0.9548536586741372,0.9717622009015554,0.984757895914213,0.9937886619441678,0.9988201506066354,
	      -0.998866404420071,-0.9940319694320907,-0.9853540840480058,-0.972864385106692,-0.9566109552428079,-0.936656618944878,-0.9130785566557919,-0.8859679795236131,-0.8554297694299461,-0.821582070859336,-0.7845558329003993,-0.7444943022260685,-0.7015524687068222,-0.6558964656854394,-0.6077029271849502,-0.5571583045146501,-0.5044581449074642,-0.4498063349740388,-0.3934143118975651,-0.3355002454194373,-0.276288193779532,-0.2160072368760418,-0.1548905899981459,-0.09317470156008614,-0.03109833832718887,0.03109833832718887,0.09317470156008614,0.1548905899981459,0.2160072368760418,0.276288193779532,0.3355002454194373,0.3934143118975651,0.4498063349740388,0.5044581449074642,0.5571583045146501,0.6077029271849502,0.6558964656854394,0.7015524687068222,0.7444943022260685,0.7845558329003993,0.821582070859336,0.8554297694299461,0.8859679795236131,0.9130785566557919,0.936656618944878,0.9566109552428079,0.972864385106692,0.9853540840480058,0.9940319694320907,0.998866404420071,
	      -0.9993050417357722,-0.9963401167719553,-0.9910133714767443,-0.983336253884626,-0.973326827789911,-0.9610087996520538,-0.9464113748584028,-0.9295691721319396,-0.9105221370785028,-0.8893154459951141,-0.8659993981540928,-0.8406292962525803,-0.8132653151227975,-0.7839723589433414,-0.7528199072605319,-0.7198818501716108,-0.6852363130542333,-0.6489654712546573,-0.6111553551723933,-0.571895646202634,-0.5312794640198946,-0.489403145707053,-0.4463660172534641,-0.4022701579639916,-0.3572201583376681,-0.311322871990211,-0.2646871622087674,-0.2174236437400071,-0.1696444204239928,-0.1214628192961206,-0.07299312178779904,-0.02435029266342443,0.02435029266342443,0.07299312178779904,0.1214628192961206,0.1696444204239928,0.2174236437400071,0.2646871622087674,0.311322871990211,0.3572201583376681,0.4022701579639916,0.4463660172534641,0.489403145707053,0.5312794640198946,0.571895646202634,0.6111553551723933,0.6489654712546573,0.6852363130542333,0.7198818501716108,0.7528199072605319,0.7839723589433414,0.8132653151227975,0.8406292962525803,0.8659993981540928,0.8893154459951141,0.9105221370785028,0.9295691721319396,0.9464113748584028,0.9610087996520538,0.973326827789911,0.983336253884626,0.9910133714767443,0.9963401167719553,0.9993050417357722,
	      -0.9998248879471319,-0.9990774599773758,-0.997733248625514,-0.9957927585349812,-0.9932571129002129,-0.9901278184917344,-0.9864067427245862,-0.9820961084357185,-0.9771984914639074,-0.9717168187471366,-0.9656543664319652,-0.9590147578536999,-0.9518019613412644,-0.9440202878302202,-0.9356743882779164,-0.9267692508789478,-0.9173101980809605,-0.9073028834017568,-0.8967532880491582,-0.8856677173453972,-0.8740527969580318,-0.8619154689395485,-0.8492629875779689,-0.8361029150609068,-0.8224431169556439,-0.8082917575079137,-0.7936572947621933,-0.7785484755064119,-0.7629743300440948,-0.746944166797062,-0.7304675667419088,-0.7135543776835874,-0.6962147083695144,-0.6784589224477192,-0.660297632272646,-0.6417416925623075,-0.6228021939105849,-0.6034904561585486,-0.5838180216287631,-0.5637966482266181,-0.5434383024128103,-0.5227551520511755,-0.5017595591361445,-0.480464072404172,-0.4588814198335522,-0.4370245010371042,-0.414906379552275,-0.3925402750332674,-0.369939555349859,-0.3471177285976355,-0.3240884350244134,-0.3008654388776772,-0.2774626201779044,-0.2538939664226943,-0.23017356422666,-0.2063155909020792,-0.1823343059853372,-0.1582440427142249,-0.1340591994611878,-0.1097942311276437,-0.08546364050451549,-0.06108196960413957,-0.0366637909687335,-0.01222369896061576,0.01222369896061576,0.0366637909687335,0.06108196960413957,0.08546364050451549,0.1097942311276437,0.1340591994611878,0.1582440427142249,0.1823343059853372,0.2063155909020792,0.23017356422666,0.2538939664226943,0.2774626201779044,0.3008654388776772,0.3240884350244134,0.3471177285976355,0.369939555349859,0.3925402750332674,0.414906379552275,0.4370245010371042,0.4588814198335522,0.480464072404172,0.5017595591361445,0.5227551520511755,0.5434383024128103,0.5637966482266181,0.5838180216287631,0.6034904561585486,0.6228021939105849,0.6417416925623075,0.660297632272646,0.6784589224477192,0.6962147083695144,0.7135543776835874,0.7304675667419088,0.746944166797062,0.7629743300440948,0.7785484755064119,0.7936572947621933,0.8082917575079137,0.8224431169556439,0.8361029150609068,0.8492629875779689,0.8619154689395485,0.8740527969580318,0.8856677173453972,0.8967532880491582,0.9073028834017568,0.9173101980809605,0.9267692508789478,0.9356743882779164,0.9440202878302202,0.9518019613412644,0.9590147578536999,0.9656543664319652,0.9717168187471366,0.9771984914639074,0.9820961084357185,0.9864067427245862,0.9901278184917344,0.9932571129002129,0.9957927585349812,0.997733248625514,0.9990774599773758,0.9998248879471319,
	      -0.9999560500189922,-0.9997684374092631,-0.9994309374662614,-0.9989435258434088,-0.9983062664730065,-0.9975192527567208,-0.9965826020233816,-0.9954964544810964,-0.9942609729224097,-0.9928763426088221,-0.9913427712075831,-0.9896604887450652,-0.9878297475648606,-0.985850822286126,-0.9837240097603155,-0.9814496290254644,-0.979028021257622,-0.9764595497192342,-0.9737445997043704,-0.970883578480743,-0.9678769152284895,-0.9647250609757064,-0.9614284885307322,-0.9579876924111781,-0.9544031887697162,-0.9506755153166283,-0.9468052312391275,-0.9427929171174625,-0.9386391748378148,-0.9343446275020031,-0.9299099193340057,-0.9253357155833162,-0.9206227024251465,-0.9157715868574904,-0.9107830965950651,-0.9056579799601446,-0.9003970057703036,-0.8950009632230845,-0.8894706617776109,-0.8838069310331583,-0.8780106206047066,-0.8720825999954883,-0.8660237584665545,-0.8598350049033764,-0.853517267679503,-0.8470714945172962,-0.8404986523457627,-0.8337997271555049,-0.8269757238508125,-0.8200276660989171,-0.8129565961764316,-0.8057635748129987,-0.7984496810321707,-0.791016011989546,-0.7834636828081838,-0.7757938264113258,-0.7680075933524456,-0.7601061516426555,-0.7520906865754921,-0.7439624005491116,-0.7357225128859178,-0.7273722596496521,-0.7189128934599714,-0.7103456833045433,-0.7016719143486851,-0.692892887742577,-0.684009920426076,-0.6750243449311628,-0.6659375091820485,-0.6567507762929732,-0.6474655243637248,-0.6380831462729114,-0.628605049469015,-0.6190326557592613,-0.6093674010963339,-0.5996107353629683,-0.5897641221544543,-0.579829038559083,-0.5698069749365687,-0.5596994346944811,-0.5495079340627186,-0.5392340018660592,-0.5288791792948223,-0.5184450196736745,-0.507933088228616,-0.4973449618521815,-0.4866822288668903,-0.4759464887869833,-0.4651393520784793,-0.45426243991759,-0.4433173839475273,-0.4323058260337413,-0.4212294180176238,-0.4100898214687165,-0.3988887074354591,-0.3876277561945156,-0.3763086569987164,-0.364933107823654,-0.35350281511297,-0.3420194935223717,-0.330484865662417,-0.3189006618401063,-0.3072686197993191,-0.2955904844601356,-0.2838680076570818,-0.2721029478763366,-0.2602970699919425,-0.2484521450010567,-0.236569949758284,-0.224652266709132,-0.212700883622626,-0.2007175933231267,-0.1887041934213888,-0.176662486044902,-0.1645942775675538,-0.1525013783386564,-0.1403856024113759,-0.1282487672706071,-0.1160926935603328,-0.1039192048105094,-0.09173012716351955,-0.07952728910023296,-0.0673125211657164,-0.05508765569463399,-0.0428545265363791,-0.03061496877997903,-0.01837081847881367,-0.00612391237518953,0.00612391237518953,0.01837081847881367,0.03061496877997903,0.0428545265363791,0.05508765569463399,0.0673125211657164,0.07952728910023296,0.09173012716351955,0.1039192048105094,0.1160926935603328,0.1282487672706071,0.1403856024113759,0.1525013783386564,0.1645942775675538,0.176662486044902,0.1887041934213888,0.2007175933231267,0.212700883622626,0.224652266709132,0.236569949758284,0.2484521450010567,0.2602970699919425,0.2721029478763366,0.2838680076570818,0.2955904844601356,0.3072686197993191,0.3189006618401063,0.330484865662417,0.3420194935223717,0.35350281511297,0.364933107823654,0.3763086569987164,0.3876277561945156,0.3988887074354591,0.4100898214687165,0.4212294180176238,0.4323058260337413,0.4433173839475273,0.45426243991759,0.4651393520784793,0.4759464887869833,0.4866822288668903,0.4973449618521815,0.507933088228616,0.5184450196736745,0.5288791792948223,0.5392340018660592,0.5495079340627186,0.5596994346944811,0.5698069749365687,0.579829038559083,0.5897641221544543,0.5996107353629683,0.6093674010963339,0.6190326557592613,0.628605049469015,0.6380831462729114,0.6474655243637248,0.6567507762929732,0.6659375091820485,0.6750243449311628,0.684009920426076,0.692892887742577,0.7016719143486851,0.7103456833045433,0.7189128934599714,0.7273722596496521,0.7357225128859178,0.7439624005491116,0.7520906865754921,0.7601061516426555,0.7680075933524456,0.7757938264113258,0.7834636828081838,0.791016011989546,0.7984496810321707,0.8057635748129987,0.8129565961764316,0.8200276660989171,0.8269757238508125,0.8337997271555049,0.8404986523457627,0.8470714945172962,0.853517267679503,0.8598350049033764,0.8660237584665545,0.8720825999954883,0.8780106206047066,0.8838069310331583,0.8894706617776109,0.8950009632230845,0.9003970057703036,0.9056579799601446,0.9107830965950651,0.9157715868574904,0.9206227024251465,0.9253357155833162,0.9299099193340057,0.9343446275020031,0.9386391748378148,0.9427929171174625,0.9468052312391275,0.9506755153166283,0.9544031887697162,0.9579876924111781,0.9614284885307322,0.9647250609757064,0.9678769152284895,0.970883578480743,0.9737445997043704,0.9764595497192342,0.979028021257622,0.9814496290254644,0.9837240097603155,0.985850822286126,0.9878297475648606,0.9896604887450652,0.9913427712075831,0.9928763426088221,0.9942609729224097,0.9954964544810964,0.9965826020233816,0.9975192527567208,0.9983062664730065,0.9989435258434088,0.9994309374662614,0.9997684374092631,0.9999560500189922,
	    };

	  unsigned order_array[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,64,128,256,};
	  unsigned order_array_length = 53;

	  // Construct the weights or knots as requested
	  if(weights)
	    construction_helper(weights_data_array, order_array, order_array_length);
	  else
	    construction_helper(knots_data_array, order_array, order_array_length);

	  break;
	}


      case 1 :
	{

	  double weights_data_array[] =
	    {
	      2,
	      1,1,
	      0.3333333333333333,1.333333333333333,0.3333333333333333,
	      0.1111111111111111,0.8888888888888888,0.8888888888888892,0.1111111111111111,
	      0.06666666666666667,0.5333333333333333,0.8,0.5333333333333334,0.06666666666666667,
	      0.03999999999999999,0.3607430412000112,0.5992569587999889,0.5992569587999887,0.3607430412000113,0.03999999999999999,
	      0.02857142857142857,0.253968253968254,0.4571428571428571,0.5206349206349207,0.4571428571428573,0.253968253968254,0.02857142857142857,
	      0.02040816326530612,0.1901410072182084,0.3522424237181591,0.4372084057983264,0.4372084057983264,0.3522424237181592,0.1901410072182084,0.02040816326530612,
	      0.01587301587301587,0.1462186492160182,0.2793650793650794,0.3617178587204897,0.3936507936507937,0.3617178587204898,0.2793650793650794,0.1462186492160182,0.01587301587301587,
	      0.01234567901234568,0.1165674565720371,0.2252843233381044,0.3019400352733685,0.3438625058041442,0.3438625058041441,0.3019400352733687,0.2252843233381044,0.1165674565720372,0.01234567901234568,
	      0.0101010101010101,0.09457905488370155,0.1856352144242478,0.2535883332836866,0.2992132704242371,0.3137662337662338,0.2992132704242371,0.2535883332836866,0.1856352144242478,0.09457905488370161,0.0101010101010101,
	      0.008264462809917354,0.0785601537462,0.1550404550825614,0.2155625460008686,0.2599173410669162,0.2826550412935365,0.2826550412935365,0.2599173410669162,0.2155625460008686,0.1550404550825614,0.0785601537462,0.008264462809917354,
	      0.006993006993006992,0.06605742495207438,0.1315425315425316,0.1847633847633847,0.2269730269730269,0.2526756937810443,0.261989861989862,0.2526756937810444,0.226973026973027,0.1847633847633848,0.1315425315425316,0.06605742495207449,0.006993006993006992,
	      0.005917159763313608,0.05646531376341445,0.1127686724898565,0.1600380261167187,0.1989924103657832,0.2259030497785645,0.239915367722349,0.239915367722349,0.2259030497785645,0.1989924103657832,0.1600380261167187,0.1127686724898566,0.05646531376341445,0.005917159763313608,
	      0.005128205128205128,0.04869938729508823,0.09782039167605217,0.1396650784956043,0.1756057890010667,0.2020514674823835,0.2188815116305734,0.2242963385820529,0.2188815116305734,0.2020514674823836,0.1756057890010668,0.1396650784956043,0.09782039167605218,0.04869938729508826,0.005128205128205128,
	      0.004444444444444444,0.04251476624752508,0.08553884025933287,0.1229401008284936,0.1557331760396737,0.1813297813297813,0.1992147813263886,0.2082841095243604,0.2082841095243604,0.1992147813263886,0.1813297813297813,0.1557331760396737,0.1229401008284936,0.08553884025933292,0.04251476624752508,0.004444444444444444,
	      0.003921568627450979,0.0373687028372056,0.07548233154315184,0.1089055525818909,0.1389564683682331,0.1631726642817033,0.1814737842364934,0.1925138646129257,0.1964101258218905,0.1925138646129257,0.1814737842364934,0.1631726642817033,0.1389564683682331,0.1089055525818909,0.07548233154315184,0.0373687028372057,0.003921568627450979,
	      0.003460207612456747,0.0331534326715737,0.06703353832804583,0.09715305543612526,0.1245714688992924,0.147430009111827,0.1654514950211523,0.1777312087202402,0.1840155841992866,0.1840155841992867,0.1777312087202401,0.1654514950211523,0.147430009111827,0.1245714688992924,0.09715305543612529,0.0670335383280458,0.0331534326715737,0.003460207612456747,
	      0.003095975232198142,0.02957099501095232,0.05996056410917944,0.0871279467770696,0.1122751964454378,0.1336369232722281,0.1511998837076236,0.1639645996014117,0.1719203928894619,0.1744950459088746,0.1719203928894619,0.1639645996014117,0.1511998837076237,0.1336369232722281,0.1122751964454378,0.0871279467770696,0.05996056410917946,0.02957099501095231,0.003095975232198142,
	      0.002770083102493074,0.0265712911416073,0.05391036208922058,0.07858988043704006,0.1016181939319847,0.1216130604542333,0.1384460983072959,0.1514060491721975,0.1602952408066014,0.1647797405573261,0.1647797405573261,0.1602952408066014,0.1514060491721975,0.138446098307296,0.1216130604542334,0.1016181939319847,0.07858988043704007,0.05391036208922061,0.02657129114160733,0.002770083102493074,
	      0.0025062656641604,0.02397864723964093,0.04875412910691163,0.07120273882096263,0.09239763039561463,0.1110236443982574,0.1271175533901578,0.1399278661926565,0.1494190324736799,0.1551202361805628,0.1571045122747909,0.1551202361805628,0.1494190324736799,0.1399278661926565,0.1271175533901578,0.1110236443982574,0.09239763039561463,0.07120273882096265,0.04875412910691165,0.02397864723964094,0.0025062656641604,
	      0.002267573696145124,0.02176910382553372,0.04427790071214075,0.06482079940835253,0.08432252492527317,0.1017224714321703,0.1169820069087314,0.1295437375048167,0.139266874070554,0.1458441499013008,0.1491828576149814,0.1491828576149814,0.1458441499013008,0.139266874070554,0.1295437375048167,0.1169820069087314,0.1017224714321703,0.08432252492527323,0.06482079940835254,0.04427790071214075,0.02176910382553377,0.002267573696145124,
	      0.0020703933747412,0.01983305181710098,0.04040669187862778,0.05923164660738842,0.07725854407420041,0.09347540554706026,0.1079500539305165,0.1201069225701775,0.1299155986669404,0.1369968405793485,0.1413635213876031,0.14278265913259,0.1413635213876031,0.1369968405793485,0.1299155986669404,0.1201069225701775,0.1079500539305165,0.09347540554706027,0.07725854407420042,0.05923164660738844,0.04040669187862778,0.01983305181710098,0.0020703933747412,
	      0.001890359168241965,0.01815905012388155,0.0370044679933253,0.05434409896562557,0.0710127642718389,0.08617432027228909,0.09984327978460976,0.1115786801022364,0.1212852892186523,0.1286990494367902,0.1337370463803083,0.1362715942822006,0.1362715942822006,0.1337370463803083,0.1286990494367902,0.1212852892186523,0.1115786801022364,0.09984327978460975,0.08617432027228915,0.07101276427183893,0.05434409896562557,0.03700446799332529,0.01815905012388157,0.001890359168241965,
	      0.001739130434782608,0.01667549216107983,0.03402583454591183,0.05001881853170057,0.06549535317337175,0.07965528767324662,0.09258363160099597,0.1038308385866787,0.1133783755617111,0.1209215189559515,0.1264522547409,0.1297676093087341,0.1309117094498707,0.1297676093087342,0.1264522547409,0.1209215189559515,0.1133783755617111,0.1038308385866787,0.09258363160099599,0.07965528767324663,0.06549535317337175,0.05001881853170057,0.03402583454591188,0.01667549216107982,0.001739130434782608,
	      0.001599999999999999,0.0153771669147882,0.03138077102269942,0.04619703562462628,0.06057509265584682,0.07384035495765288,0.08603814579548864,0.09681482453872046,0.1061089219827171,0.1136986829074585,0.119516734736279,0.1234358939129451,0.1254163749507777,0.1254163749507777,0.1234358939129451,0.1195167347362789,0.1136986829074585,0.1061089219827171,0.09681482453872051,0.08603814579548869,0.07384035495765287,0.0605750926558468,0.04619703562462631,0.03138077102269941,0.01537716691478825,0.001599999999999999,
	      0.001481481481481481,0.01421548079893025,0.02904088677824307,0.0427843154260888,0.05619092210016806,0.06861339505972898,0.08014486808181145,0.0904275436091969,0.09945432378545657,0.1069790025753655,0.1129884263378799,0.1173101238696938,0.1199583506942187,0.1208217588034729,0.1199583506942187,0.1173101238696938,0.1129884263378799,0.1069790025753655,0.09945432378545657,0.09042754360919693,0.08014486808181145,0.06861339505972898,0.05619092210016808,0.04278431542608882,0.02904088677824306,0.01421548079893026,0.001481481481481481,
	      0.001371742112482853,0.0131883954457362,0.02694455907239935,0.03974175146473947,0.05225137230908839,0.06391836544104249,0.07480524531264111,0.08462423974649817,0.09333812302591335,0.1007616931694246,0.1068427831471891,0.1114648176141417,0.1145889177859132,0.11615799435279,0.11615799435279,0.1145889177859132,0.1114648176141417,0.1068427831471892,0.1007616931694246,0.09333812302591336,0.0846242397464982,0.07480524531264113,0.0639183654410425,0.0522513723090884,0.03974175146473948,0.02694455907239936,0.01318839544573621,0.001371742112482853,
	      0.00127713920817369,0.01226186308296506,0.02507354294987612,0.03700348812448521,0.04871418451399866,0.05967172250029153,0.06997172018812491,0.07932439914181021,0.08773180864766075,0.09499323112767082,0.101096440098399,0.1058961889907875,0.1093934802640209,0.1114876766360766,0.1122062290513182,0.1114876766360766,0.1093934802640209,0.1058961889907875,0.101096440098399,0.09499323112767082,0.08773180864766074,0.07932439914181022,0.06997172018812491,0.05967172250029155,0.04871418451399866,0.03700348812448522,0.02507354294987613,0.01226186308296508,0.00127713920817369,
	      0.001189060642092746,0.01143547089105791,0.02338441509783392,0.03454314798888627,0.04551413072927685,0.05583292240186365,0.06557089935125514,0.07449054615558812,0.08257179430862692,0.08965963475345294,0.09571489906131193,0.1006343759618618,0.1043837360991885,0.1069027576183408,0.1081722089393624,0.1081722089393624,0.1069027576183408,0.1043837360991885,0.1006343759618618,0.09571489906131195,0.08965963475345297,0.08257179430862691,0.07449054615558813,0.06557089935125512,0.05583292240186363,0.04551413072927688,0.0345431479888863,0.02338441509783393,0.01143547089105793,0.001189060642092746,
	      0.001112347052280311,0.01068470435076317,0.02186518228117956,0.03231364162958044,0.04262143874002459,0.0523408741616434,0.06156665745800698,0.07006035806408085,0.07783077467610745,0.0847125838726641,0.09069650999620839,0.09566036749407249,0.09959984572594252,0.1024262301345836,0.1041510705441101,0.1047148276375041,0.1041510705441101,0.1024262301345836,0.09959984572594252,0.09566036749407252,0.0906965099962084,0.0847125838726641,0.07783077467610745,0.0700603580640809,0.06156665745800698,0.0523408741616434,0.04262143874002462,0.03231364162958046,0.02186518228117956,0.0106847043507632,0.001112347052280311,
	      0.001040582726326742,0.01000996986216579,0.02048446989741451,0.03029670497981701,0.03998861970252621,0.04916595205053239,0.05790384022127524,0.06600363286534564,0.07345730610119647,0.08013438127226709,0.08600619097229349,0.09098269668990058,0.09503508764806513,0.0981054314967449,0.1001737244866966,0.1012114090274321,0.1012114090274321,0.1001737244866967,0.09810543149674489,0.09503508764806516,0.0909826966899006,0.08600619097229351,0.08013438127226709,0.07345730610119647,0.06600363286534564,0.05790384022127525,0.0491659520505324,0.03998861970252622,0.03029670497981701,0.02048446989741453,0.0100099698621658,0.001040582726326742,
	      0.0009775171065493642,0.009393197962955013,0.01923424513268115,0.02845791667723369,0.03759434191404721,0.04626276283775175,0.05455501630398032,0.06227210954529399,0.06942757563043546,0.07588380044138848,0.0816348176549385,0.08657753844182743,0.09070611286772098,0.09394324443876872,0.09629232594548819,0.09769818820805558,0.09817857778176831,0.09769818820805558,0.09629232594548819,0.09394324443876871,0.090706112867721,0.08657753844182743,0.0816348176549385,0.0758838004413885,0.06942757563043547,0.06227210954529399,0.05455501630398031,0.04626276283775177,0.03759434191404721,0.0284579166772337,0.01923424513268119,0.009393197962955046,0.0009775171065493642,
	      0.0009182736455463723,0.00883518569074398,0.01809133857568153,0.02678439253771144,0.03540337665560393,0.04360941153622287,0.05147828712084552,0.05884182435496416,0.06570039894560709,0.07194338229704955,0.07755036084272689,0.08244320121196332,0.0865984408533037,0.08996264275683918,0.09251708033632246,0.09423029462904924,0.0950921080098187,0.0950921080098187,0.09423029462904925,0.09251708033632246,0.08996264275683918,0.08659844085330372,0.08244320121196333,0.07755036084272687,0.07194338229704955,0.06570039894560709,0.05884182435496418,0.05147828712084553,0.04360941153622287,0.03540337665560393,0.02678439253771143,0.01809133857568155,0.008835185690744009,0.0009182736455463723,
	      0.0008658008658008654,0.008322339561696946,0.01705028699490804,0.02525047534335382,0.0334001991537292,0.04117165271132335,0.04865273085293336,0.05567522316011321,0.062255837304864,0.06827883819216805,0.07374128471756851,0.07855576188428735,0.08271661439191831,0.0861568011372818,0.08887570900172628,0.09082341055245403,0.092008636283651,0.0923967957804438,0.09200863628365101,0.09082341055245403,0.08887570900172628,0.08615680113728177,0.08271661439191831,0.07855576188428735,0.07374128471756852,0.06827883819216805,0.06225583730486403,0.05567522316011323,0.04865273085293337,0.04117165271132335,0.03340019915372918,0.02525047534335384,0.01705028699490805,0.008322339561696975,0.0008658008658008654,
	      0.0008163265306122444,0.007855617388824603,0.01609364878092147,0.02384692501092013,0.03155822857402767,0.03893311572526589,0.04604615612196809,0.05275338692855699,0.05906099325581413,0.06487444199575146,0.07017971496057955,0.07490903400617871,0.07904368405938869,0.08253539095828089,0.08536753716441627,0.0875087050489908,0.08894780457950441,0.0896692889099981,0.08966928890999809,0.08894780457950439,0.0875087050489908,0.08536753716441627,0.08253539095828089,0.07904368405938871,0.0749090340061787,0.07017971496057954,0.06487444199575146,0.05906099325581415,0.05275338692855699,0.04604615612196809,0.0389331157252659,0.03155822857402767,0.02384692501092014,0.01609364878092149,0.007855617388824614,0.0008163265306122444,
	      0.0007722007722007717,0.00742460809294238,0.01521766273242866,0.02255428181550788,0.02986600309004261,0.03686759591165097,0.04364249860808216,0.05004688737975094,0.05609948383066711,0.06170199409736771,0.0668540152690916,0.07148097316851729,0.07557813444949971,0.07908738538747968,0.08200634058119245,0.08429038306856602,0.08594314004844052,0.08693199146431756,0.08726884046450803,0.08693199146431757,0.0859431400484405,0.08429038306856604,0.08200634058119244,0.07908738538747968,0.07557813444949972,0.0714809731685173,0.0668540152690916,0.06170199409736773,0.05609948383066712,0.05004688737975096,0.04364249860808216,0.03686759591165098,0.02986600309004261,0.02255428181550791,0.01521766273242867,0.00742460809294238,0.0007722007722007717,
	      0.0007304601899196489,0.007030307490142718,0.01440897402783651,0.02136582167267539,0.02830315666985924,0.03496271344610073,0.04141674583844882,0.04754078260580741,0.05334497666154384,0.05874791236371663,0.06374048780787137,0.06826377267880682,0.07230308842067565,0.07581527243397357,0.07878603528378204,0.0811854395431737,0.08300263944820513,0.0842199185369909,0.08483149488046986,0.08483149488046987,0.0842199185369909,0.08300263944820514,0.0811854395431737,0.07878603528378207,0.07581527243397357,0.07230308842067565,0.06826377267880682,0.06374048780787139,0.05874791236371666,0.05334497666154388,0.04754078260580742,0.04141674583844881,0.03496271344610073,0.02830315666985926,0.0213658216726754,0.01440897402783651,0.007030307490142751,0.0007304601899196489,
	      0.0006930006930006925,0.006664614534771629,0.01366496923560734,0.02026655410307393,0.02686111549852122,0.03319815086386309,0.03935646756226598,0.04521182719451029,0.05078412065968924,0.05598909494550464,0.06082837904862815,0.06523793952757956,0.0692141991637497,0.0727066942911031,0.07571245849575388,0.07819187368190007,0.08014574154214719,0.08154396361782708,0.08239304775413625,0.08267157517273392,0.08239304775413625,0.08154396361782705,0.08014574154214721,0.07819187368190009,0.07571245849575389,0.07270669429110312,0.06921419916374971,0.06523793952757956,0.06082837904862818,0.05598909494550464,0.05078412065968927,0.04521182719451029,0.03935646756226598,0.03319815086386307,0.02686111549852124,0.02026655410307392,0.01366496923560735,0.006664614534771622,0.0006930006930006925,
	      0.0006574621959237339,0.006328498989942889,0.01297526669525893,0.01925152226104752,0.02552405657609022,0.03156444437104093,0.03744201765981547,0.04304818211599316,0.0483957464610764,0.05341405601229049,0.05809776985320771,0.06239545923473512,0.06629579804896743,0.06976037025701511,0.07277713895609879,0.07531828155108194,0.07737376561643511,0.07892536261045552,0.07996654128321905,0.08048825925030435,0.08048825925030435,0.07996654128321905,0.0789253626104555,0.07737376561643512,0.07531828155108193,0.07277713895609879,0.06976037025701512,0.06629579804896743,0.06239545923473512,0.05809776985320771,0.05341405601229048,0.04839574646107642,0.04304818211599318,0.03744201765981549,0.03156444437104093,0.02552405657609022,0.01925152226104752,0.01297526669525896,0.006328498989942911,0.0006574621959237339,
	      0.0006253908692933079,0.006015567729891233,0.01233803907178104,0.01830901591261394,0.02428548929466838,0.03004564840197724,0.03566369984754822,0.0410313326110132,0.04616900147467561,0.05100386712698057,0.05553913518277631,0.05971946032639461,0.06354242126631117,0.06696406123651212,0.06998145072470656,0.0725595027997117,0.07469752809300673,0.0763681051519477,0.0775744602879448,0.07829613413760424,0.07854137690528259,0.07829613413760426,0.0775744602879448,0.07636810515194771,0.07469752809300674,0.07255950279971171,0.06998145072470656,0.06696406123651213,0.06354242126631117,0.05971946032639462,0.05553913518277633,0.05100386712698057,0.04616900147467561,0.0410313326110132,0.03566369984754823,0.03004564840197725,0.02428548929466839,0.01830901591261396,0.01233803907178104,0.006015567729891242,0.0006253908692933079,
	      0.0005948839976204635,0.005726732928160036,0.01174510563951802,0.01743533888979566,0.02313296372701732,0.02863445416487101,0.03400592341164632,0.03915185350905824,0.04408677229977931,0.04874891709851373,0.05313581373030693,0.05720237901185122,0.06094003375348171,0.06431462975430591,0.06731617840704303,0.06991917852575025,0.07211465084605248,0.07388476204675995,0.07522290376748784,0.07611848768324433,0.07656803680773631,0.07656803680773631,0.07611848768324433,0.07522290376748784,0.07388476204675996,0.07211465084605248,0.06991917852575025,0.06731617840704303,0.06431462975430592,0.06094003375348171,0.05720237901185122,0.05313581373030693,0.04874891709851373,0.04408677229977932,0.03915185350905824,0.03400592341164634,0.02863445416487101,0.02313296372701731,0.01743533888979567,0.01174510563951802,0.005726732928160059,0.0005948839976204635,
	      0.0005672149744753256,0.005456885128538565,0.011195184692559,0.01662122715322138,0.02206151247961364,0.02731820608452054,0.03246106345682104,0.03739487520590695,0.04214023165323361,0.04663365908478271,0.05087955602936732,0.05482977592745247,0.05848294987824097,0.06180064073556969,0.06478025694793393,0.06739073436892289,0.0696306738573883,0.07147527342712572,0.07292578774925719,0.0739630598423674,0.07459196079387197,0.07479854105765875,0.07459196079387195,0.0739630598423674,0.07292578774925719,0.0714752734271257,0.0696306738573883,0.06739073436892291,0.06478025694793395,0.06180064073556966,0.05848294987824095,0.05482977592745249,0.05087955602936732,0.04663365908478271,0.04214023165323363,0.03739487520590695,0.03246106345682105,0.02731820608452054,0.02206151247961364,0.01662122715322139,0.01119518469255902,0.005456885128538564,0.0005672149744753256,
	      0.000540832882639264,0.005206865418115447,0.01068175047398483,0.01586389650829946,0.02106122596013223,0.02609120554585279,0.03101662696848615,0.03575285265126173,0.04031545988558687,0.0446500886187375,0.04875644191968971,0.05259480857089974,0.05615884202369557,0.05941816628598924,0.06236459322930067,0.0649749267485619,0.06724130029474086,0.0691467487175752,0.07068495827392642,0.07184476694780581,0.07262213270363554,0.07301150937108311,0.07301150937108311,0.07262213270363554,0.07184476694780581,0.07068495827392643,0.0691467487175752,0.06724130029474086,0.0649749267485619,0.06236459322930066,0.05941816628598926,0.05615884202369559,0.05259480857089974,0.04875644191968971,0.04465008861873752,0.04031545988558689,0.03575285265126172,0.03101662696848615,0.0260912055458528,0.02106122596013224,0.01586389650829949,0.01068175047398483,0.005206865418115448,0.000540832882639264,
	      0.0005167958656330745,0.00497254349993739,0.01020390434338543,0.01515593365268844,0.02012826456721103,0.02494330363170869,0.02966657063252143,0.0342137468348706,0.03860525927349107,0.04278538837871101,0.04675941320234074,0.05048518675690299,0.05396225152190157,0.05715685305151434,0.0600668644312374,0.06266475184240815,0.06494885941616645,0.0668968405341451,0.06850879553603016,0.06976700683991661,0.07067417370181789,0.07121684291001318,0.07140089915089437,0.07121684291001318,0.07067417370181789,0.06976700683991661,0.06850879553603016,0.06689684053414512,0.06494885941616646,0.06266475184240815,0.0600668644312374,0.05715685305151434,0.05396225152190157,0.05048518675690301,0.04675941320234074,0.04278538837871103,0.03860525927349107,0.03421374683487061,0.02966657063252143,0.02494330363170868,0.02012826456721103,0.01515593365268845,0.01020390434338542,0.004972543499937401,0.0005167958656330745,
	      0.0004938271604938267,0.004754685789078419,0.009756382704990518,0.01449522250419855,0.01925465702951085,0.02387001749850451,0.0284007921430443,0.03277148221801521,0.0369982356967542,0.04103290105630689,0.04487682754253925,0.04849486346395685,0.05188246784051915,0.05501256882062858,0.05787852950023086,0.06045933230212313,0.06274818815385574,0.06472922113844488,0.06639660207503394,0.06773917098341606,0.06875272509457971,0.06943062355009336,0.06977067573368116,0.06977067573368116,0.06943062355009336,0.06875272509457971,0.06773917098341607,0.06639660207503394,0.0647292211384449,0.06274818815385576,0.06045933230212314,0.05787852950023088,0.0550125688206286,0.05188246784051914,0.04849486346395689,0.04487682754253926,0.04103290105630689,0.0369982356967542,0.03277148221801521,0.02840079214304431,0.02387001749850452,0.01925465702951084,0.01449522250419854,0.009756382704990534,0.00475468578907844,0.0004938271604938267,
	      0.0004728132387706851,0.004549918906835608,0.009338558447231077,0.01387576851196736,0.01843740846696191,0.02286313019525714,0.02721446543111642,0.03141623548624119,0.03548849939498459,0.03938201221537949,0.04310268433141279,0.04661339939561216,0.04991446503159971,0.05297610679954393,0.05579669470134412,0.05835177394850422,0.06063968870658444,0.06264033773601961,0.0643531782600336,0.0657619481765247,0.06686796077881772,0.06765847527591694,0.0681371765917576,0.06829459994316581,0.0681371765917576,0.06765847527591695,0.06686796077881772,0.0657619481765247,0.0643531782600336,0.06264033773601962,0.06063968870658444,0.05835177394850424,0.05579669470134413,0.05297610679954394,0.04991446503159971,0.04661339939561215,0.04310268433141281,0.03938201221537949,0.03548849939498461,0.0314162354862412,0.02721446543111642,0.02286313019525716,0.0184374084669619,0.01387576851196739,0.009338558447231086,0.00454991890683564,0.0004728132387706851,
	      0.0004526935264825709,0.004358936311639493,0.008946143006177722,0.01329595069099517,0.01767004442559142,0.02191908600409841,0.02609935151838274,0.03014308514775314,0.0340666639738143,0.03782718116347588,0.04142722290785926,0.04483552154181742,0.04804899144742766,0.05104344236766516,0.05381355737508891,0.05634032248075507,0.05861794474565336,0.0606317134993067,0.06237636560294586,0.06384107504286643,0.06502172071520999,0.06591115957503207,0.0665067845868649,0.06680504234309634,0.06680504234309634,0.06650678458686492,0.06591115957503207,0.06502172071520999,0.06384107504286643,0.06237636560294586,0.0606317134993067,0.05861794474565334,0.05634032248075508,0.0538135573750889,0.05104344236766516,0.04804899144742768,0.04483552154181743,0.04142722290785926,0.03782718116347587,0.03406666397381428,0.03014308514775317,0.02609935151838275,0.02191908600409842,0.01767004442559144,0.01329595069099517,0.008946143006177739,0.004358936311639507,0.0004526935264825709,
	      0.0004342162396873638,0.004178956876062836,0.008578703258605087,0.01275087522924873,0.01695020950305129,0.02103117305110113,0.02505156862017066,0.02894395613926466,0.03272789002857068,0.03635957496138607,0.03984535647647143,0.04315236232703318,0.04628152867824215,0.04920647366797049,0.05192604724781408,0.05441847236255171,0.05668223224153558,0.05869924887386987,0.06046865526703729,0.06197559292242948,0.06322050412135194,0.06419146399372862,0.06489068028737092,0.06530895771519693,0.06545059982049571,0.06530895771519694,0.06489068028737092,0.06419146399372862,0.06322050412135195,0.06197559292242948,0.0604686552670373,0.05869924887386988,0.05668223224153558,0.05441847236255171,0.05192604724781407,0.04920647366797049,0.04628152867824215,0.04315236232703321,0.03984535647647144,0.03635957496138605,0.03272789002857069,0.02894395613926465,0.02505156862017066,0.02103117305110114,0.01695020950305132,0.01275087522924875,0.008578703258605082,0.004178956876062839,0.0004342162396873638,
	      0.0004164931278633898,0.004010605491385264,0.008232712739549324,0.01223929186020785,0.01627259795966257,0.02019654926393445,0.02406434365584402,0.02781483820529734,0.03146439660419052,0.03497457327189035,0.03834887561628814,0.04155935855194848,0.04460410666255692,0.04746143902379512,0.05012715304034253,0.05258404852313164,0.05482722859064214,0.05684314192350467,0.05862710954665857,0.06016882083339292,0.06146437491563957,0.06250649868053845,0.06329242015374517,0.06381780458611555,0.06408121717187507,0.06408121717187507,0.06381780458611555,0.06329242015374517,0.06250649868053845,0.06146437491563959,0.06016882083339293,0.05862710954665857,0.05684314192350468,0.05482722859064216,0.05258404852313164,0.05012715304034251,0.04746143902379513,0.04460410666255692,0.04155935855194848,0.03834887561628816,0.03497457327189036,0.03146439660419052,0.02781483820529735,0.02406434365584401,0.02019654926393446,0.01627259795966258,0.01223929186020788,0.008232712739549331,0.004010605491385286,0.0004164931278633898,
	      0.0002519526329050135,0.002426825577091331,0.004985682800997247,0.007422108895215612,0.009886791733980157,0.01230107564932004,0.01470112902792591,0.01705337266361986,0.01937142575637869,0.02163508760827157,0.02384983935834243,0.0260013856529256,0.02809148349447062,0.03010909164371389,0.03205406836747507,0.03391746772918589,0.03569815436810179,0.03738870431340301,0.03898749674668905,0.04048832316267238,0.04188938934228031,0.04318553037633184,0.04437498297779119,0.04545352684357432,0.04641956917622356,0.04726977666978605,0.04800282424610153,0.04861623221720066,0.04910901048034229,0.04947951406104645,0.0497271321125967,0.04985104431403992,0.04985104431403992,0.04972713211259669,0.04947951406104646,0.04910901048034228,0.04861623221720066,0.04800282424610155,0.04726977666978605,0.04641956917622356,0.04545352684357432,0.0443749829777912,0.04318553037633184,0.04188938934228031,0.04048832316267238,0.03898749674668904,0.037388704313403,0.03569815436810179,0.03391746772918588,0.03205406836747508,0.0301090916437139,0.02809148349447064,0.02600138565292559,0.02384983935834244,0.02163508760827157,0.01937142575637868,0.01705337266361986,0.01470112902792591,0.01230107564932006,0.009886791733980153,0.007422108895215627,0.004985682800997249,0.002426825577091351,0.0002519526329050135,
	      6.200012400024783e-05,0.0005973762149118416,0.001228404325585456,0.001831569442962184,0.002445140720692403,0.003050829883647109,0.003658708757911623,0.004261545074217487,0.004863826631094891,0.005461565021654319,0.006057196343139136,0.006648123879881041,0.007235805607387576,0.007818371069374283,0.008396737678707451,0.008969463595819406,0.009537137738428422,0.01009859402781798,0.01065420834930766,0.01120300402074574,0.01174521229621522,0.01227999356617884,0.01280747767569352,0.01332692866341296,0.01383840371367647,0.01434124826411841,0.01483546677824516,0.01532047079201182,0.01579622637151619,0.01626220035154729,0.01671833100026016,0.01716413266809309,0.0175995238714267,0.01802406077249,0.01843764837928557,0.01883988042991853,0.01923065336073684,0.01960959530712628,0.01997659810049917,0.02033132186943132,0.02067365707086236,0.02100329399785076,0.0213201243926148,0.02162386731644146,0.02191441800514906,0.02219152322011034,0.02245508353487128,0.02270487259356926,0.0229407978520228,0.02316265921267624,0.02337037230693285,0.02356376282006356,0.02374275563759934,0.02390720186767724,0.02405703654135776,0.02419213591961974,0.02431244590425881,0.02441786770948752,0.02450835868263871,0.02458384484722067,0.02464429543223335,0.02468966117132483,0.02471992348106206,0.02473505774318532,0.02473505774318531,0.02471992348106206,0.02468966117132484,0.02464429543223335,0.02458384484722067,0.02450835868263871,0.02441786770948752,0.02431244590425881,0.02419213591961973,0.02405703654135777,0.02390720186767725,0.02374275563759934,0.02356376282006355,0.02337037230693285,0.02316265921267624,0.0229407978520228,0.02270487259356926,0.02245508353487128,0.02219152322011034,0.02191441800514905,0.02162386731644147,0.0213201243926148,0.02100329399785077,0.02067365707086237,0.02033132186943133,0.01997659810049917,0.01960959530712629,0.01923065336073685,0.01883988042991853,0.01843764837928556,0.01802406077249001,0.01759952387142671,0.01716413266809309,0.01671833100026017,0.01626220035154729,0.0157962263715162,0.01532047079201182,0.01483546677824517,0.01434124826411841,0.01383840371367648,0.01332692866341296,0.01280747767569352,0.01227999356617884,0.01174521229621523,0.01120300402074574,0.01065420834930767,0.01009859402781799,0.009537137738428422,0.008969463595819408,0.008396737678707456,0.007818371069374287,0.007235805607387577,0.006648123879881052,0.006057196343139142,0.005461565021654323,0.004863826631094894,0.004261545074217493,0.003658708757911625,0.003050829883647109,0.002445140720692403,0.001831569442962192,0.00122840432558546,0.0005973762149118454,6.200012400024783e-05,
	      1.537870049980768e-05,0.0001481863885321792,0.0003047906606122091,0.0004546218333176311,0.0006072451420738039,0.000758190333667881,0.0009100285731154968,0.001061031268255706,0.00121238405584298,0.001363162175987088,0.001514041635542977,0.00166444190861389,0.001814795443749086,0.001964700691793535,0.002114454869997663,0.002263761371099011,0.002412834584638254,0.002561444723686011,0.002709751797595784,0.002857571213050027,0.003005025379005977,0.003151961715725076,0.003298475571746596,0.003444437892760579,0.003589923922687423,0.003734822421001835,0.003879193304332365,0.004022939162895638,0.004166107977100523,0.004308613306836135,0.004450493671026209,0.004591671492205189,0.004732177677088058,0.004871941925800522,0.00501098894333844,0.005149254492992578,0.005286758173305948,0.005423440865349684,0.0055593179252744,0.005694334605665305,0.005828502711624328,0.005961771270897006,0.006094149097736766,0.006225588513294813,0.006356095800134502,0.006485626179865219,0.006614183783637206,0.006741726410240587,0.006868256357367531,0.006993733732978106,0.00711815926948246,0.007241495160284487,0.007363740800527748,0.007484860281145577,0.007604851855328713,0.007723681352829043,0.007841346053341023,0.00795781339072186,0.008073079817392142,0.008187114256460076,0.008299912460749601,0.008411444744305852,0.008521706272455841,0.008630668665725613,0.008738326600872674,0.008844652932122275,0.008949641935380783,0.009053267635674815,0.00915552398618168,0.009256386128238401,0.00935584776215184,0.009453885098259122,0.009550491646699823,0.009645644645658017,0.009739337471579344,0.009831548354639908,0.009922270588612357,0.01001148336438384,0.0100991799392781,0.01018534043757266,0.01026995812212519,0.0103530140267207,0.01043450145796514,0.01051440233825983,0.01059271005280798,0.01066940739434503,0.01074448785850083,0.01081793509234233,0.01088974273103286,0.01095989526196328,0.01102838648647086,0.01109520172001142,0.01116033495449131,0.01122377232270729,0.01128550802947646,0.01134552901556086,0.0114038297191431,0.01146039788076094,0.0115152281906747,0.01156830918205265,0.01161963581432855,0.01166919740707611,0.01171698920449162,0.01176300130714061,0.01180722925815999,0.01184966393440995,0.01189030119081856,0.01192913267647676,0.01196615456969921,0.01200135928830436,0.01203474334439715,0.0120662999215175,0.0120960258748269,0.01212391515102317,0.01214996495650083,0.01217416999894672,0.01219652784311505,0.01221703395586728,0.01223568626642881,0.01225248099934087,0.0122674164534251,0.01228048960969914,0.01229169914074282,0.01230104278311453,0.01230851958637052,0.01231412804192391,0.01231786757859581,0.01231973744220419,0.01231973744220419,0.01231786757859581,0.0123141280419239,0.01230851958637052,0.01230104278311453,0.01229169914074282,0.01228048960969914,0.0122674164534251,0.01225248099934087,0.01223568626642881,0.01221703395586728,0.01219652784311505,0.01217416999894671,0.01214996495650083,0.01212391515102317,0.0120960258748269,0.01206629992151749,0.01203474334439715,0.01200135928830436,0.0119661545696992,0.01192913267647676,0.01189030119081856,0.01184966393440995,0.01180722925815999,0.01176300130714061,0.01171698920449162,0.01166919740707611,0.01161963581432855,0.01156830918205266,0.0115152281906747,0.01146039788076095,0.0114038297191431,0.01134552901556086,0.01128550802947646,0.01122377232270729,0.01116033495449131,0.01109520172001142,0.01102838648647086,0.01095989526196328,0.01088974273103286,0.01081793509234233,0.01074448785850084,0.01066940739434503,0.01059271005280798,0.01051440233825983,0.01043450145796514,0.0103530140267207,0.01026995812212519,0.01018534043757266,0.0100991799392781,0.01001148336438385,0.009922270588612355,0.00983154835463991,0.009739337471579346,0.00964564464565802,0.009550491646699829,0.009453885098259129,0.009355847762151837,0.009256386128238399,0.009155523986181678,0.009053267635674813,0.008949641935380786,0.008844652932122278,0.008738326600872676,0.008630668665725611,0.008521706272455843,0.008411444744305854,0.008299912460749606,0.008187114256460081,0.00807307981739214,0.00795781339072186,0.007841346053341023,0.007723681352829043,0.007604851855328716,0.007484860281145576,0.007363740800527745,0.007241495160284487,0.007118159269482463,0.006993733732978105,0.006868256357367535,0.00674172641024059,0.006614183783637209,0.006485626179865222,0.006356095800134502,0.006225588513294813,0.006094149097736764,0.005961771270897004,0.005828502711624328,0.005694334605665304,0.005559317925274401,0.005423440865349687,0.005286758173305949,0.00514925449299258,0.00501098894333844,0.004871941925800526,0.00473217767708806,0.004591671492205188,0.004450493671026206,0.004308613306836135,0.004166107977100525,0.004022939162895638,0.003879193304332367,0.003734822421001837,0.003589923922687428,0.003444437892760581,0.003298475571746598,0.003151961715725079,0.003005025379005982,0.002857571213050024,0.002709751797595783,0.002561444723686011,0.002412834584638254,0.002263761371099012,0.002114454869997664,0.001964700691793538,0.001814795443749087,0.001664441908613892,0.00151404163554298,0.001363162175987091,0.001212384055842984,0.00106103126825571,0.000910028573115496,0.0007581903336678806,0.0006072451420738037,0.0004546218333176309,0.0003047906606122096,0.0001481863885321808,1.537870049980768e-05,
	    };

	  double knots_data_array[] =
	    {
	      0,
	      1,-1,
	      1,6.123031769111886e-17,-1,
	      1,0.5000000000000001,-0.4999999999999998,-1,
	      1,0.7071067811865476,6.123031769111886e-17,-0.7071067811865475,-1,
	      1,0.8090169943749475,0.3090169943749475,-0.3090169943749473,-0.8090169943749473,-1,
	      1,0.8660254037844387,0.5000000000000001,6.123031769111886e-17,-0.4999999999999998,-0.8660254037844387,-1,
	      1,0.9009688679024191,0.6234898018587336,0.2225209339563144,-0.2225209339563143,-0.6234898018587335,-0.900968867902419,-1,
	      1,0.9238795325112867,0.7071067811865476,0.3826834323650898,6.123031769111886e-17,-0.3826834323650897,-0.7071067811865475,-0.9238795325112867,-1,
	      1,0.9396926207859084,0.766044443118978,0.5000000000000001,0.1736481776669304,-0.1736481776669303,-0.4999999999999998,-0.7660444431189779,-0.9396926207859083,-1,
	      1,0.9510565162951535,0.8090169943749475,0.5877852522924731,0.3090169943749475,6.123031769111886e-17,-0.3090169943749473,-0.587785252292473,-0.8090169943749473,-0.9510565162951535,-1,
	      1,0.9594929736144974,0.8412535328311812,0.6548607339452851,0.4154150130018864,0.1423148382732851,-0.142314838273285,-0.4154150130018863,-0.654860733945285,-0.8412535328311811,-0.9594929736144974,-1,
	      1,0.9659258262890683,0.8660254037844387,0.7071067811865476,0.5000000000000001,0.2588190451025207,6.123031769111886e-17,-0.2588190451025206,-0.4999999999999998,-0.7071067811865475,-0.8660254037844387,-0.9659258262890682,-1,
	      1,0.970941817426052,0.8854560256532099,0.7485107481711011,0.5680647467311559,0.3546048870425356,0.120536680255323,-0.1205366802553229,-0.3546048870425355,-0.5680647467311557,-0.7485107481711012,-0.8854560256532098,-0.970941817426052,-1,
	      1,0.9749279121818236,0.9009688679024191,0.7818314824680298,0.6234898018587336,0.4338837391175582,0.2225209339563144,6.123031769111886e-17,-0.2225209339563143,-0.4338837391175581,-0.6234898018587335,-0.7818314824680298,-0.900968867902419,-0.9749279121818236,-1,
	      1,0.9781476007338057,0.9135454576426009,0.8090169943749475,0.6691306063588582,0.5000000000000001,0.3090169943749475,0.1045284632676535,-0.1045284632676533,-0.3090169943749473,-0.4999999999999998,-0.6691306063588582,-0.8090169943749473,-0.9135454576426008,-0.9781476007338057,-1,
	      1,0.9807852804032304,0.9238795325112867,0.8314696123025452,0.7071067811865476,0.5555702330196023,0.3826834323650898,0.1950903220161283,6.123031769111886e-17,-0.1950903220161282,-0.3826834323650897,-0.555570233019602,-0.7071067811865475,-0.8314696123025453,-0.9238795325112867,-0.9807852804032304,-1,
	      1,0.9829730996839018,0.9324722294043558,0.8502171357296142,0.7390089172206591,0.6026346363792564,0.4457383557765383,0.273662990072083,0.09226835946330202,-0.09226835946330189,-0.2736629900720829,-0.4457383557765382,-0.6026346363792563,-0.739008917220659,-0.850217135729614,-0.9324722294043558,-0.9829730996839018,-1,
	      1,0.984807753012208,0.9396926207859084,0.8660254037844387,0.766044443118978,0.6427876096865394,0.5000000000000001,0.3420201433256688,0.1736481776669304,6.123031769111886e-17,-0.1736481776669303,-0.3420201433256687,-0.4999999999999998,-0.6427876096865394,-0.7660444431189779,-0.8660254037844387,-0.9396926207859083,-0.984807753012208,-1,
	      1,0.9863613034027223,0.9458172417006346,0.8794737512064891,0.7891405093963936,0.6772815716257411,0.5469481581224269,0.4016954246529695,0.2454854871407992,0.08257934547233239,-0.08257934547233227,-0.2454854871407991,-0.4016954246529694,-0.5469481581224267,-0.6772815716257409,-0.7891405093963935,-0.879473751206489,-0.9458172417006346,-0.9863613034027223,-1,
	      1,0.9876883405951378,0.9510565162951535,0.8910065241883679,0.8090169943749475,0.7071067811865476,0.5877852522924731,0.4539904997395468,0.3090169943749475,0.1564344650402309,6.123031769111886e-17,-0.1564344650402308,-0.3090169943749473,-0.4539904997395467,-0.587785252292473,-0.7071067811865475,-0.8090169943749473,-0.8910065241883678,-0.9510565162951535,-0.9876883405951377,-1,
	      1,0.9888308262251285,0.9555728057861407,0.9009688679024191,0.8262387743159949,0.7330518718298263,0.6234898018587336,0.5000000000000001,0.365341024366395,0.2225209339563144,0.07473009358642439,-0.07473009358642427,-0.2225209339563143,-0.3653410243663949,-0.4999999999999998,-0.6234898018587335,-0.7330518718298263,-0.8262387743159947,-0.900968867902419,-0.9555728057861408,-0.9888308262251285,-1,
	      1,0.9898214418809327,0.9594929736144974,0.9096319953545184,0.8412535328311812,0.7557495743542583,0.6548607339452851,0.5406408174555977,0.4154150130018864,0.2817325568414298,0.1423148382732851,6.123031769111886e-17,-0.142314838273285,-0.2817325568414297,-0.4154150130018863,-0.5406408174555977,-0.654860733945285,-0.7557495743542582,-0.8412535328311811,-0.9096319953545182,-0.9594929736144974,-0.9898214418809327,-1,
	      1,0.9906859460363308,0.9629172873477992,0.917211301505453,0.8544194045464886,0.7757112907044198,0.6825531432186541,0.5766803221148672,0.4600650377311522,0.3348796121709863,0.2034560130526337,0.068242413364671,-0.06824241336467088,-0.2034560130526336,-0.3348796121709862,-0.4600650377311521,-0.5766803221148671,-0.6825531432186539,-0.7757112907044197,-0.8544194045464883,-0.917211301505453,-0.9629172873477992,-0.9906859460363308,-1,
	      1,0.9914448613738104,0.9659258262890683,0.9238795325112867,0.8660254037844387,0.7933533402912352,0.7071067811865476,0.6087614290087207,0.5000000000000001,0.3826834323650898,0.2588190451025207,0.1305261922200517,6.123031769111886e-17,-0.1305261922200516,-0.2588190451025206,-0.3826834323650897,-0.4999999999999998,-0.6087614290087207,-0.7071067811865475,-0.7933533402912351,-0.8660254037844387,-0.9238795325112867,-0.9659258262890682,-0.9914448613738104,-1,
	      1,0.9921147013144779,0.9685831611286311,0.9297764858882515,0.8763066800438636,0.8090169943749475,0.7289686274214116,0.6374239897486897,0.5358267949789965,0.4257792915650727,0.3090169943749475,0.1873813145857247,0.06279051952931353,-0.0627905195293134,-0.1873813145857246,-0.3090169943749473,-0.4257792915650727,-0.5358267949789964,-0.6374239897486897,-0.7289686274214113,-0.8090169943749473,-0.8763066800438636,-0.9297764858882513,-0.9685831611286311,-0.9921147013144778,-1,
	      1,0.992708874098054,0.970941817426052,0.9350162426854148,0.8854560256532099,0.8229838658936564,0.7485107481711011,0.6631226582407953,0.5680647467311559,0.4647231720437686,0.3546048870425356,0.2393156642875578,0.120536680255323,6.123031769111886e-17,-0.1205366802553229,-0.2393156642875577,-0.3546048870425355,-0.4647231720437685,-0.5680647467311557,-0.663122658240795,-0.7485107481711012,-0.8229838658936564,-0.8854560256532098,-0.9350162426854147,-0.970941817426052,-0.992708874098054,-1,
	      1,0.993238357741943,0.9730448705798238,0.9396926207859084,0.8936326403234123,0.8354878114129365,0.766044443118978,0.6862416378687336,0.5971585917027862,0.5000000000000001,0.3960797660391569,0.2868032327110903,0.1736481776669304,0.0581448289104759,-0.05814482891047577,-0.1736481776669303,-0.2868032327110902,-0.3960797660391568,-0.4999999999999998,-0.597158591702786,-0.6862416378687335,-0.7660444431189779,-0.8354878114129363,-0.8936326403234122,-0.9396926207859083,-0.9730448705798238,-0.993238357741943,-1,
	      1,0.9937122098932426,0.9749279121818236,0.9438833303083676,0.9009688679024191,0.8467241992282841,0.7818314824680298,0.7071067811865476,0.6234898018587336,0.5320320765153366,0.4338837391175582,0.3302790619551672,0.2225209339563144,0.1119644761033079,6.123031769111886e-17,-0.1119644761033078,-0.2225209339563143,-0.330279061955167,-0.4338837391175581,-0.5320320765153365,-0.6234898018587335,-0.7071067811865475,-0.7818314824680298,-0.8467241992282841,-0.900968867902419,-0.9438833303083676,-0.9749279121818236,-0.9937122098932426,-1,
	      1,0.9941379571543596,0.9766205557100867,0.9476531711828025,0.907575419670957,0.8568571761675893,0.7960930657056438,0.7259954919231308,0.6473862847818277,0.5611870653623824,0.4684084406997902,0.3701381553399143,0.2675283385292208,0.1617819965527648,0.05413890858541761,-0.05413890858541748,-0.1617819965527647,-0.2675283385292206,-0.3701381553399142,-0.46840844069979,-0.5611870653623823,-0.6473862847818276,-0.7259954919231308,-0.7960930657056438,-0.8568571761675893,-0.9075754196709569,-0.9476531711828023,-0.9766205557100867,-0.9941379571543596,-1,
	      1,0.9945218953682733,0.9781476007338057,0.9510565162951535,0.9135454576426009,0.8660254037844387,0.8090169943749475,0.7431448254773942,0.6691306063588582,0.5877852522924731,0.5000000000000001,0.4067366430758002,0.3090169943749475,0.2079116908177595,0.1045284632676535,6.123031769111886e-17,-0.1045284632676533,-0.2079116908177593,-0.3090169943749473,-0.4067366430758,-0.4999999999999998,-0.587785252292473,-0.6691306063588582,-0.743144825477394,-0.8090169943749473,-0.8660254037844387,-0.9135454576426008,-0.9510565162951535,-0.9781476007338057,-0.9945218953682733,-1,
	      1,0.9948693233918952,0.9795299412524945,0.9541392564000488,0.9189578116202306,0.8743466161445821,0.8207634412072763,0.7587581226927909,0.6889669190756866,0.6121059825476629,0.5289640103269624,0.4403941515576343,0.3473052528448203,0.2506525322587205,0.1514277775045767,0.05064916883871277,-0.05064916883871264,-0.1514277775045766,-0.2506525322587204,-0.3473052528448202,-0.4403941515576344,-0.5289640103269625,-0.6121059825476629,-0.6889669190756866,-0.7587581226927909,-0.8207634412072763,-0.8743466161445821,-0.9189578116202306,-0.9541392564000488,-0.9795299412524945,-0.9948693233918952,-1,
	      1,0.9951847266721969,0.9807852804032304,0.9569403357322088,0.9238795325112867,0.881921264348355,0.8314696123025452,0.773010453362737,0.7071067811865476,0.6343932841636455,0.5555702330196023,0.4713967368259978,0.3826834323650898,0.2902846772544623,0.1950903220161283,0.09801714032956077,6.123031769111886e-17,-0.09801714032956065,-0.1950903220161282,-0.2902846772544622,-0.3826834323650897,-0.4713967368259977,-0.555570233019602,-0.6343932841636454,-0.7071067811865475,-0.773010453362737,-0.8314696123025453,-0.8819212643483549,-0.9238795325112867,-0.9569403357322088,-0.9807852804032304,-0.9951847266721968,-1,
	      1,0.9954719225730846,0.9819286972627067,0.9594929736144974,0.9283679330160726,0.8888354486549235,0.8412535328311812,0.7860530947427875,0.7237340381050702,0.6548607339452851,0.5800569095711982,0.5000000000000001,0.4154150130018864,0.3270679633174218,0.2357589355094273,0.1423148382732851,0.0475819158237424,-0.04758191582374228,-0.142314838273285,-0.2357589355094272,-0.3270679633174217,-0.4154150130018863,-0.4999999999999998,-0.580056909571198,-0.654860733945285,-0.7237340381050702,-0.7860530947427873,-0.8412535328311811,-0.8888354486549234,-0.9283679330160726,-0.9594929736144974,-0.9819286972627066,-0.9954719225730846,-1,
	      1,0.9957341762950345,0.9829730996839018,0.961825643172819,0.9324722294043558,0.8951632913550623,0.8502171357296142,0.7980172272802396,0.7390089172206591,0.6736956436465572,0.6026346363792564,0.5264321628773558,0.4457383557765383,0.361241666187153,0.273662990072083,0.1837495178165703,0.09226835946330202,6.123031769111886e-17,-0.09226835946330189,-0.1837495178165702,-0.2736629900720829,-0.3612416661871529,-0.4457383557765382,-0.5264321628773559,-0.6026346363792563,-0.6736956436465572,-0.739008917220659,-0.7980172272802395,-0.850217135729614,-0.8951632913550622,-0.9324722294043558,-0.961825643172819,-0.9829730996839018,-0.9957341762950345,-1,
	      1,0.9959742939952391,0.9839295885986297,0.9639628606958532,0.9362348706397372,0.9009688679024191,0.8584487936018661,0.8090169943749475,0.753071466003611,0.6910626489868646,0.6234898018587336,0.5508969814521025,0.4738686624729987,0.3930250316539237,0.3090169943749475,0.2225209339563144,0.1342332658176555,0.04486483035051499,-0.04486483035051486,-0.1342332658176554,-0.2225209339563143,-0.3090169943749473,-0.3930250316539236,-0.4738686624729986,-0.5508969814521024,-0.6234898018587335,-0.6910626489868646,-0.7530714660036109,-0.8090169943749473,-0.8584487936018661,-0.900968867902419,-0.9362348706397372,-0.9639628606958532,-0.9839295885986297,-0.9959742939952391,-1,
	      1,0.9961946980917455,0.984807753012208,0.9659258262890683,0.9396926207859084,0.9063077870366499,0.8660254037844387,0.8191520442889918,0.766044443118978,0.7071067811865476,0.6427876096865394,0.5735764363510462,0.5000000000000001,0.4226182617406994,0.3420201433256688,0.2588190451025207,0.1736481776669304,0.08715574274765814,6.123031769111886e-17,-0.08715574274765801,-0.1736481776669303,-0.2588190451025206,-0.3420201433256687,-0.4226182617406993,-0.4999999999999998,-0.5735764363510462,-0.6427876096865394,-0.7071067811865475,-0.7660444431189779,-0.8191520442889916,-0.8660254037844387,-0.9063077870366499,-0.9396926207859083,-0.9659258262890682,-0.984807753012208,-0.9961946980917455,-1,
	      1,0.9963974885425265,0.9856159103477085,0.9677329469334989,0.9428774454610842,0.9112284903881357,0.8730141131611882,0.8285096492438422,0.7780357543184395,0.7219560939545245,0.6606747233900815,0.5946331763042867,0.5243072835572317,0.4502037448176733,0.3728564777803086,0.2928227712765504,0.2106792699957264,0.1270178197468789,0.04244120319614846,-0.04244120319614834,-0.1270178197468788,-0.2106792699957263,-0.2928227712765503,-0.3728564777803085,-0.4502037448176734,-0.5243072835572316,-0.5946331763042866,-0.6606747233900813,-0.7219560939545244,-0.7780357543184393,-0.8285096492438421,-0.8730141131611882,-0.9112284903881356,-0.9428774454610842,-0.9677329469334988,-0.9856159103477085,-0.9963974885425265,-1,
	      1,0.9965844930066698,0.9863613034027223,0.9694002659393304,0.9458172417006346,0.9157733266550574,0.8794737512064891,0.8371664782625287,0.7891405093963936,0.7357239106731317,0.6772815716257411,0.6142127126896678,0.5469481581224269,0.4759473930370736,0.4016954246529695,0.3246994692046836,0.2454854871407992,0.164594590280734,0.08257934547233239,6.123031769111886e-17,-0.08257934547233227,-0.1645945902807338,-0.2454854871407991,-0.3246994692046835,-0.4016954246529694,-0.4759473930370736,-0.5469481581224267,-0.6142127126896678,-0.6772815716257409,-0.7357239106731316,-0.7891405093963935,-0.8371664782625285,-0.879473751206489,-0.9157733266550575,-0.9458172417006346,-0.9694002659393305,-0.9863613034027223,-0.9965844930066698,-1,
	      1,0.99675730813421,0.9870502626379128,0.970941817426052,0.9485364419471455,0.9199794436588242,0.8854560256532099,0.8451900855437947,0.7994427634035012,0.7485107481711011,0.6927243535095994,0.6324453755953773,0.5680647467311559,0.5000000000000001,0.4286925614030542,0.3546048870425356,0.2782174639164527,0.2000256937760445,0.120536680255323,0.04026594010941524,-0.04026594010941512,-0.1205366802553229,-0.2000256937760443,-0.2782174639164526,-0.3546048870425355,-0.4286925614030543,-0.4999999999999998,-0.5680647467311557,-0.6324453755953772,-0.6927243535095994,-0.7485107481711012,-0.799442763403501,-0.8451900855437946,-0.8854560256532098,-0.9199794436588242,-0.9485364419471455,-0.970941817426052,-0.9870502626379128,-0.9967573081342099,-1,
	      1,0.996917333733128,0.9876883405951378,0.9723699203976766,0.9510565162951535,0.9238795325112867,0.8910065241883679,0.8526401643540922,0.8090169943749475,0.7604059656000309,0.7071067811865476,0.6494480483301837,0.5877852522924731,0.5224985647159489,0.4539904997395468,0.3826834323650898,0.3090169943749475,0.2334453638559055,0.1564344650402309,0.078459095727845,6.123031769111886e-17,-0.07845909572784487,-0.1564344650402308,-0.2334453638559053,-0.3090169943749473,-0.3826834323650897,-0.4539904997395467,-0.5224985647159488,-0.587785252292473,-0.6494480483301835,-0.7071067811865475,-0.7604059656000309,-0.8090169943749473,-0.8526401643540922,-0.8910065241883678,-0.9238795325112867,-0.9510565162951535,-0.9723699203976766,-0.9876883405951377,-0.996917333733128,-1,
	      1,0.9970658011837404,0.9882804237803485,0.9736954238777791,0.9533963920549305,0.9275024511020947,0.8961655569610556,0.8595696069872012,0.8179293607667176,0.771489179821943,0.720521593600787,0.6653257001655654,0.6062254109666381,0.5435675500012211,0.477719818512263,0.40906863717134,0.3380168784085027,0.2649815021966617,0.1903911091646684,0.1146834253984005,0.03830273369003549,-0.03830273369003537,-0.1146834253984004,-0.1903911091646683,-0.2649815021966616,-0.3380168784085026,-0.4090686371713399,-0.4777198185122627,-0.543567550001221,-0.6062254109666381,-0.6653257001655652,-0.7205215936007869,-0.771489179821943,-0.8179293607667176,-0.8595696069872012,-0.8961655569610555,-0.9275024511020946,-0.9533963920549305,-0.973695423877779,-0.9882804237803485,-0.9970658011837404,-1,
	      1,0.9972037971811801,0.9888308262251285,0.9749279121818236,0.9555728057861407,0.9308737486442042,0.9009688679024191,0.8660254037844387,0.8262387743159949,0.7818314824680298,0.7330518718298263,0.6801727377709194,0.6234898018587336,0.5633200580636221,0.5000000000000001,0.4338837391175582,0.365341024366395,0.2947551744109043,0.2225209339563144,0.1490422661761744,0.07473009358642439,6.123031769111886e-17,-0.07473009358642427,-0.1490422661761743,-0.2225209339563143,-0.2947551744109042,-0.3653410243663949,-0.4338837391175581,-0.4999999999999998,-0.5633200580636221,-0.6234898018587335,-0.6801727377709192,-0.7330518718298263,-0.7818314824680298,-0.8262387743159947,-0.8660254037844387,-0.900968867902419,-0.9308737486442041,-0.9555728057861408,-0.9749279121818236,-0.9888308262251285,-0.9972037971811801,-1,
	      1,0.9973322836635516,0.9893433680751103,0.9760758775559272,0.957600599908406,0.934016108732548,0.9054482374931466,0.8720494081438076,0.8339978178898779,0.7914964884292541,0.7447721827437819,0.694074195220634,0.6396730215588913,0.5818589155579529,0.5209403404879303,0.4572423233046386,0.391104720490156,0.3228804047714463,0.2529333823916807,0.1816368509794365,0.1093712083778745,0.03652202305765885,-0.03652202305765873,-0.1093712083778744,-0.1816368509794364,-0.2529333823916806,-0.3228804047714462,-0.3911047204901559,-0.4572423233046385,-0.5209403404879301,-0.5818589155579527,-0.6396730215588911,-0.694074195220634,-0.7447721827437819,-0.791496488429254,-0.8339978178898778,-0.8720494081438077,-0.9054482374931466,-0.9340161087325479,-0.9576005999084058,-0.9760758775559271,-0.9893433680751103,-0.9973322836635516,-1,
	      1,0.9974521146102535,0.9898214418809327,0.9771468659711595,0.9594929736144974,0.9369497249997617,0.9096319953545184,0.8776789895672557,0.8412535328311812,0.8005412409243604,0.7557495743542583,0.7071067811865476,0.6548607339452851,0.599277666511347,0.5406408174555977,0.4792489867200568,0.4154150130018864,0.3494641795990984,0.2817325568414298,0.2125652895529768,0.1423148382732851,0.07133918319923235,6.123031769111886e-17,-0.07133918319923224,-0.142314838273285,-0.2125652895529767,-0.2817325568414297,-0.3494641795990983,-0.4154150130018863,-0.4792489867200569,-0.5406408174555977,-0.599277666511347,-0.654860733945285,-0.7071067811865475,-0.7557495743542582,-0.8005412409243603,-0.8412535328311811,-0.8776789895672555,-0.9096319953545182,-0.9369497249997618,-0.9594929736144974,-0.9771468659711595,-0.9898214418809327,-0.9974521146102535,-1,
	      1,0.9975640502598242,0.9902680687415704,0.9781476007338057,0.9612616959383189,0.9396926207859084,0.9135454576426009,0.882947592858927,0.848048096156426,0.8090169943749475,0.766044443118978,0.7193398003386512,0.6691306063588582,0.6156614753256583,0.5591929034707468,0.5000000000000001,0.4383711467890775,0.3746065934159122,0.3090169943749475,0.2419218955996677,0.1736481776669304,0.1045284632676535,0.03489949670250108,-0.03489949670250096,-0.1045284632676533,-0.1736481776669303,-0.2419218955996676,-0.3090169943749473,-0.3746065934159121,-0.4383711467890775,-0.4999999999999998,-0.5591929034707467,-0.6156614753256583,-0.6691306063588582,-0.719339800338651,-0.7660444431189779,-0.8090169943749473,-0.848048096156426,-0.882947592858927,-0.9135454576426008,-0.9396926207859083,-0.9612616959383189,-0.9781476007338057,-0.9902680687415703,-0.9975640502598242,-1,
	      1,0.9976687691905392,0.9906859460363308,0.9790840876823229,0.9629172873477992,0.9422609221188205,0.917211301505453,0.8878852184023752,0.8544194045464886,0.8169698930104421,0.7757112907044198,0.7308359642781241,0.6825531432186541,0.6310879443260529,0.5766803221148672,0.5195839500354336,0.4600650377311522,0.3984010898462416,0.3348796121709863,0.2697967711570244,0.2034560130526337,0.1361666490962466,0.068242413364671,6.123031769111886e-17,-0.06824241336467088,-0.1361666490962465,-0.2034560130526336,-0.2697967711570243,-0.3348796121709862,-0.3984010898462415,-0.4600650377311521,-0.5195839500354337,-0.5766803221148671,-0.6310879443260529,-0.6825531432186539,-0.7308359642781241,-0.7757112907044197,-0.816969893010442,-0.8544194045464883,-0.8878852184023752,-0.917211301505453,-0.9422609221188204,-0.9629172873477992,-0.9790840876823228,-0.9906859460363308,-0.9976687691905392,-1,
	      1,0.9977668786231532,0.9910774881547801,0.9799617050365869,0.9644691750543766,0.9446690916079188,0.9206498866764288,0.8925188358598812,0.8604015792601394,0.8244415603417603,0.784799385278661,0.7416521056479576,0.6951924276746423,0.6456278515588024,0.5931797447293553,0.5380823531633727,0.4805817551866838,0.420934762428335,0.3594077728375128,0.2962755808856339,0.2318201502675284,0.1663293545831302,0.1000956916240984,0.03341497700767464,-0.03341497700767452,-0.1000956916240983,-0.16632935458313,-0.2318201502675283,-0.2962755808856338,-0.3594077728375127,-0.4209347624283349,-0.4805817551866837,-0.5380823531633726,-0.5931797447293552,-0.6456278515588024,-0.6951924276746423,-0.7416521056479576,-0.784799385278661,-0.8244415603417604,-0.8604015792601394,-0.8925188358598811,-0.9206498866764287,-0.9446690916079187,-0.9644691750543765,-0.9799617050365867,-0.99107748815478,-0.9977668786231532,-1,
	      1,0.9978589232386035,0.9914448613738104,0.9807852804032304,0.9659258262890683,0.9469301294951057,0.9238795325112867,0.8968727415326884,0.8660254037844387,0.8314696123025452,0.7933533402912352,0.7518398074789774,0.7071067811865476,0.6593458151000688,0.6087614290087207,0.5555702330196023,0.5000000000000001,0.4422886902190012,0.3826834323650898,0.3214394653031617,0.2588190451025207,0.1950903220161283,0.1305261922200517,0.06540312923014305,6.123031769111886e-17,-0.06540312923014292,-0.1305261922200516,-0.1950903220161282,-0.2588190451025206,-0.3214394653031616,-0.3826834323650897,-0.4422886902190011,-0.4999999999999998,-0.555570233019602,-0.6087614290087207,-0.6593458151000688,-0.7071067811865475,-0.7518398074789773,-0.7933533402912351,-0.8314696123025453,-0.8660254037844387,-0.8968727415326883,-0.9238795325112867,-0.9469301294951056,-0.9659258262890682,-0.9807852804032304,-0.9914448613738104,-0.9978589232386035,-1,
	      1,0.9979453927503363,0.9917900138232462,0.9815591569910653,0.9672948630390295,0.9490557470106686,0.9269167573460217,0.9009688679024191,0.8713187041233894,0.8380881048918406,0.8014136218679566,0.7614459583691344,0.7183493500977276,0.6723008902613168,0.6234898018587336,0.5721166601221697,0.5183925683105252,0.4625382902408354,0.4047833431223938,0.3453650544213076,0.2845275866310324,0.2225209339563144,0.1595998950333793,0.09602302590768189,0.03205157757165533,-0.03205157757165521,-0.09602302590768176,-0.1595998950333792,-0.2225209339563143,-0.2845275866310323,-0.3453650544213075,-0.4047833431223937,-0.4625382902408351,-0.518392568310525,-0.5721166601221694,-0.6234898018587335,-0.6723008902613169,-0.7183493500977275,-0.7614459583691344,-0.8014136218679565,-0.8380881048918406,-0.8713187041233892,-0.900968867902419,-0.9269167573460217,-0.9490557470106686,-0.9672948630390295,-0.9815591569910653,-0.9917900138232461,-0.9979453927503363,-1,
	      1,0.9987569212189223,0.9950307753654014,0.9888308262251285,0.9801724878485438,0.969077286229078,0.9555728057861407,0.9396926207859084,0.9214762118704076,0.9009688679024191,0.8782215733702285,0.8532908816321557,0.8262387743159949,0.7971325072229225,0.766044443118978,0.7330518718298263,0.6982368180860729,0.6616858375968594,0.6234898018587336,0.58374367223479,0.5425462638657594,0.5000000000000001,0.4562106573531631,0.4112871031306115,0.365341024366395,0.3184866502516844,0.2708404681430052,0.2225209339563144,0.1736481776669304,0.1243437046474853,0.07473009358642439,0.02493069173807303,-0.02493069173807291,-0.07473009358642427,-0.1243437046474852,-0.1736481776669303,-0.2225209339563143,-0.270840468143005,-0.3184866502516843,-0.3653410243663949,-0.4112871031306114,-0.456210657353163,-0.4999999999999998,-0.5425462638657593,-0.5837436722347896,-0.6234898018587335,-0.6616858375968595,-0.6982368180860727,-0.7330518718298263,-0.7660444431189779,-0.7971325072229225,-0.8262387743159947,-0.8532908816321556,-0.8782215733702284,-0.900968867902419,-0.9214762118704077,-0.9396926207859083,-0.9555728057861408,-0.9690772862290778,-0.9801724878485438,-0.9888308262251285,-0.9950307753654014,-0.9987569212189223,-1,
	      1,0.9996940572530831,0.9987764162142613,0.9972476383747747,0.9951086591716066,0.992360787415103,0.9890057044881307,0.9850454633172633,0.9804824871166253,0.9753195679051626,0.9695598647982465,0.9632069020746571,0.9562645670201275,0.9487371075487712,0.9406291296038439,0.9319455943394346,0.9226918150848067,0.9128734540932493,0.9024965190774262,0.8915673595333443,0.8800926628551884,0.8680794502434017,0.855535072408516,0.8424672050733575,0.8288838442763838,0.8147933014790244,0.8002041984800171,0.7851254621398549,0.7695663189185699,0.7535362892301956,0.7370451816173639,0.7201030867496006,0.7027203712489901,0.6849076713469912,0.6666758863762794,0.6480361721016052,0.6289999338937425,0.6095788197507079,0.5897847131705195,0.5696297258798572,0.5491261904230724,0.5282866526160837,0.507123863869773,0.4856507733875838,0.4638805202420894,0.4418264253353867,0.419501983248229,0.3969208539828873,0.3740968546047932,0.3510439507880777,0.3277762482701767,0.3043079842207362,0.2806535185300931,0.256827325022668,0.2328439826006414,0.2087181663233352,0.184464638427756,0.1600982392957975,0.1356338783736258,0.1110865250488047,0.08647119949074578,0.06180296346008408,0.03709691109260546,0.01236815966336301,-0.01236815966336288,-0.03709691109260534,-0.06180296346008397,-0.08647119949074565,-0.1110865250488045,-0.1356338783736257,-0.1600982392957974,-0.1844646384277558,-0.2087181663233351,-0.2328439826006413,-0.2568273250226679,-0.280653518530093,-0.3043079842207361,-0.3277762482701766,-0.3510439507880775,-0.3740968546047931,-0.3969208539828872,-0.4195019832482288,-0.4418264253353866,-0.4638805202420893,-0.4856507733875837,-0.5071238638697732,-0.5282866526160834,-0.5491261904230722,-0.5696297258798569,-0.5897847131705194,-0.6095788197507078,-0.6289999338937425,-0.6480361721016052,-0.6666758863762794,-0.6849076713469909,-0.70272037124899,-0.7201030867496003,-0.7370451816173638,-0.7535362892301954,-0.7695663189185699,-0.7851254621398549,-0.8002041984800169,-0.8147933014790243,-0.8288838442763837,-0.8424672050733574,-0.8555350724085159,-0.8680794502434017,-0.8800926628551884,-0.8915673595333444,-0.9024965190774261,-0.9128734540932492,-0.9226918150848066,-0.9319455943394345,-0.9406291296038438,-0.9487371075487712,-0.9562645670201275,-0.963206902074657,-0.9695598647982465,-0.9753195679051626,-0.9804824871166253,-0.9850454633172633,-0.9890057044881307,-0.992360787415103,-0.9951086591716066,-0.9972476383747747,-0.9987764162142613,-0.9996940572530831,-1,
	      1,0.9999241101148306,0.9996964519778716,0.9993170601430229,0.99878599219429,0.9981033287370441,0.997269173385788,0.9962836527484294,0.9951469164070644,0.9938591368952737,0.9924205096719357,0.9908312530915603,0.989091608371146,0.987201839553569,0.9851622334675065,0.9829730996839018,0.9806347704689777,0.9781476007338057,0.9755119679804366,0.9727282722446048,0.9697969360350095,0.9667184042691874,0.9634931442059831,0.9601216453746282,0.9566044195004408,0.9529420004271566,0.9491349440359013,0.9451838281608196,0.9410892525013715,0.9368518385313106,0.9324722294043558,0.9279510898565747,0.9232891061054893,0.918486985745923,0.9135454576426009,0.9084652718195237,0.9032471993461288,0.8978920322202582,0.8924005832479478,0.8867736859200619,0.8810121942857845,0.8751169828229927,0.8690889463055284,0.8629289996673897,0.8566380778638628,0.8502171357296142,0.8436671478337664,0.8369891083319778,0.8301840308155507,0.8232529481575873,0.8161969123562217,0.8090169943749475,0.8017142839800667,0.7942898895752861,0.7867449380334832,0.7790805745256705,0.7712979623471806,0.763398282741103,0.7553827347189938,0.747252534878891,0.7390089172206591,0.7306531329586932,0.7221864503320093,0.7136101544117524,0.7049255469061472,0.6961339459629267,0.6872366859692627,0.678235117349234,0.6691306063588582,0.6599245348787226,0.6506183002042422,0.6412133148335785,0.6317110062532509,0.6221128167214739,0.61242020304925,0.6026346363792564,0.5927576019625549,0.582790598933161,0.5727351400805053,0.5625927516198231,0.5523649729605059,0.5420533564724495,0.5316594672504361,0.5211848828765852,0.510631193180907,0.5000000000000001,0.4892929169339237,0.4785115691012865,0.4676575928925868,0.4567326357218406,0.4457383557765383,0.434676421765965,0.4235485126679244,0.4123563174739034,0.4011015349327188,0.3897858732926794,0.3784110500423103,0.3669787916496722,0.3554908333003182,0.3439489186339281,0.3323547994796596,0.3207102355902552,0.3090169943749475,0.2972768506312027,0.2854915862753422,0.273662990072083,0.2617928573630403,0.2498829897942308,0.2379351950426188,0.2259512865417477,0.2139330832064975,0.2018824091570104,0.1898010934418257,0.1776909697602686,0.16555387618413,0.1533916548786854,0.1412061518230915,0.1289992165302034,0.1167727017658563,0.1045284632676535,0.09226835946330202,0.07999425118854168,0.06770800140470754,0.05541147491597008,0.04310653808629573,0.03079505855617033,0.01847890495912992,0.006159946638138691,-0.006159946638138568,-0.01847890495912979,-0.0307950585561702,-0.04310653808629561,-0.05541147491596995,-0.06770800140470741,-0.07999425118854157,-0.09226835946330189,-0.1045284632676533,-0.1167727017658561,-0.1289992165302033,-0.1412061518230914,-0.1533916548786853,-0.1655538761841299,-0.1776909697602685,-0.1898010934418256,-0.2018824091570103,-0.2139330832064974,-0.2259512865417476,-0.2379351950426187,-0.2498829897942307,-0.2617928573630401,-0.2736629900720829,-0.285491586275342,-0.2972768506312026,-0.3090169943749473,-0.320710235590255,-0.3323547994796595,-0.343948918633928,-0.3554908333003181,-0.3669787916496721,-0.3784110500423102,-0.3897858732926793,-0.4011015349327187,-0.4123563174739033,-0.4235485126679243,-0.4346764217659649,-0.4457383557765382,-0.4567326357218405,-0.4676575928925867,-0.4785115691012864,-0.4892929169339235,-0.4999999999999998,-0.5106311931809068,-0.5211848828765848,-0.5316594672504362,-0.5420533564724495,-0.5523649729605058,-0.5625927516198231,-0.5727351400805052,-0.5827905989331609,-0.5927576019625548,-0.6026346363792563,-0.6124202030492498,-0.6221128167214738,-0.6317110062532507,-0.6412133148335781,-0.6506183002042422,-0.6599245348787227,-0.6691306063588582,-0.678235117349234,-0.6872366859692627,-0.6961339459629265,-0.7049255469061471,-0.7136101544117522,-0.7221864503320092,-0.7306531329586931,-0.739008917220659,-0.7472525348788908,-0.7553827347189935,-0.763398282741103,-0.7712979623471807,-0.7790805745256705,-0.7867449380334832,-0.794289889575286,-0.8017142839800666,-0.8090169943749473,-0.8161969123562216,-0.8232529481575871,-0.8301840308155505,-0.8369891083319777,-0.8436671478337662,-0.850217135729614,-0.8566380778638628,-0.8629289996673897,-0.8690889463055284,-0.8751169828229927,-0.8810121942857845,-0.8867736859200619,-0.8924005832479478,-0.897892032220258,-0.9032471993461288,-0.9084652718195236,-0.9135454576426008,-0.9184869857459229,-0.9232891061054892,-0.9279510898565747,-0.9324722294043558,-0.9368518385313106,-0.9410892525013715,-0.9451838281608195,-0.9491349440359012,-0.9529420004271565,-0.9566044195004407,-0.9601216453746281,-0.9634931442059831,-0.9667184042691874,-0.9697969360350094,-0.9727282722446048,-0.9755119679804366,-0.9781476007338057,-0.9806347704689777,-0.9829730996839018,-0.9851622334675065,-0.987201839553569,-0.989091608371146,-0.9908312530915603,-0.9924205096719357,-0.9938591368952736,-0.9951469164070644,-0.9962836527484294,-0.997269173385788,-0.9981033287370441,-0.99878599219429,-0.9993170601430229,-0.9996964519778716,-0.9999241101148306,-1,
	    };

	  unsigned order_array[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,64,128,256,};
	  unsigned order_array_length = 53;

	  // Construct the weights or knots as requested
	  if(weights)
	    construction_helper(weights_data_array, order_array, order_array_length);
	  else
	    construction_helper(knots_data_array, order_array, order_array_length);

	  break;
	}


      case 2 :
	{

	  double weights_data_array[] =
	    {
	      2,
	      0.6666666666666666,0.6666666666666666,
	      0.6666666666666667,0.6666666666666667,0.6666666666666666,
	      0.4254644007500071,0.574535599249993,0.574535599249993,0.425464400750007,
	      0.3111111111111111,0.4000000000000001,0.5777777777777778,0.4,0.3111111111111111,
	      0.2269152467244296,0.3267938603769863,0.4462908928985842,0.4462908928985841,0.3267938603769863,0.2269152467244296,
	      0.17796468096205,0.2476190476190477,0.3934638904665215,0.361904761904762,0.3934638904665215,0.2476190476190476,0.1779646809620499,
	      0.1397697435050226,0.2063696457302284,0.3142857142857144,0.3395748964790347,0.3395748964790348,0.3142857142857143,0.2063696457302284,0.1397697435050225,
	      0.1147810750857218,0.1654331942222275,0.2737903534857068,0.2790112502222169,0.333968253968254,0.279011250222217,0.2737903534857068,0.1654331942222275,0.1147810750857217,
	      0.09441954173982806,0.1411354380109716,0.2263866903636005,0.2530509772156453,0.2850073526699545,0.2850073526699544,0.2530509772156453,0.2263866903636005,0.1411354380109716,0.09441954173982806,
	      0.08004343893808848,0.1175565175565176,0.1987493987493988,0.212987012987013,0.2666617077670583,0.248003848003848,0.2666617077670583,0.212987012987013,0.1987493987493987,0.1175565175565176,0.08004343893808837,
	      0.0679557494725985,0.102289902955499,0.1688961414796902,0.1922697506411542,0.2300995573175295,0.2384888981335286,0.2384888981335286,0.2300995573175294,0.1922697506411542,0.1688961414796903,0.102289902955499,0.06795574947259848,
	      0.05895579755149852,0.08756398141964193,0.1499214887520146,0.1653493787446565,0.2123078777387938,0.2086251013741632,0.2345527488384631,0.2086251013741632,0.2123078777387938,0.1653493787446565,0.1499214887520146,0.08756398141964192,0.05895579755149849,
	      0.05120941158738113,0.07741843619139868,0.1301313630007154,0.149785348427595,0.1857742257742258,0.1964679635986112,0.2092132514200729,0.2092132514200729,0.1964679635986113,0.1857742257742258,0.1497853484275949,0.1301313630007154,0.07741843619139864,0.05120941158738113,
	      0.04521184009210766,0.06763919428824988,0.1167486898367929,0.1311133311133311,0.1710158015366053,0.1736306469815914,0.2003570018678276,0.1885669885669886,0.2003570018678276,0.1736306469815914,0.1710158015366053,0.1311133311133311,0.1167486898367929,0.06763919428824988,0.04521184009210756,
	      0.03995601467630659,0.06058044331486687,0.1030369110467109,0.1194572203372117,0.1516004910244863,0.1623668005175085,0.1796250702432303,0.183377048839679,0.183377048839679,0.1796250702432303,0.1623668005175085,0.1516004910244862,0.1194572203372117,0.1030369110467108,0.06058044331486689,0.03995601467630658,
	      0.03576294547534861,0.05376861364478319,0.09331989724146589,0.1060832459810415,0.1398288737366244,0.1450079332432274,0.170156550065808,0.1657284424250656,0.1806869963732709,0.1657284424250656,0.170156550065808,0.1450079332432274,0.1398288737366244,0.1060832459810415,0.09331989724146589,0.05376861364478316,0.03576294547534861,
	      0.03203589670062518,0.04867037737065753,0.08346231119164667,0.09724622435084129,0.1253653129286142,0.1354159146057867,0.1536315085885575,0.1589352104069294,0.1652372438563418,0.1652372438563418,0.1589352104069294,0.1536315085885574,0.1354159146057866,0.1253653129286141,0.09724622435084125,0.08346231119164664,0.04867037737065751,0.03203589670062514,
	      0.02899117856796175,0.04374159777859085,0.07621527014928345,0.08738509906729383,0.1160361757265782,0.122105022061837,0.1449403975209773,0.1444065011453591,0.1601327675088836,0.1520919809464701,0.1601327675088836,0.1444065011453591,0.1449403975209773,0.122105022061837,0.1160361757265782,0.08738509906729383,0.07621527014928343,0.04374159777859083,0.02899117856796173,
	      0.02625359736850488,0.03994423719383626,0.06890682602015491,0.08057541030252494,0.1050469697171128,0.1141543887597123,0.1318113112009618,0.1376099986766021,0.1468533151346628,0.1488439456259273,0.1488439456259273,0.1468533151346628,0.1376099986766021,0.1318113112009618,0.1141543887597122,0.1050469697171128,0.08057541030252489,0.06890682602015488,0.03994423719383625,0.02625359736850483,
	      0.02397383856658338,0.03626590512914538,0.06337243335687084,0.07311775732471802,0.09761619229654267,0.1038092671810341,0.1242477093196599,0.125774811917458,0.1411376273288309,0.1372227346381207,0.1469234458820724,0.1372227346381207,0.1411376273288309,0.125774811917458,0.1242477093196599,0.1038092671810341,0.09761619229654266,0.07311775732471801,0.06337243335687083,0.03626590512914538,0.02397383856658339,
	      0.02190455464575806,0.0333639489485321,0.05781181655165753,0.06778244516202234,0.08910706617287298,0.09726273860041823,0.1137589459703456,0.119545912894527,0.1299651349270397,0.1329678365010923,0.1365295996257343,0.1365295996257343,0.1329678365010923,0.1299651349270397,0.119545912894527,0.1137589459703456,0.09726273860041824,0.08910706617287292,0.06778244516202231,0.05781181655165753,0.0333639489485321,0.02190455464575804,
	      0.02015375303064504,0.03054757367634667,0.05349707940126579,0.06201709230380654,0.08313354854281184,0.08910537073143078,0.1073090994562439,0.1099001146921459,0.1243997798255167,0.1229739938713348,0.1332458701782994,0.1274334485803055,0.1332458701782993,0.1229739938713348,0.1243997798255167,0.1099001146921459,0.1073090994562439,0.08910537073143077,0.08313354854281184,0.06201709230380654,0.05349707940126579,0.03054757367634662,0.02015375303064505,
	      0.01855193395899458,0.02828130490708779,0.04917232037946873,0.05777091127970643,0.07642920933965271,0.08370544618774016,0.0988545813059163,0.1043942762387843,0.1150611766404668,0.1185278803542791,0.1240355141196194,0.1252154452882839,0.1252154452882839,0.1240355141196194,0.1185278803542791,0.1150611766404668,0.1043942762387843,0.09885458130591626,0.08370544618774012,0.07642920933965272,0.05777091127970645,0.04917232037946868,0.0282813049070878,0.01855193395899454,
	      0.01717844376189322,0.0260779238152801,0.04574727838905179,0.05322795913720512,0.07157635802269194,0.07718190511884848,0.09339050657215989,0.0964913608224936,0.1099419655383285,0.1100254633749169,0.1202730868326568,0.1169953877312557,0.1237847217664358,0.1169953877312557,0.1202730868326568,0.110025463374917,0.1099419655383285,0.0964913608224936,0.09339050657215987,0.07718190511884848,0.07157635802269194,0.05322795913720509,0.04574727838905176,0.0260779238152801,0.01717844376189321,
	      0.01591332921183207,0.02427502581977982,0.0423197833461823,0.04979970525744667,0.06621051307180502,0.07270361446732157,0.0865069328545057,0.09169982784977403,0.1021334352819075,0.1057561445572326,0.1122516577587538,0.1141125167497899,0.1163175137736692,0.1163175137736692,0.1141125167497899,0.1122516577587537,0.1057561445572326,0.1021334352819075,0.09169982784977403,0.08650693285450567,0.07270361446732156,0.06621051307180499,0.04979970525744665,0.04231978334618228,0.0242750258197798,0.01591332921183206,
	      0.01481614149931246,0.02251926453352875,0.0395577665408326,0.04615990609765128,0.06222600091663893,0.06741744177177753,0.0818786775581576,0.08517753023131336,0.09754750954401821,0.09854216168205163,0.1084504674071349,0.1068392018476736,0.114041955052424,0.1096519506349709,0.114041955052424,0.1068392018476736,0.1084504674071348,0.09854216168205161,0.09754750954401821,0.08517753023131337,0.08187867755815759,0.06741744177177753,0.06222600091663892,0.04615990609765127,0.0395577665408326,0.02251926453352874,0.01481614149931244,
	      0.0137996515263834,0.02106189296772671,0.036796782165302,0.04335580630675379,0.05787063269001485,0.06367769348750806,0.07621705148715324,0.08103223120569746,0.09099420565800084,0.09460096697879163,0.1015146093871649,0.1037475212632094,0.1072874948277409,0.1080434600485528,0.1080434600485528,0.1072874948277409,0.1037475212632094,0.1015146093871648,0.09460096697879161,0.09099420565800082,0.08103223120569748,0.07621705148715323,0.06367769348750807,0.05787063269001487,0.04335580630675377,0.03679678216530197,0.02106189296772671,0.01379965152638338,
	      0.01290939845532382,0.01964048817661894,0.03453833573414108,0.04039674463546401,0.05456556826620402,0.05934196335344636,0.07228505216864153,0.07560608057154682,0.08693727797722474,0.08847181589164776,0.09788506159863314,0.09737515162138188,0.1046509242391442,0.1019263764395495,0.1069395217420648,0.1019263764395495,0.1046509242391442,0.09737515162138188,0.09788506159863311,0.08847181589164775,0.08693727797722474,0.07560608057154682,0.07228505216864149,0.05934196335344637,0.05456556826620402,0.04039674463546397,0.03453833573414106,0.01964048817661894,0.0129093984553238,
	      0.01208045752791377,0.01844590602384014,0.03228242663725728,0.03807611645253615,0.05098561202169699,0.05619569570263368,0.06758273405721409,0.07202345195119504,0.08140827509650779,0.084905329348304,0.09189922978367299,0.09431228795431945,0.09862708088750186,0.09985857822758204,0.1013168183278248,0.1013168183278248,0.09985857822758204,0.09862708088750186,0.09431228795431942,0.09189922978367299,0.08490532934830398,0.0814082750965078,0.07202345195119504,0.06758273405721407,0.05619569570263367,0.05098561202169698,0.03807611645253613,0.03228242663725728,0.01844590602384013,0.01208045752791375,
	      0.01134823217605378,0.01727921091958246,0.03041295089033243,0.03563930770094848,0.0482177970508505,0.05259998209088158,0.06422714375839272,0.06747254141733675,0.07783883465448722,0.07967978344183978,0.08853257265492616,0.08875107865462227,0.09589827865186744,0.09433729173238946,0.09965322242115431,0.09622354356866958,0.09965322242115431,0.09433729173238946,0.09589827865186745,0.08875107865462226,0.08853257265492616,0.07967978344183978,0.07783883465448721,0.06747254141733673,0.06422714375839272,0.0525999820908816,0.04821779705085048,0.03563930770094848,0.03041295089033242,0.01727921091958242,0.01134823217605374,
	      0.0106634169535045,0.0162879800866775,0.02854654675922566,0.03369838504308589,0.04524179987167728,0.04993328522400221,0.06028544803676913,0.0643712271584536,0.07314606500401978,0.07648505889677426,0.08336147485750971,0.08583551153649585,0.09056331853887302,0.0920840979019617,0.09449166255976237,0.0950047215712076,0.0950047215712076,0.09449166255976234,0.09208409790196172,0.09056331853887301,0.08583551153649584,0.08336147485750968,0.07648505889677426,0.07314606500401978,0.0643712271584536,0.06028544803676909,0.04993328522400221,0.04524179987167727,0.03369838504308589,0.02854654675922566,0.01628798008667748,0.01066341695350447,
	      0.01005394129329871,0.01531868526330632,0.02698207707495557,0.03166859742212745,0.04290325444292508,0.04692112912133165,0.05740682489171495,0.06052423557326231,0.07001043992376978,0.07200968298596679,0.08028736361588908,0.08098501266031657,0.0878884028688835,0.08714410727012456,0.09255501228405577,0.09027703455204927,0.09412839751204552,0.09027703455204927,0.09255501228405577,0.08714410727012456,0.08788840286888353,0.08098501266031657,0.08028736361588908,0.07200968298596677,0.07001043992376978,0.06052423557326228,0.05740682489171494,0.04692112912133163,0.04290325444292508,0.03166859742212747,0.02698207707495555,0.01531868526330631,0.01005394129329868,
	      0.009481697868816842,0.01448723312606658,0.02542074192634194,0.03002968184645259,0.04040408530551475,0.04464460707118954,0.05407423100100589,0.05783148882152254,0.06600270754511778,0.0691617724269326,0.07580845764936581,0.07827002093698382,0.08317706447934851,0.0848630196225796,0.08787200453300111,0.08872864822714903,0.08974253761261117,0.08974253761261118,0.08872864822714906,0.08787200453300111,0.08486301962257961,0.08317706447934851,0.07827002093698379,0.07580845764936581,0.06916177242693265,0.06600270754511776,0.05783148882152252,0.05407423100100589,0.04464460707118954,0.04040408530551474,0.03002968184645259,0.02542074192634193,0.01448723312606657,0.00948169786881683,
	      0.008969009637343924,0.01367326118802712,0.02409868335990946,0.02832160154564107,0.03841199745605253,0.04209809706368062,0.0515912889241525,0.05455508228626558,0.06324639564176926,0.06530961372469007,0.07302537471291884,0.07403373290509818,0.08063178693188122,0.08046193903679091,0.08583478461296758,0.08439873850403896,0.0884763930087191,0.08572443892010649,0.08847639300871908,0.08439873850403898,0.08583478461296756,0.08046193903679093,0.08063178693188122,0.07403373290509817,0.07302537471291883,0.06530961372469007,0.06324639564176925,0.05455508228626557,0.05159128892415248,0.04209809706368062,0.03841199745605251,0.02832160154564106,0.02409868335990942,0.01367326118802712,0.008969009637343924,
	      0.008485964887575222,0.01296906765771568,0.02277960245709249,0.02692568779409434,0.03629394571839896,0.04014134172864429,0.04875116923728127,0.05220832837181642,0.05980263273472401,0.06277529464002612,0.06913248440459831,0.07153711722482871,0.07647298425985767,0.07824132165663772,0.08161323029741263,0.08269485380905874,0.08440548145846187,0.08476949166177572,0.08476949166177571,0.08440548145846186,0.08269485380905872,0.08161323029741262,0.07824132165663771,0.07647298425985768,0.07153711722482871,0.06913248440459831,0.0627752946400261,0.05980263273472396,0.05220832837181639,0.04875116923728125,0.04014134172864431,0.03629394571839895,0.02692568779409433,0.02277960245709248,0.01296906765771568,0.008485964887575189,
	      0.008050615920773009,0.01227896784960597,0.0216525554890753,0.02547511411251985,0.03458415224986446,0.0379704661762646,0.04659782858051167,0.04939811927368789,0.05737509633150602,0.05944237766262679,0.06662394091358095,0.06782819777774832,0.07409269567710451,0.0743264571097525,0.07957787506790147,0.07875974015614583,0.08292996500382845,0.08100704636813486,0.0840575765587353,0.08100704636813486,0.08292996500382847,0.07875974015614581,0.07957787506790145,0.0743264571097525,0.0740926956771045,0.0678281977777483,0.06662394091358095,0.05944237766262676,0.05737509633150602,0.04939811927368786,0.04659782858051167,0.03797046617626459,0.03458415224986448,0.02547511411251984,0.02165255548907531,0.01227896784960596,0.008050615920773016,
	      0.007639159487160807,0.01167737022893692,0.02052823733984574,0.02427680287201771,0.03277414778150618,0.03627770993377577,0.04415954317522237,0.0473445396715912,0.05439829105262065,0.05718688970395116,0.06322707708561662,0.06554883585734199,0.07041783245293887,0.07221344065050635,0.07578456016652249,0.07700793068689343,0.07918838127419786,0.07980804466223902,0.08054120591711489,0.08054120591711489,0.07980804466223899,0.07918838127419785,0.07700793068689342,0.07578456016652249,0.07221344065050635,0.07041783245293884,0.06554883585734199,0.0632270770856166,0.05718688970395116,0.05439829105262067,0.0473445396715912,0.04415954317522235,0.03627770993377575,0.03277414778150618,0.02427680287201771,0.02052823733984575,0.01167737022893689,0.007639159487160785,
	      0.007266349468477857,0.01108725733319443,0.01955979765120058,0.02303470755608177,0.03129643014056387,0.03441291810896161,0.04228211434959982,0.04491821973608899,0.05225464886556719,0.05428835344418971,0.06097024206498124,0.06229163952772455,0.06821484297509875,0.06873066898611994,0.07381028453829833,0.07344674635442013,0.07761888689053432,0.07632367854935819,0.07954691587619088,0.07729059516669598,0.07954691587619085,0.07632367854935819,0.07761888689053431,0.07344674635442011,0.07381028453829831,0.06873066898611994,0.06821484297509874,0.06229163952772455,0.06097024206498123,0.0542883534441897,0.05225464886556719,0.04491821973608899,0.04228211434959982,0.0344129181089616,0.03129643014056386,0.02303470755608176,0.01955979765120056,0.01108725733319442,0.007266349468477849,
	      0.006913009907557729,0.01056928122098102,0.018593810542238,0.02199864321297218,0.0297379668966998,0.03293969431353681,0.04017454191713343,0.0431136261239712,0.04966681023334055,0.05227856019836072,0.05799396223631948,0.06021876616161181,0.0649613490285489,0.06674780265628481,0.07040587529811353,0.0717124891822708,0.07420002855740443,0.0749963825192253,0.07625493435236795,0.07652246544106167,0.07652246544106167,0.07625493435236794,0.07499638251922532,0.0742000285574044,0.07171248918227079,0.07040587529811353,0.06674780265628484,0.06496134902854889,0.06021876616161183,0.05799396223631947,0.05227856019836073,0.04966681023334055,0.0431136261239712,0.04017454191713343,0.0329396943135368,0.02973796689669979,0.02199864321297219,0.01859381054223799,0.01056928122098102,0.006913009907557705,
	      0.006591315077489216,0.01006075474360836,0.01775565710217204,0.02092708253066298,0.02845263603347119,0.0313266335078704,0.03852930515485761,0.04100580170428297,0.04776808903373336,0.04974512608041667,0.05596420587640315,0.0573485199292903,0.06293507068452031,0.06364582699898329,0.06852516431787356,0.06849624390843766,0.07260970337607635,0.07179135780030654,0.07509748979131806,0.07345753084492131,0.07593297100660941,0.07345753084492131,0.07509748979131806,0.07179135780030654,0.07260970337607636,0.06849624390843766,0.06852516431787355,0.06364582699898327,0.06293507068452034,0.05734851992929032,0.05596420587640312,0.04974512608041667,0.04776808903373336,0.04100580170428295,0.03852930515485761,0.03132663350787039,0.02845263603347119,0.02092708253066298,0.01775565710217203,0.01006075474360834,0.006591315077489217,
	      0.006285645605961367,0.00961161162263263,0.01691968436936594,0.02002542217440113,0.02710149879488746,0.03003723460775803,0.03669611864168228,0.03941335299765841,0.04550622327360962,0.047950847346684,0.05334556486643315,0.05546692961530304,0.06004754315517039,0.06180110989724233,0.06546951011611693,0.06681825570794783,0.06949599739769576,0.07041136889329723,0.07204123731122326,0.07250382961182605,0.07305101399310329,0.07305101399310328,0.07250382961182605,0.07204123731122326,0.07041136889329723,0.06949599739769576,0.06681825570794783,0.06546951011611694,0.06180110989724233,0.06004754315517036,0.05546692961530304,0.05334556486643315,0.047950847346684,0.0455062232736096,0.03941335299765839,0.03669611864168229,0.03003723460775803,0.02710149879488746,0.02002542217440112,0.01691968436936591,0.009611611622632629,0.006285645605961365,
	      0.006006135231203551,0.009170312612119274,0.0161895253839546,0.01909467283594488,0.02597689536297483,0.02863297890125528,0.03524733856613677,0.03757166754222491,0.04381898010997718,0.0457258214710746,0.05151877848816915,0.05292865979063543,0.05819044478278049,0.05903327269997125,0.06369834357367429,0.06391526768490031,0.06793043226541126,0.067475203804764,0.07080059857118276,0.06964058197055174,0.07225043464127932,0.07036730741962821,0.07225043464127932,0.06964058197055173,0.07080059857118276,0.06747520380476402,0.06793043226541125,0.0639152676849003,0.06369834357367429,0.05903327269997125,0.05819044478278049,0.05292865979063541,0.05151877848816915,0.0457258214710746,0.04381898010997717,0.03757166754222491,0.03524733856613675,0.02863297890125528,0.02597689536297484,0.01909467283594489,0.01618952538395459,0.009170312612119281,0.006006135231203539,
	      0.005739934233779501,0.008778340167961822,0.01546129420862699,0.01830526276191744,0.02479810897582393,0.02749852502438495,0.03364352922330797,0.03616065733017996,0.04183193018655375,0.04412024043822421,0.04920532252601977,0.05122159810584373,0.05562062953699962,0.05732624021285236,0.06095315946261697,0.06231522899653321,0.06509920295910258,0.0660914001052809,0.06797810618894662,0.06858122072157533,0.06953386153850585,0.06973620709496263,0.06973620709496264,0.06953386153850587,0.06858122072157533,0.0679781061889466,0.0660914001052809,0.06509920295910257,0.0623152289965332,0.06095315946261695,0.05732624021285234,0.0556206295369996,0.05122159810584374,0.04920532252601972,0.04412024043822421,0.04183193018655375,0.03616065733017996,0.03364352922330797,0.02749852502438494,0.02479810897582392,0.01830526276191745,0.015461294208627,0.008778340167961807,0.00573993423377948,
	      0.005495545384377011,0.008392931969689715,0.01482139498950876,0.01749178198942053,0.02380875667279853,0.02626883895357505,0.03236186196378257,0.03454287291744323,0.04032763869292086,0.04215705785387144,0.04755902587315353,0.04896883855405833,0.05392173327708531,0.05485106822380276,0.05929740042604561,0.05969406222904307,0.063585964213561,0.06340755178249223,0.06670757465406607,0.06592233430127635,0.06860410175345832,0.06719155011421622,0.06924022642070718,0.06719155011421622,0.06860410175345831,0.06592233430127635,0.06670757465406607,0.06340755178249223,0.06358596421356098,0.05969406222904307,0.05929740042604559,0.05485106822380274,0.05392173327708531,0.04896883855405834,0.04755902587315353,0.04215705785387142,0.04032763869292087,0.03454287291744322,0.03236186196378257,0.02626883895357505,0.02380875667279852,0.01749178198942054,0.01482139498950873,0.008392931969689706,0.00549554538437698,
	      0.005262301525422354,0.008048834279917184,0.01418319533113695,0.01679682652151323,0.02277437716897657,0.02526580703067208,0.03095116014626824,0.03328766752360139,0.03857361987813569,0.04071667570525296,0.04550700375607543,0.04741957322409165,0.0516279854655002,0.05327650011412977,0.05682749527673823,0.05818283224661605,0.06101282238335227,0.06205096245868374,0.06410931911791001,0.06481183375254135,0.06606175201919968,0.06641615924361087,0.06683529583065421,0.0668352958306542,0.06641615924361086,0.06606175201919968,0.06481183375254135,0.06410931911791001,0.06205096245868374,0.06101282238335227,0.05818283224661608,0.05682749527673821,0.05327650011412979,0.0516279854655002,0.04741957322409163,0.04550700375607542,0.04071667570525296,0.03857361987813569,0.0332876675236014,0.03095116014626819,0.02526580703067207,0.02277437716897656,0.01679682652151322,0.01418319533113694,0.008048834279917171,0.005262301525422339,
	      0.005047389355437568,0.007710270779230353,0.01361930770862348,0.01608177702367659,0.02189960553047587,0.02418313614079593,0.02981238861863938,0.03185945754919595,0.03722800744076079,0.03897692399709671,0.04402079480640794,0.04541309619886743,0.05007490614734521,0.05105761476843934,0.05528690484192644,0.05581379976216085,0.05956768135324461,0.05960022278766258,0.06284402540180421,0.06235207164197722,0.06505989647310334,0.0640222478079962,0.06617739019457167,0.06458216734112099,0.06617739019457165,0.06402224780799619,0.06505989647310334,0.06235207164197721,0.06284402540180421,0.05960022278766256,0.0595676813532446,0.05581379976216085,0.05528690484192644,0.05105761476843935,0.05007490614734521,0.04541309619886742,0.04402079480640791,0.0389769239970967,0.0372280074407608,0.03185945754919595,0.02981238861863939,0.02418313614079592,0.02189960553047586,0.01608177702367657,0.01361930770862345,0.007710270779230359,0.005047389355437565,
	      0.004841880287512181,0.007406565289467493,0.0130569171471642,0.01546685463351595,0.02098709965711286,0.02329223473677195,0.02856533288909778,0.03073860009929806,0.03567268914435337,0.03768130908411993,0.04219363256974867,0.04400573152711525,0.04802145642509568,0.04960779460480837,0.05306061383768573,0.05439541470616858,0.05722842996202266,0.05828993058529047,0.0604565051769342,0.06122736734660083,0.06269185556013554,0.06315947563476693,0.06389779044692996,0.06405451864828351,0.06405451864828353,0.06389779044692995,0.06315947563476694,0.06269185556013554,0.06122736734660082,0.0604565051769342,0.05828993058529047,0.05722842996202265,0.05439541470616857,0.05306061383768572,0.04960779460480839,0.04802145642509567,0.04400573152711525,0.04219363256974867,0.03768130908411992,0.03567268914435336,0.03073860009929807,0.02856533288909778,0.02329223473677195,0.02098709965711285,0.01546685463351594,0.01305691714716418,0.007406565289467486,0.004841880287512158,
	      0.004651889712821409,0.007107557645826528,0.01255749036742756,0.01483502315356996,0.0202100071623087,0.02233421292731736,0.02754936662972879,0.02947186954872786,0.03446510882377789,0.03613315111990666,0.04084928745120459,0.0422124004578068,0.04660157139762206,0.04761352839864058,0.05163138243383895,0.05225126336200878,0.05585946112266511,0.05605242056147462,0.05921916049241727,0.05895702981621391,0.06165751356062751,0.06091927073124629,0.06313607562386298,0.06190819064460767,0.06363153370870064,0.06190819064460767,0.06313607562386295,0.06091927073124629,0.0616575135606275,0.0589570298162139,0.05921916049241727,0.05605242056147462,0.05585946112266513,0.0522512633620088,0.05163138243383894,0.04761352839864056,0.04660157139762206,0.04221240045780677,0.04084928745120459,0.03613315111990667,0.0344651088237779,0.02947186954872787,0.02754936662972876,0.02233421292731735,0.02021000716230869,0.01483502315356997,0.01255749036742756,0.007107557645826511,0.004651889712821389,
	      0.004469887659462549,0.006838167099685352,0.01205938276384572,0.01428835078102051,0.01940104893289329,0.02153951662309917,0.02644184520202902,0.0284673387473048,0.03308011963551493,0.03496461078497344,0.03921631813289153,0.04093225935122131,0.04475776003690031,0.04627961308215733,0.04962057676666541,0.05092552169667786,0.05373111526234587,0.05479951348644096,0.05702709655938406,0.05784283816915334,0.05945857582280339,0.06000934714413022,0.06098870636776487,0.06126618830654191,0.06159430158509293,0.06159430158509292,0.0612661883065419,0.06098870636776487,0.06000934714413022,0.05945857582280336,0.05784283816915332,0.05702709655938404,0.05479951348644096,0.05373111526234586,0.05092552169667786,0.0496205767666654,0.04627961308215735,0.0447577600369003,0.04093225935122128,0.03921631813289155,0.03496461078497343,0.0330801196355149,0.0284673387473048,0.02644184520202901,0.02153951662309918,0.01940104893289328,0.01428835078102049,0.01205938276384572,0.006838167099685343,0.00446988765946255,
	      0.002752657295060779,0.004212899622558557,0.00744237421209495,0.008826922345720096,0.01202263597568563,0.01336931863002165,0.01648649930872011,0.01778882835564131,0.02079540226874629,0.02204284787202134,0.0249097711956744,0.02609130391955806,0.02879139689416383,0.02989627105492565,0.03240411866611649,0.03342216927811872,0.03571424274924558,0.03663605183564141,0.03869088256495642,0.03950789736370321,0.04130625677136774,0.04201088394914206,0.04353595278806144,0.04412163674711546,0.04535915668202571,0.0458204448426853,0.04675884849987249,0.04709144455314446,0.04772196172653019,0.04792276712606441,0.04823950560701255,0.04830664929860381,0.04830664929860381,0.04823950560701254,0.0479227671260644,0.0477219617265302,0.04709144455314446,0.04675884849987248,0.04582044484268529,0.04535915668202571,0.04412163674711546,0.04353595278806141,0.04201088394914206,0.04130625677136772,0.03950789736370321,0.03869088256495642,0.03663605183564141,0.03571424274924555,0.03342216927811871,0.03240411866611651,0.02989627105492564,0.02879139689416383,0.02609130391955806,0.02490977119567439,0.02204284787202133,0.02079540226874629,0.01778882835564129,0.0164864993087201,0.01336931863002164,0.01202263597568563,0.008826922345720091,0.007442374212094945,0.004212899622558559,0.002752657295060761,
	      0.0006991477762184043,0.00107058193501479,0.001895131463711382,0.00225041391376509,0.003076491218610577,0.003427636737070387,0.004249507799752202,0.004597217358326011,0.005412186726809599,0.005756043545939279,0.00656193578766951,0.006901277619412231,0.00769607939128861,0.008030171724563325,0.008811947626817708,0.009140034854651766,0.009906902968440613,0.0102282280577189,0.01097835297947463,0.0112921668145019,0.01202375888338542,0.01232932560952532,0.01304064261669361,0.01333724322279695,0.01402659319001972,0.01431352822668802,0.01497927265120654,0.01525586447998761,0.01589642176231531,0.016162016523011,0.01677586543290669,0.01702983482273073,0.01761551792330198,0.01785726083673554,0.01841338781889242,0.01864233187429905,0.01916758277080687,0.01938318573776566,0.0198763139956282,0.02007806513027563,0.0205379005257456,0.02072532181763952,0.02115077320160465,0.02132342053344456,0.02171347839720571,0.02187094261748199,0.0222246814705303,0.02236658937845357,0.0226831699310421,0.0228091851727136,0.02308785631696747,0.0231976801915661,0.02343778077567671,0.02353115295038583,0.02373211334114547,0.02380881247357605,0.02397015590316382,0.0240300001701221,0.02415134386367198,0.02419419139525015,0.02427524747632949,0.02430099669445606,0.02434157286616803,0.02435016272693003,0.02435016272693004,0.02434157286616802,0.02430099669445606,0.02427524747632949,0.02419419139525016,0.02415134386367198,0.02403000017012211,0.02397015590316382,0.02380881247357605,0.02373211334114546,0.02353115295038583,0.02343778077567671,0.0231976801915661,0.02308785631696747,0.0228091851727136,0.0226831699310421,0.02236658937845358,0.02222468147053029,0.021870942617482,0.02171347839720571,0.02132342053344456,0.02115077320160464,0.02072532181763953,0.0205379005257456,0.02007806513027563,0.01987631399562819,0.01938318573776566,0.01916758277080687,0.01864233187429905,0.01841338781889241,0.01785726083673553,0.01761551792330197,0.01702983482273072,0.01677586543290668,0.016162016523011,0.01589642176231532,0.0152558644799876,0.01497927265120653,0.01431352822668802,0.01402659319001972,0.01333724322279694,0.0130406426166936,0.01232932560952533,0.01202375888338542,0.0112921668145019,0.01097835297947463,0.0102282280577189,0.009906902968440601,0.009140034854651764,0.008811947626817704,0.008030171724563325,0.007696079391288603,0.006901277619412227,0.006561935787669511,0.005756043545939281,0.005412186726809588,0.004597217358326008,0.004249507799752202,0.003427636737070386,0.003076491218610567,0.002250413913765088,0.001895131463711381,0.001070581935014791,0.0006991477762183977,
	      0.0001761673040739866,0.0002697942661301234,0.0004778353095542702,0.0005675900250460902,0.0007766669504901253,0.0008657333309864466,0.001074771264345746,0.00116348179591581,0.00137216834784244,0.001460572404503794,0.001668722069769264,0.001756805209785221,0.001964268112189649,0.002051995258620697,0.002258634883670134,0.002345962765886926,0.002551648745222864,0.002638530412456905,0.002843135734497978,0.002929522476214617,0.003132922277873141,0.003218764545819178,0.003420835552848029,0.003506083449830412,0.00370670371243854,0.003791307274861531,0.003990356049074605,0.00407426542378007,0.004271623129523452,0.004354788693068941,0.00455033691477593,0.004632709359712384,0.004826330871488348,0.004907861272855531,0.005099440078271946,0.005180079947752517,0.005369501328545211,0.005449202660630764,0.0056363532308715,0.005715068543672725,0.005899836307286121,0.005977518679625424,0.006159793089889269,0.006236396195720363,0.00641606821585129,0.00649154635668635,0.006668508520901795,0.006742816656697344,0.006916963131328751,0.006990056910133889,0.007161283554486764,0.007233119341060038,0.007401323767796065,0.007471858671332664,0.007636940306203786,0.00770613220727012,0.007867992348071957,0.007935799924814493,0.008094341799452811,0.008160724553126383,0.008315853376709465,0.008380771656555212,0.008532394687438632,0.008595809714930918,0.008743836309651408,0.008805710202125283,0.008950051869168054,0.009010347662832992,0.009150918115182946,0.009209599787524583,0.009346314993956429,0.009403347485524922,0.009536125720590649,0.00959147495617222,0.009720236848847884,0.009773869758014184,0.009898538338970081,0.009950422875999284,0.01007092362345997,0.01012102878662238,0.01023728967078481,0.01028558552098525,0.01039753704696511,0.01044399472573422,0.01055156997501173,0.01059616172183791,0.01069929639217625,0.01074199556116997,0.01084062800498038,0.01088140908086282,0.01097548034199202,0.01101431895539936,0.01110377280431598,0.01114064574641222,0.01122542871376985,0.01126031395015968,0.01134037535871613,0.01137325204265071,0.01144854403752308,0.01147939252239153,0.01154987009962875,0.01157867195072823,0.01164429298418317,0.01167103098976151,0.01173175625624645,0.01175641443781041,0.01181220764052024,0.01183477126240436,0.0118855990525931,0.0119060546307832,0.01195188662768093,0.01197022193788746,0.01201103074684525,0.01202723483182172,0.0120629960606736,0.0120770592367762,0.01210775151040818,0.01211966537339251,0.01214527034651001,0.01215502777656161,0.01217553014464717,0.01218312531064329,0.0121985128190983,0.01220394118209795,0.01221420463356264,0.01221746294952333,0.01222259620937062,0.01222368253108986,0.01222368253108986,0.01222259620937062,0.01221746294952333,0.01221420463356264,0.01220394118209795,0.0121985128190983,0.01218312531064329,0.01217553014464716,0.01215502777656161,0.01214527034651001,0.01211966537339251,0.01210775151040818,0.0120770592367762,0.0120629960606736,0.01202723483182172,0.01201103074684525,0.01197022193788746,0.01195188662768093,0.0119060546307832,0.0118855990525931,0.01183477126240436,0.01181220764052024,0.01175641443781041,0.01173175625624645,0.01167103098976151,0.01164429298418318,0.01157867195072823,0.01154987009962875,0.01147939252239153,0.01144854403752308,0.01137325204265071,0.01134037535871612,0.01126031395015968,0.01122542871376984,0.01114064574641221,0.01110377280431598,0.01101431895539936,0.01097548034199202,0.01088140908086282,0.01084062800498038,0.01074199556116997,0.01069929639217624,0.0105961617218379,0.01055156997501172,0.01044399472573422,0.01039753704696511,0.01028558552098525,0.01023728967078481,0.01012102878662238,0.01007092362345996,0.009950422875999291,0.009898538338970073,0.009773869758014185,0.009720236848847884,0.00959147495617222,0.009536125720590651,0.009403347485524923,0.00934631499395642,0.009209599787524583,0.009150918115182946,0.009010347662832995,0.008950051869168054,0.008805710202125283,0.00874383630965141,0.008595809714930918,0.008532394687438632,0.008380771656555213,0.008315853376709467,0.008160724553126374,0.008094341799452813,0.007935799924814494,0.007867992348071952,0.00770613220727012,0.007636940306203788,0.007471858671332663,0.007401323767796064,0.007233119341060035,0.007161283554486762,0.00699005691013389,0.006916963131328751,0.006742816656697343,0.006668508520901792,0.006491546356686349,0.00641606821585129,0.006236396195720363,0.00615979308988927,0.005977518679625424,0.005899836307286114,0.005715068543672723,0.005636353230871502,0.005449202660630764,0.005369501328545206,0.005180079947752517,0.005099440078271944,0.00490786127285553,0.00482633087148835,0.004632709359712381,0.004550336914775928,0.004354788693068939,0.004271623129523454,0.004074265423780073,0.003990356049074601,0.00379130727486153,0.003706703712438538,0.003506083449830412,0.003420835552848029,0.003218764545819177,0.00313292227787314,0.002929522476214615,0.002843135734497977,0.002638530412456905,0.002551648745222862,0.002345962765886924,0.002258634883670132,0.002051995258620698,0.00196426811218965,0.001756805209785217,0.001668722069769263,0.001460572404503792,0.001372168347842439,0.001163481795915811,0.001074771264345741,0.0008657333309864439,0.0007766669504901241,0.0005675900250460908,0.0004778353095542713,0.0002697942661301205,0.0001761673040739835,
	    };

	  double knots_data_array[] =
	    {
	      0,
	      -0.4999999999999998,0.5000000000000001,
	      -0.7071067811865475,6.123031769111886e-17,0.7071067811865476,
	      -0.8090169943749473,-0.3090169943749473,0.3090169943749475,0.8090169943749475,
	      -0.8660254037844387,-0.4999999999999998,6.123031769111886e-17,0.5000000000000001,0.8660254037844387,
	      -0.900968867902419,-0.6234898018587335,-0.2225209339563143,0.2225209339563144,0.6234898018587336,0.9009688679024191,
	      -0.9238795325112867,-0.7071067811865475,-0.3826834323650897,6.123031769111886e-17,0.3826834323650898,0.7071067811865476,0.9238795325112867,
	      -0.9396926207859083,-0.7660444431189779,-0.4999999999999998,-0.1736481776669303,0.1736481776669304,0.5000000000000001,0.766044443118978,0.9396926207859084,
	      -0.9510565162951535,-0.8090169943749473,-0.587785252292473,-0.3090169943749473,6.123031769111886e-17,0.3090169943749475,0.5877852522924731,0.8090169943749475,0.9510565162951535,
	      -0.9594929736144974,-0.8412535328311811,-0.654860733945285,-0.4154150130018863,-0.142314838273285,0.1423148382732851,0.4154150130018864,0.6548607339452851,0.8412535328311812,0.9594929736144974,
	      -0.9659258262890682,-0.8660254037844387,-0.7071067811865475,-0.4999999999999998,-0.2588190451025206,6.123031769111886e-17,0.2588190451025207,0.5000000000000001,0.7071067811865476,0.8660254037844387,0.9659258262890683,
	      -0.970941817426052,-0.8854560256532098,-0.7485107481711012,-0.5680647467311557,-0.3546048870425355,-0.1205366802553229,0.120536680255323,0.3546048870425356,0.5680647467311559,0.7485107481711011,0.8854560256532099,0.970941817426052,
	      -0.9749279121818236,-0.900968867902419,-0.7818314824680298,-0.6234898018587335,-0.4338837391175581,-0.2225209339563143,6.123031769111886e-17,0.2225209339563144,0.4338837391175582,0.6234898018587336,0.7818314824680298,0.9009688679024191,0.9749279121818236,
	      -0.9781476007338057,-0.9135454576426008,-0.8090169943749473,-0.6691306063588582,-0.4999999999999998,-0.3090169943749473,-0.1045284632676533,0.1045284632676535,0.3090169943749475,0.5000000000000001,0.6691306063588582,0.8090169943749475,0.9135454576426009,0.9781476007338057,
	      -0.9807852804032304,-0.9238795325112867,-0.8314696123025453,-0.7071067811865475,-0.555570233019602,-0.3826834323650897,-0.1950903220161282,6.123031769111886e-17,0.1950903220161283,0.3826834323650898,0.5555702330196023,0.7071067811865476,0.8314696123025452,0.9238795325112867,0.9807852804032304,
	      -0.9829730996839018,-0.9324722294043558,-0.850217135729614,-0.739008917220659,-0.6026346363792563,-0.4457383557765382,-0.2736629900720829,-0.09226835946330189,0.09226835946330202,0.273662990072083,0.4457383557765383,0.6026346363792564,0.7390089172206591,0.8502171357296142,0.9324722294043558,0.9829730996839018,
	      -0.984807753012208,-0.9396926207859083,-0.8660254037844387,-0.7660444431189779,-0.6427876096865394,-0.4999999999999998,-0.3420201433256687,-0.1736481776669303,6.123031769111886e-17,0.1736481776669304,0.3420201433256688,0.5000000000000001,0.6427876096865394,0.766044443118978,0.8660254037844387,0.9396926207859084,0.984807753012208,
	      -0.9863613034027223,-0.9458172417006346,-0.879473751206489,-0.7891405093963935,-0.6772815716257409,-0.5469481581224267,-0.4016954246529694,-0.2454854871407991,-0.08257934547233227,0.08257934547233239,0.2454854871407992,0.4016954246529695,0.5469481581224269,0.6772815716257411,0.7891405093963936,0.8794737512064891,0.9458172417006346,0.9863613034027223,
	      -0.9876883405951377,-0.9510565162951535,-0.8910065241883678,-0.8090169943749473,-0.7071067811865475,-0.587785252292473,-0.4539904997395467,-0.3090169943749473,-0.1564344650402308,6.123031769111886e-17,0.1564344650402309,0.3090169943749475,0.4539904997395468,0.5877852522924731,0.7071067811865476,0.8090169943749475,0.8910065241883679,0.9510565162951535,0.9876883405951378,
	      -0.9888308262251285,-0.9555728057861408,-0.900968867902419,-0.8262387743159947,-0.7330518718298263,-0.6234898018587335,-0.4999999999999998,-0.3653410243663949,-0.2225209339563143,-0.07473009358642427,0.07473009358642439,0.2225209339563144,0.365341024366395,0.5000000000000001,0.6234898018587336,0.7330518718298263,0.8262387743159949,0.9009688679024191,0.9555728057861407,0.9888308262251285,
	      -0.9898214418809327,-0.9594929736144974,-0.9096319953545182,-0.8412535328311811,-0.7557495743542582,-0.654860733945285,-0.5406408174555977,-0.4154150130018863,-0.2817325568414297,-0.142314838273285,6.123031769111886e-17,0.1423148382732851,0.2817325568414298,0.4154150130018864,0.5406408174555977,0.6548607339452851,0.7557495743542583,0.8412535328311812,0.9096319953545184,0.9594929736144974,0.9898214418809327,
	      -0.9906859460363308,-0.9629172873477992,-0.917211301505453,-0.8544194045464883,-0.7757112907044197,-0.6825531432186539,-0.5766803221148671,-0.4600650377311521,-0.3348796121709862,-0.2034560130526336,-0.06824241336467088,0.068242413364671,0.2034560130526337,0.3348796121709863,0.4600650377311522,0.5766803221148672,0.6825531432186541,0.7757112907044198,0.8544194045464886,0.917211301505453,0.9629172873477992,0.9906859460363308,
	      -0.9914448613738104,-0.9659258262890682,-0.9238795325112867,-0.8660254037844387,-0.7933533402912351,-0.7071067811865475,-0.6087614290087207,-0.4999999999999998,-0.3826834323650897,-0.2588190451025206,-0.1305261922200516,6.123031769111886e-17,0.1305261922200517,0.2588190451025207,0.3826834323650898,0.5000000000000001,0.6087614290087207,0.7071067811865476,0.7933533402912352,0.8660254037844387,0.9238795325112867,0.9659258262890683,0.9914448613738104,
	      -0.9921147013144778,-0.9685831611286311,-0.9297764858882513,-0.8763066800438636,-0.8090169943749473,-0.7289686274214113,-0.6374239897486897,-0.5358267949789964,-0.4257792915650727,-0.3090169943749473,-0.1873813145857246,-0.0627905195293134,0.06279051952931353,0.1873813145857247,0.3090169943749475,0.4257792915650727,0.5358267949789965,0.6374239897486897,0.7289686274214116,0.8090169943749475,0.8763066800438636,0.9297764858882515,0.9685831611286311,0.9921147013144779,
	      -0.992708874098054,-0.970941817426052,-0.9350162426854147,-0.8854560256532098,-0.8229838658936564,-0.7485107481711012,-0.663122658240795,-0.5680647467311557,-0.4647231720437685,-0.3546048870425355,-0.2393156642875577,-0.1205366802553229,6.123031769111886e-17,0.120536680255323,0.2393156642875578,0.3546048870425356,0.4647231720437686,0.5680647467311559,0.6631226582407953,0.7485107481711011,0.8229838658936564,0.8854560256532099,0.9350162426854148,0.970941817426052,0.992708874098054,
	      -0.993238357741943,-0.9730448705798238,-0.9396926207859083,-0.8936326403234122,-0.8354878114129363,-0.7660444431189779,-0.6862416378687335,-0.597158591702786,-0.4999999999999998,-0.3960797660391568,-0.2868032327110902,-0.1736481776669303,-0.05814482891047577,0.0581448289104759,0.1736481776669304,0.2868032327110903,0.3960797660391569,0.5000000000000001,0.5971585917027862,0.6862416378687336,0.766044443118978,0.8354878114129365,0.8936326403234123,0.9396926207859084,0.9730448705798238,0.993238357741943,
	      -0.9937122098932426,-0.9749279121818236,-0.9438833303083676,-0.900968867902419,-0.8467241992282841,-0.7818314824680298,-0.7071067811865475,-0.6234898018587335,-0.5320320765153365,-0.4338837391175581,-0.330279061955167,-0.2225209339563143,-0.1119644761033078,6.123031769111886e-17,0.1119644761033079,0.2225209339563144,0.3302790619551672,0.4338837391175582,0.5320320765153366,0.6234898018587336,0.7071067811865476,0.7818314824680298,0.8467241992282841,0.9009688679024191,0.9438833303083676,0.9749279121818236,0.9937122098932426,
	      -0.9941379571543596,-0.9766205557100867,-0.9476531711828023,-0.9075754196709569,-0.8568571761675893,-0.7960930657056438,-0.7259954919231308,-0.6473862847818276,-0.5611870653623823,-0.46840844069979,-0.3701381553399142,-0.2675283385292206,-0.1617819965527647,-0.05413890858541748,0.05413890858541761,0.1617819965527648,0.2675283385292208,0.3701381553399143,0.4684084406997902,0.5611870653623824,0.6473862847818277,0.7259954919231308,0.7960930657056438,0.8568571761675893,0.907575419670957,0.9476531711828025,0.9766205557100867,0.9941379571543596,
	      -0.9945218953682733,-0.9781476007338057,-0.9510565162951535,-0.9135454576426008,-0.8660254037844387,-0.8090169943749473,-0.743144825477394,-0.6691306063588582,-0.587785252292473,-0.4999999999999998,-0.4067366430758,-0.3090169943749473,-0.2079116908177593,-0.1045284632676533,6.123031769111886e-17,0.1045284632676535,0.2079116908177595,0.3090169943749475,0.4067366430758002,0.5000000000000001,0.5877852522924731,0.6691306063588582,0.7431448254773942,0.8090169943749475,0.8660254037844387,0.9135454576426009,0.9510565162951535,0.9781476007338057,0.9945218953682733,
	      -0.9948693233918952,-0.9795299412524945,-0.9541392564000488,-0.9189578116202306,-0.8743466161445821,-0.8207634412072763,-0.7587581226927909,-0.6889669190756866,-0.6121059825476629,-0.5289640103269625,-0.4403941515576344,-0.3473052528448202,-0.2506525322587204,-0.1514277775045766,-0.05064916883871264,0.05064916883871277,0.1514277775045767,0.2506525322587205,0.3473052528448203,0.4403941515576343,0.5289640103269624,0.6121059825476629,0.6889669190756866,0.7587581226927909,0.8207634412072763,0.8743466161445821,0.9189578116202306,0.9541392564000488,0.9795299412524945,0.9948693233918952,
	      -0.9951847266721968,-0.9807852804032304,-0.9569403357322088,-0.9238795325112867,-0.8819212643483549,-0.8314696123025453,-0.773010453362737,-0.7071067811865475,-0.6343932841636454,-0.555570233019602,-0.4713967368259977,-0.3826834323650897,-0.2902846772544622,-0.1950903220161282,-0.09801714032956065,6.123031769111886e-17,0.09801714032956077,0.1950903220161283,0.2902846772544623,0.3826834323650898,0.4713967368259978,0.5555702330196023,0.6343932841636455,0.7071067811865476,0.773010453362737,0.8314696123025452,0.881921264348355,0.9238795325112867,0.9569403357322088,0.9807852804032304,0.9951847266721969,
	      -0.9954719225730846,-0.9819286972627066,-0.9594929736144974,-0.9283679330160726,-0.8888354486549234,-0.8412535328311811,-0.7860530947427873,-0.7237340381050702,-0.654860733945285,-0.580056909571198,-0.4999999999999998,-0.4154150130018863,-0.3270679633174217,-0.2357589355094272,-0.142314838273285,-0.04758191582374228,0.0475819158237424,0.1423148382732851,0.2357589355094273,0.3270679633174218,0.4154150130018864,0.5000000000000001,0.5800569095711982,0.6548607339452851,0.7237340381050702,0.7860530947427875,0.8412535328311812,0.8888354486549235,0.9283679330160726,0.9594929736144974,0.9819286972627067,0.9954719225730846,
	      -0.9957341762950345,-0.9829730996839018,-0.961825643172819,-0.9324722294043558,-0.8951632913550622,-0.850217135729614,-0.7980172272802395,-0.739008917220659,-0.6736956436465572,-0.6026346363792563,-0.5264321628773559,-0.4457383557765382,-0.3612416661871529,-0.2736629900720829,-0.1837495178165702,-0.09226835946330189,6.123031769111886e-17,0.09226835946330202,0.1837495178165703,0.273662990072083,0.361241666187153,0.4457383557765383,0.5264321628773558,0.6026346363792564,0.6736956436465572,0.7390089172206591,0.7980172272802396,0.8502171357296142,0.8951632913550623,0.9324722294043558,0.961825643172819,0.9829730996839018,0.9957341762950345,
	      -0.9959742939952391,-0.9839295885986297,-0.9639628606958532,-0.9362348706397372,-0.900968867902419,-0.8584487936018661,-0.8090169943749473,-0.7530714660036109,-0.6910626489868646,-0.6234898018587335,-0.5508969814521024,-0.4738686624729986,-0.3930250316539236,-0.3090169943749473,-0.2225209339563143,-0.1342332658176554,-0.04486483035051486,0.04486483035051499,0.1342332658176555,0.2225209339563144,0.3090169943749475,0.3930250316539237,0.4738686624729987,0.5508969814521025,0.6234898018587336,0.6910626489868646,0.753071466003611,0.8090169943749475,0.8584487936018661,0.9009688679024191,0.9362348706397372,0.9639628606958532,0.9839295885986297,0.9959742939952391,
	      -0.9961946980917455,-0.984807753012208,-0.9659258262890682,-0.9396926207859083,-0.9063077870366499,-0.8660254037844387,-0.8191520442889916,-0.7660444431189779,-0.7071067811865475,-0.6427876096865394,-0.5735764363510462,-0.4999999999999998,-0.4226182617406993,-0.3420201433256687,-0.2588190451025206,-0.1736481776669303,-0.08715574274765801,6.123031769111886e-17,0.08715574274765814,0.1736481776669304,0.2588190451025207,0.3420201433256688,0.4226182617406994,0.5000000000000001,0.5735764363510462,0.6427876096865394,0.7071067811865476,0.766044443118978,0.8191520442889918,0.8660254037844387,0.9063077870366499,0.9396926207859084,0.9659258262890683,0.984807753012208,0.9961946980917455,
	      -0.9963974885425265,-0.9856159103477085,-0.9677329469334988,-0.9428774454610842,-0.9112284903881356,-0.8730141131611882,-0.8285096492438421,-0.7780357543184393,-0.7219560939545244,-0.6606747233900813,-0.5946331763042866,-0.5243072835572316,-0.4502037448176734,-0.3728564777803085,-0.2928227712765503,-0.2106792699957263,-0.1270178197468788,-0.04244120319614834,0.04244120319614846,0.1270178197468789,0.2106792699957264,0.2928227712765504,0.3728564777803086,0.4502037448176733,0.5243072835572317,0.5946331763042867,0.6606747233900815,0.7219560939545245,0.7780357543184395,0.8285096492438422,0.8730141131611882,0.9112284903881357,0.9428774454610842,0.9677329469334989,0.9856159103477085,0.9963974885425265,
	      -0.9965844930066698,-0.9863613034027223,-0.9694002659393305,-0.9458172417006346,-0.9157733266550575,-0.879473751206489,-0.8371664782625285,-0.7891405093963935,-0.7357239106731316,-0.6772815716257409,-0.6142127126896678,-0.5469481581224267,-0.4759473930370736,-0.4016954246529694,-0.3246994692046835,-0.2454854871407991,-0.1645945902807338,-0.08257934547233227,6.123031769111886e-17,0.08257934547233239,0.164594590280734,0.2454854871407992,0.3246994692046836,0.4016954246529695,0.4759473930370736,0.5469481581224269,0.6142127126896678,0.6772815716257411,0.7357239106731317,0.7891405093963936,0.8371664782625287,0.8794737512064891,0.9157733266550574,0.9458172417006346,0.9694002659393304,0.9863613034027223,0.9965844930066698,
	      -0.9967573081342099,-0.9870502626379128,-0.970941817426052,-0.9485364419471455,-0.9199794436588242,-0.8854560256532098,-0.8451900855437946,-0.799442763403501,-0.7485107481711012,-0.6927243535095994,-0.6324453755953772,-0.5680647467311557,-0.4999999999999998,-0.4286925614030543,-0.3546048870425355,-0.2782174639164526,-0.2000256937760443,-0.1205366802553229,-0.04026594010941512,0.04026594010941524,0.120536680255323,0.2000256937760445,0.2782174639164527,0.3546048870425356,0.4286925614030542,0.5000000000000001,0.5680647467311559,0.6324453755953773,0.6927243535095994,0.7485107481711011,0.7994427634035012,0.8451900855437947,0.8854560256532099,0.9199794436588242,0.9485364419471455,0.970941817426052,0.9870502626379128,0.99675730813421,
	      -0.996917333733128,-0.9876883405951377,-0.9723699203976766,-0.9510565162951535,-0.9238795325112867,-0.8910065241883678,-0.8526401643540922,-0.8090169943749473,-0.7604059656000309,-0.7071067811865475,-0.6494480483301835,-0.587785252292473,-0.5224985647159488,-0.4539904997395467,-0.3826834323650897,-0.3090169943749473,-0.2334453638559053,-0.1564344650402308,-0.07845909572784487,6.123031769111886e-17,0.078459095727845,0.1564344650402309,0.2334453638559055,0.3090169943749475,0.3826834323650898,0.4539904997395468,0.5224985647159489,0.5877852522924731,0.6494480483301837,0.7071067811865476,0.7604059656000309,0.8090169943749475,0.8526401643540922,0.8910065241883679,0.9238795325112867,0.9510565162951535,0.9723699203976766,0.9876883405951378,0.996917333733128,
	      -0.9970658011837404,-0.9882804237803485,-0.973695423877779,-0.9533963920549305,-0.9275024511020946,-0.8961655569610555,-0.8595696069872012,-0.8179293607667176,-0.771489179821943,-0.7205215936007869,-0.6653257001655652,-0.6062254109666381,-0.543567550001221,-0.4777198185122627,-0.4090686371713399,-0.3380168784085026,-0.2649815021966616,-0.1903911091646683,-0.1146834253984004,-0.03830273369003537,0.03830273369003549,0.1146834253984005,0.1903911091646684,0.2649815021966617,0.3380168784085027,0.40906863717134,0.477719818512263,0.5435675500012211,0.6062254109666381,0.6653257001655654,0.720521593600787,0.771489179821943,0.8179293607667176,0.8595696069872012,0.8961655569610556,0.9275024511020947,0.9533963920549305,0.9736954238777791,0.9882804237803485,0.9970658011837404,
	      -0.9972037971811801,-0.9888308262251285,-0.9749279121818236,-0.9555728057861408,-0.9308737486442041,-0.900968867902419,-0.8660254037844387,-0.8262387743159947,-0.7818314824680298,-0.7330518718298263,-0.6801727377709192,-0.6234898018587335,-0.5633200580636221,-0.4999999999999998,-0.4338837391175581,-0.3653410243663949,-0.2947551744109042,-0.2225209339563143,-0.1490422661761743,-0.07473009358642427,6.123031769111886e-17,0.07473009358642439,0.1490422661761744,0.2225209339563144,0.2947551744109043,0.365341024366395,0.4338837391175582,0.5000000000000001,0.5633200580636221,0.6234898018587336,0.6801727377709194,0.7330518718298263,0.7818314824680298,0.8262387743159949,0.8660254037844387,0.9009688679024191,0.9308737486442042,0.9555728057861407,0.9749279121818236,0.9888308262251285,0.9972037971811801,
	      -0.9973322836635516,-0.9893433680751103,-0.9760758775559271,-0.9576005999084058,-0.9340161087325479,-0.9054482374931466,-0.8720494081438077,-0.8339978178898778,-0.791496488429254,-0.7447721827437819,-0.694074195220634,-0.6396730215588911,-0.5818589155579527,-0.5209403404879301,-0.4572423233046385,-0.3911047204901559,-0.3228804047714462,-0.2529333823916806,-0.1816368509794364,-0.1093712083778744,-0.03652202305765873,0.03652202305765885,0.1093712083778745,0.1816368509794365,0.2529333823916807,0.3228804047714463,0.391104720490156,0.4572423233046386,0.5209403404879303,0.5818589155579529,0.6396730215588913,0.694074195220634,0.7447721827437819,0.7914964884292541,0.8339978178898779,0.8720494081438076,0.9054482374931466,0.934016108732548,0.957600599908406,0.9760758775559272,0.9893433680751103,0.9973322836635516,
	      -0.9974521146102535,-0.9898214418809327,-0.9771468659711595,-0.9594929736144974,-0.9369497249997618,-0.9096319953545182,-0.8776789895672555,-0.8412535328311811,-0.8005412409243603,-0.7557495743542582,-0.7071067811865475,-0.654860733945285,-0.599277666511347,-0.5406408174555977,-0.4792489867200569,-0.4154150130018863,-0.3494641795990983,-0.2817325568414297,-0.2125652895529767,-0.142314838273285,-0.07133918319923224,6.123031769111886e-17,0.07133918319923235,0.1423148382732851,0.2125652895529768,0.2817325568414298,0.3494641795990984,0.4154150130018864,0.4792489867200568,0.5406408174555977,0.599277666511347,0.6548607339452851,0.7071067811865476,0.7557495743542583,0.8005412409243604,0.8412535328311812,0.8776789895672557,0.9096319953545184,0.9369497249997617,0.9594929736144974,0.9771468659711595,0.9898214418809327,0.9974521146102535,
	      -0.9975640502598242,-0.9902680687415703,-0.9781476007338057,-0.9612616959383189,-0.9396926207859083,-0.9135454576426008,-0.882947592858927,-0.848048096156426,-0.8090169943749473,-0.7660444431189779,-0.719339800338651,-0.6691306063588582,-0.6156614753256583,-0.5591929034707467,-0.4999999999999998,-0.4383711467890775,-0.3746065934159121,-0.3090169943749473,-0.2419218955996676,-0.1736481776669303,-0.1045284632676533,-0.03489949670250096,0.03489949670250108,0.1045284632676535,0.1736481776669304,0.2419218955996677,0.3090169943749475,0.3746065934159122,0.4383711467890775,0.5000000000000001,0.5591929034707468,0.6156614753256583,0.6691306063588582,0.7193398003386512,0.766044443118978,0.8090169943749475,0.848048096156426,0.882947592858927,0.9135454576426009,0.9396926207859084,0.9612616959383189,0.9781476007338057,0.9902680687415704,0.9975640502598242,
	      -0.9976687691905392,-0.9906859460363308,-0.9790840876823228,-0.9629172873477992,-0.9422609221188204,-0.917211301505453,-0.8878852184023752,-0.8544194045464883,-0.816969893010442,-0.7757112907044197,-0.7308359642781241,-0.6825531432186539,-0.6310879443260529,-0.5766803221148671,-0.5195839500354337,-0.4600650377311521,-0.3984010898462415,-0.3348796121709862,-0.2697967711570243,-0.2034560130526336,-0.1361666490962465,-0.06824241336467088,6.123031769111886e-17,0.068242413364671,0.1361666490962466,0.2034560130526337,0.2697967711570244,0.3348796121709863,0.3984010898462416,0.4600650377311522,0.5195839500354336,0.5766803221148672,0.6310879443260529,0.6825531432186541,0.7308359642781241,0.7757112907044198,0.8169698930104421,0.8544194045464886,0.8878852184023752,0.917211301505453,0.9422609221188205,0.9629172873477992,0.9790840876823229,0.9906859460363308,0.9976687691905392,
	      -0.9977668786231532,-0.99107748815478,-0.9799617050365867,-0.9644691750543765,-0.9446690916079187,-0.9206498866764287,-0.8925188358598811,-0.8604015792601394,-0.8244415603417604,-0.784799385278661,-0.7416521056479576,-0.6951924276746423,-0.6456278515588024,-0.5931797447293552,-0.5380823531633726,-0.4805817551866837,-0.4209347624283349,-0.3594077728375127,-0.2962755808856338,-0.2318201502675283,-0.16632935458313,-0.1000956916240983,-0.03341497700767452,0.03341497700767464,0.1000956916240984,0.1663293545831302,0.2318201502675284,0.2962755808856339,0.3594077728375128,0.420934762428335,0.4805817551866838,0.5380823531633727,0.5931797447293553,0.6456278515588024,0.6951924276746423,0.7416521056479576,0.784799385278661,0.8244415603417603,0.8604015792601394,0.8925188358598812,0.9206498866764288,0.9446690916079188,0.9644691750543766,0.9799617050365869,0.9910774881547801,0.9977668786231532,
	      -0.9978589232386035,-0.9914448613738104,-0.9807852804032304,-0.9659258262890682,-0.9469301294951056,-0.9238795325112867,-0.8968727415326883,-0.8660254037844387,-0.8314696123025453,-0.7933533402912351,-0.7518398074789773,-0.7071067811865475,-0.6593458151000688,-0.6087614290087207,-0.555570233019602,-0.4999999999999998,-0.4422886902190011,-0.3826834323650897,-0.3214394653031616,-0.2588190451025206,-0.1950903220161282,-0.1305261922200516,-0.06540312923014292,6.123031769111886e-17,0.06540312923014305,0.1305261922200517,0.1950903220161283,0.2588190451025207,0.3214394653031617,0.3826834323650898,0.4422886902190012,0.5000000000000001,0.5555702330196023,0.6087614290087207,0.6593458151000688,0.7071067811865476,0.7518398074789774,0.7933533402912352,0.8314696123025452,0.8660254037844387,0.8968727415326884,0.9238795325112867,0.9469301294951057,0.9659258262890683,0.9807852804032304,0.9914448613738104,0.9978589232386035,
	      -0.9979453927503363,-0.9917900138232461,-0.9815591569910653,-0.9672948630390295,-0.9490557470106686,-0.9269167573460217,-0.900968867902419,-0.8713187041233892,-0.8380881048918406,-0.8014136218679565,-0.7614459583691344,-0.7183493500977275,-0.6723008902613169,-0.6234898018587335,-0.5721166601221694,-0.518392568310525,-0.4625382902408351,-0.4047833431223937,-0.3453650544213075,-0.2845275866310323,-0.2225209339563143,-0.1595998950333792,-0.09602302590768176,-0.03205157757165521,0.03205157757165533,0.09602302590768189,0.1595998950333793,0.2225209339563144,0.2845275866310324,0.3453650544213076,0.4047833431223938,0.4625382902408354,0.5183925683105252,0.5721166601221697,0.6234898018587336,0.6723008902613168,0.7183493500977276,0.7614459583691344,0.8014136218679566,0.8380881048918406,0.8713187041233894,0.9009688679024191,0.9269167573460217,0.9490557470106686,0.9672948630390295,0.9815591569910653,0.9917900138232462,0.9979453927503363,
	      -0.9980267284282716,-0.9921147013144778,-0.9822872507286887,-0.9685831611286311,-0.9510565162951535,-0.9297764858882513,-0.9048270524660194,-0.8763066800438636,-0.8443279255020151,-0.8090169943749473,-0.7705132427757891,-0.7289686274214113,-0.6845471059286887,-0.6374239897486897,-0.587785252292473,-0.5358267949789964,-0.481753674101715,-0.4257792915650727,-0.368124552684678,-0.3090169943749473,-0.2486898871648546,-0.1873813145857246,-0.1253332335643041,-0.0627905195293134,6.123031769111886e-17,0.06279051952931353,0.1253332335643043,0.1873813145857247,0.2486898871648547,0.3090169943749475,0.3681245526846781,0.4257792915650727,0.4817536741017153,0.5358267949789965,0.5877852522924731,0.6374239897486897,0.6845471059286886,0.7289686274214116,0.7705132427757893,0.8090169943749475,0.8443279255020151,0.8763066800438636,0.9048270524660196,0.9297764858882515,0.9510565162951535,0.9685831611286311,0.9822872507286887,0.9921147013144779,0.9980267284282716,
	      -0.9981033287370441,-0.9924205096719357,-0.9829730996839018,-0.9697969360350094,-0.9529420004271565,-0.9324722294043558,-0.9084652718195236,-0.8810121942857845,-0.850217135729614,-0.8161969123562216,-0.7790805745256705,-0.739008917220659,-0.6961339459629265,-0.6506183002042422,-0.6026346363792563,-0.5523649729605058,-0.4999999999999998,-0.4457383557765382,-0.3897858732926793,-0.3323547994796595,-0.2736629900720829,-0.2139330832064974,-0.1533916548786853,-0.09226835946330189,-0.0307950585561702,0.03079505855617033,0.09226835946330202,0.1533916548786854,0.2139330832064975,0.273662990072083,0.3323547994796596,0.3897858732926794,0.4457383557765383,0.5000000000000001,0.5523649729605059,0.6026346363792564,0.6506183002042422,0.6961339459629267,0.7390089172206591,0.7790805745256705,0.8161969123562217,0.8502171357296142,0.8810121942857845,0.9084652718195237,0.9324722294043558,0.9529420004271566,0.9697969360350095,0.9829730996839018,0.9924205096719357,0.9981033287370441,
	      -0.9988322268323265,-0.9953316347176486,-0.9895063994510511,-0.9813701261394134,-0.970941817426052,-0.9582458291091662,-0.943311813257743,-0.9261746489577763,-0.9068743608505454,-0.8854560256532098,-0.8619696668800491,-0.8364701380102267,-0.8090169943749473,-0.7796743540632223,-0.7485107481711012,-0.715598960744121,-0.681015858786797,-0.6448422127361707,-0.6071625078187112,-0.5680647467311557,-0.5276402441061325,-0.4859834132426061,-0.4431915455992412,-0.3993645835656955,-0.3546048870425355,-0.3090169943749473,-0.2627073781985869,-0.215784196767806,-0.1683570413470385,-0.1205366802553229,-0.07243480016176228,-0.02416374523613229,0.02416374523613242,0.0724348001617624,0.120536680255323,0.1683570413470386,0.2157841967678061,0.2627073781985871,0.3090169943749475,0.3546048870425356,0.3993645835656957,0.4431915455992413,0.4859834132426062,0.5276402441061327,0.5680647467311559,0.6071625078187112,0.6448422127361706,0.6810158587867972,0.7155989607441211,0.7485107481711011,0.7796743540632224,0.8090169943749475,0.8364701380102267,0.8619696668800493,0.8854560256532099,0.9068743608505454,0.9261746489577766,0.943311813257743,0.9582458291091662,0.970941817426052,0.9813701261394134,0.9895063994510511,0.9953316347176486,0.9988322268323266,
	      -0.9997034698451394,-0.998814055240823,-0.9973322836635516,-0.9952590338932358,-0.9925955354920264,-0.9893433680751103,-0.9855044603739027,-0.9810810890921943,-0.9760758775559271,-0.9704917941574054,-0.9643321505948586,-0.9576005999084058,-0.9503011343135823,-0.9424380828337144,-0.9340161087325479,-0.9250402067486519,-0.915515700133237,-0.9054482374931466,-0.8948437894408918,-0.883708645053719,-0.8720494081438077,-0.8598729933418101,-0.8471866219960603,-0.8339978178898778,-0.820314402779511,-0.8061444917553627,-0.791496488429254,-0.7763790799505742,-0.7608012318542781,-0.7447721827437819,-0.7283014388119159,-0.7113987682031776,-0.694074195220634,-0.6763379943809028,-0.6582006843207479,-0.6396730215588911,-0.6207659941167485,-0.6014908150018703,-0.5818589155579527,-0.5618819386853605,-0.5415717319361846,-0.5209403404879301,-0.4999999999999998,-0.4787631293572092,-0.4572423233046385,-0.4354503449781911,-0.413400118335283,-0.3911047204901559,-0.3685773739583617,-0.3458314388150113,-0.3228804047714462,-0.2997378831750241,-0.2764175989367714,-0.2529333823916806,-0.2292991610964899,-0.2055289515698006,-0.1816368509794364,-0.1576370287819737,-0.1335437183193978,-0.1093712083778744,-0.08513383471363556,-0.06084597155101392,-0.03652202305765873,-0.01217641480199753,0.01217641480199765,0.03652202305765885,0.06084597155101405,0.08513383471363568,0.1093712083778745,0.133543718319398,0.1576370287819738,0.1816368509794365,0.2055289515698007,0.22929916109649,0.2529333823916807,0.2764175989367715,0.2997378831750243,0.3228804047714463,0.3458314388150114,0.3685773739583618,0.391104720490156,0.4134001183352831,0.4354503449781912,0.4572423233046386,0.4787631293572091,0.5000000000000001,0.5209403404879303,0.5415717319361846,0.5618819386853604,0.5818589155579529,0.6014908150018704,0.6207659941167485,0.6396730215588913,0.6582006843207481,0.6763379943809029,0.694074195220634,0.7113987682031779,0.728301438811916,0.7447721827437819,0.7608012318542781,0.7763790799505744,0.7914964884292541,0.8061444917553627,0.820314402779511,0.8339978178898779,0.8471866219960603,0.8598729933418101,0.8720494081438076,0.8837086450537192,0.8948437894408918,0.9054482374931466,0.9155157001332371,0.9250402067486521,0.934016108732548,0.9424380828337144,0.9503011343135824,0.957600599908406,0.9643321505948586,0.9704917941574053,0.9760758775559272,0.9810810890921943,0.9855044603739027,0.9893433680751103,0.9925955354920264,0.9952590338932358,0.9973322836635516,0.998814055240823,0.9997034698451394,
	      -0.9999252866697326,-0.9997011578430937,-0.9993276470109054,-0.9988048099856441,-0.9981327248931005,-0.9973114921607055,-0.9963412345025239,-0.9952220969009173,-0.9939542465848804,-0.9925378730050517,-0.9909731878054056,-0.9892604247916259,-0.9873998398961705,-0.9853917111400267,-0.9832363385911683,-0.9809340443197178,-0.9784851723498196,-0.9758900886082343,-0.9731491808696592,-0.9702628586987844,-0.9672315533890932,-0.964055717898415,-0.960735826781242,-0.9572723761178167,-0.9536658834400057,-0.9499168876539662,-0.9460259489596189,-0.9419936487669395,-0.9378205896090794,-0.9335073950523333,-0.9290547096029599,-0.9244631986108762,-0.919733548170237,-0.9148664650169126,-0.9098626764228853,-0.9047229300875749,-0.8994479940261121,-0.8940386564545773,-0.8884957256722192,-0.8828200299406745,-0.8770124173602024,-0.8710737557429565,-0.8650049324833107,-0.8588068544252573,-0.8524804477269028,-0.8460266577220731,-0.8394464487790565,-0.8327408041565012,-0.8259107258564895,-0.8189572344748131,-0.8118813690484682,-0.8046841869003959,-0.79736676348149,-0.7899301922098955,-0.7823755843076232,-0.7747040686345038,-0.7669167915195056,-0.7590149165894429,-0.7509996245950976,-0.7428721132347862,-0.7346335969753897,-0.7262853068708802,-0.7178284903783699,-0.7092644111717052,-0.7005943489526429,-0.691819599259627,-0.6829414732742012,-0.6739612976250832,-0.6648804141899316,-0.6557001798948323,-0.6464219665115392,-0.6370471604524938,-0.6275771625636587,-0.6180133879151942,-0.6083572655900072,-0.5986102384702111,-0.5887737630215184,-0.5788493090756068,-0.5688383596104875,-0.5587424105289069,-0.5485629704348219,-0.5383015604079718,-0.5279597137765891,-0.5175389758882802,-0.5070409038791074,-0.4964670664409133,-0.4858190435869153,-0.4750984264156099,-0.4643068168730203,-0.4534458275133212,-0.4425170812578835,-0.4315222111527637,-0.4204628601246847,-0.4093406807355385,-0.398157334935449,-0.3869144938144324,-0.3756138373526925,-0.3642570541695865,-0.3528458412712995,-0.3413819037972679,-0.3298669547653846,-0.31830271481603,-0.3066909119549613,-0.2950332812951027,-0.2833315647972737,-0.2715875110098926,-0.2598028748076984,-0.2479794171295246,-0.2361189047151674,-0.2242231098413897,-0.2122938100570935,-0.2003327879177082,-0.188341830718829,-0.1763227302291461,-0.164277282422709,-0.1522072872105575,-0.140114548171769,-0.128000872283955,-0.11586806965325,-0.1037179532438342,-0.09155233860702745,-0.0793730436099985,-0.06718188816412665,-0.05498069395305778,-0.04277128416049784,-0.03055548319777877,-0.01833511643124444,-0.006112009909491491,0.006112009909491613,0.01833511643124456,0.0305554831977789,0.04277128416049795,0.0549806939530579,0.06718188816412678,0.07937304360999863,0.09155233860702758,0.1037179532438343,0.1158680696532501,0.1280008722839552,0.1401145481717691,0.1522072872105576,0.1642772824227091,0.1763227302291462,0.1883418307188291,0.2003327879177083,0.2122938100570936,0.2242231098413899,0.2361189047151675,0.2479794171295247,0.2598028748076985,0.2715875110098928,0.2833315647972738,0.2950332812951028,0.3066909119549615,0.3183027148160301,0.3298669547653847,0.3413819037972681,0.3528458412712996,0.3642570541695866,0.3756138373526926,0.3869144938144325,0.3981573349354491,0.4093406807355386,0.420462860124685,0.431522211152764,0.4425170812578836,0.4534458275133214,0.4643068168730202,0.4750984264156102,0.4858190435869156,0.4964670664409134,0.5070409038791075,0.5175389758882802,0.5279597137765895,0.5383015604079721,0.5485629704348221,0.5587424105289071,0.5688383596104875,0.5788493090756071,0.5887737630215186,0.5986102384702112,0.6083572655900074,0.6180133879151942,0.627577162563659,0.6370471604524939,0.6464219665115392,0.6557001798948323,0.6648804141899316,0.6739612976250834,0.6829414732742013,0.6918195992596271,0.7005943489526431,0.7092644111717052,0.7178284903783698,0.7262853068708804,0.7346335969753897,0.7428721132347863,0.7509996245950976,0.7590149165894428,0.7669167915195058,0.7747040686345039,0.7823755843076232,0.7899301922098955,0.7973667634814899,0.804684186900396,0.8118813690484682,0.8189572344748132,0.8259107258564895,0.8327408041565012,0.8394464487790566,0.8460266577220732,0.8524804477269028,0.8588068544252573,0.8650049324833106,0.8710737557429566,0.8770124173602025,0.8828200299406745,0.8884957256722192,0.8940386564545773,0.8994479940261122,0.9047229300875749,0.9098626764228853,0.9148664650169126,0.919733548170237,0.9244631986108763,0.9290547096029599,0.9335073950523334,0.9378205896090794,0.9419936487669394,0.946025948959619,0.9499168876539663,0.9536658834400057,0.9572723761178167,0.960735826781242,0.9640557178984152,0.9672315533890933,0.9702628586987844,0.9731491808696592,0.9758900886082343,0.9784851723498197,0.9809340443197179,0.9832363385911683,0.9853917111400267,0.9873998398961705,0.989260424791626,0.9909731878054056,0.9925378730050517,0.9939542465848804,0.9952220969009173,0.9963412345025239,0.9973114921607055,0.9981327248931006,0.9988048099856441,0.9993276470109054,0.9997011578430937,0.9999252866697326,
	    };

	  unsigned order_array[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,64,128,256,};
	  unsigned order_array_length = 53;

	  // Construct the weights or knots as requested
	  if(weights)
	    construction_helper(weights_data_array, order_array, order_array_length);
	  else
	    construction_helper(knots_data_array, order_array, order_array_length);

	  break;
	}


      default :
	throw OomphLibError("The given quadrature scheme does not exist.",
			    "weights_data_structure::weights_data_structure",
			    "OOMPHEXCEPTIONLOCATION");
      }
  }


  // Call the constructors.
  const weights_data_structure VariableOrderGaussLegendre::Weights(0,1);
  const weights_data_structure VariableOrderGaussLegendre::Knots(0,0);
  const weights_data_structure VariableOrderClenshawCurtis::Weights(1,1);
  const weights_data_structure VariableOrderClenshawCurtis::Knots(1,0);
  const weights_data_structure VariableOrderFejerSecond::Weights(2,1);
  const weights_data_structure VariableOrderFejerSecond::Knots(2,0);
}

#endif
