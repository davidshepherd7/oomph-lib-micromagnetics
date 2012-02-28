#ifndef OOMPH_VARIABLE_QUADRATURE_H
#define OOMPH_VARIABLE_QUADRATURE_H

#include "generic.h"
#include <cmath>
/* Weights generated using C++ QUADRULE by John Burkardt.

   You may want to disable line wrapping for this file as the lines
   containing the data are very long. In emacs this is done using M-x
   toggle-truncate-lines.
*/

using namespace oomph;

namespace oomph
{


  //============================================================
  /// A structure to hold the weights/knots data. We need this so we can have a
  /// constructor to initialise the const static data structures.
  //============================================================
  //??ds Had to use std::vector because our Vector class lacks some constructors.
  //??ds Ask if I can add them?
  struct weights_data_structure
  {

  private:

    /// Typedef to keep code lines reasonable and allow easy swapping between
    /// std::vector and oomph::Vector if needed.
    typedef std::vector<double> vdb;

    /// The data structure itself
    std::map<unsigned,vdb > data;

  public:

    /// Constructor - to allow construction of static constants
    weights_data_structure(const double data_array[],
			   const unsigned order_array[],
			   const unsigned &order_array_length)
    {
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
      unsigned i_y = i%order;
      unsigned i_x = i - i_y*order;
      return weight_1d(i_y,order)*weight_1d(i_x,order);
    }

    inline double knot(const unsigned &i, const unsigned &j, const unsigned &order) const
    {
#ifdef PARANOID
      error_check(i,j,order,"VariableOrderGaussLegendre::knot");
#endif
      unsigned i_y = i%order;
      if(j==1) return knot_1d(i_y,order);
      else if(j==0)
	{
	  unsigned i_x = i - i_y*order;
	  return knot_1d(i_x,order);
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
      unsigned i_z = i%(order*order);
      unsigned i_y = (i-i_z*order*order)%order;
      unsigned i_x = i - i_z*order*order - i_y*order;
      return weight_1d(i_z,order)*weight_1d(i_y,order)*weight_1d(i_x,order);
    }

    inline double knot(const unsigned &i, const unsigned &j, const unsigned &order) const
    {
#ifdef PARANOID
      error_check(i,j,order,"VariableOrderGaussLegendre::knot");
#endif
      unsigned i_z = i%(order*order);
      if(j==2) return knot_1d(i_z,order);

      unsigned i_y = (i-i_z*order*order)%order;
      if(j==1)
	{
	  return knot_1d(i_y,order);
	}
      else if(j==0)
	{
	  unsigned i_x = i - i_z*order*order - i_y*order;
	  return knot_1d(i_x,order);
	}
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

    /// Get the next order to use in an adaptive scheme.
    inline unsigned adaptive_scheme_next_order(const unsigned &order) const
    {
      if(order == 0) return 2;
      else return 2*order;
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

    /// Get the next highest order allowing reuse of all previous knots
    inline unsigned adaptive_scheme_next_order(const unsigned &order) const
    {
      if(order ==0) return 2;
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

    /// Get the next highest order allowing reuse of previous knots
    inline unsigned adaptive_scheme_next_order(const unsigned &order) const
    {
      //??ds maybe should have an error check?
      if(order==0) return 2;
      else return (2*order) + 1;
    }
  };


  //============================================================
  /// The final class for GaussLegendre on a QElement
  //============================================================
  template <unsigned DIM>
  class QVariableOrderGaussLegendre : public VariableOrderGaussLegendre,
				      public QVariableOrderQuadrature<DIM>
  {
    // Just make sure we are calling the right functions
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

  // void VariableOrderGaussLegendre::intialise_gauss_legendre_weights()
  // {
  //   double weights_data_array[] =
  //     {
  // 	2,
  // 	0.9999999999999998,0.9999999999999998,
  // 	0.5555555555555556,0.8888888888888888,0.5555555555555556,
  //     };

  //   double knots_data_array[] =
  //     {
  // 	0,
  // 	-0.5773502691896257,0.5773502691896257,
  // 	-0.7745966692414834,0,0.7745966692414834,
  //     };

  //   unsigned order_array[] = {1,2,3,};
  //   unsigned order_array_length = 3;

  //   const weights_data_structure
  //     Weights(weights_data_array,order_array,order_array_length);

  //   const weights_data_structure
  //     Knots(knots_data_array,order_array,order_array_length);
  // }

  // VariableOrderGaussLegendre<1>::intialise_gauss_legendre_weights()



  // double weights_data_arrray[] =
  //   {
  //     2,
  //     1,1,
  //     0.3333333333333333,1.333333333333333,0.3333333333333333,
  //   };

  // double knots_data_arrray[] =
  //   {
  //     0,
  //     1,-1,
  //     1,6.123031769111886e-17,-1,
  //   };

  // unsigned order_array[] = {1,2,3,};
  // unsigned order_array_length = 3;
  // const weights_data_structure VariableOrderClenshawCurtis<DIM>::
  //   Weights(weights_data_array,order_array,order_array_length);
  // const weights_data_structure VariableOrderClenshawCurtis<DIM>::
  //   Knots(knots_data_array,order_array,order_array_length);


  // // Enclose in braces to keep all these temporary arrays contained.

  // double weights_data_arrray[] =
  //   {
  //     2,
  //     0.6666666666666666,0.6666666666666666,
  //     0.6666666666666667,0.6666666666666667,0.6666666666666666,
  //   };

  // double knots_data_arrray[] =
  //   {
  //     0,
  //     -0.4999999999999998,0.5000000000000001,
  //     -0.7071067811865475,6.123031769111886e-17,0.7071067811865476,
  //   };

  // unsigned order_array[] = {1,2,3,};
  // unsigned order_array_length = 3;
  // const weights_data_structure VariableOrderFejerSecond<DIM>::
  //   Weights(weights_data_array,order_array,order_array_length);
  // const weights_data_structure VariableOrderFejerSecond<DIM>::
  //   Knots(knots_data_array,order_array,order_array_length);



}

#endif
