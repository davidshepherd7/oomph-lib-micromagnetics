
/* ??ds c++0x c++11 C++0x C++11
   This header will only compile when using the 2011 version of C++.
   We have used the new methods for initialising vectors to lists.
   To generalise to earlier versions the vector initialisers for each
   weight and knot must be replaced with the old versions, i.e.:

   Weights[0] = {1.0,1.0};
   ---->>
   a[] = {1.0,1.0};
   Weights[0](a,a+1);

   Hopefully by the time this code is being used this will not be an issue.
*/

#include "generic.h"

using namespace oomph;

namespace oomph
{

  //============================================================
  /// Class to hold all the weights for Gaussian quadrature
  /// (required to have a static list of weights).
  //============================================================
  template<unsigned DIM>
  class GaussWeights
  {
  public:
    /// Default construtor - fills in all the data
    GaussWeights();

    /// Broken copy constructor
    GaussWeights(const GaussWeights& dummy);
    //    {BrokenCopy::broken_copy("GaussWeights");}

    /// Broken assignment operator
    void operator=(const GaussWeights& dummy)
    {BrokenCopy::broken_assign("GaussWeights");}

    /// Destructor
    ~GaussWeights(){};

    /// Get function
    double get_weight(const unsigned &order, const unsigned &i) const
    {
      //??ds #ifdef PARANOID check that order >=2, i ok for given DIM
      return Weights[order-2][i];
    }

  private:
    std::vector<std::vector<double> > Weights;

  };

  //============================================================
  /// Class to hold all the knots for Gaussian quadrature
  /// (required to have a static list of knots).
  //============================================================
  template<unsigned DIM>
  class GaussKnots
  {
  public:
    /// Default construtor - fills in all the data
    GaussKnots();

    /// Broken copy constructor
    GaussKnots(const GaussKnots& dummy)
    {BrokenCopy::broken_copy("GaussKnots");}

    /// Broken assignment operator
    void operator=(const GaussKnots& dummy)
    {BrokenCopy::broken_assign("GaussKnots");}

    /// Destructor
    ~GaussKnots(){};

    /// Get function
    double get_knot(const unsigned &order, const unsigned &i, const unsigned &j) const
    {
      //??ds #ifdef PARANOID check that order >=2, i ok for given DIM
      return Knots[order-2][i][j];
    }

  private:
    std::vector<std::vector<std::vector<double> > > Knots;

  };

  //=========================================================
  /// Variable order Gaussian quadrature. The idea is to allow the
  /// order of the quadrature to be decided at run time - to allow the use of
  /// predictive schemes etc.
  /// ??ds for now this is only implemented for 1D and order
  /// from 2 to 50.
  //=========================================================
  template<unsigned DIM>
  class VariableGauss
  {
  private:

    /// A static data structure holding all the Gaussian quadrature weights
    // Another class is needed here to allow us to use static vectors.
    static const GaussWeights<DIM> Weights;

    /// A static data structure holding all the Gaussian quadrature knots
    static const GaussKnots<DIM> Knots;

  public:

    /// Default construtor
    VariableGauss(){};

    /// Broken copy constructor
    VariableGauss(const VariableGauss& dummy)
    {BrokenCopy::broken_copy("VariableGauss");}

    /// Broken assignment operator
    void operator=(const VariableGauss& dummy)
    {BrokenCopy::broken_assign("VariableGauss");}

    /// Destructor
    ~VariableGauss(){};

    /// Return the number of integration points of the order-th Gauss scheme.
    unsigned nweight(const unsigned &order) const
    {return order;}

    /// Return weight of i-th integration point in the order-th Gauss scheme.
    double weight(const unsigned &order, const unsigned &i) const
    {return Weights.get_weight(order,i);}

    /// Return local coordinate s[j] of i-th integration point in the order-th Gauss scheme.
    double knot(const unsigned &order, const unsigned &i, const unsigned &j) const
    {
      return Knots.get_knot(order,i,j);
    }

  };

}
