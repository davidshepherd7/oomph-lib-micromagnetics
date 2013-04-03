#ifndef OOMPH_POLY_INTERP_H
#define OOMPH_POLY_INTERP_H


#include "generic.h"
#include "./vector_helpers.h"

using namespace oomph;

namespace oomph
{

  /// Check if two doubles are with a tol of each other
  bool almost_equal(double a, double b, double tol=1e-10)
  {
    return std::abs(a - b) < tol;
  }

  void double_vector_add_ax(double a, const DoubleVector& x,
                             DoubleVector& y)
  {
    OOMPH_ASSERT(y.built());
    OOMPH_ASSERT(x.built());
    OOMPH_ASSERT(*x.distribution_pt() == *y.distribution_pt());

    // cache (to make sure there's no function call overhead)
   double* x_values_pt = x.values_pt();
   double* y_values_pt = y.values_pt();

   for (unsigned i = 0, ni=y.nrow_local(); i<ni; i++)
    {
      y_values_pt[i] += a * x_values_pt[i];
    }
  }

  bool almost_equal(const DoubleVector& x,
                    const DoubleVector &y,
                    const double &tol=1e-10)
  {
    if(*x.distribution_pt() != *y.distribution_pt())
      {
        return false;
      }
    else
      {
        for(unsigned i=0, ni=x.nrow_local(); i<ni; i++)
          {
            if(!almost_equal(x[i], y[i], tol))
              {
                return false;
              }
          }
      }
    return true;
  }


  // ============================================================
  ///
  // ============================================================
  class PolynomialInterpolatorBase
  {
  public:
    /// Get value at point
    virtual void eval(const double& x, DoubleVector& result) const =0;

    /// Get nth derivative at point
    virtual void eval_derivative(const double& x, const unsigned &deriv_order,
                                 DoubleVector& result) const =0;

  };


  // ============================================================
  ///
  // ============================================================
  class BarycentricLagrangeInterpolator : public PolynomialInterpolatorBase
  {
  public:

    /// Construct from locations, values lists.
    BarycentricLagrangeInterpolator(const Vector<double>& locations,
                                    const Vector<DoubleVector >& values)
      : Locations(locations), Values(values)
    { build(); }

    /// Get value at point
    void eval(const double& x, DoubleVector& result) const;

    /// Get nth derivative at point
    void eval_derivative(const double& x, const unsigned &deriv_order,
                         DoubleVector& result) const;

  private:

    /// Construct the weights
    void build()
    {
      assert(!VectorOps::contains_duplicates(Locations));
      assert(Locations.size() == Values.size());

      // Make the right sized weights vector
      Weights.clear();
      Weights.resize(Locations.size());

      // Calculate each weight according to (3.2) in Berrut2004 (O(N^2) where
      // N is the number of interpolation points).
      for(unsigned j=0, nj=Locations.size(); j<nj; j++)
        {
          double product = 1.0;
          for(unsigned k=0, nk = Locations.size(); k<nk; k++)
            {
              if(k != j) product *= (Locations[j] - Locations[k]);
            }

          Weights[j] = 1 / product;
        }
    }

    void eval_checks(const double &x) const;

    /// Data storage
    Vector<double> Locations;
    Vector<DoubleVector> Values;
    Vector<double> Weights;

  };


  /// Get value at point
  void BarycentricLagrangeInterpolator::
  eval(const double& x, DoubleVector& result) const
  {
    eval_checks(x);

    // Calculate the polynomial according to (4.2) of Berrut2004.
    DoubleVector numerator(Values[0]);
    numerator.initialise(0.0);
    double denominator = 0;
    for(unsigned i=0, ni=Locations.size(); i<ni; i++)
      {
        double temp = Weights[i]  / (x - Locations[i]);
        double_vector_add_ax(temp, Values[i], numerator);
        denominator += temp;
      }

    result = numerator;
    result /= denominator;
  }


  void BarycentricLagrangeInterpolator::
  eval_derivative(const double& x, const unsigned &deriv_order,
                  DoubleVector& result) const
  {
    eval_checks(x);

    // Only implemented for first derivatives...
    OOMPH_ASSERT(deriv_order == 1);

    DoubleVector g, dg;
    g.build(Values[0].distribution_pt(), 0.0);
    dg.build(Values[0].distribution_pt(), 0.0);
    double h = 0, dh = 0;

    for(unsigned i=0, ni=Locations.size(); i<ni; i++)
      {
        double temp = Weights[i]  / (x - Locations[i]);

        double_vector_add_ax(temp, Values[i], g);
        double_vector_add_ax(-1* temp / (x - Locations[i]), Values[i], dg);
        h += temp;
        dh += -1 * temp / (x - Locations[i]);
      }

    OOMPH_ASSERT(h != 0);
    result.build(g.distribution_pt(), 0.0);
    for(unsigned i=0, ni=g.nrow_local(); i<ni; i++)
      {
        result[i] = (dg[i] * h - g[i]*dh) / (h*h);
      }
  }


  void BarycentricLagrangeInterpolator::
  eval_checks(const double &x) const
  {
    // Check the everything has been set up right
    OOMPH_ASSERT(Locations.size() != 0);
    OOMPH_ASSERT(Locations.size() == Values.size());
    OOMPH_ASSERT(Weights.size() == Locations.size());

    // Check that x is not a given location (trivial calculation gives
    // undefined value here, should be possible to output the appropriate
    // value if this functionality is needed...).
    OOMPH_ASSERT(std::find(Locations.begin(), Locations.end(), x) == Locations.end());
  }

} // End of oomph namespace

#endif
