#ifndef OOMPH_POLY_INTERP_H
#define OOMPH_POLY_INTERP_H


#include "generic.h"
#include "./vector_helpers.h"
#include "./my_assert.h"

using namespace oomph;

namespace oomph
{

  /// Check if two doubles are with a tol of each other
  bool almost_equal(double a, double b, double tol=1e-10)
  {
    return std::abs(a - b) < tol;
  }

  void vector_add_ax(double a, const Vector<double>& x,
                             Vector<double>& y)
  {
    my_assert(y.size() == x.size());
    for (unsigned i = 0, ni=y.size(); i<ni; i++)
      {
        y[i] += a * x[i];
      }
  }

  bool almost_equal(const Vector<double>& x,
                    const Vector<double> &y,
                    const double &tol=1e-10)
  {
    if(x.size() != y.size())
      {
        return false;
      }
    else
      {
        for(unsigned i=0, ni=x.size(); i<ni; i++)
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
    virtual void eval(const double& x, Vector<double>& result) const =0;

    /// Get nth derivative at point
    virtual void eval_derivative(const double& x, const unsigned &deriv_order,
                                 Vector<double>& result) const =0;

  };


  // ============================================================
  ///
  // ============================================================
  class BarycentricLagrangeInterpolator : public PolynomialInterpolatorBase
  {
  public:

    /// Construct from locations, values lists.
    BarycentricLagrangeInterpolator(const Vector<double>& locations,
                                    const Vector<Vector<double> >& values)
      : Locations(locations), Values(values)
    { build(); }

    /// Get value at point
    void eval(const double& x, Vector<double>& result) const;

    /// Get nth derivative at point
    void eval_derivative(const double& x, const unsigned &deriv_order,
                         Vector<double>& result) const;

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
    Vector<Vector<double> > Values;
    Vector<double> Weights;

  };


  /// Get value at point
  void BarycentricLagrangeInterpolator::
  eval(const double& x, Vector<double>& result) const
  {
    eval_checks(x);

    // Calculate the polynomial according to (4.2) of Berrut2004.
    Vector<double> numerator(Values[0]);
    numerator.initialise(0.0);
    double denominator = 0;
    for(unsigned i=0, ni=Locations.size(); i<ni; i++)
      {
        double temp = Weights[i]  / (x - Locations[i]);
        vector_add_ax(temp, Values[i], numerator);
        denominator += temp;
      }

    result.reserve(numerator.size());
    for(unsigned i=0, ni=numerator.size(); i<ni; i++)
      {
        result.push_back(numerator[i] / denominator);
      }
  }


  void BarycentricLagrangeInterpolator::
  eval_derivative(const double& x, const unsigned &deriv_order,
                  Vector<double>& result) const
  {
    eval_checks(x);

    // Only implemented for first derivatives...
    my_assert(deriv_order == 1);

    Vector<double> g(Values[0].size(), 0.0), dg(Values[0].size(), 0.0);
    double h = 0, dh = 0;

    for(unsigned i=0, ni=Locations.size(); i<ni; i++)
      {
        double temp = Weights[i]  / (x - Locations[i]);

        vector_add_ax(temp, Values[i], g);
        vector_add_ax(-1* temp / (x - Locations[i]), Values[i], dg);
        h += temp;
        dh += -1 * temp / (x - Locations[i]);
      }

    my_assert(h != 0);
    result.assign(Values[0].size(), 0.0);
    for(unsigned i=0, ni=g.size(); i<ni; i++)
      {
        result[i] = (dg[i] * h - g[i]*dh) / (h*h);
      }
  }


  void BarycentricLagrangeInterpolator::
  eval_checks(const double &x) const
  {
    // Check the everything has been set up right
    my_assert(Locations.size() != 0);
    my_assert(Locations.size() == Values.size());
    my_assert(Weights.size() == Locations.size());

    // Check that x is not a given location (trivial calculation gives
    // undefined value here, should be possible to output the appropriate
    // value if this functionality is needed...).
    my_assert(std::find(Locations.begin(), Locations.end(), x) == Locations.end());
  }

} // End of oomph namespace

#endif
