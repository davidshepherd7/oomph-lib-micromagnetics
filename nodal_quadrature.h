#ifndef OOMPH_NODAL_QUADRATURE_H
#define OOMPH_NODAL_QUADRATURE_H

#include <memory>

#include "../../src/generic/Vector.h"
#include "../../src/generic/integral.h"
#include "../../src/generic/elements.h"

#include "oomph_factories.h"


namespace oomph
{

  /// Class for nodal quadrature in micromagnetics, as defined in
  /// e.g. Cimrak2008. Essentially we just use nodes as knots and \int
  /// shape(x) dx as weights. This allows some nice conservation properties.
  class NodalQuadrature : public Integral
  {
  protected:

    /// Storage for pointer to the element we are integrating over
    const FiniteElement* Ele_pt;

    /// Storage for weight values
    Vector<double> Weight;

  public:

    /// Dummy constructor
    NodalQuadrature()
    {
      Ele_pt = 0;
    }


    /// Real constructor
    NodalQuadrature(const FiniteElement* ele_pt)
    {
      Ele_pt = ele_pt;
      build();
    }


    /// Actually set up the integration scheme
    virtual void build()
    {
      // Construct integration scheme for integration of shape
      // function. Use factory so that we don't have to hard code (or
      // template by) the dimension/shape of the elements. Have to store in
      // pointer to use factory, use smart ptr so that it is deleted
      // automatically.
      std::auto_ptr<Integral> int_pt
        (Factories::gauss_integration_factory(ele_pt()->dim(),
                                              ele_pt()->nnode_1d(),
                                              ele_pt()->element_geometry()));

      // Get constants
      const unsigned nnode = ele_pt()->nnode();
      const unsigned dim = ele_pt()->dim();
      const unsigned nipt = int_pt->nweight();

      // Initialise storage for weights
      Weight.assign(nnode, 0);

      // Integrate each shape function over the element and store in
      // Weight. Loops reordered for efficiency.
      for(unsigned i=0; i<nipt; i++)
        {
          // Get integration point information
          Vector<double> s(dim);
          for(unsigned j=0; j<dim; j++)
            {s[j] = int_pt->knot(i, j);}
          const double w = int_pt->weight(i);

          // Jacobian of transformation at integration point
          const double J = J_eulerian(s);

          // Shape function values at integration point
          Shape shape(nnode);
          ele_pt()->shape(s, shape);

          // Add contribution for each node/shape function
          for(unsigned j=0; j<nnode; j++)
            {
              Weight[j] += w * J * shape(j);
            }
        }
    }

    /// Should the elemental jacobian be 1 for this quadrature?
    virtual bool unit_jacobian() const
    {
      return false;
    }

    /// Get Jacobian of transformation at point s (separate function so
    /// that it can be overloaded by other nodal quadrature schemes).
    virtual double J_eulerian(const Vector<double>& s) const
    {
      return ele_pt()->J_eulerian(s);
    }


    /// Get ele_pt, with paranoid check
    const FiniteElement* ele_pt() const
    {
#ifdef PARANOID
      if(Ele_pt == 0)
        {
          std::string err = "Element pointer not set.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return Ele_pt;
    }


    /// Set ele_pt
    void set_ele_pt(const FiniteElement* ele_pt)
    {
      Ele_pt = ele_pt;
    }


    /// Return number of integration points
    unsigned nweight() const {return ele_pt()->nnode();}


    /// Return jth s-coordinate of ith node/integration point.
    double knot(const unsigned &i, const unsigned &j) const
    {
      // This is dumb but it's the only way I can find to get the local
      // coordinate of a node.
      return knot(i)[j];
    }


    /// Return local coordinates of i-th intergration point.
    Vector<double> knot(const unsigned &i) const
    {
      Vector<double> s;
      ele_pt()->local_coordinate_of_node(i, s);
      return s;
    }


    /// Return weight of i-th node (= integral of shape dx)
    double weight(const unsigned &i) const {return Weight[i];}

  };


  /// Class to handle ??ds
  class RescaledNodalQuadrature : public NodalQuadrature
  {
  public:
    /// Constructor
    RescaledNodalQuadrature() {}

    /// Virtual destructor
    virtual ~RescaledNodalQuadrature() {}

    /// Real constructor
    RescaledNodalQuadrature(const FiniteElement* ele_pt,
                               const double& mean_element_volume)
    {
      Ele_pt = ele_pt;
      Mean_element_volume = mean_element_volume;
      build();
    }

    void build()
    {
      // Build as normal
      NodalQuadrature::build();

      // Rescale by mean element volume
      const unsigned ni = Weight.size();
      for(unsigned i=0; i<ni; i++)
        {
          Weight[i] /= Mean_element_volume;
        }
    }

  protected:

    double Mean_element_volume;

  private:
    /// Broken copy constructor
    RescaledNodalQuadrature(const RescaledNodalQuadrature& dummy)
    {BrokenCopy::broken_copy("RescaledNodalQuadrature");}

    /// Broken assignment operator
    void operator=(const RescaledNodalQuadrature& dummy)
    {BrokenCopy::broken_assign("RescaledNodalQuadrature");}

  };


  /// As nodal quadrature but do the weight calculation on a standard
  /// element (local coordinates) without mapping back to global. This is
  /// done by setting the Jacobian of the transformation (used for the
  /// weight integration only) to 1.
  class LocalNodalQuadrature : public NodalQuadrature
  {
  public:
    /// Constructor
    LocalNodalQuadrature() {}

    /// Virtual destructor
    virtual ~LocalNodalQuadrature() {}

    /// Real constructor
    LocalNodalQuadrature(const FiniteElement* ele_pt)
      {
        Ele_pt = ele_pt;
        build();
      }

    virtual double J_eulerian(const Vector<double>& s) const
    {
      return 1;
    }

  private:
    /// Broken copy constructor
    LocalNodalQuadrature(const LocalNodalQuadrature& dummy)
    {BrokenCopy::broken_copy("LocalNodalQuadrature");}

    /// Broken assignment operator
    void operator=(const LocalNodalQuadrature& dummy)
    {BrokenCopy::broken_assign("LocalNodalQuadrature");}

  };


  /// As nodal quadrature but do everything on the global level. Not sure
  /// if this works...
  class GlobalNodalQuadrature : public NodalQuadrature
  {
  public:
    /// Constructor
    GlobalNodalQuadrature() {}

    /// Virtual destructor
    virtual ~GlobalNodalQuadrature() {}

    /// Real constructor
    GlobalNodalQuadrature(const FiniteElement* ele_pt)
    {
      Ele_pt = ele_pt;
      build();
    }

    /// Should the elemental jacobian be 1 for this quadrature?
    virtual bool unit_jacobian() const
    {
      return true;
    }

  private:
    /// Broken copy constructor
    GlobalNodalQuadrature(const GlobalNodalQuadrature& dummy)
    {BrokenCopy::broken_copy("GlobalNodalQuadrature");}

    /// Broken assignment operator
    void operator=(const GlobalNodalQuadrature& dummy)
    {BrokenCopy::broken_assign("GlobalNodalQuadrature");}

  };


  /// Nodal quadrature with all weights set to 1, should work better with
  /// integration on standard elements.
  class UnitWeightNodalQuadrature : public NodalQuadrature
  {
  public:
    /// Constructor
    UnitWeightNodalQuadrature() {}

    /// Virtual destructor
    virtual ~UnitWeightNodalQuadrature() {}

    /// Real constructor
    UnitWeightNodalQuadrature(const FiniteElement* ele_pt)
      : NodalQuadrature(ele_pt) {}

    virtual void build()
    {
      // Set all weights to 1
      Weight.assign(ele_pt()->nnode(), 1);
    }

  private:
    /// Broken copy constructor
    UnitWeightNodalQuadrature(const UnitWeightNodalQuadrature& dummy)
    {BrokenCopy::broken_copy("UnitWeightNodalQuadrature");}

    /// Broken assignment operator
    void operator=(const UnitWeightNodalQuadrature& dummy)
    {BrokenCopy::broken_assign("UnitWeightNodalQuadrature");}

  };

} // End of oomph namespace

#endif
