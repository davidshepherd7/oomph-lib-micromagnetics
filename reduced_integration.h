#ifndef OOMPH_REDUCED_INTEGRATION_H
#define OOMPH_REDUCED_INTEGRATION_H

#include <memory>

#include "../../src/generic/Vector.h"
#include "../../src/generic/integral.h"
#include "../../src/generic/elements.h"

#include "oomph_factories.h"
#include "array_interpolator.h"

namespace oomph
{

/// Class for reduced integration in micromagnetics, as defined in
/// e.g. Cimrak2008. Essentially we just use nodes as knots and \int
/// shape(x) dx as weights. This allows some nice conservation properties.
class ReducedIntegration : public Integral
{
private:

  /// Storage for pointer to the element we are integrating over
  const FiniteElement* Ele_pt;

  /// Storage for weight values
  Vector<double> Weight;

public:

  /// Dummy constructor
  ReducedIntegration()
  {
    Ele_pt = 0;
  }


  /// Real constructor
  ReducedIntegration(const FiniteElement* ele_pt)
  {
    Ele_pt = ele_pt;
    build();
  }


  /// Actually set up the integration scheme
  void build()
  {
    // Crude way to check the type of element, should really add a function
    // for this to FiniteElement.
    bool is_q_element = (dynamic_cast<const QElementGeometricBase*>(ele_pt()) != 0);

    // Construct integration scheme for integration of shape function. Use
    // factory so that we don't have to hard code (or template by) the
    // dimension/shape of the elements. Have to store in pointer to use
    // factory unfortunately so use auto ptr so that it is deleted
    // automatically. In c++11 this should be replaced by unique_ptr.
    std::auto_ptr<Integral> int_pt
      (Factories::gauss_integration_factory(ele_pt()->dim(),
                                            ele_pt()->nnode_1d(),
                                            is_q_element));

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
        const double J = ele_pt()->J_eulerian(s);

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

  /// Reduced integration means we don't actually need to interpolate, so
  /// this class has the interface of an interpolator without the actual
  /// calculations.
  template<unsigned VAL>
  class ReducedIntegrationInterpolator : public GeneralArrayInterpolator<VAL>
  {
  public:
    /// Constructor
    ReducedIntegrationInterpolator(const FiniteElement* const this_element,
                                   const unsigned& time_index=0)
      : GeneralArrayInterpolator<VAL>(this_element, time_index)
    {
      Node_pt = 0;
      Node_number = 0;

#ifdef PARANOID
      if(this->This_element->has_hanging_nodes())
        {
          std::string err = "Not designed for use with hanging nodes, you would need to use non raw nodal values functions, write a new but similar class if you really need hanging nodes.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
    }

    /// Virtual destructor
    virtual ~ReducedIntegrationInterpolator() {}

    void build(const Vector<double>& s_in)
      {
        Node_pt = this->This_element->get_node_at_local_coordinate(s_in);
        Node_number = this->This_element->get_node_number(Node_pt);

        GeneralArrayInterpolator<VAL>::build(s_in);
      }

    /// Interpolate evaluation position
    void interpolate_x()
    {
      for(unsigned j=0; j<this->Dim; j++)
        {
          this->X[j] = Node_pt->x(j);
        }
    }

    void interpolate_dxdt()
      {
        for(unsigned j=0; j<this->Dim; j++)
          {
            this->Dxdt[j] = Node_pt->dx_dt(j);
          }
      }

    void interpolate_values(const unsigned &start,
                            const unsigned &end)
    {
      for(unsigned j=start; j<end; j++)
        {
          this->Values[j] = Node_pt->value(j);
        }
    }

    void interpolate_dvaluesdt(const unsigned &start,
                               const unsigned &end)
    {
      for(unsigned j=start; j<end; j++)
        {
          // still need to interpolate in time
          this->Dvaluesdt[j] = 0;
          for(unsigned i_tm=0; i_tm<this->Nprev_value_derivative; i_tm++)
            {
              this->Dvaluesdt[j] += Node_pt->value(i_tm, j)
                * (*this->Ts_weights_pt)(1, i_tm);
            }
        }
    }

    void interpolate_dvaluesdx(const unsigned &i_val)
    {
      for(unsigned i_direc=0; i_direc<this->Dim; i_direc++)
        {
          this->Dvaluesdx[i_val][i_direc] =
            Node_pt->value(this->Time_index, i_val)
            * this->Dpsidx(Node_number, i_direc);
        }
    }

  private:

    /// The node which we are faking interpolation at
    Node* Node_pt;

    /// and it's number
    unsigned Node_number;

    /// Broken copy constructor
    ReducedIntegrationInterpolator(const ReducedIntegrationInterpolator& dummy)
    {BrokenCopy::broken_copy("ReducedIntegrationInterpolator");}

    /// Broken assignment operator
    void operator=(const ReducedIntegrationInterpolator& dummy)
    {BrokenCopy::broken_assign("ReducedIntegrationInterpolator");}

  };


} // End of oomph namespace

#endif
