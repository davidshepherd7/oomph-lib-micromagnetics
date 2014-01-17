 #ifndef OOMPH_MICROMAGNETICS_BOUNDARY_ELEMENT_H
#define OOMPH_MICROMAGNETICS_BOUNDARY_ELEMENT_H


#include "./micromagnetics_element.h"
#include "./variable_order_quadrature.h"
#include <functional>
#include <algorithm>

#include "../../src/generic/mesh.h"

// magpar definitions
#include "./magpar_requirements.h"

// easier to use vector functions for now (probably not optimal...)
#include "./vector_helpers.h"

namespace oomph
{
  using namespace MathematicalConstants;
  using namespace VectorOps;

  //===========================================================================
  ///
  //===========================================================================
  class MicromagBEMElementEquations : public FaceElement
  {
  public:

    /// \short Constructor, takes the pointer to the bulk element and the
    /// index of the face to which the element is attached.
    MicromagBEMElementEquations(FiniteElement* const bulk_el_pt,
                                const int& face_index);

    ///\short  Broken empty constructor
    MicromagBEMElementEquations()
    {
      throw OomphLibError("Don't call empty constructor for MicromagBEMElementEquations",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    /// Broken copy constructor
    MicromagBEMElementEquations(const MicromagBEMElementEquations& dummy)
    {
      BrokenCopy::broken_copy("MicromagBEMElementEquations");
    }

    /// Broken assignment operator
    void operator=(const MicromagBEMElementEquations&)
    {
      BrokenCopy::broken_assign("MicromagBEMElementEquations");
    }

    unsigned self_test()
    {
      // Just pass it (we can't let this element call the normal self-tests
      // because it messes with the integration scheme).
      return 0;
    }

    /// \short Add the element's contribution to its residual vector and
    /// its Jacobian matrix. Use numerical or analytical integration as
    /// decided by second argument.
    void fill_in_contribution_to_boundary_matrix(DenseMatrix<double> &boundary_matrix,
                                                 bool use_numerical_integration) const
    {
      if(use_numerical_integration)
        {
          fill_in_be_contribution_adaptive(boundary_matrix);
        }
      else
        {
          fill_in_be_contribution_analytic(boundary_matrix);
        }
    }

    /// Function determining how to block the Jacobian. Just move boundary
    /// values of phi into their own blocks. ??ds
    void get_dof_numbers_for_unknowns(std::list<std::pair<unsigned long,unsigned> >&
                                      block_lookup_list)
    {
      // unsigned phi_1_index = phi_1_index_micromag();
      // unsigned phi_index = phi_index_micromag();

      // // We just want to move the boundary values of phi into their own blocks
      // for(unsigned nd=0; nd<nnode(); nd++)
      //  {
      //   int phi_1_local = micromag_bulk_element_pt()->nodal_local_eqn(nd,phi_1_index);
      //   if(phi_1_local >=0)
      //    {
      //     std::pair<unsigned,unsigned> block_lookup;
      //     block_lookup.first = micromag_bulk_element_pt()->eqn_number(phi_1_local);
      //     block_lookup.second = 0;
      //     block_lookup_list.push_front(block_lookup);
      //    }

      //   int phi_local = micromag_bulk_element_pt()->nodal_local_eqn(nd,phi_index);
      //   if(phi_local >=0)
      //    {
      //     std::pair<unsigned,unsigned> block_lookup;
      //     block_lookup.first = micromag_bulk_element_pt()->eqn_number(phi_local);
      //     block_lookup.second = 1;
      //     block_lookup_list.push_front(block_lookup);
      //    }
      //  }

      //??ds - removed for now!

      std::string err = "Not implemented";
      throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                          OOMPH_CURRENT_FUNCTION);
    }

    /// Number of dof types added by this element. We add two: one each for the
    /// boundary values of phi and phi_1.
    unsigned ndof_types() {return 0;} //??ds

    /// Output function - do nothing
    void output(std::ostream &outfile, const unsigned &n_plot=5) const {}

    /// \short Function giving the normal derivative of the Green's function.
    /// The parameters x and y are the positions of source and point of interest
    /// (interchangable), n is the normal unit vector out of the surface.
    double green_normal_derivative(const Vector<double>& x,
                                   const Vector<double>& y,
                                   const Vector<double>& n) const;

    /// \short Function giving the normal derivative of the Green's
    /// function.  The parameters x and y are the positions of source and
    /// point of interest (interchangable), n is the normal unit vector out
    /// of the surface. Optimised but harder to read and potentially less
    /// safe version.
    double green_normal_derivative_opt(const Vector<double>& x,
                                   const Vector<double>& y,
                                   const Vector<double>& n) const;

    /// Const access function for mesh pointer
    const Mesh* boundary_mesh_pt() const {return Boundary_mesh_pt;}

    /// Set function for mesh pointer
    void set_boundary_mesh_pt(const Mesh* const boundary_mesh_pointer)
    {Boundary_mesh_pt = boundary_mesh_pointer;}

    /// Get the max difference between two vectors relative to that element of vector1
    double max_rel_error(const Vector<double> &vector1,
                         const Vector<double> &vector2) const;

    /// dump out important values for testing
    void dump_values(const Vector< Vector<double> > &x_kn,
                     const Vector<double> &source_node_x,
                     const Vector<Vector<double> > &normal) const;


    /// ??ds Debugging function
    bool normals_match(const Vector<unsigned> &node_list) const;


  protected:

    /// The number of values to be stored at each boundary element node
    inline unsigned required_nvalue(const unsigned &n) const
    {return 0;}

  private:

    /// \short Pointer to the boundary mesh (needed to access nodes
    /// outside of this element for calculation of boundary matrix).
    const Mesh* Boundary_mesh_pt;

    /// Add the element's contribution to the boundary element matrix using adaptive quadrature.
    void fill_in_be_contribution_adaptive(DenseMatrix<double> &boundary_matrix) const;

    /// Add the element's contribution to the boundary element matrix using
    /// analytic calculations by calling analytic_integral_dgreendn with
    /// appropriate input for the element.
    void fill_in_be_contribution_analytic(DenseMatrix<double> &boundary_matrix) const;

    /// Calculate the contribution of a triangular region to the boundary
    /// element matrix using analytic calculations from Lindholm1984.
    void analytic_integral_dgreendn_triangle(const Vector<Vector<double> >& x_nds,
                                             const Vector<unsigned>& l,
                                             DenseMatrix<double>& boundary_matrix) const;

  };


  template<unsigned DIM, unsigned NNODE_1D>
  class QMicromagBEMElement : public virtual FaceGeometry<QElement<DIM, NNODE_1D> >,
                              public virtual MicromagBEMElementEquations
  {

  public:
    QMicromagBEMElement(FiniteElement* const &bulk_el_pt,
                               const int& face_index)
      : FaceGeometry<QElement<DIM, NNODE_1D> >(),
        MicromagBEMElementEquations(bulk_el_pt, face_index)
    {}

  };


  template<unsigned DIM, unsigned NNODE_1D>
  class TMicromagBEMElement : public virtual FaceGeometry<TElement<DIM, NNODE_1D> >,
                              public virtual MicromagBEMElementEquations
  {

  public:
    TMicromagBEMElement(FiniteElement* const &bulk_el_pt,
                               const int& face_index)
      : FaceGeometry<TElement<DIM, NNODE_1D> >(),
        MicromagBEMElementEquations(bulk_el_pt, face_index)
    {}

  };

}

#endif
