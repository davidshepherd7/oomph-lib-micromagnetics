#ifndef OOMPH_MICROMAGNETICS_FLUX_ELEMENT_H
#define OOMPH_MICROMAGNETICS_FLUX_ELEMENT_H

// class MicromagEquations;
#include "./micromagnetics_element.h"

// Note: we can't use a .cc file here easily because then we would need to
// instantiate the template for every possible element, of which there are
// quite a few...

namespace oomph
{


  //======================================================================
  /// A face element to impose the flux boundary condition on the potential for
  /// micromagnetics. It is not possible to use the already existing Poisson flux
  /// elements becuase they assume that:
  ///
  /// 1) The "parent" element will be called PoissonElement.
  ///
  /// 2) The Jacobian contribution due to the flux elements is always
  /// zero. This is not true for our case since the flux is fixed as m.n
  /// (where m is the magnetisation). So differentiation wrt m gives a
  /// non-zero result.
  // ======================================================================
  template <class ELEMENT>
  class MicromagFluxElement : public virtual FaceGeometry<ELEMENT>,
                              public virtual FaceElement
  {

  public:

    /// \short Constructor, takes the pointer to the "bulk" element and the
    /// index of the face to which the element is attached.
    MicromagFluxElement(FiniteElement* const &bulk_el_pt,
                        const int& face_index);

    /// Add the element's contribution to its residual vector
    void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      // Call the generic residuals function with flag set to 0 using a dummy
      // matrix argument
      fill_in_generic_residual_contribution_fluxes
        (residuals,GeneralisedElement::Dummy_matrix,0);
    }

    /// \short Add the element's contribution to its residual vector and its
    /// Jacobian matrix
    void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                          DenseMatrix<double> &jacobian)
    {
      //Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_fluxes(residuals,jacobian,1);
    }


    /// \short Output function -- forward to broken version in FiniteElement.
    void output(std::ostream &outfile, const unsigned &n_plot=5)
    {FiniteElement::output(outfile,n_plot);}

    /// \short C-style output function -- forward to broken version in
    /// FiniteElement.
    void output(FILE* file_pt, const unsigned &n_plot=5)
    {FiniteElement::output(file_pt,n_plot);}

    MicromagEquations* bulk_element_pt() const
    {
      return Bulk_element_pt;
    }

    unsigned m_index_micromag(const unsigned &i) const
    {
      return bulk_element_pt()->m_index_micromag(i);
    }

    unsigned ndof_types() {return 0;}

    /// \short Element makes no changes to dof numbering - do nothing.
    void get_dof_numbers_for_unknowns
    (std::list<std::pair<unsigned long,unsigned> >& block_lookup_list) {}

    /// Return  m . dm/dn at s.
    double interpolated_mdotdmdn_micromag(const Vector<double> &s) const
    {
      // dm_i/dn = dm_i/dx . normal

      std::string err = "Untested with midpoint, not using interpolator so prob. no good";
      throw OomphLibWarning(err, OOMPH_EXCEPTION_LOCATION,
                          OOMPH_CURRENT_FUNCTION);

      // Get shape functions
      Shape psi(nnode()); DShape dpsidx(nnode(),dim());
      dshape_local(s,psi,dpsidx);

      // Get unit normal
      Vector<double> normal(nodal_dimension(),0.0);
      outer_unit_normal(s, normal);

      // For each m direction
      Vector<double> dmdn(3,0.0), m(3,0.0);
      for(unsigned i_m=0; i_m<3; i_m++)
        {

          // Interpolate dm_idx and m_i
          Vector<double> dmidx(nodal_dimension(),0.0);
          for(unsigned l=0, nl=nnode(); l<nl; l++)
            {
              for(unsigned j=0; j<dim(); j++)
                {
                  dmidx[j] += nodal_value(l,m_index_micromag(i_m))*dpsidx(l,j);
                }

              m[i_m] += nodal_value(l,m_index_micromag(i_m)) * psi(l);
            }

          // Take dot product
          dmdn[i_m] = VectorOps::dot(dmidx, normal);
        }

      // Get m . dm/dn (because it's easy now and we probably need it)
      return VectorOps::dot(m,dmdn);
    }

  protected:

    /// \short Function to compute the shape and test functions and to return
    /// the Jacobian of mapping between local and global (Eulerian)
    /// coordinates
    double shape_test(const Vector<double> &s, Shape &psi, Shape &test)
      const
    {
      //Get the shape functions
      shape(s,psi);

      //Set the test functions to be the same as the shape functions
      test = psi;

      //Return the value of the jacobian
      return J_eulerian(s);
    }

    /// \short Add the element's contribution to its residual vector. When
    /// flag=1 (or 0): do (or don't) compute the contribution to the Jacobian as
    /// well.
    void fill_in_generic_residual_contribution_fluxes
    (Vector<double> &residuals, DenseMatrix<double> &jacobian,
     const unsigned& flag);

    void fill_in_contribution_to_mass_matrix
    (Vector<double> &residuals, DenseMatrix<double> &mmatrix);


  private:

    /// Pointer to the attached bulk element
    MicromagEquations* Bulk_element_pt;

    ///\short Broken empty constructor
    MicromagFluxElement()
    {
      throw OomphLibError
        ("Don't call empty constructor for MicromagFluxElement",
         OOMPH_CURRENT_FUNCTION,
         OOMPH_EXCEPTION_LOCATION);
    }

    /// Broken copy constructor
    MicromagFluxElement(const MicromagFluxElement& dummy)
    {
      BrokenCopy::broken_copy("MicromagFluxElement");
    }

    /// Broken assignment operator
    void operator=(const MicromagFluxElement&)
    {
      BrokenCopy::broken_assign("MicromagFluxElement");
    }
  };


  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the
  /// index of the fixed local coordinate and its value represented
  /// by an integer (+/- 1), indicating that the face is located
  /// at the max. or min. value of the "fixed" local coordinate
  /// in the bulk element.
  //===========================================================================
  template<class ELEMENT>
  MicromagFluxElement<ELEMENT>::
  MicromagFluxElement(FiniteElement* const &bulk_el_pt, const int &face_index)
    : FaceGeometry<ELEMENT>(), FaceElement(),
      Bulk_element_pt(checked_dynamic_cast<MicromagEquations*>(bulk_el_pt))
  {
    // Let the bulk element build the FaceElement, i.e. setup the pointers
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_element_pt()->build_face_element(face_index,this);
  }


  //===========================================================================
  /// Compute the element's residual vector and the Jacobian matrix.
  //===========================================================================
  template<class ELEMENT>
  void MicromagFluxElement<ELEMENT>::
  fill_in_generic_residual_contribution_fluxes
  (Vector<double> &residuals,
   DenseMatrix<double> &jacobian,
   const unsigned& flag)
  {
    // If there is some surface anisotropy then throw an error: not
    // implemented.
#ifdef PARANOID
    if(bulk_element_pt()->magnetic_parameters_pt()->surface_anisotropy_enabled())
      {
        std::ostringstream error_msg;
        error_msg << "Surface anisotropy is not implemented";
        throw OomphLibError(error_msg.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

    const unsigned n_node = nnode();
    const unsigned n_intpt = integral_pt()->nweight();
    const unsigned dim = node_pt(0)->ndim();

    Shape psi(n_node), test(n_node);

    // indicies:
    unsigned phi1_index = bulk_element_pt()->phi_1_index_micromag();
    Vector<unsigned> m_index(3, 0);
    for(unsigned j=0; j<3; j++)
      {
        m_index[j] = bulk_element_pt()->m_index_micromag(j);
      }


    //Loop over the integration points
    //--------------------------------
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {

        // Local coords
        Vector<double> s(dim-1);
        for(unsigned i=0; i<(dim-1); i++) {s[i] = integral_pt()->knot(ipt, i);}

        // Create interpolator (deals with midpoint weirdness, so still
        // very useful even for such a simple residual).
        FaceElementArrayInterpolator<5> intp(this, s);

        const double W = integral_pt()->weight(ipt) * intp.j();

        // Get normal vector, force 3 entries for ease of combination with
        // m which always has 3 entries.
        Vector<double> normal(dim,0.0);
        outer_unit_normal(s, normal);
        normal.resize(3, 0.0); // make 3d, initialise new values to zero.

        // Get mdotn at this point
        double mdotn = 0.0;
        for(unsigned j=0; j<dim; j++)
          {
            mdotn += intp.value(m_index[j]) * normal[j];
          }

        // Loop over the test functions doing residual and Jacobian
        // contributions.
        for(unsigned l=0;l<n_node;l++)
          {
            // Get indicies for phi_1 equation
            int phi_1_eqn = nodal_local_eqn (l, phi1_index);

            // Skip if a dirichlet b.c.
            if(phi_1_eqn < 0) continue;

            // Add contribution to phi residual
            residuals[phi_1_eqn] += mdotn * intp.test(l) * W;

            // Skip rest if Jacobian not requested
            if(!flag) continue;

            for(unsigned l2=0; l2<n_node; l2++)
              {
                // Loop over which m we are differentiating w.r.t
                for(unsigned j=0; j<3; j++)
                  {
                    unsigned m_unknown = nodal_local_eqn(l2, m_index[j]);

                    // Skip rest is m[j] pinned
                    if(m_unknown < 0) continue;

                    // phi_1 residual w.r.t m[j]
                    jacobian(phi_1_eqn, m_unknown)
                      += intp.psi(l2) * intp.test(l) * normal[j] * W;

                  }
              }
          }


      }

  }

  /// This might go better inside generic get jacobian etc. once I write
  /// it.
  template<class ELEMENT>
  void MicromagFluxElement<ELEMENT>::fill_in_contribution_to_mass_matrix
  (Vector<double> &residuals, DenseMatrix<double> &mmatrix)
  {
    const unsigned n_node = this->nnode();
    const unsigned eldim = this->dim();
    const unsigned n_unknowns = ndof_types();

    Shape psi(n_node), test(n_node);
    Vector<double> s(eldim);

    //Loop over the integration points
    for(unsigned ipt=0, nipt=this->integral_pt()->nweight(); ipt<nipt; ipt++)
      {
        // Get position
        for(unsigned j=0; j<eldim; j++)
          {s[j] = this->integral_pt()->knot(ipt,j);}

        // Get shape/test/coord transform Jacobian
        double J = J_eulerian(s);
        shape(s, psi);
        test = psi;

        double W = this->integral_pt()->weight(ipt) * J;

        for(unsigned l=0;l<n_node;l++)
          {
            //Loop over the unknowns
            for(unsigned i=0;i<n_unknowns;i++)
              {
                int local_eqn = this->nodal_local_eqn(l,i);
                if(local_eqn < 0) continue;

                for(unsigned l2=0;l2<n_node;l2++)
                  {
                    int local_unknown = this->nodal_local_eqn(l2,i);
                    if(local_unknown < 0) continue;

                    //??ds do we need to multiply by
                    //d_valuederivative_evaltime_by_dvalue_np1?
                    mmatrix(local_eqn, local_unknown) +=
                      psi(l2)*test(l)*W;
                  }
              }
          }
      }
  }

} // End of oomph namespace

#endif
