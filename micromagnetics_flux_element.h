#ifndef OOMPH_MICROMAGNETICS_FLUX_ELEMENT_H
#define OOMPH_MICROMAGNETICS_FLUX_ELEMENT_H

/*
  A face element to impose the flux boundary condition on the potential for
  micromagnetics. It is not possible to use the already existing Poisson flux
  elements becuase they assume that:

  1) The "parent" element will be called PoissonElement.

  2) The Jacobian contribution due to the flux elements is always zero. This is
  not true for our case since the flux is fixed as m.n (where m is the
  magnetisation). So differentiation wrt m gives a non-zero result.

*/


//??ds
#include <iomanip>

using namespace oomph;

namespace oomph
{

  //======================================================================
  ///
  //======================================================================
  template <class ELEMENT>
  class MicromagFluxElement : public virtual FaceGeometry<ELEMENT>,
                              public virtual FaceElement
  {

  public:

    /// \short Constructor, takes the pointer to the "bulk" element and the
    /// index of the face to which the element is attached.
    MicromagFluxElement(ELEMENT* const &bulk_el_pt,
                        const int& face_index);

    ///\short  Broken empty constructor
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

    // /// \short Specify the value of nodal zeta from the face geometry
    // /// The "global" intrinsic coordinate of the element when
    // /// viewed as part of a geometric object should be given by
    // /// the FaceElement representation, by default (needed to break
    // /// indeterminacy if bulk element is SolidElement)
    // //??ds don't think I need this....
    // double zeta_nodal(const unsigned &n, const unsigned &k,
    //                 const unsigned &i) const
    // {return FaceElement::zeta_nodal(n,k,i);}

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


    void fill_in_bulk_contribution_to_face_jacobian
    (DenseMatrix<double>& jacobian) const;

    /// \short Output function -- forward to broken version in FiniteElement.
    void output(std::ostream &outfile, const unsigned &n_plot=5)
    {FiniteElement::output(outfile,n_plot);}

    /// \short C-style output function -- forward to broken version in
    /// FiniteElement.
    void output(FILE* file_pt, const unsigned &n_plot=5)
    {FiniteElement::output(file_pt,n_plot);}

    /// Access function for the pointer the to bulk element
    ELEMENT* bulk_element_pt() const {return Bulk_element_pt;}

    unsigned m_index_micromag(const unsigned &i) const
    {
      return bulk_element_pt()->m_index_micromag(i);
    }

    unsigned ndof_types()
    {return 0;}

    /// \short Element makes no changes to dof numbering - do nothing.
    void get_dof_numbers_for_unknowns(std::list<std::pair<unsigned long,unsigned> >&
                                      block_lookup_list)
    {}

    /// Return  m . dm/dn at s.
    double interpolated_mdotdmdn_micromag(const Vector<double> &s) const
    {
      // dm_i/dn = dm_i/dx . normal

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

  private:

    ///The spatial dimension of the problem
    unsigned Nodal_dim;

    /// Pointer to the attached bulk element
    ELEMENT* Bulk_element_pt;
  };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////



  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the
  /// index of the fixed local coordinate and its value represented
  /// by an integer (+/- 1), indicating that the face is located
  /// at the max. or min. value of the "fixed" local coordinate
  /// in the bulk element.
  //===========================================================================
  template<class ELEMENT>
  MicromagFluxElement<ELEMENT>::
  MicromagFluxElement(ELEMENT* const &bulk_el_pt,
                      const int &face_index) :
    FaceGeometry<ELEMENT>(),
    FaceElement(),
    Bulk_element_pt(bulk_el_pt)
  {
    // Let the bulk element build the FaceElement, i.e. setup the pointers
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_element_pt()->build_face_element(face_index,this);

    // Extract the dimension of the problem from the dimension of
    // the first node
    Nodal_dim = this->node_pt(0)->ndim();
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
    // If no surface anisotropy then nothing to do here
    if(!bulk_element_pt()->magnetic_parameters_pt()->surface_anisotropy_enabled())
      {
        return;
      }

    std::ostringstream error_msg;
    error_msg << "Surface anisotropy is not implemented";
    throw OomphLibError(error_msg.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);

    const unsigned n_node = nnode();
    const unsigned n_intpt = integral_pt()->nweight();

    Shape psi(n_node), test(n_node);

    //const unsigned phi_1_index = bulk_element_pt()->phi_1_index_micromag();
    Vector<unsigned> m_index(3);
    for(unsigned i=0; i<3; i++)
      m_index[i] = bulk_element_pt()->m_index_micromag(i);

    // Get coefficients
    const double exch_c = bulk_element_pt()->exchange_coeff();
    const double llg_precess_c = bulk_element_pt()->llg_precession_coeff();

    const double boundary_exch_debug = bulk_element_pt()->magnetic_parameters_pt()
      ->boundary_exchange_debug_coeff();

    //Loop over the integration points
    //--------------------------------
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {
        // Local coords
        Vector<double> s(Nodal_dim-1);
        for(unsigned i=0;i<(Nodal_dim-1);i++)
          s[i] = integral_pt()->knot(ipt,i);

        const double J = shape_test(s,psi,test);
        const double W = integral_pt()->weight(ipt) * J;

        // Get the shape/test functions and derrivatives
        shape_test(s,psi,test);

        Vector<double> itp_x(Nodal_dim,0.0), itp_m(3,0.0);
        for(unsigned l=0; l<n_node; l++)
          {
            for(unsigned j=0; j<Nodal_dim; j++)
              itp_x[j] += nodal_position(l,j)*psi(l);
            for(unsigned j=0; j<3; j++)
              itp_m[j] += nodal_value(l,m_index[j])*psi(l);
          }

        // We have to ask the bulk element to calculate dmdx for us because the
        // derivatives of shape functions for nodes in the bulk but not in the
        // face could be non-zero on the boundary.
        Vector<double> s_bulk(3,0.0);
        get_local_coordinate_in_bulk(s,s_bulk);
        DenseDoubleMatrix itp_dmdx(3,3,0.0);
        bulk_element_pt()->interpolated_dmdx_micromag(s_bulk,itp_dmdx);

        // Vector<double> low_dim_normal(Nodal_dim,0.0), normal(3,0.0);
        Vector<double> normal(Nodal_dim,0.0);
        outer_unit_normal(s,normal);
        // for(unsigned j=0; j<Nodal_dim; j++) normal[j] = low_dim_normal[j];

        // Some pre-calculations:
        Vector<double> dmdn(3,0.0);
        for(unsigned i=0; i<3; i++)
          for(unsigned j=0; j<Nodal_dim; j++)
            dmdn[i] += normal[j] * itp_dmdx(i,j);

        // std::cout << "from residual dmdn at x = " << itp_x << " is " << dmdn << std::endl
        //           << "       m is " << itp_m << std::endl;

        Vector<double> mxdmdn(3,0.0);
        VectorOps::cross(itp_m,dmdn,mxdmdn);
        //double mdotn = VectorOps::dot(normal,itp_m);


        // Loop over the test functions doing residual and Jacobian
        // contributions.
        for(unsigned l=0;l<n_node;l++)
          {
            // Get indicies for phi_1 and m equations
            // int phi_1_eqn = nodal_local_eqn (l,phi_1_index);
            Vector<int> m_eqn(3,0);

            //??ds
            // // Phi contribution (if not a dirichlet b.c.)
            // if(phi_1_eqn >= 0)
            //   residuals[phi_1_eqn] += mdotn*test(l)*W;

            // Exchange contribution: (m x (dm/dn)) * test
            for(unsigned i=0; i<3; i++)
              {
                m_eqn[i] = nodal_local_eqn(l,m_index[i]);
                if(m_eqn[i] >=0)
                  {
                    residuals[m_eqn[i]] += llg_precess_c * exch_c
                      * boundary_exch_debug
                      * mxdmdn[i] * test(l) * W;
                  }
              }

            //??ds still need to deal with surface anisotropy...


            // // Jacobian (phi_1 w.r.t m_i)
            // if(!flag) continue;
            // for(unsigned l2=0; l2<n_node; l2++)
            //   {
            //     // Loop over which m we are differentiating w.r.t
            //     Vector<unsigned> m_unknown(3,0.0);
            //     for(unsigned j=0; j<3; j++)
            //       {
            //         m_unknown[j] = nodal_local_eqn(l2,m_index[j]);
            //         if(!(m_unknown[j] >= 0)) continue;

            //         // phi_1 residual w.r.t m_j
            //         if(phi_1_eqn >= 0)
            //           {
            //             jacobian(phi_1_eqn,m_unknown[j])
            //               += psi(l2) * test(l) * normal[j] * W;
            //           }

            //         // m residuals w.r.t m_j due to exchange boundary
            //         // effects has to be handled by bulk element since it
            //         // involves derivatives of m. So it is affected by
            //         // nodes which are not in the face element.


            //       }
            //   }
          }

        // std::cout <<  normal << std::endl;

      }

    // std::cout << "residuals are:" <<  residuals << std::endl;
  }

  template<class ELEMENT>
  void MicromagFluxElement<ELEMENT>::
  fill_in_bulk_contribution_to_face_jacobian(DenseMatrix<double>& jacobian) const
  {

    // If no surface anisotropy then nothing to do here
    if(!bulk_element_pt()->magnetic_parameters_pt()->surface_anisotropy_enabled())
      {
        return;
      }

    std::ostringstream error_msg;
    error_msg << "Surface anisotropy is not implemented";
    throw OomphLibError(error_msg.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);


    // Get some values
    const ELEMENT* be_pt = bulk_element_pt();
    const unsigned bulk_nnode = be_pt->nnode();
    const unsigned nodal_dim = this->node_pt(0)->ndim();
    const unsigned face_ele_dim = this->dim();

    // Integration using surface element points/scheme
    for(unsigned ipt=0, nipt=this->integral_pt()->nweight(); ipt<nipt;ipt++)
      {
        // Get integration data (knot, weight)
        Vector<double> face_s(face_ele_dim,0.0), bulk_s(nodal_dim,0.0);
        for(unsigned j=0; j<face_ele_dim; j++) face_s[j] = this->integral_pt()->knot(ipt,j);
        this->get_local_coordinate_in_bulk(face_s, bulk_s);

        double w = this->integral_pt()->weight(ipt);

        // Get shape functions and  for the bulk element.
        Shape psi(bulk_nnode), test(bulk_nnode);
        DShape dpsidx(bulk_nnode,nodal_dim), dtestdx(bulk_nnode,nodal_dim);
        be_pt->dshape_dtest(bulk_s,psi,dpsidx,test,dtestdx);

        // Get the Jacobian of transformation (from reference element) for
        // the face element (since that is the length/area/volume we are
        // integrating over).
        double J = J_eulerian(face_s);
        double W = w * J;

        // Get the unit normal
        Vector<double> normal(nodal_dim, 0.0);
        this->outer_unit_normal(face_s,normal);
        // std::cout << normal << std::endl;


        // Calculate n dot gradpsi for each node's shape function
        Vector<double> dpsidn(bulk_nnode,0.0);
        for(unsigned l=0; l<bulk_nnode; l++)
          {
            for(unsigned j=0; j<nodal_dim; j++)
              {
                dpsidn[l] += normal[j] * dpsidx(l,j);
              }
          }

        // std::cout << std::endl << dpsidn << std::endl << std::endl;

        // Interpolate magnetisation and dmdn
        Vector<double> itp_m(3,0.0), itp_dmdn(3,0.0);
        for(unsigned l=0; l<bulk_nnode; l++)
          {
            for(unsigned j=0; j<3; j++)
              {
                itp_m[j] += be_pt->nodal_value(l, be_pt->m_index_micromag(j))
                  * psi(l);
                itp_dmdn[j] += be_pt->nodal_value(l, be_pt->m_index_micromag(j))
                  * dpsidn[l];
              }
          }

        Vector<double> itp_x(Nodal_dim,0.0);
        for(unsigned l=0; l<bulk_nnode; l++)
          {
            for(unsigned j=0; j<Nodal_dim; j++)
              {
                itp_x[j] += be_pt->nodal_position(l,j)*psi(l);
              }
          }

        // std::cout << "from jacobian dmdn at x = " << itp_x << " is " << itp_dmdn << std::endl
        //           << "       m is " << itp_m << std::endl;

        // Get constants
        double exch_c = be_pt->exchange_coeff();
        double llg_precess_c = be_pt->llg_precession_coeff();

        const double boundary_exch_debug = bulk_element_pt()->magnetic_parameters_pt()
          ->boundary_exchange_debug_coeff();

        // Do the Jacobian calculation
        // ============================================================

        // Loop over the variable we differentiate by (nodes).
        for(unsigned l_var=0; l_var<bulk_nnode; l_var++)
          {

            // Calculate some things...
            Vector<double> stuff(3,0.0);
            for(unsigned i=0; i<3; i++)
              {
                stuff[i] = psi(l_var) * itp_dmdn[i] - itp_m[i] * dpsidn[l_var];
              } //??ds

            // std::cout << itp_dmdn << std::endl;
            // std::cout << itp_m << std::endl;

            // Loop over the variable we differentiate by (components).
            for(unsigned i_var=0; i_var<3; i_var++)
              {
                // Get the eqn number for the variable, skip if pinned.
                const int var_en =
                  be_pt->nodal_local_eqn(l_var,be_pt->m_index_micromag(i_var));
                if(var_en < 0) continue;

                // Unit vector in i_var direction
                Vector<double> jhat(3,0.0); jhat[i_var] = 1.0;

                Vector<double> jxstuff(3,0.0);
                VectorOps::cross(jhat, stuff, jxstuff); //??ds

                // for(unsigned i=0; i<3; i++)
                //   {

                //     std::cout << std::setw( 12 ) << jxstuff[i];
                //   }
                // std::cout << std::endl;

                // Loop over residuals (one for each node).
                for(unsigned l_res=0; l_res<bulk_nnode; l_res++)
                  {
                    // Check if this node is in the face element (otherwise there
                    // is no residual contribution so no Jacobian contribution).
                    if(this->get_node_number(be_pt->node_pt(l_res)) == -1) continue;
                    //??ds possibly a cleaner way of doing this? Loop over face
                    //nodes and convert number?

                    // Loop over components of residual.
                    for(unsigned i_res=0; i_res<3; i_res++)
                      {
                        // Get eqn number, skip if pinned.
                        const int res_en
                          = be_pt->nodal_local_eqn(l_res, be_pt->m_index_micromag(i_res));
                        if(res_en < 0) continue;

                        jacobian(res_en,var_en) += llg_precess_c * exch_c
                          * boundary_exch_debug
                          * jxstuff[i_res] * test(l_res) * W;

                        //??ds which value of J (=W/w) should we be using?
                        //bulk probably... but that is different to the one
                        //used for residuals...

                        // jacobian(res_en,var_en) += W
                        //   * jxstuff[i_res]; //??ds wrong!
                      }
                  }

              }
          }

      }

    //??ds

    // // Loop over residuals (one for each node).
    // for(unsigned l_res=0; l_res<bulk_nnode; l_res++)
    //   {
    //     //if(this->get_node_number(be_pt->node_pt(l_res)) == -1) continue;

    //     for(unsigned l_var=0; l_var<bulk_nnode; l_var++)
    //       {
    //         std::cout <<  "residual: " << l_res << ", variable: " << l_var << std::endl;

    //         for(unsigned j=0; j<3; j++)
    //           {
    //             for(unsigned i=0; i<3; i++)
    //               {
    //                 std::cout << std::setw( 12 ) << jacobian(l_res*3 + j, l_var*3 + i);
    //               }
    //             std::cout << std::endl;
    //           }
    //       }
    //   }

    //jacobian.output(std::cout);
  }

} // End of oomph namespace

#endif
