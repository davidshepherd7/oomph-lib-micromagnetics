#ifndef OOMPH_MICROMAGNETICS_ELEMENTS_HEADER
#define OOMPH_MICROMAGNETICS_ELEMENTS_HEADER


// Generic oomph-lib routines
#include "../../src/generic/Vector.h"
#include "../../src/generic/nodes.h"
#include "../../src/generic/Qelements.h"
#include "../../src/generic/Telements.h"
#include "../../src/generic/oomph_utilities.h"
#include "../../src/generic/oomph_definitions.h"

// Print vectors nicely
#include "./prettyprint98.hpp"

// My vector helpers
#include "./vector_helpers.h"
#include "./magnetic_materials.h"

// #include "./micromagnetics_flux_element.h"

// Magnetostatic elements are based on Poisson
#include "./template_free_poisson.h"

#include "./interpolator.h"


using namespace oomph;

namespace oomph
{

  // Forward declaration of flux element
  template <class ELEMENT> class MicromagFluxElement;

  class MMInterpolator;

  //==============================================================================
  /// A class for the maths used in solving the Landau-Lifshitz-Gilbert equations.
  //==============================================================================
  class MicromagEquations : public virtual FiniteElement
  {
  public:

    // CONSTRUCTORS ETC.
    /// Constructor (initialises the various function pointers to null).
    MicromagEquations() : Phi_source_pt(0), Phi_1_source_pt(0),
                          Magnetic_parameters_pt(0),
                          Applied_field_pt(0)
    {}

    /// Broken copy constructor
    MicromagEquations(const MicromagEquations& dummy)
    {BrokenCopy::broken_copy("MicromagEquations");}

    /// Broken assignment operator
    void operator=(const MicromagEquations&)
    {BrokenCopy::broken_assign("MicromagEquations");}

    /// Self-test: Return 0 for OK.
    unsigned self_test();

    // Equation numbering
    /// Specify nodal index for phi.
    unsigned phi_index_micromag() const {return 0;} // equation number 0

    /// Specify nodal index for phi 1
    unsigned phi_1_index_micromag() const {return 1;} //equation number 1

    /// Specify nodal index for kth component of M.
    unsigned m_index_micromag(const unsigned &k) const
    {
#ifdef PARANOID
      if(k>=3)
        throw OomphLibError("M only has 3 indices",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
#endif
      return 2 + k; // equations 2,3,4
    }

    /// Function determining how to block the Jacobian.
    void get_dof_numbers_for_unknowns(std::list<std::pair<unsigned long,unsigned> >&
                                      block_lookup_list)

    {
      // Clear list (just in case)
      block_lookup_list.clear();

      // Loop over all nodes then all unpinned values (dofs) at each node. For
      // each of these we create a pair giving the global equation number and
      // the corresponding dof type (number).
      for(unsigned nd=0; nd<nnode(); nd++)
        {
          // Put it into the block lookup list
          for(unsigned index = 0, nindex=node_pt(nd)->nvalue(); index<nindex; index++)
            {
              int local_eqn_number = this->nodal_local_eqn(nd,index);
              if(local_eqn_number >= 0)
                {
                  int global_eqn_number = eqn_number(local_eqn_number);
                  std::pair<unsigned long, unsigned> lookup;
                  lookup.first = global_eqn_number;
                  lookup.second = index;
                  block_lookup_list.push_front(lookup);
                }
            }
        }

    }

    void add_face_element_pt(FiniteElement* const face_element_pt)
    { Face_element_pts.insert(face_element_pt); }

    /// \short We need 5 values: 3 magnetisation + phi + phi1.
    unsigned required_nvalue(const unsigned &n) const
    {return 5;}

    unsigned ndof_types()
    {return required_nvalue(0);}

    typedef double (*TimeSpaceToDoubleFctPt)(const double& t, const Vector<double>&x);

    typedef Vector<double> (*TimeSpaceToDoubleVectFctPt)
    (const double& t, const Vector<double>&x);

    // typedef void (*TimeSpaceMagToDoubleVectFctPt)
    // (const double& t, const Vector<double>&x, const Vector<double>& M, Vector<double>& out);


    const MagneticParameters* magnetic_parameters_pt() const
    {return Magnetic_parameters_pt;}

    const MagneticParameters*& magnetic_parameters_pt()
    {return Magnetic_parameters_pt;}

    // SOURCE FUNCTIONS for testing
    /// Access function: Pointer to phi source function
    TimeSpaceToDoubleFctPt & phi_source_pt() {return Phi_source_pt;}

    /// Access function: Pointer to phi source function. Const version
    TimeSpaceToDoubleFctPt phi_source_pt() const {return Phi_source_pt;}

    /// Get phi source term at (Eulerian) position x.
    inline double get_phi_source(const double t,
                                 const Vector<double>& x) const
    {
      if(Phi_source_pt==0) return 0.0;
      else return (*Phi_source_pt)(t,x);
    }

    /// Access function: Pointer to phi_1 source function
    TimeSpaceToDoubleFctPt & phi_1_source_pt() {return Phi_1_source_pt;}

    /// Access function: Pointer to phi_1 source function. Const version
    TimeSpaceToDoubleFctPt phi_1_source_pt() const {return Phi_1_source_pt;}

    /// Get phi_1 source term at (Eulerian) position x.
    inline double get_phi_1_source(const double t,
                                   const Vector<double>& x) const
    {
      if(Phi_1_source_pt==0) return 0.0;
      else return (*Phi_1_source_pt)(t,x);
    }

    // /// Access function: Pointer to source function
    // TimeSpaceToDoubleVectFctPt& llg_source_pt() {return Llg_source_pt;}

    // /// Access function: Pointer to source function. Const version
    // TimeSpaceToDoubleVectFctPt llg_source_pt() const {return Llg_source_pt;}

    // /// Get LLG source term at (Eulerian) position x.
    // inline void get_source_llg(const double& t,
    //                          const Vector<double>& x,
    //                          Vector<double>& source) const
    // {
    //   if(Llg_source_pt==0) {for(unsigned j=0;j<3;j++) source[j] = 0.0;}
    //   else (*Llg_source_pt)(t,x,source);
    // }


    // APPLIED FIELD
    /// Access function: Pointer to applied field function
    TimeSpaceToDoubleVectFctPt& applied_field_pt() {return Applied_field_pt;}

    /// Access function: Pointer to applied field function. Const version
    TimeSpaceToDoubleVectFctPt applied_field_pt() const {return Applied_field_pt;}

    /// Get the applied field at Eulerian position x.
    virtual Vector<double> get_applied_field(const double& t, const Vector<double> &x,
                                             const Vector<double> &s) const
    {
      Vector<double> H_app(3,0.0);
      if(Applied_field_pt != 0)
        {
          H_app = (*Applied_field_pt)(t, x);
        }

      for(unsigned j=0; j<3; j++)
        {
          H_app[j] *= magnetic_parameters_pt()->field_normalisation_factor();
        }
      return H_app;
    }

    /// Get the crystalline anisotropy field at Eulerian position x.
    inline Vector<double> get_H_cryst_anis_field(const double& t,
                                                 const Vector<double> &x,
                                                 const Vector<double>& m) const
    {
      Vector<double> h_ca(3, 0.0);
      magnetic_parameters_pt()->
        crystalline_ansiotropy_field(t, x, m, h_ca);
      return h_ca;
    }

    void get_hca_derivative(const double& t, const Vector<double>&x,
                            const Vector<double>& m,
                            const double shape_fn_k_at_x,
                            DenseMatrix<double>& dhcadm) const
    {
      magnetic_parameters_pt()->
        crystalline_ansiotropy_field_derivative(t,x,m,shape_fn_k_at_x,dhcadm);
    }

    /// Get saturisation magnetisation at eulerian postition x.
    inline double sat_mag() const
    {
      //??ds this isn't included!
      return magnetic_parameters_pt()->normalised_saturation_magnetisation();
    }

    /// Get LLG damping coefficient.
    inline double llg_damping_coeff() const
    {
      return magnetic_parameters_pt()->normalised_gilbert_damping();
    }

    /// Get LLG precession coefficient.
    inline double llg_precession_coeff() const
    {
      return magnetic_parameters_pt()->normalised_gamma();
    }

    /// Get exchange coefficient at eulerian postition x.
    inline double exchange_coeff() const
    {
      return magnetic_parameters_pt()->normalised_hex();
    }

    /// Get magnetostatic coefficient at eulerian postition x.
    inline double magnetostatic_coeff() const
    {
      return magnetic_parameters_pt()->normalised_hms();
    }

    // // EXACT PHI FUNCTION POINTER
    // /// Access function: Pointer to exact phi function
    // TimeSpaceToDoubleFctPt& exact_phi_pt() {return Exact_phi_pt;}

    // /// Access function: Pointer to exact phi function. Const version
    // TimeSpaceToDoubleFctPt exact_phi_pt() const {return Exact_phi_pt;}

    // /// Get exact phi at eulerian postition x.
    // inline double get_exact_phi(const double& t, const Vector<double>& x) const
    // {
    //   // If no exact phi function has been set, return something crazy
    //   if(Exact_phi_pt==0) {return -1000.0;}
    //   else return (*Exact_phi_pt)(t,x);
    // }

    // // EXACT M FUNCTION POINTER
    // /// Access function: Pointer to exact M function
    // TimeSpaceToDoubleVectFctPt& exact_m_pt() {return Exact_m_pt;}

    // /// Access function: Pointer to LLG exact M. Const version
    // TimeSpaceToDoubleVectFctPt exact_m_pt() const {return Exact_m_pt;}

    // /// Get exact M at eulerian postition x.
    // inline void get_exact_m(const double& t, const Vector<double>& x,
    //                       Vector<double>& m_exact) const
    // {
    //   // If no exact M function has been set empty the vector
    //   if(Exact_m_pt==0) m_exact.clear();
    //   else (*Exact_m_pt)(t,x,m_exact);
    // }


    /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
    inline double dshape_dtest(const Vector<double> &s,
                               Shape &psi, DShape &dpsidx,
                               Shape &test, DShape &dtestdx) const
    {
      // Call the geometrical shape functions and derivatives
      const double J = this->dshape_eulerian(s,psi,dpsidx);

      // Set the test functions equal to the shape functions
      test = psi;
      dtestdx= dpsidx;

      // Return the jacobian
      return J;
    }

    /// Return FE representation of function value phi(s) at local coordinate s
    inline double interpolated_phi_micromag(const Vector<double> &s) const
    {
      //Find number of nodes
      const unsigned n_node = nnode();

      //Get the index at which the poisson unknown is stored
      const unsigned phi_nodal_index = phi_index_micromag();

      //Local shape function
      Shape psi(n_node);

      //Find values of shape function
      shape(s,psi);

      //Initialise value of u
      double itp_phi = 0.0;

      //Loop over the local nodes and sum
      for(unsigned l=0;l<n_node;l++)
        {
          itp_phi += this->nodal_value(l,phi_nodal_index)*psi[l];
        }

      return(itp_phi);
    }


    /// Return FE representation of phi_1 at local coordinate s
    inline double interpolated_phi_1_micromag(const Vector<double> &s) const
    {
      //Find number of nodes
      const unsigned n_node = nnode();

      //Get the index at which the poisson unknown is stored
      const unsigned phi_1_nodal_index = phi_1_index_micromag();

      //Local shape function
      Shape psi(n_node);

      //Find values of shape function
      shape(s,psi);

      //Initialise value of u
      double itp_phi_1 = 0.0;

      //Loop over the local nodes and sum
      for(unsigned l=0;l<n_node;l++)
        {
          itp_phi_1 += this->nodal_value(l,phi_1_nodal_index)*psi[l];
        }

      return itp_phi_1;
    }


    /// Return FE representation of M at local coordinate s and current time.
    inline void interpolated_m_micromag(const Vector<double> &s,
                                        Vector<double>& itp_m) const
    {
      // Find number of nodes
      const unsigned n_node = nnode();

      // Get local shape function
      Shape psi(n_node);
      shape(s,psi);

      // Interpolate m
      itp_m.assign(3,0.0);
      for(unsigned k=0; k<3; k++)
        {
          // For each component sum over the contributions due to each node
          // in this element.
          for(unsigned l=0;l<n_node;l++)
            {
              itp_m[k] += this->nodal_value(l, m_index_micromag(k))*psi[l];
            }
        }
    }

    /// Return FE representation of M at local coordinate s and current time.
    void interpolated_dmdx_micromag(const Vector<double> &s,
                                    Vector<Vector<double> >& itp_dmdx) const;


    /// \short Return FE representation of solution vector (phis,M,H_ex)
    /// at local coordinate s and current time.
    inline void interpolated_solution_micromag(const Vector<double> &s,
                                               Vector<double>& itp_solution) const
    {
      //Find number of nodes
      const unsigned n_node = nnode();

      // Get number of values required
      const unsigned nvalue = required_nvalue(0);

      //Local shape function
      Shape psi(n_node);

      //Find values of shape function
      shape(s,psi);

      // Initialise solution vector
      itp_solution.assign(nvalue,0.0);

      // Loop over the list of solutions
      for(unsigned i=0; i<nvalue; i++)
        {
          //Loop over the local nodes and sum
          for(unsigned l=0;l<n_node;l++)
            {
              itp_solution[i] += this->nodal_value(l,i)*psi[l];
            }
        }

    }

    // /// Get total field at position x. There is a more efficient private version
    // /// of this function for use in residual calculations etc. This version does
    // /// some extra calculations to get interpolated m and x and ensures
    // /// initilisations.
    // //??ds unchecked
    // void interpolated_ht_micromag(const Vector<double>& s,
    //                               Vector<double>& h_total)
    // {
    //   Vector<double> x(3,0.0);
    //   interpolated_x(s,x);

    //   Vector<double> m(3,0.0);
    //   interpolated_m_micromag(s,m);

    //   Vector<double> itp_dphidx(3,0.0);
    //   interpolated_dphidx_micromag(s,itp_dphidx);

    //   h_total.assign(3,0.0);
    //   interpolated_ht_micromag_efficient(x,m,itp_dphidx,h_total);
    // }

    /// Get divergence of magnetisation (for poisson source).
    double divergence_m(const unsigned &ipt)
    {
      // Get values
      const unsigned n_node = nnode();
      Shape psi(n_node);  DShape dpsidx(n_node, this->nodal_dimension());
      dshape_eulerian_at_knot(ipt,psi,dpsidx);

      // Sum up the derivatives
      double div_m = 0.0;
      for(unsigned j=0; j<nodal_dimension(); j++)
        {
          for(unsigned l=0; l<n_node; l++)
            {
              div_m += nodal_value(l,m_index_micromag(j)) * dpsidx(l,j);
            }
        }

      return div_m;
    }



    // OUTPUT FUNCTIONS
    /// Output FE representation of soln: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(std::ostream &outfile, const unsigned &n_plot=5);

    /// Output exact solution at n_plot points
    void output_fct(std::ostream &outfile, const unsigned &n_plot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt);

    /// Get error by comparing with exact solution and get norm of exact solution.
    void compute_error(std::ostream &outfile,
                       FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                       const double& time, double& error, double& norm);

    // RESIDUALS + JACOBIAN
    /// Add the element's contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      //Call the generic residuals function with flag set to 0 using a dummy matrix argument
      fill_in_generic_residual_contribution_micromag
        (residuals,GeneralisedElement::Dummy_matrix, 0);
    }

    /// \short Add the element's contribution to its residual vector and element
    /// Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                          DenseMatrix<double> &jacobian)
    {
      // Call the generic routine with the flag set to 1

      //??ds debug code
      //FiniteElement::fill_in_contribution_to_jacobian(residuals,jacobian);

      fill_in_face_element_contribution_to_jacobian(jacobian);

      fill_in_generic_residual_contribution_micromag(residuals,jacobian,1);
    }

    virtual void fill_in_face_element_contribution_to_jacobian
    (DenseMatrix<double> &jacobian) const =0;

    /// Put dM/dt at local node l (using timestepper and history values) in the vector dmdt.
    void dm_dt_micromag(const unsigned& l, Vector<double>& dmdt) const
    {
      // Get the data's timestepper
      TimeStepper* time_stepper_pt= this->node_pt(l)->time_stepper_pt();

      //Initialise dM/dt to zero
      dmdt.assign(3,0.0);

      //Loop over the timesteps, if there is a non Steady timestepper
      if (!time_stepper_pt->is_steady())
        {
          // Get number of timsteps to use (past & present)
          const unsigned n_time_steps = time_stepper_pt->ntstorage();

          // Loop over the directions
          for (unsigned j=0; j<3; j++)
            {
              // Loop over past and present times and add the contributions to
              // the time derivative
              for(unsigned t=0;t<n_time_steps;t++)
                {
                  // The "1" in weights is the derivative order (i.e. 1:
                  // dM/dt, 2: d^2M/dt^2)
                  dmdt[j] += time_stepper_pt->weight(1,t)
                    *nodal_value(t,l,m_index_micromag(j));
                }
            } // End of loop over directions
        }

    } // end of dM_dt_micromag


    /// \short Return a vector containing the magnetisation at a node.
    Vector<double> get_m(unsigned node) const
    {
      Vector<double> m(3, 0.0);
      for(unsigned j=0; j<3; j++)
        {
          m[j] = nodal_value(node, m_index_micromag(j));
        }
      return m;
    }


    /// \short Get the maximum difference in angle between the
    /// magnetisation of two nodes of the element. If this is large there
    /// is likely an error somewhere (in code or in problem setup), or
    /// there is not enough refinement.
    double max_m_angle_variation() const
    {
      double max_angle = 0;

      // Double loop over nodes: compare magnetisation angles at all nodes
      // to each other, store the maximum difference of angles.
      for(unsigned nd=0, n_nd=nnode(); nd<n_nd; nd++)
        {
          Vector<double> m1 = get_m(nd);
          for(unsigned nd2=0; nd<n_nd; nd++)
            {
              Vector<double> m2 = get_m(nd2);
              max_angle = std::max(max_angle,
                                   VectorOps::angle_diff(m1, m2));
            }
        }
      return max_angle;
    }


    /// A dummy double to hold space in Jacobian matrices for entries that are
    /// determined by the boundary element method part of the hybrid method
    /// (until it can be filled in at the problem level).
    static const double DummyBEMControlledEntry;

  protected:

    /// Fill in contribution to residuals and jacobian (if flag is set) from
    /// these equations (compatible with multiphysics)
    void fill_in_generic_residual_contribution_micromag(Vector<double> &residuals,
                                                        DenseMatrix<double> &jacobian,
                                                        const unsigned& flag) const;

    /// In PARANOID compare new and old interpolation functions to check
    /// they are the same, else do nothing.
    void check_interpolation(MMInterpolator& intp) const;

    // /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
    // virtual double dshape_dtest(const Vector<double> &s, Shape &psi, DShape &dpsidx,
    //                             Shape &test, DShape &dtestdx) const=0;

    // /// Get the field at current time (pass in x,m, dphidx for
    // /// efficiency).
    // inline void interpolated_ht_micromag_efficient(const Vector<double>& x,
    //                                                const Vector<double>& m,
    //                                                const Vector<double>& itp_dphidx,
    //                                                Vector<double>& h_total)
    //   const
    // {
    //   // Get time
    //   const double time = time_pt()->time();

    //   Vector<double> h_applied(3,0.0);
    //   get_applied_field(time, x,  h_applied);
    //   Vector<double> h_cryst_anis(3,0.0);
    //   get_H_cryst_anis_field(time, x, m, h_cryst_anis);

    //   // Magnetostatic field is -1* dphi/dx, multiply by a coeff too alow easy
    //   // switching on and off for debugging.
    //   double magnetosttaic_coeff = magnetostatic_coeff(time,x);
    //   Vector<double> h_ms(3,0.0);
    //   for(unsigned i=0; i<h_ms.size(); i++)
    //     h_ms[i] *= -1 * magnetostatic_coeff;

    //   // Take total of all fields used
    //   for(unsigned j=0; j<3; j++)
    //     {
    //       h_total[j] = h_applied[j] + h_ms[j] + h_cryst_anis[j];
    //     }
    // }

    /// Pointer to poisson source function - only for testing purposes since
    /// div(M) is our source function in calculation of the demagnetising
    /// potential.
    TimeSpaceToDoubleFctPt Phi_source_pt;

    TimeSpaceToDoubleFctPt Phi_1_source_pt;

    const MagneticParameters* Magnetic_parameters_pt;

    /// Pointer to function giving applied field.
    TimeSpaceToDoubleVectFctPt Applied_field_pt;

    // /// Pointer to the exact solution for M
    // TimeSpaceToDoubleVectFctPt Exact_m_pt;

    // /// Pointer to the exact solution for phi
    // TimeSpaceToDoubleFctPt Exact_phi_pt;

    // List of face elements attached to this element
    std::set<FiniteElement*> Face_element_pts;

  }; // End of MicromagEquations class


  //====================================================================
  /// A class combining the micromag equations with a QElement geometry
  //====================================================================
  template < unsigned DIM, unsigned NNODE_1D>
  class QMicromagElement : public QElement<DIM,NNODE_1D>, public MicromagEquations
  {
  public:
    /// Constructor
    QMicromagElement() {}

    /// Broken copy constructor
    QMicromagElement(const QMicromagElement<DIM,NNODE_1D>& dummy)
    {
      BrokenCopy::broken_copy("QMicromagElement");
    }

    /// Broken assignment operator
    void operator=(const QMicromagElement<DIM,NNODE_1D>&)
    {
      BrokenCopy::broken_assign("QMicromagElement");
    }

    /// Output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(std::ostream &outfile, const unsigned &n_plot=5)
    {MicromagEquations::output(outfile,n_plot);}

    // /// C-style output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    // void output(FILE* file_pt, const unsigned &n_plot = 5)
    // {MicromagEquations::output(file_pt,n_plot);}

    // /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
    // inline double dshape_dtest(const Vector<double> &s,
    //                            Shape &psi, DShape &dpsidx,
    //                            Shape &test, DShape &dtestdx) const
    // {
    //   // Call the geometrical shape functions and derivatives
    //   const double J = this->dshape_eulerian(s,psi,dpsidx);

    //   // Set the test functions equal to the shape functions
    //   test = psi;
    //   dtestdx= dpsidx;

    //   // Return the jacobian
    //   return J;
    // }

    void fill_in_face_element_contribution_to_jacobian
    (DenseMatrix<double> &jacobian) const
    {
      std::set<FiniteElement*>::iterator it;
      for(it=this->Face_element_pts.begin(); it!=this->Face_element_pts.end(); it++)
        {
          MicromagFluxElement<QMicromagElement<DIM,NNODE_1D> >* flux_ele_pt =
            dynamic_cast<MicromagFluxElement<QMicromagElement<DIM,NNODE_1D> >* >
            (*it);
          flux_ele_pt->fill_in_bulk_contribution_to_face_jacobian(jacobian);
        }
    }

  }; // end of QMicromagElement class declaration



  //====================================================================
  /// A class combining the micromag equations with a TElement geometry
  //====================================================================
  template < unsigned DIM, unsigned NNODE_1D>
  class TMicromagElement : public TElement<DIM,NNODE_1D>, public MicromagEquations
  {
  public:
    /// Constructor
    TMicromagElement() {}

    /// Broken copy constructor
    TMicromagElement(const TMicromagElement<DIM,NNODE_1D>& dummy)
    {
      BrokenCopy::broken_copy("TMicromagElement");
    }

    /// Broken assignment operator
    void operator=(const TMicromagElement<DIM,NNODE_1D>&)
    {
      BrokenCopy::broken_assign("TMicromagElement");
    }

    /// Output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(std::ostream &outfile, const unsigned &n_plot=5)
    {MicromagEquations::output(outfile,n_plot);}

    // /// C-style output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    // void output(FILE* file_pt, const unsigned &n_plot = 5)
    // {MicromagEquations::output(file_pt,n_plot);}

    void fill_in_face_element_contribution_to_jacobian
    (DenseMatrix<double> &jacobian) const
    {
      std::set<FiniteElement*>::iterator it;
      for(it=this->Face_element_pts.begin(); it!=this->Face_element_pts.end(); it++)
        {
          MicromagFluxElement<TMicromagElement<DIM,NNODE_1D> >* flux_ele_pt =
            dynamic_cast<MicromagFluxElement<TMicromagElement<DIM,NNODE_1D> >* >
            (*it);
          flux_ele_pt->fill_in_bulk_contribution_to_face_jacobian(jacobian);
        }
    }

  }; // end of TMicromagElement class declaration


  // =================================================================
  ///
  //??ds if you wnt to use phi_1 + phi_2 style calculation might have to
  // change things here.
  // =================================================================
  class MagnetostaticFieldEquations : public TFPoissonEquations
  {
  public:

    /// Constructor (null the pointers)
    MagnetostaticFieldEquations() : Micromag_element_pt(0) {}

    /// Destructor
    ~MagnetostaticFieldEquations() {}

    /// Get the magnetostatic field at local coordinate point s in the element.
    void magnetostatic_field(const Vector<double> &s, Vector<double> &hms) const
    {
      // Get values
      const unsigned n_node = this->nnode();
      const unsigned u_nodal_index = this->u_index_poisson();
      Shape psi(n_node); DShape dpsidx(n_node, nodal_dimension());
      this->dshape_eulerian(s,psi,dpsidx);

      // Interpolate
      hms.assign(3,0.0);
      for(unsigned j=0; j<nodal_dimension(); j++)
        {
          for(unsigned l=0; l<n_node; l++)
            {
              hms[j] -= this->nodal_value(l,u_nodal_index)*dpsidx(l,j);
            }

          // Normalise
          hms[j] *= Micromag_element_pt->magnetostatic_coeff();
        }
    }

    /// For micromagnetics the source function is the divergence of the
    /// magnetisation.
    void get_source_poisson(const unsigned& ipt,
                            const Vector<double>& x,
                            double& source) const
    {
      source = 0;
      // Get contribution from divergence of M at this integration point.
      source += Micromag_element_pt->divergence_m(ipt);

      // Get contribution from any real source functions.
      double poisson_source=0;
      TFPoissonEquations::get_source_poisson(ipt, x, poisson_source);
      source += poisson_source;
    }

    // Access functions
    // ============================================================

    /// \short Non-const access function for Micromag_element_pt.
    MicromagEquations* micromag_element_pt() const
    {return micromag_element_pt();}

    void set_micromag_element_pt(MicromagEquations* ele_pt)
    {Micromag_element_pt = ele_pt;}

    /// \short Const access function for Micromag_element_pt.
    MicromagEquations*& micromag_element_pt()
    {
      // Lots of checks because this is stupid really...
#ifdef PARANOID
      if(Micromag_element_pt == 0)
        {
          std::ostringstream error_msg;
          error_msg << "Magnetics element pointer not set.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      if(this->nnode() != Micromag_element_pt->nnode())
        {
          std::ostringstream error_msg;
          error_msg << "Elements must be the same geometry for this to "
                    << "work... sorry for the hackyness. Maybe you can fix it.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      if(this->dim() != Micromag_element_pt->dim())
        {
          std::ostringstream error_msg;
          error_msg << "Elements must be the same geometry for this to "
                    << "work... sorry for the hackyness. Maybe you can fix it.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      if(this->integral_pt() != Micromag_element_pt->integral_pt())
        {
          std::ostringstream error_msg;
          error_msg << "Elements must have the same integration scheme for this to"
                    << "work... sorry for the hackyness. Maybe you can fix it.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      // Check node positions
      for(unsigned nd=0, n_nd=this->nnode(); nd<n_nd; nd++)
        {
          for(unsigned j=0; j<this->node_pt(nd)->ndim(); j++)
            {
              if(this->node_pt(nd)->position(j)
                 != Micromag_element_pt->node_pt(nd)->position(j))
                {
                  std::ostringstream error_msg;
                  error_msg << "Mismatch in positions.";
                  throw OomphLibError(error_msg.str(),
                                      OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
                }
            }
        }
#endif
      return Micromag_element_pt;
    }

  private:

    /// Pointer to the element from which it should get the divergence of
    /// the magnetisation at a point (for use as source function).
    MicromagEquations* Micromag_element_pt;

  };

  // =================================================================
  ///
  // =================================================================
  template < unsigned DIM, unsigned NNODE_1D>
  class TMagnetostaticFieldElement : public TElement<DIM,NNODE_1D>,
                                     public MagnetostaticFieldEquations
  {
  public:

    // Overload output functions (diamond inheritance...)
    // ============================================================

    /// Output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(std::ostream &outfile, const unsigned &n_plot = 5)
    {TFPoissonEquations::output(outfile,n_plot);}

    /// Output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(std::ostream &outfile)
    {TFPoissonEquations::output(outfile,5);}

    /// C-style output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(FILE* outfile_pt, const unsigned &n_plot = 5)
    {TFPoissonEquations::output(outfile_pt, n_plot);}

    /// C-style output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(FILE* outfile_pt)
    {TFPoissonEquations::output(outfile_pt, 5);}

    // Copied from the Poisson versions of these functions:
    // ============================================================

    double dshape_and_dtest_eulerian_poisson
    (const Vector<double> &s, Shape &psi, DShape &dpsidx, Shape &test,
     DShape &dtestdx) const
    {
      const double J = this->dshape_eulerian(s,psi,dpsidx);
      test = psi;
      dtestdx= dpsidx;
      return J;
    }

    double dshape_and_dtest_eulerian_at_knot_poisson
    (const unsigned &ipt, Shape &psi, DShape &dpsidx, Shape &test,
     DShape &dtestdx) const
    {
      const double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);
      test = psi;
      dtestdx = dpsidx;
      return J;
    }

    double dshape_and_dtest_eulerian_at_knot_poisson
    (const unsigned &ipt, Shape &psi, DShape &dpsidx,
     RankFourTensor<double> &d_dpsidx_dX, Shape &test, DShape &dtestdx,
     RankFourTensor<double> &d_dtestdx_dX, DenseMatrix<double> &djacobian_dX) const
    {
      const double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx,
                                                     djacobian_dX,d_dpsidx_dX);
      test = psi;
      dtestdx = dpsidx;
      d_dtestdx_dX = d_dpsidx_dX;
      return J;
    }

  };

  // =================================================================
  ///
  // =================================================================
  template < unsigned DIM, unsigned NNODE_1D>
  class QMagnetostaticFieldElement : public QElement<DIM,NNODE_1D>,
                                     public MagnetostaticFieldEquations
  {
  public:

    // Overload output functions (diamond inheritance...)
    // ============================================================

    /// Output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(std::ostream &outfile, const unsigned &n_plot = 5)
    {TFPoissonEquations::output(outfile,n_plot);}

    /// Output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(std::ostream &outfile)
    {TFPoissonEquations::output(outfile,5);}

    /// C-style output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(FILE* outfile_pt, const unsigned &n_plot = 5)
    {TFPoissonEquations::output(outfile_pt, n_plot);}

    /// C-style output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(FILE* outfile_pt)
    {TFPoissonEquations::output(outfile_pt, 5);}

    // Copied from the Poisson versions of these functions:
    // ============================================================

    double dshape_and_dtest_eulerian_poisson
    (const Vector<double> &s, Shape &psi, DShape &dpsidx, Shape &test,
     DShape &dtestdx) const
    {
      const double J = this->dshape_eulerian(s,psi,dpsidx);
      test = psi;
      dtestdx= dpsidx;
      return J;
    }

    double dshape_and_dtest_eulerian_at_knot_poisson
    (const unsigned &ipt, Shape &psi, DShape &dpsidx, Shape &test,
     DShape &dtestdx) const
    {
      const double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);
      test = psi;
      dtestdx = dpsidx;
      return J;
    }

    double dshape_and_dtest_eulerian_at_knot_poisson
    (const unsigned &ipt, Shape &psi, DShape &dpsidx,
     RankFourTensor<double> &d_dpsidx_dX, Shape &test, DShape &dtestdx,
     RankFourTensor<double> &d_dtestdx_dX, DenseMatrix<double> &djacobian_dX) const
    {
      const double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx,
                                                     djacobian_dX,d_dpsidx_dX);
      test = psi;
      dtestdx = dpsidx;
      d_dtestdx_dX = d_dpsidx_dX;
      return J;
    }

  };



  // =================================================================
  /// Micromagnetics elements with additional coupling from magnetostatic
  /// field elements via a pointer (for use in semi-implicit methods).
  // =================================================================
  class SemiImplicitMicromagEquations : public MicromagEquations
  {
  public:

    /// For micromagnetics the source function is the divergence of the
    /// magnetisation.
    Vector<double> get_applied_field(const double& t,
                                     const Vector<double> &x,
                                     const Vector<double> &s) const
    {
      // Lots of checks because this is stupid really...
#ifdef PARANOID
      if(this->nnode() != magnetostatic_field_element_pt()->nnode())
        {
          std::ostringstream error_msg;
          error_msg << "Elements must be the same geometry for this to "
                    << "work... sorry for the hackyness. Maybe you can fix it.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      if(this->dim() != magnetostatic_field_element_pt()->dim())
        {
          std::ostringstream error_msg;
          error_msg << "Elements must be the same geometry for this to "
                    << "work... sorry for the hackyness. Maybe you can fix it.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      if(this->integral_pt() != magnetostatic_field_element_pt()->integral_pt())
        {
          std::ostringstream error_msg;
          error_msg << "Elements must have the same integration scheme for this to"
                    << "work... sorry for the hackyness. Maybe you can fix it.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      // Get contribution from field element
      Vector<double> h_ms;
      magnetostatic_field_element_pt()->magnetostatic_field(s, h_ms);

      // Get contribution from any real applied field functions.
      Vector<double> h_app = MicromagEquations::get_applied_field(t, x, s);

      //Add them up
      for(unsigned j=0; j<3; j++) h_app[j] += h_ms[j];

      return h_app;
    }

    // Access functions:
    // ============================================================

    /// \short Non-const access function for Magnetostatic_field_element_pt.
    MagnetostaticFieldEquations*& magnetostatic_field_element_pt()
    {return Magnetostatic_field_element_pt;}

    /// \short Const access function for Magnetostatic_field_element_pt.
    MagnetostaticFieldEquations* magnetostatic_field_element_pt() const
    {
#ifdef PARANOID
      if(Magnetostatic_field_element_pt == 0)
        {
          std::ostringstream error_msg;
          error_msg << "Magnetics element pointer not set.";
          throw OomphLibError(error_msg.str(), OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return Magnetostatic_field_element_pt;
    }

  private:

    MagnetostaticFieldEquations* Magnetostatic_field_element_pt;
  };


  //====================================================================
  /// A class combining the micromag equations with a QElement geometry
  //====================================================================
  template < unsigned DIM, unsigned NNODE_1D>
  class QSemiImplicitMicromagElement : public QElement<DIM,NNODE_1D>,
                                       public SemiImplicitMicromagEquations
  {
  public:

    /// Output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(std::ostream &outfile, const unsigned &n_plot=5)
    {SemiImplicitMicromagEquations::output(outfile,n_plot);}

    // /// C-style output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    // void output(FILE* file_pt, const unsigned &n_plot = 5)
    // {SemiImplicitMicromagEquations::output(file_pt,n_plot);}

    // /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
    // inline double dshape_dtest(const Vector<double> &s,
    //                            Shape &psi, DShape &dpsidx,
    //                            Shape &test, DShape &dtestdx) const
    // {
    //   // Call the geometrical shape functions and derivatives
    //   const double J = this->dshape_eulerian(s,psi,dpsidx);

    //   // Set the test functions equal to the shape functions
    //   test = psi;
    //   dtestdx= dpsidx;

    //   // Return the jacobian
    //   return J;
    // }

    //??ds
    void fill_in_face_element_contribution_to_jacobian
    (DenseMatrix<double> &jacobian) const
    {
      std::set<FiniteElement*>::iterator it;
      for(it=this->Face_element_pts.begin(); it!=this->Face_element_pts.end(); it++)
        {
          MicromagFluxElement<QSemiImplicitMicromagElement<DIM,NNODE_1D> >* flux_ele_pt =
            dynamic_cast<MicromagFluxElement<QSemiImplicitMicromagElement<DIM,NNODE_1D> >* >
            (*it);
          flux_ele_pt->fill_in_bulk_contribution_to_face_jacobian(jacobian);
        }
    }


  }; // end of QSemiImplicitMicromagElement class declaration



  //====================================================================
  /// A class combining the micromag equations with a TElement geometry
  //====================================================================
  template < unsigned DIM, unsigned NNODE_1D>
  class TSemiImplicitMicromagElement : public TElement<DIM,NNODE_1D>,
                                       public SemiImplicitMicromagEquations
  {
  public:

    /// Output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(std::ostream &outfile, const unsigned &n_plot=5)
    {SemiImplicitMicromagEquations::output(outfile,n_plot);}

    // /// C-style output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    // void output(FILE* file_pt, const unsigned &n_plot = 5)
    // {SemiImplicitMicromagEquations::output(file_pt,n_plot);}

    // /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
    // inline double dshape_dtest(const Vector<double> &s,
    //                            Shape &psi, DShape &dpsidx,
    //                            Shape &test, DShape &dtestdx) const
    // {
    //   // Call the geometrical shape functions and derivatives
    //   const double J = this->dshape_eulerian(s,psi,dpsidx);

    //   // Set the test functions equal to the shape functions
    //   test = psi;
    //   dtestdx= dpsidx;

    //   // Return the jacobian
    //   return J;
    // }

    //??ds
    void fill_in_face_element_contribution_to_jacobian
    (DenseMatrix<double> &jacobian) const
    {
      std::set<FiniteElement*>::iterator it;
      for(it=this->Face_element_pts.begin(); it!=this->Face_element_pts.end(); it++)
        {
          MicromagFluxElement<TSemiImplicitMicromagElement<DIM,NNODE_1D> >* flux_ele_pt =
            dynamic_cast<MicromagFluxElement<TSemiImplicitMicromagElement<DIM, NNODE_1D> >* >
            (*it);

#ifdef PARANOID
          if(flux_ele_pt == 0)
            {
              std::ostringstream error_msg;
              error_msg << "Failed dynamic cast.";
              throw OomphLibError(error_msg.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
#endif


          flux_ele_pt->fill_in_bulk_contribution_to_face_jacobian(jacobian);
        }
    }

  };


  class MMInterpolator : public GeneralInterpolator
  {
    // Assumption: m_index_micromag(0 - 3) are consecutive.

  private:

    // Extra storage for magnetisation values, so we can have nice vector
    // access to them.
    Vector<double> Dmdt;
    Vector<double> M;
    Vector<Vector<double> > Dmdx;
    double Div_m;

    const MicromagEquations* This_element;

  public:

    /// Default constructor
    MMInterpolator(const FiniteElement* const this_element,
                   const Vector<double> &s)
      : GeneralInterpolator(this_element, s),
        Div_m(this->NotYetCalculatedValue)
    {}

    double phi() {return this->value(This_element->phi_index_micromag());}
    const Vector<double> & dphidx()
    {return this->dvaluedx(This_element->phi_index_micromag());}

    double phi1() {return this->value(This_element->phi_1_index_micromag());}
    const Vector<double> & dphi1dx()
    {return this->dvaluedx(This_element->phi_1_index_micromag());}

    const Vector<double> &m()
    {
      if(this->uninitialised(M))
        {
          M = this->interpolate_values(This_element->m_index_micromag(0),
                                       This_element->m_index_micromag(2) + 1);
        }
      return M;
    }

    const Vector<double>& dmdt()
    {
      if(this->uninitialised(Dmdt))
        {
          Dmdt = this->interpolate_dvaluesdt(This_element->m_index_micromag(0),
                                             This_element->m_index_micromag(2) + 1);
        }
      return Dmdt;
    }

    const Vector<Vector<double> >& dmdx()
    {
      if(this->uninitialised(Dmdx))
        {
          Dmdx = this->interpolate_dvaluesdx(This_element->m_index_micromag(0),
                                             This_element->m_index_micromag(2) + 1);
        }
      return Dmdx;
    }

    double div_m()
    {
      if(this->uninitialised(Div_m))
        {
          Div_m = 0;
          for(unsigned j=0; j<Dim; j++) Div_m += dmdx()[j][j];
        }
      return Div_m;
    }
  };


  // =================================================================
  /// Face geometry elements for the elements defined above.
  // =================================================================

  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<TMagnetostaticFieldElement<DIM,NNODE_1D> >:
    public virtual TElement<DIM-1,NNODE_1D> {};

  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<QMagnetostaticFieldElement<DIM,NNODE_1D> >:
    public virtual QElement<DIM-1,NNODE_1D> {};

  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<TMicromagElement<DIM,NNODE_1D> >:
    public virtual TElement<DIM-1,NNODE_1D> {};

  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<QMicromagElement<DIM,NNODE_1D> >:
    public virtual QElement<DIM-1,NNODE_1D> {};

  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<TSemiImplicitMicromagElement<DIM,NNODE_1D> >:
    public virtual TElement<DIM-1,NNODE_1D> {};

  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<QSemiImplicitMicromagElement<DIM,NNODE_1D> >:
    public virtual QElement<DIM-1,NNODE_1D> {};
}




#endif
