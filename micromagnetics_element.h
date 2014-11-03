#ifndef OOMPH_MICROMAGNETICS_ELEMENTS_HEADER
#define OOMPH_MICROMAGNETICS_ELEMENTS_HEADER


// Generic oomph-lib routines
#include "../../src/generic/Vector.h"
#include "../../src/generic/nodes.h"
#include "../../src/generic/Qelements.h"
#include "../../src/generic/Telements.h"
#include "../../src/generic/refineable_elements.h"
#include "../../src/generic/refineable_brick_element.h"
#include "../../src/generic/refineable_quad_element.h"
#include "../../src/generic/refineable_line_element.h"
#include "../../src/generic/oomph_utilities.h"
#include "../../src/generic/oomph_definitions.h"

// This is only used for casting for an error check, get rid of it? Put in
// .cc?
#include "../../src/generic/midpoint_method.h"

// My vector helpers
#include "./vector_helpers.h"
#include "./magnetic_parameters.h"

// Magnetostatic elements are based on Poisson
#include "./template_free_poisson.h"

#include "./energy_functions.h"
#include "micromag_types.h"
#include "./residual_calculator.h"
#include "magnetostatics_calculator.h"


#include "array_interpolator.h"
#include "interpolator.h"
#include "new_interpolators.h"




namespace oomph
{

  // Forward declaration of flux element
  template <class ELEMENT> class MicromagFluxElement;

  class MMInterpolator;
  class CachingMMArrayInterpolator;

  //==============================================================================
  /// A class for the maths used in solving the Landau-Lifshitz-Gilbert equations.
  //==============================================================================
  class MicromagEquations : public virtual FiniteElement
  {
  public:

    // CONSTRUCTORS ETC.
    /// Constructor (initialises the various function pointers to null).
    MicromagEquations() : Use_fd_jacobian(false),
                          Residual_calculator_pt(0),
                          Phi_source_pt(0), Phi_1_source_pt(0),
                          Magnetic_parameters_pt(0)
    {
      Ms_calc_pt = 0;
    }

    /// Virtual destructor
    virtual ~MicromagEquations()
    {
      // easiest to just always have own calculator for now...
      delete Ms_calc_pt; Ms_calc_pt = 0;
    }

    /// Broken copy constructor
    MicromagEquations(const MicromagEquations& dummy)
    {BrokenCopy::broken_copy("MicromagEquations");}

    /// Broken assignment operator
    void operator=(const MicromagEquations&)
    {BrokenCopy::broken_assign("MicromagEquations");}

    void add_face_element_pt(FiniteElement* const face_element_pt)
    { Face_element_pts.insert(face_element_pt); }

    /// Self-test: Return 0 for OK.
    unsigned self_test()
    {
#ifdef PARANOID
      // Check that M indicies are sequential. Otherwise interpolation will
      // break...
      bool ok = ((m_index_micromag(0) +1) == m_index_micromag(1));
      ok = ok && ((m_index_micromag(1) +1) == m_index_micromag(2));
      if(!ok)
        {
          std::string error_msg = "M indicies must be sequential!";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      return 0;
    }

    // Equation numbering
    // ============================================================

    /// \short We need 5 values: 3 magnetisation + phi + phi1.
    unsigned required_nvalue(const unsigned &n) const {return 5;}

    /// Specify nodal index for phi.
    unsigned phi_index_micromag() const {return 0;} // equation number 0

    /// Specify nodal index for phi 1
    unsigned phi_1_index_micromag() const {return 1;} //equation number 1

    /// Specify nodal index for kth component of M.
    unsigned m_index_micromag(const unsigned &k) const
    {
#ifdef PARANOID
      if(k>=3) throw OomphLibError("M only has 3 indices",
                                   OOMPH_CURRENT_FUNCTION,
                                   OOMPH_EXCEPTION_LOCATION);
      if(k<0) throw OomphLibError("M index must be  >= zero.",
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
#endif
      return 2 + k; // equations 2,3,4
    }

    /// Function determining how to block the Jacobian.
    void get_dof_numbers_for_unknowns(std::list<std::pair<unsigned long,unsigned> >&
                                      block_lookup_list) const
    {
      // Temporarily use a map for storage because it's easier to overwrite
      // entries.
      std::map<unsigned long, unsigned> block_lookup_map;

      // Different dof numbers for boundary values (because BEM means they
      // sometimes have different equations), put them at the end for
      // simplicity.
      int boundary_phi_dof_number = required_nvalue(0);
      int boundary_phi_1_dof_number = required_nvalue(0) + 1;

      // Loop over all nodes then all unpinned values (dofs) at each node. For
      // each of these we create a pair giving the global equation number and
      // the corresponding dof type (number).
      for(unsigned nd=0; nd<nnode(); nd++)
        {
          Node* nd_pt = node_pt(nd);

          // Put it into the block lookup list
          for(unsigned index = 0, nindex=node_pt(nd)->nvalue(); index<nindex; index++)
            {
              int local_eqn_number = this->nodal_local_eqn(nd,index);
              if(local_eqn_number >= 0)
                {
                  int global_eqn_number = eqn_number(local_eqn_number);
                  block_lookup_map[global_eqn_number] = index;
                }
            }

          // If it's a boundary node then move the phi dofs to new blocks
          // of their own (these are BEM values).
          if(nd_pt->is_on_boundary())
            {
              int phi_local_eqn_number = this->nodal_local_eqn(nd, phi_index_micromag());
              if(phi_local_eqn_number >= 0)
                {
                  int global_eqn_number = eqn_number(phi_local_eqn_number);
                  block_lookup_map[global_eqn_number] = boundary_phi_dof_number;
                }

              int phi_1_local_eqn_number = this->nodal_local_eqn(nd, phi_1_index_micromag());
              if(phi_1_local_eqn_number >= 0)
                {
                  int global_eqn_number = eqn_number(phi_1_local_eqn_number);
                  block_lookup_map[global_eqn_number] = boundary_phi_1_dof_number;
                }
            }
        }

      // Convert to a list
      block_lookup_list.assign(block_lookup_map.begin(),
                               block_lookup_map.end());
    }

    /// \short 7 dof types for preconditioning: the 5 values and 2 more for
    /// boundary phi values.
    unsigned ndof_types() const
    {return required_nvalue(0) + 2;}


    // Pointers to things..
    // ============================================================

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

    /// Get the applied field at Eulerian position x.
    virtual Vector<double> get_applied_field(const double& t,
                                             const Vector<double> &x) const
    {
      return magnetic_parameters_pt()->h_app(t, x);
    }

    /// Get the crystalline anisotropy field at Eulerian position x.
    void get_H_cryst_anis_field(const double& t,
                                const Vector<double> &x,
                                const Vector<double>& m,
                                Vector<double> &h_ca) const
    {
      magnetic_parameters_pt()->
        crystalline_ansiotropy_field(t, x, m, h_ca);
    }

    void get_hca_derivative(const double& t, const Vector<double>&x,
                            const Vector<double>& m,
                            const double shape_fn_k_at_x,
                            double dhcadm[3][3]) const
    {
      magnetic_parameters_pt()->
        crystalline_ansiotropy_field_derivative(t,x,m,shape_fn_k_at_x,dhcadm);
    }

    /// Get LLG damping coefficient.
    inline double llg_damping_coeff() const
    {
      return magnetic_parameters_pt()->gilbert_damping();
    }

    /// Get LLG precession coefficient.
    inline double llg_precession_coeff() const
    {
      return 1.0;
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

    // Magnetostatic field
    // ============================================================

    /// Object to calculate magnetostatic field
    MagnetostaticsCalculator* Ms_calc_pt;

    /// \short Helper function for calculation of magnetostatic field.
    void get_magnetostatic_field(const Vector<double> &s,
                                 Vector<double> &h_magnetostatic) const;

    /// \short Calculation of magnetostatic field. Optimised version for
    /// calculations when we aleady have an interpolator (e.g. during
    /// residual calculations).
    virtual void get_magnetostatic_field(CachingMMArrayInterpolator* intp_pt,
                                         Vector<double> &h_magnetostatic) const;

    /// Get the time derivative of the magnetostatic field at a point.
    virtual void get_magnetostatic_field_time_derivative
    (CachingMMInterpolator* intp_pt, Vector<double> &dh_ms_dt) const;


    // OUTPUT FUNCTIONS
    // ============================================================

    /// Output FE representation of soln: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(const unsigned& t, std::ostream &outfile,
                const unsigned &n_plot) const;

    /// Output exact solution at n_plot points
    void output_fct(std::ostream &outfile, const unsigned &n_plot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt);

    void output_fct(std::ostream &outfile, const unsigned &n_plot,
                    const double& time,
                    const SolutionFunctorBase& exact_soln) const;

    /// Get error by comparing with exact solution and get norm of exact solution.
    void compute_error(std::ostream &outfile,
                       FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                       const double& time, double& error, double& norm);

    // RESIDUALS + JACOBIAN
    // ============================================================

    /// Add the element's contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      //Call the generic residuals function with flag set to 0 using a dummy matrix argument
      fill_in_generic_residual_contribution_micromag
        (residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    /// \short Add the element's contribution to its residual vector and element
    /// Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                          DenseMatrix<double> &jacobian)
    {
      fill_in_generic_residual_contribution_micromag
        (residuals, jacobian, 1);
    }

    /// This might go better inside generic get jacobian etc. once I write
    /// it for LL form.
    void fill_in_contribution_to_mass_matrix(Vector<double> &residuals,
                                             DenseMatrix<double> &mmatrix)
    {
      if(Residual_calculator_pt->use_gilbert_form())
        {
          std::string err = "Cannot do explicit time steps for Gilbert form!";
          err += " (well, it's probably possible but much easier to just use LL form)";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }

      // Get the residuals
      fill_in_contribution_to_residuals(residuals);

      const unsigned n_node = this->nnode();
      const unsigned eldim = this->dim();
      const unsigned n_unknowns = required_nvalue(0);

      Shape psi(n_node), test(n_node);
      Vector<double> s(eldim);

      //Loop over the integration points
      for(unsigned ipt=0, nipt=this->integral_pt()->nweight(); ipt<nipt; ipt++)
        {
          // Get position
          for(unsigned j=0; j<eldim; j++)
            {s[j] = this->integral_pt()->knot(ipt,j);}

          // Get shape/test
          shape(s, psi);
          test = psi;

          // Get integration weight and J of transformation
          double W = this->integral_pt()->weight(ipt) * J_eulerian(s);

          // Loop over the dofs in a node
          for(unsigned i=0;i<n_unknowns;i++)
            {
              // Double loop over nodes: equations and unknowns
              for(unsigned l=0;l<n_node;l++)
                {
                  int local_eqn = this->nodal_local_eqn(l,i);
                  if(local_eqn < 0) continue;

                  for(unsigned l2=0;l2<n_node;l2++)
                    {
                      int local_unknown = this->nodal_local_eqn(l2,i);
                      if(local_unknown < 0) continue;

                      mmatrix(local_eqn, local_unknown) +=
                        -psi(l2)*test(l)*W;
                    }
                }
            }
        }
    }


  protected:
    /// Fill in contribution to residuals and jacobian (if flag is set) from
    /// these equations (compatible with multiphysics).
    void fill_in_generic_residual_contribution_micromag(Vector<double> &residuals,
                                                        DenseMatrix<double> &jacobian,
                                                        const unsigned& flag);

  public:


    // Energy and other auxilary calculations
    // ============================================================

    /// \short Integrate a function given by func_pt over the element using
    /// the given integral_pt(). Because C++ sucks we have to do this with
    /// weird function objects. Optionally provide the quadrature to use.
    double integrate_over_element(const ElementalFunction* func_pt,
                                  const Integral* quadrature_pt=0) const;


    /// \short Return a vector containing the magnetisation at a node.
    Vector<double> get_m(const unsigned& t_hist,
                         unsigned node) const
    {
      Vector<double> m(3, 0.0);
      for(unsigned j=0; j<3; j++)
        {
          m[j] = nodal_value(t_hist, node, m_index_micromag(j));
        }
      return m;
    }
    Vector<double> get_m(unsigned node) const
    {
      return get_m(0, node);
    }

    /// \short Get the maximum difference in angle between the
    /// magnetisation of two nodes of the element. If this is large there
    /// is likely an error somewhere (in code or in problem setup), or
    /// there is not enough refinement.
    double max_m_angle_variation(const unsigned& t_hist) const
    {
      double max_angle = 0;

      // Double loop over nodes: compare magnetisation angles at all nodes
      // to each other, store the maximum difference of angles.
      for(unsigned nd=0, n_nd=nnode(); nd<n_nd; nd++)
        {
          Vector<double> m1 = get_m(t_hist, nd);
          for(unsigned nd2=0; nd2<n_nd; nd2++)
            {
              Vector<double> m2 = get_m(t_hist, nd2);
              max_angle = std::max(max_angle,
                                   VectorOps::angle_diff(m1, m2));
            }
        }
      return max_angle;
    }


    // Data storage
    // ============================================================

    // Leave these public for now because access functions are a pain and
    // we would need public write access functions anyway (so no "safer"
    // with access functions).

    /// Should the Jacobian be calcuated by (element level) finite
    /// differencing?
    bool Use_fd_jacobian;

    /// Pointer to the class (a glorified functor really) used to calculate
    /// the residual and Jacobian. Done like this rather than hard coded so
    /// that we can easily switch between LL and LLG formulations, e.g. for
    /// use in predictor step of adaptive midpoint.
    LLGResidualCalculator* Residual_calculator_pt;


  protected:

    /// Pointer to poisson source function for phi. Only for testing
    /// purposes since div(M) is our source function in calculation of the
    /// demagnetising potential.
    TimeSpaceToDoubleFctPt Phi_source_pt;

    /// Pointer to poisson source function for phi1. Only for testing
    /// purposes since div(M) is our source function in calculation of the
    /// demagnetising potential.
    TimeSpaceToDoubleFctPt Phi_1_source_pt;

    /// Pointer to class storing the magnetic parameter data.
    const MagneticParameters* Magnetic_parameters_pt;

    /// Pointer to function giving applied field.
    TimeSpaceToDoubleVectFctPt Applied_field_pt;

    // List of face elements attached to this element
    std::set<FiniteElement*> Face_element_pts;

  }; // End of MicromagEquations class

/// Class to handle ??ds
  class RefineableMicromagEquations : public virtual RefineableElement,
                                      public virtual MicromagEquations

  {
  public:
    /// Constructor
    RefineableMicromagEquations() {}

    /// Virtual destructor
    virtual ~RefineableMicromagEquations() {}

    /// \short Number of continuously interpolated values. Note: We assume
    /// that they are located at the beginning of the value_pt Vector!
    /// (Used for interpolation to son elements, for integrity check
    /// and post-processing -- we can only expect
    /// the continously interpolated values to be continous across
    /// element boundaries).
    unsigned ncont_interpolated_values() const
    {
      const unsigned dummy = 0;
      return required_nvalue(dummy);
    }

    /// \short Get all continously interpolated function values at previous
    /// timestep in this element as a Vector. (t=0: present; t>0: prev. timestep)
    /// Note: Vector sets is own size to ensure that
    /// that this function can be used in black-box fashion
    void get_interpolated_values(const unsigned& t,
                                 const Vector<double>&s,
                                 Vector<double>& values)
    {
      GeneralInterpolator intp(this, s);
      values = intp.value();
    }

    ///  Further build: Copy some pointers from father element
    void further_build()
    {
      RefineableMicromagEquations* fele_pt =
        checked_dynamic_cast<RefineableMicromagEquations*>(father_element_pt());

      this->Use_fd_jacobian = fele_pt->Use_fd_jacobian;
      this->Residual_calculator_pt = fele_pt->Residual_calculator_pt;
      this->Magnetic_parameters_pt = fele_pt->Magnetic_parameters_pt;
      this->Phi_source_pt = fele_pt->phi_source_pt();
      this->Phi_1_source_pt = fele_pt->phi_1_source_pt();
    }


  private:
    /// Broken copy constructor
    RefineableMicromagEquations(const RefineableMicromagEquations& dummy)
    {BrokenCopy::broken_copy("RefineableMicromagEquations");}

    /// Broken assignment operator
    void operator=(const RefineableMicromagEquations& dummy)
    {BrokenCopy::broken_assign("RefineableMicromagEquations");}

  };



  //====================================================================
  /// A class combining the micromag equations with a QElement geometry
  //====================================================================
  template < unsigned DIM, unsigned NNODE_1D>
  class QMicromagElement : public virtual RefineableQElement<DIM>,
                           public virtual QElement<DIM,NNODE_1D>,
                           public virtual RefineableMicromagEquations
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

    // /// Output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    // void output(std::ostream &outfile, const unsigned &n_plot=5)
    // {MicromagEquations::output(outfile,n_plot);}

    /// \short Perform additional hanging node procedures for variables
    /// that are not interpolated by all nodes. We don't have any of these
    /// in LLG so do nothing.
    void further_setup_hanging_nodes() {}

    /// \short Rebuild the element, e.g. set internal values in line with
    /// those of the sons that have now merged.
    void rebuild_from_sons(Mesh* &mesh_pt)
    {
      // ??ds don't think I need anything here... don't really know
      throw OomphLibError("Not implemented (yet?).", OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
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

    // /// Output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    // void output(std::ostream &outfile, const unsigned &n_plot=5)
    // {MicromagEquations::output(outfile,n_plot);}

    /// \short Perform additional hanging node procedures for variables
    /// that are not interpolated by all nodes. We don't have any of these
    /// in LLG so do nothing.
    void further_setup_hanging_nodes() {}

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
    void magnetostatic_field(const Vector<double> &s,
                             const TimeStepper* ts_pt,
                             Vector<double> &hms) const
    {
      GeneralInterpolator intp(this, s, ts_pt);
      hms = intp.dvaluedx(this->u_index_poisson());

      // Make sure the field has 3 dimensions (even if there are only two
      // spatial dimensions).
      hms.resize(3, 0.0);

      // Multiply by -1 and normalise
      for(unsigned j=0; j<3; j++)
        {
          hms[j] *= -1 * Micromag_element_pt->magnetostatic_coeff();
        }
    }

    /// Get the time derivative of the magnetostatic field at local
    /// coordinate point s in the element.
    void magnetostatic_field_time_derivative(const Vector<double> &s,
                                             const TimeStepper* ts_pt,
                                             Vector<double> &dh_ms_dt) const
    {
#ifdef PARANOID
      if(node_pt(0)->time_stepper_pt() == 0)
        {
          std::string error_msg = "No timestepper in poisson elements, need to assign one in mesh construction to evaluate time derivatives.";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      if(dynamic_cast<const MidpointMethod*>(ts_pt) != 0)
        {
          std::string error_msg = "Probably shouldn't be using midpoint method for time derivatives of things (unless you know what you're doing?)";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      GeneralInterpolator intp(this, s, ts_pt);
      dh_ms_dt = intp.d2valuedxdt(this->u_index_poisson());

      // Make sure the field has 3 dimensions (even if there are only two
      // spatial dimensions).
      dh_ms_dt.resize(3, 0.0);

      // Multiply by -1 and normalise
      for(unsigned j=0; j<3; j++)
        {
          dh_ms_dt[j] *= -1 * Micromag_element_pt->magnetostatic_coeff();
        }
    }

    /// For micromagnetics the source function is the divergence of the
    /// magnetisation.
    void get_source_poisson(const unsigned& ipt, const Vector<double>& x,
                            double& source) const;


    // Access functions
    // ============================================================

    /// \short Non-const access function for Micromag_element_pt.
    void set_micromag_element_pt(MicromagEquations* ele_pt)
    {Micromag_element_pt = ele_pt;}

    /// \short Const access function for Micromag_element_pt.
    MicromagEquations* micromag_element_pt() const
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



  /// \short Simpler (but slower) implementation of micromagnetics
  /// interpolator class.
  class MMInterpolator : public GeneralInterpolator
  {
    // Assumption: m_index_micromag(0 - 3) are consecutive.

  private:

    // Extra storage for magnetisation values, so we can have nice vector
    // access to them.
    Vector<double> Dmdt;
    Vector<double> M;
    double Div_m;

    const MicromagEquations* This_element;

  public:

    /// Default constructor
    MMInterpolator(const FiniteElement* const this_element,
                   const Vector<double> &s)
      : GeneralInterpolator(this_element, s),
        Div_m(InterpolatorHelpers::not_yet_calculated_value())
    {}

    /// Set different timestepper constructor
    MMInterpolator(const FiniteElement* const this_element,
                   const Vector<double> &s,
                   const TimeStepper* ts_pt)
      : GeneralInterpolator(this_element, s, ts_pt),
        Div_m(InterpolatorHelpers::not_yet_calculated_value())
    {}


    double phi() {return this->value(This_element->phi_index_micromag());}
    const Vector<double> & dphidx()
    {return this->dvaluedx(This_element->phi_index_micromag());}

    double phi1() {return this->value(This_element->phi_1_index_micromag());}
    const Vector<double> & dphi1dx()
    {return this->dvaluedx(This_element->phi_1_index_micromag());}

    const Vector<double> &m()
    {
      if(InterpolatorHelpers::uninitialised(M))
        {
          M = this->interpolate_values(This_element->m_index_micromag(0),
                                       This_element->m_index_micromag(2) + 1);
        }
      return M;
    }

    const Vector<double>& dmdt()
    {
      if(InterpolatorHelpers::uninitialised(Dmdt))
        {
          Dmdt = this->interpolate_dvaluesdt(This_element->m_index_micromag(0),
                                             This_element->m_index_micromag(2) + 1);
        }
      return Dmdt;
    }

    const Vector<double>& dmdx(const unsigned &i_val)
    {
      return this->dvaluedx(This_element->m_index_micromag(i_val));
    }

    const Vector<double>& d2mdxdt(const unsigned &i_val)
    {
      return this->d2valuedxdt(This_element->m_index_micromag(i_val));
    }

    double div_m()
    {
      if(InterpolatorHelpers::uninitialised(Div_m))
        {
          Div_m = 0;
          for(unsigned j=0; j<Dim; j++) Div_m += dmdx(j)[j];
        }
      return Div_m;
    }
  };

  class MMArrayInterpolator : public GeneralArrayInterpolator<5>
  {
    // Assumption: m_index_micromag(0 - 3) are consecutive.

  private:

    // Extra storage for magnetisation values, so we can have nicer access
    // to them.
    double* Dmdt;
    double* M;
    double Div_m;

    const MicromagEquations* This_element;

  public:

    /// Default constructor
    MMArrayInterpolator(const FiniteElement* const this_element,
                        const unsigned& time_index=0)
      : GeneralArrayInterpolator<5>(this_element, time_index),
        Dmdt(0), M(0), Div_m(InterpolatorHelpers::not_yet_calculated_value()),
        This_element(checked_dynamic_cast<const MicromagEquations*>(this_element))
    {}


    /// Overload build to also initialise our new variables
    void build(const Vector<double>& s)
    {
      // Initialise this class' storage
      Dmdt = 0;
      M = 0;
      Div_m = InterpolatorHelpers::not_yet_calculated_value();

      // Call the base version to initialise the rest
      GeneralArrayInterpolator<5>::build(s);
    }

    double phi() {return this->value()[This_element->phi_index_micromag()];}
    const double* dphidx()
    {return this->dvaluedx(This_element->phi_index_micromag());}

    double phi1() {return this->value()[This_element->phi_1_index_micromag()];}
    const double* dphi1dx()
    {return this->dvaluedx(This_element->phi_1_index_micromag());}

    const double* m()
    {
      if(InterpolatorHelpers::uninitialised(M))
        {
          this->interpolate_values(This_element->m_index_micromag(0),
                                   This_element->m_index_micromag(2) + 1);
          M = this->Values + This_element->m_index_micromag(0);
        }
      return M;
    }

    const double* dmdt()
    {
      if(InterpolatorHelpers::uninitialised(Dmdt))
        {
          this->interpolate_dvaluesdt(This_element->m_index_micromag(0),
                                      This_element->m_index_micromag(2) + 1);

          Dmdt = this->Dvaluesdt + This_element->m_index_micromag(0);
        }
      return Dmdt;
    }

    const double* dmdx(const unsigned &i_val)
    {
      return this->dvaluedx(This_element->m_index_micromag(i_val));
    }

    double div_m()
    {
      if(InterpolatorHelpers::uninitialised(Div_m))
        {
          Div_m = 0;
          for(unsigned j=0; j<this->Dim; j++) Div_m += dmdx(j)[j];
        }
      return Div_m;
    }
  };


  /// A micromagnetics specific caching interpolator. "Knows about"
  /// magnetisation and phi values, so we can have functions such as m()
  /// which gives the interpolated magnetisation.
  class CachingMMArrayInterpolator : public CachingArrayInterpolatorBase
  {

  public:
    /// Default constructor
    CachingMMArrayInterpolator() {}

    /// Destructor
    virtual ~CachingMMArrayInterpolator() {}

    virtual void build(const Vector<double>& s)
    {
      CachingArrayInterpolatorBase::build(s);

      // Clear storage
      M.reset();
      Phi.reset();
      Phi1.reset();
      Dmdt.reset();
      Dphidx.reset();
      Dphi1dx.reset();
      for(unsigned i=0; i<3; i++) {Dmdx[i].reset();}
      Div_m.reset();
    }

    const double* m()
    {
      if(!M.is_initialised())
        {
          for(unsigned j=0; j<3; j++)
            {
              M.set()[j] = Intp_pt->interpolate_value(ele_pt()->m_index_micromag(j));
            }
        }

      return M.get();
    }

    double phi()
    {
      if(!Phi.is_initialised())
        {
          Phi.set() = Intp_pt->interpolate_value(ele_pt()->phi_index_micromag());
        }

      return Phi.get();
    }

    double phi1()
    {
      if(!Phi1.is_initialised())
        {
          Phi1.set() = Intp_pt->interpolate_value(ele_pt()->phi_1_index_micromag());
        }

      return Phi1.get();
    }

    const double* dmdt()
    {
      if(!Dmdt.is_initialised())
        {
          for(unsigned j=0; j<3; j++)
            {
              Dmdt.set()[j] = Intp_pt->interpolate_dvaluedt(ele_pt()->m_index_micromag(j));
            }
        }

      return Dmdt.get();
    }

    const double* dmdx(const unsigned &i_val)
    {
      if(!Dmdx[i_val].is_initialised())
        {
          for(unsigned j=0; j<Intp_pt->Dim; j++)
            {
              Dmdx[i_val].set()[j] = Intp_pt->
                interpolate_dvaluedx(ele_pt()->m_index_micromag(i_val), j);
            }
        }

      return Dmdx[i_val].get();
    }

    const double* dphidx()
    {
      if(!Dphidx.is_initialised())
        {
          for(unsigned j=0; j<Intp_pt->Dim; j++)
            {
              Dphidx.set()[j] = Intp_pt->
                interpolate_dvaluedx(ele_pt()->phi_index_micromag(), j);
            }
        }

      return Dphidx.get();
    }

    const double* dphi1dx()
    {
      if(!Dphi1dx.is_initialised())
        {
          for(unsigned j=0; j<Intp_pt->Dim; j++)
            {
              Dphi1dx.set()[j] = Intp_pt->
                interpolate_dvaluedx(ele_pt()->phi_1_index_micromag(), j);
            }
        }

      return Dphi1dx.get();
    }

    double div_m()
    {
      if(!Div_m.is_initialised())
        {
          double div_m = 0;
          for(unsigned j=0; j<Intp_pt->Dim; j++)
            {
              div_m += this->dmdx(j)[j];
            }
          Div_m.set() = div_m;
        }
      return Div_m.get();
    }

    const MicromagEquations* ele_pt() const
    {
      return dynamic_cast<const MicromagEquations*>(Intp_pt->This_element);
    }

  protected:

    /// Storage arrays
    CachedArray<3> M;
    Cached<double> Phi;
    Cached<double> Phi1;

    CachedArray<3> Dmdt;

    CachedArray<3> Dmdx[3];
    CachedArray<3> Dphidx;
    CachedArray<3> Dphi1dx;

    Cached<double> Div_m;
  };





  class CachingMMInterpolator : public CachingVectorInterpolatorBase
  {

  public:
    /// Default constructor
    CachingMMInterpolator() {}

    /// Destructor
    virtual ~CachingMMInterpolator() {}

    virtual void build(const Vector<double>& s)
    {
      CachingVectorInterpolatorBase::build(s);

      // Clear storage
      M.reset();
      Phi.reset();
      Phi1.reset();
      Dmdt.reset();
      Dphidx.reset();
      Dphi1dx.reset();
      for(unsigned i=0; i<3; i++) {Dmdx[i].reset();}
      Div_m.reset();
    }

    const Vector<double>& m()
    {
      if(!M.is_initialised())
        {
          M.set().assign(3, 0.0);
          for(unsigned j=0; j<3; j++)
            {
              M.set()[j] = Intp_pt->interpolate_value(ele_pt()->m_index_micromag(j));
            }
        }

      return M.get();
    }

    double phi()
    {
      if(!Phi.is_initialised())
        {
          Phi.set() = Intp_pt->interpolate_value(ele_pt()->phi_index_micromag());
        }

      return Phi.get();
    }

    double phi1()
    {
      if(!Phi1.is_initialised())
        {
          Phi1.set() = Intp_pt->interpolate_value(ele_pt()->phi_1_index_micromag());
        }

      return Phi1.get();
    }

    const Vector<double>& dmdt()
    {
      if(!Dmdt.is_initialised())
        {
          Dmdt.set().assign(3, 0.0);
          for(unsigned j=0; j<3; j++)
            {
              Dmdt.set()[j] = Intp_pt->interpolate_dvaluedt(ele_pt()->m_index_micromag(j));
            }
        }

      return Dmdt.get();
    }

    const Vector<double>& dmdx(const unsigned &i_val)
    {
      if(!Dmdx[i_val].is_initialised())
        {
          Dmdx[i_val].set().assign(3, 0.0);
          for(unsigned j=0; j<Intp_pt->Dim; j++)
            {
              Dmdx[i_val].set()[j] = Intp_pt->
                interpolate_dvaluedx(ele_pt()->m_index_micromag(i_val), j);
            }
        }

      return Dmdx[i_val].get();
    }

    const Vector<double>& dphidx()
    {
      if(!Dphidx.is_initialised())
        {
          Dphidx.set().assign(3, 0.0);
          for(unsigned j=0; j<Intp_pt->Dim; j++)
            {
              Dphidx.set()[j] = Intp_pt->
                interpolate_dvaluedx(ele_pt()->phi_index_micromag(), j);
            }
        }

      return Dphidx.get();
    }

    const Vector<double>& dphi1dx()
    {
      if(!Dphi1dx.is_initialised())
        {
          Dphi1dx.set().assign(3, 0.0);
          for(unsigned j=0; j<Intp_pt->Dim; j++)
            {
              Dphi1dx.set()[j] = Intp_pt->
                interpolate_dvaluedx(ele_pt()->phi_1_index_micromag(), j);
            }
        }

      return Dphi1dx.get();
    }

    double div_m()
    {
      if(!Div_m.is_initialised())
        {
          double div_m = 0;
          for(unsigned j=0; j<Intp_pt->Dim; j++)
            {
              div_m += this->dmdx(j)[j];
            }
          Div_m.set() = div_m;
        }
      return Div_m.get();
    }

    const MicromagEquations* ele_pt() const
    {
      return dynamic_cast<const MicromagEquations*>(Intp_pt->This_element);
    }

  protected:

    /// Storage arrays
    Cached<Vector<double> > M;
    Cached<double> Phi;
    Cached<double> Phi1;

    Cached<Vector<double> > Dmdt;

    Cached<Vector<double> > Dmdx[3];
    Cached<Vector<double> > Dphidx;
    Cached<Vector<double> > Dphi1dx;

    Cached<double> Div_m;
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

  template<unsigned NNODE_1D>
  class FaceGeometry<QMagnetostaticFieldElement<1, NNODE_1D> >:
    public PointElement {};

  template<unsigned NNODE_1D>
  class FaceGeometry<TMagnetostaticFieldElement<1, NNODE_1D> >:
    public PointElement {};

  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<TMicromagElement<DIM,NNODE_1D> >:
    public virtual TElement<DIM-1,NNODE_1D> {};

  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<QMicromagElement<DIM,NNODE_1D> >:
    public virtual QElement<DIM-1,NNODE_1D> {};

  template<unsigned NNODE_1D>
  class FaceGeometry<QMicromagElement<1, NNODE_1D> >:
    public PointElement {};

  template<unsigned NNODE_1D>
  class FaceGeometry<TMicromagElement<1,NNODE_1D> >:
    public PointElement {};

}




#endif
