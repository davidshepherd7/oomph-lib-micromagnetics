#ifndef OOMPH_MICROMAGNETICS_ELEMENTS_HEADER
#define OOMPH_MICROMAGNETICS_ELEMENTS_HEADER


// Generic oomph-lib routines
#include "generic.h"

// Print vectors nicely (c++0x)
#include <iostream>
#include <vector>
#include "./prettyprint98.hpp"

// My vector helpers
#include "./vector_helpers.h"

using namespace oomph;

namespace oomph
{

  //==============================================================================
  /// A class for the maths used in solving the Landau-Lifshitz-Gilbert equations.
  //==============================================================================
  template<unsigned DIM>
  class MicromagEquations : public virtual FiniteElement
  {
  public:

    // CONSTRUCTORS ETC.
    /// Constructor (initialises the various function pointers to null).
    MicromagEquations() : Phi_source_pt(0), Phi_1_source_pt(0), Llg_source_pt(0),
			  Applied_field_pt(0), Cryst_anis_field_pt(0),
			  Sat_mag_pt(0), Llg_damp_pt(0),
			  Llg_precess_pt(0), Exchange_coeff_pt(0),
			  Magnetostatic_coeff_pt(0),
			  Exact_m_pt(0), Exact_phi_pt(0)
    {}

    /// Broken copy constructor
    MicromagEquations(const MicromagEquations& dummy)
    {BrokenCopy::broken_copy("MicromagEquations");}

    /// Broken assignment operator
    void operator=(const MicromagEquations&)
    {BrokenCopy::broken_assign("MicromagEquations");}

    /// Self-test: Return 0 for OK.
    unsigned self_test(){return 0;} //??ds write a real test sometime

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
			    "MicromagEquations::m_index_micromag",
			    OOMPH_EXCEPTION_LOCATION);
#endif
      return 2 + k; // equations 2,3,4
    }

    /// Function determining how to block the Jacobian.
    void get_dof_numbers_for_unknowns(std::list<std::pair<unsigned long,unsigned> >&
				      block_lookup_list)
    {
      // Loop over all nodes then all unpinned values (dofs) at each node. For
      // each of these we create a pair giving the global equation number and
      // the corresponding dof type (number).

      for(unsigned nd=0; nd<nnode(); nd++)
	{
	  for(unsigned dof=0; dof<node_pt(nd)->nvalue(); dof++)
	    {
	      int local_eqn_number = nodal_local_eqn(nd,dof);
	      if(local_eqn_number >= 0)
		{
		  std::pair<unsigned,unsigned> block_lookup;
		  block_lookup.first = eqn_number(local_eqn_number);
		  block_lookup.second = dof;
		  block_lookup_list.push_front(block_lookup);
		  //??ds why do we use push front?
		  //??ds why do we use lists?
		}
	    }
	}

    }

    unsigned ndof_types()
    {
      return required_nvalue(0);
    }

    typedef double (*TimeSpaceToDoubleFctPt)(const double& t, const Vector<double>&x);

    typedef void (*TimeSpaceToDoubleVectFctPt)
    (const double& t, const Vector<double>&x,Vector<double>& out);

    typedef void (*TimeSpaceMagToDoubleVectFctPt)
    (const double& t, const Vector<double>&x, const Vector<double>& M, Vector<double>& out);


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

    /// Access function: Pointer to source function
    TimeSpaceToDoubleVectFctPt& llg_source_pt() {return Llg_source_pt;}

    /// Access function: Pointer to source function. Const version
    TimeSpaceToDoubleVectFctPt llg_source_pt() const {return Llg_source_pt;}

    /// Get LLG source term at (Eulerian) position x.
    inline void get_source_llg(const double& t,
			       const Vector<double>& x,
			       Vector<double>& source) const
    {
      if(Llg_source_pt==0) {for(unsigned j=0;j<3;j++) source[j] = 0.0;}
      else (*Llg_source_pt)(t,x,source);
    }


    // APPLIED FIELD
    /// Access function: Pointer to applied field function
    TimeSpaceToDoubleVectFctPt& applied_field_pt() {return Applied_field_pt;}

    /// Access function: Pointer to applied field function. Const version
    TimeSpaceToDoubleVectFctPt applied_field_pt() const {return Applied_field_pt;}

    /// Get the applied field at Eulerian position x.
    inline void get_applied_field(const double& t, const Vector<double> &x,
				  Vector<double>& H_app) const
    {
      if(Applied_field_pt==0) H_app.assign(3,0.0);
      else (*Applied_field_pt)(t,x,H_app);
    }

    // EFFECTIVE ANISOTROPY FIELD
    /// Access function: Pointer to crystalline anisotropy field function
    TimeSpaceMagToDoubleVectFctPt& cryst_anis_field_pt() {return Cryst_anis_field_pt;}

    /// Access function: Pointer to crystalline anisotropy field function. Const version
    TimeSpaceMagToDoubleVectFctPt cryst_anis_field_pt() const {return Cryst_anis_field_pt;}

    /// Get the crystalline anisotropy field at Eulerian position x.
    inline void get_H_cryst_anis_field(const double& t,
				       const Vector<double> &x,
				       const Vector<double>& m,
				       Vector<double>& H_ca) const
    {
      if(Cryst_anis_field_pt==0) H_ca.assign(3,0.0);
      else (*Cryst_anis_field_pt)(t,x,m,H_ca);
    }

    /// Access function: Pointer to saturisation magnetisation
    double*& sat_mag_pt() {return Sat_mag_pt;}

    /// Access function: Pointer to saturisation magnetisation, const version
    double* sat_mag_pt() const {return Sat_mag_pt;}

    /// Get saturisation magnetisation at eulerian postition x.
    inline double get_sat_mag() const
    {
      if(Sat_mag_pt==0) return 1.0;
      else return *Sat_mag_pt;
    }

    /// Access function: Pointer to LLG damping coefficient
    double*& llg_damp_pt() {return Llg_damp_pt;}

    /// Access function: Pointer to LLG damping coefficient, const version
    double* llg_damp_pt() const {return Llg_damp_pt;}

    /// Get LLG damping coefficient at eulerian postition x.
    inline double get_llg_damping_coeff() const
    {
      if(Llg_damp_pt==0) {return 1.0;}
      else return *Llg_damp_pt;
    }

    // LLG PRECESSION COEFF FUNCTION POINTER
    /// Access function: Pointer to LLG precession coefficient function
    double*& llg_precess_pt() {return Llg_precess_pt;}

    /// Access function: Pointer to LLG precession coefficient function. Const version
    double* llg_precess_pt() const {return Llg_precess_pt;}

    /// Get LLG precession coefficient at eulerian postition x.
    inline double get_llg_precession_coeff() const
    {
      if(Llg_precess_pt==0) {return 1.0;}
      else return *Llg_precess_pt;
    }

    // EXCHANGE COEFF FUNCTION POINTER
    /// Access function: Pointer to exchange coefficient function
    double*& exchange_coeff_pt() {return Exchange_coeff_pt;}

    /// Access function: Pointer to exchange coefficient function. Const version
    double* exchange_coeff_pt() const {return Exchange_coeff_pt;}

    /// Get exchange coefficient at eulerian postition x.
    inline double get_exchange_coeff(const double& t, const Vector<double>& x) const
    {
      if(Exchange_coeff_pt==0) {return 1.0;}
      else return *Exchange_coeff_pt;
    }

    // MAGNETOSTATIC COEFF FUNCTION POINTER
    /// Access function: Pointer to magnetostatic coefficient function
    double*& magnetostatic_coeff_pt() {return Magnetostatic_coeff_pt;}

    /// Access function: Pointer to magnetostatic coefficient function. Const version
    double* magnetostatic_coeff_pt() const {return Magnetostatic_coeff_pt;}

    /// Get magnetostatic coefficient at eulerian postition x.
    inline double get_magnetostatic_coeff(const double& t, const Vector<double>& x) const
    {
      if(Magnetostatic_coeff_pt==0) {return 1.0;}
      else return *Magnetostatic_coeff_pt;
    }

    // EXACT PHI FUNCTION POINTER
    /// Access function: Pointer to exact phi function
    TimeSpaceToDoubleFctPt& exact_phi_pt() {return Exact_phi_pt;}

    /// Access function: Pointer to exact phi function. Const version
    TimeSpaceToDoubleFctPt exact_phi_pt() const {return Exact_phi_pt;}

    /// Get exact phi at eulerian postition x.
    inline double get_exact_phi(const double& t, const Vector<double>& x) const
    {
      // If no exact phi function has been set, return something crazy
      if(Exact_phi_pt==0) {return -1000.0;}
      else return (*Exact_phi_pt)(t,x);
    }

    // EXACT M FUNCTION POINTER
    /// Access function: Pointer to exact M function
    TimeSpaceToDoubleVectFctPt& exact_m_pt() {return Exact_m_pt;}

    /// Access function: Pointer to LLG exact M. Const version
    TimeSpaceToDoubleVectFctPt exact_m_pt() const {return Exact_m_pt;}

    /// Get exact M at eulerian postition x.
    inline void get_exact_m(const double& t, const Vector<double>& x,
			    Vector<double>& m_exact) const
    {
      // If no exact M function has been set empty the vector
      if(Exact_m_pt==0) m_exact.clear();
      else (*Exact_m_pt)(t,x,m_exact);
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
      //Find number of nodes
      const unsigned n_node = nnode();

      //Local shape function
      Shape psi(n_node);

      //Find values of shape function
      shape(s,psi);

      // Initialise m
      for(unsigned i=0; i<3; i++){itp_m[i] = 0.0;}

      // Loop over dimensions of M
      for(unsigned k=0; k<3; k++)
	{
	  //Loop over the local nodes and sum
	  for(unsigned l=0;l<n_node;l++)
	    {
	      itp_m[k] += this->nodal_value(l,m_index_micromag(k))*psi[l];
	    }
	}

    }

    /// Return FE representation of M at local coordinate s and current time.
    inline void interpolated_dmdx_micromag(const Vector<double> &s,
					   DenseDoubleMatrix& itp_dmdx) const
    {
      // Get number of nodes
      const unsigned n_node = nnode();

      // Set up memory for the shape and test functions
      Shape psi(n_node), test(n_node);
      DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);

      // Get shape and test functions
      dshape_dtest(s,psi,dpsidx,test,dtestdx);

      // Interpolate values at knot by looping over nodes adding contributions
      for(unsigned l=0;l<n_node;l++)
	for(unsigned j=0; j<DIM; j++)
	  for(unsigned k=0; k<3; k++)
	    {
	      itp_dmdx(k,j) += nodal_value(l,m_index_micromag(k))*dpsidx(l,j);
	    }

    }

    inline void interpolated_dphidx_micromag(const Vector<double>& s,
					     Vector<double>& itp_dphidx) const
    {
      //Find number of nodes
      const unsigned n_node = nnode();

      //Local shape function and derivative
      Shape psi(n_node); DShape dpsidx(n_node,DIM);
      dshape_eulerian(s,psi,dpsidx);

      // Initialise output vector
      itp_dphidx.assign(3,0.0);

      // Calculate values
      for(unsigned l=0;l<n_node;l++)
	for(unsigned j=0; j<DIM; j++)
	  itp_dphidx[j] += nodal_value(l,phi_index_micromag())*dpsidx(l,j);
    }


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


    /// Get total field at position x. There is a more efficient private version
    /// of this function for use in residual calculations etc. This version does
    /// some extra calculations to get interpolated m and x and ensures
    /// initilisations.
    //??ds unchecked
    void interpolated_ht_micromag(const Vector<double>& s,
				  Vector<double>& h_total)
    {
      Vector<double> x(3,0.0);
      interpolated_x(s,x);

      Vector<double> m(3,0.0);
      interpolated_m_micromag(s,m);

      Vector<double> itp_dphidx(3,0.0);
      interpolated_dphidx_micromag(s,itp_dphidx);

      h_total.assign(3,0.0);
      interpolated_ht_micromag_efficient(x,m,itp_dphidx,h_total);
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

    // ??ds testing if the Jacobian is correct!
    /// \short Add the element's contribution to its residual vector and element
    /// Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double> &residuals,
    					  DenseMatrix<double> &jacobian)
    {
      // Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_micromag(residuals,jacobian,1);
    }

    /// Add dM/dt at local node n (using timestepper and history values) to the vector dmdt.
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
		  // ??ds The "1" in weights is the derrivative order (i.e. 1:
		  // dM/dt, 2: d^2M/dt^2) I think...
		  dmdt[j] += time_stepper_pt->weight(1,t)
		    *nodal_value(t,l,m_index_micromag(j));
		}
	    } // End of loop over directions
	}

    } // end of dM_dt_micromag

  protected:

    /// Fill in contribution to residuals and jacobian (if flag is set) from
    /// these equations (compatible with multiphysics)
    void fill_in_generic_residual_contribution_micromag(Vector<double> &residuals,
							DenseMatrix<double> &jacobian,
							const unsigned& flag) const;

    /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
    virtual double dshape_dtest(const Vector<double> &s, Shape &psi, DShape &dpsidx,
				Shape &test, DShape &dtestdx) const=0;

    /// Get the field at current time (pass in x,m, dphidx for
    /// efficiency).
    inline void interpolated_ht_micromag_efficient(const Vector<double>& x,
						   const Vector<double>& m,
						   const Vector<double>& itp_dphidx,
						   Vector<double>& h_total)
      const
    {
      // Get time
      const double time = time_pt()->time();

      Vector<double> h_applied(3,0.0);
      get_applied_field(time, x, h_applied);
      Vector<double> h_cryst_anis(3,0.0);
      get_H_cryst_anis_field(time, x, m, h_cryst_anis);

      // Magnetostatic field is -1* dphi/dx, multiply by a coeff too alow easy
      // switching on and off for debugging.
      double magnetostatic_coeff = get_magnetostatic_coeff(time,x);
      Vector<double> h_ms(3,0.0);
      for(unsigned i=0; i<h_ms.size(); i++)
	h_ms[i] *= -1 * magnetostatic_coeff;

      // Take total of all fields used
      for(unsigned j=0; j<3; j++)
	{
	  h_total[j] = h_applied[j] + h_ms[j] + h_cryst_anis[j];
	}
    }

    /// Pointer to poisson source function - only for testing purposes since
    /// div(M) is our source function in calculation of the demagnetising
    /// potential.
    TimeSpaceToDoubleFctPt Phi_source_pt;

    TimeSpaceToDoubleFctPt Phi_1_source_pt;

    /// Pointer to LLG source function (for testing purposes)
    TimeSpaceToDoubleVectFctPt Llg_source_pt;

    /// Pointer to function giving applied field.
    TimeSpaceToDoubleVectFctPt Applied_field_pt;

    /// Pointer to function giving effective field due to the crystalline anisotropy
    TimeSpaceMagToDoubleVectFctPt Cryst_anis_field_pt;

    /// Pointer to saturisation magnetisation
    double* Sat_mag_pt;

    /// Pointer to LLG damping coefficient
    double* Llg_damp_pt;

    /// Pointer to LLG precession coefficient
    double* Llg_precess_pt;

    /// Pointer to exchange coefficient
    double* Exchange_coeff_pt;

    /// Pointer to magnetostatic coefficient
    double* Magnetostatic_coeff_pt;

    /// Pointer to the exact solution for M
    TimeSpaceToDoubleVectFctPt Exact_m_pt;

    /// Pointer to the exact solution for phi
    TimeSpaceToDoubleFctPt Exact_phi_pt;

  }; // End of MicromagEquations class

  //====================================================================
  /// A class combining the micromag equations with a QElement geometry
  //====================================================================
  template < unsigned DIM, unsigned NNODE_1D>
  class QMicromagElement : public QElement<DIM,NNODE_1D>, public MicromagEquations<DIM>
  {
  private:
    /// Static int that holds the number of variables at nodes: always the same
    static const unsigned Initial_Nvalue;

  public:
    /// Constructor: Call constructors for QElement and Micromag equations
    QMicromagElement() : QElement<DIM,NNODE_1D>(), MicromagEquations<DIM>()
    {}

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

    /// Required  # of `values' (pinned or dofs) at node n.
    inline unsigned required_nvalue(const unsigned &n) const
    {return Initial_Nvalue;}

    /// Output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(std::ostream &outfile, const unsigned &n_plot=5)
    {MicromagEquations<DIM>::output(outfile,n_plot);}

    // /// C-style output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    // void output(FILE* file_pt, const unsigned &n_plot = 5)
    // {MicromagEquations<DIM>::output(file_pt,n_plot);}

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

  }; // end of QMicromagElement class declaration

  //=======================================================================
  /// Face geometry for the QMicromagElement elements: The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<QMicromagElement<DIM,NNODE_1D> >:
    public virtual QElement<DIM-1,NNODE_1D>
  {

  public:

    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : QElement<DIM-1,NNODE_1D>() {}

  };

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //====================================================================
  /// A class combining the micromag equations with a TElement geometry
  //====================================================================
  template < unsigned DIM, unsigned NNODE_1D>
  class TMicromagElement : public TElement<DIM,NNODE_1D>, public MicromagEquations<DIM>
  {
  private:
    /// Static int that holds the number of variables at nodes: always the same
    static const unsigned Initial_Nvalue;

  public:
    /// Constructor: Call constructors for TElement and Micromag equations
    TMicromagElement() : TElement<DIM,NNODE_1D>(), MicromagEquations<DIM>()
    {}

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

    /// Required  # of `values' (pinned or dofs) at node n.
    inline unsigned required_nvalue(const unsigned &n) const
    {return Initial_Nvalue;}

    /// Output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(std::ostream &outfile, const unsigned &n_plot=5)
    {MicromagEquations<DIM>::output(outfile,n_plot);}

    // /// C-style output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    // void output(FILE* file_pt, const unsigned &n_plot = 5)
    // {MicromagEquations<DIM>::output(file_pt,n_plot);}

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

  }; // end of TMicromagElement class declaration

  //=======================================================================
  /// Face geometry for the TMicromagElement elements: The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<TMicromagElement<DIM,NNODE_1D> >:
    public virtual TElement<DIM-1,NNODE_1D>
  {

  public:

    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
    FaceGeometry() : TElement<DIM-1,NNODE_1D>() {}

  };



}




#endif
