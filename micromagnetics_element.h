#ifndef OOMPH_MICROMAGNETICS_ELEMENTS_HEADER
#define OOMPH_MICROMAGNETICS_ELEMENTS_HEADER


// Generic oomph-lib routines
#include "generic.h"

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
    MicromagEquations() : // Poisson_source_pt(0), Llg_source_pt(0),
			  Applied_field_pt(0), Cryst_anis_field_pt(0),
			  Sat_mag_pt(0), Llg_damp_pt(0),
			  Llg_precess_pt(0), Exchange_coeff_pt(0),
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
    unsigned phi_1_index_micromag() const {return 1;}

    /// Specify nodal index for kth component of M.
    unsigned m_index_micromag(const unsigned &k) const {return 2 + k;}

    /// Specify nodal index for kth component of H_ex.
    unsigned exchange_index_micromag(const unsigned &k) const {return 6 + k;}


    // SOURCE FUNCTIONS
    typedef double (*TimeSpaceToDoubleFctPt)(const double& t, const Vector<double>&x);

    typedef void (*TimeSpaceToDoubleVectFctPt)
    (const double& t, const Vector<double>&x,Vector<double>& out);

    typedef void (*TimeSpaceMagToDoubleVectFctPt)
    (const double& t, const Vector<double>&x, const Vector<double>& M, Vector<double>& out);

    // /// \short Function pointer to source function fct(x,f(x)) --
    // /// x is a Vector!
    // typedef void (*PoissonSourceFctPt)(const double& t, const Vector<double>& x, double& f);

    // /// Access function: Pointer to source function
    // PoissonSourceFctPt& poisson_source_pt() {return Poisson_source_pt;}

    // /// Access function: Pointer to source function. Const version
    // PoissonSourceFctPt poisson_source_pt() const {return Poisson_source_pt;}

    // /// Get poisson source term at (Eulerian) position x.
    // inline void get_poisson_source(const double& t,
    // 					   const unsigned& ipt,
    // 					   const Vector<double>& x,
    // 					   double& source) const
    // {
    //   if(Poisson_source_pt==0) {source = 0.0;}
    //   else (*Poisson_source_pt)(t,x,source);
    // }


    // // LLG SOURCE FUNCTION
    // /// Access function: Pointer to source function
    // TimeSpaceToDoubleVectFctPt& llg_source_pt() {return Llg_source_pt;}

    // /// Access function: Pointer to source function. Const version
    // TimeSpaceToDoubleVectFctPt llg_source_pt() const {return Llg_source_pt;}

    // /// Get LLG source term at (Eulerian) position x.
    // inline void get_source_llg(const double& t,
    // 				       const Vector<double>& x,
    // 				       Vector<double>& source) const
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
    inline void get_applied_field(const double& t, const Vector<double> &x,
				  Vector<double>& H_app) const
    {
      if(Applied_field_pt==0) std::fill(H_app.begin(), H_app.end(), 0.0);
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
      if(Cryst_anis_field_pt==0) std::fill(H_ca.begin(), H_ca.end(), 0.0);
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
    inline double get_llg_damp() const
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
    inline double get_llg_precess() const
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
    inline void get_exact_m(const double& t, const Vector<double>& x, Vector<double>& m_exact) const
    {
      // If no exact M function has been set, return something crazy
      if(Exact_m_pt==0) {for(unsigned j=0;j<3;j++) m_exact[j] = -1000.0;}
      else (*Exact_m_pt)(t,x,m_exact);
    }

    /// Calculate the cross product of vectors A and B, store the result in
    /// vector output. NOTE: the cross product is only valid for 3-dimensional
    /// vectors
    //??ds this should go somewhere else probably, maybe generic.h?
    void cross(Vector<double>& A, Vector<double>& B, Vector<double>& output) const
    {
      output[0] = A[1]*B[2] - A[2]*B[1];
      output[1] = A[2]*B[0] - A[0]*B[2];
      output[2] = A[0]*B[1] - A[1]*B[0];
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
      double interpolated_phi = 0.0;

      //Loop over the local nodes and sum
      for(unsigned l=0;l<n_node;l++)
	{
	  interpolated_phi += this->nodal_value(l,phi_nodal_index)*psi[l];
	}

      return(interpolated_phi);
    }

    /// Return FE representation of M at local coordinate s and current time.
    inline void interpolated_m_micromag(const Vector<double> &s,
					Vector<double>& interpolated_m) const
    {
      //Find number of nodes
      const unsigned n_node = nnode();

      //Local shape function
      Shape psi(n_node);

      //Find values of shape function
      shape(s,psi);

      // Initialise m
      for(unsigned i=0; i<3; i++){interpolated_m[i] = 0.0;}

      // Loop over dimensions of M
      for(unsigned k=0; k<3; k++)
	{
	  //Loop over the local nodes and sum
	  for(unsigned l=0;l<n_node;l++)
	    {
	      interpolated_m[k] += this->nodal_value(l,m_index_micromag(k))*psi[l];
	    }
	}

    }

    /// \short Return FE representation of solution vector (phi,M,H_ex)
    /// at local coordinate s and current time.
    inline void interpolated_solution_micromag(const Vector<double> &s,
					       Vector<double>& interpolated_solution) const
    {
      //Find number of nodes
      const unsigned n_node = nnode();

      //Local shape function
      Shape psi(n_node);

      //Find values of shape function
      shape(s,psi);

      // Initialise solution vector
      //??ds need to get length of solution vector really..., replace 7 with it
      for(unsigned i=0; i<7; i++){interpolated_solution[i] = 0.0;}

      // Loop over the list of solutions
      for(unsigned i=0; i<7; i++)
	{
	  //Loop over the local nodes and sum
	  for(unsigned l=0;l<n_node;l++)
	    {
	      interpolated_solution[i] += this->nodal_value(l,i)*psi[l];
	    }
	}

    }

    // OUTPUT FUNCTIONS
    /// Output with default number of plot points
    void output(std::ostream &outfile)
    {
      const unsigned n_plot=5;
      output(outfile,n_plot);
    }

    /// Output FE representation of soln: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(std::ostream &outfile, const unsigned &n_plot);

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

    // ??ds re-activate this when you're confident that the Jacobian is correct!
    // /// Add the element's contribution to its residual vector and element Jacobian matrix (wrapper)
    // void fill_in_contribution_to_jacobian(Vector<double> &residuals,
    //                                       DenseMatrix<double> &jacobian)
    // {
    //   //Call the generic routine with the flag set to 1
    //   fill_in_generic_residual_contribution_micromag(residuals,jacobian,1);
    // }

    /// Fill in contribution to residuals and jacobian (if flag is set) from these equations (compatible with multiphysics)
    void fill_in_generic_residual_contribution_micromag(Vector<double> &residuals,
							DenseMatrix<double> &jacobian,
							const unsigned& flag) const;

    /// Add dM/dt at local node n (using timestepper and history values) to the vector dmdt.
    void dm_dt_micromag(const unsigned &n, Vector<double>& dmdt) const
    {
      // Get the data's timestepper
      TimeStepper* time_stepper_pt= this->node_pt(n)->time_stepper_pt();

      //Initialise dM/dt to zero
      for (unsigned j=0; j<3; j++) {dmdt[j] = 0.0;}

      //Loop over the timesteps, if there is a non Steady timestepper
      if (!time_stepper_pt->is_steady())
	{
	  // Loop over the directions
	  for (unsigned j=0; j<3; j++)
	    {
	      // Get number of timsteps to use (past & present)
	      const unsigned n_time = time_stepper_pt->ntstorage();

	      // Loop over past and present times and add the contributions to
	      // the time derivative
	      for(unsigned t=0;t<n_time;t++)
		{
		  // ??ds The "1" in weights is the derrivative order (i.e. 1:
		  // dM/dt, 2: d^2M/dt^2) I think...
		  dmdt[j] += time_stepper_pt->weight(1,t)*nodal_value(t,n,m_index_micromag(j));
		}
	    }

	}

    } // end of dM_dt_micromag

  protected:

    /// Pointer to poisson source function - only for testing purposes since div(M) is our source function in calculation of the demagnetising potential.
    // PoissonSourceFctPt Poisson_source_pt;

    /// Pointer to LLG source function (for testing purposes)
    // TimeSpaceToDoubleVectFctPt Llg_source_pt;

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

    /// Pointer to the exact solution for M
    TimeSpaceToDoubleVectFctPt Exact_m_pt;

    /// Pointer to the exact solution for phi
    TimeSpaceToDoubleFctPt Exact_phi_pt;

    /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
    virtual double dshape_and_dtest_eulerian_micromag
    (const Vector<double> &s, Shape &psi, DShape &dpsidx,
     Shape &test, DShape &dtestdx) const=0;

    /// Shape, test functions & derivs. w.r.t. to global coords. at integration point ipt. Return Jacobian.
    virtual double dshape_and_dtest_eulerian_at_knot_micromag
    (const unsigned& ipt, Shape &psi, DShape &dpsidx,
     Shape &test, DShape &dtestdx) const=0;

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
    // We need to store phi(1 value), M(3) and H_exchange(3),
    // phi_1 and phi_2.
    //??ds clean this up - probably don't need phi or H-exchange..
    //??ds generalise somehow?
    inline unsigned required_nvalue(const unsigned &n) const
    {return 9;}

    // OUTPUT FUNCTIONS (just call from MicromagEquations class)
    /// Output function: x,y,u or x,y,z,u
    void output(std::ostream &outfile)
    {MicromagEquations<DIM>::output(outfile);}


    /// Output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    void output(std::ostream &outfile, const unsigned &n_plot)
    {MicromagEquations<DIM>::output(outfile,n_plot);}


    // /// C-style output function: x,y,u or x,y,z,u
    // void output(FILE* file_pt)
    // {MicromagEquations<DIM>::output(file_pt);}


    // /// C-style output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    // void output(FILE* file_pt, const unsigned &n_plot)
    // {MicromagEquations<DIM>::output(file_pt,n_plot);}


  protected:

    /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
    inline double dshape_and_dtest_eulerian_micromag(const Vector<double> &s,
						     Shape &psi,
						     DShape &dpsidx,
						     Shape &test,
						     DShape &dtestdx) const;

    /// Shape, test functions & derivs. w.r.t. to global coords. at integration point ipt. Return Jacobian.
    inline double dshape_and_dtest_eulerian_at_knot_micromag(const unsigned& ipt,
							     Shape &psi,
							     DShape &dpsidx,
							     Shape &test,
							     DShape &dtestdx) const;

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


}




#endif
