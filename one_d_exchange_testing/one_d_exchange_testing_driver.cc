
# include "../micromagnetics_element.cc"
# include "meshes/one_d_mesh.h"
//# include "./../one_d_demag_testing/parameters-demag4.cc"
# include "./parameters-exchange7.cc"

//============================================================
// Core parameters (others are in parameters files)
//============================================================
namespace OneDMicromagSetup
{
  using namespace MathematicalConstants;
  using namespace OneDMicromagSetup;

  // Prototypes
  // ==========================================================
  void exact_solution(const double& t, const Vector<double>& x,
		      Vector<double>& exact);
  double sat_mag(const double& t, const Vector<double>& x);


  // Constants
  //===========================================================
  double alpha = 0.7;   // Gibert damping constant
  double gamma = 0.3;   // Electromagnetic ratio

  // The coefficient of the precession term of the Landau-Lifschitz-Gilbert equation
  // (M x H)
  double llg_precession_coeff(const double& t, const Vector<double>& x)
  {return gamma/(1 + alpha*alpha);}

  // The coefficient of the damping term of the Landau-Lifschitz-Gilbert equation
  // (M x (M x H))
  double llg_damping_coeff(const double& t, const Vector<double>& x)
  { //??ds ifdef error checking then check ms is nonzero?
    //??ds need to put Ms in here instead of 1 but in these tests Ms can be zero...
    return gamma/(1 + alpha*alpha) * (alpha/1.0);}

  double exchange_coeff(const double& t, const Vector<double>& x)
  {return 1.0;}

  // Crystalline anisotropy field - set to zero if unused
  void cryst_anis_field(const double& t, const Vector<double>& x, const Vector<double>& m, Vector<double>& H_cryst_anis)
  {
    H_cryst_anis[0] = 0.0;
    H_cryst_anis[1] = 0.0;
    H_cryst_anis[2] = 0.0;
  }

  // Applied field - set to zero if unused
  void applied_field(const double& t, const Vector<double>& x, Vector<double>& H_applied)
  {
    H_applied[0] = 0.0;
    H_applied[1] = 0.0;
    H_applied[2] = 0.0;
  }

  // Poisson source, set to zero unless needed for testing
  void poisson_source_function(const double& t, const Vector<double>& x, double& source)
  {
    source = 0.0;
  }

  // Get the saturisation magnetisation at position x
  double sat_mag(const double& t, const Vector<double>& x)
  {
    Vector<double> exact(7,0.0);
    exact_solution(t,x,exact);
    return exact[1]*exact[1] + exact[2]*exact[2] + exact[3]*exact[3];
  }

  double boundary_phi(const double& t, const Vector<double>& x)
  {
    Vector<double> exact(7,0.0);
    exact_solution(t,x,exact);
    return exact[0];
  }

  void initial_m(const double& t, const Vector<double>& x, Vector<double>& m)
  {
    Vector<double> exact(7,0.0);
    exact_solution(t,x,exact);
    for(unsigned i=0; i<3; i++){m[i] = exact[i+1];}
  }

  void boundary_m(const double& t, const Vector<double>& x, Vector<double>& m)
  {
    Vector<double> exact(7,0.0);
    exact_solution(t,x,exact);
    for(unsigned i=0; i<3; i++){m[i] = exact[i+1];}
  }

  double exact_phi_solution(const double& t, const Vector<double>& x)
  {
    Vector<double> exact(7,0.0);
    exact_solution(t,x,exact);
    return exact[0];
  }

  void exact_m_solution(const double& t, const Vector<double>& x,
			Vector<double>& m)
  {
    Vector<double> exact(7,0.0);
    exact_solution(t,x,exact);
    for(unsigned i=0; i<3; i++){m[i] = exact[i+1];}
  }

}; // End of namespace

//==start_of_problem_class============================================
/// 1D Micromag problem in unit interval.
//====================================================================
template<class ELEMENT>
class OneDMicromagProblem : public Problem
{

private:

  /// Pointer to Poisson source function
  MicromagEquations<1>::PoissonSourceFctPt Poisson_source_fct_pt;

  /// Pointer to LLg source function
  MicromagEquations<1>::LlgSourceFctPt Llg_source_fct_pt;

  /// Pointer to applied field function
  MicromagEquations<1>::AppliedFieldFctPt Applied_field_fct_pt;

  /// Pointer to crystalline anisotropy effective field function
  MicromagEquations<1>::CrystAnisFieldFctPt Cryst_anis_field_fct_pt;

  /// Pointer to saturisation magnetisation function
  MicromagEquations<1>::SatMagFctPt Sat_mag_fct_pt;

  /// Pointer to LLG damping coefficient function
  MicromagEquations<1>::LlgDampFctPt Llg_damp_fct_pt;

  /// Pointer to LLG precession coefficient function
  MicromagEquations<1>::LlgPrecessFctPt Llg_precess_fct_pt;

  /// Pointer to exchange coefficient function
  MicromagEquations<1>::ExchangeCoeffFctPt Exchange_coeff_fct_pt;

  /// Pointer to exact M solution function
  MicromagEquations<1>::ExactMFctPt Exact_m_fct_pt;

  /// Pointer to exact phi solution function
  MicromagEquations<1>::ExactPhiFctPt Exact_phi_fct_pt;

  /// Pointer to control node at which the solution is documented ??ds - not sure what this is
  Node* Control_node_pt;

  // Doc info object
  DocInfo Doc_info;

  // Trace file
  std::ofstream Trace_file;

public:

  /// Constructor: Pass number of elements and pointer to source function
  OneDMicromagProblem(const unsigned& n_element,
		      MicromagEquations<1>::PoissonSourceFctPt source_fct_pt,
		      MicromagEquations<1>::LlgSourceFctPt llg_source_fct_pt,
		      MicromagEquations<1>::AppliedFieldFctPt applied_field_fct_pt,
		      MicromagEquations<1>::CrystAnisFieldFctPt cryst_anis_field_fct_pt,
		      MicromagEquations<1>::SatMagFctPt sat_mag_fct_pt,
		      MicromagEquations<1>::LlgDampFctPt llg_damp_fct_pt,
		      MicromagEquations<1>::LlgPrecessFctPt llg_precess_fct_pt,
		      MicromagEquations<1>::ExchangeCoeffFctPt exchange_coeff_fct_pt,
		      MicromagEquations<1>::ExactMFctPt exact_m_fct_pt,
		      MicromagEquations<1>::ExactPhiFctPt exact_phi_fct_pt);

  /// Destructor (empty -- all the cleanup is done in the base class)
  ~OneDMicromagProblem(){};

  /// Update the problem specs before solve: Set boundary conditions
  void actions_before_newton_solve(){};

  /// Update the problem specs after solve
  void actions_after_newton_solve(){};

  /// Doc the solution
  void doc_solution(DocInfo& doc_info, std::ofstream& trace_file);

  // TIME STEPPING FUNCTIONS
  /// Update the problem specs after solve (empty)
  void actions_after_implicit_timestep() {}

  /// Update the problem specs before next timestep
  void actions_before_implicit_timestep();

  /// Set initial condition (incl previous timesteps) according to specified function.
  void set_initial_condition();

}; // end of problem class

/// Set number of values stored at each node (7: phi, 3 M's, 3 H_ex's)
template<unsigned DIM, unsigned NNODE_1D>
const unsigned QMicromagElement<DIM,NNODE_1D>::Initial_Nvalue = 7;



//=====start_of_constructor===============================================
/// ??ds
//========================================================================
template<class ELEMENT>
OneDMicromagProblem<ELEMENT>::
OneDMicromagProblem(const unsigned& n_element,
		    MicromagEquations<1>::PoissonSourceFctPt poisson_source_fct_pt,
		    MicromagEquations<1>::LlgSourceFctPt llg_source_fct_pt,
		    MicromagEquations<1>::AppliedFieldFctPt applied_field_fct_pt,
		    MicromagEquations<1>::CrystAnisFieldFctPt cryst_anis_field_fct_pt,
		    MicromagEquations<1>::SatMagFctPt sat_mag_fct_pt,
		    MicromagEquations<1>::LlgDampFctPt llg_damp_fct_pt,
		    MicromagEquations<1>::LlgPrecessFctPt llg_precess_fct_pt,
		    MicromagEquations<1>::ExchangeCoeffFctPt exchange_coeff_fct_pt,
		    MicromagEquations<1>::ExactMFctPt exact_m_fct_pt,
		    MicromagEquations<1>::ExactPhiFctPt exact_phi_fct_pt) :

  Poisson_source_fct_pt(poisson_source_fct_pt),
  Llg_source_fct_pt(llg_source_fct_pt),
  Applied_field_fct_pt(applied_field_fct_pt),
  Cryst_anis_field_fct_pt(cryst_anis_field_fct_pt),
  Sat_mag_fct_pt(sat_mag_fct_pt),
  Llg_damp_fct_pt(llg_damp_fct_pt),
  Llg_precess_fct_pt(llg_precess_fct_pt),
  Exchange_coeff_fct_pt(exchange_coeff_fct_pt),
  Exact_m_fct_pt(exact_m_fct_pt),
  Exact_phi_fct_pt(exact_phi_fct_pt)
{
  // Allocate the timestepper -- this constructs the Problem's time object with a sufficient amount of storage to store the previous timsteps.
  add_time_stepper_pt(new BDF<2>);

  // Set domain length
  double L=1.0;

  // Build mesh and store pointer in Problem
  Problem::mesh_pt() = new OneDMesh<ELEMENT>(n_element,L,time_stepper_pt());

  // Choose a control node at which the solution is documented
  unsigned control_el = unsigned(n_element/2); // Pick a control element in the middle
  Control_node_pt=mesh_pt()->finite_element_pt(control_el)->node_pt(0);  // Choose its first node as the control node
  std::cout << "Recording trace of the solution at: " << Control_node_pt->x(0) << std::endl;


  // Set up the boundary conditions for this problem:
  // pin the phi values of the nodes at either end - int. by parts of phi
  mesh_pt()->boundary_node_pt(0,0)->pin(0);
  mesh_pt()->boundary_node_pt(1,0)->pin(0);

  // Loop over elements to set pointers to everything
  for(unsigned i=0;i<n_element;i++)
    {
      // Upcast from GeneralisedElement to the present element
      ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

      //Set the source function pointer
      elem_pt->poisson_source_fct_pt() = Poisson_source_fct_pt;

      // Set the LLG source function pointer
      elem_pt->llg_source_fct_pt() = Llg_source_fct_pt;

      // Set the applied field function pointer
      elem_pt->applied_field_fct_pt() = Applied_field_fct_pt;

      // Set the crystalline anisotropy effective field pointer
      elem_pt->cryst_anis_field_fct_pt() = Cryst_anis_field_fct_pt;

      // Set the saturisation magnetisation function pointer
      elem_pt->sat_mag_fct_pt() = Sat_mag_fct_pt;

      // Set the LLG damping coefficient function pointer
      elem_pt->llg_damp_fct_pt() = Llg_damp_fct_pt;

      // Set the LLg precession coefficient function pointer
      elem_pt->llg_precess_fct_pt() = Llg_precess_fct_pt;

      // Set the exact M solution function pointer
      elem_pt->exact_m_fct_pt() = Exact_m_fct_pt;

      // Set the exact phi function pointer
      elem_pt->exact_phi_fct_pt() = Exact_phi_fct_pt;

      // Set pointer to continous time
      elem_pt->time_pt() = time_pt();
    }


  // Setup equation numbering scheme
  assign_eqn_numbers();

} // end of constructor



//=========start of actions_before_implicit_timestep===============================
/// \short Actions before timestep: update the domain, then reset the
/// boundary conditions for the current time.
//========================================================================
template<class ELEMENT>
void OneDMicromagProblem<ELEMENT>::actions_before_implicit_timestep()
{
  // get current time
  double t = time_pt()->time();

  // Get pointer to (0th) element - exact element doesn't matter for this use.
  ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0));

  // Get index at which phi is stored (in nodal data)
  unsigned phi_nodal_index = elem_pt->phi_index_micromag();

  // Set boundary conditions (more general method that works in higher dimensions)
  // Loop over all boundaries
  unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned ibound=0;ibound<num_bound;ibound++)
    {
      // Loop over the nodes on this boundary
      unsigned num_nod=mesh_pt()->nboundary_node(ibound);
      for (unsigned inod=0;inod<num_nod;inod++)
	{
	  // Get x coordinate at this node.
	  Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
	  Vector<double> x(1,nod_pt->x(0));

	  // Get and set conditions on phi.
	  double phi_boundary_value = OneDMicromagSetup::boundary_phi(t,x);
	  nod_pt->set_value(phi_nodal_index,phi_boundary_value);
	}
    }

}


//======================start_of_set_initial_condition====================
/// \short Set initial condition: Assign previous and current values
/// from exact solution.
//========================================================================
template<class ELEMENT>
void OneDMicromagProblem<ELEMENT>::set_initial_condition()
{
  // Backup time in global Time object
  double backed_up_time=time_pt()->time();

  // Past history needs to be established for t=time0-deltat, ...
  // Then provide current values (at t=time0) which will also form
  // the initial guess for the first solve at t=time0+deltat

  //Find number of nodes in mesh
  unsigned num_nod = mesh_pt()->nnode();

  // Set continuous times at previous timesteps:
  // How many previous timesteps does the timestepper use?
  int nprev_steps=time_stepper_pt()->nprev_values();
  Vector<double> prev_time(nprev_steps+1);
  for (int t=nprev_steps;t>=0;t--)
    {
      prev_time[t]=time_pt()->time(unsigned(t));
    }

  // Loop over current & previous timesteps
  for (int t=nprev_steps;t>=0;t--)
    {
      // Continuous time
      double time = prev_time[t];
      std::cout << "setting IC at time =" << time << std::endl;

      // Loop over the nodes to set initial values everywhere
      for (unsigned n=0;n<num_nod;n++)
	{
	  // Get nodal coordinate
	  Vector<double> x(1,0.0);
	  x[0]=mesh_pt()->node_pt(n)->x(0);

	  // Get initial value of solution
	  Vector<double> initial_solution(7,0.0);
	  OneDMicromagSetup::exact_solution(time,x,initial_solution);

	  // Set all components of the solution at this node, timestep t
	  for(unsigned i=0; i<7; i++)
	    {
	      mesh_pt()->node_pt(n)->
		set_value(t,i,initial_solution[i]);
	    }
	}
    }

  // Reset backed up time for global timestepper
  time_pt()->time()=backed_up_time;

} // end of set_initial_condition


//===start_of_doc=========================================================
/// Doc the solution in tecplot format. Label files with label.
//========================================================================
template<class ELEMENT>
void OneDMicromagProblem<ELEMENT>::doc_solution(DocInfo& doc_info, std::ofstream& trace_file)
{

  char filename[100];

  // Number of plot points
  unsigned npts;
  npts=2;

  // Get the current time
  double time = time_pt()->time();

  std::cout << std::endl;
  std::cout << "=================================================" << std::endl;
  std::cout << "Doc'ing solution for t=" << time << std::endl;
  std::cout << "=================================================" << std::endl;


  // Output computed solution
  //-----------------
  std::ofstream soln_file;

  sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
	  doc_info.number());

  soln_file.open(filename);
  mesh_pt()->output(soln_file,npts);
  soln_file.close();

  if(1) //??ds some condition determining if I want full error checking output, maybe ifdef..
    {
      // Output exact solution (at many more points than the computed solution)
      //-----------------------------------------
      std::ofstream exact_file;
      sprintf(filename,"%s/exact%i.dat",doc_info.directory().c_str(),
	      doc_info.number());

      exact_file.open(filename);
      mesh_pt()->output_fct(exact_file, 10*npts, time,
			    OneDMicromagSetup::exact_solution);
      exact_file.close();


      // Output errors
      //------------------------------
      std::ofstream error_file;
      double error_norm(0.0), exact_norm(0.0);

      sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
	      doc_info.number());

      // Do the outputing
      error_file.open(filename);
      mesh_pt()->compute_error(error_file, OneDMicromagSetup::exact_solution,
			       time, error_norm, exact_norm);
      error_file.close();

      // Doc error norm:
      std::cout << "\nNorm of error   : " << sqrt(error_norm) << std::endl;
      std::cout << "Norm of solution : " << sqrt(exact_norm) << std::endl << std::endl;
      std::cout << std::endl;
    }

} // end of doc


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////



  //======start_of_main==================================================
  /// Driver for 1D Micromag problem
  //=====================================================================
int main(int argc, char *argv[])
{

  // Get t_max, dt and n_element intelligently
  double t_max, dt;
  unsigned n_element;
  if(argc == 4)
    {
      // If arguments are provided to function use them
      t_max = atof(argv[1]);
      dt = atof(argv[2]);
      n_element = atoi(argv[3]);
    }
  else
    {
      // Else use the arguments from parameters files
      std::cout << "Not enough arguments, using parameters from file" << std::endl;
      t_max = OneDMicromagSetup::t_max;
      dt = OneDMicromagSetup::dt;
      n_element = OneDMicromagSetup::n_x_elements;
    }

    OneDMicromagProblem<QMicromagElement<1,4> >
    problem(n_element,
	    OneDMicromagSetup::poisson_source_function,
	    OneDMicromagSetup::llg_source_function,
	    OneDMicromagSetup::applied_field,
	    OneDMicromagSetup::cryst_anis_field,
	    OneDMicromagSetup::sat_mag,
	    OneDMicromagSetup::llg_damping_coeff,
	    OneDMicromagSetup::llg_precession_coeff,
	    OneDMicromagSetup::exchange_coeff,
	    OneDMicromagSetup::exact_m_solution,
	    OneDMicromagSetup::exact_phi_solution);

  // SET UP OUTPUT
  // Setup labels for output
  DocInfo doc_info;

  // Output directory
  doc_info.set_directory("RESLT");

  // Output number
  doc_info.number()=0;

  // Open a trace file
  std::ofstream trace_file;
  char filename[100];
  sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
  trace_file.open(filename);
  trace_file << "VARIABLES=\"time\",\"u<SUB>FE</SUB>\","
	     << "\"u<SUB>exact</SUB>\",\"norm of error\",\"norm of solution\""
	     << std::endl;

  // Initialise timestep -- also sets the weights for all timesteppers
  // in the problem.
  problem.initialise_dt(dt);

  // Set initial condition (on M)
  problem.set_initial_condition();

  // Output initial conditions
  problem.doc_solution(doc_info,trace_file);   //Output initial condition
  doc_info.number()++; //Increment counter for solutions

  // Check whether the problem can be solved
  std::cout << "\n\n\nProblem self-test ";
  if (problem.self_test()==0)
    {
      std::cout << "passed: Problem can be solved." << std::endl;
    }
  else
    {
      throw OomphLibError("failed!",
			  "main()",
			  OOMPH_EXCEPTION_LOCATION);
    }

  // //  ??ds testing stuff
  // DenseDoubleMatrix jacobian;
  // DoubleVector res;
  // problem.get_jacobian(res,jacobian);
  // jacobian.DenseMatrix<double>::sparse_indexed_output(std::cout);
  // // ??ds hacky way to output last entry of matrix for matlab to pickup on
  // std::cout << jacobian.nrow() -1 << " " << jacobian.ncol() -1 << " 0" << std:: endl;

  // SOLVE THE PROBLEM
  // Find number of steps
  unsigned nstep = unsigned(t_max/dt);

  // Timestepping loop
  for(unsigned istep=0; istep<nstep; istep++)
    {
      std::cout << "Timestep " << istep << std::endl;

      // Take timestep
      problem.unsteady_newton_solve(dt);

      //Output solution
      problem.doc_solution(doc_info,trace_file);

      //Increment counter for solutions
      doc_info.number()++;
    }

  // Close trace file
  trace_file.close();

  return 0;

} // end of main
