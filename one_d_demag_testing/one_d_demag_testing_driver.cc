
# include "../micromagnetics_element.h"
# include "./parameters3.cc"

using namespace std;

//==start_of_problem_class============================================
/// 1D Micromag problem in unit interval.
//====================================================================
template<class ELEMENT> 
class OneDMicromagProblem : public Problem
{

private:

  /// Pointer to Poisson source function
  MicromagEquations<1>::PoissonSourceFctPt Source_fct_pt;

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

  /// Pointer to exact M solution function
  MicromagEquations<1>::ExactMFctPt Exact_m_fct_pt;

  /// Pointer to exact phi solution function
  MicromagEquations<1>::ExactPhiFctPt Exact_phi_fct_pt;

  /// Pointer to control node at which the solution is documented ??ds - not sure what this is
  Node* Control_node_pt;

  // Doc info object
  DocInfo Doc_info;
  
  // Trace file
  ofstream Trace_file;

public:

  /// Constructor: Pass number of elements and pointer to source function
  OneDMicromagProblem(const unsigned& n_element, MicromagEquations<1>::PoissonSourceFctPt source_fct_pt, MicromagEquations<1>::LlgSourceFctPt, MicromagEquations<1>::AppliedFieldFctPt applied_field_fct_pt, MicromagEquations<1>::CrystAnisFieldFctPt cryst_anis_field_fct_pt, MicromagEquations<1>::SatMagFctPt sat_mag_fct_pt, MicromagEquations<1>::LlgDampFctPt llg_damp_fct_pt, MicromagEquations<1>::LlgPrecessFctPt llg_precess_fct_pt, MicromagEquations<1>::ExactMFctPt exact_m_fct_pt, MicromagEquations<1>::ExactPhiFctPt exact_phi_fct_pt);

  /// Destructor (empty -- all the cleanup is done in the base class)
  ~OneDMicromagProblem(){};

  /// Update the problem specs before solve: Set boundary conditions
  void actions_before_newton_solve(){};

  /// Update the problem specs after solve (calculate demag field)
  void actions_after_newton_solve(const unsigned n_element);

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

/// Set number of values stored at each node (4: phi, M_x, M_y, M_z)
template<unsigned DIM, unsigned NNODE_1D>
const unsigned QMicromagElement<DIM,NNODE_1D>::Initial_Nvalue = 4;



//=====start_of_constructor===============================================
/// ??ds
//========================================================================
template<class ELEMENT> 
OneDMicromagProblem<ELEMENT>::OneDMicromagProblem(const unsigned& n_element,
						  MicromagEquations<1>::PoissonSourceFctPt source_fct_pt,
						  MicromagEquations<1>::LlgSourceFctPt llg_source_fct_pt,
						  MicromagEquations<1>::AppliedFieldFctPt applied_field_fct_pt,
						  MicromagEquations<1>::CrystAnisFieldFctPt cryst_anis_field_fct_pt,
						  MicromagEquations<1>::SatMagFctPt sat_mag_fct_pt,
						  MicromagEquations<1>::LlgDampFctPt llg_damp_fct_pt,
						  MicromagEquations<1>::LlgPrecessFctPt llg_precess_fct_pt,
						  MicromagEquations<1>::ExactMFctPt exact_m_fct_pt,
						  MicromagEquations<1>::ExactPhiFctPt exact_phi_fct_pt) :
  Source_fct_pt(source_fct_pt), Llg_source_fct_pt(llg_source_fct_pt), Applied_field_fct_pt(applied_field_fct_pt), Cryst_anis_field_fct_pt(cryst_anis_field_fct_pt), Sat_mag_fct_pt(sat_mag_fct_pt), Llg_damp_fct_pt(llg_damp_fct_pt), Llg_precess_fct_pt(llg_precess_fct_pt), Exact_m_fct_pt(exact_m_fct_pt), Exact_phi_fct_pt(exact_phi_fct_pt)
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
  cout << "Recording trace of the solution at: " << Control_node_pt->x(0) << std::endl;


  // Set up the boundary conditions for this problem: pin the nodes at either end
  mesh_pt()->boundary_node_pt(0,0)->pin(0);
  mesh_pt()->boundary_node_pt(1,0)->pin(0);


  // Loop over elements to set pointers to source function, applied field and time
  for(unsigned i=0;i<n_element;i++)
    {
      // Upcast from GeneralisedElement to the present element
      ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
      //Set the source function pointer
      elem_pt->source_fct_pt() = Source_fct_pt;

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

  // Get pointer to (0th) element - exact element doesn't matter for this use, hopefully!
  //??ds this will break if number of nodal data values changes in different elements, don't think it does change though
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
	  // Set boundary conditions at this node
	  Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
	  Vector<double> x(1,nod_pt->x(0));
	  double phi_value = OneDMicromagSetup::boundary_phi(t,x);
	  nod_pt->set_value(phi_nodal_index,phi_value);
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
 
  // Vector of exact solution value
  Vector<double> M(3);
  Vector<double> x(1);

  //Find number of nodes in mesh
  unsigned num_nod = mesh_pt()->nnode();

  // Get pointer to an element (any will do so take 0th)
  ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0)); 

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
      cout << "setting IC at time =" << time << std::endl;
   
      // Loop over the nodes to set initial values everywhere
      for (unsigned n=0;n<num_nod;n++)
	{
	  // Get nodal coordinate
	  x[0]=mesh_pt()->node_pt(n)->x(0);

	  // Get initial value of M
	  OneDMicromagSetup::initial_M(x,M);
     
	  // Assign solution
	  for(unsigned i=0; i<3; i++)
	    {
	      // Set ith direction of M on node n at time t to be M[i]
	      mesh_pt()->node_pt(n)->set_value(t,elem_pt->M_index_micromag(i),M[i]);
	    }
     
	  // Loop over coordinate directions: Mesh doesn't move, so previous position = present position
	  // ??ds presumably this is where the ALE formulation would/will/should come in
	  // for (unsigned i=0;i<1;i++)
	  //   {
	  //     mesh_pt()->node_pt(n)->x(t,i)=x[i];
	  //   }
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

  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts;
  npts=5; 

  cout << std::endl;
  cout << "=================================================" << std::endl;
  cout << "Doc'ing solution for t=" << time_pt()->time() << std::endl;
  cout << "=================================================" << std::endl;

  // Output solution 
  //-----------------
  sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
	  doc_info.number());

  some_file.open(filename);
  mesh_pt()->output(some_file,npts);

  // // Write file as a tecplot text object
  // some_file << "TEXT X=2.5,Y=93.6,F=HELV,HU=POINT,C=BLUE,H=26,T=\"time = " 
  // 	    << time_pt()->time() << "\"";
  // // ...and draw a horizontal line whose length is proportional
  // // to the elapsed time
  // some_file << "GEOMETRY X=2.5,Y=98,T=LINE,C=BLUE,LT=0.4" << std::endl;
  // some_file << "1" << std::endl;
  // some_file << "2" << std::endl;
  // some_file << " 0 0" << std::endl;
  // some_file << time_pt()->time()*20.0 << " 0" << std::endl;

  some_file.close();
} // end of doc

 

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////



  //======start_of_main==================================================
  /// Driver for 1D Micromag problem
  //=====================================================================
int main()
{

  // Set up the problem:
  unsigned n_element = OneDMicromagSetup::n_x_elements; //Number of elements
  OneDMicromagProblem<QMicromagElement<1,2> > problem(n_element, OneDMicromagSetup::source_function,
						      OneDMicromagSetup::llg_source_function,
						      OneDMicromagSetup::applied_field, 
						      OneDMicromagSetup::cryst_anis_field, 
						      OneDMicromagSetup::sat_mag, 
						      OneDMicromagSetup::llg_damping_coeff, 
						      OneDMicromagSetup::llg_precession_coeff,
						      OneDMicromagSetup::exact_M_solution,
						      OneDMicromagSetup::exact_phi_solution);


  // SET UP OUTPUT
  // Setup labels for output
  DocInfo doc_info;

  // Output directory
  doc_info.set_directory("RESLT");
 
  // Output number
  doc_info.number()=0;
 
  // Open a trace file
  ofstream trace_file;
  char filename[100];   
  sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
  trace_file.open(filename);
  trace_file << "VARIABLES=\"time\",\"u<SUB>FE</SUB>\","
	     << "\"u<SUB>exact</SUB>\",\"norm of error\",\"norm of solution\""
	     << std::endl;
  

  // SET UP TIME STEPPING
  // Choose simulation interval and timestep
  double t_max = OneDMicromagSetup::t_max;
  double dt = OneDMicromagSetup::dt;

  // Initialise timestep -- also sets the weights for all timesteppers
  // in the problem.
  problem.initialise_dt(dt);
 
  // Set initial condition (on M)
  problem.set_initial_condition();

  // Output initial conditions
  problem.doc_solution(doc_info,trace_file);   //Output initial condition
  doc_info.number()++; //Increment counter for solutions 

  // Check whether the problem can be solved
  cout << "\n\n\nProblem self-test ";
  if (problem.self_test()==0)  
    {
      cout << "passed: Problem can be solved." << std::endl;
    }
  else 
    {
      throw OomphLibError("failed!",
			  "main()",
			  OOMPH_EXCEPTION_LOCATION);
    }

  // SOLVE THE PROBLEM
  // Find number of steps
  unsigned nstep = unsigned(t_max/dt);

  // Timestepping loop
  for (unsigned istep=0;istep<nstep;istep++)
    {
      cout << "Timestep " << istep << std::endl;
   
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
