
// Include the generic routines
#include "generic.h"

// Include the micromagnetics equations
#include "../micromagnetics_element.h"

// Include the mesh
#include "meshes/simple_rectangular_quadmesh.h"

using namespace std;

//==start_of_namespace================================================
/// Namespace for various functions
//====================================================================
namespace TwoDMicromagSetup
{
  void source_function(const Vector<double>& x, double& source)
  {
    // Just set to zero unless needed for testing
    source = 0.0;
  }

  void get_initial_M(const Vector<double>& x, Vector<double>& M)
  {
    // Start with step change in M_x from 1 to -1 at x = 0.5
    // Should create a domain wall type structure
    M[0] = 1.0; //tanh(5*(x[0] - 0.5));

    // Initialise y and z components of M to zero
    M[1] = 0.0;
    M[2] = 0.0;
  }

  void get_applied_field(const Vector<double>& x, Vector<double>& H_applied)
  {
    H_applied[0] = -1.0;
    H_applied[1] = 0.0;
    H_applied[2] = 0.0;
  }

  double get_boundary_phi(const Vector<double>& x)
  {
    // Set all boundaries to zero for now
    return 0.0;
  }

  double get_llg_damping_coeff()
  {
    //??ds fill in here and pass pointers out like with source
    return 1;
  }

  double get_llg_precession_coeff()
  {
    //??ds fill in here and pass pointers out like with source
    return 0.5;
  }

  void get_cryst_anis_field(const Vector<double>& x, Vector<double>& H_cryst_anis)
  {
    H_cryst_anis[0] = 0.5;
    
    H_cryst_anis[1] = 0;

    H_cryst_anis[2] = 0;
  }

}; // End of namespace



//==start_of_problem_class============================================
/// 1D Micromag problem in unit interval.
//====================================================================
template<class ELEMENT> 
class TwoDMicromagProblem : public Problem
{

private:

  /// Dimension
  unsigned const static Element_dim = 2;

  /// Pointer to source function
  MicromagEquations<2>::PoissonSourceFctPt Source_fct_pt;

  /// Pointer to applied field function
  MicromagEquations<2>::AppliedFieldFctPt Applied_field_fct_pt;

  /// Pointer to crystalline anisotropy effective field function
  MicromagEquations<2>::CrystAnisFieldFctPt Cryst_anis_field_fct_pt;

  /// Pointer to control node at which the solution is documented ??ds - not sure what this is
  Node* Control_node_pt;

  // Doc info object
  DocInfo Doc_info;
  
  // Trace file
  ofstream Trace_file;

public:

  /// Constructor: Pass number of elements and pointer to source functions
  TwoDMicromagProblem(const unsigned& n__x,
		      const unsigned& n__y,
		      MicromagEquations<Element_dim>::PoissonSourceFctPt source_fct_pt,
		      MicromagEquations<Element_dim>::AppliedFieldFctPt applied_field_fct_pt,
		      MicromagEquations<Element_dim>::CrystAnisFieldFctPt cryst_anis_field_fct_pt);

  /// Destructor (empty -- all the cleanup is done in the base class)
  ~TwoDMicromagProblem(){};

  /// Update the problem specs before solve: Set boundary conditions
  void actions_before_newton_solve(){};

  /// Update the problem specs after solve (calculate demag field)
  void actions_after_newton_solve(const unsigned n_element_x, const unsigned n_element_y);

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
TwoDMicromagProblem<ELEMENT>::TwoDMicromagProblem(const unsigned& n_x,
						  const unsigned& n_y,
						  MicromagEquations<Element_dim>::PoissonSourceFctPt source_fct_pt,
						  MicromagEquations<Element_dim>::AppliedFieldFctPt applied_field_fct_pt,
						  MicromagEquations<Element_dim>::CrystAnisFieldFctPt cryst_anis_field_fct_pt) :
  Source_fct_pt(source_fct_pt), Applied_field_fct_pt(applied_field_fct_pt), Cryst_anis_field_fct_pt(cryst_anis_field_fct_pt)
{  
  // Allocate the timestepper -- this constructs the Problem's time object with a sufficient amount of storage to store the previous timsteps. 
  add_time_stepper_pt(new BDF<2>);

  // Set domain length 
  double l_x = 1.0;
  double l_y = 1.0;

  // Build mesh and store pointer in Problem
  Problem::mesh_pt() = new SimpleRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y,time_stepper_pt());
 
  // Choose a control node at which the solution is documented
  unsigned control_el = unsigned(n_x/2); // Pick a control element in the middle
  Control_node_pt=mesh_pt()->finite_element_pt(control_el)->node_pt(0);  // Choose its first node as the control node
  cout << "Recording trace of the solution at: " << Control_node_pt->x(0) << std::endl;


  // Set up the boundary conditions for this problem: pin the nodes at edges
  // Loop over boundarys then over nodes within each boundary
  unsigned n_bound = mesh_pt()->nboundary();
  for(unsigned i=0; i<n_bound; i++)
    {
      unsigned n_node = mesh_pt()->nboundary_node(i);
      for(unsigned j=0; j<n_node; j++)
	{
	  // Pin the jth node on the ith boundary
	  mesh_pt()->boundary_node_pt(i,j)->pin(0);
	}
    }


  // Loop over elements to set pointers to source function, applied field and time
  // (does not need extra loops for extra dimension since we loop over a list of all elements).
  unsigned n_element = mesh_pt()->nelement();
  for(unsigned i=0;i<n_element;i++)
    {
	  // Upcast from GeneralisedElement to the present element
	  ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
	  //Set the source function pointer
	  elem_pt->source_fct_pt() = Source_fct_pt;

	  // Set the applied field function pointer
	  elem_pt->applied_field_fct_pt() = Applied_field_fct_pt;

	  // Set the crystalline anisotropy effective field pointer
	  elem_pt->cryst_anis_field_fct_pt() = Cryst_anis_field_fct_pt;

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
void TwoDMicromagProblem<ELEMENT>::actions_before_implicit_timestep()
{

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
	  // Get pointer to this node
	  Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
	  // Get the Eulerian position of this node
	  Vector<double> x_node(Element_dim,0.0);
	  for(unsigned i=0;  i<Element_dim; i++){x_node[i] = nod_pt->x(i);}
	  // Get the boundary value for this position
	  double phi_value = TwoDMicromagSetup::get_boundary_phi(x_node);
	  // Set the boundary value
	  nod_pt->set_value(phi_nodal_index,phi_value);
	}
    }

}


//======================start_of_set_initial_condition====================
/// \short Set initial condition: Assign previous and current values
/// from exact solution.
//========================================================================
template<class ELEMENT>
void TwoDMicromagProblem<ELEMENT>::set_initial_condition()
{ 
  // Backup time in global Time object
  double backed_up_time=time_pt()->time();
         
  // Past history needs to be established for t=time0-deltat, ...
  // Then provide current values (at t=time0) which will also form
  // the initial guess for the first solve at t=time0+deltat
 
  // Vector of exact solution value
  Vector<double> M(3);
  Vector<double> x(Element_dim);

  //Find number of nodes in mesh
  unsigned num_nod = mesh_pt()->nnode();

  // Get pointer to an element (any will do so take 0th)
  ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0)); 

  // Set continuous times at previous timesteps:
  // How many previous timesteps does the timestepper use?
  int nprev_steps=time_stepper_pt()->nprev_values();

  //??ds why is this a seperate loop, why do we even need to store the time valuse at all?
  // - moved time getting into other loop
  // Vector<double> prev_time(nprev_steps+1);
  // for (int t=nprev_steps;t>=0;t--)
  //   {
  //     prev_time[t]=time_pt()->time(unsigned(t));
  //   }

  cout << nprev_steps << endl;

  // Loop over current & previous timesteps
  for (int t=nprev_steps; t>=0; t--)
    {
      // Continuous time
      double time = time_pt()->time(unsigned(t));
      cout << "setting IC at time =" << time << std::endl;
   
      // Loop over the nodes to set initial values everywhere
      for (unsigned n=0;n<num_nod;n++)
	{
	  // Get Eulerian nodal position
	  for(unsigned i=0; i<Element_dim; i++) 
	    {
	      x[i]=mesh_pt()->node_pt(n)->x(t,i);
	    }

	  // Get initial value of M at the nodal position x
	  TwoDMicromagSetup::get_initial_M(x,M);
     
	  // Assign solution (loop over magnetisation directions)
	  for(unsigned i=0; i<3; i++)
	    {
	      // Set ith direction of M on node n at time t to be M[i]
	      mesh_pt()->node_pt(n)->set_value(t, elem_pt->M_index_micromag(i), M[i]);
	    }
     
	  // // Loop over coordinate directions: Mesh doesn't move, so previous position = present position
	  // // ??ds presumably this is where the ALE formulation would/will/should come in
	  // for (unsigned i=0;i<Element_dim;i++)
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
void TwoDMicromagProblem<ELEMENT>::doc_solution(DocInfo& doc_info, std::ofstream& trace_file)
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
  unsigned n_element_x = 20; //Number of elements
  unsigned n_element_y = 5;

  TwoDMicromagProblem<QMicromagElement<2,2> > problem(n_element_x, n_element_y, TwoDMicromagSetup::source_function, TwoDMicromagSetup::get_applied_field, TwoDMicromagSetup::get_cryst_anis_field);
  //??ds should pass these by reference??


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
  double t_max=60;
  double dt=1;

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
