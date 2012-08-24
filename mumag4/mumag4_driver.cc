#ifndef OOMPH_COMBINED_HYBRID_DRIVER_H
#define OOMPH_COMBINED_HYBRID_DRIVER_H

/*
  description of file goes here
*/

// Include all my stuff
#include "../my_general_header.h"
#include "../hybrid_micromagnetics_problem.h"

// Mesh
#include "meshes/simple_cubic_tet_mesh.h"

using namespace oomph;
using namespace MathematicalConstants;

namespace Inputs
{

  double applied_field_coeff = 3;

  double dt = 1e-1;
  double tmax = 100;

  bool adaptive_timestepping = 1;
  bool full_jacobian_fd = 0;
  bool GMRES = 0;
  bool debug_output = 0;
  bool midpointmethod = 0;

  // Which field from mumag#4 should we use?
  unsigned field_num = 1;

  // bool debug_parameters = false;
  // double llg_damp = 1;
  // double llg_precess = 1;
  // double exch_c = 1;
  // double magstatic_coeff = 1;
  // double k1 = 0;
  // double sat_mag = 1;

  // Roughly how many elements to use
  const unsigned nx = 5;
  //const unsigned nz = 3;


  // Applied fields from mumag std problem #4, field is normalised automatically
  // within the element.
  void applied_field(const double& t, const Vector<double>& x, Vector<double>& h_app)
  {
    if(Inputs::field_num == 1)
      {
	// µ0Hx=-24.6 mT, µ0Hy= 4.3 mT, µ0Hz= 0.0 mT
	h_app.assign(3,0.0);
	h_app[0] = -24.6e-3 / mag_parameters::mu0;
	h_app[1] = 4.3e-3 / mag_parameters::mu0;
      }
    else if (Inputs::field_num == 2)
      {
	// µ0Hx=-35.5 mT, µ0Hy=-6.3 mT, µ0Hz= 0.0 mT
	h_app.assign(3,0.0);
	h_app[0] = -35.5-3 / mag_parameters::mu0;
	h_app[1] = -6.3e-3 / mag_parameters::mu0;
      }
  }
  void initial_m(const double& t, const Vector<double>& x,
		 Vector<double>& m)
  {
    m.assign(3,0.0);
    m[0] = -1;
    m[1] = -0.1;
    m[2] = -0.1;
    VectorOps::normalise(m);
  }

};

namespace oomph
{

  //======================================================================
  /// A problem class to test the combination of the hybrid BEM/FEM for
  /// magnetostatic fields and the LLG equations for micromagnetics.
  //======================================================================
  template<class BULK_ELEMENT,
	   template<class BULK_ELEMENT,unsigned DIM> class BEM_ELEMENT,
	   unsigned DIM>
  class ThreeDHybridProblem
    : public HybridMicromagneticsProblem< BULK_ELEMENT, BEM_ELEMENT, DIM >
  {

  public:

    /// Constructor
    ThreeDHybridProblem();

    /// Destructor (empty -- once the problem is done with the program is over)
    ~ThreeDHybridProblem(){};

    /// Doc the solution
    void doc_solution(DocInfo& doc_info);

    /// Set initial condition (incl previous timesteps)
    void set_initial_condition();

    /// Dummy increase timestep by a bit if successful newton solve
    double global_temporal_error_norm()
    {
      return 90;
    }

  }; // end of problem class


//======================================================================
/// Constructor
//======================================================================
template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
ThreeDHybridProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
ThreeDHybridProblem()
{

  if(Inputs::midpointmethod)
    {
      if(Inputs::adaptive_timestepping)
	{
	  std::cerr << "No adaptive midpoint yet." << std::endl;
	  exit(1);
	}
      else
	this->add_time_stepper_pt(new BDF<1>);
    }
  else // use bdf
    {
      if(Inputs::adaptive_timestepping)
	this->add_time_stepper_pt(new BDF<2>(true));
      else
	this->add_time_stepper_pt(new BDF<2>);
    }

  if(Inputs::full_jacobian_fd)
    {
      this->linear_solver_pt() = new FD_LU;
    }
  else if (Inputs::GMRES)
    {
      this->linear_solver_pt() = new GMRES<SumOfMatrices>;

      // // Set general preconditioner
      // IterativeLinearSolver* it_lin_solver_pt =
      //   dynamic_cast<IterativeLinearSolver*>(this->linear_solver_pt());
      // it_lin_solver_pt->preconditioner_pt() =
      //   new ILUZeroPreconditioner<CRDoubleMatrix>;
      //??ds can't use it because we have this sumofmatrices class
    }
  else
    {
      this->linear_solver_pt() = new SuperLUSolver;
    }

  // unsigned lx = 125, ly = 500, lz = 30;
  // unsigned nx = Inputs::nx, nz = 5, ny = Inputs::nx*unsigned((ly/lx));

  unsigned lx = 10, ly = 10, lz = 10;
  unsigned nx = Inputs::nx, nz = nx * unsigned(lz/lx), ny = Inputs::nx*unsigned(ly/lx);

  // Cubeoid
  this->bulk_mesh_pt() = new SimpleCubicTetMesh<BULK_ELEMENT>
    (nx,ny,nz,lx,ly,lz,this->time_stepper_pt());

  // For some reason we have to do this manually...
  this->bulk_mesh_pt()->setup_boundary_element_info();

  MagneticParameters* mumag4_parameters_pt = new MagneticParameters;
  mumag4_parameters_pt->set_mumag4();
  this->magnetic_parameters_pt() = mumag4_parameters_pt;

  // if(Inputs::debug_parameters)
  //   {
  //     mumag4_parameters_pt->exchange_constant() = Inputs::exch_c;
  //     mumag4_parameters_pt->k1() = Inputs::k1;
  //     mumag4_parameters_pt->gilbert_damping() = Inputs::llg_damp;
  //     mumag4_parameters_pt->gamma() = Inputs::llg_precess;
  //     mumag4_parameters_pt->saturation_magnetisation() = Inputs::sat_mag;
  //   }

  //  Create the map of sharp corner nodes and their angles for this mesh.
  // ============================================================

  std::map<Node*,double>* Corners_map_pt = new std::map<Node*,double>;
  this->corners_map_pt() = Corners_map_pt;

  for(unsigned b=0; b <this->bulk_mesh_pt()->nboundary(); b++)
    {
      for(unsigned nd=0; nd<this->bulk_mesh_pt()->nboundary_node(b); nd++)
	{
	  Node* nd_pt = this->bulk_mesh_pt()->boundary_node_pt(b,nd);
	  std::set<unsigned>* boundaries_pt;
	  nd_pt->get_boundaries_pt(boundaries_pt);
	  if(boundaries_pt->size() == DIM)
	    this->corners_map_pt()->insert(std::pair<Node*,double>(nd_pt,0.125));
	}
    }
  std::cout << *Corners_map_pt << std::endl;


  // Loop over elements in bulk mesh to set applied field function pointer
  for(unsigned i=0; i< this->bulk_mesh_pt()->nelement(); i++)
    {
      // Upcast from GeneralisedElement to the present element
      BULK_ELEMENT* elem_pt = dynamic_cast<BULK_ELEMENT*>(this->bulk_mesh_pt()->element_pt(i));

      elem_pt->applied_field_pt() = &Inputs::applied_field;
    }

  // Set up boundary element method and flux conditions on phi_1
  this->finish_building_hybrid_problem();

} // end of constructor


  //======================================================================
  /// Output function
  //======================================================================
template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
void ThreeDHybridProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
doc_solution(DocInfo& doc_info)
{
  // Number of plot points
  unsigned npts=2;

  double cts_time = this->time_pt()->time();

  // File set up
  std::ofstream some_file;
  char filename[100];
  sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
	  doc_info.number());

  // Output
  some_file.open(filename);
  this->mesh_pt()->output(some_file,npts);
  some_file.close();

  // Output means
  Vector<double> means;
  this->get_mean_bulk_values(means);

  std::ofstream means_file;
  means_file.open("./means",std::fstream::app);
  means_file << cts_time;
  for(unsigned i=0; i< means.size(); i++)
    means_file << " " << means[i];
  means_file << std::endl;
  means_file.close();

} // end of doc


  //======================================================================
  /// Set up the initial conditions
  //======================================================================
template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
void ThreeDHybridProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
set_initial_condition()
{

  // Backup time in global Time object
  double backed_up_time=this->time_pt()->time();

  // Past history needs to be established for t=time0-deltat, ...
  // Then provide current values (at t=time0) which will also form
  // the initial guess for the first solve at t=time0+deltat

  // Get M indicies
  Vector<unsigned> m_index_micromag(3,0);
  BULK_ELEMENT* elem_pt = dynamic_cast< BULK_ELEMENT* >(this->bulk_mesh_pt()->element_pt(0));
  for(unsigned i=0; i<3; i++)
    m_index_micromag[i] = elem_pt->m_index_micromag(i);

  //Find number of nodes in mesh
  unsigned num_nod = this->mesh_pt()->nnode();

  // Set continuous times at previous timesteps:
  int nprev_steps=this->time_stepper_pt()->nprev_values();
  Vector<double> prev_time(nprev_steps+1);
  for (int t=nprev_steps;t>=0;t--)
    {
      prev_time[t]=this->time_pt()->time(t);
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
	  // Get initial value of m from inputs
	  //??ds encapsulate properly
	  Vector<double> m(3,0.0), x(DIM,0.0);
	  this->mesh_pt()->node_pt(n)->position(t,x);
	  Inputs::initial_m(time,x,m);

	  // Set initial condition on m
	  for(unsigned i=0; i<3; i++)
	    this->mesh_pt()->node_pt(n)->set_value(t,m_index_micromag[i],m[i]);
	}
    }

  // Reset backed up time for global timestepper
  this->time_pt()->time()=backed_up_time;
}


} // End of oomph namespace

int main(int argc, char* argv[])
{

  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  // Inputs
  const unsigned dim = 3;
  const unsigned nnode_1d = 2;
  double dt = Inputs::dt;
  const double tmax = Inputs::tmax;

  // Dummy error for timestepper - always be ok
  const double dummy_t_eps = 100;


  // Create the problem
  ThreeDHybridProblem< TMicromagElement <dim,nnode_1d>, MicromagFaceElement, dim >
    problem;

  // Initialise timestep, initial conditions
  problem.initialise_dt(dt);
  problem.set_initial_condition();

  // Set up outputs and output initial conditions
  DocInfo doc_info;
  doc_info.set_directory("results");
  doc_info.number()=0;
  problem.doc_solution(doc_info);
  doc_info.number()++;

  // Write paramters to std::out
  problem.magnetic_parameters_pt()->output(std::cout);

  // Additional output?
  if(Inputs::debug_output)
    problem.debug_doc().enable_doc();

  /// Check problem
  if(!(problem.self_test()==0))
    throw OomphLibError("Problem self_test failed","main",
  			OOMPH_EXCEPTION_LOCATION);

  std::cout << "constructor done, everything ready" << "\n" << std::endl;

  // Open a trace file
  std::ofstream trace_file;
  char trace_filename[100];
  sprintf(trace_filename,"%s/trace.dat",doc_info.directory().c_str());
  trace_file.open(trace_filename);


  if(Inputs::adaptive_timestepping)
    {
      // Adaptive while loop

      while (problem.time_pt()->time()<tmax)
	{
	  std::cout << "Time is " << problem.time_pt()->time()<< std::endl
		    << "Current timestep is " << dt << std::endl << std::endl;


	  // Take an adaptive timestep -- the input dt is the suggested timestep.
	  // The adaptive timestepper will adjust dt until the required error
	  // tolerance is satisfied. The function returns a suggestion
	  // for the timestep that should be taken for the next step. This
	  // is either the actual timestep taken this time or a larger
	  // value if the solution was found to be "too accurate".
	  double dt_next=problem.adaptive_unsteady_newton_solve(dt,dummy_t_eps);

	  // Use dt_next as suggestion for the next timestep
	  dt=dt_next;

	  //Output solution
	  problem.doc_solution(doc_info);

	  trace_file << doc_info.number() << " " << problem.time_pt()->time()
		     << " " << dt_next << std::endl;

	  //Increment counter for solutions
	  doc_info.number()++;
	  problem.debug_doc().next_timestep();

	} // end of timestepping loop


    }
  else
    {
      unsigned nstep = int(tmax/dt);

      // Standard timestepping loop
      for(unsigned istep=0; istep<nstep; istep++)
	{
	  std::cout << "Timestep " << istep << std::endl;

	  // Take timestep
	  problem.unsteady_newton_solve(dt);

	  //Output solution
	  problem.doc_solution(doc_info);

	  //Increment counter for solutions
	  doc_info.number()++;
	  problem.debug_doc().next_timestep();
	}

    }
}

#endif
