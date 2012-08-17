#ifndef OOMPH_COMBINED_HYBRID_DRIVER_H
#define OOMPH_COMBINED_HYBRID_DRIVER_H

/*
  description of file goes here
*/

// Include all my stuff
#include "../my_general_header.h"
#include "../hybrid_micromagnetics_problem.h"

// Mesh
#include "meshes/tetgen_mesh.h"

#include "meshes/rectangular_quadmesh.h"
#include "meshes/simple_cubic_mesh.h"

using namespace oomph;
using namespace MathematicalConstants;

/*

  problems with jacobian:

  d m_y/ d m_x

  d m_z/d m_x

  with m_x at different nodes

  affected by dt


  Also flux element Jacobian is wrong, but not badly so
*/
namespace Inputs
{
  // ??ds need to do this input stuff properly..
  double magnetostatic_coeff = 1;
  double applied_field_coeff = 3;

  double dt = 1e-3;
  double tmax = 1;

  bool adaptive_timestepping = 0;
  bool full_jacobian_fd = 0;
  unsigned sum_matrix = 1;
  bool GMRES = 0;

  bool midpointmethod = 0;
  const unsigned bdforder = 2;

  //??ds temp - commented corner calculation code so only valid for sphere!

  // Roughly how many elements per side, scaled as required. The nmag cubeoid
  // example uses ~10 (although they are tets).
  const unsigned nx = 3;

  const unsigned dim = 3;
  // If changing dim you also probably need to swap element type
  // (QMicromagElement vs TMicromagElement) because the mesh changes type.


  // Remember these might not be in Jacobian yet!
  void applied_field(const double& t, const Vector<double>& x, Vector<double>& h_app)
  {
    h_app.assign(3,0.0);
    h_app[0] = +1;
    // h_app[1] = 0.01;
    // h_app[2] = 0.01;
    for(unsigned j=0; j<3; j++)
      h_app[j] *= Inputs::applied_field_coeff;
  }

  void initial_m(const double& t, const Vector<double>& x,
		 Vector<double>& m)
  {
    m.assign(3,0.0);
    m[0] = -1;
    m[1] = -0.1;
    m[2] = -0.1;
    normalise(m);

    // m[0] = sin(x[0])*sin(x[0]);
    // m[1] = cos(x[0])*cos(x[0]);
    // m[2] = x[1];
    // normalise(m);

    // double rsq = pow(x[0],2) +  pow(x[1],2) + pow(x[2],2);
    // m[0] = rsq*sin(x[0]);
    // m[1] = rsq*cos(x[0]);
    // m[2] = rsq;
    // //normalise(m);
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
    ThreeDHybridProblem(const std::string& node_file_name,
			const std::string& element_file_name,
			const std::string& face_file_name);

    /// Destructor (empty -- once the problem is done with the program is over)
    ~ThreeDHybridProblem(){};

    /// Doc the solution
    void doc_solution(DocInfo& doc_info);

    /// Set initial condition (incl previous timesteps)
    void set_initial_condition();

    /// Dummy - always reduce timestep as much as possible if fails to converge.
    double global_temporal_error_norm()
    {
      return 1e-8;
    }

    // //?? check if node is corner, nasty!
    // bool corner(Vector<double> x)
    // {
    //   bool temp(1);
    //   for(unsigned i=0; i<DIM; i++)
    // 	{
    // 	  temp = temp && ( (x[i] == 0) || (x[i] == 1) );
    // 	}
    //   return temp;
    // }

    // /// Overload get_jacobian to include the boundary matrix in sumofmatrices
    // /// format. Note that this will only work with iterative solvers since we
    // /// can only multiply when using sumofmatrices.
    // void get_jacobian(DoubleVector& residuals, SumOfMatrices& jacobian)
    // {
    //   std::cout << "Calling your get_jacobian function using SumOfMatrices." << std::endl;

    //   // Create a matrix to store the sparse part of the Jacobian and get it
    //   CRDoubleMatrix* sparse_jacobian_pt = new CRDoubleMatrix;
    //   Problem::get_jacobian(residuals,*sparse_jacobian_pt);

    //   // Set as the main (first) matrix of the sum.
    //   jacobian.main_matrix_pt() = sparse_jacobian_pt;

    //   // Set the sparse part of the Jacobian to be deleted along with the sum of
    //   // the matrices - avoid a memory leak.
    //   jacobian.set_delete_main_matrix();

    //   if(Inputs::magnetostatic_coeff != 0)
    // 	{

    // 	  // Add the boundary element matrix to the total Jacobian. It represents
    // 	  // the derivative of phi with respect to phi_1 so each entry goes in the
    // 	  // phi row and the phi_1 column of the respective element (this is done
    // 	  // via the two maps). Don't delete when done.
    // 	  jacobian.add_matrix(boundary_matrix_pt(),
    // 			      &Global_boundary_equation_num_map,
    // 			      &Global_phi_1_num_map,0);
    // 	}

    //   // //dump jacobian for tests
    //   // std::ofstream matrix_file;
    //   // matrix_file.precision(16);
    //   // char filename[100];
    //   // sprintf(filename,"matrices/sm_jacobian");
    //   // matrix_file.open(filename);
    //   // jacobian.sparse_indexed_output(matrix_file);
    //   // matrix_file.close();

    //   // residuals.output("matrices/residual");


    //   //jacobian.Matrix::sparse_indexed_output("sum_jacobian");

    //   //exit(0);
    // }






  }; // end of problem class


//======================================================================
/// Constructor
//======================================================================
template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
ThreeDHybridProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
ThreeDHybridProblem(const std::string& node_file_name,
		    const std::string& element_file_name,
		    const std::string& face_file_name)
{

  if(Inputs::midpointmethod)
    {
      if(Inputs::adaptive_timestepping)
	std::cerr << "No adaptive midpoint yet." << std::endl;

      this->add_time_stepper_pt(new BDF<1>);
    }
  else // use bdf
    {
      if(Inputs::adaptive_timestepping)
	this->add_time_stepper_pt(new BDF<Inputs::bdforder>(true));
      else
	this->add_time_stepper_pt(new BDF<Inputs::bdforder>);
    }

  if(Inputs::full_jacobian_fd)
    {
      this->linear_solver_pt() = new FD_LU;
    }
  else if (Inputs::GMRES)
    {
      if(Inputs::sum_matrix == 1)
	this->linear_solver_pt() = new GMRES<SumOfMatrices>;
      else
	this->linear_solver_pt() = new GMRES<CRDoubleMatrix>;

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

  unsigned lx = 3, ly = 3, lz = 10;
  unsigned ny = Inputs::nx, nz = unsigned(Inputs::nx * lz/lx);
  // Cubeoid
  this->bulk_mesh_pt() = new SimpleCubicMesh<BULK_ELEMENT>
    (Inputs::nx,ny,nz,lx,ly,lz,this->time_stepper_pt());



  // Bulk elements
  //------------------------------------------------------------

  // Loop over elements in bulk mesh to set function pointers
  for(unsigned i=0; i< this->bulk_mesh_pt()->nelement(); i++)
    {
      // Upcast from GeneralisedElement to the present element
      BULK_ELEMENT* elem_pt = dynamic_cast<BULK_ELEMENT*>(this->bulk_mesh_pt()->element_pt(i));

      // Set pointer to continuous time
      elem_pt->time_pt() = this->time_pt();

      // // Set the function pointers for parameters
      // //??ds fix this to use proper encapsulation asap
      // elem_pt->applied_field_pt() = &Inputs::applied_field;
      // elem_pt->cryst_anis_field_pt() = &Inputs::cryst_anis_field;
      // elem_pt->hca_derivative_pt() = &Inputs::dhcadm_k;
      // // elem_pt->sat_mag_pt() = &Inputs::sat_mag;
      // elem_pt->llg_damp_pt() = &Inputs::llg_damping;
      // elem_pt->llg_precess_pt() = &Inputs::llg_precession;
      // elem_pt->exchange_coeff_pt() = &Inputs::exchange_coeff;
      // elem_pt->magnetostatic_coeff_pt() = &Inputs::magnetostatic_coeff;
    }



  // // ??ds temp: pin the values on
  // if((Inputs::magnetostatic_coeff == 0))
  //   {
  // 	BULK_ELEMENT* some_el_pt = dynamic_cast< BULK_ELEMENT* >
  // 	  (bulk_mesh_pt()->element_pt(0));
  // 	for(unsigned b=0; b<bulk_mesh_pt()->nboundary(); b++)
  // 	  {
  // 	    for(unsigned nd=0; nd < bulk_mesh_pt()->nboundary_node(b); nd++)
  // 	      {
  // 		bulk_mesh_pt()->boundary_node_pt(b,nd)->
  // 		  pin(some_el_pt->phi_index_micromag());
  // 	      }
  // 	  }

  //pin_local_eqn_on_all_boundaries(phi_index,  bulk_mesh_pt())
  //}

  // //??temp to help with testing phi_1 pin it's value to zero at r cos(azi) sin(polar) = 0,
  // // i.e when r =0.
  // bool found = false;
  // for(unsigned nd=0; nd< bulk_mesh_pt()->nnode(); nd++)
  //   {
  // 	Vector<double> nd_x(DIM,0.0);
  // 	for(unsigned j=0; j<DIM; j++)
  // 	  nd_x[j] = bulk_mesh_pt()->node_pt(nd)->x(j);
  // 	if ( small(nd_x[0]) && small(nd_x[1]) && small(1 - nd_x[2]) )
  // 	  {
  // 	    bulk_mesh_pt()->node_pt(nd)->pin(phi_1_index());
  // 	    bulk_mesh_pt()->node_pt(nd)->set_value(phi_1_index(),0.0);
  // 	    found = true;
  // 	    break;
  // 	  }
  //   }
  // if (!(found))
  //   throw OomphLibError("No node near middle","",OOMPH_EXCEPTION_LOCATION);


  // Flux elements
  //------------------------------------------------------------

  // We want Neumann (flux) boundary condtions on phi_1 on all boundaries so
  // create the face elements needed.
  for(unsigned b=0; b < this->bulk_mesh_pt()->nboundary(); b++)
    {
      create_flux_elements(b,this->bulk_mesh_pt(),this->flux_mesh_pt());
    }

  // BEM elements
  //------------------------------------------------------------

  // Create integration scheme in case it is needed.
  QVariableOrderGaussLegendre<DIM-1>* variable_scheme_pt
    = new QVariableOrderGaussLegendre<DIM-1>;

  // Create BEM elements on all boundaries and add to BEM mesh
  build_bem_mesh(this->bem_mesh_pt());

  // Set boundary mesh pointer in element (static memeber so only do once)
  BEM_ELEMENT<BULK_ELEMENT,DIM>::set_boundary_mesh_pt(this->bem_mesh_pt());

  // Set pointers in elements
  for(unsigned i_ele = 0; i_ele < this->bem_mesh_pt()->nelement(); i_ele++)
    {
      // Cast pointer
      BEM_ELEMENT<BULK_ELEMENT,DIM>* ele_pt
	= dynamic_cast< BEM_ELEMENT< BULK_ELEMENT,DIM>* >
	(this->bem_mesh_pt()->element_pt(i_ele));

      // Set integration scheme
      ele_pt->set_integration_scheme(variable_scheme_pt);
    }

  // Build global finite element mesh
  this->add_sub_mesh(this->bulk_mesh_pt());
  this->add_sub_mesh(this->flux_mesh_pt());
  this->build_global_mesh();

  // Setup equation numbering scheme for all the finite elements
  std::cout << "FEM number of equations: " << this->assign_eqn_numbers() << std::endl;

  // Make the boundary matrix (including setting up the numbering scheme).  Note
  // that this requires the FEM numbering scheme to be already set up.
  this->build_boundary_matrix();


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

  // File set up
  std::ofstream some_file;
  char filename[100];
  sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
	  doc_info.number());

  // Output
  some_file.open(filename);
  this->mesh_pt()->output(some_file,npts);
  some_file.close();
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

  // Store command line arguments
  CommandLineArgs::setup(argc,argv);

  if(argc < 4)
    {
      throw OomphLibError("Not enough args: needs mesh inputs",
			  "main()", OOMPH_EXCEPTION_LOCATION);
    }

  // Inputs
  const unsigned dim = Inputs::dim;
  const unsigned nnode_1d = 2;
  // const double dt = 1e-3;
  // const unsigned nstep = 10;
  double dt = Inputs::dt;
  const double tmax = Inputs::tmax;

  // Dummy error for timestepper - always be ok
  const double dummy_t_eps = 100;



  // Create the problem
  ThreeDHybridProblem< TMicromagElement <dim,nnode_1d>, MicromagFaceElement, dim >
    problem(argv[1],argv[2],argv[3]);

  // problem.max_newton_iterations() = 15;
  // problem.max_residuals() = 30;

  // // dump mesh for testing
  // std::ofstream mesh_plot;
  // mesh_plot.open("./mesh_points");
  // for(unsigned nd=0; nd<problem.mesh_pt()->nnode(); nd++)
  //   {
  //     for(unsigned j=0; j<dim; j++)
  // 	mesh_plot << problem.mesh_pt()->node_pt(nd)->x(j) << " ";
  //     mesh_plot << std::endl;
  //   }
  // mesh_plot.close();


  // // dump boundary for testing
  // unsigned b = 0;
  // std::ofstream bound_plot;
  // bound_plot.open("./bound_points");
  // for(unsigned nd=0; nd<problem.bulk_mesh_pt()->nboundary_node(b); nd++)
  //   {
  //     for(unsigned j=0; j<dim; j++)
  // 	bound_plot << problem.bulk_mesh_pt()->boundary_node_pt(b,nd)->x(j) << " ";
  //     bound_plot << std::endl;
  //   }
  // bound_plot.close();


  // Initialise timestep, initial conditions
  problem.initialise_dt(dt);
  problem.set_initial_condition();

  // Set up outputs and output initial conditions
  DocInfo doc_info;
  doc_info.set_directory("results");
  doc_info.number()=0;
  problem.doc_solution(doc_info);
  doc_info.number()++;

  /// Check problem
  if(!(problem.self_test()==0))
    throw OomphLibError("Problem self_test failed","main",
  			OOMPH_EXCEPTION_LOCATION);

  std::cout << "constructor done, everything ready" << "\n" << std::endl;

  std::cout << std::endl;

  // Open a trace file
  std::ofstream trace_file;
  char trace_filename[100];
  sprintf(trace_filename,"%s/trace.dat",doc_info.directory().c_str());
  trace_file.open(trace_filename);

  // DoubleVector residual_vector;
  // problem.get_residuals(residual_vector);
  // residual_vector.output("./matrices/residuals");

  problem.boundary_matrix_pt()->Matrix::output("./matrices/bem");

  // DoubleVector dummy_res;
  // DenseMatrix<double> fd_jacobian;
  // problem.get_fd_jacobian(dummy_res,fd_jacobian);
  // fd_jacobian.Matrix::output("fd_jacobian");





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

	  // DenseDoubleMatrix jacobian;
	  // DoubleVector residuals;
	  // problem.get_jacobian(residuals,jacobian);

	  // std::ofstream residual_file;
	  // residual_file.precision(16);
	  // char filename2[100];
	  // sprintf(filename2,"results/residual");
	  // residual_file.open(filename2);
	  // residuals.output(residual_file);
	  // residual_file.close();

	  //Increment counter for solutions
	  doc_info.number()++;
	}

    }
}

#endif

