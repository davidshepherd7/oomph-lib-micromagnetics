#ifndef OOMPH_COMBINED_HYBRID_DRIVER_H
#define OOMPH_COMBINED_HYBRID_DRIVER_H

/*
  description of file goes here
*/

// Include all my stuff
#include "../my_general_header.h"
#include "../hybrid_micromagnetics_problem.h"
#include "../micromagnetics_preconditioners.h"

// Mesh
#include "meshes/simple_cubic_tet_mesh.h"


using namespace oomph;
using namespace MathematicalConstants;


namespace Inputs
{
  enum MMPrecond
    {
      identity,
      diagonal,
      block_diagonal,
      block_upper_triangular,
      block_lower_triangular,
      exact,
      block_exact,
    };

  double dt = 0.1;
  unsigned nx = 10, ny = nx, nz = nx;

  double tmax = 500;

  bool adaptive_timestepping = 0;
  bool full_jacobian_fd = 0; // fails...
  bool GMRES = 1;
  bool debug_output = 0;
  bool midpointmethod = 0;

  enum MMPrecond preconditioner_type = Inputs::block_exact;

  // Just run the solver tests?
  bool solver_tests = 1;

  // Output BEM?
  bool dump_bem = 1;
  unsigned output_timestep = 0;

  // Combine blocks?
  bool block_combine_bulk_and_boundary_phi = 0;
  bool block_combine_bulk_and_boundary_phi_1 = 0;
  bool block_combine_m_directions = 0;

  // bool debug_parameters = false;
  // double llg_damp = 1;
  // double llg_precess = 1;
  // double exch_c = 1;
  // double magstatic_coeff = 1;
  // double k1 = 0;
  // double sat_mag = 1;

  // A cube
  double lx = 10, ly = 10, lz = 10;

  // Applied fields from mumag std problem #4, field is normalised automatically
  // within the element.
  void applied_field(const double& t, const Vector<double>& x, Vector<double>& h_app)
  {
    h_app.assign(3,0.0);
    double h_strength = -1*(80e-3 / mag_parameters::mu0);
    h_app[0] = h_strength;
    h_app[1] = h_strength;
    h_app[2] = h_strength;

  }
  void initial_m(const double& t, const Vector<double>& x,
                 Vector<double>& m)
  {
    m.assign(3,0.0);
    m[0] = 1;
    m[1] = 1;
    m[2] = 1;
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
    ThreeDHybridProblem(const unsigned& a, const unsigned& b,
                        const unsigned& c);

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

    /// Construct the mapping between dofs as assigned in elements and the final
    /// blocks depending on which dofs we want lumped together.
    void build_dof_to_block_map(const bool combine_m_directions,
                                const bool combine_bulk_and_boundary_phi_1,
                                const bool combine_bulk_and_boundary_phi,
                                Vector<unsigned>& dof_to_block_map) const;

  }; // end of problem class


//======================================================================
/// Constructor
//======================================================================
template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
ThreeDHybridProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
ThreeDHybridProblem(const unsigned& a, const unsigned& b, const unsigned& c)
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

  double lx = Inputs::lx, ly = Inputs::ly, lz = Inputs::lz;
  unsigned nx = a, ny = b, nz = c;

  // Cube:
  this->bulk_mesh_pt() = new SimpleCubicTetMesh<BULK_ELEMENT>
    (nx,ny,nz,lx,ly,lz,this->time_stepper_pt());

  // For some reason we have to do this manually...
  this->bulk_mesh_pt()->setup_boundary_element_info();

  // Set up magnetic parameters - using mumag4 ones for now
  MagneticParameters* mumag4_parameters_pt = new MagneticParameters;
  mumag4_parameters_pt->set_mumag4();
  this->magnetic_parameters_pt() = mumag4_parameters_pt;


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

  // Set up the linear solver. This needs to go after all mesh construction
  // etc. because the preconditioners need the mesh to be filled in already.
  if(Inputs::full_jacobian_fd)
    {
      this->linear_solver_pt() = new FD_LU;
    }
  else if(Inputs::GMRES)
    {

      // Set solver to GMRES (on a "SumOfMatrices" because we need to add a
      // dense block for BEM) and cast the pointer ready for later.
      this->linear_solver_pt() = new GMRES<SumOfMatrices>;
      IterativeLinearSolver* it_lin_solver_pt =
        dynamic_cast<IterativeLinearSolver*>(this->linear_solver_pt());

      // Set a preconditioner according to the inputs
      bool block_preconditioner = true;
      Preconditioner* p_pt;
      switch(Inputs::preconditioner_type)
        {
        case Inputs::identity:
          p_pt = new IdentityPreconditioner;
          block_preconditioner = false;
          break;

        case Inputs::diagonal:
          p_pt = new MatrixBasedDiagPreconditioner;
          block_preconditioner = false;
          break;

        case Inputs::block_diagonal:
          p_pt = new MMBlockDiagonalPreconditioner<CRDoubleMatrix>;
          break;

        case Inputs::block_upper_triangular:
          p_pt = new MMBlockTriangularPreconditioner<CRDoubleMatrix>;
          // upper triangular by default
          break;

        case Inputs::block_lower_triangular:
          {
            p_pt = new MMBlockTriangularPreconditioner<CRDoubleMatrix>;

            // We have to convert to lower triangular (via a cast'ed pointer)
            BlockTriangularPreconditioner<CRDoubleMatrix>* btp_pt =
              dynamic_cast <BlockTriangularPreconditioner<CRDoubleMatrix>* >
              (p_pt);
            btp_pt->lower_triangular();
          }
          break;

        case Inputs::block_exact:
          p_pt = new MMExactBlockPreconditioner<CRDoubleMatrix>;
          break;

        case Inputs::exact:
          std::cout << "Preconditioner not implemented." << std::endl;
          exit(3);
          break;

        default:
          std::cout << "Preconditioner not implemented." << std::endl;
          exit(3);
        }

      it_lin_solver_pt->preconditioner_pt() = p_pt;

      if(block_preconditioner)
        {
          // Cast the pointer to a block preconditioner so we can use it's functions
          GeneralPurposeBlockPreconditioner<CRDoubleMatrix>* prec_pt =
            dynamic_cast<GeneralPurposeBlockPreconditioner<CRDoubleMatrix>* >
            (it_lin_solver_pt->preconditioner_pt());

          prec_pt->problem_pt() = this;

          // Add the main mesh (for basic blocking of the Jacobian) and the BEM mesh
          // (to block the boundary values of phi and phi_1 differently) to the
          // preconditioner.
          prec_pt->add_mesh(this->bulk_mesh_pt());
          prec_pt->add_mesh(this->bem_mesh_pt());

          // We (maybe) have to set the Newton solve tol lower because GMRES
          // tolerence isn't set high enough or because condition number bad?
          this->newton_solver_tolerance() = 1e-6;

          // Put multiple dofs into a single block depending on bools in Inputs.
          // ============================================================

          // We have to create vectors x where x[i] = the block that dof i goes
          // in.
          Vector<unsigned> dof_to_block_map;
          build_dof_to_block_map(Inputs::block_combine_m_directions,
                                 Inputs::block_combine_bulk_and_boundary_phi_1,
                                 Inputs::block_combine_bulk_and_boundary_phi,
                                 dof_to_block_map);

          // // Fill in the rest of the entries as needed to be in seperate blocks
          // for(unsigned j=0; j<dof_to_block_map.size(); j++)
          //   {
          //     if(dof_to_block_map[j] == unassigned_dof_magic_number)
          //       {
          //         dof_to_block_map[j] = dofs_so_far;
          //         dofs_so_far++;
          //       }
          //   }

          std::cout << "dof to block mapping is: " << dof_to_block_map << std::endl;

          // Put the mapping into the preconditioner.
          prec_pt->set_dof_to_block_map(dof_to_block_map);
        }

    }
  else
    {
      this->linear_solver_pt() = new SuperLUSolver;
    }

  if(Inputs::dump_bem)
    {
      this->boundary_matrix_pt()->Matrix<double, DenseMatrix<double> >::
        output("bem_matrix");
    }

} // end of constructor


  // =================================================================
  /// Construct the mapping between dofs as assigned in elements and the final
  /// blocks depending on which dofs we want lumped together.
  // =================================================================
template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
void ThreeDHybridProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
build_dof_to_block_map(const bool combine_m_directions,
                       const bool combine_bulk_and_boundary_phi_1,
                       const bool combine_bulk_and_boundary_phi,
                       Vector<unsigned>& dof_to_block_map) const
{
  // Create vector of a specific and unlikely number so that we can see if
  // the dof has been assigned a block yet.
  unsigned unassigned_dof_magic_number = 76543;
  dof_to_block_map.assign(7,unassigned_dof_magic_number);

  // Assign dofs to be in the same block in some specific cases
  unsigned dofs_so_far = 0;
  if(combine_m_directions)
    {
      dof_to_block_map[this->m_index(0)] = dofs_so_far;
      dof_to_block_map[this->m_index(1)] = dofs_so_far;
      dof_to_block_map[this->m_index(2)] = dofs_so_far;
      dofs_so_far++;
    }
  else
    {
      dof_to_block_map[this->m_index(0)] = dofs_so_far++;
      dof_to_block_map[this->m_index(1)] = dofs_so_far++;
      dof_to_block_map[this->m_index(2)] = dofs_so_far++;
    }

  if(combine_bulk_and_boundary_phi_1)
    {
      dof_to_block_map[this->phi_1_index()] = dofs_so_far;
      dof_to_block_map[5] = dofs_so_far;
      dofs_so_far++;
    }
  else
    {
      dof_to_block_map[this->phi_1_index()] = dofs_so_far++;
      dof_to_block_map[5] = dofs_so_far++;
    }

  if(combine_bulk_and_boundary_phi)
    {
      dof_to_block_map[this->phi_index()] = dofs_so_far;
      dof_to_block_map[6] = dofs_so_far;
      dofs_so_far++;
    }
  else
    {
      dof_to_block_map[this->phi_index()] = dofs_so_far++;
      dof_to_block_map[6] = dofs_so_far++;
    }

#ifdef PARANOID
  for(unsigned j=0; j<dof_to_block_map.size(); j++)
    {
      if(dof_to_block_map[j] == unassigned_dof_magic_number)
        {
          std::ostringstream error_msg;
          error_msg << "Dof number " << j << " was not assigned in the dof to block map.";
          throw OomphLibError(error_msg.str(),
                              "ThreeDHybridProblem::build_dof_to_block_map",
                              OOMPH_EXCEPTION_LOCATION);
        }
    }
#endif

  //??ds warning: if the dof numberings in MicromagneticsElement or the
  // boundary elements changes this will all be wrong! Also if you need
  // specific blocks be sure to account for this rearrangement!
}


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

namespace PreconditionerHelpers
{

  /// Generate a double in [0,1].
  double rand01()
  {
    unsigned n = 10000;
    return double( rand() % n) / double(n);
  }

  /// Fill a doublevector with random entries between zero and one.
  void randomise_vector(DoubleVector& vec)
  {
    unsigned n = vec. nrow();
    for(unsigned j=0; j<n; j++)
      vec[j] = rand01();
  }

  /// Do a superlu solve with random rhs and discard result (only useful for timings).
  double superlusolve(CRDoubleMatrix* matrix_pt)
  {
    LinearAlgebraDistribution* dist_pt = matrix_pt->distribution_pt();
    DoubleVector rhs(dist_pt), soln(dist_pt);
    randomise_vector(rhs);

    double start_time = TimingHelpers::timer();
    matrix_pt->solve(rhs,soln);

    double stop_time = TimingHelpers::timer();
    return stop_time - start_time;
  }

  double cgamgsolve(CRDoubleMatrix* matrix_pt)
  {
    LinearAlgebraDistribution* dist_pt = matrix_pt->distribution_pt();
    DoubleVector rhs(dist_pt), soln(dist_pt);
    randomise_vector(rhs);

    // Construct CG solver
    CG<CRDoubleMatrix> cg;

    double start_time = TimingHelpers::timer();

    // Solve
    HyprePreconditioner amg;
    amg.hypre_method() = HyprePreconditioner::BoomerAMG;
    cg.preconditioner_pt() = &amg;
    cg.solve(matrix_pt, rhs, soln);

    double stop_time = TimingHelpers::timer();
    return stop_time - start_time;
  }

  double cgilu0solve(CRDoubleMatrix* matrix_pt)
  {
    LinearAlgebraDistribution* dist_pt = matrix_pt->distribution_pt();
    DoubleVector rhs(dist_pt), soln(dist_pt);
    randomise_vector(rhs);

    // Construct CG solver
    CG<CRDoubleMatrix> cg;

    double start_time = TimingHelpers::timer();

    // Solve
    ILUZeroPreconditioner<CRDoubleMatrix> ilu0;
    cg.preconditioner_pt() = &ilu0;
    cg.solve(matrix_pt, rhs, soln);

    double stop_time = TimingHelpers::timer();
    return stop_time - start_time;
  }

  double gmresamgsolve(CRDoubleMatrix* matrix_pt)
  {
    LinearAlgebraDistribution* dist_pt = matrix_pt->distribution_pt();
    DoubleVector rhs(dist_pt), soln(dist_pt);
    randomise_vector(rhs);

    // Construct GMRES solver
    GMRES<CRDoubleMatrix> gmres;

    double start_time = TimingHelpers::timer();

    // Solve
    HyprePreconditioner amg;
    amg.hypre_method() = HyprePreconditioner::BoomerAMG;
    gmres.preconditioner_pt() = &amg;
    gmres.solve(matrix_pt, rhs, soln);

    double stop_time = TimingHelpers::timer();
    return stop_time - start_time;
  }

  double gmresilu0solve(CRDoubleMatrix* matrix_pt)
  {
    LinearAlgebraDistribution* dist_pt = matrix_pt->distribution_pt();
    DoubleVector rhs(dist_pt), soln(dist_pt);
    randomise_vector(rhs);

    // Construct GMRES solver
    GMRES<CRDoubleMatrix> gmres;

    double start_time = TimingHelpers::timer();

    // Solve
    ILUZeroPreconditioner<CRDoubleMatrix> ilu0;
    gmres.preconditioner_pt() = &ilu0;
    gmres.solve(matrix_pt, rhs, soln);

    double stop_time = TimingHelpers::timer();
    return stop_time - start_time;
  }

  double multiply(DoubleMatrixBase* dense_block_pt)
  {
    // Dummy dist with no communicator
    LinearAlgebraDistribution dist(0,dense_block_pt->nrow(),false);
    DoubleVector rhs(dist), soln(dist);
    randomise_vector(rhs);

    double start_time = TimingHelpers::timer();

    dense_block_pt->multiply(rhs,soln);

    double stop_time = TimingHelpers::timer();
    return stop_time - start_time;
  }

  void run_preconditioner_tests(Problem* problem_pt)
  {
    // Get a hybrid micromagnetics problem pt
    // This cast is a stupid hack but I don't care anymore... :(
    ThreeDHybridProblem< TMicromagElement <3,2>, MicromagFaceElement, 3 >*
      prob_pt =
      dynamic_cast<ThreeDHybridProblem< TMicromagElement <3,2>,
                                        MicromagFaceElement, 3 >*>
    (problem_pt);

  // Calculate the jacobian
  DoubleVector residuals;
  SumOfMatrices sumjacobian;
  prob_pt->get_jacobian(residuals, sumjacobian);

  // Get the sparse part of Jacobian and the dense sub block.
  CRDoubleMatrix* jacobian_pt
  = dynamic_cast<CRDoubleMatrix*>(sumjacobian.main_matrix_pt());
  DoubleMatrixBase* dense_block_pt = sumjacobian.added_matrix_pt(0);

  // Set up a dummy preconditioner to get blocks of the Jacobian
  DummyBlockPreconditioner<CRDoubleMatrix> blocker;
  blocker.add_mesh(prob_pt->bulk_mesh_pt());
  blocker.add_mesh(prob_pt->bem_mesh_pt());

  // Create a mapping to lump F block into one block and boundary/bulk phis into
  // one block.
  Vector<unsigned> dof_to_block_map;
  prob_pt->build_dof_to_block_map(1,1,1,dof_to_block_map);
  blocker.set_dof_to_block_map(dof_to_block_map);
  //??ds warning: if the dof numberings in MicromagneticsElement or the
  // boundary elements changes this will all be wrong!
  unsigned f_dof(0), phi1_dof(1), phi_dof(2);

  std::cout << dof_to_block_map << std::endl;

  // Finish setting up the blocking dummy preconditioner
  blocker.setup(problem_pt, jacobian_pt);

  // Get F, phi, phi1 blocks
  CRDoubleMatrix *f_block_pt(0), *phi1_block_pt(0), *phi_block_pt(0);
  blocker.get_block(f_dof, f_dof, f_block_pt);
  blocker.get_block(phi1_dof, phi1_dof, phi1_block_pt);
  blocker.get_block(phi_dof, phi_dof, phi_block_pt);

  // // Should fail if i got blocking right
  // CRDoubleMatrix *dummy(0);
  // blocker.get_block(4,4,dummy);

  std::cout << std::endl
  << "************************************************************"
  << std::endl;

  Vector<std::string> solve_types;
  solve_types.push_back("LU");
  solve_types.push_back("Krylov+amg");
  solve_types.push_back("Krylov+ilu0");

  Vector<double> f_times(solve_types.size());
  f_times[0] = superlusolve(f_block_pt);
  f_times[1] = gmresamgsolve(f_block_pt);
  f_times[2] =  -100; //gmresilu0solve(f_block_pt);

  Vector<double> phi_times(solve_types.size());
  phi_times[0] = superlusolve(phi_block_pt);
  phi_times[1] = cgamgsolve(phi_block_pt);
  phi_times[2] = cgilu0solve(phi_block_pt);

  Vector<double> phi1_times(solve_types.size());
  phi1_times[0] = superlusolve(phi1_block_pt);
  phi1_times[1] = cgamgsolve(phi1_block_pt);
  phi1_times[2] = cgilu0solve(phi1_block_pt);

  std::cout << "matrix " << solve_types << std::endl;
  std::cout << "f " << f_times << std::endl;
  std::cout << "phi " << phi_times << std::endl;
  std::cout << "phi1 " << phi1_times << std::endl;

  }

}



// takes x,y,z nelements, dt and output_step
int main(int argc, char* argv[])
{

  unsigned a, b, c, output_timestep;
  double dt;

  if(argc < 5)
    {
      std::cout << "Not enough command line inputs." << std::endl;
      // a = 30; b = 30; c = 100;
      a = Inputs::nx; b = Inputs::ny; c = Inputs::nz;
      dt = Inputs::dt;
      output_timestep = Inputs::output_timestep;
    }
  else
    {
      a = atoi(argv[1]);
      b = atoi(argv[2]);
      c = atoi(argv[3]);
      dt = atof(argv[4]);
      output_timestep = atoi(argv[5]);
    }

  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

#ifdef OOMPH_HAS_MPI
  // Initialise oomph-lib's MPI
  MPI_Helpers::init(argc,argv);
#endif

  // Inputs
  const unsigned dim = 3;
  const unsigned nnode_1d = 2;
  const double tmax = Inputs::tmax;

  // Dummy error for timestepper - always be ok
  const double dummy_t_eps = 100;


  // Create the problem
  ThreeDHybridProblem< TMicromagElement <dim,nnode_1d>, MicromagFaceElement, dim >
    problem(a,b,c);

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

  // /// Check problem
  // if(!(problem.self_test()==0))
  //   throw OomphLibError("Problem self_test failed","main",
  //                       OOMPH_EXCEPTION_LOCATION);

  std::cout << "constructor done, everything ready" << "\n" << std::endl;

  // Open a trace file
  std::ofstream trace_file;
  char trace_filename[100];
  sprintf(trace_filename,"%s/trace.dat",doc_info.directory().c_str());
  trace_file.open(trace_filename);

  // Pointer casting to get a block preconditioner pointer
  IterativeLinearSolver* it_lin_solver_pt =
    dynamic_cast<IterativeLinearSolver*>(problem.linear_solver_pt());
  GeneralPurposeBlockPreconditioner<CRDoubleMatrix>* prec_pt =
    dynamic_cast<GeneralPurposeBlockPreconditioner<CRDoubleMatrix>* >
    (it_lin_solver_pt->preconditioner_pt());

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

          // If this is the step where we wanted the blocked jacobian ouput then do it
          if(istep == output_timestep)
            {
              if(Inputs::solver_tests)
                {
                  PreconditionerHelpers::run_preconditioner_tests(&problem);
                  exit(0);
                }
              else
                prec_pt->set_block_output_to_files("blocks/jac");
            }
          else
            prec_pt->disable_block_output_to_files();

          // Take timestep
          problem.unsteady_newton_solve(dt);

          //Output solution
          problem.doc_solution(doc_info);

          //Increment counter for solutions
          doc_info.number()++;
          problem.debug_doc().next_timestep();

          // ??ds temp: Stop once output is done
          if(istep == output_timestep)
            exit(0);
        }

    }

#ifdef OOMPH_HAS_MPI
  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize();
#endif

}



#endif
