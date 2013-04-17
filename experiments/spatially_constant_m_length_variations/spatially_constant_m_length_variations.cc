/*
  description of file goes here
*/

#include "generic.h"

#include "../../my_general_header.h"
#include "../../implicit_llg_problem.h"
#include "../../midpoint_method.h"

// Floating point error checks
#include <fenv.h>

// Mesh
#include "meshes/simple_rectangular_quadmesh.h"

using namespace oomph;
using namespace MathematicalConstants;
using namespace StringConversion;

namespace Inputs
{
  void initial_m(const double& t, const Vector<double> &x, Vector<double> &m)
  {
    m.assign(3,0.0);
    m[2] = 1;
    m[0] = 0.3; // start a little bit along +x direction (otherwise we
    // have undefined things...)
    VectorOps::normalise(m);
  }

  void varying_initial_m(const double& t, const Vector<double> &x, Vector<double> &m)
  {
    m.assign(3,0.0);

    double l = 5.0;
    m[0] = sin(x[0]*2*Pi/l) + sin(x[1]*2*Pi/l);
    m[1] = cos(x[0]*2*Pi/l) + cos(x[1]*2*Pi/l);
    m[2] = 1.0 - m[0] - m[1];

    // m[0] = x[0]/2 + x[1]/2;
    // m[1] = (1 - m[0]);
    // m[2] = 0.0;

    // m[0] = x[0];
    // m[1] = 0;
    // m[2] = 1 - x[0];

    VectorOps::normalise(m);
  }

  void applied_field(const double& t, const Vector<double> &x,
                     Vector<double> &h_app)
  {
    h_app.assign(3,0.0);
    h_app[2] = -3;
  }

}

// Calculations for the switching time and phi value for switching between
// two angles in the very simple case when:
// * H is aligned with the easy axis and constant
// * There is no spatial variation
// From Mallinson2000.
namespace CompareSolutions
{

  double cart2theta(const Vector<double> &m)
  {
    double r = VectorOps::two_norm(m);
    return std::acos(m[2]/r);
  }

  double initial_theta()
  {
    Vector<double> m,x;
    Inputs::initial_m(0,x,m);
    return cart2theta(m);
  }

  double cart2phi(const Vector<double> &m)
  {
    double result_in_mpi_pi = std::atan2(m[1],m[0]);
    if (result_in_mpi_pi > 0)
      return result_in_mpi_pi - 2*Pi;
    else
      return result_in_mpi_pi;
  }

  double switching_time(const double &alpha,
                        const double &gamma,
                        const double &H,
                        const double &H_k,
                        const double &theta_start,
                        const double &theta_now)
  {
    using namespace std;
    return ( (pow(alpha,2) + 1) / (alpha * gamma) )
      * ( 1 / ( pow(H,2) - pow(H_k,2)))
      * (H * log((tan(theta_now/2))
                 / (tan(theta_start/2)))

         + H_k * log((H - H_k * cos(theta_start))
                     / (H - H_k * cos(theta_now)))

         + H_k * log( sin(theta_now)
                      / sin(theta_start)));
  }

  double switching_time_wrapper(const MagneticParameters* const parameters_pt,
                                const Vector<double> &m)
  {
    Vector<double> H, x, m_start_cartesian;
    Inputs::applied_field(0,x,H); // assuming constant field!
    Inputs::initial_m(0,x,m_start_cartesian); // assuming constant wrt x!

    double theta_start = cart2theta(m_start_cartesian);
    double theta_now = cart2theta(m);

    double analytical_time =
      switching_time(parameters_pt->normalised_gilbert_damping(),
                     parameters_pt->normalised_gamma(),
                     std::abs(H[2]),
                     parameters_pt->normalised_hk(),
                     theta_start,
                     theta_now);

    return analytical_time;
  }

  double switching_time_error(const MagneticParameters* const parameters_pt,
                              const double &time,
                              const Vector<double> &m)
  {
    double analytical_time = switching_time_wrapper(parameters_pt,m);
    double result = std::abs(analytical_time - time);
    return result;
  }

  double analytic_phi(const double &alpha,
                      const double &theta_start,
                      const double &theta_now)
  {
    double phi = (-1/alpha) * std::log((std::tan(theta_now/2))
                                       /(std::tan(theta_start/2)));
    return std::fmod(phi,2*Pi); // map into range 0:2pi
  }

  double analytic_phi_wrapper(const MagneticParameters* const parameters_pt,
                              const Vector<double> &m)
  {
    Vector<double> x, m_start_cartesian;
    Inputs::initial_m(0,x,m_start_cartesian); // assuming constant wrt x!

    double theta_start = cart2theta(m_start_cartesian);
    double theta_now = cart2theta(m);

    return analytic_phi(parameters_pt->normalised_gilbert_damping(),
                        theta_start,
                        theta_now);
  }

  double phi_error(const MagneticParameters* const parameters_pt,
                   const Vector<double> &m)
  {
    double phi_now = cart2phi(m);
    double analytical_phi = analytic_phi_wrapper(parameters_pt,m);
    double result = std::abs(analytical_phi - phi_now);
    return result;
  }

}


int main(int argc, char *argv[])
{
  double t_start = TimingHelpers::timer();

  // Start MPI
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::init(argc,argv);
#endif

  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  // Parameter defaults
  double lx(1.0), ly(1.0);
  double gilbert_damping = 0.5, hk = 0.0; // (normalised hk)

  // Read args
  MyCliArgs args(argc, argv);

  // Create problem
  const unsigned dim = 2;
  ImplicitLLGProblem<QMicromagElement<dim,2> > problem;

  // Create timestepper
  if(to_lower(args.timestepper) == "midpoint") {
    MidpointMethod* ts_pt = new MidpointMethod(args.adaptive_flag(), 2);
    problem.add_time_stepper_pt(ts_pt);
    ts_pt->Fudge_factor = 0.1;
    }
  else if(to_lower(args.timestepper) == "bdf2") {
    problem.add_time_stepper_pt(new BDF<2>(args.adaptive_flag()));
    }
  else
    throw OomphLibError("Unrecognised time stepper",
                        OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);

  // Create mesh
  problem.bulk_mesh_pt() = new SimpleRectangularQuadMesh<QMicromagElement<dim,2> >
    (args.nx, args.ny, lx, ly, problem.time_stepper_pt());
  problem.bulk_mesh_pt()->setup_boundary_element_info();

  // Set magnetic parameters
  problem.mag_parameters_pt()->set_simple_llg_parameters();
  problem.mag_parameters_pt()->gilbert_damping() = gilbert_damping;
  problem.mag_parameters_pt()->exchange_constant() = 0.0;
  problem.mag_parameters_pt()->magnetostatic_debug_coeff() = 0.0;
  problem.mag_parameters_pt()->k1() = hk * mag_parameters::mu0/2;

  // Set applied field
  problem.applied_field_fct_pt() = &Inputs::applied_field;

  // Finished setup, now we can build the problem
  problem.build();

  // Initialise problem
  problem.initialise_dt(args.dt);
  // problem.set_initial_condition(Inputs::initial_m);

  problem.set_initial_condition(Inputs::varying_initial_m); //??ds
#warning "Not actually spatiall const..."

  // Set up outputs
  DocInfo doc_info(args.outdir);
  problem.doc_solution(doc_info);

  Vector<double> mean_m; problem.mean_magnetisation(mean_m);
  std::cout << CompareSolutions::cart2theta(mean_m) << " "
            << mean_m << std::endl;
  std::cout << CompareSolutions::cart2phi(mean_m) << std::endl;
  std::cout << CompareSolutions::switching_time_error(problem.mag_parameters_pt(),
                                                      problem.time(),
                                                      mean_m) << std::endl;
  std::cout << CompareSolutions::phi_error(problem.mag_parameters_pt(),
                                           mean_m) << std::endl;

  double mz = 1;
  Vector<Vector<double> > mean_m_vec;
  Vector<double> time_vec;

  // Timestep to convergence, for 1000 steps or until tmax.
  // while((std::abs(mz + 1) > 0.05) && (doc_info.number() < 5000)
  //       && (problem.time() < tmax))
  double dt = args.dt;
  while(problem.time() < args.tmax)
    {
      // Update average mz
      mz = problem.mean_mz();

      std::cout << "step number = " << doc_info.number()
                << ", time = " << problem.time()
                << ", dt = " << dt
                << ", mz = " << mz
                << ", |m| error = " << 1 - problem.mean_nodal_magnetisation_length()
                << std::endl;

      // Step, adaptive if requested
      if(args.adaptive_flag())
        { dt = problem.adaptive_unsteady_newton_solve(dt, args.tol); }
      else
        { problem.unsteady_newton_solve(dt); }

      // Store all the data we need for error checks later
      Vector<double> mean_m; problem.mean_magnetisation(mean_m);
      mean_m_vec.push_back(mean_m);
      time_vec.push_back(problem.time());

      // Output
      problem.doc_solution(doc_info);
      //conv_data.output_this_newton_step(outdir+"/convergence_data");
    }


  // Post-processing on mean_m etc.
  // ============================================================

  unsigned N = mean_m_vec.size();

  // Get phis, thetas, length_m
  Vector<double> phi_vec(N), theta_vec(N),
    ms(N);
  std::transform(mean_m_vec.begin(), mean_m_vec.end(),
                 phi_vec.begin(), CompareSolutions::cart2phi);
  std::transform(mean_m_vec.begin(), mean_m_vec.end(),
                 theta_vec.begin(), CompareSolutions::cart2theta);
  std::transform(mean_m_vec.begin(), mean_m_vec.end(),
                 ms.begin(), VectorOps::two_norm);

  // Get errors etc.
  Vector<double> st_err_vec(N), phi_err_vec(N), analytical_phi_vec(N),
    analytical_time_vec(N);
  for(unsigned i=0; i<N; i++)
    {
      st_err_vec[i] =
        CompareSolutions::switching_time_error(problem.mag_parameters_pt(),
                                               time_vec[i],
                                               mean_m_vec[i]);
      phi_err_vec[i] =
        CompareSolutions::phi_error(problem.mag_parameters_pt(),
                                    mean_m_vec[i]);
      analytical_phi_vec[i] =
        CompareSolutions::analytic_phi_wrapper(problem.mag_parameters_pt(),
                                               mean_m_vec[i]);

      analytical_time_vec[i] =
        CompareSolutions::switching_time_wrapper(problem.mag_parameters_pt(),
                                                 mean_m_vec[i]);
    }

  // Data dump!
  std::ofstream outfile((args.outdir + "/comparison_data").c_str());
  outfile << "t mx my mz r theta phi exact_t exact_phi error_t error_phi"
          << std::endl;
  for(unsigned i=0; i<N; i++)
    {
      outfile << std::setprecision(5)
              << time_vec[i] << " " // don't want high precision for time
                                    // because we want to search for t=1
              << std::setprecision(15)
              << mean_m_vec[i][0] << " "
              << mean_m_vec[i][1] << " "
              << mean_m_vec[i][2] << " "
              << ms[i] << " "
              << theta_vec[i] << " "
              << phi_vec[i] << " "
              << analytical_time_vec[i] << " "
              << analytical_phi_vec[i] << " "
              << st_err_vec[i] << " "
              << phi_err_vec[i] << " "
              << std::endl;
    }


#ifdef OOMPH_HAS_MPI
  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize();
#endif

  std::cout << TimingHelpers::timer() - t_start << std::endl;

  return 0;
}
