/*
  description of file goes here
*/

#include "generic.h"

#include "../../my_general_header.h"
#include "../../implicit_llg_problem.h"
#include "../../midpoint_method.h"
#include "../../magnetics_helpers.h"

// Floating point error checks
#include <fenv.h>

using namespace oomph;
using namespace MathematicalConstants;
using namespace StringConversion;

// // Calculations for the switching time and phi value for switching between
// // two angles in the very simple case when:
// // * H is aligned with the easy axis and constant
// // * There is no spatial variation
// // From Mallinson2000.
// namespace CompareSolutions
// {

//   double cart2theta(const Vector<double> &m)
//   {
//     double r = VectorOps::two_norm(m);
//     return std::acos(m[2]/r);
//   }

//   double initial_theta()
//   {
//     Vector<double> m,x;
//     Inputs::initial_m(0,x,m);
//     return cart2theta(m);
//   }

//   double cart2phi(const Vector<double> &m)
//   {
//     double result_in_mpi_pi = std::atan2(m[1],m[0]);
//     if (result_in_mpi_pi > 0)
//       return result_in_mpi_pi - 2*Pi;
//     else
//       return result_in_mpi_pi;
//   }

//   double switching_time(const double &alpha,
//                         const double &gamma,
//                         const double &H,
//                         const double &H_k,
//                         const double &theta_start,
//                         const double &theta_now)
//   {
//     using namespace std;
//     return ( (pow(alpha,2) + 1) / (alpha * gamma) )
//       * ( 1 / ( pow(H,2) - pow(H_k,2)))
//       * (H * log((tan(theta_now/2))
//                  / (tan(theta_start/2)))

//          + H_k * log((H - H_k * cos(theta_start))
//                      / (H - H_k * cos(theta_now)))

//          + H_k * log( sin(theta_now)
//                       / sin(theta_start)));
//   }

//   double switching_time_wrapper(const MagneticParameters* const parameters_pt,
//                                 const Vector<double> &m)
//   {
//     Vector<double> H, x, m_start_cartesian;
//     Inputs::applied_field(0,x,H); // assuming constant field!
//     Inputs::initial_m(0,x,m_start_cartesian); // assuming constant wrt x!

//     double theta_start = cart2theta(m_start_cartesian);
//     double theta_now = cart2theta(m);

//     double analytical_time =
//       switching_time(parameters_pt->normalised_gilbert_damping(),
//                      parameters_pt->normalised_gamma(),
//                      std::abs(H[2]),
//                      parameters_pt->normalised_hk(),
//                      theta_start,
//                      theta_now);

//     return analytical_time;
//   }

//   double switching_time_error(const MagneticParameters* const parameters_pt,
//                               const double &time,
//                               const Vector<double> &m)
//   {
//     double analytical_time = switching_time_wrapper(parameters_pt,m);
//     double result = std::abs(analytical_time - time);
//     return result;
//   }

//   double analytic_phi(const double &alpha,
//                       const double &theta_start,
//                       const double &theta_now)
//   {
//     double phi = (-1/alpha) * std::log((std::tan(theta_now/2))
//                                        /(std::tan(theta_start/2)));
//     return std::fmod(phi,2*Pi); // map into range 0:2pi
//   }

//   double analytic_phi_wrapper(const MagneticParameters* const parameters_pt,
//                               const Vector<double> &m)
//   {
//     Vector<double> x, m_start_cartesian;
//     Inputs::initial_m(0,x,m_start_cartesian); // assuming constant wrt x!

//     double theta_start = cart2theta(m_start_cartesian);
//     double theta_now = cart2theta(m);

//     return analytic_phi(parameters_pt->normalised_gilbert_damping(),
//                         theta_start,
//                         theta_now);
//   }

//   double phi_error(const MagneticParameters* const parameters_pt,
//                    const Vector<double> &m)
//   {
//     double phi_now = cart2phi(m);
//     double analytical_phi = analytic_phi_wrapper(parameters_pt,m);
//     double result = std::abs(analytical_phi - phi_now);
//     return result;
//   }

// }


int main(int argc, char *argv[])
{
  // Start MPI
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::init(argc,argv);
#endif

  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  // Read and process arguments
  const unsigned dim = 2, nnode1d = 2;
  MyCliArgs<QMicromagElement<dim, nnode1d> > args(argc, argv);

  // Create problem, get timestepper and mesh from input arguments.
  ImplicitLLGProblem<QMicromagElement<dim, nnode1d> > problem;
  problem.add_time_stepper_pt(args.time_stepper_pt);
  problem.bulk_mesh_pt() = args.mesh_pt;
  problem.bulk_mesh_pt()->setup_boundary_element_info();

  // Set magnetic parameters
  double gilbert_damping = 0.5, hk = 0.0; // (normalised hk)
  problem.mag_parameters_pt()->set_simple_llg_parameters();
  problem.mag_parameters_pt()->gilbert_damping() = gilbert_damping;
  problem.mag_parameters_pt()->exchange_constant() = 0.0;
  problem.mag_parameters_pt()->magnetostatic_debug_coeff() = 0.0;
  problem.mag_parameters_pt()->k1() = hk * mag_parameters::mu0/2;

  // Set applied field
  problem.applied_field_fct_pt() = args.h_app_fct_pt;

  // Finished setup, now we can build the problem
  problem.build();

  // Initialise problem
  problem.initialise_dt(args.dt);
  problem.set_initial_condition(args.initial_m_fct_pt);

  // Set up outputs
  DocInfo doc_info(args.outdir);
  problem.doc_solution(doc_info);

  // All ready: step until completion
  double dt = args.dt;
  while(problem.time() < args.tmax)
    {
      std::cout << "step number = " << doc_info.number()
                << ", time = " << problem.time()
                << ", dt = " << dt
                << ", |m| error = " << 1 - problem.mean_nodal_magnetisation_length()
                << std::endl;

      // The Newton step itself, adaptive if requested
      if(args.adaptive_flag())
        {
          dt = problem.adaptive_unsteady_newton_solve(dt, args.tol);
        }
      else
        {
          problem.unsteady_newton_solve(dt);
        }

      // Output
      problem.doc_solution(doc_info);
    }


#ifdef OOMPH_HAS_MPI
  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize();
#endif

  return 0;
}
