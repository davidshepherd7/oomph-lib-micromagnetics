

// Floating point debugging
#include <fenv.h>

#include "../../midpoint_method.h"
#include "generic.h"
#include "./unsteady_heat_example.h"



using namespace oomph;

namespace oomph
{

}

namespace input
{

}



int main(int argc, char *argv[])
{
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  CommandLineArgs::setup(argc,argv);

  //??ds
  MidpointMethod a(false), adaptive_midpoint(true);

  BDF<1> my_bdf1;
  BDF<2> my_bdf2;
  // BDF<2>* bdf2_pt = new BDF<2>;
  // BDF<2>* bdf2_pt1 = new BDF<2>;
  // BDF<2>* bdf2_pt2 = new BDF<2>;
  // BDF<2>* bdf2_pt3 = new BDF<2>;

  // Build problem
  UnsteadyHeatProblem<QUnsteadyHeatElement<2,3> >
    problem(&ExactSolnForUnsteadyHeat::get_source, &adaptive_midpoint);

  // Setup labels for output
  DocInfo doc_info;
  doc_info.set_directory("results");
  doc_info.number()=0;

  // Open a trace file
  ofstream trace_file;
  char filename[100];
  sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
  trace_file.open(filename);
  trace_file << "\"time\", \"u FE\", \"u exact\", \"norm of error\", \"norm of solution\", \"dt\""
             << std::endl;

  // Choose simulation interval and timestep
  double t_max = 1.0;
  double tol = 5e-3;
  double dt = 1e-6;

  // Initialise timestep -- also sets the weights for all timesteppers
  // in the problem.
  problem.initialise_dt(dt);
  problem.set_initial_condition();

  //Output initial condition
  problem.doc_solution(doc_info, trace_file);

  //Increment counter for solutions
  doc_info.number()++;

  // Timestepping loop
  while(problem.time() < t_max)
    {
      cout << "Timestep " << doc_info.number() << std::endl;
      std::cout << "dt = " << dt << std::endl;

      // Take timestep
      double dt_next = problem.adaptive_unsteady_newton_solve(dt, tol);
      dt = dt_next;

      //Output solution
      double error = problem.doc_solution(doc_info,trace_file);
      OOMPH_ASSERT(error < 5e-4);

      //Increment counter for solutions
      doc_info.number()++;
    }

  // Close trace file
  trace_file.close();

  return 0;
}; // end of main
