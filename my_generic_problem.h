
#ifndef OOMPH_MY_GENERIC_PROBLEM_H
#define OOMPH_MY_GENERIC_PROBLEM_H

/*

  Possbily all output should be handled by a(nother) new class?
*/

#include "../../src/generic/problem.h"
#include "../../src/generic/assembly_handler.h"
#include "../../src/generic/oomph_utilities.h"

#include "../../src/generic/linear_solver.h"
#include "../../src/generic/iterative_linear_solver.h"
#include "../../src/generic/matrices.h"
#include "../../src/generic/block_preconditioner.h"
#include "../../src/generic/general_purpose_block_preconditioners.h"
#include "../../src/generic/general_purpose_preconditioners.h"

#include "micromag_types.h"
#include "prettyprint98.hpp"
#include "energy_functions.h"


// #include "./my_general_header.h"

#include <ctime>
#include <ostream>
#include <string>

namespace oomph
{

using namespace StringConversion;

  class ElementalFunction;

  inline std::string real_date_time()
  {
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];

    // Get time
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    // Write into cstring
    strftime(buffer, 80, "%Y-%m-%d-%H-%M-%S", timeinfo);

    // Return as "real" string
    return std::string(buffer);
  }

  struct SolverParameters
  {
    SolverParameters()
    {
      linear_solver_pt = 0;
      mass_matrix_solver_pt = 0;
    }

    LinearSolver* linear_solver_pt;

    // Newton options
    double newton_solver_tolerance;
    unsigned max_newton_iterations;
    double max_residuals;
    bool shut_up_in_newton_solve;
    bool always_take_one_newton_step;

    // Linear optimisations
    bool jacobian_reuse_is_enabled;
    bool jacobian_has_been_computed;
    bool problem_is_nonlinear;

    // Explicit solver
    LinearSolver* mass_matrix_solver_pt;
    bool mass_matrix_reuse_is_enabled;
    bool mass_matrix_has_been_computed;
    bool discontinuous_element_formulation;
  };


  class MyDocInfo : public DocInfo
  {
  public:
    /// Default constructor
    MyDocInfo() : DocInfo(), output_jacobian("never") {}

    /// Copy dump of args into args_str.
    void copy_args_string()
      {
        std::ostringstream stream;
        CommandLineArgs::doc_all_flags(stream);
        args_str.assign(stream.str());
      }

    std::string output_jacobian;
    std::string args_str;
  };


  class MyProblem : public Problem
  {
  public:
    /// Default constructor
    MyProblem() :
      Trace_filename("trace"),
      Info_filename("info"),
      Trace_seperator("; "),
      Dummy_doc_data(-1)
    {
      Dim = 0;
      Output_precision = 8;
      Error_norm_limit = -1.0;
      Solution_norm_limit = -1.0;

      // Throw a real error (not just a warning) if the output directory
      // does not exist.
      Doc_info.enable_error_if_directory_does_not_exist();

      Disable_mass_matrix_solver_optimisations = false;

      // By default output to trace file every step
      Always_write_trace = true;

      Dump = false;
      Output_ltes = false;
      Output_predictor_values = false;
      Want_doc_exact = false;

      N_steps_taken = 0;
      Total_step_time= 0;
    }

    /// Virtual destructor. Policy decision: my problem classes won't call
    /// delete on anything. That's up to the driver code or
    /// whatever. Ideally we should just start using c++11 smart pointers
    /// and not have to worry so much about memory management!
    virtual ~MyProblem() {}

    /// Do the newton solve or explicit step (different ones depending flags
    /// set).
    double smart_time_step(double dt, const double& tol)
    {
      double step_time_start = TimingHelpers::timer();

      // The Newton step itself, adaptive if requested.
      if(explicit_flag())
        {
          explicit_timestep(dt);
        }
      else if(tol != 0.0)
        {
          dt = adaptive_unsteady_newton_solve(dt, tol);
        }
      else
        {
          unsteady_newton_solve(dt);
        }

      double step_time_stop = TimingHelpers::timer();
      Total_step_time = step_time_stop - step_time_start;


      oomph_info << "Time for step " << Total_step_time
                 << std::endl;

      N_steps_taken++;

      return dt;
    }

    /// Use to specify a condition for time stepping to halt early. By
    /// default never halt early.
    virtual bool finished() const
    {
      return false;
    }

    void get_solver_parameters(SolverParameters& sp)
    {
      sp.linear_solver_pt = linear_solver_pt();

      // Newton options
      sp.newton_solver_tolerance = newton_solver_tolerance();
      sp.max_newton_iterations = max_newton_iterations();
      sp.max_residuals = max_residuals();
      sp.shut_up_in_newton_solve = Shut_up_in_newton_solve;
      sp.always_take_one_newton_step = Always_take_one_newton_step;

      // Linear optimisations
      sp.jacobian_reuse_is_enabled = jacobian_reuse_is_enabled();
      sp.jacobian_has_been_computed = Jacobian_has_been_computed;
      sp.problem_is_nonlinear = Problem_is_nonlinear;

      // Explicit solver
      sp.mass_matrix_solver_pt = mass_matrix_solver_pt();
      sp.mass_matrix_reuse_is_enabled = mass_matrix_reuse_is_enabled();
      sp.mass_matrix_has_been_computed = Mass_matrix_has_been_computed;
      sp.discontinuous_element_formulation = Discontinuous_element_formulation;
    }

    void set_solver_parameters(SolverParameters& sp)
    {
      linear_solver_pt() = sp.linear_solver_pt;

      // Newton options
      newton_solver_tolerance() = sp.newton_solver_tolerance;
      max_newton_iterations() = sp.max_newton_iterations;
      max_residuals() = sp.max_residuals;
      Shut_up_in_newton_solve = sp.shut_up_in_newton_solve;
      Always_take_one_newton_step = sp.always_take_one_newton_step;

      // Linear optimisations
      Jacobian_reuse_is_enabled = sp.jacobian_reuse_is_enabled;
      Jacobian_has_been_computed = sp.jacobian_has_been_computed;
      Problem_is_nonlinear = sp.problem_is_nonlinear;

      // Explicit solver
      mass_matrix_solver_pt() = sp.mass_matrix_solver_pt;
      Mass_matrix_reuse_is_enabled = sp.mass_matrix_reuse_is_enabled;
      Mass_matrix_has_been_computed = sp.mass_matrix_has_been_computed;
      Discontinuous_element_formulation = sp.discontinuous_element_formulation;
    }

    bool explicit_flag()
      {
        return (explicit_time_stepper_pt() != 0)
          && (time_stepper_pt()->is_steady());
      }

    bool is_steady()
      {
        return (explicit_time_stepper_pt() == 0)
          && (time_stepper_pt()->is_steady());
      }

    virtual void actions_before_newton_step()
      {
        // Output Jacobian and residuals if requested
        if(to_lower(Doc_info.output_jacobian) == "always")
          {
            // label with doc_info number and the newton step number
            std::string label = to_string(Doc_info.number())
              + "_"
              + to_string(nnewton_step_this_solve() + 1);

            dump_current_mm_or_jacobian_residuals(label);
          }
      }

    unsigned nnewton_step_this_solve() const
      {
        return Jacobian_setup_times.size();
      }

    virtual void actions_before_explicit_stage()
      {MyProblem::actions_before_newton_step();}

    virtual void actions_after_explicit_stage()
    {
      Jacobian_setup_times.push_back
        (this->mass_matrix_solver_pt()->jacobian_setup_time());
      Solver_times.push_back
        (this->mass_matrix_solver_pt()->linear_solver_solution_time());

      // No non-linear residuals to store

      const IterativeLinearSolver* its_pt
        = dynamic_cast<const IterativeLinearSolver*>(this->mass_matrix_solver_pt());
      if(its_pt != 0)
        {
          Solver_iterations.push_back(its_pt->iterations());
          Preconditioner_setup_times.push_back(its_pt->preconditioner_setup_time());
        }
      else
        {
          // Fill in dummy data
          Solver_iterations.push_back(Dummy_doc_data);
          Preconditioner_setup_times.push_back(Dummy_doc_data);
        }

      // Not quite the same as actions after newton step because we are
      // interested in what happened in the explicit solver instead of the
      // main solver.
    }

    virtual void actions_before_time_integration() {}

    virtual void actions_before_explicit_timestep()
    {
      MyProblem::actions_before_newton_solve();
      check_norm_limits();
    }

    virtual void actions_after_explicit_timestep()
      {
        MyProblem::actions_after_newton_solve();
      }

    virtual void actions_after_implicit_timestep() {}

    virtual void actions_before_implicit_timestep()
    {
      check_norm_limits();
    }


    virtual void actions_after_newton_step()
      {
        Jacobian_setup_times.push_back
          (this->linear_solver_pt()->jacobian_setup_time());
        Solver_times.push_back
          (this->linear_solver_pt()->linear_solver_solution_time());

        const IterativeLinearSolver* its_pt
          = dynamic_cast<const IterativeLinearSolver*>(this->linear_solver_pt());
        if(its_pt != 0)
          {
            Solver_iterations.push_back(its_pt->iterations());
            Preconditioner_setup_times.push_back(its_pt->preconditioner_setup_time());
          }
        else
          {
            // Fill in dummy data
            Solver_iterations.push_back(Dummy_doc_data);
            Preconditioner_setup_times.push_back(Dummy_doc_data);
          }
      }

    virtual void actions_before_newton_solve()
      {
        // Clean up times vectors
        Jacobian_setup_times.clear();
        Solver_times.clear();
        Solver_iterations.clear();
        Preconditioner_setup_times.clear();
      }



    /// Pin all dofs with index in node one of indices in all nodes. Uses a
    /// different magic pinning number so that it can be easily undone
    /// using undo_segregated_pinning(). Does not handle: global data,
    /// element data, nodes with varying nvalue.
    void segregated_pin_indices(const Vector<unsigned>& indices)
    {
      for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
        {
          Mesh* mesh_pt = this->mesh_pt(msh);
          for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
            {
              Node* nd_pt = mesh_pt->node_pt(nd);
              for(unsigned j=0; j<indices.size(); j++)
                {
                  if(!nd_pt->is_pinned(indices[j]))
                    {
                      nd_pt->eqn_number(indices[j])
                        = Data::Is_segregated_solve_pinned;
                    }
                }
            }
        }
      oomph_info << "segregated solve, without indices: " << indices
                 << " n eqn: " << assign_eqn_numbers()  << std::endl;
    }

    /// Remove pinning set up by segregated_pin_indices.
    void undo_segregated_pinning()
    {
      for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
        {
          Mesh* mesh_pt = this->mesh_pt(msh);
          for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
            {
              Node* nd_pt = mesh_pt->node_pt(nd);
              for(unsigned j=0; j<nd_pt->nvalue(); j++)
                {
                  if(nd_pt->eqn_number(j) == Data::Is_segregated_solve_pinned)
                    {
                      nd_pt->eqn_number(j) = Data::Is_unclassified;
                    }
                }
            }
        }
      oomph_info << "un-segregated n eqn " << assign_eqn_numbers() << std::endl;
    }


    /// Check that nothing is currently pinned for a segregated solve.
    void check_not_segregated(const char* function) const
    {
#ifdef PARANOID
      for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
        {
          Mesh* mesh_pt = this->mesh_pt(msh);
          for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
            {
              Node* nd_pt = mesh_pt->node_pt(nd);
              for(unsigned j=0; j<nd_pt->nvalue(); j++)
                {
                  if(nd_pt->eqn_number(j) == Data::Is_segregated_solve_pinned)
                    {
                      throw OomphLibError("Some dofs already segregated",
                                          OOMPH_EXCEPTION_LOCATION,
                                          function);
                    }
                }
            }
        }
#endif
    }


    void check_norm_limits()
      {
        // If a limit has been set
        if(Error_norm_limit != -1.0)
          {
            double error_norm = get_error_norm();

            if((error_norm != Dummy_doc_data)
               && (error_norm > Error_norm_limit))
              {
                std::string err = "Error norm " + to_string(error_norm);
                err += " exceeds the limit " + to_string(Error_norm_limit);
                throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                    OOMPH_CURRENT_FUNCTION);
              }
          }

        // Same for solution norm
        if(Solution_norm_limit != -1.0)
          {
            double solution_norm = get_solution_norm();

            if((solution_norm != Dummy_doc_data)
               && (solution_norm > Solution_norm_limit))
              {
                std::string err = "Solution norm " + to_string(solution_norm);
                err += " exceeds the limit " + to_string(Solution_norm_limit);
                throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                    OOMPH_CURRENT_FUNCTION);
              }
          }
      }

    /// ??ds
    double min_element_size();

    /// \short Write some general data about the previous time step to a
    /// trace file. Extend by overloading write_additional_trace_data(...).
    void write_trace();


    /// \short Overload to write problem specific data into trace
    /// file. Note: don't add any line endings, seperate data with the
    /// Trace_seperator string (BEFORE each data entry).
    virtual void write_additional_trace_data(std::ofstream& trace_file) const {}

    /// \short Overload to write problem specific headers into trace
    /// file. Note: don't add any line endings, seperate headers with the
    /// Trace_seperator string (BEFORE each data entry).
    virtual void write_additional_trace_headers(std::ofstream& trace_file)
      const {}

    /// Overload to write any problem specific data
    virtual void output_solution(std::ofstream& soln_file) const {}
    virtual void final_doc_additional() const {}
    virtual void initial_doc_additional() const {}


    /// \short Outputs to be done at the start of a run (i.e. outputing
    /// basic info on command line args etc, writing trace file headers,
    /// outputting initial conditions).
    void initial_doc();

    /// \short Outputs to be done at the end of a run (e.g. closing tags
    /// for XML files).
    void final_doc()
      {
        // Output Jacobian if requested
        if((to_lower(Doc_info.output_jacobian) == "at_end")
           || (to_lower(Doc_info.output_jacobian) == "always"))
          {
            dump_current_mm_or_jacobian_residuals("at_end");
          }

        // Write end of .pvd XML file
        std::ofstream pvd_file((Doc_info.directory() + "/" + "soln.pvd").c_str(),
                               std::ios::app);
        pvd_file << "</Collection>" << std::endl
                 << "</VTKFile>" << std::endl;
        pvd_file.close();

        // Write end of exact.pvd XML file
        if(doc_exact())
          {
            std::ofstream pvd_file((Doc_info.directory() + "/" + "exact.pvd").c_str(),
                                   std::ios::app);
            pvd_file << "</Collection>" << std::endl
                     << "</VTKFile>" << std::endl;
            pvd_file.close();
          }

        // Write out anything requested from an inheriting class.
        final_doc_additional();
      }

    /// \short General output function: output to trace file. Maybe output
    /// Jacobian depending on Doc_info.output_jacobian. Maybe output full
    /// solution depending on should_doc_this_step(..) function. Maybe
    /// output ltes depending on value of Output_ltes. Extend by overloading
    /// output_solution(...).
    void doc_solution();


    /// Standard output function: loop over all elements in all meshes and
    /// output.
    virtual void output_solution(const unsigned& t, std::ostream& outstream,
      const unsigned& npoints=2) const
      {
        const unsigned n_msh = nsub_mesh();
        for(unsigned msh=0; msh<n_msh; msh++)
          {
            Mesh* msh_pt = mesh_pt(msh);

            const unsigned n_ele = msh_pt->nelement();
            for(unsigned ele=0; ele<n_ele; ele++)
              {
                FiniteElement* ele_pt = msh_pt->finite_element_pt(ele);
                ele_pt->output(t, outstream, npoints);
              }
          }
      }


    /// output_solution(...) with default output time step = 0 = current
    /// time.
    void output_solution(std::ostream& outstream,
                         const unsigned& npoints=2) const
    {
      output_solution(0, outstream, npoints);
    }

    /// Standard output function: loop over all elements in all meshes and
    /// output exact solution.
    virtual void output_exact_solution(std::ostream& outstream,
                                       const unsigned& npoints=2) const
    {
      const double time = time_pt()->time();

      const unsigned n_msh = nsub_mesh();
      for(unsigned msh=0; msh<n_msh; msh++)
        {
          Mesh* msh_pt = mesh_pt(msh);

          const unsigned n_ele = msh_pt->nelement();
          for(unsigned ele=0; ele<n_ele; ele++)
            {
              FiniteElement* ele_pt = msh_pt->finite_element_pt(ele);
              ele_pt->output_fct(outstream, npoints, time,
                                 *Exact_solution_pt);
            }
        }
    }

    /// \short Error norm calculator
    virtual double get_error_norm() const
    {
      if(Exact_solution_pt != 0)
        {
          // ExactFunctionDiffSquared f;
          // f.Exact_pt = Exact_solution_pt;
          // return std::sqrt(integrate_over_problem(&f));

          // Nodal rms difference
          const double t = time_pt()->time();

          double diffsq = 0;

          const unsigned n_node = mesh_pt()->nnode();
          for(unsigned nd=0; nd<n_node; nd++)
            {
              Node* nd_pt = mesh_pt()->node_pt(nd);
              Vector<double> values(nd_pt->nvalue(), 0.0), x(dim(), 0.0);
              nd_pt->position(x);
              nd_pt->value(values);

              Vector<double> exact = exact_solution(t, x);

              const unsigned ni = values.size();
              for(unsigned i=0; i<ni; i++)
                {
                  diffsq += std::pow(values[i] - exact[i], 2);
                }
            }

          return std::sqrt(diffsq);
        }
      else
        {
          return Dummy_doc_data;
        }
    }

    /// \short Dummy solution norm calculator (overload in derived classes).
    virtual double get_solution_norm() const
    {
      DoubleVector dofs;
      get_dofs(dofs);
      return dofs.norm();
    }

    /// \short should the previous step be doc'ed? Check if we went past an
    /// entry of Doc_times in the last step. If no Doc_times have been set
    /// then always output.
    virtual bool should_doc_this_step(const double &dt, const double &time) const
      {
        // I'm sure there should be a more efficient way to do this if we
        // know that Doc_times is ordered, but it shouldn't really matter I
        // guess--Jacobian calculation and solve times will always be far
        // far larger than this brute force search.

        // If no Doc_times have been set then always output.
        if(Doc_times.empty()) return true;

        // Loop over entries of Doc_times and check if they are in the
        // range (t - dt, t].
        for(unsigned j=0; j<Doc_times.size(); j++)
          {
            if(( time >= Doc_times[j]) && ((time - dt) < Doc_times[j]))
              {
                return true;
              }
          }
        return false;
      }

    /// \short Assign a vector of times to output the full solution at.
    void set_doc_times(Vector<double> &doc_times)
      {
        Doc_times = doc_times;
      }


    /// Get an lte error norm using the same norm and values as the
    /// adaptive time stepper used.
    double lte_norm()
    {
      if(time_stepper_pt()->adaptive_flag())
        {
          // Just return the error that would be used for adaptivity.
          return global_temporal_error_norm();
        }
      else
        {
          return Dummy_doc_data;
        }
    }

    bool doc_exact() const
    {
      return Want_doc_exact && (Exact_solution_pt != 0);
    }

    virtual Vector<double> trace_values() const
    {
      unsigned nele = mesh_pt()->nelement();
      unsigned e = nele/2;

      Vector<double> values;
      if(dynamic_cast<FiniteElement*>(mesh_pt()->element_pt(e)))
        {
          // Just use an element somewhere in the middle...
          Node* nd_pt = mesh_pt()->finite_element_pt(e)->node_pt(0);
          values.assign(nd_pt->nvalue(), 0.0);
          nd_pt->value(values);
        }
      else
        {
          // Not finite elements so no idea what to use
          values.assign(1, Dummy_doc_data);
        }

      return values;
    }


    void dump_current_mm_or_jacobian_residuals(const std::string& label);


    IterativeLinearSolver* iterative_linear_solver_pt() const
    {
      return dynamic_cast<IterativeLinearSolver*>
        (this->linear_solver_pt());
    }


    /// \short Perform set up of problem.
    virtual void build(Vector<Mesh*>& bulk_mesh_pts);

    /// \short Get problem dimension (nodal dimension).
    const unsigned dim() const {return this->Dim;}

    virtual std::string problem_name() const {return "unknown";}

    /// Set all history values/dts to be the same as the present values/dt.
    virtual void set_up_impulsive_initial_condition();

    /// Assign initial conditions from function pointer
    virtual void set_initial_condition(const InitialConditionFct& ic);

    /// Hook to be overloaded with any calculations needed after setting of
    /// initial conditions.
    virtual void actions_after_set_initial_condition();

    /// Integrate a function given by func_pt over every element in a mesh
    /// and return the total. This should probably be in the mesh class but
    /// that's core oomph-lib so I'll leave it here.
    virtual double integrate_over_mesh(const ElementalFunction* func_pt,
                                       const Mesh* const mesh_pt) const;

    /// \short Integrate a function given by func_pt over every element
    /// in every bulk mesh in this problem.
    virtual double integrate_over_problem(const ElementalFunction* func_pt) const;


    virtual void dump(std::ofstream& dump_file) const
      {
        // Set very high precision to avoid any issues
        dump_file.precision(14);

        dump_file << Doc_info.number() << " # Doc_info.number()" << std::endl;
        dump_file << N_steps_taken << " # N_steps_taken" << std::endl;
        Problem::dump(dump_file);
      }

    virtual void read(std::ifstream& restart_file)
    {
      // buffer
      std::string input_string;

      // Read in Doc_info number. Read line up to termination sign then
      // ignore.
      getline(restart_file, input_string, '#');
      restart_file.ignore(80,'\n');
      Doc_info.number() = std::atoi(input_string.c_str());

      // Read in number of steps taken. Read line up to termination sign
      // then ignore.
      getline(restart_file, input_string, '#');
      restart_file.ignore(80,'\n');
      N_steps_taken = std::atoi(input_string.c_str());

      // Let base class handle the rest
      Problem::read(restart_file);

      // Decrement doc info number so that it is correct after initial doc
      // Doc_info.number()--;
    }


    /// Output lte values of each nodal value along with spatial position
    /// for plotting with paraview. To plot load the csv file, use "Table
    /// To Points", add a 3D view, set the ltes to visible and color by the
    /// lte of interest. Useful to set "Representation" to "Points" and
    /// increase point size so that things are easily visible. Not
    /// implemented for nodes with varying numbers of values (you might
    /// need something fancier than a csv file for this).
    void output_ltes(std::ostream& output) const
      {
        // Output position labels
        for(unsigned j=0; j<Dim; j++)
          {
            output << "x" << j << ", ";
          }

        // Output labels for ltes, assuming that all nodes have the same
        // number of values...
        for(unsigned j=0; j<mesh_pt()->node_pt(0)->nvalue(); j++)
          {
            output << "error" << j << ", ";
          }

        output << std::endl;

        // Output actual positions and ltes
        for(unsigned i=0, ni=mesh_pt()->nnode(); i<ni; i++)
          {
            Node* nd_pt = mesh_pt()->node_pt(i);

            // Output position of node
            for(unsigned j=0; j<Dim; j++)
              {
                output << nd_pt->x(j) << ", ";
              }

            // Output ltes of node
            for(unsigned j=0; j<nd_pt->nvalue(); j++)
              {
                // Get timestepper's error estimate for this direction of m
                // at this point.
                double error = nd_pt->time_stepper_pt()->
                  temporal_error_in_value(nd_pt, j);

                output << error << ", ";
              }

            output << std::endl;
          }
      }


    MyDocInfo Doc_info;
    unsigned Output_precision;

    unsigned N_steps_taken;
    double Total_step_time;

    std::string Trace_filename;
    std::string Info_filename;

    double Error_norm_limit;
    double Solution_norm_limit;

    /// Option to turn off optimisation of the linear solves needed for
    /// explicit timestepping (for debugging purposes).
    bool Disable_mass_matrix_solver_optimisations;

    /// Should we output to trace file every step?
    bool Always_write_trace;

    /// Should we try to output exact solution?
    bool Want_doc_exact;

    /// Should we dump ready for a restart?
    bool Dump;

    /// Should we output the local truncation error at each node as well?
    bool Output_ltes;

    /// Should we output the predicted values too?
    bool Output_predictor_values;

    /// Function pointer for exact solution
    InitialConditionFct* Exact_solution_pt;

    /// Get exact solution
    Vector<double> exact_solution(const double& t, Vector<double>& x) const
    {
#ifdef PARANOID
      if(Exact_solution_pt == 0)
        {
          std::string err = "Exact_solution_pt is null!";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif
      return (*Exact_solution_pt)(t, x);
    }

  protected:

    /// \short String to insert between fields in trace file. Use "; " by
    /// default (so that "," or " " can be used for lists if needed).
    const std::string Trace_seperator;
    double Dummy_doc_data;
    unsigned Dim;

  private:

    Vector<double> Jacobian_setup_times;
    Vector<double> Solver_times;
    Vector<double> Solver_iterations;
    Vector<double> Preconditioner_setup_times;

    /// Times at which we want to output the full solution.
    Vector<double> Doc_times;

    /// Inaccessible copy constructor
    MyProblem(const MyProblem &dummy)
    {BrokenCopy::broken_copy("MyProblem");}

    /// Inaccessible assignment operator
    void operator=(const MyProblem &dummy)
    {BrokenCopy::broken_assign("MyProblem");}
  };

} // End of oomph namespace

#endif
