
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
#include "../../src/generic/general_purpose_block_preconditioners.h"
#include "../../src/generic/hypre_solver.h"
#include "../../src/generic/general_purpose_preconditioners.h"



#include "./my_general_header.h"
#include "./git_version.h"

#include <ctime>
#include <ostream>
#include <string>

namespace oomph
{

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


  class MyDocInfo : public DocInfo
  {
  public:
    /// Default constructor
    MyDocInfo() : DocInfo(), output_jacobian("never"), Args_pt(0) {}

    MyDocInfo(const std::string& directory,
              const MyCliArgs* const args_pt,
              const std::string& output_jacobian="never")
      : DocInfo(directory), output_jacobian(output_jacobian), Args_pt(args_pt)
    {}

    std::string output_jacobian;
    const MyCliArgs* Args_pt;
  };


  class MyProblem : public Problem
  {
  public:
    /// Default constructor
    MyProblem() :
      Problem(),
      Doc_info("results", 0),
      Use_time_adaptive_newton(false),
      Trace_filename("trace"),
      Info_filename("info"),
      Trace_seperator("; "),
      Dummy_doc_data(-1)
    {
      Dim = 0;
      Output_precision = 8;
    }

    /// Destructor
    virtual ~MyProblem() {}


    double smart_newton_solve(double dt, const double& tol)
    {
      // The Newton step itself, adaptive if requested.
      if(use_explicit())
        {
          explicit_timestep(dt);
        }
      else if(Use_time_adaptive_newton)
        {
          dt = adaptive_unsteady_newton_solve(dt, tol);
        }
      else
        {
          unsteady_newton_solve(dt);
        }

      return dt;
    }

    bool use_explicit()
      {
        return (explicit_time_stepper_pt() != 0)
          && (time_stepper_pt()->is_steady());
      }

    virtual void actions_before_newton_step()
      {
        // Output Jacobian and residuals if requested
        if((to_lower(Doc_info.output_jacobian) == "always")
           && (should_doc_this_step(time_pt()->dt(), time())))
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

    virtual void actions_before_explicit_timestep()
      {
        // Output Jacobian and residuals if requested
        if((to_lower(Doc_info.output_jacobian) == "always")
           && (should_doc_this_step(time_pt()->dt(), time())))
          {
            // label with doc_info number and the newton step number
            std::string label = to_string(Doc_info.number())
              + "_"
              + to_string(nnewton_step_this_solve() + 1);

            dump_current_mm_or_jacobian_residuals(label);
          }
      }


    virtual void actions_after_explicit_timestep()
    {
      Jacobian_setup_times.push_back
        (this->explicit_solver_pt()->jacobian_setup_time());
      Solver_times.push_back
        (this->explicit_solver_pt()->linear_solver_solution_time());

      const IterativeLinearSolver* its_pt
        = dynamic_cast<const IterativeLinearSolver*>(this->explicit_solver_pt());
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

    /// \short Write some general data about the previous time step to a
    /// trace file. Extend by overloading write_additional_trace_data(...).
    void write_trace()
    {
      std::ofstream trace_file((Doc_info.directory() + "/" + Trace_filename).c_str(),
                               std::ios::app);
      trace_file.precision(Output_precision);

      // Write out data that can be done for every problem
      trace_file
        << Doc_info.number()
        << Trace_seperator << time()
        << Trace_seperator << time_pt()->dt()
        << Trace_seperator << get_error_norm()

        << Trace_seperator << Nnewton_iter_taken
        << Trace_seperator << Solver_iterations

        << Trace_seperator << Solver_times
        << Trace_seperator << Jacobian_setup_times
        << Trace_seperator << Preconditioner_setup_times

        << Trace_seperator << lte_norm()
        << Trace_seperator << trace_value()

        << Trace_seperator << std::time(0)

        // Reserved slots in case I think of more things to add later
        << Trace_seperator << Dummy_doc_data
        << Trace_seperator << Dummy_doc_data
        << Trace_seperator << Dummy_doc_data
        << Trace_seperator << Dummy_doc_data
        << Trace_seperator << Dummy_doc_data
        << Trace_seperator << Dummy_doc_data
        << Trace_seperator << Dummy_doc_data
        << Trace_seperator << Dummy_doc_data
        << Trace_seperator << Dummy_doc_data;

      //??ds residuals?
      //??ds max values? Just list all values of vectors with ; as sep?

      // Add problem specific data
      write_additional_trace_data(trace_file);

      // Finish off this line
      trace_file << std::endl;
      trace_file.close();
    }


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
    virtual void doc_solution_additional(std::ofstream& soln_file) const {}
    virtual void final_doc_additional() const {}
    virtual void initial_doc_additional() const {}


    /// \short Outputs to be done at the start of a run (i.e. outputing
    /// basic info on command line args etc, writing trace file headers,
    /// outputting initial conditions).
    void initial_doc()
      {
#ifdef PARANOID
        if(*(Doc_info.directory().end()-1) == '/')
          {
            std::string error_msg = "Don't put a / on the end of results dir";
            throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif

        // Output Jacobian if requested
        if((to_lower(Doc_info.output_jacobian) == "at_start")
           || (to_lower(Doc_info.output_jacobian) == "always"))
          {
            dump_current_mm_or_jacobian_residuals("at_start");
          }

        // pvd file
        // ============================================================
        // Write start of .pvd XML file
        std::ofstream pvd_file((Doc_info.directory() + "/" + "soln.pvd").c_str(),
                               std::ios::out);
        pvd_file << "<?xml version=\"1.0\"?>" << std::endl
                 << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">"
                 << std::endl
                 << "<Collection>" << std::endl;
        pvd_file.close();


        // Trace file
        // ============================================================

        // Clear (by overwriting) and write headers
        std::ofstream trace_file((Doc_info.directory() + "/" + Trace_filename).c_str());
        trace_file
          << "DocInfo_numbers"
          << Trace_seperator << "times"
          << Trace_seperator << "dts"
          << Trace_seperator << "error_norms"

          << Trace_seperator << "n_newton_iters"
          << Trace_seperator << "n_solver_iters"

          << Trace_seperator << "solver_times"
          << Trace_seperator << "jacobian_setup_times"
          << Trace_seperator << "preconditioner_setup_times"

          << Trace_seperator << "LTE_norms"
          << Trace_seperator << "trace_values"

          << Trace_seperator << "unix_timestamp"

          // Reserved slots in case I think of more things to add later
          << Trace_seperator << "dummy"
          << Trace_seperator << "dummy"
          << Trace_seperator << "dummy"
          << Trace_seperator << "dummy"
          << Trace_seperator << "dummy"
          << Trace_seperator << "dummy"
          << Trace_seperator << "dummy"
          << Trace_seperator << "dummy"
          << Trace_seperator << "dummy";

        // Other headers depending on the derived class
        write_additional_trace_headers(trace_file);

        // Finish the line and close
        trace_file << std::endl;
        trace_file.close();


        // Info file
        // ============================================================
        std::ofstream info_file((Doc_info.directory() + "/" + Info_filename).c_str());
        info_file
          << "real_time " << real_date_time() << std::endl
          << "unix_time " << std::time(0) << std::endl
          << "git_version " << GitVersion::VERSION << std::endl
          << "driver_name " << CommandLineArgs::Argv[0] << std::endl;
        Doc_info.Args_pt->dump_args(info_file);
        info_file.close();


        // Write initial solution and anything else problem specific
        // (e.g. more trace file headers)
        this->doc_solution();
        initial_doc_additional();
      }

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

        // Write out anything requested from an inheriting class.
        final_doc_additional();
      }

    /// \short General output function: output to trace file. Maybe output
    /// Jacobian depending on Doc_info.output_jacobian. Maybe output full
    /// solution depending on should_doc_this_step(..) function. Extend by
    /// overloading doc_solution_additional(...).
    void doc_solution()
      {
        // Always output trace file data
        write_trace();

        // Output full set of data if requested for this timestep
        if(should_doc_this_step(time_pt()->dt(), time()))
          {
            std::ofstream soln_file((Doc_info.directory() + "/" + "soln" +
                                     to_string(Doc_info.number()) + ".dat").c_str(),
                                    std::ios::out);
            soln_file.precision(Output_precision);
            doc_solution_additional(soln_file);
            soln_file.close();

            // Write the simulation time and filename to the pvd file
            std::ofstream pvd_file((Doc_info.directory() + "/" + "soln.pvd").c_str(),
                                   std::ios::app);
            pvd_file.precision(Output_precision);
            pvd_file << "<DataSet timestep=\"" << time()
              << "\" group=\"\" part=\"0\" file=\"" << "soln"
                     << Doc_info.number() << ".vtu"
              << "\"/>" << std::endl;
            pvd_file.close();

            Doc_info.number()++;
          }
      }


    /// \short Dummy error norm calculator (overload in derived classes).
    virtual double get_error_norm() const
    {return Dummy_doc_data;}

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

    virtual double trace_value() const
    {
      unsigned nele = mesh_pt()->nelement();
      unsigned e = nele/2;
      if(dynamic_cast<FiniteElement*>(mesh_pt()->element_pt(e)))
        {
          // Just use an element somewhere in the middle...
          Node* trace_nd_pt = mesh_pt()->finite_element_pt(e)->node_pt(0);
          return trace_nd_pt->value(0);
        }
      else
        {
          // Not finite elements so no idea what to use
          return Dummy_doc_data;
        }
      }

    void dump_current_mm_or_jacobian_residuals(const std::string& label)
    {
      // We actually want the mass matrix if we are doing explicit steps
      if(use_explicit())
        {
          CRDoubleMatrix M;
          DoubleVector residuals;

          // This is how you get mass matrices, ugly...
          AssemblyHandler* old_assembly_handler_pt = this->assembly_handler_pt();
          ExplicitTimeStepHandler ehandler;
          this->assembly_handler_pt() = &ehandler;
          this->get_jacobian(residuals, M);
          this->assembly_handler_pt() = old_assembly_handler_pt;

          M.sparse_indexed_output(Doc_info.directory() + "/massmatrix_" + label,
                                  Output_precision, true);
          residuals.output(Doc_info.directory() + "/explicit_residual_" + label,
                           Output_precision);
        }
      else
        {
          CRDoubleMatrix J;
          DoubleVector residuals;
          this->get_jacobian(residuals, J);

          J.sparse_indexed_output(Doc_info.directory() + "/jacobian_" + label,
                                  Output_precision, true);
          residuals.output(Doc_info.directory() + "/residual_" + label,
                           Output_precision);
        }
    }


    IterativeLinearSolver* iterative_linear_solver_pt() const
    {
      return dynamic_cast<IterativeLinearSolver*>
        (this->linear_solver_pt());
    }


    /// \short Perform set up of problem.
    virtual void build(Vector<Mesh*>& bulk_mesh_pts)
      {

        // Push all the meshes into the problem's sub mesh list
        for(unsigned j=0; j<bulk_mesh_pts.size(); j++)
          {
            add_sub_mesh(bulk_mesh_pts[j]);
          }

        // If we have an iterative solver with a block preconditioner then
        // add all the meshes to the block preconditioner as well.
        IterativeLinearSolver* its_pt = iterative_linear_solver_pt();
        if(its_pt != 0)
          {
            BlockPreconditioner<CRDoubleMatrix>* bp_pt
              = dynamic_cast<BlockPreconditioner<CRDoubleMatrix>*>
              (its_pt->preconditioner_pt());

            if(bp_pt != 0)
              {
                bp_pt->set_nmesh(nsub_mesh());
                for(unsigned i=0; i< nsub_mesh(); i++)
                  {
                    bp_pt->set_mesh(i, mesh_pt(i));
                  }
              }
          }

        // Get the problem dimension
        FiniteElement* fele_pt = dynamic_cast<FiniteElement*>
          (bulk_mesh_pts[0]->element_pt(0));
        if(fele_pt != 0)
          {
            Dim = fele_pt->nodal_dimension();
          }
        else
          {
            // Presumably if a "bulk" mesh contains non-finite elements
            // then this is not a pde, so no dimension as such.
            Dim = 0;
          }

        // Set the solver for explicit timesteps (mass matrix) to CG with a
        // diagonal predconditioner.
        IterativeLinearSolver* mm_solver_pt = new CG<CRDoubleMatrix>;
        mm_solver_pt->preconditioner_pt() =
          new MatrixBasedLumpedPreconditioner<CRDoubleMatrix>;
        explicit_solver_pt() = mm_solver_pt;

      }

    /// \short Get problem dimension (nodal dimension).
    const unsigned dim() const {return this->Dim;}

    MyDocInfo Doc_info;
    unsigned Output_precision;

    bool Use_time_adaptive_newton;

    std::string Trace_filename;
    std::string Info_filename;

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
