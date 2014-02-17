
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

      // Throw a real error (not just a warning) if the output directory
      // does not exist.
      Doc_info.enable_error_if_directory_does_not_exist();

      Disable_explicit_solver_optimisations = false;

      // By default output to trace file every step
      Always_write_trace = true;

      Dump = false;

      N_steps_taken = 0;
    }

    /// Destructor
    virtual ~MyProblem() {}

    /// Do the newton solve or explicit step (different ones depending flags
    /// set).
    double smart_time_step(double dt, const double& tol)
    {
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

      N_steps_taken++;

      return dt;
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

      // Not quite the same as actions after newton step because we are
      // interested in what happened in the explicit solver instead of the
      // main solver.
    }

    virtual void actions_before_explicit_timestep()
    {MyProblem::actions_before_newton_solve();}

    virtual void actions_after_explicit_timestep()
      {
        MyProblem::actions_after_newton_solve();
        check_error_norm_limits();
      }

    virtual void actions_after_implicit_timestep()
      {
        check_error_norm_limits();
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

    void check_error_norm_limits()
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
      }

    double min_element_size()
      {
        // Check that we have finite elements ??ds this will still go wrong
        // if there are only some none finite elements in the mesh...

        //??ds what happens with face elements?
        FiniteElement* test_pt = dynamic_cast<FiniteElement*>
          (mesh_pt(0)->element_pt(0));

        if(test_pt != 0)
          {
            double min_size = mesh_pt(0)->finite_element_pt(0)->size();

            // Loop over all meshes in problem
            for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
              {
                Mesh* mesh_pt = this->mesh_pt(msh);
                for(unsigned ele=0, nele=mesh_pt->nelement(); ele<nele; ele++)
                  {
                    FiniteElement* ele_pt = mesh_pt->finite_element_pt(ele);
                    double new_size = ele_pt->size();
                    if(new_size < min_size)
                      {
                        min_size = new_size;
                      }
                  }
              }

            return min_size;
          }
        // If it's not a finite element then we can't get a size so return
        // a dummy value.
        else
          {
            return Dummy_doc_data;
          }
      }

    /// \short Write some general data about the previous time step to a
    /// trace file. Extend by overloading write_additional_trace_data(...).
    void write_trace()
    {
      std::ofstream trace_file((Doc_info.directory() + "/" + Trace_filename).c_str(),
                               std::ios::app);
      trace_file.precision(Output_precision);

      double time = Dummy_doc_data, dt = Dummy_doc_data,
        lte_norm = Dummy_doc_data;
      if(!is_steady())
        {
          time = this->time();
          dt = this->time_pt()->dt();
          lte_norm = this->lte_norm();
        }

      // Write out data that can be done for every problem
      trace_file
        << Doc_info.number()
        << Trace_seperator << time
        << Trace_seperator << dt
        << Trace_seperator << get_error_norm()

        << Trace_seperator << Nnewton_iter_taken
        << Trace_seperator << Solver_iterations

        << Trace_seperator << Solver_times
        << Trace_seperator << Jacobian_setup_times
        << Trace_seperator << Preconditioner_setup_times

        << Trace_seperator << lte_norm
        << Trace_seperator << trace_value()

        << Trace_seperator << std::time(0)
        << Trace_seperator << min_element_size()

        // Reserved slots in case I think of more things to add later
        << Trace_seperator << Dummy_doc_data
        << Trace_seperator << Dummy_doc_data
        << Trace_seperator << Dummy_doc_data
        << Trace_seperator << Dummy_doc_data
        << Trace_seperator << Dummy_doc_data
        << Trace_seperator << Dummy_doc_data
        << Trace_seperator << Dummy_doc_data
        << Trace_seperator << Dummy_doc_data;

      //??ds residuals?

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
          << Trace_seperator << "min_element_size"

          // Reserved slots in case I think of more things to add later
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
          << "driver_name " << problem_name() << std::endl;
        info_file << Doc_info.args_str;
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
        bool doc_this_step = true;
        if(!is_steady())
          {
            doc_this_step = should_doc_this_step(time_pt()->dt(), time());
          }

        const std::string dir = Doc_info.directory();
        const std::string num = to_string(Doc_info.number());

        if(Always_write_trace || doc_this_step)
          {
            // Always output trace file data
            write_trace();
          }

        // Output full set of data if requested for this timestep
        if(doc_this_step)
          {

            // Solution itself
            std::ofstream soln_file((dir + "/" + "soln" + num + ".dat").c_str(),
                                    std::ios::out);
            soln_file.precision(Output_precision);
            doc_solution_additional(soln_file);
            soln_file.close();

            if(!is_steady())
              {
                // Write the simulation time and filename to the pvd file
                std::ofstream pvd_file((dir + "/" + "soln.pvd").c_str(),
                                       std::ios::app);
                pvd_file.precision(Output_precision);

                pvd_file << "<DataSet timestep=\"" << time()
                         << "\" group=\"\" part=\"0\" file=\"" << "soln"
                         << num << ".vtu"
                         << "\"/>" << std::endl;

                pvd_file.close();
              }


            // Maybe dump the restart data
            if(Dump)
              {
                std::ofstream dump_file((dir + "/" + "dump" + num + ".dat").c_str(),
                                        std::ios::out);
                this->dump(dump_file);
              }


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
    virtual void set_up_impulsive_initial_condition()
    {

#ifdef PARANOID
      if(nglobal_data() != 0)
        {
          std::string err = "Problem has global data which cannot be set from function pt.";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif
      unsigned nprev_steps=this->time_stepper_pt()->nprev_values();
      for(unsigned t=0; t< nprev_steps; t++)
        {
          // Loop over all nodes in all meshes in problem and set values.
          for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
            {
              Mesh* mesh_pt = this->mesh_pt(msh);

              for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
                {
                  Node* nd_pt = mesh_pt->node_pt(nd);
                  for(unsigned j=0, nj=nd_pt->nvalue(); j<nj; j++)
                    {
                      nd_pt->set_value(t, j, nd_pt->value(0, j));
                    }
                }

#ifdef PARANOID
              for(unsigned ele=0, nele=mesh_pt->nelement(); ele<nele; ele++)
                {
                  FiniteElement* ele_pt = mesh_pt->finite_element_pt(ele);
                  if(ele_pt->ninternal_data() != 0)
                    {
                      std::string err =
                        "Element with non-nodal data, cannot set via function...";
                      throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                          OOMPH_CURRENT_FUNCTION);
                    }
                }
#endif

            }
        }

      actions_after_set_initial_condition();
    }

    /// Assign initial conditions from function pointer
    virtual void set_initial_condition(InitialConditionFctPt ic_fpt)
      {
#ifdef PARANOID
        if(ic_fpt == 0)
          {
            std::string err = "Null inital condition function pointer";
            throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                OOMPH_CURRENT_FUNCTION);
          }

        if(nglobal_data() != 0)
          {
            std::string err = "Problem has global data which cannot be set from function pt.";
            throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                OOMPH_CURRENT_FUNCTION);
          }
#endif

        // Loop over current & previous timesteps
        int nprev_steps=this->time_stepper_pt()->nprev_values();
        for (int t=nprev_steps; t>=0; t--)
          {
            double time = time_pt()->time(t);

            // Loop over all nodes in all meshes in problem and set values.
            for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
              {
                Mesh* mesh_pt = this->mesh_pt(msh);

                for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
                  {
                    Node* nd_pt = mesh_pt->node_pt(nd);

                    // Get the position
                    const unsigned dim = nd_pt->ndim();
                    Vector<double> x(dim);
                    nd_pt->position(t, x);

                    // Get the values
                    Vector<double> values = ic_fpt(time, x);

                    // Copy into dofs
                    for(unsigned j=0, nj=values.size(); j<nj; j++)
                      {
                        nd_pt->set_value(t, j, values[j]);
                      }
                  }

#ifdef PARANOID
                for(unsigned ele=0, nele=mesh_pt->nelement(); ele<nele; ele++)
                  {
                    FiniteElement* ele_pt = mesh_pt->finite_element_pt(ele);
                    if(ele_pt->ninternal_data() != 0)
                      {
                        std::string err =
                          "Element with non-nodal data, cannot set via function...";
                        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                            OOMPH_CURRENT_FUNCTION);
                      }
                  }
#endif

                //??ds can't set external/internal data like this though
              }
          }

        actions_after_set_initial_condition();
      }

    /// Hook to be overloaded with any calculations needed after setting of
    /// initial conditions.
    virtual void actions_after_set_initial_condition() {}


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

    MyDocInfo Doc_info;
    unsigned Output_precision;

    unsigned N_steps_taken;

    std::string Trace_filename;
    std::string Info_filename;

    double Error_norm_limit;

    /// Option to turn off optimisation of the linear solves needed for
    /// explicit timestepping (for debugging purposes).
    bool Disable_explicit_solver_optimisations;

    // Should we output to trace file every step?
    bool Always_write_trace;

    /// Should we dump ready for a restart?
    bool Dump;

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
