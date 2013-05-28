
#ifndef OOMPH_MY_GENERIC_PROBLEM_H
#define OOMPH_MY_GENERIC_PROBLEM_H

/*

  Possbily all output should be handled by a(nother) new class?
*/

#include "../../src/generic/problem.h"
#include "../../src/generic/oomph_utilities.h"

#include "./my_general_header.h"
#include "./git_version.h"

#include <ctime>
#include <ostream>
#include <string>

using namespace oomph;

namespace oomph
{

  std::string real_date_time()
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
      Dim(0),
      Doc_info("results", 0),
      Trace_filename("trace"),
      Info_filename("info"),
      Trace_seperator("; "),
      Dummy_doc_data(-1)
    {}

    /// Destructor
    ~MyProblem() {}

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
            Preconditioner_setup_times.push_back(-1);
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

      //??ds rubbish fix..
      if(Solver_iterations.empty())
        {Solver_iterations.push_back(Dummy_doc_data);}
      if(Solver_times.empty())
        {Solver_times.push_back(Dummy_doc_data);}
      if(Jacobian_setup_times.empty())
        {Jacobian_setup_times.push_back(Dummy_doc_data);}
      if(Preconditioner_setup_times.empty())
        {Preconditioner_setup_times.push_back(Dummy_doc_data);}

      // Write out data that can be done for every problem
      trace_file
        << Doc_info.number() << Trace_seperator // 0
        << time() << Trace_seperator // 1
        << time_pt()->dt() << Trace_seperator // 2
        << get_error_norm() << Trace_seperator // 3

        << Nnewton_iter_taken << Trace_seperator // 4

        << VectorOps::mean(Solver_iterations) << Trace_seperator // 5
        << VectorOps::stddev(Solver_iterations) << Trace_seperator // 6

        << VectorOps::mean(Solver_times) << Trace_seperator // 7
        << VectorOps::stddev(Solver_times) << Trace_seperator // 8
        << VectorOps::mean(Jacobian_setup_times) << Trace_seperator // 9
        << VectorOps::stddev(Jacobian_setup_times) << Trace_seperator // 10
        << VectorOps::mean(Preconditioner_setup_times) << Trace_seperator // 11
        << VectorOps::stddev(Preconditioner_setup_times) << Trace_seperator // 12

        << lte_norm() << Trace_seperator // 13
        << trace_value() << Trace_seperator  // 14

        // Reserved slots in case I think of more things to add later
        << Dummy_doc_data << Trace_seperator // 15
        << Dummy_doc_data << Trace_seperator // 16
        << Dummy_doc_data << Trace_seperator // 17
        << Dummy_doc_data << Trace_seperator // 18
        << Dummy_doc_data << Trace_seperator // 19
        << Dummy_doc_data << Trace_seperator; // 20

      //??ds residuals?
      //??ds max values? Just list all values of vectors with ; as sep?

      // Add problem specific data
      write_additional_trace_data(trace_file);

      // Finish off this line
      trace_file << std::endl;
      trace_file.close();
    }


    /// \short Overload to write problem specific data into trace
    /// file. Note: don't add any line endings and seperate data with the
    /// Trace_seperator string.
    virtual void write_additional_trace_data(std::ofstream& trace_file) const {}


    /// ??ds
    virtual void doc_solution_additional(std::ofstream& soln_file) const = 0;
    virtual void final_doc_additional() const {};
    virtual void initial_doc_additional() const {};


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
            dump_current_jacobian("at_start");
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
          << "Doc_info number" << Trace_seperator // 0
          << "time" << Trace_seperator // 1
          << "dt" << Trace_seperator // 2
          << "error norm" << Trace_seperator // 3

          << "n newton iter" << Trace_seperator // 4
          << "mean solver iter" << Trace_seperator // 5
          << "std dev solver iter" << Trace_seperator // 6

          << "mean solver time" << Trace_seperator // 7
          << "std dev solver time" << Trace_seperator // 8
          << "mean jacobian setup time" << Trace_seperator // 9
          << "stddev jacobian setup time" << Trace_seperator // 10
          << "mean preconditioner setup time" << Trace_seperator // 11
          << "stddev preconditioner setup time" << Trace_seperator // 12

          << "LTE norm" << Trace_seperator // 13
          << "trace value" << Trace_seperator // 14

          // Reserved slots in case I think of more things to add later
          << "dummy" << Trace_seperator // 15
          << "dummy" << Trace_seperator // 16
          << "dummy" << Trace_seperator // 17
          << "dummy" << Trace_seperator // 18
          << "dummy" << Trace_seperator // 19
          << "dummy" << Trace_seperator // 20

          << std::endl;

        trace_file.close();


        // Info file
        // ============================================================
        std::ofstream info_file((Doc_info.directory() + "/" + Info_filename).c_str());
        info_file
          << "real_time " << real_date_time() << std::endl
          << "unix_time " << time() << std::endl
          << "git_version " << GitVersion::VERSION << std::endl
          << "build_time " << GitVersion::BUILD_TIME <<std::endl
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
            dump_current_jacobian("at_end");
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
            doc_solution_additional(soln_file);
            soln_file.close();

            // Write the simulation time and filename to the pvd file
            std::ofstream pvd_file((Doc_info.directory() + "/" + "soln.pvd").c_str(),
                                   std::ios::app);
            pvd_file << "<DataSet timestep=\"" << time()
              << "\" group=\"\" part=\"0\" file=\"" << "soln"
                     << Doc_info.number() << ".vtu"
              << "\"/>" << std::endl;
            pvd_file.close();

            // Output Jacobian if requested
            if(to_lower(Doc_info.output_jacobian) == "always")
              {
                dump_current_jacobian(to_string(Doc_info.number()));
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


    /// \short Calculate norm of timestepper's lte error estimate.
    double lte_norm() const
    {
      if(time_stepper_pt()->adaptive_flag())
        {
          double e = 0;
          for(unsigned nd=0, nnode=mesh_pt()->nnode(); nd<nnode; nd++)
            {
              Node* nd_pt = mesh_pt()->node_pt(nd);
              unsigned nvalue = nd_pt->nvalue();
              for(unsigned i=0; i<nvalue; i++)
                {
                  e += std::pow(nd_pt->time_stepper_pt()->temporal_error_in_value(nd_pt, i), 2);
                }
            }
          return std::sqrt(e);
        }
      else
        {
          return Dummy_doc_data;
        }
    }

    virtual double trace_value() const
      {
        // Just use an element somewhere in the middle...
        unsigned nele = mesh_pt()->nelement();
        Node* trace_nd_pt = mesh_pt()->finite_element_pt(nele/2)->node_pt(0);
        return trace_nd_pt->value(0);
      }

    void dump_current_jacobian(const std::string& label)
    {
      CRDoubleMatrix J;
      DoubleVector dummy;
      this->get_jacobian(dummy, J);
      J.sparse_indexed_output(Doc_info.directory() + "/" + "jacobian_" +
                              label);
    }

    unsigned Dim;
    MyDocInfo Doc_info;

    std::string Trace_filename;
    std::string Info_filename;

  protected:

    /// \short String to insert between fields in trace file. Use "; " by
    /// default (so that "," or " " can be used for lists if needed).
    std::string Trace_seperator;

  private:

    double Dummy_doc_data;
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
