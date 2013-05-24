
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
      Dummy_doc_data(-1)
    {}

    /// Destructor
    ~MyProblem() {}

    /// \short write out a general data file and initialise trace file
    void write_initial_data() const
    {
      // Clear the trace file (by overwriting) and write headers
      std::ofstream trace_file((Doc_info.directory() + "/" + Trace_filename).c_str());

      trace_file
        << "Doc_info.number()" << " " // 0
        << "time()" << " " // 1
        << "time_pt()->dt()" << " " // 2
        << "get_error_norm()" << " " // 3

        << "newton_iter" << " " // 4
        << "mean_solver_iter" << " " // 5
        << "std_dev_solver_iter" << " " // 6

        << "mean_solver_time" << " " // 7
        << "std_dev_solver_time" << " " // 8
        << "mean_jacobian_setup_time" << " " // 9
        << "stddev_jacobian_setup_time" << " " // 10
        << "mean_preconditioner_setup_time" << " " // 11
        << "stddev_preconditioner_setup_time" << " " // 12

        << "LTE_norm" << " " // 13
        << "trace_value" << " " // 14

        // Reserved slots in case I think of more things to add later
        << "dummy" << " " // 15
        << "dummy" << " " // 16
        << "dummy" << " " // 17
        << "dummy" << " " // 18
        << "dummy" << " " // 19
        << "dummy" << " " // 20

        << std::endl;

      trace_file.close();

      // Output other useful info
      std::ofstream info_file((Doc_info.directory() + "/" + Info_filename).c_str());

      info_file
        << "real_time " << real_date_time() << std::endl
        << "unix_time " << time() << std::endl
        << "git_version " << GitVersion::VERSION << std::endl
        << "build_time " << GitVersion::BUILD_TIME <<std::endl;

      Doc_info.Args_pt->dump_args(info_file);

      info_file.close();
    }

    /// \short Overload to write problem specific data into trace
    /// file. Note: don't add any line endings, just space seperated data.
    virtual void write_additional_trace_data(std::ofstream& trace_file) const {}

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

    /// \short Write a trace file
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
        << Doc_info.number() << " " // 0
        << time() << " " // 1
        << time_pt()->dt() << " " // 2
        << get_error_norm() << " " // 3

        << Nnewton_iter_taken << " " // 4

        << VectorOps::mean(Solver_iterations) << " " // 5
        << VectorOps::stddev(Solver_iterations) << " " // 6

        << VectorOps::mean(Solver_times) << " " // 7
        << VectorOps::stddev(Solver_times) << " " // 8
        << VectorOps::mean(Jacobian_setup_times) << " " // 9
        << VectorOps::stddev(Jacobian_setup_times) << " " // 10
        << VectorOps::mean(Preconditioner_setup_times) << " " // 11
        << VectorOps::stddev(Preconditioner_setup_times) << " " // 12

        << lte_norm() << " " // 13
        << trace_value() << " "  // 14

        // Reserved slots in case I think of more things to add later
        << Dummy_doc_data << " " // 15
        << Dummy_doc_data << " " // 16
        << Dummy_doc_data << " " // 17
        << Dummy_doc_data << " " // 18
        << Dummy_doc_data << " " // 19
        << Dummy_doc_data << " "; // 20

      //??ds residuals?
      //??ds max values? Just list all values of vectors with ; as sep?

      // Add problem specific data
      write_additional_trace_data(trace_file);

      // Finish off this line
      trace_file << std::endl;
      trace_file.close();
    }

    /// ??ds
    virtual void doc_solution_additional(std::ofstream& soln_file) const = 0;
    virtual void final_doc_additional() const {};
    virtual void initial_doc_additional() const {};

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

        write_initial_data();
        this->doc_solution();
        initial_doc_additional();
      }

    void final_doc()
      {
        // Output Jacobian if requested
        if((to_lower(Doc_info.output_jacobian) == "at_end")
           || (to_lower(Doc_info.output_jacobian) == "always"))
          {
            dump_current_jacobian("at_end");
          }

        final_doc_additional();
      }

    /// \short ??ds
    void doc_solution()
      {
        write_trace();

        // Output Jacobian if requested
        if(to_lower(Doc_info.output_jacobian) == "always")
          {
            dump_current_jacobian(to_string(Doc_info.number()));
          }

        std::ofstream soln_file((Doc_info.directory() + "/" + "soln" +
                                 to_string(Doc_info.number()) + ".dat").c_str(),
                                std::ios::out);
        doc_solution_additional(soln_file);
        soln_file.close();

        Doc_info.number()++;
      }

    /// \short Dummy error norm calculator (overload in derived classes).
    virtual double get_error_norm() const
    {
      return -1;
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

  private:

    double Dummy_doc_data;
    Vector<double> Jacobian_setup_times;
    Vector<double> Solver_times;
    Vector<double> Solver_iterations;
    Vector<double> Preconditioner_setup_times;

    /// Inaccessible copy constructor
    MyProblem(const MyProblem &dummy)
    {BrokenCopy::broken_copy("MyProblem");}

    /// Inaccessible assignment operator
    void operator=(const MyProblem &dummy)
    {BrokenCopy::broken_assign("MyProblem");}
  };

} // End of oomph namespace

#endif
