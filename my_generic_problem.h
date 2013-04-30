#ifndef OOMPH_MY_GENERIC_PROBLEM_H
#define OOMPH_MY_GENERIC_PROBLEM_H

/*

  Possbily all output should be handled by a(nother) new class?
*/

#include "../../src/generic/problem.h"
#include "../../src/generic/oomph_utilities.h"

#include "./my_general_header.h"

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
      Problem(), Dim(0),
      Doc_info("results", 0),
      Trace_filename("trace"),
      Info_filename("info")
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
        << "stddev_jacobian_setup_time" << " "; // 10

      trace_file.close();

      // Output other useful info
      std::ofstream info_file((Doc_info.directory() + "/" + Info_filename).c_str());

      info_file
        << "real_time " << real_date_time() << std::endl
        << "unix_time " << time() << std::endl;

      Doc_info.Args_pt->dump_args(info_file);

      info_file.close();
    }

    /// \short Overload to write problem specific data into trace
    /// file. Note: don't add any line endings, just space seperated data.
    virtual void write_additional_trace_data(std::ofstream& trace_file) const {}

    /// \short Write a trace file
    void write_trace() const
    {
      std::ofstream trace_file((Doc_info.directory() + "/" + Trace_filename).c_str(),
                               std::ios::app);

      Vector<double> solver_iters(1, -1);
      Vector<double> solver_times(1, -1);
      Vector<double> jacobian_setup_times(1, -1);

      double dummy = -1;

      // Write out data that can be done for every problem
      trace_file
        << Doc_info.number() << " " // 0
        << time() << " " // 1
        << time_pt()->dt() << " " // 2
        << get_error_norm() << " " // 3

        << Nnewton_iter_taken << " " // 4

        << VectorOps::mean(solver_iters) << " " // 5
        << VectorOps::stddev(solver_iters) << " " // 6

        << VectorOps::mean(solver_times) << " " // 7
        << VectorOps::stddev(solver_times) << " " // 8
        << VectorOps::mean(jacobian_setup_times) << " " // 9
        << VectorOps::stddev(jacobian_setup_times) << " " // 10

        // Reserved slots in case I think of more things to add later
        << dummy << " " // 11
        << dummy << " " // 12
        << dummy << " " // 13
        << dummy << " " // 14
        << dummy << " " // 15
        << dummy << " " // 16
        << dummy << " " // 17
        << dummy << " " // 18
        << dummy << " " // 19
        << dummy << " "; // 20

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

        write_initial_data();
        this->doc_solution();
        initial_doc_additional();
      }

    void final_doc()
      {
        // Output Jacobian if requested
        if(to_lower(Doc_info.output_jacobian) == "at_end")
          {
            CRDoubleMatrix J;
            DoubleVector dummy;
            this->get_jacobian(dummy, J);
            J.sparse_indexed_output(Doc_info.directory() + "/" + "jacobian_at_end");
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
            CRDoubleMatrix J;
            DoubleVector dummy;
            this->get_jacobian(dummy, J);
            J.sparse_indexed_output(Doc_info.directory() + "/" + "jacobian" +
                                    to_string(Doc_info.number()));
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

    unsigned Dim;
    MyDocInfo Doc_info;

    std::string Trace_filename;
    std::string Info_filename;

  private:

    /// Inaccessible copy constructor
    MyProblem(const MyProblem &dummy)
    {BrokenCopy::broken_copy("MyProblem");}

    /// Inaccessible assignment operator
    void operator=(const MyProblem &dummy)
    {BrokenCopy::broken_assign("MyProblem");}
  };

} // End of oomph namespace

#endif
