#ifndef OOMPH_MY_GENERIC_PROBLEM_H
#define OOMPH_MY_GENERIC_PROBLEM_H

/*

  Possbily all output should be handled by a(nother) new class?
*/

#include "../../src/generic/problem.h"
#include "../../src/generic/oomph_utilities.h"

#include <ctime>
#include <iostream>
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

  template<class ELEMENT>
  class MyProblem : public Problem
  {
  public:
    /// Default constructor
    MyProblem() :
      Problem(), Dim(0), Ele_pt(0),
      Doc_info("results/"),
      Trace_filename("trace"),
      Info_filename("info")
    {}

    /// Destructor
    ~MyProblem() {}

    /// \short write out a general data file and initialise trace file
    void write_initial_data() const
    {
      // Clear the trace file (just in case) and write headers
      std::ofstream trace_file((Doc_info.directory() + Trace_filename).c_str());
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


      std::ofstream info_file((Doc_info.directory() + Info_filename).c_str());

      info_file
        << "real time: " << real_date_time() << std::endl
        << "unix time: " << time() << std::endl
        << "mesh name: " << std::endl
        << "input arguments: " << std::endl;
    }

    /// \short Overload to write problem specific data into trace
    /// file. Note: don't add any line endings, just space seperated data.
    virtual void write_additional_trace_data(std::ofstream& trace_file) const {}

    /// \short Write a trace file
    void write_trace() const
    {
      std::ofstream trace_file((Doc_info.directory() + Trace_filename).c_str(),
                               std::ios::app);

      // Write out data that can be done for every problem
      trace_file
        << Doc_info.number() << " " // 0
        << time() << " " // 1
        << time_pt()->dt() << " " // 2
        << get_error_norm() << " " // 3

        << Nnewton_iter_taken << " "; // 4

        // << VectorOps::mean(solver_iters) << " " // 5
        // << VectorOps::stddev(solver_iters) << " " // 6

        // << VectorOps::mean(solver_times) << " " // 7
        // << VectorOps::std_dev(solver_times) << " " // 8
        // << VectorOps::mean(jacobian_setup_times) << " " // 9
        // << VectorOps::stddev(jacobian_setup_times) << " "; // 10


      // Add problem specific data
      write_additional_trace_data(trace_file);

      // Finish off this line
      trace_file << std::endl;
      trace_file.close();
    }

    /// ??ds
    virtual void doc_solution_additional() const = 0;
    virtual void final_doc_additional() const {};
    virtual void initial_doc_additional() const {};

    void initial_doc()
      {
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
            J.sparse_indexed_output(Doc_info.directory() + "jacobian_at_end");
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
            J.sparse_indexed_output(Doc_info.directory() + "jacobian" +
                                    to_string(Doc_info.number()));
          }

        doc_solution_additional();

        Doc_info.number()++;
      }

    /// \short Dummy error norm calculator (overload in derived classes).
    virtual double get_error_norm() const
    {
      return -1;
    }

    unsigned Dim;
    ELEMENT* Ele_pt;
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
