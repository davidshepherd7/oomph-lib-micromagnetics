#ifndef OOMPH_GENERAL_INNER_ITERATION_PRECONDTIONER_H
#define OOMPH_GENERAL_INNER_ITERATION_PRECONDTIONER_H

#include "../../src/generic/preconditioner.h"
#include "../../src/generic/iterative_linear_solver.h"

namespace oomph
{

/// A preconditioner for performing inner iteration preconditioner solves.
/// Pass in a pointer to the solver to use as a preconditioner. This is a
/// replacement to InnerIterationPreconditioner which requires you to hard
/// code the desired solver and its precondioner and makes it very hard to
/// set any parameters for them.
class GeneralInnerIterationPreconditioner : public Preconditioner
{
public:

  GeneralInnerIterationPreconditioner()
  {
    Solver_pt = 0;
    Cached_matrix_pt = 0;
  }

  virtual ~GeneralInnerIterationPreconditioner() {}

  void clean_up_memory()
  {
    delete Cached_matrix_pt; Cached_matrix_pt = 0;
    Solver_pt->clean_up_memory();
  }

  /// \short Preconditioner setup method. Setup the preconditioner for the inner
  /// iteration solver.
  void setup()
  {
    oomph_info << "Setting up inner iteration preconditioner"
               << std::endl;

#ifdef PARANOID
    if(Solver_pt == 0)
      {
        std::string err = "No underlying solver has been set.";
        throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

    // Copy the matrix (matrix_pt() is out of scope after setup), assume
    // it's a cr matrix (not sure how to avoid that). ??ds ideally
    // ownership of the original matrix should somehow be passed to the
    // preconditioner, using a smart pointer or something...
    Cached_matrix_pt = new CRDoubleMatrix(*checked_dynamic_cast<CRDoubleMatrix*>(matrix_pt()));

    // set the distribution, not really sure why the solver can't do this
    // automatically
    DistributableLinearAlgebraObject* dist_pt =
      dynamic_cast<DistributableLinearAlgebraObject*>(Cached_matrix_pt);
    if(dist_pt != 0)
      {
        this->build_distribution(dist_pt->distribution_pt());
      }
    else
      {
        LinearAlgebraDistribution dist(comm_pt(), Cached_matrix_pt->nrow(), false);
        this->build_distribution(dist);
      }

    // Awful hack #1 from InnerIterationPreconditioner: manually set up the
    // preconditioner because when resolve is enabled the preconditioner is
    // NEVER set up, not even the first time solve is called!
    Solver_pt->preconditioner_pt()->setup(Cached_matrix_pt, comm_pt());


    // Enable resolve, i.e. store the matrix pointer and the preconditioner
    // for future use.
    Solver_pt->enable_resolve();

    // Awful hack #2 from InnerIterationPreconditioner: Store the matrix
    // pointer in the solver. This is done by "solving" the matrix with a
    // zero rhs and with only one iteration because there's no automatic
    // storage of the matrix when resolve is enabled.

    // backup max_iter
    unsigned max_iter = Solver_pt->max_iter();

    // "Solve" for one iteration
    Solver_pt->max_iter() = 1;
    DoubleVector x(this->distribution_pt(),0.0);
    DoubleVector y(x);
    Solver_pt->solve(Cached_matrix_pt, x, y);

    // revert max_iter
    Solver_pt->max_iter() = max_iter;
  }

  /// \short Preconditioner solve method. Performs the specified number
  /// of Krylov iterations preconditioned with the specified preconditioner
  void preconditioner_solve(const DoubleVector &r, DoubleVector &z)
  {
    Solver_pt->resolve(r,z);
  }

  /// Read/write access to the underlying solver
  IterativeLinearSolver*& solver_pt() {return Solver_pt;}


private:

  /// Pointer to the underlying solver
  IterativeLinearSolver* Solver_pt;

  // Pointer to the matrix for the solver to solve
  DoubleMatrixBase* Cached_matrix_pt;
};

} // End of oomph namespace

#endif
