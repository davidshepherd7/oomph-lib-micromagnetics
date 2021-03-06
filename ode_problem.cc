#include "ode_problem.h"
#include "magnetics_helpers.h"


namespace oomph
{

  namespace ODEFactories
  {

    SolutionFunctorBase* exact_solutions_factory(const std::string& exact_name,
                                                 const std::string& initial_m_name,
                                                 double sc_beta)
    {
      using namespace Factories;

      if(sc_beta != 0.0 && !(exact_name == "ll" || exact_name == "mallinson"))
        {
          std::string err = "Self correcting term only implemented for LL eqn.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      if(exact_name == "damped_oscillation")
        {
          return new deriv_functions::DampedOscillation;
        }
      else if(exact_name == "simple_stiff")
        {
          return new deriv_functions::SimpleStiffTest;
        }
      else if(exact_name == "order_reduction")
        {
          return new deriv_functions::OrderReductionTest;
        }
      else if(exact_name == "strong_order_reduction")
        {
          return new deriv_functions::OrderReductionTest;
        }
      else if(exact_name == "ll")
        {
          InitialM::LLODESolution* sol_pt = new InitialM::LLODESolution;
          sol_pt->initial_m_pt = initial_m_factory(initial_m_name);
          sol_pt->sc_beta = sc_beta;
          return sol_pt;
        }
      else if(exact_name == "mallinson")
        {
          InitialM::LLGMallinsonSolution* sol_pt = new InitialM::LLGMallinsonSolution;
          sol_pt->llg_ode_solution.initial_m_pt = initial_m_factory(initial_m_name);
	  sol_pt->llg_ode_solution.sc_beta = sc_beta;
          return sol_pt;
        }

      TimeSpaceToDoubleVectFctPt fpt;
      TimeSpaceValueToDoubleVectFctPt dfpt;

      if(exact_name == "sin")
        {
          fpt = &deriv_functions::sin;
          dfpt = &deriv_functions::dsin;
        }
      else if(exact_name == "cos")
        {
          fpt = &deriv_functions::cos;
          dfpt = &deriv_functions::dcos;
        }
      else if(exact_name == "exp")
        {
          fpt = &deriv_functions::exp;
          dfpt = &deriv_functions::dexp;
        }
      else if(exact_name == "poly3")
        {
          fpt = &deriv_functions::poly3;
          dfpt = &deriv_functions::dpoly3;
        }
      else if(exact_name == "stiff_test")
        {
          fpt = &deriv_functions::stiff_test;
          dfpt = &deriv_functions::dstiff_test;
        }
      else if(exact_name == "poly2")
        {
          fpt = &deriv_functions::poly2;
          dfpt = &deriv_functions::dpoly2;
        }
      else if(exact_name == "tanh")
        {
          fpt = &deriv_functions::tanh;
          dfpt = &deriv_functions::dtanh;
        }
      else
        {
          throw OomphLibError("Unrecognised exact solution " + exact_name,
                              OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }

      return new SolutionFunctor(fpt, dfpt);
    }

  }

  double LLGODEProblem::get_switching_time_error_norm(const unsigned& t_hist) const
  {
#ifdef PARANOID
    if(!Mallinson_applicable)
      {
        std::string err = "Can't calculate switching time";
        throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

    using namespace CompareSolutions;
    double exact_time = switching_time_wrapper(Magnetic_parameters_pt,
                                               exact_solution(0.0),
                                               solution());

    double actual_time = ts_pt()->time_pt()->time(t_hist);

    return std::abs(exact_time - actual_time);
  }

  double LLGODEProblem::get_error_norm(const unsigned& t_hist) const
  {
    // Use proper solution if we have it
    if(dynamic_cast<InitialM::LLGMallinsonSolution*>(Exact_solution_pt) != 0)
      {
        return ODEProblem::get_error_norm(t_hist);
      }

    // Otherwise maybe mallinson
    else if(Mallinson_applicable)
      {
        return get_switching_time_error_norm(t_hist);
      }

    // Otherwise no solution at all
    else
      {
        return MyProblem::Dummy_doc_data;
      }
  }
}
