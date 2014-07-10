#include "ode_problem.h"

namespace oomph
{

  namespace ODEFactories
  {

    SolutionFunctorBase* exact_solutions_factory(const std::string& exact_name)
    {
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
          return new deriv_functions::LLODESolution;
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
      else
        {
          throw OomphLibError("Unrecognised exact solution " + exact_name,
                              OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }

      return new SolutionFunctor(fpt, dfpt);
    }

  }
}
