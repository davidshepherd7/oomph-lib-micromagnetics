#include "my_generic_problem.h"
#include "magnetics_helpers.h"
#include "energy_functions.h"


namespace oomph
{

  namespace MManipulation
  {

    /// \short Compute the effective damping constant (alpha) for the
    /// previous time step (see Albuquerque2001).
    double alt_effective_damping_used(MyProblem* problem_pt,
                                      std::deque<double>& previous_energies)
    {
      // Integral over all space of (dm/dt)^2 used in last step
      double dmdt_squared = integral_of_dmdt_squared(problem_pt); //??ds

      // If no change then damping is undefined
      if(dmdt_squared  == 0) return nan("");

      // Forumla from Albuquerque2001 & dAquino2005
      double dEdt = alt_dEnergydt(problem_pt, previous_energies);
      double effective_alpha = - dEdt / dmdt_squared;

      return effective_alpha;
    }


    /// \short Compute the effective damping constant (alpha) for the
    /// previous time step (see Albuquerque2001).
    double effective_damping_used(MyProblem* problem_pt)
    {
      // Integral over all space of (dm/dt)^2 used in last step
      double dmdt_squared = integral_of_dmdt_squared(problem_pt);

      // If no change then damping is undefined
      if(dmdt_squared  == 0) return nan("");

      // Forumla from Albuquerque2001 & dAquino2005
      double dEdt = dEnergydt(problem_pt);
      double effective_alpha = - dEdt / dmdt_squared;

      return effective_alpha;
    }


    double exchange_energy(MyProblem* problem_pt)
    {
      ExchangeEnergyFunction f;
      return problem_pt->integrate_over_problem(&f);
    }


    double zeeman_energy(MyProblem* problem_pt)
    {
      ZeemanEnergyFunction f;
      return problem_pt->integrate_over_problem(&f);
    }

    double crystalline_anisotropy_energy(MyProblem* problem_pt)
    {
      CrystallineAnisotropyEnergyFunction f;
      return problem_pt->integrate_over_problem(&f);
    }


    double magnetostatic_energy(MyProblem* problem_pt)
    {
      MagnetostaticEnergyFunction f;
      return problem_pt->integrate_over_problem(&f);
    }

    double integral_of_dmdt_squared(MyProblem* problem_pt)
    {
      DmdtSquaredFunction f;
      return problem_pt->integrate_over_problem(&f);
    }

    double dEnergydt(MyProblem* problem_pt)
    {
      dExchangeEnergydtFunction de_exdt;
      double I_de_exdt = problem_pt->integrate_over_problem(&de_exdt);

      dZeemanEnergydtFunction de_zeedt;
      double I_de_zeedt = problem_pt->integrate_over_problem(&de_zeedt);

      dCrystallineAnisotropydtEnergyFunction de_cadt;
      double I_de_cadt = problem_pt->integrate_over_problem(&de_cadt);

      dMagnetostaticEnergydtFunction de_ms;
      double I_de_ms = problem_pt->integrate_over_problem(&de_ms);

      return I_de_exdt + I_de_zeedt + I_de_cadt + I_de_ms;
    }

    double alt_dEnergydt(MyProblem* problem_pt,
std::deque<double>& previous_energies)
    {
      // Make a BDF2 time stepper to look up weights from (because I'm
      // lazy...)
      BDF<2> bdf;
      TimeStepper* node_ts_pt = problem_pt->mesh_pt()->
        finite_element_pt(0)->node_pt(0)->time_stepper_pt();
      bdf.time_pt() = node_ts_pt->time_pt();
      bdf.set_weights();

      // Calculate first derivative
      double deriv = 0.0;
      for(unsigned t=0;t<bdf.ntstorage();t++)
        {
          deriv += bdf.weight(1,t) * previous_energies[t];
        }

      return deriv;
    }
  }

}
