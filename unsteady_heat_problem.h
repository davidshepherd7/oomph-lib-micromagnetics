#ifndef OOMPH_UNSTEADY_HEAT_PROBLEM_H
#define OOMPH_UNSTEADY_HEAT_PROBLEM_H
//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//           Version 0.90. August 3, 2009.
//LIC//
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
//LIC//
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC//
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC//
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC//
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC//
//LIC//====================================================================
//Driver for 2D unsteady heat problem

//Generic routines
#include "generic.h"

// The unsteady heat equations
#include "unsteady_heat.h"

#include "my_generic_problem.h"
#include "vector_helpers.h"

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

namespace OscillatoryHeatEqn
{
  using namespace std;

  double k = 0.1;
  double omega1 = 0.3;
  double omega2 = 2;
  double beta = 0;

  double exact(const double &t, const Vector<double> &x)
    {
      return sin(k*x[0]) * cos(omega1 * t) * cos(omega2 * t) * exp(-beta *t);
    }

  void source(const double& t, const Vector<double> &x, double& u)
    {
      double a = sin(k*x[0]) * sin(omega1 * t) * cos(omega2 * t) * exp(-beta *t);
      double b = sin(k*x[0]) * cos(omega1 * t) * sin(omega2 * t) * exp(-beta *t);
      double u1 = exact(t,x);
      u = -1 * (u1 * (k*k - beta) - omega1 * a - omega2 * b);
    }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//=====start_of_problem_class=========================================
/// UnsteadyHeat problem
//====================================================================
class UnsteadyHeatProblem : public MyProblem
{
public:

  typedef double (*UnsteadyExactSolutionFctPt)(const double& time,
                                               const Vector<double>& x);

  /// Constructor
  UnsteadyHeatProblem() : Source_fct_pt(0), Exact_solution_fct_pt(0) {}

  /// Destructor
  ~UnsteadyHeatProblem() {}

  /// \short Update the problem specs before next timestep:
  /// Set Dirchlet boundary conditions from exact solution.
  void actions_before_implicit_timestep();

  /// \short Set initial condition (incl previous timesteps) according
  /// to specified function.
  void set_initial_condition(UnsteadyExactSolutionFctPt initial_soln_fct_pt);

  /// Doc the solution
  void doc_solution_additional(std::ofstream& some_file) const
  {mesh_pt()->output(some_file, 2);}

  /// Nothing extra (yet?)
  void write_additional_trace_data(std::ofstream& trace_file) const
  {}

  /// Error for adaptive timestepper (rms of nodal error determined by
  /// comparison with explicit timestepper result).
  double global_temporal_error_norm()
  {
    //Find out how many nodes there are in the problem
    unsigned n_node = mesh_pt()->nnode();

    // Get the error
    Vector<double> error(n_node, 0.0);
    for(unsigned i=0;i<n_node;i++)
      {
        // Get error in solution.
        error[i] = mesh_pt()->node_pt(i)->time_stepper_pt()->
          temporal_error_in_value(mesh_pt()->node_pt(i),0);
      }

    // Compute a norm
    double final_error = VectorOps::two_norm(error);

    return final_error;
  }

  double get_error_norm() const;

  /// Pointer to source function
  UnsteadyHeatEquationsBase::UnsteadyHeatSourceFctPt Source_fct_pt;

  UnsteadyExactSolutionFctPt Exact_solution_fct_pt;

  /// Pointer to control node at which the solution is documented
  Node* Control_node_pt;

void build()
  {

  // Choose a control node at which the solution is documented
  //----------------------------------------------------------
  // Total number of elements
  unsigned n_el=mesh_pt()->nelement();

  // Choose an element in the middle
  unsigned control_el=unsigned(n_el/2);

  // Choose its first node as the control node
  Control_node_pt=mesh_pt()->finite_element_pt(control_el)->node_pt(0);

  std::cout << "Recording trace of the solution at: "
            << Control_node_pt->x(0) << " "
            << Control_node_pt->x(1) << std::endl;


  // Set the boundary conditions for this problem:
  unsigned n_bound = mesh_pt()->nboundary();
  for(unsigned b=0;b<n_bound;b++)
    {
      unsigned n_node = mesh_pt()->nboundary_node(b);
      for (unsigned n=0;n<n_node;n++)
        {
          mesh_pt()->boundary_node_pt(b,n)->pin(0);
        }
    } // end of set boundary conditions


  // Complete the build of all elements so they are fully functional
  //----------------------------------------------------------------

  // Find number of elements in mesh
  unsigned n_element = mesh_pt()->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for(unsigned i=0;i<n_element;i++)
    {
      // Upcast from FiniteElement to the present element
      UnsteadyHeatEquationsBase *el_pt
        = dynamic_cast<UnsteadyHeatEquationsBase*>(mesh_pt()->element_pt(i));

      //Set the source function pointer
      el_pt->source_fct_pt() = Source_fct_pt;
    }

  // Do equation numbering
  std::cout <<"Number of equations: " << assign_eqn_numbers() << std::endl;

} // end of constructor

};

//=========start of actions_before_implicit_timestep===============================
/// \short Actions before timestep: update the domain, then reset the
/// boundary conditions for the current time.
//========================================================================
void UnsteadyHeatProblem::actions_before_implicit_timestep()
{
  // Get current time
  double time=time_pt()->time();

  //Loop over the boundaries
  unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned ibound=0;ibound<num_bound;ibound++)
    {
      // Loop over the nodes on boundary
      unsigned num_nod=mesh_pt()->nboundary_node(ibound);
      for (unsigned inod=0;inod<num_nod;inod++)
        {
          Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
          Vector<double> x(2);
          x[0]=nod_pt->x(0);
          x[1]=nod_pt->x(1);
          // Get current values of the boundary conditions from the
          // exact solution
          double u = Exact_solution_fct_pt(time, x);
          nod_pt->set_value(0, u);
        }
    }
} // end of actions_before_implicit_timestep



void UnsteadyHeatProblem::set_initial_condition
(UnsteadyExactSolutionFctPt initial_soln_fct_pt)
{
  // Backup time in global Time object
  double backed_up_time=time_pt()->time();

  // Past history needs to be established for t=time0-deltat, ...
  // Then provide current values (at t=time0) which will also form
  // the initial guess for the first solve at t=time0+deltat

  // Vector of exact solution value
  Vector<double> x(2);

  //Find number of nodes in mesh
  unsigned num_nod = mesh_pt()->nnode();

  // Set continuous times at previous timesteps:
  // How many previous timesteps does the timestepper use?
  int nprev_steps=time_stepper_pt()->nprev_values();

  Vector<double> prev_time(nprev_steps+1);
  for (int t=nprev_steps;t>=0;t--)
    {
      prev_time[t] = time_stepper_pt()->time_pt()->time(unsigned(t));
    }

  // Loop over current & previous timesteps
  for (int t=nprev_steps;t>=0;t--)
    {
      // Continuous time
      double time=prev_time[t];
      std::cout << "setting IC at time =" << time << std::endl;

      // Loop over the nodes to set initial guess everywhere
      for (unsigned n=0;n<num_nod;n++)
        {
          // Get nodal coordinates
          x[0]=mesh_pt()->node_pt(n)->x(0);
          x[1]=mesh_pt()->node_pt(n)->x(1);

          // Get exact solution at previous time
          double soln = initial_soln_fct_pt(time,x);

          // Assign solution
          mesh_pt()->node_pt(n)->set_value(t,0,soln);

          // Loop over coordinate directions: Mesh doesn't move, so
          // previous position = present position
          for (unsigned i=0;i<2;i++)
            {
              mesh_pt()->node_pt(n)->x(t,i)=x[i];
            }
        }
    }

  // Reset backed up time for global timestepper
  time_pt()->time()=backed_up_time;

} // end of set_initial_condition


double UnsteadyHeatProblem::get_error_norm() const
{
  double time = time_pt()->time();

  Vector<double> error; error.reserve(mesh_pt()->nnode());
  for(unsigned ind=0, nnd=mesh_pt()->nnode(); ind<nnd; ind++)
    {
      Node* nd_pt = mesh_pt()->node_pt(ind);

      Vector<double> approx_values(1,0.0), x(2,0.0);
      nd_pt->position(x);
      nd_pt->value(approx_values);
      double exact = Exact_solution_fct_pt(time, x);

      error.push_back(std::abs(approx_values[0] - exact));
    }

  return VectorOps::two_norm(error);
}


namespace UnsteadyHeatFactories
{
  UnsteadyHeatEquationsBase::UnsteadyHeatSourceFctPt
  source_fct_pt_factory(const std::string &source_fct_pt_name)
  {
    if(source_fct_pt_name == "oscillating")
      {
        return &OscillatoryHeatEqn::source;
      }
    else
      {
        std::string err("Unrecognised source function name ");
        err += source_fct_pt_name;
        throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
  }

  UnsteadyHeatProblem::UnsteadyExactSolutionFctPt
  exact_fct_pt_factory(const std::string &source_fct_pt_name)
  {
    if(source_fct_pt_name == "oscillating")
      {
        return &OscillatoryHeatEqn::exact;
      }
    else
      {
        std::string err("Unrecognised source function name ");
        err += source_fct_pt_name;
        throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
  }

}





class UnsteadyHeatArgs : public MyCliArgs
{
public:
  UnsteadyHeatArgs() : source_fct_pt(0), exact_fct_pt(0) {}

  virtual void set_flags()
  {
    MyCliArgs::set_flags();

    specify_command_line_flag("-source", &source_fct_pt_name);
    source_fct_pt_name = "oscillating";
  }

  virtual void run_factories()
  {
    MyCliArgs::run_factories();

    source_fct_pt_name = to_lower(source_fct_pt_name);
    source_fct_pt = UnsteadyHeatFactories::source_fct_pt_factory(source_fct_pt_name);
    exact_fct_pt = UnsteadyHeatFactories::exact_fct_pt_factory(source_fct_pt_name);
  }

  std::string source_fct_pt_name;

  UnsteadyHeatEquationsBase::UnsteadyHeatSourceFctPt source_fct_pt;
  UnsteadyHeatProblem::UnsteadyExactSolutionFctPt exact_fct_pt;
};

#endif
