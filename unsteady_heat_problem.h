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


// The unsteady heat equations
#include "../../src/unsteady_heat/unsteady_heat_elements.h"
#include "../../src/unsteady_heat/Tunsteady_heat_elements.h"
#include "../../src/unsteady_heat/unsteady_heat_flux_elements.h"



#include "my_generic_problem.h"
#include "my_general_header.h"
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

  Vector<double> exact(const double &t, const Vector<double> &x)
    {
      Vector<double> exact(1, 0.0);
      exact[0] = sin(k*x[0]) * cos(omega1 * t) * cos(omega2 * t) * exp(-beta *t);
      return exact;
    }

  void source(const double& t, const Vector<double> &x, double& u)
    {
      double a = sin(k*x[0]) * sin(omega1 * t) * cos(omega2 * t) * exp(-beta *t);
      double b = sin(k*x[0]) * cos(omega1 * t) * sin(omega2 * t) * exp(-beta *t);
      vector<double> u1 = exact(t,x);
      u = -1 * (u1[0] * (k*k - beta) - omega1 * a - omega2 * b);
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

  typedef TimeSpaceToDoubleVectFctPt UnsteadyExactSolutionFctPt;

  /// Constructor
  UnsteadyHeatProblem() : Source_fct_pt(0), Exact_solution_fct_pt(0) {}

  /// Destructor
  ~UnsteadyHeatProblem() {}

  /// \short Update the problem specs before next timestep:
  /// Set Dirchlet boundary conditions from exact solution.
  void actions_before_implicit_timestep();

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

  Vector<double> exact_solution(const double &t, const Vector<double> &x) const
    {
      #ifdef PARANOID
      if(Exact_solution_fct_pt == 0)
        {
      std::string err = "Exact_solution_fct_pt is null!";
      throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
        OOMPH_CURRENT_FUNCTION);
    }
#endif
      return Exact_solution_fct_pt(t, x);
    }

  /// Pointer to source function
  UnsteadyHeatEquationsBase::UnsteadyHeatSourceFctPt Source_fct_pt;

  UnsteadyExactSolutionFctPt Exact_solution_fct_pt;

  /// Pointer to control node at which the solution is documented
  Node* Control_node_pt;

  /// Function that does the real work of the constructors.
  void build(Vector<Mesh*>& bulk_mesh_pts);

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
          Vector<double> u = exact_solution(time, x);
          nod_pt->set_value(0, u[0]);
        }
    }
} // end of actions_before_implicit_timestep



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
      Vector<double> exact = exact_solution(time, x);

      error.push_back(std::abs(approx_values[0] - exact[0]));
    }

  return VectorOps::two_norm(error);
}


namespace UnsteadyHeatFactories
{
  inline UnsteadyHeatEquationsBase::UnsteadyHeatSourceFctPt
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

  inline UnsteadyHeatProblem::UnsteadyExactSolutionFctPt
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

  /// \short Make a mesh as specified by an input argument. Refined
  /// according to the given refinement level (in some way appropriate
  /// for that mesh type). Assumption: this will be passed into a
  /// problem, which will delete the pointer when it's done.
  Mesh* mesh_factory(const std::string& _mesh_name,
                     int refinement_level,
                     TimeStepper* time_stepper_pt,
                     double scaling_factor=1.0,
                     unsigned nnode1d = 2);
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
    using namespace UnsteadyHeatFactories;
    using namespace Factories;

    mesh_factory_pt = &mesh_factory;

    MyCliArgs::run_factories();

    source_fct_pt_name = to_lower(source_fct_pt_name);
    source_fct_pt = source_fct_pt_factory(source_fct_pt_name);
    exact_fct_pt = exact_fct_pt_factory(source_fct_pt_name);
    initial_condition_fpt = exact_fct_pt;
  }

  void assign_specific_parameters(MyProblem* problem_pt) const
  {
      UnsteadyHeatProblem* ust_pt =
        checked_dynamic_cast<UnsteadyHeatProblem*>(problem_pt);

      ust_pt->Exact_solution_fct_pt = exact_fct_pt;
      ust_pt->Source_fct_pt = source_fct_pt;
    }

  std::string source_fct_pt_name;

  UnsteadyHeatEquationsBase::UnsteadyHeatSourceFctPt source_fct_pt;
  UnsteadyHeatProblem::UnsteadyExactSolutionFctPt exact_fct_pt;
};


#endif
