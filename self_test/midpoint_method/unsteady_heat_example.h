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

// Mesh
#include "meshes/rectangular_quadmesh.h"


using namespace oomph;
using namespace MathematicalConstants;


#include "../../vector_helpers.h"

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//=====start_of_problem_class=========================================
/// UnsteadyHeat problem
//====================================================================
template<class ELEMENT>
class UnsteadyHeatProblem : public Problem
{

public:

  /// Constructor
  UnsteadyHeatProblem(
      UnsteadyHeatEquations<2>::UnsteadyHeatSourceFctPt source_fct_pt,
      FiniteElement::UnsteadyExactSolutionFctPt exact_solution_fct_pt,
      TimeStepper* ts_pt);

  /// Destructor (empty)
  ~UnsteadyHeatProblem(){}

  /// Update the problem specs after solve (empty)
  void actions_after_newton_solve() {}

  /// \short Update the problem specs before solve (empty)
  void actions_before_newton_solve() {}

  /// Update the problem specs after solve (empty)
  void actions_after_implicit_timestep() {}

  /// \short Update the problem specs before next timestep:
  /// Set Dirchlet boundary conditions from exact solution.
  void actions_before_implicit_timestep();

  /// \short Set initial condition (incl previous timesteps) according
  /// to specified function.
  void set_initial_condition();

  /// Doc the solution
  double doc_solution(DocInfo& doc_info, std::ofstream& trace_file);


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

private:

  /// Pointer to source function
  UnsteadyHeatEquations<2>::UnsteadyHeatSourceFctPt Source_fct_pt;

  FiniteElement::UnsteadyExactSolutionFctPt Exact_solution_fct_pt;

  /// Pointer to control node at which the solution is documented
  Node* Control_node_pt;

}; // end of problem class


//========start_of_constructor============================================
/// Constructor for UnsteadyHeat problem in square domain
//========================================================================
template<class ELEMENT> UnsteadyHeatProblem<ELEMENT>::
 UnsteadyHeatProblem(UnsteadyHeatEquations<2>::UnsteadyHeatSourceFctPt source_fct_pt,
                     FiniteElement::UnsteadyExactSolutionFctPt exact_solution_fct_pt,
                     TimeStepper* ts_pt) :
   Source_fct_pt(source_fct_pt), Exact_solution_fct_pt(exact_solution_fct_pt)
{


  // Allocate the timestepper -- this constructs the Problem's
  // time object with a sufficient amount of storage to store the
  // previous timsteps.
  add_time_stepper_pt(ts_pt);

  // Setup parameters for exact solution
  // -----------------------------------

  // Setup mesh
  //-----------

  // Number of elements in x and y directions
  unsigned nx=5;
  unsigned ny=5;

  // Lengths in x and y directions
  double lx=1.0;
  double ly=1.0;

  // Build mesh
  mesh_pt() = new RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly,time_stepper_pt());

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
  // ---------------------------------------------
  // All nodes are free by default -- just pin the ones that have
  // Dirichlet conditions here.
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
      ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

      //Set the source function pointer
      el_pt->source_fct_pt() = Source_fct_pt;
    }

  // Do equation numbering
  std::cout <<"Number of equations: " << assign_eqn_numbers() << std::endl;

} // end of constructor



//=========start of actions_before_implicit_timestep===============================
/// \short Actions before timestep: update the domain, then reset the
/// boundary conditions for the current time.
//========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::actions_before_implicit_timestep()
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
          Vector<double> u;
          Vector<double> x(2);
          x[0]=nod_pt->x(0);
          x[1]=nod_pt->x(1);
          // Get current values of the boundary conditions from the
          // exact solution
          Exact_solution_fct_pt(time,x,u);
          nod_pt->set_value(0,u[0]);
        }
    }
} // end of actions_before_implicit_timestep



//======================start_of_set_initial_condition====================
/// \short Set initial condition: Assign previous and current values
/// from exact solution.
//========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::set_initial_condition()
{
  // Backup time in global Time object
  double backed_up_time=time_pt()->time();

  // Past history needs to be established for t=time0-deltat, ...
  // Then provide current values (at t=time0) which will also form
  // the initial guess for the first solve at t=time0+deltat

  // Vector of exact solution value
  Vector<double> soln(1);
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
          Exact_solution_fct_pt(time,x,soln);

          // Assign solution
          mesh_pt()->node_pt(n)->set_value(t,0,soln[0]);

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



//=======start_of_doc_solution============================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
double UnsteadyHeatProblem<ELEMENT>::
doc_solution(DocInfo& doc_info, std::ofstream& trace_file)
{
  // ofstream some_file;
  // char filename[100];

  // // Number of plot points
  // unsigned npts;
  // npts=5;


  // cout << std::endl;
  // cout << "=================================================" << std::endl;
  // cout << "Docing solution for t=" << time_pt()->time() << std::endl;
  // cout << "=================================================" << std::endl;


  // // Output solution
  // //-----------------
  // sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
  //         doc_info.number());
  // some_file.open(filename);
  // mesh_pt()->output(some_file,npts);

  // // Write file as a tecplot text object
  // some_file << "TEXT X=2.5,Y=93.6,F=HELV,HU=POINT,C=BLUE,H=26,T=\"time = "
  //           << time_pt()->time() << "\"";
  // // ...and draw a horizontal line whose length is proportional
  // // to the elapsed time
  // some_file << "GEOMETRY X=2.5,Y=98,T=LINE,C=BLUE,LT=0.4" << std::endl;
  // some_file << "1" << std::endl;
  // some_file << "2" << std::endl;
  // some_file << " 0 0" << std::endl;
  // some_file << time_pt()->time()*20.0 << " 0" << std::endl;
  // some_file.close();


  // // Output exact solution
  // //----------------------
  // sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
  //         doc_info.number());
  // some_file.open(filename);
  // mesh_pt()->output_fct(some_file,npts,time_pt()->time(),
  //                       Exact_solution_fct_pt);
  // some_file.close();

  // // Doc error
  // //----------
  // double error,norm;
  // sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
  //         doc_info.number());
  // some_file.open(filename);
  // mesh_pt()->compute_error(some_file,
  //                          Exact_solution_fct_pt,
  //                          time_pt()->time(),
  //                          error,norm);
  // some_file.close();

  // // Doc solution and error
  // //-----------------------
  // cout << "error: " << error << std::endl;
  // cout << "norm : " << norm << std::endl << std::endl;

  // Get exact solution at control node
  Vector<double> x_ctrl(2);
  x_ctrl[0]=Control_node_pt->x(0);
  x_ctrl[1]=Control_node_pt->x(1);

  Vector<double> u_exact;
  Exact_solution_fct_pt(time_pt()->time(),x_ctrl,u_exact);

  trace_file << time_pt()->time() << " "
             << Control_node_pt->value(0) << " "
             << u_exact[0] << " "
             << get_error_norm() << " " // 3
             << time_pt()->dt() // 4
             << std::endl;

  // return error norm
  return 0.0;

} // end of doc_solution

template<class ELEMENT>
double UnsteadyHeatProblem<ELEMENT>::
get_error_norm() const
{
  double time = time_pt()->time();

  Vector<double> error; error.reserve(mesh_pt()->nnode());
  for(unsigned ind=0, nnd=mesh_pt()->nnode(); ind<nnd; ind++)
    {
      Node* nd_pt = mesh_pt()->node_pt(ind);

      Vector<double> approx_values(1,0.0), exact_values, x(2,0.0);
      nd_pt->position(x);
      nd_pt->value(approx_values);
      Exact_solution_fct_pt(time, x, exact_values);

      error.push_back(std::abs(approx_values[0] - exact_values[0]));
    }

  return VectorOps::two_norm(error);
}
