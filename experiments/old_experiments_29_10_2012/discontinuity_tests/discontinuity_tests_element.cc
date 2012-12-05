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
//Header file for UnsteadyHeat elements
#ifndef OOMPH_UNSTEADY_HEAT_ELEMENTS_HEADER
#define OOMPH_UNSTEADY_HEAT_ELEMENTS_HEADER
				       				       
// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


//OOMPH-LIB headers
#include "generic/nodes.h"
#include "generic/Qelements.h"
#include "generic/oomph_utilities.h"


namespace oomph
{

  //=============================================================
  /// A class for all isoparametric elements that solve the 
  /// UnsteadyHeat equations.
  /// \f[ 
  /// \frac{\partial^2 u}{\partial x_i^2}=\frac{\partial u}{\partial t}+f(t,x_j)
  /// \f] 
  /// This contains the generic maths. Shape functions, geometric
  /// mapping etc. must get implemented in derived class.
  /// Note that this class assumes an isoparametric formulation, i.e. that
  /// the scalar unknown is interpolated using the same shape funcitons
  /// as the position.
  //=============================================================
  template <unsigned DIM>
  class UnsteadyHeatEquations : public virtual FiniteElement
  {

  public:

    /// \short Function pointer to source function fct(t,x,f(x,t)) -- 
    /// x is a Vector! 
    typedef void (*UnsteadyHeatSourceFctPt)(const double& time,
					    const Vector<double>& x,
					    double& u);


    /// \short Constructor: Initialises the Source_fct_pt to null and 
    /// sets flag to use ALE formulation of the equations.
    UnsteadyHeatEquations() : Source_fct_pt(0), ALE_is_disabled(false) {}
 

    /// Broken copy constructor
    UnsteadyHeatEquations(const UnsteadyHeatEquations& dummy) 
    { 
      BrokenCopy::broken_copy("UnsteadyHeatEquations");
    } 
 
    /// Broken assignment operator
    void operator=(const UnsteadyHeatEquations&) 
    {
      BrokenCopy::broken_assign("UnsteadyHeatEquations");
    }

    /// \short Return the index at which the unknown value
    /// is stored. The default value, 0, is appropriate for single-physics
    /// problems, when there is only one variable, the value that satisfies the
    /// unsteady heat equation. 
    /// In derived multi-physics elements, this function should be overloaded
    /// to reflect the chosen storage scheme. Note that these equations require
    /// that the unknown is always stored at the same index at each node.
    virtual inline unsigned u_index_ust_heat() const {return 0;}
 
    /// \short du/dt at local node n. 
    /// Uses suitably interpolated value for hanging nodes.
    double du_dt_ust_heat(const unsigned &n) const
    {
      // Get the data's timestepper
      TimeStepper* time_stepper_pt= this->node_pt(n)->time_stepper_pt();

      //Initialise dudt
      double dudt=0.0;
   
      //Loop over the timesteps, if there is a non Steady timestepper
      if (!time_stepper_pt->is_steady())
	{
	  //Find the index at which the variable is stored
	  const unsigned u_nodal_index = u_index_ust_heat();
     
	  // Number of timsteps (past & present)
	  const unsigned n_time = time_stepper_pt->ntstorage();
     
	  //Add the contributions to the time derivative
	  for(unsigned t=0;t<n_time;t++)
	    {
	      dudt += time_stepper_pt->weight(1,t)*nodal_value(t,n,u_nodal_index);
	    }
	}
      return dudt;
    }

    /// \short Disable ALE, i.e. assert the mesh is not moving -- you do this
    /// at your own risk!
    void disable_ALE()
    {
      ALE_is_disabled=true;
    }


    /// \short (Re-)enable ALE, i.e. take possible mesh motion into account
    /// when evaluating the time-derivative. Note: By default, ALE is 
    /// enabled, at the expense of possibly creating unnecessary work 
    /// in problems where the mesh is, in fact, stationary. 
    void enable_ALE()
    {
      ALE_is_disabled=false;
    }


    /// Output with default number of plot points
    void output(std::ostream &outfile) 
    {
      unsigned nplot=5;
      output(outfile,nplot);
    }


    /// \short Output FE representation of soln: x,y,u or x,y,z,u at 
    /// n_plot^DIM plot points
    void output(std::ostream &outfile, const unsigned &nplot);

    /// Access function: Pointer to source function
    UnsteadyHeatSourceFctPt& source_fct_pt() {return Source_fct_pt;}


    /// Access function: Pointer to source function. Const version
    UnsteadyHeatSourceFctPt source_fct_pt() const {return Source_fct_pt;}


    /// \short Get source term at continous time t and (Eulerian) position x.
    /// Virtual so it can be overloaded in derived multiphysics elements. 
    virtual inline void get_source_ust_heat(const double& t,
					    const unsigned& ipt,
					    const Vector<double>& x,
					    double& source) const
    {
      //If no source function has been set, return zero
      if(Source_fct_pt==0) {source = 0.0;}
      else
	{
	  // Get source strength
	  (*Source_fct_pt)(t,x,source);
	}
    }

    /// Get flux: flux[i] = du/dx_i
    void get_flux(const Vector<double>& s, Vector<double>& flux) const
    {
      //Find out how many nodes there are in the element
      unsigned n_node = nnode();

      //Find the index at which the variable is stored
      unsigned u_nodal_index = u_index_ust_heat();

      //Set up memory for the shape and test functions
      Shape psi(n_node);
      DShape dpsidx(n_node,DIM);
 
      //Call the derivatives of the shape and test functions
      dshape_eulerian(s,psi,dpsidx);
     
      //Initialise to zero
      for(unsigned j=0;j<DIM;j++) {flux[j] = 0.0;}
   
      // Loop over nodes
      for(unsigned l=0;l<n_node;l++) 
	{
	  //Loop over derivative directions
	  for(unsigned j=0;j<DIM;j++)
	    {                               
	      flux[j] += nodal_value(l,u_nodal_index)*dpsidx(l,j);
	    }
	}
    }


    /// Compute element residual Vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      //Call the generic residuals function with flag set to 0
      //using a dummy matrix argument
      fill_in_generic_residual_contribution_ust_heat(
						     residuals,GeneralisedElement::Dummy_matrix,0);
    }


    /// Return FE representation of function value u(s) at local coordinate s
    inline double interpolated_u_ust_heat(const Vector<double> &s) const
    {
      //Find number of nodes
      unsigned n_node = nnode();

      //Find the index at which the variable is stored
      unsigned u_nodal_index = u_index_ust_heat();

      //Local shape function
      Shape psi(n_node);

      //Find values of shape function
      shape(s,psi);

      //Initialise value of u
      double interpolated_u = 0.0;

      //Loop over the local nodes and sum
      for(unsigned l=0;l<n_node;l++) 
	{
	  interpolated_u += nodal_value(l,u_nodal_index)*psi[l];
	}

      return(interpolated_u);
    }

    /// \short Self-test: Return 0 for OK
    unsigned self_test(){return 0;};


  protected:

    /// \short Shape/test functions and derivs w.r.t. to global coords at 
    /// local coord. s; return  Jacobian of mapping
    virtual double dshape_and_dtest_eulerian_ust_heat(const Vector<double> &s, 
						      Shape &psi, 
						      DShape &dpsidx, 
						      Shape &test, 
						      DShape &dtestdx) const=0;


    /// \short Shape/test functions and derivs w.r.t. to global coords at 
    /// integration point ipt; return  Jacobian of mapping
    virtual double dshape_and_dtest_eulerian_at_knot_ust_heat(
							      const unsigned &ipt, 
							      Shape &psi, 
							      DShape &dpsidx,
							      Shape &test, 
							      DShape &dtestdx)
      const=0;

    /// \short Compute element residual Vector only (if flag=and/or element 
    /// Jacobian matrix 
    virtual void fill_in_generic_residual_contribution_ust_heat(
								Vector<double> &residuals, DenseMatrix<double> &jacobian, 
								unsigned flag); 

    /// Pointer to source function:
    UnsteadyHeatSourceFctPt Source_fct_pt;

    /// \short Boolean flag to indicate if ALE formulation is disabled when 
    /// time-derivatives are computed. Only set to true if you're sure
    /// that the mesh is stationary.
    bool ALE_is_disabled;

  };






  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////



  //======================================================================
  /// QUnsteadyHeatElement elements are linear/quadrilateral/brick-shaped 
  /// UnsteadyHeat elements with isoparametric interpolation for the function.
  //======================================================================
  template <unsigned DIM, unsigned NNODE_1D>
  class QUnsteadyHeatElement : public virtual QElement<DIM,NNODE_1D>,
			       public virtual UnsteadyHeatEquations<DIM>
  {
  private:

    /// \short Static array of ints to hold number of variables at 
    /// nodes: Initial_Nvalue[n]
    static const unsigned Initial_Nvalue;
 
  public:
 
    ///\short  Constructor: Call constructors for QElement and 
    /// UnsteadyHeat equations
    QUnsteadyHeatElement() : QElement<DIM,NNODE_1D>(), 
			     UnsteadyHeatEquations<DIM>()
    { }

    /// Broken copy constructor
    QUnsteadyHeatElement(const QUnsteadyHeatElement<DIM,NNODE_1D>& dummy) 
    { 
      BrokenCopy::broken_copy("QUnsteadyHeatElement");
    } 
 
    /// Broken assignment operator
    void operator=(const QUnsteadyHeatElement<DIM,NNODE_1D>&) 
    {
      BrokenCopy::broken_assign("QUnsteadyHeatElement");
    }

    /// \short  Required  # of `values' (pinned or dofs) 
    /// at node n
    inline unsigned required_nvalue(const unsigned &n) const 
    {return Initial_Nvalue;}

    /// \short Output function:  
    ///  x,y,u   or    x,y,z,u
    void output(std::ostream &outfile)
    {UnsteadyHeatEquations<DIM>::output(outfile);}


    ///  \short Output function:  
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(std::ostream &outfile, const unsigned &n_plot)
    {UnsteadyHeatEquations<DIM>::output(outfile,n_plot);}

  protected:

    /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
    inline double dshape_and_dtest_eulerian_ust_heat(const Vector<double> &s, 
						     Shape &psi, 
						     DShape &dpsidx, 
						     Shape &test, 
						     DShape &dtestdx) const;
 

    /// \short Shape/test functions and derivs w.r.t. to global coords at 
    /// integration point ipt; return  Jacobian of mapping
    inline double dshape_and_dtest_eulerian_at_knot_ust_heat(const unsigned &ipt, 
							     Shape &psi, 
							     DShape &dpsidx,
							     Shape &test, 
							     DShape &dtestdx)
      const;

  };


  //Inline functions:


  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double QUnsteadyHeatElement<DIM,NNODE_1D>::
  dshape_and_dtest_eulerian_ust_heat(const Vector<double> &s,
				     Shape &psi, 
				     DShape &dpsidx,
				     Shape &test, 
				     DShape &dtestdx) const
  {
    //Call the geometrical shape functions and derivatives  
    double J = this->dshape_eulerian(s,psi,dpsidx);
 
    //Loop over the test functions and derivatives and set them equal to the
    //shape functions
    for(unsigned i=0;i<NNODE_1D;i++)
      {
	test[i] = psi[i]; 
	for(unsigned j=0;j<DIM;j++)
	  {
	    dtestdx(i,j) = dpsidx(i,j);
	  }
      }
 
    //Return the jacobian
    return J;
  }


  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned DIM,unsigned NNODE_1D>
  double QUnsteadyHeatElement<DIM,NNODE_1D>::
  dshape_and_dtest_eulerian_at_knot_ust_heat(
					     const unsigned &ipt,
					     Shape &psi, 
					     DShape &dpsidx,
					     Shape &test, 
					     DShape &dtestdx) const
  {
    //Call the geometrical shape functions and derivatives  
    double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);

    //Set the test functions equal to the shape functions 
    //(sets internal pointers)
    test = psi;
    dtestdx = dpsidx;

    //Return the jacobian
    return J;
  }

  //======================================================================
  // Set the data for the number of Variables at each node
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  const unsigned QUnsteadyHeatElement<DIM,NNODE_1D>::Initial_Nvalue = 1;

  //======================================================================
  /// Compute element residual Vector and/or element Jacobian matrix 
  /// 
  /// flag=1: compute both
  /// flag=0: compute only residual Vector
  ///
  /// Pure version without hanging nodes
  //======================================================================
  template <unsigned DIM>
  void  UnsteadyHeatEquations<DIM>::
  fill_in_generic_residual_contribution_ust_heat(Vector<double> &residuals, 
						 DenseMatrix<double> &jacobian, 
						 unsigned flag) 
  {
    //Find out how many nodes there are
    unsigned n_node = nnode();
  
    //Find the index at which the variable is stored
    unsigned u_nodal_index = u_index_ust_heat();

    //Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);
 
    //Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();
   
    //Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    //Integers to hold the local equation and unknowns
    int local_eqn=0; //local_unknown=0;

    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {
	//Assign values of s
	for(unsigned i=0;i<DIM;i++) s[i] = integral_pt()->knot(ipt,i);

	//Get the integral weight
	double w = integral_pt()->weight(ipt);

	//Call the derivatives of the shape and test functions
	double J = 
	  dshape_and_dtest_eulerian_at_knot_ust_heat(ipt,psi,dpsidx,test,dtestdx);

	//Premultiply the weights and the Jacobian
	double W = w*J;

	//Allocate memory for local quantities and initialise to zero
	double interpolated_u=0.0;
	double dudt=0.0;
	Vector<double> interpolated_x(DIM,0.0);
	Vector<double> interpolated_dudx(DIM,0.0);
	Vector<double> mesh_velocity(DIM,0.0);

	//Calculate function value and derivatives:
	// Loop over nodes
	for(unsigned l=0;l<n_node;l++) 
	  {
	    //Calculate the value at the nodes
	    double u_value = raw_nodal_value(l,u_nodal_index);
	    interpolated_u += u_value*psi(l);
	    dudt += du_dt_ust_heat(l)*psi(l);
	    // Loop over directions
	    for(unsigned j=0;j<DIM;j++)
	      {
		interpolated_x[j] += raw_nodal_position(l,j)*psi(l);
		interpolated_dudx[j] += u_value*dpsidx(l,j);
	      }
	  }
	
	//Get source function
	//-------------------
	double source;
	get_source_ust_heat(time(),ipt,interpolated_x,source);
	
	// Assemble residuals and Jacobian
	//--------------------------------
       
	// Loop over the test functions
	for(unsigned l=0;l<n_node;l++)
	  {
	    local_eqn = nodal_local_eqn(l,u_nodal_index);
	    /*IF it's not a boundary condition*/
	    if(local_eqn >= 0)
	      {
		
		// sign function with discontinous jump at u=0
		double sign;
		if (interpolated_u > 0.0) {sign = +1.0;}
		else {sign = -1.0;}
		
		//For equation dudt = u*f(u) + source
		// residuals[local_eqn] += (dudt - interpolated_u*sign - source)*test(l)*W;

		// for equation dudt = f(u) + source
		residuals[local_eqn] += (dudt + 10*sign - 0.3 - source)*test(l)*W;
		
	      }
	  }
      }

  }
  // End of loop over integration points
   


  //======================================================================
  /// Output function:
  ///
  ///   
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  template <unsigned DIM>
  void  UnsteadyHeatEquations<DIM>::output(std::ostream &outfile, 
					   const unsigned &nplot)
  {
    //Vector of local coordinates
    Vector<double> s(DIM);

    unsigned num_plot_points=nplot_points(nplot);
    for (unsigned iplot=0;iplot<num_plot_points;iplot++)
      {
	// Get local coordinates of plot point
	get_s_plot(iplot,nplot,s);

	// Output time and interpolated u
	outfile << time_pt()->time() << " " << interpolated_u_ust_heat(s) << std::endl;   
      }

  }

}

#endif