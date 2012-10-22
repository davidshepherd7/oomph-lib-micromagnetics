#ifndef OOMPH_IMPLICIT_LLG_EXCHANGE_PROBLEM_H
#define OOMPH_IMPLICIT_LLG_EXCHANGE_PROBLEM_H

#include "generic.h"
#include "../boundary_element_handler.h"
#include "../micromagnetics_element.h"
#include "../micromagnetics_element.cc" //??ds shouldn't need...
#include "../vector_helpers.h"

namespace oomph
{

  // ============================================================
  ///
  // ============================================================
  template<class ELEMENT>
  class ImplicitLLGExchangeProblem : public Problem
  {
  public:

    // Function pointer for initial magnetisation.
    typedef void (*InitialMFctPt)(const double& t, const Vector<double> &x,
                                  Vector<double> &m);

    // Function pointer for applied field.
    typedef typename ELEMENT::TimeSpaceToDoubleVectFctPt AppliedFieldFctPt; 
    
    /// Default constructor
    ImplicitLLGExchangeProblem(const unsigned &nx,
                               const unsigned &ny,
                               const double &lx,
                               const double &ly,
                               AppliedFieldFctPt applied_field_pt)
      
    {
      
      // Create timestepper
      BDF<2>* ts_pt = new BDF<2>;
      add_time_stepper_pt(ts_pt); 

      // Create mesh
      mesh_pt() = new SimpleRectangularTriMesh<ELEMENT>(nx,ny,lx,ly,ts_pt); 

      // Set magnetic parameters: all fields zero except exchange which has
      // coefficient 1, |M| = 1.
      ExchangeOnlyParameters.set_simple_llg_parameters();

      // Fiddle with exchange strength if needed.
      ExchangeOnlyParameters.exchange_constant() /= 5;
      //??ds complicated because of normalisation. To get |Hex| = 1 we have
      //to get the exchange constant to mu0. Here we want |Hex| = 0.2 so
      //divide by 5...
      ExchangeOnlyParameters.gilbert_damping() = 0.1;
      
      ExchangeOnlyParameters.output(std::cout);

      // Fiddle with newton paramters if needed
      // max_residuals() *= 2; 
      
      // Get a pointer to an element for use later
      Ele_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0));

      // Find out the problem dimensions
      Dim = Ele_pt->nodal_dimension();

      // //??ds boundary conditions - pin boundary m for now to zero...
      // for(unsigned b=0, nb= mesh_pt()->nboundary(); b < nb; b++)
      //   {
      //     for(unsigned nd=0, nnd=mesh_pt()->nboundary_node(b); nd<nnd; nd++)
      //       {
      //         Node* nd_pt = mesh_pt()->boundary_node_pt(b,nd);
      //         for(unsigned j=0; j<Dim; j++)
      //           {
      //             nd_pt->pin(m_index(j)); 
      //             nd_pt->set_value(m_index(j),0.0); 
      //           }
      //       }
      //   }

      // Finish off elements
      for(unsigned i=0; i< mesh_pt()->nelement(); i++)
        {
          ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

          // Set pointer to continuous time
          elem_pt->time_pt() = time_pt();

          // Set values for magnetic parameters
          elem_pt->magnetic_parameters_pt() = &ExchangeOnlyParameters;

          // Set pointer for an applied field
          elem_pt->applied_field_pt() = applied_field_pt; 
        }

      // Pin all phi dofs...??ds remove them from element?
      for(unsigned nd=0, nnode=mesh_pt()->nnode(); nd<nnode; nd++)
        {
          Node* nd_pt = mesh_pt()->node_pt(nd);
          
          nd_pt->pin(phi_index());
          nd_pt->pin(phi_1_index());
          nd_pt->set_value(phi_index(),0.0);
          nd_pt->set_value(phi_1_index(),0.0); 
        }
      
      // Do equation numbering
      std::cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
    }
    
    /// Destructor
    ~ImplicitLLGExchangeProblem()
    {
      // mesh is cleaned up by problem base class
      // timestepper is cleaned up by problem base class
    }

    /// Renormalise magnetisation to 1 (needed with BDF2)
    void renormalise_magnetisation()
    {
      for(unsigned nd=0; nd<mesh_pt()->nnode(); nd++)
        {
          Node* nd_pt = mesh_pt()->node_pt(nd);
          Vector<double> m_values(3,0.0);
          for(unsigned j=0; j<3; j++) m_values[j] = nd_pt->value(m_index(j));
          VectorOps::normalise(m_values);
          for(unsigned j=0; j<3; j++) nd_pt->set_value(m_index(j),m_values[j]);
        }
    }
    

    void actions_after_newton_solve()
    {
      // If we're using BDF we need to keep M normalised.
      if(dynamic_cast<BDF<2>* >(time_stepper_pt()) != 0)
        {
          renormalise_magnetisation(); 
        } 
    }
    

    /// Output solution
    void doc_solution(DocInfo &doc_info) const
    {
      std::ofstream some_file;
      char filename[100];

      // Number of plot points
      unsigned npts = 1;

      // Output solution with specified number of plot points per element
      sprintf(filename,"%s/soln%i.dat", doc_info.directory().c_str(),
              doc_info.number());
      some_file.open(filename);
      mesh_pt()->output(some_file,npts);
      some_file.close();
      
      // Increment number used as output label
      doc_info.number()++; 
    }
    
    /// Set up an initial M
    void set_initial_condition(const InitialMFctPt initial_m_pt);
    
    
    // Access functions
    // ============================================================
    
    unsigned m_index(const unsigned &j) const
    {return Ele_pt->m_index_micromag(j);} 

    unsigned phi_index() const
    {return Ele_pt->phi_index_micromag();} 

    unsigned phi_1_index() const
    {return Ele_pt->phi_1_index_micromag();} 
    
    
  private:
    
    /// Magnetic parameters with all but exchange field zero.
    MagneticParameters ExchangeOnlyParameters; 

    /// A pointer to an element (used to get index numbers)
    const ELEMENT* Ele_pt;
    
    /// The problem dimension
    unsigned Dim;
        
    /// Inaccessible copy constructor
    ImplicitLLGExchangeProblem(const ImplicitLLGExchangeProblem & dummy)
    {BrokenCopy::broken_copy("ImplicitLLGExchangeProblem");}
    
    /// Inaccessible assignment operator
    void operator=(const ImplicitLLGExchangeProblem &dummy)
    {BrokenCopy::broken_assign("ImplicitLLGExchangeProblem");}

  };

  //======================================================================
  /// Set up the initial conditions
  //======================================================================
  template<class ELEMENT>
  void ImplicitLLGExchangeProblem<ELEMENT>::
  set_initial_condition(const InitialMFctPt initial_m_pt)
  {
    // Backup time in global Time object
    double backed_up_time=this->time_pt()->time();

    // Past history needs to be established for t=time0-deltat, ...
    // Then provide current values (at t=time0) which will also form
    // the initial guess for the first solve at t=time0+deltat

    // Get M indicies
    Vector<unsigned> m_index_micromag(3,0);
    ELEMENT* elem_pt = dynamic_cast<ELEMENT* >(this->mesh_pt()->element_pt(0));
    for(unsigned i=0; i<3; i++)
      {
        m_index_micromag[i] = elem_pt->m_index_micromag(i);
      }

    // Find number of nodes in mesh
    unsigned num_nod = this->mesh_pt()->nnode();

    // Set continuous times at previous timesteps:
    int nprev_steps=this->time_stepper_pt()->nprev_values();
    Vector<double> prev_time(nprev_steps+1);
    for (int t=nprev_steps;t>=0;t--)
      {
        prev_time[t]=this->time_pt()->time(t);
      }

    // Loop over current & previous timesteps
    for (int t=nprev_steps;t>=0;t--)
      {
        // Continuous time
        double time = prev_time[t];
        std::cout << "setting IC at time =" << time << std::endl;

        // Loop over the nodes to set initial values everywhere
        for (unsigned n=0;n<num_nod;n++)
          {
            unsigned dim = mesh_pt()->node_pt(n)->ndim(); 

            // Get initial value of m from inputs
            Vector<double> m(3,0.0), x(dim,0.0);
            this->mesh_pt()->node_pt(n)->position(t,x);
            initial_m_pt(time,x,m);

            // Set initial condition on m
            for(unsigned i=0; i<3; i++)
              this->mesh_pt()->node_pt(n)->set_value(t,m_index_micromag[i],m[i]);
          }
      }

    // Reset backed up time for global timestepper
    this->time_pt()->time()=backed_up_time;
  }

    
}

#endif
