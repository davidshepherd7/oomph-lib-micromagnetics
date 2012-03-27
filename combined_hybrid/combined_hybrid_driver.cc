#ifndef OOMPH_COMBINED_HYBRID_DRIVER_H
#define OOMPH_COMBINED_HYBRID_DRIVER_H

/*
  description of file goes here
*/

#include "generic.h"
#include "../micromagnetics_boundary_element.h"

// Mesh
#include "meshes/rectangular_quadmesh.h"

// My variable gauss quadrature header
#include "../variable_order_quadrature.h"

using namespace oomph;
using namespace MathematicalConstants;

inline double dot3(const Vector<double>& a, const Vector<double>& b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

namespace Inputs
{
  double sat_mag = 1.0;
  double llg_damping = 0.05;
  double llg_precession = 1.0;

  void applied_field(const double& t, const Vector<double>& x, Vector<double>& H_app)
  {
    std::fill(H_app.begin(), H_app.end(), 0.0);
  }

  void cryst_anis_field(const double& t, const Vector<double>& x,
			const Vector<double>& M ,Vector<double>& H_ca)
  {
    Vector<double> Easy_axis(3,0.0); Easy_axis[0] = 1.0;
    double magnitude = dot3(Easy_axis,M);
    H_ca = Easy_axis;
    std::for_each(H_ca.begin(), H_ca.end(), [magnitude](double& elem) {elem*=magnitude;});
  }
}

namespace oomph
{

  //======================================================================
  /// A problem class to test the combination of the hybrid BEM/FEM for
  /// magnetostatic fields and the LLG equations for micromagnetics.
  //======================================================================
  template<class BULK_ELEMENT,
	   template<class BULK_ELEMENT,unsigned DIM> class FACE_ELEMENT,
	   unsigned DIM>
  class TwoDHybridProblem : public Problem
  {

  public:

    /// Constructor
    TwoDHybridProblem();

    /// Destructor (empty -- once the problem is done with the program is over)
    ~TwoDHybridProblem(){};

    /// Doc the solution
    void doc_solution(DocInfo& doc_info);


    /// Build the meshes of face elements and corner elements
    void build_face_mesh(Mesh* face_mesh_pt) const;

    /// Get the boundary element matrix (similar to problem::get_jacobian)
    void build_boundary_matrix();

    /// Set initial condition (incl previous timesteps)
    void set_initial_condition();


    /// \short Get the mapping between the global equation numbering and
    /// the boundary equation numbering.
    void create_global_boundary_equation_number_map();

    /// Get the boundary equation number from the global equation number
    unsigned convert_global_to_boundary_equation_number(const unsigned &global_num);

    /// Update the values of phi_2 on the boundary
    void update_boundary_phi_2();


    /// Access to the face mesh pointer
    Mesh* face_mesh_pt() {return Face_mesh_pt;}

    /// Access to the boundary matrix
    DenseDoubleMatrix* boundary_matrix_pt()
    {return &Boundary_matrix;}

  private:

    /// The map between the global equation numbering and the boundary equation numbering.
    std::map<unsigned,unsigned> Global_boundary_equation_num_map;

    /// Mesh containing the face elements
    Mesh* Face_mesh_pt;

    /// Doc info object
    DocInfo Doc_info;

    /// Matrix to store the relationship between phi_1 and phi_2 on the boundary
    DenseDoubleMatrix Boundary_matrix;

    // Update the problem before Newton convergence check
    void actions_before_newton_convergence_check()
    {update_boundary_phi_2();}

    /// Update the problem specs before solve
    // Nothing to do here since no dirichlet boundaries.
    void actions_before_newton_solve(){};

    /// Update the problem specs after solve (empty)
    void actions_after_newton_solve(){};

    /// Update the problem specs after solve (empty)
    void actions_after_implicit_timestep(){};

    /// Update the problem specs before next timestep (boundary conditions)
    void actions_before_implicit_timestep();

  }; // end of problem class


//======================================================================
/// Constructor
//======================================================================
template<class BULK_ELEMENT, template<class,unsigned> class FACE_ELEMENT, unsigned DIM>
TwoDHybridProblem<BULK_ELEMENT,FACE_ELEMENT,DIM>::
TwoDHybridProblem()
{
  // Allocate steady state timestepper
  add_time_stepper_pt(new BDF<2>);

  // Build rectangular mesh
  unsigned nx = 20, ny = 20;
  mesh_pt() = new RectangularQuadMesh<BULK_ELEMENT>
    (nx,ny,2.0,2.0,time_stepper_pt());

  // Create face elements on all boundaries and add to face mesh
  Face_mesh_pt = new Mesh;
  build_face_mesh(face_mesh_pt());

  // Create a variable integration scheme to do the BEM integration
  QVariableOrderGaussLegendre<DIM-1>* bem_quadrature_scheme_pt
    = new QVariableOrderGaussLegendre<DIM-1>;

  // Loop over elements in face mesh to set mesh pointer
  unsigned n_face_element = face_mesh_pt()->nelement();
  for(unsigned i=0;i<n_face_element;i++)
    {
      // Upcast from GeneralisedElement to the face element
      FACE_ELEMENT<BULK_ELEMENT,DIM>* face_elem_pt =
	dynamic_cast<FACE_ELEMENT<BULK_ELEMENT,DIM> *>(face_mesh_pt()->element_pt(i));

      // Set boundary mesh pointer in element
      face_elem_pt->set_boundary_mesh_pt(face_mesh_pt());

      // Set the integration scheme pointer in element
      face_elem_pt->set_integration_scheme(bem_quadrature_scheme_pt);
    }

  // Loop over elements in bulk mesh to set function pointers
  unsigned n_bulk_element = mesh_pt()->nelement();
  for(unsigned i=0; i<n_bulk_element; i++)
    {
      // Upcast from GeneralisedElement to the present element
      BULK_ELEMENT *elem_pt = dynamic_cast<BULK_ELEMENT*>(mesh_pt()->element_pt(i));

      // Set pointer to continous time
      elem_pt->time_pt() = time_pt();

      // Set the function pointers for parameters
      elem_pt->applied_field_pt() = &Inputs::applied_field; //??ds fix this to use proper encapsulation asap
      elem_pt->cryst_anis_field_pt() = &Inputs::cryst_anis_field;
      elem_pt->sat_mag_pt() = &Inputs::sat_mag;
      elem_pt->llg_damp_pt() = &Inputs::llg_damping;
      elem_pt->llg_precess_pt() = &Inputs::llg_precession;
    }

  // Setup equation numbering scheme
  std::cout << "FEM number of equations: " << assign_eqn_numbers() << std::endl;

} // end of constructor



  //======================================================================
  /// Output function
  //======================================================================
template<class BULK_ELEMENT, template<class,unsigned> class FACE_ELEMENT, unsigned DIM>
void TwoDHybridProblem<BULK_ELEMENT,FACE_ELEMENT,DIM>::
doc_solution(DocInfo& doc_info)
{
  // Number of plot points
  unsigned npts=5;

  // File set up
  std::ofstream some_file;
  char filename[100];
  sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
	  doc_info.number());

  // Output
  some_file.open(filename);
  mesh_pt()->output(some_file,npts);
  some_file.close();

} // end of doc


  //======================================================================
  /// Given the global equation number return the boundary equation number. Most
  /// of the function is an error check.
  //======================================================================
template<class BULK_ELEMENT, template<class,unsigned> class FACE_ELEMENT, unsigned DIM>
unsigned TwoDHybridProblem<BULK_ELEMENT,FACE_ELEMENT,DIM>::
convert_global_to_boundary_equation_number(const unsigned &global_num)
{
#ifdef PARANOID
  // Get the location of the global_num key in an iterator
  std::map<unsigned,unsigned>::iterator it
    = Global_boundary_equation_num_map.find(global_num);

  // If the iterator is placed at the end the given global equation number is
  // not in the map, so return an error.
  if(it == Global_boundary_equation_num_map.end())
    {
      std::ostringstream error_stream;
      error_stream << "Global equation number " << global_num
		   << " is not in the global to boundary map.";
      throw OomphLibError(error_stream.str(),
			  "TwoDHybridProblem::convert_global_to_boundary_equation_number",
			  OOMPH_EXCEPTION_LOCATION);
    }
#endif
  return ((*Global_boundary_equation_num_map.find(global_num)).second);
}


//======================================================================
/// Build the mesh of face elements.
//======================================================================
template<class BULK_ELEMENT, template<class,unsigned> class FACE_ELEMENT, unsigned DIM>
void TwoDHybridProblem<BULK_ELEMENT,FACE_ELEMENT,DIM>::
build_face_mesh(Mesh* face_mesh_pt) const
{
  // Create a set to temporarily store the list of boundary nodes
  // (we do it via a set because sets automatically detect duplicates)
  std::set<Node*> node_set, corner_node_set;
  std::set<Node*>::iterator it, c_it;

  // Loop over the boundaries
  unsigned n_boundary = mesh_pt()->nboundary();
  for(unsigned b=0; b<n_boundary; b++)
    {
      // Loop over the boundary nodes on boundary b making a set of nodes
      unsigned n_bound_node = mesh_pt()->nboundary_node(b);
      for(unsigned n=0;n<n_bound_node;n++)
	node_set.insert(mesh_pt()->boundary_node_pt(b,n));

      // Loop over the elements on boundary b creating face elements
      unsigned n_bound_element = mesh_pt()->nboundary_element(b);
      for(unsigned e=0;e<n_bound_element;e++)
	{
	  // Create the corresponding FaceElement
	  FACE_ELEMENT<BULK_ELEMENT,DIM>* face_element_pt = new FACE_ELEMENT<BULK_ELEMENT,DIM>
	    (mesh_pt()->boundary_element_pt(b,e),
	     mesh_pt()->face_index_at_boundary(b,e));

	  // Add the face element to the face mesh
	  face_mesh_pt->add_element_pt(face_element_pt);
	}
    }

  // Iterate over all nodes in the set and add them to the face mesh
  for(it=node_set.begin(); it!=node_set.end(); it++)
    face_mesh_pt->add_node_pt(*it);
}


//======================================================================
/// Create the map between the global equation numbering system and the
/// boundary equation numbering system.
//======================================================================
template<class BULK_ELEMENT, template<class,unsigned> class FACE_ELEMENT, unsigned DIM>
void TwoDHybridProblem<BULK_ELEMENT,FACE_ELEMENT,DIM>::
create_global_boundary_equation_number_map()
{
  // Initialise the map
  Global_boundary_equation_num_map.clear();

  // Loop over boundary nodes assigning a boundary equation number to each.
  unsigned n_boundary_node = this->face_mesh_pt()->nnode(), k=0;
  for(unsigned i_node=0; i_node<n_boundary_node; i_node++)
    {
      // Get global equation number for phi_2
      unsigned global_eqn_number = this->face_mesh_pt()->
	node_pt(i_node)->eqn_number(0);

      // Set up the pair ready to input with key="global equation number" and
      // value ="boundary equation number"=k.
      std::pair<unsigned,unsigned> input_pair = std::make_pair(global_eqn_number,k);

      // Add entry to map and store whether this was a new addition
      bool new_addition = (Global_boundary_equation_num_map.insert(input_pair)
			   ).second;

      // Increment k if this was a new addition to the map
      if(new_addition) k++;
    }
}

//=============================================================================
/// Get the fully assembled boundary matrix in dense storage.
//=============================================================================
template<class BULK_ELEMENT, template<class,unsigned> class FACE_ELEMENT, unsigned DIM>
void TwoDHybridProblem<BULK_ELEMENT,FACE_ELEMENT,DIM>::
build_boundary_matrix()
{
  // Create the mapping from global to boundary equations
  create_global_boundary_equation_number_map();

  // get the number of nodes in the boundary problem
  unsigned long n_node = face_mesh_pt()->nnode();

  // Initialise and resize the boundary matrix
  Boundary_matrix.resize(n_node,n_node);
  Boundary_matrix.initialise(0.0);

  // Loop over all elements in the face mesh
  unsigned long n_face_element = face_mesh_pt()->nelement();
  for(unsigned long e=0;e<n_face_element;e++)
    {
      // Get the pointer to the element (and cast to FiniteElement)
      FiniteElement* elem_pt =
	dynamic_cast < FiniteElement* > (face_mesh_pt()->element_pt(e));

      // Find number of nodes in the element
      unsigned long n_element_node = elem_pt->nnode();

      // Set up a matrix and dummy residual vector
      DenseMatrix<double> element_boundary_matrix(n_element_node,n_node);
      Vector<double> dummy(0);

      // Fill the matrix
      assembly_handler_pt()
	->get_jacobian(elem_pt,dummy,element_boundary_matrix);

      // Loop over the nodes in this element
      for(unsigned l=0;l<n_element_node;l++)
	{
	  // Get the boundary equation (=node) number from the global one
	  unsigned l_number =
	    this->convert_global_to_boundary_equation_number
	    (elem_pt->node_pt(l)->eqn_number(0));

	  //std::cout << "target node = " << l << std::endl;

	  // Loop over all nodes in the mesh and add contributions from this element
	  for(unsigned long source_node=0; source_node<n_node; source_node++)
	    {
	      unsigned source_number =
		this->convert_global_to_boundary_equation_number
		(face_mesh_pt()->node_pt(source_node)->eqn_number(0));

	      Boundary_matrix(l_number,source_number)
		+= element_boundary_matrix(l,source_node);
	    }
	}
    }
}

//======================================================================
/// Get the new boundary condition on phi_2 using the boundary element matrix
/// and phi_1.
//======================================================================
template<class BULK_ELEMENT, template<class,unsigned> class FACE_ELEMENT, unsigned DIM>
void TwoDHybridProblem<BULK_ELEMENT,FACE_ELEMENT,DIM>::
update_boundary_phi_2()
{
  // Get the index of phi_1 and phi_2
  BULK_ELEMENT* elem_pt = (dynamic_cast<BULK_ELEMENT*>(mesh_pt()->element_pt(0)));
  const unsigned phi_1_index = elem_pt->phi_1_index_micromag();
  const unsigned phi_index = elem_pt->phi_index_micromag();

  // Loop over all (target) nodes on the boundary
  unsigned n_boundary_node = face_mesh_pt()->nnode();
  for(unsigned target_node=0; target_node<n_boundary_node; target_node++)
    {
      // Get a pointer to the target node
      Node* target_node_pt = face_mesh_pt()->node_pt(target_node);

      // Get boundary equation number for this target node
      unsigned target_number = convert_global_to_boundary_equation_number
	(target_node_pt->eqn_number(0));

      // Double to store the value of total phi during computation
      double target_phi_value = 0;

      // Loop over all source nodes adding contribution from each
      for(unsigned source_node=0; source_node<n_boundary_node; source_node++)
	{
	  // Get a pointer to the source node
	  Node* source_node_pt = face_mesh_pt()->node_pt(source_node);

	  // Get boundary equation number for this source node
	  unsigned source_number = convert_global_to_boundary_equation_number
	    (source_node_pt->eqn_number(0));

	  // Add the contribution to total phi at the target node due to
	  // the source node (relationship is given by the boundary matrix).
	  //??ds check consistency of boundary matrix numbering
	  target_phi_value += Boundary_matrix(target_number,source_number)
	    * source_node_pt->value(phi_1_index);
	}
      // Save the total into the target node
      target_node_pt->set_value(phi_index,target_phi_value);
    }
}


//======================================================================
/// Set up the intial conditions
//======================================================================
template<class BULK_ELEMENT, template<class,unsigned> class FACE_ELEMENT, unsigned DIM>
void TwoDHybridProblem<BULK_ELEMENT,FACE_ELEMENT,DIM>::
set_initial_condition()
{
  // Set intial conditions on M
  Vector<double> initial_solution(3,0.0);
  initial_solution[0] = 1.0;


  // Backup time in global Time object
  double backed_up_time=time_pt()->time();

  // Past history needs to be established for t=time0-deltat, ...
  // Then provide current values (at t=time0) which will also form
  // the initial guess for the first solve at t=time0+deltat

  //Find number of nodes in mesh
  unsigned num_nod = mesh_pt()->nnode();

  // Set continuous times at previous timesteps:
  int nprev_steps=time_stepper_pt()->nprev_values();
  Vector<double> prev_time(nprev_steps+1);
  for (int t=nprev_steps;t>=0;t--)
    {
      prev_time[t]=time_pt()->time(t);
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
	  // Set initial condition on M, could set others here using other i values
	  for(unsigned i=1; i<4; i++)
	    {
	      mesh_pt()->node_pt(n)->
		set_value(t,i,initial_solution[i]);
	    }
	}
    }

  // Reset backed up time for global timestepper
  time_pt()->time()=backed_up_time;


  // Construct the boundary matrix
  build_boundary_matrix();
}


//======================================================================
/// Actions before timestep, we set up the boundary conditions here.
//======================================================================
template<class BULK_ELEMENT, template<class,unsigned> class FACE_ELEMENT, unsigned DIM>
void TwoDHybridProblem<BULK_ELEMENT,FACE_ELEMENT,DIM>::
actions_before_implicit_timestep()
{
  //??ds set up Neumann boundary conditions for phi_1.
}


} // End of oomph namespace

int main(int argc, char* argv[])
{

  // Inputs
  const unsigned dim = 2;
  const unsigned nnode_1d = 2;
  const double dt = 0.01;
  const unsigned nstep = 50;

  // Create the problem
  TwoDHybridProblem< QMicromagElement <dim,nnode_1d>, MicromagFaceElement, dim >
    problem;

  // Initialise timestep, initial conditions
  problem.initialise_dt(dt);
  problem.set_initial_condition();

  // Set up outputs and output initial conditions
  DocInfo doc_info;
  doc_info.set_directory("results");
  doc_info.number()=0;
  problem.doc_solution(doc_info);
  doc_info.number()++;

  /// Check problem
  if(!(problem.self_test()==0))
    throw OomphLibError("Problem self_test failed","main",
			OOMPH_EXCEPTION_LOCATION);

  // Timestepping loop
  for(unsigned istep=0; istep<nstep; istep++)
    {
      std::cout << "Timestep " << istep << std::endl;
      std::cout << ", continuous time " << problem.time_pt()->time() << std::endl;

      // Take timestep
      problem.unsteady_newton_solve(dt);

      //Output solution
      problem.doc_solution(doc_info);

      //Increment counter for solutions
      doc_info.number()++;
    }

  return 0;
}

#endif
