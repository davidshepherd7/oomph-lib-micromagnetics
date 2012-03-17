#ifndef OOMPH_HYBRID_BOUNDARY_ELEMENT_DRIVER_H
#define OOMPH_HYBRID_BOUNDARY_ELEMENT_DRIVER_H

#include <map>
#include "generic.h"
#include "../micromagnetics_boundary_element.h"

// Various meshes
#include "meshes/rectangular_quadmesh.h"
#include "meshes/collapsible_channel_mesh.h"


// My variable gauss quadrature header
#include "../variable_order_quadrature.h"

using namespace oomph;
using namespace MathematicalConstants;

namespace oomph
{
  /////////////////////////////////////// stolen from NS collapsible channel:
  //=========================================================================
  class NonOscillatingWall : public GeomObject
  {

  public:

    /// \short Constructor : It's a 2D object, parametrised by
    /// one Lagrangian coordinate. Arguments: height at ends, x-coordinate of
    /// left end and length.
    NonOscillatingWall(const double& h, const double& x_left, const double& l,
		       const double &a) :
      GeomObject(1,2), H(h), Length(l), X_left(x_left), A(a), B(0.0)
    {}

    /// Destructor:  Empty
    ~NonOscillatingWall(){}

    /// \short Position vector at Lagrangian coordinate zeta
    /// at time level t.
    void position(const unsigned& t, const Vector<double>&zeta,
		  Vector<double>& r) const
    {
      using namespace MathematicalConstants;
      // Position vector
      r[0] = zeta[0]+X_left - A*B*sin(2.0*Pi*zeta[0]/Length);
      r[1] = H + A*((Length-zeta[0])*zeta[0])/pow(0.5*Length,2);
    }

    /// \short "Current" position vector at Lagrangian coordinate zeta
    void position(const Vector<double>&zeta, Vector<double>& r) const
    {position (0, zeta, r);}

    /// Number of geometric Data in GeomObject: None.
    unsigned ngeom_data() const {return 0;}

  private:
    /// Height at ends
    double H;
    /// Length
    double Length;
    /// x-coordinate of left end
    double X_left;
    /// Amplitude of bump in wall
    double A;
    /// Relative amplitude of horizontal wall motion
    double B;

  }; // end of non-oscillating wall

  //=====start_of_TwoDMicromagProblem===============================================
  /// A problem class to solve the LLG equation using a hybrid BEM/FEM.
  //========================================================================
  template<class BULK_ELEMENT,
	   template<class BULK_ELEMENT,unsigned DIM> class FACE_ELEMENT,
	   unsigned DIM>
  class TwoDMicromagProblem : public Problem
  {

  public:

    /// Constructor
    TwoDMicromagProblem(const unsigned& n_x, const unsigned& n_y);

    /// Destructor (empty -- all the cleanup is done in the base class)
    ~TwoDMicromagProblem(){};

    /// Doc the solution
    void doc_solution(DocInfo& doc_info);

    /// Build the meshes of face elements and corner elements
    void build_face_mesh(Mesh* face_mesh_pt, Mesh* corner_mesh_pt) const;

    /// Acess to the face mesh pointer
    Mesh* face_mesh_pt() {return Face_mesh_pt;}

    /// Acess to the corner mesh pointer
    Mesh* corner_mesh_pt() {return Corner_mesh_pt;}

    /// Get the boundary element matrix (similar to problem::get_jacobian)
    void get_boundary_matrix();

    /// \short Get the mapping between the global equation numbering and
    /// the boundary equation numbering.
    void create_global_boundary_equation_number_map();

    /// Get the boundary equation number from the global equation number
    unsigned convert_global_to_boundary_equation_number(const unsigned &global_num);

    /// Update the values of phi_2 on the boundary
    void update_boundary_phi_2();

    /// Return pointer to the boundary matrix
    DenseDoubleMatrix* boundary_matrix_pt()
    {return &Boundary_matrix;}

  private:

    /// The map between the global equation numbering and the boundary equation numbering.
    std::map<unsigned,unsigned> Global_boundary_equation_num_map;

    /// Pointer to control node at which the solution is documented
    Node* Control_node_pt;

    /// Mesh containing the face elements
    Mesh* Face_mesh_pt;

    /// Mesh containing the point corner elements (basiclly just a list of elements)
    Mesh* Corner_mesh_pt;

    /// List of the point elements at (pre-discretisation) sharp corners.
    Vector<MicromagCornerAngleElement<BULK_ELEMENT,DIM>* > Corner_elements;

    /// Doc info object
    DocInfo Doc_info;

    /// Matrix to store the relationship between phi_1 and phi_2 on the boundary
    DenseDoubleMatrix Boundary_matrix;

    //  /// Trace file
    //  std::ofstream Trace_file;

    // Update the problem before Newton convergence check
    void actions_before_newton_convergence_check()
    {update_boundary_phi_2();}

    /// Update the problem specs before solve
    // Nothing to do here since no dirichlet boundaries.
    void actions_before_newton_solve(){};

    /// Update the problem specs after solve
    void actions_after_newton_solve(){};

    /// Update the problem specs after solve (empty)
    void actions_after_implicit_timestep(){};

    /// Update the problem specs before next timestep
    void actions_before_implicit_timestep(){};

    /// Set initial condition (incl previous timesteps) according to specified function.
    void set_initial_condition(){};

  }; // end of problem class

//=====start_of_constructor===============================================
///
//========================================================================
template<class BULK_ELEMENT, template<class,unsigned> class FACE_ELEMENT, unsigned DIM>
TwoDMicromagProblem<BULK_ELEMENT,FACE_ELEMENT,DIM>::
TwoDMicromagProblem(const unsigned& n_x, const unsigned& n_y)
{
  // Allocate steady state timestepper
  add_time_stepper_pt(new Steady<2>);

  // // Build rectangular mesh
  // mesh_pt() = new RectangularQuadMesh<BULK_ELEMENT>
  //   (n_x,n_y,2.0,2.0,time_stepper_pt());

  // Build mesh collapsible channel mesh
  {
    //Create the geometric object that represents the wall Parameters chosen
    // to make it very similar to the 2x2 square used before (if a = 0).
    double height = 2.0, x_left = 2.0/3, length = 2.0/3, a = 0.3;
    GeomObject* Wall_pt = new NonOscillatingWall(height, x_left, length, a);

    // Number of elements and lengths of parts of the mesh
    unsigned nup = unsigned(n_x/3), ndown = nup, ncollapsible = nup;
    double lup = x_left, ldown = 2.0/3, lcollapsible = length, ly = height;
    mesh_pt() = new
      CollapsibleChannelMesh<BULK_ELEMENT>(nup, ncollapsible, ndown, n_y,
					   lup, lcollapsible, ldown, ly,
					   Wall_pt, time_stepper_pt());
  }

  // dump mesh for testing
  unsigned n_nd = mesh_pt()->nnode();
  std::ofstream mesh_plot;
  mesh_plot.open("./mesh_points");
  for(unsigned nd=0; nd<n_nd; nd++)
    {
      mesh_plot << mesh_pt()->node_pt(nd)->x(0) << " "
		<< mesh_pt()->node_pt(nd)->x(1) << std::endl;
    }
  mesh_plot.close();

  // Create face elements on all boundaries and add to face mesh, create point
  // elements at sharp corners of the pre-discretisation object and add to
  // corner mesh.
  Face_mesh_pt = new Mesh;
  Corner_mesh_pt = new Mesh;
  build_face_mesh(face_mesh_pt(),corner_mesh_pt());

  // Loop over elements to set mesh pointer
  unsigned n_element = face_mesh_pt()->nelement();
  for(unsigned i=0;i<n_element;i++)
    {
      // Upcast from GeneralisedElement to the face element
      FACE_ELEMENT<BULK_ELEMENT,DIM>* elem_pt =
	dynamic_cast<FACE_ELEMENT<BULK_ELEMENT,DIM> *>(face_mesh_pt()->element_pt(i));

      // Set boundary mesh pointer in element
      elem_pt->set_boundary_mesh_pt(face_mesh_pt());

      //??ds might need some more pointers later
    }

  // Setup equation numbering scheme
  std::cout << "FEM number of equations: " << assign_eqn_numbers() << std::endl;

} // end of constructor

  //=====================start_of_doc=======================================
  /// ??ds write this!
  //========================================================================
template<class BULK_ELEMENT, template<class,unsigned> class FACE_ELEMENT, unsigned DIM>
void TwoDMicromagProblem<BULK_ELEMENT,FACE_ELEMENT,DIM>::
doc_solution(DocInfo& doc_info)
{

  std::ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts;
  npts=5;

  // Output solution
  sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
	  doc_info.number());
  some_file.open(filename);
  mesh_pt()->output(some_file,npts);
  some_file.close();

  // // Output exact solution
  // //----------------------
  // sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
  // 	    doc_info.number());
  // some_file.open(filename);
  // mesh_pt()->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u);
  // some_file.close();

} // end of doc

  //======start_of_convert_global_to_boundary_equation_number===============
  /// ??ds write this!
  //========================================================================
template<class BULK_ELEMENT, template<class,unsigned> class FACE_ELEMENT, unsigned DIM>
unsigned TwoDMicromagProblem<BULK_ELEMENT,FACE_ELEMENT,DIM>::
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
			  "TwoDMicromagProblem::convert_global_to_boundary_equation_number",
			  OOMPH_EXCEPTION_LOCATION);
    }
#endif
  return ((*Global_boundary_equation_num_map.find(global_num)).second);
}

//======start_of_build_face_mesh==========================================
/// \short Constuct a Mesh of FACE_ELEMENTs along the b-th boundary
/// of the mesh (which contains elements of type BULK_ELEMENT)
// ??ds also create a "corner mesh" which contains only point elements at the
// edges of the boundarys - this works in 2d and for simple-ish meshes ONLY.
//========================================================================
template<class BULK_ELEMENT, template<class,unsigned> class FACE_ELEMENT, unsigned DIM>
void TwoDMicromagProblem<BULK_ELEMENT,FACE_ELEMENT,DIM>::
build_face_mesh(Mesh* face_mesh_pt, Mesh* corner_mesh_pt) const
{
  // Create a set to temporarily store the list of boundary nodes
  // (we do it via a set because sets automatically detect duplicates)
  std::set<Node*> node_set, corner_node_set;
  std::set<Node*>::iterator it, c_it;

  // Loop over the boundaries
  unsigned n_boundary = mesh_pt()->nboundary();
  for(unsigned b=0; b<n_boundary; b++)
    {
      //Loop over the boundary nodes on boundary b making a set of nodes
      unsigned n_bound_node = mesh_pt()->nboundary_node(b);
      for(unsigned n=0;n<n_bound_node;n++)
	node_set.insert(mesh_pt()->boundary_node_pt(b,n));

      //Loop over the elements on boundary b creating face elements
      unsigned n_bound_element = mesh_pt()->nboundary_element(b);
      for(unsigned e=0;e<n_bound_element;e++)
	{
	  //Create the corresponding FaceElement
	  FACE_ELEMENT<BULK_ELEMENT,DIM>* face_element_pt = new FACE_ELEMENT<BULK_ELEMENT,DIM>
	    (mesh_pt()->boundary_element_pt(b,e),
	     mesh_pt()->face_index_at_boundary(b,e));

	  //Add the face element to the face mesh
	  face_mesh_pt->add_element_pt(face_element_pt);
	}

      // Add the first and last nodes on boundary b to the corner set
      corner_node_set.insert(mesh_pt()->boundary_node_pt(b,0));
      corner_node_set.insert(mesh_pt()->boundary_node_pt(b,n_bound_node-1));

    }

  // Iterate over all nodes in the set and add to the face mesh
  for(it=node_set.begin(); it!=node_set.end(); it++)
    face_mesh_pt->add_node_pt(*it);

  // Iterate over all nodes in the corner set and add to the corner mesh
  for(c_it=corner_node_set.begin(); c_it!=corner_node_set.end(); c_it++)
    corner_mesh_pt->add_node_pt(*c_it);
}


//======start_of_create_global_boundary_equation_number_map===============
/// Create a map from the global equation numbers to the boundary node.
/// Note since we use the global equation number as a key the map will
/// automatically reject duplicate entries.
//========================================================================
template<class BULK_ELEMENT, template<class,unsigned> class FACE_ELEMENT, unsigned DIM>
void TwoDMicromagProblem<BULK_ELEMENT,FACE_ELEMENT,DIM>::
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
void TwoDMicromagProblem<BULK_ELEMENT,FACE_ELEMENT,DIM>::
get_boundary_matrix()
{

  // get the number of nodes in the boundary problem
  unsigned long n_node = face_mesh_pt()->nnode();

  // Initialise and resize the boundary matrix
  Boundary_matrix.resize(n_node,n_node);
  Boundary_matrix.initialise(0.0);

  // Loop over all elements in the face mesh
  unsigned long n_element = face_mesh_pt()->nelement();
  for(unsigned long e=0;e<n_element;e++)
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

  // Loop over all the corner elements and get contributions from the angles
  // at sharp corners.
  for(unsigned long ce=0; ce<Corner_elements.size(); ce++)
    {
      // Find out the index of the node at which the corner element is placed
      unsigned global_j = Corner_elements[ce]->node_pt(0)->eqn_number(0);
      unsigned j = convert_global_to_boundary_equation_number(global_j);

      // Calculate the angle and add to boundary matrix entry (j,j)
      Boundary_matrix(j,j) += Corner_elements[ce]->calculate_corner_fractional_angle();
    }

}

//=============================================================================
///
// updating the boundary conditions on phi_2 from the values
// of phi_1 on the boundary goes in here. This combined with
// including it in the jacobian allows the solver to work as normal.
//=============================================================================
template<class BULK_ELEMENT, template<class,unsigned> class FACE_ELEMENT, unsigned DIM>
void TwoDMicromagProblem<BULK_ELEMENT,FACE_ELEMENT,DIM>::
update_boundary_phi_2()
{
  // Get the index of phi_2
  BULK_ELEMENT* elem_pt = (dynamic_cast<BULK_ELEMENT*>(mesh_pt()->element_pt(0)));
  unsigned phi_2_index = elem_pt->phi_2_index_micromag();

  // Loop over all (target) nodes on the boundary
  unsigned n_boundary_node = face_mesh_pt()->nnode();
  for(unsigned target_node=0; target_node<n_boundary_node; target_node++)
    {
      // Get a pointer to the target node
      Node* target_node_pt = face_mesh_pt()->node_pt(target_node);

      // Get boundary equation number for this target node
      unsigned target_number = convert_global_to_boundary_equation_number
	(target_node_pt->eqn_number(0));

      // Double to store the value of phi_2 during computation
      double target_phi_2_value = 0;

      // Loop over all source nodes adding contribution from each
      for(unsigned source_node=0; source_node<n_boundary_node; source_node++)
	{
	  // Get a pointer to the source node
	  Node* source_node_pt = face_mesh_pt()->node_pt(source_node);

	  // Get boundary equation number for this source node
	  unsigned source_number = convert_global_to_boundary_equation_number
	    (source_node_pt->eqn_number(0));

	  // Add the contribution to phi_2 at the target node due to
	  // the source node (relationship is given by the boundary matrix).
	  //??ds check consistency of boundary matrix numbering
	  target_phi_2_value += Boundary_matrix(target_number,source_number)
	    * source_node_pt->value(0); //??ds replace this with actual phi_1_index
	}
      // Save the total into the target node
      target_node_pt->set_value(phi_2_index,target_phi_2_value);
    }
}

} // End of oomph namespace



//==========start_of_main=================================================
/// ??ds
//========================================================================
int main(int argc, char* argv[])
{
  unsigned n_x, n_y;

  // Get inputs
  if(argc >= 3)
    {
      n_x = atoi(argv[1]);
      n_y = atoi(argv[2]);
    }
  else if(argc == 2)
    {
      n_x = atoi(argv[1]);
      n_y = atoi(argv[1]);
    }
  else
    {
      // Set number of elements in each direction
      n_x = 10;
      n_y = 10;
    }

  // Set dimension of problem and number of nodes along each element edge
  const unsigned dim = 2;
  const unsigned nnode_1d = 2;

  // Set up the bulk problem
  TwoDMicromagProblem<QMicromagElement<dim,nnode_1d>,
    MicromagFaceElement,
    dim>
    problem(n_x,n_y);

  std::cout << "Constructor done." << std::endl;

  // Setup doc info
  DocInfo doc_info;
  doc_info.set_directory("results");

  // // Self tests are REALLLY slow - dominates exectution time, so disable them
  // // for now.
  // // Run self tests on problem and boundary problem:
  // problem.self_test();

  // // // Dump boundary mesh positions
  // unsigned n_node = problem.face_mesh_pt()->nnode();
  // std::cout << n_node << std::endl;
  // // for(unsigned i_node=0; i_node < n_node; i_node++)
  // //   {
  // //     std::cout << boundary_problem.mesh_pt()->node_pt(i_node)->position(0)
  // // 		<< ", " << boundary_problem.mesh_pt()->node_pt(i_node)->position(1)
  // // 		<< std::endl;
  // //   }

  // Set up the boundary equation numbering system
  problem.create_global_boundary_equation_number_map();

  // // ??ds testing my variable gaussian scheme
  QVariableOrderGaussLegendre<dim-1> quadrature_scheme;

  // Set all boundary elements to use the variable order gauss integration
  unsigned n_element = problem.face_mesh_pt()->nelement();
  for(unsigned i=0; i<n_element; i++)
    {
      FiniteElement* finite_element_pt =
	dynamic_cast<FiniteElement*>(problem.face_mesh_pt()->element_pt(i));
      finite_element_pt->set_integration_scheme(&quadrature_scheme);
    }

  // Set very high precision output
  std::cout.precision(16);

  // unsigned max_order = 50;
  // for(unsigned order=2; order<=max_order; order++)
  //   {

  //   // Set the integration scheme order
  //   quadrature_scheme.set_order(order);

  //   // Create + start timer
  //   double start_time = TimingHelpers::timer();

  //   // Get the boundary matrix
  //   problem.get_boundary_matrix();

  //   // output timer result
  //   double stop_time = TimingHelpers::timer();
  //   std::cout << order << " " << stop_time - start_time << std::endl;

  //   // dump for testing
  //   std::ofstream matrix_file;
  //   matrix_file.setf(std::ios::fixed,std::ios::floatfield); // Set high precision output
  //   matrix_file.precision(16);
  //   char filename[100];
  //   sprintf(filename,"results/boundary_matrix_%u",order);
  //   matrix_file.open(filename);
  //   problem.boundary_matrix_pt()->output(matrix_file);
  //   matrix_file.close();
  // }

  // For testing adaptive schemes:
  {

    // Start timer
    double start_time = TimingHelpers::timer();

    // Get the matrix
    problem.get_boundary_matrix();

    // output timer result
    double stop_time = TimingHelpers::timer();
    std::cout << stop_time - start_time << std::endl;

    // dump for testing
    std::ofstream matrix_file;
    matrix_file.precision(16);
    char filename[100];
    sprintf(filename,"results/boundary_matrix_adaptive_nnode1d_%u",nnode_1d);
    matrix_file.open(filename);
    problem.boundary_matrix_pt()->output(matrix_file);
    matrix_file.close();
  }

  return 0;
} //end of main

#endif
