#ifndef OOMPH_BOUNDARY_ELEMENT_HANDLER_H
#define OOMPH_BOUNDARY_ELEMENT_HANDLER_H

/*
  TODO:

  * Parallelisation?

  * some way to check if the mesh has changed?

  * Fix for non-rectangle boundaries!!

  */

#include "../../src/generic/sum_of_matrices.h"

#include "prettyprint98.hpp"
#include "vector_helpers.h"

#include "micromagnetics_boundary_element.h"
#include "micromag_types.h"


namespace oomph
{

  /// Helper function to construct a bem mesh
  inline void build_bem_mesh_helper
  (const Vector<std::pair<unsigned, const Mesh*> >& bem_boundaries,
   BEMElementFactoryFctPt bem_element_factory_fpt,
   Integral* integration_scheme_pt,
   Mesh& bem_mesh)
  {
#ifdef PARANOID
    // Check list of BEM boundaries is not empty
    if(bem_boundaries.size() == 0)
      {
        std::ostringstream error_msg;
        error_msg << "No BEM boundaries are set so there is no need"
                  << " to call build_bem_mesh().";
        throw OomphLibWarning(error_msg.str(),
                              "BoundaryElementHandler::build_bem_mesh",
                              OOMPH_EXCEPTION_LOCATION);
      }


    if(bem_element_factory_fpt == 0)
      {
        throw OomphLibError("BEM factory function is null",
                            OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif

    // Create a set to temporarily store the list of boundary nodes (we use a
    // set because they automatically detect duplicates).
    std::set<Node*> node_set;

    // Loop over entries in bem_boundaries vector.
    for(unsigned i=0; i < bem_boundaries.size(); i++)
      {
        // Get mesh pointer and boundary number from vector.
        const unsigned b = bem_boundaries[i].first;
        const Mesh* mesh_pt = bem_boundaries[i].second;

#ifdef PARANOID
        if(mesh_pt->nboundary_node(b) == 0)
          {
            std::string err = "No nodes on boundary " + to_string(b);
            throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                OOMPH_CURRENT_FUNCTION);
          }
#endif

        // Loop over the nodes on boundary b adding to the set of nodes.
        for(unsigned n=0, nnd=mesh_pt->nboundary_node(b); n<nnd;n++)
          {
            node_set.insert(mesh_pt->boundary_node_pt(b,n));
          }

        // Loop over the elements on boundary b creating bem elements
        for(unsigned e=0, ne=mesh_pt->nboundary_element(b); e<ne;e++)
          {
            int face_index = mesh_pt->face_index_at_boundary(b,e);

            // Create the corresponding BEM Element
            MicromagBEMElementEquations* bem_element_pt =
              bem_element_factory_fpt(mesh_pt->boundary_element_pt(b,e),
                                      face_index);

            // Add the new BEM element to the BEM mesh
            bem_mesh.add_element_pt(bem_element_pt);

            // Set integration pointer
            bem_element_pt->set_integration_scheme(integration_scheme_pt);

            // Set the mesh pointer
            bem_element_pt->set_boundary_mesh_pt(&bem_mesh);
          }
      }

    // Iterate over all nodes in the set and add them to the BEM mesh
    std::set<Node*>::iterator it;
    for(it=node_set.begin(); it!=node_set.end(); it++)
      {
        bem_mesh.add_node_pt(*it);
      }

  }

  // =================================================================
  /// Simple class to store a list of angles associated with nodes of the
  /// boundary element mesh for assembly of the matrix.
  // Assumptions: mesh pointer provided is the boundary element mesh AND
  // the boundary element mesh numbering is being used for BEM lookup in
  // the matrix.
  // =================================================================
  class CornerAngleList
  {
  public:
    /// Default constructor
    CornerAngleList() {}

    /// Destructor
    ~CornerAngleList() {}

    /// Set up corners for a rectangular/cubeoid shape
    void set_up_rectangular_corners(const Mesh* const mesh_pt)
    {
      unsigned nnode = mesh_pt->nnode();
      Corners.assign(nnode,0.0);

#ifdef PARANOID
      if(nnode == 0)
        {
          std::ostringstream error_msg;
          error_msg << "No nodes in mesh.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      // Get dimension (from first node)
      unsigned dim = mesh_pt->node_pt(0)->ndim();

      // The corner_angle is assumed to be 1/(2^dim), since it is rectangular or
      // cubeoid.
      double corner_angle = 1 / std::pow(2.0,dim);

      // Loop over nodes
      for(unsigned nd=0; nd<nnode; nd++)
        {
          Node* nd_pt = mesh_pt->node_pt(nd);

          // Get list of boundaries that this node is on
          std::set<unsigned>* boundaries_pt;
          nd_pt->get_boundaries_pt(boundaries_pt);

#ifdef PARANOID
          if(boundaries_pt == 0)
            {
              std::ostringstream error_msg;
              error_msg << "Node is not on any boundaries, this probably "
                        << "means something has gone wrong, maybe you passed "
                        << "in the bulk mesh?";
              throw OomphLibError(error_msg.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
#endif

          // If it is on dim many boundaries then this is a corner, otherwise
          // assume it is a smooth point and so the angle is 0.5
          if(boundaries_pt->size() == dim)
            {
              Corners[nd] = corner_angle;

              std::cout << "Think I've found a corner at [";
              for(unsigned j=0; j<nd_pt->ndim(); j++)
                {
                  std::cout << nd_pt->x(j) << ", ";
                }
              std::cout << "] I gave it the angle " << corner_angle << std::endl;
            }
          else if((dim == 3) && (boundaries_pt->size() == (dim - 1)))
            {
              // in 3d we also have edges with solid angle pi/4pi
              Corners[nd] = 0.25;

              std::cout << "Think I've found an edge at [";
              for(unsigned j=0; j<nd_pt->ndim(); j++)
                {
                  std::cout << nd_pt->x(j) << ", ";
                }
              std::cout << "] I gave it the angle " << Corners[nd] << std::endl;
            }
          else
            {
              // Angle = pi/2pi or 2pi/4pi in 2 or 3 dimensions respectively.
              Corners[nd] = 0.5;
            }
        }

    }

    /// Set up corners for a smooth mesh (i.e. no sharp corners).
    void set_up_smooth_mesh(const Mesh* const mesh_pt)
    {
      Corners.assign(mesh_pt->nnode(),0.5);
    }

    /// Add the contribution due to corners to the diagonal of the boundary
    /// matrix.
    void add_corner_contributions(DoubleMatrixBase& bem_matrix) const
    {
      // Assume that the bem matrix is a densedoublematrix so that we can write
      // to it with operator().

      //??ds fix for h matrix
      DenseDoubleMatrix* bem_matrix_pt =
        checked_dynamic_cast<DenseDoubleMatrix*>(&bem_matrix);


#ifdef PARANOID
      // Check that the list has been set up
      if(!(is_set_up()))
        {
          std::ostringstream error_msg;
          error_msg << "Corner list has not been set up.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      // Check that it is the correct size
      if(bem_matrix_pt->nrow() != Corners.size())
        {
          std::ostringstream error_msg;
          error_msg << "Corners list is the wrong size for the matrix.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      // Add the fractional angles
      for(unsigned nd=0, s=Corners.size(); nd<s; nd++)
        {
          bem_matrix_pt->operator()(nd,nd) += Corners[nd];
        }
    }

    /// Create diagonal matrix just containing the corner contributions.
    void make_diagonal_corner_matrix(CRDoubleMatrix& corner_matrix) const
    {
      VectorOps::diag_cr_matrix(corner_matrix, Corners);
    }


    void set_up_from_input_data
    (const Mesh* const mesh_pt,
     const Vector<std::pair<Vector<double>, double> >* const input_data_pt)
    {
      // Initialise to default values
      Corners.assign(mesh_pt->nnode(),0.5);

      // Look through input list of corner locations + angles, find the
      // corners and add to our list.
      for(unsigned i=0; i<input_data_pt->size(); i++)
        {
          unsigned bem_node_number =
            find_node_by_position_in_mesh(mesh_pt, (*input_data_pt)[i].first);
          Corners[bem_node_number] = (*input_data_pt)[i].second;
        }
    }

    /// Check if the list has been set up.
    bool is_set_up() const
    {
      // If there is a non-zero length vector something has been set up.
      return (Corners.size() != 0);
    }

  private:

    unsigned find_node_by_position_in_mesh(const Mesh* mesh_pt,
                                           const Vector<double> &x) const
    {
      for(unsigned nd=0, nnode=mesh_pt->nnode(); nd<nnode; nd++)
        {
          Node* nd_pt = mesh_pt->node_pt(nd);

          Vector<double> node_x(x.size(), 0.0);
          nd_pt->position(node_x);

          if( VectorOps::numerically_close(x, node_x) )
            {
              return nd;
            }
        }

      std::ostringstream error_msg;
      error_msg << "Failed to find node at " << x;
      throw OomphLibError(error_msg.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
      return 0;
    }


    /// Storage for the location and angle of the corners.
    Vector<double> Corners;

    /// Inaccessible copy constructor
    CornerAngleList(const CornerAngleList& dummy)
    {BrokenCopy::broken_copy("CornerAngleList");}

    /// Inaccessible assignment operator
    void operator=(const CornerAngleList& dummy)
    {BrokenCopy::broken_assign("CornerAngleList");}
  };


  // =================================================================
  /// A class implementing all the boundary element methods stuff needed for
  /// the hybrid method in micromagnetics. Problems can then simply contain
  /// an instance of this class to gain access to hybrid method functions
  /// (see the "Prefer composition over inheritance" principle).
  ///
  /// In particular this allows the semi-implict "problem" class to not
  /// inherit from problem. This is good because many functions defined in
  /// Problem (for fully implicit methods) do not make sense for the
  /// semi-implicit "problem".
  // =================================================================
  class BoundaryElementHandler
  {

    // Note: in BEM input (output) index == jacobian col (row) == phi_1 (phi)
    // global equation numbers respectively.

  public:

    /// Default constructor
    BoundaryElementHandler() :
      Bem_element_factory_fpt(0), Input_index(0), Output_index(0)
    {
      // Make an empty mesh for the bem nodes and elements
      Bem_mesh_pt = new Mesh;

      // By default evaluate BEM integrals analytically.
      Numerical_int_bem = false;

      // By default use hlib if we have it
#ifdef OOMPH_HAS_HLIB
      Hierarchical_bem = true;
#else
      Hierarchical_bem = false;
#endif

      // Debugging flags
      Debug_disable_corner_contributions=false;

      // Null pointers
      Bem_matrix_pt = 0;
      Integration_scheme_pt = 0;
      Hmatrix_dof2idx_mapping_pt = 0;
    }

    /// Destructor
    ~BoundaryElementHandler()
    {
      // Delete bem mesh, but first clear out the node pointers so that it
      // doesn't delete them (the node pointers belong to the bulk mesh).
      Bem_mesh_pt->flush_node_storage();
      delete Bem_mesh_pt; Bem_mesh_pt = 0;

      // Delete the integrator
      delete Integration_scheme_pt;
      Integration_scheme_pt = 0;

      // Delete the matrix
      delete Bem_matrix_pt;
      Bem_matrix_pt = 0;

      delete Hmatrix_dof2idx_mapping_pt;
      Hmatrix_dof2idx_mapping_pt = 0;
    }

    /// Put the (output) values of the bem into a vector.
    void get_bem_values(DoubleVector &bem_output_values) const;

    /// Put the (output) values of the bem into a series of vectors (one
    /// per boundary).
    void get_bem_values(const Vector<DoubleVector*> &bem_output_values) const;

    /// Put the output values of the bem directly into the values of the
    /// boundary nodes.
    void get_bem_values_and_copy_into_values() const;

    /// Build the mesh, lookup schemes and matrix in that order.
    void build(const CornerDataInput& input_corner_data)
    {
      // Force use of numerical integration if needed (if nodal dimension
      // of meshes is not 3).
      if((!Numerical_int_bem) && (dimension() != 3))
        {
          OomphLibWarning("Problem appears to not be 3 dimensional. I'm forcing the use of numerical integration for boundary element calculations (analytical calculations only work for 3d).",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
          Numerical_int_bem = true;
        }

      // Can't force numerical integration in hlib
      if(Hierarchical_bem && Numerical_int_bem)
        {
          std::string err = "Can't control integration method in hlib, so we";
          err += "can't use numerical integration.";
          OomphLibWarning(err, OOMPH_EXCEPTION_LOCATION,
                          OOMPH_CURRENT_FUNCTION);
        }

      // Construct the mesh using on boundaries specified in Bem_boundaries
      build_bem_mesh();

#ifdef PARANOID
      // Try to check if equation numbering has been set up. If it hasn't
      // then all equation numbers are -10 (Data::Is_unclassified). Just
      // check the first one, don't think it can legally be -10 at this
      // point.
      if(bem_mesh_pt()->node_pt(0)->eqn_number(0) == Data::Is_unclassified)
        {
          std::string err = "An equation number is unclassified, you probably";
          err += " haven't set up the equation numbering yet.";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif

      // Construct the lookup schemes
      Input_lookup.build(bem_mesh_pt(), input_index());
      Output_lookup.build(bem_mesh_pt(), output_index());

      // Set up the corner angle data
      if(input_corner_data.size() > 0)
        {
          corner_list_pt()->set_up_from_input_data(bem_mesh_pt(),
                                                   &input_corner_data);
        }
      // If none has been provided then assume rectangular mesh ??ds
      // generalise?
      else
        {
          oomph_info << "No input corner list, assuming rectangular..."
                     << std::endl;
          corner_list_pt()->set_up_rectangular_corners(bem_mesh_pt());
        }

      if(!Hierarchical_bem)
        {
          // Say what we're doing
          oomph_info << "Building dense BEM matrix, this may take some time"
                     << std::endl;
          if(Numerical_int_bem)
            oomph_info << "Using numerical integration." << std::endl;
          else
            oomph_info << "Using analytical integration." << std::endl;

          // Construct the (dense) matrix
          build_bem_matrix();
        }
      else
        {
          oomph_info << "Building hierarchical BEM matrix" << std::endl;
          build_hierarchical_bem_matrix();
        }


    }

    /// \short Get the (nodal) dimension of the boundary meshes
    unsigned dimension()
    {
      // Use mesh pointer of the first boundary presumably all boundaries
      // must have the same nodal dimension.
      return Bem_boundaries[0].second->node_pt(0)->ndim();
    }

    /// If we are using H-matrix for bem then write out some data on
    /// it. Otherwise do nothing. (Wrapped up in a function so that
    /// llg_problem.h doesn't need to know anything about hlib, HMatrix
    /// etc.)
    void maybe_write_h_matrix_data(const std::string& outdir) const;


    /// Get the (internal, bem only) output lookup's equation number of a
    /// node.
    unsigned output_equation_number(const Node* node_pt) const
    {
      unsigned g_eqn = node_pt->eqn_number(output_index());
      return output_lookup_pt()->main_to_added(g_eqn);
    }

    /// Get the (internal, bem only) input lookup's equation number of a
    /// node.
    unsigned input_equation_number(const Node* node_pt) const
    {
      unsigned g_eqn = node_pt->eqn_number(input_index());
      return input_lookup_pt()->main_to_added(g_eqn);
    }


    // Access functions:
    // ============================================================

    /// Access to the pointer to the boundary element method mesh.
    const Mesh* bem_mesh_pt() const {return Bem_mesh_pt;}

    /// Const access to the boundary matrix
    const DoubleMatrixBase* bem_matrix_pt() const
    {return Bem_matrix_pt;}

    /// Non-const access to the boundary matrix
    DoubleMatrixBase* bem_matrix_pt() {return Bem_matrix_pt;}

    /// \short Set function for Input_index.
    void set_input_index(unsigned input_index) {Input_index = input_index;}

    /// \short Const access function for Input_index.
    unsigned input_index() const {return Input_index;}

    /// \short Set function for Output_index.
    void set_output_index(unsigned output_index) {Output_index = output_index;}

    /// \short Const access function for Output_index.
    unsigned output_index() const {return Output_index;}

    /// \short Const access function for a pointer to the map of sharp corner
    /// angles at nodes.
    const CornerAngleList* corner_list_pt() const
    {return &Corner_list;}

    /// \short Non-const access function for a pointer to the map of sharp corner
    /// angles at nodes.
    CornerAngleList* corner_list_pt() {return &Corner_list;}

    /// \short Pointer to a lookup between output value global equation
    /// numbers and node numbers within mesh.
    const AddedMainNumberingLookup* output_lookup_pt() const
    {return &Output_lookup;}

    /// \short Alias of output_lookup_pt for when we are working with Jacobians
    /// (they are the same lookup).
    const AddedMainNumberingLookup* row_lookup_pt() const
    {return output_lookup_pt();}

    /// \short Pointer to a lookup between input value global equation
    /// numbers and node numbers within mesh.
    const AddedMainNumberingLookup* input_lookup_pt() const
    {return &Input_lookup;}

    /// \short Alias of input_lookup_pt for when we are working with Jacobians
    /// (they are the same lookup).
    const AddedMainNumberingLookup* col_lookup_pt() const
    {return input_lookup_pt();}

    /// \short Get an appropriate linear algebra distribution for working
    /// with the boundary matrix.
    void get_bm_distribution(LinearAlgebraDistribution& dist) const;

    /// \short Non-const access function for Integration_scheme.
    Integral* &integration_scheme_pt() {return Integration_scheme_pt;}

    /// \short Const access function for Integration_scheme.
    const Integral* integration_scheme_pt() const {return Integration_scheme_pt;}

    MicromagBEMElementEquations* bem_element_factory(FiniteElement* const fe_pt,
                                                     const int& face) const
    {
#ifdef PARANOID
      if(Bem_element_factory_fpt == 0)
        {
          throw OomphLibError("No BEM factory function has been specified! Can't create BEM elements.",
                              OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif
      return Bem_element_factory_fpt(fe_pt, face);
    }

    BEMElementFactoryFctPt Bem_element_factory_fpt;

    /// \short Integrate BEM integrals by adaptive numerical integration or
    /// analytical integration.
    bool Numerical_int_bem;

    /// Use a hierarchical representation of the BEM matrix. Much faster
    /// but requires external library.
    bool Hierarchical_bem;

    /// Mapping from Hmatrix dofs (equivalent to the nodal numbering
    /// scheme) to Hmatrix indices (internal Hmatrix ordering used to
    /// improve hierarchical properties). Null pointer unless we are using
    /// Hmatrix bem.
    AddedMainNumberingLookup* Hmatrix_dof2idx_mapping_pt;

    /// A list of the boundaries (on various meshes) to which the boundary
    /// element method should be applied. ??ds will multiple meshes work?
    /// are they needed? probably not for anything I do...
    BemBoundaryData Bem_boundaries;

    /// Debugging flag: don't add sharp corner solid angles to bem matrix.
    bool Debug_disable_corner_contributions;

  private:

    /// \short Lookup between output value's global equation numbers and
    /// node numbers within mesh.
    AddedMainNumberingLookup Output_lookup;

    /// \short Lookup between input value's global equation numbers and
    /// node numbers within mesh.
    AddedMainNumberingLookup Input_lookup;

    /// \short Storage for the adaptive integration scheme to be used.
    Integral* Integration_scheme_pt;

    /// The pointer to the "boundary element" mesh (as in boundary element method
    /// not finite elements on the boundary).
    Mesh* Bem_mesh_pt;

    /// Pointer to storage for the list of nodal angles/solid angles.
    CornerAngleList Corner_list;

    /// The (local/elemental) index of the dof we take input values from.
    unsigned Input_index;

    /// The (local/elemental) index of the dof we are determining the
    /// boundary conditions for.
    unsigned Output_index;

    /// Matrix to store the relationship between phi_1 and phi on the boundary
    DoubleMatrixBase* Bem_matrix_pt;

    /// Construct BEM elements on boundaries listed in Bem_boundaries and add
    /// to the Bem_mesh.
    void build_bem_mesh();

    /// Construct a dense boundary matrix in Bem_matrix using pure
    /// oomph-lib code.
    void build_bem_matrix();

    /// Construct a hierarchical boundary matrix in Bem_matrix using hlib.
    void build_hierarchical_bem_matrix();

    /// \short Get the mapping between the global equation numbering and
    /// the boundary equation numbering.
    void create_global_boundary_equation_number_maps();

    /// Inaccessible copy constructor
    BoundaryElementHandler(const BoundaryElementHandler& dummy)
    {BrokenCopy::broken_copy("BoundaryElementHandler");}

    /// Inaccessible assignment operator
    void operator=(const BoundaryElementHandler& dummy)
    {BrokenCopy::broken_assign("BoundaryElementHandler");}

  };


} // End of oomph namespace

#endif
