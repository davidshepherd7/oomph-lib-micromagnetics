#ifndef OOMPH_PINNED_BOUNDARY_ELEMENT_HANDLER_H
#define OOMPH_PINNED_BOUNDARY_ELEMENT_HANDLER_H

/*
  TODO:

  * Parallelisation?

  * some way to check if the mesh has changed?

  * Fix for non-rectangle boundaries!!

  */

#include "boundary_element_handler.h"

#include "../../src/generic/sum_of_matrices.h"

#include "vector_helpers.h"

#include "micromagnetics_boundary_element.h"
#include "micromag_types.h"


namespace oomph
{

  class BemLookup
  {
  public:
    /// Constructor
    BemLookup() {}

    /// Virtual destructor
    virtual ~BemLookup() {}

    void build(const Mesh* mesh_pt,
               bool use_pinned, bool use_unpinned, const unsigned& dof_index)
    {
      Dof_index = dof_index;

      lookup.clear();

      const unsigned n_node = mesh_pt->nnode();
      for(unsigned nd=0; nd<n_node; nd++)
        {
          Node* nd_pt = mesh_pt->node_pt(nd);
          bool is_pinned = nd_pt->is_pinned(dof_index);
          if((use_pinned && is_pinned) || (use_unpinned && !is_pinned))
            {
              lookup.push_back(nd_pt);
            }
        }
    }

    unsigned node_to_bemeq(const Node* node_pt) const
    {
      std::vector<const Node*>::const_iterator it
        = std::find(lookup.begin(), lookup.end(), node_pt);
#ifdef PARANOID
      if(it == lookup.end())
        {
          std::string err = "Node not found";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return (it - lookup.begin());
    }

    const Node* bemeq_to_node(const unsigned& bemeq) const
    {
      return lookup[bemeq];
    }

    unsigned bemeq_to_global(const unsigned& bemeq) const
    {
      return bemeq_to_node(bemeq)->eqn_number(Dof_index);
    }

    unsigned size() const {return lookup.size();}

    Vector<const Node*> lookup;

    unsigned Dof_index;

  private:
    /// Broken copy constructor
    BemLookup(const BemLookup& dummy)
    {BrokenCopy::broken_copy("BemLookup");}

    /// Broken assignment operator
    void operator=(const BemLookup& dummy)
    {BrokenCopy::broken_assign("BemLookup");}

  };

  // =================================================================
  /// Simple class to store a list of angles associated with nodes of the
  /// boundary element mesh for assembly of the matrix.
  // Assumptions: mesh pointer provided is the boundary element mesh AND
  // the boundary element mesh numbering is being used for BEM lookup in
  // the matrix.
  // =================================================================
  class PinnedBoundaryCornerAngleList
  {
  public:
    /// Default constructor
    PinnedBoundaryCornerAngleList() {}

    /// Destructor
    ~PinnedBoundaryCornerAngleList() {}

    /// Set up corners for a rectangular/cubeoid shape
    void set_up_rectangular_corners(const Mesh* const mesh_pt)
    {
      unsigned nnode = mesh_pt->nnode();

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
              Corners.insert(std::make_pair(nd_pt, corner_angle));

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
              Corners.insert(std::make_pair(nd_pt, 0.25));

              std::cout << "Think I've found an edge at [";
              for(unsigned j=0; j<nd_pt->ndim(); j++)
                {
                  std::cout << nd_pt->x(j) << ", ";
                }
              std::cout << "] I gave it the angle " << 0.25 << std::endl;
            }
          else
            {
              // Angle = pi/2pi or 2pi/4pi in 2 or 3 dimensions respectively.
              Corners.insert(std::make_pair(nd_pt, 0.5));
            }
        }

    }

    /// Set up corners for a smooth mesh (i.e. no sharp corners).
    void set_up_smooth_mesh(const Mesh* const mesh_pt)
    {
      const unsigned n_node = mesh_pt->nnode();
      for(unsigned nd=0; nd<n_node; nd++)
        {
          Node* nd_pt = mesh_pt->node_pt(nd);
          Corners.insert(std::make_pair(nd_pt, 0.5));
        }
    }

    /// Add the contribution due to corners to the diagonal of the boundary
    /// matrix.
    void add_corner_contributions(DoubleMatrixBase& bem_matrix,
                                  const BemLookup& lookup_unpinned_input,
                                  const BemLookup& lookup_unpinned_output,
                                  const unsigned& dof_index) const
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
          error_msg << "Corners list is the wrong size for the matrix rows";
          error_msg << "\n bem matrix nrow: " << bem_matrix_pt->nrow();
          error_msg << "\ncorners list size " << Corners.size();
          error_msg << "\n should all be the same";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      if((bem_matrix_pt->ncol() != Corners.size())
         && (bem_matrix_pt->ncol() != Corners.size() -1))
        {
          std::string err = "Corners list is the wrong size for the matrix cols";
          err += "ncols are " + to_string(bem_matrix_pt->ncol());
          err += "corner list size is " + to_string(Corners.size());
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      // Add the fractional angles to the appropriate places
      std::map<const Node*, double>::const_iterator it;
      unsigned pinned_correction = 0;
      for(it = Corners.begin(); it != Corners.end(); ++it)
        {
          if(!it->first->is_pinned(dof_index))
            {
              const unsigned n = lookup_unpinned_input.node_to_bemeq(it->first);
              const unsigned m = lookup_unpinned_output.node_to_bemeq(it->first);
              bem_matrix_pt->operator()(m, n) += it->second;
            }
          else
            {
              pinned_correction++;
#ifdef PARANOID
              if(pinned_correction > 1)
                {
                  std::string err = "Too many pinned phi1 values!";
                  throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
                }
#endif
            }
        }
    }

    /// Create diagonal matrix just containing the corner contributions.
    void make_diagonal_corner_matrix(CRDoubleMatrix& corner_matrix) const
    {
      Vector<double> vcorners(Corners.size());

      throw OomphLibError("Not implemented (yet?).", OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);

      VectorOps::diag_cr_matrix(corner_matrix, vcorners);
    }


    void set_up_from_input_data
    (const Mesh* const mesh_pt,
     const Vector<std::pair<Vector<double>, double> >* const input_data_pt)
    {
      // Initialise to default values
      const unsigned n_node = mesh_pt->nnode();
      for(unsigned nd=0; nd<n_node; nd++)
        {
          Node* nd_pt = mesh_pt->node_pt(nd);
          Corners.insert(std::make_pair(nd_pt, 0.5));
        }

      // Look through input list of corner locations + angles, find the
      // corners and add to our list.
      for(unsigned i=0; i<input_data_pt->size(); i++)
        {
          Node* nd_pt =
            find_node_by_position_in_mesh(mesh_pt, (*input_data_pt)[i].first);
          Corners[nd_pt] = (*input_data_pt)[i].second;
        }
    }

    /// Check if the list has been set up.
    bool is_set_up() const
    {
      // If there is a non-zero length map something has been set up.
      return (Corners.size() != 0);
    }

  private:

    Node* find_node_by_position_in_mesh(const Mesh* mesh_pt,
                                        const Vector<double> &x) const
    {
      for(unsigned nd=0, nnode=mesh_pt->nnode(); nd<nnode; nd++)
        {
          Node* nd_pt = mesh_pt->node_pt(nd);

          Vector<double> node_x(x.size(), 0.0);
          nd_pt->position(node_x);

          if( VectorOps::numerically_close(x, node_x) )
            {
              return nd_pt;
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
    // Vector<double> Corners;
    std::map<const Node*, double> Corners;

    /// Inaccessible copy constructor
    PinnedBoundaryCornerAngleList(const PinnedBoundaryCornerAngleList& dummy)
    {BrokenCopy::broken_copy("PinnedBoundaryCornerAngleList");}

    /// Inaccessible assignment operator
    void operator=(const PinnedBoundaryCornerAngleList& dummy)
    {BrokenCopy::broken_assign("PinnedBoundaryCornerAngleList");}
  };


  /// A class implementing all the boundary element methods stuff needed for
  /// the hybrid method in micromagnetics. Problems can then simply contain
  /// an instance of this class to gain access to hybrid method functions
  /// (see the "Prefer composition over inheritance" principle).
  class PinnedBoundaryElementHandler : public BoundaryElementHandlerBase
  {

    // Note: in BEM input (output) index == jacobian col (row) == phi_1 (phi)
    // global equation numbers respectively.

  public:

    /// Default constructor
    PinnedBoundaryElementHandler() {}

    /// Destructor
    ~PinnedBoundaryElementHandler() {}

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
      Lookup_all_nodes.build(bem_mesh_pt(), true, true, 0);
      Lookup_unpinned_input.build(bem_mesh_pt(), false, true, input_index());
      Lookup_unpinned_output.build(bem_mesh_pt(), false, true, output_index());

      Lookup_pinned_input.build(bem_mesh_pt(), true, false, input_index());
      Lookup_pinned_output.build(bem_mesh_pt(), true, false, output_index());

      Input_lookup.build(Lookup_unpinned_input.lookup, Lookup_unpinned_input.Dof_index);
      Output_lookup.build(Lookup_unpinned_output.lookup, Lookup_unpinned_output.Dof_index);

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


    // Access functions:
    // ============================================================

    /// \short Const access function for a pointer to the map of sharp corner
    /// angles at nodes.
    const PinnedBoundaryCornerAngleList* corner_list_pt() const
    {return &Corner_list;}

    /// \short Non-const access function for a pointer to the map of sharp corner
    /// angles at nodes.
    PinnedBoundaryCornerAngleList* corner_list_pt() {return &Corner_list;}


    // Lookup schemes
    // ============================================================
    BemLookup Lookup_all_nodes;
    BemLookup Lookup_unpinned_input;
    BemLookup Lookup_unpinned_output;
    BemLookup Lookup_pinned_input;
    BemLookup Lookup_pinned_output;


  private:

    /// Pointer to storage for the list of nodal angles/solid angles.
    PinnedBoundaryCornerAngleList Corner_list;

    /// Construct a dense boundary matrix in Bem_matrix using pure
    /// oomph-lib code.
    void build_bem_matrix();

    /// Construct a hierarchical boundary matrix in Bem_matrix using hlib.
    void build_hierarchical_bem_matrix();

    /// Inaccessible copy constructor
    PinnedBoundaryElementHandler(const PinnedBoundaryElementHandler& dummy)
    {BrokenCopy::broken_copy("PinnedBoundaryElementHandler");}

    /// Inaccessible assignment operator
    void operator=(const PinnedBoundaryElementHandler& dummy)
    {BrokenCopy::broken_assign("PinnedBoundaryElementHandler");}

  };


} // End of oomph namespace

#endif
