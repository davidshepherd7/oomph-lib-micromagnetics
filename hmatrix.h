#ifndef OOMPH_H_MATRIX_BEM3D_H
#define OOMPH_H_MATRIX_BEM3D_H


// Only if we have hlib installed, otherwise we create a dummy class, see
// below.
#ifdef OOMPH_HAS_HLIB


#include <set>
#include <map>
#include <utility> // std::pair

// hlib code
#include "bem3d.h"
#include "cluster.h"
#include "laplacebem.h"
#include "supermatrix.h"
#include "hcoarsening.h"

// oomph-lib core code
#include "../../src/generic/Vector.h"
#include "../../src/generic/sum_of_matrices.h"
#include "../../src/generic/nodes.h"
#include "../../src/generic/matrices.h"

// micromag code
#include "boundary_element_handler.h"





// TODO:

// handle 2d bemgrids as well? Not so easy because would have to implement
// integration etc. myself, can't use Lindholm...

// Figure out parameters and provide access functions!


namespace hlib_helpers
{
  /// Helper function: Given a mesh create a lookup from node pointer to
  /// node numbers in the mesh. ??ds still not entirely sure how this
  /// interacts with my bem equation numbering. Might break on multiple
  /// processors when the equation numbers are not just the number of the
  /// node in the mesh... but not handling mpi yet anyway!
  inline std::map<Node*, unsigned> build_nodept2node
  (const Mesh& mesh, const BoundaryElementHandler& bem_handler,
   bool use_output=true)
  {
    std::map<Node*, unsigned> lookup;

    unsigned index;
    if(use_output)
      {
        index = bem_handler.output_index();
      }
    else
      {
        index = bem_handler.input_index();
      }

    for(unsigned nd=0, nnd=mesh.nnode(); nd<nnd; nd++)
      {
        Node* nd_pt = mesh.node_pt(nd);

        // Get bem-local-equation-number (i.e. the index of this node's
        // dof in the bem matrix).
        int g_eq = nd_pt->eqn_number(index);
        unsigned eq;
        if(use_output)
          {
            eq = bem_handler.output_lookup_pt()->main_to_added(g_eq);
          }
        else
          {
            eq = bem_handler.input_lookup_pt()->main_to_added(g_eq);
          }

        lookup[nd_pt] = eq;
      }

    return lookup;
  }


  /// Make a pair of integers which always has the first entry as the
  /// smaller of the two. Needed because pair(a, b) != pair(b, a)
  /// otherwise, so we can't use std::set to easily get a unique list.
  inline std::pair<int, int> ordered_make_pair(int a, int b)
  {
    if(a < b)
      {
        return std::make_pair(a, b);
      }
    else
      {
        return std::make_pair(b, a);
      }
  }

}


// A matrix class for hierarchical matrices for 3D BEM on triangular
// elements. Could probably be extended to work with other types of
// H-matrix fairly easily.
class HMatrix : public DoubleMatrixBase
{
public:

  HMatrix()
  {
#ifdef OOMPH_HAS_MPI
    throw OomphLibError("Not implemented for MPI",
                        OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
#endif

    // Parameters taken from Andreas Knittel's thesis
    nmin = 30;
    eps_aca = 1e-5;
    kmax = 500;
    eps = 1e-4;
    algorithm = HLIB_HCAII;
    eta = 0.25;
    quadorder = 3;
    polyorder = 4;

    // Null pointers
    This_hmatrix_pt = 0;
    Cluster_tree_pt = 0;
  }

  virtual ~HMatrix()
  {
    // Delete data in hlib data structures and null their pointers. Because
    // its C-code their delete functions aren't necessarily required to
    // handle null pointers intelligently, and they don't :( So we need to
    // wrap it all in ifs...
    if(This_hmatrix_pt != 0)
      {
        del_supermatrix(This_hmatrix_pt);
        This_hmatrix_pt = 0;
      }

    if(Cluster_tree_pt != 0)
      {
        del_clustertree(Cluster_tree_pt);
        Cluster_tree_pt = 0;
      }
  }


  /// Create an equivalent dense matrix.
  void to_dense(DenseDoubleMatrix& out) const
  {
    // Create a full (dense) hlib matrix in a temporary data structure
    fullmatrix* full_pt = new_fullmatrix(This_hmatrix_pt->rows,
                                         This_hmatrix_pt->cols);
    convertsuper2_fullmatrix(full_pt, This_hmatrix_pt);

    // Just copy the data over (might be possible to do something fancier
    // if better performance is needed, maybe just assign full_pt->
    // somewhere?).
    const unsigned cols = full_pt->cols, rows = full_pt->rows;
    const double* elements = full_pt->e;
    out.resize(rows, cols);
    for(unsigned i=0; i<rows; i++)
      {
        const unsigned ii = Cluster_tree_pt->idx2dof[i];
        for(unsigned j=0; j<cols; j++)
          {
            const unsigned jj = Cluster_tree_pt->idx2dof[j];
            out(i, j) = elements[ii + jj*rows];
          }
      }

    // Clean up the temporary data
    del_fullmatrix(full_pt);
  }

  /// Build the hmatrix for a given mesh and bem handler (needed for
  /// equation lookup at the moment...).
  void build(const Mesh& mesh,
             const BoundaryElementHandlerBase* _bem_handler_pt)
  {
    // Only works with normal bem handler (for now?)
    const BoundaryElementHandler* bem_handler_pt =
      dynamic_cast<const BoundaryElementHandler*>(_bem_handler_pt);

#ifdef PARANOID
    if(bem_handler_pt == 0)
      {
        std::string err = "Not implemented for bem handlers other than BoundaryElementHandler";
        throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

#ifdef PARANOID
    {
      FiniteElement* fe_pt = mesh.finite_element_pt(0);
      if(fe_pt->nnode() != 3 || fe_pt->dim() != 2 || fe_pt->nodal_dimension() != 3)
        {
          std::string err = "Only works for 2D triangles in 3D space ";
          err += "you appear to have " + to_string(fe_pt->dim())
            + "D elements with " + to_string(fe_pt->nnode())
            + " nodes in " + to_string(fe_pt->nodal_dimension())
            +"D space.";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
          // Implementing anything else will probably be hard, you would need
          // to write integration routines in hlib itself.
        }
    }

    if(This_hmatrix_pt != 0)
      {
        std::string err = "Already built a hmatrix in here, replacement not yet implemented";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
        // To implement this just make sure you delete all the data
        // structures (only if pointers non-null) and null their pointers.
      }
#endif

    // Get some helper functions
    using namespace hlib_helpers;


    // Copy the oomph mesh data into hlib data structures
    // ============================================================

    // Create a reverse lookup scheme which finds the node number when
    // given a node pointer.
    std::map<Node*, unsigned> nodept2node =
      build_nodept2node(mesh, *bem_handler_pt, false);

    // Get a unique list (using std::set) of all the edges in the mesh
    std::set<std::pair<int, int> > edges_set;
    for(unsigned ele=0, nele=mesh.nelement(); ele<nele; ele++)
      {
        FiniteElement* ele_pt = mesh.finite_element_pt(ele);

        unsigned nd0 = nodept2node.find(ele_pt->node_pt(0))->second;
        unsigned nd1 = nodept2node.find(ele_pt->node_pt(1))->second;
        unsigned nd2 = nodept2node.find(ele_pt->node_pt(2))->second;

        // Insert all pairs of the nodes into the set (function makes
        // sure they are in order inside the pair so that edge {0, 1} ==
        // edge {1, 0}).
        edges_set.insert(ordered_make_pair(nd0, nd1));
        edges_set.insert(ordered_make_pair(nd1, nd2));
        edges_set.insert(ordered_make_pair(nd2, nd0));
      }

    unsigned vertices = mesh.nnode();
    unsigned triangles = mesh.nelement();
    unsigned edges = edges_set.size();

    bemgrid3d* bem_grid_pt = new_bemgrid3d(vertices, edges, triangles);

    // Copy vertex data (nodal positions) into the hlib grid
    for(unsigned nd=0, nnd=mesh.nnode(); nd<nnd; nd++)
      {
        Node* nd_pt = mesh.node_pt(nd);
        bem_grid_pt->x[nd][0] = nd_pt->x(0);
        bem_grid_pt->x[nd][1] = nd_pt->x(1);
        bem_grid_pt->x[nd][2] = nd_pt->x(2);
      }

    // Copy data on which vertices make up each triangle (element) into
    // hlib grid.
    for(unsigned ele=0, nele=mesh.nelement(); ele<nele; ele++)
      {
        FaceElement* ele_pt = dynamic_cast<FaceElement*>(mesh.element_pt(ele));

        for(unsigned nd=0, nnd=ele_pt->nnode(); nd<nnd; nd++)
          {
            Node* nd_pt = ele_pt->node_pt(nd);
            bem_grid_pt->t[ele][nd] = nodept2node.find(nd_pt)->second;
          }

        // Order nodes such that Lindholm formula will get the right sign
        // of the unit normal. If the oomph-lib normal_sign() is negative
        // then the nodes are in the wrong order to get the outer unit
        // normal via cross products, so swap two of them.
        if(ele_pt->normal_sign() < 0)
          {
            std::swap(bem_grid_pt->t[ele][0],
                      bem_grid_pt->t[ele][1]);
          }
      }

    // copy to array
    std::set<std::pair<int, int> >::const_iterator it;
    unsigned i_edge = 0;
    for(it = edges_set.begin(); it != edges_set.end(); ++it)
      {
        bem_grid_pt->e[i_edge][0] = it->first;
        bem_grid_pt->e[i_edge][1] = it->second;
        i_edge++;
      }

    // Fill in some required geometrical values
    prepare_bemgrid3d(bem_grid_pt);


    // Construct the matrix itself, via some other intermediate data
    // structures.
    // ============================================================

    // Build the cluster tree
    Cluster_tree_pt = buildvertexcluster_bemgrid3d(bem_grid_pt,
                                                   HLIB_REGULAR,
                                                   nmin, 0);

    // Build the factory which will finally create the hmatrix
    surfacebemfactory* surface_bem_factory_pt = new_surfacebemfactory_dlp_collocation
      (bem_grid_pt, HLIB_LINEAR_BASIS, Cluster_tree_pt,
       bem_grid_pt, HLIB_LINEAR_BASIS, Cluster_tree_pt,
       quadorder, quadorder,
       polyorder, 0.0);


    // Finally build the matrix itself, with adaptive recompression
    This_hmatrix_pt = onthefly_hca_coarsen_supermatrix
      (Cluster_tree_pt->root, Cluster_tree_pt->root,
       surface_bem_factory_pt,
       eps_aca, kmax, eps, 1, algorithm, 0, eta, 0);


    // Clean up auxilary data structures
    // ============================================================

    // Since theses are C data types they don't have destructors to do this
    // for us, so do it manually now:
    del_bemgrid3d(bem_grid_pt);
    bem_grid_pt = 0;
    del_surfacebemfactory(surface_bem_factory_pt);
    surface_bem_factory_pt = 0;
  }

  /// Get number of rows in the matrix
  long unsigned int nrow() const {return This_hmatrix_pt->rows;}

  /// Get number of cols in the matrix
  long unsigned int ncol() const {return This_hmatrix_pt->cols;}

  /// Get an entry from the matrix (note: don't implement things using this
  /// inside a loop: that removes the entire point of using H-matrices!)
  double operator()(const long unsigned int& i,
                    const long unsigned int& j) const
  {
    return getentry_supermatrix(This_hmatrix_pt, i, j);
  }

  /// Calculate sol = M.x. x is in the geometrical order defined by the
  /// order of the nodes in the mesh that was given to build. Soln is in
  /// the same order.
  void multiply(const DoubleVector& x, DoubleVector& soln) const
  {
    matrix_vector_multiply_build_check(x, soln);

    // intermediate storage vectors
    DoubleVector permuted_x(x.distribution_pt()),
      permuted_soln(x.distribution_pt());

    // permute x into this H-matrix's ordering
    for(unsigned j=0; j<x.nrow(); j++)
      {
        permuted_x[j] = x[Cluster_tree_pt->dof2idx[j]];
      }

    // actually multiply
    this->unpermuted_multiply(permuted_x, permuted_soln);

    // permute soln from this H-matrix's ordering back to the geometrical
    // ordering
    for(unsigned j=0; j<soln.nrow(); j++)
      {
        soln[j] = permuted_soln[Cluster_tree_pt->idx2dof[j]];
      }
  }

  /// Calculate soln = M.x, but *without* doing the permutations to x and
  /// soln which bring them back to the original geometrical
  /// ordering. Useful if there are two permutations be done which can be
  /// combined or if you need to do a large number of operations before
  /// returning values to other oomph-lib code.
  void unpermuted_multiply(const DoubleVector& x, DoubleVector& soln) const
  {
    matrix_vector_multiply_build_check(x, soln);
    eval_supermatrix(This_hmatrix_pt, x.values_pt(), soln.values_pt());
  }

  /// Mutliply a vector by the transpose of the matrix. Broken function for
  /// this because it's pure virtual in DoubleMatrixBase but I don't need
  /// it...
  void multiply_transpose(const DoubleVector& x, DoubleVector& soln) const
  {
    // Probably have this function if needed, look in supermatrix.h...
    throw OomphLibError("Function not yet implemented",
                        OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
  }

  // Parameters for H matrix compression
  unsigned nmin;
  double eps_aca;
  unsigned kmax;
  double eps;
  HLIB_HCA_VARIANT algorithm;
  double eta;
  unsigned quadorder;
  unsigned polyorder;

  /// Access function to the underlying hlib matrix representation.
  supermatrix* supermatrix_pt() const
  {
    return This_hmatrix_pt;
  }

  /// Access function to cluster tree (contains the map between mesh's
  /// geometrical ordering and the H matrix's ordering).
  clustertree* cluster_tree_pt() const
  {
    return Cluster_tree_pt;
  }

private:

  /// Data storage for hlib data structures
  supermatrix* This_hmatrix_pt;
  clustertree* Cluster_tree_pt;

  /// Run pre-multiply paranoid tests and maybe build soln.
  void matrix_vector_multiply_build_check(const DoubleVector& x,
                                          DoubleVector& soln) const
  {
    // If soln is not setup then setup the distribution
    if(!soln.built())
      {
        // Resize and initialize the solution vector
        LinearAlgebraDistribution dist(0, nrow(), false);
        soln.build(&dist, 0.0);
      }
    else // Otherwise just zero it
      {
        soln.initialise(0.0);
      }

#ifdef PARANOID
    // Check to see if x.nrow and matrix.ncol are the same
    if (this->ncol() != x.distribution_pt()->nrow())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The number of rows in the x vector and the number of columns in the "
          << "matrix must be the same";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

    // Check that solution and matrix nrow are the same
    if (soln.distribution_pt()->nrow() != this->nrow())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The soln vector is setup and therefore must have the same "
          << "distribution as the matrix";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
  }

};


// If we don't have hlib
#else

// Dummy class
class HMatrix : public DoubleMatrixBase
{
  public:

  HMatrix()
  {
    std::string err = "Hlib not installed, can't create HMatrix!";
    throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                        OOMPH_CURRENT_FUNCTION);
  }

  virtual ~HMatrix() {}

};


#endif

#endif
