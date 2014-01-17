

// Oomph-lib code
#include "generic.h"
#include "micromag.h"
#include "../../src/meshes/tetgen_mesh.h"


// hlib code
#include "bem3d.h"
#include "cluster.h"
// #include "hca.h"
#include "laplacebem.h"
#include "supermatrix.h"
#include "hcoarsening.h"


// TODO

// handle 2d bemgrids as well?

// Figure out parameters and maybe provide access functions

// Play with higher optimisation levels for hlib?

// This might only be valid for one potential (Laplace), might not be quite
// the same as the one I'm using...

using namespace TimingHelpers;
using namespace Factories;

namespace hlib
{

  /// Helper function: Given a mesh create a lookup from node pointer to
  /// node numbers in the mesh.
  inline std::map<Node*, unsigned> build_nodept2node
  (const Mesh& mesh, const BoundaryElementHandler& bem_handler,
   bool use_output=true)
    {
      std::map<Node*, unsigned> lookup;

      unsigned index = 99999;
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
          // dof in the bem matrix). ??ds should we use input or output? both?
          int g_eq = nd_pt->eqn_number(index);
          unsigned eq = 9999999999;
          if(use_output)
            {
              eq = bem_handler.output_lookup_pt()->global_to_node(g_eq);
            }
          else
            {
              eq = bem_handler.input_lookup_pt()->global_to_node(g_eq);
            }

          lookup[nd_pt] = eq;
        }

      return lookup;
    }


  /// Make a pair of integers which always has the first entry as the
  /// smaller of the two. Needed because pair(a, b) != pair(b, a)
  /// otherwise.
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

  /// Convert a (dense or hierarchical) matrix in hlib format to a dense
  /// matrix in oomph-lib format.
  inline void hlib_supermatrix2densedoublematrix (supermatrix& in,
                                                  clustertree ct,
                                                  DenseDoubleMatrix& out)
    {
#warning "not sure which permutation of the supermatrix brings us back to our ordering..."

      // hlib uses pointers for everything...
      supermatrix* in_pt = &in;

      // Create a full (dense) hlib matrix in a temporary data structure
      fullmatrix* full_pt = new_fullmatrix(in_pt->rows, in_pt->cols);
      convertsuper2_fullmatrix(full_pt, in_pt);

      // Just copy the data over (might be possible to do something fancier
      // if better performance is needed, maybe just assign full_pt->
      // somewhere?).
      const unsigned cols = full_pt->cols, rows = full_pt->rows;
      const double* elements = full_pt->e;
      out.resize(rows, cols);
      for(unsigned i=0; i<rows; i++)
        {
          const unsigned ii = ct.idx2dof[i];
          for(unsigned j=0; j<cols; j++)
            {
              const unsigned jj = ct.idx2dof[j];
              out(i, j) = elements[ii + jj*rows];
            }
        }

      // Clean up the temporary data
      del_fullmatrix(full_pt);
    }


  class HMatrix : public DoubleMatrixBase
  {
  public:

    HMatrix(const unsigned& _nmin)
    {
#ifdef OOMPH_HAS_MPI
      throw OomphLibError("Not implemented for MPI, it might work for non-distributed memory though... But untested",
                          OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);

#endif

      // Random numbers for parameters, some taken from Knittel's thesis
      nmin = _nmin;
      eps_aca = 1e-7;
      kmax = 500;
      eps = 1e-3;
      algorithm = HLIB_HCAII;
      eta = 0.25;
      quadorder = 3;
      polyorder = 4;

      // Null pointers
      This_hmatrix_pt = 0;
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
    }


    /// Create an equivalent dense matrix
    void todense(DenseDoubleMatrix& out) const
      {
        hlib::hlib_supermatrix2densedoublematrix(*This_hmatrix_pt,
                                                 *Cluster_tree_pt,
                                                 out);
      }

    /// Build the hmatrix for a given mesh
    void build(const Mesh& mesh, const BoundaryElementHandler* bem_handler_pt)
    {
#ifdef PARANOID
      if(mesh.finite_element_pt(0)->nnode() != 3 ||
         mesh.finite_element_pt(0)->dim() != 2)
        {
          std::string err = "Only works for triangles";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif

      if(This_hmatrix_pt != 0)
        {
          std::string err = "Already built a hmatrix in here, replacement not yet implemented";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
          // To implement this just make sure you delete all the data
          // structures (only if pointers non-null) and null their pointers.
        }


      // Copy the oomph mesh data into hlib data structures
      // ============================================================

      // Create a reverse lookup scheme which finds the node number when
      // given a node pointer.
      std::map<Node*, unsigned> nodept2node =
        build_nodept2node(mesh, *bem_handler_pt, false);

      // Get a (unique) list of all the edges in the mesh
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
      std::vector<Node*> nodes(vertices);
      for(unsigned nd=0, nnd=mesh.nnode(); nd<nnd; nd++)
        {
          nodes[nd] = mesh.node_pt(nd);
        }

      std::random_shuffle(nodes.begin(), nodes.end());
      for(unsigned nd=0, nnd=nodes.size(); nd<nnd; nd++)
        {
          Node* nd_pt = nodes[nd];
          unsigned nd_num = nodept2node.find(nd_pt)->second;

          bem_grid_pt->x[nd_num][0] = nd_pt->x(0);
          bem_grid_pt->x[nd_num][1] = nd_pt->x(1);
          bem_grid_pt->x[nd_num][2] = nd_pt->x(2);
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
          // of the unit normal.
          if(ele_pt->normal_sign() > 0)
            {
              std::swap(bem_grid_pt->t[ele][0],
                        bem_grid_pt->t[ele][1]);
            }
        }

      std::random_shuffle(bem_grid_pt->t, bem_grid_pt->t+triangles);

      std::vector<std::pair<int, int> > edges_vec(edges_set.begin(),
                                                  edges_set.end());
      // std::set<std::pair<int, int> >::const_iterator it;
      // unsigned i_edge = 0;
      // for(it = edges_set.begin(); it != edges_set.end(); ++it)
      //   {
      //     bem_grid_pt->e[i_edge][0] = it->first;
      //     bem_grid_pt->e[i_edge][1] = it->second;
      //     i_edge++;
      //   }

      // shuffle vector
      std::random_shuffle(edges_vec.begin(), edges_vec.end());

      // copy to array
      std::vector<std::pair<int, int> >::const_iterator it;
      unsigned i_edge = 0;
      for(it = edges_vec.begin(); it != edges_vec.end(); ++it)
        {
          bem_grid_pt->e[i_edge][0] = it->first;
          bem_grid_pt->e[i_edge][1] = it->second;
          i_edge++;
        }

      // Calculate required values in the grid
      prepare_bemgrid3d(bem_grid_pt);

      // Construct the matrix itself, via some other intermediate data
      // structures.
      // ============================================================

      // Build the cluster tree ??ds what are 0 and HLIB_REGULAR?
      Cluster_tree_pt = buildvertexcluster_bemgrid3d(bem_grid_pt,
                                                     HLIB_REGULAR,
                                                     nmin, 0);

      // Build the factory which will finally create the hmatrix
      // ??ds parameter list?

// #warning "Using galerkin H-bem"
//       surfacebemfactory* surface_bem_factory_pt = new_surfacebemfactory_dlp
//         (bem_grid_pt, HLIB_LINEAR_BASIS, Cluster_tree_pt,
//          HLIB_LINEAR_BASIS, Cluster_tree_pt,
//          quadorder, quadorder,
//          polyorder, 0.0);

      surfacebemfactory* surface_bem_factory_pt = new_surfacebemfactory_dlp_collocation
        (bem_grid_pt, HLIB_LINEAR_BASIS, Cluster_tree_pt,
         bem_grid_pt, HLIB_LINEAR_BASIS, Cluster_tree_pt,
         quadorder, quadorder,
         polyorder, 0.0);


      // Finally build the matrix itself, with adaptive recompression
      // ??ds parameter list?
      This_hmatrix_pt = onthefly_hca_coarsen_supermatrix
        (Cluster_tree_pt->root, Cluster_tree_pt->root,
         surface_bem_factory_pt,
         eps_aca, kmax, eps, 1, algorithm, 0, eta, 0);

      // // Output idx2dof map
      // for(unsigned j=0; j<Cluster_tree_pt->ndof; j++)
      //   {
      //     std::cout << j << " " << Cluster_tree_pt->idx2dof[j] << std::endl;
      //   }


      // Clean up auxilary data structures
      // ============================================================

      // Since theses are all C data types they don't have destructors to do
      // this for us, so do it manually now.
      del_bemgrid3d(bem_grid_pt);
      bem_grid_pt = 0;

      // Used to get mapping between hlib-dofs and real dofs
      // del_clustertree(cluster_tree_pt);
      // cluster_tree_pt = 0;

      del_surfacebemfactory(surface_bem_factory_pt);
      surface_bem_factory_pt = 0;

      // It's plausible that some of these data structures could actually
      // be needed later if operations other than multiplication are used
      // (and I'm not going to try to test them all!)... if so feel free to
      // make these data structres class variables and to not delete them
      // until the destructor.
    }


    long unsigned int nrow() const {return This_hmatrix_pt->rows;}


    long unsigned int ncol() const {return This_hmatrix_pt->cols;}


    double operator()(const long unsigned int& i,
                      const long unsigned int& j) const
    {
      return getentry_supermatrix(This_hmatrix_pt, i, j);
    }


    void multiply(const DoubleVector& x, DoubleVector& soln) const
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

      // Do the multiplication
      eval_supermatrix(This_hmatrix_pt, x.values_pt(), soln.values_pt());

      //??ds reorder result?
    }

    void multiply_transpose(const DoubleVector& x, DoubleVector& soln) const
    {
      // Probably have this function if needed, look in supermatrix.h...
      throw OomphLibError("Function not yet implemented",
                          OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
    }

    // Parameters
    unsigned nmin;
    double eps_aca;
    unsigned kmax;
    double eps;
    HLIB_HCA_VARIANT algorithm;
    double eta;
    unsigned quadorder;
    unsigned polyorder;

    /// Access function to the underlying hlib matrix representation.
    supermatrix* supermatrix_pt()
      {
        return This_hmatrix_pt;
      }

  private:

    /// Data storage for hlib data structures
    supermatrix* This_hmatrix_pt;
    clustertree* Cluster_tree_pt;

  };

}


int main(int argc, char *argv[])
{

  // First make sure that we can create/destroy empty h matrices correctly
  // ============================================================
  {
    hlib::HMatrix hmat(50);
  }

  // Make the problem to get bem mesh etc.
  // ============================================================

  // Need a dummy problem so that we can get equation numbers :(
  MMArgs args;
  args.parse(argc, argv);
  LLGProblem problem;
  problem.Residual_calculator_pt
    = LLGFactories::residual_calculator_factory("llg");
  problem.add_time_stepper_pt(args.time_stepper_pt);
  args.assign_specific_parameters(&problem);
  // disable corner angles in bem matrx, to make matrices easier to compare
  problem.Disable_bem_corners = true;
  problem.build(args.mesh_pts);

  // Some useful pointers
  DenseDoubleMatrix* old_bem_matrix_pt = problem.Bem_handler_pt->bem_matrix_pt();
  const Mesh* surface_mesh_pt = problem.Bem_handler_pt->bem_mesh_pt();

  // Vectors for multiply test
  const unsigned nbemnodes = surface_mesh_pt->nnode();
  LinearAlgebraDistribution dist(0, nbemnodes, false);
  DoubleVector x(dist), H_soln(dist), dense_H_soln(dist), oomph_soln(dist);

  // Random entries for x vector
  Vector<double> v = VectorOps::random_vector(nbemnodes);
  x.initialise(v);


  // Multiply with pure-oomph bem matrix
  old_bem_matrix_pt->multiply(x, oomph_soln);


  // Make a H matrix and use it in a multiply
  // ============================================================

  // Use the bem mesh to make a h-matrix

  unsigned nmin = 500; // ??ds works
  // unsigned nmin = 5;
  hlib::HMatrix hmat(nmin);
  hmat.build(*surface_mesh_pt, problem.Bem_handler_pt);

  // Dump rank info
  outputrank_supermatrix(hmat.supermatrix_pt(), (const char*)"rank.ps");

  // Do the multiply
  hmat.multiply(x, H_soln);

  // Flip the sign ??ds not sure why we need to do this
  H_soln *= -1;


  // Make a dense matrix from the H matrix and try multiplying with that
  // ============================================================

  // Copy H matrix into an oomph-lib matrix
  DenseDoubleMatrix double_matrix;
  hmat.todense(double_matrix);

  // Flip sign ??ds not sure why we need this to make it match
  for(unsigned i=0; i<double_matrix.nrow(); i++)
    {
      for(unsigned j=0; j<double_matrix.ncol(); j++)
        {
          double_matrix(i,j) *= -1;
        }
    }

  // Multiply
  double_matrix.multiply(x, dense_H_soln);



  // // Output the timing and error results
  // // ============================================================
  // std::cout << "H-matrix build time " << hbuildstop - hbuildstart << std::endl;
  // std::cout << "dense build time " << dbuildstop - dbuildstart << std::endl;
  // std::cout << "oldbem build time " << oldbembuildstop - oldbembuildstart << std::endl;
  // std::cout << std::endl;
  // std::cout << "H-matrix multiply time " <<  hstop - hstart << std::endl;
  // std::cout << "dense multiply time " << dstop - dstart << std::endl;
  // std::cout << "old bem dense multiply time " << oldbemstop - oldbemstart
  //           << std::endl;


  // Dump matrices to files
  // ============================================================
  double_matrix.output("Validation/new_bem_matrix");
  old_bem_matrix_pt->output("Validation/old_bem_matrix");


  double H_soln_vs_dense = two_norm_diff(H_soln, dense_H_soln);
  double H_soln_vs_oomph = two_norm_diff(H_soln, oomph_soln);


  // Some analysis
  // ============================================================
  std::cout << std::endl;
  std::cout << "muliplication 2 norm of H-matrix vs dense version: "
            << H_soln_vs_dense << std::endl;
  std::cout << "muliplication 2 norm of H-matrix vs oomph bem: "
            << H_soln_vs_oomph << std::endl;

  //??ds need to fix ordering of H-matrix multiply output (and input?)
  // if(H_soln_vs_dense > 1e-8)
  //   {
  //     std::cout << "Test failed due to error in dense H." << std::endl;
  //     return 1;
  //   }

  // if(H_soln_vs_oomph > 1e-8)
  //   {
  //     std::cout << "Test failed due to error in oomph vs H." << std::endl;
  //     return 2;
  //   }


  // DoubleVector sold, snew, rhs(dist, 1.0);
  // double_matrix.solve(rhs, snew);
  // old_bem_matrix_pt->solve(rhs, sold);

  // // for(unsigned j=0; j<snew.nrow(); j++)
  // //   {
  // //     std::cout << j << ": " << snew[j] << " " << sold[j] << std::endl;
  // //   }

  return 0;
}
