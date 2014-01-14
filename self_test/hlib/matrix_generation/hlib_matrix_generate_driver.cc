

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

namespace hlib
{

  /// Helper function: Given a mesh create a lookup from node pointer to
  /// node numbers in the mesh.
  inline std::map<Node*, unsigned> build_nodept2node(const Mesh& mesh)
    {
      std::map<Node*, unsigned> lookup;

      for(unsigned nd=0, nnd=mesh.nnode(); nd<nnd; nd++)
        {
          Node* nd_pt = mesh.node_pt(nd);
          lookup[nd_pt] = nd;
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
                                                  DenseDoubleMatrix& out)
    {
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
          for(unsigned j=0; j<cols; j++)
            {
              out(i, j) = elements[i + j*rows];
            }
        }

      // Clean up the temporary data
      del_fullmatrix(full_pt);
    }


  class HMatrix : public DoubleMatrixBase
  {
  public:

    HMatrix()
    {
#ifdef OOMPH_HAS_MPI
      throw OomphLibError("Not implemented for MPI, it might work for non-distributed memory though... But untested",
                          OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);

#endif

      // Random numbers for parameters, some taken from Knittel's thesis
      nmin = 30;
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
        hlib::hlib_supermatrix2densedoublematrix(*This_hmatrix_pt, out);
      }

    /// Build the hmatrix for a given mesh
    void build(const Mesh& mesh)
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

      // "Copy" the oomph mesh into hlib data structures
      // ============================================================

      // Create a reverse lookup scheme which finds the node number when
      // given a node pointer.
      std::map<Node*, unsigned> nodept2node = build_nodept2node(mesh);

      // Get a (unique) list of all the edges in the mesh
      std::set<std::pair<int, int> > edges_set;
      for(unsigned ele=0, nele=mesh.nelement(); ele<nele; ele++)
        {
          FiniteElement* ele_pt = mesh.finite_element_pt(ele);

          unsigned nd0 = nodept2node.find(ele_pt->node_pt(0))->second;
          unsigned nd1 = nodept2node.find(ele_pt->node_pt(0))->second;
          unsigned nd2 = nodept2node.find(ele_pt->node_pt(0))->second;

          // Insert all pairs of the nodes into the set (function makes
          // sure they are in order inside the pair so that edge {0, 1} ==
          // edge {1, 0}.
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
          FiniteElement* ele_pt = mesh.finite_element_pt(ele);
          for(unsigned nd=0, nnd=ele_pt->nnode(); nd<nnd; nd++)
            {
              Node* nd_pt = ele_pt->node_pt(nd);
              bem_grid_pt->t[ele][nd] = nodept2node[nd_pt];
            }
        }

      std::set<std::pair<int, int> >::const_iterator it;
      unsigned i_edge = 0;
      for(it = edges_set.begin(); it != edges_set.end(); ++it)
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
      clustertree* cluster_tree_pt = buildvertexcluster_bemgrid3d(bem_grid_pt,
                                                                  HLIB_REGULAR,
                                                                  nmin, 0);

      // Build the factory which will finally create the hmatrix
      // ??ds parameter list?
      surfacebemfactory* surface_bem_factory_pt = new_surfacebemfactory_dlp_collocation
        (bem_grid_pt, HLIB_LINEAR_BASIS, cluster_tree_pt,
         bem_grid_pt, HLIB_LINEAR_BASIS, cluster_tree_pt,
         quadorder, quadorder,
         polyorder, 0.0);


      // Finally build the matrix itself, with adaptive recompression
      // ??ds parameter list?
      This_hmatrix_pt = onthefly_hca_coarsen_supermatrix
        (cluster_tree_pt->root, cluster_tree_pt->root,
         surface_bem_factory_pt,
         eps_aca, kmax, eps, 1, algorithm, 0, eta, 0);



      // Clean up auxilary data structures
      // ============================================================

      // Since theses are all C data types they don't have destructors to do
      // this for us, so do it manually now.
      del_bemgrid3d(bem_grid_pt);
      bem_grid_pt = 0;

      del_clustertree(cluster_tree_pt);
      cluster_tree_pt = 0;

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

  };

}


int main()
{

  // First make sure that we can create/destroy empty h matrices correctly
  // ============================================================
  {
    hlib::HMatrix hmat;
  }

  // Now build a real h-matrix and do some multiplication
  // ============================================================

  // Make a mesh
  const std::string refinement_level = "1";
  bool run_old_bem = false; // old version is very slow by comparison

  TetgenMesh<TMicromagElement<3, 2> > bulk_mesh
    ("./meshes/sphere." + to_string(refinement_level) + ".node",
     "./meshes/sphere." + to_string(refinement_level) + ".ele",
     "./meshes/sphere." + to_string(refinement_level) + ".face");
  bulk_mesh.setup_boundary_element_info();

  // Loop over boundaries, create a list of boundaries to include in bem
  // (all of them).
  Vector<std::pair<unsigned, const Mesh*> > bem_boundaries;
  for(unsigned b=0, nb=bulk_mesh.nboundary(); b<nb; b++)
    {
      bem_boundaries.push_back(std::make_pair(b, &bulk_mesh));
    }

  // Build a bem mesh for it
  Mesh surface_mesh;
  surface_mesh.nodes_can_be_deleted() = false;
  BEMElementFactoryFctPt bem_element_factory_fpt
    = LLGFactories::bem_element_factory_factory(bulk_mesh.finite_element_pt(0));
  build_bem_mesh_helper(bem_boundaries, bem_element_factory_fpt,
                        0, surface_mesh);

  const unsigned nbemnodes = surface_mesh.nnode();

  // Use the bem mesh to make a h-matrix
  double hbuildstart = timer();
  hlib::HMatrix hmat;
  hmat.build(surface_mesh);
  double hbuildstop = timer();

  outputrank_supermatrix(hmat.supermatrix_pt(), (const char*)"rank.ps");

  // for(unsigned j=0; j<hmat.nrow(); j++)
  //   {
  //     std::cout << hmat(j, j) << std::endl;
  //   }


  Vector<double> v = VectorOps::random_vector(nbemnodes);
  LinearAlgebraDistribution dist(0, nbemnodes, false);
  DoubleVector x(dist), soln(dist), oomphsoln(dist), oldbemsoln(dist);
  x.initialise(v);

  double hstart = timer();
  hmat.multiply(x, soln);
  double hstop = timer();


  // Copy into an oomph-lib matrix
  double dbuildstart = timer();
  DenseDoubleMatrix double_matrix;
  hmat.todense(double_matrix);
  double dbuildstop = timer();

  double dstart = timer();
  double_matrix.multiply(x, oomphsoln);
  double dstop = timer();

  double oldbembuildstart, oldbembuildstop, oldbemstart, oldbemstop;
  oldbembuildstart = timer();
  BoundaryElementHandler bem_handler;
  bem_handler.Bem_element_factory_fpt = bem_element_factory_fpt;
  bem_handler.integration_scheme_pt() = LLGFactories::
    variable_order_integrator_factory(bulk_mesh.finite_element_pt(0));
  bem_handler.set_input_index(1);
  bem_handler.set_output_index(0);
  bem_handler.set_bem_all_boundaries(&bulk_mesh);
  bem_handler.input_corner_data_pt() = 0;
  bem_handler.build();
  DenseDoubleMatrix* old_bem_matrix_pt = bem_handler.bem_matrix_pt();
  oldbembuildstop = timer();

  oldbemstart = timer();
  old_bem_matrix_pt->multiply(x, oldbemsoln);
  oldbemstop = timer();


  // Create the diagonal matrix of solid angles (since h-lib created
  // matrices don't include this we need to add it ourselves). ??ds could
  // just use a dot product...
  // ============================================================
  DenseDoubleMatrix corner_matrix(nbemnodes);
  bem_handler.corner_list_pt()->add_corner_contributions(corner_matrix);

  // Calculate contribution
  DoubleVector corner_soln_contribution(dist);
  corner_matrix.multiply(x, corner_soln_contribution);

  // Add to the solutions
  soln += corner_soln_contribution;




  std::cout << "h build time " << hbuildstop - hbuildstart << std::endl;
  std::cout << "dense build time " << dbuildstop - dbuildstart << std::endl;
  std::cout << "oldbem build time " << oldbembuildstop - oldbembuildstart << std::endl;
  std::cout << std::endl;
  std::cout << "h multiply time " <<  hstop - hstart << std::endl;
  std::cout << "dense multiply time " << dstop - dstart << std::endl;
  std::cout << "old bem dense multiply time " << oldbemstop - oldbemstart
            << std::endl;

  std::cout << std::endl;
  std::cout << "2 norm : h vs dense "
            << two_norm_diff(soln, oomphsoln) << std::endl;
  std::cout << "2 norm : h vs old bem "
            << two_norm_diff(soln, oldbemsoln) << std::endl;


  // Dump matrices to files
  bem_handler.corner_list_pt()->add_corner_contributions(double_matrix);
  double_matrix.output("Validation/new_bem_matrix");
  old_bem_matrix_pt->output("Validation/old_bem_matrix");

  return 0;
}
