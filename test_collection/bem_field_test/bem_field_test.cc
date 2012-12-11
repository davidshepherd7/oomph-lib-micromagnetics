
// Floating point debugging
#include <fenv.h>

#include "generic.h"
#include "meshes/tetgen_mesh.h"

#include "../../boundary_element_handler.h"
#include "../../generic_poisson_problem.h"
#include "../../micromagnetics_boundary_element.h"

using namespace oomph;

namespace oomph
{

  // Don't change these typedefs! They are needed for compatability with
  // poisson elements

  /// \short Function pointer to the prescribed-flux function fct(x,f(x)) --
  /// x is a Vector!
  typedef void (*PoissonPrescribedFluxFctPt)
  (const Vector<double>& x, double& flux);

  /// \short Function pointer to source function fct(x,f(x)) --
  /// x is a Vector!
  typedef void (*PoissonSourceFctPt)(const Vector<double>& x, double& f);

  // =================================================================
  /// Solve for the magnetostatic field in a mesh with known magnetisation
  /// using the hybrid FEM/BEM method.
  // =================================================================
  template<class ELEMENT>
  class ExplicitHybridMagnetostaticFieldProblem
  {
  public:

    /// Constructor
    ExplicitHybridMagnetostaticFieldProblem
    (Mesh* phi_1_mesh_pt, Mesh* phi_mesh_pt,
     const PoissonSourceFctPt divM_pt,
     const PoissonPrescribedFluxFctPt Mdotn_pt)
    {
      // Set up phi_1 problem
      Phi_1_problem.set_bulk_mesh(phi_1_mesh_pt);
      Phi_1_problem.set_source_fct_pt(divM_pt);
      for(unsigned b=0, nb=phi_1_mesh_pt->nboundary(); b < nb; b++)
        {
          // Phi_1 b.c.s are all Neumann
          Phi_1_problem.set_neumann_boundary(b, Mdotn_pt);
        }

      // Pin a single phi_1 value so that the problem is fully determined.
      // This is necessary to avoid having a free constant of integration
      // (which causes scaling problems). Just pin the first non boundary
      // node ??ds not sure this is ok...
      Node* pinned_phi_1_node_pt = phi_1_mesh_pt->get_non_boundary_node();
      pinned_phi_1_node_pt->pin(0);
      pinned_phi_1_node_pt->set_value(0,0.0);
      Phi_1_problem.build();


      // Construct the BEM (must be done before pinning phi values)
      Bem_handler.set_bem_all_boundaries(phi_1_mesh_pt);
      // both zero because they are in seperate problems
      Bem_handler.input_index() = 0;
      Bem_handler.output_index() = 0;
      Bem_handler.integration_scheme_pt() = new TVariableOrderGaussLegendre<2>; //??ds memory leak
      Bem_handler.build();


      // Now we can set up phi problem
      Phi_problem.set_bulk_mesh(phi_mesh_pt);
      unsigned nboundary = phi_mesh_pt->nboundary();
      Phi_boundary_values_pts.assign(nboundary, 0);
      for(unsigned b=0; b < nboundary; b++)
        {
          // Phi is determined by BEM
          LinearAlgebraDistribution* dist_pt =
            new LinearAlgebraDistribution(0, phi_mesh_pt->nboundary_node(b), false);

          Phi_boundary_values_pts[b] = new DoubleVector(dist_pt);
          Phi_problem.set_dirichlet_boundary_by_vector(b, Phi_boundary_values_pts[b]);
        }
      Phi_problem.build();

      // Set up linear solvers
      Phi_problem.linear_solver_pt() = new CG<CRDoubleMatrix>;
      Phi_1_problem.linear_solver_pt() = new CG<CRDoubleMatrix>;

      // Cast to iterative solver pointers
      IterativeLinearSolver* phi_it_solver_pt =
        dynamic_cast<IterativeLinearSolver*>(Phi_problem.linear_solver_pt());

      IterativeLinearSolver* phi_1_it_solver_pt =
        dynamic_cast<IterativeLinearSolver*>(Phi_1_problem.linear_solver_pt());

#ifdef OOMPH_HAS_HYPRE
      // AMG preconditioners
      HyprePreconditioner *amg_phi_pt, *amg_phi_1_pt;
      amg_phi_pt = new HyprePreconditioner;
      phi_it_solver_pt->preconditioner_pt() = amg_phi_pt;
      amg_phi_pt->hypre_method() = HyprePreconditioner::BoomerAMG;

      amg_phi_1_pt = new HyprePreconditioner;
      phi_1_it_solver_pt->preconditioner_pt() = amg_phi_1_pt;
      amg_phi_1_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
#endif

      // // ILU0 preconditioners
      // phi_it_solver_pt->preconditioner_pt()=new ILUZeroPreconditioner<CRDoubleMatrix>;
      // phi_1_it_solver_pt->preconditioner_pt()=new ILUZeroPreconditioner<CRDoubleMatrix>;


    };


    /// Blank destructor
    ~ExplicitHybridMagnetostaticFieldProblem() {}

    void solve()
    {
      // Solve for phi1
      Phi_1_problem.newton_solve();

      // Update boundary values of phi
      Bem_handler.get_bem_values(Phi_boundary_values_pts);

      // Solve for phi
      Phi_problem.newton_solve();
    }

    /// Calculate the magnetostatic field at each node (= - grad phi).
    void average_magnetostatic_field(Vector<double> &magnetostatic_field) const;

    /// Output solution (phi values only).
    void doc_solution(const unsigned &label) const;

    // Access functions
    // ============================================================

    /// \short Const access function for Phi_mesh_pt.
    Mesh* phi_mesh_pt() const {return Phi_problem.bulk_mesh_pt();}

    /// \short Const access function for Phi_mesh_pt.
    Mesh* phi1_mesh_pt() const {return Phi_1_problem.bulk_mesh_pt();}

  private:

    /// Bem handler object - provide functions and data needed for hybird
    /// BEM method (dim hardcoded to 3, probably ok).
    BoundaryElementHandler<MicromagFaceElement<ELEMENT> > Bem_handler;

    /// The problem for the preliminary poisson solve to get boundary
    /// conditions on the magnetostatic potential solve.
    GenericPoissonProblem<ELEMENT> Phi_1_problem;

    /// The problem for the "real" magnetostatic potential.
    GenericPoissonProblem<ELEMENT> Phi_problem;

    /// Intermediate storage for results of bem (ideally we would have it
    /// call a function to get the boundary values filled in but c++ member
    /// functions pointers are useless...
    Vector<DoubleVector*> Phi_boundary_values_pts;

  };

  // =================================================================
  /// Calculate the average field in the mesh (= - grad phi).
  // =================================================================
  template<class ELEMENT>
  void ExplicitHybridMagnetostaticFieldProblem<ELEMENT>::
  average_magnetostatic_field(Vector<double> &average_magnetostatic_field) const
  {
    // Pick a point in the middle of the element
    const Vector<double> s(3,0.3);
    Vector<double> total_dphidx(3,0.0);

    // Loop over all elements calculating the value in the middle of the element
    for(unsigned e=0, ne=phi_mesh_pt()->nelement(); e < ne; e++)
      {
        ELEMENT* ele_pt = dynamic_cast<ELEMENT*>
          (phi_mesh_pt()->element_pt(e));

        // Get the shape function and eulerian coordinate derivative at
        // position s.
        unsigned n_node = ele_pt->nnode();
        Shape psi(n_node); DShape dpsidx(n_node,3);
        ele_pt->dshape_eulerian(s,psi,dpsidx);

        // Interpolate grad phi
        Vector<double> interpolated_dphidx(3,0.0);
        for(unsigned l=0;l<n_node;l++)
          {
            double phi_value = ele_pt->raw_nodal_value(l,0);
            for(unsigned i=0; i<3; i++)
              {interpolated_dphidx[i] += phi_value*dpsidx(l,i);}
          }

        // Add this grad phi to the sum
        for(unsigned j=0; j<3; j++)
          {
            total_dphidx[j] += interpolated_dphidx[j];
          }
      }

    // Divide sum by number of elements to get the average. Take the
    // negative to get the field.
    double nele = double(phi_mesh_pt()->nelement());
    for(unsigned j=0; j<3; j++)
      {
        average_magnetostatic_field[j] = - total_dphidx[j] / nele;
      }
  }

  // =================================================================
  /// Output solution
  // =================================================================
  template<class ELEMENT>
  void ExplicitHybridMagnetostaticFieldProblem<ELEMENT>::
  doc_solution(const unsigned &label) const
  {
    std::ofstream some_file;
    char filename[100];

    // Number of plot points
    unsigned npts;
    npts=3;

    // Output solution with specified number of plot points per element
    sprintf(filename,"results/phi_solution%i.dat",label);
    some_file.open(filename);
    phi_mesh_pt()->output(some_file,npts);
    some_file.close();

    sprintf(filename,"results/phi1_solution%i.dat",label);
    some_file.open(filename);
    phi1_mesh_pt()->output(some_file,npts);
    some_file.close();
  }



} // End of oomph namespace

namespace Inputs
{
  //??ds try some others eventually?

  void exact_M(const Vector<double> &x, Vector<double> &M)
  {
    M.assign(3,0.0);
    M[0] = 1;
  }

  void unit_surface_normal(const Vector<double> &x, Vector<double> &n)
  {
    // n = r/|r| = normalised position since this is a sphere
    n = x;
    normalise(n);
  }

  void div_M(const Vector<double>& x, double& f)
  {
    // Finite difference it...
    double eps = 1e-10;

    Vector<double> Mxp, Mxm, Myp, Mym, Mzp, Mzm;
    Vector<double> Xxp(x), Xxm(x), Xyp(x), Xym(x), Xzp(x), Xzm(x);
    Xxp[0] += eps; Xxm[0] -= eps;
    Xyp[1] += eps; Xym[1] -= eps;
    Xzp[2] += eps; Xzm[2] -= eps;

    exact_M(Xxp,Mxp); exact_M(Xxm,Mxm);
    exact_M(Xyp,Myp); exact_M(Xym,Mym);
    exact_M(Xzp,Mzp); exact_M(Xzm,Mzm);

    f = Mxp[0] - Mxm[0]
      + Myp[1] - Mym[1]
      + Mzp[2] - Mzm[2];
  }

  void Mdotn(const Vector<double>& x, double& f)
  {
    Vector<double> M, n;
    exact_M(x,M);
    unit_surface_normal(x,n);

    f = VectorOps::dot(M,n);
  }

}

int main()
{
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  // Make a spherical mesh from the listed files generated previously.
  TetgenMesh<TPoissonElement<3,2> > sphere_mesh1("mesh.1.node",
                                                 "mesh.1.ele",
                                                 "mesh.1.face");
  // Make a spherical mesh from the listed files generated previously.
  TetgenMesh<TPoissonElement<3,2> > sphere_mesh("mesh.1.node",
                                                "mesh.1.ele",
                                                "mesh.1.face");
  //??ds - must be a more elegant way to have two identical meshes... or
  // avoid having two altogether... how can we be sure that we don't break
  // everything by applying changes (e.g. refinement) to only one mesh?


  // Set a divergence of M (M = constant so divM = 0).
  PoissonSourceFctPt divM_pt = &Inputs::div_M;

  // Set flux conditioners on phi_1 (= Mdotn).
  PoissonPrescribedFluxFctPt Mdotn_pt = &Inputs::Mdotn;

  // Make a hybrid problem
  ExplicitHybridMagnetostaticFieldProblem<TPoissonElement<3,2> >
    magnetostatic_problem(&sphere_mesh1, &sphere_mesh,
                          divM_pt, Mdotn_pt);

  // Solve it
  magnetostatic_problem.solve();

  // Output results
  magnetostatic_problem.doc_solution(0);

  // Compute field
  Vector<double> average_field(3,0.0);
  magnetostatic_problem.average_magnetostatic_field(average_field);
  double rel_err = std::abs((average_field[0] - (-1.0/3.0))/average_field[0]);
  double abs_err = std::abs(average_field[0] - (-1.0/3.0));

  std::cout << average_field << std::endl;

  // Check close/convergence to analytical value
  std::cout
    << "Field values are [" << average_field[0] << ","
    << average_field[1] << "," << average_field[2] << "]\n\n"

    << "Error compared to an exact sphere is:\n"
    << "abs error = " << abs_err << "\n"
    << "rel error = " << rel_err << std::endl;

  // Return error status. We don't expect it to be very accurate because we
  // use a pretty bad approximation to a sphere for speed reasons.
  if(rel_err > 0.25)
    {
      std::cout << std::endl
                << "Relative error is too large!" << std::endl;
      return 1;
    }

  // Also check that y,z components are *all* close to zero (not average)

  // Also check that std-dev of x is not too high?
  return 0;
}
