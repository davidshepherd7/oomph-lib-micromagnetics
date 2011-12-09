/// Boundary element classes

#include "generic.h"
#include "meshes/rectangular_quadmesh.h"


using namespace oomph;

namespace oomph
{

  //===============================================================
  /// A class for the maths used in magnetostatic boundary elements
  //===============================================================
  template<unsigned DIM>
  class MagnetostaticBEEquations : public virtual FiniteElement
  {
  private:

  public:

    /// Constructor
    MagnetostaticBEEquations(){}

    /// Broken copy constructor
    MagnetostaticBEEquations(const MagnetostaticBEEquations& dummy)
    {
      BrokenCopy::broken_copy("MagnetostaticBEEquations");
    }

    /// Broken assignment operator
    void operator=(const MagnetostaticBEEquations&)
    {
      BrokenCopy::broken_assign("MagnetostaticBEEquations");
    }

    /// Self-test: Return 0 for OK.
    unsigned self_test(){return 0;} //??ds write a real test sometime

    /// Add the element's contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      //Call the generic residuals function with flag set to 0 using a dummy matrix argument
      fill_in_generic_residual_contribution_micromag
	(residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    /// Output function at default number of plot points
    void output(std::ostream &outfile);

    /// Output function at n_plot^DIM points
    void output(std::ostream &outfile, const unsigned &n_plot);

    /// Fill in contribution to residuals and jacobian (if flag is set) from these equations (compatible with multiphysics)
    void fill_in_generic_residual_contribution_micromag(Vector<double> &residuals, DenseMatrix<double> &jacobian, const unsigned& flag);

    /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
    virtual double dshape_and_dtest_eulerian_micromag(const Vector<double> &s, Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const=0;

    /// Shape, test functions & derivs. w.r.t. to global coords. at integration point ipt. Return Jacobian.
    virtual double dshape_and_dtest_eulerian_at_knot_micromag(const unsigned& ipt, Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const=0;

  }; // End of MagnetostaticBEEquations class


  // /// Fill in contribution to residuals and jacobian (if flag is set) from these equations (compatible with multiphysics)
  // void MagnetostaticBEEquations<DIM> fill_in_generic_residual_contribution_micromag(Vector<double> &residuals, DenseMatrix<double> &jacobian, const unsigned& flag)
  // {

  // }

  // /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
  // virtual double MagnetostaticBEEquations<DIM>
  // dshape_and_dtest_eulerian_micromag(const Vector<double> &s, Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx)
  // {

  // }

  // /// Shape, test functions & derivs. w.r.t. to global coords. at integration point ipt. Return Jacobian.
  // MagnetostaticBEEquations<DIM>
  // dshape_and_dtest_eulerian_at_knot_micromag(const unsigned& ipt, Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx)
  // {

  // }












  //==================================================
  /// A class combining DIM-1 QElements with the MagnetostaticBEEquations class
  //==================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class QMagnetostaticBEElement : public QElement<DIM,NNODE_1D>, public MagnetostaticBEEquations<DIM>
  {

  private:

  public:

    /// Constructor: Call constructors for QElement and magnetostatic boundary element equations
    QMagnetostaticBEElement() : QElement<DIM,NNODE_1D>(), MagnetostaticBEEquations<DIM>()
    {}

    /// Broken copy constructor
    QMagnetostaticBEElement(const QMagnetostaticBEElement<DIM,NNODE_1D>& dummy)
    {
      BrokenCopy::broken_copy("QMagnetostaticBEElement");
    }

    /// Broken assignment operator
    void operator=(const QMagnetostaticBEElement<DIM,NNODE_1D>&)
    {
      BrokenCopy::broken_assign("QMagnetostaticBEElement");
    }

    /// Required  # of `values' (pinned or dofs) at node n.
    // We need to store phi_2 ??ds any others
    //??ds generalise somehow?
    inline unsigned required_nvalue(const unsigned &n) const
    {return 1;}

    /// Output function at default number of plot points
    void output(std::ostream &outfile)
    {MagnetostaticBEEquations<DIM>::output(outfile);}

    /// Output function at n_plot^DIM points
    void output(std::ostream &outfile, const unsigned &n_plot)
    {MagnetostaticBEEquations<DIM>::output(outfile,n_plot);}

    // /// C-style output function: x,y,u or x,y,z,u
    // void output(FILE* file_pt)
    // {MicromagEquations<DIM>::output(file_pt);}

    // /// C-style output function: x,y,u or x,y,z,u at n_plot^DIM plot points
    // void output(FILE* file_pt, const unsigned &n_plot)
    // {MicromagEquations<DIM>::output(file_pt,n_plot);}

  protected:

    //??ds just copied directly, think it's ok
    /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
    inline double dshape_and_dtest_eulerian_micromag(const Vector<double> &s,
						     Shape &psi,
						     DShape &dpsidx,
						     Shape &test,
						     DShape &dtestdx) const;

    /// Shape, test functions & derivs. w.r.t. to global coords. at integration point ipt. Return Jacobian.
    inline double dshape_and_dtest_eulerian_at_knot_micromag(const unsigned& ipt,
							     Shape &psi,
							     DShape &dpsidx,
							     Shape &test,
							     DShape &dtestdx) const;
  }; // End of QMagnetostaticBEElement class

}


int main()
{
  return 0;
}
