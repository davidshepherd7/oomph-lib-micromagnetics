#ifndef OOMPH_MICROMAG_TYPES_H
#define OOMPH_MICROMAG_TYPES_H

#include "../../src/generic/Vector.h" // Vectors
#include <utility> // pairs

#include "../../src/generic/Qelements.h" // Q face element geometry


// Header for typedefs needed in multiple places in the library

namespace oomph
{
  // Need some forward decls
  class Mesh;
  class MicromagBEMElementEquations;
  class TimeStepper;
  class FiniteElement;
  class CachingMMInterpolator;


  // Function pointer types
  // ============================================================

  /// General function of space and time which returns a double.
  typedef double (*TimeSpaceToDoubleFctPt)(const double& t,
                                           const Vector<double>&x);

  /// General function of space and time which returns a vector of doubles.
  typedef Vector<double> (*TimeSpaceToDoubleVectFctPt) (const double& t,
                                                        const Vector<double>&x);

  /// General function of time, space and a value vector which returns a
  /// vector of doubles.
  typedef Vector<double> (*TimeSpaceValueToDoubleVectFctPt)
  (const double& t, const Vector<double>&x, const Vector<double>&u);

  /// General function of time and a value which returns a double.
  typedef double (*TimeValueToDoubleFctPt) (const double& t,
                                            const double& u);

  /// Function to calculate magnetostatic field analytically
  typedef Vector<double> (*MagnetostaticFieldFctPt)
  (const double& t, const Vector<double>& x, const Vector<double>& m);

  /// Function type for applied fields
  typedef TimeSpaceToDoubleVectFctPt HAppFctPt;

  /// Function type for initial magnetisation
  typedef SolutionFunctorBase InitialMFct;

  /// Function type for use as general initial condition
  typedef SolutionFunctorBase InitialConditionFct;

  /// \short Function class for functions to integrate over elements.
  class ElementalFunction
  {
  public:

    /// \short Virtual destructor (to stop compiler complaining).
    virtual ~ElementalFunction() {}

    /// \short Function to integrate over the element.
    virtual double call(const GeneralisedElement* ele_pt,
                        CachingMMInterpolator* intp_pt) const = 0;

    // /// \short Helper function to automatically create an interpolator
    // /// object for call.
    // double call(const GeneralisedElement* ele_pt, const Vector<double> &s) const;
  };



  // Factory function pointers
  // ============================================================


  typedef FiniteElement* (*ElementFactoryFctPt)(void);

  /// Function pointer type for function to create a BEM element.
  typedef MicromagBEMElementEquations*
  (*BEMElementFactoryFctPt)(FiniteElement* const, const int&);

  /// Function pointer type for function to create flux meshes
  typedef Mesh* (*FluxMeshFactoryFctPt)(Mesh* bulk_mesh_pt,
                                        const Vector<unsigned> &boundaries);

  /// Function pointer type to create a mesh
  typedef Mesh* (*MeshFactoryFctPt)(const std::string& mesh_name,
                                    int refinement_level,
                                    TimeStepper* time_stepper_pt,
                                    double scaling_factor,
                                    unsigned nnode1d);

  /// Function pointer type to create a nodal quadrature scheme. Mean
  /// element size needed for some schemes so always needs to be passed in
  /// unfortunately.
  typedef Integral* (*NodalQuadratureFactoryFctPt)
    (const FiniteElement* ele_pt, const double& mean_size);



  // Some stl-based data structures
  // ============================================================

  /// Type to hold data on which boundaries of which meshes should be
  /// included in the bem.
  typedef Vector<std::pair<unsigned, const Mesh*> > BemBoundaryData;

  /// Type to hold data on the locations of sharp corners in bem meshes
  typedef Vector<std::pair<Vector<double>, double> > CornerDataInput;


  // QElement face geometries which for some strange reason are not in
  // oomph core.
  // ============================================================

  /// Face geometry for the 1D Qelements: Point elements
  template<unsigned NNODE_1D>
  class FaceGeometry<QElement<1,NNODE_1D> >:
    public PointElement {};

  /// Face geometry for the QElements: The spatial dimension of the face
  /// elements is one lower than that of the bulk element but they have the
  /// same number of points along their 1D edges.
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<QElement<DIM,NNODE_1D> >:
    public QElement<DIM-1,NNODE_1D>
    {};

}

#endif
