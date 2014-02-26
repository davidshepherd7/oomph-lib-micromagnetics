#ifndef OOMPH_MICROMAG_TYPES_H
#define OOMPH_MICROMAG_TYPES_H

#include "../../src/generic/Vector.h" // Vectors
#include <utility> // pairs


// Header for typedefs needed in multiple places in the library

namespace oomph
{
  // Need some forward decls
  class Mesh;
  class MicromagBEMElementEquations;
  class TimeStepper;
  class FiniteElement;



  // Function pointer types
  // ============================================================

  /// General function of space and time which returns a double.
  typedef double (*TimeSpaceToDoubleFctPt)(const double& t,
                                           const Vector<double>&x);

  /// General function of space and time which returns a vector of doubles.
  typedef Vector<double> (*TimeSpaceToDoubleVectFctPt) (const double& t,
                                                        const Vector<double>&x);

  /// General function of time and a value which returns a double.
  typedef double (*TimeValueToDoubleFctPt) (const double& t,
                                            const double& u);

  /// Function to calculate magnetostatic field analytically
  typedef Vector<double> (*MagnetostaticFieldFctPt)
  (const double& t, const Vector<double>& x, const Vector<double>& m);

  /// Function type for applied fields
  typedef TimeSpaceToDoubleVectFctPt HAppFctPt;

  /// Function type for initial magnetisation
  typedef TimeSpaceToDoubleVectFctPt InitialMFctPt;

  /// Function type for use as general initial condition. Slight overhead
  /// of vectors is worth it even in cases with one dof/space dimensions
  /// for the generality. No overhead for returning a vector due to return
  /// value optimisation.
  typedef TimeSpaceToDoubleVectFctPt InitialConditionFctPt;


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



  // Some stl-based data structures
  // ============================================================

  /// Type to hold data on which boundaries of which meshes should be
  /// included in the bem.
  typedef Vector<std::pair<unsigned, const Mesh*> > BemBoundaryData;

  /// Type to hold data on the locations of sharp corners in bem meshes
  typedef Vector<std::pair<Vector<double>, double> > CornerDataInput;

}

#endif
