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

  /// Function class for exact solutions/initial conditions/boundary
  /// conditions. This is needed so that we can have solutions that depend
  /// on problem parameters with resorting to global variables.
  class SolutionFunctor
  {
    public:
    SolutionFunctor() {}

    SolutionFunctor(TimeSpaceToDoubleVectFctPt solution_fpt)
    {
      Solution_fpt = solution_fpt;
    }

    virtual ~SolutionFunctor() {}

    SolutionFunctor(const SolutionFunctor& that)
    {
      Solution_fpt = that.Solution_fpt;
    }

    void operator=(const SolutionFunctor& that)
    {
      this->Solution_fpt = that.Solution_fpt;
    }

    /// Call the function.
    Vector<double> operator()(const double& t, const Vector<double>&x) const
    {
#ifdef PARANOID
      if(Solution_fpt == 0)
        {
          std::string err = "Solution_fpt is null!";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return Solution_fpt(t, x);
    }

    /// Overload to grab data from the problem.
    virtual void initialise_from_problem(const Problem* problem_pt) {}

    private:
    TimeSpaceToDoubleVectFctPt Solution_fpt;
  };

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
