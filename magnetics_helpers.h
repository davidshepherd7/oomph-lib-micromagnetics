#ifndef OOMPH_MAGNETICS_HELPERS_H
#define OOMPH_MAGNETICS_HELPERS_H

/*
description of file goes here
*/

#include "vector_helpers.h"
#include "micromag_types.h"

#include <deque>

#include "../../src/generic/Vector.h"
#include "../../src/generic/oomph_utilities.h"
#include "../../src/generic/oomph_definitions.h"


namespace oomph
{

  class LLGProblem;
  class MagneticParameters;

  /// Generic zero initial condition function
  template<unsigned NVAL>
  inline Vector<double> zero_initial_condition(const double& t,
                                               const Vector<double> &x)
  {
    Vector<double> v(NVAL, 0.0);
    return v;
  }

namespace HApp
{
  using namespace oomph;
  using namespace MathematicalConstants;
  using namespace StringConversion;
  using namespace VectorOps;

  inline Vector<double> zero(const double& t, const Vector<double> &x)
  {
    return zero_initial_condition<3>(t, x);
  }

  inline Vector<double> x(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[0] = 1.1;
    return h;
  }

  inline Vector<double> y(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[1] = 1.1;
    return h;
  }

  inline Vector<double> z(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[2] = 1.1;
    return h;
  }

  inline Vector<double> minus_z(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[2] = -1.1;
    return h;
  }

  inline Vector<double> nanowire(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[2] = -0.3;
    return h;
  }

  inline Vector<double> minus_x(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[0] = -1.1;
    return h;
  }

  inline Vector<double> all_directions(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 1.1);
    return h;
  }

  inline Vector<double> z_oscillating_p20(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[0] = 0.1;
    h[2] = 1.1 * std::sin((MathematicalConstants::Pi/20) * t);
    return h;
  }

  inline Vector<double> mumag4_initial(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    double t0 = 100;

    // Strong, decay "slowly" to zero
    h[0] = std::max(2.5 * (1 - (t/t0)), 0.0);
    h[1] = h[0];
    h[2] = h[0];

    return h;
  }

  inline Vector<double> mumag4_initial_strong(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    double t0 = 100;

    // Strong, decay "slowly" to zero
    h[0] = std::max(10 * (1 - (t/t0)), 0.0);
    h[1] = h[0];
    h[2] = h[0];

    return h;
  }

  inline Vector<double> mumag4_initial_crazy(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    double t0 = 500;

    // Strong, decay "slowly" to zero
    h[0] = std::max(50 * (1 - (t/t0)), 0.0);
    h[1] = h[0];
    h[2] = h[0];

    return h;
  }

  inline Vector<double> mumag4_initial_exponential(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    double t0 = 50;

    // Strong, decay to zero, last bit needs to be slower so use exponential
    h[0] = std::max(50 *std::exp(-t/t0), 0.0);
    h[1] = h[0];
    h[2] = h[0];

    return h;
  }

  inline Vector<double> mumag4_field1(const double& t, const Vector<double> &x)
  {
    double mu0 = 4*Pi*1e-7;
    double Ms = 8e5;
    double a = 1e-3/(mu0*Ms);

    Vector<double> h(3, 0.0);
    h[0] = -24.6*a;
    h[1] = 4.3*a;
    h[2] = 0.0;

    return h;
  }

  inline Vector<double> mumag4_field2(const double& t, const Vector<double> &x)
  {
    double mu0 = 4*Pi*1e-7;
    double Ms = 8e5;
    double a = 1e-3/(mu0*Ms);

    Vector<double> h(3, 0.0);
    h[0] = -35.5*a;
    h[1] = -6.3*a;
    h[2] = 0.0;

    return h;
  }

  inline Vector<double> field_smoother(const double& t, const Vector<double>& h,
                                       const double& t0)
  {
    double a = 1 - std::exp(-t/t0);

    const unsigned ni = h.size();
    Vector<double> h_final(ni, 0.0);
    for(unsigned i=0; i<ni; i++)
      {
        h_final[i] = h[i]*a;
      }

    return h_final;
  }

  // these don't seem to make any difference
  inline Vector<double> smoothed_mumag4_field1(const double& t, const Vector<double> &x)
  {
    Vector<double> h = mumag4_field1(t, x);
    return field_smoother(t, h, 1e-3);
  }

  inline Vector<double> smoothed_mumag4_field2(const double& t, const Vector<double> &x)
  {
    Vector<double> h = mumag4_field2(t, x);
    return field_smoother(t, h, 1e-3);
  }

  inline Vector<double> non_uniform_z_helper
  (const double& t, const Vector<double> &x, double l)
  {
    Vector<double> h(3, 0.0);
    h[0] = 0.1;
    h[2] = 1.1 * std::sin(MathematicalConstants::Pi * x[0] / l);
    return h;
  }

  inline Vector<double> non_uniform_z_5(const double& t, const Vector<double> &x)
  {return non_uniform_z_helper(t, x, 5);}
  inline Vector<double> non_uniform_z_50(const double& t, const Vector<double> &x)
  {return non_uniform_z_helper(t, x, 50);}
  inline Vector<double> non_uniform_z_500(const double& t, const Vector<double> &x)
  {return non_uniform_z_helper(t, x, 500);}

  inline Vector<double> tanhx_minus_z(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[2] = -1.1*tanh(5*(x[0] - 0.5));
    return h;
  }

  inline Vector<double> minus_z_above_x0(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    if(x[0] > 0)
      {
        h[2] = -1.1;
      }
    return h;
  }

  inline Vector<double> tanhx_minus_x(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[0] = -1.1*tanh(5*x[0]);
    return h;
  }

  inline Vector<double> minus_x_above_x0(const double& t, const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    if(x[0] > 0)
      {
        h[0] = -1.1;
      }
    return h;
  }


  inline Vector<double> smooth_start_z(const double& t,
                                       const Vector<double> &x)
  {
    double a = 5;
    double smoothing = (1 - std::exp(-a*t));
    Vector<double> h = HApp::z(t, x);
    for(unsigned j=0; j<3; j++)
      {
        h[j] *= smoothing;
      }
    return h;
  }

  inline Vector<double> smooth_start_minus_z(const double& t,
                                             const Vector<double> &x)
  {
    double a = 5;
    double smoothing = (1 - std::exp(-a*t));
    Vector<double> h = HApp::minus_z(t, x);
    for(unsigned j=0; j<3; j++)
      {
        h[j] *= smoothing;
      }
    return h;
  }


  inline Vector<double> linear_z(const double& t,
                                 const Vector<double> &x)
  {
    Vector<double> h(3, 0.0);
    h[2] = 1.1 * std::max(1.0 - x[0], 0.0);
    return h;
  }
}

namespace InitialM
{

  using namespace oomph;
  using namespace MathematicalConstants;
  using namespace StringConversion;
  using namespace VectorOps;

  // Note that dof storage in micromag elements is actually in the order
  // [phi, phi1, mx, my, mz], so we put mx in element 2 of the vector etc.

  inline Vector<double> x(const double& t, const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);
    m[2] = 1.0;
    m[3] = 0.2;
    normalise(m);
    return m;
  }

  inline Vector<double> y(const double& t, const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);
    m[3] = 1.0;
    m[4] = 0.2;
    normalise(m);
    return m;
  }

  inline Vector<double> z(const double& t, const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);
    m[4] = 1.0;
    m[2] = 0.2;
    normalise(m);
    return m;
  }

  inline Vector<double> ode_z(const double& t, const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);
    m[4] = 1.0;
    m[2] = 0.01;
    normalise(m);
    return m;
  }

  inline Vector<double> xyz(const double& t, const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);
    m[2] = 1.0;
    m[3] = 1.0;
    m[4] = 1.0;
    normalise(m);
    return m;
  }

  inline Vector<double> xy(const double& t, const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);
    m[2] = 1.0;
    m[3] = 1.0;
    normalise(m);
    return m;
  }

  inline Vector<double> exactly_z(const double& t, const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);
    m[4] = 1.0;
    return m;
  }

  inline Vector<double> xz(const double& t, const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);
    m[2] = 1.0;
    m[4] = 1.0;
    normalise(m);
    return m;
  }

  inline Vector<double> smoothly_varying_m_helper(double l, const double& t,
                                                  const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);

    if(x.size() == 1)
      {
        m[2] = sin(x[0]*2*Pi/l);
        m[3] = cos(x[0]*2*Pi/l);
        m[4] = 1.0 - m[2] - m[3];
      }
    else
      {
        m[2] = sin(x[0]*2*Pi/l) + sin(x[1]*2*Pi/l);
        m[3] = cos(x[0]*2*Pi/l) + cos(x[1]*2*Pi/l);
        m[4] = 1.0 - m[2] - m[3];
      }
    normalise(m);

    return m;
  }

  inline Vector<double> alt_smoothly_varying_m_helper(double l, const double& t,
                                                  const Vector<double> &x)
  {
    Vector<double> m(5, 0.0);

    m[2] = sin(x[0]*2*Pi/l)/2 + sin(x[1]*2*Pi/l)/2;
    m[3] = cos(x[0]*2*Pi/l)/2 + cos(x[1]*2*Pi/l)/2;
    m[4] = std::sqrt(1.0 - m[2]*m[2] - m[3]*m[3]); // normalise

#ifdef PARANOID
    if(std::abs(1 - m[2]*m[2] - m[3]*m[3] - m[4]*m[4]) > 1e-14)
      {
        std::string err = "Non unit length initial m!";
        throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

    return m;
  }

  inline Vector<double> smoothly_varying_5(const double& t, const Vector<double> &x)
  {return smoothly_varying_m_helper(5.0, t, x);}
  inline Vector<double> smoothly_varying_50(const double& t, const Vector<double> &x)
  {return smoothly_varying_m_helper(50.0, t, x);}
  inline Vector<double> smoothly_varying_500(const double& t, const Vector<double> &x)
  {return smoothly_varying_m_helper(500.0, t, x);}
  inline Vector<double> smoothly_varying_5000(const double& t, const Vector<double> &x)
  {return smoothly_varying_m_helper(5000.0, t, x);}

  inline Vector<double> alt_smoothly_varying_5(const double& t, const Vector<double> &x)
  {return alt_smoothly_varying_m_helper(5.0, t, x);}
  inline Vector<double> alt_smoothly_varying_50(const double& t, const Vector<double> &x)
  {return alt_smoothly_varying_m_helper(50.0, t, x);}
  inline Vector<double> alt_smoothly_varying_500(const double& t, const Vector<double> &x)
  {return alt_smoothly_varying_m_helper(500.0, t, x);}
  inline Vector<double> alt_smoothly_varying_5000(const double& t, const Vector<double> &x)
  {return alt_smoothly_varying_m_helper(5000.0, t, x);}

  inline Vector<double> smoothly_varying_50_xflipped(const double& t, const Vector<double> &x)
  {
    Vector<double> x2 = x;
    x2[0] = 1 - x2[0];
    return smoothly_varying_m_helper(50, t, x2);
  }


  class LLGWaveSolution : public SolutionFunctorBase
  {
  public:
    /// Constructor
    LLGWaveSolution(const double& c)
    {
      C = c;
      k = 2 * Pi;
    }

    /// Virtual destructor
    virtual ~LLGWaveSolution() {}
    void initialise_from_problem(const Problem* problem_pt);
    Vector<double> operator()(const double& t, const Vector<double>& x) const;

    virtual Vector<double> derivative(const double& t, const Vector<double>& x,
                                      const Vector<double>& u) const
    {
      throw OomphLibError("Not implemented (yet?).", OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    double b(const double& t) const
    {
      // Rescale time because this is a solution to the LL equation
      double t_scaled = t / (1 + damping*damping);
      return dim*k*k*damping*t_scaled;
    }

    double d(const double& t) const
    {
      using namespace std;
      using namespace MathematicalConstants;

      double a = this->C*Pi;
      return sqrt(sin(a)*sin(a)  +  exp(2*b(t))*cos(a)*cos(a));
    }

    double g(const double& t) const
    {
      using namespace std;
      using namespace MathematicalConstants;

      double a = this->C*Pi;
      return (1/damping) * log((d(t)  +  exp(b(t))*cos(a))
                               /(1  +  cos(a)));
    }

    unsigned dim;
    double damping;
    double C;
    double k;
  };


  /// Derivative function for ll equation as an ode, needs to go after
  /// llgode problem class.
  class LLODESolution : public SolutionFunctorBase
  {
  public:
    /// Virtual destructor
    virtual ~LLODESolution()
    {
      magnetic_parameters_pt = 0;
      initial_m_pt = 0;
    }

    /// Just the initial condition actually, no exact solution that can fit
    /// this framework.
    Vector<double> operator()(const double& t, const Vector<double>&x) const override
    {
      Vector<double> full_vector = initial_m_pt->operator()(t, x);

      Vector<double> m(3);
      m[0] = full_vector[2];
      m[1] = full_vector[3];
      m[2] = full_vector[4];

      return m;
    }

    /// Derivative function.
    Vector<double> derivative(const double& t, const Vector<double>& x,
                              const Vector<double>& m) const override;

    void jacobian(const double& t, const Vector<double>& x,
                  const Vector<double>& m,
                  DenseMatrix<double>& jacobian) const override;

    bool have_jacobian() const override {return true;}

    /// Get parameters from problem
    void initialise_from_problem(const Problem* problem_pt) override;

    MagneticParameters* magnetic_parameters_pt;
    InitialMFct* initial_m_pt;
  };

  /// solution from Mallinson, uses LLODESolution for LLG derivative
  /// + Jacobian
  class LLGMallinsonSolution : public SolutionFunctorBase
  {
  public:
    /// Constructor
    LLGMallinsonSolution()
    {
      mag_params_pt = 0;
    }

    /// Virtual destructor
    virtual ~LLGMallinsonSolution() {}

    void initialise_from_problem(const Problem* problem_pt);
    Vector<double> operator()(const double& t, const Vector<double>& x) const;

    Vector<double> derivative(const double& t, const Vector<double>& x,
                              const Vector<double>& u) const
    {
      return llg_ode_solution.derivative(t, x, u);
    }

    void jacobian(const double& t, const Vector<double>& x,
                  const Vector<double>& m,
                  DenseMatrix<double>& jacobian) const override
    {
      llg_ode_solution.jacobian(t, x, m, jacobian);
    }

    Vector<double> initial_m;
    const MagneticParameters* mag_params_pt;
    LLODESolution llg_ode_solution;

  };

}


  namespace MagnetostaticFieldFunctions
  {

    inline Vector<double> sphere(const double& time, const Vector<double>& x,
                                 const Vector<double>& m)
    {
      Vector<double> hms(3, 0.0);

      for(unsigned j=0; j<3; j++)
        {
          hms[j] = -m[j]/3;
        }

      return hms;
    }


    inline MagnetostaticFieldFctPt ms_factory(const std::string& ms_name)
    {
      if(ms_name == "sphere")
        {
          return &MagnetostaticFieldFunctions::sphere;
        }
      else
        {
          throw OomphLibError("Unrecognised initial magnetostatic field name "
                              + ms_name,
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }

  }

class MyProblem;

namespace MManipulation
{

  /// \short Compute the effective damping constant (alpha) for the
  /// previous time step (see Albuquerque2001).
  double alt_effective_damping_used(const LLGProblem& problem,
                                    std::deque<double>& previous_energies);


  /// \short Compute the effective damping constant (alpha) for the
  /// previous time step (see Albuquerque2001).
  double effective_damping_used(const LLGProblem& problem);

  double effective_damping_used_3(const LLGProblem& problem);


  double exchange_energy(const LLGProblem& problem,
                         const Integral* quadrature_pt=0);


  double zeeman_energy(const LLGProblem& problem,
                       const Integral* quadrature_pt=0);

  double crystalline_anisotropy_energy(const LLGProblem& problem,
                                       const Integral* quadrature_pt=0);


  double magnetostatic_energy(const LLGProblem& problem,
                              const Integral* quadrature_pt=0);

  double integral_of_dmdt_squared(const LLGProblem& problem,
                                  const Integral* quadrature_pt=0);

  double dEnergydt(const LLGProblem& problem);

  double alt_dEnergydt(const LLGProblem& problem,
                       const std::deque<double>& previous_energies);

  /// Get a vector of the nodal values of the magnetisation
  Vector<Vector<double> > nodal_magnetisations(const unsigned& t_hist,
                                               const LLGProblem& problem);
  inline Vector<Vector<double> > nodal_magnetisations(const LLGProblem& problem)
  {
    return nodal_magnetisations(0, problem);
  }

  /// Get the average of each magnetisation direction, specify either m
  /// values or problem.
  Vector<double> mean_nodal_magnetisation(const Vector<Vector<double> >& ms);
  inline Vector<double> mean_nodal_magnetisation(const unsigned& t_hist,
                                          const LLGProblem& problem)
  {
    return mean_nodal_magnetisation(nodal_magnetisations(t_hist, problem));
  }
  inline Vector<double> mean_nodal_magnetisation(const LLGProblem& problem)
  {
    return mean_nodal_magnetisation(0, problem);
  }


  /// Get list of |m| errors
  Vector<double> nodal_m_length_errors(const Vector<Vector<double> >& ms);

  /// Integrate a function given by func_pt over every element in a mesh
  /// and return the total. This should probably be in the mesh class but
  /// that's core oomph-lib so I'll leave it here.
  double integrate_over_mesh(const ElementalFunction* func_pt,
                             const Mesh* const mesh_pt,
                             const Integral* quadrature_pt=0);
}

  class TetMeshBase;

  namespace MeshCreationHelpers
  {
    /// The trivial factory function for an ELEMENT.
    template<class ELEMENT>
    inline  FiniteElement* new_element() { return new ELEMENT;}


    void brick2tet(const Mesh& brick_mesh,
                   ElementFactoryFctPt element_factory_fpt,
                   TetMeshBase& out_mesh);

    void make_boundaries_periodic(Mesh* mesh_pt, const unsigned& b1,
                                  const unsigned& b2,
                                  const unsigned& direction);


    /// Helper function to make two boundaries of a mesh periodic. This
    /// version uses a brute force search to find the appropriate nodes to
    /// link together so it should be robust but is O(N^2). In practice it
    /// seems to be un-noticably fast for meshes that I've used so far.
    void slow_make_boundaries_periodic(Mesh* mesh_pt, const unsigned& b1,
                                       const unsigned& b2,
                                       const unsigned& direction);

    /// Move all nodes of a mesh by (x, y, z). ??ds Move inside Mesh?
    void shift_mesh(const double& x, const double& y, const double& z,
                    Mesh* mesh_pt);

    /// Rotate a mesh in the x-y plane
    void rotate_mesh(const double& theta, Mesh* mesh_pt);

    /// Construct an equilateral triangle mesh using a square mesh.
    /// Replaces each square with four triangles:
    ///
    ///     o = node
    ///
    ///       o------------------------------------------------/o
    ///       | \---                                       /--- |
    ///       |     \--                                 /--     |
    ///       |        \---                         /---        |
    ///       |            \---                 /---            |
    ///       |                \--           /--                |
    ///       |                   \---   /---                   |
    ///       |                       o--                       |
    ///       |                   /---   \---                   |
    ///       |                /--           \--                |
    ///       |            /---                 \---            |
    ///       |        /---                         \---        |
    ///       |     /--                                 \--     |
    ///       | /---                                       \--- |
    ///       o------------------------------------------------\o
    ///
    Mesh* equilateral_triangle_mesh(int refinement_level,
                                    TimeStepper* time_stepper_pt,
                                    unsigned nnode1d,
                                    ElementFactoryFctPt element_factory_fpt);

    /// ??ds
    Mesh* union_jack_triangle_mesh(int refinement_level,
                                   TimeStepper* time_stepper_pt,
                                   unsigned nnode1d,
                                   ElementFactoryFctPt element_factory_fpt);

  }

  namespace ErrorNorms
  {
    /// Error norm for wave solution based on phase at x=0
    double wave_phase_error_norm(const LLGProblem& problem);

    /// m_z error norm for wave solution
    double wave_mz_error_norm(const LLGProblem& problem);
  }

}


#endif
