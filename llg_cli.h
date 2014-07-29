#ifndef OOMPH_LLG_CLI_H
#define OOMPH_LLG_CLI_H

#include "my_cli.h"
#include "oomph_factories.h"

#include "llg_factories.h"
#include "llg_problem.h"

namespace oomph
{

  /// Command line args processing class for llg problems.
  class LLGArgs : public MyCliArgs
  {
  public:
    /// Constructor: Initialise pointers to null.
    LLGArgs() : mag_params_pt(0) {}

    virtual ~LLGArgs()
    {
      delete mag_params_pt; mag_params_pt = 0;
    }

    virtual void set_flags()
    {
      MyCliArgs::set_flags();

      specify_command_line_flag("-initial-m", &initial_m_name);
      initial_m_name = "z";

      specify_command_line_flag("-h-app", &h_app_name);
      h_app_name = "minus_z";

      specify_command_line_flag("-mag-params", &mag_params_name);
      mag_params_name = "simple-llg";

      specify_command_line_flag("-renormalise", &renormalise);
      renormalise = -1;

      specify_command_line_flag("-compare-residuals-with-gauss", &compare_residuals_with_gauss,
                                "Compare residual calculations with those given by high order Gauss? default: -1.");
      compare_residuals_with_gauss = -1;

      specify_command_line_flag("-damping", &damping);
      damping = -10;

      specify_command_line_flag("-k1", &k1);
      k1 = -10;

      specify_command_line_flag("-ms-debug-coeff", &ms_debug_coeff,
                                "Debugging multiplier for magnetostatic field strength.");
      ms_debug_coeff = -10;

      specify_command_line_flag("-h-app-debug-coeff", &h_app_debug_coeff,
                                "Debugging multiplier for applied field strength.");
      h_app_debug_coeff = -10;

      // Flags automatically default to false
      specify_command_line_flag("-pin-boundary-m");

      specify_command_line_flag("-hlib-bem", &hlib_bem);
      hlib_bem = -1;

      specify_command_line_flag("-numerical-int-bem", &numerical_int_bem);
      numerical_int_bem = -1;

      specify_command_line_flag("-mallinson", &mallinson);
      mallinson = -1;

      specify_command_line_flag("-ms-method", &ms_method);
      ms_method = "implicit";

      specify_command_line_flag("-check-angles", &check_angles);
      check_angles = -1;

      specify_command_line_flag("-doc-ml-error", &doc_ml_error,
                                "Write out spatial values of error in |m|, default: -1.");
      doc_ml_error = -1;

      specify_command_line_flag("-llg-prec", &llg_prec_name,
                                "Set preconditioner for llg block, only used if overall preconditioner is a magnetostatics block preconditioner");
      llg_prec_name = "";

      specify_command_line_flag("-llg-sub-prec", &llg_sub_prec_name,
                                "Set preconditioner for llg sub block (ony if llg block prec in use)");
      llg_sub_prec_name = "";

      specify_command_line_flag("-disable-magnetostatic-solver-optimistations",
                                &disable_magnetostatic_solver_optimistations);
      disable_magnetostatic_solver_optimistations = -1;

      specify_command_line_flag("-phi1-singularity-method",
                                &phi1_singularity_method,
                                "Set to either 'pin_bulk', 'pin', 'pin_boundary', 'normalise' or 'nothing', default 'pin_bulk'.");
      phi1_singularity_method = "pin_bulk";

      specify_command_line_flag("-relax-m", &relax_m,
                                "Should the magnetisation be relaxed before starting time integration? (-1/0/1, -1 lets the class keep its default, default -1).");
      relax_m = -1;

      specify_command_line_flag("-quadrature", &quadrature_type,
                                "gauss, nodal, rnodal; default is gauss.");
      quadrature_type = "gauss";

      specify_command_line_flag("-wave-solution-c", &wave_solution_c,
                                "Constant for wave exact solution, multiplied by pi to get actual c used. Does nothing if not using it. Default 0.25 (i.e. pi/4).");
      wave_solution_c = 0.25;
    }

    bool is_decoupled(const std::string& ms_method) const
    {
      return to_lower(ms_method) == "decoupled"
        || to_lower(ms_method) == "decoupled-no-extrapolation";
    }

    virtual Preconditioner* preconditioner_factory(const std::string& name) const
    {
      // Call the factory
      return Factories::micromag_preconditioner_factory(name, llg_prec_name,
                                                        llg_sub_prec_name);
    }

    virtual void run_factories()
    {
      using namespace Factories;

      mesh_factory_pt = &llg_mesh_factory;


      MyCliArgs::run_factories();

      initial_m_name = to_lower(initial_m_name);
      h_app_name = to_lower(h_app_name);
      mag_params_name = to_lower(mag_params_name);

      initial_condition_pt = initial_m_factory(initial_m_name, wave_solution_c);

      // Build magnetic parameters object
      mag_params_pt = magnetic_parameters_factory(mag_params_name);

      if(command_line_flag_has_been_set("-damping"))
        {
          mag_params_pt->Gilbert_damping = damping;
        }

      if(command_line_flag_has_been_set("-k1"))
        {
          mag_params_pt->Anisotropy_coeff = k1;
        }

      if(command_line_flag_has_been_set("-ms-debug-coeff"))
        {
          mag_params_pt->Magnetostatic_debug_coeff = ms_debug_coeff;
        }

      if(command_line_flag_has_been_set("-h-app-debug-coeff"))
        {
          mag_params_pt->Applied_field_debug_coeff = h_app_debug_coeff;
        }

      HApp::HAppFctPt h_app_fct_pt = h_app_factory(h_app_name);
      mag_params_pt->Applied_field_fct_pt = h_app_fct_pt;


      // Copy flags into bools in this class
      pin_boundary_m = command_line_flag_has_been_set("-pin-boundary-m");
    }

    virtual void assign_specific_parameters(MyProblem* problem_pt) const
    {
      LLGProblem* llg_pt = checked_dynamic_cast<LLGProblem*>(problem_pt);

      // Set parameters
      llg_pt->set_mag_parameters_pt(mag_params_pt);


      if(renormalise == -1)
        {
          llg_pt->Renormalise_each_time_step =
            (ts_name == "bdf1" || ts_name == "bdf2" || ts_name == "tr");
        }
      else
        {
          llg_pt->Renormalise_each_time_step = renormalise;
        }

      if(compare_residuals_with_gauss != -1)
        {
          llg_pt->Compare_with_gauss_quadrature = bool(compare_residuals_with_gauss);
        }

      /// Pick the function which will be used by each element to create an
      /// quadrature scheme.
      llg_pt->Nodal_quadrature_factory_fpt =
        Factories::nodal_quadrature_factory_factory(quadrature_type);


      // Dirichlet boundries, just use same function for b.c. as initial
      // cond.
      llg_pt->Pin_boundary_m = pin_boundary_m;
      if(pin_boundary_m)
        {
          llg_pt->Boundary_solution_pt = initial_condition_pt;
        }

      if(to_lower(ms_method) == "implicit")
        {
          llg_pt->Decoupled_ms = false;
          llg_pt->Disable_ms = false;
          llg_pt->Extrapolate_decoupled_ms = false;
        }
      else if(to_lower(ms_method) == "decoupled-no-extrapolation")
        {
          llg_pt->Decoupled_ms = true;
          llg_pt->Disable_ms = false;
          llg_pt->Extrapolate_decoupled_ms = false;
        }
      else if(to_lower(ms_method) == "decoupled")
        {
          llg_pt->Decoupled_ms = true;
          llg_pt->Disable_ms = false;
          llg_pt->Extrapolate_decoupled_ms = true;
        }
      else if(to_lower(ms_method) == "disabled")
        {
          llg_pt->Decoupled_ms = false;
          llg_pt->Disable_ms = true;
          llg_pt->Extrapolate_decoupled_ms = false;
        }
      else
        {
          llg_pt->Decoupled_ms = false;
          llg_pt->Disable_ms = true;
          llg_pt->Extrapolate_decoupled_ms = false;
          llg_pt->Analytic_ms_fct_pt =
            MagnetostaticFieldFunctions::ms_factory(to_lower(ms_method));
        }

      // ??ds this should maybe be a general one?
      llg_pt->Use_fd_jacobian = use_fd_jacobian;

      // Set exact solution if we have one: only for spherical
      // nano-particles or when ms is disabled.
      if(((h_app_name == "minus_z")
          && (initial_m_name == "z")
          && (mag_params_pt->gilbert_damping() != 0.0)
          && ((to_lower(ms_method) == "disabled")
              || (mesh_name == "ut_sphere")))
         || mallinson == 1)
        {
          llg_pt->Compare_with_mallinson = true;
        }

      if(check_angles != -1)
        {
          llg_pt->Check_angles = bool(check_angles);
        }

      if(disable_magnetostatic_solver_optimistations != -1)
        {
          llg_pt->Disable_magnetostatic_solver_optimistations
            = bool(disable_magnetostatic_solver_optimistations);
        }

      llg_pt->Bem_element_factory_pt = bem_element_factory_fct_pt;


      // Assign phi1 singularity handler method
      if(phi1_singularity_method == "pin")
        {
          llg_pt->Phi_1_singularity_method = phi_1_singularity_handling::pin_any;
        }
      else if(phi1_singularity_method == "pin_bulk")
        {
          llg_pt->Phi_1_singularity_method = phi_1_singularity_handling::pin_bulk;
        }
      else if(phi1_singularity_method == "pin_boundary")
        {
          llg_pt->Phi_1_singularity_method = phi_1_singularity_handling::pin_boundary;
        }
      else if(phi1_singularity_method == "normalise")
        {
          llg_pt->Phi_1_singularity_method = phi_1_singularity_handling::normalise;
        }
      else if(phi1_singularity_method == "nothing")
        {
          llg_pt->Phi_1_singularity_method = phi_1_singularity_handling::nothing;
        }
      else if(phi1_singularity_method == "jacobian")
        {
          llg_pt->Phi_1_singularity_method = phi_1_singularity_handling::jacobian;
        }
      else
        {
          std::string err = "Unrecognised phi1 singularity handling method.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      if(relax_m == 0)
        {
          llg_pt->Relax_magnetisation = false;
        }
      else if(relax_m == 1)
        {
          llg_pt->Relax_magnetisation = true;
        }
      // else do nothing


      if(doc_ml_error != -1)
        {
          llg_pt->Doc_m_length_error = bool(doc_ml_error);
        }

    }

    MagneticParameters* mag_params_pt;

    // Strings for input to factory functions
    std::string initial_m_name;
    std::string h_app_name;
    std::string mag_params_name;
    std::string llg_prec_name;
    std::string llg_sub_prec_name;
    std::string phi1_singularity_method;

    int doc_ml_error;


    /// Flag to control renormalisation of |m| after each step. -1 =
    /// default for timestepper, 0 = off, 1 = on.
    int renormalise;

    int compare_residuals_with_gauss;

    double damping;
    double k1;
    double ms_debug_coeff;
    double h_app_debug_coeff;
    double wave_solution_c;

    int numerical_int_bem;
    int hlib_bem;
    int mallinson;
    int check_angles;
    int disable_magnetostatic_solver_optimistations;
    int relax_m;

    std::string quadrature_type;

    std::string ms_method;

    bool pin_boundary_m;

    Vector<Mesh*> phi_1_mesh_pts;
    Vector<Mesh*> phi_mesh_pts;

    BEMElementFactoryFctPt bem_element_factory_fct_pt;

  };

} // End of oomph namespace

#endif
