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
    LLGArgs() : h_app_fct_pt(0), mag_params_pt(0) {}

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

      specify_command_line_flag("-damping", &damping);
      damping = -10;

      specify_command_line_flag("-k1", &k1);
      k1 = -10;

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
                                "Set to either 'pin' or 'normalise', default 'pin'.");
      phi1_singularity_method = "pin";
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

      initial_condition_pt = initial_m_factory(initial_m_name);
      h_app_fct_pt = h_app_factory(h_app_name);
      mag_params_pt = magnetic_parameters_factory(mag_params_name);

      if(command_line_flag_has_been_set("-damping"))
        {
          mag_params_pt->Gilbert_damping = damping;
        }

      if(command_line_flag_has_been_set("-k1"))
        {
          mag_params_pt->Anisotropy_coeff = k1;
        }

      // Copy flags into bools in this class
      pin_boundary_m = command_line_flag_has_been_set("-pin-boundary-m");
    }

    virtual void assign_specific_parameters(MyProblem* problem_pt) const
    {
      LLGProblem* llg_pt = checked_dynamic_cast<LLGProblem*>(problem_pt);
      llg_pt->applied_field_fct_pt() = h_app_fct_pt;
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
          llg_pt->Phi_1_singularity_method = phi_1_singularity_handling::pin;
        }
      else if(phi1_singularity_method == "normalise")
        {
          llg_pt->Phi_1_singularity_method = phi_1_singularity_handling::normalise;
        }
      else
        {
          std::string err = "Unrecognised phi1 singularity handling method.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }

    HApp::HAppFctPt h_app_fct_pt;
    MagneticParameters* mag_params_pt;

    // Strings for input to factory functions
    std::string initial_m_name;
    std::string h_app_name;
    std::string mag_params_name;
    std::string llg_prec_name;
    std::string llg_sub_prec_name;
    std::string phi1_singularity_method;


    /// Flag to control renormalisation of |m| after each step. -1 =
    /// default for timestepper, 0 = off, 1 = on.
    int renormalise;

    double damping;
    double k1;

    int numerical_int_bem;
    int hlib_bem;
    int mallinson;
    int check_angles;
    int disable_magnetostatic_solver_optimistations;

    std::string ms_method;

    bool pin_boundary_m;

    Vector<Mesh*> phi_1_mesh_pts;
    Vector<Mesh*> phi_mesh_pts;

    BEMElementFactoryFctPt bem_element_factory_fct_pt;

  };

} // End of oomph namespace

#endif
