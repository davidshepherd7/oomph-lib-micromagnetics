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
    LLGArgs() : mag_params_pt(0)
    {
      renormalisation_handler_pt = 0;
    }

    virtual ~LLGArgs()
    {
      delete mag_params_pt; mag_params_pt = 0;
      delete renormalisation_handler_pt; renormalisation_handler_pt = 0;
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

      specify_command_line_flag("-renormalise", &renormalise_name,
                                "Choose the method to use for renormalisation of m.");
      renormalise_name = "default";

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

      specify_command_line_flag("-use-m-mallinson", &use_m_mallinson,
                                "When using mallinson exact solution use error based on M instead of time, default: -1.");
      use_m_mallinson = -1;


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

      specify_command_line_flag("-relax-m-field", &relax_field_name,
                                "Field to relax magnetisation under before time integration, default 0 = do not relax).");
      relax_m = -1;

      specify_command_line_flag("-relax-m-time", &relax_m_time,
                                "Amount of time to allow m to relax for if -relax-m-field is set, default: 300.");
      relax_m_time = 300;

      specify_command_line_flag("-do-proper-relax", &do_proper_relax,
                                "Relax m in the usual way (small incremental field reductions), default: -1.");
      do_proper_relax = -1;


      specify_command_line_flag("-quadrature", &quadrature_type,
                                "gauss, nodal, rnodal; default is gauss.");
      quadrature_type = "gauss";

      specify_command_line_flag("-wave-solution-c", &wave_solution_c,
                                "Constant for wave exact solution, multiplied by pi to get actual c used. Does nothing if not using it. Default 0.25 (i.e. pi/4).");
      wave_solution_c = 0.25;

      specify_command_line_flag("-read-cached-bem-matrix", &cached_bem_matrix_filename,
                                "Use a bem matrix from a previous run, default: \"\" (don't use any). Be careful with different methods of removing the phi singularity: different pinned nodes = different numbering!");
      cached_bem_matrix_filename = "";

      specify_command_line_flag("-write-cached-bem-matrix", &cached_bem_matrix_filename_out,
                                "Write a bem matrix to disk to reuse later, default: \"\" (don't use any). Be careful with different methods of removing the phi singularity: different pinned nodes = different numbering!");
      cached_bem_matrix_filename_out = "";

      specify_command_line_flag("-fd-jac", &fd_jac,
                                "Use finite differences to calculate the elemental Jacobians, default: -1");
      fd_jac = -1;

      specify_command_line_flag("-pin-boundary-m", &pin_boundary_m,
                                ", default: -1.");
      pin_boundary_m = -1;

      specify_command_line_flag("-gauss-quadrature-energy", &gauss_quadrature_energy,
                                "Force the use of Gaussian quadrature for energy calculations?, default: -1 = use problem's default.");
      gauss_quadrature_energy      = -1;

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

      // Construct a renormalisation handler
      if(renormalise_name == "default")
        {
          if(ts_name == "bdf1" || ts_name == "bdf2" || ts_name == "tr")
            {
              renormalisation_handler_pt = renormalisation_handler_factory("always");
            }
          else
            {
              renormalisation_handler_pt = renormalisation_handler_factory("never");
            }
        }
      else
        {
          renormalisation_handler_pt = renormalisation_handler_factory(renormalise_name);
        }

    }

    virtual void assign_specific_parameters(MyProblem* problem_pt) const
    {
      LLGProblem* llg_pt = checked_dynamic_cast<LLGProblem*>(problem_pt);

      // Set parameters
      llg_pt->set_mag_parameters_pt(mag_params_pt);

      // Set the renormalisation handler
      llg_pt->renormalisation_handler_pt = renormalisation_handler_pt;

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
      if(pin_boundary_m != -1)
        {
          llg_pt->Pin_boundary_m = bool(pin_boundary_m);
          if(bool(pin_boundary_m))
            {
              llg_pt->Boundary_solution_pt = initial_condition_pt;
            }
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

      if(fd_jac != -1)
        {
          llg_pt->Use_fd_jacobian = bool(fd_jac);
        }

      // Set exact solution if we have one: only for spherical
      // nano-particles or when ms is disabled.
      if(((h_app_name == "minus_z")
          && (initial_m_name == "z")
          && (mag_params_pt->gilbert_damping() != 0.0)
          && ((to_lower(ms_method) == "disabled")
              || (mesh_name == "ut_sphere")))
         || mallinson == 1)
        {
          if(use_m_mallinson != -1 && bool(use_m_mallinson))
            {
              llg_pt->Compare_with_mallinson_m = true;
            }
          else
            {
              llg_pt->Compare_with_mallinson = true;
            }
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

      if(relax_field_name != "")
        {
          llg_pt->H_app_relax = Factories::h_app_factory(relax_field_name);
        }
      llg_pt->Relax_m_time = relax_m_time;
      if(do_proper_relax != -1)
        {
          llg_pt->Do_proper_relax = bool(do_proper_relax);
        }


      if(doc_ml_error != -1)
        {
          llg_pt->Doc_m_length_error = bool(doc_ml_error);
        }

      if(gauss_quadrature_energy != -1)
        {
          llg_pt->Force_gaussian_quadrature_in_energy = bool(gauss_quadrature_energy);
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
    std::string relax_field_name;

    int doc_ml_error;


    /// Flag to control renormalisation of |m| after each step.
    std::string renormalise_name;
    RenormalisationHandler* renormalisation_handler_pt;

    int compare_residuals_with_gauss;

    double damping;
    double k1;
    double ms_debug_coeff;
    double h_app_debug_coeff;
    double wave_solution_c;
    double relax_m_time;
    int do_proper_relax;

    int numerical_int_bem;
    int hlib_bem;
    int mallinson;
    int check_angles;
    int disable_magnetostatic_solver_optimistations;
    int relax_m;
    int use_m_mallinson;
    int fd_jac;
    int pin_boundary_m;
    int gauss_quadrature_energy;

    std::string quadrature_type;

    std::string ms_method;


    Vector<Mesh*> phi_1_mesh_pts;
    Vector<Mesh*> phi_mesh_pts;

    BEMElementFactoryFctPt bem_element_factory_fct_pt;

    std::string cached_bem_matrix_filename;
    std::string cached_bem_matrix_filename_out;

  };

} // End of oomph namespace

#endif
