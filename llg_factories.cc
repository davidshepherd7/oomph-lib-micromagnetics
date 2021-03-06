#include "llg_factories.h"

#include "llg_problem.h"
#include "llg_preconditioners.h"
#include "boundary_element_handler.h"
#include "pinned_boundary_element_handler.h"

#include "nodal_quadrature.h"

#include "new_interpolators.h"

// Meshes for mesh factory
#include "../../src/meshes/simple_rectangular_quadmesh.h"
#include "../../src/meshes/rectangular_quadmesh.h"
#include "../../src/meshes/one_d_mesh.h"
#include "../../src/meshes/simple_rectangular_tri_mesh.h"
#include "../../src/meshes/simple_cubic_tet_mesh.h"
#include "../../src/meshes/simple_cubic_mesh.h"
#include "../../src/meshes/tetgen_mesh.h"
#include "../../src/meshes/triangle_mesh.h"
#include "./single_element_mesh.h"
#include "./simpler_cubic_mesh.h"

namespace oomph
{
  namespace Factories
  {
    BoundaryElementHandlerBase* bem_handler_factory
    (const Vector<Mesh*>& output_mesh_pts,
     const CornerDataInput* input_corner_data_pt,
     int hierarchical_bem,
     bool disable_corner_angles,
     int numerical_int_bem,
     bool allow_pinned_boundary_values,
     const std::string& bem_matrix_filename_in,
     const std::string& bem_matrix_filename_out
     )
    {
      BoundaryElementHandlerBase* bem_handler_pt = 0;
      if(allow_pinned_boundary_values)
        {
          bem_handler_pt = new PinnedBoundaryElementHandler;
        }
      else
        {
          bem_handler_pt = new BoundaryElementHandler;
        }

      // Figure out what defaults to use if any bool-like options are -1
      // ============================================================


      // Get the first finite element we can find in the a bulk mesh on
      // which we are going to construct our bem mesh. Assuimg that all
      // elements in all the meshes are the same type...
      FiniteElement* bulk_fe_pt = output_mesh_pts[0]->finite_element_pt(0);

      // Use H-lib if possible (have it and surface mesh is triangular)
      if(hierarchical_bem == -1)
        {
#ifdef OOMPH_HAS_HLIB
          if((bulk_fe_pt->nodal_dimension() == 3)
             && (bulk_fe_pt->nnode_1d() == 2)
             && (bulk_fe_pt->nnode() == 4))
            {
              hierarchical_bem = true;
            }
          else
            {
              hierarchical_bem = false;
            }
#else
          hierarchical_bem = false;
#endif
        }

      // Use analytical integation if possible, numerical otherwise
      if(numerical_int_bem == -1)
        {
          if((bulk_fe_pt->nodal_dimension() == 3)
             && (bulk_fe_pt->nnode_1d() == 2)
             && (bulk_fe_pt->nnode() == 4))
            {
              numerical_int_bem = false;
            }
          else
            {
              numerical_int_bem = true;
            }
        }

      // Next assign the parameters
      // ============================================================

      // Check that we can do hierarchical bem, if so set the parameter.
      if(hierarchical_bem)
        {
#ifndef OOMPH_HAS_HLIB
          std::string err = "Hlib library required for hierarchical bem matrix";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
#endif
        }
      bem_handler_pt->Hierarchical_bem = hierarchical_bem;

      // Figure out which element type we should use in the bem mesh
      // (based on the element type used in the bulk mesh) and store the
      // function needed to create them.
      bem_handler_pt->Bem_element_factory_fpt = Factories::
        bem_element_factory_factory(bulk_fe_pt);

      // Create an integration scheme
      bem_handler_pt->integration_scheme_pt() = Factories::
        variable_order_integrator_factory(bulk_fe_pt);

      // Figure out if we are doing phi/phi1 in separate meshes (and
      // problems) or all in one by checking what type of element we
      // have. Set the indicies accordingly.
      MicromagEquations* mele_pt = dynamic_cast<MicromagEquations*>
        (output_mesh_pts[0]->element_pt(0));
      TFPoissonEquations* pele_pt = dynamic_cast<TFPoissonEquations*>
        (output_mesh_pts[0]->element_pt(0));
      if(mele_pt != 0)
        {
          // Fully implicit/all in one mesh
          bem_handler_pt->set_input_index(mele_pt->phi_1_index_micromag());
          bem_handler_pt->set_output_index(mele_pt->phi_index_micromag());
        }
      else if(pele_pt != 0)
        {
          // Fully implicit/all in one mesh
          bem_handler_pt->set_input_index(pele_pt->u_index_poisson());
          bem_handler_pt->set_output_index(pele_pt->u_index_poisson());
        }

      // Add all boundaries of all meshes to bem boundary list. Not likely
      // to want to make use of bem on only some boundaries anytime soon so
      // just add them all. If you want it on only some boundaries write a
      // different factory (sorry).
      BemBoundaryData bem_boundaries;
      for(unsigned msh=0, nmsh=output_mesh_pts.size(); msh<nmsh; msh++)
        {
          Mesh* mesh_pt = output_mesh_pts[msh];
          for(unsigned b=0, nb=mesh_pt->nboundary(); b<nb; b++)
            {
              bem_boundaries.push_back(std::make_pair(b, mesh_pt));
            }
        }

      // Copy in the list of boundaries to operate on
      bem_handler_pt->Bem_boundaries = bem_boundaries;

      // Set cached matrix filename
      bem_handler_pt->Bem_matrix_filename_in = bem_matrix_filename_in;
      bem_handler_pt->Bem_matrix_filename_out = bem_matrix_filename_out;

      // Set debug parameters
      bem_handler_pt->Debug_disable_corner_contributions = disable_corner_angles;
      bem_handler_pt->Numerical_int_bem = numerical_int_bem;


      // Now build it
      if(input_corner_data_pt == 0)
        {
          CornerDataInput dummy;
          bem_handler_pt->build(dummy);
        }
      else
        {
          bem_handler_pt->build(*input_corner_data_pt);
        }

      return bem_handler_pt;
    }


    Preconditioner* llg_sub_preconditioner_factory(const std::string& llg_sub_prec)
    {
      Preconditioner* llg_sub_prec_pt = 0;
      if(llg_sub_prec == "block")
        {
          LLGSubBlockPreconditioner* _llg_sub_prec_pt = new LLGSubBlockPreconditioner;
          _llg_sub_prec_pt->build();
          llg_sub_prec_pt = _llg_sub_prec_pt;
        }
      else
        {
          llg_sub_prec_pt = preconditioner_factory(llg_sub_prec);
        }

      return llg_sub_prec_pt;
    }


    Preconditioner* llg_preconditioner_factory(const std::string& llg_prec,
                                               const std::string& llg_sub_prec)
    {
      Preconditioner* llg_prec_pt = 0;
      if(llg_prec == "block")
        {
          LLGBlockPreconditioner* _llg_prec_pt = new LLGBlockPreconditioner;

          Preconditioner* llg_sub_prec_pt = llg_sub_preconditioner_factory(llg_sub_prec);

          Vector<unsigned> master_to_subs_map(2);
          master_to_subs_map[0] = 0; // mx
          master_to_subs_map[1] = 1; // my

          llg_sub_prec_pt->turn_into_subsidiary_block_preconditioner(_llg_prec_pt, master_to_subs_map);
          _llg_prec_pt->J_aabb_prec_pt = llg_sub_prec_pt;

          _llg_prec_pt->build();

          llg_prec_pt = _llg_prec_pt;
        }
      else
        {
          llg_prec_pt = preconditioner_factory(llg_prec);
        }

      return llg_prec_pt;
    }


    Preconditioner* micromag_preconditioner_factory(const std::string& ms_prec,
                                                    const std::string& llg_prec,
                                                    const std::string& llg_sub_prec)
    {
      // Magnetostatics prec
      // ============================================================

      Preconditioner* ms_prec_pt = 0;
      // First check if we want a sum of matrice preconditioner:
      if(has_prefix("som-main-", ms_prec))
        {
          Preconditioner* ul_prec = micromag_preconditioner_factory
            (rest_of_name("som-main-", ms_prec), llg_prec, llg_sub_prec);

          MainMatrixOnlyPreconditioner* mm_prec_pt = new MainMatrixOnlyPreconditioner;
          mm_prec_pt->set_underlying_prec_pt(ul_prec);

          ms_prec_pt = mm_prec_pt;
        }
      // Make a preconditioner which only acts on the main matrix and
      // diagonals of added matrices of a sum of matrices.
      else if(has_prefix("som-maindiag-", ms_prec))
        {
          Preconditioner* ul_prec = micromag_preconditioner_factory
            (rest_of_name("som-maindiag-", ms_prec), llg_prec, llg_sub_prec);

          MainMatrixAndDiagsPreconditioner* mm_prec_pt
            = new MainMatrixAndDiagsPreconditioner;
          mm_prec_pt->set_underlying_prec_pt(ul_prec);

          ms_prec_pt = mm_prec_pt;
        }
      else if(ms_prec == "dummy")
        {
          // Make preconditioners
          DummyPinnedMsPreconditioner* _ms_prec_pt = new DummyPinnedMsPreconditioner;
          Preconditioner* llg_prec_pt = llg_preconditioner_factory(llg_prec,
                                                                   llg_sub_prec);

          // Set up master/subsidiary links
          Vector<unsigned> micromag_to_llg_block_map(3);
          micromag_to_llg_block_map[0] = 2; // mx
          micromag_to_llg_block_map[1] = 3; // my
          micromag_to_llg_block_map[2] = 4; // mz

          llg_prec_pt->turn_into_subsidiary_block_preconditioner
            (_ms_prec_pt, micromag_to_llg_block_map);
          _ms_prec_pt->Real_preconditioner_pt = llg_prec_pt;

          ms_prec_pt = _ms_prec_pt;
        }

      else if(ms_prec == "block")
        {
          // Make preconditioners
          MagnetostaticsBlockPreconditioner* _ms_prec_pt = new MagnetostaticsBlockPreconditioner;
          Preconditioner* llg_prec_pt = llg_preconditioner_factory(llg_prec,
                                                                   llg_sub_prec);

          // Set up master/subsidiary links (note: different to above!)
          Vector<unsigned> micromag_to_llg_block_map(3);
          micromag_to_llg_block_map[0] = 2; // mx
          micromag_to_llg_block_map[1] = 3; // my
          micromag_to_llg_block_map[2] = 4; // mz

          llg_prec_pt->turn_into_subsidiary_block_preconditioner
            (_ms_prec_pt, micromag_to_llg_block_map);
          _ms_prec_pt->Llg_preconditioner_pt = llg_prec_pt;

          _ms_prec_pt->build();

          ms_prec_pt = _ms_prec_pt;
        }

      else if(ms_prec == "block-drop-p")
        {
          // same as "block" except that we drop the P block (rather than
          // Q) of the bulk-phi+llg block to get it triangular.

          // Make preconditioners
          MagnetostaticsBlockPreconditioner* _ms_prec_pt = new MagnetostaticsBlockPreconditioner;
          Preconditioner* llg_prec_pt = llg_preconditioner_factory(llg_prec,
                                                                   llg_sub_prec);

          // Set up master/subsidiary links
          Vector<unsigned> micromag_to_llg_block_map(3);
          micromag_to_llg_block_map[0] = 2; // mx
          micromag_to_llg_block_map[1] = 3; // my
          micromag_to_llg_block_map[2] = 4; // mz

          llg_prec_pt->turn_into_subsidiary_block_preconditioner
            (_ms_prec_pt, micromag_to_llg_block_map);
          _ms_prec_pt->Llg_preconditioner_pt = llg_prec_pt;

          _ms_prec_pt->Drop_P = true;

          _ms_prec_pt->build();

          ms_prec_pt = _ms_prec_pt;
        }
      else
        {
          ms_prec_pt = preconditioner_factory(ms_prec);
        }

      return ms_prec_pt;
    }


    Vector<unsigned> dof_to_block_factory(const std::string& _name)
    {
      const std::string name = to_lower(_name);

      // Make an element to look up indicies from
      TMicromagElement<2,2> dummy_ele;

      const unsigned ndof = dummy_ele.ndof_types(); //??ds unsafe?
      Vector<unsigned> dof_to_block(ndof);

      if(name == "none")
        {
          // identity mapping
          for(unsigned j=0; j<ndof; j++)
            {
              dof_to_block[j] = j;
            }
        }

      // All m values in one block, others left alone.
      // [0, 1, 2, 2, 2, 3, 4]
      else if(name == "group-m")
        {
          unsigned k = 0;

          // Normal phi/phi1
          dof_to_block[dummy_ele.phi_index_micromag()] = k++;
          dof_to_block[dummy_ele.phi_1_index_micromag()] = k++;

          // m all in one block
          for(unsigned j=0; j<3; j++)
            {
              int index = dummy_ele.m_index_micromag(j);
              dof_to_block[index] = k;
            }
          k++;

          // boundary phi/phi1 ??ds assume they are at the end...
          dof_to_block[5] = k++;
          dof_to_block[6] = k++;
        }
      else if(name == "group-m-phi-phi-boundary")
        {
          unsigned k = 0;

          // All phi into one block
          dof_to_block[dummy_ele.phi_index_micromag()] = k;
          dof_to_block[5] = k;
          k++;

          // Simiarly for phi1
          dof_to_block[dummy_ele.phi_1_index_micromag()] = k;
          dof_to_block[6] = k;
          k++;

          // m all in one block
          for(unsigned j=0; j<3; j++)
            {
              int index = dummy_ele.m_index_micromag(j);
              dof_to_block[index] = k;
            }
          k++;
        }
      else
        {
          std::string err = "Unrecognised blocking name";
          err += name;
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }

      return dof_to_block;
    }

    /// \short Make a mesh as specified by an input argument. Refined
    /// according to the given refinement level (in some way appropriate
    /// for that mesh type). Assumption: this will be passed into a
    /// problem, which will delete the pointer when it's done.
    Mesh* llg_mesh_factory(const std::string& _mesh_name,
                           int refinement_level,
                           TimeStepper* time_stepper_pt,
                           double scaling_factor,
                           double rotate_xy_angle,
                           unsigned nnode1d)
    {
      // Ignore case in mesh names
      const std::string mesh_name = to_lower(_mesh_name);

      // Refinement always roughly the same for structured meshes
      unsigned nx = 5 * std::pow(2, refinement_level-1);

      // Make the mesh and store a pointer to it
      Mesh* mesh_pt = 0;
      if(mesh_name == "sq_square" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new SimpleRectangularQuadMesh<QMicromagElement<2,2> >
            (nx, nx, lx, lx, time_stepper_pt);
        }
      else if(mesh_name == "set_square" && nnode1d == 2)
        {
          ElementFactoryFctPt f_pt =
            &general_element_factory<TMicromagElement<2, 2> >;

          mesh_pt = MeshCreationHelpers::equilateral_triangle_mesh
            (refinement_level, time_stepper_pt, nnode1d, f_pt);

          mesh_pt->setup_boundary_element_info();

          // Turn off triangle refinement dump stuff (breaks Micromag
          // elements).
          checked_dynamic_cast<TriangleMeshBase*>(mesh_pt)->
            disable_triangulateio_restart();
        }
    else if(mesh_name == "st_union_jack_square" && nnode1d == 2)
      {
        ElementFactoryFctPt f_pt =
          &general_element_factory<TMicromagElement<2, 2> >;

        mesh_pt = MeshCreationHelpers::union_jack_triangle_mesh
          (refinement_level, time_stepper_pt, nnode1d, f_pt);

        mesh_pt->setup_boundary_element_info();

        // Turn off triangle refinement dump stuff (breaks Micromag
        // elements).
        checked_dynamic_cast<TriangleMeshBase*>(mesh_pt)->
          disable_triangulateio_restart();
      }
    else if(mesh_name == "set_square_periodic" && nnode1d == 2)
      {
        mesh_pt = llg_mesh_factory("set_square", refinement_level,
                                   time_stepper_pt, 1.0, 0.0, 2);

        // Link boundary 0 to boundary 2 and boundary 1 to boundary 3
        MeshCreationHelpers::slow_make_boundaries_periodic(mesh_pt, 1, 3, 0); // x
        MeshCreationHelpers::slow_make_boundaries_periodic(mesh_pt, 0, 2, 1); // y
      }
      else if(mesh_name == "sq_square" && nnode1d == 3)
        {
          double lx = 1.0;
          mesh_pt = new SimpleRectangularQuadMesh<QMicromagElement<2,3> >
            (nx, nx, lx, lx, time_stepper_pt);
        }
      else if(mesh_name == "sq_square_xstretch" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new SimpleRectangularQuadMesh<QMicromagElement<2,2> >
            (nx, nx, lx, lx, time_stepper_pt);

          // Stretch in the x direction (for testing jacobian of
          // transformation stuff).
          const unsigned n_node = mesh_pt->nnode();
          for(unsigned nd=0; nd<n_node; nd++)
            {
              Node* nd_pt = mesh_pt->node_pt(nd);
              nd_pt->x(0) *= 2;
            }
        }
      else if(mesh_name == "sq_square_periodic" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new RectangularQuadMesh<QMicromagElement<2,2> >
            (nx, nx, lx, lx, time_stepper_pt);

          // Link boundary 0 to boundary 2 and boundary 1 to boundary 3
          MeshCreationHelpers::slow_make_boundaries_periodic(mesh_pt, 1, 3, 0); // x
          MeshCreationHelpers::slow_make_boundaries_periodic(mesh_pt, 0, 2, 1); // y
        }
      else if(mesh_name == "st_square_periodic" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new SimpleRectangularTriMesh<TMicromagElement<2,2> >
            (nx, nx, lx, lx, time_stepper_pt);

          // Link boundary 0 to boundary 2 and boundary 1 to boundary 3
          MeshCreationHelpers::slow_make_boundaries_periodic(mesh_pt, 1, 3, 0); // x
          MeshCreationHelpers::slow_make_boundaries_periodic(mesh_pt, 0, 2, 1); // y
        }
      else if(mesh_name == "sq_line" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new OneDMesh<QMicromagElement<1,2> >
            (nx, lx, time_stepper_pt);
          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "st_line" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new OneDMesh<TMicromagElement<1,2> >
            (nx, lx, time_stepper_pt);
          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "sq_cube" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new SimplerCubicMesh<QMicromagElement<3,2> >
            (nx, nx, nx, lx, lx, lx, time_stepper_pt);
        }
      else if(mesh_name == "st_cube" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new SimpleCubicTetMesh<TMicromagElement<3,2> >
            (nx, nx, nx, lx, lx, lx, time_stepper_pt);

          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "sq_cube_periodic" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new SimplerCubicMesh<QMicromagElement<3,2> >
            (nx, nx, nx, lx, lx, lx, time_stepper_pt);

          MeshCreationHelpers::slow_make_boundaries_periodic(mesh_pt, 0, 5, 2); // x
          MeshCreationHelpers::slow_make_boundaries_periodic(mesh_pt, 1, 3, 1); // y
          MeshCreationHelpers::slow_make_boundaries_periodic(mesh_pt, 2, 4, 0); // z
        }
      else if(mesh_name == "sq_line_periodic" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new OneDMesh<QMicromagElement<1,2> >
            (nx, lx, time_stepper_pt);

          MeshCreationHelpers::slow_make_boundaries_periodic(mesh_pt, 0, 1, 0); // x
          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "st_square" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new SimpleRectangularTriMesh<TMicromagElement<2,2> >
            (nx, nx, lx, lx, time_stepper_pt);

          mesh_pt->setup_boundary_element_info();

          // Turn off triangle refinement dump stuff (breaks Micromag
          // elements).
          checked_dynamic_cast<TriangleMeshBase*>(mesh_pt)->
            disable_triangulateio_restart();
        }
      else if(mesh_name == "q_single_element" && nnode1d == 2)
        {
          mesh_pt = new QSingleElementMesh<QMicromagElement<2,2> >(time_stepper_pt);
        }
      else if(mesh_name == "t_single_element" && nnode1d == 2)
        {
          mesh_pt = new TSingleElementMesh<TMicromagElement<2,2> >(time_stepper_pt);

          mesh_pt->setup_boundary_element_info();

          // Turn off triangle refinement dump stuff (breaks Micromag
          // elements).
          checked_dynamic_cast<TriangleMeshBase*>(mesh_pt)->
            disable_triangulateio_restart();
        }
      else if(mesh_name == "ut_square" && nnode1d == 2)
        {
          mesh_pt = new TriangleMesh<TMicromagElement<2, 2> >
            ("./meshes/square." + to_string(refinement_level) + ".node",
             "./meshes/square." + to_string(refinement_level) + ".ele",
             "./meshes/square." + to_string(refinement_level) + ".poly",
             time_stepper_pt);

          // Turn off triangle refinement dump stuff (breaks Micromag
          // elements).
          checked_dynamic_cast<TriangleMeshBase*>(mesh_pt)->
            disable_triangulateio_restart();
        }
      else if(mesh_name == "st_cubeoid" && nnode1d == 2)
        {
          // nmag cubeoid
          double lx = 1, ly = lx, lz = 3*lx;
          unsigned ny = nx, nz = std::ceil(lz/lx) * nx;
          mesh_pt = new SimpleCubicTetMesh<TMicromagElement<3, 2> >
            (nx, ny, nz, lx, ly, lz, time_stepper_pt);

          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "sqt_cubeoid" && nnode1d == 2)
        {
          Mesh* qmesh_pt = llg_mesh_factory("sq_cubeoid", refinement_level,
                                            time_stepper_pt,
                                            1, rotate_xy_angle, nnode1d);
          //??ds memory leak, fix? Can't delete this mesh or nodes will
          //go...

          TetMeshBase* tmesh_pt = new TetMeshBase;
          ElementFactoryFctPt factory_fpt =
            MeshCreationHelpers::new_element<TMicromagElement<3, 2> >;

          MeshCreationHelpers::brick2tet(*qmesh_pt, factory_fpt, *tmesh_pt);

          mesh_pt = tmesh_pt;
        }
      else if(mesh_name == "ut_cubeoid" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TMicromagElement<3, 2> >
            ("./meshes/cubeoid." + to_string(refinement_level) + ".node",
             "./meshes/cubeoid." + to_string(refinement_level) + ".ele",
             "./meshes/cubeoid." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else if(mesh_name == "ut_mumag4" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TMicromagElement<3, 2> >
            ("./meshes/mumag4." + to_string(refinement_level) + ".node",
             "./meshes/mumag4." + to_string(refinement_level) + ".ele",
             "./meshes/mumag4." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else if(mesh_name == "st_mumag4" && nnode1d == 2)
        {
          mesh_pt = new SimpleCubicTetMesh<TMicromagElement<3, 2> >
            (5*nx, std::ceil(1.25*nx), 1, 500, 125, 3, time_stepper_pt);

          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "sq_mumag4" && nnode1d == 2)
        {
          unsigned this_nx = refinement_level;

          mesh_pt = new SimplerCubicMesh<QMicromagElement<3, 2> >
            (5*this_nx, std::ceil(1.25*this_nx), 1, 500, 125, 3, time_stepper_pt);

          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "sq_mumag4_normalised" && nnode1d == 2)
        {
          unsigned this_nx = refinement_level;

          mesh_pt = new SimplerCubicMesh<QMicromagElement<3, 2> >
            (5*this_nx, std::ceil(1.25*this_nx), 1, 500, 125, 3, time_stepper_pt);

          mesh_pt->scale_mesh(1/5.6858);

          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "sq_mumag4_normalised_zref" && nnode1d == 2)
        {
          unsigned this_nx = refinement_level;

          mesh_pt = new SimplerCubicMesh<QMicromagElement<3, 2> >
            (5*this_nx, std::ceil(1.25*this_nx), 2, 500, 125, 3, time_stepper_pt);
          // two elements thick instead of 1

          mesh_pt->scale_mesh(1/5.6858);

          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "sq_mumag4_2d" && nnode1d == 2)
        {
          unsigned this_nx = refinement_level;

          mesh_pt = new SimpleRectangularQuadMesh<QMicromagElement<2, 2> >
            (5*this_nx, std::ceil(1.25*this_nx), 500, 125, time_stepper_pt);

          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "sq_mumag4+" && nnode1d == 2)
        {
          unsigned this_nx = refinement_level;

          RefineableMeshBase* ref_mesh_pt
            = new SimplerCubicMesh<QMicromagElement<3, 2> >
            (5*this_nx, std::ceil(1.25*this_nx), 1, 500, 125, 3, time_stepper_pt);

          // Create an "error" vector such that only one element is refined
          double max_error =  ref_mesh_pt->max_permitted_error();
          Vector<double> fake_errors(ref_mesh_pt->nelement(), max_error/2);
          fake_errors[0] = max_error*2;

          // And adapt it
          ref_mesh_pt->adapt(fake_errors);

          ref_mesh_pt->setup_boundary_element_info();

          mesh_pt = ref_mesh_pt;
        }
      else if(mesh_name == "sqt_mumag4" && nnode1d == 2)
        {
          Mesh* qmesh_pt = llg_mesh_factory("sq_mumag4", refinement_level,
                                            time_stepper_pt,
                                            1, rotate_xy_angle, nnode1d);

          // Convert to tet mesh
          TetMeshBase* tmesh_pt = new TetMeshBase;
          ElementFactoryFctPt factory_fpt =
            MeshCreationHelpers::new_element<TMicromagElement<3, 2> >;
          MeshCreationHelpers::brick2tet(*qmesh_pt, factory_fpt, *tmesh_pt);

          // delete the Q mesh without deleting the nodes/elements
          qmesh_pt->flush_element_and_node_storage();
          delete qmesh_pt; qmesh_pt = 0;

          mesh_pt = tmesh_pt;
        }
      else if(mesh_name == "sq_cubeoid" && nnode1d == 2)
        {
          double lx = 1, ly = lx, lz = 3*lx;
          mesh_pt = new SimplerCubicMesh<QMicromagElement<3, 2> >
            (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "ut_sphere" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TMicromagElement<3, 2> >
            ("./meshes/sphere." + to_string(refinement_level) + ".node",
             "./meshes/sphere." + to_string(refinement_level) + ".ele",
             "./meshes/sphere." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else if(mesh_name == "ut_cylinder" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TMicromagElement<3, 2> >
            ("./meshes/cylinder25_40." + to_string(refinement_level) + ".node",
             "./meshes/cylinder25_40." + to_string(refinement_level) + ".ele",
             "./meshes/cylinder25_40." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else if(mesh_name == "sq_annular" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new SimpleRectangularQuadMesh<QMicromagElement<2,2> >
            (nx, nx, lx, lx, time_stepper_pt);

          const double r_min = 0.5, r_max = 1.5, phi_min = 0, phi_max = 180;

          // Warp the domain as in oomph-lib website example
          const unsigned n_node = mesh_pt->nnode();
          for(unsigned nd=0; nd<n_node; nd++)
            {
              Node* nd_pt = mesh_pt->node_pt(nd);

              // Get the nodal coordinates
              const double x_old = nd_pt->x(0);
              const double y_old = nd_pt->x(1);

              // Map from the old x/y to the new r/phi:
              const double r = r_min+(r_max-r_min)*x_old;
              const double phi = (phi_min + (phi_max-phi_min)*y_old)*Pi/180.0;

              // Set new nodal coordinates
              nd_pt->x(0)=r*std::cos(phi);
              nd_pt->x(1)=r*std::sin(phi);
            }
        }
      else
        {
          throw OomphLibError("Unrecognised mesh name " + mesh_name,
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      // Check nnode1d is correct
      if(mesh_pt->finite_element_pt(0)->nnode_1d() != nnode1d)
        {
          std::string err = "Incorrect nnode1d";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      // rotate the mesh if requested
      if(rotate_xy_angle != 0.0)
        {
          MeshCreationHelpers::rotate_mesh(rotate_xy_angle, mesh_pt);
        }

      // Scale the mesh as requested
      mesh_pt->scale_mesh(scaling_factor);

      // This should go inside an element factory but our meshes don't
      // allow that :(
      for(unsigned ele=0, nele=mesh_pt->nelement(); ele<nele; ele++)
        {
          MicromagEquations* ele_pt = checked_dynamic_cast<MicromagEquations*>
            (mesh_pt->element_pt(ele));
          ele_pt->Ms_calc_pt = new ImplicitMagnetostaticsCalculator;
        }

      // Done: pass out the mesh pointer
      return mesh_pt;
    }

    /// Pick and create a residual calculator to use
    LLGResidualCalculator* residual_calculator_factory(const std::string& residual)
    {
      if(residual == "llg")
        return new LLGResidualCalculator(true);
      else if(residual == "ll")
        return new LLGResidualCalculator(false);
      else
        throw OomphLibError("Unrecognised residual "+residual,
                            OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
    }

    /// \short Create a variable order quadrature object based on the
    /// dimension and shape of the element.
    Integral* variable_order_integrator_factory(const unsigned& dim,
                                                const unsigned& nnode_1d,
                                                const ElementGeometry::ElementGeometry& geom)
    {
      using namespace ElementGeometry;

      //??ds nnnode1d doesn't matter right??


      if((dim == 2) && (geom == T))
        {
          return new TVariableOrderGaussLegendre<1>;
        }
      else if((dim == 2) && (geom == Q))
        {
          return new QVariableOrderGaussLegendre<1>;
        }
      else if((dim == 3) && (geom == T))
        {
          return new TVariableOrderGaussLegendre<2>;
        }
      else if((dim == 3) && (geom == Q))
        {
          return new QVariableOrderGaussLegendre<2>;
        }
      else
        {
          std::string err("Cannot determine element type.\n");
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }


    /// \short Create a variable order quadrature object based on the
    /// dimension and shape of the element.
    Integral* variable_order_integrator_factory(const FiniteElement* const el_pt)
    {
      return variable_order_integrator_factory(el_pt->dim(),
                                               el_pt->nnode_1d(),
                                               el_pt->element_geometry());
    }


    /// \short Return a function which will create the appropriate BEM face
    /// element for the bulk element pointer given (should work for a
    /// pointer to any bulk element type i.e., field or llg).
    BEMElementFactoryFctPt bem_element_factory_factory
    (const FiniteElement* bulk_ele_pt)
    {
      if(dynamic_cast<const TElement<1,2>*>(bulk_ele_pt) != 0)
        {
          return &bem_element_factory<TMicromagBEMElement<1,2> >;
        }
      else if(dynamic_cast<const TElement<2, 2>*>(bulk_ele_pt) != 0)
        {
          return &bem_element_factory<TMicromagBEMElement<2,2> >;
        }
      else if(dynamic_cast<const TElement<3, 2>*>(bulk_ele_pt) != 0)
        {
          return &bem_element_factory<TMicromagBEMElement<3,2> >;
        }

      else if(dynamic_cast<const QElement<1,2>*>(bulk_ele_pt) != 0)
        {
          return &bem_element_factory<QMicromagBEMElement<1,2> >;
        }
      else if(dynamic_cast<const QElement<2,2>*>(bulk_ele_pt) != 0)
        {
          return &bem_element_factory<QMicromagBEMElement<2,2> >;
        }
      else if(dynamic_cast<const QElement<3,2>*>(bulk_ele_pt) != 0)
        {
          return &bem_element_factory<QMicromagBEMElement<3,2> >;
        }

      else if(dynamic_cast<const TElement<1,3>*>(bulk_ele_pt) != 0)
        {
          return &bem_element_factory<TMicromagBEMElement<1,3> >;
        }
      else if(dynamic_cast<const TElement<2, 3>*>(bulk_ele_pt) != 0)
        {
          return &bem_element_factory<TMicromagBEMElement<2,3> >;
        }
      else if(dynamic_cast<const TElement<3, 3>*>(bulk_ele_pt) != 0)
        {
          return &bem_element_factory<TMicromagBEMElement<3,3> >;
        }

      else if(dynamic_cast<const QElement<1,3>*>(bulk_ele_pt) != 0)
        {
          return &bem_element_factory<QMicromagBEMElement<1,3> >;
        }
      else if(dynamic_cast<const QElement<2,3>*>(bulk_ele_pt) != 0)
        {
          return &bem_element_factory<QMicromagBEMElement<2,3> >;
        }
      else if(dynamic_cast<const QElement<3,3>*>(bulk_ele_pt) != 0)
        {
          return &bem_element_factory<QMicromagBEMElement<3,3> >;
        }

      else
        {
          throw OomphLibError("Unrecognised element type",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }

    /// \short Return a factory function which will create the appropriate
    /// "flux mesh" for the bulk element pointer given.
    FluxMeshFactoryFctPt
    mm_flux_mesh_factory_factory(const FiniteElement* bulk_ele_pt)
    {
      if(dynamic_cast<const TMicromagElement<1, 2>*>(bulk_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory
            <MicromagFluxElement<TMicromagElement<1, 2> > >;
        }
      else if(dynamic_cast<const TMicromagElement<2, 2>*>(bulk_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory
            <MicromagFluxElement<TMicromagElement<2, 2> > >;
        }
      else if(dynamic_cast<const TMicromagElement<3, 2>*>(bulk_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory
            <MicromagFluxElement<TMicromagElement<3, 2> > >;
        }

      else if(dynamic_cast<const QMicromagElement<1,2>*>(bulk_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory
            <MicromagFluxElement<QMicromagElement<1,2> > >;
        }
      else if(dynamic_cast<const QMicromagElement<2,2>*>(bulk_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory
            <MicromagFluxElement<QMicromagElement<2,2> > >;
        }
      else if(dynamic_cast<const QMicromagElement<3,2>*>(bulk_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory
            <MicromagFluxElement<QMicromagElement<3,2> > >;
        }

      // nnode1d =3
      else if(dynamic_cast<const QMicromagElement<1,3>*>(bulk_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory
            <MicromagFluxElement<QMicromagElement<1,3> > >;
        }
      else if(dynamic_cast<const QMicromagElement<2,3>*>(bulk_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory
            <MicromagFluxElement<QMicromagElement<2,3> > >;
        }
      else if(dynamic_cast<const QMicromagElement<3,3>*>(bulk_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory
            <MicromagFluxElement<QMicromagElement<3,3> > >;
        }

      else
        {
          throw OomphLibError("Unrecognised element type",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

    }


    HAppFctPt h_app_factory(const std::string& field_name)
    {
      if(field_name == "zero")
        {
          return &HApp::zero;
        }
      else if(field_name == "x")
        {
          return &HApp::x;
        }
      else if(field_name == "y")
        {
          return &HApp::y;
        }
      else if(field_name == "z")
        {
          return &HApp::z;
        }
      else if(field_name == "minus_z")
        {
          return &HApp::minus_z;
        }
      else if(field_name == "nanowire")
        {
          return &HApp::nanowire;
        }
      else if(field_name == "minus_x")
        {
          return &HApp::minus_x;
        }
      else if(field_name == "all_directions")
        {
          return &HApp::all_directions;
        }
      else if(field_name == "z_oscillating_p20")
        {
          return &HApp::z_oscillating_p20;
        }
      else if(field_name == "non_uniform_z_5")
        {
          return &HApp::non_uniform_z_5;
        }
      else if(field_name == "non_uniform_z_50")
        {
          return &HApp::non_uniform_z_50;
        }
      else if(field_name == "non_uniform_z_500")
        {
          return &HApp::non_uniform_z_500;
        }
      else if(field_name == "tanhx_minus_z")
        {
          return &HApp::tanhx_minus_z;
        }
      else if(field_name == "minus_z_above_x0")
        {
          return &HApp::minus_z_above_x0;
        }
      else if(field_name == "tanhx_minus_x")
        {
          return &HApp::tanhx_minus_x;
        }
      else if(field_name == "minus_x_above_x0")
        {
          return &HApp::minus_x_above_x0;
        }
      else if(field_name == "mumag4_initial")
        {
          return &HApp::mumag4_initial;
        }
      else if(field_name == "mumag4_initial_strong")
        {
          return &HApp::mumag4_initial_strong;
        }
      else if(field_name == "mumag4_initial_crazy")
        {
          return &HApp::mumag4_initial_crazy;
        }
      else if(field_name == "mumag4_initial_exponential")
        {
          return &HApp::mumag4_initial_exponential;
        }
      else if(field_name == "mumag4_initial_proper_relax")
        {
          return &HApp::mumag4_initial_proper_relax;
        }
      else if(field_name == "mumag4_field1")
        {
          return &HApp::mumag4_field1;
        }
      else if(field_name == "mumag4_field2")
        {
          return &HApp::mumag4_field2;
        }
      else if(field_name == "smoothed_mumag4_field1")
        {
          return &HApp::smoothed_mumag4_field1;
        }
      else if(field_name == "smoothed_mumag4_field2")
        {
          return &HApp::smoothed_mumag4_field2;
        }
      else if(field_name == "smooth_start_minus_z")
        {
          return &HApp::smooth_start_minus_z;
        }
      else if(field_name == "smooth_start_z")
        {
          return &HApp::smooth_start_z;
        }
      else if(field_name == "linear_z")
        {
          return &HApp::linear_z;
        }
      else
        {
          throw OomphLibError("Unrecognised field name " + field_name,
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }

    InitialMFct* initial_m_factory(const std::string& m_name,
                                   const double& wave_solution_c)
    {
      if(m_name == "periodic_exact")
        {
          return new InitialM::LLGWaveSolution(wave_solution_c);
        }

      TimeSpaceToDoubleVectFctPt fpt;

      if(m_name == "x")
        {
          fpt = &InitialM::x;
        }
      else if(m_name == "y")
        {
          fpt = &InitialM::y;
        }
      else if(m_name == "z")
        {
          fpt = &InitialM::z;
        }
      else if(m_name == "ode_z")
        {
          fpt = &InitialM::ode_z;
        }
      else if(m_name == "xyz")
        {
          fpt = &InitialM::xyz;
        }
      else if(m_name == "xy")
        {
          fpt = &InitialM::xy;
        }
      else if(m_name == "xz")
        {
          fpt = &InitialM::xz;
        }
      else if(m_name == "exactly_z")
        {
          fpt = &InitialM::exactly_z;
        }
      else if(m_name == "smoothly_varying_5")
        {
          fpt = &InitialM::smoothly_varying_5;
        }
      else if(m_name == "smoothly_varying_50")
        {
          fpt = &InitialM::smoothly_varying_50;
        }
      else if(m_name == "smoothly_varying_50_xflipped")
        {
          fpt = &InitialM::smoothly_varying_50_xflipped;
        }
      else if(m_name == "smoothly_varying_500")
        {
          fpt = &InitialM::smoothly_varying_500;
        }
      else if(m_name == "smoothly_varying_5000")
        {
          fpt = &InitialM::smoothly_varying_5000;
        }
      else if(m_name == "alt_smoothly_varying_5")
        {
          fpt = &InitialM::alt_smoothly_varying_5;
        }
      else if(m_name == "alt_smoothly_varying_50")
        {
          fpt = &InitialM::alt_smoothly_varying_50;
        }
      else if(m_name == "alt_smoothly_varying_500")
        {
          fpt = &InitialM::alt_smoothly_varying_500;
        }
      else if(m_name == "alt_smoothly_varying_5000")
        {
          fpt = &InitialM::alt_smoothly_varying_5000;
        }
      else
        {
          throw OomphLibError("Unrecognised initial m name " + m_name,
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      return new SolutionFunctor(fpt);
    }


    Integral* rescaled_nodal_quadrature_factory(const FiniteElement* ele_pt,
                                                const double& mean_size)
    {
      return new RescaledNodalQuadrature(ele_pt, mean_size);
    }


    template<class T>
    Integral* nodal_quadrature_factory(const FiniteElement* ele_pt,
                                       const double& mean_size)
    {
      return new T(ele_pt);
    }


    NodalQuadratureFactoryFctPt
    nodal_quadrature_factory_factory(const std::string& label)
    {
      NodalQuadratureFactoryFctPt nqf_pt = 0;
      if(label == "gauss")
        {
          nqf_pt = &Factories::no_nodal_quadrature_scheme_factory;
        }
      else if(label == "nodal")
        {
          nqf_pt = &Factories::nodal_quadrature_factory<NodalQuadrature>;
        }
      else if(label == "lnodal")
        {
          nqf_pt = &Factories::nodal_quadrature_factory<LocalNodalQuadrature>;
        }
      else if(label == "gnodal")
        {
          nqf_pt = &Factories::nodal_quadrature_factory<GlobalNodalQuadrature>;
        }
      else if(label == "uwnodal")
        {
          nqf_pt = &Factories::nodal_quadrature_factory<UnitWeightNodalQuadrature>;
        }
      else if(label == "rnodal")
        {
          nqf_pt = &Factories::rescaled_nodal_quadrature_factory;
        }
      else
        {
          std::string err = "Unrecognisied integration scheme type" + label;
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      return nqf_pt;
    }

    MagneticParameters* magnetic_parameters_factory(const std::string& name)
    {
      using namespace MagParameterHelpers;
      using namespace StringConversion;

      MagneticParameters* parameters_pt = 0;

      // Set properties as used in mumag standard problem #4
      if(to_lower(name) == "mumag4")
        {
          parameters_pt = magnetic_parameters_factory("simple-llg");
          parameters_pt->Exchange_coeff = Astar_to_A(1.3e-11, 1e-9, 8e5);
          parameters_pt->Gilbert_damping = 0.02;
        }

      // Set parameters to remove as many coefficients as possible from the
      // LLG.
      else if(to_lower(name) == "simple-llg")
        {
          parameters_pt = new MagneticParameters;
        }

      // Set parameters as many coefficients as possible from the LLG. Set
      // damping = 0.5, lex = 1, hk = 1.
      else if(to_lower(name) == "simple-llg-anisotropy")
        {
          parameters_pt = magnetic_parameters_factory("simple-llg");
          parameters_pt->Anisotropy_coeff = 0.5; // This gives hk = 1
        }

      else
        {
          std::string error_msg = "Unrecognised name: "
            + to_lower(name);
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      return parameters_pt;
    }

    CachingMMArrayInterpolator* mm_array_interpolator_factory
    (const FiniteElement* ele_pt, const unsigned& time_step_index,
     const TimeStepper* alternative_ts_pt)
    {
      // Choose time stepper
      const TimeStepper* ts_pt = 0;
      if(alternative_ts_pt != 0)
        {
          ts_pt = alternative_ts_pt;
        }
      else
        {
          ts_pt = ele_pt->node_pt(0)->time_stepper_pt();
        }

      // Create a caching interpolator
      CachingMMArrayInterpolator* cintp_pt = new CachingMMArrayInterpolator;

      // Choose underlying interpolator type
      if(ts_pt->nprev_values_for_value_at_evaluation_time() != 1
         && time_step_index == 0)
        {
          cintp_pt->Intp_pt = new RawTimeInterpolatedValueInterpolator;
        }
      else if(time_step_index != 0)
        {
          HistoryInterpolator* intp_pt = new HistoryInterpolator;
          intp_pt->Time_index = time_step_index;
          cintp_pt->Intp_pt = intp_pt;
        }
      else
        {
          cintp_pt->Intp_pt = new RawInterpolator;
        }

      // Set pointers in interpolator
      cintp_pt->Intp_pt->This_element = ele_pt;
      cintp_pt->Intp_pt->Ts_pt = ts_pt;

      return cintp_pt;
    }

    CachingMMInterpolator* mm_interpolator_factory
    (const FiniteElement* ele_pt, const unsigned& time_step_index,
     const TimeStepper* alternative_ts_pt)
    {
      // Choose time stepper
      const TimeStepper* ts_pt = 0;
      if(alternative_ts_pt != 0)
        {
          ts_pt = alternative_ts_pt;
        }
      else
        {
          ts_pt = ele_pt->node_pt(0)->time_stepper_pt();
        }

      // Create a caching interpolator
      CachingMMInterpolator* cintp_pt = new CachingMMInterpolator;

      // Choose underlying interpolator type
      if(ts_pt->nprev_values_for_value_at_evaluation_time() != 1
         && time_step_index == 0)
        {
          cintp_pt->Intp_pt = new RawTimeInterpolatedValueInterpolator;
        }
      else if(time_step_index != 0)
        {
          HistoryInterpolator* intp_pt = new HistoryInterpolator;
          intp_pt->Time_index = time_step_index;
          cintp_pt->Intp_pt = intp_pt;
        }
      else
        {
          cintp_pt->Intp_pt = new RawInterpolator;
        }

      // Set pointers in interpolator
      cintp_pt->Intp_pt->This_element = ele_pt;
      cintp_pt->Intp_pt->Ts_pt = ts_pt;

      return cintp_pt;
    }

  } // end of Factories
}
