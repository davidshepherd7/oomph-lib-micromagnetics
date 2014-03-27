#include "llg_factories.h"

#include "llg_problem.h"
#include "boundary_element_handler.h"


// Meshes for mesh factory
#include "../../src/meshes/simple_rectangular_quadmesh.h"
#include "../../src/meshes/rectangular_quadmesh.h"
#include "../../src/meshes/one_d_mesh.h"
#include "../../src/meshes/simple_rectangular_tri_mesh.h"
#include "../../src/meshes/simple_cubic_tet_mesh.h"
#include "../../src/meshes/simple_cubic_mesh.h"
#include "../../src/meshes/tetgen_mesh.h"
#include "../../src/meshes/triangle_mesh.h"
#include "./multi_mesh.h"
#include "./single_element_mesh.h"

namespace oomph
{
  namespace Factories
  {


    void bem_handler_factory(BoundaryElementHandler& new_bem_handler,
                             const Vector<Mesh*>& output_mesh_pts,
                             const CornerDataInput* input_corner_data_pt,
                             int hierarchical_bem,
                             bool disable_corner_angles,
                             int numerical_int_bem)
    {
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
      new_bem_handler.Hierarchical_bem = hierarchical_bem;

      // Figure out which element type we should use in the bem mesh
      // (based on the element type used in the bulk mesh) and store the
      // function needed to create them.
      new_bem_handler.Bem_element_factory_fpt = Factories::
        bem_element_factory_factory(bulk_fe_pt);

      // Create an integration scheme
      new_bem_handler.integration_scheme_pt() = Factories::
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
          new_bem_handler.set_input_index(mele_pt->phi_1_index_micromag());
          new_bem_handler.set_output_index(mele_pt->phi_index_micromag());
        }
      else if(pele_pt != 0)
        {
          // Fully implicit/all in one mesh
          new_bem_handler.set_input_index(pele_pt->u_index_poisson());
          new_bem_handler.set_output_index(pele_pt->u_index_poisson());
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
      new_bem_handler.Bem_boundaries = bem_boundaries;

      // Set debug parameters
      new_bem_handler.Debug_disable_corner_contributions = disable_corner_angles;
      new_bem_handler.Numerical_int_bem = numerical_int_bem;

      // Now build it
      if(input_corner_data_pt == 0)
        {
          CornerDataInput dummy;
          new_bem_handler.build(dummy);
        }
      else
        {
          new_bem_handler.build(*input_corner_data_pt);
        }
    }


    BoundaryElementHandler* bem_handler_factory
    (const Vector<Mesh*>& output_mesh_pts,
     const CornerDataInput* input_corner_data_pt,
     int hierarchical_bem,
     bool disable_corner_angles,
     int numerical_int_bem)
    {
      // Create with new, fill in with factory
      BoundaryElementHandler* bem_handler_pt = new BoundaryElementHandler;
      bem_handler_factory(*bem_handler_pt, output_mesh_pts,
                          input_corner_data_pt,
                          hierarchical_bem, disable_corner_angles,
                          numerical_int_bem);
      return bem_handler_pt;
    }



    Preconditioner* block_llg_factory(const std::string &_prec_name)
    {
      // Parse the parameter string
      Vector<std::string> parameters = split_string(to_lower(_prec_name), '-');

#ifdef PARANOID
      if(parameters[0] != "blockllg")
        {
          std::string error_msg = _prec_name + " is not a block llg preconditioner";
          error_msg += "(should begin with blockllg-).";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      if(parameters.size() != 5)
        {
          std::string error_msg = "Not enough parameters in llg block string.";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      std::cout << "block llg parameters = " << parameters << std::endl;

      // Pick the basic block structure
      GeneralPurposeBlockPreconditioner<CRDoubleMatrix>* bp_pt = 0;
      if(parameters[1] == "blockexact")
        {
          bp_pt = new ExactBlockPreconditioner<CRDoubleMatrix>;
        }
      else if(parameters[1] == "uppertriangular")
        {
          BlockTriangularPreconditioner<CRDoubleMatrix>* tri_prec_pt
            = new BlockTriangularPreconditioner<CRDoubleMatrix>;
          tri_prec_pt->upper_triangular();
          bp_pt = tri_prec_pt;
        }
      else if(parameters[1] == "lowertriangular")
        {
          BlockTriangularPreconditioner<CRDoubleMatrix>* tri_prec_pt
            = new BlockTriangularPreconditioner<CRDoubleMatrix>;
          tri_prec_pt->lower_triangular();
          bp_pt = tri_prec_pt;
        }
      else if(parameters[1] == "blockdiagonal")
        {
          bp_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
        }
      else
        {
          throw OomphLibError("Unrecognised block structure setting "
                              + parameters[1], OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }


      // Create + set reordering of blocks: the two dofs given (x = m_x
      // etc.) are put into the first block, the other magnetisation dof is
      // put into the second block. Poisson dofs are in their own blocks
      // (empty for now).
      Vector<unsigned> a(5);
      if(parameters[3] == "xy")
        {
          a[0] = 2; a[1] = 3; // poisson blocks
          a[2] = 0; a[3] = 0; // first m block
          a[4] = 1; // second m block
        }
      else if(parameters[3] == "xz")
        {
          a[0] = 2; a[1] = 3; // poisson blocks
          a[2] = 0; a[4] = 0; // first m block
          a[3] = 1; // second m block
        }
      else if(parameters[3] == "yz")
        {
          a[0] = 2; a[1] = 3; // poisson blocks
          a[4] = 0; a[4] = 0; // first m block
          a[2] = 1; // second m block
        }
      else
        {
          throw OomphLibError("Unrecognised block swapping setting "
                              + parameters[3], OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      bp_pt->set_dof_to_block_map(a);


      // Pick the solver to use on the individual blocks (or on the entire
      // thing if we are using "blockexact").
      if(parameters[4] == "exact")
        {
          // Do nothing--default GeneralPurposeBlockPreconditioner solver
          // is SuperLU.
        }
      // else if(parameters[4] == "amg")
      //   {
      //     bp_pt->set_subsidiary_preconditioner_function();
      //   }
      else
        {
          throw OomphLibError("Unrecognised block structure setting "
                              + parameters[4], OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }


      // Now create the upper-left block sub preconditioner.
      // ============================================================
      GeneralPurposeBlockPreconditioner<CRDoubleMatrix>* sub_bp_pt = 0;
      if(parameters[2] == "blockexact")
        {
          sub_bp_pt = new ExactBlockPreconditioner<CRDoubleMatrix>;
        }
      else if(parameters[2] == "uppertriangular")
        {
          BlockTriangularPreconditioner<CRDoubleMatrix>* tri_prec_pt
            = new BlockTriangularPreconditioner<CRDoubleMatrix>;
          tri_prec_pt->upper_triangular();
          sub_bp_pt = tri_prec_pt;
        }
      else if(parameters[2] == "lowertriangular")
        {
          BlockTriangularPreconditioner<CRDoubleMatrix>* tri_prec_pt
            = new BlockTriangularPreconditioner<CRDoubleMatrix>;
          tri_prec_pt->lower_triangular();
          sub_bp_pt = tri_prec_pt;
        }
      else if(parameters[2] == "blockdiagonal")
        {
          sub_bp_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
        }
      else if(parameters[2] == "blockantidiagonal")
        {
          sub_bp_pt = new BlockAntiDiagonalPreconditioner<CRDoubleMatrix>;
        }
      else
        {
          std::string err = "Unrecognised sub-preconditioner block structure";
          err += " setting " + parameters[2];
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      // The blocks we want to use it on are just the ones which are in the
      // zeroth block of the master. Find out which ones these are.
      Vector<unsigned> sub_bp_mapping;
      for(unsigned j=0; j<a.size(); j++)
        {if(a[j] == 0) sub_bp_mapping.push_back(j);}

      // And set the master pointer along with this mapping
      sub_bp_pt->turn_into_subsidiary_block_preconditioner(bp_pt,
                                                           sub_bp_mapping);

      // Provide our subsidiary preconditioner to the master for use on
      // block 0, which is the first m block.
      bp_pt->set_subsidiary_preconditioner_pt(sub_bp_pt, 0);

      return bp_pt;
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
      else if(mesh_name == "sq_square_periodic" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new RectangularQuadMesh<QMicromagElement<2,2> >
            (nx, nx, lx, lx, time_stepper_pt);

          // Link boundary 0 to boundary 2 and boundary 1 to boundary 3
          MeshCreationHelpers::make_boundaries_periodic(mesh_pt, 1, 3, 0); // x
          MeshCreationHelpers::make_boundaries_periodic(mesh_pt, 0, 2, 1); // y
        }
      else if(mesh_name == "sq_line" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new OneDMesh<QMicromagElement<1,2> >
            (nx, lx, time_stepper_pt);
          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "sq_line_periodic" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new OneDMesh<QMicromagElement<1,2> >
            (nx, lx, time_stepper_pt);

          MeshCreationHelpers::make_boundaries_periodic(mesh_pt, 0, 1, 0); // x
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
      else if(mesh_name == "single-element" && nnode1d == 2)
        {
          mesh_pt = new SingleElementMesh<QMicromagElement<2,2> >(time_stepper_pt);
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
                                        1, nnode1d);
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

          mesh_pt = new SimpleCubicMesh<QMicromagElement<3, 2> >
            (5*this_nx, std::ceil(1.25*this_nx), 2, 500, 125, 3, time_stepper_pt);

          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "sqt_mumag4" && nnode1d == 2)
        {
          Mesh* qmesh_pt = llg_mesh_factory("sq_mumag4", refinement_level,
                                            time_stepper_pt,
                                            1, nnode1d);
          //??ds memory leak, fix? Can't delete this mesh or nodes will
          //go...

          // Convert to tet mesh
          TetMeshBase* tmesh_pt = new TetMeshBase;
          ElementFactoryFctPt factory_fpt =
            MeshCreationHelpers::new_element<TMicromagElement<3, 2> >;
          MeshCreationHelpers::brick2tet(*qmesh_pt, factory_fpt, *tmesh_pt);

          mesh_pt = tmesh_pt;
        }
      else if(mesh_name == "sq_cubeoid" && nnode1d == 2)
        {
          double lx = 1, ly = lx, lz = 3*lx;
          mesh_pt = new SimpleCubicMesh<QMicromagElement<3, 2> >
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
      else
        {
          throw OomphLibError("Unrecognised mesh name " + mesh_name,
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      // Scale the mesh as requested
      scale_mesh(scaling_factor, mesh_pt);

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
    /// dimension and shape of the element. Only works for
    Integral* variable_order_integrator_factory(const FiniteElement* const el_pt)
    {
      if((el_pt->nodal_dimension() == 2) && (el_pt->nvertex_node() == 3))
        {
          return new TVariableOrderGaussLegendre<1>;
        }
      else if((el_pt->nodal_dimension() == 2) && (el_pt->nvertex_node() == 4))
        {
          return new QVariableOrderGaussLegendre<1>;
        }
      else if((el_pt->nodal_dimension() == 3) && (el_pt->nvertex_node() == 4))
        {
          return new TVariableOrderGaussLegendre<2>;
        }
      else if((el_pt->nodal_dimension() == 3) && (el_pt->nvertex_node() == 8))
        {
          return new QVariableOrderGaussLegendre<2>;
        }
      else
        {
          std::string err("Cannot determine element type.\n");
          err += "Maybe it is a higher order element (NNODE_1D > 2)?\n";
          err += "Variable order quadratures are not supported for this case.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }



    /// \short Return a function which will create the appropriate BEM face
    /// element for the bulk element pointer given (should work for a
    /// pointer to any bulk element type i.e., field or llg).
    BEMElementFactoryFctPt bem_element_factory_factory
    (const FiniteElement* bulk_ele_pt)
    {
      if(dynamic_cast<const TElement<2, 2>*>(bulk_ele_pt) != 0)
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
      if(dynamic_cast<const TMicromagElement<2, 2>*>(bulk_ele_pt) != 0)
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

      else
        {
          throw OomphLibError("Unrecognised element type",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

    }
  }
}
